library(readxl)
library(tidyverse)

#import dataset

Hanson22_dataset <- read_excel("~/Downloads/jciinsight-7-157621-s149_corrected.xlsx",
                               sheet = "ScaledImpDataZeroDrug&Tobacco")

#break dataset into useful subcomponents

h22base <- Hanson22_dataset %>% pivot_longer(cols = -c(1:16)) %>% 
mutate(id = row_number()) %>%  filter(id>2828)




h22_case_control<- Hanson22_dataset %>% pivot_longer(cols = -c(1:16)) %>% 
mutate(id = row_number()) %>%  filter(id<2828) %>% filter(value %in% c("Control", "CFS")) %>% 
select( name, casecontrol =value)


h22prepost <- Hanson22_dataset %>% pivot_longer(cols = -c(1:16)) %>% mutate(id = row_number()) %>%  filter(id<2828) %>%
  filter(value %in% c("D1-PRE","D1-POST","D2-PRE","D2-POST")) %>% select(name, prepost = value)


h22sex <- Hanson22_dataset %>% pivot_longer(cols = -c(1:16)) %>% mutate(id = row_number()) %>%  filter(id<2828) %>%
  filter(value %in% c( "F", "M")) %>% select(name, sex = value)

h22bell <- Hanson22_dataset %>% pivot_longer(cols = -c(1:16)) %>% mutate(id = row_number()) %>%  filter(id<2828) %>%
  filter(value %in% c( 10,20,30,40,50,60,70,80,90,100)) %>% select(name, bell = value)

#make core dataset for working on from subcomponents

h22final<-h22base %>% left_join(h22_case_control, by = "name") %>% 
left_join(h22prepost, by = "name") %>% 
left_join(h22sex, by = "name") %>%
  left_join(h22bell, by = "name") %>% 
filter(str_detect(name, "...")) %>% 
mutate(name = str_trunc(name, side = "right", 7, ellipsis = "")) %>%
  mutate(value = as.numeric(value)) %>% rename(substance = `...2`)

#make function that permits looking at each substance in turn


substancebars <- function(x){
  h22final %>% filter(substance == x) %>%mutate(marker = paste( casecontrol,sex, name)) %>%
    mutate(marker = str_remove_all(marker, "CU-1"))%>%

    mutate(prepost = factor(prepost, levels = c("D1-PRE","D1-POST","D2-PRE","D2-POST"))) %>%
    group_by(prepost, casecontrol) %>%
    mutate(mean = mean(value)) %>%
    mutate(median = median(value)) %>%
    ungroup() %>%
    mutate(marker2 = paste0( casecontrol,prepost)) %>%
    mutate(pval = if_else(prepost == "D1-PRE", t.test(value[marker2 == 'ControlD1-PRE'],
                                                      value[marker2 == 'CFSD1-PRE'])$p.value, 1)) %>%

    mutate(pval =if_else(prepost == "D1-POST", t.test(value[marker2 == 'ControlD1-POST'],
                                                      value[marker2 == 'CFSD1-POST'])$p.value, pval)) %>%

    mutate(pval =if_else(prepost == "D2-PRE", t.test(value[marker2 == 'ControlD2-PRE'],
                                                     value[marker2 == 'CFSD2-PRE'])$p.value, pval)) %>%

    mutate(pval =if_else(prepost == "D2-POST", t.test(value[marker2 == 'ControlD2-POST'],
                                                      value[marker2 == 'CFSD2-POST'])$p.value, pval)) %>%

    ggplot()+aes(x=marker)+
    #geom_text(aes(y=value, label = bell),size = 1.6, vjust = -.2)+
    geom_col(aes (y=value, fill = casecontrol, colour = sex, size = sex), show.legend = FALSE)+
    geom_line(aes(y=mean, group = casecontrol), colour = "black", show.legend = FALSE)+
    geom_line(aes(y=median, group = casecontrol), colour = "grey", show.legend = FALSE)+
    geom_text( show.legend = FALSE, aes(x=10, y=max(value)-.5, label = paste("p=",round(pval, 3))))+
    facet_wrap(~prepost, ncol=1)+   #+scale_y_log10()+
    scale_color_manual(values = c(NA, "grey10"))+
    scale_size_manual(values = c(0, .2))+
    theme_bw()+
    theme(axis.text.x = element_text(angle=90, size =6),
          axis.title = element_blank())+
    labs(title = paste("Measurement of", x, "at the four time points."),
         subtitle = "Patients in red at left, controls in blue at right, males outlined in grey at the end of each group. \nDark line is means, light grey line is medians, p-value applies to group means.\nPatients did two cardiopulmonary exercise tests separated by 24 hours, on Day 1 (D1) and Day 2 (D2). \nScientists took blood before (PRE) and after (POST) exercise for four measurements.",
         caption = "Data: Germain et al, 2022, Plasma metabolomics reveals disrupted response and recovery following maximal exercise in \nmyalgic encephalomyelitis/chronic fatigue syndrome. JCI Insight. Supplemental data table 1, scaled data. ")
}

#Demonstrations of function
substancebars("ADP")
substancebars("AMP")
substancebars("cotinine")
substancebars("1-methylnicotinamide")
substancebars("nicotinamide")
substancebars("hypoxanthine")
substancebars("xanthine")

#a few other ways to cut up the data


#bell score vs value, four facets, one for each timepoint, plus line of best fit. could add r^2

h22final %>% filter(substance == "nicotinamide")  %>%  ggplot()+aes(x=value, y= as.numeric(bell), colour =casecontrol)+
  geom_point(size = 4,  alpha = .6)+geom_smooth(method = "lm")+scale_x_log10()+facet_wrap(~prepost)

#line graph, one line for each patient over the time points. two facets, case and control.
h22final %>% filter(substance == "arginine")  %>% mutate(prepost = factor(prepost, levels = c("D1-PRE","D1-POST","D2-PRE","D2-POST"))) %>%
  ggplot()+aes(x=prepost, y= value, group = name, colour= casecontrol)+geom_line()+geom_point()+facet_wrap(~casecontrol,ncol=2 )

# Line graph, means over 4 time points, one line for cases, one for controls
h22final%>% filter(substance == "citrulline") %>%mutate(prepost = factor(prepost, levels = c("D1-PRE","D1-POST","D2-PRE","D2-POST"))) %>%
  group_by(substance, casecontrol, prepost) %>% mutate(mean = mean(value)) %>% ungroup() %>%
  ggplot()+aes(x=prepost, y=mean, colour = casecontrol, group = casecontrol)+geom_line()+facet_wrap(~substance)
