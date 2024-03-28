library(tidyverse)
library(ggpmisc)
library(fuzzyjoin)

#find hanson 2020 data at:  https://www.mdpi.com/2218-1989/10/1/34/s1

#import supplementary file 1 and 2 and tidy data to R

hanson <- hanson_2020Supplementary_File_1 %>%
  select( -RI, -CAS, -PUBCHEM, -MASS, -CHEMSPIDER) %>%
  rbind(Supplementary_File_2)

hansonlong <- hanson %>%
  rename(median =`Ratio of medians (Control/Patient)` ) %>%
  pivot_longer(-c(1:9)) %>%
  rename(hmdb = `Group   HMDB_ID` )

#get hanson 2022 data and tidy it up

Hanson22_datasetc <- read_excel("~/Downloads/jciinsight-7-157621-s149_corrected(1).xlsx",
                               sheet = "ScaledImpDataZeroDrug&Tobacco")
h22basec <- Hanson22_datasetc %>% pivot_longer(cols = -c(1:16)) %>% mutate(id = row_number()) %>%  filter(id>2828)


h22_case_controlc<- Hanson22_datasetc %>% pivot_longer(cols = -c(1:16)) %>% mutate(id = row_number()) %>%  filter(id<2828) %>% filter(value %in% c("Control", "CFS")) %>% select( name, casecontrol =value)


h22prepostc <- Hanson22_datasetc %>% pivot_longer(cols = -c(1:16)) %>% mutate(id = row_number()) %>%  filter(id<2828) %>%
  filter(value %in% c("D1-PRE","D1-POST","D2-PRE","D2-POST")) %>% select(name, prepost = value)


h22sexc <- Hanson22_datasetc %>% pivot_longer(cols = -c(1:16)) %>% mutate(id = row_number()) %>%  filter(id<2828) %>%
  filter(value %in% c( "F", "M")) %>% select(name, sex = value)

h22bellc <- Hanson22_datasetc %>% pivot_longer(cols = -c(1:16)) %>% mutate(id = row_number()) %>%  filter(id<2828) %>%
  filter(value %in% c( 10,20,30,40,50,60,70,80,90,100)) %>% select(name, bell = value)

#make core data for working on:
h22finalc <- h22basec %>% left_join(h22_case_controlc, by = "name") %>% left_join(h22prepostc, by = "name") %>% left_join(h22sexc, by = "name") %>%
  left_join(h22bellc, by = "name") %>% filter(str_detect(name, "...")) %>% mutate(name = str_trunc(name, side = "right", 7, ellipsis = "")) %>%
  mutate(value = as.numeric(value)) %>% rename(substance = `...2`)


#manipulate 2022 data so it can be joined with 2020 data

h22forjoining <- h22finalc %>% filter(prepost == "D1-POST") %>%
  group_by(substance, casecontrol) %>%
  mutate(median = median(value)) %>%
  ungroup() %>%
  select(substance,casecontrol,SUPER_PATHWAY = `...3`, median, `...13`) %>%
  distinct() %>%
  mutate(control_div_by_patient = median/lead(median)) %>% filter(casecontrol == "Control") %>%
  select(substance, hmdb = `...13`, control_div_by_patient) %>% filter(!is.na(substance))

#join 2020 data to 2022 data and plot them.

hansonlong %>%  mutate(group = name) %>%
  mutate(group =str_remove_all(group, "[0-9]")) %>%
  group_by(group, BIOCHEMICAL) %>%
  mutate(median = median(value)) %>%
  ungroup() %>%
  select(substance = BIOCHEMICAL, hmdb,group, SUPER_PATHWAY, median) %>% distinct() %>% filter(!is.na(hmdb)) %>% filter(group == "median") %>%
stringdist_left_join(h22forjoining, by ="substance",max_dist = 0, ignore_case = TRUE,  method = "lv") %>%
  arrange(desc(control_div_by_patient)) %>%
  filter(!is.na(control_div_by_patient)) %>%
  ggplot()+aes(x=median, y= control_div_by_patient )+
  geom_point(size = 3, alpha = .2)+scale_y_continuous(limits = c(0,5))+
  geom_smooth(method = "lm", se=F)+
  xlab("Hanson 2020")+ylab("Hanson 2022")+facet_wrap(~SUPER_PATHWAY)+
  stat_poly_eq(colour = "blue")+
  labs(title = "Medians in controls divided by medians in patients")+theme_minimal()
