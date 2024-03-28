library(ggrepel)
library(tidyverse)
library(ggpmisc)
library(readxl)


#A map of metabolic phenotypes in patients with myalgic encephalomyelitis/chronic fatigue syndrome 2021
#get data
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8409979/bin/jciinsight-6-149217-s262.xlsx

  

flugefinal <-fluge %>% pivot_longer(-c(1:8)) %>%
  filter(!is.na(value)) %>%
  mutate(name = if_else(name == "MEall...9", "All_foldchange", name),
         name = if_else(name == "M1...10", "M1_foldchange", name),
         name = if_else(name == "M2...11", "M2_foldchange", name),
         name = if_else(name == "M3...12", "M3_foldchange", name),
         name = if_else(name == "MEall...14", "All_pv", name),
         name = if_else(name == "M1...15", "M1_pv", name),
         name = if_else(name == "M2...16", "M2_pv", name),
         name = if_else(name == "M3...17", "M3_pv", name),
         name = if_else(name == "MEall...19", "All_qv", name),
         name = if_else(name == "M1...20", "M1_qv", name),
         name = if_else(name == "M2...21", "M2_qv", name),
         name = if_else(name == "M3...22", "M3_qv", name)) %>%
  separate(name, into= c("subgroup", "measurement"), sep = "_")

#Get naviaux data


#get naviaux data for males and females and metadata on metabolites

#naviaux data source: https://www.metabolomicsworkbench.org/data/DRCCStudySummary.php?Mode=SetupRawDataDownload&StudyID=ST000450

#make female data long


femalelong<- female_CFS_metabolites_HMDB %>%
  pivot_longer(-c(`Sample Name`, Group, Sex)) %>%
  mutate(Sex  =if_else(Sex == "female", "Female", Sex))


#make male data long

malelong <- male_CFS_metabolites_HMDB %>%
  pivot_longer(-c(`Sample Name`, Group, Sex))

#bind together and make ready for joining.

Naviauxdata <- femalelong %>% rbind(malelong)


 navforjoining <- Naviauxdata %>% group_by(name, Group) %>%
  mutate(mean = mean(value)) %>%
  ungroup() %>%
  select(name,Group, mean) %>%
  distinct() %>%
  pivot_wider(names_from ="Group", values_from = "mean") %>%

  mutate(ratio = CFS/control) %>%
  filter(!is.na(ratio )) %>%
  arrange(desc(ratio))

#get list of naviaux metabolites to figure out which ones are in common with hanson.

fhmdb<- female_CFS_metabolites_metaID_8_19_16 %>% rbind(female_CFS_metabolites_metaID_8_19_16) %>% select(`Metabolite Name`, `HMDB ID`)


mhmdb<- male_CFS_metabolites_metaID_8_19_16 %>% rbind(male_CFS_metabolites_metaID_8_19_16) %>% select(`Metabolite Name`=`Metabolites Name`, `HMDB ID`)

navhmdb <- rbind(fhmdb, mhmdb) %>% distinct() %>% rename(name = `Metabolite Name` )

navforjoining2 <- navforjoining %>% left_join(navhmdb, by = "name") %>% rename(HMDB_ID = `HMDB ID`)


#find other metabolites in common using other technique

extras <- flugefinal %>% filter(measurement == "foldchange") %>%
  filter(subgroup == "All") %>%  select(name= Met, value) %>%
  stringdist_left_join(navforjoining2, by ="name",max_dist = 01, ignore_case = TRUE,  method = "lv") %>% #not sure if this is a fair subset of all
  filter(!is.na(ratio)) %>% filter(is.na(HMDB_ID)) 


#join and plot

flugefinal %>% filter(measurement == "foldchange") %>%
  filter(subgroup == "All") %>%  select(name= Met,HMDB_ID, value) %>%
  left_join(navforjoining2, by ="HMDB_ID") %>%
  filter(!is.na(ratio)) %>%
  rbind(extras) %>%
  filter(!name.y %in% c("15-HETE","11-HETE"  ))%>%
  #filter(!str_detect(name.x, "nicot"  ))%>%
  distinct() %>%
  mutate(consistency = value/ratio) %>%
  ggplot()+aes(x=value, y=ratio)+geom_point()+
  geom_text_repel(aes(label = name.x), size = 1.8, max.overlaps = 30)+
  scale_x_log10()+
  scale_y_log10()+
  ylab("Naviaux 2017 measurement")+
  xlab("Fluge 2021 measurement")+
  labs(title = "Fluge Mella 2021 vs Naviaux 2017",
       caption = "Fluge, Mella et al, A map of metabolic phenotypes in patients with myalgic encephalomyelitis/chronic fatigue syndrome,\n and Naviaux 2017, Metabolic features of chronic fatigue syndrome",
       subtitle = "The consistent measurements aren't interesting and the interesting ones aren't consistent")+
  geom_abline(slope= 1, intercept =0, size = 8, colour = "gold", alpha = .2)+
  geom_abline(slope= 1, intercept =0, size = 3, colour = "gold", alpha = .3)+
  theme_bw()+
  geom_smooth(method = "lm", se=F)+
  stat_poly_eq(colour = "blue")
