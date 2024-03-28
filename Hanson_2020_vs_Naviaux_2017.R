
library(tidyverse)
library(ggpmisc)
library(fuzzyjoin)

#get data for Hanson 2020
hanson <- hanson_2020Supplementary_File_1 %>%
  select( -RI, -CAS, -PUBCHEM, -MASS, -CHEMSPIDER) %>%
  rbind(Supplementary_File_2)

hansonlong <- hanson %>%
  rename(median =`Ratio of medians (Control/Patient)` ) %>%
  pivot_longer(-c(1:9)) %>%
  rename(hmdb = `Group   HMDB_ID` )

#get data for Naviaux 2017


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



#make a second list of metabolites in common using another technique

extrashansonnaviaux <-hansonlong %>% filter(name == "median" ) %>%
  mutate(value = 1/value) %>%    #inverting because her ratio is controls/patients]
  arrange(desc(value)) %>% select(name = BIOCHEMICAL, HMDB_ID = hmdb, value) %>%
  stringdist_left_join(navforjoining2, by ="name",max_dist = 0, ignore_case = TRUE,  method = "lv") %>% #not sure if this is a fair subset of all
  filter(!is.na(ratio)) %>% filter(is.na(HMDB_ID.y)| HMDB_ID.y=="N/A") %>% select(-HMDB_ID.y, HMDB_ID = HMDB_ID.x)

#combine Hanson 2020 and Naviaux and plot

hansonlong %>% filter(name == "median" ) %>%
  mutate(value = 1/value) %>%    #inverting because her ratio is controls/patients]
  arrange(desc(value)) %>% select(name = BIOCHEMICAL, HMDB_ID = hmdb, value) %>%
  filter(!is.na(HMDB_ID)) %>%
left_join(navforjoining2, by ="HMDB_ID") %>%
  filter(!is.na(ratio)) %>%
  rbind(extrashansonnaviaux) %>%
  ggplot()+aes(x=value, y=ratio)+geom_point()+
  ggrepel::geom_text_repel(aes(label = name.x), size = 1.8, max.overlaps = 10)+
  scale_x_log10()+
  scale_y_log10()+
  ylab("Naviaux 2017 measurement, patients divided by controls")+
  xlab("Hanson 2020 measurement, patients divided by controls")+
  geom_smooth(method ='lm', se=F)+
  labs(title = "Naviaux 2017 and Hanson 2020 don't line up",
       subtitle = "The consistent measurements aren't interesting and the interesting ones aren't consistent")+
  geom_abline(slope= 1, intercept =0, size = 18, colour = "gold", alpha = .2)+
  geom_abline(slope= 1, intercept =0, size = 3, colour = "gold", alpha = .3)+
  theme_bw()+
  stat_poly_eq()
