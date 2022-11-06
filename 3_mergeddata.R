globi_p_f <- FALSE
if (globi_p_f) {
  source("2_Add_LCA_CTV_parasites_family.R")
} else {
  f_parasites <- readRDS("family_parasites.Rds")
}

globi_p_g <- FALSE
if (globi_p_g) {
  source("2_Add_LCA_CTV_parasites_genus.R")
} else {
  g_parasites <- readRDS("genus_parasites.Rds")
}

globi_p_s <- FALSE
if (globi_p_s) {
  source("2_Add_LCA_CTV_parasites_species.R")
} else {
  s_parasites <- readRDS("species_parasites.Rds")
}

library(tidyverse)
library(taxonomizr)
library(ggplot2)

# f_parasites #237
# g_parasites #336
# s_parasites #328


write.csv(f_parasites, file="family_parasites.csv", row.names = FALSE)
write.csv(g_parasites, file="genus_parasites.csv", row.names = FALSE)
write.csv(s_parasites, file="species_parasites.csv", row.names = FALSE)


#-----
x <- merge(s_parasites,g_parasites,by.x="genus",by.y="genus",all.x = TRUE)

x %>%
  select(superkingdom.x,
         kingdom.x,
         phylum.x,
         subphylum.x,
         superclass.x,
         class.x,
         subclass.x,
         superorder.x,
         order.x,
         suborder.x,
         family.x,
         genus,
         species.x,
         interaction_s_parasites,
         s_parasitesLCA_v,
         s_parasitesLCA_l,
         s_parasitesLCA,
         s_parasitesCTV_v,
         s_parasitesCTV_l,
         s_parasitesCTV,
         interaction_g_parasites,
         g_parasitesLCA_v,
         g_parasitesLCA_l,
         g_parasitesLCA,
         g_parasitesCTV_v,
         g_parasitesCTV_l,
         g_parasitesCTV) %>%
  rename(superkingdom =superkingdom.x,
         kingdom=kingdom.x,
         phylum=phylum.x,
         subphylum=subphylum.x,
         superclass=superclass.x,
         class=class.x,
         subclass=subclass.x,
         superorder=superorder.x,
         order=order.x,
         suborder=suborder.x,
         family=family.x,
         species=species.x) -> species_genus_data

y <- merge(species_genus_data,f_parasites,by.x="family",by.y="family",all.x=TRUE)

y %>%
  select(superkingdom.x,
         kingdom.x,
         phylum.x,
         subphylum.x,
         superclass.x,
         class.x,
         subclass.x,
         superorder.x,
         order.x,
         suborder.x,
         genus.x,
         family,
         genus.x,
         species.x,
         interaction_s_parasites,
         s_parasitesLCA_v,
         s_parasitesLCA_l,
         s_parasitesLCA,
         s_parasitesCTV_v,
         s_parasitesCTV_l,
         s_parasitesCTV,
         interaction_g_parasites,
         g_parasitesLCA_v,
         g_parasitesLCA_l,
         g_parasitesLCA,
         g_parasitesCTV_v,
         g_parasitesCTV_l,
         g_parasitesCTV,
         interaction_f_parasites,
         f_parasitesLCA_v,
         f_parasitesLCA_l,
         f_parasitesLCA,
         f_parasitesCTV_v,
         f_parasitesCTV_l,
         f_parasitesCTV) %>%
  rename(superkingdom =superkingdom.x,
         kingdom=kingdom.x,
         phylum=phylum.x,
         subphylum=subphylum.x,
         superclass=superclass.x,
         class=class.x,
         subclass=subclass.x,
         superorder=superorder.x,
         order=order.x,
         suborder=suborder.x,
         genus=genus.x,
         species=species.x) -> mergeddata

mergeddata %>%
  transform(s_parasitesCTV=as.numeric(s_parasitesCTV)) %>%
  transform(s_parasitesLCA=as.numeric(s_parasitesLCA)) %>%
  transform(g_parasitesCTV=as.numeric(g_parasitesCTV)) %>%
  transform(g_parasitesLCA=as.numeric(g_parasitesLCA)) %>%
  transform(f_parasitesCTV=as.numeric(f_parasitesCTV)) %>%
  transform(f_parasitesLCA=as.numeric(f_parasitesLCA)) %>%
  filter(superkingdom=="Eukaryota",
         kingdom %in% c(NA,"Viridiplantae","Metazoa")) -> mergeddata2

# Einschätzung Parasit
mergeddata2 %>%
  filter(s_parasitesCTV<=4 & s_parasitesLCA <=4) %>%
  bind_cols(Parasit="Fuchsparasit") -> fuchsparasit #4

mergeddata2 %>%
  filter(s_parasitesCTV>=10 & s_parasitesLCA>=10) %>%
  bind_cols(Parasit="Kein Fuchsparasit") -> keinfuchsparasit1 #1

mergeddata2 %>%
  filter(s_parasitesCTV<=4 & s_parasitesLCA>=10) %>%
  bind_cols(Parasit="Parasit allgemein") -> parasit #12

mergeddata2 %>%
  filter(s_parasitesCTV>10 & s_parasitesLCA<=4) %>%
  bind_cols(Parasit="Parasit, aber kein Fuchsparasit") -> keinfuchsparasit2 #70

mergeddata2 %>%
  filter(s_parasitesCTV %in% c(5,6,7,8,9)) %>%
  bind_cols(Parasit="Vielleicht Parasit") -> maybeparasit1 #53
  
mergeddata2 %>%
  filter(s_parasitesLCA %in% c(5,6,7,8,9)) %>%
  bind_cols(Parasit="Vielleicht Parasit") -> maybeparasit2 #12

parasites_mergeddata <- 
  unique(data.frame(bind_rows(fuchsparasit,keinfuchsparasit1,keinfuchsparasit2,parasit,
                       maybeparasit1,maybeparasit2))) #151

parasites_mergeddata %>%
  filter(species!=s_parasitesLCA_v) -> parasites_data


# Versuche

parasites_data %>%
  filter(Parasit=="Vielleicht Parasit") -> test1







 ## -----------------------------------------------
# 1) Parasites plotten
plot_species_p <- ggplot(species_parasites,aes(s_parasitesCTV,s_parasitesLCA)) + 
  geom_point()

genus_parasites %>% filter(!is.na(g_parasitesCTV) & !is.na(g_parasitesLCA)) %>% 
  ggplot(aes(g_parasitesCTV,g_parasitesLCA,size=interaction_g_parasites)) + 
  geom_jitter()

genus_parasites %>% filter(g_parasitesCTV<4)


plot_family_p <- ggplot(family_parasites,aes(f_parasitesCTV,f_parasitesLCA)) + 
  geom_point()
