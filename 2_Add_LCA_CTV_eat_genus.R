newglobi <- FALSE
if (newglobi) {
  source("1_globidata_eat_genus.R")
} else {
  globidata_eat <- readRDS("foxglobieat.Rds")
}

library(tidyverse)
library(taxonomizr)

#Taxon Path von Vulpes Vulpes herausfinden
vulpestaxonpath <- data.frame(getTaxonomy(getId("Vulpes vulpes")))

globidata_eat %>% filter(superkingdom %in% c("Bacteria","Eukaryota",
                                             "Viruses")) -> new2

globi_eat <- rename(new2,qtax_genus=qtax2)

globi_eat %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            g_eatLCA_v=unique(species)) %>%
  filter (nspecies==1) %>%
  bind_cols(g_eatLCA_l="species",g_eatLCA="1") %>% 
  select(qtax_genus,interaction_g_eat,g_eatLCA_v,g_eatLCA_l,g_eatLCA) -> 
  g_speciesdata2 #91

globi_eat %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            g_eatLCA_v=unique(genus)) %>%
  filter (nspecies !=1 & ngenus==1) %>%
  bind_cols(g_eatLCA_l="genus",g_eatLCA="2") %>% 
  select(qtax_genus,interaction_g_eat,g_eatLCA_v,g_eatLCA_l,g_eatLCA) -> 
  g_genusdata2 #7

globi_eat %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            g_eatLCA_v=unique(family)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam==1) %>%
  bind_cols(g_eatLCA_l="family",g_eatLCA="3") %>% 
  select(qtax_genus,interaction_g_eat,g_eatLCA_v,g_eatLCA_l,g_eatLCA) -> 
  g_familydata2 #4

globi_eat %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            g_eatLCA_v=unique(order)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder==1) %>%
  bind_cols(g_eatLCA_l="order",g_eatLCA="4") %>% 
  select(qtax_genus,interaction_g_eat,g_eatLCA_v,g_eatLCA_l,g_eatLCA) -> 
  g_orderdata2 #11

globi_eat %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            g_eatLCA_v=unique(class)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass==1) %>%
  bind_cols(g_eatLCA_l="class",g_eatLCA="5") %>% 
  select(qtax_genus,interaction_g_eat,g_eatLCA_v,g_eatLCA_l,g_eatLCA) -> 
  g_classdata2 #15

globi_eat %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            g_eatLCA_v=unique(phylum)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass!=1 & 
            nphylum==1) %>%
  bind_cols(g_eatLCA_l="phylum",g_eatLCA="6") %>% 
  select(qtax_genus,interaction_g_eat,g_eatLCA_v,g_eatLCA_l,g_eatLCA) -> 
  g_phylumdata2 #20

globi_eat %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            g_eatLCA_v=unique(superkingdom)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass!=1 & 
            nphylum!=1 & nkingdom==1) %>%
  bind_cols(g_eatLCA_l="superkingdom",g_eatLCA="7") %>% 
  select(qtax_genus,interaction_g_eat,g_eatLCA_v,g_eatLCA_l,g_eatLCA) -> 
  g_superkingdomdata2 #249

globi_eat %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom))) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass!=1 & 
            nphylum!=1 & nkingdom>1) %>%
  bind_cols(g_eatLCA_l="Bac + Euk",g_eatLCA_v=">superkingdom",g_eatLCA="8") %>% 
  select(qtax_genus,interaction_g_eat,g_eatLCA_v,g_eatLCA_l,g_eatLCA) -> 
  g_nolcadata2 #27

g_eat_LCA <- 
  as.data.frame(bind_rows(g_speciesdata2, g_genusdata2,g_familydata2,
                          g_orderdata2, g_classdata2, g_phylumdata2,
                          g_superkingdomdata2,g_nolcadata2)) #424

#--- Closest to vulpes level + value finden

globi_eat %>% 
  group_by(qtax_genus) %>% 
  summarize(n_eat=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes")) %>%
  filter(nvulpesspecies>=1) %>%
  bind_cols(g_eatCTV_l="species",g_eatCTV_v="Vulpes vulpes",g_eatCTV="1") %>%
  select(qtax_genus,g_eatCTV_v,g_eatCTV_l,g_eatCTV) -> g_vspecies2 #0

globi_eat %>% 
  group_by(qtax_genus) %>% 
  summarize(n_eat=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes")) %>%
  filter(nvulpesspecies %in% c(NA,0) & nvulpesgenus>=1) %>%
  bind_cols(g_eatCTV_l="genus",g_eatCTV_v="Vulpes",g_eatCTV="2") %>%
  select(qtax_genus,g_eatCTV_v,g_eatCTV_l,g_eatCTV) -> g_vgenus2 #11

globi_eat %>% 
  group_by(qtax_genus) %>% 
  summarize(n_eat=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam>=1) %>%
  bind_cols(g_eatCTV_l="family",g_eatCTV_v="Canidae",g_eatCTV="3") %>%
  select(qtax_genus,g_eatCTV_v,g_eatCTV_l,g_eatCTV) -> g_vfam2 #5

globi_eat %>% 
  group_by(qtax_genus) %>% 
  summarize(n_eat=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpesorder>=1) %>%
  bind_cols(g_eatCTV_l="order",g_eatCTV_v="Carnivora",g_eatCTV="4") %>%
  select(qtax_genus,g_eatCTV_v,g_eatCTV_l,g_eatCTV) -> g_vorder2 #36

globi_eat %>% 
  group_by(qtax_genus) %>% 
  summarize(n_eat=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpesorder %in% c(NA,0) & 
            nvulpesclass>=1) %>%
  bind_cols(g_eatCTV_l="class",g_eatCTV_v="Mammalia",g_eatCTV="5") %>%
  select(qtax_genus,g_eatCTV_v,g_eatCTV_l,g_eatCTV) -> g_vclass2 #49

globi_eat %>% 
  group_by(qtax_genus) %>% 
  summarize(n_eat=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpesorder %in% c(NA,0) & 
            nvulpesclass %in% c(NA,0) & nvulpesphylum>=1) %>%
  bind_cols(g_eatCTV_l="phylum",g_eatCTV_v="Chordata",g_eatCTV="6") %>%
  select(qtax_genus,g_eatCTV_v,g_eatCTV_l,g_eatCTV) -> g_vphylum2 #111

globi_eat %>% 
  group_by(qtax_genus) %>% 
  summarize(n_eat=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata"),
            nvulpeskingdom = sum(superkingdom=="Eukaryota")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpesorder %in% c(NA,0) & 
            nvulpesclass %in% c(NA,0) & nvulpesphylum %in% c(NA,0) &
            nvulpeskingdom>=1) %>%
  bind_cols(g_eatCTV_l="superkingdom",g_eatCTV_v="Eukaryota",g_eatCTV="7") %>%
  select(qtax_genus,g_eatCTV_v,g_eatCTV_l,g_eatCTV) -> g_vkingdom2 #208

globi_eat %>% 
  group_by(qtax_genus) %>% 
  summarize(n_eat=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata"),
            nvulpeskingdom = sum(superkingdom=="Eukaryota")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpesorder %in% c(NA,0) & 
            nvulpesclass %in% c(NA,0) & nvulpesphylum %in% c(NA,0) &
            nvulpeskingdom %in% c(NA,0)) %>%
  bind_cols(g_eatCTV_l="No level ctv",g_eatCTV_v="No",g_eatCTV="8") %>%
  select(qtax_genus,g_eatCTV_v,g_eatCTV_l,g_eatCTV) -> g_vnomatch2 #4

# no match bedeutet hier: Bacterium hat nur Interaktionen mit anderen Bacteria

g_eat_CTV <- 
  data.frame(bind_rows(g_vspecies2,g_vgenus2,g_vfam2,g_vorder2,g_vclass2,
                       g_vphylum2,g_vkingdom2,g_vnomatch2)) #424

# LCA und closest to vulpes verbinden
g_eat <- merge(g_eat_LCA,g_eat_CTV)

saveRDS(g_eat,"g_eat.Rds")
