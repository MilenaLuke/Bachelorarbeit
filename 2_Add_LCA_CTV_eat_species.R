newglobi <- FALSE
if (newglobi) {
  source("1_globidata_eat_species.R")
} else {
  globidata_eat_s <- readRDS("globi_species_eat.Rds")
}

library(tidyverse)
library(taxonomizr)

#Taxon Path von Vulpes Vulpes herausfinden
vulpestaxonpath <- data.frame(getTaxonomy(getId("Vulpes vulpes")))

globidata_eat_s %>% filter(superkingdom %in% c("Bacteria","Eukaryota",
                                             "Viruses")) -> globi_eat_s

globi_eat_s %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            s_eatLCA_v=unique(species)) %>%
  filter (nspecies==1) %>%
  bind_cols(s_eatLCA_l="species",s_eatLCA="1") %>% 
  select(qtax_species,interaction_s_eat,s_eatLCA_v,s_eatLCA_l,s_eatLCA) -> 
  s_speciesdata2 #39

globi_eat_s %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            s_eatLCA_v=unique(genus)) %>%
  filter (nspecies !=1 & ngenus==1) %>%
  bind_cols(s_eatLCA_l="genus",s_eatLCA="2") %>% 
  select(qtax_species,interaction_s_eat,s_eatLCA_v,s_eatLCA_l,s_eatLCA) -> 
  s_genusdata2 #5

globi_eat_s %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            s_eatLCA_v=unique(family)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam==1) %>%
  bind_cols(s_eatLCA_l="family",s_eatLCA="3") %>% 
  select(qtax_species,interaction_s_eat,s_eatLCA_v,s_eatLCA_l,s_eatLCA) -> 
  s_familydata2 #3

globi_eat_s %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            s_eatLCA_v=unique(order)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder==1) %>%
  bind_cols(s_eatLCA_l="order",s_eatLCA="4") %>% 
  select(qtax_species,interaction_s_eat,s_eatLCA_v,s_eatLCA_l,s_eatLCA) -> 
  s_orderdata2 #4

globi_eat_s %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            s_eatLCA_v=unique(class)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass==1) %>%
  bind_cols(s_eatLCA_l="class",s_eatLCA="5") %>% 
  select(qtax_species,interaction_s_eat,s_eatLCA_v,s_eatLCA_l,s_eatLCA) -> 
  s_classdata2 #12

globi_eat_s %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            s_eatLCA_v=unique(phylum)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass!=1 & 
            nphylum==1) %>%
  bind_cols(s_eatLCA_l="phylum",s_eatLCA="6") %>% 
  select(qtax_species,interaction_s_eat,s_eatLCA_v,s_eatLCA_l,s_eatLCA) -> 
  s_phylumdata2 #7

globi_eat_s %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            s_eatLCA_v=unique(superkingdom)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass!=1 & 
            nphylum!=1 & nkingdom==1) %>%
  bind_cols(s_eatLCA_l="superkingdom",s_eatLCA="7") %>% 
  select(qtax_species,interaction_s_eat,s_eatLCA_v,s_eatLCA_l,s_eatLCA) -> 
  s_superkingdomdata2 #129

globi_eat_s %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom))) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass!=1 & 
            nphylum!=1 & nkingdom>1) %>%
  bind_cols(s_eatLCA_l="Bac + Euk",s_eatLCA_v=">superkingdom",s_eatLCA="8") %>% 
  select(qtax_species,interaction_s_eat,s_eatLCA_v,s_eatLCA_l,s_eatLCA) -> 
  s_nolcadata2 #23

s_eat_LCA <- 
  as.data.frame(bind_rows(s_speciesdata2, s_genusdata2,s_familydata2,
                          s_orderdata2, s_classdata2, s_phylumdata2,
                          s_superkingdomdata2,s_nolcadata2)) #222

#--- Closest to vulpes level + value finden

globi_eat_s %>% 
  group_by(qtax_species) %>% 
  summarize(s_n_eat=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes")) %>%
  filter(nvulpesspecies>=1) %>%
  bind_cols(s_eatCTV_l="species",s_eatCTV_v="Vulpes vulpes",s_eatCTV="1") %>%
  select(qtax_species,s_eatCTV_v,s_eatCTV_l,s_eatCTV) -> s_vspecies2 #0

globi_eat_s %>% 
  group_by(qtax_species) %>% 
  summarize(s_n_eat=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes")) %>%
  filter(nvulpesspecies %in% c(NA,0) & nvulpesgenus>=1) %>%
  bind_cols(s_eatCTV_l="genus",s_eatCTV_v="Vulpes",s_eatCTV="2") %>%
  select(qtax_species,s_eatCTV_v,s_eatCTV_l,s_eatCTV) -> s_vgenus2 #13

globi_eat_s %>% 
  group_by(qtax_species) %>% 
  summarize(s_n_eat=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam>=1) %>%
  bind_cols(s_eatCTV_l="family",s_eatCTV_v="Canidae",s_eatCTV="3") %>%
  select(qtax_species,s_eatCTV_v,s_eatCTV_l,s_eatCTV) -> s_vfam2 #1

globi_eat_s %>% 
  group_by(qtax_species) %>% 
  summarize(s_n_eat=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpesorder>=1) %>%
  bind_cols(s_eatCTV_l="order",s_eatCTV_v="Carnivora",s_eatCTV="4") %>%
  select(qtax_species,s_eatCTV_v,s_eatCTV_l,s_eatCTV) -> s_vorder2 #12

globi_eat_s %>% 
  group_by(qtax_species) %>% 
  summarize(s_n_eat=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpesorder %in% c(NA,0) & 
            nvulpesclass>=1) %>%
  bind_cols(s_eatCTV_l="class",s_eatCTV_v="Mammalia",s_eatCTV="5") %>%
  select(qtax_species,s_eatCTV_v,s_eatCTV_l,s_eatCTV) -> s_vclass2 #43

globi_eat_s %>% 
  group_by(qtax_species) %>% 
  summarize(s_n_eat=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpesorder %in% c(NA,0) & 
            nvulpesclass %in% c(NA,0) & nvulpesphylum>=1) %>%
  bind_cols(s_eatCTV_l="phylum",s_eatCTV_v="Chordata",s_eatCTV="6") %>%
  select(qtax_species,s_eatCTV_v,s_eatCTV_l,s_eatCTV) -> s_vphylum2 #45

globi_eat_s %>% 
  group_by(qtax_species) %>% 
  summarize(s_n_eat=n(),
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
  bind_cols(s_eatCTV_l="superkingdom",s_eatCTV_v="Eukaryota",s_eatCTV="7") %>%
  select(qtax_species,s_eatCTV_v,s_eatCTV_l,s_eatCTV) -> s_vkingdom2 #105

globi_eat_s %>% 
  group_by(qtax_species) %>% 
  summarize(s_n_eat=n(),
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
  bind_cols(s_eatCTV_l="No level ctv",s_eatCTV_v="No",s_eatCTV="8") %>%
  select(qtax_species,s_eatCTV_v,s_eatCTV_l) -> s_vnomatch2 #3

# no match bedeutet hier: Bacterium hat nur Interaktionen mit anderen Bacteria

s_eat_CTV <- 
  data.frame(bind_rows(s_vspecies2,s_vgenus2,s_vfam2,s_vorder2,s_vclass2,
                       s_vphylum2,s_vkingdom2,s_vnomatch2)) #222

# LCA und closest to vulpes verbinden
s_eat <- merge(s_eat_LCA,s_eat_CTV)

saveRDS(s_eat,"s_eat.Rds")
