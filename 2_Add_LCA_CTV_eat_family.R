newglobi <- FALSE
if (newglobi) {
  source("1_globidata_eat_family.R")
} else {
  globidata_eat_f <- readRDS("globi_family_eat.Rds")
}

library(tidyverse)
library(taxonomizr)

#Taxon Path von Vulpes Vulpes herausfinden
vulpestaxonpath <- data.frame(getTaxonomy(getId("Vulpes vulpes")))

globidata_eat_f %>% filter(superkingdom %in% c("Bacteria","Eukaryota",
                                               "Viruses")) -> globi_eat_f

globi_eat_f %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            f_eatLCA_v=unique(species)) %>%
  filter (nspecies==1) %>%
  bind_cols(f_eatLCA_l="species",f_eatLCA="1") %>% 
  select(qtax_family,interaction_f_eat,f_eatLCA_v,f_eatLCA_l,f_eatLCA) -> 
  f_speciesdata2 #62

globi_eat_f %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            f_eatLCA_v=unique(genus)) %>%
  filter (nspecies !=1 & ngenus==1) %>%
  bind_cols(f_eatLCA_l="genus",f_eatLCA="2") %>% 
  select(qtax_family,interaction_f_eat,f_eatLCA_v,f_eatLCA_l,f_eatLCA) -> 
  f_genusdata2 #4

globi_eat_f %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            f_eatLCA_v=unique(family)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam==1) %>%
  bind_cols(f_eatLCA_l="family",f_eatLCA="3") %>% 
  select(qtax_family,interaction_f_eat,f_eatLCA_v,f_eatLCA_l,f_eatLCA) -> 
  f_familydata2 #1

globi_eat_f %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            f_eatLCA_v=unique(order)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder==1) %>%
  bind_cols(f_eatLCA_l="order",f_eatLCA="4") %>% 
  select(qtax_family,interaction_f_eat,f_eatLCA_v,f_eatLCA_l,f_eatLCA) -> 
  f_orderdata2 #3

globi_eat_f %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            f_eatLCA_v=unique(class)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass==1) %>%
  bind_cols(f_eatLCA_l="class",f_eatLCA="5") %>% 
  select(qtax_family,interaction_f_eat,f_eatLCA_v,f_eatLCA_l,f_eatLCA) -> 
  f_classdata2 #1

globi_eat_f %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            f_eatLCA_v=unique(phylum)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass!=1 & 
            nphylum==1) %>%
  bind_cols(f_eatLCA_l="phylum",f_eatLCA="6") %>% 
  select(qtax_family,interaction_f_eat,f_eatLCA_v,f_eatLCA_l,f_eatLCA) -> 
  f_phylumdata2 #16

globi_eat_f %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            f_eatLCA_v=unique(superkingdom)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass!=1 & 
            nphylum!=1 & nkingdom==1) %>%
  bind_cols(f_eatLCA_l="superkingdom",f_eatLCA="7") %>% 
  select(qtax_family,interaction_f_eat,f_eatLCA_v,f_eatLCA_l,f_eatLCA) -> 
  f_superkingdomdata2 #211

globi_eat_f %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom))) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass!=1 & 
            nphylum!=1 & nkingdom>1) %>%
  bind_cols(f_eatLCA_l="Bac + Euk",f_eatLCA_v=">superkingdom",f_eatLCA="8") %>% 
  select(qtax_family,interaction_f_eat,f_eatLCA_v,f_eatLCA_l,f_eatLCA) -> 
  f_nolcadata2 #20

f_eat_LCA <- 
  as.data.frame(bind_rows(f_speciesdata2, f_genusdata2,f_familydata2,
                          f_orderdata2, f_classdata2, f_phylumdata2,
                          f_superkingdomdata2,f_nolcadata2)) #328

#--- Closest to vulpes level + value finden

globi_eat_f %>% 
  group_by(qtax_family) %>% 
  summarize(f_n_eat=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes")) %>%
  filter(nvulpesspecies>=1) %>%
  bind_cols(f_eatCTV_l="species",f_eatCTV_v="Vulpes vulpes",f_eatCTV="1") %>%
  select(qtax_family,f_eatCTV_v,f_eatCTV_l,f_eatCTV) -> f_vspecies2 #0

globi_eat_f %>% 
  group_by(qtax_family) %>% 
  summarize(f_n_eat=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes")) %>%
  filter(nvulpesspecies %in% c(NA,0) & nvulpesgenus>=1) %>%
  bind_cols(f_eatCTV_l="genus",f_eatCTV_v="Vulpes",f_eatCTV="2") %>%
  select(qtax_family,f_eatCTV_v,f_eatCTV_l,f_eatCTV) -> f_vgenus2 #3

globi_eat_f %>% 
  group_by(qtax_family) %>% 
  summarize(f_n_eat=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam>=1) %>%
  bind_cols(f_eatCTV_l="family",f_eatCTV_v="Canidae",f_eatCTV="3") %>%
  select(qtax_family,f_eatCTV_v,f_eatCTV_l,f_eatCTV) -> f_vfam2 #6

globi_eat_f %>% 
  group_by(qtax_family) %>% 
  summarize(f_n_eat=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpesorder>=1) %>%
  bind_cols(f_eatCTV_l="order",f_eatCTV_v="Carnivora",f_eatCTV="4") %>%
  select(qtax_family,f_eatCTV_v,f_eatCTV_l,f_eatCTV) -> f_vorder2 #20

globi_eat_f %>% 
  group_by(qtax_family) %>% 
  summarize(f_n_eat=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpesorder %in% c(NA,0) & 
            nvulpesclass>=1) %>%
  bind_cols(f_eatCTV_l="class",f_eatCTV_v="Mammalia",f_eatCTV="5") %>%
  select(qtax_family,f_eatCTV_v,f_eatCTV_l,f_eatCTV) -> f_vclass2 #38

globi_eat_f %>% 
  group_by(qtax_family) %>% 
  summarize(f_n_eat=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpesorder %in% c(NA,0) & 
            nvulpesclass %in% c(NA,0) & nvulpesphylum>=1) %>%
  bind_cols(f_eatCTV_l="phylum",f_eatCTV_v="Chordata",f_eatCTV="6") %>%
  select(qtax_family,f_eatCTV_v,f_eatCTV_l,f_eatCTV) -> f_vphylum2 #100

globi_eat_f %>% 
  group_by(qtax_family) %>% 
  summarize(f_n_eat=n(),
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
  bind_cols(f_eatCTV_l="superkingdom",f_eatCTV_v="Eukaryota",f_eatCTV="7") %>%
  select(qtax_family,f_eatCTV_v,f_eatCTV_l,f_eatCTV) -> f_vkingdom2 #157

globi_eat_f %>% 
  group_by(qtax_family) %>% 
  summarize(f_n_eat=n(),
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
  bind_cols(f_eatCTV_l="No level ctv",f_eatCTV_v="No",f_eatCTV="8") %>%
  select(qtax_family,f_eatCTV_v,f_eatCTV_l,f_eatCTV) -> f_vnomatch2 #4

# no match bedeutet hier: Bacterium hat nur Interaktionen mit anderen Bacteria

f_eat_CTV <- 
  data.frame(bind_rows(f_vspecies2,f_vgenus2,f_vfam2,f_vorder2,f_vclass2,
                       f_vphylum2,f_vkingdom2,f_vnomatch2)) #328

# LCA und closest to vulpes verbinden
f_eat <- merge(f_eat_LCA,f_eat_CTV)

saveRDS(f_eat,"f_eat.Rds")
