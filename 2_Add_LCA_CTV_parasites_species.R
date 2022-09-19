newglobi <- FALSE
if (newglobi) {
  source("1_globidata_parasites_species")
} else {
  globidata_parasites_s <- readRDS("globi_parasites_species.Rds")
}

library(tidyverse)
library(taxonomizr)

#Taxon Path von Vulpes Vulpes herausfinden
vulpestaxonpath <- data.frame(getTaxonomy(getId("Vulpes vulpes")))

globidata_parasites_s %>% filter(superkingdom %in% 
                                 c("Bacteria","Eukaryota","Viruses")) -> 
  globi_parasites_species

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            s_parasitesLCA_v=unique(species)) %>%
  filter (nspecies==1) %>%
  bind_cols(s_parasitesLCA_l="species",s_parasitesLCA="1") %>% 
  select(qtax_species,interaction_s_parasites,s_parasitesLCA_v,s_parasitesLCA_l,
         s_parasitesLCA) -> s_speciesdata #182

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            s_parasitesLCA_v=unique(genus)) %>%
  filter (nspecies !=1 & ngenus==1) %>%
  bind_cols(s_parasitesLCA_l="genus",s_parasitesLCA="2") %>% 
  select(qtax_species,interaction_s_parasites,s_parasitesLCA_v,s_parasitesLCA_l,
         s_parasitesLCA) -> s_genusdata #22

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            s_parasitesLCA_v=unique(family)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam==1) %>%
  bind_cols(s_parasitesLCA_l="family",s_parasitesLCA="3") %>% 
  select(qtax_species,interaction_s_parasites,s_parasitesLCA_v,s_parasitesLCA_l,
         s_parasitesLCA) -> s_familydata #32

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            s_parasitesLCA_v=unique(order)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder==1) %>%
  bind_cols(s_parasitesLCA_l="order",s_parasitesLCA="4") %>% 
  select(qtax_species,interaction_s_parasites,s_parasitesLCA_v,s_parasitesLCA_l,
         s_parasitesLCA) -> s_orderdata #18

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            s_parasitesLCA_v=unique(class)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass==1) %>%
  bind_cols(s_parasitesLCA_l="class",s_parasitesLCA="5") %>% 
  select(qtax_species,interaction_s_parasites,s_parasitesLCA_v,s_parasitesLCA_l,
         s_parasitesLCA) -> s_classdata #66

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            s_parasitesLCA_v=unique(phylum)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass!=1 & 
            nphylum==1) %>%
  bind_cols(s_parasitesLCA_l="phylum",s_parasitesLCA="6") %>% 
  select(qtax_species,interaction_s_parasites,s_parasitesLCA_v,s_parasitesLCA_l,
         s_parasitesLCA) -> s_phylumdata #37

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            s_parasitesLCA_v=unique(superkingdom)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass!=1 & 
            nphylum!=1 & nkingdom==1) %>%
  bind_cols(s_parasitesLCA_l="superkingdom",s_parasitesLCA="7") %>% 
  select(qtax_species,interaction_s_parasites,s_parasitesLCA_v,s_parasitesLCA_l,
         s_parasitesLCA) -> s_superkingdomdata #57

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom))) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass!=1 
          & nphylum!=1 & nkingdom>1) %>%
  bind_cols(s_parasitesLCA_l="Bac + Euk",
            s_parasitesLCA_v=">superkingdom",s_parasitesLCA="8") %>% 
  select(qtax_species,interaction_s_parasites,s_parasitesLCA_v,s_parasitesLCA_l,
         s_parasitesLCA) -> s_nolcadata #3

s_parasites_LCA <- 
  as.data.frame(bind_rows(s_speciesdata, s_genusdata,s_familydata,s_orderdata,
                          s_classdata,s_phylumdata,s_superkingdomdata,
                          s_nolcadata)) #417

#--- Closest to vulpes level + value finden

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes")) %>%
  filter(nvulpesspecies>=1) %>%
  bind_cols(s_parasitesCTV_l="species",s_parasitesCTV_v="Vulpes vulpes",
            s_parasitesCTV="1") %>%
  select(qtax_species,s_parasitesCTV_v,s_parasitesCTV_l,s_parasitesCTV) ->
  s_vspecies #7

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasite=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes")) %>%
  filter(nvulpesspecies %in% c(NA,0) & nvulpesgenus>=1) %>%
  bind_cols(s_parasitesCTV_l="genus",s_parasitesCTV_v="Vulpes",
            s_parasitesCTV="2") %>%
  select(qtax_species,s_parasitesCTV_v,s_parasitesCTV_l,s_parasitesCTV) -> 
  s_vgenus #10

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam>=1) %>%
  bind_cols(s_parasitesCTV_l="family",s_parasitesCTV_v="Canidae",
            s_parasitesCTV="3") %>%
  select(qtax_species,s_parasitesCTV_v,s_parasitesCTV_l,s_parasitesCTV) -> 
  s_vfam #21

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpesorder>=1) %>%
  bind_cols(s_parasitesCTV_l="order",s_parasitesCTV_v="Carnivora",
            s_parasitesCTV="4") %>%
  select(qtax_species,s_parasitesCTV_v,s_parasitesCTV_l,s_parasitesCTV) -> 
  s_vorder #47

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpesorder %in% c(NA,0) & 
            nvulpesclass>=1) %>%
  bind_cols(s_parasitesCTV_l="class",s_parasitesCTV_v="Mammalia",
            s_parasitesCTV="5") %>%
  select(qtax_species,s_parasitesCTV_v,s_parasitesCTV_l,s_parasitesCTV) -> 
  s_vclass #140

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpesorder %in% c(NA,0) & 
            nvulpesclass %in% c(NA,0) & nvulpesphylum>=1) %>%
  bind_cols(s_parasitesCTV_l="phylum",s_parasitesCTV_v="Chordata",
            s_parasitesCTV="6") %>%
  select(qtax_species,s_parasitesCTV_v,s_parasitesCTV_l,s_parasitesCTV) -> 
  s_vphylum #20

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
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
  bind_cols(s_parasitesCTV_l="superkingdom",s_parasitesCTV_v="Eukaryota",
            s_parasitesCTV="7") %>%
  select(qtax_species,s_parasitesCTV_v,s_parasitesCTV_l,s_parasitesCTV) -> 
  s_vkingdom #172

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
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
  bind_cols(s_parasitesCTV_l="No level CTV",
            s_parasitesCTV_v="No",s_parasitesCTV="8") %>%
  select(qtax_species,s_parasitesCTV_v,s_parasitesCTV_l,s_parasitesCTV) -> 
  s_vnomatch #0

# no match bedeutet hier: Bacterium hat nur Interaktionen mit anderen Bacteria

s_parasites_CTV <- 
  data.frame(bind_rows(s_vspecies,s_vgenus,s_vfam,s_vorder,s_vclass,s_vphylum,
                       s_vkingdom,s_vnomatch)) #417

# LCA und closest to vulpes verbinden

s_parasites <- merge(s_parasites_LCA,s_parasites_CTV)


saveRDS(s_parasites,"s_parasites.Rds")
