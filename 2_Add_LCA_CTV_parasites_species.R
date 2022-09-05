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
  summarize(s_nparasites=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            s_parasitesLCA_v=unique(species)) %>%
  filter (nspecies==1) %>%
  bind_cols(s_parasitesLCA_l="species") %>% 
  select(qtax_species,s_nparasites,s_parasitesLCA_v,s_parasitesLCA_l) -> 
  s_speciesdata #182

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(s_nparasites=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            s_parasitesLCA_v=unique(genus)) %>%
  filter (nspecies !=1 & ngenus==1) %>%
  bind_cols(s_parasitesLCA_l="genus") %>% 
  select(qtax_species,s_nparasites,s_parasitesLCA_v,s_parasitesLCA_l) -> 
  s_genusdata #22

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(s_nparasites=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            s_parasitesLCA_v=unique(family)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam==1) %>%
  bind_cols(s_parasitesLCA_l="family") %>% 
  select(qtax_species,s_nparasites,s_parasitesLCA_v,s_parasitesLCA_l) -> 
  s_familydata #32

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(s_nparasites=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            s_parasitesLCA_v=unique(order)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder==1) %>%
  bind_cols(s_parasitesLCA_l="order") %>% 
  select(qtax_species,s_nparasites,s_parasitesLCA_v,s_parasitesLCA_l) -> 
  s_orderdata #18

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(s_nparasites=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            s_parasitesLCA_v=unique(class)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass==1) %>%
  bind_cols(s_parasitesLCA_l="class") %>% 
  select(qtax_species,s_nparasites,s_parasitesLCA_v,s_parasitesLCA_l) -> 
  s_classdata #66

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(s_nparasites=n(),
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
  bind_cols(s_parasitesLCA_l="phylum") %>% 
  select(qtax_species,s_nparasites,s_parasitesLCA_v,s_parasitesLCA_l) -> 
  s_phylumdata #37

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(s_nparasites=n(),
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
  bind_cols(s_parasitesLCA_l="superkingdom") %>% 
  select(qtax_species,s_nparasites,s_parasitesLCA_v,s_parasitesLCA_l) -> 
  s_superkingdomdata #57

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(s_nparasites=n(),
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
            s_parasitesLCA_v=">superkingdom") %>% 
  select(qtax_species,s_nparasites,s_parasitesLCA_v,s_parasitesLCA_l) -> 
  s_nolcadata #3

s_parasites_LCA <- 
  as.data.frame(bind_rows(s_speciesdata, s_genusdata,s_familydata,s_orderdata,
                          s_classdata,s_phylumdata,s_superkingdomdata,
                          s_nolcadata)) #417

#--- Closest to vulpes level + value finden

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(s_nparasites=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes")) %>%
  filter(nvulpesspecies>=1) %>%
  bind_cols(s_parasitesCTV_l="species",s_parasitesCTV_v="Vulpes vulpes") %>%
  select(qtax_species,s_parasitesCTV_v,s_parasitesCTV_l) ->s_vspecies #7

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(s_nparasites=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes")) %>%
  filter(nvulpesspecies %in% c(NA,0) & nvulpesgenus>=1) %>%
  bind_cols(s_parasitesCTV_l="genus",s_parasitesCTV_v="Vulpes") %>%
  select(qtax_species,s_parasitesCTV_v,s_parasitesCTV_l) -> s_vgenus #10


globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(s_nparasites=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam>=1) %>%
  bind_cols(s_parasitesCTV_l="family",s_parasitesCTV_v="Canidae") %>%
  select(qtax_species,s_parasitesCTV_v,s_parasitesCTV_l) -> s_vfam #21

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(s_nparasites=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpesorder>=1) %>%
  bind_cols(s_parasitesCTV_l="order",s_parasitesCTV_v="Carnivora") %>%
  select(qtax_species,s_parasitesCTV_v,s_parasitesCTV_l) -> s_vorder #47

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(s_nparasites=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpesorder %in% c(NA,0) & 
            nvulpesclass>=1) %>%
  bind_cols(s_parasitesCTV_l="class",s_parasitesCTV_v="Mammalia") %>%
  select(qtax_species,s_parasitesCTV_v,s_parasitesCTV_l) -> s_vclass #140

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(s_nparasites=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpesorder %in% c(NA,0) & 
            nvulpesclass %in% c(NA,0) & nvulpesphylum>=1) %>%
  bind_cols(s_parasitesCTV_l="phylum",s_parasitesCTV_v="Chordata") %>%
  select(qtax_species,s_parasitesCTV_v,s_parasitesCTV_l) -> s_vphylum #20

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(s_nparasites=n(),
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
  bind_cols(s_parasitesCTV_l="superkingdom",s_parasitesCTV_v="Eukaryota") %>%
  select(qtax_species,s_parasitesCTV_v,s_parasitesCTV_l) -> s_vkingdom #172

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(s_nparasites=n(),
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
            s_parasitesCTV_v="No") %>%
  select(qtax_species,s_parasitesCTV_v,s_parasitesCTV_l) -> s_vnomatch #0

# no match bedeutet hier: Bacterium hat nur Interaktionen mit anderen Bacteria

s_parasites_CTV <- 
  data.frame(bind_rows(s_vspecies,s_vgenus,s_vfam,s_vorder,s_vclass,s_vphylum,
                       s_vkingdom,s_vnomatch)) #417

# LCA und closest to vulpes verbinden

s_parasites <- merge(s_parasites_LCA,s_parasites_CTV)


saveRDS(s_parasites,"s_parasites.Rds")
