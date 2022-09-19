newglobi <- FALSE
if (newglobi) {
  source("1_globidata_parasites_family")
} else {
  globidata_parasites_fam2 <- readRDS("globi_parasites_family.Rds")
}

library(tidyverse)
library(taxonomizr)

#Taxon Path von Vulpes Vulpes herausfinden
vulpestaxonpath <- data.frame(getTaxonomy(getId("Vulpes vulpes")))

globidata_parasites_fam2 %>% filter(superkingdom %in% 
                                   c("Bacteria","Eukaryota","Viruses")) -> 
  globi_parasites_fam

globi_parasites_fam %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasite=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            f_parasitesLCA_v=unique(species)) %>%
  filter (nspecies==1) %>%
  bind_cols(f_parasitesLCA_l="species",f_parasitesLCA="1") %>% 
  select(qtax_family,interaction_f_parasite,f_parasitesLCA_v,f_parasitesLCA_l,
         f_parasitesLCA) -> f_speciesdata #182

globi_parasites_fam %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasite=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            f_parasitesLCA_v=unique(genus)) %>%
  filter (nspecies !=1 & ngenus==1) %>%
  bind_cols(f_parasitesLCA_l="genus",f_parasitesLCA="2") %>% 
  select(qtax_family,interaction_f_parasite,f_parasitesLCA_v,f_parasitesLCA_l,
         f_parasitesLCA) -> f_genusdata #22

globi_parasites_fam %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasite=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            f_parasitesLCA_v=unique(family)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam==1) %>%
  bind_cols(f_parasitesLCA_l="family",f_parasitesLCA="3") %>% 
  select(qtax_family,interaction_f_parasite,f_parasitesLCA_v,f_parasitesLCA_l,
         f_parasitesLCA) -> f_familydata #32

globi_parasites_fam %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasite=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            f_parasitesLCA_v=unique(order)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder==1) %>%
  bind_cols(f_parasitesLCA_l="order",f_parasitesLCA="4") %>% 
  select(qtax_family,interaction_f_parasite,f_parasitesLCA_v,f_parasitesLCA_l,
         f_parasitesLCA) -> f_orderdata #18

globi_parasites_fam %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasite=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            f_parasitesLCA_v=unique(class)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass==1) %>%
  bind_cols(f_parasitesLCA_l="class",f_parasitesLCA="5") %>% 
  select(qtax_family,interaction_f_parasite,f_parasitesLCA_v,f_parasitesLCA_l,
         f_parasitesLCA) -> f_classdata #66

globi_parasites_fam %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasite=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            f_parasitesLCA_v=unique(phylum)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass!=1 & 
            nphylum==1) %>%
  bind_cols(f_parasitesLCA_l="phylum",f_parasitesLCA="6") %>% 
  select(qtax_family,interaction_f_parasite,f_parasitesLCA_v,f_parasitesLCA_l,f_parasitesLCA) -> 
  f_phylumdata #37

globi_parasites_fam %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasite=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            f_parasitesLCA_v=unique(superkingdom)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass!=1 & 
            nphylum!=1 & nkingdom==1) %>%
  bind_cols(f_parasitesLCA_l="superkingdom",f_parasitesLCA="7") %>% 
  select(qtax_family,interaction_f_parasite,f_parasitesLCA_v,f_parasitesLCA_l,f_parasitesLCA) -> 
  f_superkingdomdata #57

globi_parasites_fam %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasite=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom))) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass!=1 
          & nphylum!=1 & nkingdom>1) %>%
  bind_cols(f_parasitesLCA_l="Bac + Euk",
            f_parasitesLCA_v=">superkingdom",f_parasitesLCA="8") %>% 
  select(qtax_family,interaction_f_parasite,f_parasitesLCA_v,f_parasitesLCA_l,f_parasitesLCA) -> 
  f_nolcadata #3

f_parasites_LCA <- 
  as.data.frame(bind_rows(f_speciesdata, f_genusdata,f_familydata,f_orderdata,
                          f_classdata,f_phylumdata,f_superkingdomdata,
                          f_nolcadata)) #433

#--- Closest to vulpes level + value finden

globi_parasites_fam %>% 
  group_by(qtax_family) %>% 
  summarize(f_nparasites=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes")) %>%
  filter(nvulpesspecies>=1) %>%
  bind_cols(f_parasitesCTV_l="species",f_parasitesCTV_v="Vulpes vulpes",f_parasitesCTV="1") %>%
  select(qtax_family,f_parasitesCTV_v,f_parasitesCTV_l,f_parasitesCTV) ->f_vspecies #7

globi_parasites_fam %>% 
  group_by(qtax_family) %>% 
  summarize(f_nparasites=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes")) %>%
  filter(nvulpesspecies %in% c(NA,0) & nvulpesgenus>=1) %>%
  bind_cols(f_parasitesCTV_l="genus",f_parasitesCTV_v="Vulpes",f_parasitesCTV="2") %>%
  select(qtax_family,f_parasitesCTV_v,f_parasitesCTV_l,f_parasitesCTV) -> f_vgenus #10

globi_parasites_fam %>% 
  group_by(qtax_family) %>% 
  summarize(f_nparasites=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam>=1) %>%
  bind_cols(f_parasitesCTV_l="family",f_parasitesCTV_v="Canidae",f_parasitesCTV="3") %>%
  select(qtax_family,f_parasitesCTV_v,f_parasitesCTV_l,f_parasitesCTV) -> f_vfam #21

globi_parasites_fam %>% 
  group_by(qtax_family) %>% 
  summarize(f_nparasites=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpesorder>=1) %>%
  bind_cols(f_parasitesCTV_l="order",f_parasitesCTV_v="Carnivora",f_parasitesCTV="4") %>%
  select(qtax_family,f_parasitesCTV_v,f_parasitesCTV_l,f_parasitesCTV) -> f_vorder #47

globi_parasites_fam %>% 
  group_by(qtax_family) %>% 
  summarize(f_nparasites=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpesorder %in% c(NA,0) & 
            nvulpesclass>=1) %>%
  bind_cols(f_parasitesCTV_l="class",f_parasitesCTV_v="Mammalia",f_parasitesCTV="5") %>%
  select(qtax_family,f_parasitesCTV_v,f_parasitesCTV_l,f_parasitesCTV) -> f_vclass #140

globi_parasites_fam %>% 
  group_by(qtax_family) %>% 
  summarize(f_nparasites=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpesorder %in% c(NA,0) & 
            nvulpesclass %in% c(NA,0) & nvulpesphylum>=1) %>%
  bind_cols(f_parasitesCTV_l="phylum",f_parasitesCTV_v="Chordata",f_parasitesCTV="6") %>%
  select(qtax_family,f_parasitesCTV_v,f_parasitesCTV_l,f_parasitesCTV) -> f_vphylum #20

globi_parasites_fam %>% 
  group_by(qtax_family) %>% 
  summarize(f_nparasites=n(),
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
  bind_cols(f_parasitesCTV_l="superkingdom",f_parasitesCTV_v="Eukaryota",f_parasitesCTV="7") %>%
  select(qtax_family,f_parasitesCTV_v,f_parasitesCTV_l,f_parasitesCTV) -> f_vkingdom #172

globi_parasites_fam %>% 
  group_by(qtax_family) %>% 
  summarize(f_nparasites=n(),
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
  bind_cols(f_parasitesCTV_l="No level CTV",
            f_parasitesCTV_v="No",f_parasitesCTV="8") %>%
  select(qtax_family,f_parasitesCTV_v,f_parasitesCTV_l,f_parasitesCTV) -> f_vnomatch #0

# no match bedeutet hier: Bacterium hat nur Interaktionen mit anderen Bacteria

f_parasites_CTV <- 
  data.frame(bind_rows(f_vspecies,f_vgenus,f_vfam,f_vorder,f_vclass,f_vphylum,
                       f_vkingdom,f_vnomatch)) #433

# LCA und closest to vulpes verbinden

f_parasites <- merge(f_parasites_LCA,f_parasites_CTV)


saveRDS(f_parasites,"f_parasites.Rds")
