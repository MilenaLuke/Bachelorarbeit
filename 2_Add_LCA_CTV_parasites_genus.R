newglobi <- FALSE
if (newglobi) {
  source("1_globidata_parasites_genus")
} else {
  globidata_parasites <- readRDS("foxglobi.Rds")
}

library(tidyverse)
library(taxonomizr)

#Taxon Path von Vulpes Vulpes herausfinden
vulpestaxonpath <- data.frame(getTaxonomy(getId("Vulpes vulpes")))

globidata_parasites %>% filter(superkingdom %in% 
                                 c("Bacteria","Eukaryota","Viruses")) -> new
globi_parasites <- rename(new,qtax_genus = qtax)

globi_parasites %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            g_parasitesLCA_v=unique(species)) %>%
  filter (nspecies==1) %>%
  bind_cols(g_parasitesLCA_l="species",g_parasitesLCA="1") %>% 
  select(qtax_genus,interaction_g_parasites,g_parasitesLCA_v,g_parasitesLCA_l,
         g_parasitesLCA) -> g_speciesdata #161

globi_parasites %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            g_parasitesLCA_v=unique(genus)) %>%
  filter (nspecies !=1 & ngenus==1) %>%
  bind_cols(g_parasitesLCA_l="genus",g_parasitesLCA="2") %>% 
  select(qtax_genus,interaction_g_parasites,g_parasitesLCA_v,g_parasitesLCA_l,
         g_parasitesLCA) -> g_genusdata #24
              
globi_parasites %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            g_parasitesLCA_v=unique(family)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam==1) %>%
  bind_cols(g_parasitesLCA_l="family",g_parasitesLCA="3") %>% 
  select(qtax_genus,interaction_g_parasites,g_parasitesLCA_v,g_parasitesLCA_l,
         g_parasitesLCA) -> g_familydata #32

globi_parasites %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            g_parasitesLCA_v=unique(order)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder==1) %>%
  bind_cols(g_parasitesLCA_l="order",g_parasitesLCA="4") %>% 
  select(qtax_genus,interaction_g_parasites,g_parasitesLCA_v,g_parasitesLCA_l,
         g_parasitesLCA) -> g_orderdata #21

globi_parasites %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            g_parasitesLCA_v=unique(class)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass==1) %>%
  bind_cols(g_parasitesLCA_l="class",g_parasitesLCA="5") %>% 
  select(qtax_genus,interaction_g_parasites,g_parasitesLCA_v,g_parasitesLCA_l,
         g_parasitesLCA) -> g_classdata #102

globi_parasites %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            g_parasitesLCA_v=unique(phylum)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass!=1 & 
            nphylum==1) %>%
  bind_cols(g_parasitesLCA_l="phylum",g_parasitesLCA="6") %>% 
  select(qtax_genus,interaction_g_parasites,g_parasitesLCA_v,g_parasitesLCA_l,
         g_parasitesLCA) -> g_phylumdata #63

globi_parasites %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            g_parasitesLCA_v=unique(superkingdom)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass!=1 & 
            nphylum!=1 & nkingdom==1) %>%
  bind_cols(g_parasitesLCA_l="superkingdom",g_parasitesLCA="7") %>% 
  select(qtax_genus,interaction_g_parasites,g_parasitesLCA_v,g_parasitesLCA_l,
         g_parasitesLCA) -> g_superkingdomdata #172

globi_parasites %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom))) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass!=1 
          & nphylum!=1 & nkingdom>1) %>%
  bind_cols(g_parasitesLCA_l="Bac + Euk",
            g_parasitesLCA_v=">superkingdom",g_parasitesLCA="8") %>% 
  select(qtax_genus,interaction_g_parasites,g_parasitesLCA_v,g_parasitesLCA_l,
         g_parasitesLCA) -> g_nolcadata #57
  
g_parasites_LCA <- 
  as.data.frame(bind_rows(g_speciesdata, g_genusdata,g_familydata,g_orderdata,
                          g_classdata,g_phylumdata,g_superkingdomdata,
                          g_nolcadata)) #632

#--- Closest to vulpes level + value finden

globi_parasites %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes")) %>%
  filter(nvulpesspecies>=1) %>%
  bind_cols(g_parasitesCTV_l="species",g_parasitesCTV_v="Vulpes vulpes",
            g_parasitesCTV="1") %>%
  select(qtax_genus,g_parasitesCTV_v,g_parasitesCTV_l,g_parasitesCTV) -> 
  g_vspecies #6

globi_parasites %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes")) %>%
  filter(nvulpesspecies %in% c(NA,0) & nvulpesgenus>=1) %>%
  bind_cols(g_parasitesCTV_l="genus",g_parasitesCTV_v="Vulpes",
            g_parasitesCTV="2") %>%
  select(qtax_genus,g_parasitesCTV_v,g_parasitesCTV_l,g_parasitesCTV) -> 
  g_vgenus #11

globi_parasites %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam>=1) %>%
  bind_cols(g_parasitesCTV_l="family",g_parasitesCTV_v="Canidae",
            g_parasitesCTV="3") %>%
  select(qtax_genus,g_parasitesCTV_v,g_parasitesCTV_l,g_parasitesCTV) -> 
  g_vfam #27

globi_parasites %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpesorder>=1) %>%
  bind_cols(g_parasitesCTV_l="order",g_parasitesCTV_v="Carnivora",
            g_parasitesCTV="4") %>%
  select(qtax_genus,g_parasitesCTV_v,g_parasitesCTV_l,g_parasitesCTV) -> 
  g_vorder #60

globi_parasites %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpesorder %in% c(NA,0) & 
            nvulpesclass>=1) %>%
  bind_cols(g_parasitesCTV_l="class",g_parasitesCTV_v="Mammalia",
            g_parasitesCTV="5") %>%
  select(qtax_genus,g_parasitesCTV_v,g_parasitesCTV_l,g_parasitesCTV) -> 
  g_vclass #130

globi_parasites %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpesorder %in% c(NA,0) & 
            nvulpesclass %in% c(NA,0) & nvulpesphylum>=1) %>%
  bind_cols(g_parasitesCTV_l="phylum",g_parasitesCTV_v="Chordata",
            g_parasitesCTV="6") %>%
  select(qtax_genus,g_parasitesCTV_v,g_parasitesCTV_l,g_parasitesCTV) -> 
  g_vphylum #47

globi_parasites %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
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
  bind_cols(g_parasitesCTV_l="superkingdom",g_parasitesCTV_v="Eukaryota",
            g_parasitesCTV="7") %>%
  select(qtax_genus,g_parasitesCTV_v,g_parasitesCTV_l,g_parasitesCTV) -> 
  g_vkingdom #347

globi_parasites %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
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
  bind_cols(g_parasitesCTV_l="No level close to vulpes",g_parasitesCTV_v="No",
            g_parasitesCTV="8") %>%
  select(qtax_genus,g_parasitesCTV_v,g_parasitesCTV_l,g_parasitesCTV) -> 
  g_vnomatch #4

# no match bedeutet hier: Bacterium hat nur Interaktionen mit anderen Bacteria

g_parasites_CTV <- 
  data.frame(bind_rows(g_vspecies,g_vgenus,g_vfam,g_vorder,g_vclass,g_vphylum,
                       g_vkingdom,g_vnomatch)) #632

# LCA und closest to vulpes verbinden

g_parasites <- merge(g_parasites_LCA,g_parasites_CTV)


saveRDS(g_parasites,"g_parasites.Rds")
