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
  summarize(n_parasites=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            parasitesLCA_v=unique(species)) %>%
  filter (nspecies==1) %>%
  bind_cols(parasitesLCA_l="species") %>% 
  select(qtax_genus,n_parasites,parasitesLCA_v,parasitesLCA_l) -> 
  g_speciesdata #161

globi_parasites %>% 
  group_by(qtax_genus) %>% 
  summarize(n_parasites=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            parasitesLCA_v=unique(genus)) %>%
  filter (nspecies !=1 & ngenus==1) %>%
  bind_cols(parasitesLCA_l="genus") %>% 
  select(qtax_genus,n_parasites,parasitesLCA_v,parasitesLCA_l) -> 
  g_genusdata #24
              
globi_parasites %>% 
  group_by(qtax_genus) %>% 
  summarize(n_parasites=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            parasitesLCA_v=unique(family)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam==1) %>%
  bind_cols(parasitesLCA_l="family") %>% 
  select(qtax_genus,n_parasites,parasitesLCA_v,parasitesLCA_l) -> 
  g_familydata #32

globi_parasites %>% 
  group_by(qtax_genus) %>% 
  summarize(n_parasites=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            parasitesLCA_v=unique(order)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder==1) %>%
  bind_cols(parasitesLCA_l="order") %>% 
  select(qtax_genus,n_parasites,parasitesLCA_v,parasitesLCA_l) -> 
  g_orderdata #21

globi_parasites %>% 
  group_by(qtax_genus) %>% 
  summarize(n_parasites=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            parasitesLCA_v=unique(class)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass==1) %>%
  bind_cols(parasitesLCA_l="class") %>% 
  select(qtax_genus,n_parasites,parasitesLCA_v,parasitesLCA_l) -> 
  g_classdata #102

globi_parasites %>% 
  group_by(qtax_genus) %>% 
  summarize(n_parasites=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            parasitesLCA_v=unique(phylum)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass!=1 & 
            nphylum==1) %>%
  bind_cols(parasitesLCA_l="phylum") %>% 
  select(qtax_genus,n_parasites,parasitesLCA_v,parasitesLCA_l) -> 
  g_phylumdata #63

globi_parasites %>% 
  group_by(qtax_genus) %>% 
  summarize(n_parasites=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            parasitesLCA_v=unique(superkingdom)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass!=1 & 
            nphylum!=1 & nkingdom==1) %>%
  bind_cols(parasitesLCA_l="superkingdom") %>% 
  select(qtax_genus,n_parasites,parasitesLCA_v,parasitesLCA_l) -> 
  g_superkingdomdata #172

globi_parasites %>% 
  group_by(qtax_genus) %>% 
  summarize(n_parasites=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom))) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass!=1 
          & nphylum!=1 & nkingdom>1) %>%
  bind_cols(parasitesLCA_l="Bac + Euk",
            parasitesLCA_v=">superkingdom") %>% 
  select(qtax_genus,n_parasites,parasitesLCA_v,parasitesLCA_l) -> 
  g_nolcadata #57
  
g_parasites_LCA <- 
  as.data.frame(bind_rows(g_speciesdata, g_genusdata,g_familydata,g_orderdata,
                          g_classdata,g_phylumdata,g_superkingdomdata,
                          g_nolcadata)) #632

#--- Closest to vulpes level + value finden

globi_parasites %>% 
  group_by(qtax_genus) %>% 
  summarize(n_parasites=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes")) %>%
  filter(nvulpesspecies>=1) %>%
  bind_cols(parasitesCTV_l="species",parasitesCTV_v="Vulpes vulpes") %>%
  select(qtax_genus,parasitesCTV_v,parasitesCTV_l) ->g_vspecies #6

globi_parasites %>% 
  group_by(qtax_genus) %>% 
  summarize(n_parasites=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes")) %>%
  filter(nvulpesspecies %in% c(NA,0) & nvulpesgenus>=1) %>%
  bind_cols(parasitesCTV_l="genus",parasitesCTV_v="Vulpes") %>%
  select(qtax_genus,parasitesCTV_v,parasitesCTV_l) -> g_vgenus #11


globi_parasites %>% 
  group_by(qtax_genus) %>% 
  summarize(n_parasites=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam>=1) %>%
  bind_cols(parasitesCTV_l="family",parasitesCTV_v="Canidae") %>%
  select(qtax_genus,parasitesCTV_v,parasitesCTV_l) -> g_vfam #27

globi_parasites %>% 
  group_by(qtax_genus) %>% 
  summarize(n_parasites=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpesorder>=1) %>%
  bind_cols(parasitesCTV_l="order",parasitesCTV_v="Carnivora") %>%
  select(qtax_genus,parasitesCTV_v,parasitesCTV_l) -> g_vorder #60

globi_parasites %>% 
  group_by(qtax_genus) %>% 
  summarize(n_parasites=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpesorder %in% c(NA,0) & 
            nvulpesclass>=1) %>%
  bind_cols(parasitesCTV_l="class",parasitesCTV_v="Mammalia") %>%
  select(qtax_genus,parasitesCTV_v,parasitesCTV_l) -> g_vclass #130

globi_parasites %>% 
  group_by(qtax_genus) %>% 
  summarize(n_parasites=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes"),
            nvulpesfam = sum(family=="Canidae"),
            nvulpesorder = sum(order=="Carnivora"),
            nvulpesclass = sum(class=="Mammalia"),
            nvulpesphylum = sum(phylum=="Chordata")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpesorder %in% c(NA,0) & 
            nvulpesclass %in% c(NA,0) & nvulpesphylum>=1) %>%
  bind_cols(parasitesCTV_l="phylum",parasitesCTV_v="Chordata") %>%
  select(qtax_genus,parasitesCTV_v,parasitesCTV_l) -> g_vphylum #47

globi_parasites %>% 
  group_by(qtax_genus) %>% 
  summarize(n_parasites=n(),
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
  bind_cols(parasitesCTV_l="superkingdom",parasitesCTV_v="Eukaryota") %>%
  select(qtax_genus,parasitesCTV_v,parasitesCTV_l) -> g_vkingdom #347

globi_parasites %>% 
  group_by(qtax_genus) %>% 
  summarize(n_parasites=n(),
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
  bind_cols(parasitesCTV_l="No level close to vulpes",parasitesCTV_v="No") %>%
  select(qtax_genus,parasitesCTV_v,parasitesCTV_l) -> g_vnomatch #4

# no match bedeutet hier: Bacterium hat nur Interaktionen mit anderen Bacteria

g_parasites_CTV <- 
  data.frame(bind_rows(g_vspecies,g_vgenus,g_vfam,g_vorder,g_vclass,g_vphylum,
                       g_vkingdom,g_vnomatch)) #632

# LCA und closest to vulpes verbinden

g_parasites <- merge(g_parasites_LCA,g_parasites_CTV)


saveRDS(g_parasites,"g_parasites.Rds")

                                                                              
