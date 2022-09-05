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
  summarize(n_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            eatLCA_v=unique(species)) %>%
  filter (nspecies==1) %>%
  bind_cols(eatLCA_l="species") %>% 
  select(qtax_genus,n_eat,eatLCA_v,eatLCA_l) -> g_speciesdata2 #91

globi_eat %>% 
  group_by(qtax_genus) %>% 
  summarize(n_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            eatLCA_v=unique(genus)) %>%
  filter (nspecies !=1 & ngenus==1) %>%
  bind_cols(eatLCA_l="genus") %>% 
  select(qtax_genus,n_eat,eatLCA_v,eatLCA_l) -> g_genusdata2 #7

globi_eat %>% 
  group_by(qtax_genus) %>% 
  summarize(n_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            eatLCA_v=unique(family)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam==1) %>%
  bind_cols(eatLCA_l="family") %>% 
  select(qtax_genus,n_eat,eatLCA_v,eatLCA_l) -> g_familydata2 #4

globi_eat %>% 
  group_by(qtax_genus) %>% 
  summarize(n_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            eatLCA_v=unique(order)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder==1) %>%
  bind_cols(eatLCA_l="order") %>% 
  select(qtax_genus,n_eat,eatLCA_v,eatLCA_l) -> g_orderdata2 #11

globi_eat %>% 
  group_by(qtax_genus) %>% 
  summarize(n_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            eatLCA_v=unique(class)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass==1) %>%
  bind_cols(eatLCA_l="class") %>% 
  select(qtax_genus,n_eat,eatLCA_v,eatLCA_l) -> g_classdata2 #15

globi_eat %>% 
  group_by(qtax_genus) %>% 
  summarize(n_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            eatLCA_v=unique(phylum)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass!=1 & 
            nphylum==1) %>%
  bind_cols(eatLCA_l="phylum") %>% 
  select(qtax_genus,n_eat,eatLCA_v,eatLCA_l) -> g_phylumdata2 #20

globi_eat %>% 
  group_by(qtax_genus) %>% 
  summarize(n_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom)),
            eatLCA_v=unique(superkingdom)) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass!=1 & 
            nphylum!=1 & nkingdom==1) %>%
  bind_cols(eatLCA_l="superkingdom") %>% 
  select(qtax_genus,n_eat,eatLCA_v,eatLCA_l) -> 
  g_superkingdomdata2 #249

globi_eat %>% 
  group_by(qtax_genus) %>% 
  summarize(n_eat=n(),
            nspecies=length(unique(species)),
            ngenus=length(unique(genus)),
            nfam=length(unique(family)),
            norder=length(unique(order)),
            nclass=length(unique(class)),
            nphylum=length(unique(phylum)),
            nkingdom=length(unique(superkingdom))) %>%
  filter (nspecies !=1 & ngenus!=1 & nfam !=1 & norder!=1 & nclass!=1 & 
            nphylum!=1 & nkingdom>1) %>%
  bind_cols(eatLCA_l="Bac + Euk",eatLCA_v=">superkingdom") %>% 
  select(qtax_genus,n_eat,eatLCA_v,eatLCA_l) -> g_nolcadata2 #27

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
  bind_cols(eatCTV_l="species",eatCTV_v="Vulpes vulpes") %>%
  select(qtax_genus,eatCTV_v,eatCTV_l) -> g_vspecies2 #0

globi_eat %>% 
  group_by(qtax_genus) %>% 
  summarize(n_eat=n(),
            nvulpesspecies = sum(species=="Vulpes vulpes"),
            nvulpesgenus = sum(genus=="Vulpes")) %>%
  filter(nvulpesspecies %in% c(NA,0) & nvulpesgenus>=1) %>%
  bind_cols(eatCTV_l="genus",eatCTV_v="Vulpes") %>%
  select(qtax_genus,eatCTV_v,eatCTV_l) -> g_vgenus2 #11


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
  bind_cols(eatCTV_l="family",eatCTV_v="Canidae") %>%
  select(qtax_genus,eatCTV_v,eatCTV_l) -> g_vfam2 #5

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
  bind_cols(eatCTV_l="order",eatCTV_v="Carnivora") %>%
  select(qtax_genus,eatCTV_v,eatCTV_l) -> g_vorder2 #36

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
  bind_cols(eatCTV_l="class",eatCTV_v="Mammalia") %>%
  select(qtax_genus,eatCTV_v,eatCTV_l) -> g_vclass2 #49

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
  bind_cols(eatCTV_l="phylum",eatCTV_v="Chordata") %>%
  select(qtax_genus,eatCTV_v,eatCTV_l) -> g_vphylum2 #111

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
  bind_cols(eatCTV_l="superkingdom",eatCTV_v="Eukaryota") %>%
  select(qtax_genus,eatCTV_v,eatCTV_l) -> g_vkingdom2 #208

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
  bind_cols(eatCTV_l="No level ctv",eatCTV_v="No") %>%
  select(qtax_genus,eatCTV_v,eatCTV_l) -> g_vnomatch2 #4

# no match bedeutet hier: Bacterium hat nur Interaktionen mit anderen Bacteria

g_eat_CTV <- 
  data.frame(bind_rows(g_vspecies2,g_vgenus2,g_vfam2,g_vorder2,g_vclass2,
                       g_vphylum2,g_vkingdom2,g_vnomatch2)) #424

# LCA und closest to vulpes verbinden
g_eat <- merge(g_eat_LCA,g_eat_CTV)

saveRDS(g_eat,"g_eat.Rds")
