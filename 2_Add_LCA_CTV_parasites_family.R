newglobi <- FALSE
if (newglobi) {
  source("1_globidata_parasites_family")
} else {
  globidata_parasites_f <- readRDS("globi_parasites_family.Rds")
}

library(tidyverse)
library(taxonomizr)

#Taxon Path von Vulpes Vulpes herausfinden
desiredtaxa = c("superkingdom","kingdom","superphylum","phylum",
                "subphylum", "superclass","class","subclass", "superorder",
                "order","suborder","family","genus","species")

vulpestaxonpath <- data.frame(getTaxonomy(getId("Vulpes vulpes")))

#nur Bakterien, Viren und Eukaryoten auswählen
globidata_parasites_f %>% 
  filter(globi_superkingdom %in% c("Bacteria","Eukaryota","Viruses")) -> 
  globi_parasites_family2 #334317

# Doppelte entfernen
globi_parasites_family <- unique(globi_parasites_family2) #22238

globi_parasites_family %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasites=n()) ->x #419

#species LCA
globi_parasites_family %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasites=n(),
            nspecies=length(unique(globi_species)),
            ngenus=length(unique(globi_genus)),
            nfam=length(unique(globi_family)),
            nsuborder=length(unique(globi_suborder)),
            norder=length(unique(globi_order)),
            nsuperorder=length(unique(globi_superorder)),
            nsubclass=length(unique(globi_subclass)),
            nclass=length(unique(globi_class)),
            nsuperclass=length(unique(globi_superclass)),
            nsubphylum=length(unique(globi_subphylum)),
            nphylum=length(unique(globi_phylum)),
            nsuperphylum=length(unique(globi_superphylum)),
            nkingdom=length(unique(globi_kingdom)),
            nsuperkingdom=length(unique(globi_superkingdom)),
            f_parasitesLCA_v=unique(globi_species)) %>%
  filter (nspecies==1 & ngenus==1 & nfam==1 & nsuborder==1 & norder==1 &
            nsuperorder ==1 & nsubclass==1 & nclass==1 & nsuperclass==1 & 
            nsubphylum==1 & nphylum==1 & nkingdom==1 & nsuperkingdom==1 & 
            !is.na(f_parasitesLCA_v)) %>%
  bind_cols(f_parasitesLCA_l="species",f_parasitesLCA="1") %>%
  select(qtax_family,
         interaction_f_parasites,
         f_parasitesLCA_v,
         f_parasitesLCA_l,
         f_parasitesLCA) -> f_speciesdata #46

#genus LCA
globi_parasites_family %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasites=n(),
            nspecies=length(unique(globi_species)),
            ngenus=length(unique(globi_genus)),
            nfam=length(unique(globi_family)),
            nsuborder=length(unique(globi_suborder)),
            norder=length(unique(globi_order)),
            nsuperorder=length(unique(globi_superorder)),
            nsubclass=length(unique(globi_subclass)),
            nclass=length(unique(globi_class)),
            nsuperclass=length(unique(globi_superclass)),
            nsubphylum=length(unique(globi_subphylum)),
            nphylum=length(unique(globi_phylum)),
            nsuperphylum=length(unique(globi_superphylum)),
            nkingdom=length(unique(globi_kingdom)),
            nsuperkingdom=length(unique(globi_superkingdom)),
            f_parasitesLCA_v=unique(globi_genus)) %>%
  filter (nspecies>1 & ngenus==1 & nfam==1 & nsuborder==1 & norder==1 &
            nsuperorder ==1 & nsubclass==1 & nclass==1 & nsuperclass==1 & 
            nsubphylum==1 & nphylum==1 & nkingdom==1 & nsuperkingdom==1 & 
            !is.na(f_parasitesLCA_v)) %>%
  bind_cols(f_parasitesLCA_l="genus",f_parasitesLCA="2") %>%
  select(qtax_family,
         interaction_f_parasites,
         f_parasitesLCA_v,
         f_parasitesLCA_l,
         f_parasitesLCA) -> f_genusdata #7

#family LCA
globi_parasites_family  %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasites=n(),
            nspecies=length(unique(globi_species)),
            ngenus=length(unique(globi_genus)),
            nfam=length(unique(globi_family)),
            nsuborder=length(unique(globi_suborder)),
            norder=length(unique(globi_order)),
            nsuperorder=length(unique(globi_superorder)),
            nsubclass=length(unique(globi_subclass)),
            nclass=length(unique(globi_class)),
            nsuperclass=length(unique(globi_superclass)),
            nsubphylum=length(unique(globi_subphylum)),
            nphylum=length(unique(globi_phylum)),
            nsuperphylum=length(unique(globi_superphylum)),
            nkingdom=length(unique(globi_kingdom)),
            nsuperkingdom=length(unique(globi_superkingdom)),
            f_parasitesLCA_v=unique(globi_family)) %>%
  filter (nspecies>1 & ngenus>1 & nfam==1 & nsuborder==1 & norder==1 &
            nsuperorder ==1 & nsubclass==1 & nclass==1 & nsuperclass==1 & 
            nsubphylum==1 & nphylum==1 & nkingdom==1 & nsuperkingdom==1 & 
            !is.na(f_parasitesLCA_v))%>%
  bind_cols(f_parasitesLCA_l="family",f_parasitesLCA="3") %>%
  select(qtax_family,
         interaction_f_parasites,
         f_parasitesLCA_v,
         f_parasitesLCA_l,
         f_parasitesLCA) -> f_familydata #13

#suborder LCA 
globi_parasites_family %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasites=n(),
            nspecies=length(unique(globi_species)),
            ngenus=length(unique(globi_genus)),
            nfam=length(unique(globi_family)),
            nsuborder=length(unique(globi_suborder)),
            norder=length(unique(globi_order)),
            nsuperorder=length(unique(globi_superorder)),
            nsubclass=length(unique(globi_subclass)),
            nclass=length(unique(globi_class)),
            nsuperclass=length(unique(globi_superclass)),
            nsubphylum=length(unique(globi_subphylum)),
            nphylum=length(unique(globi_phylum)),
            nsuperphylum=length(unique(globi_superphylum)),
            nkingdom=length(unique(globi_kingdom)),
            nsuperkingdom=length(unique(globi_superkingdom)),
            f_parasitesLCA_v=unique(globi_suborder)) %>%
  filter (nspecies>1 & ngenus>1 & nfam>1 & nsuborder==1 & norder==1 &
            nsuperorder ==1 & nsubclass==1 & nclass==1 & nsuperclass==1 & 
            nsubphylum==1 & nphylum==1 & nkingdom==1 & nsuperkingdom==1 & 
            !is.na(f_parasitesLCA_v)) %>%
  bind_cols(f_parasitesLCA_l="suborder",f_parasitesLCA="4") %>%
  select(qtax_family,
         interaction_f_parasites,
         f_parasitesLCA_v,
         f_parasitesLCA_l,
         f_parasitesLCA) -> f_suborderdata #2

#order LCA
globi_parasites_family %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasites=n(),
            nspecies=length(unique(globi_species)),
            ngenus=length(unique(globi_genus)),
            nfam=length(unique(globi_family)),
            nsuborder=length(unique(globi_suborder)),
            norder=length(unique(globi_order)),
            nsuperorder=length(unique(globi_superorder)),
            nsubclass=length(unique(globi_subclass)),
            nclass=length(unique(globi_class)),
            nsuperclass=length(unique(globi_superclass)),
            nsubphylum=length(unique(globi_subphylum)),
            nphylum=length(unique(globi_phylum)),
            nsuperphylum=length(unique(globi_superphylum)),
            nkingdom=length(unique(globi_kingdom)),
            nsuperkingdom=length(unique(globi_superkingdom)),
            f_parasitesLCA_v=unique(globi_order))%>% 
  filter (nspecies>1 & ngenus>1 & nfam>1 & nsuborder>1 & norder==1 &
            nsuperorder ==1 & nsubclass==1 & nclass==1 & nsuperclass==1 & 
            nsubphylum==1 & nphylum==1 & nkingdom==1 & nsuperkingdom==1 & 
            !is.na(f_parasitesLCA_v)) %>%
  bind_cols(f_parasitesLCA_l="order",f_parasitesLCA="5") %>%
  select(qtax_family,
         interaction_f_parasites,
         f_parasitesLCA_v,
         f_parasitesLCA_l,
         f_parasitesLCA) -> f_orderdata #1

#superorder LCA
globi_parasites_family %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasites=n(),
            nspecies=length(unique(globi_species)),
            ngenus=length(unique(globi_genus)),
            nfam=length(unique(globi_family)),
            nsuborder=length(unique(globi_suborder)),
            norder=length(unique(globi_order)),
            nsuperorder=length(unique(globi_superorder)),
            nsubclass=length(unique(globi_subclass)),
            nclass=length(unique(globi_class)),
            nsuperclass=length(unique(globi_superclass)),
            nsubphylum=length(unique(globi_subphylum)),
            nphylum=length(unique(globi_phylum)),
            nsuperphylum=length(unique(globi_superphylum)),
            nkingdom=length(unique(globi_kingdom)),
            nsuperkingdom=length(unique(globi_superkingdom)),
            f_parasitesLCA_v=unique(globi_superorder)) %>%
  filter (nspecies!=1 & ngenus!=1 & nfam!=1 & nsuborder!=1 & norder!=1 &
            nsuperorder==1 & nsubclass==1 & nclass==1 & nsuperclass==1 & 
            nsubphylum==1 & nphylum==1 & nkingdom==1 & nsuperkingdom==1 &
            !is.na(f_parasitesLCA_v)) %>%
  bind_cols(f_parasitesLCA_l="superorder",f_parasitesLCA="6") %>%
  select(qtax_family,
         interaction_f_parasites,
         f_parasitesLCA_v,
         f_parasitesLCA_l,
         f_parasitesLCA) -> f_superorderdata #3

#subclass LCA
globi_parasites_family %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasites=n(),
            nspecies=length(unique(globi_species)),
            ngenus=length(unique(globi_genus)),
            nfam=length(unique(globi_family)),
            nsuborder=length(unique(globi_suborder)),
            norder=length(unique(globi_order)),
            nsuperorder=length(unique(globi_superorder)),
            nsubclass=length(unique(globi_subclass)),
            nclass=length(unique(globi_class)),
            nsuperclass=length(unique(globi_superclass)),
            nsubphylum=length(unique(globi_subphylum)),
            nphylum=length(unique(globi_phylum)),
            nsuperphylum=length(unique(globi_superphylum)),
            nkingdom=length(unique(globi_kingdom)),
            nsuperkingdom=length(unique(globi_superkingdom)),
            f_parasitesLCA_v=unique(globi_subclass)) %>%
  filter (nspecies!=1 & ngenus!=1 & nfam!=1 & nsuborder!=1 & norder!=1 &
            nsuperorder!=1 & nsubclass==1 & nclass==1 & 
            nsuperclass==1 & nsubphylum==1 & nphylum==1 & nkingdom==1 & 
            nsuperkingdom==1 & !is.na(f_parasitesLCA_v)) %>%
  bind_cols(f_parasitesLCA_l="subclass",f_parasitesLCA="7") %>%
  select(qtax_family,
         interaction_f_parasites,
         f_parasitesLCA_v,
         f_parasitesLCA_l,
         f_parasitesLCA) -> f_subclassdata #5

#class LCA
globi_parasites_family %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasites=n(),
            nspecies=length(unique(globi_species)),
            ngenus=length(unique(globi_genus)),
            nfam=length(unique(globi_family)),
            nsuborder=length(unique(globi_suborder)),
            norder=length(unique(globi_order)),
            nsuperorder=length(unique(globi_superorder)),
            nsubclass=length(unique(globi_subclass)),
            nclass=length(unique(globi_class)),
            nsuperclass=length(unique(globi_superclass)),
            nsubphylum=length(unique(globi_subphylum)),
            nphylum=length(unique(globi_phylum)),
            nsuperphylum=length(unique(globi_superphylum)),
            nkingdom=length(unique(globi_kingdom)),
            nsuperkingdom=length(unique(globi_superkingdom)),
            f_parasitesLCA_v=unique(globi_class)) %>%
  filter (nspecies!=1 & ngenus!=1 & nfam!=1 & nsuborder!=1 & norder!=1 &
            nsuperorder!=1 & nsubclass==1 & nclass==1 & 
            nsuperclass==1 & nsubphylum==1 & nphylum==1 & nkingdom==1 & 
            nsuperkingdom==1 & is.na(f_parasitesLCA_v)) %>%
  bind_cols(f_parasitesLCA_l="class",f_parasitesLCA="7") %>%
  select(qtax_family,
         interaction_f_parasites,
         f_parasitesLCA_v,
         f_parasitesLCA_l,
         f_parasitesLCA) -> f_classdata #15

#superclass LCA
globi_parasites_family %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasites=n(),
            nspecies=length(unique(globi_species)),
            ngenus=length(unique(globi_genus)),
            nfam=length(unique(globi_family)),
            nsuborder=length(unique(globi_suborder)),
            norder=length(unique(globi_order)),
            nsuperorder=length(unique(globi_superorder)),
            nsubclass=length(unique(globi_subclass)),
            nclass=length(unique(globi_class)),
            nsuperclass=length(unique(globi_superclass)),
            nsubphylum=length(unique(globi_subphylum)),
            nphylum=length(unique(globi_phylum)),
            nsuperphylum=length(unique(globi_superphylum)),
            nkingdom=length(unique(globi_kingdom)),
            nsuperkingdom=length(unique(globi_superkingdom)),
            f_parasitesLCA_v=unique(globi_superclass)) %>%
  filter (nspecies!=1 & ngenus!=1 & nfam!=1 & nsuborder!=1 & norder!=1 &
            nsuperorder!=1 & nsubclass!=1 & nclass!=1 & nsuperclass==1 & 
            nsubphylum==1 & nphylum==1 & nkingdom==1 & nsuperkingdom==1 & 
            !is.na(f_parasitesLCA_v)) %>%
  bind_cols(f_parasitesLCA_l="superclass",f_parasitesLCA="9") %>%
  select(qtax_family,
         interaction_f_parasites,
         f_parasitesLCA_v,
         f_parasitesLCA_l,
         f_parasitesLCA) -> f_superclassdata #4

#subphylum LCA
globi_parasites_family %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasites=n(),
            nspecies=length(unique(globi_species)),
            ngenus=length(unique(globi_genus)),
            nfam=length(unique(globi_family)),
            nsuborder=length(unique(globi_suborder)),
            norder=length(unique(globi_order)),
            nsuperorder=length(unique(globi_superorder)),
            nsubclass=length(unique(globi_subclass)),
            nclass=length(unique(globi_class)),
            nsuperclass=length(unique(globi_superclass)),
            nsubphylum=length(unique(globi_subphylum)),
            nphylum=length(unique(globi_phylum)),
            nsuperphylum=length(unique(globi_superphylum)),
            nkingdom=length(unique(globi_kingdom)),
            nsuperkingdom=length(unique(globi_superkingdom)),
            f_parasitesLCA_v=unique(globi_subphylum)) %>%
  filter (nspecies!=1 & ngenus!=1 & nfam!=1 & nsuborder!=1 & norder!=1 &
            nsuperorder!=1 & nsubclass!=1 & nclass!=1, nsuperclass!=1 & 
            nsubphylum==1 & nphylum==1 & nkingdom==1 & nsuperkingdom==1 &
            !is.na(f_parasitesLCA_v)) %>%
  bind_cols(f_parasitesLCA_l="subphylum",f_parasitesLCA="10") %>%
  select(qtax_family,
         interaction_f_parasites,
         f_parasitesLCA_v,
         f_parasitesLCA_l,
         f_parasitesLCA) -> f_subphylumdata #11

#phylum LCA
globi_parasites_family %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasites=n(),
            nspecies=length(unique(globi_species)),
            ngenus=length(unique(globi_genus)),
            nfam=length(unique(globi_family)),
            nsuborder=length(unique(globi_suborder)),
            norder=length(unique(globi_order)),
            nsuperorder=length(unique(globi_superorder)),
            nsubclass=length(unique(globi_subclass)),
            nclass=length(unique(globi_class)),
            nsuperclass=length(unique(globi_superclass)),
            nsubphylum=length(unique(globi_subphylum)),
            nphylum=length(unique(globi_phylum)),
            nsuperphylum=length(unique(globi_superphylum)),
            nkingdom=length(unique(globi_kingdom)),
            nsuperkingdom=length(unique(globi_superkingdom)),
            f_parasitesLCA_v=unique(globi_phylum)) %>%
  filter (nspecies!=1 & ngenus!=1 & nfam!=1 & nsuborder!=1 & norder!=1 &
            nsuperorder!=1 & nsubclass!=1 & nclass!=1 & nsuperclass!=1 & 
            nsubphylum!=1 & nphylum==1 & nkingdom==1 & nsuperkingdom==1 &
            !is.na(f_parasitesLCA_v)) %>%
  bind_cols(f_parasitesLCA_l="phylum",f_parasitesLCA="11") %>%
  select(qtax_family,
         interaction_f_parasites,
         f_parasitesLCA_v,
         f_parasitesLCA_l,
         f_parasitesLCA) -> f_phylumdata #2

#kingdom LCA
globi_parasites_family %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasites=n(),
            nspecies=length(unique(globi_species)),
            ngenus=length(unique(globi_genus)),
            nfam=length(unique(globi_family)),
            nsuborder=length(unique(globi_suborder)),
            norder=length(unique(globi_order)),
            nsuperorder=length(unique(globi_superorder)),
            nsubclass=length(unique(globi_subclass)),
            nclass=length(unique(globi_class)),
            nsuperclass=length(unique(globi_superclass)),
            nsubphylum=length(unique(globi_subphylum)),
            nphylum=length(unique(globi_phylum)),
            nsuperphylum=length(unique(globi_superphylum)),
            nkingdom=length(unique(globi_kingdom)),
            nsuperkingdom=length(unique(globi_superkingdom)),
            f_parasitesLCA_v=unique(globi_kingdom)) %>%
  filter (nspecies!=1 & ngenus!=1 & nfam!=1 & nsuborder!=1 & norder!=1 &
            nsuperorder!=1 & nsubclass!=1 & nclass!=1 & nsuperclass!=1 & 
            nsubphylum!=1 & nphylum!=1 & nkingdom==1 & nsuperkingdom==1 &
            !is.na(f_parasitesLCA_v)) %>%
  bind_cols(f_parasitesLCA_l="kingdom",f_parasitesLCA="12") %>%
  select(qtax_family,
         interaction_f_parasites,
         f_parasitesLCA_v,
         f_parasitesLCA_l,
         f_parasitesLCA) -> f_kingdomdata #63

#superkingdom LCA
globi_parasites_family %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasites=n(),
            nspecies=length(unique(globi_species)),
            ngenus=length(unique(globi_genus)),
            nfam=length(unique(globi_family)),
            nsuborder=length(unique(globi_suborder)),
            norder=length(unique(globi_order)),
            nsuperorder=length(unique(globi_superorder)),
            nsubclass=length(unique(globi_subclass)),
            nclass=length(unique(globi_class)),
            nsuperclass=length(unique(globi_superclass)),
            nsubphylum=length(unique(globi_subphylum)),
            nphylum=length(unique(globi_phylum)),
            nsuperphylum=length(unique(globi_superphylum)),
            nkingdom=length(unique(globi_kingdom)),
            nsuperkingdom=length(unique(globi_superkingdom)),
            f_parasitesLCA_v=unique(globi_superkingdom)) %>%
  filter (nspecies!=1 & ngenus!=1 & nfam!=1 & nsuborder!=1 & norder!=1 &
            nsuperorder!=1 & nsubclass!=1 & nclass!=1, nsuperclass!=1 & 
            nsubphylum!=1,nphylum!=1,nkingdom!=1,nsuperkingdom==1 & 
            !is.na(f_parasitesLCA_v)) %>%
  bind_cols(f_parasitesLCA_l="superkingdom",f_parasitesLCA="13") %>%
  select(qtax_family,
         interaction_f_parasites,
         f_parasitesLCA_v,
         f_parasitesLCA_l,
         f_parasitesLCA) -> f_superkingdomdata #58

#no LCA
globi_parasites_family %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasites=n(),
            nspecies=length(unique(globi_species)),
            ngenus=length(unique(globi_genus)),
            nfam=length(unique(globi_family)),
            nsuborder=length(unique(globi_suborder)),
            norder=length(unique(globi_order)),
            nsuperorder=length(unique(globi_superorder)),
            nsubclass=length(unique(globi_subclass)),
            nclass=length(unique(globi_class)),
            nsuperclass=length(unique(globi_superclass)),
            nsubphylum=length(unique(globi_subphylum)),
            nphylum=length(unique(globi_phylum)),
            nsuperphylum=length(unique(globi_superphylum)),
            nkingdom=length(unique(globi_kingdom)),
            nsuperkingdom=length(unique(globi_superkingdom)),
            f_parasitesLCA_v=unique(globi_superkingdom)) %>%
  filter (nspecies!=1 & ngenus!=1 & nfam!=1 & nsuborder!=1 & norder!=1 &
            nsuperorder!=1 & nsubclass!=1 & nclass!=1, nsuperclass!=1 & 
            nsubphylum!=1,nphylum!=1,nkingdom!=1,nsuperkingdom>1 &
            !is.na(f_parasitesLCA_v)) %>%
  bind_cols(f_parasitesLCA_l=">Superkingdom",f_parasitesLCA="14") %>%
  select(qtax_family,
         interaction_f_parasites,
         f_parasitesLCA_v,
         f_parasitesLCA_l,
         f_parasitesLCA) -> f_nodata #8

f_parasites_LCA <- 
  unique(as.data.frame(bind_rows(f_speciesdata,f_genusdata,f_familydata,
                                 f_suborderdata,f_orderdata,f_superorderdata,
                                 f_subclassdata,f_classdata,f_superclassdata, 
                                 f_subphylumdata,f_phylumdata,f_kingdomdata,
                                 f_superkingdomdata,f_nodata))) #238

#--- Closest to vulpes level + value finden

globi_parasites_family %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasites=n(),
            nvulpesspecies = sum(globi_species=="Vulpes vulpes")) %>%
  filter(nvulpesspecies>=1) %>%
  bind_cols(f_parasitesCTV_l="species",f_parasitesCTV_v="Vulpes vulpes",
            f_parasitesCTV="1") %>%
  select(qtax_family,f_parasitesCTV_v,f_parasitesCTV_l,f_parasitesCTV) ->
  f_vspecies #0

globi_parasites_family %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasite=n(),
            nvulpesspecies = sum(globi_species=="Vulpes vulpes"),
            nvulpesgenus = sum(globi_genus=="Vulpes"),
            nvulpesfam = sum(globi_family=="Canidae"),
            nvulpessuborder = sum(globi_suborder=="Caniformia"),
            nvulpesorder = sum(globi_order=="Carnivora"),
            nvulpessuperorder = sum(globi_superorder=="Laurasiatheria"),
            nvulpesclass = sum(globi_class=="Mammalia"),
            nvulpessuperclass = sum(globi_superclass=="Sarcopterygii"),
            nvulpessubphylum = sum(globi_subphylum=="Craniata"),
            nvulpesphylum = sum(globi_phylum=="Chordata"),
            nvulpeskingdom = sum(globi_kingdom=="Metazoa"),
            nvulpessuperkingdom = sum(globi_superkingdom=="Eukaryota")) %>%
  filter(nvulpesspecies %in% c(NA,0) & nvulpesgenus>=1) %>%
  bind_cols(f_parasitesCTV_l="genus",f_parasitesCTV_v="Vulpes",
            f_parasitesCTV="2") %>%
  select(qtax_family,f_parasitesCTV_v,f_parasitesCTV_l,f_parasitesCTV) -> 
  f_vgenus #3

globi_parasites_family %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasites=n(),
            nvulpesspecies = sum(globi_species=="Vulpes vulpes"),
            nvulpesgenus = sum(globi_genus=="Vulpes"),
            nvulpesfam = sum(globi_family=="Canidae"),
            nvulpessuborder = sum(globi_suborder=="Caniformia"),
            nvulpesorder = sum(globi_order=="Carnivora"),
            nvulpessuperorder = sum(globi_superorder=="Laurasiatheria"),
            nvulpesclass = sum(globi_class=="Mammalia"),
            nvulpessuperclass = sum(globi_superclass=="Sarcopterygii"),
            nvulpessubphylum = sum(globi_subphylum=="Craniata"),
            nvulpesphylum = sum(globi_phylum=="Chordata"),
            nvulpeskingdom = sum(globi_kingdom=="Metazoa"),
            nvulpessuperkingdom = sum(globi_superkingdom=="Eukaryota")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam>=1) %>%
  bind_cols(f_parasitesCTV_l="family",f_parasitesCTV_v="Canidae",
            f_parasitesCTV="3") %>%
  select(qtax_family,f_parasitesCTV_v,f_parasitesCTV_l,f_parasitesCTV) -> 
  f_vfam #31

globi_parasites_family %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasites=n(),
            nvulpesspecies = sum(globi_species=="Vulpes vulpes"),
            nvulpesgenus = sum(globi_genus=="Vulpes"),
            nvulpesfam = sum(globi_family=="Canidae"),
            nvulpessuborder = sum(globi_suborder=="Caniformia"),
            nvulpesorder = sum(globi_order=="Carnivora"),
            nvulpessuperorder = sum(globi_superorder=="Laurasiatheria"),
            nvulpesclass = sum(globi_class=="Mammalia"),
            nvulpessuperclass = sum(globi_superclass=="Sarcopterygii"),
            nvulpessubphylum = sum(globi_subphylum=="Craniata"),
            nvulpesphylum = sum(globi_phylum=="Chordata"),
            nvulpeskingdom = sum(globi_kingdom=="Metazoa"),
            nvulpessuperkingdom = sum(globi_superkingdom=="Eukaryota")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpessuborder >=1) %>%
  bind_cols(f_parasitesCTV_l="order",f_parasitesCTV_v="Caniformia",
            f_parasitesCTV="5") %>%
  select(qtax_family,f_parasitesCTV_v,f_parasitesCTV_l,f_parasitesCTV) -> 
  f_vsuborder #1

globi_parasites_family %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasites=n(),
            nvulpesspecies = sum(globi_species=="Vulpes vulpes"),
            nvulpesgenus = sum(globi_genus=="Vulpes"),
            nvulpesfam = sum(globi_family=="Canidae"),
            nvulpessuborder = sum(globi_suborder=="Caniformia"),
            nvulpesorder = sum(globi_order=="Carnivora"),
            nvulpessuperorder = sum(globi_superorder=="Laurasiatheria"),
            nvulpesclass = sum(globi_class=="Mammalia"),
            nvulpessuperclass = sum(globi_superclass=="Sarcopterygii"),
            nvulpessubphylum = sum(globi_subphylum=="Craniata"),
            nvulpesphylum = sum(globi_phylum=="Chordata"),
            nvulpeskingdom = sum(globi_kingdom=="Metazoa"),
            nvulpessuperkingdom = sum(globi_superkingdom=="Eukaryota")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpessuborder %in% c(NA,0) & 
            nvulpesorder>=1) %>%
  bind_cols(f_parasitesCTV_l="order",f_parasitesCTV_v="Carnivora",
            f_parasitesCTV="5") %>%
  select(qtax_family,f_parasitesCTV_v,f_parasitesCTV_l,f_parasitesCTV) -> 
  f_vorder #36

globi_parasites_family %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasites=n(),
            nvulpesspecies = sum(globi_species=="Vulpes vulpes"),
            nvulpesgenus = sum(globi_genus=="Vulpes"),
            nvulpesfam = sum(globi_family=="Canidae"),
            nvulpessuborder = sum(globi_suborder=="Caniformia"),
            nvulpesorder = sum(globi_order=="Carnivora"),
            nvulpessuperorder = sum(globi_superorder=="Laurasiatheria"),
            nvulpesclass = sum(globi_class=="Mammalia"),
            nvulpessuperclass = sum(globi_superclass=="Sarcopterygii"),
            nvulpessubphylum = sum(globi_subphylum=="Craniata"),
            nvulpesphylum = sum(globi_phylum=="Chordata"),
            nvulpeskingdom = sum(globi_kingdom=="Metazoa"),
            nvulpessuperkingdom = sum(globi_superkingdom=="Eukaryota")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpessuborder %in% c(NA,0) & 
            nvulpesorder %in% c(NA,0) & nvulpessuperorder>=1) %>%
  bind_cols(f_parasitesCTV_l="superorder",f_parasitesCTV_v="Laurasiatheria",
            f_parasitesCTV="6") %>%
  select(qtax_family,f_parasitesCTV_v,f_parasitesCTV_l,f_parasitesCTV) -> 
  f_vsuperorder #16

globi_parasites_family %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasites=n(),
            nvulpesspecies = sum(globi_species=="Vulpes vulpes"),
            nvulpesgenus = sum(globi_genus=="Vulpes"),
            nvulpesfam = sum(globi_family=="Canidae"),
            nvulpessuborder = sum(globi_suborder=="Caniformia"),
            nvulpesorder = sum(globi_order=="Carnivora"),
            nvulpessuperorder = sum(globi_superorder=="Laurasiatheria"),
            nvulpesclass = sum(globi_class=="Mammalia"),
            nvulpessuperclass = sum(globi_superclass=="Sarcopterygii"),
            nvulpessubphylum = sum(globi_subphylum=="Craniata"),
            nvulpesphylum = sum(globi_phylum=="Chordata"),
            nvulpeskingdom = sum(globi_kingdom=="Metazoa"),
            nvulpessuperkingdom = sum(globi_superkingdom=="Eukaryota")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpessuborder %in% c(NA,0) & 
            nvulpesorder %in% c(NA,0) & nvulpessuperorder %in% c(NA,0) & 
            nvulpesclass>=1) %>%
  bind_cols(f_parasitesCTV_l="class",f_parasitesCTV_v="Mammalia",
            f_parasitesCTV="7") %>%
  select(qtax_family,f_parasitesCTV_v,f_parasitesCTV_l,f_parasitesCTV) -> 
  f_vclass #70

globi_parasites_family %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasites=n(),
            nvulpesspecies = sum(globi_species=="Vulpes vulpes"),
            nvulpesgenus = sum(globi_genus=="Vulpes"),
            nvulpesfam = sum(globi_family=="Canidae"),
            nvulpessuborder = sum(globi_suborder=="Caniformia"),
            nvulpesorder = sum(globi_order=="Carnivora"),
            nvulpessuperorder = sum(globi_superorder=="Laurasiatheria"),
            nvulpesclass = sum(globi_class=="Mammalia"),
            nvulpessuperclass = sum(globi_superclass=="Sarcopterygii"),
            nvulpessubphylum = sum(globi_subphylum=="Craniata"),
            nvulpesphylum = sum(globi_phylum=="Chordata"),
            nvulpeskingdom = sum(globi_kingdom=="Metazoa"),
            nvulpessuperkingdom = sum(globi_superkingdom=="Eukaryota")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpessuborder %in% c(NA,0) & 
            nvulpesorder %in% c(NA,0) & nvulpessuperorder %in% c(NA,0) & 
            nvulpesclass%in% c(NA,0) & nvulpessuperclass >=1) %>%
  bind_cols(f_parasitesCTV_l="superclass",f_parasitesCTV_v="Sarcopterygii",
            f_parasitesCTV="8") %>%
  select(qtax_family,f_parasitesCTV_v,f_parasitesCTV_l,f_parasitesCTV) -> 
  f_vsuperclass #8

globi_parasites_family %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasites=n(),
            nvulpesspecies = sum(globi_species=="Vulpes vulpes"),
            nvulpesgenus = sum(globi_genus=="Vulpes"),
            nvulpesfam = sum(globi_family=="Canidae"),
            nvulpessuborder = sum(globi_suborder=="Caniformia"),
            nvulpesorder = sum(globi_order=="Carnivora"),
            nvulpessuperorder = sum(globi_superorder=="Laurasiatheria"),
            nvulpesclass = sum(globi_class=="Mammalia"),
            nvulpessuperclass = sum(globi_superclass=="Sarcopterygii"),
            nvulpessubphylum = sum(globi_subphylum=="Craniata"),
            nvulpesphylum = sum(globi_phylum=="Chordata"),
            nvulpeskingdom = sum(globi_kingdom=="Metazoa"),
            nvulpessuperkingdom = sum(globi_superkingdom=="Eukaryota")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpessuborder %in% c(NA,0) & 
            nvulpesorder %in% c(NA,0) & nvulpessuperorder %in% c(NA,0) & 
            nvulpesclass%in% c(NA,0) & nvulpessuperclass %in% c(NA,0) &
            nvulpessubphylum>=1) %>%
  bind_cols(f_parasitesCTV_l="subphylum",f_parasitesCTV_v="Craniata",
            f_parasitesCTV="9") %>%
  select(qtax_family,f_parasitesCTV_v,f_parasitesCTV_l,f_parasitesCTV) -> 
  f_vsubphylum #20

globi_parasites_family %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasites=n(),
            nvulpesspecies = sum(globi_species=="Vulpes vulpes"),
            nvulpesgenus = sum(globi_genus=="Vulpes"),
            nvulpesfam = sum(globi_family=="Canidae"),
            nvulpessuborder = sum(globi_suborder=="Caniformia"),
            nvulpesorder = sum(globi_order=="Carnivora"),
            nvulpessuperorder = sum(globi_superorder=="Laurasiatheria"),
            nvulpesclass = sum(globi_class=="Mammalia"),
            nvulpessuperclass = sum(globi_superclass=="Sarcopterygii"),
            nvulpessubphylum = sum(globi_subphylum=="Craniata"),
            nvulpesphylum = sum(globi_phylum=="Chordata"),
            nvulpeskingdom = sum(globi_kingdom=="Metazoa"),
            nvulpessuperkingdom = sum(globi_superkingdom=="Eukaryota")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpessuborder %in% c(NA,0) & 
            nvulpesorder %in% c(NA,0) & nvulpessuperorder %in% c(NA,0) & 
            nvulpesclass%in% c(NA,0) & nvulpessuperclass %in% c(NA,0) &
            nvulpessubphylum %in% c(NA,0) &  nvulpesphylum >=1) %>%
  bind_cols(f_parasitesCTV_l="phylum",f_parasitesCTV_v="Chordata",
            f_parasitesCTV="10") %>%
  select(qtax_family,f_parasitesCTV_v,f_parasitesCTV_l,f_parasitesCTV) -> 
  f_vphylum #19

globi_parasites_family %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasites=n(),
            nvulpesspecies = sum(globi_species=="Vulpes vulpes"),
            nvulpesgenus = sum(globi_genus=="Vulpes"),
            nvulpesfam = sum(globi_family=="Canidae"),
            nvulpessuborder = sum(globi_suborder=="Caniformia"),
            nvulpesorder = sum(globi_order=="Carnivora"),
            nvulpessuperorder = sum(globi_superorder=="Laurasiatheria"),
            nvulpesclass = sum(globi_class=="Mammalia"),
            nvulpessuperclass = sum(globi_superclass=="Sarcopterygii"),
            nvulpessubphylum = sum(globi_subphylum=="Craniata"),
            nvulpesphylum = sum(globi_phylum=="Chordata"),
            nvulpeskingdom = sum(globi_kingdom=="Metazoa"),
            nvulpessuperkingdom = sum(globi_superkingdom=="Eukaryota")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpessuborder %in% c(NA,0) & 
            nvulpesorder %in% c(NA,0) & nvulpessuperorder %in% c(NA,0) & 
            nvulpesclass%in% c(NA,0) & nvulpessuperclass %in% c(NA,0) &
            nvulpessubphylum %in% c(NA,0) &  nvulpesphylum %in% c(NA,0) & 
            nvulpeskingdom>=1) %>%
  bind_cols(f_parasitesCTV_l="kingdom",f_parasitesCTV_v="Metazoa",
            f_parasitesCTV="11") %>%
  select(qtax_family,f_parasitesCTV_v,f_parasitesCTV_l,f_parasitesCTV) -> 
  f_vkingdom #63

globi_parasites_family %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasites=n(),
            nvulpesspecies = sum(globi_species=="Vulpes vulpes"),
            nvulpesgenus = sum(globi_genus=="Vulpes"),
            nvulpesfam = sum(globi_family=="Canidae"),
            nvulpessuborder = sum(globi_suborder=="Caniformia"),
            nvulpesorder = sum(globi_order=="Carnivora"),
            nvulpessuperorder = sum(globi_superorder=="Laurasiatheria"),
            nvulpesclass = sum(globi_class=="Mammalia"),
            nvulpessuperclass = sum(globi_superclass=="Sarcopterygii"),
            nvulpessubphylum = sum(globi_subphylum=="Craniata"),
            nvulpesphylum = sum(globi_phylum=="Chordata"),
            nvulpeskingdom = sum(globi_kingdom=="Metazoa"),
            nvulpessuperkingdom = sum(globi_superkingdom=="Eukaryota")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpessuborder %in% c(NA,0) & 
            nvulpesorder %in% c(NA,0) & nvulpessuperorder %in% c(NA,0) & 
            nvulpesclass%in% c(NA,0) & nvulpessuperclass %in% c(NA,0) &
            nvulpessubphylum %in% c(NA,0) &  nvulpesphylum %in% c(NA,0) & 
            nvulpeskingdom %in% c(NA,0) & nvulpessuperkingdom>=1) %>%
  bind_cols(f_parasitesCTV_l="superkingdom",f_parasitesCTV_v="Eukaryota",
            f_parasitesCTV="12") %>%
  select(qtax_family,f_parasitesCTV_v,f_parasitesCTV_l,f_parasitesCTV) -> 
  f_vsuperkingdom #151

globi_parasites_family %>% 
  group_by(qtax_family) %>% 
  summarize(interaction_f_parasites=n(),
            nvulpesspecies = sum(globi_species=="Vulpes vulpes"),
            nvulpesgenus = sum(globi_genus=="Vulpes"),
            nvulpesfam = sum(globi_family=="Canidae"),
            nvulpessuborder = sum(globi_suborder=="Caniformia"),
            nvulpesorder = sum(globi_order=="Carnivora"),
            nvulpessuperorder = sum(globi_superorder=="Laurasiatheria"),
            nvulpesclass = sum(globi_class=="Mammalia"),
            nvulpessuperclass = sum(globi_superclass=="Sarcopterygii"),
            nvulpessubphylum = sum(globi_subphylum=="Craniata"),
            nvulpesphylum = sum(globi_phylum=="Chordata"),
            nvulpeskingdom = sum(globi_kingdom=="Metazoa"),
            nvulpessuperkingdom = sum(globi_superkingdom=="Eukaryota")) %>%
  filter (nvulpesspecies %in% c(NA,0) & nvulpesgenus %in% c(NA,0) & 
            nvulpesfam %in% c(NA,0) & nvulpessuborder %in% c(NA,0) & 
            nvulpesorder %in% c(NA,0) & nvulpessuperorder %in% c(NA,0) & 
            nvulpesclass%in% c(NA,0) & nvulpessuperclass %in% c(NA,0) &
            nvulpessubphylum %in% c(NA,0) &  nvulpesphylum %in% c(NA,0) & 
            nvulpeskingdom %in% c(NA,0) & nvulpessuperkingdom %in% c(NA,0)) %>%
  bind_cols(f_parasitesCTV_l=">Superkingdom",f_parasitesCTV_v=">Eukaryota",
            f_parasitesCTV="13") %>%
  select(qtax_family,f_parasitesCTV_v,f_parasitesCTV_l,f_parasitesCTV) -> 
  f_nomatch #1

# no match bedeutet hier: Bacterium hat nur Interaktionen mit anderen Bacteria

f_parasites_CTV <- 
  data.frame(bind_rows(f_vspecies,f_vgenus,f_vfam,f_vsuborder,f_vorder,
                       f_vsuperorder,f_vclass,f_vsuperclass,f_vsubphylum,
                       f_vphylum,f_vkingdom,f_vsuperkingdom,f_nomatch)) #419

# LCA und closest to vulpes verbinden

g_f <- merge(f_parasites_LCA,f_parasites_CTV)

qtaxonomy_f <- data.frame(getTaxonomy(getId(g_f$qtax_family),
                                      desiredTaxa = desiredtaxa))
f_parasites2 <- cbind(g_f,
                      qtaxonomy_f$superkingdom,
                      qtaxonomy_f$kingdom, 
                      qtaxonomy_f$phylum,
                      qtaxonomy_f$subphylum,
                      qtaxonomy_f$superclass,
                      qtaxonomy_f$class,
                      qtaxonomy_f$subclass,
                      qtaxonomy_f$superorder,
                      qtaxonomy_f$order,
                      qtaxonomy_f$suborder,
                      qtaxonomy_f$family,
                      qtaxonomy_f$genus,
                      qtaxonomy_f$species)

f_parasites2 %>% rename(superkingdom = "qtaxonomy_f$superkingdom", 
                        kingdom = "qtaxonomy_f$kingdom",
                        phylum="qtaxonomy_f$phylum",
                        subphylum="qtaxonomy_f$subphylum",
                        superclass="qtaxonomy_f$superclass",
                        class="qtaxonomy_f$class",
                        subclass="qtaxonomy_f$subclass",
                        superorder="qtaxonomy_f$superorder",
                        order="qtaxonomy_f$order",
                        suborder="qtaxonomy_f$suborder",
                        family="qtaxonomy_f$family",
                        genus="qtaxonomy_f$genus",
                        species="qtaxonomy_f$species") %>%
  filter(!is.na(family)) %>%
  select(superkingdom,kingdom,phylum,subphylum, superclass,class,subclass,
         superorder,order,suborder,family,genus,species,
         interaction_f_parasites,f_parasitesLCA_v,f_parasitesLCA_l,f_parasitesLCA,
         f_parasitesCTV_v,f_parasitesCTV_l,f_parasitesCTV)-> family_parasites

saveRDS(family_parasites,"family_parasites.Rds")
