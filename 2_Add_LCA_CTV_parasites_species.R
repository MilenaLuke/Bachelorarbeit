newglobi <- FALSE
if (newglobi) {
  source("1_globidata_parasites_species")
} else {
  globidata_parasites_s <- readRDS("globi_parasites_species.Rds")
}

library(tidyverse)
library(taxonomizr)

#Taxon Path von Vulpes Vulpes herausfinden
desiredtaxa = c("superkingdom","kingdom","superphylum","phylum",
                "subphylum", "superclass","class","subclass", "superorder",
                "order","suborder","family","genus","species")

vulpestaxonpath <- data.frame(getTaxonomy(getId("Vulpes vulpes"),
                                          desiredTaxa=desiredtaxa))

#nur Bakterien, Viren und Eukaryoten auswählen
globidata_parasites_s %>% 
  filter(globi_superkingdom %in% c("Bacteria","Eukaryota","Viruses")) -> 
  globi_parasites_species2 #9022

# Doppelte entfernen
globi_parasites_species <- unique(globi_parasites_species2) #3546

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n()) ->x #419

#species LCA
globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
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
            s_parasitesLCA_v=unique(globi_species)) %>%
  filter (nspecies==1 & ngenus==1 & nfam==1 & nsuborder==1 & norder==1 &
            nsuperorder ==1 & nsubclass==1 & nclass==1 & nsuperclass==1 & 
            nsubphylum==1 & nphylum==1 & nkingdom==1 & nsuperkingdom==1 & 
            !is.na(s_parasitesLCA_v)) %>%
  bind_cols(s_parasitesLCA_l="species",s_parasitesLCA="1") %>%
  select(qtax_species,
         interaction_s_parasites,
         s_parasitesLCA_v,
         s_parasitesLCA_l,
         s_parasitesLCA) -> s_speciesdata #170

#genus LCA
globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
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
            s_parasitesLCA_v=unique(globi_genus)) %>%
  filter (nspecies >1 & ngenus==1 & nfam==1 & nsuborder==1 & norder==1 &
            nsuperorder ==1 & nsubclass==1 & nclass==1 & nsuperclass==1 & 
            nsubphylum==1 & nphylum==1 & nkingdom==1 & nsuperkingdom==1 & 
            !is.na(s_parasitesLCA_v)) %>%
  bind_cols(s_parasitesLCA_l="genus",s_parasitesLCA="2") %>%
  select(qtax_species,
         interaction_s_parasites,
         s_parasitesLCA_v,
         s_parasitesLCA_l,
         s_parasitesLCA) -> s_genusdata #22

#family LCA
globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
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
            s_parasitesLCA_v=unique(globi_family)) %>%
  filter (nspecies>1 & ngenus>1 & nfam==1 & nsuborder==1 & norder==1 &
            nsuperorder ==1 & nsubclass==1 & nclass==1 & nsuperclass==1 & 
            nsubphylum==1 & nphylum==1 & nkingdom==1 & nsuperkingdom==1 & 
            !is.na(s_parasitesLCA_v))%>%
  bind_cols(s_parasitesLCA_l="family",s_parasitesLCA="3") %>%
  select(qtax_species,
         interaction_s_parasites,
         s_parasitesLCA_v,
         s_parasitesLCA_l,
         s_parasitesLCA) -> s_familydata #32

#suborder LCA 
globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
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
            s_parasitesLCA_v=unique(globi_suborder)) %>%
  filter (nspecies>1 & ngenus>1 & nfam>1 & nsuborder==1 & norder==1 &
            nsuperorder ==1 & nsubclass==1 & nclass==1 & nsuperclass==1 & 
            nsubphylum==1 & nphylum==1 & nkingdom==1 & nsuperkingdom==1 & 
            !is.na(s_parasitesLCA_v)) %>%
  bind_cols(s_parasitesLCA_l="suborder",s_parasitesLCA="4") %>%
  select(qtax_species,
         interaction_s_parasites,
         s_parasitesLCA_v,
         s_parasitesLCA_l,
         s_parasitesLCA) -> s_suborderdata #7

#order LCA
globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
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
            s_parasitesLCA_v=unique(globi_order))%>% 
  filter (nspecies>1 & ngenus>1 & nfam>1 & nsuborder>1 & norder==1 &
            nsuperorder ==1 & nsubclass==1 & nclass==1 & nsuperclass==1 & 
            nsubphylum==1 & nphylum==1 & nkingdom==1 & nsuperkingdom==1 & 
            !is.na(s_parasitesLCA_v)) %>%
  bind_cols(s_parasitesLCA_l="order",s_parasitesLCA="5") %>%
  select(qtax_species,
         interaction_s_parasites,
         s_parasitesLCA_v,
         s_parasitesLCA_l,
         s_parasitesLCA) -> s_orderdata #7

#superorder LCA
globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
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
            s_parasitesLCA_v=unique(globi_superorder)) %>%
  filter (nspecies!=1 & ngenus!=1 & nfam!=1 & nsuborder!=1 & norder!=1 &
          nsuperorder==1 & nsubclass==1 & nclass==1 & nsuperclass==1 & 
            nsubphylum==1 & nphylum==1 & nkingdom==1 & nsuperkingdom==1 &
            !is.na(s_parasitesLCA_v)) %>%
  bind_cols(s_parasitesLCA_l="superorder",s_parasitesLCA="6") %>%
  select(qtax_species,
         interaction_s_parasites,
         s_parasitesLCA_v,
         s_parasitesLCA_l,
         s_parasitesLCA) -> s_superorderdata #7

#subclass LCA
globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
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
            s_parasitesLCA_v=unique(globi_subclass)) %>%
  filter (nspecies!=1 & ngenus!=1 & nfam!=1 & nsuborder!=1 & norder!=1 &
          nsuperorder!=1 & nsubclass==1 & nclass==1 & 
            nsuperclass==1 & nsubphylum==1 & nphylum==1 & nkingdom==1 & 
            nsuperkingdom==1 & !is.na(s_parasitesLCA_v)) %>%
  bind_cols(s_parasitesLCA_l="subclass",s_parasitesLCA="7") %>%
  select(qtax_species,
         interaction_s_parasites,
         s_parasitesLCA_v,
         s_parasitesLCA_l,
         s_parasitesLCA) -> s_subclassdata #3

#class LCA
globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
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
            s_parasitesLCA_v=unique(globi_class)) %>%
  filter (nspecies!=1 & ngenus!=1 & nfam!=1 & nsuborder!=1 & norder!=1 &
            nsuperorder!=1 & nsubclass==1 & nclass==1 & 
            nsuperclass==1 & nsubphylum==1 & nphylum==1 & nkingdom==1 & 
            nsuperkingdom==1 & is.na(s_parasitesLCA_v)) %>%
  bind_cols(s_parasitesLCA_l="class",s_parasitesLCA="7") %>%
  select(qtax_species,
         interaction_s_parasites,
         s_parasitesLCA_v,
         s_parasitesLCA_l,
         s_parasitesLCA) -> s_classdata #39

#superclass LCA
globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
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
            s_parasitesLCA_v=unique(globi_superclass)) %>%
  filter (nspecies!=1 & ngenus!=1 & nfam!=1 & nsuborder!=1 & norder!=1 &
            nsuperorder!=1 & nsubclass!=1 & nclass!=1 & nsuperclass==1 & 
            nsubphylum==1 & nphylum==1 & nkingdom==1 & nsuperkingdom==1 & 
            !is.na(s_parasitesLCA_v)) %>%
  bind_cols(s_parasitesLCA_l="superclass",s_parasitesLCA="9") %>%
  select(qtax_species,
         interaction_s_parasites,
         s_parasitesLCA_v,
         s_parasitesLCA_l,
         s_parasitesLCA) -> s_superclassdata #0

#subphylum LCA
globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
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
            s_parasitesLCA_v=unique(globi_subphylum)) %>%
  filter (nspecies!=1 & ngenus!=1 & nfam!=1 & nsuborder!=1 & norder!=1 &
            nsuperorder!=1 & nsubclass!=1 & nclass!=1, nsuperclass!=1 & 
            nsubphylum==1 & nphylum==1 & nkingdom==1 & nsuperkingdom==1 &
            !is.na(s_parasitesLCA_v)) %>%
  bind_cols(s_parasitesLCA_l="subphylum",s_parasitesLCA="10") %>%
  select(qtax_species,
         interaction_s_parasites,
         s_parasitesLCA_v,
         s_parasitesLCA_l,
         s_parasitesLCA) -> s_subphylumdata #9

#phylum LCA
globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
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
            s_parasitesLCA_v=unique(globi_phylum)) %>%
  filter (nspecies!=1 & ngenus!=1 & nfam!=1 & nsuborder!=1 & norder!=1 &
            nsuperorder!=1 & nsubclass!=1 & nclass!=1 & nsuperclass!=1 & 
            nsubphylum!=1 & nphylum==1 & nkingdom==1 & nsuperkingdom==1 &
            !is.na(s_parasitesLCA_v)) %>%
  bind_cols(s_parasitesLCA_l="phylum",s_parasitesLCA="11") %>%
  select(qtax_species,
         interaction_s_parasites,
         s_parasitesLCA_v,
         s_parasitesLCA_l,
         s_parasitesLCA) -> s_phylumdata #0

#kingdom LCA
globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
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
            s_parasitesLCA_v=unique(globi_kingdom)) %>%
  filter (nspecies!=1 & ngenus!=1 & nfam!=1 & nsuborder!=1 & norder!=1 &
            nsuperorder!=1 & nsubclass!=1 & nclass!=1 & nsuperclass!=1 & 
            nsubphylum!=1 & nphylum!=1 & nkingdom==1 & nsuperkingdom==1 &
            !is.na(s_parasitesLCA_v)) %>%
  bind_cols(s_parasitesLCA_l="kingdom",s_parasitesLCA="12") %>%
  select(qtax_species,
         interaction_s_parasites,
         s_parasitesLCA_v,
         s_parasitesLCA_l,
         s_parasitesLCA) -> s_kingdomdata #27

#superkingdom LCA
globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
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
            s_parasitesLCA_v=unique(globi_superkingdom)) %>%
  filter (nspecies!=1 & ngenus!=1 & nfam!=1 & nsuborder!=1 & norder!=1 &
            nsuperorder!=1 & nsubclass!=1 & nclass!=1, nsuperclass!=1 & 
            nsubphylum!=1,nphylum!=1,nkingdom!=1,nsuperkingdom==1 & 
            !is.na(s_parasitesLCA_v)) %>%
  bind_cols(s_parasitesLCA_l="superkingdom",s_parasitesLCA="13") %>%
  select(qtax_species,
         interaction_s_parasites,
         s_parasitesLCA_v,
         s_parasitesLCA_l,
         s_parasitesLCA) -> s_superkingdomdata #11

#no LCA
globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
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
            s_parasitesLCA_v=unique(globi_superkingdom)) %>%
  filter (nspecies!=1 & ngenus!=1 & nfam!=1 & nsuborder!=1 & norder!=1 &
            nsuperorder!=1 & nsubclass!=1 & nclass!=1, nsuperclass!=1 & 
            nsubphylum!=1,nphylum!=1,nkingdom!=1,nsuperkingdom>1 &
            !is.na(s_parasitesLCA_v)) %>%
  bind_cols(s_parasitesLCA_l=">Superkingdom",s_parasitesLCA="14") %>%
  select(qtax_species,
         interaction_s_parasites,
         s_parasitesLCA_v,
         s_parasitesLCA_l,
         s_parasitesLCA) -> s_nodata #2

s_parasites_LCA <- 
  unique(as.data.frame(bind_rows(s_speciesdata,s_genusdata,s_familydata,
                                 s_suborderdata,s_orderdata,s_superorderdata,
                                 s_subclassdata,s_classdata,s_superclassdata, 
                                 s_subphylumdata,s_phylumdata,s_kingdomdata,
                                 s_superkingdomdata,s_nodata))) #336

#--- Closest to vulpes level + value finden

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
            nvulpesspecies = sum(globi_species=="Vulpes vulpes")) %>%
  filter(nvulpesspecies>=1) %>%
  bind_cols(s_parasitesCTV_l="species",s_parasitesCTV_v="Vulpes vulpes",
            s_parasitesCTV="1") %>%
  select(qtax_species,s_parasitesCTV_v,s_parasitesCTV_l,s_parasitesCTV) ->
  s_vspecies #7

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasite=n(),
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
  bind_cols(s_parasitesCTV_l="genus",s_parasitesCTV_v="Vulpes",
            s_parasitesCTV="2") %>%
  select(qtax_species,s_parasitesCTV_v,s_parasitesCTV_l,s_parasitesCTV) -> 
  s_vgenus #10

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(nvulpesspecies = sum(globi_species=="Vulpes vulpes"),
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
  bind_cols(s_parasitesCTV_l="family",s_parasitesCTV_v="Canidae",
            s_parasitesCTV="3") %>%
  select(qtax_species,s_parasitesCTV_v,s_parasitesCTV_l,s_parasitesCTV) -> 
  s_vfam #21

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
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
  bind_cols(s_parasitesCTV_l="order",s_parasitesCTV_v="Caniformia",
            s_parasitesCTV="5") %>%
  select(qtax_species,s_parasitesCTV_v,s_parasitesCTV_l,s_parasitesCTV) -> 
  s_vsuborder #9

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
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
  bind_cols(s_parasitesCTV_l="order",s_parasitesCTV_v="Carnivora",
            s_parasitesCTV="5") %>%
  select(qtax_species,s_parasitesCTV_v,s_parasitesCTV_l,s_parasitesCTV) -> 
  s_vorder #38

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
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
  bind_cols(s_parasitesCTV_l="superorder",s_parasitesCTV_v="Laurasiatheria",
            s_parasitesCTV="6") %>%
  select(qtax_species,s_parasitesCTV_v,s_parasitesCTV_l,s_parasitesCTV) -> 
  s_vsuperorder #32

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
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
  bind_cols(s_parasitesCTV_l="class",s_parasitesCTV_v="Mammalia",
            s_parasitesCTV="7") %>%
  select(qtax_species,s_parasitesCTV_v,s_parasitesCTV_l,s_parasitesCTV) -> 
  s_vclass #109

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
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
  bind_cols(s_parasitesCTV_l="superclass",s_parasitesCTV_v="Sarcopterygii",
            s_parasitesCTV="8") %>%
  select(qtax_species,s_parasitesCTV_v,s_parasitesCTV_l,s_parasitesCTV) -> 
  s_vsuperclass #10

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
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
  bind_cols(s_parasitesCTV_l="subphylum",s_parasitesCTV_v="Craniata",
            s_parasitesCTV="9") %>%
  select(qtax_species,s_parasitesCTV_v,s_parasitesCTV_l,s_parasitesCTV) -> 
  s_vsubphylum #5

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
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
  bind_cols(s_parasitesCTV_l="phylum",s_parasitesCTV_v="Chordata",
            s_parasitesCTV="10") %>%
  select(qtax_species,s_parasitesCTV_v,s_parasitesCTV_l,s_parasitesCTV) -> 
  s_vphylum #5

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
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
  bind_cols(s_parasitesCTV_l="kingdom",s_parasitesCTV_v="Metazoa",
            s_parasitesCTV="11") %>%
  select(qtax_species,s_parasitesCTV_v,s_parasitesCTV_l,s_parasitesCTV) -> 
  s_vkingdom #41

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
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
  bind_cols(s_parasitesCTV_l="superkingdom",s_parasitesCTV_v="Eukaryota",
            s_parasitesCTV="12") %>%
  select(qtax_species,s_parasitesCTV_v,s_parasitesCTV_l,s_parasitesCTV) -> 
  s_vsuperkingdom #132

globi_parasites_species %>% 
  group_by(qtax_species) %>% 
  summarize(interaction_s_parasites=n(),
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
  bind_cols(s_parasitesCTV_l=">Superkingdom",s_parasitesCTV_v=">Eukaryota",
            s_parasitesCTV="13") %>%
  select(qtax_species,s_parasitesCTV_v,s_parasitesCTV_l,s_parasitesCTV) -> 
  s_nomatch #0

# no match bedeutet hier: Bacterium hat nur Interaktionen mit anderen Bacteria

s_parasites_CTV <- 
  data.frame(bind_rows(s_vspecies,s_vgenus,s_vfam,s_vsuborder,s_vorder,
                       s_vsuperorder,s_vclass,s_vsuperclass,s_vsubphylum,
                       s_vphylum,s_vkingdom,s_vsuperkingdom,s_nomatch)) #419

# LCA und closest to vulpes verbinden

s_parasites <- merge(s_parasites_LCA,s_parasites_CTV)

qtaxonomy_s <- data.frame(getTaxonomy(getId(s_parasites$qtax_species),
                                      desiredTaxa = desiredtaxa))
s_parasites2 <- cbind(s_parasites,
                      qtaxonomy_s$superkingdom,
                      qtaxonomy_s$kingdom, 
                      qtaxonomy_s$phylum,
                      qtaxonomy_s$subphylum,
                      qtaxonomy_s$superclass,
                      qtaxonomy_s$class,
                      qtaxonomy_s$subclass,
                      qtaxonomy_s$superorder,
                      qtaxonomy_s$order,
                      qtaxonomy_s$suborder,
                      qtaxonomy_s$family,
                      qtaxonomy_s$genus,
                      qtaxonomy_s$species)

s_parasites2 %>% rename(superkingdom = "qtaxonomy_s$superkingdom", 
                        kingdom = "qtaxonomy_s$kingdom",
                        phylum="qtaxonomy_s$phylum",
                        subphylum="qtaxonomy_s$subphylum",
                        superclass="qtaxonomy_s$superclass",
                        class="qtaxonomy_s$class",
                        subclass="qtaxonomy_s$subclass",
                        superorder="qtaxonomy_s$superorder",
                        order="qtaxonomy_s$order",
                        suborder="qtaxonomy_s$suborder",
                        family="qtaxonomy_s$family",
                        genus="qtaxonomy_s$genus",
                        species="qtaxonomy_s$species") %>%
  filter(!is.na(species)) %>%
  select(superkingdom,kingdom,phylum,subphylum, superclass,class,subclass,
         superorder,order,suborder,family,genus,species,
         interaction_s_parasites,s_parasitesLCA_v,s_parasitesLCA_l,s_parasitesLCA,
         s_parasitesCTV_v,s_parasitesCTV_l,s_parasitesCTV)-> species_parasites

saveRDS(species_parasites,"species_parasites.Rds")
