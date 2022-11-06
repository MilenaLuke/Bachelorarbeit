newglobi <- FALSE
if (newglobi) {
  source("1_globidata_parasites_genus")
} else {
  globidata_parasites_g <- readRDS("globi_parasites_genus.Rds")
}

library(tidyverse)
library(taxonomizr)

#Taxon Path von Vulpes Vulpes herausfinden
desiredtaxa = c("superkingdom","kingdom","superphylum","phylum",
                "subphylum", "superclass","class","subclass", "superorder",
                "order","suborder","family","genus","species")

vulpestaxonpath <- data.frame(getTaxonomy(getId("Vulpes vulpes")))

#nur Bakterien, Viren und Eukaryoten auswählen
globidata_parasites_g %>% 
  filter(globi_superkingdom %in% c("Bacteria","Eukaryota","Viruses")) -> 
  globi_parasites_genus2 #147818

# Doppelte entfernen
globi_parasites_genus <- unique(globi_parasites_genus2) #13462

globi_parasites_genus %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n()) ->x #501

#species LCA
globi_parasites_genus %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
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
            g_parasitesLCA_v=unique(globi_genus)) %>%
  filter (nspecies==1 & ngenus==1 & nfam==1 & nsuborder==1 & norder==1 &
            nsuperorder ==1 & nsubclass==1 & nclass==1 & nsuperclass==1 & 
            nsubphylum==1 & nphylum==1 & nkingdom==1 & nsuperkingdom==1 & 
            !is.na(g_parasitesLCA_v)) %>%
  bind_cols(g_parasitesLCA_l="species",g_parasitesLCA="1") %>%
  select(qtax_genus,
         interaction_g_parasites,
         g_parasitesLCA_v,
         g_parasitesLCA_l,
         g_parasitesLCA) -> g_speciesdata #124

#genus LCA
globi_parasites_genus %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
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
            g_parasitesLCA_v=unique(globi_genus)) %>%
  filter (nspecies>1 & ngenus==1 & nfam==1 & nsuborder==1 & norder==1 &
            nsuperorder ==1 & nsubclass==1 & nclass==1 & nsuperclass==1 & 
            nsubphylum==1 & nphylum==1 & nkingdom==1 & nsuperkingdom==1 & 
            !is.na(g_parasitesLCA_v)) %>%
  bind_cols(g_parasitesLCA_l="genus",g_parasitesLCA="2") %>%
  select(qtax_genus,
         interaction_g_parasites,
         g_parasitesLCA_v,
         g_parasitesLCA_l,
         g_parasitesLCA) -> g_genusdata #18

#family LCA
globi_parasites_genus %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
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
            g_parasitesLCA_v=unique(globi_family)) %>%
  filter (nspecies>1 & ngenus>1 & nfam==1 & nsuborder==1 & norder==1 &
            nsuperorder ==1 & nsubclass==1 & nclass==1 & nsuperclass==1 & 
            nsubphylum==1 & nphylum==1 & nkingdom==1 & nsuperkingdom==1 & 
            !is.na(g_parasitesLCA_v))%>%
  bind_cols(g_parasitesLCA_l="family",g_parasitesLCA="3") %>%
  select(qtax_genus,
         interaction_g_parasites,
         g_parasitesLCA_v,
         g_parasitesLCA_l,
         g_parasitesLCA) -> g_familydata #29

#suborder LCA 
globi_parasites_genus %>% 
  group_by(qtax_genus ) %>% 
  summarize(interaction_g_parasites=n(),
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
            g_parasitesLCA_v=unique(globi_suborder)) %>%
  filter (nspecies>1 & ngenus>1 & nfam>1 & nsuborder==1 & norder==1 &
            nsuperorder ==1 & nsubclass==1 & nclass==1 & nsuperclass==1 & 
            nsubphylum==1 & nphylum==1 & nkingdom==1 & nsuperkingdom==1 & 
            !is.na(g_parasitesLCA_v)) %>%
  bind_cols(g_parasitesLCA_l="suborder",g_parasitesLCA="4") %>%
  select(qtax_genus,
         interaction_g_parasites,
         g_parasitesLCA_v,
         g_parasitesLCA_l,
         g_parasitesLCA) -> g_suborderdata #6

#order LCA
globi_parasites_genus %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
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
            g_parasitesLCA_v=unique(globi_order))%>% 
  filter (nspecies>1 & ngenus>1 & nfam>1 & nsuborder>1 & norder==1 &
            nsuperorder ==1 & nsubclass==1 & nclass==1 & nsuperclass==1 & 
            nsubphylum==1 & nphylum==1 & nkingdom==1 & nsuperkingdom==1 & 
            !is.na(g_parasitesLCA_v)) %>%
  bind_cols(g_parasitesLCA_l="order",g_parasitesLCA="5") %>%
  select(qtax_genus,
         interaction_g_parasites,
         g_parasitesLCA_v,
         g_parasitesLCA_l,
         g_parasitesLCA) -> g_orderdata #3

#superorder LCA
globi_parasites_genus %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
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
            g_parasitesLCA_v=unique(globi_superorder)) %>%
  filter (nspecies!=1 & ngenus!=1 & nfam!=1 & nsuborder!=1 & norder!=1 &
            nsuperorder==1 & nsubclass==1 & nclass==1 & nsuperclass==1 & 
            nsubphylum==1 & nphylum==1 & nkingdom==1 & nsuperkingdom==1 &
            !is.na(g_parasitesLCA_v)) %>%
  bind_cols(g_parasitesLCA_l="superorder",g_parasitesLCA="6") %>%
  select(qtax_genus,
         interaction_g_parasites,
         g_parasitesLCA_v,
         g_parasitesLCA_l,
         g_parasitesLCA) -> g_superorderdata #6

#subclass LCA
globi_parasites_genus %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
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
            g_parasitesLCA_v=unique(globi_subclass)) %>%
  filter (nspecies!=1 & ngenus!=1 & nfam!=1 & nsuborder!=1 & norder!=1 &
            nsuperorder!=1 & nsubclass==1 & nclass==1 & 
            nsuperclass==1 & nsubphylum==1 & nphylum==1 & nkingdom==1 & 
            nsuperkingdom==1 & !is.na(g_parasitesLCA_v)) %>%
  bind_cols(g_parasitesLCA_l="subclass",g_parasitesLCA="7") %>%
  select(qtax_genus,
         interaction_g_parasites,
         g_parasitesLCA_v,
         g_parasitesLCA_l,
         g_parasitesLCA) -> g_subclassdata #7

#class LCA
globi_parasites_genus %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
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
            g_parasitesLCA_v=unique(globi_class)) %>%
  filter (nspecies!=1 & ngenus!=1 & nfam!=1 & nsuborder!=1 & norder!=1 &
            nsuperorder!=1 & nsubclass==1 & nclass==1 & 
            nsuperclass==1 & nsubphylum==1 & nphylum==1 & nkingdom==1 & 
            nsuperkingdom==1 & is.na(g_parasitesLCA_v)) %>%
  bind_cols(g_parasitesLCA_l="class",g_parasitesLCA="7") %>%
  select(qtax_genus,
         interaction_g_parasites,
         g_parasitesLCA_v,
         g_parasitesLCA_l,
         g_parasitesLCA) -> g_classdata #24

#superclass LCA
globi_parasites_genus %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
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
            g_parasitesLCA_v=unique(globi_superclass)) %>%
  filter (nspecies!=1 & ngenus!=1 & nfam!=1 & nsuborder!=1 & norder!=1 &
            nsuperorder!=1 & nsubclass!=1 & nclass!=1 & nsuperclass==1 & 
            nsubphylum==1 & nphylum==1 & nkingdom==1 & nsuperkingdom==1 & 
            !is.na(g_parasitesLCA_v)) %>%
  bind_cols(g_parasitesLCA_l="superclass",g_parasitesLCA="9") %>%
  select(qtax_genus,
         interaction_g_parasites,
         g_parasitesLCA_v,
         g_parasitesLCA_l,
         g_parasitesLCA) -> g_superclassdata #2

#subphylum LCA
globi_parasites_genus %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
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
            g_parasitesLCA_v=unique(globi_subphylum)) %>%
  filter (nspecies!=1 & ngenus!=1 & nfam!=1 & nsuborder!=1 & norder!=1 &
            nsuperorder!=1 & nsubclass!=1 & nclass!=1, nsuperclass!=1 & 
            nsubphylum==1 & nphylum==1 & nkingdom==1 & nsuperkingdom==1 &
            !is.na(g_parasitesLCA_v)) %>%
  bind_cols(g_parasitesLCA_l="subphylum",g_parasitesLCA="10") %>%
  select(qtax_genus,
         interaction_g_parasites,
         g_parasitesLCA_v,
         g_parasitesLCA_l,
         g_parasitesLCA) -> g_subphylumdata #7

#phylum LCA
globi_parasites_genus %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
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
            g_parasitesLCA_v=unique(globi_phylum)) %>%
  filter (nspecies!=1 & ngenus!=1 & nfam!=1 & nsuborder!=1 & norder!=1 &
            nsuperorder!=1 & nsubclass!=1 & nclass!=1 & nsuperclass!=1 & 
            nsubphylum!=1 & nphylum==1 & nkingdom==1 & nsuperkingdom==1 &
            !is.na(g_parasitesLCA_v)) %>%
  bind_cols(g_parasitesLCA_l="phylum",g_parasitesLCA="11") %>%
  select(qtax_genus,
         interaction_g_parasites,
         g_parasitesLCA_v,
         g_parasitesLCA_l,
         g_parasitesLCA) -> g_phylumdata #3

#kingdom LCA
globi_parasites_genus %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
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
            g_parasitesLCA_v=unique(globi_kingdom)) %>%
  filter (nspecies!=1 & ngenus!=1 & nfam!=1 & nsuborder!=1 & norder!=1 &
            nsuperorder!=1 & nsubclass!=1 & nclass!=1 & nsuperclass!=1 & 
            nsubphylum!=1 & nphylum!=1 & nkingdom==1 & nsuperkingdom==1 &
            !is.na(g_parasitesLCA_v)) %>%
  bind_cols(g_parasitesLCA_l="kingdom",g_parasitesLCA="12") %>%
  select(qtax_genus,
         interaction_g_parasites,
         g_parasitesLCA_v,
         g_parasitesLCA_l,
         g_parasitesLCA) -> g_kingdomdata #39

#superkingdom LCA
globi_parasites_genus %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
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
            g_parasitesLCA_v=unique(globi_superkingdom)) %>%
  filter (nspecies!=1 & ngenus!=1 & nfam!=1 & nsuborder!=1 & norder!=1 &
            nsuperorder!=1 & nsubclass!=1 & nclass!=1, nsuperclass!=1 & 
            nsubphylum!=1,nphylum!=1,nkingdom!=1,nsuperkingdom==1 & 
            !is.na(g_parasitesLCA_v)) %>%
  bind_cols(g_parasitesLCA_l="superkingdom",g_parasitesLCA="13") %>%
  select(qtax_genus,
         interaction_g_parasites,
         g_parasitesLCA_v,
         g_parasitesLCA_l,
         g_parasitesLCA) -> g_superkingdomdata #20

#no LCA
globi_parasites_genus %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
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
            g_parasitesLCA_v=unique(globi_superkingdom)) %>%
  filter (nspecies!=1 & ngenus!=1 & nfam!=1 & nsuborder!=1 & norder!=1 &
            nsuperorder!=1 & nsubclass!=1 & nclass!=1, nsuperclass!=1 & 
            nsubphylum!=1,nphylum!=1,nkingdom!=1,nsuperkingdom>1 &
            !is.na(g_parasitesLCA_v)) %>%
  bind_cols(g_parasitesLCA_l=">Superkingdom",g_parasitesLCA="14") %>%
  select(qtax_genus,
         interaction_g_parasites,
         g_parasitesLCA_v,
         g_parasitesLCA_l,
         g_parasitesLCA) -> g_nodata #50

g_parasites_LCA <- 
  unique(as.data.frame(bind_rows(g_speciesdata,g_genusdata,g_familydata,
                                 g_suborderdata,g_orderdata,g_superorderdata,
                                 g_subclassdata,g_classdata,g_superclassdata, 
                                 g_subphylumdata,g_phylumdata,g_kingdomdata,
                                 g_superkingdomdata,g_nodata))) #338

#--- Closest to vulpes level + value finden

globi_parasites_genus %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
            nvulpesspecies = sum(globi_species=="Vulpes vulpes")) %>%
  filter(nvulpesspecies>=1) %>%
  bind_cols(g_parasitesCTV_l="species",g_parasitesCTV_v="Vulpes vulpes",
            g_parasitesCTV="1") %>%
  select(qtax_genus,g_parasitesCTV_v,g_parasitesCTV_l,g_parasitesCTV) ->
  g_vspecies #1

globi_parasites_genus %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasite=n(),
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
  bind_cols(g_parasitesCTV_l="genus",g_parasitesCTV_v="Vulpes",
            g_parasitesCTV="2") %>%
  select(qtax_genus,g_parasitesCTV_v,g_parasitesCTV_l,g_parasitesCTV) -> 
  g_vgenus #6

globi_parasites_genus %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
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
  bind_cols(g_parasitesCTV_l="family",g_parasitesCTV_v="Canidae",
            g_parasitesCTV="3") %>%
  select(qtax_genus,g_parasitesCTV_v,g_parasitesCTV_l,g_parasitesCTV) -> 
  g_vfam #26

globi_parasites_genus %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
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
  bind_cols(g_parasitesCTV_l="order",g_parasitesCTV_v="Caniformia",
            g_parasitesCTV="5") %>%
  select(qtax_genus,g_parasitesCTV_v,g_parasitesCTV_l,g_parasitesCTV) -> 
  g_vsuborder #2

globi_parasites_genus %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
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
  bind_cols(g_parasitesCTV_l="order",g_parasitesCTV_v="Carnivora",
            g_parasitesCTV="5") %>%
  select(qtax_genus,g_parasitesCTV_v,g_parasitesCTV_l,g_parasitesCTV) -> 
  g_vorder #38

globi_parasites_genus %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
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
  bind_cols(g_parasitesCTV_l="superorder",g_parasitesCTV_v="Laurasiatheria",
            g_parasitesCTV="6") %>%
  select(qtax_genus,g_parasitesCTV_v,g_parasitesCTV_l,g_parasitesCTV) -> 
  g_vsuperorder #26

globi_parasites_genus %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
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
  bind_cols(g_parasitesCTV_l="class",g_parasitesCTV_v="Mammalia",
            g_parasitesCTV="7") %>%
  select(qtax_genus,g_parasitesCTV_v,g_parasitesCTV_l,g_parasitesCTV) -> 
  g_vclass #102

globi_parasites_genus %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
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
  bind_cols(g_parasitesCTV_l="superclass",g_parasitesCTV_v="Sarcopterygii",
            g_parasitesCTV="8") %>%
  select(qtax_genus,g_parasitesCTV_v,g_parasitesCTV_l,g_parasitesCTV) -> 
  g_vsuperclass #9

globi_parasites_genus %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
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
  bind_cols(g_parasitesCTV_l="subphylum",g_parasitesCTV_v="Craniata",
            g_parasitesCTV="9") %>%
  select(qtax_genus,g_parasitesCTV_v,g_parasitesCTV_l,g_parasitesCTV) -> 
  g_vsubphylum #10

globi_parasites_genus %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
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
  bind_cols(g_parasitesCTV_l="phylum",g_parasitesCTV_v="Chordata",
            g_parasitesCTV="10") %>%
  select(qtax_genus,g_parasitesCTV_v,g_parasitesCTV_l,g_parasitesCTV) -> 
  g_vphylum #12

globi_parasites_genus %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
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
  bind_cols(g_parasitesCTV_l="kingdom",g_parasitesCTV_v="Metazoa",
            g_parasitesCTV="11") %>%
  select(qtax_genus,g_parasitesCTV_v,g_parasitesCTV_l,g_parasitesCTV) -> 
  g_vkingdom #59

globi_parasites_genus %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
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
  bind_cols(g_parasitesCTV_l="superkingdom",g_parasitesCTV_v="Eukaryota",
            g_parasitesCTV="12") %>%
  select(qtax_genus,g_parasitesCTV_v,g_parasitesCTV_l,g_parasitesCTV) -> 
  g_vsuperkingdom #207

globi_parasites_genus %>% 
  group_by(qtax_genus) %>% 
  summarize(interaction_g_parasites=n(),
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
  bind_cols(g_parasitesCTV_l=">Superkingdom",g_parasitesCTV_v=">Eukaryota",
            g_parasitesCTV="13") %>%
  select(qtax_genus,g_parasitesCTV_v,g_parasitesCTV_l,g_parasitesCTV) -> 
  g_nomatch #3

# no match bedeutet hier: Bacterium hat nur Interaktionen mit anderen Bacteria

g_parasites_CTV <- 
  data.frame(bind_rows(g_vspecies,g_vgenus,g_vfam,g_vsuborder,g_vorder,
                       g_vsuperorder,g_vclass,g_vsuperclass,g_vsubphylum,
                       g_vphylum,g_vkingdom,g_vsuperkingdom,g_nomatch)) #501

# LCA und closest to vulpes verbinden

g_p <- merge(g_parasites_LCA,g_parasites_CTV)

qtaxonomy_g <- data.frame(getTaxonomy(getId(g_p$qtax_genus),
                                      desiredTaxa = desiredtaxa))
g_parasites2 <- cbind(g_p,
                      qtaxonomy_g$superkingdom,
                      qtaxonomy_g$kingdom, 
                      qtaxonomy_g$phylum,
                      qtaxonomy_g$subphylum,
                      qtaxonomy_g$superclass,
                      qtaxonomy_g$class,
                      qtaxonomy_g$subclass,
                      qtaxonomy_g$superorder,
                      qtaxonomy_g$order,
                      qtaxonomy_g$suborder,
                      qtaxonomy_g$family,
                      qtaxonomy_g$genus,
                      qtaxonomy_g$species)

g_parasites2 %>% rename(superkingdom = "qtaxonomy_g$superkingdom", 
                        kingdom = "qtaxonomy_g$kingdom",
                        phylum="qtaxonomy_g$phylum",
                        subphylum="qtaxonomy_g$subphylum",
                        superclass="qtaxonomy_g$superclass",
                        class="qtaxonomy_g$class",
                        subclass="qtaxonomy_g$subclass",
                        superorder="qtaxonomy_g$superorder",
                        order="qtaxonomy_g$order",
                        suborder="qtaxonomy_g$suborder",
                        family="qtaxonomy_g$family",
                        genus="qtaxonomy_g$genus",
                        species="qtaxonomy_g$species") %>%
  filter(!is.na(genus)) %>%
  select(superkingdom,kingdom,phylum,subphylum, superclass,class,subclass,
         superorder,order,suborder,family,genus,species,
         interaction_g_parasites,g_parasitesLCA_v,g_parasitesLCA_l,g_parasitesLCA,
         g_parasitesCTV_v,g_parasitesCTV_l,g_parasitesCTV)-> genus_parasites

saveRDS(genus_parasites,"genus_parasites.Rds")
