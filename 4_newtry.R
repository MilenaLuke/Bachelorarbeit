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


q <- unique(globidata_parasites_s) #9022

q %>% filter(qtax_species =="Angiostrongylus costaricensis") %>%
  summarize(nclass=length(globi_class),
            numberclass=length(unique(globi_class)))-> b

# nur 1 Klasse
q %>% 
  group_by(qtax_species) %>%
  filter(!is.na(globi_class)) %>%
  summarize(n_interaction=length(globi_species),
            numberclass=length(unique(globi_class))) %>%
  filter(numberclass==1) %>%
  bind_cols(class="nur 1 Klasse")-> x

#mehr als 1 Klasse
q %>% 
  group_by(qtax_species) %>%
  summarize(n_interaction=length(globi_species),
            n_class=length(unique(globi_class)),
            class=globi_class) %>%
  filter(n_class>1 & !is.na(class)) -> y


z <-unique(y)

list_class <- as.data.frame(unique(q$globi_class))

%>%
  filter(!is.na(class)) 

%>%
  summarize(n= sum(Factor="Mammalia")) ->x
  
%>% 
  summarize(interaction_s_parasites=n())
