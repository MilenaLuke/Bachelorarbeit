globi_p_g <- FALSE
if (globi_p_g) {
  source("Add_LCA_CTV_parasites_genus.R")
} else {
  g_parasites <- readRDS("g_parasites.Rds")
}

globi_e_g <- FALSE
if (globi_e_g) {
  source("Add_LCA_CTV_eat_genus.R")
} else {
  g_eat <- readRDS("g_eat.Rds")
}

globi_p_s <- FALSE
if (globi_p_s) {
  source("Add_LCA_CTV_parasites_species.R")
} else {
  s_parasites <- readRDS("s_parasites.Rds")
}

globi_e_s <- FALSE
if (globi_e_s) {
  source("Add_LCA_CTV_eat_species.R")
} else {
  s_eat <- readRDS("s_eat.Rds")
}

library(tidyverse)
library(taxonomizr)

#Daten von Nahrung und Parasiten zusammenfügen (genus und species separat)
genus_data <- merge(g_parasites,g_eat, by.x="qtax_genus",by.y="qtax_genus", 
                    all=TRUE)

species_data <- merge(s_parasites,s_eat,by.x="qtax_species",by.y="qtax_species",
                      all=TRUE)
