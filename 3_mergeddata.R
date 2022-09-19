globi_p_f <- FALSE
if (globi_p_f) {
  source("Add_LCA_CTV_parasites_family.R")
} else {
  f_parasites <- readRDS("f_parasites.Rds")
}

globi_e_f <- FALSE
if (globi_e_f) {
  source("Add_LCA_CTV_eat_family.R")
} else {
  f_eat <- readRDS("f_eat.Rds")
}

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

family_data <- merge(f_parasites,f_eat,by.x="qtax_family",by.y="qtax_family",
                     all=TRUE)


# plotten
s <- s_parasites_numbers
x <- s$s_CTV
y <- s$s_LCA
relation <- lm(x~y)
plot(y,x,col="blue",main="Überschrift", abline(lm(x~y), cex=1.3,pch=16,
                                               xlab="CTV",ylab="LCA"))

plot(x,y,col="blue",main="Überschrift", abline(lm(x~y),
                                               xlab="CTV",ylab="LCA"))
