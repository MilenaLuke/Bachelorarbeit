library(rglobi)
library(dplyr)
library(taxonomizr)
library(data.table)

foxdata <- read.csv(file="ForMilena.csv")

# Doppelte + NA + unidentified/uncultured entfernen
unique_foxdata <- unique(foxdata)
fdata_species <- unique_foxdata[!is.na(unique_foxdata$species),]
fd_species <- fdata_species %>% filter(!grepl('unidentified|uncultured',species))

# Für jeden Eukaryoten/Bacteria Interaktionen finden 
# man erhält 1 Liste pro Eukaryot/Bacteria
inter_list_eat <- list("eatenBy","preyedUponBy")

lapplyfunction_species_eat <- function(x) {
  w_species_e <- get_interactions_by_taxa(sourcetaxon = x, 
                                          interactiontype = inter_list_eat)
  if (nrow(w_species_e)>0) w_species_e$qtax_species <- x
  data.frame(w_species_e)
}

b_species_e <- lapply(fd_species[,"species"],lapplyfunction_species_eat)

eat_species <- b_species_e

# Leere Listen (= keine Interaktionen) entfernen
eat_species2 <- eat_species[unlist(lapply(eat_species,nrow)) > 0]

fd3_species_e <- lapply(eat_species2, function(x) { 
  id <- getId(x$target_taxon_name)
  newdf2_species_e <- cbind(x,id)
  newdf2_species_e[!is.na(newdf2_species_e$id),]
})

p_species_e <- length(fd3_species_e)

mylist_species_e <- list()
z_species_e <- for (n in 1:p_species_e) {
  t_species_e <- cbind.data.frame(fd3_species_e[[n]], 
                                  data.frame(getTaxonomy(fd3_species_e[[n]]$id)))
  mylist_species_e[[n]] <- t_species_e
}

s_globi_eat <- unique(do.call(rbind,mylist_species_e))

saveRDS(s_globi_eat,"globi_species_eat.Rds")
