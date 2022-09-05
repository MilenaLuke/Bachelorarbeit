library(rglobi)
library(dplyr)
library(taxonomizr)
library(data.table)

foxdata <- read.csv(file="ForMilena.csv")

# Doppelte + NA + unidentified/uncultured in species entfernen
unique_foxdata <- unique(foxdata)
fdata_species <- unique_foxdata[!is.na(unique_foxdata$species),]
fd_species <- fdata_species %>% filter(!grepl('unidentified|uncultured', 
                                              species))

# Für jeden Eukaryoten/Bacteria Interaktionen finden 
# man erhält 1 Liste pro Eukaryot/Bacteria
inter_list <- list("parasiteOf","pathogenOf","endoparasiteOf","ectoparasiteOf")
lapplyfunction_species <- function(x) {
  w_species <- get_interactions_by_taxa(sourcetaxon = x, 
                                        interactiontype = inter_list)
  if (nrow(w_species)>0) w_species$qtax_species <- x
  data.frame(w_species)
}

a_species <- lapply(fd_species[,"species"],lapplyfunction_species)

# Leere Listen (= keine Interaktionen) entfernen
b_species <- a_species[unlist(lapply(a_species,nrow)) > 0]

fd2_species <- lapply(b_species, function(x) { 
  id_species <- getId(x$target_taxon_name)
  newdf_species <- cbind(x,id_species)
  newdf_species[!is.na(newdf_species$id_species),]
})

p_species <- length(fd2_species)

mylist_species <- list ()
z_species <- for (n in 1:p_species) {
  t_species <- cbind.data.frame(fd2_species[[n]], 
                                data.frame(getTaxonomy
                                           (fd2_species[[n]]$id_species)))
  mylist_species[[n]] <- t_species
}

globi_species <- do.call(rbind,mylist_species)

saveRDS(globi_species,"globi_parasites_species.Rds")
