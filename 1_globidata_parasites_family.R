library(rglobi)
library(dplyr)
library(taxonomizr)
library(data.table)

foxdata <- read.csv(file="ForMilena.csv")

# Doppelte + NA + unidentified/uncultured in species entfernen
unique_foxdata <- unique(foxdata)
fdata_f <- unique_foxdata[!is.na(unique_foxdata$family),]
fd_f <- fdata_f %>% filter(!grepl('unidentified|uncultured', 
                                              family))

# Für jeden Eukaryoten/Bacteria Interaktionen finden 
# man erhält 1 Liste pro Eukaryot/Bacteria
inter_list <- list("parasiteOf","pathogenOf","endoparasiteOf","ectoparasiteOf")
lapplyfunction_f <- function(x) {
  w_family <- get_interactions_by_taxa(sourcetaxon = x, 
                                        interactiontype = inter_list)
  if (nrow(w_family)>0) w_family$qtax_family <- x
  data.frame(w_family)
}

a_family <- lapply(fd_f[,"family"],lapplyfunction_f)

# Leere Listen (= keine Interaktionen) entfernen
b_family <- a_family[unlist(lapply(a_family,nrow)) > 0]

fd2_fam <- lapply(b_family, function(x) { 
  id_fam <- getId(x$target_taxon_name)
  newdf_fam <- cbind(x,id_fam)
  newdf_fam[!is.na(newdf_fam$id_fam),]
})

p_fam <- length(fd2_fam)

mylist_fam <- list ()
z_fam <- for (n in 1:p_fam) {
  t_fam <- cbind.data.frame(fd2_fam[[n]], 
                                data.frame(getTaxonomy
                                           (fd2_fam[[n]]$id_fam)))
  mylist_fam[[n]] <- t_fam
}

globi_fam <- do.call(rbind,mylist_fam)

saveRDS(globi_fam,"globi_parasites_family.Rds")
