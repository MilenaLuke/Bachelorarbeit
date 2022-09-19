library(rglobi)
library(dplyr)
library(taxonomizr)
library(data.table)

foxdata <- read.csv(file="ForMilena.csv")

# Doppelte + NA + unidentified/uncultured entfernen
unique_foxdata <- unique(foxdata)
fdata_f <- unique_foxdata[!is.na(unique_foxdata$family),]
fd_f <- fdata %>% filter(!grepl('unidentified|uncultured', family))

# Für jeden Eukaryoten/Bacteria Interaktionen finden 
# man erhält 1 Liste pro Eukaryot/Bacteria
inter_list_eat <- list("eatenBy","preyedUponBy")

lapplyfunction4 <- function(x) {
  w_family_e <- get_interactions_by_taxa(sourcetaxon = x, 
                                  interactiontype = inter_list_eat)
  if (nrow(w_family_e)>0) w_family_e$qtax_family <- x
  data.frame(w_family_e)
}

b_family_e <- lapply(fd_f[,"family"],lapplyfunction4)

# Leere Listen (= keine Interaktionen) entfernen

family_eat <- b_family_e[unlist(lapply(b_family_e,nrow)) > 0]

fd3_family_eat <- lapply(family_eat, function(x) { 
  id <- getId(x$target_taxon_name)
  newdf2_family_e <- cbind(x,id)
  newdf2_family_e[!is.na(newdf2_family_e$id),]
})

p2_family_eat <- length(fd3_family_eat)

mylist2_family_eat <- list()
z2_family_eat <- for (n in 1:p2_family_eat) {
  t2_family_eat <- cbind.data.frame(fd3_family_eat[[n]], 
                                    data.frame(
                                      getTaxonomy(fd3_family_eat[[n]]$id)))
  mylist2_family_eat[[n]] <- t2_family_eat
}

tidyfox_eat_family <- unique(do.call(rbind,mylist2_family_eat))

saveRDS(tidyfox_eat_family,"globi_family_eat.Rds")

