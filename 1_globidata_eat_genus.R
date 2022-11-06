library(rglobi)
library(dplyr)
library(taxonomizr)
library(data.table)

foxdata <- read.csv(file="ForMilena.csv")

# Doppelte + NA + unidentified/uncultured entfernen
unique_foxdata <- unique(foxdata)
fdata <- unique_foxdata[!is.na(unique_foxdata$genus),]
fd <- fdata %>% filter(!grepl('unidentified|uncultured', genus))

# Für jeden Eukaryoten/Bacteria Interaktionen finden 
# man erhält 1 Liste pro Eukaryot/Bacteria
inter_list_eat <- list("eatenBy","preyedUponBy")

lapplyfunction2 <- function(x) {
  w <- get_interactions_by_taxa(sourcetaxon = x, interactiontype = inter_list_eat)
  if (nrow(w)>0) w$qtax2 <- x
  data.frame(w)
}

b <- lapply(fd[,"genus"],lapplyfunction2)

eat <- b

# Leere Listen (= keine Interaktionen) entfernen
eat_2 <- eat[unlist(lapply(eat,nrow)) > 0]

fd3 <- lapply(eat_2, function(x) { 
  id <- getId(x$target_taxon_name)
  newdf2 <- cbind(x,id)
  newdf2[!is.na(newdf2$id),]
})

p2 <- length(fd3)

mylist2 <- list()
z2 <- for (n in 1:p2) {
  t2 <- cbind.data.frame(fd3[[n]], data.frame(getTaxonomy(fd3[[n]]$id)))
  mylist2[[n]] <- t2
}

tidyfox_eat <- unique(do.call(rbind,mylist2))

saveRDS(tidyfox_eat,"foxglobieat.Rds")

