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
inter_list <- list("parasiteOf","pathogenOf","endoparasiteOf","ectoparasiteOf")
lapplyfunction <- function(x) {
  w <- get_interactions_by_taxa(sourcetaxon = x, interactiontype = inter_list)
  if (nrow(w)>0) w$qtax <- x
  data.frame(w)
  }

a <- lapply(fd[,"genus"],lapplyfunction)

# Leere Listen (= keine Interaktionen) entfernen
b <- a[unlist(lapply(a,nrow)) > 0]

fd2 <- lapply(b, function(x) { 
  id <- getId(x$target_taxon_name)
  newdf <- cbind(x,id)
  newdf[!is.na(newdf$id),]
  })

p <- length(fd2)

mylist <- list ()
z <- for (n in 1:p) {
  t <- cbind.data.frame(fd2[[n]], data.frame(getTaxonomy(fd2[[n]]$id)))
  mylist[[n]] <- t
}

tidyfox <- do.call(rbind,mylist)

saveRDS(tidyfox,"foxglobi.Rds")
