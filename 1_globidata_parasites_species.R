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

desiredtaxa = c("superkingdom","kingdom","superphylum","phylum",
                "subphylum", "superclass","class","subclass", "superorder",
                "order","suborder","family","genus","species")

fd3_species <- do.call(rbind, fd2_species)

fd3_species %>%
  data.frame(getTaxonomy(fd3_species$id_species,desiredTaxa = desiredtaxa)) %>%
  rename(globi_superkingdom = "superkingdom",
         globi_kingdom = "kingdom",
         globi_superphylum = "superphylum",
         globi_phylum = "phylum",
         globi_subphylum="subphylum",
         globi_superclass="superclass",
         globi_class="class",
         globi_subclass="subclass",
         globi_superorder="superorder",
         globi_order="order",
         globi_suborder="suborder",
         globi_family="family",
         globi_genus="genus",
         globi_species="species") %>%
  data.frame(getTaxonomy(getId(fd3_species$qtax_species),
                         desiredTaxa = desiredtaxa)) %>%
  rename(qtax_superkingdom = "superkingdom",
         qtax_kingdom = "kingdom",
         qtax_superphylum = "superphylum",
         qtax_phylum = "phylum",
         qtax_subphylum="subphylum",
         qtax_superclass="superclass",
         qtax_class="class",
         qtax_subclass="subclass",
         qtax_superorder="superorder",
         qtax_order="order",
         qtax_suborder="suborder",
         qtax_family="family",
         qtax_genus="genus") %>%
  select(qtax_superkingdom,
         qtax_kingdom,
         qtax_superphylum,
         qtax_phylum,
         qtax_subphylum,
         qtax_superclass,
         qtax_class,
         qtax_subclass,
         qtax_superorder,
         qtax_order,
         qtax_suborder,
         qtax_family,
         qtax_genus,
         qtax_species,
         globi_superkingdom,
         globi_kingdom,
         globi_superphylum,
         globi_phylum,
         globi_subphylum,
         globi_superclass,
         globi_class,
         globi_subclass,
         globi_superorder,
         globi_order,
         globi_suborder,
         globi_family,
         globi_genus,
         globi_species) -> globi_species_data
         
saveRDS(globi_species_data,"globi_parasites_species.Rds")
