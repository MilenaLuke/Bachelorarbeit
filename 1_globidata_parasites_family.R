library(rglobi)
library(dplyr)
library(taxonomizr)
library(data.table)

foxdata <- read.csv(file="ForMilena.csv")

# Doppelte + NA + unidentified/uncultured entfernen
unique_foxdata <- unique(foxdata)
fdata_family <- unique_foxdata[!is.na(unique_foxdata$family),]
fd_family <- fdata_family %>% filter(!grepl('unidentified|uncultured', family))

# Für jeden Eukaryoten/Bacteria Interaktionen finden 
# man erhält 1 Liste pro Eukaryot/Bacteria
inter_list <- list("parasiteOf","pathogenOf","endoparasiteOf","ectoparasiteOf")
lapplyfunction_family <- function(x) {
  w_family <- get_interactions_by_taxa(sourcetaxon = x, 
                                      interactiontype = inter_list)
  if (nrow(w_family)>0) w_family$qtax <- x
  data.frame(w_family)
}

a_family <- lapply(fd_family[,"family"],lapplyfunction_family)

# Leere Listen (= keine Interaktionen) entfernen
b_family <- a_family[unlist(lapply(a_family,nrow)) > 0]

fd2_family <- lapply(b_family, function(x) { 
  id_family <- getId(x$target_taxon_name)
  newdf_family <- cbind(x,id_family)
  newdf_family[!is.na(newdf_family$id_family),]
})

p_family <- length(fd2_family)

desiredtaxa = c("superkingdom","kingdom","superphylum","phylum",
                "subphylum", "superclass","class","subclass", "superorder",
                "order","suborder","family","genus","species")

fd3_family <- do.call(rbind, fd2_family)

fd3_family %>%
  data.frame(getTaxonomy(fd3_family$id_family,desiredTaxa = desiredtaxa)) %>%
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
  data.frame(getTaxonomy(getId(fd3_family$qtax),
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
         qtax_family="family") %>%
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
         globi_species) -> globi_family_data

saveRDS(globi_family_data,"globi_parasites_family.Rds")
