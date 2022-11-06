library(rglobi)
library(dplyr)
library(taxonomizr)
library(data.table)

foxdata <- read.csv(file="ForMilena.csv")

# Doppelte + NA + unidentified/uncultured entfernen
unique_foxdata <- unique(foxdata)
fdata_genus <- unique_foxdata[!is.na(unique_foxdata$genus),]
fd_genus <- fdata %>% filter(!grepl('unidentified|uncultured', genus))

# Für jeden Eukaryoten/Bacteria Interaktionen finden 
# man erhält 1 Liste pro Eukaryot/Bacteria
inter_list <- list("parasiteOf","pathogenOf","endoparasiteOf","ectoparasiteOf")
lapplyfunction_genus <- function(x) {
  w_genus <- get_interactions_by_taxa(sourcetaxon = x, 
                                      interactiontype = inter_list)
  if (nrow(w_genus)>0) w_genus$qtax <- x
  data.frame(w_genus)
  }

a_genus <- lapply(fd_genus[,"genus"],lapplyfunction_genus)

# Leere Listen (= keine Interaktionen) entfernen
b_genus <- a_genus[unlist(lapply(a_genus,nrow)) > 0]

fd2_genus <- lapply(b_genus, function(x) { 
  id_genus <- getId(x$target_taxon_name)
  newdf_genus <- cbind(x,id_genus)
  newdf_genus[!is.na(newdf_genus$id_genus),]
  })

p_genus <- length(fd2_genus)

desiredtaxa = c("superkingdom","kingdom","superphylum","phylum",
                "subphylum", "superclass","class","subclass", "superorder",
                "order","suborder","family","genus","species")

fd3_genus <- do.call(rbind, fd2_genus)

fd3_genus %>%
  data.frame(getTaxonomy(fd3_genus$id_genus,desiredTaxa = desiredtaxa)) %>%
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
  data.frame(getTaxonomy(getId(fd3_genus$qtax),
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
         globi_species) -> globi_genus_data

saveRDS(globi_genus_data,"globi_parasites_genus.Rds")
