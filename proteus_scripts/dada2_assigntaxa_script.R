

library(dada2)
library(gridExtra)
library(Biostrings)
library(ggplot2)

# Merge sequence tables #

st1 <- readRDS("/mnt/HA/groups/rosenclassGrp/TimeSeriesGroup/sequence_tables/sequence_table_data_test14/seqtab.rds")
st2 <- readRDS("/mnt/HA/groups/rosenclassGrp/TimeSeriesGroup/sequence_tables/sequence_table_data_test15/seqtab.rds")
st3 <- readRDS("/mnt/HA/groups/rosenclassGrp/TimeSeriesGroup/sequence_tables/sequence_table_data_test16/seqtab.rds")
st4 <- readRDS("/mnt/HA/groups/rosenclassGrp/TimeSeriesGroup/sequence_tables/sequence_table_data_test17/seqtab.rds")
st5 <- readRDS("/mnt/HA/groups/rosenclassGrp/TimeSeriesGroup/sequence_tables/sequence_table_data_test18/seqtab.rds")
st6 <- readRDS("/mnt/HA/groups/rosenclassGrp/TimeSeriesGroup/sequence_tables/sequence_table_data_test19/seqtab.rds")
st7 <- readRDS("/mnt/HA/groups/rosenclassGrp/TimeSeriesGroup/sequence_tables/sequence_table_data_test20/seqtab.rds")
st8 <- readRDS("/mnt/HA/groups/rosenclassGrp/TimeSeriesGroup/sequence_tables/sequence_table_data_test21/seqtab.rds")

st <- mergeSequenceTables(st1, st2, st3, st4, st5, st6, st7, st8)
saveRDS(st, "/mnt/HA/groups/rosenclassGrp/TimeSeriesGroup/sequence_tables/seqtab_final.rds")

### Assign taxonomy (make sure the working directory is set)
seqtab <- readRDS("/mnt/HA/groups/rosenclassGrp/TimeSeriesGroup/sequence_tables/seqtab_final.rds")
seqs <- colnames(seqtab)
taxa <- assignTaxonomy(seqs, "/mnt/HA/groups/rosenclassGrp/TimeSeriesGroup/silva/silva_nr_v123_train_set.fa")

colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
unname(head(taxa))
dim(taxa)
write.csv(taxa, "/mnt/HA/groups/rosenclassGrp/TimeSeriesGroup/taxatable.csv")

# Assign species where possible
taxa.plus <- addSpecies(taxa, "/mnt/HA/groups/rosenclassGrp/TimeSeriesGroup/silva/silva_species_assignment_v123.fa", verbose=TRUE)
unname(head(taxa.plus))

write.csv(taxa, "/mnt/HA/groups/rosenclassGrp/TimeSeriesGroup/taxatable.csv")
