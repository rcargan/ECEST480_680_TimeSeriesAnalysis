

library(dada2);
library(ShortRead);
library(ggplot2);

path <- "/mnt/HA/groups/rosenclassGrp/TimeSeriesGroup/fastq_data"
fastqs <- list.files(path)
fns <- sort(fastqs)
fns <- file.path(path,fns)

#p127 <- plotQualityProfile(fns[[127]])
#ggsave("./images/p127.pdf",p127)

sample.names <- sapply(strsplit(fastqs,"[.]"),'[',1)

filt_path <- file.path(path,"../filtered_data")

for(i in seq_along(fastqs[1:2])) {
  fastq <- fastqs[[i]]
  fastqFilter(file.path(path,fastq), file.path(filt_path, fastq),
              trimLeft=10, truncLen=65,
              maxEE=1, truncQ=10, maxN=0, rm.phix=TRUE, # maybe change maxEE later? 
              compress=TRUE, verbose=TRUE) 
}
