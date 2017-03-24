# Programs to load #

library(stringr)
library(phyloseq)
library(biom)
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(vegan)
library(viridis)
library(proxy)
library(plotly)

# Set Working Directory or Set Data Directory #

data_dir <- '~/bryan'

# Create Phylseq Object #

OTU <- readRDS(file.path(data_dir,'seqtab_final.rds'))
TAX <- read.table(file.path(data_dir,'taxatable.csv'),stringsAsFactors=FALSE,sep=',',header=TRUE)
rownames(TAX) <- TAX[,1]
TAX <- TAX[,-1]
META <- read.delim(file.path(data_dir,'TimeSeries_Metadata.txt'))

OTU <- OTU[rowSums(OTU) > 1000,]


SAMPLEIDS <- intersect(rownames(OTU),META$SampleID)

META <- META %>% filter(SampleID %in% SAMPLEIDS)
OTU <- OTU[SAMPLEIDS,]

OTU <- OTU[,colSums(OTU)>0]

OTUIDS <- intersect(rownames(TAX)[!is.na(TAX$Phylum)],names(which(colSums(OTU)>0)))
OTU <- OTU[,OTUIDS]

TAX <- TAX[colnames(OTU),]


META <- as.data.frame(META)
rownames(META) <- as.character(META$SampleID)
META <- sample_data(META[rownames(OTU),])
TAX <- tax_table(as.matrix(TAX[colnames(OTU),]))
OTU <- otu_table(OTU,taxa_are_rows=FALSE)

PS <- phyloseq(OTU,META,TAX)
PS <- prune_samples(sample_sums(PS)>0,PS)

# Jensen Shannon #

jsd <- function(p,q){
  m <- (p+q)/2
  sqrt(.5*(sum(p*log(p/m)) + sum(q*log(q/m))))
}

# Heatmap DonorA Stool #

TAX2 <- as(tax_table(PS),'matrix')
META2 <- as(sample_data(PS),'data.frame')
OTU1 <- as(otu_table(PS),'matrix')
OTU2 <- OTU1 + 1
OTURA <- OTU2/rowSums(OTU2)

DON1 <- OTURA[META2 %>% 
                filter(description_s %in% 'DonorA Stool') %>% 
                select(SampleID) %>%
                unlist(),]

DIST1 <- as.matrix(proxy::dist(DON1,jsd))

DF1 <- data.frame(DIST1,SampleID=rownames(DIST1)) %>%
  gather(SampleID2,distance,-SampleID) %>%
  left_join(META2 %>% select(SampleID,description_s,collection_day_s,description_s),by='SampleID') %>%
  mutate(collection_day_s=as.integer(as.character.factor(collection_day_s))) %>%
  arrange(collection_day_s,SampleID,SampleID2) %>%
  mutate(SampleID1=factor(SampleID,levels=unique(SampleID),ordered=TRUE),
         SampleID2=factor(SampleID2,levels=unique(SampleID2),ordered=TRUE))

ss_ids <- sample(DF1$SampleID,15)

DF1 %>%
  # filter(SampleID1 %in% ss_ids,SampleID2 %in% ss_ids) %>%
  arrange(collection_day_s,SampleID,SampleID2) %>%
  mutate(SampleID1=factor(SampleID,levels=unique(SampleID),ordered=TRUE),
         SampleID2=factor(SampleID2,levels=unique(SampleID),ordered=TRUE)) %>%
  ggplot(aes(SampleID1,SampleID2,fill=distance)) +
  geom_raster() +
  geom_rug(aes(color=collection_day_s),sides='b',size=1) +
  viridis::scale_fill_viridis('JSD',direction=-1) +
  scale_color_distiller(type='div',palette=5) +
  labs(x='SampleID',y='SampleID',color='Day') +
  theme(aspect.ratio=1,
        axis.text=element_blank(),
        axis.ticks=element_blank()) +
  ggtitle(paste0(unique(as.character(DF1$description_s))))

# Heatmap DonorB Stool #

TAX2 <- as(tax_table(PS),'matrix')
META2 <- as(sample_data(PS),'data.frame')
OTU1 <- as(otu_table(PS),'matrix')
OTU2 <- OTU1 + 1
OTURA <- OTU2/rowSums(OTU2)

DON2 <- OTURA[META2 %>% 
                filter(description_s %in% 'DonorB Stool') %>% 
                select(SampleID) %>%
                unlist(),]

DIST1 <- as.matrix(proxy::dist(DON2,jsd))

DF1 <- data.frame(DIST1,SampleID=rownames(DIST1)) %>%
  gather(SampleID2,distance,-SampleID) %>%
  left_join(META2 %>% select(SampleID,description_s,collection_day_s,description_s),by='SampleID') %>%
  mutate(collection_day_s=as.integer(as.character.factor(collection_day_s))) %>%
  arrange(collection_day_s,SampleID,SampleID2) %>%
  mutate(SampleID1=factor(SampleID,levels=unique(SampleID),ordered=TRUE),
         SampleID2=factor(SampleID2,levels=unique(SampleID2),ordered=TRUE))

ss_ids <- sample(DF1$SampleID,15)

DF1 %>%
  # filter(SampleID1 %in% ss_ids,SampleID2 %in% ss_ids) %>%
  arrange(collection_day_s,SampleID,SampleID2) %>%
  mutate(SampleID1=factor(SampleID,levels=unique(SampleID),ordered=TRUE),
         SampleID2=factor(SampleID2,levels=unique(SampleID),ordered=TRUE)) %>%
  ggplot(aes(SampleID1,SampleID2,fill=distance)) +
  geom_raster() +
  geom_rug(aes(color=collection_day_s),sides='b',size=1) +
  viridis::scale_fill_viridis('JSD',direction=-1) +
  scale_color_distiller(type='div',palette=5) +
  labs(x='SampleID',y='SampleID',color='Day') +
  theme(aspect.ratio=1,
        axis.text=element_blank(),
        axis.ticks=element_blank()) +
  ggtitle(paste0(unique(as.character(DF1$description_s))))

# Heatmap DonorA Saliva #

TAX2 <- as(tax_table(PS),'matrix')
META2 <- as(sample_data(PS),'data.frame')
OTU1 <- as(otu_table(PS),'matrix')
OTU2 <- OTU1 + 1
OTURA <- OTU2/rowSums(OTU2)

DON3 <- OTURA[META2 %>% 
                filter(description_s %in% 'DonorA Saliva') %>% 
                select(SampleID) %>%
                unlist(),]

DIST1 <- as.matrix(proxy::dist(DON3,jsd))

DF1 <- data.frame(DIST1,SampleID=rownames(DIST1)) %>%
  gather(SampleID2,distance,-SampleID) %>%
  left_join(META2 %>% select(SampleID,description_s,collection_day_s,description_s),by='SampleID') %>%
  mutate(collection_day_s=as.integer(as.character.factor(collection_day_s))) %>%
  arrange(collection_day_s,SampleID,SampleID2) %>%
  mutate(SampleID1=factor(SampleID,levels=unique(SampleID),ordered=TRUE),
         SampleID2=factor(SampleID2,levels=unique(SampleID2),ordered=TRUE))

ss_ids <- sample(DF1$SampleID,15)

DF1 %>%
  # filter(SampleID1 %in% ss_ids,SampleID2 %in% ss_ids) %>%
  arrange(collection_day_s,SampleID,SampleID2) %>%
  mutate(SampleID1=factor(SampleID,levels=unique(SampleID),ordered=TRUE),
         SampleID2=factor(SampleID2,levels=unique(SampleID),ordered=TRUE)) %>%
  ggplot(aes(SampleID1,SampleID2,fill=distance)) +
  geom_raster() +
  geom_rug(aes(color=collection_day_s),sides='b',size=1) +
  viridis::scale_fill_viridis('JSD',direction=-1) +
  scale_color_distiller(type='div',palette=5) +
  labs(x='SampleID',y='SampleID',color='Day') +
  theme(aspect.ratio=1,
        axis.text=element_blank(),
        axis.ticks=element_blank()) +
  ggtitle(paste0(unique(as.character(DF1$description_s))))

#Relative Abundance vs Time #

time <- length(levels(META2$collection_day_s))
P5 <- data.frame(OTU1,SampleID=rownames(OTU1)) %>%
  gather(sequence,count,-SampleID) %>%
  left_join(META2 %>% select(SampleID,description_s,collection_day_s,description_s),by='SampleID') %>%
  mutate(collection_day_s=as.integer(as.character.factor(collection_day_s))) %>%
  left_join(data.frame(TAX2,sequence=rownames(TAX2)),by='sequence') %>%
  dplyr::mutate(Taxon=Order, # change this
                Taxon=ifelse(is.na(Taxon),'Other',as.character(Taxon))) %>%
  group_by(Taxon,description_s,collection_day_s) %>%
  dplyr::summarize(count=sum(count)) %>%
  ungroup() %>%
  group_by(Taxon) %>%
  mutate(Taxon_sum=sum(count)) %>%
  ungroup() %>%
  mutate(Taxon_rank=dense_rank(desc(Taxon_sum)),
         Taxon=ifelse(Taxon_rank > 10,'Other',Taxon)) %>% # change this
  group_by(Taxon,description_s,collection_day_s) %>%
  summarize(count=sum(count)) %>%
  ungroup() %>%
  arrange(desc(count)) %>%
  mutate(Taxon=factor(Taxon,levels=unique(Taxon))) %>%
  filter(!(Taxon %in% 'Other')) %>%
  ggplot() +
  geom_bar(aes(collection_day_s,count,fill=Taxon,color=Taxon),stat='identity',position='fill') +
  facet_grid(description_s~.) +
  theme_dark() + 
  theme(axis.text.x=element_text(angle=90,hjust=1),
        legend.position='bottom') +
  scale_fill_viridis(discrete=TRUE) +
  scale_color_viridis(discrete=TRUE) +
  scale_x_continuous(breaks=seq(0,time,by=10)) +
  labs(y='Relative Abundance',x='Day',fill='',color='')

ggplotly(P5) %>% layout(legend=list(orientation='h',x=0.0,y=100,size=1), 
                        margin=list(l=100,r=100,b=70,t=70,pad=5))
