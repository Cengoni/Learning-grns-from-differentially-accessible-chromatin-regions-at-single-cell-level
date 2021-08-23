### Peaks comparison ###
# This script was used to compare peaks stats and feature distribution of macs2 and cellranger peaks
# (3.1.1 Peak calling: comparison between MACS2 and CellRanger)

#BiocManager::install("GenomicRanges") 
#BiocManager::install("GenomicFeatures")
#BiocManager::install("ChIPseeker")
library(GenomicRanges)
library(GenomicFeatures)
library(ChIPseeker)

setwd("/home/claudia/Desktop/UGent/Thesis/peaks/comparison/")
rep1_cell <- read.table("cellranger_rep1.bed")
rep2_cell <- read.table("cellranger_rep2.bed")
rep1_macs <- read.table("rep1_macs.bed")
rep2_macs <- read.table("rep2_macs.bed")

rep1_cell$V4 <- rep1_cell$V3 - rep1_cell$V2
rep1_macs$V4 <- rep1_macs$V3 - rep1_macs$V2
rep2_cell$V4 <- rep2_cell$V3 - rep2_cell$V2
rep2_macs$V4 <- rep2_macs$V3 - rep2_macs$V2

hist(rep1_cell$V4)
hist(rep2_cell$V4)
print(mean(rep1_cell$V4))
print(mean(rep2_cell$V4))
print(median(rep1_cell$V4))
print(median(rep2_cell$V4))

hist(rep1_macs$V4)
hist(rep2_macs$V4)
print(mean(rep1_macs$V4))
print(mean(rep2_macs$V4))
print(median(rep1_macs$V4))
print(median(rep2_macs$V4))

print(var(rep1_cell$V4))
print(var(rep2_cell$V4))
print(var(rep1_macs$V4))
print(var(rep2_macs$V4))
 
# make txdb annotation
txdb <- makeTxDbFromGFF(file="~/Desktop/UGent/Thesis/enrichment/Arabidopsis_thaliana.TAIR10.48.filtered.atac.all.parsed.gtf",
                        format="gtf", dataSource="Ensembl", organism="Arabidopsis thaliana")

# annotate peaks feature distribution
peakAnno_rep1_cellranger <- annotatePeak("cellranger_rep1.bed", tssRegion = c(-5000,0), TxDb=txdb)
peakAnno_rep1_macs <- annotatePeak("rep1_macs.bed",  tssRegion = c(-5000,0), TxDb=txdb)
peakAnno_rep2_cellranger <- annotatePeak("cellranger_rep2.bed",  tssRegion = c(-5000,0), TxDb=txdb)
peakAnno_rep2_macs <- annotatePeak("rep2_macs.bed", tssRegion = c(-5000,0),TxDb=txdb)
 
# plot
peakAnnolist <-  list(cellRanger = peakAnno_rep2_cellranger, MACS2_extended=peakAnno_rep2_macs) 
plotAnnoBar(peakAnnolist)

peakAnnolist <- list(cellRanger = peakAnno_rep1_cellranger, MACS2_extended=peakAnno_rep1_macs) 
plotAnnoBar(peakAnnolist,  title = "Feature Distribution comparison - Replicate 1",)
