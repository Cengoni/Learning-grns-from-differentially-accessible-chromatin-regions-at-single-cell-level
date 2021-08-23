### Utility script
# Aim: extract TFs coordinates - 1 kb for gene-body openness filtering
library(ballgown)

### Read in annotation ###
annotation <- gffReadGR("/home/claudia/Desktop/UGent/Thesis/enrichment/Arabidopsis_thaliana.TAIR10.48.filtered.atac.all.parsed.gtf")
 
tfs <- read.csv("/home/claudia/Desktop/UGent/Thesis/enrichment/TF2fam2mot_n.txt", sep='\t', header=FALSE)
tfs <- unique(tfs$V1)

subset_annotation <- annotation[annotation$gene_id %in% tfs & annotation$type %in% 'transcript']
 

names <- subset_annotation$gene_id 
start_ <- subset_annotation@ranges@start - 1000
end_ <- subset_annotation@ranges@start + subset_annotation@ranges@width - 1
chrom_lengths <- subset_annotation@seqnames@lengths

chrom = c()
for (i in c(1:length(chrom_lengths))){
    chrom = c(chrom,rep(i, chrom_lengths[i]))
}
 

length(names)
length(start_)
length(end_)
length(chrom)

x<-data.frame(names,chrom, start_,end_)
x<-x[x$chrom %in% c(1,2,3,4,5),] 

write.csv(x, "/home/claudia/Desktop/UGent/Thesis/enrichment/tfs_pos_1kb.csv", row.names = FALSE, quote = FALSE)
