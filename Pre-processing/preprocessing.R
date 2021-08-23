# load packages
library(tidyr)
library(dplyr)
library(data.table)
library(ballgown)
library(Signac)
library(Seurat)
library(GenomeInfoDb) 
library(ggplot2)
library(patchwork)
library(ChIPpeakAnno)
library(GenomicFeatures) 
library(Matrix)
library(harmony) 
library(GenomicRanges) 
library(ChIPseeker)
 

#load(file="/ngsprojects/scgrn/data_archive/clmen/signac_scripts/farmer.RData")
#save.image(file="/ngsprojects/scgrn/data_archive/clmen/signac_scripts/farmer.RData")

### Read in annotation file ###
annotation <- gffReadGR("/ngsprojects/scgrn/data_archive/clmen/Ath/Arabidopsis_thaliana.TAIR10.48.filtered.atac.all.parsed.gtf")
seqlevelsStyle(annotation) <- 'Ensembl'
genome(annotation) <- 'TAIR10' 

# make txdb annotation
txdb <- makeTxDbFromGFF(file="/ngsprojects/scgrn/data_archive/clmen/Ath/Arabidopsis_thaliana.TAIR10.48.filtered.atac.all.parsed.gtf",
                        format="gtf", dataSource="Ensembl", organism="Arabidopsis thaliana")
ge <- genes(txdb, columns=c("tx_name", "gene_id", "tx_type")) 


### Replicate 1 ###

# read peak-by-cell matrix 
mex_dir_path <- "/ngsprojects/scgrn/data_archive/clmen/farmer_etal/cellranger/rep1_summit_py/outs/filtered_peak_bc_matrix/"
mtx_path <- paste(mex_dir_path, "matrix.mtx", sep = '/')
feature_path <- "/ngsprojects/scgrn/data_archive/clmen/farmer_etal/macs2/rep1_peaks.bed" #paste(mex_dir_path, "peaks.bed", sep = '/')
barcode_path <- paste(mex_dir_path, "barcodes.tsv", sep = '/')

m <- readMM(mtx_path)
barcodes <- fread(barcode_path, header=F)[[1]]
features <- fread(feature_path, header=F) %>% tidyr::unite(feature, sep='-')
colnames(m) <- barcodes
rownames(m) <- features$feature
m <- m*1

# feature ranges
df_features<-fread(feature_path, header=F)
colnames(df_features)<-c("chr","start","stop")
ranges_1<-makeGRangesFromDataFrame(df_features)

# metadata
metadata <- read.csv(
  file = "/ngsprojects/scgrn/data_archive/clmen/farmer_etal/cellranger/rep1_summit_py/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)

# create fragment file
frag_path<-("/ngsprojects/scgrn/data_archive/clmen/farmer_etal/cellranger/rep1_summit_py/outs/fragments.tsv.gz") 
fragments <- CreateFragmentObject(
  path = frag_path,
)

# create chromatin assay
chrom_assay <- CreateChromatinAssay(
  counts = m, 
  fragments = fragments,
  annotation = annotation,
  ranges = ranges_1,
  min.features = 1000)

# create seurat object
ath_r1 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

# assign name to assay 
ath_r1[['peaks']]


## QC ##

# add fraction of reads in peaks
ath_r1$pct_reads_in_peaks <- ath_r1$peak_region_fragments / ath_r1$passed_filters * 100 

VlnPlot(
  object = ath_r1,
  features = c('pct_reads_in_peaks', 'peak_region_fragments'),
  pt.size = 0.1,
  ncol = 2
)

ath_r1 <- ath_r1[,ath_r1$pct_reads_in_peaks>15]
ath_r1 <- ath_r1[,ath_r1$peak_region_fragments<2e4]

#compute nucleosome signal
ath_r1 <- NucleosomeSignal(object = ath_r1)

# compute TSS enrichment score per cell
ath_r1 <- TSSEnrichment(object = ath_r1, fast = FALSE)
 

VlnPlot(
  object = ath_r1,
  features = c('TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 2
)

ath_r1$high.tss <- ifelse(ath_r1$TSS.enrichment > 1.5, 'High', 'Low')
TSSPlot(ath_r1, group.by = 'high.tss') + NoLegend()
TSSPlot(ath_r1) + NoLegend()


ath_r1 <- subset(x = ath_r1,
                 subset = TSS.enrichment > 1.5 & 
                 nucleosome_signal < 5)#mean(ath_r1$nucleosome_signal)+2*sd(ath_r1$nucleosome_signal))
 
ath_r1

### Replicate 2 ###

# read peak-cell matrix 
mex_dir_path <- "/ngsprojects/scgrn/data_archive/clmen/farmer_etal/cellranger/rep2_summit_py/outs/filtered_peak_bc_matrix/"
mtx_path <- paste(mex_dir_path, "matrix.mtx", sep = '/')
feature_path <- feature_path <- "/ngsprojects/scgrn/data_archive/clmen/farmer_etal/macs2/rep2_peaks.bed" #paste(mex_dir_path, "peaks.bed", sep = '/') 
barcode_path <- paste(mex_dir_path, "barcodes.tsv", sep = '/')

m <- readMM(mtx_path)
barcodes <- fread(barcode_path, header=F)[[1]]
features <- fread(feature_path, header=F) %>% tidyr::unite(feature, sep='-')
colnames(m) <- barcodes
rownames(m) <- features$feature
m <- m*1

# feature ranges
df_features<-fread(feature_path, header=F)
colnames(df_features)<-c("chr","start","stop")
ranges_2<-makeGRangesFromDataFrame(df_features )

# metadata
metadata <- read.csv(
  file = "/ngsprojects/scgrn/data_archive/clmen/farmer_etal/cellranger/rep2_summit_py/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)

# create fragment file
frag_path<-("/ngsprojects/scgrn/data_archive/clmen/farmer_etal/cellranger/rep2_summit_py/outs/fragments.tsv.gz") 
fragments <- CreateFragmentObject(
  path = frag_path,
)

# create chromatin assay
chrom_assay <- CreateChromatinAssay(
  counts = m, 
  fragments = fragments,
  annotation = annotation,
  ranges = ranges_2,
  min.features = 1000)

# create seurat object
ath_r2 <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

# assign name to assay 
ath_r2[['peaks']]


## QC ##

# add fraction of reads in peaks
ath_r2$pct_reads_in_peaks <- ath_r2$peak_region_fragments / ath_r2$passed_filters * 100 

VlnPlot(
  object = ath_r2,
  features = c('pct_reads_in_peaks', 'peak_region_fragments'),
  pt.size = 0.1,
  ncol = 2
)

ath_r2 <- ath_r2[,ath_r2$pct_reads_in_peaks>15]
ath_r2 <- ath_r2[,ath_r2$peak_region_fragments<2e4]

#compute nucleosome signal
ath_r2 <- NucleosomeSignal(object = ath_r2)

# compute TSS enrichment score per cell
ath_r2 <- TSSEnrichment(object = ath_r2, fast = FALSE)


VlnPlot(
  object = ath_r2,
  features = c('TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 2
)

ath_r2$high.tss <- ifelse(ath_r2$TSS.enrichment > 1.5, 'High', 'Low')
TSSPlot(ath_r2, group.by = 'high.tss') + NoLegend()

#ath_r2$nucleosome_group <- ifelse(ath_r2$nucleosome_signal > 5, 'NS > 5', 'NS < 5')
#FragmentHistogram(object = ath_r2, group.by = 'nucleosome_group')

ath_r2 <- subset(x = ath_r2,
                 subset = TSS.enrichment > 1.5 & 
                 nucleosome_signal < 5)# mean(ath_r2$nucleosome_signal)+2*sd(ath_r2$nucleosome_signal))

### MERGE ###
# Create a unified set of peaks to quantify in each dataset 

ol <- findOverlapsOfPeaks(ranges_1, ranges_2, minoverlap = 20)  

makeVennDiagram(ol,
                fill=c("#009E73", "#F0E442"), 
                col=c("#D55E00", "#0072B2"),  
                cat.col=c("#D55E00", "#0072B2"))

combined.peaks<-ol$peaklist[["ranges_1///ranges_2"]]

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 5000 & peakwidths > 20]
length(combined.peaks)

# peaks location 
peakAnno <- annotatePeak(combined.peaks, tssRegion = c(-5000,0),
                         TxDb=txdb)
plotAnnoPie(peakAnno)

# rep1
# reload fragment file only for cells left from QC
frags_1 <- CreateFragmentObject(path = "/ngsprojects/scgrn/data_archive/clmen/farmer_etal/cellranger/rep1_summit_py/outs/fragments.tsv.gz",cells = colnames(ath_r1))

# quantify peaks from common set in dataset 1
counts_1 <- FeatureMatrix(fragments = frags_1,features = combined.peaks,cells = colnames(ath_r1))
 
# create objects to merge
assay_1 <- CreateChromatinAssay(counts_1, fragments = frags_1, ranges=combined.peaks)
combined_ath_1 <- CreateSeuratObject(assay_1, assay = "ATAC") 
combined_ath_1$dataset <- 'rep1' 

# rep2 
# change names to avoid replicate
new_names <- paste0(substr(colnames(ath_r2),start=1,stop=nchar(colnames(ath_r2))-1),"2")

# update object with renamed fragment file 
frags_renamed_2 <- CreateFragmentObject(path = "/ngsprojects/scgrn/data_archive/clmen/farmer_etal/cellranger/rep2_summit_py/outs/fragments_renamed.tsv.gz",cells = new_names)

# quantify peaks from common set in dataset 2
counts_renamed_2 <- FeatureMatrix(fragments = frags_renamed_2,features = combined.peaks, cells = new_names)

# create objects to merge
assay_renamed_2 <- CreateChromatinAssay(counts_renamed_2, fragments = frags_renamed_2, ranges=combined.peaks)
combined_ath_renamed_2 <- CreateSeuratObject(assay_renamed_2, assay = "ATAC")
combined_ath_renamed_2$dataset <- 'rep2'

# merge datasets
combined <- merge(
  x = combined_ath_1,
  y = combined_ath_renamed_2
)

combined[["ATAC"]] 

### Post-processing ###
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 20)
combined <- RunSVD(combined) 

# Harmony integration 
integrated <- RunHarmony(
  object = combined,
  group.by.vars = 'dataset',
  reduction = 'lsi',
  assay.use = 'ATAC',
  project.dim = FALSE
)

# re-compute the UMAP using corrected LSI embeddings
integrated <- RunUMAP(integrated, dims=1:20, reduction = 'harmony')
DimPlot(integrated, group.by = 'dataset', pt.size = 0.1)  
integrated <- FindNeighbors(object = integrated, reduction = 'harmony', dims = 1:20)
integrated <- FindClusters(object = integrated, verbose = FALSE, algorithm = 3, resolution=1.3)
DimPlot(object = integrated, label = TRUE)# + NoLegend() 
 
# Find DARs
da_peaks <- FindAllMarkers(
  object = integrated, 
  only.pos = TRUE
)

da_peaks<-da_peaks[da_peaks$p_val_adj<0.05,]
t<-table(da_peaks$gene)

# make granges from DARs
regions <- da_peaks$gene
regions.df<-as.data.frame(t(data.frame(sapply(1:length(regions), function(i){strsplit(regions[i],"[-]")}))))
colnames(regions.df)<-c("chr","start","stop") 
regions.gr<-makeGRangesFromDataFrame(regions.df)
u_regions.gr<-unique(regions.gr)

# annotate peaks 
closest_genes <- annoPeaks(
  u_regions.gr,
  genes(txdb),
  bindingType = "fullRange",
  bindingRegion = c(-5000, 1000), 
  select = "bestOne")
 
plotAnnoPie(closest_genes, main='DARs distribution')

# couple DARs and annotations into a dataframe
cg<-GRangesToString(closest_genes, sep = c("-", "-"))
da<-da_peaks[cg,c("p_val_adj","cluster")]
da<-data.frame(cbind(da, closest_genes$feature))
colnames(da)<-c("p_val_adj","cluster","feature") 

# annotate clusters 
wendrich <- read.csv("/ngsprojects/scgrn/data_archive/markers_lit_Wendrich4.txt", fill=1, header=TRUE, sep='\t') 
marker <- unique(wendrich)

nclusters<-length(levels(integrated$seurat_clusters))
identification_top<- data.frame(c(1:nclusters-1), 
                                c(1:nclusters-1)) 
colnames(identification_top)<-c("cluster","type")
 

for (i in c(1:nclusters)){
  subset<-da[da$cluster==i,]
  names<-subset$feature[order(subset$p_val_adj)]  
  names<-names[1:200]
  for (n in names){ 
    if (n %in% marker$gene){  
      print(n)
      print(i)
      print(as.character(marker[marker$gene==n,'type']))
      identification_top$type[identification_top$cluster==i] <- as.character(marker[marker$gene==n,'type'])
      break
    } 
  }
}


identification_top 

# assign annotation to clusters
new.clusters.id <- as.list(identification_top$type)
names(new.clusters.id) <- levels(integrated)
integrated_final<- RenameIdents(integrated, new.clusters.id)
DimPlot(integrated_final, pt.size = 0.1,label=TRUE)

## Other visualization
# create gene activity matrix
annotation_geneid<-annotation
annotation_geneid$gene_name<-annotation_geneid$gene_id
transcripts <- annotation[annotation$type=="transcript"]
transcripts <- transcripts[transcripts$gene_biotype == "protein_coding"] 

Annotation(integrated) <- annotation_geneid
gene.activities <- GeneActivity(integrated, features=transcripts$gene_id, extend.upstream = 5000, extend.downstream=1000, process_n=200) #50G Memory

integrated[['RNA']] <- CreateAssayObject(counts = gene.activities)
integrated <- NormalizeData(
  object = integrated,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(integrated$nCount_RNA)
)

DefaultAssay(integrated) <- 'RNA'
 
DotPlot(integrated, features = wendrich[wendrich$type=='endodermis','gene'] , col.min=0, col.max=3.5,  cols = c("white", "blue")) + RotatedAxis()
DotPlot(integrated, features = wendrich[wendrich$type=='phloem','gene'] , col.min=0, col.max=3.5,  cols = c("white", "blue")) + RotatedAxis()

# add new idenities to clusters (endodermis and phloem) 
identification_top$type[identification_top$type==0] <- 'endodermis'
identification_top$type[identification_top$type==16] <- 'endodermis'
identification_top$type[identification_top$type==17] <- 'endodermis'
identification_top$type[identification_top$type==14] <- 'phloem'
new.clusters.id <- as.list(identification_top$type)
names(new.clusters.id) <- levels(integrated)
integrated_final<- RenameIdents(integrated, new.clusters.id)
DimPlot(integrated_final, pt.size = 0.1, label=TRUE)

### Extract obtained DARs ###
# extract DAR bed
wd<-"/ngsprojects/scgrn/data_archive/clmen/cluster_DARs_second/"
peaks_cluster <- vector(mode="list", length=length(unique(Idents(integrated_final))))
names(peaks_cluster) <- unique(Idents(integrated_final))
filename<-"_DARs.bed"

for (i in identification_top$cluster){
  subset<-da_peaks[da_peaks$cluster==i,"gene"]
  if (length(subset)==0){  
    next
  } 
  else {
    peaks_cluster[identification_top$type[i+1]]<-length(subset)
  }
  id<-identification_top$type[i+1]
  id<-gsub(" ", "_", id)
  df <-  t(data.frame(strsplit(subset,'-')))
  colnames(df)<-c("chr","start","stop")
  df[,"chr"] <- paste0("Chr",df[,"chr"]) 
  write.table(df,file = paste0(wd,paste0(id,filename)), append=FALSE, quote = FALSE, 
              row.names = FALSE, col.names=FALSE, sep="\t")
}