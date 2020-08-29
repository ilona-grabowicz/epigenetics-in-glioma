## Purpose of script: Copy number variants calls from RNAseq (CaSpER tool)
## Author: Adria Roura Canalda
## Date Created: 2020-05-29
# Necessary files need to be unpacked first: retaskstoilonasmanuscript.zip

library(CaSpER)
library(GenomicRanges)
library(BiocParallel)
register(MulticoreParam(8))

#read raw count matrix
counts <- read.table("symfonia_ALL_RNAseq_counts.featureCounts", head=TRUE, row.names = 1)
exon_length <- counts$Length
counts <- counts[,-c(1,2,3,4,5,18,33)] #irrelevant info and low quality samples
counts <- counts[,-c(18,22)] #discarding samples that did not pass casper baf extract
counts <- counts[,-c(31,32,33)] #remove some NB

#convert to TPMs
tpm3 <- function(counts,len) {
  x <- counts/len
  return(t(t(x)*1e6/colSums(x)))
}

counts_tpm <- tpm3(counts, exon_length)

#cytoband information
cytoband <- read.delim("cytoBand.txt", header=F)
cytoband <- data.frame(V1=gsub("chr", "", cytoband[,1]), V2=cytoband[,2], V3=cytoband[,3], V4=substring(cytoband$V4, 1, 1), stringsAsFactors=F)
start <- do.call(rbind, lapply(split(cytoband$V2, paste0(cytoband$V1, cytoband$V4)), min))
end <- do.call(rbind, lapply(split(cytoband$V3, paste0(cytoband$V1, cytoband$V4)), max))
cytoband <- data.frame(V1=gsub("p", "", gsub("q", "", rownames(start))), V2=start, V3=end, V4=rownames(start), stringsAsFactors=F)
cytoband <- cytoband [as.vector(unlist(sapply(c(1:22, "X"), function(x) which(cytoband$V1 %in% x)))), ]
cytoband$V4[grep("q", cytoband$V4)] <- "q"
cytoband$V4[grep("p", cytoband$V4)] <- "p"
rownames(cytoband) <- NULL

#annotations, important in hg38
annotation <- generateAnnotation(id_type="ensembl_gene_id", genes=rownames(counts_tpm), ishg19 = FALSE, centromere_hg38)

#error: different dimensions counts and annotations throws and error
keep <- rownames(counts_tpm) %in% annotation$Gene
counts_tpm <- counts_tpm[keep,]

#read baf files
loh <- readBAFExtractOutput(path="./baf_files/", sequencing.type="bulk", suffix="baf")
names(loh) <- gsub(".baf", "", names(loh))
colnames(counts_tpm) <- names(loh)

loh.name.mapping <- data.frame(loh.name = names(loh), sample.name = names(loh))

#1) Generate casper object
object <- CreateCasperObject(raw.data=counts_tpm,loh.name.mapping=loh.name.mapping, sequencing.type="bulk", 
                             cnv.scale=3, loh.scale=3, matrix.type="normalized", expr.cutoff=4.5,
                             annotation=annotation, method="iterative", loh=loh, filter="median",  
                             control.sample.ids=c("RNAseq_NB1","RNAseq_NB2","RNAseq_NB3","RNAseq_NB5") , cytoband=cytoband, genomeVersion = "hg38", log.transformed = FALSE)

#2) Pairwise comparison of scales from BAF and expression signals
final.objects <- runCaSpER(object, removeCentromere=TRUE, cytoband=cytoband, method="iterative")

#3) Harmonization and Summarization of CNV calls from multiple scales and from multiple pairwise comparison of BAF and Expression Signals (Large-Scale, segment-based, gene-based)

# a) Large scale CNV (each of the chr arms = amp, del, neutral)
finalChrMat <- extractLargeScaleEvents(final.objects, thr=0.75) 

# b) Segment-based CNV
gamma <- 6
all.segments <- do.call(rbind, lapply(final.objects, function(x) x@segments))
segment.summary <- extractSegmentSummary(final.objects)
loss <- segment.summary$all.summary.loss
gain <- segment.summary$all.summary.gain
loh <- segment.summary$all.summary.loh
loss.final <- loss[loss$count>gamma, ]
gain.final <- gain[gain$count>gamma, ]
loh.final <- loh[loh$count>gamma, ]

# c) Gene-based CNV (1 amp, -1 del, 0 neutral)
all.summary<- rbind(loss.final, gain.final)
colnames(all.summary) [2:4] <- c("Chromosome", "Start",   "End")
rna <-  GRanges(seqnames = Rle(gsub("q", "", gsub("p", "", all.summary$Chromosome))), 
                IRanges(all.summary$Start, all.summary$End))   
ann.gr <- makeGRangesFromDataFrame(final.objects[[1]]@annotation.filt, keep.extra.columns = TRUE, seqnames.field="Chr")
hits <- findOverlaps(rna, ann.gr)
genes <- splitByOverlap(ann.gr, rna, "GeneSymbol")
genes.ann <- lapply(genes, function(x) x[!(x=="")])
all.genes <- unique(final.objects[[1]]@annotation.filt[,2])
all.samples <- unique(as.character(final.objects[[1]]@segments$ID))
rna.matrix <- gene.matrix(seg=all.summary, all.genes=all.genes, all.samples=all.samples, genes.ann=genes.ann)

#save segment-based CNV
write.table(all.summary, file = "segment_based_CNV_CaSpER_all_samples.txt",row.names=F, col.names=T, quote=F, sep="\t")

#save segment-based CNV
write.table(rna.matrix, file = "gene_based_CNV_CaSpER_all_samples.txt",row.names=T, col.names=T, quote=F, sep="\t")


#Visualization
library(GenomeGraphs)

# a) plotHeatmap: gene expression from different smoothing scales (row=samples, column=chromosomes)
obj <- final.objects[[9]]
plotHeatmap(object=obj, fileName="heatmap.png",cnv.scale= 3, cluster_cols = T, cluster_rows = T, show_rownames = T, only_soi = T)

# b) plot large scale events among all samples (large amplifications and deletions)
plotLargeScaleEvent(object=obj, fileName="large.scale.events2.pdf") 

# c) plot large scale events using event summary matrix 1: amplification, -1:deletion, 0: neutral (75% at least 7 out of 9 scales, of consistent CNV call)
plotLargeScaleEvent2(finalChrMat, fileName="large.scale.events.summarized.pdf")

# d) plotBAFAllSamples: plot BAF deviation for all samples
plotBAFAllSamples(loh = obj@loh.median.filtered.data,  fileName="LOHAllSamples.png")

# e) plot BAF deviation for each sample in seperate pages 
plotBAFInSeperatePages(loh =obj@loh.median.filtered.data, folderName="LOHPlotsSeperate") 

# f) plot BAF signal in different scales for all samples
plotBAFOneSample(object, fileName="LohPlotsAllScales.pdf") 

# g) plotGEAndBAFOneSample: Gene expression and BAF signal for one sample in one plot
plotGEAndBAFOneSample(object=obj, cnv.scale=3, loh.scale=3, sample= "RNAseq_GB01")

#Intersect Ilona's TADs with CNV segmend-based
library(tidygenomics)
library(dplyr)
segs<- rbind(loss.final, gain.final)
colnames(segs) [2:4] <- c("Chromosome", "Start",   "End")
segs$chr <- gsub(segs$Chromosome,pattern = c("q|p"), replacement = "")
segs$chr <- paste("chr", segs$chr, sep = "")
segs$chr <- as.factor(segs$chr)
segs <- arrange(segs,chr, Start, End)
segs <- segs[,c(1,2,9,3,4,5,6,7,8)]
write.table(segs[c(3,4,5)], file = "segment_based_CNV_CaSpER_all_samples.bed",row.names=F, col.names=T, quote=F, sep="\t")

#Ilona's TADS
TADS <- read.table("TADS_coordinates_rm_chrX.csv", head=TRUE, row.names = 1, sep = ","); TADS <- TADS[,c(1:3)]
TADS$TAD_nr <- rownames(TADS)
colnames(TADS) <- c("chr", "Start", "End", "TAD_nr")
TADS <- arrange(TADS, chr, Start,End)
#write.table(TADS, file = "TADS_IGV.bed",row.names=F, col.names=T, quote=F, sep="\t")

segs_TADS <- genome_intersect(TADS, segs, by=c("chr", "Start", "End"), mode = "both")
#write.table(segs_TADS, file = "segment_based_CNV_CaSpER_all_samples_intersected_TADs.bed",row.names=F, col.names=T, quote=F, sep="\t")

#filtering based on patients (at least 2) type of structural variation 
segs_TADS_filtered <- segs_TADS
segs_TADS_filtered$merge <- paste(segs_TADS_filtered$chr, segs_TADS_filtered$Start, segs_TADS_filtered$type, sep = ":")

#select only duplicated cases
segs_TADS_filtered <- segs_TADS_filtered %>% group_by(merge) %>% filter(n() > 1)
segs_TADS_filtered <- as.data.frame(segs_TADS_filtered)
segs_TADS_filtered <- segs_TADS_filtered[duplicated(segs_TADS_filtered$merge),] 

#write.table(segs_TADS_filtered, file = "segment_based_CNV_CaSpER_intersected_TADs_rm_1sample.bed",row.names=F, col.names=T, quote=F, sep="\t")
#write.table(segs_TADS_filtered[,c(1,9,10)], file = "segment_based_CNV_CaSpER_intersected_TADs_rm_1sample_IGV.bed",row.names=F, col.names=T, quote=F, sep="\t")

