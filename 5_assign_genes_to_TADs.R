# Author: Ilona E. Grabowicz
# Date: 2019-2020

library(tidyr)
library(dplyr)
library(data.table)

############ This script allows to assign genes to TADs by their chromosomal coordinates and save the assignments.

# dataset <- 'RNASeq'
# comparison <-  'grade1_grade4'
# part <- 'promoters'
#############
funkcja <- function(dataset, comparison, part) {
  TADs <- read.csv('./input_files/TADborders-Fetal_brain_D_cortical_plates_Won_et_al_2016.bed', stringsAsFactors = F)
  genes <- read.csv(paste0('DeSeq2_', comparison, '.csv'), sep='\t')
  genes <- separate(genes,col = gene, into=c('ENSEMBL', 'symbol', 'location'), sep = '_') 
  genes <- separate(genes, col=location, into=c('chromosome', 'rest'), sep=':')
  genes <- separate(genes, col=rest, into=c('start', 'stop'))
  genes <- separate(genes, col=chromosome, into=c('chromosome', 'x'))
  genes$x<-NULL
  genes$start<-as.integer(genes$start)
  genes$stop<-as.integer(genes$stop)
  # Subsetting to the same genes in all datasets:
  DNAme_proms <- read.csv(paste0('DeSeq2_DNAme_', comparison, '_', part, '.csv'), sep='\t', stringsAsFactors = F)
  DNAme_genes <- read.csv(paste0('DeSeq2_DNAme_', comparison, '_', part, '.csv'), sep='\t', stringsAsFactors = F)
  if (dataset == 'RNASeq'){
    DE <- genes }
  if (dataset!='RNASeq') {
    DE <- read.csv(paste0('../2018_04_05_DESeq2_glioma_grades_', dataset, '_', part, '/DeSeq2_', comparison, '.csv'), sep='\t', stringsAsFactors = F)
    genes_coords <- genes[,c('ENSEMBL', 'symbol', 'chromosome', 'start', 'stop')]
    names(DE)[names(DE) == 'gene'] <- 'ENSEMBL'
    DE <- inner_join(x=DE, y=genes_coords, by = c("ENSEMBL" = "ENSEMBL"))
  }
  # Which genes are in which TADs?
  list_of_TADs <-list()
  TADs_with_symbols <- list()
  for (row in 1:nrow(TADs)) { #nrow(TADs)
   # cat('.')
    my_DE <- DE[DE$chromosome == TADs[row,"chr"],]
    is.included <- mapply(function(de_start,de_stop,tad_start,tad_stop){de_start >= tad_start &  de_stop <= tad_stop}, my_DE$start, my_DE$stop, TADs[row,]$start, TADs[row,]$stop)
    list_of_TADs[[row]] <- my_DE[is.included,'log2_FC']
    TADs_with_symbols[[row]] <- my_DE[is.included,'ENSEMBL'] #symbol - then instead of ENSEMBL symbols are used
  }
  # Remove empty TADs and Kruskal-Wallis test
  list_of_TADs_noNA <- lapply(list_of_TADs, function(x) x[!is.na(x)])
  list_of_TADs_no_empty <- list_of_TADs_noNA[lapply(list_of_TADs_noNA,length)>0] # Limiting to TADs containing desired minimum number of genes
  a <- kruskal.test(list_of_TADs_no_empty)
  print(a$p.value)
  print(a)
  
  TADs_with_symbols <- lapply(TADs_with_symbols, function(x) {as.data.frame(t(x))})
  TADs_with_symbols <- lapply(TADs_with_symbols, function(x) { if(length(x)==0) data.frame(V1 = "xxx") else x} )
  TADs_with_symbols <- rbindlist(TADs_with_symbols, fill = T)
  TADs_with_symbols$V1[TADs_with_symbols$V1 == "xxx"] <- NA
  names(TADs_with_symbols) <- paste0("gene_",1:ncol(TADs_with_symbols))
  dim(TADs_with_symbols)
  means <- lapply(list_of_TADs, function(x) mean(x, na.rm = T)) #result is a list
  medians <- lapply(list_of_TADs, function(x) median(x, na.rm = T))
  TADs_df <- data.frame(TAD_nr=1:nrow(TADs), chr=TADs$chr, start=TADs$start, stop=TADs$stop,means=unlist(means), medians=unlist(medians), TADs_with_symbols)
  sorted_TADs <- arrange(TADs_df, desc(means))
  if (dataset=='RNASeq'){
  write.csv(sorted_TADs, paste0('./TADs/TADs_with_genes_ENSEMBL_and_FC_', dataset, '_', comparison, '.csv'))
  }
  else{
    write.csv(sorted_TADs, paste0('./TADs/TADs_with_genes_ENSEMBL_and_FC_', dataset, '_', comparison,'_', part, '.csv'))
  }
}

rnaseq <- funkcja('RNASeq', 'grade1_grade23', 'genes')
rnaseq <- funkcja('RNASeq', 'grade23_grade4', 'genes')
rnaseq <- funkcja('RNASeq', 'grade1_grade4', 'promoters')
h3k4me3 <- funkcja('H3K4me3', 'grade1_grade4', 'promoters')
h3k4me3 <- funkcja('H3K4me3', 'grade1_grade23', 'promoters')
h3k4me3 <- funkcja('H3K4me3', 'grade23_grade4', 'promoters')
dnase <- funkcja('DNase', 'grade1_grade4', 'genes')
dnase <- funkcja('DNase', 'grade1_grade23', 'genes')
dnase <- funkcja('DNase', 'grade23_grade4', 'genes')
dnase <- funkcja('DNase', 'grade1_grade4', 'promoters')
dnase <- funkcja('DNase', 'grade1_grade23', 'promoters')
dnase <- funkcja('DNase', 'grade23_grade4', 'promoters')
atac <- funkcja('ATACseq', 'grade1_grade4', 'genes')
atac <- funkcja('ATACseq', 'grade1_grade23', 'genes')
atac <- funkcja('ATACseq', 'grade23_grade4', 'genes')
atac <- funkcja('ATACseq', 'grade1_grade4', 'promoters')
atac <- funkcja('ATACseq', 'grade1_grade23', 'promoters')
atac <- funkcja('ATACseq', 'grade23_grade4', 'promoters')
h3k27ac<- funkcja('H3K27ac', 'grade1_grade4', 'genes')
h3k27ac <- funkcja('H3K27ac', 'grade1_grade23', 'genes')
h3k27ac <- funkcja('H3K27ac', 'grade23_grade4', 'genes')
h3k27ac<- funkcja('H3K27ac', 'grade1_grade4', 'promoters')
h3k27ac <- funkcja('H3K27ac', 'grade1_grade23', 'promoters')
h3k27ac <- funkcja('H3K27ac', 'grade23_grade4', 'promoters')
dname <- funkcja('DNAme', 'grade1_grade4', 'genes')
dname <- funkcja('DNAme', 'grade1_grade23', 'genes')
dname <- funkcja('DNAme', 'grade23_grade4', 'genes')
dname <- funkcja('DNAme', 'grade1_grade4', 'promoters')
dname <- funkcja('DNAme', 'grade1_grade23', 'promoters')
dname <- funkcja('DNAme', 'grade23_grade4', 'promoters')
