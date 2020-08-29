# Author: Ilona E. Grabowicz
# Date: 2019-2020

# Identifying genes with DEMs using DeSeq2 tool.

library("DESeq2")
library("dplyr")

data_type <- 'H3K27ac' 
gene_part <- 'promoters'
run_deseq <- function(subtype_1, subtype_2, data_type, gene_part){
  if (data_type=='DNAmeth'){
  data <- read.csv(paste0('./input_files/ChIPseq_ATAC_DNase/methyl_table_', gene_part, '_cov_od80.tsv'), header=T, stringsAsFactors = F, row.names = 1, sep='\t')
  }
  else{
    data <- read.csv(paste0('./input_files/ChIPseq_ATAC_DNase/', data_type, '_coverage_on_', gene_part, '.tsv'), header=T, stringsAsFactors = F, sep='\t')
    data$X <- NULL
    data <- distinct(data, gene_ID, .keep_all = TRUE)
    rownames(data) <- data[,1]
    data$gene_ID <- NULL
    }
  my_pattern_list<-list(subtype_1, subtype_2)
  cts <- data
  comparison <- grepl(paste(my_pattern_list, collapse='|'), colnames(cts))
  cts <- cts[,comparison]
  coldata <- data.frame(matrix(nrow=ncol(cts),ncol=1)) # Metadata table
  rownames(coldata) <- colnames(cts)
  colnames(coldata) <- c('grade')
  coldata$grade[grepl('PA', rownames(coldata))] <- 'grade1'
  coldata$grade[grepl('DA', rownames(coldata))] <- 'grade23'
  coldata$grade[grepl('GB|PG', rownames(coldata))] <- 'grade4'
  coldata$grade<-factor(coldata$grade)
  dds <- DESeqDataSetFromMatrix(countData = cts, # DeSeq object
                                colData = coldata,
                                design = ~ grade)
  levels<-c()
  if (subtype_1=='PA'|subtype_2=='PA'){
    levels <- c(levels, 'grade1')
  }
  if (subtype_1=='DA'|subtype_2=='DA'){
    levels <- c(levels, 'grade23')
  }
  if (subtype_1=='GB|PG'|subtype_2=='GB|PG'){
    levels <- c(levels, 'grade4')
  }
  
  dds$grade <- factor(dds$grade, levels = levels) 
  dds <- DESeq(dds)
  res <- results(dds)
  resOrdered <- data.frame(res[order(res$pvalue),]) # Result
  write.csv(resOrdered, paste0('DeSeq2_', data_type, '_', levels[1],'_', levels[2], '_', gene_part, '.csv')) 
}

run_deseq('PA', 'DA', 'H3K27ac', 'promoters')
run_deseq('PA', 'DA', 'H3K27ac', 'genes')
run_deseq('PA', 'DA', 'H3K4me3', 'promoters')
run_deseq('PA', 'DA', 'H3K4me3', 'genes')
run_deseq('PA', 'DA', 'DNase', 'promoters')
run_deseq('PA', 'DA', 'DNase', 'genes')
run_deseq('PA', 'DA', 'ATACseq', 'promoters')
run_deseq('PA', 'DA', 'ATACseq', 'genes')
run_deseq('PA', 'DA', 'DNAmeth', 'promoters')
run_deseq('PA', 'DA', 'DNAmeth', 'genes')

run_deseq('PA', 'GB|PG', 'H3K27ac', 'promoters')
run_deseq('PA', 'GB|PG', 'H3K27ac', 'genes')
run_deseq('PA', 'GB|PG', 'H3K4me3', 'promoters')
run_deseq('PA', 'GB|PG', 'H3K4me3', 'genes')
run_deseq('PA', 'GB|PG', 'DNase', 'promoters')
run_deseq('PA', 'GB|PG', 'DNase', 'genes')
run_deseq('PA', 'GB|PG', 'ATACseq', 'promoters')
run_deseq('PA', 'GB|PG', 'ATACseq', 'genes')
run_deseq('PA', 'GB|PG', 'DNAmeth', 'promoters')
run_deseq('PA', 'GB|PG', 'DNAmeth', 'genes')

run_deseq('DA', 'GB|PG', 'H3K27ac', 'promoters')
run_deseq('DA', 'GB|PG', 'H3K27ac', 'genes')
run_deseq('DA', 'GB|PG', 'H3K4me3', 'promoters')
run_deseq('DA', 'GB|PG', 'H3K4me3', 'genes')
run_deseq('DA', 'GB|PG', 'DNase', 'promoters')
run_deseq('DA', 'GB|PG', 'DNase', 'genes')
run_deseq('DA', 'GB|PG', 'ATACseq', 'promoters')
run_deseq('DA', 'GB|PG', 'ATACseq', 'genes')
run_deseq('DA', 'GB|PG', 'DNAmeth', 'promoters')
run_deseq('DA', 'GB|PG', 'DNAmeth', 'genes')

