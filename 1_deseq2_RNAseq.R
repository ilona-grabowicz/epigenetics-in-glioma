# Author: Ilona E. Grabowicz
# Date: 2019-2020

# Identifying DEGs using DeSeq2 tool.

library("DESeq2")
# DE genes (RNAseq data)
data <- read.csv('./input_files/RNAseq/combined_RNAs_non_normalized_no_underscores.csv', header=T, stringsAsFactors = F, row.names = 1, sep='\t')

run_deseq <- function(subtype_1, subtype_2){
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
  write.csv(resOrdered, paste0('DeSeq2_', levels[1],'_', levels[2],'.csv'), sep='\t') 
}
run_deseq('PA', 'DA')
run_deseq('DA', 'GB|PG')
run_deseq('PA', 'GB|PG')  
