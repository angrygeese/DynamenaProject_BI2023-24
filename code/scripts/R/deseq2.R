# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")
# install.packages("dplyr")
# install.packages("gplots")

library("DESeq2")
library("dplyr")

deseq2_workflow_ctrl <- function(file){
  # Set up the conditions based on the experimental setup.
  cond_1 = rep("DMSO", 3)
  cond_2 = rep("AZK", 3)
  
  # Read the data from the standard input.
  countData = read.table(file, header=TRUE, sep="\t", row.names=1)
  
  # Build the dataframe from the conditions
  samples = names(countData)
  condition = factor(c(cond_1, cond_2))
  colData = data.frame(samples=samples, condition=condition)
  
  # Create DESEq2 dataset.
  dds = DESeqDataSetFromMatrix(countData=countData, colData=colData, design = ~condition)
  
  #Set the reference to be compared
  dds$condition = relevel(dds$condition,"DMSO")
  
  # # Control features for estimating size factors (normalize counts)
  ctrlGenes <- c(which(rownames(countData) == "NODE_35809_length_1886_cov_615.371677_g1855_i7"))
  dds <- estimateSizeFactors(dds, controlGenes=ctrlGenes)
  
  # Run deseq
  dds = DESeq(dds)
  
  # Format the results.
  res = results(dds)
  
  # Sort the results data frame by the padj and foldChange columns.
  sorted = res[with(res, order(padj, -log2FoldChange)), ]
  
  # Turn it into a dataframe to have proper column names.
  sorted.df = data.frame("id"=rownames(sorted),sorted)
  
  # Write the table out.
  write.table(sorted.df, file="../data_folder/processed_data/4.Expression_count/DESeq2/result_ctrl.txt", sep="\t", col.names=NA, quote=FALSE)
  
  # Get normalized counts and write this to a file
  nc = counts(dds,normalized=TRUE)
  
  # Turn it into a dataframe to have proper column names.
  dt = data.frame("id"=rownames(nc),nc)
  
  # Save the normalize data matrix.
  write.table(dt, file="../data_folder/processed_data/4.Expression_count/DESeq2/norm-matrix-deseq2_ctrl.txt", sep="\t",  row.name=FALSE, col.names=TRUE,quote=FALSE)
}