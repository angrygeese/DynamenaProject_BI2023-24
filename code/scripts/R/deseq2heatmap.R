# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# install.packages("dplyr")
# install.packages("gplots")

library("dplyr")
library("gplots")

# draw heatmap
draw_heatmap <- function(file){
  # Read normalized counts
  data = read.table(file, header=T, sep="\t", as.is=TRUE)
  
  gene = data[,1]
  vals = as.matrix(data[,2:ncol(data)])
  
  # Adds a little noise to each element
  # To avoid the clusteing function failing on zero
  # variance datalines.
  vals = jitter(vals, factor = 1, amount=0.00001)
  
  
  # Calculate zscore
  score = NULL
  for (i in 1:nrow(vals)) {
    row=vals[i,]
    zscore=(row-mean(row))/sd(row)
    score =rbind(score,zscore)
  }
  
  row.names(score) = gene
  zscore=score
  
  # Generate heatmap
  mat = as.matrix(zscore)
  
  # Opent the drawing device.
  pdf('DMSO-AZK_heatmap.pdf')
  # png('../data_folder/processed_data/4. Expression count/DESeq2/DMSO-AZK_heatmap.png', width = 1920, height = 1920, units='px', res = 300)
  
  colors = colorRampPalette(c("green","black","red"),space="rgb")(256)
  heatmap.2(mat,col=colors,density.info="none",trace="none", margins=c(14,14),labCol=c("DMSO_1", "DMSO_2", "DMSO_3", "AZK_1", "AZK_2", "AZK_3"),lhei=c(1,5))
  
  invisible(dev.off())
}