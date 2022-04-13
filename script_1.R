# ██████╗ ███████╗███████╗██╗    ██████╗ 
# ██╔══██╗██╔════╝██╔════╝██║    ██╔══██╗
# ██║  ██║█████╗  █████╗  ██║    ██████╔╝
# ██║  ██║██╔══╝  ██╔══╝  ██║    ██╔══██╗
# ██████╔╝███████╗██║     ██║    ██║  ██║
# ╚═════╝ ╚══════╝╚═╝     ╚═╝    ╚═╝  ╚═╝
# Trenton Dailey-Chwalibóg, M.P.H., Ph.D.
# Universiteit Gent

setwd("C:/Users/tdaileyc/Google Drive/9 ugent/training/duo/defi_R")

library(ggplot2)
library(grid)
library(gridExtra)
library(cowplot)

################################################################################
# Import data ##################################################################
################################################################################

expData <- read.table("cell-cycle_SCERE_DUO.txt", 
                      row.names = 1, 
                      sep = "\t", 
                      header = T)

################################################################################
# Use Gaëlle's function ########################################################
################################################################################

plotGenes <- function(expData, title = "", yMin = 0, yMax = NULL, meanProfile = TRUE){
  
  # Check function parameters
  if(is.null(yMax)){
    
    print("You must specify a maximal value for Y axis")
    
  }else{
    
    # Representation of the first expression profile
    plot(1:ncol(expData), expData[1,], col = "grey", type = "l",
         ylim = c(floor(yMin), ceiling(yMax)),
         xlab = "Time point", ylab = "Gene expression level",
         main = title)
    
    # Add expression profile for other genes
    for(i in 2:nrow(expData)){
      
      lines(1:ncol(expData), expData[i,], col = "grey")
      
      # end of for()  
    }
    
    # Average expression profile
    if(meanProfile == TRUE){
      expMean = apply(expData, 2, mean)
      lines(1:ncol(expData), expMean, col = "red", 
            lwd = 1.5, lty = "dashed")
    }
    
    # end of else()   
  }
  
  # end of function plotGenes()  
}

################################################################################
# set cluster number ###########################################################
################################################################################

N = 7

################################################################################
################################################################################
# k-means ######################################################################
################################################################################
################################################################################

################################################################################
# Euclidean ####################################################################
################################################################################

resKmeans <- kmeans(expData, centers = N)

for(i in 1:N){
  clus_eucl_kmea <- expData[which(resKmeans$cluster == i),]
  plotGenes(clus_eucl_kmea, yMax = 100000)
  if(nrow(clus_eucl_kmea)>1){
    heatmap(as.matrix(clus_eucl_kmea))}
}

################################################################################
# Correlation ##################################################################
################################################################################

cor(expData)

cor(t(expData))

matDist <- as.dist(1 - cor(t(expData)))

resKmeans2 <- kmeans(matDist, centers = N)

for(i in 1:N){
  clu1_corr_kmea <- expData[which(resKmeans2$cluster == i),]
  plotGenes(clu1_corr_kmea, yMax = 100000)
  heatmap(as.matrix(clu1_corr_kmea))
}

################################################################################
################################################################################
# HCL ##########################################################################
################################################################################
################################################################################

################################################################################
# Euclidean #####################################################################
################################################################################

d <- dist(expData)

resHCL <- hclust(d)

for(i in 1:N){
  clus_eucl_hcl <- expData[which(cutree(resHCL, k = N) == i),]
  plotGenes(clus_eucl_hcl, yMax = 100000)
  heatmap(as.matrix(clus_eucl_hcl))
}

################################################################################
# Correlation matrix ###########################################################
################################################################################

resHCL2 <- hclust(matDist)

for(i in 1:N){
  clus_corr_hcl <- expData[which(cutree(resHCL2, k = N) == i),]
  plotGenes(clus_corr_hcl, yMax = 100000)
  heatmap(as.matrix(clus_corr_hcl))
}
