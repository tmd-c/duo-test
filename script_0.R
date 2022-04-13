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
################################################################################
# k-means ######################################################################
################################################################################
################################################################################

################################################################################
# Euclidean ####################################################################
################################################################################

N = 4

resKmeans <- kmeans(expData, centers = N)

View(resKmeans)

# Clusters #####################################################################

clu1_eucl_kmea <- expData[which(resKmeans$cluster == 1),]
clu2_eucl_kmea <- expData[which(resKmeans$cluster == 2),]
clu3_eucl_kmea <- expData[which(resKmeans$cluster == 3),]
clu4_eucl_kmea <- expData[which(resKmeans$cluster == 4),]

# plotGenes ####################################################################

plotGenes(expData, yMax = 100000)

plotGenes(clu1_eucl_kmea, yMax = 100000)
plotGenes(clu2_eucl_kmea, yMax = 100000)
plotGenes(clu3_eucl_kmea, yMax = 100000)
plotGenes(clu4_eucl_kmea, yMax = 100000)

# Heatmaps #####################################################################

heatmap(as.matrix(expData))

heatmap(as.matrix(clu1_eucl_kmea))
heatmap(as.matrix(clu2_eucl_kmea))
heatmap(as.matrix(clu3_eucl_kmea))
heatmap(as.matrix(clu4_eucl_kmea))

################################################################################
# Correlation ##################################################################
################################################################################

cor(expData)

cor(t(expData))

matDist <- as.dist(1 - cor(t(expData)))

resKmeans2 <- kmeans(matDist, centers = 4)

# Clusters #####################################################################

clu1_corr_kmea <- expData[which(resKmeans2$cluster == 1),]
clu2_corr_kmea <- expData[which(resKmeans2$cluster == 2),]
clu3_corr_kmea <- expData[which(resKmeans2$cluster == 3),]
clu4_corr_kmea <- expData[which(resKmeans2$cluster == 4),]

# plotGenes ####################################################################

plotGenes(clu1_corr_kmea, yMax = 100000)
plotGenes(clu2_corr_kmea, yMax = 100000)
plotGenes(clu3_corr_kmea, yMax = 100000)
plotGenes(clu4_corr_kmea, yMax = 100000)

# heatmaps #####################################################################

heatmap(as.matrix(clu1_corr_kmea))
heatmap(as.matrix(clu2_corr_kmea))
heatmap(as.matrix(clu3_corr_kmea))
heatmap(as.matrix(clu4_corr_kmea))

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

plot(resHCL)

# Clusters #####################################################################

clu1_eucl_hcl <- expData[which(cutree(resHCL, k = N) == 1),]
clu2_eucl_hcl <- expData[which(cutree(resHCL, k = N) == 2),]
clu3_eucl_hcl <- expData[which(cutree(resHCL, k = N) == 3),]
clu4_eucl_hcl <- expData[which(cutree(resHCL, k = N) == 4),]

# plotGenes ####################################################################

plotGenes(clu1_eucl_hcl, yMax = 100000)
plotGenes(clu2_eucl_hcl, yMax = 100000)
plotGenes(clu3_eucl_hcl, yMax = 100000)
plotGenes(clu4_eucl_hcl, yMax = 100000)

# Heatmaps #####################################################################

heatmap(as.matrix(clu1_eucl_hcl))
heatmap(as.matrix(clu2_eucl_hcl))
heatmap(as.matrix(clu3_eucl_hcl))
heatmap(as.matrix(clu4_eucl_hcl))

################################################################################
# Correlation matrix ###########################################################
################################################################################

resHCL2 <- hclust(matDist)

plot(resHCL2)

# Clusters #####################################################################

clu1_corr_hcl <- expData[which(cutree(resHCL2, k = N) == 1),]
clu2_corr_hcl <- expData[which(cutree(resHCL2, k = N) == 2),]
clu3_corr_hcl <- expData[which(cutree(resHCL2, k = N) == 3),]
clu4_corr_hcl <- expData[which(cutree(resHCL2, k = N) == 4),]

# plotGenes ####################################################################

plotGenes(clu1_corr_hcl, yMax = 100000)
plotGenes(clu2_corr_hcl, yMax = 100000)
plotGenes(clu3_corr_hcl, yMax = 100000)
plotGenes(clu4_corr_hcl, yMax = 100000)

# Heatmaps #####################################################################

heatmap(as.matrix(clu1_corr_hcl))
heatmap(as.matrix(clu2_corr_hcl))
heatmap(as.matrix(clu3_corr_hcl))
heatmap(as.matrix(clu4_corr_hcl))