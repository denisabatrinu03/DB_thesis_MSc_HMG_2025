########## Install packages ########################################

install.packages(c("ggplot2", "factoextra"))
library(ggplot2)
library(factoextra) 


########### Input data #########################################
# Clear environment
rm(list=ls())
# Load data 
Proteomics.Data = readRDS("/Users/denisabatrinu/Desktop/Research_Project_25_PCD_PG/MS_analysis_22APR/data.rds")  
Proteomics.Data = as.data.frame(Proteomics.Data) # make sure data is in data.frame format
# Replace missing gene names in the "genes_lung" column with "Unknown"
Proteomics.Data$genes_lung[is.na(Proteomics.Data$genes_lung)] <- "Unknown"
# Set unique row names based on the first column ("genes_lung" column)
rownames(Proteomics.Data) = make.unique(Proteomics.Data[, 1])
# Remove the first column after assigning it as row names (avoid duplication)
Proteomics.Data = Proteomics.Data[, -1]    
# Load metadata
Meta.Data = readRDS("/Users/denisabatrinu/Desktop/Research_Project_25_PCD_PG/MS_analysis_22APR/meta_with_pcs.rds") 



######## Run PCA on transposed proteomics data (samples as rows, proteins as columns) ############
pca_result = prcomp(as.matrix(t(Proteomics.Data)), center = TRUE, scale. = TRUE)
# Scree plot
fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 50)) +
  labs(title = "PCA Scree Plot") +
  theme_minimal()
# Extract PC scores and match to Sample ID to create a merged data frame
pc_scores = as.data.frame(pca_result$x)
pc_scores$SampleID = rownames(pc_scores)
merged_data = merge(pc_scores, Meta.Data, by = "SampleID")
# Plot PC1 vs PC2 
ggplot(merged_data, aes(x = PC1.x, y = PC2.x, color = Genotype, shape = Day)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "PC1 vs PC2",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)")) +
  theme_minimal() +
  theme(text = element_text(size = 14))
