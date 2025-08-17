############## Install ##################
install.packages(c("emmeans", "dplyr", "ggplot2"))
library(emmeans)
library(dplyr)
library(ggplot2)  



######## Reload clean data ##################
# Clear workspace again
# Reload proteomics and metadata
Proteomics.Data = readRDS("data.rds")    
Meta.Data = readRDS("meta_with_pcs.rds")        
# Replace missing gene names with "Unknown"
if (any(is.na(Proteomics.Data$genes_lung))) {
  Proteomics.Data$genes_lung[is.na(Proteomics.Data$genes_lung)] <- "Unknown"
}
# Replace semicolons in gene names with underscores 
Proteomics.Data$genes_lung = gsub(";", "_", Proteomics.Data$genes_lung)
# Convert to data.frame, set gene names as row names, remove name column
Proteomics.Data = as.data.frame(Proteomics.Data)
rownames(Proteomics.Data) = make.unique(Proteomics.Data$genes_lung)
Proteomics.Data$genes_lung = NULL

# Align metadata to expression data
Meta.Data = Meta.Data[match(colnames(Proteomics.Data), Meta.Data$SampleID), ]
# This line verifies that after reordering, every sample in Meta.Data is in the correct order.
# If there's a mismatch, the script stops with an error, preventing incorrect downstream analysis.
if (!all(Meta.Data$SampleID == colnames(Proteomics.Data))) stop("Mismatch in sample order.")

# Extract experimental covariates
Group = as.factor(Meta.Data$Genotype)           
PC1 = as.numeric(Meta.Data$PC1)
PC2 = as.numeric(Meta.Data$PC2)
PC3 = as.numeric(Meta.Data$PC3)
Day = as.factor(Meta.Data$Day)           

# Initialize empty results dataframe
gene_names = rownames(Proteomics.Data)

library(emmeans)

results.lm = data.frame()
skipped_genes = 0

for (gene in gene_names) {
  gene_expression = as.numeric(Proteomics.Data[gene, ])
  
  dat <- data.frame(
    expr = gene_expression,
    Group = Group,
    PC1 = PC1,
    PC2 = PC2,
    Day = Day
  )
  
  dat$Group <- factor(dat$Group, levels = c("ctrl", "wt", "mut"))
  
  tryCatch({
    model <- lm(expr ~ Group + PC1, data = dat)
    
    # All pairwise comparisons between groups
    em <- emmeans(model, pairwise ~ Group)
    contrast_df <- as.data.frame(em$contrasts)
    
    results.lm <- rbind(results.lm, data.frame(
      Gene = gene,
      Contrast = contrast_df$contrast,
      Estimate = contrast_df$estimate,
      Std_Error = contrast_df$SE,
      t_value = contrast_df$t.ratio,
      p_value = contrast_df$p.value
    ))
    
  }, error = function(e) {
    message(paste("Skipping gene:", gene, "Error:", e$message))
    skipped_genes <<- skipped_genes + 1
  })
}


# check t dist
message(paste("Total genes skipped due to model errors:", skipped_genes))
head(results.lm)
results.lm <- results.lm[abs(results.lm$t_value) <10,]

hist(results.lm$t_value)
hist(results.lm$t_value,
     main = "Distribution of t-values from linear regression model",
     xlab = "t-value (Genotype effect size / SE)",
     ylab = "Number of proteins",
     breaks = 50,
     col = "gray20",
     border = "white")

sum(abs(results.lm$t_value<10))
length(results.lm$t_value)



# FDR correction per contrast
results.lm = results.lm %>%
  group_by(Contrast) %>%
  mutate(FDR = p.adjust(p_value, method = "fdr")) %>%
  arrange(FDR) %>%
  ungroup()

head(results.lm)


