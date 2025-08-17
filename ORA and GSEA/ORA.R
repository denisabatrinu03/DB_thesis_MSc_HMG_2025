###### Install packages #######
# Install
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "enrichplot", "org.Mm.eg.db", "DOSE"))
install.packages(c("ggplot2", "stringr"))

# Load
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(ggplot2)
# library(stringr)   # optional since call stringr::str_wrap()


# results.lm <- read_excel()


# Filter for significantly differentially expressed genes by Group (mut or wt)
results.lm.Groupmut.sig = results.lm[results.lm$Contrast == "ctrl - mut", ]
results.lm.Groupwt.sig  = results.lm[results.lm$Contrast == "ctrl - wt", ]
results.lm.Groupwtmut.sig  = results.lm[results.lm$Contrast == "wt - mut", ]


# ---- Please read the following comments carefully before running the code below.

############## Define background gene universe ###############################
# Use all protein-coding genes from this dataset background
background_genes = unique(rownames(Proteomics.Data))
# Convert background gene symbols to ENTREZ IDs
gene2entrez_bg = bitr(background_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
background_entrez = unique(gene2entrez_bg$ENTREZID) # Returns error that 20% of genes are missed...
# Check % mapped
mapped_percent = nrow(gene2entrez_bg) / length(background_genes) * 100
cat(sprintf("Mapped %.1f%% of gene symbols to Entrez IDs.\n", mapped_percent))
# Optional: see failed ones
failed = setdiff(background_genes, gene2entrez_bg$SYMBOL)
head(failed)
print(failed) # Vector is too long, R will truncate putput. Use: 
cat(failed, sep = "\n") # instead. OR, save to a file to analyse later:
writeLines(failed, "unmapped_gene_symbols.txt")


# ############## Select significant gene lists for ORA ######################

# Filter genes with significant differential expression (FDR < 0.05)
results.lm.Groupmut.sig = results.lm[results.lm$Contrast == "ctrl - mut" & results.lm$FDR < 0.05, ]
results.lm.Groupwt.sig  = results.lm[results.lm$Contrast == "ctrl - wt"  & results.lm$FDR < 0.05, ]
results.lm.Groupwtmut.sig  = results.lm[results.lm$Contrast == "wt - mut"  & results.lm$FDR < 0.05, ]


# Convert significant gene symbols to ENTREZ IDs
entrez_mut = bitr(results.lm.Groupmut.sig$Gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
entrez_wt  = bitr(results.lm.Groupwt.sig$Gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID


############## Perform GO Enrichment Analysis (Biological Process) ##################
go_mut = enrichGO(gene = entrez_mut,
                  OrgDb = org.Mm.eg.db,
                  #  universe = background_entrez,  # Optional: set background genes
                  keyType = "ENTREZID",
                  ont = "CC",           # Biological Process category
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2)

go_wt = enrichGO(gene = entrez_wt,
                 OrgDb = org.Mm.eg.db,
                 # universe = background_entrez,
                 keyType = "ENTREZID",
                 ont = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2)


################ Review Top Enriched Terms ##################################
head(as.data.frame(go_mut), 10)
MT=as.data.frame(go_mut)

head(as.data.frame(go_wt), 10)
WT=as.data.frame(go_wt)


# Check how many terms are enriched
length(WT$Description)
head(WT$Description)

length(MT$Description)
head(MT$Description)


########### List the proteins associated to a specific pathway of interest ############
# Example: inspect genes enriched in “protein modification by small protein removal”
First.Path=go_wt[go_wt$Description %in% c("protein modification by small protein removal"),]

# Extract gene IDs from the selected term and convert to gene symbols
A = gsub("/", ",", unique(First.Path$geneID)) # Replace '/' with ',' to allow splitting
A_vector = unlist(strsplit(A, ",")) # Split into individual ENTREZ IDs

# Convert ENTREZ IDs back to gene symbols for interpretation
bitr_result = bitr(A_vector, fromType = "ENTREZID", 
                   toType = "SYMBOL", 
                   OrgDb = org.Mm.eg.db)

bitr_result


################### Visualising GSEA results ##########################
library(enrichplot) # For enrichment result visualization

# Bar plots
# For top GO terms in go_mut
barplot(go_mut, showCategory = 20, title = "Top 20 Enriched GO Terms (MT vs CTRL)")
p1 <- barplot(go_mut, showCategory = 20, title = "Top 20 Enriched GO Terms (MT vs CTRL)")
p1 + theme(axis.text.y = element_text(size = 10)) +
  scale_y_discrete(labels = function(labels) stringr::str_wrap(labels, width = 50))
# For WT
barplot(go_wt, showCategory = 20, title = "Enriched GO Terms (WT vs CTRL)")
p2 <- barplot(go_wt, showCategory = 20, title = "Enriched GO Terms (WT vs CTRL)")
p2 + theme(axis.text.y = element_text(size = 10)) +
  scale_y_discrete(labels = function(labels) stringr::str_wrap(labels, width = 50))
