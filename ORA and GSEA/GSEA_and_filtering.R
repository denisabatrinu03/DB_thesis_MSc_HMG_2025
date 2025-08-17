############## Install packages ###############
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "DOSE"))  # DOSE is a dependency of clusterProfiler
install.packages(c("dplyr", "tidyr"))

library(dplyr)
library(tidyr)
library(clusterProfiler)

# results.lm <- read_excel()


# Define genes of interest (cilia-related or PAQosome interactors)
aA = c("Ruvbl1", "Ruvbl2", "Rpap3", "Pih1d1", "Hsp90", 
       "Uri1", "Uxt", "Pdrg1", "Pfdn2", "Pfdn6", "Asdurf", 
       "Wdr92", "Polr2e", "Rpb5", "Telo2", "Tti1", "Tti2", 
       "Nufip1", "Znhit3", "Znhit6", "Znhit2", "Ecd", "Aar2", 
       "Naf1", "Shq1", "Nopchap1", "Hop", "p23") 

# Filter for significantly differentially expressed genes by Group (mut or wt)
results.lm.Groupmut.sig = results.lm[results.lm$Contrast == "ctrl - mut", ]
results.lm.Groupmut.sig[results.lm.Groupmut.sig$Gene %in% aA, ]

results.lm.Groupwt.sig  = results.lm[results.lm$Contrast == "ctrl - wt", ]
results.lm.Groupwt.sig[results.lm.Groupwt.sig$Gene %in% aA, ]

results.lm.Groupwtmut.sig  = results.lm[results.lm$Contrast == "wt - mut", ]
results.lm.Groupwtmut.sig[results.lm.Groupwtmut.sig$Gene %in% aA, ]

# Save the filtered significant genes
install.packages("writexl")   
library(writexl)

ctrl_mut_aA   <- dplyr::filter(results.lm.Groupmut.sig,  Gene %in% aA)
ctrl_wt_aA    <- dplyr::filter(results.lm.Groupwt.sig,   Gene %in% aA)
wtmut_aA <- dplyr::filter(results.lm.Groupwtmut.sig, Gene %in% aA)

write_xlsx(
  list(
    ctrl_vs_mut_aA = ctrl_mut_aA,
    ctrl_vs_wt_aA  = ctrl_wt_aA,
    wt_vs_mut_aA   = wtmut_aA
  ),
  path = "filter_PAQ_lm_results_by_contrast_after_tdisttresh10.xlsx"
)


######## Pathway analysis #############

results.lm.Groupwt.sig.New=results.lm.Groupwt.sig[abs(results.lm.Groupwt.sig$t_value) < 4,]
results.lm.Groupmut.sig.New=results.lm.Groupmut.sig[abs(results.lm.Groupmut.sig$t_value) < 4,]
results.lm.Groupwtmut.sig.New=results.lm.Groupwtmut.sig[abs(results.lm.Groupwtmut.sig$t_value) < 4,]


## Also repeat the following code for the other 2 contrasts:
# This run is ctrl vs mut:
Genes_Ranked_mut = results.lm.Groupmut.sig.New %>%
  drop_na(Estimate, p_value) %>%
  arrange(desc(Estimate)) %>%
  mutate(rank_metric = -log10(p_value+1) * sign(Estimate)) %>%
  arrange(desc(rank_metric)) %>%
  filter(!is.na(Gene)) %>%
  distinct(Gene, .keep_all = TRUE) 

# This run is ctrl vs wt:
Genes_Ranked_wt = results.lm.Groupwt.sig.New %>%
  drop_na(Estimate, p_value) %>%
  arrange(desc(Estimate)) %>%
  mutate(rank_metric = -log10(p_value+1) * sign(Estimate)) %>%
  arrange(desc(rank_metric)) %>%
  filter(!is.na(Gene)) %>%
  distinct(Gene, .keep_all = TRUE) 

# This run is wt vs mut:
Genes_Ranked_wtmut = results.lm.Groupwtmut.sig.New %>%
  drop_na(Estimate, p_value) %>%
  arrange(desc(Estimate)) %>%
  mutate(rank_metric = -(-log10(p_value + 1) * sign(Estimate))) %>% # WT > MUT
  arrange(desc(rank_metric)) %>%
  filter(!is.na(Gene)) %>%
  distinct(Gene, .keep_all = TRUE) 


#========================= Combine GMT files ==========================#
## This section combines multiple GMT files into a single file for GSEA analysis.

gmt_directory <- "/Users/denisabatrinu/Desktop/Research_Project_25_PCD_PG/MS_analysis_22APR"

gmt_files <- list.files(gmt_directory, pattern = "\\.gmt$", full.names = TRUE)

read_gmt <- function(file) {
  readLines(file)
}

all_gmt_data <- unlist(lapply(gmt_files, read_gmt))

output_file <- "combined_gmt_file.gmt"
writeLines(all_gmt_data, con = output_file)

Ranked_Gene_List_wt = Genes_Ranked_wt$rank_metric
names(Ranked_Gene_List_wt) = Genes_Ranked_wt$Gene

Ranked_Gene_List_wtmut = Genes_Ranked_wtmut$rank_metric
names(Ranked_Gene_List_wt) = Genes_Ranked_wt$Gene

GMT_FILE = "combined_gmt_file.gmt"
GMT_FILE_GENE_SETS = read.gmt(GMT_FILE)

set.seed(123) 
Resulta_GSEA_mut = GSEA(
  geneList = Ranked_Gene_List_mut,
  TERM2GENE = GMT_FILE_GENE_SETS,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  minGSSize = 5,
  maxGSSize = 800,
  verbose = TRUE
)

set.seed(123) 
Resulta_GSEA_wt = GSEA(
  geneList = Ranked_Gene_List_wt,
  TERM2GENE = GMT_FILE_GENE_SETS,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  minGSSize = 5,
  maxGSSize = 800,
  verbose = TRUE
)

set.seed(123) 
Resulta_GSEA_wtmut = GSEA(
  geneList = Ranked_Gene_List_wtmut,
  TERM2GENE = GMT_FILE_GENE_SETS,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  minGSSize = 5,
  maxGSSize = 800,
  verbose = TRUE
)

dim(Resulta_GSEA_wt)
head(Resulta_GSEA_wt)
str(Resulta_GSEA_wt)

# Extract signif pathways:
# Note that only 25 pathways have been able to run from the 41 specified
# presumably due to the proteomic dataset including only ~7000 proteins
Res_wt_v_ctrl_pathways <- Resulta_GSEA_wt@result$Description
class(Res_wt_v_ctrl_pathways)
length(Res_wt_v_ctrl_pathways)
Res_wt_v_ctrl_pathways