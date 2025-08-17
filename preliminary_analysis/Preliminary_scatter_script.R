########## Install packages ########################################

install.packages(c("readxl", "dplyr", "ggplot2", "ggrepel"))
library(readxl)
library(dplyr)    
library(ggplot2)
library(ggrepel)



############### load and edit data #####################################
data <- read_excel("/Users/denisabatrinu/Desktop/DB_thesis_2025_MSc_HMG/Mass_Spec_raw_data/raw.xlsx", sheet = 2)
# CAREFUL! data contains numbers written with commas as decimal separators instead of periods so readxl reads that as text...

# Fix warnings:
# Identify character columns.
char_cols <- names(data)[sapply(data, is.character)]
char_cols
str(data)

# Columns to convert
cols_to_convert <- c("mut_wt_sigA", "mut_ctrl_sigA", "tstatistic_mut_ctrl_lung", "minus_log_pval_ttest_mut_ctrl_lung")

# Conversion
data[cols_to_convert] <- lapply(data[cols_to_convert], function(x) {
  as.numeric(gsub(",", ".", x))
})

str(data)

# create new column that is TRUE only when protein has "+" in "Ttest_wt_over_ct" and "Asig_wt_over_ct" are <0.05  
data$signif_highlight <- (data$ttest_sig_mut_wt_lung == "+" & data$asig_mut_wt_lung == "+")


# create gene list
categories <- list(
  "PAQosome network" = c("Ruvbl1", "Ruvbl2", "Rpap3", "Pih1d1", "Hsp90", "Uri1", "Uxt", "Pdrg1", "Pfdn2", "Pfdn6", "Asdurf", 
                         "Wdr92", "Polr2e", "Rpb5", "Telo2", "Tti1", "Tti2", "Nufip1", "Znhit3", "Znhit6", "Znhit2", "Ecd", 
                         "Aar2", "Naf1", "Shq1", "Nopchap1", "Hop", "p23"),
  "PCD genes" = c("Ccdc103", "Ccdc114", "Ccdc151", "Ccdc164", "Ccdc40", "Ccdc65", "Ccno", "Cfap54", "Cfap57", "Cfap74", 
                  "Cfap221", "Cfap298", "Cfap300", "Clxn", "Dnaaf1", "Dnaaf2", "Dnaaf3", "Dnaaf4", "Dnaaf5", "Dnaaf6", 
                  "Dnaaf11", "Dnah1", "Dnah5", "Dnah7", "Dnah9", "Dnah10", "Dnah11", "Dnai1", "Dnai2", "Dnajb13", 
                  "Dnal1", "Drc1", "Foxj1", "Gas2l2", "Gas8", "Hydin", "Kiaa0586", "Lrrc56", "Mcidas", "Nek10", 
                  "Nme5", "Nme8", "Odad1", "Odad2", "Odad3", "Odad4", "Ofd1", "Rpgr", "Rsph1", "Rsph3", "Rsph4a", "Rsph9", 
                  "Spag1","Spef2", "Stk36", "Tp73", "Ttc12", "Zmynd10"),
  "IFT" = c("Ift20", "Tubb6", "Ift46", "Ift43", "Ift22", "Dynll1", "Tnpo1", "Dynlt2b", "Cluap1", "Tubb1", "Ift27", "Ift52",
            "Kif3a", "Dynll2", "Kif3c", "Trip11", "Tuba1b", "Ift140", "Dync2li1", "Tuba4a", "Kifap3", "Kif3b", "Ift80", "Dynlt5",
            "Ift74", "Kif17", "Ift81", "Ift122", "Ift57", "Traf3ip1", "Dynlrb2", "Ttc21b", "Tubb4b", "Wdr19", "Ift172", "Dync2i2", 
            "Dync2i1", "Tuba1c", "Tubb2b", "Dync2h1", "Dynlrb1", "Ift56", "Tubb2a", "Tubb3", "Tubb4a", "Ift25", "Wdr35", "Tuba3b", 
            "Tuba3a", "Tuba1a", "Ift70a1", "Ift70a2", "Ift70b", "Dynlt2a1", "Dynlt2a3")
)


data <- data %>%
  mutate(category = case_when(
    genes_lung %in% categories[["PAQosome network"]] ~ "PAQosome network",
    genes_lung %in% categories[["PCD genes"]] ~ "PCD genes",
    genes_lung %in% categories[["IFT"]] ~ "IFT",
    signif_highlight ~ "Other significant",
    TRUE ~ "Not significant"
  ))


# Define min/max for shaded regions
x_min <- min(data$mut_wt_lung, na.rm = TRUE)
x_max <- max(data$mut_wt_lung, na.rm = TRUE)
x_mid <- 0


# Scatter plot with shaded regions
ggplot(data, aes(x = mut_wt_lung, y = sum_intensity_lung)) +
  annotate("rect", xmin = x_mid, xmax = x_max, ymin = -Inf, ymax = Inf,
           fill = "#F7C6B1", alpha = 0.2) +
  annotate("rect", xmin = x_min, xmax = x_mid, ymin = -Inf, ymax = Inf,
           fill = "lightblue", alpha = 0.2) +
  geom_point(aes(color = category), size = 1, alpha = 0.8) +
  scale_color_manual(values = c(
    "PAQosome network" = "#7F6244",
    "PCD genes" = "black",
    "IFT" = "red",
    "Other significant" = "lightpink",
    "Not significant" = "grey70"
  )) +
  geom_text_repel(
    data = subset(data, signif_highlight & genes_lung %in% unlist(categories)),
    aes(label = genes_lung, color = category),
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = "grey50",
    show.legend = FALSE
  ) +
  annotate("text", x = x_max, y = max(data$sum_intensity_lung, na.rm = TRUE),
           label = "Mutant Dnaaf19 enriched", color = "red",
           hjust = 1, vjust = 1.2, size = 5) +
  annotate("text", x = x_min, y = max(data$sum_intensity_lung, na.rm = TRUE),
           label = "Wildtype Dnaaf19 enriched", color = "blue4",
           hjust = 0, vjust = 1.2, size = 5) +
  theme_minimal() +
  labs(
    x = "log2 median MUT/WT",
    y = "Intensity",
    color = "Category"
  ) +
  theme(
    legend.position = "right",
    text = element_text(size = 12)
  )

