#Comparison of proteins present in the four samples

#clear workspace
rm(list=ls())


#load libraries
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(writexl)

#loaddataset
DE_WUDEA_PASNvsTSB1 <- read_excel("~/BIO253 Correlation 6850 v JE2/DE_WUDEA_PASNvsTSB1.xlsx",
                                  sheet = "diff_exp_analysis")
JE2_TSBvPASN_data <- read_excel("~/BIO253 Correlation 6850 v JE2/JE2_TSBvPASN_data.xlsx")
SAstrainSpecificIDs_to_Uniprot_2_ <- read_excel("~/BIO253 Correlation 6850 v JE2/SAstrainSpecificIDs_to_Uniprot (2).xlsx",
                                                sheet = "Sheet1", skip = 6)
SAstrainSpecificIDs_to_Uniprot_onlyJE2 <- read_excel("~/BIO253 Correlation 6850 v JE2/SAstrainSpecificIDs_to_Uniprot_onlyJE2.xlsx",
                                                     sheet = "Sheet1", skip = 6)
X6850_2020_data <- read_excel("~/BIO253 Correlation 6850 v JE2/6850_2020_data.xlsx",
                              skip = 1)
JE2_2020_data <- read_excel("~/BIO253 Correlation 6850 v JE2/JE2_2020_data.xlsx",
                            skip = 2)


SA6850_2024 <- DE_WUDEA_PASNvsTSB1
SA6850_2020 <- X6850_2020_data
JE2_2024 <- JE2_TSBvPASN_data
JE2_2020 <- JE2_2020_data
Uniprot_all <- SAstrainSpecificIDs_to_Uniprot_2_
Uniprot_JE2 <- SAstrainSpecificIDs_to_Uniprot_onlyJE2


rm(SAstrainSpecificIDs_to_Uniprot_2_)
rm(SAstrainSpecificIDs_to_Uniprot_onlyJE2)
rm(DE_WUDEA_PASNvsTSB1)
rm(JE2_TSBvPASN_data)
rm(X6850_2020_data)
rm(JE2_2020_data)


#Add Uniprot ID to all (new column with same name for all)
JE2_2020$uniprot_ID <- JE2_2020$SAuniprotID
SA6850_2020$uniprot_ID <- SA6850_2020$SAuniprotID
SA6850_2024 <- SA6850_2024 %>%
    mutate(locus_tag = str_extract(description, "(?<=\\[locus_tag=)[^\\]]+"))
JE2_2024 <- JE2_2024 %>%
    mutate(locus_tag = str_extract(description, "(?<=\\[locus_tag=)[^\\]]+"))
Uniprot_all <- Uniprot_all %>%
    mutate(uniprot_ID = str_extract(BlastOrthologue, "(?<=\\|)[^\\|]+"))
SA6850_2024 <- left_join(SA6850_2024, Uniprot_all[, c("myLocTag", "uniprot_ID")],
    by = c("locus_tag" = "myLocTag"))
JE2_2024 <- left_join(JE2_2024, Uniprot_all[, c("myLocTag", "uniprot_ID")],
                         by = c("locus_tag" = "myLocTag"))

#make venn diagram - ChatGPT Code
############################################################
# Packages
############################################################
#install.packages("ggVennDiagram")   # only once if not installed
library(ggVennDiagram)
library(dplyr)

############################################################
# Build sets: proteins present per strain/year
############################################################

set_JE2_2020    <- JE2_2020$uniprot_ID    %>% unique() %>% na.omit()
set_SA6850_2020 <- SA6850_2020$uniprot_ID %>% unique() %>% na.omit()
set_SA6850_2024 <- SA6850_2024$uniprot_ID %>% unique() %>% na.omit()
set_JE2_2024    <- JE2_2024$uniprot_ID    %>% unique() %>% na.omit()

protein_sets <- list(
    `JE2 2020`     = set_JE2_2020,
    `6850 2020`    = set_SA6850_2020,
    `6850 2024`    = set_SA6850_2024,
    `JE2 2024`     = set_JE2_2024
)

############################################################
# Venn diagram (ggplot-based)
############################################################

p_venn <- ggVennDiagram(protein_sets) +
    scale_fill_gradient(
        low  = "#deebf7",
        high = "#08519c",
        limits = c(0, 400),   # <-- THIS caps the color scale
        oob = scales::squish   # <-- values >250 get squished to darkest blue
    ) +
    theme_void() +
    theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    ggtitle("Overlap of detected proteins (Uniprot IDs)")


p_venn

#Extract all uniprotIDs to use as background in String
uniprot_ID_all <- unique(c(JE2_2020$uniprot_ID, JE2_2024$uniprot_ID, SA6850_2020$uniprot_ID, SA6850_2024$uniprot_ID))
write.table(uniprot_ID_all, "allUniprotID_unique.txt", quote = FALSE, row.names = FALSE)
getwd()

############################################################
############################################################
# Filter by FDR < 5 and extract Uniprot IDs
# (if your FDR is 0–1, change 5 to 0.05)
############################################################

JE2_2020_FDR    <- JE2_2020    %>% filter(FDR < 0.05)
SA6850_2020_FDR <- SA6850_2020 %>% filter(FDR < 0.05)
SA6850_2024_FDR <- SA6850_2024 %>% filter(FDR < 0.05)
JE2_2024_FDR    <- JE2_2024    %>% filter(FDR < 0.05)

set_JE2_2020    <- JE2_2020_FDR$uniprot_ID    %>% unique() %>% na.omit()
set_SA6850_2020 <- SA6850_2020_FDR$uniprot_ID %>% unique() %>% na.omit()
set_SA6850_2024 <- SA6850_2024_FDR$uniprot_ID %>% unique() %>% na.omit()
set_JE2_2024    <- JE2_2024_FDR$uniprot_ID    %>% unique() %>% na.omit()

protein_sets_FDR <- list(
    `JE2 2020 (FDR<5)`     = set_JE2_2020,
    `6850 2020 (FDR<5)`    = set_SA6850_2020,
    `6850 2024 (FDR<5)`    = set_SA6850_2024,
    `JE2 2024 (FDR<5)`     = set_JE2_2024
)

############################################################
# Venn diagram with capped blue scale (max colour at 250)
############################################################

p_venn_FDR <- ggVennDiagram(protein_sets_FDR, label_alpha = 0) +
    scale_fill_gradient(
        low  = "#deebf7",
        high = "#08519c",
        limits = c(0, 400),          # cap colour scale at 250
        oob = scales::squish         # values >250 use darkest blue
    ) +
    theme_void() +
    theme(
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16)
    ) +
    ggtitle("Overlap of detected proteins (Uniprot IDs, FDR < 0.05)")

p_venn_FDR
############################################################
# SAVE BOTH PLOTS WITH ggsave()
############################################################

outdir <- "C:/Users/Lea Sonderegger/Documents/BIO253 Correlation Output/Proteins detected"

# Full plot
ggsave(
    filename = file.path(outdir, "Venn_Proteins_all.png"),
    plot = p_venn,
    width = 10,
    height = 8,
    dpi = 300
)

# FDR < 5 plot
ggsave(
    filename = file.path(outdir, "Venn_Proteins_FDR_less_0.05.png"),
    plot = p_venn_FDR,
    width = 10,
    height = 8,
    dpi = 300
)
############################################################
