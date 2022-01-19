# -------------------------------------------------------------------------
# OG4953 UpSetR Plots -----------------------------------------------------
# Sam Old, 30th September 2021 --------------------------------------------
# -------------------------------------------------------------------------


# Load required packages --------------------------------------------------

library(tidyverse) # tidyverse_1.3.1
library(UpSetR)    # UpSetR_1.4.0


# Load datasets and set initial variables ---------------------------------
# Adapted from an in-house Shiny app used for data interrogation ----------

github        <- "https://raw.githubusercontent.com/SameOldSamOld/RNA-Seq_OG4953/master/data/"
all_data      <- read_csv(paste0(github, "RNASeq_raw_counts_gene_STAR_FR_2019-11-01.csv"))
gene_info     <- read_csv(paste0(github, "gene_info.csv"))
all_pval      <- read_csv(paste0(github, "all_pval.csv"))
all_l2fc      <- read_csv(paste0(github, "all_l2fc.csv"))

colnames(all_data)[1] <- "gene_name"

mean_cutoff   <- 20
l2fc_cutoff   <- 1
pval_cutoff   <- 0.05
biotypeFilter <- "protein_coding"

dataSelected  <- c(
  "KOm24pvsB6m24p",
  "KOm24nvsB6m24n",
  "KOs103pvsB6s103p",
  "KOs103nvsB6s103n",
  "KOe11bpvsB6e11bp",
  "KOeTNvsB6eTN"
)

B6m24p  <- c("B61.m.24p",  "B62.m.24p",  "B63.m.24p" )
B6m24n  <- c("B61.m.24n",                "B63.m.24n" ) # The one sample removed
B6s103p <- c("B61.s.103p", "B62.s.103p", "B63.s.103p")
B6s103n <- c("B61.s.103n", "B62.s.103n", "B63.s.103n")
B6e11bp <- c("B61.e.11bp", "B62.e.11bp", "B63.e.11bp")
B6eTN   <- c("B61.e.TN",   "B62.e.TN",   "B63.e.TN"  )
KOm24p  <- c("KO1.m.24p",  "KO2.m.24p",  "KO3.m.24p" )
KOm24n  <- c("KO1.m.24n",  "KO2.m.24n",  "KO3.m.24n" )
KOs103p <- c("KO1.s.103p", "KO2.s.103p", "KO3.s.103p")
KOs103n <- c("KO1.s.103n", "KO2.s.103n", "KO3.s.103n")
KOe11bp <- c("KO1.e.11bp", "KO2.e.11bp", "KO3.e.11bp")
KOeTN   <- c("KO1.e.TN",   "KO2.e.TN",   "KO3.e.TN"  )


# Final DEGs anad prepare data for UpSetR ---------------------------------

genes_of_biotype <- gene_info %>%
  dplyr::filter(gene_biotype %in% biotypeFilter)

pval_data <- all_pval %>%
  dplyr::filter(gene_name %in% genes_of_biotype$external_gene_name) %>%
  dplyr::select(all_of(dataSelected))

l2fc_data <- all_l2fc %>%
  dplyr::filter(gene_name %in% genes_of_biotype$external_gene_name) %>%
  dplyr::select(gene_name, all_of(dataSelected))

## apply mean count cutoff
# Need to filter counts PER set in a loop -----------------------------

for (i in dataSelected) {
  ss <- unlist(strsplit(i, "vs"))
  mean_pass <- all_data %>%
    dplyr::select(gene_name, get(ss[1]), get(ss[2])) %>%
    dplyr::filter(rowMeans(dplyr::select(., c(get(ss[1]),get(ss[2])))) >= mean_cutoff) %>%
    dplyr::pull(gene_name)
  l2fc_data[!l2fc_data$gene_name %in% mean_pass, i] <- 0
}

l2fc_data <- dplyr::select(l2fc_data, -gene_name)

pval_data_cut <- pval_data < pval_cutoff
l2fc_data_neg <- l2fc_data < -l2fc_cutoff
l2fc_data_pos <- l2fc_data >  l2fc_cutoff
l2fc_data_com <- l2fc_data_neg + l2fc_data_pos

final_data    <- l2fc_data_com + pval_data_cut
final_data    <- final_data == 2
final_data    <- final_data[rowSums(final_data) > 0,]
final_data    <- as.data.frame(final_data + 1 - 1)



# Plot Supplemental Figure 2.h --------------------------------------------

upset(data = final_data,
      nintersects = NA,
      nsets = length(dataSelected),
      sets = dataSelected,
      point.size = 5,
      group.by = "degree",
      matrix.color = "black",
      main.bar.color = "black",
      sets.bar.color = "black",
      shade.color = "black",
      text.scale = 1.6)
