
# -------------------------------------------------------------------------
# OG4953 Principal Component Analysis Plots -------------------------------
# Sam Old, 30th September 2021 --------------------------------------------
# -------------------------------------------------------------------------


# Load required packages --------------------------------------------------

library(biomaRt)   # biomaRt_2.46.3
library(DESeq2)    # DESeq2_1.30.1
library(tidyverse) # tidyverse_1.3.1



# -------------------------------------------------------------------------
# Load required datasets --------------------------------------------------

metadata <- read.csv(
  file = "https://raw.githubusercontent.com/SameOldSamOld/RNA-Seq_OG4953/master/data/metadata.csv",
  row.names = 1, header = TRUE)
metadata$Sample.name <- gsub("-", ".", metadata$Sample.name)

raw_counts <- read.csv(
  file = "https://raw.githubusercontent.com/SameOldSamOld/RNA-Seq_OG4953/master/data/RNASeq_raw_counts_gene_STAR_FR_2019-11-01.csv",
  row.names = 1, header = TRUE)



# -------------------------------------------------------------------------
# Convert raw_counts into vstpk_counts ------------------------------------
# Estimated run time: 5 minutes -------------------------------------------

mart <- useEnsembl(biomart = "ensembl",
                   dataset = "mmusculus_gene_ensembl",
                   version = 91)
transcripts <- getBM(attributes = c("external_gene_name", "ensembl_gene_id_version", "transcript_length", "description", "gene_biotype"),
                     filters = "external_gene_name",
                     values = rownames(raw_counts),
                     bmHeader = FALSE,
                     mart = mart)

# Check correct version of ensembl database was used: should return all TRUE
table(rownames(raw_counts) %in% transcripts$external_gene_name)

# Keep only the longest transcript for VSTpk gene normalisation
longest_transcripts <- transcripts %>%
  dplyr::group_by(external_gene_name) %>%
  dplyr::slice(which.max(transcript_length)) %>%
  dplyr::arrange(external_gene_name)

# Prepare VST Matrix
counts <- raw_counts %>%
  arrange(rownames(.)) %>%
  dplyr::select(metadata$Sample.name)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = metadata,
                              ~ Marker) %>%
  DESeq()
dds.vst <- varianceStabilizingTransformation(dds, fitType = "local")
dds.vst.mat <- assays(dds.vst)[[1]];

# log rules to calculate VSTpk normalisation. log1p is used to estimate log2 while accounting for zero counts
vst.logklengths <- log1p(longest_transcripts[,"transcript_length"]/1000)/log1p(2);
rownames(vst.logklengths) <- longest_transcripts$external_gene_name

# Create and fill final VSTpk matrix
dds.vstpk.mat <- matrix( data = NA, nrow = nrow(dds.vst.mat), ncol = ncol(dds.vst.mat), dimnames = dimnames(dds.vst.mat))
for (i in rownames(dds.vstpk.mat)) {
  dds.vstpk.mat[i,] <- dds.vst.mat[i,] - vst.logklengths[i,]
}

vstpk_counts <- as_tibble(dds.vstpk.mat) %>%
  add_column(gene_name = rownames(dds.vstpk.mat)) %>%
  dplyr::select(gene_name, all_of(colnames(dds.vstpk.mat))) %>%
  arrange(gene_name)




# -------------------------------------------------------------------------
# Plot PCA function -------------------------------------------------------

plotPCA <- function(z = NULL, importance = NULL, pcx = NULL, pcy = NULL, topGenes = NULL) {

  name1 = paste0("PC", pcx)
  name2 = paste0("PC", pcy)

  dat <- z %>%
    dplyr::rename(!!name1 := x, !!name2 := y) %>%
    tibble::as_tibble()
  dat

  pca_labs <- c(
    "103n-s" = bquote(paste("CD103"^"-")),
    "103p-s" = bquote(paste("CD103"^"+")),
    "11bp-e" = bquote(paste("CD11b"^"hi")),
    "TN-e"   = bquote(paste("CD11b"^"lo")),
    "24n-m"  = bquote(paste("CD24"^"-")),
    "24p-m"  = bquote(paste("CD24"^"+"))
  )

  ggplot(dat, aes(x = get(name1), y = get(name2), shape = Cell, fill = Cell, alpha = Immunisation)) +
    geom_hline(yintercept = 0, colour = "grey80") +
    geom_vline(xintercept = 0, colour = "grey80") +
    geom_point(size = 4) +
    scale_alpha_discrete(name = "Mouse Strain",
                         range = c(1, 0.45)) +
    scale_fill_manual(name = "Cell Subsets",
                      labels = pca_labs,
                      values = c("103n-s" = "#00BFC4",
                                 "103p-s" = "#F8766D",
                                 "11bp-e" = "#619CFF",
                                 "TN-e"   = "#B79F00",
                                 "24n-m"  = "#F564E3",
                                 "24p-m"  = "#00BA38")) +
    scale_shape_manual(name = "Cell Subsets",
                       labels = pca_labs,
                       values = c("103n-s" = 23,
                                  "103p-s" = 23,
                                  "11bp-e" = 22,
                                  "TN-e"   = 22,
                                  "24n-m"  = 21,
                                  "24p-m"  = 21)) +
    labs(title = paste("Principal Component Analysis -", topGenes, "most variable genes"),
         x = paste0(name1, " - ", importance[pcx], "%"),
         y = paste0(name2, " - ", importance[pcy], "%")) +
    theme_bw()
}


# Perform PCA Function ----------------------------------------------------

performPCA <- function(pcx = 1, pcy = 2, columns = "", vstpk_data = NULL, metadata = NULL) {

  data <- vstpk_data %>%
    dplyr::select_if(is_double) %>%
    t()

  # Organise PCA data
  data <- data[columns, order(apply(data,2,var), decreasing = T)]            # Order PCA data for subsetting
  data <- data[,apply(data, 2, var, na.rm = TRUE) != 0]                      # remove 0 Variation, NA rows

  # Compute PCA
  data.pca <- prcomp(data, center = T)

  # Generate PCA data
  x <- as.numeric(data.pca$x[,pcx])
  y <- as.numeric(data.pca$x[,pcy])

  z <- as.data.frame(cbind(x,y), row.names = columns)

  ## CUSTOM metadata columns
  meta <- metadata %>%
    dplyr::filter(Sample.name %in% rownames(z)) %>%                             # filter only selected rows
    slice(match(rownames(z), Sample.name))                                # order by pca data rownames
  z$Cell         <- as.factor(paste(meta$Cell, meta$Tissue, sep = "-"))
  z$Immunisation <- as.factor(meta$B6)
  z$label        <- as.character(meta$Marker)
  z$Marker       <- as.character(as.factor(meta$Sample.name))
  importance     <- as.numeric(summary(data.pca)$importance["Proportion of Variance",])*100

  return(list(z,importance))
}




# -------------------------------------------------------------------------
# Plotting of PCA figures -------------------------------------------------

# Figure 1.f. -------------------------------------------------------------


# Load Variables

selected_columns   <- metadata %>%
  dplyr::filter(`B6` == "B6") %>%
  dplyr::filter(!`Sample.name` %in% "B62.m.24n") %>%  # Outlier sample
  dplyr::pull(`Sample.name`)

pca <- performPCA(
  pcx = 1, pcy = 2,
  columns = selected_columns,
  vstpk_data = vstpk_counts,
  metadata = metadata)

plotPCA(z = pca[[1]], importance = pca[[2]], pcx = 1, pcy = 2, topGenes = ncol(data))

pca <- performPCA(
  pcx = 1, pcy = 3,
  columns = selected_columns,
  vstpk_data = vstpk_counts,
  metadata = metadata)

plotPCA(z = pca[[1]], importance = pca[[2]], pcx = 1, pcy = 3, topGenes = ncol(data))


# Supplemental Figure 2.f -------------------------------------------------

# S2.f: Ear draining lymph node Sample PCA
selected_columns   <- metadata %>%
  dplyr::filter(`Tissue` == "e") %>%
  dplyr::pull(`Sample.name`)

pca <- performPCA(
  pcx = 1, pcy = 2,
  columns = selected_columns,
  vstpk_data = vstpk_counts,
  metadata = metadata)

plotPCA(z = pca[[1]], importance = pca[[2]], pcx = 1, pcy = 2, topGenes = ncol(data))

pca <- performPCA(
  pcx = 1, pcy = 3,
  columns = selected_columns,
  vstpk_data = vstpk_counts,
  metadata = metadata)

plotPCA(z = pca[[1]], importance = pca[[2]], pcx = 1, pcy = 3, topGenes = ncol(data))


# S2.f: Mediastinal lymph node Sample PCA
selected_columns   <- metadata %>%
  dplyr::filter(`Tissue` == "m") %>%
  dplyr::filter(!`Sample.name` %in% "B62.m.24n") %>%  # Outlier sample
  dplyr::pull(`Sample.name`)

pca <- performPCA(
  pcx = 1, pcy = 2,
  columns = selected_columns,
  vstpk_data = vstpk_counts,
  metadata = metadata)

plotPCA(z = pca[[1]], importance = pca[[2]], pcx = 1, pcy = 2, topGenes = ncol(data))


# S2.f: Skin Sample PCA
selected_columns   <- metadata %>%
  dplyr::filter(`Tissue` == "s") %>%
  dplyr::pull(`Sample.name`)

pca <- performPCA(
  pcx = 1, pcy = 2,
  columns = selected_columns,
  vstpk_data = vstpk_counts,
  metadata = metadata)

plotPCA(z = pca[[1]], importance = pca[[2]], pcx = 1, pcy = 2, topGenes = ncol(data))
