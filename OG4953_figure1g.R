### Sam Old
### 11th August 2020
### Create some figures for STAT6KO publication


# Load Packages and Data --------------------------------------------------

library(DESeq2)       # v1.28.1
library(ggdendro)     # v0.1.21
library(ggplot2)      # v3.3.2
library(ggpubr)       # v0.4.0
library(gridExtra)    # v2.3
library(pheatmap)     # v1.0.12
library(RColorBrewer) # v1.1-2
library(tidyverse)    # v1.3.0

raw_counts   <- read.csv("/Users/siold/Proj/FR_RNA-Seq/data/RNASeq_raw_counts_gene_STAR_FR_2019-11-01.csv",
                       row.names = 1, header = T)
vstpk_counts <- read_csv('/Users/siold/Proj/FR_RNA-Seq/shiny_august2020/FR19/data/vstpk_data.csv')
metadata     <- read_csv('/Users/siold/Proj/FR_RNA-Seq/shiny_august2020/FR19/data/metadata.csv')


# Load gene lists ---------------------------------------------------------

## New Gene Lists
# CD103- and CD24- SIGNATURE
group1_sigraw <- "Sirpb1a Sirpb1b Sirpb1c Tmem132c Plaur Zfp521 Ryr3 Vegfa Rgl1 Esrrg Usp2 H2-Ob Cyp4f16 Myo1f Bmp2 Dab2ip Lrrc18 Lag3 Angpt1 Nid1 Wdfy4 Clec9a Pde7b Plbd1 Sema4d Il1r2 Tm6sf1 Ncald Ccl4 Havcr2 Cass4 Slamf8 Plk3 Tns3 Tspan9 Cd4 Eef2k Tlr7 Lrp8 Rnf26 Alox5ap Aspm Ms4a4c Sla Mef2c Dtx1 Tg Tmem156 Atp8b4 Csf1r Dpp10 Elmsan1 Mapk14 Tpd52 Fshr Pik3ap1 Wdfy2 Cyth4 March1 Itgal Fmnl3 Erdr1 Rora Tbc1d30 Egr2 Fcrl1 Ifitm3 Egr3 Tgm2 Maml3 Stk38l Usp6nl H60b Rin2 Slc16a6 Dnal1 Tle1 Zfp36l2 Gpr35 Vmn2r90 Tagap Arhgap25 Evi2b Btla Cln8 Raet1e Ly6e Ms4a6b Kmo Trit1 Arhgap15 Chn2 Best1 A530032D15Rik Ms4a6c Stap1 9930111J21Rik2 9930111J21Rik1 Ifi207 Bach2 Rem1"
# CD103+ SIGNATURE
group2_sigraw <- "Fndc5 Car2 Pla2g7 Rnase6 Syngr3 Cmtm8 Ppp1r1a Ido1 Cfb Ccnd2 Ptk6 Gm20547 Gpr55 Arl14 Mpzl2 Cd1d1 4933407L21Rik Dnah2 Palm Pygl Zc3h12d Glycam1 Ahrr Apol10b Mylk Plb1 Slc9a3r2 Epcam Exoc3l4 B3gnt3 Glrx Tgtp1 Pcbp4 Pbld2 Gpr62 Tmem86a Nme9 Tpm2 Slc43a3 Apol7e Tm4sf5 Gm12216 Rasl11a Ulk3 Apol7b 6430531B16Rik Large1 Esam Insl6 Rgs3 Mpi Usp18 Apol7c Zbtb12 Mospd1 Coro1b Tmem140 Nudt22 Cmc2 Sept12 Lima1 Tm2d2 Apol7a Nelfe Uhrf1bp1 3110062M04Rik Ccdc120 P3h1 Mtmr4 Tap2 Gm15821 Atmin Hsh2d Akt1s1 Rragd Cwf19l1 Uba7 Ocln Cptp Fam234b Myom1 Lad1 Gm20661 Tmem8 Tbc1d17 Tnnt2"
# CD24+ and CD11b+ SIGNATURE
group3a_sigraw <- "Ephx1 S100a10 P2ry1 Cyp2ab1 Sod3 P3h2 Slc27a3 Dkkl1 Serpinb1a Capg Tmem159 Gm2115 Apbb2 Pros1 Phlda3 Fxyd3 Mt2 Ccl17 Lgalsl Ehf Cavin1 Cxcl1 Gpr85 5730409E04Rik Pdlim1 Itgb8 Anxa1 Plppr4 Tsku Tpi1 Mt1 Lgi4 Mfge8 Macc1 Ddias Cpt1c Slamf9"
# TN SIGNATURE
group3b_sigraw <- "Pi16 Myrip Khdc1c Efna5 Khdc1a Sned1 Lamc2 Sim1 Gm5796 Gm3468 Eya1 Ptprf Lhfp Cd8a Gm3488 Il9r Gm3264 Gm3739 Gm3411 Gm3500 Gm3373 Gm8206 Arhgap29 Gm3252 Klf12 Ube2ql1 Gm3383 Gm3558 Gm3512 Gm2956 Gm3696 Ranbp17 Gm3194 Mterf4 Scd1 Gm8281 Prox2 Calm4 Clmn Gm3164 Gm2897 Dennd5b Gm3636 Gm3020 Gm10409 Gm3173 Tlx2 Gm3002 Gm3752 Gm10406 Dclk1 Tgfbr3 Gm8108 Gm3667 Gm3005 2610528A11Rik Trim43b Myb Gm2974 Gm2237 Cdhr1 Antxr1 Rhbdf1 Snrnp25 Snap91 Bmp1 1700025B11Rik Shisa9 Msrb3 Izumo4 Tnfsf4 Clcn5 Slc41a2 Neto2 Jcad Fndc10 Ankrd42 Pde5a S1pr1"

to_vector <- function(group = NULL) return(unique(unlist(strsplit(group, " "))))
# JOHANNES REQUESTEDL: swap of group1 and group 2 as labels
group1_sig  <- sort(to_vector(group1_sigraw))
group2_sig  <- sort(to_vector(group2_sigraw))
group3a_sig <- sort(to_vector(group3a_sigraw))
group3b_sig <- sort(to_vector(group3b_sigraw))

genes <- unique(c(group1_sig, group2_sig, group3a_sig, group3b_sig))
length(genes) == sum(length(group1_sig),length(group2_sig), length(group3a_sig), length(group3b_sig))



# Filter samples, and perform DESeq2 normalisation --------------------------------------------

# Filter raw counts on metadata to remove sample #14
counts <- raw_counts %>%
  dplyr::select(metadata$Sample.name) %>%
  dplyr::arrange(rownames(.))

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = metadata,
                              ~ Tissue) %>%
  DESeq()

all_normalised <- dds %>%
  counts(normalized = TRUE)

counts.cut <- all_normalised %>%
  as.data.frame() %>%
  dplyr::select(metadata %>%
                  dplyr::filter(B6 == "B6") %>%
                  dplyr::pull(Sample.name))


# Final Dendrogram --------------------------------------------------------

counts.cut.log10 <- log10(counts.cut + 1)
phm              <- pheatmap(mat = counts.cut.log10[genes,],
                             cluster_rows = FALSE,
                             cluster_cols = TRUE)

hc   <- phm$tree_col
dend <- as.dendrogram(hc)


# Twist dendrogram branches to group technical replicates (vis) -----------

order.dendrogram(dend)
neworder  <- as.integer(c(11:9,5:4,8:6,3:1,14:12,17:15))
dend.2    <- reorder(dend,
                     wts = order(match(colnames(counts.cut[,neworder]), colnames(counts.cut))), mean)
dend_data <- dendro_data(dend.2, type = "rectangle")

classify      <- dend_data$labels;
classify$labs <- dend_data$labels$label;
classify$labs <- sub('.*\\.', '', classify$labs)
classify$cell <- c("103","103","103","24","24","103","103","103","24","24","24","11b","11b","11b","11b","11b","11b")

# Format labels for final figure
dend_data$labels$final_label <- c(
  "CD103^-{}*1","CD103^-{}*3","CD103^-{}*2",
  "CD24^-{}*3", "CD24^-{}*1",
  "CD103^+{}*2","CD103^+{}*3","CD103^+{}*1",
  "CD24^+{}*3","CD24^+{}*2","CD24^+{}*1",
  "CD11b^hi* 3","CD11b^hi* 2","CD11b^hi* 1",
  "CD11b^lo* 3","CD11b^lo* 2","CD11b^lo* 1")

gg <- list()
gg[[1]] <- ggplot(segment(dend_data)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  # geom_text(data = dend_data$labels, aes(x, y, label = final_label),
  #           check_overlap = TRUE, hjust = 0, angle = 0, size = 3) +
  annotate("text", x = dend_data$labels$x, y = dend_data$labels$y - 20,
           label = dend_data$labels$final_label, parse = TRUE) +
  geom_point(data = classify, aes(x = x, y = y - 4, fill = labs, shape = cell),
             color = "black", size = 4.3) +
  scale_fill_manual(name = "Cell Subsets",
                    # labels = pca_labs,
                    values = c("#00BFC4","#F8766D","#619CFF","#F564E3","#00BA38","#B79F00")) +
  scale_shape_manual(values = c(23, 22, 21)) +
  ylim(25, 100) +
  coord_flip() +
  theme_classic() +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.line = element_blank(),
        legend.position = "none") +
  scale_y_reverse(expand = c(0.2, 0)) +
  ylab("") +
  xlab("")
gg[[1]]
ggsave(
  filename = "B6_dendrogram_leaf-rearranged_5x8.pdf",
  path = "/Users/siold/Proj/FR_RNA-Seq/plots/dendrogram/",
  plot = gg[[1]],
  device = "pdf",
  width = 5,
  height = 8,
  units = "in",
  dpi = 300
)
ggsave(
  filename = "B6_dendrogram_leaf-rearranged_5x8.svg",
  path = "/Users/siold/Proj/FR_RNA-Seq/plots/dendrogram/",
  plot = gg[[1]],
  device = "svg",
  width = 5,
  height = 8,
  units = "in",
  dpi = 300
)

ggsave(
  filename = "B6_dendrogram_leaf-rearranged_5x8_HIGH-DPI.pdf",
  path = "/Users/siold/Proj/FR_RNA-Seq/plots/dendrogram/",
  plot = gg[[1]],
  device = "pdf",
  width = 5,
  height = 8,
  units = "in",
  dpi = 320
)
ggsave(
  filename = "B6_dendrogram_leaf-rearranged_5x8_HIGH-DPI.svg",
  path = "/Users/siold/Proj/FR_RNA-Seq/plots/dendrogram/",
  plot = gg[[1]],
  device = "svg",
  width = 5,
  height = 8,
  units = "in",
  dpi = 320
)

# Figure 1. Put the heatmap together with z score scaling  ----------------

vstpk_selected <- vstpk_counts %>%
  dplyr::filter(gene_name %in% unique(genes)) %>%
  dplyr::select(gene_name, starts_with("B6")) %>%
  as.data.frame() #%>%
rownames(vstpk_selected) <- vstpk_selected$gene_name
vstpk_selected <- vstpk_selected %>%
  dplyr::select(starts_with("B6"))


order <- c("B61.s.103p","B62.s.103p","B63.s.103p","B61.s.103n","B62.s.103n","B63.s.103n",
           "B61.m.24n", "B63.m.24n","B61.m.24p" , "B62.m.24p" , "B63.m.24p",
           "B61.e.11bp","B62.e.11bp","B63.e.11bp","B61.e.TN","B62.e.TN","B63.e.TN")
res_selected <- vstpk_selected[unique(c(group1_sig, group2_sig, rev(group3a_sig), group3b_sig)), order]
final_dat <- t(scale(t(res_selected),  center = TRUE))


# annotation_row <- matrix(data = c(genes, rep("", length(genes))), ncol = length(genes), byrow = TRUE)
annotation_row <- data.frame(Group = factor(c(rep("1   (101)", length(group1_sig)),
                                              rep("2   (86)", length(group2_sig)),
                                              rep("3a (37)",length(group3a_sig)),
                                              # rep("3 overlap (36)",length(group3_over)),
                                              rep("3b (78)",length(group3b_sig)))))
rownames(annotation_row) <- genes
# skip to line 111 to get classify label
labels <- as.character(classify$label)

pheatmap::pheatmap(final_dat[,labels],
                   color = viridisLite::viridis(512),
                   cluster_cols = FALSE, cluster_rows = FALSE,
                   show_rownames = FALSE, annotation_row = annotation_row,
                   # gaps_row = c(59,272,360, 396), gaps_col = c(6,11),
                   main = paste("Scaled VSTpk expression,", length(genes), "genes, B6 mice"),
                   width = 6, height = 5)

pheatmap::pheatmap(final_dat[,labels],
                   color = viridisLite::viridis(512),
                   cluster_cols = FALSE, cluster_rows = FALSE,
                   show_rownames = FALSE, annotation_row = annotation_row,
                   # gaps_row = c(59,272,360, 396), gaps_col = c(6,11),
                   main = paste("Scaled VSTpk expression,", length(genes), "genes, B6 mice"),
                   width = 6, height = 5,
                   filename = "/Users/siold/Proj/FR_RNA-Seq/plots/heatmap/6x5_heatmap.pdf")

pheatmap::pheatmap(final_dat[,labels],
                   color = viridisLite::viridis(512),
                   cluster_cols = FALSE, cluster_rows = FALSE,
                   show_rownames = FALSE, annotation_row = annotation_row,
                   # gaps_row = c(59,272,360, 396), gaps_col = c(6,11),
                   main = paste("Scaled VSTpk expression,", length(genes), "genes, B6 mice"),
                   width = 8, height = 8,
                   filename = "/Users/siold/Proj/FR_RNA-Seq/plots/heatmap/8x8_heatmap.pdf")



# Format the plots together -----------------------------------------------

gg[[1]]
x <- pheatmap::pheatmap(t(final_dat[,rev(labels)]),
                        color = viridisLite::viridis(512),
                        cluster_cols = FALSE, cluster_rows = FALSE,
                        show_colnames = FALSE, show_rownames = FALSE,
                        annotation_col = annotation_row,
                        cellheight = 30,
                        # gaps_row = c(59,272,360, 396), gaps_col = c(6,11),
                        main = paste("Scaled VSTpk expression",length(rownames(final_dat)),"genes, B6 mice"))
gg[[2]] <- x[[4]]
dev.off()
lay <- rbind(c(NA,NA,2, 2, 2),
             c(1, 1, 2, 2, 2),
             c(1, 1, 2, 2, 2),
             c(1, 1, 2, 2, 2),
             c(1, 1, 2, 2, 2),
             c(1, 1, 2, 2, 2),
             c(1, 1, 2, 2, 2),
             c(1, 1, 2, 2, 2),
             c(1, 1, 2, 2, 2),
             c(1, 1, 2, 2, 2),
             c(1, 1, 2, 2, 2),
             c(1, 1, 2, 2, 2),
             c(1, 1, 2, 2, 2),
             c(1, 1, 2, 2, 2),
             c(1, 1, 2, 2, 2),
             c(1, 1, 2, 2, 2))
g <- grid.arrange(arrangeGrob(grobs = gg, ncol = 2, layout_matrix = lay))
setwd('/Users/siold/Proj/FR_RNA-Seq/plots/')
ggsave("horizontal_combined.pdf", g, height = 8.2, width = 11.1, units = "in")

# PCA Plot Images for publications ----------------------------------------

data <- t(log10(counts.cut+1))
data <- data[,apply(data, 2, var, na.rm = TRUE) != 0]
# Remove top1%?
# data <- data[,round(ncol(data)/100):ncol(data)]
# Remove Lowly Expressed Data ?
# data <- data[,colMeans(data) > log2(5)]
# Select top x genes
topGenes <- ncol(data)
topGenes <- as.numeric(topGenes)
data <- data[,1:topGenes]
# compute PCA
data.pca <- prcomp(data, center = T)

# Generate PCA data
x = data.pca$x[,1]
y = data.pca$x[,2]
z = data.pca$x[,3]
gg_data <- x %>%
  tibble::as_tibble() %>%
  tibble::add_column("y" = y,
                     "z" = z,
                     "Label" = names(y),
                     "Tissue" = substr(names(y), 5,5),
                     "Cell" = substr(names(y), 7, 1000000L)) %>%
  dplyr::rename("x" = "value")

pca_labs <- c(
  bquote(paste("CD103"^"-")),
  bquote(paste("CD103"^"+")),
  bquote(paste("CD11b"^"hi")),
  bquote(paste("CD24"^"-")),
  bquote(paste("CD24"^"+")),
  bquote(paste("CD11b"^"lo"))
)
gg[[3]] <- ggplot(gg_data, aes(x = x, y = y, shape = Cell, fill = Cell)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(size = 4) +
  scale_fill_manual(name = "Cell Subsets",
                    labels = pca_labs,
                    values = c("#00BFC4","#F8766D","#619CFF","#F564E3","#00BA38","#B79F00")) +
  scale_shape_manual(name = "Cell Subsets",
                     labels = pca_labs,
                     values = c(23, 23, 22, 21, 21, 22)) +
  labs(title = paste("Principal Component Analysis - Top",topGenes, "most variable genes"),
       y = paste("PC2 -", round(as.numeric(summary(data.pca)$importance["Proportion of Variance","PC2"])*100,2), "%"),
       x = paste("PC1 -", round(as.numeric(summary(data.pca)$importance["Proportion of Variance","PC1"])*100,2), "%")) +
  theme_minimal()
gg[[3]]
gg[[4]] <- ggplot(gg_data, aes(x = x, y = z, shape = Cell, fill = Cell)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(size = 4) +
  scale_fill_manual(name = "Cell Subsets",
                    labels = pca_labs,
                    values = c("#00BFC4","#F8766D","#619CFF","#F564E3","#00BA38","#B79F00")) +
  scale_shape_manual(name = "Cell Subsets",
                     labels = pca_labs,
                     values = c(23, 23, 22, 21, 21, 22)) +
  labs(title = paste("Principal Component Analysis - Top",topGenes, "most variable genes"),
       y = paste("PC3 -", round(as.numeric(summary(data.pca)$importance["Proportion of Variance","PC3"])*100,2), "%"),
       x = paste("PC1 -", round(as.numeric(summary(data.pca)$importance["Proportion of Variance","PC1"])*100,2), "%")) +
  theme_minimal()

ggsave(
  filename = "B6_PCA-1v2_4x5v2.pdf",
  path = "/Users/siold/Proj/FR_RNA-Seq/plots/pca/",
  plot = gg[[3]],
  device = "pdf",
  width = 5,
  height = 4,
  units = "in",
  dpi = 300
)

ggsave(
  filename = "B6_PCA-1v3_4x5v2.pdf",
  path = "/Users/siold/Proj/FR_RNA-Seq/plots/pca/",
  plot = gg[[4]],
  device = "pdf",
  width = 5,
  height = 4,
  units = "in",
  dpi = 300
)
ggsave(
  filename = "B6_PCA-1v2_4x5v2.svg",
  path = "/Users/siold/Proj/FR_RNA-Seq/plots/pca/",
  plot = gg[[3]],
  device = "svg",
  width = 5,
  height = 4,
  units = "in",
  dpi = 300
)

ggsave(
  filename = "B6_PCA-1v3_4x5v2.svg",
  path = "/Users/siold/Proj/FR_RNA-Seq/plots/pca/",
  plot = gg[[4]],
  device = "svg",
  width = 5,
  height = 4,
  units = "in",
  dpi = 300
)
ggsave(
  filename = "B6_PCA-1v2_4x5v2.tiff",
  path = "/Users/siold/Proj/FR_RNA-Seq/plots/pca/",
  plot = gg[[3]],
  device = "tiff",
  width = 5,
  height = 4,
  units = "in",
  dpi = 300
)

ggsave(
  filename = "B6_PCA-1v3_4x5v2.tiff",
  path = "/Users/siold/Proj/FR_RNA-Seq/plots/pca/",
  plot = gg[[4]],
  device = "tiff",
  width = 5,
  height = 4,
  units = "in",
  dpi = 300
)



lay2 <- rbind(c(3, 3, 3, 4, 4, 4),
              c(3, 3, 3, 4, 4, 4),
              c(3, 3, 3, 4, 4, 4),
              c(3, 3, 3, 4, 4, 4),
              c(3, 3, 3, 4, 4, 4),
              c(3, 3, 3, 4, 4, 4),
              c(3, 3, 3, 4, 4, 4),
              c(3, 3, 3, 4, 4, 4),
              c(3, 3, 3, 4, 4, 4),
              c(NA,NA,2, 2, 2, 2),
              c(1, 1, 2, 2, 2, 2),
              c(1, 1, 2, 2, 2, 2),
              c(1, 1, 2, 2, 2, 2),
              c(1, 1, 2, 2, 2, 2),
              c(1, 1, 2, 2, 2, 2),
              c(1, 1, 2, 2, 2, 2),
              c(1, 1, 2, 2, 2, 2),
              c(1, 1, 2, 2, 2, 2),
              c(1, 1, 2, 2, 2, 2),
              c(1, 1, 2, 2, 2, 2),
              c(1, 1, 2, 2, 2, 2),
              c(1, 1, 2, 2, 2, 2),
              c(1, 1, 2, 2, 2, 2),
              c(1, 1, 2, 2, 2, 2),
              c(1, 1, 2, 2, 2, 2))
g2 <- grid.arrange(arrangeGrob(grobs = gg, ncol = 2, layout_matrix = lay2))
setwd('/Users/siold/Proj/FR_RNA-Seq/plots/combined')
ggsave("combined.pdf", g2, height = 13, width = 11.1, units = "in")
ggsave("combined.svg", g2, height = 13, width = 11.1, units = "in", device = "svg")
ggsave("combined.eps", g2, height = 13, width = 11.1, units = "in", device = "eps")
ggsave("combined.tiff", g2, height = 13, width = 11.1, units = "in", device = "tiff")
ggsave("combined.tex", g2, height = 13, width = 11.1, units = "in", device = "tex")
ggsave("combined.ps", g2, height = 13, width = 11.1, units = "in", device = "ps")



# requested to spin the dendrgram and heatmap -----------------------------

setwd('/Users/siold/Proj/FR_RNA-Seq/plots/dend+hm/')

dd <- list()
dd[[1]] <- ggplot(segment(dend_data)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  # geom_text(data = dend_data$labels, aes(x, y, label = final_label),
  #           check_overlap = TRUE, hjust = 0, angle = 0, size = 3) +
  annotate("text", x = dend_data$labels$x, y = dend_data$labels$y,
           label = dend_data$labels$final_label, parse = TRUE, angle = 90, hjust = 1) +
  geom_point(data = classify, aes(x = x, y = y , fill = labs, shape = cell),
             color = "black", size = 4.3) +
  # scale_fill_manual(values = brewer.pal(6, "Set2")) +
  scale_fill_manual(name = "Cell Subsets",
                    # labels = pca_labs,
                    values = c("#00BFC4","#F8766D","#619CFF","#F564E3","#00BA38","#B79F00")) +
  scale_shape_manual(values = c(23, 22, 21)) +
  # ylim(25, 100) +
  theme_classic() +
  theme(axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.line = element_blank(),
        legend.position = "none",
        # plotmargin(top, right, bot, left )
        plot.margin = margin(0, 2.2, 0, 0.5, "in")) +
  # scale_y_reverse(expand = c(0.2, 0)) +
  ylab("") +
  xlab("")
dd[[1]]
annotation_row$DEGs <- NULL#rep(" ", 497)
xx <- pheatmap::pheatmap(final_dat[,labels],
                         color = viridisLite::turbo(512),
                         cluster_cols = FALSE, cluster_rows = FALSE,
                         show_colnames = FALSE, show_rownames = FALSE,
                         annotation_row = annotation_row,
                         cellheight = 1, cellwidth = 29.23529)
# gaps_row = c(59,272,360, 396), gaps_col = c(6,11),
# main = "Scaled VSTpk expression, 497 genes, B6 mice")
dd[[2]] <- xx[[4]]
lay3 = rbind(1,1,2,2,2)
g3 <- grid.arrange(arrangeGrob(grobs = dd, ncol = 1, layout_matrix = lay3))

ggsave("combined2.pdf", g3, height = 12, width = 10, units = "in")
ggsave("combined2.svg", g3, height = 12, width = 10, units = "in", device = "svg")
ggsave("combined2.eps", g3, height = 12, width = 10, units = "in", device = "eps")
ggsave("combined2.tiff", g3, height = 12, width = 10, units = "in", device = "tiff")
# ggsave("combined.tex", g3, height = 30.48, width = 25.4, units = "cm", device = "tex") #metric not avail?
ggsave("combined2.ps", g3, height = 12, width = 10, units = "in", device = "ps")

xx <- pheatmap::pheatmap(final_dat[,labels],
                         color = viridisLite::viridis(512),
                         cluster_cols = FALSE, cluster_rows = FALSE,
                         show_colnames = FALSE, show_rownames = FALSE,
                         annotation_row = annotation_row,
                         cellheight = 0.5, cellwidth = 14.61764)
dd[[1]] <- dd[[1]] + theme(plot.margin = margin(0, 1.1, 0.3, 0.25, "in"))
dd[[2]] <- xx[[4]]
g3 <- grid.arrange(arrangeGrob(grobs = dd, ncol = 1, layout_matrix = lay3))


ggsave("combined_small.pdf", g3, height = 6, width = 6, units = "in")
ggsave("combined_small.svg", g3, height = 6, width = 6, units = "in", device = "svg")
ggsave("combined_small.eps", g3, height = 6, width = 6, units = "in", device = "eps")
ggsave("combined_small.tiff", g3, height = 6, width = 6, units = "in", device = "tiff")
# ggsave("combined.tex", g3, height = 30.48, width = 25.4, units = "cm", device = "tex") #metric not avail?
ggsave("combined_small.ps", g3, height = 6, width = 6, units = "in", device = "ps")

