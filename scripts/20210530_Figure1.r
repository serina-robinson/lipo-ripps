# Install packages
pacman::p_load("tidyverse", "Biostrings", "randomcoloR","umap",
               "data.table", "vegan", "RColorBrewer", "ggrepel",
               "scales", "ggpubr", "patchwork", "caret", "viridis")

# Set random seed
set.seed(1234)

# Read in the dataset
# dat <- read_csv("data/combined_Nif11_NHLP_cleavage_motifs/1921_full_Nif11_NHLP_genome_nbs_with_clans_precursors.csv") %>%
dat <- read_csv("data/combined_Nif11_NHLP_cleavage_motifs/1922_full_Nif11_NHLP_genome_nbs_with_clans_precursors.csv") %>%  
  janitor::clean_names() %>%
  dplyr::mutate(grouping_id = paste0(query, "_", genus_species)) %>%
  group_by(grouping_id) %>%
  dplyr::add_count(name = "pfam_count", pfam_long) %>%
  dplyr::add_count(name = "clan_count", clan_keep) %>%
  ungroup() %>%
  dplyr::mutate(precursor = case_when(
    grepl("PF07862", pfam_id1) ~ "Nif11",
    grepl("PF07862", pfam_id2) ~ "Nif11",
    TRUE ~ "NHLP")) %>%
  dplyr::mutate(nif11_queries = case_when(
    precursor == "Nif11" ~ query,
    TRUE ~ NA_character_)) %>%
  dplyr::mutate(precursor_type = case_when(
        query %in% nif11_queries ~ "Nif11",
        TRUE ~ "NHLP")) %>%
  dplyr::mutate(known_product = case_when(
    grepl("WP_007357382.1", query) ~ "landornamide",
    # grepl("WP_019503884.1", query) ~ "spliceotides", # removed bc product not published
    grepl("aerodenitrificans", genus_species) ~ "aeronamide",
    grepl("polytheonamide", genus_species) ~ "polytheonamide",
    grepl("WP_036784545", query) ~ "pha",
    grepl("WP_007358513.1", query) ~ "ksp",
    grepl("WP_012412980.1", query) ~ "npu",
    TRUE ~ "unknown"))


nost <- dat[grep("WP_012412980.1", dat$query),]

# Find nbs for characterized NPs
char <- read_csv("data/characterized_nps/proteusin_key_nbs_pull/main_co_occur.csv")

# Optional: checking characterized products
# proc <- dat %>%
#   dplyr::filter(grepl("Prochlorococcus marinus", genus_species))
# 
# kampto <- dat %>%
#   dplyr::filter(grepl("Kamptonema", genus_species))
# 
# nosto <- dat %>%
#   dplyr::filter(grepl("WP_012412980.1", query))
# 
# # pleuro <- dat %>%
# #   dplyr::filter(grepl("WP_019503884.1", query))
# 
# microvirg <- dat %>%
#   dplyr::filter(grepl("aerodenitrificans", genus_species))
# 
# ento <- dat %>%
#   dplyr::filter(grepl("polytheonamide", genus_species))


# Reshape data frome long to wide
dat_trim <- dat
mat <- reshape2::dcast(dat_trim, grouping_id ~ pfam_long, value.var = "pfam_count") 
mat_bool <- mat
mat_bool[mat > 1] <- 1

# Barplot of most common co-occurring PFAMs
bpldf <- data.frame(sort(colSums(mat_bool[,2:ncol(mat_bool)]), decreasing = T)[1:40]) %>%
  rownames_to_column()
colnames(bpldf) <- c("PFAM", "Count")

bpldf <- bpldf %>% 
  dplyr::filter(PFAM != "NA") %>%
  dplyr::filter(PFAM != "NA_NA") %>%
  dplyr::mutate(coloring = ifelse(PFAM %in% "PF00583_Acetyltransf_1", "forestgreen", "gray60")) %>%
  dplyr::mutate("% BGCs containing PFAM" = round((Count/2347) * 100, 2))

bpldf$PFAM <- factor(bpldf$PFAM, levels = rev(as.character(bpldf$PFAM)))

pdf("output/Supplementary_fig_S1.pdf", width = 8.5, height = 12)
bpl_s1 <- ggplot(data = bpldf, aes(x = PFAM, y = `% BGCs containing PFAM`)) + 
  geom_bar(stat = "identity", fill = bpldf$coloring) + 
  theme_pubr() +
 # theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
  #scale_x_discrete(expand = expansion(mult = c(0, 0))) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  coord_flip()
bpl_s1
dev.off()

bpl_for_tab <- bpldf %>%
  dplyr::select(-coloring)
write_csv(bpl_for_tab, "output/Supplementary_table_S1.csv")

bpl_maturases <- bpldf %>%
  dplyr::filter(PFAM %in% c("PF04055_Radical_SAM",
                            "PF13575_DUF4135",
                            "PF00583_Acetyltransf_1",
                            "PF02624_YcaO",
                            "PF00881_Nitroreductase",
                            "PF13561_adh_short_C2",
                            "PF13649_Methyltransf_25",
                            "PF00106_adh_short",
                            "PF01593_Amino_oxidase")) 

bpl_maturases$PFAM <- gsub("_", " ", bpl_maturases$PFAM)
bpl_maturases$PFAM[grep(" ", bpl_maturases$PFAM)]
bpl_maturases$`Protein family` <- c("Radical SAM (PF04055)",
               "DUF4135 (PF13575)",
               "Acyltransferase (PF00583)",
               "YcaO (PF02624)",
               "Nitroreductase (PF00881)",
               "Enoyl-ACP reductase (PF13561)",
               "Methyltransferase (PF13649)",
               "Short-chain dehydrogenase (PF00106)",
               "Amine oxidoreductase (PF01593)")
   
bpl_maturases$Family<- factor(bpl_maturases$`Protein family`, 
                              levels = rev(c("Radical SAM (PF04055)",
                                        "DUF4135 (PF13575)",
                                        "Acyltransferase (PF00583)",
                                        "YcaO (PF02624)",
                                        "Nitroreductase (PF00881)",
                                        "Enoyl-ACP reductase (PF13561)",
                                        "Methyltransferase (PF13649)",
                                        "Short-chain dehydrogenase (PF00106)",
                                        "Amine oxidoreductase (PF01593)")))

bpl <- ggplot(data = bpl_maturases, aes(x = Family, y = `% BGCs containing PFAM`)) + 
  geom_bar(stat = "identity", fill = bpl_maturases$coloring) + 
  theme_pubr() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)))
bpl

# Matrix of PFAM clans
#  clanmat <- reshape2::dcast(dat_trim, grouping_id ~ clan_keep, value.var = "clan_count")
# sort(colSums(clanmat[,2:ncol(clanmat)]), decreasing = T)[1:15]

keep_tax_info <- dat_trim %>%
  dplyr::select(grouping_id, known_product, precursor_type, genvar, genus_species, 
                nucleotide_acc, gtdb_tax, phylum, order, class, 
                family, genus_1) %>%
  distinct()

mat_df <- mat %>%
  right_join(keep_tax_info, ., by = "grouping_id") %>%
 #  left_join(., clanmat, by = "grouping_id") %>%
  distinct(grouping_id, .keep_all = T) %>%
  mutate(kingdom = "Bacteria")

mat_num <- data.frame(mat_df[, grep("^PF", colnames(mat_df))]) 
rownames(mat_num) <- mat_df$grouping_id

# Optional: Compare PCoA to UMAP
# mat_num_umap <- umap(mat_num)
# umap_plot_df <- data.frame(mat_num_umap$layout) %>%
#   mutate(grouping_id = rownames(mat_num_umap$data)) %>%
#   inner_join(., mat_df, by = "grouping_id") 
# 
# head(umap_plot_df)[1:5, 1:5]
# 
# umap_plot <- ggplot(
#   umap_plot_df,
#   aes(
#     x = X1,
#     y = X2
#   )
# ) +
#   geom_point() +
#   theme_classic() 
# umap_plot

# Remove columns without variance across clusters
nozdat <- nearZeroVar(mat_num, 
                      freqCut = 98.9/1.1,
                      saveMetrics = TRUE)
which_rem <- rownames(nozdat)[nozdat[,"nzv"] == TRUE]

findat <- mat_num %>%
  dplyr::select(-all_of(which_rem))

#  "PF02979", "PF14407", # NHLP precursors
#  "PF07862", # Nif11 precursor

# Normalize
findat[1:10, 1:10]
mat_reg <- sweep(findat, 1, rowSums(findat), '/') 
mat_fin <- mat_reg[!is.na(rowSums(mat_reg)),]

# Run Principal Coordinates Analysis using vegan package
d.jaccard <- vegdist(mat_fin, method = "jaccard") # jaccard distance metric
pc.bray <- cmdscale(d.jaccard, eig = TRUE, k = 2)
pc1 <- pc.bray$GOF[1]
pc2 <- (pc.bray$GOF[2]-pc.bray$GOF[1])

xlab <- paste0("PC 1 (",round(pc1*100,2),"% explained var.)")
ylab <- paste0("PC 2 (",round(pc2*100,2),"% explained var.)")

plotdf <- data.frame(grouping_id = rownames(pc.bray$points),
                     PC1 = pc.bray$points[,1],
                     PC2 = pc.bray$points[,2]) %>%
  left_join(mat_df, by = "grouping_id")

# Set custom color palette for plotting
pal <- colorRampPalette(brewer.pal(8, "Set1"))(8)
pal <- colorRampPalette(brewer.pal(8, "Set1"))(8)
pal

pal2 <- c(pal[c(1:2, 4:5, 8)], "seagreen2", "gray80",  "black")
pal2

# PCoA of taxonomy assignments (phylum level)
prod <- ggplot() +
  geom_point(data = plotdf,
             alpha = ifelse(plotdf$known_product == "unknown", 0.4, 1),
             size = ifelse(plotdf$known_product == "unknown", 1, 3),
             aes(
               x = PC1,
               y = PC2,
               color = known_product)) +
  theme_pubr() +
  theme(legend.text=element_text(size=14),
    legend.title = element_text(size = 16),
    legend.position = "none") +
  scale_color_manual(values = pal2) +
  geom_text_repel(data = subset(plotdf, known_product != "unknown"),
    size = 5,
    aes(x = PC1,
        y = PC2,
        label = known_product)) +
  guides(colour = guide_legend(override.aes = list(size=5),
                               ncol = 2))
prod


# Calculate weighted averages of BGCs
# Normalize by dividing by column sums
pfams.norm <- sweep(mat_fin, 2, colSums(mat_fin), '/')

# Use matrix multiplication to calculated weighted average 
# of each PFAM Clan along each axis
wa <- t(pfams.norm) %*% pc.bray$points[,1:2] # weighted average
wa_df <- data.frame(CL = rownames(wa),
                    WA_1 = wa[,1],
                    WA_2 = wa[,2])

# GNAT-family acyltransferases PF00583
acetyl <- ggplot() +
  geom_point( data = plotdf,
              alpha = 0.8,
              aes(
                x = PC1,
                y = PC2,
                color = PF00583_Acetyltransf_1)) +
  theme_pubr() +
  scale_color_viridis(name = "GNAT-acyltransferases\n(PF00583)\nper BGC",
                      direction = -1, 
                      option = "magma",
                      end = 0.9,
                      na.value = "gray80") 
  #theme(axis.title.x = element_blank())
acetyl

# C39 Peptidases
c39 <- ggplot() +
  geom_point(data = plotdf,
             alpha = 0.8,
             aes(
               x = PC1,
               y = PC2,
               color = PF03412_Peptidase_C39)) +
  theme_pubr() +
  scale_color_viridis(name = "Peptidase C39\n(PF03412)\nper BGC",
                      direction = -1, 
                      option = "magma",
                      end = 0.9,
                      na.value = "gray80")  #+
c39

abc <- ggplot() +
  geom_point(data = plotdf,
             alpha = 0.8,
             aes(
               x = PC1,
               y = PC2,
               color = PF00005_ABC_tran)) +
  theme_pubr() +
  scale_color_viridis(name = "ABC trans.\n(PF00005)\n per BGC",
                      direction = -1, 
                      option = "magma",
                      end = 0.9,
                      na.value = "gray80")  #+
abc

# Radical SAMs
sam <- ggplot() +
  geom_point(data = plotdf,
             alpha = 0.6,
             aes(
               x = PC1,
               y = PC2,
               color = PF04055_Radical_SAM)) +
  theme_pubr() +
  scale_color_viridis(name ="Radical SAMs\n (PF04055)\n per BGC",
                      direction = -1,
                      option = "magma",
                      end = 0.9,
                      na.value = "gray80")
sam


# NHLP vs. Nif11 precursors
colnames(plotdf)[grep("PF07862", colnames(plotdf))]
table(plotdf$PF07862_Nif11 > 0)

nif11 <- ggplot() +
  geom_point(data = plotdf,
             alpha = 0.6,
             aes(
               x = PC1,
               y = PC2,
               color = precursor_type)) +
  theme_pubr() +
  #theme(axis.title.x = element_blank()) +
  scale_color_manual(values = magma(8)[c(5, 7)],
                     name = "RiPP\nprecursor\ntype")
nif11


# DUF4135
duf4135 <- ggplot() +
  geom_point(data = plotdf,
             alpha = 0.6,
             aes(
               x = PC1,
               y = PC2,
               color = PF13575_DUF4135)) +
  theme_classic() +
  #theme(axis.title.x = element_blank()) +
  scale_color_viridis(name ="DUF4135 \n  (PF13575) \n per BGC",
                      direction = -1,
                      option = "magma",
                      end = 0.9,
                      na.value = "gray80")  #+
duf4135

# Set additional colors in custom palette
pal_distinct <- palette(colorRampPalette(brewer.pal(8,"Set2"))(8))
colorblind <- c(pal[1], "#CC6677",  "#999933",  "steelblue4", "#AA4499",  
                 pal_distinct[5:7], 
                "#882255", "#6699CC",  pal_distinct[2:3], "gray40",
                "navyblue", "mistyrose1", "#88CCEE", pal_distinct[4], #"#661100", 
                "orange", "goldenrod", "purple4", "lightseagreen")
# "#117733"
show_col(colorblind)
plotdf$phylum[is.na(plotdf$phylum)] <- "Other"

cr_df <- data.frame(sort(table(plotdf$phylum), decreasing = T), stringsAsFactors = F) 
colnames(cr_df) <- c("Phylum", "Freq")
sum(cr_ord$Freq)
cr_df$Phylum <- as.character(cr_df$Phylum)
cr_ord <- cr_df[order(cr_df$Phylum),] # order alphabetically
cr_ord$lab <- round((cr_ord$Freq/sum(cr_ord$Freq) * 100), 1)
cr_ord$lab[cr_ord$lab < 3.0] <- ""
cr_ord$lab <- ifelse(cr_ord$lab == "", "", paste0(cr_ord$lab, "%"))

donut <- ggdonutchart(data = cr_ord,
                      x = "Freq",
                      label = "lab",
                      fill = "Phylum",
                      lab.font = c(6, "plain", "black"),
                      color = "white",
                      palette = colorblind)
donut

phygraph <- ggplot(data = plotdf) +
  geom_point(data = plotdf,
             alpha = 0.8,
             size = 2,
             aes(
               x = PC1,
               y = PC2,
               color = phylum)) +
  scale_color_manual(values = colorblind,
                     name = "Phylum") +
  theme_pubr()  +
  guides(colour = guide_legend(override.aes = list(size= 3),
                               ncol = 1))
  #guides(colour = guide_legend(override.aes = list(size=2),))
phygraph


phygraph <- phygraph +
  ggtitle('A') +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        legend.text=element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.title = element_text(size = 24),
        legend.position = "none")
donut <- donut +
  ggtitle('Z') +
  theme(axis.text = element_text(size = 16),
        axis.title = element_blank(),
        plot.title = element_blank(),
        legend.position = "none")
prod <- prod +
  ggtitle('B') +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 24))
bpl <- bpl +
  ggtitle('C') +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 24),
        axis.title.y = element_blank()) +
  coord_flip()
bpl
# sam <- sam + 
#   ggtitle('C') +
#   theme(axis.text = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         plot.title = element_text(size = 20),
#         legend.position = "right")
c39 <- c39 + 
  ggtitle('F') +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 24),
        legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        legend.position = "right")
abc <- abc +
  ggtitle('E') +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 24),
        legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        legend.position = "right")
nif11 <- nif11 +
  ggtitle('D') +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        plot.title = element_text(size =  24),
        legend.text=element_text(size=16),
        legend.title = element_text(size = 14),
        legend.position = "right")
acetyl <- acetyl +
  ggtitle('G') +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 24),
        legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        legend.position = "right")

#ZZZZZZ
# layout <- '
# AAAAAAAABBB
# CCCCCCCEDFG
# CCCCCCCHHHH'
#CCDDHH
#EEGGHH

layout <- '
AAAAAABBBCCCCCCCC
DDDDDEEEHHHHHHHHH
DDDDDEEEHHHHHHHHH
FFFFFGGGHHHHHHHHH
FFFFFGGGHHHHHHHHH'

# pdf("output/Figure1.pdf", width = 21, height = 14)
# wrap_plots(
#   A = prod,
#   B = phygraph,
#   C = donut,
#   D = nif11,
#   #D = abc,
#   E = abc,
#   F = c39,
#   G = acetyl,
#   H = bpl, 
#   heights = c(8,1.5,1.5,1.5,1.5),
#  # widths = c(3,1,1,
#   #           1,1,3),
#   # widths = c(1, 1, 2, 1),
#   # 
#   design = layout)
# dev.off()


layout1 <- '
AAAABBBB
DDDDCCC#'

# New figure (as per Flo's specifications)
pdf("output/Figure1_part1.pdf", width = 15, height = 12)
wrap_plots(
  A = phygraph,
  B = donut,
  C = bpl,
  D = prod,
  heights = c(3,3),
  widths = c(2,1),
  design = layout1)
dev.off()

# pdf("output/Figure1_part1.pdf", width = 14, height = 10)
# phygraph + donut + prod + bpl + plot_layout(widths = 1)
# dev.off()

layout2 <-'
EEEFFF
GGGHHH'

pdf("output/Figure1_part2.pdf", width = 15, height = 10)
wrap_plots(  
  E = nif11,
  F = abc,
  G = c39,
  H = acetyl,
#  heights = c(3,4,2,2),
  widths = c(1,1,1,1),
 # widths = c(3,1,1,
  #           1,1,3),
  # widths = c(1, 1, 2, 1),
  #
  design = layout2)
dev.off()


c39 <- c39 + 
  ggtitle('F') +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 24),
        legend.text=element_text(size=12),
        legend.title = element_text(size = 14),
        legend.position = "top")
abc <- abc +
  ggtitle('E') +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 24),
        legend.text=element_text(size=12),
        legend.title = element_text(size = 14),
        legend.position = "top")
nif11 <- nif11 +
  ggtitle('D') +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        plot.title = element_text(size =  24),
        legend.text=element_text(size=14),
        legend.title = element_text(size = 14),
        legend.position = "top")
acetyl <- acetyl +
  ggtitle('G') +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        plot.title = element_text(size = 24),
        legend.text=element_text(size=12),
        legend.title = element_text(size = 14),
        legend.position = "top")

layout3 <-'
EEEFFFGGGHHH'
pdf("output/Figure1_part2_horiz.pdf", width = 15, height = 5)
wrap_plots(  
  E = nif11,
  F = abc,
  G = c39,
  H = acetyl,
  #  heights = c(3,4,2,2),
  widths = c(1,1,1,1),
  # widths = c(3,1,1,
  #           1,1,3),
  # widths = c(1, 1, 2, 1),
  #
  design = layout3)
dev.off()


