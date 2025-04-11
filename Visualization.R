#Conditional Package Installation
packages <- c("dplyr", "ggplot2", "tximport", "DESeq2")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}
# For Bioconductor packages like DESeq2 and tximport:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_packages <- c("tximport", "DESeq2")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

#Load Packages
library(dplyr)
library(ggplot2)
library(tximport)
library(DESeq2)


#Import Quant.sf Data
metadata <- read.csv("sample_metadata.csv", stringsAsFactors = FALSE)
metadata$folder <- paste0("salmon_output/", metadata$accession, "_quant")
files <- file.path(metadata$folder, "quant.sf")
names(files) <- metadata$accession


txi <- tximport(files, type = "salmon", txOut = TRUE)
rownames(metadata) <- metadata$accession
names(files)

#Differential Gene Expression
metadata$species <- relevel(factor(metadata$species), ref = "virginalis")
dds <- DESeqDataSetFromTximport(txi, colData = metadata, design = ~ species)
dds <- DESeq(dds)
resultsNames(dds) 

res_vir_fall <- results(dds, contrast = c("species", "virginalis", "fallax"))
res_vir_clar <- results(dds, contrast = c("species", "virginalis", "clarkii"))


# Combine results
df1 <- as.data.frame(res_vir_fall) %>%
  mutate(gene = rownames(.), contrast = "Virginalis vs Fallax")

df2 <- as.data.frame(res_vir_clar) %>%
  mutate(gene = rownames(.), contrast = "Virginalis vs Clarkii")

combined_df <- bind_rows(df1, df2) %>%
  filter(!is.na(padj)) %>%
  mutate(
    signif = padj < 0.05,
    label = ifelse(signif, "*", "")
  )


# Change Na to 0
combined_df$gene <- factor(combined_df$gene, levels = unique(combined_df$gene))
combined_df$log2FoldChange[is.na(combined_df$log2FoldChange)] <- 0
combined_df$padj[is.na(combined_df$padj)] <- 1  # mark as non-significant
combined_df$label[is.na(combined_df$label)] <- ""  # no asterisk



# Plot
ggplot(combined_df, aes(x = log2FoldChange, y = gene, fill = signif)) +
  geom_col(width = 0.6, color = "black") +  # Adds outline for 0-height bars
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_text(aes(label = label),
            hjust = ifelse(combined_df$log2FoldChange > 0, -0.3, 1.3),
            size = 6,
            na.rm = TRUE) +
  facet_wrap(~ contrast) +
  scale_fill_manual(values = c("grey80", "firebrick")) +
  scale_x_continuous(
    limits = c(-10, 2),
    breaks = seq(-10, 2, by = 2)
  ) +
  labs(
    x = "Log2 Fold Change", y = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 14),
    strip.background = element_rect(fill = "grey90", color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.major.y = element_blank()
  )
