# Marble crayfish parthenogenesis
# Date: April 25th

## What is the plan
# summry
what genes in marble crayfish are homologous with other species that aren't and are parthenogenisis using known genes in Meiosis to find homologues genes between parthogenesis speices and non-parthogensis species. 
# genes 
  1. K10318 (EMI2,FBXO43)
  2. K03363 (CDC20)
  3. K05867 (CDC25C)
  4. K02087 (CDC2,CDK1)
  5. Cyclin B
  6. Mos
  7. PLCZ
  8. IP3R
  9. Calm
# organisms
 - Marbled Crayfish, Procambarus Virginalis
 - Slough Crayfish, Procambarus fallax, no morphilogical difference from the marbled crayfish
 - Redswamp crayfish, Procambarus clarkii
# Data to down load genes
1. Genes
   reference.fasta
2. organism genomes
   -
## DOWNLOAD DATA
```
module load anaconda3
conda create -n sra-tools -c bioconda -c conda-forge sra-tools
```


```
vi download_sra.sh
```
- Type I
```
#!/bin/bash
module load anaconda3
conda activate sra-tools
 while read -r SRR; do
   echo "Downloading $SRR..."
   prefetch --max-size 100G $SRR 
   fasterq-dump --gzip $SRR --split-files -O fastq_files/
done < SRR_Accessions.txt
conda deactivate
```
-script may give you an error if anaconda3 was loaded before script was run, exit the HPC and run it again.
-Run script using
```
bash download_sra.sh
```
# Quantify Gene Expression in Salmon
```
module load anaconda3
conda create -n salmon_env salmon -c bioconda -c conda-forge
conda activate salmon_env
```
- Building salmon index
```
salmon index -t reference.fasta -i salmon_index
```
- Write a slurm script to run salmon using all the sequences at the same time
```
vi salmon.slurm
```
-Type `I` and paste the following:
```
#!/bin/bash
#SBATCH --job-name=salmon_quant
#SBATCH --output=logs/salmon_%A_%a.out
#SBATCH --error=logs/salmon_%A_%a.err
#SBATCH --array=1-34
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

module load anaconda3
conda activate salmon_env

# Get accession number for this task
ACCESSION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" SRR_Accessions.txt)

# Run Salmon
salmon quant -i salmon_index \
             -l A \
             -1 fastq_files/${ACCESSION}_1.fastq \
             -2 fastq_files/${ACCESSION}_2.fastq \
             -p 8 \
             -o salmon_output/${ACCESSION}_quant
```
- Escape and type `:wq`
- Run the slurm scrip
```
sbatch salmon.slurm

### We are relaxing the salmon mapping, since results suggest mapping parameters were too stringent.
- The following script was run


# Findings
The old quant.sf files were used to visualize gene expression in R. The following R script was used:
- Note that folder `salmon_output` and `sample_metadata.csv` must be in the same R directory. The R script was also added to this repository
```
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
  ```


