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


# Findings
