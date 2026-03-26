# Homework 4: pipelines, plotting, genome assembly
```
wget https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/dmel_r6.66_FB2025_05/fasta/dmel-all-chromosome-r6.66.fasta.gz
```

# Part ONE: Calculate the following for all sequencess up to 100kb and over 100kb
## 1. Total number of nucleotides <= 100kb is 6178042 and > 100kb is 137547960
```
< dmel-all-chromosome-r6.66.fasta.gz \
bioawk -c fastx \
'length($seq) <= 100000 { count += length($seq) } END { print count+0 }'

< dmel-all-chromosome-r6.66.fasta.gz \
bioawk -c fastx \
'length($seq) > 100000 { count += length($seq) } END { print count+0 }'
```
## 2. Total number of N's <= 100kb is 662593 and > 100kb is 490385
```
< dmel-all-chromosome-r6.66.fasta.gz \
bioawk -c fastx \
'length($seq) <= 100000 { count += gsub(/[Nn]/, "", $seq) } END { print count+0 }' 
< dmel-all-chromosome-r6.66.fasta.gz bioawk -c fastx 'length($seq) > 100000 { count += gsub(/[Nn]/, "", $seq) } END { print count+0 }'
```
## 3. Total number of sequences <= 100kb is 1863 and > 100kb is 7
```
< dmel-all-chromosome-r6.66.fasta.gz \
bioawk -c fastx \
'{ if (length($seq) <= 100000) count++ } END { print count+0 }'
< dmel-all-chromosome-r6.66.fasta.gz \
bioawk -c fastx \
'{ if (length($seq) > 100000) count++ } END { print count+0 }'
```

# Part TWO: Plots of the following for all sequences <= 100kb and > 100kb
## 1. Sequence length distribution
```
# Create separate files for seq <= 100kb and >100kb
< dmel-all-chromosome-r6.66.fasta.gz \
bioawk -c fastx \
'{ if (length($seq) <= 100000)print $name "\t" length($seq)}' >  dmel_100kb_or_less.txt
< dmel-all-chromosome-r6.66.fasta.gz \
bioawk -c fastx \
'{ if (length($seq) > 100000)print $name "\t" length($seq)}' >  dmel_over_100kb.txt
```
```
#Moved files to desktop, continue in R:
#For seq 100kb or less:
library(ggplot2)
#Define files
input_file <- "~/Desktop/dmel_100kb_or_less.txt"
output_file <- "~/Desktop/dmel_100kb_or_less_histogram.png"
#Define histogram
data <- read.table(input_file, sep = "\t", header = FALSE, col.names = c("Name", "Length"))
p <- ggplot(data, aes(x = Length)) + geom_histogram(fill = "mediumseagreen", color = "black", bins = 50) + theme_minimal() + labs(title = "Sequence Lengths (100kb or less)", x = "Sequence Length (bp)", y = "Frequency")
ggsave(output_file, plot = p, width = 8, height = 6, dpi = 300)
```
<img width="2400" height="1800" alt="dmel_100kb_or_less_histogram" src="https://github.com/user-attachments/assets/2fa290bb-3710-408e-8793-58873b87f664" />


```
#For seq over 100kb:
input_file <- "~/Desktop/dmel_over_100kb.txt"
output_file <- "~/Desktop/dmel_over_100kb_histogram.png"
#Define histogram
data <- read.table(input_file, sep = "\t", header = FALSE, col.names = c("Name", "Length"))
p <- ggplot(data, aes(x = Length)) + geom_histogram(fill = "mediumseagreen", color = "black", bins = 50) + theme_minimal() +
labs(title = "Sequence Lengths (over 100kb)", x = "Sequence Length (bp)", y = "Frequency")
ggsave(output_file, plot = p, width = 8, height = 6, dpi = 300)
```
<img width="2400" height="1800" alt="dmel_over_100kb_histogram" src="https://github.com/user-attachments/assets/2e68dc15-2027-488b-a560-b2107d97e467" />

## 2. Sequence GC% distribution
```
# Create separate files for seq <= 100kb and >100kb
< dmel-all-chromosome-r6.66.fasta.gz \
bioawk -c fastx \
'{ if (length($seq) <= 100000) print $name "\t" length($seq) "\t" gc($seq) }' > dmel_gc_100kb_or_less.txt
< dmel-all-chromosome-r6.66.fasta.gz \
bioawk -c fastx \
'{ if (length($seq) > 100000) print $name "\t" length($seq) "\t" gc($seq) }' > dmel_gc_over_100kb.txt

#Moved files to desktop, continue in R:
#For seq 100kb or less:
library(ggplot2)
#Define files

input_file <- "~/Desktop/dmel_gc_100kb_or_less.txt"
output_file <- "~/Desktop/dmel_gc_100kb_or_less_histogram.png"

#Define histogram
data <- read.table(input_file, sep = "\t", header = FALSE, col.names = c("Name", "Length", "GC"))
p <- ggplot(data, aes(x = GC)) + geom_histogram(fill = "mediumseagreen", color = "black", bins = 50) + theme_minimal() +
labs(title = "GC Distribution (100kb or less)", x = "GC content", y = "Frequency")
ggsave(output_file, plot = p, width = 8, height = 6, dpi = 300)
```
<img width="2400" height="1800" alt="dmel_gc_100kb_or_less_histogram" src="https://github.com/user-attachments/assets/41b125bf-d6a3-4d18-b1dd-393882b5015a" />

```
#For seq over 100kb:
library(ggplot2)

#Define files
input_file <- "~/Desktop/dmel_gc_over_100kb.txt"
output_file <- "~/Desktop/dmel_gc_over_100kb_histogram.png"

#Define histogram
data <- read.table(input_file, sep = "\t", header = FALSE, col.names = c("Name", "Length", "GC"))
p <- ggplot(data, aes(x = GC)) + geom_histogram(fill = "mediumseagreen", color = "black", bins = 50) + theme_minimal() +
labs(title = "GC Distribution (over 100kb)", x = "GC content", y = "Frequency")
ggsave(output_file, plot = p, width = 8, height = 6, dpi = 300)
```
<img width="2400" height="1800" alt="dmel_gc_over_100kb_histogram" src="https://github.com/user-attachments/assets/ace4ef7a-455c-42c0-b242-283d80fec0f2" />


## 3. Cumulative sequence size sorted from largest to smallest
```
# For seq 100kb or less
gunzip dmel-all-chromosome-r6.66.fasta.gz \
| bioawk -c fastx \
'{ if (length($seq) <= 100000) print length($seq) }' \
| sort -rn \
| gawk 'BEGIN {print "Length\tAssembly"} {print $1 "\tDmel_100kb_or_less"}' > dmel_100kb_or_less_lengths.txt \
| plotCDF2 dmel_100kb_or_less_lengths.txt cumulative_sizes_100kb_or_less.png
```
<img width="640" height="480" alt="cumulative_sizes_100kb_or_less" src="https://github.com/user-attachments/assets/2198553f-04af-41b5-9df1-1c19bc70902d" />

```
#For seq over 100kb
cat dmel-all-chromosome-r6.66.fasta \
| bioawk -c fastx \
'{ if (length($seq) > 100000) print length($seq) }' \
| sort -rn \
| gawk 'BEGIN {print "Length\tAssembly"} {print $1 "\tDmel_over_100kb"}' > dmel_over_100kb_lengths.txt \
| plotCDF2 dmel_over_100kb_lengths.txt cumulative_sizes_over_100kb.png
```
<img width="640" height="480" alt="cumulative_sizes_over_100kb" src="https://github.com/user-attachments/assets/d91a0ad2-6344-4e8a-a611-051b24203745" />

# Part THREEa: Assemble a genome using Pacbio HiFi reads
```
wget https://hpc.oit.uci.edu/~solarese/ee282/iso1_onp_a2_1kb.fastq.gz

#!/usr/bin/env bash
#SBATCH --job-name=hifiasm
#SBATCH --partition=standard
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=hifiasm_%j.out
#SBATCH --error=hifiasm_%j.err

source /data/homezvol1/mnahmou/miniforge3/etc/profile.d/conda.sh
conda activate ee282

# Run hifiasm
hifiasm \
-o hifi_fly_assembly \
-t 16 \
ISO_HiFi_Shukla2025.fasta.gz
```
# Part THREEb: Assembly assessment
## 1. Calculate the N50 of the above assembly (~21Mbp)
```
awk '/^S/{print ">"$2"\n"$3}' hifi_fly_assembly.bp.p_ctg.gfa > hifi_fly_assembly.fasta
bioawk -c fastx '{print length($seq)}' hifi_fly_assembly.fasta | sort -rn | awk '{a[NR]=$1; s+=$1} END {for(i=1;i<=NR;i++){c+=a[i]; if(c>=s/2){print "N50: " a[i]; exit}}}'
N50: 21715751
```
## 2. Compare assembly to the contig assembly and scaffold
Reference contig N50: 21.5 Mbp

Assembly contig N50: 21.7 Mbp 

(good!!)
## 3. Calculate BUSCO scores of both assemblies and compare them
For Shukla assembly: 
```
#!/usr/bin/env bash
#SBATCH --job-name=busco
#SBATCH --partition=standard
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=busco_%j.out
#SBATCH --error=busco_%j.err

source /data/homezvol1/mnahmou/miniforge3/etc/profile.d/conda.sh
conda activate ee282

# Run busco on fly genome
busco -i hifi_fly_assembly.fasta -o busco_fly_eval -m genome -l diptera_odb10 -c 16
```
***** Results: *****

	C:99.8%[S:99.6%,D:0.2%],F:0.0%,M:0.2%,n:3285,E:10.7%	   
	3280	Complete BUSCOs (C)	(of which 351 contain internal stop codons)		   
	3273	Complete and single-copy BUSCOs (S)	   
	7	Complete and duplicated BUSCOs (D)	   
	0	Fragmented BUSCOs (F)			   
	5	Missing BUSCOs (M)			   
	3285	Total BUSCO groups searched	

For GenBank reference:
```
#!/usr/bin/env bash
#SBATCH --job-name=busco
#SBATCH --partition=standard
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=busco_%j.out
#SBATCH --error=busco_%j.err

source /data/homezvol1/mnahmou/miniforge3/etc/profile.d/conda.sh
conda activate ee282

# Run busco on fly genome
busco -i GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz -o busco_ref_fly_eval -m genome -l diptera_odb10 -c 16
```
***** Results: *****

	C:99.9%[S:99.7%,D:0.3%],F:0.0%,M:0.1%,n:3285,E:10.7%	   
	3283	Complete BUSCOs (C)	(of which 351 contain internal stop codons)		   
	3274	Complete and single-copy BUSCOs (S)	   
	9	Complete and duplicated BUSCOs (D)	   
	0	Fragmented BUSCOs (F)			   
	2	Missing BUSCOs (M)			   
	3285	Total BUSC<img width="5000" height="5000" alt="map_dmel-all-chromosome-r6 67_to_hifi_fly_assembly" src="https://github.com/user-attachments/assets/54ce55de-d789-4f7a-84cf-a30b16120258" />

The Shukla assembly is pretty darn good!

# Extra Credit: assembly comparison to contig assembly 
Uploaded target sequence: dmel-all-chromosome-r6.67.fasta
and Query sequence: hifi_fly_assembly.fasta
to the web browser page at: https://dgenies.toulouse.inra.fr/result/eHfe8_20260318205519
<img width="5000" height="5000" alt="map_dmel-all-chromosome-r6 67_to_hifi_fly_assembly" src="https://github.com/user-attachments/assets/72362171-578f-4c44-a6bc-ebe415f03a38" />
