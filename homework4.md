# Homework 4: pipelines, plotting, genome assembly
```
wget https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/dmel_r6.66_FB2025_05/fasta/dmel-all-chromosome-r6.66.fasta.gz
```

# Calculate the following for all sequencess up to 100kb and over 100kb
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

# Plots of the following for all sequences <= 100kb and > 100kb
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
p <- ggplot(data, aes(x = Length)) + geom_histogram(fill = "mediumseagreen", color = "black", bins = 50) + theme_minimal() +
labs(title = "Sequence Lengths (100kb or less)", x = "Sequence Length (bp)", y = "Frequency")
ggsave(output_file, plot = p, width = 8, height = 6, dpi = 300)

#For seq over 100kb:
input_file <- "~/Desktop/dmel_over_100kb.txt"
output_file <- "~/Desktop/dmel_over_100kb_histogram.png"
#Define histogram
data <- read.table(input_file, sep = "\t", header = FALSE, col.names = c("Name", "Length"))
p <- ggplot(data, aes(x = Length)) + geom_histogram(fill = "mediumseagreen", color = "black", bins = 50) + theme_minimal() +
labs(title = "Sequence Lengths (100kb or less)", x = "Sequence Length (bp)", y = "Frequency")
ggsave(output_file, plot = p, width = 8, height = 6, dpi = 300)
```
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

#For seq over 100kb:
library(ggplot2)
#Define files
input_file <- "~/Desktop/dmel_gc_over_100kb.txt"
output_file <- "~/Desktop/dmel_gc_over_100kb_histogram.png"
#Define histogram
data <- read.table(input_file, sep = "\t", header = FALSE, col.names = c("Name", "Length", "GC"))
p <- ggplot(data, aes(x = GC)) + geom_histogram(fill = "mediumseagreen", color = "black", bins = 50) + theme_minimal() +
labs(title = "GC Distribution (100kb or less)", x = "GC content", y = "Frequency")
ggsave(output_file, plot = p, width = 8, height = 6, dpi = 300)
```

## 3. Cumulative sequence size sorted from largest to smallest
```
< dmel-all-chromosome-r6.66.fasta.gz \
bioawk -c fastx \
'{ if (length($seq) <= 100000) print length($seq) }' \
| sort -rn \
| gawk 'BEGIN {print "Length\tAssembly"} {print $1 "\tDmel_100kb_or_less"}' > dmel_100kb_or_less_lengths.txt \
| plotCDF2 dmel_100kb_or_less_lengths.txt cumulative_sizes.png
```
```
< dmel-all-chromosome-r6.66.fasta.gz \
bioawk -c fastx \
'{ if (length($seq) <= 100000) print length($seq) }' \
| sort -rn \
| gawk 'BEGIN {print "Length\tAssembly"} {print $1 "\tDmel_100kb_or_less"}' > dmel_100kb_or_less_lengths.txt \
| plotCDF2 dmel_100kb_or_less_lengths.txt cumulative_sizes.png
```

# Assemble a genome using Pacbio HiFi reads
```
cd ~/hw4/fasta
wget https://hpc.oit.uci.edu/~solarese/ee282/iso1_onp_a2_1kb.fastq.gz

#!/usr/bin/env bash
#SBATCH --job-name=hifiasm
#SBATCH --partition=standard

#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=hifiasm_%j.out
#SBATCH --error=hifiasm_%j.err

# Initialize mamba
source /data/homezvol1/mnahmou/miniforge3/etc/profile.d/conda.sh
# Activate your environment
conda activate ee282
# Run hifiasm
hifiasm \
-o hifi_fly_assembly \
-t 16 \
ISO_HiFi_Shukla2025.fasta.gz
```
# Assembly assessment
## 1. Calculate the N50 of the above assembly (~21Mbp)
```
awk '/^S/{print ">"$2"\n"$3}' hifi_fly_assembly.bp.p_ctg.gfa > hifi_fly_assembly.fasta
bioawk -c fastx '{print length($seq)}' hifi_fly_assembly.fasta | sort -rn | awk '{a[NR]=$1; s+=$1} END {for(i=1;i<=NR;i++){c+=a[i]; if(c>=s/2){print "N50: " a[i]; exit}}}'
N50: 21715751
```
## 2. Compare assembly to the contig assembly and scaffold
Reference contig N50: 21.5 Mbp
Assembly contig N50: 21.7 Mbp (good!!)
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

# Initialize mamba
source /data/homezvol1/mnahmou/miniforge3/etc/profile.d/conda.sh
# Activate your environment
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

# Initialize mamba
source /data/homezvol1/mnahmou/miniforge3/etc/profile.d/conda.sh
# Activate your environment
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
	3285	Total BUSCO groups searched	

 The Shukla assembly is pretty darn good!
# Extra Credit: assembly comparison to contig assembly 
Uploaded target sequence: dmel-all-chromosome-r6.67.fasta
and Query sequence: hifi_fly_assembly.fasta
to the web browser page at: https://dgenies.toulouse.inra.fr/result/eHfe8_20260318205519
![dotplot] map_dmel-all-chromosome-r6.67_to_hifi_fly_assembly.png


create a dataframe with long and short sequences
color=length on a histogram

plotCD utilityF? made by jj


hifiasm available on the cluster?

what is the N50?
use fasta splitby N

whats a busco score?


extra credit: d-genies


bioawk -c fastx ` [print length($seq) } ' dmel-all-chromosome-r6... | less

Calculating N50

< dmel-all-chromosome-r6.66.fasta.gz \
bioawk -c fastx \
` { tot=length($seq) + tot ; print length($seq) } END { print tot } ' \
| sort -k1,1rn \ 
| gawk ' NR == 1 { tot = $1 } NR > 1 { cs=$1+cs; print $1 "\t" cs / tot } ' \
| gawk ' $2 >= 0.5 { print $12; exit } ' \
>out.txt

< dmel-all-chromosome-r6.66.fasta.gz \
bioawk -c fastx \
` { tot=length($seq) + tot ; print length($seq) } END { print tot } ' \
| sort -k1,1rn | gawk ' NR == 1 { tot = $1 } NR > 1 { cs=$1+cs; if (cs / tot >= 0.5) { print $1; exit } } 
>out.txt

less `which plotCDF2`

faSplitByN \
< dmel...... /dev/stdout 10 \
| bioawk -c fastx ' { print length($seq) } ' \
| sort k1,1rn | gawk ' BEGIN  { print "Length\tAssembly" } { print $1 "t\dmel_r6_scaff" } ' \
| less
 mkfifo, background it


mkfifo scaff_fifo & 

contig_fifo &

faSplitByN \
< dmel...... /dev/stdout 10 \
| bioawk -c fastx ' { print length($seq) } ' \
| sort k1,1rn | gawk ' BEGIN  { print "Length\tAssembly" } { print $1 "t\dmel_r6_scaff"	} ' \
| less


\tAssembly\n0    the zero starts the plot from zero

