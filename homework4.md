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

# Assembly assessment
```
# Get assembly
cd ~/hw4/fasta
wget https://hpc.oit.uci.edu/~solarese/ee282/iso1_onp_a2_1kb.fastq.gz

## 1. Calculate the N50 of the above assembly
## 2. Compare assembly to the contig assembly and scaffold
## 3. Calculate BUSCO scores of both assemblies and compare them

# Extra Credit: assembly comparison to contig assembly 

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

