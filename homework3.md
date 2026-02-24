# Homework 3: Pipelines

**MD5 checkums**
dmel-all-chromosome-r6.66.fasta.gz: OK

**faSize**
143726002 bases (1152978 N's 142573024 real 142573024 upper 0 lower) in 1870 sequences in 1 files
Total size: mean 76858.8 sd 1382100.2 min 544 (211000022279089) max 32079331 (3R) median 1577
N count: mean 616.6 sd 6960.7
U count: mean 76242.3 sd 1379508.4
L count: mean 0.0 sd 0.0
%0.00 masked total, %0.00 masked real

**MD5 checkums**
dmel-all-r6.66.gtf.gz: OK

**Features**
 190176 exon
 163377 CDS
  46856 5UTR
  33778 3UTR
  30922 start_codon
  30862 stop_codon
  30836 mRNA
  17872 gene
   3059 ncRNA
    485 miRNA
    365 pseudogene
    312 tRNA
    270 snoRNA
    262 pre_miRNA
    115 rRNA
     32 snRNA
   4226 3R
   3649 2R
   3508 2L
   3481 3L
   2704 X
    114 4
    113 Y

script used "run.sh":

```
#!/bin/bash

#set environment
eval "$(mamba shell hook --shell bash)"
mamba activate ee282
report="$(pwd)/homework3.txt"

#fasta summary
mkdir -p fasta
cd fasta
wget -q https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/dmel_r6.66_FB2025_05/fasta/dmel-all-chromosome-r6.66>
wget -q https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/dmel_r6.66_FB2025_05/fasta/md5sum.txt
echo "xxxxxMD5 checkumsxxxxx" >> "$report"
grep "dmel-all-chromosome-r6.66.fasta.gz" md5sum.txt | md5sum -c - | tee md5check.txt >> "$report"
echo "xxxxxfaSize outputxxxxx" >> "$report"
faSize dmel-all-chromosome-r6.66.fasta.gz | tee faSize_output.txt >> "$report"

cd ../

#gtf summary
mkdir -p gtf
cd gtf
wget -q https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/dmel_r6.66_FB2025_05/gtf/dmel-all-r6.66.gtf.gz
wget -q https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/dmel_r6.66_FB2025_05/gtf/md5sum.txt
echo "xxxxxMD5 checkumsxxxxx" >> "$report"
grep "dmel-all-r6.66.gtf.gz" md5sum.txt | md5sum -c - | tee md5check.txt >> "$report"
echo "xxxxxFEATURESxxxxx" >> "$report"
bioawk -c gff '{print $feature}' dmel-all-r6.66.gtf.gz | sort | uniq -c | sort -rn | tee summary.txt >> "$report"
bioawk -c gff '$feature == "gene" && $seqname ~ /^(X|Y|2L|2R|3L|3R|4)$/ {print $seqname}' dmel-all-r6.66.gtf.gz | sort>
```
