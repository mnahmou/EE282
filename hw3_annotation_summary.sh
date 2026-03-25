#!/bin/bash

#set environment
eval "$(mamba shell hook --shell bash)"
mamba activate ee282
report="$(pwd)/homework3.txt"

#gtf summary
mkdir -p gtf
cd gtf
wget -q https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/dmel_r6.66_FB2025_05/gtf/dmel-all-r6.66.gtf.gz
wget -q https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/dmel_r6.66_FB2025_05/gtf/md5sum.txt
echo "xxxxxMD5 checkumsxxxxx" >> "$report"
grep "dmel-all-r6.66.gtf.gz" md5sum.txt | md5sum -c - | tee md5check.txt >> "$report"
echo "xxxxxFEATURESxxxxx" >> "$report"
bioawk -c gff '{print $feature}' dmel-all-r6.66.gtf.gz | sort | uniq -c | sort -rn | tee summary.txt >> "$report"
bioawk -c gff '$feature == "gene" && $seqname ~ /^(X|Y|2L|2R|3L|3R|4)$/ {print $seqname}' dmel-all-r6.66.gtf.gz | sort | uniq -c | sort -rn | tee -a summary.txt >> "$report"
