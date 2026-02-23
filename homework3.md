# Homework 3: Pipelines

## Activate mamba environment
```mamb a activate ee282```

## Download .gz and gff  files from FlyBase
```wget https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/dmel_r6.66_FB2025_05/fasta/index.html \
wget https://s3ftp.flybase.org/genomes/Drosophila_melanogaster/dmel_r6.66_FB2025_05/gtf/index.html```


## faSize to calculate total number of nucelotides, number of Ns, sequences

```faSize *.fasta.gz > faSize_output.txt```

## Verifying file integrity: md5sum

```md5sum * fasta.gz > md5sum.txt```
