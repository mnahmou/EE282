bioawk to get all seqs up to 100kb and over 100kb

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

