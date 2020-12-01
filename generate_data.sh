# firstly, sort and remove duplicates for the bam file

for i in *.bam ; 
do bedtools coverage -abam $i -b /DATABASE/refSeq_gene_TSS_1k.bed > ${i##*/}.bed ;
done

# Script for generate data for Fig.2a-d and Fig.3b
# Using analyze_TSS_coverage.py from Nucleosome_ctDNA (https://github.com/PeterUlz/Nucleosome_ctDNA)
for i in *.bam; 
do python ~/Nucleosome_ctDNA-master/scripts/analyze_TSS_coverage.py -gl top.txt -rg ~/Nucleosome_ctDNA-master/ref/refSeq_extended_names_strand.bed -m 0 -b $i -t 20 > $i.Top.txt;
done
