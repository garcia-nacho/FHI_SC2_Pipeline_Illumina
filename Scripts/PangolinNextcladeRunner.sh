#!/bin/bash
cat *.fa* > multifasta.fasta
source activate pangolin
pangolin --update
pangolin /home/docker/Fastq/multifasta.fasta --outfile /home/docker/Fastq/multifasta_pangolin_out.csv
conda deactivate

nextclade --input-fasta /home/docker/Fastq/multifasta.fasta --output-csv /home/docker/Fastq/multifasta_Nextclade.results.csv
nextclade_output_converter.py /home/docker/Fastq/multifasta_Nextclade.results.csv >> /home/docker/Fastq/multifasta_Nextclade.results2.csv

awk -F ',' '{print $1 "," $2 "," $4}' /home/docker/Fastq/multifasta_pangolin_out.csv > pangolin_out.csv
awk -F ';' '{print $1 "," $2}' /home/docker/Fastq/multifasta_Nextclade.results.csv > nextclade_out2.csv

cat nextclade_out2.csv | sed "s/, /\//g" > nextclade_out3.csv && mv nextclade_out3.csv nextclade_out2.csv #ny fra 22.06.21 Kamilla&Nacho

(head -n 1 pangolin_out.csv && tail -n +2 pangolin_out.csv | sort) > pangolin_out_sorted.csv
(head -n 1 nextclade_out2.csv && tail -n +2 nextclade_out2.csv | sort) > nextclade.out2_sorted.csv
(head -n 1 /home/docker/Fastq/multifasta_Nextclade.results2.csv && tail -n +2 /home/docker/Fastq/multifasta_Nextclade.results2.csv | sort) > /home/docker/Fastq/multifasta_Nextclade.results2_sorted.csv

paste -d, /home/docker/Fastq/multifasta_Nextclade.results2_sorted.csv pangolin_out_sorted.csv > NextcladeAndPangolin.out.csv
paste -d, NextcladeAndPangolin.out.csv nextclade.out2_sorted.csv > NextcladeAndPangolin.out2.csv

sed 's/,/\t/g' NextcladeAndPangolin.out2.csv | sed 's/ORF10/ORF10\t/g' > /home/docker/Fastq/multifasta_NextcladeAndPangolin.csv