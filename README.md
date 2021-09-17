# FHI SARS-CoV-2 Illumina Pipeline
Bioinformatic pipeline for SARS-CoV-2 sequence analysis used at the [Folkehelseinstituttet](https://www.fhi.no)

## Description
Docker-based solution for sequence analysis of SARS-CoV-2 Illumina samples 

## Protocols supported
ArticV3   
ArticV4   

## Outputs
-Summary including: Mutations found, pangolin lineage, number of reads, covergae, depth, etc...   
-Bam files   
-Consensus sequences   
-Aligned consensus sequences   
-Consensus nucleotide sequence for protein S   
-Indel and frameshift identification   
-Quality-control plot for the plate to detect possible contaminations   
-Phylogenetic-tree plot of the samples analyzed   
-Noise during variant calling accross the genome   
-Quality-control for contaminations/low-quality samples   
-Amplicon efficacy of the selected protocol for all the samples   

## Image building 
<code>git clone https://github.com/garcia-nacho/FHI_SC2_Pipeline_Illumina \   
&& cd FHI_SC2_Pipeline_Illumina \   
&& docker build -t garcianacho/fhisc2:Illumina . </code>
   
Alternativetly, it is posible to *pull* updated builds from dockerhub:

<code>docker pull garcianacho/fhisc2:Illumina</code>

# Running the pipeline


