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
<code>git clone https://github.com/garcia-nacho/FHI_SC2_Pipeline_Illumina </code>  
<code> cd FHI_SC2_Pipeline_Illumina </code>   
<code> docker build -t garcianacho/fhisc2:Illumina . </code>
   
Alternativetly, it is posible to *pull* updated builds from [dockerhub](https://hub.docker.com/repository/docker/garcianacho/fhisc2):

<code>docker pull garcianacho/fhisc2:Illumina</code>

## Running the pipeline

The script expects this folder structure:
./_    
  |-ExperimentXX.xlsx      
  |-Sample1     
      |-Sample1_SX_LXXXX_R1.fastq.gz       
      |-Sample1_SX_LXXXX_R2.fastq.gz      
  |-Sample2      
      |-Sample2_SX_LXXXX_R1.fastq.gz   
      |-Sample2_SX_LXXXX_R2.fastq.gz   
  |-Sample3   
      |-Sample2_SX_LXXXX_R1.fastq.gz   
      |-Sample2_SX_LXXXX_R2.fastq.gz
  |-...   


Where the fastq.gz files er inside independent folders for each Sample

The scripts expects a .xlsx file, that contains information about the position of the samples on a 96-well-plate and the DNA concentration (alternatively this column can be used for the Ct-values).
If the file is not properly formated the script will run without errors but the QC-plot will not be generated or it will contain errors. 
Note that the scripts takes the name of the experiment from the name of the xlsx file. If the file is not found the names of the output files might be incorrect. 
It is possible to douwload a template of the xlsx file [here](https://github.com/garcia-nacho/FHI_SC2_Pipeline_Illumina/blob/master/Template_FHISC2_Illumina.xlsx?raw=true)

