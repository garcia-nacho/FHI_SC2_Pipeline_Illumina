# FHI's SARS-CoV-2 Illumina Pipeline
Bioinformatic pipeline for SARS-CoV-2 sequence analysis used at the [Folkehelseinstituttet](https://www.fhi.no)

### Warning!!!!!!!!!
**This repository is longer mantained here.**    
To get the last updated version of it please visit:
[Folkehelseinstuttet's github](https://github.com/folkehelseinstituttet/FHI_SC2_Pipeline_Illumina/)

## Description
Docker-based solution for sequence analysis of SARS-CoV-2 Illumina samples 

## Primer schemes supported
[ArticV3](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V3)   
[ArticV4](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V4)   

## Installation
<code>git clone https://github.com/garcia-nacho/FHI_SC2_Pipeline_Illumina </code>  
<code> cd FHI_SC2_Pipeline_Illumina </code>   
<code> docker build -t garcianacho/fhisc2:Illumina . </code>
 
*Note that building the image for the first time can take up to two hours.* 
 
Alternativetly, it is posible to *pull* updated builds from [Dockerhub](https://hub.docker.com/repository/docker/garcianacho/fhisc2):

<code>docker pull garcianacho/fhisc2:Illumina</code>

## Running the pipeline
*ArticV4:*   
<code>docker run -it --rm -v $(pwd):/home/docker/Fastq garcianacho/fhisc2:Illumina SARS-CoV-2_Illumina_Docker_V12.sh ArticV4</code>    
   
*ArticV3:*   
<code>docker run -it --rm -v $(pwd):/home/docker/Fastq garcianacho/fhisc2:Illumina SARS-CoV-2_Illumina_Docker_V12.sh ArticV3</code>

*Note that older versions of docker might require the flag --privileged and that multiuser systems might require the flag -u 1000 to run*

The script expects the following folder structure where the *fastq.gz* files are placed inside independent folders for each Sample
   
<pre>
./ExpXX    
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

</pre>
   

The script also expects a *.xlsx* file, that contains information about the position of the samples on a 96-well-plate and the DNA concentration (alternatively this column can be used for the Ct-values).
If the file is not properly formated the script will run without errors but the Quality-control plot will not be generated or it will contain errors. 
Note that the script takes the name of the experiment from the name of the xlsx file. If the file is not found the names of the output files might be incorrect. 
It is possible to download a template of the xlsx file [here](https://github.com/garcia-nacho/FHI_SC2_Pipeline_Illumina/blob/master/Template_FHISC2_Illumina.xlsx?raw=true)

## Outputs
-Summary including mutations found, pangolin lineage, number of reads, coverage, depth, etc...   
-Bam files   
-Consensus sequences   
-Aligned consensus sequences   
-Consensus nucleotide sequence for gene *S*   
-*Indels* and frameshift identification   
-Quality-control plot for the plate to detect possible contaminations   
-Phylogenetic-tree plot of the samples   
-Noise during variant calling across the genome   
-Quality-control for contaminations/low-quality samples   
-Amplicon efficacy of the selected primer-set for all the samples   

## Contributors
 
