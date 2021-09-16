#Corona Swiss-Army-Knife Docker Image
#Nacho Garcia 2021 / iggl@fhi.no

library(msa)
library(seqinr)

#Docker

results <-"/home/docker/Fastq/"   #Docker
results<-list.dirs(results)
results<-results[grep("summaries$",results)]
results<-paste(results,"/",sep = "")

result.folder <-results   #Docker
spike <-"/home/docker/CommonFiles/spike.fa"

# # #Testing
# result.folder <-"/home/nacho/DockerImages/DockerBundling/SpikeExtractor/"
# spike <-"/home/nacho/DockerImages/CommonFiles/spike.fa"

args=commandArgs(TRUE)

if(length(args)==0){
  print("Error. No fasta file found!")
}else{
  
  sequences<-readDNAStringSet(c(paste(result.folder,"/fasta/",args[1],sep = ""),spike))
  
  sequences.aln<-list()
  pb<-txtProgressBar(min = 1, max = length(sequences), initial = 1)
  for (i in 1:(length(sequences)-1)) {
    setTxtProgressBar(pb, i)    
    seqs.to.aln<-sequences[c(i, grep("Spike", names(sequences)))]
    alignment<-msa(seqs.to.aln, "Muscle")
    x<-DNAMultipleAlignment(alignment)
    DNAStr = as(x, "DNAStringSet")
    Aligned.samples <-  unlist(base::strsplit(as.character(DNAStr[grep("Spike", names(seqs.to.aln))]),"") )
    
    start.read<-which(Aligned.samples!="-")
    start.read<-start.read[which(start.read==min(start.read))]
    
    end.read<-which(Aligned.samples!="-")
    end.read<-end.read[which(end.read==max(end.read))]
    
    Aligned.samples<-as.character(Aligned.samples[start.read:end.read])
    
    sequences.aln[[i]]<-Aligned.samples
    names(sequences.aln)[i]<-names(DNAStr)[-grep("Spike", names(DNAStr))]
    
  }
  write.fasta(sequences = sequences.aln, names = paste(names(sequences.aln),"_Spike",sep = ""),
              paste(result.folder,args[1],"_Spike.fa",sep = ""))
  print(" ")            
  print("Thanks for using Nacho's Corona Swiss-Army-Knife Docker Image")
}

