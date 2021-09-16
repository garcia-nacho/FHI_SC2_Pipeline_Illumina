#Corona Swiss-Army-Knife Docker Image
#Nacho Garcia 2021 / iggl@fhi.no

library(seqinr)
library(writexl)

args=commandArgs(TRUE)
 cores<-10
if(length(args)==0){ 
  cores<-10
}else{
  cores<-as.numeric(gsub("c","",args[1]))
  }

#Docker
 results <-"/home/docker/Fastq/"   #Docker
 results<-list.dirs(results)
 results<-results[grep("Frameshift",results)]
 results<-paste(results,"/",sep = "")
result.folder <-results   #Docker
common.folder <-"/home/docker/CommonFiles/"
source("/home/docker/Scripts/CSAK_DeletionFinder_v05.R")

# # #Testing
# result.folder <-"/home/nacho/Documents/Corona/CSAK_docker_image/DockerBundling/multisubmit-2021-05-10/"
# common.folder <-"/home/nacho/Documents/Corona/CSAK_docker_image/CommonFiles/"
# source("/home/nacho/Documents/Corona/CSAK_docker_image/Scripts/CSAK_DeletionFinder_v05.R")


multifasta<-list.files(result.folder, pattern = ".*\\.fa.*$", full.names = TRUE)

deletion_results<-DeletionFinder(total.fasta = multifasta,
               results.folder = result.folder,
               cores=cores,
               reference = paste(common.folder,"nCoV-2019.reference.fasta",sep = ""),
               genelist=paste(common.folder,"corona genemap.csv",sep = ""))

genes<-read.csv(paste(common.folder,"corona genemap.csv",sep = ""))
non.codding<-c(1:29903)
for (i in 1:nrow(genes)) {
  non.codding<-non.codding[-which(non.codding %in% c(genes$start[i]:genes$end[i]))]
}


#Cleaning

deletion_results$Frameshift[which(deletion_results$Deletions=="1[28271] / 3[21991] / 6[21765] / 9[11288]" & deletion_results$Insertions=="NO")]<-"NO"

#Check insertions and deletion region of genes
positions.to.test<-list()
for (i in 1:nrow(deletion_results)) {
  if(deletion_results$Insertions[i]=="NO" & deletion_results$Frameshift[i]=="YES"){
    dummy<-unlist(base::strsplit(deletion_results$Deletions[i],"/"))
    del.check<-FALSE
    for (j in 1:length(dummy)) {
      size<-as.numeric(gsub("\\[.*", "",dummy[j]))
      
      if(size%%3 !=0){
        positions.del<-gsub("\\].*","",gsub(".*\\[","",dummy[j]))
        positions.del<-as.numeric(unlist(base::strsplit(positions.del,";")))
        if(length(positions.del)==length(positions.del[which(positions.del %in% non.codding)]) & deletion_results$Insertions[i]=="NO"){
          if(!del.check) deletion_results$Frameshift[i]<-"NO"
          rm(positions.del)
        }else{
          if(length(positions.del[which(positions.del %in% non.codding)])>0) positions.del<-positions.del[-which(positions.del %in% non.codding)]
          if(length(positions.to.test)==0) positions.to.test[[1]]<-positions.del
          if(length(positions.to.test)>0)positions.to.test[[length(positions.to.test)]]<-c(positions.to.test[[length(positions.to.test)]],positions.del)
          deletion_results$Frameshift[i]<-"YES"
          del.check<-TRUE
        }
      }
      
      
    }
    if(deletion_results$Insertions[i]=="NO" & deletion_results$Frameshift[i]=="YES"){
      names(positions.to.test)[length(positions.to.test)]<-deletion_results$Sample[i]
    }  
  }
  
}

deletion_results<-deletion_results[order(deletion_results$Frameshift, decreasing = TRUE),]
write_xlsx(deletion_results[,c(1:4)],paste(result.folder,"FrameShift_", gsub("\\.fa.*","",gsub(".*/","", multifasta)),".xlsx",sep = ""))

print("Thanks for using Nacho's Corona Swiss-Army-Knife Docker Image")
