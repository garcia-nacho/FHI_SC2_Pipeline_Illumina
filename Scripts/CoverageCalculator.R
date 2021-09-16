library("parallel")
library("foreach")
library(doSNOW)
library(progress)
library(Biostrings)
library("ggplot2")
library("pdftools")

#NoiseRunner2

cores<-8
input.folder<-"/home/docker/Fastq/"
#input.folder<-"/media/nacho/Data/DockerImages/Tests/Illumina/Run640kopi/"

results <-paste(input.folder, "temp/",sep = "")   
if(!dir.exists(results)) dir.create(results)

bamfiles<-list.files(input.folder, pattern = ".bam$", full.names = TRUE, recursive = TRUE)
bamfiles<-bamfiles[grep("summaries/bam/", bamfiles)]

samp<-c(1:length(bamfiles))

pb <- progress_bar$new(
  format = "Index: :samp.pb [:bar] :elapsed | eta: :eta",
  total = length(bamfiles),    # 100 
  width = 60)


progress <- function(n){
  pb$tick(tokens = list(samp.pb = samp[n]))
} 

opts <- list(progress = progress)


cluster.cores<-makeCluster(cores)
registerDoSNOW(cluster.cores)

out.par<-foreach(i=1:length(bamfiles), .verbose=FALSE, .options.snow = opts) %dopar%{
  
  try(system(paste("FINex -f ",bamfiles[i], " > ",
                   gsub(".*/",results,gsub("\\.bam", "_NoisExtractorResult2.tsv",bamfiles[i])), sep = "")))
  
}
stopCluster(cluster.cores)


results<-list.files(results,full.names = TRUE)

primers<-read.csv("/home/docker/CommonFiles/artic-ncov2019/primer_schemes/nCoV-2019/V4/SARS-CoV-2.primer.bed",
                  sep = "\t", header = FALSE)
summary<-list.files(input.folder,full.names = TRUE, pattern = "_NextcladeAndPangolin.csv",recursive = TRUE)

summary<-read.csv(summary, sep = "\t")


pb<-txtProgressBar(min = 1, max = length(results), initial = 1)
try(rm(out))
for (i in 1:length(results)) {
  setTxtProgressBar(pb, i)
  dummy<-read.csv(results[i],sep = "\t", header = FALSE)
  dummy<-dummy[,c(1:3)]
  colnames(dummy)<-c("Base","Noise","Reads")
  genome.position<-as.data.frame(c(1:29903))
  colnames(genome.position)<-"Base"
  
  dummy<-merge(genome.position, dummy, by="Base", all=TRUE)
  dummy$Reads[which(is.na(dummy$Reads))]<-0
  
  dummy$Sample<-gsub("_.*","",gsub(".*/","",results[i]))
  if(!exists("out")){
    out<-dummy
  }else{
    out<-rbind(out,dummy)
  }
  
}

primers$Amplicon<-NA
primers$Amplicon[grep("RIGHT", primers$V4)]<-gsub("_RIGHT","",primers$V4[grep("RIGHT", primers$V4)])
primers$Amplicon[grep("LEFT", primers$V4)]<-gsub("_LEFT","",primers$V4[grep("LEFT", primers$V4)])
primers$Side<-gsub(".*_","",primers$V4)

Amplicons<-as.data.frame(unique(primers$Amplicon))
colnames(Amplicons)<-"Amplicon"
Amplicons$Start<-NA
Amplicons$End<-NA
out$Amplicon<-NA

for (i in 1:nrow(Amplicons)) {
  
  Amplicons$Start[i]<-primers$V2[which(primers$Side=="LEFT" & primers$Amplicon==Amplicons$Amplicon[i])]
  Amplicons$End[i]<-primers$V3[which(primers$Side=="RIGHT" & primers$Amplicon==Amplicons$Amplicon[i])]
  
  out$Amplicon[which(out$Base>=Amplicons$Start[i] & out$Base<=Amplicons$End[i])]<-Amplicons$Amplicon[i]
  
}

out.agg<-aggregate(Reads~Amplicon, out, FUN=mean)
out.agg.sd<-aggregate(Reads~Amplicon, out, FUN=sd)
out.agg$SD<-out.agg.sd$Reads
out.agg$Amplicon<-factor(out.agg$Amplicon, levels = 
                           unique(as.character(out.agg$Amplicon))[order(as.numeric(gsub(".*_","",unique(as.character(out.agg$Amplicon)))))])

runname<-gsub("_.*","",gsub(".*/","",list.files(input.folder,full.names = TRUE, pattern = "_NextcladeAndPangolin.csv",recursive = TRUE)))

            ggplot(out.agg)+
            geom_bar(aes(Amplicon, Reads), fill="red", stat = "identity")+
            geom_errorbar(aes(y=Reads, x=Amplicon, ymin=Reads, ymax=Reads+SD ))+
            theme_minimal()+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
            ylab("Raw Number of Reads (mean)")+
              ggtitle(paste("ArticV4 efficiency on", runname))
runid<-list.files(input.folder,full.names = TRUE, pattern = "_NextcladeAndPangolin.csv",recursive = TRUE)
            
ggsave(gsub("_NextcladeAndPangolin.csv","_AmpliconV4.pdf", runid), width = 12, height = 6)
            
            
            
out.agg<-aggregate(Reads~Base, out, FUN=mean)
out.agg.sd<-aggregate(Reads~Base, out, FUN=sd)
out.agg$SD<-out.agg.sd$Reads
out.agg$ymin<-out.agg$Reads-out.agg$SD
out.agg$ymin[which(out.agg$ymin<0)]<-0

ggplot(out.agg)+
  geom_line(aes(Base, Reads), colour="red")+
  geom_ribbon(aes(Base, Reads, ymax=Reads+SD, ymin=ymin), alpha=0.3)+
  theme_minimal()+
  ylab("Raw Number of Reads (mean)")+
  ggtitle(paste("Depth", runname))
ggsave(gsub("_NextcladeAndPangolin.csv","_Depth.pdf", runid), width = 12, height = 6)


colnames(summary)[1]<-"Sample"
out<-merge(out, summary[,c("Sample","lineage")], by="Sample")
out$SampleLineage<-paste(out$Sample, out$lineage, sep = " / ")

out$Amplicon<-factor(out$Amplicon, levels = 
                       unique(as.character(out$Amplicon))[order(as.numeric(gsub(".*_","",unique(as.character(out$Amplicon)))))])

  ggplot(out)+
  geom_bar(aes(Amplicon, Reads), fill="red", stat = "identity")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("Raw Number of Reads (mean)")+
  facet_wrap(~SampleLineage)

ggsave(gsub("_NextcladeAndPangolin.csv","_Amplicon_Sample.pdf", runid), width = 22, height = 16)
  
 ggplot(out)+
   geom_line(aes(Base, Reads), colour="red")+
   theme_minimal()+
   ylab("Raw Number of Reads (mean)")+
   facet_wrap(~SampleLineage)
 
 ggsave(gsub("_NextcladeAndPangolin.csv","_Depth_Sample.pdf", runid), width = 22, height = 16)

 