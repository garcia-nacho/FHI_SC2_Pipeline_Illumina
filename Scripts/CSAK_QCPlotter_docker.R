#Corona Swiss-Army-Knife Docker Image
#Nacho Garcia 2021 / iggl@fhi.no
#V2
#07.Jun.2021

library(ggplot2)
library(readxl)
library(geomnet)
#library(gridExtra)
#library(scales)

#input.folder<-"/home/nacho/Documents/Corona/CSAK_docker_image/Tests/QCPlotter/Run615_Corona__summaries/"
#input.folder<-"C:/Users/IGGL/Desktop/QCPLot/QCPlotter/"
input.folder<-"/home/docker/Fastq/"   #Docker

summaries<-list.files(input.folder, pattern = ".*_summaries_and_Pangolin.csv", full.names = TRUE, recursive = TRUE)
summaries.file<-summaries
name<-gsub("_.*","",gsub(".*/","", summaries))
Plate<- list.files(input.folder, pattern = ".*\\.xlsx", full.names = TRUE, recursive = TRUE)
if(length(Plate[grep("FrameShift",Plate)])==1) Plate<-Plate[-grep("FrameShift",Plate)]
if(length(Plate[grep("NoisExtractor",Plate)])==1) Plate<-Plate[-grep("NoisExtractor",Plate)]
 

Plate<-read_xlsx(Plate, skip = 2)
Plate<-Plate[-which(is.na(Plate$`Well number`)),c("Well number", "Sample ID", names(Plate)[7])]
colnames(Plate)<-c("Position","name","ctvalue")
positions<-expand.grid(toupper(letters[1:8]), c(1:12))
positions<-paste(as.character(positions$Var1), as.character(positions$Var2),sep = "")
pad<-as.data.frame(positions[-which(positions %in% Plate$Position)])
colnames(pad)<-"Position"
pad$name<-NA
pad$ctvalue<-NA
Plate<-rbind(Plate, pad)

summaries<-read.csv(summaries, sep = "\t")

total<- merge(summaries, Plate, all=TRUE, by="name")
total$PX<-as.numeric(gsub("[A-Z]","", total$Position))
total$PY<-gsub("[0-9]","", total$Position)

#reordering factors for plotting
total$PX<-factor(as.character(total$PX), levels = c(1:12))
total$PY<-factor(as.character(total$PY), toupper(letters[c(8:1)]))

#Plot
total$name[which(is.na(total$name))]<-"Empty"

total$coverage.bin<-NA
total$coverage.bin[which(total$Percent.covered.above.depth.9.>=97)]<-">97"
total$coverage.bin[which(total$Percent.covered.above.depth.9.<97)]<-"<97"
total$lineage[which(is.na(total$lineage))]<-"--"
total$ctvalue[which(is.na(total$ctvalue))]<-"--"

ggplot(total, aes(PX, PY, label=paste(gsub("Artic","",name),"\n", lineage, "\n", ctvalue))) +
  geom_circle(aes(x = PX, y=PY, fill=as.numeric(gsub(",",".",Percent.covered.above.depth.9.))), radius=0.05, alpha=0.4)+
  scale_fill_steps(low = "red", high = "green",  breaks=c(97),limits=c(0,100), show.limits=TRUE)+
  #labs(fill = "Coverage")+
  #geom_circle(aes(x = PX, y = PY, colour=as.numeric(gsub(",",".",Percent.covered.above.depth.9.))), radius=0.051)+
  geom_circle(aes(x = PX, y = PY, colour=as.numeric(gsub(",",".",Average.depth.))), radius=0.051)+
  scale_color_gradient2(low = "red", high = "blue", mid = "blue", midpoint = median(as.numeric(gsub(",",".",total$Average.depth.)),na.rm = TRUE))+
  labs(fill = c("Coverage"), color="Average Depth")+
  geom_text(size=3)+
  scale_x_discrete(position = "top")+
    theme_minimal()+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.minor = element_blank(),
        panel.grid.major =element_blank(),axis.title.x = element_blank(), axis.title.y = element_blank(),
        plot.title = element_text(size=22), 
        axis.text =element_text(size=22) )+
  ggtitle(name)
  
ggsave(gsub("_summaries_and_Pangolin.csv","_QCPlot.pdf",summaries.file), width = 18, height = 10)
print("Thanks for using Nacho's Corona Swiss-Army-Knife Docker Image")
