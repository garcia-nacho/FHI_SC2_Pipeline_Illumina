#if(!require("readxl"))install.packages("readxl")
#if(!require("stringr"))install.packages("stringr")
library("stringr")
library("readxl")
  
  args<-"/home/docker/Fastq"
  
  files <- list.files(args[1],full.names = TRUE,recursive = FALSE)
  dirs <- list.dirs(args[1], full.names = TRUE,recursive = FALSE) 
  dirs <- list.dirs(dirs, full.names = TRUE,recursive = FALSE) 
  dirs <- dirs[grep("Oppsett", dirs)]
  df<- read_excel(files[grep(".xlsx", files)])
  colnames(df)<-df[1,]
  df<-df[-1,]
  letter<-unlist(str_split(dirs,""))
  letter<-letter[-grep("/",letter)]
  letter<-toupper(letter[length(letter)])
  oppsettname<-gsub(".*/","",dirs)
  if(!"Kommentar" %in% names(df)){
  
  if(letter=="A") row.to.get<-c(which(df[,1]=="A1"):which(df[,1]=="H3"))
  if(letter=="B") row.to.get<-c(which(df[,1]=="A4"):which(df[,1]=="H6"))
  if(letter=="C") row.to.get<-c(which(df[,1]=="A7"):which(df[,1]=="H9"))
  if(letter=="D") row.to.get<-c(which(df[,1]=="A10"):which(df[,1]=="H12"))
  
  }else{ #48Cells
    
    if(letter=="A") row.to.get<-grep("Flowcell 1", df$Kommentar)
    if(letter=="B") row.to.get<-grep("Flowcell 2", df$Kommentar)
    if(letter=="C") row.to.get<-grep("Flowcell 3", df$Kommentar)
    if(letter=="D") row.to.get<-grep("Flowcell 4", df$Kommentar)
    
  }
  
  #Changes to use Barcode instead Labware#
  df<-df[row.to.get,c("Barcode","SequenceID")]
  
  if(length(which(is.na(df$Barcode)))>0) df<-df[-which(is.na(df$Barcode)),]
  df<-df[,c(2,1)]
  
  write.table(df, paste(dirs,"/",oppsettname,".tsv",sep = ""), col.names = FALSE, sep = "\t", row.names = FALSE, quote = FALSE)


