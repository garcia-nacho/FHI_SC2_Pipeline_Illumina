#Corona Swiss-Army-Knife Docker Image
#Nacho Garcia 2021 / iggl@fhi.no

library("rvest")

#pangolintable.path<-"/home/nacho/DockerImages/CommonFiles/PangolinTable.bk.csv"
pangolintable.path<-"/home/docker/CommonFiles/PangolinTable.bk.csv"

results <-"/home/docker/Fastq/"   #Docker
results<-list.dirs(results)
results<-results[grep("summaries$",results)]
results<-paste(results,"/",sep = "")
input.folder<-results   #Docker

#PangoScrapper

try(content <- read_html("https://cov-lineages.org/lineages.html"))
try(tables <- content %>% html_table(fill = TRUE))
try(pangolin.table<-tables[[1]])

warning.pango<-FALSE
if(!exists("pangolin.table")){ 
  pangolin.table<-read.csv(pangolintable.path)
warning.pango<-TRUE
}

names(pangolin.table)<-tolower(names(pangolin.table))

csvs <- list.files(input.folder, full.names = TRUE, pattern = "\\.csv$")
  
  df.pango<-read.csv(csvs[grep(".*NextcladeAndPangolin.csv$", csvs)], sep="\t")
  df.pango[is.na(df.pango)]<-"" #Added 20April 2021 Nacho
  colnames(pangolin.table)[which(colnames(pangolin.table)=="Lineage")]<-"lineage" #Fix Nacho 01Jun2021
  df.pango<-merge(df.pango, pangolin.table, by="lineage",all.x=TRUE) #Nacho 19April2021
  if(warning.pango) df.pango$`most common countries`<-paste("OUTDATED:", df.pango$`most common countries`)
  
  df.summ<-read.csv(csvs[grep(".*summaries.csv$", csvs)], sep = "\t")
  
  df.pango$name<-gsub("/.*","",df.pango$name) #Nanopore fix
  df.summ<- df.summ[-grep("i,a,", df.summ[,2]),] #Nanopore fix
  
  df.summ<-t(df.summ)
  colnames(df.summ)<-df.summ[1,]
  df.summ<-as.data.frame(df.summ)
  df.summ<-df.summ[-1,]
  df.summ$name<-rownames(df.summ)
  
  df.summ$name<-gsub("^X", "", df.summ$name) #Nanopore fix
  #names(df.pango)[which(names(df.pango)=="taxon")]<-"name" #Updated 20May2021

  
  df.summ<-merge(df.summ, df.pango, by="name",all=TRUE) #April 13 Nacho. Added all=TRUE to keep all records
  write.table(df.summ, gsub(".csv","_and_Pangolin.csv",csvs[grep(".*summaries.csv$", csvs)],), sep = "\t", quote = FALSE, row.names = FALSE)
  print("Thanks for using Nacho's Corona Swiss-Army-Knife Docker Image")


