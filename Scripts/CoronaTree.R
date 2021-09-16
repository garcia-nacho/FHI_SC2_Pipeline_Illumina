library(ggtree)
library(seqinr)
library(phangorn)
library(ape)

fasta.file<-list.files("/home/docker/Fastq/", pattern = ".*aligned.fasta", full.names = TRUE )
dis.matrix<-as.data.frame(as.matrix(dist.alignment(read.alignment(file = fasta.file, format= 'fasta'))))

my_wpgma <- wpgma(dis.matrix)

ggtree(my_wpgma)+   theme_tree()+
  geom_tiplab( size=3)+
  hexpand(.4, direction = 1)

ggsave("/home/docker/Fastq/Tree.pdf", width = 8.27, height = 12)
