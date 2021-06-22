library(ggplot2)
library(gggenes)
library(RColorBrewer)
library(gridExtra)

setwd("/Users/ben/Desktop/")

genes <- read.csv("mat_loci.csv",stringsAsFactors = F)
genes$gene <- as.character(genes$gene)

genes$new_start[genes$strand == "plus"] <- genes$start[genes$strand == "plus"]
genes$new_end[genes$strand == "plus"] <- genes$end[genes$strand == "plus"]
genes$new_start[genes$strand == "minus"] <- genes$end[genes$strand == "minus"]
genes$new_end[genes$strand == "minus"] <- genes$start[genes$strand == "minus"]

genes <- genes[genes$gene != "2401",]
genes <- genes[genes$gene != "2402",]

genes <- genes[!is.na(genes$gene),]

hdA <- genes[genes$locus == "hd" & genes$start < 5500000,]
hdB <- genes[genes$locus == "hd" & genes$start > 5530000,]

#set x limits to have 7kb sequence
p_hdA <- ggplot(hdA,aes(xmin=new_start,xmax=new_end,y=locus,fill=alt_name)) + theme_genes() + 
  geom_gene_arrow() +
  facet_wrap(~ locus, scales = "free", ncol = 1) + 
  scale_fill_brewer(palette="Set3") +
  scale_x_continuous(limits=c(5455500,5462500))

p_hdB <- ggplot(hdB,aes(xmin=new_start,xmax=new_end,y=locus,fill=alt_name)) + theme_genes() + 
  geom_gene_arrow() +
  facet_wrap(~ locus, scales = "free", ncol = 1) + 
  scale_fill_brewer(palette="Set3") +
  scale_x_continuous(limits=c(5534500,5541500))

p_total <- ggplot(genes,aes(xmin=new_start,xmax=new_end,y=locus,fill=alt_name)) + theme_genes() + 
  geom_gene_arrow() +
  facet_wrap(~locus, scales = "free", ncol = 1) + 
  scale_fill_brewer(palette="Set3")

pdf("P5_mat_loci.pdf")
grid.arrange(p_hdA,p_hdB,p_total,layout_matrix=rbind(c(1,2),c(3,3)),nrow=2)
dev.off()

#this prints linkage group 2 with and without the mating type
#set working directory
setwd("/Users/ben/Downloads/")

library(LinkageMapView)

data1417 <- read.csv("outmap_1417markers_forpicture.csv")
#add mating type
data1417_addMT <- read.csv("outmap_1417markers_addMT_forpicture.csv")

#only LG1
lmv.linkage.plot(data1417, "plot_whole_1417_dupnr_LG2.pdf", mapthese = c(2), dupnbr = TRUE)

#only LG1 MT formated
MTlist <- list()
locus <- c("MatingType")
col <- c("red")
MTlist [[1]] <- list(locus=locus, col=col)
lmv.linkage.plot(data1417_addMT, "plot_whole_1417_addMTformat_dupnr_LG2.pdf", mapthese = c(1), revthese=c(1),dupnbr = TRUE, markerformatlist = MTlist)
