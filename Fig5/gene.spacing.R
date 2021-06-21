setwd("/Users/ben/Downloads/")

library(ggplot2)
library(gridExtra)

##############################
term_cds_lengths <- read.table("Termitomyces_P5_cds.lengths.txt")
  starts <- sort(floor(runif(lengths(term_cds_lengths$V1),1,72e6)))
  ends <- starts + term_cds_lengths$V1
  intergenic_simterm <- c()
  for (i in seq(2,length(ends))){
    intergenic_simterm <- c(starts[i]-ends[i-1],intergenic_simterm)}
intergenic_simterm <- intergenic_simterm[intergenic_simterm>0]

term_spaces <- read.table("Termitomyces_P5_intergenic.spaces.txt")
intergenic_simterm <- sample(intergenic_simterm,length(term_spaces$V1),replace=T)
term <- ggplot() + 
  geom_density(aes(x=intergenic_simterm),lwd=0.3)+
  geom_histogram(aes(x=term_spaces$V1,y=stat(density)),bins=200,fill="tan4",alpha=0.7) + 
  labs(x="")+scale_x_log10(limits=c(1,200000),breaks=c(1e0,1e1,1e2,1e3,1e4,1e5),labels=c()) + 
  theme_classic() + 
  theme(plot.margin=margin(0,0,0,0))+
  #geom_vline(aes(xintercept=mean(term_spaces$V1))) + 
  scale_y_continuous(expand=c(0,0),limits=c(0,1))+
  geom_vline(aes(xintercept=median(term_spaces$V1)),lty=1,lwd=0.5)
##############################

##############################
schi_cds_lengths <- read.table("Schizophyllum_cds.lengths.txt")
starts <- sort(floor(runif(lengths(schi_cds_lengths$V1),1,32e6)))
ends <- starts + schi_cds_lengths$V1
intergenic_simschi <- c()
for (i in seq(2,length(ends))){
  intergenic_simschi <- c(starts[i]-ends[i-1],intergenic_simschi)}
intergenic_simschi <- intergenic_simschi[intergenic_simschi>0]

schi_spaces <- read.table("/Users/ben/Downloads/Schizophyllum_intergenic.spaces.txt")
intergenic_simschi <- sample(intergenic_simschi,length(schi_spaces$V1),replace=T)
schi <- ggplot() + 
  geom_density(aes(x=intergenic_simschi),lwd=0.3)+
  geom_histogram(aes(x=schi_spaces$V1,y=stat(density)),bins=200,fill="orange2",alpha=0.7) + 
  labs(x="Intergenic Space (bp)") + scale_x_log10(limits=c(1,200000),breaks=c(1e0,1e1,1e2,1e3,1e4,1e5),labels=c("1","10","100","1000","10,000","100,000")) + 
  theme_classic() + 
  theme(plot.margin=margin(0,0,0,0))+
  #geom_vline(aes(xintercept=mean(schi_spaces$V1))) + 
  scale_y_continuous(expand=c(0,0),limits=c(0,1))+
  geom_vline(aes(xintercept=median(schi_spaces$V1)),lty=1,lwd=0.5)+
  theme(axis.text.x=element_text(size=12),
        axis.title.x=element_text(size=12),
        plot.margin=margin(0,0,0,0))
##############################

##############################
copr_cds_lengths <- read.table("Coprinopsis_cds.lengths.txt")
starts <- sort(floor(runif(lengths(copr_cds_lengths$V1),1,37e6)))
ends <- starts + copr_cds_lengths$V1
intergenic_simcopr <- c()
for (i in seq(2,length(ends))){
  intergenic_simcopr <- c(starts[i]-ends[i-1],intergenic_simcopr)}
intergenic_simcopr <- intergenic_simcopr[intergenic_simcopr>0]

copr_spaces <- read.table("/Users/ben/Downloads/Coprinopsis_intergenic.spaces.txt")
intergenic_simcopr <- sample(intergenic_simcopr,length(copr_spaces$V1),replace=T)
copr <- ggplot() + 
  geom_density(aes(x=intergenic_simcopr),lwd=0.3)+
  geom_histogram(aes(x=copr_spaces$V1,y=stat(density)),bins=200,fill="cornflowerblue",alpha=0.7) + 
  labs(x="") + scale_x_log10(limits=c(1,200000),breaks=c(1e0,1e1,1e2,1e3,1e4,1e5),labels=c()) + 
  theme_classic() + 
  #geom_vline(aes(xintercept=mean(copr_spaces$V1))) +
  scale_y_continuous(expand=c(0,0),limits=c(0,1))+
  theme(plot.margin=margin(0,0,0,0))+
  geom_vline(aes(xintercept=median(copr_spaces$V1)),lty=1,lwd=0.5)
##############################

##############################
agar_cds_lengths <- read.table("Agaricus_cds.lengths.txt")
starts <- sort(floor(runif(lengths(agar_cds_lengths$V1),1,37e6)))
ends <- starts + agar_cds_lengths$V1
intergenic_simagar <- c()
for (i in seq(2,length(ends))){
  intergenic_simagar <- c(starts[i]-ends[i-1],intergenic_simagar)}
intergenic_simagar <- intergenic_simagar[intergenic_simagar>0]

agar_spaces <- read.table("/Users/ben/Downloads/Agaricus_intergenic.spaces.txt")
intergenic_simagar <- sample(intergenic_simagar,length(agar_spaces$V1),replace=T)
agar <- ggplot() + 
  geom_density(aes(x=intergenic_simagar),lwd=0.3)+
  geom_histogram(aes(x=agar_spaces$V1,y=stat(density)),bins=200,fill="green2",alpha=0.7) + 
  labs(x="") + scale_x_log10(limits=c(1,200000),breaks=c(1e0,1e1,1e2,1e3,1e4,1e5),labels=c()) + 
  theme_classic() + 
  theme(plot.margin=margin(0,0,0,0))+
  #geom_vline(aes(xintercept=mean(agar_spaces$V1))) + 
  scale_y_continuous(expand=c(0,0),limits=c(0,1))+
  geom_vline(aes(xintercept=median(agar_spaces$V1)),lty=1,lwd=0.5)
##############################

##############################
armi_cds_lengths <- read.table("armillaria_cds.lengths.txt")
starts <- sort(floor(runif(lengths(armi_cds_lengths$V1),1,37e6)))
ends <- starts + armi_cds_lengths$V1
intergenic_simarmi <- c()
for (i in seq(2,length(ends))){
  intergenic_simarmi <- c(starts[i]-ends[i-1],intergenic_simarmi)}
intergenic_simarmi <- intergenic_simarmi[intergenic_simarmi>0]

armi_spaces <- read.table("/Users/ben/Downloads/armillaria_intergenic.spaces.txt")
intergenic_simarmi <- sample(intergenic_simarmi,length(armi_spaces$V1),replace=T)
armi <- ggplot() + 
  geom_density(aes(x=intergenic_simarmi),lwd=0.3)+
  geom_histogram(aes(x=armi_spaces$V1,y=stat(density)),bins=200,fill="purple",alpha=0.7) + 
  labs(x="") + scale_x_log10(limits=c(1,200000),breaks=c(1e0,1e1,1e2,1e3,1e4,1e5),labels=c()) + 
  theme_classic() + 
  theme(plot.margin=margin(0,0,0,0))+
  #geom_vline(aes(xintercept=mean(armi_spaces$V1))) + 
  scale_y_continuous(expand=c(0,0),limits=c(0,1))+
  geom_vline(aes(xintercept=median(armi_spaces$V1)),lty=1,lwd=0.5)
##############################

grid.arrange(agar,copr,term,armi,schi,nrow=5)

svg("gene_spacing.svg",width=7,height=4.5)
grid.arrange(agar,copr,term,armi,schi,nrow=5,heights=c(1,1,1,1,1.1))
dev.off()
