library(ggplot2)

setwd("/Users/ben/Downloads/")

gene_density <- read.table("Termitomyces_P5.genesinwindows.bed")
gene_density$V1 <- factor(gene_density$V1,levels=c("scaffold_1","scaffold_2","scaffold_3","scaffold_4","scaffold_5",
                                                   "scaffold_6","scaffold_7","scaffold_8","scaffold_9","scaffold_10",
                                                   "scaffold_11","scaffold_12","scaffold_13","scaffold_14","scaffold_15"))
at <- read.table("Termitomyces_P5.atpercent.bed")
at$V4 <- 1-at$V4
#now need to rescale the AT percentage
#the genes range from 0 to 400, so we want the AT to range approximately from that
at$scaled <- -2500*at$V4
at$scaled <- at$scaled+850
med_at <- median(at$V4)

head(gene_density)
head(at)

range(gene_density$V4)
range(at$scaled)

#plot the entire genomes-----
svg("scaffolds.svg",width=6,height=8)
ggplot(data=gene_density,aes(x=(V2+V3)/2)) + theme_classic() +
  geom_line(aes(y=V4),lwd=0.4)+
  geom_area(aes(y=V4),fill="cornflowerblue") +
  geom_line(aes(y=-75),lwd=0.5) +
  geom_line(aes(y=at$scaled),lwd=0.4) +
  geom_ribbon(aes(ymin=(med_at*-2500)+850,ymax=at$scaled),fill="orange2")+
  scale_y_continuous(breaks=c((0.40*-2500)+850,(0.55*-2500)+850,(0.69*-2500)+850,30,400),labels=c("60","45","30","0","400"))+
  scale_x_continuous(breaks=c(0e6,2e6,4e6,6e6),labels=c(0,2,4,6))+
  labs(x="Position (Mb)",y="GC Content (%)          Gene Density")+
  theme(strip.text=element_text(angle=90),
        panel.spacing.y = unit(0.6,"lines"),
        strip.text.y = element_text(angle=0),
        axis.text.y=element_text(size=5))+
  facet_grid(rows=vars(gene_density$V1))
dev.off()

scaf6_density <- gene_density[gene_density$V1 == "scaffold_6",]
scaf6_at      <- at[at$V1 == "scaffold_6",]

#svg("scaffold_6.svg",width=6,height=2)
#ggplot(data=scaf6_density,aes(x=(V2+V3)/2)) + theme_classic() +
  geom_line(aes(y=V4),lwd=0.4)+
  geom_area(aes(y=V4),fill="cornflowerblue") +
  geom_line(aes(y=-75),lwd=1.5) +
  geom_line(aes(y=scaf6_at$scaled),lwd=0.4) +
  geom_ribbon(aes(ymin=(med_at*-2500)+850,ymax=scaf6_at$scaled),fill="orange2")+
  scale_y_continuous(breaks=c((0.40*-2500)+850,(0.55*-2500)+850,(0.70*-2500)+850,0,400),labels=c("60","45","30","0","400"))+
  scale_x_continuous(breaks=c(0e6,2e6,4e6,6e6),labels=c(0,2,4,6))+
  labs(x="Position (Mb)",y="GC Content (%)          Gene Density")+
  theme(strip.text=element_text(angle=90),
        panel.spacing.y = unit(1.2,"lines"))
#dev.off()
