setwd("/Users/ben/Downloads/")
setwd("/Users/ben/Library/CloudStorage/OneDrive-WageningenUniversity&Research/sabine_genome/other\ files")
#open file filtered data file : 2239 markers, 88 individuals
data <- read.csv("compare.manual_complete.csv")

#open file with contig lengths
TIGlength <- read.table("TIGlength.txt", header = TRUE)

library(ggplot2)
library(gridExtra)
#plot with all data, colored by LG
ggplot(data, aes(Mbp, position)) + geom_point(aes(color = LG)) 

#make all 14 LGs in one plot
#first make the base plot
base <- ggplot() + 
  scale_color_brewer(palette = "Dark2") + 
  theme_classic() +
  guides(size = "none") +
  xlim(0, 8.022900) +
  ylim(0, 158) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        panel.border = element_rect(fill = "NA"),
        legend.position = c(0.8, 0.3),
        legend.text = element_text(size = 7),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, "cm"),
        plot.title = element_text(hjust = 0.5))

LG1 <- base +
  labs(x = "physical position (Mb)", y = "genetic position (cM)", title = "LG1") +
  geom_point(data=subset(data,LG=="LG1"),aes(Mbp,position,color=TIG),shape=1)
LG1

LG2 <- base + 
  labs(x = "physical position (Mb)", title = "LG2") +
  geom_point(data=subset(data,LG=="LG2"),aes(Mbp,position,color=TIG),shape=1)
LG2

#LG3a
LG3a <- base +
  labs(x = "physical position (Mbp)", title = "LG3a") +
  geom_point(data=subset(data,LG=="LG3a"),aes(Mbp,position,color=TIG),shape=1)
LG3a
#LG3b
LG3b <- data[data$LG=="LG3b",]
p4 <- ggplot(LG3b, aes(Mbp, position)) + 
  geom_point(aes(colour = TIG), shape = 1) + 
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "physical position (Mbp)", title = "LG3b") +
  theme_classic() +theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.title.x=element_blank(),axis.text.x=element_blank())+
  guides(size = "none") +
  xlim(0, 8.022900) +
  ylim(0, 158) +
  theme(panel.border = element_rect(fill = "NA"), legend.position = c(0.8, 0.3), legend.text = element_text(size = 7), legend.title = element_blank(), legend.key.size = unit(0.3, "cm"), plot.title = element_text(hjust = 0.5))

#LG4
LG4 <- data[data$LG=="LG4",]
p5 <- ggplot(LG4, aes(Mbp, position)) + 
  geom_point(aes(colour = TIG), shape = 1) + 
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "physical position (Mbp)", title = "LG4") +
  theme_classic() +theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.title.x=element_blank(),axis.text.x=element_blank())+
  guides(size = "none") +
  xlim(0, 8.022900) +
  ylim(0, 158) +
  theme(panel.border = element_rect(fill = "NA"), legend.position = c(0.8, 0.3), legend.text = element_text(size = 7), legend.title = element_blank(), legend.key.size = unit(0.3, "cm"), plot.title = element_text(hjust = 0.5))

#LG5+LG11
LG5_LG11 <- data[data$LG=="LG5_LG11",]
p6 <- ggplot(LG5_LG11, aes(Mbp, position)) + 
  geom_point(aes(colour = TIG), shape = 1) + 
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "physical position (Mbp)", title = "LG5_LG11") +
  theme_classic() +theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.title.x=element_blank(),axis.text.x=element_blank())+
  guides(size = "none") +
  xlim(0, 8.022900) +
  ylim(0, 158) +
  theme(panel.border = element_rect(fill = "NA"), legend.position = c(0.8, 0.3), legend.text = element_text(size = 7), legend.title = element_blank(), legend.key.size = unit(0.3, "cm"), plot.title = element_text(hjust = 0.5))

#LG6
LG6 <- data[data$LG=="LG6",]
p7 <- ggplot(LG6, aes(Mbp, position)) + 
  geom_point(aes(colour = TIG), shape = 1) + 
  scale_color_brewer(palette = "Dark2") + 
  labs(x = " ", y = "genetic position (cM)", title = "LG6") +
  theme_classic() +
  guides(size = "none") +
  xlim(0, 8.022900) +
  ylim(0, 158) +
  theme(panel.border = element_rect(fill = "NA"), legend.position = c(0.8, 0.3), legend.text = element_text(size = 7), legend.title = element_blank(), legend.key.size = unit(0.3, "cm"), plot.title = element_text(hjust = 0.5))

#LG7_LG11
LG7_LG11 <- data[data$LG=="LG7_LG11",]
p8 <- ggplot(LG7_LG11, aes(Mbp, position)) + 
  geom_point(aes(colour = TIG), shape = 1) + 
  scale_color_brewer(palette = "Dark2") + 
  labs(x = " ", title = "LG7_LG8") +
  theme_classic() +theme(axis.title.y=element_blank(),axis.text.y=element_blank())+
  guides(size = "none") +
  xlim(0, 8.022900) +
  ylim(0, 158) +
  theme(panel.border = element_rect(fill = "NA"), legend.position = c(0.8, 0.3), legend.text = element_text(size = 7), legend.title = element_blank(), legend.key.size = unit(0.3, "cm"), plot.title = element_text(hjust = 0.5))

#LG9
LG9 <- data[data$LG=="LG9",]
p9 <- ggplot(LG9, aes(Mbp, position)) + 
  geom_point(aes(colour = TIG), shape = 1) + 
  scale_color_brewer(palette = "Dark2") + 
  labs(x = " ", title = "LG9") +
  theme_classic() +theme(axis.title.y=element_blank(),axis.text.y=element_blank())+
  guides(size = "none") +
  xlim(0, 8.022900) +
  ylim(0, 158) +
  theme(panel.border = element_rect(fill = "NA"), legend.position = c(0.8, 0.3), legend.text = element_text(size = 7), legend.title = element_blank(), legend.key.size = unit(0.3, "cm"), plot.title = element_text(hjust = 0.5))

#LG10_TIG058
LG10_TIG058 <- data[data$LG=="LG10_TIG058",]
p10 <- ggplot(LG10_TIG058, aes(Mbp, position)) + 
  geom_point(aes(colour = TIG), shape = 1) + 
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "physical position (Mb)", title = "LG10_TIG058") +
  theme_classic() +theme(axis.title.y=element_blank(),axis.text.y=element_blank())+
  guides(size = "none") +
  xlim(0, 8.022900) +
  ylim(0, 158) +
  theme(panel.border = element_rect(fill = "NA"), legend.position = c(0.8, 0.3), legend.text = element_text(size = 7), legend.title = element_blank(), legend.key.size = unit(0.3, "cm"), plot.title = element_text(hjust = 0.5))

#LG12
LG12 <- data[data$LG=="LG12",]
p11 <- ggplot(LG12, aes(Mbp, position)) + 
  geom_point(aes(colour = TIG), shape = 1) + 
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "physical position (Mb)", title = "LG12") +
  theme_classic() +theme(axis.title.y=element_blank(),axis.text.y=element_blank())+
  guides(size = "none") +
  xlim(0, 8.022900) +
  ylim(0, 158) +
  theme(panel.border = element_rect(fill = "NA"), legend.position = c(0.8, 0.3), legend.text = element_text(size = 7), legend.title = element_blank(), legend.key.size = unit(0.3, "cm"), plot.title = element_text(hjust = 0.5))

#LG13
LG13 <- data[data$LG=="LG13",]
p12 <- ggplot(LG13, aes(Mbp, position)) + 
  geom_point(aes(colour = TIG), shape = 1) + 
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "physical position (Mb)", y = "genetic position (cM)", title = "LG13") +
  theme_classic() +
  guides(size = "none") +
  xlim(0, 8.022900) +
  ylim(0, 158) +
  theme(panel.border = element_rect(fill = "NA"), legend.position = c(0.8, 0.3), legend.text = element_text(size = 7), legend.title = element_blank(), legend.key.size = unit(0.3, "cm"), plot.title = element_text(hjust = 0.5))

#LG14
LG14 <- data[data$LG=="LG14",]
p13 <- ggplot(LG14, aes(Mbp, position)) + 
  geom_point(aes(colour = TIG), shape = 1) + 
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "physical position (Mb)", title = "LG14") +
  theme_classic() +theme(axis.title.y=element_blank(),axis.text.y=element_blank())+
  guides(size = "none") +
  xlim(0, 8.022900) +
  ylim(0, 158) +
  theme(panel.border = element_rect(fill = "NA"), legend.position = c(0.8, 0.3), legend.text = element_text(size = 7), legend.title = element_blank(), legend.key.size = unit(0.3, "cm"), plot.title = element_text(hjust = 0.5))

#LG15
p14 <- ggplot(subset(data,LG=="LG15"), aes(Mbp, position)) + 
  geom_point(aes(colour = TIG), shape = 1) + 
  scale_color_brewer(palette = "Dark2") + 
  labs(x = "physical position (Mb)", title = "LG15") +
  theme_classic() +theme(axis.title.y=element_blank(),axis.text.y=element_blank())+
  guides(size = "none") +
  xlim(0, 8.022900) +
  ylim(0, 158) +
  theme(panel.border = element_rect(fill = "NA"), legend.position = c(0.8, 0.3), legend.text = element_text(size = 7), legend.title = element_blank(), legend.key.size = unit(0.3, "cm"), plot.title = element_text(hjust = 0.5))

library(gridExtra)
svg("genetic.v.physical.svg",width=10,height=6)
grid.arrange(p1, p2, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, nrow = 3,heights=c(1,1.2,1.2),widths=c(1.2,1,1,1,1))
dev.off()
