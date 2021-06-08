setwd("/Users/ben/Downloads/")

library(ggplot)
library(gridExtra)

dat<-read.csv("distance_comparison.csv",header=T)
head(dat)
#size<-5767440 ## will be relevant

L<-length(dat$N) ## number of windows in which cross-overs were detected

### we will first calculate marker distance. 
### as we can see in tail, the last marker is not listed in position.in.bp
### therefore we extract the positions from the marker info

marker_distance<-c(dat$markerdist)

plot(marker_distance)
hist(log(marker_distance))

#### first lets plot if there is a relationship between marker distance
#### and number of cross-overs

plot(dat$N~log(marker_distance))
### indeed we can check that the high number of cross overs co-occurs
### with larger distances
dat[which(dat$N>10),]

### this is highly correlated
cor.test(dat$N,marker_distance, method="spearman")

### now we assume (H0) that the relative marker distance
### determines the outcome of recombination and want to test this

### we construct a probability distribution, given the marker distance

p_dist<- marker_distance / sum(marker_distance)
### the sum of this should be one
sum(p_dist)

### now we can draw a multinominal random set of cross-overs
sum(dat$N)
sample<-rmultinom(1,sum(dat$N),prob=p_dist)

### now we can test how likely a random sample is and how 
### likely our real sample is
dmultinom(sample,prob=p_dist) # is 0
dmultinom(dat$N,prob=p_dist) # is also 0
### both are very unlikely, but for sure, our cross-over pattern
### seem less likely given the marker distance

### lets make a nice loop


Preal<-dmultinom(dat$N,prob=p_dist) 

Prandom<-c()
for (i in 1:1000){
sample<-rmultinom(1,sum(dat$N)+1,prob=p_dist)
Prandom<-c(Prandom,log10(dmultinom(sample,prob=p_dist)))
hist(Prandom, xlim=c(-125,0))
lines(c(log10(Preal),log10(Preal)), c(0,i),col="red")
}

#### 1000 random samples is enough to see that our
### real distribution is very much non random

### now we have seen this we can take not 1 but 1000 random samples,
### and test for outlier genome positions
### we do thousand subsamples
K<-10000

randomset<-rmultinom(K,sum(dat$N),prob=p_dist)

### for the first number of cross overs is less or more than expected,
### and get the minimum p value for this
min((K-sum(dat$N[1]<randomset[1,]))/K,
1-sum(dat$N[1]>randomset[1,])/K)

p_outlier<-c()
for (i in 1:L)
{
p_outlier<-c(p_outlier, min((K-sum(dat$N[i]<randomset[i,]))/K,
1-sum(dat$N[i]>randomset[i,])/K))
}

p.adj<-p.adjust(p_outlier, method="BH", n=L)

min(p.adj)

### some windows are outliers
plot(p.adj,xlab="intermarker space",ylab="adjusted p value")
### if we zoom in 
plot(p.adj,xlab="intermarker space",ylab="adjusted p value",ylim=c(0,0.05))
### we find which
which(p.adj<0.05)

exp_N<- rowMeans(randomset)

single_cols<-rep("black", L)
single_cols[which(p.adj<0.05 & dat$N>exp_N)]<-"red"
single_cols[which(p.adj<0.05 & dat$N<exp_N)]<-"blue"

MINS<-c()
MAXS<-c()
MEANS<-c()
for (i in 1:L)
{
MAXS<-c(MAXS, max(randomset[i,]))
MINS<-c(MINS, min(randomset[i,]))
MEANS<-c(MEANS,mean(randomset[i,]))
}

single_dat_max<-cbind(MINS,MEANS,MAXS,marker_distance)[order(marker_distance),]
single_dat <- dat
single_marker_distance <- marker_distance
singles_dots <- ggplot() + theme_minimal() +
  labs(x="Distance (Mb)",y="Recombination Rate (cM)",title="A) Marker Pairs") +
  scale_color_manual(values=c("black","blue","red"))+
  scale_x_continuous(labels=c("0.0","0.5","1.0","1.5"),breaks=c(0,0.5e6,1e6,1.5e6))+
  geom_point(aes(y=single_dat$N,x=single_marker_distance,col=single_cols),alpha=0.5)+
  geom_line(aes(y=single_dat_max[,1],x=single_dat_max[,4]),lty=3)+
  geom_line(aes(y=single_dat_max[,2],x=single_dat_max[,4]),lty=1)+
  geom_line(aes(y=single_dat_max[,3],x=single_dat_max[,4]),lty=3)+
  geom_text(aes(x=1e6,y=1,label="Minimum"))+
  geom_text(aes(x=1e6,y=22,label="Average"))+
  geom_text(aes(x=1e6,y=43,label="Maximum"))+
  theme(legend.position="none", axis.line=element_line())

singles_dots

hotspots <- c(which(p.adj<0.05 & dat$N>exp_N))
coldspots <- c(which(p.adj<0.05 & dat$N<exp_N))
hot_and_coldspots <- which(p.adj<0.05)

hotspots_data <- dat[hotspots,]
coldspots_data <- dat[coldspots,]

write.csv(hotspots_data, 'hotspots.csv')
write.csv(coldspots_data, 'coldspots.csv')


#####Je zou ook nog het relatieve aantal cross overs (in proportie per chromosoom) 
####kunnen berekenen en vergelijken met relatieve proportie van marker distance. 
###Voor de maximale cross overs per LG is dit
UNILG <-unique(dat$LG)
ENR<-c()
for (j in 1:length(UNILG)){
  LL<-which(dat$LG==UNILG[j])
  WL<-which(dat$N[LL]==max(dat$N[LL]))
  ENR<-c(ENR,(dat$N[LL[WL]]/ sum(dat$N[LL]))/(dat$markerdist[LL[WL]]/ sum(dat$markerdist[LL])))
}

ratio_CO <-c()
for (j in 1:length(UNILG)){
  LL<-which(dat$LG==UNILG[j])
  WL<-which(dat$N[LL]==max(dat$N[LL]))
  ratio_CO<-c(ratio_CO,(dat$N[LL[WL]]/ sum(dat$N[LL])))
}

ratio_dist<-c()
for (j in 1:length(UNILG)){
  LL<-which(dat$LG==UNILG[j])
  WL<-which(dat$N[LL]==max(dat$N[LL]))
  ratio_dist<-c(ratio_dist,(dat$markerdist[LL[WL]]/ sum(dat$markerdist[LL])))
}

#now for grouped distances----
rm(dat)
rm(dat_max)
dat<-read.csv("distance_comparison_grouped.csv",header=T)

head(dat)

L<-length(dat$N) ## number of windows in which cross-overs were detected

tail(dat)

### we will first calculate marker distance. 
### as we can see in tail, the last marker is not listed in position.in.bp
### therefore we extract the positions from the marker info


marker_distance<-c(dat$markerdist)

plot(marker_distance)
hist(log(marker_distance))

#### first lets plot if there is a relationship between marker distance
#### and number of cross-overs

plot(dat$N~log(marker_distance))
### indeed we can check that the high number of cross overs co-occurs
### with larger distances
dat[which(dat$N>10),]

### this is highly correlated
cor.test(dat$N,marker_distance, method="spearman")

### now we assume (H0) that the relative marker distance
### determines the outcome of recombination and want to test this

### we construct a probability distribution, given the marker distance

p_dist<- marker_distance / sum(marker_distance)
### the sum of this should be one
sum(p_dist)

### now we can draw a multinominal random set of cross-overs
sum(dat$N)
sample<-rmultinom(1,sum(dat$N),prob=p_dist)

### now we can test how likely a random sample is and how 
### likely our real sample is
dmultinom(sample,prob=p_dist) # is 0
dmultinom(dat$N,prob=p_dist) # is also 0
### both are very unlikely, but for sure, our cross-over pattern
### seem less likely given the marker distance

### lets make a nice loop


Preal<-dmultinom(dat$N,prob=p_dist) 

Prandom<-c()
for (i in 1:1000){
  sample<-rmultinom(1,sum(dat$N)+1,prob=p_dist)
  Prandom<-c(Prandom,log10(dmultinom(sample,prob=p_dist)))
  hist(Prandom, xlim=c(-300,0))
  lines(c(log10(Preal),log10(Preal)), c(0,i),col="red")
}

#### 1000 random samples is enough to see that our
### real distribution is very much non random

### now we have seen this we can take not 1 but 1000 random samples,
### and test for outlier genome positions
### we do thousand subsamples
K<-10000

randomset<-rmultinom(K,sum(dat$N),prob=p_dist)

### for the first number of cross overs is less or more than expected,
### and get the minimum p value for this
min((K-sum(dat$N[1]<randomset[1,]))/K,
    1-sum(dat$N[1]>randomset[1,])/K)

p_outlier<-c()
for (i in 1:L)
{
  p_outlier<-c(p_outlier, min((K-sum(dat$N[i]<randomset[i,]))/K,
                              1-sum(dat$N[i]>randomset[i,])/K))
}

p.adj<-p.adjust(p_outlier, method="BH", n=L)

min(p.adj)

### some windows are outliers
plot(p.adj,xlab="intermarker space",ylab="adjusted p value")
### if we zoom in 
plot(p.adj,xlab="intermarker space",ylab="adjusted p value",ylim=c(0,0.05))
### we find which
which(p.adj<0.05)

exp_N<- rowMeans(randomset)



cols<-rep("black", L)
cols[which(p.adj<0.05 & dat$N>exp_N)]<-"red"
cols[which(p.adj<0.05 & dat$N<exp_N)]<-"blue"
plot(dat$N~ log(marker_distance),col=cols,ylab="number of cross-overs",
     xlab="marker distance",main="cold and hot spots")

MINS<-c()
MAXS<-c()
MEANS<-c()
for (i in 1:L)
{
  MAXS<-c(MAXS, max(randomset[i,]))
  MINS<-c(MINS, min(randomset[i,]))
  MEANS<-c(MEANS,mean(randomset[i,]))
}

dat_max<-cbind(MINS,MEANS,MAXS,marker_distance)[order(marker_distance),]

points(dat_max[,1]~log(dat_max[,4]), lty=2, type='l')
points(dat_max[,2]~log(dat_max[,4]), lty=1, type='l')
points(dat_max[,3]~log(dat_max[,4]), lty=2, type='l')

hotspots <- c(which(p.adj<0.05 & dat$N>exp_N))
coldspots <- c(which(p.adj<0.05 & dat$N<exp_N))
hot_and_coldspots <- which(p.adj<0.05)


hotorcold <- rep("hot", length(hotspots))
hotspots_data <- dat[hotspots,]
hotspots_data <- cbind(hotspots_data, hotorcold)
coldspots_data <- dat[coldspots,]
hotorcold <- rep("cold", length(coldspots))
coldspots_data <- cbind(coldspots_data, hotorcold)

hot_and_coldspots_data <- rbind(hotspots_data, coldspots_data)
hot_and_coldspots_data_sort <- hot_and_coldspots_data[order(hot_and_coldspots_data$LG, hot_and_coldspots_data$start),]
write.csv(hot_and_coldspots_data_sort, "hotandcold_MANUAL_new_grouped.csv")

library(pander)
pander(hot_and_coldspots_data_sort)

write.csv(hotspots_data, 'hotspots.csv')
write.csv(coldspots_data, 'coldspots.csv')

#########################plot non log
markerdistMb <- marker_distance/1000000
cols<-rep("Average", L)
cols[which(p.adj<0.05 & dat$N>exp_N)]<-"Hot"
cols[which(p.adj<0.05 & dat$N<exp_N)]<-"Cold"
grouped_dots <- ggplot() + theme_minimal() +
  labs(x="Distance (Mb)",y="",title="B) Grouped Markers",color="Recombination\nRate") +
  scale_color_manual(values=c("black","blue","red"))+
  scale_x_continuous(labels=c("0.0","0.5","1.0","1.5"),breaks=c(0,0.5,1.0,1.5))+ expand_limits(x=0)+
  geom_point(aes(y=dat$N,x=markerdistMb,col=cols),alpha=0.5)+
  geom_line(aes(y=dat_max[,1],x=dat_max[,4]),lty=3)+
  geom_line(aes(y=dat_max[,2],x=dat_max[,4]),lty=1)+
  geom_line(aes(y=dat_max[,3],x=dat_max[,4]),lty=3)+
  geom_text(aes(x=1.2,y=3,label="Minimum"))+
  geom_text(aes(x=1.2,y=26,label="Average"))+
  geom_text(aes(x=1.2,y=48,label="Maximum"))+
  theme(axis.line=element_line())

grouped_dots
MINS<-c()
MAXS<-c()
MEANS<-c()
for (i in 1:L)
{
  MAXS<-c(MAXS, max(randomset[i,]))
  MINS<-c(MINS, min(randomset[i,]))
  MEANS<-c(MEANS,mean(randomset[i,]))
}

dat_max<-cbind(MINS,MEANS,MAXS,markerdistMb)[order(markerdistMb),]

points(dat_max[,1]~dat_max[,4], lty=2, type='l')
points(dat_max[,2]~dat_max[,4], lty=1, type='l')
points(dat_max[,3]~dat_max[,4], lty=2, type='l')
title("B) Grouped Regions")



#####Je zou ook nog het relatieve aantal cross overs (in proportie per chromosoom) 
####kunnen berekenen en vergelijken met relatieve proportie van marker distance. 
###Voor de maximale cross overs per LG is dit
UNILG <-unique(dat$LG)
ENR<-c()
for (j in 1:length(UNILG)){
  LL<-which(dat$LG==UNILG[j])
  WL<-which(dat$N[LL]==max(dat$N[LL]))
  ENR<-c(ENR,(dat$N[LL[WL]]/ sum(dat$N[LL]))/(dat$markerdist[LL[WL]]/ sum(dat$markerdist[LL])))
}

ratio_CO <-c()
for (j in 1:length(UNILG)){
  LL<-which(dat$LG==UNILG[j])
  WL<-which(dat$N[LL]==max(dat$N[LL]))
  ratio_CO<-c(ratio_CO,(dat$N[LL[WL]]/ sum(dat$N[LL])))
}

ratio_dist<-c()
for (j in 1:length(UNILG)){
  LL<-which(dat$LG==UNILG[j])
  WL<-which(dat$N[LL]==max(dat$N[LL]))
  ratio_dist<-c(ratio_dist,(dat$markerdist[LL[WL]]/ sum(dat$markerdist[LL])))
}
svg("Fig3.distance.comparison.svg",width=10,height=4)
grid.arrange(singles_dots,grouped_dots,widths=c(1,1.3),ncol=2)
dev.off()
