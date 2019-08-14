setwd("B:/Affiliated people/Sabrina_2019/Degree Distributions/")
data<-read.csv("Degree Distribution by Phylostratum yan polar.csv")
# data<-read.csv("ICM network for final bridgeness dist.csv")
# data<-read.csv("ICM network for final centrality dist.csv")

library(ggplot2)
library(plyr)
library(reshape2)
colnames(data)<-c("gene","degree","Phylostratum")
data$Phylostratum<-factor(data$Phylostratum)
data$phylo_bin<-c(rep("1",length(which(data$Phylostratum==1))),rep("2:3",length(na.omit(match(data$Phylostratum,c(2:3))))),
                  rep("4:7",length(na.omit(match(data$Phylostratum,c(4:7))))),rep("8:19",length(na.omit(match(data$Phylostratum,c(8:19))))))
data$phylo_bin<-factor(data$phylo_bin)
dls<-list()
for(i in 1:length(levels(data$phylo_bin))){
  dls[[i]]<-data[which(data$phylo_bin==levels(data$phylo_bin)[i]),]}

ct<-list()
for(i in 1:length(levels(data$phylo_bin))){
  ct[[i]]<-count(dls[[i]]$degree)
  ct[[i]][,2]<-ct[[i]][,2]/sum(ct[[i]]$freq)
  ct[[i]]<-ct[[i]][order(-ct[[i]]$x),]
  ct[[i]]$cumsum<-cumsum(ct[[i]]$freq)
}

dat<-do.call(rbind.data.frame,ct)
phylo_bin<-rep(levels(data$phylo_bin)[1],length(ct[[1]]$x))

for(i in 2:length(levels(data$phylo_bin))){
  phylo_bin<-c(phylo_bin,rep(levels(data$phylo_bin)[i],length(ct[[i]]$x)))
}
dat$phylo_bin<-phylo_bin

# png(file="Degree Distribution.png",height=860,width=920)
a<-ggplot(dat,aes(x=x,y=cumsum))+
  geom_point(stat="identity",aes(col=phylo_bin),size=2)+
  geom_smooth(se = FALSE, method = "lm",fullrange=TRUE,aes(col=phylo_bin))+
  scale_y_continuous(trans = "log10", breaks = c(0.01, 1.00))+
  scale_x_continuous(trans = "log10", breaks = c(1, 10, 100, 1000))+
  expand_limits(x = c(0,1300), y = c(0.00015,2))+
  annotation_logticks()+
  labs(title="Degree distribution of genes associated to phylostratum within polar",
       x="Degree",y="Cumulative Frequency",fill="Cell line")+
  theme(plot.title = element_text(size=22),axis.title = element_text(size=18),axis.text = element_text(size=18),legend.text = element_text(size=12),legend.title = element_text(size=12))
# scale_color_manual(values=c("#56B4E9", "#E69F00", "#999999","#458B00"))+
# geom_vline(xintercept = mean(dls[[1]]$degree),col="#56B4E9",size=1.5)+
# geom_vline(xintercept = mean(dls[[2]]$degree),col="#E69F00",size=1.5)+
# geom_vline(xintercept = mean(dls[[3]]$degree),col="#999999",size=1.5)+
# geom_vline(xintercept = mean(dls[[4]]$degree),col="#458B00",size=1.5)
a
# dev.off()       


setwd("B:/Affiliated people/Sabrina_2019/Bridgeness Distributions/")
data<-read.csv("Bridgeness Distribution by Phylostratum Yan polar control.csv")
# data<-read.csv("ICM network for final bridgeness dist.csv")
# data<-read.csv("ICM network for final centrality dist.csv")

library(ggplot2)
library(plyr)
library(reshape2)
colnames(data)<-c("Gene","Bridgness","Phylostratum")
data$Phylostratum<-factor(data$Phylostratum)
data$phylo_bin<-c(rep("1",length(which(data$Phylostratum==1))),rep("2:3",length(na.omit(match(data$Phylostratum,c(2:3))))),
                  rep("4:7",length(na.omit(match(data$Phylostratum,c(4:7))))),rep("8:19",length(na.omit(match(data$Phylostratum,c(8:19))))))
data$phylo_bin<-factor(data$phylo_bin)
dls<-list()
for(i in 1:length(levels(data$phylo_bin))){
  dls[[i]]<-data[which(data$phylo_bin==levels(data$phylo_bin)[i]),]}

ct<-list()
for(i in 1:length(levels(data$phylo_bin))){
  ct[[i]]<-count(dls[[i]]$Bridgness)
  ct[[i]][,2]<-ct[[i]][,2]/sum(ct[[i]]$freq)
  ct[[i]]<-ct[[i]][order(-ct[[i]]$x),]
  ct[[i]]$cumsum<-cumsum(ct[[i]]$freq)
}

dat<-do.call(rbind.data.frame,ct)
phylo_bin<-rep(levels(data$phylo_bin)[1],length(ct[[1]]$x))

for(i in 2:length(levels(data$phylo_bin))){
  phylo_bin<-c(phylo_bin,rep(levels(data$phylo_bin)[i],length(ct[[i]]$x)))
}
dat$phylo_bin<-phylo_bin


# png(file="Bridgeness Distribution.png",height=860,width=920)
a<-ggplot(dat,aes(x=x,y=cumsum))+
  geom_point(stat="identity",aes(col=phylo_bin),size=2)+
  geom_smooth(se = FALSE, method = "lm",fullrange=TRUE,aes(col=phylo_bin))+
  scale_y_continuous(trans = "log10")+
  scale_x_continuous(trans = "log10")+
  annotation_logticks()+
  labs(title="Bridgeness distribution of genes associated to phylostratum within the polar control",
       x="Bridgeness",y="Cumulative Frequency",fill="Cell line")+
  theme(plot.title = element_text(size=22),axis.title = element_text(size=18),axis.text = element_text(size=18),legend.text = element_text(size=12),legend.title = element_text(size=12))
# scale_color_manual(values=c("#56B4E9", "#E69F00", "#999999","#458B00"))+
# geom_vline(xintercept = mean(dls[[1]]$degree),col="#56B4E9",size=1.5)+
# geom_vline(xintercept = mean(dls[[2]]$degree),col="#E69F00",size=1.5)+
# geom_vline(xintercept = mean(dls[[3]]$degree),col="#999999",size=1.5)+
# geom_vline(xintercept = mean(dls[[4]]$degree),col="#458B00",size=1.5)
a
# dev.off() ,

mat<-data.frame(matrix(NA,4,1))
colnames(mat)<-c("gradient")
rownames(mat)<-c("1","2:3","4:7","8:19")

dat2<-dat
dat2$x<-log10(dat2$x)
dat2$cumsum<-log10(dat2$cumsum)
test<-lm(cumsum~x,data=dat2[which(dat2$phylo_bin=="1"),])
mat[1,1]<-coef(test)[2]
test<-lm(cumsum~x,data=dat2[which(dat2$phylo_bin=="2:3"),])
mat[2,1]<-coef(test)[2]
test<-lm(cumsum~x,data=dat2[which(dat2$phylo_bin=="4:7"),])
mat[3,1]<-coef(test)[2]
test<-lm(cumsum~x,data=dat2[which(dat2$phylo_bin=="8:19"),])
mat[4,1]<-coef(test)[2]
mat

setwd("B:/Affiliated people/Sabrina_2019/Community Centrality Distribution/")
data<-read.csv("Community Centrality Distribution by Phylostratum Yan polar control.csv")
# data<-read.csv("ICM network for final bridgeness dist.csv")
# data<-read.csv("ICM network for final centrality dist.csv")

library(ggplot2)
library(plyr)
library(reshape2)
colnames(data)<-c("Gene","Community","Phylostratum")
data$Phylostratum<-factor(data$Phylostratum)
data$phylo_bin<-c(rep("1",length(which(data$Phylostratum==1))),rep("2:3",length(na.omit(match(data$Phylostratum,c(2:3))))),
                  rep("4:7",length(na.omit(match(data$Phylostratum,c(4:7))))),rep("8:19",length(na.omit(match(data$Phylostratum,c(8:19))))))
data$phylo_bin<-factor(data$phylo_bin)
dls<-list()
for(i in 1:length(levels(data$phylo_bin))){
  dls[[i]]<-data[which(data$phylo_bin==levels(data$phylo_bin)[i]),]}

ct<-list()
for(i in 1:length(levels(data$phylo_bin))){
  ct[[i]]<-count(dls[[i]]$Community)
  ct[[i]][,2]<-ct[[i]][,2]/sum(ct[[i]]$freq)
  ct[[i]]<-ct[[i]][order(-ct[[i]]$x),]
  ct[[i]]$cumsum<-cumsum(ct[[i]]$freq)
}

dat<-do.call(rbind.data.frame,ct)
phylo_bin<-rep(levels(data$phylo_bin)[1],length(ct[[1]]$x))

for(i in 2:length(levels(data$phylo_bin))){
  phylo_bin<-c(phylo_bin,rep(levels(data$phylo_bin)[i],length(ct[[i]]$x)))
}
dat$phylo_bin<-phylo_bin


# png(file="Community Centrality Distribution.png",height=860,width=920)
a<-ggplot(dat,aes(x=x,y=cumsum))+
  geom_point(stat="identity",aes(col=phylo_bin),size=2)+
  geom_smooth(se = FALSE, method = "lm",fullrange=TRUE,aes(col=phylo_bin))+
  scale_y_continuous(trans = "log10")+
  scale_x_continuous(trans = "log10")+
  annotation_logticks()+
  labs(title="Community Centrality distribution of genes associated to phylostratum in the polar control",
       x="Centrality",y="Cumulative Frequency",fill="Cell line")+
  theme(plot.title = element_text(size=22),axis.title = element_text(size=18),axis.text = element_text(size=18),legend.text = element_text(size=12),legend.title = element_text(size=12))
# scale_color_manual(values=c("#56B4E9", "#E69F00", "#999999","#458B00"))+
# geom_vline(xintercept = mean(dls[[1]]$degree),col="#56B4E9",size=1.5)+
# geom_vline(xintercept = mean(dls[[2]]$degree),col="#E69F00",size=1.5)+
# geom_vline(xintercept = mean(dls[[3]]$degree),col="#999999",size=1.5)+
# geom_vline(xintercept = mean(dls[[4]]$degree),col="#458B00",size=1.5)
a
# dev.off() ,

wilcox.test(x~phylo_bin,data=dat[which(!is.na(match(as.character(dat$phylo_bin),c("4:7","8:19")))),])

library(readxl)
Bridgeness_vs_Centrality_ICM <- read_excel("Bridgeness vs Centrality ICM stiparo.xlsx")

ggplot(Bridgeness_vs_Centrality_ICM, aes(x=Centrality, y=Bridgeness))+
  geom_point()+
  geom_smooth(method = lm)+
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans = 'log10')+
  ggtitle("ICM stiparo")

library(readxl)                                                                              
setwd("B:/Affiliated people/Sabrina_2019/Bridgeness vs Centrality/")
data<-read_excel("Bridgeness vs Centrality yan polar.xlsx")


library(ggplot2)
library(plyr)
library(reshape2)

colnames(data)<-c("Bridgeness","Centrality","Phylostratum", "gene")
data$Phylostratum<-factor(data$Phylostratum)
data$phylo_bin<-c(rep("1",length(which(data$Phylostratum==1))),rep("2:3",length(na.omit(match(data$Phylostratum,c(2:3))))),
                  rep("4:7",length(na.omit(match(data$Phylostratum,c(4:7))))),rep("8:19",length(na.omit(match(data$Phylostratum,c(8:19))))))
data$phylo_bin<-factor(data$phylo_bin)
dls<-list()
for(i in 1:length(levels(data$phylo_bin))){
  dls[[i]]<-data[which(data$phylo_bin==levels(data$phylo_bin)[i]),]}

highlight.gene <- "GATA3"
data$highlight <-ifelse(data$gene == highlight.gene, "highlight", "normal")
textdf <- data[data$gene == highlight.gene]
mycolours <- c("highligh" = "red", "normal" = "grey50")

ggplot(data, aes(x=Centrality, y=Bridgeness))+
  geom_point(stat = "identity", aes(col=highlight))+
  geom_smooth(method = lm)+
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans = 'log10')+
  ggtitle("Polar Yan control")

data<-data.frame(data)
rownames(data)<-as.character(1:1289)
data2<-data
data2$Bridgeness<-log10(data2$Bridgeness)
data2$Centrality<-log10(data2$Centrality)
data2<-data2[which(data2$Bridgeness!="-Inf"),]
test<-lm(Bridgeness~Centrality,data=data2)
inter<-coef(test)[1]
grad<-coef(test)[2]
data2$expected<-NA
for(x in 1:nrow(data2)){
  data2$expected[x]<-grad*(data2$Centrality[x])+inter
}
data2$date<-NA
data2$date[which(data2$Bridgeness>data2$expected)]<-1
data2$date[which(is.na(data2$date))]<-0
data$date<-data2$date[match(rownames(data),rownames(data2))]

ggplot(data, aes(x=Centrality, y=Bridgeness))+
  geom_point(stat = "identity", aes(col=data$date))+
  geom_smooth(method = lm)+
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans = 'log10')+
  ggtitle("Yan polar")
write.table(data, file = "Polar yan surface date party.txt",sep = "\t",row.names = FALSE)


setwd("B:/Affiliated people/Sabrina_2019/Degree Distributions/")
data<-read.csv("Degree Distribution by HAR yan polar surface.csv")
# data<-read.csv("ICM network for final bridgeness dist.csv")
# data<-read.csv("ICM network for final centrality dist.csv")

library(ggplot2)
library(plyr)
library(reshape2)
colnames(data)<-c("present","Gene","Degree","Phylostratum")

data$present<-factor(data$present)
data$HAR<-c(rep("yes",length(which(data$present=="yes"))),rep("no",length(which(data$present=="null"))))
data$HAR<-factor(data$HAR)
dls<-list()
for(i in 1:length(levels(data$HAR))){
  dls[[i]]<-data[which(data$HAR==levels(data$HAR)[i]),]}

ct<-list()
for(i in 1:length(levels(data$HAR))){
  ct[[i]]<-count(dls[[i]]$Degree)
  ct[[i]][,2]<-ct[[i]][,2]/sum(ct[[i]]$freq)
  ct[[i]]<-ct[[i]][order(-ct[[i]]$x),]
  ct[[i]]$cumsum<-cumsum(ct[[i]]$freq)
}

dat<-do.call(rbind.data.frame,ct)
HAR<-rep(levels(data$HAR)[1],length(ct[[1]]$x))

for(i in 2:length(levels(data$HAR))){
  HAR<-c(HAR,rep(levels(data$HAR)[i],length(ct[[i]]$x)))
}
dat$HAR<-HAR



# png(file="Degree Distribution.png",height=860,width=920)
a<-ggplot(dat,aes(x=x,y=cumsum))+
  geom_point(stat="identity",aes(col=HAR),size=2)+
  geom_smooth(se = FALSE, method = "lm",fullrange=TRUE,aes(col=HAR))+
  scale_y_continuous(trans = "log10", breaks = c(0.01, 1.00))+
  scale_x_continuous(trans = "log10", breaks = c(1, 10, 100, 1000))+
  expand_limits(x = c(0,1000), y = c(0.0001,1))+
  annotation_logticks()+
  labs(title="Degree distribution of human accelerated genes within the polar surface",
       x="Degree",y="Cumulative Frequency",fill="Cell line")+
  theme(plot.title = element_text(size=22),axis.title = element_text(size=18),axis.text = element_text(size=18),legend.text = element_text(size=12),legend.title = element_text(size=12))
# scale_color_manual(values=c("#56B4E9", "#E69F00", "#999999","#458B00"))+
# geom_vline(xintercept = mean(dls[[1]]$degree),col="#56B4E9",size=1.5)+
# geom_vline(xintercept = mean(dls[[2]]$degree),col="#E69F00",size=1.5)+
# geom_vline(xintercept = mean(dls[[3]]$degree),col="#999999",size=1.5)+
# geom_vline(xintercept = mean(dls[[4]]$degree),col="#458B00",size=1.5)
a
# dev.off() ,

setwd("B:/Affiliated people/Sabrina_2019/Bridgeness Distributions/")
data<-read.csv("Bridgeness Distribution by HAR yan polar control.csv")
# data<-read.csv("ICM network for final bridgeness dist.csv")
# data<-read.csv("ICM network for final centrality dist.csv")

library(ggplot2)
library(plyr)
library(reshape2)
colnames(data)<-c("present","Gene","Bridgeness","Phylostratum")

data$present<-factor(data$present)
data$HAR<-c(rep("yes",length(which(data$present=="yes"))),rep("no",length(which(data$present=="null"))))
data$HAR<-factor(data$HAR)
dls<-list()
for(i in 1:length(levels(data$HAR))){
  dls[[i]]<-data[which(data$HAR==levels(data$HAR)[i]),]}

ct<-list()
for(i in 1:length(levels(data$HAR))){
  ct[[i]]<-count(dls[[i]]$Bridgeness)
  ct[[i]][,2]<-ct[[i]][,2]/sum(ct[[i]]$freq)
  ct[[i]]<-ct[[i]][order(-ct[[i]]$x),]
  ct[[i]]$cumsum<-cumsum(ct[[i]]$freq)
}

dat<-do.call(rbind.data.frame,ct)
HAR<-rep(levels(data$HAR)[1],length(ct[[1]]$x))

for(i in 2:length(levels(data$HAR))){
  HAR<-c(HAR,rep(levels(data$HAR)[i],length(ct[[i]]$x)))
}
dat$HAR<-HAR



# png(file="Bridgeness Distribution.png",height=860,width=920)
a<-ggplot(dat,aes(x=x,y=cumsum))+
  geom_point(stat="identity",aes(col=HAR),size=2)+
  scale_y_continuous(trans = "log10")+
  scale_x_continuous(trans = "log10")+
  annotation_logticks()+
  labs(title="Bridgeness distribution of human accelerated genes within the polar control",
       x="Bridgeness",y="Cumulative Frequency",fill="Cell line")+
  theme(plot.title = element_text(size=22),axis.title = element_text(size=18),axis.text = element_text(size=18),legend.text = element_text(size=12),legend.title = element_text(size=12))
# scale_color_manual(values=c("#56B4E9", "#E69F00", "#999999","#458B00"))+
# geom_vline(xintercept = mean(dls[[1]]$degree),col="#56B4E9",size=1.5)+
# geom_vline(xintercept = mean(dls[[2]]$degree),col="#E69F00",size=1.5)+
# geom_vline(xintercept = mean(dls[[3]]$degree),col="#999999",size=1.5)+
# geom_vline(xintercept = mean(dls[[4]]$degree),col="#458B00",size=1.5)
a
# dev.off() ,


setwd("B:/Affiliated people/Sabrina_2019/Community Centrality Distribution")
data<-read.csv("Community Centrality Distribution by HAR yan mural surface.csv")
# data<-read.csv("ICM network for final bridgeness dist.csv")
# data<-read.csv("ICM network for final centrality dist.csv")

library(ggplot2)
library(plyr)
library(reshape2)
colnames(data)<-c("present","Gene","Community","Phylostratum")

data$present<-factor(data$present)
data$HAR<-c(rep("yes",length(which(data$present=="yes"))),rep("no",length(which(data$present=="null"))))
data$HAR<-factor(data$HAR)
dls<-list()
for(i in 1:length(levels(data$HAR))){
  dls[[i]]<-data[which(data$HAR==levels(data$HAR)[i]),]}

ct<-list()
for(i in 1:length(levels(data$HAR))){
  ct[[i]]<-count(dls[[i]]$Community)
  ct[[i]][,2]<-ct[[i]][,2]/sum(ct[[i]]$freq)
  ct[[i]]<-ct[[i]][order(-ct[[i]]$x),]
  ct[[i]]$cumsum<-cumsum(ct[[i]]$freq)
}

dat<-do.call(rbind.data.frame,ct)
HAR<-rep(levels(data$HAR)[1],length(ct[[1]]$x))

for(i in 2:length(levels(data$HAR))){
  HAR<-c(HAR,rep(levels(data$HAR)[i],length(ct[[i]]$x)))
}
dat$HAR<-HAR



# png(file="Centrality Distribution.png",height=860,width=920)
a<-ggplot(dat,aes(x=x,y=cumsum))+
  geom_point(stat="identity",aes(col=HAR),size=2)+
  geom_smooth(se = FALSE, method = "lm",fullrange=TRUE,aes(col=HAR))+
  scale_y_continuous(trans = "log10")+
  scale_x_continuous(trans = "log10")+
  annotation_logticks()+
  labs(title="Community Centrality distribution of HAR genes within the mural surface",
       x="Community",y="Cumulative Frequency",fill="Cell line")+
  theme(plot.title = element_text(size=22),axis.title = element_text(size=18),axis.text = element_text(size=18),legend.text = element_text(size=12),legend.title = element_text(size=12))
# scale_color_manual(values=c("#56B4E9", "#E69F00", "#999999","#458B00"))+
# geom_vline(xintercept = mean(dls[[1]]$degree),col="#56B4E9",size=1.5)+
# geom_vline(xintercept = mean(dls[[2]]$degree),col="#E69F00",size=1.5)+
# geom_vline(xintercept = mean(dls[[3]]$degree),col="#999999",size=1.5)+
# geom_vline(xintercept = mean(dls[[4]]$degree),col="#458B00",size=1.5)
a
# dev.off() ,

library(readxl)                                                                              
setwd("B:/Affiliated people/Sabrina_2019/Bridgeness vs Centrality")
data<-read_excel("Bridgeness vs Centrality HAR yan mural surface.xlsx")


library(ggplot2)
library(plyr)
library(reshape2)

colnames(data)<-c("present","Bridgeness","Centrality","Phylostratum","gene")
data$present<-factor(data$present)
data$HAR<-c(rep("yes",length(which(data$present=="yes"))),rep("no",length(which(data$present=="null"))))
data$HAR<-factor(data$HAR)
dls<-list()
for(i in 1:length(levels(data$HAR))){
  dls[[i]]<-data[which(data$HAR==levels(data$HAR)[i]),]}

data<-data.frame(data)
rownames(data)<-as.character(1:1289)
data2<-data
data2$Bridgeness<-log10(data2$Bridgeness)
data2$Centrality<-log10(data2$Centrality)
data2<-data2[which(data2$Bridgeness!="-Inf"),]
test<-lm(Bridgeness~Centrality,data=data2)
inter<-coef(test)[1]
grad<-coef(test)[2]
data2$expected<-NA
for(x in 1:nrow(data2)){
  data2$expected[x]<-grad*(data2$Centrality[x])+inter
}
data2$date<-NA
data2$date[which(data2$Bridgeness>data2$expected)]<-1
data2$date[which(is.na(data2$date))]<-0
data$date<-data2$date[match(rownames(data),rownames(data2))]


ggplot(data, aes(x=Centrality, y=Bridgeness))+
  geom_point(stat = "identity", aes(col=data$date))+
  geom_smooth(method = lm)+
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans = 'log10')+
  ggtitle("mural HAR")

writeClipboard(names(data))
write.table(data, file = "Yan polar surface HAR date party.txt",sep = "\t",row.names = FALSE)

library(readxl)                                                                              
setwd("B:/Affiliated people/Sabrina_2019/Bridgeness vs Centrality")
data<-read_excel("Bridgeness vs Centrality HAR yan polar surface.xlsx")


library(ggplot2)
library(plyr)
library(reshape2)

colnames(data)<-c("present","Bridgeness","Centrality","Phylostratum","gene")
data$present<-factor(data$present)
data$HAR<-c(rep("yes",length(which(data$present=="yes"))),rep("no",length(which(data$present=="null"))))
data$HAR<-factor(data$HAR)
dls<-list()
for(i in 1:length(levels(data$HAR))){
  dls[[i]]<-data[which(data$HAR==levels(data$HAR)[i]),]}


ggplot(data, aes(x=Centrality, y=Bridgeness))+
  geom_point(stat = "identity", aes(col=HAR))+
  geom_smooth(method = lm)+
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans = 'log10')+
  ggtitle("Har genes polar surface")

wilcox.test(x~HAR,data=dat[which(!is.na(match(as.character(dat$HAR),c("no","yes")))),])

wilcox.test(x~phylo_bin,data=dat[which(!is.na(match(as.character(dat$phylo_bin),c("4:7","8:19")))),])
wilcox.test(x~phylo_bin,data=dat[which(!is.na(match(as.character(dat$phylo_bin),c("1","2:3")))),])
wilcox.test(x~phylo_bin,data=dat[which(!is.na(match(as.character(dat$phylo_bin),c("8:19","2:3")))),])
wilcox.test(x~phylo_bin,data=dat[which(!is.na(match(as.character(dat$phylo_bin),c("4:7","1")))),])
wilcox.test(x~phylo_bin,data=dat[which(!is.na(match(as.character(dat$phylo_bin),c("1","8:19")))),])
wilcox.test(x~phylo_bin,data=dat[which(!is.na(match(as.character(dat$phylo_bin),c("4:7","2:3")))),])

younggenes<-matrix(1:4, ncol = 2, nrow = 2, byrow = TRUE)
younggenes
Polar<- c(11,199)
mural<- c(3,216)
trophectoderm<- rbind(Polar, mural)

trophectoderm<- data.frame(Polar, mural)

library(Rcmdr)

setwd("B:/Affiliated people/Sabrina_2019/degree distributions/")
data<-read.csv("degree distribution by phylostratum 8-19 mural vs polar.csv")
# data<-read.csv("ICM network for final bridgeness dist.csv")
# data<-read.csv("ICM network for final centrality dist.csv")

library(ggplot2)
library(plyr)
library(reshape2)
colnames(data)<-c("gene","degree","network")
data$network<-factor(data$network)

dls<-list()
for(i in 1:length(levels(data$network))){
  dls[[i]]<-data[which(data$network==levels(data$network)[i]),]}

ct<-list()
for(i in 1:length(levels(data$network))){
  ct[[i]]<-count(dls[[i]]$degree)
  ct[[i]][,2]<-ct[[i]][,2]/sum(ct[[i]]$freq)
  ct[[i]]<-ct[[i]][order(-ct[[i]]$x),]
  ct[[i]]$cumsum<-cumsum(ct[[i]]$freq)
}

dat<-do.call(rbind.data.frame,ct)
network<-rep(levels(data$network)[1],length(ct[[1]]$x))

for(i in 2:length(levels(data$network))){
  network<-c(network,rep(levels(data$network)[i],length(ct[[i]]$x)))
}
dat$network<-network

# png(file="Degree Distribution.png",height=860,width=920)
a<-ggplot(dat,aes(x=x,y=cumsum))+
  geom_point(stat="identity",aes(col=network),size=2)+
  geom_smooth(se = FALSE, method = "lm",fullrange=TRUE,aes(col=network))+
  scale_y_continuous(trans = "log10", breaks = c(1e-08, 1e-06, 1e-04, 1e-02, 1))+
  scale_x_continuous(trans = "log10", breaks = c(1, 10, 100, 1000))+
  expand_limits(x = c(0,2500), y = c(0.0000000001, 1))+
  annotation_logticks()+
  labs(title="Degree distribution of genes associated to phylostratum within the Surfaceome",
       x="Degree",y="Cumulative Frequency",fill="Cell line")+
  theme(plot.title = element_text(size=22),axis.title = element_text(size=18),axis.text = element_text(size=18),legend.text = element_text(size=12),legend.title = element_text(size=12))
# scale_color_manual(values=c("#56B4E9", "#E69F00", "#999999","#458B00"))+
# geom_vline(xintercept = mean(dls[[1]]$degree),col="#56B4E9",size=1.5)+
# geom_vline(xintercept = mean(dls[[2]]$degree),col="#E69F00",size=1.5)+
# geom_vline(xintercept = mean(dls[[3]]$degree),col="#999999",size=1.5)+
# geom_vline(xintercept = mean(dls[[4]]$degree),col="#458B00",size=1.5)
a
# dev.off() 


wilcox.test(x~network,data=dat[which(!is.na(match(as.character(dat$network),c("MURAL","POLAR")))),])

