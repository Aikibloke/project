setwd("B:/Sabrina_2019/Degree Distributions/")
data<-read.csv("Degree Distribution by Phylostratum Ishi transcriptome.csv")
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
  scale_y_continuous(trans = "log10")+
  scale_x_continuous(trans = "log10")+
  annotation_logticks()+
  labs(title="Degree distribution of genes associated to phylostratum within the Ishi transcriptome",
       x="Degree",y="Cumulative Frequency",fill="Cell line")+
  theme(plot.title = element_text(size=22),axis.title = element_text(size=18),axis.text = element_text(size=18),legend.text = element_text(size=12),legend.title = element_text(size=12))
# scale_color_manual(values=c("#56B4E9", "#E69F00", "#999999","#458B00"))+
# geom_vline(xintercept = mean(dls[[1]]$degree),col="#56B4E9",size=1.5)+
# geom_vline(xintercept = mean(dls[[2]]$degree),col="#E69F00",size=1.5)+
# geom_vline(xintercept = mean(dls[[3]]$degree),col="#999999",size=1.5)+
# geom_vline(xintercept = mean(dls[[4]]$degree),col="#458B00",size=1.5)
a
# dev.off() ,


wilcox.test(x~phylo_bin,data=dat[which(!is.na(match(as.character(dat$phylo_bin),c("4:7","8:19")))),])
wilcox.test(x~phylo_bin,data=dat[which(!is.na(match(as.character(dat$phylo_bin),c("1","2:3")))),])
wilcox.test(x~phylo_bin,data=dat[which(!is.na(match(as.character(dat$phylo_bin),c("8:19","2:3")))),])
wilcox.test(x~phylo_bin,data=dat[which(!is.na(match(as.character(dat$phylo_bin),c("4:7","1")))),])
wilcox.test(x~phylo_bin,data=dat[which(!is.na(match(as.character(dat$phylo_bin),c("1","8:19")))),])
wilcox.test(x~phylo_bin,data=dat[which(!is.na(match(as.character(dat$phylo_bin),c("4:7","2:3")))),])

setwd("B:/Sabrina_2019/Bridgeness Distributions/")
data<-read.csv("Bridgeness Distribution by Phylostratum Ishi transcriptome.csv")
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
  labs(title="Bridgeness distribution of genes associated to phylostratum within the Ishi transcriptome",
       x="Bridgeness",y="Cumulative Frequency",fill="Cell line")+
  theme(plot.title = element_text(size=22),axis.title = element_text(size=18),axis.text = element_text(size=18),legend.text = element_text(size=12),legend.title = element_text(size=12))
# scale_color_manual(values=c("#56B4E9", "#E69F00", "#999999","#458B00"))+
# geom_vline(xintercept = mean(dls[[1]]$degree),col="#56B4E9",size=1.5)+
# geom_vline(xintercept = mean(dls[[2]]$degree),col="#E69F00",size=1.5)+
# geom_vline(xintercept = mean(dls[[3]]$degree),col="#999999",size=1.5)+
# geom_vline(xintercept = mean(dls[[4]]$degree),col="#458B00",size=1.5)
a
# dev.off() ,


ggplot(data, aes(x=phylo_bin,y=Community))+
  geom_boxplot(aes(col=phylo_bin),notch = TRUE)+
  scale_y_continuous(trans = "log10")+
  ggtitle("Community centrality score of the surfaceome across phylogentic levels")+xlab("Phylostratum")+ylab("Community Centrality Score")+
  theme_classic()

kruskal.test(Bridgness ~ phylo_bin, data = data)   
pairwise.wilcox.test(data$Bridgness, data$phylo_bin)

setwd("B:/Sabrina_2019/Community Centrality Distribution/")
data<-read.csv("Community Centrality Distribution by Phylostratum Ishi transcriptome.csv")
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
  labs(title="Community Centrality distribution of genes associated to phylostratum in the Ishi transcriptome",
       x="Centrality",y="Cumulative Frequency",fill="Cell line")+
  theme(plot.title = element_text(size=22),axis.title = element_text(size=18),axis.text = element_text(size=18),legend.text = element_text(size=12),legend.title = element_text(size=12))
# scale_color_manual(values=c("#56B4E9", "#E69F00", "#999999","#458B00"))+
# geom_vline(xintercept = mean(dls[[1]]$degree),col="#56B4E9",size=1.5)+
# geom_vline(xintercept = mean(dls[[2]]$degree),col="#E69F00",size=1.5)+
# geom_vline(xintercept = mean(dls[[3]]$degree),col="#999999",size=1.5)+
# geom_vline(xintercept = mean(dls[[4]]$degree),col="#458B00",size=1.5)
a
# dev.off() ,


ggplot(data, aes(x=phylo_bin,y=Community))+
  geom_boxplot(aes(col=phylo_bin),notch = TRUE)+
  scale_y_continuous(trans = "log10")+
  ggtitle("Community centrality score of the surfaceome across phylogentic levels")+xlab("Phylostratum")+ylab("Community Centrality Score")+
  theme_classic()

kruskal.test(Community ~ phylo_bin, data = data)   
pairwise.wilcox.test(data$Community, data$phylo_bin, p.adjust.method = "BH")

library(readxl)                                                                              
setwd("B:/Sabrina_2019/")
data<-read_excel("Bridgeness vs Centrality Ishi transcriptome.xlsx")


library(ggplot2)
library(plyr)
library(reshape2)

colnames(data)<-c("Bridgeness","Centrality","Phylostratum")
data$Phylostratum<-factor(data$Phylostratum)
data$phylo_bin<-c(rep("1",length(which(data$Phylostratum==1))),rep("2:3",length(na.omit(match(data$Phylostratum,c(2:3))))),
                  rep("4:7",length(na.omit(match(data$Phylostratum,c(4:7))))),rep("8:19",length(na.omit(match(data$Phylostratum,c(8:19))))))
data$phylo_bin<-factor(data$phylo_bin)
dls<-list()
for(i in 1:length(levels(data$phylo_bin))){
  dls[[i]]<-data[which(data$phylo_bin==levels(data$phylo_bin)[i]),]}

ggplot(data, aes(x=Centrality, y=Bridgeness))+
  geom_point(stat = "identity", aes(col=phylo_bin))+
  geom_smooth(method = lm)+
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans = 'log10')+
  ggtitle("Ishi transcriptome")

data<-data.frame(data)
rownames(data)<-as.character(1:4853)
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
  ggtitle("Ishi transcriptome")
write.table(data, file = "Ishi transcriptome date party.txt",sep = "\t",row.names = FALSE)

library(readxl)
library(ggplot2)
library(plyr)
library(reshape2)

setwd("B:/Sabrina_2019/")
data<-read_excel("Bridgeness vs Centrality Ishi transcriptome.xlsx")

colnames(data)<-c("Bridgeness","Centrality","Phylostratum")
data$Phylostratum<-factor(data$Phylostratum)
data$phylo_bin<-c(rep("1",length(which(data$Phylostratum==1))),rep("2:3",length(na.omit(match(data$Phylostratum,c(2:3))))),
                  rep("4:7",length(na.omit(match(data$Phylostratum,c(4:7))))),rep("8:19",length(na.omit(match(data$Phylostratum,c(8:19))))))
data$phylo_bin<-factor(data$phylo_bin)
dls<-list()
for(i in 1:length(levels(data$phylo_bin))){
  dls[[i]]<-data[which(data$phylo_bin==levels(data$phylo_bin)[i]),]}
newdata <- data[ !(data$phylo_bin %in% c("8:19","4:7","2:3")), ]

ggplot(newdata, aes(x=Centrality, y=Bridgeness))+
  geom_point(stat = "identity", aes(col=phylo_bin))+
  geom_smooth(method = lm)+
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans = 'log10')+
  ggtitle("Ishi transcriptome")

  
  
  setwd("B:/Sabrina_2019/")
data<-read.csv("Degree Distribution by Phylostratum Ishi e-surf.csv")
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
  scale_y_continuous(trans = "log10")+
  scale_x_continuous(trans = "log10")+
  annotation_logticks()+
  labs(title="Degree distribution of genes associated to phylostratum within the Ishi e-surf",
       x="Degree",y="Cumulative Frequency",fill="Cell line")+
  theme(plot.title = element_text(size=22),axis.title = element_text(size=18),axis.text = element_text(size=18),legend.text = element_text(size=12),legend.title = element_text(size=12))
# scale_color_manual(values=c("#56B4E9", "#E69F00", "#999999","#458B00"))+
# geom_vline(xintercept = mean(dls[[1]]$degree),col="#56B4E9",size=1.5)+
# geom_vline(xintercept = mean(dls[[2]]$degree),col="#E69F00",size=1.5)+
# geom_vline(xintercept = mean(dls[[3]]$degree),col="#999999",size=1.5)+
# geom_vline(xintercept = mean(dls[[4]]$degree),col="#458B00",size=1.5)
a
# dev.off() ,


wilcox.test(x~phylo_bin,data=dat[which(!is.na(match(as.character(dat$phylo_bin),c("4:7","8:19")))),])
wilcox.test(x~phylo_bin,data=dat[which(!is.na(match(as.character(dat$phylo_bin),c("1","2:3")))),])
wilcox.test(x~phylo_bin,data=dat[which(!is.na(match(as.character(dat$phylo_bin),c("8:19","2:3")))),])
wilcox.test(x~phylo_bin,data=dat[which(!is.na(match(as.character(dat$phylo_bin),c("4:7","1")))),])
wilcox.test(x~phylo_bin,data=dat[which(!is.na(match(as.character(dat$phylo_bin),c("1","8:19")))),])
wilcox.test(x~phylo_bin,data=dat[which(!is.na(match(as.character(dat$phylo_bin),c("4:7","2:3")))),])

setwd("B:/Sabrina_2019/")
data<-read.csv("Bridgeness Distribution by Phylostratum Ishi e-surf.csv")
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
  labs(title="Bridgeness distribution of genes associated to phylostratum within the Ishi e-surf",
       x="Bridgeness",y="Cumulative Frequency",fill="Cell line")+
  theme(plot.title = element_text(size=22),axis.title = element_text(size=18),axis.text = element_text(size=18),legend.text = element_text(size=12),legend.title = element_text(size=12))
# scale_color_manual(values=c("#56B4E9", "#E69F00", "#999999","#458B00"))+
# geom_vline(xintercept = mean(dls[[1]]$degree),col="#56B4E9",size=1.5)+
# geom_vline(xintercept = mean(dls[[2]]$degree),col="#E69F00",size=1.5)+
# geom_vline(xintercept = mean(dls[[3]]$degree),col="#999999",size=1.5)+
# geom_vline(xintercept = mean(dls[[4]]$degree),col="#458B00",size=1.5)
a
# dev.off() ,

setwd("B:/Sabrina_2019/Community Centrality Distribution/")
data<-read.csv("Community Centrality Distribution by Phylostratum Ishi e-surf.csv")
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
  labs(title="Community Centrality distribution of genes associated to phylostratum in the Ishi e-surf",
       x="Centrality",y="Cumulative Frequency",fill="Cell line")+
  theme(plot.title = element_text(size=22),axis.title = element_text(size=18),axis.text = element_text(size=18),legend.text = element_text(size=12),legend.title = element_text(size=12))
# scale_color_manual(values=c("#56B4E9", "#E69F00", "#999999","#458B00"))+
# geom_vline(xintercept = mean(dls[[1]]$degree),col="#56B4E9",size=1.5)+
# geom_vline(xintercept = mean(dls[[2]]$degree),col="#E69F00",size=1.5)+
# geom_vline(xintercept = mean(dls[[3]]$degree),col="#999999",size=1.5)+
# geom_vline(xintercept = mean(dls[[4]]$degree),col="#458B00",size=1.5)
a
# dev.off() ,

ggplot(data, aes(x=phylo_bin,y=Community))+
  geom_boxplot(aes(col=phylo_bin))+
  scale_y_continuous(trans = "log10")+
  ggtitle("Community centrality score of the e-surfaceome across phylogentic levels")+xlab("Phylostratum")+ylab("Community Centrality Score")+
  theme_classic()

kruskal.test(Community ~ phylo_bin, data = data)   


library(readxl)                                                                              
setwd("B:/Sabrina_2019/")
data<-read_excel("Bridgeness vs Centrality Ishi e-surf.xlsx")


library(ggplot2)
library(plyr)
library(reshape2)

colnames(data)<-c("Bridgeness","Centrality","Phylostratum")
data$Phylostratum<-factor(data$Phylostratum)
data$phylo_bin<-c(rep("1",length(which(data$Phylostratum==1))),rep("2:3",length(na.omit(match(data$Phylostratum,c(2:3))))),
                  rep("4:7",length(na.omit(match(data$Phylostratum,c(4:7))))),rep("8:19",length(na.omit(match(data$Phylostratum,c(8:19))))))
data$phylo_bin<-factor(data$phylo_bin)
dls<-list()
for(i in 1:length(levels(data$phylo_bin))){
  dls[[i]]<-data[which(data$phylo_bin==levels(data$phylo_bin)[i]),]}

ggplot(data, aes(x=Centrality, y=Bridgeness))+
  geom_point(stat = "identity", aes(col=phylo_bin))+
  geom_smooth(method = lm)+
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans = 'log10')+
  ggtitle("Ishi e-surf")

data<-data.frame(data)
rownames(data)<-as.character(1:4853)
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
  ggtitle("Ishi e-surf")
write.table(data, file = "Ishi e-surf date party.txt",sep = "\t",row.names = FALSE)

library(readxl)
library(ggplot2)
library(plyr)
library(reshape2)

setwd("B:/Sabrina_2019/")
data<-read_excel("Bridgeness vs Centrality Ishi e-surf.xlsx")

colnames(data)<-c("Bridgeness","Centrality","Phylostratum")
data$Phylostratum<-factor(data$Phylostratum)
data$phylo_bin<-c(rep("1",length(which(data$Phylostratum==1))),rep("2:3",length(na.omit(match(data$Phylostratum,c(2:3))))),
                  rep("4:7",length(na.omit(match(data$Phylostratum,c(4:7))))),rep("8:19",length(na.omit(match(data$Phylostratum,c(8:19))))))
data$phylo_bin<-factor(data$phylo_bin)
dls<-list()
for(i in 1:length(levels(data$phylo_bin))){
  dls[[i]]<-data[which(data$phylo_bin==levels(data$phylo_bin)[i]),]}
newdata <- data[ !(data$phylo_bin %in% c("4:7","2:3","1")), ]

ggplot(newdata, aes(x=Centrality, y=Bridgeness))+
  geom_point(stat = "identity", aes(col=phylo_bin))+
  geom_smooth(method = lm)+
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans = 'log10')+
  ggtitle("Ishi e-surf")

setwd("B:/Sabrina_2019/")
data<-read.csv("Degree Distribution by Phylostratum Ishi control.csv")
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
  scale_y_continuous(trans = "log10")+
  scale_x_continuous(trans = "log10")+
  annotation_logticks()+
  labs(title="Degree distribution of genes associated to phylostratum within the Ishi control",
       x="Degree",y="Cumulative Frequency",fill="Cell line")+
  theme(plot.title = element_text(size=22),axis.title = element_text(size=18),axis.text = element_text(size=18),legend.text = element_text(size=12),legend.title = element_text(size=12))
# scale_color_manual(values=c("#56B4E9", "#E69F00", "#999999","#458B00"))+
# geom_vline(xintercept = mean(dls[[1]]$degree),col="#56B4E9",size=1.5)+
# geom_vline(xintercept = mean(dls[[2]]$degree),col="#E69F00",size=1.5)+
# geom_vline(xintercept = mean(dls[[3]]$degree),col="#999999",size=1.5)+
# geom_vline(xintercept = mean(dls[[4]]$degree),col="#458B00",size=1.5)
a
# dev.off() ,

wilcox.test(x~phylo_bin,data=dat[which(!is.na(match(as.character(dat$phylo_bin),c("4:7","8:19")))),])
wilcox.test(x~phylo_bin,data=dat[which(!is.na(match(as.character(dat$phylo_bin),c("1","2:3")))),])
wilcox.test(x~phylo_bin,data=dat[which(!is.na(match(as.character(dat$phylo_bin),c("8:19","2:3")))),])
wilcox.test(x~phylo_bin,data=dat[which(!is.na(match(as.character(dat$phylo_bin),c("4:7","1")))),])
wilcox.test(x~phylo_bin,data=dat[which(!is.na(match(as.character(dat$phylo_bin),c("1","8:19")))),])
wilcox.test(x~phylo_bin,data=dat[which(!is.na(match(as.character(dat$phylo_bin),c("4:7","2:3")))),])



setwd("B:/Sabrina_2019/Bridgeness Distributions")
data<-read.csv("Bridgeness Distribution by HAR Ishi e-surf .csv")
# data<-read.csv("ICM network for final bridgeness dist.csv")
# data<-read.csv("ICM network for final centrality dist.csv")

library(ggplot2)
library(plyr)
library(reshape2)
colnames(data)<-c("present","Gene","Bridgeness","list")

data$present<-factor(data$present)
data$HAR<-c(rep("yes",length(which(data$present==1))),rep("no",length(which(data$present==0))))
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
  labs(title="Bridgeness distribution of human accelerated genes within the Ishikawa e-surfaceome",
       x="Bridgeness",y="Cumulative Frequency",fill="Cell line")+
  theme(plot.title = element_text(size=22),axis.title = element_text(size=18),axis.text = element_text(size=18),legend.text = element_text(size=12),legend.title = element_text(size=12))
# scale_color_manual(values=c("#56B4E9", "#E69F00", "#999999","#458B00"))+
# geom_vline(xintercept = mean(dls[[1]]$degree),col="#56B4E9",size=1.5)+
# geom_vline(xintercept = mean(dls[[2]]$degree),col="#E69F00",size=1.5)+
# geom_vline(xintercept = mean(dls[[3]]$degree),col="#999999",size=1.5)+
# geom_vline(xintercept = mean(dls[[4]]$degree),col="#458B00",size=1.5)
a
# dev.off() ,


setwd("B:/Sabrina_2019/Community Centrality Distribution")
data<-read.csv("Community Centrality Distribution by HAR Ishi e-surf.csv")
# data<-read.csv("ICM network for final bridgeness dist.csv")
# data<-read.csv("ICM network for final centrality dist.csv")

library(ggplot2)
library(plyr)
library(reshape2)
colnames(data)<-c("present","Gene","Community","Phylostratum","List")

data$present<-factor(data$present)
data$HAR<-c(rep("yes",length(which(data$present==1))),rep("no",length(which(data$present==0))))
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
  labs(title="Community Centrality distribution of human accelerated genes within the Ishi e-surf",
       x="Community",y="Cumulative Frequency",fill="Cell line")+
  theme(plot.title = element_text(size=22),axis.title = element_text(size=18),axis.text = element_text(size=18),legend.text = element_text(size=12),legend.title = element_text(size=12))
# scale_color_manual(values=c("#56B4E9", "#E69F00", "#999999","#458B00"))+
# geom_vline(xintercept = mean(dls[[1]]$degree),col="#56B4E9",size=1.5)+
# geom_vline(xintercept = mean(dls[[2]]$degree),col="#E69F00",size=1.5)+
# geom_vline(xintercept = mean(dls[[3]]$degree),col="#999999",size=1.5)+
# geom_vline(xintercept = mean(dls[[4]]$degree),col="#458B00",size=1.5)
a
# dev.off() ,

setwd("B:/Sabrina_2019/Degree Distributions")
data<-read.csv("Degree Distribution by HAR Ishi e-surf.csv")
# data<-read.csv("ICM network for final bridgeness dist.csv")
# data<-read.csv("ICM network for final centrality dist.csv")

library(ggplot2)
library(plyr)
library(reshape2)
colnames(data)<-c("present","Gene","Degree","Phylostratum","List")

data$present<-factor(data$present)
data$HAR<-c(rep("yes",length(which(data$present==1))),rep("no",length(which(data$present==0))))
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
  scale_y_continuous(trans = "log10")+
  scale_x_continuous(trans = "log10")+
  annotation_logticks()+
  labs(title="Degree distribution of human accelerated genes within the Ishi e-surf",
       x="Degree",y="Cumulative Frequency",fill="Cell line")+
  theme(plot.title = element_text(size=22),axis.title = element_text(size=18),axis.text = element_text(size=18),legend.text = element_text(size=12),legend.title = element_text(size=12))
# scale_color_manual(values=c("#56B4E9", "#E69F00", "#999999","#458B00"))+
# geom_vline(xintercept = mean(dls[[1]]$degree),col="#56B4E9",size=1.5)+
# geom_vline(xintercept = mean(dls[[2]]$degree),col="#E69F00",size=1.5)+
# geom_vline(xintercept = mean(dls[[3]]$degree),col="#999999",size=1.5)+
# geom_vline(xintercept = mean(dls[[4]]$degree),col="#458B00",size=1.5)
a
# dev.off() ,

library(readxl)                                                                              
setwd("B:/Sabrina_2019/Bridgeness vs Centrality")
data<-read_excel("Bridgeness vs Centrality HAR Ishi e-surf.xlsx")


library(ggplot2)
library(plyr)
library(reshape2)

colnames(data)<-c("present","Bridgeness","Centrality","Phylostratum","gene","list")
data$present<-factor(data$present)
data$HAR<-c(rep("yes",length(which(data$present==1))),rep("no",length(which(data$present==0))))
data$HAR<-factor(data$HAR)
dls<-list()
for(i in 1:length(levels(data$HAR))){
  dls[[i]]<-data[which(data$HAR==levels(data$HAR)[i]),]}


ggplot(data, aes(x=Centrality, y=Bridgeness))+
  geom_point(stat = "identity", aes(col=HAR))+
  geom_smooth(method = lm)+
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans = 'log10')+
  ggtitle("Har genes Ishikawa e-surfaceome")


setwd("B:/Sabrina_2019/Bridgeness Distributions")
data<-read.csv("Bridgeness Distribution by HAR Ishi transcriptome.csv")
# data<-read.csv("ICM network for final bridgeness dist.csv")
# data<-read.csv("ICM network for final centrality dist.csv")

library(ggplot2)
library(plyr)
library(reshape2)
colnames(data)<-c("present","Gene","Bridgeness","list")

data$present<-factor(data$present)
data$HAR<-c(rep("yes",length(which(data$present==1))),rep("no",length(which(data$present==0))))
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
  labs(title="Bridgeness distribution of human accelerated genes within the Ishikawa transcriptome",
       x="Bridgeness",y="Cumulative Frequency",fill="Cell line")+
  theme(plot.title = element_text(size=22),axis.title = element_text(size=18),axis.text = element_text(size=18),legend.text = element_text(size=12),legend.title = element_text(size=12))
# scale_color_manual(values=c("#56B4E9", "#E69F00", "#999999","#458B00"))+
# geom_vline(xintercept = mean(dls[[1]]$degree),col="#56B4E9",size=1.5)+
# geom_vline(xintercept = mean(dls[[2]]$degree),col="#E69F00",size=1.5)+
# geom_vline(xintercept = mean(dls[[3]]$degree),col="#999999",size=1.5)+
# geom_vline(xintercept = mean(dls[[4]]$degree),col="#458B00",size=1.5)
a
# dev.off() ,


setwd("B:/Sabrina_2019/Community Centrality Distribution")
data<-read.csv("Community Centrality Distribution by HAR Ishi transcriptome.csv")
# data<-read.csv("ICM network for final bridgeness dist.csv")
# data<-read.csv("ICM network for final centrality dist.csv")

library(ggplot2)
library(plyr)
library(reshape2)
colnames(data)<-c("present","Gene","Community","Phylostratum","List")

data$present<-factor(data$present)
data$HAR<-c(rep("yes",length(which(data$present==1))),rep("no",length(which(data$present==0))))
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
  labs(title="Community Centrality distribution of human accelerated genes within the Ishi transcriptome",
       x="Community",y="Cumulative Frequency",fill="Cell line")+
  theme(plot.title = element_text(size=22),axis.title = element_text(size=18),axis.text = element_text(size=18),legend.text = element_text(size=12),legend.title = element_text(size=12))
# scale_color_manual(values=c("#56B4E9", "#E69F00", "#999999","#458B00"))+
# geom_vline(xintercept = mean(dls[[1]]$degree),col="#56B4E9",size=1.5)+
# geom_vline(xintercept = mean(dls[[2]]$degree),col="#E69F00",size=1.5)+
# geom_vline(xintercept = mean(dls[[3]]$degree),col="#999999",size=1.5)+
# geom_vline(xintercept = mean(dls[[4]]$degree),col="#458B00",size=1.5)
a
# dev.off() ,

setwd("B:/Sabrina_2019/Degree Distributions")
data<-read.csv("Degree Distribution by HAR Ishi transcriptome.csv")
# data<-read.csv("ICM network for final bridgeness dist.csv")
# data<-read.csv("ICM network for final centrality dist.csv")

library(ggplot2)
library(plyr)
library(reshape2)
colnames(data)<-c("present","Gene","Degree","Phylostratum","List")

data$present<-factor(data$present)
data$HAR<-c(rep("yes",length(which(data$present==1))),rep("no",length(which(data$present==0))))
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
  scale_y_continuous(trans = "log10")+
  scale_x_continuous(trans = "log10")+
  annotation_logticks()+
  labs(title="Degree distribution of human accelerated genes within the Ishi trancriptome",
       x="Degree",y="Cumulative Frequency",fill="Cell line")+
  theme(plot.title = element_text(size=22),axis.title = element_text(size=18),axis.text = element_text(size=18),legend.text = element_text(size=12),legend.title = element_text(size=12))
# scale_color_manual(values=c("#56B4E9", "#E69F00", "#999999","#458B00"))+
# geom_vline(xintercept = mean(dls[[1]]$degree),col="#56B4E9",size=1.5)+
# geom_vline(xintercept = mean(dls[[2]]$degree),col="#E69F00",size=1.5)+
# geom_vline(xintercept = mean(dls[[3]]$degree),col="#999999",size=1.5)+
# geom_vline(xintercept = mean(dls[[4]]$degree),col="#458B00",size=1.5)
a
# dev.off() ,

library(readxl)                                                                              
setwd("B:/Sabrina_2019/Bridgeness vs Centrality")
data<-read_excel("Bridgeness vs Centrality HAR Ishi transcriptome.xlsx")


library(ggplot2)
library(plyr)
library(reshape2)

colnames(data)<-c("present","Bridgeness","Centrality","Phylostratum","gene","list")
data$present<-factor(data$present)
data$HAR<-c(rep("yes",length(which(data$present==1))),rep("no",length(which(data$present==0))))
data$HAR<-factor(data$HAR)
dls<-list()
for(i in 1:length(levels(data$HAR))){
  dls[[i]]<-data[which(data$HAR==levels(data$HAR)[i]),]}


ggplot(data, aes(x=Centrality, y=Bridgeness))+
  geom_point(stat = "identity", aes(col=HAR))+
  geom_smooth(method = lm)+
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans = 'log10')+
  ggtitle("Har genes Ishikawa transcriptome")

setwd("B:/Sabrina_2019/Degree Distributions")
data<-read.csv("Degree Distribution by HAR Ishi control.csv")
# data<-read.csv("ICM network for final bridgeness dist.csv")
# data<-read.csv("ICM network for final centrality dist.csv")

library(ggplot2)
library(plyr)
library(reshape2)
colnames(data)<-c("present","Gene","Degree","Phylostratum","List")

data$present<-factor(data$present)
data$HAR<-c(rep("yes",length(which(data$present==1))),rep("no",length(which(data$present==0))))
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
  scale_y_continuous(trans = "log10")+
  scale_x_continuous(trans = "log10")+
  annotation_logticks()+
  labs(title="Degree distribution of human accelerated genes within the Ishi control",
       x="Degree",y="Cumulative Frequency",fill="Cell line")+
  theme(plot.title = element_text(size=22),axis.title = element_text(size=18),axis.text = element_text(size=18),legend.text = element_text(size=12),legend.title = element_text(size=12))
# scale_color_manual(values=c("#56B4E9", "#E69F00", "#999999","#458B00"))+
# geom_vline(xintercept = mean(dls[[1]]$degree),col="#56B4E9",size=1.5)+
# geom_vline(xintercept = mean(dls[[2]]$degree),col="#E69F00",size=1.5)+
# geom_vline(xintercept = mean(dls[[3]]$degree),col="#999999",size=1.5)+
# geom_vline(xintercept = mean(dls[[4]]$degree),col="#458B00",size=1.5)
a
# dev.off() ,

wilcox.test(x~HAR,data=dat[which(!is.na(match(as.character(dat$HAR),c("no","yes")))),])
