setwd("B:/Affiliated people/Sabrina_2019/Degree distributions/")
data<-read.csv("Degree Distribution by Phylostratum Proteome.csv")
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
  scale_y_continuous(trans = "log10", breaks = c(1e-08, 1e-06, 1e-04, 1e-02, 1))+
  scale_x_continuous(trans = "log10", breaks = c(1, 10, 100, 1000))+
  expand_limits(x = c(0,2500), y = c(0.0000000001, 1))+
  annotation_logticks()+
  labs(title="Degree distribution of genes associated to phylostratum within the Proteome",
       x="Degree",y="Cumulative Frequency",fill="Cell line")+
  theme(plot.title = element_text(size=22),axis.title = element_text(size=18),axis.text = element_text(size=18),legend.text = element_text(size=12),legend.title = element_text(size=12))
# scale_color_manual(values=c("#56B4E9", "#E69F00", "#999999","#458B00"))+
# geom_vline(xintercept = mean(dls[[1]]$degree),col="#56B4E9",size=1.5)+
# geom_vline(xintercept = mean(dls[[2]]$degree),col="#E69F00",size=1.5)+
# geom_vline(xintercept = mean(dls[[3]]$degree),col="#999999",size=1.5)+
# geom_vline(xintercept = mean(dls[[4]]$degree),col="#458B00",size=1.5)
a
# dev.off()       

wilcox.test(x~phylo_bin,data=dat[which(!is.na(match(as.character(dat$phylo_bin),c("1","8:19")))),])

setwd("B:/Sabrina_2019/")
data<-read.csv("Betweeness Distribution by Phylostratum Proteome.csv")
# data<-read.csv("ICM network for final bridgeness dist.csv")
# data<-read.csv("ICM network for final centrality dist.csv")

library(ggplot2)
library(plyr)
library(reshape2)
colnames(data)<-c("gene","betweeness","Phylostratum")
data$Phylostratum<-factor(data$Phylostratum)
data$phylo_bin<-c(rep("1",length(which(data$Phylostratum==1))),rep("2:3",length(na.omit(match(data$Phylostratum,c(2:3))))),
                  rep("4:7",length(na.omit(match(data$Phylostratum,c(4:7))))),rep("8:19",length(na.omit(match(data$Phylostratum,c(8:19))))))
data$phylo_bin<-factor(data$phylo_bin)
dls<-list()
for(i in 1:length(levels(data$phylo_bin))){
  dls[[i]]<-data[which(data$phylo_bin==levels(data$phylo_bin)[i]),]}

ct<-list()
for(i in 1:length(levels(data$phylo_bin))){
  ct[[i]]<-count(dls[[i]]$betweeness)
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

# png(file="Betweeness Distribution.png",height=860,width=920)
a<-ggplot(dat,aes(x=x,y=cumsum))+
  geom_point(stat="identity",aes(col=phylo_bin),size=2)+
  geom_smooth(se = FALSE, method = "lm",fullrange=TRUE,aes(col=phylo_bin))+
  scale_y_continuous(trans = "log10")+
  scale_x_continuous(trans = "log10")+
  annotation_logticks()+
  labs(title="Betweeness distribution of genes associated to phylostratum within the Proteome",
       x="Betweeness Centrality",y="Cumulative Frequency",fill="Cell line")+
  theme(plot.title = element_text(size=22),axis.title = element_text(size=18),axis.text = element_text(size=18),legend.text = element_text(size=12),legend.title = element_text(size=12))
# scale_color_manual(values=c("#56B4E9", "#E69F00", "#999999","#458B00"))+
# geom_vline(xintercept = mean(dls[[1]]$degree),col="#56B4E9",size=1.5)+
# geom_vline(xintercept = mean(dls[[2]]$degree),col="#E69F00",size=1.5)+
# geom_vline(xintercept = mean(dls[[3]]$degree),col="#999999",size=1.5)+
# geom_vline(xintercept = mean(dls[[4]]$degree),col="#458B00",size=1.5)
a
# dev.off()     

setwd("B:/Sabrina_2019/")
data<-read.csv("Bridgeness Distribution by Phylostratum Proteome.csv")
# data<-read.csv("ICM network for final bridgeness dist.csv")
# data<-read.csv("ICM network for final centrality dist.csv")

library(ggplot2)
library(plyr)
library(reshape2)
colnames(data)<-c("Gene","Bridgeness","Phylostratum")
data$Phylostratum<-factor(data$Phylostratum)
data$phylo_bin<-c(rep("1",length(which(data$Phylostratum==1))),rep("2:3",length(na.omit(match(data$Phylostratum,c(2:3))))),
                  rep("4:7",length(na.omit(match(data$Phylostratum,c(4:7))))),rep("8:19",length(na.omit(match(data$Phylostratum,c(8:19))))))
data$phylo_bin<-factor(data$phylo_bin)
dls<-list()
for(i in 1:length(levels(data$phylo_bin))){
  dls[[i]]<-data[which(data$phylo_bin==levels(data$phylo_bin)[i]),]}

ct<-list()
for(i in 1:length(levels(data$phylo_bin))){
  ct[[i]]<-count(dls[[i]]$Bridgeness)
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
  labs(title="Bridgeness distribution of genes associated to phylostratum within the Proteome",
       x="Bridgeness",y="Cumulative Frequency",fill="Cell line")+
  theme(plot.title = element_text(size=22),axis.title = element_text(size=18),axis.text = element_text(size=18),legend.text = element_text(size=12),legend.title = element_text(size=12))
# scale_color_manual(values=c("#56B4E9", "#E69F00", "#999999","#458B00"))+
# geom_vline(xintercept = mean(dls[[1]]$degree),col="#56B4E9",size=1.5)+
# geom_vline(xintercept = mean(dls[[2]]$degree),col="#E69F00",size=1.5)+
# geom_vline(xintercept = mean(dls[[3]]$degree),col="#999999",size=1.5)+
# geom_vline(xintercept = mean(dls[[4]]$degree),col="#458B00",size=1.5)
a
# dev.off() ,

setwd("B:/Sabrina_2019/")
data<-read.csv("Community Centrality Distribution by Phylostratum Proteome.csv")
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
  labs(title="Community Centrality distribution of genes associated to phylostratum in the Proteome",
       x="Centrality",y="Cumulative Frequency",fill="Cell line")+
  theme(plot.title = element_text(size=22),axis.title = element_text(size=18),axis.text = element_text(size=18),legend.text = element_text(size=12),legend.title = element_text(size=12))
# scale_color_manual(values=c("#56B4E9", "#E69F00", "#999999","#458B00"))+
# geom_vline(xintercept = mean(dls[[1]]$degree),col="#56B4E9",size=1.5)+
# geom_vline(xintercept = mean(dls[[2]]$degree),col="#E69F00",size=1.5)+
# geom_vline(xintercept = mean(dls[[3]]$degree),col="#999999",size=1.5)+
# geom_vline(xintercept = mean(dls[[4]]$degree),col="#458B00",size=1.5)
a
# dev.off() ,

wilcox.test(x~phylo_bin,data=dat[which(!is.na(match(as.character(dat$phylo_bin),c("2:3","4:7")))),])

library(readxl)
Bridgeness_vs_Centrality_Proteome <- read_excel("Bridgeness vs Centrality Proteome.xlsx")

ggplot(Bridgeness_vs_Centrality_Proteome, aes(x=Centrality, y=Bridgeness))+
  geom_point()+
  geom_smooth(method = lm)+
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans = 'log10')+
  ggtitle("Proteome")

library(readxl)
library(ggplot2)
library(plyr)
library(reshape2)

setwd("B:/Sabrina_2019/")
data<-read_excel("Bridgeness vs Centrality Proteome.xlsx")

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
  ggtitle("Proteome")

data<-data.frame(data)
rownames(data)<-as.character(1:2211)
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
  ggtitle("Proteome")

writeClipboard(names(data))
write.table(data, file = "Proteome date party.txt",sep = "\t",row.names = FALSE)

library(readxl)
library(ggplot2)
library(plyr)
library(reshape2)

setwd("B:/Sabrina_2019/")
data<-read_excel("Bridgeness vs Centrality Proteome.xlsx")

colnames(data)<-c("Bridgeness","Centrality","Phylostratum")
data$Phylostratum<-factor(data$Phylostratum)
data$phylo_bin<-c(rep("1",length(which(data$Phylostratum==1))),rep("2:3",length(na.omit(match(data$Phylostratum,c(2:3))))),
                  rep("4:7",length(na.omit(match(data$Phylostratum,c(4:7))))),rep("8:19",length(na.omit(match(data$Phylostratum,c(8:19))))))
data$phylo_bin<-factor(data$phylo_bin)
dls<-list()
for(i in 1:length(levels(data$phylo_bin))){
  dls[[i]]<-data[which(data$phylo_bin==levels(data$phylo_bin)[i]),]}
newdata <- data[ !(data$phylo_bin %in% c("8:19","1","2:3")), ]

ggplot(newdata, aes(x=Centrality, y=Bridgeness))+
  geom_point(stat = "identity", aes(col=phylo_bin))+
  geom_smooth(method = lm)+
  scale_y_continuous(trans='log10')+
  scale_x_continuous(trans = 'log10')+
  ggtitle("Proteome")

