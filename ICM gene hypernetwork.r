library("gplots")
setwd("B:\\Terry\\TE_ICM_SEPT2017\\TEvICM_FASEB RESPONSE_2019")
data <- read.table("./Our HESCs/HESC data RMA from Qlucore.txt",sep="\t")
genes<-read.table("./New Analysis and Papers/TE and ICM specific _ defined by stirparo z scores.csv",header=TRUE,sep=",")
te_genes<-as.character(unlist(genes[-which(genes[,2]==""),2]))
icm_genes<-as.character(unlist(genes[-which(genes[,1]==""),1]))

hesc<-data[4:nrow(data),4:ncol(data)]
colnames(hesc)<-as.character(unlist(data[3,4:ncol(data)]))
rownames(hesc)<-as.character(unlist(data[4:nrow(data),1]))

genannot<-data[4:nrow(data),1:2]
samannot<-t(data[1:3,4:ncol(data)])

rn<-rownames(hesc)
hesc<-apply(hesc,2,as.numeric)
rownames(hesc)<-rn

#ICM

MAN1<-hesc[,which(samannot[,1]=="MAN1")]
HUES3<-hesc[,which(samannot[,1]=="HUES3")]
HUES7<-hesc[,which(samannot[,1]=="HUES7")]

MAN1<-cor(t(MAN1[na.omit(match(icm_genes,genannot$V2)),]),t(MAN1[-na.omit(match(icm_genes,genannot$V2)),]))
HUES3<-cor(t(HUES3[na.omit(match(icm_genes,genannot$V2)),]),t(HUES3[-na.omit(match(icm_genes,genannot$V2)),]))
HUES7<-cor(t(HUES7[na.omit(match(icm_genes,genannot$V2)),]),t(HUES7[-na.omit(match(icm_genes,genannot$V2)),]))

MAN1[which(MAN1>=0.8 |MAN1<(-0.8))]<-1
MAN1[which(MAN1!=1)]<-0
MAN1[which(is.na(MAN1))]<-0
HUES3[which(HUES3>=0.8 |HUES3<(-0.8))]<-1
HUES3[which(HUES3!=1)]<-0
HUES3[which(is.na(HUES3))]<-0
HUES7[which(HUES7>=0.8 |HUES7<(-0.8))]<-1
HUES7[which(HUES7!=1)]<-0
HUES7[which(is.na(HUES7))]<-0

man1_mxmt<-(MAN1 %*% t(MAN1))-1
hues3_mxmt<-(HUES3 %*% t(HUES3))-1
hues7_mxmt<-(HUES7 %*% t(HUES7))-1

setwd("R:\\Terry\\TE_ICM_SEPT2017\\TEvICM_FASEB RESPONSE_2019\\ICM gene hypernetwork\\")

png("MAN1_MxMt_Heatmap (ICM).png",height=1000,width=1000)
heatmap.2(man1_mxmt,trace="none")
dev.off()
# hc<-hclust(dist(t(man1_mxmt)))
# plot(hc,cex=0.5)
# png("MAN1_MxMt_Dendrogram (ICM).png",height=1000,width=1000)
# plot(hc,cex=0.5)
# dev.off()
# k<-3
# ct<- cutree(hc, k)
# rect.hclust(hc, k)
# mat<-as.matrix(ct[which(ct==3)])
# mat<-cbind(mat,as.character(unlist(genannot[match(rownames(mat),genannot$V1),2])))
# write.csv(mat,file="MAN1_MxMt_Centralcluster (ICM).csv")

png("HUES3_MxMt_Heatmap (ICM).png",height=1000,width=1000)
heatmap.2(hues3_mxmt,trace="none")
dev.off()
# hc<-hclust(dist(t(hues3_mxmt)))
# plot(hc,cex=0.5)
# png("HUES3_MxMt_Dendrogram (ICM).png",height=1000,width=1000)
# plot(hc,cex=0.5)
# dev.off()
# k<-2
# ct<- cutree(hc, k)
# rect.hclust(hc, k)
# mat<-as.matrix(ct[which(ct==2)])
# mat<-cbind(mat,as.character(unlist(genannot[match(rownames(mat),genannot$V1),2])))
# write.csv(mat,file="HUES3_MxMt_Centralcluster (ICM).csv")

png("HUES7_MxMt_Heatmap (ICM).png",height=1000,width=1000)
heatmap.2(hues7_mxmt,trace="none")
dev.off()
# hc<-hclust(dist(t(hues7_mxmt)))
# plot(hc,cex=0.5)
# png("HUES7_MxMt_Dendrogram (ICM).png",height=1000,width=1000)
# plot(hc,cex=0.5)
# dev.off()
# k<-3
# ct<- cutree(hc, k)
# rect.hclust(hc, k)
# mat<-as.matrix(ct[which(ct==3)])
# mat<-cbind(mat,as.character(unlist(genannot[match(rownames(mat),genannot$V1),2])))
# write.csv(mat,file="HUES7_MxMt_Centralcluster (ICM).csv")

#TE

MAN1<-hesc[,which(samannot[,1]=="MAN1")]
HUES3<-hesc[,which(samannot[,1]=="HUES3")]
HUES7<-hesc[,which(samannot[,1]=="HUES7")]

MAN1<-cor(t(MAN1[na.omit(match(te_genes,genannot$V2)),]),t(MAN1[-na.omit(match(te_genes,genannot$V2)),]))
HUES3<-cor(t(HUES3[na.omit(match(te_genes,genannot$V2)),]),t(HUES3[-na.omit(match(te_genes,genannot$V2)),]))
HUES7<-cor(t(HUES7[na.omit(match(te_genes,genannot$V2)),]),t(HUES7[-na.omit(match(te_genes,genannot$V2)),]))

MAN1[which(MAN1>=0.8 |MAN1<(-0.8))]<-1
MAN1[which(MAN1!=1)]<-0
MAN1[which(is.na(MAN1))]<-0
HUES3[which(HUES3>=0.8 |HUES3<(-0.8))]<-1
HUES3[which(HUES3!=1)]<-0
HUES3[which(is.na(HUES3))]<-0
HUES7[which(HUES7>=0.8 |HUES7<(-0.8))]<-1
HUES7[which(HUES7!=1)]<-0
HUES7[which(is.na(HUES7))]<-0

man1_mxmt<-(MAN1 %*% t(MAN1))-1
hues3_mxmt<-(HUES3 %*% t(HUES3))-1
hues7_mxmt<-(HUES7 %*% t(HUES7))-1

setwd("R:\\Terry\\TE_ICM_SEPT2017\\TEvICM_FASEB RESPONSE_2019\\ICM gene hypernetwork\\")

png("MAN1_MxMt_Heatmap (TE).png",height=1000,width=1000)
heatmap.2(man1_mxmt,trace="none")
dev.off()
# hc<-hclust(dist(t(man1_mxmt)))
# png("MAN1_MxMt_Dendrogram (TE).png",height=1000,width=1000)
# plot(hc,cex=0.5)
# dev.off()
# k<-3
# ct<- cutree(hc, k)
# rect.hclust(hc, k)
# mat<-as.matrix(ct[which(ct==3)])
# mat<-cbind(mat,as.character(unlist(genannot[match(rownames(mat),genannot$V1),2])))
# write.csv(mat,file="MAN1_MxMt_Centralcluster (TE).csv")

png("HUES3_MxMt_Heatmap (TE).png",height=1000,width=1000)
heatmap.2(hues3_mxmt,trace="none")
dev.off()
# hc<-hclust(dist(t(hues3_mxmt)))
# png("HUES3_MxMt_Dendrogram (TE).png",height=1000,width=1000)
# plot(hc,cex=0.5)
# dev.off()
# k<-2
# ct<- cutree(hc, k)
# rect.hclust(hc, k)
# mat<-as.matrix(ct[which(ct==2)])
# mat<-cbind(mat,as.character(unlist(genannot[match(rownames(mat),genannot$V1),2])))
# write.csv(mat,file="HUES3_MxMt_Centralcluster (TE).csv")

png("HUES7_MxMt_Heatmap (TE).png",height=1000,width=1000)
heatmap.2(hues7_mxmt,trace="none")
dev.off()
# hc<-hclust(dist(t(hues7_mxmt)))
# png("HUES7_MxMt_Dendrogram (TE).png",height=1000,width=1000)
# plot(hc,cex=0.5)
# dev.off()
# k<-3
# ct<- cutree(hc, k)
# # rect.hclust(hc, k)
# mat<-as.matrix(ct[which(ct==3)])
# mat<-cbind(mat,as.character(unlist(genannot[match(rownames(mat),genannot$V1),2])))
# write.csv(mat,file="HUES7_MxMt_Centralcluster (TE).csv")