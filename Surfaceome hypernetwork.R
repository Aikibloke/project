library("gplots")
setwd("B://Affiliated people/Sabrina_2019/stiparo/")
genes <-read.csv("Stiparo TE gene direct partners.csv")
setwd("B://Affiliated people/Sabrina_2019/hypernetwork/")
data <-read.csv("Surfaceome list for hypernetwork.csv")

TE_genes<-as.character(unlist(genes))

hesc2<-data[1:nrow(data),3:ncol(data)]
rownames(hesc2)<-as.character(unlist(data[1:nrow(data),2]))

genannot<-data[1:nrow(data),1:2]

rn<-rownames(hesc2)
hesc2<-apply(hesc2,2,as.numeric)
rownames(hesc2)<-rn

surfaceome<-hesc2[1:nrow(hesc2),28:ncol(hesc2)]

surfaceome<-cor(t(surfaceome[na.omit(match(TE_genes,genannot$Gene.Name)),]),t(surfaceome[-na.omit(match(TE_genes,genannot$Gene.Name)),]))

Test<-t(surfaceome[-na.omit(match(TE_genes,genannot$Gene.Name))])
View(Test)[1, 1:10]
hist(surfaceome)

surfaceome[which(surfaceome>=0.8 |surfaceome<(-0.8))]<-1
surfaceome[which(surfaceome!=1)]<-0
surfaceome[which(is.na(surfaceome))]<-0

surfaceome_mxmt<-(surfaceome %*% t(surfaceome))-1

png("surfaceome_mxmt2.png",height=1000,width=1000)
heatmap.2(surfaceome_mxmt,trace="none")
dev.off()


# cluster_6<-read_excel("R:/Affiliated people/Sabrina_2019/Hypernetwork/surfaceome_mxmt_clustergenes2.xlsx",sheet = "Cluster 6 Surfaceome ")
# cluster_6<-as.character(unlist(cluster_6[,2]))
# 
# galois<-surfaceome[match(cluster_6,rownames(surfaceome)),]
# galois<-galois[,which(colSums(galois)==31)]
# 
# write.csv(galois,"surfaceome galois full data.csv",row.names = TRUE)



hc<-hclust(dist(t(surfaceome_mxmt)))
plot(hc,cex=0.5)
# png("surfaceome_mxmt_Dendrogram.png",height=1000,width=1000)
# plot(hc,cex=0.5)
# dev.off()
k<-7
ct<- cutree(hc, k)
rect.hclust(hc, k)
mat<-as.matrix(ct[which(ct=1)])
mat<-cbind(mat,as.character(unlist(genannot[match(rownames(mat),genannot$Accession.Number),1])))
# write.csv(mat,file="surfaceome_mxmt_clustergenes.csv")
writeClipboard(mat)
write.table(mat, file = "surfaceome_mxmt_clustergenes.txt",sep = "\t",row.names = FALSE)


hc<-hclust(dist(t(surfaceome_mxmt)))
plot(hc,cex=0.5)
# png("surfaceome_mxmt_Dendrogram.png",height=1000,width=1000)
# plot(hc,cex=0.5)
# dev.off()
k<-11
ct<- cutree(hc, k)
rect.hclust(hc, k)
mat<-as.matrix(ct)
mat<-cbind(mat,as.character(unlist(genannot[match(rownames(mat),genannot$Accession.Number),1])))
# write.csv(mat,file="surfaceome_mxmt_clustergenes.csv")
writeClipboard(mat)
write.table(mat, file = "surfaceome_mxmt_clustergenes2.txt",sep = "\t",row.names = FALSE)

#Organising String protein names

setwd("B://Affiliated people/Sabrina_2019/Hypernetwork/")
names <-read.csv("STRING protein names.csv")
interactions <-read.csv("STRING binding interactions.csv")

interact2 <-interactions[1:nrow(interactions),1:2]
colnames(interact2)<-c("ID A","ID B")

merge(names, interact2, by.x ='protein_external_id', by.y = 'ID A')
merge(names, interact2, by.x = 'protein_external_id', by.y = 'ID B')

output<-merge(interact2, names, by.x = 'ID A', by.y = 'protein_external_id')
colnames(output)[3]<- 'geneA'

output2<-merge(output, names, by.x = 'ID B', by.y= 'protein_external_id')
colnames(output2)[4]<- 'geneB'

library(xlsx)
write.xlsx(output2, file= "binding interactions with names.xlsx")
library(csv)
write.csv(output2,file="binding_interactions_genes.csv", sep = "/B")



