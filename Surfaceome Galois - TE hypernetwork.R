library(mixOmics)
library(factoextra)
library(gplots)
library(ggplot2)

setwd("B:/Terry/TE_ICM_SEPT2017/TEvICM_FASEB RESPONSE_2019/")
Petr<-read.table("./Petropolous/rpkm.txt")
Petr_samples<-read.table("./Petropolous/E-MTAB-3929.sdrf.txt",sep="\t",header=T)

rn<-rownames(Petr)
Petr<-apply(Petr,2,as.numeric)
rownames(Petr)<-rn
Petr<-Petr[,c(which(Petr_samples$Characteristics.developmental.stage.=="embryonic day 6") ,which(Petr_samples$Characteristics.developmental.stage.=="embryonic day 7"))]
Petr_samples<-Petr_samples[c(which(Petr_samples$Characteristics.developmental.stage.=="embryonic day 6") ,which(Petr_samples$Characteristics.developmental.stage.=="embryonic day 7")),]
Petr<-Petr[,which(Petr_samples$Characteristics.inferred.lineage.=="trophectoderm")]
Petr_samples<-Petr_samples[which(Petr_samples$Characteristics.inferred.lineage.=="trophectoderm"),]

setwd("B://Affiliated people/Sabrina_2019/hypernetwork/")
surfaceome_galois<-read.csv("surfaceome galois full data.csv",header=T,row.names=1,check.names = FALSE)
annot <-read.csv("Surfaceome list for hypernetwork.csv")
annot<-annot[match(colnames(surfaceome_galois),annot$Accession.Number),]
genes<-as.character(annot$Gene.Name)
genes<-genes[!is.na(genes)]
interactions<-read.csv("STRING binding interactions.csv",header=T)
proteins<-read.csv("STRING protein names.csv",header=T)

genes_alt_names<-as.character(proteins$protein_external_id[match(genes,proteins$preferred_name)])
bound_proteins_alt_names<-as.character(interactions$item_id_b[na.omit(match(interactions$ï..item_id_a,genes_alt_names))])
genes_names<-as.character(proteins$preferred_name[match(bound_proteins_alt_names,proteins$protein_external_id)])
genes_names_a<-unique(genes_names)
bound_proteins_alt_names<-as.character(interactions$ï..item_id_a[na.omit(match(interactions$item_id_b,genes_alt_names))])
genes_names<-as.character(proteins$preferred_name[match(bound_proteins_alt_names,proteins$protein_external_id)])
genes_names_b<-c(unique(genes_names),genes_names_a)

cordat<-cor(t(Petr[na.omit(match(genes_names_b,rownames(Petr))),]),t(Petr[-na.omit(match(genes_names_b,rownames(Petr))),]))

# Check distribution of R values using random selection of probes
# rand<-round(rnorm(n=1000,sd = 10000))
# rand[which(rand<0)]<-rand[which(rand<0)]*(-1)
# cordat<-cor(t(Petr[rand,]),t(Petr[-rand,]))
# hist(cordat)

cordat[which(cordat>=0.2 |cordat<(-0.2))]<-1
cordat[which(cordat!=1)]<-0
cordat[which(is.na(cordat))]<-0

mxmt<-cordat%*%t(cordat)

png("Surfaceome Galois vs TE transcriptome - mxmt.png",height=1000,width=1000)
heatmap.2(mxmt,trace="none")
dev.off()

hc<-hclust(dist(t(mxmt)))
plot(hc,cex=0.5)
k<-2
ct<- cutree(hc, k)
rect.hclust(hc, k)
names<-names(ct[which(ct==2)])

galois<-cordat[match(names,rownames(cordat)),]
galois<-galois[,which(colSums(galois)==nrow(galois))]
write.csv(galois,"Surfaceome Galois vs TE transcriptome _ Galois.csv",row.names=T)
