library(mixOmics)
library(factoextra)
library(gplots)
library(ggplot2)

setwd("R:/Terry/TE_ICM_SEPT2017/TEvICM_FASEB RESPONSE_2019/")
Petr<-read.table("./Petropolous/rpkm.txt")
Petr_samples<-read.table("./Petropolous/E-MTAB-3929.sdrf.txt",sep="\t",header=T)

rn<-rownames(Petr)
Petr<-apply(Petr,2,as.numeric)
rownames(Petr)<-rn
Petr<-Petr[,c(which(Petr_samples$Characteristics.developmental.stage.=="embryonic day 6") ,which(Petr_samples$Characteristics.developmental.stage.=="embryonic day 7"))]
Petr_samples<-Petr_samples[c(which(Petr_samples$Characteristics.developmental.stage.=="embryonic day 6") ,which(Petr_samples$Characteristics.developmental.stage.=="embryonic day 7")),]
Petr<-Petr[,which(Petr_samples$Characteristics.inferred.lineage.=="epiblast")]
Petr_samples<-Petr_samples[which(Petr_samples$Characteristics.inferred.lineage.=="epiblast"),]

setwd("R://Affiliated people/Sabrina_2019/hypernetwork/")
surfaceome_ligand_galois<-read.csv("Surfaceome Galois vs TE transcriptome _ Galois.csv",header=T,row.names=1,check.names = FALSE)
annot <-read.csv("Surfaceome list for hypernetwork.csv")
annot<-annot[na.omit(match(colnames(surfaceome_ligand_galois),annot$Gene.Name)),]
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

cordat[which(cordat>=0.2 |cordat<(-0.2))]<-1
cordat[which(cordat!=1)]<-0
cordat[which(is.na(cordat))]<-0

mxmt<-cordat%*%t(cordat)

png("TE vs ICM transcriptome - mxmt.png",height=1000,width=1000)
heatmap.2(mxmt,trace="none")
dev.off()