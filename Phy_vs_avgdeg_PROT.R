library(readxl)
library(ggplot2)
phylostratum_vs_degree <- read_excel("B:/Sabrina_2019/phylostratum_vs_degree.xlsx", 
                                     sheet = "ProtAvg")
View(phylostratum_vs_degree)
plot.default(phylostratum_vs_degree)
ggplot(phylostratum_vs_degree, sheet=ProtAvg, aes(x=Phylostratum, y=Degree_Average))+
  geom_point()+
  scale_x_discrete(limits=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19"))+
  geom_smooth()+
  ggtitle("Proteome")

phylostratum_vs_degree <- read_excel("B:/Sabrina_2019/phylostratum_vs_degree.xlsx", 
                                     sheet = "SurAvg")
View(phylostratum_vs_degree)
plot.default(phylostratum_vs_degree)
ggplot(phylostratum_vs_degree, sheet=SurAvges, aes(x=Phylostratum, y=Degree_Average))+
  geom_point()+
  scale_x_discrete(limits=c("1":"19"))+
  geom_smooth()+
  ggtitle("Surfaceome")

phylostratum_vs_degree <- read_excel("B:/Sabrina_2019/phylostratum_vs_degree.xlsx", 
                                     sheet = "GlcAvg")

View(phylostratum_vs_degree)
plot.default(phylostratum_vs_degree)
ggplot(phylostratum_vs_degree, sheet=GlcAvg, aes(x=Phylostratum, y=Degree_Average))+
  geom_point()+
  scale_x_discrete(limits=c("1":"19"))+
  geom_smooth()+
  ggtitle("GlcNAcome")
