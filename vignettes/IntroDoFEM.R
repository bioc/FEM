### R code from vignette source 'IntroDoFEM.Rnw'

###################################################
### code chunk number 1: IntroDoFEM.Rnw:44-47 (eval = FALSE)
###################################################
## if (!requireNamespace("BiocManager", quietly=TRUE))
##     install.packages("BiocManager")
## BiocManager::install("marray")


###################################################
### code chunk number 2: IntroDoFEM.Rnw:51-55
###################################################
library("igraph");
library("marray");
library("corrplot");
library("graph");


###################################################
### code chunk number 3: IntroDoFEM.Rnw:60-61
###################################################
library("FEM");	


###################################################
### code chunk number 4: IntroDoFEM.Rnw:64-65
###################################################
data(Toydata)


###################################################
### code chunk number 5: IntroDoFEM.Rnw:69-70
###################################################
names(Toydata)


###################################################
### code chunk number 6: IntroDoFEM.Rnw:74-75
###################################################
head(Toydata$statM)


###################################################
### code chunk number 7: IntroDoFEM.Rnw:81-83
###################################################
rownames(Toydata$statM)[which(Toydata$statM[,1]>2)]->tennodes;
tennodes;


###################################################
### code chunk number 8: IntroDoFEM.Rnw:86-87
###################################################
head(Toydata$statR)


###################################################
### code chunk number 9: IntroDoFEM.Rnw:90-91
###################################################
rownames(Toydata$statM)[which(Toydata$statR[,1]< -2)];


###################################################
### code chunk number 10: IntroDoFEM.Rnw:98-101
###################################################
mod.idx <- match(Toydata$tennodes,rownames(Toydata$adj));
plot.igraph(graph.adjacency(Toydata$adjacency,mod="undirected"),
vertex.size=8,mark.groups=mod.idx,mark.col="yellow")


###################################################
### code chunk number 11: IntroDoFEM.Rnw:108-109
###################################################
intFEM.o <- list(statM=Toydata$statM,statR=Toydata$statR,adj=Toydata$adj);


###################################################
### code chunk number 12: IntroDoFEM.Rnw:111-113
###################################################
DoFEMtoy.o  <- DoFEMbi(intFEM.o,nseeds=1,gamma=1.1,nMC=1000,
minSize=5,maxSize=50,writeOUT=TRUE,nameSTUDY="TOY");


###################################################
### code chunk number 13: IntroDoFEM.Rnw:131-132
###################################################
DoFEMtoy.o$fem


###################################################
### code chunk number 14: IntroDoFEM.Rnw:137-138
###################################################
head(DoFEMtoy.o$topmod$SGMS2,n=5L)


###################################################
### code chunk number 15: IntroDoFEM.Rnw:142-145
###################################################
mod.idx<- match(DoFEMtoy.o$topmod$SGMS2[,1],rownames(Toydata$adj));
plot.igraph(graph.adjacency(Toydata$adjacency,mod="undirected"),
vertex.size=8,mark.groups=mod.idx,mark.col="yellow")


###################################################
### code chunk number 16: IntroDoFEM.Rnw:150-152
###################################################
sensitivity=length(intersect(tennodes,rownames(Toydata$adj)[mod.idx]))/length(tennodes);
sensitivity


###################################################
### code chunk number 17: IntroDoFEM.Rnw:161-163
###################################################
data(Realdata);
attributes(Realdata);


###################################################
### code chunk number 18: IntroDoFEM.Rnw:167-168
###################################################
intFEM.o <- list(statM=Realdata$statM,statR=Realdata$statR,adj=Realdata$adj)


###################################################
### code chunk number 19: IntroDoFEM.Rnw:171-174
###################################################
#Realdata$fembi.o <- DoFEMbi(intFEM.o,
#                            nseeds=100,gamma=0.5,nMC=1000,
#                            minSize=5,maxSize=50,writeOUT=TRUE,nameSTUDY="TCGA-EC");


###################################################
### code chunk number 20: IntroDoFEM.Rnw:179-180
###################################################
Realdata$fembi.o$fem


###################################################
### code chunk number 21: IntroDoFEM.Rnw:184-185
###################################################
Realdata$fembi.o$topmod


###################################################
### code chunk number 22: IntroDoFEM.Rnw:189-193
###################################################
library("marray");
library("corrplot");
HAND2.mod<-Realdata$fembi.o$topmod$HAND2;
HAND2.graphNEL.o=FemModShow(Realdata$fembi.o$topmod$HAND2,name="HAND2",Realdata$fembi.o)


###################################################
### code chunk number 23: IntroDoFEM.Rnw:209-212
###################################################
#for(m in 1:length(names(Realdata$fembi.o$topmod))){
#FemModShow(Realdata$fembi.o$topmod[[m]],
#name=names(Realdata$fembi.o$topmod)[m],Realdata$fembi.o)}


###################################################
### code chunk number 24: IntroDoFEM.Rnw:227-228
###################################################
#DoIntFEM450k(statM.o,statR.o,adj.m,cM,cR)


###################################################
### code chunk number 25: IntroDoFEM.Rnw:247-253
###################################################
#statM.o <- GenStatM(dnaM.m,phenoM.v,"450k");
#intEpi.o=DoIntEpi450k(statM.o,adj.m,c=1)
#EpiMod.o=DoEpiMod(intEpi.o,
#                          nseeds=100,gamma=0.5,nMC=1000,sizeR.
#                          v=c(1,100), minsizeOUT=10,writeOUT=TRUE,
#                          nameSTUDY="TCGA-EC",ew.v=NULL);


###################################################
### code chunk number 26: IntroDoFEM.Rnw:262-269
###################################################
#statR.o <- GenStatR(exp.m,pheno.v);
#intExp.o=DoIntExp(statR.o,adj.m)
#ExpMod.o=DoExpMod(intExp.o,
#                  nseeds=100,gamma=0.5,nMC=1000,
#                  sizeR.v=c(1,100),minsizeOUT=10,
#                  writeOUT=TRUE,nameSTUDY="TCGA-EC",ew.v=NULL)
#


###################################################
### code chunk number 27: IntroDoFEM.Rnw:277-290
###################################################
#library(minfi);
#require(IlluminaHumanMethylation450kmanifest);

#baseDIR <- getwd();# the base dir of the Rawdata
#setwd(baseDIR);
#targets <- read.450k.sheet(baseDIR);#read the csv file.
#RGset <- read.450k.exp(baseDIR); #Reads an entire 450k experiment 
#                                  using a sample sheet
#MSet.raw <- preprocessRaw(RGset);#Converts the Red/Green channel for an Illumina 
#                                  methylation array into methylation signal,
#                                  without using any normalization
#beta.m <- getBeta(MSet.raw,type = "Illumina");# get normalized beta 
#pval.m <- detectionP(RGset,type ="m+u")


