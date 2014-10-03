DoIntFEM450k <-
function(dnaM.m,exp.m,phenoM.v,phenoR.v,adj.m){

data("map450kEID.v");
data("probeInfoALL.lv");
extractFn <- function(tmp.v,ext.idx){
    return(tmp.v[ext.idx]);
}
avVEC <- function(data.m,idx.v){
 av.v <- apply(matrix(data.m[idx.v,],nrow=length(idx.v),ncol=ncol(data.m)),2,mean,na.rm=T);
 return(av.v);
}


MeanRep <- function(tmp.m){

 tmp.v <- union(rownames(tmp.m),rownames(tmp.m));
 ng <- length(tmp.v);
 newtmp.m <- matrix(nrow=ng,ncol=ncol(tmp.m));
 for( g in 1:ng){
  which(rownames(tmp.m)==tmp.v[g]) -> tmp.idx;
  newtmp.m[g,] <- avVEC(tmp.m,tmp.idx);
  #print(paste(g," done out of ",ng,sep=""));
 }
 colnames(newtmp.m) <- colnames(tmp.m);
 rownames(newtmp.m) <- tmp.v;

 return(newtmp.m);
}

require(igraph);
require(limma);

eid450k.v <- as.numeric(map450kEID.v[match(rownames(dnaM.m),names(map450kEID.v))]);
map.idx <- match(rownames(dnaM.m),probeInfoALL.lv[[5]]);
probeInfo.lv <- lapply(probeInfoALL.lv,extractFn,map.idx);

### loop through gene groups
beta.lm <- list();
for(g in 1:6){
group.idx <- which(probeInfo.lv[[3]]==g);
tmp.m <- dnaM.m[group.idx,];
rownames(tmp.m) <- eid450k.v[group.idx];
sel.idx <- which(is.na(rownames(tmp.m))==FALSE);
beta.lm[[g]] <- MeanRep(tmp.m[sel.idx,]);
print(paste("Done for regional gene group ",g,sep=""));
}

unqEID.v <- unique(c(rownames(beta.lm[[2]]),rownames(beta.lm[[4]]),rownames(beta.lm[[1]])));
### construct average DNA methylation matrix
avbeta.m <- matrix(nrow=length(unqEID.v),ncol=ncol(dnaM.m));
colnames(avbeta.m) <- colnames(dnaM.m);
rownames(avbeta.m) <- unqEID.v;	
# 1:TSS1500 4:1exon 2:TSS200, so 1:TSS200 will overwrite the 1exon and tss1500. 
for(gr in c(1,4,2)){
 avbeta.m[match(rownames(beta.lm[[gr]]),rownames(avbeta.m)),] <- beta.lm[[gr]];#
}



#### average gene expression if needed
avexp.m <- exp.m;
if(length(unique(rownames(exp.m)))!=length(rownames(exp.m))){
 avexp.m <- MeanRep(exp.m);
}


commonEID.v <- intersect(intersect(rownames(adj.m),rownames(avbeta.m)),rownames(avexp.m));
mapA.idx <- match(commonEID.v,rownames(adj.m));
tmpA.m <- adj.m[mapA.idx,mapA.idx];

mapM.idx <- match(commonEID.v,rownames(avbeta.m));
tmpM.m <- avbeta.m[mapM.idx,];
mapR.idx <- match(commonEID.v,rownames(avexp.m));
tmpR.m <- avexp.m[mapR.idx,];

gr.o <- graph.adjacency(tmpA.m,mode="undirected");
comp.l <- clusters(gr.o);
ngpc.v <- summary(factor(comp.l$member));
maxCLid <- as.numeric(names(ngpc.v)[which.max(ngpc.v)]);
maxc.idx <- which(comp.l$member==maxCLid);

tmpA.m <- tmpA.m[maxc.idx,maxc.idx];
gr.o <-  graph.adjacency(tmpA.m,mode="undirected");
tmpE.m <- get.edgelist(gr.o);
tmpM.m <- tmpM.m[maxc.idx,];
tmpR.m <- tmpR.m[maxc.idx,];

#### now generate statistics in Limma
out.o <- DoLimma(tmpM.m,phenoM.v);
selcol.idx <- match(c("t","P.Value"),colnames(out.o$top[[1]]));
ordrow.idx <- match(rownames(tmpM.m),rownames(out.o$top[[1]]));
idINCL <- FALSE;
if(length(intersect(c("ID"),colnames(out.o$top[[1]])))==1){
 idINCL <- TRUE;
 ordrow.idx <- match(rownames(tmpM.m),out.o$top[[1]][,1]);
}
statM.m <- out.o$top[[1]][ordrow.idx,selcol.idx];
rownames(statM.m) <- rownames(tmpM.m);



out.o <- DoLimma(tmpR.m,phenoR.v);
selcol.idx <- match(c("t","P.Value"),colnames(out.o$top[[1]]));
ordrow.idx <- match(rownames(tmpR.m),rownames(out.o$top[[1]]));
if(length(intersect(c("ID"),colnames(out.o$top[[1]])))==1){
 ordrow.idx <- match(rownames(tmpR.m),out.o$top[[1]][,1]);
}
statR.m <- out.o$top[[1]][ordrow.idx,selcol.idx];
rownames(statR.m) <- rownames(tmpR.m);

######################################################################################
#annotation.m
######################################################################################

annotation.m=matrix(nrow=nrow(tmpM.m),ncol=2)
rownames(annotation.m)=rownames(tmpM.m)
colnames(annotation.m)=c("EntrezID","GeneSymbol")
annotation.m[,1]=rownames(tmpM.m)
#library(org.Hs.eg.db);
#use  org.Hs.egSYMBOL it is ok when run the source("DoInt450k.R"),however when in hte pacakge it failed with errori: Error in as.list.default(x[mapped_genes]) :  no method for coercing this S4 class to a vector:
#x <- org.Hs.egSYMBOL
# Get the gene symbol that are mapped to an entrez gene identifiers
#mapped_genes <- mappedkeys(x)
# Convert to a list
#xx <- as.list(x[mapped_genes])
data(Entrez.GeneSybo.list);
genesymbol=Entrez.GeneSybo.list[annotation.m[,1]]#extract the genesymbols
length(genesymbol);
v=c();
for(i in 1:length(genesymbol)){v=c(v,genesymbol[[i]])}
annotation.m[,2]=v
###################################################################
return(list(statR=statR.m,statM=statM.m,adj=tmpA.m,avexp=avexp.m,avbeta=avbeta.m,annotation=annotation.m));
}
