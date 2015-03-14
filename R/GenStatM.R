### GenStatM.R

GenStatM <- function(dnaM.m,pheno.v){
dnaM.m=as.matrix(dnaM.m)
#####
##generate average gene's beta matrix:avbeta.m
data("probe450kfemanno");
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

eid450k.v <- as.numeric(probe450kfemanno$eid[match(rownames(dnaM.m),names(probe450kfemanno$eid))]);
map.idx <- match(rownames(dnaM.m),probe450kfemanno[[5]]);
probeInfo.lv <- lapply(probe450kfemanno,extractFn,map.idx);

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



### now derive statistics



data.m <- avbeta.m;

##do lima
### construct model matrix
sampletype.f <- as.factor(pheno.v);
design.sample <- model.matrix(~0 + sampletype.f);
colnames(design.sample) <- levels(sampletype.f);
sampletypes.v <- levels(sampletype.f);

lmf.o <- lmFit(data.m,design.sample);

### construct contrast matrix
ntypes <- length(levels(sampletype.f));
ncomp <- 0.5*ntypes*(ntypes-1);
cont.m <- matrix(0,nrow=ncol(design.sample),ncol=ncomp);
tmp.v <- vector();
c <- 1;
for(i1 in 1:(ntypes-1)){
 for(i2 in (i1+1):ntypes){
   cont.m[i1,c] <- -1;
   cont.m[i2,c] <- 1;
   tmp.v[c] <- paste(sampletypes.v[i2],"--",sampletypes.v[i1],sep="");
   c <- c+1;
 }
}
rownames(cont.m) <- sampletypes.v; # sampletype.v determined separately
colnames(cont.m) <- tmp.v;

### do linear model to contrasts
lmf2.o <- contrasts.fit(lmf.o,cont.m);
### empirical Bayesian estimation of differentially expressed genes (DEGs)
bay.o <- eBayes(lmf2.o);
### build ranked list of DEGs for each comparison
top.lm <- list();
for(c in 1:ncol(cont.m)){
top.lm[[c]] <- topTable(bay.o,coef=c,adjust="fdr",number=nrow(data.m));
}

return(list(top=top.lm,cont=cont.m,avbeta=avbeta.m));

}
