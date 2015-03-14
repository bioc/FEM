### GenStatR.R

GenStatR <- function(exp.m,pheno.v){
if(length(grep("[a-zA-Z]",rownames(exp.m)))!=0){print("ERROR: The rownames of exp.m should be EntrezID");break}

### Library needed
require(limma);
### average probes mapping to the same gene
#nspg.v <- summary(factor(rownames(exp.m)),maxsum = 10000000);## maxsum : integer, indicating how many levels should be shown for ‘factor’s
#avexp.m <- rowsum(exp.m,group=rownames(exp.m))/nspg.v;#the m=rowsum(tmp.m,group=rownames(tmp.m)), the rownames(m) is same to names(nspg.v)

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
avexp.m=MeanRep(exp.m)

### construct model matrix
sampletype.f <- as.factor(pheno.v);
design.sample <- model.matrix(~0 + sampletype.f);
colnames(design.sample) <- levels(sampletype.f);
sampletypes.v <- levels(sampletype.f);


### do linear model fit
data.m <- avexp.m;
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

return(list(top=top.lm,cont=cont.m,avexp=avexp.m));

}
