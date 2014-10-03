DoEpiMod <-
function(statM.m,adj.m,nseeds=100,gamma=0.5,nMC=1000,sizeR.v=c(1,100),minsizeOUT=10,writeOUT=TRUE,nameSTUDY="X",ew.v=NULL){

PasteVector <- function(v){
  vt <- v[1];
  if(length(v) > 1){
   for(g in 2:length(v)){
    vt <- paste(vt,v[g],sep=" ")

   }
  }
  vt <- paste(vt," EnD",sep="");
  out.v <- sub(" EnD","",vt);
  out.v <- sub("NA , ","",out.v);
  out.v <- sub(" , NA","",out.v);
  out.v <- sub(" , NA , "," , ",out.v);
  return(out.v);
}

Heaviside <- function(v){
  out.v <- v;
  out.v[which(v>=0)] <- 1;
  out.v[which(v<0)] <- 0;  
  return(out.v);
}

WriteOutPval <- function(pv.v,round.min=3,round.max=50){
  round.v <- round.min:round.max
  th.v <- 10^(-round.v);
  outpv.v <- vector(length=length(pv.v));
  done.idx <- which(pv.v >= th.v[1]);
  outpv.v[done.idx] <- round(pv.v[done.idx],round.min);
  todo.idx <- setdiff(1:length(pv.v),done.idx);
  for(i in todo.idx){
    if(length(which(th.v <= pv.v[i]))>0){
     outpv.v[i] <- round(pv.v[i],round.v[min(which(th.v <= pv.v[i]))]);
    }
    else{
     outpv.v[i] <- 0;
    }
  }
  return(outpv.v);
}

#require(igraph);
####################################
###3generatte the annotaion matrix
####################################
annotation.m=matrix(nrow=nrow(statM.m),ncol=2)
rownames(annotation.m)=rownames(statM.m)
colnames(annotation.m)=c("EntrezID","GeneSymbol")
annotation.m[,1]=rownames(statM.m)
data(Entrez.GeneSybo.list);
genesymbol=Entrez.GeneSybo.list[annotation.m[,1]]#extract the genesymbols
length(genesymbol);
v=c();
#for(i in 1:length(genesymbol)){v=c(v,genesymbol[[i]])}
for(i in 1:length(genesymbol)){if(length(genesymbol[[i]])==0){gene="None"}else{gene=genesymbol[[i]]};v=c(v,gene)}
annotation.m[,2]=v
anno.m=annotation.m
if(length(statM.m[,1])==84){#determine whether it is toydata.
data(toydata);
anno.m=toydata$annotation;
}
#########################end of generation annoataion matrix
#if(!identical(rownames(statM.m),rownames(statR.m))){
#    print("Please provide gene identifiers as rownames to statM.v and statR.v, and these should be identical");
#}
### ensure that variance of statistics are scaled so they are comparable
statM.v <- statM.m[,1];
#statR.v <- statR.m[,1];
#sdD <- sqrt(var(statM.v));
#sdR <- sqrt(var(statR.v));
#statR.v <- statR.v*sdD/sdR ;

### integrate statistics 
nameSTUDY <- paste("Epi-",nameSTUDY,sep="");

statI.v <- abs(statM.v);
#statI.v[which(sign(statM.v)==sign(statR.v))] <- 0; ### inconsistent ones set to zero
names(statI.v) <- rownames(statM.m);

### find subnetworks around seeds
ntop <- nseeds;

### use greedy approach: rank nodes to define seeds
rankN.s <- sort(statI.v,decreasing=TRUE,index.return=TRUE);
seedsN.v <- names(statI.v)[rankN.s$ix[1:ntop]];

### now define weights on network
print("Constructing weighted network");
tmpA.m <- adj.m;
gr.o <-  graph.adjacency(tmpA.m,mode="undirected");
tmpE.m <- get.edgelist(gr.o);
if(is.null(ew.v)){
tmpW.v <- vector(length=nrow(tmpE.m));
for(e in 1:nrow(tmpE.m)){
  match(tmpE.m[e,],rownames(tmpA.m)) -> map.idx;
  tmpW.v[e] <- mean(statI.v[map.idx]);
  print(e);
}
}
else{
  tmpW.v <- ew.v
}

### a number of edges might have a weight of zero which would later alter the topology of network. this is not desired, hence we replace 0 by the minimum non-zero value.
minval <- min(setdiff(tmpW.v,0))
tmpW.v[tmpW.v==0] <- minval;

### defined weighted graph and adjacency matrix
grW.o <- set.edge.attribute(gr.o,"weight",value=tmpW.v);
adjW.m <- get.adjacency(grW.o,attr="weight")

### Run Spin-Glass algorithm
print("Running Spin-Glass algorithm");
sizeN.v <- vector();
sgcN.lo <- list();
for(v in 1:ntop){
 sgcN.o <- spinglass.community(gr.o,weights=tmpW.v,spins=25,start.temp=1,stop.temp=0.1,cool.fact=0.99,update.rule=c("config"),gamma=gamma,vertex=rankN.s$ix[v]);
 sizeN.v[v] <- length(sgcN.o$comm);
 sgcN.lo[[v]] <- sgcN.o;
 print(paste("Done for seed ",v,sep=""));
}
names(sizeN.v) <- seedsN.v;
print(sizeN.v);
### compute modularities
modN.v <- vector();
for(v in 1:ntop){
 subgr.o <- induced.subgraph(grW.o,sgcN.lo[[v]]$comm);
 modN.v[v] <- mean(get.edge.attribute(subgr.o,name="weight"))
}
names(modN.v) <- seedsN.v;
print(modN.v);

### now determine significance against randomisation of profiles
print("Starting Monte Carlo Runs");
modNmc.m <- matrix(nrow=ntop,ncol=nMC);
for(m in 1:ntop){
  subgr.o <- induced.subgraph(gr.o,sgcN.lo[[m]]$comm);
  nN <- sizeN.v[m];
  if( (nN> sizeR.v[1]) && (nN< sizeR.v[2])){
  tmpEL.m <- get.edgelist(subgr.o);
  for(run in 1:nMC){
   permN.idx <- sample(1:nrow(tmpA.m),nrow(tmpA.m),replace=FALSE);
   tmpEW.v <- vector();
   for(e in 1:nrow(tmpEL.m)){
     match(tmpEL.m[e,],rownames(tmpA.m)[permN.idx]) -> map.idx;
     tmpEW.v[e] <- mean(statI.v[map.idx]);
   }
   subgrW.o <- set.edge.attribute(subgr.o,"weight",value=tmpEW.v)
   modNmc.m[m,run] <- mean(get.edge.attribute(subgrW.o,name="weight"));
  }
  }
  print(paste("Done for seed/module ",m,sep=""));
}


modNpv.v <- rep(1,ntop);
for(v in 1:ntop){
  if( (sizeN.v[v] > sizeR.v[1]) && (sizeN.v[v]< sizeR.v[2])){
    modNpv.v[v] <- length(which(modNmc.m[v,] > modN.v[v]))/nMC;
  }
}
names(modNpv.v) <- seedsN.v;
print(modNpv.v);

### summarize hits
print("Summarising and generating output");
selpvN.idx <- which(modNpv.v < 0.05);
selSize.idx <- which(sizeN.v >= minsizeOUT);
selMod.idx <- intersect(selpvN.idx,selSize.idx);
print(selMod.idx);
print(seedsN.v);
topmodN.m <- matrix(nrow=length(selMod.idx),ncol=6);
match(seedsN.v[selMod.idx],anno.m[,1]) -> map.idx;
seedsSYM.v <- anno.m[map.idx,2];

topmodN.m[,1] <- seedsN.v[selMod.idx];
topmodN.m[,2] <- seedsSYM.v;
topmodN.m[,3:5] <- cbind(sizeN.v[selMod.idx],modN.v[selMod.idx],modNpv.v[selMod.idx]);
mi <- 1;
for(m in selMod.idx){
  tmpEID.v <- rownames(tmpA.m)[sgcN.lo[[m]]$comm];
  genes.v <- anno.m[match(tmpEID.v,anno.m[,1]),2];
  topmodN.m[mi,6] <- PasteVector(genes.v);
  mi <- mi+1;
}
colnames(topmodN.m) <- c("EntrezID(Seed)","Symbol(Seed)","Size","Mod","P","Genes");

if(writeOUT){
write.table(topmodN.m,file=paste("topFEM-",nameSTUDY,".txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE);
}

seltopmodN.lm <- list();
for(m in 1:length(selMod.idx)){
  tmpEID.v <- rownames(tmpA.m)[sgcN.lo[[selMod.idx[m]]]$comm]
  match(tmpEID.v,anno.m[,1]) -> map.idx;
  match(tmpEID.v,rownames(tmpA.m)) -> map1.idx;
  seltopmodN.m <- cbind(anno.m[map.idx,1:2],statM.m[map1.idx,],statI.v[map1.idx]);
  seltopmodN.lm[[m]] <- seltopmodN.m;
  colnames(seltopmodN.lm[[m]]) <- c("EntrezID","Symbol","stat(DNAm)","P(DNAm)","stat(Int)");
}
names(seltopmodN.lm) <- seedsSYM.v

if(writeOUT){

for(m in 1:length(selMod.idx)){
  out.m <- seltopmodN.lm[[m]];
  out.m[,3] <- round(as.numeric(out.m[,3]),2);
  out.m[,4] <- WriteOutPval(as.numeric(out.m[,4]),round.max=100);  
  out.m[,5] <- round(as.numeric(out.m[,5]),2);  
  write(paste(seedsSYM.v[m]," (",nrow(seltopmodN.lm[[m]]), " genes)",sep=""),file=paste("topFEMLists-",nameSTUDY,".txt",sep=""),ncolumns=1,append=TRUE);
  write.table(out.m,file=paste("topFEMLists-",nameSTUDY,".txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE,append=TRUE);
}

}




return(list(size=sizeN.v,mod=modN.v,pvalue=modNpv.v,selmod=selMod.idx,fem=topmodN.m,topmod=seltopmodN.lm,sgc=sgcN.lo,ew=tmpW.v));


}
