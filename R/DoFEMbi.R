DoFEMbi <- function(intFEM.o,nseeds=10,gamma=1.1,nMC=100,minSize=5,maxSize=50,writeOUT=TRUE,nameSTUDY="X",pvth=0.05,ncores=4){

### auxilliary functions
Heaviside <- function(v) {
        out.v <- rep(0,length(v));
        out.v[which(v >= 0)] <- 1
        return(out.v)
}

PasteVector <- function(v){
    if(length(v)>1){
      out.v <- v[1];
        for(i in 2:length(v)){
            out.v <- paste(out.v,v[i],sep=",");
        }
    }
    else if(length(v)==1){
        out.v <- v;
    }
    return(out.v);
}

### WriteOutPval.R
### function that rounds p-values for nice table generation

WriteOutPval <- function(pv.v,round.min=3,round.max=20){

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
   

EstEdgeWeight <- function(idx,tmpE.m,statI.v){
   map.idx <- match(tmpE.m[idx, ], names(statI.v));
   weight <- mean(statI.v[map.idx])
   return(weight);
}

ConstWeightedNet <- function(adj.m,statI.v,ncores=8){
 gr.o <- graph.adjacency(adj.m, mode = "undirected")
 tmpE.m <- get.edgelist(gr.o)
 idx.l <- as.list(1:nrow(tmpE.m));
 mcl.o <- mclapply(idx.l,EstEdgeWeight,tmpE.m,statI.v,mc.cores=ncores);
 tmpW.v <-  unlist(mcl.o)
 minval <- min(setdiff(tmpW.v, 0))
 tmpW.v[tmpW.v == 0] <- minval
 grW.o <- set.edge.attribute(gr.o, "weight", value = tmpW.v);
 return(list(grW=grW.o,ew=tmpW.v));
}
 
ComputeADH <- function(mod.idx,adjW.m,k.v,gamma){
    nE <- sum(adjW.m)/2;
    Mrs <- sum(apply(adjW.m[mod.idx,-mod.idx],1,sum));
    Ks <- sum(k.v[mod.idx]);
    Kr <- sum(k.v[-mod.idx]);
    EMrs <- Ks*Kr/(2*nE);
    a12 <- Mrs - gamma*EMrs;
    return(a12);
}

FindNN <- function(idx,adjW.m){
    nn.idx <- which(adjW.m[idx,]!=0);
    return(nn.idx);
}

FindNNofMod <- function(mod.idx,adjW.m){
    out <- lapply(as.list(mod.idx),FindNN,adjW.m);
    nn.idx <- setdiff(unlist(out),mod.idx);
    return(nn.idx);
}

AssignNN <- function(mod.idx,nn.idx){
  modNN.l <- lapply(as.list(nn.idx),union,mod.idx);
  return(modNN.l);
}


FindModuleSeed <- function(seed,adjW.m,gamma,maxSize=maxSize){
  k.v <- apply(adjW.m,1,sum); ### generalized degree distribution for weighted network
  nE <- sum(adjW.m)/2; #### generalized mass edge number for weighted network
  ### compute initial adhesion of seed to rest of network
  Mrs <- sum(adjW.m[seed,-seed]);
  Ks <- sum(k.v[seed]);
  Kr <- sum(k.v[-seed]);
  EMrs <- Ks*Kr/(2*nE);
  aS <- Mrs - gamma*EMrs;
  print(paste("Initial Adhesion=",aS,sep=""));
  mod.idx <- seed;
  modSize <- 1;
  adhOld <- aS;
  dADH <- -1;
  while ( (dADH < 0) && (modSize <= maxSize)){
    print(adhOld);
    print(paste("Current Module Size=",modSize,sep=""));
    nn.idx <- FindNNofMod(mod.idx,adjW.m);
    modNN.l <- AssignNN(mod.idx,nn.idx);
    adhNN.v <- unlist(lapply(modNN.l,ComputeADH,adjW.m,k.v,gamma));
    min.idx <- which.min(adhNN.v);
    optNN <- nn.idx[min.idx];
    adhNew <- min(adhNN.v);
    dADH <- adhNew - adhOld;
    if(dADH < 0){
        mod.idx <- c(mod.idx,optNN);
        modSize <- modSize + 1;
    }
    adhOld <- adhNew;
  }
  return(mod.idx);
}

CompModularity <- function(mod.idx,adjW.m){
    tmpAW.m <- as.matrix(adjW.m[mod.idx,mod.idx]);
    tmpA.m <- sign(tmpAW.m);
    return(mean(tmpAW.m[tmpA.m==1]));
}

doMC <- function(mcIDX,selMod.l,adj.m,statI.v){
    permN.idx <- sample(1:nrow(adj.m), nrow(adj.m),replace = FALSE);
    gr.o <- graph.adjacency(adj.m,diag=FALSE);
    nmod <- length(selMod.l);
    modMC.v <- vector(length=nmod);
    for (m in 1:nmod) {
      subgr.o <- induced.subgraph(gr.o,selMod.l[[m]]);
      tmpEL.m <- get.edgelist(subgr.o);
      map.m <- apply(tmpEL.m,1,match,rownames(adj.m[permN.idx,]));
      tmpEW.v <- apply(matrix(statI.v[map.m],nrow=nrow(tmpEL.m),ncol=2,byrow=TRUE),1,mean);
      modMC.v[m] <- mean(tmpEW.v);
    }
    return(modMC.v);
}

######## load data ###########
adj.m <- intFEM.o$adj;
statM.m <- intFEM.o$statM;
statR.m <- intFEM.o$statR;
statM.v <- statM.m[, 1]
statR.v <- statR.m[, 1]
sdD <- sqrt(var(statM.v))
sdR <- sqrt(var(statR.v))
statR.v <- statR.v * sdD/sdR
statI.v <- (Heaviside(statM.v) * Heaviside(-statR.v) + Heaviside(-statM.v) * 
        Heaviside(statR.v)) * abs(statM.v - statR.v)
statI.v[which(sign(statM.v) == sign(statR.v))] <- 0
names(statI.v) <- rownames(statM.m)
rankN.s <- sort(statI.v, decreasing = TRUE, index.return = TRUE)
rankedN.v <- names(statI.v)[rankN.s$ix];
seedsN.v <- rankedN.v[1:nseeds];
    
#### construct weighted network
print("Construct Weighted net");
wnet.o <- ConstWeightedNet(adj.m,statI.v,ncores=ncores);
grW.o <- wnet.o$grW;
adjW.m <- get.adjacency(grW.o, attr = "weight")
print("Done");
seed.idx <- match(seedsN.v,rownames(adjW.m));
seedIDX.l <- as.list(seed.idx);
print("Finding Modules around seeds");
mod.l <- mclapply(seedIDX.l,FindModuleSeed,adjW.m,gamma,maxSize=maxSize,mc.cores=ncores);
print("Done");
modSize.v <- unlist(lapply(mod.l,length));
selMod.idx <- which(modSize.v >= minSize);
selMod.l <- mod.l[selMod.idx];
#### compute modularity of modules and significance
print("Computing modularity");
modularity.v <- unlist(mclapply(selMod.l,CompModularity,adjW.m,mc.cores=ncores));
print("Starting Monte Carlo Runs")
idxMC.l <- as.list(1:nMC);
mcl.o <- mclapply(idxMC.l,doMC,selMod.l,adj.m,statI.v,mc.cores=ncores);
modMC.m <- matrix(unlist(mcl.o),nrow=nMC,ncol=length(selMod.l),byrow=TRUE);
modPV.v <- rep(1,length(selMod.l));
for (v in 1:length(modPV.v)) {
  modPV.v[v] <- length(which(modMC.m[,v] > modularity.v[v]))/nMC;
}
print("Summarising and generating output")
x <- org.Hs.egSYMBOL
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])
mapEIDtoSYM.v <- unlist(xx)
map.idx <- match(rownames(adj.m), names(mapEIDtoSYM.v))
anno.m <- cbind(rownames(adj.m), mapEIDtoSYM.v[map.idx])
colnames(anno.m) <- c("EntrezID", "Symbol")

selPV.idx <- which(modPV.v < pvth)
if(length(selPV.idx)==0){
  print("No significant modules found. Try relaxing significance threshold, change gamma or increase seed number!");
}
else if (length(selPV.idx)>=1){
  topmodN.m <- matrix(nrow = length(selPV.idx), ncol = 6)
  map.idx <- match(seedsN.v[selMod.idx[selPV.idx]], anno.m[, 1])
  seedsSYM.v <- anno.m[map.idx, 2]
  topmodN.m[, 1] <- seedsN.v[selMod.idx[selPV.idx]];
  topmodN.m[, 2] <- seedsSYM.v
  topmodN.m[, 3:5] <- cbind(modSize.v[selMod.idx[selPV.idx]], modularity.v[selPV.idx], modPV.v[selPV.idx]);
  mi <- 1
  for (m in selPV.idx){
        tmpEID.v <- rownames(adj.m)[selMod.l[[m]]];
        genes.v <- anno.m[match(tmpEID.v, anno.m[, 1]), 2]
        topmodN.m[mi, 6] <- PasteVector(genes.v)
        mi <- mi + 1
   }
   colnames(topmodN.m) <- c("EntrezID(Seed)", "Symbol(Seed)", 
        "Size", "Mod", "P", "Genes")
    if (writeOUT) {
        write.table(topmodN.m, file = paste("topFEM-", nameSTUDY, 
            ".txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
    }
    seltopmodN.lm <- list()
    for (m in 1:length(selPV.idx)) {
        tmpEID.v <- rownames(adj.m)[selMod.l[[selPV.idx[m]]]];
        map.idx <- match(tmpEID.v, anno.m[, 1])
        map1.idx <- match(tmpEID.v, rownames(adj.m))
        seltopmodN.m <- cbind(anno.m[map.idx, 1:2], statM.m[map1.idx, ],statR.m[map1.idx, ], statI.v[map1.idx])
        seltopmodN.lm[[m]] <- seltopmodN.m
        colnames(seltopmodN.lm[[m]]) <- c("EntrezID", "Symbol", 
            "stat(DNAm)", "P(DNAm)", "stat(mRNA)", "P(mRNA)", 
            "stat(Int)")
    }
    names(seltopmodN.lm) <- seedsSYM.v
    if (writeOUT) {
        for (m in 1:length(selPV.idx)) {
            out.m <- seltopmodN.lm[[m]]
            out.m[, 3] <- round(as.numeric(out.m[, 3]), 2)
            out.m[, 4] <- WriteOutPval(as.numeric(out.m[, 4]), 
                round.max = 100)
            out.m[, 5] <- round(as.numeric(out.m[, 5]), 2)
            out.m[, 6] <- WriteOutPval(as.numeric(out.m[, 6]), 
                round.max = 100)
            out.m[, 7] <- round(as.numeric(out.m[, 7]), 2)
            write(paste(seedsSYM.v[m], " (", nrow(seltopmodN.lm[[m]]), 
                " genes)", sep = ""), file = paste("topFEMLists-", 
                nameSTUDY, ".txt", sep = ""), ncolumns = 1, append = TRUE)
            write.table(out.m, file = paste("topFEMLists-", nameSTUDY, 
                ".txt", sep = ""), quote = FALSE, sep = "\t", 
                row.names = FALSE, append = TRUE)
        }
    }
   }
    return(list(size = modSize.v[selMod.idx[selPV.idx]], mod = modularity.v[selPV.idx], pv = modPV.v[selPV.idx],sigmod = selMod.idx[selPV.idx], fem = topmodN.m, topmod = seltopmodN.lm, sgc=selMod.l, ew=wnet.o$ew , adj = intFEM.o$adj));


} ### EOF


