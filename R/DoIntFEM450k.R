DoIntFEM450k <-
function(statM.o,statR.o,adj.m,cM,cR){
if(length(grep("[a-zA-Z]",rownames(adj.m)))!=0){print("ERROR: The rownames of adj.m should be EntrezID");break}

require(igraph);
avbeta.m <- statM.o$avbeta;
avexp.m <- statR.o$avexp;
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

#### now extract statistics
selcol.idx <- match(c("t","P.Value"),colnames(statM.o$top[[cM]]));
ordrow.idx <- match(rownames(tmpM.m),rownames(statM.o$top[[cM]]));
idINCL <- FALSE;
if(length(intersect(c("ID"),colnames(statM.o$top[[cM]])))==1){
 idINCL <- TRUE;
 ordrow.idx <- match(rownames(tmpM.m),statM.o$top[[cM]][,1]);
}
statM.m <- statM.o$top[[cM]][ordrow.idx,selcol.idx];
rownames(statM.m) <- rownames(tmpM.m);

selcol.idx <- match(c("t","P.Value"),colnames(statR.o$top[[cR]]));
ordrow.idx <- match(rownames(tmpR.m),rownames(statR.o$top[[cR]]));
idINCL <- FALSE;
if(length(intersect(c("ID"),colnames(statR.o$top[[cR]])))==1){
 idINCL <- TRUE;
 ordrow.idx <- match(rownames(tmpR.m),statM.o$top[[cR]][,1]);
}
statR.m <- statR.o$top[[cR]][ordrow.idx,selcol.idx];
rownames(statR.m) <- rownames(tmpR.m);

return(list(statR=statR.m,statM=statM.m,adj=tmpA.m,avexp=avexp.m,avbeta=avbeta.m));
}
