DoIntEpi450k <-
function(statM.o,adj.m,c){

if(length(grep("[a-zA-Z]",rownames(adj.m)))!=0){print("ERROR: The rownames of adj.m should be EntrezID");break}

avbeta.m <- statM.o$avbeta;
commonEID.v <- intersect(rownames(adj.m),rownames(avbeta.m));
mapA.idx <- match(commonEID.v,rownames(adj.m));
tmpA.m <- adj.m[mapA.idx,mapA.idx];

mapM.idx <- match(commonEID.v,rownames(avbeta.m));
tmpM.m <- avbeta.m[mapM.idx,];

gr.o <- graph.adjacency(tmpA.m,mode="undirected");
comp.l <- clusters(gr.o);
ngpc.v <- summary(factor(comp.l$member));
maxCLid <- as.numeric(names(ngpc.v)[which.max(ngpc.v)]);
maxc.idx <- which(comp.l$member==maxCLid);
#get the max connected network
tmpA.m <- tmpA.m[maxc.idx,maxc.idx];
gr.o <-  graph.adjacency(tmpA.m,mode="undirected");
tmpE.m <- get.edgelist(gr.o);
tmpM.m <- tmpM.m[maxc.idx,];

#### now extract statistics
selcol.idx <- match(c("t","P.Value"),colnames(statM.o$top[[c]]));
ordrow.idx <- match(rownames(tmpM.m),rownames(statM.o$top[[c]]));
if(length(intersect(c("ID"),colnames(statM.o$top[[c]])))==1){
 ordrow.idx <- match(rownames(tmpM.m),statM.o$top[[c]][,1]);
}
statM.m <- statM.o$top[[c]][ordrow.idx,selcol.idx];
rownames(statM.m) <- rownames(tmpM.m);

return(list(statM=statM.m,adj=tmpA.m,avbeta=avbeta.m));
}
