DoIntExp <-
function(statR.o,adj.m,c){

if(length(grep("[a-zA-Z]",rownames(adj.m)))!=0){print("ERROR: The rownames of adj.m should be EntrezID");break}

avexp.m <- statR.o$avexp;
commonEID.v <- intersect(rownames(adj.m),rownames(avexp.m));
mapA.idx <- match(commonEID.v,rownames(adj.m));
tmpA.m <- adj.m[mapA.idx,mapA.idx];

mapR.idx <- match(commonEID.v,rownames(avexp.m));
tmpR.m <- avexp.m[mapR.idx,];

gr.o <- graph.adjacency(tmpA.m,mode="undirected");
comp.l <- clusters(gr.o);
ngpc.v <- summary(factor(comp.l$member));
maxCLid <- as.numeric(names(ngpc.v)[which.max(ngpc.v)]);
maxc.idx <- which(comp.l$member==maxCLid);
#get the max connected network
tmpA.m <- tmpA.m[maxc.idx,maxc.idx];
gr.o <-  graph.adjacency(tmpA.m,mode="undirected");
tmpE.m <- get.edgelist(gr.o);
tmpR.m <- tmpR.m[maxc.idx,];

#### now extract statistics
selcol.idx <- match(c("t","P.Value"),colnames(statR.o$top[[c]]));
ordrow.idx <- match(rownames(tmpR.m),rownames(statR.o$top[[c]]));
if(length(intersect(c("ID"),colnames(statR.o$top[[c]])))==1){
 ordrow.idx <- match(rownames(tmpR.m),statR.o$top[[c]][,1]);
}
statR.m <- statR.o$top[[c]][ordrow.idx,selcol.idx];
rownames(statR.m) <- rownames(tmpR.m);

return(list(statR=statR.m,adj=tmpA.m,avexp=avexp.m));
}
