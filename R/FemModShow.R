FemModShow <-
function(mod,name,fem.o,mode= "Integration"){
edgeweight=fem.o$ew
adjacency=fem.o$adj
###############add shape to the vertex of igraph

mycircle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }

  mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
         vertex.size, vertex.frame.width,
         FUN=function(x, y, bg, fg, size, lwd) {
           symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                   circles=size, add=TRUE, inches=FALSE)
         })
}

#library("igraph");
#library("marray");
#library("corrplot");

#mod is from fem result such as HAND2, hand2<-fembi.o$topmod$HAND2. 
#overall ideas: give the mtval, rtval, vmcolor, vrcolor, label to vertex; give weight, edgewidth to the edge of realgraph. then subgraph the HAND2 mod

realgraph=graph.adjacency(adjacency,mode="undirected")

#add the values
E(realgraph)$weight= edgeweight#give the weight to the edges

#V(realgraph)$mtval=realdata$statM[,1] #give the tvalue of Methylation on the vertex
#V(realgraph)$rtval=realdata$statR[,1] #give the tvalue of Methylation on the vertex

#there are vertex  1 edge width, 2 methylation colar, 3 rna expression, 4 gene symbol vertex label, 4 label size to provide to the igraph plot the other arguments need to set  are 5 label cex (size)
#vertex.v=as.vector(V(realgraph))
#edgevectro=as.vector(E(realgraph)))
#head(fembi.o$ew)
#################################################################################
# 1 edge size
E(realgraph)$edgewidth<-ifelse(E(realgraph)$weight<=1,0.5,E(realgraph)$weight/2)


#################################################################################
# 2 and 3  transform methylation value to node color, rna expression
mod=as.data.frame(mod);
if(mode=="Epi"){
mod[,"stat(mRNA)"]=mod[,"stat(DNAm)"];
}else if(mode=="Exp"){
mod[,"stat(DNAm)"]=mod[,"stat(mRNA)"];
}
# if  there is mode is epi we add the stat(mRNA) colum and set them as 0;
# sub graph mod and inhebited the edgewidth add the mval and rval
mod.graph=igraph::induced.subgraph(realgraph,v=as.vector(mod[,1]))
print(mod[,1])
#the fembi.o$topmod$HAND2 is a dataframe and it has "Symbol",  "stat(DNAm)" "stat(mRNA)" which is useful
mtval.v=vector()
for(i in V(mod.graph)$name){mtval.v=c(mtval.v,(as.vector(mod[i,"stat(DNAm)"])))} #add the mtval to the mod graph one by one according to the squence of mod graph name
V(mod.graph)$mtval=mtval.v;


rtval.v=vector()
for(i in V(mod.graph)$name){rtval.v=c(rtval.v,(as.vector(mod[i,"stat(mRNA)"])))} #add the 
V(mod.graph)$rtval=rtval.v;
print(rtval.v)
#print(rtval.v)
#V(realgraph)$mtval=realdata$statM[,1] #give the tvalue of Methylation on the vertex
#V(realgraph)$rtval=realdata$statR[,1] #give the tvalue of Methylation on the vertex
# add the vm.color, vr.color
vm.color=rep(0,length(V(mod.graph)));
vr.color=rep(0,length(V(mod.graph)));
# color scheme generate 100 colors
tmcolor.scheme<-maPalette(low = "yellow", high ="blue", mid="grey", k =100);
trcolor.scheme<-maPalette(low = "green", high ="red", mid="grey", k =100);

tmcolor.scheme[13:88]="#BEBEBE";#( floor(-1.5/0.04)+51,floor(1.5/0.04)+51) which is [13,88] is grey "#BEBEBE". 1.5 is the thresh hold for the t values 
tmcolor.scheme[1:12]=maPalette(low = "yellow", high="lightyellow", k =12);
tmcolor.scheme[89:100]=maPalette(low = "lightblue", high="blue", k =12);
trcolor.scheme[13:88]="#BEBEBE";#[ floor(-1.5/0.04)+51,floor(1.5/0.04)+51] which is [13,88] is grey "#BEBEBE". 1.5 is the thresh hold for the t values 
trcolor.scheme[1:12]=maPalette(low = "green", high="lightgreen", k =12);
trcolor.scheme[89:100]=maPalette(low = "#DC6868", high="red", k =12);

#give the color according the mtval. (-2,2),floor get the integer mod on 0.04 + 51 then get the according color.
tmcolor.position=floor(as.numeric(V(mod.graph)$mtval)/0.04)+51;
tmcolor.position[which(tmcolor.position<1)]<-1;
tmcolor.position[which(tmcolor.position>100)]<-100;
vm.color=tmcolor.scheme[tmcolor.position];
V(mod.graph)$vmcolor<-vm.color#add he vm.color to the vertex value
print(vm.color)

#add the frame color idea: get the tr color position then get the color from trcolor.scheme
trcolor.position=floor(as.numeric(V(mod.graph)$rtval)/0.04)+51;
print(trcolor.position)
trcolor.position[which(trcolor.position<1)]<-1;
trcolor.position[which(trcolor.position>100)]<-100;
vr.color=trcolor.scheme[trcolor.position];
#vr.color[which(V(mod.graph)$rtval>=2)]<-trcolor.scheme[100];
#vr.color[which(V(mod.graph)$rtval<=(-2))]<-trcolor.scheme[1];
#vr.color[which(V(mod.graph)$rtval>-1 & V(realgraph)$rtval<1)]<-"lightgrey"
#vr.color[which(V(mod.graph)$rtval>=1)]<-"red"
#vr.color[which(V(mod.graph)$rtval<(-1))]<-"green"
V(mod.graph)$vrcolor<-vr.color# the rna expression color

if(mode=="Exp"){
	V(mod.graph)$color<-V(mod.graph)$vrcolor;
}else{
	V(mod.graph)$color<-V(mod.graph)$vmcolor #use the vmcolor as the vertex color but if the mod is Exp them the vertex color is vrcolor
}	
print(vr.color)
#################################################################################
#add the mod label value

#!!!!!!!!!!!!!!!V(mod.graph)$label<-realdata$annotation[,2]


#subnet the hand 2
#V(mod.graph)$labels=label.v
label.v=vector()
for(i in V(mod.graph)$name){label.v=c(label.v,(as.vector(mod[i,"Symbol"])))} #add the V(mod.graph)$name's labels one by one from mod["$name","Symbol"]


#################################################################################
#create subgraph label.cex value
V(mod.graph)$label.cex=rep(0.5,length(as.vector(V(mod.graph)))); #all the cex first set as 0.7
V(mod.graph)$label.cex[which(as.vector(V(mod.graph)$name)==as.vector(mod[1,1]))]=0.8 #only the firt mod name was set as 1


#################################################################################
# generate the plot 
# when you want to plot the vertex shape, and its frame width first you should load the api script api bellow
add.vertex.shape("fcircle", clip=igraph.shape.noclip,plot=mycircle, parameters=list(vertex.frame.color=1,vertex.frame.width=1))


pdf(paste(name,".mod.pdf",sep=""))
if(mode =="Integration"){
#plot(mod.graph,vertex.shape="fcircle", vertex.frame.color=V(mod.graph)$vrcolor,vertex.frame.width=6,vertex.label=E(mod.graph)$label,vertex.label.dist=1,vertex.label.cex=0.7,vertex.label.font=3,edge.color="grey",edge.width=E(mod.graph)$edgewidth,main=paste(name,"-centred module",sep=""))
#plot(HAND2.graph,vertex.shape="fcircle", vertex.frame.color=V(HAND2.graph)$vrcolor,vertex.frame.width=6,vertex.label=E(HAND2.graph)$label,vertex.label.dist=1,vertex.label.cex=label.cex,vertex.label.font=3,edge.color="grey",edge.width=E(HAND2.graph)$edgewidth)
	plot(mod.graph,layout=layout.fruchterman.reingold,vertex.shape="fcircle", vertex.frame.color=V(mod.graph)$vrcolor,vertex.frame.width=4,vertex.size=10,vertex.label=label.v,vertex.label.dist=0.6,vertex.label.cex=V(mod.graph)$label.cex,vertex.label.font=3,edge.color="grey",edge.width=E(mod.graph)$edgewidth)
#plot(mod.graph,vertex.shape="fcircle", vertex.frame.color=V(mod.graph)$vrcolor,vertex.frame.width=6,vertex.label=paste(sep="\n",V(mod.graph)$label,V(mod.graph)$mtval,V(mod.graph)$rtval),vertex.label.dist=0.8,vertex.label.cex=label.cex,vertex.label.font=3,edge.color="grey",edge.width=E(mod.graph)$edgewidth)

	colorlegend(trcolor.scheme,seq(-2,2,0.5),ratio.colbar=0.3 ,xlim=c(-1.55,-1.4),ylim=c(-0.5,0),align="r",cex=0.5)
	colorlegend(tmcolor.scheme,seq(-2,2,0.5),ratio.colbar=0.3 ,xlim=c(-1.55,-1.4),ylim=c(0.5,1),align="r",cex=0.5)
#colorlegend(tmcolor.scheme,c("-2","-1.5"," ","0"," ","1.5","2"),ratio.colbar=0.3 ,xlim=c(-1.6,-1.4),ylim=c(0.5,1),align="r",cex=0.7)

	text(-1.50, 0.43, c("t(DNAm)\nCore"),cex=0.6)
	text(-1.50, -0.57,c("t(mRNA)\nBorder"),cex=0.6)
#legend(50, 150, legend=c("assignTo","submitBy"), col=1:2, lwd=3, lty=c(1,2));
#plot(mod.graph,vertex.shape="fcircle", vertex.frame.color=V(mod.graph)$vrcolor,vertex.frame.width=6,vertex.label=paste(sep="\n",E(mod.graph)$label,E(mod.graph)$mtval,E(mod.graph)$rtval),vertex.label.dist=1,edge.color="grey",edge.width=E(mod.graph)$edgewidth)

#plot(mod.graph,ertex.shape="fcircle", vertex.frame.color=V(mod.graph)$vrcolor,vertex.frame.width=6,vertex.label=E(mod.graph)$name,vertex.label.dist=1,edge.color="grey",edge.width=E(mod.graph)$edgewidth)a
}
else if(mode =="Epi"){
	# if the mode is Epi the frame need not to show
	plot(mod.graph,layout=layout.fruchterman.reingold,vertex.frame.color=NA,vertex.size=10,vertex.label=label.v,vertex.label.dist=0.6,vertex.label.cex=V(mod.graph)$label.cex,vertex.label.font=3,edge.color="grey",edge.width=E(mod.graph)$edgewidth)
#If you don't want vertices to have a frame, supply NA as the color name.
#plot(mod.graph,vertex.shape="fcircle", vertex.frame.color=V(mod.graph)$vrcolor,vertex.frame.width=6,vertex.label=paste(sep="\n",V(mod.graph)$label,V(mod.graph)$mtval,V(mod.graph)$rtval),vertex.label.dist=0.8,vertex.label.cex=label.cex,vertex.label.font=3,edge.color="grey",edge.width=E(mod.graph)$edgewidth)

	#colorlegend(trcolor.scheme,seq(-2,2,0.5),ratio.colbar=0.3 ,xlim=c(-1.55,-1.4),ylim=c(-0.5,0),align="r",cex=0.5)
	colorlegend(tmcolor.scheme,seq(-2,2,0.5),ratio.colbar=0.3 ,xlim=c(-1.55,-1.4),ylim=c(0.5,1),align="r",cex=0.5)
#colorlegend(tmcolor.scheme,c("-2","-1.5"," ","0"," ","1.5","2"),ratio.colbar=0.3 ,xlim=c(-1.6,-1.4),ylim=c(0.5,1),align="r",cex=0.7)

	text(-1.50, 0.43, c("t(DNAm)"),cex=0.6)
	#text(-1.50, -0.57,c("t(mRNA)\nBorder"),cex=0.6)


}
else if(mode =="Exp"){
        #plot(mod.graph,layout=layout.fruchterman.reingold,vertex.shape="fcircle", vertex.frame.color=V(mod.graph)$vrcolor,vertex.frame.width=4,vertex.size=10,vertex.label=label.v,vertex.label.dist=0.6,vertex.label.cex=V(mod.graph)$label.cex,vertex.label.font=3,edge.color="grey",edge.width=E(mod.graph)$edgewidth)
	plot(mod.graph,layout=layout.fruchterman.reingold,vertex.frame.color=NA,vertex.size=10,vertex.label=label.v,vertex.label.dist=0.6,vertex.label.cex=V(mod.graph)$label.cex,vertex.label.font=3,edge.color="grey",edge.width=E(mod.graph)$edgewidth)
#plot(mod.graph,vertex.shape="fcircle", vertex.frame.color=V(mod.graph)$vrcolor,vertex.frame.width=6,vertex.label=paste(sep="\n",V(mod.graph)$label,V(mod.graph)$mtval,V(mod.graph)$rtval),vertex.label.dist=0.8,vertex.label.cex=label.cex,vertex.label.font=3,edge.color="grey",edge.width=E(mod.graph)$edgewidth)

        colorlegend(trcolor.scheme,seq(-2,2,0.5),ratio.colbar=0.3 ,xlim=c(-1.55,-1.4),ylim=c(-0.5,0),align="r",cex=0.5)
        #colorlegend(tmcolor.scheme,seq(-2,2,0.5),ratio.colbar=0.3 ,xlim=c(-1.55,-1.4),ylim=c(0.5,1),align="r",cex=0.5)
#colorlegend(tmcolor.scheme,c("-2","-1.5"," ","0"," ","1.5","2"),ratio.colbar=0.3 ,xlim=c(-1.6,-1.4),ylim=c(0.5,1),align="r",cex=0.7)

        #text(-1.50, 0.43, c("t(DNAm)\nCore"),cex=0.6)
        text(-1.50, -0.57,c("t(mRNA)"),cex=0.6)


}
dev.off()
return(igraph.to.graphNEL(mod.graph));

}
