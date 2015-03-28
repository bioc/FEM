\name{Toydata}
\alias{Toydata}
\docType{data}
\title{
Artificial data set. 
}

\description{A list object containing artificial data to test the FEM
  algorithm. It contains artifically created statistics of differential
  DNA methylation, mRNA expression, as well as an artificial adjacency
  matrix to illustrate the application of FEM.}

\value{
\item{statM}{a matrix of simulated statistics and p-values of differential DNA Methylation with 84 rows (genes).}

\item{statR}{a matrix of simulated statistics and p-values of differential RNA Expression with 84 rows (genes).}

\item{adjacency}{artificial adjacency matrix of 84 rows and columns,
  ordered according to the rows of previous two matrices. Adjacency
  matrix is a random graph except for containing a clique defined by 10 nodes.}

\item{tennodes}{a vector of 10 integer indices, labeling the nodes/genes which form the clique in the network.}


}
\usage{data(Toydata)}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
A systems-level integrative framework for genome-wide DNA methylation and gene expression data identifies differential gene expression modules under epigenetic control. Jiao Y, Widschwendter M, Teschendorff AE. Bioinformatics. 2014;30(16):2360-2366
}
\author{
"Yinming Jiao"<20907099@zju.edu.cn>, "Andrew E Teschendorff"<andrew@picb.ac.cn>
}
\examples{
data(Toydata);
intFEM.o <- list(statM=Toydata$statM,statR=Toydata$statR,adj=Toydata$adjacency);
DoFEMbi(intFEM.o,nseeds=1,gamma=0.5,nMC=1000,sizeR.v=c(1,100),minsizeOUT=10,
 writeOUT=TRUE,nameSTUDY="TEST",ew.v=NULL);
}
\keyword{datasets}
