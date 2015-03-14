\name{Realdata}
\alias{Realdata}
\docType{data}
\title{
A dataset object derived from real DNA methylation and gene expression
data from normal endometrial and endometrial cancer tissue}
\description{
Realdata is a list containing data derived from mRNA Expression and DNA
methylation from an endometrial cancer study. See Jones et al PLoS
Med.2013 for details about the data.

The list has five entries:

Realdata$statM: a matrix of statistics and p-values of differential DNA
Methylation between endometrial cancer and normal tissue (cancer compared to normal) with rownames annotated with entrez gene IDs.

Realdata$statR: a matrix of statistics and p-values of differential RNA
Expression between endometrial cancer and normal tissue (cancer compared to normal) with rownames annotated with entrez gene IDs.

Realdata$adjacency: adjacency matrix of PPI network, with number of rows
equal to the number of rows of Realdata$statM, ordered in same way and with same gene identifier. The resulting graph is connected. 

Realdata$fembi.o: this entry represents the output of the function
DoFEMbi() (see example below). This is included so as to avoid having to
rerun the (lengthy) example from scratch.

}


\usage{data(Realdata)}
\format{
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Jones A, Teschendorff AE, Li Q, Hayward JD, Kannan A, et al. (2013) Role of dna methylation and epigenetic silencing of hand2 in endometrial cancer development. PLoS Med 10:e1001551.
}
\examples{
#data(Realdata);
#intFEM.o <- list(statM=Realdata$statM,statR=Realdata$statR,adj=Realdata$adjacency);
#DoFEMbi(intFEM.o,nseeds=100,gamma=0.5,nMC=1000,sizeR.v=c(1,100),minsizeOUT=10,writeOUT=TRUE,nameSTUDY="TEST",ew.v=NULL);
}
\keyword{datasets}
