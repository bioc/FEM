\name{DoFEMbi}
\alias{DoFEMbi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identifies interactome hotspots of differential promoter DNAm
methylation and differential mRNA expression. 
}
\description{
DoFEMbi identifies interactome hotspots of differential promoter
methylation and differential expression in relation to a phenotype of interest, where an inverse association between methylation and gene expression is assumed.
}
\usage{
DoFEMbi(intFEM.o, nseeds = 10, gamma = 1.1, nMC = 100, 
  minSize=5, maxSize=50, writeOUT = TRUE, nameSTUDY = "X", pvth=0.05, ncores=4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

  \item{intFEM.o}{The output of the DoIntFEM450k function.}
  \item{nseeds}{An integer specifying the number of seeds and therefore
  modules to search for. By default this number is 10.}

  \item{gamma}{A parameter of the spin-glass algorithm, which determines
  the average module size. Default value generally leads to modules in
  the desired size range (10-100 genes).}

  \item{nMC}{Number of Monte Carlo runs for establishing statistical significance of modularity values under randomisation of the molecular profiles on the network.}
  
  \item{minSize}{
 Minimum size of modules to report as interesting.
}
  \item{maxSize}{
 Maximum size of modules to report as interesting.
}

  \item{writeOUT}{A logical to indicate whether to write out tables in text format.}

  \item{nameSTUDY}{A name for the study, to be used as label in the
    output files.}

  \item{pvth}{P value threshold for detecting significant modules (default is 0.05).}

  \item{ncores}{No. of parallel processing computer cores (default value is 2)}


}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with following entries:
\item{size}{A vector of inferred module sizes for each of the ntop seeds.}
\item{mod}{A vector of associated modularities.}
\item{pv}{A vector of associated significance P-values with resolution of nMC}
\item{sigmod}{Index positions of significant modules of size at least minSize}
\item{fem}{A summary matrix of the selected modules.}
\item{topmod}{A list of summary matrices for each of the selected
  modules.}
\item{sgc}{A list of the spin-glass module detection algorithm for each
  seed.}
\item{ew}{The edge-weight vector of the integrated network.}
\item{adj}{ adjacency matrix of the maximally connected integrated
  network (at present only maximally connected subnetwork is used).It is same to intFEM.o$adj, and wil be used for FemModShow function}
}

\references{
A systems-level integrative framework for genome-wide DNA methylation and gene expression data identifies differential gene expression modules under epigenetic control. Jiao Y, Widschwendter M, Teschendorff AE. Bioinformatics. 2014;30(16):2360-2366
}
\author{
"Yinming Jiao"<20907099@zju.edu.cn>, "Andrew E Teschendorff"<andrew@picb.ac.cn>
}

\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(Toydata);
intFEM.o <- list(statM=Toydata$statM,statR=Toydata$statR,adj=Toydata$adj);
DoFEMbi(intFEM.o,nseeds=1,gamma=1.1,nMC=1000,
 minSize=10,maxSize=50,writeOUT=TRUE,nameSTUDY="TEST");
#You can also test on the Realdata contains matched DNA methylation and RNA Expression of 17 normal and 118 cancer samples. Since running on the realdata is time-consuming, we comment it out.   
#data(Realdata);
#intFEM.o <- list(statM=Realdata$statM,statR=Realdata$statR,adj=Realdata$adjacency);
#DoFEMbi(intFEM.o,nseeds=100,gamma=0.5,nMC=1000,nameSTUDY="TEST");
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
