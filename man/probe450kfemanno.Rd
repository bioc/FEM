\name{probe450kfemanno}
\alias{probe450kfemanno}
\docType{data}
\title{
Processed annotation of Illumina 450k probes for use in the FEM algorithm.
}
\description{A list of entries giving various annotation information of
  the 450k probes}

\value{
\item{typeC}{a vector of 0's and 1's indicating whether a probe maps to a CpG (1) or not (0).}

\item{CGI}{a vector of 0's and 1's indicating whether a probe maps to a CGI (1) or not (0).}


\item{design}{a vector of 1's and 2's, specifying design type of probe (1 = type I, 2 = type II).}

\item{probeID}{the 450k probe IDs.}

\item{GeneGroup}{an integer vector specifying which gene region the
probe maps to. Regions are annotated as follows: 1=TSS1500, 2=TSS200, 3=5'UTR, 4=1stExon, 5=gene body, 6=3'UTR. Probes with ambiguous mappings are assigned an
NA.}

\item{eid}{mapping of 450k probes to Entrez Gene ID.}
}
\usage{data(probe450kfemanno)}
\format{
}
\details{
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
}
\keyword{datasets}
