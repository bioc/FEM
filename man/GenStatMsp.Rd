\name{GenStatMsp}
\alias{GenStatMsp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generation of differential DNA methylation statistics for SingleProbe Mode
}
\description{
Given an Illumina 450k or EPIC (850k) data matrix and a phenotype vector, similiar to 
the function of GenStatM, it will generate statistics of differential DNA methylation 
using limma, at the gene-level, for all pairwise comparisons of phenotype levels. However,
instead of the average beta values for the probes in each region, it will select the most 
significant differentially methylated CpGs among TSS200 and 1st Exon regions, if not present
then select the the most significant differentially methylated CpGs among TSS1500
}
\usage{
GenStatMsp(dnaM.m,pheno.v,"chiptype")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dnaM.m}{
The DNA methylation beta valued data matrix with rownames annotated to
Illumina 450k or EPIC(850k) probe IDs
}
  \item{pheno.v}{
The phenotype vector. 
}
  \item{chiptype}{
A parameter specifying the the input data matrix, it should be either "450k" for
Illumina 450k matrix or "EPIC" for Illumina EPIC matrix. default for "450k"
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{top}{A list of matrices, from the limma output, ranking genes
  according to differential methylation between two phenotypes as
  specified in the contrasts matrix, i.e. top[[i]] will contain the
  ranking for the i'th contrast.}
\item{cont}{The contrasts matrix, with columns labeling the contrasts.}
\item{probeID}{A list of vectors, with each entry a vector of probeIDs that present 
the most differentially methylated cpg associated with each gene, i.e. probeID[[i]] 
will contain the vector for the i'th contrast.}
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

}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
