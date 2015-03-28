\name{GenStatM}
\alias{GenStatM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generation of differential DNA methylation statistics
}
\description{
Given an Illumina 450k data matrix and a phenotype vector, it will
generate statistics of differential DNA methylation using limma, at the
gene-level, for all pairwise comparisons of phenotype levels.
}
\usage{
GenStatM(dnaM.m,pheno.v)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dnaM.m}{
The DNA methylation beta valued data matrix with rownames annotated to
Illumina 450k probe IDs
}
  \item{pheno.v}{
The phenotype vector. 
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
\item{avbeta}{The beta-valued DNAm data matrix at the gene-level,
  following the procedure as described in Jiao et al Bioinformatics 2014.}
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
