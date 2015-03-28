\name{GenStatR}
\alias{GenStatR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generation of differential expression statistics
}
\description{
Given a gene expression data matrix annotated to entrez gene IDs and a
phenotype vector, it will generate statistics of differential expression
using limma, at the gene-level, for all pairwise comparisons of phenotype levels.
}
\usage{
GenStatR(exp.m,pheno.v)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{exp.m}{
The gene expression data matrix with rownames annotated to entrez gene
IDs. If there are multiple rows with the same gene ID, function will
average these.
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
  according to differential expression between two phenotypes as
  specified in the contrasts matrix, i.e. top[[i]] will contain the
  ranking for the i'th contrast.}
\item{cont}{The contrasts matrix, with columns labeling the contrasts.}
\item{avexp}{The averaged expression matrix at the gene-level.}
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
