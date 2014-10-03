\name{DoIntEpi450k}
\alias{DoIntEpi450k}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
DoIntEpi450k
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Generate differentcial methylation statistics using 450K methylation matrix.
}
\usage{
DoIntEpi450k(dnaM.m, phenoM.v, adj.m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

Arguments:

  \item{dnaM.m}{
normalised DNA methylation 450k data matrix, with rownames annotated to 450k probe IDs.
}
  \item{phenoM.v}{
phenotype vector corresponding to dnaM.m
%%     ~~Describe \code{phenoM.v} here~~
}
  \item{adj.m}{
adjacency matrix of a network of relations (e.g. PPI network) with rownames/colnames annotated to NCBI Entrez gene IDs.
Note: The PPI network can be derived from the Pathway Commons resource\cite{Cerami2011} and follows the procedure described in \cite{West2013}. The PIN used in previous papers  is available at http://sourceforge.net/projects/funepimod. The PPI network consists of 8434 genes annotated to NCBI Entrez identifiers, and is sparse containing 303600 documented interactions (edges).  If the user wishes they can use a different PPI network or generate statR and statM using different method.
%%     ~~Describe \code{adj.m} here~~
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
