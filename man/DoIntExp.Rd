\name{DoIntExp}
\alias{DoIntExp}
\title{
DoIntExp
}
\description{
generate the statR, adjacency, annotation matrix for the DoExpMod.
}
\usage{
DoIntExp(exp.m, phenoR.v, adj.m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
Arguments:
  
  \item{exp.m}{
normalized gene expression data matrix with rownames annotated to NCBI Entrez gene IDs. If the mapped Entrez gene IDs are not unique,
 we use the average value of the same Entrez gene ID as the expresssion value.
%%     ~~Describe \code{exp.m} here~~
}
  
  \item{phenoR.v}{
 phenotype vector corresponding to dnaM.m
%%     ~~Describe \code{phenoR.v} here~~
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
