\name{DoEpiMod}
\alias{DoEpiMod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
DoEpiMod
}
\description{
Indentify  differiential methylation hotspots in the network. Edge weights in the interactome network reflect the combined differiential methylation statistics(absolute values) of the genes making up the edge. 
}
\usage{
DoEpiMod(statM.m, adj.m, nseeds = 100, gamma = 0.5, nMC = 1000, sizeR.v = c(1, 100), minsizeOUT = 10, writeOUT = TRUE, nameSTUDY = "X", ew.v = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{

Arguments:

	\item{statM.m}{
a matrix of statistics and p-values of differential methylation (one row for each gene promoter) with rownames annotated with entrez gene IDs.
}
  
	\item{adj.m}{
adjacency matrix with number of rows and columns equal to length of statM.v and statR.v, ordered in same way and with same gene identifier. The resulting graph is assumed to be connected.
}
  \item{nseeds}{
number of seeds/modules to search for. This should be a number such that P-values of significance after multiple testing is less than some reasonable FDR threshold, i.e.  0.3.
}
  \item{gamma}{
tuning parameter of spin-glass algorithm. Default value generally leads to modules in the desired size range (10-100).
}
  \item{nMC}{
number of Monte Carlo runs for establishing statistical significance of modularity values under randomisation of the molecular profiles on the network.
}
  \item{sizeR.v}{
desired size range for modules
}
  \item{minsizeOUT}{
 minimum size of modules to report as interesting
}

  \item{writeOUT}{
 a logical to indicate whether to write out tables in text format
}
  \item{nameSTUDY}{
a name for the study.
}
  \item{ew.v}{
The adjacency edge weight vector
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
