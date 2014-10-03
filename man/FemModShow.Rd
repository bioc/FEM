\name{FemModShow}
\alias{FemModShow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
FemModShow
}
\description{
generate particular module net which is from FEM result object such as fembi.o which can be loaded by "data(fembi.o)". and also it will return an igraph object.
}
\usage{
FemModShow(mod, name = "mod", edgeweight, adjacency,mode="integration")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mod}{
particular module of the FEM result object
}
  \item{name}{
the name of the module
}
  \item{edgeweight}{
FEM result object's edgeweight
}
  \item{adjacency}{
the whole net adjacency matrix
%%     ~~Describe \code{realdata} here~~
}
\item{mode}{
There are three mode, "integration", "Epi","Epx". "integration" means the module is from DoFEMbi, "Epi" means the module is from DoEpiMod, "Exp" means the module is from DoExpMod. 
%%     ~~Describe \code{realdata} here~~
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
data(fembi.o)
data(realdata)
FemModShow(fembi.o$topmod$HAND2,name="HAND2",fembi.o$ew,realdata$adjacency)
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
