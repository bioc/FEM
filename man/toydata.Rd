\name{toydata}
\alias{toydata}
\docType{data}
\title{
toydata
}
\description{
Artifical created statitics of Methylation, RNA Expression, and also the adjacency matrix  and annotation matrix. Thees data are used to test and prove that FEM's ability to find hotspot or module based on inverse association between methylation and gene expression.
}
\usage{data(toydata)}
\format{
  The format is:
List of 4
 $ statM     : num [1:84, 1:2] -0.06511 0.00116 0.19583 3.93402 -0.0254 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:84] "1" "2" "3" "4" ...
  .. ..$ : NULL
 $ statR     : num [1:84, 1:2] -0.0959 -0.033 0.1779 -2.5759 -0.1286 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:84] "1" "2" "3" "4" ...
  .. ..$ : NULL
 $ adjacency :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
  .. ..@ i       : int [1:300] 26 30 79 5 40 3 10 18 27 29 ...
  .. ..@ p       : int [1:85] 0 2 3 5 17 20 21 23 25 27 ...
  .. ..@ Dim     : int [1:2] 84 84
  .. ..@ Dimnames:List of 2
  .. .. ..$ : chr [1:84] "1" "2" "3" "4" ...
  .. .. ..$ : NULL
  .. ..@ x       : num [1:300] 1 1 1 1 1 1 1 1 1 1 ...
  .. ..@ factors : list()
 $ annotation: chr [1:84, 1:2] "1" "2" "3" "4" ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:2] "EentrezID" "GeneSymbol"
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(toydata)
## maybe str(toydata) ; plot(toydata) ...
}
\keyword{datasets}
