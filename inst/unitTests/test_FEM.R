## This is a very trivial demo  of
## the RUnit test case execution system:
## ---------------------------------

## functions to be tested (usually defined in a different
## file from where the test cases are located):




## test functions:
## ---------------------


.setUp <- function() {  ## called before each test case, see also .tearDown()
  print(".setUp")
}


test.FEM <- function() {
  #load Toydata 
data(Toydata);
intFEM.o <- list(statM=Toydata$statM,statR=Toydata$statR,adj=Toydata$adj);
DoFEMbi(intFEM.o,nseeds=1,gamma=0.5,nMC=1000,sizeR.v=c(1,100),minsizeOUT=10,writeOUT=TRUE,nameSTUDY="TEST",ew.v=NULL);
  #test DoFEMbi function
}




## How to run the tests (do not uncomment in this file,
## but execute the commands at the R prompt):
## All you have to do is to adapt the directory locations.
## ------------------------------------------------

## define the test suite:
#testsuite.cf <- defineTestSuite("cfConversion", dirs="directoryOfThisFile")

## run test suite:
#testResult <- runTestSuite(testsuite.cf)

## print text protocol to console:
#printTextProtocol(testResult)

## print HTML version to a file:
#printHTMLProtocol(testResult, fileName="someFileName.html")

## In this case we also have a shortcut
#runTestFile("directoryOfThisFile/runitcfConversion.r")
