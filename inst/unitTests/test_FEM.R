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
  #load toydata 
  data(toydata);
  #test DoFEMbi function
  DoFEMbi_test.o  <- with(toydata,
                             DoFEMbi(statM,statR,adjacency,
                             nseeds=1,gamma=0.5,nMC=1000,
                             sizeR.v=c(1,100),minsizeOUT=10,
                             writeOUT=TRUE,
                             nameSTUDY="TEST_DoFEMbi",ew.v=NULL));
  #test DoEpiMod function
  DoEpiMod_test.o<- with(toydata,
			DoEpiMod(statM,adjacency,nseeds=1,gamma=0.5,nMC=1000,
				sizeR.v=c(1,100),minsizeOUT=10,writeOUT=TRUE,
				nameSTUDY="TEST_DoEpiMod",ew.v=NULL));
  #test DoExpMod function
  DoExpMod_test.o<- with(toydata,
			DoEpiMod(statR,adjacency,nseeds=1,gamma=0.5,nMC=1000,
                                sizeR.v=c(1,100),minsizeOUT=10,writeOUT=TRUE,
                                nameSTUDY="TEST_DoExpMod",ew.v=NULL));
  #load realdata and DoFEMbi returned object fembi.o which is drived from realdata, we provide it here because it takes a little long time to generate it
  data(realdata);
  data(fembi.o);
  HAND2.mod<-fembi.o$topmod$HAND2;
  #test FemModShow function 
  HAND2.graphNEL.o=FemModShow(fembi.o$topmod$HAND2,name="HAND2",fembi.o$ew,realdata$adjacency)
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
