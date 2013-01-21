CC = g++
CFLAGS = -Wall -c -g

ROOFIT_INCLUDE := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep INCLUDE= | sed 's|INCLUDE=||')
ROOFIT_LIBDIR := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIBDIR= | sed 's|LIBDIR=||')

INCLUDES = -I. -I$(ROOTSYS)/include  -I$(ROOFIT_INCLUDE)/ -I$(CMSSW_BASE)/src -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -I$(CMSSW_BASE)/src/UserCode/pandolf/ -I$(CMSSW_BASE)/src/CondFormats/JetMETObjects/interface

ROOTSYS  ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs` 

EXTRALIBS  :=  -L$(ROOTSYS)/lib -L$(ROOFIT_LIBDIR)/ -L$(CMSSW_BASE)/lib/slc5_amd64_gcc462/ -lHtml -lMathCore -lGenVector -lMinuit -lEG -lRooFitCore -lRooFit



provaQG: QGLikelihoodCalculator.o Bins.o provaQG.cpp
	$(CC) -Wall $(INCLUDES) -o provaQG provaQG.cpp QGLikelihoodCalculator.o Bins.o $(ROOTFLAG) $(EXTRALIBS) 

QGLikelihoodCalculator.o: $(CMSSW_BASE)/src/QGLikelihood/src/QGLikelihoodCalculator.cc
	$(CC) $(CFLAGS) $(INCLUDES) $(CMSSW_BASE)/src/QGLikelihood/src/QGLikelihoodCalculator.cc $(ROOTFLAG)


Bins.o: $(CMSSW_BASE)/src/QGLikelihood/src/Bins.cc
	$(CC) $(CFLAGS) $(INCLUDES) $(CMSSW_BASE)/src/QGLikelihood/src/Bins.cc $(ROOTFLAG)





clean:
	rm *.o
