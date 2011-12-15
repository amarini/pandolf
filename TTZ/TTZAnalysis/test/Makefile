CC = g++
CFLAGS = -Wall -c -g

ROOFIT_INCLUDE := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep INCLUDE= | sed 's|INCLUDE=||')
ROOFIT_LIBDIR := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIBDIR= | sed 's|LIBDIR=||')

INCLUDES = -I. -I$(ROOTSYS)/include  -I$(ROOFIT_INCLUDE)/ -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -I$(CMSSW_BASE)/src/UserCode/pandolf/

ROOTSYS  ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs` 

EXTRALIBS  :=  -L$(ROOTSYS)/lib -L$(ROOFIT_LIBDIR)/ -lHtml -lMathCore -lGenVector -lMinuit -lEG -lRooFitCore -lRooFit


merge_and_setWeights: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/merge_and_setWeights.cpp
	$(CC) -Wall $(INCLUDES) -o merge_and_setWeights $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/merge_and_setWeights.cpp $(ROOTFLAG) $(EXTRALIBS)


do2ndLevel_TTZ: Ntp1Analyzer.o Ntp1Analyzer_TTZ.o AnalysisElectron.o AnalysisMuon.o AnalysisJet.o PUWeight.o do2ndLevel_TTZ.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o do2ndLevel_TTZ do2ndLevel_TTZ.cpp Ntp1Analyzer.o PUWeight.o Ntp1Analyzer_TTZ.o AnalysisElectron.o AnalysisMuon.o AnalysisJet.o $(ROOTFLAG) 

do2ndLevel_TTZ_DATA: Ntp1Analyzer.o Ntp1Analyzer_TTZ.o AnalysisElectron.o AnalysisMuon.o AnalysisJet.o PUWeight.o do2ndLevel_TTZ_DATA.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o do2ndLevel_TTZ_DATA do2ndLevel_TTZ_DATA.cpp Ntp1Analyzer.o Ntp1Analyzer_TTZ.o AnalysisElectron.o AnalysisMuon.o AnalysisJet.o PUWeight.o $(ROOTFLAG)




finalize_TTZTrilepton: Ntp1Finalizer.o Ntp1Finalizer_TTZTrilepton.o finalize_TTZTrilepton.cpp  QGLikelihoodCalculator.o AnalysisJet.o BTagSFUtil.o PUWeight.o
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -I$(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood -o finalize_TTZTrilepton finalize_TTZTrilepton.cpp Ntp1Finalizer.o Ntp1Finalizer_TTZTrilepton.o QGLikelihoodCalculator.o AnalysisJet.o BTagSFUtil.o PUWeight.o $(ROOTFLAG) $(EXTRALIBS)


drawTTZTrilepton: DrawBase.o fitTools.o drawTTZTrilepton.cpp
	$(CC) -Wall $(INCLUDES) -o drawTTZTrilepton drawTTZTrilepton.cpp DrawBase.o fitTools.o $(ROOTFLAG) $(EXTRALIBS)

drawTTZTrilepton_DATA: DrawBase.o fitTools.o drawTTZTrilepton_DATA.cpp
	$(CC) -Wall $(INCLUDES) -o drawTTZTrilepton_DATA drawTTZTrilepton_DATA.cpp DrawBase.o fitTools.o $(ROOTFLAG) $(EXTRALIBS)




Ntp1Analyzer.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/Ntp1Analyzer.C
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/Ntp1Analyzer.C $(ROOTFLAG)

Ntp1Analyzer_TTZ.o: Ntp1Analyzer_TTZ.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -I$(CMSSW_BASE)/src/UserCode/emanuele/CommonTools Ntp1Analyzer_TTZ.cc $(ROOTFLAG)


Ntp1Finalizer.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/Ntp1Finalizer.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/ $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/Ntp1Finalizer.C $(ROOTFLAG)

Ntp1Finalizer_TTZTrilepton.o: Ntp1Finalizer_TTZTrilepton.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/ -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -I$(CMSSW_BASE)/src/UserCode/emanuele/CommonTools Ntp1Finalizer_TTZTrilepton.cc $(ROOTFLAG)




TFitParticleEtEtaPhi.o: $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TFitParticleEtEtaPhi.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TFitParticleEtEtaPhi.cc $(ROOTFLAG)

TAbsFitParticle.o: $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TAbsFitParticle.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TAbsFitParticle.cc $(ROOTFLAG)

TFitConstraintM.o: $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TFitConstraintM.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TFitConstraintM.cc $(ROOTFLAG)

TFitConstraintMGaus.o: $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TFitConstraintMGaus.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TFitConstraintMGaus.cc $(ROOTFLAG)

TFitConstraintEp.o: $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TFitConstraintEp.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TFitConstraintEp.cc $(ROOTFLAG)

TAbsFitConstraint.o: $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TAbsFitConstraint.cc 
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TAbsFitConstraint.cc $(ROOTFLAG)

TKinFitter.o: $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TKinFitter.cc 
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/TKinFitter.cc $(ROOTFLAG)



DrawBase.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/DrawBase.C fitTools.o
	$(CC) $(CFLAGS) $(INCLUDES) fitTools.o $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/DrawBase.C $(ROOTFLAG) $(EXTRALIBS) 

fitTools.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/fitTools.C
	$(CC) $(CFLAGS) $(INCLUDES) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/fitTools.C $(ROOTFLAG) $(EXTRALIBS)

QGLikelihoodCalculator.o: $(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood/QGLikelihoodCalculator.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood/ $(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood/QGLikelihoodCalculator.C $(ROOTFLAG)


BTagSFUtil.o: $(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/BTagSFUtil.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/ $(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/BTagSFUtil.cc $(ROOTFLAG)

PUWeight.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/PUWeight.C
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/PUWeight.C $(ROOTFLAG)


AnalysisJet.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisJet.cc 
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisJet.cc $(ROOTFLAG)

AnalysisElectron.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisElectron.cc AnalysisLepton.o
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisElectron.cc $(ROOTFLAG)

AnalysisMuon.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisMuon.cc AnalysisLepton.o
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisMuon.cc $(ROOTFLAG)

AnalysisLepton.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisLepton.cc
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisLepton.cc $(ROOTFLAG)

AnalysisNeutrino.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisNeutrino.h
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisNeutrino.h $(ROOTFLAG)



clean:
	rm *.o
