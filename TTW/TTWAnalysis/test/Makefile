CC = g++
CFLAGS = -Wall -c -g

ROOFIT_INCLUDE := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep INCLUDE= | sed 's|INCLUDE=||')
ROOFIT_LIBDIR := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIBDIR= | sed 's|LIBDIR=||')

INCLUDES = -I. -I$(ROOTSYS)/include  -I$(ROOFIT_INCLUDE)/ -I$(CMSSW_BASE)/src/ -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -I$(CMSSW_BASE)/src/UserCode/pandolf/

ROOTSYS  ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs` 

EXTRALIBS  :=  -L$(ROOTSYS)/lib -L$(ROOFIT_LIBDIR)/ -lHtml -lMathCore -lGenVector -lMinuit -lEG -lRooFitCore -lRooFit


merge_and_setWeights: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/merge_and_setWeights.cpp
	$(CC) -Wall $(INCLUDES) -o merge_and_setWeights $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/merge_and_setWeights.cpp $(ROOTFLAG) $(EXTRALIBS)


do2ndLevel_TTW: Ntp1Analyzer.o Ntp1Analyzer_TTW.o AnalysisElectron.o AnalysisMuon.o AnalysisJet.o PUWeight.o QGLikelihoodCalculator.o do2ndLevel_TTW.cpp
	$(CC) -Wall $(INCLUDES) -o do2ndLevel_TTW do2ndLevel_TTW.cpp Ntp1Analyzer.o PUWeight.o QGLikelihoodCalculator.o Ntp1Analyzer_TTW.o AnalysisElectron.o AnalysisMuon.o AnalysisJet.o $(ROOTFLAG) 

do2ndLevel_TTW_DATA: Ntp1Analyzer.o Ntp1Analyzer_TTW.o AnalysisElectron.o AnalysisMuon.o AnalysisJet.o PUWeight.o do2ndLevel_TTW_DATA.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o do2ndLevel_TTW_DATA do2ndLevel_TTW_DATA.cpp Ntp1Analyzer.o Ntp1Analyzer_TTW.o AnalysisElectron.o AnalysisMuon.o AnalysisJet.o PUWeight.o $(ROOTFLAG)




finalize_TTW: Ntp1Finalizer.o Ntp1Finalizer_TTW.o finalize_TTW.cpp  QGLikelihoodCalculator.o AnalysisJet.o PUWeight.o
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -I$(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood -o finalize_TTW finalize_TTW.cpp Ntp1Finalizer.o Ntp1Finalizer_TTW.o QGLikelihoodCalculator.o AnalysisJet.o PUWeight.o $(ROOTFLAG) $(EXTRALIBS)

finalize_TTWTrilepton: Ntp1Finalizer.o Ntp1Finalizer_TTWTrilepton.o finalize_TTWTrilepton.cpp  QGLikelihoodCalculator.o AnalysisJet.o BTagSFUtil.o PUWeight.o
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -I$(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood -o finalize_TTWTrilepton finalize_TTWTrilepton.cpp Ntp1Finalizer.o Ntp1Finalizer_TTWTrilepton.o QGLikelihoodCalculator.o AnalysisJet.o BTagSFUtil.o PUWeight.o $(ROOTFLAG) $(EXTRALIBS)

finalize_TTWDilepton: Ntp1Finalizer.o Ntp1Finalizer_TTWDilepton.o finalize_TTWDilepton.cpp  QGLikelihoodCalculator.o AnalysisJet.o BTagSFUtil.o PUWeight.o
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -I$(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood -o finalize_TTWDilepton finalize_TTWDilepton.cpp Ntp1Finalizer.o Ntp1Finalizer_TTWDilepton.o QGLikelihoodCalculator.o AnalysisJet.o BTagSFUtil.o PUWeight.o $(ROOTFLAG) $(EXTRALIBS)


drawTTWTrilepton: DrawBase.o fitTools.o drawTTWTrilepton.cpp
	$(CC) -Wall $(INCLUDES) -o drawTTWTrilepton drawTTWTrilepton.cpp DrawBase.o fitTools.o $(ROOTFLAG) $(EXTRALIBS)

drawTTWTrilepton_DATA: DrawBase.o fitTools.o drawTTWTrilepton_DATA.cpp
	$(CC) -Wall $(INCLUDES) -o drawTTWTrilepton_DATA drawTTWTrilepton_DATA.cpp DrawBase.o fitTools.o $(ROOTFLAG) $(EXTRALIBS)

drawTTW_leptIso: DrawBase.o fitTools.o drawTTW_leptIso.cpp
	$(CC) -Wall $(INCLUDES) -o drawTTW_leptIso drawTTW_leptIso.cpp DrawBase.o fitTools.o $(ROOTFLAG) $(EXTRALIBS)

drawTTH: DrawBase.o fitTools.o drawTTH.cpp
	$(CC) -Wall $(INCLUDES) -o drawTTH drawTTH.cpp DrawBase.o fitTools.o $(ROOTFLAG) $(EXTRALIBS)




Ntp1Analyzer.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/Ntp1Analyzer.C
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/Ntp1Analyzer.C $(ROOTFLAG)

Ntp1Analyzer_TTW.o: Ntp1Analyzer_TTW.cc QGLikelihoodCalculator.o
	$(CC) $(CFLAGS) $(INCLUDES) QGLikelihoodCalculator.o Ntp1Analyzer_TTW.cc $(ROOTFLAG)


Ntp1Finalizer.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/Ntp1Finalizer.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/ $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/Ntp1Finalizer.C $(ROOTFLAG)

Ntp1Finalizer_TTW.o: Ntp1Finalizer_TTW.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/ -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -I$(CMSSW_BASE)/src/UserCode/emanuele/CommonTools Ntp1Finalizer_TTW.cc $(ROOTFLAG)

Ntp1Finalizer_TTWTrilepton.o: Ntp1Finalizer_TTWTrilepton.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/ -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -I$(CMSSW_BASE)/src/UserCode/emanuele/CommonTools Ntp1Finalizer_TTWTrilepton.cc $(ROOTFLAG)

Ntp1Finalizer_TTWDilepton.o: Ntp1Finalizer_TTWDilepton.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/ -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -I$(CMSSW_BASE)/src/UserCode/emanuele/CommonTools Ntp1Finalizer_TTWDilepton.cc $(ROOTFLAG)




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

QGLikelihoodCalculator.o: $(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood/src/QGLikelihoodCalculator.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood/ $(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood/src/QGLikelihoodCalculator.cc $(ROOTFLAG)


BTagSFUtil.o: $(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/src/BTagSFUtil.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/ $(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/src/BTagSFUtil.cc $(ROOTFLAG)

PUWeight.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/PUWeight.C
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/PUWeight.C $(ROOTFLAG)


AnalysisJet.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisJet.cc 
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisJet.cc $(ROOTFLAG)

AnalysisElectron.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisElectron.cc AnalysisLepton.o
	$(CC) $(CFLAGS) $(INCLUDES) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisElectron.cc $(ROOTFLAG)

AnalysisMuon.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisMuon.cc AnalysisLepton.o
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisMuon.cc $(ROOTFLAG)

AnalysisLepton.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisLepton.cc
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisLepton.cc $(ROOTFLAG)

AnalysisNeutrino.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisNeutrino.h
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/AnalysisNeutrino.h $(ROOTFLAG)



clean:
	rm *.o
