CC = g++
CFLAGS = -Wall -c -g

ROOFIT_INCLUDE := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep INCLUDE= | sed 's|INCLUDE=||')
ROOFIT_LIBDIR := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIBDIR= | sed 's|LIBDIR=||')

INCLUDES = -I. -I$(ROOTSYS)/include  -I$(ROOFIT_INCLUDE)/ -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -I$(CMSSW_BASE)/src/UserCode/pandolf/

ROOTSYS  ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs` 

EXTRALIBS  :=  -L$(ROOTSYS)/lib -L$(ROOFIT_LIBDIR)/ -lHtml -lMathCore -lGenVector -lMinuit -lEG -lRooFitCore -lRooFit


merge_and_setWeights: merge_and_setWeights.cpp
	$(CC) -Wall $(INCLUDES) -o merge_and_setWeights merge_and_setWeights.cpp $(ROOTFLAG) $(EXTRALIBS)


do2ndLevel_TTZ: Ntp1Analyzer.o Ntp1Analyzer_TTZ.o AnalysisElectron.o AnalysisMuon.o PUWeight.o do2ndLevel_TTZ.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o do2ndLevel_TTZ do2ndLevel_TTZ.cpp Ntp1Analyzer.o PUWeight.o Ntp1Analyzer_TTZ.o AnalysisElectron.o AnalysisMuon.o $(ROOTFLAG) 

do2ndLevel_TTZ_DATA: Ntp1Analyzer.o Ntp1Analyzer_TTZ.o AnalysisElectron.o AnalysisMuon.o PUWeight.o do2ndLevel_TTZ_DATA.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o do2ndLevel_TTZ_DATA do2ndLevel_TTZ_DATA.cpp Ntp1Analyzer.o Ntp1Analyzer_TTZ.o AnalysisElectron.o AnalysisMuon.o PUWeight.o $(ROOTFLAG)




finalize_TTZ: Ntp1Finalizer.o Ntp1Finalizer_TTZ.o finalize_TTZ.cpp fitTools.o TAbsFitParticle.o TAbsFitConstraint.o TFitConstraintM.o TFitParticleEtEtaPhi.o TKinFitter.o QGLikelihoodCalculator.o HelicityLikelihoodDiscriminant.o RooPentaSpinTwo.o RooBphione.o RooBtheta.o RooBkg2l2jHL.o RooSpinZero5DV2.o RooBkgd2L2JV2.o DiJetKinFitter.o AnalysisJet.o BTagSFUtil.o
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -I$(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood -o finalize_TTZ finalize_TTZ.cpp Ntp1Finalizer.o Ntp1Finalizer_TTZ.o fitTools.o TAbsFitParticle.o TAbsFitConstraint.o TFitConstraintM.o TFitParticleEtEtaPhi.o TKinFitter.o QGLikelihoodCalculator.o HelicityLikelihoodDiscriminant.o RooPentaSpinTwo.o RooBphione.o RooBtheta.o RooBkg2l2jHL.o RooSpinZero5DV2.o RooBkgd2L2JV2.o DiJetKinFitter.o AnalysisJet.o BTagSFUtil.o $(ROOTFLAG) -lMinuit


drawTTZ: DrawBase.o fitTools.o drawTTZ.cpp
	$(CC) -Wall $(INCLUDES) -o drawTTZ drawTTZ.cpp DrawBase.o fitTools.o $(ROOTFLAG) $(EXTRALIBS)

drawTTZ_DATA: DrawBase.o fitTools.o drawTTZ_DATA.cpp
	$(CC) -Wall $(INCLUDES) -o drawTTZ_DATA drawTTZ_DATA.cpp DrawBase.o fitTools.o $(ROOTFLAG) $(EXTRALIBS)




Ntp1Analyzer.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/Ntp1Analyzer.C
	$(CC) $(CFLAGS) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/Ntp1Analyzer.C $(ROOTFLAG)

Ntp1Analyzer_TTZ.o: Ntp1Analyzer_TTZ.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -I$(CMSSW_BASE)/src/UserCode/emanuele/CommonTools Ntp1Analyzer_TTZ.C $(ROOTFLAG)


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

DiJetKinFitter.o: $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/DiJetKinFitter.cc 
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/DiJetKinFitter.cc $(ROOTFLAG)

LeptonNeutrinoKinFitter.o: $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/LeptonNeutrinoKinFitter.cc 
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/LeptonNeutrinoKinFitter.cc $(ROOTFLAG)

GlobalFitter.o: $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/GlobalFitter.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit $(CMSSW_BASE)/src/UserCode/pandolf/KinematicFit/GlobalFitter.cc $(ROOTFLAG)


DrawBase.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/DrawBase.C fitTools.o
	$(CC) $(CFLAGS) $(INCLUDES) fitTools.o $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/DrawBase.C $(ROOTFLAG) $(EXTRALIBS) 

fitTools.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/fitTools.C
	$(CC) $(CFLAGS) $(INCLUDES) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/fitTools.C $(ROOTFLAG) $(EXTRALIBS)

QGLikelihoodCalculator.o: $(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood/QGLikelihoodCalculator.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood/ $(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood/QGLikelihoodCalculator.C $(ROOTFLAG)

HelicityLikelihoodDiscriminant.o: $(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/HelicityLikelihoodDiscriminant.cc RooPentaSpinTwo.o RooBphione.o RooBtheta.o RooBkg2l2jHL.o RooSpinZero5DV2.o RooBkgd2L2JV2.o 
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/ RooPentaSpinTwo.o RooBphione.o RooBtheta.o RooBkg2l2jHL.o RooSpinZero5DV2.o RooBkgd2L2JV2.o $(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/HelicityLikelihoodDiscriminant.cc $(ROOTFLAG)

RooBphione.o: $(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs/RooBphione.cc 
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs $(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs/RooBphione.cc $(ROOTFLAG)

RooBtheta.o: $(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs/RooBtheta.cc 
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs $(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs/RooBtheta.cc $(ROOTFLAG)

RooPentaSpinTwo.o: $(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs/RooPentaSpinTwo.cc 
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs $(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs/RooPentaSpinTwo.cc $(ROOTFLAG)

RooBkg2l2jHL.o: $(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs/RooBkg2l2jHL.cc 
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs $(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs/RooBkg2l2jHL.cc $(ROOTFLAG)

RooSpinZero5DV2.o: $(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs/RooSpinZero5DV2.cc 
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs $(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs/RooSpinZero5DV2.cc $(ROOTFLAG)

RooBkgd2L2JV2.o: $(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs/RooBkgd2L2JV2.cc 
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs $(CMSSW_BASE)/src/UserCode/pandolf/HelicityLikelihoodDiscriminant/PDFs/RooBkgd2L2JV2.cc $(ROOTFLAG)

RooFermi.o: RooFermi.cc
	$(CC) $(CFLAGS) $(INCLUDES) RooFermi.cc $(ROOTFLAG) $(EXTRALIBS)

RooCruijffPdf.o: RooCruijffPdf.cc
	$(CC) $(CFLAGS) $(INCLUDES) RooCruijffPdf.cc $(ROOTFLAG) $(EXTRALIBS)

BTagSFUtil.o: $(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/BTagSFUtil.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/ $(CMSSW_BASE)/src/UserCode/pandolf/BTagSFUtil/BTagSFUtil.cc $(ROOTFLAG)

PUWeight.o: PUWeight.C
	$(CC) $(CFLAGS) PUWeight.C $(ROOTFLAG)


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
