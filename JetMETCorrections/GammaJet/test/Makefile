CC = g++
CFLAGS = -Wall -c

ROOFIT_INCLUDE := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep INCLUDE= | sed 's|INCLUDE=||')
ROOFIT_LIBDIR := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIBDIR= | sed 's|LIBDIR=||')

INCLUDES = -I. -I$(ROOTSYS)/include  -I$(ROOFIT_INCLUDE)/ -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -I$(CMSSW_BASE)/src/UserCode/pandolf/

ROOTSYS  ?= ERROR_RootSysIsNotDefined

ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs` -lMinuit -lRooFitCore -lRooFit

EXTRALIBS  :=  -L$(ROOTSYS)/lib -L$(ROOFIT_LIBDIR)/ -lHtml -lMathCore -lGenVector -lMinuit -lEG -lRooFitCore -lRooFit



do2ndLevel_DiJet: AnalysisJet.o TreeAnalyzer.o TreeAnalyzer_DiJet.o
	$(CC) -Wall -o do2ndLevel_DiJet do2ndLevel_DiJet.cpp AnalysisJet.o TreeAnalyzer.o TreeAnalyzer_DiJet.o $(ROOTFLAG)

do2ndLevel_PhotonJet: AnalysisJet.o AnalysisPhoton.o TreeAnalyzer.o TreeAnalyzer_PhotonJet.o do2ndLevel_PhotonJet.cpp fitTools.o
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o do2ndLevel_PhotonJet do2ndLevel_PhotonJet.cpp AnalysisJet.o AnalysisPhoton.o TreeAnalyzer.o TreeAnalyzer_PhotonJet.o fitTools.o $(ROOTFLAG) -lMinuit

#do2ndLevel_PhotonJet_batch: AnalysisJet.o AnalysisPhoton.o TreeAnalyzer.o TreeAnalyzer_PhotonJet.o Utilities.o JetCorrectorParameters.o SimpleJetCorrector.o FactorizedJetCorrector.o fitTools.o do2ndLevel_PhotonJet_batch.cpp
#	$(CC) -Wall -o do2ndLevel_PhotonJet_batch do2ndLevel_PhotonJet_batch.cpp AnalysisJet.o AnalysisPhoton.o TreeAnalyzer.o TreeAnalyzer_PhotonJet.o Utilities.o JetCorrectorParameters.o SimpleJetCorrector.o FactorizedJetCorrector.o fitTools.o $(ROOTFLAG) -lMinuit

do2ndLevel_PhotonJet_batch: AnalysisJet.o AnalysisPhoton.o TreeAnalyzer.o TreeAnalyzer_PhotonJet.o  fitTools.o QGLikelihoodCalculator.o do2ndLevel_PhotonJet_batch.cpp
	$(CC) -Wall $(INCLUDES) -o do2ndLevel_PhotonJet_batch do2ndLevel_PhotonJet_batch.cpp AnalysisJet.o AnalysisPhoton.o TreeAnalyzer.o TreeAnalyzer_PhotonJet.o fitTools.o QGLikelihoodCalculator.o $(ROOTFLAG) $(EXTRALIBS)

do2ndLevel_VHgg_batch: AnalysisJet.o AnalysisPhoton.o TreeAnalyzer.o TreeAnalyzer_VHgg.o  do2ndLevel_VHgg_batch.cpp
	$(CC) -Wall -o do2ndLevel_VHgg_batch do2ndLevel_VHgg_batch.cpp AnalysisJet.o AnalysisPhoton.o TreeAnalyzer.o TreeAnalyzer_VHgg.o $(ROOTFLAG) $(EXTRALIBS)

do2ndLevel_MultiJet_batch: AnalysisJet.o TreeAnalyzer.o TreeAnalyzer_MultiJet.o  QGLikelihoodCalculator.o do2ndLevel_MultiJet_batch.cpp
	$(CC) -Wall $(INCLUDES) -o do2ndLevel_MultiJet_batch do2ndLevel_MultiJet_batch.cpp AnalysisJet.o TreeAnalyzer.o TreeAnalyzer_MultiJet.o QGLikelihoodCalculator.o $(ROOTFLAG) $(EXTRALIBS)

do2ndLevel_PhotonJetInclusive_batch: AnalysisJet.o AnalysisPhoton.o TreeAnalyzer.o TreeAnalyzer_PhotonJetInclusive.o fitTools.o do2ndLevel_PhotonJetInclusive_batch.cpp
	$(CC) -Wall -o do2ndLevel_PhotonJetInclusive_batch do2ndLevel_PhotonJetInclusive_batch.cpp AnalysisJet.o AnalysisPhoton.o TreeAnalyzer.o TreeAnalyzer_PhotonJetInclusive.o fitTools.o $(ROOTFLAG) -lMinuit


do2ndLevel_Calo_vs_PF_batch: AnalysisJet.o AnalysisPhoton.o TreeAnalyzer.o TreeAnalyzer_Calo_vs_PF.o do2ndLevel_Calo_vs_PF_batch.cpp fitTools.o
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o do2ndLevel_Calo_vs_PF_batch do2ndLevel_Calo_vs_PF_batch.cpp AnalysisJet.o AnalysisPhoton.o TreeAnalyzer.o TreeAnalyzer_Calo_vs_PF.o fitTools.o $(ROOTFLAG) -lMinuit

do2ndLevel_ZGamma:  TreeAnalyzer.o TreeAnalyzer_ZGamma.o
	$(CC) -Wall -o do2ndLevel_ZGamma do2ndLevel_ZGamma.cpp TreeAnalyzer.o TreeAnalyzer_ZGamma.o $(ROOTFLAG)


prova: AnalysisJet.o AnalysisPhoton.o TreeAnalyzer.o TreeAnalyzer_PhotonJet.o
	$(CC) -Wall -o prova prova_TreeAnalyzer_Photon.cpp AnalysisJet.o AnalysisPhoton.o TreeAnalyzer.o TreeAnalyzer_PhotonJet.o $(ROOTFLAG)


AnalysisJet.o: AnalysisJet.C
	$(CC) $(CFLAGS) AnalysisJet.C $(ROOTFLAG)

AnalysisPhoton.o: AnalysisPhoton.C
	$(CC) $(CFLAGS) AnalysisPhoton.C $(ROOTFLAG)

TreeAnalyzer.o: TreeAnalyzer.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools TreeAnalyzer.C $(ROOTFLAG)

TreeAnalyzer_PhotonJet.o: TreeAnalyzer_PhotonJet.C fitTools.o QGLikelihoodCalculator.o
	$(CC) $(CFLAGS) $(INCLUDES) fitTools.o QGLikelihoodCalculator.o TreeAnalyzer_PhotonJet.C $(ROOTFLAG)

TreeAnalyzer_VHgg.o: TreeAnalyzer_VHgg.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools TreeAnalyzer_VHgg.C $(ROOTFLAG)

TreeAnalyzer_MultiJet.o: TreeAnalyzer_MultiJet.C QGLikelihoodCalculator.o
	$(CC) $(CFLAGS) $(INCLUDES) QGLikelihoodCalculator.o TreeAnalyzer_MultiJet.C $(ROOTFLAG)

TreeAnalyzer_PhotonJetInclusive.o: TreeAnalyzer_PhotonJetInclusive.C fitTools.o
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools fitTools.o TreeAnalyzer_PhotonJetInclusive.C $(ROOTFLAG)

TreeAnalyzer_Calo_vs_PF.o: TreeAnalyzer_Calo_vs_PF.C fitTools.o
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools fitTools.o TreeAnalyzer_Calo_vs_PF.C $(ROOTFLAG)

TreeAnalyzer_ZGamma.o: TreeAnalyzer_ZGamma.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools TreeAnalyzer_ZGamma.C $(ROOTFLAG)

TreeAnalyzer_DiJet.o: TreeAnalyzer_DiJet.C
	$(CC) $(CFLAGS) TreeAnalyzer_DiJet.C $(ROOTFLAG)



TreeFinalizer.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/TreeFinalizer.C
	$(CC) $(CFLAGS) $(INCLUDES) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/TreeFinalizer.C $(ROOTFLAG)

TreeFinalizerC_MultiJet.o: TreeFinalizerC_MultiJet.C TreeFinalizer.o
	$(CC) $(CFLAGS) $(INCLUDES) TreeFinalizer.o TreeFinalizerC_MultiJet.C $(ROOTFLAG)

finalize_MultiJet: finalize_MultiJet.cpp TreeFinalizer.o TreeFinalizerC_MultiJet.o 
	$(CC)  $(INCLUDES) -Wall -o finalize_MultiJet TreeFinalizer.o TreeFinalizerC_MultiJet.o finalize_MultiJet.cpp `${ROOTSYS}/bin/root-config --cflags --libs`


make_omogeneizzato: make_omogeneizzato.cpp
	$(CC)  $(INCLUDES) -Wall -o make_omogeneizzato make_omogeneizzato.cpp `${ROOTSYS}/bin/root-config --cflags --libs`


drawPhotonJet_2bkg: fitTools.o DrawBase.o drawPhotonJet_2bkg.C
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o drawPhotonJet_2bkg drawPhotonJet_2bkg.C fitTools.o DrawBase.o $(ROOTFLAG) $(EXTRALIBS)

drawPF_vs_Calo: fitTools.o DrawBase.o drawPF_vs_Calo.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o drawPF_vs_Calo drawPF_vs_Calo.cpp fitTools.o DrawBase.o $(ROOTFLAG) $(EXTRALIBS)

drawEfficiencies_vs_QGCut: DrawBase.o fitTools.o drawEfficiencies_vs_QGCut.C
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -o drawEfficiencies_vs_QGCut drawEfficiencies_vs_QGCut.C DrawBase.o fitTools.o $(ROOTFLAG) -lMinuit

drawQGStudies: fitTools.o DrawBase.o drawQGStudies.C
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -o drawQGStudies drawQGStudies.C fitTools.o DrawBase.o $(ROOTFLAG) $(EXTRALIBS)

drawMultiJetQG: fitTools.o DrawBase.o drawMultiJetQG.C
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -o drawMultiJetQG drawMultiJetQG.C fitTools.o DrawBase.o $(ROOTFLAG) $(EXTRALIBS)

drawDiMultiJetQG: fitTools.o DrawBase.o drawDiMultiJetQG.C
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -o drawDiMultiJetQG drawDiMultiJetQG.C fitTools.o DrawBase.o $(ROOTFLAG) $(EXTRALIBS)

drawQGFits: fitTools.o DrawBase.o drawQGFits.cpp
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -o drawQGFits drawQGFits.cpp fitTools.o DrawBase.o $(ROOTFLAG) $(EXTRALIBS)

drawVHgg: fitTools.o DrawBase.o drawVHgg.C
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/ -o drawVHgg drawVHgg.C fitTools.o DrawBase.o $(ROOTFLAG) $(EXTRALIBS)

drawPhotonJet_ONLYDATA: fitTools.o DrawBase.o drawPhotonJet_ONLYDATA.C
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o drawPhotonJet_ONLYDATA drawPhotonJet_ONLYDATA.C fitTools.o DrawBase.o $(ROOTFLAG) -lMinuit

drawPhotonJetExtrap: fitTools.o DrawBase.o DrawExtrap.o drawPhotonJetExtrap.C
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o drawPhotonJetExtrap drawPhotonJetExtrap.C fitTools.o DrawBase.o DrawExtrap.o $(ROOTFLAG) $(EXTRALIBS)

draw_all_methods_vs_pt: fitTools.o DrawBase.o draw_all_methods_vs_pt.C
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o draw_all_methods_vs_pt draw_all_methods_vs_pt.C fitTools.o DrawBase.o $(ROOTFLAG) $(EXTRALIBS)

draw_PhotonJet_vs_QCD_vs_pt: fitTools.o DrawBase.o draw_PhotonJet_vs_QCD_vs_pt.C
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o draw_PhotonJet_vs_QCD_vs_pt draw_PhotonJet_vs_QCD_vs_pt.C fitTools.o DrawBase.o $(ROOTFLAG) -lMinuit

drawSystPhotonJetExtrap: fitTools.o drawSystPhotonJetExtrap.C
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o drawSystPhotonJetExtrap drawSystPhotonJetExtrap.C fitTools.o $(ROOTFLAG) -lMinuit

draw_dijet_gammajet: fitTools.o DrawBase.o draw_dijet_gammajet.C
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o draw_dijet_gammajet draw_dijet_gammajet.C fitTools.o DrawBase.o $(ROOTFLAG) -lMinuit

drawResoBreakdown: fitTools.o DrawBase.o drawResoBreakdown.C
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o drawResoBreakdown drawResoBreakdown.C fitTools.o DrawBase.o $(ROOTFLAG) -lMinuit

drawMinBias: fitTools.o DrawBase.o
	$(CC) -Wall -o drawMinBias drawMinBias.cpp fitTools.o DrawBase.o $(ROOTFLAG) -lMinuit

drawZGamma: fitTools.o DrawBase.o drawZGamma.cpp
	$(CC) -Wall $(INCLUDES) -o drawZGamma drawZGamma.cpp fitTools.o DrawBase.o $(ROOTFLAG) $(EXTRALIBS)


QGLikelihoodCalculator.o: $(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood/QGLikelihoodCalculator.C 
	$(CC) $(CFLAGS) $(INCLUDES) $(CMSSW_BASE)/src/UserCode/pandolf/QGLikelihood/QGLikelihoodCalculator.C $(ROOTFLAG) $(EXTRALIBS)

DrawBase.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/DrawBase.C fitTools.o
	$(CC) $(CFLAGS) $(INCLUDES) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/DrawBase.C fitTools.o $(ROOTFLAG) $(EXTRALIBS)

fitTools.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/fitTools.C
	$(CC) $(CFLAGS) $(INCLUDES) $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/fitTools.C $(ROOTFLAG) $(EXTRALIBS)

DrawExtrap.o: DrawExtrap.C DrawBase.o
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools DrawExtrap.C DrawBase.o $(ROOTFLAG)

merge_and_setWeights_PhotonJet: merge_and_setWeights_PhotonJet.cpp
	$(CC) -Wall -o merge_and_setWeights_PhotonJet merge_and_setWeights_PhotonJet.cpp $(ROOTFLAG) $(EXTRALIBS)

merge_and_setWeights_Calo_vs_PF: merge_and_setWeights_Calo_vs_PF.cpp
	$(CC) -Wall -o merge_and_setWeights_Calo_vs_PF merge_and_setWeights_Calo_vs_PF.cpp $(ROOTFLAG)


Utilities.o: CondFormats/JetMETObjects/src/Utilities.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/JetMETCorrections/GammaJet/test/ CondFormats/JetMETObjects/src/Utilities.cc $(ROOTFLAG)

JetCorrectorParameters.o: CondFormats/JetMETObjects/src/JetCorrectorParameters.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/JetMETCorrections/GammaJet/test/ CondFormats/JetMETObjects/src/JetCorrectorParameters.cc $(ROOTFLAG)

FactorizedJetCorrector.o: CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/JetMETCorrections/GammaJet/test/ CondFormats/JetMETObjects/src/FactorizedJetCorrector.cc $(ROOTFLAG)

SimpleJetCorrector.o: CondFormats/JetMETObjects/src/SimpleJetCorrector.cc
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/JetMETCorrections/GammaJet/test/ CondFormats/JetMETObjects/src/SimpleJetCorrector.cc $(ROOTFLAG)


clean:
	rm *.o
