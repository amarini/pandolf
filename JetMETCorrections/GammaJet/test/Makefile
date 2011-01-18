CC = g++
CFLAGS = -Wall -c
ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs`


do2ndLevel_DiJet: AnalysisJet.o TreeAnalyzer.o TreeAnalyzer_DiJet.o
	$(CC) -Wall -o do2ndLevel_DiJet do2ndLevel_DiJet.cpp AnalysisJet.o TreeAnalyzer.o TreeAnalyzer_DiJet.o $(ROOTFLAG)

do2ndLevel_PhotonJet: AnalysisJet.o AnalysisPhoton.o TreeAnalyzer.o TreeAnalyzer_PhotonJet.o
	$(CC) -Wall -o do2ndLevel_PhotonJet do2ndLevel_PhotonJet.cpp AnalysisJet.o AnalysisPhoton.o TreeAnalyzer.o TreeAnalyzer_PhotonJet.o $(ROOTFLAG) -lMinuit

do2ndLevel_ZGamma:  TreeAnalyzer.o TreeAnalyzer_ZGamma.o
	$(CC) -Wall -o do2ndLevel_ZGamma do2ndLevel_ZGamma.cpp TreeAnalyzer.o TreeAnalyzer_ZGamma.o $(ROOTFLAG)

do2ndLevel_Calo_vs_PF:  TreeAnalyzer_Calo_vs_PF.o
	$(CC) -Wall -o do2ndLevel_Calo_vs_PF do2ndLevel_Calo_vs_PF.cpp TreeAnalyzer_Calo_vs_PF.o $(ROOTFLAG)

prova: AnalysisJet.o AnalysisPhoton.o TreeAnalyzer.o TreeAnalyzer_PhotonJet.o
	$(CC) -Wall -o prova prova_TreeAnalyzer_Photon.cpp AnalysisJet.o AnalysisPhoton.o TreeAnalyzer.o TreeAnalyzer_PhotonJet.o $(ROOTFLAG)

AnalysisJet.o: AnalysisJet.C
	$(CC) $(CFLAGS) AnalysisJet.C $(ROOTFLAG)

AnalysisPhoton.o: AnalysisPhoton.C
	$(CC) $(CFLAGS) AnalysisPhoton.C $(ROOTFLAG)

TreeAnalyzer.o: TreeAnalyzer.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools TreeAnalyzer.C $(ROOTFLAG)

TreeAnalyzer_PhotonJet.o: TreeAnalyzer_PhotonJet.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools TreeAnalyzer_PhotonJet.C $(ROOTFLAG)

TreeAnalyzer_Calo_vs_PF.o: TreeAnalyzer_Calo_vs_PF.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools TreeAnalyzer_Calo_vs_PF.C $(ROOTFLAG)

TreeAnalyzer_ZGamma.o: TreeAnalyzer_ZGamma.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools TreeAnalyzer_ZGamma.C $(ROOTFLAG)

TreeAnalyzer_DiJet.o: TreeAnalyzer_DiJet.C
	$(CC) $(CFLAGS) TreeAnalyzer_DiJet.C $(ROOTFLAG)

drawPhotonJet: fitTools.o DrawBase.o
	$(CC) -Wall -o drawPhotonJet drawPhotonJet.C fitTools.o DrawBase.o $(ROOTFLAG) -lMinuit

drawPhotonJet_2bkg: fitTools.o DrawBase.o drawPhotonJet_2bkg.C
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o drawPhotonJet_2bkg drawPhotonJet_2bkg.C fitTools.o DrawBase.o $(ROOTFLAG) -lMinuit

drawPhotonJet_ONLYDATA: fitTools.o DrawBase.o drawPhotonJet_ONLYDATA.C
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o drawPhotonJet_ONLYDATA drawPhotonJet_ONLYDATA.C fitTools.o DrawBase.o $(ROOTFLAG) -lMinuit

draw_Calo_vs_PF: fitTools.o DrawBase.o draw_Calo_vs_PF.C
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o draw_Calo_vs_PF draw_Calo_vs_PF.C fitTools.o DrawBase.o $(ROOTFLAG) -lMinuit

drawPhotonJetExtrap: fitTools.o DrawBase.o DrawExtrap.o drawPhotonJetExtrap.C
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o drawPhotonJetExtrap drawPhotonJetExtrap.C fitTools.o DrawBase.o DrawExtrap.o $(ROOTFLAG) -lMinuit

draw_all_methods_vs_pt: fitTools.o DrawBase.o draw_all_methods_vs_pt.C
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o draw_all_methods_vs_pt draw_all_methods_vs_pt.C fitTools.o DrawBase.o $(ROOTFLAG) -lMinuit

draw_PhotonJet_vs_QCD_vs_pt: fitTools.o DrawBase.o draw_PhotonJet_vs_QCD_vs_pt.C
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o draw_PhotonJet_vs_QCD_vs_pt draw_PhotonJet_vs_QCD_vs_pt.C fitTools.o DrawBase.o $(ROOTFLAG) -lMinuit

drawSystPhotonJetExtrap: fitTools.o drawSystPhotonJetExtrap.C
	$(CC) -Wall -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools -o drawSystPhotonJetExtrap drawSystPhotonJetExtrap.C fitTools.o $(ROOTFLAG) -lMinuit

drawMinBias: fitTools.o DrawBase.o
	$(CC) -Wall -o drawMinBias drawMinBias.cpp fitTools.o DrawBase.o $(ROOTFLAG) -lMinuit

DrawBase.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/DrawBase.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/DrawBase.C $(ROOTFLAG)

fitTools.o: $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/fitTools.C
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools $(CMSSW_BASE)/src/UserCode/pandolf/CommonTools/fitTools.C $(ROOTFLAG) -lMinuit

DrawExtrap.o: DrawExtrap.C DrawBase.o
	$(CC) $(CFLAGS) -I$(CMSSW_BASE)/src/UserCode/pandolf/CommonTools DrawExtrap.C DrawBase.o $(ROOTFLAG)

merge_and_setWeights_PhotonJet: merge_and_setWeights_PhotonJet.cpp
	$(CC) -Wall -o merge_and_setWeights_PhotonJet merge_and_setWeights_PhotonJet.cpp $(ROOTFLAG)

merge_and_setWeights_Calo_vs_PF: merge_and_setWeights_Calo_vs_PF.cpp
	$(CC) -Wall -o merge_and_setWeights_Calo_vs_PF merge_and_setWeights_Calo_vs_PF.cpp $(ROOTFLAG)


clean:
	rm *.o
