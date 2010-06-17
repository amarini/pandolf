CC = g++
CFLAGS = -Wall -c
ROOTFLAG = `${ROOTSYS}/bin/root-config --cflags --libs`


do2ndLevel_Inclusive: AnalysisJet.o TreeAnalyzer.o TreeAnalyzer_Inclusive.o
	$(CC) -Wall -o do2ndLevel_Inclusive do2ndLevel_Inclusive.cpp AnalysisJet.o TreeAnalyzer.o TreeAnalyzer_Inclusive.o $(ROOTFLAG)

do2ndLevel_DiJet: AnalysisJet.o TreeAnalyzer.o TreeAnalyzer_DiJet.o
	$(CC) -Wall -o do2ndLevel_DiJet do2ndLevel_DiJet.cpp AnalysisJet.o TreeAnalyzer.o TreeAnalyzer_DiJet.o $(ROOTFLAG)

do2ndLevel_PhotonJet: AnalysisJet.o AnalysisPhoton.o TreeAnalyzer.o TreeAnalyzer_PhotonJet.o
	$(CC) -Wall -o do2ndLevel_PhotonJet do2ndLevel_PhotonJet.cpp AnalysisJet.o AnalysisPhoton.o TreeAnalyzer.o TreeAnalyzer_PhotonJet.o $(ROOTFLAG)

do2ndLevel_HZZlljj: AnalysisJet.o TreeAnalyzer.o TreeAnalyzer_HZZlljj.o
	$(CC) -Wall -o do2ndLevel_HZZlljj do2ndLevel_HZZlljj.cpp AnalysisJet.o TreeAnalyzer.o TreeAnalyzer_HZZlljj.o $(ROOTFLAG)

prova: AnalysisJet.o AnalysisPhoton.o TreeAnalyzer.o TreeAnalyzer_PhotonJet.o
	$(CC) -Wall -o prova prova_TreeAnalyzer_Photon.cpp AnalysisJet.o AnalysisPhoton.o TreeAnalyzer.o TreeAnalyzer_PhotonJet.o $(ROOTFLAG)

AnalysisJet.o: AnalysisJet.C
	$(CC) $(CFLAGS) AnalysisJet.C $(ROOTFLAG)

AnalysisPhoton.o: AnalysisPhoton.C
	$(CC) $(CFLAGS) AnalysisPhoton.C $(ROOTFLAG)

TreeAnalyzer.o: TreeAnalyzer.C
	$(CC) $(CFLAGS) TreeAnalyzer.C $(ROOTFLAG)

TreeAnalyzer_PhotonJet.o: TreeAnalyzer_PhotonJet.C
	$(CC) $(CFLAGS) TreeAnalyzer_PhotonJet.C $(ROOTFLAG)

TreeAnalyzer_HZZlljj.o: TreeAnalyzer_HZZlljj.C
	$(CC) $(CFLAGS) TreeAnalyzer_HZZlljj.C $(ROOTFLAG)

TreeAnalyzer_DiJet.o: TreeAnalyzer_DiJet.C
	$(CC) $(CFLAGS) TreeAnalyzer_DiJet.C $(ROOTFLAG)

TreeAnalyzer_Inclusive.o: TreeAnalyzer_Inclusive.C
	$(CC) $(CFLAGS) TreeAnalyzer_Inclusive.C $(ROOTFLAG)


clean:
	rm TreeAnalyzer.o
	rm AnalysisPhoton.o
	rm AnalysisJet.o
	rm TreeAnalyzer_PhotonJet.o
