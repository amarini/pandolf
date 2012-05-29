# ---------------------------------------------------------------------------------------
# USAGE: source finalize_TTZTrilepton_all.sh [selType] [bTaggerType] 
# ---------------------------------------------------------------------------------------



#./finalize_TTZTrilepton DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1 $1 $2 $3
./finalize_TTZTrilepton DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Fall11 $1 $2
./finalize_TTZTrilepton WJetsToLNu_TuneZ2_7TeV-madgraph-tauola_Summer11 $1 $2
./finalize_TTZTrilepton TTJ_Fall11_highstat $1 $2
./finalize_TTZTrilepton TTW_TuneZ2_7TeV-madgraphCMSSW42xPUv2_spadhi $1 $2
./finalize_TTZTrilepton TTZ_TuneZ2_7TeV-madgraphCMSSW42xPUv3_spadhi $1 $2
./finalize_TTZTrilepton VV_Summer11 $1 $2
./finalize_TTZTrilepton GVJets_7TeV-madgraph_Fall11 $1 $2

./finalize_TTZTrilepton DATA_Run2011_FULL $1 $2
