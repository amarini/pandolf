# ---------------------------------------------------------------------------------------
# USAGE: source finalize_TTZTrilepton_all.sh [selType] [bTaggerType] [PUType] [leptType]
# ---------------------------------------------------------------------------------------



./finalize_TTZTrilepton DYJetsToLL_TuneZ2_M-50_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1 $1 $2 $3
./finalize_TTZTrilepton TTJets_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1 $1 $2 $3
./finalize_TTZTrilepton TTW_TuneZ2_7TeV-madgraphCMSSW42xPUv2_spadhi $1 $2 $3
./finalize_TTZTrilepton TTZ_TuneZ2_7TeV-madgraphCMSSW42xPUv3_spadhi $1 $2 $3
./finalize_TTZTrilepton WZJetsTo3LNu_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1 $1 $2 $3
