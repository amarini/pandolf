# usage: source finalize_all.sh [redntpProdVersion] [selectionType]


 cp finalize_VHgg finalize_VHgg_tmp
 ./finalize_VHgg_tmp  $1 HToGG_M-125_8TeV-pythia6 $2
 ./finalize_VHgg_tmp  $1 DiPhoton_8TeV-pythia6 $2
 ./finalize_VHgg_tmp  $1 V_8TeV $2
 ./finalize_VHgg_tmp  $1 GJet_doubleEMEnriched_TuneZ2star_8TeV-pythia6 $2
 ./finalize_VHgg_tmp  $1 QCD_doubleEMEnriched_TuneZ2star_8TeV-pythia6 $2
 ./finalize_VHgg_tmp  $1 TT_8TeV $2
 ./finalize_VHgg_tmp  $1 VV_8TeV $2
 ./finalize_VHgg_tmp  $1 VGG_8TeV $2
 rm finalize_VHgg_tmp
