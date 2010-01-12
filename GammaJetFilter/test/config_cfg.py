import FWCore.ParameterSet.Config as cms

process = cms.Process("photonSkim")

process.load("FWCore.MessageLogger.MessageLogger_cfi")

#process.MessageLogger = cms.Service("MessageLogger",
#    cout = cms.untracked.PSet(
#        threshold = cms.untracked.string('WARNING'),
#        noLineBreaks = cms.untracked.bool(True),
#        noTimeStamps = cms.untracked.bool(True),
#        default = cms.untracked.PSet(
#            limit = cms.untracked.int32(0)
#        ),
#        EcalPositionFromTrack = cms.untracked.PSet(
#            limit = cms.untracked.int32(0)
#        )
#    ),
#    categories = cms.untracked.vstring('EcalPositionFromTrack'),
#    destinations = cms.untracked.vstring('cout')
#)


process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.load("PhysicsTools.HepMCCandAlgos.genParticleCandidates_cfi")

process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")

process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")

#process.load("MagneticField.Engine.uniformMagneticField_cfi")

process.load("Geometry.TrackerGeometryBuilder.trackerGeometry_cfi")

process.load("Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi")

process.load("Geometry.CaloEventSetup.CaloTopology_cfi")


process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck")

process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),
# fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/cms/store/relval/CMSSW_3_1_0_pre11/RelValZmumuJets_Pt_20_300_GEN/GEN-SIM-RECO/MC_31X_V1_LowLumiPileUp-v1/0001/FE549F7D-2C65-DE11-9B82-001D09F2AF1E.root')
#    fileNames = cms.untracked.vstring('file:/cmsrm/pc17/pandolf/eventi_PhotonJetPt170_Summer09_10TeV.root')
#    fileNames = cms.untracked.vstring('file:/cmsrm/pc17/delre/7A62C541-37A0-DE11-BF0C-00215E221098.root')
# fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/cms/store/relval/CMSSW_3_1_0_pre11/RelValZmumuJets_Pt_20_300_GEN/GEN-SIM-RECO/MC_31X_V1_LowLumiPileUp-v1/0001/FE549F7D-2C65-DE11-9B82-001D09F2AF1E.root')
#    fileNames = cms.untracked.vstring('file:/tmp/voutila/Cern/data/summer09/raw/PhotonJet_Pt80to120_Summer09-MC_31X_V3-v1_x100.root')
#    fileNames = cms.untracked.vstring('file:/cmsrm/pc18/pandolf/CMSSW_3_3_4/src/JetMETCorrections/GammaJet/test/bit40or41skim.root')
#    fileNames = cms.untracked.vstring('file:Run124120_display_Dec14.root')
#    fileNames = cms.untracked.vstring('file:/cmsrm/pc18/pandolf/data/DiJetFilter5/DiJetSkim_124120.root')
     fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/cms/store/caf/user/meridian/MinimumBias/BeamCommissioning09_BSCFilter_v6/c32fbdf401a82bf088a335dc690defa6/bscFilter_124030_1.root',
'rfio:/castor/cern.ch/cms/store/caf/user/meridian/MinimumBias/BeamCommissioning09_BSCFilter_v6/c32fbdf401a82bf088a335dc690defa6/bscFilter_124030_10.root',
'rfio:/castor/cern.ch/cms/store/caf/user/meridian/MinimumBias/BeamCommissioning09_BSCFilter_v6/c32fbdf401a82bf088a335dc690defa6/bscFilter_124030_11.root',
'rfio:/castor/cern.ch/cms/store/caf/user/meridian/MinimumBias/BeamCommissioning09_BSCFilter_v6/c32fbdf401a82bf088a335dc690defa6/bscFilter_124030_12.root',
'rfio:/castor/cern.ch/cms/store/caf/user/meridian/MinimumBias/BeamCommissioning09_BSCFilter_v6/c32fbdf401a82bf088a335dc690defa6/bscFilter_124030_13.root',
'rfio:/castor/cern.ch/cms/store/caf/user/meridian/MinimumBias/BeamCommissioning09_BSCFilter_v6/c32fbdf401a82bf088a335dc690defa6/bscFilter_124030_14.root',
'rfio:/castor/cern.ch/cms/store/caf/user/meridian/MinimumBias/BeamCommissioning09_BSCFilter_v6/c32fbdf401a82bf088a335dc690defa6/bscFilter_124030_15.root',
'rfio:/castor/cern.ch/cms/store/caf/user/meridian/MinimumBias/BeamCommissioning09_BSCFilter_v6/c32fbdf401a82bf088a335dc690defa6/bscFilter_124030_16.root',
'rfio:/castor/cern.ch/cms/store/caf/user/meridian/MinimumBias/BeamCommissioning09_BSCFilter_v6/c32fbdf401a82bf088a335dc690defa6/bscFilter_124030_17.root',
'rfio:/castor/cern.ch/cms/store/caf/user/meridian/MinimumBias/BeamCommissioning09_BSCFilter_v6/c32fbdf401a82bf088a335dc690defa6/bscFilter_124030_18.root',
'rfio:/castor/cern.ch/cms/store/caf/user/meridian/MinimumBias/BeamCommissioning09_BSCFilter_v6/c32fbdf401a82bf088a335dc690defa6/bscFilter_124030_2.root',
'rfio:/castor/cern.ch/cms/store/caf/user/meridian/MinimumBias/BeamCommissioning09_BSCFilter_v6/c32fbdf401a82bf088a335dc690defa6/bscFilter_124030_3.root',
'rfio:/castor/cern.ch/cms/store/caf/user/meridian/MinimumBias/BeamCommissioning09_BSCFilter_v6/c32fbdf401a82bf088a335dc690defa6/bscFilter_124030_4.root',
'rfio:/castor/cern.ch/cms/store/caf/user/meridian/MinimumBias/BeamCommissioning09_BSCFilter_v6/c32fbdf401a82bf088a335dc690defa6/bscFilter_124030_5.root',
'rfio:/castor/cern.ch/cms/store/caf/user/meridian/MinimumBias/BeamCommissioning09_BSCFilter_v6/c32fbdf401a82bf088a335dc690defa6/bscFilter_124030_6.root',
'rfio:/castor/cern.ch/cms/store/caf/user/meridian/MinimumBias/BeamCommissioning09_BSCFilter_v6/c32fbdf401a82bf088a335dc690defa6/bscFilter_124030_7.root',
'rfio:/castor/cern.ch/cms/store/caf/user/meridian/MinimumBias/BeamCommissioning09_BSCFilter_v6/c32fbdf401a82bf088a335dc690defa6/bscFilter_124030_8.root',
'rfio:/castor/cern.ch/cms/store/caf/user/meridian/MinimumBias/BeamCommissioning09_BSCFilter_v6/c32fbdf401a82bf088a335dc690defa6/bscFilter_124030_9.root'
)

)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.photonJetSkim = cms.EDFilter("GammaJetFilter",
    Photonsrc = cms.untracked.InputTag("photons"),
    jetspfakt5 = cms.untracked.InputTag("ak5PFJets"),
    jetspfakt7 = cms.untracked.InputTag("ak7PFJets")

)

#### the path
process.p = cms.Path(process.photonJetSkim)

#### output
process.outputSkim = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *','drop *_MEtoEDMConverter_*_*'),
    fileName = cms.untracked.string("photonJetSkim.root"),
    dataset = cms.untracked.PSet(
      dataTier = cms.untracked.string('RAW-RECO'),
      filterName = cms.untracked.string('photonJetSkim')
    ),

    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p'))
)

process.outpath = cms.EndPath(process.outputSkim)
