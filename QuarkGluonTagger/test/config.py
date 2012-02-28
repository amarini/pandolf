import FWCore.ParameterSet.Config as cms
import RecoJets.JetProducers.JetIDParams_cfi
theJetIDParams = RecoJets.JetProducers.JetIDParams_cfi.JetIDParams.clone()
from RecoJets.Configuration.RecoGenJets_cff import *
from RecoMET.Configuration.RecoGenMET_cff import *
from RecoMET.Configuration.GenMETParticles_cff import *

process = cms.Process("myprocess")

process.load("FWCore.MessageLogger.MessageLogger_cfi")


process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("PhysicsTools.HepMCCandAlgos.genParticleCandidates_cfi")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.TrackerGeometryBuilder.trackerGeometry_cfi")
process.load("Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load('Configuration/StandardSequences/Reconstruction_cff')


process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START52_V2A::All') 
#process.GlobalTag.globaltag = cms.string('START42_V16::All') 


#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck")

process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(
    'file:/tmp/pandolf/EC9EFACC-C34C-E111-B207-003048FFCBA8.root'
)

)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.options = cms.untracked.PSet(
#    SkipEvent = cms.untracked.vstring('ProductNotFound')
    wantSummary = cms.untracked.bool(True)
)


#############   Include the corrections ##########
process.load("JetMETCorrections.Configuration.DefaultJEC_cff")
process.load("JetMETCorrections.Type1MET.MetType1Corrections_cff")
process.load('RecoJets.Configuration.RecoPFJets_cff')


process.ak5PFJets.doAreaFastjet = True


from RecoJets.JetProducers.kt4PFJets_cfi import *

process.kt6PFJets = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )

process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)



##--------- value map for qg likelihood variable -----
process.qglAK5PF   = cms.EDProducer("QuarkGluonTagger",
          jets     = cms.InputTag("ak5PFJets"),
          rho      = cms.InputTag('kt6PFJetsForIsolation','rho'),
          jec      = cms.string('ak5PFL1FastL2L3'),
          #jec      = cms.string('ak5PFL1FastL2L3Residual'),
          filename = cms.string('/afs/cern.ch/user/p/pandolf/scratch1/CMSSW_4_2_8_patch7/src/pandolf/QGLikelihood/QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2.root')
)



process.p = cms.Path( process.ak5PFJets * process.kt6PFJets * process.kt6PFJetsForIsolation * process.qglAK5PF )
#process.p = cms.Path( process.ak5PFJets * process.kt6PFJets * process.ak5PFJetsL1FastL2L3 * process.kt6PFJetsForIsolation * process.qglAK5PF )
