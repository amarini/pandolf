# $Id: common_dump_config.py,v 1.7 2011/11/20 21:45:20 meridian Exp $
#
#  common configuration to dump ntuples in MC and data
#    all changes affecting the path and additional modules msut be done here
#

import FWCore.ParameterSet.Config as cms
import RecoJets.JetProducers.JetIDParams_cfi
theJetIDParams = RecoJets.JetProducers.JetIDParams_cfi.JetIDParams.clone()
from RecoJets.Configuration.RecoGenJets_cff import *

process = cms.Process("myprocess")

process.load("FWCore.MessageLogger.MessageLogger_cfi")

is41X=False

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("PhysicsTools.HepMCCandAlgos.genParticleCandidates_cfi")
process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.load("Geometry.CaloEventSetup.CaloGeometry_cfi")
process.load("Geometry.TrackerGeometryBuilder.trackerGeometry_cfi")
process.load("Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi")
process.load("Geometry.CaloEventSetup.CaloTopology_cfi")
process.load('Configuration/StandardSequences/Reconstruction_cff')


process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),
    fileNames = cms.untracked.vstring(
    #'file:/tmp/pandolf/events_GluGluToHToZZTo2L2Q_M-550_7TeV-powheg-pythia6_Summer11_PROVA.root'
    'file:/afs/cern.ch/work/p/pandolf/public/DoubleElectron_RunA_53X_Jul13rereco.root'
)

)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

process.options = cms.untracked.PSet(
#    SkipEvent = cms.untracked.vstring('ProductNotFound')
    wantSummary = cms.untracked.bool(True)
)

process.MessageLogger.cerr.FwkReport.reportEvery = 1




process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('START52_V2A::All')
#process.GlobalTag.globaltag = cms.string('START42_V16::All')



#############   Include the corrections ##########
process.load("JetMETCorrections.Configuration.DefaultJEC_cff")

process.ak5PFJets.doAreaFastjet = True



process.myanalysis = cms.EDAnalyzer("ProvaQGAnalyzer",
    debug = cms.bool(False),
#    PFParameters = PhotonFixParameters,
    outFileName = cms.untracked.string("output.root"),                                    
    JetCorrectionService_pfakt5 = cms.string('ak5PFL1FastL2L3'),
    jetspfakt5 = cms.untracked.InputTag("ak5PFJets"),
    pfjetptthr = cms.double(4.),
    pfjetnmin = cms.int32(10),
    JetIDParams = theJetIDParams,
    Xsec = cms.double(1.)
                                    
)

# compute rho with PF candidates
process.load('RecoJets.Configuration.RecoPFJets_cff')
#process.kt6PFJets = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True,  doAreaFastjet = cms.bool(True), voronoiRfact = cms.double(0.9) )
process.kt6PFJets = process.kt6PFJets.clone( rParam = 0.6, doRhoFastjet = True )

process.kt6PFJetsForIso = process.kt6PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIso.Rho_EtaMax = cms.double(2.5)




##--------- value map for qg likelihood variable -----
process.qglAK5PF   = cms.EDProducer("QuarkGluonTagger2012",
          jets     = cms.InputTag("ak5PFJets"),
          rho      = cms.InputTag('kt6PFJetsForIso','rho'),
          jec      = cms.string('ak5PFL1FastL2L3'),
          #jec      = cms.string('ak5PFL1FastL2L3Residual'),
          #filename = cms.string('/afs/cern.ch/user/p/pandolf/scratch1/CMSSW_4_2_8_patch7/src/pandolf/QGLikelihood/QG_QCD_Pt-15to3000_TuneZ2_Flat_7TeV_pythia6_Summer11-PU_S3_START42_V11-v2.root')
)



process.analysisSequence = cms.Sequence(  )
    
process.analysisSequence *=  (process.kt6PFJets * process.ak5PFJets *process.kt6PFJetsForIso * process.ak5PFJetsL1FastL2L3 * process.qglAK5PF * process.myanalysis)

process.p = cms.Path(process.analysisSequence)

