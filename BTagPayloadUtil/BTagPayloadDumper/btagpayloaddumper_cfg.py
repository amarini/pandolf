import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("EmptySource")

#process.source = cms.Source("PoolSource",
#    # replace 'myfile.root' with the source file you want to use
#    fileNames = cms.untracked.vstring(
#        'file:events_relval_TTbar_423.root'
#    )
#)

process.load("CondCore.DBCommon.CondDBCommon_cfi")


#Data measurements from Summer11
process.load ("RecoBTag.PerformanceDB.BTagPerformanceDB1107")
process.load ("RecoBTag.PerformanceDB.PoolBTagPerformanceDB1107")


process.myanalysis = cms.EDAnalyzer('BTagPayloadDumper'
)


process.p = cms.Path(process.myanalysis)
