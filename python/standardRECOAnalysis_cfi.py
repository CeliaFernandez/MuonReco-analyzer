import FWCore.ParameterSet.Config as cms

recoAnalysis = cms.EDAnalyzer('standardRECOAnalyzer',
    nameOfOutput = cms.string('output.root'),
    EventInfo = cms.InputTag("generator"),
    RunInfo = cms.InputTag("generator"),
    BeamSpot = cms.InputTag("offlineBeamSpot"),
    muonCollection = cms.InputTag("muons"),
    GenCollection = cms.InputTag("genParticles"),
    theGenEventInfoProduct = cms.InputTag("generator"),
)


