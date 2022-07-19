import FWCore.ParameterSet.Config as cms

standaloneDuplicates = cms.EDAnalyzer('standaloneDuplicates',
    nameOfOutput = cms.string('output.root'),
    EventInfo = cms.InputTag("generator"),
    RunInfo = cms.InputTag("generator"),
    BeamSpot = cms.InputTag("offlineBeamSpot"),
    muonCollection = cms.InputTag("muons")
)


