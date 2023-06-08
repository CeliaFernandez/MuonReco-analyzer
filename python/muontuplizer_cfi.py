import FWCore.ParameterSet.Config as cms

ntuples = cms.EDAnalyzer('recoMuonNtuplizer',
    nameOfOutput = cms.string('output.root'),
    isData = cms.bool(True),
    EventInfo = cms.InputTag("generator"),
    RunInfo = cms.InputTag("generator"),
    BeamSpot = cms.InputTag("offlineBeamSpot"),
    muonCollection = cms.InputTag("muons"),
    trackCollection = cms.InputTag("generalTracks")
)


