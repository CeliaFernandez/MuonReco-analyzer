import FWCore.ParameterSet.Config as cms

displacedMuonsPATAnalysis = cms.EDAnalyzer('displacedMuonsPATAnalyzer',
    nameOfOutput = cms.string('output.root'),
    EventInfo = cms.InputTag("generator"),
    RunInfo = cms.InputTag("generator"),
    BeamSpot = cms.InputTag("offlineBeamSpot"),
    muonCollection = cms.InputTag("slimmedDisplacedMuons")
)


