import FWCore.ParameterSet.Config as cms

displacedMuonsRECOAnalysis = cms.EDAnalyzer('displacedMuonsRECOAnalyzer',
    nameOfOutput = cms.string('output.root'),
    EventInfo = cms.InputTag("generator"),
    RunInfo = cms.InputTag("generator"),
    BeamSpot = cms.InputTag("offlineBeamSpot"),
    muonCollection = cms.InputTag("displacedMuons")
)


