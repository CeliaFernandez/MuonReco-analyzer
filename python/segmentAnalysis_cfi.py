import FWCore.ParameterSet.Config as cms

segmentAnalysis = cms.EDAnalyzer('segmentAnalyzer',
    nameOfOutput = cms.string('output.root'),
    EventInfo = cms.InputTag("generator"),
    RunInfo = cms.InputTag("generator"),

    ## Segments
    DTRecSegment4DCollectionLabel = cms.InputTag("dt4DSegments"),
    CSCSegmentCollectionLabel     = cms.InputTag("cscSegments"),
    GEMSegmentCollectionLabel     = cms.InputTag("gemSegments"),
#    ME0SegmentCollectionLabel     = cms.InputTag("me0Segments"),

    ## Hits
    RPCHitCollectionLabel         = cms.InputTag("rpcRecHits"),
    GEMHitCollectionLabel         = cms.InputTag("gemRecHits"),

    BeamSpot = cms.InputTag("offlineBeamSpot")
)

