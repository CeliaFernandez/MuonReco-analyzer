import FWCore.ParameterSet.Config as cms


standAloneMuonsInMuons = cms.EDFilter("CandViewSelector",
                                      src = cms.InputTag("muons"),
                                      cut = cms.string( "isStandAloneMuon" )
                                     )

standAlonePairFilter = cms.EDFilter("CandCountFilter",
                           src = cms.InputTag("standAloneMuonsInMuons"),
                           minNumber = cms.uint32(2),
                          )
