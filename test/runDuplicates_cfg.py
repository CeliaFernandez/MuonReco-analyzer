import FWCore.ParameterSet.Config as cms
import os

# Run standaloneDuplicates plugin
#
# @ Celia Fernandez Madrazo 


process = cms.Process("demo")
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.Services_cff')

# Debug printout & summaries.
process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport.limit = cms.untracked.int32(50)
#process.MessageLogger.cerr.default.limit = cms.untracked.int32(50)

process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True),
  # Set up multi-threaded run. Must be consistent with config.JobType.numCores in crab_cfg.py.
  #numberOfThreads=cms.untracked.uint32(8)
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

# Select number of events to be processed
nEvents = -1
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(nEvents) )

# Read events
isData = True # Always true by default (running on MC is useless)
inputdir = "/eos/user/f/fernance/DYJetsToMuMu_M-50_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/skimmed-duplicatedmuons_0eta/220601_141623/"
inputdir = "/eos/user/f/fernance/DYJetsToMuMu_M-50_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/crab_pickEvents/220718_103333/0000/"
listOfFiles = ['file:' + os.path.join(inputdir, f) for f in os.listdir(inputdir) if '.root' in f]

#listOfFiles = ['file:/afs/cern.ch/work/f/fernance/private/MuonPOG/L3-RECO/MuonReco-analysis/DuplicatedMuons-study/CMSSW_10_6_20/src/RecoMuon/MuonIdentification/test/redo_muons.root']

if (isData):
  process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( listOfFiles ),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents = cms.untracked.uint32(0)
  )
  process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mcRun2_asymptotic_v13')

## Re-do muons 
## Clone the muons1stStep process and exclude generalTracks as info saved in AOD is not
## enough to re-run the Tracker Muon reconstruction
##
## (links collection is made from existing 'muons' collection)
## IMPORTANT: For this setup to run, timing computation must be dropped from MuonIdProducer

process.load("RecoMuon.MuonIdentification.links_cfi")
process.load("RecoMuon.MuonIdentification.muonIdProducerSequence_cff")
process.load("RecoMuon.Configuration.RecoMuonPPonly_cff")
process.load("RecoMuon.Configuration.RecoMuon_cff")

process.fixedMuons1stStep = process.muons1stStep.clone()
process.fixedMuons1stStep.inputCollectionLabels = cms.VInputTag(cms.InputTag("globalMuonLinks"),
                                                           cms.InputTag("standAloneMuons","UpdatedAtVtx"))
process.fixedMuons1stStep.inputCollectionTypes = cms.vstring('links',
                                                        'outer tracks')

process.fixedMuons1stStep.fillIsolation = False
process.fixedMuons1stStep.writeIsoDeposits = False
process.fixedMuons1stStep.fillEnergy = False
process.fixedMuons1stStep.storeCrossedHcalRecHits = False
process.fixedMuons1stStep.fillGlobalTrackRefits = False
process.fixedMuons1stStep.fillCaloCompatibility = False
process.fixedMuons1stStep.fillMatching = False
process.fixedMuons1stStep.fillShowerDigis = False
process.fixedMuons1stStep.fillTrackerKink = False
process.fixedMuons1stStep.runArbitrationCleaner = False
process.fixedMuons1stStep.arbitrateTrackerMuons = False
process.fixedMuons1stStep.fillShowerDigis = False

process.fixedMuons = process.muons.clone(InputMuons = cms.InputTag("fixedMuons1stStep"))
process.fixedMuons.PFCandidates = cms.InputTag("particleFlow")
process.fixedMuons.FillPFIsolation = False
process.fixedMuons.FillSelectorMaps = False
process.fixedMuons.FillCosmicsIdMap =  False
process.fixedMuons.FillTimingInfo = False
process.fixedMuons.FillDetectorBasedIsolation = False
process.fixedMuons.FillShoweringInfo = False



## Define two sequences, one for the duplicates and one to evaluate the fixing
process.load("MuonReco-Analysis.MuonReco-analyzer.standaloneDuplicates_cfi")
process.standaloneDuplicates.nameOfOutput = 'outputRAW.root'
process.standaloneDuplicates.muonCollection = 'muons'
process.fixedDuplicates = process.standaloneDuplicates.clone()
process.fixedDuplicates.nameOfOutput = 'outputFIXED.root'
process.fixedDuplicates.muonCollection = 'fixedMuons'

## Running sequence
process.duplicates = cms.Sequence(process.globalMuonLinks*
                                  process.fixedMuons1stStep*
                                  process.fixedMuons*
                                  process.standaloneDuplicates*
                                  process.fixedDuplicates)

#process.duplicates = cms.Sequence(process.standaloneDuplicates*process.fixedDuplicates)

process.p = cms.EndPath(process.duplicates)

