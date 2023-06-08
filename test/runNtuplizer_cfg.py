import FWCore.ParameterSet.Config as cms
import os


process = cms.Process("demo")
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.Services_cff')

# Debug printout and summary.
process.load("FWCore.MessageService.MessageLogger_cfi")

process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True),
  # Set up multi-threaded run. Must be consistent with config.JobType.numCores in crab_cfg.py.
  #numberOfThreads=cms.untracked.uint32(8)
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag

# Select number of events to be processed
nEvents = 1000
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(nEvents) )

# Read events
#listOfFiles = ['file:/eos/user/f/fernance/Muon-POG/TrackAssociatorIssue/CRAB_PrivateMC/step3_iter_fixed/step3_RAW2DIGI_L1Reco_RECO_RECOSIM_PU_400.root']
#listOfFiles = ['file:/eos/user/f/fernance/PR-testing/TrackAssociator/Jun23Iter/ultimateTest/modified2/CMSSW_13_0_0_pre3/src/step3_RAW2DIGI_L1Reco_RECO_RECOSIM_PU.root']
listOfFiles = ['file:/eos/user/f/fernance/PR-testing/TrackAssociator/Jun23Iter/ultimateTest/original/CMSSW_13_0_0_pre3/src/step3_RAW2DIGI_L1Reco_RECO_RECOSIM_PU.root']
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( listOfFiles ),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents = cms.untracked.uint32(0)
  )
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2022_realistic')

## Define the process to run 
## 
process.load("MuonReco-Analysis.MuonReco-analyzer.muontuplizer_cfi")

process.p = cms.EndPath(process.ntuples)

