import FWCore.ParameterSet.Config as cms
import os

#
# Test file to evaluate duplicated muons in reco::Muon 'muons' collection
# 2 steps:
#
#  - Redo the muon merging (MuonIdProducer) with standalone and global muons
#    to create a new collection 'fixedMuons' where fix can be applied
#    (Tracker muons cannot be redone from AOD)
#
#  - Run the standaloneDuplicates.cc plugin on both raw and fixed collections
#    to compare
#
# @ Celia Fernandez Madrazo (Wed 19/07/2022)
#

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
nEvents = -1
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(nEvents) )

# Read events
listOfFiles = ['file:/eos/user/f/fernance/Muon-POG/TrackAssociatorIssue/CRAB_PrivateMC/step3_iter_fixed/step3_RAW2DIGI_L1Reco_RECO_RECOSIM_PU_400.root']
listOfFiles = ['/store/data/Run2022D/Muon/AOD/PromptReco-v2/000/357/734/00000/011a59b3-42bc-493c-bc2b-73cbddcb0e79.root']

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( listOfFiles ),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents = cms.untracked.uint32(0)
  )
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2022_realistic')
process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_Prompt_v4')

## Define the process to run 
## 

process.load("MuonReco-Analysis.MuonReco-analyzer.segmentAnalysis_cfi")
process.segmentAnalysis.nameOfOutput = 'output_segmentAnalysis.root'
process.p = cms.EndPath(process.segmentAnalysis)

