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
#listOfFiles = ['/store/relval/CMSSW_12_6_0_pre2/RelValSingleMuPt100/GEN-SIM-RECO/125X_mcRun3_2022_realistic_v3-v1/2580000/1d014bdd-abc3-41b9-b4f2-2334f009ff7f.root']
listOfFiles = ['/store/relval/CMSSW_12_6_0_pre2/RelValTTbar_14TeV/GEN-SIM-RECO/125X_mcRun3_2022_realistic_v3-v1/2580000/4dd813aa-b6fd-4480-a38c-4cd7709e6a9d.root']
listOfFiles.append('/store/relval/CMSSW_12_6_0_pre2/RelValTTbar_14TeV/GEN-SIM-RECO/125X_mcRun3_2022_realistic_v3-v1/2580000/b0c132e7-3224-4344-897f-a5c558c4f8ce.root')
#listOfFiles = ['']

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( listOfFiles ),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents = cms.untracked.uint32(0)
  )
process.GlobalTag = GlobalTag(process.GlobalTag, '124X_mcRun3_2022_realistic_v5')

## Define the process to run 
## 

process.load("MuonReco-Analysis.MuonReco-analyzer.segmentAnalysis_cfi")
process.segmentAnalysis.nameOfOutput = 'output_segmentAnalysis.root'
process.p = cms.EndPath(process.segmentAnalysis)

