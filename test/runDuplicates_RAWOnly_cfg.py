import FWCore.ParameterSet.Config as cms
import os

#
# Test file to run the duplicatedMuons analysis in RAW mode i.e. without 
# rerunning the MuonIdProducer sequences.
#
# @ Celia Fernandez Madrazo (Wed 19/09/2022)
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
isData = True # Always true by default (running on MC is useless)
#listOfFiles = ['/store/data/Run2022G/Muon/AOD/PromptReco-v1/000/362/433/00000/96cf0f3c-8aa0-45db-b898-a3820a80ecba.root']
listOfFiles = ['/store/data/Run2022G/Muon/AOD/PromptReco-v1/000/362/438/00000/79a4613f-bbde-4e41-aefc-d590b2a53f82.root']

if (isData):
  process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( listOfFiles ),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents = cms.untracked.uint32(0)
  )
  process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_Prompt_v10')


## Load the analysis sequence
process.load("MuonReco-Analysis.MuonReco-analyzer.standaloneDuplicates_cfi")
process.standaloneDuplicates.nameOfOutput = 'outputRAW.root'
process.standaloneDuplicates.muonCollection = 'muons'

## Running sequence
process.p = cms.EndPath(process.standaloneDuplicates)

