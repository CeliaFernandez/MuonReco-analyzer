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

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag

# Select number of events to be processed
nEvents = -1
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(nEvents) )

# Read events
isData = True # Always true by default (running on MC is useless)
listOfFiles = ['/store/mc/RunIISummer20UL16RECO/DYJetsToMuMu_M-50_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/AODSIM/106X_mcRun2_asymptotic_v13-v2/00000/0001EFA7-97DF-4A45-84A9-79448BFD814B.root']

if (isData):
  process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( listOfFiles ),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents = cms.untracked.uint32(0)
  )
  process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mcRun2_asymptotic_v13')


## Load the analysis sequence
process.load("MuonReco-Analysis.MuonReco-analyzer.standaloneDuplicates_cfi")
process.standaloneDuplicates.nameOfOutput = 'outputRAW.root'
process.standaloneDuplicates.muonCollection = 'muons'

## Running sequence
process.p = cms.EndPath(process.standaloneDuplicates)

