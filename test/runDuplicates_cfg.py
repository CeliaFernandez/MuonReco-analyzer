import FWCore.ParameterSet.Config as cms
import os

# Run standaloneDuplicates plugin
#
# @ Celia Fernandez Madrazo 


process = cms.Process("demo")

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
listOfFiles = ['file:' + os.path.join(inputdir, f) for f in os.listdir(inputdir) if '.root' in f]
if (isData):
  process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( listOfFiles ),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents = cms.untracked.uint32(0)
  )
  process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data')

process.load("MuonReco-Analysis.MuonReco-analyzer.standaloneDuplicates_cfi")

#process.options = cms.untracked.PSet()
#process.options.numberOfThreads = cms.untracked.uint32(8)

process.p = cms.EndPath(process.standaloneDuplicates)

