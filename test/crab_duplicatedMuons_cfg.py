import CRABClient
from CRABClient.UserUtilities import config

config = config()

# General
config.General.workArea = 'crab_projects'
config.General.requestName = 'duplicatedMuons-results'
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.instance = 'prod'

# JobType
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runDuplicates_cfg.py'
config.JobType.maxMemoryMB = 2500
config.JobType.allowUndistributedCMSSW = True
config.JobType.outputFiles = ['outputRAW.root', 'outputFIXED.root']

# Data
config.Data.inputDataset = '/DYJetsToMuMu_M-50_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 12 # splitting to get ~500 output files
config.Data.outLFNDirBase = '/store/user/fernance/' # modify accordingly
config.Data.publication = False
config.Data.outputDatasetTag = 'duplicatedMuons-Jul22'

# Site
config.Site.storageSite = 'T3_CH_CERNBOX'
