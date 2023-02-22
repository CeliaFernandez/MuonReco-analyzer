import CRABClient
from CRABClient.UserUtilities import config

config = config()

# General
config.General.workArea = 'crab_projects'
config.General.requestName = 'segmentAnalysis-results'
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.instance = 'prod'

# JobType
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runRECOPlusSegmentAnalysis.py'
config.JobType.maxMemoryMB = 4000
config.JobType.allowUndistributedCMSSW = True
config.JobType.outputFiles = ['output_segmentAnalysis.root']

# Data
config.Data.inputDataset = '/TTbar_14TeV_TuneCP5_fixTrackAssociator/phys_muon-fernance-step2_221020-b77b279da2031f1030dd70c6bd2efe6c/USER'
config.Data.inputDBS = 'phys03'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1 # splitting to get ~995 output files
config.Data.outLFNDirBase = '/store/user/fernance/' # modify accordingly
config.Data.publication = False
config.Data.outputDatasetTag = 'segmentAnalysis-Fall22-uncorrected'

# Site
config.Site.storageSite = 'T3_CH_CERNBOX'
