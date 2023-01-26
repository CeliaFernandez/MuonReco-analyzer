import CRABClient
from CRABClient.UserUtilities import config

config = config()

# General
config.General.workArea = 'crab_projects'
config.General.requestName = 'muonAnalysis-PU70'
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.instance = 'prod'

# JobType
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runRECOAnalysis_cfg.py'
config.JobType.maxMemoryMB = 2500
config.JobType.allowUndistributedCMSSW = True
config.JobType.outputFiles = ['output.root']

# Data
config.Data.inputDataset = '/Mu_FlatPt-1to1000-gun/Run3Summer22EEDR-Poisson70KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v1/AODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 100 # splitting to get ~500 output files
config.Data.outLFNDirBase = '/store/user/fernance/MuonPU' # modify accordingly
config.Data.publication = False
config.Data.outputDatasetTag = 'muonAnalysis-PU70'

# Site
config.Site.storageSite = 'T3_CH_CERNBOX'
