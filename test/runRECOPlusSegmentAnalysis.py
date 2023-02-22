# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 -s RAW2DIGI,L1Reco,RECO,RECOSIM --conditions auto:phase1_2022_realistic --datatier AODSIM --eventcontent AODSIM --geometry DB:Extended --era Run3 --pileup Run3_Flat55To75_PoissonOOTPU --pileup_input das:/RelValMinBias_14TeV/CMSSW_12_5_0_pre4-124X_mcRun3_2022_realistic_v10_BS2022-v1/GEN-SIM --filein file:step2.root --fileout file:step3.root -n 50 --customise=Validation/Performance/TimeMemoryInfo.py
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process('RECO',Run3)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mix_Run3_Flat55To75_PoissonOOTPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.RecoSim_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(50),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:step2.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    makeTriggerResults = cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:50'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)


# Other statements
process.mix.input.fileNames = cms.untracked.vstring(['/store/relval/CMSSW_12_5_0_pre4/RelValMinBias_14TeV/GEN-SIM/124X_mcRun3_2022_realistic_v10_BS2022-v1/10000/2e6c6adc-61d1-464a-b320-105dfa3849f7.root', '/store/relval/CMSSW_12_5_0_pre4/RelValMinBias_14TeV/GEN-SIM/124X_mcRun3_2022_realistic_v10_BS2022-v1/10000/499ded67-cabc-4cfc-aa48-5de71bf20d2f.root', '/store/relval/CMSSW_12_5_0_pre4/RelValMinBias_14TeV/GEN-SIM/124X_mcRun3_2022_realistic_v10_BS2022-v1/10000/63fea605-cb9e-4822-9e85-84ecf6bd28d8.root', '/store/relval/CMSSW_12_5_0_pre4/RelValMinBias_14TeV/GEN-SIM/124X_mcRun3_2022_realistic_v10_BS2022-v1/10000/6e724aec-d658-4230-9e9a-861bd9b9f70e.root', '/store/relval/CMSSW_12_5_0_pre4/RelValMinBias_14TeV/GEN-SIM/124X_mcRun3_2022_realistic_v10_BS2022-v1/10000/921e189b-c678-4c67-abd4-d5947770e7ee.root', '/store/relval/CMSSW_12_5_0_pre4/RelValMinBias_14TeV/GEN-SIM/124X_mcRun3_2022_realistic_v10_BS2022-v1/10000/c212d1f0-678d-4e2e-b5d7-bba59b7ec74e.root'])
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2022_realistic', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.recosim_step = cms.Path(process.recosim)

# Load analyzer
process.load("MuonReco-Analysis.MuonReco-analyzer.segmentAnalysis_cfi")
process.segmentAnalysis.nameOfOutput = 'output_segmentAnalysis.root'
process.endjob_step = cms.EndPath(process.segmentAnalysis)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.recosim_step,process.endjob_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)


# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
