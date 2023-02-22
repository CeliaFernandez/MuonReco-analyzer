import FWCore.ParameterSet.Config as cms
import os
from FWCore.ParameterSet.VarParsing import VarParsing

## VarParsing object
#options = VarParsing('python')
#options.register('input', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Set input dir')
#options.register('output', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Set input dir')
#options.parseArguments()


## Process
process = cms.Process("Analysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi")

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.GlobalTag.globaltag = cms.string("102X_upgrade2018_realistic_v15")

process.load("MuonReco-Analysis.MuonReco-analyzer.standardRECOAnalysis_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


inputdirs = ['/eos/cms/store/mc/Run3Summer22EEDR/Mu_FlatPt-1to1000-gun/AODSIM/Poisson60KeepRAW_124X_mcRun3_2022_realistic_postEE_v1-v2/2820000/']
#inputdirs= options.input.split(',')
listOfFiles = []
for inputdir in inputdirs:
    listOfFiles += ['file:' + os.path.join(inputdir, f) for f in os.listdir(inputdir) if '.root' in f]
listOfFiles = listOfFiles[:3]
#listOfFiles = ['file:/eos/user/f/fernance/Muon-POG/displacedCollection-implementation/11834.21/step3_220411/output_0.root']

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
        listOfFiles
    )
)

#process.options = cms.untracked.PSet(
#    numberOfThreads = cms.untracked.uint32(8)
#)


#process.recoValidator.nameOfOutput = options.output + '.root'
#process.recoValidator.DisplacedCollection = 'muons'



## path definitions
process.p      = cms.Path(
    process.recoAnalysis

)


