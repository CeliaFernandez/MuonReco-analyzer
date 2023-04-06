![alt text](https://github.com/CeliaFernandez/MuonReco-analyzer/blob/main/.header.png?raw=true)
# MuonReco Analyzer

## How to install

Run the following set of commands:
```
cmsrel CMSSW_X_Y_Z
cd CMSSW_X_Y_Z/src
mkdir MuonReco-Analysis
cd MuonReco-Analysis
git clone git@github.com:CeliaFernandez/MuonReco-analyzer.git
scram b -j 8
```

## Plugins

### Duplicated muons
(Working release is ```CMSSW_10_6_20```)
Code developed to study the presence of standalone muons that are made from the same set of hits (duplicates). This shouldn't be allowed by CMSSW as cleaning is applied on standalone tracks.

Plugin is available in ```plugins/standaloneDuplicates.cc``` and config parameters in ```python/standaloneDuplicates_cfi.py```.

It can be run iteratively through:
```
cmsRun test/runDuplicates_cfg.py
```
or with crab:
```
crab submit crab_duplicatedMuons_cfg.py
```


Note that input files and some absolute paths may need to be modified.

## Skims

Skims devoted to muon studies can be found in test/skims.
Instructions to produce an skim are detailed in this twiki: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookPickEvents

### Duplicated muons
Event list with duplicated muons reported at eta~0 are available in test/skims/eventlist_duplicates_eta0.txt file. Then the skim is produced throught the following command:
```
edmPickEvents.py "/DYJetsToMuMu_M-50_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM" eventlist_duplicates_eta0.txt --crab
```

## Scripts

Set of useful scripts to use when studying the reco::Muon object.

#### printDisplacedSize_AOD.py and printDisplacedSize_MINI.py

Used to measure the size taken by the displacedMuons and slimmedDisplacedMuons collections.

