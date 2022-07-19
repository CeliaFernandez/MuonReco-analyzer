# MuonReco Analyzer

## How to install

Run the following set of commands:
```
cmsrel CMSSW_10_6_20
cd CMSSW_10_6_20/src
mkdir MuonReco-Analysis
cd MuonReco-Analysis
git clone git@github.com:CeliaFernandez/MuonReco-analyzer.git
scram b -j 8
```

## Skims

Skims devoted to muon studies can be found in test/skims.
Instructions to produce an skim are detailed in this twiki: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookPickEvents

### Duplicated muons
Event list with duplicated muons reported at eta~0 are available in test/skims/eventlist_duplicates_eta0.txt file. Then the skim is produced throught the following command:
```
edmPickEvents.py "/DYJetsToMuMu_M-50_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM" eventlist_duplicates_eta0.txt --crab
```
