import ROOT as r
import numpy as np

referenceTime = '/afs/cern.ch/work/f/fernance/private/MuonPOG/L3-RECO/MuonReco-analysis/TrackAssociator-study/CMSSW_12_6_0_pre2/src/timing/step3_TimeMemoryInfo_sta.log'
testTime = '/afs/cern.ch/work/f/fernance/private/MuonPOG/L3-RECO/MuonReco-analysis/TrackAssociator-study/CMSSW_12_6_0_pre2/src/timing/step3_TimeMemoryInfo_mod.log'

with open(referenceTime) as rf:
    refInfo = rf.readlines()

with open(testTime) as tf:
    testInfo = tf.readlines()

report = {}

### Loop over reference
skip1 = False
skip2 = False
for line in refInfo:
    if skip1 and skip2:
        info = line.split(' ')
        while '' in info:
            info.remove('')
        if len(info) == 5:
            collection = info[-1]
            collection.replace('\n', '')
            report[collection] = []
            timeref = info[-2]
            report[collection].append(float(timeref))
        else:
            skip1 = False
            skip2 = False
    if 'TimeReport ---------- Module Summary ---[Real sec]----' in line:
        skip1 = True
    if 'TimeReport  per event     per exec    per visit  Name' in line:
        skip2 = True

### Loop over testing
skip1 = False
skip2 = False
for line in testInfo:
    if skip1 and skip2:
        info = line.split(' ')
        while '' in info:
            info.remove('')
        if len(info) == 5:
            collection = info[-1]
            collection.replace('\n', '')
            timetest = info[-2]
            report[collection].append(float(timetest))
        else:
            skip1 = False
            skip2 = False
    if 'TimeReport ---------- Module Summary ---[Real sec]----' in line:
        skip1 = True
    if 'TimeReport  per event     per exec    per visit  Name' in line:
        skip2 = True

print(report)

r.gROOT.LoadMacro('include/tdrstyle.C')
r.gROOT.SetBatch(1)
r.setTDRStyle()


#### 2D time comparison plot
time_2d = r.TH2F('time_2d', '', 100, np.logspace(-5, 1, 101), 100, np.logspace(-5, 1, 101))
time_1d = r.TH1F('time_1d', 'Time report ; (t_{mod} - t_{ref})/t_{ref} ; Number of sequences', 200, -50, 50)
time_1d_ref = r.TH1F('time_1d_ref', '', 200, -40, 40)
observation = 'muons1stStep'
for key in report.keys():
    coll = report[key]
    time_2d.Fill(coll[0], coll[1])
    if coll[0] > 0.0001:
        merit = (coll[1] - coll[0])/coll[0]*100
        if observation in key:
            print('hola')
            time_1d_ref.Fill(merit)
        else:
            time_1d.Fill(merit)


c1 = r.TCanvas('c1', '', 500, 500)

time_1d.SetMarkerSize(3)
time_1d.SetLineColor(r.kBlue)
time_1d.Draw('HIST')
time_1d_ref.SetMarkerSize(3)
time_1d_ref.SetLineColor(r.kRed)
time_1d_ref.Draw('HIST,SAME')
c1.Print('Prueba.png')

"""
time_2d.SetMarkerSize(1)
time_2d.Draw('P')
c1.SetLogy(1)
c1.SetLogx(1)
c1.Print('Prueba.png')
"""






