import ROOT as r
from ROOT import gDirectory, gROOT, gStyle, gPad
import optparse
import os
from include.plotTools import *

if __name__ == "__main__":

    r.gROOT.LoadMacro('include/tdrstyle.C')
    r.gROOT.SetBatch(1)
    r.setTDRStyle()
    r.gStyle.SetPadRightMargin(0.12)

    filenameFIXED = '/eos/user/f/fernance/DYJetsToMuMu_M-50_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/duplicatedMuons-Jul22/220719_150942/outputFIXED.root'
    filenameRAW = '/eos/user/f/fernance/DYJetsToMuMu_M-50_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/duplicatedMuons-Jul22/220719_150942/outputRAW.root'

    ## Plot histograms
    h_RAW_nSTA = getObject(filenameRAW, 'h_nSTA')
    h_RAW_nGLB = getObject(filenameRAW, 'h_nGLB')
    h_RAW_STA_eta = getObject(filenameRAW, 'h_STA_eta')
    h_RAW_GLB_eta = getObject(filenameRAW, 'h_GLB_eta')
    h_RAW_eta1_eta2 = getObject(filenameRAW, 'h_etaGLB_etaSTA')
    h_RAW_eta1_eta2_low = getObject(filenameRAW, 'h_etaGLB_etaSTA_loweta')
    h_RAW_nDT1_nDT2 = getObject(filenameRAW, 'h_nDTGLB_nDTSTA')
    h_RAW_nCSC1_nCSC2 = getObject(filenameRAW, 'h_nCSCGLB_nCSCSTA')
    h_RAW_nHits1_nHits2 = getObject(filenameRAW, 'h_nHitsGLB_nHitsSTA')

    h_FIXED_nSTA = getObject(filenameFIXED, 'h_nSTA')
    h_FIXED_nGLB = getObject(filenameFIXED, 'h_nGLB')
    h_FIXED_STA_eta = getObject(filenameFIXED, 'h_STA_eta')
    h_FIXED_GLB_eta = getObject(filenameFIXED, 'h_GLB_eta')
    h_FIXED_eta1_eta2 = getObject(filenameFIXED, 'h_etaGLB_etaSTA')
    h_FIXED_eta1_eta2_low = getObject(filenameFIXED, 'h_etaGLB_etaSTA_loweta')
    h_FIXED_nDT1_nDT2 = getObject(filenameFIXED, 'h_nDTGLB_nDTSTA')
    h_FIXED_nCSC1_nCSC2 = getObject(filenameFIXED, 'h_nCSCGLB_nCSCSTA')
    h_FIXED_nHits1_nHits2 = getObject(filenameFIXED, 'h_nHitsGLB_nHitsSTA')

    ## Tune
    h_RAW_STA_eta.Rebin(2)
    h_RAW_GLB_eta.Rebin(2)
    h_FIXED_STA_eta.Rebin(2)
    h_FIXED_GLB_eta.Rebin(2)


    plot2D(h_RAW_eta1_eta2, '/eos/user/f/fernance/www/MuonPOG/MuonDuplicates-220622/', zlog = True)
    plot2D(h_RAW_eta1_eta2_low, '/eos/user/f/fernance/www/MuonPOG/MuonDuplicates-220622/', maxDigits = 3)
    plot2D(h_RAW_nDT1_nDT2, '/eos/user/f/fernance/www/MuonPOG/MuonDuplicates-220622/')
    plot2D(h_RAW_nCSC1_nCSC2, '/eos/user/f/fernance/www/MuonPOG/MuonDuplicates-220622/')
    plot2D(h_RAW_nHits1_nHits2, '/eos/user/f/fernance/www/MuonPOG/MuonDuplicates-220622/')

    plotValidation(h_FIXED_nSTA, h_RAW_nSTA, '/eos/user/f/fernance/www/MuonPOG/MuonDuplicates-220622/', 'Fixed', 'Original', "", ylog = False)
    plotValidation(h_FIXED_STA_eta, h_RAW_STA_eta, '/eos/user/f/fernance/www/MuonPOG/MuonDuplicates-220622/', 'Fixed', 'Original', "", ylog = False)
    plotValidation(h_FIXED_nGLB, h_RAW_nGLB, '/eos/user/f/fernance/www/MuonPOG/MuonDuplicates-220622/', 'Fixed', 'Original', "", ylog = False)
    plotValidation(h_FIXED_GLB_eta, h_RAW_GLB_eta, '/eos/user/f/fernance/www/MuonPOG/MuonDuplicates-220622/', 'Fixed', 'Original', "", ylog = False)
