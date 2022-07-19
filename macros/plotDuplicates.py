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

    filename = '/afs/cern.ch/work/f/fernance/private/MuonPOG/L3-RECO/MuonReco-analysis/Analysis/CMSSW_10_6_20/src/MuonReco-Analysis/MuonReco-analyzer/output.root'

    h_pt1_pt2 = getObject(filename, 'h_ptGLB_ptSTA')
    h_eta1_eta2 = getObject(filename, 'h_etaGLB_etaSTA')
    h_eta1_eta2_low = getObject(filename, 'h_etaGLB_etaSTA_loweta')
    h_phi1_phi2 = getObject(filename, 'h_phiGLB_phiSTA')
    h_nDT1_nDT2 = getObject(filename, 'h_nDTGLB_nDTSTA')
    h_nCSC1_nCSC2 = getObject(filename, 'h_nCSCGLB_nCSCSTA')
    h_nHits1_nHits2 = getObject(filename, 'h_nHitsGLB_nHitsSTA')

    plot2D(h_pt1_pt2, '/eos/user/f/fernance/www/MuonPOG/MuonDuplicates-220622/')
    plot2D(h_eta1_eta2, '/eos/user/f/fernance/www/MuonPOG/MuonDuplicates-220622/')
    plot2D(h_eta1_eta2_low, '/eos/user/f/fernance/www/MuonPOG/MuonDuplicates-220622/')
    plot2D(h_phi1_phi2, '/eos/user/f/fernance/www/MuonPOG/MuonDuplicates-220622/')
    plot2D(h_nDT1_nDT2, '/eos/user/f/fernance/www/MuonPOG/MuonDuplicates-220622/')
    plot2D(h_nCSC1_nCSC2, '/eos/user/f/fernance/www/MuonPOG/MuonDuplicates-220622/')
    plot2D(h_nHits1_nHits2, '/eos/user/f/fernance/www/MuonPOG/MuonDuplicates-220622/')

