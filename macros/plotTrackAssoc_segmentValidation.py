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

    filenameFIXED = '/afs/cern.ch/work/f/fernance/private/MuonPOG/L3-RECO/MuonReco-analysis/Analysis/CMSSW_12_6_0_pre2/src/MuonReco-Analysis/MuonReco-analyzer/output_segmentAnalysis.root'
    filenameRAW = '/afs/cern.ch/work/f/fernance/private/MuonPOG/L3-RECO/MuonReco-analysis/Analysis/CMSSW_12_6_0_pre2/src/MuonReco-Analysis/MuonReco-analyzer/output_segmentAnalysis_original.root'

    ## Plot histograms
    RAW_trackerMuons_size             = getObject(filenameRAW, 'h_trackerMuons_size')
    RAW_trackerMuons_numberOfMatches  = getObject(filenameRAW, 'h_trackerMuons_numberOfMatches')
    RAW_trackerMuons_numberOfSegments = getObject(filenameRAW, 'h_trackerMuons_numberOfSegments')
    RAW_trackerMuons_segmentX         = getObject(filenameRAW, 'h_trackerMuons_segmentX')
    RAW_trackerMuons_segmentY         = getObject(filenameRAW, 'h_trackerMuons_segmentY')

    FIXED_trackerMuons_size             = getObject(filenameFIXED, 'h_trackerMuons_size')
    FIXED_trackerMuons_numberOfMatches  = getObject(filenameFIXED, 'h_trackerMuons_numberOfMatches')
    FIXED_trackerMuons_numberOfSegments = getObject(filenameFIXED, 'h_trackerMuons_numberOfSegments')
    FIXED_trackerMuons_segmentX         = getObject(filenameFIXED, 'h_trackerMuons_segmentX')
    FIXED_trackerMuons_segmentY         = getObject(filenameFIXED, 'h_trackerMuons_segmentY')

    plotValidation(FIXED_trackerMuons_size, RAW_trackerMuons_size, '/eos/user/f/fernance/www/MuonPOG/TrackAssociator-Validation221016/', 'Modified', 'Original', "", ylog = False)
    plotValidation(FIXED_trackerMuons_numberOfMatches, RAW_trackerMuons_numberOfMatches, '/eos/user/f/fernance/www/MuonPOG/TrackAssociator-Validation221016/', 'Modified', 'Original', "", ylog = False)
    plotValidation(FIXED_trackerMuons_numberOfSegments, RAW_trackerMuons_numberOfSegments, '/eos/user/f/fernance/www/MuonPOG/TrackAssociator-Validation221016/', 'Modified', 'Original', "", ylog = False)
    plotValidation(FIXED_trackerMuons_segmentX, RAW_trackerMuons_segmentX, '/eos/user/f/fernance/www/MuonPOG/TrackAssociator-Validation221016/', 'Modified', 'Original', "", ylog = False)
    plotValidation(FIXED_trackerMuons_segmentY, RAW_trackerMuons_segmentY, '/eos/user/f/fernance/www/MuonPOG/TrackAssociator-Validation221016/', 'Modified', 'Original', "", ylog = False)
