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

    filename = '/afs/cern.ch/work/f/fernance/private/MuonPOG/L3-RECO/MuonReco-analysis/TrackAssociator-study/CMSSW_12_6_0_pre2/src/Results.root'

    ## Plot histograms
    h_nExtrapolatedTracks_init             = getObject(filename, 'h_nExtrapolatedTracks_init')
    h_nExtrapolatedTracks_mod             = getObject(filename, 'h_nExtrapolatedTracks_mod')
    h_nExtrapolations_init             = getObject(filename, 'h_nExtrapolations_init')
    h_nExtrapolations_mod             = getObject(filename, 'h_nExtrapolations_mod')
    h_EventSummary_init             = getObject(filename, 'h_EventSummary_init')
    h_EventSummary_mod             = getObject(filename, 'h_EventSummary_mod')


    plotComparison(h_nExtrapolatedTracks_mod, h_nExtrapolatedTracks_init, '/eos/user/f/fernance/www/MuonPOG/TrackAssociator-Validation221016/', 'Modified', 'Original', "", ylog = False)
    plotComparison(h_nExtrapolations_mod, h_nExtrapolations_init, '/eos/user/f/fernance/www/MuonPOG/TrackAssociator-Validation221016/', 'Modified', 'Original', "", ylog = False)
    plotComparison(h_EventSummary_mod, h_EventSummary_init, '/eos/user/f/fernance/www/MuonPOG/TrackAssociator-Validation221016/', 'Modified', 'Original', "", ylog = False, ww = 700)



