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

    filenamePU70 = '/eos/user/f/fernance/MuonPU/Mu_FlatPt-1to1000-gun/muonAnalysis-PU70/230125_113712/output.root'
    filenamePU60 = '/eos/user/f/fernance/MuonPU/Mu_FlatPt-1to1000-gun/muonAnalysis-PU60/230125_121142/output.root'

    ## Plot histograms
    h_mu_pt_70 = getObject(filenamePU70, 'h_mu_pt')
    h_mu_eta_70 = getObject(filenamePU70, 'h_mu_eta')
    h_gen_pt_70 = getObject(filenamePU70, 'h_gen_pt')
    h_gen_eta_70 = getObject(filenamePU70, 'h_gen_eta')
    e_mu_pt_70 = getObject(filenamePU70, 'e_mu_pt')
    e_mu_eta_70 = getObject(filenamePU70, 'e_mu_eta')

    h_mu_pt_60 = getObject(filenamePU60, 'h_mu_pt')
    h_mu_eta_60 = getObject(filenamePU60, 'h_mu_eta')
    h_gen_pt_60 = getObject(filenamePU60, 'h_gen_pt')
    h_gen_eta_60 = getObject(filenamePU60, 'h_gen_eta')
    e_mu_pt_60 = getObject(filenamePU60, 'e_mu_pt')
    e_mu_eta_60 = getObject(filenamePU60, 'e_mu_eta')

    weight = h_gen_pt_60.GetEntries()/h_gen_pt_70.GetEntries()   
    h_mu_pt_70.Scale(weight)
    h_mu_eta_70.Scale(weight)
    h_gen_pt_70.Scale(weight)
    h_gen_eta_70.Scale(weight)

    plotValidation('h_mu_pt', h_mu_pt_70, h_mu_pt_60, '/eos/user/f/fernance/www/MuonPOG/MuonPerformancePU/', 'Poisson70KeepRAW', 'Poisson60KeepRAW', "", ylog = False)
    plotValidation('h_mu_eta', h_mu_eta_70, h_mu_eta_60, '/eos/user/f/fernance/www/MuonPOG/MuonPerformancePU/', 'Poisson70KeepRAW', 'Poisson60KeepRAW', "", ylog = False)
    #plotValidation('h_gen_pt', h_gen_pt_70, h_gen_pt_60, '/eos/user/f/fernance/www/MuonPOG/MuonPerformancePU/', 'Poisson70KeepRAW', 'Poisson60KeepRAW', "", ylog = False)
    plotValidation('h_gen_pt_w', h_gen_pt_70, h_gen_pt_60, '/eos/user/f/fernance/www/MuonPOG/MuonPerformancePU/', 'Poisson70KeepRAW', 'Poisson60KeepRAW', "", ylog = False)
    plotValidation('h_gen_eta_w', h_gen_eta_70, h_gen_eta_60, '/eos/user/f/fernance/www/MuonPOG/MuonPerformancePU/', 'Poisson70KeepRAW', 'Poisson60KeepRAW', "", ylog = False)
    plotEfficiency('e_mu_pt', e_mu_pt_70, e_mu_pt_60, '/eos/user/f/fernance/www/MuonPOG/MuonPerformancePU/', 'Poisson70KeepRAW', 'Poisson60KeepRAW', "", ylog = False)
    plotEfficiency('e_mu_eta', e_mu_eta_70, e_mu_eta_60, '/eos/user/f/fernance/www/MuonPOG/MuonPerformancePU/', 'Poisson70KeepRAW', 'Poisson60KeepRAW', "", ylog = False)
