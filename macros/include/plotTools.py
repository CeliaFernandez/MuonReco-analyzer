import ROOT as r
from ROOT import gDirectory, gROOT, gStyle, TLatex
import optparse
import os
import copy

def getObject(filename, key):

    _f = r.TFile(filename)
    _h = _f.Get(key)
    _hcopy = copy.deepcopy(_h)
    _f.Close()

    return _hcopy

def plot2D(histo, output, zlog = False, rebin = False):

    histo.Sumw2()

    # Make the plot bonito
    histo.GetXaxis().SetTitleSize(0.045)
    histo.GetYaxis().SetTitleSize(0.045)

    c1 = r.TCanvas("c1", "", 600, 600)
    c1.cd()

    if zlog:
        c1.SetLogz(1)

    histo.Draw('COLZ')



    # CMS logo
    latex = TLatex()
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextColor(r.kBlack);
    latex.SetTextFont(42);
    latex.SetTextAlign(11);
    latex.SetTextSize(0.055);
    latex.DrawLatex(0.13, 0.93, "#bf{CMS}")

    latexb = TLatex()
    latexb.SetNDC();
    latexb.SetTextAngle(0);
    latexb.SetTextColor(r.kBlack);
    latexb.SetTextFont(42);
    latexb.SetTextAlign(11);
    latexb.SetTextSize(0.04);
    latexb.DrawLatex(0.25, 0.93, "#it{Internal}")
    


    if output[-1] != '/': output = output + '/'
    c1.SaveAs(output + histo.GetName()+'.png')

