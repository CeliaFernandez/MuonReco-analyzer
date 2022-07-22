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

def plot2D(histo, output, zlog = False, rebin = False, maxDigits = False):

    histo.Sumw2()

    # Make the plot bonito
    histo.GetXaxis().SetTitleSize(0.045)
    histo.GetYaxis().SetTitleSize(0.045)
    if maxDigits:
        histo.GetXaxis().SetMaxDigits(maxDigits)
        histo.GetYaxis().SetMaxDigits(maxDigits)


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
    if maxDigits:
        latex.DrawLatex(0.23, 0.93, "#bf{CMS}")
    else:
        latex.DrawLatex(0.13, 0.93, "#bf{CMS}")

    latexb = TLatex()
    latexb.SetNDC();
    latexb.SetTextAngle(0);
    latexb.SetTextColor(r.kBlack);
    latexb.SetTextFont(42);
    latexb.SetTextAlign(11);
    latexb.SetTextSize(0.04);
    if maxDigits:
        latexb.DrawLatex(0.35, 0.93, "#it{Internal}")
    else:
        latexb.DrawLatex(0.25, 0.93, "#it{Internal}")
    


    if output[-1] != '/': output = output + '/'
    c1.SaveAs(output + histo.GetName()+'.png')


def plotValidation(target, reference, output, tlabel, rlabel, relval, ylog = False, rebin = False):

    target.Sumw2()
    reference.Sumw2()

    if rebin:
        target.Rebin(rebin)
        reference.Rebin(rebin)

    target.SetMarkerSize(1)
    target.SetMarkerStyle(24)
    target.SetMarkerColor(r.kBlue)
    target.SetLineColor(r.kBlue)
    target.SetLineWidth(2)
    target.GetYaxis().SetLabelSize(0.045)

    reference.SetMarkerSize(1)
    reference.SetMarkerStyle(20)
    reference.SetLineWidth(2)
    reference.SetLineColor(r.kBlack)
    reference.SetMarkerColor(r.kBlack)
    reference.GetYaxis().SetTitleSize(0.045)
    reference.GetYaxis().SetLabelSize(0.045)
    reference.GetXaxis().SetLabelSize(0)
    reference.GetXaxis().SetTitle(reference.GetTitle())

    ratio = target.Clone(target.GetName() + '_ratio')
    ratio.Reset()
    for n in range(1, target.GetNbinsX() + 1):
        tv = target.GetBinContent(n)
        rv = reference.GetBinContent(n)
        te = target.GetBinError(n)
        re = reference.GetBinError(n)
        if tv != 0.0 and rv != 0.0:
            value = tv/rv
            error = value * ( (te/tv)**2 + (re/rv)**2 )**0.5
        else:
            value = 0.0
            error = 0.0
        ratio.SetBinContent(n, value)
        ratio.SetBinError(n, error)


    ratio.SetTitle(";"+reference.GetXaxis().GetTitle()+";Ratio")
    ratio.GetYaxis().CenterTitle()
    ratio.GetYaxis().SetTitleOffset(0.4)
    ratio.GetYaxis().SetTitleSize(0.14)
    ratio.GetYaxis().SetLabelSize(0.14)
    ratio.GetXaxis().SetLabelSize(0.14)
    ratio.GetXaxis().SetTitleSize(0.14)
    ratio.SetMarkerColor(r.kRed)
    ratio.SetLineColor(r.kRed)
    ratio.SetLineWidth(2)
    #ratio.SetFillColor(r.kRed)
    ratio.SetMarkerStyle(20)
    ratio.Sumw2()

    target.SetTitle(';;')
    reference.SetTitle(';;Counts')
    ymax = max([target.GetMaximum(), reference.GetMaximum()])
    if not ylog:
        reference.SetMaximum(1.4*ymax)
        reference.SetMinimum(0.0)
    else:
        reference.SetMaximum(10.*ymax)
        reference.SetMinimum(0.1)

    c1 = r.TCanvas("c1", "", 550, 600)
    c1.cd()

    pad1 = r.TPad("pad1", "pad1", 0, 0.25, 1, 1.0)
    pad1.SetBottomMargin(0.03)
    if ylog:
        pad1.SetLogy(1)
    pad1.Draw()

    ### pad 2 drawing
    r.gStyle.SetOptStat(0)
    pad2 = r.TPad("pad2", "pad2", 0, 0.0, 1, 0.25)
    pad2.SetTopMargin(0.0);
    pad2.SetBottomMargin(0.30);
    pad2.Draw();

    ### pad 1 drawing
    pad1.cd()
    reference.Draw('HIST')
    target.Draw('P SAMES')

    legend = r.TLegend(0.5, 0.76, 0.85, 0.87)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.035)
    legend.SetLineWidth(0)
    legend.SetBorderSize(0)
    legend.AddEntry(reference, rlabel + " ({0})".format(int(reference.GetEntries())), 'l')
    legend.AddEntry(target, tlabel + " ({0})".format(int(target.GetEntries())), 'pl')
    legend.Draw()


    ### pad2 drawing
    pad2.cd()
    ratio.Draw("P,SAME")
    line = r.TLine(ratio.GetBinLowEdge(1), 1, ratio.GetBinLowEdge(ratio.GetNbinsX()+1), 1)
    line.Draw("Same")


    ### RelVal text
    pad1.cd()
    rvlabel = r.TLatex()
    rvlabel.SetNDC()
    rvlabel.SetTextAngle(0)
    rvlabel.SetTextColor(r.kBlack)
    rvlabel.SetTextFont(42)
    rvlabel.SetTextAlign(31)
    rvlabel.SetTextSize(0.035)
    rvlabel.DrawLatex(0.85, 0.935, relval)

    ## CMS logo
    latex = TLatex()
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextColor(r.kBlack);
    latex.SetTextFont(42);
    latex.SetTextAlign(11);
    latex.SetTextSize(0.065);
    latex.DrawLatex(0.17, 0.83, "#bf{CMS}")

    latexb = TLatex()
    latexb.SetNDC();
    latexb.SetTextAngle(0);
    latexb.SetTextColor(r.kBlack);
    latexb.SetTextFont(42);
    latexb.SetTextAlign(11);
    latexb.SetTextSize(0.042);
    latexb.DrawLatex(0.17, 0.78, "#it{Internal}")

    ## Save the plot
    if output[-1] != '/': output = output + '/'
    c1.SaveAs(output + target.GetName()+'.png')








