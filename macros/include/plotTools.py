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

def plotComparison(name, tree_new, tree_ref, var, cut, nbin, xmin, xmax, label = ''):

    histo_new = r.TH1F(name + "_new", "", nbin, xmin, xmax)
    histo_ref = r.TH1F(name + "_ref", "", nbin, xmin, xmax)
    histo_new.Sumw2()
    histo_ref.Sumw2()
    tree_ref.Project(histo_ref.GetName(), var, cut, "")
    tree_new.Project(histo_new.GetName(), var, cut, "")
    histo_sub = histo_new.Clone(name + "_sub")
    histo_sub.Add(histo_ref, -1)

    ## Tune histos
    histo_ref.SetMarkerSize(1)
    histo_ref.SetMarkerStyle(20)
    histo_ref.SetLineWidth(2)
    histo_ref.SetLineColor(r.kBlack)
    histo_ref.SetMarkerColor(r.kBlack)
    histo_ref.GetYaxis().SetTitleSize(0.045)
    histo_ref.GetYaxis().SetLabelSize(0.045)
    histo_ref.GetXaxis().SetLabelSize(0.045)
    histo_ref.GetXaxis().SetTitle(var)
    histo_ref.GetYaxis().SetTitle('Counts')
    histo_ref.SetMinimum(histo_sub.GetMinimum())

    histo_new.SetMarkerSize(1)
    histo_new.SetMarkerStyle(24)
    histo_new.SetMarkerColor(r.kRed)
    histo_new.SetLineColor(r.kRed)
    histo_new.SetLineWidth(2)
    histo_new.GetYaxis().SetLabelSize(0.045)

    histo_sub.SetMarkerSize(1)
    histo_sub.SetMarkerStyle(20)
    histo_sub.SetLineWidth(2)
    histo_sub.SetLineColor(r.kBlue)
    histo_sub.SetMarkerColor(r.kBlue)
    histo_sub.GetYaxis().SetTitleSize(0.045)

    ## Legend
    legend = r.TLegend(0.55, 0.76, 0.8, 0.87)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.035)
    legend.SetLineWidth(0)
    legend.SetBorderSize(0)
    legend.AddEntry(histo_ref, 'Reference', 'l')
    legend.AddEntry(histo_new, 'New', 'pl')
    legend.AddEntry(histo_sub, 'New - Reference', 'l')

    ## Zero line
    line = r.TLine(xmin, 0, xmax, 0)
    line.SetLineColor(r.kGray+2)
    line.SetLineWidth(2)
    line.SetLineStyle(7)

    ## Plot histos
    c1 = r.TCanvas("c1", "", 700, 600)
    c1.cd()
    #c1.SetLogy(1)

    histo_ref.Draw("HIST")
    line.Draw("SAME")
    histo_new.Draw("P, SAMES")
    histo_sub.Draw("HIST,SAME")
    legend.Draw()

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

    uplabel = TLatex()
    uplabel.SetNDC();
    uplabel.SetTextAngle(0);
    uplabel.SetTextColor(r.kBlack);
    uplabel.SetTextFont(42);
    uplabel.SetTextAlign(11);
    uplabel.SetTextSize(0.032);
    uplabel.DrawLatex(0.13, 0.93, label)


    c1.SaveAs('tupleplots/' + name + '.png')


def plotComparisonRatio(name, tree_new, tree_ref, var, cut, nbin, xmin, xmax, label = '', output = 'ratiotuples', ylog = False):

    histo_new = r.TH1F(name + "_new", "", nbin, xmin, xmax)
    histo_ref = r.TH1F(name + "_ref", "", nbin, xmin, xmax)
    histo_new.Sumw2()
    histo_ref.Sumw2()
    tree_ref.Project(histo_ref.GetName(), var, cut, "")
    tree_new.Project(histo_new.GetName(), var, cut, "")
    histo_sub = histo_new.Clone(name + "_sub")
    histo_sub.Add(histo_ref, -1)

    ## Tune histos
    histo_ref.SetMarkerSize(1)
    histo_ref.SetMarkerStyle(20)
    histo_ref.SetLineWidth(2)
    histo_ref.SetLineColor(r.kBlack)
    histo_ref.SetMarkerColor(r.kBlack)
    histo_ref.GetYaxis().SetTitleSize(0.045)
    histo_ref.GetYaxis().SetLabelSize(0.045)
    histo_ref.GetXaxis().SetLabelSize(0)
    histo_ref.GetXaxis().SetTitle(var)
    histo_ref.GetYaxis().SetTitle('Counts')

    histo_new.SetMarkerSize(1)
    histo_new.SetMarkerStyle(24)
    histo_new.SetMarkerColor(r.kRed)
    histo_new.SetLineColor(r.kRed)
    histo_new.SetLineWidth(2)
    histo_new.GetYaxis().SetLabelSize(0.045)

    histo_sub.SetMarkerSize(1)
    histo_sub.SetMarkerStyle(20)
    histo_sub.SetLineWidth(2)
    histo_sub.SetLineColor(r.kBlue)
    histo_sub.SetMarkerColor(r.kBlue)
    histo_sub.GetYaxis().SetTitleSize(0.14)
    histo_sub.GetXaxis().SetTitleSize(0.14)
    histo_sub.GetYaxis().SetLabelSize(0.14)
    histo_sub.GetXaxis().SetLabelSize(0.14)
    histo_sub.GetXaxis().SetTitle(var)
    histo_sub.GetYaxis().SetTitle('New - Ref')
    histo_sub.GetYaxis().CenterTitle()
    histo_sub.GetYaxis().SetTitleOffset(0.45)
    histo_sub.GetXaxis().SetTitleOffset(1.0)

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
    histo_ref.Draw("HIST")
    histo_new.Draw("P, SAMES")
     

    legend = r.TLegend(0.5, 0.76, 0.8, 0.87)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.035)
    legend.SetLineWidth(0)
    legend.SetBorderSize(0)
    legend.AddEntry(histo_ref, 'Reference', 'l')
    legend.AddEntry(histo_new, 'New', 'pl')
    legend.AddEntry(histo_sub, 'New - Reference', 'l')
    legend.Draw()


    ### pad2 drawing
    pad2.cd()
    histo_sub.Draw("HIST,SAME")
    line = r.TLine(xmin, 0, xmax, 0)
    line.SetLineColor(r.kGray+2)
    line.SetLineWidth(2)
    line.SetLineStyle(7)
    line.Draw("Same")


    ## CMS logo
    pad1.cd()
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

    uplabel = TLatex()
    uplabel.SetNDC();
    uplabel.SetTextAngle(0);
    uplabel.SetTextColor(r.kBlack);
    uplabel.SetTextFont(42);
    uplabel.SetTextAlign(11);
    uplabel.SetTextSize(0.035);
    uplabel.DrawLatex(0.13, 0.93, label)

    ## Save the plot
    if output[-1] != '/': output = output + '/'
    c1.SaveAs(output + name +'.png')



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


def plotValidation(name, target, reference, output, tlabel, rlabel, relval, ylog = False, rebin = False):

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
    #reference.GetXaxis().SetTitle(reference.GetTitle())

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
    ratio.SetFillColor(r.kRed)
    ratio.SetMarkerStyle(20)
    ratio.Sumw2()

    target.SetTitle(';;')
    reference.SetTitle(';;Counts')
    ymax = max([target.GetMaximum(), reference.GetMaximum()])
    if not ylog:
        reference.SetMaximum(1.4*ymax)
        reference.SetMinimum(0.0)
    else:
        reference.SetMaximum(100.*ymax)
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
    reference.Draw('P')
    target.Draw('P SAMES')

    legend = r.TLegend(0.35, 0.76, 0.7, 0.87)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.035)
    legend.SetLineWidth(0)
    legend.SetBorderSize(0)
    legend.AddEntry(reference, rlabel + " ({0})".format(int(reference.GetEntries())), 'pl')
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
    if relval:
        rvlabel.DrawLatex(0.88, 0.935, relval)

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
    c1.SaveAs(output + name +'.png')


def plotEfficiency(name, target, reference, output, tlabel, rlabel, relval, ylog = False, rebin = False):

    target.SetMarkerSize(1)
    target.SetMarkerStyle(24)
    target.SetMarkerColor(r.kBlue)
    target.SetLineColor(r.kBlue)
    target.SetLineWidth(2)
    #target.GetYaxis().SetLabelSize(0.045)

    reference.SetMarkerSize(1)
    reference.SetMarkerStyle(20)
    reference.SetLineWidth(2)
    reference.SetLineColor(r.kBlack)
    reference.SetMarkerColor(r.kBlack)
    #reference.GetYaxis().SetTitleSize(0.045)
    #reference.GetYaxis().SetLabelSize(0.045)
    #reference.GetXaxis().SetLabelSize(0)
    #reference.GetXaxis().SetTitle(reference.GetTitle())

    tmp_num = target.GetTotalHistogram().Clone()
    tmp_den = reference.GetTotalHistogram().Clone()
    for n in range(0,tmp_num.GetNbinsX()):
        tmp_den.SetBinContent(n+1, reference.GetEfficiency(n+1))
        tmp_num.SetBinContent(n+1, target.GetEfficiency(n+1))
        tmp_den.SetBinError(n+1, reference.GetEfficiencyErrorLow(n+1))
        tmp_num.SetBinError(n+1, target.GetEfficiencyErrorLow(n+1))
    ratio = tmp_num.Clone(name+'_ratio')
    ratio.Divide(tmp_den) 

    #ratio.SetTitle(";"+reference.GetXaxis().GetTitle()+";Ratio")
    ratio.GetYaxis().SetTitle("Ratio")
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
    _h = r.TH1F("h", ";;Efficiency", 1, tmp_num.GetBinLowEdge(1), tmp_num.GetBinLowEdge(tmp_num.GetNbinsX()) + tmp_num.GetBinWidth(tmp_num.GetNbinsX()))
    _h.SetMinimum(0.0)
    _h.SetMaximum(1.25)
    _h.GetYaxis().SetTitleSize(0.045)
    _h.GetYaxis().SetTitle(tmp_num.GetYaxis().GetTitle())
    _h.GetYaxis().SetLabelSize(0.045)
    _h.GetXaxis().SetLabelSize(0)
    pad1.cd()
    _h.Draw("AXIS")
    reference.Draw('P, SAME')
    target.Draw('P, SAME')

    legend = r.TLegend(0.35, 0.76, 0.7, 0.87)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.035)
    legend.SetLineWidth(0)
    legend.SetBorderSize(0)
    legend.AddEntry(reference, rlabel, 'pl')
    legend.AddEntry(target, tlabel, 'pl')
    legend.Draw()


    ### pad2 drawing
    pad2.cd()
    #ratio.SetMinimum(0.9)
    #ratio.SetMaximum(1.1)
    ratio.Draw("P")
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
    c1.SaveAs(output + name +'.png')






