from typing import ValuesView
import ROOT
import os
import array

from . import functions as func

def plot_FFs(SRlike_hist, ARlike_hist, config, process, njets):
    c = ROOT.TCanvas("c", "", 700, 700)
    c.SetRightMargin(0.05)
    c.SetLeftMargin(0.16)
    c.SetBottomMargin(0.12)

    # pad = ROOT.TPad( "pad", "",  0.03, 0.03, 0.97, 0.97)
    # pad.Draw()
    # pad.cd()

    ROOT.gROOT.SetStyle("CMS")
    ROOT.gStyle.SetOptStat(0) # set off of the histogram statistics box
    ROOT.gStyle.SetTextFont(42) # chosing font, see https://root.cern/root/html534/TAttText.html

    ratio = SRlike_hist.Clone()
    ratio.Divide(ARlike_hist)

    ratio.SetAxisRange(0, 0.5, "Y")
    ratio.SetMarkerStyle(20)
    ratio.SetMarkerSize(1.2)
    ratio.SetLineWidth(2)
    ratio.SetLineColor(ROOT.kBlack)
    ratio.SetTitle("channel: et, {} jets".format(njets.replace(" ","")))
    ratio.GetXaxis().SetMoreLogLabels()
    ratio.GetXaxis().SetNoExponent()
    ratio.GetYaxis().SetTitle("FF({})".format(process))
    ratio.GetYaxis().SetLabelSize(0.045)
    ratio.GetXaxis().SetLabelSize(0.045)
    ratio.GetYaxis().SetTitleSize(0.05)
    ratio.GetXaxis().SetTitleSize(0.05)
    # ratio.SetAxisRange(0.,700.,axis="X")
    
    fit = ratio.Fit("pol1", "SF")  # https://root.cern/doc/master/classTH1.html#a7e7d34c91d5ebab4fc9bba3ca47dabdd
    cl68 = fit.GetConfidenceIntervals(cl=0.68)
    chi2 = fit.Chi2()
    dof = fit.Ndf()
    print("Confidence intervals:", cl68)


    # bindata = ROOT.Fit.BinData()
    # print(bindata)
    # ROOT.Fit.FillData(bindata, ratio)
    # print(bindata.GetValue())
    # cl68 = fit.GetConfidenceIntervals(bindata, cl=0.68)
    # print(cl68)
    
    print("-" * 50)

    if config["target_process"][process]["tau_pt_logx"]:
        c.SetLogx()
    
    ratio.Draw()

    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextFont(42)
    text.SetTextSize(0.035)
    text.DrawLatex(0.6, 0.8, "#chi^{} / N(dof) = {} / {}".format(2, round(chi2, 2),dof))

    # fit = ROOT.TH1D("fit", "Fit with .95 conf.band", 4, array.array("d", config["tau_pt_bins"]))
    # ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(fit, 0.68)
    # # fit.SetStats(False)
    # # fit.SetFillStyle(1001)
    # fit.SetFillColorAlpha(2, 0.2)
    # fit.SetAxisRange(0, 0.5, "Y")
    
    # #cl68.Draw("E3 SAME")
    # fit.SetTitle("channel: et, {} jets".format(njets.replace(" ","")))
    # fit.Draw("E3 SAME")
    
    func.check_output_path(os.getcwd() + "/workdir/" + config["workdir_name"])
    c.SaveAs("workdir/{}/ff_{}_{}jets.png".format(config["workdir_name"], process, njets.replace(" ","")))
    c.SaveAs("workdir/{}/ff_{}_{}jets.pdf".format(config["workdir_name"], process, njets.replace(" ","")))
    c.Close()

def plot_histogram(hists, config, region, process, njets):
    c = ROOT.TCanvas("c", "", 850, 700)
   
    ROOT.gStyle.SetOptStat(0) # set off of the histogram statistics box
    ROOT.gStyle.SetTextFont(42) # chosing font, see https://root.cern/root/html534/TAttText.html

    mc = hists[process]
    mc.SetLineWidth(2)
    mc.SetFillStyle(1001)
    mc.SetLineColor(ROOT.kBlack)
    mc.SetFillColor(ROOT.kAzure - 9)
    mc.SetTitle("channel: et, {} jets".format(njets.replace(" ","")))
    mc.SetYTitle("Events")
    mc.SetMaximum(mc.GetMaximum()+ 0.5* mc.GetMaximum())
    mc.Draw("HIST")
    if config["target_process"][process]["tau_pt_logx"]:
        c.SetLogx()

    data = hists["data_substracted"]
    data.SetMarkerStyle(20)
    data.SetMarkerSize(1.2)
    data.SetLineWidth(2)
    data.SetLineColor(ROOT.kBlack)
    
    data.Draw("E SAME")

    legend = ROOT.TLegend(0.55, 0.60, 0.88, 0.88)
    legend.SetTextFont(42)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.04)
    legend.SetTextAlign(32)
    legend.AddEntry(data, "Data")
    legend.AddEntry(mc, process)
    legend.Draw("SAME")

    func.check_output_path(os.getcwd() + "/workdir/" + config["workdir_name"])
    c.SaveAs("workdir/{}/hist_{}_{}_{}jets.png".format(config["workdir_name"], process, region, njets.replace(" ","")))
    c.SaveAs("workdir/{}/hist_{}_{}_{}jets.pdf".format(config["workdir_name"], process, region, njets.replace(" ","")))
    c.Close()

def plot_full_histogram(hists, config, region, process, njets):
    c = ROOT.TCanvas("c", "", 850, 700)
   
    ROOT.gStyle.SetOptStat(0) # set off of the histogram statistics box
    ROOT.gStyle.SetTextFont(42) # chosing font, see https://root.cern/root/html534/TAttText.html

    stack = ROOT.THStack()
    ttbar_J = hists["ttbar_J"]
    ttbar_L = hists["ttbar_L"]
    dib_J = hists["diboson_J"]
    dib_L = hists["diboson_L"]
    DYjets_J = hists["DYjets_J"]
    DYjets_L = hists["DYjets_L"]
    Wjets = hists["Wjets"]
    emb = hists["embedding"]
    qcd = hists["QCD"]
    for h, color in zip([qcd, dib_J, dib_L, Wjets, ttbar_J, ttbar_L, DYjets_J, DYjets_L, emb], [(204,204,204), (100, 192, 232), (100, 212, 232), (137,160,44), (155, 152, 204), (175, 152, 204), (191, 34, 41), (211, 34, 41), (255,204,0)]): 
        h.SetLineWidth(1)
        h.SetFillStyle(1001)
        h.SetLineColor(ROOT.kBlack)
        h.SetFillColor(ROOT.TColor.GetColor(*color))
        stack.Add(h)
    stack.SetTitle("channel: et, {} jets".format(njets.replace(" ","")))
    stack.SetMaximum(stack.GetMaximum()+ 0.5* stack.GetMaximum())
    stack.Draw("HIST")
    if config["target_process"][process]["tau_pt_logx"]:
        c.SetLogx()
    #c.SetLogy()

    data = hists["data"]
    data.SetMarkerStyle(20)
    data.SetMarkerSize(1.2)
    data.SetLineWidth(2)
    data.SetLineColor(ROOT.kBlack)
    data.Draw("E SAME")

    legend = ROOT.TLegend(0.55, 0.60, 0.88, 0.88)
    legend.SetTextFont(42)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.04)
    legend.SetTextAlign(32)
    legend.AddEntry(data, "Data")
    legend.AddEntry(emb, "Embedded")
    legend.AddEntry(ttbar_J, "ttbar J")
    legend.AddEntry(ttbar_L, "ttbar L")
    legend.AddEntry(Wjets, "Wjets")
    legend.AddEntry(DYjets_J, "DYjets J")
    legend.AddEntry(DYjets_L, "DYjets L")
    legend.AddEntry(dib_J, "diboson J")
    legend.AddEntry(dib_L, "diboson L")
    legend.AddEntry(qcd, "QCD")
    legend.Draw("SAME")

    func.check_output_path(os.getcwd() + "/workdir/" + config["workdir_name"])
    c.SaveAs("workdir/{}/full_hist_{}_{}_{}jets.png".format(config["workdir_name"], process, region, njets.replace(" ","")))
    c.SaveAs("workdir/{}/full_hist_{}_{}_{}jets.pdf".format(config["workdir_name"], process, region, njets.replace(" ","")))
    c.Close()

def plot_ratio_histogram(hists, config, region, process, njets):
    # create required parts
    mc = hists[process]
    mc.SetLineColor(ROOT.kBlack)
    mc.SetFillColor(ROOT.kAzure - 9)
    mc.SetLineWidth(2)
    mc.GetYaxis().SetTitleSize(20)
    mc.GetYaxis().SetTitleFont(43)
    mc.GetYaxis().SetTitleOffset(1.55)
    mc.SetStats(0)

    data = hists["data_substracted"]
    data.SetLineColor(ROOT.kBlack)
    data.SetLineWidth(2)

    ratio = data.Clone("ratio")
    ratio.SetLineColor(ROOT.kBlack)
    ratio.SetMarkerStyle(21)
    ratio.SetTitle("")
    ratio.SetMinimum(0.75)
    ratio.SetMaximum(1.25)
    # Set up plot for markers and errors
    ratio.Sumw2()
    ratio.SetStats(0)
    ratio.Divide(mc)
 
    # Adjust y-axis settings
    y = ratio.GetYaxis()
    y.SetTitle("data/mc")
    y.SetNdivisions(505)
    y.SetTitleSize(20)
    y.SetTitleFont(43)
    y.SetTitleOffset(1.55)
    y.SetLabelFont(43)
    y.SetLabelSize(15)
 
    # Adjust x-axis settings
    x = ratio.GetXaxis()
    x.SetTitleSize(20)
    x.SetTitleFont(43)
    x.SetTitleOffset(4.0)
    x.SetLabelFont(43)
    x.SetLabelSize(15)

    c = ROOT.TCanvas("c", "canvas", 800, 800)
    # Upper histogram plot is pad1
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(1)  # joins upper and lower plot
    pad1.Draw()
    # Lower ratio plot is pad2
    c.cd()  # returns to main canvas before defining pad2
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
    pad2.SetTopMargin(0)  # joins upper and lower plot
    pad2.SetBottomMargin(0.2)
    pad2.SetGridy()
    pad2.Draw()
 
    # draw everything
    pad1.cd()
    mc.Draw("HIST")
    data.Draw("SAME")
    # to avoid clipping the bottom zero, redraw a small axis
    #mc.GetYaxis().SetLabelSize(0.0)
    # axis = ROOT.TGaxis(-5, 20, -5, 220, 20, 220, 510, "")
    # axis.SetLabelFont(43)
    # axis.SetLabelSize(15)
    # axis.Draw()
    pad2.cd()
    ratio.Draw("ep")

    func.check_output_path(os.getcwd() + "/workdir/" + config["workdir_name"])
    c.SaveAs("workdir/{}/ratio_hist_{}_{}_{}jets.png".format(config["workdir_name"], process, region, njets.replace(" ","")))
    c.SaveAs("workdir/{}/ratio_hist_{}_{}_{}jets.pdf".format(config["workdir_name"], process, region, njets.replace(" ","")))
    c.Close()