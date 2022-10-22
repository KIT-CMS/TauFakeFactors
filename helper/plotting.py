import ROOT
import os

from . import functions as func

def plot_FFs(SRlike_hist, ARlike_hist, dir_name, process, njets):
    c = ROOT.TCanvas("c", "", 850, 700)
   
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
    ratio.Fit("pol1")
    print("-" * 50)

    c.SetLogx()
    ratio.Draw()

    func.check_output_path(os.getcwd() + "/workdir/" + dir_name)
    c.SaveAs("workdir/{}/ff_{}_{}jets.png".format(dir_name, process, njets.replace(" ","")))
    c.SaveAs("workdir/{}/ff_{}_{}jets.pdf".format(dir_name, process, njets.replace(" ","")))
    c.Close()

def plot_histogram(hists, dir_name, process, region, njets):
    c = ROOT.TCanvas("c", "", 850, 700)
   
    ROOT.gStyle.SetOptStat(0) # set off of the histogram statistics box
    ROOT.gStyle.SetTextFont(42) # chosing font, see https://root.cern/root/html534/TAttText.html

    data = hists["data_substracted"]
    data.SetMarkerStyle(20)
    data.SetMarkerSize(1.2)
    data.SetLineWidth(2)
    data.SetLineColor(ROOT.kBlack)
    data.Draw("E")

    mc = hists[process]
    mc.SetLineWidth(2)
    mc.SetFillStyle(1001)
    mc.SetLineColor(ROOT.kBlack)
    mc.SetFillColor(ROOT.kAzure - 9)
    mc.Draw("HIST SAME")
    c.SetLogx()

    func.check_output_path(os.getcwd() + "/workdir/" + dir_name)
    c.SaveAs("workdir/{}/hist_{}_{}_{}jets.png".format(dir_name, process, region, njets.replace(" ","")))
    c.SaveAs("workdir/{}/hist_{}_{}_{}jets.pdf".format(dir_name, process, region, njets.replace(" ","")))
    c.Close()

