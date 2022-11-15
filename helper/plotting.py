import os
import ROOT

from . import functions as func

FF_YAxis = {"ttbar": "FF_{t#bar{t}}", "Wjets": "FF_{Wjets}", "QCD": "FF_{QCD}"}
channel_dict = {"et": "e#tau_{h}", "mt": "#mu#tau_{h}", "tt": "#tau_{h}#tau_{h}"}
era_dict = {
    "2017": "41.48 fb^{-1} (2017, 13 TeV)",
    "2018": "59.83 fb^{-1} (2018, 13 TeV)",
}
color_dict = {
    "QCD": (204, 204, 204),
    "diboson_J": (100, 192, 232),
    "diboson_L": (100, 212, 232),
    "Wjets": (137, 160, 44),
    "ttbar_J": (155, 152, 204),
    "ttbar_L": (175, 152, 204),
    "DYjets_J": (191, 34, 41),
    "DYjets_L": (211, 34, 41),
    "embedding": (255, 204, 0),
}
label_dict = {
    "QCD": "QCD multijet",
    "diboson_J": "Diboson (jet#rightarrow#tau_{h})",
    "diboson_L": "Diboson",
    "Wjets": "W+jets",
    "ttbar_J": "t#bar{t} (jet#rightarrow#tau_{h})",
    "ttbar_L": "t#bar{t}",
    "DYjets_J": "Z#rightarrow ll (jet#rightarrow#tau_{h})",
    "DYjets_L": "Z#rightarrow ll",
    "embedding": "#tau embedded",
    "data": "Data",
    "data_subtracted": "reduced Data",
}
var_dict = {"tau_pt": "p_{T}(#tau_{h}) (GeV)", "tau_phi": "#phi(#tau_{h})"}


def plot_FFs(ratio, config, process, njets):
    c = ROOT.TCanvas("c", "", 700, 700)
    c.SetRightMargin(0.05)
    c.SetLeftMargin(0.16)
    c.SetBottomMargin(0.12)

    ROOT.gStyle.SetOptStat(0)  # set off of the histogram statistics box
    ROOT.gStyle.SetTextFont(
        42
    )  # chosing font, see https://root.cern/root/html534/TAttText.html

    ratio.SetAxisRange(0, 0.5, "Y")
    ratio.SetMarkerStyle(20)
    ratio.SetMarkerSize(1.2)
    ratio.SetLineWidth(2)
    ratio.SetLineColor(ROOT.kBlack)
    ratio.SetTitle("")
    ratio.GetXaxis().SetMoreLogLabels()
    ratio.GetXaxis().SetNoExponent()
    ratio.GetYaxis().SetTitle(FF_YAxis[process])
    ratio.GetXaxis().SetTitle("p_{T}(#tau_{h}) (GeV)")
    ratio.GetYaxis().SetLabelSize(0.04)
    ratio.GetXaxis().SetLabelSize(0.04)
    ratio.GetYaxis().SetTitleSize(0.05)
    ratio.GetXaxis().SetTitleSize(0.05)

    fit = ratio.Fit(
        "pol1", "SF"
    )  # https://root.cern/doc/master/classTH1.html#a7e7d34c91d5ebab4fc9bba3ca47dabdd
    cl68 = fit.GetConfidenceIntervals(cl=0.68)
    chi2 = fit.Chi2()
    dof = fit.Ndf()
    print("Confidence intervals:", cl68)
    print("-" * 50)

    if config["target_process"][process]["tau_pt_logx"]:
        c.SetLogx()

    ratio.Draw()

    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextSize(0.04)
    text.DrawLatex(
        0.165,
        0.915,
        "channel: {}, {} jets".format(
            channel_dict[config["channel"]], njets.replace(" ", "")
        ),
    )
    text.SetTextSize(0.035)
    text.DrawLatex(0.6, 0.915, "{}".format(era_dict[config["era"]]))
    text.SetTextSize(0.035)
    text.DrawLatex(
        0.62, 0.8, "{} = {} / {}".format("#chi^{2} / N_{dof}", round(chi2, 2), dof)
    )

    func.check_output_path(os.getcwd() + "/workdir/" + config["workdir_name"])
    c.SaveAs(
        "workdir/{}/ff_{}_{}jets.png".format(
            config["workdir_name"], process, njets.replace(" ", "")
        )
    )
    c.SaveAs(
        "workdir/{}/ff_{}_{}jets.pdf".format(
            config["workdir_name"], process, njets.replace(" ", "")
        )
    )
    c.Close()


def plot_data_mc(hists, config, region, var, process, njets, sig, bkg):
    c = ROOT.TCanvas("c", "", 850, 700)
    c.SetRightMargin(0.05)
    c.SetLeftMargin(0.16)
    c.SetBottomMargin(0.12)

    ROOT.gStyle.SetOptStat(0)  # set off of the histogram statistics box
    ROOT.gStyle.SetTextFont(
        42
    )  # chosing font, see https://root.cern/root/html534/TAttText.html

    stack = ROOT.THStack()
    mc = hists[bkg[0]].Clone()
    for sample in bkg:
        hists[sample].SetLineWidth(1)
        hists[sample].SetFillStyle(1001)
        hists[sample].SetLineColor(ROOT.kBlack)
        color = color_dict[sample]
        hists[sample].SetFillColor(ROOT.TColor.GetColor(*color))
        stack.Add(hists[sample])
        if sample != bkg[0]:
            mc.Add(hists[sample])
    stack.SetTitle("")
    mc.SetFillStyle(3004)
    mc.SetFillColor(ROOT.kGray + 2)
    # mc.SetFillColorAlpha(ROOT.kGray+2, 0.35)
    stack.Draw("HIST")
    stack.GetYaxis().SetTitle("N_{Events}")
    stack.GetXaxis().SetTitle(var_dict[var])
    stack.GetYaxis().SetLabelSize(0.04)
    stack.GetXaxis().SetLabelSize(0.04)
    stack.GetYaxis().SetTitleSize(0.05)
    stack.GetXaxis().SetTitleSize(0.05)
    mc.Draw("E2 SAME")

    if config["target_process"][process]["tau_pt_logx"]:
        c.SetLogx()
    # c.SetLogy()

    data = hists[sig]
    stack.SetMaximum(data.GetMaximum() + 0.4 * data.GetMaximum())
    data.SetMarkerStyle(20)
    data.SetMarkerSize(1.2)
    data.SetLineWidth(2)
    data.SetLineColor(ROOT.kBlack)
    data.Draw("E SAME")

    legend = ROOT.TLegend(0.65, 0.47, 0.93, 0.88)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)
    legend.SetTextAlign(32)
    legend.AddEntry(data, label_dict[sig], "lep")
    for sample in bkg:
        legend.AddEntry(hists[sample], label_dict[sample], "f")
    legend.Draw("SAME")

    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextSize(0.04)
    text.DrawLatex(
        0.23,
        0.915,
        "channel: {}, {} jets".format(
            channel_dict[config["channel"]], njets.replace(" ", "")
        ),
    )
    text.SetTextSize(0.035)
    text.DrawLatex(0.66, 0.915, "{}".format(era_dict[config["era"]]))

    func.check_output_path(os.getcwd() + "/workdir/" + config["workdir_name"])
    c.SaveAs(
        "workdir/{}/hist_{}_{}_{}_{}jets.png".format(
            config["workdir_name"], sig, process, region, njets.replace(" ", "")
        )
    )
    c.SaveAs(
        "workdir/{}/hist_{}_{}_{}_{}jets.pdf".format(
            config["workdir_name"], sig, process, region, njets.replace(" ", "")
        )
    )
    c.Close()


def plot_data_mc_ratio(hists, config, region, var, process, njets, sig, bkg):

    ROOT.gStyle.SetOptStat(0)  # set off of the histogram statistics box
    ROOT.gStyle.SetTextFont(
        42
    )  # chosing font, see https://root.cern/root/html534/TAttText.html

    stack = ROOT.THStack()
    mc = hists[bkg[0]].Clone()
    for sample in bkg:
        hists[sample].SetLineWidth(1)
        hists[sample].SetFillStyle(1001)
        hists[sample].SetLineColor(ROOT.kBlack)
        color = color_dict[sample]
        hists[sample].SetFillColor(ROOT.TColor.GetColor(*color))
        stack.Add(hists[sample])
        if sample != bkg[0]:
            mc.Add(hists[sample])
    stack.SetTitle("")
    mc.SetFillStyle(3004)
    mc.SetFillColor(ROOT.kGray + 2)

    data = hists[sig]
    data.SetMarkerStyle(20)
    data.SetMarkerSize(1.2)
    data.SetLineWidth(2)
    data.SetLineColor(ROOT.kBlack)

    ratio = data.Clone("ratio")
    ratio.SetLineColor(ROOT.kBlack)
    ratio.SetMarkerStyle(20)
    ratio.SetTitle("")
    ratio.SetMinimum(0.75)
    ratio.SetMaximum(1.25)
    ratio.Divide(mc)

    # Adjust y-axis settings
    y = ratio.GetYaxis()
    y.SetTitle("#frac{data}{simulation}")
    y.SetNdivisions(505)
    y.SetTitleSize(0.10)
    y.SetTitleOffset(0.5)
    y.SetLabelSize(0.09)

    # Adjust x-axis settings
    x = ratio.GetXaxis()
    x.SetTitle(var_dict[var])
    x.SetTitleSize(0.12)
    x.SetLabelSize(0.1)

    c = ROOT.TCanvas("c", "", 700, 800)

    if config["target_process"][process]["tau_pt_logx"]:
        c.SetLogx()
    # c.SetLogy()

    # Upper histogram plot is pad1
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.25, 1, 1.0)
    pad1.SetRightMargin(0.05)
    pad1.SetLeftMargin(0.16)
    pad1.Draw()
    # Lower ratio plot is pad2
    c.cd()
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.01, 1, 0.3)
    pad2.SetBottomMargin(0.3)
    pad2.SetRightMargin(0.05)
    pad2.SetLeftMargin(0.16)
    pad2.SetGridy()
    pad2.Draw()

    # draw everything
    pad1.cd()
    stack.Draw("HIST")
    stack.GetYaxis().SetTitle("N_{Events}")
    stack.SetMaximum(data.GetMaximum() + 0.4 * data.GetMaximum())
    stack.GetYaxis().SetLabelSize(0.04)
    stack.GetXaxis().SetLabelSize(0.0)
    stack.GetYaxis().SetTitleSize(0.05)
    stack.GetXaxis().SetTitleSize(0.05)
    mc.Draw("E2 SAME")

    data.Draw("E SAME")

    legend = ROOT.TLegend(0.65, 0.47, 0.93, 0.88)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)
    legend.SetTextAlign(32)
    legend.AddEntry(data, label_dict[sig], "lep")
    for sample in bkg:
        legend.AddEntry(hists[sample], label_dict[sample], "f")
    legend.Draw("SAME")

    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextSize(0.04)
    text.DrawLatex(
        0.23,
        0.915,
        "channel: {}, {} jets".format(
            channel_dict[config["channel"]], njets.replace(" ", "")
        ),
    )
    text.SetTextSize(0.035)
    text.DrawLatex(0.65, 0.915, "{}".format(era_dict[config["era"]]))

    pad2.cd()
    ratio.Draw("ep")

    func.check_output_path(os.getcwd() + "/workdir/" + config["workdir_name"])
    c.SaveAs(
        "workdir/{}/hist_ratio_{}_{}_{}_{}jets.png".format(
            config["workdir_name"], sig, process, region, njets.replace(" ", "")
        )
    )
    c.SaveAs(
        "workdir/{}/hist_ratio_{}_{}_{}_{}jets.pdf".format(
            config["workdir_name"], sig, process, region, njets.replace(" ", "")
        )
    )
    c.Close()
