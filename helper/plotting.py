import os
import ROOT
from io import StringIO
from wurlitzer import pipes, STDOUT

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
    "diboson_L": (100, 232, 232),
    "diboson_T": (100, 232, 202),
    "Wjets": (137, 160, 44),
    "ttbar_J": (155, 152, 204),
    "ttbar_L": (195, 152, 204),
    "ttbar_T": (195, 152, 174),
    "DYjets_J": (191, 34, 41),
    "DYjets_L": (231, 34, 41),
    "DYjets_T": (231, 34, 11),
    "ST_J": (184, 207, 92),
    "ST_L": (184, 237, 92),
    "ST_T": (214, 237, 92),
    "embedding": (255, 204, 0),
    "tau_fakes": (137, 160, 44),
    "data_ff": (255, 204, 0),
}
unc_color_dict = {
    "fit_graph_slope": ROOT.kBlue,
    "fit_graph_norm": ROOT.kRed,
    "fit_graph_mc_sub": ROOT.kGreen,
}

label_dict = {
    "QCD": "QCD multijet",
    "diboson_J": "Diboson (jet#rightarrow#tau_{h})",
    "diboson_L": "Diboson (lep#rightarrow#tau_{h})",
    "diboson_T": "Diboson (genuine #tau_{h})",
    "Wjets": "W+jets",
    "ttbar_J": "t#bar{t} (jet#rightarrow#tau_{h})",
    "ttbar_L": "t#bar{t} (lep#rightarrow#tau_{h})",
    "ttbar_T": "t#bar{t} (genuine #tau_{h})",
    "DYjets_J": "Z#rightarrow ll (jet#rightarrow#tau_{h})",
    "DYjets_L": "Z#rightarrow ll (lep#rightarrow#tau_{h})",
    "DYjets_T": "Z#rightarrow ll (genuine #tau_{h})",
    "ST_J": "ST (jet#rightarrow#tau_{h})",
    "ST_L": "ST (lep#rightarrow#tau_{h})",
    "ST_T": "ST (genuine #tau_{h})",
    "embedding": "#tau embedded",
    "data": "Data",
    "data_subtracted": "reduced Data",
    "data_ff": "Data with FFs",
    "tau_fakes": "jet#rightarrow#tau_{h}",
    "fit_graph_slope": "best fit (slope unc.)",
    "fit_graph_norm": "best fit (normalization unc.)",
    "fit_graph_mc_sub": "best fit (MC subtraction unc.)",
}
var_dict = {
    "pt_1": "p_{T}(e/#mu) (GeV)",
    "pt_2": "p_{T}(#tau_{h}) (GeV)",
    "eta_1": "#eta(e/#mu)",
    "eta_2": "#eta(#tau_{h})",
    "phi_1": "#phi(e/#mu)",
    "phi_2": "#phi(#tau_{h})",
    "iso_1": "isolation (e/#mu)",
    "mt_1": "m_{T}(e/#mu) (GeV)",
    "m_vis": "m_{vis} (GeV)",
    "mjj": "m_{jj} (GeV)",
    "met": "MET (GeV)",
    "metphi": "MET #phi",
    "bpt_1": "p_{T}(leading b) (GeV)",
    "beta_1": "#eta (leading b)",
    "bphi_1": "#phi (leading b)",
    "bpt_2": "p_{T}(second b) (GeV)",
    "njets": "number of jets",
    "nbtag": "number of b-tagged jets",
}

cat_dict = {
    "incl": "incl.",
    "njets": "jets",
    "nbtag": "b-jets",
    "deltaR_ditaupair": "#Delta" + "R(l#tau_{h})",
}


def plot_FFs(ff_ratio, uncertainties, process, config, split, save_path):
    ROOT.PyConfig.IgnoreCommandLineOptions = True
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    c = ROOT.TCanvas("c", "", 700, 700)
    c.SetRightMargin(0.05)
    c.SetLeftMargin(0.16)
    c.SetBottomMargin(0.12)

    ROOT.gStyle.SetOptStat(0)  # set off of the histogram statistics box
    ROOT.gStyle.SetTextFont(
        42
    )  # chosing font, see https://root.cern/root/html534/TAttText.html

    ff_ratio.SetAxisRange(0, 0.5, "Y")
    ff_ratio.SetMarkerStyle(20)
    ff_ratio.SetMarkerSize(1.2)
    ff_ratio.SetLineWidth(2)
    ff_ratio.SetLineColor(ROOT.kBlack)
    ff_ratio.SetTitle("")
    ff_ratio.GetXaxis().SetMoreLogLabels()
    ff_ratio.GetXaxis().SetNoExponent()
    ff_ratio.GetYaxis().SetTitle(FF_YAxis[process])
    ff_ratio.GetXaxis().SetTitle(
        var_dict[config["target_process"][process]["var_dependence"]]
    )
    ff_ratio.GetYaxis().SetLabelSize(0.04)
    ff_ratio.GetXaxis().SetLabelSize(0.04)
    ff_ratio.GetYaxis().SetTitleSize(0.05)
    ff_ratio.GetXaxis().SetTitleSize(0.05)

    ff_ratio.Draw()
    for unc in uncertainties:
        uncertainties[unc].SetLineWidth(2)
        uncertainties[unc].SetLineColor(unc_color_dict[unc])
        uncertainties[unc].SetFillColorAlpha(unc_color_dict[unc], 0.35)
        uncertainties[unc].Draw("E3 L SAME")

    legend = ROOT.TLegend(0.2, 0.68, 0.5, 0.88)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.03)
    legend.SetTextAlign(12)
    legend.AddEntry(ff_ratio, "measured", "lep")
    for unc in uncertainties:
        legend.AddEntry(uncertainties[unc], label_dict[unc], "fl")
    legend.Draw("SAME")

    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextSize(0.03)
    text.DrawLatex(
        0.165,
        0.917,
        "channel: {}, {}".format(
            channel_dict[config["channel"]],
            ", ".join(
                ["{} {}".format(split[var], cat_dict[var]) for var in split.keys()]
            ),
        ),
    )
    text.SetTextSize(0.035)
    text.DrawLatex(0.6, 0.915, "{}".format(era_dict[config["era"]]))
    # text.SetTextSize(0.035)
    # text.DrawLatex(
    #     0.62, 0.8, "{} = {} / {}".format("#chi^{2} / N_{dof}", round(chi2, 2), dof)
    # )

    out = StringIO()
    with pipes(stdout=out, stderr=STDOUT):
        c.SaveAs(
            "{}/ff_{}_{}.png".format(
                save_path,
                process,
                "_".join(["{}_{}".format(split[var], var) for var in split.keys()]),
            )
        )
        c.SaveAs(
            "{}/ff_{}_{}.pdf".format(
                save_path,
                process,
                "_".join(["{}_{}".format(split[var], var) for var in split.keys()]),
            )
        )
    print(out.getvalue())
    c.Close()


def plot_data_mc(hists, config, var, process, region, data, samples, split, save_path):
    ROOT.PyConfig.IgnoreCommandLineOptions = True
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    c = ROOT.TCanvas("c", "", 850, 700)
    c.SetRightMargin(0.05)
    c.SetLeftMargin(0.16)
    c.SetBottomMargin(0.12)

    ROOT.gStyle.SetOptStat(0)  # set off of the histogram statistics box
    ROOT.gStyle.SetTextFont(
        42
    )  # chosing font, see https://root.cern/root/html534/TAttText.html

    stack = ROOT.THStack()
    mc = hists[samples[0]].Clone()
    for sample in samples:
        hists[sample].SetLineWidth(1)
        hists[sample].SetFillStyle(1001)
        hists[sample].SetLineColor(ROOT.kBlack)
        color = color_dict[sample]
        hists[sample].SetFillColor(ROOT.TColor.GetColor(*color))
        stack.Add(hists[sample])
        if sample != samples[0]:
            mc.Add(hists[sample])
    stack.SetTitle("")
    mc.SetFillStyle(3004)
    mc.SetFillColor(ROOT.kGray + 2)
    stack.Draw("HIST")
    stack.GetYaxis().SetTitle("N_{Events}")
    stack.GetXaxis().SetTitle(var_dict[var])
    stack.GetYaxis().SetLabelSize(0.04)
    stack.GetXaxis().SetLabelSize(0.04)
    stack.GetYaxis().SetTitleSize(0.05)
    stack.GetXaxis().SetTitleSize(0.05)
    mc.Draw("E2 SAME")

    obs = hists[data]
    stack.SetMaximum(
        max(
            obs.GetMaximum() + 0.4 * obs.GetMaximum(),
            stack.GetMaximum() + 0.4 * stack.GetMaximum(),
        )
    )
    obs.SetMarkerStyle(20)
    obs.SetMarkerSize(1.2)
    obs.SetLineWidth(2)
    obs.SetLineColor(ROOT.kBlack)
    obs.Draw("E SAME")

    legend = ROOT.TLegend(0.65, 0.47, 0.93, 0.88)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)
    legend.SetTextAlign(12)
    legend.AddEntry(obs, label_dict[data], "lep")
    for sample in samples:
        legend.AddEntry(hists[sample], label_dict[sample], "f")
    legend.Draw("SAME")

    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextSize(0.03)
    text.DrawLatex(
        0.23,
        0.917,
        "channel: {}, {}".format(
            channel_dict[config["channel"]],
            ", ".join(
                ["{} {}".format(split[var], cat_dict[var]) for var in split.keys()]
            ),
        ),
    )
    text.SetTextSize(0.035)
    text.DrawLatex(0.66, 0.915, "{}".format(era_dict[config["era"]]))

    if data == "data_subtracted":
        hist_str = "reduced"
    else:
        hist_str = "full"

    out = StringIO()
    with pipes(stdout=out, stderr=STDOUT):
        c.SaveAs(
            "{}/hist_{}_{}_{}_{}_{}.png".format(
                save_path,
                hist_str,
                var,
                process,
                region,
                "_".join(["{}_{}".format(split[var], var) for var in split.keys()]),
            )
        )
        c.SaveAs(
            "{}/hist_{}_{}_{}_{}_{}.pdf".format(
                save_path,
                hist_str,
                var,
                process,
                region,
                "_".join(["{}_{}".format(split[var], var) for var in split.keys()]),
            )
        )
    print(out.getvalue())
    c.Close()


def plot_data_mc_ratio(
    hists, config, var, process, region, data, samples, split, save_path
):
    ROOT.PyConfig.IgnoreCommandLineOptions = True
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    ROOT.gStyle.SetOptStat(0)  # set off of the histogram statistics box
    ROOT.gStyle.SetTextFont(
        42
    )  # chosing font, see https://root.cern/root/html534/TAttText.html
    stack = ROOT.THStack()
    mc = hists[samples[0]].Clone()
    for sample in samples:
        hists[sample].SetLineWidth(1)
        hists[sample].SetFillStyle(1001)
        hists[sample].SetLineColor(ROOT.kBlack)
        color = color_dict[sample]
        hists[sample].SetFillColor(ROOT.TColor.GetColor(*color))
        stack.Add(hists[sample])
        if sample != samples[0]:
            mc.Add(hists[sample])
    stack.SetTitle("")
    mc.SetFillStyle(3004)
    mc.SetFillColor(ROOT.kGray + 2)

    obs = hists[data]
    obs.SetMarkerStyle(20)
    obs.SetMarkerSize(1.2)
    obs.SetLineWidth(2)
    obs.SetLineColor(ROOT.kBlack)

    ratio = obs.Clone("ratio")
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
    stack.SetMaximum(obs.GetMaximum() + 0.4 * obs.GetMaximum())
    stack.GetYaxis().SetLabelSize(0.04)
    stack.GetXaxis().SetLabelSize(0.0)
    stack.GetYaxis().SetTitleSize(0.05)
    stack.GetXaxis().SetTitleSize(0.05)
    mc.Draw("E2 SAME")

    obs.Draw("E SAME")

    legend = ROOT.TLegend(0.65, 0.47, 0.93, 0.88)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)
    legend.SetTextAlign(12)
    legend.AddEntry(obs, label_dict[data], "lep")
    for sample in samples:
        legend.AddEntry(hists[sample], label_dict[sample], "f")
    legend.Draw("SAME")

    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextSize(0.03)
    text.DrawLatex(
        0.23,
        0.917,
        "channel: {}, {}".format(
            channel_dict[config["channel"]],
            ", ".join(
                ["{} {}".format(split[var], cat_dict[var]) for var in split.keys()]
            ),
        ),
    )
    text.SetTextSize(0.035)
    text.DrawLatex(0.65, 0.915, "{}".format(era_dict[config["era"]]))

    pad2.cd()
    ratio.Draw("ep")

    if data == "data_subtracted":
        hist_str = "reduced"
    else:
        hist_str = "full"

    out = StringIO()
    with pipes(stdout=out, stderr=STDOUT):
        c.SaveAs(
            "{}/hist_ratio_{}_{}_{}_{}_{}.png".format(
                save_path,
                hist_str,
                var,
                process,
                region,
                "_".join(["{}_{}".format(split[var], var) for var in split.keys()]),
            )
        )
        c.SaveAs(
            "{}/hist_ratio_{}_{}_{}_{}_{}.pdf".format(
                save_path,
                hist_str,
                var,
                process,
                region,
                "_".join(["{}_{}".format(split[var], var) for var in split.keys()]),
            )
        )
    print(out.getvalue())
    c.Close()


def fraction_plot(hists, config, var, region, samples, split, save_path):
    ROOT.PyConfig.IgnoreCommandLineOptions = True
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    c = ROOT.TCanvas("c", "", 750, 700)
    c.SetRightMargin(0.05)
    c.SetLeftMargin(0.16)
    c.SetBottomMargin(0.12)

    ROOT.gStyle.SetOptStat(0)  # set off of the histogram statistics box
    ROOT.gStyle.SetTextFont(
        42
    )  # chosing font, see https://root.cern/root/html534/TAttText.html

    stack = ROOT.THStack()
    for sample in samples:
        hists[sample].SetLineWidth(1)
        hists[sample].SetFillStyle(1001)
        hists[sample].SetLineColor(ROOT.kBlack)
        color = color_dict[sample]
        hists[sample].SetFillColor(ROOT.TColor.GetColor(*color))
        stack.Add(hists[sample])
    stack.SetTitle("")
    stack.Draw("HIST")
    stack.GetYaxis().SetTitle("Fraction")
    stack.GetXaxis().SetTitle(var_dict[var])
    stack.GetYaxis().SetLabelSize(0.04)
    stack.GetXaxis().SetLabelSize(0.04)
    stack.GetYaxis().SetTitleSize(0.05)
    stack.GetXaxis().SetTitleSize(0.05)

    legend = ROOT.TLegend(0.65, 0.47, 0.93, 0.88)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)
    legend.SetTextAlign(32)
    for sample in samples:
        legend.AddEntry(hists[sample], label_dict[sample], "f")
    legend.Draw("SAME")

    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextSize(0.03)
    text.DrawLatex(
        0.23,
        0.915,
        "channel: {}, {}".format(
            channel_dict[config["channel"]],
            ", ".join(
                ["{} {}".format(split[var], cat_dict[var]) for var in split.keys()]
            ),
        ),
    )
    text.SetTextSize(0.035)
    text.DrawLatex(0.66, 0.915, "{}".format(era_dict[config["era"]]))

    out = StringIO()
    with pipes(stdout=out, stderr=STDOUT):
        c.SaveAs(
            "{}/fraction_{}_{}_{}.png".format(
                save_path,
                var,
                region,
                "_".join(["{}_{}".format(split[var], var) for var in split.keys()]),
            )
        )
        c.SaveAs(
            "{}/fraction_{}_{}_{}.pdf".format(
                save_path,
                var,
                region,
                "_".join(["{}_{}".format(split[var], var) for var in split.keys()]),
            )
        )
    print(out.getvalue())
    c.Close()


def plot_correction(
    corr_ratio, uncertainty, var, process, correction, config, save_path
):
    ROOT.PyConfig.IgnoreCommandLineOptions = True
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    c = ROOT.TCanvas("can", "", 700, 700)
    c.SetRightMargin(0.05)
    c.SetLeftMargin(0.16)
    c.SetBottomMargin(0.12)

    ROOT.gStyle.SetOptStat(0)  # set off of the histogram statistics box
    ROOT.gStyle.SetTextFont(
        42
    )  # chosing font, see https://root.cern/root/html534/TAttText.html

    corr_ratio.SetAxisRange(0, 2, "Y")
    corr_ratio.SetMarkerStyle(20)
    corr_ratio.SetMarkerSize(1.2)
    corr_ratio.SetLineWidth(2)
    corr_ratio.SetLineColor(ROOT.kBlack)
    corr_ratio.SetTitle("")
    corr_ratio.GetXaxis().SetMoreLogLabels()
    corr_ratio.GetXaxis().SetNoExponent()
    corr_ratio.GetYaxis().SetTitle("Correction")
    corr_ratio.GetXaxis().SetTitle(var_dict[var])
    corr_ratio.GetYaxis().SetLabelSize(0.04)
    corr_ratio.GetXaxis().SetLabelSize(0.04)
    corr_ratio.GetYaxis().SetTitleSize(0.05)
    corr_ratio.GetXaxis().SetTitleSize(0.05)

    corr_ratio.Draw()

    uncertainty.SetLineWidth(2)
    uncertainty.SetLineColor(ROOT.kOrange)
    uncertainty.SetFillColorAlpha(ROOT.kOrange, 0.35)

    uncertainty.Draw("E3 L SAME")

    legend = ROOT.TLegend(0.65, 0.68, 0.9, 0.88)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.03)
    legend.SetTextAlign(12)
    legend.AddEntry(corr_ratio, "measured", "lep")
    legend.AddEntry(uncertainty, "smoothed curve", "fl")

    legend.Draw("SAME")

    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextSize(0.03)
    text.DrawLatex(
        0.165,
        0.917,
        "channel: {}, {}".format(channel_dict[config["channel"]], process),
    )
    text.SetTextSize(0.035)
    text.DrawLatex(0.6, 0.915, "{}".format(era_dict[config["era"]]))
    text.SetTextSize(0.035)
    if "non_closure" in correction:
        text.DrawLatex(0.2, 0.8, "non closure correction")
    elif "DR_SR" in correction:
        text.DrawLatex(0.2, 0.8, "DR to SR correction")

    out = StringIO()
    with pipes(stdout=out, stderr=STDOUT):
        c.SaveAs("{}/corr_{}_{}.png".format(save_path, process, correction))
        c.SaveAs("{}/corr_{}_{}.pdf".format(save_path, process, correction))
    print(out.getvalue())
    c.Close()
