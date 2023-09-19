import ROOT
from io import StringIO
from wurlitzer import pipes, STDOUT
import logging
from typing import Any, Dict, List

import configs.general_definitions as gd


def plot_FFs(variable: str, ff_ratio: Any, uncertainties: Dict[str, Any], era: str, channel: str, process: str, category: Dict[str, str], output_path: str) -> None:
    '''
    Function which produces a fake factor plot.

    Args:
        variable: Name of the variable the fake factor is measured in
        ff_ratio: Histogram of the ratio (fake factor)
        uncertainties: Dictionary with graphs for each uncertainty/variation corresponding to "ff_ratio"
        era: Information about the era is added to the plot
        channel: Information about the channel is added to the plot
        process: Information about the process is added to the plot 
        category: Information about the category split is added to the plot 
        output_path: Path where the plot should be stored

    Return:
        None
    '''
    log = logging.getLogger("ff_calculation")

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
    ff_ratio.GetYaxis().SetTitle(gd.FF_YAxis[process])
    ff_ratio.GetXaxis().SetTitle(gd.variable_dict[channel][variable])
    ff_ratio.GetYaxis().SetLabelSize(0.04)
    ff_ratio.GetXaxis().SetLabelSize(0.04)
    ff_ratio.GetYaxis().SetTitleSize(0.05)
    ff_ratio.GetXaxis().SetTitleSize(0.05)

    ff_ratio.Draw()
    for unc in uncertainties:
        uncertainties[unc].SetLineWidth(2)
        uncertainties[unc].SetLineColor(gd.color_dict[unc])
        uncertainties[unc].SetFillColorAlpha(gd.color_dict[unc], 0.35)
        uncertainties[unc].Draw("E3 L SAME")

    legend = ROOT.TLegend(0.2, 0.68, 0.5, 0.88)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.03)
    legend.SetTextAlign(12)
    legend.AddEntry(ff_ratio, "measured", "lep")
    for unc in uncertainties:
        legend.AddEntry(uncertainties[unc], gd.label_dict[unc], "fl")
    legend.Draw("SAME")

    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextSize(0.03)
    text.DrawLatex(
        0.165,
        0.917,
        f"channel: {gd.channel_dict[channel]}, {', '.join([f'{gd.category_dict[var]} {category[var]}' for var in category.keys()])}",
    )
    text.SetTextSize(0.035)
    text.DrawLatex(0.6, 0.915, f"{gd.era_dict[era]}")

    out = StringIO()
    with pipes(stdout=out, stderr=STDOUT):
        c.SaveAs(
            f"{output_path}/ff_{process}_{'_'.join([f'{var}_{category[var]}' for var in category.keys()])}.png"
        )
        c.SaveAs(
            f"{output_path}/ff_{process}_{'_'.join([f'{var}_{category[var]}' for var in category.keys()])}.pdf"
        )
    log.info(out.getvalue())
    c.Close()


def plot_data_mc(variable: str, hists: Dict[str, Any], era: str, channel: str, process: str, region: str, data: str, samples: List[str], category: Dict[str, str], output_path: str) -> None:
    '''
    Function which produces a data to MC control plot.

    Args:
        variable: Name of the variable the fake factor is measured in
        hists: Dictionary with histogram for all processes and data
        era: Information about the era is added to the plot
        channel: Information about the channel is added to the plot
        process: Information about the process is added to the plot 
        region: Information about the fake factor calculation region is added to the plot name
        data: Name of the data process in "hists"
        samples: List of processes to be considered from "hists" for MC
        category: Information about the category split is added to the plot 
        output_path: Path where the plot should be stored
    
    Return:
        None
    '''
    log = logging.getLogger("ff_calculation")

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
        color = gd.color_dict[sample]
        hists[sample].SetFillColor(ROOT.TColor.GetColor(*color))
        stack.Add(hists[sample])
        if sample != samples[0]:
            mc.Add(hists[sample])
    stack.SetTitle("")
    mc.SetFillStyle(3004)
    mc.SetFillColor(ROOT.kGray + 2)
    stack.Draw("HIST")
    stack.GetYaxis().SetTitle("N_{Events}")
    stack.GetXaxis().SetTitle(gd.variable_dict[channel][variable])
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
    legend.AddEntry(obs, gd.label_dict[data], "lep")
    for sample in samples:
        legend.AddEntry(hists[sample], gd.label_dict[sample], "f")
    legend.Draw("SAME")

    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextSize(0.03)
    text.DrawLatex(
        0.23,
        0.917,
        f"channel: {gd.channel_dict[channel]}, {', '.join([f'{gd.category_dict[var]} {category[var]}' for var in category.keys()])}",
    )
    text.SetTextSize(0.035)
    text.DrawLatex(0.66, 0.915, "{}".format(gd.era_dict[era]))

    if data == "data_subtracted":
        hist_str = "reduced"
    else:
        hist_str = "full"

    out = StringIO()
    with pipes(stdout=out, stderr=STDOUT):
        c.SaveAs(
            f"{output_path}/hist_{hist_str}_{variable}_{process}_{region}_{'_'.join([f'{var}_{category[var]}' for var in category.keys()])}.png"
        )
        c.SaveAs(
            f"{output_path}/hist_{hist_str}_{variable}_{process}_{region}_{'_'.join([f'{var}_{category[var]}' for var in category.keys()])}.pdf"
        )
    log.info(out.getvalue())
    c.Close()


def plot_data_mc_ratio(variable: str, hists: Dict[str, Any], era: str, channel: str, process: str, region: str, data: str, samples: List[str], category: Dict[str, str], output_path: str) -> None:
    '''
    Function which produces a data to MC control plot with a ratio plot.

    Args:
        variable: Name of the variable the fake factor is measured in
        hists: Dictionary with histogram for all processes and data
        era: Information about the era is added to the plot
        channel: Information about the channel is added to the plot
        process: Information about the process is added to the plot 
        region: Information about the fake factor calculation region is added to the plot name
        data: Name of the data process in "hists"
        samples: List of processes to be considered from "hists" for MC
        category: Information about the category split is added to the plot 
        output_path: Path where the plot should be stored
    
    Return:
        None
    '''
    log = logging.getLogger("ff_calculation")

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
        color = gd.color_dict[sample]
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
    x.SetTitle(gd.variable_dict[channel][variable])
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
    legend.AddEntry(obs, gd.label_dict[data], "lep")
    for sample in samples:
        legend.AddEntry(hists[sample], gd.label_dict[sample], "f")
    legend.Draw("SAME")

    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextSize(0.03)
    text.DrawLatex(
        0.23,
        0.917,
        f"channel: {gd.channel_dict[channel]}, {', '.join([f'{gd.category_dict[var]} {category[var]}' for var in category.keys()])}",
    )
    text.SetTextSize(0.035)
    text.DrawLatex(0.65, 0.915, "{}".format(gd.era_dict[era]))

    pad2.cd()
    ratio.Draw("ep")

    if data == "data_subtracted":
        hist_str = "reduced"
    else:
        hist_str = "full"

    out = StringIO()
    with pipes(stdout=out, stderr=STDOUT):
        c.SaveAs(
            f"{output_path}/hist_ratio_{hist_str}_{variable}_{process}_{region}_{'_'.join([f'{var}_{category[var]}' for var in category.keys()])}.png"
        )
        c.SaveAs(
            f"{output_path}/hist_ratio_{hist_str}_{variable}_{process}_{region}_{'_'.join([f'{var}_{category[var]}' for var in category.keys()])}.pdf"
        )
    log.info(out.getvalue())
    c.Close()


def plot_fractions(variable: str, hists: Dict[str, Any], era: str, channel: str, region: str, processes: List[str], category: Dict[str, str], output_path: str) -> None:
    '''
    Function which produces a fraction plot where the sum of all considered process contributions 
    results to 1 for each bin in the fraction histogram.

    Args:
        variable: Name of the variable the fraction is measured in
        hists: Dictionary with fraction histograms for all considered "processes"
        era: Information about the era is added to the plot
        channel: Information about the channel is added to the plot
        region: Information about the fraction calculation region is added to the plot name
        processes: List of processes which are considered for the fraction calculation
        category: Information about the category split is added to the plot 
        output_path: Path where the plot should be stored
    
    Return:
        None
    '''
    log = logging.getLogger("ff_calculation")

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
    for sample in processes:
        hists[sample].SetLineWidth(1)
        hists[sample].SetFillStyle(1001)
        hists[sample].SetLineColor(ROOT.kBlack)
        color = gd.color_dict[sample]
        hists[sample].SetFillColor(ROOT.TColor.GetColor(*color))
        stack.Add(hists[sample])
    stack.SetTitle("")
    stack.Draw("HIST")
    stack.GetYaxis().SetTitle("Fraction")
    stack.GetXaxis().SetTitle(gd.variable_dict[channel][variable])
    stack.GetYaxis().SetLabelSize(0.04)
    stack.GetXaxis().SetLabelSize(0.04)
    stack.GetYaxis().SetTitleSize(0.05)
    stack.GetXaxis().SetTitleSize(0.05)

    legend = ROOT.TLegend(0.65, 0.47, 0.93, 0.88)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.035)
    legend.SetTextAlign(32)
    for sample in processes:
        legend.AddEntry(hists[sample], gd.label_dict[sample], "f")
    legend.Draw("SAME")

    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextSize(0.03)
    text.DrawLatex(
        0.23,
        0.915,
        f"channel: {gd.channel_dict[channel]}, {', '.join([f'{gd.category_dict[var]} {category[var]}' for var in category.keys()])}",
    )
    text.SetTextSize(0.035)
    text.DrawLatex(0.66, 0.915, f"{gd.era_dict[era]}")

    out = StringIO()
    with pipes(stdout=out, stderr=STDOUT):
        c.SaveAs(
            f"{output_path}/fraction_{variable}_{region}_{'_'.join([f'{var}_{category[var]}' for var in category.keys()])}.png"
        )
        c.SaveAs(
            f"{output_path}/fraction_{variable}_{region}_{'_'.join([f'{var}_{category[var]}' for var in category.keys()])}.pdf"
        )
    log.info(out.getvalue())
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
    if var not in ["pt_1", "pt_2", "mt_1", "boosted_pt_1", "boosted_pt_2", "boosted_mt_1", "iso_1", "boosted_iso_1"]:
        corr_ratio.GetXaxis().SetTitle(var_dict[var])
    else:
        corr_ratio.GetXaxis().SetTitle(var_dict[var][config["channel"]])
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
