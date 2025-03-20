import array
import logging
import os
from io import StringIO
from typing import Any, Dict, Iterable, List, Tuple, Union

import ROOT
from wurlitzer import STDOUT, pipes

import configs.general_definitions as gd


def save_plots_and_data(
    output_path: str,
    plot_filename_and_obj: Union[Tuple[str, Any], None] = None,
    data_filename_and_obj: Union[Tuple[str, Any], None] = None,
    plot_extensions: Iterable[str] = ("png", "pdf"),
    data_extensions: Iterable[str] = ("root",),
    logger: Union[logging.Logger, None] = None,
) -> None:
    """
    Save plot and/or data ROOT objects to disk in designated subdirectories within the output path.

    This function creates subdirectories for each file extension (e.g. "png", "pdf", "root") in the
    provided output path if they do not already exist and saves the corresponding object there.
    Plot objects that currently work are TCanvas objects, while data objects can be either a ROOT object
    (e.g. a histogram) or a dictionary of ROOT objects (e.g. a dictionary of ROOT histograms).

    Args:
        output_path (str): Base directory where output subfolders (e.g. "png", "pdf", "root") will be created.
        plot_filename_and_obj (Optional[Tuple[str, Any]]):
            A tuple containing the base filename (without extension) and the plot object (e.g. a TCanvas) to be saved.
            If None, no plot will be saved.
        data_filename_and_obj (Optional[Tuple[str, Any]]):
            A tuple containing the base filename (without extension) and the data object
        plot_extensions (Iterable[str]): Iterable of file extensions for plot outputs. Default is ("png", "pdf").
        data_extensions (Iterable[str]): Iterable of file extensions for data outputs. Default is ("root",).
        logger (Optional[logging.Logger]): Logger object to which status messages will be written; if None, logging is skipped.

    Returns:
        None
    """
    base_filename_plot, plot_obj = "", None
    base_filename_data, data_obj = "", None

    def save_root_canvas(ext):
        plot_obj.SaveAs(os.path.join(output_path, ext, f"{base_filename_plot}.{ext}"))

    def save_root_hists(ext):
        root_filepath = os.path.join(output_path, ext, f"{base_filename_data}.{ext}")
        root_file = ROOT.TFile(root_filepath, "RECREATE")
        if isinstance(data_obj, dict):
            for key, hist in data_obj.items():
                hist.Write(key)
        elif isinstance(data_obj, ROOT.TH1D):
            data_obj.Write(base_filename_data)
        root_file.Close()

    tasks = []
    if plot_filename_and_obj is not None:
        base_filename_plot, plot_obj = plot_filename_and_obj
        tasks.append((plot_extensions, save_root_canvas))
    if data_filename_and_obj is not None:
        base_filename_data, data_obj = data_filename_and_obj
        tasks.append((data_extensions, save_root_hists))

    for extensions, _func in tasks:
        for extension in extensions:
            os.makedirs(os.path.join(output_path, extension), exist_ok=True)
            out = StringIO()
            with pipes(stdout=out, stderr=STDOUT):
                _func(extension)
            if logger:
                logger.info(out.getvalue())


def set_optional_TGraph_style_and_draw(obj: ROOT.TGraphAsymmErrors) -> None:
    """
    Function to set the style of a TGraphAsymmErrors object and draw it.

    Args:
        obj: TGraphAsymmErrors object to be styled and drawn.
    """
    n = obj.GetN()
    if n > 0:
        x0, y0 = array.array("d", [0.]), array.array("d", [0.])
        xn, yn = array.array("d", [0.]), array.array("d", [0.])
        obj.GetPoint(0, x0, y0)
        obj.GetPoint(n - 1, xn, yn)
        # Set x-axis limits: first point minus lower error, last point plus upper error.
        x_low = x0[0] - obj.GetErrorXlow(0)
        x_high = xn[0] + obj.GetErrorXhigh(n - 1)
        obj.GetXaxis().SetLimits(x_low, x_high)
    obj.SetLineStyle(0)  # Disable connecting lines
    obj.Draw("AP")  # Plot the axis and markers (and errors)


def plot_FFs(
    variable: str,
    ff_ratio: Any,
    uncertainties: Dict[str, Any],
    era: str,
    channel: str,
    process: str,
    category: Dict[str, str],
    output_path: str,
    logger: str,
    draw_option: str,
    save_data: bool = False,
) -> None:
    """
    Function which produces a fake factor plot.

    Args:
        variable: Name of the variable the fake factor is measured in
        ff_ratio: Histogram of the ratio (fake factor) or TGraph
        uncertainties: Dictionary with graphs for each uncertainty/variation corresponding to "ff_ratio"
        era: Information about the era is added to the plot
        channel: Information about the channel is added to the plot
        process: Information about the process is added to the plot
        category: Information about the category split is added to the plot
        output_path: Path where the plot should be stored
        logger: Name of the logger that should be used
        draw_option: Based on chosen fit_option to correctly plot a fitted function or a histogram as measurement
        save_data: Boolean to additionally save the histograms to a root file

    Return:
        None
    """
    log = logging.getLogger(logger)

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

    ff_ratio.SetMaximum(max(ff_ratio.GetMaximum() + 0.4 * ff_ratio.GetMaximum(), 0.5))
    ff_ratio.SetMinimum(0.0)
    ff_ratio.SetMarkerStyle(20)
    ff_ratio.SetMarkerSize(1.2)
    ff_ratio.SetLineWidth(2)
    ff_ratio.SetLineColor(ROOT.kBlack)
    ff_ratio.SetTitle("")
    ff_ratio.GetXaxis().SetMoreLogLabels()
    ff_ratio.GetXaxis().SetNoExponent()
    ff_ratio.GetYaxis().SetTitle(gd.FF_YAxis[process])
    ff_ratio.GetXaxis().SetTitle(
        gd.variable_dict[channel].get(variable, variable)
    )
    ff_ratio.GetYaxis().SetLabelSize(0.04)
    ff_ratio.GetXaxis().SetLabelSize(0.04)
    ff_ratio.GetYaxis().SetTitleSize(0.05)
    ff_ratio.GetXaxis().SetTitleSize(0.05)

    if isinstance(ff_ratio, ROOT.TGraphAsymmErrors):
        set_optional_TGraph_style_and_draw(ff_ratio)
    else:
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
        legend.AddEntry(uncertainties[unc], gd.label_dict[unc][draw_option], "fl")
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

    selection = '_'.join([f'{var}_{category[var]}' for var in category.keys()])
    base_filename = f"ff_{variable}_{process}_{selection}"
    save_plots_and_data(
        output_path=output_path,
        plot_filename_and_obj=(f"{base_filename}", c),
        data_filename_and_obj=(f"data_{base_filename}", {"ff_ratio": ff_ratio, **uncertainties}) if save_data else None,
        logger=log,
    )
    c.Close()


def plot_data_mc(
    variable: str,
    hists: Dict[str, Any],
    era: str,
    channel: str,
    process: str,
    region: str,
    data: str,
    samples: List[str],
    category: Dict[str, str],
    output_path: str,
    logger: str,
    yscale: str = "linear",
    save_data: bool = False,
) -> None:
    """
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
        logger: Name of the logger that should be used
        yscale: String do set the scaling of the y-axis, options are "linear" or "log"
        save_data: Boolean to additionally save the histograms to a root file

    Return:
        None
    """
    log = logging.getLogger(logger)

    ROOT.PyConfig.IgnoreCommandLineOptions = True
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    c = ROOT.TCanvas("c", "", 850, 700)
    c.SetRightMargin(0.05)
    c.SetLeftMargin(0.16)
    c.SetBottomMargin(0.12)

    if yscale == "log":
        c.SetLogy()

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
    stack.GetXaxis().SetTitle(
        gd.variable_dict[channel].get(variable, variable)
    )
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

    hist_str = "reduced" if data == "data_subtracted" else "full"
    selection = '_'.join([f'{var}_{category[var]}' for var in category.keys()])
    base_filename = f"hist_{hist_str}_{variable}_{process}_{region}_{selection}"
    save_plots_and_data(
        output_path=output_path,
        plot_filename_and_obj=(f"{base_filename}_{yscale}", c),
        data_filename_and_obj=(f"data_{base_filename}", hists) if save_data else None,
        logger=log,
    )
    c.Close()


def plot_data_mc_ratio(
    variable: str,
    hists: Dict[str, Any],
    era: str,
    channel: str,
    process: str,
    region: str,
    data: str,
    samples: List[str],
    category: Dict[str, str],
    output_path: str,
    logger: str,
    yscale: str = "linear",
    save_data: bool = False,
) -> None:
    """
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
        logger: Name of the logger that should be used
        yscale: String do set the scaling of the y-axis, options are "linear" or "log"
        save_data: Boolean to additionally save the histograms to a root file

    Return:
        None
    """
    log = logging.getLogger(logger)

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
    x.SetTitle(
        gd.variable_dict[channel].get(variable, variable)
    )
    x.SetTitleSize(0.12)
    x.SetLabelSize(0.1)

    c = ROOT.TCanvas("c", "", 700, 800)

    # Upper histogram plot is pad1
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.25, 1, 1.0)
    pad1.SetRightMargin(0.05)
    pad1.SetLeftMargin(0.16)

    if yscale == "log":
        pad1.SetLogy()

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

    if yscale == "log":
        c.SetLogy()

    hist_str = "reduced" if data == "data_subtracted" else "full"
    selection = '_'.join([f'{var}_{category[var]}' for var in category.keys()])
    base_filename = f"hist_ratio_{hist_str}_{variable}_{process}_{region}_{selection}"
    save_plots_and_data(
        output_path=output_path,
        plot_filename_and_obj=(f"{base_filename}_{yscale}", c),
        data_filename_and_obj=(f"data_{base_filename}", hists) if save_data else None,
        logger=log,
    )
    c.Close()


def plot_fractions(
    variable: str,
    hists: Dict[str, Any],
    era: str,
    channel: str,
    region: str,
    fraction_name: str,
    processes: List[str],
    category: Dict[str, str],
    output_path: str,
    logger: str,
    save_data: bool = False,
) -> None:
    """
    Function which produces a fraction plot where the sum of all considered process contributions
    results to 1 for each bin in the fraction histogram.

    Args:
        variable: Name of the variable the fraction is measured in
        hists: Dictionary with fraction histograms for all considered "processes"
        era: Information about the era is added to the plot
        channel: Information about the channel is added to the plot
        region: Information about the fraction calculation region is added to the plot name
        fraction_name: Name of the fraction, especially for tt channel relevant where 2 different fractions can be calculated
        processes: List of processes which are considered for the fraction calculation
        category: Information about the category split is added to the plot
        output_path: Path where the plot should be stored
        logger: Name of the logger that should be used
        save_data: Boolean to additionally save the histograms to a root file

    Return:
        None
    """
    log = logging.getLogger(logger)

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
    stack.GetXaxis().SetTitle(
        gd.variable_dict[channel].get(variable, variable)
    )
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

    selection = '_'.join([f'{var}_{category[var]}' for var in category.keys()])
    base_filename = f"fraction_{fraction_name}_{variable}_{region}_{selection}"
    save_plots_and_data(
        output_path=output_path,
        plot_filename_and_obj=(f"{base_filename}", c),
        data_filename_and_obj=(f"data_{base_filename}", hists) if save_data else None,
        logger=log,
    )
    c.Close()


def plot_correction(
    variable: str,
    corr_hist: Any,
    corr_graph: Any,
    corr_name: str,
    era: str,
    channel: str,
    process: str,
    output_path: str,
    logger: str,
    category: Dict[str, str] = None,
    save_data: bool = False,
) -> None:
    """
    Function which produces a plot of the correction including the correction histogram
    and the TGraph (with variations) from a smoothed fit which is derived from the histogram.

    Args:
        variable: Name of the variable the correction is measured in
        corr_hist: Histogram of the measured correction
        corr_graph: TGraph of the measured correction which is a smoothed function of "corr_hist" (includes variation)
        corr_name: Name of the correction
        era: Information about the era is added to the plot
        channel: Information about the channel is added to the plot
        process: Name of the target process the correction is calculated for
        output_path: Path where the plot should be stored
        logger: Name of the logger that should be used
        category: Information about the category split is added to the plot
        save_data: Boolean to additionally save the histograms to a root file

    Return:
        None
    """
    log = logging.getLogger(logger)

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

    # corr_hist.SetAxisRange(0, 2, "Y")
    corr_hist.SetMaximum(
        max(corr_hist.GetMaximum() + 0.4 * corr_hist.GetMaximum(), 2.0)
    )
    corr_hist.SetMinimum(0.0)
    corr_hist.SetMarkerStyle(20)
    corr_hist.SetMarkerSize(1.2)
    corr_hist.SetLineWidth(2)
    corr_hist.SetLineColor(ROOT.kBlack)
    corr_hist.SetTitle("")
    corr_hist.GetXaxis().SetMoreLogLabels()
    corr_hist.GetXaxis().SetNoExponent()
    corr_hist.GetYaxis().SetTitle("Correction")
    corr_hist.GetXaxis().SetTitle(
        gd.variable_dict[channel].get(variable, variable)
    )
    corr_hist.GetYaxis().SetLabelSize(0.04)
    corr_hist.GetXaxis().SetLabelSize(0.04)
    corr_hist.GetYaxis().SetTitleSize(0.05)
    corr_hist.GetXaxis().SetTitleSize(0.05)

    if isinstance(corr_hist, ROOT.TGraphAsymmErrors):
        set_optional_TGraph_style_and_draw(corr_hist)
    else:
        corr_hist.Draw()

    corr_graph.SetLineWidth(2)
    corr_graph.SetLineColor(ROOT.kOrange)
    corr_graph.SetFillColorAlpha(ROOT.kOrange, 0.35)

    corr_graph.Draw("E3 L SAME")

    legend = ROOT.TLegend(0.65, 0.68, 0.9, 0.88)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.03)
    legend.SetTextAlign(12)
    legend.AddEntry(corr_hist, "measured", "lep")
    legend.AddEntry(corr_graph, "smoothed curve", "fl")

    legend.Draw("SAME")

    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextSize(0.03)
    if category is not None:
        plot_text = f"channel: {gd.channel_dict[channel]}, {process}, {', '.join([f'{gd.category_dict[var]} {category[var]}' for var in category.keys()])}"
    else:
        plot_text = f"channel: {gd.channel_dict[channel]}, {process}"
    text.DrawLatex(
        0.165,
        0.917,
        plot_text,
    )
    text.SetTextSize(0.035)
    text.DrawLatex(0.6, 0.915, "{}".format(gd.era_dict[era]))
    text.SetTextSize(0.035)
    if "non_closure" in corr_name:
        text.DrawLatex(0.2, 0.8, "non closure correction")
    elif "DR_SR" in corr_name:
        text.DrawLatex(0.2, 0.8, "DR to SR correction")

    selection = ""
    if category is not None:
        selection = '_'.join([""] + [f'{var}_{category[var]}' for var in category.keys()])

    save_plots_and_data(
        output_path=output_path,
        plot_filename_and_obj=(f"corr_{process}_{corr_name}{selection}", c),
        data_filename_and_obj=(
            f"data_corr_{process}_{corr_name}{selection}",
            {"corr_hist": corr_hist, "corr_graph": corr_graph},
        ) if save_data else None,
        logger=log,
    )
    c.Close()
