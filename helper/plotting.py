import array
import copy
import logging
import os
import pickle
from copy import deepcopy
from functools import partial
from typing import Any, Callable, Dict, Iterable, List, Tuple, Union

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import mplhep as hep
import numpy as np
import ROOT
from matplotlib.legend_handler import HandlerBase, HandlerTuple
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import configs.general_definitions as gd
import helper.ff_functions as func

hep.style.use(hep.style.CMS)


class VTuple(tuple):
    """Vertical tuple: render items stacked (two lines like '=')."""
    pass


class MergedTuple(tuple):
    """Merged tuple: for combined handles like patch + line, using HandlerTuple."""
    pass


class HandlerEquals(HandlerBase):
    """
    Construct handler to render "=" inside legends with adjustable properties for each line.
    Supported Propreties are:
        lw: linewidth
        ls: linestyle
        c: color
    """
    def __init__(self, *args, spacing_factor=0.5, **kwargs):
        super().__init__(*args, **kwargs)
        self.spacing_factor = spacing_factor

    def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
        def get_props(item):
            if isinstance(item, dict):
                return item
            elif isinstance(item, Line2D):
                return {
                    'c': item.get_color(),
                    'lw': item.get_linewidth(),
                    'ls': item.get_linestyle()
                }
            else:
                raise ValueError(f"Unsupported type in handle: {type(item)}")

        if len(orig_handle) != 2:
            raise ValueError("HandlerEquals expects a tuple of exactly two items (dicts or Line2D).")

        props1, props2 = [get_props(p) for p in orig_handle]

        mid, offset = height / 2, (height * self.spacing_factor) / 2
        y_top, y_bottom = ydescent + mid + offset, ydescent + mid - offset

        line1 = Line2D([xdescent, xdescent + width], [y_top, y_top])
        line1.set(**props1)

        line2 = Line2D([xdescent, xdescent + width], [y_bottom, y_bottom])
        line2.set(**props2)

        return [line1, line2]


class SmartSciFormatter(ticker.Formatter):
    """
    A custom formatter for axis tick labels to format numbers using scientific notation.
    Adapts based on axis scale (logarithmic or linear) and tick values range.

    Converts into a LaTeX implementation with a mantissa and exponent: mantissa · 10^(exponent)
    """
    def __init__(self, axis, use_dot=True):
        self.axis = axis
        self.use_dot = use_dot

    def __call__(self, x, pos=None):
        if np.isclose(float(x), 0.0):
            return "0"

        ticks = [tick for tick in self.axis.get_ticklocs() if tick != 0.0]
        mantissa, exponent = f"{x:.1e}".split("e")
        mantissa = mantissa.rstrip("0").rstrip(".")
        dot = r"\cdot" if self.use_dot else r"\times"

        if self.axis.get_scale() == "log":
            mantissas = [float(f"{tick:.1e}".split("e")[0]) for tick in ticks if tick != 0.0]
            if all(np.isclose(m, 1.0) for m in mantissas):
                return rf"$10^{{{int(exponent)}}}$"
            else:
                return rf"${mantissa}{dot}10^{{{int(exponent)}}}$"

        # this is for linear scale
        if all(abs(t) > 1e3 or abs(t) < 1e-2 for t in ticks):
            return rf"${mantissa}{dot}10^{{{int(exponent)}}}$"
        else:
            return f"{x}".rstrip("0").rstrip(".")


HANDLER_MAP = {MergedTuple: HandlerTuple(ndivide=1, pad=0.0), VTuple: HandlerEquals()}


def save_plots_and_data(
    output_path: str,
    plot_filename_and_obj: Union[Tuple[str, Any], None] = None,
    data_filename_and_obj: Union[Tuple[str, Any], None] = None,
    plot_extensions: Iterable[str] = ("png", "pdf"),
    data_extensions: Iterable[str] = ("pickle",),
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

    def save_canvas(ext):
        plot_obj.savefig(os.path.join(output_path, ext, f"{base_filename_plot}.{ext}"))
        if logger:
            logger.info(f"Saved {os.path.join(output_path, extension, f'{base_filename_plot}.{ext}')}")

    def save_data(ext):
        filepath = os.path.join(output_path, ext, f"{base_filename_data}.{ext}")
        if isinstance(data_obj, dict):
            with open(filepath, "wb") as f:
                pickle.dump(data_obj, f, pickle.HIGHEST_PROTOCOL)
        else:
            raise NotImplementedError(f"{ext} type not supported (yet)")
        if logger:
            logger.info(f"Saved {os.path.join(output_path, extension, f'{base_filename_data}.{ext}')}")

    tasks = []
    if plot_filename_and_obj is not None:
        base_filename_plot, plot_obj = plot_filename_and_obj
        tasks.append((plot_extensions, save_canvas))
    if data_filename_and_obj is not None:
        base_filename_data, data_obj = data_filename_and_obj
        tasks.append((data_extensions, save_data))

    for extensions, _func in tasks:
        for extension in extensions:
            os.makedirs(os.path.join(output_path, extension), exist_ok=True)
            _func(extension)


def TH1_to_numpy(
    hist: ROOT.TH1,
    is_data: bool = False,
) -> Tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    Union[List[np.ndarray], None],
    Union[List[np.ndarray], None],
]:
    """
    Convert a ROOT TH1 histpgram to numpy arrays.

    Args:
        hist: The TH1 object to convert.
        is_data: Flag to indicate that it is a data histogram.

    Returns:
        Tuple containging:
            - edges: bin edge values as numpy histogram
            - x: x values as numpy array
            - y: y values as numpy array
            - x_errors: List of numpy arrays containing the x errors (low and high)
            - y_errors: List of numpy arrays containing the y errors (low and high)
    """
    n_bins = hist.GetNbinsX()
    edges = np.array(
        [
            hist.GetBinLowEdge(i + 1)
            for i in range(n_bins)
        ] + [hist.GetBinLowEdge(n_bins) + hist.GetBinWidth(n_bins)]
    )
    y_errors = [
        np.array([hist.GetBinErrorLow(i + 1) for i in range(n_bins)]),
        np.array([hist.GetBinErrorUp(i + 1) for i in range(n_bins)]),
    ]

    x_errors, x = None, None

    if is_data:  # if the histogram has center of mass x: extract during TGraph creation
        (_, x, y, x_err_down, x_err_up, _, _) = func.build_TGraph(
            hist,
            return_components=True,
            add_xerrors_in_graph=True,
        )
        x_errors = [np.array(x_err_down), np.array(x_err_up)]
    else:
        y = [hist.GetBinContent(i + 1) for i in range(n_bins)]

    return (edges, np.array(x), np.array(y), x_errors, y_errors)


def TGraphAsymmErrors_to_numpy(
    graph: ROOT.TGraphAsymmErrors,
) -> Tuple[
    np.ndarray,
    np.ndarray,
    List[np.ndarray],
    List[np.ndarray],
]:
    """
    Convert a ROOT TGraphAsymmErrors to numpy arrays.

    Args:
        graph: The TGraphAsymmErrors object to convert.

    Returns:
        Tuple containging:
            - edges: bin edge values as numpy histogram
            - x: x values as numpy array
            - y: y values as numpy array
            - x_errors: List of numpy arrays containing the x errors (low and high)
            - y_errors: List of numpy arrays containing the y errors (low and high)
    """
    n = graph.GetN()
    bin_edges = []

    x, y, xe_low, xe_high, ye_low, ye_high = (np.empty(n) for _ in range(6))
    for i in range(n):
        _x, _y = array.array("d", [0.]), array.array("d", [0.])
        graph.GetPoint(i, _x, _y)
        x[i], y[i] = _x[0], _y[0]
        xe_low[i], xe_high[i] = graph.GetErrorXlow(i), graph.GetErrorXhigh(i)
        ye_low[i], ye_high[i] = graph.GetErrorYlow(i), graph.GetErrorYhigh(i)
        if i == 0:
            bin_edges.append(x[i] - xe_low[i])
        bin_edges.append(x[i] + xe_high[i])

    return (
        np.array(bin_edges),
        np.array(x),
        np.array(y),
        [np.array(xe_low), np.array(xe_high)],
        [np.array(ye_low), np.array(ye_high)],
    )


def interval_splitted_category_label(category, var):
    """
    Function to split the category label into lower and upper bounds for LaTeX formatting
    if the category is defined as an interval.
    Args:
        category: The category dictionary containing the variable and its corresponding category.
        var: The variable for which the category is defined.
    Returns:
        str: The formatted string for the category label.
    """
    unit = r"\ GeV" if any(it in var for it in {'pt_'}) else ""
    lower, upper = category[var].split("#&&#")

    upper = f"{gd.category_dict[var].replace('$', '')} {upper}" + unit
    lower = f"{gd.category_dict[var].replace('$', '')} {lower}" + unit
    return f"$^{{{upper}}}_{{{lower}}}$"


def default_category_label(category, var):
    """
    Function to format the category label for LaTeX formatting

    Args:
        category: The category dictionary containing the variable and its corresponding category.
        var: The variable for which the category is defined.
    Returns:
        str: The formatted string for the category label.
    """
    unit = r"\ GeV" if any(it in var for it in {'pt_'}) else ""
    if var and "incl" not in var.lower():
        return f"{gd.category_dict[var]} ${category[var]}{unit}$"
    else:
        return f"{gd.category_dict[var]}"


def latex_adjust_selection_string(string):
    """
    Function to adjust the selection operators in the string to be used in a LaTeX formating.

    Args:
        string: The string to be adjusted.

    Returns:
        str: The adjusted string.
    """
    for a, b in [
        ("==", r"=\,"),
        (">=", r"\geq\,"),
        ("<=", r"\leq\,"),
        (">", r">\,"),
        ("<", r"<\,"),
        ("&&", r"\wedge\,"),
        ("||", r"\vee\,"),
        ("!=", r"\neq\,"),
    ]:
        string = string.replace(a, b if "$" in string else "$" + b + "$")
    return string


def plot_FFs(
    variable: str,
    ff_ratio: Any,
    uncertainties: Dict[str, Callable],
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

    ff_ratio = deepcopy(ff_ratio)
    uncertainties = deepcopy(uncertainties)

    recreation_summary = dict(
        variable=variable,
        era=era,
        channel=channel,
        process=process,
        category=category,
        draw_option=draw_option,
    )

    nominal_key, unc_up_key, unc_down_key = "nominal", "unc_up", "unc_down"
    mc_up_key, mc_down_key = "mc_subtraction_unc_up", "mc_subtraction_unc_down"

    if isinstance(ff_ratio, ROOT.TH1):
        _, x, y, x_err, y_err = TH1_to_numpy(ff_ratio, is_data=True)
    elif isinstance(ff_ratio, ROOT.TGraphAsymmErrors):
        _, x, y, x_err, y_err = TGraphAsymmErrors_to_numpy(ff_ratio)
    elif isinstance(ff_ratio, tuple):  # for recreating the plot
        x, y, x_err, y_err = ff_ratio

    _xmin, _xmax = x[0] - x_err[0][0], x[-1] + x_err[1][-1]
    x_sample = np.linspace(_xmin, _xmax, 1000)
    for k in uncertainties:  # for recreating the plot
        if not isinstance(uncertainties[k], np.ndarray):
            uncertainties[k] = np.vectorize(uncertainties[k])(x_sample)

    fig, ax = plt.subplots(1, 1, figsize=(11, 8))
    ax.set(
        xlim=(_xmin, _xmax),
        xlabel=gd.variable_dict[channel].get(variable, variable),
        ylim=(0, max(0.3, 1.7 * max(y))),
        ylabel=gd.FF_YAxis[process],
    )

    if category is not None:
        plot_text = f"{gd.channel_dict[channel]}, " + ", ".join(
            interval_splitted_category_label(category, var)
            if "#" in category[var]
            else default_category_label(category, var)
            for var in category.keys()
        )
    else:
        plot_text = gd.channel_dict[channel]
    plot_text = latex_adjust_selection_string(plot_text)

    ax.set_title(plot_text, loc="left", fontsize=25)
    ax.set_title(gd.era_dict[era], loc="right", fontsize=25)

    hep.cms.text(text=gd.default_CMS_text, ax=ax, loc=2, fontsize=20)

    data_handler = ax.errorbar(
        x,
        y,
        xerr=x_err,
        yerr=y_err,
        fmt="o",
        color="black",
        capsize=3,
        elinewidth=2,
    )

    handlers, labels = [data_handler], ["measured"]

    color = gd.color_dict["fit_graph_unc"]
    nominal, up, down = (
        uncertainties[nominal_key],
        uncertainties[unc_up_key],
        uncertainties[unc_down_key],
    )

    ax.plot(x_sample, nominal, lw=3, color=color)
    ax.fill_between(x_sample, up, down, alpha=0.25, lw=0, color=color)

    handlers.append((Patch(color=color, alpha=0.25), Line2D([-1], [-1], color=color, lw=3)))
    labels.append(gd.label_dict["fit_graph_unc"][draw_option])

    if mc_up_key in uncertainties and mc_down_key in uncertainties:
        up_mc, down_mc = uncertainties[mc_up_key], uncertainties[mc_down_key]
        mc_color = gd.color_dict["fit_graph_mc_sub"]
        ax.fill_between(x_sample, up_mc, down_mc, alpha=0.25, lw=0, color=mc_color)

        handlers.append((Patch(color=mc_color, alpha=0.25), Line2D([-1], [-1], color=mc_color, lw=3)))
        labels.append(gd.label_dict["fit_graph_mc_sub"][draw_option])

    # Add legend with both the combined entry and the histogram
    ax.legend(
        handles=handlers,
        labels=labels,
        loc="upper right",
        fontsize=20,
        labelspacing=0.3,
    )

    fig.tight_layout()

    selection = "_".join([f"{var}_{category[var]}" for var in category.keys()])
    base_filename = f"ff_{variable}_{process}_{selection}"
    save_plots_and_data(
        output_path=output_path,
        plot_filename_and_obj=(f"{base_filename}", fig),
        data_filename_and_obj=(
            f"data_{base_filename}",
            {
                "ff_ratio": (x, y, x_err, y_err),
                "uncertainties": uncertainties,
                **recreation_summary,
            },
        ) if save_data else None,
        logger=log,
    )
    plt.close()


def plot_data_mc_ratio(
    variable: str,
    hists: Dict[str, Any],  # TH1F-style ROOT histograms
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

    hists = {k: v.Clone() for k, v in hists.items()}  # clone method accounts for center-of-mass x values

    recreation_summary = dict(
        variable=variable,
        era=era,
        channel=channel,
        process=process,
        region=region,
        data=data,
        samples=samples,
        category=category,
        yscale=yscale,
    )

    if "edges" not in hists:  # for recreating the plot
        hists["edges"] = TH1_to_numpy(hists[samples[0]])[0]
    edges = hists["edges"]
    bin_centers, bin_widths = 0.5 * (edges[:-1] + edges[1:]), np.diff(edges)

    mc_stack, mc_labels, mc_colors = [], [], []
    total_mc = np.zeros_like(bin_centers)
    total_mc_err_up2 = np.zeros_like(bin_centers)
    total_mc_err_down2 = np.zeros_like(bin_centers)

    for sample in samples:
        if isinstance(hists[sample], ROOT.TH1):  # for recreating the plot
            _, _, values, _, y_errors = TH1_to_numpy(hists[sample])
            hists[sample] = (values, y_errors)
        values, y_errors = hists[sample]

        mc_stack.append(values)
        mc_labels.append(gd.label_dict[sample])
        mc_colors.append(gd.color_dict[sample])
        total_mc += values
        total_mc_err_down2 += y_errors[0]**2
        total_mc_err_up2 += y_errors[1]**2

    total_mc_err_up = np.sqrt(total_mc_err_up2)
    total_mc_err_down = np.sqrt(total_mc_err_down2)

    if isinstance(hists[data], ROOT.TH1):  # for recreating the plot
        _, data_x_values, data_y_values, data_x_errors, data_y_errors = TH1_to_numpy(hists[data], is_data=True)
        hists[data] = (data_x_values, data_y_values, data_x_errors, data_y_errors)
    data_x_values, data_y_values, data_x_errors, data_y_errors = hists[data]

    fig, (ax_top, ax_ratio) = plt.subplots(
        2, 1, figsize=(13, 12), gridspec_kw={"height_ratios": [0.75, 0.25]}, sharex=True
    )
    fig.subplots_adjust(hspace=0.05)

    ax_top.hist(
        [bin_centers] * len(mc_stack),
        bins=edges,
        weights=mc_stack,
        stacked=True,
        histtype="stepfilled",
        label=mc_labels,
        color=mc_colors,
        edgecolor="black",
    )

    ax_top.bar(
        bin_centers,
        height=total_mc_err_up + total_mc_err_down,
        width=bin_widths,
        bottom=total_mc - total_mc_err_down,
        color="gray",
        alpha=0.5,
        label="Stat. unc.",
        linewidth=0,
    )

    ax_top.errorbar(
        data_x_values,
        data_y_values,
        xerr=data_x_errors,
        yerr=data_y_errors,
        fmt="o",
        color="black",
        label=gd.label_dict[data],
        capsize=3,
        elinewidth=3,
    )

    ax_top.set_ylabel(r"$N_{Events}$")
    ax_top.set_yscale(yscale)
    if yscale == "linear":
        ax_top.set_ylim(0., 1.9 * max(max(data_y_values), max(total_mc)))
    else:
        ax_top.set_ylim(None, 1000.0 * max(max(data_y_values), max(total_mc)))
    ax_top.set_xlim(edges[0], edges[-1])
    ax_top.legend(ncol=2, loc="upper right", fontsize=20, labelspacing=0.1)

    formatter = SmartSciFormatter(ax_top.yaxis)
    ax_top.yaxis.set_major_formatter(formatter)

    ratio = np.divide(
        data_y_values,
        total_mc,
        out=np.ones_like(data_y_values),
        where=total_mc > 0,
    )
    ratio_err_down = np.divide(
        data_y_errors[0],
        total_mc,
        out=np.zeros_like(data_y_errors[0]),
        where=total_mc > 0,
    )
    ratio_err_up = np.divide(
        data_y_errors[1],
        total_mc,
        out=np.zeros_like(data_y_errors[1]),
        where=total_mc > 0,
    )

    ax_ratio.axhline(1.0, color="gray", linestyle="--")
    ax_ratio.bar(
        bin_centers,
        height=((total_mc_err_up + total_mc_err_down) / total_mc),
        width=bin_widths,
        bottom=(1 - total_mc_err_down / total_mc),
        color="gray",
        alpha=0.5,
        label="Stat. unc.",
        linewidth=0,
    )

    ax_ratio.errorbar(
        data_x_values,
        ratio,
        xerr=data_x_errors,
        yerr=(ratio_err_down, ratio_err_up),
        fmt="o",
        color="black",
        capsize=3,
        elinewidth=3,
    )

    ax_ratio.set_ylim(0.75, 1.25)
    ax_ratio.set_ylabel("Data / MC")
    ax_ratio.set_xlabel(gd.variable_dict[channel].get(variable, variable))
    ax_ratio.grid(axis="y")

    hep.cms.text(text=gd.default_CMS_text, ax=ax_top, loc=2, fontsize=20)

    if category is not None:
        plot_text = f"{gd.channel_dict[channel]}, {process.replace('_', ' ')}, " + ", ".join(
            interval_splitted_category_label(category, var)
            if "#" in category[var]
            else default_category_label(category, var)
            for var in category.keys()
        )
    else:
        plot_text = f"{gd.channel_dict[channel]}, {process}"
    plot_text = latex_adjust_selection_string(plot_text)

    ax_top.set_title(plot_text, loc="left", fontsize=25)
    ax_top.set_title(gd.era_dict[era], loc="right", fontsize=25)

    fig.tight_layout()

    hist_str = "reduced" if data == "data_subtracted" else "full"
    selection = '_'.join([f'{var}_{category[var]}' for var in category.keys()])
    base_filename = f"hist_ratio_{hist_str}_{variable}_{process}_{region}_{selection}"
    save_plots_and_data(
        output_path=output_path,
        plot_filename_and_obj=(f"{base_filename}_{yscale}", fig),
        data_filename_and_obj=(
            f"data_{base_filename}",
            {
                "hists": hists,
                **recreation_summary,
            },
        ) if save_data else None,
        logger=log,
    )
    plt.close(fig)


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

    hists = {k: v.Clone() for k, v in hists.items()}

    recreation_summary = dict(
        variable=variable,
        era=era,
        channel=channel,
        region=region,
        fraction_name=fraction_name,
        processes=processes,
        category=category,
    )

    if "edges" not in hists:  # for recreating the plot
        hists["edges"] = TH1_to_numpy(hists[processes[0]])[0]
    edges = hists["edges"]
    bin_centers = 0.5 * (edges[:-1] + edges[1:])

    # Stack the fractions
    mc_stack, mc_labels, mc_colors = [], [], []
    total_mc = np.zeros_like(bin_centers)

    for sample in processes:
        if isinstance(hists[sample], ROOT.TH1):  # for recreating the plot
            _, _, values, _, _ = TH1_to_numpy(hists[sample])
            hists[sample] = values
        values = hists[sample]
        mc_stack.append(values)
        mc_labels.append(gd.label_dict[sample])
        mc_colors.append(gd.color_dict[sample])
        total_mc += values

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.hist(
        [bin_centers] * len(processes),
        bins=edges,
        weights=mc_stack,
        stacked=True,
        color=mc_colors,
        label=mc_labels,
        histtype="stepfilled",
        edgecolor="black",
    )

    ax.set_ylabel("Fraction")
    ax.set_xlabel(gd.variable_dict[channel].get(variable, variable))
    ax.set_ylim(0, 1.3)
    ax.set_xlim(edges[0], edges[-1])

    ax.legend(loc="upper right", fontsize=20, labelspacing=0.2)

    hep.cms.text(text=gd.default_CMS_text, ax=ax, loc=2, fontsize=20)

    if category is not None:
        plot_text = f"{gd.channel_dict[channel]}, " + ", ".join(
            interval_splitted_category_label(category, var)
            if "#" in category[var]
            else default_category_label(category, var)
            for var in category.keys()
        )
    else:
        plot_text = f"{gd.channel_dict[channel]}"
    plot_text = latex_adjust_selection_string(plot_text)

    ax.set_title(plot_text, loc="left", fontsize=25)
    ax.set_title(gd.era_dict[era], loc="right", fontsize=25)

    fig.tight_layout()

    selection = '_'.join([f'{var}_{category[var]}' for var in category.keys()])
    base_filename = f"fraction_{fraction_name}_{variable}_{region}_{selection}"
    save_plots_and_data(
        output_path=output_path,
        plot_filename_and_obj=(f"{base_filename}", fig),
        data_filename_and_obj=(
            f"data_{base_filename}",
            {
                "hists": hists,
                **recreation_summary,
            },
        ) if save_data else None,
        logger=log,
    )
    plt.close(fig)


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
        corr_graph: Dict of the measured correction which is a smoothed function of "corr_hist" (includes variation)
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

    corr_hist = deepcopy(corr_hist)
    corr_graph = deepcopy(corr_graph)

    recreation_summary = dict(
        variable=variable,
        corr_graph=corr_graph,
        corr_name=corr_name,
        era=era,
        channel=channel,
        process=process,
        category=category,
    )

    lw = 3

    (
        color,
        color_stat_unct,
        color_bandwidth_unct,
        color_mc_sub_unct,
    ) = (
        gd.color_dict["correction_graph"],
        gd.color_dict["correction_graph_stat_unct"],
        gd.color_dict["correction_graph_bandwidth_unct"],
        gd.color_dict["correction_graph_mc_sub_unct"],
    )

    def pad(item):
        return np.pad(item, (0, 1), "edge")

    def step(_ax, _y, _c=color, _lw=lw, _ls="-"):
        return _ax.step(
            corr_graph["edges"],
            pad(corr_graph[_y]) if isinstance(_y, str) else _y,
            where="post",
            color=_c,
            lw=_lw,
            ls=_ls,
        )

    def bin_lookup(_values, edges, values):
        _idx = np.searchsorted(edges, _values, side="right") - 1
        _idx = np.clip(_idx, 0, len(values) - 1)
        return values[_idx]

    if isinstance(corr_hist, ROOT.TGraphAsymmErrors):
        _, x, y, x_err, y_err = TGraphAsymmErrors_to_numpy(corr_hist)
    elif isinstance(corr_hist, tuple):  # in case of recreating the plot
        x, y, x_err, y_err = corr_hist

    fig, ax = plt.subplots(2, 1, figsize=(14, 12), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    ax[0].set(
        xlim=(x[0] - x_err[0][0], x[-1] + x_err[1][-1]),
        ylim=(0, max(2.2, 1.5 * max(y))),
        ylabel="Correction",
    )
    ax[1].set(
        xlabel=gd.variable_dict[channel].get(variable, variable),
        ylabel="Ratio",
        ylim=(0.5, 1.5),
    )

    ax[0].hlines(1, *ax[0].get_xlim(), colors="black", lw=0.8)
    ax[1].hlines(1, *ax[1].get_xlim(), colors=color, lw=lw)

    if category is not None:
        plot_text = f"{gd.channel_dict[channel]}, {process}, " + ", ".join(
            interval_splitted_category_label(category, var)
            if "#" in category[var]
            else default_category_label(category, var)
            for var in category.keys()
        )
    else:
        plot_text = f"{gd.channel_dict[channel]}, {process}"
    plot_text = latex_adjust_selection_string(plot_text)

    ax[0].set_title(plot_text, loc="left", fontsize=20)
    ax[0].set_title(gd.era_dict[era], loc="right", fontsize=20)

    if "non_closure" in corr_name:
        ax[0].text(0.08, 0.8, "non closure correction", transform=ax[0].transAxes, fontsize=20)
    elif "DR_SR" in corr_name:
        ax[0].text(0.08, 0.8, "DR to SR correction", transform=ax[0].transAxes, fontsize=20)

    hep.cms.text(text=gd.default_CMS_text, ax=ax[0], loc=1, fontsize=20)

    total_unct_up = np.sqrt(
        (
            corr_graph["StatUp"] - corr_graph["nominal"]
        )**2 + (
            corr_graph["BandAsymUp"] - corr_graph["nominal"]
        )**2 + (
            corr_graph["MCShiftUp"] - corr_graph["nominal"]
        )**2
    )
    total_unct_down = np.sqrt(
        (
            corr_graph["StatDown"] - corr_graph["nominal"]
        ) ** 2 + (
            corr_graph["BandAsymDown"] - corr_graph["nominal"]
        ) ** 2 + (
            corr_graph["MCShiftDown"] - corr_graph["nominal"]
        ) ** 2
    )

    ax[0].fill_between(
        corr_graph["edges"],
        pad(corr_graph["nominal"] - total_unct_down),
        pad(corr_graph["nominal"] + total_unct_up),
        alpha=0.15,
        step="post",
        color=color,
        lw=0,
    )
    ax[1].fill_between(
        corr_graph["edges"],
        pad((corr_graph["nominal"] - total_unct_down) / corr_graph["nominal"]),
        pad((corr_graph["nominal"] + total_unct_up) / corr_graph["nominal"]),
        alpha=0.15,
        step="post",
        color=color,
        lw=0,
    )

    _step = partial(step, _lw=lw - 2, _c=color_stat_unct)
    stat_unct_up, *_ = _step(ax[0], "StatUp")
    stat_unct_down, *_ = _step(ax[0], "StatDown")
    _step(ax[1], pad(corr_graph["StatUp"] / corr_graph["nominal"]))
    _step(ax[1], pad(corr_graph["StatDown"] / corr_graph["nominal"]))

    _step = partial(step, _lw=lw - 2, _c=color_bandwidth_unct)
    bandwidth_unct_up, *_ = _step(ax[0], "BandAsymUp", _ls="dashed")
    bandwidth_unct_down, *_ = _step(ax[0], "BandAsymDown", _ls="dashdot")
    _step(ax[1], pad(corr_graph["BandAsymUp"] / corr_graph["nominal"]), _ls="dashed")
    _step(ax[1], pad(corr_graph["BandAsymDown"] / corr_graph["nominal"]), _ls="dashdot")

    _step = partial(step, _lw=lw - 2, _c=color_mc_sub_unct)
    mc_sub_unct_up, *_ = _step(ax[0], "MCShiftUp", _ls="dashed")
    mc_sub_unct_down, *_ = _step(ax[0], "MCShiftDown", _ls="dashdot")
    _step(ax[1], pad(corr_graph["MCShiftUp"] / corr_graph["nominal"]), _ls="dashed")
    _step(ax[1], pad(corr_graph["MCShiftDown"] / corr_graph["nominal"]), _ls="dashdot")

    step(ax[0], "nominal")

    data_handle = ax[0].errorbar(x, y, xerr=x_err, yerr=y_err, fmt="o", color="black", capsize=3, elinewidth=2)
    ax[1].errorbar(
        x,
        y / bin_lookup(x, corr_graph["edges"], corr_graph["nominal"]),
        xerr=x_err,
        yerr=y_err / np.interp(x, corr_graph["edges"][:-1] + np.diff(corr_graph["edges"]) / 2, corr_graph["nominal"]),
        fmt="o",
        color="black",
        capsize=3,
        elinewidth=2,
    )

    ax[0].legend(
        [
            data_handle,
            MergedTuple((Patch(color=color, alpha=0.15), Line2D([], [], color=color, linewidth=lw))),
            VTuple((stat_unct_up, stat_unct_down)),
            VTuple((bandwidth_unct_up, bandwidth_unct_down)),
            VTuple((mc_sub_unct_up, mc_sub_unct_down)),
        ],
        [
            "measured",
            "Smoothed Curve + Total unct.",
            "Stat. unct.",
            "Bandwidth unct.",
            "MC subtraction unct.",
        ],
        loc="upper right",
        fontsize=20,
        labelspacing=0.3,
        handler_map=HANDLER_MAP,
    )

    fig.tight_layout()

    selection = ""
    if category is not None:
        selection = '_'.join([""] + [f'{var}_{category[var]}' for var in category.keys()])

    save_plots_and_data(
        output_path=output_path,
        plot_filename_and_obj=(f"corr_{process}_{corr_name}{selection}", fig),
        data_filename_and_obj=(
            f"data_corr_{process}_{corr_name}{selection}",
            {
                "corr_hist": (x, y, x_err, y_err),
                **recreation_summary,
            },
        ) if save_data else None,
        logger=log,
    )
    plt.close()
