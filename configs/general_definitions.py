from collections import defaultdict
from copy import deepcopy
from typing import Any, List, Union

import ROOT

from helper.hooks_and_patches import _EXTRA_PARAM_FLAG, _EXTRA_PARAM_MEANS

#### For fake factor fits ####

random_seed = 19
default_fit_option = "poly_1"
default_correction_option = "smoothed"

default_CMS_text = "Own work (Data/Simulation)"


class AutoGetDict(dict):
    def __getitem__(self, key: Any) -> Any:
        if key not in self:
            print(f"WARNING: Key '{key}' not found in dictionary, returning the key itself.")
        return self.get(key, key)


def get_default_fit_function_limit_kwargs(
    binning: List[float],
    hist: Union[ROOT.TH1, None] = None,
) -> dict:
    """
    Function to determine the default fit function limit kwargs for the fit function.
    Default is set to be the binning range. If the histogram has the extra parameter flag
    the limits are set to first bin (left edge) and the center of mass of the last bin
    plus the minimum of the distance to the next bin edge and the distance to the previous
    bin edge.

    Args:
        binning: List of floats defining the binning range
        hist: Histogram object to be used for the fit function

    Returns:
        dict: Dictionary with the fit function limit kwargs
    """
    params = {
        "limit_x": {
            "nominal": (binning[0], binning[-1]),
            "up": (-float("inf"), float("inf")),
            "down": (-float("inf"), float("inf")),
        },
    }
    if hist is not None and hasattr(hist, _EXTRA_PARAM_FLAG) and getattr(hist, _EXTRA_PARAM_FLAG):
        x, n = getattr(hist, _EXTRA_PARAM_MEANS)[-1], hist.GetNbinsX()
        a, b = hist.GetBinLowEdge(n), hist.GetBinLowEdge(n + 1)
        params["limit_x"]["nominal"] = (binning[0], x + min(b - x, x - a))
    return params


def get_default_bandwidth(
    binning: List[float],
) -> float:
    """
    Function to determine the smoothness factor for the fit function, default is set to be
    1/5 of the binning range.

    Args:
        binning: List of floats defining the binning range
    Returns:
        float: Smoothness factor for the fit function
    """

    factor = (binning[-1] - binning[0]) / 5.0
    return factor


#### For plotting ####
OFFICIAL_CMS_COLOR_PALLET = {
    6: ["#5790fc", "#f89c20", "#e42536", "#964a8b", "#9c9ca1", "#7a21dd"],
    10: ["#3f90da", "#ffa90e", "#bd1f01", "#94a4a2", "#832db6", "#a96b59", "#e76300", "#b9ac70", "#717581", "#92dadd"],
}


# label definitions for y-axis
FF_YAxis = AutoGetDict(
    {
        "ttbar": r"$F_{F}^{t\bar{t}}$",
        "ttbar_subleading": r"$F_{F}^{t\bar{t}}$",
        "Wjets": r"$F_{F}^{Wjets}$",
        "QCD": r"$F_{F}^{QCD}$",
        "QCD_subleading": r"$F_{F}^{QCD}$",
    }
)

# definitions for channels
channel_dict = AutoGetDict({"et": r"$e\tau_{h}$", "mt": r"$\mu\tau_{h}$", "tt": r"$\tau_{h}\tau_{h}$"})

# definitions for era and luminosity
era_dict = AutoGetDict(
    {
        "2016preVFP": r"$19.5\,fb^{-1}$ (2016preVFP, 13 TeV)",
        "2016postVFP": r"$16.8\,fb^{-1}$ (2016postVFP, 13 TeV)",
        "2017": r"$41.5\,fb^{-1}$ (2017, 13 TeV)",
        "2018": r"$59.8\,fb^{-1}$ (2018, 13 TeV)",
    }
)

# definitions for process color is the histograms + colors for the fitted graphs
color_dict = {
    "QCD": OFFICIAL_CMS_COLOR_PALLET[10][7],  # (185, 172, 112)
    "diboson_J": "#94A484",  # (148, 164, 132)
    "diboson_L": OFFICIAL_CMS_COLOR_PALLET[10][3],  # (148, 164, 162)
    "diboson_T": "#94A4C0",  # (148, 164, 192)
    "Wjets": OFFICIAL_CMS_COLOR_PALLET[10][6],  # (231, 99, 0)
    "ttbar_J": "#652DB6",  # (101, 45, 182)
    "ttbar_L": OFFICIAL_CMS_COLOR_PALLET[10][4],  # (131, 45, 182)
    "ttbar_T": "#A12DB6",  # (161, 45, 182)
    "DYjets_J": "#3F72DA",  # (63, 114, 218)
    "DYjets_L": OFFICIAL_CMS_COLOR_PALLET[10][0],  # (63, 144, 218)
    "DYjets_T": "#3FAEDA",  # (63, 174, 218)
    "ST_J": "#717563",  # (113, 117, 99)
    "ST_L": OFFICIAL_CMS_COLOR_PALLET[10][8],  # (113, 117, 129)
    "ST_T": "#71759F",  # (113, 117, 159)
    "embedding": OFFICIAL_CMS_COLOR_PALLET[10][1],  # (255, 169, 14)
    "tau_fakes": "#B9AC70",  # (185, 172, 112)
    "data_ff": OFFICIAL_CMS_COLOR_PALLET[10][1],  # (255, 169, 14)
    "fit_graph_mc_sub": OFFICIAL_CMS_COLOR_PALLET[6][0],
    "fit_graph_unc": OFFICIAL_CMS_COLOR_PALLET[6][2],
    "correction_graph": OFFICIAL_CMS_COLOR_PALLET[10][0],
    "correction_graph_stat_unct": OFFICIAL_CMS_COLOR_PALLET[10][0],
    "correction_graph_bandwidth_unct": OFFICIAL_CMS_COLOR_PALLET[10][1],
    "correction_graph_mc_sub_unct": OFFICIAL_CMS_COLOR_PALLET[10][5],
}

# definitions for process labels on the plots
label_dict = AutoGetDict(
    {
        "QCD": "QCD multijet",
        "diboson_J": r"diboson (jet$\rightarrow\tau_{h}$)",
        "diboson_L": r"diboson ($\ell\rightarrow\tau_{h}$)",
        "diboson_T": r"diboson ($\tau\rightarrow\tau_{h}$)",
        "Wjets": r"W$\rightarrow\ell\nu$",
        "ttbar_J": r"$t\bar{t}$ (jet$\rightarrow\tau_{h}$)",
        "ttbar_L": r"$t\bar{t}$ ($\ell\rightarrow\tau_{h}$)",
        "ttbar_T": r"$t\bar{t}$ ($\tau\rightarrow\tau_{h}$)",
        "DYjets_J": r"$Z\rightarrow\ell\ell$ (jet$\rightarrow\tau_{h}$)",
        "DYjets_L": r"$Z\rightarrow\ell\ell$ ($\ell\rightarrow\tau_{h}$)",
        "DYjets_T": r"$Z\rightarrow\ell\ell$ ($\tau\rightarrow\tau_{h}$)",
        "ST_J": r"t (jet$\rightarrow\tau_{h}$)",
        "ST_L": r"t ($\ell\rightarrow\tau_{h}$)",
        "ST_T": r"t ($\tau\rightarrow\tau_{h}$)",
        "embedding": r"$\tau$ embedded",
        "data": "data",
        "data_subtracted": "reduced data",
        "data_ff": r"data with $F_F$'s",
        "tau_fakes": r"jet$\rightarrow\tau_{h}$",
        "fit_graph_unc": {
            **{f"poly_{i}": f"poly({i}) best fit unc." for i in range(7)},
            "binwise": "plain histogram unc.",
        },
        "fit_graph_mc_sub": {
            **{f"poly_{i}": f"poly({i}) best fit (MC subtr. unc.)" for i in range(7)},
            "binwise": "plain histogram (MC subtr. unc.)",
        },
    }
)

# definitions to translate variable to readable language, channel dependent
channel_indipendent_variable_dict = AutoGetDict(
    {
        "njets": r"$N_{jets}$",
        "metphi": r"$\phi(p_T^{miss})$",
        "met": r"$p_T^{miss}$ (GeV)",
        "m_vis": r"$m_{vis}$ (GeV)",
        "nbtag": r"$N_{b-jets}$",
    }
)
variable_dict = {
    "mt": AutoGetDict(
        {
            "pt_1": r"$p_{T}^{\mu}$ (GeV)",
            "eta_1": r"$\eta^{\mu}$",
            "phi_1": r"$\phi^{\mu}$",
            "iso_1": r"$iso^{\mu}$",
            "mt_1": r"$m_{T}(\mu,p_T^{miss})$ (GeV)",
            "pt_2": r"$p_{T}^{\tau_{h}}$ (GeV)",
            "eta_2": r"$\eta^{\tau_{h}}$",
            "phi_2": r"$\phi^{\tau_{h}}$",
            "mass_2": r"$\tau_{h}$ mass",
            "deltaR_ditaupair": r"$\Delta R(\mu,\tau_{h})$",
            "tau_decaymode_2": r"$\tau_{h}^{DM}$",
            **channel_indipendent_variable_dict,
        }
    ),
    "et": AutoGetDict(
        {
            "pt_1": r"$p_{T}^{e}$ (GeV)",
            "eta_1": r"$\eta^{e}$",
            "phi_1": r"$\phi^{e}$",
            "iso_1": r"$iso^{e}$",
            "mt_1": r"$m_{T}(e,p_T^{miss})$ (GeV)",
            "pt_2": r"$p_{T}^{\tau_{h}}$ (GeV)",
            "eta_2": r"$\eta^{\tau_{h}}$",
            "phi_2": r"$\phi^{\tau_{h}}$",
            "mass_2": r"$\tau_{h}$ mass",
            "deltaR_ditaupair": r"$\Delta R(e,\tau_{h})$",
            "tau_decaymode_2": r"$\tau_{h}^{DM}$",
            **channel_indipendent_variable_dict,
        },
    ),
    "tt": AutoGetDict(
        {
            "pt_1": r"$p_{T}^{\tau_{h,1}}$ (GeV)",
            "eta_1": r"$\eta^{\tau_{h,1}}$",
            "phi_1": r"$\phi^{\tau_{h,1}}$",
            "pt_2": r"$p_{T}^{\tau_{h,2}}$ (GeV)",
            "eta_2": r"$\eta^{\tau_{h,2}}$",
            "phi_2": r"$\phi^{\tau_{h,2}}$",
            "mass_1": r"$\tau_{h,1}$ mass",
            "mass_2": r"$\tau_{h,2}$ mass",
            "deltaR_ditaupair": r"$\Delta R(\tau_{h, 1},\tau_{h, 2})$",
            "tau_decaymode_1": r"$\tau_{h,1}^{DM}$",
            "tau_decaymode_2": r"$\tau_{h,2}^{DM}$",
            **channel_indipendent_variable_dict,
        }
    ),
}

# definitions to translate category cuts to readable language
category_dict = AutoGetDict(
    {
        "incl": r"incl.",
        "njets": r"$N_{jets}$",
        "nbtag": r"$N_{b-jets}$",
        "deltaR_ditaupair": r"$\Delta R(\ell,\tau_{h})$",
        "tau_decaymode_1": r"$\tau_{h,1}^{DM}$",
        "tau_decaymode_2": r"$\tau_{h,2}^{DM}$",
        "pt_1": r"$p_{T}^{\mu}$",  # for mt channel
        "pt_2": r"$p_{T}^{\tau_{h}}$",
    }
)

### For correctionlib ###

# intern naming translation helper for variables
variable_translator = {
    "QCD": "QCD",
    "Wjets": "Wjets",
    "ttbar_J": "ttbar",
}
# definitions for the variable type, needed for correctionlib
variable_type = defaultdict(
    lambda: "real",
    {
        "pt_2": "real",
        "pt_1": "real",
        "mt_1": "real",
        "mass_1": "real",
        "mass_2": "real",
        "iso_1": "real",
        "m_vis": "real",
        "bpt_1": "real",
        "njets": "real",
        "nbtag": "real",
        "deltaR_ditaupair": "real",
        "tau_decaymode_1": "real",
        "tau_decaymode_2": "real",
    }
)
# definitions for variable descriptions, needed for correctionlib, #var_max and #var_min are replaced later by using the variable binning
variable_description = AutoGetDict(
    {
        "pt_2": "transverse momentum of the subleading hadronic tau in the tau pair; measured between #var_min and #var_max GeV; for higher/lower pt's the edge values are used",
        "pt_1": "transverse momentum of the leading leptonic/hadronic tau in the tau pair; measured between #var_min and #var_max GeV; for higher/lower pt's the edge values are used",
        "iso_1": "isolation of the lepton in the tau pair; measured between #var_min and #var_max GeV; for higher/lower isolation values the edge values are used",
        "mt_1": "transverse mass of the lepton and MET in the tau pair; measured between #var_min and #var_max GeV; for higher/lower mt's the edge values are used",
        "mass_1": "mass of the leading hadronic tau in the tau pair; measured between #var_min and #var_max GeV; for higher/lower masses the edge values are used",
        "mass_2": "mass of the subleading hadronic tau in the tau pair; measured between #var_min and #var_max GeV; for higher/lower masses the edge values are used",
        "m_vis": "invariant mass of the visible di-tau decay products; measured between #var_min and #var_max GeV; for higher/lower m_vis's the edge values are used",
        "bpair_pt_1": "transverse momentum of the hardest b-tagged jet; measured between #var_min and #var_max GeV; for higher/lower pt's the edge values are used",
        "njets": "number of jets in an event; the defined categories are ",
        "nbtag": "number of b-tagged jets in an event; the defined categories are ",
        "deltaR_ditaupair": "spatial distance between the tau pair with deltaR ",
        "tau_decaymode_1": "decay mode of the leading tau in the tau pair; the defined categories are ",
        "tau_decaymode_2": "decay mode of the subleading tau in the tau pair; the defined categories are ",
    }
)
# intern naming translation helper for fit uncertainties
# naming has to match the one used in helper/ff_functions.py -> fit_function()
ff_variation_dict = {
    **{
        k: {
            "unc_up": "FFuncUp",
            "unc_down": "FFuncDown",
            "mc_subtraction_unc_up": "FFmcSubUncUp",
            "mc_subtraction_unc_down": "FFmcSubUncDown",
        }
        for k in ["QCD", "QCD_subleading", "Wjets"]
    },
    **{
        k: {
            "unc_up": "FFuncUp",
            "unc_down": "FFuncDown",
        }
        for k in ["ttbar", "ttbar_subleading"]
    },
}
# intern naming translation helper for fraction uncertainties
# naming has to match the one used in helper/ff_functions.py -> add_fraction_variations()
frac_variation_dict = {
    "QCD": {
        "frac_QCD_up": "fracQCDUncUp",
        "frac_QCD_down": "fracQCDUncDown",
    },
    "Wjets": {
        "frac_Wjets_up": "fracWjetsUncUp",
        "frac_Wjets_down": "fracWjetsUncDown",
    },
    "ttbar_J": {
        "frac_ttbar_J_up": "fracTTbarUncUp",
        "frac_ttbar_J_down": "fracTTbarUncDown",
    },
}
