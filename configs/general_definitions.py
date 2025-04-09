import ROOT
from typing import Union, List
from helper.hooks_and_patches import _EXTRA_PARAM_FLAG, _EXTRA_PARAM_MEANS
from copy import deepcopy

#### For fake factor fits ####

random_seed = 19
default_fit_option = "poly_1"
default_correction_option = "smoothed"

default_CMS_text = "Own work (Data/Simulation)"


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
FF_YAxis_root = {
    "ttbar": "FF_{t#bar{t}}",
    "ttbar_subleading": "FF_{t#bar{t}}",
    "Wjets": "FF_{Wjets}",
    "QCD": "FF_{QCD}",
    "QCD_subleading": "FF_{QCD}",
}
FF_YAxis_root_mpl = {
    "ttbar": r"$F_{F}^{t\bar{t}}$",
    "ttbar_subleading": r"$F_{F}^{t\bar{t}}$",
    "Wjets": r"$F_{F}^{Wjets}$",
    "QCD": r"$F_{F}^{QCD}$",
    "QCD_subleading": r"$F_{F}^{QCD}$",
}
FF_YAxis = FF_YAxis_root  # current default: ROOT plotting

# definitions for channels
channel_dict_root = {"et": "e#tau_{h}", "mt": "#mu#tau_{h}", "tt": "#tau_{h}#tau_{h}"}
channel_dict_mpl = {"et": r"$e\tau_{h}$", "mt": r"$\mu\tau_{h}$", "tt": r"$\tau_{h}\tau_{h}$"}
channel_dict = channel_dict_root  # current default: ROOT plotting

# definitions for era and luminosity TODO: 2016
era_dict_root = {
    "2016preVFP": "19.5 fb^{-1} (2016preVFP, 13 TeV)",
    "2016postVFP": "16.8 fb^{-1} (2016postVFP, 13 TeV)",
    "2017": "41.5 fb^{-1} (2017, 13 TeV)",
    "2018": "59.8 fb^{-1} (2018, 13 TeV)",
}
era_dict_mpl = {
    "2016preVFP": r"$19.5\,fb^{-1}$ (2016preVFP, 13 TeV)",
    "2016postVFP": r"$16.8\,fb^{-1}$ (2016postVFP, 13 TeV)",
    "2017": r"$41.5\,fb^{-1}$ (2017, 13 TeV)",
    "2018": r"$59.8\,fb^{-1}$ (2018, 13 TeV)",
}
era_dict = era_dict_root  # current default: ROOT plotting

# definitions for process color is the histograms + colors for the fitted graphs
color_dict = {
    "QCD": (185, 172, 112),
    "diboson_J": (148, 164, 132),
    "diboson_L": (148, 164, 162),
    "diboson_T": (148, 164, 192),
    "Wjets": (231, 99, 0),
    "ttbar_J": (101, 45, 182),
    "ttbar_L": (131, 45, 182),
    "ttbar_T": (161, 45, 182),
    "DYjets_J": (63, 114, 218),
    "DYjets_L": (63, 144, 218),
    "DYjets_T": (63, 174, 218),
    "ST_J": (113, 117, 99),
    "ST_L": (113, 117, 129),
    "ST_T": (113, 117, 159),
    "embedding": (255, 169, 14),
    "tau_fakes": (185, 172, 112),
    "data_ff": (255, 169, 14),
    "fit_graph_mc_sub": ROOT.kGreen,
    "fit_graph_unc": ROOT.kRed,
}
# definitions for process labels on the plots
label_dict_root = {
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
    "fit_graph_unc": {
        **{f"poly_{i}": f"poly({i}) best fit uncertainty" for i in range(7)},
        "binwise": "plain histogram unc.",
    },
    "fit_graph_mc_sub": {
        **{f"poly_{i}": f"poly({i}) best fit (MC subtraction unc.)" for i in range(7)},
        "binwise": "plain histogram (MC subtraction unc.)",
    },
}
label_dict_mpl = deepcopy(label_dict_root).update(
    {
        "diboson_J": r"Diboson (jet$\rightarrow\tau_{h}$)",
        "diboson_L": r"Diboson ($\ell\rightarrow\tau_{h}$)",
        "diboson_T": r"Diboson (genuine $\tau_{h}$)",
        "ttbar_J": r"$t\bar{t}$ (jet$\rightarrow\tau_{h}$)",
        "ttbar_L": r"$t\bar{t}$ ($\ell\rightarrow\tau_{h}$)",
        "ttbar_T": r"$t\bar{t}$ (genuine $\tau_{h}$)",
        "DYjets_J": r"$Z\rightarrow\ell\ell$ (jet$\rightarrow\tau_{h}$)",
        "DYjets_L": r"$Z\rightarrow\ell\ell$ ($\ell\rightarrow\tau_{h}$)",
        "DYjets_T": r"$Z\rightarrow\ell\ell$ (genuine $\tau_{h}$)",
        "ST_J": r"ST (jet$\rightarrow\tau_{h}$)",
        "ST_L": r"ST ($\ell\rightarrow\tau_{h}$)",
        "ST_T": r"ST (genuine $\tau_{h}$)",
        "embedding": r"$\tau$ embedded",
        "data_ff": "Data with $F_F$s",
        "tau_fakes": r"jet$\rightarrow\tau_{h}$",
    },
)
label_dict = label_dict_root  # current default: ROOT plotting

# definitions to translate variable to readable language, channel dependent
variable_dict_root = {
    "et": {
        "pt_1": "p_{T}(e) (GeV)",
        "eta_1": "#eta(e)",
        "phi_1": "#phi(e)",
        "iso_1": "isolation (e)",
        "mt_1": "m_{T}(e, #slash{E}_{T}) (GeV)",
        "pt_2": "p_{T}(#tau_{h}) (GeV)",
        "eta_2": "#eta(#tau_{h})",
        "phi_2": "#phi(#tau_{h})",
        "mass_2": "#tau_{h} mass",
        "m_vis": "m_{vis} (GeV)",
        "met": "MET (GeV)",
        "metphi": "MET #phi",
        "njets": "number of jets",
        "nbtag": "number of b-tagged jets",
        "deltaR_ditaupair": "#DeltaR(e#tau_{h})",
        "tau_decaymode_2": "#tau_{h}^{DM} ",
    },
    "mt": {
        "pt_1": "p_{T}(#mu) (GeV)",
        "eta_1": "#eta(#mu)",
        "phi_1": "#phi(#mu)",
        "iso_1": "isolation (#mu)",
        "mt_1": "m_{T}(#mu, #slash{E}_{T}) (GeV)",
        "pt_2": "p_{T}(#tau_{h}) (GeV)",
        "eta_2": "#eta(#tau_{h})",
        "phi_2": "#phi(#tau_{h})",
        "mass_2": "#tau_{h} mass",
        "m_vis": "m_{vis} (GeV)",
        "met": "MET (GeV)",
        "metphi": "MET #phi",
        "njets": "number of jets",
        "nbtag": "number of b-tagged jets",
        "deltaR_ditaupair": "#DeltaR(#mu#tau_{h})",
        "tau_decaymode_2": "#tau_{h}^{DM}",
    },
    "tt": {
        "pt_1": "leading p_{T}(#tau_{h}) (GeV)",
        "eta_1": "leading #eta(#tau_{h})",
        "phi_1": "leading #phi(#tau_{h})",
        "pt_2": "subleading p_{T}(#tau_{h}) (GeV)",
        "eta_2": "subleading #eta(#tau_{h})",
        "phi_2": "subleading #phi(#tau_{h})",
        "mass_1": "leading #tau_{h} mass",
        "mass_2": "subleading #tau_{h} mass",
        "m_vis": "m_{vis} (GeV)",
        "met": "MET (GeV)",
        "metphi": "MET #phi",
        "njets": "number of jets",
        "nbtag": "number of b-tagged jets",
        "deltaR_ditaupair": "#DeltaR(#tau_{h}#tau_{h})",
        "tau_decaymode_1": "#tau_{h, 1}^{DM}",
        "tau_decaymode_2": "#tau_{h, 2}^{DM}",
    },
}
variable_dict_mpl = {
    "mt": {
        "pt_1": r"$p_{T}(\mu)$ (GeV)",
        "eta_1": r"$\eta(\mu)$",
        "phi_1": r"$\phi(\mu)$",
        "iso_1": r"$iso(\mu)$",
        "mt_1": r"$m_{T}(\mu, \slash{E}_{T})$ (GeV)",
        "pt_2": r"$p_{T}(\tau_{h})$ (GeV)",
        "eta_2": r"$\eta(\tau_{h})$",
        "phi_2": r"$\phi(\tau_{h})$",
        "mass_2": r"$\tau_{h}$ mass",
        "m_vis": r"$m_{vis}$ (GeV)",
        "met": r"$MET$ (GeV)",
        "metphi": r"$MET\ \phi$",
        "njets": r"$N_{jets}$",
        "nbtag": r"$N_{b-jets}$",
        "deltaR_ditaupair": r"$\Delta R(\mu \tau_{h})$",
        "tau_decaymode_2": r"$\tau_{h}^{DM}$",
    },
    "et": {
        "pt_1": r"$p_{T}(e)$ (GeV)",
        "eta_1": r"$\eta(e)$",
        "phi_1": r"$\phi(e)$",
        "iso_1": r"$iso(e)$",
        "mt_1": r"$m_{T}(e, \slash{E}_{T})$ (GeV)",
        "pt_2": r"$p_{T}(\tau_{h})$ (GeV)",
        "eta_2": r"$\eta(\tau_{h})$",
        "phi_2": r"$\phi(\tau_{h})$",
        "mass_2": r"$\tau_{h}$ mass",
        "m_vis": r"$m_{vis}$ (GeV)",
        "met": r"$MET$ (GeV)",
        "metphi": r"$MET \phi$",
        "njets": r"$N_{jets}$",
        "nbtag": r"$N_{b-jets}$",
        "deltaR_ditaupair": r"$\Delta R(e \tau_{h})$",
        "tau_decaymode_2": r"$\tau_{h}^{DM}$",
    },
    "tt": {
        "pt_1": r"$p_{T}(\tau_{h, 1})$ (GeV)",
        "eta_1": r"$\eta(\tau_{h, 1})$",
        "phi_1": r"$\phi(\tau_{h, 1})$",
        "pt_2": r"$p_{T}(\tau_{h, 2})$ (GeV)",
        "eta_2": r"$\eta(\tau_{h, 2})$",
        "phi_2": r"$\phi(\tau_{h, 2})$",
        "mass_1": r"$\tau_{h, 1}$ mass",
        "mass_2": r"$\tau_{h, 2}$ mass",
        "m_vis": r"$m_{vis}$ (GeV)",
        "met": r"$MET$ (GeV)",
        "metphi": r"$MET \phi$",
        "njets": r"$N_{jets}$",
        "nbtag": r"$N_{b-jets}$",
        "deltaR_ditaupair": r"$\Delta R(\tau_{h} \tau_{h})$",
        "tau_decaymode_1": r"$\tau_{h, 1}^{DM}$",
        "tau_decaymode_2": r"$\tau_{h, 2}^{DM}$",
    },
}
variable_dict = variable_dict_root  # current default: ROOT plotting

# definitions to translate category cuts to readable language
category_dict_root = {
    "incl": "incl.",
    "njets": "N_{jets}",
    "nbtag": "N_{b-jets}",
    "deltaR_ditaupair": "#Delta" + "R(l#tau_{h})",
    "tau_decaymode_1": "#tau_{h, 1}^{DM}",
    "tau_decaymode_2": "#tau_{h, 2}^{DM}",
    "pt_1": "p_{T}(#mu)",
}
category_dict_mpl = {
    "incl": r"incl.",
    "njets": r"$N_{jets}$",
    "nbtag": r"$N_{b-jets}$",
    "deltaR_ditaupair": r"$\Delta R(\ell\tau_{h})$",
    "tau_decaymode_1": r"$\tau_{h, 1}^{DM}$",
    "tau_decaymode_2": r"$\tau_{h, 2}^{DM}$",
    "pt_1": r"$p_{T}(e)$",
    "pt_2": r"$p_{T}(\tau_{h})$",
}
category_dict = category_dict_root  # current default: ROOT plotting

### For correctionlib ###

# intern naming translation helper for variables
variable_translator = {
    "QCD": "QCD",
    "Wjets": "Wjets",
    "ttbar_J": "ttbar",
}
# definitions for the variable type, needed for correctionlib
variable_type = {
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
# definitions for variable descriptions, needed for correctionlib, #var_max and #var_min are replaced later by using the variable binning
variable_description = {
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
