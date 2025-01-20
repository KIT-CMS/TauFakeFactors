import ROOT

### for preselection ###

# list of output features, they can change depending on the analysis or channel (tau vs jet ID is added internally)
output_features = {
    "nmssm": {
        "et": ["fj_Xbb_hadflavor", "fj_Xbb_nBhad", "fj_Xbb_nChad", "weight", "btag_weight", "pNet_Xbb_weight", "njets", "nbtag", "q_1", "pt_2", "q_2", "gen_match_2", "m_vis", "mt_1", "deltaR_ditaupair", "pt_1", "iso_1", "metphi", "extramuon_veto", "extraelec_veto", "dilepton_veto", "fj_Xbb_pt", "fj_Xbb_eta", "fj_Xbb_particleNet_XbbvsQCD", "bpair_pt_1", "bpair_pt_2", "bpair_btag_value_2", "bpair_eta_2", "met", "mass_2", "tau_decaymode_2"],
        "mt": ["fj_Xbb_hadflavor", "fj_Xbb_nBhad", "fj_Xbb_nChad", "weight", "btag_weight", "pNet_Xbb_weight", "njets", "nbtag", "q_1", "pt_2", "q_2", "gen_match_2", "m_vis", "mt_1", "deltaR_ditaupair", "pt_1", "iso_1", "metphi", "extramuon_veto", "extraelec_veto", "dilepton_veto", "fj_Xbb_pt", "fj_Xbb_eta", "fj_Xbb_particleNet_XbbvsQCD", "bpair_pt_1", "bpair_pt_2", "bpair_btag_value_2", "bpair_eta_2", "met", "mass_2", "tau_decaymode_2"],
        "tt": ["fj_Xbb_hadflavor", "fj_Xbb_nBhad", "fj_Xbb_nChad", "weight", "btag_weight", "pNet_Xbb_weight", "njets", "nbtag", "q_1", "pt_2", "q_2", "gen_match_1", "gen_match_2", "m_vis", "mt_1", "deltaR_ditaupair", "pt_1", "iso_1", "metphi", "extramuon_veto", "extraelec_veto", "fj_Xbb_pt", "fj_Xbb_eta", "fj_Xbb_particleNet_XbbvsQCD", "bpair_pt_1", "bpair_pt_2", "bpair_btag_value_2", "bpair_eta_2", "met", "mass_1", "tau_decaymode_1", "mass_2", "tau_decaymode_2"],
    },
    "nmssm_boosted": {
        "et": ["fj_Xbb_hadflavor", "fj_Xbb_nBhad", "fj_Xbb_nChad", "weight", "btag_weight", "pNet_Xbb_weight", "njets", "nbtag", "q_1", "pt_2", "q_2", "gen_match_2", "m_vis", "mt_1", "deltaR_ditaupair", "pt_1", "iso_1", "metphi", "extramuon_veto", "extraelec_veto", "dilepton_veto", "fj_Xbb_pt", "fj_Xbb_eta", "fj_Xbb_particleNet_XbbvsQCD", "bpair_pt_1", "bpair_pt_2", "bpair_btag_value_2", "bpair_eta_2", "met", "mass_2", "tau_decaymode_2"],
        "mt": ["fj_Xbb_hadflavor", "fj_Xbb_nBhad", "fj_Xbb_nChad", "weight", "btag_weight", "pNet_Xbb_weight", "njets", "nbtag", "q_1", "pt_2", "q_2", "gen_match_2", "m_vis", "mt_1", "deltaR_ditaupair", "pt_1", "iso_1", "metphi", "extramuon_veto", "extraelec_veto", "dilepton_veto", "fj_Xbb_pt", "fj_Xbb_eta", "fj_Xbb_particleNet_XbbvsQCD", "bpair_pt_1", "bpair_pt_2", "bpair_btag_value_2", "bpair_eta_2", "met", "mass_2", "tau_decaymode_2"],
        "tt": ["fj_Xbb_hadflavor", "fj_Xbb_nBhad", "fj_Xbb_nChad", "weight", "btag_weight", "pNet_Xbb_weight", "njets", "nbtag", "q_1", "pt_2", "q_2", "gen_match_1", "gen_match_2", "m_vis", "mt_1", "deltaR_ditaupair", "pt_1", "iso_1", "metphi", "extramuon_veto", "extraelec_veto", "fj_Xbb_pt", "fj_Xbb_eta", "fj_Xbb_particleNet_XbbvsQCD", "bpair_pt_1", "bpair_pt_2", "bpair_btag_value_2", "bpair_eta_2", "met", "mass_1", "tau_decaymode_1", "mass_2", "tau_decaymode_2"],
    },
    "smhtt_ul": {
        "et": ["weight", "btag_weight", "njets", "nbtag", "q_1", "pt_2", "q_2", "gen_match_2", "m_vis", "mt_1", "deltaR_ditaupair", "pt_1", "iso_1", "metphi", "extramuon_veto", "extraelec_veto", "tau_decaymode_2"],
        "mt": ["weight", "btag_weight", "njets", "nbtag", "q_1", "pt_2", "q_2", "gen_match_2", "m_vis", "mt_1", "deltaR_ditaupair", "pt_1", "iso_1", "metphi", "extramuon_veto", "extraelec_veto", "dimuon_veto", "tau_decaymode_2"],
        "tt": ["weight", "btag_weight", "njets", "nbtag", "q_1", "pt_2", "q_2", "gen_match_1", "gen_match_2", "m_vis", "mt_1", "deltaR_ditaupair", "pt_1", "iso_1", "metphi", "extramuon_veto", "extraelec_veto", "tau_decaymode_2", "tau_decaymode_1"],
    },
}

### For plotting ###

# label definitions for y-axis
FF_YAxis = {
    "ttbar": "FF_{t#bar{t}}",
    "ttbar_subleading": "FF_{t#bar{t}}",
    "Wjets": "FF_{Wjets}",
    "QCD": "FF_{QCD}",
    "QCD_subleading": "FF_{QCD}",
}
# definitions for channels
channel_dict = {"et": "e#tau_{h}", "mt": "#mu#tau_{h}", "tt": "#tau_{h}#tau_{h}"}
# definitions for era and luminosity TODO: 2016
era_dict = {
    "2016preVFP": "19.52 fb^{-1} (2016preVFP, 13 TeV)",
    "2016postVFP": "16.81 fb^{-1} (2016postVFP, 13 TeV)",
    "2017": "41.48 fb^{-1} (2017, 13 TeV)",
    "2018": "59.83 fb^{-1} (2018, 13 TeV)",
}
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
    "fit_graph_unct": ROOT.kRed,
}
# definitions for process labels on the plots
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
    "fit_graph_unct": "best fit uncertainty",
    "fit_graph_mc_sub": "best fit (MC subtraction unc.)",
}
# definitions to translate variable to readable language, channel dependent
variable_dict = {
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
# definitions to translate category cuts to readable language
category_dict = {
    "incl": "incl.",
    "njets": "N_{jets}",
    "nbtag": "N_{b-jets}",
    "deltaR_ditaupair": "#Delta" + "R(l#tau_{h})",
    "tau_decaymode_1": "#tau_{h, 1}^{DM}",
    "tau_decaymode_2": "#tau_{h, 2}^{DM}",
}

### For correctionlib ###

# intern naming translation helper for variables
variable_translator = {
    "pt_2": "subleading_lep_pt",
    "pt_1": "leading_lep_pt",
    "mass_1": "leading_lep_mass",
    "mass_2": "subleading_lep_mass",
    "mt_1": "lep_mt",
    "iso_1": "lep_iso",
    "m_vis": "m_vis",
    "bpt_1": "b_pt",
    "njets": "njets",
    "nbtag": "nbtags",
    "deltaR_ditaupair": "dR_tau_pair",
    "QCD": "QCD",
    "Wjets": "Wjets",
    "ttbar_J": "ttbar",
    "tau_decaymode_1": "leading_tau_decaymode",
    "tau_decaymode_2": "subleading_tau_decaymode",
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
    "bpt_1": "transverse momentum of the hardest b-tagged jet; measured between #var_min and #var_max GeV; for higher/lower pt's the edge values are used",
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
# intern naming translation helper for correction uncertainties
corr_variation_dict = {
    "non_closure_subleading_lep_pt": "nonClosureSubleadingLepPt",
    "non_closure_leading_lep_pt": "nonClosureLeadingLepPt",
    "non_closure_m_vis": "nonClosureMvis",
    "non_closure_b_pt": "nonClosureBPt",
    "non_closure_lep_iso": "nonClosureLepIso",
    "non_closure_subleading_lep_mass": "nonClosureSubleadingLepMass",
    "non_closure_leading_lep_mass": "nonClosureLeadingLepMass",
    "non_closure_lep_mt": "nonClosureLepMT",
    "DR_SR": "DRtoSR",
}
