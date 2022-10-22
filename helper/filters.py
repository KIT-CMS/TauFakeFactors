import sys


def had_tau_pt_cut(rdf, channel):
    if channel == "et":
        rdf = rdf.Filter("(pt_2 > 30)", "had. tau pT cut")
    elif channel == "mt":
        rdf = rdf.Filter("(pt_2 > 30)", "had. tau pT cut")
    elif channel == "tt":
        rdf = rdf.Filter("(pt_1 > 40) && (pt_2 > 40)", "had. tau pT cuts")
    else:
        sys.exit("Eventfilter: tau pt: Such a channel is not defined: {}".format(channel))
    return rdf

def had_tau_eta_cut(rdf, channel):
    if channel == "et":
        rdf = rdf.Filter("(abs(eta_2) < 2.3)", "had. tau eta cut")
    elif channel == "mt":
        rdf = rdf.Filter("(abs(eta_2) < 2.3)", "had. tau eta cut")
    elif channel == "tt":
        rdf = rdf.Filter(
            "(abs(eta_1) < 2.1) && (abs(eta_2) < 2.1)", "had. tau eta cuts"
        )
    else:
        sys.exit("Eventfilter: tau eta: Such a channel is not defined: {}".format(channel))
    return rdf

def had_tau_decay_mode_cut(rdf, channel):
    if channel == "et":
        rdf = rdf.Filter(
            "(decaymode_2 == 0) || (decaymode_2 == 1) || (decaymode_2 == 10) || (decaymode_2 == 11)",
            "had. tau decay mode cuts",
        )
    elif channel == "mt":
        rdf = rdf.Filter(
            "(decaymode_2 == 0) || (decaymode_2 == 1) || (decaymode_2 == 10) || (decaymode_2 == 11)",
            "had. tau decay mode cuts",
        )
    elif channel == "tt":
        rdf = rdf.Filter(
            "((decaymode_1 == 0) || (decaymode_1 == 1) || (decaymode_1 == 10) || (decaymode_1 == 11)) && ((decaymode_2 == 0) || (decaymode_2 == 1) || (decaymode_2 == 10) || (decaymode_2 == 11))",
            "had. tau decay mode cuts",
        )
    else:
        sys.exit("Eventfilter: tau dm: Such a channel is not defined: {}".format(channel))
    return rdf

def had_tau_vsLep_id_cut(rdf, channel):
    if channel == "et":
        rdf = rdf.Filter(
            "(id_tau_vsMu_VLoose_2 > 0.5) && (id_tau_vsEle_Tight_2 > 0.5)",
            "tau vs leptons id cuts",
        )
    elif channel == "mt":
        rdf = rdf.Filter(
            "(id_tau_vsMu_Tight_2 > 0.5) && (id_tau_vsEle_VVLoose_2 > 0.5)",
            "tau vs leptons id cuts",
        )
    elif channel == "tt":
        rdf = rdf.Filter(
            "(id_tau_vsMu_VLoose_1 > 0.5) && (id_tau_vsEle_VVLoose_1 > 0.5) && (id_tau_vsMu_VLoose_2 > 0.5) && (id_tau_vsEle_VVLoose_2 > 0.5)",
            "tau vs leptons id cuts",
        )
    else:
        sys.exit("Eventfilter: tau id: Such a channel is not defined: {}".format(channel))
    return rdf

def lep_iso_cut(rdf, channel):
    if channel == "et":
        rdf = rdf.Filter(
            "iso_1 < 0.15",
            "electron isolation cut",
        )
    elif channel == "mt":
        rdf = rdf.Filter(
            "iso_1 < 0.15",
            "muon isolation cut",
        )
    elif channel == "tt":
        pass
    else:
        sys.exit("Eventfilter: extra lep veto: Such a channel is not defined: {}".format(channel))
    return rdf

def extra_leptons_veto(rdf, channel):
    if channel == "et":
        rdf = rdf.Filter(
            "(extramuon_veto < 0.5) && (extraelec_veto < 0.5) && (dimuon_veto < 0.5)",
            "extra leptons veto",
        )
    elif channel == "mt":
        rdf = rdf.Filter(
            "(extramuon_veto < 0.5) && (extraelec_veto < 0.5) && (dimuon_veto < 0.5)",
            "extra leptons veto",
        )
    elif channel == "tt":
        rdf = rdf.Filter(
            "(extramuon_veto < 0.5) && (extraelec_veto < 0.5) && (dimuon_veto < 0.5)",
            "extra leptons veto",
        )
    else:
        sys.exit("Eventfilter: extra lep veto: Such a channel is not defined: {}".format(channel))
    return rdf

def trigger_cut(rdf, channel):
    if channel == "et":
        rdf = rdf.Filter(
            "((pt_1 > 33) && ((trg_single_ele32 > 0.5) || (trg_single_ele35 > 0.5))) || ((pt_2 > 35) && (abs(eta_2) < 2.1) && (pt_1 <= 33) && (pt_1 > 25) && ((trg_cross_ele24tau30 > 0.5) || (trg_cross_ele24tau30_hps > 0.5)))",
            "single electron or cross trigger + electron pT cuts",
        )
    elif channel == "mt":
        rdf = rdf.Filter(
            "((pt_1 > 25) && ((trg_single_mu24 > 0.5) || (trg_single_mu27 > 0.5))) || ((pt_2 > 32) && (pt_1 <= 25) && (pt_1 > 21) && ((trg_cross_mu20tau27 > 0.5) || (trg_cross_mu20tau27_hps > 0.5)))",
            "single muon or cross trigger + muon pT cuts",
        )
    elif channel == "tt":
        rdf = rdf.Filter(
            "(trg_double_tau35_tightiso_tightid > 0.5) || (trg_double_tau35_mediumiso_hps > 0.5) || (trg_double_tau40_mediumiso_tightid > 0.5) || (trg_double_tau40_tightiso > 0.5)",
            "tau double trigger cuts",
        )
    else:
        sys.exit("Eventfilter: trigger: Such a channel is not defined: {}".format(channel))
    return rdf

def tau_gen_match_split(rdf, channel, tau_gen_mode):
    if channel == "et":
        if tau_gen_mode == "all":
            pass
        elif tau_gen_mode == "T":
            rdf = rdf.Filter(
                "(gen_match_1 == 3) && (gen_match_2 == 5)",
                "tau gen. match split cuts",
            )
        elif tau_gen_mode == "J":
            rdf = rdf.Filter(
                "(gen_match_2 == 6)",
                "tau gen. match split cuts",
            )
        elif tau_gen_mode == "L":
            rdf = rdf.Filter(
                "((gen_match_1 != 3) && (gen_match_2 != 5)) && (gen_match_2 != 6)",
                "tau gen. match split cuts",
            )
        else:
            print("Eventfilter: tau gen match: et: Such a tau gen. match is not defined: {}".format(tau_gen_mode))

    elif channel == "mt":
        if tau_gen_mode == "all":
            pass
        elif tau_gen_mode == "T":
            rdf = rdf.Filter(
                "(gen_match_1 == 4) && (gen_match_2 == 5)",
                "tau gen. match split cuts",
            )
        elif tau_gen_mode == "J":
            rdf = rdf.Filter(
                "(gen_match_2 == 6)",
                "tau gen. match split cuts",
            )
        elif tau_gen_mode == "L":
            rdf = rdf.Filter(
                "((gen_match_1 != 4) && (gen_match_2 != 5)) && (gen_match_2 != 6)",
                "tau gen. match split cuts",
            )
        else:
            print("Eventfilter: tau gen match: mt: Such a tau gen. match is not defined: {}".format(tau_gen_mode))

    elif channel == "tt":
        if tau_gen_mode == "all":
            pass
        elif tau_gen_mode == "T":
            rdf = rdf.Filter(
                "(gen_match_1 == 5) && (gen_match_2 == 5)",
                "tau gen. match split cuts",
            )
        elif tau_gen_mode == "J":
            rdf = rdf.Filter(
                "(gen_match_1 == 6) || (gen_match_2 == 6)",
                "tau gen. match split cuts",
            )
        elif tau_gen_mode == "L":
            rdf = rdf.Filter(
                "!((gen_match_1 == 5) && (gen_match_2 == 5)) && !((gen_match_1 == 6) || (gen_match_2 == 6))",
                "tau gen. match split cuts",
            )
        else:
            print("Eventfilter: tau gen match: tt: Such a tau gen. match is not defined: {}".format(tau_gen_mode))

    else:
        sys.exit("Eventfilter: tau gen match: Such a channel is not defined: {}".format(channel))

    return rdf

def emb_tau_gen_match(rdf, channel):
    if channel == "et":
        rdf = rdf.Filter(
            "(gen_match_1 == 3) && (gen_match_2 == 5)",
            "embedding tau gen. matching",
        )
    elif channel == "mt":
        rdf = rdf.Filter(
            "(gen_match_1 == 4) && (gen_match_2 == 5)",
            "embedding tau gen. matching",
        )
    elif channel == "tt":
        rdf = rdf.Filter(
            "(gen_match_1 == 5) && (gen_match_2 == 5)",
            "embedding tau gen. matching",
        )
    else:
        sys.exit("Eventfilter: emb tau gen match: Such a channel is not defined: {}".format(channel))
    return rdf