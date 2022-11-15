import sys


def had_tau_pt_cut(rdf, channel, config):
    if channel == "et":
        rdf = rdf.Filter(
            "(pt_2 {})".format(config["had_tau_pt"]),
            "cut had. tau pT {}".format(config["had_tau_pt"]),
        )
    elif channel == "mt":
        rdf = rdf.Filter(
            "(pt_2 {})".format(config["had_tau_pt"]),
            "cut had. tau pT {}".format(config["had_tau_pt"]),
        )
    elif channel == "tt":
        rdf = rdf.Filter(
            "(pt_1 {}) && (pt_2 {})".format(config["had_tau_pt"], config["had_tau_pt"]),
            "cut had. tau pT {}".format(config["had_tau_pt"]),
        )
    else:
        sys.exit(
            "Eventfilter: tau pt: Such a channel is not defined: {}".format(channel)
        )
    return rdf


def had_tau_eta_cut(rdf, channel, config):
    if channel == "et":
        rdf = rdf.Filter(
            "(abs(eta_2) {})".format(config["had_tau_eta"]),
            "cut had. tau abs(eta) {}".format(config["had_tau_eta"]),
        )
    elif channel == "mt":
        rdf = rdf.Filter(
            "(abs(eta_2) {})".format(config["had_tau_eta"]),
            "cut had. tau abs(eta) {}".format(config["had_tau_eta"]),
        )
    elif channel == "tt":
        rdf = rdf.Filter(
            "(abs(eta_1) {}) && (abs(eta_2) {})".format(
                config["had_tau_eta"], config["had_tau_eta"]
            ),
            "cut had. tau abs(eta) {}".format(config["had_tau_eta"]),
        )
    else:
        sys.exit(
            "Eventfilter: tau eta: Such a channel is not defined: {}".format(channel)
        )
    return rdf


def had_tau_decay_mode_cut(rdf, channel, config):
    cut_string = ""
    for i, dm in enumerate(config["had_tau_decay_mode"]):
        if i == 0:
            cut_string += "(decaymode_2 == {})".format(dm)
        else:
            cut_string += "|| (decaymode_2 == {})".format(dm)

    if channel == "et":
        rdf = rdf.Filter(
            cut_string,
            "cut had. tau decay modes: {}".format(
                " ".join(dm for dm in config["had_tau_decay_mode"])
            ),
        )
    elif channel == "mt":
        rdf = rdf.Filter(
            cut_string,
            "cut had. tau decay modes: {}".format(
                " ".join(dm for dm in config["had_tau_decay_mode"])
            ),
        )
    elif channel == "tt":
        rdf = rdf.Filter(
            cut_string,
            "cut second had. tau decay modes: {}".format(
                " ".join(dm for dm in config["had_tau_decay_mode"])
            ),
        )
        rdf = rdf.Filter(
            cut_string.replace("decaymode_2", "decaymode_1"),
            "cut first had. tau decay modes: {}".format(
                " ".join(dm for dm in config["had_tau_decay_mode"])
            ),
        )
    else:
        sys.exit(
            "Eventfilter: tau decay mode: Such a channel is not defined: {}".format(
                channel
            )
        )
    return rdf


def had_tau_vsLep_id_cut(rdf, channel, config):
    if channel == "et":
        rdf = rdf.Filter(
            "(id_tau_vsMu_{}_2 > 0.5) && (id_tau_vsEle_{}_2 > 0.5)".format(
                config["had_tau_id_vs_mu"], config["had_tau_id_vs_ele"]
            ),
            "tau {} vsEle and {} vsMu id cuts".format(
                config["had_tau_id_vs_ele"], config["had_tau_id_vs_mu"]
            ),
        )
    elif channel == "mt":
        rdf = rdf.Filter(
            "(id_tau_vsMu_{}_2 > 0.5) && (id_tau_vsEle_{}_2 > 0.5)".format(
                config["had_tau_id_vs_mu"], config["had_tau_id_vs_ele"]
            ),
            "tau {} vsEle and {} vsMu id cuts".format(
                config["had_tau_id_vs_ele"], config["had_tau_id_vs_mu"]
            ),
        )
    elif channel == "tt":
        rdf = rdf.Filter(
            "(id_tau_vsMu_{}_1 > 0.5) && (id_tau_vsEle_{}_1 > 0.5) && (id_tau_vsMu_{}_2 > 0.5) && (id_tau_vsEle_{}_2 > 0.5)".format(
                config["had_tau_id_vs_mu"],
                config["had_tau_id_vs_ele"],
                config["had_tau_id_vs_mu"],
                config["had_tau_id_vs_ele"],
            ),
            "tau {} vsEle and {} vsMu id cuts".format(
                config["had_tau_id_vs_ele"], config["had_tau_id_vs_mu"]
            ),
        )
    else:
        sys.exit(
            "Eventfilter: tau id vs leptons: Such a channel is not defined: {}".format(
                channel
            )
        )
    return rdf


def lep_iso_cut(rdf, channel, config):
    if channel == "et":
        rdf = rdf.Filter(
            "iso_1 {}".format(config["lep_iso"]),
            "cut electron isolation {}".format(config["lep_iso"]),
        )
    elif channel == "mt":
        rdf = rdf.Filter(
            "iso_1 {}".format(config["lep_iso"]),
            "cut muon isolation {}".format(config["lep_iso"]),
        )
    elif channel == "tt":
        pass
    else:
        sys.exit(
            "Eventfilter: lepton isolation: Such a channel is not defined: {}".format(
                channel
            )
        )
    return rdf


def trigger_cut(rdf, channel):
    if channel == "et":
        rdf = rdf.Filter(
            "((pt_1 > 33) && ((trg_single_ele32 > 0.5) || (trg_single_ele35 > 0.5)))",  # || ((pt_2 > 35) && (abs(eta_2) < 2.1) && (pt_1 <= 33) && (pt_1 > 25) && ((trg_cross_ele24tau30 > 0.5) || (trg_cross_ele24tau30_hps > 0.5)))",
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
        sys.exit(
            "Eventfilter: trigger: Such a channel is not defined: {}".format(channel)
        )
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
                "(!(gen_match_1==3 && gen_match_2==5)) && (!(gen_match_2 == 6))",
                "tau gen. match split cuts",
            )
        else:
            print(
                "Eventfilter: tau gen match: et: Such a tau gen. match is not defined: {}".format(
                    tau_gen_mode
                )
            )

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
            print(
                "Eventfilter: tau gen match: mt: Such a tau gen. match is not defined: {}".format(
                    tau_gen_mode
                )
            )

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
            print(
                "Eventfilter: tau gen match: tt: Such a tau gen. match is not defined: {}".format(
                    tau_gen_mode
                )
            )

    else:
        sys.exit(
            "Eventfilter: tau gen match: Such a channel is not defined: {}".format(
                channel
            )
        )

    return rdf


def emb_tau_gen_match(rdf, channel):
    if channel == "et":
        rdf = rdf.Filter(
            # "(gen_match_1 == 3) && (gen_match_2 == 5)",
            "((gen_match_1 > 2 && gen_match_1 < 6) && (gen_match_2 > 2 && gen_match_2 < 6))",
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
        sys.exit(
            "Eventfilter: emb tau gen match: Such a channel is not defined: {}".format(
                channel
            )
        )
    return rdf
