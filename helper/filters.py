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


def had_tau_id_vsEle_cut(rdf, channel, config):
    if channel == "et":
        rdf = rdf.Filter(
            "(id_tau_vsEle_{WP}_2 > 0.5)".format(WP=config["had_tau_id_vs_ele"]),
            "cut on {WP} tau vs electrons id".format(WP=config["had_tau_id_vs_ele"]),
        )
    elif channel == "mt":
        rdf = rdf.Filter(
            "(id_tau_vsEle_{WP}_2 > 0.5)".format(WP=config["had_tau_id_vs_ele"]),
            "cut on {WP} tau vs electrons id".format(WP=config["had_tau_id_vs_ele"]),
        )
    elif channel == "tt":
        rdf = rdf.Filter(
            "(id_tau_vsEle_{WP}_1 > 0.5) && (id_tau_vsEle_{WP}_2 > 0.5)".format(
                WP=config["had_tau_id_vs_ele"]
            ),
            "cut on {WP} tau vs electrons id".format(
                WP=config["had_tau_id_vs_ele"]
            ),
        )
    else:
        sys.exit(
            "Eventfilter: tau id vs electron: Such a channel is not defined: {}".format(
                channel
            )
        )
    return rdf

def had_tau_id_vsMu_cut(rdf, channel, config):
    if channel == "et":
        rdf = rdf.Filter(
            "(id_tau_vsMu_{WP}_2 > 0.5)".format(WP=config["had_tau_id_vs_mu"]),
            "cut on {WP} tau vs muons id".format(WP=config["had_tau_id_vs_mu"]),
        )
    elif channel == "mt":
        rdf = rdf.Filter(
            "(id_tau_vsMu_{WP}_2 > 0.5)".format(WP=config["had_tau_id_vs_mu"]),
            "cut on {WP} tau vs muons id".format(WP=config["had_tau_id_vs_mu"]),
        )
    elif channel == "tt":
        rdf = rdf.Filter(
            "(id_tau_vsMu_{WP}_1 > 0.5) && (id_tau_vsMu_{WP}_2 > 0.5)".format(
                WP=config["had_tau_id_vs_mu"],
            ),
            "cut on {WP} tau vs muons id".format(
                WP=config["had_tau_id_vs_mu"]
            ),
        )
    else:
        sys.exit(
            "Eventfilter: tau id vs muon: Such a channel is not defined: {}".format(
                channel
            )
        )
    return rdf


def had_tau_id_vsJet_cut(rdf, channel, config):
    if isinstance(config["had_tau_id_vs_jet"], str):
        if channel == "et":
            rdf = rdf.Filter(
                "(id_tau_vsJet_{WP}_2 > 0.5)".format(WP=config["had_tau_id_vs_jet"]),
                "cut on {WP} tau vs jets id".format(WP=config["had_tau_id_vs_jet"]),
            )
        elif channel == "mt":
            rdf = rdf.Filter(
                "(id_tau_vsJet_{WP}_2 > 0.5)".format(WP=config["had_tau_id_vs_jet"]),
                "cut on {WP} tau vs jets id".format(WP=config["had_tau_id_vs_jet"]),
            )
        elif channel == "tt":
            rdf = rdf.Filter(
                "(id_tau_vsJet_{WP}_2 > 0.5)".format(WP=config["had_tau_id_vs_jet"]),
                "cut on {WP} tau vs jets id".format(WP=config["had_tau_id_vs_jet"]),
            )
        else:
            sys.exit(
                "Eventfilter: tau id vs jet: Such a channel is not defined: {}".format(
                    channel
                )
            )

    elif isinstance(config["had_tau_id_vs_jet"], list):
        if channel == "et":
            rdf = rdf.Filter(
                "(id_tau_vsJet_{lower_WP}_2 > 0.5) && (id_tau_vsJet_{upper_WP}_2 < 0.5)".format(
                    lower_WP=config["had_tau_id_vs_jet"][0], 
                    upper_WP=config["had_tau_id_vs_jet"][1],
                ),
                "cut on tau vs jets id between {lower_WP} and {upper_WP}".format(
                    lower_WP=config["had_tau_id_vs_jet"][0], 
                    upper_WP=config["had_tau_id_vs_jet"][1],
                ),
            )
        elif channel == "mt":
            rdf = rdf.Filter(
                "(id_tau_vsJet_{upper_WP}_2 > 0.5) && (id_tau_vsJet_{upper_WP}_2 < 0.5)".format(
                    lower_WP=config["had_tau_id_vs_jet"][0], 
                    upper_WP=config["had_tau_id_vs_jet"][1],
                ),
                "cut on tau vs jets id between {lower_WP} and {upper_WP}".format(
                    lower_WP=config["had_tau_id_vs_jet"][0], 
                    upper_WP=config["had_tau_id_vs_jet"][1],
                ),
            )
        elif channel == "tt":
            rdf = rdf.Filter(
                "(id_tau_vsJet_{lower_WP}_2 > 0.5) && (id_tau_vsJet_{upper_WP}_2 < 0.5)".format(
                    lower_WP=config["had_tau_id_vs_jet"][0], 
                    upper_WP=config["had_tau_id_vs_jet"][1],
                ),
                "cut on tau vs jets id between {lower_WP} and {upper_WP}".format(
                    lower_WP=config["had_tau_id_vs_jet"][0], 
                    upper_WP=config["had_tau_id_vs_jet"][1],
                ),
            )
        else:
            sys.exit(
                "Eventfilter: tau id vs jet: Such a channel is not defined: {}".format(
                    channel
                )
            )
    else:
        sys.exit(
                "Eventfilter: tau id vs jet: Such a type is not defined: {}".format(
                    config["had_tau_id_vs_jet"]
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


def lep_mt_cut(rdf, channel, config):
    if channel == "et":
        rdf = rdf.Filter(
            "(mt_1 {})".format(config["lep_mt"]),
            "W boson origin cut on lepton mt {}".format(config["lep_mt"]),
        )
    elif channel == "mt":
        rdf = rdf.Filter(
            "(mt_1 {})".format(config["lep_mt"]),
            "W boson origin cut on lepton mT {}".format(config["lep_mt"]),
        )
    elif channel == "tt":
        rdf = rdf.Filter(
            "(mt_1 {})".format(config["lep_mt"]),
            "W boson origin cut on lepton mT {}".format(config["lep_mt"]),
        )
    else:
        sys.exit(
            "Eventfilter: lepton transverse mass: Such a channel is not defined: {}".format(
                channel
            )
        )
    return rdf


def no_extra_lep_cut(rdf, channel, config):
    no_lep = "> 0.5" if config["no_extra_lep"] else "< 0.5"
    if channel == "et":
        rdf = rdf.Filter(
            "(no_extra_lep {})".format(no_lep),
            "exclude additional leptons: {}".format(config["no_extra_lep"]),
        )
    elif channel == "mt":
        rdf = rdf.Filter(
            "(no_extra_lep {})".format(no_lep),
            "exclude additional leptons: {}".format(config["no_extra_lep"]),
        )
    elif channel == "tt":
        rdf = rdf.Filter(
            "(no_extra_lep {})".format(no_lep),
            "exclude additional leptons: {}".format(config["no_extra_lep"]),
        )
    else:
        sys.exit(
            "Eventfilter: no extra lep: Such a channel is not defined: {}".format(
                channel
            )
        )
    return rdf


def tau_pair_sign_cut(rdf, channel, sign):
    if sign == "same":
        sign_cut = "> 0"
    elif sign == "opposite":
        sign_cut = "< 0"
    else:
        sys.exit("Eventfilter: SS/OS: Such a sign is not defined: {}".format(sign))

    if channel == "et":
        rdf = rdf.Filter(
            "((q_1 * q_2) {})".format(sign_cut), "{} sign cut".format(sign)
        )
    elif channel == "mt":
        rdf = rdf.Filter(
            "((q_1 * q_2) {})".format(sign_cut), "{} sign cut".format(sign)
        )
    elif channel == "tt":
        rdf = rdf.Filter(
            "((q_1 * q_2) {})".format(sign_cut), "{} sign cut".format(sign)
        )
    else:
        sys.exit(
            "Eventfilter: SS/OS: Such a channel is not defined: {}".format(channel)
        )
    return rdf


def trigger_cut(rdf, channel):
    if channel == "et":
        rdf = rdf.Filter(
            "((pt_1 > 33) && ((trg_single_ele32 > 0.5) || (trg_single_ele35 > 0.5)))",  # || ((pt_2 > 35) && (abs(eta_2) < 2.1) && (pt_1 <= 33) && (pt_1 > 25) && ((trg_cross_ele24tau30 > 0.5) || (trg_cross_ele24tau30_hps > 0.5)))",
            "single electron trigger + electron pT cuts",
        )
    elif channel == "mt":
        rdf = rdf.Filter(
            "((pt_1 > 25) && ((trg_single_mu24 > 0.5) || (trg_single_mu27 > 0.5)))",  # || ((pt_2 > 32) && (pt_1 <= 25) && (pt_1 > 21) && ((trg_cross_mu20tau27 > 0.5) || (trg_cross_mu20tau27_hps > 0.5)))",
            "single muon trigger + muon pT cuts",
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


def bjet_number_cut(rdf, channel, nbtags):
    if channel == "et":
        rdf = rdf.Filter(
            "(nbtag {})".format(nbtags),
            "cut on {} b-tagged jets".format(nbtags),
        )
    elif channel == "mt":
        rdf = rdf.Filter(
            "(nbtag {})".format(nbtags),
            "cut on {} b-tagged jets".format(nbtags),
        )
    elif channel == "tt":
        rdf = rdf.Filter(
            "(nbtag {})".format(nbtags),
            "cut on {} b-tagged jets".format(nbtags),
        )
    else:
        sys.exit(
            "Eventfilter: bjet number: Such a channel is not defined: {}".format(
                channel
            )
        )
    return rdf


def jet_number_cut(rdf, channel, njets):
    if channel == "et":
        rdf = rdf.Filter("(njets {})".format(njets), "cut on {} jets".format(njets))
    elif channel == "mt":
        rdf = rdf.Filter("(njets {})".format(njets), "cut on {} jets".format(njets))
    elif channel == "tt":
        rdf = rdf.Filter("(njets {})".format(njets), "cut on {} jets".format(njets))
    else:
        sys.exit(
            "Eventfilter: jet number: Such a channel is not defined: {}".format(
                channel
            )
        )
    return rdf


def tau_origin_split(rdf, channel, tau_gen_mode):
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
