import sys


def emb_tau_gen_match(rdf, channel: str):
    """
    Function to make sure the tau pair final state objects originated for genuine taus for embedded events.

    Args:
        rdf: root DataFrame object
        channel: Analysis channel of the tau analysis e.g. "et", "mt" or "tt"

    Return:
        Filtered root DataFrame object
    """
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
        raise ValueError(
            f"Eventfilter: emb tau gen match: Such a channel is not defined: {channel}"
        )
    return rdf


def emb_boostedtau_gen_match(rdf, channel: str):
    """
    Function to make sure the boostedtau pair final state objects originated for genuine taus for embedded events.

    Args:
        rdf: root DataFrame object
        channel: Analysis channel of the tau analysis e.g. "et", "mt" or "tt"

    Return:
        Filtered root DataFrame object
    """
    if channel == "et":
        rdf = rdf.Filter(
            "(boosted_gen_match_1 == 3) && (boosted_gen_match_2 == 5)",
            "embedding boostedtau gen. matching",
        )
    elif channel == "mt":
        rdf = rdf.Filter(
            "(boosted_gen_match_1 == 4) && (boosted_gen_match_2 == 5)",
            "embedding boostedtau gen. matching",
        )
    elif channel == "tt":
        rdf = rdf.Filter(
            "(boosted_gen_match_1 == 5) && (boosted_gen_match_2 == 5)",
            "embedding boostedtau gen. matching",
        )
    else:
        raise ValueError(
            f"Eventfilter: emb boostedtau gen match: Such a channel is not defined: {channel}"
        )
    return rdf


def tau_origin_split(rdf, channel: str, tau_gen_mode: str):
    """
    Function to apply a cut based on the origin of the tau pair. This is needed to prevent an overestimation because
    genuine tau pairs are estimated with embedded events and events with jets faking hadronic taus are the targed the
    fake factor measurement. There are diffrent possibilities where a lepton (electron, muon, hadronic tau) originate from.
    Electrons and muons are either genuine or decay from a tau. Hadtronic taus are either decays from a genuine tau
    or jets/leptons misidentified as a hadronic tau.

    Args:
        rdf: root DataFrame object
        channel: Analysis channel of the tau analysis e.g. "et", "mt" or "tt"
        tau_gen_mode: There are 4 options:
                      "T" both taus in the tau pair originated from genuine taus,
                      "J" a hadronic tau is faked by a jet,
                      "L" everything else where mainly the muons or electrons are genuine,
                      "all" when no split should be applied

    Return:
        Filtered root DataFrame object
    """
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
            raise ValueError(
                f"Eventfilter: tau gen match: et: Such a tau gen. match is not defined: {tau_gen_mode}"
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
                "(!(gen_match_1==4 && gen_match_2==5)) && (!(gen_match_2 == 6))",
                "tau gen. match split cuts",
            )
        else:
            raise ValueError(
                f"Eventfilter: tau gen match: mt: Such a tau gen. match is not defined: {tau_gen_mode}"
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
                "(!((gen_match_1 == 5) && (gen_match_2 == 5)) && !((gen_match_1 == 6) || (gen_match_2 == 6)))",
                "tau gen. match split cuts",
            )
        else:
            raise ValueError(
                f"Eventfilter: tau gen match: tt: Such a tau gen. match is not defined: {tau_gen_mode}"
            )

    else:
        raise ValueError(
            f"Eventfilter: tau gen match: Such a channel is not defined: {channel}"
        )

    return rdf


def boostedtau_origin_split(rdf, channel: str, tau_gen_mode: str):
    """
    Function to apply a cut based on the origin of the boosted tau pair. This is needed to prevent an overestimation because
    genuine tau pairs are estimated with embedded events and events with jets faking hadronic taus are the targed the
    fake factor measurement. There are diffrent possibilities where a lepton (electron, muon, hadronic tau) originate from.
    Electrons and muons are either genuine or decay from a tau. Hadtronic taus are either decays from a genuine tau
    or jets/leptons misidentified as a hadronic tau.

    Args:
        rdf: root DataFrame object
        channel: Analysis channel of the tau analysis e.g. "et", "mt" or "tt"
        tau_gen_mode: There are 4 options:
                      "T" both taus in the boosted tau pair originated from genuine taus,
                      "J" a hadronic tau is faked by a jet,
                      "L" everything else where mainly the muons or electrons are genuine,
                      "all" when no split should be applied

    Return:
        Filtered root DataFrame object
    """
    if channel == "et":
        if tau_gen_mode == "all":
            pass
        elif tau_gen_mode == "T":
            rdf = rdf.Filter(
                "(boosted_gen_match_1 == 3) && (boosted_gen_match_2 == 5)",
                "boosted tau gen. match split cuts",
            )
        elif tau_gen_mode == "J":
            rdf = rdf.Filter(
                "(boosted_gen_match_2 == 6)",
                "boosted tau gen. match split cuts",
            )
        elif tau_gen_mode == "L":
            rdf = rdf.Filter(
                "(!(boosted_gen_match_1==3 && boosted_gen_match_2==5)) && (!(boosted_gen_match_2 == 6))",
                "boosted tau gen. match split cuts",
            )
        else:
            raise ValueError(
                f"Eventfilter: boosted tau gen match: et: Such a tau gen. match is not defined: {tau_gen_mode}"
            )

    elif channel == "mt":
        if tau_gen_mode == "all":
            pass
        elif tau_gen_mode == "T":
            rdf = rdf.Filter(
                "(boosted_gen_match_1 == 4) && (boosted_gen_match_2 == 5)",
                "boosted tau gen. match split cuts",
            )
        elif tau_gen_mode == "J":
            rdf = rdf.Filter(
                "(boosted_gen_match_2 == 6)",
                "boosted tau gen. match split cuts",
            )
        elif tau_gen_mode == "L":
            rdf = rdf.Filter(
                "(!(boosted_gen_match_1==4 && boosted_gen_match_2==5)) && (!(boosted_gen_match_2 == 6))",
                "boosted tau gen. match split cuts",
            )
        else:
            raise ValueError(
                f"Eventfilter: boosted tau gen match: mt: Such a tau gen. match is not defined: {tau_gen_mode}"
            )

    elif channel == "tt":
        if tau_gen_mode == "all":
            pass
        elif tau_gen_mode == "T":
            rdf = rdf.Filter(
                "(boosted_gen_match_1 == 5) && (boosted_gen_match_2 == 5)",
                "boosted tau gen. match split cuts",
            )
        elif tau_gen_mode == "J":
            rdf = rdf.Filter(
                "(boosted_gen_match_1 == 6) || (boosted_gen_match_2 == 6)",
                "boosted tau gen. match split cuts",
            )
        elif tau_gen_mode == "L":
            rdf = rdf.Filter(
                "(!((boosted_gen_match_1 == 5) && (boosted_gen_match_2 == 5)) && !((boosted_gen_match_1 == 6) || (boosted_gen_match_2 == 6)))",
                "boosted tau gen. match split cuts",
            )
        else:
            raise ValueError(
                f"Eventfilter: boosted tau gen match: tt: Such a tau gen. match is not defined: {tau_gen_mode}"
            )

    else:
        raise ValueError(
            f"Eventfilter: tau gen match: Such a channel is not defined: {channel}"
        )

    return rdf


# def had_tau_pt_cut(rdf, channel, config):
#     if channel == "et":
#         rdf = rdf.Filter(
#             "(pt_2 {})".format(config["had_tau_pt"]),
#             "cut had. tau pT {}".format(config["had_tau_pt"]),
#         )
#     elif channel == "mt":
#         rdf = rdf.Filter(
#             "(pt_2 {})".format(config["had_tau_pt"]),
#             "cut had. tau pT {}".format(config["had_tau_pt"]),
#         )
#     elif channel == "tt":
#         rdf = rdf.Filter(
#             "(pt_1 {}) && (pt_2 {})".format(config["had_tau_pt"], config["had_tau_pt"]),
#             "cut had. tau pT {}".format(config["had_tau_pt"]),
#         )
#     else:
#         sys.exit(
#             "Eventfilter: tau pt: Such a channel is not defined: {}".format(channel)
#         )
#     return rdf


# def had_boostedtau_pt_cut(rdf, channel, config):
#     if channel == "et":
#         rdf = rdf.Filter(
#             "(boosted_pt_2 {})".format(config["had_boostedtau_pt"]),
#             "cut had. boosted tau pT {}".format(config["had_boostedtau_pt"]),
#         )
#     elif channel == "mt":
#         rdf = rdf.Filter(
#             "(boosted_pt_2 {})".format(config["had_boostedtau_pt"]),
#             "cut had. boosted tau pT {}".format(config["had_boostedtau_pt"]),
#         )
#     elif channel == "tt":
#         rdf = rdf.Filter(
#             "(boosted_pt_1 {}) && (boosted_pt_2 {})".format(config["had_boostedtau_pt"], config["had_boostedtau_pt"]),
#             "cut had. boosted tau pT {}".format(config["had_boostedtau_pt"]),
#         )
#     else:
#         sys.exit(
#             "Eventfilter: boosted tau pt: Such a channel is not defined: {}".format(channel)
#         )
#     return rdf


# def had_tau_eta_cut(rdf, channel, config):
#     if channel == "et":
#         rdf = rdf.Filter(
#             "(abs(eta_2) {})".format(config["had_tau_eta"]),
#             "cut had. tau abs(eta) {}".format(config["had_tau_eta"]),
#         )
#     elif channel == "mt":
#         rdf = rdf.Filter(
#             "(abs(eta_2) {})".format(config["had_tau_eta"]),
#             "cut had. tau abs(eta) {}".format(config["had_tau_eta"]),
#         )
#     elif channel == "tt":
#         rdf = rdf.Filter(
#             "(abs(eta_1) {}) && (abs(eta_2) {})".format(
#                 config["had_tau_eta"], config["had_tau_eta"]
#             ),
#             "cut had. tau abs(eta) {}".format(config["had_tau_eta"]),
#         )
#     else:
#         sys.exit(
#             "Eventfilter: tau eta: Such a channel is not defined: {}".format(channel)
#         )
#     return rdf


# def had_boostedtau_eta_cut(rdf, channel, config):
#     if channel == "et":
#         rdf = rdf.Filter(
#             "(abs(boosted_eta_2) {})".format(config["had_boostedtau_eta"]),
#             "cut had. boosted tau abs(eta) {}".format(config["had_boostedtau_eta"]),
#         )
#     elif channel == "mt":
#         rdf = rdf.Filter(
#             "(abs(boosted_eta_2) {})".format(config["had_boostedtau_eta"]),
#             "cut had. boosted tau abs(eta) {}".format(config["had_boostedtau_eta"]),
#         )
#     elif channel == "tt":
#         rdf = rdf.Filter(
#             "(abs(boosted_eta_1) {}) && (abs(boosted_eta_2) {})".format(
#                 config["had_boostedtau_eta"], config["had_boostedtau_eta"]
#             ),
#             "cut had. boosted tau abs(eta) {}".format(config["had_boostedtau_eta"]),
#         )
#     else:
#         sys.exit(
#             "Eventfilter: boosted tau eta: Such a channel is not defined: {}".format(channel)
#         )
#     return rdf


# def had_tau_decay_mode_cut(rdf, channel, config):
#     cut_string = ""
#     for i, dm in enumerate(config["had_tau_decay_mode"]):
#         if i == 0:
#             cut_string += "(tau_decaymode_2 == {})".format(dm)
#         else:
#             cut_string += "|| (tau_decaymode_2 == {})".format(dm)

#     if channel == "et":
#         rdf = rdf.Filter(
#             cut_string,
#             "cut had. tau decay modes: {}".format(
#                 " ".join(dm for dm in config["had_tau_decay_mode"])
#             ),
#         )
#     elif channel == "mt":
#         rdf = rdf.Filter(
#             cut_string,
#             "cut had. tau decay modes: {}".format(
#                 " ".join(dm for dm in config["had_tau_decay_mode"])
#             ),
#         )
#     elif channel == "tt":
#         rdf = rdf.Filter(
#             cut_string,
#             "cut second had. tau decay modes: {}".format(
#                 " ".join(dm for dm in config["had_tau_decay_mode"])
#             ),
#         )
#         rdf = rdf.Filter(
#             cut_string.replace("tau_decaymode_2", "tau_decaymode_1"),
#             "cut first had. tau decay modes: {}".format(
#                 " ".join(dm for dm in config["had_tau_decay_mode"])
#             ),
#         )
#     else:
#         sys.exit(
#             "Eventfilter: tau decay mode: Such a channel is not defined: {}".format(
#                 channel
#             )
#         )
#     return rdf


# def had_boostedtau_decay_mode_cut(rdf, channel, config):
#     cut_string = ""
#     for i, dm in enumerate(config["had_boostedtau_decay_mode"]):
#         if i == 0:
#             cut_string += "(boosted_tau_decaymode_2 == {})".format(dm)
#         else:
#             cut_string += "|| (boosted_tau_decaymode_2 == {})".format(dm)

#     if channel == "et":
#         rdf = rdf.Filter(
#             cut_string,
#             "cut had. boosted tau decay modes: {}".format(
#                 " ".join(dm for dm in config["had_boostedtau_decay_mode"])
#             ),
#         )
#     elif channel == "mt":
#         rdf = rdf.Filter(
#             cut_string,
#             "cut had. boosted tau decay modes: {}".format(
#                 " ".join(dm for dm in config["had_boostedtau_decay_mode"])
#             ),
#         )
#     elif channel == "tt":
#         rdf = rdf.Filter(
#             cut_string,
#             "cut second had. boosted tau decay modes: {}".format(
#                 " ".join(dm for dm in config["had_boostedtau_decay_mode"])
#             ),
#         )
#         rdf = rdf.Filter(
#             cut_string.replace("boosted_tau_decaymode_2", "boosted_tau_decaymode_1"),
#             "cut first had. boosted tau decay modes: {}".format(
#                 " ".join(dm for dm in config["had_boostedtau_decay_mode"])
#             ),
#         )
#     else:
#         sys.exit(
#             "Eventfilter: boosted tau decay mode: Such a channel is not defined: {}".format(
#                 channel
#             )
#         )
#     return rdf


# def had_tau_id_vsEle_cut(rdf, channel, config):
#     if channel == "et":
#         rdf = rdf.Filter(
#             "(id_tau_vsEle_{WP}_2 > 0.5)".format(WP=config["had_tau_id_vs_ele"]),
#             "cut on {WP} tau vs electrons id".format(WP=config["had_tau_id_vs_ele"]),
#         )
#     elif channel == "mt":
#         rdf = rdf.Filter(
#             "(id_tau_vsEle_{WP}_2 > 0.5)".format(WP=config["had_tau_id_vs_ele"]),
#             "cut on {WP} tau vs electrons id".format(WP=config["had_tau_id_vs_ele"]),
#         )
#     elif channel == "tt":
#         rdf = rdf.Filter(
#             "(id_tau_vsEle_{WP}_1 > 0.5) && (id_tau_vsEle_{WP}_2 > 0.5)".format(
#                 WP=config["had_tau_id_vs_ele"]
#             ),
#             "cut on {WP} tau vs electrons id".format(WP=config["had_tau_id_vs_ele"]),
#         )
#     else:
#         sys.exit(
#             "Eventfilter: tau id vs electron: Such a channel is not defined: {}".format(
#                 channel
#             )
#         )
#     return rdf


# def had_boostedtau_id_antiEle_cut(rdf, channel, config):
#     if channel == "et":
#         rdf = rdf.Filter(
#             "(id_boostedtau_antiEle_{WP}_2 > 0.5)".format(WP=config["had_boostedtau_id_antiEle"]),
#             "cut on {WP} tau old anti electrons id".format(WP=config["had_boostedtau_id_antiEle"]),
#         )
#     elif channel == "mt":
#         rdf = rdf.Filter(
#             "(id_boostedtau_antiEle_{WP}_2 > 0.5)".format(WP=config["had_boostedtau_id_antiEle"]),
#             "cut on {WP} tau old anti electrons id".format(WP=config["had_boostedtau_id_antiEle"]),
#         )
#     elif channel == "tt":
#         rdf = rdf.Filter(
#             "(id_boostedtau_antiEle_{WP}_1 > 0.5) && (id_boostedtau_antiEle_{WP}_2 > 0.5)".format(
#                 WP=config["had_boostedtau_id_antiEle"]
#             ),
#             "cut on {WP} tau old anti electrons id".format(WP=config["had_boostedtau_id_antiEle"]),
#         )
#     else:
#         sys.exit(
#             "Eventfilter: boosted tau id anti electron: Such a channel is not defined: {}".format(
#                 channel
#             )
#         )
#     return rdf


# def had_tau_id_vsMu_cut(rdf, channel, config):
#     if channel == "et":
#         rdf = rdf.Filter(
#             "(id_tau_vsMu_{WP}_2 > 0.5)".format(WP=config["had_tau_id_vs_mu"]),
#             "cut on {WP} tau vs muons id".format(WP=config["had_tau_id_vs_mu"]),
#         )
#     elif channel == "mt":
#         rdf = rdf.Filter(
#             "(id_tau_vsMu_{WP}_2 > 0.5)".format(WP=config["had_tau_id_vs_mu"]),
#             "cut on {WP} tau vs muons id".format(WP=config["had_tau_id_vs_mu"]),
#         )
#     elif channel == "tt":
#         rdf = rdf.Filter(
#             "(id_tau_vsMu_{WP}_1 > 0.5) && (id_tau_vsMu_{WP}_2 > 0.5)".format(
#                 WP=config["had_tau_id_vs_mu"],
#             ),
#             "cut on {WP} tau vs muons id".format(WP=config["had_tau_id_vs_mu"]),
#         )
#     else:
#         sys.exit(
#             "Eventfilter: tau id vs muon: Such a channel is not defined: {}".format(
#                 channel
#             )
#         )
#     return rdf


# def had_boostedtau_id_antiMu_cut(rdf, channel, config):
#     if channel == "et":
#         rdf = rdf.Filter(
#             "(id_boostedtau_antiMu_{WP}_2 > 0.5)".format(WP=config["had_boostedtau_id_antiMu"]),
#             "cut on {WP} tau old anti muons id".format(WP=config["had_boostedtau_id_antiMu"]),
#         )
#     elif channel == "mt":
#         rdf = rdf.Filter(
#             "(id_boostedtau_antiMu_{WP}_2 > 0.5)".format(WP=config["had_boostedtau_id_antiMu"]),
#             "cut on {WP} tau old anti muons id".format(WP=config["had_boostedtau_id_antiMu"]),
#         )
#     elif channel == "tt":
#         rdf = rdf.Filter(
#             "(id_boostedtau_antiMu_{WP}_1 > 0.5) && (id_boostedtau_antiMu_{WP}_2 > 0.5)".format(
#                 WP=config["had_boostedtau_id_antiMu"],
#             ),
#             "cut on {WP} tau old anti muons id".format(WP=config["had_boostedtau_id_antiMu"]),
#         )
#     else:
#         sys.exit(
#             "Eventfilter: boosted tau id anti muon: Such a channel is not defined: {}".format(
#                 channel
#             )
#         )
#     return rdf


def had_tau_id_vsJet_cut(rdf, channel, config):
    if channel in ["et", "mt"]:
        if isinstance(config["had_tau_id_vs_jet"], str):
            rdf = rdf.Filter(
                "(id_tau_vsJet_{WP}_2 > 0.5)".format(WP=config["had_tau_id_vs_jet"]),
                "cut on {WP} tau vs jets id".format(WP=config["had_tau_id_vs_jet"]),
            )
        elif isinstance(config["had_tau_id_vs_jet"], list):
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

    elif channel == "tt":
        if isinstance(config["had_tau_id_vs_jet_1"], str):
            rdf = rdf.Filter(
                "(id_tau_vsJet_{WP}_1 > 0.5)".format(WP=config["had_tau_id_vs_jet_1"]),
                "cut on {WP} tau vs jets id 1".format(WP=config["had_tau_id_vs_jet_1"]),
            )
        elif isinstance(config["had_tau_id_vs_jet_1"], list):
            rdf = rdf.Filter(
                "(id_tau_vsJet_{lower_WP}_1 > 0.5) && (id_tau_vsJet_{upper_WP}_1 < 0.5)".format(
                    lower_WP=config["had_tau_id_vs_jet_1"][0],
                    upper_WP=config["had_tau_id_vs_jet_1"][1],
                ),
                "cut on tau vs jets id 1 between {lower_WP} and {upper_WP}".format(
                    lower_WP=config["had_tau_id_vs_jet_1"][0],
                    upper_WP=config["had_tau_id_vs_jet_1"][1],
                ),
            )

        if isinstance(config["had_tau_id_vs_jet_2"], str):
            rdf = rdf.Filter(
                "(id_tau_vsJet_{WP}_2 > 0.5)".format(WP=config["had_tau_id_vs_jet_2"]),
                "cut on {WP} tau vs jets id 2".format(WP=config["had_tau_id_vs_jet_2"]),
            )
        elif isinstance(config["had_tau_id_vs_jet_2"], list):
            rdf = rdf.Filter(
                "(id_tau_vsJet_{lower_WP}_2 > 0.5) && (id_tau_vsJet_{upper_WP}_2 < 0.5)".format(
                    lower_WP=config["had_tau_id_vs_jet_2"][0],
                    upper_WP=config["had_tau_id_vs_jet_2"][1],
                ),
                "cut on tau vs jets id 2 between {lower_WP} and {upper_WP}".format(
                    lower_WP=config["had_tau_id_vs_jet_2"][0],
                    upper_WP=config["had_tau_id_vs_jet_2"][1],
                ),
            )
    else:
        sys.exit(
            "Eventfilter: tau id vs jet: Such a channel is not defined: {}".format(
                channel
            )
        )
    return rdf


def had_boostedtau_id_iso_cut(rdf, channel, config):
    if channel in ["et", "mt"]:
        if isinstance(config["had_boostedtau_id_iso"], str):
            rdf = rdf.Filter(
                "(id_boostedtau_iso_{WP}_2 > 0.5)".format(
                    WP=config["had_boostedtau_id_iso"]
                ),
                "cut on {WP} boosted tau iso id".format(
                    WP=config["had_boostedtau_id_iso"]
                ),
            )
        elif isinstance(config["had_boostedtau_id_iso"], list):
            if config["had_boostedtau_id_iso"][0] != "":
                rdf = rdf.Filter(
                    "(id_boostedtau_iso_{lower_WP}_2 > 0.5) && (id_boostedtau_iso_{upper_WP}_2 < 0.5)".format(
                        lower_WP=config["had_boostedtau_id_iso"][0],
                        upper_WP=config["had_boostedtau_id_iso"][1],
                    ),
                    "cut on boosted tau iso id between {lower_WP} and {upper_WP}".format(
                        lower_WP=config["had_boostedtau_id_iso"][0],
                        upper_WP=config["had_boostedtau_id_iso"][1],
                    ),
                )
            else:
                rdf = rdf.Filter(
                    "(id_boostedtau_iso_{upper_WP}_2 < 0.5)".format(
                        upper_WP=config["had_boostedtau_id_iso"][1],
                    ),
                    "cut on boosted tau iso id between {lower_WP} and {upper_WP}".format(
                        lower_WP=config["had_boostedtau_id_iso"][0],
                        upper_WP=config["had_boostedtau_id_iso"][1],
                    ),
                )

    elif channel == "tt":
        if isinstance(config["had_boostedtau_id_iso_1"], str):
            rdf = rdf.Filter(
                "(id_boostedtau_iso_{WP}_1 > 0.5)".format(
                    WP=config["had_boostedtau_id_iso_1"]
                ),
                "cut on {WP} boosted tau iso id 1".format(
                    WP=config["had_boostedtau_id_iso_1"]
                ),
            )
        elif isinstance(config["had_boostedtau_id_iso_1"], list):
            if config["had_boostedtau_id_iso"][0] != "":
                rdf = rdf.Filter(
                    "(id_boostedtau_iso_{lower_WP}_1 > 0.5) && (id_boostedtau_iso_{upper_WP}_1 < 0.5)".format(
                        lower_WP=config["had_boostedtau_id_iso_1"][0],
                        upper_WP=config["had_boostedtau_id_iso_1"][1],
                    ),
                    "cut on boosted tau iso id 1 between {lower_WP} and {upper_WP}".format(
                        lower_WP=config["had_boostedtau_id_iso_1"][0],
                        upper_WP=config["had_boostedtau_id_iso_1"][1],
                    ),
                )
            else:
                rdf = rdf.Filter(
                    "(id_boostedtau_iso_{upper_WP}_1 < 0.5)".format(
                        upper_WP=config["had_boostedtau_id_iso_1"][1],
                    ),
                    "cut on boosted tau iso id 1 between {lower_WP} and {upper_WP}".format(
                        lower_WP=config["had_boostedtau_id_iso_1"][0],
                        upper_WP=config["had_boostedtau_id_iso_1"][1],
                    ),
                )

        if isinstance(config["had_boostedtau_id_iso_2"], str):
            rdf = rdf.Filter(
                "(id_boostedtau_iso_{WP}_2 > 0.5)".format(
                    WP=config["had_boostedtau_id_iso_2"]
                ),
                "cut on {WP} boosted tau iso id 2".format(
                    WP=config["had_boostedtau_id_iso_2"]
                ),
            )
        elif isinstance(config["had_boostedtau_id_iso_2"], list):
            if config["had_boostedtau_id_iso"][0] != "":
                rdf = rdf.Filter(
                    "(id_boostedtau_iso_{lower_WP}_2 > 0.5) && (id_boostedtau_iso_{upper_WP}_2 < 0.5)".format(
                        lower_WP=config["had_boostedtau_id_iso_2"][0],
                        upper_WP=config["had_boostedtau_id_iso_2"][1],
                    ),
                    "cut on boosted tau iso id 2 between {lower_WP} and {upper_WP}".format(
                        lower_WP=config["had_boostedtau_id_iso_2"][0],
                        upper_WP=config["had_boostedtau_id_iso_2"][1],
                    ),
                )
            else:
                rdf = rdf.Filter(
                    "(id_boostedtau_iso_{upper_WP}_2 < 0.5)".format(
                        upper_WP=config["had_boostedtau_id_iso_2"][1],
                    ),
                    "cut on boosted tau iso id 2 between {lower_WP} and {upper_WP}".format(
                        lower_WP=config["had_boostedtau_id_iso_2"][0],
                        upper_WP=config["had_boostedtau_id_iso_2"][1],
                    ),
                )
    else:
        sys.exit(
            "Eventfilter: boosted tau id iso: Such a channel is not defined: {}".format(
                channel
            )
        )
    return rdf


def lep_iso_cut(rdf, channel, config):
    if channel == "et":
        # check if lep_iso is a list
        if isinstance(config["lep_iso"], list):
            rdf = rdf.Filter(
                "(iso_1 {}) && (iso_1 {})".format(
                    config["lep_iso"][0], config["lep_iso"][1]
                ),
                "cut muon isolation {} and {}".format(
                    config["lep_iso"][0], config["lep_iso"][1]
                ),
            )
        else:
            rdf = rdf.Filter(
                "iso_1 {}".format(config["lep_iso"]),
                "cut muon isolation {}".format(config["lep_iso"]),
            )
    elif channel == "mt":
        # check if lep_iso is a list
        if isinstance(config["lep_iso"], list):
            rdf = rdf.Filter(
                "(iso_1 {}) && (iso_1 {})".format(
                    config["lep_iso"][0], config["lep_iso"][1]
                ),
                "cut muon isolation {} and {}".format(
                    config["lep_iso"][0], config["lep_iso"][1]
                ),
            )
        else:
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


def boosted_lep_iso_cut(rdf, channel, config):
    if channel == "et":
        # check if lep_iso is a list
        if isinstance(config["boosted_lep_iso"], list):
            rdf = rdf.Filter(
                "(boosted_iso_1 {}) && (boosted_iso_1 {})".format(
                    config["boosted_lep_iso"][0], config["boosted_lep_iso"][1]
                ),
                "cut boosted electron isolation {} and {}".format(
                    config["boosted_lep_iso"][0], config["boosted_lep_iso"][1]
                ),
            )
        else:
            rdf = rdf.Filter(
                "boosted_iso_1 {}".format(config["boosted_lep_iso"]),
                "cut boosted electron isolation {}".format(config["boosted_lep_iso"]),
            )
    elif channel == "mt":
        # check if lep_iso is a list
        if isinstance(config["boosted_lep_iso"], list):
            rdf = rdf.Filter(
                "(boosted_iso_1 {}) && (boosted_iso_1 {})".format(
                    config["boosted_lep_iso"][0], config["boosted_lep_iso"][1]
                ),
                "cut boosted muon isolation {} and {}".format(
                    config["boosted_lep_iso"][0], config["boosted_lep_iso"][1]
                ),
            )
        else:
            rdf = rdf.Filter(
                "boosted_iso_1 {}".format(config["boosted_lep_iso"]),
                "cut boosted muon isolation {}".format(config["boosted_lep_iso"]),
            )
    elif channel == "tt":
        pass
    else:
        sys.exit(
            "Eventfilter: boosted lepton isolation: Such a channel is not defined: {}".format(
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
        pass
    else:
        sys.exit(
            "Eventfilter: transverse mass of lepton + MET: Such a channel is not defined: {}".format(
                channel
            )
        )
    return rdf


def boosted_lep_mt_cut(rdf, channel, config):
    if channel == "et":
        rdf = rdf.Filter(
            "(boosted_mt_1 {})".format(config["boosted_lep_mt"]),
            "W boson origin cut on boosted lepton mt {}".format(
                config["boosted_lep_mt"]
            ),
        )
    elif channel == "mt":
        rdf = rdf.Filter(
            "(boosted_mt_1 {})".format(config["boosted_lep_mt"]),
            "W boson origin cut on boosted lepton mT {}".format(
                config["boosted_lep_mt"]
            ),
        )
    elif channel == "tt":
        pass
    else:
        sys.exit(
            "Eventfilter: boosted transverse mass of lepton + MET: Such a channel is not defined: {}".format(
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


def boostedtau_pair_sign_cut(rdf, channel, sign):
    if sign == "same":
        sign_cut = "> 0"
    elif sign == "opposite":
        sign_cut = "< 0"
    else:
        sys.exit(
            "Eventfilter: boosted SS/OS: Such a sign is not defined: {}".format(sign)
        )

    if channel == "et":
        rdf = rdf.Filter(
            "((boosted_q_1 * boosted_q_2) {})".format(sign_cut),
            "boosted {} sign cut".format(sign),
        )
    elif channel == "mt":
        rdf = rdf.Filter(
            "((boosted_q_1 * boosted_q_2) {})".format(sign_cut),
            "boosted {} sign cut".format(sign),
        )
    elif channel == "tt":
        rdf = rdf.Filter(
            "((boosted_q_1 * boosted_q_2) {})".format(sign_cut),
            "boosted {} sign cut".format(sign),
        )
    else:
        sys.exit(
            "Eventfilter: boosted SS/OS: Such a channel is not defined: {}".format(
                channel
            )
        )
    return rdf


def tau_pair_dR_cut(rdf, channel, dR):
    if channel == "et":
        rdf = rdf.Filter(
            "(deltaR_ditaupair {})".format(dR), "deltaR(tau pair) {} cut".format(dR)
        )
    elif channel == "mt":
        rdf = rdf.Filter(
            "(deltaR_ditaupair {})".format(dR), "deltaR(tau pair) {} cut".format(dR)
        )
    elif channel == "tt":
        rdf = rdf.Filter(
            "(deltaR_ditaupair {})".format(dR), "deltaR(tau pair) {} cut".format(dR)
        )
    else:
        sys.exit(
            "Eventfilter: deltaR(tautau): Such a channel is not defined: {}".format(
                channel
            )
        )
    return rdf


def boostedtau_pair_dR_cut(rdf, channel, dR):
    if channel == "et":
        rdf = rdf.Filter(
            "(boosted_deltaR_ditaupair {})".format(dR),
            "boosted deltaR(tau pair) {} cut".format(dR),
        )
    elif channel == "mt":
        rdf = rdf.Filter(
            "(boosted_deltaR_ditaupair {})".format(dR),
            "boosted deltaR(tau pair) {} cut".format(dR),
        )
    elif channel == "tt":
        rdf = rdf.Filter(
            "(boosted_deltaR_ditaupair {})".format(dR),
            "boosted deltaR(tau pair) {} cut".format(dR),
        )
    else:
        sys.exit(
            "Eventfilter: boosted deltaR(tautau): Such a channel is not defined: {}".format(
                channel
            )
        )
    return rdf


# def trigger_cut(rdf, channel):
#     if channel == "et":
#         rdf = rdf.Filter(
#             "((pt_1 > 33) && ((trg_single_ele32 > 0.5) || (trg_single_ele35 > 0.5)))",  # || ((pt_2 > 35) && (abs(eta_2) < 2.1) && (pt_1 <= 33) && (pt_1 > 25) && ((trg_cross_ele24tau30 > 0.5) || (trg_cross_ele24tau30_hps > 0.5)))",
#             "single electron trigger + electron pT cuts",
#         )
#     elif channel == "mt":
#         rdf = rdf.Filter(
#             "((pt_1 > 25) && ((trg_single_mu24 > 0.5) || (trg_single_mu27 > 0.5)))",  # || ((pt_2 > 32) && (pt_1 <= 25) && (pt_1 > 21) && ((trg_cross_mu20tau27 > 0.5) || (trg_cross_mu20tau27_hps > 0.5)))",
#             "single muon trigger + muon pT cuts",
#         )
#     elif channel == "tt":
#         rdf = rdf.Filter(
#             "(trg_double_tau35_tightiso_tightid > 0.5) || (trg_double_tau35_mediumiso_hps > 0.5) || (trg_double_tau40_mediumiso_tightid > 0.5) || (trg_double_tau40_tightiso > 0.5)",
#             "tau double trigger cuts",
#         )
#     else:
#         sys.exit(
#             "Eventfilter: trigger: Such a channel is not defined: {}".format(channel)
#         )
#     return rdf

# def boosted_trigger_cut(rdf, channel):
#     if channel == "et":
#         rdf = rdf.Filter(
#             "((boosted_pt_1 > 120) && ((trg_single_ele115_boosted > 0.5)))",
#             "single electron trigger (no iso) + boosted electron pT cuts",
#         )
#     elif channel == "mt":
#         rdf = rdf.Filter(
#             "((boosted_pt_1 > 55) && (trg_single_mu50_boosted > 0.5))",
#             "single muon trigger (no iso) + boosted muon pT cuts",
#         )
#     elif channel == "tt":
#         rdf = rdf.Filter(
#             "(trg_double_tau35_tightiso_tightid > 0.5) || (trg_double_tau35_mediumiso_hps > 0.5) || (trg_double_tau40_mediumiso_tightid > 0.5) || (trg_double_tau40_tightiso > 0.5)",
#             "tau double trigger cuts",
#         )
#     else:
#         sys.exit(
#             "Eventfilter: boosted trigger: Such a channel is not defined: {}".format(channel)
#         )
#     return rdf


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
            "Eventfilter: jet number: Such a channel is not defined: {}".format(channel)
        )
    return rdf



