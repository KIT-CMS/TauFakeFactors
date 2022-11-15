import sys


def same_opposite_sign(rdf, channel, sign):
    if sign == "same":
        sign_cut = "> 0"
    elif sign == "opposite":
        sign_cut = "< 0"
    else:
        sys.exit("Regionfilter: SS/OS: Such a sign is not defined: {}".format(sign))

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
            "Regionfilter: SS/OS: Such a channel is not defined: {}".format(channel)
        )
    return rdf


def btag_number(rdf, channel, config):
    if channel == "et":
        rdf = rdf.Filter(
            "(nbtag {})".format(config["nbtag"]),
            "cut on {} b-tagged jets".format(config["nbtag"]),
        )
    elif channel == "mt":
        rdf = rdf.Filter(
            "(nbtag {})".format(config["nbtag"]),
            "cut on {} b-tagged jets".format(config["nbtag"]),
        )
    elif channel == "tt":
        rdf = rdf.Filter(
            "(nbtag {})".format(config["nbtag"]),
            "cut on {} b-tagged jets".format(config["nbtag"]),
        )
    else:
        sys.exit(
            "Regionfilter: btag number: Such a channel is not defined: {}".format(
                channel
            )
        )
    return rdf


def jet_number(rdf, channel, njets):
    if channel == "et":
        rdf = rdf.Filter("njets {}".format(njets), "cut on {} jets".format(njets))
    elif channel == "mt":
        rdf = rdf.Filter("njets {}".format(njets), "cut on {} jets".format(njets))
    elif channel == "tt":
        rdf = rdf.Filter("njets {}".format(njets), "cut on {} jets".format(njets))
    else:
        sys.exit(
            "Regionfilter: jet number: Such a channel is not defined: {}".format(
                channel
            )
        )
    return rdf


def no_extra_leptons(rdf, channel, config):
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
            "Regionfilter: no extra lep: Such a channel is not defined: {}".format(
                channel
            )
        )
    return rdf


def tau_id_vs_jets_WP(rdf, channel, config):
    if channel == "et":
        rdf = rdf.Filter(
            "(id_tau_vsJet_{}_2 > 0.5)".format(config["tau_id_vsJet"]),
            "cut on {} tau vs jets id".format(config["tau_id_vsJet"]),
        )
    elif channel == "mt":
        rdf = rdf.Filter(
            "(id_tau_vsJet_{}_2 > 0.5)".format(config["tau_id_vsJet"]),
            "cut on {} tau vs jets id".format(config["tau_id_vsJet"]),
        )
    elif channel == "tt":
        rdf = rdf.Filter(
            "(id_tau_vsJet_{}_2 > 0.5)".format(config["tau_id_vsJet"]),
            "cut on {} tau vs jets id".format(config["tau_id_vsJet"]),
        )
    else:
        sys.exit(
            "Regionfilter: tau vs jet ID: Such a channel is not defined: {}".format(
                channel
            )
        )
    return rdf


def tau_id_vs_jets_between_WPs(rdf, channel, config):
    if channel == "et":
        rdf = rdf.Filter(
            "(id_tau_vsJet_{}_2 > 0.5) && (id_tau_vsJet_{}_2 < 0.5)".format(
                config["tau_id_vsJet"][0], config["tau_id_vsJet"][1]
            ),
            "cut on tau vs jets id between {} and {}".format(
                config["tau_id_vsJet"][0], config["tau_id_vsJet"][1]
            ),
        )
    elif channel == "mt":
        rdf = rdf.Filter(
            "(id_tau_vsJet_{}_2 > 0.5) && (id_tau_vsJet_{}_2 < 0.5)".format(
                config["tau_id_vsJet"][0], config["tau_id_vsJet"][1]
            ),
            "cut on tau vs jets id between {} and {}".format(
                config["tau_id_vsJet"][0], config["tau_id_vsJet"][1]
            ),
        )
    elif channel == "tt":
        rdf = rdf.Filter(
            "(id_tau_vsJet_{}_2 > 0.5) && (id_tau_vsJet_{}_2 < 0.5)".format(
                config["tau_id_vsJet"][0], config["tau_id_vsJet"][1]
            ),
            "cut on tau vs jets id between {} and {}".format(
                config["tau_id_vsJet"][0], config["tau_id_vsJet"][1]
            ),
        )
    else:
        sys.exit(
            "Regionfilter: tau vs jet ID between WPs: Such a channel is not defined: {}".format(
                channel
            )
        )
    return rdf


def lepton_mT(rdf, channel, config):
    if channel == "et":
        rdf = rdf.Filter(
            "(mt_1 {})".format(config["lep_mT"]),
            "W boson origin cut on lepton mT {}".format(config["lep_mT"]),
        )
    elif channel == "mt":
        rdf = rdf.Filter(
            "(mt_1 {})".format(config["lep_mT"]),
            "W boson origin cut on lepton mT {}".format(config["lep_mT"]),
        )
    elif channel == "tt":
        rdf = rdf.Filter(
            "(mt_1 {})".format(config["lep_mT"]),
            "W boson origin cut on lepton mT {}".format(config["lep_mT"]),
        )
    else:
        sys.exit(
            "Regionfilter: lepton transverse mass: Such a channel is not defined: {}".format(
                channel
            )
        )
    return rdf
