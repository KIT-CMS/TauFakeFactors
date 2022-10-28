import sys

def signal_region_QCD(rdf, channel):
    if channel == "et":
        rdf = rdf.Filter("((q_1 * q_2) < 0)", "opposite sign cut")
        rdf = rdf.Filter("(id_tau_vsJet_Medium_2 > 0.5)", "tau vs jets id cuts")
    elif channel == "mt":
        rdf = rdf.Filter("((q_1 * q_2) < 0)", "opposite sign cut")
        rdf = rdf.Filter("(id_tau_vsJet_Medium_2 > 0.5)", "tau vs jets id cuts")
    elif channel == "tt":
        rdf = rdf.Filter("((q_1 * q_2) < 0)", "opposite sign cut")
        rdf = rdf.Filter(
            "(id_tau_vsJet_Medium_1 > 0.5) && (id_tau_vsJet_Medium_2 > 0.5)",
            "tau vs jets id cuts",
        )
    else:
        sys.exit("Regionfilter: signal: QCD: Such a channel is not defined: {}".format(channel))
    return rdf

def application_region_QCD(rdf, channel):
    if channel == "et":
        rdf = rdf.Filter("((q_1 * q_2) < 0)", "opposite sign cut")
        rdf = rdf.Filter(
            "(id_tau_vsJet_VVVLoose_2 > 0.5) && (id_tau_vsJet_Medium_2 < 0.5)",
            "tau vs jets id cuts",
        )
    elif channel == "mt":
        rdf = rdf.Filter("((q_1 * q_2) < 0)", "opposite sign cut")
        rdf = rdf.Filter(
            "(id_tau_vsJet_VVVLoose_2 > 0.5) && (id_tau_vsJet_Medium_2 < 0.5)",
            "tau vs jets id cuts",
        )
    elif channel == "tt":
        rdf = rdf.Filter("((q_1 * q_2) < 0)", "opposite sign cut")
        rdf = rdf.Filter(
            "(id_tau_vsJet_VVVLoose_1 > 0.5) && (id_tau_vsJet_Medium_1 < 0.5) && (id_tau_vsJet_VVVLoose_2 > 0.5) && (id_tau_vsJet_Medium_2 < 0.5)",
            "tau vs jets id cuts",
        )
    else:
        sys.exit("Regionfilter:  application: QCD: Such a channel is not defined: {}".format(channel))
    return rdf

def signal_like_region_QCD(rdf, channel, njets):
    if channel == "et":
        rdf = rdf.Filter("((q_1 * q_2) > 0)", "same sign cut")
        rdf = rdf.Filter("(mt_1 < 50)", "W boson origin exclusion cut")
        rdf = rdf.Filter("(id_tau_vsJet_Medium_2 > 0.5)", "tau vs jets id cuts")
        rdf = rdf.Filter("njets {}".format(njets), "cut on {} jets".format(njets))
    elif channel == "mt":
        rdf = rdf.Filter("((q_1 * q_2) > 0)", "same sign cut")
        rdf = rdf.Filter("(mt_1 < 50)", "W boson origin exclusion cut")
        rdf = rdf.Filter("(id_tau_vsJet_Medium_2 > 0.5)", "tau vs jets id cuts")
        rdf = rdf.Filter("njets == {}".format(njets), "cut on {} jets".format(njets))
    elif channel == "tt":
        rdf = rdf.Filter("((q_1 * q_2) > 0)", "same sign cut")
        rdf = rdf.Filter(
            "(id_tau_vsJet_Medium_1 > 0.5) && (id_tau_vsJet_Medium_2 > 0.5)",
            "tau vs jets id cuts",
        )
        rdf = rdf.Filter("njets == {}".format(njets), "cut on {} jets".format(njets))
    else:
        sys.exit("Regionfilter: signal-like: QCD: Such a channel is not defined: {}".format(channel))
    return rdf

def application_like_region_QCD(rdf, channel, njets):
    if channel == "et":
        rdf = rdf.Filter("((q_1 * q_2) > 0)", "same sign cut")
        rdf = rdf.Filter("(mt_1 < 50)", "W boson origin exclusion cut")
        rdf = rdf.Filter(
            "(id_tau_vsJet_VVVLoose_2 > 0.5) && (id_tau_vsJet_Medium_2 < 0.5)",
            "tau vs jets id cuts",
        )
        rdf = rdf.Filter("njets {}".format(njets), "cut on {} jets".format(njets))
    elif channel == "mt":
        rdf = rdf.Filter("((q_1 * q_2) > 0)", "same sign cut")
        rdf = rdf.Filter("(mt_1 < 50)", "W boson origin exclusion cut")
        rdf = rdf.Filter(
            "(id_tau_vsJet_VVVLoose_2 > 0.5) && (id_tau_vsJet_Medium_2 < 0.5)",
            "tau vs jets id cuts",
        )
        rdf = rdf.Filter("njets == {}".format(njets), "cut on {} jets".format(njets))
    elif channel == "tt":
        rdf = rdf.Filter("((q_1 * q_2) > 0)", "same sign cut")
        rdf = rdf.Filter(
            "(id_tau_vsJet_VVVLoose_1 > 0.5) && (id_tau_vsJet_Medium_1 < 0.5) && (id_tau_vsJet_VVVLoose_2 > 0.5) && (id_tau_vsJet_Medium_2 < 0.5)",
            "tau vs jets id cuts",
        )
        rdf = rdf.Filter("njets == {}".format(njets), "cut on {} jets".format(njets))
    else:
        sys.exit("Regionfilter: application-like: QCD: Such a channel is not defined: {}".format(channel))
    return rdf


def signal_like_region_Wjets(rdf, channel, njets):
    if channel == "et":
        rdf = rdf.Filter("(nbtag == 0)", "no b-tagged jets cut")
        rdf = rdf.Filter("(mt_1 > 70)", "W boson origin cut")
        rdf = rdf.Filter("(id_tau_vsJet_Medium_2 > 0.5)", "tau vs jets id cuts")
        rdf = rdf.Filter("njets {}".format(njets), "cut on {} jets".format(njets))
    elif channel == "mt":
        rdf = rdf.Filter("(nbtag == 0)", "no b-tagged jets cut")
        rdf = rdf.Filter("(mt_1 > 70)", "W boson origin cut")
        rdf = rdf.Filter("(id_tau_vsJet_Medium_2 > 0.5)", "tau vs jets id cuts")
        rdf = rdf.Filter("njets {}".format(njets), "cut on {} jets".format(njets))
    elif channel == "tt":
        rdf = rdf.Filter("(nbtag == 0)", "no b-tagged jets cut")
        rdf = rdf.Filter(
            "(id_tau_vsJet_Medium_1 > 0.5) && (id_tau_vsJet_Medium_2 > 0.5)",
            "tau vs jets id cuts",
        )
        rdf = rdf.Filter("njets {}".format(njets), "cut on {} jets".format(njets))
    else:
        sys.exit("Regionfilter: signal-like: Wjets: Such a channel is not defined: {}".format(channel))
    return rdf

def application_like_region_Wjets(rdf, channel, njets):
    if channel == "et":
        rdf = rdf.Filter("(nbtag == 0)", "no b-tagged jets cut")
        rdf = rdf.Filter("(mt_1 > 70)", "W boson origin cut")
        rdf = rdf.Filter(
            "(id_tau_vsJet_VVVLoose_2 > 0.5) && (id_tau_vsJet_Medium_2 < 0.5)",
            "tau vs jets id cuts",
        )
        rdf = rdf.Filter("njets {}".format(njets), "cut on {} jets".format(njets))
    elif channel == "mt":
        rdf = rdf.Filter("(nbtag == 0)", "no b-tagged jets cut")
        rdf = rdf.Filter("(mt_1 > 70)", "W boson origin cut")
        rdf = rdf.Filter(
            "(id_tau_vsJet_VVVLoose_2 > 0.5) && (id_tau_vsJet_Medium_2 < 0.5)",
            "tau vs jets id cuts",
        )
        rdf = rdf.Filter("njets {}".format(njets), "cut on {} jets".format(njets))
    elif channel == "tt":
        rdf = rdf.Filter("(nbtag == 0)", "no b-tagged jets cut")
        rdf = rdf.Filter(
            "(id_tau_vsJet_VVVLoose_1 > 0.5) && (id_tau_vsJet_Medium_1 < 0.5) && (id_tau_vsJet_VVVLoose_2 > 0.5) && (id_tau_vsJet_Medium_2 < 0.5)",
            "tau vs jets id cuts",
        )
        rdf = rdf.Filter("njets {}".format(njets), "cut on {} jets".format(njets))
    else:
        sys.exit("Regionfilter: application-like: Wjets: Such a channel is not defined: {}".format(channel))
    return rdf

def same_sign(rdf, channel):
    if channel == "et":
        rdf = rdf.Filter("((q_1 * q_2) > 0)", "same sign cut")
    elif channel == "mt":
        rdf = rdf.Filter("((q_1 * q_2) > 0)", "same sign cut")
    elif channel == "tt":
        rdf = rdf.Filter("((q_1 * q_2) > 0)", "same sign cut")
    else:
        sys.exit("Regionfilter: QCD SS: Such a channel is not defined: {}".format(channel))
    return rdf

def opposite_sign(rdf, channel):
    if channel == "et":
        rdf = rdf.Filter("((q_1 * q_2) < 0)", "opposite sign cut")
    elif channel == "mt":
        rdf = rdf.Filter("((q_1 * q_2) < 0)", "opposite sign cut")
    elif channel == "tt":
        rdf = rdf.Filter("((q_1 * q_2) < 0)", "opposite sign cut")
    else:
        sys.exit("Regionfilter: QCD SS: Such a channel is not defined: {}".format(channel))
    return rdf