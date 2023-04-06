import sys
import array


def lep_iso_weight(rdf, channel):
    if channel == "et":
        rdf = rdf.Redefine("weight", "weight * iso_wgt_ele_1")
    elif channel == "mt":
        rdf = rdf.Redefine("weight", "weight * iso_wgt_mu_1")
    elif channel == "tt":
        rdf = rdf.Redefine("weight", "weight")
    else:
        sys.exit("Weight calc: iso: Such a channel is not defined: {}".format(channel))

    return rdf


def lep_id_weight(rdf, channel):
    if channel == "et":
        rdf = rdf.Redefine("weight", "weight * id_wgt_ele_1")
    elif channel == "mt":
        rdf = rdf.Redefine("weight", "weight * id_wgt_mu_1")
    elif channel == "tt":
        rdf = rdf.Redefine("weight", "weight")
    else:
        sys.exit("Weight calc: id: Such a channel is not defined: {}".format(channel))

    return rdf


def apply_btag_weight(rdf):
    # Calculating correction ratio for b-tagging SFs https://twiki.cern.ch/twiki/bin/view/CMS/BTagShapeCalibration#Effect_on_event_yields
    rdf = rdf.Define("wgt_with_btag", "weight * btag_weight")

    # measure corr. ratio for N jets (0 to 8); the highest N jet here is an arbitrary choice
    xbinning = array.array("d", [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5])
    nbinsx = len(xbinning) - 1

    # histogram without b-tagging SFs
    h = rdf.Histo1D(("", "", nbinsx, xbinning), "njets", "weight")
    h = h.GetValue()
    # histogram with b-tagging SFs
    h_btag = rdf.Histo1D(("", "", nbinsx, xbinning), "njets", "wgt_with_btag")
    h_btag = h_btag.GetValue()

    ratio = h.Clone()
    ratio.Divide(h_btag)

    # generating expression for b-tagging weights with the N jets dependent corr. ratio
    btag_wgt = "btag_weight*("
    for n in range(nbinsx):
        if ratio.GetBinContent(n + 1) != 0.0:
            btag_wgt += "(njets=={})*{}+".format(n, ratio.GetBinContent(n + 1))
        else:
            btag_wgt += "(njets=={})*{}+".format(n, 1.0)
    btag_wgt += "(njets>8)*1.)"

    # applying the b-tagging SFs
    rdf = rdf.Redefine("weight", "weight*{}".format(btag_wgt))
    return rdf


def apply_tau_id_vsJet_weight(rdf, channel, config):
    if isinstance(config["had_tau_id_vs_jet"], str):
        if channel == "et":
            rdf = rdf.Redefine(
                "weight",
                "weight * ((gen_match_2==5) * ((id_tau_vsJet_{WP}_2>0.5)*id_wgt_tau_vsJet_{WP}_2 + (id_tau_vsJet_{WP}_2<0.5)) + (gen_match_2!=5))".format(
                    WP=config["had_tau_id_vs_jet"]
                ),
            )
        elif channel == "mt":
            rdf = rdf.Redefine(
                "weight",
                "weight * ((gen_match_2==5) * ((id_tau_vsJet_{WP}_2>0.5)*id_wgt_tau_vsJet_{WP}_2 + (id_tau_vsJet_{WP}_2<0.5)) + (gen_match_2!=5))".format(
                    WP=config["had_tau_id_vs_jet"]
                ),
            )
        elif channel == "tt":
            rdf = rdf.Redefine(
                "weight",
                "weight * ((gen_match_1==5) * ((id_tau_vsJet_{WP}_1>0.5)*id_wgt_tau_vsJet_{WP}_1 + (id_tau_vsJet_{WP}_1<0.5)) + (gen_match_1!=5))".format(
                    WP=config["had_tau_id_vs_jet"]
                ),
            )
            rdf = rdf.Redefine(
                "weight",
                "weight * ((gen_match_2==5) * ((id_tau_vsJet_{WP}_2>0.5)*id_wgt_tau_vsJet_{WP}_2 + (id_tau_vsJet_{WP}_2<0.5)) + (gen_match_2!=5))".format(
                    WP=config["had_tau_id_vs_jet"]
                ),
            )
        else:
            sys.exit(
                "Weight calc: tau id vs jet: Such a channel is not defined: {}".format(
                    channel
                )
            )
    elif isinstance(config["had_tau_id_vs_jet"], list):
        if channel == "et":
            rdf = rdf.Redefine(
                "weight",
                "weight * ((gen_match_2==5) * ((id_tau_vsJet_{upper_WP}_2>0.5) + (id_tau_vsJet_{upper_WP}_2<0.5)*(id_tau_vsJet_{lower_WP}_2>0.5)*id_wgt_tau_vsJet_{lower_WP}_2 + (id_tau_vsJet_{lower_WP}_2<0.5)) + (gen_match_2!=5))".format(
                    lower_WP=config["had_tau_id_vs_jet"][0],
                    upper_WP=config["had_tau_id_vs_jet"][1],
                ),
            )
        elif channel == "mt":
            rdf = rdf.Redefine(
                "weight",
                "weight * ((gen_match_2==5) * ((id_tau_vsJet_{upper_WP}_2>0.5) + (id_tau_vsJet_{upper_WP}_2<0.5)*(id_tau_vsJet_{lower_WP}_2>0.5)*id_wgt_tau_vsJet_{lower_WP}_2 + (id_tau_vsJet_{lower_WP}_2<0.5)) + (gen_match_2!=5))".format(
                    lower_WP=config["had_tau_id_vs_jet"][0],
                    upper_WP=config["had_tau_id_vs_jet"][1],
                ),
            )
        elif channel == "tt":
            rdf = rdf.Redefine(
                "weight",
                "weight * ((gen_match_1==5) * ((id_tau_vsJet_{upper_WP}_1>0.5) + (id_tau_vsJet_{upper_WP}_1<0.5)*(id_tau_vsJet_{lower_WP}_1>0.5)*id_wgt_tau_vsJet_{lower_WP}_1 + (id_tau_vsJet_{lower_WP}_1<0.5)) + (gen_match_1!=5))".format(
                    lower_WP=config["had_tau_id_vs_jet"][0],
                    upper_WP=config["had_tau_id_vs_jet"][1],
                ),
            )
            rdf = rdf.Redefine(
                "weight",
                "weight * ((gen_match_2==5) * ((id_tau_vsJet_{upper_WP}_2>0.5) + (id_tau_vsJet_{upper_WP}_2<0.5)*(id_tau_vsJet_{lower_WP}_2>0.5)*id_wgt_tau_vsJet_{lower_WP}_2 + (id_tau_vsJet_{lower_WP}_2<0.5)) + (gen_match_2!=5))".format(
                    lower_WP=config["had_tau_id_vs_jet"][0],
                    upper_WP=config["had_tau_id_vs_jet"][1],
                ),
            )
        else:
            sys.exit(
                "Weight calc: tau id vs jet: Such a channel is not defined: {}".format(
                    channel
                )
            )
    else:
        sys.exit(
            "Weight calc: tau id vs jet: Such a type is not defined: {}".format(
                config["had_tau_id_vs_jet"]
            )
        )

    return rdf


def tau_id_vsMu_weight(rdf, channel, config):
    if channel == "et":
        rdf = rdf.Redefine(
            "weight",
            "weight * ((gen_match_2==5) * ((id_tau_vsMu_{WP}_2>0.5)*id_wgt_tau_vsMu_{WP}_2 + (id_tau_vsMu_{WP}_2<0.5)) + (gen_match_2!=5))".format(
                WP=config["had_tau_id_vs_mu"]
            ),
        )
    elif channel == "mt":
        rdf = rdf.Redefine(
            "weight",
            "weight * ((gen_match_2==5) * ((id_tau_vsMu_{WP}_2>0.5)*id_wgt_tau_vsMu_{WP}_2 + (id_tau_vsMu_{WP}_2<0.5)) + (gen_match_2!=5))".format(
                WP=config["had_tau_id_vs_mu"]
            ),
        )
    elif channel == "tt":
        rdf = rdf.Redefine(
            "weight",
            "weight * ((gen_match_1==5) * ((id_tau_vsMu_{WP}_1>0.5)*id_wgt_tau_vsMu_{WP}_1 + (id_tau_vsMu_{WP}_1<0.5)) + (gen_match_1!=5))".format(
                WP=config["had_tau_id_vs_mu"]
            ),
        )
        rdf = rdf.Redefine(
            "weight",
            "weight * ((gen_match_2==5) * ((id_tau_vsMu_{WP}_2>0.5)*id_wgt_tau_vsMu_{WP}_2 + (id_tau_vsMu_{WP}_2<0.5)) + (gen_match_2!=5))".format(
                WP=config["had_tau_id_vs_mu"]
            ),
        )
    else:
        sys.exit(
            "Weight calc: tau id vs muon: Such a channel is not defined: {}".format(
                channel
            )
        )

    return rdf


def tau_id_vsEle_weight(rdf, channel, config):
    if channel == "et":
        rdf = rdf.Redefine(
            "weight",
            "weight * ((gen_match_2==5) * ((id_tau_vsEle_{WP}_2>0.5)*id_wgt_tau_vsEle_{WP}_2 + (id_tau_vsEle_{WP}_2<0.5)) + (gen_match_2!=5))".format(
                WP=config["had_tau_id_vs_ele"]
            ),
        )
    elif channel == "mt":
        rdf = rdf.Redefine(
            "weight",
            "weight * ((gen_match_2==5) * ((id_tau_vsEle_{WP}_2>0.5)*id_wgt_tau_vsEle_{WP}_2 + (id_tau_vsEle_{WP}_2<0.5)) + (gen_match_2!=5))".format(
                WP=config["had_tau_id_vs_ele"]
            ),
        )
    elif channel == "tt":
        rdf = rdf.Redefine(
            "weight",
            "weight * ((gen_match_1==5) * ((id_tau_vsEle_{WP}_1>0.5)*id_wgt_tau_vsEle_{WP}_1 + (id_tau_vsEle_{WP}_1<0.5)) + (gen_match_1!=5))".format(
                WP=config["had_tau_id_vs_ele"]
            ),
        )
        rdf = rdf.Redefine(
            "weight",
            "weight * ((gen_match_2==5) * ((id_tau_vsEle_{WP}_2>0.5)*id_wgt_tau_vsEle_{WP}_2 + (id_tau_vsEle_{WP}_2<0.5)) + (gen_match_2!=5))".format(
                WP=config["had_tau_id_vs_ele"]
            ),
        )
    else:
        sys.exit(
            "Weight calc: tau id vs electron: Such a channel is not defined: {}".format(
                channel
            )
        )

    return rdf


def lumi_weight(rdf, era):
    if era == "2016preVFP":
        rdf = rdf.Redefine("weight", "weight * 19.52 * 1000.")
    elif era == "2016postVFP":
        rdf = rdf.Redefine("weight", "weight * 16.81 * 1000.")
    elif era == "2017":
        rdf = rdf.Redefine("weight", "weight * 41.48 * 1000.")
    elif era == "2018":
        rdf = rdf.Redefine("weight", "weight * 59.83 * 1000.")
    else:
        sys.exit("Weight calc: lumi: Such an era is not defined: {}".format(era))

    return rdf


def pileup_weight(rdf):
    return rdf.Redefine("weight", "weight * puweight")


def trigger_weight(rdf, channel, process):
    if channel == "et":
        rdf = rdf.Redefine("weight", "weight * trg_wgt_single_ele32orele35")
    elif channel == "mt":
        rdf = rdf.Redefine("weight", "weight * trg_wgt_single_mu24ormu27")
    elif channel == "tt":
        rdf = rdf.Redefine("weight", "weight")
    else:
        sys.exit(
            "Weight calc: trigger: Such a channel is not defined: {}".format(channel)
        )

    return rdf


def gen_weight(rdf, sample):
    number_generated_events_weight = 1.0 / float(sample["nevents"])
    cross_section_per_event_weight = float(sample["xsec"])
    negative_events_fraction = float(sample["generator_weight"])
    rdf = rdf.Define(
        "numberGeneratedEventsWeight",
        "(float){}".format(number_generated_events_weight),
    )
    rdf = rdf.Define(
        "crossSectionPerEventWeight", "(float){}".format(cross_section_per_event_weight)
    )
    rdf = rdf.Define(
        "negativeEventsFraction", "(float){}".format(negative_events_fraction)
    )

    return rdf.Redefine(
        "weight",
        "weight * numberGeneratedEventsWeight * crossSectionPerEventWeight * (( 1.0 / negativeEventsFraction) * ( ((genWeight<0) * -1) + ((genWeight>=0) * 1)))",
    )


def Z_pt_reweight(rdf, process):
    if process == "DYjets":
        return rdf.Redefine("weight", "weight * ZPtMassReweightWeight")
    else:
        return rdf


def Top_pt_reweight(rdf, process):
    if process == "ttbar":
        return rdf.Redefine("weight", "weight * topPtReweightWeight")
    else:
        return rdf


def emb_gen_weight(rdf):
    return rdf.Redefine(
        "weight",
        "weight * emb_genweight * emb_idsel_wgt_1 * emb_idsel_wgt_2 * emb_triggersel_wgt",
    )
