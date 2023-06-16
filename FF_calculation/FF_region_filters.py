import helper.filters as filters
import helper.weights as weights


def region_filter(rdf, channel, cut_config, sample):

    if "tau_pair_sign" in cut_config:
        rdf = filters.tau_pair_sign_cut(rdf, channel, cut_config["tau_pair_sign"])
    if "deltaR_ditaupair" in cut_config:
        rdf = filters.tau_pair_dR_cut(rdf, channel, cut_config["deltaR_ditaupair"])
    if "lep_mt" in cut_config:
        rdf = filters.lep_mt_cut(rdf, channel, cut_config)
    if "lep_iso" in cut_config:
        rdf = filters.lep_iso_cut(rdf, channel, cut_config)
    if "no_extra_lep" in cut_config:
        rdf = filters.no_extra_lep_cut(rdf, channel, cut_config)
    if "had_tau_id_vs_jet" in cut_config or "had_tau_id_vs_jet_1" in cut_config:
        if sample not in ["data"]:
            rdf = weights.apply_tau_id_vsJet_weight(rdf, channel, cut_config)
        rdf = filters.had_tau_id_vsJet_cut(rdf, channel, cut_config)
    if "njets" in cut_config:
        rdf = filters.jet_number_cut(rdf, channel, cut_config["njets"])
    if "nbtag" in cut_config:
        if sample not in ["data", "embedding"]:
            rdf = weights.apply_btag_weight(rdf)
        rdf = filters.bjet_number_cut(rdf, channel, cut_config["nbtag"])

    return rdf
