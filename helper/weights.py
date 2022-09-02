import sys

def iso_weight(rdf, channel):
    if channel == "et":
        rdf = rdf.Redefine("weight", "weight * iso_wgt_ele_1")
    elif channel == "mt":
        rdf = rdf.Redefine("weight", "weight * iso_wgt_mu_1")
    elif channel == "tt":
        rdf = rdf.Redefine("weight", "weight")
    else:
        sys.exit("Weight calc: iso: Such a channel is not defined: {}".format(channel))
    
    return rdf

def id_weight(rdf, channel):
    if channel == "et":
        rdf = rdf.Redefine("weight", "weight * id_wgt_ele_1")
    elif channel == "mt":
        rdf = rdf.Redefine("weight", "weight * id_wgt_mu_1")
    elif channel == "tt":
        rdf = rdf.Redefine("weight", "weight")
    else:
        sys.exit("Weight calc: id: Such a channel is not defined: {}".format(channel))
    
    return rdf

def tau_id_vsJet_weight(rdf, channel):
    if channel == "et":
        rdf = rdf.Redefine("weight", "weight * ((gen_match_2==5) * ((id_tau_vsJet_Medium_2>0.5)*id_wgt_tau_vsJet_Medium_2 + (id_tau_vsJet_Medium_2<0.5)*(id_tau_vsJet_VVVLoose_2>0.5)*id_wgt_tau_vsJet_VVVLoose_2 + (id_tau_vsJet_VVVLoose_2<0.5)) + (gen_match_2!=5))")
    elif channel == "mt":
        rdf = rdf.Redefine("weight", "weight * ((gen_match_2==5) * ((id_tau_vsJet_Medium_2>0.5)*id_wgt_tau_vsJet_Medium_2 + (id_tau_vsJet_Medium_2<0.5)*(id_tau_vsJet_VVVLoose_2>0.5)*id_wgt_tau_vsJet_VVVLoose_2 + (id_tau_vsJet_VVVLoose_2<0.5)) + (gen_match_2!=5))")
    elif channel == "tt":
        rdf = rdf.Redefine("weight", "weight * ((gen_match_1==5) * ((id_tau_vsJet_Medium_1>0.5)*id_wgt_tau_vsJet_Medium_1 + (id_tau_vsJet_Medium_1<0.5)*(id_tau_vsJet_VVVLoose_1>0.5)*id_wgt_tau_vsJet_VVVLoose_1 + (id_tau_vsJet_VVVLoose_1<0.5)) + (gen_match_1!=5))")
        rdf = rdf.Redefine("weight", "weight * ((gen_match_2==5) * ((id_tau_vsJet_Medium_2>0.5)*id_wgt_tau_vsJet_Medium_2 + (id_tau_vsJet_Medium_2<0.5)*(id_tau_vsJet_VVVLoose_2>0.5)*id_wgt_tau_vsJet_VVVLoose_2 + (id_tau_vsJet_VVVLoose_2<0.5)) + (gen_match_2!=5))")
    else:
        sys.exit("Weight calc: tau id vsJet: Such a channel is not defined: {}".format(channel))

    return rdf

def tau_id_vsMu_weight(rdf, channel):
    if channel == "et":
        rdf = rdf.Redefine("weight", "weight * ((gen_match_2==5) * ((id_tau_vsMu_VLoose_2>0.5)*id_wgt_tau_vsMu_VLoose_2 + (id_tau_vsMu_VLoose_2<0.5)) + (gen_match_2!=5))")
    elif channel == "mt":
        rdf = rdf.Redefine("weight", "weight * ((gen_match_2==5) * ((id_tau_vsMu_Tight_2>0.5)*id_wgt_tau_vsMu_Tight_2 + (id_tau_vsMu_Tight_2<0.5)) + (gen_match_2!=5))")
    elif channel == "tt":
        rdf = rdf.Redefine("weight", "weight * ((gen_match_1==5) * ((id_tau_vsMu_VLoose_1>0.5)*id_wgt_tau_vsMu_VLoose_1 + (id_tau_vsMu_VLoose_1<0.5)) + (gen_match_1!=5))")
        rdf = rdf.Redefine("weight", "weight * ((gen_match_2==5) * ((id_tau_vsMu_VLoose_2>0.5)*id_wgt_tau_vsMu_VLoose_2 + (id_tau_vsMu_VLoose_2<0.5)) + (gen_match_2!=5))")
    else:
        sys.exit("Weight calc: tau id vsMu: Such a channel is not defined: {}".format(channel))

    return rdf

def tau_id_vsEle_weight(rdf, channel):
    if channel == "et":
        rdf = rdf.Redefine("weight", "weight * ((gen_match_2==5) * ((id_tau_vsEle_Tight_2>0.5)*id_wgt_tau_vsEle_Tight_2 + (id_tau_vsEle_Tight_2<0.5)) + (gen_match_2!=5))")
    elif channel == "mt":
        rdf = rdf.Redefine("weight", "weight * ((gen_match_2==5) * ((id_tau_vsEle_VVLoose_2>0.5)*id_wgt_tau_vsEle_VVLoose_2 + (id_tau_vsEle_VVLoose_2<0.5)) + (gen_match_2!=5))")
    elif channel == "tt":
        rdf = rdf.Redefine("weight", "weight * ((gen_match_1==5) * ((id_tau_vsEle_VVLoose_1>0.5)*id_wgt_tau_vsEle_VVLoose_1 + (id_tau_vsEle_VVLoose_1<0.5)) + (gen_match_1!=5))")
        rdf = rdf.Redefine("weight", "weight * ((gen_match_2==5) * ((id_tau_vsEle_VVLoose_2>0.5)*id_wgt_tau_vsEle_VVLoose_2 + (id_tau_vsEle_VVLoose_2<0.5)) + (gen_match_2!=5))")
    else:
        sys.exit("Weight calc: tau id vsEle: Such a channel is not defined: {}".format(channel))

    return rdf

def lumi_weight(rdf, era):
    if era == "2016":
        rdf = rdf.Redefine("weight", "weight * 36.33 * 1000.")
    elif era == "2017":
        rdf = rdf.Redefine("weight", "weight * 41.529 * 1000.")
    elif era == "2018":
        rdf = rdf.Redefine("weight", "weight * 59.74 * 1000.")
    else:
        sys.exit("Weight calc: lumi: Such an era is not defined: {}".format(era))

    return rdf

def pileup_weight(rdf):
    return rdf.Redefine("weight", "weight * puweight")

def trigger_weight(rdf, channel, process):
    if channel == "et":
        rdf = rdf.Redefine("weight", "weight * trg_wgt_single_ele32orele35")
    elif channel == "mt":
        if process == "embedding":
            rdf = rdf.Redefine("weight", "weight * trg_wgtsingle_mu24Ormu27")
        else:
            rdf = rdf.Redefine("weight", "weight * trg_wgt_single_mu24ormu27")
    elif channel == "tt":
        rdf = rdf.Redefine("weight", "weight")
    else:
        sys.exit("Weight calc: trigger: Such a channel is not defined: {}".format(channel))

    return rdf

def gen_weight(rdf, sample):
    number_generated_events_weight = 1. / float(sample["nevents"])
    cross_section_per_event_weight = float(sample["xsec"])
    rdf = rdf.Define("numberGeneratedEventsWeight", "(float){}".format(number_generated_events_weight))
    rdf = rdf.Define("crossSectionPerEventWeight", "(float){}".format(cross_section_per_event_weight))

    return rdf.Redefine("weight", "weight * numberGeneratedEventsWeight * crossSectionPerEventWeight")

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
    return rdf.Redefine("weight", "weight * emb_genweight * emb_idsel_wgt_1 * emb_idsel_wgt_2 * emb_triggersel_wgt")