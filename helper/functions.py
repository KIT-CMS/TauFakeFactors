"""
Collection of helpful functions for other scripts
"""

import sys
import os
import glob
import array
import ROOT


def check_for_empty_file(path, tree):
    f = ROOT.TFile.Open(path)
    return bool(f.Get(tree))


def get_ntuples(config, sample):
    sample_paths = (
        config["ntuple_path"]
        + "/"
        + config["era"]
        + "/"
        + sample
        + "/"
        + config["channel"]
        + "/*.root"
    )
    print(
        "The following files are loaded for era: {}, channel: {}".format(
            config["era"], config["channel"]
        )
    )
    print("-" * 50)

    sample_path_list = glob.glob(sample_paths)
    cleaned_sample_path_list = glob.glob(sample_paths)

    for f in sample_path_list:
        if check_for_empty_file(f, config["tree"]):
            print(f)
        else:
            print(
                "File {} has no {} tree. Removed from file list.".format(
                    f, config["tree"]
                )
            )
            cleaned_sample_path_list.remove(f)
    print("-" * 50)

    if cleaned_sample_path_list == []:
        sys.exit("Input files: No files found for {}".format(sample))

    return cleaned_sample_path_list

def get_samples(config):
    sample_paths = (
        config["file_path"]
        + "/preselection/"
        + config["era"]
        + "/"
        + config["channel"]
        + "/*.root"
    )
    print(
        "The following files are loaded for era: {}, channel: {}".format(
            config["era"], config["channel"]
        )
    )
    print("-" * 50)
    sample_path_list = glob.glob(sample_paths)
    tmp_list = glob.glob(sample_paths)

    for f in tmp_list:
        sample = f.rsplit("/")[-1].rsplit(".")[0]
        if config["use_embedding"] and "_T" in sample:
            sample_path_list.remove(f)
        elif not config["use_embedding"] and sample == "embedding":
            sample_path_list.remove(f)
    
    for f in sample_path_list:    
        print(f)
    print("-" * 50)

    return sample_path_list


def check_output_path(output_path):
    if not os.path.exists(output_path):
        print(
            "Output directory does not exist! Making directory {}".format(output_path)
        )
        print("-" * 50)
        os.makedirs(output_path, exist_ok=True)


def get_output_name(output_path, process, tau_gen_mode, idx=None):
    if tau_gen_mode == "all":
        tau_gen_mode = ""
    else:
        tau_gen_mode = "_" + tau_gen_mode

    if idx is not None:
        return output_path + "/" + process + tau_gen_mode + "_" + str(idx) + ".root"
    else:
        return output_path + "/" + process + tau_gen_mode + ".root"


def get_split_combinations(categories):
    combinations = list()
    split_vars = list(categories.keys())
    
    if len(split_vars) == 1:
        for n in categories[split_vars[0]]:
            combinations.append({split_vars[0]: n})
    elif len(split_vars) == 2:
        for n in categories[split_vars[0]]:
            for m in categories[split_vars[1]]:
                combinations.append({split_vars[0]: n, split_vars[1]: m})
    else:
        sys.exit("Category splitting is only defined up to 2 dimensions.")

    return tuple(split_vars), combinations


def QCD_SS_estimate(hists):
    qcd = hists["data"].Clone()

    for sample in hists:
        if sample not in ["data", "data_subtracted", "QCD"]:
            qcd.Add(hists[sample], -1)
    # check for negative bins
    for i in range(qcd.GetNbinsX()):
        if qcd.GetBinContent(i) < 0.0:
            qcd.SetBinContent(i, 0.0)
    return qcd

def calc_fraction(hists, target, processes):
    mc = hists[processes[0]].Clone()
    for p in processes:
        if p != processes[0]:
            mc.Add(hists[p])
    frac = hists[target].Clone()
    frac.Divide(mc)
    return frac

def calculate_QCD_FF(SRlike, ARlike):
    ratio = SRlike["data_subtracted"].Clone()
    ratio.Divide(ARlike["data_subtracted"])
    return ratio


def calculate_Wjets_FF(SRlike, ARlike):
    ratio = SRlike["data_subtracted"].Clone()
    ratio.Divide(ARlike["data_subtracted"])
    return ratio


def calculate_ttbar_FF(SR, AR, SRlike, ARlike):
    ratio_mc = SR["ttbar_J"].Clone()
    ratio_mc.Divide(AR["ttbar_J"])

    ratio_DR_data = SRlike["data_subtracted"].Clone()
    ratio_DR_data.Divide(ARlike["data_subtracted"])

    ratio_DR_mc = SRlike["ttbar_J"].Clone()
    ratio_DR_mc.Divide(ARlike["ttbar_J"])

    sf = ratio_DR_data.GetMaximum() / ratio_DR_mc.GetMaximum()
    ratio_mc.Scale(sf)

    return ratio_mc


def fit_function(ff_hist):
    nbins = ff_hist.GetNbinsX() 

    x = list()
    y = list()
    error_y_up = list()
    error_y_down = list()
    for nbin in range(nbins):
        x.append(ff_hist.GetBinCenter(nbin+1))
        y.append(ff_hist.GetBinContent(nbin+1))
        error_y_up.append(ff_hist.GetBinErrorUp(nbin+1))
        error_y_down.append(ff_hist.GetBinErrorLow(nbin+1))

    x = array.array("d", x)
    y = array.array("d", y)
    error_y_up = array.array("d", error_y_up)
    error_y_down = array.array("d", error_y_down)

    graph = ROOT.TGraphAsymmErrors(nbins, x, y, 0, 0, error_y_down, error_y_up)

    fit = graph.Fit("pol1", "SFN")
    print("-" * 50)
    fitted_func = lambda pt: fit.Parameter(1) * pt + fit.Parameter(0)
    cs_expression = "{}*x+{}".format(fit.Parameter(1), fit.Parameter(0))

    cl68 = fit.GetConfidenceIntervals(cl=0.68)
    chi2 = fit.Chi2()
    dof = fit.Ndf()

    y_fit = list()
    error_y_fit_up = list()
    error_y_fit_down = list()

    for nbin in range(nbins):
        y_fit.append(fitted_func(x[nbin]))
        error_y_fit_up.append(cl68[nbin])
        error_y_fit_down.append(cl68[nbin])

    y_fit = array.array("d", y_fit)
    error_y_fit_up = array.array("d", error_y_fit_up)
    error_y_fit_down = array.array("d", error_y_fit_down)

    fit_graph = ROOT.TGraphAsymmErrors(nbins, x, y_fit, 0, 0, error_y_fit_down, error_y_fit_up)
    
    return fit_graph, chi2, dof, cs_expression

def get_yields_from_hists(hists, processes):
    fracs = dict()
    categories = hists.keys()
    for cat in categories:
        nfracs = dict()
        for p in processes:
            h = hists[cat][p]
            l = list()
            for b in range(h.GetNbinsX()):
                l.append(h.GetBinContent(b+1)) 
            nfracs[p] = l
        fracs[cat] = nfracs
    
    return fracs


def calculating_FF(rdf):
    import correctionlib
    correctionlib.register_pyroot_binding()

    ROOT.gInterpreter.ProcessLine('''
        float calc_ff(float pt, int njets, float mt, int nbtag) {
        if (njets > 2) {
            njets = 2;
        }
        if (nbtag > 2) {
            nbtag = 2;
        }
        auto qcd = correction::CorrectionSet::from_file("workdir/json_test_v3/et/fake_factors.json")->at("QCD_fake_factors");
        auto wjets = correction::CorrectionSet::from_file("workdir/json_test_v3/et/fake_factors.json")->at("Wjets_fake_factors");
        auto ttbar = correction::CorrectionSet::from_file("workdir/json_test_v3/et/fake_factors.json")->at("ttbar_fake_factors");
        auto frac = correction::CorrectionSet::from_file("workdir/json_test_v3/et/fake_factors.json")->at("process_fractions");
        float qcd_ff = 1.;
        float wjets_ff = 1.;
        float ttbar_ff = 1.;
        qcd_ff = qcd->evaluate({pt, njets});
        wjets_ff = wjets->evaluate({pt, njets});
        if (njets == 0) {
            njets = 1;
        }
        ttbar_ff = ttbar->evaluate({pt, njets});
        float qcd_frac = 1.;
        float wjets_frac = 1.;
        float ttbar_frac = 1.;
        qcd_frac = frac->evaluate({"QCD", mt, nbtag});
        wjets_frac = frac->evaluate({"Wjets", mt, nbtag});
        ttbar_frac = frac->evaluate({"ttbar", mt, nbtag});
        float full_ff = qcd_frac*qcd_ff+wjets_frac*wjets_ff+ttbar_frac*ttbar_ff;
        return full_ff;
        }''')

    rdf = rdf.Define(
        "fake_factor",
        "calc_ff(pt_2,njets,mt_1,nbtag)",
        )
    
    return rdf