"""
Collection of helpful functions for other scripts
"""

import sys
import os
import glob
import ROOT


def check_for_empty_file(path, tree):
    f = ROOT.TFile.Open(path)
    return bool(f.Get(tree))


def get_samples(config, sample):
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


def QCD_SS_estimate(hists):
    qcd = hists["data"].Clone()

    for sample in hists:
        if sample not in ["data", "data_subtracted", "QCD"] and "_T" not in sample:
            qcd.Add(hists[sample], -1)
    # check for negative bins
    for i in range(qcd.GetNbinsX()):
        if qcd.GetBinContent(i) < 0.0:
            qcd.SetBinContent(i, 0.0)
    return qcd


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


def fit_function():
    # bindata = ROOT.Fit.BinData()
    # print(bindata)
    # ROOT.Fit.FillData(bindata, ratio)
    # print(bindata.GetValue())
    # cl68 = fit.GetConfidenceIntervals(bindata, cl=0.68)
    # print(cl68)

    # fit = ROOT.TH1D("fit", "Fit with .95 conf.band", 4, array.array("d", config["tau_pt_bins"]))
    # ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(fit, 0.68)
    # # fit.SetStats(False)
    # # fit.SetFillStyle(1001)
    # fit.SetFillColorAlpha(2, 0.2)
    # fit.SetAxisRange(0, 0.5, "Y")

    # #cl68.Draw("E3 SAME")
    # fit.SetTitle("channel: et, {} jets".format(njets.replace(" ","")))
    # fit.Draw("E3 SAME")
    pass
