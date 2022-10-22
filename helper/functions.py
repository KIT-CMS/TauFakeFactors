import ROOT
import glob
import sys
import os


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
            print("File {} has no {} tree. Removed from file list.".format(f, config["tree"]))
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