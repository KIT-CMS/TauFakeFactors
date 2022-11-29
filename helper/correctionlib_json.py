import os
import sys
import gzip
import correctionlib.schemav2 as cs
import rich

import helper.functions as func

cat_def_dict = {
    "njets": {"==0": 0, "==1": 1, "<=1": 1, ">=2": 2},
    "nbtag": {"==0": 0, "==1": 1, "<=1": 1, ">=2": 2},
}

var_dict = {
    "pt_2": "tau_pt", 
    "mt_1": "lep_mt",
    "njets": "njets",
    "nbtag": "nbtags",
    "QCD": "QCD",
    "Wjets": "Wjets",
    "ttbar_J": "ttbar",
}
var_type = {
    "pt_2": "real", 
    "mt_1": "real",
    "njets": "int",
    "nbtag": "int",
}
var_discription = {
    "pt_2": "transverse momentum of the hadronic tau in the tau pair; the fake factors are measured between #var_min and #var_max GeV; for higher/lower pt's the edge values are used",
    "mt_1": "transverse mass of the lepton in the tau pair; the fractions are measured between #var_min and #var_max GeV; for higher/lower mt's the edge values are used",
    "njets": "number of jets in an event; the defined categories are ",
    "nbtag": "number of b-tagged jets in an event; the defined categories are ",
}

def cat_transalator(c):
    strings = c.split("#")
    if len(strings) == 2:
        return (strings[0], strings[1])
    elif len(strings) == 4:
        return (strings[0], cat_def_dict[strings[0]][strings[1]], strings[2], cat_def_dict[strings[2]][strings[3]])
    else:
        sys.exit("Category splitting is only defined up to 2 dimensions.")

def generate_cs_json(config, ff_functions, fractions):    
    var = config["target_process"]["QCD"]["var_dependence"]
    binning = config["target_process"]["QCD"]["var_bins"]
    QCD_ff = make_1D_ff("QCD", (var, binning), ff_functions["QCD"])

    var = config["target_process"]["Wjets"]["var_dependence"]
    binning = config["target_process"]["Wjets"]["var_bins"]
    Wjets_ff = make_1D_ff("Wjets", (var, binning), ff_functions["Wjets"])
    
    var = config["target_process"]["ttbar"]["var_dependence"]
    binning = config["target_process"]["ttbar"]["var_bins"]
    ttbar_ff = make_1D_ff("ttbar", (var, binning), ff_functions["ttbar"])
    
    var = config["process_fractions"]["var_dependence"]
    binning = config["process_fractions"]["var_bins"]
    fraction = make_1D_fractions((var, binning), fractions)

    cset = cs.CorrectionSet(
        schema_version=2,
        description="Fake factors for tau analysis",
        corrections=[
            QCD_ff,
            Wjets_ff,
            ttbar_ff,
            fraction,
        ],
    )
    

    func.check_output_path(os.getcwd() + "/workdir/" + config["workdir_name"] + "/" + config["channel"])
    save_path = "workdir/{}/{}/".format(config["workdir_name"], config["channel"])

    with open(save_path + "fake_factors.json", "w") as fout:
        fout.write(cset.json(exclude_unset=True, indent=4))

    with gzip.open(save_path + "fake_factors.json.gz", "wt") as fout:
        fout.write(cset.json(exclude_unset=True, indent=4))


def make_1D_ff(process, variable, ff_functions):
    # read out first category from the fraction dict to get the category variable 
    cat_input = cat_transalator(list(ff_functions.keys())[0])[1]
    # get the cuts on the categories
    cat_list = list()
    for cat in ff_functions:
        cat_list.append(cat_transalator(cat)[0])

    ff = cs.Correction(
        name="{}_fake_factors".format(process),
        description="Calculation of the {} part the for data-driven background estimation (fake factors) for misindentified jets as tau leptons in H->tautau analysis.".format(process),
        version=1,
        inputs=[
            cs.Variable(name=var_dict[variable[0]], type=var_type[variable[0]], description=var_discription[variable[0]].replace("#var_min", str(min(variable[1]))).replace("#var_max", str(max(variable[1])))),
            cs.Variable(name=var_dict[cat_input], type=var_type[cat_input], description=var_discription[cat_input]+", ".join(cat_list)),
        ],
        output=cs.Variable(name="{}_ff".format(process), type="real", description="{} part of the fake factor".format(process)),
        data=cs.Category(
            nodetype="category",
            input=var_dict[cat_input],
            content=[
                cs.CategoryItem(
                    key=cat_def_dict[cat_input][cat_transalator(cat)[0]],
                    value=cs.Binning(
                        nodetype="binning",
                        input=var_dict[variable[0]],
                        edges=[0,min(variable[1]), max(variable[1]), (max(variable[1])+1)],
                        content=[
                            eval(ff_functions[cat].replace("x", str(min(variable[1])))),
                            cs.Formula(
                                nodetype="formula",
                                variables=[var_dict[variable[0]]],
                                parser="TFormula",
                                expression=ff_functions[cat],
                            ),
                            eval(ff_functions[cat].replace("x", str(max(variable[1])))),
                        ],
                        flow="clamp",
                    ),
                )
                for cat in ff_functions
            ]
        ),
    )
    rich.print(ff)

    return ff

def make_1D_fractions(variable, fractions):
    # read out first category from the fraction dict to get the category variable 
    cat_input = cat_transalator(list(fractions.keys())[0])[1]
    # get the cuts on the categories
    cat_list = list()
    for cat in fractions:
        cat_list.append(cat_transalator(cat)[0])
    process_list = list(fractions[list(fractions.keys())[0]].keys())

    frac = cs.Correction(
        name="process_fractions",
        description="Calculation of process contributions (fractions) for the fake factor calculation.",
        version=1,
        inputs=[
        cs.Variable(name="process", type="string", description="name of the process"),
        cs.Variable(name=var_dict[variable[0]], type=var_type[variable[0]], description=var_discription[variable[0]].replace("#var_min", str(min(variable[1]))).replace("#var_max", str(max(variable[1])))),
        cs.Variable(name=var_dict[cat_input], type=var_type[cat_input], description=var_discription[cat_input]+", ".join(cat_list)),
    ],
        output=cs.Variable(name="fraction", type="real", description="process fraction"),
        data=cs.Category(
            nodetype="category",
            input="process",
            content=[
                cs.CategoryItem(
                    key=var_dict[p],
                    value=cs.Category(
                        nodetype="category",
                        input=var_dict[cat_input],
                        content=[
                            cs.CategoryItem(
                                key=cat_def_dict[cat_input][cat_transalator(cat)[0]],
                                value=cs.Binning(
                                    nodetype="binning",
                                    input=var_dict[variable[0]],                                         
                                    edges=variable[1],
                                    content=fractions[cat][p],
                                    flow="clamp",
                                ),
                            )
                            for cat in fractions
                        ]
                    ),
                )
                for p in process_list
                # cs.CategoryItem(
                #     key="Wjets",
                #     value=cs.Category(
                #         nodetype="category",
                #         input=var_dict[cat_input],
                #         content=[
                #             cs.CategoryItem(
                #                 key=cat_def_dict[cat_input][cat_transalator(cat)[0]],
                #                 value=cs.Binning(
                #                     nodetype="binning",
                #                     input=var_dict[variable[0]],                                         
                #                     edges=variable[1],
                #                     content=fractions[cat]["Wjets"],
                #                     flow="clamp",
                #                 ),
                #             )
                #             for cat in fractions
                #         ]
                #     ),
                # ),
                # cs.CategoryItem(
                #     key="ttbar",
                #     value=cs.Category(
                #         nodetype="category",
                #         input=var_dict[cat_input],
                #         content=[
                #             cs.CategoryItem(
                #                 key=cat_def_dict[cat_input][cat_transalator(cat)[0]],
                #                 value=cs.Binning(
                #                     nodetype="binning",
                #                     input=var_dict[variable[0]],                                         
                #                     edges=variable[1],
                #                     content=fractions[cat]["ttbar_J"],
                #                     flow="clamp",
                #                 ),
                #             )
                #             for cat in fractions
                #         ]
                #     ),
                # ),
            ]
        ),
    )
    rich.print(frac)

    return frac