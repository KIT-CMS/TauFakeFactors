import os
import sys
import gzip
import correctionlib.schemav2 as cs
import rich


var_dict = {
    "pt_2": "tau_pt",
    "pt_1": "lep_pt",
    "mt_1": "lep_mt",
    "m_vis": "m_vis",
    "njets": "njets",
    "nbtag": "nbtags",
    "deltaR_ditaupair": "dR_tau_pair",
    "QCD": "QCD",
    "Wjets": "Wjets",
    "ttbar_J": "ttbar",
}
var_type = {
    "pt_2": "real",
    "pt_1": "real",
    "mt_1": "real",
    "m_vis": "real",
    "njets": "real",
    "nbtag": "real",
    "deltaR_ditaupair": "real",
}
var_discription = {
    "pt_2": "transverse momentum of the hadronic tau in the tau pair; measured between #var_min and #var_max GeV; for higher/lower pt's the edge values are used",
    "pt_1": "transverse momentum of the leptonic tau in the tau pair; measured between #var_min and #var_max GeV; for higher/lower pt's the edge values are used",
    "mt_1": "transverse mass of the lepton in the tau pair; measured between #var_min and #var_max GeV; for higher/lower mt's the edge values are used",
    "m_vis": "invariant mass of the visible tau decay products; measured between #var_min and #var_max GeV; for higher/lower m_vis's the edge values are used",
    "njets": "number of jets in an event; the defined categories are ",
    "nbtag": "number of b-tagged jets in an event; the defined categories are ",
    "deltaR_ditaupair": "spatial distance between the tau pair with deltaR ",
}


def cat_transalator(c):
    strings = c.split("#")
    if len(strings) == 2:
        return (strings[0], strings[1])
    elif len(strings) == 4:
        return ((strings[0], strings[2]), (strings[1], strings[3]))
    else:
        sys.exit("Category splitting is only defined up to 2 dimensions.")


def generate_ff_cs_json(config, ff_functions, fractions, save_path):
    cs_corrections = list()

    if "QCD" in config["target_process"]:
        var = config["target_process"]["QCD"]["var_dependence"]
        binning = config["target_process"]["QCD"]["var_bins"]
        ff_unc_order = {  # the order has to match the one used in helper/functions.py -> fit_function()
            "FFslopeUncUp": 1,
            "FFslopeUncDown": 2,
            "FFnormUncUp": 3,
            "FFnormUncDown": 4,
            "FFmcSubUncUp": 5,
            "FFmcSubUncDown": 6,
        }
        if len(config["target_process"]["QCD"]["split_categories"]) == 1:
            QCD_ff = make_1D_ff(
                "QCD",
                (var, binning),
                ff_functions["QCD"],
                config["target_process"]["QCD"],
                ff_unc_order,
            )
        elif len(config["target_process"]["QCD"]["split_categories"]) == 2:
            QCD_ff = make_2D_ff(
                "QCD",
                (var, binning),
                ff_functions["QCD"],
                config["target_process"]["QCD"],
                ff_unc_order,
            )
        cs_corrections.append(QCD_ff)

    if "Wjets" in config["target_process"]:
        var = config["target_process"]["Wjets"]["var_dependence"]
        binning = config["target_process"]["Wjets"]["var_bins"]
        ff_unc_order = {  # the order has to match the one used in helper/functions.py -> fit_function()
            "FFslopeUncUp": 1,
            "FFslopeUncDown": 2,
            "FFnormUncUp": 3,
            "FFnormUncDown": 4,
            "FFmcSubUncUp": 5,
            "FFmcSubUncDown": 6,
        }
        if len(config["target_process"]["Wjets"]["split_categories"]) == 1:
            Wjets_ff = make_1D_ff(
                "Wjets",
                (var, binning),
                ff_functions["Wjets"],
                config["target_process"]["Wjets"],
                ff_unc_order,
            )
        elif len(config["target_process"]["Wjets"]["split_categories"]) == 2:
            Wjets_ff = make_2D_ff(
                "Wjets",
                (var, binning),
                ff_functions["Wjets"],
                config["target_process"]["Wjets"],
                ff_unc_order,
            )
        cs_corrections.append(Wjets_ff)

    if "ttbar" in config["target_process"]:
        var = config["target_process"]["ttbar"]["var_dependence"]
        binning = config["target_process"]["ttbar"]["var_bins"]
        ff_unc_order = {  # the order has to match the one used in helper/functions.py -> fit_function()
            "FFslopeUncUp": 1,
            "FFslopeUncDown": 2,
            "FFnormUncUp": 3,
            "FFnormUncDown": 4,
        }
        if len(config["target_process"]["ttbar"]["split_categories"]) == 1:
            ttbar_ff = make_1D_ff(
                "ttbar",
                (var, binning),
                ff_functions["ttbar"],
                config["target_process"]["ttbar"],
                ff_unc_order,
            )
        elif len(config["target_process"]["ttbar"]["split_categories"]) == 2:
            ttbar_ff = make_2D_ff(
                "ttbar",
                (var, binning),
                ff_functions["ttbar"],
                config["target_process"]["ttbar"],
                ff_unc_order,
            )
        cs_corrections.append(ttbar_ff)

    if "process_fractions" in config:
        var = config["process_fractions"]["var_dependence"]
        binning = config["process_fractions"]["var_bins"]
        frac_unc = {  # the naming has to match the one used in helper/functions.py -> add_fraction_variations()
            "fracQCDUncUp": "frac_QCD_up",
            "fracQCDUncDown": "frac_QCD_down",
            "fracWjetsUncUp": "frac_Wjets_up",
            "fracWjetsUncDown": "frac_Wjets_down",
            "fracTTbarUncUp": "frac_ttbar_J_up",
            "fracTTbarUncDown": "frac_ttbar_J_down",
        }
        fraction = make_1D_fractions(
            (var, binning), fractions, config["process_fractions"], frac_unc
        )
        cs_corrections.append(fraction)

    cset = cs.CorrectionSet(
        schema_version=2,
        description="Fake factors for tau analysis",
        corrections=cs_corrections,
    )

    with open(save_path + "fake_factors.json", "w") as fout:
        fout.write(cset.json(exclude_unset=True, indent=4))

    with gzip.open(save_path + "fake_factors.json.gz", "wt") as fout:
        fout.write(cset.json(exclude_unset=True, indent=4))


def make_1D_ff(process, variable, ff_functions, config, uncertainties):
    # get categories from config
    cat_inputs = list(config["split_categories"].keys())
    cat_values = [config["split_categories"][cat] for cat in config["split_categories"]]

    ff = cs.Correction(
        name="{}_fake_factors".format(process),
        description="Calculation of the {} part the for data-driven background estimation (fake factors) for misindentified jets as tau leptons in H->tautau analysis.".format(
            process
        ),
        version=1,
        inputs=[
            cs.Variable(
                name=var_dict[variable[0]],
                type=var_type[variable[0]],
                description=var_discription[variable[0]]
                .replace("#var_min", str(min(variable[1])))
                .replace("#var_max", str(max(variable[1]))),
            ),
            cs.Variable(
                name=var_dict[cat_inputs[0]],
                type=var_type[cat_inputs[0]],
                description=var_discription[cat_inputs[0]] + ", ".join(cat_values[0]),
            ),
            cs.Variable(
                name="syst",
                type="string",
                description="Uncertainties from the best fit for the fake factor measurement.",
            ),
        ],
        output=cs.Variable(
            name="{}_ff".format(process),
            type="real",
            description="{} part of the fake factor".format(process),
        ),
        data=cs.Category(
            nodetype="category",
            input="syst",
            content=[
                cs.CategoryItem(
                    key=unc,
                    value=cs.Binning(
                        nodetype="binning",
                        input=var_dict[cat_inputs[0]],
                        edges=config["split_categories_binedges"][cat_inputs[0]],
                        content=[
                            cs.Binning(
                                nodetype="binning",
                                input=var_dict[variable[0]],
                                edges=[
                                    0,
                                    min(variable[1]),
                                    max(variable[1]),
                                    (max(variable[1]) + 1),
                                ],
                                content=[
                                    eval(
                                        ff_functions[cat1][idx].replace(
                                            "x", str(min(variable[1]))
                                        )
                                    ),
                                    cs.Formula(
                                        nodetype="formula",
                                        variables=[var_dict[variable[0]]],
                                        parser="TFormula",
                                        expression=ff_functions[cat1][idx],
                                    ),
                                    eval(
                                        ff_functions[cat1][idx].replace(
                                            "x", str(max(variable[1]))
                                        )
                                    ),
                                ],
                                flow="clamp",
                            )
                            for cat1 in ff_functions
                        ],
                        flow="clamp",
                    ),
                )
                for unc, idx in uncertainties.items()
            ],
            default=cs.Binning(
                nodetype="binning",
                input=var_dict[cat_inputs[0]],
                edges=config["split_categories_binedges"][cat_inputs[0]],
                content=[
                    cs.Binning(
                        nodetype="binning",
                        input=var_dict[variable[0]],
                        edges=[
                            0,
                            min(variable[1]),
                            max(variable[1]),
                            (max(variable[1]) + 1),
                        ],
                        content=[
                            eval(
                                ff_functions[cat1][0].replace(
                                    "x", str(min(variable[1]))
                                )
                            ),
                            cs.Formula(
                                nodetype="formula",
                                variables=[var_dict[variable[0]]],
                                parser="TFormula",
                                expression=ff_functions[cat1][0],
                            ),
                            eval(
                                ff_functions[cat1][0].replace(
                                    "x", str(max(variable[1]))
                                )
                            ),
                        ],
                        flow="clamp",
                    )
                    for cat1 in ff_functions
                ],
                flow="clamp",
            ),
        ),
    )
    rich.print(ff)

    return ff


def make_2D_ff(process, variable, ff_functions, config, uncertainties):
    # get categories from config
    cat_inputs = list(config["split_categories"].keys())
    cat_values = [config["split_categories"][cat] for cat in config["split_categories"]]

    ff = cs.Correction(
        name="{}_fake_factors".format(process),
        description="Calculation of the {} part the for data-driven background estimation (fake factors) for misindentified jets as tau leptons in H->tautau analysis.".format(
            process
        ),
        version=1,
        inputs=[
            cs.Variable(
                name=var_dict[variable[0]],
                type=var_type[variable[0]],
                description=var_discription[variable[0]]
                .replace("#var_min", str(min(variable[1])))
                .replace("#var_max", str(max(variable[1]))),
            ),
            cs.Variable(
                name=var_dict[cat_inputs[0]],
                type=var_type[cat_inputs[0]],
                description=var_discription[cat_inputs[0]] + ", ".join(cat_values[0]),
            ),
            cs.Variable(
                name=var_dict[cat_inputs[1]],
                type=var_type[cat_inputs[1]],
                description=var_discription[cat_inputs[1]] + ", ".join(cat_values[1]),
            ),
            cs.Variable(
                name="syst",
                type="string",
                description="Uncertainties from the best fit for the fake factor measurement.",
            ),
        ],
        output=cs.Variable(
            name="{}_ff".format(process),
            type="real",
            description="{} part of the fake factor".format(process),
        ),
        data=cs.Category(
            nodetype="category",
            input="syst",
            content=[
                cs.CategoryItem(
                    key=unc,
                    value=cs.Binning(
                        nodetype="binning",
                        input=var_dict[cat_inputs[0]],
                        edges=config["split_categories_binedges"][cat_inputs[0]],
                        content=[
                            cs.Binning(
                                nodetype="binning",
                                input=var_dict[cat_inputs[1]],
                                edges=config["split_categories_binedges"][
                                    cat_inputs[1]
                                ],
                                content=[
                                    cs.Binning(
                                        nodetype="binning",
                                        input=var_dict[variable[0]],
                                        edges=[
                                            0,
                                            min(variable[1]),
                                            max(variable[1]),
                                            (max(variable[1]) + 1),
                                        ],
                                        content=[
                                            eval(
                                                ff_functions[cat1][cat2][idx].replace(
                                                    "x", str(min(variable[1]))
                                                )
                                            ),
                                            cs.Formula(
                                                nodetype="formula",
                                                variables=[var_dict[variable[0]]],
                                                parser="TFormula",
                                                expression=ff_functions[cat1][cat2][
                                                    idx
                                                ],
                                            ),
                                            eval(
                                                ff_functions[cat1][cat2][idx].replace(
                                                    "x", str(max(variable[1]))
                                                )
                                            ),
                                        ],
                                        flow="clamp",
                                    )
                                    for cat2 in ff_functions[cat1]
                                ],
                                flow="clamp",
                            )
                            for cat1 in ff_functions
                        ],
                        flow="clamp",
                    ),
                )
                for unc, idx in uncertainties.items()
            ],
            default=cs.Binning(
                nodetype="binning",
                input=var_dict[cat_inputs[0]],
                edges=config["split_categories_binedges"][cat_inputs[0]],
                content=[
                    cs.Binning(
                        nodetype="binning",
                        input=var_dict[cat_inputs[1]],
                        edges=config["split_categories_binedges"][cat_inputs[1]],
                        content=[
                            cs.Binning(
                                nodetype="binning",
                                input=var_dict[variable[0]],
                                edges=[
                                    0,
                                    min(variable[1]),
                                    max(variable[1]),
                                    (max(variable[1]) + 1),
                                ],
                                content=[
                                    eval(
                                        ff_functions[cat1][cat2][0].replace(
                                            "x", str(min(variable[1]))
                                        )
                                    ),
                                    cs.Formula(
                                        nodetype="formula",
                                        variables=[var_dict[variable[0]]],
                                        parser="TFormula",
                                        expression=ff_functions[cat1][cat2][0],
                                    ),
                                    eval(
                                        ff_functions[cat1][cat2][0].replace(
                                            "x", str(max(variable[1]))
                                        )
                                    ),
                                ],
                                flow="clamp",
                            )
                            for cat2 in ff_functions[cat1]
                        ],
                        flow="clamp",
                    )
                    for cat1 in ff_functions
                ],
                flow="clamp",
            ),
        ),
    )
    rich.print(ff)

    return ff


def make_1D_fractions(variable, fractions, config, uncertainties):
    # get categories from config
    cat_inputs = list(config["split_categories"].keys())
    cat_values = [config["split_categories"][cat] for cat in config["split_categories"]]

    frac = cs.Correction(
        name="process_fractions",
        description="Calculation of process contributions (fractions) for the fake factor calculation.",
        version=1,
        inputs=[
            cs.Variable(
                name="process", type="string", description="name of the process"
            ),
            cs.Variable(
                name=var_dict[variable[0]],
                type=var_type[variable[0]],
                description=var_discription[variable[0]]
                .replace("#var_min", str(min(variable[1])))
                .replace("#var_max", str(max(variable[1]))),
            ),
            cs.Variable(
                name=var_dict[cat_inputs[0]],
                type=var_type[cat_inputs[0]],
                description=var_discription[cat_inputs[0]] + ", ".join(cat_values[0]),
            ),
            cs.Variable(
                name="syst",
                type="string",
                description="Uncertainties from the +/-7% variation of the process fractions.",
            ),
        ],
        output=cs.Variable(
            name="fraction", type="real", description="process fraction"
        ),
        data=cs.Category(
            nodetype="category",
            input="syst",
            content=[
                cs.CategoryItem(
                    key=unc,
                    value=cs.Category(
                        nodetype="category",
                        input="process",
                        content=[
                            cs.CategoryItem(
                                key=var_dict[p],
                                value=cs.Binning(
                                    nodetype="binning",
                                    input=var_dict[cat_inputs[0]],
                                    edges=config["split_categories_binedges"][
                                        cat_inputs[0]
                                    ],
                                    content=[
                                        cs.Binning(
                                            nodetype="binning",
                                            input=var_dict[variable[0]],
                                            edges=variable[1],
                                            content=fractions[cat][unc_name][p],
                                            flow="clamp",
                                        )
                                        for cat in fractions
                                    ],
                                    flow="clamp",
                                ),
                            )
                            for p in config["processes"]
                        ],
                    ),
                )
                for unc, unc_name in uncertainties.items()
            ],
            default=cs.Category(
                nodetype="category",
                input="process",
                content=[
                    cs.CategoryItem(
                        key=var_dict[p],
                        value=cs.Binning(
                            nodetype="binning",
                            input=var_dict[cat_inputs[0]],
                            edges=config["split_categories_binedges"][cat_inputs[0]],
                            content=[
                                cs.Binning(
                                    nodetype="binning",
                                    input=var_dict[variable[0]],
                                    edges=variable[1],
                                    content=fractions[cat]["nominal"][p],
                                    flow="clamp",
                                )
                                for cat in fractions
                            ],
                            flow="clamp",
                        ),
                    )
                    for p in config["processes"]
                ],
            ),
        ),
    )
    rich.print(frac)

    return frac


corr_unc_dict = {
    "non_closure": "nonClosure",
}


def generate_corr_cs_json(config, corrections, save_path):

    cs_ff_corrections = list()

    for process in corrections.keys():
        for correction in corrections[process].keys():
            var = config["target_process"][process]["corrections"][correction][
                "var_dependence"
            ]
            corr_dict = corrections[process][correction]
            corr = make_1D_correction(process, var, correction, corr_dict)
            cs_ff_corrections.append(corr)

    cset = cs.CorrectionSet(
        schema_version=2,
        description="Corrections for fake factors for tau analysis",
        corrections=cs_ff_corrections,
    )

    with open(save_path + "FF_corrections.json", "w") as fout:
        fout.write(cset.json(exclude_unset=True, indent=4))

    with gzip.open(save_path + "FF_corrections.json.gz", "wt") as fout:
        fout.write(cset.json(exclude_unset=True, indent=4))


def make_1D_correction(process, variable, correction, corr_dict):

    cs_corr = cs.Correction(
        name="{}_{}_correction".format(process, correction),
        description="Calculation of the {} correction for {} for data-driven background estimation (fake factors) for misindentified jets as tau leptons in H->tautau analysis.".format(
            correction,
            process,
        ),
        version=1,
        inputs=[
            cs.Variable(
                name=var_dict[variable],
                type=var_type[variable],
                description=var_discription[variable]
                .replace("#var_min", str(min(corr_dict["edges"])))
                .replace("#var_max", str(max(corr_dict["edges"]))),
            ),
            cs.Variable(
                name="syst",
                type="string",
                description="Uncertainty from the correction measurement.",
            ),
        ],
        output=cs.Variable(
            name="{}_{}_correction".format(process, correction),
            type="real",
            description="{} correction for {}".format(correction, process),
        ),
        data=cs.Category(
            nodetype="category",
            input="syst",
            content=[
                cs.CategoryItem(
                    key="{}CorrUp".format(corr_unc_dict[correction]),
                    value=cs.Binning(
                        nodetype="binning",
                        input=var_dict[variable],
                        edges=list(corr_dict["edges"]),
                        content=list(corr_dict["up"]),
                        flow="clamp",
                    ),
                ),
                cs.CategoryItem(
                    key="{}CorrDown".format(corr_unc_dict[correction]),
                    value=cs.Binning(
                        nodetype="binning",
                        input=var_dict[variable],
                        edges=list(corr_dict["edges"]),
                        content=list(corr_dict["down"]),
                        flow="clamp",
                    ),
                ),
            ],
            default=cs.Binning(
                nodetype="binning",
                input=var_dict[variable],
                edges=list(corr_dict["edges"]),
                content=list(corr_dict["down"]),
                flow="clamp",
            ),
        ),
    )
    rich.print(cs_corr)

    return cs_corr
