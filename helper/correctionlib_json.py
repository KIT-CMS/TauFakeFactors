import os
import gzip
import correctionlib.schemav2 as cs
import rich
import logging
from typing import List, Dict, Union, Tuple

import configs.general_definitions as gd


def generate_ff_corrlib_json(
    config: Dict[str, Union[str, Dict, List]], ff_functions: Dict[str, Dict[str, Union[Dict[str, str], Dict[str, Dict[str, str]]]]], fractions: Dict[str, Dict[str, Dict[str, List[float]]]], output_path: str, for_corrections: bool = False
) -> None:
    '''
    Function which produces a correctionlib file based on the measured fake factors and fractions (including variations).

    Args:
        config: A dictionary with all the relevant information for the fake factors
        ff_functions: Dictionary of fake factor functions as strings, e.g. ff_function[PROCESS][CATEGORY_1][CATEGORY_2][VARIATION] if dimension of categories is 2
        fractions: Dictionary of fraction values, e.g. fractions[CATEGORY][VARIATION][PROCESS] 
        output_path: Path where the generated correctionlib files should be stored
        for_corrections: Boolean which decides where the produced correctionlib file is stored, this is needed because additional fake factors are needed to calculate some corrections and therefore they are stored in a different place (default: False)

    Return:
        None
    '''
    log = logging.getLogger("ff_calculation")

    corrlib_corrections = list()

    for process in config["target_processes"]:
        var = config["target_processes"][process]["var_dependence"]
        binning = config["target_processes"][process]["var_bins"]

        if len(config["target_processes"][process]["split_categories"]) == 1:
            process_ff = make_1D_ff(
                process=process,
                process_conf=config["target_processes"][process],
                variable_info=(var, binning),
                ff_functions=ff_functions[process],
                # ff_unc_order,
            )
        elif len(config["target_processes"][process]["split_categories"]) == 2:
            process_ff = make_2D_ff(
                process=process,
                process_conf=config["target_processes"][process],
                variable_info=(var, binning),
                ff_functions=ff_functions[process],
                # ff_unc_order,
            )
        corrlib_corrections.append(process_ff)

    if "process_fractions" in config:
        var = config["process_fractions"]["var_dependence"]
        binning = config["process_fractions"]["var_bins"]
        
        frac_unc = dict()  
        for process in config["process_fractions"]["processes"]:
            frac_unc.update(gd.frac_variation_dict[process])

        fraction = make_1D_fractions(
            fraction_conf=config["process_fractions"],
            variable_info=(var, binning), 
            fractions=fractions, 
            uncertainties=frac_unc,
        )
        corrlib_corrections.append(fraction)

    cset = cs.CorrectionSet(
        schema_version=2,
        description="Fake factors for tau analysis",
        corrections=corrlib_corrections,
    )

    if not for_corrections:
        with open(os.path.join(output_path, f"fake_factors_{config['channel']}.json"), "w") as fout:
            fout.write(cset.json(exclude_unset=True, indent=4))

        with gzip.open(os.path.join(output_path, f"fake_factors_{config['channel']}.json.gz"), "wt") as fout:
            fout.write(cset.json(exclude_unset=True, indent=4))
        
        log.info("Correctionlib file saved to " + os.path.join(output_path, f"fake_factors_{config['channel']}.json"))

    elif for_corrections:
        with open(os.path.join(output_path, f"fake_factors_{config['channel']}_for_corrections.json"), "w") as fout:
            fout.write(cset.json(exclude_unset=True, indent=4))

        with gzip.open(os.path.join(output_path, f"fake_factors_{config['channel']}_for_corrections.json.gz"), "wt") as fout:
            fout.write(cset.json(exclude_unset=True, indent=4))

        log.info("Correctionlib file saved to " + os.path.join(output_path, f"fake_factors_{config['channel']}_for_corrections.json"))


def make_1D_ff(process: str, process_conf: Dict[str, Union[Dict, List, str]], variable_info: Tuple[str, List[float]], ff_functions: Dict[str, Dict[str, str]]) -> cs.Correction:
    '''
    Function which produces a correctionlib Correction based on the measured fake factors (including variations) for the case of 1D categories.

    Args:
        process: Name of the process for which the Correction is produced 
        process_conf: A dictionary with all the relevant information for the fake factors of a specific "process"
        variable_info: Tuple with information (name and binning) about the variable the fake factors depends on
        ff_functions: Dictionary of fake factor functions as strings, e.g. ff_function[PROCESS][CATEGORY_1][VARIATION] 

    Return:
        Correction object from correcionlib
    '''
    # get categories from config
    cat_inputs = list(process_conf["split_categories"].keys())
    cat_values = [process_conf["split_categories"][cat] for cat in process_conf["split_categories"]]

    ff = cs.Correction(
        name=f"{process}_fake_factors",
        description=f"Calculation of the {process} part the for data-driven background estimation (fake factors) for misindentified jets as tau leptons in H->tautau analysis.",
        version=1,
        inputs=[
            cs.Variable(
                name=gd.variable_translator[variable_info[0]],
                type=gd.variable_type[variable_info[0]],
                description=gd.variable_discription[variable_info[0]]
                .replace("#var_min", str(min(variable_info[1])))
                .replace("#var_max", str(max(variable_info[1]))),
            ),
            cs.Variable(
                name=gd.variable_translator[cat_inputs[0]],
                type=gd.variable_type[cat_inputs[0]],
                description=gd.variable_discription[cat_inputs[0]] + ", ".join(cat_values[0]),
            ),
            cs.Variable(
                name="syst",
                type="string",
                description="Uncertainties from the best fit for the fake factor measurement.",
            ),
        ],
        output=cs.Variable(
            name=f"{process}_ff",
            type="real",
            description=f"{process} part of the fake factor",
        ),
        data=cs.Category(
            nodetype="category",
            input="syst",
            content=[
                cs.CategoryItem(
                    key=process+unc_name,
                    value=cs.Binning(
                        nodetype="binning",
                        input=gd.variable_translator[cat_inputs[0]],
                        edges=process_conf["split_categories_binedges"][cat_inputs[0]],
                        content=[
                            cs.Binning(
                                nodetype="binning",
                                input=gd.variable_translator[variable_info[0]],
                                edges=[
                                    0,
                                    min(variable_info[1]),
                                    max(variable_info[1]),
                                    (max(variable_info[1]) + 1),
                                ],
                                content=[
                                    eval(
                                        ff_functions[cat1][unc].replace(
                                            "x", str(min(variable_info[1]))
                                        )
                                    ),
                                    cs.Formula(
                                        nodetype="formula",
                                        variables=[gd.variable_translator[variable_info[0]]],
                                        parser="TFormula",
                                        expression=ff_functions[cat1][unc],
                                    ),
                                    eval(
                                        ff_functions[cat1][unc].replace(
                                            "x", str(max(variable_info[1]))
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
                for unc, unc_name in gd.ff_variation_dict[process].items()
            ],
            default=cs.Binning(
                nodetype="binning",
                input=gd.variable_translator[cat_inputs[0]],
                edges=process_conf["split_categories_binedges"][cat_inputs[0]],
                content=[
                    cs.Binning(
                        nodetype="binning",
                        input=gd.variable_translator[variable_info[0]],
                        edges=[
                            0,
                            min(variable_info[1]),
                            max(variable_info[1]),
                            (max(variable_info[1]) + 1),
                        ],
                        content=[
                            eval(
                                ff_functions[cat1]["nominal"].replace(
                                    "x", str(min(variable_info[1]))
                                )
                            ),
                            cs.Formula(
                                nodetype="formula",
                                variables=[gd.variable_translator[variable_info[0]]],
                                parser="TFormula",
                                expression=ff_functions[cat1]["nominal"],
                            ),
                            eval(
                                ff_functions[cat1]["nominal"].replace(
                                    "x", str(max(variable_info[1]))
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


def make_2D_ff(process: str, process_conf: Dict[str, Union[Dict, List, str]], variable_info: Tuple[str, List[float]], ff_functions: Dict[str, Dict[str, Dict[str, str]]]) -> cs.Correction:
    '''
    Function which produces a correctionlib Correction based on the measured fake factors (including variations) for the case of 2D categories.

    Args:
        process: Name of the process for which the Correction is produced 
        process_conf: A dictionary with all the relevant information for the fake factors of a specific "process"
        variable_info: Tuple with information (name and binning) about the variable the fake factors depends on
        ff_functions: Dictionary of fake factor functions as strings, e.g. ff_function[PROCESS][CATEGORY_1][CATEGORY_2][VARIATION] 

    Return:
        Correction object from correcionlib
    '''
    # get categories from config
    cat_inputs = list(process_conf["split_categories"].keys())
    cat_values = [process_conf["split_categories"][cat] for cat in process_conf["split_categories"]]

    ff = cs.Correction(
        name=f"{process}_fake_factors",
        description="Calculation of the {process} part the for data-driven background estimation (fake factors) for misindentified jets as tau leptons in H->tautau analysis.",
        version=1,
        inputs=[
            cs.Variable(
                name=gd.variable_translator[variable_info[0]],
                type=gd.variable_type[variable_info[0]],
                description=gd.variable_discription[variable_info[0]]
                .replace("#var_min", str(min(variable_info[1])))
                .replace("#var_max", str(max(variable_info[1]))),
            ),
            cs.Variable(
                name=gd.variable_translator[cat_inputs[0]],
                type=gd.variable_type[cat_inputs[0]],
                description=gd.variable_discription[cat_inputs[0]] + ", ".join(cat_values[0]),
            ),
            cs.Variable(
                name=gd.variable_translator[cat_inputs[1]],
                type=gd.variable_type[cat_inputs[1]],
                description=gd.variable_discription[cat_inputs[1]] + ", ".join(cat_values[1]),
            ),
            cs.Variable(
                name="syst",
                type="string",
                description="Uncertainties from the best fit for the fake factor measurement.",
            ),
        ],
        output=cs.Variable(
            name=f"{process}_ff",
            type="real",
            description=f"{process} part of the fake factor",
        ),
        data=cs.Category(
            nodetype="category",
            input="syst",
            content=[
                cs.CategoryItem(
                    key=process+unc_name,
                    value=cs.Binning(
                        nodetype="binning",
                        input=gd.variable_translator[cat_inputs[0]],
                        edges=process_conf["split_categories_binedges"][cat_inputs[0]],
                        content=[
                            cs.Binning(
                                nodetype="binning",
                                input=gd.variable_translator[cat_inputs[1]],
                                edges=process_conf["split_categories_binedges"][
                                    cat_inputs[1]
                                ],
                                content=[
                                    cs.Binning(
                                        nodetype="binning",
                                        input=gd.variable_translator[variable_info[0]],
                                        edges=[
                                            0,
                                            min(variable_info[1]),
                                            max(variable_info[1]),
                                            (max(variable_info[1]) + 1),
                                        ],
                                        content=[
                                            eval(
                                                ff_functions[cat1][cat2][unc].replace(
                                                    "x", str(min(variable_info[1]))
                                                )
                                            ),
                                            cs.Formula(
                                                nodetype="formula",
                                                variables=[gd.variable_translator[variable_info[0]]],
                                                parser="TFormula",
                                                expression=ff_functions[cat1][cat2][
                                                    unc
                                                ],
                                            ),
                                            eval(
                                                ff_functions[cat1][cat2][unc].replace(
                                                    "x", str(max(variable_info[1]))
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
                for unc, unc_name in gd.ff_variation_dict[process].items()
            ],
            default=cs.Binning(
                nodetype="binning",
                input=gd.variable_translator[cat_inputs[0]],
                edges=process_conf["split_categories_binedges"][cat_inputs[0]],
                content=[
                    cs.Binning(
                        nodetype="binning",
                        input=gd.variable_translator[cat_inputs[1]],
                        edges=process_conf["split_categories_binedges"][cat_inputs[1]],
                        content=[
                            cs.Binning(
                                nodetype="binning",
                                input=gd.variable_translator[variable_info[0]],
                                edges=[
                                    0,
                                    min(variable_info[1]),
                                    max(variable_info[1]),
                                    (max(variable_info[1]) + 1),
                                ],
                                content=[
                                    eval(
                                        ff_functions[cat1][cat2]["nominal"].replace(
                                            "x", str(min(variable_info[1]))
                                        )
                                    ),
                                    cs.Formula(
                                        nodetype="formula",
                                        variables=[gd.variable_translator[variable_info[0]]],
                                        parser="TFormula",
                                        expression=ff_functions[cat1][cat2]["nominal"],
                                    ),
                                    eval(
                                        ff_functions[cat1][cat2]["nominal"].replace(
                                            "x", str(max(variable_info[1]))
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


def make_1D_fractions(fraction_conf: Dict[str, Union[Dict, List, str]], variable_info: Tuple[str, List[float]], fractions: Dict[str, Dict[str, Dict[str, List[float]]]], uncertainties: Dict[str, str]) -> cs.Correction:
    '''
    Function which produces a correctionlib Correction based on the measured fractions (including variations).

    Args:
        fraction_conf: A dictionary with all the relevant information for the fraction calculation
        variable_info: Tuple with information (name and binning) about the variable the fractions depends on
        fractions: Dictionary of fraction values, e.g. fractions[CATEGORY][VARIATION][PROCESS]

    Return:
        Correction object from correcionlib
    '''
    # get categories from config
    cat_inputs = list(fraction_conf["split_categories"].keys())
    cat_values = [fraction_conf["split_categories"][cat] for cat in fraction_conf["split_categories"]]

    frac = cs.Correction(
        name="process_fractions",
        description="Calculation of process contributions (fractions) for the fake factor calculation.",
        version=1,
        inputs=[
            cs.Variable(
                name="process", type="string", description="name of the process"
            ),
            cs.Variable(
                name=gd.variable_translator[variable_info[0]],
                type=gd.variable_type[variable_info[0]],
                description=gd.variable_discription[variable_info[0]]
                .replace("#var_min", str(min(variable_info[1])))
                .replace("#var_max", str(max(variable_info[1]))),
            ),
            cs.Variable(
                name=gd.variable_translator[cat_inputs[0]],
                type=gd.variable_type[cat_inputs[0]],
                description=gd.variable_discription[cat_inputs[0]] + ", ".join(cat_values[0]),
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
                    key=unc_name,
                    value=cs.Category(
                        nodetype="category",
                        input="process",
                        content=[
                            cs.CategoryItem(
                                key=gd.variable_translator[p],
                                value=cs.Binning(
                                    nodetype="binning",
                                    input=gd.variable_translator[cat_inputs[0]],
                                    edges=fraction_conf["split_categories_binedges"][
                                        cat_inputs[0]
                                    ],
                                    content=[
                                        cs.Binning(
                                            nodetype="binning",
                                            input=gd.variable_translator[variable_info[0]],
                                            edges=variable_info[1],
                                            content=fractions[cat][unc][p],
                                            flow="clamp",
                                        )
                                        for cat in fractions
                                    ],
                                    flow="clamp",
                                ),
                            )
                            for p in fraction_conf["processes"]
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
                        key=gd.variable_translator[p],
                        value=cs.Binning(
                            nodetype="binning",
                            input=gd.variable_translator[cat_inputs[0]],
                            edges=fraction_conf["split_categories_binedges"][cat_inputs[0]],
                            content=[
                                cs.Binning(
                                    nodetype="binning",
                                    input=gd.variable_translator[variable_info[0]],
                                    edges=variable_info[1],
                                    content=fractions[cat]["nominal"][p],
                                    flow="clamp",
                                )
                                for cat in fractions
                            ],
                            flow="clamp",
                        ),
                    )
                    for p in fraction_conf["processes"]
                ],
            ),
        ),
    )
    rich.print(frac)

    return frac


def generate_corr_cs_json(config, save_path, corrections, for_DRtoSR=False):

    cs_ff_corrections = list()

    for process in corrections.keys():
        for correction in corrections[process].keys():
            if "non_closure" in correction and not for_DRtoSR:
                corr_var = correction.split("non_closure_")[1]
                var = config["target_process"][process]["non_closure"][corr_var][
                    "var_dependence"
                ]
            elif "non_closure" in correction and for_DRtoSR:
                corr_var = correction.split("non_closure_")[1]
                var = config["target_process"][process]["DR_SR"]["non_closure"][
                    corr_var
                ]["var_dependence"]
            else:
                var = config["target_process"][process][correction]["var_dependence"]
            corr_dict = corrections[process][correction]
            corr = make_1D_correction(process, var, correction, corr_dict)
            cs_ff_corrections.append(corr)

    cset = cs.CorrectionSet(
        schema_version=2,
        description="Corrections for fake factors for tau analysis",
        corrections=cs_ff_corrections,
    )

    if not for_DRtoSR:
        with open(
            save_path + "/FF_corrections_{}.json".format(config["channel"]), "w"
        ) as fout:
            fout.write(cset.json(exclude_unset=True, indent=4))

        with gzip.open(
            save_path + "/FF_corrections_{}.json.gz".format(config["channel"]), "wt"
        ) as fout:
            fout.write(cset.json(exclude_unset=True, indent=4))

    elif for_DRtoSR:
        with open(
            save_path + "/FF_corrections_{}_for_DRtoSR.json".format(config["channel"]),
            "w",
        ) as fout:
            fout.write(cset.json(exclude_unset=True, indent=4))

        with gzip.open(
            save_path
            + "/FF_corrections_{}_for_DRtoSR.json.gz".format(config["channel"]),
            "wt",
        ) as fout:
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
                name=gd.variable_translator[variable],
                type=gd.variable_type[variable],
                description=gd.variable_discription[variable]
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
                    key=process + "{}CorrUp".format(gd.corr_variation_dict[correction]),
                    value=cs.Binning(
                        nodetype="binning",
                        input=gd.variable_translator[variable],
                        edges=list(corr_dict["edges"]),
                        content=list(corr_dict["up"]),
                        flow="clamp",
                    ),
                ),
                cs.CategoryItem(
                    key=process + "{}CorrDown".format(gd.corr_variation_dict[correction]),
                    value=cs.Binning(
                        nodetype="binning",
                        input=gd.variable_translator[variable],
                        edges=list(corr_dict["edges"]),
                        content=list(corr_dict["down"]),
                        flow="clamp",
                    ),
                ),
            ],
            default=cs.Binning(
                nodetype="binning",
                input=gd.variable_translator[variable],
                edges=list(corr_dict["edges"]),
                content=list(corr_dict["down"]),
                flow="clamp",
            ),
        ),
    )
    rich.print(cs_corr)

    return cs_corr
