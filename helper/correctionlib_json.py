import gzip
import json
import os
from typing import Any, Dict, List, Tuple, Union

import correctionlib.schemav2 as cs
import numpy as np
import rich

import configs.general_definitions as gd


def get_edges_and_content(
    item: Union[str, Any],
    variable_info: Tuple[str, list],
):
    # import ipdb; ipdb.set_trace()
    if isinstance(item, str):
        return {
            "edges": [
                min(variable_info[1]),
                max(variable_info[1]),
            ],
            "content": [
                cs.Formula(
                    nodetype="formula",
                    variables=[variable_info[0]],
                    parser="TFormula",
                    expression=item,
                    parameters=None,
                ),
            ],
        }
    return {
        "edges": variable_info[1],
        "content": item,
    }


def write_json(path: str, item: Union[cs.CorrectionSet, cs.Correction]) -> None:
    _open, _mode = open, "w"
    if path.endswith(".gz"):
        _open, _mode = gzip.open, "wt"
    with _open(path, _mode) as fout:
        fout.write(json.dumps(item.model_dump(exclude_unset=True), indent=4))


def generate_ff_corrlib_json(
    config: Dict[str, Union[str, Dict, List]],
    ff_functions: Dict[
        str, Dict[str, Union[Dict[str, str], Dict[str, Dict[str, str]]]]
    ],
    fractions: Union[Dict[str, Dict[str, Dict[str, List[float]]]], None],
    fractions_subleading: Union[Dict[str, Dict[str, Dict[str, List[float]]]], None],
    output_path: str,
    for_corrections: bool = False,
) -> None:
    """
    Function which produces a correctionlib file based on the measured fake factors and fractions (including variations).

    Args:
        config: A dictionary with all the relevant information for the fake factors
        ff_functions: Dictionary of fake factor functions as strings, e.g. ff_function[PROCESS][CATEGORY_1][CATEGORY_2][VARIATION] if dimension of categories is 2
        fractions: Dictionary of fraction values, e.g. fractions[CATEGORY][VARIATION][PROCESS]
        fractions_subleading: Second dictionary of fraction values, e.g. fractions[CATEGORY][VARIATION][PROCESS] (can be relvant for tt channel)
        output_path: Path where the generated correctionlib files should be stored
        for_corrections: Boolean which decides where the produced correctionlib file is stored, this is needed because additional fake factors are needed to calculate some corrections and therefore they are stored in a different place (default: False)

    Return:
        None
    """
    corrlib_corrections = list()

    for process in config["target_processes"]:
        if process in ff_functions:
            var = config["target_processes"][process]["var_dependence"]
            binning = config["target_processes"][process]["var_bins"]

            if len(config["target_processes"][process]["split_categories"]) == 1:
                process_ff = make_1D_ff(
                    process=process,
                    process_conf=config["target_processes"][process],
                    variable_info=(var, binning),
                    ff_functions=ff_functions[process],
                )
            elif len(config["target_processes"][process]["split_categories"]) == 2:
                process_ff = make_2D_ff(
                    process=process,
                    process_conf=config["target_processes"][process],
                    variable_info=(var, binning),
                    ff_functions=ff_functions[process],
                )
            corrlib_corrections.append(process_ff)

    if "process_fractions" in config and fractions is not None:
        var = config["process_fractions"]["var_dependence"]
        binning = config["process_fractions"]["var_bins"]

        frac_unc = dict()
        for process in config["process_fractions"]["processes"]:
            frac_unc.update(gd.frac_variation_dict[process])

        fraction = make_1D_fractions(
            fraction_conf=config["process_fractions"],
            variable_info=(var, binning),
            fractions=fractions,
            fraction_name="process_fractions",
            uncertainties=frac_unc,
        )
        corrlib_corrections.append(fraction)

    if "process_fractions_subleading" in config and fractions_subleading is not None:
        var = config["process_fractions_subleading"]["var_dependence"]
        binning = config["process_fractions_subleading"]["var_bins"]

        frac_unc = dict()
        for process in config["process_fractions_subleading"]["processes"]:
            frac_unc.update(gd.frac_variation_dict[process])

        fraction = make_1D_fractions(
            fraction_conf=config["process_fractions_subleading"],
            variable_info=(var, binning),
            fractions=fractions_subleading,
            fraction_name="process_fractions_subleading",
            uncertainties=frac_unc,
        )
        corrlib_corrections.append(fraction)

    cset = cs.CorrectionSet(
        schema_version=2,
        description="Fake factors for tau analysis",
        corrections=corrlib_corrections,
        compound_corrections=None,
    )

    base_filepath = os.path.join(output_path, f"fake_factors_{config['channel']}")
    suffix = "_for_corrections" if for_corrections else ""
    write_json(f"{base_filepath}{suffix}.json", cset)
    write_json(f"{base_filepath}{suffix}.json.gz", cset)


def make_1D_ff(
    process: str,
    process_conf: Dict[str, Union[Dict, List, str]],
    variable_info: Tuple[str, List[float]],
    ff_functions: Dict[str, Dict[str, str]],
) -> cs.Correction:
    """
    Function which produces a correctionlib Correction based on the measured fake factors (including variations) for the case of 1D categories.

    Args:
        process: Name of the process for which the Correction is produced
        process_conf: A dictionary with all the relevant information for the fake factors of a specific "process"
        variable_info: Tuple with information (name and binning) about the variable the fake factors depends on
        ff_functions: Dictionary of fake factor functions as strings, e.g. ff_function[PROCESS][CATEGORY_1][VARIATION]

    Return:
        Correction object from correcionlib
    """
    # get categories from config
    cat_inputs = list(process_conf["split_categories"].keys())
    cat_values = [
        process_conf["split_categories"][cat]
        for cat in process_conf["split_categories"]
    ]

    ff = cs.Correction(
        name=f"{process}_fake_factors",
        description=f"Calculation of the {process} part the for data-driven background estimation (fake factors) for misindentified jets as tau leptons in H->tautau analysis.",
        version=1,
        generic_formulas=None,
        inputs=[
            cs.Variable(
                name=variable_info[0],
                type=gd.variable_type[variable_info[0]],
                description=gd.variable_description[variable_info[0]]
                .replace("#var_min", str(min(variable_info[1])))
                .replace("#var_max", str(max(variable_info[1]))),
            ),
            cs.Variable(
                name=cat_inputs[0],
                type=gd.variable_type[cat_inputs[0]],
                description=gd.variable_description[cat_inputs[0]]
                + ", ".join(cat_values[0]),
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
                    key=process + unc_name,
                    value=cs.Binning(
                        nodetype="binning",
                        input=cat_inputs[0],
                        edges=process_conf["split_categories_binedges"][cat_inputs[0]],
                        content=[
                            cs.Binning(
                                nodetype="binning",
                                input=variable_info[0],
                                **get_edges_and_content(
                                    ff_functions[cat1][unc], variable_info
                                ),
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
                input=cat_inputs[0],
                edges=process_conf["split_categories_binedges"][cat_inputs[0]],
                content=[
                    cs.Binning(
                        nodetype="binning",
                        input=variable_info[0],
                        **get_edges_and_content(
                            ff_functions[cat1]["nominal"], variable_info
                        ),
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


def make_2D_ff(
    process: str,
    process_conf: Dict[str, Union[Dict, List, str]],
    variable_info: Tuple[str, List[float]],
    ff_functions: Dict[str, Dict[str, Dict[str, str]]],
) -> cs.Correction:
    """
    Function which produces a correctionlib Correction based on the measured fake factors (including variations) for the case of 2D categories.

    Args:
        process: Name of the process for which the Correction is produced
        process_conf: A dictionary with all the relevant information for the fake factors of a specific "process"
        variable_info: Tuple with information (name and binning) about the variable the fake factors depends on
        ff_functions: Dictionary of fake factor functions as strings, e.g. ff_function[PROCESS][CATEGORY_1][CATEGORY_2][VARIATION]

    Return:
        Correction object from correcionlib
    """
    # get categories from config
    cat_inputs = list(process_conf["split_categories"].keys())
    cat_values = [
        process_conf["split_categories"][cat]
        for cat in process_conf["split_categories"]
    ]

    ff = cs.Correction(
        name=f"{process}_fake_factors",
        description=f"Calculation of the {process} part the for data-driven background estimation (fake factors) for misindentified jets as tau leptons in H->tautau analysis.",
        version=1,
        generic_formulas=None,
        inputs=[
            cs.Variable(
                name=variable_info[0],
                type=gd.variable_type[variable_info[0]],
                description=gd.variable_description[variable_info[0]]
                .replace("#var_min", str(min(variable_info[1])))
                .replace("#var_max", str(max(variable_info[1]))),
            ),
            cs.Variable(
                name=cat_inputs[0],
                type=gd.variable_type[cat_inputs[0]],
                description=gd.variable_description[cat_inputs[0]]
                + ", ".join(cat_values[0]),
            ),
            cs.Variable(
                name=cat_inputs[1],
                type=gd.variable_type[cat_inputs[1]],
                description=gd.variable_description[cat_inputs[1]]
                + ", ".join(cat_values[1]),
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
                    key=process + unc_name,
                    value=cs.Binning(
                        nodetype="binning",
                        input=cat_inputs[0],
                        edges=process_conf["split_categories_binedges"][cat_inputs[0]],
                        content=[
                            cs.Binning(
                                nodetype="binning",
                                input=cat_inputs[1],
                                edges=process_conf["split_categories_binedges"][
                                    cat_inputs[1]
                                ],
                                content=[
                                    cs.Binning(
                                        nodetype="binning",
                                        input=variable_info[0],
                                        **get_edges_and_content(
                                            ff_functions[cat1][cat2][unc], variable_info
                                        ),
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
                input=cat_inputs[0],
                edges=process_conf["split_categories_binedges"][cat_inputs[0]],
                content=[
                    cs.Binning(
                        nodetype="binning",
                        input=cat_inputs[1],
                        edges=process_conf["split_categories_binedges"][cat_inputs[1]],
                        content=[
                            cs.Binning(
                                nodetype="binning",
                                input=variable_info[0],
                                **get_edges_and_content(
                                    ff_functions[cat1][cat2]["nominal"], variable_info
                                ),
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


def make_1D_fractions(
    fraction_conf: Dict[str, Union[Dict, List, str]],
    variable_info: Tuple[str, List[float]],
    fractions: Dict[str, Dict[str, Dict[str, List[float]]]],
    fraction_name: str,
    uncertainties: Dict[str, str],
) -> cs.Correction:
    """
    Function which produces a correctionlib Correction based on the measured fractions (including variations).

    Args:
        fraction_conf: A dictionary with all the relevant information for the fraction calculation
        variable_info: Tuple with information (name and binning) about the variable the fractions depends on
        fractions: Dictionary of fraction values, e.g. fractions[CATEGORY][VARIATION][PROCESS]
        fraction_name: Name of the calculated fraction, relevant if more than one fraction should be added
        uncertainties: Dictionary of uncertainty names which should be added

    Return:
        Correction object from correcionlib
    """
    # get categories from config
    cat_inputs = list(fraction_conf["split_categories"].keys())
    cat_values = [
        fraction_conf["split_categories"][cat]
        for cat in fraction_conf["split_categories"]
    ]

    frac = cs.Correction(
        name=fraction_name,
        description=f"Calculation of process contributions ({fraction_name}) for the fake factor calculation.",
        version=1,
        generic_formulas=None,
        inputs=[
            cs.Variable(
                name="process",
                type="string",
                description="name of the process",
            ),
            cs.Variable(
                name=variable_info[0],
                type=gd.variable_type[variable_info[0]],
                description=gd.variable_description[variable_info[0]]
                .replace("#var_min", str(min(variable_info[1])))
                .replace("#var_max", str(max(variable_info[1]))),
            ),
            cs.Variable(
                name=cat_inputs[0],
                type=gd.variable_type[cat_inputs[0]],
                description=gd.variable_description[cat_inputs[0]]
                + ", ".join(cat_values[0]),
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
                    key=fraction_name + unc_name,
                    value=cs.Category(
                        nodetype="category",
                        input="process",
                        content=[
                            cs.CategoryItem(
                                key=gd.variable_translator[p],
                                value=cs.Binning(
                                    nodetype="binning",
                                    input=cat_inputs[0],
                                    edges=fraction_conf["split_categories_binedges"][
                                        cat_inputs[0]
                                    ],
                                    content=[
                                        cs.Binning(
                                            nodetype="binning",
                                            input=variable_info[0],
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
                        default=None,
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
                            input=cat_inputs[0],
                            edges=fraction_conf["split_categories_binedges"][
                                cat_inputs[0]
                            ],
                            content=[
                                cs.Binning(
                                    nodetype="binning",
                                    input=variable_info[0],
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
                default=None,
            ),
        ),
    )
    rich.print(frac)

    return frac


def generate_correction_corrlib_json(
    config: Dict[str, Union[str, Dict, List]],
    corrections: Dict[str, Dict[str, Dict[str, np.ndarray]]],
    output_path: str,
    for_DRtoSR: bool = False,
) -> None:
    """
    Function which produces a correctionlib file based on the measured fake factor corrections (including variations).

    Args:
        config: A dictionary with all the relevant information for the fake factors
        corrections: Dictionary of correction arrays, e.g. corrections[PROCESS][CORRECTION][VARIATION]
        output_path: Path where the generated correctionlib files should be stored
        for_DRtoSR: Boolean which decides where the produced correctionlib file is stored, this is needed because additional
                    fake factor corrections are needed to calculate the DR to SR correction and therefore they are stored in
                    a different place (default: False)

    Return:
        None
    """
    corrlib_ff_corrections = list()

    compound_correction_names = dict()
    compound_inputs = dict()

    for process in corrections:
        compound_correction_names[process] = list()
        compound_inputs[process] = list()

        for correction in corrections[process]:
            if "non_closure" in correction and not for_DRtoSR:
                tmp_var = correction.split("non_closure_")[1]
                variable = config["target_processes"][process]["non_closure"][tmp_var]["var_dependence"]
            elif "non_closure" in correction and for_DRtoSR:
                tmp_var = correction.split("non_closure_")[1]
                variable = config["target_processes"][process]["DR_SR"]["non_closure"][tmp_var]["var_dependence"]
            else:
                variable = config["target_processes"][process][correction]["var_dependence"]
            correction_dict = corrections[process][correction]
            corr = make_1D_correction(
                process,
                variable,
                correction,
                correction_dict,
            )
            corrlib_ff_corrections.append(corr)

            if "non_closure" in correction:
                compound_correction_names[process].append(f"{process}_{correction}_correction")
                compound_inputs[process].append(
                    cs.Variable(
                        name=variable,
                        type=gd.variable_type[variable],
                        description=gd.variable_description[variable]
                        .replace("#var_min", str(min(correction_dict["edges"])))
                        .replace("#var_max", str(max(correction_dict["edges"]))),
                    )
                )

    cset = cs.CorrectionSet(
        schema_version=2,
        description="Corrections for fake factors for tau analysis",
        corrections=corrlib_ff_corrections,
        compound_corrections=[
            cs.CompoundCorrection(
                name=f"{process}_compound_correction",
                inputs=compound_inputs[process] + [
                    cs.Variable(
                        name="syst",
                        type="string",
                        description="Uncertainty from the correction measurement.",
                    )
                ],
                output=cs.Variable(
                    name=f"{process}_compound_correction",
                    type="real",
                    description=f"compound correction for {process}",
                ),
                inputs_update=[],
                input_op="*",
                output_op="*",
                stack=compound_correction_names[process],
            )
            for process in compound_correction_names
        ],
    )

    base_filepath = os.path.join(output_path, f"FF_corrections_{config['channel']}")
    suffix = "_for_DRtoSR" if for_DRtoSR else ""
    write_json(f"{base_filepath}{suffix}.json", cset)
    write_json(f"{base_filepath}{suffix}.json.gz", cset)


def make_1D_correction(
    process: str, variable: str, correction_name: str, correction: Dict[str, np.ndarray]
) -> cs.Correction:
    """
    Function which produces a correctionlib Correction based on the measured fake factor corrections (including variations).

    Args:
        process: Name of the process for which the Correction is produced
        variable: Name of the variable the correction depends on
        correction_name: Name of the correction
        correction: Dictionary with the correction information for edge values and nominal, up and down variations

    Return:
        Correction object from correcionlib
    """
    corr = cs.Correction(
        name=f"{process}_{correction_name}_correction",
        description=f"Calculation of the {correction_name} correction for {process} for data-driven background estimation (fake factors) for misindentified jets as tau leptons in H->tautau analysis.",
        version=1,
        generic_formulas=None,
        inputs=[
            cs.Variable(
                name=variable,
                type=gd.variable_type[variable],
                description=gd.variable_description[variable]
                .replace("#var_min", str(min(correction["edges"])))
                .replace("#var_max", str(max(correction["edges"]))),
            ),
            cs.Variable(
                name="syst",
                type="string",
                description="Uncertainty from the correction measurement.",
            ),
        ],
        output=cs.Variable(
            name=f"{process}_{correction_name}_correction",
            type="real",
            description=f"{correction_name} correction for {process}",
        ),
        data=cs.Category(
            nodetype="category",
            input="syst",
            content=[
                cs.CategoryItem(
                    key=process + f"_{correction_name}_CorrUp",
                    value=cs.Binning(
                        nodetype="binning",
                        input=variable,
                        edges=list(correction["edges"]),
                        content=list(correction["up"]),
                        flow="clamp",
                    ),
                ),
                cs.CategoryItem(
                    key=process + f"_{correction_name}_CorrDown",
                    value=cs.Binning(
                        nodetype="binning",
                        input=variable,
                        edges=list(correction["edges"]),
                        content=list(correction["down"]),
                        flow="clamp",
                    ),
                ),
            ],
            default=cs.Binning(
                nodetype="binning",
                input=variable,
                edges=list(correction["edges"]),
                content=list(correction["nominal"]),
                flow="clamp",
            ),
        ),
    )
    rich.print(corr)

    return corr
