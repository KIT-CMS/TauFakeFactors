import gzip
import json
import os
from itertools import product
from typing import Any, Dict, List, Tuple, Union

import correctionlib.schemav2 as cs
import numpy as np
import rich

import configs.general_definitions as gd
import helper.ff_functions as ff_func


def replace_varmin_varmax(string: str, obj: Dict[str, Any]) -> str:
    """
    This function replaces the placeholders for the variable minimum and maximum with the actual values.

    Args:
        string: String where the placeholders should be replaced
        obj: Dictionary with the actual values for the placeholders

    Return:
        String with the replaced placeholders

    """
    # 1D FF, 1D fractions
    if isinstance(obj, dict) and all(isinstance(it, list) for it in obj.values()):
        string = string.replace(
            "#var_min and #var_max GeV",
            '; '.join(f"({k.replace('#', '')}): ({min(obj[k])}, {max(obj[k])}) GeV" for k in obj),
        )
    # 2D FF
    elif isinstance(obj, dict) and all(
        isinstance(subdict, dict) and all(
            isinstance(lst, list) for lst in subdict.values()
        )
        for subdict in obj.values()
    ):
        string = string.replace(
            "#var_min and #var_max GeV",
            '; '.join(
                f"({outer.replace('#', '')} && {inner.replace('#', '')}): ({min(lst)}, {max(lst)}) GeV"
                for outer, subdict in obj.items()
                for inner, lst in subdict.items()
            )
        )
    # 2D corrections v1
    elif isinstance(obj, dict) and all("edges" in subdict for subdict in obj.values()):
        string = string.replace(
            "#var_min and #var_max GeV",
            '; '.join(f"({k.replace('#', '')}): ({min(obj[k]['edges'])}, {max(obj[k]['edges'])}) GeV" for k in obj),
        )
    # 2D corrections v2 (TODO: might be generalized)
    elif isinstance(obj, dict) and all(isinstance(it, list) for it in obj.values()):
        string = string.replace(
            "#var_min and #var_max GeV",
            '; '.join(f"({k.replace('#', '')}): ({min(obj[k])}, {max(obj[k])}) GeV" for k in obj),
        )
    # 1D corrections categories string
    elif isinstance(obj, dict) and "edges" in obj:
        if "GeV" not in string and "#var_min" not in string and "#var_max" not in string:
            string += str(list(obj["edges"])).replace("[", "(").replace("]", ")")
    # 1D corrections categories without splitting
    elif isinstance(obj, list):
        string = string.replace(
            "#var_min and #var_max GeV",
            f"({min(obj)}, {max(obj)}) GeV",
        )

    return string


def get_edges_and_content(
    item: Union[str, Any],
    variable_name: str,
    variable_bins: list,
):
    if isinstance(item, str):
        return {
            "edges": [
                min(variable_bins),
                max(variable_bins),
            ],
            "content": [
                cs.Formula(
                    nodetype="formula",
                    variables=[variable_name],
                    parser="TFormula",
                    expression=item,
                    parameters=None,
                ),
            ],
        }
    return {
        "edges": variable_bins,
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
    output_path: Union[str, None] = None,
    for_corrections: bool = False,
) -> cs.CorrectionSet:
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
            proc_conf = config["target_processes"][process]
            var = proc_conf["var_dependence"]
            binning = ff_func.SplitQuantities(proc_conf).to_dict("var_bins")

            if len(proc_conf["split_categories"]) == 1:
                process_ff = make_1D_ff(
                    process=process,
                    process_conf=proc_conf,
                    variable_info=(var, binning),
                    ff_functions=ff_functions[process],
                )
            elif len(proc_conf["split_categories"]) == 2:
                process_ff = make_2D_ff(
                    process=process,
                    process_conf=proc_conf,
                    variable_info=(var, binning),
                    ff_functions=ff_functions[process],
                )
            corrlib_corrections.append(process_ff)

    if "process_fractions" in config and fractions is not None:
        frac_conf = config["process_fractions"]
        var = frac_conf["var_dependence"]
        binning = ff_func.SplitQuantities(frac_conf).to_dict("var_bins")

        frac_unc = dict()
        for proc in frac_conf["processes"]:
            frac_unc.update(gd.frac_variation_dict[proc])

        fraction = make_1D_fractions(
            fraction_conf=frac_conf,
            variable_info=(var, binning),
            fractions=fractions,
            fraction_name="process_fractions",
            uncertainties=frac_unc,
        )
        corrlib_corrections.append(fraction)

    if "process_fractions_subleading" in config and fractions_subleading is not None:
        frac_conf = config["process_fractions_subleading"]
        var = frac_conf["var_dependence"]
        binning = ff_func.SplitQuantities(frac_conf).to_dict("var_bins")

        frac_unc = dict()
        for proc in frac_conf["processes"]:
            frac_unc.update(gd.frac_variation_dict[proc])

        fraction = make_1D_fractions(
            fraction_conf=frac_conf,
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

    if output_path is not None:
        base_filepath = os.path.join(output_path, f"fake_factors_{config['channel']}")
        suffix = "_for_corrections" if for_corrections else ""
        write_json(f"{base_filepath}{suffix}.json", cset)
        write_json(f"{base_filepath}{suffix}.json.gz", cset)
    return cset


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
    cat_values = [process_conf["split_categories"][cat] for cat in process_conf["split_categories"]]

    ff = cs.Correction(
        name=f"{process}_fake_factors",
        description=f"Calculation of the {process} part the for data-driven background estimation (fake factors) for misindentified jets as tau leptons in H->tautau analysis.",
        version=1,
        generic_formulas=None,
        inputs=[
            cs.Variable(
                name=variable_info[0],
                type=gd.variable_type[variable_info[0]],
                description=replace_varmin_varmax(
                    gd.variable_description[variable_info[0]],
                    variable_info[1],
                )
            ),
            cs.Variable(
                name=cat_inputs[0],
                type=gd.variable_type[cat_inputs[0]],
                description=f"{gd.variable_description[cat_inputs[0]]} {', '.join(cat_values[0])}",
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
                                    ff_functions[cat1][unc],
                                    variable_info[0],
                                    variable_info[1][cat1],
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
                            ff_functions[cat1]["nominal"],
                            variable_info[0],
                            variable_info[1][cat1],
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
    cat_values_conf = process_conf["split_categories"]
    binedges_conf = process_conf["split_categories_binedges"]
    cat2_values_for_desc = []
    cat2_conf = cat_values_conf[cat_inputs[1]]
    if isinstance(cat2_conf, dict):  # New format
        for cat_list in cat2_conf.values():
            cat2_values_for_desc.extend(cat_list)
        cat2_values_for_desc = list(dict.fromkeys(cat2_values_for_desc))
    else:  # Old format
        cat2_values_for_desc = cat2_conf

    ff = cs.Correction(
        name=f"{process}_fake_factors",
        description=f"Calculation of the {process} part the for data-driven background estimation (fake factors) for misindentified jets as tau leptons in H->tautau analysis.",
        version=1,
        generic_formulas=None,
        inputs=[
            cs.Variable(
                name=variable_info[0],
                type=gd.variable_type[variable_info[0]],
                description=replace_varmin_varmax(
                    gd.variable_description[variable_info[0]],
                    variable_info[1],
                ),
            ),
            cs.Variable(
                name=cat_inputs[0],
                type=gd.variable_type[cat_inputs[0]],
                description=f"{gd.variable_description[cat_inputs[0]]} {', '.join(cat_values_conf[cat_inputs[0]])}",
            ),
            cs.Variable(
                name=cat_inputs[1],
                type=gd.variable_type[cat_inputs[1]],
                description=f"{gd.variable_description[cat_inputs[1]]} {', '.join(cat2_values_for_desc)}",
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
                        edges=binedges_conf[cat_inputs[0]],
                        content=[
                            cs.Binning(
                                nodetype="binning",
                                input=cat_inputs[1],
                                edges=(
                                    binedges_conf[cat_inputs[1]][cat1.split("#")[-1]]
                                    if isinstance(binedges_conf[cat_inputs[1]], dict)
                                    else binedges_conf[cat_inputs[1]]
                                ),
                                content=[
                                    cs.Binning(
                                        nodetype="binning",
                                        input=variable_info[0],
                                        **get_edges_and_content(
                                            ff_functions[cat1][cat2][unc],
                                            variable_info[0],
                                            variable_info[1][cat1][cat2],
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
                edges=binedges_conf[cat_inputs[0]],
                content=[
                    cs.Binning(
                        nodetype="binning",
                        input=cat_inputs[1],
                        edges=(
                            binedges_conf[cat_inputs[1]][cat1.split("#")[-1]]
                            if isinstance(binedges_conf[cat_inputs[1]], dict)
                            else binedges_conf[cat_inputs[1]]
                        ),
                        content=[
                            cs.Binning(
                                nodetype="binning",
                                input=variable_info[0],
                                **get_edges_and_content(
                                    ff_functions[cat1][cat2]["nominal"],
                                    variable_info[0],
                                    variable_info[1][cat1][cat2],
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
    cat_values = [fraction_conf["split_categories"][cat] for cat in fraction_conf["split_categories"]]
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
                description=replace_varmin_varmax(
                    gd.variable_description[variable_info[0]],
                    variable_info[1],
                )
            ),
            cs.Variable(
                name=cat_inputs[0],
                type=gd.variable_type[cat_inputs[0]],
                description=f"{gd.variable_description[cat_inputs[0]]} {', '.join(cat_values[0])}",
            ),
            cs.Variable(
                name="syst",
                type="string",
                description="Uncertainties from the +/-7% variation of the process fractions.",
            ),
        ],
        output=cs.Variable(
            name="fraction",
            type="real",
            description="process fraction",
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
                                    edges=fraction_conf["split_categories_binedges"][cat_inputs[0]],
                                    content=[
                                        cs.Binning(
                                            nodetype="binning",
                                            input=variable_info[0],
                                            edges=variable_info[1][cat],
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
                            edges=fraction_conf["split_categories_binedges"][cat_inputs[0]],
                            content=[
                                cs.Binning(
                                    nodetype="binning",
                                    input=variable_info[0],
                                    edges=variable_info[1][cat],
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


def generate_correction_corrlib(
    config: Dict[str, Union[str, Dict, List]],
    corrections: Dict[str, Dict[str, Dict[str, np.ndarray]]],
    output_path: Union[str, None] = None,
    for_DRtoSR: bool = False,
) -> cs.CorrectionSet:
    """
    Function which produces a correctionlib file based on the measured fake factor corrections (including variations).

    Args:
        config: A dictionary with all the relevant information for the fake factors
        corrections: Dictionary of correction arrays, e.g. corrections[PROCESS][CORRECTION][VARIATION]
        output_path: Path where the generated correctionlib files should be stored if provided
        for_DRtoSR: Boolean which decides where the produced correctionlib file is stored, this is needed because additional
                    fake factor corrections are needed to calculate the DR to SR correction and therefore they are stored in
                    a different place (default: False)

    Return:
        None
    """
    corrlib_ff_corrections = []

    compound_correction_names = {}
    compound_inputs = {}
    compound_inputs_splits = {}

    class UniqueList(list):
        def unique_append(self, item):
            if item not in self:
                self.append(item)

    for process in corrections:
        compound_correction_names[process] = []
        compound_inputs[process] = []
        compound_inputs_splits[process] = UniqueList()

        for correction in corrections[process]:
            if "non_closure" in correction and not for_DRtoSR:
                tmp_var = correction.split("non_closure_")[1]
                variable = config["target_processes"][process]["non_closure"][tmp_var]["var_dependence"]
                correction_conf = config["target_processes"][process]["non_closure"][tmp_var]
            elif "non_closure" in correction and for_DRtoSR:
                tmp_var = correction.split("non_closure_")[1]
                variable = config["target_processes"][process]["DR_SR"]["non_closure"][tmp_var]["var_dependence"]
                correction_conf = config["target_processes"][process]["DR_SR"]["non_closure"][tmp_var]
            else:
                variable = config["target_processes"][process][correction]["var_dependence"]
                correction_conf = config["target_processes"][process][correction]

            correction_dict = corrections[process][correction]
            is_2D = "split_categories" in correction_conf

            binning = correction_conf["var_bins"]

            correction_variations = correction_conf.get("correction_variations", gd.default_correction_variations)
            correction_variations = ["".join(it) for it in product(correction_variations, ["Up", "Down"])]

            if is_2D:
                _obj = ff_func.SplitQuantities(correction_conf)
                split_variables, binning = _obj.split_variables, _obj.to_dict("var_bins")
                corr = make_2D_correction(
                    process=process,
                    corr_conf=correction_conf,
                    correction_name=correction,
                    corrections=correction_dict,
                    variable=variable,
                    variations=correction_variations,
                )
            else:
                corr = make_1D_correction(
                    process=process,
                    variable=variable,
                    correction_name=correction,
                    correction=correction_dict,
                    variations=correction_variations,
                )

            corrlib_ff_corrections.append(corr)

            if "non_closure" in correction:
                compound_correction_names[process].append(f"{process}_{correction}_correction")
                if is_2D:
                    compound_inputs_splits[process].unique_append(
                        cs.Variable(
                            name=split_variables[0],
                            type=gd.variable_type[split_variables[0]],
                            description=gd.variable_description[split_variables[0]],
                        ),
                    )
                compound_inputs[process].append(
                    cs.Variable(
                        name=variable,
                        type=gd.variable_type[variable],
                        description=replace_varmin_varmax(
                            gd.variable_description[variable],
                            binning,
                        ),
                    )
                )
        compound_inputs[process].extend(compound_inputs_splits[process])

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

    if output_path is not None:
        base_filepath = os.path.join(output_path, f"FF_corrections_{config['channel']}")
        suffix = "_for_DRtoSR" if for_DRtoSR else ""
        write_json(f"{base_filepath}{suffix}.json", cset)
        write_json(f"{base_filepath}{suffix}.json.gz", cset)

    return cset


def make_1D_correction(
    process: str,
    variable: str,
    correction_name: str,
    correction: Dict[str, np.ndarray],
    variations: Tuple[str],
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
                description=replace_varmin_varmax(
                    gd.variable_description[variable],
                    correction,
                ),
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
                    key=f"{process}_{correction_name}_Corr{variation}",
                    value=cs.Binning(
                        nodetype="binning",
                        input=variable,
                        edges=list(correction["edges"]),
                        content=list(correction[variation]),
                        flow="clamp",
                    ),
                )
                for variation in variations
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


def make_2D_correction(
    process: str,
    corr_conf: Dict[str, Union[Dict, List, str]],
    corrections: Dict[str, Dict[str, np.ndarray]],
    correction_name: str,
    variable: str,
    variations: Tuple[str],
) -> cs.Correction:
    cat_inputs = list(corr_conf["split_categories"].keys())
    cat_values = [corr_conf["split_categories"][cat] for cat in corr_conf["split_categories"]]

    corr = cs.Correction(
        name=f"{process}_{correction_name}_correction",
        description=f"Calculation of the {correction_name} correction for {process} for data-driven background estimation (fake factors) for misindentified jets as tau leptons in H->tautau analysis.",
        version=1,
        generic_formulas=None,
        inputs=[
            cs.Variable(
                name=variable,
                type=gd.variable_type[variable],
                description=replace_varmin_varmax(
                    gd.variable_description[variable],
                    corrections,
                ),
            ),
            cs.Variable(
                name=cat_inputs[0],
                type=gd.variable_type[cat_inputs[0]],
                description=f"{gd.variable_description[cat_inputs[0]]} {', '.join(cat_values[0])}",
            ),
            cs.Variable(
                name="syst",
                type="string",
                description="Uncertainty from the correction measurement",
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
                    key=f"{process}_{correction_name}_Corr{variation}",
                    value=cs.Binning(
                        nodetype="binning",
                        input=cat_inputs[0],
                        edges=corr_conf["split_categories_binedges"][cat_inputs[0]],
                        content=[
                            cs.Binning(
                                nodetype="binning",
                                input=variable,
                                edges=list(corrections[cat1]["edges"]),
                                content=list(corrections[cat1][variation]),
                                flow="clamp",
                            )
                            for cat1 in corrections
                        ],
                        flow="clamp",
                    ),
                )
                for variation in variations
            ],
            default=cs.Binning(
                nodetype="binning",
                input=cat_inputs[0],
                edges=corr_conf["split_categories_binedges"][cat_inputs[0]],
                content=[
                    cs.Binning(
                        nodetype="binning",
                        input=variable,
                        edges=list(corrections[cat1]["edges"]),
                        content=list(corrections[cat1]["nominal"]),
                        flow="clamp",
                    )
                    for cat1 in corrections
                ],
                flow="clamp",
            ),
        ),
    )
    rich.print(corr)

    return corr
