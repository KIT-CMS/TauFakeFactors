import argparse
from collections import OrderedDict
from itertools import pairwise
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd
import uproot
from rich.console import Console
from rich.table import Table
from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap, CommentedSeq
from ruamel.yaml.scalarstring import DoubleQuotedScalarString

EQUIPOPULATED_BINNING_OPTIONS_KEY = "equipopulated_binning_options"

parser = argparse.ArgumentParser(description="Equipopulated binning adjustments of configuration.")
parser.add_argument("--config", type=str, required=True, help="Path to the YAML configuration.")
parser.add_argument("--cuts-config", type=str, default=None, help="Path to a YAML correction containing cuts for processes (only for corrections).")
parser.add_argument("--processes", type=str, nargs="+", required=True, help="List of processes (QCD Wjets ttbar process_fractions) to adjust binning for.")
parser.add_argument("--cut-region", type=str, default="SRlike", help="Cut region to apply (e.g., SRlike, ARlike), process_fractions will always use AR_cuts.")
parser.add_argument("--dry-run", action="store_true", help="If set, do not modify the config file.")


def to_commented_map(items: dict) -> CommentedMap:
    """
    Recursively converts a python dict to a ruamel.yaml CommentedMap here leaf lists are 
    formatted using flow style and keys are always quoted.

    Args:
        items (dict): The dictionary to convert.

    Returns:
        CommentedMap: A ruamel.yaml CommentedMap with flow style lists and quoted keys
    """
    if not isinstance(items, dict):
        return items
    commented_map = CommentedMap()
    for key, value in items.items():
        if any(op in str(key) for op in ["<=", ">=", "<", ">", "==", "#"]):
            key = DoubleQuotedScalarString(key)

        if isinstance(value, dict):
            commented_map[key] = to_commented_map(value)
        elif isinstance(value, list):
            is_leaf_list = not any(isinstance(item, (dict, list)) for item in value)

            if is_leaf_list:  # create a CommentedSeq and set flow style.
                comment_sequence = CommentedSeq(value)
                comment_sequence.fa.set_flow_style()
                commented_map[key] = comment_sequence
            else:  # otherwise: recursively process each item.
                commented_map[key] = [to_commented_map(item) for item in value]
        else:
            commented_map[key] = value
    return commented_map


def _equipopulated_binned_variable(item: Union[np.ndarray, pd.Series], n_bins: int) -> np.ndarray:
    """
    Helper function to calculate equipopulated bin edges.

    Args:
        item (pd.Series or np.ndarray): The data to bin.
        n_bins (int): The number of bins to create.

    Returns:
        list: A list of bin edges.
    """
    if len(item) == 0 or n_bins == 0:
        return []
    return np.quantile(item, np.linspace(0, 1, n_bins + 1))


def get_binning(
    df: pd.DataFrame,
    binning_config: dict,
    process_cuts: str,
    category_splits: dict,
    original_category_splits: dict,
    var_dependence: str,
    prepend_space: int = 0,
    parent_cat_keys: Union[list, None] = None,
) -> tuple:
    """
    Recursively calculates equipopulated binning for a variable based on category splits.

    Args:
        df (pd.DataFrame): The input DataFrame containing the data.
        binning_config (dict): Configuration for binning, including n_bins and variable_config.
        process_cuts (str): Cuts to apply to the DataFrame.
        category_splits (dict): Dictionary defining how to split categories.
        original_category_splits (dict): Original category splits for reference.
        var_dependence (str): The variable to be binned.
        prepend_space (int): Number of spaces to prepend for formatting.
        parent_cat_keys (list, optional): List of parent category keys for recursive calls.

    Returns:
        tuple: A tuple containing:
    """
    parent_cat_keys = parent_cat_keys or []
    current_split_var = next(iter(category_splits))

    table = Table(
        show_header=True,
        header_style="bold magenta",
        box=None,
        show_lines=False,
        title=f"Binning for {current_split_var}",
    )
    table.add_column("Category", style="dim", width=30)
    table.add_column("Total Events", justify="right")
    table.add_column("Bins")
    table.add_column("Events/Bin", justify="left")

    generated_edges = {}
    generated_categories = {}
    remaining_splits = {k: v for k, v in category_splits.items() if k != current_split_var}

    if isinstance(category_splits[current_split_var], (list, CommentedSeq)) and category_splits[current_split_var]:
        categories = category_splits[current_split_var]
        generated_categories[current_split_var] = categories
    elif all(
        [
            current_split_var in binning_config.get("n_bins", {}),
            isinstance(binning_config["n_bins"][current_split_var], int),
        ]
    ):  # Get categories for current split level
        n_split_bins = binning_config["n_bins"][current_split_var]
        var_conf = binning_config["variable_config"][current_split_var]
        current_df = df.query(process_cuts)
        split_edges = _equipopulated_binned_variable(current_df[current_split_var], n_split_bins)

        if len(split_edges) > 0:
            split_edges[0] = var_conf.get("min", split_edges[0])
            split_edges[-1] = var_conf.get("max", split_edges[-1])

        rounding = var_conf.get("rounding", 2)
        generated_edges[current_split_var] = np.round(split_edges, rounding).tolist()

        categories = []
        for i in range(len(split_edges) - 1):
            lower = round(split_edges[i], rounding)
            upper = round(split_edges[i + 1], rounding)

            if i == 0:  # syntax setting
                categories.append(f"<={upper}")
            elif i == len(split_edges) - 2:
                categories.append(f">{lower}")
            else:
                categories.append(f">{lower}#&&#<={upper}")

        generated_categories[current_split_var] = categories
    else:
        # Categories are explicitly provided in the config: i.e. njets or tau_decaymode_2
        # where the discrete categories are defined.
        if isinstance(category_splits[current_split_var], (list, CommentedSeq)):
            categories = category_splits[current_split_var]
            generated_categories[current_split_var] = categories
        else:  # Categories need to be generated based on n_bins.
            temp_config = binning_config.get("n_bins", {})  # Config for current split variable setting
            for key in parent_cat_keys:  # parent_cat_keys: path to the current category, e.g., ['==0']
                if isinstance(temp_config, dict) and key in temp_config:
                    temp_config = temp_config[key]
                else:
                    break  # Path does not exist (global setting for variable)

            if isinstance(temp_config, dict) and current_split_var in temp_config:
                n_bins = temp_config[current_split_var]
                var_conf = binning_config.get("variable_config", {}).get(current_split_var, {})
            else:
                raise ValueError(
                    f"Unable to determine categories for '{current_split_var}'"
                    "neither trough equipopulated binning option nor explicit categories."
                )

            min_val, max_val = var_conf.get("min"), var_conf.get("max")
            range_val = max_val - min_val if min_val is not None and max_val is not None else 1
            categories = [
                f">{min_val + i * range_val / n_bins}#&&#<={min_val + (i + 1) * range_val / n_bins}"
                for i in range(n_bins)
            ]  # Default categories if none are provided            
            generated_categories[current_split_var] = categories

    output_bins = OrderedDict()
    table = Table(show_header=True, header_style="bold magenta")

    category_header = f"Category: {current_split_var}"
    if parent_cat_keys:
        original_split_vars = list(original_category_splits.keys())
        parent_var_name = original_split_vars[len(parent_cat_keys) - 1]
        parent_key = parent_cat_keys[-1]
        category_header += f" (in {parent_var_name} {parent_key})"

    table.add_column(category_header)
    table.add_column("Num Events")
    table.add_column(f"Bin Edges ({var_dependence})")
    table.add_column("Num Events/Bin")

    for cat_key in categories:
        query_cat_cut = (
            f"{current_split_var} {cat_key}"
            if not any(op in cat_key for op in ["<", ">", "=="])
            else cat_key
        )
        if "#&&#" in query_cat_cut:
            sub_cuts = query_cat_cut.split("#&&#")
            query_cat_cut = f"({current_split_var} {sub_cuts[0]} & {current_split_var} {sub_cuts[1]})"
        elif "#||#" in query_cat_cut:
            sub_cuts = query_cat_cut.split("#||#")
            query_cat_cut = f"({current_split_var} {sub_cuts[0]} | {current_split_var} {sub_cuts[1]})"
        else:
            query_cat_cut = f"({current_split_var} {cat_key})"

        combined_cut = f"{process_cuts} & {query_cat_cut}"
        sub_df = df.query(combined_cut)
        current_cat_keys = parent_cat_keys + [cat_key]

        if remaining_splits:
            table.add_row(f"{' ' * prepend_space}{cat_key}", str(len(sub_df)), "", "")
            nested_bins, _, nested_edges, nested_cats = get_binning(
                df=df,
                binning_config=binning_config,
                process_cuts=combined_cut,
                category_splits=remaining_splits,
                original_category_splits=original_category_splits,
                var_dependence=var_dependence,
                prepend_space=prepend_space + 2,
                parent_cat_keys=current_cat_keys,
            )
            output_bins[cat_key] = nested_bins
            for var, edges in nested_edges.items():
                if var not in generated_edges:
                    generated_edges[var] = {}
                generated_edges[var][cat_key] = edges
            for var, cats in nested_cats.items():
                if var not in generated_categories:
                    generated_categories[var] = {}
                generated_categories[var][cat_key] = cats
        else:  # final level of splitting
            target_variable = var_dependence
            temp_config, n_bins = binning_config.get("n_bins", {}), None
            original_split_vars = list(original_category_splits.keys())

            cat_keys_map = dict(zip(original_split_vars, current_cat_keys))
            if (search_order := [var for var in original_split_vars if var in temp_config]):
                temp_config = temp_config[search_order[0]]
                for split_var in search_order:
                    if isinstance(temp_config, int):
                        n_bins = temp_config
                        break
                    if cat_keys_map[split_var] in temp_config:
                        temp_config = temp_config[cat_keys_map[split_var]]
                    else:
                        break

            if isinstance(temp_config, int):
                n_bins = temp_config

            if n_bins is None:  # Fallback for flat n_bins definition
                for split_var, key in zip(original_split_vars, current_cat_keys):
                    if not (isinstance(temp_config, dict) and split_var in temp_config):
                        break
                    temp_config = temp_config[split_var]
                    if isinstance(temp_config, int):
                        n_bins = temp_config
                        break
                    if key in temp_config:
                        temp_config = temp_config[key]
                    else:
                        if isinstance(temp_config, int):
                            n_bins = temp_config
                        break  # Stop searching deeper
                if isinstance(temp_config, int):
                    n_bins = temp_config

            if n_bins and n_bins > 0:
                var_conf = binning_config.get("variable_config", {}).get(target_variable, {})

                filtered_df = sub_df
                if (min_val := var_conf.get("min")) is not None:
                    filtered_df = filtered_df.query(f"{target_variable} >= {min_val}")
                if (max_val := var_conf.get("max")) is not None:
                    filtered_df = filtered_df.query(f"{target_variable} <= {max_val}")

                bins = _equipopulated_binned_variable(filtered_df[target_variable], n_bins)
                if len(bins) > 0:
                    bins[0] = var_conf.get("min", bins[0])
                    bins[-1] = var_conf.get("max", bins[-1])

                if (rounding := var_conf.get("rounding", 2)) is not None:
                    bins = np.round(bins, rounding).tolist()

                for option in ["add_left", "add_right"]:
                    if option not in binning_config:
                        continue

                    temp_config = binning_config[option]
                    split_vars = list(original_category_splits.keys())[:len(parent_cat_keys)]

                    path_found = True
                    for split_var, key in zip(split_vars, parent_cat_keys):
                        if not (isinstance(temp_config, dict) and split_var in temp_config):
                            path_found = False
                            break
                        temp_config = temp_config[split_var]
                        if not (isinstance(temp_config, dict) and key in temp_config):
                            break
                        temp_config = temp_config[key]

                    if not (path_found and isinstance(temp_config, (int, float, list, tuple))):
                        continue  # no valid value

                    value_to_add = temp_config

                    if not isinstance(value_to_add, (list, tuple)):
                        value_to_add = [value_to_add]

                    bins = value_to_add + bins if option == "add_left" else bins + value_to_add

                events_per_bin = []
                if len(bins) > 1:
                    pairs = list(pairwise(bins))
                    for lower_edge, upper_edge in pairs[:-1]:
                        count = sub_df.query(f"{lower_edge} <= {target_variable} < {upper_edge}").shape[0]
                        events_per_bin.append(count)

                    lower_edge, upper_edge = pairs[-1]  # exclusive lower edge, inclusive upper edge
                    count = sub_df.query(f"{lower_edge} <= {target_variable} <= {upper_edge}").shape[0]
                    events_per_bin.append(count)

                table.add_row(f"{' ' * prepend_space}{cat_key}", str(len(sub_df)), str(bins), str(events_per_bin))
                output_bins[cat_key] = bins
            else:
                table.add_row(f"{' ' * prepend_space}{cat_key}", str(len(sub_df)), "No binning", "")

    console.print(table)

    return output_bins, parent_cat_keys, generated_edges, generated_categories


if __name__ == "__main__":
    args = parser.parse_args()
    args.config = Path(args.config).resolve()
    console = Console()
    console.print(f"Loading configuration from [cyan]{args.config}[/cyan]")
    console.print(f"Processes to adjust binning for: [cyan]{args.processes}[/cyan]")
    console.print(f"Cut region: [cyan]{args.cut_region}[/cyan]")

    yaml = YAML()
    yaml.preserve_quotes = True
    yaml.indent(mapping=2, sequence=4, offset=2)
    yaml.width = 4096

    cuts_config = None
    if args.cuts_config:
        console.print(f"Loading cuts from [cyan]{args.cuts_config}[/cyan]")
        with open(Path(args.cuts_config).resolve(), "r") as f:
            cuts_config = yaml.load(f)

    with open((args.config.parent / "common_settings.yaml").resolve(), "r") as f:
        common_settings = yaml.load(f)

    with open(args.config, "r") as f:
        config = yaml.load(f)

    base_directory, directory, era, channel, tree, file = (
        Path(common_settings["output_path"]),
        "preselection",
        common_settings["era"],
        config["channel"],
        common_settings["tree"],
        "data.root",
    )
    root_file_path = base_directory / directory / era / channel / file
    console.print(f"Loading data from [cyan]{root_file_path}[/cyan]")
    dataframe = uproot.open(root_file_path)[tree].arrays(library="pd")

    if len((bool_cols := dataframe.select_dtypes(include=['bool']).columns)) > 0:
        dataframe[bool_cols] = dataframe[bool_cols].astype(int)
        console.print(f"Converted {len(bool_cols)} boolean columns to integers.")

    for process in args.processes:
        console.print(f"\n[bold green]Processing: {process}[/bold green]")

        process_config, cuts = None, None

        if process == "process_fractions":
            if "process_fractions" not in config:
                console.print("[red]Block 'process_fractions' not found in config. Skipping.[/red]")
                continue
            process_config = config["process_fractions"]
            cuts = " & ".join(f"({c})" for c in process_config["AR_cuts"].values())
            console.print("Using AR cuts for process_fractions.")

        elif process in config.get("target_processes", {}):
            process_config = config["target_processes"][process]

            cuts_source = cuts_config["target_processes"][process] if cuts_config else process_config
            cuts = " & ".join(f"({c})" for c in cuts_source[f"{args.cut_region}_cuts"].values())

        else:
            console.print(f"[red]Process '{process}' not found in target_processes or as process_fractions. Skipping.[/red]")
            continue

        if EQUIPOPULATED_BINNING_OPTIONS_KEY in process_config:
            binning_config = process_config[EQUIPOPULATED_BINNING_OPTIONS_KEY]
            temp_cuts = cuts.replace('&&', '&').replace('!', '~').replace('||', '|')

            console.print(f"Variable to bin: [bold magenta]{process_config['var_dependence']}[/bold magenta]")
            console.print(f"Base cuts: [dim]{temp_cuts}[/dim]")

            new_bins, _, new_edges, new_cats = get_binning(
                df=dataframe,
                binning_config=binning_config,
                process_cuts=temp_cuts,
                category_splits=process_config["split_categories"],
                original_category_splits=process_config["split_categories"],
                var_dependence=process_config["var_dependence"],
            )

            process_config["var_bins"] = to_commented_map(new_bins)
            process_config["split_categories"].update(to_commented_map(new_cats))
            if "split_categories_binedges" not in process_config or not process_config["split_categories_binedges"]:
                process_config["split_categories_binedges"] = CommentedMap()
            process_config["split_categories_binedges"].update(to_commented_map(new_edges))

        def get_merged_cuts(
            var_config: Union[dict, CommentedMap],
            parent_config: Union[dict, CommentedMap, None] = None,
        ) -> str:
            """
            Function to merge cuts for corrections, taking into account parent configurations.
            var_config represents the current correction block configuration,
            parent_config is the parent correction block configuration (if any), i.e. DR_SR 
            if provided.

            Args:
                var_config (dict): Configuration for the current correction block.
                parent_config (dict, optional): Parent configuration to override cuts if available.

            Returns:
                str: Merged cuts as a string.
            """

            merged_cuts = {}
            if cuts_config:  # start from external file if provided
                cuts_source = cuts_config["target_processes"][process]
                merged_cuts = cuts_source.get(f"{args.cut_region}_cuts", {}).copy()

            if parent_config:  # overwrite with parent cuts if available i.e. DR_SR
                parent_cuts = parent_config.get(f"{args.cut_region}_cuts", {})
                merged_cuts.update(parent_cuts)

            local_cuts = var_config.get(f"{args.cut_region}_cuts", {})
            merged_cuts.update(local_cuts)

            return " & ".join(f"({c})" for c in merged_cuts.values())

        if "non_closure" in process_config:  # top level non_closure corrections
            for var, var_config in process_config["non_closure"].items():
                if EQUIPOPULATED_BINNING_OPTIONS_KEY not in var_config:
                    continue

                final_cuts = get_merged_cuts(var_config)
                temp_cuts = final_cuts.replace('&&', '&').replace('!', '~').replace('||', '|')

                console.print(f"  Processing correction: [bold magenta]non_closure/{var}[/bold magenta]")
                console.print(f"    Variable to bin: [bold magenta]{var_config['var_dependence']}[/bold magenta]")
                console.print(f"    Base cuts: [dim]{temp_cuts}[/dim]")

                new_bins, _, new_edges, new_cats = get_binning(
                    df=dataframe,
                    binning_config=var_config[EQUIPOPULATED_BINNING_OPTIONS_KEY],
                    process_cuts=temp_cuts,
                    category_splits=var_config["split_categories"],
                    original_category_splits=var_config["split_categories"],
                    var_dependence=var_config["var_dependence"],
                )

                var_config["var_bins"] = to_commented_map(new_bins)
                var_config["split_categories"].update(to_commented_map(new_cats))
                if "split_categories_binedges" not in var_config or not var_config["split_categories_binedges"]:
                    var_config["split_categories_binedges"] = CommentedMap()
                var_config["split_categories_binedges"].update(to_commented_map(new_edges))

        if "DR_SR" in process_config:
            if EQUIPOPULATED_BINNING_OPTIONS_KEY not in process_config["DR_SR"]:
                console.print("[red]DR_SR does not have equipopulated_binning_options. Skipping.[/red]")
                continue

            final_cuts = get_merged_cuts(process_config["DR_SR"])
            temp_cuts = final_cuts.replace('&&', '&').replace('!', '~').replace('||', '|')

            console.print(f"  Processing correction: [bold magenta]DR_SR[/bold magenta]")
            console.print(f"    Variable to bin: [bold magenta]{process_config['DR_SR']['var_dependence']}[/bold magenta]")
            console.print(f"    Base cuts: [dim]{temp_cuts}[/dim]")

            new_bins, _, new_edges, new_cats = get_binning(
                df=dataframe,
                binning_config=process_config["DR_SR"][EQUIPOPULATED_BINNING_OPTIONS_KEY],
                process_cuts=temp_cuts,
                category_splits=process_config["DR_SR"]["split_categories"],
                original_category_splits=process_config["DR_SR"]["split_categories"],
                var_dependence=process_config["DR_SR"]["var_dependence"],
            )

            process_config["DR_SR"]["var_bins"] = to_commented_map(new_bins)
            process_config["DR_SR"]["split_categories"].update(to_commented_map(new_cats))
            if "split_categories_binedges" not in process_config["DR_SR"] or not process_config["DR_SR"]["split_categories_binedges"]:
                process_config["DR_SR"]["split_categories_binedges"] = CommentedMap()
            process_config["DR_SR"]["split_categories_binedges"].update(to_commented_map(new_edges))

            if "non_closure" in process_config["DR_SR"]:
                for var, var_config in process_config["DR_SR"]["non_closure"].items():
                    if EQUIPOPULATED_BINNING_OPTIONS_KEY not in var_config:
                        continue

                    final_cuts = get_merged_cuts(var_config, parent_config=process_config["DR_SR"])
                    temp_cuts = final_cuts.replace('&&', '&').replace('!', '~').replace('||', '|')

                    console.print(f"  Processing correction: [bold magenta]DR_SR/non_closure/{var}[/bold magenta]")
                    console.print(f"    Variable to bin: [bold magenta]{var_config['var_dependence']}[/bold magenta]")
                    console.print(f"    Base cuts: [dim]{temp_cuts}[/dim]")

                    new_bins, _, new_edges, new_cats = get_binning(
                        df=dataframe,
                        binning_config=var_config[EQUIPOPULATED_BINNING_OPTIONS_KEY],
                        process_cuts=temp_cuts,
                        category_splits=var_config["split_categories"],
                        original_category_splits=var_config["split_categories"],
                        var_dependence=var_config["var_dependence"],
                    )

                    var_config["var_bins"] = to_commented_map(new_bins)
                    var_config["split_categories"].update(to_commented_map(new_cats))
                    if "split_categories_binedges" not in var_config or not var_config["split_categories_binedges"]:
                        var_config["split_categories_binedges"] = CommentedMap()
                    var_config["split_categories_binedges"].update(to_commented_map(new_edges))

    console.print(f"\n[bold]Writing updated configuration to [cyan]{args.config}[/cyan][/bold]")

    with open(args.config, "w") as f:
        yaml.dump(config, f)
