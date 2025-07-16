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
import helper.functions as func

EQUIPOPULATED_BINNING_OPTIONS_KEY = "equipopulated_binning_options"

parser = argparse.ArgumentParser(description="Equipopulated binning adjustments of configuration.")
parser.add_argument(
    "--config",
    type=str,
    required=True,
    help="Path to the YAML configuration that will be adjusted.",
)
parser.add_argument(
    "--cuts-config",
    type=str,
    default=None,
    help="Path to a YAML config that contains cuts for processes (can be the fake factors config). Will not be adjusted.",
)
parser.add_argument(
    "--processes",
    type=str,
    nargs="+",
    required=False,
    help="List of processes to adjust binning for.",
    default=["QCD", "Wjets", "ttbar", "process_fractions"],
)
parser.add_argument(
    "--cut-region",
    type=str,
    default="SRlike",
    help="Cut region to apply (e.g., SRlike, ARlike), process_fractions will always use AR_cut."
)
parser.add_argument(
    "--dry-run",
    action="store_true",
    help="If set, do not modify the config file.",
)


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


def get_n_bins(
    n_bins_config: Union[int, list, CommentedSeq, dict, CommentedMap],
    category_splits: dict,
    full_cat_path: list,
) -> Union[int, None]:
    """
    Traverse the nested n_bins_config (lists of lists or dicts) according to full_cat_path and
    category_splits and search for the number of bins to use for the given category path.

    Args:
        n_bins_config (int, list, CommentedSeq, dict, CommentedMap): Configuration for the number of bins.
        category_splits (dict): Dictionary defining how to split categories.
        full_cat_path (list): List of category keys representing the path to traverse.

    Returns:
        int or None: The number of bins for the given category path, or None if not found.
    """
    value, prev_key, splits = n_bins_config, None, list(category_splits)

    for split, cat_key in zip(splits, full_cat_path):
        if isinstance(value, int):
            return value

        if isinstance(value, (dict, CommentedMap)):
            value = value.get(cat_key)
            if value is None:
                return None
            continue

        cats = category_splits[split]
        if isinstance(cats, dict):
            cats = cats.get(prev_key, [])

        try:
            idx = cats.index(cat_key)
            if not isinstance(value, (list, CommentedSeq)) or idx >= len(value):
                return None
            value, prev_key = value[idx], cat_key
        except ValueError:
            return None

    return value if isinstance(value, int) else None


def get_binning(
    df: pd.DataFrame,
    binning_config: dict,
    process_cuts: str,
    category_splits: dict,
    var_dependence: str,
    original_category_splits: Union[dict, None] = None,
    prepend_space: int = 0,
    parent_cat_keys: Union[list, None] = None,
    generated_categories_for_splits: Union[dict, None] = None,
) -> tuple:
    """
    Recursively calculates equipopulated binning for a variable based on category splits.

    Args:
        df (pd.DataFrame): The input DataFrame containing the data.
        binning_config (dict): Configuration for binning, including n_bins and variable_config.
        process_cuts (str): Cuts to apply to the DataFrame.
        category_splits (dict): Dictionary defining how to split categories.
        var_dependence (str): The variable to be binned.
        original_category_splits (dict): Original category splits for reference.
        prepend_space (int): Number of spaces to prepend for formatting.
        parent_cat_keys (list, optional): List of parent category keys for recursive calls.
        generated_categories_for_splits (dict, optional): Dynamically generated categories.

    Returns:
        tuple: A tuple containing:
    """
    original_category_splits = original_category_splits or category_splits
    parent_cat_keys = parent_cat_keys or []
    generated_categories_for_splits = generated_categories_for_splits or {}

    current_split_var = next(iter(category_splits))
    generated_edges = {}
    generated_categories = {}
    remaining_splits = {k: v for k, v in category_splits.items() if k != current_split_var}

    current_split_value = category_splits[current_split_var]
    if isinstance(current_split_value, (list, CommentedSeq)) and current_split_value:
        categories = current_split_value
        generated_categories[current_split_var] = categories
        generated_categories_for_splits[current_split_var] = categories
    elif (current_split_var in binning_config.get("n_bins", {}) and isinstance(binning_config["n_bins"][current_split_var], int)):
        n_split_bins = binning_config["n_bins"][current_split_var]
        var_conf = binning_config["variable_config"][current_split_var]
        current_df = df.query(process_cuts)
        split_edges = _equipopulated_binned_variable(current_df[current_split_var], n_split_bins)

        if split_edges.size > 0:
            split_edges[0] = var_conf.get("min", split_edges[0])
            split_edges[-1] = var_conf.get("max", split_edges[-1])

        rounding = var_conf.get("rounding", 2)
        generated_edges[current_split_var] = np.round(split_edges, rounding).tolist()

        categories = []
        for i, (lower, upper) in enumerate(pairwise(split_edges)):
            lower, upper = round(lower, rounding), round(upper, rounding)

            if i == 0:  # syntax setting
                categories.append(f"<={upper}")
            elif i == len(split_edges) - 2:
                categories.append(f">{lower}")
            else:
                categories.append(f">{lower}#&&#<={upper}")

        generated_categories[current_split_var] = categories
        generated_categories_for_splits[current_split_var] = categories
    else:
        raise ValueError(
            f"Unable to determine categories for '{current_split_var}' "
            "neither through equipopulated binning option nor explicit categories."
        )

    table = Table(show_header=True, header_style="bold magenta")
    category_header = f"Category: {current_split_var}"
    if parent_cat_keys:
        original_split_vars = list(original_category_splits.keys())
        parent_var_name = original_split_vars[len(parent_cat_keys) - 1]
        parent_key = parent_cat_keys[-1]
        category_header += f" (in {parent_var_name} {parent_key})"

    for col_name in [category_header, "Num Events", f"Bin Edges ({var_dependence})", "Num Events/Bin"]:
        table.add_column(col_name)

    output_bins = OrderedDict()

    for cat_key in categories:
        query_cat_cut = (cat_key if any(op in cat_key for op in ["<", ">", "=="]) else f"{current_split_var} {cat_key}")

        for _sep, _join in [("#&&#", " & "), ("#||#", " | ")]:
            if _sep in query_cat_cut:
                sub_cuts = query_cat_cut.split(_sep)
                query_cat_cut = f"({current_split_var} {sub_cuts[0]} {_join} {current_split_var} {sub_cuts[1]})"
                break
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
                generated_categories_for_splits=generated_categories_for_splits,
            )
            output_bins[cat_key] = nested_bins
            for var, edges in nested_edges.items():
                generated_edges.setdefault(var, {})[cat_key] = edges
            for var, cats in nested_cats.items():
                generated_categories.setdefault(var, {})[cat_key] = cats
        else:  # final level of splitting
            target_variable = var_dependence
            n_bins = get_n_bins(
                n_bins_config=binning_config.get("var_dependence_n_bins"),
                category_splits=generated_categories_for_splits,
                full_cat_path=current_cat_keys,
            )

            if n_bins and n_bins > 0:
                var_conf = binning_config.get("variable_config", {}).get(target_variable, {})

                filtered_df = sub_df
                if (min_val := var_conf.get("min")) is not None:
                    filtered_df = filtered_df.query(f"{target_variable} >= {min_val}")
                if (max_val := var_conf.get("max")) is not None:
                    filtered_df = filtered_df.query(f"{target_variable} <= {max_val}")

                bins = _equipopulated_binned_variable(filtered_df[target_variable], n_bins)
                if bins.size > 0:
                    bins[0] = var_conf.get("min", bins[0])
                    bins[-1] = var_conf.get("max", bins[-1])

                if (rounding := var_conf.get("rounding", 2)):
                    bins = np.round(bins, rounding).tolist()

                for option in ["add_left", "add_right"]:
                    if option not in binning_config:
                        continue

                    temp_config = binning_config[option]
                    split_vars = list(original_category_splits.keys())[:len(parent_cat_keys)]

                    try:
                        for split_var, key in zip(split_vars, parent_cat_keys):
                            temp_config = temp_config[split_var][key]

                        if isinstance(temp_config, (int, float, list, tuple)):
                            value_to_add = temp_config if isinstance(temp_config, (list, tuple)) else [temp_config]
                            bins = value_to_add + bins if option == "add_left" else bins + value_to_add
                    except (KeyError, TypeError):
                        continue

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

    cuts_config = None
    if args.cuts_config:
        console.print(f"Loading cuts from [cyan]{args.cuts_config}[/cyan]")
        with open(Path(args.cuts_config).resolve(), "r") as f:
            cuts_config = func.configured_yaml.load(f)

    with open((args.config.parent / "common_settings.yaml").resolve(), "r") as f:
        common_settings = func.configured_yaml.load(f)

    with open(args.config, "r") as f:
        config = func.configured_yaml.load(f)

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

    if (bool_cols := dataframe.select_dtypes(include=['bool']).columns).any():
        dataframe[bool_cols] = dataframe[bool_cols].astype(int)
        console.print(f"Converted {list(bool_cols)} boolean columns to integers.")

    def get_merged_cuts(
        var_config: Union[dict, CommentedMap],
        parent_config: Union[dict, CommentedMap, None] = None,
    ) -> str:
        """Merge cuts for corrections, taking into account parent configurations."""
        merged_cuts = {}

        # Start from external file if provided
        if cuts_config:
            cuts_source = cuts_config["target_processes"][process]
            merged_cuts = cuts_source.get(f"{args.cut_region}_cuts", {}).copy()

        # Overwrite with parent cuts if available
        if parent_config:
            merged_cuts.update(parent_config.get(f"{args.cut_region}_cuts", {}))

        # Apply local cuts
        merged_cuts.update(var_config.get(f"{args.cut_region}_cuts", {}))

        return " & ".join(f"({c})" for c in merged_cuts.values())

    for process in args.processes:
        console.print(f"\n[bold green]Processing: {process}[/bold green]")

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
                var_dependence=process_config["var_dependence"],
            )

            process_config["var_bins"] = func.to_commented_map(new_bins)
            process_config["split_categories"].update(func.to_commented_map(new_cats))
            process_config.setdefault("split_categories_binedges", CommentedMap()).update(func.to_commented_map(new_edges))
        else:
            console.print(f"[red]Process '{process}' does not have equipopulated_binning_options. Skipping.[/red]")

        if "non_closure" in process_config:
            for var, var_config in process_config["non_closure"].items():
                if EQUIPOPULATED_BINNING_OPTIONS_KEY not in var_config:
                    console.print(f"[red]Correction '{var}' does not have equipopulated_binning_options. Skipping.[/red]")
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
                    var_dependence=var_config["var_dependence"],
                )

                var_config["var_bins"] = func.to_commented_map(new_bins)
                var_config["split_categories"].update(func.to_commented_map(new_cats))
                var_config.setdefault("split_categories_binedges", CommentedMap()).update(func.to_commented_map(new_edges))

        if "DR_SR" in process_config:
            dr_sr_config = process_config["DR_SR"]
            if EQUIPOPULATED_BINNING_OPTIONS_KEY not in dr_sr_config:
                console.print("[red]DR_SR does not have equipopulated_binning_options. Skipping.[/red]")
                continue

            final_cuts = get_merged_cuts(dr_sr_config)
            temp_cuts = final_cuts.replace('&&', '&').replace('!', '~').replace('||', '|')

            console.print("  Processing correction: [bold magenta]DR_SR[/bold magenta]")
            console.print(f"    Variable to bin: [bold magenta]{dr_sr_config['var_dependence']}[/bold magenta]")
            console.print(f"    Base cuts: [dim]{temp_cuts}[/dim]")

            new_bins, _, new_edges, new_cats = get_binning(
                df=dataframe,
                binning_config=dr_sr_config[EQUIPOPULATED_BINNING_OPTIONS_KEY],
                process_cuts=temp_cuts,
                category_splits=dr_sr_config["split_categories"],
                var_dependence=dr_sr_config["var_dependence"],
            )

            dr_sr_config["var_bins"] = func.to_commented_map(new_bins)
            dr_sr_config["split_categories"].update(func.to_commented_map(new_cats))
            dr_sr_config.setdefault("split_categories_binedges", CommentedMap()).update(func.to_commented_map(new_edges))

            if "non_closure" in dr_sr_config:
                for var, var_config in dr_sr_config["non_closure"].items():
                    if EQUIPOPULATED_BINNING_OPTIONS_KEY not in var_config:
                        console.print(f"[red]Correction '{var}' does not have equipopulated_binning_options. Skipping.[/red]")
                        continue

                    final_cuts = get_merged_cuts(var_config, parent_config=dr_sr_config)
                    temp_cuts = final_cuts.replace('&&', '&').replace('!', '~').replace('||', '|')

                    console.print(f"  Processing correction: [bold magenta]DR_SR/non_closure/{var}[/bold magenta]")
                    console.print(f"    Variable to bin: [bold magenta]{var_config['var_dependence']}[/bold magenta]")
                    console.print(f"    Base cuts: [dim]{temp_cuts}[/dim]")

                    new_bins, _, new_edges, new_cats = get_binning(
                        df=dataframe,
                        binning_config=var_config[EQUIPOPULATED_BINNING_OPTIONS_KEY],
                        process_cuts=temp_cuts,
                        category_splits=var_config["split_categories"],
                        var_dependence=var_config["var_dependence"],
                    )

                    var_config["var_bins"] = func.to_commented_map(new_bins)
                    var_config["split_categories"].update(func.to_commented_map(new_cats))
                    var_config.setdefault("split_categories_binedges", CommentedMap()).update(func.to_commented_map(new_edges))

    if not args.dry_run:
        console.print(f"\n[bold]Writing updated configuration to [cyan]{args.config}[/cyan][/bold]")
        with open(args.config, "w") as f:
            func.configured_yaml.dump(config, f)
    else:
        console.print("[bold red]Dry run mode: configuration not modified.[/bold red]")
