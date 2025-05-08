#!/usr/bin/env python3
"""
# Lambda Decay Analyzer

This script analyzes EDM4hep ROOT files containing Monte Carlo particle data to extract
Lambda hyperon (PDG: 3122) decay information. It categorizes Lambda particles into three cases:

1. Lambda particles with no daughters (stable/undecayed)
2. Lambda → proton + π⁻ decay channel
3. Lambda → neutron + π⁰ decay channel

For each category, the script saves a CSV dataset containing:
- Event and particle indices
- Lambda properties (momentum, endpoint)
- Decay product properties (when applicable)

The script supports processing multiple input files, printing progress updates, and
saving comprehensive statistics about the Lambda particles found.

Usage:
    python lambda_decay_analyzer.py -i FILE1 [FILE2 ...] [-o OUTPUT_PREFIX] [--debug]
"""

import argparse
import uproot
import awkward as ak
import pandas as pd
import numpy as np
import time
import os
from rich.table import Table
from rich.console import Console
from rich.progress import Progress


def fill_particle_properties(record, prefix, evt_idx, part_idx, arrays):
    """Fill particle properties into the record with a given prefix.

    Args:
        record: Dictionary to update with particle properties
        prefix: Prefix to use for the property names (e.g., 'lam', 'prot')
        evt_idx: Event index
        part_idx: Particle index
        arrays: Dictionary of arrays from the ROOT file

    Returns:
        Updated record dictionary
    """
    record.update({
        f"{prefix}_pdg": int(arrays["MCParticles.PDG"][evt_idx][part_idx]),
        f"{prefix}_mom_x": float(arrays["MCParticles.momentum.x"][evt_idx][part_idx]),
        f"{prefix}_mom_y": float(arrays["MCParticles.momentum.y"][evt_idx][part_idx]),
        f"{prefix}_mom_z": float(arrays["MCParticles.momentum.z"][evt_idx][part_idx]),
        f"{prefix}_endpoint_x": float(arrays["MCParticles.endpoint.x"][evt_idx][part_idx]),
        f"{prefix}_endpoint_y": float(arrays["MCParticles.endpoint.y"][evt_idx][part_idx]),
        f"{prefix}_endpoint_z": float(arrays["MCParticles.endpoint.z"][evt_idx][part_idx]),
    })
    return record


def process_tree(tree, max_events=None, debug=False):
    """Process the ROOT tree and extract Lambda decay information.

    Args:
        tree: The uproot tree object to process
        max_events: Maximum number of events to process (None = all)
        debug: Whether to print debug tables

    Returns:
        Tuple of (DataFrames dictionary, statistics dictionary)
    """
    # Read event data with comprehensive particle properties
    arrays = tree.arrays(
        expressions=[
            # Basic MC info
            "MCParticles.PDG",
            # Full momentum vector
            "MCParticles.momentum.x",
            "MCParticles.momentum.y",
            "MCParticles.momentum.z",
            # Full endpoint vector
            "MCParticles.endpoint.x",
            "MCParticles.endpoint.y",
            "MCParticles.endpoint.z",
            # Daughter indices
            "MCParticles.daughters_begin",
            "MCParticles.daughters_end",
            # Parent indices
            "MCParticles.parents_begin",
            "MCParticles.parents_end",
            # The actual Podio::ObjectID vectors for daughters & parents
            "_MCParticles_daughters.index",
            "_MCParticles_parents.index",
        ],
        entry_stop=max_events,
        library="ak"
    )

    console = Console()

    # Prepare lists to collect data for each decay category
    no_decay_data = []
    proton_piminus_data = []
    neutron_pizero_data = []

    # Get the total number of events
    n_events = len(arrays["MCParticles.PDG"])

    # Statistics tracking
    stats = {
        "total_events": n_events,
        "events_with_lambda": 0,
        "events_without_lambda": 0,
        "lambdas_found": 0,
        "no_decay": 0,
        "proton_piminus": 0,
        "neutron_pizero": 0,
        "other_decay": 0
    }

    # Set up progress tracking
    last_update_time = time.time()
    start_time = time.time()

    # Loop over each event
    for i_evt in range(n_events):
        # Print progress update every 5 seconds
        current_time = time.time()
        if current_time - last_update_time >= 5:
            elapsed = current_time - start_time
            events_per_sec = (i_evt + 1) / elapsed if elapsed > 0 else 0
            percent_complete = (i_evt + 1) / n_events * 100
            print(f"Processing: {i_evt + 1}/{n_events} events ({percent_complete:.1f}%) - "
                  f"{events_per_sec:.1f} events/sec - Lambdas found: {stats['lambdas_found']}")
            last_update_time = current_time

        # Get the PDG codes for this event
        event_pdg = arrays["MCParticles.PDG"][i_evt]

        # Find all Lambda particles in this event
        lambda_mask = (event_pdg == 3122)
        lambda_indices = np.where(lambda_mask)[0]

        # Count events with/without lambdas
        if len(lambda_indices) > 0:
            stats["events_with_lambda"] += 1
        else:
            stats["events_without_lambda"] += 1

        # Update total lambdas count
        stats["lambdas_found"] += len(lambda_indices)

        # For each Lambda, check its decay mode
        for row_idx, lambda_idx in enumerate(lambda_indices):
            # Get daughter info for this Lambda
            dbeg = arrays["MCParticles.daughters_begin"][i_evt][lambda_idx]
            dend = arrays["MCParticles.daughters_end"][i_evt][lambda_idx]

            # Create basic record with indices
            record = {
                "event": i_evt,
                "lam": row_idx,  # Lambda index within event's lambdas
                "lambda_index": int(lambda_idx)  # Absolute particle index
            }

            # Add Lambda properties to the record
            fill_particle_properties(record, "lam", i_evt, lambda_idx, arrays)

            # Case 1: Lambda with no daughters
            if (dbeg < 0) or (dend <= dbeg):
                no_decay_data.append(record)
                stats["no_decay"] += 1
                continue

            # Get the daughter indices for this Lambda
            daughter_indices = arrays["_MCParticles_daughters.index"][i_evt][dbeg:dend]

            # Create dictionary mapping PDG codes to particle indices for easy lookup
            daughter_by_pdg = {}
            for idx in daughter_indices:
                pdg = event_pdg[idx]
                daughter_by_pdg[pdg] = idx

            # Case 2: Lambda → proton + π-
            if 2212 in daughter_by_pdg and -211 in daughter_by_pdg:
                proton_idx = daughter_by_pdg[2212]
                piminus_idx = daughter_by_pdg[-211]

                # Add daughter properties to the record
                fill_particle_properties(record, "prot", i_evt, proton_idx, arrays)
                fill_particle_properties(record, "piminus", i_evt, piminus_idx, arrays)

                proton_piminus_data.append(record)
                stats["proton_piminus"] += 1

            # Case 3: Lambda → neutron + π0
            elif 2112 in daughter_by_pdg and 111 in daughter_by_pdg:
                neutron_idx = daughter_by_pdg[2112]
                pizero_idx = daughter_by_pdg[111]

                # Add daughter properties to the record
                fill_particle_properties(record, "neut", i_evt, neutron_idx, arrays)
                fill_particle_properties(record, "pizero", i_evt, pizero_idx, arrays)

                neutron_pizero_data.append(record)
                stats["neutron_pizero"] += 1

            else:
                # Other decay modes
                stats["other_decay"] += 1

        # Print debug tables if requested
        if debug:
            print_debug_table(i_evt, arrays, console)

    # Convert lists to pandas DataFrames
    df_no_decay = pd.DataFrame(no_decay_data)
    df_proton_piminus = pd.DataFrame(proton_piminus_data)
    df_neutron_pizero = pd.DataFrame(neutron_pizero_data)

    # Create dictionary of DataFrames
    df_dict = {
        "no_decay": df_no_decay,
        "proton_piminus": df_proton_piminus,
        "neutron_pizero": df_neutron_pizero
    }

    return df_dict, stats


def print_debug_table(i_evt, arrays, console):
    """Print debug tables for an event.

    Args:
        i_evt: Event index
        arrays: Dictionary of arrays from the ROOT file
        console: Rich console object for printing
    """
    table = Table(title=f"Event {i_evt}")
    table.expand = True  # allow full width
    table.add_column("Idx", no_wrap=True)
    table.add_column("PDG", no_wrap=True)
    table.add_column("pz", no_wrap=True)
    table.add_column("endZ", no_wrap=True)
    table.add_column("Parents (idx)", no_wrap=True)
    table.add_column("Daughters (idx)", no_wrap=True)

    event_pdg = arrays["MCParticles.PDG"][i_evt]
    event_pz = arrays["MCParticles.momentum.z"][i_evt]
    event_endz = arrays["MCParticles.endpoint.z"][i_evt]
    event_daughters_begin = arrays["MCParticles.daughters_begin"][i_evt]
    event_daughters_end = arrays["MCParticles.daughters_end"][i_evt]
    event_parents_begin = arrays["MCParticles.parents_begin"][i_evt]
    event_parents_end = arrays["MCParticles.parents_end"][i_evt]
    event_daughters_idx = arrays["_MCParticles_daughters.index"][i_evt]
    event_parents_idx = arrays["_MCParticles_parents.index"][i_evt]

    n_parts = len(event_pdg)

    for i_part in range(n_parts):
        dbeg = event_daughters_begin[i_part]
        dend = event_daughters_end[i_part]
        pbeg = event_parents_begin[i_part]
        pend = event_parents_end[i_part]

        # Build a short string for the parent's indexes
        if (pbeg >= 0) and (pend > pbeg):
            these_par_idxs = event_parents_idx[pbeg:pend]
            par_str = ", ".join(str(x) for x in these_par_idxs)
        else:
            par_str = "None"

        # Build a short string for the daughter's indexes
        if (dbeg >= 0) and (dend > dbeg):
            these_dau_idxs = event_daughters_idx[dbeg:dend]
            dau_str = ", ".join(str(x) for x in these_dau_idxs)
        else:
            dau_str = "None"

        table.add_row(
            str(i_part),
            str(event_pdg[i_part]),
            f"{event_pz[i_part]:.2f}",
            f"{event_endz[i_part]:.2f}",
            par_str,
            dau_str
        )

    console.print(table)


def save_dataframes(df_dict, output_prefix, descriptions):
    """Save DataFrames to CSV format.

    Args:
        df_dict: Dictionary of DataFrames
        output_prefix: Prefix for output files
        descriptions: Dictionary of descriptions for each DataFrame
    """
    for name, df in df_dict.items():
        if df.empty:
            print(f"No data for {name}, skipping...")
            continue

        output_file = f"{output_prefix}_{name}.csv"
        df.to_csv(output_file, index=False)
        print(f"Saved {len(df)} {descriptions[name]} to {output_file}")


def generate_statistics_markdown(stats, input_files, elapsed_time, output_prefix):
    """Generate statistics in Markdown format.

    Args:
        stats: Dictionary of statistics
        input_files: List of input files processed
        elapsed_time: Total processing time in seconds
        output_prefix: Prefix for output files

    Returns:
        Markdown formatted string with statistics
    """
    # Calculate percentages of decay modes
    total_decays = stats["lambdas_found"]
    pct_no_decay = (stats["no_decay"] / total_decays * 100) if total_decays > 0 else 0
    pct_proton_piminus = (stats["proton_piminus"] / total_decays * 100) if total_decays > 0 else 0
    pct_neutron_pizero = (stats["neutron_pizero"] / total_decays * 100) if total_decays > 0 else 0
    pct_other = (stats["other_decay"] / total_decays * 100) if total_decays > 0 else 0

    # Format input files list for markdown
    files_list = "\n".join([f"- {os.path.basename(f)}" for f in input_files])

    # Create markdown content - use ASCII instead of Unicode for better compatibility
    markdown = f"""# Lambda Decay Analysis Results

## Processing Summary

- **Files Processed**: {len(input_files)}
- **Total Processing Time**: {elapsed_time:.2f} seconds

### Input Files:
{files_list}

## Statistics

| Statistic | Value | Percentage |
|-----------|-------|------------|
| Total Events | {stats["total_events"]} | 100% |
| Events with Lambda | {stats["events_with_lambda"]} | {stats["events_with_lambda"] / stats["total_events"] * 100:.1f}% |
| Events without Lambda | {stats["events_without_lambda"]} | {stats["events_without_lambda"] / stats["total_events"] * 100:.1f}% |
| **Total Lambda particles** | **{stats["lambdas_found"]}** | **100%** |
| Lambda with no decay | {stats["no_decay"]} | {pct_no_decay:.1f}% |
| Lambda -> p + pi- | {stats["proton_piminus"]} | {pct_proton_piminus:.1f}% |
| Lambda -> n + pi0 | {stats["neutron_pizero"]} | {pct_neutron_pizero:.1f}% |
| Lambda with other decay | {stats["other_decay"]} | {pct_other:.1f}% |

## Outputs

The following CSV files were created:

- `{output_prefix}_no_decay.csv`: Lambda particles with no decay
- `{output_prefix}_proton_piminus.csv`: Lambda -> p + pi- decays
- `{output_prefix}_neutron_pizero.csv`: Lambda -> n + pi0 decays

Each file contains detailed kinematic information about the Lambda particles and their decay products.
"""
    return markdown


def print_statistics(stats, input_files, elapsed_time, output_prefix):
    """Print statistics about the processed data."""
    markdown = generate_statistics_markdown(stats, input_files, elapsed_time, output_prefix)
    print("\n" + markdown)


def save_statistics_markdown(stats, output_prefix, input_files, elapsed_time):
    """Save statistics to a Markdown file."""
    # Generate markdown content
    markdown = generate_statistics_markdown(stats, input_files, elapsed_time, output_prefix)

    # Save markdown to file with explicit UTF-8 encoding
    output_file = f"{output_prefix}_statistics.md"
    with open(output_file, "w", encoding="utf-8") as f:
        f.write(markdown)

    print(f"Saved statistics to {output_file}")


def process_input_file(input_file, tree_name, max_events, debug):
    """Process a single input file.

    Args:
        input_file: Path to the input file
        tree_name: Name of the tree in the ROOT file
        max_events: Maximum number of events to process
        debug: Whether to print debug tables

    Returns:
        Tuple of (DataFrames dictionary, statistics dictionary)
    """
    print(f"Opening file: {input_file}")
    with uproot.open(input_file) as f:
        tree = f[tree_name]
        return process_tree(tree, max_events=max_events, debug=debug)


def merge_dataframes(df_dicts):
    """Merge multiple DataFrame dictionaries.

    Args:
        df_dicts: List of DataFrame dictionaries

    Returns:
        Merged DataFrame dictionary
    """
    # Initialize with empty DataFrames
    merged = {
        "no_decay": pd.DataFrame(),
        "proton_piminus": pd.DataFrame(),
        "neutron_pizero": pd.DataFrame()
    }

    # Concatenate each DataFrame
    for df_dict in df_dicts:
        for name, df in df_dict.items():
            if not df.empty:
                if merged[name].empty:
                    merged[name] = df
                else:
                    merged[name] = pd.concat([merged[name], df], ignore_index=True)

    return merged


def merge_statistics(stats_list):
    """Merge multiple statistics dictionaries.

    Args:
        stats_list: List of statistics dictionaries

    Returns:
        Merged statistics dictionary
    """
    merged = {
        "total_events": 0,
        "events_with_lambda": 0,
        "events_without_lambda": 0,
        "lambdas_found": 0,
        "no_decay": 0,
        "proton_piminus": 0,
        "neutron_pizero": 0,
        "other_decay": 0
    }

    # Sum all statistics
    for stats in stats_list:
        for key in merged:
            merged[key] += stats[key]

    return merged


def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Lambda decay analyzer: Extracts Lambda hyperon decays into pandas DataFrames.")
    parser.add_argument("-i", "--input-files", required=True, nargs='+',
                        help="Path(s) to one or more EDM4hep ROOT files with MCParticles.")
    parser.add_argument("-t", "--tree-name", default="events",
                        help="Name of the TTree (default 'events').")
    parser.add_argument("-o", "--output-prefix", default="lambda_decays",
                        help="Prefix for output files (default 'lambda_decays').")
    parser.add_argument("--max-events", type=int, default=None,
                        help="Max number of events to read per file (default=all).")
    parser.add_argument("--debug", action="store_true",
                        help="Print debug tables for events.")
    args = parser.parse_args()

    # Create descriptions for the DataFrames
    descriptions = {
        "no_decay": "Lambda particles with no decay",
        "proton_piminus": "Lambda → p + π- decays",
        "neutron_pizero": "Lambda → n + π0 decays"
    }

    # Lists to collect results from each file
    all_df_dicts = []
    all_stats = []

    # Process each input file
    start_time = time.time()
    for input_file in args.input_files:
        df_dict, stats = process_input_file(
            input_file,
            args.tree_name,
            args.max_events,
            args.debug
        )
        all_df_dicts.append(df_dict)
        all_stats.append(stats)

    # Merge results from all files
    merged_df_dict = merge_dataframes(all_df_dicts)
    merged_stats = merge_statistics(all_stats)

    # Calculate total processing time
    total_time = time.time() - start_time

    # Print overall statistics
    print("\nProcessing complete!")
    print(f"Total processing time: {total_time:.2f} seconds")
    print_statistics(merged_stats)

    # Save the merged DataFrames
    save_dataframes(merged_df_dict, args.output_prefix, descriptions)

    # Save statistics to markdown file
    save_statistics_markdown(merged_stats, args.output_prefix, args.input_files, total_time)


if __name__ == "__main__":
    main()