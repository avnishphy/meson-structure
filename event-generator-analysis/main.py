#!/usr/bin/env python3
"""
Main script to analyze k_lambda data using uproot and scikit-hep/hist.
This version processes multiple trees simultaneously (Evnts and Process) and
adds Lambda kinematic histograms.
"""

import argparse
import uproot

# Local modules
from histograms import create_histograms, fill_histograms
from plotting import make_plots


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Analyze k_lambda ROOT TTree (multiple trees) using scikit-hep hist.")
    parser.add_argument("--input-file", "-i", required=True, help="Path to the input ROOT file.")
    parser.add_argument("--chunk-size", "-c", type=int, default=100_000, help="Number of events to process per chunk (default: %(default)s).")
    parser.add_argument("--events", "-e", type=int, default=None, help="Maximum number of events to process (default: all).")
    parser.add_argument("--outdir", "-o", default="plots", help="Directory to save output plots (default: %(default)s).")
    return parser.parse_args()


def main():
    args = parse_args()

    # Open the file and get the trees we need
    file = uproot.open(args.input_file)
    tree_evnts = file["Evnts"]
    print(f"TREE Evnts has: {tree_evnts.num_entries} events")
    tree_process = file["Process"]
    # (If you need Meta, add: tree_meta = file["Meta"])

    # Define the branches to read from each tree.
    branches_evnts = [
        "TDIS_Q2", "TDIS_xbj", "TDIS_y", "TDIS_t",
        "pkx_Lab", "pky_Lab",
        "lamb_scat"  # assumed TLorentzVector with fields like "px", "py", "pz", "E"
    ]
    branches_process = [
        "ElambE_Lab"  # Lambda energy from Process tree
    ]

    # Create histograms (including new lambda parameters)
    hists = create_histograms()

    # Create iterators for each tree. The events are assumed to be in the same order.
    iter_evnts = tree_evnts.iterate(
        expressions=branches_evnts,
        step_size=args.chunk_size
    )
    iter_process = tree_process.iterate(
        expressions=branches_process,
        step_size=args.chunk_size
    )

    events_processed = 0
    for chunk_evnts, chunk_process in zip(iter_evnts, iter_process):
        # Fill histograms using data from both trees
        fill_histograms(hists, chunk_evnts, chunk_process)

        events_processed += len(chunk_evnts["TDIS_Q2"])
        if args.max_events and events_processed >= args.max_events:
            break

    print(f"Total events processed: {events_processed}")

    # Create and save plots
    make_plots(hists, outdir=args.outdir)


if __name__ == "__main__":
    main()
