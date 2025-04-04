#!/usr/bin/env python3

"""
analyze_kinematics.py

Reads the 'Evnts' TTree from a large ROOT file. For each of the branches:
    - "P_Inc."
    - "e_Inc."
    - "e_Scat."
    - "k."
    - "lamb_scat."

we build 1D histograms of:
  1) pT = sqrt(px^2 + py^2)
  2) pz
  3) theta in milliradians
  4) phi in milliradians

We store them in scikit-hep/hist objects and produce PNG plots in 'plots/' by default.

Usage example:
  python analyze_kinematics.py --input-file k_lambda_crossing_0_10.0on100.0.root --max-events 100000
"""

import argparse
import os
import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
import hist
from hist import Hist

branches  = ["P_Inc.", "e_Inc.", "e_Scat.", "k.",   "lamb_scat."]
particles = ["inc_p" , "inc_e",  "scat_e",  "kaon", "lambda"]


def create_histograms():
    """
    Returns a dict of hist.Hist objects for each particle's:
      - pt
      - pz
      - theta_mrad
      - phi_mrad

    We'll index them by (particle_name, var_name).
    """
    # Binning choices
    pt_edges       = hist.axis.Regular(100, 0, 10, name="pt", label="pT (GeV/c)")
    pz_edges       = hist.axis.Regular(200, -100, 100, name="pz", label="pz (GeV/c)")
    p_hist         = hist.axis.Regular(200, -100, 100, name="P", label="P-total (GeV/c)")
    theta_mrad_ax  = hist.axis.Regular(180, 0, 180, name="theta_mrad", label=r"$\theta$ (mrad)")
    # phi from -180 to 180, in mrad
    phi_mrad_ax    = hist.axis.Regular(360, -180, 180, name="phi_mrad", label=r"$\phi$ (mrad)")

    # We have 5 particle branches


    # We'll build a dictionary: hists[(part, "pt")] = Hist(...)
    hists = {}
    for part in particles:
        # 1) pT
        h_pt = Hist(pt_edges, storage=hist.storage.Double())
        # 2) pz
        h_pz = Hist(pz_edges, storage=hist.storage.Double())
        # 3) theta in mrad
        h_th = Hist(theta_mrad_ax, storage=hist.storage.Double())
        # 4) phi in mrad
        h_phi = Hist(phi_mrad_ax, storage=hist.storage.Double())
        # Total momentum
        h_p = Hist(p_hist, storage=hist.storage.Double())

        hists[(part, "pt")]         = h_pt
        hists[(part, "pz")]         = h_pz
        hists[(part, "theta_mrad")] = h_th
        hists[(part, "phi_mrad")]   = h_phi
        hists[(part, "p")]   = h_p

    return hists


def fill_histograms(hists, chunk):
    """
    Extract the relevant variables from each particle branch, compute:
      - pT = sqrt(px^2 + py^2)
      - pz
      - theta [mrad]
      - phi [mrad]

    Then fill the corresponding histograms in the hists dict.
    """
    # List of particle branches
    particles = ["P_Inc.", "e_Inc.", "e_Scat.", "k.", "lamb_scat."]

    # For each branch, we have awkward arrays where each entry is a record:
    # { fP: { fX, fY, fZ }, fE }
    # We'll get px = fP.fX, py = fP.fY, pz = fP.fZ
    for part in particles:
        arr = chunk[part]
        # arr has length = # events in chunk
        px = arr["fP"]["fX"]
        py = arr["fP"]["fY"]
        pz = arr["fP"]["fZ"]

        # compute derived
        pt = np.sqrt(px**2 + py**2)

        # total momentum magnitude
        p_mag = np.sqrt(px**2 + py**2 + pz**2 + 1e-30)  # avoid zero in divisions

        # angle definitions
        # polar angle: theta = arctan2(pt, pz) in [0..pi], in radians
        theta = np.arctan2(pt, pz)
        # or eqv. theta = np.arccos(pz / p_mag) => same result
        # we want in milliradians => multiply by 1e3
        theta_mrad = theta * 1e3

        # azimuth: phi = arctan2(py, px), in [-pi, pi], convert to mrad
        phi = np.arctan2(py, px)
        phi_mrad = phi * 1e3 * (180 / np.pi) / 180  # Actually let's do consistent approach

        # Actually, let's do direct: phi (rad) * 1e3 => mrad
        # But maybe we want phi from -180..180 in degrees, then *1e3?
        # The user specifically asked for "Phi (in milirads)", so let's do phi (rad)->(rad*mrad)
        phi_mrad = phi * 1e3

        # fill
        hists[(part, "pt")].fill(pt=pt)
        hists[(part, "pz")].fill(pz=pz)
        hists[(part, "theta_mrad")].fill(theta_mrad=theta_mrad)
        hists[(part, "phi_mrad")].fill(phi_mrad=phi_mrad)
        hists[(part, "p")].fill(p_mag)


def plot_histograms(hists, outdir="plots"):
    """
    For each (particle, var) in hists, produce a 1D matplotlib plot and save it.
    """
    os.makedirs(outdir, exist_ok=True)

    for key, hist_obj in hists.items():
        part, var = key  # e.g. ("P_Inc.", "pt")

        # We'll create a simple figure
        fig, ax = plt.subplots(figsize=(6,5))
        hist_obj.plot1d(ax=ax)
        ax.set_title(f"{part} - {var}")
        # The Hist should have an axis label from the hist creation steps
        # but we can do a fallback
        if hist_obj.axes[0].label is not None:
            ax.set_xlabel(hist_obj.axes[0].label)
        else:
            ax.set_xlabel(var)
        ax.set_ylabel("Counts")
        fig.tight_layout()

        # sanitize the part name for filename usage
        part_sanitized = part.replace(".", "")  # remove dot
        plt_filename = f"{part_sanitized}_{var}.png"
        plt_path = os.path.join(outdir, plt_filename)
        plt.savefig(plt_path)
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Analyze kinematics of final state particles from Evnts TTree.")
    parser.add_argument("--input-file", "-i", required=True, help="Path to the ROOT file with the Evnts TTree.")
    parser.add_argument("--outdir", "-o", default="plots", help="Directory to save output plots (default: plots).")
    parser.add_argument("--max-events", "-m", type=int, default=None, help="Maximum number of events to process (default: all).")
    parser.add_argument("--chunk-size", type=int, default=100_000, help="How many TTree entries to read per chunk (default: 100000).")
    args = parser.parse_args()

    # create histograms
    hists = create_histograms()

    # read TTree in chunks
    tree_name = "Evnts"


    events_processed = 0
    for data_chunk in uproot.iterate(
            f"{args.input_file}:{tree_name}",
            expressions=branches,
            library="ak",
            step_size=args.chunk_size
    ):
        fill_histograms(hists, data_chunk)
        # count how many events in this chunk
        chunk_len = len(data_chunk[branches[0]])  # or any branch
        events_processed += chunk_len
        if args.max_events and events_processed >= args.max_events:
            break

    print(f"Processed {events_processed} events in total.")
    # now produce plots
    plot_histograms(hists, outdir=args.outdir)


if __name__ == "__main__":
    main()
