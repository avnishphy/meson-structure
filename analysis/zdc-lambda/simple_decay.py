#!/usr/bin/env python3

import argparse
import re
import os
import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot Lambda (PDG=3122) decay endpoint.z and decay product angular distributions."
    )
    parser.add_argument("-i", "--input-files", nargs="+", required=True,
                        help="List of EDM4hep ROOT files containing MCParticles.")
    parser.add_argument("-t", "--tree-name", default="events",
                        help="Name of the TTree (default 'events').")
    parser.add_argument("-o", "--outdir", default="plots",
                        help="Output directory for PDF/PNG plots.")
    parser.add_argument("--step-size", default="1000",
                        help="Number of events to read per chunk (integer).")
    return parser.parse_args()

def extract_beam_label(filename):
    """
    Attempt to parse something like '_5x41_', '_10x100_', '_18x275_'
    from the filename. Return that substring if found, else 'unknown'.
    """
    base = os.path.basename(filename)
    match = re.search(r'_(\d+x\d+)_', base)
    if match:
        return match.group(1)
    return "unknown"

# PDG codes
PDG_LAMBDA = 3122
PDG_PROTON = 2212
PDG_PION_MINUS = -211

def calculate_theta_phi(px, py, pz):
    """Calculate theta and phi angles from momentum components."""
    p = np.sqrt(px**2 + py**2 + pz**2)
    # Guard against p=0 just in case
    if p == 0:
        return 0.0, 0.0
    theta = np.arccos(pz / p) * 180.0 / np.pi  # Convert to degrees
    phi = np.arctan2(py, px) * 180.0 / np.pi   # Convert to degrees
    return theta, phi

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Dictionary to store Lambda decay data
    data_dict = {}

    # Color map for the overlays
    color_map = {
        "5x41": "red",
        "10x100": "blue",
        "18x275": "green",
    }

    # 1) Group input files by beam label
    beams_map = {}
    for f in args.input_files:
        beam = extract_beam_label(f)
        beams_map.setdefault(beam, []).append(f)

    # 2) Read in chunks of step_size events, gather Lambda endpoints and decay products
    for beam_label, file_list in beams_map.items():
        print(f"\n=== Beam: {beam_label} ===")

        # Initialize data structures for this beam if needed
        if beam_label not in data_dict:
            data_dict[beam_label] = {
                'endpoint_z': [],
                'proton_theta': [],
                'proton_phi': [],
                'pion_theta': [],
                'pion_phi': []
            }

        # Convert step_size to integer
        step_n = int(args.step_size)

        # Loop over each file for this beam
        for fname in file_list:
            print(f"   Reading: {fname}")

            # We'll make a small dict for uproot's iterate
            file_dict = {fname: args.tree_name}

            # Use uproot.iterate with axis="entry" so that partial events are NEVER split
            for chunk in uproot.iterate(
                    file_dict,
                    expressions=[
                        "MCParticles.PDG",
                        "MCParticles.endpoint.z",
                        "MCParticles.momentum.x",
                        "MCParticles.momentum.y",
                        "MCParticles.momentum.z",
                        "MCParticles.daughters_begin",
                        "MCParticles.daughters_end"
                    ],
                    step_size=step_n,    # chunk size in number of events
                    axis="entry",        # preserves event boundaries
                    library="ak",
            ):
                pdg = chunk["MCParticles.PDG"]
                endpointz = chunk["MCParticles.endpoint.z"]
                px = chunk["MCParticles.momentum.x"]
                py = chunk["MCParticles.momentum.y"]
                pz = chunk["MCParticles.momentum.z"]
                daughters_begin = chunk["MCParticles.daughters_begin"]
                daughters_end = chunk["MCParticles.daughters_end"]

                n_events = len(pdg)  # number of events in this chunk
                for i_evt in range(n_events):
                    # Find all Lambdas in this event
                    is_lambda_evt = (pdg[i_evt] == PDG_LAMBDA)
                    lambda_indices = np.where(is_lambda_evt)[0]

                    # Process each Lambda
                    for lambda_idx in lambda_indices:
                        # Store endpoint.z
                        data_dict[beam_label]['endpoint_z'].append(endpointz[i_evt][lambda_idx])

                        # Check if this Lambda has daughters
                        d_begin = daughters_begin[i_evt][lambda_idx]
                        d_end = daughters_end[i_evt][lambda_idx]

                        # If d_begin == -1 or d_end == -1, there's no daughter info
                        if d_end > d_begin:
                            # Get PDG codes of daughters
                            d_pdgs = pdg[i_evt][d_begin:d_end]

                            # Check for proton + pion- among daughters
                            has_proton = PDG_PROTON in d_pdgs
                            has_pion_minus = PDG_PION_MINUS in d_pdgs

                            if has_proton and has_pion_minus:
                                # Indices of the proton and pion
                                proton_idx = d_begin + np.where(d_pdgs == PDG_PROTON)[0][0]
                                pion_idx = d_begin + np.where(d_pdgs == PDG_PION_MINUS)[0][0]

                                # Get momentum components
                                p_px = px[i_evt][proton_idx]
                                p_py = py[i_evt][proton_idx]
                                p_pz = pz[i_evt][proton_idx]

                                pi_px = px[i_evt][pion_idx]
                                pi_py = py[i_evt][pion_idx]
                                pi_pz = pz[i_evt][pion_idx]

                                # Calculate angles
                                p_theta, p_phi = calculate_theta_phi(p_px, p_py, p_pz)
                                pi_theta, pi_phi = calculate_theta_phi(pi_px, pi_py, pi_pz)

                                # Store them
                                data_dict[beam_label]['proton_theta'].append(p_theta)
                                data_dict[beam_label]['proton_phi'].append(p_phi)
                                data_dict[beam_label]['pion_theta'].append(pi_theta)
                                data_dict[beam_label]['pion_phi'].append(pi_phi)

    # 3) Original plots for Lambda endpoints

    # We'll define range=0..40m
    zrange_m = (0, 40.0)
    nbins = 50

    ###############################
    # (a) Single overlay
    ###############################
    plt.figure(figsize=(8, 6))
    for beam_label in ["5x41", "10x100", "18x275"]:
        if beam_label not in data_dict:
            continue
        zvals_mm = data_dict[beam_label]['endpoint_z']
        if len(zvals_mm) == 0:
            continue
        # convert mm->m
        zvals_m = [z/1000.0 for z in zvals_mm]

        plt.hist(zvals_m,
                 bins=nbins,
                 range=zrange_m,
                 histtype="stepfilled",
                 alpha=0.5,
                 facecolor=color_map.get(beam_label, "gray"),
                 label=beam_label)

    plt.xlabel("Lambda decay endpoint z [m]")
    plt.ylabel("Events")
    plt.title("All Lambdas (PDG=3122), 0-40 m range (overlay)")
    plt.legend()

    out_pdf = os.path.join(args.outdir, "Lambda_endpoint_overlay.pdf")
    out_png = os.path.join(args.outdir, "Lambda_endpoint_overlay.png")
    plt.savefig(out_pdf)
    plt.savefig(out_png)
    plt.close()
    print(f"Saved overlay => {out_pdf} / {out_png}")

    ###############################
    # (b) Vertical subplots for each beam
    ###############################
    beams_order = ["5x41", "10x100", "18x275"]
    fig, axs = plt.subplots(len(beams_order), 1, figsize=(8, 12), sharex=True)

    for i, beam_label in enumerate(beams_order):
        ax = axs[i]
        ax.set_title(f"Beam: {beam_label}")
        ax.set_xlim(zrange_m)
        ax.set_ylabel("Events")
        zvals_mm = data_dict.get(beam_label, {}).get('endpoint_z', [])
        if len(zvals_mm) == 0:
            ax.text(0.5, 0.5, "No data", ha='center', va='center', transform=ax.transAxes)
            continue
        zvals_m = [z/1000.0 for z in zvals_mm]
        ax.hist(zvals_m, bins=nbins, range=zrange_m,
                histtype="stepfilled", alpha=0.5, facecolor=color_map.get(beam_label, "gray"))

    axs[-1].set_xlabel("Lambda decay endpoint z [m]")
    fig.tight_layout()

    out_pdf = os.path.join(args.outdir, "Lambda_endpoint_subplots.pdf")
    out_png = os.path.join(args.outdir, "Lambda_endpoint_subplots.png")
    plt.savefig(out_pdf)
    plt.savefig(out_png)
    plt.close()
    print(f"Saved vertical subplots => {out_pdf} / {out_png}")

    ###############################
    # (c) Individual single-plot per beam
    ###############################
    for beam_label in beams_order:
        zvals_mm = data_dict.get(beam_label, {}).get('endpoint_z', [])
        if len(zvals_mm) == 0:
            print(f"No data for beam={beam_label}, skipping single-plot.")
            continue

        zvals_m = [z/1000.0 for z in zvals_mm]
        plt.figure(figsize=(8,6))
        plt.hist(zvals_m,
                 bins=nbins,
                 range=zrange_m,
                 histtype="stepfilled",
                 alpha=0.5,
                 facecolor=color_map.get(beam_label,"gray"))
        plt.xlabel("Lambda decay endpoint z [m]")
        plt.ylabel("Events")
        plt.title(f"All Lambdas (PDG=3122), 0-40 m range\nBeam: {beam_label}")

        out_pdf = os.path.join(args.outdir, f"Lambda_endpoint_{beam_label}.pdf")
        out_png = os.path.join(args.outdir, f"Lambda_endpoint_{beam_label}.png")
        plt.savefig(out_pdf)
        plt.savefig(out_png)
        plt.close()
        print(f"Saved single-plot => {out_pdf} / {out_png}")

    ###############################
    # 4) Theta-Phi distributions
    ###############################

    # We'll create the angular distributions in a grid
    n_rows = len(beams_order)
    n_cols = 4  # proton theta, proton phi, pion theta, pion phi
    fig = plt.figure(figsize=(16, 4 * n_rows))
    gs = GridSpec(n_rows, n_cols, figure=fig)

    for row, beam_label in enumerate(beams_order):
        if beam_label not in data_dict:
            continue
        beam_data = data_dict[beam_label]

        # Check if we have any proton/pion angles for this beam
        if len(beam_data['proton_theta']) == 0:
            continue

        # 1. Proton theta
        ax1 = fig.add_subplot(gs[row, 0])
        ax1.hist(beam_data['proton_theta'], bins=36, range=(0, 180),
                 histtype='stepfilled', alpha=0.7, color=color_map.get(beam_label, 'gray'))
        ax1.set_xlabel(r"Proton $\theta$ [degrees]")
        ax1.set_ylabel("Count")
        if row == 0:
            ax1.set_title(r"Proton $\theta$ Distribution")
        ax1.text(0.05, 0.95, f"Beam: {beam_label}", transform=ax1.transAxes,
                 verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

        # 2. Proton phi
        ax2 = fig.add_subplot(gs[row, 1])
        ax2.hist(beam_data['proton_phi'], bins=36, range=(-180, 180),
                 histtype='stepfilled', alpha=0.7, color=color_map.get(beam_label, 'gray'))
        ax2.set_xlabel(r"Proton $\phi$ [degrees]")
        if row == 0:
            ax2.set_title(r"Proton $\phi$ Distribution")
        ax2.text(0.05, 0.95, f"Beam: {beam_label}", transform=ax2.transAxes,
                 verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

        # 3. Pion theta
        ax3 = fig.add_subplot(gs[row, 2])
        ax3.hist(beam_data['pion_theta'], bins=36, range=(0, 180),
                 histtype='stepfilled', alpha=0.7, color=color_map.get(beam_label, 'gray'))
        ax3.set_xlabel(r"Pion $\theta$ [degrees]")
        if row == 0:
            ax3.set_title(r"Pion $\theta$ Distribution")
        ax3.text(0.05, 0.95, f"Beam: {beam_label}", transform=ax3.transAxes,
                 verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

        # 4. Pion phi
        ax4 = fig.add_subplot(gs[row, 3])
        ax4.hist(beam_data['pion_phi'], bins=36, range=(-180, 180),
                 histtype='stepfilled', alpha=0.7, color=color_map.get(beam_label, 'gray'))
        ax4.set_xlabel(r"Pion $\phi$ [degrees]")
        if row == 0:
            ax4.set_title(r"Pion $\phi$ Distribution")
        ax4.text(0.05, 0.95, f"Beam: {beam_label}", transform=ax4.transAxes,
                 verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

    fig.tight_layout()
    angular_pdf = os.path.join(args.outdir, "Lambda_decay_angular_distributions.pdf")
    angular_png = os.path.join(args.outdir, "Lambda_decay_angular_distributions.png")
    plt.savefig(angular_pdf)
    plt.savefig(angular_png)
    plt.close()
    print(f"Saved angular distributions => {angular_pdf} / {angular_png}")

    ###############################
    # 5) 2D Theta-Phi scatter plots
    ###############################

    fig = plt.figure(figsize=(12, 4 * n_rows))
    gs = GridSpec(n_rows, 2, figure=fig)

    for row, beam_label in enumerate(beams_order):
        if beam_label not in data_dict:
            continue

        beam_data = data_dict[beam_label]

        if len(beam_data['proton_theta']) == 0:
            continue

        # Proton (theta vs phi)
        ax1 = fig.add_subplot(gs[row, 0])
        ax1.scatter(beam_data['proton_phi'], beam_data['proton_theta'],
                    alpha=0.5, c=color_map.get(beam_label, 'gray'), s=5)
        ax1.set_xlabel(r"Proton $\phi$ [degrees]")
        ax1.set_ylabel(r"Proton $\theta$ [degrees]")
        ax1.set_xlim(-180, 180)
        ax1.set_ylim(0, 180)
        if row == 0:
            ax1.set_title(r"Proton $\theta$ vs $\phi$")
        ax1.text(0.05, 0.95, f"Beam: {beam_label}", transform=ax1.transAxes,
                 verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

        # Pion (theta vs phi)
        ax2 = fig.add_subplot(gs[row, 1])
        ax2.scatter(beam_data['pion_phi'], beam_data['pion_theta'],
                    alpha=0.5, c=color_map.get(beam_label, 'gray'), s=5)
        ax2.set_xlabel(r"Pion $\phi$ [degrees]")
        ax2.set_ylabel(r"Pion $\theta$ [degrees]")
        ax2.set_xlim(-180, 180)
        ax2.set_ylim(0, 180)
        if row == 0:
            ax2.set_title(r"Pion $\theta$ vs $\phi$")
        ax2.text(0.05, 0.95, f"Beam: {beam_label}", transform=ax2.transAxes,
                 verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

    fig.tight_layout()
    angular_2d_pdf = os.path.join(args.outdir, "Lambda_decay_angular_2d_scatter.pdf")
    angular_2d_png = os.path.join(args.outdir, "Lambda_decay_angular_2d_scatter.png")
    plt.savefig(angular_2d_pdf)
    plt.savefig(angular_2d_png)
    plt.close()
    print(f"Saved 2D angular scatter plots => {angular_2d_pdf} / {angular_2d_png}")

    # Print some summary statistics
    for beam_label in beams_order:
        if beam_label not in data_dict:
            continue

        beam_data = data_dict[beam_label]
        n_lambdas = len(beam_data['endpoint_z'])
        n_decay_products = len(beam_data['proton_theta'])

        if n_lambdas > 0:
            fraction = 100.0 * n_decay_products / n_lambdas
        else:
            fraction = 0.0

        print(f"\nBeam {beam_label} statistics:")
        print(f"  Total Lambdas: {n_lambdas}")
        print(f"  Lambdas with (p + pi-) decay: {n_decay_products} ({fraction:.1f}% of total)")

if __name__ == "__main__":
    main()
