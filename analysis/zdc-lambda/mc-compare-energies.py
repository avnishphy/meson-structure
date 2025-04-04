#!/usr/bin/env python3

import argparse
import os
import re

import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt

# PDG codes
PDG_LAMBDA = 3122
PDG_PROTON = 2212
PDG_NEUTRON = 2112
PDG_PION_PLUS = 211
PDG_PION_MINUS = -211


def parse_args():
    """ Parses command line arguments"""
    parser = argparse.ArgumentParser(description="All-in-one script that correctly uses _MCParticles_daughters.index to find child PDGs, chunk by chunk.")
    parser.add_argument("-i", "--input-files", nargs="+", required=True, help="List of EDM4hep ROOT files containing MCParticles.")
    parser.add_argument("-t", "--tree-name", default="events", help="Name of the TTree (default 'events').")
    parser.add_argument("-o", "--outdir", default="plots", help="Output directory for PDF/PNG plots.")
    parser.add_argument("--step-size", default="1000", help="Number of events to read per chunk (integer).")
    parser.add_argument("--beam-labels", nargs="+",
                        default=["5x41", "10x100", "18x275"],
                        help="Which beam labels to compare in final overlay.")
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


def calculate_theta_phi(px, py, pz):
    """
    Calculate theta (0..180 deg) and phi (-180..180 deg)
    from momentum components.
    """
    p = np.sqrt(px**2 + py**2 + pz**2)
    if p == 0:
        return 0.0, 0.0
    theta = np.arccos(pz / p) * 180.0 / np.pi
    phi   = np.arctan2(py, px) * 180.0 / np.pi
    return theta, phi


def process_chunk(
        pdg, px, py, pz, endpoint_z,
        d_beg, d_end,
        daughters_idx_array,
        # Arrays for final results:
        lambda_z_vals_m,
        proton_theta_vals, proton_phi_vals,
        pion_theta_vals,   pion_phi_vals,
        pz_lambdas,        pt_lambdas,
        pz_protons,        pt_protons,
        pz_neutrons,       pt_neutrons,
        pz_pions,          pt_pions,
        debug_counters
):
    """
    Process one chunk of data, event-by-event, using the correct approach:
      - d_beg, d_end are offsets into the _MCParticles_daughters.index array
      - For each particle, we slice daughters_idx_array[ d_beg[i] : d_end[i] ]
        to get the actual child indices in the main array
    """

    n_events = len(pdg)  # number of events in this chunk

    for i_evt in range(n_events):
        evt_pdg = pdg[i_evt]       # shape (nParticlesInEvent,)
        evt_px  = px[i_evt]
        evt_py  = py[i_evt]
        evt_pz  = pz[i_evt]
        evt_z   = endpoint_z[i_evt]

        # The arrays for daughters_begin, daughters_end in this event
        evt_db = d_beg[i_evt]      # shape (nParticlesInEvent,)
        evt_de = d_end[i_evt]      # shape (nParticlesInEvent,)

        # The event-level _MCParticles_daughters.index array
        # shape: (nDaughtersForThisEvent,)
        evt_dau_idx = daughters_idx_array[i_evt]

        # First, fill pz, pT for all relevant species in this event
        n_parts = len(evt_pdg)
        for i_part in range(n_parts):
            pdg_ = evt_pdg[i_part]
            px_  = evt_px[i_part]
            py_  = evt_py[i_part]
            pz_  = evt_pz[i_part]
            pt_  = np.sqrt(px_**2 + py_**2)

            if pdg_ == PDG_LAMBDA:
                pz_lambdas.append(pz_)
                pt_lambdas.append(pt_)
            elif pdg_ == PDG_PROTON:
                pz_protons.append(pz_)
                pt_protons.append(pt_)
            elif pdg_ == PDG_NEUTRON:
                pz_neutrons.append(pz_)
                pt_neutrons.append(pt_)
            elif pdg_ in [PDG_PION_PLUS, PDG_PION_MINUS]:
                pz_pions.append(pz_)
                pt_pions.append(pt_)

        # Identify Lambdas in this event
        lambda_indices = np.where(evt_pdg == PDG_LAMBDA)[0]
        for lam_idx in lambda_indices:
            debug_counters["all_lambdas"] += 1

            lam_z_mm = evt_z[lam_idx]
            lambda_z_vals_m.append(lam_z_mm / 1000.0)  # mm->m

            d_b = evt_db[lam_idx]
            d_e = evt_de[lam_idx]
            # If invalid, skip
            if (d_b < 0) or (d_e <= d_b):
                continue

            debug_counters["valid_dtrs"] += 1

            # Correct approach: slice the event_dau_idx array to get child indices
            child_indices = evt_dau_idx[d_b : d_e]
            # child_indices is an array of ints referencing the main MCParticles array
            # Now get the child PDGs
            child_pdgs = evt_pdg[child_indices]

            # If we want to find p + pi- among children, let's check child_pdgs
            if (PDG_PROTON in child_pdgs) and (PDG_PION_MINUS in child_pdgs):
                debug_counters["has_p_pim"] += 1
                # find which child is the proton, which is the pi-
                local_proton_idxs = np.where(child_pdgs == PDG_PROTON)[0]
                local_pionm_idxs  = np.where(child_pdgs == PDG_PION_MINUS)[0]

                if len(local_proton_idxs) > 0 and len(local_pionm_idxs) > 0:
                    # get the actual index in the main arrays
                    p_idx   = child_indices[ local_proton_idxs[0] ]
                    pim_idx = child_indices[ local_pionm_idxs[0] ]

                    # retrieve momentum from the main arrays
                    p_px = evt_px[p_idx]
                    p_py = evt_py[p_idx]
                    p_pz = evt_pz[p_idx]

                    pi_px = evt_px[pim_idx]
                    pi_py = evt_py[pim_idx]
                    pi_pz = evt_pz[pim_idx]

                    # angles
                    p_theta, p_phi = calculate_theta_phi(p_px, p_py, p_pz)
                    pi_theta, pi_phi = calculate_theta_phi(pi_px, pi_py, pi_pz)

                    proton_theta_vals.append(p_theta)
                    proton_phi_vals.append(p_phi)
                    pion_theta_vals.append(pi_theta)
                    pion_phi_vals.append(pi_phi)

                    debug_counters["filled_angles"] += 1

def analyze_beam(beam_label, file_list, tree_name, step_size, outdir):
    """
    For the given beam_label and list of files:
      1) create final arrays
      2) read each file chunk-by-chunk
      3) fill arrays & debug counters
      4) after reading, produce plots + debug output
      5) return the Lambda endpoints for final overlay
    """

    # Arrays that accumulate data from all chunks
    lambda_z_vals_m = []
    proton_theta_vals = []
    proton_phi_vals   = []
    pion_theta_vals   = []
    pion_phi_vals     = []

    pz_lambdas  = []
    pt_lambdas  = []
    pz_protons  = []
    pt_protons  = []
    pz_neutrons = []
    pt_neutrons = []
    pz_pions    = []
    pt_pions    = []

    # Debug counters
    debug_counters = {
        "all_lambdas": 0,
        "valid_dtrs": 0,
        "has_p_pim": 0,
        "filled_angles": 0
    }

    # ---- chunk-based reading & immediate processing ----
    for fname in file_list:
        file_dict = {fname: tree_name}
        # We must also read the _MCParticles_daughters.index array
        for chunk in uproot.iterate(
                file_dict,
                expressions=[
                    "MCParticles.PDG",
                    "MCParticles.momentum.x",
                    "MCParticles.momentum.y",
                    "MCParticles.momentum.z",
                    "MCParticles.endpoint.z",
                    "MCParticles.daughters_begin",
                    "MCParticles.daughters_end",
                    # The array that stores which index in MCParticles is each daughter
                    "_MCParticles_daughters.index",
                ],
                step_size=step_size,
                axis="entry",
                library="ak"):
            # chunk is an Awkward array dict
            # We'll pass chunk["_MCParticles_daughters.index"] into the process function
            # so it can do the correct indexing
            process_chunk(
                pdg=chunk["MCParticles.PDG"],
                px=chunk["MCParticles.momentum.x"],
                py=chunk["MCParticles.momentum.y"],
                pz=chunk["MCParticles.momentum.z"],
                endpoint_z=chunk["MCParticles.endpoint.z"],
                d_beg=chunk["MCParticles.daughters_begin"],
                d_end=chunk["MCParticles.daughters_end"],
                daughters_idx_array=chunk["_MCParticles_daughters.index"],

                # pass the final arrays and debug counters to fill
                lambda_z_vals_m=lambda_z_vals_m,
                proton_theta_vals=proton_theta_vals,
                proton_phi_vals=proton_phi_vals,
                pion_theta_vals=pion_theta_vals,
                pion_phi_vals=pion_phi_vals,
                pz_lambdas=pz_lambdas,
                pt_lambdas=pt_lambdas,
                pz_protons=pz_protons,
                pt_protons=pt_protons,
                pz_neutrons=pz_neutrons,
                pt_neutrons=pt_neutrons,
                pz_pions=pz_pions,
                pt_pions=pt_pions,
                debug_counters=debug_counters
            )

    # === after processing all chunks for this beam, produce plots ===
    os.makedirs(outdir, exist_ok=True)

    # 1) Lambda endpoint
    zrange = (0, 40)
    plt.figure(figsize=(8,6))
    plt.hist(lambda_z_vals_m, bins=50, range=zrange,
             histtype='stepfilled', alpha=0.5, color='blue')
    plt.xlabel("Lambda decay endpoint z [m]")
    plt.ylabel("Count")
    plt.title(f"Lambda endpoints (0-40 m)\nBeam: {beam_label}")
    out_pdf = os.path.join(outdir, f"Lambda_endpoint_{beam_label}.pdf")
    out_png = os.path.join(outdir, f"Lambda_endpoint_{beam_label}.png")
    plt.savefig(out_pdf)
    plt.savefig(out_png)
    plt.close()

    # 2) Angular distributions
    if len(proton_theta_vals) > 0:
        fig, axs = plt.subplots(2,2, figsize=(10, 8))

        axs[0,0].hist(proton_theta_vals, bins=36, range=(0, 180), histtype='stepfilled', alpha=0.7)
        axs[0,0].set_xlabel(r"Proton $\theta$ [deg]")
        axs[0,0].set_ylabel("Count")

        axs[0,1].hist(proton_phi_vals, bins=36, range=(-180,180), histtype='stepfilled', alpha=0.7)
        axs[0,1].set_xlabel(r"Proton $\phi$ [deg]")
        axs[0,1].set_ylabel("Count")

        axs[1,0].hist(pion_theta_vals, bins=36, range=(0,180), histtype='stepfilled', alpha=0.7)
        axs[1,0].set_xlabel(r"Pion$^-$ $\theta$ [deg]")
        axs[1,0].set_ylabel("Count")

        axs[1,1].hist(pion_phi_vals, bins=36, range=(-180,180), histtype='stepfilled', alpha=0.7)
        axs[1,1].set_xlabel(r"Pion$^-$ $\phi$ [deg]")
        axs[1,1].set_ylabel("Count")

        fig.suptitle(f"Lambda->p + pi- angles\nBeam: {beam_label}")
        fig.tight_layout()

        out_pdf = os.path.join(outdir, f"Lambda_decay_angles_{beam_label}.pdf")
        out_png = os.path.join(outdir, f"Lambda_decay_angles_{beam_label}.png")
        plt.savefig(out_pdf)
        plt.savefig(out_png)
        plt.close()
    else:
        print(f"  [INFO] No (p+pi-) decays found for beam={beam_label}, skipping angle plots.")

    # 3) pz vs pT for Lambdas, p, n, pi
    fig, axs = plt.subplots(2,2, figsize=(10,8))
    axs = axs.flatten()

    axs[0].scatter(pz_lambdas, pt_lambdas, s=5, alpha=0.5, color='green')
    axs[0].set_xlabel(r"$p_z$ [GeV/c]")
    axs[0].set_ylabel(r"$p_T$ [GeV/c]")
    axs[0].set_title("Lambdas")

    axs[1].scatter(pz_protons, pt_protons, s=5, alpha=0.5, color='red')
    axs[1].set_xlabel(r"$p_z$ [GeV/c]")
    axs[1].set_ylabel(r"$p_T$ [GeV/c]")
    axs[1].set_title("Protons")

    axs[2].scatter(pz_neutrons, pt_neutrons, s=5, alpha=0.5, color='blue')
    axs[2].set_xlabel(r"$p_z$ [GeV/c]")
    axs[2].set_ylabel(r"$p_T$ [GeV/c]")
    axs[2].set_title("Neutrons")

    axs[3].scatter(pz_pions, pt_pions, s=5, alpha=0.5, color='magenta')
    axs[3].set_xlabel(r"$p_z$ [GeV/c]")
    axs[3].set_ylabel(r"$p_T$ [GeV/c]")
    axs[3].set_title("Pions (±)")

    fig.suptitle(f"pZ vs pT\nBeam: {beam_label}")
    fig.tight_layout()

    out_pdf = os.path.join(outdir, f"pz_pt_{beam_label}.pdf")
    out_png = os.path.join(outdir, f"pz_pt_{beam_label}.png")
    plt.savefig(out_pdf)
    plt.savefig(out_png)
    plt.close()

    # ========== DEBUG OUTPUT ==========
    n_lambdas         = debug_counters["all_lambdas"]
    n_valid_dtrs      = debug_counters["valid_dtrs"]
    n_has_p_pim       = debug_counters["has_p_pim"]
    n_filled_angles   = debug_counters["filled_angles"]
    print(f"\n[DEBUG: {beam_label}]")
    print(f"   Found Lambdas: {n_lambdas}")
    print(f"   Lambdas with valid daughter indices: {n_valid_dtrs}")
    print(f"   Lambdas with (p + pi-): {n_has_p_pim}")
    print(f"   # times we actually filled angles: {n_filled_angles}")

    return lambda_z_vals_m


def make_final_comparison(lambda_endpoints_dict, beam_labels, outdir):
    """
    Create an overlay / comparison histogram of Lambda endpoints for multiple beams.
    """
    os.makedirs(outdir, exist_ok=True)
    colors = ["red","blue","green","magenta","cyan"]

    plt.figure(figsize=(8,6))
    nbins  = 50
    zrange = (0, 40)

    for i, beam_label in enumerate(beam_labels):
        if beam_label not in lambda_endpoints_dict:
            print(f"  [WARN] No data for beam={beam_label} => skipping overlay.")
            continue

        z_m = lambda_endpoints_dict[beam_label]
        plt.hist(z_m, bins=nbins, range=zrange,
                 histtype='stepfilled',
                 alpha=0.5, color=colors[i % len(colors)],
                 label=beam_label)

    plt.xlabel("Lambda decay endpoint z [m]")
    plt.ylabel("Count")
    plt.title("Overlay: Lambda endpoints across beams")
    plt.legend()

    out_pdf = os.path.join(outdir, "Lambda_endpoint_overlay.pdf")
    out_png = os.path.join(outdir, "Lambda_endpoint_overlay.png")
    plt.savefig(out_pdf)
    plt.savefig(out_png)
    plt.close()

    print(f"[INFO] Final endpoint overlay saved to:\n  {out_pdf}\n  {out_png}")


def main():
    args = parse_args()
    step_size = int(args.step_size)
    os.makedirs(args.outdir, exist_ok=True)

    # Group files by beam
    beams_map = {}
    for fname in args.input_files:
        beam = extract_beam_label(fname)
        beams_map.setdefault(beam, []).append(fname)

    # We'll store final endpoints for each beam in a dict
    lambda_endpoints_dict = {}

    # Analyze beam by beam
    for beam_label, file_list in beams_map.items():
        print(f"\n=== Analyzing beam: {beam_label} ===")
        z_vals_m = analyze_beam(beam_label, file_list, args.tree_name, step_size, args.outdir)
        lambda_endpoints_dict[beam_label] = z_vals_m

    # Make final overlay across user-specified beam labels
    print("\n=== Final comparison overlay of Lambda endpoints ===")
    make_final_comparison(lambda_endpoints_dict, args.beam_labels, args.outdir)

if __name__ == "__main__":
    main()
