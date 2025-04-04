"""
This analysis focuses mainly on MCParticles,
which are basically true generator particles + whatever Geant4 adds: decay products, scattering etc.

"""


import argparse
import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Any

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
    parser = argparse.ArgumentParser(description="All-in-one script that correctly uses _MCParticles_daughters.index to find child PDGs.")
    parser.add_argument("input-files", nargs="+", help="List of EDM4hep ROOT files containing MCParticles.")
    parser.add_argument("-t", "--tree-name", default="events", help="Name of the TTree (default 'events').")
    parser.add_argument("-o", "--outdir", default="plots", help="Output directory for PDF/PNG plots.")
    parser.add_argument("--step-size", default="1000", type=int, help="Number of events to read per chunk (integer).")
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

    Note: Fixed phi calculation to avoid 180-degree rotation.
    """
    p = np.sqrt(px**2 + py**2 + pz**2)
    if p == 0:
        return 0.0, 0.0

    theta = np.arccos(pz / p) * 180.0 / np.pi
    # Original: phi = np.arctan2(py, px) * 180.0 / np.pi
    # Using arctan2 already gives the correct quadrant, no need to rotate
    phi = np.arctan2(py, px) * 180.0 / np.pi

    return theta, phi


@dataclass
class ParticleData:
    """Data container for particle information collected from all chunks"""
    # Lambda decay vertex
    lambda_z: List[float] = field(default_factory=list)
    # Decay product angles
    proton_theta: List[float] = field(default_factory=list)
    proton_phi: List[float] = field(default_factory=list)
    pion_theta: List[float] = field(default_factory=list)
    pion_phi: List[float] = field(default_factory=list)
    # Momentum components for different particle types
    pz_lambda: List[float] = field(default_factory=list)
    pt_lambda: List[float] = field(default_factory=list)
    pz_proton: List[float] = field(default_factory=list)
    pt_proton: List[float] = field(default_factory=list)
    pz_neutron: List[float] = field(default_factory=list)
    pt_neutron: List[float] = field(default_factory=list)
    pz_pion: List[float] = field(default_factory=list)
    pt_pion: List[float] = field(default_factory=list)
    # Debug counters
    counters: Dict[str, int] = field(default_factory=lambda: {
        "all_lambdas": 0,
        "valid_dtrs": 0,
        "has_p_pim": 0,
        "filled_angles": 0
    })


def process_chunk(chunk, data: ParticleData):
    """
    Process one chunk of data, event-by-event

    Args:
        chunk: Dict of awkward arrays for this chunk of events
        data: ParticleData instance to populate
    """
    pdg = chunk["MCParticles.PDG"]
    px = chunk["MCParticles.momentum.x"]
    py = chunk["MCParticles.momentum.y"]
    pz = chunk["MCParticles.momentum.z"]
    endpoint_z = chunk["MCParticles.endpoint.z"]
    d_beg = chunk["MCParticles.daughters_begin"]
    d_end = chunk["MCParticles.daughters_end"]
    daughters_idx_array = chunk["_MCParticles_daughters.index"]

    n_events = len(pdg)  # number of events in this chunk

    for i_evt in range(n_events):
        evt_pdg = pdg[i_evt]       # shape (nParticlesInEvent,)
        evt_px = px[i_evt]
        evt_py = py[i_evt]
        evt_pz = pz[i_evt]
        evt_z = endpoint_z[i_evt]

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
            px_ = evt_px[i_part]
            py_ = evt_py[i_part]
            pz_ = evt_pz[i_part]
            pt_ = np.sqrt(px_**2 + py_**2)

            if pdg_ == PDG_LAMBDA:
                data.pz_lambda.append(pz_)
                data.pt_lambda.append(pt_)
            elif pdg_ == PDG_PROTON:
                data.pz_proton.append(pz_)
                data.pt_proton.append(pt_)
            elif pdg_ == PDG_NEUTRON:
                data.pz_neutron.append(pz_)
                data.pt_neutron.append(pt_)
            elif pdg_ in [PDG_PION_PLUS, PDG_PION_MINUS]:
                data.pz_pion.append(pz_)
                data.pt_pion.append(pt_)

        # Identify Lambdas in this event
        lambda_indices = np.where(evt_pdg == PDG_LAMBDA)[0]
        for lam_idx in lambda_indices:
            data.counters["all_lambdas"] += 1

            lam_z_mm = evt_z[lam_idx]
            data.lambda_z.append(lam_z_mm / 1000.0)  # mm->m

            d_b = evt_db[lam_idx]
            d_e = evt_de[lam_idx]
            # If invalid, skip
            if (d_b < 0) or (d_e <= d_b):
                continue

            data.counters["valid_dtrs"] += 1

            # Correct approach: slice the event_dau_idx array to get child indices
            child_indices = evt_dau_idx[d_b: d_e]
            # child_indices is an array of ints referencing the main MCParticles array
            # Now get the child PDGs
            child_pdgs = evt_pdg[child_indices]

            # If we want to find p + pi- among children, let's check child_pdgs
            if (PDG_PROTON in child_pdgs) and (PDG_PION_MINUS in child_pdgs):
                data.counters["has_p_pim"] += 1
                # find which child is the proton, which is the pi-
                local_proton_idxs = np.where(child_pdgs == PDG_PROTON)[0]
                local_pionm_idxs = np.where(child_pdgs == PDG_PION_MINUS)[0]

                if len(local_proton_idxs) > 0 and len(local_pionm_idxs) > 0:
                    # get the actual index in the main arrays
                    p_idx = child_indices[local_proton_idxs[0]]
                    pim_idx = child_indices[local_pionm_idxs[0]]

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

                    data.proton_theta.append(p_theta)
                    data.proton_phi.append(p_phi)
                    data.pion_theta.append(pi_theta)
                    data.pion_phi.append(pi_phi)

                    data.counters["filled_angles"] += 1


def analyze_beam(beam_label, file_list, tree_name, step_size, outdir):
    """
    For the given beam_label and list of files:
      1) create a ParticleData container
      2) read each file chunk-by-chunk
      3) fill the container with data
      4) after reading, produce plots + debug output
      5) return the Lambda endpoints for final overlay
    """
    # Create data container for this beam
    data = ParticleData()

    # Read files chunk-by-chunk and process
    for fname in file_list:
        file_dict = {fname: tree_name}
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
                    "_MCParticles_daughters.index",
                ],
                step_size=step_size):
            # Process this chunk and update our data container
            process_chunk(chunk, data)

    # === after processing all chunks for this beam, produce plots ===
    os.makedirs(outdir, exist_ok=True)

    # 1) Lambda endpoint
    zrange = (0, 40)
    plt.figure(figsize=(8, 6))
    plt.hist(data.lambda_z, bins=50, range=zrange,
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
    if len(data.proton_theta) > 0:
        fig, axs = plt.subplots(2, 2, figsize=(10, 8))

        axs[0, 0].hist(data.proton_theta, bins=36, range=(0, 180), histtype='stepfilled', alpha=0.7)
        axs[0, 0].set_xlabel(r"Proton $\theta$ [deg]")
        axs[0, 0].set_ylabel("Count")

        axs[0, 1].hist(data.proton_phi, bins=36, range=(-180, 180), histtype='stepfilled', alpha=0.7)
        axs[0, 1].set_xlabel(r"Proton $\phi$ [deg]")
        axs[0, 1].set_ylabel("Count")

        axs[1, 0].hist(data.pion_theta, bins=36, range=(0, 180), histtype='stepfilled', alpha=0.7)
        axs[1, 0].set_xlabel(r"Pion$^-$ $\theta$ [deg]")
        axs[1, 0].set_ylabel("Count")

        axs[1, 1].hist(data.pion_phi, bins=36, range=(-180, 180), histtype='stepfilled', alpha=0.7)
        axs[1, 1].set_xlabel(r"Pion$^-$ $\phi$ [deg]")
        axs[1, 1].set_ylabel("Count")

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
    fig, axs = plt.subplots(2, 2, figsize=(10, 8))
    axs = axs.flatten()

    axs[0].scatter(data.pz_lambda, data.pt_lambda, s=5, alpha=0.5, color='green')
    axs[0].set_xlabel(r"$p_z$ [GeV/c]")
    axs[0].set_ylabel(r"$p_T$ [GeV/c]")
    axs[0].set_title("Lambdas")

    axs[1].scatter(data.pz_proton, data.pt_proton, s=5, alpha=0.5, color='red')
    axs[1].set_xlabel(r"$p_z$ [GeV/c]")
    axs[1].set_ylabel(r"$p_T$ [GeV/c]")
    axs[1].set_title("Protons")

    axs[2].scatter(data.pz_neutron, data.pt_neutron, s=5, alpha=0.5, color='blue')
    axs[2].set_xlabel(r"$p_z$ [GeV/c]")
    axs[2].set_ylabel(r"$p_T$ [GeV/c]")
    axs[2].set_title("Neutrons")

    axs[3].scatter(data.pz_pion, data.pt_pion, s=5, alpha=0.5, color='magenta')
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
    print(f"\n[DEBUG: {beam_label}]")
    print(f"   Found Lambdas: {data.counters['all_lambdas']}")
    print(f"   Lambdas with valid daughter indices: {data.counters['valid_dtrs']}")
    print(f"   Lambdas with (p + pi-): {data.counters['has_p_pim']}")
    print(f"   # times we actually filled angles: {data.counters['filled_angles']}")

    return data.lambda_z


def make_final_comparison(lambda_endpoints_dict, beam_labels, outdir):
    """
    Create an overlay / comparison histogram of Lambda endpoints for multiple beams.
    """
    os.makedirs(outdir, exist_ok=True)
    colors = ["red", "blue", "green", "magenta", "cyan"]

    plt.figure(figsize=(8, 6))
    nbins = 50
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
        z_vals_m = analyze_beam(beam_label, file_list, args.tree_name, args.step_size, args.outdir)
        lambda_endpoints_dict[beam_label] = z_vals_m

    # Make final overlay across user-specified beam labels
    print("\n=== Final comparison overlay of Lambda endpoints ===")
    make_final_comparison(lambda_endpoints_dict, args.beam_labels, args.outdir)


if __name__ == "__main__":
    main()