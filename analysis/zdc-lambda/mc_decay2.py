#!/usr/bin/env python3

import argparse
import re
import os
import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
import rich

def parse_args():
    parser = argparse.ArgumentParser(
        description="Identify Lambda decays (p+pi- or n+pi0) in EDM4hep, iterate with uproot, "
                    "and debug print how many are skipped or accepted."
    )
    parser.add_argument("-i", "--input-files", nargs="+", required=True,
                        help="ROOT files with MCParticles TTree (EDM4hep format).")
    parser.add_argument("-t", "--tree-name", default="events",
                        help="Name of the TTree (default: 'events').")
    parser.add_argument("-o", "--outdir", default="plots",
                        help="Directory for output PDF plots.")
    parser.add_argument("--step-size", default="100",
                        help="Chunk size for uproot.iterate (e.g. '100MB').")
    return parser.parse_args()


def extract_beam_label(filename):
    base = os.path.basename(filename)
    match = re.search(r'_(\d+x\d+)_', base)
    if match:
        return match.group(1)
    return "unknown"


# PDG codes
PDG_LAMBDA  = 3122
PDG_PROTON  = 2212
PDG_PIMINUS = -211
PDG_NEUTRON = 2112
PDG_PI0     = 111

def ensure_beam(data_dict, beam):
    if beam not in data_dict:
        data_dict[beam] = {"ppi-": [], "npi0": []}
    return data_dict[beam]

def analyze_chunk(chunk, data_dict, beam_label):
    """
    Process one chunk and accumulate endpoint.z for Lambda->p pi- or n pi0.
    Includes debug counters to show how many Lambdas we skip/accept.
    """

    pdg            = chunk["MCParticles.PDG"]
    daughters_begin= chunk["MCParticles.daughters_begin"]
    daughters_end  = chunk["MCParticles.daughters_end"]
    endpoint_z     = chunk["MCParticles.endpoint.z"]


    beam_dict = ensure_beam(data_dict, beam_label)

    # Debug counters
    total_lambdas     = 0
    skip_out_of_range = 0
    skip_not_2dau     = 0
    matched_ppi       = 0
    matched_npi0      = 0
    skip_other        = 0
    is_lambda = pdg == PDG_LAMBDA
    rich.print(is_lambda)
    fd = daughters_begin[is_lambda]
    ld = daughters_end[is_lambda]
    fdpdg = pdg[fd]
    ldpdg = pdg[ld]
    rich.print(fd)


    lambda_indices = ak.where(is_lambda)
    rich.print(is_lambda)
    exit(1)

    #
    # n_events = len(pdg)
    # for i_evt in range(n_events):
    #     # Which entries are Lambdas?
    #     is_lambda = (pdg[i_evt] == PDG_LAMBDA)
    #     lambda_indices = np.where(is_lambda)[0]
    #     total_lambdas += len(lambda_indices)
    #
    #     n_particles = len(pdg[i_evt])
    #
    #     for lam_idx in lambda_indices:
    #         d_begin = daughters_begin[i_evt][lam_idx]
    #         d_end   = daughters_end[i_evt][lam_idx]
    #
    #         # If negative or inverted range => skip
    #         if d_begin < 0 or (d_end < d_begin):
    #             skip_out_of_range += 1
    #             continue
    #
    #         if d_end > n_particles:
    #             # Out-of-bounds reference => skip
    #             skip_out_of_range += 1
    #             continue
    #
    #         n_daughters = d_end - d_begin
    #         if n_daughters != 2:
    #             skip_not_2dau += 1
    #             continue
    #
    #         d_idx1 = d_begin
    #         d_idx2 = d_begin + 1
    #         if d_idx2 >= n_particles:
    #             skip_out_of_range += 1
    #             continue
    #
    #         pdg_d1 = pdg[i_evt][d_idx1]
    #         pdg_d2 = pdg[i_evt][d_idx2]
    #
    #         daughter_set = frozenset({pdg_d1, pdg_d2})
    #
    #         if daughter_set == {PDG_PROTON, PDG_PIMINUS}:
    #             beam_dict["ppi-"].append(endpoint_z[i_evt][lam_idx])
    #             matched_ppi += 1
    #         elif daughter_set == {PDG_NEUTRON, PDG_PI0}:
    #             beam_dict["npi0"].append(endpoint_z[i_evt][lam_idx])
    #             matched_npi0 += 1
    #         else:
    #             skip_other += 1
    #
    # # Print debug for this chunk
    # print(f"[DEBUG] Beam={beam_label}  chunk stats:")
    # print(f"   total Lambdas found: {total_lambdas}")
    # print(f"   skip_out_of_range:   {skip_out_of_range}")
    # print(f"   skip_not_2dau:       {skip_not_2dau}")
    # print(f"   skip_other_daughters:{skip_other}")
    # print(f"   matched p+pi-:       {matched_ppi}")
    # print(f"   matched n+pi0:       {matched_npi0}")


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    data_dict = {}
    # Color map
    color_map = {
        "5x41": "red",
        "10x100": "blue",
        "18x275": "green",
    }

    # 1) Group files by beam
    beams_map = {}
    for f in args.input_files:
        beam = extract_beam_label(f)
        beams_map.setdefault(beam, []).append(f)

    # 2) Iterate
    for beam_label, file_list in beams_map.items():
        print(f"\n=== Beam: {beam_label} ===")
        for fname in file_list:
            print(f"   Reading file: {fname}")
            file_dict = {fname: args.tree_name}

            for chunk in uproot.iterate(
                    file_dict,
                    expressions=[
                        "MCParticles.PDG",
                        "MCParticles.daughters_begin",
                        "MCParticles.daughters_end",
                        "MCParticles.endpoint.z",
                    ],
                    step_size=int(args.step_size),
                    library="ak",
            ):
                analyze_chunk(chunk, data_dict, beam_label)

    # Now produce two overlaid histograms:
    #   1) p+pi- decays
    #   2) n+pi0 decays
    plt.figure(figsize=(8,6))
    for beam in ["5x41", "10x100", "18x275"]:
        if beam not in data_dict:
            continue
        zvals = data_dict[beam]["ppi-"]
        if len(zvals) == 0:
            continue
        plt.hist(zvals,
                 bins=50,
                 range=(0,2000),
                 histtype="stepfilled",
                 alpha=0.5,
                 facecolor=color_map.get(beam,"gray"),
                 label=f"{beam} (p+pi-)")
    plt.xlabel("Lambda decay endpoint.z [mm]")
    plt.ylabel("Events")
    plt.title("Lambda -> p + pi-  (Endpoint Z)")
    plt.legend()
    out_ppi = os.path.join(args.outdir, "decayZ_lambda_ppiminus.pdf")
    plt.savefig(out_ppi)
    print(f"Saved {out_ppi}")
    plt.close()

    plt.figure(figsize=(8,6))
    for beam in ["5x41", "10x100", "18x275"]:
        if beam not in data_dict:
            continue
        zvals = data_dict[beam]["npi0"]
        if len(zvals) == 0:
            continue
        plt.hist(zvals,
                 bins=50,
                 range=(0,2000),
                 histtype="stepfilled",
                 alpha=0.5,
                 facecolor=color_map.get(beam,"gray"),
                 label=f"{beam} (n+pi0)")
    plt.xlabel("Lambda decay endpoint.z [mm]")
    plt.ylabel("Events")
    plt.title("Lambda -> n + pi0  (Endpoint Z)")
    plt.legend()
    out_npi = os.path.join(args.outdir, "decayZ_lambda_npi0.pdf")
    plt.savefig(out_npi)
    print(f"Saved {out_npi}")
    plt.close()


if __name__ == "__main__":
    main()
