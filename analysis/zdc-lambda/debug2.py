#!/usr/bin/env python3

import uproot
import awkward as ak
import numpy as np
import vector
import matplotlib.pyplot as plt
import argparse
import os
import rich

def parse_args():
    parser = argparse.ArgumentParser(description="MC Truth Analysis of Lambda Particles")
    parser.add_argument("-i", "--input-file", required=True, help="Input ROOT file")
    parser.add_argument("-o", "--output", default="lambda_mc_analysis", help="Output directory")
    parser.add_argument("--tree", default="events", help="Tree name")
    return parser.parse_args()

def find_lambda_decays(data):
    """
    Find Lambda particles that decay to proton + pion in the MC truth.

    Parameters:
    -----------
    data : awkward array
        Data containing MCParticle information

    Returns:
    --------
    dict
        Dictionary containing Lambda info and their decay products
    """
    # Get PDG codes
    pdg = data["MCParticles.PDG"]

    # Find all Lambdas (PDG code 3122)
    lambda_mask = (pdg == 3122)

    # Get daughter indices
    daughters_begin = data["MCParticles.daughters_begin"]
    daughters_end = data["MCParticles.daughters_end"]

    # Initialize lists to store results
    lambda_info = {
        "idx": [],
        "mass": [],
        "px": [],
        "py": [],
        "pz": [],
        "proton_idx": [],
        "pion_idx": [],
        "proton_px": [],
        "proton_py": [],
        "proton_pz": [],
        "proton_mass": [],
        "pion_px": [],
        "pion_py": [],
        "pion_pz": [],
        "pion_mass": [],
        "reconstructed_mass": []
    }

    # Loop through events
    for evt_idx in range(len(pdg)):
        # Find all Lambdas in this event
        evt_lambda_indices = ak.where(lambda_mask[evt_idx])[0]

        for lambda_idx in evt_lambda_indices:
            # Get daughter range
            start_idx = daughters_begin[evt_idx][lambda_idx]
            end_idx = daughters_end[evt_idx][lambda_idx]

            # Skip if no daughters
            if end_idx <= start_idx:
                continue

            # Get daughter PDG codes
            daughter_pdgs = pdg[evt_idx][start_idx:end_idx]

            # Check if Lambda decays to proton + pion
            has_proton = (2212 in daughter_pdgs)
            has_pion = (-211 in daughter_pdgs)

            if has_proton and has_pion:
                # This Lambda decays to proton + pion
                proton_daughter_idx = ak.where(daughter_pdgs == 2212)[0][0] + start_idx
                pion_daughter_idx = ak.where(daughter_pdgs == -211)[0][0] + start_idx

                # Get Lambda properties
                lambda_mass = data["MCParticles.mass"][evt_idx][lambda_idx]
                lambda_px = data["MCParticles.momentum.x"][evt_idx][lambda_idx]
                lambda_py = data["MCParticles.momentum.y"][evt_idx][lambda_idx]
                lambda_pz = data["MCParticles.momentum.z"][evt_idx][lambda_idx]

                # Get proton properties
                proton_px = data["MCParticles.momentum.x"][evt_idx][proton_daughter_idx]
                proton_py = data["MCParticles.momentum.y"][evt_idx][proton_daughter_idx]
                proton_pz = data["MCParticles.momentum.z"][evt_idx][proton_daughter_idx]
                proton_mass = data["MCParticles.mass"][evt_idx][proton_daughter_idx]

                # Get pion properties
                pion_px = data["MCParticles.momentum.x"][evt_idx][pion_daughter_idx]
                pion_py = data["MCParticles.momentum.y"][evt_idx][pion_daughter_idx]
                pion_pz = data["MCParticles.momentum.z"][evt_idx][pion_daughter_idx]
                pion_mass = data["MCParticles.mass"][evt_idx][pion_daughter_idx]

                # Calculate reconstructed mass from decay products
                try:
                    # Calculate proton energy
                    proton_e = np.sqrt(proton_px**2 + proton_py**2 + proton_pz**2 + proton_mass**2)

                    # Calculate pion energy
                    pion_e = np.sqrt(pion_px**2 + pion_py**2 + pion_pz**2 + pion_mass**2)

                    # Create 4-vectors
                    proton_vec = vector.obj(px=proton_px, py=proton_py, pz=proton_pz, E=proton_e)
                    pion_vec = vector.obj(px=pion_px, py=pion_py, pz=pion_pz, E=pion_e)

                    # Add 4-vectors to get Lambda
                    reconstructed_vec = proton_vec + pion_vec
                    reconstructed_mass = reconstructed_vec.mass

                    # Store results
                    lambda_info["idx"].append(lambda_idx)
                    lambda_info["mass"].append(lambda_mass)
                    lambda_info["px"].append(lambda_px)
                    lambda_info["py"].append(lambda_py)
                    lambda_info["pz"].append(lambda_pz)
                    lambda_info["proton_idx"].append(proton_daughter_idx)
                    lambda_info["pion_idx"].append(pion_daughter_idx)
                    lambda_info["proton_px"].append(proton_px)
                    lambda_info["proton_py"].append(proton_py)
                    lambda_info["proton_pz"].append(proton_pz)
                    lambda_info["proton_mass"].append(proton_mass)
                    lambda_info["pion_px"].append(pion_px)
                    lambda_info["pion_py"].append(pion_py)
                    lambda_info["pion_pz"].append(pion_pz)
                    lambda_info["pion_mass"].append(pion_mass)
                    lambda_info["reconstructed_mass"].append(reconstructed_mass)

                except Exception as e:
                    print(f"Error reconstructing Lambda mass: {e}")

    return lambda_info


def find_reconstructed_matches(reco_data,reco_main_branch, mc_data, lambda_mc_info):
    """
    Try to find matches between reconstructed particles and MC truth.

    Parameters:
    -----------
    data : awkward array
        Data containing ReconstructedParticles
    lambda_mc_info : dict
        Dictionary of MC truth Lambda information

    Returns:
    --------
    dict
        Dictionary with matching information
    """
    # Create a dictionary to store matching information
    matches = {
        "lambda_idx": [],
        "mc_mass": [],
        "reco_mass": [],
        "proton_pt_mc": [],
        "proton_pt_reco": [],
        "pion_pt_mc": [],
        "pion_pt_reco": []
    }

    # Extract reconstructed particles
    reco_pdg =  reco_data[f"{reco_main_branch}.PDG"]
    reco_px =   reco_data[f"{reco_main_branch}.momentum.x"]
    reco_py =   reco_data[f"{reco_main_branch}.momentum.y"]
    reco_pz =   reco_data[f"{reco_main_branch}.momentum.z"]
    reco_mass = reco_data[f"{reco_main_branch}.mass"]

    # Get MC information arrays
    proton_px_mc = np.array(lambda_mc_info["proton_px"])
    proton_py_mc = np.array(lambda_mc_info["proton_py"])
    proton_pz_mc = np.array(lambda_mc_info["proton_pz"])
    pion_px_mc = np.array(lambda_mc_info["pion_px"])
    pion_py_mc = np.array(lambda_mc_info["pion_py"])
    pion_pz_mc = np.array(lambda_mc_info["pion_pz"])

    mc_pdgs = mc_data["MCParticles.PDG"]

    # Loop through events to find matches
    for evt_idx in range(len(reco_pdg)):

        evt_reco_pdgs = reco_pdg[evt_idx]
        evt_mc_pdgs = mc_pdgs[evt_idx]

        print(evt_reco_pdgs)
        print(evt_mc_pdgs)


        # Skip if no protons or pions in this event
        if 2212 not in reco_pdg[evt_idx] or -211 not in reco_pdg[evt_idx]:
            continue

        # Find protons and pions
        proton_mask = (reco_pdg[evt_idx] == 2212)
        pion_mask = (reco_pdg[evt_idx] == -211)

        proton_indices = ak.where(proton_mask)[0]
        pion_indices = ak.where(pion_mask)[0]

        # Loop through MC Lambdas that match this event
        for lambda_idx, mc_lambda_mass in enumerate(lambda_mc_info["mass"]):
            # Get MC decay product momenta
            p_px_mc = proton_px_mc[lambda_idx]
            p_py_mc = proton_py_mc[lambda_idx]
            p_pz_mc = proton_pz_mc[lambda_idx]
            pi_px_mc = pion_px_mc[lambda_idx]
            pi_py_mc = pion_py_mc[lambda_idx]
            pi_pz_mc = pion_pz_mc[lambda_idx]

            best_match_score = float('inf')
            best_match_mass = None
            best_p_pt_mc = None
            best_p_pt_reco = None
            best_pi_pt_mc = None
            best_pi_pt_reco = None

            # Calculate MC pT values
            p_pt_mc = np.sqrt(p_px_mc**2 + p_py_mc**2)
            pi_pt_mc = np.sqrt(pi_px_mc**2 + pi_py_mc**2)

            # Try all combinations of reconstructed protons and pions
            for p_idx in proton_indices:
                p_px_reco = reco_px[evt_idx][p_idx]
                p_py_reco = reco_py[evt_idx][p_idx]
                p_pz_reco = reco_pz[evt_idx][p_idx]
                p_m_reco = reco_mass[evt_idx][p_idx]
                p_pt_reco = np.sqrt(p_px_reco**2 + p_py_reco**2)

                for pi_idx in pion_indices:
                    pi_px_reco = reco_px[evt_idx][pi_idx]
                    pi_py_reco = reco_py[evt_idx][pi_idx]
                    pi_pz_reco = reco_pz[evt_idx][pi_idx]
                    pi_m_reco = reco_mass[evt_idx][pi_idx]
                    pi_pt_reco = np.sqrt(pi_px_reco**2 + pi_py_reco**2)

                    # Calculate match score (simpler version - just pT difference)
                    match_score = abs(p_pt_reco - p_pt_mc) + abs(pi_pt_reco - pi_pt_mc)

                    if match_score < best_match_score:
                        best_match_score = match_score

                        # Calculate corrected energy values
                        p_e_reco = np.sqrt(p_px_reco**2 + p_py_reco**2 + p_pz_reco**2 + p_m_reco**2)
                        pi_e_reco = np.sqrt(pi_px_reco**2 + pi_py_reco**2 + pi_pz_reco**2 + pi_m_reco**2)

                        # Calculate invariant mass
                        p_vec = vector.obj(px=p_px_reco, py=p_py_reco, pz=p_pz_reco, E=p_e_reco)
                        pi_vec = vector.obj(px=pi_px_reco, py=pi_py_reco, pz=pi_pz_reco, E=pi_e_reco)
                        reco_vec = p_vec + pi_vec

                        best_match_mass = reco_vec.mass
                        best_p_pt_mc = p_pt_mc
                        best_p_pt_reco = p_pt_reco
                        best_pi_pt_mc = pi_pt_mc
                        best_pi_pt_reco = pi_pt_reco

            # Save the best match if we found one
            if best_match_score < 1.0:  # Arbitrary threshold
                matches["lambda_idx"].append(lambda_idx)
                matches["mc_mass"].append(mc_lambda_mass)
                matches["reco_mass"].append(best_match_mass)
                matches["proton_pt_mc"].append(best_p_pt_mc)
                matches["proton_pt_reco"].append(best_p_pt_reco)
                matches["pion_pt_mc"].append(best_pi_pt_mc)
                matches["pion_pt_reco"].append(best_pi_pt_reco)

    return matches


def main():
    args = parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(args.output) if os.path.dirname(args.output) else '.', exist_ok=True)

    # Open the file and get the tree
    filename = args.input_file
    treename = args.tree

    print(f"Processing file: {filename}")
    print(f"Tree: {treename}")

    # Read the file
    file = uproot.open(filename)
    tree = file[treename]

    # Check for MCParticles collection
    if "MCParticles" not in tree.keys():
        print("Error: MCParticles collection not found in the tree")
        return

    # Get needed branches for MC truth
    mc_branches = [
        "MCParticles.PDG",
        "MCParticles.mass",
        "MCParticles.momentum.x",
        "MCParticles.momentum.y",
        "MCParticles.momentum.z",
        "MCParticles.daughters_begin",
        "MCParticles.daughters_end"
    ]

    # Get needed branches for reconstructed particles
    reco_main_branch = "ReconstructedChargedParticles"
    reco_branches = [
        f"{reco_main_branch}.PDG",
        f"{reco_main_branch}.momentum.x",
        f"{reco_main_branch}.momentum.y",
        f"{reco_main_branch}.momentum.z",
        f"{reco_main_branch}.energy",
        f"{reco_main_branch}.mass"
    ]

    # Read MC data
    print("Reading MC truth data...")
    mc_data = tree.arrays(mc_branches)
    print(ak.type(mc_data[0]))

    # Read reconstructed data
    print("Reading reconstructed particle data...")
    reco_data = tree.arrays(reco_branches)
    print(ak.type(reco_data[0]))

    # Find Lambda decays in MC truth
    print("Finding Lambda decays in MC truth...")
    lambda_info = find_lambda_decays(mc_data)

    num_lambdas = len(lambda_info["mass"])
    print(f"Found {num_lambdas} Lambda->p+π⁻ decays in MC truth")

    if num_lambdas == 0:
        print("No Lambda decays found in MC truth. Exiting.")
        return

    # Print some Lambda properties
    print("\nMC Lambda properties (first 5):")
    for i in range(min(5, num_lambdas)):
        print(f"  Lambda {i}:")
        print(f"    True mass: {lambda_info['mass'][i]:.6f} GeV")
        print(f"    Reconstructed mass from daughters: {lambda_info['reconstructed_mass'][i]:.6f} GeV")
        print(f"    Momentum: px={lambda_info['px'][i]:.3f}, py={lambda_info['py'][i]:.3f}, pz={lambda_info['pz'][i]:.3f} GeV")
        print(f"    Proton momentum: px={lambda_info['proton_px'][i]:.3f}, py={lambda_info['proton_py'][i]:.3f}, pz={lambda_info['proton_pz'][i]:.3f} GeV")
        print(f"    Pion momentum: px={lambda_info['pion_px'][i]:.3f}, py={lambda_info['pion_py'][i]:.3f}, pz={lambda_info['pion_pz'][i]:.3f} GeV")

        # Calculate pT
        lambda_pt = np.sqrt(lambda_info['px'][i]**2 + lambda_info['py'][i]**2)
        proton_pt = np.sqrt(lambda_info['proton_px'][i]**2 + lambda_info['proton_py'][i]**2)
        pion_pt = np.sqrt(lambda_info['pion_px'][i]**2 + lambda_info['pion_py'][i]**2)

        print(f"    Lambda pT: {lambda_pt:.3f} GeV")
        print(f"    Proton pT: {proton_pt:.3f} GeV")
        print(f"    Pion pT: {pion_pt:.3f} GeV")

    # Plot histograms of MC Lambda properties

    # 1. Lambda mass
    plt.figure(figsize=(10, 6))
    plt.hist(lambda_info["mass"], bins=50, range=(1.0, 1.2))
    plt.axvline(1.11568, color='r', linestyle='--', label="PDG Lambda mass (1.11568 GeV)")
    plt.xlabel("Lambda mass [GeV]")
    plt.ylabel("Count")
    plt.title("MC Lambda Mass")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.savefig(f"{args.output}_mc_lambda_mass.png")

    # 2. Lambda pT
    lambda_pt = np.sqrt(np.array(lambda_info["px"])**2 + np.array(lambda_info["py"])**2)
    plt.figure(figsize=(10, 6))
    plt.hist(lambda_pt, bins=50)
    plt.xlabel("Lambda pT [GeV]")
    plt.ylabel("Count")
    plt.title("MC Lambda pT")
    plt.grid(True, alpha=0.3)
    plt.savefig(f"{args.output}_mc_lambda_pt.png")

    # 3. Proton pT
    proton_pt = np.sqrt(np.array(lambda_info["proton_px"])**2 + np.array(lambda_info["proton_py"])**2)
    plt.figure(figsize=(10, 6))
    plt.hist(proton_pt, bins=50)
    plt.xlabel("Proton pT [GeV]")
    plt.ylabel("Count")
    plt.title("Proton pT from Lambda Decay")
    plt.grid(True, alpha=0.3)
    plt.savefig(f"{args.output}_mc_proton_pt.png")

    # 4. Pion pT
    pion_pt = np.sqrt(np.array(lambda_info["pion_px"])**2 + np.array(lambda_info["pion_py"])**2)
    plt.figure(figsize=(10, 6))
    plt.hist(pion_pt, bins=50)
    plt.xlabel("Pion pT [GeV]")
    plt.ylabel("Count")
    plt.title("Pion pT from Lambda Decay")
    plt.grid(True, alpha=0.3)
    plt.savefig(f"{args.output}_mc_pion_pt.png")

    # 5. Comparison of MC mass vs reconstructed mass
    plt.figure(figsize=(10, 6))
    plt.scatter(lambda_info["mass"], lambda_info["reconstructed_mass"], alpha=0.5)
    plt.plot([1.0, 1.2], [1.0, 1.2], 'r--')  # Diagonal line
    plt.xlabel("MC Lambda mass [GeV]")
    plt.ylabel("Reconstructed mass from MC daughters [GeV]")
    plt.title("MC Mass vs Reconstructed Mass")
    plt.grid(True, alpha=0.3)
    plt.savefig(f"{args.output}_mc_vs_reco_mass.png")

    # Try to find matching reconstructed particles
    print("\nLooking for matches between MC truth and reconstruction...")
    matches = find_reconstructed_matches(reco_data, reco_main_branch, mc_data, lambda_info)

    num_matches = len(matches["lambda_idx"])
    print(f"Found {num_matches} potential matches between MC and reconstructed particles")

    if num_matches > 0:
        # Plot matching results
        plt.figure(figsize=(10, 6))
        plt.scatter(matches["mc_mass"], matches["reco_mass"], alpha=0.5)
        plt.plot([1.0, 1.2], [1.0, 1.2], 'r--')  # Diagonal line
        plt.xlabel("MC Lambda mass [GeV]")
        plt.ylabel("Reconstructed Lambda mass [GeV]")
        plt.title("MC vs Reconstructed Lambda Mass")
        plt.grid(True, alpha=0.3)
        plt.savefig(f"{args.output}_match_mc_vs_reco_mass.png")

        # Plot pT comparisons
        plt.figure(figsize=(10, 6))
        plt.scatter(matches["proton_pt_mc"], matches["proton_pt_reco"], alpha=0.5, label="Protons")
        plt.scatter(matches["pion_pt_mc"], matches["pion_pt_reco"], alpha=0.5, label="Pions")
        plt.plot([0, 10], [0, 10], 'r--')  # Diagonal line
        plt.xlabel("MC particle pT [GeV]")
        plt.ylabel("Reconstructed particle pT [GeV]")
        plt.title("MC vs Reconstructed Particle pT")
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.savefig(f"{args.output}_match_pt_comparison.png")

    print(f"\nPlots saved with prefix: {args.output}")

if __name__ == "__main__":
    main()