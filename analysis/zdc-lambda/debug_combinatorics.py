#!/usr/bin/env python3

import uproot
import awkward as ak
import numpy as np
import vector
import matplotlib.pyplot as plt
import argparse
import os


def fix_four_vectors(data, branch_prefix):
    """
    Fix the energy component of 4-vectors using stored mass values.
    This corrects for issues where E = |p| incorrectly.

    Parameters:
    -----------
    data : awkward array
        The data from uproot containing particle information
    branch_prefix : str
        Prefix for the branch names (e.g., "ReconstructedParticles.")

    Returns:
    --------
    awkward array
        Data with corrected energy values
    """
    # Get momentum components
    px = data[f"{branch_prefix}momentum.x"]
    py = data[f"{branch_prefix}momentum.y"]
    pz = data[f"{branch_prefix}momentum.z"]

    # Calculate momentum magnitude
    p_squared = px**2 + py**2 + pz**2
    p_mag = np.sqrt(p_squared)

    # Get mass from branch (this appears to be correct)
    mass = data[f"{branch_prefix}mass"]

    # Calculate correct energy: E^2 = p^2 + m^2
    energy_squared = p_squared + mass**2
    corrected_energy = np.sqrt(energy_squared)

    # Create a copy of the data
    fixed_data = ak.copy(data)

    # Replace the energy values
    fixed_data[f"{branch_prefix}energy"] = corrected_energy

    return fixed_data


def check_energy_momentum_consistency(data, branch_prefix, output_dir):
    """
    Check if the energy and momentum values are consistent with the mass.

    Parameters:
    -----------
    data : awkward array
        The data containing particle information
    branch_prefix : str
        Prefix for the branch names
    """
    # Get momentum components and energy
    px = data[f"{branch_prefix}momentum.x"]
    py = data[f"{branch_prefix}momentum.y"]
    pz = data[f"{branch_prefix}momentum.z"]
    energy = data[f"{branch_prefix}energy"]

    # Calculate momentum magnitude squared
    p_squared = px**2 + py**2 + pz**2

    # Calculate E^2 - p^2 (should equal m^2)
    m_squared = energy**2 - p_squared

    # Flatten arrays for histogram
    flat_m_squared = ak.flatten(m_squared)

    # Plot the distribution of m^2
    plt.figure(figsize=(10, 6))
    plt.hist(ak.to_numpy(flat_m_squared), bins=100, range=(-0.1, 2.0))
    plt.xlabel("E^2 - p^2 (Should equal m^2) [GeV^2]")
    plt.ylabel("Count")
    plt.title("Consistency Check: E^2 - p^2 Distribution")
    plt.axvline(0.938**2, color='r', linestyle='--',
                label=f"Proton m^2 ({0.938**2:.3f} GeV^2)")
    plt.axvline(0.140**2, color='g', linestyle='--',
                label=f"Pion m^2 ({0.140**2:.3f} GeV^2)")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.savefig(f"{output_dir}_energy_momentum_check.png")

    return

def main():
    parser = argparse.ArgumentParser(description="Diagnostic tool for Lambda reconstruction")
    parser.add_argument("-i", "--input-file", required=True, help="Input ROOT file")
    parser.add_argument("-o", "--output", default="lambda_diagnostics", help="Output directory")
    parser.add_argument("--collection", default="ReconstructedParticles", help="Particle collection name")
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(args.output) if os.path.dirname(args.output) else '.', exist_ok=True)

    # Open the file and get the tree
    filename = args.input_file
    collection = args.collection

    print(f"Processing file: {filename}")
    print(f"Collection: {collection}")

    # Build the branch names
    branch_prefix = f"{collection}."
    needed_branches = [
        f"{branch_prefix}PDG",
        f"{branch_prefix}charge",
        f"{branch_prefix}momentum.x",
        f"{branch_prefix}momentum.y",
        f"{branch_prefix}momentum.z",
        f"{branch_prefix}energy",
        f"{branch_prefix}mass",  # This might exist and help diagnose issues
    ]

    # Read the file
    tree = uproot.open(filename)["events"]

    # Read the data
    data = tree.arrays(needed_branches)

    # Check consistency before fixing
    print("\nChecking energy-momentum consistency before fixing...")
    check_energy_momentum_consistency(data, branch_prefix)

    # Fix the energy values using the correct mass
    print("Fixing energy values using the mass information...")
    fixed_data = fix_four_vectors(data, branch_prefix)

    # Check consistency after fixing
    print("Checking energy-momentum consistency after fixing...")
    check_energy_momentum_consistency(fixed_data, branch_prefix)

    # Continue with analysis using fixed_data
    data = fixed_data

    # Check PDG code distribution
    pdg = data[f"{branch_prefix}PDG"]
    all_pdgs = ak.flatten(pdg)
    unique_pdgs, counts = np.unique(ak.to_numpy(all_pdgs), return_counts=True)

    print("\nPDG code distribution:")
    for pdg_code, count in zip(unique_pdgs, counts):
        print(f"  PDG {pdg_code}: {count} particles")

    # Check for protons (2212) and pions (-211)
    proton_mask = (pdg == 2212)
    pion_mask = (pdg == -211)

    n_protons = ak.sum(proton_mask, axis=1)
    n_pions = ak.sum(pion_mask, axis=1)

    print(f"\nEvents with protons: {ak.sum(n_protons > 0)} / {len(n_protons)}")
    print(f"Events with pions: {ak.sum(n_pions > 0)} / {len(n_pions)}")
    print(f"Events with both: {ak.sum((n_protons > 0) & (n_pions > 0))} / {len(n_pions)}")

    # Get momentum components
    px = data[f"{branch_prefix}momentum.x"]
    py = data[f"{branch_prefix}momentum.y"]
    pz = data[f"{branch_prefix}momentum.z"]
    energy = data[f"{branch_prefix}energy"]

    # Check if mass is available
    mass_branch = f"{branch_prefix}mass"
    if mass_branch in data.fields:
        mass = data[mass_branch]
        has_mass = True
        print("\nMass branch is available")

        # Check some proton masses
        proton_mass = mass[proton_mask]
        if len(proton_mass) > 0:
            flat_proton_mass = ak.flatten(proton_mass)
            if len(flat_proton_mass) > 0:
                print(f"First 5 proton masses: {ak.to_numpy(flat_proton_mass[:5])}")
                print(f"Mean proton mass: {ak.mean(flat_proton_mass)}")

        # Check some pion masses
        pion_mass = mass[pion_mask]
        if len(pion_mass) > 0:
            flat_pion_mass = ak.flatten(pion_mass)
            if len(flat_pion_mass) > 0:
                print(f"First 5 pion masses: {ak.to_numpy(flat_pion_mass[:5])}")
                print(f"Mean pion mass: {ak.mean(flat_pion_mass)}")
    else:
        has_mass = False
        print("\nMass branch is not available")

    # Print the first proton and pion found
    for evt_idx in range(len(pdg)):
        if ak.sum(proton_mask[evt_idx]) > 0 and ak.sum(pion_mask[evt_idx]) > 0:
            # Get proton and pion indices
            p_idx = ak.argmax(proton_mask[evt_idx], axis=0, keepdims=False)
            pi_idx = ak.argmax(pion_mask[evt_idx], axis=0, keepdims=False)

            # Get proton properties
            p_px = px[evt_idx][p_idx]
            p_py = py[evt_idx][p_idx]
            p_pz = pz[evt_idx][p_idx]
            p_E = energy[evt_idx][p_idx]

            # Get pion properties
            pi_px = px[evt_idx][pi_idx]
            pi_py = py[evt_idx][pi_idx]
            pi_pz = pz[evt_idx][pi_idx]
            pi_E = energy[evt_idx][pi_idx]

            print(f"\nExample from event {evt_idx}:")
            print(f"  Proton (PDG={pdg[evt_idx][p_idx]}):")
            print(f"    px = {p_px}, py = {p_py}, pz = {p_pz}, E = {p_E}")
            if has_mass:
                print(f"    mass from branch = {mass[evt_idx][p_idx]}")

            print(f"  Pion (PDG={pdg[evt_idx][pi_idx]}):")
            print(f"    px = {pi_px}, py = {pi_py}, pz = {pi_pz}, E = {pi_E}")
            if has_mass:
                print(f"    mass from branch = {mass[evt_idx][pi_idx]}")

            # Calculate pT, E^2 - p^2, etc.
            p_pt = np.sqrt(p_px**2 + p_py**2)
            p_p = np.sqrt(p_px**2 + p_py**2 + p_pz**2)
            p_m2 = p_E**2 - p_p**2

            pi_pt = np.sqrt(pi_px**2 + pi_py**2)
            pi_p = np.sqrt(pi_px**2 + pi_py**2 + pi_pz**2)
            pi_m2 = pi_E**2 - pi_p**2

            print("  Derived quantities:")
            print(f"    Proton: pT = {p_pt}, p = {p_p}, E^2-p^2 = {p_m2}, sqrt(E^2-p^2) = {np.sqrt(abs(p_m2)) if p_m2 > 0 else 'imaginary'}")
            print(f"    Pion: pT = {pi_pt}, p = {pi_p}, E^2-p^2 = {pi_m2}, sqrt(E^2-p^2) = {np.sqrt(abs(pi_m2)) if pi_m2 > 0 else 'imaginary'}")

            # Try different ways to create 4-vectors

            # Method 1: Direct from px,py,pz,E
            try:
                p_vec1 = vector.obj(px=p_px, py=p_py, pz=p_pz, E=p_E)
                pi_vec1 = vector.obj(px=pi_px, py=pi_py, pz=pi_pz, E=pi_E)
                lambda_vec1 = p_vec1 + pi_vec1
                print("\n  Vector calculations method 1 (px,py,pz,E):")
                print(f"    Proton mass = {p_vec1.mass}")
                print(f"    Pion mass = {pi_vec1.mass}")
                print(f"    Lambda mass = {lambda_vec1.mass}")
                print(f"    Lambda pT = {lambda_vec1.pt}")
            except Exception as e:
                print(f"\n  Vector calculations method 1 failed: {e}")

            # Method 2: Try with PDG masses
            try:
                proton_pdg_mass = 0.93827  # GeV
                pion_pdg_mass = 0.13957  # GeV

                # Recalculate energy using PDG mass
                p_E2 = np.sqrt(p_p**2 + proton_pdg_mass**2)
                pi_E2 = np.sqrt(pi_p**2 + pion_pdg_mass**2)

                p_vec2 = vector.obj(px=p_px, py=p_py, pz=p_pz, E=p_E2)
                pi_vec2 = vector.obj(px=pi_px, py=pi_py, pz=pi_pz, E=pi_E2)
                lambda_vec2 = p_vec2 + pi_vec2
                print("\n  Vector calculations method 2 (with PDG masses):")
                print(f"    Proton recalculated E = {p_E2} (original = {p_E})")
                print(f"    Pion recalculated E = {pi_E2} (original = {pi_E})")
                print(f"    Proton mass = {p_vec2.mass}")
                print(f"    Pion mass = {pi_vec2.mass}")
                print(f"    Lambda mass = {lambda_vec2.mass}")
                print(f"    Lambda pT = {lambda_vec2.pt}")
            except Exception as e:
                print(f"\n  Vector calculations method 2 failed: {e}")

            # Method 3: Try with pt, eta, phi, mass
            try:
                p_eta = 0.5 * np.log((p_p + p_pz) / (p_p - p_pz)) if p_p != p_pz else 0
                p_phi = np.arctan2(p_py, p_px)

                pi_eta = 0.5 * np.log((pi_p + pi_pz) / (pi_p - pi_pz)) if pi_p != pi_pz else 0
                pi_phi = np.arctan2(pi_py, pi_px)

                p_vec3 = vector.obj(pt=p_pt, eta=p_eta, phi=p_phi, mass=proton_pdg_mass)
                pi_vec3 = vector.obj(pt=pi_pt, eta=pi_eta, phi=pi_phi, mass=pion_pdg_mass)
                lambda_vec3 = p_vec3 + pi_vec3
                print("\n  Vector calculations method 3 (pt,eta,phi,mass):")
                print(f"    Proton: pt={p_pt}, eta={p_eta}, phi={p_phi}")
                print(f"    Pion: pt={pi_pt}, eta={pi_eta}, phi={pi_phi}")
                print(f"    Proton mass = {p_vec3.mass}")
                print(f"    Pion mass = {pi_vec3.mass}")
                print(f"    Lambda mass = {lambda_vec3.mass}")
                print(f"    Lambda pT = {lambda_vec3.pt}")
            except Exception as e:
                print(f"\n  Vector calculations method 3 failed: {e}")

            # Just try the pure calculation directly
            try:
                # Convert to array for easier calculation
                p_4vec = np.array([p_E, p_px, p_py, p_pz])
                pi_4vec = np.array([pi_E, pi_px, pi_py, pi_pz])
                lambda_4vec = p_4vec + pi_4vec

                # Calculate mass from E^2 - p^2
                lambda_m2 = lambda_4vec[0]**2 - (lambda_4vec[1]**2 + lambda_4vec[2]**2 + lambda_4vec[3]**2)
                lambda_m = np.sqrt(lambda_m2) if lambda_m2 > 0 else np.sqrt(abs(lambda_m2))

                print("\n  Pure numpy 4-vector calculation:")
                print(f"    Lambda 4-vector = [{lambda_4vec[0]}, {lambda_4vec[1]}, {lambda_4vec[2]}, {lambda_4vec[3]}]")
                print(f"    Lambda m^2 = {lambda_m2}")
                print(f"    Lambda mass = {lambda_m}" + (" (imaginary!)" if lambda_m2 < 0 else ""))
            except Exception as e:
                print(f"\n  Pure numpy calculation failed: {e}")

            # Exit after finding the first example
            break

    # Create mass histograms
    lambda_masses = []
    lambda_masses_method2 = []
    lambda_masses_method3 = []

    # Try to accumulate masses using different methods
    for evt_idx in range(len(pdg)):
        # Get event masks
        evt_proton_mask = proton_mask[evt_idx]
        evt_pion_mask = pion_mask[evt_idx]

        if ak.sum(evt_proton_mask) > 0 and ak.sum(evt_pion_mask) > 0:
            # Get proton and pion indices
            proton_indices = ak.where(evt_proton_mask)[0]
            pion_indices = ak.where(evt_pion_mask)[0]

            # Process each combination
            for p_idx in proton_indices:
                p_px = px[evt_idx][p_idx]
                p_py = py[evt_idx][p_idx]
                p_pz = pz[evt_idx][p_idx]
                p_E = energy[evt_idx][p_idx]
                p_p = np.sqrt(p_px**2 + p_py**2 + p_pz**2)

                for pi_idx in pion_indices:
                    pi_px = px[evt_idx][pi_idx]
                    pi_py = py[evt_idx][pi_idx]
                    pi_pz = pz[evt_idx][pi_idx]
                    pi_E = energy[evt_idx][pi_idx]
                    pi_p = np.sqrt(pi_px**2 + pi_py**2 + pi_pz**2)

                    # Method 1: Direct from px,py,pz,E (with fixed energies)
                    try:
                        p_vec = vector.obj(px=p_px, py=p_py, pz=p_pz, E=p_E)
                        pi_vec = vector.obj(px=pi_px, py=pi_py, pz=pi_pz, E=pi_E)
                        lambda_vec = p_vec + pi_vec
                        lambda_masses.append(lambda_vec.mass)
                    except Exception as e:
                        continue

                    # Method 2: Verify with PDG masses
                    try:
                        proton_pdg_mass = 0.93827  # GeV
                        pion_pdg_mass = 0.13957  # GeV

                        # Recalculate energy using PDG mass
                        p_E2 = np.sqrt(p_p**2 + proton_pdg_mass**2)
                        pi_E2 = np.sqrt(pi_p**2 + pion_pdg_mass**2)

                        p_vec2 = vector.obj(px=p_px, py=p_py, pz=p_pz, E=p_E2)
                        pi_vec2 = vector.obj(px=pi_px, py=pi_py, pz=pi_pz, E=pi_E2)
                        lambda_vec2 = p_vec2 + pi_vec2
                        lambda_masses_method2.append(lambda_vec2.mass)
                    except Exception as e:
                        continue

                    # Method 3: Using pt, eta, phi, mass
                    try:
                        p_pt = np.sqrt(p_px**2 + p_py**2)
                        p_eta = 0.5 * np.log((p_p + p_pz) / (p_p - p_pz)) if abs(p_p - p_pz) > 1e-10 else 0
                        p_phi = np.arctan2(p_py, p_px)

                        pi_pt = np.sqrt(pi_px**2 + pi_py**2)
                        pi_eta = 0.5 * np.log((pi_p + pi_pz) / (pi_p - pi_pz)) if abs(pi_p - pi_pz) > 1e-10 else 0
                        pi_phi = np.arctan2(pi_py, pi_px)

                        p_vec3 = vector.obj(pt=p_pt, eta=p_eta, phi=p_phi, mass=proton_pdg_mass)
                        pi_vec3 = vector.obj(pt=pi_pt, eta=pi_eta, phi=pi_phi, mass=pion_pdg_mass)
                        lambda_vec3 = p_vec3 + pi_vec3
                        lambda_masses_method3.append(lambda_vec3.mass)
                    except Exception as e:
                        continue

    # Print summary statistics
    print("\nSummary of Lambda reconstruction:")
    print(f"  Method 1 (Fixed energies): {len(lambda_masses)} candidates")
    print(f"  Method 2 (PDG masses): {len(lambda_masses_method2)} candidates")
    print(f"  Method 3 (pt,eta,phi,mass): {len(lambda_masses_method3)} candidates")

    if len(lambda_masses) > 0:
        print(f"  Method 1 mean mass: {np.mean(lambda_masses):.4f} GeV")
    if len(lambda_masses_method2) > 0:
        print(f"  Method 2 mean mass: {np.mean(lambda_masses_method2):.4f} GeV")
    if len(lambda_masses_method3) > 0:
        print(f"  Method 3 mean mass: {np.mean(lambda_masses_method3):.4f} GeV")

    # Plot mass histogram - using a very wide range for diagnostic purposes
    plt.figure(figsize=(12, 6))

    if len(lambda_masses) > 0:
        plt.hist(lambda_masses, bins=100, range=(0, 10),
                 label=f"Method 1: Fixed E (n={len(lambda_masses)})", alpha=0.7)

    if len(lambda_masses_method2) > 0:
        plt.hist(lambda_masses_method2, bins=100, range=(0, 10), alpha=0.7,
                 label=f"Method 2: PDG masses (n={len(lambda_masses_method2)})")

    if len(lambda_masses_method3) > 0:
        plt.hist(lambda_masses_method3, bins=100, range=(0, 10), alpha=0.7,
                 label=f"Method 3: pt,eta,phi (n={len(lambda_masses_method3)})")

    plt.axvline(1.11568, color='r', linestyle='--', label="PDG Lambda mass (1.11568 GeV)")
    plt.legend()
    plt.xlabel("Lambda mass [GeV]")
    plt.ylabel("Counts")
    plt.title("Lambda mass distribution (wide range)")
    plt.grid(True, alpha=0.3)
    plt.savefig(f"{args.output}_wide_mass.png")

    # Plot with a narrower range focusing on the PDG mass
    plt.figure(figsize=(12, 6))

    if len(lambda_masses) > 0:
        plt.hist(lambda_masses, bins=100, range=(0.5, 2.0),
                 label=f"Method 1: Fixed E (n={len(lambda_masses)})", alpha=0.7)

    if len(lambda_masses_method2) > 0:
        plt.hist(lambda_masses_method2, bins=100, range=(0.5, 2.0), alpha=0.7,
                 label=f"Method 2: PDG masses (n={len(lambda_masses_method2)})")

    if len(lambda_masses_method3) > 0:
        plt.hist(lambda_masses_method3, bins=100, range=(0.5, 2.0), alpha=0.7,
                 label=f"Method 3: pt,eta,phi (n={len(lambda_masses_method3)})")

    plt.axvline(1.11568, color='r', linestyle='--', label="PDG Lambda mass (1.11568 GeV)")
    plt.legend()
    plt.xlabel("Lambda mass [GeV]")
    plt.ylabel("Counts")
    plt.title("Lambda mass distribution (focused range)")
    plt.grid(True, alpha=0.3)
    plt.savefig(f"{args.output}_narrow_mass.png")

    print(f"\nPlots saved to:")
    print(f"  {args.output}_energy_momentum_check.png")
    print(f"  {args.output}_wide_mass.png")
    print(f"  {args.output}_narrow_mass.png")

if __name__ == "__main__":
    main()