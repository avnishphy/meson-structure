#!/usr/bin/env python3

import argparse
import os
import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(
        description="Analyze reconstructed particles and plot angular distributions."
    )
    parser.add_argument("-i", "--input-file", required=True,
                        help="EDM4hep ROOT file containing ReconstructedParticles.")
    parser.add_argument("-t", "--tree-name", default="events",
                        help="Name of the TTree (default 'events').")
    parser.add_argument("-o", "--outdir", default="particle_plots",
                        help="Output directory for plots.")
    parser.add_argument("-n", "--num-events", type=int, default=None,
                        help="Number of events to process (default: all events).")
    return parser.parse_args()

def get_particle_name(pdg_code):
    """Map PDG code to particle name."""
    pdg_names = {
        # Leptons
        11: "e-",
        -11: "e+",
        13: "μ-",
        -13: "μ+",
        15: "τ-",
        -15: "τ+",
        # Neutrinos
        12: "νe",
        -12: "ν̄e",
        14: "νμ",
        -14: "ν̄μ",
        16: "ντ",
        -16: "ν̄τ",
        # Light mesons
        111: "π⁰",
        211: "π+",
        -211: "π-",
        221: "η",
        321: "K+",
        -321: "K-",
        310: "K⁰S",
        130: "K⁰L",
        # Light baryons
        2212: "p",
        -2212: "p̄",
        2112: "n",
        -2112: "n̄",
        # Strange baryons
        3122: "Λ",
        -3122: "Λ̄",
        3222: "Σ+",
        3212: "Σ⁰",
        3112: "Σ-",
        # Photon and other bosons
        22: "γ",
        23: "Z",
        24: "W+",
        -24: "W-",
        25: "h",
        # Zero code
        0: "unknown"
    }
    return pdg_names.get(pdg_code, f"PDG {pdg_code}")

def calculate_theta_phi(px, py, pz):
    """
    Calculate theta and phi angles from momentum components.
    Returns angles in degrees.

    Parameters:
    -----------
    px, py, pz : arrays
        Momentum components

    Returns:
    --------
    tuple
        (theta, phi) in degrees
    """
    p = np.sqrt(px**2 + py**2 + pz**2)
    # Handle zero momentum case
    mask = (p > 0)

    theta = np.zeros_like(p)
    theta[mask] = np.arccos(pz[mask] / p[mask]) * 180.0 / np.pi  # Convert to degrees

    phi = np.arctan2(py, px) * 180.0 / np.pi   # Convert to degrees

    return theta, phi

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print(f"Processing file: {args.input_file}")

    # Open ROOT file
    file = uproot.open(args.input_file)
    tree = file[args.tree_name]

    # Define branches to read
    branches = [
        "ReconstructedParticles.PDG",
        "ReconstructedParticles.momentum.x",
        "ReconstructedParticles.momentum.y",
        "ReconstructedParticles.momentum.z"
    ]

    # Read data
    print("Reading reconstructed particle data...")
    data = tree.arrays(branches, library="ak", entry_stop=args.num_events)

    # Get PDG codes and momentum
    pdg_codes = data["ReconstructedParticles.PDG"]
    px = data["ReconstructedParticles.momentum.x"]
    py = data["ReconstructedParticles.momentum.y"]
    pz = data["ReconstructedParticles.momentum.z"]

    # Flatten PDG codes to count occurrences
    all_pdgs = ak.flatten(pdg_codes)

    # Count occurrences of each particle type
    pdg_counts = {}
    for pdg in np.unique(ak.to_numpy(all_pdgs)):
        count = np.sum(all_pdgs == pdg)
        pdg_counts[pdg] = count

    # Sort particles by occurrence (descending)
    sorted_pdgs = sorted(pdg_counts.items(), key=lambda x: x[1], reverse=True)

    # Create a list of (PDG code, name, count)
    particles_data = [(pdg, get_particle_name(pdg), count) for pdg, count in sorted_pdgs]

    # Create a DataFrame for easier plotting
    df = pd.DataFrame(particles_data, columns=['PDG', 'Particle', 'Count'])
    print("\nParticle counts:")
    print(df)

    # Plot particle counts (horizontal bar chart)
    plt.figure(figsize=(12, max(8, len(df) * 0.4)))  # Adjust height based on number of particles
    bars = plt.barh(df['Particle'], df['Count'])

    # Add count numbers at the end of each bar
    for bar in bars:
        width = bar.get_width()
        plt.text(width*1.01, bar.get_y() + bar.get_height()/2,
                 f'{int(width)}', va='center')

    plt.xlabel('Count')
    plt.title('Reconstructed Particle Types')
    plt.grid(axis='x', alpha=0.3)
    plt.tight_layout()

    # Save the plot
    particle_counts_path = os.path.join(args.outdir, "particle_counts.png")
    plt.savefig(particle_counts_path)
    print(f"Saved particle counts plot to: {particle_counts_path}")
    plt.close()

    # Process angular distributions for protons and pions
    particle_types = {
        'proton': 2212,
        'antiproton': -2212,
        'piplus': 211,
        'piminus': -211
    }

    # Dictionary to store angular distributions
    angular_data = {}

    print("\nCollecting angular distributions...")
    for particle_name, pdg_code in particle_types.items():
        print(f"Processing {particle_name} (PDG={pdg_code})...")

        # Initialize arrays for this particle type
        theta_values = []
        phi_values = []

        # Loop through all events
        num_events = len(pdg_codes)
        event_counter = 0  # For debug purposes

        for evt_idx in range(num_events):
            # Find all particles of this type in the event
            particle_mask = (pdg_codes[evt_idx] == pdg_code)

            # Skip if no particles of this type
            if ak.sum(particle_mask) == 0:
                continue

            event_counter += 1

            # Get indices of the particles
            indices = ak.where(particle_mask)[0]

            # Get momentum components
            event_px = px[evt_idx][indices]
            event_py = py[evt_idx][indices]
            event_pz = pz[evt_idx][indices]

            # Calculate theta and phi
            theta, phi = calculate_theta_phi(event_px, event_py, event_pz)

            # Add to our lists
            theta_values.extend(ak.to_list(theta))
            phi_values.extend(ak.to_list(phi))

        # Store the results
        angular_data[particle_name] = {
            'theta': theta_values,
            'phi': phi_values
        }

        print(f"  Found in {event_counter} events")
        print(f"  Collected {len(theta_values)} particles")

    # Create plots for angular distributions
    print("\nCreating angular distribution plots...")

    # Set up the figure with 2 rows (theta, phi) and 4 columns (particle types)
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))

    # Plot theta distributions in the first row
    particle_labels = {'proton': 'Proton', 'antiproton': 'Antiproton',
                       'piplus': 'π+', 'piminus': 'π-'}

    row_labels = ['θ distribution (degrees)', 'φ distribution (degrees)']

    for col, (particle_name, data) in enumerate(angular_data.items()):
        # Plot theta (top row)
        if len(data['theta']) > 0:
            axes[0, col].hist(data['theta'], bins=36, range=(0, 180))
            axes[0, col].set_title(f"{particle_labels[particle_name]} θ")
            axes[0, col].set_xlabel('θ (degrees)')
            axes[0, col].text(0.05, 0.95, f"n={len(data['theta'])}",
                              transform=axes[0, col].transAxes,
                              verticalalignment='top',
                              bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
        else:
            axes[0, col].text(0.5, 0.5, "No data", ha='center', va='center',
                              transform=axes[0, col].transAxes)

        # Plot phi (bottom row)
        if len(data['phi']) > 0:
            axes[1, col].hist(data['phi'], bins=36, range=(-180, 180))
            axes[1, col].set_title(f"{particle_labels[particle_name]} φ")
            axes[1, col].set_xlabel('φ (degrees)')
            axes[1, col].text(0.05, 0.95, f"n={len(data['phi'])}",
                              transform=axes[1, col].transAxes,
                              verticalalignment='top',
                              bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
        else:
            axes[1, col].text(0.5, 0.5, "No data", ha='center', va='center',
                              transform=axes[1, col].transAxes)

    # Add row labels
    for ax, row in zip(axes[:,0], row_labels):
        ax.set_ylabel(row, size='large')

    plt.tight_layout()
    angular_dist_path = os.path.join(args.outdir, "angular_distributions.png")
    plt.savefig(angular_dist_path)
    print(f"Saved angular distributions plot to: {angular_dist_path}")
    plt.close()

    # Create 2D scatter plots for theta vs phi
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()

    for i, (particle_name, data) in enumerate(angular_data.items()):
        if len(data['theta']) > 0 and len(data['phi']) > 0:
            scatter = axes[i].scatter(data['phi'], data['theta'], alpha=0.5, s=5)
            axes[i].set_title(f"{particle_labels[particle_name]}: θ vs φ")
            axes[i].set_xlabel('φ (degrees)')
            axes[i].set_ylabel('θ (degrees)')
            axes[i].set_xlim(-180, 180)
            axes[i].set_ylim(0, 180)
            axes[i].grid(True, alpha=0.3)
            axes[i].text(0.05, 0.95, f"n={len(data['theta'])}",
                         transform=axes[i].transAxes,
                         verticalalignment='top',
                         bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
        else:
            axes[i].text(0.5, 0.5, "No data", ha='center', va='center',
                         transform=axes[i].transAxes)

    plt.tight_layout()
    scatter_plot_path = os.path.join(args.outdir, "theta_phi_scatter.png")
    plt.savefig(scatter_plot_path)
    print(f"Saved theta-phi scatter plots to: {scatter_plot_path}")
    plt.close()

    # Create 2D histograms for theta vs phi
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()

    for i, (particle_name, data) in enumerate(angular_data.items()):
        if len(data['theta']) > 0 and len(data['phi']) > 0:
            hist2d = axes[i].hist2d(data['phi'], data['theta'],
                                    bins=[36, 18], range=[[-180, 180], [0, 180]],
                                    cmap='viridis')
            cbar = plt.colorbar(hist2d[3], ax=axes[i])
            cbar.set_label('Count')
            axes[i].set_title(f"{particle_labels[particle_name]}: θ vs φ")
            axes[i].set_xlabel('φ (degrees)')
            axes[i].set_ylabel('θ (degrees)')
            axes[i].text(0.05, 0.95, f"n={len(data['theta'])}",
                         transform=axes[i].transAxes,
                         verticalalignment='top',
                         bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))
        else:
            axes[i].text(0.5, 0.5, "No data", ha='center', va='center',
                         transform=axes[i].transAxes)

    plt.tight_layout()
    heatmap_path = os.path.join(args.outdir, "theta_phi_heatmap.png")
    plt.savefig(heatmap_path)
    print(f"Saved theta-phi heatmaps to: {heatmap_path}")

    print(f"\nAll plots saved to directory: {args.outdir}")

    # If there are a large number of "unknown" particles (PDG=0), print a warning
    if 0 in pdg_counts and pdg_counts[0] > 0.1 * sum(pdg_counts.values()):
        percentage = pdg_counts[0] * 100.0 / sum(pdg_counts.values())
        print(f"\nWARNING: {pdg_counts[0]} particles ({percentage:.1f}%) have PDG=0 (unknown)")

    # Print out momenta statistics for the main particle types
    print("\nMomentum statistics:")
    for particle_name, pdg_code in particle_types.items():
        particle_mask = (all_pdgs == pdg_code)
        if ak.sum(particle_mask) == 0:
            continue

        # Flatten all momentum arrays to get values for this particle type
        flat_px = ak.flatten(px)[particle_mask]
        flat_py = ak.flatten(py)[particle_mask]
        flat_pz = ak.flatten(pz)[particle_mask]

        # Calculate p and pT
        p = np.sqrt(flat_px**2 + flat_py**2 + flat_pz**2)
        pt = np.sqrt(flat_px**2 + flat_py**2)

        print(f"\n{particle_labels[particle_name]} (PDG={pdg_code}):")
        print(f"  Count: {ak.sum(particle_mask)}")
        print(f"  Mean p: {ak.mean(p):.3f} GeV")
        print(f"  Mean pT: {ak.mean(pt):.3f} GeV")
        print(f"  Max pT: {ak.max(pt):.3f} GeV")

if __name__ == "__main__":
    main()