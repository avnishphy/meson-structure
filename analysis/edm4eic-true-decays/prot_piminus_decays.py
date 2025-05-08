#!/usr/bin/env python3
"""
# Lambda Decay Analyzer - Visualization

This script analyzes the CSV files created by the Lambda decay extractor and generates
histograms and 2D plots of particle endpoint distributions.

It creates:
1. 1D histograms of endpoint_z for Lambda, proton, and pion
2. 2D histograms of endpoint_z vs endpoint_y for each particle type

Usage:
    python lambda_decay_visualizer.py -i INPUT_PREFIX -o OUTPUT_FOLDER
"""

import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from hist import Hist


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Analyze Lambda decay data and create histograms.")
    parser.add_argument("-i", "--input-prefix", required=True,
                        help="Prefix of input CSV files (e.g., 'lambda_decays')")
    parser.add_argument("-o", "--output-folder", required=True,
                        help="Folder to save output plots")
    return parser.parse_args()


def load_data(input_prefix):
    """Load the CSV data files.

    Args:
        input_prefix: Prefix of the CSV files

    Returns:
        Dictionary of DataFrames
    """
    data = {}

    # Try to load proton_piminus data (Lambda → p + π-)
    proton_piminus_file = f"{input_prefix}_proton_piminus.csv"
    if os.path.exists(proton_piminus_file):
        data['proton_piminus'] = pd.read_csv(proton_piminus_file)
        print(f"Loaded {len(data['proton_piminus'])} Lambda → p + π- decays")

    # Try to load neutron_pizero data (Lambda → n + π0)
    neutron_pizero_file = f"{input_prefix}_neutron_pizero.csv"
    if os.path.exists(neutron_pizero_file):
        data['neutron_pizero'] = pd.read_csv(neutron_pizero_file)
        print(f"Loaded {len(data['neutron_pizero'])} Lambda → n + π0 decays")

    # Try to load no_decay data
    no_decay_file = f"{input_prefix}_no_decay.csv"
    if os.path.exists(no_decay_file):
        data['no_decay'] = pd.read_csv(no_decay_file)
        print(f"Loaded {len(data['no_decay'])} Lambda with no decay")

    return data


def create_1d_histograms(data, output_folder):
    """Create 1D histograms of endpoint_z for each particle.

    Args:
        data: Dictionary of DataFrames
        output_folder: Folder to save plots
    """
    os.makedirs(output_folder, exist_ok=True)

    # Process Lambda → p + π- data
    if 'proton_piminus' in data:
        df = data['proton_piminus']

        # Lambda endpoint_z
        lambda_endz = df['lam_endpoint_z'].values / 1000.0  # Convert mm to meters
        lambda_endz = lambda_endz[(lambda_endz >= -5) & (lambda_endz <= 40)]  # Filter to requested range

        # Create a histogram
        h_lambda = Hist(
            np.linspace(-5, 40, 100, name="endpoint_z (m)")
        )

        # Fill the histogram
        h_lambda.fill(lambda_endz)

        # Plot the histogram
        fig, ax = plt.subplots(figsize=(10, 6))
        h_lambda.plot(ax=ax)
        ax.set_title("Lambda Endpoint Z Distribution")
        ax.set_xlabel("Endpoint Z [m]")
        ax.set_ylabel("Counts")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(os.path.join(output_folder, "lambda_endpoint_z.png"), dpi=300)
        plt.close(fig)

        # Proton endpoint_z
        proton_endz = df['prot_endpoint_z'].values / 1000.0  # Convert mm to meters
        proton_endz = proton_endz[(proton_endz >= -5) & (proton_endz <= 40)]  # Filter to requested range

        # Create a histogram
        h_proton = Hist(
            np.linspace(-5, 40, 100, name="endpoint_z (m)")
        )

        # Fill the histogram
        h_proton.fill(proton_endz)

        # Plot the histogram
        fig, ax = plt.subplots(figsize=(10, 6))
        h_proton.plot(ax=ax)
        ax.set_title("Proton Endpoint Z Distribution")
        ax.set_xlabel("Endpoint Z [m]")
        ax.set_ylabel("Counts")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(os.path.join(output_folder, "proton_endpoint_z.png"), dpi=300)
        plt.close(fig)

        # Pion endpoint_z
        pion_endz = df['piminus_endpoint_z'].values / 1000.0  # Convert mm to meters
        pion_endz = pion_endz[(pion_endz >= -5) & (pion_endz <= 40)]  # Filter to requested range

        # Create a histogram
        h_pion = Hist(
            np.linspace(-5, 40, 100, name="endpoint_z (m)")
        )

        # Fill the histogram
        h_pion.fill(pion_endz)

        # Plot the histogram
        fig, ax = plt.subplots(figsize=(10, 6))
        h_pion.plot(ax=ax)
        ax.set_title("π- Endpoint Z Distribution")
        ax.set_xlabel("Endpoint Z [m]")
        ax.set_ylabel("Counts")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(os.path.join(output_folder, "pion_endpoint_z.png"), dpi=300)
        plt.close(fig)

        # Compare all three on one plot
        fig, ax = plt.subplots(figsize=(12, 7))

        # Create histograms with the same binning
        bins = np.linspace(-5, 40, 100)

        # Plot histograms
        ax.hist(lambda_endz, bins=bins, alpha=0.5, label="Lambda")
        ax.hist(proton_endz, bins=bins, alpha=0.5, label="Proton")
        ax.hist(pion_endz, bins=bins, alpha=0.5, label="π-")

        ax.set_title("Comparison of Endpoint Z Distributions")
        ax.set_xlabel("Endpoint Z [m]")
        ax.set_ylabel("Counts")
        ax.grid(True, alpha=0.3)
        ax.legend()

        fig.tight_layout()
        fig.savefig(os.path.join(output_folder, "endpoint_z_comparison.png"), dpi=300)
        plt.close(fig)


def create_2d_histograms(data, output_folder):
    """Create 2D histograms of endpoint_z vs endpoint_y for each particle.

    Args:
        data: Dictionary of DataFrames
        output_folder: Folder to save plots
    """
    os.makedirs(output_folder, exist_ok=True)

    # Process Lambda → p + π- data
    if 'proton_piminus' in data:
        df = data['proton_piminus']

        # Lambda endpoint_z vs endpoint_y (convert mm to meters)
        lambda_endz = df['lam_endpoint_z'].values / 1000.0
        lambda_endy = df['lam_endpoint_y'].values / 1000.0

        # Filter to requested z range and remove NaN values
        mask = (lambda_endz >= -5) & (lambda_endz <= 40) & np.isfinite(lambda_endz) & np.isfinite(lambda_endy)
        lambda_endz = lambda_endz[mask]
        lambda_endy = lambda_endy[mask]

        # Create a 2D histogram
        h_lambda_2d = Hist(
            np.linspace(-5, 40, 100, name="endpoint_z (m)"),
            np.linspace(np.min(lambda_endy), np.max(lambda_endy), 100, name="endpoint_y (m)")
        )

        # Fill the histogram
        h_lambda_2d.fill(lambda_endz, lambda_endy)

        # Plot the 2D histogram
        fig, ax = plt.subplots(figsize=(10, 8))
        mesh = ax.pcolormesh(
            h_lambda_2d.axes[0].edges,
            h_lambda_2d.axes[1].edges,
            h_lambda_2d.values().T,
            cmap='viridis',
            norm='log'
        )
        fig.colorbar(mesh, ax=ax, label="Counts")
        ax.set_title("Lambda Endpoint Z vs Y")
        ax.set_xlabel("Endpoint Z [m]")
        ax.set_ylabel("Endpoint Y [m]")
        fig.tight_layout()
        fig.savefig(os.path.join(output_folder, "lambda_endpoint_zy.png"), dpi=300)
        plt.close(fig)

        # Proton endpoint_z vs endpoint_y (convert mm to meters)
        proton_endz = df['prot_endpoint_z'].values / 1000.0
        proton_endy = df['prot_endpoint_y'].values / 1000.0

        # Filter to requested z range and remove NaN values
        mask = (proton_endz >= -5) & (proton_endz <= 40) & np.isfinite(proton_endz) & np.isfinite(proton_endy)
        proton_endz = proton_endz[mask]
        proton_endy = proton_endy[mask]

        # Create a 2D histogram
        h_proton_2d = Hist(
            np.linspace(-5, 40, 100, name="endpoint_z (m)"),
            np.linspace(np.min(proton_endy), np.max(proton_endy), 100, name="endpoint_y (m)")
        )

        # Fill the histogram
        h_proton_2d.fill(proton_endz, proton_endy)

        # Plot the 2D histogram
        fig, ax = plt.subplots(figsize=(10, 8))
        mesh = ax.pcolormesh(
            h_proton_2d.axes[0].edges,
            h_proton_2d.axes[1].edges,
            h_proton_2d.values().T,
            cmap='viridis',
            norm='log'
        )
        fig.colorbar(mesh, ax=ax, label="Counts")
        ax.set_title("Proton Endpoint Z vs Y")
        ax.set_xlabel("Endpoint Z [m]")
        ax.set_ylabel("Endpoint Y [m]")
        fig.tight_layout()
        fig.savefig(os.path.join(output_folder, "proton_endpoint_zy.png"), dpi=300)
        plt.close(fig)

        # Pion endpoint_z vs endpoint_y (convert mm to meters)
        pion_endz = df['piminus_endpoint_z'].values / 1000.0
        pion_endy = df['piminus_endpoint_y'].values / 1000.0

        # Filter to requested z range and remove NaN values
        mask = (pion_endz >= -5) & (pion_endz <= 40) & np.isfinite(pion_endz) & np.isfinite(pion_endy)
        pion_endz = pion_endz[mask]
        pion_endy = pion_endy[mask]

        # Create a 2D histogram
        h_pion_2d = Hist(
            np.linspace(-5, 40, 100, name="endpoint_z (m)"),
            np.linspace(np.min(pion_endy), np.max(pion_endy), 100, name="endpoint_y (m)")
        )

        # Fill the histogram
        h_pion_2d.fill(pion_endz, pion_endy)

        # Plot the 2D histogram
        fig, ax = plt.subplots(figsize=(10, 8))
        mesh = ax.pcolormesh(
            h_pion_2d.axes[0].edges,
            h_pion_2d.axes[1].edges,
            h_pion_2d.values().T,
            cmap='viridis',
            norm='log'
        )
        fig.colorbar(mesh, ax=ax, label="Counts")
        ax.set_title("π- Endpoint Z vs Y")
        ax.set_xlabel("Endpoint Z [m]")
        ax.set_ylabel("Endpoint Y [m]")
        fig.tight_layout()
        fig.savefig(os.path.join(output_folder, "pion_endpoint_zy.png"), dpi=300)
        plt.close(fig)


def main():
    """Main function."""
    args = parse_args()

    # Load data
    data = load_data(args.input_prefix)

    if not data:
        print("No data files found. Please check the input prefix.")
        return

    # Create output folder if it doesn't exist
    os.makedirs(args.output_folder, exist_ok=True)

    # Create histograms
    create_1d_histograms(data, args.output_folder)
    create_2d_histograms(data, args.output_folder)

    print(f"Analysis complete. Plots saved to {args.output_folder}")


if __name__ == "__main__":
    main()