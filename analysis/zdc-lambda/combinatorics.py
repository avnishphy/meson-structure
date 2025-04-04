#!/usr/bin/env python3

import argparse
import numpy as np
import awkward as ak
import uproot
import hist
import vector
import matplotlib.pyplot as plt
import mplhep as hep
import re, os
from scipy.optimize import curve_fit
import warnings
from scipy.optimize import OptimizeWarning

warnings.simplefilter("ignore", OptimizeWarning)

hep.style.use("CMS")
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['savefig.facecolor'] = 'white'
plt.rcParams['savefig.bbox'] = 'tight'

############################
# 1) Command-line arguments
############################
def parse_args():
    parser = argparse.ArgumentParser(
        description="Lambda mass reconstruction analysis: combines proton+pion candidates using cartesian products."
    )
    parser.add_argument(
        "-i", "--input-files",
        nargs="+",
        required=True,
        help="List of .root files for one or more beam energies."
    )
    parser.add_argument(
        "-o", "--outdir",
        default="lambda_mass_plots",
        help="Output directory for plots"
    )
    parser.add_argument(
        "--tree",
        default="events",
        help="Name of the TTree (default: 'events')"
    )
    parser.add_argument(
        "--reco-collection",
        default="ReconstructedParticles",
        help="Name of the ReconstructedParticles collection (default: 'ReconstructedParticles')"
    )
    parser.add_argument(
        "--step-size",
        default="100MB",
        help="Chunk size for uproot.iterate (default: '100MB')"
    )
    parser.add_argument(
        "--mass-window",
        default=[1.08, 1.15],
        nargs=2,
        type=float,
        help="Mass window for Lambda signal region (default: 1.08 1.15)"
    )
    parser.add_argument(
        "--pt-cut",
        default=0.0,
        type=float,
        help="pT cut for proton and pion candidates in GeV (default: 0.0)"
    )
    parser.add_argument(
        "--vertex-cut",
        default=10.0,
        type=float,
        help="Maximum vertex distance for particle combinations in mm (default: 10.0)"
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug output"
    )
    return parser.parse_args()

############################
# 2) Utility functions
############################
def gauss(x, A, mu, sigma):
    """Simple Gaussian function for curve_fit."""
    return A * np.exp(-0.5 * ((x - mu)/sigma)**2)

def gaus_plus_pol1(x, A, mu, sigma, p0, p1):
    """Gaussian + linear background for curve_fit."""
    return A * np.exp(-0.5 * ((x - mu)/sigma)**2) + p0 + p1*x

def extract_beam_label(filename):
    """
    Parse a substring like '5x41', '10x100', '18x275' from filename,
    or return 'unknownBeam' if not found.
    """
    match = re.search(r'_(\d+x\d+)_', os.path.basename(filename))
    if match:
        return match.group(1)
    return "unknownBeam"

def vertex_distance(point1_x, point1_y, point1_z, point2_x, point2_y, point2_z):
    """Calculate 3D distance between two points."""
    return np.sqrt((point1_x - point2_x)**2 +
                   (point1_y - point2_y)**2 +
                   (point1_z - point2_z)**2)

############################
# 3) Analysis class for each beam
############################
class LambdaMassAnalysis:
    """
    Holds histograms for Lambda mass reconstruction analysis.
    Uses proton+pion combinations to reconstruct Lambda candidates.
    """
    def __init__(self, beam_label, mass_window=(1.08, 1.15), debug=False):
        self.beam_label = beam_label
        self.mass_window = mass_window
        self.lambda_pdg_mass = 1.115683  # GeV
        self.debug = debug

        # 1D histograms
        self.h_proton_count = hist.Hist(
            hist.axis.Regular(20, 0, 20, name="nproton", label="Number of proton candidates per event")
        )
        self.h_pion_count = hist.Hist(
            hist.axis.Regular(20, 0, 20, name="npion", label="Number of pion candidates per event")
        )
        self.h_lambda_count = hist.Hist(
            hist.axis.Regular(20, 0, 20, name="nlambda", label="Number of Lambda candidates per event")
        )

        # Mass histograms
        self.h_mass_all = hist.Hist(
            hist.axis.Regular(100, 1.0, 1.3, name="mass", label=r"$m_{p\pi^-}$ [GeV]")
        )
        self.h_mass_selected = hist.Hist(
            hist.axis.Regular(100, 1.0, 1.3, name="mass", label=r"$m_{p\pi^-}$ [GeV]")
        )

        # Kinematic distributions for selected candidates
        self.h_pt = hist.Hist(
            hist.axis.Regular(50, 0, 5, name="pt", label=r"$p_T^{\Lambda}$ [GeV]")
        )
        self.h_eta = hist.Hist(
            hist.axis.Regular(50, -5, 5, name="eta", label=r"$\eta^{\Lambda}$")
        )
        self.h_phi = hist.Hist(
            hist.axis.Regular(50, -np.pi, np.pi, name="phi", label=r"$\phi^{\Lambda}$ [rad]")
        )

        # 2D mass vs. kinematic variable
        self.h_mass_vs_pt = hist.Hist(
            hist.axis.Regular(50, 0, 5, name="pt", label=r"$p_T^{\Lambda}$ [GeV]"),
            hist.axis.Regular(50, 1.0, 1.3, name="mass", label=r"$m_{p\pi^-}$ [GeV]")
        )

        # For scatter plots of selected candidates
        self.proton_pt_list = []
        self.pion_pt_list = []
        self.mass_list = []
        self.pt_list = []
        self.vertex_dist_list = []

        # Counters for debugging
        self.total_events = 0
        self.events_with_protons = 0
        self.events_with_pions = 0
        self.events_with_both = 0
        self.total_protons = 0
        self.total_pions = 0
        self.total_combinations_possible = 0
        self.combinations_after_cuts = 0

    def fill_candidate_counts(self, n_protons, n_pions):
        """
        Fill histograms of proton and pion candidate counts.
        """
        # Update counters
        self.total_events += len(n_protons)
        self.events_with_protons += np.sum(n_protons > 0)
        self.events_with_pions += np.sum(n_pions > 0)
        self.events_with_both += np.sum((n_protons > 0) & (n_pions > 0))
        self.total_protons += np.sum(n_protons)
        self.total_pions += np.sum(n_pions)

        # Fill histograms
        self.h_proton_count.fill(nproton=n_protons)
        self.h_pion_count.fill(npion=n_pions)

        # Debug information
        if self.debug:
            print(f"Events processed: {len(n_protons)}")
            print(f"Events with protons: {np.sum(n_protons > 0)}")
            print(f"Events with pions: {np.sum(n_pions > 0)}")
            print(f"Events with both: {np.sum((n_protons > 0) & (n_pions > 0))}")
            print(f"Total protons: {np.sum(n_protons)}")
            print(f"Total pions: {np.sum(n_pions)}")
            print(f"Proton counts: {n_protons[:20]}...")
            print(f"Pion counts: {n_pions[:20]}...")

    def process_combinations(self,
                             pdg, px, py, pz, energy,
                             ref_x, ref_y, ref_z,
                             pt_cut=0.0, vertex_cut=float('inf')):
        """
        Process combinations directly from the full event arrays.
        This avoids issues with conversion to Python lists.
        """
        # Create masks for proton and pion
        proton_mask = (pdg == 2212)  # Proton PDG code
        pion_mask = (pdg == -211)    # Pi- PDG code

        n_protons = ak.sum(proton_mask, axis=1)
        n_pions = ak.sum(pion_mask, axis=1)

        # Fill histograms for candidate counts
        self.fill_candidate_counts(ak.to_numpy(n_protons), ak.to_numpy(n_pions))

        # Print debug info about a few protons and pions
        if self.debug:
            # Find an event with at least one proton and one pion
            for i in range(min(10, len(proton_mask))):
                n_p = ak.sum(proton_mask[i])
                n_pi = ak.sum(pion_mask[i])

                if n_p > 0 and n_pi > 0:
                    print(f"\nExample event {i}:")
                    print(f"  {n_p} protons, {n_pi} pions")

                    # Get first proton
                    p_idx = ak.argmax(proton_mask[i], axis=0, keepdims=False)
                    print(f"  First proton (idx={p_idx}):")
                    print(f"    PDG = {pdg[i][p_idx]}")
                    print(f"    px = {px[i][p_idx]}, py = {py[i][p_idx]}, pz = {pz[i][p_idx]}, E = {energy[i][p_idx]}")
                    print(f"    ref_x = {ref_x[i][p_idx]}, ref_y = {ref_y[i][p_idx]}, ref_z = {ref_z[i][p_idx]}")

                    # Get first pion
                    pi_idx = ak.argmax(pion_mask[i], axis=0, keepdims=False)
                    print(f"  First pion (idx={pi_idx}):")
                    print(f"    PDG = {pdg[i][pi_idx]}")
                    print(f"    px = {px[i][pi_idx]}, py = {py[i][pi_idx]}, pz = {pz[i][pi_idx]}, E = {energy[i][pi_idx]}")
                    print(f"    ref_x = {ref_x[i][pi_idx]}, ref_y = {ref_y[i][pi_idx]}, ref_z = {ref_z[i][pi_idx]}")

                    break

        # Create a selection that finds events with at least one proton and one pion
        events_with_both = (n_protons > 0) & (n_pions > 0)

        # Filter to only events with both particles
        if not ak.any(events_with_both):
            if self.debug:
                print("No events with both protons and pions found.")
            return 0

        # Extract proton and pion data
        proton_px = px[proton_mask]
        proton_py = py[proton_mask]
        proton_pz = pz[proton_mask]
        proton_E = energy[proton_mask]
        proton_ref_x = ref_x[proton_mask]
        proton_ref_y = ref_y[proton_mask]
        proton_ref_z = ref_z[proton_mask]

        pion_px = px[pion_mask]
        pion_py = py[pion_mask]
        pion_pz = pz[pion_mask]
        pion_E = energy[pion_mask]
        pion_ref_x = ref_x[pion_mask]
        pion_ref_y = ref_y[pion_mask]
        pion_ref_z = ref_z[pion_mask]

        # Debug info about extracted particles
        if self.debug:
            print("\nExtracted particles:")
            print(f"  Protons: {len(proton_px)} with shape {proton_px.type if hasattr(proton_px, 'type') else 'unknown'}")
            print(f"  Pions: {len(pion_px)} with shape {pion_px.type if hasattr(pion_px, 'type') else 'unknown'}")

            # Show a few values
            if len(proton_px) > 0:
                print(f"  First few proton px values: {proton_px[:5] if len(proton_px) >= 5 else proton_px}")
            if len(pion_px) > 0:
                print(f"  First few pion px values: {pion_px[:5] if len(pion_px) >= 5 else pion_px}")

        # Try a simpler direct approach - iterate through events and process manually
        total_combinations = 0
        mass_list = []
        pt_list = []
        eta_list = []
        phi_list = []
        vertex_dist_list = []
        proton_pt_list = []
        pion_pt_list = []

        for evt_idx in range(len(pdg)):
            # Skip events without both particles
            if not events_with_both[evt_idx]:
                continue

            # Get protons and pions for this event
            evt_proton_mask = proton_mask[evt_idx]
            evt_pion_mask = pion_mask[evt_idx]

            # Skip if no matches (redundant but safe)
            if not ak.any(evt_proton_mask) or not ak.any(evt_pion_mask):
                continue

            # Get proton indices and properties
            proton_indices = ak.where(evt_proton_mask)[0]
            pion_indices = ak.where(evt_pion_mask)[0]

            # Debug info
            if self.debug and evt_idx < 3:
                print(f"\nEvent {evt_idx}:")
                print(f"  Proton indices: {proton_indices}")
                print(f"  Pion indices: {pion_indices}")

            # Iterate through all combinations
            for p_idx in proton_indices:
                p_px = px[evt_idx][p_idx]
                p_py = py[evt_idx][p_idx]
                p_pz = pz[evt_idx][p_idx]
                p_E = energy[evt_idx][p_idx]
                p_ref_x = ref_x[evt_idx][p_idx]
                p_ref_y = ref_y[evt_idx][p_idx]
                p_ref_z = ref_z[evt_idx][p_idx]

                # Calculate proton pT
                p_pt = np.sqrt(p_px**2 + p_py**2)

                # Apply pT cut
                if p_pt < pt_cut:
                    continue

                for pi_idx in pion_indices:
                    pi_px = px[evt_idx][pi_idx]
                    pi_py = py[evt_idx][pi_idx]
                    pi_pz = pz[evt_idx][pi_idx]
                    pi_E = energy[evt_idx][pi_idx]
                    pi_ref_x = ref_x[evt_idx][pi_idx]
                    pi_ref_y = ref_y[evt_idx][pi_idx]
                    pi_ref_z = ref_z[evt_idx][pi_idx]

                    # Debug the first combination
                    if self.debug and evt_idx < 3 and len(mass_list) == 0:
                        print(f"\nFirst combination (event {evt_idx}, proton {p_idx}, pion {pi_idx}):")
                        print(f"  Proton: px={p_px}, py={p_py}, pz={p_pz}, E={p_E}")
                        print(f"  Pion: px={pi_px}, py={pi_py}, pz={pi_pz}, E={pi_E}")

                    # Calculate pion pT
                    pi_pt = np.sqrt(pi_px**2 + pi_py**2)

                    # Apply pT cut
                    if pi_pt < pt_cut:
                        continue

                    # Calculate vertex distance
                    vtx_dist = np.sqrt((p_ref_x - pi_ref_x)**2 +
                                       (p_ref_y - pi_ref_y)**2 +
                                       (p_ref_z - pi_ref_z)**2)

                    # Apply vertex cut
                    if vtx_dist > vertex_cut:
                        continue

                    # Create 4-vectors and calculate mass
                    try:
                        p_vec = vector.obj(px=p_px, py=p_py, pz=p_pz, E=p_E)
                        pi_vec = vector.obj(px=pi_px, py=pi_py, pz=pi_pz, E=pi_E)
                        lambda_vec = p_vec + pi_vec

                        # Get properties
                        mass = lambda_vec.mass
                        pt = lambda_vec.pt
                        eta = lambda_vec.eta
                        phi = lambda_vec.phi

                        # Debug the first successful combination
                        if self.debug and total_combinations == 0:
                            print(f"\nFirst successful combination:")
                            print(f"  Lambda mass: {mass}")
                            print(f"  Lambda pT: {pt}")
                            print(f"  Lambda eta: {eta}")
                            print(f"  Lambda phi: {phi}")
                            print(f"  Vertex distance: {vtx_dist}")

                        # Store properties
                        mass_list.append(mass)
                        pt_list.append(pt)
                        eta_list.append(eta)
                        phi_list.append(phi)
                        vertex_dist_list.append(vtx_dist)
                        proton_pt_list.append(p_pt)
                        pion_pt_list.append(pi_pt)

                        # Count this combination
                        total_combinations += 1

                    except Exception as e:
                        if self.debug:
                            print(f"Error creating vectors: {e}")
                        continue

        # Fill histograms with the collected data
        if total_combinations > 0:
            # Convert lists to numpy arrays
            masses = np.array(mass_list)
            pts = np.array(pt_list)
            etas = np.array(eta_list)
            phis = np.array(phi_list)

            # Update histograms
            self.h_mass_all.fill(mass=masses)

            # Signal region for kinematic histograms
            signal_mask = (masses > self.mass_window[0]) & (masses < self.mass_window[1])

            if np.any(signal_mask):
                signal_masses = masses[signal_mask]
                signal_pts = pts[signal_mask]
                signal_etas = etas[signal_mask]
                signal_phis = phis[signal_mask]

                self.h_mass_selected.fill(mass=signal_masses)
                self.h_pt.fill(pt=signal_pts)
                self.h_eta.fill(eta=signal_etas)
                self.h_phi.fill(phi=signal_phis)

                # Fill 2D histogram
                for m, p in zip(signal_masses, signal_pts):
                    self.h_mass_vs_pt.fill(pt=p, mass=m)

                # Store for scatter plots
                self.mass_list.extend(signal_masses)
                self.pt_list.extend(signal_pts)
                self.proton_pt_list.extend(np.array(proton_pt_list)[signal_mask])
                self.pion_pt_list.extend(np.array(pion_pt_list)[signal_mask])
                self.vertex_dist_list.extend(np.array(vertex_dist_list)[signal_mask])

        # Update counters
        self.total_combinations_possible += total_combinations
        self.combinations_after_cuts += total_combinations

        if self.debug:
            print(f"\nProcessed {total_combinations} combinations in this chunk")

        return total_combinations

    def print_summary(self):
        """Print summary statistics"""
        print("\n=== Analysis Summary ===")
        print(f"Total events processed: {self.total_events}")
        print(f"Events with protons: {self.events_with_protons} ({self.events_with_protons/max(1,self.total_events)*100:.1f}%)")
        print(f"Events with pions: {self.events_with_pions} ({self.events_with_pions/max(1,self.total_events)*100:.1f}%)")
        print(f"Events with both: {self.events_with_both} ({self.events_with_both/max(1,self.total_events)*100:.1f}%)")
        print(f"Total protons found: {self.total_protons}")
        print(f"Total pions found: {self.total_pions}")
        print(f"Average protons per event: {self.total_protons/max(1,self.total_events):.2f}")
        print(f"Average pions per event: {self.total_pions/max(1,self.total_events):.2f}")
        print(f"Total combinations processed: {self.total_combinations_possible}")
        print(f"Combinations after cuts: {self.combinations_after_cuts} ({self.combinations_after_cuts/max(1,self.total_combinations_possible)*100:.2f}%)")
        print(f"Histogram entries (all mass): {sum(self.h_mass_all.values())}")
        print(f"Histogram entries (selected mass): {sum(self.h_mass_selected.values())}")

    def finalize_and_plot(self, outdir):
        """Create plots and save them in the specified output directory."""
        # Print summary statistics
        self.print_summary()

        beam_outdir = os.path.join(outdir, self.beam_label)
        os.makedirs(beam_outdir, exist_ok=True)

        # 1) Candidate counts
        fig, axs = plt.subplots(1, 3, figsize=(18, 6))
        fig.suptitle(f"Beam: {self.beam_label} - Particle Candidate Counts")

        # Proton count
        plt.sca(axs[0])
        vals = self.h_proton_count.values()
        edges = self.h_proton_count.axes[0].edges
        plt.hist(edges[:-1], bins=edges, weights=vals, histtype='step', color='blue', lw=2)
        plt.yscale('log')
        plt.ylim(0.5)
        plt.xlabel("Number of Proton Candidates")
        plt.ylabel("Events")

        # Pion count
        plt.sca(axs[1])
        vals = self.h_pion_count.values()
        edges = self.h_pion_count.axes[0].edges
        plt.hist(edges[:-1], bins=edges, weights=vals, histtype='step', color='red', lw=2)
        plt.yscale('log')
        plt.ylim(0.5)
        plt.xlabel("Number of Pion Candidates")
        plt.ylabel("Events")

        # Lambda candidates
        plt.sca(axs[2])
        vals = self.h_lambda_count.values()
        edges = self.h_lambda_count.axes[0].edges
        plt.hist(edges[:-1], bins=edges, weights=vals, histtype='step', color='green', lw=2)
        plt.yscale('log')
        plt.ylim(0.5)
        plt.xlabel("Number of Lambda Candidates")
        plt.ylabel("Events")

        plt.tight_layout()
        base_fname = os.path.join(beam_outdir, f"candidate_counts_{self.beam_label}")
        plt.savefig(f"{base_fname}.pdf")
        plt.savefig(f"{base_fname}.png", dpi=300)
        plt.close()

        # 2) Mass distributions with fit
        fig, axs = plt.subplots(1, 2, figsize=(16, 8))
        fig.suptitle(f"Beam: {self.beam_label} - Lambda Mass Reconstruction")

        # All combinations
        plt.sca(axs[0])
        vals_all = self.h_mass_all.values()
        edges = self.h_mass_all.axes[0].edges
        plt.hist(edges[:-1], bins=edges, weights=vals_all,
                 histtype='step', color='black', lw=2, label="All combinations")
        plt.axvline(self.lambda_pdg_mass, ls='--', color='green',
                    label=f"PDG mass: {self.lambda_pdg_mass:.6f} GeV")
        plt.legend()
        plt.xlabel(r"$m_{p\pi^-}$ [GeV]")
        plt.ylabel("Combinations / bin")

        # Selected combinations with fit
        plt.sca(axs[1])
        vals_sel = self.h_mass_selected.values()
        edges = self.h_mass_selected.axes[0].edges
        plt.hist(edges[:-1], bins=edges, weights=vals_sel,
                 histtype='stepfilled', color='skyblue', edgecolor='blue', label="Selected")
        plt.axvline(self.lambda_pdg_mass, ls='--', color='green',
                    label=f"PDG mass: {self.lambda_pdg_mass:.6f} GeV")

        # Fit with Gaussian + linear background
        bins = 0.5 * (edges[:-1] + edges[1:])
        try:
            # Initial guesses based on histogram properties
            max_bin = np.argmax(vals_sel)
            peak_pos = bins[max_bin]
            peak_height = vals_sel[max_bin]

            # Reasonable initial params for the fit
            p0 = [peak_height, peak_pos, 0.01, 0.1*peak_height, 0]

            # Fit only the region around the peak
            mask = (bins > peak_pos - 0.1) & (bins < peak_pos + 0.1)
            if np.sum(mask) > 5:  # Ensure enough points for fitting
                popt, _ = curve_fit(gaus_plus_pol1, bins[mask], vals_sel[mask], p0=p0)

                # Plot the fit
                x_fit = np.linspace(1.0, 1.3, 1000)
                plt.plot(x_fit, gaus_plus_pol1(x_fit, *popt), 'r-',
                         label=f"Fit: $\mu={popt[1]:.4f}$, $\sigma={popt[2]:.4f}$ GeV")
        except Exception as e:
            print(f"Fitting failed: {e}")

        plt.legend()
        plt.xlabel(r"$m_{p\pi^-}$ [GeV]")
        plt.ylabel("Combinations / bin")

        plt.tight_layout()
        base_fname = os.path.join(beam_outdir, f"lambda_mass_{self.beam_label}")
        plt.savefig(f"{base_fname}.pdf")
        plt.savefig(f"{base_fname}.png", dpi=300)
        plt.close()

        # Only proceed with other plots if we have combinations
        if self.total_combinations_possible == 0:
            return

        # 3) Kinematic distributions
        fig, axs = plt.subplots(1, 3, figsize=(18, 6))
        fig.suptitle(f"Beam: {self.beam_label} - Lambda Candidate Kinematics")

        # pT distribution
        plt.sca(axs[0])
        vals = self.h_pt.values()
        edges = self.h_pt.axes[0].edges
        plt.hist(edges[:-1], bins=edges, weights=vals, histtype='step', color='blue', lw=2)
        plt.xlabel(r"$p_T^{\Lambda}$ [GeV]")
        plt.ylabel("Candidates / bin")

        # eta distribution
        plt.sca(axs[1])
        vals = self.h_eta.values()
        edges = self.h_eta.axes[0].edges
        plt.hist(edges[:-1], bins=edges, weights=vals, histtype='step', color='red', lw=2)
        plt.xlabel(r"$\eta^{\Lambda}$")
        plt.ylabel("Candidates / bin")

        # phi distribution
        plt.sca(axs[2])
        vals = self.h_phi.values()
        edges = self.h_phi.axes[0].edges
        plt.hist(edges[:-1], bins=edges, weights=vals, histtype='step', color='green', lw=2)
        plt.xlabel(r"$\phi^{\Lambda}$ [rad]")
        plt.ylabel("Candidates / bin")

        plt.tight_layout()
        base_fname = os.path.join(beam_outdir, f"lambda_kinematics_{self.beam_label}")
        plt.savefig(f"{base_fname}.pdf")
        plt.savefig(f"{base_fname}.png", dpi=300)
        plt.close()

        # 4) 2D mass vs. pT
        plt.figure(figsize=(10, 8))
        h = self.h_mass_vs_pt
        hist_array = h.values()

        plt.pcolormesh(h.axes[0].edges, h.axes[1].edges, hist_array.T,
                       cmap='viridis', shading='auto')
        plt.colorbar(label='Candidates')
        plt.xlabel(r"$p_T^{\Lambda}$ [GeV]")
        plt.ylabel(r"$m_{p\pi^-}$ [GeV]")
        plt.axhline(self.lambda_pdg_mass, ls='--', color='red', alpha=0.7)
        plt.title(f"Beam: {self.beam_label} - Mass vs. pT")

        base_fname = os.path.join(beam_outdir, f"mass_vs_pt_{self.beam_label}")
        plt.savefig(f"{base_fname}.pdf")
        plt.savefig(f"{base_fname}.png", dpi=300)
        plt.close()

        # 5) Scatter plots
        if self.proton_pt_list and self.pion_pt_list:
            proton_pt = np.array(self.proton_pt_list)
            pion_pt = np.array(self.pion_pt_list)

            if len(proton_pt) > 0 and len(pion_pt) > 0:
                plt.figure(figsize=(10, 8))
                plt.scatter(proton_pt, pion_pt, s=2, alpha=0.5)
                plt.xlabel(r"Proton $p_T$ [GeV]")
                plt.ylabel(r"Pion $p_T$ [GeV]")
                plt.title(f"Beam: {self.beam_label} - p vs. π pT for Λ candidates")

                base_fname = os.path.join(beam_outdir, f"proton_pion_pt_{self.beam_label}")
                plt.savefig(f"{base_fname}.pdf")
                plt.savefig(f"{base_fname}.png", dpi=300)
                plt.close()

        # 6) Vertex distance for selected candidates
        if self.vertex_dist_list:
            vertex_dist = np.array(self.vertex_dist_list)

            if len(vertex_dist) > 0:
                plt.figure(figsize=(10, 8))
                plt.hist(vertex_dist, bins=50, range=(0, 20))
                plt.xlabel("Vertex Distance [mm]")
                plt.ylabel("Candidates / bin")
                plt.title(f"Beam: {self.beam_label} - Vertex Distance for Λ candidates")

                base_fname = os.path.join(beam_outdir, f"vertex_distance_{self.beam_label}")
                plt.savefig(f"{base_fname}.pdf")
                plt.savefig(f"{base_fname}.png", dpi=300)
                plt.close()

############################
# 4) Main function
############################
def main():
    args = parse_args()
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    # Group input files by beam label
    beams_map = {}
    for f in args.input_files:
        beam = extract_beam_label(f)
        beams_map.setdefault(beam, []).append(f)

    # Build the branch names for the specific collection
    collection = args.reco_collection
    branch_prefix = f"{collection}."

    needed_branches = [
        f"{branch_prefix}PDG",
        f"{branch_prefix}charge",
        f"{branch_prefix}momentum.x",
        f"{branch_prefix}momentum.y",
        f"{branch_prefix}momentum.z",
        f"{branch_prefix}energy",
        f"{branch_prefix}referencePoint.x",
        f"{branch_prefix}referencePoint.y",
        f"{branch_prefix}referencePoint.z",
    ]

    beam_analyses = {}

    for beam_label, filelist in beams_map.items():
        print(f"\n=== Beam: {beam_label} ===")
        print("   Files:", filelist)
        beam_analyses[beam_label] = LambdaMassAnalysis(
            beam_label,
            mass_window=args.mass_window,
            debug=args.debug
        )
        analysis_obj = beam_analyses[beam_label]

        file_dict = {fname: args.tree for fname in filelist}

        # Print all available branches in first file if debugging
        if args.debug:
            print("\nAvailable branches:")
            f = uproot.open(filelist[0])
            tree = f[args.tree]
            for branch in sorted(tree.keys()):
                print(f"  {branch}")

        chunk_counter = 0

        for chunk in uproot.iterate(
                file_dict,
                expressions=needed_branches,
                step_size=args.step_size,
                library="ak"
        ):
            chunk_counter += 1
            if chunk_counter == 1 or args.debug:
                print(f"[Debug] Processing chunk #{chunk_counter} with {len(chunk)} events")

            # Extract relevant branches
            pdg = chunk[f"{branch_prefix}PDG"]
            charge = chunk[f"{branch_prefix}charge"]
            px = chunk[f"{branch_prefix}momentum.x"]
            py = chunk[f"{branch_prefix}momentum.y"]
            pz = chunk[f"{branch_prefix}momentum.z"]
            energy = chunk[f"{branch_prefix}energy"]
            ref_x = chunk[f"{branch_prefix}referencePoint.x"]
            ref_y = chunk[f"{branch_prefix}referencePoint.y"]
            ref_z = chunk[f"{branch_prefix}referencePoint.z"]

            # Print PDG codes distribution if debugging
            if args.debug and chunk_counter == 1:
                all_pdgs = ak.flatten(pdg)
                unique_pdgs, counts = np.unique(ak.to_numpy(all_pdgs), return_counts=True)
                print("\nPDG code distribution:")
                for pdg_code, count in zip(unique_pdgs, counts):
                    print(f"  PDG {pdg_code}: {count} particles")

            # Process this chunk
            n_combinations = analysis_obj.process_combinations(
                pdg, px, py, pz, energy,
                ref_x, ref_y, ref_z,
                pt_cut=args.pt_cut,
                vertex_cut=args.vertex_cut
            )

            print(f"[Debug] Processed {n_combinations} combinations in chunk #{chunk_counter}")

    # Create final plots
    for beam_label, analysis_obj in beam_analyses.items():
        analysis_obj.finalize_and_plot(outdir)

    print("\nDone! Output saved to:", outdir)

if __name__ == "__main__":
    main()