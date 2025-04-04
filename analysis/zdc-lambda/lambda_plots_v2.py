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
from   scipy.optimize import curve_fit
import warnings
from   scipy.optimize import OptimizeWarning

warnings.simplefilter("ignore", OptimizeWarning)

hep.style.use("CMS")
plt.rcParams['figure.facecolor']    = 'white'
plt.rcParams['savefig.facecolor']   = 'white'
plt.rcParams['savefig.bbox']        = 'tight'


############################
# 1) Command-line arguments
############################
def parse_args():
    parser = argparse.ArgumentParser(
        description="Lambda recon analysis: chunked uproot + Scikit-HEP hist, PDG-based selection, vectorized CM angles."
    )
    parser.add_argument(
        "-i", "--input-files",
        nargs="+",
        required=True,
        help="List of .root files for one or more beam energies."
    )
    parser.add_argument(
        "-o", "--outdir",
        default="chunked_plots",
        help="Output directory for plots"
    )
    parser.add_argument(
        "--tree",
        default="events",
        help="Name of the TTree (default: 'events')"
    )
    parser.add_argument(
        "--step-size",
        default="100MB",
        help="Chunk size for uproot.iterate (default: '100MB')"
    )
    parser.add_argument(
        "--no-cm-angles",
        action="store_true",
        help="Skip the neutron CM-angle resolution if you don't have that collection."
    )
    return parser.parse_args()


############################
# 2) Utility functions
############################
def gauss(x, A, mu, sigma):
    """Simple Gaussian function for curve_fit."""
    return A * np.exp(-0.5 * ((x - mu)/sigma)**2)

def extract_beam_label(filename):
    """
    Parse a substring like '5x41', '10x100', '18x275' from filename,
    or return 'unknownBeam' if not found.
    """
    match = re.search(r'_(\d+x\d+)_', os.path.basename(filename))
    if match:
        return match.group(1)
    return "unknownBeam"


############################
# 3) Analysis class for each beam
############################
class LambdaAnalysis:
    """
    Holds histograms (h_nclusters, h_dtheta, h_dzvtx, h_mass) and arrays
    for final scatter plots (theta, z). Also stores cm_dphi, cm_dtheta
    for neutron resolution.
    """
    def __init__(self, beam_label):
        self.beam_label = beam_label

        # 1D hist objects
        self.h_nclusters = hist.Hist(
            hist.axis.Regular(20, 0, 20, name="nclust", label="ZDC cluster count")
        )
        self.h_dtheta = hist.Hist(
            hist.axis.Regular(50, -1, 1, name="dtheta",
                              label=r"$\theta_{\mathrm{rec}}-\theta_{\mathrm{truth}}$ [mrad]")
        )
        self.h_dzvtx = hist.Hist(
            hist.axis.Regular(50, -10, 10, name="dzvtx",
                              label=r"$z^{\mathrm{rec}}-\;z^{\mathrm{truth}}$ [m]")
        )
        self.h_mass = hist.Hist(
            hist.axis.Regular(100, 1.0, 1.25, name="mass", label=r"$m_{\Lambda}^{\mathrm{rec}}$ [GeV]")
        )

        # For scatter plots (angles, z-vertex)
        self.truth_theta_list = []
        self.recon_theta_list = []
        self.z_truth_list = []
        self.z_recon_list = []

        # CM angle resolution: store final arrays
        self.cm_dphi   = []
        self.cm_dtheta = []

    #######################################
    # Fill methods
    #######################################
    def fill_cluster_count(self, cluster_x):
        """Histogram of ZDC cluster count."""
        nclust = ak.num(cluster_x, axis=1)
        self.h_nclusters.fill(nclust=nclust)

    def fill_truth_recon_angles(self,
                                px_lambda_truth, py_lambda_truth, pz_lambda_truth,
                                px_lambda_recon, py_lambda_recon, pz_lambda_recon,
                                tilt=-0.025):
        """
        Fill angle scatter & residual for events with a truth-lambda and at least 1 recon-lambda.
        """
        valid_truth = ~ak.is_none(px_lambda_truth)
        has_rec = (ak.num(px_lambda_recon, axis=1) > 0)
        fill_mask = valid_truth & has_rec
        if not ak.any(fill_mask):
            return

        px_t = px_lambda_truth[fill_mask]
        py_t = py_lambda_truth[fill_mask]
        pz_t = pz_lambda_truth[fill_mask]

        px_r = px_lambda_recon[fill_mask][:, 0]
        py_r = py_lambda_recon[fill_mask][:, 0]
        pz_r = pz_lambda_recon[fill_mask][:, 0]

        # calculate truth angle
        pt_t = np.hypot(px_t*np.cos(tilt) - pz_t*np.sin(tilt), py_t)
        th_t = np.arctan2(pt_t, pz_t*np.cos(tilt) + px_t*np.sin(tilt))

        # calculate recon angle
        pt_r = np.hypot(px_r*np.cos(tilt) - pz_r*np.sin(tilt), py_r)
        th_r = np.arctan2(pt_r, pz_r*np.cos(tilt) + px_r*np.sin(tilt))

        # store for scatter
        self.truth_theta_list.append(ak.to_numpy(th_t))
        self.recon_theta_list.append(ak.to_numpy(th_r))

        # fill residual histogram
        dtheta_mrad = (th_r - th_t)*1e3
        self.h_dtheta.fill(dtheta=dtheta_mrad)

    def fill_zvertex(self,
                     z_lambda_truth,
                     px_lambda_recon, pz_lambda_recon,
                     x_ref, z_ref,
                     tilt=-0.025):
        """
        Fill z-vertex difference and store scatter arrays.
        Note: z_lambda_truth is in mm, z_ref is also in mm
        """
        valid_truth = ~ak.is_none(z_lambda_truth)
        has_rec = (ak.num(px_lambda_recon, axis=1) > 0)
        fill_mask = valid_truth & has_rec
        if not ak.any(fill_mask):
            return

        z_t = z_lambda_truth[fill_mask]
        x_r = x_ref[fill_mask][:, 0]
        z_r = z_ref[fill_mask][:, 0]

        # Calculate reconstructed z-vertex with tilt correction
        # Note: keeping all values in mm for consistency
        z_vtx_rec = z_r*np.cos(tilt) + x_r*np.sin(tilt)

        # Store arrays in mm for scatter plot
        self.z_truth_list.append(ak.to_numpy(z_t))
        self.z_recon_list.append(ak.to_numpy(z_vtx_rec))

        # Fill residual histogram in [m] for display
        # Convert from mm to m for the histogram
        dz_m = (z_vtx_rec - z_t)/1000.0
        self.h_dzvtx.fill(dzvtx=ak.to_numpy(dz_m))


    def fill_mass(self, mass_r):
        """Fill the mass of the first reconstructed Lambda (if present)."""
        has_lambda = (ak.num(mass_r, axis=1) > 0)
        if not ak.any(has_lambda):
            return
        m_sel = mass_r[has_lambda][:, 0]
        self.h_mass.fill(mass=ak.to_numpy(m_sel))


    def fill_cm_angles_vectorized(
            self,
            pdg_mc,        # [nEvents, nParticles]
            px_mc, py_mc, pz_mc, E_mc,
            pdg_decay,     # [nEvents, nDecay]
            px_decay, py_decay, pz_decay
    ):
        """
        An efficient approach to fill CM neutron angles
        for each event, picking the first Lambda (3122)
        and the first neutron (2112). Then store results
        in self.cm_dphi, self.cm_dtheta.

        This uses Awkward Arrays for the particle selection
        but processes the boost/angle calculations in batches
        to avoid limitations with the vector library.
        """
        # Get the first lambda and neutron indices for each event - this is vectorized
        lambda_mask = (pdg_mc == 3122)
        neutron_mask = (pdg_mc == 2112)
        rec_neutron_mask = (pdg_decay == 2112)

        # Use local_index to get position indices within each event
        idx_local = ak.local_index(pdg_mc, axis=1)
        lambda_indices = idx_local[lambda_mask]
        neutron_indices = idx_local[neutron_mask]

        rec_idx_local = ak.local_index(pdg_decay, axis=1)
        rec_neutron_indices = rec_idx_local[rec_neutron_mask]

        # Get the first lambda and neutron for each event
        first_lambda_idx = ak.firsts(lambda_indices, axis=1)
        first_neutron_idx = ak.firsts(neutron_indices, axis=1)
        first_rec_neutron_idx = ak.firsts(rec_neutron_indices, axis=1)

        # Filter for events with all particles
        valid_events = (~ak.is_none(first_lambda_idx)) & (~ak.is_none(first_neutron_idx)) & (~ak.is_none(first_rec_neutron_idx))
        if not ak.any(valid_events):
            return

        # Convert to numpy arrays for safe processing
        first_lambda_idx_np = ak.to_numpy(first_lambda_idx[valid_events])
        first_neutron_idx_np = ak.to_numpy(first_neutron_idx[valid_events])
        first_rec_neutron_idx_np = ak.to_numpy(first_rec_neutron_idx[valid_events])

        # Get filtered event data
        px_mc_valid = px_mc[valid_events]
        py_mc_valid = py_mc[valid_events]
        pz_mc_valid = pz_mc[valid_events]
        E_mc_valid = E_mc[valid_events]

        px_decay_valid = px_decay[valid_events]
        py_decay_valid = py_decay[valid_events]
        pz_decay_valid = pz_decay[valid_events]

        # Process in NumPy for maximum efficiency
        px_mc_np = ak.to_numpy(ak.flatten(px_mc_valid))
        py_mc_np = ak.to_numpy(ak.flatten(py_mc_valid))
        pz_mc_np = ak.to_numpy(ak.flatten(pz_mc_valid))
        E_mc_np = ak.to_numpy(ak.flatten(E_mc_valid))

        px_decay_np = ak.to_numpy(ak.flatten(px_decay_valid))
        py_decay_np = ak.to_numpy(ak.flatten(py_decay_valid))
        pz_decay_np = ak.to_numpy(ak.flatten(pz_decay_valid))

        # Get cumulative sums of particle counts to map flat indices
        mc_counts = ak.num(px_mc_valid, axis=1)
        decay_counts = ak.num(px_decay_valid, axis=1)

        mc_offsets = np.insert(np.cumsum(ak.to_numpy(mc_counts)), 0, 0)
        decay_offsets = np.insert(np.cumsum(ak.to_numpy(decay_counts)), 0, 0)

        # Process each event using fast NumPy indexing
        lambda_px = []
        lambda_py = []
        lambda_pz = []
        lambda_E = []

        neutron_px = []
        neutron_py = []
        neutron_pz = []
        neutron_E = []

        rec_neutron_px = []
        rec_neutron_py = []
        rec_neutron_pz = []
        rec_neutron_E = []

        # Extract components with fast NumPy indexing
        for i in range(len(first_lambda_idx_np)):
            # Lambda indices in flattened arrays
            lambda_flat_idx = mc_offsets[i] + first_lambda_idx_np[i]
            lambda_px.append(px_mc_np[lambda_flat_idx])
            lambda_py.append(py_mc_np[lambda_flat_idx])
            lambda_pz.append(pz_mc_np[lambda_flat_idx])
            lambda_E.append(E_mc_np[lambda_flat_idx])

            # Neutron indices in flattened arrays
            neutron_flat_idx = mc_offsets[i] + first_neutron_idx_np[i]
            neutron_px.append(px_mc_np[neutron_flat_idx])
            neutron_py.append(py_mc_np[neutron_flat_idx])
            neutron_pz.append(pz_mc_np[neutron_flat_idx])
            neutron_E.append(E_mc_np[neutron_flat_idx])

            # Reconstructed neutron indices
            rec_neutron_flat_idx = decay_offsets[i] + first_rec_neutron_idx_np[i]
            rec_neutron_px.append(px_decay_np[rec_neutron_flat_idx])
            rec_neutron_py.append(py_decay_np[rec_neutron_flat_idx])
            rec_neutron_pz.append(pz_decay_np[rec_neutron_flat_idx])

            # Calculate energy
            mn = 0.9396
            rec_E = np.sqrt(
                px_decay_np[rec_neutron_flat_idx]**2 +
                py_decay_np[rec_neutron_flat_idx]**2 +
                pz_decay_np[rec_neutron_flat_idx]**2 +
                mn**2
            )
            rec_neutron_E.append(rec_E)

        # Make vector calculations in batches (100 at a time)
        batch_size = 100
        for batch_start in range(0, len(lambda_px), batch_size):
            batch_end = min(batch_start + batch_size, len(lambda_px))

            # Process batch
            for i in range(batch_start, batch_end):
                # Create vector objects
                lambda_vec = vector.obj(
                    px=lambda_px[i],
                    py=lambda_py[i],
                    pz=lambda_pz[i],
                    E=lambda_E[i]
                )

                neutron_vec = vector.obj(
                    px=neutron_px[i],
                    py=neutron_py[i],
                    pz=neutron_pz[i],
                    E=neutron_E[i]
                )

                rec_neutron_vec = vector.obj(
                    px=rec_neutron_px[i],
                    py=rec_neutron_py[i],
                    pz=rec_neutron_pz[i],
                    E=rec_neutron_E[i]
                )

                # Boost to CM frame
                n_cm_truth = neutron_vec.boost_p4(lambda_vec)
                n_cm_rec = rec_neutron_vec.boost_p4(lambda_vec)

                # Calculate angles
                phi_truth = n_cm_truth.phi
                theta_truth = n_cm_truth.theta
                phi_rec = n_cm_rec.phi
                theta_rec = n_cm_rec.theta

                # Handle phi wrapping for proper difference calculation
                # This ensures phi differences are always in [-π, π]
                dphi_raw = phi_rec - phi_truth
                dphi_wrapped = np.arctan2(np.sin(dphi_raw), np.cos(dphi_raw))

                # Calculate residuals with correct wrapping
                dphi = dphi_wrapped * np.sin(theta_truth) * 1000.0  # Convert to mrad
                dtheta = (theta_rec - theta_truth) * 1000.0  # Convert to mrad

                # Add to results
                self.cm_dphi.append(float(dphi))
                self.cm_dtheta.append(float(dtheta))


    #######################################
    # 5) Final plotting
    #######################################
    def finalize_and_plot(self, outdir):
        beam_outdir = os.path.join(outdir, self.beam_label)
        os.makedirs(beam_outdir, exist_ok=True)

        # 1) nclusters => single panel
        plt.figure()
        vals = self.h_nclusters.values()
        edges= self.h_nclusters.axes[0].edges
        plt.hist(edges[:-1], bins=edges, weights=vals, histtype='step', color='black')
        plt.yscale('log')
        plt.ylim(1)
        plt.xlabel("Number of ZDC clusters")
        plt.title(f"Beam: {self.beam_label}")

        # Save both PDF and PNG
        base_fname = os.path.join(beam_outdir, f"nclust_{self.beam_label}")
        plt.savefig(f"{base_fname}.pdf")
        plt.savefig(f"{base_fname}.png", dpi=300)
        plt.close()

        # 2) Angles => 2-panel (scatter, residual)
        truth_theta = np.concatenate(self.truth_theta_list) if self.truth_theta_list else np.array([])
        recon_theta = np.concatenate(self.recon_theta_list) if self.recon_theta_list else np.array([])

        fig, axs = plt.subplots(1, 2, figsize=(16, 8))
        fig.suptitle(f"Beam: {self.beam_label} — Reconstructed vs Truth angles")

        # Left: scatter
        plt.sca(axs[0])
        if len(truth_theta) > 0 and len(recon_theta) > 0:
            plt.scatter(truth_theta*1e3, recon_theta*1e3, s=2)
            plt.xlabel(r"$\theta^{*\mathrm{truth}}_{\Lambda}$ [mrad]")
            plt.ylabel(r"$\theta^{*\mathrm{recon}}_{\Lambda}$ [mrad]")
            plt.xlim(0, 3.2)
            plt.ylim(0, 3.2)
        else:
            plt.text(0.5, 0.5, "No angle data", ha='center')

        # Right: residual hist + fit
        plt.sca(axs[1])
        vals = self.h_dtheta.values()
        edges= self.h_dtheta.axes[0].edges
        plt.hist(edges[:-1], bins=edges, weights=vals, color='gray', edgecolor='black')
        bc = 0.5*(edges[:-1]+edges[1:])
        mask = np.abs(bc) < 0.3
        yvals= vals[mask]
        xvals= bc[mask]
        sigma= np.sqrt(yvals) + (yvals==0)
        if np.sum(yvals)>0:
            try:
                popt, pcov = curve_fit(gauss, xvals, yvals, p0=[max(yvals), 0, 0.05], sigma=sigma)
                xx = np.linspace(-1, 1, 200)
                plt.plot(xx, gauss(xx, *popt), 'r-', label=f"$\\sigma={popt[2]:.3f}$ mrad")
                plt.legend()
            except:
                pass
        plt.xlabel(r"$\theta^{*\mathrm{rec}}_{\Lambda} - \theta^{*\mathrm{truth}}_{\Lambda}$ [mrad]")
        plt.ylabel("Counts")

        plt.tight_layout()
        base_fname = os.path.join(beam_outdir, f"thetastar_{self.beam_label}")
        plt.savefig(f"{base_fname}.pdf")
        plt.savefig(f"{base_fname}.png", dpi=300)
        plt.close()

        # 3) z-vtx => 2-panel (scatter, residual)
        z_truth_flat = np.concatenate(self.z_truth_list) if self.z_truth_list else np.array([])
        z_recon_flat = (np.concatenate(self.z_recon_list) if self.z_recon_list else np.array([]))/1000.0

        fig, axs = plt.subplots(1, 2, figsize=(16, 8))
        fig.suptitle(f"Beam: {self.beam_label} — z-vtx reconstruction")

        # Left: scatter truth vs recon (values in mm)
        plt.sca(axs[0])
        if len(z_truth_flat) > 0 and len(z_recon_flat) > 0:
            # Check the range to set appropriate limits
            z_min = min(np.min(z_truth_flat), np.min(z_recon_flat))
            z_max = max(np.max(z_truth_flat), np.max(z_recon_flat))

            # Plot without dividing by 1000 (keep in mm)
            plt.scatter(z_truth_flat, z_recon_flat, s=2)
            plt.xlabel(r"$z^{\mathrm{truth}}_{vtx}$ [mm]")
            plt.ylabel(r"$z^{\mathrm{recon}}_{vtx}$ [mm]")

            # Set range with padding
            pad = (z_max - z_min) * 0.1
            plt.xlim(z_min - pad, z_max + pad)
            plt.ylim(z_min - pad, z_max + pad)

            # Add diagonal line
            plt.plot([z_min - pad, z_max + pad], [z_min - pad, z_max + pad], 'r--', alpha=0.7)
        else:
            plt.text(0.5, 0.5, "No z-vertex data", ha='center')

        # Right: residual histogram + fit (values already in m)
        plt.sca(axs[1])
        vals = self.h_dzvtx.values()
        edges= self.h_dzvtx.axes[0].edges
        plt.hist(edges[:-1], bins=edges, weights=vals, color='gray', edgecolor='black')
        bc = 0.5*(edges[:-1]+edges[1:])

        # Determine appropriate range for the fit based on the data
        if len(vals) > 0:
            non_zero_vals = vals[vals > 0]
            if len(non_zero_vals) > 0:
                mask = np.abs(bc) < np.percentile(np.abs(bc[vals > 0]), 95)
                yvals = vals[mask]
                xvals = bc[mask]
                sigma = np.sqrt(yvals) + (yvals==0)

                if np.sum(yvals) > 0:
                    try:
                        popt, pcov = curve_fit(gauss, xvals, yvals, p0=[max(yvals), 0, 1], sigma=sigma)
                        xx = np.linspace(min(xvals), max(xvals), 200)
                        plt.plot(xx, gauss(xx, *popt), 'r-', label=f"$\\sigma={popt[2]:.2f}$ m")
                        plt.legend()
                    except:
                        pass

        plt.xlabel(r"$(z^{\mathrm{rec}}_{vtx}-z^{\mathrm{truth}}_{vtx})$ [m]")
        plt.ylabel("Counts")

        plt.tight_layout()
        base_fname = os.path.join(beam_outdir, f"zvtx_{self.beam_label}")
        plt.savefig(f"{base_fname}.pdf")
        plt.savefig(f"{base_fname}.png", dpi=300)
        plt.close()

        # 4) Mass => single-panel
        plt.figure()
        vals = self.h_mass.values()
        edges= self.h_mass.axes[0].edges
        plt.hist(edges[:-1], bins=edges, weights=vals, color='gray', edgecolor='black')
        pdg_mass = 1.115683
        plt.axvline(pdg_mass, ls='--', color='g', lw=2)
        bc = 0.5*(edges[:-1]+edges[1:])
        mask = (bc>pdg_mass-0.05) & (bc<pdg_mass+0.05)
        yvals= vals[mask]
        xvals= bc[mask]
        sigma= np.sqrt(yvals)+(yvals==0)
        if np.sum(yvals)>0:
            try:
                popt, pcov = curve_fit(gauss, xvals, yvals,
                                       p0=[max(yvals), pdg_mass, 0.04],
                                       sigma=sigma)
                xx = np.linspace(0.8,1.3,200)
                plt.plot(xx, gauss(xx, *popt), 'r-',
                         label=f"$\\sigma={popt[2]:.3f}$ GeV")
                plt.legend()
            except:
                pass
        plt.xlabel(r"$m_{\Lambda}^{\mathrm{rec}}$ [GeV]")
        plt.ylabel("Counts")
        plt.title(f"Beam: {self.beam_label} — Reconstructed Lambda mass")
        base_fname = os.path.join(beam_outdir, f"lambda_mass_{self.beam_label}")
        plt.savefig(f"{base_fname}.pdf")
        plt.savefig(f"{base_fname}.png", dpi=300)
        plt.close()

        # 5) CM neutron angles => if any were filled
        if len(self.cm_dphi) > 0 or len(self.cm_dtheta) > 0:
            # Calculate appropriate range for phi plot
            phi_array = np.array(self.cm_dphi)
            phi_range = np.percentile(np.abs(phi_array), 99)
            phi_range = min(max(phi_range * 1.2, 50), 300)  # Reasonable range

            plt.figure()
            plt.hist(self.cm_dphi, bins=100, range=(-phi_range, phi_range))
            plt.xlabel(r"$(\phi^n_{\mathrm{rec}}-\phi^n_{\mathrm{truth}})\times \sin(\theta^n_{\mathrm{truth}})$ [mrad]")
            base_fname = os.path.join(beam_outdir, f"neutron_phi_cm_res_{self.beam_label}")
            plt.savefig(f"{base_fname}.pdf")
            plt.savefig(f"{base_fname}.png", dpi=300)
            plt.close()

            # Calculate appropriate range for theta plot
            theta_array = np.array(self.cm_dtheta)
            theta_range = np.percentile(np.abs(theta_array), 99)
            theta_range = min(max(theta_range * 1.2, 20), 100)  # Reasonable range

            plt.figure()
            plt.hist(self.cm_dtheta, bins=100, range=(-theta_range, theta_range))
            plt.xlabel(r"$\theta^n_{\mathrm{rec}}-\theta^n_{\mathrm{truth}}$ [mrad]")
            base_fname = os.path.join(beam_outdir, f"neutron_theta_cm_res_{self.beam_label}")
            plt.savefig(f"{base_fname}.pdf")
            plt.savefig(f"{base_fname}.png", dpi=300)
            plt.close()


############################
# 4) Main
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

    # Minimal set of branches we need
    needed_branches = [
        "HcalFarForwardZDCClusters.position.x",
        "MCParticles.PDG",
        "MCParticles.momentum.x",
        "MCParticles.momentum.y",
        "MCParticles.momentum.z",
        "MCParticles.vertex.z",
        "MCParticles.mass",
        "ReconstructedFarForwardZDCLambdas.momentum.x",
        "ReconstructedFarForwardZDCLambdas.momentum.y",
        "ReconstructedFarForwardZDCLambdas.momentum.z",
        "ReconstructedFarForwardZDCLambdas.mass",
        "ReconstructedFarForwardZDCLambdas.referencePoint.x",
        "ReconstructedFarForwardZDCLambdas.referencePoint.z"
    ]
    if not args.no_cm_angles:
        needed_branches += [
            "ReconstructedFarForwardZDCLambdaDecayProductsCM.PDG",
            "ReconstructedFarForwardZDCLambdaDecayProductsCM.momentum.x",
            "ReconstructedFarForwardZDCLambdaDecayProductsCM.momentum.y",
            "ReconstructedFarForwardZDCLambdaDecayProductsCM.momentum.z",
        ]

    beam_analyses = {}

    for beam_label, filelist in beams_map.items():
        print(f"\n=== Beam: {beam_label} ===")
        print("   Files:", filelist)
        beam_analyses[beam_label] = LambdaAnalysis(beam_label)
        analysis_obj = beam_analyses[beam_label]

        file_dict = {fname: args.tree for fname in filelist}

        chunk_counter = 0
        for chunk in uproot.iterate(
                file_dict,
                expressions=needed_branches,
                step_size=args.step_size,
                library="ak"
        ):
            chunk_counter += 1
            if chunk_counter == 1:
                print(f"[Debug] Fields in first chunk for beam={beam_label}:")
                for k in chunk.fields:
                    print("   ", k)
            print(f"[Debug] #Events in this chunk = {len(chunk)}")

            # 1) fill cluster count
            if "HcalFarForwardZDCClusters.position.x" in chunk.fields:
                cluster_x = chunk["HcalFarForwardZDCClusters.position.x"]
                analysis_obj.fill_cluster_count(cluster_x)

            # 2) truth vs recon angles by PDG=3122
            have_angle_fields = {
                "MCParticles.PDG",
                "MCParticles.momentum.x",
                "MCParticles.momentum.y",
                "MCParticles.momentum.z",
                "ReconstructedFarForwardZDCLambdas.momentum.x",
                "ReconstructedFarForwardZDCLambdas.momentum.y",
                "ReconstructedFarForwardZDCLambdas.momentum.z"
            }.issubset(chunk.fields)
            if have_angle_fields:
                pdg_mc = chunk["MCParticles.PDG"]
                px_mc  = chunk["MCParticles.momentum.x"]
                py_mc  = chunk["MCParticles.momentum.y"]
                pz_mc  = chunk["MCParticles.momentum.z"]

                # PDG=3122 => truth-lambda
                is_lambda = (pdg_mc == 3122)
                px_lambda = px_mc[is_lambda]
                py_lambda = py_mc[is_lambda]
                pz_lambda = pz_mc[is_lambda]
                has_truth_lambda = (ak.num(px_lambda, axis=1) > 0)

                px_lambda_first = ak.mask(px_lambda[:,0], has_truth_lambda)
                py_lambda_first = ak.mask(py_lambda[:,0], has_truth_lambda)
                pz_lambda_first = ak.mask(pz_lambda[:,0], has_truth_lambda)

                # Reconstructed
                px_rec = chunk["ReconstructedFarForwardZDCLambdas.momentum.x"]
                py_rec = chunk["ReconstructedFarForwardZDCLambdas.momentum.y"]
                pz_rec = chunk["ReconstructedFarForwardZDCLambdas.momentum.z"]

                analysis_obj.fill_truth_recon_angles(px_lambda_first,
                                                     py_lambda_first,
                                                     pz_lambda_first,
                                                     px_rec, py_rec, pz_rec)

            # 3) z-vtx
            have_zvtx_fields = {
                "MCParticles.PDG",
                "MCParticles.vertex.z",
                "ReconstructedFarForwardZDCLambdas.referencePoint.x",
                "ReconstructedFarForwardZDCLambdas.referencePoint.z",
                "ReconstructedFarForwardZDCLambdas.momentum.x"
            }.issubset(chunk.fields)
            if have_zvtx_fields:
                pdg_mc = chunk["MCParticles.PDG"]
                z_mc   = chunk["MCParticles.vertex.z"]

                is_lambda = (pdg_mc == 3122)
                z_lambda  = z_mc[is_lambda]
                has_lam   = (ak.num(z_lambda, axis=1) > 0)
                z_lambda_first = ak.mask(z_lambda[:,0], has_lam)

                px_rec = chunk["ReconstructedFarForwardZDCLambdas.momentum.x"]
                pz_rec = chunk["ReconstructedFarForwardZDCLambdas.momentum.z"]
                x_ref  = chunk["ReconstructedFarForwardZDCLambdas.referencePoint.x"]
                z_ref  = chunk["ReconstructedFarForwardZDCLambdas.referencePoint.z"]

                analysis_obj.fill_zvertex(z_lambda_first, px_rec, pz_rec, x_ref, z_ref)

            # 4) mass
            if "ReconstructedFarForwardZDCLambdas.mass" in chunk.fields:
                mass_r = chunk["ReconstructedFarForwardZDCLambdas.mass"]
                analysis_obj.fill_mass(mass_r)

            # 5) CM angles (vectorized)
            if not args.no_cm_angles:
                needed_cms = {
                    "ReconstructedFarForwardZDCLambdaDecayProductsCM.PDG",
                    "ReconstructedFarForwardZDCLambdaDecayProductsCM.momentum.x",
                    "ReconstructedFarForwardZDCLambdaDecayProductsCM.momentum.y",
                    "ReconstructedFarForwardZDCLambdaDecayProductsCM.momentum.z",
                    "MCParticles.PDG",
                    "MCParticles.momentum.x",
                    "MCParticles.momentum.y",
                    "MCParticles.momentum.z",
                    "MCParticles.mass"
                }
                if needed_cms.issubset(chunk.fields):
                    pdg_mc  = chunk["MCParticles.PDG"]
                    px_mc   = chunk["MCParticles.momentum.x"]
                    py_mc   = chunk["MCParticles.momentum.y"]
                    pz_mc   = chunk["MCParticles.momentum.z"]
                    mass_mc = chunk["MCParticles.mass"]
                    E_mc    = np.sqrt(mass_mc**2 + px_mc**2 + py_mc**2 + pz_mc**2)

                    pdg_cm = chunk["ReconstructedFarForwardZDCLambdaDecayProductsCM.PDG"]
                    px_cm  = chunk["ReconstructedFarForwardZDCLambdaDecayProductsCM.momentum.x"]
                    py_cm  = chunk["ReconstructedFarForwardZDCLambdaDecayProductsCM.momentum.y"]
                    pz_cm  = chunk["ReconstructedFarForwardZDCLambdaDecayProductsCM.momentum.z"]

                    # Vectorized fill => appends to self.cm_dphi, self.cm_dtheta
                    analysis_obj.fill_cm_angles_vectorized(
                        pdg_mc, px_mc, py_mc, pz_mc, E_mc,
                        pdg_cm, px_cm, py_cm, pz_cm
                    )

    # Produce final plots
    for beam_label, analysis_obj in beam_analyses.items():
        analysis_obj.finalize_and_plot(outdir)


if __name__ == "__main__":
    main()
