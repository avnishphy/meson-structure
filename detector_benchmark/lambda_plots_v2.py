#!/usr/bin/env python3

import argparse
import numpy as np
import awkward as ak
import uproot
import hist
import matplotlib.pyplot as plt
import mplhep as hep
import re, os
from   scipy.optimize import curve_fit

hep.style.use("CMS")
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['savefig.facecolor'] = 'white'
plt.rcParams['savefig.bbox']     = 'tight'

############################
# 1) Command-line arguments
############################
def parse_args():
    parser = argparse.ArgumentParser(
        description="Lambda recon analysis, chunked uproot + Scikit-HEP hist, using PDG-based selection."
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
    return A * np.exp(-0.5 * ((x - mu)/sigma)**2)

def extract_beam_label(filename):
    """
    Heuristic to parse e.g. '5x41', '10x100', '18x275' from the filename.
    Adjust to your naming convention as needed.
    """
    match = re.search(r'_(\d+x\d+)_', os.path.basename(filename))
    if match:
        return match.group(1)
    return "unknownBeam"

############################
# 3) Analysis class for each beam
############################
class LambdaAnalysis:
    def __init__(self, beam_label):
        self.beam_label = beam_label
        # 1D hist objects
        self.h_nclusters = hist.Hist(
            hist.axis.Regular(20, 0, 20, name="nclust", label="ZDC cluster count")
        )
        self.h_dtheta = hist.Hist(
            hist.axis.Regular(50, -1, 1, name="dtheta",
                              label=r"$\theta_{rec}-\theta_{truth}$ [mrad]")
        )
        self.h_dzvtx = hist.Hist(
            hist.axis.Regular(50, -10, 10, name="dzvtx",
                              label=r"$z^{rec}_{vtx}-z^{truth}_{vtx}$ [m]")
        )
        self.h_mass = hist.Hist(
            hist.axis.Regular(100, 1.0, 1.25, name="mass", label=r"$m_{\Lambda}^{rec}$ [GeV]")
        )

        # Arrays for scatter plots
        self.truth_theta_list = []
        self.recon_theta_list = []

        # Arrays for CM angle resolution
        self.cm_dphi   = []
        self.cm_dtheta = []

    #######################################
    # Fill methods for each type of plot
    #######################################
    def fill_cluster_count(self, cluster_x):
        nclust = ak.num(cluster_x, axis=1)
        self.h_nclusters.fill(nclust=nclust)

    def fill_truth_recon_angles(self,
                                px_lambda_truth, py_lambda_truth, pz_lambda_truth,
                                px_lambda_recon, py_lambda_recon, pz_lambda_recon,
                                tilt=-0.025):
        """
        We now have the truth-level Lambda from PDG=3122, plus the first
        reconstructed Lambda. Compute the angle difference in the same
        way the old code did.

        - px_lambda_truth, py_lambda_truth, pz_lambda_truth: (Nevents) => first Lambda per event
        - px_lambda_recon, py_lambda_recon, pz_lambda_recon: jagged => we pick the first recon if multiple
        """
        # We'll only fill events that have at least one recon-lambda.
        has_rec = (ak.num(px_lambda_recon, axis=1) > 0)
        # Corresponding truth must also exist => if px_lambda_truth is a float array,
        # we rely on shape matching. Or we specifically do `has_truth = px_lambda_truth != None`.
        # But here, if an event has no truth-lambda, we may have a masked value or empty.
        # Let's check that the lengths match. We'll assume we only call this if a truth-lambda was found.

        # flatten truth arrays
        # We have shape (Nevents, ) for px_lambda_truth, etc. They might be None for events with no lambda.
        # Let's create a mask that identifies events with valid truth-lambda
        valid_truth = ~ak.is_none(px_lambda_truth)  # or ~ak.is_none
        # combined mask
        fill_mask = has_rec & valid_truth

        # select only events that have both rec-lambda and truth-lambda
        px_t = px_lambda_truth[fill_mask]
        py_t = py_lambda_truth[fill_mask]
        pz_t = pz_lambda_truth[fill_mask]

        # recon => pick only the first
        px_r_first = px_lambda_recon[fill_mask][:, 0]
        py_r_first = py_lambda_recon[fill_mask][:, 0]
        pz_r_first = pz_lambda_recon[fill_mask][:, 0]

        # angles
        pt_truth = np.hypot(px_t*np.cos(tilt) - pz_t*np.sin(tilt), py_t)
        th_truth = np.arctan2(pt_truth, pz_t*np.cos(tilt) + px_t*np.sin(tilt))

        pt_recon = np.hypot(px_r_first*np.cos(tilt) - pz_r_first*np.sin(tilt), py_r_first)
        th_recon = np.arctan2(pt_recon, pz_r_first*np.cos(tilt) + px_r_first*np.sin(tilt))

        # store in lists for final scatter
        self.truth_theta_list.append(ak.to_numpy(th_truth))
        self.recon_theta_list.append(ak.to_numpy(th_recon))

        # fill dtheta histogram (in mrad)
        dtheta_mrad = (th_recon - th_truth)*1e3
        self.h_dtheta.fill(dtheta=dtheta_mrad)

    def fill_zvertex(self,
                     z_lambda_truth,
                     px_lambda_recon, pz_lambda_recon,
                     x_ref, z_ref,
                     tilt=-0.025):
        """
        Fill the difference in z-vertex: (recon - truth).
        Now that we have z_lambda_truth from PDG=3122, we only fill for events
        that have a recon-lambda as well.
        """
        has_rec = (ak.num(px_lambda_recon, axis=1) > 0)
        valid_truth = ~ak.is_none(z_lambda_truth)
        fill_mask = has_rec & valid_truth

        z_truth_sel = z_lambda_truth[fill_mask]

        x_ref_sel = x_ref[fill_mask][:, 0]
        z_ref_sel = z_ref[fill_mask][:, 0]
        z_vtx_rec = z_ref_sel*np.cos(tilt) + x_ref_sel*np.sin(tilt)

        dz_m = (z_vtx_rec - z_truth_sel)/1000.0  # mm -> m, if that's your unit
        self.h_dzvtx.fill(dzvtx=ak.to_numpy(dz_m))

    def fill_mass(self, mass_r):
        """
        Fill the mass distribution for the first Reconstructed Lambda
        (if multiple exist in an event).
        """
        has_lambda = (ak.num(mass_r, axis=1) > 0)
        m_sel = mass_r[has_lambda][:, 0]
        self.h_mass.fill(mass=ak.to_numpy(m_sel))

    def fill_cm_angles(self,
                       pdg_mc, px_mc, py_mc, pz_mc, E_mc,
                       pdg_decay, px_decay, py_decay, pz_decay):
        """
        If we also want to identify the Lambda at truth level by PDG=3122
        and the neutron by PDG=2112, we can do that *inside* this method
        or outside. For now, let's replicate the old code's partial approach:
        - We find the event's *first* Lambda (PDG=3122).
        - Then we find the event's first neutron (PDG=2112).
        - Boost the truth-neutron into the Lambda rest frame => get angles.
        - Then find the reco neutron in ReconstructedFarForwardZDCLambdaDecayProductsCM (PDG=2112)
          and do the same boost, measure the difference.

        The old code used fixed indices [2],[3]. Now we do PDG-based.
        """
        import ROOT

        for i_evt in range(len(pdg_mc)):
            # find Lambdas
            lambda_mask_evt = (pdg_mc[i_evt] == 3122)
            n_lambda = np.count_nonzero(lambda_mask_evt)
            if n_lambda < 1:
                continue
            # pick the first
            lambda_idx = np.nonzero(lambda_mask_evt)[0][0]
            # build lambda TLorentzVector
            px_l = px_mc[i_evt][lambda_idx]
            py_l = py_mc[i_evt][lambda_idx]
            pz_l = pz_mc[i_evt][lambda_idx]
            E_l  = E_mc[i_evt][lambda_idx]

            l_vec = ROOT.TLorentzVector(px_l, py_l, pz_l, E_l)

            # find neutrons
            n_mask_evt = (pdg_mc[i_evt] == 2112)
            if np.count_nonzero(n_mask_evt) < 1:
                continue
            n_idx = np.nonzero(n_mask_evt)[0][0]
            px_n  = px_mc[i_evt][n_idx]
            py_n  = py_mc[i_evt][n_idx]
            pz_n  = pz_mc[i_evt][n_idx]
            E_n   = E_mc[i_evt][n_idx]

            n_vec_truth = ROOT.TLorentzVector(px_n, py_n, pz_n, E_n)

            # Boost truth neutron to Lambda rest frame
            boost = -l_vec.BoostVector()
            n_cm_truth = n_vec_truth.Clone()
            n_cm_truth.Boost(boost)
            phi_truth  = n_cm_truth.Phi()
            theta_truth= n_cm_truth.Theta()

            # Now find the *reconstructed* neutron in the decay products
            if len(pdg_decay[i_evt]) == 0:
                continue
            # pick the first PDG=2112
            neutron_matches = np.where(pdg_decay[i_evt] == 2112)[0]
            if len(neutron_matches) == 0:
                continue
            n_idx_rec = neutron_matches[0]
            pxn = px_decay[i_evt][n_idx_rec]
            pyn = py_decay[i_evt][n_idx_rec]
            pzn = pz_decay[i_evt][n_idx_rec]
            mn  = 0.9396  # approximate neutron mass
            E_rec_n = np.sqrt(pxn**2 + pyn**2 + pzn**2 + mn**2)
            n_vec_rec = ROOT.TLorentzVector(pxn, pyn, pzn, E_rec_n)
            n_vec_rec.Boost(boost)
            phi_rec   = n_vec_rec.Phi()
            theta_rec = n_vec_rec.Theta()

            dphi   = (phi_rec - phi_truth)*np.sin(theta_truth)*1000.0
            dtheta = (theta_rec - theta_truth)*1000.0
            self.cm_dphi.append(dphi)
            self.cm_dtheta.append(dtheta)

    #######################################
    # Final plotting
    #######################################
    def finalize_and_plot(self, outdir):
        beam_outdir = os.path.join(outdir, self.beam_label)
        os.makedirs(beam_outdir, exist_ok=True)

        # 1) nclusters
        plt.figure()
        vals = self.h_nclusters.values()
        edges= self.h_nclusters.axes[0].edges
        plt.hist(edges[:-1], bins=edges, weights=vals, histtype='step', color='black')
        plt.yscale('log')
        plt.ylim(1)
        plt.xlabel("Number of ZDC clusters")
        plt.title(f"Beam: {self.beam_label}")
        fname = os.path.join(beam_outdir, f"nclust_{self.beam_label}.pdf")
        plt.savefig(fname)
        plt.close()

        # 2) \theta^* scatter + residual
        truth_theta = np.concatenate(self.truth_theta_list) if self.truth_theta_list else np.array([])
        recon_theta = np.concatenate(self.recon_theta_list) if self.recon_theta_list else np.array([])

        fig, axs = plt.subplots(1, 3, figsize=(24, 8))
        fig.suptitle(f"Beam: {self.beam_label} — Reconstructed vs Truth angles")

        # Panel 1: scatter
        plt.sca(axs[0])
        if len(truth_theta) > 0 and len(recon_theta) > 0:
            plt.scatter(truth_theta*1e3, recon_theta*1e3, s=2)
            plt.xlabel(r"$\theta^{*\mathrm{truth}}_{\Lambda}$ [mrad]")
            plt.ylabel(r"$\theta^{*\mathrm{recon}}_{\Lambda}$ [mrad]")
            plt.xlim(0, 3.2)
            plt.ylim(0, 3.2)
        else:
            plt.text(0.5, 0.5, "No data for scatter", ha='center')

        # Panel 2: dtheta hist + fit
        plt.sca(axs[1])
        vals = self.h_dtheta.values()
        edges= self.h_dtheta.axes[0].edges
        plt.hist(edges[:-1], bins=edges, weights=vals, color='gray', edgecolor='black')
        bc = 0.5*(edges[:-1] + edges[1:])
        mask = np.abs(bc) < 0.3
        yvals= vals[mask]
        xvals= bc[mask]
        sigma= np.sqrt(yvals) + (yvals==0)
        if np.sum(yvals)>0:
            try:
                popt, pcov = curve_fit(gauss, xvals, yvals, p0=[max(yvals), 0, 0.05], sigma=sigma)
                xx = np.linspace(-1, 1, 200)
                plt.plot(xx, gauss(xx, *popt), 'r-', label=f"sigma={popt[2]:.3f} mrad")
                plt.legend()
            except:
                pass
        plt.xlabel(r"$\theta^{*\mathrm{rec}}_{\Lambda}-\theta^{*\mathrm{truth}}_{\Lambda}$ [mrad]")
        plt.ylabel("Counts")

        # Panel 3: single beam
        plt.sca(axs[2])
        plt.text(0.1, 0.5,
                 "Single beam => no momentum loop.\n"
                 "Check combined if desired.",
                 ha='left')
        plt.axis('off')

        plt.tight_layout()
        outfn = os.path.join(beam_outdir, f"thetastar_{self.beam_label}.pdf")
        plt.savefig(outfn)
        plt.close()

        # 3) z-vtx
        fig, axs = plt.subplots(1, 3, figsize=(24,8))
        fig.suptitle(f"Beam: {self.beam_label} — z-vtx resolution")

        # Panel 1: no scatter stored chunkwise
        plt.sca(axs[0])
        plt.text(0.5, 0.5, "No chunkwise scatter stored", ha='center')
        plt.axis('off')

        # Panel 2: dzvtx hist + fit
        plt.sca(axs[1])
        vals = self.h_dzvtx.values()
        edges= self.h_dzvtx.axes[0].edges
        plt.hist(edges[:-1], bins=edges, weights=vals, color='gray', edgecolor='black')
        bc   = 0.5*(edges[:-1]+edges[1:])
        mask = np.abs(bc) < 5
        yvals= vals[mask]
        xvals= bc[mask]
        sigma= np.sqrt(yvals) + (yvals==0)
        if np.sum(yvals)>0:
            try:
                popt, pcov = curve_fit(gauss, xvals, yvals, p0=[max(yvals), 0, 1], sigma=sigma)
                xx = np.linspace(-5,5,200)
                plt.plot(xx, gauss(xx, *popt), 'r-', label=f"sigma={popt[2]:.2f} m")
                plt.legend()
            except:
                pass
        plt.xlabel(r"$z^{rec}_{vtx}-z^{truth}_{vtx}$ [m]")
        plt.ylabel("Counts")

        # Panel 3
        plt.sca(axs[2])
        plt.text(0.1, 0.5,
                 "Single beam => no momentum loop.",
                 ha='left')
        plt.axis('off')
        plt.tight_layout()
        outfn = os.path.join(beam_outdir, f"zvtx_{self.beam_label}.pdf")
        plt.savefig(outfn)
        plt.close()

        # 4) Lambda mass
        fig, axs = plt.subplots(1,2, figsize=(16,8))
        fig.suptitle(f"Beam: {self.beam_label} — Reconstructed Lambda mass")

        plt.sca(axs[0])
        vals = self.h_mass.values()
        edges= self.h_mass.axes[0].edges
        plt.hist(edges[:-1], bins=edges, weights=vals, color='gray', edgecolor='black')
        pdg_mass = 1.115683
        plt.axvline(pdg_mass, ls='--', color='g', lw=2)
        bc   = 0.5*(edges[:-1]+edges[1:])
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
                         label=f"sigma={popt[2]:.3f} GeV")
                plt.legend()
            except:
                pass
        plt.xlabel(r"$m_{\Lambda}^{rec}$ [GeV]")
        plt.ylabel("Counts")

        plt.sca(axs[1])
        plt.text(0.1, 0.5,
                 "Single beam => no momentum loop.",
                 ha='left')
        plt.axis('off')
        plt.tight_layout()
        outfn = os.path.join(beam_outdir, f"lambda_mass_{self.beam_label}.pdf")
        plt.savefig(outfn)
        plt.close()

        # 5) CM neutron angles
        if len(self.cm_dphi) > 0 or len(self.cm_dtheta) > 0:
            plt.figure()
            plt.hist(self.cm_dphi, bins=100, range=(-300,300))
            plt.xlabel(r"$(\phi^n_{cm,rec}-\phi^n_{cm,truth})\times\sin(\theta^n_{cm,truth})$ [mrad]")
            outfn = os.path.join(beam_outdir, f"neutron_phi_cm_res_{self.beam_label}.pdf")
            plt.savefig(outfn)
            plt.close()

            plt.figure()
            plt.hist(self.cm_dtheta, bins=100, range=(-1000,1000))
            plt.xlabel(r"$\theta^n_{cm,rec}-\theta^n_{cm,truth}$ [mrad]")
            outfn = os.path.join(beam_outdir, f"neutron_theta_cm_res_{self.beam_label}.pdf")
            plt.savefig(outfn)
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

    # We only read the branches we need.
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
        "ReconstructedFarForwardZDCLambdas.referencePoint.z",
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
                library="ak",
        ):
            chunk_counter += 1
            if chunk_counter == 1:
                print(f"[Debug] Fields in first chunk for beam={beam_label}:")
                for k in chunk.fields:
                    print("   ", k)
                print(f"[Debug] #Events in this chunk = {len(chunk)}")

            fields = chunk.fields

            # 1) fill cluster count if present
            if "HcalFarForwardZDCClusters.position.x" in fields:
                cluster_x = chunk["HcalFarForwardZDCClusters.position.x"]
                analysis_obj.fill_cluster_count(cluster_x)

            # 2) truth vs recon angles by PDG=3122
            #   find the *first* truth-lambda (if any) in each event
            if {"MCParticles.PDG",
                "MCParticles.momentum.x",
                "MCParticles.momentum.y",
                "MCParticles.momentum.z",
                "ReconstructedFarForwardZDCLambdas.momentum.x",
                "ReconstructedFarForwardZDCLambdas.momentum.y",
                "ReconstructedFarForwardZDCLambdas.momentum.z"}.issubset(fields):
                pdg_mc = chunk["MCParticles.PDG"]
                px_mc  = chunk["MCParticles.momentum.x"]
                py_mc  = chunk["MCParticles.momentum.y"]
                pz_mc  = chunk["MCParticles.momentum.z"]

                # select PDG=3122
                is_lambda = (pdg_mc == 3122)
                # px_lambda, py_lambda, pz_lambda => jagged arrays of shape (Nevents, #lambdas)
                px_lambda = px_mc[is_lambda]
                py_lambda = py_mc[is_lambda]
                pz_lambda = pz_mc[is_lambda]

                # if an event has >=1 lambda, pick the first
                has_truth_lambda = (ak.num(px_lambda, axis=1)>0)

                # For events with no lambda, we'll store None or masked
                # Let's do an option: fill with None if no lambda
                # e.g.:
                px_lambda_first = ak.where(has_truth_lambda,
                                           px_lambda[:,0],
                                           None)
                py_lambda_first = ak.where(has_truth_lambda,
                                           py_lambda[:,0],
                                           None)
                pz_lambda_first = ak.where(has_truth_lambda,
                                           pz_lambda[:,0],
                                           None)

                # Reconstructed
                px_rec = chunk["ReconstructedFarForwardZDCLambdas.momentum.x"]
                py_rec = chunk["ReconstructedFarForwardZDCLambdas.momentum.y"]
                pz_rec = chunk["ReconstructedFarForwardZDCLambdas.momentum.z"]

                analysis_obj.fill_truth_recon_angles(px_lambda_first,
                                                     py_lambda_first,
                                                     pz_lambda_first,
                                                     px_rec, py_rec, pz_rec)

            # 3) z-vtx difference => also use the same "first lambda" approach
            if {"MCParticles.PDG",
                "MCParticles.vertex.z",
                "ReconstructedFarForwardZDCLambdas.referencePoint.x",
                "ReconstructedFarForwardZDCLambdas.referencePoint.z",
                "ReconstructedFarForwardZDCLambdas.momentum.x"}.issubset(fields):
                pdg_mc = chunk["MCParticles.PDG"]
                z_mc   = chunk["MCParticles.vertex.z"]

                is_lambda = (pdg_mc == 3122)
                z_lambda  = z_mc[is_lambda]
                has_lam   = (ak.num(z_lambda, axis=1)>0)
                z_lambda_first = ak.where(has_lam, z_lambda[:,0], None)

                px_rec = chunk["ReconstructedFarForwardZDCLambdas.momentum.x"]
                pz_rec = chunk["ReconstructedFarForwardZDCLambdas.momentum.z"]
                x_ref  = chunk["ReconstructedFarForwardZDCLambdas.referencePoint.x"]
                z_ref  = chunk["ReconstructedFarForwardZDCLambdas.referencePoint.z"]

                analysis_obj.fill_zvertex(z_lambda_first, px_rec, pz_rec, x_ref, z_ref)

            # 4) fill mass
            if "ReconstructedFarForwardZDCLambdas.mass" in fields:
                mass_r = chunk["ReconstructedFarForwardZDCLambdas.mass"]
                analysis_obj.fill_mass(mass_r)

            # 5) CM angles
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
                if needed_cms.issubset(fields):
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

                    analysis_obj.fill_cm_angles(
                        pdg_mc, px_mc, py_mc, pz_mc, E_mc,
                        pdg_cm, px_cm, py_cm, pz_cm
                    )

    # Final plots
    for beam_label, analysis_obj in beam_analyses.items():
        analysis_obj.finalize_and_plot(outdir)

if __name__ == "__main__":
    main()
