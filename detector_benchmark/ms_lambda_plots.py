#!/usr/bin/env python3
"""
python ms_lambda_plots.py \
    -i data/k_lambda_5x41_5000evt_001.edm4eic.root \
       data/k_lambda_10x100_5000evt_001.edm4eic.root \
       data/k_lambda_18x275_5000evt_001.edm4eic.root \
    -o my_plots
"""

import argparse
import numpy as np
import awkward as ak
import uproot
import matplotlib.pyplot as plt
import mplhep as hep
import sys, os, re

# Use a CMS-like style
hep.style.use("CMS")
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['savefig.facecolor'] = 'white'
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams["figure.figsize"] = (7, 7)


#######################
# 1) Parse arguments
#######################
def parse_args():
    parser = argparse.ArgumentParser(description="Refactored Lambda reconstruction analysis script, reproducing the same plots as the old code.")
    parser.add_argument("-i", "--input-files", nargs="+", required=True, help="List of .root files to analyze (e.g. k_lambda_5x41_XXXX.edm4eic.root, etc.)")
    parser.add_argument("-o", "--outdir", default="output_plots", help="Directory where plots will be saved.")
    parser.add_argument("--tree", default="events", help="Name of the TTree/TTrees to read (default: 'events').")
    parser.add_argument("--no-combine", action="store_true", help="If set, skip any combined plots over all beam energies.")
    return parser.parse_args()


#######################
# 2) Utility functions
#######################
def gauss(x, A, mu, sigma):
    """Simple Gaussian for curve_fit."""
    return A * np.exp(-0.5 * ((x - mu)/sigma)**2)


def extract_beam_energy(filename):
    """
    Attempt to parse a beam-energy label (e.g. '5x41', '10x100', '18x275')
    from the filename. Adjust the regex/pattern as needed.
    """
    # Example filenames: k_lambda_5x41_5000evt_200.edm4eic.root
    #                   k_lambda_10x100_5000evt_001.edm4eic.root
    # We look for a pattern like: '_5x41_' or '_10x100_' or '_18x275_' ...
    match = re.search(r'_(\d+x\d+)_', filename)
    if match:
        return match.group(1)  # e.g. "5x41"
    return "unknownBeam"


#######################
# 3) Beam-by-beam analysis
#######################
def analyze_beam(beam_label, events, outdir):
    """
    Perform the same histograms, fits, and plots
    that the old script did for one beam's dataset.
    Saves beam-specific plots, e.g. ncluster_5x41.pdf, etc.
    """

    print(f"\n=== Analyzing beam: {beam_label} ===")
    # Make a subdirectory for each beam to avoid name collisions
    beam_outdir = os.path.join(outdir, beam_label)
    os.makedirs(beam_outdir, exist_ok=True)

    #########################
    # (a) Number of clusters
    #########################
    if "HcalFarForwardZDCClusters.position.x" in events.fields:
        nclusters = []
        cluster_x = events["HcalFarForwardZDCClusters.position.x"]
        # Count the cluster array length per event
        for evt_i in range(len(events)):
            nclusters.append(len(cluster_x[evt_i]))
        nclusters = np.array(nclusters)

        plt.figure()
        plt.hist(nclusters, bins=20, range=(0, 20))
        plt.xlabel("Number of ZDC clusters")
        plt.yscale("log")
        plt.title(f"Beam: {beam_label}")
        plt.ylim(1)
        outfile = os.path.join(beam_outdir, f"nclust_{beam_label}.pdf")
        plt.savefig(outfile)
        plt.close()
        print("Saved cluster-count plot:", outfile)
    else:
        print("WARNING: 'HcalFarForwardZDCClusters.position.x' not in branches; skipping cluster plot.")

    #########################
    # (b) Truth-level \theta*
    #     The old script used MCParticles[2] as the Lambda (?), be sure your data matches
    #########################
    tilt = -0.025  # beam tilt from old code
    px_truth = None
    py_truth = None
    pz_truth = None

    # Check if MCParticles are available
    if all(x in events.fields for x in (
            "MCParticles.momentum.x",
            "MCParticles.momentum.y",
            "MCParticles.momentum.z"
    )):
        px_truth = events["MCParticles.momentum.x"][:, 2]
        py_truth = events["MCParticles.momentum.y"][:, 2]
        pz_truth = events["MCParticles.momentum.z"][:, 2]
        pt_truth = np.hypot(px_truth*np.cos(tilt) - pz_truth*np.sin(tilt), py_truth)
        theta_truth = np.arctan2(pt_truth, pz_truth*np.cos(tilt) + px_truth*np.sin(tilt))
    else:
        print("WARNING: MCParticles branches not fully present; skipping truth-level angle computation.")
        px_truth, py_truth, pz_truth = None, None, None
        pt_truth, theta_truth = None, None

    #########################
    # (c) Reconstructed Lambda angle, mass, vertex, etc.
    #########################
    recon_fields = [
        "ReconstructedFarForwardZDCLambdas.momentum.x",
        "ReconstructedFarForwardZDCLambdas.momentum.y",
        "ReconstructedFarForwardZDCLambdas.momentum.z",
        "ReconstructedFarForwardZDCLambdas.mass",
        "ReconstructedFarForwardZDCLambdas.referencePoint.x",
        "ReconstructedFarForwardZDCLambdas.referencePoint.z",
    ]
    have_recon = all(rf in events.fields for rf in recon_fields)
    if not have_recon:
        print("WARNING: ReconstructedFarForwardZDCLambdas.* not found; skipping recon-based plots.")
        return

    px_rec = events["ReconstructedFarForwardZDCLambdas.momentum.x"]
    py_rec = events["ReconstructedFarForwardZDCLambdas.momentum.y"]
    pz_rec = events["ReconstructedFarForwardZDCLambdas.momentum.z"]
    m_rec  = events["ReconstructedFarForwardZDCLambdas.mass"]
    # Energy from sqrt(px^2 + py^2 + pz^2 + m^2)
    E_rec = np.sqrt(px_rec**2 + py_rec**2 + pz_rec**2 + m_rec**2)
    theta_rec = np.arctan2(
        np.hypot(px_rec*np.cos(tilt) - pz_rec*np.sin(tilt), py_rec),
        pz_rec*np.cos(tilt) + px_rec*np.sin(tilt)
    )

    x_ref = events["ReconstructedFarForwardZDCLambdas.referencePoint.x"]
    z_ref = events["ReconstructedFarForwardZDCLambdas.referencePoint.z"]
    z_vtx = z_ref*np.cos(tilt) + x_ref*np.sin(tilt)

    ###############
    # \theta^* plots
    ###############
    # We replicate the old code that made a triple-panel figure
    fig, axs = plt.subplots(1, 3, figsize=(24, 8))
    fig.suptitle(f"Beam: {beam_label} — Reconstructed vs Truth angles")

    # Panel 1: scatter of truth vs recon
    plt.sca(axs[0])
    if (theta_truth is not None) and (theta_rec is not None):
        plt.scatter(
            ak.to_numpy(ak.flatten(theta_truth)) * 1000.0,
            ak.to_numpy(ak.flatten(theta_rec)) * 1000.0,
            s=2
        )
        plt.xlabel(r"$\theta^{*\mathrm{truth}}_{\Lambda}$ [mrad]")
        plt.ylabel(r"$\theta^{*\mathrm{recon}}_{\Lambda}$ [mrad]")
        plt.xlim(0, 3.2)
        plt.ylim(0, 3.2)
    else:
        plt.text(0.5, 0.5, "MC truth or recon angles not available", ha='center')

    # Panel 2: histogram of residual (theta_rec - theta_truth)
    plt.sca(axs[1])
    if (theta_truth is not None) and (theta_rec is not None):
        dtheta = ak.to_numpy(ak.flatten(theta_rec - theta_truth)) * 1000.0
        yvals, xedges, _ = plt.hist(dtheta, bins=50, range=(-1, 1))
        bc = 0.5 * (xedges[:-1] + xedges[1:])

        # Gaussian fit
        from scipy.optimize import curve_fit
        mask = np.abs(bc) < 0.3
        p0 = [np.max(yvals), 0.0, 0.05]
        sigma = np.sqrt(yvals[mask]) + (yvals[mask] == 0)
        try:
            coeff, cov = curve_fit(gauss, bc[mask], yvals[mask], p0=p0, sigma=sigma, maxfev=10000)
            xx = np.linspace(-1, 1, 200)
            plt.plot(xx, gauss(xx, *coeff), color='tab:orange', label=f"Fit: sigma={coeff[2]:.3f} mrad")
            plt.legend()
        except:
            pass

        plt.xlabel(r"$\theta^{*\mathrm{rec}}_{\Lambda} - \theta^{*\mathrm{truth}}_{\Lambda}$ [mrad]")
        plt.ylabel("Counts")
    else:
        plt.text(0.5, 0.5, "No MC truth or recon angles", ha='center')

    # Panel 3: resolution vs. ??? (In the old code, they looped over momenta. We have only 1 beam here.)
    # We can store the fitted sigma for this beam:
    plt.sca(axs[2])
    plt.text(0.1, 0.5,
             "In the old script, this panel summarized\n"
             "the sigma vs momentum or beam.\n"
             "Here, we have only one beam in isolation.\n"
             "Check combined plot for multi-beam comparison.",
             ha='left')
    plt.axis('off')

    fig.tight_layout()
    outfile = os.path.join(beam_outdir, f"thetastar_{beam_label}.pdf")
    plt.savefig(outfile)
    plt.close()
    print("Saved:", outfile)

    ######################
    # (d) z-vtx resolution
    ######################
    # Old code again used a triple-panel figure
    fig, axs = plt.subplots(1, 3, figsize=(24, 8))
    fig.suptitle(f"Beam: {beam_label} — Lambda vertex reconstruction")

    # Panel 1: scatter of truth vs recon for z-vtx
    plt.sca(axs[0])
    if "MCParticles.vertex.z" in events.fields:
        z_truth = events["MCParticles.vertex.z"][:, 3]  # index 3? (the old script used 3 for the 4th particle)
        plt.scatter(
            ak.to_numpy(ak.flatten(z_truth / 1000.0)),
            ak.to_numpy(ak.flatten(z_vtx / 1000.0)),
            s=2
        )
        plt.xlabel(r"$z^{\mathrm{truth}}_{\mathrm{vtx}}$ [m]")
        plt.ylabel(r"$z^{\mathrm{recon}}_{\mathrm{vtx}}$ [m]")
        plt.xlim(0, 40)
        plt.ylim(0, 40)
    else:
        plt.text(0.5, 0.5, "MCParticles.vertex.z not found", ha='center')

    # Panel 2: distribution of residual
    plt.sca(axs[1])
    if "MCParticles.vertex.z" in events.fields:
        z_truth = events["MCParticles.vertex.z"][:, 3]
        dz = (z_vtx - z_truth) / 1000.0
        dz_flat = ak.to_numpy(ak.flatten(dz))
        yvals, xedges, _ = plt.hist(dz_flat, bins=50, range=(-10, 10))
        bc = 0.5*(xedges[:-1] + xedges[1:])

        from scipy.optimize import curve_fit
        mask = np.abs(bc) < 5
        p0 = [np.max(yvals), 0.0, 1.0]
        sigma = np.sqrt(yvals[mask]) + (yvals[mask] == 0)
        try:
            coeff, cov = curve_fit(gauss, bc[mask], yvals[mask], p0=p0, sigma=sigma, maxfev=10000)
            xx = np.linspace(-5, 5, 200)
            plt.plot(xx, gauss(xx, *coeff), color='tab:orange', label=f"Fit: sigma={coeff[2]:.2f} m")
            plt.legend()
        except:
            pass
        plt.xlabel(r"$z^{\mathrm{recon}}_{\mathrm{vtx}} - z^{\mathrm{truth}}_{\mathrm{vtx}}$ [m]")
        plt.ylabel("Counts")
    else:
        plt.text(0.5, 0.5, "No vertex info", ha='center')

    # Panel 3: In the old script, they did sigma vs momentum again. We'll just say single beam.
    plt.sca(axs[2])
    plt.text(0.1, 0.5,
             "Single-beam case:\n"
             "No momentum dependence here.\n"
             "Check combined step if needed.",
             ha='left')
    plt.axis('off')

    fig.tight_layout()
    outfile = os.path.join(beam_outdir, f"zvtx_{beam_label}.pdf")
    plt.savefig(outfile)
    plt.close()
    print("Saved:", outfile)

    ######################
    # (e) Reconstructed mass
    ######################
    # Similar approach to the old code
    lambda_mass_pdg = 1.115683
    fig, axs = plt.subplots(1, 2, figsize=(16, 8))
    fig.suptitle(f"Beam: {beam_label} — Reconstructed Lambda mass")

    plt.sca(axs[0])
    mass_flat = ak.to_numpy(ak.flatten(m_rec))
    yvals, xedges, _ = plt.hist(mass_flat, bins=100, range=(1.0, 1.25))
    plt.axvline(lambda_mass_pdg, ls='--', color='tab:green', lw=2)
    plt.text(lambda_mass_pdg + 0.005, max(yvals)*1.05, "PDG mass", color='tab:green')
    plt.xlabel(r"$m_{\Lambda}^{\mathrm{recon}}$ [GeV]")
    plt.ylabel("Counts")

    bc = 0.5*(xedges[:-1] + xedges[1:])
    from scipy.optimize import curve_fit
    # Fit around ±0.05 GeV from PDG
    mask = (bc > lambda_mass_pdg-0.05) & (bc < lambda_mass_pdg+0.05)
    p0 = [np.max(yvals), lambda_mass_pdg, 0.04]
    sigma = np.sqrt(yvals[mask]) + (yvals[mask] == 0)
    try:
        coeff, cov = curve_fit(gauss, bc[mask], yvals[mask], p0=p0, sigma=sigma, maxfev=10000)
        xx = np.linspace(0.8, 1.3, 200)
        plt.plot(xx, gauss(xx, *coeff), color='tab:orange',
                 label=fr"Gauss fit: $\sigma={coeff[2]:.3f}$ GeV")
        plt.legend()
        # debug: print("mass fit sigma, error:", coeff[2], np.sqrt(cov[2][2]))
    except:
        pass

    plt.sca(axs[1])
    # For the old code, they looped over p in 'momenta' and plotted sigmas. We only have one beam here.
    plt.text(0.1, 0.5,
             "Single beam: see left panel for mass distribution & fit.\n"
             "Combined or multiple beam energies can be placed here.\n",
             ha='left')
    plt.axis('off')

    fig.tight_layout()
    outfile = os.path.join(beam_outdir, f"lambda_mass_{beam_label}.pdf")
    plt.savefig(outfile)
    plt.close()
    print("Saved:", outfile)

    ######################
    # (f) CM angle resolution (neutron)
    ######################
    # Old code: uses TPyROOT + TLorentzVector. We can replicate in Python with vector, or just keep it awkward-based.
    # For brevity, we will do exactly what the old code does, using python's "boost" approach is possible,
    # or we keep the same method of partial usage of ROOT. We'll do something simple here to preserve the logic.

    try:
        import ROOT
        have_mc = all(x in events.fields for x in (
            "MCParticles.momentum.x",
            "MCParticles.momentum.y",
            "MCParticles.momentum.z",
            "MCParticles.mass",
        ))
        have_recon_cm = all(x in events.fields for x in (
            "ReconstructedFarForwardZDCLambdaDecayProductsCM.PDG",
            "ReconstructedFarForwardZDCLambdaDecayProductsCM.momentum.x",
            "ReconstructedFarForwardZDCLambdaDecayProductsCM.momentum.y",
            "ReconstructedFarForwardZDCLambdaDecayProductsCM.momentum.z"
        ))
    except ImportError:
        have_mc = False
        have_recon_cm = False
        print("Cannot import ROOT => skipping CM angle resolution plots.")

    if have_mc and have_recon_cm:
        pdg_cm = events["ReconstructedFarForwardZDCLambdaDecayProductsCM.PDG"]
        pxcm   = events["ReconstructedFarForwardZDCLambdaDecayProductsCM.momentum.x"]
        pycm   = events["ReconstructedFarForwardZDCLambdaDecayProductsCM.momentum.y"]
        pzcm   = events["ReconstructedFarForwardZDCLambdaDecayProductsCM.momentum.z"]

        # Full MC:
        px_mc = events["MCParticles.momentum.x"]
        py_mc = events["MCParticles.momentum.y"]
        pz_mc = events["MCParticles.momentum.z"]
        m_mc  = events["MCParticles.mass"]
        E_mc  = np.sqrt(px_mc**2 + py_mc**2 + pz_mc**2 + m_mc**2)

        # We'll replicate the old indexing:
        # L = MCParticles index=2, n = index=3.
        # Then boost to Lambda rest frame to get the neutron's angles.
        phi_resid = []
        theta_resid = []

        for i_evt in range(len(events)):
            # If event has the right shape
            if (len(px_mc[i_evt]) < 4) or (len(pxcm[i_evt]) < 2):
                continue
            l_vec = ROOT.TLorentzVector(
                px_mc[i_evt][2],
                py_mc[i_evt][2],
                pz_mc[i_evt][2],
                E_mc[i_evt][2]
            )
            n_vec = ROOT.TLorentzVector(
                px_mc[i_evt][3],
                py_mc[i_evt][3],
                pz_mc[i_evt][3],
                E_mc[i_evt][3]
            )
            # Boost n -> Lambda rest frame
            n_cm_truth = n_vec.Clone()
            l_boost = -l_vec.BoostVector()
            n_cm_truth.Boost(l_boost)
            phi_truth  = n_cm_truth.Phi()
            theta_truth = n_cm_truth.Theta()

            # Reconstructed neutron in the CM
            # PDG=2112 => neutron
            # Indices might vary if there's multiple decay products.
            # We'll find the neutron row
            pdg_evt = pdg_cm[i_evt]
            px_evt  = pxcm[i_evt]
            py_evt  = pycm[i_evt]
            pz_evt  = pzcm[i_evt]
            # We only look at the entry where PDG=2112
            # If none found => skip
            matches = np.where(pdg_evt == 2112)[0]
            if len(matches) == 0:
                continue
            # For simplicity, assume 1 neutron => take first match
            idx_n = matches[0]
            n_vec_rec = ROOT.TLorentzVector(
                px_evt[idx_n], py_evt[idx_n], pz_evt[idx_n],
                np.sqrt(px_evt[idx_n]**2 + py_evt[idx_n]**2 + pz_evt[idx_n]**2 + 0.9396**2) # approximate neutron mass
            )
            # Boost to the same Lambda rest frame as the old code did
            n_vec_rec.Boost(l_boost)
            phi_rec   = n_vec_rec.Phi()
            theta_rec = n_vec_rec.Theta()

            # In the old code, they do: (phi_rec - phi_truth)* sin(thetacmtruth)
            # and (theta_rec - theta_truth).
            # We'll replicate exactly:
            dphi = (phi_rec - phi_truth)*np.sin(theta_truth)
            dtheta = (theta_rec - theta_truth)

            phi_resid.append(dphi)
            theta_resid.append(dtheta)

        # final hist
        phi_resid = np.array(phi_resid)*1000.0
        theta_resid = np.array(theta_resid)*1000.0

        plt.figure()
        plt.hist(phi_resid, bins=100, range=(-300,300))
        plt.xlabel(r"$(\phi^n_{cm,rec} - \phi^n_{cm,truth}) \times \sin(\theta^n_{cm,truth})$ [mrad]")
        outfile = os.path.join(beam_outdir, f"neutron_phi_cm_res_{beam_label}.pdf")
        plt.savefig(outfile)
        plt.close()
        print("Saved:", outfile)

        plt.figure()
        plt.hist(theta_resid, bins=100, range=(-1000,1000))
        plt.xlabel(r"$\theta^n_{cm,rec} - \theta^n_{cm,truth}$ [mrad]")
        outfile = os.path.join(beam_outdir, f"neutron_theta_cm_res_{beam_label}.pdf")
        plt.savefig(outfile)
        plt.close()
        print("Saved:", outfile)
    else:
        print(f"Skipping CM neutron-angles for beam {beam_label} (missing branches or ROOT).")


#######################
# 4) Combined plots
#######################
def plot_combined_theta(arrays_dict, outdir):
    """
    Reproduce the "combined scatter" of all beams in one figure,
    similar to how the old code combined momenta=100..275.
    """
    # Just replicate the relevant piece from old script:
    #   x += list(theta_truth)
    #   y += list(theta_rec)
    #   then scatter
    #   then histogram residual
    #   then do sigma fit vs. beam ?

    all_theta_truth = []
    all_theta_recon = []
    beam_labels = sorted(arrays_dict.keys())

    for beam in beam_labels:
        events = arrays_dict[beam]
        # We need the same checks
        if not all(x in events.fields for x in [
            "MCParticles.momentum.x",
            "MCParticles.momentum.y",
            "MCParticles.momentum.z",
            "ReconstructedFarForwardZDCLambdas.momentum.x",
            "ReconstructedFarForwardZDCLambdas.momentum.y",
            "ReconstructedFarForwardZDCLambdas.momentum.z"
        ]):
            print(f"Skipping combined theta for {beam} - missing branches.")
            continue

        tilt = -0.025
        px_truth = events["MCParticles.momentum.x"][:, 2]
        py_truth = events["MCParticles.momentum.y"][:, 2]
        pz_truth = events["MCParticles.momentum.z"][:, 2]
        pt_truth = np.hypot(px_truth*np.cos(tilt) - pz_truth*np.sin(tilt), py_truth)
        th_truth = np.arctan2(pt_truth, pz_truth*np.cos(tilt) + px_truth*np.sin(tilt))

        px_rec = events["ReconstructedFarForwardZDCLambdas.momentum.x"]
        py_rec = events["ReconstructedFarForwardZDCLambdas.momentum.y"]
        pz_rec = events["ReconstructedFarForwardZDCLambdas.momentum.z"]
        th_rec = np.arctan2(
            np.hypot(px_rec*np.cos(tilt) - pz_rec*np.sin(tilt), py_rec),
            pz_rec*np.cos(tilt) + px_rec*np.sin(tilt)
        )

        # Flatten
        all_theta_truth.append(ak.to_numpy(ak.flatten(th_truth)))
        all_theta_recon.append(ak.to_numpy(ak.flatten(th_rec)))

    # Combine everything
    if len(all_theta_truth) == 0:
        print("No combined data for multi-beam theta star!")
        return
    all_theta_truth = np.concatenate(all_theta_truth)
    all_theta_recon = np.concatenate(all_theta_recon)

    fig, axs = plt.subplots(1, 3, figsize=(24, 8))
    fig.suptitle("Combined \u03B8* (all beams)")

    # Panel 1: scatter
    plt.sca(axs[0])
    plt.scatter(all_theta_truth*1e3, all_theta_recon*1e3, s=2)
    plt.xlabel(r"$\theta^{*\mathrm{truth}}_{\Lambda}$ [mrad]")
    plt.ylabel(r"$\theta^{*\mathrm{recon}}_{\Lambda}$ [mrad]")
    plt.xlim(0, 3.2)
    plt.ylim(0, 3.2)

    # Panel 2: residual hist + fit
    dtheta = (all_theta_recon - all_theta_truth)*1e3
    plt.sca(axs[1])
    yvals, xedges, _ = plt.hist(dtheta, bins=50, range=(-1,1))
    bc = 0.5*(xedges[:-1]+xedges[1:])
    from scipy.optimize import curve_fit
    mask = np.abs(bc) < 0.3
    p0 = [np.max(yvals), 0., 0.05]
    sigma = np.sqrt(yvals[mask]) + (yvals[mask] == 0)
    try:
        coeff, cov = curve_fit(gauss, bc[mask], yvals[mask], p0=p0, sigma=sigma, maxfev=10000)
        xx = np.linspace(-1, 1, 200)
        plt.plot(xx, gauss(xx, *coeff), color='tab:orange',
                 label=fr"Fit $\sigma={coeff[2]:.3f}$ mrad")
        plt.legend()
    except:
        pass
    plt.xlabel(r"$\theta^{*\mathrm{rec}}_{\Lambda} - \theta^{*\mathrm{truth}}_{\Lambda}$ [mrad]")
    plt.ylabel("Counts")

    # Panel 3: Usually old code did sigma vs momentum. We'll do beam-labeled:
    # We'll do a quick loop over each beam, fit the local dtheta, store the sigma
    plt.sca(axs[2])
    beam_sigmas = []
    beam_errs   = []
    for beam in beam_labels:
        events = arrays_dict[beam]
        # same checks
        if not all(x in events.fields for x in [
            "MCParticles.momentum.x",
            "MCParticles.momentum.y",
            "MCParticles.momentum.z",
            "ReconstructedFarForwardZDCLambdas.momentum.x",
            "ReconstructedFarForwardZDCLambdas.momentum.y",
            "ReconstructedFarForwardZDCLambdas.momentum.z"
        ]):
            beam_sigmas.append((beam, 0., 0.))
            continue
        tilt = -0.025
        px_t = events["MCParticles.momentum.x"][:, 2]
        py_t = events["MCParticles.momentum.y"][:, 2]
        pz_t = events["MCParticles.momentum.z"][:, 2]
        pt_t = np.hypot(px_t*np.cos(tilt) - pz_t*np.sin(tilt), py_t)
        th_t = np.arctan2(pt_t, pz_t*np.cos(tilt) + px_t*np.sin(tilt))

        px_r = events["ReconstructedFarForwardZDCLambdas.momentum.x"]
        py_r = events["ReconstructedFarForwardZDCLambdas.momentum.y"]
        pz_r = events["ReconstructedFarForwardZDCLambdas.momentum.z"]
        th_r = np.arctan2(
            np.hypot(px_r*np.cos(tilt) - pz_r*np.sin(tilt), py_r),
            pz_r*np.cos(tilt) + px_r*np.sin(tilt)
        )
        dth = ak.to_numpy(ak.flatten(th_r - th_t))*1e3
        hh, xx = np.histogram(dth, bins=100, range=(-1,1))
        bc = 0.5*(xx[:-1]+xx[1:])
        mask = np.abs(bc)<0.3
        p0 = [np.max(hh), 0., 0.06]
        from scipy.optimize import curve_fit
        sig = np.sqrt(hh[mask]) + (hh[mask] == 0)
        try:
            c, v = curve_fit(gauss, bc[mask], hh[mask], p0=p0, sigma=sig, maxfev=10000)
            beam_sigmas.append((beam, c[2], np.sqrt(v[2][2])))
        except:
            beam_sigmas.append((beam, 0., 0.))

    # Make a simple errorbar plot
    beams_ = [b[0] for b in beam_sigmas]
    svals_ = [b[1] for b in beam_sigmas]
    serrs_ = [b[2] for b in beam_sigmas]
    xvals = np.arange(len(beams_))  # or define your numeric scale
    plt.errorbar(xvals, svals_, yerr=serrs_, fmt='o', color='k')
    plt.xticks(xvals, beams_, rotation=45)
    plt.ylabel(r"$\sigma[\theta^*_{\Lambda}]$ [mrad]")
    plt.title("Resolution by beam energy")

    fig.tight_layout()
    out_fig = os.path.join(outdir, "thetastar_recon_combined.pdf")
    plt.savefig(out_fig)
    plt.close()
    print("Saved combined plot:", out_fig)


#######################
# 5) Main
#######################
def main():
    args = parse_args()
    outdir = args.outdir
    os.makedirs(outdir, exist_ok=True)

    # 5a) Group input files by beam label
    beams_map = {}  # beam_label -> list of root files
    for f in args.input_files:
        label = extract_beam_energy(f)
        beams_map.setdefault(label, []).append(f)

    # 5b) Read them into memory or do chunked approach
    #     The old script used uproot.concatenate. For large samples,
    #     you might prefer iterate. For simplicity we do the same as the old code:
    arrays_sim = {}
    for beam_label, flist in beams_map.items():
        print(f"Concatenating files for beam={beam_label}:")
        for ff in flist:
            print("   ", ff)
        # We pass a dict {filename: treename} as in the old script
        file_dict = {fname: args.tree for fname in flist}
        # Combine all events for that beam
        arrays_sim[beam_label] = uproot.concatenate(file_dict, library="ak")

    # 5c) Analyze each beam individually
    for beam_label in sorted(arrays_sim.keys()):
        analyze_beam(beam_label, arrays_sim[beam_label], outdir)

    # 5d) Optionally do a combined \theta^* plot across all beams
    if not args.no_combine and len(arrays_sim) > 1:
        plot_combined_theta(arrays_sim, outdir)


if __name__ == "__main__":
    main()
