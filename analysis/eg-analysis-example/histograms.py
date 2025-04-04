#!/usr/bin/env python3
"""
Defines histogram containers and helper functions for filling them
using scikit-hep/hist.

In addition to the original kinematic histograms (TDIS_Q2, TDIS_xbj, etc.),
we now add Lambda kinematic histograms.
"""

import hist
from hist import Hist
import numpy as np
import awkward as ak  # For vectorized operations on awkward arrays


def create_histograms():
    """
    Create and return a dictionary of histograms using scikit-hep/hist.
    """
    # 1D histograms for DIS kinematics
    h_Q2 = Hist.new.Reg(100, 0, 100, name="Q2", label="Q^2 (GeV^2)").Double()
    h_xbj = Hist.new.Reg(100, 0, 1, name="xbj", label="x_Bj").Double()
    h_y   = Hist.new.Reg(100, 0, 1, name="y", label="y").Double()
    h_t   = Hist.new.Reg(75, -1.5, 0, name="t", label="-t (GeV^2)").Double()

    # 2D histograms for DIS kinematics
    h_xbj_Q2_2d = (
        Hist.new
        .Reg(100, 0, 1, name="xbj", label="x_Bj")
        .Reg(100, 0, 100, name="Q2", label="Q^2 (GeV^2)")
        .Double()
    )
    h_pkx_pky_2d = (
        Hist.new
        .Reg(100, -5, 5, name="pkx", label="pK_x (GeV/c)")
        .Reg(100, -5, 5, name="pky", label="pK_y (GeV/c)")
        .Double()
    )

    # New histograms for Lambda kinematics
    h_Elamb = Hist.new.Reg(100, 0, 100, name="Elamb", label="Lambda Energy (GeV)").Double()
    h_lambda_pt = Hist.new.Reg(100, 0, 10, name="lambda_pt", label="Lambda pT (GeV/c)").Double()
    h_lambda_rapidity = Hist.new.Reg(100, -5, 5, name="lambda_rapidity", label="Lambda Rapidity").Double()

    return {
        "Q2": h_Q2,
        "xbj": h_xbj,
        "y": h_y,
        "t": h_t,
        "xbj_Q2_2d": h_xbj_Q2_2d,
        "pkx_pky_2d": h_pkx_pky_2d,
        # Lambda histograms:
        "Elamb": h_Elamb,
        "lambda_pt": h_lambda_pt,
        "lambda_rapidity": h_lambda_rapidity,
    }


def fill_histograms(hists, chunk_evnts, chunk_process):
    """
    Fill each histogram with data from the chunks.

    Parameters:
      hists: dictionary of histograms.
      chunk_evnts: dictionary of numpy/awkward arrays from the Evnts tree.
      chunk_process: dictionary of arrays from the Process tree.
    """
    # Fill 1D DIS histograms from the Evnts tree
    Q2_vals = chunk_process["TDIS_Q2"]
    xbj_vals = chunk_process["TDIS_xbj"]
    y_vals = chunk_process["TDIS_y"]
    t_vals = chunk_process["TDIS_t"]

    hists["Q2"].fill(Q2=Q2_vals)
    hists["xbj"].fill(xbj=xbj_vals)
    hists["y"].fill(y=y_vals)
    hists["t"].fill(t=t_vals)

    # Fill 2D DIS histograms
    hists["xbj_Q2_2d"].fill(xbj=xbj_vals, Q2=Q2_vals)
    hists["pkx_pky_2d"].fill(pkx=chunk_evnts["pkx_Lab"], pky=chunk_evnts["pky_Lab"])

    # --- New: Fill Lambda kinematics histograms ---

    # 1. Lambda energy from the Process tree
    # Use the branch "ElambE_Lab" from the Process tree.
    Elamb_vals = chunk_process["ElambE_Lab"]
    hists["Elamb"].fill(Elamb=Elamb_vals)

    # 2. Lambda pT and rapidity from the Evnts tree "lamb_scat" branch.
    # We assume "lamb_scat" is an awkward array of TLorentzVectors with fields "px", "py", "pz", "E".
    lamb = chunk_evnts["lamb_scat"]

    # If lamb is not already an awkward array with named fields, you might need to adjust this.
    # Compute transverse momentum: pT = sqrt(px^2 + py^2)
    # Compute rapidity: y = 0.5*ln((E+pz)/(E-pz))
    # (Be sure to handle division by zero if necessary.)

    # Using awkward arrays for vectorized calculations:
    lambda_px = lamb["px"]
    lambda_py = lamb["py"]
    lambda_pz = lamb["pz"]
    lambda_E  = lamb["E"]

    lambda_pt = np.sqrt(lambda_px**2 + lambda_py**2)
    # Protect against division by zero (avoid log of zero)
    # Here we add a small epsilon.
    eps = 1e-9
    lambda_rapidity = 0.5 * np.log((lambda_E + lambda_pz + eps) / (lambda_E - lambda_pz + eps))

    hists["lambda_pt"].fill(lambda_pt=lambda_pt)
    hists["lambda_rapidity"].fill(lambda_rapidity=lambda_rapidity)
