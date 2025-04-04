#!/usr/bin/env python3
"""
Produces and saves 1D & 2D plots from scikit-hep/hist objects.
Now includes plots for Lambda kinematics.
"""

import os
import matplotlib.pyplot as plt


def make_plots(hists, outdir="plots"):
    """
    Create and save all plots.
    """
    os.makedirs(outdir, exist_ok=True)

    # 1) DIS histograms
    _plot_and_save(hists["Q2"], "TDIS_Q2 distribution", "h_Q2.png", xlabel="Q^2 (GeV^2)")
    _plot_and_save(hists["xbj"], "TDIS_xbj distribution", "h_xbj.png", xlabel="x_Bj")
    _plot_and_save(hists["y"], "TDIS_y distribution", "h_y.png", xlabel="y")
    _plot_and_save(hists["t"], "TDIS_t distribution", "h_t.png", xlabel="-t (GeV^2)")

    _plot_2d_and_save(
        hists["xbj_Q2_2d"], "TDIS_xbj vs TDIS_Q2", "h_xbj_Q2_2d.png",
        xlabel="x_Bj", ylabel="Q^2 (GeV^2)", log=True
    )
    _plot_2d_and_save(
        hists["pkx_pky_2d"], "Kaon momentum: pkx vs pky", "h_pkx_pky_2d.png",
        xlabel="pK_x (GeV/c)", ylabel="pK_y (GeV/c)", log=True
    )

    # 2) Lambda kinematic histograms
    _plot_and_save(
        hists["Elamb"], "Lambda Energy distribution", "h_Elamb.png",
        xlabel="Lambda Energy (GeV)"
    )
    _plot_and_save(
        hists["lambda_pt"], "Lambda pT distribution", "h_lambda_pt.png",
        xlabel="Lambda pT (GeV/c)"
    )
    _plot_and_save(
        hists["lambda_rapidity"], "Lambda Rapidity distribution", "h_lambda_rapidity.png",
        xlabel="Lambda Rapidity"
    )


def _plot_and_save(hist_obj, title, filename, xlabel="Value"):
    """Helper to plot a 1D histogram and save it."""
    fig, ax = plt.subplots(figsize=(6, 5))
    hist_obj.plot1d(ax=ax)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    fig.tight_layout()
    fig.savefig(os.path.join("plots", filename))
    plt.close(fig)


def _plot_2d_and_save(hist_obj, title, filename, xlabel="X", ylabel="Y", log=False):
    """Helper to plot a 2D histogram and save it."""
    fig, ax = plt.subplots(figsize=(6, 5))
    hist_obj.plot2d(ax=ax, norm="log" if log else None, cmap="viridis")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    fig.savefig(os.path.join("plots", filename))
    plt.close(fig)
