# Analysis Scripts Documentation

## analysis/zdc-lambda/

This directory contains analysis scripts related to Lambda hyperons in the Zero Degree Calorimeter (ZDC).

- [analysis/zdc-lambda/](https://github.com/JeffersonLab/meson-structure/blob/main/analysis/zdc-lambda/) - 
  Lambda reconstruction analysis using EICrecon full sim-recon output 



### Files:


- [analysis/zdc-lambda/original_lambda_plots.py](https://github.com/JeffersonLab/meson-structure/blob/main/analysis/zdc-lambda/original_lambda_plots.py) - Original ZDC Lambda analysis
  Original analysis file for ZDC reconstructed lambda hyperons. This script doesn't work with meson-structure sim-recon files but is kept for reference

- [analysis/zdc-lambda/ms_lambda_plots.py](https://github.com/JeffersonLab/meson-structure/blob/main/analysis/zdc-lambda/ms_lambda_plots.py) - Refactored Lambda reconstruction analysis
  Refactored Lambda reconstruction analysis script, reproducing the same plots as the old code. This script can handle multiple input files with different beam energies.

- [analysis/zdc-lambda/reco-data.py](https://github.com/JeffersonLab/meson-structure/blob/main/analysis/zdc-lambda/reco-data.py) - Analysis of reconstructed particles
  Analyzes reconstructed particles and plots angular distributions. It creates histograms and scatter plots showing the angular distributions of various particles.


## analysis/edm4eic-metadata/

This directory contains scripts for extracting and analyzing metadata from EDM4eic files.

### Files:

- [analysis/edm4eic-metadata/true_value_analysis.py](https://github.com/JeffersonLab/meson-structure/blob/main/analysis/edm4eic-metadata/true_value_analysis.py) - Processes event-level metadata from EDM4eic files and builds histograms

  Shows event-level metadata from EDM4eic files and builds 1D histograms of all numeric key-values. Metadata from the original event generator files are copied through the simulation chain. Important values such as true Q2, Bjorken x, etc. are extracted from special branches of the 'event' tree ("GPStringKeys" and "GPStringValues").

## analysis/eg-analysis-example/

This directory contains example analysis scripts for studying the kinematics of Lambda hyperons.

### Files:

- [analysis/eg-analysis-example/histograms.py](https://github.com/JeffersonLab/meson-structure/blob/main/analysis/eg-analysis-example/histograms.py) - Histogram definition for Lambda kinematics

  Defines histogram containers and helper functions for filling them using scikit-hep/hist. In addition to the original kinematic histograms (TDIS_Q2, TDIS_xbj, etc.), it includes Lambda kinematic histograms.

- [analysis/eg-analysis-example/main.py](https://github.com/JeffersonLab/meson-structure/blob/main/analysis/eg-analysis-example/main.py) - Main analysis script for k_lambda data

  Main script to analyze k_lambda data using uproot and scikit-hep/hist. This version processes multiple trees simultaneously (Evnts and Process) and adds Lambda kinematic histograms.

- [analysis/eg-analysis-example/plotting.py](https://github.com/JeffersonLab/meson-structure/blob/main/analysis/eg-analysis-example/plotting.py) - Plotting utilities for analysis results

  Produces and saves 1D & 2D plots from scikit-hep/hist objects. Now includes plots for Lambda kinematics.



## tools/

This directory contains utility scripts for working with ROOT files and metadata.

### Files:

- [tools/events-metadata.py](https://github.com/JeffersonLab/meson-structure/blob/main/tools/events-metadata.py) - Event metadata extractor

  Shows event level metadata from the first 10 events in a ROOT file, processing GP*Keys and GP*Values branches.

- [tools/plot_metadata.py](https://github.com/JeffersonLab/meson-structure/blob/main/tools/plot_metadata.py) - Metadata plotting tool

  Extracts event-level metadata (e.g., dis_xbj, dis_q2) from a ROOT file and creates plots of the extracted values.

- [tools/plot_metadata2.py](https://github.com/JeffersonLab/meson-structure/blob/main/tools/plot_metadata2.py) - Enhanced metadata plotting

  Shows event-level metadata from EDM4eic files and builds 1D histograms of all numeric key-values. Processes data in chunks for efficient memory usage.

- [tools/podio_metadata.py](https://github.com/JeffersonLab/meson-structure/blob/main/tools/podio_metadata.py) - PODIO metadata analyzer

  Extracts and analyzes podio_metadata from a ROOT file, including class definitions, collection relationships, and ID information.

- [tools/root-metadata.py](https://github.com/JeffersonLab/meson-structure/blob/main/tools/root-metadata.py) - ROOT metadata extractor

  Reads and processes various metadata trees (metadata, runs, podio_metadata) from a ROOT file and outputs to a JSON file.

- [tools/root-metadata2.py](https://github.com/JeffersonLab/meson-structure/blob/main/tools/root-metadata2.py) - Enhanced ROOT metadata extractor

  Similar to root-metadata.py but with additional processing for complex structures in the podio_metadata tree.

## Overall Summary

The repository contains scripts for analyzing data related to Lambda hyperons in particle physics experiments, particularly those reconstructed in Zero Degree Calorimeters (ZDC). The scripts use various Python packages like uproot, awkward arrays, and scikit-hep/hist to process ROOT files from simulations and experiments.

Key functionality includes:
1. Extraction and analysis of event metadata
2. Reconstruction of Lambda particles from proton-pion combinations
3. Comparison between MC truth and reconstructed particles
4. Analysis of Lambda decay modes and kinematics
5. Visualization of angular distributions and other kinematic variables
6. Tools for working with ROOT files and metadata

The scripts span from low-level debugging tools to complete analysis chains that produce publication-quality plots.

::: info
This page is autogenerated generated with this promt:
Based on contents of uproot analysis files in analysis folder
and its subfolders make analysis.md with a summary, what each analysis script does.

1. When describing a file with analysis do, start with full link to the file in JeffersonLab/meson-structure and title - their relative path. Here is an example:
   ``` - [analysis/zdc-lambda/original_lambda_plots.py] https://github.com/JeffersonLab/meson-structure/blob/main/analysis/zdc-lambda/original_lambda_plots.py - Original analysis file for ZDC reconstructed lambda hyperons ```
2. If .py file has a description in the beginning, e.g. the first long commend in """ - this is the best description of the analysis.
   you can generate title out of it and use it as the description.
3. Provide summary for each directory too. If directory has README.md - this is the best description for the directory, use it
4. Provide overall summary after
:::