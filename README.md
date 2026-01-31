# Example code for arXiv:2510.02033

[![arXiv](https://img.shields.io/badge/arXiv-2510.02033-b31b1b.svg)](https://arxiv.org/abs/2510.02033)

This repository contains the source code and data for the paper:

**"QNM families: classification and competition"** *Zhen-Hao Yang, Liang-Bi Wu, Xiao-Mei Kuang, and Wei-Liang Qian*

The codebase implements four different methods (Pseudospectral, Time-Domain/Matrix Pencil, WKB, and Direct Integration) to analyze the Quasinormal Modes (QNMs) of gravitational decoupled hairy black holes.

## Repository Structure

The repository is organized into four main directories, each corresponding to a specific computational method discussed in the paper:

```text
.
├── Pseudospectral/             # Mathematica notebook for hyperboloidal pseudospectral method
├── Timedomain & Matrix pencil/ # Julia scripts for time evolution and QNM extraction
├── WKB/                        # WKB approximation codes
└── Direct integration/         # Shooting method implementation

```
# Detailed Usage Guide

## 1. Pseudospectral Method
Located in - `Pseudospectral/`.

Content: Contains a single Mathematica notebook (`.nb`).

Description: Implements the pseudospectral method based on the hyperboloidal coordinate regularization scheme for the master equation of the gravitational decoupled black hole.

Reproducibility: Running this notebook reproduces the QNM data presented in Table I of the paper.

## 2. Time-Domain Evolution & Matrix Pencil Method
Located in `Timedomain & Matrix pencil/`.

Prerequisites: Julia (v1.10 recommended).

Workflow:

Time Evolution: Run `FD_Evolution_ky.jl`.

Input: Requires a pre-calculated discretized tortoise coordinate inverse function CSV file (e.g., `r(rs)_grid_M1_L6_a700_lo952.csv` provided in the folder).

Output: Generates the time-domain waveform plot (`TD_M1_L6_a700_lo952.pdf`) and the waveform data (`TD_data_M1L6S1_Echo.csv`).

QNM Extraction: Run `MP_qnmExtractor_ky.jl`.

Input: Reads the waveform data generated in step 1. Note: Ensure all input/output files are in the same directory.

Output: Extracts QNM frequencies, complex amplitudes, and the Energy Fractions.

## 3. WKB Approximation
Located in `WKB/`.

Content:

WKB.m: The WKB calculation package developed by Konoplya et al.

`WKB_ky.nb`: Mathematica notebook applying WKB.m to the hairy black hole metric.

Reproducibility: Reproduces the WKB approximation data in Table I.

## 4. Direct Integration (Shooting Method)
Located in `Direct integration/`.

Content:

`general_seriesHa_13` & `general_seriesINFb_14`: Correspond to the series expansion formulas (Eq. D6 and D7) in the Appendix.

`general_asolall_13` & `general_bsolall_14`: The solved coefficients for the series expansions.

`Shooting_ky.nb`: The main notebook.

Description: The notebook loads the pre-computed series solutions to construct high-precision boundary values and performs the shooting method integration.

Reproducibility: Reproduces the eigenfunctions and potentials shown in Fig. 3.
