# Supplemental Material for arXiv:2510.02033

[![arXiv](https://img.shields.io/badge/arXiv-2510.02033-b31b1b.svg)](https://arxiv.org/abs/2510.02033)

This repository contains the source code, data, and supplementary materials for the paper:

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
