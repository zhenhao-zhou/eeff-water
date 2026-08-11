# EeFF-water LAMMPS source files

This directory contains the LAMMPS source files that were modified or added for EeFF-water. It is intentionally not a complete LAMMPS distribution.

Use these files with the base release recorded in `../lammps_patch/BASE_VERSION.txt`. Copy them into the corresponding LAMMPS `src/` directory before compilation. `MANIFEST.tsv` lists the target installation path and classifies each file as a modified upstream file or a new EeFF-water file.

The source code in this release preserves the production model implementation. Repository cleanup was limited to documentation, comments, input-file consistency, and a strict argument-count check in `PairEffCut::settings()`; the fitted numerical parameter values and interaction formulas were not changed.

Implementation-level fixed restraints and numerical guards are documented in `../IMPLEMENTATION_NOTES.md`.
