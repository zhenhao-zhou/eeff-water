# EeFF-water supplementary code and data

This archive provides the modified LAMMPS source files, EeFF-water parameter information, figure source data, and representative simulation inputs associated with the manuscript.

## Contents

- `src/`: modified LAMMPS source files required for the EeFF-water simulations.
- `lammps_patch/`: base LAMMPS version information and installation notes.
- `parameters/`: final EeFF-water parameter file and parameter-order information.
- `examples/`: representative LAMMPS input and configuration files for benchmark simulations.
- `Source_Data/`: numerical source data used for figure generation.

## Required software

The simulations require a LAMMPS build containing the modified eFF source files in `src/`. The base LAMMPS version and build notes are provided in `lammps_patch/`.

## Parameters

The final EeFF-water parameter set should be provided in `parameters/` together with the parameter order and units used in the manuscript. Keep the atom-type convention consistent with the LAMMPS input files:

- O nuclei: type 1
- H nuclei: type 2
- electron pairs: type 3

## Example simulations

The `examples/` directory should contain representative production inputs and initial configurations for the main benchmark categories reported in the manuscript:

- ambient liquid water
- ambient proton-defect water
- compressed liquid water
- superionic water
- mixed close-packed oxygen-lattice water, if reported as a main benchmark

Detailed notes on the minimum required case-study inputs are provided in `examples/REQUIRED_INPUTS.md`.

## Source data

The `Source_Data/` directory contains numerical data used to generate the manuscript and supplementary figures. Where full trajectories are too large to include, processed source data and representative input files are provided, and the complete large data should be deposited in a public archive.

## Citation

If you use these files, please cite the associated manuscript. A machine-readable citation template is provided in `CITATION.cff` and should be updated after publication.
