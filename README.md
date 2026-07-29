# EeFF-water: code, parameters and data

This repository provides the modified LAMMPS source files, final force-field parameters, parameterization dataset, numerical source data and representative simulation inputs associated with the EeFF water model.

EeFF is an explicit-electron reactive force field in which effective oxygen cores, hydrogen nuclei and valence electron-pair wave packets are propagated as interacting particles.

## Repository contents

- `src/`: modified LAMMPS source files used for EeFF-water simulations.
- `lammps_patch/`: base LAMMPS version information and compilation notes.
- `parameters/`: final EeFF-water parameter set and parameter-order information.
- `training_data/`: energy and force reference data used for parameterization.
- `examples/`: representative LAMMPS inputs and initial configurations.
- `Source_Data/`: numerical source data used to generate the manuscript and supplementary figures.

## Parameterization dataset

The parameterization dataset contains 1,772 configurations, including 1,762 isolated molecular and ionic cluster configurations and 10 condensed-phase configurations containing 64 water molecules.

The dataset includes distorted water monomers, neutral water clusters, protonated and deprotonated water clusters, double-defect clusters and ambient liquid-water configurations.

Cluster reference energies and nuclear forces were calculated at the M06-2X/def2-TZVP level. The condensed-phase configurations were selected from a published revPBE0-D3 liquid- and solid-water dataset.

Energies were expressed relative to the lowest-energy configuration within each compositional class. Compressed and superionic water configurations were excluded from the parameterization dataset and were used to evaluate model transferability.

Additional documentation is provided in `training_data/README.md`.

## Required software

The simulations require a LAMMPS build containing the modified eFF source files provided in `src/`. Information on the base LAMMPS version and compilation procedure is provided in `lammps_patch/`.

## EeFF-water parameters

The final EeFF-water parameter set is provided in `parameters/`.

The atom-type convention used in the simulations is:

- oxygen cores: type 1
- hydrogen nuclei: type 2
- electron pairs: type 3

## Source data

The `Source_Data/` directory contains numerical values used to generate the manuscript and supplementary figures.

Processed source data are provided when complete molecular-dynamics trajectories are too large for inclusion in the repository.

## Citation

Please cite the associated manuscript when using the EeFF-water source code, parameters or data:

"A transferable explicit electron force field for molecular, ionic and superionic water"

A machine-readable citation record is provided in `CITATION.cff` and will be updated with the final publication information.

## License

Licensing information for the modified LAMMPS source files, EeFF-water parameters, parameterization dataset and processed source data is provided in `LICENSE.md`.
