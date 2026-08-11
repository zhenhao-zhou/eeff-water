# EeFF-water: code, parameters and data

This repository contains the modified LAMMPS source files, final EeFF-water parameter set, parameterization dataset, processed numerical source data included with this release, and representative simulation inputs associated with the manuscript:

"A transferable explicit electron force field for molecular, ionic and superionic water"

EeFF is an explicit-electron reactive force field in which effective oxygen cores, hydrogen nuclei, and valence electron-pair wave packets are propagated as interacting particles.

## Repository contents

- `src/`: modified LAMMPS source files used by EeFF-water.
- `lammps_patch/`: base LAMMPS version and build instructions.
- `parameters/`: final EeFF-water parameters and the LAMMPS argument order.
- `training_data/`: energy and nuclear-force reference data used for parameterization.
- `examples/`: representative LAMMPS inputs and initial configurations.
- `Source_Data/`: processed numerical data included in the supplied archive.
- `scripts/validate_repository.py`: static consistency checks for the release package.
- `IMPLEMENTATION_NOTES.md`: implementation details that are not fitted parameters.

## Parameterization dataset

The parameterization dataset contains 1,772 configurations:

- 1,762 isolated molecular and ionic cluster configurations.
- 10 condensed-phase configurations containing 64 water molecules.

The dataset includes distorted water monomers, neutral water clusters, protonated and deprotonated water clusters, double-defect clusters, and ambient liquid-water configurations.

Reference energies and nuclear forces for isolated molecular and ionic clusters were calculated at the M06-2X/def2-TZVP level. The 10 condensed-phase configurations were selected from a published revPBE0-D3 water dataset. Energies are expressed relative to the lowest-energy configuration within each compositional class.

Compressed liquid water and superionic water configurations are not included in the parameterization dataset and are used only for transferability tests.

See `training_data/README.md` for file organization and units.

## EeFF-water parameters

The final LAMMPS `pair_style eff/cut` command is provided in:

```text
parameters/eeff_water_pair_style.inc
```

All representative inputs read this same include file so that a single parameter set is used throughout the repository.

The modified pair style currently accepts 59 numeric arguments. Arguments 0-56 correspond to the EeFF-water parameters reported for the model. Arguments 57-58 are retained only for compatibility with legacy spin-3 electronic particles. No spin-3 particles occur in the EeFF-water training data or representative simulation data supplied here, so these final two compatibility bounds are inactive for the reported water model.

The atom-type convention is:

- oxygen cores: type 1
- hydrogen nuclei: type 2
- electron-pair wave packets: type 3

See `parameters/README.md` and `parameters/manuscript_parameter_mapping.csv` for details.

## Building LAMMPS

The source modifications were used with LAMMPS 2 Aug 2023 - Update 3. Build instructions and the recorded compiler/MPI environment are provided in `lammps_patch/`.

## Representative examples

The `examples/` directory contains:

- fixed-nuclei electron minimization for a six-water cluster;
- 512-water ambient liquid water at 300 K and 1 atm;
- proton-defect water for H3O+ and OH- systems;
- compressed liquid water at 1020 K and 11 GPa;
- bcc, fcc, and hcp superionic-water examples at 100 GPa using a representative 1000 K production temperature.

These are representative inputs, not a complete archive of every independent trajectory used in the manuscript. See `examples/README.md` for run instructions and simulation lengths.

## Source data

`Source_Data/` contains the processed numerical datasets present in the supplied archive. The files have been normalized so that files carrying a `.csv` extension are valid comma-separated files. Complete molecular-dynamics trajectories are not included because of their size.

The supplied archive does not contain every processed dataset used for every manuscript and Supporting Information figure. No missing numerical values were reconstructed or invented during repository cleanup. See `Source_Data/README.md` for the exact coverage of this release.

## Validation

Run the repository consistency checks from the repository root:

```bash
python scripts/validate_repository.py
```

The validator checks the parameter count, training-data counts and label dimensions, absence of spin-3 particles in the supplied EeFF-water data, example file references, the ambient-water production length, and CSV formatting.

## Citation

Please cite the associated manuscript when using the EeFF-water source code, parameters, or data. A machine-readable citation record is provided in `CITATION.cff`.

## License

See `LICENSE.md` for the license scope of the modified LAMMPS source and the EeFF-water-specific files.
