# EeFF-water parameterization dataset

This directory contains the energy and force reference data used to parameterize the EeFF water model.

## Dataset composition

The dataset contains 1,772 configurations:

- 1,762 isolated molecular and ionic cluster configurations
- 10 condensed-phase configurations containing 64 water molecules

The dataset includes:

- distorted water monomers
- neutral water clusters, (H2O)n with n = 2-6
- protonated water clusters containing H3O+ defects
- deprotonated water clusters containing OH- defects
- clusters containing two hydronium defects
- clusters containing two hydroxide defects
- ambient liquid-water configurations

## File organization

Each subdirectory represents one compositional or configurational class.

Files named `data.*` are LAMMPS data files containing the nuclear and electron-pair coordinates used to initialize electronic energy minimization.

The corresponding reference energies and nuclear forces are stored in:

- `energy.txt`
- `force.txt`

The entries in these files follow the indexed order of the corresponding `data.*` files.

## Reference calculations

Reference energies and nuclear forces for isolated molecular and ionic cluster configurations were calculated at the M06-2X/def2-TZVP level.

The 10 condensed-phase liquid-water configurations were selected from a published dataset calculated at the revPBE0-D3 level.

Energies were expressed relative to the lowest-energy configuration within each compositional class. Absolute energies obtained at the two electronic-structure levels were therefore not directly compared.

## Units

- Energy: kcal mol-1
- Nuclear force: kcal mol-1 Angstrom-1

Only forces acting on oxygen and hydrogen nuclei are included as reference force labels.

## Electron-pair coordinates

The electron-pair coordinates and wave-packet radii in the LAMMPS data files are initial values used for electronic energy minimization and are not quantum-mechanical reference labels.

For each candidate EeFF parameter set, the nuclear coordinates were held fixed while the electron-pair coordinates and radii were optimized by energy minimization before evaluating the EeFF energy and nuclear forces.

## Training-domain scope

Compressed liquid water and superionic water configurations were not included in the parameterization dataset. These conditions were reserved for evaluation of model transferability.

## Related files

The final EeFF-water parameters are provided in `../parameters/`.

The modified LAMMPS source files are provided in `../src/`.
