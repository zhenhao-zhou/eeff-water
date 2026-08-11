# EeFF-water parameterization dataset

This directory contains the energy and nuclear-force reference data used to parameterize the EeFF-water model.

## Dataset composition

The dataset contains 1,772 configurations:

- distorted H2O monomers: 100
- neutral water clusters, (H2O)n with n = 2-6: 462
- protonated clusters containing one H3O+ defect: 300
- deprotonated clusters containing one OH- defect: 300
- clusters containing two H3O+ defects: 300
- clusters containing two OH- defects: 300
- ambient condensed-phase water configurations: 10

The 1,762 isolated molecular/cluster configurations plus the 10 condensed-phase configurations give the total of 1,772.

## File organization

Each subdirectory represents one compositional or configurational class.

Files named `data.*` are LAMMPS data files containing the nuclear and electron-pair coordinates used to initialize electronic energy minimization.

The corresponding reference energies and nuclear forces are stored in:

- `energy.txt`
- `force.txt`

Entries follow the natural numeric order of the corresponding `data.*` index within each directory.

## Reference calculations

Reference energies and nuclear forces for isolated molecular and ionic cluster configurations were calculated at the M06-2X/def2-TZVP level.

The 10 condensed-phase liquid-water configurations were selected from a published revPBE0-D3 water dataset.

Energies are expressed relative to the lowest-energy configuration within each compositional class. Absolute energies from the two electronic-structure data sources are therefore not directly mixed.

## Units

- energy: kcal mol-1
- nuclear force: kcal mol-1 Angstrom-1

Only forces acting on oxygen and hydrogen nuclei are included as reference labels.

## Electron-pair coordinates

The electron-pair coordinates and wave-packet radii in the LAMMPS data files are initialization values for electronic minimization and are not quantum-mechanical reference labels.

For each candidate EeFF parameter set, nuclear coordinates are held fixed while electron-pair center coordinates and radii are minimized before evaluating EeFF energies and nuclear forces.

## Training-domain scope

Compressed liquid water and superionic water configurations are not included in the parameterization dataset. These conditions are reserved for transferability evaluation.

## Related files

- final parameters: `../parameters/`
- modified LAMMPS source: `../src/`
