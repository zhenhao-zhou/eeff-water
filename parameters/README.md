# EeFF-water parameters

`eeff_water_pair_style.inc` contains the complete `pair_style eff/cut` command used by the representative examples, followed by `pair_coeff * *`.

The modified `PairEffCut::settings()` interface accepts exactly 59 numeric arguments.

## Manuscript parameters versus compatibility arguments

Arguments 0-56 are the 57 EeFF-water values that define the model reported in the manuscript and Supporting Information, including the 8 Angstrom cutoff and the optimized radius bounds for the `abs(spin) = 1` electron-pair particles.

Arguments 57-58 are legacy radius bounds retained for compatibility with `abs(spin) = 3` electronic particles in the modified pEFF code. The supplied EeFF-water training data and all representative EeFF-water simulation data files contain only spin values 0, 1, and 2. Therefore, arguments 57-58 are inactive for EeFF-water.

`parameter_order.csv` records the exact argument order read by the source code. `manuscript_parameter_mapping.csv` maps these arguments to the parameter labels used in the manuscript and marks the two inactive compatibility values explicitly.

## Use in a LAMMPS input

From an example directory two levels below the repository root:

```lammps
include ../../parameters/eeff_water_pair_style.inc
```

The include file also executes:

```lammps
pair_coeff * *
```

All representative inputs in this release read the same include file instead of duplicating the 59 numeric arguments.
