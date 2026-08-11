# EeFF-water implementation notes

This file documents implementation details that are fixed in the supplied source code and are not part of the optimized EeFF-water `pair_style` parameter vector.

## 1. Nuclear-overlap guards

`src/pair_eff_cut.cpp` contains short-range harmonic guards for extremely small nucleus-nucleus separations. They contribute to the diagnostic `eadd` energy component when the corresponding distance is below the fixed threshold:

- H-H: 0.4 Angstrom
- O-H: 0.6 Angstrom
- O-O: 1.8 Angstrom

The fixed coefficient used by these guards is 200 in LAMMPS `real` units. These values are hard-coded implementation safeguards and are not optimized EeFF-water parameters.

## 2. O-electron-pair localization restraint

The modified `bond_harmonic` implementation is a flat-bottom restraint rather than the standard LAMMPS harmonic bond. For the EeFF-water data files, `bond_coeff 1 10000.0 1.0 -0.1` gives zero bond energy inside the allowed interval and applies a quadratic penalty outside it. Because physical distances are non-negative, the effective upper localization bound for the supplied O-electron-pair bonds is 1.0 Angstrom.

This restraint is specified independently of the `pair_style eff/cut` parameter vector. The LAMMPS data files explicitly list the O-electron-pair bonds used by the restraint.

## 3. Electron-radius bounds and the 59 pair-style arguments

The water model uses electron-pair particles with `abs(spin) = 1`. Their optimized radius bounds are:

- maximum radius: 0.933901288659 Angstrom
- minimum radius: 0.340662334454 Angstrom

The modified source also retains two legacy radius bounds for `abs(spin) = 3` particles. No spin-3 particles are present in any EeFF-water training configuration or representative simulation data file in this repository. Consequently, the last two numeric arguments of the current pair-style interface are compatibility values and are inactive for the reported water model.

## 4. Constant-pressure simulations

For the supplied constant-pressure examples, the pressure passed to the barostat uses the nuclear translational temperature together with the total interaction virial. The fictitious electron translational kinetic energy is therefore not counted as a thermodynamic kinetic contribution to pressure, while interactions involving explicit electrons contribute through the virial.

## 5. Energy diagnostics

The modified pair style exposes five diagnostic components through `compute ... pair eff/cut`:

1. electron kinetic energy
2. Pauli energy
3. Coulomb energy
4. electron-radius restraint energy
5. additional short-range nuclear-overlap guard energy (`eadd`)

These diagnostics are printed in the representative inputs to make implementation-level contributions visible during testing.
