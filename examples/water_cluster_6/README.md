# Six-water-cluster electronic minimization example

This directory contains a minimal EeFF-water test for a perturbed six-water cluster.

## Files

- `data.perturbed_6_1`: six-water initial configuration.
- `in.minimize`: fixed-nuclei electronic minimization input.

## Run

```bash
/path/to/lmp -in in.minimize
```

The nuclear coordinates are fixed while the electron-pair center and radius degrees of freedom are minimized. This is the smallest supplied runtime test of the modified EeFF-water LAMMPS build.
