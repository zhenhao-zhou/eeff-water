# Six-water-cluster minimization example

This folder contains a minimal EeFF-water example for a perturbed six-water cluster.

## Files

- `data.perturbed_6_1`: initial cluster configuration with 6 water formula units, represented by 6 O cores, 12 H nuclei and 24 explicit electron-pair particles.
- `in.minimize`: LAMMPS input file for fixed-nuclei electron minimization.
- `log.perturbed_6_1`: reference output log from an earlier run.

## Requirements

The input requires a LAMMPS build with the EeFF-water `eff/cut` implementation. The final EeFF-water parameters are read from:

```bash
../../parameters/eeff_water_pair_style.inc
```

## Run

Run from this directory with:

```bash
/path/to/lmp -in in.minimize
```

The input keeps the nuclear coordinates fixed and minimizes the electron-center and electron-radius degrees of freedom.

## Notes

This example is intended as a small initialization and minimization test. The supplied log was generated with an earlier parameter set and is retained only for provenance, so it should not be used as a strict numerical regression reference for the final parameter set.