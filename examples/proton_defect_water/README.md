# Proton-defect water example

This folder contains representative EeFF-water inputs for ambient proton-defect water.

## Files

* `data.511h2o\\\_1h3o`: initial configuration for 511 H2O + 1 H3O+.
* `data.511h2o\\\_1oh`: initial configuration for 511 H2O + 1 OH-.
* `in.defect`: LAMMPS input file.

## Run

Run from this directory with:

```bash
/path/to/lmp -in in.defect
```

The input reads the final EeFF-water parameters from `../../parameters/eeff\\\_water\\\_pair\\\_style.inc`.

## Notes

This example is provided as a representative input file for proton-defect water simulations.



