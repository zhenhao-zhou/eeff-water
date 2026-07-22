# Ambient water example

This folder contains a representative EeFF-water input for ambient liquid water with 512 water formula units.

## Files

* `data.512h2o\\\\\\\_300k\\\\\\\_1atm`: initial configuration.
* `in.ambient`: LAMMPS input file.

## Run

Run from this directory with:

```bash
/path/to/lmp -in in.ambient
```

The input reads the final EeFF-water parameters from `../../parameters/eeff\\\\\\\_water\\\\\\\_pair\\\\\\\_style.inc`.

## Notes

This example is provided as a representative input file for the ambient-water benchmark.

