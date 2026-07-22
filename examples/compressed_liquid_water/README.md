# Compressed liquid water example

This folder contains a representative EeFF-water input for compressed liquid water.

## Files

* `data.512h2o`: initial configuration with 512 water formula units.
* `in.lewis`: LAMMPS input file.

## Run

Run from this directory with:

```bash
/path/to/lmp -in in.compressed
```

The input reads the final EeFF-water parameters from `../../parameters/eeff\_water\_pair\_style.inc`.

## Notes

This example is provided as a representative input file for compressed-liquid-water simulations outside the training range.

