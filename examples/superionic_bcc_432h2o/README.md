# 432-H2O bcc superionic-water example

This folder contains a representative EeFF-water input for bcc superionic water.

## Files

* `data.bcc\_432H2O`: initial configuration with 432 water formula units, represented by 432 O cores, 864 H nuclei and 1728 explicit electron-pair particles.
* `in.bcc`: LAMMPS input file for the bcc simulation.
* `log.bcc\_432H2O`: reference output log from the supplied run.

## Requirements

The input requires a LAMMPS build with the EeFF-water `eff/cut` implementation and support for `fix temp/csvr`. The final EeFF-water parameters are read from:

```bash
../../parameters/eeff\_water\_pair\_style.inc
```

## Run

Run from this directory with:

```bash
/path/to/lmp -in in.bcc
```

The input uses a 0.1 fs time step, an initial temperature-ramp stage, and a 200000-step production stage. The barostat pressure is evaluated using the nuclear kinetic contribution and total virial.

## Notes

This example is provided as a representative simulation input, not as a fast unit test. The reference log is included only to help check that the input, data file and parameter file are read consistently.

