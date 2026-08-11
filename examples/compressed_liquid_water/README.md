# Compressed liquid-water example

This directory contains a representative EeFF-water input for 512 water formula units at 1020 K and 11 GPa.

## Files

- `data.512h2o`: initial configuration.
- `in.compressed`: LAMMPS input.

## Run

```bash
/path/to/lmp -in in.compressed
```

## Protocol encoded in the input

- time step: 0.2 fs
- nuclear translational thermostat: 1020 K, 10 fs damping
- electron-pair translational thermostat: 1020 K, 100 fs damping
- electron-radius thermostat: 1 K, 1 fs damping
- isotropic barostat: 108561 atm (approximately 11 GPa), 5000 fs damping
- production segment: 2,000,000 steps = 400 ps

This file is a representative 11 GPa case. The complete benchmark study also includes other pressures; this directory does not claim to archive every production trajectory.
