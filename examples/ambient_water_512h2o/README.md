# Ambient liquid-water example

This directory contains a representative EeFF-water input for 512 water formula units at 300 K and 1 atm.

## Files

- `data.512h2o_300k_1atm`: initial configuration.
- `in.ambient`: LAMMPS input.

## Run

```bash
/path/to/lmp -in in.ambient
```

## Protocol encoded in the input

- time step: 0.2 fs
- nuclear translational thermostat: 300 K, 10 fs damping
- electron-pair translational thermostat: 300 K, 100 fs damping
- electron-radius thermostat: 1 K, 1 fs damping
- isotropic barostat: 1 atm, 5000 fs damping
- production segment: 5,000,000 steps = 1 ns

The pressure used by the barostat is evaluated from the nuclear kinetic contribution and the total interaction virial.
