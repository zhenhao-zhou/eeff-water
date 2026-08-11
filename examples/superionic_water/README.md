# Superionic-water examples

This directory contains separate representative EeFF-water inputs for bcc, fcc, and hcp oxygen sublattices at 100 GPa.

## Files

- `data.bcc_432h2o` and `in.bcc`: 432 H2O bcc case.
- `data.fcc_500h2o` and `in.fcc`: 500 H2O fcc case.
- `data.hcp_512h2o` and `in.hcp`: 512 H2O hcp case.

## Run

For example:

```bash
/path/to/lmp -in in.bcc
```

Use `in.fcc` or `in.hcp` for the other lattices.

## Protocol encoded in the representative inputs

- pressure: 986900 atm (approximately 100 GPa)
- time step: 0.1 fs
- electron-radius thermostat: 1 K, 1 fs damping
- nuclear/electron translational damping: 10/100 fs
- initial temperature ramp: nuclei 400 to 1000 K, electrons 200 to 500 K, 6000 steps = 0.6 ps
- representative production: nuclei 1000 K, electrons 500 K, 200000 steps = 20 ps

The manuscript temperature scan covers multiple temperatures. The three inputs here provide a single representative 1000 K case for each oxygen lattice rather than duplicating the full scan.
