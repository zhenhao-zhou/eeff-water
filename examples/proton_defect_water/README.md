# Proton-defect water examples

This directory contains separate runnable EeFF-water inputs for hydronium and hydroxide defects at 300 K and 1 atm.

## Files

- `data.511h2o_1h3o`: 511 H2O + 1 H3O+ initial configuration.
- `data.511h2o_1oh`: 511 H2O + 1 OH- initial configuration.
- `in.h3o`: hydronium-defect input.
- `in.oh`: hydroxide-defect input.

## Run

Hydronium:

```bash
/path/to/lmp -in in.h3o
```

Hydroxide:

```bash
/path/to/lmp -in in.oh
```

Both inputs use a 0.2 fs time step and a 5,000,000-step production segment, corresponding to 1 ns.
