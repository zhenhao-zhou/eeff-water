# Representative EeFF-water examples

Run each input from its own directory so that relative paths to `../../parameters/eeff_water_pair_style.inc` resolve correctly.

| Directory | System | Representative conditions | Production segment |
|---|---|---|---|
| `water_cluster_6` | perturbed (H2O)6 cluster | fixed nuclei, electronic minimization | minimization only |
| `ambient_water_512h2o` | 512 H2O | 300 K, 1 atm | 1 ns |
| `proton_defect_water` | 511 H2O + H3O+ or 511 H2O + OH- | 300 K, 1 atm | 1 ns |
| `compressed_liquid_water` | 512 H2O | 1020 K, 11 GPa | 400 ps |
| `superionic_water` | bcc/fcc/hcp oxygen lattices | 100 GPa, representative 1000 K case | 20 ps after a 0.6 ps ramp |

These inputs are representative production-style examples. They do not enumerate every independent random seed, temperature, or pressure point used in the complete benchmark study.
