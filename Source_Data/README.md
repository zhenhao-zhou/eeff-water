# Processed source data included in this release

This directory contains the processed numerical datasets that were present in the supplied EeFF-water archive. No missing data were reconstructed during repository cleanup.

## Current coverage

### `RDF/`

Ambient liquid-water O-O, O-H, and H-H radial distribution functions for EeFF and comparison models. Files with a `.csv` extension are normalized as two-column comma-separated files with the header:

```text
r_A,g_r
```

### `RDF_H3O/`

H3O+-O radial distribution data for EeFF, AIMD/CPMD references, and comparison models. These files retain their original text-column formats where the extension is `.txt`.

### `RDF_OH/`

OH--O radial distribution data for EeFF, AIMD/CPMD references, and comparison models.

### `density_250_340K_three_trajectories.xlsx`

Processed density data over 250-340 K, including the three-trajectory EeFF results and comparison/reference values contained in the supplied workbook.

### `diffusion_H2O_H3O_OH_1020K.xlsx`

Processed diffusion data for water and proton-defect systems, including compressed-water values contained in the supplied workbook.

### `Superionic_water/BCC/D_T_400_4000_dph2opd.txt`

The superionic-water comparison-model dataset that was present in the supplied archive.

## Scope limitation

The supplied archive does not contain every processed dataset needed to reconstruct every manuscript and Supporting Information figure, and it does not contain full molecular-dynamics trajectories. In particular, complete EeFF bcc/fcc/hcp temperature-scan source tables, vibrational-spectrum source data, ionic-population sensitivity/time-series data, and all XRD source arrays were not present in the supplied archive. These values are therefore not fabricated or inferred in this cleaned release.

For a final archival release, any additional processed datasets that are intended to support those figures should be added here with explicit column definitions and units.
