# EeFF-water parameters

`eeff_water_pair_style.inc` contains 59 full-precision numeric arguments in the order read by `PairEffCut::settings()`.

Before release, compare this file byte-for-byte/numerically with the final production input and the final SI parameter table. The current source code should also be updated to reject an incorrect argument count explicitly rather than checking only `narg < 1`.

Use from a LAMMPS input file:

```lammps
include ../../parameters/eeff_water_pair_style.inc
```

The include file also executes `pair_coeff * *`. Remove that final line if a particular example requires type-specific cutoffs.
