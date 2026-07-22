# Building the EeFF-water LAMMPS executable

This archive contains the modified LAMMPS source files required for the EeFF-water simulations. The base LAMMPS version is recorded in `BASE_VERSION.txt`.

Recommended installation procedure:

1. Download the base LAMMPS version specified in `BASE_VERSION.txt`.
2. Copy the files in `src/` from this archive into the corresponding LAMMPS `src/` directory, replacing the original files where names overlap.
3. Build LAMMPS using the packages and compiler settings used in the manuscript.
4. Verify the build using the smoke-test inputs in `examples/`.

The exact executable name, compiler version, MPI version and optional acceleration packages should be recorded here before public release.
