# Building the EeFF-water LAMMPS executable

The supplied source modifications target LAMMPS 2 Aug 2023 - Update 3 and the legacy make build workflow.

## 1. Obtain the base source

Obtain a clean copy of the official LAMMPS 2 Aug 2023 - Update 3 source release.

## 2. Copy the EeFF-water source files

From the repository root, copy the files under `src/` into the base LAMMPS `src/` directory, replacing files with the same names:

```bash
cp src/* /path/to/lammps/src/
```

`src/MANIFEST.tsv` records whether each file is a modified upstream file or a new EeFF-water file.

## 3. Enable required packages

Using the legacy LAMMPS make workflow:

```bash
cd /path/to/lammps/src
make yes-eff
make yes-molecule
make yes-extra-fix
```

The original production build also enabled EXTRA-DUMP, KSPACE, MANYBODY, REAXFF, and RIGID, but these packages are not required by every EeFF-water example.

## 4. Compile

Build with an MPI makefile appropriate for the local system, for example:

```bash
make -j 8 mpi
```

The archived production environment used GCC/G++ 8.5.0 and OpenMPI 4.1.6. Other compatible compilers and MPI implementations may also work but should be validated locally.

## 5. Static repository validation

Before running LAMMPS, validate the repository package:

```bash
python scripts/validate_repository.py
```

## 6. Minimal runtime test

The smallest supplied test is the fixed-nuclei six-water-cluster minimization:

```bash
cd examples/water_cluster_6
/path/to/lmp -in in.minimize
```

The other example inputs are production-scale demonstrations and can be computationally expensive.
