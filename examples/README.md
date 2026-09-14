# MCcol / gpMC Examples Directory

This directory contains standalone, self-contained examples demonstrating the core capabilities of `gpMC` (General Purpose Monte Carlo).

Each subdirectory is configured with its own input files, initial configuration data, and an executable runner script (`run.sh`) that checks for the compiled binary and executes the simulation.

---

## Directory Index

| Directory | Ensemble | Potential | Initial Format | Key Features Illustrated |
| :--- | :--- | :--- | :--- | :--- |
| [`01_nvt_silica_morse_dlpoly/`](./01_nvt_silica_morse_dlpoly/) | $NVT$ (Canonical) | Morse + Ewald | DL_POLY `CONFIG` | Silicalite crystal (4608 atoms), Ewald reciprocal sum, link-cell domain decomposition, spline potential table. |
| [`02_npt_silica_morse_dlpoly/`](./02_npt_silica_morse_dlpoly/) | $NpT$ (Isobaric) | Morse + Ewald | DL_POLY `CONFIG` | Isotropic box volume trial moves (`"isotr"`), pressure coupling at $P = 166.054\ \text{bar}$. |
| [`03_nvt_silica_lammps_input/`](./03_nvt_silica_lammps_input/) | $NVT$ (Canonical) | Morse + Ewald | LAMMPS `data.atoms` | Native LAMMPS data parser (flexible column support), outputting `last.conf` and `gpMC.lammpstrj`. |
| [`04_restart_checkpoint/`](./04_restart_checkpoint/) | $NVT$ (Canonical) | Morse + Ewald | DL_POLY `CONFIG` | Hot-restart workflow: Phase 1 creates `dump*.dmp`, Phase 2 resumes seamlessly using `restart.dmp`. |
| [`05_lammps_md_comparison/`](./05_lammps_md_comparison/) | Reference MD | Morse / LJ | LAMMPS `data.atoms` | Benchmark LAMMPS scripts and reference profiles ($g(r)$, MSD, thermodynamics) for cross-validation. |

---

## Quickstart Guide

### 1. Build the Executable
Ensure the binary `bin/gpMC.exe` is compiled from the root directory:
```bash
# Compile with GNU Fortran (recommended)
make -f Makefile.gfortran clean all

# Or compile with Intel Fortran
make clean all
```

### 2. Run an Example
Navigate to any example directory and execute its runner script:
```bash
cd examples/01_nvt_silica_morse_dlpoly
./run.sh
```

---

## Common Output Files

All examples produce standard simulation outputs in the working directory:

- **`thermoins.dat`**: Instantaneous thermodynamic values (sweep number, acceptance %, $E_{\text{tot}}$, $E_{\text{sr}}$, $E_{\text{vdw}}$, $E_{\text{fourier}}$, $E_{\text{self}}$, $E_{\text{coul}}$, volume).
- **`thermoaver.dat`**: Block-averaged values over blocks of size `nb`, with accumulated error statistics.
- **`gmix.dat`**: Multicomponent partial pair correlation functions $g_{ij}(r)$ computed on a radial grid of resolution `deltagr`.
- **`CONFIG.last`** / **`last.conf`**: Final atomic coordinates formatted for DL_POLY 2 or LAMMPS respectively.
- **`gpMC.lammpstrj`**: Trajectory in LAMMPS custom dump format, compatible with visualization tools like **VMD** and **OVITO**.
- **`dump<YYYYMMDDHHMM>.dmp`**: Complete binary state checkpoint for resuming calculations.
