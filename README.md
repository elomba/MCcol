# MCcol: General Purpose Monte Carlo Code (`gpMC`)

[![Fortran](https://img.shields.io/badge/Language-Fortran_90%2F95-734f96.svg)](https://fortran-lang.org/)
[![C](https://img.shields.io/badge/Language-C99-555555.svg)](https://en.wikipedia.org/wiki/C_(programming_language))
[![Compiler](https://img.shields.io/badge/Compiler-GFortran%20%7C%20Intel_ifx-blue.svg)](https://gcc.gnu.org/wiki/GFortran)
[![License](https://img.shields.io/badge/License-Academic_Research-green.svg)](https://www.csic.es)

**`MCcol`** is an advanced, high-performance atomistic Monte Carlo simulation suite designed for classical condensed-matter physics, materials science, and physical chemistry. The core simulation engine, **`gpMC`** (*General Purpose Monte Carlo*), simulates bulk multicomponent, molecular, and ionic systems with arbitrary stoichiometry in periodic orthorhombic domains.

---

## Authors & Affiliations

- **Enrique Lomba** ([enrique.lomba@csic.es](mailto:enrique.lomba@csic.es))
- **Eva G. Noya** ([eva.noya@iqf.csic.es](mailto:eva.noya@iqf.csic.es))

*Instituto de Química Física Blas Cabrera (formerly Rocasolano)*  
*Consejo Superior de Investigaciones Científicas (IQFR-CSIC)*  
*C/ Serrano 119, 28006 Madrid, Spain*

---

## Table of Contents

1. [Key Features & Theoretical Background](#key-features--theoretical-background)
   - [Thermodynamic Ensembles](#1-thermodynamic-ensembles)
   - [Interatomic Potentials](#2-interatomic-potentials)
   - [Long-Range Electrostatics (Ewald Summation)](#3-long-range-electrostatics-ewald-summation)
   - [Performance Optimizations](#4-performance-optimizations)
   - [Fault Tolerance & Emergency Checkpointing](#5-fault-tolerance--emergency-checkpointing)
   - [Interoperability (DL_POLY & LAMMPS)](#6-interoperability-dl_poly--lammps)
2. [Source Code Architecture](#source-code-architecture)
3. [Compilation & Building](#compilation--building)
   - [Prerequisites](#prerequisites)
   - [Building with GNU Fortran (gfortran)](#building-with-gnu-fortran-gfortran)
   - [Building with Intel Fortran (ifx / ifort)](#building-with-intel-fortran-ifx--ifort)
   - [Build Targets](#build-targets)
4. [Input File Specifications](#input-file-specifications)
   - [`system.dat` (System Topology & Interactions)](#systemdat-system-topology--interactions)
   - [`runMC.dat` (Run Parameters & Ensemble Control)](#runmcdat-run-parameters--ensemble-control)
   - [Coordinate Files (`CONFIG` and `data.atoms`)](#coordinate-files-config-and-dataatoms)
5. [Output Files & Data Analysis](#output-files--data-analysis)
6. [Examples Directory](#examples-directory)
7. [References](#references)

---

## Key Features & Theoretical Background

### 1. Thermodynamic Ensembles

- **Canonical Ensemble ($NVT$)**: Fixed atom count $N$, simulation box volume $V$, and temperature $T$. Single-particle trial moves are generated within an anisotropic displacement window $[-\Delta r_x, +\Delta r_x] \times [-\Delta r_y, +\Delta r_y] \times [-\Delta r_z, +\Delta r_z]$ and evaluated via the standard Metropolis acceptance criterion:
  $$\mathcal{P}_{\text{acc}}(o \to n) = \min\left(1, e^{-\beta \Delta U}\right), \quad \beta = \frac{1}{k_B T}$$

- **Isobaric-Isothermal Ensemble ($NpT$)**: Fixed particle count $N$, pressure $P$, and temperature $T$. In addition to particle displacements, the code performs trial volume changes with two scaling modes:
  - **Isotropic (`"isotr"`)**: Scales all three box vectors uniformly by $\alpha = (V'/V)^{1/3}$.
  - **Anisotropic / Orthorhombic (`"ortho"`)**: Allows independent length fluctuations along $x$, $y$, and $z$, enabling anisotropic stress relaxation.
  - Metropolis acceptance for an isotropic volume change:
    $$\mathcal{P}_{\text{acc}}(V \to V') = \min\left(1, \exp\left[-\beta \left(\Delta U + P(V'-V) - (N+1) k_B T \ln\frac{V'}{V}\right)\right]\right)$$

### 2. Interatomic Potentials

- **Morse Potential (`keyp = 1`)**:
  $$u_{ij}(r) = D_{e,ij} \left[ e^{-2\gamma_{ij}(r - r_{0,ij})} - 2 e^{-\gamma_{ij}(r - r_{0,ij})} \right]$$
  where $D_{e}$ is well depth, $\gamma$ governs potential stiffness, and $r_0$ is equilibrium bond separation.
- **Lennard-Jones 12-6 Potential (`keyp = 2`)**:
  $$u_{ij}(r) = 4\varepsilon_{ij} \left[ \left(\frac{\sigma_{ij}}{r}\right)^{12} - \left(\frac{\sigma_{ij}}{r}\right)^6 \right]$$
- **Potential Truncation & Shifting (`pshift = .true.`)**: Shifts potential energy to zero at cutoff $r_c$:
  $$u_{\text{shifted}}(r) = u(r) - u(r_c)$$
- **Energy Unit Systems**: Native conversion across `eV`, `K` ($k_B T$), and `kcal/mol`.

### 3. Long-Range Electrostatics (Ewald Summation)

Periodic boundary conditions with Coulombic charges require an Ewald decomposition into real-space, reciprocal-space, and self-interaction terms:

1. **Real-Space Screened Sum**:
   $$U_{\text{real}} = \frac{1}{2} \sum_{i=1}^N \sum_{j \neq i}^N \frac{q_i q_j \text{erfc}(\kappa r_{ij})}{4\pi\varepsilon_0 r_{ij}}$$
   where $\kappa$ is the Gaussian charge screening parameter.
2. **Reciprocal-Space Fourier Sum**:
   $$U_{\text{Fourier}} = \frac{1}{2\varepsilon_0 V} \sum_{\mathbf{k} > 0} \frac{\exp\left(-\frac{k^2}{4\kappa^2}\right)}{k^2} |\rho(\mathbf{k})|^2$$
   evaluated using half-space symmetry ($\mathbf{k} > 0$) with precomputed trilinear factor tables ($e^{i k_x x}$, $e^{i k_y y}$, $e^{i k_z z}$) for maximum efficiency.
3. **Electrostatic Self-Energy**:
   $$U_{\text{self}} = -\frac{\kappa}{4\pi^{3/2}\varepsilon_0} \sum_{i=1}^N q_i^2$$
4. **$O(K)$ Incremental Reciprocal Updates**:
   When displacing a single particle $i \to i'$, the charge density structure factor $\rho(\mathbf{k}) = \sum_{j} q_j e^{i \mathbf{k}\cdot\mathbf{r}_j}$ updates in $O(K)$ operations instead of recalculating over all $N$ atoms:
   $$\Delta \rho(\mathbf{k}) = q_i \left( e^{i \mathbf{k}\cdot\mathbf{r}_i'} - e^{i \mathbf{k}\cdot\mathbf{r}_i} \right)$$

### 4. Performance Optimizations

- **3D Link-Cell Domain Decomposition**: Subdivides the simulation cell into a 3D grid of size $\ge r_c + \Delta r_{\max}$. Neighbor evaluations check only the 27 neighboring subcells, reducing pairwise complexity from $O(N^2)$ to $O(N)$ for large systems.
- **Paul Breeuwsma Cubic Spline Interpolation**: Precalculates interatomic pair potentials on a dense radial lookup grid with smooth $C^2$-continuous 4-point cubic spline interpolation, eliminating expensive transcendental operations (`exp`, `sqrt`) during Monte Carlo loops.

### 5. Fault Tolerance & Emergency Checkpointing

- **Signal Interception**: A POSIX signal handler (`src/libutil.c`) intercepts termination signals (`SIGTERM`, `SIGINT`).
- **Cluster Walltime Protection**: When terminated by cluster resource managers (e.g., SLURM, PBS, Grid Engine), `gpMC` catches the signal, performs an orderly serialization of all coordinates, RNG states, metric tensors, and accumulators to `dump<TIMESTAMP>.dmp`, and exits safely without corrupting data.
- **Hot-Restart (`restart = .true.`)**: Resumes from `restart.dmp` seamlessly, setting output files to append mode.

### 6. Interoperability (DL_POLY & LAMMPS)

- **Input Formats**:
  - `initcf = "dlp"`: DL_POLY 2 coordinate file (`CONFIG`).
  - `initcf = "lmp"`: LAMMPS data file (`data.atoms`), auto-detecting `full`, `molecular`, and `atomic` atom styles.
- **Output Formats**:
  - `CONFIG.last`: DL_POLY 2 formatted final structure.
  - `last.conf`: LAMMPS data file formatted final structure.
  - `gpMC.lammpstrj`: LAMMPS custom trajectory dump, directly visualizable in **VMD** and **OVITO**.

---

## Source Code Architecture

All source code resides in [`src/`](./src/). Every module, subroutine, and interface includes complete scientific and algorithmic documentation:

| Source File | Primary Purpose | Key Subroutines / Data |
| :--- | :--- | :--- |
| [`src/set_precision.f90`](file:///home/elomba/MCcol/src/set_precision.f90) | Working precision constants | Single (`skind`), double (`dkind`, `wp = selected_real_kind(15, 307)`) precision. |
| [`src/Definitions.f90`](file:///home/elomba/MCcol/src/Definitions.f90) | Global state & module declarations | `configuration`, `potential`, `rundata`, `properties`, `linkcell`, `interp`, `interfaces`. |
| [`src/Distance.f90`](file:///home/elomba/MCcol/src/Distance.f90) | Distance calculations & MIC | `dist2`: Minimum image convention, coordinate unscaling, squared distance $r_{ij}^2$. |
| [`src/Interactions.f90`](file:///home/elomba/MCcol/src/Interactions.f90) | Pair interaction kernels | `fpot_Morse`, `fpot_LJ`, `fpot_elecMorse`, `fpot_elecLJ`, potential shifting (`ucut`). |
| [`src/Cells.f90`](file:///home/elomba/MCcol/src/Cells.f90) | 3D link-cell decomposition | `Init_cell`, `build_cells`, `update_cell_list`: 27-neighbor list stencil and $O(1)$ linked lists. |
| [`src/Energy.f90`](file:///home/elomba/MCcol/src/Energy.f90) | System energy calculation | `energ` (all-pairs), `energ_cell` (link-cell), `fourier` (reciprocal Ewald sum), `Eshort_r`. |
| [`src/Move.f90`](file:///home/elomba/MCcol/src/Move.f90) | Metropolis trial displacements | `move_natoms`: link-cell trial moves, cubic spline energy delta, $O(K)$ Fourier update. |
| [`src/VolumeMove.f90`](file:///home/elomba/MCcol/src/VolumeMove.f90) | $NpT$ volume trial moves | `move_volume`: isotropic (`"isotr"`) and anisotropic (`"ortho"`) box rescaling. |
| [`src/Init.f90`](file:///home/elomba/MCcol/src/Init.f90) | System initialization | `Init_conf`, `Init_rundata`, `read_potpars`, `Init_pot`, `Init_selfe`, `Init_fourier`, `Init_interp`. |
| [`src/Readconf.f90`](file:///home/elomba/MCcol/src/Readconf.f90) | Structure input parsing | `dlplmp_readconf`: DL_POLY `CONFIG` and LAMMPS `data.atoms` readers. |
| [`src/Structure.f90`](file:///home/elomba/MCcol/src/Structure.f90) | Structural correlations | `gr`: Multicomponent partial radial distribution functions $g_{ij}(r)$ with volume shell normalization. |
| [`src/Thermo.f90`](file:///home/elomba/MCcol/src/Thermo.f90) | Thermodynamic statistics | `averages`: Block accumulators, instantaneous sampling, equilibration tracking. |
| [`src/Dump.f90`](file:///home/elomba/MCcol/src/Dump.f90) | State serialization & restarts | `cierra`: Unformatted binary stream dump; `load`: Seamless restart restoration. |
| [`src/Output.f90`](file:///home/elomba/MCcol/src/Output.f90) | Formatted logging & tables | `initout`, `run_info`, `printout`, `printgr`, `end_printout`. |
| [`src/WriteCfg.f90`](file:///home/elomba/MCcol/src/WriteCfg.f90) | Configuration and trajectory output | `writecfg_dlp`, `writecfg_lmp`, `dump_trj` (LAMMPS custom trajectory format). |
| [`src/Util.f90`](file:///home/elomba/MCcol/src/Util.f90) | Timing utilities | `cputime`: High-precision POSIX CPU time accounting. |
| [`src/libutil.c`](file:///home/elomba/MCcol/src/libutil.c) | C system-level interop | `catch_`: Signal handler (`SIGTERM`, `SIGINT`); `cputime_`: POSIX `times()`; `gethost_`: Hostname lookup. |
| [`src/Main.f90`](file:///home/elomba/MCcol/src/Main.f90) | Main program driver | `gpMC`: Simulation lifecycle, phase transitions, main sweep loop. |

---

## Compilation & Building

### Prerequisites
- **Fortran Compiler**: GNU Fortran (`gfortran` $\ge$ 9.0) or Intel Fortran (`ifx` / `ifort` $\ge$ 13.0)
- **C Compiler**: `gcc` or `icx` / `icc`
- **Build System**: GNU `make`

### Building with GNU Fortran (gfortran)
```bash
# Clean and compile optimized executable (bin/gpMC.exe)
make -f Makefile.gfortran clean all

# Or compile with full debugging symbols and runtime bounds checking
make -f Makefile.gfortran clean debug
```

### Building with Intel Fortran (ifx / ifort)
```bash
# Clean and compile optimized executable
make clean all

# Or compile with debugging symbols
make clean debug
```

### Build Targets
- `make all`: Compiles the binary to `bin/gpMC.exe` and creates the `results/` directory.
- `make clean`: Removes all compiled object files (`obj/*.o`) and module definitions (`modules/*.mod`).
- `make debug`: Compiles with `-O0 -g -fcheck=all` (gfortran) or `-debug -O0 -CB` (Intel).

---

## Input File Specifications

Every simulation requires two configuration files in the working directory: `system.dat` and `runMC.dat`, along with the initial structure file (`CONFIG` or `data.atoms`).

### `system.dat` (System Topology & Interactions)

Line-by-line syntax for `system.dat`:

```text
.false.                 ! Line 1:  restart (.true. to load restart.dmp, .false. for new run)
"dlp"                   ! Line 2:  initcf: initial configuration format ('dlp' or 'lmp')
2 4608                  ! Line 3:  nsp (number of species), natoms (total atom count)
1 "O"  -0.65            ! Line 4:  species 1 ID, label (in quotes), partial charge (e)
2 "Si"  1.3             ! Line 5:  species 2 ID, label, partial charge
"eV"                    ! Line 6:  units ('eV', 'K', or 'kcal/mol')
1 1 "mors"              ! Line 7:  pair (i, j) and potential type ('mors' or 'lj')
0.023272 1.3731 3.791 9.0 ! Line 8:  Morse parameters: De (eV), gamma (1/A), r0 (A), cutoff rc (A)
1 2 "morse"             ! Line 9:  pair (1, 2)
1.99597  2.6518 1.628 9.0 ! Line 10: Morse parameters for O-Si
2 2 "morse"             ! Line 11: pair (2, 2)
0.007695 2.0446 3.7598 9.0! Line 12: Morse parameters for Si-Si
.true. .false. 1        ! Line 13: elect (.true./.false.), pshift (.true./.false.), fou_type (1 = standard Ewald)
0.287184 15 15 15       ! Line 14: kappa (1/A), kx_max, ky_max, kz_max (Ewald parameters)
```

> **Interaction Parameters Note**:
> - For **Morse** (`"mors"` or `"morse"`): Specify `D_e  gamma  r_0  r_c`.
> - For **Lennard-Jones** (`"lj"`): Specify `epsilon  sigma  r_c`.
> - Interactions must be supplied for all unique pairs in upper triangular order: $(1,1), (1,2), \dots, (1,S), (2,2), \dots, (S,S)$.

---

### `runMC.dat` (Run Parameters & Ensemble Control)

#### Canonical Ensemble ($NVT$)
```text
"nvt"                   ! Line 1: ensemble ('nvt')
100 20 10 20 20         ! Line 2: nstep, nequil, nb (averaging block), npgr (g(r) dump), ntraj (trajectory dump)
0.05 0.05 0.05          ! Line 3: rdmax(1:3): max trial displacement along x, y, z (fractional coords)
4000.0                  ! Line 4: temp: temperature (Kelvin)
0.05                    ! Line 5: deltagr: histogram radial bin width (Angstroms)
12345 67890 54321 98765 11223 44556 77889 99001  ! Line 6: random number seed (8 integers for gfortran)
```

#### Isobaric-Isothermal Ensemble ($NpT$)
```text
"npt"                   ! Line 1: ensemble ('npt')
100 20 10 20 20         ! Line 2: nstep, nequil, nb, npgr, ntraj
0.05 0.05 0.05 0.10 "isotr" ! Line 3: rdmax(1:3), vdmax (max volume delta), scaling ('isotr' or 'ortho')
4000.0 166.054          ! Line 4: temp (Kelvin), pres (bars)
0.05                    ! Line 5: deltagr (Angstroms)
12345 67890 54321 98765 11223 44556 77889 99001  ! Line 6: random number seed
```

---

### Coordinate Files (`CONFIG` and `data.atoms`)

- **DL_POLY 2 `CONFIG`**:
  - Line 1: Title / description.
  - Line 2: `levcfg` (0, 1, or 2), `imcon` (boundary condition, 1-3 for periodic orthorhombic).
  - Lines 3–5: Lattice box vectors $\mathbf{a}, \mathbf{b}, \mathbf{c}$.
  - Atom records: species name, index, followed by Cartesian coordinates $x, y, z$.
- **LAMMPS `data.atoms`**:
  - Standard LAMMPS data format with headers specifying atom counts, atom types, and orthogonal box bounds (`xlo xhi`, `ylo yhi`, `zlo zhi`).
  - Supported atom styles: `full` (7 columns: `id mol type q x y z`), `molecular` (6 columns), and `atomic` (5 columns).

---

## Output Files & Data Analysis

| Output File | Format | Description |
| :--- | :--- | :--- |
| **`thermoins.dat`** | ASCII table | Instantaneous values recorded at every sweep: step, acceptance %, $E_{\text{tot}}$, $E_{\text{sr}}$, $E_{\text{vdw}}$, $E_{\text{fourier}}$, $E_{\text{self}}$, $E_{\text{coul}}$, volume $V$. |
| **`thermoaver.dat`** | ASCII table | Block-averaged thermodynamic values computed every `nb` sweeps with cumulative averages and statistical errors. |
| **`gmix.dat`** | Multi-column ASCII | Partial pair distribution functions $g_{ij}(r)$ for all species pairs across radial distance $r$. |
| **`CONFIG.last`** | DL_POLY 2 format | Final atomic coordinates and lattice vectors (generated when `initcf = "dlp"`). |
| **`last.conf`** | LAMMPS data format | Final atomic configuration in LAMMPS data format (generated when `initcf = "lmp"`). |
| **`gpMC.lammpstrj`** | LAMMPS custom dump | Molecular trajectory snapshots recorded every `ntraj` sweeps, visualizable in VMD or OVITO. |
| **`dump<YYYYMMDDHHMM>.dmp`** | Unformatted binary stream | Checkpoint file containing complete microscopic and statistical state for restart operations. |

---

## Examples Directory

A dedicated [`examples/`](./examples/) directory provides fully configured, self-contained test cases:

- **[`01_nvt_silica_morse_dlpoly/`](./examples/01_nvt_silica_morse_dlpoly/)**: Canonical ($NVT$) simulation of silicalite zeolite (4608 atoms) with Morse potentials and Ewald electrostatics using a DL_POLY `CONFIG` file.
- **[`02_npt_silica_morse_dlpoly/`](./examples/02_npt_silica_morse_dlpoly/)**: Isobaric-Isothermal ($NpT$) simulation of silicalite with isotropic box volume moves.
- **[`03_nvt_silica_lammps_input/`](./examples/03_nvt_silica_lammps_input/)**: Canonical ($NVT$) simulation initialized from LAMMPS `data.atoms`.
- **[`04_restart_checkpoint/`](./examples/04_restart_checkpoint/)**: Two-phase simulation demonstrating binary checkpoint generation and hot restart.
- **[`05_lammps_md_comparison/`](./examples/05_lammps_md_comparison/)**: Cross-validation benchmark scripts and reference profiles ($g(r)$, MSD, thermodynamics) comparing Monte Carlo with Molecular Dynamics (LAMMPS).

To run any example:
```bash
cd examples/01_nvt_silica_morse_dlpoly
./run.sh
```

---

## References

1. **Allen, M. P., & Tildesley, D. J.** (2017). *Computer Simulation of Liquids* (2nd ed.). Oxford University Press.
2. **Frenkel, D., & Smit, B.** (2002). *Understanding Molecular Simulation: From Algorithms to Applications* (2nd ed.). Academic Press.
3. **Breeuwsma, P.** (2003). *Cubic Spline Interpolation Techniques in Physical Computations*.
4. **Smith, W., & Forester, T. R.** (1996). *DL_POLY_2.0: A general-purpose parallel molecular dynamics simulation package*. Journal of Molecular Graphics, 14(3), 136-141.
5. **Thompson, A. P. et al.** (2022). *LAMMPS - a flexible simulation tool for particle-based materials modeling at the atomic, meso, and continuum scales*. Computer Physics Communications, 271, 108171.
