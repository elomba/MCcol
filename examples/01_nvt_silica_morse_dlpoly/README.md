# Example 01: NVT Simulation of Silicalite Zeolite (Morse + Ewald)

## Overview
This example demonstrates a Canonical ensemble ($NVT$) Monte Carlo simulation of a bulk silicalite zeolite ($\text{SiO}_2$) crystal containing 4,608 atoms (3,072 Oxygen atoms and 1,536 Silicon atoms) in an orthorhombic unit cell of dimensions approximately $40.15 \times 39.72 \times 40.54\ \text{Å}^3$.

The interaction model combines:
- **Morse interatomic pair potentials** for bonded and short-range non-bonded interactions.
- **Full Ewald summation** for long-range Coulombic electrostatics with formal charges ($q_\text{O} = -0.65e$, $q_\text{Si} = +1.30e$).
- **3D Link-Cell domain decomposition** for accelerated $O(N)$ neighbor searching.
- **Paul Breeuwsma cubic spline interpolation** for rapid potential evaluations.

Initial atomic coordinates are supplied via a standard DL_POLY 2 `CONFIG` file.

---

## File Contents
- `CONFIG`: Initial DL_POLY 2 coordinate file (4,608 atoms).
- `system.dat`: System topology, species charges, Morse potential parameters, energy units (`eV`), and Ewald electrostatics settings.
- `runMC.dat`: Run parameters for the $NVT$ ensemble ($T = 4000\ \text{K}$, 100 steps, equilibration, averaging blocks, random seed).
- `run.sh`: Automated execution script.

---

## Input File Breakdown

### `system.dat`
```text
.false.                 ! restart: .false. for a fresh simulation from CONFIG
"dlp"                   ! initcf: 'dlp' indicates DL_POLY 2 format ('CONFIG')
2 4608                  ! nsp (number of species), natoms (total atom count)
1 "O"  -0.65            ! species index 1, symbol 'O', partial charge -0.65e
2 "Si"  1.3             ! species index 2, symbol 'Si', partial charge +1.30e
"eV"                    ! energy units: 'eV' (also supports 'K' and 'kcal/mol')
1 1 "mors"              ! interaction between species 1-1 (O-O): Morse
0.023272  1.3731  3.791  9.0   ! D (eV), gamma (1/A), r0 (A), cutoff rc (A)
1 2 "morse"             ! interaction between species 1-2 (O-Si): Morse
1.99597   2.6518  1.628  9.0   ! D (eV), gamma (1/A), r0 (A), rc (A)
2 2 "morse"             ! interaction between species 2-2 (Si-Si): Morse
0.007695  2.0446  3.7598 9.0   ! D (eV), gamma (1/A), r0 (A), rc (A)
.true.  .false.  1      ! elect (.true. enables Ewald), pshift (potential shifting), fou_type (1 = standard Ewald)
0.287184 15 15 15       ! kappa (screening parameter 1/A), kx, ky, kz (reciprocal wavevector cutoffs)
```

### `runMC.dat`
```text
"nvt"                   ! simulation ensemble ('nvt' or 'npt')
100 20 10 20 20         ! nstep, nequil, nb (averaging block), npgr (g(r) interval), ntraj (trajectory interval)
0.05 0.05 0.05          ! rdmax(1:3): max trial displacement along x, y, z (fractional coordinates)
4000.0                  ! temp: temperature in Kelvin
0.05                    ! deltagr: bin width for radial distribution functions (A)
12345 67890 54321 98765 11223 44556 77889 99001  ! random number generator seed (8 integers for gfortran)
```

---

## How to Run
Execute the runner script or call the binary directly:
```bash
# Option 1: Using the runner script
./run.sh

# Option 2: Direct execution
../../bin/gpMC.exe
```

---

## Expected Output Files
Upon successful completion, the following files are produced:
| Output File | Description |
| :--- | :--- |
| `thermoins.dat` | Instantaneous total energy, short-range energy, Coulomb energy, Fourier reciprocal energy, self-energy, and acceptance ratio recorded at each sweep. |
| `thermoaver.dat` | Block-averaged thermodynamic values across `nb` sweep intervals with accumulated statistics. |
| `gmix.dat` | Multicomponent partial radial distribution functions $g_{\text{O-O}}(r)$, $g_{\text{O-Si}}(r)$, and $g_{\text{Si-Si}}(r)$. |
| `CONFIG.last` | Final configuration in DL_POLY 2 format. |
| `gpMC.lammpstrj` | Atom trajectory snapshots in LAMMPS custom format (plottable in VMD/OVITO). |
| `dump<YYYYMMDDHHMM>.dmp` | Binary unformatted checkpoint stream for restart operations. |
