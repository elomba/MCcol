# Example 03: Canonical (NVT) Simulation from LAMMPS Data File

## Overview
This example demonstrates how to initialize a Monte Carlo simulation directly from a standard LAMMPS data file (`data.atoms`).

`gpMC` features a built-in parser for LAMMPS structure files supporting:
- Header sections: `atoms`, `atom types`, and orthogonal simulation box bounds `xlo xhi`, `ylo yhi`, `zlo zhi`.
- Atom data sections with flexible column formats:
  - **LAMMPS `full` style**: `id molecule-tag atom-type q x y z` (7 columns)
  - **LAMMPS `molecular` style**: `id molecule-tag atom-type x y z` (6 columns)
  - **LAMMPS `atomic` style**: `id atom-type x y z` (5 columns)

When initialized with `initcf = "lmp"`, the simulation reads `data.atoms` and produces final configuration dumps in LAMMPS data format (`last.conf`) as well as trajectory snapshots in LAMMPS custom dump format (`gpMC.lammpstrj`).

---

## File Contents
- `data.atoms`: Initial atomic structure in LAMMPS data format (4,608 atoms).
- `system.dat`: System configuration with `initcf = "lmp"`.
- `runMC.dat`: $NVT$ simulation parameters.
- `run.sh`: Automated execution script.

---

## Input File Breakdown

### `system.dat`
```text
.false.                 ! restart flag
"lmp"                   ! initcf: 'lmp' specifies reading from 'data.atoms'
2 4608                  ! 2 species (O and Si), 4608 atoms
1 "O"  -0.65            ! species 1
2 "Si"  1.3             ! species 2
"eV"                    ! energy units
1 1 "mors"              ! O-O Morse interaction
0.023272  1.3731  3.791  9.0
1 2 "morse"             ! O-Si Morse interaction
1.99597   2.6518  1.628  9.0
2 2 "morse"             ! Si-Si Morse interaction
0.007695  2.0446  3.7598 9.0
.true.  .false.  1      ! electrostatics enabled (Ewald)
0.287184 15 15 15       ! Ewald parameters
```

---

## How to Run
```bash
./run.sh
# or
../../bin/gpMC.exe
```

---

## Expected Output Files
| Output File | Description |
| :--- | :--- |
| `thermoins.dat` | Instantaneous energetics and acceptance ratios. |
| `thermoaver.dat` | Block-averaged thermodynamic observables. |
| `gmix.dat` | Partial pair distribution functions $g_{\text{O-O}}$, $g_{\text{O-Si}}$, $g_{\text{Si-Si}}$. |
| `last.conf` | Final atom coordinates written in standard LAMMPS data file format. |
| `gpMC.lammpstrj` | Atom trajectory snapshots in LAMMPS custom format (loadable into VMD or OVITO). |
| `dump<YYYYMMDDHHMM>.dmp` | Binary checkpoint file. |
