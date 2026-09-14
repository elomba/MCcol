# Example 05: Comparison with Molecular Dynamics (LAMMPS)

## Overview
This directory provides reference input scripts and benchmark simulation profiles from **LAMMPS** (Large-scale Atomic/Molecular Massively Parallel Simulator) that serve as cross-validation benchmarks for `gpMC`.

Because Monte Carlo ($MC$) and Molecular Dynamics ($MD$) sample the same equilibrium thermodynamic ensembles (e.g., canonical $NVT$ or isobaric-isothermal $NpT$), macroscopic equilibrium averages and structural correlation functions computed by `gpMC` must agree within statistical uncertainty with those produced by LAMMPS under identical thermodynamic state points and Hamiltonian parameters.

---

## File Contents
| File | Description |
| :--- | :--- |
| `data.atoms` | Initial atomic configuration (orthorhombic simulation box with 4,608 atoms). |
| `morse.lmp` | LAMMPS Molecular Dynamics input script for Morse potential + PPPM long-range electrostatics. |
| `lj.lmp` | LAMMPS input script for Lennard-Jones (12-6) + PPPM electrostatics. |
| `forcefield.morse` | Pair style definitions and potential parameters for the Morse force field. |
| `forcefield.lj` | Pair style definitions and potential parameters for the Lennard-Jones force field. |
| `gr_P_1_T_400.rdf` | Reference radial distribution function $g(r)$ calculated via LAMMPS MD at $T = 400\ \text{K}$ and $P = 1\ \text{bar}$. |
| `msd_P_1_T_400.profile` | Reference Mean Squared Displacement (MSD) profile. |
| `outvars_P_1_T_400.profile` | Reference thermodynamic output trajectory (temperature, pressure, potential energy, density). |

---

## Comparing MC (`gpMC`) with MD (LAMMPS)

### 1. Radial Distribution Functions $g_{ij}(r)$
- In `gpMC`, partial pair distribution functions are accumulated across simulation sweeps and written to `gmix.dat` (column 1 is distance $r$ in Å, subsequent columns are $g_{ij}(r)$).
- In LAMMPS, `fix ave/time` with `compute rdf` writes the reference curve to `gr_P_1_T_400.rdf`.
- Overlaying column 2 ($g_{\text{O-O}}$) and column 3 ($g_{\text{O-Si}}$) from `gmix.dat` onto `gr_P_1_T_400.rdf` directly confirms that the equilibrium structural coordination and peak positions match.

### 2. Energetics and Pressure
- Compare the block-averaged potential energy per atom $\langle E_{\text{tot}} \rangle / N$ from `thermoaver.dat` in `gpMC` with the average `PotEng/atoms` in `outvars_P_1_T_400.profile`.
- Ensure energy units match: `gpMC` supports `eV`, `K`, and `kcal/mol`. In LAMMPS `units real`, energy is in $\text{kcal/mol}$; in `units metal`, energy is in $\text{eV}$.

### 3. Running LAMMPS (if installed)
```bash
lmp -in morse.lmp
# or
lmp -in lj.lmp
```
