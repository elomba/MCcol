# Example 02: NpT Simulation of Silicalite Zeolite (Isotropic Volume Fluctuations)

## Overview
This example demonstrates an Isobaric-Isothermal ($NpT$) Monte Carlo simulation of silicalite zeolite ($\text{SiO}_2$, 4,608 atoms) at high temperature ($T = 4000\ \text{K}$) and pressure ($P = 166.054\ \text{bar} = 0.166\ \text{katm}$).

In the $NpT$ ensemble, trial moves comprise:
1. **Single-particle trial displacements** distributed uniformly in $[-\Delta r_x, +\Delta r_x] \times [-\Delta r_y, +\Delta r_y] \times [-\Delta r_z, +\Delta r_z]$.
2. **Volume trial moves** attempting isotropic box rescaling ($\ln V \to \ln V + \Delta$), updating cell dimensions $L_x, L_y, L_z$, atom coordinates, and reciprocal-space wavevectors.

The acceptance probability for an isotropic volume change $V \to V'$ is governed by the standard Metropolis criterion:
$$\mathcal{P}_{\text{acc}} = \min\left(1, \exp\left[-\beta \left(\Delta U + P(V'-V) - (N+1) k_B T \ln(V'/V)\right)\right]\right)$$

---

## File Contents
- `CONFIG`: Initial DL_POLY 2 configuration file (4,608 atoms).
- `system.dat`: Morse potential parameters, species charges, and Ewald parameters.
- `runMC.dat`: $NpT$ ensemble parameters, target pressure, maximum volume displacement $\Delta V_{\max}$, and scaling mode (`"isotr"`).
- `run.sh`: Automated execution script.

---

## Input File Breakdown

### `runMC.dat`
```text
"npt"                           ! simulation ensemble ('npt')
50 10 10 10 10                  ! nstep, nequil, nb, npgr, ntraj
0.05 0.05 0.05 0.10 "isotr"     ! rdmax(1:3), vdmax, scaling mode ('isotr' or 'ortho')
4000.0 166.054                  ! temp (K), pres (bar)
0.05                            ! deltagr (A)
12345 67890 54321 98765 11223 44556 77889 99001  ! random seed
```

> **Note on Scaling Mode:**
> - `"isotr"`: Scales all box dimensions proportionally ($L_x, L_y, L_z$ scaled by $\alpha = (V'/V)^{1/3}$), preserving the cubic/orthorhombic aspect ratio.
> - `"ortho"`: Allows independent fluctuations of each box dimension, suitable for systems with anisotropic compressibility or stress relaxation.

---

## How to Run
```bash
# Using runner script
./run.sh

# Or direct execution
../../bin/gpMC.exe
```

---

## Expected Output Files
| Output File | Description |
| :--- | :--- |
| `thermoins.dat` | Instantaneous total energy, Coulomb components, and instantaneous simulation box volume $V(t)$ at each sweep. |
| `thermoaver.dat` | Block-averaged energy, volume $\langle V \rangle$, and average box lengths $\langle L_x \rangle, \langle L_y \rangle, \langle L_z \rangle$. |
| `gmix.dat` | Partial pair correlation functions $g_{\text{O-O}}(r)$, $g_{\text{O-Si}}(r)$, and $g_{\text{Si-Si}}(r)$. |
| `CONFIG.last` | Final configuration in DL_POLY 2 format with updated box vector headers. |
| `gpMC.lammpstrj` | Trajectory file with updated box boundaries recorded at `ntraj` intervals. |
| `dump<YYYYMMDDHHMM>.dmp` | Binary restart checkpoint. |
