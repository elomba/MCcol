# Example 04: Checkpoint and Restart Workflow

## Overview
This example illustrates `gpMC`'s state serialization and hot-restart capabilities.
During a simulation, `gpMC` automatically writes an unformatted binary stream checkpoint named `dump<YYYYMMDDHHMM>.dmp` upon normal completion or when receiving an external termination signal (`SIGTERM` or `SIGINT`).

The checkpoint records the complete simulation state:
- Exact random number generator state (reproducible Markov chain).
- Scaled atomic positions $\mathbf{s}_i$ and partial charges $q_i$.
- Unit cell dimensions and metric tensors ($\mathbf{a}, \mathbf{b}, \mathbf{c}$, volume $V$).
- Precomputed reciprocal-space wavevectors and dynamic structure factors ($\rho(\mathbf{k})$, $e^{i \mathbf{k}\cdot\mathbf{r}_i}$).
- Cubic spline potential lookup tables ($u(r)$, $r_{\min}^2$).
- Thermodynamic accumulators and block averages ($\langle E_{\text{tot}} \rangle, \langle V \rangle$).
- Partial pair distribution function histograms ($h_{ij}(r)$).

---

## Restart Protocol

### 1. Identify the Dump File
Locate the generated binary dump file (e.g., `dump202609140940.dmp`) and copy or link it to `restart.dmp`:
```bash
cp dump202609140940.dmp restart.dmp
```

### 2. Enable Restart in `system.dat`
Set the first line of `system.dat` from `.false.` to `.true.`:
```text
.true.                  ! restart flag (.true. tells gpMC to load restart.dmp)
"dlp"
...
```

### 3. Update Target Sweeps in `runMC.dat`
Adjust `nstep` in `runMC.dat` to the new target total sweep number (e.g., from 20 to 40).
When `gpMC` restarts, it loads the step counter from `restart.dmp` and proceeds seamlessly up to the new `nstep`.

Outputs (`thermoins.dat`, `thermoaver.dat`) are automatically opened in **append mode**, preserving historical trajectory data.

---

## Automated Demonstration
Run the provided runner script to execute both phases sequentially:
```bash
./run.sh
```

---

## Signal Handling and Emergency Checkpointing
`gpMC` installs a POSIX signal handler via C interop (`libutil.c`). If a simulation running under a cluster queue manager (SLURM, PBS) is about to exceed its walltime and receives `SIGTERM`, `gpMC` intercepts the signal, performs an orderly unformatted stream dump, and exits gracefully, preventing data loss.
