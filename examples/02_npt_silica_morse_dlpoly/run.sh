#!/bin/bash
# Runner script for Example 02: NpT Silicalite (Morse + Ewald, Isotropic Volume Moves)
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

BIN="../../bin/gpMC.exe"

if [ ! -f "$BIN" ]; then
    echo "Error: $BIN executable not found! Please build gpMC first:"
    echo "  make -f Makefile.gfortran clean all"
    exit 1
fi

echo "=========================================================="
echo " Running Example 02: Silicalite NpT Simulation (gpMC)     "
echo "=========================================================="
"$BIN"
echo "=========================================================="
echo " Example 02 completed successfully!                       "
echo " Generated outputs:"
echo "   - thermoins.dat    : Instantaneous thermodynamic properties & volume"
echo "   - thermoaver.dat   : Block averages of energy, volume, and density"
echo "   - gmix.dat         : Partial pair correlation functions g_ij(r)"
echo "   - CONFIG.last      : Final atom coordinates & updated box vectors"
echo "   - gpMC.lammpstrj   : LAMMPS trajectory file"
echo "   - dump*.dmp        : Binary checkpoint restart file"
echo "=========================================================="
