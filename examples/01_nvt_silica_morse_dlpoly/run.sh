#!/bin/bash
# Runner script for Example 01: NVT Silicalite (Morse + Ewald)
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
echo " Running Example 01: Silicalite NVT Simulation (gpMC)     "
echo "=========================================================="
"$BIN"
echo "=========================================================="
echo " Example 01 completed successfully!                       "
echo " Generated outputs:"
echo "   - thermoins.dat    : Instantaneous thermodynamic properties"
echo "   - thermoaver.dat   : Block averages of energy and pressure"
echo "   - gmix.dat         : Partial pair correlation functions g_ij(r)"
echo "   - CONFIG.last      : Final atom coordinates (DL_POLY format)"
echo "   - gpMC.lammpstrj   : LAMMPS trajectory file"
echo "   - dump*.dmp        : Binary checkpoint restart file"
echo "=========================================================="
