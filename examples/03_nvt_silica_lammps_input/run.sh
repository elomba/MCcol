#!/bin/bash
# Runner script for Example 03: NVT Silicalite from LAMMPS data.atoms
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
echo " Running Example 03: LAMMPS data.atoms NVT Simulation     "
echo "=========================================================="
"$BIN"
echo "=========================================================="
echo " Example 03 completed successfully!                       "
echo " Generated outputs:"
echo "   - thermoins.dat    : Instantaneous thermodynamic properties"
echo "   - thermoaver.dat   : Block averages of energy and pressure"
echo "   - gmix.dat         : Partial pair correlation functions g_ij(r)"
echo "   - last.conf        : Final configuration in LAMMPS data format"
echo "   - gpMC.lammpstrj   : LAMMPS trajectory file"
echo "   - dump*.dmp        : Binary checkpoint restart file"
echo "=========================================================="
