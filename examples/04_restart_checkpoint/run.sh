#!/bin/bash
# Runner script demonstrating checkpoint save and restart in gpMC
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
echo " Step 1: Running Phase 1 (Fresh start, 20 sweeps)         "
echo "=========================================================="
cp system.dat.initial system.dat
cat << 'PHASE1' > runMC.dat
"nvt"
20 5 5 10 10
0.05 0.05 0.05
4000.0
0.05
12345 67890 54321 98765 11223 44556 77889 99001
PHASE1

"$BIN"

echo "Phase 1 completed. Locating binary checkpoint dump..."
LATEST_DMP=$(ls -t dump*.dmp 2>/dev/null | head -n 1)
if [ -z "$LATEST_DMP" ]; then
    echo "Error: No dump file produced!"
    exit 1
fi
echo "Found checkpoint dump: $LATEST_DMP"

echo "=========================================================="
echo " Step 2: Preparing for Phase 2 Restart                   "
echo "=========================================================="
cp "$LATEST_DMP" restart.dmp
cp system.dat.restart system.dat

cat << 'PHASE2' > runMC.dat
"nvt"
40 5 5 10 10
0.05 0.05 0.05
4000.0
0.05
12345 67890 54321 98765 11223 44556 77889 99001
PHASE2

echo "=========================================================="
echo " Step 3: Running Phase 2 (Resuming from sweep 20 to 40)   "
echo "=========================================================="
"$BIN"

echo "=========================================================="
echo " Restart test completed successfully!                     "
echo "=========================================================="
