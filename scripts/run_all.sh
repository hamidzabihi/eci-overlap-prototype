#!/usr/bin/env bash
set -euo pipefail
: "${SHARC:?Set SHARC}" "${ECI_PROTO:?Set ECI_PROTO}"
cd "$ECI_PROTO/scripts"

echo "########## 00 check_env ##########"
./00_check_env.sh

echo "########## 01 build_interface ##########"
./01_build_interface.sh

echo "########## 02 setup_demo ##########"
./02_setup_demo.sh

echo "########## 03 run_demo ##########"
if command -v g16 >/dev/null 2>&1 || command -v gaussian >/dev/null 2>&1 || \
   command -v orca >/dev/null 2>&1; then
    ./03_run_demo.sh || echo "  demo run failed; continuing to verify"
else
    echo "  skip: no Gaussian/ORCA on PATH; demo will run on Tomislav's cluster"
fi

echo "########## 04 verify ##########"
./04_verify.sh

echo "########## 05 package ##########"
./05_package.sh

echo "########## all done ##########"
