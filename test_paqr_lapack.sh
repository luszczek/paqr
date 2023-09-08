#!/bin/bash
# Script that compares the execution time of QR, PAQR, and QRCP
# on matrices where the non-deficient columns are placed either at 
# the beginning, the middle, or the end

if [ $# -lt 2 ]; then
  echo "Usage $0 <nb_rows> <nb_columns>"
  exit 1
else
  nrow=$1
  ncol=$2
fi

# Case where A is full rank
echo 'Case where A is full rank'
OMP_NUM_THREADS=1 ./lapack/test_dgepaqrf $nrow $ncol 1

# Case where the non-deficient columns are at the beginning
echo 'Case where A has non-deficient columns are at the beginning'
OMP_NUM_THREADS=1 ./lapack/test_dgepaqrf $nrow $ncol 2

# Case where the non-deficient columns are in the middle
echo 'Case where A has non-deficient columns are in the middle'
OMP_NUM_THREADS=1 ./lapack/test_dgepaqrf $nrow $ncol 4

# Case where the non-deficient columns are at the end
echo 'Case where A has non-deficient columns are at the end'
OMP_NUM_THREADS=1 ./lapack/test_dgepaqrf $nrow $ncol 5
