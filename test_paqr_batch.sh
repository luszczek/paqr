#!/bin/bash

############################################
# change test parameters
paqr_tolerance=1       # specify multiples of unit roundoff
threads_27x20=27       # choose a value between 20 and 27
threads_125x56=96      # choose a value between 56 and 125
batch=1000             # number of matrices
iter=10                # repititions of each test
check=-c               # keep it for error checking
############################################

source /home/.bashrc
cd /home/paqr-gpu/testing

echo "###########################################"
echo "#               TESTING qr_gpu            #"
echo "# INFO: performance is under MAGMA column #"
echo "###########################################"
./testing_dpaqrf_batched --align 1 --niter ${iter} --batch ${batch} -N 1 -c --tol ${paqr_tolerance} --nb ${threads_27x20}  --version 3
./testing_dpaqrf_batched --align 1 --niter ${iter} --batch ${batch} -N 2 -c --tol ${paqr_tolerance} --nb ${threads_125x56} --version 3

echo "###########################################"
echo "#               TESTING paqr_gpu          #"
echo "# INFO: performance is under MAGMA column #"
echo "###########################################"
./testing_dpaqrf_batched --align 1 --niter ${iter} --batch ${batch} -N 1 -c --tol ${paqr_tolerance} --nb ${threads_27x20}  --version 1
./testing_dpaqrf_batched --align 1 --niter ${iter} --batch ${batch} -N 2 -c --tol ${paqr_tolerance} --nb ${threads_125x56} --version 1

cd /home
