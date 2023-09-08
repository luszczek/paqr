#make blaslib 
make lapacklib
#gcc -g -o test_dgepoqrf test_dgepoqrf.c -lm -lgfortran -L . -llapack -lrefblas
gcc -g -o test_dgepoqrf test_dgepoqrf.c -lm -lgfortran -L . -llapack -L$MKLROOT/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -lm -lgomp
echo 'Afull'
./test_dgepoqrf $1 $2 1
echo 'Abeg'
./test_dgepoqrf $1 $2 2
echo 'Amid'
./test_dgepoqrf $1 $2 4
echo 'Aend'
./test_dgepoqrf $1 $2 5
#mex poqr.c -R2018a -lm -lgfortran -L. -llapack -lrefblas
mex poqr.c -R2018a -lm -lgfortran -L. -llapack -L$MKLROOT/lib/intel64 -lmkl_gf_lp64
