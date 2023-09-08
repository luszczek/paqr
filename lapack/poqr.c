#include <stdlib.h>
#include <stdio.h>
#include <time.h>
//#include "mex.h"
#include "matrix.h"

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))

/* Interface to the Fortran routine */
void dgepoqrf_
(
    int const* m,
    int const* n,
    double* A,
    int const* lda,
    double* tau,
    int* nd,
    int* deff,
    double* work,
    int const* lwork,
    int* info
);

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* variable declarations here */

    int m, n, lda, lwork, info;
    double *A, *RV, *tau, *work;
    int nd, *deff;
    struct timespec start, end;
    double runtime;

    /* code here */

    /* get dimensions of the input matrix */
    m = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);

    lda = m; //XXX
    lwork = min(m,n) * min(m,n);
    info = 0;

    nd = 0;

    /* create a pointer to the real data in the input matrix  */
    A = mxGetPr(prhs[0]);

    /* create the output matrix */
    plhs[0] = mxDuplicateArray(prhs[0]);
    plhs[1] = mxCreateNumericMatrix(1, min(m,n), mxDOUBLE_CLASS, mxREAL);
    plhs[2] = mxCreateNumericMatrix(1, min(m,n), mxINT32_CLASS, mxREAL);

    /* get a pointer to the real data in the output matrix */
    RV = mxGetPr(plhs[0]);
    tau = mxGetPr(plhs[1]);
    deff = mxGetPr(plhs[2]);

    work = (double*) malloc(lwork * sizeof(double));

    /* call the computational routine */
    clock_gettime(CLOCK_REALTIME, &start);
    dgepoqrf_
    (
        &m,
        &n,
        RV,
        &lda,
        tau,
        &nd,
        deff,
        work,
        &lwork,
        &info
    );
    clock_gettime(CLOCK_REALTIME, &end);
    runtime = (end.tv_sec - start.tv_sec) * 1e9;
    runtime = (runtime + (end.tv_nsec - start.tv_nsec)) * 1e-9;
    printf("Elapsed time is %.9e seconds.\n", runtime);

    free(work);

    return;
}

