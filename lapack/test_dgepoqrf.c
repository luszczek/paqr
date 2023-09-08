#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define min(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))
#define COND(type,n,j) (cond##type((n),(j)))

//#define LAPACK_dgepoqrf LAPACK_GLOBAL(dgepoqrf,DGEPOQRF)
void dgeqrf_
(
    int const* m,
    int const* n,
    double* A,
    int const* lda,
    double* tau,
    double* work,
    int const* lwork,
    int* info
);

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

void dgeqp3_
(
    int const* m,
    int const* n,
    double* A,
    int const* lda,
    int* jpvt,
    double* tau,
    double* work,
    int const* lwork,
    int* info
);

double drand(double amin, double amax)
{
    double range = (amax - amin); 
    double div = RAND_MAX / range;
    return amin + (rand() / div);
}

int cond0(n, j) {return 0;}
int cond1(n, j) {return 1;}
int cond2(n, j) {return j < n/2;}
int cond3(n, j) {return j / 2 * 2 == j;}
int cond4(n, j) {return j < n/4 || j >= n*3/4;}
int cond5(n, j) {return j >= n/2;}
int cond6(n, j) {return (((double) abs(rand()))/((double) RAND_MAX) < 0.5);}

void create_matrix(int m, int n, int type, double* A)
{
    int i, j, ND = 0;
    srand(0);
    int (*cond)(int, int);
    switch (type)
    {
        case 0: cond = &cond0; break;
        case 1: cond = &cond1; break;
        case 2: cond = &cond2; break;
        case 3: cond = &cond3; break;
        case 4: cond = &cond4; break;
        case 5: cond = &cond5; break;
        case 6: cond = &cond6; break;
    }
    for (j = 0; j<n; j++)
    {
        if (cond(n,j))
        {
            for (i = 0; i<m; i++)
            {
                A[i+j*m] = drand(-1., 1.);
            }
        }
        else
        {
            ND++;
            for (i = 0; i<m; i++)
            {
                A[i+j*m] = 0.;
            }
        }
//        A[i*n+i] = ((i / 2) * 2 == i) ? 1. : 0. ;
    }
//    printf(" ND orig    %d\n", ND);
}

int main(int argc, char** argv)
{
    int m = atoi(argv[1]);
    int n = atoi(argv[2]);
    int type = atoi(argv[3]);
    int j;
    double* A = (double*) malloc(m * n * sizeof(double));
    int lda = m; //XXX
    double* tau = (double*) malloc(min(m,n) * sizeof(double));
    int nd;
    int* deff = (int*) malloc(n * sizeof(int));
    int lwork = min(m,n) * min(m,n);
    double* work = (double*) malloc(lwork * sizeof(double));
    int info = 0;
    int* jpvt = (int*) malloc(n * sizeof(int));

    clock_t tic, toc;

    // WarmupQR

    create_matrix(m, n, type, A);

    tic = clock();
    dgeqrf_ ( &m, &n, A, &lda, tau, work, &lwork, &info);
    toc = clock();
    printf("WarmupQR : %.2f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

    // QR

    create_matrix(m, n, type, A);

    tic = clock();
    dgeqrf_ ( &m, &n, A, &lda, tau, work, &lwork, &info);
    toc = clock();
    printf("QR       : %.2f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

    // PAQR

    create_matrix(m, n, type, A);
//    for (j = 0; j<n; j++)
//    {
//        deff[j] = -1;
//    }

    tic = clock();
    dgepoqrf_
    ( &m, &n, A, &lda, tau, &nd, deff, work, &lwork, &info);
    toc = clock();

//    for (j = 0; j<n; j++)
//    {
//        printf("%d ", deff[j]);
//    }
//    printf("\n");

    printf("PAQR     : %.2f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

    // QRCP

    create_matrix(m, n, type, A);

    tic = clock();
    dgeqp3_ ( &m, &n, A, &lda, jpvt, tau, work, &lwork, &info);
    toc = clock();
    printf("QRCP     : %.2f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

    free(A);
    free(tau);
    free(deff);
    free(work);

    return 0;
}

