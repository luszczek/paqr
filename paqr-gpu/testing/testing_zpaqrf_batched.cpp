/*
    -- MAGMA (version 2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date

       @author Tingxing Dong
       @author Azzam Haidar

       @precisions normal z -> s d c

*/

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include "flops.h"
#include "magma_v2.h"
#include "magma_lapack.h"
#include "testings.h"

#if defined(_OPENMP)
#include <omp.h>
#include "../control/magma_threadsetting.h"  // internal header
#endif

#define PRECISION_z

////////////////////////////////////////////////////////////////////////////////
#ifdef PRECISION_d
// Redefine printf to only be on when compiled with DEBUG defined
#ifdef DEBUG
#define dbprintf(...) printf(__VA_ARGS__)
#else
#define dbprintf(...)
#endif

// Multinomial choose(p,q) - probably a better way to do this in C
inline int polychoose(int* p, int* q)
{
  int ret = 1;
  // tgamma(p+1) = p factorial
  for (int d=0; d < 3; d++)
    ret *= tgamma(p[d]+1) / (tgamma(q[d]+1)*tgamma(p[d]-q[d]+1));
  return ret;
}

// Create matrices A for Pth-order interpolation in 3D
// --> A * x \approx b , A is Nrow x Ncol, Nrow > Ncol
// inputs:
// radius - determines order of polynomial, values 1-5 are reasonable in 3D
// batch - number of matrices to generate
// outputs:
// nrow, ncol - size of each matrix
// outmat - pointer to memory for matrix entries, col-major, batch index:
//    row1, col1, matrix1
//    row2, col1, matrix1
//    ...
//    rowN-1, colN, matrix1
//    rowN, colN, matrix1
//    row1, col1, matrix2
//    row2, col1, matrix2
//    ...
int genmat(int radius, int batch, int* outrow, int* outcol, double** outmat)
{
    int R = radius; // radius of interpolation stencil
    int P = 2*radius+1; // order of polynomial interpolation
    // int N = 1 + (R*(8 + R*(6 + R*4)))/3; // star-shaped interpolation points
    int N = (2*R+1)*(2*R+1)*(2*R+1); // box-shaped interpolation points
    // Size of the matrix, Nrow x Ncol, Nrow > Ncol for WLS
    int Nrow = N;
    *outrow = Nrow;
    int Ncol = (P+1)*(P+2)*(P+3)/6;
    *outcol = Ncol;
    dbprintf("Radius = %d --> order = %d\n", R, P);
    dbprintf("(Nrow, Ncol) = (%d, %d)\n", Nrow, Ncol);
    int Nmat = batch; // number of matrices to generate
    dbprintf("batch size Nmat = %d \n", Nmat);
    // Allocate memory for return variable
    double* mat = (double*) malloc(Nmat * Nrow * Ncol * sizeof(double));
    *outmat = mat;

    // Determine what size problem you'll need to generate the batch
    // int Nx = (2*R+1); // size of each cube we'll loop over
    int Nx = R+1; // size of each cube we'll loop over
    int Npts = Nx*Nx*Nx; // total points of each cube we'll loop over
    dbprintf("size of cube Nx = %d \n", Nx);
    int Ncube = 1+Nmat/Npts; // number of cubes
    int Ntot = Npts*Ncube; // number of total points to generate moments for
    dbprintf("num cubes Ncube = %d, Npts = %d, Ntot = %d \n", Ncube, Npts, Ntot);

    int row, col, i, j, k, x, y, z;
    double x0, y0, z0;
    int* powers = (int*) malloc(3 * Ncol * sizeof(int));
    double* mom = (double*) malloc(Ntot * Ncol * sizeof(double));

    // set the power values in 2D array, [Ncol][3]
    int** pvals = (int**) malloc(Ncol * sizeof(int*));
    int ix=0;
    for (int pk=0; pk <= P; pk++)
    for (int pj=0; pj <= P; pj++)
    for (int pi=0; pi <= P; pi++)
    {
      if (pi+pj+pk <= P) // ignore otherwise
      {
        pvals[ix] = (powers+3*ix);
        pvals[ix][0]=pi;
        pvals[ix][1]=pj;
        pvals[ix][2]=pk;
        ix++;
      }
    }

    dbprintf("Power list:\n");
    for (ix=0; ix < Ncol; ix++)
      dbprintf("  pvals[%d]=(%d,%d,%d)\n", ix,
        pvals[ix][0], pvals[ix][1], pvals[ix][2]);
    dbprintf("\n");

    // Moment indexes for centroids, treated differently
    int ploc[3] = {1, P+1, (P+1)*(P+2)/2};

    // Loop through the pts, creating random center and calculating moments
    double** mvals = (double**) malloc(Ntot * sizeof(double*));
    for (int ip=0; ip < Ntot; ip++)
    {
      mvals[ip] = (mom+Ncol*ip);
      int ix = ip % Npts; // loop around after each unit of cube Npts
      i = ix % Nx;
      j = ix/Nx % Nx;
      k = ix/(Nx*Nx);

      //  Calculate random centroid location in [0,1]
      double xc[3], dx[3];
      for (int d=0; d < 3; d++)
      {
        xc[d] = (double) rand() / (double) RAND_MAX;
        // For debugging use this line instead for regular cell moments
        // xc[d] = 0.5;
        //  Assume dx is the min of distance to nearest grid line
        dx[d] = 2.*fabs( (xc[d] > (1-xc[d])) ? 1-xc[d] : xc[d] );
      }
      xc[0] += i; // shift centroid to this index
      xc[1] += j; // y centroid
      xc[2] += k; // z centroid
      // Copy into their correct moment indexes
      mvals[ip][0] = 1.;
      mvals[ip][ploc[0]] = xc[0]; // x centroid
      mvals[ip][ploc[1]] = xc[1]; // y centroid
      mvals[ip][ploc[2]] = xc[2]; // z centroid
      // Fill in the others - relative to centroid
      for (int ic=2; ic < Ncol; ic++) // skip first 2 moments we have set
      {
        if (ic==ploc[1] || ic==ploc[2]) // skip others we have set
          continue;
        // Set all the other moments as if a slab dimensions dx[d]
        double tmp = 1;
        int even = 1;
        for (int d=0; d < 3; d++)
        {
          int pd = pvals[ic][d];
          even = even && !(pd % 2);
          tmp *= pow(dx[d]/2., pd) / (double) (pd+1);
        }
        mvals[ip][ic] = (even) ? tmp : 0; // Any odd moments are 0 about xc
      }
      double tmp;
      tmp = xc[0]+1;
      dbprintf("  xc[%4d]=( %1.2e, %1.2e, %1.2e)\n", ip, xc[0], xc[1], xc[2]);
      dbprintf("  dx[%4d]=( %1.2e, %1.2e, %1.2e)\n", ip, dx[0], dx[1], dx[2]);
      dbprintf("  mom[%4d] for (%d,%d,%d) = ", ip, i, j, k);
      for (int ic=0; ic < Ncol; ic++)
        dbprintf(" %1.2e, ", mvals[ip][ic]);
      dbprintf("\n");
    }
    dbprintf("\n");

    // For the star-shaped stencil, we'll need offsets
    int** offset = (int**) malloc(Nrow * sizeof(int*));
    int* offvals = (int*) malloc(3*Nrow * sizeof(int));
    ix=0;
    int offsetix = 0;
    dbprintf("Calculating stencil offsets:\n");
    for (int k=-R; k <= R; k++)
    for (int j=-R; j <= R; j++)
    for (int i=-R; i <= R; i++)
    {
      // comment out for box-shaped
      // if (abs(i)+abs(j)+abs(k) <= R) // ignore otherwise for star-shaped
      {
        offset[ix] = (offvals+3*ix);
        offset[ix][0]=i;
        offset[ix][1]=j;
        offset[ix][2]=k;
        if (i==0 && j==0 && k==0)
          offsetix = ix;
        ix++;
      }
    }
    for (int ix=0; ix < Nrow; ix++)
      dbprintf("  offset[%d]=(%d,%d,%d)\n", ix,
        offset[ix][0], offset[ix][1], offset[ix][2]);
    dbprintf("\n");

    // Build up the moment matrices
    double** avals = (double**) malloc(Nmat * sizeof(double*));
    for (int ip=0; ip < Nmat; ip++)
    {
      avals[ip] = (mat+Ncol*Nrow*ip);
      int ix = ip % Npts; // loop around after each unit of cube Npts
      int ic = ip / Npts; // which "cube" we're on
      // Indices of this point location
      i = ix % Nx;
      j = ix/Nx % Nx;
      k = ix/(Nx*Nx);

      // Matrix entries are shifted to this cell's xc
      double xc0[3], dx0[3];
      for (int d=0; d < 3; d++)
        xc0[d] = mvals[ip][ploc[d]]; // save centroid to shift to

      // Copy other cell's moments into the matrix
      dbprintf("For matrix[%4d]:\n", ip);
      for (int r=0; r < Nrow; r++)
      {
        if (r == offsetix) // Do special thing for this cell's moments
        {
          // Copy center cell's moments into the matrix
          for (int ic=0; ic < Ncol; ic++)
            avals[ip][r + ic*Nrow] = mvals[ip][ic];
          // Set the first moments (centroids) to be 0
          for (int d=0; d < 3; d++)
            avals[ip][r + ploc[d]*Nrow] = 0;
          dbprintf("  row[%d] (%d,%d,%d) center\n",r,i,j,k);
          continue;
        }

        // 1. Loop over the offset
        int ir = offset[r][0] + i;
        int jr = offset[r][1] + j;
        int kr = offset[r][2] + k;

        // 2. If any of the offset indicess are <0 or >Nx, zero-fill
        if ((ir < 0 || ir >= Nx) ||
            (jr < 0 || jr >= Nx) ||
            (kr < 0 || kr >= Nx))
        {
          for (int ic=0; ic < Ncol; ic++)
            avals[ip][r + ic*Nrow] = 0;
          dbprintf("  row[%d] (%d,%d,%d) offset, =0\n",r,ir,jr,kr);
          continue;
        }

        // 3. Otherwise, calculate its index and copy moments into
        int ixr = ic*Npts + ir + Nx*(jr + kr*Nx);

        // Calculate the shift then set original centroids to 0
        for (int d=0; d < 3; d++)
          dx0[d] = mvals[ixr][ploc[d]] - xc0[d]; // d-th centroid

        dbprintf("  row[%d] index = %d, offset=(%1.2e,%1.2e,%1.2e)\n",
            r,ixr,dx0[0],dx0[1],dx0[2]);

        // Calculate the shift-based weight, dist^(-P-1)
        double dist = sqrt(dx0[0]*dx0[0]+dx0[1]*dx0[1]+dx0[2]*dx0[2]);
        double weight = pow(dist,-P-1);

        // 4. shift it to the current centroid using dx0 = old - new shift
        for (int ic=0; ic < Ncol; ic++)
        {
          avals[ip][r + ic*Nrow] = 0; // zero to accumulate

          for (int is=0; is < Ncol; is++) // shifted location
          {
            int* p = pvals[ic]; // new moment
            int* q = pvals[is]; // shifted moment
            if ((q[0] > p[0]) || (q[1] > p[1]) || (q[2] > p[2]))
              continue; // skip if not a valid q

            int pCq = polychoose(p,q);
            /*
            dbprintf("  (%d, %d, %d) choose (%d, %d, %d) = %d \n",
                p[0],p[1],p[2],q[0],q[1],q[2],pCq);
            */
            double shift = pCq*pow(dx0[0],p[0]-q[0])
                              *pow(dx0[1],p[1]-q[1])
                              *pow(dx0[2],p[2]-q[2]);
            // Fix for centroid entries which should be 0
            double momq = ((q[0]+q[1]+q[2])==1) ? 0 : mvals[ixr][is];
            // accumulate the shifted moment info
            avals[ip][r + ic*Nrow] += shift*momq;
          }
        }

        dbprintf("  row[%d] index = %d, shifted moments=\n",r,ixr);
        for (int ic=0; ic < Ncol; ic++)
          dbprintf("    %1.2e, ", avals[ip][r + ic*Nrow]);
        dbprintf("\n");

        // multiply by the weights
        for (int ic=0; ic < Ncol; ic++)
          avals[ip][r + ic*Nrow] *= weight;
      }
    }
    dbprintf("\n");

    // Free local memory
    free(avals);
    free(offset);
    free(offvals);
    free(mvals);
    free(pvals);
    free(mom);
    free(powers);

    return 0;
}
#endif

////////////////////////////////////////////////////////////////////////////////
void get_QR_error(magma_int_t M, magma_int_t N, magma_int_t min_mn,
                    magmaDoubleComplex *h_R,  magmaDoubleComplex *h_A, magma_int_t lda,
                    magmaDoubleComplex *tau,
                    magmaDoubleComplex *Q,  magma_int_t ldq,
                    magmaDoubleComplex *R,  magma_int_t ldr,
                    magmaDoubleComplex *h_work,  magma_int_t lwork,
                    double *work, double *error, double *error2)
{
    /* h_R:input the factorized matrix by lapack QR,
       h_A:input the original matrix copy
       tau: input
    */

    const double             d_neg_one = MAGMA_D_NEG_ONE;
    const double             d_one     = MAGMA_D_ONE;
    const magmaDoubleComplex c_neg_one = MAGMA_Z_NEG_ONE;
    const magmaDoubleComplex c_one     = MAGMA_Z_ONE;
    const magmaDoubleComplex c_zero    = MAGMA_Z_ZERO;
    double           Anorm;

    magma_int_t info;

    // generate M by K matrix Q, where K = min(M,N)
    lapackf77_zlacpy( "Lower", &M, &min_mn, h_R, &lda, Q, &ldq );
    lapackf77_zungqr( &M, &min_mn, &min_mn, Q, &ldq, tau, h_work, &lwork, &info );
    assert( info == 0 );

    // copy K by N matrix R
    lapackf77_zlaset( "Lower", &min_mn, &N, &c_zero, &c_zero, R, &ldr );
    lapackf77_zlacpy( "Upper", &min_mn, &N, h_R, &lda,        R, &ldr );

    // error = || R - Q^H*A || / (N * ||A||)
    blasf77_zgemm( "Conj", "NoTrans", &min_mn, &N, &M,
    &c_neg_one, Q, &ldq, h_A, &lda, &c_one, R, &ldr );

    Anorm = lapackf77_zlange( "1", &M,      &N, h_A, &lda, work );
    *error = lapackf77_zlange( "1", &min_mn, &N, R,   &ldr, work );

    if ( N > 0 && Anorm > 0 )
        *error /= (N*Anorm);

    // set R = I (K by K identity), then R = I - Q^H*Q
    // error = || I - Q^H*Q || / N
    lapackf77_zlaset( "Upper", &min_mn, &min_mn, &c_zero, &c_one, R, &ldr );
    blasf77_zherk( "Upper", "Conj", &min_mn, &M, &d_neg_one, Q, &ldq, &d_one, R, &ldr );
    *error2 = safe_lapackf77_zlanhe( "1", "Upper", &min_mn, R, &ldr, work );
    if ( N > 0 )
        *error2 /= N;
}


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing zgeqrf_batched
*/
int main( int argc, char** argv)
{
    TESTING_CHECK( magma_init() );
    magma_print_environment();

    real_Double_t    gflops, magma_perf, magma_time, cublas_perf=0, cublas_time=0, cpu_perf, cpu_time;
    double           magma_error, cublas_error, magma_error2, cublas_error2;

    magmaDoubleComplex *h_A, *h_R, *h_Amagma, *h_R2, *tau, *h_work, tmp[1], unused[1];
    magmaDoubleComplex *d_A, *dtau_magma, *dtau_cublas;

    magmaDoubleComplex **dA_array = NULL;
    magmaDoubleComplex **dtau_array = NULL;

    magma_int_t  *dinfo_magma, *dinfo_cublas, *dflag, *hflag;
    magma_int_t **dflag_array;

    magma_int_t M, N, lda, ldda, lwork, n2, info, min_mn;
    magma_int_t ione     = 1;
    magma_int_t ISEED[4] = {0,0,0,1};
    int status = 0;

    magma_int_t batchCount;
    magma_int_t column;

    magma_opts opts( MagmaOptsBatched );
    opts.parse_opts( argc, argv );
    batchCount = opts.batchcount;

    double tol = 30. * lapackf77_dlamch("E");

    printf("%% BatchCount   M     N   MAGMA Gflop/s (ms)   %s Gflop/s (ms)    CPU Gflop/s (ms)   |R - Q^H*A|_mag   |I - Q^H*Q|_mag   |R - Q^H*A|_cub   |I - Q^H*Q|_cub\n", g_platform_str);
    printf("%%============================================================================================================================================================\n");
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            #ifdef PRECISION_d
            // get just the sizes to proceed with allocation
            double* dtmp;
            int outrow, outcol;
            int radius = (int)(opts.nsize[itest]);
            radius = max(radius, 1);
            radius = min(radius, 5);
            genmat(radius, 1, &outrow, &outcol, &dtmp);
            M     = (int)outrow;
            N     = (int)outcol;
            magma_free_cpu( dtmp );
            #else
            M     = opts.msize[itest];
            N     = opts.nsize[itest];
            #endif
            min_mn = min(M, N);
            lda    = M;
            n2     = lda*N * batchCount;
            ldda = M;
            ldda   = magma_roundup( M, opts.align );  // multiple of 32 by default

            gflops = (FLOPS_ZGEQRF( M, N ) + FLOPS_ZGEQRT( M, N )) / 1e9 * batchCount;

            /* Allocate memory for the matrix */

            // in double precision, the genmat function allocates & initializes the space for h_A
            #ifdef PRECISION_d
            if(opts.version  == 3) {
                TESTING_CHECK( magma_zmalloc_cpu( &h_A,       n2 ));
                lapackf77_zlarnv( &ione, ISEED, &n2, h_A );
            }
            else{
                genmat(radius, batchCount, &outrow, &outcol, &h_A);
            }
            #else
            TESTING_CHECK( magma_zmalloc_cpu( &h_A,       n2 ));
            lapackf77_zlarnv( &ione, ISEED, &n2, h_A );
            #endif

            TESTING_CHECK( magma_zmalloc_cpu( &tau,   min_mn * batchCount ));
            TESTING_CHECK( magma_zmalloc_cpu( &h_Amagma,  n2 ));
            TESTING_CHECK( magma_zmalloc_pinned( &h_R,    n2 ));
            TESTING_CHECK( magma_zmalloc_cpu( &h_R2, n2 ));
            TESTING_CHECK( magma_imalloc_cpu( &hflag, min_mn * batchCount ));

            TESTING_CHECK( magma_zmalloc( &d_A,   ldda*N * batchCount ));
            TESTING_CHECK( magma_zmalloc( &dtau_magma,  min_mn * batchCount ));
            TESTING_CHECK( magma_zmalloc( &dtau_cublas, min_mn * batchCount ));

            TESTING_CHECK( magma_imalloc( &dflag, min_mn * batchCount ));
            TESTING_CHECK( magma_malloc( (void**) &dflag_array, batchCount * sizeof(magma_int_t*) ));

            TESTING_CHECK( magma_imalloc( &dinfo_magma,  batchCount ));
            TESTING_CHECK( magma_imalloc( &dinfo_cublas, batchCount ));

            TESTING_CHECK( magma_malloc( (void**) &dA_array,   batchCount * sizeof(magmaDoubleComplex*) ));
            TESTING_CHECK( magma_malloc( (void**) &dtau_array, batchCount * sizeof(magmaDoubleComplex*) ));

            // to determine the size of lwork
            lwork = -1;
            lapackf77_zgeqrf( &M, &N, unused, &M, unused, tmp, &lwork, &info );
            lwork = (magma_int_t)MAGMA_Z_REAL( tmp[0] );
            lwork = max(lwork, N*N);

            TESTING_CHECK( magma_zmalloc_cpu( &h_work, lwork * batchCount ));

            column = N * batchCount;

            /* normalize all columns */
            for(magma_int_t s = 0; s < N*batchCount; s++) {
                magmaDoubleComplex  beta = MAGMA_Z_ZERO;
                magmaDoubleComplex *hcol = h_A + s * lda;
                double colnorm   = lapackf77_zlange( "F", &M, &ione, hcol, &lda, (double*)h_work );
                if( colnorm == MAGMA_D_ZERO ) colnorm = MAGMA_D_ONE;

                beta = MAGMA_Z_MAKE(1./colnorm, MAGMA_D_ZERO);
                blasf77_zscal( &M, &beta, hcol, &ione );
            }

            /* introduce zero columns */
            #ifndef PRECISION_d
            if(opts.nrhs > 1) {
                for(magma_int_t s = 0; s < batchCount; s++) {
                    magmaDoubleComplex *hTmp = h_A + s * lda * N;
                    for(magma_int_t icol = 0; icol < N; icol++) {
                        if(icol % opts.nrhs == 0) {
                            magmaDoubleComplex *hcol = hTmp + icol * lda;
                            memset(hcol, 0, M*sizeof(magmaDoubleComplex));
                        }
                    }
                }
            }
            #endif

            lapackf77_zlacpy( MagmaFullStr, &M, &column, h_A, &lda, h_R, &lda );

            if(batchCount == 1) {
                magma_zprint(M, N, h_A, lda);
            }

            /* ====================================================================
               Performs operation using MAGMA
               =================================================================== */
            magma_zsetmatrix( M, column, h_R, lda,  d_A, ldda, opts.queue );
            magma_zset_pointer( dA_array, d_A, 1, 0, 0, ldda*N, batchCount, opts.queue );
            magma_zset_pointer( dtau_array, dtau_magma, 1, 0, 0, min_mn, batchCount, opts.queue );
            magma_iset_pointer( dflag_array, dflag, 1, 0, 0, min_mn, batchCount, opts.queue );
            magma_memset_async(dflag, 0, min_mn*batchCount*sizeof(magma_int_t), opts.queue);

            magma_time = magma_sync_wtime( opts.queue );

            if(opts.version == 1 || opts.version == 3) {
                magma_int_t nthreads = (opts.nb == 0 || opts.nb > M) ? N : opts.nb;
                double paqr_tol = lapackf77_dlamch("E") * opts.tolerance;
                info = magma_zpaqr2_fused_sm_batched(
                        M, N,
                        dA_array, 0, 0, ldda,
                        dtau_array, 0,
                        dflag_array, 0,
                        dinfo_magma, nthreads, 0,
                        paqr_tol, batchCount, opts.queue);
            }
            else if(opts.version == 2){
                magma_int_t nthreads = (opts.nb == 0 || opts.nb > M) ? N : opts.nb;
                info = magma_zgeqr2_fused_sm_batched(
                        M, N,
                        dA_array, 0, 0, ldda,
                        dtau_array, 0,
                        dinfo_magma, nthreads, 0,
                        batchCount, opts.queue);
            }

            magma_time = magma_sync_wtime( opts.queue ) - magma_time;
            magma_perf = gflops / magma_time;

            magma_zgetmatrix( M, column, d_A, ldda, h_Amagma, lda, opts.queue );
            magma_zgetvector(min_mn*batchCount, dtau_magma, 1, tau, 1, opts.queue );
            magma_getvector( min_mn * batchCount, sizeof(magma_int_t), dflag, 1, hflag, 1, opts.queue );

            if(batchCount == 1) {
                magma_zprint(M, N, h_Amagma, lda);
                printf("flags = [ ");
                for(int ii = 0; ii < min_mn; ii++) {
                    printf("%d ", hflag[ii]);
                }
                printf("];\n");
                printf("tau = ");
                magma_zprint(1, min_mn, tau, 1);
            }

            if(opts.ngpu == 2) {
                for(magma_int_t s = 0; s < batchCount; s++) {
                    magma_int_t *iflag = hflag + s * min_mn;
                    magma_int_t isum = 0;
                    for(magma_int_t ii = 0; ii < min_mn; ii++) {
                        isum += iflag[ii];
                    }

                    if(isum != 0) {
                        printf("matrix %d found rank deficient -- rank = %d\n", s, min_mn-isum);
                    }
                }
            }

            if (info != 0) {
                printf("magma_zgeqrf_batched returned error %lld: %s.\n",
                       (long long) info, magma_strerror( info ));
            }

            /* ====================================================================
               Performs operation using CUBLAS
               =================================================================== */

            magma_zsetmatrix( M, column, h_R, lda,  d_A, ldda, opts.queue );
            magma_zset_pointer( dA_array, d_A, 1, 0, 0, ldda*N, batchCount, opts.queue );
            magma_zset_pointer( dtau_array, dtau_cublas, 1, 0, 0, min_mn, batchCount, opts.queue );

            cublas_time = magma_sync_wtime( opts.queue );

            int device_info;  // not magma_int_t
            #ifdef MAGMA_HAVE_CUDA
            /* cublasZgeqrfBatched is only available from CUBLAS v6.5 */
            #if CUDA_VERSION >= 6050
            cublasZgeqrfBatched( opts.handle, int(M), int(N),
                                 dA_array, int(ldda), dtau_array,
                                 &device_info, int(batchCount) );
            #endif
            #else
            hipblasZgeqrfBatched( opts.handle, int(M), int(N),
                                 (hipblasDoubleComplex**)dA_array, int(ldda),
                                 (hipblasDoubleComplex**)dtau_array,
                                 &device_info, int(batchCount) );
            #endif
            cublas_time = magma_sync_wtime( opts.queue ) - cublas_time;
            cublas_perf = gflops / cublas_time;

            if (device_info != 0) {
                printf("cublasZgeqrfBatched returned error %lld: %s.\n",
                       (long long) device_info, magma_strerror( device_info ));
            }

            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_wtime();
                // #define BATCHED_DISABLE_PARCPU
                #if !defined (BATCHED_DISABLE_PARCPU) && defined(_OPENMP)
                magma_int_t nthreads = magma_get_lapack_numthreads();
                magma_set_lapack_numthreads(1);
                magma_set_omp_numthreads(nthreads);
                #pragma omp parallel for schedule(dynamic)
                #endif
                for (magma_int_t s=0; s < batchCount; s++)
                {
                    magma_int_t locinfo;
                    lapackf77_zgeqrf(&M, &N, h_A + s * lda * N, &lda, tau + s * min_mn, h_work + s * lwork, &lwork, &locinfo);
                    if (locinfo != 0) {
                        printf("lapackf77_zgeqrf matrix %lld returned error %lld: %s.\n",
                               (long long) s, (long long) locinfo, magma_strerror( locinfo ));
                    }
                }

                #if !defined (BATCHED_DISABLE_PARCPU) && defined(_OPENMP)
                    magma_set_lapack_numthreads(nthreads);
                #endif

                cpu_time = magma_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0) {
                    printf("lapackf77_zgeqrf returned error %lld: %s.\n",
                           (long long) info, magma_strerror( info ));
                }
                printf("%10lld %5lld %5lld    %7.2f (%7.2f)     %7.2f (%7.2f)   %7.2f (%7.2f)",
                       (long long) batchCount, (long long) M, (long long) N,
                       magma_perf,  1000.*magma_time,
                       cublas_perf, 1000.*cublas_time,
                       cpu_perf,    1000.*cpu_time );
            }
            else {
                printf("%10lld %5lld %5lld    %7.2f (%7.2f)     %7.2f (%7.2f)     ---  (  ---  )",
                       (long long) batchCount, (long long) M, (long long) N,
                       magma_perf,  1000.*magma_time,
                       cublas_perf, 1000.*cublas_time );
            }

            if (opts.check) {
                /* =====================================================================
                   Check the MAGMA and CUBLAS results by computing residuals
                   =================================================================== */
                magma_int_t ldq = M;
                magma_int_t ldr = min_mn;
                magmaDoubleComplex *Q, *R;
                double *work;

                TESTING_CHECK( magma_zmalloc_cpu( &Q,    batchCount*ldq*min_mn ));  // M by K
                TESTING_CHECK( magma_zmalloc_cpu( &R,    batchCount*ldr*N ));       // K by N
                TESTING_CHECK( magma_dmalloc_cpu( &work, batchCount*min_mn ));

                /* check magma result */
                magma_error  = 0;
                magma_error2 = 0;
                #pragma omp parallel for reduction(max:magma_error,magma_error2)
                for (int i=0; i < batchCount; i++) {
                    // copy h_Amagma to h_R2 and account for bad column marked by hflag
                    magmaDoubleComplex *hR1 = h_R  + i*lda*N;
                    magmaDoubleComplex *hR2 = h_R2 + i*lda*N;
                    magma_int_t *iflag = hflag + i * min_mn;
                    magma_int_t new_N  = 0;
                    for(magma_int_t ic = 0; ic < N; ic++) {
                        if(iflag[ic] == 0) {
                            memcpy(hR2 + new_N * lda, hR1 + ic * lda, lda * sizeof(magmaDoubleComplex));
                            new_N++;
                        }
                    }

                    if(batchCount == 1) {
                        magma_zprint(M,     N, hR1, lda);
                        magma_zprint(M, new_N, hR2, lda);
                    }

                    // error checking
                    double err, err2;
                    get_QR_error(M, new_N, min(M, new_N),
                             h_Amagma + i*lda*N, hR2, lda, tau + i*min_mn,
                             Q + i*ldq*min_mn, ldq, R + i*ldr*N, ldr, h_work + i*lwork, lwork,
                             work + i*min_mn, &err, &err2);
                    magma_error  = magma_max_nan( err,  magma_error  );
                    magma_error2 = magma_max_nan( err2, magma_error2 );
                }

                /* check cublas result */
                cublas_error  = 0;
                cublas_error2 = 0;
                #if ((defined(MAGMA_HAVE_CUDA) && CUDA_VERSION >= 6050) || defined(MAGMA_HAVE_HIP))
                magma_zgetvector(min_mn*batchCount, dtau_cublas, 1, tau, 1, opts.queue );
                magma_zgetmatrix( M, column, d_A, ldda, h_A, lda, opts.queue );
                #pragma omp parallel for reduction(max:cublas_error,cublas_error2)
                for (int i=0; i < batchCount; i++) {
                    double err, err2;
                    get_QR_error(M, N, min_mn,
                             h_A + i*lda*N, h_R + i*lda*N, lda, tau + i*min_mn,
                             Q + i*ldq*min_mn, ldq, R + i*ldr*N, ldr, h_work + i*lwork, lwork,
                             work + i*min_mn, &err, &err2);
                    cublas_error  = magma_max_nan( err,  cublas_error  );
                    cublas_error2 = magma_max_nan( err2, cublas_error2 );
                }
                #endif

                magma_free_cpu( Q    );  Q    = NULL;
                magma_free_cpu( R    );  R    = NULL;
                magma_free_cpu( work );  work = NULL;

                bool okay = (magma_error < tol && magma_error2 < tol);
                status += ! okay;

                printf("   %15.2e   %15.2e   %15.2e   %15.2e   %s\n",
                       magma_error, magma_error2,
                       cublas_error, cublas_error2,
                       (okay ? "ok" : "failed") );
            }
            else {
                printf("\n");
            }

            magma_free_cpu( tau    );
            magma_free_cpu( h_A    );
            magma_free_cpu( h_Amagma );
            magma_free_cpu( h_R2 );
            magma_free_cpu( h_work );
            magma_free_cpu( hflag );
            magma_free_pinned( h_R    );

            magma_free( d_A   );
            magma_free( dtau_magma  );
            magma_free( dtau_cublas );

            magma_free( dinfo_magma );
            magma_free( dinfo_cublas );

            magma_free( dflag );
            magma_free( dflag_array );

            magma_free( dA_array   );
            magma_free( dtau_array  );

            fflush( stdout );
        }
        if ( opts.niter > 1 ) {
            printf( "\n" );
        }
    }

    opts.cleanup();
    TESTING_CHECK( magma_finalize() );
    return status;
}
