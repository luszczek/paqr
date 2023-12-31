/*
    -- MAGMA (version 2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date

       @generated from testing/testing_zgeqrf_native_cusolver.cpp, normal z -> d, Tue Apr 12 21:49:53 2022
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

// 
void
get_QR_error( 
    magma_opts* opts, 
    magma_int_t M, magma_int_t N, 
    magmaDouble_ptr h_A,
    magmaDouble_ptr h_R,
    magma_int_t lda, 
    magmaDouble_ptr d_A, 
    magma_int_t ldda, 
    magmaDouble_ptr tau, 
    magmaDouble_ptr dT, 
    magmaDouble_ptr h_work, 
    magma_int_t lwork, 
    double *error, 
    double *error2)
{
    magma_int_t min_mn = min(M, N);
    magma_int_t nb     = magma_get_dgeqrf_nb( M, N );
    magma_int_t info   = 0;
    double Anorm;
    const double             d_neg_one = MAGMA_D_NEG_ONE;
    const double             d_one     = MAGMA_D_ONE;
    const double c_neg_one = MAGMA_D_NEG_ONE;
    const double c_one     = MAGMA_D_ONE;
    const double c_zero    = MAGMA_D_ZERO;
    
    if ( opts->version == 3 ) {
        // copy diagonal blocks of R back to A
        for( int i=0; i < min_mn-nb; i += nb ) {
            magma_int_t ib = min( min_mn-i, nb );
            magmablas_dlacpy( MagmaUpper, ib, ib, &dT[min_mn*nb + i*nb], nb, &d_A[ i + i*ldda ], ldda, opts->queue );
        }
    }

    /* =====================================================================
       Check the result, following zqrt01 except using the reduced Q.
       This works for any M,N (square, tall, wide).
       Only for version 2, which has LAPACK complaint output.
       Or   for version 3, after restoring diagonal blocks of A above.
       =================================================================== */
    magma_dgetmatrix( M, N, d_A, ldda, h_R, lda, opts->queue );
    
    magma_int_t ldq = M;
    magma_int_t ldr = min_mn;
    double *Q, *R;
    double *work;
    TESTING_CHECK( magma_dmalloc_cpu( &Q,    ldq*min_mn ));  // M by K
    TESTING_CHECK( magma_dmalloc_cpu( &R,    ldr*N ));       // K by N
    TESTING_CHECK( magma_dmalloc_cpu( &work, min_mn ));
    
    // generate M by K matrix Q, where K = min(M,N)
    lapackf77_dlacpy( "Lower", &M, &min_mn, h_R, &lda, Q, &ldq );
    lapackf77_dorgqr( &M, &min_mn, &min_mn, Q, &ldq, tau, h_work, &lwork, &info );
    assert( info == 0 );
    
    // copy K by N matrix R
    lapackf77_dlaset( "Lower", &min_mn, &N, &c_zero, &c_zero, R, &ldr );
    lapackf77_dlacpy( "Upper", &min_mn, &N, h_R, &lda,        R, &ldr );
    
    // error = || R - Q^H*A || / (N * ||A||)
    blasf77_dgemm( "Conj", "NoTrans", &min_mn, &N, &M,
                   &c_neg_one, Q, &ldq, h_A, &lda, &c_one, R, &ldr );
    Anorm = lapackf77_dlange( "1", &M,      &N, h_A, &lda, work );
    (*error) = lapackf77_dlange( "1", &min_mn, &N, R,   &ldr, work );
    if ( N > 0 && Anorm > 0 )
        (*error) /= (N*Anorm);
    
    // set R = I (K by K identity), then R = I - Q^H*Q
    // error = || I - Q^H*Q || / N
    lapackf77_dlaset( "Upper", &min_mn, &min_mn, &c_zero, &c_one, R, &ldr );
    blasf77_dsyrk( "Upper", "Conj", &min_mn, &M, &d_neg_one, Q, &ldq, &d_one, R, &ldr );
    (*error2) = safe_lapackf77_dlansy( "1", "Upper", &min_mn, R, &ldr, work );
    if ( N > 0 )
        (*error2) /= N;
    
    magma_free_cpu( Q    );  Q    = NULL;
    magma_free_cpu( R    );  R    = NULL;
    magma_free_cpu( work );  work = NULL;    
}

void get_QR_error_solve( 
    magma_opts* opts, 
    magma_int_t M, magma_int_t N, 
    magmaDouble_ptr h_A, magma_int_t lda, 
    magmaDouble_ptr d_A, magma_int_t ldda, 
    magmaDouble_ptr tau, 
    magmaDouble_ptr dT, 
    double *error)
{
    const double c_zero    = MAGMA_D_ZERO;
    const double c_one     = MAGMA_D_ONE;
    const double c_neg_one = MAGMA_D_NEG_ONE;
    magma_int_t ISEED[4] = {0,0,0,1};
    const magma_int_t ione = 1;
    magma_int_t info   = 0;
    double tmp[1];
    
    /* =====================================================================
       Check the result by solving consistent linear system, A*x = b.
       Only for versions 1 & 3 with M >= N.
       =================================================================== */
    magma_int_t lwork2;
    double *x, *b, *hwork;
    magmaDouble_ptr d_B;

    // initialize RHS, b = A*random
    TESTING_CHECK( magma_dmalloc_cpu( &x, N ));
    TESTING_CHECK( magma_dmalloc_cpu( &b, M ));
    lapackf77_dlarnv( &ione, ISEED, &N, x );
    blasf77_dgemv( "Notrans", &M, &N, &c_one, h_A, &lda, x, &ione, &c_zero, b, &ione );
    // copy to GPU
    TESTING_CHECK( magma_dmalloc( &d_B, M ));
    magma_dsetvector( M, b, 1, d_B, 1, opts->queue );

    if ( opts->version == 1 ) {
        // allocate hwork
        magma_dgeqrs_gpu( M, N, 1,
                          d_A, ldda, tau, dT,
                          d_B, M, tmp, -1, &info );
        lwork2 = (magma_int_t)MAGMA_D_REAL( tmp[0] );
        TESTING_CHECK( magma_dmalloc_cpu( &hwork, lwork2 ));

        // solve linear system
        magma_dgeqrs_gpu( M, N, 1,
                          d_A, ldda, tau, dT,
                          d_B, M, hwork, lwork2, &info );
        if (info != 0) {
            printf("magma_dgeqrs returned magma_error %lld: %s.\n",
                   (long long) info, magma_strerror( info ));
        }
        magma_free_cpu( hwork );
    }
    #ifdef HAVE_CUBLAS
    else if ( opts->version == 3 ) {
        // allocate hwork
        magma_dgeqrs3_gpu( M, N, 1,
                           d_A, ldda, tau, dT,
                           d_B, M, tmp, -1, &info );
        lwork2 = (magma_int_t)MAGMA_D_REAL( tmp[0] );
        TESTING_CHECK( magma_dmalloc_cpu( &hwork, lwork2 ));

        // solve linear system
        magma_dgeqrs3_gpu( M, N, 1,
                           d_A, ldda, tau, dT,
                           d_B, M, hwork, lwork2, &info );
        if (info != 0) {
            printf("magma_dgeqrs3 returned magma_error %lld: %s.\n",
                   (long long) info, magma_strerror( info ));
        }
        magma_free_cpu( hwork );
    }
    #endif
    else {
        printf( "Unknown version %lld\n", (long long) opts->version );
        return;
    }
    magma_dgetvector( N, d_B, 1, x, 1, opts->queue );

    // compute r = Ax - b, saved in b
    blasf77_dgemv( "Notrans", &M, &N, &c_one, h_A, &lda, x, &ione, &c_neg_one, b, &ione );

    // compute residual |Ax - b| / (max(m,n)*|A|*|x|)
    double norm_x, norm_A, norm_r, work[1];
    norm_A = lapackf77_dlange( "F", &M, &N, h_A, &lda, work );
    norm_r = lapackf77_dlange( "F", &M, &ione, b, &M, work );
    norm_x = lapackf77_dlange( "F", &N, &ione, x, &N, work );

    magma_free_cpu( x );
    magma_free_cpu( b );
    magma_free( d_B );

    (*error) = norm_r / (max(M,N) * norm_A * norm_x);
        
}
/* ////////////////////////////////////////////////////////////////////////////
   -- Testing dgeqrf
*/
int main( int argc, char** argv)
{
    TESTING_CHECK( magma_init() );
    magma_print_environment();

    const double c_zero    = MAGMA_D_ZERO;
    const magma_int_t        ione      = 1;
    real_Double_t    gflops, magma_perf, magma_time, cusolver_perf=0, cusolver_time=0, cpu_perf=0, cpu_time=0;
    double           magma_error=0, magma_error2=0, cusolver_error=0, cusolver_error2=0;
    double *h_A, *h_R, *tau, *h_work, tmp[1];
    magmaDouble_ptr d_A, dT;
    magma_int_t M, N, n2, lda, ldda, lwork, info, min_mn, nb, size;
    magma_int_t ISEED[4] = {0,0,0,1};
    
    magma_opts opts;
    opts.parse_opts( argc, argv );
    
    int status = 0;
    double tol = opts.tolerance * lapackf77_dlamch("E");
    
    // version 3 can do either check
    if (opts.check == 1 && opts.version == 1) {
        opts.check = 2;
        printf( "%% version 1 requires check 2 (solve A*x=b)\n" );
    }
    if (opts.check == 2 && opts.version == 2) {
        opts.check = 1;
        printf( "%% version 2 requires check 1 (R - Q^H*A)\n" );
    }
    
    printf( "%% version %lld\n", (long long) opts.version );
    if ( opts.check == 1 ) {
        printf("%%                                                                                ----- MAGMA ERROR -------   ----- CUSOLVER ERROR ----\n");
        printf("%%   M     N   CPU Gflop/s (sec)   MAGMA Gflop/s (sec)   CUSOLVER Gflop/s (sec)   |R - Q^H*A|   |I - Q^H*Q|   |R - Q^H*A|   |I - Q^H*Q|\n");
        printf("%%=====================================================================================================================================\n");
    }
    else {
        printf("%%                                                                                 MAGMA ERROR   CUSOLVER ERROR\n");
        printf("%%   M     N   CPU Gflop/s (sec)   MAGMA Gflop/s (sec)   CUSOLVER Gflop/s (sec)    |b - A*x|      |b - A*x|    \n");
        printf("%%=============================================================================================================\n");
    }
    for( int itest = 0; itest < opts.ntest; ++itest ) {
        for( int iter = 0; iter < opts.niter; ++iter ) {
            M = opts.msize[itest];
            N = opts.nsize[itest];
            min_mn = min( M, N );
            lda    = M;
            n2     = lda*N;
            ldda   = magma_roundup( M, opts.align );  // multiple of 32 by default
            nb     = magma_get_dgeqrf_nb( M, N );
            gflops = FLOPS_DGEQRF( M, N ) / 1e9;
            
            // query for workspace size
            lwork = -1;
            lapackf77_dgeqrf( &M, &N, NULL, &M, NULL, tmp, &lwork, &info );
            lwork = (magma_int_t)MAGMA_D_REAL( tmp[0] );
            
            TESTING_CHECK( magma_dmalloc_cpu( &tau,    min_mn ));
            TESTING_CHECK( magma_dmalloc_cpu( &h_A,    n2     ));
            TESTING_CHECK( magma_dmalloc_cpu( &h_work, lwork  ));
            
            TESTING_CHECK( magma_dmalloc_pinned( &h_R,    n2     ));
            
            TESTING_CHECK( magma_dmalloc( &d_A,    ldda*N ));
            
            if ( opts.version == 1 || opts.version == 3 ) {
                size = (2*min(M, N) + magma_roundup( N, 32 ) )*nb;
                TESTING_CHECK( magma_dmalloc( &dT, size ));
                magmablas_dlaset( MagmaFull, size, 1, c_zero, c_zero, dT, size, opts.queue );
            }
            
            /* Initialize the matrix */
            lapackf77_dlarnv( &ione, ISEED, &n2, h_A );
            lapackf77_dlacpy( MagmaFullStr, &M, &N, h_A, &lda, h_R, &lda );
            
            /* ====================================================================
               Performs operation using MAGMA
               =================================================================== */
            magma_dsetmatrix( M, N, h_R, lda, d_A, ldda, opts.queue );
            nb = magma_get_dgeqrf_nb( M, N );
            
            magma_time = magma_wtime();
            if ( opts.version == 1 ) {
                // stores dT, V blocks have zeros, R blocks inverted & stored in dT
                magma_dgeqrf_gpu( M, N, d_A, ldda, tau, dT, &info );
            }
            else if ( opts.version == 2 ) {
                // LAPACK complaint arguments
                magma_dgeqrf2_gpu( M, N, d_A, ldda, tau, &info );
            }
            #ifdef HAVE_CUBLAS
            else if ( opts.version == 3 ) {
                // stores dT, V blocks have zeros, R blocks stored in dT
                magma_dgeqrf3_gpu( M, N, d_A, ldda, tau, dT, &info );
            }
            #endif
            else {
                printf( "Unknown version %lld\n", (long long) opts.version );
                return -1;
            }
            magma_time = magma_wtime() - magma_time;
            magma_perf = gflops / magma_time;
            if (info != 0) {
                printf("magma_dgeqrf returned magma_error %lld: %s.\n",
                       (long long) info, magma_strerror( info ));
            }
            
            
            if ( opts.check == 1 && (opts.version == 2 || opts.version == 3) ) {
                get_QR_error( &opts, M, N, 
                              h_A, h_R, lda, 
                              d_A, ldda, 
                              tau, dT, 
                              h_work, lwork, 
                              &magma_error, &magma_error2);
            }
            else if ( opts.check == 2 && M >= N && (opts.version == 1 || opts.version == 3) ) {
                get_QR_error_solve( &opts, M, N, 
                                    h_A, lda, 
                                    d_A, ldda, 
                                    tau, dT, &magma_error);
            }
            
            /* ====================================================================
               Performs operation using CUSOLVER
               =================================================================== */
            if( opts.version == 2){
                magma_dsetmatrix( M, N, h_A, lda, d_A, ldda, opts.queue );
            
                cusolver_time = magma_wtime();
                magma_dgeqrf_cusolver_gpu(M, N, d_A, ldda, tau, &info );
                magma_queue_sync( NULL );
                cusolver_time = magma_wtime() - cusolver_time;
                cusolver_perf = gflops / cusolver_time;
                if (info != 0) {
                    printf("magma_dgeqrf_cusolver_gpu returned error %lld: %s.\n",
                           (long long) info, magma_strerror( info ));
                    cusolver_perf = 0.0;
                    cusolver_time = 0.0;
                }
                else{
                    if ( opts.check == 1 && (opts.version == 2 || opts.version == 3) ) {
                        get_QR_error( &opts, M, N, 
                                      h_A, h_R, lda, 
                                      d_A, ldda, 
                                      tau, dT, 
                                      h_work, lwork, 
                                      &cusolver_error, &cusolver_error2);
                    }
                    else if ( opts.check == 2 && M >= N && (opts.version == 1 || opts.version == 3) ) {
                            get_QR_error_solve( &opts, M, N, 
                                            h_A, lda, 
                                            d_A, ldda, 
                                            tau, dT, &cusolver_error);
                    }
                }
            }
            /* =====================================================================
               Performs operation using LAPACK
               =================================================================== */
            if ( opts.lapack ) {
                cpu_time = magma_wtime();
                lapackf77_dgeqrf( &M, &N, h_A, &lda, tau, h_work, &lwork, &info );
                cpu_time = magma_wtime() - cpu_time;
                cpu_perf = gflops / cpu_time;
                if (info != 0) {
                    printf("lapackf77_dgeqrf returned magma_error %lld: %s.\n",
                           (long long) info, magma_strerror( info ));
                }
            }
            
            /* =====================================================================
               Print performance and magma_error.
               =================================================================== */
            printf("%5lld %5lld   ", (long long) M, (long long) N );
            if ( opts.lapack ) {
                printf( "%7.2f ( %7.2f )", cpu_perf, cpu_time );
            }
            else {
                printf("  ---   (  ---  )" );
            }
            // print magma
            printf( "   %7.2f ( %7.2f )   %7.2f ( %7.2f )", magma_perf, magma_time, cusolver_perf, cusolver_time );
            if ( opts.check == 1 ) {
                bool okay = (magma_error < tol && magma_error2 < tol && cusolver_error < tol && cusolver_error2 < tol);
                status += ! okay;
                printf( "%11.2e   %11.2e   %11.2e   %11.2e   %s\n", magma_error, magma_error2, cusolver_error, cusolver_error2, (okay ? "ok" : "failed") );
            }
            else if ( opts.check == 2 ) {
                if ( M >= N ) {
                    bool okay = (magma_error < tol && cusolver_error < tol);
                    status += ! okay;
                    printf( "%10.2e   %10.2e   %s\n", magma_error, cusolver_error, (okay ? "ok" : "failed") );
                }
                else {
                    printf( "(magma_error check only for M >= N)\n" );
                }
            }
            else {
                printf( "    ---\n" );
            }
            
            magma_free_cpu( tau    );
            magma_free_cpu( h_A    );
            magma_free_cpu( h_work );
            
            magma_free_pinned( h_R );
            
            magma_free( d_A );
            
            if ( opts.version == 1 || opts.version == 3 ) {
                magma_free( dT );
            }
            
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
