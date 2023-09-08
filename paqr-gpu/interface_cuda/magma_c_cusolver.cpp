/*
    -- MAGMA (version 2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date

       @author Mark Gates
       @author Ahmad Abdelfattah
       @generated from interface_cuda/magma_z_cusolver.cpp, normal z -> c, Tue Apr 12 21:49:48 2022
*/
#include "cusolverDn.h"
#include "magma_internal.h"
#include "error.h"

#define COMPLEX

#ifdef MAGMA_HAVE_CUDA

#define PRECISION_c

#if   defined(PRECISION_z)
#define magma_cgetrf_cusolver               cusolverDnZgetrf
#define magma_cgetrf_cusolver_bufferSize    cusolverDnZgetrf_bufferSize

#define magma_cgeqrf_cusolver               cusolverDnZgeqrf
#define magma_cgeqrf_cusolver_bufferSize    cusolverDnZgeqrf_bufferSize

#define magma_cpotrf_cusolver               cusolverDnZpotrf
#define magma_cpotrf_cusolver_bufferSize    cusolverDnZpotrf_bufferSize

#define magma_cpotrf_batched_cusolver_w     cusolverDnZpotrfBatched
#elif defined(PRECISION_c)
#define magma_cgetrf_cusolver               cusolverDnCgetrf
#define magma_cgetrf_cusolver_bufferSize    cusolverDnCgetrf_bufferSize

#define magma_cgeqrf_cusolver               cusolverDnCgeqrf
#define magma_cgeqrf_cusolver_bufferSize    cusolverDnCgeqrf_bufferSize

#define magma_cpotrf_cusolver    cusolverDnCpotrf
#define magma_cpotrf_cusolver_bufferSize    cusolverDnCpotrf_bufferSize

#define magma_cpotrf_batched_cusolver_w     cusolverDnCpotrfBatched
#elif defined(PRECISION_d)
#define magma_sgetrf_cusolver               cusolverDnDgetrf
#define magma_sgetrf_cusolver_bufferSize    cusolverDnDgetrf_bufferSize

#define magma_sgeqrf_cusolver               cusolverDnDgeqrf
#define magma_sgeqrf_cusolver_bufferSize    cusolverDnDgeqrf_bufferSize

#define magma_spotrf_cusolver               cusolverDnDpotrf
#define magma_spotrf_cusolver_bufferSize    cusolverDnDpotrf_bufferSize

#define magma_spotrf_batched_cusolver_w     cusolverDnDpotrfBatched
#elif defined(PRECISION_s)
#define magma_sgetrf_cusolver               cusolverDnSgetrf
#define magma_sgetrf_cusolver_bufferSize    cusolverDnSgetrf_bufferSize

#define magma_sgeqrf_cusolver               cusolverDnSgeqrf
#define magma_sgeqrf_cusolver_bufferSize    cusolverDnSgeqrf_bufferSize

#define magma_spotrf_cusolver               cusolverDnSpotrf
#define magma_spotrf_cusolver_bufferSize    cusolverDnSpotrf_bufferSize

#define magma_spotrf_batched_cusolver_w     cusolverDnSpotrfBatched
#else
#error "One of PRECISION_{s,d,c,z} must be defined."
#endif


// =============================================================================
// LU
extern "C" magma_int_t
magma_cgetrf_cusolver_gpu(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex_ptr dA, magma_int_t ldda,
    magma_int_t *ipiv,
    magma_int_t *info )
{
    int lwork, hinfo = 0;    /* not magma_int_t */
    int *dipiv, *dinfo;    /* not magma_int_t */
    magma_int_t min_mn = min(m, n);
    magmaFloatComplex_ptr dwork;
    magma_queue_t queue;
    magma_device_t cdev;
    magma_getdevice( &cdev );
    magma_queue_create( cdev, &queue );


    cusolverDnHandle_t handle;
    cusolverDnCreate(&handle);

    magma_cgetrf_cusolver_bufferSize(handle, (int)m, (int)n, (cuFloatComplex*)dA, (int)ldda, &lwork);
    if( MAGMA_SUCCESS != magma_cmalloc(        &dwork, (magma_int_t)lwork  ) ||
        MAGMA_SUCCESS != magma_malloc( (void**)&dinfo, 1      * sizeof(int)) ||
        MAGMA_SUCCESS != magma_malloc( (void**)&dipiv, min_mn * sizeof(int))
      ){
          (*info) = MAGMA_ERR_DEVICE_ALLOC;
          goto cleanup;
    }

    magma_cgetrf_cusolver( handle,
                           (int)m, (int)n,
                           (cuFloatComplex*)dA,    (int)ldda,
                           (cuFloatComplex*)dwork, (int*)dipiv,
                           (int*)dinfo );

    #ifdef MAGMA_ILP64
    int *hipiv;
    magma_malloc_cpu((void**)&hpiv, min_mn * sizeof(int));
    magma_getvector(min_mn, sizeof(int), dipiv, 1, hipiv, 1, queue);
    magma_queue_sync( queue );
    for(int i = 0; i < min_mn; i++){
        ipiv[i] = (magma_int_t)hipiv[i];
    }
    magma_free_cpu( hipiv );
    #else
    magma_getvector(min_mn, sizeof(int), dipiv, 1, ipiv, 1, queue);
    #endif

    magma_getvector(1, sizeof(int), dinfo, 1, &hinfo, 1, queue);    /* copy int, not magma_int_t*/
    (*info) = (magma_int_t)hinfo;


    magma_queue_sync( queue );
cleanup:
    magma_free( dwork );
    magma_free( dinfo );
    magma_free( dipiv );
    cusolverDnDestroy(handle);
    magma_queue_destroy( queue );
    return (*info);
}

// =============================================================================
// QR
extern "C" magma_int_t
magma_cgeqrf_cusolver_gpu(
    magma_int_t m, magma_int_t n,
    magmaFloatComplex_ptr dA, magma_int_t ldda,
    magmaFloatComplex_ptr tau,
    magma_int_t *info )
{
    int lwork, hinfo = 0;    /* not magma_int_t */
    int *dinfo = NULL;      /* not magma_int_t */
    magma_int_t min_mn = min(m, n);
    magmaFloatComplex_ptr dwork=NULL, dtau=NULL;
    magma_queue_t queue;
    magma_device_t cdev;
    magma_getdevice( &cdev );
    magma_queue_create( cdev, &queue );

    cusolverDnHandle_t handle;
    cusolverDnCreate(&handle);

    magma_cgeqrf_cusolver_bufferSize(handle, (int)m, (int)n, (cuFloatComplex*)dA, (int)ldda, &lwork);
    if( ( MAGMA_SUCCESS != magma_cmalloc( &dwork, (magma_int_t)lwork ) ) ||
        ( MAGMA_SUCCESS != magma_cmalloc( &dtau, min_mn )              ) ||
        ( MAGMA_SUCCESS != magma_malloc((void**)&dinfo, 1*sizeof(int)) )
      ){
          (*info) = MAGMA_ERR_DEVICE_ALLOC;
          goto cleanup;
    }

    magma_cgeqrf_cusolver( handle,
                           (int)m, (int)n,
                           (cuFloatComplex*)dA,    (int)ldda,
                           (cuFloatComplex*)dtau,
                           (cuFloatComplex*)dwork, lwork,
                           (int*)dinfo );

    magma_cgetvector(min_mn, dtau, 1, tau, 1, queue);

    magma_getvector(1, sizeof(int), dinfo, 1, &hinfo, 1, queue);    /* copy int, not magma_int_t*/
    (*info) = (magma_int_t)hinfo;

    magma_queue_sync( queue );
cleanup:
    magma_free( dwork );
    magma_free( dtau );
    magma_free( dinfo );
    cusolverDnDestroy(handle);
    magma_queue_destroy( queue );
    return (*info);
}
// =============================================================================
// Cholesky
extern "C" magma_int_t
magma_cpotrf_cusolver_gpu(
    magma_uplo_t uplo, magma_int_t n,
    magmaFloatComplex_ptr dA, magma_int_t ldda,
    magma_int_t *info )
{
    int lwork, hinfo = 0;    /* not magma_int_t */
    int *dinfo;              /* not magma_int_t */
    magmaFloatComplex_ptr dwork;
    cublasFillMode_t cusolver_uplo;

    cusolver_uplo = CUBLAS_FILL_MODE_LOWER;
    if(uplo == MagmaLower)
        cusolver_uplo = CUBLAS_FILL_MODE_LOWER;
    else if(uplo == MagmaUpper)
        cusolver_uplo = CUBLAS_FILL_MODE_UPPER;

    magma_queue_t queue;
    magma_device_t cdev;
    magma_getdevice( &cdev );
    magma_queue_create( cdev, &queue );

    cusolverDnHandle_t handle;
    cusolverDnCreate(&handle);

    magma_cpotrf_cusolver_bufferSize(handle, cusolver_uplo, (int)n, (cuFloatComplex*)dA, (int)ldda, &lwork);

    if( MAGMA_SUCCESS != magma_cmalloc(        &dwork, (magma_int_t)lwork  ) ||
        MAGMA_SUCCESS != magma_malloc( (void**)&dinfo, 1      * sizeof(int))
      ){
          (*info) = MAGMA_ERR_DEVICE_ALLOC;
          goto cleanup;
    }

    magma_cpotrf_cusolver( handle, cusolver_uplo, (int)n,
                           (cuFloatComplex*)dA, (int)ldda, dwork, (int)lwork, (int*)dinfo );


    magma_getvector(1, sizeof(int), dinfo, 1, &hinfo, 1, queue);    /* copy int, not magma_int_t*/
    (*info) = (magma_int_t)hinfo;
    magma_queue_sync( queue );

cleanup:
    magma_free( dwork );
    magma_free( dinfo );
    cusolverDnDestroy(handle);
    magma_queue_destroy( queue );
    return (*info);
}

// =============================================================================
// Batch Cholesky
extern "C" magma_int_t
magma_cpotrf_batched_cusolver(
    magma_uplo_t uplo, magma_int_t n,
    magmaFloatComplex **dA_array, magma_int_t ldda,
    magma_int_t *info_array,  magma_int_t batchCount,
    magma_queue_t queue)
{
    magma_queue_sync( queue );
    cusolverDnHandle_t handle;
    cusolverDnCreate(&handle);
    cublasFillMode_t cusolver_uplo = CUBLAS_FILL_MODE_LOWER;
    if(uplo == MagmaLower)
        cusolver_uplo = CUBLAS_FILL_MODE_LOWER;
    else if(uplo == MagmaUpper)
        cusolver_uplo = CUBLAS_FILL_MODE_UPPER;

    magma_cpotrf_batched_cusolver_w(handle, cusolver_uplo, n, dA_array, ldda, info_array, batchCount);
    cusolverDnDestroy(handle);
    magma_queue_sync( queue );
    return 0;
}


#endif // MAGMA_HAVE_CUDA

#undef COMPLEX
