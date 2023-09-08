/*
    -- MAGMA (version 2.0) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date

       @author Ahmad Abdelfattah

       @precisions normal z -> s d c
*/

////////////////////////////////////////////////////////////////////////////////
#define SLDA(n)              ( (((n)+1)%4) == 0 ? (n) : (n+1) )
#define sA(i,j)              sA[(j) * slda + (i)]
#define NTCOL(M)             ((M > 32) ? 1 : 2)
#define _TPC_                (16)

////////////////////////////////////////////////////////////////////////////////
//             For sm kernel
////////////////////////////////////////////////////////////////////////////////
static __device__ __inline__
void zgeqr2_compute_vtA_device(
        int &m, int &n, int &j,
        magmaDoubleComplex *sA, const int &slda,
        magmaDoubleComplex *sY, magmaDoubleComplex &tau,
        magmaDoubleComplex *sTmp,
        const int &tx, const int &ntx)
{
    magmaDoubleComplex zsum = MAGMA_Z_ZERO;

    const int ncols= n-j-1;
    const int tpc  = ntx / ncols; // threads-per-column
    const int nath = ncols * tpc; // # active threads
    const int tx_  = tx % tpc;
    const int ty_  = tx / tpc;
    sTmp += tpc * ty_;
    if( tx < nath ) {
        for(int i = tx_+j; i < m; i+=tpc) {
            zsum += MAGMA_Z_CONJ( sA(i,j) ) * sA(i,ty_+j+1);
        }
        sTmp[ tx_ ] = zsum;
    }
    __syncthreads();
    // reduce
    if( tx < nath && tx_ == 0) {
        zsum = MAGMA_Z_ZERO;
        for(int i = 0; i < tpc; i++) {
            zsum += sTmp[i];
        }

        sY[ty_+j+1] = zsum * MAGMA_Z_CONJ( tau );; // sTmp differs based on the value of ty_;
    }
}

////////////////////////////////////////////////////////////////////////////////
//             For sm kernel
////////////////////////////////////////////////////////////////////////////////
static __device__ __inline__
void zgeqr2_compute_norm(
        int n,
        magmaDoubleComplex* x, double* dx,
        const int &tx, const int &ntx)
{
    double sum = MAGMA_D_ZERO;
    for(int itx = tx; itx < n; itx+=ntx) {
        sum += MAGMA_Z_REAL( x[itx] ) * MAGMA_Z_REAL( x[itx] ) +
               MAGMA_Z_IMAG( x[itx] ) * MAGMA_Z_IMAG( x[itx] ) ;
    }
    dx[ tx ] = sum;
    // there is a sync at the beginning & end of magma_sum_reduce_n
    __syncthreads();
    // at this point the length of dx is <= ntx (which is 1024 max.)
    if ( ntx >  512 ) { if ( tx <  512 && tx +  512 < ntx ) { dx[tx] += dx[tx+ 512]; }  __syncthreads(); }
    if ( ntx >  256 ) { if ( tx <  256 && tx +  256 < ntx ) { dx[tx] += dx[tx+ 256]; }  __syncthreads(); }
    if ( ntx >  128 ) { if ( tx <  128 && tx +  128 < ntx ) { dx[tx] += dx[tx+ 128]; }  __syncthreads(); }
    if ( ntx >   64 ) { if ( tx <   64 && tx +   64 < ntx ) { dx[tx] += dx[tx+  64]; }  __syncthreads(); }
    // continue with serial sum
    sum = MAGMA_D_ZERO;
    if( tx == 0 ) {
        for( int i = 0; i < min(ntx,64); i++ ) {
            sum += dx[i];
        }
        dx[0] = sum;
    }
    __syncthreads();
}
