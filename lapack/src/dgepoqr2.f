*> \brief \b DGEQR2 computes the QR factorization of a general rectangular matrix using an unblocked algorithm.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DGEQR2 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeqr2.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeqr2.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeqr2.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGEQR2( M, N, A, LDA, TAU, WORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, M, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGEQR2 computes a QR factorization of a real m-by-n matrix A:
*>
*>    A = Q * ( R ),
*>            ( 0 )
*>
*> where:
*>
*>    Q is a m-by-m orthogonal matrix;
*>    R is an upper-triangular n-by-n matrix;
*>    0 is a (m-n)-by-n zero matrix, if m > n.
*>
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the m by n matrix A.
*>          On exit, the elements on and above the diagonal of the array
*>          contain the min(m,n) by n upper trapezoidal matrix R (R is
*>          upper triangular if m >= n); the elements below the diagonal,
*>          with the array TAU, represent the orthogonal matrix Q as a
*>          product of elementary reflectors (see Further Details).
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] TAU
*> \verbatim
*>          TAU is DOUBLE PRECISION array, dimension (min(M,N))
*>          The scalar factors of the elementary reflectors (see Further
*>          Details).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup doubleGEcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The matrix Q is represented as a product of elementary reflectors
*>
*>     Q = H(1) H(2) . . . H(k), where k = min(m,n).
*>
*>  Each H(i) has the form
*>
*>     H(i) = I - tau * v * v**T
*>
*>  where tau is a real scalar, and v is a real vector with
*>  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
*>  and tau in TAU(i).
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DGEPOQR2( M, N, A, LDA, TAU, ND, DEFF, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
      INTEGER            ND
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
      INTEGER            DEFF( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, K
      DOUBLE PRECISION   AII
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARF, DLARFG, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGEQR2', -INFO )
         RETURN
      END IF
*
      K = MIN( M, N )
      ND = 0
*
      DO 10 I = 1, K
*
*        Generate elementary reflector H(i) to annihilate A(i+1-ND:m,i)
*
         CALL DLARFGPO( M-I+ND+1,
     $                  A( I-ND, I    ), A( MIN( I+1-ND, M ), I    ), 1,
     $                  A( I-ND, I-ND ), A( MIN( I+1-ND, M ), I-ND ), 1,
     $                  TAU( I ) )
         IF( abs(A( I-ND, I-ND )) < EPSILON(0.) ) THEN
             DEFF(I) = 1
             ND = ND + 1
         ELSE
             DEFF(I) = 0
             IF( I.LT.N ) THEN
*    
*               Apply H(i) to A(i:m,i+1:n) from the left
*    
                AII = A( I-ND, I-ND )
                A( I-ND, I-ND ) = ONE
                CALL DLARF( 'Left', M-I+ND+1, N-I, A( I-ND, I-ND ), 1,
     $                  TAU( I ), A( I-ND, I+1 ), LDA, WORK )
                A( I-ND, I-ND ) = AII
             END IF
         END IF
   10 CONTINUE
      RETURN
*
*     End of DGEQR2
*
      END