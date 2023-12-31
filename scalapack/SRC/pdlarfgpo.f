      SUBROUTINE PDLARFGPO( N, ALPHA, IAX, JAX, X, IX, JX, DESCX, INCX,
     $                      TAU , ND, DEFF )
*
*  -- ScaLAPACK auxiliary routine (version 1.7) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     May 1, 1997
*
*     .. Scalar Arguments ..
      INTEGER            IAX, INCX, IX, JAX, JX, N, ND, DEFF
      DOUBLE PRECISION   ALPHA
*     ..
*     .. Array Arguments ..
      INTEGER            DESCX( * )
      DOUBLE PRECISION   TAU( * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  PDLARFG generates a real elementary reflector H of order n, such
*  that
*
*     H * sub( X ) = H * ( x(iax,jax) ) = ( alpha ),   H' * H = I.
*                        (      x     )   (   0   )
*
*  where alpha is a scalar, and sub( X ) is an (N-1)-element real
*  distributed vector X(IX:IX+N-2,JX) if INCX = 1 and X(IX,JX:JX+N-2) if
*  INCX = DESCX(M_).  H is represented in the form
*
*        H = I - tau * ( 1 ) * ( 1 v' ) ,
*                      ( v )
*
*  where tau is a real scalar and v is a real (N-1)-element
*  vector.
*
*  If the elements of sub( X ) are all zero, then tau = 0 and H is
*  taken to be the unit matrix.
*
*  Otherwise  1 <= tau <= 2.
*
*  Notes
*  =====
*
*  Each global data object is described by an associated description
*  vector.  This vector stores the information required to establish
*  the mapping between an object element and its corresponding process
*  and memory location.
*
*  Let A be a generic term for any 2D block cyclicly distributed array.
*  Such a global array has an associated description vector DESCA.
*  In the following comments, the character _ should be read as
*  "of the global array".
*
*  NOTATION        STORED IN      EXPLANATION
*  --------------- -------------- --------------------------------------
*  DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
*                                 DTYPE_A = 1.
*  CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
*                                 the BLACS process grid A is distribu-
*                                 ted over. The context itself is glo-
*                                 bal, but the handle (the integer
*                                 value) may vary.
*  M_A    (global) DESCA( M_ )    The number of rows in the global
*                                 array A.
*  N_A    (global) DESCA( N_ )    The number of columns in the global
*                                 array A.
*  MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
*                                 the rows of the array.
*  NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
*                                 the columns of the array.
*  RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
*                                 row of the array A is distributed.
*  CSRC_A (global) DESCA( CSRC_ ) The process column over which the
*                                 first column of the array A is
*                                 distributed.
*  LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
*                                 array.  LLD_A >= MAX(1,LOCr(M_A)).
*
*  Let K be the number of rows or columns of a distributed matrix,
*  and assume that its process grid has dimension p x q.
*  LOCr( K ) denotes the number of elements of K that a process
*  would receive if K were distributed over the p processes of its
*  process column.
*  Similarly, LOCc( K ) denotes the number of elements of K that a
*  process would receive if K were distributed over the q processes of
*  its process row.
*  The values of LOCr() and LOCc() may be determined via a call to the
*  ScaLAPACK tool function, NUMROC:
*          LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
*          LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
*  An upper bound for these quantities may be computed by:
*          LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
*          LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
*
*  Because vectors may be viewed as a subclass of matrices, a
*  distributed vector is considered to be a distributed matrix.
*
*  Arguments
*  =========
*
*  N       (global input) INTEGER
*          The global order of the elementary reflector. N >= 0.
*
*  ALPHA   (local output) DOUBLE PRECISION
*          On exit, alpha is computed in the process scope having the
*          vector sub( X ).
*
*  IAX     (global input) INTEGER
*          The global row index in X of X(IAX,JAX).
*
*  JAX     (global input) INTEGER
*          The global column index in X of X(IAX,JAX).
*
*  X       (local input/local output) DOUBLE PRECISION, pointer into the
*          local memory to an array of dimension (LLD_X,*). This array
*          contains the local pieces of the distributed vector sub( X ).
*          Before entry, the incremented array sub( X ) must contain
*          the vector x. On exit, it is overwritten with the vector v.
*
*  IX      (global input) INTEGER
*          The row index in the global array X indicating the first
*          row of sub( X ).
*
*  JX      (global input) INTEGER
*          The column index in the global array X indicating the
*          first column of sub( X ).
*
*  DESCX   (global and local input) INTEGER array of dimension DLEN_.
*          The array descriptor for the distributed matrix X.
*
*  INCX    (global input) INTEGER
*          The global increment for the elements of X. Only two values
*          of INCX are supported in this version, namely 1 and M_X.
*          INCX must not be zero.
*
*  TAU     (local output) DOUBLE PRECISION array, dimension  LOCc(JX)
*          if INCX = 1, and LOCr(IX) otherwise. This array contains the
*          Householder scalars related to the Householder vectors.
*          TAU is tied to the distributed matrix X.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,
     $                   LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,
     $                     CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,
     $                     RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            ICTXT, IIAX, INDXTAU, IXCOL, IXROW, J, JJAX,
     $                   KNT, MYCOL, MYROW, NPCOL, NPROW
      DOUBLE PRECISION   BETA, RSAFMN, SAFMIN, XNORM, EPS
*     ..
*     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, DGEBR2D, DGEBS2D, PDSCAL,
     $                   PDSCALCOPY, INFOG2L, PDNRM2
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY2
      EXTERNAL           DLAMCH, DLAPY2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN
*     ..
*     .. Executable Statements ..
*
      DEFF = 0
*
*     Get grid parameters.
*
      ICTXT = DESCX( CTXT_ )
      CALL BLACS_GRIDINFO( ICTXT, NPROW, NPCOL, MYROW, MYCOL )
*
      IF( INCX.EQ.DESCX( M_ ) ) THEN
*
*        sub( X ) is distributed across a process row.
*
         CALL INFOG2L( IX, JAX, DESCX, NPROW, NPCOL, MYROW, MYCOL,
     $                 IIAX, JJAX, IXROW, IXCOL )
*
         IF( MYROW.NE.IXROW )
     $      RETURN
*
*        Broadcast X(IAX,JAX) across the process row.
*
         IF( MYCOL.EQ.IXCOL ) THEN
            J = IIAX+(JJAX-1)*DESCX( LLD_ )
            CALL DGEBS2D( ICTXT, 'Rowwise', ' ', 1, 1, X( J ), 1 )
            ALPHA = X( J )
         ELSE
            CALL DGEBR2D( ICTXT, 'Rowwise', ' ', 1, 1, ALPHA, 1,
     $                    MYROW, IXCOL )
         END IF
*
         INDXTAU = IIAX - ND !FIXME
*
      ELSE
*
*        sub( X ) is distributed across a process column.
*
         CALL INFOG2L( IAX, JX, DESCX, NPROW, NPCOL, MYROW, MYCOL,
     $                 IIAX, JJAX, IXROW, IXCOL )
*
         IF( MYCOL.NE.IXCOL )
     $      RETURN
*
*        Broadcast X(IAX,JAX) across the process column.
*
         IF( MYROW.EQ.IXROW ) THEN
            J = IIAX+(JJAX-1)*DESCX( LLD_ )
            CALL DGEBS2D( ICTXT, 'Columnwise', ' ', 1, 1, X( J ), 1 )
            ALPHA = X( J )
         ELSE
            CALL DGEBR2D( ICTXT, 'Columnwise', ' ', 1, 1, ALPHA, 1,
     $                    IXROW, MYCOL )
         END IF
*
         INDXTAU = JJAX - ND !FIXME
*
      END IF
*
      IF( N.LE.0 ) THEN
******         TAU( INDXTAU ) = ZERO
         DEFF = 1
         RETURN
      END IF
*
      CALL PDNRM2( N-1, XNORM, X, IX, JX, DESCX, INCX )
*
      IF( XNORM.EQ.ZERO ) THEN
*
*        H = I
*
******         TAU( INDXTAU ) = ZERO
         DEFF = 1
*
      ELSE
*
*        General case
*
         BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
         EPS  = DLAMCH( 'E' )
*
         IF( ABS( BETA ).LT. 0.1 * EPS ) THEN
*
******             TAU( INDXTAU ) = ZERO
             DEFF = 1
*
         ELSE
*
             SAFMIN = DLAMCH( 'S' )
             RSAFMN = ONE / SAFMIN
*
             IF( ABS( BETA ).LT.SAFMIN ) THEN
*         
*               XNORM, BETA may be inaccurate; scale X and recompute them
*         
                KNT = 0
   10           CONTINUE
                KNT = KNT + 1
                CALL PDSCAL( N-1, RSAFMN, X, IX, JX, DESCX, INCX )
                BETA = BETA*RSAFMN
                ALPHA = ALPHA*RSAFMN
                IF( ABS( BETA ).LT.SAFMIN )
     $             GO TO 10
*         
*               New BETA is at most 1, at least SAFMIN
*         
                CALL PDNRM2( N-1, XNORM, X, IX, JX, DESCX, INCX )
                BETA = -SIGN( DLAPY2( ALPHA, XNORM ), ALPHA )
                TAU( INDXTAU ) = ( BETA-ALPHA ) / BETA
*FIXME use lines 27[7-8] instead of 2[79-82]
*                CALL PDSCALCOPY( N-1, ONE/(ALPHA-BETA), X, IX, JX,
*     $                          DESCX, INCX, X, IX, JX-ND, DESCX, INCX )
                CALL PDSCAL( N-1, ONE/(ALPHA-BETA), X, IX, JX,
     $                       DESCX, INCX )
                CALL PDCOPY( N-1, X, IX, JX,
     $                       DESCX, INCX, X, IX, JX-ND, DESCX, INCX )
*         
*               If ALPHA is subnormal, it may lose relative accuracy
*         
                ALPHA = BETA
                DO 20 J = 1, KNT
                   ALPHA = ALPHA*SAFMIN
   20           CONTINUE
             ELSE
                TAU( INDXTAU ) = ( BETA-ALPHA ) / BETA
*FIXME use lines 29[3-4] instead of 29[5-8]
*                CALL PDSCALCOPY( N-1, ONE/(ALPHA-BETA), X, IX, JX,
*     $                          DESCX, INCX, X, IX, JX-ND, DESCX, INCX )
                CALL PDSCAL( N-1, ONE/(ALPHA-BETA), X, IX, JX,
     $                       DESCX, INCX )
                CALL PDCOPY( N-1, X, IX, JX,
     $                       DESCX, INCX, X, IX, JX-ND, DESCX, INCX )
                ALPHA = BETA
             END IF
         END IF
      END IF
*
      RETURN
*
*     End of PDLARFG
*
      END
