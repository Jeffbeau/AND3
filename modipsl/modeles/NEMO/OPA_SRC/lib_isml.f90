!
! subroutines for PCG or SOR solvers
! (used if the ISML library is not available)
!
! linrg
!   gauss
!   vmov
!   desremopt
!   dtrsv
!   dger
!   xerbla
!   lsame
! folr (empty)
!
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/lib_isml.f90,v 1.4 2006/03/09 17:21:50 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
!---------------------------------------------------------
   SUBROUTINE linrg(kn,pa,klda,painv,kldainv)

      !! compute inverse matrix

      IMPLICIT NONE 
      INTEGER kn,klda,kldainv
      REAL (kind=8) ::   pa(kn,kn),painv(kn,kn)
      REAL (kind=8) ::   zb(kn,kn)
      REAL (kind=8) ::   zv(kn)
      INTEGER iplin(kn)
      INTEGER ji

      IF( kn /= klda .OR. kn /= kldainv ) THEN 
          write(0,*)'change your parameters'
          STOP 
      ENDIF 

      CALL vmov( kn*kn, pa, painv )

      CALL gauss( kn, painv, iplin, zv )

      zb(:,:) = 0.e0   
      DO ji = 1, kn
        zb(ji,ji) = 1.e0
        CALL desremopt( kn, painv, iplin, zb(1,ji), zb(1,ji), zv )
      END DO 
      CALL vmov( kn*kn, zb, painv )

   END SUBROUTINE linrg
!---------------------------------------------------------
   SUBROUTINE gauss(kn,pa,kplin,pv)

      IMPLICIT NONE 
      INTEGER kn
      INTEGER ji,jj,jk
      INTEGER ik,ipp
      REAL (kind=8) ::   pa(kn,kn),pv(kn)
!!    REAL (kind=8) ::   zpivmax,zalpha
      REAL (kind=8) ::   zalpha
      INTEGER kplin(kn)
      INTEGER isamax
      EXTERNAL isamax

!  factorisation de Gauss de la matrice a avec pivot partiel . 
!  initialisation des pointeurs .
      DO ji=1,kn
        kplin(ji)=ji
      END DO 
      DO jk=1,kn-1
!  recherche du pivot maximal .
!!      ik=jk
!!      zpivmax=dabs(pa(jk,jk))
!!        DO ji=jk,kn
!!          IF(dabs(pa(ji,jk)) > zpivmax) THEN 
!!            zpivmax=dabs(pa(ji,jk))
!!            ik=ji
!!          ENDIF 
!!        END DO 
        ik=isamax( kn-jk+1, pa(jk,jk) )+jk-1
!  permutation de la ligne jk et de la ligne ik .
        IF(jk == 58) THEN
            PRINT *,'matrix ',(pa(jk,ji),ji=1,kn)
            PRINT *,' pivot ',ik,kplin(ik),kplin(jk)
        ENDIF 
        ipp=kplin(ik)
        kplin(ik)=kplin(jk)
        kplin(jk)=ipp
        DO jj=1,kn
          pv(jj)=pa(ik,jj)
          pa(ik,jj)=pa(jk,jj)
          pa(jk,jj)=pv(jj)
        END DO 
        IF(jk == 58) THEN
            PRINT *,'matrix ',(pa(jk,ji),ji=1,kn)
            PRINT *,' pivot ',ik,kplin(ik),kplin(jk)
        ENDIF 
!  calcul des coefficients de la colonne k ligne a ligne .
        DO ji=jk+1,kn
          IF(pa(jk,jk) == 0) THEN 
              PRINT *,'probleme diagonale nulle',jk,pa(jk,jk)
              pa(ji,jk)=pa(ji,jk)/1.E-20
          ENDIF 
          IF(pa(jk,jk) /= 0) pa(ji,jk)=pa(ji,jk)/pa(jk,jk)
        END DO 
!!        DO ji=jk+1,kn
!!          DO jj=jk+1,kn
!!            pa(ji,jj)=pa(ji,jj)-pa(ji,jk)*pa(jk,jj)
!!          END DO 
!!        END DO 
        zalpha=-1.
        CALL dger(kn-jk,kn-jk,zalpha,pa(jk+1,jk),1,pa(jk,jk+1),kn,   &
            pa(jk+1,jk+1),kn)
      END DO 

   END SUBROUTINE gauss
!---------------------------------------------------------
   FUNCTION isamax( I, X )
      DIMENSION X(I)
      ISAMAX = 0
      XMIN = -huge(1.)
      DO N = 1, I
         IF(ABS(X(N)) > XMIN ) THEN
            XMIN = X(N)
            ISAMAX = N
         ENDIF
      END DO
   END FUNCTION isamax
!---------------------------------------------------------
   SUBROUTINE vmov(kn,px,py)

      IMPLICIT NONE 
      INTEGER kn
      REAL (kind=8) ::   px(kn),py(kn)
      INTEGER ji

      DO ji=1,kn
         py(ji)=px(ji)
      END DO 

   END SUBROUTINE vmov
!---------------------------------------------------------
   subroutine desremopt(n,a,plin,y,x,v)
      implicit none
      integer n,i,  j0
!!    integer n,i,j,j0
      real (kind=8) ::   a(n,n),x(n),y(n),v(n)
      integer plin(n)
!  descente remontee du systeme .
!  initialisation du vecteur resultat .
!  prise en compte de la permutation des lignes .
      do i=1,n
        v(i)=y(plin(i))
      end do
      do i=1,n
        if(v(i) /= 0.) then
          j0=i-1
          goto 1
        endif
      end do
1     continue
!  descente du systeme L v = v , L est a diagonale unitaire .
!!      do j=j0+1,n
!!        do i=j+1,n
!!          v(i)=v(i)-a(i,j)*v(j)
!!        end do
!!      end do
      call dtrsv('L','N','U',n-j0,a(j0+1,j0+1),n,v(j0+1),1)
!  remontee du systeme U v = v . 
!!      do j=n,1,-1
!!        v(j)=v(j)/a(j,j)
!!        do i=1,j-1
!!          v(i)=v(i)-a(i,j)*v(j)
!!        end do
!!      end do          
      call dtrsv('U','N','N',n,a,n,v,1)
!  prise en compte de la permutation des colonnes .
      do i=1,n
        x(i)=v(i)
      end do

   end SUBROUTINE desremopt
!---------------------------------------------------------
   SUBROUTINE DTRSV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
!!    .. Scalar Arguments ..
      INTEGER            INCX, LDA, N
      CHARACTER (len=1) ::   DIAG, TRANS, UPLO
!!    .. Array Arguments ..
!     DOUBLE PRECISION   A( LDA, * ), X( * )
      REAL (kind=8) ::     A( LDA, * ), X( * )
!!    ..

!! Purpose
!! =======

!! DTRSV  solves one of the systems of equations

!!    A*x = b,   or   A'*x = b,

!! where b and x are n element vectors and A is an n by n unit, or
!! non-unit, upper or lower triangular matrix.

!! No test for singularity or near-singularity is included in this
!! routine. Such tests must be performed before calling this routine.

!! Parameters
!! ==========

!! UPLO   - CHARACTER*1.
!!          On entry, UPLO specifies whether the matrix is an upper or
!!          lower triangular matrix as follows:

!!             UPLO = 'U' or 'u'   A is an upper triangular matrix.

!!             UPLO = 'L' or 'l'   A is a lower triangular matrix.

!!          Unchanged on exit.

!! TRANS  - CHARACTER*1.
!!          On entry, TRANS specifies the equations to be solved as
!!          follows:

!!             TRANS = 'N' or 'n'   A*x = b.

!!             TRANS = 'T' or 't'   A'*x = b.

!!             TRANS = 'C' or 'c'   A'*x = b.

!!          Unchanged on exit.

!! DIAG   - CHARACTER*1.
!!          On entry, DIAG specifies whether or not A is unit
!!          triangular as follows:

!!             DIAG = 'U' or 'u'   A is assumed to be unit triangular.

!!             DIAG = 'N' or 'n'   A is not assumed to be unit
!!                                 triangular.

!!          Unchanged on exit.

!! N      - INTEGER.
!!          On entry, N specifies the order of the matrix A.
!!          N must be at least zero.
!!          Unchanged on exit.

!! A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!!          Before entry with  UPLO = 'U' or 'u', the leading n by n
!!          upper triangular part of the array A must contain the upper
!!          triangular matrix and the strictly lower triangular part of
!!          A is not referenced.
!!          Before entry with UPLO = 'L' or 'l', the leading n by n
!!          lower triangular part of the array A must contain the lower
!!          triangular matrix and the strictly upper triangular part of
!!          A is not referenced.
!!          Note that when  DIAG = 'U' or 'u', the diagonal elements of
!!          A are not referenced either, but are assumed to be unity.
!!          Unchanged on exit.

!! LDA    - INTEGER.
!!          On entry, LDA specifies the first dimension of A as declared
!!          in the calling (sub) program. LDA must be at least
!!          max( 1, n ).
!!          Unchanged on exit.

!! X      - DOUBLE PRECISION array of dimension at least
!!          ( 1 + ( n - 1 )*abs( INCX ) ).
!!          Before entry, the incremented array X must contain the n
!!          element right-hand side vector b. On exit, X is overwritten
!!          with the solution vector x.

!! INCX   - INTEGER.
!!          On entry, INCX specifies the increment for the elements of
!!          X. INCX must not be zero.
!!          Unchanged on exit.


!! Level 2 Blas routine.

!! -- Written on 22-October-1986.
!!    Jack Dongarra, Argonne National Lab.
!!    Jeremy Du Croz, Nag Central Office.
!!    Sven Hammarling, Nag Central Office.
!!    Richard Hanson, Sandia National Labs.


!!    .. Parameters ..
!     DOUBLE PRECISION   ZERO
!     PARAMETER        ( ZERO = 0.0D+0 )
      REAL (kind=8) ::                 ZERO
      PARAMETER        ( ZERO = 0.0 )
!!    .. Local Scalars ..
!     DOUBLE PRECISION   TEMP
      REAL (kind=8) ::                 TEMP
      INTEGER            I, INFO, IX, J, JX, KX
      LOGICAL            NOUNIT
!!    .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
!!    .. External Subroutines ..
      EXTERNAL           XERBLA
!!    .. Intrinsic Functions ..
      INTRINSIC          MAX
!!    ..
!!    .. Executable Statements ..

!!    Test the input parameters.

      INFO = 0
      IF     ( .NOT.LSAME( UPLO , 'U' ).AND.   &
               .NOT.LSAME( UPLO , 'L' )      )THEN
         INFO = 1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ).AND.   &
               .NOT.LSAME( TRANS, 'T' ).AND.   &
               .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 2
      ELSE IF( .NOT.LSAME( DIAG , 'U' ).AND.   &
               .NOT.LSAME( DIAG , 'N' )      )THEN
         INFO = 3
      ELSE IF( N < 0 )THEN
         INFO = 4
      ELSE IF( LDA < MAX( 1, N ) )THEN
         INFO = 6
      ELSE IF( INCX == 0 )THEN
         INFO = 8
      END IF
      IF( INFO /= 0 )THEN
         CALL XERBLA( 'DTRSV ', INFO )
         RETURN
      END IF

!!    Quick return if possible.

      IF( N == 0 )   RETURN

      NOUNIT = LSAME( DIAG, 'N' )

!!    Set up the start point in X if the increment is not unity. This
!!    will be  ( N - 1 )*INCX  too small for descending loops.

      IF( INCX <= 0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX /= 1 )THEN
         KX = 1
      END IF

!!    Start the operations. In this version the elements of A are
!!    accessed sequentially with one pass through A.

      IF( LSAME( TRANS, 'N' ) )THEN

!!       Form  x := inv( A )*x.

         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX == 1 )THEN
               DO 20, J = N, 1, -1
                  IF( X( J ) /= ZERO )THEN
                     IF( NOUNIT )   X( J ) = X( J )/A( J, J )
                     TEMP = X( J )
                     DO 10, I = J - 1, 1, -1
                        X( I ) = X( I ) - TEMP*A( I, J )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO 40, J = N, 1, -1
                  IF( X( JX ) /= ZERO )THEN
                     IF( NOUNIT )  X( JX ) = X( JX )/A( J, J )
                     TEMP = X( JX )
                     IX   = JX
                     DO 30, I = J - 1, 1, -1
                        IX      = IX      - INCX
                        X( IX ) = X( IX ) - TEMP*A( I, J )
   30                CONTINUE
                  END IF
                  JX = JX - INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX == 1 )THEN
               DO 60, J = 1, N
                  IF( X( J ) /= ZERO )THEN
                     IF( NOUNIT )   X( J ) = X( J )/A( J, J )
                     TEMP = X( J )
                     DO 50, I = J + 1, N
                        X( I ) = X( I ) - TEMP*A( I, J )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE
               JX = KX
               DO 80, J = 1, N
                  IF( X( JX ) /= ZERO )THEN
                     IF( NOUNIT )   X( JX ) = X( JX )/A( J, J )
                     TEMP = X( JX )
                     IX   = JX
                     DO 70, I = J + 1, N
                        IX      = IX      + INCX
                        X( IX ) = X( IX ) - TEMP*A( I, J )
   70                CONTINUE
                  END IF
                  JX = JX + INCX
   80          CONTINUE
            END IF
         END IF
      ELSE

!!       Form  x := inv( A' )*x.

         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX == 1 )THEN
               DO 100, J = 1, N
                  TEMP = X( J )
                  DO 90, I = 1, J - 1
                     TEMP = TEMP - A( I, J )*X( I )
   90             CONTINUE
                  IF( NOUNIT )   TEMP = TEMP/A( J, J )
                  X( J ) = TEMP
  100          CONTINUE
            ELSE
               JX = KX
               DO 120, J = 1, N
                  TEMP = X( JX )
                  IX   = KX
                  DO 110, I = 1, J - 1
                     TEMP = TEMP - A( I, J )*X( IX )
                     IX   = IX   + INCX
  110             CONTINUE
                  IF( NOUNIT )   TEMP = TEMP/A( J, J )
                  X( JX ) = TEMP
                  JX      = JX   + INCX
  120          CONTINUE
            END IF
         ELSE
            IF( INCX == 1 )THEN
               DO 140, J = N, 1, -1
                  TEMP = X( J )
                  DO 130, I = N, J + 1, -1
                     TEMP = TEMP - A( I, J )*X( I )
  130             CONTINUE
                  IF( NOUNIT )    TEMP = TEMP/A( J, J )
                  X( J ) = TEMP
  140          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 160, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = KX
                  DO 150, I = N, J + 1, -1
                     TEMP = TEMP - A( I, J )*X( IX )
                     IX   = IX   - INCX
  150             CONTINUE
                  IF( NOUNIT )    TEMP = TEMP/A( J, J )
                  X( JX ) = TEMP
                  JX      = JX   - INCX
  160          CONTINUE
            END IF
         END IF
      END IF

   END SUBROUTINE DTRSV
!---------------------------------------------------------
   SUBROUTINE DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
!!    .. Scalar Arguments ..
!     DOUBLE PRECISION   ALPHA
      REAL (kind=8) ::                 ALPHA
      INTEGER            INCX, INCY, LDA, M, N
!!    .. Array Arguments ..
!     DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
      REAL (kind=8) ::   A( LDA, * ), X( * ), Y( * )
!!    ..

!! Purpose
!! =======

!! DGER   performs the rank 1 operation

!!    A := alpha*x*y' + A,

!! where alpha is a scalar, x is an m element vector, y is an n element
!! vector and A is an m by n matrix.

!! Parameters
!! ==========

!! M      - INTEGER.
!!          On entry, M specifies the number of rows of the matrix A.
!!          M must be at least zero.
!!          Unchanged on exit.

!! N      - INTEGER.
!!          On entry, N specifies the number of columns of the matrix A.
!!          N must be at least zero.
!!          Unchanged on exit.

!! ALPHA  - DOUBLE PRECISION.
!!          On entry, ALPHA specifies the scalar alpha.
!!          Unchanged on exit.

!! X      - DOUBLE PRECISION array of dimension at least
!!          ( 1 + ( m - 1 )*abs( INCX ) ).
!!          Before entry, the incremented array X must contain the m
!!          element vector x.
!!          Unchanged on exit.

!! INCX   - INTEGER.
!!          On entry, INCX specifies the increment for the elements of
!!          X. INCX must not be zero.
!!          Unchanged on exit.

!! Y      - DOUBLE PRECISION array of dimension at least
!!          ( 1 + ( n - 1 )*abs( INCY ) ).
!!          Before entry, the incremented array Y must contain the n
!!          element vector y.
!!          Unchanged on exit.

!! INCY   - INTEGER.
!!          On entry, INCY specifies the increment for the elements of
!!          Y. INCY must not be zero.
!!          Unchanged on exit.

!! A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!!          Before entry, the leading m by n part of the array A must
!!          contain the matrix of coefficients. On exit, A is
!!          overwritten by the updated matrix.

!! LDA    - INTEGER.
!!          On entry, LDA specifies the first dimension of A as declared
!!          in the calling (sub) program. LDA must be at least
!!          max( 1, m ).
!!          Unchanged on exit.


!! Level 2 Blas routine.

!! -- Written on 22-October-1986.
!!    Jack Dongarra, Argonne National Lab.
!!    Jeremy Du Croz, Nag Central Office.
!!    Sven Hammarling, Nag Central Office.
!!    Richard Hanson, Sandia National Labs.


!!    .. Parameters ..
!     DOUBLE PRECISION   ZERO
!     PARAMETER        ( ZERO = 0.0D+0 )
      REAL (kind=8) ::                 ZERO
      PARAMETER        ( ZERO = 0.0 )
!!    .. Local Scalars ..
!     DOUBLE PRECISION   TEMP
      REAL (kind=8) ::                 TEMP
      INTEGER            I, INFO, IX, J, JY, KX
!!    .. External Subroutines ..
      EXTERNAL           XERBLA
!!    .. Intrinsic Functions ..
      INTRINSIC          MAX
!!    ..
!!    .. Executable Statements ..

!!    Test the input parameters.

      INFO = 0
      IF     ( M < 0 )THEN
         INFO = 1
      ELSE IF( N < 0 )THEN
         INFO = 2
      ELSE IF( INCX == 0 )THEN
         INFO = 5
      ELSE IF( INCY == 0 )THEN
         INFO = 7
      ELSE IF( LDA < MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO /= 0 )THEN
         CALL XERBLA( 'DGER  ', INFO )
         RETURN
      END IF

!!    Quick return if possible.

      IF( ( M == 0 ).OR.( N == 0 ).OR.( ALPHA == ZERO ) )   &
         RETURN

!!    Start the operations. In this version the elements of A are
!!    accessed sequentially with one pass through A.

      IF( INCY > 0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX == 1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ) /= ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX > 0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ) /= ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF

   END SUBROUTINE DGER
!---------------------------------------------------------
   SUBROUTINE XERBLA ( SRNAME, INFO )
!!    ..    Scalar Arguments ..
      INTEGER            INFO
      CHARACTER (len=6) ::    SRNAME
!!    ..

!! Purpose
!! =======

!! XERBLA  is an error handler for the Level 2 BLAS routines.

!! It is called by the Level 2 BLAS routines if an input parameter is
!! invalid.

!! Installers should consider modifying the STOP statement in order to
!! call system-specific exception-handling facilities.

!! Parameters
!! ==========

!! SRNAME - CHARACTER*6.
!!          On entry, SRNAME specifies the name of the routine which
!!          called XERBLA.

!! INFO   - INTEGER.
!!          On entry, INFO specifies the position of the invalid
!!          parameter in the parameter-list of the calling routine.


!! Auxiliary routine for Level 2 Blas.

!! Written on 20-July-1986.

!!    .. Executable Statements ..

      WRITE (*,99999) SRNAME, INFO

      STOP

99999 FORMAT ( ' ** On entry to ', A6, ' parameter number ', I2,   &
               ' had an illegal value' )

   END SUBROUTINE XERBLA
!-----------------------------------------------------------
   FUNCTION lsame( c1, c2 )
      logical lsame
      CHARACTER (len=*), INTENT(in) ::   c1, c2
      IF( c1 == c2 ) THEN
          lsame=.TRUE.
      ELSE
          lsame=.FALSE.
      ENDIF
   END FUNCTION lsame
