


PGF90 (Version     10.5)          09/23/2010  18:36:34      page 1

Switches: -noasm -nodclchk -nodebug -nodlines -noline -list
          -idir ../../../lib
          -idir ../../../lib/oce
          -idir /distrib/local/netcdf/pgf/include/
          -inform severe -opt 3 -nosave -object -noonetrip
          -depchk on -nostandard     
          -nosymbol -noupcase    

Filename: lib_isml.f90

(    1) !
(    2) ! subroutines for PCG or SOR solvers
(    3) ! (used if the ISML library is not available)
(    4) !
(    5) ! linrg
(    6) !   gauss
(    7) !   vmov
(    8) !   desremopt
(    9) !   dtrsv
(   10) !   dger
(   11) !   xerbla
(   12) !   lsame
(   13) ! folr (empty)
(   14) !
(   15)    !!----------------------------------------------------------------------
(   16)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(   17)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/lib_isml.f90,v 1.4 2006/03/09 17:21:50 opalod Exp $ 
(   18)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   19)    !!----------------------------------------------------------------------
(   20) !---------------------------------------------------------
(   21)    SUBROUTINE linrg(kn,pa,klda,painv,kldainv)
(   22) 
(   23)       !! compute inverse matrix
(   24) 
(   25)       IMPLICIT NONE 
(   26)       INTEGER kn,klda,kldainv
(   27)       REAL (kind=8) ::   pa(kn,kn),painv(kn,kn)
(   28)       REAL (kind=8) ::   zb(kn,kn)
(   29)       REAL (kind=8) ::   zv(kn)
(   30)       INTEGER iplin(kn)
(   31)       INTEGER ji
(   32) 
(   33)       IF( kn /= klda .OR. kn /= kldainv ) THEN 
(   34)           write(0,*)'change your parameters'
(   35)           STOP 
(   36)       ENDIF 
(   37) 
(   38)       CALL vmov( kn*kn, pa, painv )
(   39) 
(   40)       CALL gauss( kn, painv, iplin, zv )
(   41) 
(   42)       zb(:,:) = 0.e0   
(   43)       DO ji = 1, kn
(   44)         zb(ji,ji) = 1.e0
(   45)         CALL desremopt( kn, painv, iplin, zb(1,ji), zb(1,ji), zv )
(   46)       END DO 
(   47)       CALL vmov( kn*kn, zb, painv )
(   48) 
(   49)    END SUBROUTINE linrg






PGF90 (Version     10.5)          09/23/2010  18:36:34      page 2

(   50) !---------------------------------------------------------
(   51)    SUBROUTINE gauss(kn,pa,kplin,pv)
(   52) 
(   53)       IMPLICIT NONE 
(   54)       INTEGER kn
(   55)       INTEGER ji,jj,jk
(   56)       INTEGER ik,ipp
(   57)       REAL (kind=8) ::   pa(kn,kn),pv(kn)
(   58) !!    REAL (kind=8) ::   zpivmax,zalpha
(   59)       REAL (kind=8) ::   zalpha
(   60)       INTEGER kplin(kn)
(   61)       INTEGER isamax
(   62)       EXTERNAL isamax
(   63) 
(   64) !  factorisation de Gauss de la matrice a avec pivot partiel . 
(   65) !  initialisation des pointeurs .
(   66)       DO ji=1,kn
(   67)         kplin(ji)=ji
(   68)       END DO 
(   69)       DO jk=1,kn-1
(   70) !  recherche du pivot maximal .
(   71) !!      ik=jk
(   72) !!      zpivmax=dabs(pa(jk,jk))
(   73) !!        DO ji=jk,kn
(   74) !!          IF(dabs(pa(ji,jk)) > zpivmax) THEN 
(   75) !!            zpivmax=dabs(pa(ji,jk))
(   76) !!            ik=ji
(   77) !!          ENDIF 
(   78) !!        END DO 
(   79)         ik=isamax( kn-jk+1, pa(jk,jk) )+jk-1
(   80) !  permutation de la ligne jk et de la ligne ik .
(   81)         IF(jk == 58) THEN
(   82)             PRINT *,'matrix ',(pa(jk,ji),ji=1,kn)
(   83)             PRINT *,' pivot ',ik,kplin(ik),kplin(jk)
(   84)         ENDIF 
(   85)         ipp=kplin(ik)
(   86)         kplin(ik)=kplin(jk)
(   87)         kplin(jk)=ipp
(   88)         DO jj=1,kn
(   89)           pv(jj)=pa(ik,jj)
(   90)           pa(ik,jj)=pa(jk,jj)
(   91)           pa(jk,jj)=pv(jj)
(   92)         END DO 
(   93)         IF(jk == 58) THEN
(   94)             PRINT *,'matrix ',(pa(jk,ji),ji=1,kn)
(   95)             PRINT *,' pivot ',ik,kplin(ik),kplin(jk)
(   96)         ENDIF 
(   97) !  calcul des coefficients de la colonne k ligne a ligne .
(   98)         DO ji=jk+1,kn
(   99)           IF(pa(jk,jk) == 0) THEN 
(  100)               PRINT *,'probleme diagonale nulle',jk,pa(jk,jk)
(  101)               pa(ji,jk)=pa(ji,jk)/1.E-20
(  102)           ENDIF 
(  103)           IF(pa(jk,jk) /= 0) pa(ji,jk)=pa(ji,jk)/pa(jk,jk)
(  104)         END DO 
(  105) !!        DO ji=jk+1,kn
(  106) !!          DO jj=jk+1,kn
(  107) !!            pa(ji,jj)=pa(ji,jj)-pa(ji,jk)*pa(jk,jj)






PGF90 (Version     10.5)          09/23/2010  18:36:34      page 3

(  108) !!          END DO 
(  109) !!        END DO 
(  110)         zalpha=-1.
(  111)         CALL dger(kn-jk,kn-jk,zalpha,pa(jk+1,jk),1,pa(jk,jk+1),kn,   &
(  112)             pa(jk+1,jk+1),kn)
(  113)       END DO 
(  114) 
(  115)    END SUBROUTINE gauss
(  116) !---------------------------------------------------------
(  117)    FUNCTION isamax( I, X )
(  118)       DIMENSION X(I)
(  119)       ISAMAX = 0
(  120)       XMIN = -huge(1.)
(  121)       DO N = 1, I
(  122)          IF(ABS(X(N)) > XMIN ) THEN
(  123)             XMIN = X(N)
(  124)             ISAMAX = N
(  125)          ENDIF
(  126)       END DO
(  127)    END FUNCTION isamax
(  128) !---------------------------------------------------------
(  129)    SUBROUTINE vmov(kn,px,py)
(  130) 
(  131)       IMPLICIT NONE 
(  132)       INTEGER kn
(  133)       REAL (kind=8) ::   px(kn),py(kn)
(  134)       INTEGER ji
(  135) 
(  136)       DO ji=1,kn
(  137)          py(ji)=px(ji)
(  138)       END DO 
(  139) 
(  140)    END SUBROUTINE vmov
(  141) !---------------------------------------------------------
(  142)    subroutine desremopt(n,a,plin,y,x,v)
(  143)       implicit none
(  144)       integer n,i,  j0
(  145) !!    integer n,i,j,j0
(  146)       real (kind=8) ::   a(n,n),x(n),y(n),v(n)
(  147)       integer plin(n)
(  148) !  descente remontee du systeme .
(  149) !  initialisation du vecteur resultat .
(  150) !  prise en compte de la permutation des lignes .
(  151)       do i=1,n
(  152)         v(i)=y(plin(i))
(  153)       end do
(  154)       do i=1,n
(  155)         if(v(i) /= 0.) then
(  156)           j0=i-1
(  157)           goto 1
(  158)         endif
(  159)       end do
(  160) 1     continue
(  161) !  descente du systeme L v = v , L est a diagonale unitaire .
(  162) !!      do j=j0+1,n
(  163) !!        do i=j+1,n
(  164) !!          v(i)=v(i)-a(i,j)*v(j)
(  165) !!        end do






PGF90 (Version     10.5)          09/23/2010  18:36:34      page 4

(  166) !!      end do
(  167)       call dtrsv('L','N','U',n-j0,a(j0+1,j0+1),n,v(j0+1),1)
(  168) !  remontee du systeme U v = v . 
(  169) !!      do j=n,1,-1
(  170) !!        v(j)=v(j)/a(j,j)
(  171) !!        do i=1,j-1
(  172) !!          v(i)=v(i)-a(i,j)*v(j)
(  173) !!        end do
(  174) !!      end do          
(  175)       call dtrsv('U','N','N',n,a,n,v,1)
(  176) !  prise en compte de la permutation des colonnes .
(  177)       do i=1,n
(  178)         x(i)=v(i)
(  179)       end do
(  180) 
(  181)    end SUBROUTINE desremopt
(  182) !---------------------------------------------------------
(  183)    SUBROUTINE DTRSV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
(  184) !!    .. Scalar Arguments ..
(  185)       INTEGER            INCX, LDA, N
(  186)       CHARACTER (len=1) ::   DIAG, TRANS, UPLO
(  187) !!    .. Array Arguments ..
(  188) !     DOUBLE PRECISION   A( LDA, * ), X( * )
(  189)       REAL (kind=8) ::     A( LDA, * ), X( * )
(  190) !!    ..
(  191) 
(  192) !! Purpose
(  193) !! =======
(  194) 
(  195) !! DTRSV  solves one of the systems of equations
(  196) 
(  197) !!    A*x = b,   or   A'*x = b,
(  198) 
(  199) !! where b and x are n element vectors and A is an n by n unit, or
(  200) !! non-unit, upper or lower triangular matrix.
(  201) 
(  202) !! No test for singularity or near-singularity is included in this
(  203) !! routine. Such tests must be performed before calling this routine.
(  204) 
(  205) !! Parameters
(  206) !! ==========
(  207) 
(  208) !! UPLO   - CHARACTER*1.
(  209) !!          On entry, UPLO specifies whether the matrix is an upper or
(  210) !!          lower triangular matrix as follows:
(  211) 
(  212) !!             UPLO = 'U' or 'u'   A is an upper triangular matrix.
(  213) 
(  214) !!             UPLO = 'L' or 'l'   A is a lower triangular matrix.
(  215) 
(  216) !!          Unchanged on exit.
(  217) 
(  218) !! TRANS  - CHARACTER*1.
(  219) !!          On entry, TRANS specifies the equations to be solved as
(  220) !!          follows:
(  221) 
(  222) !!             TRANS = 'N' or 'n'   A*x = b.
(  223) 






PGF90 (Version     10.5)          09/23/2010  18:36:34      page 5

(  224) !!             TRANS = 'T' or 't'   A'*x = b.
(  225) 
(  226) !!             TRANS = 'C' or 'c'   A'*x = b.
(  227) 
(  228) !!          Unchanged on exit.
(  229) 
(  230) !! DIAG   - CHARACTER*1.
(  231) !!          On entry, DIAG specifies whether or not A is unit
(  232) !!          triangular as follows:
(  233) 
(  234) !!             DIAG = 'U' or 'u'   A is assumed to be unit triangular.
(  235) 
(  236) !!             DIAG = 'N' or 'n'   A is not assumed to be unit
(  237) !!                                 triangular.
(  238) 
(  239) !!          Unchanged on exit.
(  240) 
(  241) !! N      - INTEGER.
(  242) !!          On entry, N specifies the order of the matrix A.
(  243) !!          N must be at least zero.
(  244) !!          Unchanged on exit.
(  245) 
(  246) !! A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
(  247) !!          Before entry with  UPLO = 'U' or 'u', the leading n by n
(  248) !!          upper triangular part of the array A must contain the upper
(  249) !!          triangular matrix and the strictly lower triangular part of
(  250) !!          A is not referenced.
(  251) !!          Before entry with UPLO = 'L' or 'l', the leading n by n
(  252) !!          lower triangular part of the array A must contain the lower
(  253) !!          triangular matrix and the strictly upper triangular part of
(  254) !!          A is not referenced.
(  255) !!          Note that when  DIAG = 'U' or 'u', the diagonal elements of
(  256) !!          A are not referenced either, but are assumed to be unity.
(  257) !!          Unchanged on exit.
(  258) 
(  259) !! LDA    - INTEGER.
(  260) !!          On entry, LDA specifies the first dimension of A as declared
(  261) !!          in the calling (sub) program. LDA must be at least
(  262) !!          max( 1, n ).
(  263) !!          Unchanged on exit.
(  264) 
(  265) !! X      - DOUBLE PRECISION array of dimension at least
(  266) !!          ( 1 + ( n - 1 )*abs( INCX ) ).
(  267) !!          Before entry, the incremented array X must contain the n
(  268) !!          element right-hand side vector b. On exit, X is overwritten
(  269) !!          with the solution vector x.
(  270) 
(  271) !! INCX   - INTEGER.
(  272) !!          On entry, INCX specifies the increment for the elements of
(  273) !!          X. INCX must not be zero.
(  274) !!          Unchanged on exit.
(  275) 
(  276) 
(  277) !! Level 2 Blas routine.
(  278) 
(  279) !! -- Written on 22-October-1986.
(  280) !!    Jack Dongarra, Argonne National Lab.
(  281) !!    Jeremy Du Croz, Nag Central Office.






PGF90 (Version     10.5)          09/23/2010  18:36:34      page 6

(  282) !!    Sven Hammarling, Nag Central Office.
(  283) !!    Richard Hanson, Sandia National Labs.
(  284) 
(  285) 
(  286) !!    .. Parameters ..
(  287) !     DOUBLE PRECISION   ZERO
(  288) !     PARAMETER        ( ZERO = 0.0D+0 )
(  289)       REAL (kind=8) ::                 ZERO
(  290)       PARAMETER        ( ZERO = 0.0 )
(  291) !!    .. Local Scalars ..
(  292) !     DOUBLE PRECISION   TEMP
(  293)       REAL (kind=8) ::                 TEMP
(  294)       INTEGER            I, INFO, IX, J, JX, KX
(  295)       LOGICAL            NOUNIT
(  296) !!    .. External Functions ..
(  297)       LOGICAL            LSAME
(  298)       EXTERNAL           LSAME
(  299) !!    .. External Subroutines ..
(  300)       EXTERNAL           XERBLA
(  301) !!    .. Intrinsic Functions ..
(  302)       INTRINSIC          MAX
(  303) !!    ..
(  304) !!    .. Executable Statements ..
(  305) 
(  306) !!    Test the input parameters.
(  307) 
(  308)       INFO = 0
(  309)       IF     ( .NOT.LSAME( UPLO , 'U' ).AND.   &
(  310)                .NOT.LSAME( UPLO , 'L' )      )THEN
(  311)          INFO = 1
(  312)       ELSE IF( .NOT.LSAME( TRANS, 'N' ).AND.   &
(  313)                .NOT.LSAME( TRANS, 'T' ).AND.   &
(  314)                .NOT.LSAME( TRANS, 'C' )      )THEN
(  315)          INFO = 2
(  316)       ELSE IF( .NOT.LSAME( DIAG , 'U' ).AND.   &
(  317)                .NOT.LSAME( DIAG , 'N' )      )THEN
(  318)          INFO = 3
(  319)       ELSE IF( N < 0 )THEN
(  320)          INFO = 4
(  321)       ELSE IF( LDA < MAX( 1, N ) )THEN
(  322)          INFO = 6
(  323)       ELSE IF( INCX == 0 )THEN
(  324)          INFO = 8
(  325)       END IF
(  326)       IF( INFO /= 0 )THEN
(  327)          CALL XERBLA( 'DTRSV ', INFO )
(  328)          RETURN
(  329)       END IF
(  330) 
(  331) !!    Quick return if possible.
(  332) 
(  333)       IF( N == 0 )   RETURN
(  334) 
(  335)       NOUNIT = LSAME( DIAG, 'N' )
(  336) 
(  337) !!    Set up the start point in X if the increment is not unity. This
(  338) !!    will be  ( N - 1 )*INCX  too small for descending loops.
(  339) 






PGF90 (Version     10.5)          09/23/2010  18:36:34      page 7

(  340)       IF( INCX <= 0 )THEN
(  341)          KX = 1 - ( N - 1 )*INCX
(  342)       ELSE IF( INCX /= 1 )THEN
(  343)          KX = 1
(  344)       END IF
(  345) 
(  346) !!    Start the operations. In this version the elements of A are
(  347) !!    accessed sequentially with one pass through A.
(  348) 
(  349)       IF( LSAME( TRANS, 'N' ) )THEN
(  350) 
(  351) !!       Form  x := inv( A )*x.
(  352) 
(  353)          IF( LSAME( UPLO, 'U' ) )THEN
(  354)             IF( INCX == 1 )THEN
(  355)                DO 20, J = N, 1, -1
(  356)                   IF( X( J ) /= ZERO )THEN
(  357)                      IF( NOUNIT )   X( J ) = X( J )/A( J, J )
(  358)                      TEMP = X( J )
(  359)                      DO 10, I = J - 1, 1, -1
(  360)                         X( I ) = X( I ) - TEMP*A( I, J )
(  361)    10                CONTINUE
(  362)                   END IF
(  363)    20          CONTINUE
(  364)             ELSE
(  365)                JX = KX + ( N - 1 )*INCX
(  366)                DO 40, J = N, 1, -1
(  367)                   IF( X( JX ) /= ZERO )THEN
(  368)                      IF( NOUNIT )  X( JX ) = X( JX )/A( J, J )
(  369)                      TEMP = X( JX )
(  370)                      IX   = JX
(  371)                      DO 30, I = J - 1, 1, -1
(  372)                         IX      = IX      - INCX
(  373)                         X( IX ) = X( IX ) - TEMP*A( I, J )
(  374)    30                CONTINUE
(  375)                   END IF
(  376)                   JX = JX - INCX
(  377)    40          CONTINUE
(  378)             END IF
(  379)          ELSE
(  380)             IF( INCX == 1 )THEN
(  381)                DO 60, J = 1, N
(  382)                   IF( X( J ) /= ZERO )THEN
(  383)                      IF( NOUNIT )   X( J ) = X( J )/A( J, J )
(  384)                      TEMP = X( J )
(  385)                      DO 50, I = J + 1, N
(  386)                         X( I ) = X( I ) - TEMP*A( I, J )
(  387)    50                CONTINUE
(  388)                   END IF
(  389)    60          CONTINUE
(  390)             ELSE
(  391)                JX = KX
(  392)                DO 80, J = 1, N
(  393)                   IF( X( JX ) /= ZERO )THEN
(  394)                      IF( NOUNIT )   X( JX ) = X( JX )/A( J, J )
(  395)                      TEMP = X( JX )
(  396)                      IX   = JX
(  397)                      DO 70, I = J + 1, N






PGF90 (Version     10.5)          09/23/2010  18:36:34      page 8

(  398)                         IX      = IX      + INCX
(  399)                         X( IX ) = X( IX ) - TEMP*A( I, J )
(  400)    70                CONTINUE
(  401)                   END IF
(  402)                   JX = JX + INCX
(  403)    80          CONTINUE
(  404)             END IF
(  405)          END IF
(  406)       ELSE
(  407) 
(  408) !!       Form  x := inv( A' )*x.
(  409) 
(  410)          IF( LSAME( UPLO, 'U' ) )THEN
(  411)             IF( INCX == 1 )THEN
(  412)                DO 100, J = 1, N
(  413)                   TEMP = X( J )
(  414)                   DO 90, I = 1, J - 1
(  415)                      TEMP = TEMP - A( I, J )*X( I )
(  416)    90             CONTINUE
(  417)                   IF( NOUNIT )   TEMP = TEMP/A( J, J )
(  418)                   X( J ) = TEMP
(  419)   100          CONTINUE
(  420)             ELSE
(  421)                JX = KX
(  422)                DO 120, J = 1, N
(  423)                   TEMP = X( JX )
(  424)                   IX   = KX
(  425)                   DO 110, I = 1, J - 1
(  426)                      TEMP = TEMP - A( I, J )*X( IX )
(  427)                      IX   = IX   + INCX
(  428)   110             CONTINUE
(  429)                   IF( NOUNIT )   TEMP = TEMP/A( J, J )
(  430)                   X( JX ) = TEMP
(  431)                   JX      = JX   + INCX
(  432)   120          CONTINUE
(  433)             END IF
(  434)          ELSE
(  435)             IF( INCX == 1 )THEN
(  436)                DO 140, J = N, 1, -1
(  437)                   TEMP = X( J )
(  438)                   DO 130, I = N, J + 1, -1
(  439)                      TEMP = TEMP - A( I, J )*X( I )
(  440)   130             CONTINUE
(  441)                   IF( NOUNIT )    TEMP = TEMP/A( J, J )
(  442)                   X( J ) = TEMP
(  443)   140          CONTINUE
(  444)             ELSE
(  445)                KX = KX + ( N - 1 )*INCX
(  446)                JX = KX
(  447)                DO 160, J = N, 1, -1
(  448)                   TEMP = X( JX )
(  449)                   IX   = KX
(  450)                   DO 150, I = N, J + 1, -1
(  451)                      TEMP = TEMP - A( I, J )*X( IX )
(  452)                      IX   = IX   - INCX
(  453)   150             CONTINUE
(  454)                   IF( NOUNIT )    TEMP = TEMP/A( J, J )
(  455)                   X( JX ) = TEMP






PGF90 (Version     10.5)          09/23/2010  18:36:34      page 9

(  456)                   JX      = JX   - INCX
(  457)   160          CONTINUE
(  458)             END IF
(  459)          END IF
(  460)       END IF
(  461) 
(  462)    END SUBROUTINE DTRSV
(  463) !---------------------------------------------------------
(  464)    SUBROUTINE DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
(  465) !!    .. Scalar Arguments ..
(  466) !     DOUBLE PRECISION   ALPHA
(  467)       REAL (kind=8) ::                 ALPHA
(  468)       INTEGER            INCX, INCY, LDA, M, N
(  469) !!    .. Array Arguments ..
(  470) !     DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
(  471)       REAL (kind=8) ::   A( LDA, * ), X( * ), Y( * )
(  472) !!    ..
(  473) 
(  474) !! Purpose
(  475) !! =======
(  476) 
(  477) !! DGER   performs the rank 1 operation
(  478) 
(  479) !!    A := alpha*x*y' + A,
(  480) 
(  481) !! where alpha is a scalar, x is an m element vector, y is an n element
(  482) !! vector and A is an m by n matrix.
(  483) 
(  484) !! Parameters
(  485) !! ==========
(  486) 
(  487) !! M      - INTEGER.
(  488) !!          On entry, M specifies the number of rows of the matrix A.
(  489) !!          M must be at least zero.
(  490) !!          Unchanged on exit.
(  491) 
(  492) !! N      - INTEGER.
(  493) !!          On entry, N specifies the number of columns of the matrix A.
(  494) !!          N must be at least zero.
(  495) !!          Unchanged on exit.
(  496) 
(  497) !! ALPHA  - DOUBLE PRECISION.
(  498) !!          On entry, ALPHA specifies the scalar alpha.
(  499) !!          Unchanged on exit.
(  500) 
(  501) !! X      - DOUBLE PRECISION array of dimension at least
(  502) !!          ( 1 + ( m - 1 )*abs( INCX ) ).
(  503) !!          Before entry, the incremented array X must contain the m
(  504) !!          element vector x.
(  505) !!          Unchanged on exit.
(  506) 
(  507) !! INCX   - INTEGER.
(  508) !!          On entry, INCX specifies the increment for the elements of
(  509) !!          X. INCX must not be zero.
(  510) !!          Unchanged on exit.
(  511) 
(  512) !! Y      - DOUBLE PRECISION array of dimension at least
(  513) !!          ( 1 + ( n - 1 )*abs( INCY ) ).






PGF90 (Version     10.5)          09/23/2010  18:36:34      page 10

(  514) !!          Before entry, the incremented array Y must contain the n
(  515) !!          element vector y.
(  516) !!          Unchanged on exit.
(  517) 
(  518) !! INCY   - INTEGER.
(  519) !!          On entry, INCY specifies the increment for the elements of
(  520) !!          Y. INCY must not be zero.
(  521) !!          Unchanged on exit.
(  522) 
(  523) !! A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
(  524) !!          Before entry, the leading m by n part of the array A must
(  525) !!          contain the matrix of coefficients. On exit, A is
(  526) !!          overwritten by the updated matrix.
(  527) 
(  528) !! LDA    - INTEGER.
(  529) !!          On entry, LDA specifies the first dimension of A as declared
(  530) !!          in the calling (sub) program. LDA must be at least
(  531) !!          max( 1, m ).
(  532) !!          Unchanged on exit.
(  533) 
(  534) 
(  535) !! Level 2 Blas routine.
(  536) 
(  537) !! -- Written on 22-October-1986.
(  538) !!    Jack Dongarra, Argonne National Lab.
(  539) !!    Jeremy Du Croz, Nag Central Office.
(  540) !!    Sven Hammarling, Nag Central Office.
(  541) !!    Richard Hanson, Sandia National Labs.
(  542) 
(  543) 
(  544) !!    .. Parameters ..
(  545) !     DOUBLE PRECISION   ZERO
(  546) !     PARAMETER        ( ZERO = 0.0D+0 )
(  547)       REAL (kind=8) ::                 ZERO
(  548)       PARAMETER        ( ZERO = 0.0 )
(  549) !!    .. Local Scalars ..
(  550) !     DOUBLE PRECISION   TEMP
(  551)       REAL (kind=8) ::                 TEMP
(  552)       INTEGER            I, INFO, IX, J, JY, KX
(  553) !!    .. External Subroutines ..
(  554)       EXTERNAL           XERBLA
(  555) !!    .. Intrinsic Functions ..
(  556)       INTRINSIC          MAX
(  557) !!    ..
(  558) !!    .. Executable Statements ..
(  559) 
(  560) !!    Test the input parameters.
(  561) 
(  562)       INFO = 0
(  563)       IF     ( M < 0 )THEN
(  564)          INFO = 1
(  565)       ELSE IF( N < 0 )THEN
(  566)          INFO = 2
(  567)       ELSE IF( INCX == 0 )THEN
(  568)          INFO = 5
(  569)       ELSE IF( INCY == 0 )THEN
(  570)          INFO = 7
(  571)       ELSE IF( LDA < MAX( 1, M ) )THEN






PGF90 (Version     10.5)          09/23/2010  18:36:34      page 11

(  572)          INFO = 9
(  573)       END IF
(  574)       IF( INFO /= 0 )THEN
(  575)          CALL XERBLA( 'DGER  ', INFO )
(  576)          RETURN
(  577)       END IF
(  578) 
(  579) !!    Quick return if possible.
(  580) 
(  581)       IF( ( M == 0 ).OR.( N == 0 ).OR.( ALPHA == ZERO ) )   &
(  582)          RETURN
(  583) 
(  584) !!    Start the operations. In this version the elements of A are
(  585) !!    accessed sequentially with one pass through A.
(  586) 
(  587)       IF( INCY > 0 )THEN
(  588)          JY = 1
(  589)       ELSE
(  590)          JY = 1 - ( N - 1 )*INCY
(  591)       END IF
(  592)       IF( INCX == 1 )THEN
(  593)          DO 20, J = 1, N
(  594)             IF( Y( JY ) /= ZERO )THEN
(  595)                TEMP = ALPHA*Y( JY )
(  596)                DO 10, I = 1, M
(  597)                   A( I, J ) = A( I, J ) + X( I )*TEMP
(  598)    10          CONTINUE
(  599)             END IF
(  600)             JY = JY + INCY
(  601)    20    CONTINUE
(  602)       ELSE
(  603)          IF( INCX > 0 )THEN
(  604)             KX = 1
(  605)          ELSE
(  606)             KX = 1 - ( M - 1 )*INCX
(  607)          END IF
(  608)          DO 40, J = 1, N
(  609)             IF( Y( JY ) /= ZERO )THEN
(  610)                TEMP = ALPHA*Y( JY )
(  611)                IX   = KX
(  612)                DO 30, I = 1, M
(  613)                   A( I, J ) = A( I, J ) + X( IX )*TEMP
(  614)                   IX        = IX        + INCX
(  615)    30          CONTINUE
(  616)             END IF
(  617)             JY = JY + INCY
(  618)    40    CONTINUE
(  619)       END IF
(  620) 
(  621)    END SUBROUTINE DGER
(  622) !---------------------------------------------------------
(  623)    SUBROUTINE XERBLA ( SRNAME, INFO )
(  624) !!    ..    Scalar Arguments ..
(  625)       INTEGER            INFO
(  626)       CHARACTER (len=6) ::    SRNAME
(  627) !!    ..
(  628) 
(  629) !! Purpose






PGF90 (Version     10.5)          09/23/2010  18:36:34      page 12

(  630) !! =======
(  631) 
(  632) !! XERBLA  is an error handler for the Level 2 BLAS routines.
(  633) 
(  634) !! It is called by the Level 2 BLAS routines if an input parameter is
(  635) !! invalid.
(  636) 
(  637) !! Installers should consider modifying the STOP statement in order to
(  638) !! call system-specific exception-handling facilities.
(  639) 
(  640) !! Parameters
(  641) !! ==========
(  642) 
(  643) !! SRNAME - CHARACTER*6.
(  644) !!          On entry, SRNAME specifies the name of the routine which
(  645) !!          called XERBLA.
(  646) 
(  647) !! INFO   - INTEGER.
(  648) !!          On entry, INFO specifies the position of the invalid
(  649) !!          parameter in the parameter-list of the calling routine.
(  650) 
(  651) 
(  652) !! Auxiliary routine for Level 2 Blas.
(  653) 
(  654) !! Written on 20-July-1986.
(  655) 
(  656) !!    .. Executable Statements ..
(  657) 
(  658)       WRITE (*,99999) SRNAME, INFO
(  659) 
(  660)       STOP
(  661) 
(  662) 99999 FORMAT ( ' ** On entry to ', A6, ' parameter number ', I2,   &
(  663)                ' had an illegal value' )
(  664) 
(  665)    END SUBROUTINE XERBLA
(  666) !-----------------------------------------------------------
(  667)    FUNCTION lsame( c1, c2 )
(  668)       logical lsame
(  669)       CHARACTER (len=*), INTENT(in) ::   c1, c2
(  670)       IF( c1 == c2 ) THEN
(  671)           lsame=.TRUE.
(  672)       ELSE
(  673)           lsame=.FALSE.
(  674)       ENDIF
(  675)    END FUNCTION lsame
