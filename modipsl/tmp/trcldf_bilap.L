


PGF90 (Version     10.5)          09/23/2010  09:49:20      page 1

Switches: -noasm -nodclchk -nodebug -nodlines -noline -list
          -idir /usr/local/include
          -idir ../../../lib
          -idir ../../../lib/oce
          -idir /usr/local/include
          -idir /usr/include/mpich2-x86_64
          -idir /usr/include/mpich2-x86_64
          -inform severe -opt 3 -nosave -object -noonetrip
          -depchk on -nostandard     
          -nosymbol -noupcase    

Filename: trcldf_bilap.F90

(    1) MODULE trcldf_bilap
(    2)    !!==============================================================================
(    3)    !!                   ***  MODULE  trcldf_bilap  ***
(    4)    !! Ocean passive tracers:  horizontal component of the lateral tracer mixing trend
(    5)    !!==============================================================================
(    6) #if defined key_passivetrc
(    7)    !!----------------------------------------------------------------------
(    8)    !!   trc_ldf_bilap : update the tracer trend with the horizontal diffusion
(    9)    !!                   using a iso-level biharmonic operator
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce_trc         ! ocean dynamics and active tracers variables
(   13)    USE trc             ! ocean passive tracers variables
(   14)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   15)    USE prtctl_trc      ! Print control for debbuging
(   16) 
(   17)    IMPLICIT NONE
(   18)    PRIVATE
(   19) 
(   20)    !! * Routine accessibility
(   21)    PUBLIC trc_ldf_bilap   ! routine called by step.F90
(   22) 
(   23)    !! * Substitutions
(   24) #  include "passivetrc_substitute.h90"
(   25)    !!----------------------------------------------------------------------
(   26)    !!   TOP 1.0 , LOCEAN-IPSL (2005) 
(   27)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trcldf_bilap.F90,v 1.9 2006/04/10 15:38:54 opalod Exp $ 
(   28)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   29)    !!----------------------------------------------------------------------
(   30) 
(   31) CONTAINS
(   32)    
(   33)    SUBROUTINE trc_ldf_bilap( kt )
(   34)       !!----------------------------------------------------------------------
(   35)       !!                  ***  ROUTINE trc_ldf_bilap  ***
(   36)       !!
(   37)       !! ** Purpose :   Compute the before horizontal tracer tra diffusive 
(   38)       !!      trend and add it to the general trend of tracer equation.
(   39)       !!
(   40)       !! ** Method  :   4th order diffusive operator along model level surfaces 
(   41)       !!      evaluated using before fields (forward time scheme). The hor.
(   42)       !!      diffusive trends of passive tracer is given by:
(   43)       !!       * s-coordinate ('key_s_coord' defined), the vertical scale 
(   44)       !!      factors e3. are inside the derivatives:
(   45)       !!      Laplacian of trb:
(   46)       !!         zlt   = 1/(e1t*e2t*e3t) {  di-1[ e2u*e3u/e1u di(trb) ]






PGF90 (Version     10.5)          09/23/2010  09:49:20      page 2

(   47)       !!                                  + dj-1[ e1v*e3v/e2v dj(trb) ]  }
(   48)       !!      Multiply by the eddy diffusivity coef. and insure lateral bc:
(   49)       !!        zlt   = ahtt * zlt
(   50)       !!        call to lbc_lnk
(   51)       !!      Bilaplacian (laplacian of zlt):
(   52)       !!         difft = 1/(e1t*e2t*e3t) {  di-1[ e2u*e3u/e1u di(zlt) ]
(   53)       !!                                  + dj-1[ e1v*e3v/e2v dj(zlt) ]  }
(   54)       !!       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
(   55)       !!      Laplacian of trb:
(   56)       !!         zlt   = 1/(e1t*e2t) {  di-1[ e2u/e1u di(trb) ]
(   57)       !!                              + dj-1[ e1v/e2v dj(trb) ] }
(   58)       !!      Multiply by the eddy diffusivity coef. and insure lateral bc:
(   59)       !!        zlt   = ahtt * zlt
(   60)       !!        call to lbc_lnk
(   61)       !!      Bilaplacian (laplacian of zlt):
(   62)       !!         difft = 1/(e1t*e2t) {  di-1[ e2u/e1u di(zlt) ]
(   63)       !!                              + dj-1[ e1v/e2v dj(zlt) ]  }
(   64)       !!
(   65)       !!      Add this trend to the general trend tra :
(   66)       !!         tra = tra + difft 
(   67)       !!
(   68)       !! ** Action : - Update tra arrays with the before iso-level
(   69)       !!               biharmonic mixing trend.
(   70)       !!             - Save the trends in trtrd ('key_trc_diatrd')
(   71)       !!
(   72)       !! History :
(   73)       !!        !  91-11  (G. Madec)  Original code
(   74)       !!        !  93-03  (M. Guyon)  symetrical conditions
(   75)       !!        !  95-11  (G. Madec)  suppress volumetric scale factors
(   76)       !!        !  96-01  (G. Madec)  statement function for e3
(   77)       !!        !  96-01  (M. Imbard)  mpp exchange
(   78)       !!        !  97-07  (G. Madec)  optimization, and ahtt
(   79)       !!        !  00-05  (MA Foujols) add lbc for tracer trends
(   80)       !!        !  00-10  (MA Foujols E. Kestenare) use passive tracer coefficient
(   81)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   82)       !!   9.0  !  04-03  (C. Ethe )  F90: Free form and module
(   83)       !!----------------------------------------------------------------------
(   84)       !! * Arguments
(   85)       INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
(   86) 
(   87)       !! * Local declarations
(   88)       INTEGER ::   ji, jj, jk, jn             ! dummy loop indices
(   89) #if defined key_partial_steps
(   90)       INTEGER ::   iku, ikv                   ! temporary integers
(   91) #endif
(   92)       REAL(wp) ::   ztra     ! temporary scalars
(   93) 
(   94)       REAL(wp), DIMENSION(jpi,jpj) ::   & 
(   95)          zeeu, zeev, zbtr, zlt                 ! workspace
(   96)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   & 
(   97)          ztu, ztv                              ! workspace
(   98)       CHARACTER (len=22) :: charout
(   99)       !!----------------------------------------------------------------------
(  100) 
(  101)       IF( kt == nittrc000 ) THEN
(  102)          IF(lwp) WRITE(numout,*)
(  103)          IF(lwp) WRITE(numout,*) 'trc_ldf_bilap : iso-level biharmonic operator'
(  104)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'






PGF90 (Version     10.5)          09/23/2010  09:49:20      page 3

(  105)       ENDIF
(  106)       ! 
(  107) 
(  108)       DO jn = 1, jptra
(  109)                                                           ! ===============
(  110)          DO jk = 1, jpkm1                                 ! Horizontal slab
(  111)             !                                             ! ===============
(  112) 
(  113)             ! 0. Initialization of metric arrays (for z- or s-coordinates)
(  114)             ! ----------------------------------
(  115) 
(  116)             DO jj = 1, jpjm1
(  117)                DO ji = 1, fs_jpim1   ! vector opt.
(  118) #if defined key_s_coord || defined key_partial_steps
(  119)                   ! s-coordinates, vertical scale factor are used
(  120)                   zbtr(ji,jj) = 1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
(  121)                   zeeu(ji,jj) = e2u(ji,jj) * fse3u(ji,jj,jk) / e1u(ji,jj) * umask(ji,jj,jk)
(  122)                   zeev(ji,jj) = e1v(ji,jj) * fse3v(ji,jj,jk) / e2v(ji,jj) * vmask(ji,jj,jk)
(  123) #else
(  124)                   ! z-coordinates, no vertical scale factors
(  125)                   zbtr(ji,jj) = 1. / ( e1t(ji,jj)*e2t(ji,jj) )
(  126)                   zeeu(ji,jj) = e2u(ji,jj) / e1u(ji,jj) * umask(ji,jj,jk)
(  127)                   zeev(ji,jj) = e1v(ji,jj) / e2v(ji,jj) * vmask(ji,jj,jk)
(  128) #endif
(  129)                END DO
(  130)             END DO
(  131) 
(  132) 
(  133)             ! 1. Laplacian
(  134)             ! ------------
(  135) 
(  136)             ! First derivative (gradient)
(  137)             DO jj = 1, jpjm1
(  138)                DO ji = 1, fs_jpim1   ! vector opt.
(  139)                   ztu(ji,jj,jk) = zeeu(ji,jj) * ( trb(ji+1,jj  ,jk,jn) - trb(ji,jj,jk,jn) )
(  140)                   ztv(ji,jj,jk) = zeev(ji,jj) * ( trb(ji  ,jj+1,jk,jn) - trb(ji,jj,jk,jn) )
(  141)                END DO
(  142)             END DO
(  143) #if defined key_partial_steps
(  144) 
(  145)             DO jj = 1, jpj-1
(  146)                DO ji = 1, jpi-1
(  147)                   ! last level
(  148)                   iku = MIN ( mbathy(ji,jj), mbathy(ji+1,jj  ) ) - 1
(  149)                   ikv = MIN ( mbathy(ji,jj), mbathy(ji  ,jj+1) ) - 1
(  150)                   IF( iku == jk ) THEN
(  151)                      ztu(ji,jj,jk) = zeeu(ji,jj) * gtu(ji,jj)
(  152)                   ENDIF
(  153)                   IF( ikv == jk ) THEN
(  154)                      ztv(ji,jj,jk) = zeev(ji,jj) * gtv(ji,jj)
(  155)                   ENDIF
(  156)                END DO
(  157)             END DO
(  158) #endif
(  159) 
(  160)             ! Second derivative (divergence)
(  161)             DO jj = 2, jpjm1
(  162)                DO ji = fs_2, fs_jpim1   ! vector opt.






PGF90 (Version     10.5)          09/23/2010  09:49:20      page 4

(  163)                   zlt(ji,jj) = zbtr(ji,jj) * (  ztu(ji,jj,jk) - ztu(ji-1,jj,jk) + ztv(ji,jj,jk) - ztv(ji,jj-1,jk)  )
(  164)                END DO
(  165)             END DO
(  166) 
(  167)             ! Multiply by the eddy diffusivity coefficient
(  168)             DO jj = 2, jpjm1
(  169)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  170)                   zlt(ji,jj) = fsahtrt(ji,jj,jk) * zlt(ji,jj)
(  171)                END DO
(  172)             END DO
(  173) 
(  174)             ! Lateral boundary conditions on the laplacian zlt   (unchanged sgn)
(  175)             CALL lbc_lnk( zlt, 'T', 1. )  
(  176) 
(  177)             ! 2. Bilaplacian
(  178)             ! --------------
(  179) 
(  180)             ! third derivative (gradient)
(  181)             DO jj = 1, jpjm1
(  182)                DO ji = 1, fs_jpim1   ! vector opt.
(  183)                   ztu(ji,jj,jk) = zeeu(ji,jj) * ( zlt(ji+1,jj  ) - zlt(ji,jj) )
(  184)                   ztv(ji,jj,jk) = zeev(ji,jj) * ( zlt(ji  ,jj+1) - zlt(ji,jj) )
(  185)                END DO
(  186)             END DO
(  187) 
(  188)             ! fourth derivative (divergence) and add to the general tracer trend
(  189)             DO jj = 2, jpjm1
(  190)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  191)                   ! horizontal diffusive trends
(  192)                   ztra = zbtr(ji,jj) * (  ztu(ji,jj,jk) - ztu(ji-1,jj,jk) + ztv(ji,jj,jk) - ztv(ji,jj-1,jk)  )
(  193)                   ! add it to the general tracer trends
(  194)                   tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + ztra
(  195) #if defined key_trc_diatrd
(  196)                   ! save the horizontal diffusive trends
(  197)                   IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),4) = (  ztu(ji,jj,jk) - ztu(ji-1,jj,jk) ) * zbtr(ji,jj)
(  198)                   IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),5) = (  ztv(ji,jj,jk) - ztv(ji-1,jj,jk) ) * zbtr(ji,jj)
(  199) #endif
(  200)                END DO
(  201)             END DO
(  202)             !                                             ! ===============
(  203)          END DO                                           ! Horizontal slab
(  204)          !                                                ! ===============
(  205) #if defined key_trc_diatrd
(  206)          ! Lateral boundary conditions on the laplacian zlt   (unchanged sgn)
(  207)          IF (luttrd(jn)) CALL lbc_lnk( trtrd(:,:,:,ikeep(jn),5), 'T', 1. )  
(  208) #endif
(  209)       END DO
(  210) 
(  211)      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
(  212)          WRITE(charout, FMT="('ldf - bilap')")
(  213)          CALL prt_ctl_trc_info(charout)
(  214)          CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
(  215)       ENDIF
(  216) 
(  217)    END SUBROUTINE trc_ldf_bilap
(  218) 
(  219) #else
(  220)    !!----------------------------------------------------------------------






PGF90 (Version     10.5)          09/23/2010  09:49:20      page 5

(  221)    !!   Default option                                         Empty module
(  222)    !!----------------------------------------------------------------------
(  223) CONTAINS
(  224)    SUBROUTINE trc_ldf_bilap( kt )  
(  225)       INTEGER, INTENT(in) :: kt
(  226) !      WRITE(*,*) 'trc_ldf_bilap: You should not have seen this print! error?', kt
(  227)    END SUBROUTINE trc_ldf_bilap
(  228) #endif
(  229)    !!==============================================================================
(  230) END MODULE trcldf_bilap






















































PGF90 (Version     10.5)          09/23/2010  09:49:20      page 6

(    1) # 1 "trcldf_bilap.F90"
(    1) MODULE trcldf_bilap
(    2)    !!==============================================================================
(    3)    !!                   ***  MODULE  trcldf_bilap  ***
(    4)    !! Ocean passive tracers:  horizontal component of the lateral tracer mixing trend
(    5)    !!==============================================================================
(    6) # 220
(  220)    !!----------------------------------------------------------------------
(  221)    !!   Default option                                         Empty module
(  222)    !!----------------------------------------------------------------------
(  223) CONTAINS
(  224)    SUBROUTINE trc_ldf_bilap( kt )  
(  225)       INTEGER, INTENT(in) :: kt
(  226) !      WRITE(*,*) 'trc_ldf_bilap: You should not have seen this print! error?', kt
(  227)    END SUBROUTINE trc_ldf_bilap
(  228) # 229
(  229)    !!==============================================================================
(  230) END MODULE trcldf_bilap
