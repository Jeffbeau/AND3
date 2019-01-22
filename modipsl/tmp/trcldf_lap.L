


PGF90 (Version     10.5)          09/23/2010  09:49:21      page 1

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

Filename: trcldf_lap.F90

(    1) MODULE trcldf_lap
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  trcldf_lap  ***
(    4)    !! Ocean passive tracers:  horizontal component of the lateral tracer mixing trend
(    5)    !!==============================================================================
(    6) #if defined key_passivetrc
(    7)    !!----------------------------------------------------------------------
(    8)    !!   trc_ldf_lap  : update the tracer trend with the horizontal diffusion
(    9)    !!                 using a iso-level harmonic (laplacien) operator.
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce_trc             ! ocean dynamics and active tracers variables
(   13)    USE trc                 ! ocean passive tracers variables
(   14)    USE prtctl_trc          ! Print control for debbuging
(   15) 
(   16)    IMPLICIT NONE
(   17)    PRIVATE
(   18) 
(   19)    !! * Routine accessibility
(   20)    PUBLIC trc_ldf_lap  ! routine called by step.F90
(   21) 
(   22)    !! * Substitutions
(   23) #  include "passivetrc_substitute.h90"
(   24)    !!----------------------------------------------------------------------
(   25)    !!   TOP 1.0 , LOCEAN-IPSL (2005) 
(   26)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trcldf_lap.F90,v 1.9 2006/04/10 15:38:55 opalod Exp $ 
(   27)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   28)    !!----------------------------------------------------------------------
(   29)    
(   30) CONTAINS
(   31) 
(   32)    SUBROUTINE trc_ldf_lap( kt )
(   33)       !!----------------------------------------------------------------------
(   34)       !!                  ***  ROUTINE trc_ldf_lap  ***
(   35)       !!                   
(   36)       !! ** Purpose :   Compute the before horizontal tracer diffusive 
(   37)       !!      trend and add it to the general trend of tracer equation.
(   38)       !!
(   39)       !! ** Method  :   Second order diffusive operator evaluated using before
(   40)       !!      fields (forward time scheme). The horizontal diffusive trends of 
(   41)       !!      the passive tracer is given by:
(   42)       !!       * s-coordinate ('key_s_coord' defined), the vertical scale 
(   43)       !!      factors e3. are inside the derivatives:
(   44)       !!          difft = 1/(e1t*e2t*e3t) {  di-1[ aht e2u*e3u/e1u di(trb) ]
(   45)       !!                                   + dj-1[ aht e1v*e3v/e2v dj(trb) ] }
(   46)       !!       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:






PGF90 (Version     10.5)          09/23/2010  09:49:21      page 2

(   47)       !!          difft = 1/(e1t*e2t) {  di-1[ aht e2u/e1u di(trb) ]
(   48)       !!                               + dj-1[ aht e1v/e2v dj(trb) ] }
(   49)       !!      Add this trend to the general tracer trend tra :
(   50)       !!          tra = tra + difft
(   51)       !!
(   52)       !! ** Action  : - Update tra arrays with the before iso-level 
(   53)       !!                harmonic mixing trend.
(   54)       !!              - Save the trends in trtrd ('key_trc_diatrd')
(   55)       !!
(   56)       !! History :
(   57)       !!   1.0  !  87-06  (P. Andrich, D. L Hostis)  Original code
(   58)       !!        !  91-11  (G. Madec)
(   59)       !!        !  95-02  (M. Levy)    passive tracers
(   60)       !!        !  95-11  (G. Madec)  suppress volumetric scale factors
(   61)       !!        !  96-01  (G. Madec)  statement function for e3
(   62)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   63)       !!   9.0  !  04-03  (C. Ethe)   passive tracer
(   64)       !!----------------------------------------------------------------------
(   65)       USE oce_trc          , ztu => ua,  &  ! use ua as workspace
(   66)          &                   ztv => va      ! use va as workspace
(   67) 
(   68)       !! * Arguments
(   69)       INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
(   70)       
(   71)       !! * Local save
(   72)       REAL(wp), DIMENSION(jpi,jpj), SAVE ::   &
(   73)          ze1ur, ze2vr, zbtr2              ! scale factor coefficients
(   74)       
(   75)       !! * Local declarations
(   76)       INTEGER ::   ji, jj, jk,jn         ! dummy loop indices
(   77)       REAL(wp) ::   &
(   78)          zabe1, zabe2, zbtr              ! temporary scalars
(   79) 
(   80)       REAL(wp) ::   &
(   81)          ztra, ztrax, ztray              ! workspace
(   82)       CHARACTER (len=22) :: charout
(   83)       !!----------------------------------------------------------------------
(   84)       
(   85)       IF( kt == nittrc000 ) THEN
(   86)          IF(lwp) WRITE(numout,*)
(   87)          IF(lwp) WRITE(numout,*) 'trc_ldf_lap : iso-level laplacian diffusion'
(   88)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
(   89)          ze1ur(:,:) = e2u(:,:) / e1u(:,:)
(   90)          ze2vr(:,:) = e1v(:,:) / e2v(:,:)
(   91)          zbtr2(:,:) = 1. / ( e1t(:,:) * e2t(:,:) )
(   92)       ENDIF
(   93)  
(   94)       DO jn = 1, jptra
(   95)          
(   96)          !                                                  ! =============
(   97)          DO jk = 1, jpkm1                                   ! Vertical slab
(   98)             !                                               ! =============
(   99)             ! 1. First derivative (gradient)
(  100)             ! -------------------
(  101)             DO jj = 1, jpjm1
(  102)                DO ji = 1, fs_jpim1   ! vector opt.
(  103) #if defined key_s_coord
(  104)                   zabe1 = fsahtru(ji,jj,jk) * umask(ji,jj,jk) * ze1ur(ji,jj) * fse3u(ji,jj,jk)






PGF90 (Version     10.5)          09/23/2010  09:49:21      page 3

(  105)                   zabe2 = fsahtrv(ji,jj,jk) * vmask(ji,jj,jk) * ze2vr(ji,jj) * fse3v(ji,jj,jk)
(  106) #else
(  107)                   zabe1 = fsahtru(ji,jj,jk) * umask(ji,jj,jk) * ze1ur(ji,jj)
(  108)                   zabe2 = fsahtrv(ji,jj,jk) * vmask(ji,jj,jk) * ze2vr(ji,jj)
(  109) #endif
(  110)                   ztu(ji,jj,jk) = zabe1 * ( trb(ji+1,jj  ,jk,jn) - trb(ji,jj,jk,jn) )
(  111)                   ztv(ji,jj,jk) = zabe2 * ( trb(ji  ,jj+1,jk,jn) - trb(ji,jj,jk,jn) )
(  112)                END DO
(  113)             END DO
(  114) 
(  115) 
(  116)             ! 2. Second derivative (divergence)
(  117)             ! --------------------
(  118)             DO jj = 2, jpjm1
(  119)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  120) #if defined key_s_coord
(  121)                   zbtr = zbtr2(ji,jj) / fse3t(ji,jj,jk)
(  122) #else
(  123)                   zbtr = zbtr2(ji,jj)
(  124) #endif
(  125)                   ! horizontal diffusive trends
(  126)                   ztrax = zbtr * ( ztu(ji,jj,jk) - ztu(ji-1,jj,jk) )
(  127)                   ztray = zbtr * ( ztv(ji,jj,jk) - ztv(ji,jj-1,jk) )
(  128) 
(  129)                   ! add it to the general tracer trends
(  130)                   tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + ztrax + ztray
(  131) 
(  132) #if defined key_trc_diatrd
(  133)                   ! save the horizontal diffusive trends
(  134)                   IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),4) = ztrax
(  135)                   IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),5) = ztray
(  136) #endif
(  137)                END DO
(  138)             END DO
(  139)             !                                               ! =============
(  140)          END DO                                             !  End of slab  
(  141)          !                                                  ! =============
(  142) 
(  143)       END DO
(  144) 
(  145)      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
(  146)          WRITE(charout, FMT="('ldf - lap')")
(  147)          CALL prt_ctl_trc_info(charout)
(  148)          CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
(  149)       ENDIF
(  150) 
(  151)    END SUBROUTINE trc_ldf_lap
(  152) 
(  153) #else
(  154)    !!----------------------------------------------------------------------
(  155)    !!   Default option                                         Empty module
(  156)    !!----------------------------------------------------------------------
(  157) CONTAINS
(  158)    SUBROUTINE trc_ldf_lap( kt )  
(  159)       INTEGER, INTENT(in) :: kt
(  160) !      WRITE(*,*) 'trc_ldf_lap: You should not have seen this print! error?', kt
(  161)    END SUBROUTINE trc_ldf_lap
(  162) #endif






PGF90 (Version     10.5)          09/23/2010  09:49:21      page 4

(  163) 
(  164)    !!==============================================================================
(  165) END MODULE trcldf_lap





























































PGF90 (Version     10.5)          09/23/2010  09:49:21      page 5

(    1) # 1 "trcldf_lap.F90"
(    1) MODULE trcldf_lap
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  trcldf_lap  ***
(    4)    !! Ocean passive tracers:  horizontal component of the lateral tracer mixing trend
(    5)    !!==============================================================================
(    6) # 154
(  154)    !!----------------------------------------------------------------------
(  155)    !!   Default option                                         Empty module
(  156)    !!----------------------------------------------------------------------
(  157) CONTAINS
(  158)    SUBROUTINE trc_ldf_lap( kt )  
(  159)       INTEGER, INTENT(in) :: kt
(  160) !      WRITE(*,*) 'trc_ldf_lap: You should not have seen this print! error?', kt
(  161)    END SUBROUTINE trc_ldf_lap
(  162) 
(  163) # 164
(  164)    !!==============================================================================
(  165) END MODULE trcldf_lap
