


PGF90 (Version     10.5)          09/23/2010  09:49:22      page 1

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

Filename: trczdf_exp.F90

(    1) MODULE trczdf_exp
(    2)    !!==============================================================================
(    3)    !!                    ***  MODULE  trczdf_exp  ***
(    4)    !! Ocean passive tracers:  vertical component of the tracer mixing trend using
(    5)    !!                        an explicit time-stepping (time spllitting scheme)
(    6)    !!==============================================================================
(    7) #if defined key_passivetrc
(    8)    !!----------------------------------------------------------------------
(    9)    !!   trc_zdf_exp  : update the tracer trend with the vertical diffusion
(   10)    !!                  using an explicit time stepping
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE oce_trc          ! ocean dynamics and active tracers variables
(   14)    USE trc              ! ocean passive tracers variables
(   15)    USE trctrp_lec       ! passive tracers transport
(   16)    USE prtctl_trc          ! Print control for debbuging
(   17) 
(   18)    IMPLICIT NONE
(   19)    PRIVATE
(   20) 
(   21)    !! * Routine accessibility
(   22)    PUBLIC trc_zdf_exp          ! routine called by step.F90
(   23) 
(   24)    !! * Module variable
(   25)    REAL(wp), DIMENSION(jpk) ::   &
(   26)       rdttrc                     ! vertical profile of 2 x tracer time-step
(   27) 
(   28)    !! * Substitutions
(   29) #  include "passivetrc_substitute.h90"
(   30)    !!----------------------------------------------------------------------
(   31)    !!   TOP 1.0 , LOCEAN-IPSL (2005) 
(   32)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trczdf_exp.F90,v 1.8 2005/12/07 10:30:00 opalod Exp $ 
(   33)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   34)    !!----------------------------------------------------------------------
(   35) 
(   36) CONTAINS
(   37) 
(   38)    SUBROUTINE trc_zdf_exp( kt )
(   39)       !!----------------------------------------------------------------------
(   40)       !!                  ***  ROUTINE trc_zdf_exp  ***
(   41)       !!                   
(   42)       !! ** Purpose :   Compute the trend due to the vertical tracer mixing 
(   43)       !!      using an explicit time stepping and add it to the general trend 
(   44)       !!      of the tracer equations.
(   45)       !!
(   46)       !! ** Method  :   The vertical diffusion of tracers  is given by:






PGF90 (Version     10.5)          09/23/2010  09:49:22      page 2

(   47)       !!         difft = dz( avt dz(trb) ) = 1/e3t dk+1( avt/e3w dk(trb) )
(   48)       !!      It is evaluated with an Euler scheme, using a time splitting
(   49)       !!      technique.
(   50)       !!      Surface and bottom boundary conditions: no diffusive flux on
(   51)       !!      both tracers (bottom, applied through the masked field avt).
(   52)       !!      Add this trend to the general trend tra :
(   53)       !!          tra = tra + dz( avt dz(t) ) if lk_zdfddm= T)
(   54)       !!
(   55)       !! ** Action : - Update tra with the before vertical diffusion trend
(   56)       !!             - Save the trends  in trtrd ('key_trc_diatrd')
(   57)       !!
(   58)       !! History :
(   59)       !!   6.0  !  90-10  (B. Blanke)  Original code
(   60)       !!   7.0  !  91-11  (G. Madec)
(   61)       !!        !  92-06  (M. Imbard)  correction on tracer trend loops
(   62)       !!        !  96-01  (G. Madec)  statement function for e3
(   63)       !!        !  97-05  (G. Madec)  vertical component of isopycnal
(   64)       !!        !  97-07  (G. Madec)  geopotential diffusion in s-coord
(   65)       !!        !  98-03  (L. Bopp MA Foujols) passive tracer generalisation
(   66)       !!        !  00-05  (MA Foujols) add lbc for tracer trends
(   67)       !!        !  00-06  (O Aumont)  correct isopycnal scheme suppress
(   68)       !!        !                     avt multiple correction
(   69)       !!        !  00-08  (G. Madec)  double diffusive mixing
(   70)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   71)       !!   9.0  !  04-03  (C. Ethe )  adapted for passive tracers
(   72)       !!---------------------------------------------------------------------
(   73)       !! * Arguments
(   74)       INTEGER, INTENT( in ) ::   kt           ! ocean time-step index
(   75)       
(   76)       !! * Local declarations
(   77)       INTEGER ::   ji, jj, jk, jl, jn             ! dummy loop indices
(   78)       REAL(wp) ::   &
(   79)          zlavmr,                 &  ! ???
(   80)          zave3r, ze3tr,          &  ! ???
(   81)          ztra                  !
(   82)       REAL(wp), DIMENSION(jpi,jpk) ::   &
(   83)          zwx, zwy
(   84)       CHARACTER (len=22) :: charout
(   85)       !!---------------------------------------------------------------------
(   86) 
(   87)       IF(lwp .AND. kt == nittrc000 ) THEN
(   88)          WRITE(numout,*)
(   89)          WRITE(numout,*) 'trc_zdf_exp : vertical tracer mixing'
(   90)          WRITE(numout,*) '~~~~~~~~~~~~~~~'
(   91)       ENDIF
(   92) 
(   93)       ! 0. Local constant initialization
(   94)       ! --------------------------------
(   95)       IF( ln_trcadv_cen2 .OR. ln_trcadv_tvd ) THEN
(   96)          ! time step = 2 rdttra with Arakawa or TVD advection scheme
(   97)          IF( neuler == 0 .AND. kt == nittrc000 ) THEN
(   98)             rdttrc(:) =  rdttra(:) * FLOAT(ndttrc)             ! restarting with Euler time stepping
(   99)          ELSEIF( kt <= nittrc000 + 1 ) THEN
(  100)             rdttrc(:) = 2. * rdttra(:) * FLOAT(ndttrc)         ! leapfrog
(  101)          ENDIF
(  102)       ELSE
(  103)          rdttrc(:) =  rdttra(:) * FLOAT(ndttrc)      
(  104)       ENDIF






PGF90 (Version     10.5)          09/23/2010  09:49:22      page 3

(  105) 
(  106) 
(  107)       zlavmr = 1. / FLOAT( n_trczdf_exp )
(  108) 
(  109)       DO jn = 1, jptra
(  110) 
(  111)          !                                                ! ===============
(  112)          DO jj = 2, jpjm1                                 !  Vertical slab
(  113)             !                                             ! ===============
(  114)             ! 1. Initializations
(  115)             ! ------------------
(  116) 
(  117)             ! Surface & bottom boundary conditions: no flux
(  118)             DO ji = 2, jpim1
(  119)                zwy(ji, 1 ) = 0.e0
(  120)                zwy(ji,jpk) = 0.e0
(  121)             END DO
(  122) 
(  123)             ! zwx and zwz arrays set to before tracer values
(  124)             DO jk = 1, jpk
(  125)                DO ji = 2, jpim1
(  126)                   zwx(ji,jk) = trb(ji,jj,jk,jn)
(  127)                END DO
(  128)             END DO
(  129) 
(  130) 
(  131)             ! 2. Time splitting loop
(  132)             ! ----------------------
(  133) 
(  134)             DO jl = 1, n_trczdf_exp
(  135) 
(  136)                ! first vertical derivative
(  137)                ! double diffusion: fstravs(ji,jj,jk) = avt(ji,jj,jk) /= avs (key_trc_zdfddm) 
(  138)                !                   fstravs(ji,jj,jk) = avs(ji,jj,jk) = avt
(  139)                DO jk = 2, jpk
(  140)                   DO ji = 2, jpim1
(  141)                      zave3r = 1.e0 / fse3w(ji,jj,jk) 
(  142)                      zwy(ji,jk) = fstravs(ji,jj,jk) * ( zwx(ji,jk-1) - zwx(ji,jk) ) * zave3r
(  143)                   END DO
(  144)                END DO
(  145) 
(  146) 
(  147)                ! trend estimation at kt+l*2*rdt/n_zdfexp
(  148)                DO jk = 1, jpkm1
(  149)                   DO ji = 2, jpim1
(  150)                      ze3tr = zlavmr / fse3t(ji,jj,jk)
(  151)                      ! 2nd vertical derivative
(  152)                      ztra = ( zwy(ji,jk) - zwy(ji,jk+1) ) * ze3tr
(  153)                      ! update the tracer trends
(  154)                      tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + ztra
(  155)                      ! update tracer fields at kt+l*2*rdt/n_trczdf_exp
(  156)                      zwx(ji,jk) = zwx(ji,jk) + rdttrc(jk) * ztra * tmask(ji,jj,jk)
(  157)                   END DO
(  158)                END DO
(  159)             END DO
(  160)             !                                             ! ===============
(  161)          END DO                                           !   End of slab
(  162)          !                                                ! ===============






PGF90 (Version     10.5)          09/23/2010  09:49:22      page 4

(  163)       END DO
(  164) 
(  165)       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
(  166)          WRITE(charout, FMT="('zdf - exp')")
(  167)          CALL prt_ctl_trc_info(charout)
(  168)          CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
(  169)       ENDIF
(  170) 
(  171)    END SUBROUTINE trc_zdf_exp
(  172) 
(  173) #else
(  174)    !!----------------------------------------------------------------------
(  175)    !!   Dummy module :                      NO passive tracer
(  176)    !!----------------------------------------------------------------------
(  177) CONTAINS
(  178)    SUBROUTINE trc_zdf_exp (kt )              ! Empty routine
(  179)       INTEGER, INTENT(in) :: kt
(  180) !      WRITE(*,*) 'trc_zdf_exp: You should not have seen this print! error?', kt
(  181)    END SUBROUTINE trc_zdf_exp
(  182) #endif
(  183)    
(  184)    !!==============================================================================
(  185) END MODULE trczdf_exp









































PGF90 (Version     10.5)          09/23/2010  09:49:22      page 5

(    1) # 1 "trczdf_exp.F90"
(    1) MODULE trczdf_exp
(    2)    !!==============================================================================
(    3)    !!                    ***  MODULE  trczdf_exp  ***
(    4)    !! Ocean passive tracers:  vertical component of the tracer mixing trend using
(    5)    !!                        an explicit time-stepping (time spllitting scheme)
(    6)    !!==============================================================================
(    7) # 174
(  174)    !!----------------------------------------------------------------------
(  175)    !!   Dummy module :                      NO passive tracer
(  176)    !!----------------------------------------------------------------------
(  177) CONTAINS
(  178)    SUBROUTINE trc_zdf_exp (kt )              ! Empty routine
(  179)       INTEGER, INTENT(in) :: kt
(  180) !      WRITE(*,*) 'trc_zdf_exp: You should not have seen this print! error?', kt
(  181)    END SUBROUTINE trc_zdf_exp
(  182) # 183
(  183)    
(  184)    !!==============================================================================
(  185) END MODULE trczdf_exp
