


PGF90 (Version     10.5)          09/23/2010  09:49:23      page 1

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

Filename: trczdf_imp.F90

(    1) MODULE trczdf_imp
(    2)    !!==============================================================================
(    3)    !!                    ***  MODULE  trczdf_imp  ***
(    4)    !! Ocean passive tracers:  vertical component of the tracer mixing trend
(    5)    !!==============================================================================
(    6) #if defined key_passivetrc
(    7)    !!----------------------------------------------------------------------
(    8)    !!   trc_zdf_imp  : update the tracer trend with the vertical diffusion
(    9)    !!                  using an implicit time-stepping.
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce_trc             ! ocean dynamics and active tracers variables
(   13)    USE trc                 ! ocean passive tracers variables
(   14) !!DB
(   15) #ifdef key_BGCM_01
(   16) !   USE bgcm_01_initrc
(   17)    USE lib_bgcm_01
(   18) #else
(   19)    USE trctrp_lec      ! passive tracers transport
(   20) #endif
(   21)    USE prtctl_trc
(   22) 
(   23)    IMPLICIT NONE
(   24)    PRIVATE
(   25) 
(   26)    !! * Routine accessibility
(   27)    PUBLIC trc_zdf_imp          ! routine called by step.F90
(   28) 
(   29)    !! * Module variable
(   30)    REAL(wp), DIMENSION(jpk) ::   &
(   31)       rdttrc                     ! vertical profile of 2 x tracer time-step
(   32) 
(   33)    !! * Substitutions
(   34) #  include "passivetrc_substitute.h90"
(   35)    !!----------------------------------------------------------------------
(   36)    !!   TOP 1.0 , LOCEAN-IPSL (2005) 
(   37)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trczdf_imp.F90,v 1.9 2006/04/10 15:38:55 opalod Exp $ 
(   38)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   39)    !!----------------------------------------------------------------------
(   40) 
(   41) CONTAINS
(   42) 
(   43)    SUBROUTINE trc_zdf_imp( kt )
(   44)       !!----------------------------------------------------------------------
(   45)       !!                  ***  ROUTINE trc_zdf_imp  ***
(   46)       !!






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 2

(   47)       !! ** Purpose :   Compute the trend due to the vertical tracer mixing 
(   48)       !!      using an implicit time stepping and add it to the general trend
(   49)       !!      of the tracer equations.
(   50)       !!
(   51)       !! ** Method  :   The vertical diffusion of tracers tra is given by:
(   52)       !!          difft = dz( avt dz(t) ) = 1/e3t dk+1( avt/e3w dk(tra) )
(   53)       !!      It is thus evaluated using a backward time scheme
(   54)       !!      Surface and bottom boundary conditions: no diffusive flux on
(   55)       !!      both tracers (bottom, applied through the masked field avt).
(   56)       !!      Add this trend to the general trend tra :
(   57)       !!          tra = tra + dz( avt dz(t) )
(   58)       !!         (tra = tra + dz( avs dz(t) ) if lk_zdfddmtrc=T)
(   59)       !!
(   60)       !! ** Action  : - Update tra with the before vertical diffusion trend
(   61)       !!              - save the trends in trtrd ('key_trc_diatrd')
(   62)       !!
(   63)       !! History :
(   64)       !!   6.0  !  90-10  (B. Blanke)  Original code
(   65)       !!   7.0  !  91-11  (G. Madec)
(   66)       !!        !  92-06  (M. Imbard)  correction on tracer trend loops
(   67)       !!        !  96-01  (G. Madec)  statement function for e3
(   68)       !!        !  97-05  (G. Madec)  vertical component of isopycnal
(   69)       !!        !  97-07  (G. Madec)  geopotential diffusion in s-coord
(   70)       !!        !  98-03  (L. Bopp MA Foujols) passive tracer generalisation
(   71)       !!        !  00-05  (MA Foujols) add lbc for tracer trends
(   72)       !!        !  00-06  (O Aumont)  correct isopycnal scheme suppress
(   73)       !!        !                     avt multiple correction
(   74)       !!        !  00-08  (G. Madec)  double diffusive mixing
(   75)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   76)       !!   9.0  !  04-03  (C. Ethe )  adapted for passive tracers
(   77)       !!---------------------------------------------------------------------
(   78)       !! * Arguments
(   79)       INTEGER, INTENT( in ) ::   kt           ! ocean time-step index
(   80)       INTEGER ::   ikst, ikenm2, ikstp1
(   81)       !! * Local declarations
(   82)       INTEGER ::   ji, jj, jk, jn             ! dummy loop indices
(   83)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(   84)          zwd, zws, zwi,          &  ! ???
(   85)          zwx, zwy, zwt              ! ???
(   86)       REAL(wp) ::  ztra      ! temporary scalars
(   87) 
(   88)       REAL(wp), DIMENSION(jpi,jpj,jpk,jptra) ::   &
(   89)          ztrd
(   90)       CHARACTER (len=22) :: charout
(   91)       !!---------------------------------------------------------------------
(   92) 
(   93)       IF(lwp .AND. kt == nittrc000 ) THEN
(   94)          WRITE(numout,*)
(   95)          WRITE(numout,*) 'trc_zdf_implicit : vertical tracer mixing'
(   96)          WRITE(numout,*) '~~~~~~~~~~~~~~~'
(   97)       ENDIF
(   98) 
(   99)       ! 0. Local constant initialization
(  100)       ! --------------------------------
(  101)       IF( ln_trcadv_cen2 .OR. ln_trcadv_tvd ) THEN
(  102)          ! time step = 2 rdttra with Arakawa or TVD advection scheme
(  103)          IF( neuler == 0 .AND. kt == nittrc000 ) THEN
(  104)             rdttrc(:) =  rdttra(:) * FLOAT(ndttrc)             ! restarting with Euler time stepping






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 3

(  105)          ELSEIF( kt <= nittrc000 + 1 ) THEN
(  106)             rdttrc(:) = 2. * rdttra(:) * FLOAT(ndttrc)         ! leapfrog
(  107)          ENDIF
(  108)       ELSE
(  109)          rdttrc(:) =  rdttra(:) * FLOAT(ndttrc)      
(  110)       ENDIF
(  111) 
(  112)       DO jn = 1 , jptra
(  113) 
(  114) 	 ! Initialisation     
(  115) 	 zwd( 1 ,:,:) = 0.e0     ;     zwd(jpi,:,:) = 0.e0
(  116) 	 zws( 1 ,:,:) = 0.e0     ;     zws(jpi,:,:) = 0.e0
(  117) 	 zwi( 1 ,:,:) = 0.e0     ;     zwi(jpi,:,:) = 0.e0
(  118) 	 zwt( 1 ,:,:) = 0.e0     ;     zwt(jpi,:,:) = 0.e0     
(  119)          zwt(  :,:,1) = 0.e0     ;     zwt(  :,:,jpk) = 0.e0
(  120)          !                                          
(  121)          ! 0. Matrix construction
(  122)          ! ----------------------
(  123) 
(  124)          ! Diagonal, inferior, superior
(  125)          ! (including the bottom boundary condition via avs masked
(  126)          DO jk = 1, jpkm1                                                     
(  127)             DO jj = 2, jpjm1                                      
(  128)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  129)                   zwi(ji,jj,jk) = - rdttrc(jk) * fstravs(ji,jj,jk  ) /( fse3t(ji,jj,jk) * fse3w(ji,jj,jk  ) )
(  130)                   zws(ji,jj,jk) = - rdttrc(jk) * fstravs(ji,jj,jk+1) /( fse3t(ji,jj,jk) * fse3w(ji,jj,jk+1) )
(  131)                   zwd(ji,jj,jk) = 1. - zwi(ji,jj,jk) - zws(ji,jj,jk)
(  132)                END DO
(  133)             END DO
(  134)          END DO
(  135) 
(  136)          ! Surface boudary conditions
(  137)          DO jj = 2, jpjm1        
(  138)             DO ji = fs_2, fs_jpim1
(  139)                zwi(ji,jj,1) = 0.e0
(  140)                zwd(ji,jj,1) = 1. - zws(ji,jj,1)
(  141)             END DO
(  142)          END DO
(  143)          
(  144)          ! Second member construction
(  145)          DO jk = 1, jpkm1
(  146)             DO jj = 2, jpjm1     
(  147)                DO ji = fs_2, fs_jpim1
(  148)                   zwy(ji,jj,jk) = trb(ji,jj,jk,jn) + rdttrc(jk) * tra(ji,jj,jk,jn)
(  149)                END DO
(  150)             END DO
(  151)          END DO
(  152)          
(  153)  
(  154) 	! Matrix inversion from the first level
(  155) 	ikst = 1
(  156) 
(  157) #   include "zdf.matrixsolver.vopt.h90"        
(  158)  
(  159)          
(  160) #if defined key_trc_diatrd
(  161)          ! Compute and save the vertical diffusive of tracers trends
(  162) #  if defined key_trc_ldfiso






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 4

(  163)          DO jk = 1, jpkm1
(  164)             DO jj = 2, jpjm1
(  165)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  166)                   ztra = ( zwx(ji,jj,jk) - trb(ji,jj,jk,jn) ) / rdttrc(jk)
(  167)                   IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),6) = ztra - tra(ji,jj,jk,jn) + trtrd(ji,jj,jk,ikeep(jn),6)
(  168)                END DO
(  169)             END DO
(  170)          END DO
(  171) #  else
(  172)          DO jk = 1, jpkm1
(  173)             DO jj = 2, jpjm1
(  174)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  175)                   ztra = ( zwx(ji,jj,jk) - trb(ji,jj,jk,jn) ) / rdttrc(jk)
(  176)                   IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),6) = ztra - tra(ji,jj,jk,jn)
(  177)                END DO
(  178)             END DO
(  179)          END DO
(  180) #  endif
(  181) #endif  
(  182)          ! Save the masked passive tracer after in tra
(  183)          ! (c a u t i o n:  tracer not its trend, Leap-frog scheme done
(  184)          !                  it will not be done in tranxt)
(  185)          DO jk = 1, jpkm1
(  186)             DO jj = 2, jpjm1
(  187)                DO ji = fs_2, fs_jpim1
(  188)                   tra(ji,jj,jk,jn) = zwx(ji,jj,jk) * tmask(ji,jj,jk)
(  189)                END DO
(  190)             END DO
(  191)          END DO
(  192) 
(  193)          IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
(  194)             ztrd(:,:,:,:) = 0.
(  195)             DO jk = 1, jpkm1
(  196)                DO jj = 2, jpjm1
(  197)                   DO ji = fs_2, fs_jpim1
(  198)                      ztrd(ji,jj,jk,jn) = ( zwx(ji,jj,jk) - trb(ji,jj,jk,jn) ) / rdttrc(jk)
(  199)                   END DO
(  200)                END DO
(  201)             END DO
(  202)          ENDIF
(  203) 
(  204)       END DO
(  205) 
(  206)       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
(  207)          WRITE(charout, FMT="('zdf - imp')")
(  208)          CALL prt_ctl_trc_info(charout)
(  209)          CALL prt_ctl_trc(tab4d=ztrd, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
(  210)       ENDIF
(  211) 
(  212)    END SUBROUTINE trc_zdf_imp
(  213) 
(  214) #else
(  215)    !!----------------------------------------------------------------------
(  216)    !!   Dummy module :                      NO passive tracer
(  217)    !!----------------------------------------------------------------------
(  218) CONTAINS
(  219)    SUBROUTINE trc_zdf_imp (kt )              ! Empty routine
(  220)       INTEGER, INTENT(in) :: kt






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 5

(  221) !      WRITE(*,*) 'trc_zdf_imp: You should not have seen this print! error?', kt
(  222)    END SUBROUTINE trc_zdf_imp
(  223) #endif
(  224)    
(  225) !!==============================================================================
(  226) END MODULE trczdf_imp


























































PGF90 (Version     10.5)          09/23/2010  09:49:23      page 6

(    1) # 1 "trczdf_imp.F90"
(    1) MODULE trczdf_imp
(    2)    !!==============================================================================
(    3)    !!                    ***  MODULE  trczdf_imp  ***
(    4)    !! Ocean passive tracers:  vertical component of the tracer mixing trend
(    5)    !!==============================================================================
(    6) # 215
(  215)    !!----------------------------------------------------------------------
(  216)    !!   Dummy module :                      NO passive tracer
(  217)    !!----------------------------------------------------------------------
(  218) CONTAINS
(  219)    SUBROUTINE trc_zdf_imp (kt )              ! Empty routine
(  220)       INTEGER, INTENT(in) :: kt
(  221) !      WRITE(*,*) 'trc_zdf_imp: You should not have seen this print! error?', kt
(  222)    END SUBROUTINE trc_zdf_imp
(  223) # 224
(  224)    
(  225) !!==============================================================================
(  226) END MODULE trczdf_imp
