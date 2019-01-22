


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

Filename: trczdf_iso_vopt.F90

(    1) MODULE trczdf_iso_vopt
(    2)    !!==============================================================================
(    3)    !!                 ***  MODULE  trczdf_iso_vopt  ***
(    4)    !! Ocean passive tracers:  vertical component of the tracer mixing trend
(    5)    !!==============================================================================
(    6) #if defined key_passivetrc && ( defined key_ldfslp   ||   defined key_esopa )
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ldfslp'                  rotation of the lateral mixing tensor
(    9)    !!----------------------------------------------------------------------
(   10)    !!   trc_zdf_iso_vopt : Update the tracer trend with the vertical part of 
(   11)    !!                  the isopycnal or geopotential s-coord. operator and
(   12)    !!                  the vertical diffusion. vector optimization, use
(   13)    !!                  k-j-i loops.
(   14)    !!   trc_zdf_iso  :
(   15)    !!   trc_zdf_zdf  :
(   16)    !!----------------------------------------------------------------------
(   17)    !! * Modules used
(   18)    USE oce_trc         ! ocean dynamics and tracers variables
(   19)    USE trc             ! ocean passive tracers variables 
(   20)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   21)    USE trctrp_lec      ! passive tracers transport
(   22)    USE prtctl_trc          ! Print control for debbuging
(   23) 
(   24)    IMPLICIT NONE
(   25)    PRIVATE
(   26) 
(   27)    !! * Routine accessibility
(   28)    PUBLIC trc_zdf_iso_vopt   !  routine called by step.F90
(   29) 
(   30)    !! * Module variables
(   31)    REAL(wp), DIMENSION(jpk) ::  &
(   32)       rdttrc                          ! vertical profile of 2 x time-step
(   33) 
(   34)    !! * Substitutions
(   35) #  include "passivetrc_substitute.h90"
(   36)    !!----------------------------------------------------------------------
(   37)    !!   TOP 1.0 , LOCEAN-IPSL (2005) 
(   38)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trczdf_iso_vopt.F90,v 1.9 2006/04/10 15:38:55 opalod Exp $ 
(   39)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   40)    !!----------------------------------------------------------------------
(   41) 
(   42) CONTAINS
(   43)    
(   44)    SUBROUTINE trc_zdf_iso_vopt( kt )
(   45)       !!----------------------------------------------------------------------
(   46)       !!                  ***  ROUTINE trc_zdf_iso_vopt  ***






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 2

(   47)       !!
(   48)       !! ** Purpose :
(   49)       !! ** Method  :
(   50)       !! ** Action  :
(   51)       !!
(   52)       !! History :
(   53)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(   54)       !!   9.0  !  04-03  (C. Ethe)   adapted for passive tracers
(   55)       !!---------------------------------------------------------------------
(   56)       !! * Arguments
(   57)       INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
(   58)       CHARACTER (len=22) :: charout
(   59)       !!---------------------------------------------------------------------
(   60) 
(   61)       IF( kt == nittrc000 ) THEN
(   62)          IF(lwp)WRITE(numout,*)
(   63)          IF(lwp)WRITE(numout,*) 'trc_zdf_iso_vopt : vertical mixing computation'
(   64)          IF(lwp)WRITE(numout,*) '~~~~~~~~~~~~~~~~   is  iso-neutral diffusion : implicit vertical time stepping'
(   65) #if defined key_trcldf_eiv && defined key_diaeiv 
(   66)          w_trc_eiv(:,:,:) = 0.e0
(   67) #endif
(   68)       ENDIF
(   69) 
(   70) 
(   71)       ! I. vertical extra-diagonal part of the rotated tensor
(   72)       ! -----------------------------------------------------
(   73) 
(   74)       CALL trc_zdf_iso
(   75) 
(   76)       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
(   77)          WRITE(charout, FMT="('zdf - 1')")
(   78)          CALL prt_ctl_trc_info(charout)
(   79)          CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
(   80)       ENDIF
(   81) 
(   82)       ! II. vertical diffusion (including the vertical diagonal part of the rotated tensor)
(   83)       ! ----------------------
(   84) 
(   85)       CALL trc_zdf_zdf( kt )
(   86) 
(   87)       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
(   88)          WRITE(charout, FMT="('zdf - 2')")
(   89)          CALL prt_ctl_trc_info(charout)
(   90)          CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
(   91)       ENDIF
(   92) 
(   93)    END SUBROUTINE trc_zdf_iso_vopt
(   94) 
(   95) 
(   96)    SUBROUTINE trc_zdf_zdf( kt )
(   97)       !!----------------------------------------------------------------------
(   98)       !!                 ***  ROUTINE trc_zdf_zdf  ***
(   99)       !!                    
(  100)       !! ** Purpose :   Compute the trend due to the vertical tracer diffusion
(  101)       !!     including the vertical component of lateral mixing (only for 2nd
(  102)       !!     order operator, for fourth order it is already computed and add
(  103)       !!     to the general trend in traldf.F) and add it to the general trend
(  104)       !!     of the tracer equations.






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 3

(  105)       !!
(  106)       !! ** Method  :   The vertical component of the lateral diffusive trends
(  107)       !!      is provided by a 2nd order operator rotated along neural or geo-
(  108)       !!      potential surfaces to which an eddy induced advection can be 
(  109)       !!      added. It is computed using before fields (forward in time) and 
(  110)       !!      isopycnal or geopotential slopes computed in routine ldfslp.
(  111)       !!
(  112)       !!      Second part: vertical trend associated with the vertical physics
(  113)       !!      ===========  (including the vertical flux proportional to dk[t]
(  114)       !!                  associated with the lateral mixing, through the
(  115)       !!                  update of avt)
(  116)       !!      The vertical diffusion of tracers  is given by:
(  117)       !!             difft = dz( avt dz(t) ) = 1/e3t dk+1( avt/e3w dk(t) )
(  118)       !!      It is computed using a backward time scheme (t=tra).
(  119)       !!      Surface and bottom boundary conditions: no diffusive flux on
(  120)       !!      both tracers (bottom, applied through the masked field avt).
(  121)       !!      Add this trend to the general trend tra :
(  122)       !!         tra = tra + dz( avt dz(t) )
(  123)       !!         (tra = tra + dz( avs dz(t) ) if lk_trc_zdfddm=T )
(  124)       !!
(  125)       !!      Third part: recover avt resulting from the vertical physics
(  126)       !!      ==========  alone, for further diagnostics (for example to
(  127)       !!                  compute the turbocline depth in diamld).
(  128)       !!         avt = zavt
(  129)       !!         (avs = zavs if lk_trc_zdfddm=T )
(  130)       !!
(  131)       !!      'key_trdtra' defined: trend saved for futher diagnostics.
(  132)       !!
(  133)       !!      macro-tasked on vertical slab (jj-loop)
(  134)       !!
(  135)       !! ** Action  : - Update tra with before vertical diffusion trend
(  136)       !!              - Save the trend in trtrd  ('key_trc_diatrd')
(  137)       !!
(  138)       !! History :
(  139)       !!   6.0  !  90-10  (B. Blanke)  Original code
(  140)       !!   7.0  !  91-11 (G. Madec)
(  141)       !!        !  92-06 (M. Imbard) correction on tracer trend loops
(  142)       !!        !  96-01 (G. Madec) statement function for e3
(  143)       !!        !  97-05 (G. Madec) vertical component of isopycnal
(  144)       !!        !  97-07 (G. Madec) geopotential diffusion in s-coord
(  145)       !!        !  98-03  (L. Bopp MA Foujols) passive tracer generalisation
(  146)       !!        !  00-05  (MA Foujols) add lbc for tracer trends
(  147)       !!        !  00-06  (O Aumont)  correct isopycnal scheme suppress
(  148)       !!        !                     avt multiple correction
(  149)       !!        !  00-08  (G. Madec)  double diffusive mixing
(  150)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(  151)       !!   9.0  !  04-03  (C. Ethe )  adapted for passive tracers
(  152)       !!---------------------------------------------------------------------
(  153)       !! * Modules used
(  154)       USE oce_trc, ONLY :   zwd   => ua,  &  ! ua, va used as
(  155)                             zws   => va      ! workspace
(  156)       !! * Arguments
(  157)       INTEGER, INTENT( in ) ::   kt          ! ocean time-step index
(  158) 
(  159)       !! * Local declarations
(  160)       INTEGER ::   ji, jj, jk,jn                ! dummy loop indices
(  161)       REAL(wp) ::   &
(  162)          zavi, zrhs                          ! temporary scalars






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 4

(  163)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  164)          zwi, zwt, zavsi                     ! temporary workspace arrays
(  165)       REAL(wp) ::    ztra              !temporary scalars
(  166) #  if defined key_trc_diatrd
(  167)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   ztrd
(  168) #  endif
(  169)       !!---------------------------------------------------------------------
(  170) 
(  171) 
(  172)       ! I. Local constant initialization
(  173)       ! --------------------------------
(  174)       ! ... time step = 2 rdttra ex
(  175)       IF( ln_trcadv_cen2 .OR. ln_trcadv_tvd ) THEN
(  176)          ! time step = 2 rdttra with Arakawa or TVD advection scheme
(  177)          IF( neuler == 0 .AND. kt == nittrc000 ) THEN
(  178)             rdttrc(:) =  rdttra(:) * FLOAT(ndttrc)             ! restarting with Euler time stepping
(  179)          ELSEIF( kt <= nittrc000 + 1 ) THEN
(  180)             rdttrc(:) = 2. * rdttra(:) * FLOAT(ndttrc)         ! leapfrog
(  181)          ENDIF
(  182)       ELSE
(  183)          rdttrc(:) =  rdttra(:) * FLOAT(ndttrc)      
(  184)       ENDIF
(  185) 
(  186)       DO jn = 1, jptra
(  187)          
(  188)          zwd( 1 ,:,:)=0.e0     ;     zwd(jpi,:,:)=0.e0
(  189)          zws( 1 ,:,:)=0.e0     ;     zws(jpi,:,:)=0.e0
(  190)          zwi( 1 ,:,:)=0.e0     ;     zwi(jpi,:,:)=0.e0
(  191) 
(  192)          zwt( 1 ,:,:)=0.e0     ;     zwt(jpi,:,:)=0.e0
(  193)          zwt(  :,:,1)=0.e0     ;     zwt(:,:,jpk)= 0.e0
(  194)          zavsi( 1 ,:,:)=0.e0   ;     zavsi(jpi,:,:)=0.e0 
(  195)          zavsi(  :,:,1)=0.e0   ;     zavsi(:,:,jpk)=0.e0
(  196) 
(  197) #  if defined key_trc_diatrd
(  198)          ! save the tra trend
(  199)          ztrd(:,:,:) = tra(:,:,:,jn)
(  200) #  endif
(  201) 
(  202)          ! II. Vertical trend associated with the vertical physics
(  203)          ! =======================================================
(  204)          !     (including the vertical flux proportional to dk[t] associated
(  205)          !      with the lateral mixing, through the avt update)
(  206)          !     dk[ avt dk[ (t,s) ] ] diffusive trends
(  207) 
(  208) 
(  209)          ! II.0 Matrix construction
(  210)          ! ------------------------        
(  211)          ! update and save of avt (and avs if double diffusive mixing)
(  212)          DO jk = 2, jpkm1
(  213)             DO jj = 2, jpjm1
(  214)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  215)                   zavi = fsahtw(ji,jj,jk) * (                 &   ! vertical mixing coef. due to lateral mixing
(  216)                      &                           wslpi(ji,jj,jk) * wslpi(ji,jj,jk)      &
(  217)                      &                         + wslpj(ji,jj,jk) * wslpj(ji,jj,jk)  )
(  218)                   zavsi(ji,jj,jk) = fstravs(ji,jj,jk) + zavi        ! dd mixing: zavsi = total vertical mixing coef. on trac
(  219) 
(  220)                END DO






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 5

(  221)             END DO
(  222)          END DO
(  223) 
(  224) 
(  225)          ! II.2 Vertical diffusion on tracer
(  226)          ! ---------------------------========
(  227) 
(  228)          ! Rebuild the Matrix as avt /= avs
(  229) 
(  230)          ! Diagonal, inferior, superior  (including the bottom boundary condition via avs masked)
(  231)          DO jk = 1, jpkm1
(  232)             DO jj = 2, jpjm1
(  233)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  234)                   zwi(ji,jj,jk) = - rdttrc(jk) * zavsi(ji,jj,jk  ) / ( fse3t(ji,jj,jk) * fse3w(ji,jj,jk  ) )
(  235)                   zws(ji,jj,jk) = - rdttrc(jk) * zavsi(ji,jj,jk+1) / ( fse3t(ji,jj,jk) * fse3w(ji,jj,jk+1) )
(  236)                   zwd(ji,jj,jk) = 1. - zwi(ji,jj,jk) - zws(ji,jj,jk)
(  237)                END DO
(  238)             END DO
(  239)          END DO
(  240) 
(  241)          ! Surface boudary conditions
(  242)          DO jj = 2, jpjm1
(  243)             DO ji = fs_2, fs_jpim1   ! vector opt.
(  244)                zwi(ji,jj,1) = 0.e0
(  245)                zwd(ji,jj,1) = 1. - zws(ji,jj,1)
(  246)             END DO
(  247)          END DO
(  248) 
(  249)          !! Matrix inversion from the first level
(  250)          !!----------------------------------------------------------------------
(  251)          !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
(  252)          !
(  253)          !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
(  254)          !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
(  255)          !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
(  256)          !        (        ...               )( ...  ) ( ...  )
(  257)          !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
(  258)          !
(  259)          !   m is decomposed in the product of an upper and lower triangular
(  260)          !   matrix
(  261)          !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
(  262)          !   The second member is in 2d array zwy
(  263)          !   The solution is in 2d array zwx
(  264)          !   The 3d arry zwt is a work space array
(  265)          !   zwy is used and then used as a work space array : its value is modified!
(  266) 
(  267)          ! first recurrence:   Tk = Dk - Ik Sk-1 / Tk-1   (increasing k)
(  268)          DO jj = 2, jpjm1
(  269)             DO ji = fs_2, fs_jpim1
(  270)                zwt(ji,jj,1) = zwd(ji,jj,1)
(  271)             END DO
(  272)          END DO
(  273)          DO jk = 2, jpkm1
(  274)             DO jj = 2, jpjm1
(  275)                DO ji = fs_2, fs_jpim1
(  276)                   zwt(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1)  /zwt(ji,jj,jk-1)
(  277)                END DO
(  278)             END DO






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 6

(  279)          END DO
(  280) 
(  281)          ! second recurrence:    Zk = Yk - Ik / Tk-1  Zk-1
(  282)          DO jj = 2, jpjm1
(  283)             DO ji = fs_2, fs_jpim1
(  284)                tra(ji,jj,1,jn) = trb(ji,jj,1,jn) + rdttrc(1) * tra(ji,jj,1,jn)
(  285)             END DO
(  286)          END DO
(  287)          DO jk = 2, jpkm1
(  288)             DO jj = 2, jpjm1
(  289)                DO ji = fs_2, fs_jpim1
(  290)                   zrhs = trb(ji,jj,jk,jn) + rdttrc(jk) * tra(ji,jj,jk,jn)   ! zrhs=right hand side
(  291)                   tra(ji,jj,jk,jn) = zrhs - zwi(ji,jj,jk) / zwt(ji,jj,jk-1) * tra(ji,jj,jk-1,jn)
(  292)                END DO
(  293)             END DO
(  294)          END DO
(  295) 
(  296)          ! third recurrence: Xk = (Zk - Sk Xk+1 ) / Tk
(  297)          ! Save the masked passive tracer after in tra
(  298)          ! (c a u t i o n: passive tracer not its trend, Leap-frog scheme done it will not be done in tranxt)
(  299)          DO jj = 2, jpjm1
(  300)             DO ji = fs_2, fs_jpim1
(  301)                tra(ji,jj,jpkm1,jn) = tra(ji,jj,jpkm1,jn) / zwt(ji,jj,jpkm1) * tmask(ji,jj,jpkm1)
(  302)             END DO
(  303)          END DO
(  304)          DO jk = jpk-2, 1, -1
(  305)             DO jj = 2, jpjm1
(  306)                DO ji = fs_2, fs_jpim1
(  307)                   tra(ji,jj,jk,jn) = ( tra(ji,jj,jk,jn) - zws(ji,jj,jk) * tra(ji,jj,jk+1,jn) ) / zwt(ji,jj,jk) * tmask(ji,jj
(  308)                END DO
(  309)             END DO
(  310)          END DO
(  311) 
(  312) #if defined key_trc_diatrd
(  313)          ! Compute and save the vertical diffusive passive tracer trends
(  314) #  if defined key_trc_ldfiso 
(  315)          DO jk = 1, jpkm1
(  316)             DO jj = 2, jpjm1
(  317)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  318)                   ztra = ( tra(ji,jj,jk,jn) - trb(ji,jj,jk,jn) ) / rdttrc(jk)
(  319)                   IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),6) = ztra - ztrd(ji,jj,jk) + trtrd(ji,jj,jk,ikeep(jn),6)
(  320)                END DO
(  321)             END DO
(  322)          END DO
(  323) #  else
(  324)          DO jk = 1, jpkm1
(  325)             DO jj = 2, jpjm1
(  326)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  327)                   ztra = ( tra(ji,jj,jk,jn) - trb(ji,jj,jk,jn) ) / rdttrc(jk)
(  328)                   IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),6) = ztra - ztrd(ji,jj,jk)
(  329)                END DO
(  330)             END DO
(  331)          END DO
(  332) #  endif
(  333) #endif
(  334) 
(  335)       END DO
(  336) 






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 7

(  337)    END SUBROUTINE trc_zdf_zdf
(  338) 
(  339) 
(  340)    SUBROUTINE trc_zdf_iso
(  341)       !!----------------------------------------------------------------------
(  342)       !!                  ***  ROUTINE trc_zdf_iso  ***
(  343)       !!
(  344)       !! ** Purpose :
(  345)       !!     Compute the trend due to the vertical tracer diffusion inclu-
(  346)       !!     ding the vertical component of lateral mixing (only for second
(  347)       !!     order operator, for fourth order it is already computed and
(  348)       !!     add to the general trend in traldf.F) and add it to the general
(  349)       !!     trend of the tracer equations.
(  350)       !!
(  351)       !! ** Method :
(  352)       !!         The vertical component of the lateral diffusive trends is
(  353)       !!      provided by a 2nd order operator rotated along neural or geopo-
(  354)       !!      tential surfaces to which an eddy induced advection can be added
(  355)       !!      It is computed using before fields (forward in time) and isopyc-
(  356)       !!      nal or geopotential slopes computed in routine ldfslp.
(  357)       !!
(  358)       !!      First part: vertical trends associated with the lateral mixing
(  359)       !!      ==========  (excluding the vertical flux proportional to dk[t] )
(  360)       !!      vertical fluxes associated with the rotated lateral mixing:
(  361)       !!         zftw =-aht {  e2t*wslpi di[ mi(mk(trb)) ]
(  362)       !!                     + e1t*wslpj dj[ mj(mk(trb)) ]  }
(  363)       !!      save avt coef. resulting from vertical physics alone in zavt:
(  364)       !!         zavt = avt
(  365)       !!      update and save in zavt the vertical eddy viscosity coefficient:
(  366)       !!         avt = avt + wslpi^2+wslj^2
(  367)       !!      add vertical Eddy Induced advective fluxes (lk_traldf_eiv=T):
(  368)       !!         zftw = zftw + { di[aht e2u mi(wslpi)]
(  369)       !!                    +dj[aht e1v mj(wslpj)] } mk(trb)
(  370)       !!      take the horizontal divergence of the fluxes:
(  371)       !!         difft = 1/(e1t*e2t*e3t) dk[ zftw ] 
(  372)       !!      Add this trend to the general trend tra :
(  373)       !!         tra = tra + difft
(  374)       !!
(  375)       !! ** Action :
(  376)       !!         Update tra arrays with the before vertical diffusion trend
(  377)       !!         Save in trtrd arrays the trends if 'key_trc_diatrd' defined
(  378)       !!
(  379)       !! History :
(  380)       !!   6.0  !  90-10  (B. Blanke)  Original code
(  381)       !!   7.0  !  91-11  (G. Madec)
(  382)       !!        !  92-06  (M. Imbard) correction on tracer trend loops
(  383)       !!        !  96-01  (G. Madec) statement function for e3
(  384)       !!        !  97-05  (G. Madec) vertical component of isopycnal
(  385)       !!        !  97-07  (G. Madec) geopotential diffusion in s-coord
(  386)       !!        !  98-03  (L. Bopp MA Foujols) passive tracer generalisation
(  387)       !!        !  00-05  (MA Foujols) add lbc for tracer trends
(  388)       !!        !  00-06  (O Aumont)  correct isopycnal scheme suppress
(  389)       !!        !                     avt multiple correction
(  390)       !!        !  00-08  (G. Madec)  double diffusive mixing
(  391)       !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
(  392)       !!   9.0  !  04-03  (C. Ethe )  adapted for passive tracers
(  393)       !!---------------------------------------------------------------------
(  394)       !! * Modules used






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 8

(  395)       USE oce_trc, ONLY :   zwx => ua,  &  ! use ua, va as
(  396)                             zwy => va      ! workspace arrays
(  397) 
(  398)       !! * Local declarations
(  399)       INTEGER ::   ji, jj, jk,jn       ! dummy loop indices
(  400) #if defined key_partial_steps
(  401)       INTEGER ::   iku, ikv
(  402) #endif
(  403)       REAL(wp) ::   &
(  404)          ztavg,                  &  ! temporary scalars
(  405)          zcoef0, zcoef3,         &  !    "         "
(  406)          zcoef4,                 &  !    "         "
(  407)          zbtr, zmku, zmkv,       &  !    "         "
(  408) #if defined key_trcldf_eiv
(  409)          zcoeg3,                 &  !    "         "
(  410)          zuwki, zvwki,           &  !    "         "
(  411)          zuwk, zvwk,             &  !    "         "
(  412) #endif
(  413)          ztav
(  414)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  415)          zwz, zwt, ztfw             ! temporary workspace arrays
(  416)       !!---------------------------------------------------------------------
(  417) 
(  418)       DO jn = 1, jptra
(  419) 
(  420)          ! 0. Local constant initialization
(  421)          ! --------------------------------
(  422)          ztavg = 0.e0
(  423) 
(  424)          zwx( 1 ,:,:)=0.e0     ;     zwx(jpi,:,:)=0.e0
(  425)          zwy( 1 ,:,:)=0.e0     ;     zwy(jpi,:,:)=0.e0
(  426)          zwz( 1 ,:,:)=0.e0     ;     zwz(jpi,:,:)=0.e0
(  427)          zwt( 1 ,:,:)=0.e0     ;     zwt(jpi,:,:)=0.e0
(  428)          ztfw( 1 ,:,:)=0.e0    ;     ztfw(jpi,:,:)=0.e0
(  429) 
(  430)          ! I. Vertical trends associated with lateral mixing
(  431)          ! -------------------------------------------------
(  432)          !    (excluding the vertical flux proportional to dk[t] )
(  433) 
(  434) 
(  435)          ! I.1 horizontal tracer gradient
(  436)          ! ------------------------------
(  437) 
(  438)          DO jk = 1, jpkm1
(  439)             DO jj = 1, jpjm1
(  440)                DO ji = 1, fs_jpim1   ! vector opt.
(  441)                   ! i-gradient of passive tracer at ji
(  442)                   zwx (ji,jj,jk) = ( trb(ji+1,jj,jk,jn)-trb(ji,jj,jk,jn) ) * umask(ji,jj,jk)
(  443)                   ! j-gradient of passive tracer at jj
(  444)                   zwy (ji,jj,jk) = ( trb(ji,jj+1,jk,jn)-trb(ji,jj,jk,jn) ) * vmask(ji,jj,jk)
(  445)                END DO
(  446)             END DO
(  447)          END DO
(  448) #  if defined key_partial_steps
(  449)          ! partial steps correction at the bottom ocean level 
(  450)          DO jj = 1, jpjm1
(  451)             DO ji = 1, fs_jpim1   ! vector opt.
(  452)                ! last ocean level






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 9

(  453)                iku  = MIN( mbathy(ji,jj), mbathy(ji+1,jj  ) ) - 1
(  454)                ikv  = MIN( mbathy(ji,jj), mbathy(ji  ,jj+1) ) - 1
(  455)                ! i-gradient of passive tracer
(  456)                zwx (ji,jj,iku) = gtru(ji,jj,jn)
(  457)                ! j-gradient of passive tracer
(  458)                zwy (ji,jj,ikv) = gtrv(ji,jj,jn)  
(  459)             END DO
(  460)          END DO
(  461) #endif
(  462) 
(  463) 
(  464)          ! I.2 Vertical fluxes
(  465)          ! -------------------
(  466) 
(  467)          ! Surface and bottom vertical fluxes set to zero
(  468)          ztfw(:,:, 1 ) = 0.e0
(  469)          ztfw(:,:,jpk) = 0.e0
(  470) 
(  471)          ! interior (2=<jk=<jpk-1)
(  472)          DO jk = 2, jpkm1
(  473)             DO jj = 2, jpjm1
(  474)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  475)                   zcoef0 = - fsahtw(ji,jj,jk) * tmask(ji,jj,jk)
(  476) 
(  477)                   zmku = 1./MAX(   umask(ji  ,jj,jk-1) + umask(ji-1,jj,jk)      &
(  478)                      &           + umask(ji-1,jj,jk-1) + umask(ji  ,jj,jk), 1.  )
(  479) 
(  480)                   zmkv = 1./MAX(   vmask(ji,jj  ,jk-1) + vmask(ji,jj-1,jk)      &
(  481)                      &           + vmask(ji,jj-1,jk-1) + vmask(ji,jj  ,jk), 1.  )
(  482) 
(  483)                   zcoef3 = zcoef0 * e2t(ji,jj) * zmku * wslpi (ji,jj,jk)
(  484)                   zcoef4 = zcoef0 * e1t(ji,jj) * zmkv * wslpj (ji,jj,jk)
(  485) 
(  486)                   ztfw(ji,jj,jk) = zcoef3 * (   zwx(ji  ,jj  ,jk-1) + zwx(ji-1,jj  ,jk)      &
(  487)                      &                        + zwx(ji-1,jj  ,jk-1) + zwx(ji  ,jj  ,jk)  )   &
(  488)                      &           + zcoef4 * (   zwy(ji  ,jj  ,jk-1) + zwy(ji  ,jj-1,jk)      &
(  489)                      &                        + zwy(ji  ,jj-1,jk-1) + zwy(ji  ,jj  ,jk)  )
(  490)                END DO
(  491)             END DO
(  492)          END DO
(  493) 
(  494) #if defined key_trcldf_eiv
(  495)          !                              ! ---------------------------------------!
(  496)          !                              ! Eddy induced vertical advective fluxes !
(  497)          !                              ! ---------------------------------------!
(  498)          zwx(:,:, 1 ) = 0.e0
(  499)          zwx(:,:,jpk) = 0.e0
(  500) 
(  501)          DO jk = 2, jpkm1
(  502)             DO jj = 2, jpjm1
(  503)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  504) #   if defined key_traldf_c2d || defined key_traldf_c3d
(  505)                   zuwki = ( wslpi(ji,jj,jk) + wslpi(ji-1,jj,jk) )   &
(  506)                      &  * fsaeitru(ji-1,jj,jk) * e2u(ji-1,jj) * umask(ji-1,jj,jk)
(  507)                   zuwk  = ( wslpi(ji,jj,jk) + wslpi(ji+1,jj,jk) )   &
(  508)                      &  * fsaeitru(ji  ,jj,jk) * e2u(ji  ,jj) * umask(ji  ,jj,jk)
(  509)                   zvwki = ( wslpj(ji,jj,jk) + wslpj(ji,jj-1,jk) )   &
(  510)                      &  * fsaeitrv(ji,jj-1,jk) * e1v(ji,jj-1) * vmask(ji,jj-1,jk)






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 10

(  511)                   zvwk  = ( wslpj(ji,jj,jk) + wslpj(ji,jj+1,jk) )   &
(  512)                      &  * fsaeitrv(ji,jj  ,jk) * e1v(ji  ,jj) * vmask(ji  ,jj,jk)
(  513) 
(  514)                   zcoeg3 = + 0.25 * tmask(ji,jj,jk) * ( zuwk - zuwki + zvwk - zvwki )
(  515) #   else
(  516)                   zuwki = ( wslpi(ji,jj,jk) + wslpi(ji-1,jj,jk) )   &
(  517)                      &  * e2u(ji-1,jj) * umask(ji-1,jj,jk)
(  518)                   zuwk  = ( wslpi(ji,jj,jk) + wslpi(ji+1,jj,jk) )   &
(  519)                      &  * e2u(ji  ,jj) * umask(ji  ,jj,jk)
(  520)                   zvwki = ( wslpj(ji,jj,jk) + wslpj(ji,jj-1,jk) )   &
(  521)                      &  * e1v(ji,jj-1) * vmask(ji,jj-1,jk)
(  522)                   zvwk  = ( wslpj(ji,jj,jk) + wslpj(ji,jj+1,jk) )   &
(  523)                      &  * e1v(ji  ,jj) * vmask(ji  ,jj,jk)
(  524) 
(  525)                   zcoeg3 = + 0.25 * tmask(ji,jj,jk) * fsaeiw(ji,jj,jk)   &
(  526)                      &            * ( zuwk - zuwki + zvwk - zvwki )
(  527) #   endif
(  528)                   zwx(ji,jj,jk) = + zcoeg3 * ( trb(ji,jj,jk,jn) + trb(ji,jj,jk-1,jn) )
(  529) 
(  530)                   ztfw(ji,jj,jk) = ztfw(ji,jj,jk) + zwx(ji,jj,jk)
(  531) #   if defined key_diaeiv
(  532)                   w_trc_eiv(ji,jj,jk) = -2. * zcoeg3 / ( e1t(ji,jj)*e2t(ji,jj) )
(  533) #   endif
(  534)                END DO
(  535)             END DO
(  536)          END DO
(  537) #endif
(  538) 
(  539)          ! I.5 Divergence of vertical fluxes added to the general tracer trend
(  540)          ! -------------------------------------------------------------------
(  541) 
(  542)          DO jk = 1, jpkm1
(  543)             DO jj = 2, jpjm1
(  544)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  545)                   zbtr =  1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
(  546)                   ztav = (  ztfw(ji,jj,jk) - ztfw(ji,jj,jk+1)  ) * zbtr
(  547)                   tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + ztav
(  548) #if defined key_trc_diatrd
(  549) #   if defined key_trcldf_eiv
(  550)                   ztavg = ( zwx(ji,jj,jk) - zwx(ji,jj,jk+1) ) * zbtr
(  551)                   !  WARNING trtrd(ji,jj,jk,7) used for vertical gent velocity trend  not for damping !!!
(  552)                   IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),7) = ztavg
(  553) #   endif
(  554)                   IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),6) = ztav - ztavg
(  555) #endif
(  556)                END DO
(  557)             END DO
(  558)          END DO
(  559) 
(  560)       END DO
(  561) 
(  562)    END SUBROUTINE trc_zdf_iso
(  563) 
(  564) #else
(  565)    !!----------------------------------------------------------------------
(  566)    !!   Dummy module :             NO rotation of the lateral mixing tensor
(  567)    !!----------------------------------------------------------------------
(  568) CONTAINS






PGF90 (Version     10.5)          09/23/2010  09:49:23      page 11

(  569)    SUBROUTINE trc_zdf_iso_vopt( kt )              ! empty routine
(  570) !      WRITE(*,*) 'trc_zdf_iso_vopt: You should not have seen this print! error?', kt
(  571)    END SUBROUTINE trc_zdf_iso_vopt
(  572) #endif
(  573) 
(  574)    !!==============================================================================
(  575) END MODULE trczdf_iso_vopt

























































PGF90 (Version     10.5)          09/23/2010  09:49:23      page 12

(    1) # 1 "trczdf_iso_vopt.F90"
(    1) MODULE trczdf_iso_vopt
(    2)    !!==============================================================================
(    3)    !!                 ***  MODULE  trczdf_iso_vopt  ***
(    4)    !! Ocean passive tracers:  vertical component of the tracer mixing trend
(    5)    !!==============================================================================
(    6) # 565
(  565)    !!----------------------------------------------------------------------
(  566)    !!   Dummy module :             NO rotation of the lateral mixing tensor
(  567)    !!----------------------------------------------------------------------
(  568) CONTAINS
(  569)    SUBROUTINE trc_zdf_iso_vopt( kt )              ! empty routine
(  570) !      WRITE(*,*) 'trc_zdf_iso_vopt: You should not have seen this print! error?', kt
(  571)    END SUBROUTINE trc_zdf_iso_vopt
(  572) 
(  573) # 574
(  574)    !!==============================================================================
(  575) END MODULE trczdf_iso_vopt
