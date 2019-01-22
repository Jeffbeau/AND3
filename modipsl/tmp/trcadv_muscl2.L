


PGF90 (Version     10.5)          09/23/2010  09:49:18      page 1

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

Filename: trcadv_muscl2.F90

(    1) MODULE trcadv_muscl2
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  trcadv_muscl2  ***
(    4)    !! Ocean passive tracers:  horizontal & vertical advective trend
(    5)    !!==============================================================================
(    6) #if defined key_passivetrc
(    7)    !!----------------------------------------------------------------------
(    8)    !!   tra_adv_muscl2 : update the tracer trend with the horizontal
(    9)    !!                    and vertical advection trends using MUSCL2 scheme
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce_trc         ! ocean dynamics and active tracers variables
(   13)    USE trc             ! ocean passive tracers variables
(   14)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   15)    USE trcbbl          ! advective passive tracers in the BBL
(   16)    USE prtctl_trc          ! Print control for debbuging
(   17) 
(   18)    IMPLICIT NONE
(   19)    PRIVATE
(   20) 
(   21)    !! * Accessibility
(   22)    PUBLIC trc_adv_muscl2        ! routine called by trcstp.F90
(   23) 
(   24)    !! * Module variable
(   25)    REAL(wp), DIMENSION(jpk) ::   &
(   26)       rdttrc                     ! vertical profile of tracer time-step
(   27) 
(   28)    !! * Substitutions
(   29) #  include "passivetrc_substitute.h90"
(   30)    !!----------------------------------------------------------------------
(   31)    !!   TOP 1.0 , LOCEAN-IPSL (2005) 
(   32)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trcadv_muscl2.F90,v 1.11 2006/04/10 15:38:54 opalod Exp $ 
(   33)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   34)    !!----------------------------------------------------------------------
(   35) 
(   36) CONTAINS
(   37) 
(   38)    SUBROUTINE trc_adv_muscl2( kt )
(   39)       !!----------------------------------------------------------------------
(   40)       !!                   ***  ROUTINE trc_adv_muscl2  ***
(   41)       !!
(   42)       !! ** Purpose :   Compute the now trend due to total advection of passi-
(   43)       !!      ve tracer using a MUSCL scheme (Monotone Upstream-
(   44)       !!      Centered Scheme for Conservation Laws) and add it to the general
(   45)       !!      tracer trend.
(   46)       !!






PGF90 (Version     10.5)          09/23/2010  09:49:18      page 2

(   47)       !! ** Method : MUSCL scheme plus centered scheme at ocean boundaries
(   48)       !!
(   49)       !! ** Action : - update tra with the now advective tracer trends
(   50)       !!             - save trends in trtrd ('key_trc_diatrd')
(   51)       !!
(   52)       !! References :                
(   53)       !!      Estubier, A., and M. Levy, Notes Techn. Pole de Modelisation
(   54)       !!	IPSL, Sept. 2000 (http://www.lodyc.jussieu.fr/opa)
(   55)       !!
(   56)       !! History :
(   57)       !!        !  06-00  (A.Estublier)  for passive tracers
(   58)       !!   9.0  !  03-04  (C. Ethe, G. Madec)  F90: Free form and module
(   59)       !!----------------------------------------------------------------------
(   60)       !! * modules used
(   61) #if defined key_trcbbl_adv
(   62)       USE oce_trc            , zun => ua,  &  ! use ua as workspace
(   63)          &                     zvn => va      ! use va as workspace
(   64)       REAL(wp), DIMENSION(jpi,jpj,jpk) :: zwn
(   65) #else
(   66)       USE oce_trc            , zun => un,  &  ! When no bbl, zun == un
(   67)                                zvn => vn,  &  !              zvn == vn
(   68)                                zwn => wn      !              zwn == wn
(   69) #endif
(   70)       !! * Arguments
(   71)       INTEGER, INTENT( in ) ::   kt         ! ocean time-step
(   72) 
(   73)       !! * Local declarations
(   74)       INTEGER ::   ji, jj, jk,jn               ! dummy loop indices
(   75)       REAL(wp), DIMENSION (jpi,jpj,jpk) ::   &
(   76)          zt1, zt2, ztp1, ztp2
(   77)       REAL(wp) ::   zu, zv, zw, zeu, zev, zew, zbtr, ztra
(   78)       REAL(wp) ::   z0u, z0v, z0w
(   79)       REAL(wp) ::   zzt1, zzt2, zalpha
(   80) 
(   81) #if defined key_trc_diatrd
(   82)       REAL(wp) ::   ztai, ztaj
(   83)       REAL(wp) ::   zfui, zfvj
(   84) #endif
(   85)       CHARACTER (len=22) :: charout
(   86)       !!----------------------------------------------------------------------
(   87) 
(   88)       IF( kt == nittrc000 .AND. lwp ) THEN
(   89)          WRITE(numout,*)
(   90)          WRITE(numout,*) 'trc_adv_muscl2 : MUSCL2 advection scheme'
(   91)          WRITE(numout,*) '~~~~~~~~~~~~~~~'
(   92)          rdttrc(:) =  rdttra(:) * FLOAT(ndttrc)
(   93)       ENDIF
(   94)       
(   95) #if defined key_trcbbl_adv        
(   96)       ! Advective bottom boundary layer
(   97)       ! -------------------------------
(   98)       zun(:,:,:) = un (:,:,:) - u_trc_bbl(:,:,:)
(   99)       zvn(:,:,:) = vn (:,:,:) - v_trc_bbl(:,:,:)
(  100)       zwn(:,:,:) = wn (:,:,:) + w_trc_bbl( :,:,:)
(  101) #endif
(  102) 
(  103) 
(  104)       DO jn = 1, jptra






PGF90 (Version     10.5)          09/23/2010  09:49:18      page 3

(  105) 
(  106)          ! I. Horizontal advective fluxes
(  107)          ! ------------------------------
(  108) 
(  109)          ! first guess of the slopes
(  110)          ! interior values
(  111)          DO jk = 1, jpkm1
(  112)             DO jj = 1, jpjm1      
(  113)                DO ji = 1, fs_jpim1   ! vector opt.
(  114)                   zt1(ji,jj,jk) = umask(ji,jj,jk) * ( trb(ji+1,jj,jk,jn) - trb(ji,jj,jk,jn) )
(  115)                   zt2(ji,jj,jk) = vmask(ji,jj,jk) * ( trb(ji,jj+1,jk,jn) - trb(ji,jj,jk,jn) )
(  116)                END DO
(  117)             END DO
(  118)          END DO
(  119)          ! bottom values
(  120)          zt1(:,:,jpk) = 0.e0
(  121)          zt2(:,:,jpk) = 0.e0
(  122) 
(  123)          ! lateral boundary conditions on zt1, zt2  (changed sign)
(  124)          CALL lbc_lnk( zt1, 'U', -1. )
(  125)          CALL lbc_lnk( zt2, 'V', -1. )
(  126) 
(  127)          ! Slopes
(  128)          ! interior values
(  129)          DO jk = 1, jpkm1
(  130)             DO jj = 2, jpj
(  131)                DO ji = fs_2, jpi   ! vector opt.
(  132)                   ztp1(ji,jj,jk) =                    ( zt1(ji,jj,jk) + zt1(ji-1,jj  ,jk) )   &
(  133)                      &           * ( 0.25 + SIGN( 0.25, zt1(ji,jj,jk) * zt1(ji-1,jj  ,jk) ) )
(  134)                   ztp2(ji,jj,jk) =                    ( zt2(ji,jj,jk) + zt2(ji  ,jj-1,jk) )   &
(  135)                      &           * ( 0.25 + SIGN( 0.25, zt2(ji,jj,jk) * zt2(ji  ,jj-1,jk) ) )
(  136)                END DO
(  137)             END DO
(  138)          END DO
(  139)          ! bottom values
(  140)          ztp1(:,:,jpk) = 0.e0 
(  141)          ztp2(:,:,jpk) = 0.e0
(  142) 
(  143)          ! Slopes limitation
(  144)          DO jk = 1, jpkm1
(  145)             DO jj = 2, jpj
(  146)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  147)                   ztp1(ji,jj,jk) = SIGN( 1., ztp1(ji,jj,jk) )   &
(  148)                      &           * MIN(    ABS( ztp1(ji  ,jj,jk) ),   &
(  149)                      &                  2.*ABS( zt1 (ji-1,jj,jk) ),   &
(  150)                      &                  2.*ABS( zt1 (ji  ,jj,jk) ) )
(  151) 
(  152)                   ztp2(ji,jj,jk) = SIGN( 1., ztp2(ji,jj,jk) )   &
(  153)                      &           * MIN(    ABS( ztp2(ji,jj  ,jk) ),   &
(  154)                      &                  2.*ABS( zt2 (ji,jj-1,jk) ),   &
(  155)                      &                  2.*ABS( zt2 (ji,jj  ,jk) ) )
(  156)                END DO
(  157)             END DO
(  158)          END DO
(  159) 
(  160)          ! Advection terms
(  161)          ! interior values
(  162)          DO jk = 1, jpkm1






PGF90 (Version     10.5)          09/23/2010  09:49:18      page 4

(  163)             DO jj = 2, jpjm1      
(  164)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  165)                   ! volume fluxes
(  166) #if defined key_s_coord || defined key_partial_steps
(  167)                   zeu = e2u(ji,jj) * fse3u(ji,jj,jk) * zun(ji,jj,jk)
(  168)                   zev = e1v(ji,jj) * fse3v(ji,jj,jk) * zvn(ji,jj,jk)
(  169) #else
(  170)                   zeu = e2u(ji,jj) * zun(ji,jj,jk)
(  171)                   zev = e1v(ji,jj) * zvn(ji,jj,jk)
(  172) #endif
(  173)                   ! MUSCL fluxes
(  174)                   z0u = SIGN( 0.5, zun(ji,jj,jk) )            
(  175)                   zalpha = 0.5 - z0u
(  176)                   zu  = z0u - 0.5 * zun(ji,jj,jk) * rdttrc(jk) / e1u(ji,jj)
(  177)                   zzt1 = trb(ji+1,jj,jk,jn) + zu*ztp1(ji+1,jj,jk)
(  178)                   zzt2 = trb(ji  ,jj,jk,jn) + zu*ztp1(ji  ,jj,jk)
(  179)                   zt1(ji,jj,jk) = zeu * ( zalpha * zzt1 + (1.-zalpha) * zzt2 )
(  180) 
(  181)                   z0v = SIGN( 0.5, zvn(ji,jj,jk) )            
(  182)                   zalpha = 0.5 - z0v
(  183)                   zv  = z0v - 0.5 * zvn(ji,jj,jk) * rdttrc(jk) / e2v(ji,jj)
(  184)                   zzt1 = trb(ji,jj+1,jk,jn) + zv*ztp2(ji,jj+1,jk)
(  185)                   zzt2 = trb(ji,jj  ,jk,jn) + zv*ztp2(ji,jj  ,jk)
(  186)                   zt2(ji,jj,jk) = zev * ( zalpha * zzt1 + (1.-zalpha) * zzt2 )
(  187) 
(  188)                END DO
(  189)             END DO
(  190)          END DO
(  191) 
(  192) 
(  193)          DO jk = 1, jpkm1
(  194)             DO jj = 2, jpjm1
(  195)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  196) #if defined key_s_coord || defined key_partial_steps
(  197)                   zev = e1v(ji,jj) * fse3v(ji,jj,jk)
(  198)                   IF( umask(ji,jj,jk) == 0. ) THEN
(  199)                      IF( zun(ji+1,jj,jk) > 0. .AND. ji /= jpi ) THEN
(  200)                         zt1(ji+1,jj,jk) = e2u(ji+1,jj)* fse3u(ji+1,jj,jk)   &
(  201)                            &            * zun(ji+1,jj,jk) * ( trb(ji+1,jj,jk,jn) + trb(ji+2,jj,jk,jn) ) * 0.5
(  202)                      ENDIF
(  203)                      IF( zun(ji-1,jj,jk) < 0. ) THEN
(  204)                         zt1(ji-1,jj,jk) = e2u(ji-1,jj)* fse3u(ji-1,jj,jk)   &
(  205)                            &            * zun(ji-1,jj,jk) * ( trb(ji-1,jj,jk,jn) + trb(ji  ,jj,jk,jn) ) * 0.5
(  206)                      ENDIF
(  207)                   ENDIF
(  208)                   IF( vmask(ji,jj,jk) == 0. ) THEN
(  209)                      IF( zvn(ji,jj+1,jk) > 0. .AND. jj /= jpj ) THEN
(  210)                         zt2(ji,jj+1,jk) = e1v(ji,jj+1) * fse3v(ji,jj+1,jk)   &
(  211)                            &            * zvn(ji,jj+1,jk) * ( trb(ji,jj+1,jk,jn) + trb(ji,jj+2,jk,jn) ) * 0.5
(  212)                      ENDIF
(  213)                      IF( zvn(ji,jj-1,jk) < 0. ) THEN
(  214)                         zt2(ji,jj-1,jk) = e1v(ji,jj-1)* fse3v(ji,jj-1,jk)   &
(  215)                            &            * zvn(ji,jj-1,jk) * ( trb(ji,jj-1,jk,jn) + trb(ji  ,jj,jk,jn) ) * 0.5
(  216)                      ENDIF
(  217)                   ENDIF
(  218) 
(  219) #else
(  220)                   IF( umask(ji,jj,jk) == 0. ) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:18      page 5

(  221)                      IF( zun(ji+1,jj,jk) > 0. .AND. ji /= jpi ) THEN
(  222)                         zt1(ji+1,jj,jk) = e2u(ji+1,jj) * zun(ji+1,jj,jk) * ( trb(ji+1,jj,jk,jn) + trb(ji+2,jj,jk,jn) ) * 0.5
(  223)                      ENDIF
(  224)                      IF( zun(ji-1,jj,jk) < 0. ) THEN
(  225)                         zt1(ji-1,jj,jk) = e2u(ji-1,jj) * zun(ji-1,jj,jk) * ( trb(ji-1,jj,jk,jn) + trb(ji  ,jj,jk,jn) ) * 0.5
(  226)                      ENDIF
(  227)                   ENDIF
(  228)                   IF( vmask(ji,jj,jk) == 0. ) THEN
(  229)                      IF( zvn(ji,jj+1,jk) > 0. .AND. jj /= jpj ) THEN
(  230)                         zt2(ji,jj+1,jk) = e1v(ji,jj+1) * zvn(ji,jj+1,jk) * ( trb(ji,jj+1,jk,jn) + trb(ji,jj+2,jk,jn) ) * 0.5
(  231)                      ENDIF
(  232)                      IF( zvn(ji,jj-1,jk) < 0. ) THEN
(  233)                         zt2(ji,jj-1,jk) = e1v(ji,jj-1) * zvn(ji,jj-1,jk) * ( trb(ji,jj-1,jk,jn) + trb(ji  ,jj,jk,jn) ) * 0.5
(  234)                      ENDIF
(  235)                   ENDIF
(  236) #endif
(  237)                END DO
(  238)             END DO
(  239)          END DO
(  240) 
(  241)          ! lateral boundary conditions on zt1, zt2   (changed sign)
(  242)          CALL lbc_lnk( zt1, 'U', -1. )
(  243)          CALL lbc_lnk( zt2, 'V', -1. )
(  244) 
(  245)          ! Compute and add the horizontal advective trend
(  246) 
(  247)          DO jk = 1, jpkm1
(  248)             DO jj = 2, jpjm1      
(  249)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  250) #if defined key_s_coord || defined key_partial_steps
(  251)                   zbtr = 1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
(  252) #else
(  253)                   zbtr = 1. / ( e1t(ji,jj)*e2t(ji,jj) )
(  254) #endif
(  255)                   ! horizontal advective trends
(  256)                   ztra = - zbtr * ( zt1(ji,jj,jk) - zt1(ji-1,jj  ,jk  )   &
(  257)                      &            + zt2(ji,jj,jk) - zt2(ji  ,jj-1,jk  ) )
(  258)                   ! add it to the general tracer trends
(  259)                   tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + ztra
(  260) #if defined key_trc_diatrd
(  261)                   ! recompute the trends in i- and j-direction as Uh gradh(T)
(  262) #   if defined key_s_coord || defined key_partial_steps
(  263)                   zfui =  e2u(ji  ,jj) * fse3u(ji,  jj,jk) * un(ji,  jj,jk)   &
(  264)                      & -  e2u(ji-1,jj) * fse3u(ji-1,jj,jk) * un(ji-1,jj,jk)
(  265)                   zfvj =  e1v(ji,jj  ) * fse3v(ji,jj  ,jk) * vn(ji,jj  ,jk)   &
(  266)                      & -  e1v(ji,jj-1) * fse3v(ji,jj-1,jk) * vn(ji,jj-1,jk)
(  267) #   else
(  268)                   zfui = e2u(ji  ,jj) * un(ji,  jj,jk)   &
(  269)                      & - e2u(ji-1,jj) * un(ji-1,jj,jk)
(  270)                   zfvj = e1v(ji,jj  ) * vn(ji,jj  ,jk)   &
(  271)                      & - e1v(ji,jj-1) * vn(ji,jj-1,jk)
(  272) #   endif
(  273)                   ztai =-zbtr * (  zt1(ji,jj,jk) - zt1(ji-1,jj  ,jk) - trn(ji,jj,jk,jn) * zfui  )
(  274)                   ztaj =-zbtr * (  zt2(ji,jj,jk) - zt2(ji  ,jj-1,jk) - trn(ji,jj,jk,jn) * zfvj  )
(  275)                   ! save i- and j- advective trends computed as Uh gradh(T)
(  276)                   IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),1) = ztai
(  277)                   IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),2) = ztaj
(  278) 






PGF90 (Version     10.5)          09/23/2010  09:49:18      page 6

(  279) #endif
(  280)                END DO
(  281)             END DO
(  282)          END DO
(  283)       ENDDO
(  284) 
(  285)       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
(  286)          WRITE(charout, FMT="('muscl2 - had')")
(  287)          CALL prt_ctl_trc_info(charout)
(  288)          CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
(  289)       ENDIF
(  290) 
(  291)       ! II. Vertical advective fluxes
(  292)       ! -----------------------------
(  293) 
(  294)       DO jn = 1, jptra
(  295) 
(  296)          ! First guess of the slope
(  297)          ! interior values
(  298)          DO jk = 2, jpkm1
(  299)             zt1(:,:,jk) = tmask(:,:,jk) * ( trb(:,:,jk-1,jn) - trb(:,:,jk,jn) )
(  300)          END DO
(  301)          ! surface and bottom boundary conditions
(  302)          zt1 (:,:, 1 ) = 0.e0 
(  303)          zt1 (:,:,jpk) = 0.e0
(  304) 
(  305)          ! Slopes
(  306)          DO jk = 2, jpkm1
(  307)             DO jj = 1, jpj
(  308)                DO ji = 1, jpi
(  309)                   ztp1(ji,jj,jk) =                    ( zt1(ji,jj,jk) + zt1(ji,jj,jk+1) )   &
(  310)                      &           * ( 0.25 + SIGN( 0.25, zt1(ji,jj,jk) * zt1(ji,jj,jk+1) ) )
(  311)                END DO
(  312)             END DO
(  313)          END DO
(  314) 
(  315)          ! Slopes limitation
(  316)          ! interior values
(  317)          DO jk = 2, jpkm1
(  318)             DO jj = 1, jpj
(  319)                DO ji = 1, jpi
(  320)                   ztp1(ji,jj,jk) = SIGN( 1., ztp1(ji,jj,jk) )   &
(  321)                      &           * MIN(    ABS( ztp1(ji,jj,jk  ) ),   &
(  322)                      &                  2.*ABS( zt1 (ji,jj,jk+1) ),   &
(  323)                      &                  2.*ABS( zt1 (ji,jj,jk  ) ) )
(  324)                END DO
(  325)             END DO
(  326)          END DO
(  327) 
(  328)          ! surface values
(  329)          ztp1(:,:,1) = 0. 
(  330) 
(  331)          ! vertical advective flux
(  332)          ! interior values
(  333)          DO jk = 1, jpkm1
(  334)             DO jj = 2, jpjm1      
(  335)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  336)                   zew = zwn(ji,jj,jk+1)






PGF90 (Version     10.5)          09/23/2010  09:49:18      page 7

(  337)                   z0w = SIGN( 0.5, zwn(ji,jj,jk+1) )
(  338)                   zalpha = 0.5 + z0w
(  339)                   zw  = z0w - 0.5 * zwn(ji,jj,jk+1)* rdttrc(jk)/ fse3w(ji,jj,jk+1)
(  340)                   zzt1 = trb(ji,jj,jk+1,jn) + zw*ztp1(ji,jj,jk+1)
(  341)                   zzt2 = trb(ji,jj,jk  ,jn) + zw*ztp1(ji,jj,jk  )
(  342)                   zt1(ji,jj,jk+1) = zew * ( zalpha * zzt1 + (1.-zalpha)*zzt2 )
(  343)                END DO
(  344)             END DO
(  345)          END DO
(  346)          DO jk = 2, jpkm1
(  347)             DO jj = 2, jpjm1
(  348)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  349)                   IF( tmask(ji,jj,jk+1) == 0. ) THEN
(  350)                      IF( zwn(ji,jj,jk) > 0. ) THEN
(  351)                         zt1(ji,jj,jk) = zwn(ji,jj,jk) * ( trb(ji,jj,jk-1,jn) + trb(ji,jj,jk,jn) ) * 0.5
(  352)                      ENDIF
(  353)                   ENDIF
(  354)                END DO
(  355)             END DO
(  356)          END DO
(  357) 
(  358)          ! surface values
(  359)          IF( lk_dynspg_rl ) THEN        ! rigid lid : flux set to zero
(  360)             zt1(:,:, 1 ) = 0.e0
(  361)          ELSE                           ! free surface
(  362)             zt1(:,:, 1 ) = zwn(:,:,1) * trb(:,:,1,jn)
(  363)          ENDIF
(  364) 
(  365)          ! bottom values
(  366)          zt1(:,:,jpk) = 0.e0
(  367) 
(  368)          ! Compute & add the vertical advective trend
(  369) 
(  370)          DO jk = 1, jpkm1
(  371)             DO jj = 2, jpjm1      
(  372)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  373)                   zbtr = 1. / fse3t(ji,jj,jk)
(  374)                   ! horizontal advective trends
(  375)                   ztra = - zbtr * ( zt1(ji,jj,jk) - zt1(ji,jj,jk+1) )
(  376)                   ! add it to the general tracer trends
(  377)                   tra(ji,jj,jk,jn) =  tra(ji,jj,jk,jn) + ztra
(  378) #if defined key_trc_diatrd
(  379)                   ! save the vertical advective trends computed as w gradz(T)
(  380)                   IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),3) = ztra - trn(ji,jj,jk,jn) * hdivn(ji,jj,jk)
(  381) #endif
(  382)                END DO
(  383)             END DO
(  384)          END DO
(  385) 
(  386)       END DO
(  387) 
(  388)       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
(  389)          WRITE(charout, FMT="('muscl2 - zad')")
(  390)          CALL prt_ctl_trc_info(charout)
(  391)          CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
(  392)       ENDIF
(  393) 
(  394)    END SUBROUTINE trc_adv_muscl2






PGF90 (Version     10.5)          09/23/2010  09:49:18      page 8

(  395) 
(  396) #else
(  397) 
(  398)    !!----------------------------------------------------------------------
(  399)    !!   Default option                                         Empty module
(  400)    !!----------------------------------------------------------------------
(  401) CONTAINS
(  402)    SUBROUTINE trc_adv_muscl2( kt ) 
(  403)       INTEGER, INTENT(in) :: kt
(  404) !      WRITE(*,*) 'trc_adv_muscl2: You should not have seen this print! error?', kt
(  405)    END SUBROUTINE trc_adv_muscl2
(  406) #endif
(  407) 
(  408)    !!======================================================================
(  409) END MODULE trcadv_muscl2

















































PGF90 (Version     10.5)          09/23/2010  09:49:18      page 9

(    1) # 1 "trcadv_muscl2.F90"
(    1) MODULE trcadv_muscl2
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  trcadv_muscl2  ***
(    4)    !! Ocean passive tracers:  horizontal & vertical advective trend
(    5)    !!==============================================================================
(    6) 
(    7) # 398
(  398)    !!----------------------------------------------------------------------
(  399)    !!   Default option                                         Empty module
(  400)    !!----------------------------------------------------------------------
(  401) CONTAINS
(  402)    SUBROUTINE trc_adv_muscl2( kt ) 
(  403)       INTEGER, INTENT(in) :: kt
(  404) !      WRITE(*,*) 'trc_adv_muscl2: You should not have seen this print! error?', kt
(  405)    END SUBROUTINE trc_adv_muscl2
(  406) 
(  407) # 408
(  408)    !!======================================================================
(  409) END MODULE trcadv_muscl2
