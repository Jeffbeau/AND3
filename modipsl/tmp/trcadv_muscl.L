


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

Filename: trcadv_muscl.F90

(    1) MODULE trcadv_muscl
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  trcadv_muscl  ***
(    4)    !! Ocean passive tracers:  horizontal & vertical advective trend
(    5)    !!==============================================================================
(    6) #if defined key_passivetrc
(    7)    !!----------------------------------------------------------------------
(    8)    !!   trc_adv_muscl : update the tracer trend with the horizontal
(    9)    !!                   and vertical advection trends using MUSCL scheme
(   10)    !!----------------------------------------------------------------------
(   11)    !! * Modules used
(   12)    USE oce_trc         ! ocean dynamics and active tracers variables
(   13)    USE trc             ! ocean passive tracers variables
(   14)    USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
(   15)    USE trcbbl          ! advective passive tracers in the BBL
(   16)    USE lib_mpp
(   17)    USE prtctl_trc      ! Print control for debbuging
(   18) 
(   19)    IMPLICIT NONE
(   20)    PRIVATE
(   21) 
(   22)    !! * Accessibility
(   23)    PUBLIC trc_adv_muscl  ! routine called by trcstp.F90
(   24) 
(   25)    !! * Substitutions
(   26) #  include "passivetrc_substitute.h90"
(   27)    !!----------------------------------------------------------------------
(   28)    !!   TOP 1.0 , LOCEAN-IPSL (2005) 
(   29)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trcadv_muscl.F90,v 1.11 2006/04/10 15:38:54 opalod Exp $ 
(   30)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   31)    !!----------------------------------------------------------------------
(   32) 
(   33) CONTAINS
(   34) 
(   35)    SUBROUTINE trc_adv_muscl( kt )
(   36)       !!----------------------------------------------------------------------
(   37)       !!                    ***  ROUTINE trc_adv_muscl  ***
(   38)       !!          
(   39)       !! ** Purpose :   Compute the now trend due to total advection of any pas-
(   40)       !!      sive tracer using a MUSCL scheme (Monotone Upstream-centered Scheme
(   41)       !!      for Conservation Laws) and add it to the general tracer trend.
(   42)       !!
(   43)       !! ** Method  :
(   44)       !!
(   45)       !! ** Action  : - update tra with the now advective tracer trends
(   46)       !!              - save trends in trtrd ('key_trc_diatrd')






PGF90 (Version     10.5)          09/23/2010  09:49:18      page 2

(   47)       !!
(   48)       !! References :                
(   49)       !!      Estubier, A., and M. Levy, Notes Techn. Pole de Modelisation
(   50)       !!	IPSL, Sept. 2000 (http://www.lodyc.jussieu.fr/opa)
(   51)       !!
(   52)       !! History :
(   53)       !!        !  06-00  (A.Estublier)  for passive tracers
(   54)       !!   9.0  !  03-04  (C. Ethe, G. Madec)  F90: Free form and module
(   55)       !!----------------------------------------------------------------------
(   56)       !! * modules used
(   57) #if defined key_trcbbl_adv
(   58)       USE oce_trc            , zun => ua,  &  ! use ua as workspace
(   59)          &                     zvn => va      ! use va as workspace
(   60)       REAL(wp), DIMENSION(jpi,jpj,jpk) :: zwn
(   61) #else
(   62)       USE oce_trc            , zun => un,  &  ! When no bbl, zun == un
(   63)                                zvn => vn,  &  !              zvn == vn
(   64)                                zwn => wn      !              zwn == wn
(   65) #endif
(   66) 
(   67)       !! * Arguments
(   68)       INTEGER, INTENT( in ) ::   kt         ! ocean time-step
(   69) 
(   70)       !! * Local declarations
(   71)       INTEGER ::   ji, jj, jk,jn            ! dummy loop indices
(   72)       REAL(wp), DIMENSION (jpi,jpj,jpk) ::   &
(   73)          zt1, zt2, ztp1, ztp2
(   74) 
(   75)       REAL(wp) ::   zu, zv, zw, zeu, zev, zew, zbtr, ztra
(   76)       REAL(wp) ::   z0u, z0v, z0w
(   77)       REAL(wp) ::   zzt1, zzt2, zalpha, z2dtt
(   78) #if defined key_trc_diatrd
(   79)       REAL(wp) ::   ztai, ztaj
(   80)       REAL(wp) ::   zfui, zfvj
(   81) #endif
(   82)       CHARACTER (len=22) :: charout
(   83)       !!----------------------------------------------------------------------
(   84) 
(   85) 
(   86)       IF( kt == nittrc000 .AND. lwp ) THEN
(   87)          WRITE(numout,*)
(   88)          WRITE(numout,*) 'trc_adv : MUSCL advection scheme'
(   89)          WRITE(numout,*) '~~~~~~~'
(   90)       ENDIF
(   91) 
(   92)  
(   93) 
(   94) #if defined key_trcbbl_adv
(   95)       ! Advective bottom boundary layer
(   96)       ! -------------------------------
(   97)       zun(:,:,:) = un (:,:,:) - u_trc_bbl(:,:,:)
(   98)       zvn(:,:,:) = vn (:,:,:) - v_trc_bbl(:,:,:)
(   99)       zwn(:,:,:) = wn (:,:,:) + w_trc_bbl(:,:,:)
(  100) #endif
(  101) 
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
(  123)          ! lateral boundary conditions on zt1, zt2
(  124)          CALL lbc_lnk( zt1, 'U', -1. )   
(  125)          CALL lbc_lnk( zt2, 'V', -1. ) 
(  126) 
(  127) 
(  128)          ! Slopes
(  129)          ! interior values
(  130)          DO jk = 1, jpkm1
(  131)             DO jj = 2, jpj
(  132)                DO ji = fs_2, jpi   ! vector opt.
(  133)                   ztp1(ji,jj,jk) =                    ( zt1(ji,jj,jk) + zt1(ji-1,jj  ,jk) )   &
(  134)                      &           * ( 0.25 + SIGN( 0.25, zt1(ji,jj,jk) * zt1(ji-1,jj  ,jk) ) )
(  135)                   ztp2(ji,jj,jk) =                    ( zt2(ji,jj,jk) + zt2(ji  ,jj-1,jk) )   &
(  136)                      &           * ( 0.25 + SIGN( 0.25, zt2(ji,jj,jk) * zt2(ji  ,jj-1,jk) ) )
(  137)                END DO
(  138)             END DO
(  139)          END DO
(  140)          ! bottom values
(  141)          ztp1(:,:,jpk) = 0.e0
(  142)          ztp2(:,:,jpk) = 0.e0
(  143) 
(  144)          ! Slopes limitation
(  145)          DO jk = 1, jpkm1
(  146)             DO jj = 2, jpj
(  147)                DO ji = fs_2, jpi   ! vector opt.
(  148)                   ztp1(ji,jj,jk) = SIGN( 1., ztp1(ji,jj,jk) )   &
(  149)                      &           * MIN(    ABS( ztp1(ji  ,jj,jk) ),   &
(  150)                      &                  2.*ABS( zt1 (ji-1,jj,jk) ),   &
(  151)                      &                  2.*ABS( zt1 (ji  ,jj,jk) ) )
(  152) 
(  153)                   ztp2(ji,jj,jk) = SIGN( 1., ztp2(ji,jj,jk) )   &
(  154)                      &           * MIN(    ABS( ztp2(ji,jj  ,jk) ),   &
(  155)                      &                  2.*ABS( zt2 (ji,jj-1,jk) ),   &
(  156)                      &                  2.*ABS( zt2 (ji,jj  ,jk) ) )
(  157) 
(  158)                END DO
(  159)             END DO
(  160)          END DO
(  161) 
(  162)          ! Advection terms






PGF90 (Version     10.5)          09/23/2010  09:49:18      page 4

(  163)          ! interior values
(  164)          DO jk = 1, jpkm1
(  165)             DO jj = 2, jpjm1      
(  166)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  167)                   ! volume fluxes
(  168) #if defined key_s_coord || defined key_partial_steps
(  169)                   zeu = e2u(ji,jj) * fse3u(ji,jj,jk) * zun(ji,jj,jk)
(  170)                   zev = e1v(ji,jj) * fse3v(ji,jj,jk) * zvn(ji,jj,jk)
(  171) #else
(  172)                   zeu = e2u(ji,jj) * zun(ji,jj,jk)
(  173)                   zev = e1v(ji,jj) * zvn(ji,jj,jk)
(  174) #endif
(  175)                   ! MUSCL fluxes
(  176)                   z2dtt = rdttra(jk) * FLOAT(ndttrc)
(  177)                   z0u = SIGN( 0.5, zun(ji,jj,jk) )            
(  178)                   zalpha = 0.5 - z0u
(  179)                   zu  = z0u - 0.5 * zun(ji,jj,jk) * z2dtt / e1u(ji,jj)
(  180)                   zzt1 = trb(ji+1,jj,jk,jn) + zu*ztp1(ji+1,jj,jk)
(  181)                   zzt2 = trb(ji  ,jj,jk,jn) + zu*ztp1(ji  ,jj,jk)
(  182)                   zt1(ji,jj,jk) = zeu * ( zalpha * zzt1 + (1.-zalpha) * zzt2 )
(  183)                   z0v = SIGN( 0.5, zvn(ji,jj,jk) )            
(  184)                   zalpha = 0.5 - z0v
(  185)                   zv  = z0v - 0.5 * zvn(ji,jj,jk) * z2dtt / e2v(ji,jj)
(  186)                   zzt1 = trb(ji,jj+1,jk,jn) + zv*ztp2(ji,jj+1,jk)
(  187)                   zzt2 = trb(ji,jj  ,jk,jn) + zv*ztp2(ji,jj  ,jk)
(  188)                   zt2(ji,jj,jk) = zev * ( zalpha * zzt1 + (1.-zalpha) * zzt2 )
(  189)                END DO
(  190)             END DO
(  191)          END DO
(  192) 
(  193)          ! lateral boundary conditions on zt1, zt2 (changed sign)
(  194)          CALL lbc_lnk( zt1, 'U', -1. ) 
(  195)          CALL lbc_lnk( zt2, 'V', -1. ) 
(  196) 
(  197)          ! Compute and add the horizontal advective trend
(  198) 
(  199)          DO jk = 1, jpkm1
(  200)             DO jj = 2, jpjm1      
(  201)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  202) #if defined key_s_coord || defined key_partial_steps
(  203)                   zbtr = 1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
(  204) #else
(  205)                   zbtr = 1. / ( e1t(ji,jj)*e2t(ji,jj) )
(  206) #endif
(  207)                   ! horizontal advective trends
(  208)                   ztra = - zbtr * ( zt1(ji,jj,jk) - zt1(ji-1,jj  ,jk  )   &
(  209)                      &            + zt2(ji,jj,jk) - zt2(ji  ,jj-1,jk  ) )
(  210)                   ! add it to the general tracer trends
(  211)                   tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + ztra
(  212) #if defined key_trc_diatrd
(  213)                   ! recompute the trends in i- and j-direction as Uh gradh(T)
(  214) #   if defined key_s_coord || defined key_partial_steps
(  215)                   zfui =  e2u(ji  ,jj) * fse3u(ji,  jj,jk) * un(ji,  jj,jk)   &
(  216)                      & -  e2u(ji-1,jj) * fse3u(ji-1,jj,jk) * un(ji-1,jj,jk)
(  217)                   zfvj =  e1v(ji,jj  ) * fse3v(ji,jj  ,jk) * vn(ji,jj  ,jk)   &
(  218)                      & -  e1v(ji,jj-1) * fse3v(ji,jj-1,jk) * vn(ji,jj-1,jk)
(  219) #   else
(  220)                   zfui = e2u(ji  ,jj) * un(ji,  jj,jk)   &






PGF90 (Version     10.5)          09/23/2010  09:49:18      page 5

(  221)                      & - e2u(ji-1,jj) * un(ji-1,jj,jk)
(  222)                   zfvj = e1v(ji,jj  ) * vn(ji,jj  ,jk)   &
(  223)                      & - e1v(ji,jj-1) * vn(ji,jj-1,jk)
(  224) #   endif
(  225)                   ztai =-zbtr * (  zt1(ji,jj,jk) - zt1(ji-1,jj  ,jk) - trn(ji,jj,jk,jn) * zfui  )
(  226)                   ztaj =-zbtr * (  zt2(ji,jj,jk) - zt2(ji  ,jj-1,jk) - trn(ji,jj,jk,jn) * zfvj  )
(  227)                   ! save i- and j- advective trends computed as Uh gradh(T)
(  228)                   IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),1) = ztai
(  229)                   IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),2) = ztaj
(  230) #endif
(  231)                END DO
(  232)             END DO
(  233)          END DO
(  234)       ENDDO
(  235) 
(  236)       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
(  237)          WRITE(charout, FMT="('muscl - had')")
(  238)          CALL prt_ctl_trc_info(charout)
(  239)          CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
(  240)       ENDIF
(  241) 
(  242)          ! II. Vertical advective fluxes
(  243)          ! -----------------------------
(  244) 
(  245)       DO jn = 1, jptra
(  246)          ! First guess of the slope
(  247)          ! interior values
(  248)          DO jk = 2, jpkm1
(  249)             zt1(:,:,jk) = tmask(:,:,jk) * ( trb(:,:,jk-1,jn) - trb(:,:,jk,jn) )
(  250)          END DO
(  251)          ! surface and bottom boundary conditions
(  252)          zt1 (:,:, 1 ) = 0.e0 
(  253)          zt1 (:,:,jpk) = 0.e0
(  254)          ! Slopes
(  255)          DO jk = 2, jpkm1
(  256)             DO jj = 1, jpj
(  257)                DO ji = 1, jpi
(  258)                   ztp1(ji,jj,jk) =                    ( zt1(ji,jj,jk) + zt1(ji,jj,jk+1) )   &
(  259)                      &           * ( 0.25 + SIGN( 0.25, zt1(ji,jj,jk) * zt1(ji,jj,jk+1) ) )
(  260)                END DO
(  261)             END DO
(  262)          END DO
(  263) 
(  264)          ! Slopes limitation
(  265)          ! interior values
(  266)          DO jk = 2, jpkm1
(  267)             DO jj = 1, jpj
(  268)                DO ji = 1, jpi
(  269)                   ztp1(ji,jj,jk) = SIGN( 1., ztp1(ji,jj,jk) )   &
(  270)                      &           * MIN(    ABS( ztp1(ji,jj,jk  ) ),   &
(  271)                      &                  2.*ABS( zt1 (ji,jj,jk+1) ),   &
(  272)                      &                  2.*ABS( zt1 (ji,jj,jk  ) ) )
(  273)                END DO
(  274)             END DO
(  275)          END DO
(  276)          ! surface values
(  277)          ztp1(:,:,1) = 0. 
(  278)          ! vertical advective flux






PGF90 (Version     10.5)          09/23/2010  09:49:18      page 6

(  279)          ! interior values
(  280)          DO jk = 1, jpkm1
(  281)             DO jj = 2, jpjm1      
(  282)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  283)                   z2dtt = rdttra(jk) * FLOAT(ndttrc)
(  284)                   zew = zwn(ji,jj,jk+1)
(  285)                   z0w = SIGN( 0.5, zwn(ji,jj,jk+1) )
(  286)                   zalpha = 0.5 + z0w
(  287)                   zw  = z0w - 0.5 * zwn(ji,jj,jk+1)*z2dtt / fse3w(ji,jj,jk+1)
(  288)                   zzt1 = trb(ji,jj,jk+1,jn) + zw*ztp1(ji,jj,jk+1)
(  289)                   zzt2 = trb(ji,jj,jk  ,jn) + zw*ztp1(ji,jj,jk  )
(  290)                   zt1(ji,jj,jk+1) = zew * ( zalpha * zzt1 + (1.-zalpha)*zzt2 )
(  291)                END DO
(  292)             END DO
(  293)          END DO
(  294)          ! surface values
(  295)          IF( lk_dynspg_rl ) THEN        ! rigid lid : flux set to zero
(  296)             zt1(:,:, 1 ) = 0.e0
(  297)          ELSE                           ! free surface
(  298)             zt1(:,:, 1 ) = zwn(:,:,1) * trb(:,:,1,jn)
(  299)          ENDIF
(  300) 
(  301)          ! bottom values
(  302)          zt1(:,:,jpk) = 0.e0
(  303) 
(  304)          ! Compute & add the vertical advective trend
(  305) 
(  306)          DO jk = 1, jpkm1
(  307)             DO jj = 2, jpjm1      
(  308)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  309)                   zbtr = 1. / fse3t(ji,jj,jk)
(  310)                   ! horizontal advective trends
(  311)                   ztra = - zbtr * ( zt1(ji,jj,jk) - zt1(ji,jj,jk+1) )
(  312)                   ! add it to the general tracer trends
(  313)                   tra(ji,jj,jk,jn) =  tra(ji,jj,jk,jn) + ztra
(  314) #if defined key_trc_diatrd
(  315)                   ! save the vertical advective trends computed as w gradz(T)
(  316)                   IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),3) = ztra - trn(ji,jj,jk,jn) * hdivn(ji,jj,jk)
(  317) #endif
(  318)                END DO
(  319)             END DO
(  320)          END DO
(  321) 
(  322)       END DO
(  323) 
(  324)       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
(  325)          WRITE(charout, FMT="('muscl - zad')")
(  326)          CALL prt_ctl_trc_info(charout)
(  327)          CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
(  328)       ENDIF
(  329) 
(  330) END SUBROUTINE trc_adv_muscl
(  331) 
(  332) #else
(  333)    !!----------------------------------------------------------------------
(  334)    !!   Default option                                         Empty module
(  335)    !!----------------------------------------------------------------------
(  336) CONTAINS






PGF90 (Version     10.5)          09/23/2010  09:49:18      page 7

(  337)    SUBROUTINE trc_adv_muscl( kt )  
(  338)       INTEGER, INTENT(in) :: kt
(  339) !      WRITE(*,*) 'trc_adv_muscl: You should not have seen this print! error?', kt
(  340)    END SUBROUTINE trc_adv_muscl
(  341) #endif
(  342) 
(  343)    !!======================================================================
(  344) END MODULE trcadv_muscl
























































PGF90 (Version     10.5)          09/23/2010  09:49:18      page 8

(    1) # 1 "trcadv_muscl.F90"
(    1) MODULE trcadv_muscl
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  trcadv_muscl  ***
(    4)    !! Ocean passive tracers:  horizontal & vertical advective trend
(    5)    !!==============================================================================
(    6) # 333
(  333)    !!----------------------------------------------------------------------
(  334)    !!   Default option                                         Empty module
(  335)    !!----------------------------------------------------------------------
(  336) CONTAINS
(  337)    SUBROUTINE trc_adv_muscl( kt )  
(  338)       INTEGER, INTENT(in) :: kt
(  339) !      WRITE(*,*) 'trc_adv_muscl: You should not have seen this print! error?', kt
(  340)    END SUBROUTINE trc_adv_muscl
(  341) 
(  342) # 343
(  343)    !!======================================================================
(  344) END MODULE trcadv_muscl
