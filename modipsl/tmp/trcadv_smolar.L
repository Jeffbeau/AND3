


PGF90 (Version     10.5)          09/23/2010  09:49:19      page 1

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

Filename: trcadv_smolar.F90

(    1) MODULE trcadv_smolar
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  trcadv_smolar  ***
(    4)    !! Ocean passive tracers:  horizontal & vertical advective trend
(    5)    !!==============================================================================
(    6) #if defined key_passivetrc
(    7)    !!----------------------------------------------------------------------
(    8)    !!   trc_adv_smolar : update the passive tracer trend with the horizontal
(    9)    !!                  and vertical advection trends using a Smolarkiewicz 
(   10)    !!                  FCT scheme
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE oce_trc             ! ocean dynamics and active tracers variables
(   14)    USE trc                 ! ocean passive tracers variables
(   15)    USE lbclnk              ! ocean lateral boundary conditions (or mpp link)
(   16)    USE trcbbl              ! advective passive tracers in the BBL
(   17)    USE prtctl_trc      ! Print control for debbuging
(   18) 
(   19)    IMPLICIT NONE
(   20)    PRIVATE
(   21) 
(   22)    !! * Accessibility
(   23)    PUBLIC trc_adv_smolar    ! routine called by trcstp.F90
(   24) 
(   25)    !! * Module variable
(   26)    REAL(wp), DIMENSION(jpk) ::   &
(   27)       rdttrc                     ! vertical profile of tracer time-step
(   28)  
(   29)    !! * Substitutions
(   30) #  include "passivetrc_substitute.h90"
(   31)    !!----------------------------------------------------------------------
(   32)    !!   TOP 1.0 , LOCEAN-IPSL (2005) 
(   33)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trcadv_smolar.F90,v 1.11 2006/04/10 15:38:54 opalod Exp $ 
(   34)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   35)    !!----------------------------------------------------------------------
(   36) CONTAINS
(   37) 
(   38)    SUBROUTINE trc_adv_smolar( kt )
(   39)       !!----------------------------------------------------------------------
(   40)       !!                   ***  ROUTINE trc_adv_smolar  ***
(   41)       !!
(   42)       !! ** Purpose :   Compute the now trend due to total advection of passi-
(   43)       !!      ve tracer using a Smolarkiewicz FCT (Flux Corrected Transport ) 
(   44)       !!      scheme and add it to the general tracer trend.
(   45)       !!
(   46)       !! ** Method : Computation of not exactly the advection but the






PGF90 (Version     10.5)          09/23/2010  09:49:19      page 2

(   47)       !!             transport term, i.e.  div(u*tra). 
(   48)       !!             Computes the now horizontal and vertical advection with
(   49)       !!             the complete 3d method.
(   50)       !! 
(   51)       !!       note: - sc is an empirical factor to be used with care
(   52)       !!             - this advection scheme needs an euler-forward time scheme
(   53)       !!
(   54)       !! ** Action : - update tra with the now advective tracer trends
(   55)       !!             - save trends in trtrd ('key_trc_diatrd')
(   56)       !!
(   57)       !! References :                
(   58)       !!     Piotr K. Smolarkiewicz, 1983,
(   59)       !!       "A simple positive definit advection
(   60)       !!        scheme with small IMPLICIT diffusion"
(   61)       !!        Monthly Weather Review, pp 479-486
(   62)       !!
(   63)       !! History :
(   64)       !!         !  87-06 (pa-dl) Original
(   65)       !!         !  91-11 (G. Madec)
(   66)       !!         !  94-08 (A. Czaja)
(   67)       !!         !  95-09 (M. Levy) passive tracers
(   68)       !!         !  98-03 (M.A. Foujols) lateral boundary conditions
(   69)       !!         !  99-02 (M.A. Foujols) lbc in conjonction with ORCA
(   70)       !!         !  00-05 (MA Foujols) add lbc for tracer trends
(   71)       !!         !  00-10 (MA Foujols and E.Kestenare) INCLUDE instead of routine
(   72)       !!         !  01-05 (E.Kestenare) fix bug in trtrd indexes
(   73)       !!         !  02-05 (M-A Filiberti, and M.Levy) correction in trtrd computation
(   74)       !!   9.0   !  03-04  (C. Ethe)  F90: Free form and module
(   75)       !!----------------------------------------------------------------------
(   76)       !! * modules used
(   77) #if defined key_trcbbl_adv
(   78)       USE oce_trc            , zun => ua,  &  ! use ua as workspace
(   79)          &                     zvn => va      ! use va as workspace
(   80)       REAL(wp), DIMENSION(jpi,jpj,jpk) :: zwn
(   81) #else
(   82)       USE oce_trc            , zun => un,  &  ! When no bbl, zun == un
(   83)                                zvn => vn,  &  !              zvn == vn
(   84)                                zwn => wn      !              zwn == wn
(   85) #endif
(   86)       !! * Arguments
(   87)       INTEGER, INTENT( in ) ::   kt         ! ocean time-step
(   88) 
(   89)       !! * Local declarations
(   90)       INTEGER :: ji, jj, jk,jt, jn            ! dummy loop indices
(   91) 
(   92)       REAL(wp), DIMENSION (jpi,jpj,jpk) ::   &      
(   93)          zti, ztj,              &
(   94)          zaa, zbb, zcc,         &
(   95)          zx , zy , zz ,         &
(   96)          zkx, zky, zkz,         &
(   97)          zbuf
(   98) 
(   99) #if defined key_trc_diatrd
(  100)       REAL(wp) :: zgm, zgz
(  101) #endif
(  102) 
(  103)       REAL(wp) :: zbtr, ztra
(  104)       REAL(wp) :: zfp_ui, zfp_vj, zfm_ui, zfm_vj, zfp_w, zfm_w






PGF90 (Version     10.5)          09/23/2010  09:49:19      page 3

(  105)       CHARACTER (len=22) :: charout
(  106)       !!----------------------------------------------------------------------
(  107) 
(  108) 
(  109)       IF( kt == nittrc000  .AND. lwp ) THEN
(  110)          WRITE(numout,*)
(  111)          WRITE(numout,*) 'trc_adv_smolar : SMOLARKIEWICZ advection scheme'
(  112)          WRITE(numout,*) '~~~~~~~~~~~~~~~'
(  113)          rdttrc(:) = rdttra(:) * FLOAT(ndttrc)
(  114)       ENDIF
(  115) 
(  116) 
(  117) #if defined key_trcbbl_adv        
(  118)       ! Advective bottom boundary layer
(  119)       ! -------------------------------
(  120)       zun(:,:,:) = un (:,:,:) - u_trc_bbl(:,:,:)
(  121)       zvn(:,:,:) = vn (:,:,:) - v_trc_bbl(:,:,:)
(  122)       zwn(:,:,:) = wn (:,:,:) + w_trc_bbl( :,:,:)
(  123) #endif
(  124) 
(  125)       ! tracer loop parallelized (macrotasking)
(  126)       ! =======================================
(  127)       
(  128)       DO jn = 1, jptra
(  129)          
(  130)          ! 1. tracer flux in the 3 directions
(  131)          ! ----------------------------------
(  132)          
(  133)          ! 1.1 mass flux at u v and t-points and initialization
(  134) 
(  135)         DO jk = 1,jpk
(  136) 
(  137)            DO jj = 1,jpj
(  138)               DO ji = 1,jpi
(  139)                  zaa(ji,jj,jk) = e2u(ji,jj)*fse3u(ji,jj,jk) * zun(ji,jj,jk)
(  140)                  zbb(ji,jj,jk) = e1v(ji,jj)*fse3v(ji,jj,jk) * zvn(ji,jj,jk)
(  141)                  zcc(ji,jj,jk) = e1t(ji,jj)*e2t(ji,jj)      * zwn(ji,jj,jk)
(  142)                  zbuf(ji,jj,jk) = 0.
(  143)                  ztj(ji,jj,jk) = 0.
(  144)                  zx(ji,jj,jk) = 0.
(  145)                  zy(ji,jj,jk) = 0.
(  146)                  zz(ji,jj,jk) = 0.
(  147)                  zti(ji,jj,jk) = trn(ji,jj,jk,jn)
(  148) #if defined key_trc_diatrd
(  149)                  IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),1) = 0.
(  150)                  IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),2) = 0.
(  151)                  IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),3) = 0.
(  152) #endif
(  153)               END DO
(  154)            END DO
(  155)            
(  156)            ! 1.2 calcul of intermediate field with an upstream advection scheme
(  157)            !     and mass fluxes calculated above
(  158)            
(  159)            ! calcul of tracer flux in the i and j direction
(  160)            
(  161)            DO jj=1,jpj
(  162)               zkx(  1,jj,jk)=0.






PGF90 (Version     10.5)          09/23/2010  09:49:19      page 4

(  163)               zkx(jpi,jj,jk)=0.
(  164)            END DO
(  165)            
(  166)            DO ji=1,jpi
(  167)               zky(ji,  1,jk)=0.
(  168)               zky(ji,jpj,jk)=0.
(  169)            END DO
(  170)            
(  171)            DO jj = 2,jpjm1
(  172)               DO ji = 2,jpim1
(  173)                  zfp_ui = 0.5 * ( zaa(ji,jj,jk) + ABS( zaa(ji,jj,jk) ) )
(  174)                  zfp_vj = 0.5 * ( zbb(ji,jj,jk) + ABS( zbb(ji,jj,jk) ) )
(  175)                  zfm_ui = 0.5 * ( zaa(ji,jj,jk) - ABS( zaa(ji,jj,jk) ) )
(  176)                  zfm_vj = 0.5 * ( zbb(ji,jj,jk) - ABS( zbb(ji,jj,jk) ) )            
(  177)                  zkx(ji,jj,jk) = zfp_ui * zti(ji,jj,jk) + zfm_ui * zti(ji+1,jj  ,jk) 
(  178)                  zky(ji,jj,jk) = zfp_vj * zti(ji,jj,jk) + zfm_vj * zti(ji    ,jj+1,jk)             
(  179)               END DO
(  180)            END DO
(  181) 
(  182)         END DO
(  183) 
(  184)          ! II. Vertical advection
(  185)          ! ----------------------
(  186) 
(  187)          ! Surface value
(  188)          IF( lk_dynspg_rl ) THEN        ! rigid lid : flux set to zero
(  189)             zkz(:,:, 1 ) = 0.e0  
(  190)          ELSE                           ! free surface
(  191)             zkz(:,:, 1 ) = zwn(:,:,1) * trn(:,:,1,jn) * tmask(ji,jj,1)
(  192)          ENDIF
(  193) 
(  194)         DO jk = 2,jpk
(  195)           DO jj = 1,jpj
(  196)             DO ji = 1,jpi
(  197)                zfp_w = 0.5 * ( zcc(ji,jj,jk) + ABS( zcc(ji,jj,jk) ) )
(  198)                zfm_w = 0.5 * ( zcc(ji,jj,jk) - ABS( zcc(ji,jj,jk) ) )        
(  199)                zkz(ji,jj,jk) = zfp_w * zti(ji,jj,jk) + zfm_w * zti(ji,jj,jk-1)     
(  200)             END DO
(  201)           END DO
(  202)         END DO
(  203) 
(  204) ! ... Lateral boundary conditions on zk[xy]
(  205)       CALL lbc_lnk( zkx, 'U', -1. )
(  206)       CALL lbc_lnk( zky, 'V', -1. )
(  207) 
(  208) 
(  209) ! 2. calcul of after field using an upstream advection scheme
(  210) ! -----------------------------------------------------------
(  211) 
(  212)         DO jk = 1,jpkm1
(  213)           DO jj = 2,jpjm1
(  214)             DO ji = 2,jpim1
(  215)               zbtr = 1./(e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk))
(  216)               ztj(ji,jj,jk) = -zbtr*    &
(  217)      &            ( zkx(ji,jj,jk) - zkx(ji - 1,jj,jk)  &
(  218)      &            + zky(ji,jj,jk) - zky(ji,jj - 1,jk)  &
(  219)      &            + zkz(ji,jj,jk) - zkz(ji,jj,jk + 1) )
(  220) #if defined key_trc_diatrd






PGF90 (Version     10.5)          09/23/2010  09:49:19      page 5

(  221)               IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),1) = trtrd(ji,jj,jk,ikeep(jn),1) -  &
(  222)      &                       zbtr*( zkx(ji,jj,jk) - zkx(ji - 1,jj,jk) )
(  223) 
(  224)               IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),2) = trtrd(ji,jj,jk,ikeep(jn),2) -  &
(  225)      &            zbtr*( zky(ji,jj,jk) - zky(ji,jj - 1,jk) )
(  226) 
(  227)               IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),3) = trtrd(ji,jj,jk,ikeep(jn),3) -  &
(  228)      &            zbtr*( zkz(ji,jj,jk) - zkz(ji,jj,jk + 1) )
(  229) #endif
(  230)             END DO
(  231)           END DO
(  232)         END DO
(  233) 
(  234) ! 2.1 start of antidiffusive correction loop
(  235) 
(  236)         DO jt = 1,ncortrc
(  237) 
(  238) ! 2.2 calcul of intermediary field zti
(  239) 
(  240)           DO jk = 1,jpkm1
(  241)             DO jj = 2,jpjm1
(  242)               DO ji = 2,jpim1
(  243)                 zti(ji,jj,jk) = zti(ji,jj,jk)+rdttrc(jk)*ztj(ji,jj,jk)
(  244)                 zbuf(ji,jj,jk) = zbuf(ji,jj,jk) + ztj(ji,jj,jk)
(  245)               END DO
(  246)             END DO
(  247)           END DO
(  248) 
(  249) ! ... Lateral boundary conditions on zti
(  250)       CALL lbc_lnk( zti, 'T', 1. )
(  251) 
(  252) 
(  253) ! 2.3 calcul of the antidiffusive flux
(  254) 
(  255)           DO jk = 1,jpkm1
(  256)             DO jj = 2,jpjm1
(  257)               DO ji = 2,jpim1
(  258)                 zx(ji,jj,jk) = ( abs(zaa(ji,jj,jk)) - rdttrc(jk)       &
(  259)      &              *zaa(ji,jj,jk)**2/                          &
(  260)      &              (e1u(ji,jj)*e2u(ji,jj)*fse3u(ji,jj,jk) ) )  &
(  261)      &              *(zti(ji + 1,jj,jk) - zti( ji ,jj,jk))      &
(  262)      &              /(zti( ji ,jj,jk) + zti(ji + 1,jj,jk) + rtrn)    &
(  263)      &              * rsc
(  264) 
(  265)                 zy(ji,jj,jk) = ( abs(zbb(ji,jj,jk)) - rdttrc(jk)       &
(  266)      &              *zbb(ji,jj,jk)**2/                          &
(  267)      &              (e1v(ji,jj)*e2v(ji,jj)*fse3v(ji,jj,jk) ) )  &
(  268)      &              *(zti(ji,jj + 1,jk) - zti(ji, jj ,jk))      &
(  269)      &              /(zti(ji, jj ,jk) + zti(ji,jj + 1,jk) + rtrn)    &
(  270)      &              * rsc
(  271)               END DO
(  272)             END DO
(  273)           END DO
(  274) 
(  275)           DO jk = 2,jpkm1
(  276)             DO jj = 2,jpjm1
(  277)               DO ji = 2,jpim1
(  278)                 zz(ji,jj,jk) = ( abs(zcc(ji,jj,jk)) - rdttrc(jk)*zcc(ji,jj,jk)**2  &  






PGF90 (Version     10.5)          09/23/2010  09:49:19      page 6

(  279)      &              /( e1t(ji,jj)*e2t(ji,jj)*fse3w(ji,jj,jk) ) ) &
(  280)      &               *( zti(ji,jj,jk) - zti(ji,jj,jk - 1) )/  &
(  281)      &                ( zti(ji,jj,jk) + zti(ji,jj,jk - 1) + rtrn )* rsc*( -1.)
(  282)               END DO
(  283)             END DO
(  284)           END DO
(  285) 
(  286) ! 2.4 cross terms
(  287) 
(  288)           IF (crosster) THEN
(  289)               DO jk = 2,jpkm1
(  290)                 DO jj = 2,jpjm1
(  291)                   DO ji = 2,jpim1
(  292)                     zx(ji,jj,jk) = zx(ji,jj,jk) &
(  293)      &                  - 0.5*rdttrc(jk)*rsc*zaa(ji,jj,jk)*0.25* &
(  294)      &                  (    (zbb(ji  ,jj - 1,jk  ) + zbb(ji + 1,jj - 1 &
(  295)      &                  ,jk  ) + zbb(ji + 1,jj  ,jk  ) + zbb(ji  ,jj  &
(  296)      &                  ,jk))* (zti(ji  ,jj + 1,jk  ) + zti(ji + 1,jj + &
(  297)      &                  1,jk  ) - zti(ji + 1,jj - 1,jk  ) - zti(ji  ,jj &
(  298)      &                  - 1,jk  ))/ (zti(ji  ,jj + 1,jk  ) + zti(ji + 1 &
(  299)      &                  ,jj + 1,jk  ) + zti(ji + 1,jj - 1,jk  ) + zti(ji &
(  300)      &                  ,jj - 1,jk  ) + rtrn) + (zcc(ji  ,jj  ,jk  ) + &
(  301)      &                  zcc(ji + 1,jj  ,jk  ) + zcc(ji  ,jj  ,jk + 1) + &
(  302)      &                  zcc(ji + 1,jj  ,jk + 1))* (zti(ji  ,jj  ,jk - 1) &
(  303)      &                  + zti(ji + 1,jj  ,jk - 1) - zti(ji  ,jj  ,jk + 1 &
(  304)      &                  )- zti(ji + 1,jj  ,jk + 1))/ (zti(ji  ,jj  ,jk - &
(  305)      &                  1) + zti(ji + 1,jj  ,jk - 1) + zti(ji  ,jj  ,jk &
(  306)      &                  +1) + zti(ji + 1,jj  ,jk + 1) + rtrn))/(e1u(ji &
(  307)      &                  ,jj)*e2u(ji,jj)*fse3u(ji,jj,jk))*vmask(ji  ,jj - &
(  308)      &                  1,jk  )*vmask(ji + 1,jj - 1,jk  )*vmask(ji + 1 &
(  309)      &                  ,jj,jk)*vmask(ji  ,jj  ,jk  )*tmask(ji  ,jj  ,jk &
(  310)      &                  )*tmask(ji + 1,jj  ,jk  )*tmask(ji  ,jj  ,jk + 1 &
(  311)      &                  )*tmask(ji + 1,jj  ,jk + 1)
(  312) 
(  313)                     zy(ji,jj,jk) = zy(ji,jj,jk)    &   
(  314)      &                  - 0.5*rdttrc(jk)*rsc*zbb(ji,jj,jk)*0.25*    &   
(  315)      &                  (    (zaa(ji - 1,jj  ,jk  ) + zaa(ji - 1,jj + 1    &   
(  316)      &                  ,jk  ) + zaa(ji  ,jj  ,jk  ) + zaa(ji  ,jj + 1    &   
(  317)      &                  ,jk))* (zti(ji + 1,jj + 1,jk  ) + zti(ji + 1,jj    &   
(  318)      &                  ,jk  ) - zti(ji - 1,jj + 1,jk  ) - zti(ji - 1,jj    &   
(  319)      &                  ,jk  ))/ (zti(ji + 1,jj + 1,jk  ) + zti(ji + 1    &   
(  320)      &                  ,jj  ,jk  ) + zti(ji - 1,jj + 1,jk  ) + zti(ji    &   
(  321)      &                  - 1,jj  ,jk  ) + rtrn) + (zcc(ji  ,jj  ,jk  )    &   
(  322)      &                  + zcc(ji  ,jj  ,jk + 1) + zcc(ji  ,jj + 1,jk  )    &   
(  323)      &                  + zcc(ji  ,jj + 1,jk + 1))* (zti(ji  ,jj  ,jk -    &   
(  324)      &                  1) + zti(ji  ,jj + 1,jk - 1) - zti(ji  ,jj  ,jk    &   
(  325)      &                  +1) - zti(ji  ,jj + 1,jk + 1))/ (zti(ji  ,jj    &   
(  326)      &                  ,jk- 1) + zti(ji  ,jj + 1,jk - 1) + zti(ji  ,jj    &   
(  327)      &                  ,jk+ 1) + zti(ji  ,jj + 1,jk + 1) + rtrn))    &   
(  328)      &                  /(e1v(ji,jj)*e2v(ji,jj)*fse3v(ji,jj,jk))    &   
(  329)      &                  *umask(ji - 1,jj,jk  )*umask(ji - 1,jj + 1,jk  )    &   
(  330)      &                  *umask(ji  ,jj,jk  )*umask(ji  ,jj + 1,jk  )    &   
(  331)      &                  *tmask(ji  ,jj,jk)*tmask(ji  ,jj  ,jk + 1)    &   
(  332)      &                  *tmask(ji  ,jj + 1,jk)*tmask(ji  ,jj + 1,jk + 1)      
(  333) 
(  334)                     zz(ji,jj,jk) = zz(ji,jj,jk)    &   
(  335)      &                  - 0.5*rdttrc(jk)*rsc*zcc(ji,jj,jk)*0.25*    &   
(  336)      &                  (    (zaa(ji - 1,jj  ,jk  ) + zaa(ji  ,jj  ,jk    &   






PGF90 (Version     10.5)          09/23/2010  09:49:19      page 7

(  337)      &                  ) + zaa(ji  ,jj  ,jk - 1) + zaa(ji - 1,jj  ,jk -    &   
(  338)      &                  1))*(zti(ji + 1,jj  ,jk - 1) + zti(ji + 1,jj    &   
(  339)      &                  ,jk  ) - zti(ji - 1,jj  ,jk  ) - zti(ji - 1,jj    &   
(  340)      &                  ,jk - 1))/(zti(ji + 1,jj  ,jk - 1) + zti(ji + 1    &   
(  341)      &                  ,jj,jk  ) + zti(ji - 1,jj  ,jk  ) + zti(ji - 1    &   
(  342)      &                  ,jj,jk - 1) + rtrn) + (zbb(ji  ,jj - 1,jk  )    &   
(  343)      &                  + zbb(ji  ,jj  ,jk  ) + zbb(ji  ,jj  ,jk - 1)    &   
(  344)      &                  + zbb(ji  ,jj - 1,jk - 1))*(zti(ji  ,jj + 1,jk -    &   
(  345)      &                  1) + zti(ji  ,jj + 1,jk  ) - zti(ji  ,jj - 1,jk    &   
(  346)      &                  ) - zti(ji  ,jj - 1,jk - 1))/(zti(ji  ,jj + 1,jk    &   
(  347)      &                  - 1) + zti(ji  ,jj + 1,jk  ) + zti(ji  ,jj - 1    &   
(  348)      &                  ,jk  ) + zti(ji  ,jj - 1,jk - 1) + rtrn))    &   
(  349)      &                  /(e1t(ji,jj)*e2t(ji,jj)*fse3w(ji,jj,jk))    &   
(  350)      &                  *umask(ji - 1,jj,jk  )*umask(ji  ,jj  ,jk  )    &   
(  351)      &                  *umask(ji  ,jj,jk- 1)*umask(ji - 1,jj  ,jk - 1)    &   
(  352)      &                  *vmask(ji  ,jj- 1,jk)*vmask(ji  ,jj  ,jk  )    &   
(  353)      &                  *vmask(ji  ,jj  ,jk-1)*vmask(ji  ,jj - 1,jk - 1)       
(  354)                   END DO
(  355)                 END DO
(  356)               END DO
(  357) 
(  358)               DO jj = 2,jpjm1
(  359)                 DO ji = 2,jpim1
(  360)                   zx(ji,jj,1) = zx(ji,jj,1)    &   
(  361)      &                - 0.5*rdttrc(jk)*rsc*zaa(ji,jj,1)*0.25*    &   
(  362)      &                ( (zbb(ji  ,jj - 1,1  ) + zbb(ji + 1,jj - 1,1  )    &   
(  363)      &                + zbb(ji + 1,jj  ,1  ) + zbb(ji  ,jj  ,1  ))    &   
(  364)      &                *(zti(ji  ,jj + 1,1  ) + zti(ji + 1,jj + 1,1  )    &   
(  365)      &                - zti(ji + 1,jj - 1,1  ) - zti(ji  ,jj - 1,1  ))    &   
(  366)      &                /(zti(ji  ,jj + 1,1  ) + zti(ji + 1,jj + 1,1  )    &   
(  367)      &                + zti(ji + 1,jj - 1,1  ) + zti(ji  ,jj - 1,1  ) +    &   
(  368)      &                rtrn))/(e1u(ji,jj)*e2u(ji,jj)*fse3u(ji,jj,1))    &   
(  369)      &                *vmask(ji  ,jj - 1,1  )*vmask(ji + 1,jj - 1,1  )    &   
(  370)      &                *vmask(ji + 1,jj  ,1  )*vmask(ji  ,jj  ,1  )    
(  371) 
(  372)                  zy(ji,jj,1) = zy(ji,jj,1)    &   
(  373)      &                - 0.5*rdttrc(jk)*rsc*zbb(ji,jj,1)*0.25*    &   
(  374)      &                ( (zaa(ji-1  ,jj ,1  ) + zaa(ji - 1,jj + 1,1  )    &   
(  375)      &                + zaa(ji ,jj  ,1  ) + zaa(ji  ,jj + 1  ,1  ))    &   
(  376)      &                *(zti(ji + 1,jj + 1,1  ) + zti(ji + 1,jj ,1  )    &   
(  377)      &                - zti(ji - 1,jj + 1,1  ) - zti(ji - 1,jj ,1  ))    &   
(  378)      &                /(zti(ji + 1,jj + 1,1  ) + zti(ji + 1,jj ,1  )    &   
(  379)      &                + zti(ji - 1,jj + 1,1  ) + zti(ji - 1,jj ,1  ) +    &   
(  380)      &                rtrn))/(e1v(ji,jj)*e2v(ji,jj)*fse3v(ji,jj,1))    &   
(  381)      &                *umask(ji - 1,jj,1  )*umask(ji - 1,jj + 1,1  )    &   
(  382)      &                *umask(ji    ,jj,1  )*umask(ji  ,jj + 1 ,1  )    
(  383) 
(  384)                 END DO
(  385)               END DO
(  386)           ENDIF
(  387) 
(  388)           ! ... Lateral boundary conditions on z[xyz]
(  389)           CALL lbc_lnk( zx, 'U', -1. )
(  390)           CALL lbc_lnk( zy, 'V', -1. )
(  391)           CALL lbc_lnk( zz, 'W',  1. )
(  392) 
(  393) ! 2.4 reinitialization
(  394) 






PGF90 (Version     10.5)          09/23/2010  09:49:19      page 8

(  395)           DO jk = 1,jpk
(  396)             DO jj = 1,jpj
(  397)               DO ji = 1,jpi
(  398)                 zaa(ji,jj,jk) = zx(ji,jj,jk)
(  399)                 zbb(ji,jj,jk) = zy(ji,jj,jk)
(  400)                 zcc(ji,jj,jk) = zz(ji,jj,jk)
(  401)               END DO
(  402)             END DO
(  403)           END DO
(  404) 
(  405) ! 2.5 calcul of the final field:
(  406) !    advection by antidiffusive mass fluxes and an upstream scheme
(  407) 
(  408)           DO jk = 1,jpk
(  409)              DO jj = 2,jpjm1
(  410)                 DO ji = 2,jpim1
(  411)                    zfp_ui = 0.5 * ( zaa(ji,jj,jk) + ABS( zaa(ji,jj,jk) ) )
(  412)                    zfp_vj = 0.5 * ( zbb(ji,jj,jk) + ABS( zbb(ji,jj,jk) ) )
(  413)                    zfm_ui = 0.5 * ( zaa(ji,jj,jk) - ABS( zaa(ji,jj,jk) ) )
(  414)                    zfm_vj = 0.5 * ( zbb(ji,jj,jk) - ABS( zbb(ji,jj,jk) ) )            
(  415)                    zkx(ji,jj,jk) = zfp_ui * zti(ji,jj,jk) + zfm_ui * zti(ji+1,jj  ,jk) 
(  416)                    zky(ji,jj,jk) = zfp_vj * zti(ji,jj,jk) + zfm_vj * zti(ji    ,jj+1,jk)             
(  417)                 END DO
(  418)              END DO
(  419)           END DO
(  420) 
(  421)           DO jk = 2,jpk
(  422)              DO jj = 1,jpj
(  423)                 DO ji = 1,jpi
(  424)                    zfp_w = 0.5 * ( zcc(ji,jj,jk) + ABS( zcc(ji,jj,jk) ) )
(  425)                    zfm_w = 0.5 * ( zcc(ji,jj,jk) - ABS( zcc(ji,jj,jk) ) )       
(  426)                    zkz(ji,jj,jk) = zfp_w * zti(ji,jj,jk) + zfm_w * zti(ji,jj,jk-1)     
(  427)                 END DO
(  428)              END DO
(  429)           END DO
(  430) 
(  431) 
(  432) ! ... Lateral boundary conditions on zk[xy]
(  433)       CALL lbc_lnk( zkx, 'U', -1. )
(  434)       CALL lbc_lnk( zky, 'V', -1. )
(  435) 
(  436) 
(  437) ! 2.6. calcul of after field using an upstream advection scheme
(  438) 
(  439)           DO jk = 1,jpkm1
(  440)             DO jj = 2,jpjm1
(  441)               DO ji = 2,jpim1
(  442)                 zbtr = 1./(e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk))
(  443)                 ztj(ji,jj,jk) = -zbtr*     &  
(  444)      &              ( zkx(ji,jj,jk) - zkx(ji - 1,jj,jk)    &  
(  445)      &              + zky(ji,jj,jk) - zky(ji,jj - 1,jk)    &  
(  446)      &              + zkz(ji,jj,jk) - zkz(ji,jj,jk + 1) )
(  447) #if defined key_trc_diatrd
(  448)                 IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),1) = trtrd(ji,jj,jk,ikeep(jn),1) -    &  
(  449)      &              zbtr*( zkx(ji,jj,jk) - zkx(ji - 1,jj,jk) )   
(  450) 
(  451)                 IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),2) = trtrd(ji,jj,jk,ikeep(jn),2) -    &  
(  452)      &              zbtr*( zky(ji,jj,jk) - zky(ji,jj - 1,jk) )






PGF90 (Version     10.5)          09/23/2010  09:49:19      page 9

(  453) 
(  454)                 IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),3) = trtrd(ji,jj,jk,ikeep(jn),3) -    &  
(  455)      &              zbtr*( zkz(ji,jj,jk) - zkz(ji,jj,jk + 1) )
(  456) #endif
(  457)               END DO
(  458)             END DO
(  459)           END DO
(  460) 
(  461) ! 2.6 END of antidiffusive correction loop
(  462) 
(  463)         END DO
(  464) 
(  465) ! 3. trend due to horizontal and vertical advection of tracer jn
(  466) ! --------------------------------------------------------------
(  467) 
(  468)         DO jk = 1,jpk
(  469)           DO jj = 2,jpjm1
(  470)             DO ji = 2,jpim1
(  471)               ztra = ( zbuf(ji,jj,jk) + ztj(ji,jj,jk) ) * tmask(ji,jj,jk)
(  472)               tra(ji,jj,jk,jn) = tra(ji,jj,jk,jn) + ztra
(  473)             END DO
(  474)           END DO
(  475)         END DO
(  476) 
(  477) ! 4.0 convert the transport trend into advection trend
(  478) ! ----------------------------------------------------
(  479) 
(  480) #if defined key_trc_diatrd
(  481)         DO jk = 1,jpk
(  482)           DO jj = 2,jpjm1
(  483)             DO  ji = 2,jpim1
(  484)               zbtr = 1./(e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk))
(  485)               zgm = zbtr * trn(ji,jj,jk,jn) *     &  
(  486)      &            ( zun(ji  ,jj,jk) * e2u(ji  ,jj) * fse3u(ji  ,jj,jk)    &  
(  487)      &            -zun(ji-1,jj,jk) * e2u(ji-1,jj) * fse3u(ji-1,jj,jk))
(  488) 
(  489)               zgz = zbtr * trn(ji,jj,jk,jn) *     &  
(  490)      &            ( zvn(ji,jj  ,jk) * e1v(ji,jj  ) * fse3v(ji,jj  ,jk)    &  
(  491)      &            -zvn(ji,jj-1,jk) * e1v(ji,jj-1) * fse3v(ji,jj-1,jk))
(  492) 
(  493)               IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),1) = trtrd(ji,jj,jk,ikeep(jn),1) + zgm
(  494)               IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),2) = trtrd(ji,jj,jk,ikeep(jn),2) + zgz
(  495)               IF (luttrd(jn)) trtrd(ji,jj,jk,ikeep(jn),3) = trtrd(ji,jj,jk,ikeep(jn),3)    &  
(  496)      &            - trn(ji,jj,jk,jn) * hdivn(ji,jj,jk)
(  497)             END DO
(  498)           END DO
(  499)         END DO
(  500) 
(  501)         ! Lateral boundary conditions on trtrd:
(  502) 
(  503)         IF (luttrd(jn)) CALL lbc_lnk( trtrd(:,:,:,ikeep(jn),1), 'T', 1. )
(  504)         IF (luttrd(jn)) CALL lbc_lnk( trtrd(:,:,:,ikeep(jn),2), 'T', 1. )
(  505)         IF (luttrd(jn)) CALL lbc_lnk( trtrd(:,:,:,ikeep(jn),3), 'T', 1. )
(  506) #endif
(  507) 
(  508)  
(  509)         ! END of tracer loop
(  510)         ! ==================






PGF90 (Version     10.5)          09/23/2010  09:49:19      page 10

(  511)      ENDDO
(  512) 
(  513)       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
(  514)          WRITE(charout, FMT="('smolar - adv')")
(  515)          CALL prt_ctl_trc_info(charout)
(  516)          CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
(  517)       ENDIF
(  518)      
(  519)   END SUBROUTINE trc_adv_smolar
(  520) 
(  521) #else
(  522)    !!----------------------------------------------------------------------
(  523)    !!   Default option                                         Empty module
(  524)    !!----------------------------------------------------------------------
(  525) CONTAINS
(  526)    SUBROUTINE trc_adv_smolar( kt ) 
(  527)       INTEGER, INTENT(in) :: kt
(  528) !      WRITE(*,*) 'trc_adv_smolar: You should not have seen this print! error?', kt
(  529)    END SUBROUTINE trc_adv_smolar
(  530) #endif
(  531) 
(  532)    !!======================================================================
(  533) END MODULE trcadv_smolar









































PGF90 (Version     10.5)          09/23/2010  09:49:19      page 11

(    1) # 1 "trcadv_smolar.F90"
(    1) MODULE trcadv_smolar
(    2)    !!==============================================================================
(    3)    !!                       ***  MODULE  trcadv_smolar  ***
(    4)    !! Ocean passive tracers:  horizontal & vertical advective trend
(    5)    !!==============================================================================
(    6) # 522
(  522)    !!----------------------------------------------------------------------
(  523)    !!   Default option                                         Empty module
(  524)    !!----------------------------------------------------------------------
(  525) CONTAINS
(  526)    SUBROUTINE trc_adv_smolar( kt ) 
(  527)       INTEGER, INTENT(in) :: kt
(  528) !      WRITE(*,*) 'trc_adv_smolar: You should not have seen this print! error?', kt
(  529)    END SUBROUTINE trc_adv_smolar
(  530) 
(  531) # 532
(  532)    !!======================================================================
(  533) END MODULE trcadv_smolar
