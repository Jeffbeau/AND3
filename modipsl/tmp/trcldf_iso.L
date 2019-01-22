


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

Filename: trcldf_iso.F90

(    1) MODULE trcldf_iso
(    2)    !!==============================================================================
(    3)    !!                    ***  MODULE  trcldf_iso  ***
(    4)    !! Ocean passive tracers:  horizontal component of the lateral tracer mixing trend
(    5)    !!==============================================================================
(    6) #if key_passivetrc && defined key_ldfslp 
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_ldfslp'                  rotation of the lateral mixing tensor
(    9)    !!----------------------------------------------------------------------
(   10)    !!   trc_ldf_iso : update the tracer trend with the horizontal component
(   11)    !!                 of iso neutral laplacian operator or horizontal 
(   12)    !!                 laplacian operator in s-coordinate
(   13)    !!----------------------------------------------------------------------
(   14)    !! * Modules used
(   15)    USE oce_trc      ! ocean dynamics and tracers variables
(   16)    USE trc          ! ocean passive tracers variables
(   17)    USE prtctl_trc   ! Print control for debbuging
(   18) 
(   19)    IMPLICIT NONE
(   20)    PRIVATE
(   21) 
(   22)    !! * Routine accessibility
(   23)    PUBLIC trc_ldf_iso  ! routine called by step.F90
(   24) 
(   25)    !! * Substitutions
(   26) #  include "passivetrc_substitute.h90"
(   27)    !!----------------------------------------------------------------------
(   28)    !!   TOP 1.0 , LOCEAN-IPSL (2005) 
(   29)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/TRP/trcldf_iso.F90,v 1.9 2006/04/10 15:38:54 opalod Exp $ 
(   30)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   31)    !!----------------------------------------------------------------------
(   32) 
(   33) CONTAINS
(   34) 
(   35)    SUBROUTINE trc_ldf_iso( kt )
(   36)       !!----------------------------------------------------------------------
(   37)       !!                  ***  ROUTINE trc_ldf_iso  ***
(   38)       !! 
(   39)       !! ** Purpose :   Compute the before horizontal tracer  diffusive 
(   40)       !!      trend and add it to the general trend of tracer equation.
(   41)       !!
(   42)       !! ** Method  :   The horizontal component of the lateral diffusive trends 
(   43)       !!      is provided by a 2nd order operator rotated along neural or geopo-
(   44)       !!      tential surfaces to which an eddy induced advection can be added
(   45)       !!      It is computed using before fields (forward in time) and isopyc-
(   46)       !!      nal or geopotential slopes computed in routine ldfslp.






PGF90 (Version     10.5)          09/23/2010  09:49:20      page 2

(   47)       !!
(   48)       !!      horizontal fluxes associated with the rotated lateral mixing:
(   49)       !!         zftu = (aht+ahtb0) e2u*e3u/e1u di[ tb ]
(   50)       !!               - aht       e2u*uslp    dk[ mi(mk(tb)) ]
(   51)       !!         zftv = (aht+ahtb0) e1v*e3v/e2v dj[ tb ]
(   52)       !!               - aht       e2u*vslp    dk[ mj(mk(tb)) ]
(   53)       !!      add horizontal Eddy Induced advective fluxes (lk_traldf_eiv=T):
(   54)       !!         zftu = zftu - dk-1[ aht e2u mi(wslpi) ] mi( tb ) 
(   55)       !!         zftv = zftv - dk-1[ aht e1v mj(wslpj) ] mj( tb ) 
(   56)       !!      take the horizontal divergence of the fluxes:
(   57)       !!         difft = 1/(e1t*e2t*e3t) {  di-1[ zftu ] +  dj-1[ zftv ]  }
(   58)       !!      Add this trend to the general trend tra :
(   59)       !!         tra = tra + difft
(   60)       !!
(   61)       !! ** Action  : - Update tra arrays with the before isopycnal or
(   62)       !!                geopotential s-coord harmonic mixing trend.
(   63)       !!              - Save the trends in trtrd ('key_trc_diatrd')
(   64)       !!
(   65)       !! History :
(   66)       !!        !  94-08  (G. Madec, M. Imbard)
(   67)       !!        !  97-05  (G. Madec)  split into traldf and trazdf
(   68)       !!        !  98-03  (L. Bopp, MA Foujols) passive tracer generalisation
(   69)       !!        !  00-10  (MA Foujols E Kestenare) USE passive tracer coefficient
(   70)       !!   8.5  !  02-08  (G. Madec)  Free form, F90
(   71)       !!   9.0  !  04-03  (C. Ethe)  Free form, F90
(   72)       !!----------------------------------------------------------------------
(   73)       !! * Modules used
(   74)       USE oce_trc       , zftu => ua,  &  ! use ua as workspace
(   75)          &                zfsu => va      ! use va as workspace
(   76) 
(   77)       !! * Arguments
(   78)       INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
(   79) 
(   80)       !! * Local declarations
(   81)       INTEGER ::   ji, jj, jk,jn             ! dummy loop indices
(   82)       REAL(wp) ::   &
(   83)          zabe1, zabe2, zcof1, zcof2,   &  ! temporary scalars
(   84)          zmsku, zmskv, zbtr,           &
(   85) #if defined key_trcldf_eiv
(   86)          zcg1, zcg2, zuwk, zvwk,       &
(   87)          zuwk1, zvwk1,                 &
(   88) #endif
(   89)          ztra
(   90) 
(   91)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(   92)          zdkt, zdk1t            ! workspace
(   93) 
(   94) #if defined key_trcldf_eiv
(   95)       REAL(wp), DIMENSION(jpi,jpj) ::   &
(   96)          zftug, zftvg
(   97) #endif
(   98) 
(   99)       REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(  100)          zftv                       ! workspace
(  101)       CHARACTER (len=22) :: charout
(  102)       !!----------------------------------------------------------------------
(  103) 
(  104)       IF( kt == nittrc000 ) THEN






PGF90 (Version     10.5)          09/23/2010  09:49:20      page 3

(  105)          IF(lwp) WRITE(numout,*)
(  106)          IF(lwp) WRITE(numout,*) 'trc_ldf_iso : iso neutral lateral diffusion or'
(  107)          IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   horizontal laplacian diffusion in s-coordinate'
(  108) #if defined key_trcldf_eiv && defined key_diaeiv
(  109)          u_trc_eiv(:,:,:) = 0.e0
(  110)          v_trc_eiv(:,:,:) = 0.e0
(  111) #endif
(  112)       ENDIF
(  113) 
(  114) 
(  115)       DO jn = 1, jptra
(  116) 
(  117)          !                                                ! ===============
(  118)          DO jk = 1, jpkm1                                 ! Horizontal slab
(  119)             !                                             ! ===============
(  120)             ! 1. Vertical tracer gradient at level jk and jk+1
(  121)             ! ------------------------------------------------
(  122)             ! surface boundary condition: zdkt(jk=1)=zdkt(jk=2)
(  123) 
(  124)             zdk1t(:,:) = ( trb(:,:,jk,jn) - trb(:,:,jk+1,jn) ) * tmask(:,:,jk+1)
(  125) 
(  126)             IF( jk == 1 ) THEN
(  127)                zdkt(:,:) = zdk1t(:,:)
(  128)             ELSE
(  129)                zdkt(:,:) = ( trb(:,:,jk-1,jn) - trb(:,:,jk,jn) ) * tmask(:,:,jk)
(  130)             ENDIF
(  131) 
(  132) 
(  133)             ! 2. Horizontal fluxes
(  134)             ! --------------------
(  135) 
(  136)             DO jj = 1 , jpjm1
(  137)                DO ji = 1, fs_jpim1   ! vector opt.
(  138)                   zabe1 = ( fsahtru(ji,jj,jk) + ahtrb0 ) * e2u(ji,jj) * fse3u(ji,jj,jk) / e1u(ji,jj)
(  139)                   zabe2 = ( fsahtrv(ji,jj,jk) + ahtrb0 ) * e1v(ji,jj) * fse3v(ji,jj,jk) / e2v(ji,jj)
(  140) 
(  141)                   zmsku = 1. / MAX(   tmask(ji+1,jj,jk  ) + tmask(ji,jj,jk+1)   &
(  142)                      + tmask(ji+1,jj,jk+1) + tmask(ji,jj,jk  ), 1. )
(  143) 
(  144)                   zmskv = 1. / MAX(   tmask(ji,jj+1,jk  ) + tmask(ji,jj,jk+1)   &
(  145)                      + tmask(ji,jj+1,jk+1) + tmask(ji,jj,jk  ), 1. )
(  146) 
(  147)                   zcof1 = -fsahtru(ji,jj,jk) * e2u(ji,jj) * uslp(ji,jj,jk) * zmsku
(  148)                   zcof2 = -fsahtrv(ji,jj,jk) * e1v(ji,jj) * vslp(ji,jj,jk) * zmskv
(  149) 
(  150)                   zftu(ji,jj,jk) = umask(ji,jj,jk) * (   zabe1 * (   trb(ji+1,jj,jk,jn) - trb(ji,jj,jk,jn)  )   &
(  151)                      &                              + zcof1 * (   zdkt (ji+1,jj) + zdk1t(ji,jj)      &
(  152)                      &                                          + zdk1t(ji+1,jj) + zdkt (ji,jj)  )  )
(  153) 
(  154)                   zftv(ji,jj,jk) = vmask(ji,jj,jk) * (   zabe2 * (   trb(ji,jj+1,jk,jn) - trb(ji,jj,jk,jn)  )   &
(  155)                      &                              + zcof2 * (   zdkt (ji,jj+1) + zdk1t(ji,jj)      &
(  156)                      &                                          + zdk1t(ji,jj+1) + zdkt (ji,jj)  )  )
(  157) 
(  158)                END DO
(  159)             END DO
(  160) 
(  161) #   if defined key_trcldf_eiv
(  162)             !                              ! ---------------------------------------!






PGF90 (Version     10.5)          09/23/2010  09:49:20      page 4

(  163)             !                              ! Eddy induced vertical advective fluxes !
(  164)             !                              ! ---------------------------------------!
(  165)             DO jj = 1, jpjm1
(  166)                DO ji = 1, fs_jpim1   ! vector opt.
(  167)                   zuwk = ( wslpi(ji,jj,jk  ) + wslpi(ji+1,jj,jk  ) ) * fsaeitru(ji,jj,jk  ) * umask(ji,jj,jk  )
(  168)                   zuwk1= ( wslpi(ji,jj,jk+1) + wslpi(ji+1,jj,jk+1) ) * fsaeitru(ji,jj,jk+1) * umask(ji,jj,jk+1)
(  169)                   zvwk = ( wslpj(ji,jj,jk  ) + wslpj(ji,jj+1,jk  ) ) * fsaeitrv(ji,jj,jk  ) * vmask(ji,jj,jk  )
(  170)                   zvwk1= ( wslpj(ji,jj,jk+1) + wslpj(ji,jj+1,jk+1) ) * fsaeitrv(ji,jj,jk+1) * vmask(ji,jj,jk+1)
(  171) 
(  172)                   zcg1= -0.25 * e2u(ji,jj) * umask(ji,jj,jk) * ( zuwk-zuwk1 )
(  173)                   zcg2= -0.25 * e1v(ji,jj) * vmask(ji,jj,jk) * ( zvwk-zvwk1 )
(  174) 
(  175)                   zftug(ji,jj) = zcg1 * ( trb(ji+1,jj,jk,jn) + trb(ji,jj,jk,jn) )
(  176)                   zftvg(ji,jj) = zcg2 * ( trb(ji,jj+1,jk,jn) + trb(ji,jj,jk,jn) )
(  177) 
(  178)                   zftu(ji,jj,jk) = zftu(ji,jj,jk) + zftug(ji,jj)
(  179)                   zftv(ji,jj,jk) = zftv(ji,jj,jk) + zftvg(ji,jj)
(  180) 
(  181) #   if defined key_diaeiv
(  182)                   u_trc_eiv(ji,jj,jk) = -2. * zcg1 / ( e2u(ji,jj) * fse3u(ji,jj,jk) )
(  183)                   v_trc_eiv(ji,jj,jk) = -2. * zcg2 / ( e1v(ji,jj) * fse3v(ji,jj,jk) )
(  184) #   endif
(  185)                END DO
(  186)             END DO
(  187) #   endif
(  188) 
(  189)             ! II.4 Second derivative (divergence) and add to the general trend
(  190)             ! ----------------------------------------------------------------
(  191) 
(  192)             DO jj = 2 , jpjm1
(  193)                DO ji = fs_2, fs_jpim1   ! vector opt.
(  194)                   zbtr= 1. / ( e1t(ji,jj)*e2t(ji,jj)*fse3t(ji,jj,jk) )
(  195)                   ztra = zbtr * (  zftu(ji,jj,jk) - zftu(ji-1,jj  ,jk)   &
(  196)                      &          + zftv(ji,jj,jk) - zftv(ji  ,jj-1,jk)  )
(  197)                   tra (ji,jj,jk,jn) = tra (ji,jj,jk,jn) + ztra
(  198) #if defined key_trc_diatrd
(  199)                   IF (luttrd(jn)) trtrd (ji,jj,jk,ikeep(jn),4) = ( zftu(ji,jj,jk) - zftu(ji-1,jj,jk  ) ) * zbtr
(  200)                   IF (luttrd(jn)) trtrd (ji,jj,jk,ikeep(jn),5) = ( zftv(ji,jj,jk) - zftv(ji,jj-1,jk  ) ) * zbtr
(  201) #endif
(  202)                END DO
(  203)             END DO
(  204)             !                                          ! ===============
(  205)          END DO                                        !   End of slab  
(  206)          !                                             ! ===============
(  207) 
(  208)       END DO
(  209) 
(  210)       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
(  211)          WRITE(charout, FMT="('ldf - iso')")
(  212)          CALL prt_ctl_trc_info(charout)
(  213)          CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm,clinfo2='trd')
(  214)       ENDIF
(  215) 
(  216)    END SUBROUTINE trc_ldf_iso
(  217) 
(  218) #else
(  219)    !!----------------------------------------------------------------------
(  220)    !!   Dummy module :             No rotation of the lateral mixing tensor






PGF90 (Version     10.5)          09/23/2010  09:49:20      page 5

(  221)    !!----------------------------------------------------------------------
(  222) CONTAINS
(  223)    SUBROUTINE trc_ldf_iso( kt )               ! Empty routine
(  224)       INTEGER, INTENT(in) :: kt
(  225) !      WRITE(*,*) 'trc_ldf_iso: You should not have seen this print! error?', kt
(  226)    END SUBROUTINE trc_ldf_iso
(  227) #endif
(  228) 
(  229)    !!==============================================================================
(  230) END MODULE trcldf_iso






















































PGF90 (Version     10.5)          09/23/2010  09:49:20      page 6

(    1) # 1 "trcldf_iso.F90"
(    1) MODULE trcldf_iso
(    2)    !!==============================================================================
(    3)    !!                    ***  MODULE  trcldf_iso  ***
(    4)    !! Ocean passive tracers:  horizontal component of the lateral tracer mixing trend
(    5)    !!==============================================================================
(    6) # 219
(  219)    !!----------------------------------------------------------------------
(  220)    !!   Dummy module :             No rotation of the lateral mixing tensor
(  221)    !!----------------------------------------------------------------------
(  222) CONTAINS
(  223)    SUBROUTINE trc_ldf_iso( kt )               ! Empty routine
(  224)       INTEGER, INTENT(in) :: kt
(  225) !      WRITE(*,*) 'trc_ldf_iso: You should not have seen this print! error?', kt
(  226)    END SUBROUTINE trc_ldf_iso
(  227) 
(  228) # 229
(  229)    !!==============================================================================
(  230) END MODULE trcldf_iso
