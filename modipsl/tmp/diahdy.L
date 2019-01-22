


PGF90 (Version     10.5)          09/23/2010  09:49:27      page 1

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

Filename: diahdy.F90

(    1) MODULE diahdy
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  diahdy  ***
(    4)    !! Ocean diagnostics : computation the dynamical heigh
(    5)    !!======================================================================
(    6) #if   defined key_diahdy   ||   defined key_esopa
(    7)    !!----------------------------------------------------------------------
(    8)    !!   'key_diahdy' :                          dynamical heigh diagnostics
(    9)    !!----------------------------------------------------------------------
(   10)    !!   dia_hdy      : dynamical heigh computation
(   11)    !!----------------------------------------------------------------------
(   12)    !! * Modules used
(   13)    USE oce             ! ocean dynamics and tracers
(   14)    USE dom_oce         ! ocean space and time domain
(   15)    USE phycst          ! physical constants
(   16)    USE in_out_manager  ! I/O manager
(   17) 
(   18)    IMPLICIT NONE
(   19)    PRIVATE
(   20) 
(   21)    !! * Routine accessibility
(   22)    PUBLIC dia_hdy     ! called in step.F90 module
(   23) 
(   24)    !! * Shared module variables
(   25)    LOGICAL, PUBLIC, PARAMETER ::   lk_diahdy = .TRUE.   !: dynamical heigh flag
(   26) 
(   27)    !! * Module variables
(   28)    REAL(wp), DIMENSION(jpk) ::   &
(   29)       rhosp         ! ???
(   30) 
(   31)    REAL(wp), DIMENSION(jpi,jpj,jpk) ::   &
(   32)       hdy           ! dynamical heigh
(   33) 
(   34)    !! * Substitutions
(   35) #  include "domzgr_substitute.h90"
(   36)    !!----------------------------------------------------------------------
(   37)    !!   OPA 9.0 , LOCEAN-IPSL (2005) 
(   38)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DIA/diahdy.F90,v 1.3 2005/03/27 18:34:55 opalod Exp $ 
(   39)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   40)    !!----------------------------------------------------------------------
(   41) 
(   42) CONTAINS
(   43)    
(   44)    SUBROUTINE dia_hdy ( kt )
(   45)       !!---------------------------------------------------------------------
(   46)       !!                  ***  ROUTINE dia_hdy  ***






PGF90 (Version     10.5)          09/23/2010  09:49:27      page 2

(   47)       !!
(   48)       !! ** Purpose :   Computes the dynamical heigh
(   49)       !!
(   50)       !! ** Method  : Millero + Poisson
(   51)       !!
(   52)       !! References : 
(   53)       !!	A. E. Gill, atmosphere-ocean dynamics 7.7 pp 215
(   54)       !!
(   55)       !! History :
(   56)       !!        !  9x-xx (P. Delecluse, C. Perigaud)  Original code
(   57)       !!        !  93-10  (C. Perigaud)  a trapezoidal vertical integration 
(   58)       !!                                 consistent WITH the code
(   59)       !!        !  93-12  (G. Madec M. Imbard)
(   60)       !!        !  96-03  (N. Ferry)  integration at t-points
(   61)       !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
(   62)       !!----------------------------------------------------------------------
(   63)       !! * Arguments
(   64)       INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
(   65) 
(   66)       !! * Local declarations
(   67)       INTEGER :: ji, jj, jk
(   68)       INTEGER :: ihdsup, ik
(   69) 
(   70)       REAL(wp) :: zgdsup, za, zb, zciint, zfacto, zhd
(   71)       REAL(wp) :: zp, zh, zt, zs, zxk, zq, zsr, zr1, zr2, zr3, zr4
(   72)       REAL(wp) :: ze, zbw, zc, zd, zaw, zb1, za1, zkw, zk0
(   73)       REAL(wp), DIMENSION(jpi,jpj,jpk) :: zsva
(   74)       REAL(wp), DIMENSION(jpk)         :: zwkx, zwky, zwkz
(   75)       REAL(wp) :: fsatg
(   76)       REAL(wp) :: pfps, pfpt, pfphp  
(   77) 
(   78)       ! Adiabatic laspse rate fsatg, defined as the change of temperature
(   79)       ! per unit pressure for adiabatic change of pressure of an element
(   80)       ! of seawater (bryden,h.,1973,deep-sea res.,20,401-408).
(   81)       ! units:
(   82)       !      pressure        pfphp    decibars
(   83)       !      temperature     pfpt     deg celsius (ipts-68)
(   84)       !      salinity        pfps     (ipss-78)
(   85)       !      adiabatic       fsatg    deg. c/decibar
(   86)       ! checkvalue: atg=3.255976e-4 c/dbar for pfps=40 (ipss-78),
(   87)       ! pfpt=40 deg c, pfphp=10000 decibars
(   88)       
(   89)       fsatg(pfps,pfpt,pfphp)   &
(   90)          = (((-2.1687e-16*pfpt+1.8676e-14)*pfpt-4.6206e-13)*pfphp    &
(   91)          +((2.7759e-12*pfpt-1.1351e-10)*(pfps-35.)+((-5.4481e-14*pfpt    &
(   92)          +8.733e-12)*pfpt-6.7795e-10)*pfpt+1.8741e-8))*pfphp    &
(   93)          +(-4.2393e-8*pfpt+1.8932e-6)*(pfps-35.)    &
(   94)          +((6.6228e-10*pfpt-6.836e-8)*pfpt+8.5258e-6)*pfpt+3.5803e-5
(   95)       !!----------------------------------------------------------------------
(   96) 
(   97)       ! 1. height dynamic
(   98)       ! -----------------
(   99)       ! depth for reference
(  100) 
(  101)       zgdsup = 1500.
(  102)       
(  103)       ! below for hdyn levitus
(  104)       






PGF90 (Version     10.5)          09/23/2010  09:49:27      page 3

(  105)       IF( kt == nit000 ) THEN
(  106)          IF(lwp) WRITE(numout,*)
(  107)          IF(lwp) WRITE(numout,*) 'dia_hdy : computation of dynamical heigh'
(  108)          IF(lwp) WRITE(numout,*) '~~~~~~~'
(  109) # if defined key_s_coord || defined key_partial_steps
(  110)          ! Dynamic height diagnostics  not yet implemented
(  111)          IF(lwp) WRITE(numout,cform_err)
(  112)          IF(lwp) WRITE(numout,*) '          key_s_coord or key_partial_steps used'
(  113)          IF(lwp) WRITE(numout,*) '          Dynamical height diagnostics not yet implemented'
(  114)          nstop = nstop + 1
(  115) # endif
(  116) 
(  117)          DO jk = 1, jpk
(  118)             IF( fsdepw(1,1,jk) > zgdsup ) GOTO 110
(  119)          END DO
(  120)          IF(lwp) WRITE(numout,*)'problem zgdsup greater than gdepw(jpk)'
(  121)          STOP 'dia_hdy'
(  122) 110      CONTINUE
(  123)          ihdsup = jk - 1
(  124)          IF(lwp) WRITE(numout,*)' ihdsup = ', ihdsup
(  125) 
(  126)          ! Interpolation coefficients for zgdsup-gdepw(ihdsup) layer
(  127) 
(  128)          za = fsdepw(1,1,ihdsup  )
(  129)          zb = fsdepw(1,1,ihdsup+1)
(  130)          IF( za > zgdsup .OR. zb < zgdsup ) THEN
(  131)             IF(lwp) WRITE(numout,*) za, zb, ihdsup, zgdsup
(  132)             IF(lwp) WRITE(numout,*) ' bad ihdsup'
(  133)             STOP
(  134)          ENDIF
(  135)          
(  136)          zciint = (zgdsup - za) / (zb - za)
(  137) 
(  138)          ! Computes the specific volume reference in situ temperature
(  139)          
(  140)          DO jk = 1, jpk
(  141)             zp = 0.e0
(  142)             zh = fsdept(1,1,jk)
(  143)             zt = 0.e0
(  144)             zs = 35.
(  145)             zxk= zh * fsatg( zs, zt, zp )
(  146)             zt = zt + 0.5 * zxk
(  147)             zq = zxk
(  148)             zp = zp + 0.5 * zh
(  149)             zxk= zh*fsatg( zs, zt, zp )
(  150)             zt = zt + 0.29289322 * ( zxk - zq )
(  151)             zq = 0.58578644 * zxk + 0.121320344 * zq
(  152)             zxk= zh * fsatg( zs, zt, zp )
(  153)             zt = zt + 1.707106781 * ( zxk - zq )
(  154)             zq = 3.414213562 * zxk - 4.121320344 * zq
(  155)             zp = zp + 0.5 * zh
(  156)             zxk= zh * fsatg( zs, zt, zp )
(  157)             zwkx(jk) = zt + ( zxk - 2.0 * zq ) / 6.0
(  158)          END DO
(  159) 
(  160)          ! In situ density (add the compression terms)
(  161) 
(  162)          DO jk = 1, jpk






PGF90 (Version     10.5)          09/23/2010  09:49:27      page 4

(  163)             zt = zwkx(jk)
(  164)             zs = 35.
(  165)             ! square root salinity
(  166)             zsr = sqrt( abs( zs ) )
(  167)             zwky(jk) = zsr
(  168)             ! compute density pure water at atm pressure
(  169)             zr1= ((((6.536332e-9*zt-1.120083e-6)*zt+1.001685e-4)*zt   &
(  170)                -9.095290e-3)*zt+6.793952e-2)*zt+999.842594
(  171)             ! seawater density atm pressure
(  172)             zr2= (((5.3875e-9*zt-8.2467e-7)*zt+7.6438e-5)*zt   &
(  173)                -4.0899e-3)*zt+8.24493e-1
(  174)             zr3= (-1.6546e-6*zt+1.0227e-4)*zt-5.72466e-3
(  175)             zr4= 4.8314e-4
(  176)             zwkz(jk)= (zr4*zs + zr3*zsr + zr2)*zs + zr1
(  177)          END DO
(  178) 
(  179)          DO jk = 1, jpk
(  180)             zt = zwkx(jk)
(  181)             zs = 35.
(  182)             zsr= zwky(jk)
(  183)             zh = fsdept(1,1,jk)
(  184) 
(  185)             ze = ( 9.1697e-11*zt+2.0816e-9 ) *zt-9.9348e-8
(  186)             zbw= ( 5.2787e-9*zt-6.12293e-7 ) * zt+8.50935e-6
(  187)             zb = zbw + ze * zs
(  188) 
(  189)             zd = 1.91075e-4
(  190)             zc = (-1.6078e-6*zt-1.0981e-5)*zt+2.2838e-3
(  191)             zaw= ((-5.77905e-7*zt+1.16092e-4)*zt+1.43713e-3)*zt+3.239908
(  192)             za = ( zd*zsr + zc)*zs + zaw
(  193) 
(  194)             zb1= (-5.3009e-3*zt+1.6483e-1)*zt+7.944e-1
(  195)             za1= ((-6.1670e-4*zt+1.09987e-1)*zt-6.03459)*zt+546.746
(  196)             zkw= (((-5.155288e-4*zt+1.360477e-1)*zt-23.27105)*zt   &
(  197)                 +1484.206)*zt+196522.1
(  198)             zk0= (zb1*zsr + za1)*zs + zkw
(  199)             ! evaluate pressure polynomial
(  200)             zwkz(jk) = zwkz(jk) / ( 1.0 - zh / ( zk0+zh*(za+zb*zh) ) )
(  201)          END DO
(  202) 
(  203)          DO jk = 1, jpk
(  204)             rhosp(jk) = zwkz(jk)
(  205)          END DO
(  206)       ENDIF
(  207) 
(  208)       ! Computes the specific volume anomaly
(  209) 
(  210)       DO jk = 1, jpkm1
(  211)          DO jj = 1, jpj
(  212)             DO ji = 1, jpi
(  213)                IF( tmask(ji,jj,jk) /= 0. ) THEN
(  214)                   zsva(ji,jj,jk) = ( rau0*rhd(ji,jj,jk)+rau0 -rhosp(jk) ) / rhosp(jk)
(  215)                ELSE
(  216)                   zsva(ji,jj,jk)=0.
(  217)                ENDIF
(  218)             END DO
(  219)          END DO
(  220)       END DO






PGF90 (Version     10.5)          09/23/2010  09:49:27      page 5

(  221) 
(  222)       ! zfacto coefficient to cmg
(  223)       
(  224)       ! zfacto= 1.  e+2
(  225)       !           mg->cmg
(  226)       zfacto = 1.0 * 1.e2
(  227)       
(  228)       ! Fisrt compute at depth ik=ihdsup
(  229)       
(  230)       ik = ihdsup
(  231)       DO jj = 1, jpj
(  232)          DO ji = 1, jpi
(  233)             zhd = zfacto * zciint * fse3t(ji,jj,ik) * zsva(ji,jj,ik)
(  234)             hdy(ji,jj,ik) = zhd * tmask(ji,jj,ik) * tmask(ji,jj,ik-1)
(  235)          END DO
(  236)       END DO
(  237)       
(  238)       ! Then compute other terms except level jk=1
(  239)       
(  240)       DO jk = ihdsup-1, 2, -1
(  241)          DO jj = 1, jpj
(  242)             DO ji = 1, jpi
(  243)                zhd = hdy(ji,jj,jk+1) + zfacto * fse3t(ji,jj,jk) * zsva(ji,jj,jk)
(  244)                hdy(ji,jj,jk) = zhd * tmask(ji,jj,jk) * tmask(ji,jj,jk-1)
(  245)             END DO
(  246)          END DO
(  247)       END DO
(  248)       
(  249)       ! Then compute other the last layer term jk=1
(  250)       
(  251)       ik = 1
(  252)       DO jj = 1, jpj
(  253)          DO ji = 1, jpi
(  254)             zhd = hdy(ji,jj,ik+1) + zfacto * fse3t(ji,jj,ik) * zsva(ji,jj,ik)
(  255)             hdy(ji,jj,ik) = zhd * tmask(ji,jj,ik)
(  256)          END DO
(  257)       END DO
(  258) 
(  259)    END SUBROUTINE dia_hdy
(  260) 
(  261) #else
(  262)    !!----------------------------------------------------------------------
(  263)    !!   Default option :                       NO dynamic heigh diagnostics
(  264)    !!----------------------------------------------------------------------
(  265)    USE in_out_manager
(  266)    LOGICAL, PUBLIC, PARAMETER ::   lk_diahdy = .FALSE.   !: dynamical heigh flag
(  267) CONTAINS
(  268)    SUBROUTINE dia_hdy( kt )               ! Empty routine
(  269)       if(lwp) WRITE(numout,*) 'diahdy: You should not have seen this print! error?', kt
(  270)    END SUBROUTINE dia_hdy
(  271) #endif
(  272) 
(  273)    !!======================================================================
(  274) END MODULE diahdy










PGF90 (Version     10.5)          09/23/2010  09:49:27      page 6

(    1) # 1 "diahdy.F90"
(    1) MODULE diahdy
(    2)    !!======================================================================
(    3)    !!                       ***  MODULE  diahdy  ***
(    4)    !! Ocean diagnostics : computation the dynamical heigh
(    5)    !!======================================================================
(    6) # 262
(  262)    !!----------------------------------------------------------------------
(  263)    !!   Default option :                       NO dynamic heigh diagnostics
(  264)    !!----------------------------------------------------------------------
(  265)    USE in_out_manager
(  266)    LOGICAL, PUBLIC, PARAMETER ::   lk_diahdy = .FALSE.   !: dynamical heigh flag
(  267) CONTAINS
(  268)    SUBROUTINE dia_hdy( kt )               ! Empty routine
(  269)       if(lwp) WRITE(numout,*) 'diahdy: You should not have seen this print! error?', kt
(  270)    END SUBROUTINE dia_hdy
(  271) 
(  272) # 273
(  273)    !!======================================================================
(  274) END MODULE diahdy
