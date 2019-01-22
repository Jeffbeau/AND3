


PGF90 (Version     10.5)          09/23/2010  09:49:41      page 1

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

Filename: rivers.F90

(    1) !!DB:2009.07.30
(    2) !!Interior river input routines
(    3) !!Idea: put all routines in here and call them from the various 
(    4) !!routines that need them.
(    5) !!For SOPA-MC these are: traadv_tvd advection scheme, and dynspg_ts
(    6) 
(    7) !!NB: SLE = St.Lawrence Estuary is the only river input included
(    8) !!    as of 2009.07.30
(    9) !!NB: a time-dependent ramp factor is used (standard for SOPA-MC code)
(   10) !!    This must be changed if used outside of SOPA-MC code
(   11) 
(   12) #ifdef key_RIVER_INPUT
(   13) MODULE rivers
(   14) 
(   15)    !! * Modules used
(   16)    USE oce             ! ocean dynamics and active tracers
(   17)    USE dom_oce         ! ocean space and time domain
(   18)    USE in_out_manager  ! I/O manager
(   19)    USE lib_mpp
(   20)    USE lbclnk          ! ocean lateral boundary condition (or mpp link) 
(   21)    USE daymod
(   22)    USE phycst
(   23) 
(   24)    IMPLICIT NONE
(   25)    PRIVATE
(   26) 
(   27) !!DB: 2009.07.30
(   28)    REAL(wp), PUBLIC :: rivsum, vel_rivsum, da_rivsum, dz_rivsum
(   29) 
(   30)    PUBLIC river_SLE_01, river_SLE_02, river_SLE_03
(   31) 
(   32) 
(   33)    !! * Accessibility
(   34) !   PUBLIC tra_adv_tvd    ! routine called by step.F90
(   35) 
(   36)    !! * Substitutions
(   37) #  include "domzgr_substitute.h90"
(   38) #  include "vectopt_loop_substitute.h90"
(   39) 
(   40) CONTAINS
(   41) 
(   42) !!DB:
(   43) !!First routine, to be called by advection scheme
(   44) !! (NB: SOPA on MC-domain only works with traadv_tvd scheme)
(   45)    SUBROUTINE river_SLE_01( kt )
(   46) 






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 2

(   47)       !! * Arguments
(   48)       INTEGER, INTENT( in ) ::   kt         ! ocean time-step
(   49) 
(   50)       !! * Local declarations
(   51)       INTEGER  ::   ji, jj, jk              ! dummy loop indices
(   52) 
(   53) !!DB: 2009.07.24 -- SLE runoff related (see also dynspg_ts.)
(   54) !!    Put a u,v +/-  S,T at head of SLE
(   55) !!(1) Crude fit to 1914-2007 monthly rivsum data -->
(   56) !!    If rivsum.dat exists routine returns data value;
(   57) !!    else routine returns fit-to-climatological-rivsum-data
(   58) !!    i.e.: rivsum = 10. + 7.*exp(-(float(nday_year)-121.5)**2/37.5**2) !!mSv
(   59)       call rivsum_val(kt,rivsum)
(   60) 
(   61) !!(2) Convert to vel in m/s based on hardwired entry cell configuration
(   62)       da_rivsum = 0.0
(   63)       do jj = mj0(153),mj1(155)      
(   64)          do ji = mi0(14), mi1(14)       
(   65)             da_rivsum = da_rivsum + e2u(ji,jj)
(   66)          enddo
(   67)       enddo
(   68)       call mpp_sum(da_rivsum)     !!to be safe, broadcast this to all processors
(   69)       dz_rivsum = e3t(1)+e3t(2)  !!NB: inaccurate if dz=dz(i,j,k)
(   70)       vel_rivsum = (rivsum*1000.)/(da_rivsum*dz_rivsum) 
(   71) !!DBG: output a value 
(   72)       if(lwp .AND. mod(kt-nit000,int(rday/rdt)) == 0) write(2900+narea,*)kt, nyear, nday_year, rivsum, vel_rivsum
(   73) 
(   74)       do jj = mj0(153),mj1(155)      
(   75)          do ji = mi0(14), mi1(14)       !!i=15 is t-point ===> set u(14) 
(   76)             do jk = 1, 2
(   77)                un(ji,jj,jk) = ramp * vel_rivsum
(   78)                tn(ji,jj,jk) = tn(ji+1,jj,jk)  !no flux for T
(   79)                tb(ji,jj,jk) = tb(ji+1,jj,jk)  !no flux for T
(   80)                sn(ji,jj,jk) = 1.0   ! 1.0 is usual val; 5.0 ; 0.0 (JC recommendation)
(   81)                sb(ji,jj,jk) = 1.0   ! 1.0
(   82)             enddo
(   83)          enddo
(   84)       enddo
(   85) 
(   86) 
(   87)    END SUBROUTINE river_SLE_01
(   88) 
(   89) 
(   90) !!2nd routine; called by dynspg_ts
(   91)    SUBROUTINE river_SLE_02(zua,zub,zwx)
(   92) 
(   93)       !! * Arguments
(   94)      REAL(wp), DIMENSION(jpi,jpj) ::  &
(   95)           zua, zub, zwx
(   96) 
(   97)           !! * Local declarations
(   98)      INTEGER  ::   ji, jj, jk              ! dummy loop indices
(   99) 
(  100) !!DB: Adjust the above quantities to reflect runoff at head of SLE
(  101)       do jj = mj0(153),mj1(155)
(  102)          do ji = mi0(14), mi1(14)       !!i=15 is t-point ===> set u(14) 
(  103)             do jk = 1, 2
(  104)                zua(ji,jj) = zua(ji,jj) + fse3u(ji,jj,jk) * ramp * vel_rivsum






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 3

(  105)                zub(ji,jj) = zub(ji,jj) + fse3u(ji,jj,jk) * ramp * vel_rivsum
(  106)                zwx(ji,jj) = zwx(ji,jj) + e2u(ji,jj) * fse3u(ji,jj,jk)  * vel_rivsum
(  107)             enddo
(  108)          enddo
(  109)       enddo
(  110) 
(  111)    END SUBROUTINE river_SLE_02
(  112) 
(  113) 
(  114) !!3rd routine; called by dynspg_ts
(  115)    SUBROUTINE river_SLE_03(ua_e)
(  116)       !! * Arguments
(  117)      REAL(wp), DIMENSION(jpi,jpj) ::  &
(  118)           ua_e
(  119) 
(  120)           !! * Local declarations
(  121)      INTEGER  ::   ji, jj, jk              ! dummy loop indices
(  122) 
(  123)       do jj = mj0(153),mj1(155)
(  124)          do ji = mi0(14), mi1(14)       !!i=15 is t-point ===> set u(14) 
(  125)             ua_e(ji,jj) = ramp * vel_rivsum * (fse3u(ji,jj,1)+fse3u(ji,jj,2))
(  126)          enddo
(  127)       enddo
(  128) 
(  129)    END SUBROUTINE river_SLE_03
(  130) 
(  131) 
(  132) !!DB
(  133) !!Return rivsum value
(  134) !!Routine looks for datafile: rivsum.dat which is of the form
(  135) !! year month day real_time rivsum(m3/s)
(  136) !! (NB: 1914...2007 = 94 years x 12 months = 1128 values)
(  137) !!and extracts the value based on nyear nmonth.
(  138) !!TO DO: interpolate monthly values to nday_year
(  139) !!If the file does not exist (or key_CORE_ANNUAL = off), then it
(  140) !!returns the climatological fit
(  141)    SUBROUTINE rivsum_val(kt, rivsum )
(  142)      
(  143)      !! * Arguments
(  144)      INTEGER, INTENT( in ) ::  kt
(  145)      REAL(wp), INTENT( out ) ::  rivsum
(  146) 
(  147) #if defined key_CORE_ANNUAL
(  148) !     REAL(wp), SAVE :: rivsum_data(94*12)
(  149)      REAL(wp), SAVE :: rivsum_data(120*12)  !Choose a big number
(  150)      REAL(wp) :: c1, c2, c3, c4
(  151)      INTEGER  :: i, recnum
(  152)      INTEGER, SAVE  ::  lastrec
(  153) #endif
(  154) 
(  155) 
(  156) !!(1) Crude fit to 1914-2007 monthly rivsum data
(  157) !!    Guarantees a value for rivsum
(  158)      rivsum = 10. + 7.*exp(-(float(nday_year)-121.5)**2/37.5**2) !!mSv
(  159) 
(  160) #if defined key_CORE_ANNUAL
(  161)      if(kt == nit000) then
(  162)         rivsum_data(:) = 0.0






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 4

(  163)         lastrec = 120*12
(  164)         open(66,file='rivsum.dat',STATUS='OLD',ERR=666)
(  165)         do i = 1, lastrec
(  166)            read(66,*,END=555)c1,c2,c3,c4,rivsum_data(i)
(  167)         enddo
(  168) 555     close(66)
(  169)         lastrec = i-1
(  170)         if(lwp) write(numout2,*)'DBG: RIVERS: Read rivsum.dat; lastrec = ', lastrec
(  171) 666     continue
(  172)         if(rivsum_data(1) == 0.0) then
(  173)            if(lwp) write(numout2,*)'DBG: RIVERS: Cannot find rivsum.dat ===> using climatology'
(  174)         endif
(  175)      endif  !!kt==nit000
(  176)      
(  177) !!DB: Check if there is data; YES ===> assign to rivsum; NO ===> use climatology above
(  178) !!NB: if current year.mon not_in_file ===> use climatology above
(  179)      if(rivsum_data(1) /= 0.0) then
(  180)         recnum = (nyear-1914)*12 + nmonth
(  181)         if(recnum <= lastrec) then
(  182)            rivsum = rivsum_data(recnum)/1000.0 !!mSv required by other code
(  183)         endif
(  184)      endif
(  185) 
(  186) #endif     
(  187)      
(  188)    END SUBROUTINE rivsum_val
(  189) 
(  190) 
(  191)    !!======================================================================
(  192) END MODULE rivers
(  193) 
(  194) #endif
































PGF90 (Version     10.5)          09/23/2010  09:49:41      page 5

(    1) # 1 "rivers.F90"
(    1) !!DB:2009.07.30
(    2) !!Interior river input routines
(    3) !!Idea: put all routines in here and call them from the various 
(    4) !!routines that need them.
(    5) !!For SOPA-MC these are: traadv_tvd advection scheme, and dynspg_ts
(    6) 
(    7) !!NB: SLE = St.Lawrence Estuary is the only river input included
(    8) !!    as of 2009.07.30
(    9) !!NB: a time-dependent ramp factor is used (standard for SOPA-MC code)
(   10) !!    This must be changed if used outside of SOPA-MC code
(   11) 
(   12) # 13
(   13) MODULE rivers
(   14) 
(   15)    !! * Modules used
(   16)    USE oce             ! ocean dynamics and active tracers
(   17)    USE dom_oce         ! ocean space and time domain
(   18)    USE in_out_manager  ! I/O manager
(   19)    USE lib_mpp
(   20)    USE lbclnk          ! ocean lateral boundary condition (or mpp link) 
(   21)    USE daymod
(   22)    USE phycst
(   23) 
(   24)    IMPLICIT NONE
(   25)    PRIVATE
(   26) 
(   27) !!DB: 2009.07.30
(   28)    REAL(wp), PUBLIC :: rivsum, vel_rivsum, da_rivsum, dz_rivsum
(   29) 
(   30)    PUBLIC river_SLE_01, river_SLE_02, river_SLE_03
(   31) 
(   32) 
(   33)    !! * Accessibility
(   34) !   PUBLIC tra_adv_tvd    ! routine called by step.F90
(   35) 
(   36)    !! * Substitutions
(   37) # 1 "./domzgr_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                    ***  domzgr_substitute.h90   ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute fsdep. and fse.., the vert. depth and scale
(    5)    !!      factors depending on the vertical coord. used, using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/domzgr_substitute.h90,v 1.2 2005/03/27 18:34:57 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 46
(   46)    !! z-coord:  substitution  fsdep.(,,) ==>  gdep()
(   47)    !!                         fse3.(,,)  ==>  e3.()
(   48) # 38 "rivers.F90"
(   38) # 1 "./vectopt_loop_substitute.h90"
(    1)    !!----------------------------------------------------------------------
(    2)    !!                   ***  vectopt_loop_substitute  ***
(    3)    !!----------------------------------------------------------------------
(    4)    !! ** purpose :   substitute the inner loop starting and inding indices 






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 6

(    5)    !!      to allow unrolling of do-loop using CPP macro.
(    6)    !!----------------------------------------------------------------------
(    7)    !!----------------------------------------------------------------------
(    8)    !!  OPA 9.0 , LOCEAN-IPSL (2005) 
(    9)    !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/vectopt_loop_substitute.h90,v 1.2 2005/03/27 18:34:49 opalod Exp $ 
(   10)    !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
(   11)    !!----------------------------------------------------------------------
(   12) # 39 "rivers.F90"
(   39) 
(   40) # 40
(   40) CONTAINS
(   41) 
(   42) !!DB:
(   43) !!First routine, to be called by advection scheme
(   44) !! (NB: SOPA on MC-domain only works with traadv_tvd scheme)
(   45)    SUBROUTINE river_SLE_01( kt )
(   46) 
(   47)       !! * Arguments
(   48)       INTEGER, INTENT( in ) ::   kt         ! ocean time-step
(   49) 
(   50)       !! * Local declarations
(   51)       INTEGER  ::   ji, jj, jk              ! dummy loop indices
(   52) 
(   53) !!DB: 2009.07.24 -- SLE runoff related (see also dynspg_ts.)
(   54) !!    Put a u,v +/-  S,T at head of SLE
(   55) !!(1) Crude fit to 1914-2007 monthly rivsum data -->
(   56) !!    If rivsum.dat exists routine returns data value;
(   57) !!    else routine returns fit-to-climatological-rivsum-data
(   58) !!    i.e.: rivsum = 10. + 7.*exp(-(float(nday_year)-121.5)**2/37.5**2) !!mSv
(   59)       call rivsum_val(kt,rivsum)
(   60) 
(   61) !!(2) Convert to vel in m/s based on hardwired entry cell configuration
(   62)       da_rivsum = 0.0
(   63)       do jj = mj0(153),mj1(155)      
(   64)          do ji = mi0(14), mi1(14)       
(   65)             da_rivsum = da_rivsum + e2u(ji,jj)
(   66)          enddo
(   67)       enddo
(   68)       call mpp_sum(da_rivsum)     !!to be safe, broadcast this to all processors
(   69)       dz_rivsum = e3t(1)+e3t(2)  !!NB: inaccurate if dz=dz(i,j,k)
(   70)       vel_rivsum = (rivsum*1000.)/(da_rivsum*dz_rivsum) 
(   71) !!DBG: output a value 
(   72)       if(lwp .AND. mod(kt-nit000,int(rday/rdt)) == 0) write(2900+narea,*)kt, nyear, nday_year, rivsum, vel_rivsum
(   73) 
(   74)       do jj = mj0(153),mj1(155)      
(   75)          do ji = mi0(14), mi1(14)       !!i=15 is t-point ===> set u(14) 
(   76)             do jk = 1, 2
(   77)                un(ji,jj,jk) = ramp * vel_rivsum
(   78)                tn(ji,jj,jk) = tn(ji+1,jj,jk)  !no flux for T
(   79)                tb(ji,jj,jk) = tb(ji+1,jj,jk)  !no flux for T
(   80)                sn(ji,jj,jk) = 1.0   ! 1.0 is usual val; 5.0 ; 0.0 (JC recommendation)
(   81)                sb(ji,jj,jk) = 1.0   ! 1.0
(   82)             enddo
(   83)          enddo
(   84)       enddo
(   85) 
(   86) 
(   87)    END SUBROUTINE river_SLE_01






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 7

(   88) 
(   89) 
(   90) !!2nd routine; called by dynspg_ts
(   91)    SUBROUTINE river_SLE_02(zua,zub,zwx)
(   92) 
(   93)       !! * Arguments
(   94)      REAL(wp), DIMENSION(jpi,jpj) ::  &
(   95)           zua, zub, zwx
(   96) 
(   97)           !! * Local declarations
(   98)      INTEGER  ::   ji, jj, jk              ! dummy loop indices
(   99) 
(  100) !!DB: Adjust the above quantities to reflect runoff at head of SLE
(  101)       do jj = mj0(153),mj1(155)
(  102)          do ji = mi0(14), mi1(14)       !!i=15 is t-point ===> set u(14) 
(  103)             do jk = 1, 2
(  104)                zua(ji,jj) = zua(ji,jj) + e3t(jk) * ramp * vel_rivsum
(  105)                zub(ji,jj) = zub(ji,jj) + e3t(jk) * ramp * vel_rivsum
(  106)                zwx(ji,jj) = zwx(ji,jj) + e2u(ji,jj) * e3t(jk)  * vel_rivsum
(  107)             enddo
(  108)          enddo
(  109)       enddo
(  110) 
(  111)    END SUBROUTINE river_SLE_02
(  112) 
(  113) 
(  114) !!3rd routine; called by dynspg_ts
(  115)    SUBROUTINE river_SLE_03(ua_e)
(  116)       !! * Arguments
(  117)      REAL(wp), DIMENSION(jpi,jpj) ::  &
(  118)           ua_e
(  119) 
(  120)           !! * Local declarations
(  121)      INTEGER  ::   ji, jj, jk              ! dummy loop indices
(  122) 
(  123)       do jj = mj0(153),mj1(155)
(  124)          do ji = mi0(14), mi1(14)       !!i=15 is t-point ===> set u(14) 
(  125)             ua_e(ji,jj) = ramp * vel_rivsum * (e3t(1)+e3t(2))
(  126)          enddo
(  127)       enddo
(  128) 
(  129)    END SUBROUTINE river_SLE_03
(  130) 
(  131) 
(  132) !!DB
(  133) !!Return rivsum value
(  134) !!Routine looks for datafile: rivsum.dat which is of the form
(  135) !! year month day real_time rivsum(m3/s)
(  136) !! (NB: 1914...2007 = 94 years x 12 months = 1128 values)
(  137) !!and extracts the value based on nyear nmonth.
(  138) !!TO DO: interpolate monthly values to nday_year
(  139) !!If the file does not exist (or key_CORE_ANNUAL = off), then it
(  140) !!returns the climatological fit
(  141)    SUBROUTINE rivsum_val(kt, rivsum )
(  142)      
(  143)      !! * Arguments
(  144)      INTEGER, INTENT( in ) ::  kt
(  145)      REAL(wp), INTENT( out ) ::  rivsum






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 8

(  146) 
(  147) 
(  148) 
(  149) # 156
(  156) !!(1) Crude fit to 1914-2007 monthly rivsum data
(  157) !!    Guarantees a value for rivsum
(  158)      rivsum = 10. + 7.*exp(-(float(nday_year)-121.5)**2/37.5**2) !!mSv
(  159) 
(  160) # 187
(  187)      
(  188)    END SUBROUTINE rivsum_val
(  189) 
(  190) 
(  191)    !!======================================================================
(  192) END MODULE rivers
(  193) 
PGF90-W-0155-The type of FLOAT is now double precision with -r8  (rivers.F90: 158)
