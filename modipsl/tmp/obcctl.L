


PGF90 (Version     10.5)          09/23/2010  09:49:38      page 1

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

Filename: obcctl.F90

(    1) !!DB
(    2) !!2008.05 -- Added tmask_i as in (latest) obcvol to account for inter CPU
(    3) !!overlap problems. This fixed a volume conservation problem which affected
(    4) !!my results. 
(    5) !! **(NB: should check shape factors e2t etc in obc transport calcs) 
(    6) !!I have added a lot of DBG code. Instead of deleting it in this version, I
(    7) !!just comment it out, as it is useful to retain for future debugging
(    8) Module obcctl
(    9) #ifdef	key_obc
(   10)   !!=================================================================================
(   11)   !!                       ***  MODULE  obcctl  ***
(   12)   !! Ocean dynamics:    Balancing  the water mass in the domain to avoid potential
(   13)   !!                    unbalance due to OBC 
(   14)   !!=================================================================================
(   15) 
(   16) 
(   17)   !!----------------------------------------------------------------------------------
(   18)   !! * Modules used
(   19)   USE oce             ! ocean dynamics and tracers
(   20)   USE dom_oce         ! ocean space and time domain
(   21)   USE dynspg_oce      ! surface pressure gradient variables
(   22)   USE obc_oce         ! ocean open boundary conditions
(   23)   USE phycst          ! physical constants
(   24)   USE lib_mpp
(   25)   USE in_out_manager  ! I/O logical units
(   26) 
(   27)   IMPLICIT NONE
(   28)   PRIVATE
(   29) # include "domzgr_substitute.h90"
(   30)   !! * Accessibility
(   31)   PUBLIC obc_ctl  ! 
(   32) 
(   33)   !!---------------------------------------------------------------------------------
(   34)   !!  OPA XX , BIO 
(   35)   !!---------------------------------------------------------------------------------
(   36) 
(   37) CONTAINS
(   38) 
(   39)   SUBROUTINE obc_ctl(kt)
(   40)      !!------------------------------------------------------------------------------
(   41)      !!                      SUBROUTINE obc_ctl
(   42)      !!                     **********************
(   43)      !! ** Purpose :
(   44)      !!
(   45)      !! History : 
(   46)      !! xxx !  07-07  (Z. Wang) original






PGF90 (Version     10.5)          09/23/2010  09:49:38      page 2

(   47)      !! xxx !  07-27  DB modified
(   48)      !! DB: 
(   49)      !! (1) Adjusted sum_obc to mirror divergence calc, i.e. div = (east-west)+(north-south)
(   50)      !! (2) Changed adjustment to vels to reflect (1) 
(   51)      !!------------------------------------------------------------------------------
(   52) 
(   53) !!DB
(   54)     INTEGER, INTENT( in ) ::   kt          ! ocean time-step index
(   55)     REAL(wp):: transE,transW,transS,transN,areaE,areaW,areaS,areaN
(   56) 
(   57)     REAL(wp):: sum_obc, area, resid  ! temporary var.
(   58)     INTEGER :: ji,jj,jk,ij,ii
(   59)     COMPLEX(wp) :: cobcd,carea
(   60)     INTEGER :: itype
(   61)     cobcd = cmplx(0.e0, 0.e0, wp)
(   62)     carea = cmplx(0.e0, 0.e0, wp)
(   63) 
(   64)     sum_obc=0.0
(   65)     transE=0.0;    transW=0.0;    transS=0.0;    transN=0.0;
(   66)     areaE=0.0;     areaW=0.0;     areaS=0.0;     areaN=0.0; 
(   67) 
(   68)     IF( lp_obc_east  ) THEN 
(   69)        DO ji = nie0, nie1 
(   70)           DO  jk = 1, jpkm1
(   71)              DO jj = nje0p1, nje1m1
(   72)                 ij = jj -1 + njmpp
(   73)                 CALL DDPDD( cmplx(ufoe(jj,jk)*e2t(ji,jj)*fse3t(ji,jj,jk)*uemsk(jj,jk)*tmask_i(ji+1,jj) ,0.e0, wp), cobcd, 1,
(   74)                 CALL DDPDD( cmplx(e2t(ji,jj)*fse3t(ji,jj,jk)*uemsk(jj,jk)*tmask_i(ji+1,jj) ,0.e0, wp), carea, 1, itype) 
(   75)              END DO
(   76)           END DO
(   77)        END DO
(   78) 
(   79)     ENDIF
(   80) 
(   81)     IF( lp_obc_west  )  THEN 
(   82)        DO ji=niw0,niw1
(   83)           DO jk = 1, jpkm1
(   84)              DO jj = njw0p1, njw1m1
(   85)                 ij = jj -1 + njmpp
(   86)                 CALL DDPDD( cmplx(-ufow(jj,jk)*e2t(ji,jj)*fse3t(ji,jj,jk)*uwmsk(jj,jk)*tmask_i(ji,jj) ,0.e0, wp), cobcd, 1, 
(   87)                 CALL DDPDD( cmplx(e2t(ji,jj)*fse3t(ji,jj,jk)*uwmsk(jj,jk)*tmask_i(ji,jj),0.e0, wp), carea, 1, itype) 
(   88)              END DO
(   89)           END DO
(   90)        END DO
(   91) 
(   92)     ENDIF
(   93)     
(   94)     IF( lp_obc_north )  THEN 
(   95)        DO jj = njn0,njn1
(   96)           DO jk = 1, jpkm1
(   97)              DO ji = nin0p1, nin1m1
(   98)                 ii = ji -1 + nimpp
(   99)                 CALL DDPDD( cmplx(vfon(ji,jk)*e2t(ji,jj)*fse3t(ji,jj,jk)*vnmsk(ji,jk) ,0.e0, wp), cobcd, 1, itype)
(  100)                 CALL DDPDD( cmplx(e2t(ji,jj)*fse3t(ji,jj,jk)*vnmsk(ji,jk) ,0.e0, wp), carea, 1, itype)              
(  101)              END DO
(  102)           END DO
(  103)        END DO
(  104)        






PGF90 (Version     10.5)          09/23/2010  09:49:38      page 3

(  105)     ENDIF
(  106)     
(  107)     IF( lp_obc_south )  THEN 
(  108)        DO  jj = njs0,njs1
(  109)           DO jk = 1, jpkm1
(  110)              DO ji = nis0p1, nis1m1
(  111)                 ii = ji -1 + nimpp
(  112)                 CALL DDPDD( cmplx(-vfos(ji,jk)*e2t(ji,jj)*fse3t(ji,jj,jk)*vsmsk(ji,jk)* tmask_i(ji,jj) ,0.e0, wp), cobcd, 1,
(  113)                 CALL DDPDD( cmplx(e2t(ji,jj)*fse3t(ji,jj,jk)*vsmsk(ji,jk)* tmask_i(ji,jj),0.e0, wp), carea, 1, itype)       
(  114) 
(  115)              END DO
(  116)           END DO
(  117)        END DO
(  118) 
(  119)     ENDIF
(  120) 
(  121) 
(  122)     IF(lk_mpp) call  mpp_sum(cobcd)
(  123)     IF(lk_mpp) call  mpp_sum(carea)
(  124)     
(  125)     sum_obc=real(cobcd)
(  126)     area=real(carea)
(  127)     
(  128)     resid=sum_obc/area
(  129) 
(  130) !!DBG: Output once-per-day
(  131)     if(lwp .AND. mod(kt-nit000,int(rday/rdt)) == 0) then
(  132)        write(numout2,'(A35,2x,i10,2x,3(e12.6,1x))')  &
(  133)             'DB -- obc_vol_ctl: kt, residual: ',kt,resid
(  134)     endif
(  135)     
(  136) !!DB: changed adjustments to reflect div calc
(  137)     IF( lp_obc_east  ) THEN 
(  138)        DO ji = nie0, nie1 
(  139)           DO  jk = 1, jpkm1
(  140)              DO jj = nje0p1, nje1m1
(  141)                 ij = jj -1 + njmpp
(  142)                  ufoe(jj,jk)=(ufoe(jj,jk)-resid)*uemsk(jj,jk) 
(  143)              END DO
(  144)           END DO
(  145) 
(  146)        END DO
(  147) 
(  148)     ENDIF
(  149) 
(  150)     IF( lp_obc_west  )  THEN 
(  151)        DO ji=niw0,niw1
(  152)           DO jk = 1, jpkm1
(  153)              DO jj = njw0p1, njw1m1
(  154)                 ij = jj -1 + njmpp
(  155)                  ufow(jj,jk)=(ufow(jj,jk)+resid)*uwmsk(jj,jk)
(  156)              END DO
(  157)           END DO
(  158) 
(  159)        END DO
(  160)     ENDIF
(  161) 
(  162)     IF( lp_obc_north )  THEN 






PGF90 (Version     10.5)          09/23/2010  09:49:38      page 4

(  163)        DO jj = njn0,njn1
(  164)           DO jk = 1, jpkm1
(  165)              DO ji = nin0p1, nin1m1
(  166)                 ii = ji -1 + nimpp
(  167)                  vfon(ji,jk)=(vfon(ji,jk)-resid)*vnmsk(ji,jk)
(  168)              END DO
(  169)           END DO
(  170)        END DO
(  171) 
(  172)     ENDIF
(  173) 
(  174)     IF( lp_obc_south )  THEN 
(  175)        DO  jj = njs0,njs1
(  176)           DO jk = 1, jpkm1
(  177)              DO ji = nis0p1, nis1m1
(  178)                 ii = ji -1 + nimpp
(  179)                  vfos(ji,jk)=(vfos(ji,jk)+resid)*vsmsk(ji,jk)
(  180)              END DO
(  181)           END DO
(  182) 
(  183)        END DO
(  184)     ENDIF
(  185) 
(  186) 
(  187)      
(  188)    END SUBROUTINE obc_ctl
(  189)    
(  190)    
(  191) #else
(  192)    !!=================================================================================
(  193)    !!                       ***  MODULE  obcfla  ***
(  194)    !! Ocean dynamics:   Flather's algorithm at open boundaries for the time-splitting
(  195)    !!=================================================================================
(  196)  CONTAINS
(  197)    
(  198)    SUBROUTINE obc_ctl
(  199)    END SUBROUTINE obc_ctl
(  200) #endif
(  201)    
(  202)  END MODULE obcctl
(  203)  























PGF90 (Version     10.5)          09/23/2010  09:49:38      page 5

(    1) # 1 "obcctl.F90"
(    1) !!DB
(    2) !!2008.05 -- Added tmask_i as in (latest) obcvol to account for inter CPU
(    3) !!overlap problems. This fixed a volume conservation problem which affected
(    4) !!my results. 
(    5) !! **(NB: should check shape factors e2t etc in obc transport calcs) 
(    6) !!I have added a lot of DBG code. Instead of deleting it in this version, I
(    7) !!just comment it out, as it is useful to retain for future debugging
(    8) Module obcctl
(    9) # 10
(   10)   !!=================================================================================
(   11)   !!                       ***  MODULE  obcctl  ***
(   12)   !! Ocean dynamics:    Balancing  the water mass in the domain to avoid potential
(   13)   !!                    unbalance due to OBC 
(   14)   !!=================================================================================
(   15) 
(   16) 
(   17)   !!----------------------------------------------------------------------------------
(   18)   !! * Modules used
(   19)   USE oce             ! ocean dynamics and tracers
(   20)   USE dom_oce         ! ocean space and time domain
(   21)   USE dynspg_oce      ! surface pressure gradient variables
(   22)   USE obc_oce         ! ocean open boundary conditions
(   23)   USE phycst          ! physical constants
(   24)   USE lib_mpp
(   25)   USE in_out_manager  ! I/O logical units
(   26) 
(   27)   IMPLICIT NONE
(   28)   PRIVATE
(   29) # 1 "./domzgr_substitute.h90"
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
(   48) # 30 "obcctl.F90"
(   30) # 30
(   30)   !! * Accessibility
(   31)   PUBLIC obc_ctl  ! 
(   32) 
(   33)   !!---------------------------------------------------------------------------------
(   34)   !!  OPA XX , BIO 
(   35)   !!---------------------------------------------------------------------------------
(   36) 
(   37) CONTAINS
(   38) 
(   39)   SUBROUTINE obc_ctl(kt)
(   40)      !!------------------------------------------------------------------------------
(   41)      !!                      SUBROUTINE obc_ctl






PGF90 (Version     10.5)          09/23/2010  09:49:38      page 6

(   42)      !!                     **********************
(   43)      !! ** Purpose :
(   44)      !!
(   45)      !! History : 
(   46)      !! xxx !  07-07  (Z. Wang) original
(   47)      !! xxx !  07-27  DB modified
(   48)      !! DB: 
(   49)      !! (1) Adjusted sum_obc to mirror divergence calc, i.e. div = (east-west)+(north-south)
(   50)      !! (2) Changed adjustment to vels to reflect (1) 
(   51)      !!------------------------------------------------------------------------------
(   52) 
(   53) !!DB
(   54)     INTEGER, INTENT( in ) ::   kt          ! ocean time-step index
(   55)     REAL(wp):: transE,transW,transS,transN,areaE,areaW,areaS,areaN
(   56) 
(   57)     REAL(wp):: sum_obc, area, resid  ! temporary var.
(   58)     INTEGER :: ji,jj,jk,ij,ii
(   59)     COMPLEX(wp) :: cobcd,carea
(   60)     INTEGER :: itype
(   61)     cobcd = cmplx(0.e0, 0.e0, wp)
(   62)     carea = cmplx(0.e0, 0.e0, wp)
(   63) 
(   64)     sum_obc=0.0
(   65)     transE=0.0;    transW=0.0;    transS=0.0;    transN=0.0;
(   66)     areaE=0.0;     areaW=0.0;     areaS=0.0;     areaN=0.0; 
(   67) 
(   68)     IF( lp_obc_east  ) THEN 
(   69)        DO ji = nie0, nie1 
(   70)           DO  jk = 1, jpkm1
(   71)              DO jj = nje0p1, nje1m1
(   72)                 ij = jj -1 + njmpp
(   73)                 CALL DDPDD( cmplx(ufoe(jj,jk)*e2t(ji,jj)*e3t(jk)*uemsk(jj,jk)*tmask_i(ji+1,jj) ,0.e0, wp), cobcd, 1, itype)
(   74)                 CALL DDPDD( cmplx(e2t(ji,jj)*e3t(jk)*uemsk(jj,jk)*tmask_i(ji+1,jj) ,0.e0, wp), carea, 1, itype) 
(   75)              END DO
(   76)           END DO
(   77)        END DO
(   78) 
(   79)     ENDIF
(   80) 
(   81)     IF( lp_obc_west  )  THEN 
(   82)        DO ji=niw0,niw1
(   83)           DO jk = 1, jpkm1
(   84)              DO jj = njw0p1, njw1m1
(   85)                 ij = jj -1 + njmpp
(   86)                 CALL DDPDD( cmplx(-ufow(jj,jk)*e2t(ji,jj)*e3t(jk)*uwmsk(jj,jk)*tmask_i(ji,jj) ,0.e0, wp), cobcd, 1, itype)
(   87)                 CALL DDPDD( cmplx(e2t(ji,jj)*e3t(jk)*uwmsk(jj,jk)*tmask_i(ji,jj),0.e0, wp), carea, 1, itype) 
(   88)              END DO
(   89)           END DO
(   90)        END DO
(   91) 
(   92)     ENDIF
(   93)     
(   94)     IF( lp_obc_north )  THEN 
(   95)        DO jj = njn0,njn1
(   96)           DO jk = 1, jpkm1
(   97)              DO ji = nin0p1, nin1m1
(   98)                 ii = ji -1 + nimpp
(   99)                 CALL DDPDD( cmplx(vfon(ji,jk)*e2t(ji,jj)*e3t(jk)*vnmsk(ji,jk) ,0.e0, wp), cobcd, 1, itype)






PGF90 (Version     10.5)          09/23/2010  09:49:38      page 7

(  100)                 CALL DDPDD( cmplx(e2t(ji,jj)*e3t(jk)*vnmsk(ji,jk) ,0.e0, wp), carea, 1, itype)              
(  101)              END DO
(  102)           END DO
(  103)        END DO
(  104)        
(  105)     ENDIF
(  106)     
(  107)     IF( lp_obc_south )  THEN 
(  108)        DO  jj = njs0,njs1
(  109)           DO jk = 1, jpkm1
(  110)              DO ji = nis0p1, nis1m1
(  111)                 ii = ji -1 + nimpp
(  112)                 CALL DDPDD( cmplx(-vfos(ji,jk)*e2t(ji,jj)*e3t(jk)*vsmsk(ji,jk)* tmask_i(ji,jj) ,0.e0, wp), cobcd, 1, itype)
(  113)                 CALL DDPDD( cmplx(e2t(ji,jj)*e3t(jk)*vsmsk(ji,jk)* tmask_i(ji,jj),0.e0, wp), carea, 1, itype)              
(  114) 
(  115)              END DO
(  116)           END DO
(  117)        END DO
(  118) 
(  119)     ENDIF
(  120) 
(  121) 
(  122)     IF(lk_mpp) call  mpp_sum(cobcd)
(  123)     IF(lk_mpp) call  mpp_sum(carea)
(  124)     
(  125)     sum_obc=real(cobcd)
(  126)     area=real(carea)
(  127)     
(  128)     resid=sum_obc/area
(  129) 
(  130) !!DBG: Output once-per-day
(  131)     if(lwp .AND. mod(kt-nit000,int(rday/rdt)) == 0) then
(  132)        write(numout2,'(A35,2x,i10,2x,3(e12.6,1x))')  &
(  133)             'DB -- obc_vol_ctl: kt, residual: ',kt,resid
(  134)     endif
(  135)     
(  136) !!DB: changed adjustments to reflect div calc
(  137)     IF( lp_obc_east  ) THEN 
(  138)        DO ji = nie0, nie1 
(  139)           DO  jk = 1, jpkm1
(  140)              DO jj = nje0p1, nje1m1
(  141)                 ij = jj -1 + njmpp
(  142)                  ufoe(jj,jk)=(ufoe(jj,jk)-resid)*uemsk(jj,jk) 
(  143)              END DO
(  144)           END DO
(  145) 
(  146)        END DO
(  147) 
(  148)     ENDIF
(  149) 
(  150)     IF( lp_obc_west  )  THEN 
(  151)        DO ji=niw0,niw1
(  152)           DO jk = 1, jpkm1
(  153)              DO jj = njw0p1, njw1m1
(  154)                 ij = jj -1 + njmpp
(  155)                  ufow(jj,jk)=(ufow(jj,jk)+resid)*uwmsk(jj,jk)
(  156)              END DO
(  157)           END DO






PGF90 (Version     10.5)          09/23/2010  09:49:38      page 8

(  158) 
(  159)        END DO
(  160)     ENDIF
(  161) 
(  162)     IF( lp_obc_north )  THEN 
(  163)        DO jj = njn0,njn1
(  164)           DO jk = 1, jpkm1
(  165)              DO ji = nin0p1, nin1m1
(  166)                 ii = ji -1 + nimpp
(  167)                  vfon(ji,jk)=(vfon(ji,jk)-resid)*vnmsk(ji,jk)
(  168)              END DO
(  169)           END DO
(  170)        END DO
(  171) 
(  172)     ENDIF
(  173) 
(  174)     IF( lp_obc_south )  THEN 
(  175)        DO  jj = njs0,njs1
(  176)           DO jk = 1, jpkm1
(  177)              DO ji = nis0p1, nis1m1
(  178)                 ii = ji -1 + nimpp
(  179)                  vfos(ji,jk)=(vfos(ji,jk)+resid)*vsmsk(ji,jk)
(  180)              END DO
(  181)           END DO
(  182) 
(  183)        END DO
(  184)     ENDIF
(  185) 
(  186) 
(  187)      
(  188)    END SUBROUTINE obc_ctl
(  189)    
(  190)    
(  191) # 201
(  201)    
(  202)  END MODULE obcctl
(  203)  
