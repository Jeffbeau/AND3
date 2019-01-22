


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

Filename: obcfla.F90

(    1) MODULE obcfla
(    2) #if defined key_obc && defined key_dynspg_ts
(    3)   !!=================================================================================
(    4)   !!                       ***  MODULE  obcfla  ***
(    5)   !! Ocean dynamics:   Flather's algorithm at open boundaries for the time-splitting
(    6)   !!=================================================================================
(    7) 
(    8)   !!---------------------------------------------------------------------------------
(    9)   !!   obc_fla_ts        : call the subroutine for each open boundary
(   10)   !!   obc_fla_ts_east   : Flather on the east  open boundary velocities & ssh
(   11)   !!   obc_fla_ts_west   : Flather on the west  open boundary velocities & ssh
(   12)   !!   obc_fla_ts_north  : Flather on the north open boundary velocities & ssh
(   13)   !!   obc_fla_ts_south  : Flather on the south open boundary velocities & ssh
(   14)   !!----------------------------------------------------------------------------------
(   15) 
(   16)   !!----------------------------------------------------------------------------------
(   17)   !! * Modules used
(   18)   USE oce             ! ocean dynamics and tracers
(   19)   USE dom_oce         ! ocean space and time domain
(   20)   USE dynspg_oce      ! surface pressure gradient variables
(   21)   USE phycst          ! physical constants
(   22)   USE obc_oce         ! ocean open boundary conditions
(   23)   USE obcdta          ! ocean open boundary conditions: climatology
(   24) 
(   25)   IMPLICIT NONE
(   26)   PRIVATE
(   27) 
(   28)   !! * Accessibility
(   29)   PUBLIC obc_fla_ts  ! routine called in dynspg_ts (free surface time splitting case)
(   30) 
(   31)   !!---------------------------------------------------------------------------------
(   32)   !!  OPA 9.0 , LOCEAN-IPSL (2005)
(   33)   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obcfla.F90,v 1.2 2006/01/03 15:04:15 opalod Exp $
(   34)   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   35)   !!---------------------------------------------------------------------------------
(   36) 
(   37) CONTAINS
(   38) 
(   39)   SUBROUTINE obc_fla_ts
(   40)      !!------------------------------------------------------------------------------
(   41)      !!                      SUBROUTINE obc_fla_ts
(   42)      !!                     **********************
(   43)      !! ** Purpose :
(   44)      !!      Apply Flather's algorithm at open boundaries for the time-splitting
(   45)      !!      free surface case (barotropic variables)
(   46)      !!






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 2

(   47)      !!      This routine is called in dynspg_ts.F90 routine
(   48)      !!
(   49)      !!      The logical variable lp_obc_east, and/or lp_obc_west, and/or lp_obc_north,
(   50)      !!      and/or lp_obc_south allow the user to determine which boundary is an
(   51)      !!      open one (must be done in the obc_par.F90 file).
(   52)      !!
(   53)      !! ** Reference :
(   54)      !!         Flather, R. A., 1976, Mem. Soc. R. Sci. Liege, Ser. 6, 10, 141-164
(   55)      !!
(   56)      !! History :
(   57)      !!   9.0  !  05-12  (V. Garnier) original
(   58)      !!------------------------------------------------------------------------------
(   59) 
(   60)      IF( lp_obc_east  )   CALL obc_fla_ts_east
(   61)      IF( lp_obc_west  )   CALL obc_fla_ts_west
(   62)      IF( lp_obc_north )   CALL obc_fla_ts_north
(   63)      IF( lp_obc_south )   CALL obc_fla_ts_south
(   64) 
(   65)   END SUBROUTINE obc_fla_ts
(   66) 
(   67) 
(   68)   SUBROUTINE obc_fla_ts_east
(   69)      !!------------------------------------------------------------------------------
(   70)      !!                  ***  SUBROUTINE obc_fla_ts_east  ***
(   71)      !!
(   72)      !! ** Purpose :
(   73)      !!      Apply Flather's algorithm on east OBC velocities ua, va
(   74)      !!      Fix sea surface height (sshn_e) on east open boundary
(   75)      !!
(   76)      !!  History :
(   77)      !!   9.0  !  05-12  (V. Garnier) original
(   78)      !!------------------------------------------------------------------------------
(   79)      !! * Local declaration
(   80)      INTEGER ::   ji, jj, jk ! dummy loop indices
(   81)      !!------------------------------------------------------------------------------
(   82) 
(   83)      DO ji = nie0, nie1
(   84) !ylu         DO jk = 1, jpkm1
(   85)            DO jj = 1, jpj
(   86)               ua_e(ji,jj) = (  ubtfoe(jj) + sqrt( grav*hu(ji,jj) )           &
(   87)                  &            * (  sshn_e(ji,jj)  &
(   88)                  &            - sshfoe(jj) )  ) * uemsk(jj,1)
(   89) !ylu                  &            * ( ( sshn_e(ji,jj) + sshn_e(ji+1,jj) ) * 0.5  &
(   90) !ylu                  &            - sshfoe(jj) )  ) * uemsk(jj,jk)
(   91) !ylu               va_e(ji,jj) = vbtfoe(jj) * uemsk(jj,jk)
(   92)            END DO
(   93) !ylu         END DO
(   94)      END DO
(   95) 
(   96)      DO ji = nie0p1, nie1p1
(   97)         DO jj = 1, jpj
(   98)            sshfoe_b(ji,jj) = sshfoe_b(ji,jj) + sqrt( grav*hur(ji,jj) )     &
(   99)               &             * ( ( sshn_e(ji,jj) + sshn_e(ji+1,jj) ) * 0.5  &
(  100)               &                 - sshfoe(jj) ) * uemsk(jj,1)
(  101)            ssha_e(ji,jj) = ssha_e(ji,jj) * ( 1. - temsk(jj,1) ) &
(  102)               &            + temsk(jj,1) * sshfoe(jj)
(  103)         END DO
(  104)      END DO






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 3

(  105) 
(  106)   END SUBROUTINE obc_fla_ts_east
(  107) 
(  108) 
(  109)   SUBROUTINE obc_fla_ts_west
(  110)      !!------------------------------------------------------------------------------
(  111)      !!                  ***  SUBROUTINE obc_fla_ts_west  ***
(  112)      !!
(  113)      !! ** Purpose :
(  114)      !!      Apply Flather's algorithm on west OBC velocities ua, va
(  115)      !!      Fix sea surface height (sshn_e) on west open boundary
(  116)      !!
(  117)      !!  History :
(  118)      !!   9.0  !  05-12  (V. Garnier) original
(  119)      !!------------------------------------------------------------------------------
(  120)      !! * Local declaration
(  121)      INTEGER ::   ji, jj, jk ! dummy loop indices
(  122)      !!------------------------------------------------------------------------------
(  123) 
(  124)      DO ji = niw0, niw1
(  125) !ylu         DO jk = 1, jpkm1
(  126)            DO jj = 1, jpj
(  127)               ua_e(ji,jj) = ( ubtfow(jj) - sqrt( grav * hu(ji,jj) )          &
(  128)                  &            * (  sshn_e(ji+1,jj)  &
(  129)                  &                - sshfow(jj) ) ) * uwmsk(jj,1)
(  130) !!DBG- severe
(  131) !              ua_e(ji,jj) = 0.0
(  132) 
(  133) !ylu                  &            * ( ( sshn_e(ji,jj) + sshn_e(ji+1,jj) ) * 0.5  &
(  134) !ylu                  &                - sshfow(jj) ) ) * uwmsk(jj,jk)
(  135) !ylu               va_e(ji,jj) = vbtfow(jj) * uwmsk(jj,jk)
(  136)            END DO
(  137) !ylu         END DO
(  138) 
(  139)         DO jj = 1, jpj
(  140)            sshfow_b(ji,jj) = sshfow_b(ji,jj) - sqrt( grav * hur(ji,jj) )     &
(  141)                              * ( ( sshn_e(ji,jj) + sshn_e(ji+1,jj) ) * 0.5   &
(  142)                                 - sshfow(jj) ) * uwmsk(jj,1)
(  143)            ssha_e(ji,jj) = ssha_e(ji,jj) * ( 1. - twmsk(jj,1) ) &
(  144)               &            + twmsk(jj,1)*sshfow(jj)
(  145)         END DO
(  146)      END DO
(  147) 
(  148) 
(  149)   END SUBROUTINE obc_fla_ts_west
(  150) 
(  151)   SUBROUTINE obc_fla_ts_north
(  152)      !!------------------------------------------------------------------------------
(  153)      !!                     SUBROUTINE obc_fla_ts_north
(  154)      !!                    *************************
(  155)      !! ** Purpose :
(  156)      !!      Apply Flather's algorithm on north OBC velocities ua, va
(  157)      !!      Fix sea surface height (sshn_e) on north open boundary
(  158)      !!
(  159)      !!  History :
(  160)      !!   9.0  !  05-12  (V. Garnier) original
(  161)      !!------------------------------------------------------------------------------
(  162)      !! * Local declaration






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 4

(  163)      INTEGER ::   ji, jj, jk ! dummy loop indices
(  164)      !!------------------------------------------------------------------------------
(  165) 
(  166)      DO jj = njn0, njn1
(  167) !ylu         DO jk = 1, jpkm1
(  168)            DO ji = 1, jpi
(  169)               va_e(ji,jj) = ( vbtfon(ji) + sqrt( grav * hv(ji,jj) )           &
(  170)                  &            * (  sshn_e(ji,jj)    &
(  171)                  &                - sshfon(ji) ) ) * vnmsk(ji,1) 
(  172) !ylu                  &            * ( ( sshn_e(ji,jj) + sshn_e(ji,jj+1) ) * 0.5   &
(  173) !ylu                  &                - sshfon(ji) ) ) * vnmsk(ji,jk)
(  174) !ylu               ua_e(ji,jj) = ubtfon(ji) * vnmsk(ji,jk)
(  175)            END DO
(  176) !ylu         END DO
(  177)      END DO
(  178)      DO jj = njn0p1, njn1p1
(  179)         DO ji = 1, jpi
(  180)            sshfon_b(ji,jj) = sshfon_b(ji,jj) + sqrt( grav * hvr(ji,jj) )  &
(  181)               &              * ( ( sshn_e(ji,jj) + sshn_e(ji,jj+1) ) * 0.5    &
(  182)               &                  - sshfon(ji) ) * vnmsk(ji,1)
(  183)            ssha_e(ji,jj) = ssha_e(ji,jj) * ( 1. - tnmsk(ji,1) ) &
(  184)               &            + sshfon(ji) * tnmsk(ji,1)
(  185)         END DO
(  186)      END DO
(  187) 
(  188)   END SUBROUTINE obc_fla_ts_north
(  189) 
(  190)   SUBROUTINE obc_fla_ts_south
(  191)      !!------------------------------------------------------------------------------
(  192)      !!                     SUBROUTINE obc_fla_ts_south
(  193)      !!                    *************************
(  194)      !! ** Purpose :
(  195)      !!      Apply Flather's algorithm on south OBC velocities ua, va
(  196)      !!      Fix sea surface height (sshn_e) on south open boundary
(  197)      !!
(  198)      !!  History :
(  199)      !!   9.0  !  05-12  (V. Garnier) original
(  200)      !!------------------------------------------------------------------------------
(  201)      !! * Local declaration
(  202)      INTEGER ::   ji, jj, jk ! dummy loop indices
(  203) 
(  204)      !!------------------------------------------------------------------------------
(  205) 
(  206) 
(  207)      DO jj = njs0, njs1
(  208) !ylu         DO jk = 1, jpkm1
(  209)            DO ji = 1, jpi
(  210)               va_e(ji,jj) = ( vbtfos(ji) - sqrt( grav * hv(ji,jj) )            &
(  211)                  &            * ( sshn_e(ji,jj+1)   &
(  212)                  &                - sshfos(ji) ) ) * vsmsk(ji,1)
(  213) !ylu                  &            * ( ( sshn_e(ji,jj) + sshn_e(ji,jj+1) ) * 0.5    &
(  214) !ylu                  &                - sshfos(ji) ) ) * vsmsk(ji,jk)
(  215) !ylu               ua_e(ji,jj) = ubtfos(ji) * vsmsk(ji,jk)
(  216)            END DO
(  217) !ylu         END DO
(  218)         DO ji = 1, jpi
(  219)            sshfos_b(ji,jj) = sshfos_b(ji,jj) - sqrt( grav * hvr(ji,jj) )      &
(  220)               &              * ( ( sshn_e(ji,jj) + sshn_e(ji,jj+1) ) * 0.5    &






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 5

(  221)               &                  - sshfos(ji) ) * vsmsk(ji,1)
(  222)            ssha_e(ji,jj) = ssha_e(ji,jj) * (1. - tsmsk(ji,1) ) &
(  223)               &            + tsmsk(ji,1) * sshfos(ji)
(  224)         END DO
(  225)      END DO
(  226) 
(  227)   END SUBROUTINE obc_fla_ts_south
(  228) #else
(  229)   !!=================================================================================
(  230)   !!                       ***  MODULE  obcfla  ***
(  231)   !! Ocean dynamics:   Flather's algorithm at open boundaries for the time-splitting
(  232)   !!=================================================================================
(  233) CONTAINS
(  234)   USE in_out_manager
(  235)   SUBROUTINE obc_fla_ts
(  236)      if(lwp) WRITE(numout,*) 'obc_fla_ts: You should not have seen this print! error?'
(  237)   END SUBROUTINE obc_fla_ts
(  238) #endif
(  239) 
(  240) END MODULE obcfla












































PGF90 (Version     10.5)          09/23/2010  09:49:41      page 6

(    1) # 1 "obcfla.F90"
(    1) MODULE obcfla
(    2) # 3
(    3)   !!=================================================================================
(    4)   !!                       ***  MODULE  obcfla  ***
(    5)   !! Ocean dynamics:   Flather's algorithm at open boundaries for the time-splitting
(    6)   !!=================================================================================
(    7) 
(    8)   !!---------------------------------------------------------------------------------
(    9)   !!   obc_fla_ts        : call the subroutine for each open boundary
(   10)   !!   obc_fla_ts_east   : Flather on the east  open boundary velocities & ssh
(   11)   !!   obc_fla_ts_west   : Flather on the west  open boundary velocities & ssh
(   12)   !!   obc_fla_ts_north  : Flather on the north open boundary velocities & ssh
(   13)   !!   obc_fla_ts_south  : Flather on the south open boundary velocities & ssh
(   14)   !!----------------------------------------------------------------------------------
(   15) 
(   16)   !!----------------------------------------------------------------------------------
(   17)   !! * Modules used
(   18)   USE oce             ! ocean dynamics and tracers
(   19)   USE dom_oce         ! ocean space and time domain
(   20)   USE dynspg_oce      ! surface pressure gradient variables
(   21)   USE phycst          ! physical constants
(   22)   USE obc_oce         ! ocean open boundary conditions
(   23)   USE obcdta          ! ocean open boundary conditions: climatology
(   24) 
(   25)   IMPLICIT NONE
(   26)   PRIVATE
(   27) 
(   28)   !! * Accessibility
(   29)   PUBLIC obc_fla_ts  ! routine called in dynspg_ts (free surface time splitting case)
(   30) 
(   31)   !!---------------------------------------------------------------------------------
(   32)   !!  OPA 9.0 , LOCEAN-IPSL (2005)
(   33)   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/OBC/obcfla.F90,v 1.2 2006/01/03 15:04:15 opalod Exp $
(   34)   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
(   35)   !!---------------------------------------------------------------------------------
(   36) 
(   37) CONTAINS
(   38) 
(   39)   SUBROUTINE obc_fla_ts
(   40)      !!------------------------------------------------------------------------------
(   41)      !!                      SUBROUTINE obc_fla_ts
(   42)      !!                     **********************
(   43)      !! ** Purpose :
(   44)      !!      Apply Flather's algorithm at open boundaries for the time-splitting
(   45)      !!      free surface case (barotropic variables)
(   46)      !!
(   47)      !!      This routine is called in dynspg_ts.F90 routine
(   48)      !!
(   49)      !!      The logical variable lp_obc_east, and/or lp_obc_west, and/or lp_obc_north,
(   50)      !!      and/or lp_obc_south allow the user to determine which boundary is an
(   51)      !!      open one (must be done in the obc_par.F90 file).
(   52)      !!
(   53)      !! ** Reference :
(   54)      !!         Flather, R. A., 1976, Mem. Soc. R. Sci. Liege, Ser. 6, 10, 141-164
(   55)      !!
(   56)      !! History :
(   57)      !!   9.0  !  05-12  (V. Garnier) original






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 7

(   58)      !!------------------------------------------------------------------------------
(   59) 
(   60)      IF( lp_obc_east  )   CALL obc_fla_ts_east
(   61)      IF( lp_obc_west  )   CALL obc_fla_ts_west
(   62)      IF( lp_obc_north )   CALL obc_fla_ts_north
(   63)      IF( lp_obc_south )   CALL obc_fla_ts_south
(   64) 
(   65)   END SUBROUTINE obc_fla_ts
(   66) 
(   67) 
(   68)   SUBROUTINE obc_fla_ts_east
(   69)      !!------------------------------------------------------------------------------
(   70)      !!                  ***  SUBROUTINE obc_fla_ts_east  ***
(   71)      !!
(   72)      !! ** Purpose :
(   73)      !!      Apply Flather's algorithm on east OBC velocities ua, va
(   74)      !!      Fix sea surface height (sshn_e) on east open boundary
(   75)      !!
(   76)      !!  History :
(   77)      !!   9.0  !  05-12  (V. Garnier) original
(   78)      !!------------------------------------------------------------------------------
(   79)      !! * Local declaration
(   80)      INTEGER ::   ji, jj, jk ! dummy loop indices
(   81)      !!------------------------------------------------------------------------------
(   82) 
(   83)      DO ji = nie0, nie1
(   84) !ylu         DO jk = 1, jpkm1
(   85)            DO jj = 1, jpj
(   86)               ua_e(ji,jj) = (  ubtfoe(jj) + sqrt( grav*hu(ji,jj) )           &
(   87)                  &            * (  sshn_e(ji,jj)  &
(   88)                  &            - sshfoe(jj) )  ) * uemsk(jj,1)
(   89) !ylu                  &            * ( ( sshn_e(ji,jj) + sshn_e(ji+1,jj) ) * 0.5  &
(   90) !ylu                  &            - sshfoe(jj) )  ) * uemsk(jj,jk)
(   91) !ylu               va_e(ji,jj) = vbtfoe(jj) * uemsk(jj,jk)
(   92)            END DO
(   93) !ylu         END DO
(   94)      END DO
(   95) 
(   96)      DO ji = nie0p1, nie1p1
(   97)         DO jj = 1, jpj
(   98)            sshfoe_b(ji,jj) = sshfoe_b(ji,jj) + sqrt( grav*hur(ji,jj) )     &
(   99)               &             * ( ( sshn_e(ji,jj) + sshn_e(ji+1,jj) ) * 0.5  &
(  100)               &                 - sshfoe(jj) ) * uemsk(jj,1)
(  101)            ssha_e(ji,jj) = ssha_e(ji,jj) * ( 1. - temsk(jj,1) ) &
(  102)               &            + temsk(jj,1) * sshfoe(jj)
(  103)         END DO
(  104)      END DO
(  105) 
(  106)   END SUBROUTINE obc_fla_ts_east
(  107) 
(  108) 
(  109)   SUBROUTINE obc_fla_ts_west
(  110)      !!------------------------------------------------------------------------------
(  111)      !!                  ***  SUBROUTINE obc_fla_ts_west  ***
(  112)      !!
(  113)      !! ** Purpose :
(  114)      !!      Apply Flather's algorithm on west OBC velocities ua, va
(  115)      !!      Fix sea surface height (sshn_e) on west open boundary






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 8

(  116)      !!
(  117)      !!  History :
(  118)      !!   9.0  !  05-12  (V. Garnier) original
(  119)      !!------------------------------------------------------------------------------
(  120)      !! * Local declaration
(  121)      INTEGER ::   ji, jj, jk ! dummy loop indices
(  122)      !!------------------------------------------------------------------------------
(  123) 
(  124)      DO ji = niw0, niw1
(  125) !ylu         DO jk = 1, jpkm1
(  126)            DO jj = 1, jpj
(  127)               ua_e(ji,jj) = ( ubtfow(jj) - sqrt( grav * hu(ji,jj) )          &
(  128)                  &            * (  sshn_e(ji+1,jj)  &
(  129)                  &                - sshfow(jj) ) ) * uwmsk(jj,1)
(  130) !!DBG- severe
(  131) !              ua_e(ji,jj) = 0.0
(  132) 
(  133) !ylu                  &            * ( ( sshn_e(ji,jj) + sshn_e(ji+1,jj) ) * 0.5  &
(  134) !ylu                  &                - sshfow(jj) ) ) * uwmsk(jj,jk)
(  135) !ylu               va_e(ji,jj) = vbtfow(jj) * uwmsk(jj,jk)
(  136)            END DO
(  137) !ylu         END DO
(  138) 
(  139)         DO jj = 1, jpj
(  140)            sshfow_b(ji,jj) = sshfow_b(ji,jj) - sqrt( grav * hur(ji,jj) )     &
(  141)                              * ( ( sshn_e(ji,jj) + sshn_e(ji+1,jj) ) * 0.5   &
(  142)                                 - sshfow(jj) ) * uwmsk(jj,1)
(  143)            ssha_e(ji,jj) = ssha_e(ji,jj) * ( 1. - twmsk(jj,1) ) &
(  144)               &            + twmsk(jj,1)*sshfow(jj)
(  145)         END DO
(  146)      END DO
(  147) 
(  148) 
(  149)   END SUBROUTINE obc_fla_ts_west
(  150) 
(  151)   SUBROUTINE obc_fla_ts_north
(  152)      !!------------------------------------------------------------------------------
(  153)      !!                     SUBROUTINE obc_fla_ts_north
(  154)      !!                    *************************
(  155)      !! ** Purpose :
(  156)      !!      Apply Flather's algorithm on north OBC velocities ua, va
(  157)      !!      Fix sea surface height (sshn_e) on north open boundary
(  158)      !!
(  159)      !!  History :
(  160)      !!   9.0  !  05-12  (V. Garnier) original
(  161)      !!------------------------------------------------------------------------------
(  162)      !! * Local declaration
(  163)      INTEGER ::   ji, jj, jk ! dummy loop indices
(  164)      !!------------------------------------------------------------------------------
(  165) 
(  166)      DO jj = njn0, njn1
(  167) !ylu         DO jk = 1, jpkm1
(  168)            DO ji = 1, jpi
(  169)               va_e(ji,jj) = ( vbtfon(ji) + sqrt( grav * hv(ji,jj) )           &
(  170)                  &            * (  sshn_e(ji,jj)    &
(  171)                  &                - sshfon(ji) ) ) * vnmsk(ji,1) 
(  172) !ylu                  &            * ( ( sshn_e(ji,jj) + sshn_e(ji,jj+1) ) * 0.5   &
(  173) !ylu                  &                - sshfon(ji) ) ) * vnmsk(ji,jk)






PGF90 (Version     10.5)          09/23/2010  09:49:41      page 9

(  174) !ylu               ua_e(ji,jj) = ubtfon(ji) * vnmsk(ji,jk)
(  175)            END DO
(  176) !ylu         END DO
(  177)      END DO
(  178)      DO jj = njn0p1, njn1p1
(  179)         DO ji = 1, jpi
(  180)            sshfon_b(ji,jj) = sshfon_b(ji,jj) + sqrt( grav * hvr(ji,jj) )  &
(  181)               &              * ( ( sshn_e(ji,jj) + sshn_e(ji,jj+1) ) * 0.5    &
(  182)               &                  - sshfon(ji) ) * vnmsk(ji,1)
(  183)            ssha_e(ji,jj) = ssha_e(ji,jj) * ( 1. - tnmsk(ji,1) ) &
(  184)               &            + sshfon(ji) * tnmsk(ji,1)
(  185)         END DO
(  186)      END DO
(  187) 
(  188)   END SUBROUTINE obc_fla_ts_north
(  189) 
(  190)   SUBROUTINE obc_fla_ts_south
(  191)      !!------------------------------------------------------------------------------
(  192)      !!                     SUBROUTINE obc_fla_ts_south
(  193)      !!                    *************************
(  194)      !! ** Purpose :
(  195)      !!      Apply Flather's algorithm on south OBC velocities ua, va
(  196)      !!      Fix sea surface height (sshn_e) on south open boundary
(  197)      !!
(  198)      !!  History :
(  199)      !!   9.0  !  05-12  (V. Garnier) original
(  200)      !!------------------------------------------------------------------------------
(  201)      !! * Local declaration
(  202)      INTEGER ::   ji, jj, jk ! dummy loop indices
(  203) 
(  204)      !!------------------------------------------------------------------------------
(  205) 
(  206) 
(  207)      DO jj = njs0, njs1
(  208) !ylu         DO jk = 1, jpkm1
(  209)            DO ji = 1, jpi
(  210)               va_e(ji,jj) = ( vbtfos(ji) - sqrt( grav * hv(ji,jj) )            &
(  211)                  &            * ( sshn_e(ji,jj+1)   &
(  212)                  &                - sshfos(ji) ) ) * vsmsk(ji,1)
(  213) !ylu                  &            * ( ( sshn_e(ji,jj) + sshn_e(ji,jj+1) ) * 0.5    &
(  214) !ylu                  &                - sshfos(ji) ) ) * vsmsk(ji,jk)
(  215) !ylu               ua_e(ji,jj) = ubtfos(ji) * vsmsk(ji,jk)
(  216)            END DO
(  217) !ylu         END DO
(  218)         DO ji = 1, jpi
(  219)            sshfos_b(ji,jj) = sshfos_b(ji,jj) - sqrt( grav * hvr(ji,jj) )      &
(  220)               &              * ( ( sshn_e(ji,jj) + sshn_e(ji,jj+1) ) * 0.5    &
(  221)               &                  - sshfos(ji) ) * vsmsk(ji,1)
(  222)            ssha_e(ji,jj) = ssha_e(ji,jj) * (1. - tsmsk(ji,1) ) &
(  223)               &            + tsmsk(ji,1) * sshfos(ji)
(  224)         END DO
(  225)      END DO
(  226) 
(  227)   END SUBROUTINE obc_fla_ts_south
(  228) 
(  229) # 240
(  240) END MODULE obcfla
