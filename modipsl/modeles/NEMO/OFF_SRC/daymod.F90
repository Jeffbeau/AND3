MODULE daymod
   !!======================================================================
   !!                       ***  MODULE  daymod  ***
   !! Ocean        :  calendar 
   !!=====================================================================

   !!----------------------------------------------------------------------
   !!   day        : calendar
   !!----------------------------------------------------------------------
   !! * Modules used
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC day        ! called by step.F90

   !! * Shared module variables
   INTEGER , PUBLIC ::   &  !:
      nyear     ,   &  !: current year
      nmonth    ,   &  !: current month
      nday      ,   &  !: current day of the month
      nday_year ,   &  !: curent day counted from jan 1st of the current year
      ndastp           !: time step date in year/month/day aammjj
   REAL(wp), PUBLIC ::   &  !:
       adatrj   ,   &  !: number of elapsed days since the begining of the run
       adatrj0         !: value of adatrj at nit000-1 (before the present run).
       !               !  it is the accumulated duration of previous runs
       !               !  that may have been run with different time steps.
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OFF_SRC/daymod.F90,v 1.1.1.1 2005/11/14 10:41:07 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE day( kt )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE day  ***
      !! 
      !! ** Purpose :   Compute the date with a day iteration IF necessary.
      !!
      !! ** Method  : - ???
      !!
      !! ** Action  : - nyear     : current year
      !!              - nmonth    : current month of the year nyear
      !!              - nday      : current day of the month nmonth
      !!              - nday_year : current day of the year nyear
      !!              - ndastp    : =nyear*10000+nmonth*100+nday
      !!              - adatrj    : date in days since the beginning of the run
      !!
      !! History :
      !!        !  94-09  (M. Pontaud M. Imbard)  Original code
      !!        !  97-03  (O. Marti)
      !!        !  97-05  (G. Madec) 
      !!        !  97-08  (M. Imbard)
      !!   9.0  !  03-09  (G. Madec)  F90 + nyear, nmonth, nday
      !!        !  04-01  (A.M. Treguier) new calculation based on adatrj
      !!----------------------------------------------------------------------      
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step indices

      !! * Local declarations
      INTEGER  ::   js                   ! dummy loop indice
      INTEGER  ::   iend, iday0, iday1   ! temporary integers
      REAL(wp) :: zadatrjn, zadatrjb     ! adatrj at timestep kt-1 and kt-2 
      CHARACTER (len=25) :: charout
      !!----------------------------------------------------------------------

      ! 0.  initialization of adatrj0 and nday, nmonth,nyear, nday_year.
      !     ndastp has been initialized in domain.F90 or restart.F90
      !-----------------------------------------------------------------

      IF( kt == nit000 ) THEN

         IF( .NOT.ln_rstart )   adatrj0 = 0.e0      ! adatrj0 initialized in rst_read when restart 

         adatrj  = adatrj0
         nyear   =   ndastp / 10000
         nmonth  = ( ndastp - (nyear * 10000) ) / 100
         nday    =   ndastp - (nyear * 10000) - ( nmonth * 100 ) 

         ! Calculates nday_year, day since january 1st (useful to read  daily forcing fields)
         nday_year =  nday
         !                               ! accumulates days of previous months of this year
         DO js = 1, nmonth-1
            IF( nleapy == 1 .AND. MOD( nyear, 4 ) == 0 ) THEN
               nday_year = nday_year + nbiss(js)
            ELSE
               nday_year = nday_year + nobis(js)
            ENDIF
         END DO

      ENDIF

      ! I.  calculates adatrj, zadatrjn, zadatrjb.
      ! ------------------------------------------------------------------

      adatrj    = adatrj0 + ( kt - nit000 + 1 ) * rdttra(1) / rday
      zadatrjn  = adatrj0 + ( kt - nit000     ) * rdttra(1) / rday
      zadatrjb  = adatrj0 + ( kt - nit000 - 1 ) * rdttra(1) / rday


      ! II.  increment the date.  The date corresponds to 'now' variables (kt-1),
      !      which is the time step of forcing fields. 
      !      Do not do this at nit000  unless nrstdt= 2
      !      In that case ndastp (read in restart) was for step nit000-2
      ! -------------------------------------------------------------------

      iday0 = INT( zadatrjb )
      iday1 = INT( zadatrjn )

      IF( iday1 - iday0 >= 1 .AND. ( kt /= nit000 .OR. nrstdt == 2 ) ) THEN

         ! increase calendar
         nyear  =   ndastp / 10000
         nmonth = ( ndastp - (nyear * 10000) ) / 100
         nday   =   ndastp - (nyear * 10000) - ( nmonth * 100 ) 
         nday = nday + 1
         IF( nleapy == 1 .AND. MOD( nyear, 4 ) == 0 ) THEN
            iend = nbiss(nmonth)
         ELSEIF( nleapy > 1 ) THEN 
            iend = nleapy
         ELSE 
            iend = nobis(nmonth)
         ENDIF
         IF( nday == iend + 1 ) THEN
            nday  = 1
            nmonth = nmonth + 1
            IF( nmonth == 13 ) THEN
               nmonth  = 1
               nyear = nyear + 1
            ENDIF
         ENDIF
         ndastp = nyear * 10000 + nmonth * 100 + nday

         ! Calculates nday_year, day since january 1st (useful to read  daily forcing fields)
         nday_year =  nday
         !                                ! accumulates days of previous months of this year
         DO js = 1, nmonth-1
            IF( nleapy == 1 .AND. MOD( nyear, 4 ) == 0 ) THEN
               nday_year = nday_year + nbiss(js)
            ELSE
               nday_year = nday_year + nobis(js)
            ENDIF
         END DO

         IF(lwp) WRITE(numout,*)' ==============>> time-step =', kt, ' New day, DATE= ',   &
            &                   nyear, '/', nmonth, '/', nday, 'nday_year:', nday_year
      ENDIF

   END SUBROUTINE day

   !!======================================================================
END MODULE daymod
