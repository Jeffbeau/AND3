#ifdef key_RIVER_INPUT
MODULE rivers
!!=========================================================================================
!!Interior rivers input routines derived from DB's versionn improved from FD's code
!!Idea: put all routines in here and call them from the various routines that need them.
!!For SOPA-MC these are: traadv_tvd advection scheme, and dynspg_ts

!!NB: a time-dependent ramp factor is used (standard for SOPA-MC code)
!!    This must be changed if used outside of SOPA-MC code
!-------------------------------------------------------------------------------------------
! calls originating from:
!
!  step.F90
!  traadv_tvd.F90
!  dynspg_ts.F90
!  dynspg_flt.F90    !Only if  key_dynspg_flt is defined
!  dynnxt.F90
!
!-------------------------------------------------------------------------------------------
!  History:
   !!DB:2009.07.30
   !!Interior river input routines
   !!Idea: put all routines in here and call them from the various 
   !!routines that need them.
   !!For SOPA-MC these are: traadv_tvd advection scheme, and dynspg_ts

   !!NB: SLE = St.Lawrence Estuary is the only river input included
   !!    as of 2009.07.30
   !!NB: a time-dependent ramp factor is used (standard for SOPA-MC code)
   !!    This must be changed if used outside of SOPA-MC code

   !!2010.10.27 -- DBG'ing this routine, using FD code
   !!2010.11.01 -- DBG finished in this and associated routines
   !!JC: October-November 2010: Finishing multiple river input based on DB and FD approach
   !!2010.11.09 -- JC: Time interpolation (days) for runoff input

!-------------------------------------------------------------------------------------------

   !! * Modules used
   USE oce             ! ocean dynamics and active tracers
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O manager
   USE lib_mpp
   USE lbclnk          ! ocean lateral boundary condition (or mpp link) 
   USE daymod
   USE phycst

   IMPLICIT NONE
   PRIVATE
   
   INTEGER,PUBLIC, PARAMETER :: nrivermax=80,first_year_runoff=1914,last_year_runoff=2020
                                                       !!JC: a river can be spread accros
                                                       !! more that one cell...
   INTEGER,PUBLIC :: iriv_a(nrivermax),jriv_a(nrivermax)      !!first cell where river falls
   INTEGER,PUBLIC :: iriv_b(nrivermax),jriv_b(nrivermax)      !!last cell where river falls
   REAL(wp) :: RivFlow(nrivermax), vel_riv(nrivermax), da_riv(nrivermax)
   REAL(wp) :: width,hight                             !cell geometry
   INTEGER  :: ncell(nrivermax)                !number of cells for each river (e.g. St. Lawrence is 3)
   REAL(wp) :: river_monthly_flow(nrivermax,first_year_runoff:last_year_runoff,12)
   INTEGER,PUBLIC :: nriver,iside(nrivermax),ishift(nrivermax),jshift(nrivermax)
   INTEGER,PUBLIC :: mbathy_riv(nrivermax) !Probably useless for now, as min # layer=2
   LOGICAL :: first
   LOGICAL, PUBLIC :: do_river(nrivermax)
   REAL, PARAMETER :: FW_Salt=1.0 ! Freshwater salinity
!
!-------------------------------------------------------------------------------------------
!JC: Feb, 14, 2001.
!
   INTEGER,PUBLIC, DIMENSION( jpidta ) ::  mi0b, mi1b   !JC: global ==> local domain i-indice
   INTEGER,PUBLIC, DIMENSION( jpjdta ) ::  mj0b, mj1b   !JC: global ==> local domain j-indice
                                                 !JC: But without the first/last rows and colums
!-------------------------------------------------------------------------------------------

   PUBLIC Riv_init, riv_tra, riv_dyn, riv_dyna, riv_ts !JC: Fred's convention: Subroutines called by other modules
   PUBLIC river_set_mask,river_unset_mask

   data first /.true./  !Initilisation flag

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
!
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!   River flow input from dry cell    |
!   into wet cell.                    |       Velocity vector components distribution
!  iside(=1-4): From where river is   |
!           entering wet cell         |                    (N)
!                                     |                     ^
!              (N)                    |                     ^ 
!        ------ 1 --------            |              ----v(i,j,k)-----
!        |               |            |              |               |
!        |   River       |            |              |               |
!        |    cell       |            |              |               |
!    (W) 4   (on wet     2 (E)        |       (W) u(i-1,j,k)       u(i,j,k) >> (E)
!        |    point)     |            |              |               |
!        |               |            |              |               |
!        ------ 3 --------            |              ----v(i,j-1,k)---
!              (S)                    |                     (S)
!                                     |
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------


CONTAINS

     SUBROUTINE Riv_init(kt)
!!
!!DB,JC
!!Read river flow values
     
     INTEGER, INTENT( in ) ::  kt

     integer read_year,read_month
     INTEGER  :: ji, jj, jk              ! dummy loop indices
     INTEGER  :: ir                      ! dummy loop index for rivers
      

!-------------------------------------------------------------------------------------------
!
      !JC:Feb 14, 2011: Patch for rivers on tiles borders
      !
      ! data domain indices ==> local domain indices
      ! (return (m.0,m.1)=(2,1) if data domain gridpoint is to the west/south of the 
      ! local domain, or (m.0,m.1)=(jp.,jp.-1) to the east/north of local domain. 
      !excluding first/last row/column
      !July 20, 2011: Now using nlei and nlej to limit the do loops
!
      DO ji = 1, jpidta
        mi0b(ji) = MAX( 1, MIN( ji - jpizoom + 1 - nimpp + 1, nlei+1 ) )
        mi1b(ji) = MAX( 0, MIN( ji - jpizoom + 1 - nimpp + 1, nlei   ) )
        !mi0b(ji) = MAX( nldi, MIN( ji - jpizoom + 1 - nimpp + 1, nlei+1 ) )
        !mi1b(ji) = MAX( nldi-1, MIN( ji - jpizoom + 1 - nimpp + 1, nlei   ) )
      END DO
      DO jj = 1, jpjdta
        mj0b(jj) = MAX( 1, MIN( jj - jpjzoom + 1 - njmpp + 1, nlej+1 ) )
        mj1b(jj) = MAX( 0, MIN( jj - jpjzoom + 1 - njmpp + 1, nlej   ) )
        !mj0b(jj) = MAX( nldj, MIN( jj - jpjzoom + 1 - njmpp + 1, nlej+1 ) )
        !mj1b(jj) = MAX( nldj-1, MIN( jj - jpjzoom + 1 - njmpp + 1, nlej   ) )
      END DO
!
!-------------------------------------------------------------------------------------------
        river_monthly_flow(:,:,:)=-1.  !ARRAY! Initialization for check later.
        !open(66,file='Rivers.dat',STATUS='OLD')
        open(66,file='Damed_Rivers.dat',STATUS='OLD')
!        open(66,file='Rivers_RoPM.dat',STATUS='OLD')
!        open(66,file='Rivers_Rom.dat',STATUS='OLD')
!        open(66,file='Rivers_SMOB.dat',STATUS='OLD')
!        open(66,file='Rivers_alldamed.dat',STATUS='OLD')
!        open(66,file='Rivers_SL.dat',STATUS='OLD')
!        open(66,file='Rivers_natural.dat',STATUS='OLD')
!        print*,'Using Rivers natural, i.e. hydrological model with corrected SL obs'
!        print*,'Using Damed Rivers'
!DL        print*,'Using Rivers RoPM, i.e. damed Romaine and Petit-Mecatina'
!DL        print*,'Using Rivers SL'
!
        read(66,*) nriver  !number of rivers read in
        if(nriver > nrivermax) stop 'nrivermax not big enough'

        if(lwp) print*,'Number of rivers =',nriver
        do ir=1,nriver
            read(66,*)iriv_a(ir),jriv_a(ir),iriv_b(ir),jriv_b(ir),iside(ir) !global index of rivers
           if(lwp) print*, "river at:",iriv_a(ir),jriv_a(ir)," to ",iriv_b(ir),jriv_b(ir)
!           print *,'mj0b(jriv_a(ir)),mj1b(jriv_b(ir))=',mj0b(jriv_a(ir)),mj1b(jriv_b(ir))
!           print *,'mi0b(iriv_a(ir)), mi1b(iriv_b(ir))=',mi0b(iriv_a(ir)), mi1b(iriv_b(ir))
           if(lwp) print*,'iside(ir)=',iside(ir)
        enddo
        
        do
             read(66,*,end=555) read_year,read_month,(river_monthly_flow(ir,read_year,read_month),ir=1,nriver)
             if(read_year.lt.first_year_runoff) stop 'Runoff problem (1)'
             if(read_year.gt.last_year_runoff) stop 'Runoff problem (2)'
        enddo
555     close(66)
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!
!        if(lwp) print*, 'Warning: reducing river runoff by 30% (except St. Lawrence)'
!        if(lwp) print*, 'Warning: reducing river runoff by 30% (except St. Lawrence)'
!        if(lwp) print*, 'Warning: reducing river runoff by 30% (except St. Lawrence)'
!        if(lwp) print*, 'Warning: reducing river runoff by 30% (except St. Lawrence)'
!        if(lwp) print*, 'Warning: reducing river runoff by 30% (except St. Lawrence)'
!        river_monthly_flow(2:nrivermax,:,:)= 0.7*river_monthly_flow(2:nrivermax,:,:)
!
!        if(lwp) print*, 'Warning: increasing river runoff by 30% (except St. Lawrence)'
!        if(lwp) print*, 'Warning: increasing river runoff by 30% (except St. Lawrence)'
!        if(lwp) print*, 'Warning: increasing river runoff by 30% (except St. Lawrence)'
!        if(lwp) print*, 'Warning: increasing river runoff by 30% (except St. Lawrence)'
!        river_monthly_flow(2:nrivermax,:,:)= 1.3*river_monthly_flow(2:nrivermax,:,:)
!
!DL, june 2013, decreasing river flow of damne river to bring back to observed mean
        if(lwp) print*, 'Warning: decreasing river runoff of damned rivers'
        river_monthly_flow(25,:,:)= river_monthly_flow(25,:,:)*0.79    !Manicouagan
        river_monthly_flow(26,:,:)= river_monthly_flow(26,:,:)*0.76    !Outardes
        river_monthly_flow(29,:,:)= river_monthly_flow(29,:,:)*0.77    !Betsiamites
        river_monthly_flow(38,:,:)= river_monthly_flow(38,:,:)*0.95   !Saguenay
!
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!JC: To BE erased:
!JC:          !To check if there is input errors by comparing with orginal file (check last river)
!JC:              do read_year=first_year_runoff,last_year_runoff
!JC:              do read_month=1,12  
!JC:               if(lwp) write(2900,*) read_year,read_month,river_monthly_flow(78,read_year,read_month)
!JC:               enddo 
!JC:               enddo 
!JC:        flush(2900)
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!
          ncell(:)=0
          do_river=.false. !initialisation,array
          DO ir=1,nriver
!            Determine (i,j) of upstream cell
!
             if(iside(ir).eq.1) then
                ishift(ir)= 0
                jshift(ir)= 1
             endif
             if(iside(ir).eq.2) then
                ishift(ir)= 1
                jshift(ir)= 0
             endif
             if(iside(ir).eq.3) then
                ishift(ir)= 0
                jshift(ir)= -1
             endif
             if(iside(ir).eq.4) then
                ishift(ir)= -1
                jshift(ir)= 0
             endif

!
!            Determine number of cells per river and check river geometry.
!            Calculate simple vertical extension
!
                 do jj = mj0b(jriv_a(ir)),mj1b(jriv_b(ir))      
                    do ji = mi0b(iriv_a(ir)), mi1b(iriv_b(ir))
                    do_river(ir)=.true.
                    ncell(ir)=ncell(ir)+1  !summing number of cell for each river
                    !if(ji+ishift(ir).lt.1) then !boundary check
                    if(ji+ishift(ir).lt.nldi) then !boundary check
                         print*,'River problem: upstream point index lower than 1'
                         print*,'for river #',ir
                         print*,'ji,ji+ishift(ir)=',ji,ji+ishift(ir)
                    !no     stop ' stopping model'
                         do_river(ir)=.false.
                         ncell(ir)=ncell(ir)-1
                    endif
                    !if(ji+ishift(ir).gt.jpi) then !boundary check
                    if(ji+ishift(ir).gt.nlei) then !boundary check
                         print*,'River problem: upstream point greater index than nlei'
                         print*,'for river #',ir
                         print*,'ji,ji+ishift(ir)=',ji,ji+ishift(ir)
                     !no    stop ' stopping model'
                         do_river(ir)=.false.
                         ncell(ir)=ncell(ir)-1
                    endif
                    !if(jj+jshift(ir).lt.1) then !boundary check
                    if(jj+jshift(ir).lt.nldj) then !boundary check
                         print*,'River problem: upstream point index lower than 1'
                         print*,'for river #',ir
                         print*,'jj,jj+jshift(ir)=',jj,jj+jshift(ir)
                      !no   stop ' stopping model'
                         do_river(ir)=.false.
                         ncell(ir)=ncell(ir)-1
                    endif
                    !if(jj+jshift(ir).gt.jpj) then !boundary check
                    if(jj+jshift(ir).gt.nlej) then !boundary check
                         print*,'River problem: upstream point index greater than nlej'
                         print*,'for river #',ir
                         print*,'jj,jj+jshift(ir)=',jj,jj+jshift(ir)
                      !no   stop ' stopping model'
                         do_river(ir)=.false.
                         ncell(ir)=ncell(ir)-1
                    endif
                    if(tmask(ji,jj,1).ne.1.) then
                         print*,'River point is not wet'
                         print*,'river no: ',ir
                         stop ' stopping model'
                    endif
                    if(tmask(ji+ishift(ir),jj+jshift(ir),1).ne.0.) then
                         print*,'Upstream point for river is not dry'
                         print*,'river no: ',ir
                         stop ' stopping model'
                    endif
                    mbathy_riv(ir)=min(mbathy(ji,jj),2)  !JC: that could be useful later if one wants to remove
                                                         !    min number of layers set to 2...
                    !print*,'ir,mbathy_riv(ir)=', ir, mbathy_riv(ir)
                    if(mbathy_riv(ir).le.0) then
                       print*,'Problem with number of layers at river point'
                         print*,'river no: ',ir
                         stop ' stopping model'
                    endif
                    print*,'River #', ir, "is in... ji,jj,ishift,jshift,mi0b,mj0b,mi1b,mj1b,jpi,jpj",&
                        ji,jj,ishift(ir),jshift(ir),mi0b(iriv_a(ir)),mj0b(jriv_a(ir)),mi1b(iriv_b(ir)),mj1b(jriv_b(ir)),jpi,jpj
                 enddo
                 enddo
                 call mpp_sum(ncell(ir))     !!JC: to be safe, broadcast ncell to all processors
                                              !JC: In case the domain split happens at river points
           endDO
          DO ir=1,nriver !Doing another loop to avoid messed up printing. Done only one time
           if(lwp) then
                if(ncell(ir).eq.0) then
                 print*,"Warning river #",ir,"isnt in the model"
                 print*,"this will cause a division by zero while calculating vel_riv in riv_dyn"
                endif
                 print*,'ir,ncell(ir)=',ir,ncell(ir)
            endif
           endDO
!-------------------------------------------------------------------------------------------
!
        first=.false.  !Initilisation done only on first use
!
   END SUBROUTINE Riv_init

   SUBROUTINE riv_dyn (kt) !formerly river_SLE_01( kt )
!!
!!First routine, to be called before advection scheme
!! (NB: SOPA on MC-domain only works with traadv_tvd scheme)
!!

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step

      !! * Local declarations
      INTEGER  ::   ji, jj, jk              ! dummy loop indices
      INTEGER  ::   ir               ! dummy loop index for rivers

      if(first) call riv_init(kt) !initialization if first call
!!
      call RivFlow_val(kt)  !JC: no need to pass RivFlow we stay in module

!!    Convert to velocity in m/s
      DO ir=1,nriver
      da_riv(ir) = 0.0
      if(do_river(ir)) then
      do jj = mj0b(jriv_a(ir)),mj1b(jriv_b(ir))      
         do ji = mi0b(iriv_a(ir)), mi1b(iriv_b(ir))
         width=0.; hight=0.
           if(iside(ir).eq.1) width = e1v(ji,jj) !sum accross grid cell
           if(iside(ir).eq.2) width = e2u(ji,jj)  !where river is set
           if(iside(ir).eq.3) width = e1v(ji,jj-1)
           if(iside(ir).eq.4) width = e2u(ji-1,jj)
           do jk = 1, mbathy_riv(ir)
               hight= hight+e3t(jk)  !!NB: inaccurate if dz=dz(i,j,k)
           enddo
!           hight= hight+sshn(ji,jj)  !JC: add water level 
           da_riv(ir)=da_riv(ir)+ width*hight
         enddo
      enddo
      endif
      ENDDO
!
!-------------------------------------------------------------------------------------------
      DO ir=1,nriver  !New loop nescessary to sump da_riv(:) over all the processoer
         call mpp_sum(da_riv(ir))     !!to be safe, broadcast this to all processors 
                                      !!in case the domain split happens at river points
          vel_riv(ir) = RivFlow(ir)/da_riv(ir) 
      ENDDO
!-------------------------------------------------------------------------------------------
!
!
      DO ir=1,nriver
      if(do_river(ir)) then
      do jj = mj0b(jriv_a(ir)),mj1b(jriv_b(ir))      
         do ji = mi0b(iriv_a(ir)), mi1b(iriv_b(ir))
            do jk = 1, mbathy_riv(ir)
               if(iside(ir).eq.1) then
                  vn(ji,jj,jk) = -vel_riv(ir)
               endif
               if(iside(ir).eq.2) then
                  un(ji,jj,jk) = -vel_riv(ir)
               endif
               if(iside(ir).eq.3) then
                  vn(ji,jj-1,jk) = vel_riv(ir)
               endif
               if(iside(ir).eq.4) then
                  un(ji-1,jj,jk) = vel_riv(ir)
               endif
            enddo
         enddo
      enddo
      endif
      ENDDO
      !call lbc_lnk( un, 'U', 1. )       ! lateral boundary conditions 
      !call lbc_lnk( vn, 'V', 1. )

   END SUBROUTINE riv_dyn  

   SUBROUTINE riv_dyna (kt) !formerly none
!!
!!DB:
!! (NB: SOPA on MC-domain only works with traadv_tvd scheme)
!!
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step

      !! * Local declarations
      INTEGER  ::   ji, jj, jk              ! dummy loop indices
      INTEGER  ::   ir                      ! dummy loop index for rivers

      DO ir=1,nriver
      if(do_river(ir)) then
      do jj = mj0b(jriv_a(ir)),mj1b(jriv_b(ir))      
         do ji = mi0b(iriv_a(ir)), mi1b(iriv_b(ir))
            do jk = 1, mbathy_riv(ir)
               if(iside(ir).eq.1) then
                  va(ji,jj,jk) = -vel_riv(ir)
               endif
               if(iside(ir).eq.2) then
                  ua(ji,jj,jk) = -vel_riv(ir)
               endif
               if(iside(ir).eq.3) then
                  va(ji,jj-1,jk) = vel_riv(ir)
               endif
               if(iside(ir).eq.4) then
                  ua(ji-1,jj,jk) = vel_riv(ir)
               endif
            enddo
         enddo
      enddo
      endif
      ENDDO
      !call lbc_lnk( ua, 'U', -1. )       ! lateral boundary conditions 
      !call lbc_lnk( va, 'V', -1. )

   END SUBROUTINE riv_dyna  

   SUBROUTINE riv_tra (kt) ! Formerly river_SLE_01b( kt )
!!
!!DB:
!! (NB: SOPA on MC-domain only works with traadv_tvd scheme)
!!

      INTEGER, INTENT( in ) ::   kt         ! ocean time-step
      INTEGER  ::   ji, jj, jk              ! dummy loop indices
      INTEGER  ::   ir                      ! dummy loop index for rivers

      if(first) call riv_init(kt) !initialization if first call
      DO ir=1,nriver
      if(do_river(ir)) then
      do jj = mj0b(jriv_a(ir)),mj1b(jriv_b(ir))      
         do ji = mi0b(iriv_a(ir)), mi1b(iriv_b(ir))
            do jk = 1, mbathy_riv(ir)
                tb(ji+ishift(ir),jj+jshift(ir),jk) = max(0.,tb(ji,jj,jk))  !no T change due to river
                tn(ji+ishift(ir),jj+jshift(ir),jk) = max(0.,tn(ji,jj,jk))  !no T change due to river
                sb(ji+ishift(ir),jj+jshift(ir),jk) = FW_Salt
                sn(ji+ishift(ir),jj+jshift(ir),jk) = FW_Salt
            enddo
!-------------------------------------------------------------------------------------------
            ! if(ir.eq.12) then
            !    print*,"ji,jj=",ji,jj
            !    print*,"tn,sn",tn(ji,jj,1),sn(ji,jj,1)
            !    print*,"tn+shift,sn+shift",tn(ji+ishift(ir),jj+jshift(ir),1),sn(ji+ishift(ir),jj+jshift(ir),1)
            ! endif
!-------------------------------------------------------------------------------------------
         enddo
      enddo
      endif
      ENDDO
!      call lbc_lnk( tn, 'T', 1. )       ! lateral boundary conditions 
!      call lbc_lnk( tb, 'T', 1. )
!      call lbc_lnk( sn, 'T', 1. )       ! lateral boundary conditions 
!      call lbc_lnk( sb, 'T', 1. )


   END SUBROUTINE riv_tra


   SUBROUTINE riv_ts (ua_e,va_e)   !JC: formerly river_SLE_03(ua_e,va_e) !JC: General case for many rivers
!!
!!routine called by dynspg_ts
!!
     !! * Arguments
     REAL(wp), DIMENSION(jpi,jpj) ::  ua_e,va_e

     !! * Local declarations
     INTEGER  ::   ji, jj, jk              ! dummy loop indices
     INTEGER  ::   ir                      ! dummy loop index for rivers
      DO ir=1,nriver
      if(do_river(ir)) then
      do jj = mj0b(jriv_a(ir)),mj1b(jriv_b(ir))      
         do ji = mi0b(iriv_a(ir)), mi1b(iriv_b(ir))
            if(iside(ir).eq.1) va_e(ji,jj) = -RivFlow(ir)/ncell(ir)/e1v(ji,jj)
            if(iside(ir).eq.2) ua_e(ji,jj) = -RivFlow(ir)/ncell(ir)/e2u(ji,jj)
            if(iside(ir).eq.3) va_e(ji,jj-1) = RivFlow(ir)/ncell(ir)/e1v(ji,jj-1)
            if(iside(ir).eq.4) ua_e(ji-1,jj) = RivFlow(ir)/ncell(ir)/e2u(ji-1,jj)
!-------------------------------------------------------------------------------------------
      !       if(ir.eq.12) then
      !          print*,"ji,jj=",ji,jj
      !          print*,"ua_e,va_e",ua_e(ji,jj),va_e(ji,jj)
      !          print*,"tn,sn",tn(ji,jj,1),sn(ji,jj,1)
      !          print*,"tn+shift,sn+shift",tn(ji+ishift(ir),jj+jshift(ir),1),sn(ji+ishift(ir),jj+jshift(ir),1)
      !       endif
!-------------------------------------------------------------------------------------------
         enddo
      enddo
      endif
      ENDDO
      !call lbc_lnk( ua_e  , 'U', -1. )
      !call lbc_lnk( va_e  , 'V', -1. )
      !JC:no need to call lbc_lnk since it is done right after the call in dynspg_ts

   END SUBROUTINE riv_ts


   SUBROUTINE RivFlow_val(kt)
!!
!!DB
!!TO DO: interpolate monthly values to nday_year
!!JC: Done
     
     INTEGER, INTENT( in ) ::  kt
     INTEGER  :: i,ji, jj, jk              ! dummy loop indices
     INTEGER  :: ir                        ! dummy loop index for rivers
   !! * Local declarations
      integer ndm(12),month_before,month_after,year_before,year_after
      data ndm/31,28,31,30,31,30,31,31,30,31,30,31/ !JC: Don't care about leap year for now
      integer :: y1,y2,m1,m2     !temporary variables for weights
      real :: weight

!-------------------------------------------------------------------------------------------
!JC:     Weights for time interpolation. Stop at the day level. Runoff assumed to be monthly mean
!        That means that runoff will be constant on specific day.
!
      year_before=nyear
      year_after=nyear
      month_before=nmonth-1
      month_after=nmonth+1

      if(nday.ge.real(ndm(nmonth))/2.) then
         if(nmonth.eq.12)then
             year_after=year_after+1
             month_after=1
         endif
!
         if(year_after.gt.last_year_runoff) stop 'Runoff problem(3)'
!
         weight = (real(nday)-real(ndm(nmonth))/2.)/ (real(ndm(nmonth)+ndm(month_after))/2.)
!
         y2=year_after
         m2=month_after
         y1=nyear
         m1=nmonth
      else
         if(nmonth.eq.1)then
             year_before=year_before-1
             month_before=12
         endif
         y2=nyear
         m2=nmonth
         y1=year_before
         m1=month_before
!
         if(year_before.lt.first_year_runoff) stop 'Runoff problem(5)'
!
          weight = (real(ndm(month_before))/2.+real(nday))/(real(ndm(nmonth)+ndm(month_before))/2.)
      end if
!-------------------------------------------------------------------------------------------
!JC:
!Return river flow.

     DO ir=1,nriver
        if(river_monthly_flow(ir,y2,m2).lt.0.) stop 'Model wants out of range runoff: stop'
!       !JC: Runoff now in  m/s, not mSv
        RivFlow(ir) = ramp*(river_monthly_flow(ir,y2,m2) *weight + river_monthly_flow(ir,y1,m1)*(1.-weight))
        if(lwp.and.ir.eq.1) write(2901,*) kt,RivFlow(ir)
        if(lwp.and.ir.eq.78) write(2902,*) kt,RivFlow(ir)
     ENDDO
     
   END SUBROUTINE RivFlow_val

   Subroutine  river_set_mask
     !! * Local declarations
     INTEGER  ::   ji, jj, jk              ! dummy loop indices
     INTEGER  ::   ir               ! dummy loop index for rivers

      DO ir=1,nriver
      if(do_river(ir)) then
      do jj = mj0b(jriv_a(ir)),mj1b(jriv_b(ir))      
         do ji = mi0b(iriv_a(ir)), mi1b(iriv_b(ir))
            do jk = 1, mbathy_riv(ir)
                  tmask(ji+ishift(ir),jj+jshift(ir),jk) = 1.0
            enddo
         enddo
      enddo
      endif
      ENDDO
      !call lbc_lnk( tmask, 'T', 1. )       ! lateral boundary conditions 

   End Subroutine  river_set_mask

   Subroutine  river_unset_mask
          !! * Local declarations

     INTEGER  ::   ji, jj, jk              ! dummy loop indices
     INTEGER  ::   ir               ! dummy loop index for rivers
      DO ir=1,nriver
      if(do_river(ir)) then
      do jj = mj0b(jriv_a(ir)),mj1b(jriv_b(ir))      
         do ji = mi0b(iriv_a(ir)), mi1b(iriv_b(ir))
            do jk = 1, mbathy_riv(ir)
                  tmask(ji+ishift(ir),jj+jshift(ir),jk) = 0.0
            enddo
         enddo
      enddo
      endif
      ENDDO
      !call lbc_lnk( tmask, 'T', 1. )       ! lateral boundary conditions 

    End  subroutine river_unset_mask

!!=========================================================================================
END MODULE rivers

#endif
