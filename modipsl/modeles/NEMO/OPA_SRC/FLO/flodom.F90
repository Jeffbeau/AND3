!!DB 2008.05.16
!!Modified routine for more specific DB particle tracking needs
!!For older code, see OLD_CODE/ or look elsewhere for an older version


MODULE flodom
   !!======================================================================
   !!                       ***  MODULE  flodom  ***
   !! Ocean floats :   domain
   !!======================================================================
#if   defined key_floats   ||   defined key_esopa
   !!----------------------------------------------------------------------
   !!   'key_floats'                                     float trajectories
   !!----------------------------------------------------------------------
   !!   flo_dom        : initialization of floats
   !!   findmesh       : compute index of position 
   !!   dstnce         : compute distance between face mesh and floats 
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE flo_oce         ! ocean drifting floats
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distribued memory computing library

   IMPLICIT NONE

   !! * Accessibility
   PRIVATE  dstnce
   PUBLIC flo_dom     ! routine called by floats.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/FLO/flodom.F90,v 1.5 2005/09/22 10:24:53 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

  SUBROUTINE flo_dom
      !! ---------------------------------------------------------------------
      !!                  ***  ROUTINE flo_dom  ***
      !!                 
      !!  ** Purpose :   Initialisation of floats
      !!
      !!  ** Method  :   We put the floats  in the domain with the latitude,
      !!       the longitude (degree) and the depth (m).
      !!
      !!----------------------------------------------------------------------      
      !! * Local declarations
     LOGICAL  :: llinmesh
     CHARACTER (len=21) ::  clname
     INTEGER  :: ji, jj, jk               ! DO loop index on 3 directions
     INTEGER  :: jfl, jfl1                ! number of floats   
     INTEGER  :: inum = 11                ! logical unit for file read
!!DB
     INTEGER, PARAMETER :: max_num_floats = 100000
!     INTEGER, DIMENSION (jpnfl    )  ::   &
     INTEGER, DIMENSION (:),ALLOCATABLE  ::   &
          iimfl, ijmfl, ikmfl,    &          ! index mesh of floats
          idomfl,  ivtest, ihtest
     REAL(wp) :: zdxab, zdyad
!!DB
     INTEGER :: num_floats
     INTEGER, DIMENSION(:),ALLOCATABLE  ::    &
          tmp_nisobfl,    &  ! 0 for a isobar float
          !              ! 1 for a float following the w velocity
          tmp_ngrpfl         ! number to identify searcher group
     REAL(wp), DIMENSION(:),ALLOCATABLE ::    &
          tmp_flxx,       &  ! longitude of float (decimal degree)
          tmp_flyy,       &  ! latitude of float (decimal degree)
          tmp_flzz 

     !!---------------------------------------------------------------------
     
     ! Initialisation with the geographical position or restart
     
     IF(lwp) WRITE(numout,*) 'flo_dom : compute initial position of floats'
     IF(lwp) WRITE(numout,*) '~~~~~~~~'
     IF(lwp) WRITE(numout,*) '           jpnfl determined from init_floats'
     
     
     IF(lwp) WRITE(numout,*) '                     init_float read '

! First initialisation of floats
! the initials positions of floats are written in a file
! with a variable to know if it is a isobar float a number 
! to identified who want the trajectories of this float and 
! an index for the number of the float         
! open the init file 
!!DB: I keep the above, although it is not clear why
!!Change so that float file is: flxx  flyy  flzz  nisobfl  ngrpfl -- 1 for each jpnfl floats
     ALLOCATE(tmp_nisobfl(max_num_floats)); ALLOCATE(tmp_ngrpfl(max_num_floats))
     ALLOCATE(tmp_flxx(max_num_floats));ALLOCATE(tmp_flyy(max_num_floats))
     ALLOCATE(tmp_flzz(max_num_floats))
     
     clname='init_float'
     num_floats = 0
     open(inum,FILE=clname,STATUS='old',ERR=222)   !!Allow for no-file-found
     do while(num_floats >= 0 .AND. num_floats <= max_num_floats)
        num_floats = num_floats + 1
        jfl = num_floats      !lazy
        read(inum,*,END=223) tmp_flxx(jfl),tmp_flyy(jfl),tmp_flzz(jfl),tmp_nisobfl(jfl),tmp_ngrpfl(jfl)
     enddo
223  close(inum)
     num_floats = num_floats - 1
     if(lwp)write(numout,*)'Found ',num_floats,' floats'
     goto 224
!!DB file-not-found
222  if(lwp)write(numout,*)'NB: Float file init_floats not found'

!!DB Allocate arrays here
224  if(num_floats == 0) then
!!DB adjust allocation -- or not as it seems like allocating zero-size arrays is OK
        jpnfl = 0 ! Force correct value for jpnfl as it is needed elsewhere
     else
!!DB allocate arrays to num_floats elements and set jpnfl = num_floats
        jpnfl = num_floats
        ALLOCATE(iimfl(jpnfl));   ALLOCATE(ijmfl(jpnfl));   ALLOCATE(ikmfl(jpnfl))
        ALLOCATE(idomfl(jpnfl));  ALLOCATE(ivtest(jpnfl));  ALLOCATE(ihtest(jpnfl))
        ALLOCATE(nisobfl(jpnfl)); ALLOCATE(ngrpfl(jpnfl))
        ALLOCATE(flxx(jpnfl));    ALLOCATE(flyy(jpnfl));    ALLOCATE(flzz(jpnfl))
        ALLOCATE(tpifl(jpnfl));   ALLOCATE(tpjfl(jpnfl));   ALLOCATE(tpkfl(jpnfl))
!!Assign working arrays
        nisobfl(1:jpnfl) = tmp_nisobfl(1:jpnfl); ngrpfl(1:jpnfl) = tmp_ngrpfl(1:jpnfl); 
        flxx(1:jpnfl) = tmp_flxx(1:jpnfl);flyy(1:jpnfl) = tmp_flyy(1:jpnfl);flzz(1:jpnfl) = tmp_flzz(1:jpnfl);
     endif
     DEALLOCATE(tmp_nisobfl); DEALLOCATE(tmp_ngrpfl);
     DEALLOCATE(tmp_flxx);DEALLOCATE(tmp_flyy);DEALLOCATE(tmp_flzz)

     
     ! Test to find the grid point coordonate with the geographical position         
     DO jfl = 1, jpnfl

        ihtest(jfl) = 0
        ivtest(jfl) = 0
        ikmfl(jfl) = 0
# if   defined key_mpp_mpi   ||   defined key_mpp_shmem
        DO ji = MAX(nldi,2), nlei
           DO jj = MAX(nldj,2), nlej
# else
        DO ji = 2, jpi
           DO jj = 2, jpj
# endif                  
                    ! for each float we find the indexes of the mesh 

                    !!DB: Recalling:
                    !                               Vij
                    !                       ---------------Fij
                    !                       |              |
                    !                       |              |
                    !                U i-1,j|       Tij    |Uij
                    !                       |              |
                    !                       |              |
                    !                       ----------------
                    !                               Vi,j-1
                    !
                    ! ====> if .true. then float is in T-cell(ji,jj); ji,jj local indices


              CALL findmesh(glamf(ji-1,jj-1),gphif(ji-1,jj-1),   &
                   glamf(ji-1,jj  ),gphif(ji-1,jj  ),   &
                   glamf(ji  ,jj  ),gphif(ji  ,jj  ),   &
                   glamf(ji  ,jj-1),gphif(ji  ,jj-1),   &
                   flxx(jfl)       ,flyy(jfl)       ,   &
                   glamt(ji  ,jj  ),gphit(ji  ,jj  ), llinmesh)
              IF(llinmesh) THEN
                 iimfl(jfl)  = ji
                 ijmfl(jfl)  = jj
                 ihtest(jfl) = ihtest(jfl)+1
                 DO jk = 1, jpk-1
                    IF( (fsdepw(ji,jj,jk) <= flzz(jfl)) .AND. (fsdepw(ji,jj,jk+1) >  flzz(jfl)) ) THEN
                       ikmfl(jfl)  = jk
                       ivtest(jfl) = ivtest(jfl) + 1
                    ENDIF
                 END DO
              ENDIF
           END DO
        END DO
        
        ! If the float is in a mesh computed by an other processor we put iimfl=ijmfl=-1            
        IF( ihtest(jfl) == 0 ) THEN
           iimfl(jfl) = -1
           ijmfl(jfl) = -1
        ENDIF
     END DO
     
     ! A zero in the sum of the arrays "ihtest" and "ivtest"          
     IF( lk_mpp )   CALL mpp_sum(ihtest,jpnfl)   ! sums over the global domain
     IF( lk_mpp )   CALL mpp_sum(ivtest,jpnfl)
     
     DO jfl = 1, jpnfl
        IF( (ihtest(jfl) > 1 ) .OR. ( ivtest(jfl) > 1 )) THEN
           IF(lwp) WRITE(numout,*) 'THE FLOAT',jfl,' IS NOT IN ONLY ONE MESH'
        ENDIF
        IF( ihtest(jfl) == 0 ) THEN 
           IF(lwp) WRITE(numout,*)'THE FLOAT',jfl,' IS IN NO MESH'
        ENDIF
     END DO
     
     ! We compute the distance between the float and the face of  the mesh         
     DO jfl = 1, jpnfl
        ! Made only if the float is in the domain of the processor
        IF( (iimfl(jfl) >= 0 ) .AND. ( ijmfl(jfl) >= 0 ) ) THEN
           
           ! TEST TO KNOW IF THE FLOAT IS NOT INITIALISED IN THE COAST
           
           idomfl(jfl) = 0
           IF( tmask(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)) == 0. ) idomfl(jfl)=1
           
           ! Computation of the distance between the float
           ! and the faces of the mesh
           !            zdxab
           !             .
           !        B----.---------C
           !        |    .         |
           !        |<------>flo   |
           !        |        ^     |
           !        |        |.....|....zdyad
           !        |        |     |
           !        A--------|-----D
           
           zdxab = dstnce(flxx(jfl),flyy(jfl),glamf(iimfl(jfl)-1,ijmfl(jfl)-1),flyy(jfl))                
           zdyad = dstnce(flxx(jfl),flyy(jfl),flxx(jfl),gphif(iimfl(jfl)-1,ijmfl(jfl)-1))
           
           ! Translation of this distances (in meter) in indexes
           
           !!DB: If I take tp?fl to denote T-cell, then T-cell(i,j) goes from
           !!     (i-1/2 ... i+1/2 , j-1/2 ... j+1/2) 
           
           tpifl(jfl) = (iimfl(jfl)-0.5)+zdxab/ e1u(iimfl(jfl)-1,ijmfl(jfl))+(mig(1)-jpizoom)
           tpjfl(jfl) = (ijmfl(jfl)-0.5)+zdyad/ e2v(iimfl(jfl),ijmfl(jfl)-1)+(mjg(1)-jpjzoom)
           tpkfl(jfl) = (fsdepw(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)+1) - flzz(jfl))*(ikmfl(jfl))                     &
                / (fsdepw(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)+1) - fsdepw(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)))   &
                + (flzz(jfl) - fsdepw(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)))*(ikmfl(jfl)+1)                     &
                / (fsdepw(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)+1) - fsdepw(iimfl(jfl),ijmfl(jfl),ikmfl(jfl)))
        ELSE
           tpifl (jfl) = 0.e0
           tpjfl (jfl) = 0.e0
           tpkfl (jfl) = 0.e0
           idomfl(jfl) = 0
        ENDIF
     END DO
     
     ! The sum of all the arrays tpifl, tpjfl, tpkfl give 3 arrays with the positions of all the floats. 
     IF( lk_mpp )   CALL mpp_sum( tpifl , jpnfl )   ! sums over the global domain
     IF( lk_mpp )   CALL mpp_sum( tpjfl , jpnfl )
     IF( lk_mpp )   CALL mpp_sum( tpkfl , jpnfl )
     IF( lk_mpp )   CALL mpp_sum( idomfl, jpnfl )
     
     !!DB
     ! Print the initial positions of the floats
     if( .NOT. ln_rstflo ) THEN 
        !!DB         ! WARNING : initial position not in the sea .OR. it's ijk position
        DO jfl = 1, jpnfl
           IF( idomfl(jfl) == 1 ) THEN
              IF(lwp) WRITE(numout,*)'*****************************'
              IF(lwp) WRITE(numout,*)'!!!!!!!  WARNING   !!!!!!!!!!'
              IF(lwp) WRITE(numout,*)'*****************************'
              IF(lwp) WRITE(numout,*)'The float number',jfl,'is out of the sea.'
              IF(lwp) WRITE(numout,*)'geographical position',flxx(jfl),flyy(jfl),flzz(jfl)
              IF(lwp) WRITE(numout,*)'index position',tpifl(jfl),tpjfl(jfl),tpkfl(jfl)
           else
              IF(lwp) WRITE(numout,*)'DBG: Float number ',jfl
              IF(lwp) WRITE(numout,*)'Geographical position',flxx(jfl),flyy(jfl),flzz(jfl)
              IF(lwp) WRITE(numout,*)'Index position',tpifl(jfl),tpjfl(jfl),tpkfl(jfl)
           endif
           
        END DO
     endif
     
   END SUBROUTINE flo_dom
         

   SUBROUTINE findmesh( pax, pay, pbx, pby,   &
                        pcx, pcy, pdx, pdy,   &
                        px  ,py  ,ptx, pty, ldinmesh )
      !! -------------------------------------------------------------
      !!                ***  ROUTINE findmesh  ***
      !!      
      !! ** Purpose :   Find the index of mesh for the point spx spy
      !!
      !! ** Method  : 
      !!
      !! History :
      !!   8.0  !  98-07 (Y.Drillet)  Original code
      !!----------------------------------------------------------------------
      !! * Arguments
      REAL(wp) ::   &
         pax, pay, pbx, pby,    &     ! ???
         pcx, pcy, pdx, pdy,    &     ! ???
         px, py,                &     ! longitude and latitude
         ptx, pty                     ! ???
      LOGICAL ::  ldinmesh            ! ???

      !! * local declarations
      REAL(wp) ::   &
         zabt, zbct, zcdt, zdat, zabpt, zbcpt, zcdpt, zdapt,  &
         psax,psay,psbx,psby,psx,psy
      REAL(wp) ::  fsline                ! Statement function

      !! * Substitutions
      fsline(psax, psay, psbx, psby, psx, psy) = psy  * ( psbx - psax )   &
                                               - psx  * ( psby - psay )   &
                                               + psax *   psby - psay * psbx
      !!---------------------------------------------------------------------
      
      ! 4 semi plane defined by the 4 points and including the T point
      zabt = fsline(pax,pay,pbx,pby,ptx,pty)
      zbct = fsline(pbx,pby,pcx,pcy,ptx,pty)
      zcdt = fsline(pcx,pcy,pdx,pdy,ptx,pty)
      zdat = fsline(pdx,pdy,pax,pay,ptx,pty)
      
      ! 4 semi plane defined by the 4 points and including the extrememity
      zabpt = fsline(pax,pay,pbx,pby,px,py)
      zbcpt = fsline(pbx,pby,pcx,pcy,px,py)
      zcdpt = fsline(pcx,pcy,pdx,pdy,px,py)
      zdapt = fsline(pdx,pdy,pax,pay,px,py)
       
      ! We compare the semi plane T with the semi plane including the point
      ! to know if it is in this  mesh.
      ! For numerical reasons it is possible that for a point which is on
      ! the line we don't have exactly zero with fsline function. We want 
      ! that a point can't be in 2 mesh in the same time, so we put the 
      ! coefficient to zero if it is smaller than 1.E-12
      
      IF( ABS(zabpt) <= 1.E-12 ) zabpt = 0.
      IF( ABS(zbcpt) <= 1.E-12 ) zbcpt = 0.
      IF( ABS(zcdpt) <= 1.E-12 ) zcdpt = 0.
      IF( ABS(zdapt) <= 1.E-12 ) zdapt = 0.
      IF( (zabt*zabpt >  0.) .AND. (zbct*zbcpt >= 0. ) .AND. ( zcdt*zcdpt >= 0. ) .AND. ( zdat*zdapt > 0. )   &
         .AND. ( px <= MAX(pcx,pdx) ) .AND. ( px >= MIN(pax,pbx) )    &
         .AND. ( py <= MAX(pby,pcy) ) .AND. ( py >= MIN(pay,pdy) ) ) THEN
         ldinmesh=.TRUE.
      ELSE
         ldinmesh=.FALSE.
      ENDIF

   END SUBROUTINE findmesh


   FUNCTION dstnce( pla1, phi1, pla2, phi2 )
      !! -------------------------------------------------------------
      !!                 ***  Function dstnce  ***
      !!          
      !! ** Purpose :   returns distance (in m) between two geographical
      !!                points
      !! ** Method  : 
      !!         
      !!----------------------------------------------------------------------
      !! * Arguments
      REAL(wp), INTENT(in) ::   pla1, phi1, pla2, phi2   ! ???

      !! * Local variables
      REAL(wp) ::   dly1, dly2, dlx1, dlx2, dlx, dls, dld, dpi
      REAL(wp) ::   dstnce
      !!---------------------------------------------------------------------
      
      dpi  = 2.* ASIN(1.)
      dls  = dpi / 180.
      dly1 = phi1 * dls
      dly2 = phi2 * dls
      dlx1 = pla1 * dls
      dlx2 = pla2 * dls

      dlx = SIN(dly1) * SIN(dly2) + COS(dly1) * COS(dly2) * COS(dlx2-dlx1)
 
      IF( ABS(dlx) > 1.0 ) dlx = 1.0

      dld = ATAN(DSQRT( ( 1-dlx )/( 1+dlx ) )) * 222.24 / dls
      dstnce = dld * 1000.

   END FUNCTION dstnce

#  else
   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE flo_dom                 ! Empty routine
   END SUBROUTINE flo_dom
#endif

   !!======================================================================
END MODULE flodom
