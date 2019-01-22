MODULE solver
   !!======================================================================
   !!                     ***  MODULE  solver  ***
   !! Ocean solver :  initialization of ocean solver
   !!=====================================================================

   !!----------------------------------------------------------------------
   !!   solver_init: solver initialization
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE zdf_oce         ! ocean vertical physics variables
   USE sol_oce         ! solver variables
   USE solmat          ! ???
   USE solisl          ! ???
   USE obc_oce         ! Lateral open boundary condition
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp
   USE dynspg_oce      ! choice/control of key cpp for surface pressure gradient

   IMPLICIT NONE

   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/SOL/solver.F90,v 1.13 2006/03/20 17:27:15 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE solver_init( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE solver_init  ***
      !!                   
      !! ** Purpose :   Initialization for the solver of the elliptic equation:
      !!       * default option: barotropic stream function system
      !!         and islands initialization (if lk_isl=T)
      !!       * lk_dynspg_flt = T : transport divergence system. No specific
      !!         treatment of islands.
      !!      
      !! ** Method :
      !!       - Compute the local depth of the water column at u- and v-point
      !!      (lk_dynspg_flt = T) or its inverse (lk_dynspg_rl = T).
      !!      The local depth of the water column is computed by summing 
      !!      the vertical scale factors. For its inverse, the thickness of
      !!      the first model level is imposed as lower bound. The inverse of
      !!      this depth is THEN taken and masked, so that the inverse of the
      !!      local depth is zero when the local depth is zero.
      !!       - Construct the matrix of the elliptic system by a call to
      !!      solmat.F routine.
      !!       - island (if lk_isl=T)
      !!            isl_dom: find islands from the bathymetry file
      !!            isl_bsf: compute the island barotropic stream function
      !!            isl_mat: compute the inverse island matrix
      !!            set mbathy to the number of non-zero w-levels of a water
      !!            column (the minimum value of mbathy is 2):
      !!                  mbathy = min( mbathy, 1 ) + 1
      !!
      !! ** Action : - hur, hvr : masked inverse of the local depth at
      !!                                u- and v-point. (lk_dynspg_rl = T)
      !!             - hu, hv   : masked local depth at u- and v- points
      !!                                (lk_dynspg_flt = T)
      !!             - c_solver_pt : nature of the gridpoint at which the
      !!                                solver is applied
      !! References :
      !!      Jensen, 1986: adv. phys. oceanogr. num. mod.,ed. o brien,87-110.
      !!      Madec & Marti, 1990: internal rep. LODYC, 90/03., 29pp.
      !!
      !! History :
      !!        !  90-10  (G. Madec)  Original code           
      !!        !  93-02  (O. Marti)                         
      !!        !  97-02  (G. Madec)  local depth inverse computation
      !!        !  98-10  (G. Roullet, G. Madec)  free surface 
      !!   9.0  !  03-07  (G. Madec)  free form, F90
      !!    "   !  05-11  (V. Garnier) Surface pressure gradient organization
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT(in) :: kt

      !! * Local declarations
      INTEGER :: ji, jj   ! dummy loop indices
      CHARACTER(len=80) :: clname

      NAMELIST/namsol/ nsolv, nsol_arp, nmin, nmax, nmod, eps, resmax, sor, epsisl, nmisl, rnu
      !!----------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'solver_init : solver to compute the surface pressure gradient'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'

      ! open elliptic solver statistics file
      clname = 'solver.stat'
      CALL ctlopn( numsol, clname, 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL',   &
                   1, numout, lwp, 1 )


      ! 0. Define the solver parameters
      !    ----------------------------
      ! Namelist namsol : elliptic solver / islands / free surface
      REWIND( numnam )
      READ  ( numnam, namsol )

#if defined key_feti
      ! FETI algorithm, we force nsolv at 3
      nsolv = 3
#endif


      ! 0. Parameter control and print
      !    ---------------------------

      ! Control print
      IF(lwp) WRITE(numout,*) '          Namelist namsol : set solver parameters'

      IF(lwp) THEN
         WRITE(numout,*) '             type of elliptic solver            nsolv    = ', nsolv
         WRITE(numout,*) '             absolute/relative (0/1) precision  nsol_arp = ', nsol_arp
         WRITE(numout,*) '             minimum iterations for solver      nmin     = ', nmin
         WRITE(numout,*) '             maximum iterations for solver      nmax     = ', nmax
         WRITE(numout,*) '             frequency for test                 nmod     = ', nmod
         WRITE(numout,*) '             absolute precision of solver       eps      = ', eps
         WRITE(numout,*) '             absolute precision for SOR solver  resmax   = ', resmax
         WRITE(numout,*) '             optimal coefficient of sor         sor      = ', sor
         IF(lk_isl) WRITE(numout,*) '             absolute precision stream fct    epsisl   = ', epsisl
         IF(lk_isl) WRITE(numout,*) '             maximum pcg iterations island    nmisl    = ', nmisl
         WRITE(numout,*) '             free surface parameter         rnu    = ', rnu
         WRITE(numout,*)
      ENDIF

      IF( lk_dynspg_flt ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '          free surface formulation'
         IF( lk_isl ) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) ' key_islands inconsistent with key_dynspg_flt'
            nstop = nstop + 1
         ENDIF
      ELSEIF( lk_dynspg_rl ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '          Rigid lid formulation'
      ELSE
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          Choose only one surface pressure gradient calculation: filtered or rigid-lid'
         IF(lwp) WRITE(numout,*) '          Should not call this routine if dynspg_exp or dynspg_ts has been chosen'
         nstop = nstop + 1
      ENDIF
      IF( lk_dynspg_flt .AND. lk_dynspg_rl ) THEN
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          Chose between free surface or rigid-lid, not both'
         nstop = nstop + 1
      ENDIF

      SELECT CASE ( nsolv )

      CASE ( 1 )                ! preconditioned conjugate gradient solver
         IF(lwp) WRITE(numout,*) '          a preconditioned conjugate gradient solver is used'
         IF( jpr2di /= 0 .AND. jpr2dj /= 0 ) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) ' jpr2di and jpr2dj should be equal to zero'
            nstop = nstop + 1
         ENDIF

      CASE ( 2 )                ! successive-over-relaxation solver
         IF(lwp) WRITE(numout,*) '          a successive-over-relaxation solver is used'
         IF( jpr2di /= 0 .AND. jpr2dj /= 0 ) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) ' jpr2di and jpr2dj should be equal to zero'
            nstop = nstop + 1
         ENDIF

      CASE ( 3 )                ! FETI solver
         IF(lwp) WRITE(numout,*) '          the FETI solver is used'
         IF( jpr2di /= 0 .AND. jpr2dj /= 0 ) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) ' jpr2di and jpr2dj should be equal to zero'
            nstop = nstop + 1
         ENDIF
         IF( .NOT.lk_mpp ) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) ' The FETI algorithm is used only with the key_mpp_... option'
            nstop = nstop + 1
         ELSE
            IF( jpnij == 1 ) THEN
               IF(lwp) WRITE(numout,cform_err)
               IF(lwp) WRITE(numout,*) ' The FETI algorithm needs more than one processor'
               nstop = nstop + 1
            ENDIF
         ENDIF
         
      CASE ( 4 )                ! successive-over-relaxation solver with extra outer halo
         IF(lwp) WRITE(numout,*) '          a successive-over-relaxation solver with extra outer halo is used'
         IF(lwp) WRITE(numout,*) '          with jpr2di =', jpr2di, ' and  jpr2dj =', jpr2dj
         IF( .NOT. lk_mpp .AND. jpr2di /= 0 .AND. jpr2dj /= 0 ) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) ' jpr2di and jpr2dj are not equal to zero'
            IF(lwp) WRITE(numout,*) ' In this case this algorithm should be used only with the key_mpp_... option'
            nstop = nstop + 1
         ELSE
            IF( ( ( jperio == 1 .OR. jperio == 4 .OR. jperio == 6 ) .OR. ( jpni /= 1 ) ) &
              &  .AND. ( jpr2di /= jpr2dj ) ) THEN  
               IF(lwp) WRITE(numout,cform_err)
               IF(lwp) WRITE(numout,*) '          jpr2di should be equal to jpr2dj'
               nstop = nstop + 1
            ENDIF
         ENDIF

      CASE DEFAULT
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) '          bad flag value for nsolv = ', nsolv
         nstop = nstop + 1
         
      END SELECT

      ! Grid-point at which the solver is applied
      ! -----------------------------------------

      IF( lk_dynspg_rl ) THEN       ! rigid-lid
         IF( lk_mpp ) THEN
            c_solver_pt = 'G'   ! G= F with special staff ??? which one?
         ELSE
            c_solver_pt = 'F'
         ENDIF
      ELSE                          ! free surface T-point
         IF( lk_mpp ) THEN
            c_solver_pt = 'S'   ! S=T with special staff ??? which one?
         ELSE
            c_solver_pt = 'T'
         ENDIF
      ENDIF


      ! Construction of the elliptic system matrix
      ! ------------------------------------------

      CALL sol_mat( kt )


      IF( lk_isl ) THEN
      
         ! Islands in the domain
         ! ---------------------

         IF ( jpisl == 0 ) THEN
             IF(lwp)WRITE(numout,cform_err)
             IF(lwp)WRITE(numout,*) ' bad islands parameter jpisl =', jpisl
             nstop = nstop + 1
         ENDIF

         ! open Island streamfunction statistic file
         CALL ctlopn( numisp, 'islands.stat', 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL',   &
            &         1     , numout        , lwp      , 1                         )
   
         CALL isl_dom       ! Island identification

         CALL isl_bsf       ! Island barotropic stream function

         CALL isl_mat       ! Comput and invert the island matrix

         ! mbathy set to the number of w-level (minimum value 2)
         DO jj = 1, jpj
            DO ji = 1, jpi
               mbathy(ji,jj) = MAX( 1, mbathy(ji,jj) ) + 1
            END DO
         END DO

      ENDIF

   END SUBROUTINE solver_init

   !!======================================================================
END MODULE solver
