MODULE domain
   !!==============================================================================
   !!                       ***  MODULE domain   ***
   !! Ocean initialization : domain initialization
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!   dom_init       : initialize the space and time domain
   !!   dom_nam        : read and contral domain namelists
   !!   dom_ctl        : control print for the ocean domain
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce             ! 
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager
   USE daymod          ! calendar
   USE lib_mpp         ! distributed memory computing library
   USE flxrnf          ! runoffs

   USE domstp          ! domain: set the time-step
   USE domrea          ! domain: write the meshmask file
   USE dommsk          ! domain : mask

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC dom_init       ! called by opa.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL  (2005)
   !!   $Header: /home/opalod/NEMOCVSROOT/NEMO/OFF_SRC/DOM/domain.F90,v 1.2 2005/11/16 16:12:12 opalod Exp $
   !!   This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dom_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dom_init  ***
      !!                    
      !! ** Purpose :   Domain initialization. Call the routines that are 
      !!      required to create the arrays which define the space and time
      !!      domain of the ocean model.
      !!
      !! ** Method  :
      !!      - dom_stp: defined the model time step
      !!      - dom_rea: read the meshmask file if nmsh=1
      !!
      !! History :
      !!        !  90-10  (C. Levy - G. Madec)  Original code
      !!        !  91-11  (G. Madec)
      !!        !  92-01  (M. Imbard) insert time step initialization
      !!        !  96-06  (G. Madec) generalized vertical coordinate 
      !!        !  97-02  (G. Madec) creation of domwri.F
      !!        !  01-05  (E.Durand - G. Madec) insert closed sea
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER ::   iconf = 0         ! temporary integers
      !!----------------------------------------------------------------------

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dom_init : domain initialization'
         WRITE(numout,*) '~~~~~~~~'
      ENDIF

      CALL dom_nam                        ! read namelist ( namrun, namdom, namcla )

      CALL dom_stp                        ! Time step

      CALL dom_rea      ! Create a domain file

      CALL dom_msk      ! Masks

      CALL dom_ctl    ! Domain control

   END SUBROUTINE dom_init


   SUBROUTINE dom_nam
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dom_nam  ***
      !!                    
      !! ** Purpose :   read domaine namelists and print the variables.
      !!
      !! ** input   : - namrun namelist
      !!              - namdom namelist
      !!              - namcla namelist
      !!
      !! History :
      !!   9.0  !  03-08  (G. Madec)  Original code
      !!----------------------------------------------------------------------
      !! * Modules used
      USE ioipsl
      NAMELIST/namrun/ no    , cexper   , ln_rstart , nrstdt , nit000,          &
         &             nitend, ndate0   , nleapy   , ninist , nstock,           &
         &             nprint, nwrite   , nrunoff  , ln_ctl , nictls, nictle,   &
         &             njctls, njctle   , nbench   , isplt  , jsplt

      NAMELIST/namdom/ e3zps_min, e3zps_rat, nmsh  ,   &
         &             nacc  , atfp     , rdt      , rdtmin , rdtmax,   &
         &             rdth 

      NAMELIST/namcla/ n_cla
      !!----------------------------------------------------------------------

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dom_nam  : domain initialization through namelist read'
         WRITE(numout,*) '~~~~~~~ '
      ENDIF

      ! Namelist namrun : parameters of the run
      REWIND( numnam )
      READ  ( numnam, namrun )

      IF(lwp) THEN
         WRITE(numout,*) '        Namelist namrun'
         WRITE(numout,*) '           job number                      no        = ', no
         WRITE(numout,*) '           experiment name for output      cexper    = ', cexper
         WRITE(numout,*) '           restart logical                 ln_rstart = ', ln_rstart
         WRITE(numout,*) '           control of time step            nrstdt    = ', nrstdt
         WRITE(numout,*) '           number of the first time step   nit000    = ', nit000
         WRITE(numout,*) '           number of the last time step    nitend    = ', nitend
         WRITE(numout,*) '           initial calendar date aammjj    ndate0    = ', ndate0
         WRITE(numout,*) '           leap year calendar (0/1)        nleapy    = ', nleapy
         WRITE(numout,*) '           initial state output            ninist    = ', ninist
         WRITE(numout,*) '           level of print                  nprint    = ', nprint
         WRITE(numout,*) '           frequency of restart file       nstock    = ', nstock
         WRITE(numout,*) '           frequency of output file        nwrite    = ', nwrite
         WRITE(numout,*) '           runoff option                   nrunoff   = ', nrunoff
         WRITE(numout,*) '           run control (for debugging)     ln_ctl    = ', ln_ctl
         WRITE(numout,*) '           Start i indice for SUM control  nictls    = ', nictls
         WRITE(numout,*) '           End i indice for SUM control    nictle    = ', nictle
         WRITE(numout,*) '           Start j indice for SUM control  njctls    = ', njctls
         WRITE(numout,*) '           End j indice for SUM control    njctle    = ', njctle
         WRITE(numout,*) '           number of proc. following i     isplt     = ', isplt
         WRITE(numout,*) '           number of proc. following j     jsplt     = ', jsplt
         WRITE(numout,*) '           benchmark parameter (0/1)       nbench    = ', nbench
      ENDIF

      ndastp = ndate0                ! Assign initial date to current date

! ... Control the sub-domain area indices for the print control
      IF(ln_ctl)   THEN
         IF( lk_mpp ) THEN
            ! the domain is forced to the real splitted domain in MPI
            isplt = jpni ; jsplt = jpnj ; ijsplt = jpni*jpnj
         ELSE
            IF( isplt == 1 .AND. jsplt == 1  ) THEN
               IF(lwp) WRITE(numout,cform_war)
               IF(lwp) WRITE(numout,*)'          - isplt & jsplt are equal to 1'
               IF(lwp) WRITE(numout,*)'          - the print control will be done over the whole domain'
               IF(lwp) WRITE(numout,*)
            ENDIF

            ! compute the total number of processors ijsplt
            ijsplt = isplt*jsplt
         ENDIF

         IF(lwp) WRITE(numout,*)'          - The total number of processors over which the'
         IF(lwp) WRITE(numout,*)'            print control will be done is ijsplt : ', ijsplt

         ! Control the indices used for the SUM control
         IF( nictls+nictle+njctls+njctle == 0 )   THEN
            ! the print control is done over the default area
            lsp_area = .FALSE.
         ELSE
            ! the print control is done over a specific  area
            lsp_area = .TRUE.
            IF( nictls < 1 .OR. nictls > jpiglo )   THEN
               IF(lwp) WRITE(numout,cform_war)
               IF(lwp) WRITE(numout,*)'          - nictls must be 1<=nictls>=jpiglo, it is forced to 1'
               IF(lwp) WRITE(numout,*)
               nwarn = nwarn + 1
               nictls = 1
            ENDIF

            IF( nictle < 1 .OR. nictle > jpiglo )   THEN
               IF(lwp) WRITE(numout,cform_war)
               IF(lwp) WRITE(numout,*)'          - nictle must be 1<=nictle>=jpiglo, it is forced to jpiglo'
               IF(lwp) WRITE(numout,*)
               nwarn = nwarn + 1
               nictle = jpjglo
            ENDIF

            IF( njctls < 1 .OR. njctls > jpjglo )   THEN
               IF(lwp) WRITE(numout,cform_war)
               IF(lwp) WRITE(numout,*)'          - njctls must be 1<=njctls>=jpjglo, it is forced to 1'
               IF(lwp) WRITE(numout,*)
               nwarn = nwarn + 1
               njctls = 1
            ENDIF

            IF( njctle < 1 .OR. njctle > jpjglo )   THEN
               IF(lwp) WRITE(numout,cform_war)
               IF(lwp) WRITE(numout,*)'          - njctle must be 1<=njctle>= jpjglo, it is forced to jpjglo'
               IF(lwp) WRITE(numout,*)
               nwarn = nwarn + 1
               njctle = jpjglo
            ENDIF

         ENDIF          ! IF( nictls+nictle+njctls+njctle == 0 )
       ENDIF            ! IF(ln_ctl)

! ... Control of output frequency
      IF ( nstock == 0 ) THEN
          IF(lwp)WRITE(numout,cform_war)
          IF(lwp)WRITE(numout,*) '           nstock = ', nstock, ' it is forced to ', nitend
          nstock = nitend
          nwarn = nwarn + 1
      ENDIF
      IF ( nwrite == 0 ) THEN
          IF(lwp)WRITE(numout,cform_war)
          IF(lwp)WRITE(numout,*) '           nwrite = ', nwrite, ' it is forced to ', nitend
          nwrite = nitend
          nwarn = nwarn + 1
      ENDIF

      SELECT CASE ( nleapy )   ! Choose calendar for IOIPSL
      CASE (  1 ) 
         CALL ioconf_calendar('gregorian')
         IF(lwp) WRITE(numout,*) '           The IOIPSL calendar is "gregorian", i.e. leap year'
      CASE (  0 )
         CALL ioconf_calendar('noleap')
         IF(lwp) WRITE(numout,*) '           The IOIPSL calendar is "noleap", i.e. no leap year'
      CASE ( 30 )
         CALL ioconf_calendar('360d')
         IF(lwp) WRITE(numout,*) '           The IOIPSL calendar is "360d", i.e. 360 days in a year'
      END SELECT

      SELECT CASE ( nleapy )   ! year=raajj*days day=rjjhh*hours hour=rhhmm*minutes etc ...
      CASE ( 1 )
         raajj = 365.25
         raass = raajj * rjjss
         rmoss = raass/raamo
      CASE ( 0 )
         raajj = 365.
         raass = raajj * rjjss
         rmoss = raass/raamo
      CASE DEFAULT
         raajj = FLOAT( nleapy ) * raamo
         raass =        raajj    * rjjss
         rmoss = FLOAT( nleapy ) * rjjss
      END SELECT
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) '           nb of days per year      raajj = ', raajj,' days'
         WRITE(numout,*) '           nb of seconds per year   raass = ', raass, ' s'
         WRITE(numout,*) '           nb of seconds per month  rmoss = ', rmoss, ' s'
      ENDIF

      ! Namelist namdom : space/time domain (bathymetry, mesh, timestep)
      REWIND( numnam )
      READ  ( numnam, namdom )

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) '        Namelist namdom'
         WRITE(numout,*) '           minimum thickness of partial   e3zps_min = ', e3zps_min, ' (m)'
         WRITE(numout,*) '              step level                  e3zps_rat = ', e3zps_rat
         WRITE(numout,*) '           flag write mesh/mask file(s)   nmsh      = ', nmsh
         WRITE(numout,*) '                = 0   no file created                 '
         WRITE(numout,*) '                = 1   mesh_mask                       '
         WRITE(numout,*) '                = 2   mesh and mask                   '
         WRITE(numout,*) '                = 3   mesh_hgr, msh_zgr and mask      '
         WRITE(numout,*) '           acceleration of converge       nacc      = ', nacc
         WRITE(numout,*) '           asselin time filter parameter  atfp      = ', atfp
         WRITE(numout,*) '           time step                      rdt       = ', rdt
         WRITE(numout,*) '           minimum time step on tracers   rdtmin    = ', rdtmin
         WRITE(numout,*) '           maximum time step on tracers   rdtmax    = ', rdtmax
         WRITE(numout,*) '           depth variation tracer step    rdth      = ', rdth
      ENDIF



      ! Default values
      n_cla = 0

      ! Namelist cross land advection
      REWIND( numnam )
      READ  ( numnam, namcla )
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) '        Namelist namcla'
         WRITE(numout,*) '           cross land advection        n_cla        = ',n_cla
      ENDIF

   END SUBROUTINE dom_nam


   SUBROUTINE dom_ctl
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dom_ctl  ***
      !!
      !! ** Purpose :   Domain control.
      !!
      !! ** Method  :   compute and print extrema of masked scale factors
      !!
      !! History :
      !!   8.5  !  02-08  (G. Madec)    Original code
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER ::   iimi1, ijmi1, iimi2, ijmi2, iima1, ijma1, iima2, ijma2
      INTEGER, DIMENSION(2) ::   iloc      ! 
      REAL(wp) ::   ze1min, ze1max, ze2min, ze2max
      !!----------------------------------------------------------------------

      ! Extrema of the scale factors

      IF(lwp)WRITE(numout,*)
      IF(lwp)WRITE(numout,*) 'dom_ctl : extrema of the masked scale factors'
      IF(lwp)WRITE(numout,*) '~~~~~~~'

      IF (lk_mpp) THEN
         CALL mpp_minloc( e1t(:,:), tmask(:,:,1), ze1min, iimi1,ijmi1 )
         CALL mpp_minloc( e2t(:,:), tmask(:,:,1), ze2min, iimi2,ijmi2 )
         CALL mpp_maxloc( e1t(:,:), tmask(:,:,1), ze1max, iima1,ijma1 )
         CALL mpp_maxloc( e2t(:,:), tmask(:,:,1), ze2max, iima2,ijma2 )
      ELSE
         ze1min = MINVAL( e1t(:,:), mask = tmask(:,:,1) == 1.e0 )    
         ze2min = MINVAL( e2t(:,:), mask = tmask(:,:,1) == 1.e0 )    
         ze1max = MAXVAL( e1t(:,:), mask = tmask(:,:,1) == 1.e0 )    
         ze2max = MAXVAL( e2t(:,:), mask = tmask(:,:,1) == 1.e0 )    

         iloc  = MINLOC( e1t(:,:), mask = tmask(:,:,1) == 1.e0 )
         iimi1 = iloc(1) + nimpp - 1
         ijmi1 = iloc(2) + njmpp - 1
         iloc  = MINLOC( e2t(:,:), mask = tmask(:,:,1) == 1.e0 )
         iimi2 = iloc(1) + nimpp - 1
         ijmi2 = iloc(2) + njmpp - 1
         iloc  = MAXLOC( e1t(:,:), mask = tmask(:,:,1) == 1.e0 )
         iima1 = iloc(1) + nimpp - 1
         ijma1 = iloc(2) + njmpp - 1
         iloc  = MAXLOC( e2t(:,:), mask = tmask(:,:,1) == 1.e0 )
         iima2 = iloc(1) + nimpp - 1
         ijma2 = iloc(2) + njmpp - 1
      ENDIF

      IF(lwp) THEN
         WRITE(numout,"(14x,'e1t maxi: ',1f10.2,' at i = ',i5,' j= ',i5)") ze1max, iima1, ijma1
         WRITE(numout,"(14x,'e1t mini: ',1f10.2,' at i = ',i5,' j= ',i5)") ze1min, iimi1, ijmi1
         WRITE(numout,"(14x,'e2t maxi: ',1f10.2,' at i = ',i5,' j= ',i5)") ze2max, iima2, ijma2
         WRITE(numout,"(14x,'e2t mini: ',1f10.2,' at i = ',i5,' j= ',i5)") ze2min, iimi2, ijmi2
      ENDIF

   END SUBROUTINE dom_ctl

   !!======================================================================
END MODULE domain
