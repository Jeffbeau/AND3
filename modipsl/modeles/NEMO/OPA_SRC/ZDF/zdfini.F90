MODULE zdfini
   !!======================================================================
   !!              ***  MODULE  zdfini  ***
   !! Ocean physics : define vertical mixing variables
   !!=====================================================================

   !!----------------------------------------------------------------------
   !!   zdf_init    : initialization, namelist read, and parameters control
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_oce         ! mesh and scale factors
   USE ldftra_oce      ! ocean active tracers: lateral physics
   USE ldfdyn_oce      ! ocean dynamics lateral physics
   USE zdf_oce         ! TKE vertical mixing          
   USE lib_mpp         ! distribued memory computing
   USE zdftke          ! TKE vertical mixing  
   USE zdfkpp          ! KPP vertical mixing          
   USE zdfddm          ! double diffusion mixing      
   USE zdfevd          ! enhanced vertical diffusion  
   USE zdfric          ! Richardson vertical mixing   
   USE tranpc          ! convection: non penetrative adjustment
   USE ldfslp          ! iso-neutral slopes

   USE in_out_manager  ! I/O manager

   IMPLICIT NONE
   PRIVATE

   !! *  Routine accessibility
   PUBLIC zdf_init          ! routine called by opa.F90
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/ZDF/zdfini.F90,v 1.4 2005/09/22 10:57:33 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------
   
CONTAINS

   SUBROUTINE zdf_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_init  ***
      !! 
      !! ** Purpose :   initializations of the vertical ocean physics
      !!
      !! ** Method  :   Read namelist namzdf, control logicals 
      !!
      !! History :
      !!        !  97-06  (G. Madec)  Original code from inimix
      !!   8.5  !  02-08  (G. Madec)  F90 : free form
      !!   9.0  !  05-06  (C. Ethe) KPP parameterization
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER ::   ioptio       ! temporary scalar
      LOGICAL ::           & !!! namzdf: vertical diffusion
         ln_zdfexp = .FALSE.   ! explicit vertical diffusion scheme flag

      !! * Namelist
      NAMELIST/namzdf/ ln_zdfevd, ln_zdfnpc,   &
         &             avm0, avt0, avevd, nevdm, ln_zdfexp, n_zdfexp
      !!----------------------------------------------------------------------
      !!  OPA 9.0, LODYC-IPSL (2003)
      !!----------------------------------------------------------------------

      ! Read namzdf namelist : vertical mixing parameters
      ! --------------------
      REWIND( numnam )
      READ  ( numnam, namzdf )

      ! Parameter print
      ! ---------------
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'zdf_init: vertical physics'
         WRITE(numout,*) '~~~~~~~~'
         WRITE(numout,*) '          Namelist namzdf : set vertical mixing mixing parameters'
         WRITE(numout,*) '             enhanced vertical diffusion      ln_zdfevd = ', ln_zdfevd
         WRITE(numout,*) '             non-penetrative convection       ln_zdfnpc = ', ln_zdfnpc
         WRITE(numout,*) '             vertical eddy viscosity             avm0   = ', avm0
         WRITE(numout,*) '             vertical eddy diffusivity           avt0   = ', avt0
         WRITE(numout,*) '             vertical coefficient for evd        avevd  = ', avevd
         WRITE(numout,*) '                applied on momentum (=1/0)       nevdm  = ', nevdm
         WRITE(numout,*) '             time splitting / backward scheme ln_zdfexp = ', ln_zdfexp
         WRITE(numout,*) '             number of time step               n_zdfexp = ', n_zdfexp
      ENDIF

      ! Parameter & logicals controls
      ! -----------------------------
      ! ... vertical mixing
      ! time stepping scheme (N.B. TKE or KPP schemes => force the use of implicit scheme)
      IF( ( ln_zdfexp .AND. .NOT.lk_zdftke ) .OR. ( ln_zdfexp .AND. .NOT.lk_zdfkpp ) ) THEN  
         l_trazdf_exp = .TRUE.           ! use explicit scheme
         l_trazdf_imp = .FALSE.
         l_dynzdf_exp = .TRUE.           ! use explicit scheme
         l_dynzdf_imp = .FALSE.
      ELSE
         l_trazdf_exp = .FALSE.          ! use implicit scheme
         l_trazdf_imp = .TRUE. 
         l_dynzdf_exp = .FALSE.          ! use implicit scheme
         l_dynzdf_imp = .TRUE. 
      ENDIF
      IF( l_trazdf_iso .OR. l_trazdf_iso_vo ) THEN  
         l_trazdf_exp = .FALSE.          ! iso-neutral diffusion : 
         l_trazdf_imp = .FALSE.          ! implicit scheme included in iso-neutral routine
      ENDIF
      IF( l_dynldf_iso ) THEN  
         l_dynzdf_exp = .FALSE.          ! iso-neutral diffusion :
         l_dynzdf_imp = .FALSE.          ! implicit scheme included in iso-neutral routine
      ENDIF
#if defined key_autotasking
      IF( l_dynzdf_imp ) THEN
         l_dynzdf_imp     = .FALSE.
         l_dynzdf_imp_tsk = .TRUE.
      ENDIF
#else
      l_dynzdf_imp_tsk = .FALSE.
#endif
      IF( lk_esopa  ) THEN
         l_trazdf_exp = .TRUE.           ! esopa: use all options
         l_trazdf_imp = .TRUE.
         l_dynzdf_exp     = .TRUE.           ! esopa: use all options
         l_dynzdf_imp     = .TRUE.
         l_dynzdf_imp_tsk = .TRUE.
      ENDIF

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '          vertical mixing option :'
      ioptio = 0
      IF( lk_zdfcst ) THEN
         IF(lwp) WRITE(numout,*) '             constant eddy diffusion coef.'
         ioptio = ioptio+1
      ENDIF
      IF( lk_zdfric ) THEN
         IF(lwp) WRITE(numout,*) '             Richardson dependent eddy coef.'
         ioptio = ioptio+1
      ENDIF
      IF( lk_zdftke ) THEN
         IF(lwp) WRITE(numout,*) '             TKE dependent eddy coef.'
         ioptio = ioptio+1
      ENDIF
      IF( lk_zdfkpp ) THEN
         IF(lwp) WRITE(numout,*) '             KPP dependent eddy coef.'
         ioptio = ioptio+1
      ENDIF
      IF( ioptio == 0 .OR. ioptio > 1 .AND. .NOT. lk_esopa ) THEN
          IF(lwp) WRITE(numout,cform_err)
          IF(lwp) WRITE(numout,*) ' one and only one vertical diffusion option has to be defined '
          nstop = nstop + 1
      ENDIF

      ! ... Convection
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '          convection :'
      ioptio = 0
      IF( ln_zdfnpc ) THEN
         IF(lwp) WRITE(numout,*) '             use non penetrative convective scheme'
         ioptio = ioptio+1
      ENDIF
      IF( ln_zdfevd ) THEN
         IF(lwp) WRITE(numout,*) '             use enhanced vertical dif. scheme'
         ioptio = ioptio+1
      ENDIF
      IF( lk_zdftke ) THEN
         IF(lwp) WRITE(numout,*) '             use the 1.5 turbulent closure'
      ENDIF
      IF( lk_zdfkpp ) THEN
         IF(lwp) WRITE(numout,*) '             use the KPP closure scheme'
         IF(lk_mpp) THEN
            IF(lwp) WRITE(numout,cform_err)
            IF(lwp) WRITE(numout,*) '             The KPP scheme is not ready to run in MPI'
         ENDIF
      ENDIF
      IF ( ioptio > 1 .AND. .NOT. lk_esopa ) THEN
          IF(lwp) WRITE(numout,cform_err)
          IF(lwp) WRITE(numout,*) ' chose between ln_zdfnpc'
          IF(lwp) WRITE(numout,*) '           and ln_zdfevd'
          nstop = nstop + 1
      ENDIF
      IF( ioptio == 0 .AND. .NOT. lk_zdftke ) THEN
         IF(lwp) WRITE(numout,cform_err)
         IF(lwp) WRITE(numout,*) ' except for TKE scheme, a convection scheme is'
         IF(lwp) WRITE(numout,*) ' required: ln_zdfevd or ln_zdfnpc logicals'
         nstop = nstop + 1
      ENDIF

   END SUBROUTINE zdf_init

   !!======================================================================
END MODULE zdfini
