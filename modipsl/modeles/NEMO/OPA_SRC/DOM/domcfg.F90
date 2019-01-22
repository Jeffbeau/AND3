MODULE domcfg
   !!==============================================================================
   !!                       ***  MODULE domcfg   ***
   !! Ocean initialization : domain configuration initialization
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!   dom_cfg        : initialize the domain configuration
   !!----------------------------------------------------------------------
   !! * Modules used
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distributed memory computing library
   USE solisl          ! ???

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC dom_cfg        ! called by opa.F90
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DOM/domcfg.F90,v 1.4 2006/04/10 15:46:07 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dom_cfg
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dom_cfg  ***
      !!                    
      !! ** Purpose :   set the domain configuration
      !!
      !! ** Method  :
      !!
      !! History :
      !!   9.0  !  03-09  (G. Madec)  Original code
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER ::   iconf = 0         ! temporary integers
      !!----------------------------------------------------------------------

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dom_cfg : set the ocean configuration'
         WRITE(numout,*) '~~~~~~~      ocean model configuration used :',   &
            &                             ' cp_cfg = ', cp_cfg, ' jp_cfg = ', jp_cfg
      ENDIF

      ! Global domain boundary conditions
      ! ---------------------------------
      IF(lwp) THEN
         WRITE(numout,*) '          global domain lateral boundaries'

         IF( jperio == 0 ) WRITE(numout,*) '             jperio= 0, closed'
         IF( jperio == 1 ) WRITE(numout,*) '             jperio= 1, cyclic east-west'
         IF( jperio == 2 ) WRITE(numout,*) '             jperio= 2, equatorial symmetric'
         IF( jperio == 3 ) WRITE(numout,*) '             jperio= 3, north fold with T-point pivot'
         IF( jperio == 4 ) WRITE(numout,*) '             jperio= 4, cyclic east-west and',   &
                                                                  ' north fold with T-point pivot'
         IF( jperio == 5 ) WRITE(numout,*) '             jperio= 5, north fold with F-point pivot'
         IF( jperio == 6 ) WRITE(numout,*) '             jperio= 6, cyclic east-west and',   &
                                                                  ' north fold with F-point pivot'
      ENDIF
      IF( jperio <  0 .OR. jperio > 6 ) THEN
          IF(lwp) WRITE(numout,cform_err)
          IF(lwp) WRITE(numout,*) 'jperio is out of range'
          nstop = nstop + 1
      ENDIF


      ! global domain versus zoom and/or local domain
      ! ---------------------------------------------

      CALL dom_glo 

   END SUBROUTINE dom_cfg


   SUBROUTINE dom_glo
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dom_glo  ***
      !!
      !! ** Purpose :   initialization for global domain, zoom and local domain
      !!
      !! ** Method  :   
      !!
      !! ** Action  : - mig  , mjg : 
      !!              - mi0  , mi1   :
      !!              - mj0, , mj1   :
      !!
      !! History :
      !!   8.5  !  02-08  (G. Madec)    Original code
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER ::   ji, jj            ! dummy loop argument
      !!----------------------------------------------------------------------

      ! Local domain 
      ! ============

      ! local domain indices ==> data domain indices
      DO ji = 1, jpi
        mig(ji) = ji + jpizoom - 1 + nimpp - 1
      END DO
      DO jj = 1, jpj
        mjg(jj) = jj + jpjzoom - 1 + njmpp - 1
      END DO

      ! data domain indices ==> local domain indices
      ! (return (m.0,m.1)=(1,0) if data domain gridpoint is to the west/south of the 
      ! local domain, or (m.0,m.1)=(jp.+1,jp.) to the east/north of local domain. 
      DO ji = 1, jpidta
        mi0(ji) = MAX( 1, MIN( ji - jpizoom + 1 - nimpp + 1, jpi+1 ) )
        mi1(ji) = MAX( 0, MIN( ji - jpizoom + 1 - nimpp + 1, jpi   ) )
      END DO
      DO jj = 1, jpjdta
        mj0(jj) = MAX( 1, MIN( jj - jpjzoom + 1 - njmpp + 1, jpj+1 ) )
        mj1(jj) = MAX( 0, MIN( jj - jpjzoom + 1 - njmpp + 1, jpj   ) )
      END DO

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dom_glo : domain: data / local '
         WRITE(numout,*) '~~~~~~~ '
         WRITE(numout,*) '          data input domain    : jpidta = ', jpidta,   &
            &                                            ' jpjdta = ', jpjdta, ' jpkdta = ', jpkdta
         WRITE(numout,*) '          global or zoom domain: jpiglo = ', jpiglo,   &
            &                                            ' jpjglo = ', jpjglo, ' jpk    = ', jpk
         WRITE(numout,*) '          local domain         : jpi    = ', jpi   ,   &
            &                                            ' jpj    = ', jpj   , ' jpk    = ', jpk
         WRITE(numout,*)
         WRITE(numout,*) '          south-west indices    jpizoom = ', jpizoom,   &
            &                                           ' jpjzoom = ', jpjzoom
         WRITE(numout,*)
         WRITE(numout,*) '          conversion local  ==> data i-index domain'
         WRITE(numout,25)              (mig(ji),ji = 1,jpi)
         WRITE(numout,*)
         WRITE(numout,*) '          conversion data   ==> local  i-index domain'
         WRITE(numout,*) '             starting index'
         WRITE(numout,25)              (mi0(ji),ji = 1,jpidta)
         WRITE(numout,*) '             ending index'
         WRITE(numout,25)              (mi1(ji),ji = 1,jpidta)
         WRITE(numout,*)
         WRITE(numout,*) '          conversion local  ==> data j-index domain'
         WRITE(numout,25)              (mjg(jj),jj = 1,jpj)
         WRITE(numout,*)
         WRITE(numout,*) '          conversion data  ==> local j-index domain'
         WRITE(numout,*) '             starting index'
         WRITE(numout,25)              (mj0(jj),jj = 1,jpjdta)
         WRITE(numout,*) '             ending index'
         WRITE(numout,25)              (mj1(jj),jj = 1,jpjdta)
      ENDIF
 25   FORMAT( 100(10x,19i4,/) )

      ! Zoom domain
      ! ===========

      ! zoom control
      IF( jpiglo + jpizoom - 1  >  jpidta .OR.   &
          jpjglo + jpjzoom - 1  >  jpjdta      ) THEN
         IF(lwp)WRITE(numout,cform_err)
         IF(lwp)WRITE(numout,*)' global or zoom domain exceed the data domain ! '
         nstop = nstop + 1
      ENDIF

      ! set zoom flag
      IF ( jpiglo < jpidta .OR. jpjglo < jpjdta )   lzoom = .TRUE.

      ! set zoom type flags
      IF( lzoom .AND. jpizoom /= 1 )   lzoom_w = .TRUE.                     ! 
      IF( lzoom .AND. jpjzoom /= 1 )   lzoom_s = .TRUE.
      IF( lzoom .AND. jpiglo + jpizoom -1 /= jpidta )   lzoom_e = .TRUE.
      IF( lzoom .AND. jpjglo + jpjzoom -1 /= jpjdta )   lzoom_n = .TRUE.

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) '          zoom flags : '
         WRITE(numout,*) '             lzoom   = ', lzoom  , ' (T = zoom, F = global )'
         WRITE(numout,*) '             lzoom_e = ', lzoom_e, ' (T = forced closed east  boundary)'
         WRITE(numout,*) '             lzoom_w = ', lzoom_w, ' (T = forced closed west  boundary)'
         WRITE(numout,*) '             lzoom_s = ', lzoom_s, ' (T = forced closed South boundary)'
         WRITE(numout,*) '             lzoom_n = ', lzoom_n, ' (T = forced closed North boundary)'
      ENDIF
      IF(  ( lzoom_e .OR. lzoom_w )  .AND.  ( jperio == 1 .OR. jperio == 4 .OR. jperio == 6 )  ) THEN
         IF(lwp)WRITE(numout,cform_err)
         IF(lwp)WRITE(numout,*)' Your zoom choice is inconsistent with east-west cyclic boundary condition'
         nstop = nstop + 1
      ENDIF
      IF(  lzoom_n  .AND.  ( 3 <= jperio .AND. jperio <= 6 )  ) THEN
         IF(lwp)WRITE(numout,cform_err)
         IF(lwp)WRITE(numout,*)' Your zoom choice is inconsistent with North fold boundary condition'
         nstop = nstop + 1
      ENDIF
      IF(  lzoom  .AND.  lk_isl  ) THEN
         IF(lwp)WRITE(numout,cform_err)
         IF(lwp)WRITE(numout,*)' key_islands and zoom are not allowed'
         nstop = nstop + 1
      ENDIF

!!DB: ORCA-related
      ! Pre-defined arctic/antarctic zoom of ORCA configuration flag
!      IF( cp_cfg == "orca" ) THEN
!         SELECT CASE ( jp_cfg )
!         !                                        ! =======================
!         CASE ( 2 )                               !  ORCA_R2 configuration
!            !                                     ! =======================
!            IF(  jpiglo  == 142    .AND. jpjglo  ==  53 .AND.   &
!               & jpizoom ==  21    .AND. jpjzoom ==  97         )   lzoom_arct = .TRUE.
!            IF(  jpiglo  == jpidta .AND. jpjglo  ==  50 .AND.   &
!               & jpizoom ==   1    .AND. jpjzoom ==   1         )   lzoom_anta = .TRUE.
!            !                                     ! =======================
!         CASE ( 05 )                              !  ORCA_R05 configuration
!            !                                     ! =======================
!            IF(  jpiglo  == 562    .AND. jpjglo  == 202 .AND.   &
!               & jpizoom ==  81    .AND. jpjzoom == 301         )   lzoom_arct = .TRUE.
!            IF(  jpiglo  == jpidta .AND. jpjglo  == 187 .AND.   &
!               & jpizoom ==   1    .AND. jpjzoom ==   1         )   lzoom_anta = .TRUE.
!         END SELECT
!         !
!         IF(lwp) WRITE(numout,*) '          ORCA configuration: antarctic/arctic zoom flags : '
!         IF(lwp) WRITE(numout,*) '             lzoom_arct = ', lzoom_arct, ' (T=   arctic zoom, F=global)'
!         IF(lwp) WRITE(numout,*) '             lzoom_anta = ', lzoom_anta, ' (T=antarctic zoom, F=global)'
!         !
!      ENDIF
         
   END SUBROUTINE dom_glo

   !!======================================================================
END MODULE domcfg
