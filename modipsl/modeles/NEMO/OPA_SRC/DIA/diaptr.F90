MODULE diaptr
   !!======================================================================
   !!                       ***  MODULE  diaptr  ***
   !! Ocean physics:  brief description of the purpose of the module
   !!                 (please no more than 2 lines)
   !!=====================================================================
   !!----------------------------------------------------------------------
   !!   dia_ptr      : Poleward Transport Diagnostics module
   !!   dia_ptr_init : Initialization, namelist read
   !!   dia_ptr_wri  : Output of poleward fluxes
   !!   ptr_vjk      : "zonal" sum computation of a "meridional" flux array
   !!   ptr_vtjk     : "zonal" mean computation of a tracer field
   !!   ptr_vj       : "zonal" and vertical sum computation of a "meridional"
   !!                : flux array; Generic interface: ptr_vj_3d, ptr_vj_2d
   !!----------------------------------------------------------------------
   !! * Modules used
   USE oce           ! ocean dynamics and active tracers
   USE dom_oce       ! ocean space and time domain
   USE ldftra_oce    ! ???
   USE lib_mpp
   USE in_out_manager
   USE dianam
   USE phycst

   IMPLICIT NONE
   PRIVATE

   INTERFACE ptr_vj
      MODULE PROCEDURE ptr_vj_3d, ptr_vj_2d
   END INTERFACE

   !! *  Routine accessibility
   PUBLIC dia_ptr_init   ! call in opa module
   PUBLIC dia_ptr        ! call in step module
   PUBLIC ptr_vj         ! call by tra_ldf & tra_adv routines
   PUBLIC ptr_vjk        ! call by tra_ldf & tra_adv routines

   !! * Share Module variables
   LOGICAL, PUBLIC ::       & !!! ** init namelist (namptr) **
      ln_diaptr = .FALSE.,  &  !: Poleward transport flag (T) or not (F)
      ln_subbas = .FALSE.      !: Atlantic/Pacific/Indian basins calculation
   INTEGER, PUBLIC ::       & !!: ** ptr namelist (namptr) **
      nf_ptr = 15              !: frequency of ptr computation
   REAL(wp), PUBLIC, DIMENSION(jpj) ::   &   !!: poleward transport
      pht_adv, pst_adv,     &  !: heat and salt: advection
      pht_ove, pst_ove,     &  !: heat and salt: overturning
      pht_ldf, pst_ldf,     &  !: heat and salt: lateral diffusion
#if defined key_diaeiv
      pht_eiv, pst_eiv,     &  !: heat and salt: bolus advection
#endif
      ht_atl,ht_ind,ht_pac, &  !: heat
      st_atl,st_ind,st_pac     !: salt
   REAL(wp),DIMENSION(jpi,jpj) ::   &
      abasin,pbasin,ibasin     !: return function value
     

   !! Module variables
   REAL(wp), DIMENSION(jpj,jpk) ::   &  
      tn_jk  , sn_jk  ,  &  !: "zonal" mean temperature and salinity
      v_msf_atl       ,  &  !: "meridional" Stream-Function
      v_msf_glo       ,  &  !: "meridional" Stream-Function
      v_msf_ipc       ,  &  !: "meridional" Stream-Function
#if defined key_diaeiv
      v_msf_eiv       ,  &  !: bolus "meridional" Stream-Function
#endif
      surf_jk_r             !: inverse of the ocean "zonal" section surface

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/OPA_SRC/DIA/diaptr.F90,v 1.10 2006/03/20 16:05:53 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------

CONTAINS

   FUNCTION ptr_vj_3d( pva )   RESULT ( p_fval )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE ptr_vj_3d  ***
      !!
      !! ** Purpose :   "zonal" and vertical sum computation of a "meridional"
      !!      flux array
      !!
      !! ** Method  : - i-k sum of pva using the interior 2D vmask (vmask_i).
      !!      pva is supposed to be a masked flux (i.e. * vmask*e1v*e3v)
      !!
      !! ** Action  : - p_fval: i-k-mean poleward flux of pva
      !!
      !! History :
      !!   9.0  !  03-09  (G. Madec)  Original code
      !!----------------------------------------------------------------------
      !! * arguments
      REAL(wp) , INTENT(in), DIMENSION(jpi,jpj,jpk) ::   &
         pva                         ! mask flux array at V-point

      !! * local declarations
      INTEGER  ::   ji, jj, jk        ! dummy loop arguments
#if ! defined key_agrif
      INTEGER  ::   ijpj = jpj        ! ???
#else
      INTEGER  ::   ijpj             ! ???
#endif      
      REAL(wp),DIMENSION(jpj) ::   &
         p_fval                       ! function value
      !!--------------------------------------------------------------------
#if defined key_agrif
      ijpj = jpj
#endif      

      p_fval(:) = 0.e0
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! Vector opt.
               p_fval(jj) = p_fval(jj) + pva(ji,jj,jk) * tmask_i(ji,jj+1) * tmask_i(ji,jj) 
            END DO
         END DO
      END DO

      IF( lk_mpp )   CALL mpp_sum( p_fval, ijpj )     !!bug  I presume

   END FUNCTION ptr_vj_3d



   FUNCTION ptr_vj_2d( pva )   RESULT ( p_fval )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE ptr_vj_2d  ***
      !!
      !! ** Purpose :   "zonal" and vertical sum computation of a "meridional"
      !!      flux array
      !!
      !! ** Method  : - i-k sum of pva using the interior 2D vmask (vmask_i).
      !!      pva is supposed to be a masked flux (i.e. * vmask*e1v*e3v)
      !!
      !! ** Action  : - p_fval: i-k-mean poleward flux of pva
      !!
      !! History :
      !!   9.0  !  03-09  (G. Madec)  Original code
      !!----------------------------------------------------------------------
      !! * arguments
      REAL(wp) , INTENT(in), DIMENSION(jpi,jpj) ::   &
         pva                         ! mask flux array at V-point

      !! * local declarations
      INTEGER  ::   ji,jj             ! dummy loop arguments
#if ! defined key_agrif
      INTEGER  ::   ijpj = jpj        ! ???
#else
      INTEGER  ::   ijpj             ! ???
#endif      
      REAL(wp),DIMENSION(jpj) ::   &
         p_fval                       ! function value
      !!--------------------------------------------------------------------
#if defined key_agrif
      ijpj = jpj
#endif      
 
      p_fval(:) = 0.e0
      DO jj = 2, jpjm1
         DO ji = fs_2, fs_jpim1   ! Vector opt.
            p_fval(jj) = p_fval(jj) + pva(ji,jj) * tmask_i(ji,jj+1) * tmask_i(ji,jj)
         END DO
      END DO

      IF( lk_mpp )   CALL mpp_sum( p_fval, ijpj )     !!bug  I presume
 
    END FUNCTION ptr_vj_2d



   FUNCTION ptr_vjk( pva )   RESULT ( p_fval )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE ptr_vjk  ***
      !!
      !! ** Purpose :   "zonal" sum computation of a "meridional" flux array
      !!
      !! ** Method  : - i-sum of pva using the interior 2D vmask (vmask_i).
      !!      pva is supposed to be a masked flux (i.e. * vmask*e1v*e3v)
      !!
      !! ** Action  : - p_fval: i-k-mean poleward flux of pva
      !!
      !! History :
      !!   9.0  !  03-09  (G. Madec)  Original code
      !!----------------------------------------------------------------------
      !! * arguments
      REAL(wp) , INTENT(in), DIMENSION(jpi,jpj,jpk) ::   &
         pva                         ! mask flux array at V-point

      !! * local declarations
      INTEGER  ::   ji, jj, jk        ! dummy loop arguments
      INTEGER, DIMENSION (1) :: ish
      INTEGER, DIMENSION (2) :: ish2
      REAL(wp),DIMENSION(jpj*jpk) ::   &
         zwork                        ! temporary vector for mpp_sum
      REAL(wp),DIMENSION(jpj,jpk) ::   &
         p_fval                       ! return function value
      !!--------------------------------------------------------------------
 
      p_fval(:,:) = 0.e0

      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
           DO ji = fs_2, fs_jpim1
            p_fval(jj,jk) = p_fval(jj,jk) + pva(ji,jj,jk) * e1v(ji,jj) * fse3v(ji,jj,jk) &
               &            * tmask_i(ji,jj+1) * tmask_i(ji,jj)
           END DO
         END DO
      END DO

      IF(lk_mpp)   THEN
         ish(1) = jpj*jpk ; ish2(1)=jpj ; ish2(2)=jpk
         zwork(:)= RESHAPE(p_fval, ish )
         CALL mpp_sum(zwork, jpj*jpk )
         p_fval(:,:)= RESHAPE(zwork,ish2)
      END IF

   END FUNCTION ptr_vjk

   FUNCTION ptr_vtjk( pva )   RESULT ( p_fval )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE ptr_vtjk  ***
      !!
      !! ** Purpose :   "zonal" mean computation of a tracer field
      !!
      !! ** Method  : - i-sum of mj(pva) using the interior 2D vmask (vmask_i)
      !!      multiplied by the inverse of the surface of the "zonal" ocean
      !!      section
      !!
      !! ** Action  : - p_fval: i-k-mean poleward flux of pva
      !!
      !! History :
      !!   9.0  !  03-09  (G. Madec)  Original code
      !!   9.0  !  06-01  (A. Biastoch)  Allow sub-basins computation
      !!----------------------------------------------------------------------
      !! * arguments
      REAL(wp) , INTENT(in), DIMENSION(jpi,jpj,jpk) ::   &
         pva                         ! mask flux array at V-point
 
      !! * local declarations
      INTEGER  ::   ji, jj, jk        ! dummy loop arguments
      INTEGER, DIMENSION (1) :: ish
      INTEGER, DIMENSION (2) :: ish2
      REAL(wp),DIMENSION(jpj*jpk) ::   &
         zwork                        ! temporary vector for mpp_sum
      REAL(wp),DIMENSION(jpj,jpk) ::   &
         p_fval                       ! return function value
      !!-------------------------------------------------------------------- 

      p_fval(:,:) = 0.e0
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = fs_2, fs_jpim1   ! Vector opt.
               p_fval(jj,jk) = p_fval(jj,jk) + ( pva(ji,jj,jk) + pva(ji,jj+1,jk) )              &
                  &                          * e1v(ji,jj) * fse3v(ji,jj,jk) * vmask(ji,jj,jk)   &
                  &                          * tmask_i(ji,jj+1) * tmask_i(ji,jj)
            END DO
         END DO
      END DO
      p_fval(:,:) = p_fval(:,:) * 0.5
      IF(lk_mpp)   THEN
         ish(1) = jpj*jpk ; ish2(1)=jpj ; ish2(2)=jpk
         zwork(:)= RESHAPE(p_fval, ish )
         CALL mpp_sum(zwork, jpj*jpk )
         p_fval(:,:)= RESHAPE(zwork,ish2)
      END IF

   END FUNCTION ptr_vtjk


   SUBROUTINE dia_ptr( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dia_ptr  ***
      !!----------------------------------------------------------------------
      !! * Moudules used
      USE ioipsl

      !! * Argument
      INTEGER, INTENT(in) ::   kt   ! ocean time step index

      !! * Local variables
      INTEGER ::   jk,jj,ji               ! dummy loop
      REAL(wp) ::    &
         zsverdrup,  &              ! conversion from m3/s to Sverdrup
         zpwatt,     &              ! conversion from W    to PW
         zggram                     ! conversion from g    to Pg

      REAL(wp), DIMENSION(jpi,jpj,jpk) ::  &
         v_atl , v_ipc,                    &
         vt_atl, vt_pac, vt_ind,           &
         vs_atl, vs_pac, vs_ind,           &
         zv_eiv
      CHARACTER (len=32) ::   &
         clnam = 'subbasins.nc'                
      INTEGER ::  itime,inum,ipi,ipj,ipk       ! temporary integer
      INTEGER, DIMENSION (1) ::   istep
      REAL(wp) ::    zdate0,zsecond,zdt        ! temporary scalars
      REAL(wp), DIMENSION(jpidta,jpjdta) ::   &
         zlamt, zphit, zdta             ! temporary workspace (NetCDF read)
      REAL(wp), DIMENSION(jpk) ::   &
         zdept                          ! temporary workspace (NetCDF read)
      !!----------------------------------------------------------------------

      IF( kt == nit000 .OR. MOD( kt, nf_ptr ) == 0 )   THEN

         zsverdrup = 1.e-6
         zpwatt    = 1.e-15
         zggram    = 1.e-6
         ipi       = jpidta
         ipj       = jpjdta
         ipk       = 1
         itime     = 1
         zsecond   = 0.e0
         zdate0    = 0.e0
   
# if defined key_diaeiv
         zv_eiv(:,:,:) = v_eiv(:,:,:)
# else
         zv_eiv(:,:,:) = 0.e0
# endif

         ! "zonal" mean temperature and salinity at V-points
         tn_jk(:,:) = ptr_vtjk( tn(:,:,:) ) * surf_jk_r(:,:)
         sn_jk(:,:) = ptr_vtjk( sn(:,:,:) ) * surf_jk_r(:,:)

         !--------------------------------------------------------
         ! overturning calculation:
 
         IF( ln_subbas ) THEN              ! Basins computation

            IF( kt == nit000 ) THEN                ! load basin mask
               itime = 1
               ipi   = jpidta
               ipj   = jpjdta
               ipk   = 1
               zdt   = 0.e0
               istep = 0
               clnam = 'subbasins.nc'

               CALL flinopen(clnam,1,jpidta,1,jpjdta,.FALSE.,ipi,ipj, &
                  &          ipk,zlamt,zphit,zdept,itime,istep,zdate0,zdt,inum)

               ! get basins:
               abasin (:,:) = 0.e0
               pbasin (:,:) = 0.e0
               ibasin (:,:) = 0.e0

               ! Atlantic basin
               CALL flinget(inum,'atlmsk',jpidta,jpjdta,1,itime,1,   &
                  &         0,1,jpidta,1,jpjdta,zdta(:,:))
               DO jj = 1, nlcj                                 ! interior values
                  DO ji = 1, nlci
                     abasin (ji,jj) = zdta( mig(ji), mjg(jj) )
                  END DO
               END DO

               ! Pacific basin
               CALL flinget(inum,'pacmsk',jpidta,jpjdta,1,itime,1,   &
                  &         0,1,jpidta,1,jpjdta,zdta(:,:))
               DO jj = 1, nlcj                                 ! interior values
                  DO ji = 1, nlci
                     pbasin (ji,jj) = zdta( mig(ji), mjg(jj) )
                  END DO
               END DO

               ! Indian basin
               CALL flinget(inum,'indmsk',jpidta,jpjdta,1,itime,1,   &
                  &         0,1,jpidta,1,jpjdta,zdta(:,:))
               DO jj = 1, nlcj                                 ! interior values
                  DO ji = 1, nlci
                     ibasin (ji,jj) = zdta( mig(ji), mjg(jj) )
                  END DO
               END DO

               CALL flinclo(inum)

            ENDIF

            ! basin separation:
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ! basin separated velocity
                  v_atl(ji,jj,:) = (vn(ji,jj,:)+zv_eiv(ji,jj,:))*abasin(ji,jj)
                  v_ipc(ji,jj,:) = (vn(ji,jj,:)+zv_eiv(ji,jj,:))*(pbasin(ji,jj)+ibasin(ji,jj))

                  ! basin separated T times V on T points
                  vt_ind(ji,jj,:) = tn(ji,jj,:) *                                 &
                     &              ( (vn    (ji,jj,:) + vn    (ji,jj-1,:))*0.5   &
                     &              + (zv_eiv(ji,jj,:) + zv_eiv(ji,jj-1,:))*0.5 ) 
                  vt_atl(ji,jj,:) = vt_ind(ji,jj,:) * abasin(ji,jj)
                  vt_pac(ji,jj,:) = vt_ind(ji,jj,:) * pbasin(ji,jj)
                  vt_ind(ji,jj,:) = vt_ind(ji,jj,:) * ibasin(ji,jj)

                  ! basin separated S times V on T points
                  vs_ind(ji,jj,:) = sn(ji,jj,:) *                                 &
                     &              ( (vn    (ji,jj,:) + vn    (ji,jj-1,:))*0.5   &
                     &              + (zv_eiv(ji,jj,:) + zv_eiv(ji,jj-1,:))*0.5 ) 
                  vs_atl(ji,jj,:) = vs_ind(ji,jj,:) * abasin(ji,jj)
                  vs_pac(ji,jj,:) = vs_ind(ji,jj,:) * pbasin(ji,jj)
                  vs_ind(ji,jj,:) = vs_ind(ji,jj,:) * ibasin(ji,jj)
               END DO
            END DO

         ENDIF

         ! horizontal integral and vertical dz 
         v_msf_glo(:,:) = ptr_vjk( vn(:,:,:) ) 
#if defined key_diaeiv
         v_msf_eiv(:,:) = ptr_vjk( v_eiv(:,:,:) ) 
#endif
         IF( ln_subbas ) THEN
            v_msf_atl(:,:) = ptr_vjk( v_atl(:,:,:) ) 
            v_msf_ipc(:,:) = ptr_vjk( v_ipc(:,:,:) ) 
            ht_atl(:) = SUM(ptr_vjk( vt_atl(:,:,:)),2 )
            ht_pac(:) = SUM(ptr_vjk( vt_pac(:,:,:)),2 )
            ht_ind(:) = SUM(ptr_vjk( vt_ind(:,:,:)),2 )
            st_atl(:) = SUM(ptr_vjk( vs_atl(:,:,:)),2 )
            st_pac(:) = SUM(ptr_vjk( vs_pac(:,:,:)),2 )
            st_ind(:) = SUM(ptr_vjk( vs_ind(:,:,:)),2 )
         ENDIF

         ! poleward tracer transports: 
         ! overturning components:
         pht_ove(:) = SUM( v_msf_glo(:,:) * tn_jk(:,:), 2 )   ! SUM over jk
         pst_ove(:) = SUM( v_msf_glo(:,:) * sn_jk(:,:), 2 )   ! SUM over jk
#if defined key_diaeiv
         pht_eiv(:) = SUM( v_msf_eiv(:,:) * tn_jk(:,:), 2 )   ! SUM over jk
         pst_eiv(:) = SUM( v_msf_eiv(:,:) * sn_jk(:,:), 2 )   ! SUM over jk
#endif
      
         ! conversion in PW and G g
         zpwatt = zpwatt * rau0 * rcp
         pht_adv(:) = pht_adv(:) * zpwatt  
         pht_ove(:) = pht_ove(:) * zpwatt
         pht_ldf(:) = pht_ldf(:) * zpwatt
         pst_adv(:) = pst_adv(:) * zggram
         pst_ove(:) = pst_ove(:) * zggram
         pst_ldf(:) = pst_ldf(:) * zggram
#if defined key_diaeiv
         pht_eiv(:) = pht_eiv(:) * zpwatt
         pst_eiv(:) = pst_eiv(:) * zggram
#endif
         IF( ln_subbas ) THEN
            ht_atl(:) = ht_atl(:) * zpwatt
            ht_pac(:) = ht_pac(:) * zpwatt
            ht_ind(:) = ht_ind(:) * zpwatt
            st_atl(:) = st_atl(:) * zggram 
            st_pac(:) = st_pac(:) * zggram
            st_ind(:) = st_ind(:) * zggram
         ENDIF

         ! "Meridional" Stream-Function
         DO jk = 2,jpk 
            v_msf_glo(:,jk) = v_msf_glo(:,jk-1) + v_msf_glo(:,jk)
         END DO
         v_msf_glo(:,:) = v_msf_glo(:,:) * zsverdrup

#if defined key_diaeiv
         ! Bolus "Meridional" Stream-Function
         DO jk = 2,jpk 
            v_msf_eiv(:,jk) = v_msf_eiv(:,jk-1) + v_msf_eiv(:,jk)
         END DO
         v_msf_eiv(:,:) = v_msf_eiv(:,:) * zsverdrup
#endif

         IF( ln_subbas ) THEN
            DO jk = 2,jpk 
               v_msf_atl(:,jk) = v_msf_atl(:,jk-1) + v_msf_atl(:,jk)
               v_msf_ipc(:,jk) = v_msf_ipc(:,jk-1) + v_msf_ipc(:,jk)
            END DO
            v_msf_atl(:,:) = v_msf_atl(:,:) * zsverdrup
            v_msf_ipc(:,:) = v_msf_ipc(:,:) * zsverdrup
         ENDIF

         ! outputs
         CALL dia_ptr_wri( kt )

      ENDIF

      ! Close the file
      IF( kt == nitend ) CALL histclo( numptr )

   END SUBROUTINE dia_ptr


   SUBROUTINE dia_ptr_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dia_ptr_init  ***
      !!                   
      !! ** Purpose :   Initialization, namelist read
      !!
      !! ** Method  :   
      !!
      !! ** input   :   Namlist namptr
      !!
      !! ** Action  :  
      !!
      !! history :
      !!   9.0  !  03-08  (Autor Names)  Original code
      !!----------------------------------------------------------------------
      !! * local declarations
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   z_1         ! temporary workspace

      NAMELIST/namptr/ ln_diaptr, ln_subbas, nf_ptr
      !!----------------------------------------------------------------------

      ! Read Namelist namptr : poleward transport parameters
      REWIND ( numnam )
      READ   ( numnam, namptr )


      ! Control print
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dia_ptr_init : poleward transport and msf initialization'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '          Namelist namptr : set ptr parameters'
         WRITE(numout,*) '             Switch for ptr diagnostic (T) or not (F) ln_diaptr = ', ln_diaptr
         WRITE(numout,*) '             Atla/Paci/Ind basins computation         ln_subbas = ', ln_subbas
         WRITE(numout,*) '             Frequency of computation                    nf_ptr = ', nf_ptr
      ENDIF

      ! inverse of the ocean "zonal" v-point section
      z_1(:,:,:) = 1.e0
      surf_jk_r(:,:) = ptr_vtjk( z_1(:,:,:) )
      WHERE( surf_jk_r(:,:) /= 0.e0 )   surf_jk_r(:,:) = 1.e0 / surf_jk_r(:,:)

   END SUBROUTINE dia_ptr_init

   !!---------------------------------------------------------------------
   !!   Default option :                                       NetCDF file
   !!---------------------------------------------------------------------

   SUBROUTINE dia_ptr_wri( kt )
      !!---------------------------------------------------------------------
      !!                ***  ROUTINE dia_ptr_wri  ***
      !!
      !! ** Purpose :   output of poleward fluxes
      !!
      !! ** Method  :   NetCDF file
      !!
      !! History :
      !!   9.0  !  03-09  (G. Madec)  Original code
      !!----------------------------------------------------------------------
      USE ioipsl          ! NetCDF IPSL library
      USE daymod

      !! * Arguments
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index

      !! * Save variables   
      INTEGER, SAVE ::   nhoridz, ndepidzt, ndepidzw, ndex(1)

      !! * Local variables
      CHARACTER (len=40) ::   &
         clhstnam, clop             ! temporary names
      INTEGER ::   iline, it, ji    !
      REAL(wp) ::   &
         zsto, zout, zdt, zmax, &   ! temporary scalars
         zjulian
      REAL(wp), DIMENSION(jpj) ::   zphi, zfoo
      !!----------------------------------------------------------------------
      
      ! Define frequency of output and means
      zdt = rdt
      IF( nacc == 1 ) zdt = rdtmin
#if defined key_diainstant
         zsto = nf_ptr * zdt
         clop = "inst(x)"               ! no use of the mask value (require less cpu time)
         !!! clop="inst(only(x))"       ! put 1.e+20 on land (very expensive!!)
#else
         zsto = zdt
         clop = "ave(x)"                ! no use of the mask value (require less cpu time)
         !!! clop="ave(only(x))"        ! put 1.e+20 on land (very expensive!!)
#endif
      zout = nf_ptr * zdt
      zmax = ( nitend - nit000 + 1 ) * zdt
      
         
      ! define time axis
      it = kt - nit000 + 1

      ! Initialization
      ! --------------
      IF( kt == nit000 ) THEN
      
      zdt = rdt
      IF( nacc == 1 ) zdt = rdtmin

         ! Reference latitude
         ! ------------------
         !                                           ! =======================
!!DB: Orca-related
!         IF( cp_cfg == "orca" ) THEN                 !   ORCA configurations
!            !                                        ! =======================!
!
!            IF( jp_cfg == 05  )   iline = 192   ! i-line that passes near the North Pole
!            IF( jp_cfg == 025 )   iline = 384   ! i-line that passes near the North Pole
!            IF( jp_cfg == 2   )   iline =  48   ! i-line that passes near the North Pole
!            IF( jp_cfg == 4   )   iline =  24   ! i-line that passes near the North Pole
!            zphi(:) = 0.e0
!            DO ji = mi0(iline), mi1(iline) 
!               zphi(:) = gphiv(ji,:)         ! if iline is in the local domain
!               ! correct highest latitude for ORCA05
!               IF( jp_cfg == 05  ) zphi(jpj) = zphi(jpjm1) + (zphi(jpjm1)-zphi(jpj-2))/2.
!               IF( jp_cfg == 05  ) zphi(jpj) = MIN( zphi(jpj), 90.)!!!
!
!            END DO
!            ! provide the correct zphi to all local domains
!            IF( lk_mpp )   CALL mpp_sum( zphi, jpj )        !
!
!            !                                        ! =======================
!         ELSE                                        !   OTHER configurations!
            !                                        ! =======================!
            zphi(:) = gphiv(1,:)             ! assume lat/lon coordinate, select the first i-line
            !
!         ENDIF

         ! OPEN netcdf file 
         ! ----------------
         ! Define frequency of output and means
         zsto = nf_ptr * zdt
         clop = "ave(x)"
         zout = nf_ptr * zdt
         zfoo(:) = 0.e0

         ! Compute julian date from starting date of the run

         CALL ymds2ju( nyear, nmonth, nday, 0.e0, zjulian )

         CALL dia_nam( clhstnam, nf_ptr, 'diaptr' )
         IF(lwp)WRITE( numout,*)" Name of diaptr NETCDF file ",clhstnam

         ! Horizontal grid : zphi()
         CALL histbeg(clhstnam, 1, zfoo, jpj, zphi,   &
            1, 1, 1, jpj, 0, zjulian, zdt, nhoridz, numptr, domain_id=nidom )
         ! Vertical grids : gdept, gdepw
         CALL histvert( numptr, "deptht", "Vertical T levels",   &
            "m", jpk, gdept, ndepidzt )
         CALL histvert( numptr, "depthw", "Vertical W levels",   &
            "m", jpk, gdepw, ndepidzw )
         
         !  Zonal mean T and S
         
         CALL histdef( numptr, "zotemglo", "Zonal Mean Temperature","C" ,   &
            1, jpj, nhoridz, jpk, 1, jpk, ndepidzt, 32, clop, zsto, zout )
         CALL histdef( numptr, "zosalglo", "Zonal Mean Salinity","PSU"  ,   &
            1, jpj, nhoridz, jpk, 1, jpk, ndepidzt, 32, clop, zsto, zout )

         !  Meridional Stream-Function (eulerian and bolus)
         
         CALL histdef( numptr, "zomsfglo", "Meridional Stream-Function: Global","Sv" ,   &
            1, jpj, nhoridz, jpk, 1, jpk, ndepidzw, 32, clop, zsto, zout )
         IF( ln_subbas ) THEN
            CALL histdef( numptr, "zomsfatl", "Meridional Stream-Function: Atlantic","Sv" ,   &
               1, jpj, nhoridz, jpk, 1, jpk, ndepidzw, 32, clop, zsto, zout )
            CALL histdef( numptr, "zomsfipc", "Meridional Stream-Function: Indo-Pacific","Sv" ,&
               1, jpj, nhoridz, jpk, 1, jpk, ndepidzw, 32, clop, zsto, zout )
         ENDIF

         !  Heat transport 

         CALL histdef( numptr, "sophtadv", "Advective Heat Transport"      ,   &
            "PW", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
         CALL histdef( numptr, "sophtldf", "Diffusive Heat Transport"      ,   &
            "PW",1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
         CALL histdef( numptr, "sophtove", "Overturning Heat Transport"    ,   &
            "PW",1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
         IF( ln_subbas ) THEN
            CALL histdef( numptr, "sohtatl", "Heat Transport Atlantic"      ,  &
               "PW", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
            CALL histdef( numptr, "sohtpac", "Heat Transport Pacific"      ,   &
               "PW", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
            CALL histdef( numptr, "sohtind", "Heat Transport Indic"      ,     &
               "PW", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
         ENDIF


         !  Salt transport 

         CALL histdef( numptr, "sopstadv", "Advective Salt Transport"      ,   &
            "Giga g/s", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
         CALL histdef( numptr, "sopstldf", "Diffusive Salt Transport"      ,   &
            "Giga g/s", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
         CALL histdef( numptr, "sopstove", "Overturning Salt Transport"    ,   &
            "Giga g/s", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )

#if defined key_diaeiv
         ! Eddy induced velocity
         CALL histdef( numptr, "zomsfeiv", "Bolus Meridional Stream-Function: global",   &
            "Sv"      , 1, jpj, nhoridz, jpk, 1, jpk, ndepidzw, 32, clop, zsto, zout )
         CALL histdef( numptr, "sophteiv", "Bolus Advective Heat Transport",   &
            "PW"      , 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
         CALL histdef( numptr, "sopsteiv", "Bolus Advective Salt Transport",   &
            "Giga g/s", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
#endif
         IF( ln_subbas ) THEN
            CALL histdef( numptr, "sostatl", "Salt Transport Atlantic"      ,    &
               "Giga g/s", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
            CALL histdef( numptr, "sostpac", "Salt Transport Pacific"      ,     &
               "Giga g/s", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
            CALL histdef( numptr, "sostind", "Salt Transport Indic"      ,       &
               "Giga g/s", 1, jpj, nhoridz, 1, 1, 1, -99, 32, clop, zsto, zout )
         ENDIF
         

         CALL histend( numptr )

      ENDIF

      IF( MOD( kt, nf_ptr ) == 0 ) THEN

         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'dia_ptr : write Poleward Transports at time-step : ', kt
            WRITE(numout,*) '~~~~~~~~'
            WRITE(numout,*)
         ENDIF

         ! define time axis
         it= kt - nit000 + 1
         ndex(1) = 0
         CALL histwrite( numptr, "zotemglo", it, tn_jk    , jpj*jpk, ndex )
         CALL histwrite( numptr, "zosalglo", it, sn_jk    , jpj*jpk, ndex )
         ! overturning outputs:
         CALL histwrite( numptr, "zomsfglo", it, v_msf_glo , jpj*jpk, ndex )
         IF( ln_subbas ) THEN
            CALL histwrite( numptr, "zomsfatl", it, v_msf_atl , jpj*jpk, ndex )
            CALL histwrite( numptr, "zomsfipc", it, v_msf_ipc , jpj*jpk, ndex )
         ENDIF
         ! heat transport outputs:
         IF( ln_subbas ) THEN
            CALL histwrite( numptr, "sohtatl", it, ht_atl  , jpj, ndex )
            CALL histwrite( numptr, "sohtpac", it, ht_pac  , jpj, ndex )
            CALL histwrite( numptr, "sohtind", it, ht_ind  , jpj, ndex )
            CALL histwrite( numptr, "sostatl", it, st_atl  , jpj, ndex )
            CALL histwrite( numptr, "sostpac", it, st_pac  , jpj, ndex )
            CALL histwrite( numptr, "sostind", it, st_ind  , jpj, ndex )
         ENDIF

         CALL histwrite( numptr, "sophtadv", it, pht_adv  , jpj, ndex )
         CALL histwrite( numptr, "sophtldf", it, pht_ldf  , jpj, ndex )
         CALL histwrite( numptr, "sophtove", it, pht_ove  , jpj, ndex )
         CALL histwrite( numptr, "sopstadv", it, pst_adv  , jpj, ndex )
         CALL histwrite( numptr, "sopstldf", it, pst_ldf  , jpj, ndex )
         CALL histwrite( numptr, "sopstove", it, pst_ove  , jpj, ndex )
#if defined key_diaeiv
         CALL histwrite( numptr, "zomsfeiv", it, v_msf_eiv, jpj*jpk, ndex )
         CALL histwrite( numptr, "sophteiv", it, pht_eiv  , jpj    , ndex )
         CALL histwrite( numptr, "sopsteiv", it, pst_eiv  , jpj    , ndex )
#endif
 
      ENDIF

   END SUBROUTINE dia_ptr_wri

   !!======================================================================
END MODULE diaptr
