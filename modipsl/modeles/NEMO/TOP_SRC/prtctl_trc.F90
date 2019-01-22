MODULE prtctl_trc
   !!==============================================================================
   !!                       ***  MODULE prtctl   ***
   !! Ocean system   : print all SUM trends for each processor domain
   !!==============================================================================
#if defined key_passivetrc

   USE par_trc_trp
   USE oce_trc          ! ocean space and time domain variables
   USE in_out_manager   ! I/O manager
   USE lib_mpp          ! distributed memory computing

   IMPLICIT NONE
   PRIVATE

   !! * Module declaration
   INTEGER, DIMENSION(:), ALLOCATABLE, SAVE ::   &  !:
      nlditl , nldjtl ,   &  !: first, last indoor index for each i-domain
      nleitl , nlejtl ,   &  !: first, last indoor index for each j-domain
      nimpptl, njmpptl,   &  !: i-, j-indexes for each processor
      nlcitl , nlcjtl ,   &  !: dimensions of every subdomain
      ibonitl, ibonjtl

   REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   &  !:
      tra_ctl                   !: previous trend values

   !! * Routine accessibility
   PUBLIC prt_ctl_trc         ! called by all subroutines
   PUBLIC prt_ctl_trc_info    !
   PUBLIC prt_ctl_trc_init    ! called by opa.F90
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2005) 
   !! $Header: /home/opalod/NEMOCVSROOT/NEMO/TOP_SRC/prtctl_trc.F90,v 1.2 2006/03/21 15:53:52 opalod Exp $ 
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt 
   !!----------------------------------------------------------------------


CONTAINS

   SUBROUTINE prt_ctl_trc (tab4d, mask, clinfo, ovlap, kdim, clinfo2)
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE prt_ctl  ***
      !!
      !! ** Purpose : - print sum control 3D arrays over the same area 
      !!                in mono and mpp case. This way can be usefull when
      !!                debugging a new parametrization in mono or mpp. 
      !!
      !! ** Method  : 2 possibilities exist when setting the ln_ctl parameter to
      !!                .true. in the ocean namelist:
      !!              - to debug a MPI run .vs. a mono-processor one; 
      !!                the control print will be done over each sub-domain.
      !!                The nictl[se] and njctl[se] parameters in the namelist must 
      !!                be set to zero and [ij]splt to the corresponding splitted
      !!                domain in MPI along respectively i-, j- directions.
      !!              - to debug a mono-processor run over the whole domain/a specific area; 
      !!                in the first case the nictl[se] and njctl[se] parameters must be set
      !!                to zero else to the indices of the area to be controled. In both cases
      !!                isplt and jsplt must be set to 1.
      !!              - All arguments of the above calling sequence are optional so their
      !!                name must be explicitly typed if used. For instance if the mask
      !!                array tmask(:,:,:) must be passed through the prt_ctl subroutine, 
      !!                it must looks like: CALL prt_ctl(mask=tmask).
      !!
      !!                    tab4d   : 4D array
      !!                    mask    : mask (3D) to apply to the tab4d array
      !!                    clinfo  : information about the tab3d array
      !!                    ovlap   : overlap value
      !!                    kdim    : k- direction for 4D arrays 
      !!
      !! History :
      !!   9.0  !  05-07  (C. Talandier) original code
      !!        !  05-10  (C. Ethe     ) adapted to passive tracer
      !!----------------------------------------------------------------------
      !! * Arguments
      REAL(wp), DIMENSION(:,:,:,:), INTENT(in), OPTIONAL :: tab4d
      REAL(wp), DIMENSION(:,:,:), INTENT(in), OPTIONAL :: mask
      CHARACTER (len=*), DIMENSION(:), INTENT(in), OPTIONAL :: clinfo
      CHARACTER (len=*), INTENT(in), OPTIONAL :: clinfo2
      INTEGER, INTENT(in), OPTIONAL :: ovlap
      INTEGER, INTENT(in), OPTIONAL :: kdim

      !! * Local declarations
      INTEGER  :: overlap, numid, jn, js, sind, eind, kdir
      REAL(wp) :: zsum, zvctl
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zmask, ztab3d
      CHARACTER (len=20), DIMENSION(jptra) :: cl
      CHARACTER (len=10) :: cl2
      !!----------------------------------------------------------------------

      ! Arrays, scalars initialization 
      overlap       = 0
      kdir          = jpkm1
      zsum          = 0.e0
      zvctl         = 0.e0
      cl(:)         = ''
      cl2           = ''
      ztab3d(:,:,:) = 0.e0
      zmask (:,:,:) = 1.e0

      ! Control of optional arguments

      IF( PRESENT(ovlap)   )  overlap       = ovlap
      IF( PRESENT(kdim)    )  kdir          = kdim
      IF( PRESENT(clinfo ) )  cl(:)         = clinfo(:)
      IF( PRESENT(clinfo2) )  cl2           = clinfo2
      IF( PRESENT(mask)    )  zmask (:,:,:) = mask(:,:,:)

      IF( lk_mpp )   THEN
         ! processor number
         sind = narea
         eind = narea
      ELSE
         ! processors total number
         sind = 1
         eind = ijsplt
      ENDIF

      ! Loop over each sub-domain, i.e. the total number of processors ijsplt
      DO js = sind, eind

         numid = 90 + js

         ! Set indices for the SUM control
         IF( .NOT. lsp_area ) THEN
            IF (lk_mpp )   THEN
               nictls = MAX( 1, nlditl(js) - overlap )
               nictle = nleitl(js) + overlap * MIN( 1, nlcitl(js) - nleitl(js)) 
               njctls = MAX( 1, nldjtl(js) - overlap )
               njctle = nlejtl(js) + overlap * MIN( 1, nlcjtl(js) - nlejtl(js))
               ! Do not take into account the bound of the domain
               IF( ibonitl(js) == -1 .OR. ibonitl(js) == 2 ) nictls = MAX(2, nictls)
               IF( ibonitl(js) ==  1 .OR. ibonitl(js) == 2 ) nictle = MIN(nictle, nleitl(js) - 1)
               IF( ibonjtl(js) == -1 .OR. ibonjtl(js) == 2 ) njctls = MAX(2, njctls)
               IF( ibonjtl(js) ==  1 .OR. ibonjtl(js) == 2 ) njctle = MIN(njctle, nlejtl(js) - 1)
            ELSE
               nictls = MAX( 1, nimpptl(js) + nlditl(js) - 1 - overlap )
               nictle = nimpptl(js) + nleitl(js) - 1 + overlap * MIN( 1, nlcitl(js) - nleitl(js) ) 
               njctls = MAX( 1, njmpptl(js) + nldjtl(js) - 1 - overlap )
               njctle = njmpptl(js) + nlejtl(js) - 1 + overlap * MIN( 1, nlcjtl(js) - nlejtl(js) ) 
               ! Do not take into account the bound of the domain
               IF( ibonitl(js) == -1 .OR. ibonitl(js) == 2 ) nictls = MAX(2, nictls)
               IF( ibonjtl(js) == -1 .OR. ibonjtl(js) == 2 ) njctls = MAX(2, njctls)
               IF( ibonitl(js) ==  1 .OR. ibonitl(js) == 2 ) nictle = MIN(nictle, nimpptl(js) + nleitl(js) - 2)
               IF( ibonjtl(js) ==  1 .OR. ibonjtl(js) == 2 ) njctle = MIN(njctle, njmpptl(js) + nlejtl(js) - 2)
            ENDIF
         ENDIF
         
         IF( PRESENT(clinfo2) ) THEN
            DO jn = 1, jptra
               zvctl  = tra_ctl(jn,js)
               ztab3d(:,:,:) = tab4d(:,:,:,jn)
               zsum          = SUM( ztab3d(nictls:nictle,njctls:njctle,1:kdir) &
                  &                 *zmask(nictls:nictle,njctls:njctle,1:kdir) )
               WRITE(numid,FMT="(3x,a,' : ',D23.16)") cl(jn), zsum-zvctl
               tra_ctl(jn,js) = zsum
            ENDDO
         ELSE
            DO jn = 1, jptra
               ztab3d(:,:,:) = tab4d(:,:,:,jn)
               zsum          = SUM( ztab3d(nictls:nictle,njctls:njctle,1:kdir) &
                  &               * zmask(nictls:nictle,njctls:njctle,1:kdir) )
               WRITE(numid,FMT="(3x,a,' : ',D23.16)") cl(jn), zsum
            END DO
         ENDIF
         

      ENDDO

   END SUBROUTINE prt_ctl_trc

   SUBROUTINE prt_ctl_trc_info (clinfo)
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE prt_ctl_trc_info  ***
      !!
      !! ** Purpose : - print information without any computation
      !!
      !! ** Action  : - input arguments
      !!                    clinfo : information to print
      !!
      !! History :
      !!   9.0  !  05-07  (C. Talandier) original code
      !!----------------------------------------------------------------------
      !! * Arguments
      CHARACTER (len=*), INTENT(in) ::   clinfo

      !! * Local declarations
      INTEGER ::  numid, js, sind, eind
      !!----------------------------------------------------------------------

      IF( lk_mpp )   THEN
         ! processor number
         sind = narea
         eind = narea
      ELSE
         ! total number of processors
         sind = 1
         eind = ijsplt
      ENDIF

      ! Loop over each sub-domain, i.e. number of processors ijsplt
      DO js = sind, eind
         numid = 90 + js
         WRITE(numid,*)clinfo
      ENDDO


   END SUBROUTINE prt_ctl_trc_info

   SUBROUTINE prt_ctl_trc_init
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE prt_ctl_trc_init  ***
      !!
      !! ** Purpose :   open ASCII files & compute indices
      !!
      !! History :
      !!   9.0  !  05-07  (C. Talandier) original code
      !!        !  05-10  (C. Ethe     ) adapted to passive tracer
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER ::   js, numid, sind, eind
      CHARACTER (len=31) :: clfile_out
      CHARACTER (len=27) :: clb_name
      CHARACTER (len=19) :: cl_run
      !!----------------------------------------------------------------------

      ! Allocate arrays
      ALLOCATE(nlditl (ijsplt))
      ALLOCATE(nldjtl (ijsplt))
      ALLOCATE(nleitl (ijsplt))
      ALLOCATE(nlejtl (ijsplt))
      ALLOCATE(nimpptl(ijsplt))
      ALLOCATE(njmpptl(ijsplt))
      ALLOCATE(nlcitl (ijsplt))
      ALLOCATE(nlcjtl (ijsplt))
      ALLOCATE(tra_ctl(jptra,ijsplt))
      ALLOCATE(ibonitl(ijsplt))
      ALLOCATE(ibonjtl(ijsplt))

      ! Initialization 
      tra_ctl (:,:)=0.e0

      IF( lk_mpp ) THEN
         sind = narea
         eind = narea
         clb_name = "('mpp.top.output_',I3.3)"
         cl_run = 'MULTI processor run'
         ! use indices for each area computed by mpp_init subroutine
         nlditl(:) = nldit(:) 
         nleitl(:) = nleit(:) 
         nldjtl(:) = nldjt(:) 
         nlejtl(:) = nlejt(:) 
         !
         nimpptl(:) = nimppt(:)
         njmpptl(:) = njmppt(:)
         !
         nlcitl(:) = nlcit(:)
         nlcjtl(:) = nlcjt(:)
         !
         ibonitl(:) = ibonit(:)
         ibonjtl(:) = ibonjt(:)
      ELSE
         sind = 1
         eind = ijsplt
         clb_name = "('mono.top.output_',I3.3)"
         cl_run = 'MONO processor run '
         ! compute indices for each area as done in mpp_init subroutine
         CALL sub_dom
      ENDIF

      DO js = sind, eind
         numid = 90 + js
         WRITE(clfile_out,FMT=clb_name) js-1
         OPEN ( UNIT=numid, FILE=TRIM(clfile_out),FORM='FORMATTED' )
         WRITE(numid,*)
         WRITE(numid,*) '                 L O D Y C - I P S L'
         WRITE(numid,*) '                     O P A model'
         WRITE(numid,*) '            Ocean General Circulation Model'
         WRITE(numid,*) '               version OPA 9.0  (2005) '
         WRITE(numid,*)
         WRITE(numid,*) '                   PROC number: ', js
         WRITE(numid,*)
         WRITE(numid,FMT="(19x,a20)")cl_run

         ! Print the SUM control indices
         IF( .NOT. lsp_area )   THEN
            IF ( lk_mpp )   THEN
               nictls = nlditl(js) 
               nictle = nleitl(js)
               njctls = nldjtl(js)
               njctle = nlejtl(js)
            ELSE
               nictls = nimpptl(js) + nlditl(js) - 1
               nictle = nimpptl(js) + nleitl(js) - 1
               njctls = njmpptl(js) + nldjtl(js) - 1
               njctle = njmpptl(js) + nlejtl(js) - 1
            ENDIF
         ENDIF
         WRITE(numid,*) 
         WRITE(numid,*) 'prt_tra_ctl :  Sum control indices'
         WRITE(numid,*) '~~~~~~~'
         WRITE(numid,*)
         WRITE(numid,9000)'                                nlej   = ', nlejtl(js), '              '
         WRITE(numid,9000)'                  ------------- njctle = ', njctle, ' -------------'
         WRITE(numid,9001)'                  |                                       |'
         WRITE(numid,9001)'                  |                                       |'
         WRITE(numid,9001)'                  |                                       |'
         WRITE(numid,9002)'           nictls = ', nictls,  '                           nictle = ', nictle
         WRITE(numid,9002)'           nldi   = ', nlditl(js),  '                           nlei   = ', nleitl(js)
         WRITE(numid,9001)'                  |                                       |'
         WRITE(numid,9001)'                  |                                       |'
         WRITE(numid,9001)'                  |                                       |'
         WRITE(numid,9004)'  njmpp  = ',njmpptl(js),'   ------------- njctls = ', njctls, ' -------------'
         WRITE(numid,9003)'           nimpp  = ', nimpptl(js), '        nldj   = ', nldjtl(js), '              '
         WRITE(numid,*)
         WRITE(numid,*)

9000     FORMAT(a41,i4.4,a14)
9001     FORMAT(a59)
9002     FORMAT(a20,i4.4,a36,i3.3)
9003     FORMAT(a20,i4.4,a17,i4.4)
9004     FORMAT(a11,i4.4,a26,i4.4,a14)
      ENDDO

   END SUBROUTINE prt_ctl_trc_init


   SUBROUTINE sub_dom
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sub_dom  ***
      !!                    
      !! ** Purpose :   Lay out the global domain over processors. 
      !!                CAUTION: 
      !!                This part has been extracted from the mpp_init
      !!                subroutine and names of variables/arrays have been 
      !!                slightly changed to avoid confusion but the computation
      !!                is exactly the same. Any modification about indices of
      !!                each sub-domain in the mppini.F90 module should be reported 
      !!                here.
      !!
      !! ** Method  :   Global domain is distributed in smaller local domains.
      !!                Periodic condition is a function of the local domain position
      !!                (global boundary or neighbouring domain) and of the global
      !!                periodic
      !!                Type :         jperio global periodic condition
      !!                               nperio local  periodic condition
      !!
      !! ** Action  : - set domain parameters
      !!                    nimpp     : longitudinal index 
      !!                    njmpp     : latitudinal  index
      !!                    nperio    : lateral condition type 
      !!                    narea     : number for local area
      !!                    nlcil      : first dimension
      !!                    nlcjl      : second dimension
      !!                    nbondil    : mark for "east-west local boundary"
      !!                    nbondjl    : mark for "north-south local boundary"
      !!
      !! History :
      !!        !  94-11  (M. Guyon)  Original code
      !!        !  95-04  (J. Escobar, M. Imbard)
      !!        !  98-02  (M. Guyon)  FETI method
      !!        !  98-05  (M. Imbard, J. Escobar, L. Colombet )  SHMEM and MPI versions
      !!   8.5  !  02-08  (G. Madec)  F90 : free form
      !!----------------------------------------------------------------------
      !! * Local variables
      INTEGER ::   ji, jj, js               ! dummy loop indices
      INTEGER ::   &
         ii, ij,                         &  ! temporary integers
         irestil, irestjl,               &  !    "          "
         ijpi  , ijpj, nlcil,            &  ! temporary logical unit
         nlcjl , nbondil, nbondjl,       &
         nrecil, nrecjl, nldil, nleil, nldjl, nlejl

      INTEGER, DIMENSION(:,:), ALLOCATABLE ::   &
         iimpptl, ijmpptl, ilcitl, ilcjtl       ! temporary workspace
      REAL(wp) ::   zidom, zjdom            ! temporary scalars
      !!----------------------------------------------------------------------

      !  1. Dimension arrays for subdomains
      ! -----------------------------------
      !  Computation of local domain sizes ilcitl() ilcjtl()
      !  These dimensions depend on global sizes isplt,jsplt and jpiglo,jpjglo
      !  The subdomains are squares leeser than or equal to the global
      !  dimensions divided by the number of processors minus the overlap
      !  array (cf. par_oce.F90).

      ijpi = ( jpiglo-2*jpreci + (isplt-1) ) / isplt + 2*jpreci
      ijpj = ( jpjglo-2*jprecj + (jsplt-1) ) / jsplt + 2*jprecj

      ALLOCATE(ilcitl (isplt,jsplt))
      ALLOCATE(ilcjtl (isplt,jsplt))

      nrecil  = 2 * jpreci
      nrecjl  = 2 * jprecj
      irestil = MOD( jpiglo - nrecil , isplt )
      irestjl = MOD( jpjglo - nrecjl , jsplt )

      IF(  irestil == 0 )   irestil = isplt
      DO jj = 1, jsplt
         DO ji = 1, irestil
            ilcitl(ji,jj) = ijpi
         END DO
         DO ji = irestil+1, isplt
            ilcitl(ji,jj) = ijpi -1
         END DO
      END DO
      
      IF( irestjl == 0 )   irestjl = jsplt
      DO ji = 1, isplt
         DO jj = 1, irestjl
            ilcjtl(ji,jj) = ijpj
         END DO
         DO jj = irestjl+1, jsplt
            ilcjtl(ji,jj) = ijpj -1
         END DO
      END DO
      
      zidom = nrecil
      DO ji = 1, isplt
         zidom = zidom + ilcitl(ji,1) - nrecil
      END DO
      
      zjdom = nrecjl
      DO jj = 1, jsplt
         zjdom = zjdom + ilcjtl(1,jj) - nrecjl
      END DO

      !  2. Index arrays for subdomains
      ! -------------------------------

      ALLOCATE(iimpptl(isplt,jsplt))
      ALLOCATE(ijmpptl(isplt,jsplt))
      
      iimpptl(:,:) = 1
      ijmpptl(:,:) = 1
      
      IF( isplt > 1 ) THEN
         DO jj = 1, jsplt
            DO ji = 2, isplt
               iimpptl(ji,jj) = iimpptl(ji-1,jj) + ilcitl(ji-1,jj) - nrecil
            END DO
         END DO
      ENDIF

      IF( jsplt > 1 ) THEN
         DO jj = 2, jsplt
            DO ji = 1, isplt
               ijmpptl(ji,jj) = ijmpptl(ji,jj-1)+ilcjtl(ji,jj-1)-nrecjl
            END DO
         END DO
      ENDIF
      
      ! 3. Subdomain description
      ! ------------------------

      DO js = 1, ijsplt
         ii = 1 + MOD( js-1, isplt )
         ij = 1 + (js-1) / isplt
         nimpptl(js) = iimpptl(ii,ij)
         njmpptl(js) = ijmpptl(ii,ij)
         nlcitl (js) = ilcitl (ii,ij)     
         nlcil       = nlcitl (js)     
         nlcjtl (js) = ilcjtl (ii,ij)     
         nlcjl       = nlcjtl (js)
         nbondjl = -1                                    ! general case
         IF( js   >  isplt          )   nbondjl = 0      ! first row of processor
         IF( js   >  (jsplt-1)*isplt )  nbondjl = 1     ! last  row of processor
         IF( jsplt == 1             )   nbondjl = 2      ! one processor only in j-direction
         ibonjtl(js) = nbondjl
         
         nbondil = 0                                     ! 
         IF( MOD( js, isplt ) == 1 )   nbondil = -1      !
         IF( MOD( js, isplt ) == 0 )   nbondil =  1      !
         IF( isplt            == 1 )   nbondil =  2      ! one processor only in i-direction
         ibonitl(js) = nbondil
         
         nldil =  1   + jpreci
         nleil = nlcil - jpreci
         IF( nbondil == -1 .OR. nbondil == 2 )   nldil = 1
         IF( nbondil ==  1 .OR. nbondil == 2 )   nleil = nlcil
         nldjl =  1   + jprecj
         nlejl = nlcjl - jprecj
         IF( nbondjl == -1 .OR. nbondjl == 2 )   nldjl = 1
         IF( nbondjl ==  1 .OR. nbondjl == 2 )   nlejl = nlcjl
         nlditl(js) = nldil
         nleitl(js) = nleil
         nldjtl(js) = nldjl
         nlejtl(js) = nlejl
      END DO

      DEALLOCATE(iimpptl)
      DEALLOCATE(ijmpptl)
      DEALLOCATE(ilcitl)
      DEALLOCATE(ilcjtl)

   END SUBROUTINE sub_dom
 
#else
   !!----------------------------------------------------------------------
   !!   Dummy module :                      NO passive tracer
   !!----------------------------------------------------------------------
#endif
    
   !!======================================================================

END MODULE prtctl_trc
