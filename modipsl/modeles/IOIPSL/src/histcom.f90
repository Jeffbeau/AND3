MODULE histcom
!-
!$Header: /home/ioipsl/CVSROOT/IOIPSL/src/histcom.f90,v 2.3 2005/10/10 08:02:57 adm Exp $
!-
  USE netcdf
!-
  USE stringop,  ONLY : nocomma,cmpblank,findpos,find_str,strlowercase
  USE mathelp,   ONLY : mathop,moycum,trans_buff,buildop
  USE fliocom,   ONLY : flio_dom_file,flio_dom_att
  USE calendar
  USE errioipsl, ONLY : ipslerr
!-
  IMPLICIT NONE
!-
  PRIVATE
  PUBLIC :: histbeg, histdef, histhori, histvert, histend, &
 &          histwrite, histclo, histsync, ioconf_modname
!---------------------------------------------------------------------
!- Some confusing vocabulary in this code !
!- =========================================
!-
!- A REGULAR grid is a grid which is i,j indexes
!- and thus it is stored in a 2D matrix.
!- This is opposed to a IRREGULAR grid which is only in a vector
!- and where we do not know which neighbors we have.
!- As a consequence we need the bounds for each grid-cell.
!-
!- A RECTILINEAR grid is a special case of a regular grid
!- in which all longitudes for i constant are equal
!- and all latitudes for j constant.
!- In other words we do not need the full 2D matrix
!- to describe the grid, just two vectors.
!---------------------------------------------------------------------
  INTERFACE histwrite
!---------------------------------------------------------------------
!- The "histwrite" routines will give the data to the I/O system.
!- It will trigger the operations to be performed,
!- and the writting to the file if needed
!-
!- We test for the work to be done at this time here so that at a
!- later stage we can call different operation and write subroutine
!- for the REAL and INTEGER interfaces
!-
!- INPUT
!- pfileid  : The ID of the file on which this variable is to be,
!-            written. The variable should have been defined in
!-            this file before.
!- pvarname : The short name of the variable
!- pitau    : Current timestep
!- pdata    : The variable, I mean the real data !
!- nbindex  : The number of indexes provided. If it is equal to
!-            the size of the full field as provided in histdef
!-            then nothing is done.
!- nindex   : The indices used to expand the variable (pdata)
!-            onto the full field.
!---------------------------------------------------------------------
!- histwrite - we have to prepare different type of fields :
!-             real and integer, 1,2 or 3D
    MODULE PROCEDURE histwrite_r1d,histwrite_r2d,histwrite_r3d
  END INTERFACE
!-
  INTERFACE histbeg
!!  MODULE PROCEDURE histbeg_regular,histbeg_irregular
    MODULE PROCEDURE histbeg_totreg,histbeg_regular,histbeg_irregular
  END INTERFACE
!-
  INTERFACE histhori
    MODULE PROCEDURE histhori_regular,histhori_irregular
  END INTERFACE
!-
! Fixed parameter
!-
  INTEGER,PARAMETER :: nb_files_max=20, nb_var_max=400, &
 &                     nb_hax_max=5, nb_zax_max=10, nbopp_max=10
  REAL,PARAMETER :: missing_val=1.e20
!-                 or HUGE(1.0) maximum real number
!-
  INTEGER :: bufftmp_max(nb_files_max) = 1
!-
! Time variables
!-
  INTEGER,SAVE :: itau0(nb_files_max)=0
  REAL,DIMENSION(nb_files_max),SAVE ::date0,deltat
!-
! Counter of elements
!-
  INTEGER,SAVE :: nb_files=0
  INTEGER,DIMENSION(nb_files_max),SAVE :: nb_var=0, nb_tax=0
!-
! NETCDF IDs for files and axes
!-
  INTEGER,DIMENSION(nb_files_max),SAVE :: ncdf_ids,xid,yid,tid
  CHARACTER(LEN=500),SAVE :: assc_file=''
!-
! General definitions in the NETCDF file
!-
  INTEGER,DIMENSION(nb_files_max,2),SAVE :: &
 &   full_size=0,slab_ori,slab_sz
!-
! The horizontal axes
!-
  INTEGER,SAVE :: nb_hax(nb_files_max)=0
  CHARACTER(LEN=25),SAVE :: hax_name(nb_files_max,nb_hax_max,2)
!-
! The vertical axes
!-
  INTEGER,SAVE :: nb_zax(nb_files_max)=0
  INTEGER,DIMENSION(nb_files_max,nb_zax_max),SAVE :: &
  &  zax_size,zax_ids,zax_name_length
  CHARACTER(LEN=20),SAVE :: zax_name(nb_files_max,nb_zax_max)
!-
! Informations on each variable
!-
  INTEGER,DIMENSION(nb_files_max,nb_var_max),SAVE :: &
 &  name_length,nbopp
  CHARACTER(LEN=20),DIMENSION(nb_files_max,nb_var_max),SAVE :: &
 &  name,unit_name
  CHARACTER(LEN=80),DIMENSION(nb_files_max,nb_var_max),SAVE :: &
 &  title,fullop
  CHARACTER(LEN=7),SAVE :: topp(nb_files_max,nb_var_max)
  CHARACTER(LEN=7),SAVE :: sopps(nb_files_max,nb_var_max,nbopp_max)
  REAL,SAVE :: scal(nb_files_max,nb_var_max,nbopp_max)
!- Sizes of the associated grid and zommed area
  INTEGER,DIMENSION(nb_files_max,nb_var_max,3),SAVE :: &
 &   scsize,zorig,zsize
!- Sizes for the data as it goes through the various math operations
  INTEGER,SAVE :: datasz_in(nb_files_max,nb_var_max,3) = -1
  INTEGER,SAVE :: datasz_max(nb_files_max,nb_var_max) = -1
!-
  INTEGER,DIMENSION(nb_files_max,nb_var_max),SAVE :: &
 &  var_haxid,var_zaxid,var_axid,ncvar_ids
!-
  REAL,SAVE :: minmax(nb_files_max,nb_var_max,2)
!-
  REAL,DIMENSION(nb_files_max,nb_var_max),SAVE :: &
 &  freq_opp,freq_wrt
  INTEGER,DIMENSION(nb_files_max,nb_var_max),SAVE :: &
 &  last_opp,last_wrt,last_opp_chk,last_wrt_chk,nb_opp,nb_wrt,point
!-
! Book keeping for the buffers
!-
  INTEGER,SAVE :: buff_pos=0
  REAL,ALLOCATABLE,SAVE :: buffer(:)
  LOGICAL,SAVE :: &
 &  zoom(nb_files_max)=.FALSE., regular(nb_files_max)=.TRUE.
!-
! Book keeping of the axes
!-
  INTEGER,DIMENSION(nb_files_max,nb_var_max),SAVE :: &
 &  tdimid,tax_last,tax_name_length
  CHARACTER(LEN=40),DIMENSION(nb_files_max,nb_var_max),SAVE :: &
 &  tax_name
!-
! A list of functions which require special action
! (Needs to be updated when functions are added
!  but they are well located here)
!-
  CHARACTER(LEN=120),SAVE :: &
 &  indchfun = 'scatter, fill, gather, coll', &
 &  fuchnbout = 'scatter, fill'
!- Some configurable variables with locks
  CHARACTER(LEN=80),SAVE :: model_name='An IPSL model'
  LOGICAL,SAVE :: lock_modname=.FALSE.
!-
!===
CONTAINS
!===
!-
SUBROUTINE histbeg_totreg                     &
 & (pfilename, pim, plon, pjm, plat,          &
 &  par_orix, par_szx, par_oriy, par_szy,     &
 &  pitau0, pdate0, pdeltat, phoriid, pfileid, domain_id)
!---------------------------------------------------------------------
!- This is just an interface for histbeg_regular in case when
!- the user provides plon and plat as vectors. Obviously this can only
!- be used for very regular grids.
!-
!- INPUT
!-
!- pfilename : Name of the netcdf file to be created
!- pim       : Size of arrays in longitude direction
!- plon      : Coordinates of points in longitude
!- pjm       : Size of arrays in latitude direction
!- plat      : Coordinates of points in latitude
!-
!- The next 4 arguments allow to define a horizontal zoom
!- for this file. It is assumed that all variables to come
!- have the same index space. This can not be assumed for
!- the z axis and thus we define the zoom in histdef.
!-
!- par_orix  : Origin of the slab of data within the X axis (pim)
!- par_szx   : Size of the slab of data in X
!- par_oriy  : Origin of the slab of data within the Y axis (pjm)
!- par_szy   : Size of the slab of data in Y
!-
!- pitau0    : time step at which the history tape starts
!- pdate0    : The Julian date at which the itau was equal to 0
!- pdeltat   : Time step in seconds. Time step of the counter itau
!-             used in histwrt for instance
!-
!- OUTPUT
!-
!- phoriid   : ID of the horizontal grid
!- pfileid   : ID of the netcdf file
!-
!- Optional INPUT arguments
!-
!- domain_id       : Domain identifier
!-
!- TO DO
!-
!- This package should be written in f90
!- and use the following features :
!-    - structures for the meta-data of the files and variables
!-    - memory allocation as needed
!-    - Pointers
!-
!- VERSION
!-
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*) :: pfilename
  INTEGER,INTENT(IN) :: pim,pjm
  REAL,DIMENSION(pim),INTENT(IN) :: plon
  REAL,DIMENSION(pjm),INTENT(IN) :: plat
  INTEGER,INTENT(IN):: par_orix,par_szx,par_oriy,par_szy
  INTEGER,INTENT(IN) :: pitau0
  REAL,INTENT(IN) :: pdate0,pdeltat
  INTEGER,INTENT(OUT) :: pfileid,phoriid
  INTEGER,INTENT(IN),OPTIONAL :: domain_id
!-
  LOGICAL :: check = .FALSE.
!-
  REAL,ALLOCATABLE,DIMENSION(:,:) :: lon_tmp,lat_tmp
!---------------------------------------------------------------------
  IF (check) WRITE(*,*) "histbeg_totreg"
!-
  ALLOCATE (lon_tmp(pim,pjm),lat_tmp(pim,pjm))
!-
  lon_tmp(:,:) = SPREAD(plon(:),2,pjm)
  lat_tmp(:,:) = SPREAD(plat(:),1,pim)
!-
  CALL histbeg_regular &
 &  (pfilename,pim,lon_tmp,pjm,lat_tmp, &
 &   par_orix,par_szx,par_oriy,par_szy, &
 &   pitau0,pdate0,pdeltat,phoriid,pfileid, &
 &   .TRUE.,domain_id)
!-
  DEALLOCATE (lon_tmp, lat_tmp)
!----------------------------
END SUBROUTINE histbeg_totreg
!===
SUBROUTINE histbeg_regular &
 & (pfilename, pim, plon, pjm, plat,         &
 &  par_orix, par_szx, par_oriy, par_szy,    &
 &  pitau0, pdate0, pdeltat, phoriid, pfileid, &
 &  opt_rectilinear, domain_id)
!---------------------------------------------------------------------
!- This subroutine initializes a netcdf file and returns the ID.
!- It will set up the geographical space on which the data will be
!- stored and offers the possibility of seting a zoom.
!- It also gets the global parameters into the I/O subsystem.
!-
!- INPUT
!-
!- pfilename : Name of the netcdf file to be created
!- pim       : Size of arrays in longitude direction
!- plon      : Coordinates of points in longitude
!- pjm       : Size of arrays in latitude direction
!- plat      : Coordinates of points in latitude
!-
!- The next 4 arguments allow to define a horizontal zoom
!- for this file. It is assumed that all variables to come
!- have the same index space. This can not be assumed for
!- the z axis and thus we define the zoom in histdef.
!-
!- par_orix  : Origin of the slab of data within the X axis (pim)
!- par_szx   : Size of the slab of data in X
!- par_oriy  : Origin of the slab of data within the Y axis (pjm)
!- par_szy   : Size of the slab of data in Y
!-
!- pitau0    : time step at which the history tape starts
!- pdate0    : The Julian date at which the itau was equal to 0
!- pdeltat   : Time step in seconds. Time step of the counter itau
!-             used in histwrt for instance
!-
!- OUTPUT
!-
!- phoriid   : ID of the horizontal grid
!- pfileid   : ID of the netcdf file
!-
!- Optional INPUT arguments
!-
!- opt_rectilinear : If true we know the grid is rectilinear
!- domain_id       : Domain identifier
!-
!- TO DO
!-
!- This package should be written in F90 and use the following
!- feature :
!-   - structures for the meta-data of the files and variables
!-   - memory allocation as needed
!-   - Pointers
!-
!- VERSION
!-
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*) :: pfilename
  INTEGER,INTENT(IN) :: pim,pjm
  REAL,DIMENSION(pim,pjm),INTENT(IN) :: plon,plat
  INTEGER,INTENT(IN):: par_orix, par_szx, par_oriy, par_szy
  INTEGER,INTENT(IN) :: pitau0
  REAL,INTENT(IN) :: pdate0, pdeltat
  INTEGER,INTENT(OUT) :: pfileid, phoriid
  LOGICAL,INTENT(IN),OPTIONAL :: opt_rectilinear
  INTEGER,INTENT(IN),OPTIONAL :: domain_id
!-
  INTEGER :: ncid, iret
  INTEGER :: lengf, lenga
  CHARACTER(LEN=120) :: file
  CHARACTER(LEN=30) :: timenow
  LOGICAL :: rectilinear
!-
  LOGICAL :: check = .FALSE.
!---------------------------------------------------------------------
  nb_files = nb_files+1
  pfileid  = nb_files
!-
! 1.0 Transfering into the common for future use
!-
  IF (check) WRITE(*,*) "histbeg_regular 1.0"
!-
  itau0(pfileid) = pitau0
  date0(pfileid) = pdate0
  deltat(pfileid) = pdeltat
!-
  IF (PRESENT(opt_rectilinear)) THEN
    rectilinear = opt_rectilinear
  ELSE
    rectilinear = .FALSE.
  ENDIF
!-
! 2.0 Initializes all variables for this file
!-
  IF (check) WRITE(*,*) "histbeg_regular 2.0"
!-
  IF (nb_files > nb_files_max) THEN
    CALL ipslerr (3,"histbeg", &
   &  'Table of files too small. You should increase nb_files_max', &
   &  'in M_HISTCOM.f90 in order to accomodate all these files', ' ')
  ENDIF
!-
  nb_var(pfileid) = 0
  nb_tax(pfileid) = 0
  nb_hax(pfileid) = 0
  nb_zax(pfileid) = 0
!-
  slab_ori(pfileid,1:2) = (/ par_orix, par_oriy /)
  slab_sz(pfileid,1:2)  = (/ par_szx,  par_szy  /)
!-
! 3.0 Opening netcdf file and defining dimensions
!-
  IF (check) WRITE(*,*) "histbeg_regular 3.0"
!-
! Add DOMAIN number and ".nc" suffix in file name if needed
!-
  file  = pfilename
  CALL flio_dom_file (file,domain_id)
!-
! Keep track of the name of the files opened
!-
  lengf=LEN_TRIM(file)
  lenga=LEN_TRIM(assc_file)
  IF (nb_files == 1) THEN
    assc_file=file(1:lengf)
  ELSE IF ( (lenga+lengf) < 500) THEN
    assc_file = assc_file(1:lenga)//' '//file(1:lengf)
  ELSE IF (     ((lenga+7) < 500) &
         & .AND.(INDEX(assc_file(1:lenga),'et.al.') < 1) ) THEN
    assc_file = assc_file(1:lenga)//' et.al.'
  ELSE
    CALL ipslerr (2,"histbeg", &
   & 'The file names do not fit into the associate_file attribute.', &
   & 'Use shorter names if you wish to keep the information.',' ')
  ENDIF
!-
  iret = NF90_CREATE (file, NF90_CLOBBER, ncid)
!-
  IF (rectilinear) THEN
    iret = NF90_DEF_DIM (ncid, 'lon', par_szx, xid(nb_files))
    iret = NF90_DEF_DIM (ncid, 'lat', par_szy, yid(nb_files))
  ELSE
    iret = NF90_DEF_DIM (ncid, 'x', par_szx, xid(nb_files))
    iret = NF90_DEF_DIM (ncid, 'y', par_szy, yid(nb_files))
  ENDIF
!-
! 4.0 Declaring the geographical coordinates and other attributes
!-
  IF (check) WRITE(*,*) "histbeg_regular 4.0"
!-
! 4.3 Global attributes
!-
  iret = NF90_PUT_ATT (ncid,NF90_GLOBAL,'Conventions','GDT 1.3')
  iret = NF90_PUT_ATT (ncid,NF90_GLOBAL,'file_name',TRIM(file))
  iret = NF90_PUT_ATT (ncid,NF90_GLOBAL,'production',TRIM(model_name))
  lock_modname = .TRUE.
  CALL ioget_timestamp (timenow)
  iret = NF90_PUT_ATT (ncid,NF90_GLOBAL,'TimeStamp',TRIM(timenow))
!-
! Add DOMAIN attributes if needed
!-
  CALL flio_dom_att (ncid,domain_id)
!-
! 5.0 Saving some important information on this file in the common
!-
  IF (check) WRITE(*,*) "histbeg_regular 5.0"
!-
  ncdf_ids(pfileid) = ncid
  full_size(pfileid,1:2) = (/ pim, pjm /)
!-
! 6.0 storing the geographical coordinates
!-
  IF ( (pim /= par_szx).OR.(pjm /= par_szy) )   zoom(pfileid)=.TRUE.
  regular(pfileid)=.TRUE.
!-
  CALL histhori_regular (pfileid, pim, plon, pjm, plat, &
 &  ' ' , 'Default grid', phoriid, rectilinear)
!-----------------------------
END SUBROUTINE histbeg_regular
!===
SUBROUTINE histbeg_irregular &
 &  (pfilename, pim, plon, plon_bounds, plat, plat_bounds,   &
 &   pitau0, pdate0, pdeltat, phoriid, pfileid, domain_id)
!---------------------------------------------------------------------
!- This subroutine initializes a netcdf file and returns the ID.
!- This version is for totaly irregular grids. In this case all
!- all the data comes in as vectors and for the grid we have
!- the coordinates of the 4 corners.
!- It also gets the global parameters into the I/O subsystem.
!-
!- INPUT
!-
!- pfilename   : Name of the netcdf file to be created
!- pim         : Size of arrays in longitude direction
!- plon        : Coordinates of points in longitude
!- plon_bounds : The 2 corners of the grid in longitude
!- plat        : Coordinates of points in latitude
!- plat_bounds : The 2 corners of the grid in latitude
!-
!- pitau0    : time step at which the history tape starts
!- pdate0    : The Julian date at which the itau was equal to 0
!- pdeltat   : Time step in seconds. Time step of the counter itau
!-             used in histwrt for instance
!-
!- OUTPUT
!-
!- phoriid   : ID of the horizontal grid
!- pfileid   : ID of the netcdf file
!-
!- Optional INPUT arguments
!-
!- domain_id       : Domain identifier
!-
!- TO DO
!-
!- This package should be written in F90 and use the following
!- feature :
!-   - structures for the meta-data of the files and variables
!-   - memory allocation as needed
!-   - Pointers
!-
!- VERSION
!-
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*) :: pfilename
  INTEGER,INTENT(IN) :: pim
  REAL,DIMENSION(pim),INTENT(IN) :: plon,plat
  REAL,DIMENSION(:,:),INTENT(IN) :: plon_bounds,plat_bounds
  INTEGER,INTENT(IN) :: pitau0
  REAL,INTENT(IN) :: pdate0, pdeltat
  INTEGER,INTENT(OUT) :: pfileid, phoriid
  INTEGER,INTENT(IN),OPTIONAL :: domain_id
!-
  INTEGER :: ncid, iret
  INTEGER :: lengf, lenga
  CHARACTER(LEN=120) :: file
  CHARACTER(LEN=30) :: timenow
!-
  LOGICAL :: check = .FALSE.
!---------------------------------------------------------------------
  nb_files = nb_files+1
  pfileid  = nb_files
!-
! 1.0 Transfering into the common for future use
!-
  IF (check) WRITE(*,*) "histbeg_irregular 1.0"
!-
  itau0(pfileid) = pitau0
  date0(pfileid) = pdate0
  deltat(pfileid) = pdeltat
!-
! 2.0 Initializes all variables for this file
!-
  IF (check) WRITE(*,*) "histbeg_irregular 2.0"
!-
  IF (nb_files > nb_files_max) THEN
    CALL ipslerr (3,"histbeg", &
   &  'Table of files too small. You should increase nb_files_max', &
   &  'in M_HISTCOM.f90 in order to accomodate all these files', ' ')
  ENDIF
!-
  nb_var(pfileid) = 0
  nb_tax(pfileid) = 0
  nb_hax(pfileid) = 0
  nb_zax(pfileid) = 0
!-
  slab_ori(pfileid,1:2) = (/ 1, 1 /)
  slab_sz(pfileid,1:2)  = (/ pim, 1 /)
!-
! 3.0 Opening netcdf file and defining dimensions
!-
  IF (check) WRITE(*,*) "histbeg_irregular 3.0"
!-
! Add DOMAIN number and ".nc" suffix in file name if needed
!-
  file  = pfilename
  CALL flio_dom_file (file,domain_id)
!-
! Keep track of the name of the files opened
!-
  lengf=LEN_TRIM(file)
  lenga=LEN_TRIM(assc_file)
  IF (nb_files == 1) THEN
    assc_file=file(1:lengf)
  ELSE IF ( (lenga+lengf) < 500) THEN
    assc_file = assc_file(1:lenga)//' '//file(1:lengf)
  ELSE IF (     ((lenga+7) < 500) &
         & .AND.(INDEX(assc_file(1:lenga),'et.al.') < 1) ) THEN
    assc_file = assc_file(1:lenga)//' et.al.'
  ELSE
    CALL ipslerr (2,"histbeg", &
   & 'The file names do not fit into the associate_file attribute.', &
   & 'Use shorter names if you wish to keep the information.',' ')
  ENDIF
!-
  iret = NF90_CREATE (file, NF90_CLOBBER, ncid)
!-
  iret = NF90_DEF_DIM (ncid, 'x', pim, xid(nb_files))
  yid(nb_files) = 0
!-
!- 4.0 Declaring the geographical coordinates and other attributes
!-
   IF (check) WRITE(*,*) "histbeg_irregular 4.0"
!-
! 4.3 Global attributes
!-
  iret = NF90_PUT_ATT (ncid,NF90_GLOBAL,'Conventions','GDT 1.3')
  iret = NF90_PUT_ATT (ncid,NF90_GLOBAL,'file_name',TRIM(file))
  iret = NF90_PUT_ATT (ncid,NF90_GLOBAL,'production',TRIM(model_name))
  lock_modname = .TRUE.
  CALL ioget_timestamp (timenow)
  iret = NF90_PUT_ATT (ncid,NF90_GLOBAL,'TimeStamp',TRIM(timenow))
!-
! Add DOMAIN attributes if needed
!-
  CALL flio_dom_att (ncid,domain_id)
!-
! 5.0 Saving some important information on this file in the common
!-
  IF (check) WRITE(*,*) "histbeg_irregular 5.0"
!-
  ncdf_ids(pfileid) = ncid
  full_size(pfileid,1:2) = (/ pim, 1 /)
!-
! 6.0 storing the geographical coordinates
!-
  zoom(pfileid)=.FALSE.
  regular(pfileid)=.FALSE.
!-
  CALL histhori_irregular &
 &  (pfileid, pim, plon, plon_bounds, plat, plat_bounds, &
 &   ' ' , 'Default grid', phoriid)
!-------------------------------
END SUBROUTINE histbeg_irregular
!===
SUBROUTINE histhori_regular &
 &  (pfileid,pim,plon,pjm,plat,phname,phtitle,phid,opt_rectilinear)
!---------------------------------------------------------------------
!- This subroutine is made to declare a new horizontale grid.
!- It has to have the same number of points as
!- the original and thus in this routine we will only
!- add two variable (longitude and latitude).
!- Any variable in the file can thus point to this pair
!- through an attribute. This routine is very usefull
!- to allow staggered grids.
!-
!- INPUT
!-
!- pfileid : The id of the file to which the grid should be added
!- pim     : Size in the longitude direction
!- plon    : The longitudes
!- pjm     : Size in the latitude direction
!- plat    : The latitudes
!- phname  : The name of grid
!- phtitle : The title of the grid
!-
!- OUTPUT
!-
!- phid    : Id of the created grid
!-
!- OPTIONAL
!-
!- opt_rectilinear : If true we know the grid is rectilinear.
!-
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: pfileid, pim, pjm
  REAL,INTENT(IN),DIMENSION(pim,pjm) :: plon, plat
  CHARACTER(LEN=*),INTENT(IN) :: phname, phtitle
  INTEGER,INTENT(OUT) :: phid
  LOGICAL,INTENT(IN),OPTIONAL :: opt_rectilinear
!-
  CHARACTER(LEN=25) :: lon_name, lat_name
  CHARACTER(LEN=80) :: tmp_title, tmp_name
  INTEGER :: ndim
  INTEGER,DIMENSION(2) :: dims(2)
  INTEGER :: nlonid, nlatid
  INTEGER :: orix, oriy, par_szx, par_szy
  INTEGER :: iret, ncid
  LOGICAL :: rectilinear
!-
  LOGICAL :: check = .FALSE.
!---------------------------------------------------------------------
!-
! 1.0 Check that all fits in the buffers
!-
  IF (    (pim /= full_size(pfileid,1)) &
    & .OR.(pjm /= full_size(pfileid,2)) ) THEN
    CALL ipslerr (3,"histhori", &
   &  'The new horizontal grid does not have the same size', &
   &  'as the one provided to histbeg. This is not yet ', &
   &  'possible in the hist package.')
  ENDIF
!-
  IF (PRESENT(opt_rectilinear)) THEN
    rectilinear = opt_rectilinear
  ELSE
    rectilinear = .FALSE.
  ENDIF
!-
! 1.1 Create all the variables needed
!-
  IF (check) WRITE(*,*) "histhori_regular 1.0"
!-
  ncid = ncdf_ids(pfileid)
!-
  ndim = 2
  dims(1:2) = (/ xid(pfileid), yid(pfileid) /)
!-
  tmp_name = phname
  IF (rectilinear) THEN
    IF (nb_hax(pfileid) == 0) THEN
      lon_name = 'lon'
      lat_name = 'lat'
    ELSE
      lon_name = 'lon_'//TRIM(tmp_name)
      lat_name = 'lat_'//TRIM(tmp_name)
    ENDIF
  ELSE
    IF (nb_hax(pfileid) == 0) THEN
      lon_name = 'nav_lon'
      lat_name = 'nav_lat'
    ELSE
      lon_name = 'nav_lon_'//TRIM(tmp_name)
      lat_name = 'nav_lat_'//TRIM(tmp_name)
    ENDIF
  ENDIF
!-
! 1.2 Save the informations
!-
  phid =  nb_hax(pfileid)+1
  nb_hax(pfileid) = phid
!-
  hax_name(pfileid,phid,1:2) = (/ lon_name, lat_name /)
  tmp_title = phtitle
!-
! 2.0 Longitude
!-
  IF (check) WRITE(*,*) "histhori_regular 2.0"
!-
  IF (rectilinear) THEN
    ndim = 1
    dims(1:1) = (/ xid(pfileid) /)
  ENDIF
  iret = NF90_DEF_VAR (ncid,lon_name,NF90_FLOAT,dims(1:ndim),nlonid)
  iret = NF90_PUT_ATT (ncid,nlonid,'units',"degrees_east")
  iret = NF90_PUT_ATT (ncid,nlonid,'valid_min', &
 &                     REAL(MINVAL(plon),KIND=4))
  iret = NF90_PUT_ATT (ncid,nlonid,'valid_max', &
 &                     REAL(MAXVAL(plon),KIND=4))
  iret = NF90_PUT_ATT (ncid,nlonid,'long_name',"Longitude")
  iret = NF90_PUT_ATT (ncid,nlonid,'nav_model',TRIM(tmp_title))
!-
! 3.0 Latitude
!-
  IF (check) WRITE(*,*) "histhori_regular 3.0"
!-
  IF (rectilinear) THEN
    ndim = 1
    dims(1:1) = (/ yid(pfileid) /)
  ENDIF
  iret = NF90_DEF_VAR (ncid,lat_name,NF90_FLOAT,dims(1:ndim),nlatid)
  iret = NF90_PUT_ATT (ncid,nlatid,'units',"degrees_north")
  iret = NF90_PUT_ATT (ncid,nlatid,'valid_min', &
 &                     REAL(MINVAL(plat),KIND=4))
  iret = NF90_PUT_ATT (ncid,nlatid,'valid_max', &
 &                     REAL(MAXVAL(plat),KIND=4))
  iret = NF90_PUT_ATT (ncid,nlatid,'long_name',"Latitude")
  iret = NF90_PUT_ATT (ncid,nlatid,'nav_model',TRIM(tmp_title))
!-
  iret = NF90_ENDDEF (ncid)
!-
! 4.0 storing the geographical coordinates
!-
  IF (check) WRITE(*,*) "histhori_regular 4.0"
!-
  orix = slab_ori(pfileid,1)
  oriy = slab_ori(pfileid,2)
  par_szx = slab_sz(pfileid,1)
  par_szy = slab_sz(pfileid,2)
!-
! Transfer the longitude
!-
  IF (rectilinear) THEN
    iret = NF90_PUT_VAR (ncid,nlonid,plon(1:par_szx,1))
  ELSE
    iret = NF90_PUT_VAR (ncid,nlonid,plon, &
 &           start=(/orix,oriy/),count=(/par_szx,par_szy/))
  ENDIF
!-
! Transfer the latitude
!-
  IF ( rectilinear ) THEN
    iret = NF90_PUT_VAR (ncid,nlatid,plat(1,1:par_szy))
  ELSE
    iret = NF90_PUT_VAR (ncid,nlatid,plat, &
 &           start=(/orix,oriy/),count=(/par_szx,par_szy/))
  ENDIF
!-
  iret = NF90_REDEF (ncid)
!------------------------------
END SUBROUTINE histhori_regular
!===
SUBROUTINE histhori_irregular &
 &  (pfileid, pim, plon, plon_bounds, plat, plat_bounds, &
 &   phname, phtitle, phid)
!---------------------------------------------------------------------
!- This subroutine is made to declare a new horizontale grid.
!- It has to have the same number of points as
!- the original and thus in this routine we will only
!- add two variable (longitude and latitude).
!- Any variable in the file can thus point to this pair
!- through an attribute. This routine is very usefull
!- to allow staggered grids.
!-
!- INPUT
!-
!- pfileid     : The id of the file to which the grid should be added
!- pim         : Size in the longitude direction
!- plon        : The longitudes
!- plon_bounds : The boundaries of the grid in longitude
!- plat        : The latitudes
!- plat_bounds : Boundaries of the grid in latitude
!- phname      : The name of grid
!- phtitle     : The title of the grid
!-
!- OUTPUT
!-
!- phid    : Id of the created grid
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: pfileid, pim
  REAL,DIMENSION(pim),INTENT(IN) :: plon, plat
  REAL,DIMENSION(:,:),INTENT(IN) :: plon_bounds, plat_bounds
  CHARACTER(LEN=*), INTENT(IN) :: phname, phtitle
  INTEGER,INTENT(OUT) :: phid
!-
  CHARACTER(LEN=25) :: lon_name, lat_name
  CHARACTER(LEN=30) :: lonbound_name, latbound_name
  CHARACTER(LEN=80) :: tmp_title, tmp_name, longname
  INTEGER :: ndim, dims(2)
  INTEGER :: ndimb, dimsb(2)
  INTEGER :: nbbounds
  INTEGER :: nlonid, nlatid, nlonidb, nlatidb
  INTEGER :: iret, ncid, twoid
  LOGICAL :: transp = .FALSE.
  REAL, ALLOCATABLE, DIMENSION(:,:) :: bounds_trans
!-
  LOGICAL :: check = .FALSE.
!---------------------------------------------------------------------
!-
! 1.0 Check that all fits in the buffers
!-
  IF (    (pim /= full_size(pfileid,1)) &
    & .OR.(full_size(pfileid,2) /= 1) ) THEN
    CALL ipslerr (3,"histhori", &
   &  'The new horizontal grid does not have the same size', &
   &  'as the one provided to histbeg. This is not yet ', &
   &  'possible in the hist package.')
  ENDIF
!-
! 1.1 Create all the variables needed
!-
  IF (check) WRITE(*,*) 'histhori_irregular 1.0'
!-
  ncid = ncdf_ids(pfileid)
!-
  IF ( SIZE(plon_bounds, DIM=1) == pim ) THEN
    nbbounds = SIZE(plon_bounds, DIM=2)
    transp = .TRUE.
  ELSEIF ( SIZE(plon_bounds, DIM=2) == pim ) THEN
    nbbounds = SIZE(plon_bounds, DIM=1)
    transp = .FALSE.
  ELSE
    CALL ipslerr (3,"histhori", &
   &  'The boundary variable does not have any axis corresponding', &
   &  'to the size of the longitude or latitude variable', &
   &  '.')
  ENDIF
!-
  IF (.NOT.ALLOCATED(bounds_trans)) THEN
    ALLOCATE(bounds_trans(nbbounds,pim))
  ENDIF
!-
  iret = NF90_DEF_DIM (ncid, 'nbnd', nbbounds, twoid)
  ndim = 1
  dims(1) = xid(pfileid)
  ndimb = 2
  dimsb(1:2) = (/ twoid, xid(pfileid) /)
!-
  tmp_name = phname
  IF (nb_hax(pfileid) == 0) THEN
    lon_name = 'nav_lon'
    lat_name = 'nav_lat'
  ELSE
    lon_name = 'nav_lon_'//TRIM(tmp_name)
    lat_name = 'nav_lat_'//TRIM(tmp_name)
  ENDIF
  lonbound_name = TRIM(lon_name)//'_bounds'
  latbound_name = TRIM(lat_name)//'_bounds'
!-
! 1.2 Save the informations
!-
  phid =  nb_hax(pfileid)+1
  nb_hax(pfileid) = phid
!-
  hax_name(pfileid,phid,1:2) = (/ lon_name, lat_name /)
  tmp_title = phtitle
!-
! 2.0 Longitude
!-
  IF (check) WRITE(*,*) "histhori_irregular 2.0"
!-
  iret = NF90_DEF_VAR (ncid,lon_name,NF90_FLOAT,dims(1:ndim),nlonid)
  iret = NF90_PUT_ATT (ncid,nlonid,'units',"degrees_east")
  iret = NF90_PUT_ATT (ncid,nlonid,'valid_min', &
 &                     REAL(MINVAL(plon),KIND=4))
  iret = NF90_PUT_ATT (ncid,nlonid,'valid_max', &
 &                     REAL(MAXVAL(plon),KIND=4))
  iret = NF90_PUT_ATT (ncid,nlonid,'long_name',"Longitude")
  iret = NF90_PUT_ATT (ncid,nlonid,'nav_model',TRIM(tmp_title))
!-
! 2.1 Longitude bounds
!-
  iret = NF90_PUT_ATT (ncid,nlonid,'bounds',TRIM(lonbound_name))
  iret = NF90_DEF_VAR (ncid,lonbound_name,NF90_FLOAT, &
 &                     dimsb(1:ndimb),nlonidb)
  longname = 'Boundaries for coordinate variable '//TRIM(lon_name)
  iret = NF90_PUT_ATT (ncid,nlonidb,'long_name',TRIM(longname))
!-
! 3.0 Latitude
!-
  IF (check) WRITE(*,*) "histhori_irregular 3.0"
!-
  iret = NF90_DEF_VAR (ncid,lat_name,NF90_FLOAT,dims(1:ndim),nlatid)
  iret = NF90_PUT_ATT (ncid,nlatid,'units',"degrees_north")
  iret = NF90_PUT_ATT (ncid,nlatid,'valid_min', &
 &                     REAL(MINVAL(plat),KIND=4))
  iret = NF90_PUT_ATT (ncid,nlatid,'valid_max', &
 &                     REAL(MAXVAL(plat),KIND=4))
  iret = NF90_PUT_ATT (ncid,nlatid,'long_name',"Latitude")
  iret = NF90_PUT_ATT (ncid,nlatid,'nav_model',TRIM(tmp_title))
!-
! 3.1 Latitude bounds
!-
  iret = NF90_PUT_ATT (ncid,nlatid,'bounds',TRIM(latbound_name))
  iret = NF90_DEF_VAR (ncid,latbound_name,NF90_FLOAT, &
 &                     dimsb(1:ndimb),nlatidb)
  longname = 'Boundaries for coordinate variable '//TRIM(lat_name)
  iret = NF90_PUT_ATT (ncid,nlatidb,'long_name',TRIM(longname))
!-
  iret = NF90_ENDDEF (ncid)
!-
! 4.0 storing the geographical coordinates
!-
  IF (check) WRITE(*,*) "histhori_irregular 4.0"
!-
! 4.1 Write the longitude
!-
  iret = NF90_PUT_VAR (ncid, nlonid, plon(1:pim))
!-
! 4.2 Write the longitude bounds
!-
  IF (transp) THEN
    bounds_trans = TRANSPOSE(plon_bounds)
  ELSE
    bounds_trans = plon_bounds
  ENDIF
  iret = NF90_PUT_VAR (ncid, nlonidb, bounds_trans(1:nbbounds,1:pim))
!-
! 4.3 Write the latitude
!-
  iret = NF90_PUT_VAR (ncid, nlatid, plat(1:pim))
!-
! 4.4 Write the latitude bounds
!-
  IF (transp) THEN
    bounds_trans = TRANSPOSE(plat_bounds)
  ELSE
    bounds_trans = plat_bounds
  ENDIF
  iret = NF90_PUT_VAR (ncid,nlatidb,bounds_trans(1:nbbounds,1:pim))
!-
  iret = NF90_REDEF (ncid)
!--------------------------------
END SUBROUTINE histhori_irregular
!===
SUBROUTINE histvert (pfileid, pzaxname, pzaxtitle, &
 &                   pzaxunit, pzsize, pzvalues, pzaxid, pdirect)
!---------------------------------------------------------------------
!- This subroutine defines a vertical axis and returns it s id.
!- It gives the user the possibility to the user to define many
!- different vertical axes. For each variable defined with histdef a
!- vertical axis can be specified with by it s ID.
!-
!- INPUT
!-
!- pfileid  : ID of the file the variable should be archived in
!- pzaxname : Name of the vertical axis
!- pzaxtitle: title of the vertical axis
!- pzaxunit : Units of the vertical axis
!- pzsize   : size of the vertical axis
!- pzvalues : Coordinate values of the vetical axis
!-
!- pdirect  : is an optional argument which allows to specify the
!-            the positive direction of the axis : up or down.
!- OUTPUT
!-
!- pzaxid   : Returns the ID of the axis.
!-            Note that this is not the netCDF ID !
!-
!- VERSION
!-
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: pfileid,pzsize
  CHARACTER(LEN=*),INTENT(IN) :: pzaxname,pzaxunit,pzaxtitle
  REAL,INTENT(IN) :: pzvalues(pzsize)
  INTEGER, INTENT(OUT) :: pzaxid
  CHARACTER(LEN=*),INTENT(IN), OPTIONAL :: pdirect
!-
  INTEGER :: pos, iv, nb, zdimid, zaxid_tmp
  CHARACTER(LEN=20) :: str20, tab_str20(nb_zax_max)
  INTEGER tab_str20_length(nb_zax_max)
  CHARACTER(LEN=70) :: str70, str71, str72
  CHARACTER(LEN=80) :: str80
  CHARACTER(LEN=20) :: direction
  INTEGER :: iret, leng, ncid
  LOGICAL :: check = .FALSE.
!---------------------------------------------------------------------
!-
! 1.0 Verifications :
!    Do we have enough space for an extra axis ?
!    Is the name already in use ?
!-
  IF (check) WRITE(*,*) "histvert : 1.0 Verifications", &
 &                      pzaxname,'---',pzaxunit,'---',pzaxtitle
!-
! - Direction of axis. Can we get if from the user.
!   If not we put unknown.
!-
  IF (PRESENT(pdirect)) THEN
    direction = TRIM(pdirect)
    CALL strlowercase (direction)
  ELSE
    direction = 'unknown'
  ENDIF
!-
! Check the consistency of the attribute
!-
  IF (     (direction /= 'unknown') &
 &    .AND.(direction /= 'up')      &
 &    .AND.(direction /= 'down')   ) THEN
    direction = 'unknown'
    str80 = 'The specified axis was : '//TRIM(direction)
    CALL ipslerr (2,"histvert",&
   & "The specified direction for the vertical axis is not possible.",&
   & "it is replaced by : unknown", str80)
  ENDIF
!-
  IF ( nb_zax(pfileid)+1 > nb_zax_max) THEN
    CALL ipslerr (3,"histvert", &
   &  'Table of vertical axes too small. You should increase ',&
   &  'nb_zax_max in M_HISTCOM.f90 in order to accomodate all ', &
   &  'these variables ')
  ENDIF
!-
  iv = nb_zax(pfileid)
  IF ( iv > 1) THEN
    str20 = pzaxname
    nb = iv-1
    tab_str20(1:nb) = zax_name(pfileid,1:nb)
    tab_str20_length(1:nb) = zax_name_length(pfileid,1:nb)
    CALL find_str(nb, tab_str20, tab_str20_length, str20, pos)
  ELSE
    pos = 0
  ENDIF
!-
  IF ( pos > 0) THEN
    str70 = "Vertical axis already exists"
    WRITE(str71,'("Check variable ",a," in file",I3)') str20,pfileid
    str72 = "Can also be a wrong file ID in another declaration"
    CALL ipslerr (3,"histvert", str70, str71, str72)
  ENDIF
!-
  iv = nb_zax(pfileid)+1
!-
! 2.0 Add the information to the file
!-
  IF (check) &
 &  WRITE(*,*) "histvert : 2.0 Add the information to the file"
!-
  ncid = ncdf_ids(pfileid)
!-
  leng = MIN(LEN_TRIM(pzaxname),20)
  iret = NF90_DEF_DIM (ncid,pzaxname(1:leng),pzsize,zaxid_tmp)
  iret = NF90_DEF_VAR (ncid,pzaxname(1:leng),NF90_FLOAT, &
 &                     zaxid_tmp,zdimid)
!-
  leng = MIN(LEN_TRIM(pzaxunit),20)
  iret = NF90_PUT_ATT (ncid, zdimid, 'units', pzaxunit(1:leng))
  iret = NF90_PUT_ATT (ncid, zdimid, 'positive', TRIM(direction))
!-
  iret = NF90_PUT_ATT (ncid,zdimid,'valid_min', &
 &                     REAL(MINVAL(pzvalues(1:pzsize)),KIND=4))
  iret = NF90_PUT_ATT (ncid,zdimid,'valid_max', &
 &                     REAL(MAXVAL(pzvalues(1:pzsize)),KIND=4))
!-
  leng = MIN(LEN_TRIM(pzaxname),20)
  iret = NF90_PUT_ATT (ncid, zdimid, 'title', pzaxname(1:leng))
  leng = MIN(LEN_TRIM(pzaxtitle),80)
  iret = NF90_PUT_ATT (ncid, zdimid, 'long_name', pzaxtitle(1:leng))
!-
  iret = NF90_ENDDEF (ncid)
!-
  iret = NF90_PUT_VAR (ncid, zdimid, pzvalues(1:pzsize))
!-
  iret = NF90_REDEF (ncid)
!-
!- 3.0 add the information to the common
!-
  IF ( check) &
  &  WRITE(*,*) "histvert : 3.0 add the information to the common"
!-
  nb_zax(pfileid) = iv
  zax_size(pfileid, iv) = pzsize
  zax_name(pfileid, iv) = pzaxname
  zax_name_length(pfileid, iv) = LEN_TRIM(pzaxname)
  zax_ids(pfileid, iv) = zaxid_tmp
  pzaxid =  iv
!----------------------
END SUBROUTINE histvert
!===
SUBROUTINE histdef (pfileid, pvarname, ptitle, punit, &
 &                  pxsize, pysize, phoriid, pzsize,  &
 &                  par_oriz, par_szz, pzid,          &
 &                  pnbbyt, popp, pfreq_opp, pfreq_wrt)
!---------------------------------------------------------------------
!- With this subroutine each variable to be archived on the history
!- tape should be declared.
!-
!- It gives the user the choise of operation
!- to be performed on the variables, the frequency of this operation
!- and finaly the frequency of the archiving.
!-
!- INPUT
!-
!- pfileid  : ID of the file the variable should be archived in
!- pvarname : Name of the variable, short and easy to remember
!- ptitle   : Full name of the variable
!- punit    : Units of the variable
!-
!- The next 3 arguments give the size of that data
!- that will be passed to histwrite. The zoom will be
!- done there with the horizontal information obtained
!- in histbeg and the vertical information to follow.
!-
!- pxsize   : Size in X direction (size of the data that will be
!-            given to histwrite)
!- pysize   : Size in Y direction
!- phoriid  : ID of the horizontal axis
!-
!- The next two arguments give the vertical zoom to use.
!-
!- pzsize   : Size in Z direction (If 1 then no axis is declared
!-            for this variable and pzid is not used)
!- par_oriz : Off set of the zoom
!- par_szz  : Size of the zoom
!-
!- pzid     : ID of the vertical axis to use. It has to have
!-            the size of the zoom.
!- pnbbyt   : Number of bytes on which to store in netCDF (Not opp.)
!- popp     : Operation to be performed. The following options
!-            exist today :
!- inst : keeps instantaneous values for writting
!- ave  : Computes the average from call between writes
!- pfreq_opp: Frequency of this operation (in seconds)
!- pfreq_wrt: Frequency at which the variable should be
!-            written (in seconds)
!-
!- VERSION
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: pfileid, pxsize, pysize, pzsize, pzid
  INTEGER,INTENT(IN) :: par_oriz, par_szz, pnbbyt, phoriid
  CHARACTER(LEN=*),INTENT(IN) :: pvarname, punit, popp
  CHARACTER(LEN=*),INTENT(IN) :: ptitle
  REAL,INTENT(IN) :: pfreq_opp, pfreq_wrt
!-
  INTEGER :: iv, i, nb
  CHARACTER(LEN=70) :: str70, str71, str72
  CHARACTER(LEN=20) :: tmp_name
  CHARACTER(LEN=20) :: str20, tab_str20(nb_var_max)
  INTEGER :: tab_str20_length(nb_var_max)
  CHARACTER(LEN=40) :: str40, tab_str40(nb_var_max)
  INTEGER :: tab_str40_length(nb_var_max)
  CHARACTER(LEN=10) :: str10
  CHARACTER(LEN=80) :: tmp_str80
  CHARACTER(LEN=7) :: tmp_topp, tmp_sopp(nbopp_max)
  CHARACTER(LEN=120) :: ex_topps
  REAL :: tmp_scal(nbopp_max), un_an, un_jour, test_fopp, test_fwrt
  INTEGER :: pos, buff_sz
!-
  LOGICAL :: check = .FALSE.
!---------------------------------------------------------------------
  ex_topps = 'ave, inst, t_min, t_max, t_sum, once, never, l_max, l_min'
!-
  nb_var(pfileid) = nb_var(pfileid)+1
  iv = nb_var(pfileid)
!-
  IF ( iv > nb_var_max) THEN
    CALL ipslerr (3,"histdef", &
   &  'Table of variables too small. You should increase nb_var_max',&
   &  'in M_HISTCOM.f90 in order to accomodate all these variables', &
   &  ' ')
  ENDIF
!-
! 1.0 Transfer informations on the variable to the common
!     and verify that it does not already exist
!-
  IF (check) WRITE(*,*) "histdef : 1.0"
!-
  IF (iv > 1) THEN
    str20 = pvarname
    nb = iv-1
    tab_str20(1:nb) = name(pfileid,1:nb)
    tab_str20_length(1:nb) = name_length(pfileid,1:nb)
    CALL find_str (nb, tab_str20, tab_str20_length, str20, pos)
  ELSE
    pos = 0
  ENDIF
!-
  IF (pos > 0) THEN
    str70 = "Variable already exists"
    WRITE(str71,'("Check variable  ",a," in file",I3)') str20,pfileid
    str72 = "Can also be a wrong file ID in another declaration"
    CALL ipslerr (3,"histdef", str70, str71, str72)
  ENDIF
!-
  name(pfileid,iv) = pvarname
  name_length(pfileid,iv) = LEN_TRIM(name(pfileid,iv))
  title(pfileid,iv) = ptitle
  unit_name(pfileid,iv) = punit
  tmp_name =  name(pfileid,iv)
!-
! 1.1 decode the operations
!-
  fullop(pfileid,iv) = popp
  tmp_str80 = popp
  CALL buildop &
 &  (tmp_str80, ex_topps, tmp_topp, nbopp_max, missing_val, &
 &   tmp_sopp, tmp_scal, nbopp(pfileid,iv))
!-
  topp(pfileid,iv) = tmp_topp
  DO i=1,nbopp(pfileid,iv)
    sopps(pfileid,iv,i) = tmp_sopp(i)
    scal(pfileid,iv,i) = tmp_scal(i)
  ENDDO
!-
! 1.2 If we have an even number of operations
!     then we need to add identity
!-
  IF (2*INT(nbopp(pfileid,iv)/2.0) == nbopp(pfileid,iv)) THEN
    nbopp(pfileid,iv) = nbopp(pfileid,iv)+1
    sopps(pfileid,iv,nbopp(pfileid,iv)) = 'ident'
    scal(pfileid,iv,nbopp(pfileid,iv)) = missing_val
  ENDIF
!-
! 2.0 Put the size of the variable in the common and check
!-
  IF (check) &
 &  WRITE(*,*) "histdef : 2.0", pfileid,iv,nbopp(pfileid,iv), &
 &    sopps(pfileid,iv,1:nbopp(pfileid,iv)), &
 &    scal(pfileid,iv,1:nbopp(pfileid,iv))
!-
  scsize(pfileid,iv,1:3) = (/ pxsize, pysize, pzsize /)
!-
  zorig(pfileid,iv,1:3) = &
 &  (/ slab_ori(pfileid,1), slab_ori(pfileid,2), par_oriz /)
!-
  zsize(pfileid,iv,1:3) = &
 &  (/ slab_sz(pfileid,1), slab_sz(pfileid,2), par_szz /)
!-
! Is the size of the full array the same as that of the coordinates  ?
!-
  IF (    (pxsize > full_size(pfileid,1)) &
 &    .OR.(pysize > full_size(pfileid,2)) ) THEN
!-
    str70 = "The size of the variable is different "// &
 &          "from the one of the coordinates"
    WRITE(str71,'("Size of coordinates :", 2I4)') &
 &   full_size(pfileid,1), full_size(pfileid,2)
    WRITE(str72,'("Size declared for variable ",a," :",2I4)') &
 &   TRIM(tmp_name), pxsize, pysize
    CALL ipslerr (3,"histdef", str70, str71, str72)
  ENDIF
!-
! Is the size of the zoom smaler than the coordinates ?
!-
  IF (    (full_size(pfileid,1) < slab_sz(pfileid,1)) &
 &    .OR.(full_size(pfileid,2) < slab_sz(pfileid,2)) ) THEN
    str70 = &
 &   "Size of variable should be greater or equal to those of the zoom"
    WRITE(str71,'("Size of XY zoom :", 2I4)') &
 &   slab_sz(pfileid,1),slab_sz(pfileid,1)
    WRITE(str72,'("Size declared for variable ",a," :",2I4)') &
 &   TRIM(tmp_name), pxsize, pysize
    CALL ipslerr (3,"histdef", str70, str71, str72)
  ENDIF
!-
! 2.1 We store the horizontal grid information with minimal
!     and a fall back onto the default grid
!-
  IF ( phoriid > 0 .AND. phoriid <= nb_hax(pfileid)) THEN
    var_haxid(pfileid,iv) = phoriid
  ELSE
    var_haxid(pfileid,iv) = 1
    CALL ipslerr (2,"histdef", &
   &  'We use the default grid for variable as an invalide',&
   &  'ID was provided for variable : ', pvarname)
  ENDIF
!-
! 2.2 Check the vertical coordinates if needed
!-
  IF (par_szz > 1) THEN
!-
!-- Does the vertical coordinate exist ?
!-
    IF ( pzid > nb_zax(pfileid)) THEN
      WRITE(str70, &
 &    '("The vertical coordinate chosen for variable ",a)') &
 &     TRIM(tmp_name)
      str71 = " Does not exist."
      CALL ipslerr (3,"histdef",str70,str71, " ")
    ENDIF
!-
!-- Is the vertical size of the variable equal to that of the axis ?
!-
    IF (par_szz /= zax_size(pfileid,pzid)) THEN
      str20 = zax_name(pfileid,pzid)
      str70 = "The size of the zoom does not correspond "// &
 &            "to the size of the chosen vertical axis"
      WRITE(str71,'("Size of zoom in z :", I4)') par_szz
      WRITE(str72,'("Size declared for axis ",a," :",I4)') &
 &     TRIM(str20), zax_size(pfileid,pzid)
      CALL ipslerr (3,"histdef", str70, str71, str72)
    ENDIF
!-
!-- Is the zoom smaler that the total size of the variable ?
!-
    IF ( pzsize < par_szz ) THEN
      str20 = zax_name(pfileid,pzid)
      str70 = "The vertical size of variable "// &
 &            "is smaller than that of the zoom."
      WRITE(str71,'("Declared vertical size of data :", I5)') pzsize
      WRITE(str72,'("Size of zoom for variable ",a," = ",I5)') &
 &     TRIM(tmp_name),par_szz
      CALL ipslerr (3,"histdef", str70, str71, str72)
    ENDIF
    var_zaxid(pfileid,iv) = pzid
  ELSE
    var_zaxid(pfileid,iv) = -99
  ENDIF
!-
! 3.0 Determine the position of the variable in the buffer
!     If it is instantaneous output then we do not use the buffer
!-
  IF (check) WRITE(*,*) "histdef : 3.0"
!-
! 3.1 We get the size of the arrays histwrite will get and check
!     that they fit into the tmp_buffer
!-
  buff_sz = zsize(pfileid,iv,1)*zsize(pfileid,iv,2)*zsize(pfileid,iv,3)
!-
! 3.2 move the pointer of the buffer array for operation
!     which need bufferisation
!-
  IF (     (TRIM(tmp_topp) /= "inst") &
 &    .AND.(TRIM(tmp_topp) /= "once") &
 &    .AND.(TRIM(tmp_topp) /= "never") )THEN
    point(pfileid,iv) = buff_pos+1
    buff_pos = buff_pos+buff_sz
    IF (check) THEN
      WRITE(*,*) "histdef : 3.2 bufpos for iv = ",iv, &
 &               " pfileid = ",pfileid," is = ",point(pfileid,iv)
    ENDIF
  ENDIF
!-
! 4.0 Transfer the frequency of the operations and check
!     for validity. We have to pay attention to negative values
!     of the frequency which indicate monthly time-steps.
!     The strategy is to bring it back to seconds for the tests
!-
  IF (check) WRITE(*,*) "histdef : 4.0"
!-
  freq_opp(pfileid,iv) = pfreq_opp
  freq_wrt(pfileid,iv) = pfreq_wrt
!-
  CALL ioget_calendar(un_an, un_jour)
  IF ( pfreq_opp < 0) THEN
    CALL ioget_calendar(un_an)
    test_fopp = pfreq_opp*(-1.)*un_an/12.*un_jour
  ELSE
    test_fopp = pfreq_opp
  ENDIF
  IF ( pfreq_wrt < 0) THEN
    CALL ioget_calendar(un_an)
    test_fwrt = pfreq_wrt*(-1.)*un_an/12.*un_jour
  ELSE
    test_fwrt = pfreq_wrt
  ENDIF
!-
! 4.1 Frequency of operations and output should be larger than deltat !
!-
  IF (test_fopp < deltat(pfileid)) THEN
    str70 = 'Frequency of operations should be larger than deltat'
    WRITE(str71,'("It is not the case for variable ",a," :",F10.4)') &
 &   TRIM(tmp_name),pfreq_opp
    str72 = "PATCH : frequency set to deltat"
!-
    CALL ipslerr (2,"histdef", str70, str71, str72)
!-
    freq_opp(pfileid,iv) = deltat(pfileid)
  ENDIF
!-
  IF (test_fwrt < deltat(pfileid)) THEN
    str70 = 'Frequency of output should be larger than deltat'
    WRITE(str71,'("It is not the case for variable ",a," :",F10.4)') &
 &   TRIM(tmp_name),pfreq_wrt
    str72 = "PATCH : frequency set to deltat"
!-
    CALL ipslerr (2,"histdef", str70, str71, str72)
!-
    freq_wrt(pfileid,iv) = deltat(pfileid)
  ENDIF
!-
! 4.2 First the existence of the operation is tested and then
!     its compaticility with the choice of frequencies
!-
  IF (TRIM(tmp_topp) == "inst") THEN
    IF (test_fopp /= test_fwrt) THEN
      str70 = 'For instantaneous output the frequency '// &
 &            'of operations and output'
      WRITE(str71, &
 &     '("should be the same, this was not case for variable ",a)') &
 &      TRIM(tmp_name)
      str72 = "PATCH : The smalest frequency of both is used"
      CALL ipslerr (2,"histdef", str70, str71, str72)
      IF ( test_fopp < test_fwrt) THEN
        freq_opp(pfileid,iv) = pfreq_opp
        freq_wrt(pfileid,iv) = pfreq_opp
      ELSE
        freq_opp(pfileid,iv) = pfreq_wrt
        freq_wrt(pfileid,iv) = pfreq_wrt
      ENDIF
    ENDIF
  ELSE IF (INDEX(ex_topps,TRIM(tmp_topp)) > 0) THEN
    IF (test_fopp > test_fwrt) THEN
      str70 = 'For averages the frequency of operations '// &
&             'should be smaller or equal'
      WRITE(str71, &
 &     '("to that of output. It is not the case for variable ",a)') &
 &     TRIM(tmp_name)
      str72 = 'PATCH : The output frequency is used for both'
      CALL ipslerr (2,"histdef", str70, str71, str72)
      freq_opp(pfileid,iv) = pfreq_wrt
    ENDIF
  ELSE
    WRITE (str70,'("Operation on variable ",a," is unknown")') &
&    TRIM(tmp_name)
    WRITE (str71, '("operation requested is :",a)') tmp_topp
    WRITE (str72, '("File ID :",I3)') pfileid
    CALL ipslerr (3,"histdef", str70, str71, str72)
  ENDIF
!-
! 5.0 Initialize other variables of the common
!-
  IF (check) WRITE(*,*) "histdef : 5.0"
!-
  minmax(pfileid,iv,1) =  ABS(missing_val)
  minmax(pfileid,iv,2) = -ABS(missing_val)
!-
  last_opp(pfileid,iv) = itau0(pfileid)       ! - freq_opp(pfileid,iv)/2./deltat(pfileid)
  last_wrt(pfileid,iv) = itau0(pfileid)       ! - freq_wrt(pfileid,iv)/2./deltat(pfileid)
  last_opp_chk(pfileid,iv) = itau0(pfileid)   ! - freq_opp(pfileid,iv)/2./deltat(pfileid)
  last_wrt_chk(pfileid,iv) = itau0(pfileid)   ! - freq_wrt(pfileid,iv)/2./deltat(pfileid)
  nb_opp(pfileid,iv) = 0
  nb_wrt(pfileid,iv) = 0
!-
! 6.0 Get the time axis for this variable
!-
  IF (check) WRITE(*,*) "histdef : 6.0"
!-
  IF ( freq_wrt(pfileid,iv) > 0 ) THEN
    WRITE(str10,'(I8.8)') INT(freq_wrt(pfileid,iv))
    str40 = TRIM(tmp_topp)//"_"//TRIM(str10)
  ELSE
    WRITE(str10,'(I2.2,"month")') ABS(INT(freq_wrt(pfileid,iv)))
    str40 = TRIM(tmp_topp)//"_"//TRIM(str10)
  ENDIF
!-
  DO i=1,nb_tax(pfileid)
    tab_str40(i) = tax_name(pfileid,i)
    tab_str40_length(i) = tax_name_length(pfileid,i)
  ENDDO
!-
  CALL find_str (nb_tax(pfileid),tab_str40,tab_str40_length,str40,pos)
!-
! No time axis for once, l_max, l_min or never operation
!-
  IF (     (TRIM(tmp_topp) /= 'once')  &
 &    .AND.(TRIM(tmp_topp) /= 'never') &
 &    .AND.(TRIM(tmp_topp) /= 'l_max') &
 &    .AND.(TRIM(tmp_topp) /= 'l_min') ) THEN
    IF ( pos < 0) THEN
      nb_tax(pfileid) = nb_tax(pfileid)+1
      tax_name(pfileid,nb_tax(pfileid)) = str40
      tax_name_length(pfileid, nb_tax(pfileid)) = LEN_TRIM(str40)
      tax_last(pfileid,nb_tax(pfileid)) = 0
      var_axid(pfileid,iv) = nb_tax(pfileid)
    ELSE
      var_axid(pfileid,iv) = pos
    ENDIF
  ELSE
    IF (check)   WRITE(*,*) "histdef : 7.0 ",TRIM(tmp_topp),'----'
    var_axid(pfileid,iv) = -99
  ENDIF
!-
! 7.0 prepare frequence of writing and operation
!     for never or once operation
!-
  IF (    (TRIM(tmp_topp) == 'once')  &
 &    .OR.(TRIM(tmp_topp) == 'never') ) THEN
    freq_opp(pfileid,iv) = 0.
    freq_wrt(pfileid,iv) = 0.
  ENDIF
!---------------------
END SUBROUTINE histdef
!===
SUBROUTINE histend (pfileid)
!---------------------------------------------------------------------
!- This subroutine end the decalaration of variables and sets the
!- time axes in the netcdf file and puts it into the write mode.
!-
!- INPUT
!-
!- pfileid : ID of the file to be worked on
!-
!- VERSION
!-
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER, INTENT(IN) :: pfileid
!-
  INTEGER :: ncid, ncvarid
  INTEGER :: iret, ndim, iv, itx, ziv
  INTEGER :: itax
  INTEGER :: dims(4), dim_cnt
  INTEGER :: year, month, day, hours, minutes
  REAL :: sec
  REAL :: rtime0
  CHARACTER(LEN=20) :: tname, tunit
  CHARACTER(LEN=30) :: str30
  CHARACTER(LEN=80) :: ttitle
  CHARACTER(LEN=120) :: assoc
  CHARACTER(LEN=70) :: str70
  CHARACTER(LEN=3),DIMENSION(12) :: cal =   &
 &  (/ 'JAN','FEB','MAR','APR','MAY','JUN', &
 &     'JUL','AUG','SEP','OCT','NOV','DEC' /)
  CHARACTER(LEN=7) :: tmp_opp
!-
  LOGICAL :: check = .FALSE.
!---------------------------------------------------------------------
  ncid = ncdf_ids(pfileid)
!-
! 1.0 Create the time axes
!-
  IF (check) WRITE(*,*) "histend : 1.0"
!---
  iret = NF90_DEF_DIM (ncid,'time_counter',NF90_UNLIMITED,tid(pfileid))
!-
! 1.1 Define all the time axes needed for this file
!-
  DO itx=1,nb_tax(pfileid)
    dims(1) = tid(pfileid)
    IF (nb_tax(pfileid) > 1) THEN
      str30 = "t_"//tax_name(pfileid,itx)
    ELSE
      str30 = "time_counter"
    ENDIF
    iret = NF90_DEF_VAR (ncid,str30,NF90_FLOAT, &
 &                       dims(1),tdimid(pfileid,itx))
!---
!   To transform the current itau into a real date and take it
!   as the origin of the file requires the time counter to change.
!   Thus it is an operation the user has to ask for.
!   This function should thus only be re-instated
!   if there is a ioconf routine to control it.
!---
!-- rtime0 = itau2date(itau0(pfileid), date0(pfileid), deltat(pfileid))
    rtime0 = date0(pfileid)
!-
    CALL ju2ymds(rtime0, year, month, day, sec)
!---
!   Catch any error induced by a change in calendar !
!---
    IF (year < 0) THEN
      year = 2000+year
    ENDIF
!-
    hours = INT(sec/(60.*60.))
    minutes = INT((sec-hours*60.*60.)/60.)
    sec = sec-(hours*60.*60.+minutes*60.)
!-
    WRITE(str70,7000) year, month, day, hours, minutes, INT(sec)
    iret = NF90_PUT_ATT (ncid,tdimid(pfileid,itx),'units',TRIM(str70))
!-
    CALL ioget_calendar (str30)
    iret = NF90_PUT_ATT (ncid,tdimid(pfileid,itx), &
 &                       'calendar',TRIM(str30))
!-
    iret = NF90_PUT_ATT (ncid,tdimid(pfileid,itx),'title','Time')
!-
    iret = NF90_PUT_ATT (ncid,tdimid(pfileid,itx), &
 &                       'long_name','Time axis')
!-
    WRITE(str70,7001) year, cal(month), day, hours, minutes, INT(sec)
    iret = NF90_PUT_ATT (ncid,tdimid(pfileid,itx), &
 &                       'time_origin',TRIM(str70))
  ENDDO
!-
! The formats we need
!-
7000 FORMAT('seconds since ', I4.4,'-',I2.2,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)
7001 FORMAT(' ', I4.4,'-',A3,'-',I2.2,' ',I2.2,':',I2.2,':',I2.2)
!-
! 2.0 declare the variables
!-
  IF (check) WRITE(*,*) "histend : 2.0"
!-
  DO iv=1,nb_var(pfileid)
!---
    itax = var_axid(pfileid,iv)
!---
    tname = name(pfileid,iv)
    tunit = unit_name(pfileid,iv)
    ttitle = title(pfileid,iv)
!---
    IF ( regular(pfileid) ) THEN
      dims(1:2) = (/ xid(pfileid), yid(pfileid) /)
      dim_cnt = 2
    ELSE
      dims(1) = xid(pfileid)
      dim_cnt = 1
    ENDIF
!---
    tmp_opp = topp(pfileid,iv)
    ziv = var_zaxid(pfileid,iv)
!---
!   2.1 dimension of field
!---
    IF ( (TRIM(tmp_opp) /= 'never')) THEN
      IF (     (TRIM(tmp_opp) /= 'once')  &
     &    .AND.(TRIM(tmp_opp) /= 'l_max') &
     &    .AND.(TRIM(tmp_opp) /= 'l_min') ) THEN
        IF (ziv == -99) THEN
          ndim = dim_cnt+1
          dims(dim_cnt+1:dim_cnt+2) = (/ tid(pfileid), 0 /)
        ELSE
          ndim = dim_cnt+2
          dims(dim_cnt+1:dim_cnt+2) = (/ zax_ids(pfileid,ziv), tid(pfileid) /)
        ENDIF
      ELSE
        IF (ziv == -99) THEN
          ndim = dim_cnt
          dims(dim_cnt+1:dim_cnt+2) = (/ 0, 0 /)
        ELSE
          ndim = dim_cnt+1
          dims(dim_cnt+1:dim_cnt+2) = (/ zax_ids(pfileid,ziv), 0 /)
        ENDIF
      ENDIF
!-
      iret = NF90_DEF_VAR (ncid,TRIM(tname),NF90_FLOAT, &
 &                         dims(1:ABS(ndim)),ncvarid)
!-
      ncvar_ids(pfileid,iv) = ncvarid
!-
      iret = NF90_PUT_ATT (ncid,ncvarid,'units',TRIM(tunit))
!-
      iret = NF90_PUT_ATT (ncid,ncvarid,'missing_value', &
 &                         REAL(missing_val,KIND=4))
      iret = NF90_PUT_ATT (ncid,ncvarid,'valid_min', &
 &                         REAL(minmax(pfileid,iv,1),KIND=4))
      iret = NF90_PUT_ATT (ncid,ncvarid,'valid_max', &
 &                         REAL(minmax(pfileid,iv,2),KIND=4))
!-
      iret = NF90_PUT_ATT (ncid,ncvarid,'long_name',TRIM(ttitle))
!-
      iret = NF90_PUT_ATT (ncid,ncvarid,'short_name',TRIM(tname))
!-
      iret = NF90_PUT_ATT (ncid,ncvarid,'online_operation', &
 &                         TRIM(fullop(pfileid,iv)))
!-
      SELECT CASE(ndim)
      CASE(-3)
        str30 = 'ZYX'
      CASE(2)
        str30 = 'YX'
      CASE(3)
        str30 = 'TYX'
      CASE(4)
        str30 = 'TZYX'
      CASE DEFAULT
        CALL ipslerr (3,"histend", &
       &  'less than 2 or more than 4 dimensions are not', &
       &  'allowed at this stage',' ')
      END SELECT
!-
      iret = NF90_PUT_ATT (ncid,ncvarid,'axis',TRIM(str30))
!-
      assoc='nav_lat nav_lon'
      ziv = var_zaxid(pfileid, iv)
      IF (ziv > 0) THEN
        str30 = zax_name(pfileid,ziv)
        assoc = TRIM(str30)//' '//TRIM(assoc)
      ENDIF
!-
      IF (itax > 0) THEN
        IF ( nb_tax(pfileid) > 1) THEN
          str30 = "t_"//tax_name(pfileid,itax)
        ELSE
          str30 = "time_counter"
        ENDIF
        assoc = TRIM(str30)//' '//TRIM(assoc)
!-
        IF (check) THEN
          WRITE(*,*) "histend : 2.0.n, freq_opp, freq_wrt", &
 &                   freq_opp(pfileid,iv), freq_wrt(pfileid,iv)
        ENDIF
!-
        iret = NF90_PUT_ATT (ncid,ncvarid,'interval_operation', &
 &                           REAL(freq_opp(pfileid,iv),KIND=4))
        iret = NF90_PUT_ATT (ncid,ncvarid,'interval_write', &
 &                           REAL(freq_wrt(pfileid,iv),KIND=4))
      ENDIF
      iret = NF90_PUT_ATT (ncid,ncvarid,'associate',TRIM(assoc))
    ENDIF
  ENDDO
!-
! 3.0 Put the netcdf file into wrte mode
!-
  IF (check) WRITE(*,*) "histend : 3.0"
!-
  iret = NF90_ENDDEF (ncid)
!-
! 4.0 Give some informations to the user
!-
  IF (check) WRITE(*,*) "histend : 4.0"
!-
  WRITE(str70,'("All variables have been initialized on file :",I3)') pfileid
  CALL ipslerr (1,'histend',str70,'',' ')
!---------------------
END SUBROUTINE histend
!===
SUBROUTINE histwrite_r1d (pfileid,pvarname,pitau,pdata,nbindex,nindex)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: pfileid, pitau, nbindex, nindex(nbindex)
  REAL,DIMENSION(:),INTENT(IN) :: pdata
  CHARACTER(LEN=*),INTENT(IN) :: pvarname
!-
  LOGICAL :: do_oper, do_write, largebuf
  INTEGER :: varid, io, nbpt_in, nbpt_out
  REAL, ALLOCATABLE,SAVE :: buff_tmp(:)
  INTEGER,SAVE :: buff_tmp_sz
  CHARACTER(LEN=7) :: tmp_opp
!-
  LOGICAL :: check = .FALSE.
!---------------------------------------------------------------------
!-
! 1.0 Try to catch errors like specifying the wrong file ID.
!     Thanks Marine for showing us what errors users can make !
!-
  IF ( (pfileid < 1).OR.(pfileid > nb_files) ) THEN
    CALL ipslerr (3,"histwrite", &
 &    'Illegal file ID in the histwrite of variable',pvarname,' ')
  ENDIF
!-
! 1.1 Find the id of the variable to be written and the real time
!-
  CALL histvar_seq (pfileid,pvarname,varid)
!-
! 2.0 do nothing for never operation
!-
  tmp_opp = topp(pfileid,varid)
!-
  IF (TRIM(tmp_opp) == "never") THEN
    last_opp_chk(pfileid,varid) = -99
    last_wrt_chk(pfileid,varid) = -99
  ENDIF
!-
! 3.0 We check if we need to do an operation
!-
  IF (last_opp_chk(pfileid,varid) == pitau) THEN
    CALL ipslerr (3,"histwrite", &
 &    'This variable as already been analysed at the present', &
 &    'time step',' ')
  ENDIF
!-
  CALL isittime &
 &  (pitau, date0(pfileid), deltat(pfileid), freq_opp(pfileid,varid), &
 &   last_opp(pfileid,varid), last_opp_chk(pfileid,varid), do_oper)
!-
! 4.0 We check if we need to write the data
!-
  IF (last_wrt_chk(pfileid,varid) == pitau) THEN
    CALL ipslerr (3,"histwrite", &
 &    'This variable as already been written for the present', &
 &    'time step',' ')
  ENDIF
!-
  CALL isittime &
 &  (pitau, date0(pfileid), deltat(pfileid), freq_wrt(pfileid,varid), &
 &   last_wrt(pfileid,varid), last_wrt_chk(pfileid,varid), do_write)
!-
! 5.0 histwrite called
!-
  IF (do_oper.OR.do_write) THEN
!-
!-- 5.1 Get the sizes of the data we will handle
!-
    IF (datasz_in(pfileid,varid,1) <= 0) THEN
!---- There is the risk here that the user has over-sized the array.
!---- But how can we catch this ?
!---- In the worst case we will do impossible operations
!---- on part of the data !
      datasz_in(pfileid,varid,1) = SIZE(pdata)
      datasz_in(pfileid,varid,2) = -1
      datasz_in(pfileid,varid,3) = -1
    ENDIF
!-
!-- 5.2 The maximum size of the data will give the size of the buffer
!-
    IF (datasz_max(pfileid,varid) <= 0) THEN
      largebuf = .FALSE.
      DO io=1,nbopp(pfileid,varid)
        IF (INDEX(fuchnbout,sopps(pfileid,varid,io)) > 0) THEN
          largebuf = .TRUE.
        ENDIF
      ENDDO
      IF (largebuf) THEN
        datasz_max(pfileid,varid) = &
 &        scsize(pfileid,varid,1) &
 &       *scsize(pfileid,varid,2) &
 &       *scsize(pfileid,varid,3)
      ELSE
        datasz_max(pfileid,varid) = &
 &        datasz_in(pfileid,varid,1)
      ENDIF
    ENDIF
!-
    IF (.NOT.ALLOCATED(buff_tmp)) THEN
      IF (check) THEN
        WRITE(*,*) &
 &       "histwrite_r1d : allocate buff_tmp for buff_sz = ", &
 &       datasz_max(pfileid,varid)
      ENDIF
      ALLOCATE (buff_tmp(datasz_max(pfileid,varid)))
      buff_tmp_sz = datasz_max(pfileid,varid)
    ELSE IF (datasz_max(pfileid,varid) > buff_tmp_sz) THEN
      IF (check) THEN
        WRITE(*,*) &
 &       "histwrite_r1d : re-allocate buff_tmp for buff_sz = ", &
 &       datasz_max(pfileid,varid)
      ENDIF
      DEALLOCATE (buff_tmp)
      ALLOCATE (buff_tmp(datasz_max(pfileid,varid)))
      buff_tmp_sz = datasz_max(pfileid,varid)
    ENDIF
!-
!-- We have to do the first operation anyway.
!-- Thus we do it here and change the ranke
!-- of the data at the same time. This should speed up things.
!-
    nbpt_in = datasz_in(pfileid,varid,1)
    nbpt_out = datasz_max(pfileid,varid)
    CALL mathop (sopps(pfileid,varid,1), nbpt_in, pdata, &
 &               missing_val, nbindex, nindex, &
 &               scal(pfileid,varid,1), nbpt_out, buff_tmp)
    CALL histwrite_real (pfileid, varid, pitau, nbpt_out, &
 &            buff_tmp, nbindex, nindex, do_oper, do_write)
  ENDIF
!-
! 6.0 Manage time steps
!-
  IF ((TRIM(tmp_opp) /= "once").AND.(TRIM(tmp_opp) /= "never")) THEN
    last_opp_chk(pfileid,varid) = pitau
    last_wrt_chk(pfileid,varid) = pitau
  ELSE
    last_opp_chk(pfileid,varid) = -99
    last_wrt_chk(pfileid,varid) = -99
  ENDIF
!---------------------------
END SUBROUTINE histwrite_r1d
!===
SUBROUTINE histwrite_r2d (pfileid,pvarname,pitau,pdata,nbindex,nindex)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: pfileid, pitau, nbindex, nindex(nbindex)
  REAL,DIMENSION(:,:),INTENT(IN) :: pdata
  CHARACTER(LEN=*),INTENT(IN) :: pvarname
!-
  LOGICAL :: do_oper, do_write, largebuf
  INTEGER :: varid, io, nbpt_in(1:2), nbpt_out
  REAL, ALLOCATABLE,SAVE :: buff_tmp(:)
  INTEGER,SAVE :: buff_tmp_sz
  CHARACTER(LEN=7) :: tmp_opp
!-
  LOGICAL :: check = .FALSE.
!---------------------------------------------------------------------
!-
! 1.0 Try to catch errors like specifying the wrong file ID.
!     Thanks Marine for showing us what errors users can make !
!-
  IF ( (pfileid < 1).OR.(pfileid > nb_files) ) THEN
    CALL ipslerr (3,"histwrite", &
 &    'Illegal file ID in the histwrite of variable',pvarname,' ')
  ENDIF
!-
! 1.1 Find the id of the variable to be written and the real time
!-
  CALL histvar_seq (pfileid,pvarname,varid)
!-
! 2.0 do nothing for never operation
!-
  tmp_opp = topp(pfileid,varid)
!-
  IF (TRIM(tmp_opp) == "never") THEN
    last_opp_chk(pfileid,varid) = -99
    last_wrt_chk(pfileid,varid) = -99
  ENDIF
!-
! 3.0 We check if we need to do an operation
!-
  IF (last_opp_chk(pfileid,varid) == pitau) THEN
    CALL ipslerr (3,"histwrite", &
 &    'This variable as already been analysed at the present', &
 &    'time step',' ')
  ENDIF
!-
  CALL isittime &
 &  (pitau, date0(pfileid), deltat(pfileid), freq_opp(pfileid,varid), &
 &   last_opp(pfileid,varid), last_opp_chk(pfileid,varid), do_oper)
!-
! 4.0 We check if we need to write the data
!-
  IF (last_wrt_chk(pfileid,varid) == pitau) THEN
    CALL ipslerr (3,"histwrite", &
 &    'This variable as already been written for the present', &
 &    'time step',' ')
  ENDIF
!-
  CALL isittime &
 &  (pitau, date0(pfileid), deltat(pfileid), freq_wrt(pfileid,varid), &
 &   last_wrt(pfileid,varid), last_wrt_chk(pfileid,varid), do_write)
!-
! 5.0 histwrite called
!-
  IF (do_oper.OR.do_write) THEN
!-
!-- 5.1 Get the sizes of the data we will handle
!-
    IF (datasz_in(pfileid,varid,1) <= 0) THEN
!---- There is the risk here that the user has over-sized the array.
!---- But how can we catch this ?
!---- In the worst case we will do impossible operations
!---- on part of the data !
      datasz_in(pfileid,varid,1) = SIZE(pdata, DIM=1)
      datasz_in(pfileid,varid,2) = SIZE(pdata, DIM=2)
      datasz_in(pfileid,varid,3) = -1
    ENDIF
!-
!-- 5.2 The maximum size of the data will give the size of the buffer
!-
    IF (datasz_max(pfileid,varid) <= 0) THEN
      largebuf = .FALSE.
      DO io=1,nbopp(pfileid,varid)
        IF (INDEX(fuchnbout,sopps(pfileid,varid,io)) > 0) THEN
          largebuf = .TRUE.
        ENDIF
      ENDDO
      IF (largebuf) THEN
        datasz_max(pfileid,varid) = &
 &        scsize(pfileid,varid,1) &
 &       *scsize(pfileid,varid,2) &
 &       *scsize(pfileid,varid,3)
      ELSE
        datasz_max(pfileid,varid) = &
 &        datasz_in(pfileid,varid,1) &
 &       *datasz_in(pfileid,varid,2)
      ENDIF
    ENDIF
!-
    IF (.NOT.ALLOCATED(buff_tmp)) THEN
      IF (check) THEN
        WRITE(*,*) &
 &       "histwrite_r2d : allocate buff_tmp for buff_sz = ", &
 &       datasz_max(pfileid,varid)
      ENDIF
      ALLOCATE (buff_tmp(datasz_max(pfileid,varid)))
      buff_tmp_sz = datasz_max(pfileid,varid)
    ELSE IF (datasz_max(pfileid,varid) > buff_tmp_sz) THEN
      IF (check) THEN
        WRITE(*,*) &
 &       "histwrite_r2d : re-allocate buff_tmp for buff_sz = ", &
 &       datasz_max(pfileid,varid)
      ENDIF
      DEALLOCATE (buff_tmp)
      ALLOCATE (buff_tmp(datasz_max(pfileid,varid)))
      buff_tmp_sz = datasz_max(pfileid,varid)
    ENDIF
!-
!-- We have to do the first operation anyway.
!-- Thus we do it here and change the ranke
!-- of the data at the same time. This should speed up things.
!-
    nbpt_in(1:2) = datasz_in(pfileid,varid,1:2)
    nbpt_out = datasz_max(pfileid,varid)
    CALL mathop (sopps(pfileid,varid,1), nbpt_in, pdata, &
 &               missing_val, nbindex, nindex, &
 &               scal(pfileid,varid,1), nbpt_out, buff_tmp)
    CALL histwrite_real (pfileid, varid, pitau, nbpt_out, &
 &            buff_tmp, nbindex, nindex, do_oper, do_write)
  ENDIF
!-
! 6.0 Manage time steps
!-
  IF ((TRIM(tmp_opp) /= "once").AND.(TRIM(tmp_opp) /= "never")) THEN
    last_opp_chk(pfileid,varid) = pitau
    last_wrt_chk(pfileid,varid) = pitau
  ELSE
    last_opp_chk(pfileid,varid) = -99
    last_wrt_chk(pfileid,varid) = -99
  ENDIF
!---------------------------
END SUBROUTINE histwrite_r2d
!===
SUBROUTINE histwrite_r3d (pfileid,pvarname,pitau,pdata,nbindex,nindex)
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: pfileid, pitau, nbindex, nindex(nbindex)
  REAL,DIMENSION(:,:,:),INTENT(IN) :: pdata
  CHARACTER(LEN=*),INTENT(IN) :: pvarname
!-
  LOGICAL :: do_oper, do_write, largebuf
  INTEGER :: varid, io, nbpt_in(1:3), nbpt_out
  REAL, ALLOCATABLE,SAVE :: buff_tmp(:)
  INTEGER,SAVE :: buff_tmp_sz
  CHARACTER(LEN=7) :: tmp_opp
!-
  LOGICAL :: check = .FALSE.
!---------------------------------------------------------------------
!-
! 1.0 Try to catch errors like specifying the wrong file ID.
!     Thanks Marine for showing us what errors users can make !
!-
  IF ( (pfileid < 1).OR.(pfileid > nb_files) ) THEN
    CALL ipslerr (3,"histwrite", &
 &    'Illegal file ID in the histwrite of variable',pvarname,' ')
  ENDIF
!-
! 1.1 Find the id of the variable to be written and the real time
!-
  CALL histvar_seq (pfileid,pvarname,varid)
!-
! 2.0 do nothing for never operation
!-
  tmp_opp = topp(pfileid,varid)
!-
  IF (TRIM(tmp_opp) == "never") THEN
    last_opp_chk(pfileid,varid) = -99
    last_wrt_chk(pfileid,varid) = -99
  ENDIF
!-
! 3.0 We check if we need to do an operation
!-
  IF (last_opp_chk(pfileid,varid) == pitau) THEN
    CALL ipslerr (3,"histwrite", &
 &    'This variable as already been analysed at the present', &
 &    'time step',' ')
  ENDIF
!-
  CALL isittime &
 &  (pitau, date0(pfileid), deltat(pfileid), freq_opp(pfileid,varid), &
 &   last_opp(pfileid,varid), last_opp_chk(pfileid,varid), do_oper)
!-
! 4.0 We check if we need to write the data
!-
  IF (last_wrt_chk(pfileid,varid) == pitau) THEN
    CALL ipslerr (3,"histwrite", &
 &    'This variable as already been written for the present', &
 &    'time step',' ')
  ENDIF
!-
  CALL isittime &
 &  (pitau, date0(pfileid), deltat(pfileid), freq_wrt(pfileid,varid), &
 &   last_wrt(pfileid,varid), last_wrt_chk(pfileid,varid), do_write)
!-
! 5.0 histwrite called
!-
  IF (do_oper.OR.do_write) THEN
!-
!-- 5.1 Get the sizes of the data we will handle
!-
    IF (datasz_in(pfileid,varid,1) <= 0) THEN
!---- There is the risk here that the user has over-sized the array.
!---- But how can we catch this ?
!---- In the worst case we will do impossible operations
!---- on part of the data !
      datasz_in(pfileid,varid,1) = SIZE(pdata, DIM=1)
      datasz_in(pfileid,varid,2) = SIZE(pdata, DIM=2)
      datasz_in(pfileid,varid,3) = SIZE(pdata, DIM=3)
    ENDIF
!-
!-- 5.2 The maximum size of the data will give the size of the buffer
!-
    IF (datasz_max(pfileid,varid) <= 0) THEN
      largebuf = .FALSE.
      DO io =1,nbopp(pfileid,varid)
        IF (INDEX(fuchnbout,sopps(pfileid,varid,io)) > 0) THEN
          largebuf = .TRUE.
        ENDIF
      ENDDO
      IF (largebuf) THEN
        datasz_max(pfileid,varid) = &
 &        scsize(pfileid,varid,1) &
 &       *scsize(pfileid,varid,2) &
 &       *scsize(pfileid,varid,3)
      ELSE
        datasz_max(pfileid,varid) = &
 &        datasz_in(pfileid,varid,1) &
 &       *datasz_in(pfileid,varid,2) &
 &       *datasz_in(pfileid,varid,3)
      ENDIF
    ENDIF
!-
    IF (.NOT.ALLOCATED(buff_tmp)) THEN
      IF (check) THEN
        WRITE(*,*) &
 &       "histwrite_r1d : allocate buff_tmp for buff_sz = ", &
 &       datasz_max(pfileid,varid)
      ENDIF
      ALLOCATE (buff_tmp(datasz_max(pfileid,varid)))
      buff_tmp_sz = datasz_max(pfileid,varid)
    ELSE IF (datasz_max(pfileid,varid) > buff_tmp_sz) THEN
      IF (check) THEN
        WRITE(*,*) &
 &       "histwrite_r1d : re-allocate buff_tmp for buff_sz = ", &
 &       datasz_max(pfileid,varid)
      ENDIF
      DEALLOCATE (buff_tmp)
      ALLOCATE (buff_tmp(datasz_max(pfileid,varid)))
      buff_tmp_sz = datasz_max(pfileid,varid)
    ENDIF
!-
!-- We have to do the first operation anyway.
!-- Thus we do it here and change the ranke
!-- of the data at the same time. This should speed up things.
!-
    nbpt_in(1:3) = datasz_in(pfileid,varid,1:3)
    nbpt_out = datasz_max(pfileid,varid)
    CALL mathop (sopps(pfileid,varid,1), nbpt_in, pdata, &
 &               missing_val, nbindex, nindex, &
 &               scal(pfileid,varid,1), nbpt_out, buff_tmp)
    CALL histwrite_real (pfileid, varid, pitau, nbpt_out, &
 &            buff_tmp, nbindex, nindex, do_oper, do_write)
  ENDIF
!-
! 6.0 Manage time steps
!-
  IF ((TRIM(tmp_opp) /= "once").AND.(TRIM(tmp_opp) /= "never")) THEN
    last_opp_chk(pfileid,varid) = pitau
    last_wrt_chk(pfileid,varid) = pitau
  ELSE
    last_opp_chk(pfileid,varid) = -99
    last_wrt_chk(pfileid,varid) = -99
  ENDIF
!---------------------------
END SUBROUTINE histwrite_r3d
!===
SUBROUTINE histwrite_real &
 & (pfileid,varid,pitau,nbdpt,buff_tmp,nbindex,nindex,do_oper,do_write)
!---------------------------------------------------------------------
!- This subroutine is internal and does the calculations and writing
!- if needed. At a later stage it should be split into an operation
!- and writing subroutines.
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(IN) :: pfileid,pitau,varid, &
 &                      nbindex,nindex(nbindex),nbdpt
  REAL,DIMENSION(:)  :: buff_tmp
  LOGICAL,INTENT(IN) :: do_oper,do_write
!-
  INTEGER :: tsz, ncid, ncvarid
  INTEGER :: i, iret, ipt, itax
  INTEGER :: io, nbin, nbout
  INTEGER,DIMENSION(4) :: corner, edges
  INTEGER :: itime
!-
  REAL :: rtime
  CHARACTER(LEN=7) :: tmp_opp
  REAL :: mindata, maxdata
  REAL,ALLOCATABLE,SAVE :: buff_tmp2(:)
  INTEGER,SAVE          :: buff_tmp2_sz
  REAL,ALLOCATABLE,SAVE :: buffer_used(:)
  INTEGER,SAVE          :: buffer_sz
!-xx  INTEGER :: ji
!-
  LOGICAL :: check = .FALSE.
!---------------------------------------------------------------------
  IF (check) THEN
    WRITE(*,*) "histwrite 0.0 :  VAR : ", name(pfileid,varid)
    WRITE(*,*) "histwrite 0.0 : nbindex, nindex :", &
    &  nbindex,nindex(1:MIN(3,nbindex)),'...',nindex(MAX(1,nbindex-3):nbindex)
  ENDIF
!-
! The sizes which can be encoutered
!-
  tsz = zsize(pfileid,varid,1)*zsize(pfileid,varid,2)*zsize(pfileid,varid,3)
!-
! 1.0 We allocate the memory needed to store the data between write
!     and the temporary space needed for operations.
!     We have to keep precedent buffer if needed
!-
  IF (.NOT. ALLOCATED(buffer)) THEN
    IF (check) WRITE(*,*) "histwrite_real 1.0 allocate buffer ",buff_pos
    ALLOCATE(buffer(buff_pos))
    buffer_sz = buff_pos
    buffer(:)=0.0
  ELSE IF (buffer_sz < buff_pos) THEN
    IF (check) WRITE(*,*) "histwrite_real 1.0.1 re-allocate buffer for ",buff_pos," instead of ",SIZE(buffer)
    IF (SUM(buffer)/=0.0) THEN
      IF (check) WRITE (*,*) ' histwrite : buffer has been used. We have to save it before re-allocating it '
      ALLOCATE (buffer_used(buffer_sz))
      buffer_used(:)=buffer(:)
      DEALLOCATE (buffer)
      ALLOCATE (buffer(buff_pos))
      buffer_sz = buff_pos
      buffer(:SIZE(buffer_used))=buffer_used
      DEALLOCATE (buffer_used)
    ELSE
      IF (check) WRITE (*,*) ' histwrite : buffer has not been used. We have just to re-allocating it '
      DEALLOCATE (buffer)
      ALLOCATE (buffer(buff_pos))
      buffer_sz = buff_pos
      buffer(:)=0.0
    ENDIF
  ENDIF
!-
! The buffers are only deallocated when more space is needed. This
! reduces the umber of allocates but increases memory needs.
!-
  IF (.NOT.ALLOCATED(buff_tmp2)) THEN
    IF (check) THEN
      WRITE(*,*) "histwrite_real 1.1 allocate buff_tmp2 ",SIZE(buff_tmp)
    ENDIF
    ALLOCATE (buff_tmp2(datasz_max(pfileid,varid)))
    buff_tmp2_sz = datasz_max(pfileid,varid)
  ELSE IF ( datasz_max(pfileid,varid) > buff_tmp2_sz) THEN
    IF (check) THEN
      WRITE(*,*) "histwrite_real 1.2 re-allocate buff_tmp2 : ", &
     & SIZE(buff_tmp)," instead of ",SIZE(buff_tmp2)
    ENDIF
    DEALLOCATE (buff_tmp2)
    ALLOCATE (buff_tmp2(datasz_max(pfileid,varid)))
    buff_tmp2_sz = datasz_max(pfileid,varid)
  ENDIF
!-
  rtime = pitau * deltat(pfileid)
  tmp_opp = topp(pfileid,varid)
!-
! 3.0 Do the operations or transfer the slab of data into buff_tmp
!-
  IF (check) WRITE(*,*) "histwrite: 3.0", pfileid
!-
! 3.1 DO the Operations only if needed
!-
  IF ( do_oper ) THEN
    i = pfileid
    nbout = nbdpt
!-
!-- 3.4 We continue the sequence of operations
!--     we started in the interface routine
!-
    DO io = 2, nbopp(i,varid),2
      nbin = nbout
      nbout = datasz_max(i,varid)
      CALL mathop(sopps(i,varid,io),nbin,buff_tmp,missing_val, &
 &      nbindex,nindex,scal(i,varid,io),nbout,buff_tmp2)
      IF (check) THEN
        WRITE(*,*) &
 &       "histwrite: 3.4a nbout : ",nbin,nbout,sopps(i,varid,io)
      ENDIF
!-
      nbin = nbout
      nbout = datasz_max(i,varid)
      CALL mathop(sopps(i,varid,io+1),nbin,buff_tmp2,missing_val, &
 &      nbindex,nindex,scal(i,varid,io+1),nbout,buff_tmp)
      IF (check) THEN
        WRITE(*,*) &
 &       "histwrite: 3.4b nbout : ",nbin,nbout,sopps(i,varid,io+1)
      ENDIF
    ENDDO
!-
!   3.5 Zoom into the data
!-
    IF (check) THEN
      WRITE(*,*) &
 &     "histwrite: 3.5 size(buff_tmp) : ",SIZE(buff_tmp)
      WRITE(*,*) &
 &     "histwrite: 3.5 slab in X :",zorig(i,varid,1),zsize(i,varid,1)
      WRITE(*,*) &
 &     "histwrite: 3.5 slab in Y :",zorig(i,varid,2),zsize(i,varid,2)
      WRITE(*,*) &
 &     "histwrite: 3.5 slab in Z :",zorig(i,varid,3),zsize(i,varid,3)
      WRITE(*,*) &
 &     "histwrite: 3.5 slab of input:", &
 &     scsize(i,varid,1),scsize(i,varid,2),scsize(i,varid,3)
    ENDIF
    CALL trans_buff &
 &      (zorig(i,varid,1),zsize(i,varid,1), &
 &       zorig(i,varid,2),zsize(i,varid,2), &
 &       zorig(i,varid,3),zsize(i,varid,3), &
 &       scsize(i,varid,1),scsize(i,varid,2),scsize(i,varid,3), &
 &       buff_tmp, buff_tmp2_sz,buff_tmp2)
!-
!-- 4.0 Get the min and max of the field (buff_tmp)
!-
    IF (check) WRITE(*,*) "histwrite: 4.0 buff_tmp",pfileid,varid, &
 &    TRIM(tmp_opp),'----',LEN_TRIM(tmp_opp),nbindex
!-
    mindata =  ABS(missing_val)
    maxdata = -ABS(missing_val)
!-xx DO ji=1,tsz
!-xx   IF (buff_tmp2(ji)/= missing_val) THEN
!-xx     mindata = MIN(mindata, buff_tmp2(ji))
!-xx     maxdata = MAX(maxdata, buff_tmp2(ji))
!-xx   ENDIF
!-xx ENDDO
!-??
!-??   mindata = MINVAL(buff_tmp2(1:tsz),
!-?? &                  MASK = buff_tmp2(1:tsz) /= missing_val)
!-??   maxdata = MAXVAL(buff_tmp2(1:tsz),
!-?? &                  MASK = buff_tmp2(1:tsz) /= missing_val)
!-??
    minmax(pfileid,varid,1) = mindata
    minmax(pfileid,varid,2) = maxdata
!-
!-- 5.0 Do the operations if needed. In the case of instantaneous
!--     output we do not transfer to the buffer.
!-
    IF (check) WRITE(*,*) "histwrite: 5.0", pfileid, "tsz :", tsz
!-
    ipt = point(pfileid,varid)
!-
!   WRITE(*,*) 'OPE ipt, buffer :', pvarname, ipt, varid
!-
    IF (     (TRIM(tmp_opp) /= "inst") &
   &    .AND.(TRIM(tmp_opp) /= "once") ) THEN
      CALL moycum(tmp_opp,tsz,buffer(ipt:), &
     &       buff_tmp2,nb_opp(pfileid,varid))
    ENDIF
!-
    last_opp(pfileid,varid) = pitau
    nb_opp(pfileid,varid) = nb_opp(pfileid,varid)+1
!-
   ENDIF
!-
! 6.0 Write to file if needed
!-
  IF (check) WRITE(*,*) "histwrite: 6.0", pfileid
!-
  IF ( do_write ) THEN
!-
    ncvarid = ncvar_ids(pfileid,varid)
    ncid = ncdf_ids(pfileid)
!-
!-- 6.1 Do the operations that are needed before writting
!-
    IF (check) WRITE(*,*) "histwrite: 6.1", pfileid
!-
    IF (     (TRIM(tmp_opp) /= "inst") &
   &    .AND.(TRIM(tmp_opp) /= "once") ) THEN
      rtime = (rtime+last_wrt(pfileid,varid)*deltat(pfileid))/2.0
    ENDIF
!-
!-- 6.2 Add a value to the time axis of this variable if needed
!-
    IF (     (TRIM(tmp_opp) /= "l_max") &
   &    .AND.(TRIM(tmp_opp) /= "l_min") &
   &    .AND.(TRIM(tmp_opp) /= "once") ) THEN
!-
      IF (check) WRITE(*,*) "histwrite: 6.2", pfileid
!-
      itax = var_axid(pfileid, varid)
      itime = nb_wrt(pfileid, varid)+1
!-
      IF (tax_last(pfileid, itax) < itime) THEN
        iret = NF90_PUT_VAR (ncid,tdimid(pfileid,itax),(/ rtime /), &
 &                            start=(/ itime /),count=(/ 1 /))
        tax_last(pfileid, itax) = itime
      ENDIF
    ELSE
      itime=1
    ENDIF
!-
!-- 6.3 Write the data. Only in the case of instantaneous output
!       we do not write the buffer.
!-
    IF (check) THEN
      WRITE(*,*) "histwrite: 6.3",pfileid,ncid,ncvarid,varid,itime
    ENDIF
!-
    IF (scsize(pfileid,varid,3) == 1) THEN
      IF (regular(pfileid)) THEN
        corner(1:4) = (/ 1, 1, itime, 0 /)
        edges(1:4) = (/ zsize(pfileid,varid,1), &
 &                      zsize(pfileid,varid,2), &
 &                       1, 0 /)
      ELSE
        corner(1:4) = (/ 1, itime, 0, 0 /)
        edges(1:4) = (/ zsize(pfileid,varid,1), 1, 0, 0 /)
      ENDIF
    ELSE
      IF ( regular(pfileid) ) THEN
        corner(1:4) = (/ 1, 1, 1, itime /)
        edges(1:4) = (/ zsize(pfileid,varid,1), &
 &                      zsize(pfileid,varid,2), &
 &                      zsize(pfileid,varid,3), 1 /)
      ELSE
        corner(1:4) = (/ 1, 1, itime, 0 /)
        edges(1:4) = (/ zsize(pfileid,varid,1), &
 &                      zsize(pfileid,varid,3), 1, 0 /)
      ENDIF
    ENDIF
!-
    ipt = point(pfileid,varid)
!-
    IF (     (TRIM(tmp_opp) /= "inst") &
 &      .AND.(TRIM(tmp_opp) /= "once") ) THEN
      iret = NF90_PUT_VAR (ncid,ncvarid,buffer(ipt:), &
 &                       start=corner(1:4),count=edges(1:4))
    ELSE
      iret = NF90_PUT_VAR (ncid,ncvarid,buff_tmp2, &
 &                       start=corner(1:4),count=edges(1:4))
    ENDIF
!-
    last_wrt(pfileid,varid) = pitau
    nb_wrt(pfileid,varid) = nb_wrt(pfileid,varid)+1
    nb_opp(pfileid,varid) = 0
!---
!   After the write the file can be synchronized so that no data is
!   lost in case of a crash. This feature gives up on the benefits of
!   buffering and should only be used in debuging mode. A flag is
!   needed here to switch to this mode.
!---
!   iret = NF90_SYNC (ncid)
!-
  ENDIF
!----------------------------
END SUBROUTINE histwrite_real
!===
SUBROUTINE histvar_seq (pfid,pvarname,pvid)
!---------------------------------------------------------------------
!- This subroutine optimized the search for the variable in the table.
!- In a first phase it will learn the succession of the variables
!- called and then it will use the table to guess what comes next.
!- It is the best solution to avoid lengthy searches through array
!- vectors.
!-
!- ARGUMENTS :
!-
!- pfid  : id of the file on which we work
!- pvarname : The name of the variable we are looking for
!- pvid     : The var id we found
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  INTEGER,INTENT(in)  :: pfid
  CHARACTER(LEN=*),INTENT(IN) :: pvarname
  INTEGER,INTENT(out) :: pvid
!-
  LOGICAL,SAVE :: learning(nb_files_max)=.TRUE.
  INTEGER,SAVE :: overlap(nb_files_max) = -1
  INTEGER,SAVE :: varseq(nb_files_max, nb_var_max*3)
  INTEGER,SAVE :: varseq_len(nb_files_max) = 0
  INTEGER,SAVE :: varseq_pos(nb_files_max)
  INTEGER,SAVE :: varseq_err(nb_files_max) = 0
  INTEGER      :: ib, nb, sp, nx, pos
  CHARACTER(LEN=20),DIMENSION(nb_var_max) :: tab_str20
  CHARACTER(LEN=20) :: str20
  CHARACTER(LEN=70) :: str70
  INTEGER      :: tab_str20_length(nb_var_max)
!-
  LOGICAL :: check = .FALSE.
!---------------------------------------------------------------------
  nb = nb_var(pfid)
!-
  IF (check) THEN
    WRITE(*,*) 'histvar_seq, start of the subroutine :',learning(pfid)
  ENDIF
!-
  IF (learning(pfid)) THEN
!-
!-- 1.0 We compute the length over which we are going
!--     to check the overlap
!-
    IF (overlap(pfid) <= 0) THEN
      IF (nb_var(pfid) > 6) THEN
        overlap(pfid) = nb_var(pfid)/3*2
      ELSE
        overlap(pfid) = nb_var(pfid)
      ENDIF
    ENDIF
!-
!-- 1.1 Find the position of this string
!-
    str20 = pvarname
    tab_str20(1:nb) = name(pfid,1:nb)
    tab_str20_length(1:nb) = name_length(pfid,1:nb)
!-
    CALL find_str (nb, tab_str20, tab_str20_length, str20, pos)
!-
    IF (pos > 0) THEN
      pvid = pos
    ELSE
      CALL ipslerr (3,"histvar_seq", &
 &      'The name of the variable you gave has not been declared', &
 &      'You should use subroutine histdef for declaring variable', &
 &      TRIM(str20))
    ENDIF
!-
!-- 1.2 If we have not given up we store the position
!--     in the sequence of calls
!-
    IF ( varseq_err(pfid) .GE. 0 ) THEN
      sp = varseq_len(pfid)+1
      IF (sp <= nb_var_max*3) THEN
        varseq(pfid,sp) = pvid
        varseq_len(pfid) = sp
      ELSE
        CALL ipslerr (2,"histvar_seq",&
 &       'The learning process has failed and we give up. '// &
 &       'Either you sequence is',&
 &       'too complex or I am too dumb. '// &
 &       'This will only affect the efficiency',&
 &       'of your code. Thus if you wish to save time'// &
 &       ' contact the IOIPSL team. ')
        WRITE(*,*) 'The sequence we have found up to now :'
        WRITE(*,*) varseq(pfid,1:sp-1)
        varseq_err(pfid) = -1
      ENDIF
!-
!---- 1.3 Check if we have found the right overlap
!-
      IF (varseq_len(pfid) .GE. overlap(pfid)*2) THEN
!-
!------ We skip a few variables if needed as they could come
!------ from the initialisation of the model.
!-
        DO ib = 0, sp-overlap(pfid)*2
          IF ( learning(pfid) .AND.&
            & SUM(ABS(varseq(pfid,ib+1:ib+overlap(pfid)) -&
            & varseq(pfid,sp-overlap(pfid)+1:sp))) == 0 ) THEN
            learning(pfid) = .FALSE.
            varseq_len(pfid) = sp-overlap(pfid)-ib
            varseq_pos(pfid) = overlap(pfid)+ib
            varseq(pfid,1:varseq_len(pfid)) = &
 &            varseq(pfid,ib+1:ib+varseq_len(pfid))
          ENDIF
        ENDDO
      ENDIF
    ENDIF
  ELSE
!-
!-- 2.0 Now we know how the calls to histwrite are sequenced
!--     and we can get a guess at the var ID
!-
    nx = varseq_pos(pfid)+1
    IF (nx > varseq_len(pfid)) nx = 1
!-
    pvid = varseq(pfid, nx)
!-
    IF (    (INDEX(name(pfid,pvid),pvarname) <= 0)         &
   &    .OR.(name_length(pfid,pvid) /= len_trim(pvarname)) ) THEN
      str20 = pvarname
      tab_str20(1:nb) = name(pfid,1:nb)
      tab_str20_length(1:nb) = name_length(pfid,1:nb)
      CALL find_str (nb,tab_str20,tab_str20_length,str20,pos)
      IF (pos > 0) THEN
        pvid = pos
      ELSE
        CALL ipslerr (3,"histvar_seq", &
 &  'The name of the variable you gave has not been declared',&
 &  'You should use subroutine histdef for declaring variable',str20)
      ENDIF
      varseq_err(pfid) = varseq_err(pfid)+1
    ELSE
!-
!---- We only keep the new position if we have found the variable
!---- this way. This way an out of sequence call to histwrite does
!---- not defeat the process.
!-
      varseq_pos(pfid) = nx
    ENDIF
!-
    IF (varseq_err(pfid) .GE. 10) THEN
      WRITE(str70,'("for file ",I3)') pfid
      CALL ipslerr (2,"histvar_seq", &
 &  'There were 10 errors in the learned sequence of variables',&
 &  str70,'This looks like a bug, please report it.')
         varseq_err(pfid) = 0
    ENDIF
  ENDIF
!-
  IF (check) THEN
    WRITE(*,*) &
 &   'histvar_seq, end of the subroutine :',TRIM(pvarname),pvid
  ENDIF
!-------------------------
END SUBROUTINE histvar_seq
!===
SUBROUTINE histsync (file)
!---------------------------------------------------------------------
!- This subroutine will synchronise all
!- (or one if defined) opened files.
!-
!- VERSION
!-
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
! file  : optional argument for fileid
  INTEGER,INTENT(in),OPTIONAL :: file
!-
  INTEGER :: ifile,ncid,iret
!-
  LOGICAL :: file_exists
  LOGICAL :: check = .FALSE.
!---------------------------------------------------------------------
  IF (check) WRITE(*,*) 'Entering loop on files :', nb_files
!-
! 1.The loop on files to synchronise
!-
  DO ifile = 1,nb_files
!-
    IF (PRESENT(file)) THEN
      file_exists = (ifile == file)
    ELSE
      file_exists = .TRUE.
    ENDIF
!-
    IF ( file_exists ) THEN
      IF (check) THEN
        WRITE(*,*) 'Synchronising specified file number :', file
      ENDIF
      ncid = ncdf_ids(ifile)
      iret = NF90_SYNC (ncid)
    ENDIF
!-
  ENDDO
!----------------------
END SUBROUTINE histsync
!===
SUBROUTINE histclo (fid)
!---------------------------------------------------------------------
!- This subroutine will close all (or one if defined) opened files
!-
!- VERSION
!-
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
! fid  : optional argument for fileid
  INTEGER,INTENT(in),OPTIONAL :: fid
!-
  INTEGER :: ifile,ncid,iret,iv,ncvarid
  INTEGER :: start_loop,end_loop
  CHARACTER(LEN=70) :: str70
!-
  LOGICAL :: check=.FALSE.
!---------------------------------------------------------------------
  IF (check) WRITE(*,*) 'Entering loop on files :', nb_files
!-
  IF (PRESENT(fid)) THEN
    start_loop = fid
    end_loop = fid
  ELSE
    start_loop = 1
    end_loop = nb_files
  ENDIF
!-
  DO ifile=start_loop,end_loop
    IF (check) WRITE(*,*) 'Closing specified file number :', ifile
    ncid = ncdf_ids(ifile)
    iret = NF90_REDEF (ncid)
!-
!-- 1. The loop on the number of variables to add
!-     some final information
!-
    IF ( check ) WRITE(*,*) 'Entering loop on vars :', nb_var(ifile)
    DO iv = 1,nb_var(ifile)
      ncvarid = ncvar_ids(ifile,iv)
      IF (check) THEN
        WRITE(*,*) 'min value for file :',ifile,' var n. :',iv, &
       &           ' is : ',minmax(ifile,iv,1)
        WRITE(*,*) 'max value for file :',ifile,' var n. :',iv, &
       &           ' is : ',minmax(ifile,iv,2)
      ENDIF
!-
!---- 1.1 Put the min and max values on the file
!-
      iret = NF90_PUT_ATT (ncid,ncvarid,'valid_min', &
 &                         REAL(minmax(ifile,iv,1),KIND=4))
      iret = NF90_PUT_ATT (ncid,ncvarid,'valid_max', &
 &                         REAL(minmax(ifile,iv,2),KIND=4))
    ENDDO
!-
!-- 2.0 We list the names of the other files
!--     in the associated_file attribute
!-
    IF (nb_files > 1 ) THEN
      iret = NF90_PUT_ATT (ncid,NF90_GLOBAL,'associate_file', &
 &                         TRIM(assc_file))
    ENDIF
    IF ( check ) WRITE(*,*) 'close file :', ncid
    iret = NF90_CLOSE (ncid)
    IF (iret /= NF90_NOERR) THEN
      WRITE(str70,'("This file has been already closed :",I3)') ifile
      CALL ipslerr (2,'histclo',str70,'',' ')
    ENDIF
  ENDDO
!---------------------
END SUBROUTINE histclo
!===
SUBROUTINE ioconf_modname (str)
!---------------------------------------------------------------------
!- This subroutine allows to configure the name
!- of the model written into the file
!---------------------------------------------------------------------
  IMPLICIT NONE
!-
  CHARACTER(LEN=*),INTENT(IN) :: str
!---------------------------------------------------------------------
  IF (.NOT.lock_modname) THEN
    model_name = str(1:MIN(LEN_TRIM(str),80))
    lock_modname = .TRUE.
  ELSE
    CALL ipslerr (2,"ioconf_modname", &
   &  'The model name can only be changed once and only', &
   &  'before it is used. It is now set to :',model_name)
  ENDIF
!----------------------------
END SUBROUTINE ioconf_modname
!-
!===
!-
!-----------------
END MODULE histcom
