#if defined key_agrif
      SUBROUTINE Agrif_InitWorkspace
!
!     Modules used:
!
      Use par_oce
      Use dom_oce
      USE Agrif_Util
!
!     Declarations:
!      
      IMPLICIT NONE
!
!     Variables      
!

!
!     Begin
!
      if ( .NOT. Agrif_Root() ) then
         jpiglo = nbcellsx + 2 + 2*nbghostcells
         jpjglo = nbcellsy + 2 + 2*nbghostcells
         jpi = ( jpiglo-2*jpreci + (jpni-1+0) ) / jpni + 2*jpreci
         jpj = ( jpjglo-2*jprecj + (jpnj-1+0) ) / jpnj + 2*jprecj
         jpim1 = jpi-1
         jpjm1 = jpj-1
         jpkm1 = jpk-1                                        
         jpij  = jpi*jpj
         jpidta = jpiglo
         jpjdta = jpjglo
         jpizoom = 1
         jpjzoom = 1
         nperio = 0
         jperio = 0
      endif


      Return
      End Subroutine Agrif_InitWorkspace

!
      SUBROUTINE Agrif_InitValues
!     ------------------------------------------------------------------
!     You should declare the variable which has to be interpolated here
!     -----------------------------------------------------------------
!
!     Modules used:
!  
      USE Agrif_Util
      USE oce
      USE dom_oce
      USE opa
#if   defined key_tradmp   ||   defined key_esopa
      USE tradmp
#endif
      USE sol_oce
      USE in_out_manager
#if defined key_ice_lim
      USE ice_oce
#endif
#if defined key_agrif
     USE agrif_opa_update
     USE agrif_opa_interp
     USE agrif_opa_sponge
#endif
!
!     Declarations:
!      
      Implicit none
!
!     Variables
!
      REAL(wp) tabtemp(jpi,jpj,jpk)
! 
      LOGICAL check_namelist
!
!
!     Begin
!
#if defined key_orca_r025 || defined key_orca_r05 || defined key_orca_r2 || defined key_orca_r4
      jp_cfg = -1  ! set special value for jp_cfg on fine grids
      cp_cfg = "default"
#endif

      Call opa_init  ! Initializations of each fine grid
!
!     Specific fine grid Initializations
!
#if defined key_tradmp || defined key_esopa
! no tracer damping on fine grids
      lk_tradmp = .FALSE.
#endif
!      
!     Declaration of the type of variable which have to be interpolated
!
      Call Agrif_Set_type(un,(/1,2,0/),(/2,3,0/))
      Call Agrif_Set_type(vn,(/2,1,0/),(/3,2,0/))

      Call Agrif_Set_type(ua,(/1,2,0/),(/2,3,0/))
      Call Agrif_Set_type(va,(/2,1,0/),(/3,2,0/))

      Call Agrif_Set_type(e1u,(/1,2/),(/2,3/))
      Call Agrif_Set_type(e2v,(/2,1/),(/3,2/))
            
      Call Agrif_Set_type(tn,(/2,2,0/),(/3,3,0/))
      Call Agrif_Set_type(sn,(/2,2,0/),(/3,3,0/)) 

      Call Agrif_Set_type(tb,(/2,2,0/),(/3,3,0/))
      Call Agrif_Set_type(sb,(/2,2,0/),(/3,3,0/)) 
      
      Call Agrif_Set_type(ta,(/2,2,0/),(/3,3,0/))
      Call Agrif_Set_type(sa,(/2,2,0/),(/3,3,0/))       
            
      Call Agrif_Set_type(sshn,(/2,2/),(/3,3/))
      Call Agrif_Set_type(gcb,(/2,2/),(/3,3/))

!
!     Space directions for each variables
!
      Call Agrif_Set_raf(un,(/'x','y','N'/))
      Call Agrif_Set_raf(vn,(/'x','y','N'/))
      
      Call Agrif_Set_raf(ua,(/'x','y','N'/))
      Call Agrif_Set_raf(va,(/'x','y','N'/))

      Call Agrif_Set_raf(e1u,(/'x','y'/))
      Call Agrif_Set_raf(e2v,(/'x','y'/))

      Call Agrif_Set_raf(tn,(/'x','y','N'/))
      Call Agrif_Set_raf(sn,(/'x','y','N'/))
      
      Call Agrif_Set_raf(tb,(/'x','y','N'/))
      Call Agrif_Set_raf(sb,(/'x','y','N'/))
      
      Call Agrif_Set_raf(ta,(/'x','y','N'/))
      Call Agrif_Set_raf(sa,(/'x','y','N'/))      
            
      Call Agrif_Set_raf(sshn,(/'x','y'/))
      Call Agrif_Set_raf(gcb,(/'x','y'/))

!
!     type of interpolation

      Call Agrif_Set_bcinterp(tn,interp=AGRIF_linear)
      Call Agrif_Set_bcinterp(sn,interp=AGRIF_linear)
      
      Call Agrif_Set_bcinterp(ta,interp=AGRIF_linear)
      Call Agrif_Set_bcinterp(sa,interp=AGRIF_linear)
               
      Call Agrif_Set_bcinterp(un,interp1=Agrif_linear,interp2=AGRIF_ppm)
      Call Agrif_Set_bcinterp(vn,interp1=AGRIF_ppm,interp2=Agrif_linear)

      Call Agrif_Set_bcinterp(ua,interp1=Agrif_linear,interp2=AGRIF_ppm)
      Call Agrif_Set_bcinterp(va,interp1=AGRIF_ppm,interp2=Agrif_linear)
      
      Call Agrif_Set_bcinterp(e1u,interp1=Agrif_linear,interp2=AGRIF_ppm)
      Call Agrif_Set_bcinterp(e2v,interp1=AGRIF_ppm,interp2=Agrif_linear)

!
!     Location of interpolation
!
      Call Agrif_Set_bc(un,(/0,1/))
      Call Agrif_Set_bc(vn,(/0,1/))
      
      Call Agrif_Set_bc(e1u,(/0,0/))
      Call Agrif_Set_bc(e2v,(/0,0/))

      Call Agrif_Set_bc(tn,(/0,1/))
      Call Agrif_Set_bc(sn,(/0,1/))

      Call Agrif_Set_bc(ta,(/-3*Agrif_irhox(),0/))
      Call Agrif_Set_bc(sa,(/-3*Agrif_irhox(),0/))

      Call Agrif_Set_bc(ua,(/-2*Agrif_irhox(),0/))
      Call Agrif_Set_bc(va,(/-2*Agrif_irhox(),0/))

!     Update type
      
      Call Agrif_Set_Updatetype(tn, update = AGRIF_Update_Average)
      Call Agrif_Set_Updatetype(sn, update = AGRIF_Update_Average)
      
      Call Agrif_Set_Updatetype(tb, update = AGRIF_Update_Average)
      Call Agrif_Set_Updatetype(sb, update = AGRIF_Update_Average)

      Call Agrif_Set_Updatetype(sshn, update = AGRIF_Update_Average)
      Call Agrif_Set_Updatetype(gcb,update = AGRIF_Update_Average)

      Call Agrif_Set_Updatetype(un,update1 = Agrif_Update_Copy, update2 = Agrif_Update_Average)
      Call Agrif_Set_Updatetype(vn,update1 = Agrif_Update_Average, update2 = Agrif_Update_Copy)

      Call Agrif_Set_Updatetype(e1u,update1 = Agrif_Update_Copy, update2=Agrif_Update_Average)
      Call Agrif_Set_Updatetype(e2v,update1 = Agrif_Update_Average, update2=Agrif_Update_Copy)

! First interpolations of potentially non zero fields

       Agrif_SpecialValue=0.
       Agrif_UseSpecialValue = .TRUE.
       Call Agrif_Bc_variable(tabtemp,tn,calledweight=1.)
       Call Agrif_Bc_variable(tabtemp,sn,calledweight=1.)
       Call Agrif_Bc_variable(tabtemp,un,calledweight=1.,procname=interpu)
       Call Agrif_Bc_variable(tabtemp,vn,calledweight=1.,procname=interpv)

       Call Agrif_Bc_variable(tabtemp,ta,calledweight=1.,procname=interptn)
       Call Agrif_Bc_variable(tabtemp,sa,calledweight=1.,procname=interpsn)

       Call Agrif_Bc_variable(tabtemp,ua,calledweight=1.,procname=interpun)
       Call Agrif_Bc_variable(tabtemp,va,calledweight=1.,procname=interpvn)
       Agrif_UseSpecialValue = .FALSE.
!

!
      check_namelist = .true.
!      
      IF( check_namelist ) then     
!
! check time steps           
!
       If( nint(Agrif_Rhot()) * nint(rdt) .ne. Agrif_Parent(rdt) ) then
              Write(*,*) 'incompatible time step between grids'
              Write(*,*) 'parent grid value : ',Agrif_Parent(rdt)
              Write(*,*) 'child  grid value : ',nint(rdt)
              Write(*,*) 'value on parent grid should be : ',rdt*Agrif_Rhot()
              stop
       Endif
           
       If( Agrif_IRhot() * (Agrif_Parent(nitend)- &
       Agrif_Parent(nit000)+1) .ne. (nitend-nit000+1) ) then
            Write(*,*) 'incompatible run length between grids'
            Write(*,*) 'parent grid value : ',(Agrif_Parent(nitend)- &
            Agrif_Parent(nit000)+1),' time step'
            Write(*,*) 'child  grid value : ', &
            (nitend-nit000+1),' time step'
            Write(*,*) 'value on child grid should be : ', &
            Agrif_IRhot() * (Agrif_Parent(nitend)- &
            Agrif_Parent(nit000)+1)
           stop
       Endif           
!
!
#if defined key_partial_steps                                  
!
! check parameters for partial steps 
!
       If( Agrif_Parent(e3zps_min) .ne. e3zps_min ) then
            Write(*,*) 'incompatible e3zps_min between grids'
            Write(*,*) 'parent grid :',Agrif_Parent(e3zps_min)
            Write(*,*) 'child grid  :',e3zps_min
            Write(*,*) 'those values should be identical'
            stop
       Endif          
!          
       If( Agrif_Parent(e3zps_rat) .ne. e3zps_rat ) then
            Write(*,*) 'incompatible e3zps_rat between grids'
            Write(*,*) 'parent grid :',Agrif_Parent(e3zps_rat)
            Write(*,*) 'child grid  :',e3zps_rat
            Write(*,*) 'those values should be identical'                  
            stop
       Endif                  
#endif        
!            
      ENDIF
!
!

      Call Agrif_Update_tra(0)
      Call Agrif_Update_dyn(0)
      
      nbcline = 0

      Return
      End Subroutine Agrif_InitValues
!
      SUBROUTINE Agrif_detect(g,sizex)
!
!     Modules used:
!  
      Use Agrif_Types
!
!
!     Declarations:
!      
!
!     Variables      
!
      Integer, Dimension(2) :: sizex
      Integer, Dimension(sizex(1),sizex(2))   :: g 
!
!     Begin
!
!

!
      Return
      End Subroutine Agrif_detect
      
#if defined key_mpp_mpi
!
!     **************************************************************************
!!!   Subroutine Agrif_InvLoc
!     **************************************************************************
!
      Subroutine Agrif_InvLoc(indloc,nprocloc,i,indglob)
!  
!     Description:
!
      USE dom_oce

!     Declarations:
!  
!!      Implicit none
!
      Integer :: indglob,indloc,nprocloc,i
!
!
      SELECT CASE(i)

      CASE(1)
        indglob = indloc + nimppt(nprocloc+1) - 1

      CASE(2)
        indglob = indloc + njmppt(nprocloc+1) - 1 

      CASE(3)
        indglob = indloc

      CASE(4)
        indglob = indloc

      END SELECT
!
!
      End Subroutine Agrif_InvLoc
#endif

             
#else
      subroutine Subcalledbyagrif
         write(*,*) 'Impossible to bet here'
      end subroutine Subcalledbyagrif
#endif
