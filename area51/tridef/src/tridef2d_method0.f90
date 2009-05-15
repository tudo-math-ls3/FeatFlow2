!##############################################################################
!# ****************************************************************************
!# <name> tridef2d_method0 </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstration program how to deform a grid
!# 
!# </purpose>
!##############################################################################

module tridef2d_method0

  use fsystem
  use genoutput
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use ucd
  use pprocerror
  use genoutput
  use griddeform
    
  use tridef2d_callback
  
  implicit none

contains

  ! ***************************************************************************

!<subroutine>

  subroutine tridef2d_simple
  
!<description>
  ! In this program we are only concerned with a
  ! grid deformation problem. We prescribe a monitor function
  ! and adjust the grid accordingly.
  !
  ! 1.) Read in parametrisation
  ! 2.) Read in triangulation
  ! 3.) Create a discretisation
  ! 4.) Init the Deformation structures
  ! 5.) Start the deformation procedure
  ! 7.) Write the original and the deformed to GMV files
  ! 8.) Release all variables, finish
!</description>

!</subroutine>

    ! Definitions of variables.
    !
    ! We need a couple of variables for this problem. Let's see...
    !
    ! An object for saving the domain:
    type(t_boundary) :: rboundary
    
    type(t_griddefWork) :: rgriddefWork
    
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation.
    ! This contains also information about trial/test functions,...
    type(t_blockDiscretisation) :: rdiscretisation
    
    ! In this vector we will save the area of the deformed grid
    type(t_vectorScalar) :: rvectorAreaQ0 

    ! NLMAX receives the level where we want to solve.
    integer :: NLMAX,NLMIN
    
    ! dummy
    integer :: idummy,iloop    
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    
    ! To these arrays we write the monitorfunction,
    ! elementArea, rhs, etc... in order to
    ! check these values in gmv
    real(dp), dimension(:), pointer :: p_Ddata
    real(dp), dimension(:), pointer :: p_Ddata1
    real(dp), dimension(:), pointer :: p_Ddata2        
    real(dp), dimension(:), pointer :: p_Ddata3 
    real(dp), dimension(:), pointer :: p_Ddata4 
    real(dp), dimension(:), pointer :: p_Ddata5 
    real(dp), dimension(:), pointer :: p_Ddata6 
    real(dp), dimension(:), pointer :: p_Ddata7 
    type(t_vectorBlock) :: rvectorAreaBlockQ0
    type(t_vectorBlock) :: rvectorAreaBlockQ1                 
    
    ! grid deform structure
    type(t_griddefInfo) :: rgriddefInfo

    ! Ok, let's start. 
    !
    ! Here we set the level of the
    ! grid we want to deform
    NLMAX = 6
    NLMIN = 6
    
    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    call boundary_read_prm(rboundary, './pre/QUAD.prm')
        
    ! Now read in the basic triangulation.
    call tria_readTriFile2D (rtriangulation, './pre/QUAD.tri', rboundary)
     
    ! Refine it.
    call tria_quickRefine2LevelOrdering (NLMAX-1,rtriangulation,rboundary)
    
    ! And create information about adjacencies and everything one needs from
    ! a triangulation.
    call tria_initStandardMeshFromRaw (rtriangulation,rboundary)

    ! griddeformation setup
    call griddef_deformationInit(rgriddefInfo,rtriangulation,NLMIN,NLMAX,rboundary,1)                 
    
    do iloop=NLMIN,NLMAX
      call griddef_buildHGrid(rgriddefInfo,rtriangulation,iloop)
    end do
    
    
    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. In this simple problem, we only have one block.
    call spdiscr_initBlockDiscr (rdiscretisation,1,&
                                   rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation,&
                                   rboundary)

    call spdiscr_initBlockDiscr (rdiscretisation,1,&
                                   rtriangulation,&
                                   rboundary)

    
    ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
    ! structures for every component of the solution vector.
    ! Initialise the first element of the list to specify the element
    ! and cubature rule for this solution component:
    call spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
                                   EL_E011,CUB_G2X2,&
                                   rtriangulation,&
                                   rboundary)
                                   
    
    ! Deform
    call griddef_performDeformation(rgriddefInfo, rgriddefWork,idummy,&
                                    .TRUE., .FALSE., .FALSE., &
                                    .FALSE., NLMAX, 0, 0,&
                                    tridef2d_monitorfct)    
                 
    ! That's it, rvectorBlock now contains our solution. We can now
    ! start the postprocessing. 
    ! Start UCD export to GMV file:

    call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       'gmv/u2d_0_simple1.gmv')
    
    call lsyssc_getbase_double (rgriddefWork%rvecGradBlock%RvectorBlock(1),p_Ddata)
    call ucd_addVariableVertexBased (rexport,'Ux',UCD_VAR_STANDARD, p_Ddata)
    
    call lsyssc_getbase_double (rgriddefWork%rvecGradBlock%RvectorBlock(2),p_Ddata1)
    call ucd_addVariableVertexBased (rexport,'Uy',UCD_VAR_STANDARD, p_Ddata1)
    
    call lsyssc_getbase_double (rgriddefWork%rSolBlock%RvectorBlock(1),p_Ddata2)
    call ucd_addVariableVertexBased (rexport,'Sol',UCD_VAR_STANDARD, p_Ddata2)
    
    call lsyssc_getbase_double (rgriddefWork%rrhsBlock%RvectorBlock(1),p_Ddata3)
    call ucd_addVariableVertexBased (rexport,'rhs',UCD_VAR_STANDARD, p_Ddata3) 

    call lsyssc_getbase_double (rgriddefWork%rvectorAreaBlockQ1%RvectorBlock(1),p_Ddata4)
    call ucd_addVariableVertexBased (rexport,'area',UCD_VAR_STANDARD, p_Ddata4) 

    call lsyssc_getbase_double (rgriddefWork%rvectorMonFuncQ1%RvectorBlock(1),p_Ddata5)
    call ucd_addVariableVertexBased (rexport,'Mon',UCD_VAR_STANDARD, p_Ddata5) 
    
    ! write the file to disc, that's it.
    call ucd_write (rexport)
    call ucd_release (rexport)

    call griddef_getAreaDeformed(rgriddefInfo,&
             rvectorAreaBlockQ0,rvectorAreaBlockQ1,rvectorAreaQ0)
    
    call ucd_startGMV (rexport,UCD_FLAG_STANDARD,rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation,&
                       'gmv/u2d_0_simple_grid1.gmv')

    ! Area in den Zellen und Area in den Knoten    
    
    call lsysbl_getbase_double (rvectorAreaBlockQ0,p_Ddata6)
    call ucd_addVariableElementBased (rexport,'AreaCells',UCD_VAR_STANDARD, p_Ddata6) 
    
    call lsysbl_getbase_double (rvectorAreaBlockQ1,p_Ddata7)
    call ucd_addVariableVertexBased (rexport,'AreaNode',UCD_VAR_STANDARD, p_Ddata7) 
    
    call ucd_write (rexport)
    call ucd_release (rexport)

    ! Release the triangulation. 
    call tria_done (rgriddefInfo%p_rhLevels(NLMAX)%rtriangulation)
    call tria_done (rtriangulation)
    
    ! Release Deformation structures
    call griddef_DeformationDone(rgriddefInfo,rgriddefWork)

    
    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
!    call spdiscr_releaseBlockDiscr(rdiscretisation)
    
    
    call lsysbl_releaseVector(rvectorAreaBlockQ0)
    call lsysbl_releaseVector(rvectorAreaBlockQ1) 
    call lsyssc_releaseVector(rvectorAreaQ0) 
    
    ! Finally release the domain, that's it.
    call boundary_release (rboundary)
    
  end subroutine

!******************************************************************************
  
  
  !<subroutine>  
    subroutine tridef2d_monitorfct(DvertexCoords,Dentries)
  
  
    !<description>
      ! In this function we build the nodewise area distribution out 
      ! of an elementwise distribution
    !</description>

    !<inputoutput>
     real(DP), dimension(:,:) :: DvertexCoords
     real(DP), dimension(:) :: Dentries
    !</inputoutput>

    !</subroutine>
    ! local variables
     real(dp),dimension(:,:),allocatable :: Dpoints
     integer(PREC_VERTEXIDX) :: ive,i1,ipoints
     integer :: iMethod
     real(DP) :: Dist,t,dt,dmin
     iMethod = 0
      
     ipoints = ubound(Dentries,1) 
      
      
     select case(iMethod)
       case(0)
         ! loop over all vertices and compute the monitor function
        do ive=1, ipoints
          Dentries(ive) = 0.5_dp + DvertexCoords(1,ive)
         end do
       case(1)
         ! loop over all vertices and compute the monitor function
         do ive=1,ipoints
           Dist = sqrt((0.5_dp - DvertexCoords(1,ive))**2 + (0.5_dp - DvertexCoords(2,ive))**2)
           ! Good now define the monitor function
           Dist = abs(Dist - 0.2_dp)/0.2_dp
           Dist=max(dist,0.1_dp)
           Dist=min(1.0_dp,dist)
           Dentries(ive)=Dist
         end do
       case(2)
        do ive=1,ipoints
          Dentries(ive) = 1.0_dp
        end do
       case(3)
       
        allocate(Dpoints(2,10000))
        dt = 6.28_dp/real(10000)
        t  = 0.0_dp
        do i1=1,10000
         Dpoints(1,i1) = 0.5_dp + 0.1_dp * cos(t)
         Dpoints(2,i1) = 0.5_dp + 0.2_dp * sin(t)
         t = t + dt
        end do
       
        
       
        do ive=1,ipoints
          
          dmin = 10000.0_dp
            
          do i1=1,10000
            Dist = sqrt((Dpoints(1,i1)-DvertexCoords(1,ive))**2 + (Dpoints(2,i1)-DvertexCoords(2,ive))**2)
            dmin =min(Dist,dmin)
          end do
          
          Dentries(ive) = dmin
          
        end do
       case default
     end select
      
  end subroutine

end module
