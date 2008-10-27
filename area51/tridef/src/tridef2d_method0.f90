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

MODULE tridef2d_method0

  USE fsystem
  USE genoutput
  USE storage
  USE linearsolver
  USE boundary
  USE bilinearformevaluation
  USE linearformevaluation
  USE cubature
  USE matrixfilters
  USE vectorfilters
  USE bcassembly
  USE triangulation
  USE spatialdiscretisation
  USE ucd
  USE pprocerror
  USE genoutput
  USE griddeform
    
  USE poisson2d_callback
  
  IMPLICIT NONE

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE tridef2d_simple
  
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
    TYPE(t_boundary) :: rboundary
    
    TYPE(t_griddefWork) :: rgriddefWork
    
    ! An object for saving the triangulation on the domain
    TYPE(t_triangulation) :: rtriangulation

    ! An object specifying the discretisation.
    ! This contains also information about trial/test functions,...
    TYPE(t_blockDiscretisation) :: rdiscretisation
    
    ! In this vector we will save the area of the deformed grid
    TYPE(t_vectorScalar) :: rvectorAreaQ0 

    ! NLMAX receives the level where we want to solve.
    INTEGER :: NLMAX
    
    ! dummy
    INTEGER :: idummy    
    
    ! Output block for UCD output to GMV file
    TYPE(t_ucdExport) :: rexport
    
    ! To these arrays we write the monitorfunction,
    ! elementArea, rhs, etc... in order to
    ! check these values in gmv
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata1
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata2        
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata3 
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata4 
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata5 
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata6 
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata7 
    TYPE(t_vectorBlock) :: rvectorAreaBlockQ0
    TYPE(t_vectorBlock) :: rvectorAreaBlockQ1                 
    
    ! grid deform structure
    TYPE(t_griddefInfo) :: rgriddefInfo

    ! Ok, let's start. 
    !
    ! Here we set the level of the
    ! grid we want to deform
    NLMAX = 6
    
    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    CALL boundary_read_prm(rboundary, './pre/QUAD.prm')
        
    ! Now read in the basic triangulation.
    CALL tria_readTriFile2D (rtriangulation, './pre/QUAD.tri', rboundary)
     
    ! Refine it.
    CALL tria_quickRefine2LevelOrdering (NLMAX-1,rtriangulation,rboundary)
    
    ! And create information about adjacencies and everything one needs from
    ! a triangulation.
    CALL tria_initStandardMeshFromRaw (rtriangulation,rboundary)
    
    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies the blocks in the
    ! solution vector. In this simple problem, we only have one block.
    CALL spdiscr_initBlockDiscr2D (rdiscretisation,1,&
                                   rtriangulation, rboundary)
    
    ! rdiscretisation%Rdiscretisations is a list of scalar discretisation
    ! structures for every component of the solution vector.
    ! Initialise the first element of the list to specify the element
    ! and cubature rule for this solution component:
    CALL spdiscr_initDiscr_simple (rdiscretisation%RspatialDiscr(1), &
                                   EL_E011,CUB_G2X2,rtriangulation, rboundary)
    
    ! Discretisation is set up, now go for the griddeformation
    CALL griddef_deformationInit(rgriddefInfo,rtriangulation,NLMAX,rboundary)                 
    
    ! Deform
    CALL griddef_performDeformation(rgriddefInfo, rgriddefWork,idummy,&
                                     rdiscretisation,&
                                    .TRUE., .FALSE., .FALSE., &
                                    .FALSE., NLMAX, 0, 0,rboundary)    
                 
    ! That's it, rvectorBlock now contains our solution. We can now
    ! start the postprocessing. 
    ! Start UCD export to GMV file:

    CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,rtriangulation,&
                       'gmv/u2d_0_simple1.gmv')
    
    CALL lsyssc_getbase_double (rgriddefWork%rvecGradBlock%RvectorBlock(1),p_Ddata)
    CALL ucd_addVariableVertexBased (rexport,'Ux',UCD_VAR_STANDARD, p_Ddata)
    
    CALL lsyssc_getbase_double (rgriddefWork%rvecGradBlock%RvectorBlock(2),p_Ddata1)
    CALL ucd_addVariableVertexBased (rexport,'Uy',UCD_VAR_STANDARD, p_Ddata1)
    
    CALL lsyssc_getbase_double (rgriddefWork%rSolBlock%RvectorBlock(1),p_Ddata2)
    CALL ucd_addVariableVertexBased (rexport,'Sol',UCD_VAR_STANDARD, p_Ddata2)
    
    CALL lsyssc_getbase_double (rgriddefWork%rrhsBlock%RvectorBlock(1),p_Ddata3)
    CALL ucd_addVariableVertexBased (rexport,'rhs',UCD_VAR_STANDARD, p_Ddata3) 

    CALL lsyssc_getbase_double (rgriddefWork%rvectorAreaBlockQ1%RvectorBlock(1),p_Ddata4)
    CALL ucd_addVariableVertexBased (rexport,'area',UCD_VAR_STANDARD, p_Ddata4) 

    CALL lsyssc_getbase_double (rgriddefWork%rvectorMonFuncQ1%RvectorBlock(1),p_Ddata5)
    CALL ucd_addVariableVertexBased (rexport,'Mon',UCD_VAR_STANDARD, p_Ddata5) 
    
    ! Write the file to disc, that's it.
    CALL ucd_write (rexport)
    CALL ucd_release (rexport)

    CALL griddef_getAreaDeformed(rgriddefInfo,&
             rvectorAreaBlockQ0,rvectorAreaBlockQ1,rvectorAreaQ0)
    
    CALL ucd_startGMV (rexport,UCD_FLAG_STANDARD,rgriddefInfo%rDeftriangulation,&
                       'gmv/u2d_0_simple_grid1.gmv')

    ! Area in den Zellen und Area in den Knoten    
    
    CALL lsysbl_getbase_double (rvectorAreaBlockQ0,p_Ddata6)
    CALL ucd_addVariableElementBased (rexport,'AreaCells',UCD_VAR_STANDARD, p_Ddata6) 
    
    CALL lsysbl_getbase_double (rvectorAreaBlockQ1,p_Ddata7)
    CALL ucd_addVariableVertexBased (rexport,'AreaNode',UCD_VAR_STANDARD, p_Ddata7) 
    
    CALL ucd_write (rexport)
    CALL ucd_release (rexport)
    
    ! Release Deformation structures
    CALL griddef_DeformationDone(rgriddefInfo,rgriddefWork)

    
    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    CALL spdiscr_releaseBlockDiscr(rdiscretisation)
    
    ! Release the triangulation. 
    CALL tria_done (rgriddefInfo%rDeftriangulation)    
    CALL tria_done (rtriangulation)
    
    CALL lsysbl_releaseVector(rvectorAreaBlockQ0)
    CALL lsysbl_releaseVector(rvectorAreaBlockQ1) 
    CALL lsyssc_releaseVector(rvectorAreaQ0) 
    
    ! Finally release the domain, that's it.
    CALL boundary_release (rboundary)
    
  END SUBROUTINE

END MODULE
