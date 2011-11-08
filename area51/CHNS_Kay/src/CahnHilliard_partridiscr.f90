!##############################################################################
!# ****************************************************************************
!# <name> CahnHilliard_partridiscr </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines to read the parametrisation, create the
!# triangulation and set up the discretisation for the heat conduction problem.
!# The following routines can be found here:
!#
!# 0.) CH_initSolution
!#     -> Give the initial solution of Cahn-Hilliard equation
!#
!# 1.) CH_initParamTriang
!#     -> Read .PRM/.TRI files. Generate meshes on all levels.
!#
!# 2.) CH_doneParamTriang
!#     Clean up parametrisation/triangulation, release memory.
!#
!# 3.) CH_initDiscretisation
!#     -> Initialise the spatial discretisation.
!#
!# 4.) CH_doneDiscretisation
!#     -> Clean up the discretisation, release memory.
!# </purpose>
!##############################################################################

module CahnHilliard_partridiscr

  use fsystem
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
  use sortstrategy
  use coarsegridcorrection
  use ucd
  use timestepping
  use genoutput
  use element
  
  use collection
  use paramlist
    
  use CahnHilliard_callback
  use CahnHilliard_basic
  
  IMPLICIT NONE
  
CONTAINS

  ! ***************************************************************************

!<subroutine>
  subroutine CH_initSolution(rCHproblem, rCHvector)
!<description>
  ! Initialises the initial solution vector into rvector. Depending on the settings
  ! in the DAT file this is either zero or read from a file.
  !
  ! The routine assumes that basic mass matrices have already been constructed.
!</description>

!<input>
  ! A problem structure saving problem-dependent information.
  type(t_CHproblem), intent(INOUT) :: rCHproblem
!</input>

!<inputoutput>
  ! The solution vector to be initialised. Must be set up according to the
  ! maximum level NLMAX in rproblem!
  type(t_vectorBlock), intent(INOUT) :: rCHvector
!</inputoutput>

! local variables
  integer :: i
  real(DP), dimension(:,:), pointer ::  p_DvertexCoords
  real(DP), dimension(:), pointer ::  p_vectordata
  real(DP), dimension(:), pointer ::  p_data

    call lsyssc_getbase_double(rCHvector%rvectorBlock(1), p_vectordata)
    call storage_getbase_double2D(rCHproblem%RlevelInfo(&
             rCHproblem%NLMAX)%rtriangulation%h_DvertexCoords,p_DvertexCoords)

    do i=1,rCHproblem%RlevelInfo(rCHproblem%NLMAX)%rtriangulation%NVT
      call CH_iniconPhi(p_DvertexCoords(1,i),p_DvertexCoords(2,i), p_vectordata(i))
    end do

    ! for initial solution of chemical potential
	call lsyssc_getbase_double(rCHvector%rvectorBlock(2), p_vectordata)
    call storage_getbase_double2D(rCHproblem%RlevelInfo(&
             rCHproblem%NLMAX)%rtriangulation%h_DvertexCoords,p_DvertexCoords)

    do i=1,rCHproblem%RlevelInfo(rCHproblem%NLMAX)%rtriangulation%NVT
      call CH_iniconChemP(p_DvertexCoords(1,i),p_DvertexCoords(2,i), p_vectordata(i))
    end do

 
  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine CH_initParamTriang (NLMIN,NLMAX,rCHproblem)
  
!<description>
  ! This routine initialises the parametrisation and triangulation of the
  ! domain. The corresponding .prm/.tri files are read from disc and
  ! the triangulation is refined as described by the parameter ilv.
!</description>

!<input>
  ! Minimum refinement level of the mesh; = coarse grid = level 1
  integer, intent(IN) :: NLMIN
  
  ! Maximum refinement level
  integer, intent(IN) :: NLMAX
!</input>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_CHproblem), intent(INOUT) :: rCHproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i
  
    ! Initialise the level in the problem structure

    rCHproblem%NLMIN = NLMIN
    rCHproblem%NLMAX = NLMAX

    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    call boundary_read_prm(rCHproblem%rboundary, './pre/QUAD.prm')
        

    ! Now read in the basic triangulation.
    call tria_readTriFile2D (rCHproblem%RlevelInfo(rCHproblem%NLMIN)%rtriangulation, &
        './pre/QUAD.tri', rCHproblem%rboundary)

    ! Refine the mesh up to the minimum level
    call tria_quickRefine2LevelOrdering(rCHproblem%NLMIN-1,&
        rCHproblem%RlevelInfo(rCHproblem%NLMIN)%rtriangulation,rCHproblem%rboundary)
    
    ! Create information about adjacencies and everything one needs from
    ! a triangulation. Afterwards, we have the coarse mesh.
    call tria_initStandardMeshFromRaw (&
        rCHproblem%RlevelInfo(rCHproblem%NLMIN)%rtriangulation,rCHproblem%rboundary)
    
    ! Now, refine to level up to nlmax.
    do i=rCHproblem%NLMIN+1,rCHproblem%NLMAX
      call tria_refine2LevelOrdering (rCHproblem%RlevelInfo(i-1)%rtriangulation,&
          rCHproblem%RlevelInfo(i)%rtriangulation, rCHproblem%rboundary)
      call tria_initStandardMeshFromRaw (rCHproblem%RlevelInfo(i)%rtriangulation,&
          rCHproblem%rboundary)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine CH_initDiscretisation (rCHproblem)
  
!<description>
  ! This routine initialises the discretisation structure of the underlying
  ! problem and saves it to the problem structure.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_CHproblem), intent(INOUT), TARGET :: rCHproblem
!</inputoutput>

!</subroutine>

    ! local variables
  
    ! An object for saving the domain:
    type(t_boundary), POINTER :: p_rboundary
    
    ! An object for saving the triangulation on the domain
    type(t_triangulation), POINTER :: p_rtriangulation

    ! MCai
    ! An object for the block discretisation on one level
	! MCai,
	! In CH problem, we have two blocks. If we apply Ciarlet-Raviart type mixed
	! finite element, two blocks have same discretisation. But, it is better to
	! use two pointers to denote two blocks
	
    ! An object for the block discretisation on one level
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
	! Specifically, p_rdiscretisation(1)=p_rdiscretisation_phase
    ! p_rdiscretisation(2)=p_rdiscretisation_chemPoten
    type(t_spatialDiscretisation), pointer :: p_rdiscretisationLaplace, p_rdiscretisationMass
	character(LEN=SYS_NAMELEN) :: sstr
    integer :: i, k, ielementType,icubtemp
    integer(i32) :: icubA,icubB,icubF
    integer(i32) ::  ieltype_A, ieltype_B


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Which discretisation is to use?
    ! Which cubature formula should be used?
    call parlst_getvalue_int (rCHproblem%rparamList,'CH-DISCRETISATION',&
                             'iElementType',ielementType,0)  ! Default Q1 element

    call parlst_getvalue_string (rCHproblem%rparamList,'CH-DISCRETISATION',&
                                 'scubA',sstr,'')
    if (sstr .eq. '') then
    	icubtemp = CUB_G2X2
      call parlst_getvalue_int (rCHproblem%rparamList,'CH-DISCRETISATION',&
                              'icubA',icubtemp,icubtemp)
      icubA = icubtemp
    else
      icubA = cub_igetID(sstr)
    end if

    call parlst_getvalue_string (rCHproblem%rparamList,'CH-DISCRETISATION',&
                              'scubB',sstr,'')
    if (sstr .eq. '') then
    	icubtemp = CUB_G2X2
      call parlst_getvalue_int (rCHproblem%rparamList,'CH-DISCRETISATION',&
                                'icubB',icubtemp,icubtemp)
      icubB = icubtemp
    else
      icubB = cub_igetID(sstr)
    end if

    call parlst_getvalue_string (rCHproblem%rparamList,'CH-DISCRETISATION',&
                               'scubF',sstr,'')
    if (sstr .eq. '') then
    	icubtemp = CUB_G2X2
      call parlst_getvalue_int (rCHproblem%rparamList,'CH-DISCRETISATION',&
                              'icubF',icubtemp,icubtemp)
      icubF = icubtemp
    else
      icubF = cub_igetID(sstr)
    end if


    select case (ielementType)
      case (0)
        ieltype_A = EL_Q1
        ieltype_B = EL_Q1
      case (1)
        ieltype_A = EL_Q2
        ieltype_B = EL_Q2
      case default
      call output_line (&
        'Unknown discretisation: iElementType = '//sys_siL(ielementType,10), &
        OU_CLASS_ERROR,OU_MODE_STD,'CH_initDiscretisation')
      call sys_halt()
    end select

    ! Now set up discrezisation structures on all levels:
    do i=rCHproblem%NLMIN,rCHproblem%NLMAX
      ! Ask the problem structure to give us the boundary and triangulation.
      ! We need it for the discretisation.
      p_rboundary => rCHproblem%rboundary
      p_rtriangulation => rCHproblem%RlevelInfo(i)%rtriangulation
     
      allocate(p_rdiscretisation)
      allocate(p_rdiscretisationLaplace)
      allocate(p_rdiscretisationMass)
      ! Now we can start to initialise the discretisation. At first, set up
      ! a block discretisation structure that specifies the blocks in the
      ! solution vector.
	  
      ! MCai,
      ! In CH problem, we have two blocks. If we apply Ciarlet-Raviart type mixed
      ! finite element, two blocks have same discretisation.
	
      p_rdiscretisation => rCHproblem%RlevelInfo(i)%rdiscretisation

      ! Now we can start to initialise the discretisation. At first, set up
      ! a block discretisation structure that specifies 2 blocks in the
      ! solution vector.
      call spdiscr_initBlockDiscr(&
        p_rdiscretisation,2,p_rtriangulation,p_rboundary)

      ! rdiscretisation%RspatialDiscr is a list of scalar
      ! discretisation structures for every component of the solution vector.
      ! We have a solution vector with two components:
      !  Component 1 = Phase variable
      !  Component 2 = Chemical potential
      ! For simplicity, we set up one discretisation structure for the phase var
      ! then copy the discretisation to chemical potential
      call spdiscr_initDiscr_simple ( &
        p_rdiscretisation%RspatialDiscr(1), &
        ieltype_A,icubA,p_rtriangulation, p_rboundary)
 
      ! Manually set the cubature formula for the RHS as the above routine
      ! uses the same for matrix and vectors.
      p_rdiscretisation%RspatialDiscr(1)% &
        RelementDistr(1)%ccubTypeLinForm = icubF

      ! 2nd discretisation, icubB should be equal to icubA, we do not need it now.
      call spdiscr_initDiscr_simple ( &
        p_rdiscretisation%RspatialDiscr(2), &
        ieltype_B,icubB,p_rtriangulation, p_rboundary)

      ! Manually set the cubature formula for the RHS as the above routine
      ! uses the same for matrix and vectors.
      p_rdiscretisation%RspatialDiscr(2)% &
        RelementDistr(1)%ccubTypeLinForm = icubF

      ! Save the discretisation structure to our local LevelInfo structure
      ! for later use.
      rCHproblem%RlevelInfo(i)%p_rdiscretisation => p_rdiscretisation

      call spdiscr_duplicateDiscrSc (p_rdiscretisation%RspatialDiscr(1), &
          rCHproblem%RlevelInfo(i)%rdiscretisationLaplace,.true.)
      call spdiscr_duplicateDiscrSc (p_rdiscretisation%RspatialDiscr(1), &
          rCHproblem%RlevelInfo(i)%rdiscretisationMass,.true.)

      p_rdiscretisationLaplace => rCHproblem%RlevelInfo(i)%rdiscretisationLaplace
      p_rdiscretisationMass => rCHproblem%RlevelInfo(i)%rdiscretisationMass

      ! Initialise the cubature formula appropriately.
      do k = 1,p_rdiscretisationMass%inumFESpaces
         p_rdiscretisationLaplace%RelementDistr(k)%ccubTypeBilForm = icubA
         p_rdiscretisationMass%RelementDistr(k)%ccubTypeBilForm = icubA
      end do

    end do
                                   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine CH_doneDiscretisation (rCHproblem)
  
!<description>
  ! Releases the discretisation from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_CHproblem), intent(INOUT), TARGET :: rCHproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i

    do i=rCHproblem%NLMAX,rCHproblem%NLMIN,-1
      ! Delete the block discretisation together with the associated
      ! scalar spatial discretisations....
      call spdiscr_releaseBlockDiscr(rCHproblem%RlevelInfo(i)%p_rdiscretisation)
      
      ! and remove the allocated block discretisation structure from the heap.

      ! and remove the allocated block discretisation structure from the heap.
      ! why we do not need to deallocate?
      !  deallocate(rCHproblem%RlevelInfo(i)%p_rdiscretisation)

 !     deallocate(rCHproblem%RlevelInfo(i)%p_rdiscretisation)
      
    end do
    
  end subroutine
    
  ! ***************************************************************************

!<subroutine>

  subroutine CH_doneParamTriang (rCHproblem)
  
!<description>
  ! Releases the triangulation and parametrisation from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_CHproblem), intent(INOUT), TARGET :: rCHproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i

    do i=rCHproblem%NLMAX,rCHproblem%NLMIN,-1
      ! Release the triangulation
      call tria_done (rCHproblem%RlevelInfo(i)%rtriangulation)
    end do
    
    ! Finally release the domain.
    call boundary_release (rCHproblem%rboundary)
    
  end subroutine

end module
