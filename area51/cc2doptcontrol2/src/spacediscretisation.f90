!##############################################################################
!# ****************************************************************************
!# <name> spacediscretisation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the basic spatial discretisation related routines for
!# CC2D. Here, matrix and RHS creation routines can be found as well as
!# routines to initialise/clean up discretisation structures.
!#
!# The following routines can be found here:
!#
!# 1.) fget1LevelDiscretisationNavSt2D
!#     -> Callback routine for spdsc_get1LevelDiscretisation to be used in
!#        routines from fespacehierarchy.
!#
!# 2.) spdsc_get1LevelDiscrNavSt2D
!#     -> Creates a FE discretisation for one spatial level
!#
!# #####
!# 1.) cc_initDiscretisation
!#     -> Initialise the discretisation structure inside of the problem
!#        structure using the parameters from the INI/DAT files.
!#
!# 2.) cc_allocMatVec / cc_generateBasicMat / cc_doneMatVec
!#     -> Allocates/ Generates / Releases memory for vectors/matrices on all levels.
!#
!# 3.) cc_allocStaticMatrices / cc_generateStaticMatrices /
!#     cc_releaseStaticMatrices
!#     -> Allocates/Assembles/Releases matrix entries of static matrices
!#        (Stokes, B) on one level
!#
!# 5.) cc_generateBasicRHS
!#     -> Generates a general RHS vector without any boundary conditions
!#        implemented
!#
!# 6.) cc_doneDiscretisation
!#     -> Cleanup of the underlying discretisation structures
!#
!# 7.) cc_initInitialSolution
!#     -> Init solution vector according to parameters in the DAT file
!#
!# 8.) cc_writeSolution
!#     -> Write solution vector as configured in the DAT file.
!#
!#
!# Furthermore:
!#
!# </purpose>
!##############################################################################

module spacediscretisation

  use fsystem
  use storage
  use genoutput
  use cubature
  use element
  use boundary
  use triangulation
  use spatialdiscretisation
  use collection
  
!  use linearsolver
!  use boundary
!  use bilinearformevaluation
!  use linearformevaluation
!  use cubature
!  use matrixfilters
!  use vectorfilters
!  use bcassembly
!  use triangulation
!  use spatialdiscretisation
!  use coarsegridcorrection
!  use spdiscprojection
!  use nonlinearsolver
!  use paramlist
!  use stdoperators
!  use vectorio
!  use multilevelprojection
!
!  use collection
!  use convection
!  use analyticprojection
!
!  use basicstructures
!  use user_callback
!  use spacematvecassembly
!
!  use analyticsolution
  
  implicit none
  
  private
  
  public :: fget1LevelDiscretisationNavSt2D
  public :: spdsc_get1LevelDiscrNavSt2D
  
contains

  ! ***************************************************************************

!<subroutine>
  
  subroutine fget1LevelDiscretisationNavSt2D(rtriangulation,rdiscr,rboundary,rcollection)

!<description>
  ! Callback routine wrapper for spdsc_get1LevelDiscrNavSt2D. Allows to create
  ! a discretisation with the routines from fespacehierarchy.
  ! Expects the following information in rcollection:
  !   rcollection%IquickAccess(1) = ieltype
  !   rcollection%IquickAccess(2) = nequations (=3:primal space. =6: primal+dual space)
  !   rcollection%IquickAccess(3) = SPDISC_CUB_AUTOMATIC or icubA
  !   rcollection%IquickAccess(4) = SPDISC_CUB_AUTOMATIC or icubB
  !   rcollection%IquickAccess(5) = SPDISC_CUB_AUTOMATIC or icubF
!</description>

  type(t_triangulation), intent(in) :: rtriangulation
  type(t_blockDiscretisation), intent(out) :: rdiscr
  type(t_collection), intent(inout), optional :: rcollection
  type(t_boundary), intent(in), optional :: rboundary
    
!</subroutine>
    
    ! Create the discretisation corresponding to the triangulation
    ! and element type ieltype.
    call spdsc_get1LevelDiscrNavSt2D (rboundary,rtriangulation,&
      rcollection%IquickAccess(2),rdiscr,&
      rcollection%IquickAccess(1),&
      int(rcollection%IquickAccess(3),I32),&
      int(rcollection%IquickAccess(4),I32),&
      int(rcollection%IquickAccess(5),I32))
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdsc_get1LevelDiscrNavSt2D (rboundary,rtriangulation,&
      nequations,rdiscretisation,ieltype,ccubA,ccubB,ccubF)
  
!<description>
  ! Sets up the discretisation structures for the (primal) velocity and pressure
  ! as specified by the parameters in the parameter list rparList.
!</description>

!<input>
  ! Number of equations.
  ! =3 for primal or dual discretisation structure.
  ! =6 for coupled primal/dual discretisation.
  integer, intent(IN) :: nequations
  
  ! Underlying boundary.
  type(t_boundary), intent(IN), target :: rboundary
  
  ! Underlying triangulation
  type(t_triangulation), intent(IN), target :: rtriangulation
  
  ! Element type.
  ! 0 = Q1~(E031) / Q1~(E031) / Q0
  ! 1 = Q1~(E030) / Q1~(E030) / Q0
  ! 2 = Q1~(EM31) / Q1~(EM31) / Q0
  ! 3 = Q1~(EM30) / Q1~(EM30) / Q0 = standard
  ! 4 = Q2 (E013) / Q2 (E013) / QP1
  ! 5 = Q1~(EM30) / Q1~(EM30) / Q0 unpivoted (much faster than 3 but less stable)
  ! 6 = Q1~(EM30) / Q1~(EM30) / Q0 unscaled (slightly faster than 3 but less stable)
  integer, intent(in) :: ieltype
  
  ! OPTIONAL: Cubature type for constructing bilinear forms.
  ! =SPDISC_CUB_AUTOMATIC: use default.
  integer(I32), intent(in), optional :: ccubA

  ! OPTIONAL: Cubature type for constructing bilinear forms on mixed matrices.
  ! =SPDISC_CUB_AUTOMATIC: use default.
  integer(I32), intent(in), optional :: ccubB
  
  ! OPTIONAL: Cubature type for constructing linear forms.
  ! =SPDISC_CUB_AUTOMATIC: use default.
  integer(I32), intent(in), optional :: ccubF
!</input>

!<inputoutput>
  ! The discretisation structure to be set up.
  type(t_blockDiscretisation), intent(INOUT) :: rdiscretisation
!</inputoutput>

!</subroutine>

    integer(I32) :: icubA,icubB,icubF
    integer :: i,j

    icubA = SPDISC_CUB_AUTOMATIC
    icubB = SPDISC_CUB_AUTOMATIC
    icubF = SPDISC_CUB_AUTOMATIC

    if (present(ccubA)) icubA = ccubA
    if (present(ccubB)) icubB = ccubB
    if (present(ccubF)) icubF = ccubF

!    if (present(rparlist)) then
!      ! Get parameters about the discretisation from the parameter list
!      call parlst_getvalue_int (rparlist,'CC-DISCRETISATION',&
!                                'iElementType',ielementType,3)
!
!      if (present(ieltype)) ielementType = ieltype
!
!      call parlst_getvalue_string (rparlist,'CC-DISCRETISATION',&
!                                  'scubStokes',sstr,'')
!      if (sstr .eq. '') then
!        call parlst_getvalue_int (rparlist,'CC-DISCRETISATION',&
!                                  'icubStokes',icubA,CUB_G2X2)
!      else
!        icubA = cub_igetID(sstr)
!      end if
!
!      call parlst_getvalue_string (rparlist,'CC-DISCRETISATION',&
!                                  'scubB',sstr,'')
!      if (sstr .eq. '') then
!        call parlst_getvalue_int (rparlist,'CC-DISCRETISATION',&
!                                  'icubB',icubB,CUB_G2X2)
!      else
!        icubB = cub_igetID(sstr)
!      end if
!
!      call parlst_getvalue_string (rparlist,'CC-DISCRETISATION',&
!                                  'scubF',sstr,'')
!      if (sstr .eq. '') then
!        call parlst_getvalue_int (rparlist,'CC-DISCRETISATION',&
!                                  'icubF',icubF,CUB_G2X2)
!      else
!        icubF = cub_igetID(sstr)
!      end if
!    else
!      icubA = SPDISC_CUB_AUTOMATIC
!      icubB = SPDISC_CUB_AUTOMATIC
!      icubF = SPDISC_CUB_AUTOMATIC
!    end if

    call spdiscr_initBlockDiscr (rdiscretisation,nequations,&
                                 rtriangulation, rboundary)

    select case (ieltype)
    case (0)
      ! rdiscretisation%Rdiscretisations is a list of scalar
      ! discretisation structures for every component of the solution vector.
      ! We have a solution vector with three components:
      !  Component 1 = X-velocity
      !  Component 2 = Y-velocity
      !  Component 3 = Pressure
      ! For simplicity, we set up one discretisation structure for the
      ! velocity...
      call spdiscr_initDiscr_simple ( &
                  rdiscretisation%RspatialDiscr(1), &
                  EL_E031,icubA, &
                  rtriangulation, rboundary)
                  
      ! Manually set the cubature formula for the RHS as the above routine
      ! uses the same for matrix and vectors.
      rdiscretisation%RspatialDiscr(1)% &
        RelementDistr(1)%ccubTypeLinForm = icubF
                  
      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd component (Y-velocity). This needs no additional memory,
      ! as both structures will share the same dynamic information afterwards,
      ! but we have to be careful when releasing the discretisation structures
      ! at the end of the program!
      ! rdiscretisation%RspatialDiscr(2) = &
      !   rdiscretisation%RspatialDiscr(1)
      ! New implementation: Use spdiscr_duplicateDiscrSc!
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1), &
          rdiscretisation%RspatialDiscr(2), .true.)
  
      ! For the pressure (3rd component), we set up a separate discretisation
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(2), &
          EL_Q0, icubB, rdiscretisation%RspatialDiscr(3))

    case (1)
      ! rdiscretisation%Rdiscretisations is a list of scalar
      ! discretisation structures for every component of the solution vector.
      ! We have a solution vector with three components:
      !  Component 1 = X-velocity
      !  Component 2 = Y-velocity
      !  Component 3 = Pressure
      ! For simplicity, we set up one discretisation structure for the
      ! velocity...
      call spdiscr_initDiscr_simple ( &
                  rdiscretisation%RspatialDiscr(1), &
                  EL_E030,icubA, &
                  rtriangulation, rboundary)
                  
      ! Manually set the cubature formula for the RHS as the above routine
      ! uses the same for matrix and vectors.
      rdiscretisation%RspatialDiscr(1)% &
        RelementDistr(1)%ccubTypeLinForm = icubF
                  
      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd component (Y-velocity). This needs no additional memory,
      ! as both structures will share the same dynamic information afterwards,
      ! but we have to be careful when releasing the discretisation structures
      ! at the end of the program!
      ! rdiscretisation%RspatialDiscr(2) = &
      !   rdiscretisation%RspatialDiscr(1)
      ! New implementation: Use spdiscr_duplicateDiscrSc!
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1), &
          rdiscretisation%RspatialDiscr(2), .true.)
  
      ! For the pressure (3rd component), we set up a separate discretisation
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(2), &
          EL_Q0, icubB, rdiscretisation%RspatialDiscr(3))

    case (2)
      ! rdiscretisation%Rdiscretisations is a list of scalar
      ! discretisation structures for every component of the solution vector.
      ! We have a solution vector with three components:
      !  Component 1 = X-velocity
      !  Component 2 = Y-velocity
      !  Component 3 = Pressure
      ! For simplicity, we set up one discretisation structure for the
      ! velocity...
      call spdiscr_initDiscr_simple ( &
                  rdiscretisation%RspatialDiscr(1), &
                  EL_EM31,icubA, &
                  rtriangulation, rboundary)
                  
      ! Manually set the cubature formula for the RHS as the above routine
      ! uses the same for matrix and vectors.
      rdiscretisation%RspatialDiscr(1)% &
        RelementDistr(1)%ccubTypeLinForm = icubF
                  
      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd component (Y-velocity). This needs no additional memory,
      ! as both structures will share the same dynamic information afterwards,
      ! but we have to be careful when releasing the discretisation structures
      ! at the end of the program!
      ! rdiscretisation%RspatialDiscr(2) = &
      !   rdiscretisation%RspatialDiscr(1)
      ! New implementation: Use spdiscr_duplicateDiscrSc!
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1), &
          rdiscretisation%RspatialDiscr(2), .true.)
  
      ! For the pressure (3rd component), we set up a separate discretisation
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(2), &
          EL_Q0, icubB, rdiscretisation%RspatialDiscr(3))

    case (3)
      ! rdiscretisation%Rdiscretisations is a list of scalar
      ! discretisation structures for every component of the solution vector.
      ! We have a solution vector with three components:
      !  Component 1 = X-velocity
      !  Component 2 = Y-velocity
      !  Component 3 = Pressure
      ! For simplicity, we set up one discretisation structure for the
      ! velocity...
      call spdiscr_initDiscr_simple ( &
                  rdiscretisation%RspatialDiscr(1), &
                  EL_EM30,icubA, &
                  rtriangulation, rboundary)
                  
      ! Manually set the cubature formula for the RHS as the above routine
      ! uses the same for matrix and vectors.
      rdiscretisation%RspatialDiscr(1)% &
        RelementDistr(1)%ccubTypeLinForm = icubF
                  
      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd component (Y-velocity). This needs no additional memory,
      ! as both structures will share the same dynamic information afterwards,
      ! but we have to be careful when releasing the discretisation structures
      ! at the end of the program!
      ! rdiscretisation%RspatialDiscr(2) = &
      !   rdiscretisation%RspatialDiscr(1)
      ! New implementation: Use spdiscr_duplicateDiscrSc!
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1), &
          rdiscretisation%RspatialDiscr(2), .true.)
  
      ! For the pressure (3rd component), we set up a separate discretisation
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(2), &
          EL_Q0, icubB, rdiscretisation%RspatialDiscr(3))

    case (4)
      ! rdiscretisation%Rdiscretisations is a list of scalar
      ! discretisation structures for every component of the solution vector.
      ! We have a solution vector with three components:
      !  Component 1 = X-velocity
      !  Component 2 = Y-velocity
      !  Component 3 = Pressure
      ! For simplicity, we set up one discretisation structure for the
      ! velocity...
      call spdiscr_initDiscr_simple ( &
                  rdiscretisation%RspatialDiscr(1), &
                  EL_Q2,icubA, &
                  rtriangulation, rboundary)
                  
      ! Manually set the cubature formula for the RHS as the above routine
      ! uses the same for matrix and vectors.
      rdiscretisation%RspatialDiscr(1)% &
        RelementDistr(1)%ccubTypeLinForm = icubF
                  
      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd component (Y-velocity). This needs no additional memory,
      ! as both structures will share the same dynamic information afterwards,
      ! but we have to be careful when releasing the discretisation structures
      ! at the end of the program!
      ! rdiscretisation%RspatialDiscr(2) = &
      !   rdiscretisation%RspatialDiscr(1)
      ! New implementation: Use spdiscr_duplicateDiscrSc!
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1), &
          rdiscretisation%RspatialDiscr(2), .true.)
  
      ! For the pressure (3rd component), we set up a separate discretisation
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(2), &
          EL_QP1NP, icubB, rdiscretisation%RspatialDiscr(3))
                  
    case (5)
      ! rdiscretisation%Rdiscretisations is a list of scalar
      ! discretisation structures for every component of the solution vector.
      ! We have a solution vector with three components:
      !  Component 1 = X-velocity
      !  Component 2 = Y-velocity
      !  Component 3 = Pressure
      ! For simplicity, we set up one discretisation structure for the
      ! velocity...
      call spdiscr_initDiscr_simple ( &
                  rdiscretisation%RspatialDiscr(1), &
                  EL_EM30_UNPIVOTED,icubA, &
                  rtriangulation, rboundary)
                  
      ! Manually set the cubature formula for the RHS as the above routine
      ! uses the same for matrix and vectors.
      rdiscretisation%RspatialDiscr(1)% &
        RelementDistr(1)%ccubTypeLinForm = icubF
                  
      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd component (Y-velocity). This needs no additional memory,
      ! as both structures will share the same dynamic information afterwards,
      ! but we have to be careful when releasing the discretisation structures
      ! at the end of the program!
      ! rdiscretisation%RspatialDiscr(2) = &
      !   rdiscretisation%RspatialDiscr(1)
      ! New implementation: Use spdiscr_duplicateDiscrSc!
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1), &
          rdiscretisation%RspatialDiscr(2), .true.)
  
      ! For the pressure (3rd component), we set up a separate discretisation
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(2), &
          EL_Q0, icubB, rdiscretisation%RspatialDiscr(3))

    case (6)
      ! rdiscretisation%Rdiscretisations is a list of scalar
      ! discretisation structures for every component of the solution vector.
      ! We have a solution vector with three components:
      !  Component 1 = X-velocity
      !  Component 2 = Y-velocity
      !  Component 3 = Pressure
      ! For simplicity, we set up one discretisation structure for the
      ! velocity...
      call spdiscr_initDiscr_simple ( &
                  rdiscretisation%RspatialDiscr(1), &
                  EL_EM30_UNSCALED,icubA, &
                  rtriangulation, rboundary)
                  
      ! Manually set the cubature formula for the RHS as the above routine
      ! uses the same for matrix and vectors.
      rdiscretisation%RspatialDiscr(1)% &
        RelementDistr(1)%ccubTypeLinForm = icubF
                  
      ! ...and copy this structure also to the discretisation structure
      ! of the 2nd component (Y-velocity). This needs no additional memory,
      ! as both structures will share the same dynamic information afterwards,
      ! but we have to be careful when releasing the discretisation structures
      ! at the end of the program!
      ! rdiscretisation%RspatialDiscr(2) = &
      !   rdiscretisation%RspatialDiscr(1)
      ! New implementation: Use spdiscr_duplicateDiscrSc!
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1), &
          rdiscretisation%RspatialDiscr(2), .true.)
  
      ! For the pressure (3rd component), we set up a separate discretisation
      ! structure, as this uses different finite elements for trial and test
      ! functions.
      call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(2), &
          EL_Q0, icubB, rdiscretisation%RspatialDiscr(3))

    case default
      print *,'Unknown discretisation: ieltype = ',ieltype
      stop
    end select
    
    if (nequations .eq. 6) then

      ! -----------------------------------------------------------------------
      ! Optimal control extension
      ! -----------------------------------------------------------------------
      ! The variables 4,5,6 are discretised like the variables 1,2,3
      ! in the discretisation structure.
      ! Variable 1 = x-velocity  -->  Variable 4 = dual x-velocity
      ! Variable 2 = y-velocity  -->  Variable 5 = dual y-velocity
      ! Variable 3 = pressure    -->  Variable 6 = dual pressure
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(1),&
                                     rdiscretisation%RspatialDiscr(4))
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(2),&
                                     rdiscretisation%RspatialDiscr(5))
      call spdiscr_duplicateDiscrSc (rdiscretisation%RspatialDiscr(3),&
                                     rdiscretisation%RspatialDiscr(6))
          
    end if
    
    ! Post correction of the cubature formulas: If icubF=automatic,
    ! put the standard cubature formulas there.
    if (icubF .eq. SPDISC_CUB_AUTOMATIC) then
      do i=1,size(rdiscretisation%RspatialDiscr)
        do j=1,size(rdiscretisation%RspatialDiscr(i)%RelementDistr)
          rdiscretisation%RspatialDiscr(i)%RelementDistr(j)%ccubTypeLinForm = &
              spdiscr_getStdCubature(&
              rdiscretisation%RspatialDiscr(i)%RelementDistr(j)%celement)
        end do
      end do
    end if
    
  end subroutine

!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine cc_initDiscretisation (rproblem)
!
!!<description>
!  ! This routine initialises the discretisation structure of the underlying
!  ! problem and saves it to the problem structure.
!!</description>
!
!!<inputoutput>
!  ! A problem structure saving problem-dependent information.
!  type(t_problem), intent(INOUT), target :: rproblem
!!</inputoutput>
!
!!</subroutine>
!
!  ! local variables
!  integer :: I,k
!
!  ! Number of equations in our problem.
!  ! velocity+velocity+pressure + dual velocity+dual velocity+dual pressure = 6
!  integer, parameter :: nequations = 6
!
!    ! An object for saving the triangulation on the domain
!    type(t_triangulation), pointer :: p_rtriangulation
!
!    ! An object for the block discretisation on one level
!    type(t_blockDiscretisation), pointer :: p_rdiscretisation
!    type(t_spatialDiscretisation), pointer :: p_rdiscretisationMass
!    type(t_spatialDiscretisation), pointer :: p_rdiscretisationMassPressure
!
!    integer :: icubM
!    character(LEN=SYS_NAMELEN) :: sstr
!
!    ! Set up discrezisation structures on all levels:
!
!    do i=rproblem%NLMIN,rproblem%NLMAX
!      ! Ask the problem structure to give us the boundary and triangulation.
!      ! We need it for the discretisation.
!      p_rtriangulation => rproblem%RlevelInfo(i)%p_rtriangulation
!
!      ! Now we can start to initialise the discretisation. At first, set up
!      ! a block discretisation structure that specifies the blocks in the
!      ! solution vector.
!      call spdsc_get1LevelDiscrNavSt2D (&
!          rproblem%rboundary,p_rtriangulation,nequations,&
!          rproblem%RlevelInfo(i)%rdiscretisation,rparlist=rproblem%rparamList)
!
!      ! Save the discretisation structure to our local LevelInfo structure
!      ! for later use.
!      p_rdiscretisation => rproblem%RlevelInfo(i)%rdiscretisation
!
!      ! -----------------------------------------------------------------------
!      ! Separated discretisation structures for primal and dual problem
!      ! -----------------------------------------------------------------------
!      ! Create a separate block discretisation structure, only for the primal
!      ! space. Make a copy of the one we have, i.e. reuse it.
!      ! Mark the substructures as being a copy from another discretisation structure,
!      ! sharing all information with part 1..3 of the global problem.
!      call spdiscr_deriveBlockDiscr (rproblem%RlevelInfo(i)%rdiscretisation, &
!          rproblem%RlevelInfo(i)%rdiscretisationPrimal, 1,3)
!
!      ! -----------------------------------------------------------------------
!      ! Mass matrices
!      ! -----------------------------------------------------------------------
!
!      ! Initialise a discretisation structure for the mass matrix.
!      ! Copy the discretisation structure of the first (Stokes) block
!      ! and replace the cubature-formula identifier by that which is to be
!      ! used for the mass matrix.
!      call spdiscr_duplicateDiscrSc (p_rdiscretisation%RspatialDiscr(1),&
!          rproblem%RlevelInfo(i)%rstaticInfo%rdiscretisationMass,.true.)
!      call spdiscr_duplicateDiscrSc (p_rdiscretisation%RspatialDiscr(3), &
!          rproblem%RlevelInfo(i)%rstaticInfo%rdiscretisationMassPressure,.true.)
!
!      p_rdiscretisationMass => rproblem%RlevelInfo(i)%rstaticInfo%rdiscretisationMass
!      p_rdiscretisationMassPressure => &
!          rproblem%RlevelInfo(i)%rstaticInfo%rdiscretisationMassPressure
!
!      call parlst_getvalue_string (rproblem%rparamList,'CC-DISCRETISATION',&
!                                  'scubStokes',sstr,'')
!      if (sstr .eq. '') then
!        call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
!                                  'icubStokes',icubM,CUB_G2X2)
!      else
!        icubM = cub_igetID(sstr)
!      end if
!
!      ! Initialise the cubature formula appropriately.
!      do k = 1,p_rdiscretisationMass%inumFESpaces
!        p_rdiscretisationMass%RelementDistr(k)%ccubTypeBilForm = icubM
!        p_rdiscretisationMassPressure%RelementDistr(k)%ccubTypeBilForm = icubM
!      end do
!
!    end do
!
!  end subroutine
!
!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine cc_doneDiscretisation (rproblem)
!
!!<description>
!  ! Releases the discretisation from the heap.
!!</description>
!
!!<inputoutput>
!  ! A problem structure saving problem-dependent information.
!  type(t_problem), intent(INOUT), target :: rproblem
!!</inputoutput>
!
!!</subroutine>
!
!  ! local variables
!  integer :: i
!
!    do i=rproblem%NLMAX,rproblem%NLMIN,-1
!      ! Remove the main block discretisation structure.
!      call spdiscr_releaseBlockDiscr(rproblem%RlevelInfo(i)%rdiscretisation,.true.)
!
!      ! Release the block discretisation structures of the primal and dual
!      ! space.
!      call spdiscr_releaseBlockDiscr(rproblem%RlevelInfo(i)%rdiscretisationPrimal,.true.)
!
!      ! -----------------------------------------------------------------------
!      ! Mass matrix problem
!      ! -----------------------------------------------------------------------
!
!      ! Release the mass matrix discretisation.
!      call spdiscr_releaseDiscr(rproblem%RlevelInfo(i)%rstaticInfo%rdiscretisationMass)
!
!    end do
!
!  end subroutine
!
!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine cc_allocMatVec (rproblem)
!
!!<description>
!  ! Allocates memory for all matrices and vectors of the problem on the heap
!  ! by evaluating the parameters in the problem structure.
!  ! Matrices/vectors of global importance are added to the collection
!  ! structure of the problem, given in rproblem. Matrix/vector entries
!  ! are not calculated, they are left uninitialised.\\\\
!!</description>
!
!!<inputoutput>
!  ! A problem structure saving problem-dependent information.
!  type(t_problem), intent(INOUT), target :: rproblem
!
!  ! A vector structure for the solution vector. The structure is initialised,
!  ! memory is allocated for the data entries.
!  type(t_vectorBlock), intent(INOUT) :: rvector
!
!  ! A vector structure for the RHS vector. The structure is initialised,
!  ! memory is allocated for the data entries.
!  type(t_vectorBlock), intent(INOUT) :: rrhs
!!</inputoutput>
!
!!</subroutine>
!
!  ! local variables
!  integer :: i
!
!    ! Initialise all levels...
!    do i=rproblem%NLMIN,rproblem%NLMAX
!
!      ! Generate the matrices on this level.
!      if (i .eq. rproblem%NLMIN) then
!        call cc_allocStaticMatrices (rproblem,&
!            rproblem%RlevelInfo(i)%rdiscretisation,&
!            rproblem%RlevelInfo(i)%rstaticInfo)
!      else
!        call cc_allocStaticMatrices (rproblem,&
!            rproblem%RlevelInfo(i)%rdiscretisation,&
!            rproblem%RlevelInfo(i)%rstaticInfo,&
!            rproblem%RlevelInfo(i-1)%rdiscretisation)
!      end if
!
!    end do
!
!  end subroutine
!
!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine cc_generateBasicMat (rproblem)
!
!!<description>
!  ! Calculates the entries of all static matrices (Mass, B,...) on all levels.
!  !
!  ! Memory for those matrices must have been allocated before with
!  ! allocMatVec!
!!</description>
!
!!<inputoutput>
!  ! A problem structure saving problem-dependent information.
!  type(t_problem), intent(INOUT) :: rproblem
!!</inputoutput>
!
!!</subroutine>
!
!    ! local variables
!    integer :: i
!
!    do i=rproblem%NLMIN,rproblem%NLMAX
!      call cc_generateStaticMatrices (rproblem,&
!          rproblem%RlevelInfo(i)%rdiscretisation,&
!          rproblem%RlevelInfo(i)%rstaticInfo)
!    end do
!
!  end subroutine
!
!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine cc_doneMatVec (rproblem,rvector,rrhs)
!
!!<description>
!  ! Releases system matrix and vectors.
!!</description>
!
!!<inputoutput>
!  ! A problem structure saving problem-dependent information.
!  type(t_problem), intent(INOUT), target :: rproblem
!
!  ! A vector structure for the solution vector. The structure is cleaned up,
!  ! memory is released.
!  type(t_vectorBlock), intent(INOUT) :: rvector
!
!  ! A vector structure for the RHS vector. The structure is cleaned up,
!  ! memory is released.
!  type(t_vectorBlock), intent(INOUT) :: rrhs
!!</inputoutput>
!
!!</subroutine>
!
!    integer :: i
!
!    ! Release matrices and vectors on all levels
!    do i=rproblem%NLMAX,rproblem%NLMIN,-1
!
!      ! Release the static matrices.
!      call cc_releaseStaticMatrices (rproblem%RlevelInfo(i)%rstaticInfo)
!
!    end do
!
!    ! Delete solution/RHS vector
!    call lsysbl_releaseVector (rvector)
!    call lsysbl_releaseVector (rrhs)
!
!  end subroutine
!
!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine cc_allocStaticMatrices (rproblem,rdiscretisation,rstaticInfo,&
!      rdiscretisationCoarse)
!
!!<description>
!  ! Allocates memory and generates the structure of all static matrices
!  ! in rstaticInfo.
!!</description>
!
!!<inputoutput>
!  ! A problem structure saving problem-dependent information.
!  type(t_problem), intent(inout) :: rproblem
!
!  ! Discretisation structure that defines how to discretise the different
!  ! operators.
!  type(t_blockDiscretisation), intent(in), target :: rdiscretisation
!
!  ! A t_staticLevelInfo structure. The static matrices in this structure are generated.
!  type(t_staticLevelInfo), intent(inout), target :: rstaticInfo
!
!  ! OPTIONAL: Discretisation structure of the level below the level
!  ! identified by rdiscretisation. Must be specified on all levels except
!  ! for the coarse mesh.
!  type(t_blockDiscretisation), intent(in), optional :: rdiscretisationCoarse
!!</inputoutput>
!
!!</subroutine>
!
!  ! local variables
!  integer :: cmatBuildType
!  integer :: i
!
!    ! When the jump stabilisation is used, we have to create an extended
!    ! matrix stencil!
!    cmatBuildType = BILF_MATC_ELEMENTBASED
!
!    call parlst_getvalue_int (rproblem%rparamList, 'CC-DISCRETISATION', &
!                              'IUPWIND1', i)
!    if ((i .eq. CCMASM_STAB_EDGEORIENTED) .or.&
!        (i .eq. CCMASM_STAB_EDGEORIENTED2) .or. &
!        (i .eq. CCMASM_STAB_EDGEORIENTED3)) &
!       cmatBuildType = BILF_MATC_EDGEBASED
!
!    call parlst_getvalue_int (rproblem%rparamList, 'CC-DISCRETISATION', &
!                              'IUPWIND2', i)
!    if ((i .eq. CCMASM_STAB_EDGEORIENTED) .or.&
!        (i .eq. CCMASM_STAB_EDGEORIENTED2) .or. &
!        (i .eq. CCMASM_STAB_EDGEORIENTED3)) &
!       cmatBuildType = BILF_MATC_EDGEBASED
!
!    ! -----------------------------------------------------------------------
!    ! Basic (Navier-) Stokes problem
!    ! -----------------------------------------------------------------------
!
!    ! The global system looks as follows:
!    !
!    !    ( A11       B1  M/a          )  ( u )  =  ( f1)
!    !    (      A22  B2       M/a     )  ( v )     ( f2)
!    !    ( B1^T B2^T                  )  ( p )     ( fp)
!    !    ( -M            A44  A45  B1 )  ( l1)     ( z1)
!    !    (      -M       A54  A55  B2 )  ( l2)     ( z2)
!    !    (               B1^T B2^T    )  ( xi)     ( 0 )
!    !
!    ! with A = L + nonlinear Convection. We compute in advance
!    ! a standard Stokes matrix L which can be added later to the
!    ! convection matrix, resulting in the nonlinear system matrix.
!    !
!    ! At first, we create 'template' matrices that define the structure
!    ! of each of the submatrices in that global system. These matrices
!    ! are later used to 'derive' the actual Laplace/Stokes matrices
!    ! by sharing the structure (KCOL/KLD).
!    !
!    ! Get a pointer to the template FEM matrix. This is used for the
!    ! Laplace/Stokes matrix and probably for the mass matrix.
!
!    ! Create the matrix structure
!    call bilf_createMatrixStructure (&
!              rdiscretisation%RspatialDiscr(1),LSYSSC_MATRIX9,&
!              rstaticInfo%rmatrixTemplateFEM,cconstrType=cmatBuildType)
!
!    ! Create the matrices structure of the pressure using the 3rd
!    ! spatial discretisation structure in rdiscretisation%RspatialDiscr.
!    call bilf_createMatrixStructure (&
!              rdiscretisation%RspatialDiscr(3),LSYSSC_MATRIX9,&
!              rstaticInfo%rmatrixTemplateFEMPressure)
!
!    ! Create the matrices structure of the pressure using the 3rd
!    ! spatial discretisation structure in rdiscretisation%RspatialDiscr.
!    call bilf_createMatrixStructure (&
!              rdiscretisation%RspatialDiscr(3),LSYSSC_MATRIX9,&
!              rstaticInfo%rmatrixTemplateGradient,&
!              rdiscretisation%RspatialDiscr(1))
!
!    ! Transpose the B-structure to get the matrix template for the
!    ! divergence matrices.
!    call lsyssc_transposeMatrix (rstaticInfo%rmatrixTemplateGradient,&
!        rstaticInfo%rmatrixTemplateDivergence,LSYSSC_TR_STRUCTURE)
!
!    ! Ok, now we use the matrices from above to create the actual submatrices
!    ! that are used in the global system.
!    !
!    ! Connect the Stokes matrix to the template FEM matrix such that they
!    ! use the same structure.
!    !
!    ! Don't create a content array yet, it will be created by
!    ! the assembly routines later.
!    call lsyssc_duplicateMatrix (rstaticInfo%rmatrixTemplateFEM,&
!                rstaticInfo%rmatrixStokes,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
!
!    ! Allocate memory for the entries; don't initialise the memory.
!    call lsyssc_allocEmptyMatrix (rstaticInfo%rmatrixStokes,LSYSSC_SETM_UNDEFINED)
!
!    ! In the global system, there are two coupling matrices B1 and B2.
!    ! Both have the structure of the template gradient matrix.
!    ! So connect the two B-matrices to the template gradient matrix
!    ! such that they share the same structure.
!    ! Create the matrices structure of the pressure using the 3rd
!    ! spatial discretisation structure in rdiscretisation%RspatialDiscr.
!    !
!    ! Don't create a content array yet, it will be created by
!    ! the assembly routines later.
!    ! Allocate memory for the entries; don't initialise the memory.
!
!    call lsyssc_duplicateMatrix (rstaticInfo%rmatrixTemplateGradient,&
!        rstaticInfo%rmatrixB1,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
!
!    call lsyssc_duplicateMatrix (rstaticInfo%rmatrixTemplateGradient,&
!        rstaticInfo%rmatrixB2,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
!
!    call lsyssc_allocEmptyMatrix (rstaticInfo%rmatrixB1,LSYSSC_SETM_UNDEFINED)
!    call lsyssc_allocEmptyMatrix (rstaticInfo%rmatrixB2,LSYSSC_SETM_UNDEFINED)
!
!    ! Set up memory for the divergence matrices D1 and D2.
!    call lsyssc_duplicateMatrix (rstaticInfo%rmatrixTemplateDivergence,&
!        rstaticInfo%rmatrixD1,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
!
!    call lsyssc_duplicateMatrix (rstaticInfo%rmatrixTemplateDivergence,&
!        rstaticInfo%rmatrixD2,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
!
!    call lsyssc_allocEmptyMatrix (rstaticInfo%rmatrixD1,LSYSSC_SETM_UNDEFINED)
!    call lsyssc_allocEmptyMatrix (rstaticInfo%rmatrixD2,LSYSSC_SETM_UNDEFINED)
!
!    ! The D1^T and D2^T matrices are by default the same as B1 and B2.
!    ! These matrices may be different for special VANCA variants if
!    ! B1 and B2 is different from D1 and D2 (which is actually a rare case).
!    call lsyssc_duplicateMatrix (rstaticInfo%rmatrixB1,&
!        rstaticInfo%rmatrixD1T,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!
!    call lsyssc_duplicateMatrix (rstaticInfo%rmatrixB2,&
!        rstaticInfo%rmatrixD2T,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!
!
!  end subroutine
!
!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine cc_generateStaticMatrices (rproblem,rdiscretisation,rstaticInfo)
!
!!<description>
!  ! Calculates entries of all static matrices (Stokes, B-matrices,...)
!  ! in the specified problem structure, i.e. the entries of all matrices
!  ! that do not change during the computation or which serve as a template for
!  ! generating other matrices.
!  !
!  ! Memory for those matrices must have been allocated before with
!  ! allocMatVec!
!!</description>
!
!!<inputoutput>
!  ! A problem structure saving problem-dependent information.
!  type(t_problem), intent(inout) :: rproblem
!
!  ! Discretisation structure that defines how to discretise the different
!  ! operators.
!  type(t_blockDiscretisation), intent(in), target :: rdiscretisation
!
!  ! A t_staticLevelInfo structure. The static matrices in this structure are generated.
!  type(t_staticLevelInfo), intent(inout), target :: rstaticInfo
!!</inputoutput>
!
!!</subroutine>
!
!    ! Stabilisation stuff
!    type(t_jumpStabilisation) :: rjumpStabil
!    integer :: iupwind1,iupwind2
!    real(DP) :: dupsam1,dupsam2
!
!    ! Initialise the collection for the assembly process with callback routines.
!    ! Basically, this stores the simulation time in the collection if the
!    ! simulation is nonstationary.
!    call user_initCollectForAssembly (rproblem,0.0_DP,rproblem%rcollection)
!
!    ! -----------------------------------------------------------------------
!    ! Basic (Navier-) Stokes problem
!    ! -----------------------------------------------------------------------
!
!    ! The global system looks as follows:
!    !
!    !    ( A11       B1  M/a          )  ( u )  =  ( f1)
!    !    (      A22  B2       M/a     )  ( v )     ( f2)
!    !    ( B1^T B2^T                  )  ( p )     ( fp)
!    !    ( -M            A44  A45  B1 )  ( l1)     ( z1)
!    !    (      -M       A54  A55  B2 )  ( l2)     ( z2)
!    !    (               B1^T B2^T    )  ( xi)     ( 0 )
!    !
!    ! with A = L + nonlinear Convection. We compute in advance
!    ! a standard Stokes matrix L which can be added later to the
!    ! convection matrix, resulting in the nonlinear system matrix,
!    ! as well as both B-matrices.
!
!    ! Assemble the Stokes operator:
!    call stdop_assembleLaplaceMatrix (rstaticInfo%rmatrixStokes,.true.,&
!        rproblem%rphysicsPrimal%dnu)
!
!    ! Build the first pressure matrix B1.
!    call stdop_assembleSimpleMatrix (rstaticInfo%rmatrixB1,&
!        DER_FUNC,DER_DERIV_X,-1.0_DP)
!
!    ! Build the second pressure matrix B2.
!    call stdop_assembleSimpleMatrix (rstaticInfo%rmatrixB2,&
!        DER_FUNC,DER_DERIV_Y,-1.0_DP)
!
!    ! Set up the matrices D1 and D2 by transposing B1 and B2.
!    ! For that purpose, virtually transpose B1/B2.
!    call lsyssc_transposeMatrix (rstaticInfo%rmatrixB1,&
!        rstaticInfo%rmatrixD1,LSYSSC_TR_CONTENT)
!
!    call lsyssc_transposeMatrix (rstaticInfo%rmatrixB2,&
!        rstaticInfo%rmatrixD2,LSYSSC_TR_CONTENT)
!
!    ! -----------------------------------------------------------------------
!    ! EOJ matrix, if the edge oriented stabilisation is active
!    ! -----------------------------------------------------------------------
!
!    call parlst_getvalue_int (rproblem%rparamList, 'CC-DISCRETISATION', &
!                              'IUPWIND1', iupwind1)
!    call parlst_getvalue_int (rproblem%rparamList, 'CC-DISCRETISATION', &
!                              'IUPWIND2', iupwind2)
!    call parlst_getvalue_double (rproblem%rparamList, 'CC-DISCRETISATION', &
!                              'DUPSAM1', dupsam1)
!    call parlst_getvalue_double (rproblem%rparamList, 'CC-DISCRETISATION', &
!                              'DUPSAM2', dupsam2)
!
!    if ((iupwind1 .eq. CCMASM_STAB_EDGEORIENTED3)) then
!      ! Create an empty matrix
!      call lsyssc_duplicateMatrix (rstaticInfo%rmatrixTemplateFEM,&
!          rstaticInfo%rmatrixEOJ1,LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!      call lsyssc_clearMatrix (rstaticInfo%rmatrixEOJ1)
!
!      ! Set up the jump stabilisation structure.
!      ! There's not much to do, only initialise the viscosity...
!      rjumpStabil%dnu = rproblem%rphysicsPrimal%dnu
!
!      ! Set stabilisation parameter
!      rjumpStabil%dgamma = dupsam1
!
!      ! Matrix weight
!      rjumpStabil%dtheta = 1.0_DP
!
!      ! Call the jump stabilisation technique to stabilise that stuff.
!      ! We can assemble the jump part any time as it's independent of any
!      ! convective parts...
!      call conv_jumpStabilisation2d (&
!                          rjumpStabil, CONV_MODMATRIX, &
!                          rstaticInfo%rmatrixEOJ1)
!    end if
!
!    if ((iupwind2 .eq. CCMASM_STAB_EDGEORIENTED3)) then
!
!      ! Perhaps the second matrix is like the first, then we don't
!      ! have to assemble it.
!      if ((dupsam1 .eq. dupsam2) .and. (iupwind1 .eq. iupwind2)) then
!        ! Take that from the primal space.
!        call lsyssc_duplicateMatrix (rstaticInfo%rmatrixEOJ1,rstaticInfo%rmatrixEOJ2,&
!            LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
!      else
!        ! Create an empty matrix
!        call lsyssc_duplicateMatrix (rstaticInfo%rmatrixTemplateFEM,&
!            rstaticInfo%rmatrixEOJ2,LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
!        call lsyssc_clearMatrix (rstaticInfo%rmatrixEOJ2)
!
!        ! Set up the jump stabilisation structure.
!        ! There's not much to do, only initialise the viscosity...
!        rjumpStabil%dnu = rproblem%rphysicsPrimal%dnu
!
!        ! Set stabilisation parameter
!        rjumpStabil%dgamma = dupsam2
!
!        ! Matrix weight
!        rjumpStabil%dtheta = 1.0_DP
!
!        ! Call the jump stabilisation technique to stabilise that stuff.
!        ! We can assemble the jump part any time as it's independent of any
!        ! convective parts...
!        call conv_jumpStabilisation2d (&
!                            rjumpStabil, CONV_MODMATRIX, &
!                            rstaticInfo%rmatrixEOJ2)
!      end if
!    end if
!
!    ! -----------------------------------------------------------------------
!    ! Mass matrices. They are used in so many cases, it's better we always
!    ! have them available.
!    ! -----------------------------------------------------------------------
!
!    ! If there is an existing mass matrix, release it.
!    call lsyssc_releaseMatrix (rstaticInfo%rmatrixMass)
!    call lsyssc_releaseMatrix (rstaticInfo%rmatrixMassPressure)
!
!    ! Generate mass matrix. The matrix has basically the same structure as
!    ! our template FEM matrix, so we can take that.
!    call lsyssc_duplicateMatrix (rstaticInfo%rmatrixTemplateFEM,&
!                rstaticInfo%rmatrixMass,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
!    call lsyssc_duplicateMatrix (rstaticInfo%rmatrixTemplateFEMPressure,&
!                rstaticInfo%rmatrixMassPressure,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
!
!    ! Change the discretisation structure of the mass matrix to the
!    ! correct one; at the moment it points to the discretisation structure
!    ! of the Stokes matrix...
!    call lsyssc_assignDiscrDirectMat (rstaticInfo%rmatrixMass,&
!        rstaticInfo%rdiscretisationMass)
!    call lsyssc_assignDiscrDirectMat (rstaticInfo%rmatrixMassPressure,&
!        rstaticInfo%rdiscretisationMassPressure)
!
!    ! Call the standard matrix setup routine to build the matrix.
!    call stdop_assembleSimpleMatrix (rstaticInfo%rmatrixMass,DER_FUNC,DER_FUNC)
!    call stdop_assembleSimpleMatrix (rstaticInfo%rmatrixMassPressure,DER_FUNC,DER_FUNC)
!
!    ! Clean up the collection (as we are done with the assembly, that's it.
!    call user_doneCollectForAssembly (rproblem,rproblem%rcollection)
!
!  end subroutine
!
!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine cc_releaseStaticMatrices (rstaticInfo)
!
!!<description>
!  ! Releases all static matrices from a t_staticLevelInfo structure.
!!</description>
!
!!<inputoutput>
!  ! A t_staticLevelInfo structure to be cleaned up.
!  type(t_staticLevelInfo), intent(inout), target :: rstaticInfo
!!</inputoutput>
!
!!</subroutine>
!
!    ! If there is an existing mass matrix, release it.
!    call lsyssc_releaseMatrix (rstaticInfo%rmatrixMass)
!    call lsyssc_releaseMatrix (rstaticInfo%rmatrixMassPressure)
!
!    ! Release Stokes, B1, B2,... matrices
!    call lsyssc_releaseMatrix (rstaticInfo%rmatrixD2T)
!    call lsyssc_releaseMatrix (rstaticInfo%rmatrixD1T)
!    call lsyssc_releaseMatrix (rstaticInfo%rmatrixD2)
!    call lsyssc_releaseMatrix (rstaticInfo%rmatrixD1)
!    call lsyssc_releaseMatrix (rstaticInfo%rmatrixB2)
!    call lsyssc_releaseMatrix (rstaticInfo%rmatrixB1)
!    call lsyssc_releaseMatrix (rstaticInfo%rmatrixStokes)
!    if (lsyssc_hasMatrixStructure(rstaticInfo%rmatrixEOJ2)) then
!      call lsyssc_releaseMatrix (rstaticInfo%rmatrixEOJ2)
!    end if
!    if (lsyssc_hasMatrixStructure(rstaticInfo%rmatrixEOJ1)) then
!      call lsyssc_releaseMatrix (rstaticInfo%rmatrixEOJ1)
!    end if
!
!    ! Release the template matrices. This is the point, where the
!    ! memory of the matrix structure is released.
!    call lsyssc_releaseMatrix (rstaticInfo%rmatrixTemplateDivergence)
!    call lsyssc_releaseMatrix (rstaticInfo%rmatrixTemplateGradient)
!    call lsyssc_releaseMatrix (rstaticInfo%rmatrixTemplateFEM)
!    call lsyssc_releaseMatrix (rstaticInfo%rmatrixTemplateFEMPressure)
!
!  end subroutine
!
!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine cc_generateBasicRHS (rproblem,dtime,rrhs)
!
!!<description>
!  ! Calculates the entries of the basic right-hand-side vector on the finest
!  ! level. Boundary conditions or similar things are not implemented into
!  ! the vector.
!  ! Memory for the RHS vector must have been allocated in advance.
!!</description>
!
!!<input>
!  ! Current simulation time.
!  real(DP), intent(in) :: dtime
!!</input>
!
!!<inputoutput>
!  ! A problem structure saving problem-dependent information.
!  type(t_problem), intent(INOUT), target :: rproblem
!
!  ! The RHS vector which is to be filled with data.
!  type(t_vectorBlock), intent(INOUT) :: rrhs
!!</inputoutput>
!
!!</subroutine>
!
!  ! local variables
!
!    ! A bilinear and linear form describing the analytic problem to solve
!    type(t_linearForm) :: rlinform
!
!    ! A pointer to the discretisation structure with the data.
!    type(t_blockDiscretisation), pointer :: p_rdiscretisation
!
!    ! DEBUG!!!
!    real(DP), dimension(:), pointer :: p_Drhs
!
!    ! DEBUG!!!
!    call lsysbl_getbase_double (rrhs,p_Drhs)
!
!    ! Get a pointer to the RHS on the finest level as well as to the
!    ! block discretisation structure:
!    p_rdiscretisation => rrhs%p_rspaceDiscr
!
!    ! The vector structure is already prepared, but the entries are missing.
!    !
!    ! At first set up the corresponding linear form (f,Phi_j):
!    rlinform%itermCount = 1
!    rlinform%Idescriptors(1) = DER_FUNC
!
!    ! ... and then discretise the RHS to the first subvector of
!    ! the block vector using the discretisation structure of the
!    ! first block.
!    !
!    ! We pass our collection structure as well to this routine,
!    ! so the callback routine has access to everything what is
!    ! in the collection.
!    !
!    ! Note that the vector is unsorted after calling this routine!
!    !
!    ! Initialise the collection for the assembly process with callback routines.
!    ! Basically, this stores the simulation time in the collection if the
!    ! simulation is nonstationary.
!    call user_initCollectForAssembly (rproblem,dtime,rproblem%rcollection)
!
!    ! Discretise the X-velocity part:
!    call linf_buildVectorScalar (&
!              p_rdiscretisation%RspatialDiscr(1),rlinform,.true.,&
!              rrhs%RvectorBlock(1),user_coeff_RHS_x,&
!              rproblem%rcollection)
!
!    ! And the Y-velocity part:
!    call linf_buildVectorScalar (&
!              p_rdiscretisation%RspatialDiscr(2),rlinform,.true.,&
!              rrhs%RvectorBlock(2),user_coeff_RHS_y,&
!              rproblem%rcollection)
!
!    ! The third subvector must be zero initially - as it represents the RHS of
!    ! the equation "div(u) = 0".
!    call lsyssc_clearVector(rrhs%RvectorBlock(3))
!
!    ! The RHS terms for the dual equation are calculated similarly using
!    ! the desired 'target' flow field.
!    !
!    ! Discretise the X-velocity part:
!    call linf_buildVectorScalar (&
!              p_rdiscretisation%RspatialDiscr(1),rlinform,.true.,&
!              rrhs%RvectorBlock(4),user_coeff_TARGET_x,&
!              rproblem%rcollection)
!
!    ! And the Y-velocity part:
!    call linf_buildVectorScalar (&
!              p_rdiscretisation%RspatialDiscr(2),rlinform,.true.,&
!              rrhs%RvectorBlock(5),user_coeff_TARGET_y,&
!              rproblem%rcollection)
!
!    ! Depending on the formulation, to get a reference dual velocity,
!    ! it might be necessary to switch the sign of the target velocity field
!    ! because the RHS of the dual equation is '-z'!
!    ! Remember that it this case the signs of the mass matrices that couple
!    ! primal and dual velocity must be changed, too!
!
!    if (rproblem%roptcontrol%ispaceTimeFormulation .eq. 0) then
!      call lsyssc_scaleVector (rrhs%RvectorBlock(4),-1.0_DP)
!      call lsyssc_scaleVector (rrhs%RvectorBlock(5),-1.0_DP)
!    end if
!
!    ! Dual pressure RHS is =0.
!    call lsyssc_clearVector(rrhs%RvectorBlock(6))
!
!    ! Clean up the collection (as we are done with the assembly, that's it.
!    call user_doneCollectForAssembly (rproblem,rproblem%rcollection)
!
!  end subroutine
!
!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine ffunction_initSol (cderivative,rdiscretisation, &
!                nelements,npointsPerElement,Dpoints, &
!                IdofsTest,rdomainIntSubset,&
!                Dvalues,rcollection)
!
!  use basicgeometry
!  use triangulation
!  use collection
!  use scalarpde
!  use domainintegration
!
!!<description>
!  ! Callback routine for the calculation of the initial solution vector
!  ! from an analytical expression.
!!</description>
!
!!<input>
!  ! This is a DER_xxxx derivative identifier (from derivative.f90) that
!  ! specifies what to compute: DER_FUNC=function value, DER_DERIV_X=x-derivative,...
!  ! The result must be written to the Dvalue-array below.
!  integer, intent(IN)                                         :: cderivative
!
!  ! The discretisation structure that defines the basic shape of the
!  ! triangulation with references to the underlying triangulation,
!  ! analytic boundary boundary description etc.
!  type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
!
!  ! Number of elements, where the coefficients must be computed.
!  integer, intent(IN)                                         :: nelements
!
!  ! Number of points per element, where the coefficients must be computed
!  integer, intent(IN)                                         :: npointsPerElement
!
!  ! This is an array of all points on all the elements where coefficients
!  ! are needed.
!  ! DIMENSION(NDIM2D,npointsPerElement,nelements)
!  ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
!  real(DP), dimension(:,:,:), intent(IN)  :: Dpoints
!
!  ! An array accepting the DOF's on all elements trial in the trial space.
!  ! DIMENSION(\#local DOF's in trial space,Number of elements)
!  integer, dimension(:,:), intent(IN) :: IdofsTest
!
!  ! This is a t_domainIntSubset structure specifying more detailed information
!  ! about the element set that is currently being integrated.
!  ! It's usually used in more complex situations (e.g. nonlinear matrices).
!  type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset
!
!  ! Optional: A collection structure to provide additional
!  ! information to the coefficient routine.
!  type(t_collection), intent(INOUT), optional      :: rcollection
!
!!</input>
!
!!<output>
!  ! This array has to receive the values of the (analytical) function
!  ! in all the points specified in Dpoints, or the appropriate derivative
!  ! of the function, respectively, according to cderivative.
!  !   DIMENSION(npointsPerElement,nelements)
!  real(DP), dimension(:,:), intent(OUT)                      :: Dvalues
!!</output>
!
!!</subroutine>
!
!    real(DP) :: dtime
!    integer :: i,j,icomponent
!
!    real(DP), dimension(:), allocatable :: DvaluesAct
!
!    type(t_fparser), pointer :: p_rparser
!    real(dp), dimension(:,:), allocatable :: p_Dval
!
!    ! Get the parser object with the RHS expressions from the collection
!    p_rparser => collct_getvalue_pars (rcollection, 'SOLPARSER')
!
!    ! Current time
!    dtime = rcollection%DquickAccess(1)
!
!    ! X/Y component
!    icomponent = rcollection%IquickAccess(1)
!
!    ! Prepare the array with the values for the function.
!    ! X-coordinate, Y-coordinate, time.
!    allocate(p_Dval(3,npointsPerElement*nelements))
!    do i=1,nelements
!      do j=1,npointsPerElement
!        p_Dval (1,(i-1)*npointsPerElement+j) = Dpoints(1,j,i)
!        p_Dval (2,(i-1)*npointsPerElement+j) = Dpoints(2,j,i)
!        p_Dval (3,(i-1)*npointsPerElement+j) = dtime
!      end do
!    end do
!
!    ! Evaluate the 1st expression for the X-rhs
!    allocate(DvaluesAct(npointsPerElement*nelements))
!    call fparser_evalFunction (p_rparser, icomponent, 2, p_Dval, DvaluesAct)
!
!    ! Reshape the data, that's it.
!    do i=0,nelements-1
!      do j=1,npointsPerElement
!        Dvalues(j,i+1) = DvaluesAct(i*npointsPerElement+j)
!      end do
!    end do
!
!    deallocate(DvaluesAct)
!    deallocate(p_Dval)
!
!  end subroutine
!
!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine cc_initRHS (rproblem,rrhs)
!
!!<description>
!  ! Initialises the basic RHS based on the information in the DAT file.
!!</description>
!
!!<inputoutput>
!  ! A problem structure saving problem-dependent information.
!  type(t_problem), intent(INOUT) :: rproblem
!
!  ! An analytic solution structure receiving the abstract RHS.
!  type(t_anSolution) :: rrhs
!!</inputoutput>
!
!!</subroutine>
!
!    integer :: irhs
!    character(len=SYS_STRLEN) :: sstring
!    character(len=SYS_STRLEN), dimension(2) :: SrhsExpressions
!
!    ! Initialise a basic RHS
!    call ansol_init_meshless (rrhs)
!
!    ! Specification of the RHS
!    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
!                              'iRHS',irhs,0)
!
!    ! If the RHS is given as expression, create a parser object for the
!    ! expression.
!    if (irhs .eq. -1) then
!
!      ! Create an analytic solution structure based on the
!      ! RHS expressions.
!      call parlst_getvalue_string (rproblem%rparamList,'CC-DISCRETISATION',&
!                                  'srhsExpressionX',sstring,'''''')
!      read(sstring,*) SrhsExpressions(1)
!
!      call parlst_getvalue_string (rproblem%rparamList,'CC-DISCRETISATION',&
!                                  'srhsExpressionY',sstring,'''''')
!      read(sstring,*) SrhsExpressions(2)
!
!      call ansol_configAnalytical (rrhs,SrhsExpressions)
!
!    else if (irhs .eq. 0) then
!
!      ! Initialise a zero flow.
!      call ansol_configAnalytical (rrhs)
!
!    end if
!
!    ! The id int the RHS structure is irhs.
!    rrhs%iid = irhs
!
!  end subroutine
!
!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine cc_doneRHS (rrhs)
!
!!<description>
!  ! Cleans up the RHS.
!!</description>
!
!!<inputoutput>
!  ! An analytic solution structure receiving the abstract RHS.
!  type(t_anSolution) :: rrhs
!!</inputoutput>
!
!!</subroutine>
!
!    ! Release the RHS in case it was initialised.
!    call ansol_done(rrhs)
!
!  end subroutine
!
!  ! ***************************************************************************
!
!!<subroutine>
!
!  subroutine cc_writeSolution (rproblem,rvector)
!
!!<description>
!  ! Writes a solution vector rvector to a file as configured in the parameters
!  ! in the DAT file.
!!</description>
!
!!<input>
!  ! A problem structure saving problem-dependent information.
!  type(t_problem), intent(IN) :: rproblem
!
!  ! The solution vector to be written out. Must be set up according to the
!  ! maximum level NLMAX in rproblem!
!  type(t_vectorBlock), intent(IN) :: rvector
!!</input>
!
!!</subroutine>
!
!    ! local variables
!    integer(I32) :: idestLevel
!    type(t_vectorBlock) :: rvector1,rvector2
!    type(t_vectorScalar) :: rvectorTemp
!    character(LEN=SYS_STRLEN) :: sfile,sfileString
!    integer :: ilev
!    integer :: NEQ
!    type(t_interlevelProjectionBlock) :: rprojection
!    logical :: bformatted
!
!    ! Get the parameter what to do with rvector
!    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
!                              'iSolutionWrite',idestLevel,0)
!    call parlst_getvalue_string (rproblem%rparamList,'CC-DISCRETISATION',&
!                                 'sSolutionWrite',sfileString,'')
!
!    if (idestLevel .eq. 0) return ! nothing to do.
!
!    ! Remove possible ''-characters
!    read(sfileString,*) sfile
!
!    bformatted = idestLevel .gt. 0
!    idestLevel = rproblem%NLMAX-abs(idestLevel)+1 ! level where to write out
!
!    if (idestLevel .lt. rproblem%NLMIN) then
!      print *,'Warning: Level for solution vector is < NLMIN! &
!              &Writing out at level NLMIN!'
!      idestLevel = rproblem%NLMIN
!    end if
!
!    ! Interpolate the solution down to level istart.
!    call lsysbl_copyVector (rvector,rvector1)   ! creates new rvector1!
!
!    do ilev = rproblem%NLMAX,idestLevel+1,-1
!
!      ! Initialise a vector for the lower level and a prolongation structure.
!      call lsysbl_createVectorBlock (&
!          rproblem%RlevelInfo(ilev-1)%rdiscretisation,rvector2,.false.)
!
!      call mlprj_initProjectionVec (rprojection,rvector2)
!
!      ! Interpolate to the next higher level.
!      ! (Don't 'restrict'! Restriction would be for the dual space = RHS vectors!)
!
!      NEQ = mlprj_getTempMemoryVec (rprojection,rvector2,rvector1)
!      if (NEQ .ne. 0) call lsyssc_createVector (rvectorTemp,NEQ,.false.)
!      call mlprj_performInterpolation (rprojection,rvector2,rvector1, &
!                                       rvectorTemp)
!      if (NEQ .ne. 0) call lsyssc_releaseVector (rvectorTemp)
!
!      ! Swap rvector1 and rvector2. Release the fine grid vector.
!      call lsysbl_swapVectors (rvector1,rvector2)
!      call lsysbl_releaseVector (rvector2)
!
!      call mlprj_doneProjection (rprojection)
!
!    end do
!
!    ! Write out the solution.
!    if (bformatted) then
!      call vecio_writeBlockVectorHR (rvector1, 'SOLUTION', .true.,&
!                                    0, sfile, '(E22.15)')
!    else
!      call vecio_writeBlockVectorHR (rvector1, 'SOLUTION', .true.,0, sfile)
!    end if
!
!    ! Release temp memory.
!    call lsysbl_releaseVector (rvector1)
!
!  end subroutine

end module
