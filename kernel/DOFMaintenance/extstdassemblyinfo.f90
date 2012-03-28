!##############################################################################
!# ****************************************************************************
!# <name> extstdassemblyinfo </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains additional assembly information that allows to
!# configure the assembly routines. The basic structure t_scalarCubatureInfo
!# controls e.g. the cubature rule used during the assembly.
!#
!# The t_stdCubatureData structure contains basic dynamic information
!# about a cubature formula. It can be initialised using a cubature
!# formula ID.
!#
!# The t_stdFEBasisEvalData structure contains basic information
!# for the evaluatiuon of finite element basis functions on a set of
!# points on a set of elements.
!#
!# Routines in this module:
!#
!# 1.) easminfo_initStdCubature
!#      -> Create a structure with all necessary information about
!#         a cubature formula.
!#
!# 2.) easminfo_doneStdCubature
!#     -> Release the structure with the cubature information.
!#
!# 3.) easminfo_initStdFEBasisEval
!#     -> Create an evaluation structure to evaluate a finite element
!#        on a set of points for a set of elements.
!#
!# 4.) easminfo_doneStdFEBasisEval
!#     -> Release an initialised evaluation structure.
!#
!#
!# </purpose>
!##############################################################################

module extstdassemblyinfo

!$use omp_lib
  use fsystem
  use genoutput
  use storage
  use basicgeometry
  use boundary
  use boundaryaux
  use cubature
  use scalarpde
  use element
  use transformation
  use spatialdiscretisation

  implicit none

  private

!<types>

!<typeblock>

  ! Standard information for cubature.
  type t_stdCubatureData

    ! ID of the cubature rule.
    integer(I32) :: ccubature = CUB_UNDEFINED

    ! Number of cubature points per element
    integer :: ncubp = 0

    ! Cubature weights
    real(DP), dimension(:), pointer :: p_Domega

    ! An array that takes coordinates of the cubature formula on the reference element.
    real(DP), dimension(:,:), pointer :: p_DcubPtsRef

  end type

!</typeblock>

!<typeblock>

  ! Standard information for evaluating FE basis functions in
  ! a set of cubature points on a set of elements.
  type t_stdFEBasisEvalData

    ! Element ID
    integer(I32) :: celement = EL_UNDEFINED    

    ! Number of local DOF`s
    integer :: ndofLocal = 0

    ! Number of vertices per element
    integer :: NVE = 0

    ! Highest supported derivative ID
    integer :: nmaxderivative = 0

    ! Number of points on each element
    integer :: npoints = 0

    ! Number of elements simultaneously supported
    ! by this structure
    integer :: nelements = 0

    ! Type of transformation
    integer(I32) :: ctrafoType = TRAFO_ID_UNKNOWN

    ! Arrays saving the DOF`s in the elements
    integer, dimension(:,:), pointer :: p_Idofs

    ! Arrays for the basis function values in the points.
    ! dimension(ndofLocal,nmaxderivative,npoints,nelements)
    real(DP), dimension(:,:,:,:), pointer :: p_Dbas

  end type

!</typeblock>

!</types>

!<constants>

!</constants>

  public :: t_stdCubatureData
  public :: easminfo_initStdCubature
  public :: easminfo_doneStdCubature

  public :: t_stdFEBasisEvalData
  public :: easminfo_initStdFEBasisEval
  public :: easminfo_doneStdFEBasisEval

contains

  ! ***************************************************************************

!<subroutine>

  subroutine easminfo_initStdCubature (ccubature,rcubatureInfo)

!<description>
  ! Initialises the standard cubature structure.
!</description>

!<input>
  ! ID of the cubature rule to be used.
  integer(I32), intent(in) :: ccubature
!</input>

!<inputoutput>
  ! Structure to be set up.
  type(t_stdCubatureData), intent(out) :: rcubatureInfo
!</inputoutput>

!</subroutine>

    ! Cubature rule ID.
    rcubatureInfo%ccubature = ccubature

    ! Get the number of cubature points for the cubature formula
    rcubatureInfo%ncubp = cub_igetNumPts(ccubature)

    ! Allocate two arrays for the points and the weights
    allocate(rcubatureInfo%p_Domega(rcubatureInfo%ncubp))
    allocate(rcubatureInfo%p_DcubPtsRef(&
        cub_igetCoordDim(ccubature),rcubatureInfo%ncubp))

    ! Get the cubature formula
    call cub_getCubature(ccubature,rcubatureInfo%p_DcubPtsRef,rcubatureInfo%p_Domega)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine easminfo_doneStdCubature (rcubatureInfo)

!<description>
  ! Cleans up a standard cubature structure.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_stdCubatureData), intent(inout) :: rcubatureInfo
!</inputoutput>

!</subroutine>

    ! Cubature rule ID.
    rcubatureInfo%ccubature = CUB_UNDEFINED

    ! Get the number of cubature points for the cubature formula
    rcubatureInfo%ncubp = 0

    ! Allocate two arrays for the points and the weights
    deallocate(rcubatureInfo%p_Domega)
    deallocate(rcubatureInfo%p_DcubPtsRef)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine easminfo_initStdFEBasisEval (celement,&
      nmaxderivative,npoints,nelements,rfeBasisEvalData)

!<description>
  ! Initialises the standard cubature structure.
!</description>

!<input>
  ! ID of the element to be used.
  integer(I32), intent(in) :: celement

  ! Highest derivative ID which will occur during the evaluation
  integer, intent(in) :: nmaxderivative

  ! Number of points on each element.
  integer, intent(in) :: npoints

  ! Maximum supported number of elements.
  integer, intent(in) :: nelements
!</input>

!<output>
  ! Structure to be set up.
  type(t_stdFEBasisEvalData), intent(out) :: rfeBasisEvalData
!</output>

!</subroutine>

    ! Element ID.
    rfeBasisEvalData%celement = celement

    ! Get the number of local DOF`s for trial and test functions
    rfeBasisEvalData%ndofLocal = elem_igetNDofLoc(celement)

    ! Get the number of vertices of the element, specifying the transformation
    ! form the reference to the real element.
    rfeBasisEvalData%NVE = elem_igetNVE(celement)

    ! Get from the element space the type of coordinate system
    ! that is used there:
    rfeBasisEvalData%ctrafoType = elem_igetTrafoType(celement)

    ! Number of points and number of elements.
    rfeBasisEvalData%npoints = npoints
    rfeBasisEvalData%nelements = nelements

    ! Maximum supported derivative
    rfeBasisEvalData%nmaxDerivative = nmaxderivative

    ! Allocate memory for the DOF`s of all the elements.
    allocate(rfeBasisEvalData%p_Idofs(rfeBasisEvalData%ndofLocal,nelements))

    ! Allocate arrays for the values of the test functions.
    allocate(rfeBasisEvalData%p_Dbas(rfeBasisEvalData%ndofLocal,&
             rfeBasisEvalData%nmaxDerivative,npoints,nelements))

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine easminfo_doneStdFEBasisEval (rfeBasisEvalData)

!<description>
  ! Cleans up a standard structure for the evaluation of FE basis functions.
!</description>

!<inputoutput>
  ! Structure to be cleaned up.
  type(t_stdFEBasisEvalData), intent(inout) :: rfeBasisEvalData
!</inputoutput>

!</subroutine>

    rfeBasisEvalData%celement = EL_UNDEFINED
    rfeBasisEvalData%ndofLocal = 0
    rfeBasisEvalData%NVE = 0
    rfeBasisEvalData%ctrafoType = TRAFO_ID_UNKNOWN
    rfeBasisEvalData%npoints = 0
    rfeBasisEvalData%nelements = 0
    rfeBasisEvalData%nmaxDerivative = 0

    deallocate(rfeBasisEvalData%p_Dbas)
    deallocate(rfeBasisEvalData%p_Idofs)

  end subroutine


end module
