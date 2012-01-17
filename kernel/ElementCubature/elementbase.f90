!##############################################################################
!# ****************************************************************************
!# <name> elementbase </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains a set of constant definitions and structures which
!# are used by the different element modules.
!#
!# </purpose>
!##############################################################################

module elementbase

!$use omp_lib
  use basicgeometry
  use derivatives
  use fsystem
  use perfconfig
  use transformation
  use triangulation

  implicit none

  private
  
!<constants>

!<constantblock description="Element evaluation tags. Defines the basis information that is \
!  necessary to evaluate an element. All tags define bits of a bitfield and can be \
!  combined by an OR command.">

  ! Calculate the coordinates of teh points that form the element.
  ! (Usually the corners of the element.)
  integer(I32), parameter, public :: EL_EVLTAG_COORDS       = 2**0

  ! Prepare the coordinates on the reference element of the points where to evaluate.
  ! A set of coordinates on one element must be specified in the preparation routine.
  ! These coordinates are transferred to all elements where to evaluate.
  ! For the preparation routine on one element: Transfer the source coordinates
  ! on the reference element to the evaluation structure.
  integer(I32), parameter, public :: EL_EVLTAG_REFPOINTS    = 2**1

  ! Calculate the real coordinates of the points where to evaluate.
  ! The real coordinates are calculated based on the corresponding mapping from the
  ! reference to the real element.
  integer(I32), parameter, public :: EL_EVLTAG_REALPOINTS   = 2**2
  
  ! Calculate the Jacobian matrix of the points where to evaluate
  integer(I32), parameter, public :: EL_EVLTAG_JAC          = 2**3

  ! Calculate the Jacobian determinant in the points where to evaluate
  integer(I32), parameter, public :: EL_EVLTAG_DETJ         = 2**4
  
  ! Calculate the twist indices
  integer(I32), parameter, public :: EL_EVLTAG_TWISTIDX     = 2**5
  
!</constantblock>

!</constants>


!<types>

!<typeblock>

  ! This structure collects information that is necessary by an element
  ! to be evaluated on a single point on a single cell.
  type t_evalElement
  
    ! Array with coordinates of the corners that form the real element.
    ! DIMENSION(#space dimensions)
    !  Dcoords(1)=x-coordinate,
    !  Dcoords(2)=y-coordinate.
    !  Dcoords(3)=z-coordinate (only in 3D).
    real(DP), dimension(NDIM3D,TRIA_MAXNVE) :: Dcoords

    ! Values of the Jacobian matrix that defines the mapping between the
    ! reference element and the real elements. For every point i (e.g. in 2D):
    !  Djac(1) = J_i(1,1)
    !  Djac(2) = J_i(2,1)
    !  Djac(3) = J_i(1,2)
    !  Djac(4) = J_i(2,2)
    ! ...
    ! Remark: Only used for calculating derivatives; can be set to 0.0
    ! when derivatives are not used.
    real(DP), dimension(NDIM3D*NDIM3D) :: Djac
    
    ! Determinant of the mapping from the reference element to the real
    ! elements in the point to be evaluated.
    ! Remark: Only used for calculating derivatives; can be set to 1.0
    ! when derivatives are not needed. Must not be set to 0.0!
    real(DP) :: ddetj
    
    ! Array with coordinates of the points where to evaluate.
    ! The coordinates are expected on the reference element in the corresponding
    ! coordinate system. E.g. for real coordinates, this means:
    !  Dpoints(1)=x-coordinates,
    !  Dpoints(2)=y-coordinates,
    !  Dpoints(3)=z-coordinates (in 3D).
    ! but the entries here can also be e.g. in barycentric coordinates,
    ! if the coordinate system on the reference element prescribes that.
    real(DP), dimension(TRAFO_MAXDIMREFCOORD) :: DpointRef

    ! Array with coordinates of the points where to evaluate.
    ! The coordinates are expected on the real element.
    ! It is assumed that:
    !  Dpoints(1)=x-coordinates,
    !  Dpoints(2)=y-coordinates.
    !  Dpoints(3)=z-coordinates (only 3d).
    real(DP), dimension(NDIM3D) :: DpointReal
    
    ! Twist index entry of the element.
    integer(I32) :: itwistIndex

  end type

  public :: t_evalElement
!</typeblock>

!<typeblock>

  ! This structure collects information that is necessary by an element
  ! to be evaluated on a set of points on a set of element cells.
  type t_evalElementSet
  
    ! Number of points on every element where to evalate the basis functions.
    integer :: npointsPerElement = 0
    
    ! Number of elements, this element set consists of
    integer  :: nelements = 0
    
    ! Array with coordinates of the corners that form the real element.
    ! DIMENSION(#space dimensions,NVE,nelements)
    !  Dcoords(1,.,.)=x-coordinates,
    !  Dcoords(2,.,.)=y-coordinates.
    !  Dcoords(3,.,.)=z-coordinates (only in 3D).
    ! furthermore:
    !  Dcoords(:,i,.) = Coordinates of vertex i
    ! furthermore:
    !  Dcoords(:,:,j) = Coordinates of all corner vertices of element j
    real(DP), dimension(:,:,:), pointer :: p_Dcoords => null()

    ! Values of the Jacobian matrix that defines the mapping between the
    ! reference element and the real elements. For every point i (e.g. in 2D):
    !  Djac(1,i,.) = J_i(1,1,.)
    !  Djac(2,i,.) = J_i(2,1,.)
    !  Djac(3,i,.) = J_i(1,2,.)
    !  Djac(4,i,.) = J_i(2,2,.)
    ! Remark: Only used for calculating derivatives; can be set to 0.0
    ! when derivatives are not used.
    !  Djac(:,:,j) refers to the determinants of the points of element j.
    real(DP), dimension(:,:,:), pointer :: p_Djac => null()
    
    ! Determinant of the mapping from the reference element to the real
    ! elements for every of the npointsPerElement points on all the elements.
    !  Ddetj(i,.) = Determinant of point i
    !  Ddetj(:,j) = determinants of all points on element j
    ! Remark: Only used for calculating derivatives; can be set to 1.0
    ! when derivatives are not needed. Must not be set to 0.0!
    real(DP), dimension(:,:), pointer :: p_Ddetj => null()
    
    ! Array with coordinates of the points where to evaluate.
    ! The coordinates are expected on the reference element.
    ! It is assumed that:
    !  Dpoints(1,.)=x-coordinates,
    !  Dpoints(2,.)=y-coordinates.
    !  Dpoints(3,.)=y-coordinates.
    ! furthermore:
    !  Dpoints(:,i,.) = Coordinates of point i
    ! furthermore:
    !  Dpoints(:,:,j) = Coordinates of all points on element j
    real(DP), dimension(:,:,:), pointer :: p_DpointsRef => null()

    ! Array with coordinates of the points where to evaluate.
    ! The coordinates are expected on the real element.
    ! It is assumed that:
    !  Dpoints(1,.)=x-coordinates,
    !  Dpoints(2,.)=y-coordinates.
    !  Dpoints(3,.)=z-coordinates (only 3d).
    ! furthermore:
    !  Dpoints(:,i,.) = Coordinates of point i
    ! furthermore:
    !  Dpoints(:,:,j) = Coordinates of all points on element j
    real(DP), dimension(:,:,:), pointer :: p_DpointsReal => null()
    
    ! Twist index array of the elements.
    ! Array with DIMENSION(nelements)
    integer(I32), dimension(:), pointer :: p_ItwistIndex => null()
  
    ! .TRUE., if the array with coordinates on the reference element
    ! is maintained by the caller; prevents release of memory in the
    ! cleanup routines.
    logical :: bforeignPointsRef = .false.
    
    ! .TRUE., if the array with coordinates on the real element
    ! is maintained by the caller; prevents release of memory in
    ! cleanup routines.
    logical :: bforeignPointsReal = .false.
  
    ! .TRUE., if the array p_Dcoords is maintained by the caller; prevents
    ! release of memory in cleanup routines.
    logical :: bforeignCoords = .false.
  
    ! Pointer to a performance pointer
    type(t_perfconfig), pointer :: p_rperfconfig => null()

  end type

  public :: t_evalElementSet

!</typeblock>

!</types>

end module
