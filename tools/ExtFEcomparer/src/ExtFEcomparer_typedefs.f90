module ExtFEcomparer_typedefs

  use fsystem
  use storage

  use boundary
  use cubature

  use triangulation
  use spatialdiscretisation
  use paramlist
  use linearsystemblock

  use collection

  implicit none


 !<constant block>
 ! This block defines all constants used in the ExtFEcomparer
 ! Used for dimension/type of input/...

 ! Dimension
  integer, parameter, public :: ExtFE_NDIM1 = 1
  integer, parameter, public :: ExtFE_NDIM2 = 2
  integer, parameter, public :: ExtFE_NDIM3 = 3

  ! Formatted/unformatted vector file?
  integer, parameter, public :: ExtFE_formatted_unsorted = 1
  integer, parameter, public :: ExtFE_unformatted_unsorted = 2
  integer, parameter, public :: ExtFE_formatted_sorted = 3
  integer, parameter, public :: ExtFE_unformatted_sorted = 4

  integer, parameter, public :: ExtFE_formatted = 5
  integer, parameter, public :: ExtFE_unformatted = 6

  ! "Vector type": like i.e. cc2d (u,v,p) where we read in the
  ! solution of the system without "postprocessing"
  ! or flagship-style - which calculates in the conservative variables
  ! but can write out variables that are not-conservative
  integer, parameter, public :: ExtFE_isPostprocessed = 1
  integer, parameter, public :: ExtFE_isSolutionBlockvector = 2

  ! What kind of elements did you use?
  ! An element pair (like in cc2d) -
  ! where you have one element for the speed
  ! and a different one for the pressure
  ! or did you use one element for everything?
  integer, parameter, public :: ExtFE_OneElement = 1
  integer, parameter, public :: ExtFE_ElementPair = 2


  ! For the sake of readability
  integer, parameter, public :: ExtFE_DO = 1
  integer, parameter, public :: ExtFE_DONT = 0

  integer, parameter, public :: ExtFE_calcPointValue = 1
  integer, parameter, public :: ExtFE_calcL2Norm = 2


  integer, parameter, public :: ExtFE_DER_FUNC_1D = 0
  integer, parameter, public :: ExtFE_DER_DERIV_1D_X = 1
  integer, parameter, public :: ExtFE_DER_FUNC_2D = 0
  integer, parameter, public :: ExtFE_DER_DERIV_2D_X = 1
  integer, parameter, public :: ExtFE_DER_DERIV_2D_Y = 2
  integer, parameter, public :: ExtFE_DER_FUNC_3D = 0
  integer, parameter, public :: ExtFE_DER_DERIV_3D_X = 1
  integer, parameter, public :: ExtFE_DER_DERIV_3D_Y = 2
  integer, parameter, public :: ExtFE_DER_DERIV_3D_Z = 3

  ! We have 3 derivatives + the function itself, so
  integer, parameter, public :: ExtFE_NumDerivs = 4

  ! Some day (if used in 3d) we might find out
  ! that the SYS_STRLEN of 256 is not long enough
  ! to read in expressions for ChiOmega
  ! and store them.
  ! So we define an own constant here. At the moment
  ! we keep it at the SYS_STRLEN, but if we find out
  ! that this is not long enough we only need to change
  ! this constant here and still all dependencies
  ! and arithmetics work.
  ! It might become relevant in the output/postprocessing
  ! routines at first
  integer, parameter, public :: ExtFE_STRLEN = SYS_STRLEN

  !</constant Block>



  type t_problem

    ! We need a parameterlist for sure
   type(t_parlist) :: rparamlist

   integer:: NLMIN, NLMAX

   ! We need the triangulation
   type(t_triangulation) :: rtriangulation

   ! and the discretisation
   type(t_blockDiscretisation) :: rdiscretisation

   ! The Type of elements
   integer:: ielemtype

   ! An object for saving the domain:
   type(t_boundary) :: rboundary

   !We also save the vector in here
   type(t_vectorBlock) :: coeffVector

    ! A collection object that saves structural data and some
    ! problem-dependent information which is e.g. passed to
    ! callback routines.
    type(t_collection) :: rcollection

    ! We need the path of the vector and the parametrisation
    character(LEN=ExtFE_STRLEN) :: sPRMFile, sTRIFile
    character(LEN=ExtFE_STRLEN) :: sVectorFile

    ! Dimension of the problem
    integer :: iDimension

    ! Formatted/unformatted/sorted/unsorted?
    integer :: vectorFileFormat

    ! i.e. Speed/Pressure from cc2d
    ! or one  containing only 1 variable
    integer :: vectorType

    ! Element pair or one Element?
    integer :: elementSetting

    ! How many Variables are in the Vector?
    ! This is not the number of equations - in cc2d
    ! it would be 3 as we have u,v and p as variables!
    integer :: NVAR

  end type



  type t_postprocessing
    ! This structure stores every information we get
    ! to hand it over to the output-routine
    ! it might be overkill at the moment,
    ! but it gives more structure to the code
    ! and allows to easily extend the code

    ! Everything regarding L2-Output
    logical :: writeOutL2results = .false.
    character(LEN=ExtFE_STRLEN) :: L2filepath

    ! Results + Which components where involved
    ! All handles shall become arrays:
    ! The result array in double, the other in integer
    integer :: h_L2Results = ST_NOHANDLE
    integer :: h_L2CompFunc = ST_NOHANDLE
    ! What domain was used?
    integer :: h_L2ChiOmega = ST_NOHANDLE
    integer :: h_L2TriFile = ST_NOHANDLE


    ! Everything regarding the pointvalues
    logical :: writeOutPointCalucations = .false.
    character(LEN=ExtFE_STRLEN) :: PointFilepath
    ! PointResults will become a 1D-array (DP)
    ! PointCoordinates will become a 2D-Array (DP)
    ! PointFuncComponents will become a 2D-Array(int)
    ! with the format (component1 derivative1 comp2 deriv2)
    integer :: h_PointResults = ST_NOHANDLE
    integer :: h_PointCoordinates = ST_NOHANDLE
    integer :: h_PointFuncComponents = ST_NOHANDLE



  end type

end module
