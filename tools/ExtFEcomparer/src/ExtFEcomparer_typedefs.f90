module ExtFEcomparer_typedefs

  use fsystem
  use storage

  use boundary
  use cubature

  use triangulation
  use spatialdiscretisation
  use paramlist
  use linearsystemblock

  use fparser

  implicit none


 !<constant block>
 ! This block defines all constants used in the ExtFEcomparer
 ! Used for dimension/type of input/...

 ! Dimension
  integer, parameter, public :: ExtFE_NDIM1 = 1
  integer, parameter, public :: ExtFE_NDIM2 = 2
  integer, parameter, public :: ExtFE_NDIM3 = 3

  ! These are the possibilities to write out with
  ! vecio_writeBlockVectorHR, we want to read them in
  ! so we need the types.
  integer, parameter, public :: ExtFE_formatted_unsorted = 1
  integer, parameter, public :: ExtFE_unformatted_unsorted = 2
  integer, parameter, public :: ExtFE_formatted_sorted = 3
  integer, parameter, public :: ExtFE_unformatted_sorted = 4

  ! Another input way is without the routine. For that we need
  ! to know if the file is formatted or unformatted
  integer, parameter, public :: ExtFE_formatted = 5
  integer, parameter, public :: ExtFE_unformatted = 6

  ! "Vector type": like i.e. cc2d (u,v,p) where we read in the
  ! solution of the system without "postprocessing"
  ! or flagship-style - which calculates in the conservative variables
  ! but can write out variables that are not-conservative
  ! (or, in other words: ExtFE_isSolutionBlockvector means
  ! that it is written out with vecio_writeBlockVectorHR)
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

  ! We need to do point evaluations, but also in the
  ! Derivatives. These ones are to translate the
  ! input to the global constants.
  ! Attention: Function is always zero, x-deriv always 1,
  ! y-deriv always 2 and z-deriv always 3
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

  ! Constants for the UCD-Output
  ! UCD-Types supported
  integer, parameter, public :: ExtFE_UCD_VTK = 1

  ! which way to add a component?
  ! Vertex-Based?
  integer, parameter, public :: ExtFE_UCD_VERT_BASED = 1
  ! Element based?
  integer, parameter, public :: ExtFE_UCD_ELEM_BASED = 2

  ! Translation of the input of ucd-flags
  integer, parameter, public :: ExtFE_UCD_OUT_STANDARD = 1
  integer, parameter, public :: ExtFE_UCD_OUT_DISCONTINUOUS = 2

  ! To generate UCD-Output of the functions, we have to
  ! project the solution vector on a Q1 or Q0 space as
  ! most programs don't understand other formats.
  ! These constants are used during the input and calculation
  ! for this sake
  integer, parameter, public :: ExtFE_UCD_POLY_CONST = 0
  integer, parameter, public :: ExtFE_UCD_POLY_LINEAR = 1
  !</constant Block>


! This is our global "Problem-Structure". However, this name
! is not really good, better would be "wannabe-FEM-Function":
! It collects all information that is necessary to reconstruct
! a FEM-Function, and then the coefficent vector in the end
! will be a FEM-function

  type t_problem

    ! We need a parameterlist for sure
   type(t_parlist) :: rparamlist

   ! How often was the mesh refined?
   ! nRefinements would be a better name, but
   ! it is called NLMAX in most applications I saw so
   ! I kept the name
   integer:: NLMAX

   ! We need the triangulation
   type(t_triangulation) :: rtriangulation

   ! and the discretisation
   type(t_blockDiscretisation) :: rdiscretisation

   ! The Type of elements
   ! Either an ID telling which element pair/triple/... was used...
   integer:: iElemPairID
   ! Or the ID of an element that was used for all variables
   ! of the FEM-Function
   integer(I32) :: ielemType

   ! An object for saving the domain:
   type(t_boundary) :: rboundary

   !We also save the vector in here
   type(t_vectorBlock) :: coeffVector

   ! We need the path of the vector and the parametrisation
   character(LEN=ExtFE_STRLEN) :: sPRMFile, sTRIFile
   character(LEN=ExtFE_STRLEN) :: sVectorFile

   ! Dimension of the problem
   integer :: iDimension

   ! We need to know what kind of file we will get to be able
   ! to read it in the correct way
   ! One of these constants:
   ! ExtFE_formatted_unsorted, ExtFE_unformatted_unsorted,
   ! ExtFE_formatted_sorted,   ExtFE_unformatted_sorted,
   ! ExtFE_formatted,          ExtFE_unformatted
   integer :: vectorFileFormat

   ! We need more info to read it correct:
   ! Is it written out with vecio_writeBlockVectorHR
   ! or in the other format? One of these constants:
   ! ExtFE_isPostprocessed, ExtFE_isSolutionBlockvector
   integer :: vectorType

   ! Element pair or one Element?
   ! One of these constants:
   ! ExtFE_OneElement, ExtFE_ElementPair
   integer :: elementSetting

   ! How many Variables are in the Vector?
   ! This is not the number of equations - in cc2d
   ! it would be 3 as we have u,v and p as variables!
   ! Because we always create first the discretisation and
   ! then the vector according to it (as else we cannot
   ! guarantee that all variables in the vector are set correct)
   ! we need this information in advance
   integer :: NVAR

  end type


 ! This is our "postprocessing structure".
 ! It got the name because it saves all information what
 ! we actually want to do with our 2 FEM-Functions.
 ! During the init, it shall get everything that is
 ! possible, to be precise it shall get everything
 ! but the results.
 ! It is also handed over to the calculation routine
 ! as it has all information what we want to calculate.
 ! After that, it goes to the postprocessing routine
 ! where the results are postprocessed
 ! Postprocessed includes every bit of file i/o,
 ! so even writing out the original FEM-Functions in
 ! an other format is considered postprocessing and thus
 ! done in the postprocessing routines
 ! We make heavily use of the handle-system provided by
 ! the Featflow2 storage system. Storing an integer
 ! and "creating a pointer from it" where needed
 ! is easier and cheaper than sending the whole structure
 ! (with the according arrays instead of handles)
 ! through the code as it gets rather big

  ! Intial all values are those that most likely
  ! a) produce an error
  ! b) do not trigger any computation
  ! c) are empty
  type t_postprocessing

    ! General information we need everywhere:
    ! Dimension
    integer :: nDim = -1
    ! How many variables are in first/second function?
    integer :: nVarFirst = -1
    integer :: nVarSecond = -1

    ! Everything regarding vector output
    ! Do we have to write them out?
    logical :: writeOutOrigVector1 = .false.
    logical :: writeOutOrigVector2 = .false.
    ! If yes, where to write?
    character(LEN=ExtFE_STRLEN) :: sOrigVecPathOutFirst = ''
    character(LEN=ExtFE_STRLEN) :: sOrigVecPathOutSecond = ''
    ! Which format? Every fortran format specifier is allowed,
    ! but not every makes sense. Standard will be '(E22.15)'
    character(LEN=ExtFE_STRLEN) :: sOrigVec1OutFMT = ''
    character(LEN=ExtFE_STRLEN) :: sOrigVec2OutFMT = ''
    ! We need a pointer to point to the vectors we want
    ! to write out. Else we do not know what to write out
    type(t_vectorBlock), pointer :: OrigVecFirst => NULL()
    type(t_vectorBlock), pointer :: OrigVecSecond => NULL()


    ! Everything regarding L2-Output and calculations
    ! How many L2-Calculations?
    integer :: nL2Calculations = -1
    ! Write the results out in a file?
    logical :: writeOutL2results = .false.
    ! If we write them out in a file, we need a path
    ! where to save it
    character(LEN=ExtFE_STRLEN) :: L2filepath = ''
    ! h_L2Results shall save the handle to a
    ! an array that stores nL2Calculations doubles
    ! - entry i: result of calculation i
    integer :: h_L2Results = ST_NOHANDLE
    ! h_L2CompFunc shall store the handle to an
    ! array that saves 2*nL2Calculations INTs
    ! entry (1,i) is the component of function
    ! one in calculation i, (2,i) the component of
    ! function 2 in calculation i.
    ! if one comp. is -1 then this component is not used,
    ! i.e we calc ||f|| and not ||f-g||
    integer :: h_L2CompFunc = ST_NOHANDLE
    ! For the Output we want print on which domain we
    ! calculated. h_L2ChiOmega stores the handle to a
    ! char-array with the characteristic function of
    ! that area and h_L2TriFile stores a handle with the
    ! tri-File (aka mesh) that was used for the calculation
    ! Index aritmetic: entries (i-1)*ExtFE_STRLEN+1 to i*ExtFE_STRLEN
    ! belong to calculation i.
    integer :: h_L2ChiOmega = ST_NOHANDLE
    integer :: h_L2TriFile = ST_NOHANDLE
    ! For each integral we need a cubature rule. h_L2CubRule
    ! stores the handle to an INT32-array which stores in
    ! component i the ID of the cub-rule we use
    integer :: h_L2CubRule = ST_NOHANDLE
    ! We need to pass the characteristic function
    ! somehow to the callback routine that actually
    ! calculates the Integral. We could just parse in place,
    ! but for performance reasons we parse it in the init,
    ! in component i of the parser is the function for
    ! the calculation i.
    ! It still needs some init
    type(t_fparser) :: pL2ChiOmegaParser

    ! Everything regarding the pointvalues
    ! How many PointCalculations?
    integer :: nPointCalculations = -1
    ! Do we want to write them out to a file?
    logical :: writeOutPointCalucations = .false.
    ! If yes, we need a path for the file where to write
    ! it out.
    character(LEN=ExtFE_STRLEN) :: PointFilepath = ''
    ! In h_PointResults we store the handle to an array
    ! where we write out the results of the point
    ! evaluations. It is a double array, size: nPointCalculations.
    ! In entry i we have the result of the calculation i
    integer :: h_PointResults = ST_NOHANDLE
    ! We need the coordinates where we shall calculate.
    ! This int saves the handle to a 2D-double-array that stores
    ! these information.
    ! entry (1,i): x-component for calculation i
    ! entry (2,i): y-component for calculation i
    ! entry (3,i): z-component for calculation i
    integer :: h_PointCoordinates = ST_NOHANDLE
    ! We also need to know which components
    ! and which derivatives of each
    ! function to take. h_PointFuncComponents stores
    ! a handle to an 2D-Int-Array with this information:
    ! Format:
    ! (1,i): Which component of function 1 for calculation i
    ! (2,i): Which derivative of the component of function 1 in calc i
    ! (3,i): Which component of function 2 for calculation i
    ! (4,i): Which derivative of the component of function 2 in calc i
    ! If you want to evaluate only function 1, set the component of
    ! function 1 to -1. The derivative identifier must be present
    ! even if you do not want to evaluate something from the first function!
    ! For the second function goes the same.
    ! Supported derivative Identifiers are the constants from above!
    integer :: h_PointFuncComponents = ST_NOHANDLE

    ! We want to be able to do UCD-Output of the mesh
    ! and of the functions itself.
    ! Which UCD-Format? One of these constants:
    ! ExtFE_UCD_VTK
    integer :: UCD_Format
    ! Do we want to write out the meshes?
    logical :: ucd_OUT_meshes = .false.
    ! If yes, where to save them?
    character(LEN=ExtFE_STRLEN) :: UCD_meshOneOutPath = ''
    character(LEN=ExtFE_STRLEN) :: UCD_meshTwoOutPath = ''
    ! We actually need a pointer to the triangulations.
    ! The UCD_Start routine creates the UCD-Mesh-Output from it:
    type(t_triangulation), pointer :: UCD_MeshOnePointer => NULL()
    type(t_triangulation), pointer :: UCD_MeshTwoPointer => NULL()
    ! Write out the functions?
    logical :: ucd_OUT_orig_functions_one = .false.
    logical :: ucd_OUT_orig_functions_two = .false.
    ! If yes, where to save?
    character(LEN=ExtFE_STRLEN) :: UCD_FEfunctionOneOrigOutPath = ''
    character(LEN=ExtFE_STRLEN) :: UCD_FEfunctionTwoOrigOutPath = ''
    ! We also need some space to save the FE-Function
    type(t_vectorBlock), pointer :: UCD_feFunction_first_orig => NULL()
    type(t_vectorBlock), pointer :: UCD_feFunction_second_orig => NULL()
    ! For the UCD-Output, we need to create a Linear/Constant
    ! discretisation and project out FE-Function on that one as some
    ! UCD-writers only understand Q1/Q0 vectors:
    type(t_blockDiscretisation), pointer :: UCDBlockDiscrFirst => NULL()
    type(t_blockDiscretisation), pointer :: UCDBlockDiscrSecond => NULL()
    ! We can create this discretisation for each variable. These int
    ! stores the handle of an Array (size: nVar) which tells us if we use
    ! Q1 or Q0 discretisation for this variable. Entry i belongs to variable i
    ! One of the constants ExtFE_UCD_POLY_CONST, ExtFE_UCD_POLY_LINEAR
    integer :: h_UCD_AddElemProjectFirst = ST_NOHANDLE
    integer :: h_UCD_AddElemProjectSecond = ST_NOHANDLE
    ! We can add each variable vertex or element based.
    ! This int stores the handle of an array which tells us
    ! how to do. Entry i belongs to variable i
    ! One of the constants ExtFE_UCD_VERT_BASED, ExtFE_UCD_ELEM_BASED
    integer :: h_UCD_AddTypeOrigFirst = ST_NOHANDLE
    integer :: h_UCD_AddTypeOrigSecond = ST_NOHANDLE
    ! We can treat vectors different in the output:
    ! This int stores a handle to an array of Ints, size:
    ! (/nVecs,nDimension/), entry (i,1)
    ! tells us the x-comp of vector i, (i,2) the y and (i,3) the z component.
    ! To be exact, it tells us which variable of the input function
    ! is the x/y/z-component of that vector
    ! Note: For every vector, all components have to be set! That is, if
    ! you have a 3D-Problem you must specify a 3D-Vector.
    ! ASSUMTION: EVERY Variable in a vector is written out
    ! with the same Add-type!
    integer :: h_UCD_VecsFirstOrig = ST_NOHANDLE
    integer :: h_UCD_VecsSecondOrig = ST_NOHANDLE
    ! And at last, these int stores a handle to
    ! an array where all scalar variables (that are not
    ! part of a vector as specified above) are listed.
    ! Size: iScalarVars
    ! for CC2D, 1 and 2 would be in the vector list
    ! (as 1 is speed in x, 2 is speed in y-direction)
    ! and 3 would end here (since it is the pressure and scalar)
    integer :: h_UCD_ScalarFirstOrig = ST_NOHANDLE
    integer :: h_UCD_ScalarSecondOrig = ST_NOHANDLE
    ! Later, we might add the possibility to write out
    ! discontinious solutions. Atm this can only store
    ! the global UCD-Variable`UCD_FLAG_STANDARD
    integer(I32) :: UCD_Style = -1

  end type

end module
