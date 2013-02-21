!##############################################################################
!# ****************************************************************************
!# <name> boundary </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module realises the handling of parametrised boundary components.
!# It provides support for lines, circles and arcs in 2D. The parametrisation
!# can be read from a .prm file.
!#
!# The following routines are provided by this module:
!#
!# 1.) boundary_read_prm
!#     -> Read a .prm file, create a boundary structure
!#
!# 2.) boundary_release
!#     -> Releases a boundary structure from the heap
!#
!# 3.) boundary_igetNBoundComp
!#     -> Get the number of boundary components
!#
!# 4.) boundary_igetNsegments
!#     -> Get the number of boundary segments in a boundary component
!#
!# 5.) boundary_dgetMaxParVal
!#     -> Return the maximum parameter value of a boundary component
!#
!# 6.) boundary_getCoords = boundary_getCoords /
!#                          boundary_getCoords_mult /
!#                          boundary_getCoords_sim
!#     -> Calculate the coordinates of a point given by its parameter value
!#
!# 7.) boundary_createRegion
!#     -> Get the characteristics of a boundary segment and create
!#        a boundary region structure from it.
!#
!# 8.) boundary convertRegion
!#     -> Converts the parameter values of a boundary region from 0-1
!#        parametrisation to length parametrisation and back
!#
!# 9.) boundary_isInRegion
!#     -> Tests whether a node with a specific parameter value
!#        is inside of a given boundary region or not.
!#
!# 10.) boundary_convertParameter
!#      -> Allows to convert a parameter value from 0-1 parametrisation to
!#         length parametrisation and back
!#
!# 11.) boundary_convertParameterList = boundary_convertParameter_mult /
!#                                      boundary_convertParameter_sim
!#      -> Converts a list of parameter values from 0-1 parametrisation to
!#         length parametrisation and back
!#
!# 12.) boundary_getRegionLength
!#      -> Calculates the length of a boundary region.
!#
!# 13.) boundary_getNormalVec2D = boundary_getNormalVec2D /
!#                                boundary_getNormalVec2D_mult /
!#                                boundary_getNormalVec2D_sim
!#      -> Calculate the outward unit normal vector
!#         of a boundary component in 2D
!#
!# 14.) boundary_calcNormalVec2D = boundary_calcNormalVec2D /
!#                                 boundary_calcNormalVec2D_mult /
!#                                 boundary_calcNormalVec2D_sim
!#     -> Calculates the normal to a curve given by a set
!#        of sampling points at a point on the curve.
!#
!# It contains the following set of auxiliary routines:
!#
!# 1.) boundary_getNormal2D
!#     -> Calculates the normal vector at a point on a specific boundary segment.
!#
!# </purpose>
!##############################################################################

module boundary

!$use omp_lib
  use storage
  use fsystem
  use error
  use io
  use genoutput
  use fparser

  implicit none

  private

!<constants>
!<constantblock description="Global constants for boundary">
  ! maximal degree of NURBS
  integer, parameter, public :: PAR_MAXNURBSDEGREE = 35
!</constantblock>

!<constantblock description="types of boundary segments">
  ! boundary segment type line
  integer, parameter, public :: BOUNDARY_TYPE_LINE = 1

  ! boundary segment type circle
  integer, parameter, public :: BOUNDARY_TYPE_CIRCLE = 2

  ! boundary segment type open nurbs
  integer, parameter, public :: BOUNDARY_TYPE_OPENNURBS = 3

  ! boundary segment type closed nurbs
  integer, parameter, public :: BOUNDARY_TYPE_CLOSEDNURBS = 4

  ! boundary segment analytic
  integer, parameter, public :: BOUNDARY_TYPE_ANALYTIC = 5
  
  ! boundary segment by expression
  integer, parameter, public :: BOUNDARY_TYPE_EXPRESSION = 6
!</constantblock>

!<constantblock description="kinds of boundary segments">
  ! boundary kind fictitious
  integer, parameter, public :: BOUNDARY_KIND_FICTITIOUS = 0

  ! boundary kind geometric
  integer, parameter, public :: BOUNDARY_KIND_GEOMETRIC = 1
!</constantblock>

!<constantblock description="boundary segment header definition">

  ! boundary segment header offset for type
  integer, parameter, public :: BOUNDARY_SEGHEADER_TYPE = 1

  ! boundary segment header offset for offset in the data vector
  integer, parameter, public :: BOUNDARY_SEGHEADER_OFFSET = 2

  ! boundary segment header offset for nurbs degree
  integer, parameter, public :: BOUNDARY_SEGHEADER_NURBSDEGREE = 3

  ! boundary segment header offset for analytic expression
  integer, parameter, public :: BOUNDARY_SEGHEADER_EXPRESSION = 3

  ! boundary segment header offset for number of control points
  integer, parameter, public :: BOUNDARY_SEGHEADER_NCNTRLPNTS = 4

  ! boundary segment header length
  integer, parameter, public :: BOUNDARY_SEGHEADER_LENGTH = 4

!</constantblock>

!<constantblock description="Boundary region type qualifier">

  ! The boundary region is a usual curved region specified
  ! by min/max. parameter value
  integer, parameter, public :: BDR_TP_CURVE = 0

  ! The boundary region is a point region, i.e. consisting of
  ! only one point. Min/Max. parameter values are identical.
  integer, parameter, public :: BDR_TP_POINT = 1

!</constantblock>

!<constantblock description="Bitfield constants for t\_boundaryRegion%iproperties">

  ! The startpoint (point with min. parameter value) belongs
  ! to the boundary segment. If not set, the point with the
  ! min. parameter value does not belong to the boundary segment.
  integer, parameter, public :: BDR_PROP_WITHSTART = 2**0

  ! The endpoint (point with min. parameter value) belongs
  ! to the boundary segment. If not set, the point with the
  ! max. parameter value does not belong to the boundary segment.
  integer, parameter, public :: BDR_PROP_WITHEND = 2**1

!</constantblock>

!<constantblock description="Type constants for type of parametrisation to use.">

  ! Use 0-1 parametrisation
  integer, parameter, public :: BDR_PAR_01       = 0

  ! Use parametrisation for the arc length
  integer, parameter, public :: BDR_PAR_LENGTH   = 1

!</constantblock>

!<constantblock description="Type constants for cnormalMean in the calculation of normal vectors.">

  ! Calculate the mean of the left and right normal.
  integer, parameter, public :: BDR_NORMAL_MEAN         = 0

  ! Calculate the right normal (i.e. from the point to the interval with
  ! increasing parameter value).
  integer, parameter, public :: BDR_NORMAL_RIGHT        = 1

  ! Calculate the left normal (i.e. from the point to the interval with
  ! decreasing parameter value).
  integer, parameter, public :: BDR_NORMAL_LEFT         = 2

!</constantblock>

!</constants>

!<types>

!<typeblock>

  ! Structure defining a boundary region with minimum/maximum parameter
  ! value in 2D.
  type t_boundaryRegion

    ! Type of parametrisation, the parameters refer to.
    ! One of the BDR_PAR_xxxx constants.
    integer  :: cparType  = BDR_PAR_01

    ! Minimum parameter value of the region
    real(DP) :: dminParam = 0.0_DP

    ! Maximum parameter value of the region
    real(DP) :: dmaxParam = 0.0_DP

    ! Type of the region. One of the BDR_TP_xxxx constants
    integer :: ctype = BDR_TP_CURVE

    ! Type of the boundary segment corresponding to this boundary
    ! region. One of the BOUNDARY_TYPE_xxxx constants.
    integer :: isegmentType = BOUNDARY_TYPE_LINE

    ! Number of the boundary component that contains the segment
    integer :: iboundCompIdx = 0

    ! Maximum parameter value of the boundary component iboundCompIdx
    real(DP) :: dmaxParamBC = 0.0_DP

    ! Number of the boundary segment that corresponds to this region
    ! =0 if the boundary region does not correspond to a specific
    ! boundary segment.
    integer :: iboundSegIdx = 0

    ! Bitfield specifying properties of the region. A combination
    ! of BDR_PROP_xxxx constants.
    integer :: iproperties = BDR_PROP_WITHSTART

  end type t_boundaryRegion

  public :: t_boundaryRegion

!</typeblock>

!<typeblock>

  ! Boundary structure of the domain
  type t_boundary

    private

    ! number of geometric boundary components
    integer :: iboundarycount_g = -1

    ! total number of boundary components
    integer :: iboundarycount   = -1

    ! handle to double precision array: For every boundary component, maximum
    ! parameter value in length-parametrisation.
    integer :: h_DmaxPar = ST_NOHANDLE

    ! handle for a vector containing the number of segments per boundary component
    integer :: h_IsegCount = ST_NOHANDLE

    ! contains handles for data vectors of boundary components
    integer :: h_Idbldatavec_handles = ST_NOHANDLE

    ! contains handles for offset vectors of boundary components
    integer :: h_Iintdatavec_handles = ST_NOHANDLE

    ! contains offsets in function parser object
    integer :: h_Iintdatavec_fparser = ST_NOHANDLE

    ! Function parser storing analytic expressions
    type(t_fparser), pointer :: p_rfparser => null()

  end type t_boundary

  public :: t_boundary

!</typeblock>
! </types>

  interface boundary_getCoords
    module procedure boundary_getCoords
    module procedure boundary_getCoords_mult
  end interface

  interface boundary_convertParameterList
    module procedure boundary_convertParameter_mult
    module procedure boundary_convertParameter_sim
  end interface

  interface boundary_getNormalVec2D
    module procedure boundary_getNormalVec2D
    module procedure boundary_getNormalVec2D_mult
    module procedure boundary_getNormalVec2D_sim
  end interface

  interface boundary_calcNormalVec2D
    module procedure boundary_calcNormalVec2D
    module procedure boundary_calcNormalVec2D_mult
    module procedure boundary_calcNormalVec2D_sim
  end interface

  public :: boundary_read_prm
  public :: boundary_release
  public :: boundary_igetNBoundComp
  public :: boundary_igetNsegments
  public :: boundary_dgetMaxParVal
  public :: boundary_getCoords
  public :: boundary_createRegion
  public :: boundary_isInRegion
  public :: boundary_convertRegion
  public :: boundary_convertParameter
  public :: boundary_convertParameterList
  public :: boundary_convertParameter_mult
  public :: boundary_convertParameter_sim
  public :: boundary_getRegionLength
  public :: boundary_getNormalVec2D
  public :: boundary_getNormalVec2D_mult
  public :: boundary_getNormalVec2D_sim
  public :: boundary_calcNormalVec2D
  public :: boundary_calcNormalVec2D_mult
  public :: boundary_calcNormalVec2D_sim

contains

  !************************************************************************

!<function>

  integer function boundary_igetNBoundComp(rboundary) result (iresult)

!<description>
  ! This function returns the total number of boundary components.
!</description>

!<result>
  ! Total number of boundary components
!</result>

!<input>

  ! boundary structure
  type(t_boundary), intent(in) :: rboundary

!</input>

!</function>

    iresult = rboundary%iboundarycount

  end function boundary_igetNBoundComp

  !************************************************************************

!<function>

  recursive real(DP) function boundary_dgetLength (rboundary, &
                     iboundCompIdx, dt1, dt2) result (dresult)

  !<description>
  ! This function returns the euclidian length between the point
  ! dt1 on the boundary curve and the point dt2.
  !</description>

  !<result>
  ! Real euclidian length
  !</result>

  !<input>

  ! boundary structure
  type(t_boundary), intent(in) :: rboundary

  ! boundary index
  integer :: iboundCompIdx

  ! start parameter value
  real(DP), intent(in) :: dt1

  ! end parameter value
  real(DP), intent(in) :: dt2

  !</input>

!</function>

    real(DP) :: dx1,dx2,dy1,dy2,dbl,dtmax
    integer :: ip2,ip1,i

    if ((iboundCompIdx .gt. rboundary%iboundarycount) .or.&
        (iboundCompIdx.lt.0)) then
      call output_line ('iboundCompIdx out of bounds!', &
                        OU_CLASS_WARNING,OU_MODE_STD,'boundary_dgetLength')
      dresult = -1
      return
    endif

    ip1=ceiling(dt1)
    ip2=floor(dt2)

    if (ip2.lt.ip1) then

      call boundary_getCoords(rboundary,iboundCompIdx,dt1,dx1,dy1)
      call boundary_getCoords(rboundary,iboundCompIdx,dt2,dx2,dy2)

      dbl=sqrt((dx2-dx1)*(dx2-dx1)+(dy2-dy1)*(dy2-dy1))

    else

      if (dt2.lt.dt1) then

        ! Special case: Endpoint wrapped from the end of the boundary
        ! component back to the beginning.

        dtmax=boundary_dgetMaxParVal(rboundary,iboundCompIdx)

        dbl=boundary_dgetLength(rboundary,iboundCompIdx,dt1,dtmax)
        dbl=dbl+boundary_dgetLength(rboundary,iboundCompIdx,dtmax,dt2)

      else

        if (real(ip1,DP).eq.dt1) then
          call boundary_getCoords(rboundary,iboundCompIdx,dt1,dx2,dy2)
          dbl=0.0_DP
        else
          call boundary_getCoords(rboundary,iboundCompIdx,dt1,dx1,dy1)
          call boundary_getCoords(rboundary,iboundCompIdx,real(ip1,DP),dx2,dy2)
          dbl=sqrt((dx2-dx1)*(dx2-dx1)+(dy2-dy1)*(dy2-dy1))
        endif

        do i=ip1+1,ip2

          dx1=dx2
          dy1=dy2

          call boundary_getCoords(rboundary,iboundCompIdx,real(i,DP),dx2,dy2)
          dbl=dbl+sqrt((dx2-dx1)*(dx2-dx1)+(dy2-dy1)*(dy2-dy1))

        enddo

        if (real(ip2,DP).ne.dt2) then
          call boundary_getCoords(rboundary,iboundCompIdx,dt2,dx1,dy1)
          dbl=dbl+sqrt((dx1-dx2)*(dx1-dx2)+(dy1-dy2)*(dy1-dy2))
        endif

      endif

    endif

    dresult=dbl

  end function boundary_dgetLength


  !************************************************************************

!<function>

  real(DP) function boundary_dgetMaxParVal(rboundary, iboundCompIdx, cparType)

!<description>
  ! This function returns the parametric length of the boundary component iboundCompIdx
!</description>

!<result>
  ! Parametric length of boundary component iboundCompIdx.
!</result>

!<input>

  ! boundary structure
  type(t_boundary), intent(in) :: rboundary

  ! index of boundary component
  integer, intent(in) :: iboundCompIdx

  ! OPTIONAL: Type of parametrisation to use.
  ! One of the BDR_PAR_xxxx constants. If not given, BDR_PAR_01 is assumed.
  integer, intent(in), optional :: cparType

!</input>

!</function>

    real(DP),dimension(:),pointer :: p_DmaxPar
    integer,dimension(:),pointer :: p_IsegCount

    ! if iboundCompIdx exceeds the total number of boundary components or is negative, abort
    if ((iboundCompIdx .gt. rboundary%iboundarycount) .or. (iboundCompIdx.lt.0)) then
      call output_line ('iboundCompIdx out of bounds!', &
          OU_CLASS_WARNING,OU_MODE_STD,'boundary_dgetMaxParVal')
      boundary_dgetMaxParVal = -1
      return
    endif

    if (present(cparType)) then
      select case (cparType)
      case (BDR_PAR_LENGTH)
        ! Length-parametrisation
        ! Get vector with component length - depending on the parametrisation
        ! to use.
        call storage_getbase_double(rboundary%h_DmaxPar,p_DmaxPar)

        ! Get the maximum parameter value
        boundary_dgetMaxParVal = p_DmaxPar(iboundCompIdx)

        return
      end select
    end if

    ! 0-1 parametrisation. Maximum parameter value is just the number of
    ! segments.
    call storage_getbase_int(rboundary%h_IsegCount,p_IsegCount)
    boundary_dgetMaxParVal = real(p_IsegCount(iboundCompIdx),DP)

  end function boundary_dgetMaxParVal

  !************************************************************************

!<function>

  integer function boundary_igetNsegments(rboundary, iboundCompIdx)

!<description>
  ! This function returns the number of boundary segments in
  ! the boundary component iboundCompIdx.
!</description>

!<result>
  ! Number of boundary segments on boundary component iboundCompIdx.
!</result>

!<input>

  ! boundary structure
  type(t_boundary), intent(in) :: rboundary

  ! index of boundary component
  integer, intent(in) :: iboundCompIdx

!</input>

!</function>

    integer,dimension(:),pointer :: p_IsegCount

    ! if iboundCompIdx exceeds the total number of boundary components or is negative, abort
    if ((iboundCompIdx .gt. rboundary%iboundarycount) .or. (iboundCompIdx.lt.0)) then
      call output_line ('iboundCompIdx out of bounds!', &
          OU_CLASS_WARNING,OU_MODE_STD,'boundary_dgetNsegments')
      boundary_igetNsegments = -1
      return
    endif

    ! get vector with component length
    call storage_getbase_int(rboundary%h_IsegCount,p_IsegCount)

    ! Get the maximum parameter value
    boundary_igetNsegments = p_IsegCount(iboundCompIdx)

  end function boundary_igetNsegments

  !************************************************************************

!<subroutine>

  subroutine boundary_read_prm(rboundary, sfilename)

!<description>
  ! This routine reads a .PRM file into memory. The boundary structure
  ! rboundary is initialised with the data from the file.
  ! The parameter sfilename gives the name of the .prm file to read.
  ! If p_rboundary is NULL(), a new structure will be created.
  ! Otherwise, the existing structure is recreated/updated.
!</description>

!<input>
  ! The name of the .prm file to read.
  character(LEN=*), intent(in) :: sfilename
! </input>

!<output>
  ! Boundary structure, to be filled with data
  type(t_boundary), intent(out) :: rboundary
!</output>

!</subroutine>

    ! local variables

    ! Input channel for reading
    integer :: iunit
    integer :: ibcomponent, isegment, ibct,ihandle
    integer, dimension(:), pointer :: p_IsegInfo, p_IsegCount
    integer, dimension(:), pointer :: p_IdbleSegInfo_handles,p_IintSegInfo_handles
    integer, dimension(:), pointer :: p_IintSegInfo_fparser
    real(DP), dimension(:), pointer :: p_DsegInfo, p_DmaxPar
    integer :: ityp,nspline,npar
    integer :: i,isegrel,icomp,nurbsdegree,ncntrlpnts,nknots
    integer :: idblemem  ! Counts the memory we need
    integer :: iexpression  ! Counts the number of expressions
    real(DP) :: dl,dtemp
    character(len=1024) :: cbuffer
    character(len=1), dimension(1) :: Cpar
    real(DP), dimension(1) :: D1,D2,D3

    ! Current parameter value in length-parametrisation
    real(DP) :: dmaxpar

    ! Open the file
    call io_openFileForReading(sfilename, iunit)

    ! Read "NBCT"
    read (iunit,*)

    ! Read NBCT - Number of boundary components
    read (iunit,*) rboundary%iboundarycount_g

    rboundary%iboundarycount = rboundary%iboundarycount_g

    ! Allocate an array containing handles. Each handle refers
    ! to integer data for a boundary component.
    call storage_new("boundary_read_prm", "h_Idbldatavec_handles", &
                       rboundary%iboundarycount, ST_INT, &
                       rboundary%h_Idbldatavec_handles, ST_NEWBLOCK_ZERO)
    call storage_getbase_int(rboundary%h_Idbldatavec_handles, p_IdbleSegInfo_handles)

    ! Allocate an array containing of handles. Each handle refers
    ! to integer data for a boundary component.
    call storage_new("boundary_read_prm", "h_Iintdatavec_handles", &
                     rboundary%iboundarycount, ST_INT, &
                     rboundary%h_Iintdatavec_handles, ST_NEWBLOCK_ZERO)
    call storage_getbase_int(rboundary%h_Iintdatavec_handles, p_IintSegInfo_handles)

    ! Allocate an array containing the maximum parameter values for each
    ! boundary component in length-parametrisation
    call storage_new("boundary_read_prm", "h_DmaxPar", &
                     rboundary%iboundarycount, ST_DOUBLE, &
                     rboundary%h_DmaxPar, ST_NEWBLOCK_ZERO)
    call storage_getbase_double(rboundary%h_DmaxPar, p_DmaxPar)

    ! Allocate an array containing the number of boundary segments in each
    ! boundary component
    call storage_new("boundary_read_prm", "h_IsegCount", &
                     rboundary%iboundarycount, ST_INT, &
                     rboundary%h_IsegCount, ST_NEWBLOCK_ZERO)
    call storage_getbase_int(rboundary%h_IsegCount, p_IsegCount)

    ! Allocate an array containing the offsets of boundary segments in
    ! the global function parser object
    call storage_new("boundary_read_prm", "h_Iintdatavec_fparser", &
                     rboundary%iboundarycount, ST_INT, &
                     rboundary%h_Iintdatavec_fparser, ST_NEWBLOCK_ZERO)
    call storage_getbase_int(rboundary%h_Iintdatavec_fparser, p_IintSegInfo_fparser)

    ! No we have to gather information about boundary segments.
    ! This makes it necessary to read the boundary definition file
    ! multiple times, but who cares :)

    ! iexpression counts how many analytic expressions we need in total
    iexpression = 0
    
    ! Initialise the boundary components, loop through them

    do ibcomponent = 1,rboundary%iboundarycount_g
      ! Read "IBCT"
      read (iunit,*)

      ! Read IBCT
      read (iunit,*) ibct

      ! Read "NCOMP"
      read (iunit,*)

      ! Read NCOMP = number of boundary segments in this boundary component
      read (iunit,*) p_IsegCount(ibct)

      ! Allocate an integer-array which saves the segment type and
      ! the start positions (0-based!) of each segment information
      ! block in the segment-array. It is as long as the number
      ! of segments indicates * 2.

      call storage_new("boundary_read_prm", "h_Isegcount", &
                       BOUNDARY_SEGHEADER_LENGTH*p_IsegCount(ibct), ST_INT, &
                       ihandle, ST_NEWBLOCK_ZERO)
      p_IintSegInfo_handles(ibct) = ihandle
      call storage_getbase_int(ihandle, p_IsegInfo)

      ! Evaluate the boundary segment definition to get the memory
      ! needed for it:
      read (iunit,*)         ! "ITYP NSPLINE NPAR"

      ! idblemem counts how much memory we need.
      idblemem = 0

      ! evaluate the "ITYP NSPLINE NPAR" block
      do isegment = 0,p_IsegCount(ibct)-1

        ! read ITYP NSPLINE NPAR of that segment
        read (iunit,*) ityp, nspline, npar

        ! What do we have here?
        select case (ityp)
        case (BOUNDARY_TYPE_LINE)
          ! Type 1: Line.
          ! Save the segment type into the first element of each
          ! 2-tuple in the integer array:

          p_IsegInfo (BOUNDARY_SEGHEADER_LENGTH*isegment+&
                      BOUNDARY_SEGHEADER_TYPE) = BOUNDARY_TYPE_LINE

          ! Save the start position of this segment to the segment-start array.

          p_IsegInfo (BOUNDARY_SEGHEADER_LENGTH*isegment+&
                      BOUNDARY_SEGHEADER_OFFSET) = idblemem

          ! A line consists of
          ! - Startpoint
          ! - Endpoint, relative to the startpoint
          ! Furthermore, in the first element of the segment we save the
          ! minimum parmeter value and in the second one the length.
          ! So we need 2*2+2 doubles.
          idblemem = idblemem + 6

        case (BOUNDARY_TYPE_CIRCLE)
          ! Type 2: Circle / arc.
          ! Save the segment type into the first element of each
          ! 2-tuple in the integer array:

          p_IsegInfo (BOUNDARY_SEGHEADER_LENGTH*isegment+&
                      BOUNDARY_SEGHEADER_TYPE) = BOUNDARY_TYPE_CIRCLE

          ! Save the start position of this segment to the segment-start array.

          p_IsegInfo (BOUNDARY_SEGHEADER_LENGTH*isegment+&
                      BOUNDARY_SEGHEADER_OFFSET) = idblemem

          ! A circle/arc consists of
          ! - Center
          ! - Radius
          ! - Parameter value of the start point of the arc.
          !   Refers to the unit circle.
          ! - Parameter value of the end point of the arc.
          !   Refers to the unit circle.
          ! Furthermore, in the first element of the segment we save the
          ! minimum parameter value and in the second one the length.
          ! So we need 3*2+2 doubles.
          idblemem = idblemem + 8

        case (BOUNDARY_TYPE_OPENNURBS)
          ! Type 3: Open NURBS
          ! Save the segment type into the first element of each
          ! 2-tuple in the integer array:

          p_IsegInfo (BOUNDARY_SEGHEADER_LENGTH*isegment+&
                      BOUNDARY_SEGHEADER_TYPE) = BOUNDARY_TYPE_OPENNURBS

          ! Save the start position of this segment to the segment-start array.
          
          p_IsegInfo (BOUNDARY_SEGHEADER_LENGTH*isegment+&
                      BOUNDARY_SEGHEADER_OFFSET) = idblemem

          ! An open NURBS consists of
          ! - coordinates (x,y) of control points and weights w
          ! - knot parameter values
          ! Furthermore, in the first element of the segment we save the
          ! minimum parameter value and in the second one the length.
          !
          ! nspline = degree of NURBS
          ! npar    = number of control points
          !
          ! length of open NURBS segment:
          ! three values per control point (x,y,weight) 
          ! + 2 values for parameter value and length 
          ! + length of (uncompressed) node vector = nspline + npar + 1
          !
          ! So we need 3*npar  + (npar+1 + nspline) + 2 doubles
          idblemem = idblemem + 3*npar + (npar+1 + nspline) + 2

        case (BOUNDARY_TYPE_CLOSEDNURBS)
          ! Type 4: Closed NURBS
          ! Save the segment type into the first element of each
          ! 2-tuple in the integer array:

          p_IsegInfo (BOUNDARY_SEGHEADER_LENGTH*isegment+&
                      BOUNDARY_SEGHEADER_TYPE) = BOUNDARY_TYPE_CLOSEDNURBS

          ! Save the start position of this segment to the segment-start array.
          
          p_IsegInfo (BOUNDARY_SEGHEADER_LENGTH*isegment+&
                      BOUNDARY_SEGHEADER_OFFSET) = idblemem

          ! A closed NURBS consists of
          ! - coordinates (x,y) of control points and weights w
          ! - knot parameter values
          ! Furthermore, in the first element of the segment we save the
          ! minimum parameter value and in the second one the length.
          !
          ! nspline = degree of NURBS
          ! npar    = number of control points
          !
          ! length of closed NURBS segment:
          ! three values per control point (x,y,weight) 
          ! + 2 values for parameter value and length 
          ! + length of (uncompressed) node vector = nspline + npar + 1
          !
          ! So we need  3*npar + (npar+1) + 2 doubles
          idblemem = idblemem + 3*npar + (npar+1) + 2

        case (BOUNDARY_TYPE_EXPRESSION)
          ! Type 6: Analytic expression 
          ! Save the segment type into the first element of each
          ! 2-tuple in the integer array:

          p_IsegInfo (BOUNDARY_SEGHEADER_LENGTH*isegment+&
                      BOUNDARY_SEGHEADER_TYPE) = BOUNDARY_TYPE_EXPRESSION

          ! Save the start position of this segment to the segment-start array.

          p_IsegInfo (BOUNDARY_SEGHEADER_LENGTH*isegment+&
                      BOUNDARY_SEGHEADER_OFFSET) = idblemem

          ! Save the start position of this segment in the function parser.
          p_IsegInfo (BOUNDARY_SEGHEADER_LENGTH*isegment+&
                      BOUNDARY_SEGHEADER_EXPRESSION) = iexpression

          ! An analytic expression consists of
          ! - Expression for x(t)
          ! - Expression for y(t)
          ! - Expression for x'(t)
          ! - Expression for y'(t)
          ! So we need 4 expressions
          iexpression = iexpression + 4

          ! Furthermore, in the first element of the segment we save the
          ! minimum parmeter value and in the second one the length.
          ! So we need 2 doubles.
          idblemem = idblemem + 2
        end select
      end do

      ! Allocate memory for all the information that defines our
      ! boundary segments:
      call storage_new("boundary_read_prm", "h_IdbleSegInfo_handles", &
                       idblemem, ST_DOUBLE, &
                       ihandle, ST_NEWBLOCK_ZERO)
      p_IdbleSegInfo_handles(ibct) = ihandle

    end do

    ! Do we need to allocate a funciton parser?
    if (iexpression .gt. 0) then
      allocate(rboundary%p_rfparser)
      call fparser_create(rboundary%p_rfparser, iexpression)
    end if

    ! Now we build the segment information.

    ! Ignore "PARAMETER"
    read (iunit,*)  ! "PARAMETER"

    ! Read the boundary information
    do ibcomponent = 1,rboundary%iboundarycount_g

      ! dmaxPar counts for every boundary component the length - and
      ! thus the maximum parameter value.
      dmaxPar = 0.0_DP

      ! Get a pointer to the boundary component info.
      ! Here, all double-precision information of the current boundary
      ! component is saved.
      call storage_getbase_double(&
          int(p_IdbleSegInfo_handles(ibcomponent)), p_DsegInfo)

      ! Get a pointer to the start position of the segments in the
      ! current boundary component.
      call storage_getbase_int(&
          int(p_IintSegInfo_handles(ibcomponent)), p_ISegInfo)

      ! Build up the segments in the block:
      do isegment = 0,p_IsegCount(ibcomponent)-1

        ! Get the relative position of the segment information in the
        ! segment information array.
        isegrel = p_IsegInfo(BOUNDARY_SEGHEADER_LENGTH*isegment+BOUNDARY_SEGHEADER_OFFSET)

        ! What do we have here? The type identifier is in the
        ! first entry of the 2-tupel in the integer-array:
        select case (p_IsegInfo(BOUNDARY_SEGHEADER_LENGTH*isegment+BOUNDARY_SEGHEADER_TYPE))
        case (BOUNDARY_TYPE_LINE)
          ! Type 1: Line. Consists of
          ! - Startpoint
          ! - Endpoint, relative to the startpoint
          ! We read the startpoint and the relative endpoint.
          ! The line is then saved as
          ! - startpoint
          ! - direction vector
          ! - length of the line.
          !
          ! Startpoint
          read(iunit,*) p_DsegInfo(isegrel+3),p_DsegInfo(isegrel+4)

          ! Relative endpoint
          read(iunit,*) p_DsegInfo(isegrel+5),p_DsegInfo(isegrel+6)

          ! Save the initial parameter value (in of the arc length-parametrisation)
          ! to the first entry
          p_DsegInfo(isegrel+1) = dmaxpar

          ! This is just necessary for the length-parametrisation. The
          ! 0-1 parametrisation defines the starting point just as the
          ! integer number isegment.

          ! Calculate the length and save it in the 2nd position
          dl = sqrt(p_DsegInfo(isegrel+5)**2 + p_DsegInfo(isegrel+6)**2)
          p_DsegInfo(isegrel+2) = dl

          ! Normalise the direction vector
          !p_DsegInfo(isegrel+4:isegrel+6) = p_DsegInfo(isegrel+4:isegrel+7)/dl

          ! Increase the maximum parameter value
          dmaxPar = dmaxPar + dl

        case (BOUNDARY_TYPE_CIRCLE)
          ! Type 2: Circle / arc. Consists of
          ! - Center
          ! - Radius
          ! - Parameter value of the start point of the arc.
          !   Refers to the unit circle.
          ! - Parameter value of the end point of the arc.
          !   Refers to the unit circle.
          !
          ! The first entry in the segment information array is again
          ! the length of the circle - to be computed later.
          ! First save the midpoint of the circle
          read(iunit,*) p_DsegInfo(isegrel+3),p_DsegInfo(isegrel+4)

          ! Then the radius - the second entry is a dummy;
          read(iunit,*) p_DsegInfo(isegrel+5),p_DsegInfo(isegrel+6)

          ! Finally get the "arc positions" of the startpoint and the endpoint
          ! of the arc:
          read(iunit,*) p_DsegInfo(isegrel+7),p_DsegInfo(isegrel+8)

          ! Save the initial parameter value (in length-parametrisation)
          ! to the first entry
          p_DsegInfo(isegrel+1) = dmaxpar

          ! Now compute the real length of the arc.
          dl = p_DsegInfo(isegrel+5) * &
              abs(p_DsegInfo(isegrel+8)-p_DsegInfo(isegrel+7))
          p_DsegInfo(isegrel+2) = dl

          ! Increase the maximum parameter value
          dmaxPar = dmaxPar + dl

        case (BOUNDARY_TYPE_OPENNURBS)
          ! Type 3: Open NURBS Consists of
          ! - coordinates (x,y) of control points and weights w
          ! - knot parameter values
          !
          ! nspline = number of control points
          ! npar    = number of lines in parameter section
          !           structure of PARAMETERS
          !           degree   #control pts   #knots
          !           x(i)     y(i)           w(i)      i=1..#control pts
          !           p(j)                              j=1..#knots
          !
          ! The first entry in the segment information array is again
          ! the length of the NURBS - to be computed later.
          read (iunit,*) nurbsdegree,ncntrlpnts,nknots

          ! First save the degree of the open NURBS
          p_IsegInfo (BOUNDARY_SEGHEADER_LENGTH*isegment+&
                      BOUNDARY_SEGHEADER_NURBSDEGREE) = nurbsdegree

          ! Then save the number of control points
          p_IsegInfo (BOUNDARY_SEGHEADER_LENGTH*isegment+&
                      BOUNDARY_SEGHEADER_NCNTRLPNTS) = ncntrlpnts

          ! Then save ncntrlpnts triples (x,y,w)
          do i = 0, ncntrlpnts-1
            read (iunit,*) p_DsegInfo(isegrel+3*i+1),&
                           p_DsegInfo(isegrel+3*i+2),&
                           p_DsegInfo(isegrel+3*i+3)
          end do

          ! The first entry of the node vector comes up ndegree+1
          ! times, therefore we have to repeat the writing ndegree times
          read(iunit,*) dtemp
          do i = 1, nurbsdegree + 1
            p_DsegInfo(isegrel+3*ncntrlpnts+i) = dtemp
          end do

          ! Read the interior entries of the node vector. The node
          ! vector has the (uncompressed) length ndegree + itemp
          ! +1. We have to subtract the start and end point, which
          ! both occur (ndegree+1) times.
          do i = 1, ncntrlpnts + nurbsdegree + 1 - 2*(nurbsdegree+1)
            read(iunit,*) p_DsegInfo(isegrel + 3*ncntrlpnts + i + nurbsdegree +1)
          end do

          ! Also the last control point parameter comes up (ndegree+1) times
          read(iunit,*) dtemp
          do i = 1, nurbsdegree + 1
            p_DsegInfo(isegrel+4*ncntrlpnts+i) = dtemp
          end do

          ! Save the initial parameter value (in length-parametrisation)
          ! to the first entry
          p_DsegInfo(isegrel+1) = dmaxpar

        case (BOUNDARY_TYPE_EXPRESSION)
          ! Type 6: Analytic expression. Consists of
          ! - Expression for x(t)
          ! - Expression for y(t)
          ! - Expression for x'(t)
          ! - Expression for y'(t)
          ! The four expressions must be 0-1 parametrised, with a parameter
          ! value t in [0,1].
          ! We read these four expressions into the function parser object.

          ! Get the absolute position of the analytic expression in the
          ! function parser object
          icomp = p_IsegInfo(BOUNDARY_SEGHEADER_LENGTH*isegment+BOUNDARY_SEGHEADER_EXPRESSION)
          
          Cpar = (/"t"/)
          
          read(iunit,fmt='(A)') cbuffer ! x(t)
          call sys_dequote(cbuffer)
          call fparser_parseFunction(rboundary%p_rfparser, icomp+1, cbuffer, cpar)
          read(iunit,fmt='(A)') cbuffer ! y(t)
          call sys_dequote(cbuffer)
          call fparser_parseFunction(rboundary%p_rfparser, icomp+2, cbuffer, cpar)
          read(iunit,fmt='(A)') cbuffer ! x'(t)
          call sys_dequote(cbuffer)
          call fparser_parseFunction(rboundary%p_rfparser, icomp+3, cbuffer, cpar)
          read(iunit,fmt='(A)') cbuffer ! y'(t)
          call sys_dequote(cbuffer)
          call fparser_parseFunction(rboundary%p_rfparser, icomp+4, cbuffer, cpar)

          ! Save the initial parameter value (in length-parametrisation)
          ! to the first entry
          p_DsegInfo(isegrel+1) = dmaxpar

          ! Calculate the real length of the analytically given path
          ! by adaptive simpson, applied to the derivative of the path.
          dl = adaptiveSimpson(rboundary%p_rfparser, icomp+3, icomp+4, 0.0_DP, 1.0_DP, 1e-8_DP, 10)
          p_DsegInfo(isegrel+2) = dl
          
          ! Increase the maximum parameter value
          dmaxPar = dmaxPar + dl

        end select
      end do

      ! Save the maximum parameter value for that component
      p_DmaxPar(ibcomponent) = dmaxPar

    end do

    ! Close the file, finish
    close(iunit)

  contains

    ! Here, some auxiliary functions follow

    ! Compute the length of the analytically given parameter curve
    ! by the adaptive Simpson quadrature rule
    function adaptiveSimpson(rfparser, icompx, icompy,&
        dparmin, dparmax, dtolerance, nrecursion) result(dlength)

      ! input parameters
      type(t_fparser), intent(in) :: rfparser
      integer, intent(in) :: icompx,icompy,nrecursion
      real(DP), intent(in) :: dparmin,dparmax,dtolerance
      
      ! result
      real(DP) :: dlength
      
      ! local variables
      real(DP) :: dparmid,dfmin,dfmax,dfmid,dl
      real(DP) :: daux1,daux2,daux3,daux4,daux5,daux6

      dparmid = (dparmin+dparmax)*0.5_DP
      
      D1 = (/dparmin/)
      D2 = (/dparmid/)
      D3 = (/dparmax/)
      
      call fparser_evalFunction(rfparser, icompx, D1, daux1)
      call fparser_evalFunction(rfparser, icompx, D2, daux2)
      call fparser_evalFunction(rfparser, icompx, D3, daux3)

      call fparser_evalFunction(rfparser, icompy, D1, daux4)
      call fparser_evalFunction(rfparser, icompy, D2, daux5)
      call fparser_evalFunction(rfparser, icompy, D3, daux6)

      dfmin = sqrt(daux1*daux1 + daux4*daux4)
      dfmid = sqrt(daux2*daux2 + daux5*daux5)
      dfmax = sqrt(daux3*daux3 + daux6*daux6)

      dl = (dparmax-dparmin)*(dfmin+4.0_DP*dfmid+dfmax)/6.0_DP

      dlength = adaptiveSimpsonAux(rfparser, icompx, icompy,&
          dparmin, dparmax, dtolerance, dl, dfmin, dfmid, dfmax, nrecursion)
    end function adaptiveSimpson

    ! Auxiliary routine for adaptive Simpson rule
    recursive function adaptiveSimpsonAux(rfparser, icompx, icompy, dparmin, dparmax,&
        dtolerance, dl, dfmin, dfmid, dfmax, nrecursion) result(dlength)

      ! input parameters
      type(t_fparser), intent(in) :: rfparser
      integer, intent(in) :: icompx,icompy,nrecursion
      real(DP), intent(in) :: dparmin,dparmax,dtolerance
      real(DP), intent(in) :: dl,dfmin,dfmid,dfmax
      
      ! result
      real(DP) :: dlength

      ! local variables
      real(DP) :: dparmid,dparleft,dparright,dfleft,dfright,dlleft,dlright
      real(DP) :: daux1,daux2,daux3,daux4

      dparmid   = (dparmin+dparmax)*0.5_DP
      dparleft  = (dparmin+dparmid)*0.5_DP
      dparright = (dparmid+dparmax)*0.5_DP
      
      D1 = (/dparleft/) 
      D2 = (/dparright/)

      call fparser_evalFunction(rfparser, icompx, D1, daux1)
      call fparser_evalFunction(rfparser, icompx, D2, daux2)

      call fparser_evalFunction(rfparser, icompy, D1, daux3)
      call fparser_evalFunction(rfparser, icompy, D2, daux4)

      dfleft  = sqrt(daux1*daux1 + daux3*daux3)
      dfright = sqrt(daux2*daux2 + daux4*daux4)

      dlleft  = (dparmax-dparmin)*(dfmin+4.0_DP*dfleft+dfmid)/12.0_DP
      dlright = (dparmax-dparmin)*(dfmid+4.0_DP*dfright+dfmax)/12.0_DP

      dlength = dlleft + dlright

      if ((nrecursion .le. 0) .or. (abs(dlength -dl) .le. 15*dtolerance)) then
        dlength = dlength + (dlength-dl)/15.0_DP
      else
        dlength = adaptiveSimpsonAux(rfparser, icompx, icompy, dparmin, dparmid,&
            0.5_DP*dtolerance, dlleft, dfmin, dfleft, dfmid, nrecursion-1) +&
                  adaptiveSimpsonAux(rfparser, icompx, icompy, dparmid, dparmax,&
            0.5_DP*dtolerance, dlright, dfmid, dfright, dfmax, nrecursion-1)
      end if

    end function adaptiveSimpsonAux

  end subroutine boundary_read_prm

  !************************************************************************

!<subroutine>

  subroutine boundary_release(rboundary)

!<description>
  ! This routine releases a boundary object from memory.
!</description>

!<inputoutput>
  ! Boundary structure, to be released.
  type(t_boundary), intent(inout) :: rboundary
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,ihandle
    integer, dimension(:), pointer :: p_IdbleSegInfo_handles,p_IintSegInfo_handles

    ! Check if data arrays are allocated
    if ((rboundary%h_Iintdatavec_handles .ne. ST_NOHANDLE) .and.&
        (rboundary%h_Idbldatavec_handles .ne. ST_NOHANDLE)) then

      ! Get the pointers to the segment information arrays for the current
      ! boundary component:
      call storage_getbase_int(rboundary%h_Iintdatavec_handles,p_IintSegInfo_handles)
      call storage_getbase_int(rboundary%h_Idbldatavec_handles,p_IdbleSegInfo_handles)

      ! Release the handles of the integer- and double-precision
      ! data blocks:
      do i=1,rboundary%iboundarycount
        ihandle = p_IintSegInfo_handles(i)
        call storage_free (ihandle)
        p_IintSegInfo_handles(i) = ihandle

        ihandle = p_IdbleSegInfo_handles(i)
        call storage_free (ihandle)
        p_IdbleSegInfo_handles(i) = ihandle
      end do

      ! Release data arrays in the structure
      call storage_free (rboundary%h_Iintdatavec_handles)
      call storage_free (rboundary%h_Idbldatavec_handles)
    end if

    ! Release all arrays in the structure
    if (rboundary%h_IsegCount .ne. ST_NOHANDLE)&
        call storage_free (rboundary%h_IsegCount)
    if (rboundary%h_DmaxPar .ne. ST_NOHANDLE)&
        call storage_free (rboundary%h_DmaxPar)
    if (rboundary%h_Iintdatavec_fparser .ne. ST_NOHANDLE)&
        call storage_free (rboundary%h_Iintdatavec_fparser)

    ! Release function parser if required
    if (associated(rboundary%p_rfparser)) then
      call fparser_release(rboundary%p_rfparser)
      deallocate(rboundary%p_rfparser)
      rboundary%p_rfparser => null()
    end if

    rboundary%iboundarycount_g = 0
    rboundary%iboundarycount = 0

    ! That is it...

  end subroutine boundary_release

  !************************************************************************

!<subroutine>

  pure subroutine boundary_normaliseParValue2D(IsegCount,DmaxPar,iboundCompIdx,cpar,dpar)

!<description>
  ! INTERNAL SUBROUTINE.
  ! Normalises the parameter value dpar into the range [0,max.par. value).
!</description>

!<input>
  ! Segment-count array
  integer, dimension(:), intent(in) :: IsegCount

  ! Array wirth maximum parameter values for all BC`s
  real(DP), dimension(:), intent(in) :: DmaxPar

  ! Number of the boundary component that is under consideration
  integer, intent(in) :: iboundCompIdx

  ! Type of parametrisation, format of dpar (0-1, length par.,...)
  integer, intent(in) :: cpar
!</input>

!<inputoutput>
  ! Parameter value to be normalised. Is replaced by the normalised
  ! parameter value.
  real(DP), intent(inout) :: dpar
!</inputoutput>

!</subroutine>

    select case (cpar)
    case (BDR_PAR_01)
      dpar = mod(dpar,real(IsegCount(iboundCompIdx),DP))
      !IF (dt .GE. REAL(p_IsegCount(iboundCompIdx),DP) ) THEN
      !  dpar = 0.0_DP
      !ELSE
      !  dpar = dt
      !ENDIF
    case (BDR_PAR_LENGTH)
      dpar = mod(dpar,DmaxPar(iboundCompIdx))
      !IF (dt .GE. p_DmaxPar(iboundCompIdx) ) THEN
      !  dpar = 0.0_DP
      !ELSE
      !  dpar = dt
      !ENDIF
    end select

  end subroutine boundary_normaliseParValue2D

  !************************************************************************

!<subroutine>

  pure subroutine boundary_getSegmentInfo2D(&
      IsegCount,DmaxPar,IsegInfo,DsegInfo,iboundCompIdx,dpar,cparType,iorientation,&
      iseg,istartidx,dcurrentpar,dendpar,dseglength,dparloc,isegtype)

!<description>
  ! INTERNAL SUBROUTINE.
  ! From a given parameter value dpar in a given parametrisation cpar,
  ! this routine calculates information about the corresponding boundary
  ! segment that contains the point.
!</description>

!<input>
    ! Integer segment array for segment counter
    integer, dimension(:), intent(in) :: IsegCount

    ! Double precision array defining the maximum parameter values
    ! of each boundary component
    real(DP), dimension(:), intent(in) :: DmaxPar

    ! Integer segment info array
    integer, dimension(:), intent(in) :: IsegInfo

    ! Double precision segment info array
    real(DP), dimension(:), intent(in) :: DsegInfo

    ! Number of the boundary component
    integer, intent(in) :: iboundCompIdx

    ! Parameter value of the point of interest
    real(DP), intent(in) :: dpar

    ! Type of parametrisation of dpar (BDR_PAR_01, BDR_PAR_LENGTH,...)
    integer, intent(in) :: cparType

    ! How to orient if the boundary segment is not unique (e.g. on
    ! corner points of the discretisation).
    ! =0: return the segment following the point
    ! =1: return the segment previous to the point
    integer, intent(in) :: iorientation
!</input>

!<output>
    ! Segment number (0,1,2,...) of the segment that contains dpar
    integer, intent(out) :: iseg

    ! Start index of the segment in IsegInfo that contains dpar
    integer, intent(out) :: istartidx

    ! Start parameter value of the segment that contains dpar
    real(DP), intent(out) :: dcurrentpar

    ! End parameter value of the segment that contains dpar
    real(DP), intent(out) :: dendpar

    ! Segment length of the segment that contains dpar
    real(DP), intent(out) :: dseglength

    ! Local parameter value of dpar inside of the segment
    real(DP), intent(out) :: dparloc

    ! Type of the segment
    integer, intent(out) :: isegtype
!</output>

!</subroutine>

    ! Determine the segment
    select case (cparType)
    case (BDR_PAR_01)

      ! Easy case: 0-1 parametrisation
      !
      ! Standard case: point somewhere in the inner of the segment
      iseg = aint(dpar)
      dcurrentpar = real(iseg,DP)

      ! Subtract the start position of the current boundary component
      ! from the parameter value to get the 'local' parameter value
      ! (0 .le. dparloc .le. length(segment))
      dparloc = dpar - dcurrentpar

      if ((iorientation .ne. 0) .and. (dparloc .eq. 0.0_DP)) then
        ! We are in a point where the segment is not unique and
        ! the caller wants to have the 'previous' segment. Find it!

        iseg = aint(dpar)-1
        if (iseg .eq. -1) then
          ! Get the last segment!
          iseg = IsegCount(iboundCompIdx)-1
          dcurrentpar = real(iseg,DP)
          dparloc = dpar + real(IsegCount(iboundCompIdx),DP) - dcurrentpar
        else
          dcurrentpar = real(iseg,DP)
          dparloc = dpar - dcurrentpar
        end if

      end if

      ! Determine Start index of the segment in the double-prec. block
      istartidx = IsegInfo(BOUNDARY_SEGHEADER_LENGTH*iseg+BOUNDARY_SEGHEADER_OFFSET)

      dendpar     = iseg + 1.0_DP
      dseglength  = DsegInfo(2+istartidx)

    case (BDR_PAR_LENGTH)

      ! In the length-parametrisation, we have to search.
      ! The orientation flag tells us whether to search for "<" or
      ! "<="!
      if (iorientation .ne. 0) then
        do iseg = 0,IsegCount(iboundCompIdx)-1

          ! Determine Start index of the segment in the double-prec. block
          istartidx = IsegInfo(BOUNDARY_SEGHEADER_LENGTH*iseg+BOUNDARY_SEGHEADER_OFFSET)

          ! Get the start and end parameter value
          dcurrentpar = DsegInfo(1+istartidx)
          dendpar = dcurrentpar + DsegInfo(2+istartidx)
          dseglength = DsegInfo(2+istartidx)

          ! At least one of the IF-commands in the loop will activate
          ! the exit - because of the 'dt' check above!
          if (dpar .le. dendpar) exit

        end do
      else
        do iseg = 0,IsegCount(iboundCompIdx)-1

          ! Determine Start index of the segment in the double-prec. block
          istartidx = IsegInfo(BOUNDARY_SEGHEADER_LENGTH*iseg+BOUNDARY_SEGHEADER_OFFSET)

          ! Get the start and end parameter value
          dcurrentpar = DsegInfo(1+istartidx)
          dendpar = dcurrentpar + DsegInfo(2+istartidx)
          dseglength = DsegInfo(2+istartidx)

          ! At least one of the IF-commands in the loop will activate
          ! the exit - because of the 'dt' check above!
          if (dpar .lt. dendpar) exit

        end do
      end if

      ! Subtract the start position of the current boundary component
      ! from the parameter value to get the 'local' parameter value
      ! (0 .le. dparloc .le. length(segment))
      dparloc = dpar - dcurrentpar

    case DEFAULT
      iseg = 0
      istartidx = 0
      dcurrentpar = 0.0_DP
      dendpar = 0.0_DP
      dseglength = 0.0_DP
      dparloc = 0.0_DP
      isegtype = 0
      return
    end select

    ! Use the segment type to determine how to calculate
    ! the coordinate. Remember that the segment type is noted
    ! in the first element of the integer block of each segment!
    isegtype = IsegInfo(BOUNDARY_SEGHEADER_LENGTH*iseg+BOUNDARY_SEGHEADER_TYPE)

  end subroutine boundary_getSegmentInfo2D

  !************************************************************************

!<subroutine>

  subroutine boundary_getCoords(rboundary, iboundCompIdx, dt, dx, dy, cparType)

!<description>
  ! This routine returns for a given parameter value dt the
  ! cartesian coordinates of the point on the boundary component
  ! iboundCompIdx.
!</description>

!<input>

  ! boundary structure
  type(t_boundary), intent(in) :: rboundary

  ! index of boundary component
  integer, intent(in) :: iboundCompIdx

  ! parametric value of boundary point
  real(DP), intent(in) :: dt

  ! OPTIONAL: Type of parametrisation to use.
  ! One of the BDR_PAR_xxxx constants. If not given, BDR_PAR_01 is assumed.
  integer, intent(in), optional :: cparType

!</input>

!<output>

  ! x-coordinate of boundary point
  real(DP), intent(out) :: dx

  ! y-coordinate of boundary point
  real(DP), intent(out) :: dy

!</output>

!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_IdbleSegInfo_handles,p_IintSegInfo_handles
    integer, dimension(:), pointer :: p_IsegInfo, p_IsegCount
    real(DP), dimension(:), pointer     :: p_DsegInfo, p_DmaxPar
    integer :: cpar ! local copy of cparType
    real(DP), dimension(1) :: D1

    real(DP) :: dpar, dcurrentpar, dparloc, dphi, dendpar, dseglength
    integer :: iseg,isegtype,istartidx,icomp

    cpar = BDR_PAR_01
    if (present(cparType)) cpar = cparType

    if ((iboundCompIdx.gt.rboundary%iboundarycount).or.(iboundCompIdx.lt.0)) then
      call output_line ('iboundCompIdx out of bounds!', &
          OU_CLASS_ERROR,OU_MODE_STD,'boundary_getCoords')
      call sys_halt()
    endif

    if (dt .lt. 0.0_DP) then
      call output_line ('Negative parameter values invalid!', &
          OU_CLASS_ERROR,OU_MODE_STD,'boundary_getCoords')
      call sys_halt()
    end if

    ! Get the pointers to the segment information arrays for the current
    ! boundary component:
    call storage_getbase_int(rboundary%h_Iintdatavec_handles,p_IintSegInfo_handles)
    call storage_getbase_int(int(p_IintSegInfo_handles(iboundCompIdx)),p_IsegInfo)

    call storage_getbase_int(rboundary%h_Idbldatavec_handles,p_IdbleSegInfo_handles)
    call storage_getbase_double(int(p_IdbleSegInfo_handles(iboundCompIdx)),p_DsegInfo)

    ! Get the segment-count array and the maximum-parameter array
    call storage_getbase_int(rboundary%h_IsegCount,p_IsegCount)
    call storage_getbase_double(rboundary%h_DmaxPar,p_DmaxPar)

    ! Normalise the parameter value to the range [0,TMAX)
    dpar = dt
    call boundary_normaliseParValue2D(p_IsegCount,p_DmaxPar,iboundCompIdx,cpar,dpar)

    ! Get information about the segment that contains the point
    call boundary_getSegmentInfo2D(&
        p_IsegCount,p_DmaxPar,p_IsegInfo,p_DsegInfo,iboundCompIdx,dpar,cpar,0,&
        iseg,istartidx,dcurrentpar,dendpar,dseglength,dparloc,isegtype)

    if (dseglength .eq. 0.0_DP) dseglength = 1.0_DP ! trick to avoid div/0

    ! Use the segment type to determine how to calculate
    ! the coordinate. Remember that the segment type is noted
    ! in the first element of the integer block of each segment!

    isegtype = p_IsegInfo(BOUNDARY_SEGHEADER_LENGTH*iseg+BOUNDARY_SEGHEADER_TYPE)

    select case (isegType)

      ! case of line
    case (BOUNDARY_TYPE_LINE)

      ! As we save the parametrisation in 0-1 parametrisation,
      ! when we have length-parametrisation, we have to normalise
      ! dparloc to 0 <= dparloc <= 1.
      if (cpar .eq. BDR_PAR_LENGTH) dparloc = dparloc / dseglength

      ! Calculate the x/y coordinates from the startpoint and
      ! the unit direction vector.
      dx = p_DsegInfo(istartidx+3) + dparloc*p_DsegInfo(istartidx+5)
      dy = p_DsegInfo(istartidx+4) + dparloc*p_DsegInfo(istartidx+6)

      ! case of circle segment
    case (BOUNDARY_TYPE_CIRCLE)

      ! Rescale dparloc with the length of the arc to get a value
      ! between 0 and 1; important for sin/cos functions later.
      ! In the 0-1 parametrisation, this is already the case.
      if (cpar .eq. BDR_PAR_LENGTH) dparloc = dparloc / dseglength

      ! Get the rotation angle.
      ! Use the initial rotation angle, saved at position 7 of the double
      ! precision data block
      dphi = p_DsegInfo(istartidx+7) &
           + dparloc * (p_DsegInfo(istartidx+8)-p_DsegInfo(istartidx+7))

      ! And calculate the x/y coordinate with sin/cos; the radius is
      ! to be found in element 5 of the double precision data block!
      ! The center of the circle is at position 3/4.
      dx = p_DsegInfo(istartidx+3) + p_DsegInfo(istartidx+5)*cos(dphi)
      dy = p_DsegInfo(istartidx+4) + p_DsegInfo(istartidx+5)*sin(dphi)

      ! case of analytic expression
    case (BOUNDARY_TYPE_EXPRESSION)

      ! Rescale dparloc with the length of the arc to get a value
      ! between 0 and 1; important for sin/cos functions later.
      ! In the 0-1 parametrisation, this is already the case.
      if (cpar .eq. BDR_PAR_LENGTH) dparloc = dparloc / dseglength

      ! Get absolute position of segment
      icomp = p_IsegInfo(BOUNDARY_SEGHEADER_LENGTH*iseg+BOUNDARY_SEGHEADER_EXPRESSION)

      ! Calculate the x/y coordinates from the function parser.
      D1 = (/dparloc/)
      call fparser_evalFunction(rboundary%p_rfparser, icomp+1, D1, dx)
      call fparser_evalFunction(rboundary%p_rfparser, icomp+2, D1, dy)

    case DEFAULT
      call output_line ('Wrong segment type: isegType='//sys_siL(isegType,10), &
                        OU_CLASS_ERROR,OU_MODE_STD,'boundary_getCoords')
      call sys_halt()
    end select

  end subroutine boundary_getCoords

  !************************************************************************

!<subroutine>

  subroutine boundary_getCoords_mult(rboundary, iboundCompIdx, Dt,&
                                     Dx, Dy, cparType)

!<description>
  ! This routine returns for a given attay of parameter value dt the
  ! cartesian coordinates of the points on the boundary component
  ! iboundCompIdx.
!</description>

!<input>

  ! boundary structure
  type(t_boundary), intent(in) :: rboundary

  ! index of boundary component
  integer, intent(in) :: iboundCompIdx

  ! parametric values of boundary points
  real(DP), dimension(:), intent(in) :: Dt

  ! OPTIONAL: Type of parametrisation to use.
  ! One of the BDR_PAR_xxxx constants. If not given, BDR_PAR_01 is assumed.
  integer, intent(in), optional :: cparType

!</input>

!<output>

  ! x-coordinates of boundary points
  real(DP), dimension(:), intent(out) :: Dx

  ! y-coordinates of boundary points
  real(DP), dimension(:), intent(out) :: Dy

!</output>

!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_IdbleSegInfo_handles,p_IintSegInfo_handles
    integer, dimension(:), pointer :: p_IsegInfo, p_IsegCount
    real(DP), dimension(:), pointer     :: p_DsegInfo, p_DmaxPar
    integer :: cpar ! local copy of cparType
    real(DP), dimension(1) :: D1

    real(DP) :: dpar, dcurrentpar, dparloc, dphi, dendpar, dseglength
    integer :: iseg,isegtype,istartidx,ipoint,icomp

    if ((size(Dx) .ne. size(Dt)) .or. (size(Dy) .ne. size(Dt))) then
      call output_line ('size(Dt) /= size(Dnx) /= size(Dny)!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'boundary_getCoords_mult')
      call sys_halt()
    end if

    cpar = BDR_PAR_01
    if (present(cparType)) cpar = cparType

    if ((iboundCompIdx.gt.rboundary%iboundarycount).or.(iboundCompIdx.lt.0)) then
      call output_line ('iboundCompIdx out of bounds!', &
          OU_CLASS_ERROR,OU_MODE_STD,'boundary_getCoords_mult')
      call sys_halt()
    endif

    ! Get the pointers to the segment information arrays for the current
    ! boundary component:
    call storage_getbase_int(rboundary%h_Iintdatavec_handles,p_IintSegInfo_handles)
    call storage_getbase_int(int(p_IintSegInfo_handles(iboundCompIdx)),p_IsegInfo)

    call storage_getbase_int(rboundary%h_Idbldatavec_handles,p_IdbleSegInfo_handles)
    call storage_getbase_double(int(p_IdbleSegInfo_handles(iboundCompIdx)),p_DsegInfo)

    ! Get the segment-count array and the maximum-parameter array
    call storage_getbase_int(rboundary%h_IsegCount,p_IsegCount)
    call storage_getbase_double(rboundary%h_DmaxPar,p_DmaxPar)

    ! Process each parameter value individually
    do ipoint = 1, size(Dt)

      dpar = Dt(ipoint)
      if (dpar .lt. 0.0_DP) then
        call output_line ('Negative parameter values invalid!', &
            OU_CLASS_ERROR,OU_MODE_STD,'boundary_getCoords_mult')
        call sys_halt()
      end if

      ! Normalise the parameter value to the range [0,TMAX)
      call boundary_normaliseParValue2D(p_IsegCount,p_DmaxPar,iboundCompIdx,cpar,dpar)

      ! Get information about the segment that contains the point
      call boundary_getSegmentInfo2D(&
          p_IsegCount,p_DmaxPar,p_IsegInfo,p_DsegInfo,iboundCompIdx,dpar,cpar,0,&
          iseg,istartidx,dcurrentpar,dendpar,dseglength,dparloc,isegtype)

      if (dseglength .eq. 0.0_DP) dseglength = 1.0_DP ! trick to avoid div/0

      ! Use the segment type to determine how to calculate
      ! the coordinate. Remember that the segment type is noted
      ! in the first element of the integer block of each segment!

      isegtype = p_IsegInfo(BOUNDARY_SEGHEADER_LENGTH*iseg+BOUNDARY_SEGHEADER_TYPE)

      select case (isegType)

        ! case of line
      case (BOUNDARY_TYPE_LINE)

        ! As we save the parametrisation in 0-1 parametrisation,
        ! when we have length-parametrisation, we have to normalise
        ! dparloc to 0 <= dparloc <= 1.
        if (cpar .eq. BDR_PAR_LENGTH) dparloc = dparloc / dseglength

        ! Calculate the x/y coordinates from the startpoint and
        ! the unit direction vector.
        Dx(ipoint) = p_DsegInfo(istartidx+3) + dparloc*p_DsegInfo(istartidx+5)
        Dy(ipoint) = p_DsegInfo(istartidx+4) + dparloc*p_DsegInfo(istartidx+6)

        ! case of circle segment
      case (BOUNDARY_TYPE_CIRCLE)

        ! Rescale dparloc with the length of the arc to get a value
        ! between 0 and 1; important for sin/cos functions later.
        ! In the 0-1 parametrisation, this is already the case.
        if (cpar .eq. BDR_PAR_LENGTH) dparloc = dparloc / dseglength

        ! Get the rotation angle.
        ! Use the initial rotation angle, saved at position 7 of the double
        ! precision data block
        dphi = p_DsegInfo(istartidx+7) &
             + dparloc * (p_DsegInfo(istartidx+8)-p_DsegInfo(istartidx+7))

        ! And calculate the x/y coordinate with sin/cos; the radius is
        ! to be found in element 5 of the double precision data block!
        ! The center of the circle is at position 3/4.
        Dx(ipoint) = p_DsegInfo(istartidx+3) + p_DsegInfo(istartidx+5)*cos(dphi)
        Dy(ipoint) = p_DsegInfo(istartidx+4) + p_DsegInfo(istartidx+5)*sin(dphi)

        ! case of analytic expression
      case (BOUNDARY_TYPE_EXPRESSION)

        ! Rescale dparloc with the length of the arc to get a value
        ! between 0 and 1; important for sin/cos functions later.
        ! In the 0-1 parametrisation, this is already the case.
        if (cpar .eq. BDR_PAR_LENGTH) dparloc = dparloc / dseglength
        
        ! Get absolute position of segment
        icomp = p_IsegInfo(BOUNDARY_SEGHEADER_LENGTH*iseg+BOUNDARY_SEGHEADER_EXPRESSION)
        
        ! Calculate the x/y coordinates from the function parser.
        D1 = (/dparloc/)
        call fparser_evalFunction(rboundary%p_rfparser, icomp+1, D1, Dx(ipoint))
        call fparser_evalFunction(rboundary%p_rfparser, icomp+2, D1, Dy(ipoint))

      case DEFAULT
        call output_line ('Wrong segment type: isegType='//sys_siL(isegType,10), &
                          OU_CLASS_ERROR,OU_MODE_STD,'boundary_getCoords_mult')
        call sys_halt()
      end select
    end do

  end subroutine boundary_getCoords_mult

  !************************************************************************

!<subroutine>

  subroutine boundary_getCoords_sim(rboundary, iboundCompIdx, Dt,&
                                    Dx, Dy, cparType)

!<description>
  ! This routine returns for a given attay of parameter value dt the
  ! cartesian coordinates of the points on the boundary component
  ! iboundCompIdx.
!</description>

!<input>

  ! boundary structure
  type(t_boundary), intent(in) :: rboundary

  ! index of boundary component
  integer, intent(in) :: iboundCompIdx

  ! parametric values of boundary points
  real(DP), dimension(:,:), intent(in) :: Dt

  ! OPTIONAL: Type of parametrisation to use.
  ! One of the BDR_PAR_xxxx constants. If not given, BDR_PAR_01 is assumed.
  integer, intent(in), optional :: cparType

!</input>

!<output>

  ! x-coordinates of boundary points
  real(DP), dimension(:,:), intent(out) :: Dx

  ! y-coordinates of boundary points
  real(DP), dimension(:,:), intent(out) :: Dy

!</output>

!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_IdbleSegInfo_handles,p_IintSegInfo_handles
    integer, dimension(:), pointer :: p_IsegInfo, p_IsegCount
    real(DP), dimension(:), pointer     :: p_DsegInfo, p_DmaxPar
    integer :: cpar ! local copy of cparType
    real(DP), dimension(1) :: D1

    real(DP) :: dpar, dcurrentpar, dparloc, dphi, dendpar, dseglength
    integer :: iseg,isegtype,istartidx,ipoint,iel,icomp

    if (any(shape(Dx) .ne. shape(Dt)) .or. any(shape(Dy) .ne. shape(Dt))) then
      call output_line ('size(Dt) /= size(Dnx) /= size(Dny)!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'boundary_getCoords_sim')
      call sys_halt()
    end if

    cpar = BDR_PAR_01
    if (present(cparType)) cpar = cparType

    if ((iboundCompIdx.gt.rboundary%iboundarycount).or.(iboundCompIdx.lt.0)) then
      call output_line ('iboundCompIdx out of bounds!', &
          OU_CLASS_ERROR,OU_MODE_STD,'boundary_getCoords_sim')
      call sys_halt()
    endif

    ! Get the pointers to the segment information arrays for the current
    ! boundary component:
    call storage_getbase_int(rboundary%h_Iintdatavec_handles,p_IintSegInfo_handles)
    call storage_getbase_int(int(p_IintSegInfo_handles(iboundCompIdx)),p_IsegInfo)

    call storage_getbase_int(rboundary%h_Idbldatavec_handles,p_IdbleSegInfo_handles)
    call storage_getbase_double(int(p_IdbleSegInfo_handles(iboundCompIdx)),p_DsegInfo)

    ! Get the segment-count array and the maximum-parameter array
    call storage_getbase_int(rboundary%h_IsegCount,p_IsegCount)
    call storage_getbase_double(rboundary%h_DmaxPar,p_DmaxPar)

    ! Process each parameter value individually
    do iel = 1, size(Dt,2)
      do ipoint = 1, size(Dt,1)

        dpar = Dt(ipoint,iel)
        if (dpar .lt. 0.0_DP) then
          call output_line ('Negative parameter values invalid!', &
              OU_CLASS_ERROR,OU_MODE_STD,'boundary_getCoords_sim')
          call sys_halt()
        end if

        ! Normalise the parameter value to the range [0,TMAX)
        call boundary_normaliseParValue2D(p_IsegCount,p_DmaxPar,iboundCompIdx,cpar,dpar)

        ! Get information about the segment that contains the point
        call boundary_getSegmentInfo2D(&
            p_IsegCount,p_DmaxPar,p_IsegInfo,p_DsegInfo,iboundCompIdx,dpar,cpar,0,&
            iseg,istartidx,dcurrentpar,dendpar,dseglength,dparloc,isegtype)

        if (dseglength .eq. 0.0_DP) dseglength = 1.0_DP ! trick to avoid div/0

        ! Use the segment type to determine how to calculate
        ! the coordinate. Remember that the segment type is noted
        ! in the first element of the integer block of each segment!

        isegtype = p_IsegInfo(BOUNDARY_SEGHEADER_LENGTH*iseg+BOUNDARY_SEGHEADER_TYPE)

        select case (isegType)

          ! case of line
        case (BOUNDARY_TYPE_LINE)

          ! As we save the parametrisation in 0-1 parametrisation,
          ! when we have length-parametrisation, we have to normalise
          ! dparloc to 0 <= dparloc <= 1.
          if (cpar .eq. BDR_PAR_LENGTH) dparloc = dparloc / dseglength

          ! Calculate the x/y coordinates from the startpoint and
          ! the unit direction vector.
          Dx(ipoint,iel) = p_DsegInfo(istartidx+3) + dparloc*p_DsegInfo(istartidx+5)
          Dy(ipoint,iel) = p_DsegInfo(istartidx+4) + dparloc*p_DsegInfo(istartidx+6)

          ! case of circle segment
        case (BOUNDARY_TYPE_CIRCLE)

          ! Rescale dparloc with the length of the arc to get a value
          ! between 0 and 1; important for sin/cos functions later.
          ! In the 0-1 parametrisation, this is already the case.
          if (cpar .eq. BDR_PAR_LENGTH) dparloc = dparloc / dseglength

          ! Get the rotation angle.
          ! Use the initial rotation angle, saved at position 7 of the double
          ! precision data block
          dphi = p_DsegInfo(istartidx+7) &
               + dparloc * (p_DsegInfo(istartidx+8)-p_DsegInfo(istartidx+7))

          ! And calculate the x/y coordinate with sin/cos; the radius is
          ! to be found in element 5 of the double precision data block!
          ! The center of the circle is at position 3/4.
          Dx(ipoint,iel) = p_DsegInfo(istartidx+3) + p_DsegInfo(istartidx+5)*cos(dphi)
          Dy(ipoint,iel) = p_DsegInfo(istartidx+4) + p_DsegInfo(istartidx+5)*sin(dphi)

          ! case of analytic expression
        case (BOUNDARY_TYPE_EXPRESSION)
          
          ! Rescale dparloc with the length of the arc to get a value
          ! between 0 and 1; important for sin/cos functions later.
          ! In the 0-1 parametrisation, this is already the case.
          if (cpar .eq. BDR_PAR_LENGTH) dparloc = dparloc / dseglength
          
          ! Get absolute position of segment
          icomp = p_IsegInfo(BOUNDARY_SEGHEADER_LENGTH*iseg+BOUNDARY_SEGHEADER_EXPRESSION)
          
          ! Calculate the x/y coordinates from the function parser.
          D1 = (/dparloc/)
          call fparser_evalFunction(rboundary%p_rfparser, icomp+1, D1, Dx(ipoint,iel))
          call fparser_evalFunction(rboundary%p_rfparser, icomp+2, D1, Dy(ipoint,iel))
          
        case DEFAULT
          call output_line ('Wrong segment type: isegType='//sys_siL(isegType,10), &
                            OU_CLASS_ERROR,OU_MODE_STD,'boundary_getCoords_sim')
          call sys_halt()
        end select
      end do
    end do

  end subroutine boundary_getCoords_sim

  !************************************************************************

!<function>

  function boundary_convertParameter(rboundary, iboundCompIdx, dt, &
                                     cparTypeSource, cparTypeDest) result(dresult)

!<description>
  ! This function allows to convert a parameter value dt from 0-1
  ! parametrisation to length parametrisation and back.
  ! cparTypeSource specifies the type of parametrisation of dt.
  ! cparTypeDest specifies the destination type, dt should be converted
  ! to. The return value is the converted parameter value.
!</description>

!<input>

  ! boundary structure
  type(t_boundary), intent(in) :: rboundary

  ! index of boundary component
  integer, intent(in) :: iboundCompIdx

  ! parameter value of boundary point
  real(DP), intent(in) :: dt

  ! Type of parametrisation of DT.
  ! One of the BDR_PAR_xxxx constants.
  integer, intent(in) :: cparTypeSource

  ! Type of parametrisation, DT should be converted to.
  ! One of the BDR_PAR_xxxx constants.
  integer, intent(in) :: cparTypeDest

!</input>

!<result>
  ! The converted parameter value.
  real(DP) dresult
!</result>

!</function>

    ! local variables
    integer, dimension(:), pointer :: p_IdbleSegInfo_handles,p_IintSegInfo_handles
    integer, dimension(:), pointer :: p_IsegInfo, p_IsegCount
    real(DP), dimension(:), pointer :: p_DsegInfo, p_DmaxPar

    real(DP) :: dpar, dcurrentpar, dendpar, dparloc, dseglength, dtmax
    integer :: iseg,istartidx

    ! Small check
    if ((iboundCompIdx.gt.rboundary%iboundarycount).or.(iboundCompIdx.lt.0)) then
      call output_line ('iboundCompIdx out of bounds!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'boundary_convertParameter')
      call sys_halt()
    endif

    ! In case, source and destination type is the same, it is easy.
    if (cparTypeSource .eq. cparTypeDest) then
      dresult = dt
      return
    end if

    ! Get the pointers to the segment information arrays for the current
    ! boundary component:
    call storage_getbase_int(rboundary%h_Iintdatavec_handles,p_IintSegInfo_handles)
    call storage_getbase_int(int(p_IintSegInfo_handles(iboundCompIdx)),p_IsegInfo)

    call storage_getbase_int(rboundary%h_Idbldatavec_handles,p_IdbleSegInfo_handles)
    call storage_getbase_double(int(p_IdbleSegInfo_handles(iboundCompIdx)),p_DsegInfo)

    ! Get the segment-count array and the maximum-parameter array
    call storage_getbase_int(rboundary%h_IsegCount,p_IsegCount)
    call storage_getbase_double(rboundary%h_DmaxPar,p_DmaxPar)

    ! If the parameter value exceeds the parameter interval on the boundary
    ! component, truncate the parameter value!
    select case (cparTypeSource)
    case (BDR_PAR_01)
      dtmax = real(p_IsegCount(iboundCompIdx),DP)
    case (BDR_PAR_LENGTH)
      dtmax = p_DmaxPar(iboundCompIdx)
    end select
    dpar = mod(dt,dtmax)

    ! Find segment iseg the parameter value belongs to.
    ! Remember that in the first element in the double precision block of
    ! each segment, the length of the segment is noted!

    dcurrentpar = 0.0_DP

    ! Determine the segment
    select case (cparTypeSource)
    case (BDR_PAR_01)

      ! Easy case: 0-1 parametrisation
      iseg = aint(dpar)

      ! Determine Start index of the segment in the double-prec. block
      istartidx = p_IsegInfo(BOUNDARY_SEGHEADER_LENGTH*iseg+BOUNDARY_SEGHEADER_OFFSET)

      ! Get the segment length for later use
      dseglength  = p_DsegInfo(2+istartidx)

      ! Start parameter value in length parametrisation
      dcurrentpar = p_DsegInfo(1+istartidx)

      ! Subtract the start position of the current boundary component
      ! from the parameter value to get the 'local' parameter value
      ! (0 .le. dparloc .le. 1.0)
      dparloc = dpar - real(iseg,DP)

    case (BDR_PAR_LENGTH)

      ! In the length-parametrisation, we have to search.
      do iseg = 0,p_IsegCount(iboundCompIdx)-1

        ! Determine Start index of the segment in the double-prec. block
        istartidx = p_IsegInfo(BOUNDARY_SEGHEADER_LENGTH*iseg+BOUNDARY_SEGHEADER_OFFSET)

        ! Get the start and end parameter value
        dcurrentpar = p_DsegInfo(1+istartidx)
        dseglength = p_DsegInfo(2+istartidx)
        dendpar = dcurrentpar + dseglength

        ! At least one of the IF-commands in the loop will activate
        ! the exit - because of the 'dt' check above!
        if (dpar .lt. dendpar) exit

      end do

      ! Subtract the start position of the current boundary component
      ! from the parameter value to get the 'local' parameter value
      ! (0 .le. dparloc .le. length(segment))
      dparloc = dpar - dcurrentpar

      if (dseglength .eq. 0.0_DP) dseglength = 1.0_DP ! trick to avoid div/0

      ! Divide by the segment length to get the local parameter value
      ! in 0-1 parametrisation.
      dparloc = dparloc / dseglength

    end select

    ! The local parameter value dparloc is now always in the range 0..1.
    !
    ! How shoule we convert?
    select case (cparTypeDest)
    case (BDR_PAR_01)
      ! Convert to 0-1.
      ! Take the number of teh segment as basis and add
      ! the local parameter value to get the 0-1 parameter value.
      dresult = real(iseg,DP) + dparloc

    case (BDR_PAR_LENGTH)
      ! Convert to length parametrisation. dparloc gives us the local
      ! parameter value in 0-1 parametrisation. Interpolate
      ! linearly.
      dresult = dcurrentpar + dparloc*dseglength

    end select

  end function boundary_convertParameter

  !************************************************************************

!<subroutine>

  subroutine boundary_convertParameter_mult(rboundary, iboundCompIdx, DparSource,&
                                            DparDest, cparTypeSource, cparTypeDest)

!<description>
  ! This function allows to convert an array of parameter values dt from 0-1
  ! parametrisation to length parametrisation and back.
  ! cparTypeSource specifies the type of parametrisation of DparSource.
  ! cparTypeDest specifies the destination type, dt should be converted
  ! to. The return value is the converted array.
!</description>

!<input>

  ! boundary structure
  type(t_boundary), intent(in) :: rboundary

  ! index of boundary component
  integer, intent(in) :: iboundCompIdx

  ! An array with parameter values of boundary points that should be converted.
  ! Parameter values < 0 are not converted. Parameter values > maximum parameter
  ! value are mapped to the range 0..maximum parameter value.
  real(DP), dimension(:), intent(in) :: DparSource

  ! Type of parametrisation of DparSource.
  ! One of the BDR_PAR_xxxx constants.
  integer, intent(in) :: cparTypeSource

  ! Type of parametrisation, DT should be converted to.
  ! One of the BDR_PAR_xxxx constants.
  integer, intent(in) :: cparTypeDest

!</input>

!<inputoutput>
  ! Destinatino array where to write the converted parameter values to.
  ! May coincide with DparSource.
  real(DP), dimension(:), intent(inout) :: DparDest
!</inputoutput>

!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_IdbleSegInfo_handles,p_IintSegInfo_handles
    integer, dimension(:), pointer :: p_IsegInfo, p_IsegCount
    real(DP), dimension(:), pointer     :: p_DsegInfo, p_DmaxPar

    real(DP) :: dpar, dcurrentpar, dendpar, dparloc, dseglength, dtmax
    integer :: iseg,istartidx,ipoint

    ! Small check
    if ((iboundCompIdx.gt.rboundary%iboundarycount).or.(iboundCompIdx.lt.0)) then
      call output_line ('iboundCompIdx out of bounds!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'boundary_convertParameter_mult')
      call sys_halt()
    endif

    ! In case, source and destination type is the same, it is easy.
    if (cparTypeSource .eq. cparTypeDest) then
      DparDest(:) = DparSource(:)
      return
    end if

    ! Get the pointers to the segment information arrays for the current
    ! boundary component:
    call storage_getbase_int(rboundary%h_Iintdatavec_handles,p_IintSegInfo_handles)
    call storage_getbase_int(int(p_IintSegInfo_handles(iboundCompIdx)),p_IsegInfo)

    call storage_getbase_int(rboundary%h_Idbldatavec_handles,p_IdbleSegInfo_handles)
    call storage_getbase_double(int(p_IdbleSegInfo_handles(iboundCompIdx)),p_DsegInfo)

    ! Get the segment-count array and the maximum-parameter array
    call storage_getbase_int(rboundary%h_IsegCount,p_IsegCount)
    call storage_getbase_double(rboundary%h_DmaxPar,p_DmaxPar)

    ! If the parameter value exceeds the parameter interval on the boundary
    ! component, truncate the parameter value!
    select case (cparTypeSource)
    case (BDR_PAR_01)
      dtmax = real(p_IsegCount(iboundCompIdx),DP)
    case (BDR_PAR_LENGTH)
      dtmax = p_DmaxPar(iboundCompIdx)
    end select

    do ipoint = 1,size(DparSource)

      ! Do not convert parameter values < 0.
      if (DparSource(ipoint) .ge. 0.0_DP) then

        dpar = mod(DparSource(ipoint),dtmax)

        ! Find segment iseg the parameter values belong to.
        ! Remember that in the first element in the double precision block of
        ! each segment, the length of the segment is noted!

        dcurrentpar = 0.0_DP

        ! Determine the segment
        select case (cparTypeSource)
        case (BDR_PAR_01)

          ! Easy case: 0-1 parametrisation
          iseg = aint(dpar)

          ! Determine Start index of the segment in the double-prec. block
          istartidx = p_IsegInfo(BOUNDARY_SEGHEADER_LENGTH*iseg+BOUNDARY_SEGHEADER_OFFSET)

          ! Get the segment length for later use
          dseglength  = p_DsegInfo(2+istartidx)

          ! Start parameter value in length parametrisation
          dcurrentpar = p_DsegInfo(1+istartidx)

          ! Subtract the start position of the current boundary component
          ! from the parameter value to get the 'local' parameter value
          ! (0 .le. dparloc .le. 1.0)
          dparloc = dpar - real(iseg,DP)

        case (BDR_PAR_LENGTH)

          ! In the length-parametrisation, we have to search.
          do iseg = 0,p_IsegCount(iboundCompIdx)-1

            ! Determine Start index of the segment in the double-prec. block
            istartidx = p_IsegInfo(BOUNDARY_SEGHEADER_LENGTH*iseg+BOUNDARY_SEGHEADER_OFFSET)

            ! Get the start and end parameter value
            dcurrentpar = p_DsegInfo(1+istartidx)
            dseglength = p_DsegInfo(2+istartidx)
            dendpar = dcurrentpar + dseglength

            ! At least one of the IF-commands in the loop will activate
            ! the exit - because of the 'dt' check above!
            if (dpar .lt. dendpar) exit

          end do

          ! Subtract the start position of the current boundary component
          ! from the parameter value to get the 'local' parameter value
          ! (0 .le. dparloc .le. length(segment))
          dparloc = dpar - dcurrentpar

          if (dseglength .eq. 0.0_DP) dseglength = 1.0_DP ! trick to avoid div/0

          ! Divide by the segment length to get the local parameter value
          ! in 0-1 parametrisation.
          dparloc = dparloc / dseglength

        end select

        ! The local parameter value dparloc is now always in the range 0..1.
        !
        ! How shoule we convert?
        select case (cparTypeDest)
        case (BDR_PAR_01)
          ! Convert to 0-1.
          ! Take the number of teh segment as basis and add
          ! the local parameter value to get the 0-1 parameter value.
          DparDest(ipoint) = real(iseg,DP) + dparloc

        case (BDR_PAR_LENGTH)
          ! Convert to length parametrisation. dparloc gives us the local
          ! parameter value in 0-1 parametrisation. Interpolate
          ! linearly.
          DparDest(ipoint) = dcurrentpar + dparloc*dseglength

        end select

      end if

    end do !ipoint

  end subroutine boundary_convertParameter_mult

  !************************************************************************

!<subroutine>

  subroutine boundary_convertParameter_sim(rboundary, iboundCompIdx, DparSource,&
                                           DparDest, cparTypeSource, cparTypeDest)

!<description>
  ! This function allows to convert an array of parameter values dt from 0-1
  ! parametrisation to length parametrisation and back.
  ! cparTypeSource specifies the type of parametrisation of DparSource.
  ! cparTypeDest specifies the destination type, dt should be converted
  ! to. The return value is the converted array.
!</description>

!<input>

  ! boundary structure
  type(t_boundary), intent(in) :: rboundary

  ! index of boundary component
  integer, intent(in) :: iboundCompIdx

  ! An array with parameter values of boundary points that should be converted.
  ! Parameter values < 0 are not converted. Parameter values > maximum parameter
  ! value are mapped to the range 0..maximum parameter value.
  real(DP), dimension(:,:), intent(in) :: DparSource

  ! Type of parametrisation of DparSource.
  ! One of the BDR_PAR_xxxx constants.
  integer, intent(in) :: cparTypeSource

  ! Type of parametrisation, DT should be converted to.
  ! One of the BDR_PAR_xxxx constants.
  integer, intent(in) :: cparTypeDest

!</input>

!<inputoutput>
  ! Destinatino array where to write the converted parameter values to.
  ! May coincide with DparSource.
  real(DP), dimension(:,:), intent(inout) :: DparDest
!</inputoutput>

!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_IdbleSegInfo_handles,p_IintSegInfo_handles
    integer, dimension(:), pointer :: p_IsegInfo, p_IsegCount
    real(DP), dimension(:), pointer     :: p_DsegInfo, p_DmaxPar

    real(DP) :: dpar, dcurrentpar, dendpar, dparloc, dseglength, dtmax
    integer :: iseg,istartidx,ipoint,iel

    ! Small check
    if ((iboundCompIdx.gt.rboundary%iboundarycount).or.(iboundCompIdx.lt.0)) then
      call output_line ('iboundCompIdx out of bounds!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'boundary_convertParameter_sim')
      call sys_halt()
    endif

    ! In case, source and destination type is the same, it is easy.
    if (cparTypeSource .eq. cparTypeDest) then
      DparDest(:,:) = DparSource(:,:)
      return
    end if

    ! Get the pointers to the segment information arrays for the current
    ! boundary component:
    call storage_getbase_int(rboundary%h_Iintdatavec_handles,p_IintSegInfo_handles)
    call storage_getbase_int(int(p_IintSegInfo_handles(iboundCompIdx)),p_IsegInfo)

    call storage_getbase_int(rboundary%h_Idbldatavec_handles,p_IdbleSegInfo_handles)
    call storage_getbase_double(int(p_IdbleSegInfo_handles(iboundCompIdx)),p_DsegInfo)

    ! Get the segment-count array and the maximum-parameter array
    call storage_getbase_int(rboundary%h_IsegCount,p_IsegCount)
    call storage_getbase_double(rboundary%h_DmaxPar,p_DmaxPar)

    ! If the parameter value exceeds the parameter interval on the boundary
    ! component, truncate the parameter value!
    select case (cparTypeSource)
    case (BDR_PAR_01)
      dtmax = real(p_IsegCount(iboundCompIdx),DP)
    case (BDR_PAR_LENGTH)
      dtmax = p_DmaxPar(iboundCompIdx)
    end select

    do iel = 1, size(DparSource,2)
      do ipoint = 1, size(DparSource,1)

        ! Do not convert parameter values < 0.
        if (DparSource(ipoint,iel) .ge. 0.0_DP) then

          dpar = mod(DparSource(ipoint,iel),dtmax)

          ! Find segment iseg the parameter values belong to.
          ! Remember that in the first element in the double precision block of
          ! each segment, the length of the segment is noted!

          dcurrentpar = 0.0_DP

          ! Determine the segment
          select case (cparTypeSource)
          case (BDR_PAR_01)

            ! Easy case: 0-1 parametrisation
            iseg = aint(dpar)

            ! Determine Start index of the segment in the double-prec. block
            istartidx = p_IsegInfo(BOUNDARY_SEGHEADER_LENGTH*iseg+BOUNDARY_SEGHEADER_OFFSET)

            ! Get the segment length for later use
            dseglength  = p_DsegInfo(2+istartidx)

            ! Start parameter value in length parametrisation
            dcurrentpar = p_DsegInfo(1+istartidx)

            ! Subtract the start position of the current boundary component
            ! from the parameter value to get the 'local' parameter value
            ! (0 .le. dparloc .le. 1.0)
            dparloc = dpar - real(iseg,DP)

          case (BDR_PAR_LENGTH)

            ! In the length-parametrisation, we have to search.
            do iseg = 0,p_IsegCount(iboundCompIdx)-1

              ! Determine Start index of the segment in the double-prec. block
              istartidx = p_IsegInfo(BOUNDARY_SEGHEADER_LENGTH*iseg+BOUNDARY_SEGHEADER_OFFSET)

              ! Get the start and end parameter value
              dcurrentpar = p_DsegInfo(1+istartidx)
              dseglength = p_DsegInfo(2+istartidx)
              dendpar = dcurrentpar + dseglength

              ! At least one of the IF-commands in the loop will activate
              ! the exit - because of the 'dt' check above!
              if (dpar .lt. dendpar) exit

            end do

            ! Subtract the start position of the current boundary component
            ! from the parameter value to get the 'local' parameter value
            ! (0 .le. dparloc .le. length(segment))
            dparloc = dpar - dcurrentpar

            if (dseglength .eq. 0.0_DP) dseglength = 1.0_DP ! trick to avoid div/0

            ! Divide by the segment length to get the local parameter value
            ! in 0-1 parametrisation.
            dparloc = dparloc / dseglength

          end select

          ! The local parameter value dparloc is now always in the range 0..1.
          !
          ! How shoule we convert?
          select case (cparTypeDest)
          case (BDR_PAR_01)
            ! Convert to 0-1.
            ! Take the number of teh segment as basis and add
            ! the local parameter value to get the 0-1 parameter value.
            DparDest(ipoint,iel) = real(iseg,DP) + dparloc

          case (BDR_PAR_LENGTH)
            ! Convert to length parametrisation. dparloc gives us the local
            ! parameter value in 0-1 parametrisation. Interpolate
            ! linearly.
            DparDest(ipoint,iel) = dcurrentpar + dparloc*dseglength

          end select

        end if

      end do ! ipoint
    end do ! iel

  end subroutine boundary_convertParameter_sim

  !************************************************************************

!<subroutine>

  subroutine boundary_createRegion (rboundary, iboundCompIdx, iboundSegIdx, &
                                    rregion, cpartype)

!<description>
  ! This routine creates a boundary region from a boundary segment.
  ! Where the boundary segment is saved internally, the boundary region
  ! is a parametrised description of the boundary segment, i.e.
  ! giving the caller information about the minimum and maximum parameter
  ! value etc.
  ! In the standard setting, the startpoint belongs to the boundary
  ! segment, the endpoint not.
!</description>

!<input>

  ! boundary structure
  type(t_boundary), intent(in) :: rboundary

  ! Index of boundary component.
  integer, intent(in) :: iboundCompIdx

  ! Index of the boundary segment.
  ! =0: Create a boundary region that covers the whole boundary component.
  integer, intent(in) :: iboundSegIdx

  ! OPTIONAL: Type of parametrisation to use.
  ! One of the BDR_PAR_xxxx constants. If not given, BDR_PAR_01 is assumed.
  integer, intent(in), optional :: cparType

!</input>

!<output>

  ! Boundary region that is characterised by the boundary segment
  type(t_boundaryRegion), intent(out) :: rregion

!</output>

!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_IdbleSegInfo_handles,p_IintSegInfo_handles
    integer, dimension(:), pointer :: p_IsegInfo, p_IsegCount
    real(DP), dimension(:), pointer :: p_DsegInfo, p_DmaxPar

    real(DP) :: dcurrentpar, dendpar, dmaxpar
    integer :: istartidx

    integer :: cpar ! local copy of cparType

    cpar = BDR_PAR_01
    if (present(cparType)) cpar = cparType

    if ((iboundCompIdx .gt. rboundary%iboundarycount) .or. (iboundCompIdx.lt.0)) then
      call output_line ('iboundCompIdx out of bounds!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'boundary_createRegion')
      call sys_halt()
    endif

    ! Get the pointers to the segment information arrays for the current
    ! boundary component:
    call storage_getbase_int(rboundary%h_Iintdatavec_handles,p_IintSegInfo_handles)
    call storage_getbase_int(int(p_IintSegInfo_handles(iboundCompIdx)),p_IsegInfo)

    call storage_getbase_int(rboundary%h_Idbldatavec_handles,p_IdbleSegInfo_handles)
    call storage_getbase_double(int(p_IdbleSegInfo_handles(iboundCompIdx)),p_DsegInfo)

    ! Get the segment-count array and the maximum-parameter array
    call storage_getbase_int(rboundary%h_IsegCount,p_IsegCount)
    call storage_getbase_double(rboundary%h_DmaxPar,p_DmaxPar)

    if (iboundSegIdx .ne. 0) then

      if ((iboundSegIdx .gt. p_IsegCount(iboundCompIdx)) .or. (iboundSegIdx.lt.0)) then
        call output_line ('iboundSegIdx out of bounds!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'boundary_createRegion')
        call sys_halt()
      endif

      ! Find segment iseg the parameter value belongs to.
      ! Remember that in the first element in the double precision block of
      ! each segment, the length of the segment is noted!

      istartidx = p_IsegInfo(BOUNDARY_SEGHEADER_LENGTH*0+BOUNDARY_SEGHEADER_OFFSET)
      dcurrentpar = 0.0_DP

      ! Determine Start index of the segment in the double-prec. block
      istartidx = p_IsegInfo(BOUNDARY_SEGHEADER_LENGTH*(iboundSegIdx-1)+BOUNDARY_SEGHEADER_OFFSET)

      ! Get the start and end parameter value - depending on the parametrisation
      select case (cpar)
      case (BDR_PAR_01)
        dcurrentpar = real(iboundSegIdx-1,DP)
        dendpar     = real(iboundSegIdx-1+1,DP)
        dmaxpar     = real(p_IsegCount(iboundCompIdx),DP)
      case (BDR_PAR_LENGTH)
        dcurrentpar = p_DsegInfo(1+istartidx)
        dendpar     = dcurrentpar + p_DsegInfo(2+istartidx)
        dmaxpar     = p_DmaxPar(iboundCompIdx)
      end select

      ! Set the segment type in the boundary region structure
      rregion%isegmentType = p_IsegInfo(BOUNDARY_SEGHEADER_LENGTH*(iboundSegIdx-1)+BOUNDARY_SEGHEADER_TYPE)

    else

      ! Create a boundary region that covers the whole boundary component.
      select case (cpar)
      case (BDR_PAR_01)
        dcurrentpar = real(iboundSegIdx-1,DP)
        dmaxpar     = real(p_IsegCount(iboundCompIdx),DP)
      case (BDR_PAR_LENGTH)
        dcurrentpar = p_DsegInfo(1+istartidx)
        dmaxpar     = p_DmaxPar(iboundCompIdx)
      end select

      dcurrentpar = 0.0_DP
      dendpar     = dmaxpar

      ! We have an unspecified boundary segment
      rregion%isegmentType = BOUNDARY_TYPE_ANALYTIC

    end if

    ! Create the boundary region structure
    rregion%cparType = cpar
    rregion%dminParam = dcurrentpar
    rregion%dmaxParam = dendpar
    rregion%ctype = BDR_TP_CURVE
    rregion%iproperties = BDR_PROP_WITHSTART
    rregion%iboundCompIdx = iboundCompIdx
    rregion%iboundSegIdx = iboundSegIdx
    rregion%dmaxParamBC = dmaxpar

  end subroutine boundary_createRegion

  ! ***************************************************************************

!<subroutine>

  subroutine boundary_convertRegion(rboundary, rregion, cparType)

!<description>
  ! This subroutine converts the parameter values of a boundary region
  ! from 0-1 parametrisation to length parametrisation and back.
!</description>

!<input>

  ! boundary structure
  type(t_boundary), intent(in) :: rboundary

  ! Type of parametrisation, the boundary region should be converted to
  ! One of the BDR_PAR_xxxx constants.
  integer, intent(in) :: cparType

!</input>

!<inputoutput>

  ! boundary region to be converted
  type(t_boundaryRegion), intent(inout) :: rregion

!</inputoutput>

!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_IsegCount
    real(DP), dimension(:), pointer :: p_DmaxPar

    real(DP) :: dminParam,dmaxParam,dmaxParamBC
    integer :: iboundCompIdx

    ! Check of we have to convert the region at all
    if (rregion%cparType .eq. cparType) return

    ! Get the segment-count array and the maximum-parameter array
    call storage_getbase_int(rboundary%h_IsegCount,p_IsegCount)
    call storage_getbase_double(rboundary%h_DmaxPar,p_DmaxPar)

    ! Get boundary component from boundary region
    iboundCompIdx = rregion%iboundCompIdx

    ! If the parameter value exceeds the parameter interval on the
    ! boundary component, truncate the parameter value!
    select case (cparType)
    case (BDR_PAR_01)
      dmaxParamBC = real(p_IsegCount(iboundCompIdx),DP)
    case (BDR_PAR_LENGTH)
      dmaxParamBC = p_DmaxPar(iboundCompIdx)
    end select

    ! Convert parametrisation of minimum parameter value
    dminParam = boundary_convertParameter(rboundary, iboundCompIdx,&
        rregion%dminParam, rregion%cparType, cparType)

    ! Convert parametrisation of maximum parameter value
    dmaxParam = boundary_convertParameter(rboundary, iboundCompIdx,&
        rregion%dmaxParam, rregion%cparType, cparType)

    ! Check if the maximum parameter valeu is smaller than the minimum
    ! parameter value and set maximum parameter value accordingly
    if (dmaxParam .lt. dminParam) dmaxParam=dmaxParamBC

    ! Convert the boundary region structure
    rregion%cparType = cparType
    rregion%dminParam = dminParam
    rregion%dmaxParam = dmaxParam
    rregion%dmaxParamBC = DmaxParamBC

  end subroutine boundary_convertRegion

  ! ***************************************************************************

!<function>

  logical function boundary_isInRegion (rregion,iboundCompIdx,dparam)

!<description>
  ! Checks whether a point given by a parameter value of a point is in a
  ! boundary region.
!</description>

!<input>
  ! A boundary region structure describing a part of a boundary
  type(t_boundaryRegion), intent(in) :: rregion

  ! The number of the boundary component of the point.
  integer, intent(in) :: iboundCompIdx

  ! The parameter value of the point to be checked.
  ! Must be in the range 0..max. par. value.
  ! Points with negative parameter values are not in the region by definition.
  ! The parametrisation type (0-1 or length parametrisation) must match
  ! the parametrisation in the boundary region rregion!
  real(DP), intent(in) :: dparam
!</input>

!<result>
  ! TRUE, if the point is inside the region,
  ! FALSE otherwise.
!</result>

!</function>

    ! local variables
    real(DP) :: dminpar,dmaxpar, dpar1,dpar2

    ! Default setting: not inside.
    boundary_isInRegion = .false.

    ! Correct boundary component?
    if (iboundCompIdx .ne. rregion%iboundCompIdx) return

    ! If the parameter value is negative, the point is surely not in the
    ! region. Parameter values must be positive!
    ! (Note: Do not change this behaviour to 'rounding' into 0..tmax,
    ! since there may be routines (like in the triangulation) thay
    ! depend on this to figure out which vertices are on the physical
    ! boundary and which not!)
    if (dparam .lt. 0.0_DP) return

    ! Does the region cover the complete boundary?
    if (rregion%dmaxParam-rregion%dminParam .gt. rregion%dmaxParamBC) then
      boundary_isInRegion = .true.
      return
    end if

    ! Region does not cover the complete boundary, but a part of it.
    ! Which part?

    ! Get the real bounds...
    dminpar = mod(rregion%dminParam,rregion%dmaxParamBC)
    dmaxpar = mod(rregion%dmaxParam,rregion%dmaxParamBC)

    ! And make sure, dmaxpar >= dminpar!
    if ((dmaxpar .le. dminpar) .and. &
        (rregion%dminParam .ne. rregion%dmaxParam)) then
      dmaxpar = dmaxpar + rregion%dmaxParamBC
    end if

    ! Get the parameter value on the boundary - in the set [0..TMAX)
    ! and [TMAX..2*TMAX).
    dpar1 = mod(dparam,rregion%dmaxParamBC)
    dpar2 = dpar1 + rregion%dmaxParamBC

    ! Check if dpar1 or dpar2 is in that region.
    ! For 'normal' boundary regions, dpar1 should be inside, dpar2 not.
    ! For boundary regions crossing the maximum parameter value,
    ! dpar2 will be inside and dpar1 not.

    if ((dpar1 .ge. dminpar) .and. &
        (dpar1 .le. dmaxpar)) then

      ! What is up with the endpoints?
      if ( (dpar1 .eq. dminpar) .and. &
          (iand(rregion%iproperties,BDR_PROP_WITHSTART) .eq. 0)) return
      if ( (dpar1 .eq. dmaxpar) .and. &
          (iand(rregion%iproperties,BDR_PROP_WITHEND) .eq. 0)) return

      ! It is inside.
      boundary_isInRegion = .true.

    else if ((dpar2 .ge. dminpar) .and. &
             (dpar2 .le. dmaxpar)) then

      ! What is up with the endpoints?
      if ( (dpar2 .eq. dminpar) .and. &
           (iand(rregion%iproperties,BDR_PROP_WITHSTART) .eq. 0)) return
      if ( (dpar2 .eq. dmaxpar) .and. &
           (iand(rregion%iproperties,BDR_PROP_WITHEND) .eq. 0)) return

      ! It is inside.
      boundary_isInRegion = .true.

    end if

  end function boundary_isInRegion

  ! ***************************************************************************

!<function>

  real(DP) function boundary_getRegionLength (rboundary,rregion) result (dlength)

!<description>
  ! Calculates the length of a part of the boundary identified by rregion
  ! on boundary rboundary.
!</description>

!<input>
  ! Boundary structure, the boundary region should refer to
  type(t_boundary), intent(in) :: rboundary

  ! The boundary reg8ion structure which length is to be computed
  type(t_boundaryRegion), intent(in) :: rregion
!</input>

!<result>
  ! The length of the boundary region rregion on the boundary rboundary.
!</result>

!</function>

    ! local variables
    real(DP) :: dlen1,dlen2

    ! Is the boundary region parametrised for the length? Then it is easy...
    if (rregion%cparType .eq. BDR_PAR_LENGTH) then
      dlength = rregion%dmaxParam - rregion%dminParam
      return
    end if

    ! Otherwise, compute the parameter values in length-parametrisation
    ! and subtract them to get the length.
    dlen1 = boundary_convertParameter(rboundary, rregion%iboundCompIdx, &
                                      rregion%dminParam, &
                                      rregion%cparType, BDR_PAR_LENGTH)

    dlen2 = boundary_convertParameter(rboundary, rregion%iboundCompIdx, &
                                      rregion%dmaxParam, &
                                      rregion%cparType, BDR_PAR_LENGTH)

    dlength = dlen2-dlen1

  end function boundary_getRegionLength

  !************************************************************************

!<subroutine>

  subroutine boundary_getNormalVec2D(rboundary, iboundCompIdx, dt, &
                                     dnx, dny, cnormalMean, cparType)

!<description>
  ! This routine returns for a given parameter value dt the
  ! the outward unit normal vector of the point on the boundary component
  ! iboundCompIdx.\\
  ! dt is truncated to the interval [0,max. par. value of iboundCompIdx).
  !
  ! The parameters cnormalMean are optional. This parameters
  ! define the behaviour of the function if dt specifies the parameter
  ! value of a corner point on the boundary where the normal vector
  ! is usually not unique. It allows to choose whether to calculate the
  ! 'right', 'left' or averaged normal vector.
!</description>

!<input>

    ! boundary structure
    type(t_boundary), intent(in) :: rboundary

    ! index of boundary component
    integer, intent(in) :: iboundCompIdx

    ! parametric value of boundary point
    real(DP), intent(in) :: dt

    ! OPTIONAL: For points with non-unique normal vectors, this decides
    ! how to calculate the normal vector. One of the BDR_NORMAL_xxxx
    ! constants.
    ! BDR_NORMAL_MEAN calculates the mean of the 'right' and 'left'
    !   normal vector. This is the standard setting
    ! BDR_NORMAL_LEFT calculates the 'left' normal vector (which arises
    !   in the limit when appoximating dt by 0->dt).
    ! BDR_NORMAL_RIGHT calculates the 'right' normal vector (which arises
    !   in the limit when approximating dt by dtmax->dt).
    integer, intent(in), optional :: cnormalMean

    ! OPTIONAL: Type of parametrisation to use.
    ! One of the BDR_PAR_xxxx constants. If not given, BDR_PAR_01 is assumed.
    integer, intent(in), optional :: cparType

!</input>

!<output>

    ! x-coordinate of normal vector
    real(DP), intent(out) :: dnx

    ! y-coordinate of normal vector
    real(DP), intent(out) :: dny

!</output>

!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_IdbleSegInfo_handles,p_IintSegInfo_handles
    integer, dimension(:), pointer :: p_IsegInfo, p_IsegCount
    real(DP), dimension(:), pointer     :: p_DsegInfo, p_DmaxPar
    integer :: cpar ! local copy of cparType
    integer :: cnormalMeanCalc

    real(DP) :: dpar, dcurrentpar, dparloc, dendpar, dseglength, dnorm, dnx0, dny0
    integer :: iseg,isegtype,istartidx

    cpar = BDR_PAR_01
    if (present(cparType)) cpar = cparType

    cnormalMeanCalc = BDR_NORMAL_MEAN
    if (present(cnormalMean)) cnormalMeanCalc = cnormalMean

    if ((iboundCompIdx.gt.rboundary%iboundarycount).or.(iboundCompIdx.lt.0)) then
      call output_line ('iboundCompIdx out of bounds!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'boundary_getNormalVec2D')
      call sys_halt()
    endif

    ! Get the pointers to the segment information arrays for the current
    ! boundary component:
    call storage_getbase_int(rboundary%h_Iintdatavec_handles,p_IintSegInfo_handles)
    call storage_getbase_int(int(p_IintSegInfo_handles(iboundCompIdx)),p_IsegInfo)

    call storage_getbase_int(rboundary%h_Idbldatavec_handles,p_IdbleSegInfo_handles)
    call storage_getbase_double(int(p_IdbleSegInfo_handles(iboundCompIdx)),p_DsegInfo)

    ! Get the segment-count array and the maximum-parameter array
    call storage_getbase_int(rboundary%h_IsegCount,p_IsegCount)
    call storage_getbase_double(rboundary%h_DmaxPar,p_DmaxPar)

    ! Normalise the parameter value to the range [0,TMAX)
    dpar = dt
    call boundary_normaliseParValue2D(p_IsegCount,p_DmaxPar,iboundCompIdx,cpar,dpar)

    ! Find segment iseg the parameter value belongs to.
    ! Remember that in the first element in the double precision block of
    ! each segment, the length of the segment is noted!
    call boundary_getSegmentInfo2D(&
      p_IsegCount,p_DmaxPar,p_IsegInfo,p_DsegInfo,iboundCompIdx,dpar,cpar,0,&
      iseg,istartidx,dcurrentpar,dendpar,dseglength,dparloc,isegtype)

    if (dseglength .eq. 0.0_DP) dseglength = 1.0_DP ! trick to avoid div/0

    ! If the local parameter value of dt is <> 0, then the normal vector can
    ! be computed uniquely at boundary point dt.
    ! Otherwise, we have to average the (two) normal vector(s) of the adjacent
    ! boundary components. Check if we are in the simple case.
    if (dparloc .ne. 0.0_DP) then

      call boundary_getNormal2D (rboundary%p_rfparser,iseg,isegType,cpar,&
          dparLoc,dseglength,istartidx,p_IsegInfo,p_DsegInfo,dnx,dny)

    else

      ! Ok, we are at the endpoint of the interval. Now, cnormalMean decides on
      ! what to do.
      dnx = 0.0_DP
      dny = 0.0_DP
      if ((cnormalMeanCalc .eq. BDR_NORMAL_MEAN) .or. &
          (cnormalMeanCalc .eq. BDR_NORMAL_RIGHT)) then

        ! Calculate the 'right' normal into dnx/dny.
        !
        ! We already have the segment 'right' to the point.
        ! Get the corresponding normal vector.

        call boundary_getNormal2D (rboundary%p_rfparser,iseg,isegType,cpar,&
            dparLoc,dseglength,istartidx,p_IsegInfo,p_DsegInfo,dnx,dny)

      end if

      if ((cnormalMeanCalc .eq. BDR_NORMAL_MEAN) .or. &
          (cnormalMeanCalc .eq. BDR_NORMAL_LEFT)) then

        ! Calculate the 'left' normal into dnx/dny.
        !
        ! Find segment iseg previous to the parameter value belongs to.
        ! Remember that in the first element in the double precision block of
        ! each segment, the length of the segment is noted!

        call boundary_getSegmentInfo2D(&
          p_IsegCount,p_DmaxPar,p_IsegInfo,p_DsegInfo,iboundCompIdx,dpar,cpar,1,&
          iseg,istartidx,dcurrentpar,dendpar,dseglength,dparloc,isegtype)

        ! Calculate the normal into dnx0/dny0
        call boundary_getNormal2D (rboundary%p_rfparser,iseg,isegType,cpar,&
            dparLoc,dseglength,istartidx,p_IsegInfo,p_DsegInfo,dnx0,dny0)

        ! and add it to dnx/dny
        dnx = dnx+dnx0
        dny = dny+dny0
      end if

      ! Normalise the vector -- in case two vectors are summed up.
      dnorm = sqrt(dnx*dnx+dny*dny)
      dnx = dnx/dnorm
      dny = dny/dnorm
    end if

  end subroutine boundary_getNormalVec2D

  !************************************************************************

!<subroutine>

  subroutine boundary_getNormalVec2D_mult(rboundary, iboundCompIdx, Dt, &
                                          Dnx, Dny, cnormalMean, cparType)

!<description>
  ! This routine returns for a given array of parameter values dt the
  ! the outward unit normal vectors of the points on the boundary component
  ! iboundCompIdx.\\
  ! Dt is truncated to the interval [0,max. par. value of iboundCompIdx).
  !
  ! The parameters cnormalMean are optional. This parameters
  ! define the behaviour of the function if Dt specifies the parameter
  ! values of a corner point on the boundary where the normal vector
  ! is usually not unique. It allows to choose whether to calculate the
  ! 'right', 'left' or averaged normal vector.
!</description>

!<input>

    ! boundary structure
    type(t_boundary), intent(in) :: rboundary

    ! index of boundary component
    integer, intent(in) :: iboundCompIdx

    ! parametric values of boundary points
    real(DP), dimension(:), intent(in) :: Dt

    ! OPTIONAL: For points with non-unique normal vectors, this decides
    ! how to calculate the normal vector. One of the BDR_NORMAL_xxxx
    ! constants.
    ! BDR_NORMAL_MEAN calculates the mean of the 'right' and 'left'
    !   normal vector. This is the standard setting
    ! BDR_NORMAL_LEFT calculates the 'left' normal vector (which arises
    !   in the limit when appoximating dt by 0->dt).
    ! BDR_NORMAL_RIGHT calculates the 'right' normal vector (which arises
    !   in the limit when approximating dt by dtmax->dt).
    integer, intent(in), optional :: cnormalMean

    ! OPTIONAL: Type of parametrisation to use.
    ! One of the BDR_PAR_xxxx constants. If not given, BDR_PAR_01 is assumed.
    integer, intent(in), optional :: cparType

!</input>

!<output>

    ! x-coordinates of normal vectors
    real(DP), dimension(:), intent(out) :: Dnx

    ! y-coordinates of normal vectors
    real(DP), dimension(:), intent(out) :: Dny

!</output>

!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_IdbleSegInfo_handles,p_IintSegInfo_handles
    integer, dimension(:), pointer :: p_IsegInfo, p_IsegCount
    real(DP), dimension(:), pointer     :: p_DsegInfo, p_DmaxPar
    integer :: cpar ! local copy of cparType
    integer :: cnormalMeanCalc

    real(DP) :: dpar, dcurrentpar, dparloc, dendpar, dseglength, dnorm, dnx0, dny0
    integer :: iseg,isegtype,istartidx,ipoint

    if ((size(Dnx) .ne. size(Dt)) .or. (size(Dny) .ne. size(Dt))) then
      call output_line ('size(Dt) /= size(Dnx) /= size(Dny)!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'boundary_getNormalVec2D_mult')
      call sys_halt()
    end if

    cpar = BDR_PAR_01
    if (present(cparType)) cpar = cparType

    cnormalMeanCalc = BDR_NORMAL_MEAN
    if (present(cnormalMean)) cnormalMeanCalc = cnormalMean

    if ((iboundCompIdx.gt.rboundary%iboundarycount).or.(iboundCompIdx.lt.0)) then
      call output_line ('iboundCompIdx out of bounds!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'boundary_getNormalVec2D_mult')
      call sys_halt()
    endif

    ! Get the pointers to the segment information arrays for the current
    ! boundary component:
    call storage_getbase_int(rboundary%h_Iintdatavec_handles,p_IintSegInfo_handles)
    call storage_getbase_int(int(p_IintSegInfo_handles(iboundCompIdx)),p_IsegInfo)

    call storage_getbase_int(rboundary%h_Idbldatavec_handles,p_IdbleSegInfo_handles)
    call storage_getbase_double(int(p_IdbleSegInfo_handles(iboundCompIdx)),p_DsegInfo)

    ! Get the segment-count array and the maximum-parameter array
    call storage_getbase_int(rboundary%h_IsegCount,p_IsegCount)
    call storage_getbase_double(rboundary%h_DmaxPar,p_DmaxPar)

    ! Process each parameter value individually
    do ipoint = 1, size(Dt)

      ! Normalise the parameter value to the range [0,TMAX)
      dpar = Dt(ipoint)
      call boundary_normaliseParValue2D(p_IsegCount,p_DmaxPar,iboundCompIdx,cpar,dpar)

      ! Find segment iseg the parameter value belongs to.
      ! Remember that in the first element in the double precision block of
      ! each segment, the length of the segment is noted!
      call boundary_getSegmentInfo2D(&
          p_IsegCount,p_DmaxPar,p_IsegInfo,p_DsegInfo,iboundCompIdx,dpar,cpar,0,&
          iseg,istartidx,dcurrentpar,dendpar,dseglength,dparloc,isegtype)

      if (dseglength .eq. 0.0_DP) dseglength = 1.0_DP ! trick to avoid div/0

      ! If the local parameter value of dt is <> 0, then the normal vector can
      ! be computed uniquely at boundary point dt.
      ! Otherwise, we have to average the (two) normal vector(s) of the adjacent
      ! boundary components. Check if we are in the simple case.
      if (dparloc .ne. 0.0_DP) then

        call boundary_getNormal2D (rboundary%p_rfparser,iseg, isegType,cpar,&
            dparLoc,dseglength,istartidx,p_IsegInfo,p_DsegInfo,Dnx(ipoint),Dny(ipoint))

      else

        ! Ok, we are at the endpoint of the interval. Now, cnormalMean decides on
        ! what to do.
        Dnx(ipoint) = 0.0_DP
        Dny(ipoint) = 0.0_DP
        if ((cnormalMeanCalc .eq. BDR_NORMAL_MEAN) .or. &
            (cnormalMeanCalc .eq. BDR_NORMAL_RIGHT)) then

          ! Calculate the 'right' normal into dnx/dny.
          !
          ! We already have the segment 'right' to the point.
          ! Get the corresponding normal vector.

          call boundary_getNormal2D (rboundary%p_rfparser,iseg,isegType,cpar,&
              dparLoc,dseglength,istartidx,p_IsegInfo,p_DsegInfo,Dnx(ipoint),Dny(ipoint))

        end if

        if ((cnormalMeanCalc .eq. BDR_NORMAL_MEAN) .or. &
            (cnormalMeanCalc .eq. BDR_NORMAL_LEFT)) then

          ! Calculate the 'left' normal into dnx/dny.
          !
          ! Find segment iseg previous to the parameter value belongs to.
          ! Remember that in the first element in the double precision block of
          ! each segment, the length of the segment is noted!

          call boundary_getSegmentInfo2D(&
              p_IsegCount,p_DmaxPar,p_IsegInfo,p_DsegInfo,iboundCompIdx,dpar,cpar,1,&
              iseg,istartidx,dcurrentpar,dendpar,dseglength,dparloc,isegtype)

          ! Calculate the normal into dnx0/dny0
          call boundary_getNormal2D (rboundary%p_rfparser,iseg,isegType,cpar,&
              dparLoc,dseglength,istartidx,p_IsegInfo,p_DsegInfo,dnx0,dny0)

          ! and add it to dnx/dny
          Dnx(ipoint) = Dnx(ipoint)+dnx0
          Dny(ipoint) = Dny(ipoint)+dny0
        end if

        ! Normalise the vector -- in case two vectors are summed up.
        dnorm = sqrt(Dnx(ipoint)*Dnx(ipoint)+Dny(ipoint)*Dny(ipoint))
        Dnx(ipoint) = Dnx(ipoint)/dnorm
        Dny(ipoint) = Dny(ipoint)/dnorm
      end if

    end do

  end subroutine boundary_getNormalVec2D_mult

  !************************************************************************

!<subroutine>

  subroutine boundary_getNormalVec2D_sim(rboundary, iboundCompIdx, Dt, &
                                         Dnx, Dny, cnormalMean, cparType)

!<description>
  ! This routine returns for a given array of parameter values dt the
  ! the outward unit normal vectors of the points on the boundary component
  ! iboundCompIdx.\\
  ! Dt is truncated to the interval [0,max. par. value of iboundCompIdx).
  !
  ! The parameters cnormalMean are optional. This parameters
  ! define the behaviour of the function if Dt specifies the parameter
  ! values of a corner point on the boundary where the normal vector
  ! is usually not unique. It allows to choose whether to calculate the
  ! 'right', 'left' or averaged normal vector.
!</description>

!<input>

    ! boundary structure
    type(t_boundary), intent(in) :: rboundary

    ! index of boundary component
    integer, intent(in) :: iboundCompIdx

    ! parametric values of boundary points
    real(DP), dimension(:,:), intent(in) :: Dt

    ! OPTIONAL: For points with non-unique normal vectors, this decides
    ! how to calculate the normal vector. One of the BDR_NORMAL_xxxx
    ! constants.
    ! BDR_NORMAL_MEAN calculates the mean of the 'right' and 'left'
    !   normal vector. This is the standard setting
    ! BDR_NORMAL_LEFT calculates the 'left' normal vector (which arises
    !   in the limit when appoximating dt by 0->dt).
    ! BDR_NORMAL_RIGHT calculates the 'right' normal vector (which arises
    !   in the limit when approximating dt by dtmax->dt).
    integer, intent(in), optional :: cnormalMean

    ! OPTIONAL: Type of parametrisation to use.
    ! One of the BDR_PAR_xxxx constants. If not given, BDR_PAR_01 is assumed.
    integer, intent(in), optional :: cparType

!</input>

!<output>

    ! x-coordinates of normal vectors
    real(DP), dimension(:,:), intent(out) :: Dnx

    ! y-coordinates of normal vectors
    real(DP), dimension(:,:), intent(out) :: Dny

!</output>

!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_IdbleSegInfo_handles,p_IintSegInfo_handles
    integer, dimension(:), pointer :: p_IsegInfo, p_IsegCount
    real(DP), dimension(:), pointer     :: p_DsegInfo, p_DmaxPar
    integer :: cpar ! local copy of cparType
    integer :: cnormalMeanCalc

    real(DP) :: dpar, dcurrentpar, dparloc, dendpar, dseglength, dnorm, dnx0, dny0
    integer :: iseg,isegtype,istartidx,ipoint,iel

    if (any(shape(Dnx) .ne. shape(Dt)) .or. any(shape(Dny) .ne. shape(Dt))) then
      call output_line ('size(Dt) /= size(Dnx) /= size(Dny)!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'boundary_getNormalVec2D_sim')
      call sys_halt()
    end if

    cpar = BDR_PAR_01
    if (present(cparType)) cpar = cparType

    cnormalMeanCalc = BDR_NORMAL_MEAN
    if (present(cnormalMean)) cnormalMeanCalc = cnormalMean

    if ((iboundCompIdx.gt.rboundary%iboundarycount).or.(iboundCompIdx.lt.0)) then
      call output_line ('iboundCompIdx out of bounds!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'boundary_getNormalVec2D_sim')
      call sys_halt()
    endif

    ! Get the pointers to the segment information arrays for the current
    ! boundary component:
    call storage_getbase_int(rboundary%h_Iintdatavec_handles,p_IintSegInfo_handles)
    call storage_getbase_int(int(p_IintSegInfo_handles(iboundCompIdx)),p_IsegInfo)

    call storage_getbase_int(rboundary%h_Idbldatavec_handles,p_IdbleSegInfo_handles)
    call storage_getbase_double(int(p_IdbleSegInfo_handles(iboundCompIdx)),p_DsegInfo)

    ! Get the segment-count array and the maximum-parameter array
    call storage_getbase_int(rboundary%h_IsegCount,p_IsegCount)
    call storage_getbase_double(rboundary%h_DmaxPar,p_DmaxPar)

    ! Process each parameter value individually
    do iel = 1, size(Dt,2)
      do ipoint = 1, size(Dt,1)

        ! Normalise the parameter value to the range [0,TMAX)
        dpar = Dt(ipoint,iel)
        call boundary_normaliseParValue2D(p_IsegCount,p_DmaxPar,iboundCompIdx,cpar,dpar)

        ! Find segment iseg the parameter value belongs to.
        ! Remember that in the first element in the double precision block of
        ! each segment, the length of the segment is noted!
        call boundary_getSegmentInfo2D(&
            p_IsegCount,p_DmaxPar,p_IsegInfo,p_DsegInfo,iboundCompIdx,dpar,cpar,0,&
            iseg,istartidx,dcurrentpar,dendpar,dseglength,dparloc,isegtype)

        if (dseglength .eq. 0.0_DP) dseglength = 1.0_DP ! trick to avoid div/0

        ! If the local parameter value of dt is <> 0, then the normal vector can
        ! be computed uniquely at boundary point dt.
        ! Otherwise, we have to average the (two) normal vector(s) of the adjacent
        ! boundary components. Check if we are in the simple case.
        if (dparloc .ne. 0.0_DP) then

          call boundary_getNormal2D (rboundary%p_rfparser,iseg,isegType,cpar,&
              dparLoc,dseglength,istartidx,p_IsegInfo,p_DsegInfo,Dnx(ipoint,iel),Dny(ipoint,iel))

        else

          ! Ok, we are at the endpoint of the interval. Now, cnormalMean decides on
          ! what to do.
          Dnx(ipoint,iel) = 0.0_DP
          Dny(ipoint,iel) = 0.0_DP
          if ((cnormalMeanCalc .eq. BDR_NORMAL_MEAN) .or. &
              (cnormalMeanCalc .eq. BDR_NORMAL_RIGHT)) then

            ! Calculate the 'right' normal into dnx/dny.
            !
            ! We already have the segment 'right' to the point.
            ! Get the corresponding normal vector.

            call boundary_getNormal2D (rboundary%p_rfparser,iseg,isegType,cpar,&
                dparLoc,dseglength,istartidx,p_IsegInfo,p_DsegInfo,Dnx(ipoint,iel),Dny(ipoint,iel))

          end if

          if ((cnormalMeanCalc .eq. BDR_NORMAL_MEAN) .or. &
              (cnormalMeanCalc .eq. BDR_NORMAL_LEFT)) then

            ! Calculate the 'left' normal into dnx/dny.
            !
            ! Find segment iseg previous to the parameter value belongs to.
            ! Remember that in the first element in the double precision block of
            ! each segment, the length of the segment is noted!

            call boundary_getSegmentInfo2D(&
                p_IsegCount,p_DmaxPar,p_IsegInfo,p_DsegInfo,iboundCompIdx,dpar,cpar,1,&
                iseg,istartidx,dcurrentpar,dendpar,dseglength,dparloc,isegtype)

            ! Calculate the normal into dnx0/dny0
            call boundary_getNormal2D (rboundary%p_rfparser,iseg,isegType,cpar,&
                dparLoc,dseglength,istartidx,p_IsegInfo,p_DsegInfo,dnx0,dny0)

            ! and add it to dnx/dny
            Dnx(ipoint,iel) = Dnx(ipoint,iel)+dnx0
            Dny(ipoint,iel) = Dny(ipoint,iel)+dny0
          end if

          ! Normalise the vector -- in case two vectors are summed up.
          dnorm = sqrt(Dnx(ipoint,iel)*Dnx(ipoint,iel)+&
                       Dny(ipoint,iel)*Dny(ipoint,iel))
          Dnx(ipoint,iel) = Dnx(ipoint,iel)/dnorm
          Dny(ipoint,iel) = Dny(ipoint,iel)/dnorm
        end if

      end do
    end do

  end subroutine boundary_getNormalVec2D_sim

  !************************************************************************

!<subroutine>

  subroutine boundary_getNormal2D (rfparser,iseg,isegType,cpar,&
      dparLoc,dseglength,istartidx,IsegInfo,DsegInfo,dnx,dny)

!<description>
    ! Calculates the normal of a point on a specific boundary segment.
!</description>

    !<input>

    ! Funciton parser
    type(t_fparser), intent(in) :: rfparser

    ! Segment number (0,1,2,...) of the segment
    integer, intent(in) :: iseg

    ! Segment type
    integer, intent(in) :: isegType

    ! Type of the parameter value (0-1, length par.,...)
    integer, intent(in) :: cpar

    ! Local parameter value of the point in the boundary segment
    real(DP), intent(in) :: dparLoc

    ! Length of the segment
    real(DP), intent(in) :: dseglength

    ! Start index of the segment in DsegInfo
    integer, intent(in) :: istartIdx

    ! Integer segment info array
    integer, dimension(:), intent(in) :: IsegInfo

    ! Double precision segment info array
    real(DP), dimension(:), intent(in) :: DsegInfo

    !</input>

    !<output>

    ! Normal vector
    real(DP), intent(out) :: dnx,dny

    !</output>
!</subroutine>

    ! local variables
    real(DP) :: dploc,dphi,dnorm,dnx0,dny0
    integer :: icomp
    real(DP), dimension(1) :: D1

    dploc = dparloc

    select case (isegType)

      ! case of line
    case (BOUNDARY_TYPE_LINE)

      ! Calculate the x/y components of the normal vector from the
      ! startpoint and the uni direction vector.
      dnx0 =  DsegInfo(istartidx+6)
      dny0 = -DsegInfo(istartidx+5)

      ! case of circle segment
    case (BOUNDARY_TYPE_CIRCLE)

      ! Rescale dparloc with the length of the arc to get a value
      ! between 0 and 1; important for sin/cos functions later.
      ! In the 0-1 parametrisation, this is already the case.
      if (cpar .eq. BDR_PAR_LENGTH) dploc = dploc / dseglength

      ! Get the rotation angle.
      ! Use the initial rotation angle, saved at position 6 of the double
      ! precision data block
      dphi = DsegInfo(istartidx+7) &
           + dploc * (DsegInfo(istartidx+8)-DsegInfo(istartidx+7))

      ! And calculate the x/y components of the normal vector with sin/cos;
      ! the radius is to be found in element 5 of the double precision data block
      dnx0 = DsegInfo(istartidx+5)*cos(dphi)
      dny0 = DsegInfo(istartidx+5)*sin(dphi)

      if (DsegInfo(istartidx+7) .gt. DsegInfo(istartidx+8)) then
        dnx0 = -dnx0
        dny0 = -dny0
      end if

      ! case of analytic expression
    case (BOUNDARY_TYPE_EXPRESSION)
      
      ! Rescale dparloc with the length of the arc to get a value
      ! between 0 and 1; important for sin/cos functions later.
      ! In the 0-1 parametrisation, this is already the case.
      if (cpar .eq. BDR_PAR_LENGTH) dploc = dploc / dseglength
      
      ! Get absolute position of segment
      icomp = IsegInfo(BOUNDARY_SEGHEADER_LENGTH*iseg+BOUNDARY_SEGHEADER_EXPRESSION)

      ! Calculate the x/y coordinates from the function parser.
      D1 = (/dploc/)
      call fparser_evalFunction(rfparser, icomp+3, D1, dny0)
      dny0 = -dny0
      call fparser_evalFunction(rfparser, icomp+4, D1, dnx0)

    case DEFAULT
      call output_line ('Wrong segment type: isegType='//sys_siL(isegType,10), &
                        OU_CLASS_ERROR,OU_MODE_STD,'boundary_getNormal2D')
      call sys_halt()
    end select

    ! Normalise the vector
    dnorm = sqrt(dnx0*dnx0+dny0*dny0)
    dnx = dnx0/dnorm
    dny = dny0/dnorm

  end subroutine boundary_getNormal2D

  !************************************************************************

!<subroutine>

  subroutine boundary_calcNormalVec2D(DsamplePoints, Dpoint, dnx, dny, cpolynomial)

!<description>
  ! This routine calculates the normal vector at the coordinates given
  ! by Dpoints throught the polynomial best-fit curve through the
  ! points DsamplePoints. If the optional parameter cpolynomial is
  ! given, then this values is used for the least-squares
  ! fit. Otherwise, the number of sampling points minus one is used
  ! for the least-squares fit.
!</description>

!<input>

    ! Coordinates of the sampling points
    ! Dimension: DsamplePoints(nim,nsamplePoints)
    real(DP), dimension(:,:), intent(in) :: DsamplePoints

    ! Coordinate of the point where to evaluate the normal vector
    ! Dimension: Dpoint(ndim)
    real(DP), dimension(:), intent(in) :: Dpoint

    ! OPTIONAL: Polynomial order of the least-sqaures fit
    integer, intent(in), optional :: cpolynomial

!</input>

!<output>

    ! x-coordinate of normal vector
    real(DP), intent(out) :: dnx

    ! y-coordinate of normal vector
    real(DP), intent(out) :: dny

!</output>

!</subroutine>

    ! local variables
    real(DP), dimension(:,:), allocatable :: Da
    real(DP), dimension(:), allocatable :: Db,Dc
    real(DP) :: dscale,det1,det2
    integer :: ipolynomial,nsamplepoints

    ipolynomial = size(DsamplePoints,2)-1
    if (present(cpolynomial)) ipolynomial = cpolynomial

    ! Get number of sampling points
    nsamplepoints = size(DsamplePoints,2)

    ! What polynomial degree should be used?
    select case(ipolynomial)

    case (1)
      ! linear interpolation - straight line
      !
      ! The normal to a straight line is constant and does not
      ! depend on the point on the line where it is evaluated

      ! Do we have two sampling points? - would make a least squares
      ! fit unnecessary
      if (nsamplepoints .eq. 2) then

        ! Compute outward normal vector from the fomula
        !   $ n_x = y_2-y-1 $
        !   $ n_y = x_1-x_2 $
        dnx = DsamplePoints(2,2) - DsamplePoints(2,1)
        dny = DsamplePoints(1,1) - DsamplePoints(1,2)

        ! ... and normalise the vector to unity
        dscale = sqrt(dnx**2+dny**2)
        dnx = dnx/dscale
        dny = dny/dscale

      elseif (nsamplepoints .gt. 2) then
        ! Do a linear least squares fit of the sampling points

        ! Allocate temporal memory
        allocate(Da(2,2), Db(2), Dc(2))

        ! Compute 2x2-matrix
        Da(1,1) = nsamplepoints
        Da(1,2) = sum(DsamplePoints(1,:))
        Da(2,1) = Da(1,2)
        Da(2,2) = sum(DsamplePoints(1,:)**2)

        ! Compute 2x1-vector
        Db(1) = sum(DsamplePoints(2,:))
        Db(2) = sum(DsamplePoints(1,:)*DsamplePoints(2,:))

        ! Solve linear 2x2-system via Cramers rule
        det1 = Da(1,1)*Da(2,2)-Da(1,2)*Da(2,1)

        ! Check if system has unique solution
        if (abs(det1) .gt. SYS_EPSREAL_DP) then
          ! Linear system has a unique solution. Thus, the
          ! tangential vector is given by the slope of the line
          det2 = Db(1)*Da(2,2)-Db(2)*Da(2,1)
          Dc(1) = det2/det1
          det2 = Da(1,1)*Db(2)-Da(1,2)*Db(1)
          Dc(2) = det2/det1

          ! Check if the x-value of the first sampling point is
          ! smaller than the x-value of the last sampling point
          if (DsamplePoints(1,1) .lt.&
              DsamplePoints(1,nsamplepoints)) then
            ! Adopt slope of the line as tangential vector
            dscale = sqrt(Dc(2)**2+1)
            dnx =  Dc(2)/dscale
            dny = -1.0_DP/dscale
          else
            ! Adopt negative slope of the line as tangential vector
            dscale = sqrt(Dc(2)**2+1)
            dnx = -Dc(2)/dscale
            dny =  1.0_DP/dscale
          end if
        else
          ! Linear system has no unique solution because the
          ! straight line is parallel to the y-axis. Therefore, we
          ! determine the direction of the line and return the
          ! outward unit normal vector accordingly.
          if (DsamplePoints(2,1) .lt.&
              DsamplePoints(2,nsamplepoints)) then
            dnx = 1.0_DP
            dny = 0.0_DP
          else
            dnx = -1.0_DP
            dny =  0.0_DP
          end if
        end if

        ! Deallocate temporal memory
        deallocate(Da, Db, Dc)

      else
        call output_line ('Insufficient number of sampling points!', &
            OU_CLASS_WARNING,OU_MODE_STD,'boundary_calcNormalVec2D')
        call sys_halt()
      end if

    case default
      call output_line ('Unsupported polynomial order!', &
          OU_CLASS_WARNING,OU_MODE_STD,'boundary_calcNormalVec2D')
      call sys_halt()
    end select

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine boundary_calcNormalVec2D_mult(DsamplePoints, Dpoints, Dnx, Dny, cpolynomial)

!<description>
  ! This routine calculates the normal vector at the coordinates given
  ! by Dpoints throught the polynomial best-fit curve through the
  ! points DsamplePoints. If the optional parameter cpolynomial is
  ! given, then this values is used for the least-squares
  ! fit. Otherwise, the number of sampling points minus one is used
  ! for the least-squares fit.
!</description>

!<input>

    ! Coordinates of the sampling points
    ! Dimension: DsamplePoints(nim,nsamplePoints)
    real(DP), dimension(:,:), intent(in) :: DsamplePoints

    ! Coordinate of the point where to evaluate the normal vector
    ! Dimension: Dpoint(ndim,npoints)
    real(DP), dimension(:,:), intent(in) :: Dpoints

    ! OPTIONAL: Polynomial order of the least-sqaures fit
    integer, intent(in), optional :: cpolynomial

!</input>

!<output>

    ! x-coordinates of normal vector
    ! Dimension: Dnx(npoints)
    real(DP), dimension(:), intent(out) :: Dnx

    ! y-coordinates of normal vector
    ! Dimension: Dny(npoints)
    real(DP), dimension(:), intent(out) :: Dny

!</output>

!</subroutine>

    ! local variables
    real(DP), dimension(:,:), allocatable :: Da
    real(DP), dimension(:), allocatable :: Db,Dc
    real(DP) :: dscale,det1,det2
    integer :: ipolynomial,ipoint,nsamplepoints,npoints

    ipolynomial = size(DsamplePoints,2)-1
    if (present(cpolynomial)) ipolynomial = cpolynomial

    if ((size(Dnx) .ne. size(Dpoints,2)) .or. (size(Dny) .ne. size(Dpoints,2))) then
      call output_line ('size(Dpoints,2) /= size(Dnx) /= size(Dny)!', &
          OU_CLASS_ERROR,OU_MODE_STD,'boundary_calcNormalVec2D_mult')
      call sys_halt()
    end if

    ! Get number of points
    npoints = size(Dpoints,2)
    nsamplepoints = size(DsamplePoints,2)

    ! What polynomial degree should be used?
    select case(ipolynomial)

    case (1)
      ! linear interpolation - straight line
      !
      ! The normal to a straight line is constant and does not
      ! depend on the point on the line where it is evaluated

      ! Do we have two sampling points? - would make a least squares
      ! fit unnecessary
      if (nsamplepoints .eq. 2) then

        do ipoint = 1, npoints
          ! Compute outward normal vector from the fomula
          !   $ n_x = y_2-y-1 $
          !   $ n_y = x_1-x_2 $
          Dnx(ipoint) = DsamplePoints(2,2) - DsamplePoints(2,1)
          Dny(ipoint) = DsamplePoints(1,1) - DsamplePoints(1,2)

          ! ... and normalise the vector to unity
          dscale = sqrt(Dnx(ipoint)**2+Dny(ipoint)**2)
          Dnx(ipoint) = Dnx(ipoint)/dscale
          Dny(ipoint) = Dny(ipoint)/dscale
        end do

      elseif (nsamplepoints .gt. 2) then
        ! Do a linear least squares fit of the sampling points

        ! Allocate temporal memory
        allocate(Da(2,2), Db(2), Dc(2))

        ! Compute 2x2-matrix
        Da(1,1) = nsamplepoints
        Da(1,2) = sum(DsamplePoints(1,:))
        Da(2,1) = Da(1,2)
        Da(2,2) = sum(DsamplePoints(1,:)**2)

        ! Compute 2x1-vector
        Db(1) = sum(DsamplePoints(2,:))
        Db(2) = sum(DsamplePoints(1,:)*DsamplePoints(2,:))

        ! Solve linear 2x2-system via Cramers rule
        det1 = Da(1,1)*Da(2,2)-Da(1,2)*Da(2,1)

        ! Check if system has unique solution
        if (abs(det1) .gt. SYS_EPSREAL_DP) then
          ! Linear system has a unique solution. Thus, the
          ! tangential vector is given by the slope of the line
          det2 = Db(1)*Da(2,2)-Db(2)*Da(2,1)
          Dc(1) = det2/det1
          det2 = Da(1,1)*Db(2)-Da(1,2)*Db(1)
          Dc(2) = det2/det1
          ! Check if the x-value of the first sampling point is
          ! smaller than the x-value of the last sampling point
          if (DsamplePoints(1,1) .lt.&
              DsamplePoints(1,nsamplepoints)) then
            ! Adopt slope of the line as tangential vector
            do ipoint = 1, npoints
              dscale = sqrt(Dc(2)**2+1)
              Dnx(ipoint) =  Dc(2)/dscale
              Dny(ipoint) = -1.0_DP/dscale
            end do
          else
            ! Adopt negative slope of the line as tangential vector
            do ipoint = 1, npoints
              dscale = sqrt(Dc(2)**2+1)
              Dnx(ipoint) = -Dc(2)/dscale
              Dny(ipoint) =  1.0_DP/dscale
            end do
          end if
        else
          ! Linear system has no unique solution because the
          ! straight line is parallel to the y-axis. Therefore, we
          ! determine the direction of the line and return the
          ! outward unit normal vector accordingly.
          if (DsamplePoints(2,1) .lt.&
              DsamplePoints(2,nsamplepoints)) then
            do ipoint = 1, npoints
              Dnx(ipoint) = 1.0_DP
              Dny(ipoint) = 0.0_DP
            end do
          else
            do ipoint = 1, npoints
              Dnx(ipoint) = -1.0_DP
              Dny(ipoint) =  0.0_DP
            end do
          end if
        end if

        ! Deallocate temporal memory
        deallocate(Da, Db, Dc)

      else
        call output_line ('Insufficient number of sampling points!', &
            OU_CLASS_WARNING,OU_MODE_STD,'boundary_calcNormalVec2D_mult')
        call sys_halt()
      end if

    case default
      call output_line ('Unsupported polynomial order!', &
          OU_CLASS_WARNING,OU_MODE_STD,'boundary_calcNormalVec2D_mult')
      call sys_halt()
    end select

  end subroutine

  !************************************************************************

!<subroutine>

  subroutine boundary_calcNormalVec2D_sim(DsamplePoints, Dpoints, Dnx, Dny, cpolynomial)

!<description>
  ! This routine calculates the normal vector at the coordinates given
  ! by Dpoints throught the polynomial best-fit curve through the
  ! points DsamplePoints. If the optional parameter cpolynomial is
  ! given, then this values is used for the least-squares
  ! fit. Otherwise, the number of sampling points minus one is used
  ! for the least-squares fit.
!</description>

!<input>

    ! Coordinates of the sampling points
    ! Dimension: DsamplePoints(nim,nsamplePoints,nelements)
    real(DP), dimension(:,:,:), intent(in) :: DsamplePoints

    ! Coordinate of the point where to evaluate the normal vector
    ! Dimension: Dpoint(ndim,npoints,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Dpoints

    ! OPTIONAL: Polynomial order of the least-sqaures fit
    integer, intent(in), optional :: cpolynomial

!</input>

!<output>

    ! x-coordinates of normal vector
    ! Dimension: Dnx(npoints,nelements)
    real(DP), dimension(:,:), intent(out) :: Dnx

    ! y-coordinates of normal vector
    ! Dimension: Dny(npoints,nelements)
    real(DP), dimension(:,:), intent(out) :: Dny

!</output>

!</subroutine>

    ! local variables
    real(DP), dimension(:,:), allocatable :: Da
    real(DP), dimension(:), allocatable :: Db,Dc
    real(DP) :: dscale,det1,det2
    integer :: ipolynomial,ipoint,iel,nsamplepoints,npoints,nelements

    ipolynomial = size(DsamplePoints,2)-1
    if (present(cpolynomial)) ipolynomial = cpolynomial

    if ((size(Dnx,1) .ne. size(Dpoints,2)) .or.&
        (size(Dnx,2) .ne. size(Dpoints,3)) .or.&
        (size(Dny,1) .ne. size(Dpoints,2)) .or.&
        (size(Dny,2) .ne. size(Dpoints,3))) then
      call output_line ('size(Dpoints,2:3) /= size(Dnx) /= size(Dny)!', &
          OU_CLASS_ERROR,OU_MODE_STD,'boundary_calcNormalVec2D_sim')
      call sys_halt()
    end if

    ! Get number of points and elements
    npoints = size(Dpoints,2)
    nelements = size(Dpoints,3)
    nsamplepoints = size(DsamplePoints,2)

    ! What polynomial degree should be used?
    select case(ipolynomial)

    case (1)
      ! linear interpolation - straight line
      !
      ! The normal to a straight line is constant and does not
      ! depend on the point on the line where it is evaluated

      ! Do we have two sampling points? - would make a least squares
      ! fit unnecessary
      if (nsamplepoints .eq. 2) then

        do ipoint = 1, npoints
          do iel = 1, nelements
            ! Compute outward normal vector from the fomula
            !   $ n_x = y_2-y-1 $
            !   $ n_y = x_1-x_2 $
            Dnx(ipoint,iel) = DsamplePoints(2,2,iel) - DsamplePoints(2,1,iel)
            Dny(ipoint,iel) = DsamplePoints(1,1,iel) - DsamplePoints(1,2,iel)

            ! ... and normalise the vector to unity
            dscale = sqrt(Dnx(ipoint,iel)**2+Dny(ipoint,iel)**2)
            Dnx(ipoint,iel) = Dnx(ipoint,iel)/dscale
            Dny(ipoint,iel) = Dny(ipoint,iel)/dscale
          end do
        end do

      elseif (nsamplepoints .gt. 2) then
        ! Do a linear least squares fit of the sampling points

        ! Allocate temporal memory
        allocate(Da(2,2), Db(2), Dc(2))

        do iel = 1, nelements
          ! Compute 2x2-matrix
          Da(1,1) = nsamplepoints
          Da(1,2) = sum(DsamplePoints(1,:,iel))
          Da(2,1) = Da(1,2)
          Da(2,2) = sum(DsamplePoints(1,:,iel)**2)

          ! Compute 2x1-vector
          Db(1) = sum(DsamplePoints(2,:,iel))
          Db(2) = sum(DsamplePoints(1,:,iel)*DsamplePoints(2,:,iel))

          ! Solve linear 2x2-system via Cramers rule
          det1 = Da(1,1)*Da(2,2)-Da(1,2)*Da(2,1)

          ! Check if system has unique solution
          if (abs(det1) .gt. SYS_EPSREAL_DP) then
            ! Linear system has a unique solution. Thus, the
            ! tangential vector is given by the slope of the line
            det2 = Db(1)*Da(2,2)-Db(2)*Da(2,1)
            Dc(1) = det2/det1
            det2 = Da(1,1)*Db(2)-Da(1,2)*Db(1)
            Dc(2) = det2/det1
            ! Check if the x-value of the first sampling point is
            ! smaller than the x-value of the last sampling point
            if (DsamplePoints(1,1,iel) .lt.&
                DsamplePoints(1,nsamplepoints,iel)) then
              ! Adopt slope of the line as tangential vector
              do ipoint = 1, npoints
                dscale = sqrt(Dc(2)**2+1)
                Dnx(ipoint,iel) =  Dc(2)/dscale
                Dny(ipoint,iel) = -1.0_DP/dscale
              end do
            else
              ! Adopt negative slope of the line as tangential vector
              do ipoint = 1, npoints
                dscale = sqrt(Dc(2)**2+1)
                Dnx(ipoint,iel) = -Dc(2)/dscale
                Dny(ipoint,iel) =  1.0_DP/dscale
              end do
            end if
          else
            ! Linear system has no unique solution because the
            ! straight line is parallel to the y-axis. Therefore, we
            ! determine the direction of the line and return the
            ! outward unit normal vector accordingly.
            if (DsamplePoints(2,1,iel) .lt.&
                DsamplePoints(2,nsamplepoints,iel)) then
              do ipoint = 1, npoints
                Dnx(ipoint,iel) = 1.0_DP
                Dny(ipoint,iel) = 0.0_DP
              end do
            else
              do ipoint = 1, npoints
                Dnx(ipoint,iel) = -1.0_DP
                Dny(ipoint,iel) =  0.0_DP
              end do
            end if
          end if
        end do

        ! Deallocate temporal memory
        deallocate(Da, Db, Dc)

      else
        call output_line ('Insufficient number of sampling points!', &
            OU_CLASS_WARNING,OU_MODE_STD,'boundary_calcNormalVec2D_sim')
        call sys_halt()
      end if

    case default
      call output_line ('Unsupported polynomial order!', &
          OU_CLASS_WARNING,OU_MODE_STD,'boundary_calcNormalVec2D_sim')
      call sys_halt()
    end select

  end subroutine

end module boundary
