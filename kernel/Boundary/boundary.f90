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
!#     -> Get the number of bonudary components 
!#
!# 4.) boundary_igetNsegments
!#     -> Get the number of boundary segments in a boundary component
!#
!# 5.) boundary_dgetMaxParVal
!#     -> Return the maximum parameter value of a boundary component
!#
!# 6.) boundary_getCoords
!#     -> Calculate the coordinates of a point given by its parameter value
!#
!# 7.) boundary_createRegion
!#     -> Get the characteristics of a boundary segment and create
!#        a boundary region structure from it.
!#
!# 8.) boundary_isInRegion
!#     -> Tests whether a node with a specific parameter value
!#        is inside of a given boundary region or not.
!#
!# 9.) boundary_convertParameter
!#     -> Allows to convert a parameter value from 0-1 parametrisation to
!#        length parametrisation and back
!#
!# 10.) boundary_convertParameterList
!#      -> Converts a list of parameter values from 0-1 parametrisation to
!#         length parametrisation and back
!#
!# 11.) boundary_getRegionLength
!#      -> Calculates the length of a boundary region.
!#
!# 12.) boundary_getNormalVec2D
!#      -> Calculate the outward unit normal vector 
!#         of a boundary component in 2D
!#
!# </purpose>
!##############################################################################

module boundary

  use storage
  use fsystem
  use error
  use io
  use genoutput

  implicit none

!<constants>
!<constantblock>
  !maximal degree of NURBS
  integer, parameter :: PAR_MAXNURBSDEGREE=35
!</constantblock>

!<constantblock description="types of boundary segments">
  ! boundary segment type line
  integer, parameter :: BOUNDARY_TYPE_LINE = 1

  ! boundary segment type circle
  integer, parameter :: BOUNDARY_TYPE_CIRCLE = 2

  ! boundary segment type open nurbs
  integer, parameter :: BOUNDARY_TYPE_OPENNURBS = 3

  ! boundary segment type closed nurbs
  integer, parameter :: BOUNDARY_TYPE_CLOSEDNURBS = 4

  ! boundary segment analytic
  integer, parameter :: BOUNDARY_TYPE_ANALYTIC = 5
!</constantblock>

!<constantblock description="kinds of boundary segments">
  ! boundary kind fictitious
  integer, parameter :: BOUNDARY_KIND_FICTITIOUS = 0

  ! boundary kind geometric
  integer, parameter :: BOUNDARY_KIND_GEOMETRIC = 1
!</constantblock>

!<constantblock description="boundary segment header definition">

  ! boundary segment header offset for type
  integer, parameter :: BOUNDARY_SEGHEADER_TYPE = 1

  ! boundary segment header offset for offset in the data vector
  integer, parameter :: BOUNDARY_SEGHEADER_OFFSET = 2

  ! boundary segment header offset for nurbs degree
  integer, parameter :: BOUNDARY_SEGHEADER_NURBSDEGREE = 3

  ! boundary segment header offset for number of control points
  integer, parameter :: BOUNDARY_SEGHEADER_NCNTRLPNTS = 4

  ! boundary segment header length
  integer, parameter :: BOUNDARY_SEGHEADER_LENGTH = 4

!</constantblock>

!<constantblock description="Boundary region type qualifier">

  ! The boundary region is a usual curved region specified
  ! by min/max. parameter value
  integer, parameter :: BDR_TP_CURVE = 0

  ! The boundary region is a point region, i.e. consisting of
  ! only one point. Min/Max. parameter values are identical.
  integer, parameter :: BDR_TP_POINT = 1

!</constantblock>

!<constantblock description="Bitfield constants for t\_boundaryRegion%iproperties">

  ! The startpoint (point with min. parameter value) belongs
  ! to the boundary segment. If not set, the point with the
  ! min. parameter value does not belong to the boundary segment.
  integer(I32), parameter :: BDR_PROP_WITHSTART = 2**0

  ! The endpoint (point with min. parameter value) belongs
  ! to the boundary segment. If not set, the point with the
  ! max. parameter value does not belong to the boundary segment.
  integer(I32), parameter :: BDR_PROP_WITHEND = 2**1

!</constantblock>

!<constantblock description="Type constants for type of parametrisation to use.">

  ! Use 0-1 parametrisation
  integer, parameter :: BDR_PAR_01       = 0

  ! Use parametrisation for the arc length
  integer, parameter :: BDR_PAR_LENGTH   = 1

!</constantblock>

!<constantblock description="Type constants for cnormalMean in the calculation of normal vectors.">

  ! Calculate the mean of the left and right normal.
  integer, parameter :: BDR_NORMAL_MEAN         = 0

  ! Calculate the right normal (i.e. from the point to the interval with
  ! increasing parameter value).
  integer, parameter :: BDR_NORMAL_RIGHT        = 1

  ! Calculate the left normal (i.e. from the point to the interval with
  ! decreasing parameter value).
  integer, parameter :: BDR_NORMAL_LEFT         = 2

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
    integer :: iboundSegIdx = 0

    ! Bitfield specifying properties of the region. A combination
    ! of BDR_PROP_xxxx constants.
    integer(I32) :: iproperties = BDR_PROP_WITHSTART
  
  end type t_boundaryRegion
  
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

  end type t_boundary
  
!</typeblock>
! </types>

  private :: boundary_getSegmentInfo2D

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

  !boundary structure
  type(t_boundary), intent(IN) :: rboundary
  
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
  type(t_boundary), intent(IN) :: rboundary

  ! boundary index
  integer :: iboundCompIdx

  !start parameter value
  real(DP), intent(IN) :: dt1

  !end parameter value
  real(DP), intent(IN) :: dt2
  
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
  type(t_boundary), intent(IN) :: rboundary

  ! index of boundary component
  integer, intent(IN) :: iboundCompIdx
  
  ! OPTIONAL: Type of parametrisation to use.
  ! One of the BDR_PAR_xxxx constants. If not given, BDR_PAR_01 is assumed.
  integer, intent(IN), optional :: cparType
  
!</input>
  
!</function>

    real(DP),dimension(:),pointer :: p_DmaxPar
    integer(I32),dimension(:),pointer :: p_IsegCount

    !if iboundCompIdx exceeds the total number of boundary components or is negative, abort
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

  !boundary structure
  type(t_boundary), intent(IN) :: rboundary

  !index of boundary component
  integer, intent(IN) :: iboundCompIdx
  
!</input>
  
!</function>

    integer(I32),dimension(:),pointer :: p_IsegCount

    !if iboundCompIdx exceeds the total number of boundary components or is negative, abort
    if ((iboundCompIdx .gt. rboundary%iboundarycount) .or. (iboundCompIdx.lt.0)) then
      call output_line ('iboundCompIdx out of bounds!', &
          OU_CLASS_WARNING,OU_MODE_STD,'boundary_dgetNsegments')
      boundary_igetNsegments = -1
      return
    endif
    
    !get vector with component length
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
  character(LEN=*), intent(IN) :: sfilename
! </input>
  
!<output>
  ! Boundary structure, to be filled with data
  type(t_boundary), intent(OUT) :: rboundary
!</output>
  
!</subroutine>

    ! local variables
  
    ! Input channel for reading
    integer :: iunit
    integer :: ibcomponent, isegment, ibct,ihandle
    integer(I32), dimension(:), pointer :: p_IsegInfo, p_IsegCount
    integer(I32), dimension(:), pointer :: p_IdbleSegInfo_handles,p_IintSegInfo_handles
    real(DP), dimension(:), pointer :: p_DsegInfo, p_DmaxPar
    integer :: ityp, nspline, npar
    integer :: isegrel
    integer :: idblemem  ! Counts the memory we need
    real(DP) :: dl
    
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
    call storage_new1D("boundary_read_prm", "h_Idbldatavec_handles", &
                       int(rboundary%iboundarycount,I32), ST_INT, &
                       rboundary%h_Idbldatavec_handles, ST_NEWBLOCK_ZERO)
    call storage_getbase_int(rboundary%h_Idbldatavec_handles, p_IdbleSegInfo_handles)
    
    ! Allocate an array containing of handles. Each handle refers
    ! to integer data for a boundary component.
    call storage_new("boundary_read_prm", "h_Iintdatavec_handles", &
                     int(rboundary%iboundarycount,I32), ST_INT, &
                     rboundary%h_Iintdatavec_handles, ST_NEWBLOCK_ZERO)
    call storage_getbase_int(rboundary%h_Iintdatavec_handles, p_IintSegInfo_handles)

    ! Allocate an array containing the maximum parameter values for each
    ! boundary component in length-parametrisation
    call storage_new("boundary_read_prm", "h_DmaxPar", &
                     int(rboundary%iboundarycount,I32), ST_DOUBLE, &
                     rboundary%h_DmaxPar, ST_NEWBLOCK_ZERO)
    call storage_getbase_double(rboundary%h_DmaxPar, p_DmaxPar)

    ! Allocate an array containing the number of boundary segments in each
    ! boundary component
    call storage_new("boundary_read_prm", "h_IsegCount", &
                     int(rboundary%iboundarycount,I32), ST_INT, &
                     rboundary%h_IsegCount, ST_NEWBLOCK_ZERO)
    call storage_getbase_int(rboundary%h_IsegCount, p_IsegCount)

    ! No we have to gather information about boundary segments.
    ! This makes it necessary to read the boundary definition file
    ! multiple times, but who cares :)
    
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
      
      ! Allocate an integer-array which saves the segtment type and
      ! the start positions (0-based!) of each segment information
      ! block in the segment-array. It's as long as the number
      ! of segments indicates * 2.
      
      call storage_new("boundary_read_prm", "h_Isegcount", &
                       int(2*p_IsegCount(ibct),I32), ST_INT, &
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
        case (1)
          ! Type 1: Line. 
          ! Save the segment type into the first element of each
          ! 2-tuple in the integer array:
          
          p_IsegInfo (1+2*isegment) = BOUNDARY_TYPE_LINE
          
          ! Save the start position of this segment to the segment-start array.
          
          p_IsegInfo (1+2*isegment+1) = idblemem 
          
          ! A line consists of
          ! - Startpoint
          ! - Endpoint, relative to the startpoint
          ! Furthermore, in the first element of the segment we save the
          ! minimum parmeter value and in the second one the length.
          ! So we need 2*2+2 doubles.
          idblemem = idblemem + 6
          
        case (2)
          ! Type 2: Circle / arc. 
          ! Save the segment type into the first element of each
          ! 2-tuple in the integer array:
          
          p_IsegInfo (1+2*isegment) = BOUNDARY_TYPE_CIRCLE
          
          ! Save the start position of this segment to the segment-start array.
          
          p_IsegInfo (1+2*isegment+1) = idblemem 
          
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
        end select
      end do
      
      ! Allocate memory for all the information that defines our
      ! boundary segments:
      
      call storage_new("boundary_read_prm", "h_IdbleSegInfo_handles", &
                       int(idblemem,I32), ST_DOUBLE, &
                       ihandle, ST_NEWBLOCK_ZERO)
      p_IdbleSegInfo_handles(ibct) = ihandle
      
    end do
    
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
        isegrel = p_IsegInfo(1+2*isegment+1)
        
        ! What do we have here? The type identifier is in the
        ! first entry of the 2-tupel in the integer-array:
        select case (p_IsegInfo(1+2*isegment))
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
          
        end select
      end do
      
      ! Save the maximum parameter value for that component
      p_DmaxPar(ibcomponent) = dmaxPar
      
    end do
    
    ! Close the file, finish
    close(iunit)
    
  end subroutine boundary_read_prm

!************************************************************************

!<subroutine>

  subroutine boundary_release(rboundary)

!<description>
  ! This routine releases a boundary object from memory.
!</description>

!<inputoutput>
  ! Boundary structure, to be released.
  type(t_boundary), intent(INOUT) :: rboundary
!</inputoutput>
  
!</subroutine>

    ! local variables
    integer :: i,ihandle
    integer(I32), dimension(:), pointer :: p_IdbleSegInfo_handles,p_IintSegInfo_handles

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
    
    rboundary%iboundarycount_g = 0
    rboundary%iboundarycount = 0
    
    ! That's it...
    
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
  integer(I32), dimension(:), intent(IN) :: IsegCount

  ! Array wirth maximum parameter values for all BC's
  real(DP), dimension(:), intent(IN) :: DmaxPar
  
  ! Number of the boundary component that is under consideration
  integer, intent(IN) :: iboundCompIdx
  
  ! Type of parametrisation, format of dpar (0-1, length par.,...)
  integer(I32), intent(IN) :: cpar
!</input>

!<inputoutput>
  ! Parameter value to be normalised. Is replaced by the normalised
  ! parameter value.
  real(DP), intent(INOUT) :: dpar
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
    integer(I32), dimension(:), intent(IN) :: IsegCount

    ! Double precision array defining the maximum parameter values
    ! of each boundary component
    real(DP), dimension(:), intent(IN) :: DmaxPar

    ! Integer segment info array
    integer(I32), dimension(:), intent(IN) :: IsegInfo
    
    ! Double precision segment info array
    real(DP), dimension(:), intent(IN) :: DsegInfo
    
    ! Number of the boundary component
    integer, intent(IN) :: iboundCompIdx
    
    ! Parameter value of the point of interest
    real(DP), intent(IN) :: dpar
    
    ! Type of parametrisation of dpar (BDR_PAR_01, BDR_PAR_LENGTH,...)
    integer, intent(IN) :: cparType
    
    ! How to orient if the boundary segment is not unique (e.g. on
    ! corner points of the discretisation).
    ! =0: return the segment following the point
    ! =1: return the segment previous to the point
    integer, intent(IN) :: iorientation
!</input>

!<output>
    ! Segment number (0,1,2,...) of the segment that contains dpar
    integer, intent(OUT) :: iseg

    ! Start index of the segment in IsegInfo that contains dpar
    integer, intent(OUT) :: istartidx
    
    ! Start parameter value of the segment that contains dpar
    real(DP), intent(OUT) :: dcurrentpar

    ! End parameter value of the segment that contains dpar
    real(DP), intent(OUT) :: dendpar
    
    ! Segment length of the segment that contains dpar
    real(DP), intent(OUT) :: dseglength
    
    ! Local parameter value of dpar inside of the segment
    real(DP), intent(OUT) :: dparloc

    ! Type of the segment
    integer, intent(OUT) :: isegtype
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
      istartidx = IsegInfo(1+2*iseg+1) 

      dendpar     = iseg + 1.0_DP
      dseglength  = DsegInfo(2+istartidx)
      
    case (BDR_PAR_LENGTH)

      ! In the length-parametrisation, we have to search.
      ! The orientation flag tells us whether to search for "<" or
      ! "<="!
      if (iorientation .ne. 0) then
        do iseg = 0,IsegCount(iboundCompIdx)-1

          ! Determine Start index of the segment in the double-prec. block
          istartidx = IsegInfo(1+2*iseg+1) 

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
          istartidx = IsegInfo(1+2*iseg+1) 

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
    isegtype = IsegInfo(1+2*iseg)

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

  !boundary structure
  type(t_boundary), intent(IN) :: rboundary

  !index of boundary component
  integer, intent(IN) :: iboundCompIdx

  !parametric value of boundary point
  real(DP), intent(IN) :: dt
  
  ! OPTIONAL: Type of parametrisation to use.
  ! One of the BDR_PAR_xxxx constants. If not given, BDR_PAR_01 is assumed.
  integer, intent(IN), optional :: cparType
  
!</input>

!<output>

  !x-coordinate of boundary point
  real(DP), intent(OUT) :: dx

  !y-coordinate of boundary point
  real(DP), intent(OUT) :: dy
    
!</output>

!</subroutine>

    ! local variables
    integer(I32), dimension(:), pointer :: p_IdbleSegInfo_handles,p_IintSegInfo_handles
    integer(I32), dimension(:), pointer :: p_IsegInfo, p_IsegCount
    real(DP), dimension(:), pointer     :: p_DsegInfo, p_DmaxPar
    integer :: cpar ! local copy of cparType
    
    real(DP) :: dpar, dcurrentpar, dparloc, dphi, dendpar, dseglength
    integer :: iseg,isegtype,istartidx
    
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
    
    isegtype = p_IsegInfo(1+2*iseg)
    
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
      
    case DEFAULT
      call output_line ('Wrong segment type: isegType='//sys_siL(isegType,10), &
                        OU_CLASS_ERROR,OU_MODE_STD,'boundary_getCoords')
      call sys_halt()
    end select

  end subroutine boundary_getCoords

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
  type(t_boundary), intent(IN) :: rboundary

  ! index of boundary component
  integer, intent(IN) :: iboundCompIdx

  ! parameter value of boundary point
  real(DP), intent(IN) :: dt
  
  ! Type of parametrisation of DT.
  ! One of the BDR_PAR_xxxx constants. 
  integer, intent(IN) :: cparTypeSource

  ! Type of parametrisation, DT should be converted to.
  ! One of the BDR_PAR_xxxx constants. 
  integer, intent(IN) :: cparTypeDest
  
!</input>

!<result>
  ! The converted parameter value.
  real(DP) dresult
!</result>

!</function>

    ! local variables
    integer(I32), dimension(:), pointer :: p_IdbleSegInfo_handles,p_IintSegInfo_handles
    integer(I32), dimension(:), pointer :: p_IsegInfo, p_IsegCount
    real(DP), dimension(:), pointer     :: p_DsegInfo, p_DmaxPar
    
    real(DP) :: dpar, dcurrentpar, dendpar, dparloc, dseglength, dtmax
    integer :: iseg,istartidx

    ! Small check
    if ((iboundCompIdx.gt.rboundary%iboundarycount).or.(iboundCompIdx.lt.0)) then
      call output_line ('iboundCompIdx out of bounds!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'boundary_convertParameter')
      call sys_halt()
    endif

    ! In case, source and destination type is the same, it's easy.
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
      istartidx = p_IsegInfo(1+2*iseg+1) 

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
        istartidx = p_IsegInfo(1+2*iseg+1) 
        
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

  subroutine boundary_convertParameterList(rboundary, iboundCompIdx, DparSource,&
                                           DparDest, cparTypeSource, cparTypeDest)

!<description>
  ! This function allows to convert an arra of parameter values dt from 0-1 
  ! parametrisation to length parametrisation and back.
  ! cparTypeSource specifies the type of parametrisation of DparSource.
  ! cparTypeDest specifies the destination type, dt should be converted 
  ! to. The return value is the converted array.
!</description>

!<input>

  ! boundary structure
  type(t_boundary), intent(IN) :: rboundary

  ! index of boundary component
  integer, intent(IN) :: iboundCompIdx

  ! An array with parameter values of boundary points that should be converted.
  ! Parameter values < 0 are not converted. Parameter values > maximum parameter
  ! value are mapped to the range 0..maximum parameter value.
  real(DP), dimension(:), intent(IN) :: DparSource

  ! Type of parametrisation of DparSource.
  ! One of the BDR_PAR_xxxx constants. 
  integer, intent(IN) :: cparTypeSource

  ! Type of parametrisation, DT should be converted to.
  ! One of the BDR_PAR_xxxx constants. 
  integer, intent(IN) :: cparTypeDest
  
!</input>

!<inputoutput>
  ! Destinatino array where to write the converted parameter values to.
  ! May coincide with DparSource.
  real(DP), dimension(:), intent(INOUT) :: DparDest
!</inputoutput>

!</subroutine>

    ! local variables
    integer(I32), dimension(:), pointer :: p_IdbleSegInfo_handles,p_IintSegInfo_handles
    integer(I32), dimension(:), pointer :: p_IsegInfo, p_IsegCount
    real(DP), dimension(:), pointer     :: p_DsegInfo, p_DmaxPar
    
    real(DP) :: dpar, dcurrentpar, dendpar, dparloc, dseglength, dtmax
    integer :: iseg,istartidx,iidx

    ! Small check
    if ((iboundCompIdx.gt.rboundary%iboundarycount).or.(iboundCompIdx.lt.0)) then
      call output_line ('iboundCompIdx out of bounds!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'boundary_convertParameterList')
      call sys_halt()
    endif

    ! In case, source and destination type is the same, it's easy.
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

    do iidx = 1,size(DparSource)
    
      ! Don't convert parameter values < 0.
      if (DparSource(iidx) .ge. 0.0_DP) then

        dpar = mod(DparSource(iidx),dtmax)
      
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
          istartidx = p_IsegInfo(1+2*iseg+1) 

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
            istartidx = p_IsegInfo(1+2*iseg+1) 
            
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
          DparDest(iidx) = real(iseg,DP) + dparloc
          
        case (BDR_PAR_LENGTH)
          ! Convert to length parametrisation. dparloc gives us the local
          ! parameter value in 0-1 parametrisation. Interpolate
          ! linearly.
          DparDest(iidx) = dcurrentpar + dparloc*dseglength
        
        end select
        
      end if
     
    end do ! iidx
    
  end subroutine boundary_convertParameterList

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
  type(t_boundary), intent(IN) :: rboundary

  ! Index of boundary component.
  integer, intent(IN) :: iboundCompIdx

  ! Index of the boundary segment.
  ! =0: Create a boundary region that covers the whole boundary component.
  integer, intent(IN) :: iboundSegIdx
  
!</input>

!<output>

  ! Boundary region that is characterised by the boundary segment
  type(t_boundaryRegion), intent(OUT) :: rregion
  
  ! OPTIONAL: Type of parametrisation to use.
  ! One of the BDR_PAR_xxxx constants. If not given, BDR_PAR_01 is assumed.
  integer, intent(IN), optional :: cparType

!</output>
  
!</subroutine>

    ! local variables
    integer(I32), dimension(:), pointer :: p_IdbleSegInfo_handles,p_IintSegInfo_handles
    integer(I32), dimension(:), pointer :: p_IsegInfo, p_IsegCount
    real(DP), dimension(:), pointer     :: p_DsegInfo, p_DmaxPar
    
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
      
      istartidx = p_IsegInfo(1+2*0+1) 
      dcurrentpar = 0.0_DP
      
      ! Determine Start index of the segment in the double-prec. block
      istartidx = p_IsegInfo(1+2*(iboundSegIdx-1)+1) 
      
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
      rregion%isegmentType = p_IsegInfo(1+2*(iboundSegIdx-1))
      
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
  
!<function>

  logical function boundary_isInRegion (rregion,iboundCompIdx,dparam)
  
!<description>
  ! Checks whether a point given by a parameter value of a point is in a 
  ! boundary region.
!</description>

!<input>
  ! A boundary region structure describing a part of a boundary
  type(t_boundaryRegion), intent(IN) :: rregion

  ! The number of the boundary component of the point.
  integer, intent(IN) :: iboundCompIdx
  
  ! The parameter value of the point to be checked.
  ! Must be in the range 0..max. par. value.
  ! Points with negative parameer values are not in the region by definition.
  ! The parametrisation type (0-1 or length parametrisation) must match 
  ! the parametrisation in the boundary region rregion!
  real(DP), intent(IN) :: dparam
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
    ! (Note: Don't change this behaviour to 'rounding' into 0..tmax,
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
      
      ! What's up with the endpoints?
      if ( (dpar1 .eq. dminpar) .and. &
          (iand(rregion%iproperties,BDR_PROP_WITHSTART) .eq. 0)) return
      if ( (dpar1 .eq. dmaxpar) .and. &
          (iand(rregion%iproperties,BDR_PROP_WITHEND) .eq. 0)) return
      
      ! It's inside.
      boundary_isInRegion = .true.
      
    else if ((dpar2 .ge. dminpar) .and. &
             (dpar2 .le. dmaxpar)) then
      
      ! What's up with the endpoints?
      if ( (dpar2 .eq. dminpar) .and. &
           (iand(rregion%iproperties,BDR_PROP_WITHSTART) .eq. 0)) return
      if ( (dpar2 .eq. dmaxpar) .and. &
           (iand(rregion%iproperties,BDR_PROP_WITHEND) .eq. 0)) return
      
      ! It's inside.
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
  type(t_boundary), intent(IN) :: rboundary
  
  ! The boundary reg8ion structure which length is to be computed
  type(t_boundaryRegion), intent(IN) :: rregion
!</input>

!<result>
  ! The length of the boundary region rregion on the boundary rboundary.
!</result>

!</function>

    ! local variables
    real(DP) :: dlen1,dlen2

    ! Is the boundary region parametrised for the length? Then it's easy...
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
    type(t_boundary), intent(IN) :: rboundary

    ! index of boundary component
    integer, intent(IN) :: iboundCompIdx
    
    ! parametric value of boundary point
    real(DP), intent(IN) :: dt
    
    ! OPTIONAL: For points with non-unique normal vectors, this decides
    ! how to calculate the normal vector. One of the BDR_NORMAL_xxxx
    ! constants.
    ! BDR_NORMAL_MEAN calculates the mean of the 'right' and 'left'
    !   normal vector. This is the standard setting
    ! BDR_NORMAL_LEFT calculates the 'left' normal vector (which arises
    !   in the limit when appoximating dt by 0->dt).
    ! BDR_NORMAL_RIGHT calculates the 'right' normal vector (which arises
    !   in the limit when approximating dt by dtmax->dt).
    integer, intent(IN), optional :: cnormalMean
    
    ! OPTIONAL: Type of parametrisation to use.
    ! One of the BDR_PAR_xxxx constants. If not given, BDR_PAR_01 is assumed.
    integer, intent(IN), optional :: cparType
  
!</input>

!<output>

    !x-coordinate of normal vector
    real(DP), intent(OUT) :: dnx
    
    !y-coordinate of normal vector
    real(DP), intent(OUT) :: dny
    
!</output>

!</subroutine>

    ! local variables
    integer(I32), dimension(:), pointer :: p_IdbleSegInfo_handles,p_IintSegInfo_handles
    integer(I32), dimension(:), pointer :: p_IsegInfo, p_IsegCount
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
    
      call getNormal (isegType,cpar,dparLoc,dseglength,istartidx,p_DsegInfo,dnx,dny)

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
        
        call getNormal (isegType,cpar,dparLoc,dseglength,istartidx,p_DsegInfo,dnx,dny)
        
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
        call getNormal (isegType,cpar,dparLoc,dseglength,istartidx,p_DsegInfo,dnx0,dny0)
        
        ! and add it to dnx/dny
        dnx = dnx+dnx0
        dny = dny+dny0
      end if

      ! Normalise the vector -- in case two vectors are summed up.      
      dnorm = sqrt(dnx*dnx+dny*dny)
      dnx = dnx/dnorm
      dny = dny/dnorm
    end if

  contains
  
    subroutine getNormal (isegType,cpar,dparLoc,dseglength,istartidx,DsegInfo,dnx,dny)
      
    ! Calculates the normal of a point on a specific boundary segment.
    
    !<input>
    
    ! Segment type
    integer, intent(IN) :: isegType
    
    ! Type of the parameter value (0-1, length par.,...)
    integer, intent(IN) :: cpar
    
    ! Local parameter value of the point in the boundary segment
    real(DP), intent(IN) :: dparLoc
    
    ! Length of the segment
    real(DP), intent(IN) :: dseglength
    
    ! Start index of the segment in DsegInfo
    integer, intent(IN) :: istartIdx
    
    ! Double precision segment info array
    real(DP), dimension(:), intent(IN) :: DsegInfo
    
    !</input>
    
    !<output>
    
    ! Normal vector
    real(DP), intent(OUT) :: dnx,dny
    
    !</output>
    
      ! local variables
      real(DP) :: dploc,dphi,dnorm
      
      dploc = dparloc
  
      select case (isegType)

      ! case of line
      case (BOUNDARY_TYPE_LINE)

        ! Calculate the x/y components of the normal vector from the 
        ! startpoint and the uni direction vector.
        dnx0 =  p_DsegInfo(istartidx+6)
        dny0 = -p_DsegInfo(istartidx+5)

      ! case of circle segment
      case (BOUNDARY_TYPE_CIRCLE)

        ! Rescale dparloc with the length of the arc to get a value
        ! between 0 and 1; important for sin/cos functions later.
        ! In the 0-1 parametrisation, this is already the case.
        if (cpar .eq. BDR_PAR_LENGTH) dploc = dploc / dseglength

        ! Get the rotation angle.
        ! Use the initial rotation angle, saved at position 6 of the double
        ! precision data block
        dphi = p_DsegInfo(istartidx+7) &
            + dploc * (p_DsegInfo(istartidx+8)-p_DsegInfo(istartidx+7))

        ! And calculate the x/y components of the normal vector with sin/cos;
        ! the radius is to be found in element 5 of the double precision data block
        dnx0 = -p_DsegInfo(istartidx+5)*cos(dphi)
        dny0 = -p_DsegInfo(istartidx+5)*sin(dphi)

      case DEFAULT
        call output_line ('Wrong segment type: isegType='//sys_siL(isegType,10), &
                          OU_CLASS_ERROR,OU_MODE_STD,'getNormal')
        call sys_halt()
      end select

      ! Normalize the vector
      dnorm = sqrt(dnx0*dnx0+dny0*dny0)
      dnx = dnx0/dnorm
      dny = dny0/dnorm
      
    end subroutine getNormal

  end subroutine boundary_getNormalVec2D

end module boundary
