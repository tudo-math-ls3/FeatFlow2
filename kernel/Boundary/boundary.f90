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
!#     -> Calculate the coordinatex of a point given by its parameter value
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
!# 12.) boundary_getNormalVec
!#      -> Calculate the outward unit normal vector of a boundary component
!#
!# </purpose>
!##############################################################################

MODULE boundary

  USE storage
  USE fsystem
  USE error
  USE io
  USE genoutput

  IMPLICIT NONE

!<constants>
!<constantblock>
  !maximal degree of NURBS
  INTEGER, PARAMETER :: PAR_MAXNURBSDEGREE=35
!</constantblock>

!<constantblock description="types of boundary segments">
  ! boundary segment type line
  INTEGER, PARAMETER :: BOUNDARY_TYPE_LINE = 0

  ! boundary segment type circle
  INTEGER, PARAMETER :: BOUNDARY_TYPE_CIRCLE = 1

  ! boundary segment type open nurbs
  INTEGER, PARAMETER :: BOUNDARY_TYPE_OPENNURBS = 2

  ! boundary segment type closed nurbs
  INTEGER, PARAMETER :: BOUNDARY_TYPE_CLOSEDNURBS = 3

  ! boundary segment analytic
  INTEGER, PARAMETER :: BOUNDARY_TYPE_ANALYTIC = 4
!</constantblock>

!<constantblock description="kinds of boundary segments">
  ! boundary kind fictitious
  INTEGER, PARAMETER :: BOUNDARY_KIND_FICTITIOUS = 0

  ! boundary kind geometric
  INTEGER, PARAMETER :: BOUNDARY_KIND_GEOMETRIC = 1
!</constantblock>

!<constantblock description="boundary segment header definition">

  ! boundary segment header offset for type
  INTEGER, PARAMETER :: BOUNDARY_SEGHEADER_TYPE = 1

  ! boundary segment header offset for offset in the data vector
  INTEGER, PARAMETER :: BOUNDARY_SEGHEADER_OFFSET = 2

  ! boundary segment header offset for nurbs degree
  INTEGER, PARAMETER :: BOUNDARY_SEGHEADER_NURBSDEGREE = 3

  ! boundary segment header offset for number of control points
  INTEGER, PARAMETER :: BOUNDARY_SEGHEADER_NCNTRLPNTS = 4

  ! boundary segment header length
  INTEGER, PARAMETER :: BOUNDARY_SEGHEADER_LENGTH = 4

!</constantblock>

!<constantblock description="Boundary region type qualifier">

  ! The boundary region is a usual curved region specified
  ! by min/max. parameter value
  INTEGER, PARAMETER :: BDR_TP_CURVE = 0

  ! The boundary region is a point region, i.e. consisting of
  ! only one point. Min/Max. parameter values are identical.
  INTEGER, PARAMETER :: BDR_TP_POINT = 1

!</constantblock>

!<constantblock description="Bitfield constants for t\_boundaryRegion%iproperties">

  ! The startpoint (point with min. parameter value) belongs
  ! to the boundary segment. If not set, the point with the
  ! min. parameter value does not belong to the boundary segment.
  INTEGER(I32), PARAMETER :: BDR_PROP_WITHSTART = 2**0

  ! The endpoint (point with min. parameter value) belongs
  ! to the boundary segment. If not set, the point with the
  ! max. parameter value does not belong to the boundary segment.
  INTEGER(I32), PARAMETER :: BDR_PROP_WITHEND = 2**1

!</constantblock>

!<constantblock description="Type constants for type of parametrisation to use.">

  ! Use 0-1 parametrisation
  INTEGER, PARAMETER :: BDR_PAR_01       = 0

  ! Use parametrisation for the arc length
  INTEGER, PARAMETER :: BDR_PAR_LENGTH   = 1

!</constantblock>

!<constantblock description="Type constants for cnormalMean in the calculation of normal vectors.">

  ! Calculate the mean of the left and right normal.
  INTEGER, PARAMETER :: BDR_NORMAL_MEAN         = 0

  ! Calculate the right normal (i.e. from the point to the interval with
  ! increasing parameter value).
  INTEGER, PARAMETER :: BDR_NORMAL_RIGHT        = 1

  ! Calculate the left normal (i.e. from the point to the interval with
  ! decreasing parameter value).
  INTEGER, PARAMETER :: BDR_NORMAL_LEFT         = 2

!</constantblock>

!</constants>

!<types>

!<typeblock>

  ! Structure defining a boundary region with minimum/maximum parameter
  ! value in 2D.
  TYPE t_boundaryRegion
  
    ! Type of parametrisation, the parameters refer to.
    ! One of the BDR_PAR_xxxx constants.
    INTEGER  :: cparType  = BDR_PAR_01
  
    ! Minimum parameter value of the region
    REAL(DP) :: dminParam = 0.0_DP
    
    ! Maximum parameter value of the region
    REAL(DP) :: dmaxParam = 0.0_DP
    
    ! Type of the region. One of the BDR_TP_xxxx constants
    INTEGER :: ctype = BDR_TP_CURVE
    
    ! Type of the boundary segment corresponding to this boundary
    ! region. One of the BOUNDARY_TYPE_xxxx constants.
    INTEGER :: isegmentType = BOUNDARY_TYPE_LINE
    
    ! Number of the boundary component that contains the segment
    INTEGER :: iboundCompIdx = 0
    
    ! Maximum parameter value of the boundary component iboundCompIdx
    REAL(DP) :: dmaxParamBC = 0.0_DP
    
    ! Number of the boundary segment that corresponds to this region
    INTEGER :: iboundSegIdx = 0

    ! Bitfield specifying properties of the region. A combination
    ! of BDR_PROP_xxxx constants.
    INTEGER(I32) :: iproperties = BDR_PROP_WITHSTART
  
  END TYPE
  
!</typeblock>

!<typeblock>

  ! Boundary structure of the domain
  TYPE t_boundary
  
    PRIVATE 

    ! number of geometric boundary components
    INTEGER :: iboundarycount_g = -1

    ! total number of boundary components
    INTEGER :: iboundarycount   = -1

    ! handle to double precision array: For every boundary component, maximum
    ! parameter value in length-parametrisation.
    INTEGER :: h_DmaxPar = ST_NOHANDLE

    ! handle for a vector containing the number of segments per boundary component
    INTEGER :: h_IsegCount = ST_NOHANDLE

    ! contains handles for data vectors of boundary components
    INTEGER :: h_Idbldatavec_handles = ST_NOHANDLE

    ! contains handles for offset vectors of boundary components
    INTEGER :: h_Iintdatavec_handles = ST_NOHANDLE

  END TYPE
  
!</typeblock>
! </types>

  PRIVATE :: boundary_getSegmentInfo2D

  CONTAINS

!************************************************************************

!<function>

  INTEGER FUNCTION boundary_igetNBoundComp(rboundary) RESULT (iresult)

!<description>
  ! This function returns the total number of boundary components.
!</description>

!<result>
  ! Total number of boundary components
!</result>

!<input>

  !boundary structure
  TYPE(t_boundary), INTENT(IN) :: rboundary
  
!</input>
  
!</function>

    iresult = rboundary%iboundarycount

  END FUNCTION 

!************************************************************************

!<function>

  RECURSIVE REAL(DP) FUNCTION boundary_dgetLength (rboundary, &
                     iboundCompIdx, dt1, dt2) RESULT (dresult)
    
  !<description>
  ! This function returns the euclidian length between the point
  ! dt1 on the boundary curve and the point dt2.
  !</description>

  !<result>
  ! Real euclidian length
  !</result>

  !<input>

  ! boundary structure
  TYPE(t_boundary), INTENT(IN) :: rboundary

  ! boundary index
  INTEGER :: iboundCompIdx

  !start parameter value
  REAL(DP), INTENT(IN) :: dt1

  !end parameter value
  REAL(DP), INTENT(IN) :: dt2
  
  !</input>
  
!</function>

    REAL(DP) :: dx1,dx2,dy1,dy2,dbl,dtmax
    INTEGER :: ip2,ip1,i

    IF ((iboundCompIdx .gt. rboundary%iboundarycount) .or.&
        (iboundCompIdx.lt.0)) THEN
      CALL output_line ('iboundCompIdx out of bounds!', &
                        OU_CLASS_WARNING,OU_MODE_STD,'boundary_dgetLength')
      dresult = -1
      RETURN
    ENDIF

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

  END FUNCTION 


!************************************************************************

!<function>

  REAL(DP) FUNCTION boundary_dgetMaxParVal(rboundary, iboundCompIdx, cparType)

!<description>
  ! This function returns the parametric length of the boundary component iboundCompIdx
!</description>

!<result>
  ! Parametric length of boundary component iboundCompIdx.
!</result>

!<input>

  ! boundary structure
  TYPE(t_boundary), INTENT(IN) :: rboundary

  ! index of boundary component
  INTEGER, INTENT(IN) :: iboundCompIdx
  
  ! OPTIONAL: Type of parametrisation to use.
  ! One of the BDR_PAR_xxxx constants. If not given, BDR_PAR_01 is assumed.
  INTEGER, INTENT(IN), OPTIONAL :: cparType
  
!</input>
  
!</function>

  REAL(DP),DIMENSION(:),POINTER :: p_DmaxPar
  INTEGER(I32),DIMENSION(:),POINTER :: p_IsegCount

  !if iboundCompIdx exceeds the total number of boundary components or is negative, abort
  if ((iboundCompIdx .gt. rboundary%iboundarycount) .or. (iboundCompIdx.lt.0)) then
    CALL output_line ('iboundCompIdx out of bounds!', &
                      OU_CLASS_WARNING,OU_MODE_STD,'boundary_dgetMaxParVal')
    boundary_dgetMaxParVal = -1
    RETURN
  ENDIF

  IF (PRESENT(cparType)) THEN
    SELECT CASE (cparType)
    CASE (BDR_PAR_LENGTH)
      ! Length-parametrisation
      ! Get vector with component length - depending on the parametrisation
      ! to use.
      CALL storage_getbase_double(rboundary%h_DmaxPar,p_DmaxPar)

      ! Get the maximum parameter value
      boundary_dgetMaxParVal = p_DmaxPar(iboundCompIdx)
      
      RETURN
    END SELECT
  END IF  
  
  ! 0-1 parametrisation. Maximum parameter value is just the number of
  ! segments.
  CALL storage_getbase_int(rboundary%h_IsegCount,p_IsegCount)
  boundary_dgetMaxParVal = REAL(p_IsegCount(iboundCompIdx),DP)

  END FUNCTION 

!************************************************************************

!<function>

  INTEGER FUNCTION boundary_igetNsegments(rboundary, iboundCompIdx)

!<description>
  ! This function returns the number of boundary segments in
  ! the boundary component iboundCompIdx.
!</description>

!<result>
  ! Number of boundary segments on boundary component iboundCompIdx.
!</result>

!<input>

  !boundary structure
  TYPE(t_boundary), INTENT(IN) :: rboundary

  !index of boundary component
  INTEGER, INTENT(IN) :: iboundCompIdx
  
!</input>
  
!</function>

  INTEGER(I32),DIMENSION(:),POINTER :: p_IsegCount

  !if iboundCompIdx exceeds the total number of boundary components or is negative, abort
  if ((iboundCompIdx .gt. rboundary%iboundarycount) .or. (iboundCompIdx.lt.0)) then
    CALL output_line ('iboundCompIdx out of bounds!', &
                      OU_CLASS_WARNING,OU_MODE_STD,'boundary_dgetNsegments')
    boundary_igetNsegments = -1
    RETURN
  ENDIF

  !get vector with component length
  CALL storage_getbase_int(rboundary%h_IsegCount,p_IsegCount)
  
  ! Get the maximum parameter value
  boundary_igetNsegments = p_IsegCount(iboundCompIdx)

  END FUNCTION 

!************************************************************************

!<subroutine>

  SUBROUTINE boundary_read_prm(rboundary, sfilename)

!<description>
  ! This routine reads a .PRM file into memory. The boundary structure
  ! rboundary is initialised with the data from the file.
  ! The parameter sfilename gives the name of the .prm file to read.
  ! If p_rboundary is NULL(), a new structure will be created. 
  ! Otherwise, the existing structure is recreated/updated.
!</description>

!<input>
  ! The name of the .prm file to read.
  CHARACTER(LEN=*), INTENT(IN) :: sfilename
! </input>
  
!<output>
  ! Boundary structure, to be filled with data
  TYPE(t_boundary), INTENT(OUT) :: rboundary
!</output>
  
!</subroutine>

  ! local variables
  
  ! Input channel for reading
  INTEGER :: iunit
  INTEGER :: ibcomponent, isegment, ibct,ihandle
  INTEGER(I32), DIMENSION(:), POINTER :: p_IsegInfo, p_IsegCount
  INTEGER(I32), DIMENSION(:), POINTER :: p_IdbleSegInfo_handles,p_IintSegInfo_handles
  REAL(DP), DIMENSION(:), POINTER :: p_DsegInfo, p_DmaxPar
  INTEGER :: ityp, nspline, npar
  INTEGER :: isegrel
  INTEGER :: idblemem  ! Counts the memory we need
  REAL(DP) :: dl
  
  ! Current parameter value in length-parametrisation
  REAL(DP) :: dmaxpar
  
  ! Open the file
  CALL io_openFileForReading(sfilename, iunit)
  
  ! Read "NBCT"
  READ (iunit,*)
  
  ! Read NBCT - Number of boundary components
  READ (iunit,*) rboundary%iboundarycount_g
  
  rboundary%iboundarycount = rboundary%iboundarycount_g
  
  ! Allocate an array containing handles. Each handle refers
  ! to integer data for a boundary component.
  CALL storage_new1D("boundary_read", "h_Idbldatavec_handles", &
                  INT(rboundary%iboundarycount,I32), ST_INT, &
                  rboundary%h_Idbldatavec_handles, ST_NEWBLOCK_ZERO)
  CALL storage_getbase_int(rboundary%h_Idbldatavec_handles, p_IdbleSegInfo_handles)
  
  ! Allocate an array containing of handles. Each handle refers
  ! to integer data for a boundary component.
  CALL storage_new("boundary_read", "h_Iintdatavec_handles", &
                  INT(rboundary%iboundarycount,I32), ST_INT, &
                  rboundary%h_Iintdatavec_handles, ST_NEWBLOCK_ZERO)
  CALL storage_getbase_int(rboundary%h_Iintdatavec_handles, p_IintSegInfo_handles)

  ! Allocate an array containing the maximum parameter values for each
  ! boundary component in length-parametrisation
  CALL storage_new("boundary_read", "h_DmaxPar", &
                  INT(rboundary%iboundarycount,I32), ST_DOUBLE, &
                  rboundary%h_DmaxPar, ST_NEWBLOCK_ZERO)
  CALL storage_getbase_double(rboundary%h_DmaxPar, p_DmaxPar)

  ! Allocate an array containing the number of boundary segments in each
  ! boundary component
  CALL storage_new("boundary_read", "h_IsegCount", &
                  INT(rboundary%iboundarycount,I32), ST_INT, &
                  rboundary%h_IsegCount, ST_NEWBLOCK_ZERO)
  CALL storage_getbase_int(rboundary%h_IsegCount, p_IsegCount)

  ! No we have to gather information about boundary segments.
  ! This makes it necessary to read the boundary definition file
  ! multiple times, but who cares :)

  ! Initialise the boundary components, loop through them
  
  DO ibcomponent = 1,rboundary%iboundarycount_g
    ! Read "IBCT"
    READ (iunit,*)
    
    ! Read IBCT
    READ (iunit,*) ibct
    
    ! Read "NCOMP"
    READ (iunit,*)
    
    ! Read NCOMP = number of boundary segments in this boundary component
    READ (iunit,*) p_IsegCount(ibct)
    
    ! Allocate an integer-array which saves the segtment type and
    ! the start positions (0-based!) of each segment information
    ! block in the segment-array. It's as long as the number
    ! of segments indicates * 2.
    
    CALL storage_new("boundary_read", "h_Isegcount", &
                    INT(2*p_IsegCount(ibct),I32), ST_INT, &
                    ihandle, ST_NEWBLOCK_ZERO)
    p_IintSegInfo_handles(ibct) = ihandle
    CALL storage_getbase_int(ihandle, p_IsegInfo)

    ! Evaluate the boundary segment definition to get the memory
    ! needed for it:
    READ (iunit,*)         ! "ITYP NSPLINE NPAR"
    
    ! idblemem counts how much memory we need.
    idblemem = 0
    
    ! evaluate the "ITYP NSPLINE NPAR" block
    DO isegment = 0,p_IsegCount(ibct)-1
    
      ! read ITYP NSPLINE NPAR of that segment
      READ (iunit,*) ityp, nspline, npar
      
      ! What do we have here?
      SELECT CASE (ityp)
      CASE (1)
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
      CASE (2)
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
      END SELECT
    END DO
    
    ! Allocate memory for all the information that defines our
    ! boundary segments:
    
    CALL storage_new("boundary_read", "h_IdbleSegInfo_handles", &
                    INT(idblemem,I32), ST_DOUBLE, &
                    ihandle, ST_NEWBLOCK_ZERO)
    p_IdbleSegInfo_handles(ibct) = ihandle
    
  END DO
    
  ! Now we build the segment information.
  
  ! Ignore "PARAMETER"
  READ (iunit,*)  ! "PARAMETER"
  
  ! Read the boundary information
  DO ibcomponent = 1,rboundary%iboundarycount_g
  
    ! dmaxPar counts for every bonudary component the length - and
    ! thus the maximum parameter value.
    dmaxPar = 0.0_DP
  
    ! Get a pointer to the boundary component info.
    ! Here, all double-precision information of the current boundary
    ! component is saved.
    CALL storage_getbase_double(&
         INT(p_IdbleSegInfo_handles(ibcomponent)), p_DsegInfo)
    
    ! Get a pointer to the start position of the segments in the
    ! current boundary component.
    CALL storage_getbase_int(&
         INT(p_IintSegInfo_handles(ibcomponent)), p_ISegInfo)

    ! Build up the segments in the block:
    DO isegment = 0,p_IsegCount(ibcomponent)-1
    
      ! Get the relative position of the segment information in the 
      ! segment information array.
      isegrel = p_IsegInfo(1+2*isegment+1)
      
      ! What do we have here? The type identifier is in the
      ! first entry of the 2-tupel in the integer-array:
      SELECT CASE (p_IsegInfo(1+2*isegment))
      CASE (BOUNDARY_TYPE_LINE)
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
        READ(iunit,*) p_DsegInfo(isegrel+3),p_DsegInfo(isegrel+4)
        
        ! Relative endpoint
        READ(iunit,*) p_DsegInfo(isegrel+5),p_DsegInfo(isegrel+6)
        
        ! Save the initial parameter value (in of the arc length-parametrisation)
        ! to the first entry
        p_DsegInfo(isegrel+1) = dmaxpar
        
        ! This is just necessary for the length-parametrisation. The
        ! 0-1 parametrisation defines the starting point just as the
        ! integer number isegment.
        
        ! Calculate the length and save it in the 2nd position
        dl = SQRT(p_DsegInfo(isegrel+5)**2 + p_DsegInfo(isegrel+6)**2)
        p_DsegInfo(isegrel+2) = dl
        
        ! Normalise the direction vector
        !p_DsegInfo(isegrel+4:isegrel+6) = p_DsegInfo(isegrel+4:isegrel+7)/dl
        
        ! Increase the maximum parameter value
        dmaxPar = dmaxPar + dl
        
      CASE (BOUNDARY_TYPE_CIRCLE)
      
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
        READ(iunit,*) p_DsegInfo(isegrel+3),p_DsegInfo(isegrel+4)
        
        ! Then the radius - the second entry is a dummy; 
        READ(iunit,*) p_DsegInfo(isegrel+5),p_DsegInfo(isegrel+6)
        
        ! Finally get the "arc positions" of the startpoint and the endpoint
        ! of the arc:
        READ(iunit,*) p_DsegInfo(isegrel+7),p_DsegInfo(isegrel+8)
        
        ! Save the initial parameter value (in length-parametrisation)
        ! to the first entry
        p_DsegInfo(isegrel+1) = dmaxpar
        
        ! Now compute the real length of the arc.
        dl = p_DsegInfo(isegrel+5) * &
             ABS(p_DsegInfo(isegrel+8)-p_DsegInfo(isegrel+7))
        p_DsegInfo(isegrel+2) = dl
        
        ! Increase the maximum parameter value
        dmaxPar = dmaxPar + dl
        
      END SELECT
    END DO
    
    ! Save the maximum parameter value for that component
    p_DmaxPar(ibcomponent) = dmaxPar
    
  END DO
  
  ! Close the file, finish
  CLOSE(iunit)
    
  END SUBROUTINE 

!************************************************************************

!<subroutine>

  SUBROUTINE boundary_release(rboundary)

!<description>
  ! This routine releases a boundary object from memory.
!</description>

!<inputoutput>
  ! Boundary structure, to be released.
  TYPE(t_boundary), INTENT(INOUT) :: rboundary
!</inputoutput>
  
!</subroutine>

  ! local variables
  INTEGER :: i,ihandle
  INTEGER(I32), DIMENSION(:), POINTER :: p_IdbleSegInfo_handles,p_IintSegInfo_handles
  
  ! Get the pointers to the segment information arrays for the current
  ! boundary component:
  CALL storage_getbase_int(rboundary%h_Iintdatavec_handles,p_IintSegInfo_handles)
  CALL storage_getbase_int(rboundary%h_Idbldatavec_handles,p_IdbleSegInfo_handles)

  ! Release the handles of the integer- and double-precision
  ! data blocks:
  DO i=1,rboundary%iboundarycount
    ihandle = p_IintSegInfo_handles(i)
    CALL storage_free (ihandle)
    p_IintSegInfo_handles(i) = ihandle
    
    ihandle = p_IdbleSegInfo_handles(i)
    CALL storage_free (ihandle)
    p_IdbleSegInfo_handles(i) = ihandle
  END DO
  
  ! Release all arrays in the structure
  CALL storage_free (rboundary%h_Iintdatavec_handles)
  CALL storage_free (rboundary%h_Idbldatavec_handles)
  CALL storage_free (rboundary%h_IsegCount)
  CALL storage_free (rboundary%h_DmaxPar)
  
  rboundary%iboundarycount_g = 0
  rboundary%iboundarycount = 0

  ! That's it...

  END SUBROUTINE

  !************************************************************************

!<subroutine>

  PURE SUBROUTINE boundary_normaliseParValue2D(IsegCount,DmaxPar,iboundCompIdx,cpar,dpar)

!<description>
  ! INTERNAL SUBROUTINE.
  ! Normalises the parameter value dpar into the range [0,max.par. value).
!</description>

!<input>
  ! Segment-count array
  INTEGER(I32), DIMENSION(:), INTENT(IN) :: IsegCount

  ! Array wirth maximum parameter values for all BC's
  REAL(DP), DIMENSION(:), INTENT(IN) :: DmaxPar
  
  ! Number of the boundary component that is under consideration
  INTEGER, INTENT(IN) :: iboundCompIdx
  
  ! Type of parametrisation, format of dpar (0-1, length par.,...)
  INTEGER(I32), INTENT(IN) :: cpar
!</input>

!<inputoutput>
  ! Parameter value to be normalised. Is replaced by the normalised
  ! parameter value.
  REAL(DP), INTENT(INOUT) :: dpar
!</inputoutput>

!</subroutine>

    SELECT CASE (cpar)
    CASE (BDR_PAR_01)
      dpar = MOD(dpar,REAL(IsegCount(iboundCompIdx),DP))
      !IF (dt .GE. REAL(p_IsegCount(iboundCompIdx),DP) ) THEN
      !  dpar = 0.0_DP
      !ELSE
      !  dpar = dt
      !ENDIF
    CASE (BDR_PAR_LENGTH)
      dpar = MOD(dpar,DmaxPar(iboundCompIdx))
      !IF (dt .GE. p_DmaxPar(iboundCompIdx) ) THEN
      !  dpar = 0.0_DP
      !ELSE
      !  dpar = dt
      !ENDIF
    END SELECT
    
  END SUBROUTINE

  !************************************************************************

!<subroutine>

  PURE SUBROUTINE boundary_getSegmentInfo2D(&
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
    INTEGER(I32), DIMENSION(:), INTENT(IN) :: IsegCount

    ! Double precision array defining the maximum parameter values
    ! of each boundary component
    REAL(DP), DIMENSION(:), INTENT(IN) :: DmaxPar

    ! Integer segment info array
    INTEGER(I32), DIMENSION(:), INTENT(IN) :: IsegInfo
    
    ! Double precision segment info array
    REAL(DP), DIMENSION(:), INTENT(IN) :: DsegInfo
    
    ! Number of the boundary component
    INTEGER, INTENT(IN) :: iboundCompIdx
    
    ! Parameter value of the point of interest
    REAL(DP), INTENT(IN) :: dpar
    
    ! Type of parametrisation of dpar (BDR_PAR_01, BDR_PAR_LENGTH,...)
    INTEGER, INTENT(IN) :: cparType
    
    ! How to orient if the boundary segment is not unique (e.g. on
    ! corner points of the discretisation).
    ! =0: return the segment following the point
    ! =1: return the segment previous to the point
    INTEGER, INTENT(IN) :: iorientation
!</input>

!<output>
    ! Segment number (0,1,2,...) of the segment that contains dpar
    INTEGER, INTENT(OUT) :: iseg

    ! Start index of the segment in IsegInfo that contains dpar
    INTEGER, INTENT(OUT) :: istartidx
    
    ! Start parameter value of the segment that contains dpar
    REAL(DP), INTENT(OUT) :: dcurrentpar

    ! End parameter value of the segment that contains dpar
    REAL(DP), INTENT(OUT) :: dendpar
    
    ! Segment length of the segment that contains dpar
    REAL(DP), INTENT(OUT) :: dseglength
    
    ! Local parameter value of dpar inside of the segment
    REAL(DP), INTENT(OUT) :: dparloc

    ! Type of the segment
    INTEGER, INTENT(OUT) :: isegtype
!</output>

!</subroutine>

    ! Determine the segment
    SELECT CASE (cparType)
    CASE (BDR_PAR_01)

      ! Easy case: 0-1 parametrisation
      !
      ! Standard case: point somewhere in the inner of the segment
      iseg = AINT(dpar)
      dcurrentpar = REAL(iseg,DP)

      ! Subtract the start position of the current boundary component
      ! from the parameter value to get the 'local' parameter value
      ! (0 .le. dparloc .le. length(segment))
      dparloc = dpar - dcurrentpar
      
      IF ((iorientation .NE. 0) .AND. (dparloc .EQ. 0.0_DP)) THEN
        ! We are in a point where the segment is not unique and
        ! the caller wants to have the 'previous' segment. Find it!

        iseg = AINT(dpar)-1
        IF (iseg .EQ. -1) THEN
          ! Get the last segment!
          iseg = IsegCount(iboundCompIdx)-1
          dcurrentpar = REAL(iseg,DP)
          dparloc = dpar + REAL(IsegCount(iboundCompIdx),DP) - dcurrentpar
        ELSE
          dcurrentpar = REAL(iseg,DP)
          dparloc = dpar - dcurrentpar
        END IF
        
      END IF

      ! Determine Start index of the segment in the double-prec. block
      istartidx = IsegInfo(1+2*iseg+1) 

      dendpar     = iseg + 1.0_DP
      dseglength  = DsegInfo(2+istartidx)
      
    CASE (BDR_PAR_LENGTH)

      ! In the length-parametrisation, we have to search.
      ! The orientation flag tells us whether to search for "<" or
      ! "<="!
      IF (iorientation .NE. 0) THEN
        DO iseg = 0,IsegCount(iboundCompIdx)-1

          ! Determine Start index of the segment in the double-prec. block
          istartidx = IsegInfo(1+2*iseg+1) 

          ! Get the start and end parameter value
          dcurrentpar = DsegInfo(1+istartidx)
          dendpar = dcurrentpar + DsegInfo(2+istartidx)
          dseglength = DsegInfo(2+istartidx)

          ! At least one of the IF-commands in the loop will activate
          ! the exit - because of the 'dt' check above!
          IF (dpar .LE. dendpar) EXIT

        END DO
      ELSE
        DO iseg = 0,IsegCount(iboundCompIdx)-1

          ! Determine Start index of the segment in the double-prec. block
          istartidx = IsegInfo(1+2*iseg+1) 

          ! Get the start and end parameter value
          dcurrentpar = DsegInfo(1+istartidx)
          dendpar = dcurrentpar + DsegInfo(2+istartidx)
          dseglength = DsegInfo(2+istartidx)

          ! At least one of the IF-commands in the loop will activate
          ! the exit - because of the 'dt' check above!
          IF (dpar .LT. dendpar) EXIT

        END DO
      END IF

      ! Subtract the start position of the current boundary component
      ! from the parameter value to get the 'local' parameter value
      ! (0 .le. dparloc .le. length(segment))
      dparloc = dpar - dcurrentpar

    CASE DEFAULT
      iseg = 0
      istartidx = 0
      dcurrentpar = 0.0_DP
      dendpar = 0.0_DP
      dseglength = 0.0_DP
      dparloc = 0.0_DP
      isegtype = 0
      RETURN
    END SELECT

    ! Use the segment type to determine how to calculate
    ! the coordinate. Remember that the segment type is noted
    ! in the first element of the integer block of each segment!
    isegtype = IsegInfo(1+2*iseg)

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE boundary_getCoords(rboundary, iboundCompIdx, dt, dx, dy, cparType)

!<description>
  ! This routine returns for a given parameter value dt the
  ! cartesian coordinates of the point on the boundary component 
  ! iboundCompIdx.
!</description>

!<input>

  !boundary structure
  TYPE(t_boundary), INTENT(IN) :: rboundary

  !index of boundary component
  INTEGER, INTENT(IN) :: iboundCompIdx

  !parametric value of boundary point
  REAL(DP), INTENT(IN) :: dt
  
  ! OPTIONAL: Type of parametrisation to use.
  ! One of the BDR_PAR_xxxx constants. If not given, BDR_PAR_01 is assumed.
  INTEGER, INTENT(IN), OPTIONAL :: cparType
  
!</input>

!<output>

  !x-coordinate of boundary point
  REAL(DP), INTENT(OUT) :: dx

  !y-coordinate of boundary point
  REAL(DP), INTENT(OUT) :: dy
    
!</output>

!</subroutine>

  ! local variables
  INTEGER(I32), DIMENSION(:), POINTER :: p_IdbleSegInfo_handles,p_IintSegInfo_handles
  INTEGER(I32), DIMENSION(:), POINTER :: p_IsegInfo, p_IsegCount
  REAL(DP), DIMENSION(:), POINTER     :: p_DsegInfo, p_DmaxPar
  INTEGER :: cpar ! local copy of cparType
  
  REAL(DP) :: dpar, dcurrentpar, dparloc, dphi, dendpar, dseglength
  INTEGER :: iseg,isegtype,istartidx

  cpar = BDR_PAR_01
  IF (PRESENT(cparType)) cpar = cparType

  if ((iboundCompIdx.gt.rboundary%iboundarycount).or.(iboundCompIdx.lt.0)) then
    CALL output_line ('iboundCompIdx out of bounds!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'boundary_getCoords')
    CALL sys_halt()
  ENDIF

  ! Get the pointers to the segment information arrays for the current
  ! boundary component:
  CALL storage_getbase_int(rboundary%h_Iintdatavec_handles,p_IintSegInfo_handles)
  CALL storage_getbase_int(INT(p_IintSegInfo_handles(iboundCompIdx)),p_IsegInfo)
  
  CALL storage_getbase_int(rboundary%h_Idbldatavec_handles,p_IdbleSegInfo_handles)
  CALL storage_getbase_double(INT(p_IdbleSegInfo_handles(iboundCompIdx)),p_DsegInfo)
  
  ! Get the segment-count array and the maximum-parameter array
  CALL storage_getbase_int(rboundary%h_IsegCount,p_IsegCount)
  CALL storage_getbase_double(rboundary%h_DmaxPar,p_DmaxPar)

  ! Normalise the parameter value to the range [0,TMAX)
  dpar = dt
  CALL boundary_normaliseParValue2D(p_IsegCount,p_DmaxPar,iboundCompIdx,cpar,dpar)

  ! Get information about the segment that contains the point
  CALL boundary_getSegmentInfo2D(&
      p_IsegCount,p_DmaxPar,p_IsegInfo,p_DsegInfo,iboundCompIdx,dpar,cpar,0,&
      iseg,istartidx,dcurrentpar,dendpar,dseglength,dparloc,isegtype)
      
  IF (dseglength .EQ. 0.0_DP) dseglength = 1.0_DP ! trick to avoid div/0

  ! Use the segment type to determine how to calculate
  ! the coordinate. Remember that the segment type is noted
  ! in the first element of the integer block of each segment!

  isegtype = p_IsegInfo(1+2*iseg)

  SELECT CASE (isegType)

  ! case of line
  CASE (BOUNDARY_TYPE_LINE)
  
    ! As we save the parametrisation in 0-1 parametrisation,
    ! when we have length-parametrisation, we have to normalise
    ! dparloc to 0 <= dparloc <= 1.
    IF (cpar .EQ. BDR_PAR_LENGTH) dparloc = dparloc / dseglength
  
    ! Calculate the x/y coordinates from the startpoint and
    ! the unit direction vector.
    dx = p_DsegInfo(istartidx+3) + dparloc*p_DsegInfo(istartidx+5)
    dy = p_DsegInfo(istartidx+4) + dparloc*p_DsegInfo(istartidx+6)

  ! case of circle segment
  CASE (BOUNDARY_TYPE_CIRCLE)

    ! Rescale dparloc with the length of the arc to get a value
    ! between 0 and 1; important for sin/cos functions later.
    ! In the 0-1 parametrisation, this is already the case.
    IF (cpar .EQ. BDR_PAR_LENGTH) dparloc = dparloc / dseglength

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

  CASE DEFAULT
    CALL output_line ('Wrong segment type: isegType='//sys_siL(isegType,10), &
                      OU_CLASS_ERROR,OU_MODE_STD,'boundary_getCoords')
    CALL sys_halt()
  END SELECT

  END SUBROUTINE 

!************************************************************************

!<function>

  REAL(DP) FUNCTION boundary_convertParameter(rboundary, iboundCompIdx, dt, &
                                              cparTypeSource, cparTypeDest) &
           RESULT(dresult)

!<description>
  ! This function allows to convert a parameter value dt from 0-1 
  ! parametrisation to length parametrisation and back.
  ! cparTypeSource specifies the type of parametrisation of dt.
  ! cparTypeDest specifies the destination type, dt should be converted 
  ! to. The return value is the converted parameter value.
!</description>

!<input>

  ! boundary structure
  TYPE(t_boundary), INTENT(IN) :: rboundary

  ! index of boundary component
  INTEGER, INTENT(IN) :: iboundCompIdx

  ! parameter value of boundary point
  REAL(DP), INTENT(IN) :: dt
  
  ! Type of parametrisation of DT.
  ! One of the BDR_PAR_xxxx constants. 
  INTEGER, INTENT(IN) :: cparTypeSource

  ! Type of parametrisation, DT should be converted to.
  ! One of the BDR_PAR_xxxx constants. 
  INTEGER, INTENT(IN) :: cparTypeDest
  
!</input>

!<result>
  ! The converted parameter value.
!</result>

!</function>

    ! local variables
    INTEGER(I32), DIMENSION(:), POINTER :: p_IdbleSegInfo_handles,p_IintSegInfo_handles
    INTEGER(I32), DIMENSION(:), POINTER :: p_IsegInfo, p_IsegCount
    REAL(DP), DIMENSION(:), POINTER     :: p_DsegInfo, p_DmaxPar
    
    REAL(DP) :: dpar, dcurrentpar, dendpar, dparloc, dseglength, dtmax
    INTEGER :: iseg,istartidx

    ! Small check
    IF ((iboundCompIdx.gt.rboundary%iboundarycount).or.(iboundCompIdx.lt.0)) then
      CALL output_line ('iboundCompIdx out of bounds!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'boundary_convertParameter')
      CALL sys_halt()
    ENDIF

    ! In case, source and destination type is the same, it's easy.
    IF (cparTypeSource .EQ. cparTypeDest) THEN
      dresult = dt
      RETURN
    END IF

    ! Get the pointers to the segment information arrays for the current
    ! boundary component:
    CALL storage_getbase_int(rboundary%h_Iintdatavec_handles,p_IintSegInfo_handles)
    CALL storage_getbase_int(INT(p_IintSegInfo_handles(iboundCompIdx)),p_IsegInfo)
    
    CALL storage_getbase_int(rboundary%h_Idbldatavec_handles,p_IdbleSegInfo_handles)
    CALL storage_getbase_double(INT(p_IdbleSegInfo_handles(iboundCompIdx)),p_DsegInfo)
    
    ! Get the segment-count array and the maximum-parameter array
    CALL storage_getbase_int(rboundary%h_IsegCount,p_IsegCount)
    CALL storage_getbase_double(rboundary%h_DmaxPar,p_DmaxPar)

    ! If the parameter value exceeds the parameter interval on the boundary 
    ! component, truncate the parameter value!
    SELECT CASE (cparTypeSource)
    CASE (BDR_PAR_01)
      dtmax = REAL(p_IsegCount(iboundCompIdx),DP)
    CASE (BDR_PAR_LENGTH)
      dtmax = p_DmaxPar(iboundCompIdx)
    END SELECT
    dpar = MOD(dt,dtmax)

    ! Find segment iseg the parameter value belongs to.
    ! Remember that in the first element in the double precision block of
    ! each segment, the length of the segment is noted!
    
    dcurrentpar = 0.0_DP
    
    ! Determine the segment 
    SELECT CASE (cparTypeSource)
    CASE (BDR_PAR_01)
    
      ! Easy case: 0-1 parametrisation
      iseg = AINT(dpar)

      ! Determine Start index of the segment in the double-prec. block
      istartidx = p_IsegInfo(1+2*iseg+1) 

      ! Get the segment length for later use
      dseglength  = p_DsegInfo(2+istartidx)
      
      ! Start parameter value in length parametrisation
      dcurrentpar = p_DsegInfo(1+istartidx)
      
      ! Subtract the start position of the current boundary component
      ! from the parameter value to get the 'local' parameter value
      ! (0 .le. dparloc .le. 1.0)
      dparloc = dpar - REAL(iseg,DP)

    CASE (BDR_PAR_LENGTH)
    
      ! In the length-parametrisation, we have to search.
      DO iseg = 0,p_IsegCount(iboundCompIdx)-1
        
        ! Determine Start index of the segment in the double-prec. block
        istartidx = p_IsegInfo(1+2*iseg+1) 
        
        ! Get the start and end parameter value
        dcurrentpar = p_DsegInfo(1+istartidx)
        dseglength = p_DsegInfo(2+istartidx)
        dendpar = dcurrentpar + dseglength
        
        ! At least one of the IF-commands in the loop will activate
        ! the exit - because of the 'dt' check above!
        IF (dpar .LT. dendpar) EXIT
        
      END DO

      ! Subtract the start position of the current boundary component
      ! from the parameter value to get the 'local' parameter value
      ! (0 .le. dparloc .le. length(segment))
      dparloc = dpar - dcurrentpar

      IF (dseglength .EQ. 0.0_DP) dseglength = 1.0_DP ! trick to avoid div/0
      
      ! Divide by the segment length to get the local parameter value
      ! in 0-1 parametrisation.
      dparloc = dparloc / dseglength

    END SELECT
    
    ! The local parameter value dparloc is now always in the range 0..1.
    !    
    ! How shoule we convert?
    SELECT CASE (cparTypeDest)
    CASE (BDR_PAR_01)
      ! Convert to 0-1. 
      ! Take the number of teh segment as basis and add
      ! the local parameter value to get the 0-1 parameter value.
      dresult = REAL(iseg,DP) + dparloc
      
    CASE (BDR_PAR_LENGTH)
      ! Convert to length parametrisation. dparloc gives us the local
      ! parameter value in 0-1 parametrisation. Interpolate
      ! linearly.
      dresult = dcurrentpar + dparloc*dseglength
    
    END SELECT
    
  END FUNCTION

!************************************************************************

!<subroutine>

  SUBROUTINE boundary_convertParameterList(rboundary, iboundCompIdx, &
      DparSource, DparDest, cparTypeSource, cparTypeDest)

!<description>
  ! This function allows to convert an arra of parameter values dt from 0-1 
  ! parametrisation to length parametrisation and back.
  ! cparTypeSource specifies the type of parametrisation of DparSource.
  ! cparTypeDest specifies the destination type, dt should be converted 
  ! to. The return value is the converted array.
!</description>

!<input>

  ! boundary structure
  TYPE(t_boundary), INTENT(IN) :: rboundary

  ! index of boundary component
  INTEGER, INTENT(IN) :: iboundCompIdx

  ! An array with parameter values of boundary points that should be converted.
  REAL(DP), DIMENSION(:), INTENT(IN) :: DparSource

  ! Type of parametrisation of DparSource.
  ! One of the BDR_PAR_xxxx constants. 
  INTEGER, INTENT(IN) :: cparTypeSource

  ! Type of parametrisation, DT should be converted to.
  ! One of the BDR_PAR_xxxx constants. 
  INTEGER, INTENT(IN) :: cparTypeDest
  
!</input>

!<inputoutput>
  ! Destinatino array where to write the converted parameter values to.
  ! May coincide with DparSource.
  REAL(DP), DIMENSION(:), INTENT(INOUT) :: DparDest
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER(I32), DIMENSION(:), POINTER :: p_IdbleSegInfo_handles,p_IintSegInfo_handles
    INTEGER(I32), DIMENSION(:), POINTER :: p_IsegInfo, p_IsegCount
    REAL(DP), DIMENSION(:), POINTER     :: p_DsegInfo, p_DmaxPar
    
    REAL(DP) :: dpar, dcurrentpar, dendpar, dparloc, dseglength, dtmax
    INTEGER :: iseg,istartidx,iidx

    ! Small check
    IF ((iboundCompIdx.gt.rboundary%iboundarycount).or.(iboundCompIdx.lt.0)) then
      CALL output_line ('iboundCompIdx out of bounds!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'boundary_convertParameterList')
      CALL sys_halt()
    ENDIF

    ! In case, source and destination type is the same, it's easy.
    IF (cparTypeSource .EQ. cparTypeDest) THEN
      DparDest(:) = DparSource(:)
      RETURN
    END IF

    ! Get the pointers to the segment information arrays for the current
    ! boundary component:
    CALL storage_getbase_int(rboundary%h_Iintdatavec_handles,p_IintSegInfo_handles)
    CALL storage_getbase_int(INT(p_IintSegInfo_handles(iboundCompIdx)),p_IsegInfo)
    
    CALL storage_getbase_int(rboundary%h_Idbldatavec_handles,p_IdbleSegInfo_handles)
    CALL storage_getbase_double(INT(p_IdbleSegInfo_handles(iboundCompIdx)),p_DsegInfo)
    
    ! Get the segment-count array and the maximum-parameter array
    CALL storage_getbase_int(rboundary%h_IsegCount,p_IsegCount)
    CALL storage_getbase_double(rboundary%h_DmaxPar,p_DmaxPar)

    ! If the parameter value exceeds the parameter interval on the boundary 
    ! component, truncate the parameter value!
    SELECT CASE (cparTypeSource)
    CASE (BDR_PAR_01)
      dtmax = REAL(p_IsegCount(iboundCompIdx),DP)
    CASE (BDR_PAR_LENGTH)
      dtmax = p_DmaxPar(iboundCompIdx)
    END SELECT

    DO iidx = 1,SIZE(DparSource)
    
      dpar = MOD(DparSource(iidx),dtmax)

      ! Find segment iseg the parameter values belong to.
      ! Remember that in the first element in the double precision block of
      ! each segment, the length of the segment is noted!
      
      dcurrentpar = 0.0_DP
      
      ! Determine the segment 
      SELECT CASE (cparTypeSource)
      CASE (BDR_PAR_01)
      
        ! Easy case: 0-1 parametrisation
        iseg = AINT(dpar)

        ! Determine Start index of the segment in the double-prec. block
        istartidx = p_IsegInfo(1+2*iseg+1) 

        ! Get the segment length for later use
        dseglength  = p_DsegInfo(2+istartidx)
        
        ! Start parameter value in length parametrisation
        dcurrentpar = p_DsegInfo(1+istartidx)
        
        ! Subtract the start position of the current boundary component
        ! from the parameter value to get the 'local' parameter value
        ! (0 .le. dparloc .le. 1.0)
        dparloc = dpar - REAL(iseg,DP)

      CASE (BDR_PAR_LENGTH)
      
        ! In the length-parametrisation, we have to search.
        DO iseg = 0,p_IsegCount(iboundCompIdx)-1
          
          ! Determine Start index of the segment in the double-prec. block
          istartidx = p_IsegInfo(1+2*iseg+1) 
          
          ! Get the start and end parameter value
          dcurrentpar = p_DsegInfo(1+istartidx)
          dseglength = p_DsegInfo(2+istartidx)
          dendpar = dcurrentpar + dseglength
          
          ! At least one of the IF-commands in the loop will activate
          ! the exit - because of the 'dt' check above!
          IF (dpar .LT. dendpar) EXIT
          
        END DO

        ! Subtract the start position of the current boundary component
        ! from the parameter value to get the 'local' parameter value
        ! (0 .le. dparloc .le. length(segment))
        dparloc = dpar - dcurrentpar

        IF (dseglength .EQ. 0.0_DP) dseglength = 1.0_DP ! trick to avoid div/0
        
        ! Divide by the segment length to get the local parameter value
        ! in 0-1 parametrisation.
        dparloc = dparloc / dseglength

      END SELECT
      
      ! The local parameter value dparloc is now always in the range 0..1.
      !    
      ! How shoule we convert?
      SELECT CASE (cparTypeDest)
      CASE (BDR_PAR_01)
        ! Convert to 0-1. 
        ! Take the number of teh segment as basis and add
        ! the local parameter value to get the 0-1 parameter value.
        DparDest(iidx) = REAL(iseg,DP) + dparloc
        
      CASE (BDR_PAR_LENGTH)
        ! Convert to length parametrisation. dparloc gives us the local
        ! parameter value in 0-1 parametrisation. Interpolate
        ! linearly.
        DparDest(iidx) = dcurrentpar + dparloc*dseglength
      
      END SELECT
     
    END DO ! iidx
    
  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE boundary_createRegion (rboundary, iboundCompIdx, iboundSegIdx, &
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
  TYPE(t_boundary), INTENT(IN) :: rboundary

  ! Index of boundary component.
  INTEGER, INTENT(IN) :: iboundCompIdx

  ! Index of the boundary segment.
  ! =0: Create a boundary region that covers the whole boundary component.
  INTEGER, INTENT(IN) :: iboundSegIdx
  
!</input>

!<output>

  ! Boundary region that is characterised by the boundary segment
  TYPE(t_boundaryRegion), INTENT(OUT) :: rregion
  
  ! OPTIONAL: Type of parametrisation to use.
  ! One of the BDR_PAR_xxxx constants. If not given, BDR_PAR_01 is assumed.
  INTEGER, INTENT(IN), OPTIONAL :: cparType

!</output>
  
!</subroutine>

  ! local variables
  INTEGER(I32), DIMENSION(:), POINTER :: p_IdbleSegInfo_handles,p_IintSegInfo_handles
  INTEGER(I32), DIMENSION(:), POINTER :: p_IsegInfo, p_IsegCount
  REAL(DP), DIMENSION(:), POINTER     :: p_DsegInfo, p_DmaxPar
  
  REAL(DP) :: dcurrentpar, dendpar, dmaxpar
  INTEGER :: istartidx

  INTEGER :: cpar ! local copy of cparType
  
  cpar = BDR_PAR_01
  IF (PRESENT(cparType)) cpar = cparType

  IF ((iboundCompIdx .GT. rboundary%iboundarycount) .OR. (iboundCompIdx.lt.0)) THEN
    CALL output_line ('iboundCompIdx out of bounds!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'boundary_createregion')
    CALL sys_halt()
  ENDIF

  ! Get the pointers to the segment information arrays for the current
  ! boundary component:
  CALL storage_getbase_int(rboundary%h_Iintdatavec_handles,p_IintSegInfo_handles)
  CALL storage_getbase_int(INT(p_IintSegInfo_handles(iboundCompIdx)),p_IsegInfo)
  
  CALL storage_getbase_int(rboundary%h_Idbldatavec_handles,p_IdbleSegInfo_handles)
  CALL storage_getbase_double(INT(p_IdbleSegInfo_handles(iboundCompIdx)),p_DsegInfo)
  
  ! Get the segment-count array and the maximum-parameter array
  CALL storage_getbase_int(rboundary%h_IsegCount,p_IsegCount)
  CALL storage_getbase_double(rboundary%h_DmaxPar,p_DmaxPar)

  IF (iboundSegIdx .NE. 0) THEN

    IF ((iboundSegIdx .GT. p_IsegCount(iboundCompIdx)) .OR. (iboundSegIdx.lt.0)) THEN
      CALL output_line ('iboundSegIdx out of bounds!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'boundary_createregion')
      CALL sys_halt()
    ENDIF

    ! Find segment iseg the parameter value belongs to.
    ! Remember that in the first element in the double precision block of
    ! each segment, the length of the segment is noted!
    
    istartidx = p_IsegInfo(1+2*0+1) 
    dcurrentpar = 0.0_DP
    
    ! Determine Start index of the segment in the double-prec. block
    istartidx = p_IsegInfo(1+2*(iboundSegIdx-1)+1) 
    
    ! Get the start and end parameter value - depending on the parametrisation
    SELECT CASE (cpar)
    CASE (BDR_PAR_01)
      dcurrentpar = REAL(iboundSegIdx-1,DP)
      dendpar     = REAL(iboundSegIdx-1+1,DP)
      dmaxpar     = REAL(p_IsegCount(iboundCompIdx),DP)
    CASE (BDR_PAR_LENGTH)
      dcurrentpar = p_DsegInfo(1+istartidx)
      dendpar     = dcurrentpar + p_DsegInfo(2+istartidx)
      dmaxpar     = p_DmaxPar(iboundCompIdx)
    END SELECT
  
    ! Set the segment type in the boundary region structure
    rregion%isegmentType = p_IsegInfo(1+2*(iboundSegIdx-1))

  ELSE
  
    ! Create a boundary region that covers the whole boundary component.
    SELECT CASE (cpar)
    CASE (BDR_PAR_01)
      dcurrentpar = REAL(iboundSegIdx-1,DP)
      dmaxpar     = REAL(p_IsegCount(iboundCompIdx),DP)
    CASE (BDR_PAR_LENGTH)
      dcurrentpar = p_DsegInfo(1+istartidx)
      dmaxpar     = p_DmaxPar(iboundCompIdx)
    END SELECT
    
    dcurrentpar = 0.0_DP
    dendpar     = dmaxpar
  
    ! We have an unspecified boundary segment
    rregion%isegmentType = BOUNDARY_TYPE_ANALYTIC
  
  END IF

  ! Create the boundary region structure
  rregion%cparType = cpar
  rregion%dminParam = dcurrentpar
  rregion%dmaxParam = dendpar
  rregion%ctype = BDR_TP_CURVE
  rregion%iproperties = BDR_PROP_WITHSTART
  rregion%iboundCompIdx = iboundCompIdx
  rregion%iboundSegIdx = iboundSegIdx
  rregion%dmaxParamBC = dmaxpar

  END SUBROUTINE
                                    
  ! ***************************************************************************
  
!<function>

  LOGICAL FUNCTION boundary_isInRegion (rregion,iboundCompIdx,dparam)
  
!<description>
  ! Checks whether a point given by a parameter value of a point is in a 
  ! boundary region.
!</description>

!<input>
  ! A boundary region structure describing a part of a boundary
  TYPE(t_boundaryRegion), INTENT(IN) :: rregion

  ! The number of the boundary component of the point.
  INTEGER, INTENT(IN) :: iboundCompIdx
  
  ! The parameter value of the point to be checked.
  ! Must be in the range 0..max. par. value.
  ! The parametrisation type (0-1 or length parametrisation) must match 
  ! the parametrisation in the boundary region rregion!
  REAL(DP), INTENT(IN) :: dparam
!</input>

!<result>
  ! TRUE, if the point is inside the region,
  ! FALSE otherwise.
!</result>

!</function>

  ! local variables
  REAL(DP) :: dminpar,dmaxpar, dpar1,dpar2
  
  ! Default setting: not inside.
  boundary_isInRegion = .FALSE.
  
  ! Correct boundary component?
  IF (iboundCompIdx .NE. rregion%iboundCompIdx) RETURN
  
  ! Does the region cover the complete boundary?
  IF (rregion%dmaxParam-rregion%dminParam .GT. rregion%dmaxParamBC) THEN
    boundary_isInRegion = .TRUE.
    RETURN
  END IF
  
  ! Region does not cover the complete boundary, but a part of it.
  ! Which part?
  
  ! Get the real bounds...
  dminpar = MOD(rregion%dminParam,rregion%dmaxParamBC)
  dmaxpar = MOD(rregion%dmaxParam,rregion%dmaxParamBC)
  
  ! And make sure, dmaxpar >= dminpar!
  IF ((dmaxpar .LE. dminpar) .AND. &
      (rregion%dminParam .NE. rregion%dmaxParam)) THEN
    dmaxpar = dmaxpar + rregion%dmaxParamBC
  END IF
  
  ! Get the parameter value on the boundary - in the set [0..TMAX)
  ! and [TMAX..2*TMAX).
  dpar1 = MOD(dparam,rregion%dmaxParamBC)
  dpar2 = dpar1 + rregion%dmaxParamBC
  
  ! Check if dpar1 or dpar2 is in that region.
  ! For 'normal' boundary regions, dpar1 should be inside, dpar2 not.
  ! For boundary regions crossing the maximum parameter value,
  ! dpar2 will be inside and dpar1 not.
  
  IF ((dpar1 .GE. dminpar) .AND. &
      (dpar1 .LE. dmaxpar)) THEN

    ! What's up with the endpoints?
    IF ( (dpar1 .EQ. dminpar) .AND. &
        (IAND(rregion%iproperties,BDR_PROP_WITHSTART) .EQ. 0)) RETURN
    IF ( (dpar1 .EQ. dmaxpar) .AND. &
        (IAND(rregion%iproperties,BDR_PROP_WITHEND) .EQ. 0)) RETURN
        
    ! It's inside.
    boundary_isInRegion = .TRUE.
    
  ELSE IF ((dpar2 .GE. dminpar) .AND. &
           (dpar2 .LE. dmaxpar)) THEN

    ! What's up with the endpoints?
    IF ( (dpar2 .EQ. dminpar) .AND. &
        (IAND(rregion%iproperties,BDR_PROP_WITHSTART) .EQ. 0)) RETURN
    IF ( (dpar2 .EQ. dmaxpar) .AND. &
        (IAND(rregion%iproperties,BDR_PROP_WITHEND) .EQ. 0)) RETURN
        
    ! It's inside.
    boundary_isInRegion = .TRUE.
    
  END IF

  END FUNCTION
 
  ! ***************************************************************************
  
!<function>

  REAL(DP) FUNCTION boundary_getRegionLength (rboundary,rregion) RESULT (dlength)
  
!<description>
  ! Calculates the length of a part of the boundary identified by rregion
  ! on boundary rboundary.
!</description>

!<input>
  ! Boundary structure, the boundary region should refer to
  TYPE(t_boundary), INTENT(IN) :: rboundary
  
  ! The boundary reg8ion structure which length is to be computed
  TYPE(t_boundaryRegion), INTENT(IN) :: rregion
!</input>

!<result>
  ! The length of the boundary region rregion on the boundary rboundary.
!</result>

!</function>

    ! local variables
    REAL(DP) :: dlen1,dlen2

    ! Is the boundary region parametrised for the length? Then it's easy...
    IF (rregion%cparType .EQ. BDR_PAR_LENGTH) THEN
      dlength = rregion%dmaxParam - rregion%dminParam
      RETURN
    END IF
    
    ! Otherwise, compute the parameter values in length-parametrisation
    ! and subtract them to get the length.
    dlen1 = boundary_convertParameter(rboundary, rregion%iboundCompIdx, &
                                      rregion%dminParam, &
                                      rregion%cparType, BDR_PAR_LENGTH)

    dlen2 = boundary_convertParameter(rboundary, rregion%iboundCompIdx, &
                                      rregion%dmaxParam, &
                                      rregion%cparType, BDR_PAR_LENGTH)
                                      
    dlength = dlen2-dlen1
    
  END FUNCTION

  !************************************************************************

!<subroutine>

  SUBROUTINE boundary_getNormalVec(rboundary, iboundCompIdx, dt, &
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
    TYPE(t_boundary), INTENT(IN) :: rboundary

    ! index of boundary component
    INTEGER, INTENT(IN) :: iboundCompIdx
    
    ! parametric value of boundary point
    REAL(DP), INTENT(IN) :: dt
    
    ! OPTIONAL: For points with non-unique normal vectors, this decides
    ! how to calculate the normal vector. One of the BDR_NORMAL_xxxx
    ! constants.
    ! BDR_NORMAL_MEAN calculates the mean of the 'right' and 'left'
    !   normal vector. This is the standard setting
    ! BDR_NORMAL_LEFT calculates the 'left' normal vector (which arises
    !   in the limit when appoximating dt by 0->dt).
    ! BDR_NORMAL_RIGHT calculates the 'right' normal vector (which arises
    !   in the limit when approximating dt by dtmax->dt).
    INTEGER, INTENT(IN), OPTIONAL :: cnormalMean
    
    ! OPTIONAL: Type of parametrisation to use.
    ! One of the BDR_PAR_xxxx constants. If not given, BDR_PAR_01 is assumed.
    INTEGER, INTENT(IN), OPTIONAL :: cparType
  
!</input>

!<output>

    !x-coordinate of normal vector
    REAL(DP), INTENT(OUT) :: dnx
    
    !y-coordinate of normal vector
    REAL(DP), INTENT(OUT) :: dny
    
!</output>

!</subroutine>

    ! local variables
    INTEGER(I32), DIMENSION(:), POINTER :: p_IdbleSegInfo_handles,p_IintSegInfo_handles
    INTEGER(I32), DIMENSION(:), POINTER :: p_IsegInfo, p_IsegCount
    REAL(DP), DIMENSION(:), POINTER     :: p_DsegInfo, p_DmaxPar
    INTEGER :: cpar ! local copy of cparType
    INTEGER :: cnormalMeanCalc

    REAL(DP) :: dpar, dcurrentpar, dparloc, dendpar, dseglength, dnorm, dnx0, dny0
    INTEGER :: iseg,isegtype,istartidx

    cpar = BDR_PAR_01
    IF (PRESENT(cparType)) cpar = cparType

    cnormalMeanCalc = BDR_NORMAL_MEAN
    IF (PRESENT(cnormalMean)) cnormalMeanCalc = cnormalMean

    if ((iboundCompIdx.gt.rboundary%iboundarycount).or.(iboundCompIdx.lt.0)) then
      CALL output_line ('iboundCompIdx out of bounds!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'boundary_getNormalVec')
      CALL sys_halt()
    ENDIF

    ! Get the pointers to the segment information arrays for the current
    ! boundary component:
    CALL storage_getbase_int(rboundary%h_Iintdatavec_handles,p_IintSegInfo_handles)
    CALL storage_getbase_int(INT(p_IintSegInfo_handles(iboundCompIdx)),p_IsegInfo)

    CALL storage_getbase_int(rboundary%h_Idbldatavec_handles,p_IdbleSegInfo_handles)
    CALL storage_getbase_double(INT(p_IdbleSegInfo_handles(iboundCompIdx)),p_DsegInfo)

    ! Get the segment-count array and the maximum-parameter array
    CALL storage_getbase_int(rboundary%h_IsegCount,p_IsegCount)
    CALL storage_getbase_double(rboundary%h_DmaxPar,p_DmaxPar)

    ! Normalise the parameter value to the range [0,TMAX)
    dpar = dt
    CALL boundary_normaliseParValue2D(p_IsegCount,p_DmaxPar,iboundCompIdx,cpar,dpar)

    ! Find segment iseg the parameter value belongs to.
    ! Remember that in the first element in the double precision block of
    ! each segment, the length of the segment is noted!
    CALL boundary_getSegmentInfo2D(&
      p_IsegCount,p_DmaxPar,p_IsegInfo,p_DsegInfo,iboundCompIdx,dpar,cpar,0,&
      iseg,istartidx,dcurrentpar,dendpar,dseglength,dparloc,isegtype)
    
    IF (dseglength .EQ. 0.0_DP) dseglength = 1.0_DP ! trick to avoid div/0

    ! If the local parameter value of dt is <> 0, then the normal vector can 
    ! be computed uniquely at boundary point dt.
    ! Otherwise, we have to average the (two) normal vector(s) of the adjacent
    ! boundary components. Check if we are in the simple case.
    IF (dparloc .NE. 0.0_DP) THEN
    
      CALL getNormal (isegType,cpar,dparLoc,dseglength,istartidx,p_DsegInfo,dnx,dny)

    ELSE

      ! Ok, we are at the endpoint of the interval. Now, cnormalMean decides on
      ! what to do.
      dnx = 0.0_DP
      dny = 0.0_DP
      IF ((cnormalMeanCalc .EQ. BDR_NORMAL_MEAN) .OR. &
          (cnormalMeanCalc .EQ. BDR_NORMAL_RIGHT)) THEN
          
        ! Calculate the 'right' normal into dnx/dny.
        !
        ! We already have the segment 'right' to the point.
        ! Get the corresponding normal vector.
        
        CALL getNormal (isegType,cpar,dparLoc,dseglength,istartidx,p_DsegInfo,dnx,dny)
        
      END IF

      IF ((cnormalMeanCalc .EQ. BDR_NORMAL_MEAN) .OR. &
          (cnormalMeanCalc .EQ. BDR_NORMAL_LEFT)) THEN
          
        ! Calculate the 'left' normal into dnx/dny.
        !
        ! Find segment iseg previous to the parameter value belongs to.
        ! Remember that in the first element in the double precision block of
        ! each segment, the length of the segment is noted!

        CALL boundary_getSegmentInfo2D(&
          p_IsegCount,p_DmaxPar,p_IsegInfo,p_DsegInfo,iboundCompIdx,dpar,cpar,1,&
          iseg,istartidx,dcurrentpar,dendpar,dseglength,dparloc,isegtype)
        
        ! Calculate the normal into dnx0/dny0
        CALL getNormal (isegType,cpar,dparLoc,dseglength,istartidx,p_DsegInfo,dnx0,dny0)
        
        ! and add it to dnx/dny
        dnx = dnx+dnx0
        dny = dny+dny0
      END IF

      ! Normalise the vector -- in case two vectors are summed up.      
      dnorm = SQRT(dnx*dnx+dny*dny)
      dnx = dnx/dnorm
      dny = dny/dnorm
    END IF

  CONTAINS
  
    SUBROUTINE getNormal (isegType,cpar,dparLoc,dseglength,istartidx,DsegInfo,dnx,dny)
      
    ! Calculates the normal of a point on a specific boundary segment.
    
    !<input>
    
    ! Segment type
    INTEGER, INTENT(IN) :: isegType
    
    ! Type of the parameter value (0-1, length par.,...)
    INTEGER, INTENT(IN) :: cpar
    
    ! Local parameter value of the point in the boundary segment
    REAL(DP), INTENT(IN) :: dparLoc
    
    ! Length of the segment
    REAL(DP), INTENT(IN) :: dseglength
    
    ! Start index of the segment in DsegInfo
    INTEGER, INTENT(IN) :: istartIdx
    
    ! Double precision segment info array
    REAL(DP), DIMENSION(:), INTENT(IN) :: DsegInfo
    
    !</input>
    
    !<output>
    
    ! Normal vector
    REAL(DP), INTENT(OUT) :: dnx,dny
    
    !</output>
    
      ! local variables
      REAL(DP) :: dploc,dphi,dnorm
      
      dploc = dparloc
  
      SELECT CASE (isegType)

      ! case of line
      CASE (BOUNDARY_TYPE_LINE)

        ! Calculate the x/y components of the normal vector from the 
        ! startpoint and the uni direction vector.
        dnx0 =  p_DsegInfo(istartidx+6)
        dny0 = -p_DsegInfo(istartidx+5)

      ! case of circle segment
      CASE (BOUNDARY_TYPE_CIRCLE)

        ! Rescale dparloc with the length of the arc to get a value
        ! between 0 and 1; important for sin/cos functions later.
        ! In the 0-1 parametrisation, this is already the case.
        IF (cpar .EQ. BDR_PAR_LENGTH) dploc = dploc / dseglength

        ! Get the rotation angle.
        ! Use the initial rotation angle, saved at position 6 of the double
        ! precision data block
        dphi = p_DsegInfo(istartidx+7) &
            + dploc * (p_DsegInfo(istartidx+8)-p_DsegInfo(istartidx+7))

        ! And calculate the x/y components of the normal vector with sin/cos;
        ! the radius is to be found in element 5 of the double precision data block
        dnx0 = -p_DsegInfo(istartidx+5)*cos(dphi)
        dny0 = -p_DsegInfo(istartidx+5)*sin(dphi)

      CASE DEFAULT
        CALL output_line ('Wrong segment type: isegType='//sys_siL(isegType,10), &
                          OU_CLASS_ERROR,OU_MODE_STD,'boundary_getNormalVec')
        CALL sys_halt()
      END SELECT

      ! Normalize the vector
      dnorm = SQRT(dnx0*dnx0+dny0*dny0)
      dnx = dnx0/dnorm
      dny = dny0/dnorm
      
    END SUBROUTINE  

  END SUBROUTINE boundary_getNormalVec

END MODULE 
