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
!# 4.) boundary_dgetNsegments
!#     -> Get the number of boundary segments in a boundary component
!#
!# 5.) boundary_createRegion
!#     -> Get the characteristics of a boundary segment and create
!#        a boundary region structure from it.
!#
!# 6.) boundary_isInRegion
!#     -> Tests whether a node with a specific parameter value
!#        is inside of a given boundary region or not.
!# </purpose>
!##############################################################################

MODULE boundary

  USE genoutput
  USE storage
  USE fsystem
  USE error
  USE io

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

!<constantblock description="Bitfield constants for t_boundaryRegion%iproperties">

  ! The startpoint (point with min. parameter value) belongs
  ! to the boundary segment. If not set, the point with the
  ! min. parameter value does not belong to the boundary segment.
  INTEGER, PARAMETER :: BDR_PROP_WITHSTART = 2**0

  ! The endpoint (point with min. parameter value) belongs
  ! to the boundary segment. If not set, the point with the
  ! max. parameter value does not belong to the boundary segment.
  INTEGER, PARAMETER :: BDR_PROP_WITHEND = 2**1

!</constantblock>

!</constants>

!<types>

!<typeblock>

  ! Structure defining a boundary region with minimum/maximum parameter
  ! value in 2D.
  TYPE t_boundaryRegion
  
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
    INTEGER(I32) :: iproperties = BDR_PROP_WITHSTART + BDR_PROP_WITHEND
  
  END TYPE
  
!</typeblock>

!<typeblock>

  ! Boundary structure of the domain
  TYPE t_boundary
  
    PRIVATE 

    ! number of fictitious boundary components
    INTEGER :: iboundarycount_f = -1

    ! number of geometric boundary components
    INTEGER :: iboundarycount_g = -1

    ! total number of boundary components
    INTEGER :: iboundarycount   = -1

    ! handle to double precision array: For every boundary component, maximum
    ! parameter value.
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

    IF ((iboundCompIdx.gt.rboundary%iboundarycount).or.(iboundCompIdx.lt.0)) THEN
      PRINT *,'Warning in boundary_dgetLength'
      dresult = -1
      RETURN
    ENDIF

    ip1=ceiling(dt1)
    ip2=floor(dt2)

    if (ip2.lt.ip1) then

      call boundary_getcoords(rboundary,iboundCompIdx,dt1,dx1,dy1)
      call boundary_getcoords(rboundary,iboundCompIdx,dt2,dx2,dy2)

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
          call boundary_getcoords(rboundary,iboundCompIdx,dt1,dx2,dy2)
          dbl=0.0_DP
        else
          call boundary_getcoords(rboundary,iboundCompIdx,dt1,dx1,dy1)
          call boundary_getcoords(rboundary,iboundCompIdx,real(ip1,DP),dx2,dy2)
          dbl=sqrt((dx2-dx1)*(dx2-dx1)+(dy2-dy1)*(dy2-dy1))
        endif

        do i=ip1+1,ip2

          dx1=dx2
          dy1=dy2

          call boundary_getcoords(rboundary,iboundCompIdx,real(i,DP),dx2,dy2)
          dbl=dbl+sqrt((dx2-dx1)*(dx2-dx1)+(dy2-dy1)*(dy2-dy1))

        enddo

        if (real(ip2,DP).ne.dt2) then
          call boundary_getcoords(rboundary,iboundCompIdx,dt2,dx1,dy1)
          dbl=dbl+sqrt((dx1-dx2)*(dx1-dx2)+(dy1-dy2)*(dy1-dy2))
        endif

      endif

    endif

    dresult=dbl

  END FUNCTION 


!************************************************************************

!<function>

  REAL(DP) FUNCTION boundary_dgetMaxParVal(rboundary, iboundCompIdx)

  !<description>
  ! This function returns the parametric length of the boundary component iboundCompIdx
  !</description>

  !<result>
  ! Parametric length of boundary component iboundCompIdx.
  !</result>

  !<input>

  !boundary structure
  TYPE(t_boundary), INTENT(IN) :: rboundary

  !index of boundary component
  INTEGER, INTENT(IN) :: iboundCompIdx
  
  !</input>
  
!</function>

  REAL(DP),DIMENSION(:),POINTER :: p_DmaxPar

  !if iboundCompIdx exceeds the total number of boundary components or is negative, abort
  if ((iboundCompIdx .gt. rboundary%iboundarycount) .or. (iboundCompIdx.lt.0)) then
    PRINT *,'Warning in boundary_dgetMaxParVal'
    boundary_dgetMaxParVal = -1
    RETURN
  ENDIF

  !get vector with component length
  CALL storage_getbase_double(rboundary%h_DmaxPar,p_DmaxPar)
  
  ! Get the maximum parameter value
  boundary_dgetMaxParVal = p_DmaxPar(iboundCompIdx)

  END FUNCTION 

!************************************************************************

!<function>

  INTEGER FUNCTION boundary_dgetNsegments(rboundary, iboundCompIdx)

  !<description>
  ! This function returns the number of boundary segments in
  ! the boundary component iboundCompIdx.
  !</description>

  !<result>
  ! Parametric length of boundary component iboundCompIdx.
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
    PRINT *,'Warning in boundary_dgetNsegments'
    boundary_dgetNsegments = -1
    RETURN
  ENDIF

  !get vector with component length
  CALL storage_getbase_int(rboundary%h_IsegCount,p_IsegCount)
  
  ! Get the maximum parameter value
  boundary_dgetNsegments = p_IsegCount(iboundCompIdx)

  END FUNCTION 

!************************************************************************

!<subroutine>

  SUBROUTINE boundary_read_prm(rboundary, sfilename)

  !<description>
  ! This routine reads a .PRM file into memory. The boundary structure
  ! rboundary is initialised with the data from the file.
  ! The parameter sfilename gives the name of the .prm file to read.
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
  INTEGER :: ibcomponent, isegment, ibct
  INTEGER(I32), DIMENSION(:), POINTER :: p_IsegInfo, p_IsegCount
  INTEGER(I32), DIMENSION(:), POINTER :: p_IdbleSegInfo_handles,p_IintSegInfo_handles
  REAL(DP), DIMENSION(:), POINTER :: p_DsegInfo, p_DmaxPar
  INTEGER :: ityp, nspline, npar
  INTEGER :: isegrel
  INTEGER :: idblemem  ! Counts the memory we need
  REAL(DP) :: dl,dmaxpar
  
  ! Open the file
  CALL io_openFileForReading(sfilename, iunit)
  
  ! Read "NBCT"
  READ (iunit,*)
  
  ! Read NBCT - Number of boundary components
  READ (iunit,*) rboundary%iboundarycount_g
  
  ! Fictitious boundary components not supported at the moment
  rboundary%iboundarycount_f = 0
  
  rboundary%iboundarycount = rboundary%iboundarycount_g
  
  ! Allocate an array containing handles. Each handle Each handle refers
  ! to integer data for a boundary component.
  CALL storage_new1D("boundary_read", "h_Idbldatavec_handles", &
                  rboundary%iboundarycount, ST_INT, &
                  rboundary%h_Idbldatavec_handles, ST_NEWBLOCK_ZERO)
  CALL storage_getbase_int(rboundary%h_Idbldatavec_handles, p_IdbleSegInfo_handles)
  
  ! Allocate an array containing of handles. Each handle refers
  ! to integer data for a boundary component.
  CALL storage_new("boundary_read", "h_Iintdatavec_handles", &
                  rboundary%iboundarycount, ST_INT, &
                  rboundary%h_Iintdatavec_handles, ST_NEWBLOCK_ZERO)
  CALL storage_getbase_int(rboundary%h_Iintdatavec_handles, p_IintSegInfo_handles)

  ! Allocate an array containing the maximum parameter values for each
  ! boundary component.
  CALL storage_new("boundary_read", "h_DmaxPar", &
                  rboundary%iboundarycount, ST_DOUBLE, &
                  rboundary%h_DmaxPar, ST_NEWBLOCK_ZERO)
  CALL storage_getbase_double(rboundary%h_DmaxPar, p_DmaxPar)

  ! Allocate an array containing the number of boundary segments in each
  ! boundary component
  CALL storage_new("boundary_read", "h_IsegCount", &
                  rboundary%iboundarycount, ST_INT, &
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
                    2*p_IsegCount(ibct), ST_INT, &
                    p_IintSegInfo_handles(ibct), ST_NEWBLOCK_ZERO)
    CALL storage_getbase_int(p_IintSegInfo_handles(ibct), p_IsegInfo)

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
                    idblemem, ST_DOUBLE, &
                    p_IdbleSegInfo_handles(ibct), ST_NEWBLOCK_ZERO)
    
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
    CALL storage_getbase_double(p_IdbleSegInfo_handles(ibcomponent), p_DsegInfo)
    
    ! Get a pointer to the start position of the segments in the
    ! current boundary component.
    CALL storage_getbase_int(p_IintSegInfo_handles(ibcomponent), p_ISegInfo)

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
        ! - normalised direction vector
        ! - length of the line.
        !
        ! Startpoint
        READ(iunit,*) p_DsegInfo(isegrel+3),p_DsegInfo(isegrel+4)
        
        ! Relative endpoint
        READ(iunit,*) p_DsegInfo(isegrel+5),p_DsegInfo(isegrel+6)
        
        ! Save the initial parameter value in the first entry
        p_DsegInfo(isegrel+1) = dmaxpar
        
        ! Calculate the length and save it in the first position
        dl = SQRT(p_DsegInfo(isegrel+5)**2 + p_DsegInfo(isegrel+6)**2)
        p_DsegInfo(isegrel+2) = dl
        
        ! Normalise the direction vector
        p_DsegInfo(isegrel+4:isegrel+6) = p_DsegInfo(isegrel+4:isegrel+7)/dl
        
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
        
        ! Save the initial parameter value in the first entry
        p_DsegInfo(isegrel+1) = dmaxpar
        
        ! Now compute the real length of the arc.
        dl = p_DsegInfo(isegrel+4) * &
             ABS(p_DsegInfo(isegrel+7)-p_DsegInfo(isegrel+6))
        p_DsegInfo(isegrel+2) = dl
        
        ! Increase the maximum parameter value
        dmaxPar = dmaxPar + dl
        
      END SELECT
    END DO
    
    ! Save the maximum parameter value for that component
    p_DmaxPar(ibcomponent) = dmaxPar
    
  END DO
    
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
  INTEGER :: i
  INTEGER(I32), DIMENSION(:), POINTER :: p_IdbleSegInfo_handles,p_IintSegInfo_handles
  
  ! Get the pointers to the segment information arrays for the current
  ! boundary component:
  CALL storage_getbase_int(rboundary%h_Iintdatavec_handles,p_IintSegInfo_handles)
  CALL storage_getbase_int(rboundary%h_Idbldatavec_handles,p_IdbleSegInfo_handles)

  ! Release the handles of the integer- and double-precision
  ! data blocks:
  DO i=1,rboundary%iboundarycount
    CALL storage_free (p_IintSegInfo_handles(i))
    CALL storage_free (p_IdbleSegInfo_handles(i))
  END DO
  
  ! Release all arrays in the structure
  CALL storage_free (rboundary%h_Iintdatavec_handles)
  CALL storage_free (rboundary%h_Idbldatavec_handles)
  CALL storage_free (rboundary%h_IsegCount)
  CALL storage_free (rboundary%h_DmaxPar)
  
  rboundary%iboundarycount_f = 0
  rboundary%iboundarycount_g = 0
  rboundary%iboundarycount = 0

  ! That's it...

  END SUBROUTINE

!************************************************************************

!<subroutine>

  SUBROUTINE boundary_getcoords(rboundary, iboundCompIdx, dt, dx, dy)

  !<description>
  ! This routine returns for a given parameter value dt the
  ! cartesian coordinates of the point on the boundary component 
  ! iboundCompIdx.\\
  ! dt is bounded by the maximum parameter value, i.e. when dt is
  ! larger than the maximum parameter value, dt=0 is taken.
  !</description>

  !<input>

  !boundary structure
  TYPE(t_boundary), INTENT(IN) :: rboundary

  !index of boundary component
  INTEGER, INTENT(IN) :: iboundCompIdx

  !parametric value of boundary point
  REAL(DP), INTENT(IN) :: dt
  
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
  
  REAL(DP) :: dpar, dcurrentpar, dparloc, dphi, dendpar
  INTEGER :: iseg,isegtype,istartidx

  if ((iboundCompIdx.gt.rboundary%iboundarycount).or.(iboundCompIdx.lt.0)) then
    PRINT *,'Error in boundary_getcoords'
    STOP
  ENDIF

  ! Get the pointers to the segment information arrays for the current
  ! boundary component:
  CALL storage_getbase_int(rboundary%h_Iintdatavec_handles,p_IintSegInfo_handles)
  CALL storage_getbase_int(p_IintSegInfo_handles(iboundCompIdx),p_IsegInfo)
  
  CALL storage_getbase_int(rboundary%h_Idbldatavec_handles,p_IdbleSegInfo_handles)
  CALL storage_getbase_double(p_IdbleSegInfo_handles(iboundCompIdx),p_DsegInfo)
  
  ! Get the segment-count array and the maximum-parameter array
  CALL storage_getbase_int(rboundary%h_IsegCount,p_IsegCount)
  CALL storage_getbase_double(rboundary%h_DmaxPar,p_DmaxPar)

  !if the parameter value exceeds the parameter interval on the boundary component
  
  IF (dt .GE. p_DmaxPar(iboundCompIdx) ) THEN
    dpar = 0.0 
  ELSE
    dpar = dt
  ENDIF

  ! Find segment iseg the parameter value belongs to.
  ! Remember that in the first element in the double precision block of
  ! each segment, the length of the segment is noted!
  
  dcurrentpar = 0.0_DP
  
  DO iseg = 0,p_IsegCount(iboundCompIdx)-1
    
    ! Determine Start index of the segment in the double-prec. block
    istartidx = p_IsegInfo(1+2*iseg+1) 
    
    ! Get the start and end parameter value
    dcurrentpar = p_DsegInfo(1+istartidx)
    dendpar = dcurrentpar + p_DsegInfo(2+istartidx)
    
    ! At least one of the IF-commands in the loop will activate
    ! the exit - because of the 'dt' check above!
    IF (dpar .LT. dendpar) EXIT
    
  END DO

  ! Subtract the start position of the current boundary component
  ! from the parameter value to get the 'local' parameter value
  ! (0 .le. dparloc .le. length(segment))
  dparloc = dpar - dcurrentpar

  ! Use the segment type to determine how to calculate
  ! the coordinate. Remember that the segment type is noted
  ! in the first element of the integer block of each segment!

  isegtype = p_IsegInfo(1+2*iseg)

  SELECT CASE (isegType)

  ! case of line
  CASE (BOUNDARY_TYPE_LINE)
  
    ! Calculate the x/y coordinates from the startpoint and
    ! the unid direction vector.
    dx = p_DsegInfo(istartidx+3) + dparloc*p_DsegInfo(istartidx+5)
    dy = p_DsegInfo(istartidx+4) + dparloc*p_DsegInfo(istartidx+6)

  ! case of circle segment
  CASE (BOUNDARY_TYPE_CIRCLE)

    ! Rescale dparloc with the length of the arc to get a value
    ! between 0 and 1; important for sin/cos functions later
    dparloc = dparloc / p_DsegInfo(istartidx+2)

    ! Get the rotation angle.
    ! Use the initial rotation angle, saved at position 6 of the double
    ! precision data block
    dphi = p_DsegInfo(istartidx+7) &
         + dparloc * (p_DsegInfo(istartidx+8)-p_DsegInfo(istartidx+7))

    ! And calculate the x/y coordinate with sin/cos; the radius is
    ! to be found in element 4 of the double precision data block!
    ! The center of the circle is at position 2/3.
    dx = p_DsegInfo(istartidx+3) + p_DsegInfo(istartidx+5)*cos(dphi)
    dy = p_DsegInfo(istartidx+4) + p_DsegInfo(istartidx+5)*sin(dphi)

  CASE DEFAULT
    PRINT *,'boundary_getcoords: Wrong segment type'
    STOP
  END SELECT

  END SUBROUTINE 

!************************************************************************

!<subroutine>

  SUBROUTINE boundary_createRegion (rboundary, iboundCompIdx, iboundSegIdx, &
                                    rregion)
                                    
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

  ! index of boundary component
  INTEGER, INTENT(IN) :: iboundCompIdx

  ! index of the boundary segment
  INTEGER, INTENT(IN) :: iboundSegIdx
  
  !</input>

  !<output>

  ! Boundary region that is characterised by the boundary segment
  TYPE(t_boundaryRegion), INTENT(OUT) :: rregion
    
  !</output>
  
!</subroutine>

  ! local variables
  INTEGER(I32), DIMENSION(:), POINTER :: p_IdbleSegInfo_handles,p_IintSegInfo_handles
  INTEGER(I32), DIMENSION(:), POINTER :: p_IsegInfo, p_IsegCount
  REAL(DP), DIMENSION(:), POINTER     :: p_DsegInfo, p_DmaxPar
  
  REAL(DP) :: dcurrentpar, dendpar
  INTEGER :: istartidx

  IF ((iboundCompIdx.gt.rboundary%iboundarycount).or.(iboundCompIdx.lt.0)) THEN
    PRINT *,'Error in boundary_createregion'
    STOP
  ENDIF

  ! Get the pointers to the segment information arrays for the current
  ! boundary component:
  CALL storage_getbase_int(rboundary%h_Iintdatavec_handles,p_IintSegInfo_handles)
  CALL storage_getbase_int(p_IintSegInfo_handles(iboundCompIdx),p_IsegInfo)
  
  CALL storage_getbase_int(rboundary%h_Idbldatavec_handles,p_IdbleSegInfo_handles)
  CALL storage_getbase_double(p_IdbleSegInfo_handles(iboundCompIdx),p_DsegInfo)
  
  ! Get the segment-count array and the maximum-parameter array
  CALL storage_getbase_int(rboundary%h_IsegCount,p_IsegCount)
  CALL storage_getbase_double(rboundary%h_DmaxPar,p_DmaxPar)

  IF ((iboundCompIdx.gt.p_IsegCount(iboundCompIdx)).or.(iboundSegIdx.lt.0)) THEN
    PRINT *,'Error in boundary_createregion: iboundSegIdx out of bounds.'
    STOP
  ENDIF

  ! Find segment iseg the parameter value belongs to.
  ! Remember that in the first element in the double precision block of
  ! each segment, the length of the segment is noted!
  
  istartidx = p_IsegInfo(1+2*0+1) 
  dcurrentpar = 0.0_DP
  
  ! Determine Start index of the segment in the double-prec. block
  istartidx = p_IsegInfo(1+2*(iboundSegIdx-1)+1) 
  
  ! Get the start and end parameter value
  dcurrentpar = p_DsegInfo(1+istartidx)
  dendpar = dcurrentpar + p_DsegInfo(2+istartidx)
  
  ! Create the boundary region structure
  rregion%dminParam = dcurrentpar
  rregion%dmaxParam = dendpar
  rregion%ctype = BDR_TP_CURVE
  rregion%iproperties = BDR_PROP_WITHSTART
  rregion%isegmentType = p_IsegInfo(1+2*(iboundSegIdx-1))
  rregion%iboundCompIdx = iboundCompIdx
  rregion%iboundSegIdx = iboundSegIdx
  rregion%dmaxParamBC = p_DmaxPar(iboundCompIdx)

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
  ! Must be in the range 0..max. par. value
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
 
  END MODULE 

