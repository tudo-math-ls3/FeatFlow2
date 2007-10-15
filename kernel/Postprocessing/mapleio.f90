!#########################################################################
!# ***********************************************************************
!# <name> mapleio </name>
!# ***********************************************************************
!#
!# <purpose>
!# This module contains different routines for MAPLE export.
!# This feature can be used e.g. for debugging purposes to quickly analyse
!# situaltions inside of the program.
!#
!# The following routines can be found in this module:
!#
!# 1.) mapleio_writePointArray
!#     -> Writes a MAPLE command to an output channel for creating a 
!#        point array
!#
!# 2.) mapleio_writePolygonPlot
!#     -> Writes a MAPLE PLOT command to an output channel for creating a 
!#        polygon
!#
!# 3.) mapleio_writePointPlot
!#     -> Writes a MAPLE PLOT command to an output channel for creating a 
!#        point set
!'
!# </purpose>
!#########################################################################

MODULE mapleio

  USE fsystem
  USE storage
  USE io
  USE basicgeometry
  
  IMPLICIT NONE

  INTERFACE mapleio_writePointPlot
    MODULE PROCEDURE mapleio_writePointPlotSingle
    MODULE PROCEDURE mapleio_writePointPlotMult
  END INTERFACE

CONTAINS

  ! ***************************************************************************

!<subroutine>
  SUBROUTINE mapleio_writePointArray (ichannel,Dpoints,sname,iindex)
  
  !<description>
    ! Writes a MAPLE command to output channel ichannel which creates
    ! a MAPLE array containing the points in Dpoints. The MAPLE variable
    ! name will get the name 'sname'. If
    ! iindex is specified, the array gets the name 'sname(iindex)'.
  !</description>
    
  !<input>
    ! Output channel number where to write the MAPLE output to.
    INTEGER, INTENT(IN) :: ichannel
    
    ! Point set that should be written out. This is a list of (2D or 2D)
    ! point coordinates that form the polygon.
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dpoints
    
    ! Name of the variable to be used by MAPLE.
    CHARACTER(len=*), INTENT(IN) :: sname
    
    ! OPTIONAL: Index of the variable (if sname should be treated as
    ! an array of polygons). If present, the variable gets the name
    ! 'sname(iindex)'.
    INTEGER, INTENT(IN), OPTIONAL :: iindex
  !</input>
    
!</subroutine>

    INTEGER :: ipoint,idim
    
    ! Create the first line
    IF (.NOT. PRESENT(iindex)) THEN
      WRITE (ichannel,'(A)') sname//':=['
    ELSE
      WRITE (ichannel,'(A)') sname//'('//TRIM(sys_siL(iindex,10))//'):=['
    END IF
    
    ! Write the points
    DO ipoint=1,UBOUND(Dpoints,2)
    
      WRITE (ichannel,'(A)',ADVANCE='NO') &
        '     ['//TRIM(sys_sdL(Dpoints(1,ipoint),10))
        
      DO idim=2,UBOUND(Dpoints,1)
        WRITE (ichannel,'(A)',ADVANCE='NO') &
          ','//TRIM(sys_sdL(Dpoints(idim,ipoint),10))
      END DO
      
      IF (ipoint .EQ. UBOUND(Dpoints,2)) THEN
        WRITE (ichannel,'(A)',ADVANCE='NO') ']'
      ELSE
        WRITE (ichannel,'(A)',ADVANCE='YES') '],'
      END IF
      
    END DO

    WRITE (ichannel,'(A)') '];';
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  SUBROUTINE mapleio_writePolygonPlot (ichannel,Dpoints,sname,iindex)
  
  !<description>
    ! Writes a MAPLE command to output channel ichannel which creates
    ! a MAPLE PLOT command for a polygon defined  by the points in Dpoints. 
    ! The MAPLE variable name will get the name 'sname'. If
    ! iindex is specified, the variable gets the name 'sname(iindex)'.
  !</description>
    
  !<input>
    ! Output channel number where to write the MAPLE output to.
    INTEGER, INTENT(IN) :: ichannel
    
    ! Polygon that should be created. This is a list of (2D or 2D)
    ! point coordinates that form the polygon.
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dpoints
    
    ! Name of the variable to be used by MAPLE.
    CHARACTER(len=*), INTENT(IN) :: sname
    
    ! OPTIONAL: Index of the variable (if sname should be treated as
    ! an array of polygons). If present, the variable gets the name
    ! 'sname(iindex)'.
    INTEGER, INTENT(IN), OPTIONAL :: iindex
  !</input>
    
!</subroutine>

    INTEGER :: ipoint,idim
    
    ! Create the first line
    IF (.NOT. PRESENT(iindex)) THEN
      WRITE (ichannel,'(A)') sname//':=polygonplot(['
    ELSE
      WRITE (ichannel,'(A)') sname//'('//TRIM(sys_siL(iindex,10))//&
        '):=polygonplot(['
    END IF
    
    ! Write the points
    DO ipoint=1,UBOUND(Dpoints,2)
    
      WRITE (ichannel,'(A)',ADVANCE='NO') &
        '     ['//TRIM(sys_sdL(Dpoints(1,ipoint),10))
        
      DO idim=2,UBOUND(Dpoints,1)
        WRITE (ichannel,'(A)',ADVANCE='NO') &
          ','//TRIM(sys_sdL(Dpoints(idim,ipoint),10))
      END DO
      
      IF (ipoint .EQ. UBOUND(Dpoints,2)) THEN
        WRITE (ichannel,'(A)',ADVANCE='NO') ']'
      ELSE
        WRITE (ichannel,'(A)',ADVANCE='YES') '],'
      END IF
      
    END DO

    WRITE (ichannel,'(A)')']);';
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  SUBROUTINE mapleio_writePointPlotSingle (ichannel,Dpoint,sname,iindex)
  
  !<description>
    ! Writes a MAPLE command to output channel ichannel which creates
    ! a MAPLE PLOT command for a point defined by Dpoint. 
    ! The MAPLE variable name will get the name 'sname'. If
    ! iindex is specified, the variable gets the name 'sname(iindex)'.
  !</description>
    
  !<input>
    ! Output channel number where to write the MAPLE output to.
    INTEGER, INTENT(IN) :: ichannel
    
    ! Point that should be written out.
    REAL(DP), DIMENSION(:), INTENT(IN) :: Dpoint
    
    ! Name of the variable to be used by MAPLE.
    CHARACTER(len=*), INTENT(IN) :: sname
    
    ! OPTIONAL: Index of the variable (if sname should be treated as
    ! an array of polygons). If present, the variable gets the name
    ! 'sname(iindex)'.
    INTEGER, INTENT(IN), OPTIONAL :: iindex
  !</input>
    
!</subroutine>

    INTEGER :: idim
    
    ! Create the first line
    IF (.NOT. PRESENT(iindex)) THEN
      WRITE (ichannel,'(A)') sname//':=pointplot('
    ELSE
      WRITE (ichannel,'(A)') sname//'('//TRIM(sys_siL(iindex,10))//&
        '):=polygonplot('
    END IF
    
    ! Write the points
    WRITE (ichannel,'(A)',ADVANCE='NO') &
      '     ['//TRIM(sys_sdL(Dpoint(1),10))
      
    DO idim=2,SIZE(Dpoint)
      WRITE (ichannel,'(A)',ADVANCE='NO') &
        ','//TRIM(sys_sdL(Dpoint(idim),10))
    END DO
    
    WRITE (ichannel,'(A)')']);';
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  SUBROUTINE mapleio_writePointPlotMult (ichannel,Dpoints,sname,iindex)
  
  !<description>
    ! Writes a MAPLE command to output channel ichannel which creates
    ! a MAPLE PLOT command for a point set defined  by the points in Dpoints. 
    ! The MAPLE variable name will get the name 'sname'. If
    ! iindex is specified, the variable gets the name 'sname(iindex)'.
  !</description>
    
  !<input>
    ! Output channel number where to write the MAPLE output to.
    INTEGER, INTENT(IN) :: ichannel
    
    ! Polygon that should be created. This is a list of (2D or 2D)
    ! point coordinates that form the polygon.
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dpoints
    
    ! Name of the variable to be used by MAPLE.
    CHARACTER(len=*), INTENT(IN) :: sname
    
    ! OPTIONAL: Index of the variable (if sname should be treated as
    ! an array of polygons). If present, the variable gets the name
    ! 'sname(iindex)'.
    INTEGER, INTENT(IN), OPTIONAL :: iindex
  !</input>
    
!</subroutine>

    INTEGER :: ipoint,idim
    
    ! Create the first line
    IF (.NOT. PRESENT(iindex)) THEN
      WRITE (ichannel,'(A)') sname//':=pointplot('
    ELSE
      WRITE (ichannel,'(A)') sname//'('//TRIM(sys_siL(iindex,10))//&
        '):=polygonplot(['
    END IF
    
    ! Write the point
    WRITE (ichannel,'(A)',ADVANCE='NO') &
      '     ['//sys_sd(Dpoints(1,ipoint),10)
      
    DO idim=2,UBOUND(Dpoints,1)
      WRITE (ichannel,'(A)',ADVANCE='NO') &
        ','//sys_sd(Dpoints(idim,ipoint),10)
    END DO
    
    IF (ipoint .EQ. UBOUND(Dpoints,2)) THEN
      WRITE (ichannel,'(A)',ADVANCE='YES') '],'
    ELSE
      WRITE (ichannel,'(A)',ADVANCE='NO') ']'
    END IF

    WRITE (ichannel,'(A)')');';
  
  END SUBROUTINE

END MODULE
