!##############################################################################
!# ****************************************************************************
!# <name> transformation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines for the transformation between a reference
!# element and a real element, for triangular and quadrilateral 2D elements.
!# (3D support will come later...)
!#
!# 1.) trafo_calcJacPrepare
!#     -> calculates auxiliary Jacobian factors for the transformation
!#        from a reference quadrilateral to a real quadrilateral
!#
!# 2.) trafo_calcJac
!#     -> calculates the Jacobian matrix + Jacobian determinant of the mapping
!#        from  the reference to a real quadrilateral element
!#
!# 3.) trafo_calcRealCoords
!#     -> maps a point from the reference element to the real element
!#
!# 4.) trafo_calcTrafo
!#     -> calculates the Jacobian matrix + Jacobian determinant of the mapping
!#        from  the reference to a real quadrilateral element
!#      -> maps a point from the reference element to the real element
!#     (so performing the same task as elem_calcJac and elem_calcRealCoords
!#      in one routine)
!#
!# </purpose>
!##############################################################################

MODULE transformation

  USE fsystem
  USE basicgeometry

  IMPLICIT NONE
  
!<constants>
!<constantblock description="maximal values">

  ! (maximum) number of vertices per element, supported by the transformation
  INTEGER, PARAMETER :: TRAFO_MAXNVE = 4

  ! Number of entries in the Jacobian matrix, defining the mapping between
  ! the reference element and the real element. 2x2-matrix=4 elements
  INTEGER, PARAMETER :: TRAFO_NJACENTRIES = 4

  ! Number of entries in the array with the auxiliary Jacobian factors
  INTEGER, PARAMETER :: TRAFO_NAUXJAC = 4

!</constantblock>


!<constantblock description="id values for the different types of transformation">

  ! Unspecified transformation
  INTEGER, PARAMETER :: TRAFO_IDUNKNOWN   = 0

  ! Linear triangular transformation
  INTEGER, PARAMETER :: TRAFO_IDLINTRI    = 1

  ! Bilinear quadrilateral transformation
  INTEGER, PARAMETER :: TRAFO_IDBILINQUAD = 2
!</constantblock>


!</constants>

CONTAINS

! **********************************************************************
! General trnasformation support for multiple points on multiple
! elements. Wrapper routines for all types of transformations.

!<subroutine>

  SUBROUTINE trafo_calctrafo_sim (ctrafoType,nelements,npointsPerEl,Dcoords,&
                                  DpointsRef,Djac,Ddetj,DpointsReal)

!<input>

  ! Type of transformation to calculate
  INTEGER, INTENT(IN) :: ctrafoType

  ! Number of elements where to calculate the transformation
  INTEGER, INTENT(IN) :: nelements
  
  ! Number of points in each element where to calculate the transformation
  INTEGER, INTENT(IN) :: npointsPerEl

  ! Coordinates of the corners of all the elements
  !  Dcoord(1,i,.) = x-coordinates of corner i on an element, 
  !  Dcoord(2,i,.) = y-coordinates of corner i on an element.
  REAL(DP), DIMENSION(NDIM2D,TRAFO_MAXNVE,nelements), INTENT(IN) :: Dcoords

  ! Coordinates of the points on the reference element for each element 
  ! where to calculate the mapping.
  !  DpointsRef(1,i,.) = x-coordinates of point i on an element, 
  !  DpointsRef(2,i,.) = y-coordinates of point i on an element.
  REAL(DP), DIMENSION(NDIM2D,npointsPerEl,nelements), INTENT(IN) :: DpointsRef
  
!</inbut>

!<output>
  
  ! The Jacobian matrix of the mapping for each point.
  REAL(DP), DIMENSION(TRAFO_NJACENTRIES,npointsPerEl,nelements), INTENT(OUT) :: Djac
  
  ! Jacobian determinants of the mapping for all the points from the
  ! reference element to the real element.
  REAL(DP), DIMENSION(npointsPerEl,nelements), INTENT(OUT) :: Ddetj
  
  ! OPTIONAL: Array receiving the coordinates of the points in DpointsRef,
  ! mapped from the reference element to the real element.
  ! If not specified, they are not computed.
  REAL(DP), DIMENSION(NDIM2D,npointsPerEl,nelements), INTENT(OUT), OPTIONAL :: DpointsReal

!</output>

!</subroutine>

  ! local variables
  INTEGER :: iel, ipt
  
  ! auxiliary factors for the bilinear quad mapping
  REAL(DP), DIMENSION(TRAFO_NAUXJAC) :: DjacPrep
  
  ! What type of transformation do we have?
  
  SELECT CASE (ctrafoType)
  
  CASE (TRAFO_IDLINTRI)
    PRINT *,'Triangular transformation currently not supported!'
    STOP
  
  CASE (TRAFO_IDBILINQUAD)
  
    ! Calculate with or without coordinates?
    IF (.NOT. PRESENT(DpointsReal)) THEN
    
      ! Loop over the elements
      DO iel = 1,nelements
        ! Prepare the calculation of the Jacobi determinants
        CALL trafo_calcJacPrepare2(Dcoords(:,:,iel), DjacPrep)
        
        ! Loop over the points
        DO ipt=1,npointsPerEl
          ! Calculate the Jacobian matrix and determinant
          CALL trafo_calcJac2 (Dcoords(:,:,iel),DjacPrep,Djac(:,ipt,iel),Ddetj(ipt,iel), &
                               DpointsRef(1,ipt,iel),DpointsRef(2,ipt,iel))
        END DO ! ipt
        
      END DO ! iel
    
    ELSE

      ! Loop over the elements
      DO iel = 1,nelements
        ! Prepare the calculation of the Jacobi determinants
        CALL trafo_calcJacPrepare2(Dcoords(:,:,iel), DjacPrep)
        
        ! Loop over the points
        DO ipt=1,npointsPerEl
          ! Calculate the Jacobian matrix and determinant
          CALL trafo_calcTrafo2 (Dcoords(:,:,iel),DjacPrep,Djac(:,ipt,iel),Ddetj(ipt,iel), &
                                 DpointsRef(1,ipt,iel),DpointsRef(2,ipt,iel), &
                                 DpointsReal(1,ipt,iel),DpointsReal(2,ipt,iel))
        END DO ! ipt
        
      END DO ! iel

    END IF
    
  END SELECT

  END SUBROUTINE

! **********************************************************************
! General transformation support, multiple points. Wrapper routines 
! for all types of transformations.

!<subroutine>

  SUBROUTINE trafo_calctrafo_mult (ctrafoType,npointsPerEl,Dcoords,&
                                   DpointsRef,Djac,Ddetj,DpointsReal)

!<input>

  ! Type of transformation to calculate
  INTEGER, INTENT(IN) :: ctrafoType

  ! Number of points in each element where to calculate the transformation
  INTEGER, INTENT(IN) :: npointsPerEl

  ! Coordinates of the corners of all the elements
  !  Dcoord(1,i) = x-coordinates of corner i on an element, 
  !  Dcoord(2,i) = y-coordinates of corner i on an element.
  REAL(DP), DIMENSION(NDIM2D,TRAFO_MAXNVE), INTENT(IN) :: Dcoords

  ! Coordinates of the points on the reference element for each element 
  ! where to calculate the mapping.
  !  DpointsRef(1,i) = x-coordinates of point i on an element, 
  !  DpointsRef(2,i) = y-coordinates of point i on an element.
  REAL(DP), DIMENSION(NDIM2D,npointsPerEl), INTENT(IN) :: DpointsRef
  
!</inbut>

!<output>
  
  ! The Jacobian matrix of the mapping for each point.
  REAL(DP), DIMENSION(TRAFO_NJACENTRIES,npointsPerEl), INTENT(OUT) :: Djac
  
  ! Jacobian determinants of the mapping for all the points from the
  ! reference element to the real element.
  REAL(DP), DIMENSION(npointsPerEl), INTENT(OUT) :: Ddetj
  
  ! OPTIONAL: Array receiving the coordinates of the points in DpointsRef,
  ! mapped from the reference element to the real element.
  ! If not specified, they are not computed.
  REAL(DP), DIMENSION(NDIM2D,npointsPerEl), INTENT(OUT), OPTIONAL :: DpointsReal

!</output>

!</subroutine>

  ! local variables
  INTEGER :: iel, ipt
  
  ! auxiliary factors for the bilinear quad mapping
  REAL(DP), DIMENSION(TRAFO_NAUXJAC) :: DjacPrep
  
  ! What type of transformation do we have?
  
  SELECT CASE (ctrafoType)
  
  CASE (TRAFO_IDLINTRI)
    PRINT *,'Triangular transformation currently not supported!'
    STOP
  
  CASE (TRAFO_IDBILINQUAD)
  
    ! Calculate with or without coordinates?
    IF (.NOT. PRESENT(DpointsReal)) THEN
    
      ! Prepare the calculation of the Jacobi determinants
      CALL trafo_calcJacPrepare2(Dcoords, DjacPrep)
      
      ! Loop over the points
      DO ipt=1,npointsPerEl
        ! Calculate the Jacobian matrix and determinant
        CALL trafo_calcJac2 (Dcoords,DjacPrep,Djac(:,ipt),Ddetj(ipt), &
                             DpointsRef(1,ipt),DpointsRef(2,ipt))
      END DO ! ipt
    
    ELSE

      ! Prepare the calculation of the Jacobi determinants
      CALL trafo_calcJacPrepare2(Dcoords, DjacPrep)
      
      ! Loop over the points
      DO ipt=1,npointsPerEl
        ! Calculate the Jacobian matrix and determinant
        CALL trafo_calcTrafo2 (Dcoords,DjacPrep,Djac(:,ipt),Ddetj(ipt), &
                                DpointsRef(1,ipt),DpointsRef(2,ipt), &
                                DpointsReal(1,ipt),DpointsReal(2,ipt))
      END DO ! ipt

    END IF
    
  END SELECT

  END SUBROUTINE

! **********************************************************************
! Explaination of the quadrilateral transformation:
!
! We want to perform a transformation from the reference quadrilateral
! [-1,1]x[-1,1] onto a "real" quadrilaterl:
!
!   (-1,1) ___ (1,1)         (x4,y4)   ___  (x3,y3)
!         |  |                        /  \
!         |__|           =>          /____\
!  (-1,-1)    (1,-1)         (x1,y1)        (x2,y2)   
!
! By theory this can by done with a bilinear mapping, i.e. a mapping
! Psi:R^2->R^2  of the form:
!
!      Psi (xi1) = ( a1 + a2*xi1 + a3*xi2 + a4*xi1*xi2 )
!          (xi2)   ( b1 + b2*xi1 + b3*xi2 + b4*xi1*xi2 )
!
! Our special transformation has to map:
!
!  Psi(-1,-1) = (x1,y1)
!  Psi( 1,-1) = (x2,y2)
!  Psi( 1, 1) = (x3,y3)
!  Psi(-1, 1) = (x4,y4)
!
! This gives the linear system:
!
!  a1 - a2 - a3 - a4 = x1       b1 - b2 - b3 - b4 = y1
!  a1 + a2 - a3 - a4 = x2       b1 + b2 - b3 - b4 = y2
!  a1 + a2 + a3 + a4 = x3       b1 + b2 + b3 + b4 = y3
!  a1 - a2 + a3 - a4 = x4       b1 - b2 + b3 - b4 = y4
!
! Reorder this to calculate the ai:
!
!  a1 = 1/4 * ( x1 + x2 + x3 + x4)      b1 = 1/4 * ( y1 + y2 + y3 + y4)
!  a2 = 1/4 * (-x1 + x2 + x3 - x4)      b2 = 1/4 * (-y1 + y2 + y3 - y4)
!  a3 = 1/4 * (-x1 - x2 + x3 + x4)      b3 = 1/4 * (-y1 - y2 + y3 + y4)
!  a4 = 1/4 * ( x1 - x2 + x3 - x4)      b4 = 1/4 * ( y1 - y2 + y3 - y4)
!
! The factors in the brackets in these equations are only dependent on 
! the corners of the quadrilateral, not of the current point. So they
! are constant for all points we want to map from the reference
! element to the real one. We call them here "auxiliary Jacobian factors".
! They can be calculated with elem_calcJacPrepare in advance for all points that
! have to be mapped. To be more exact, elem_calcJacPrepare calculates:
!
!  J1 = 1/2 * (-x1 - x2 + x3 - x4)
!  J2 = 1/2 * ( x1 - x2 + x3 - x4)
!
!  J3 = 1/2 * (-y1 + y2 - y3 + y4)
!  J4 = 1/2 * (-y1 + y2 + y3 - y4)
!
! Using these factors, one can write:
!
!  a1  =  1/4 * ( x1 + x2 + x3 + x4)  =  1/2 * (x1 + x2 + J1)
!  a2  =  1/4 * (-x1 + x2 + x3 - x4)  =  1/2 * (x2 - x1 + J3)
!  a3  =  1/4 * (-x1 - x2 + x3 + x4)  =  1/2 * J1 
!  a4  =  1/4 * ( x1 - x2 + x3 - x4)  =  1/2 * J2
!
!  b1  =  1/4 * ( y1 + y2 + y3 + y4)  =  1/2 * (y1 + y3 + J2)
!  b2  =  1/4 * (-y1 + y2 + y3 - y4)  =  1/2 * J4
!  b3  =  1/4 * (-y1 - y2 + y3 + y4)  =  1/2 * (y3 - y1 - J4)
!  b4  =  1/4 * ( y1 - y2 + y3 - y4)  =  1/2 * J3
!
! The Jacobian matrix of the bilinear transformation is now 
! calculated as usual by partial differentiation. Thing above 
! coefficients ai and bi one can write:
!
! DPhi ( xi1 )  =  ( a2 + a4*xi2    a3 + a4*xi1 )
!      ( xi2 )     ( b2 + b4*xi2    b3 + b4*xi1 )
!
!               =  ( 1/2*(x2-x1+J2) + 1/2*J2*xi2           1/2*J1 + 1/2*J2*xi1 )
!                  (         1/2*J4 + 1/2*J3*xi2   1/2*(y3-y1-J4) + 1/2*J3*xi1 )
!
! which gives the Jacobian determinant of the 2x2-matrix:
!
!   det DPhi (xi1,xi2)  =  (DPhi[1,1]*DPhi[2,2] - DPhi[1,2]*DPhi[2,1]) (xi1,xi2)
!
! Using these information makes it possible to map a point (XI1,XI2)
! on the reference element to coordinates (XX,YY) on the real element by:
!
!   (XX,YY) := Psi(XI1,XI2) 
! **********************************************************************

!<subroutine>

  PURE SUBROUTINE trafo_calcJacPrepare2(Dcoords, DjacPrep)
  
!<description>
  ! Calculate auxiliary Jacobian factors for standard quadrilateral.
  ! This routine builds up constant factors that are later used during
  ! the transformation from the reference element to the real element.
!</description>

  !<input>
  
  ! Coordinates of the four corners of the real quadrilateral.
  ! Dcoord(1,.) = x-coordinates, 
  ! Dcoord(2,.) = y-coordinates.
  REAL(DP), DIMENSION(NDIM2D,TRAFO_MAXNVE), INTENT(IN) :: Dcoords
  
  !</input>
  
  !<output>
  
  ! Arrays with auxiliary jacobian factors for later computation
  REAL(DP), DIMENSION(TRAFO_NAUXJAC), INTENT(OUT) :: DjacPrep
  
  !</output>

!</subroutine>

  DjacPrep(1) = 0.5_DP * (-Dcoords(1,1) - Dcoords(1,2) + Dcoords(1,3) + Dcoords(1,4))
  DjacPrep(2) = 0.5_DP * ( Dcoords(1,1) - Dcoords(1,2) + Dcoords(1,3) - Dcoords(1,4))
  DjacPrep(3) = 0.5_DP * (-Dcoords(2,1) + Dcoords(2,2) - Dcoords(2,3) + Dcoords(2,4))
  DjacPrep(4) = 0.5_DP * (-Dcoords(2,1) + Dcoords(2,2) + Dcoords(2,3) - Dcoords(2,4))

  END SUBROUTINE 

  !************************************************************************

  PURE SUBROUTINE trafo_calcTrafo2 (Dcoord,DjacPrep,Djac,ddetj, &
                                   dparx,dpary,dxreal,dyreal)

  !<description>

  ! This subroutine performs two tasks:
  ! -Initialisation of a a given 2x2 matrix with the
  !  mapping information from the reference element to the "real"
  !  quadrilateral. Calculation of the Jacobian determinant
  ! -Transformation of a given point on the reference element onto
  !  the "real" element
  ! Both things are performed simultaneously because the jacobian
  ! determinant is dependent of the point.
  !
  ! Before this routine can be called, the auxiliary factors DjacPrep
  ! have to be calculated with elem_calcJacPrepare for the 
  ! considered element.
  
  !</description>  

  !<input>
  
  ! Coordinates of the fout points forming the element.
  ! Dcoord(1,.) = x-coordinates,
  ! Dcoord(2,.) = y-coordinates.
  REAL(DP), DIMENSION(NDIM2D,TRAFO_MAXNVE), INTENT(IN) :: Dcoord (2,4)
  
  ! Auxiliary constants for the considered element with
  ! coordinates in Dcoord; have to be computed previously
  ! by elem_calcJacPrepare.
  REAL(DP), DIMENSION(TRAFO_NAUXJAC), INTENT(IN) :: DjacPrep
  
  ! Coordinates of a point on the reference element
  REAL(DP), INTENT(IN) :: dparx,dpary
  
  !</input>
  
  !<output>
  
  ! The Jacobian matrix of the mapping from the reference to the 
  ! real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  REAL(DP), DIMENSION(TRAFO_NJACENTRIES), INTENT(OUT) :: Djac
  
  ! The determinant of the mapping.
  REAL(DP), INTENT(OUT) :: ddetj
  
  ! X-Coordinates of a point on the real element
  REAL(DP), INTENT(OUT) :: dxreal
  
  ! X-Coordinates of a point on the real element
  REAL(DP), INTENT(OUT) :: dyreal

  !</output>
  
!</subroutine>

  Djac(1) = 0.5E0_DP * ((Dcoord(1,2)-Dcoord(1,1)+DjacPrep(2)) &
                          + DjacPrep(2)*dpary )
  Djac(2) = 0.5E0_DP * ( DjacPrep(4) - DjacPrep(3)*dpary )
  Djac(3) = 0.5E0_DP * ( DjacPrep(1) + DjacPrep(2)*dparx )
  Djac(4) = 0.5E0_DP * ((Dcoord(2,3)-Dcoord(2,1)-DjacPrep(4)) &
                          - DjacPrep(3)*dparx )

  ! Determinant of the mapping

  ddetj = Djac(1)*Djac(4) - Djac(3)*Djac(2)
  
  ! Map the point to the real element
  
  dxreal = 0.5E0_DP*((Dcoord(1,1)+Dcoord(1,2)+DjacPrep(1)) + &
                     (Dcoord(1,2)-Dcoord(1,1)+DjacPrep(2))*dparx + &
                     DjacPrep(1)*dpary + &
                     DjacPrep(2)*dparx*dpary )
  dyreal = 0.5E0_DP*((Dcoord(2,1)+Dcoord(2,3)+DjacPrep(3)) + &
                     DjacPrep(4)*dparx + &
                     (Dcoord(2,3)-Dcoord(2,1)-DjacPrep(4))*dpary - &
                     DjacPrep(3)*dparx*dpary )
  
  END SUBROUTINE

  !************************************************************************

  PURE SUBROUTINE trafo_calcJac2 (Dcoord,DjacPrep,Djac,ddetj,dparx,dpary)

  !<description>

  ! Calculate Jacobian determinant of mapping from reference- to
  ! real element.
  !
  ! This routine only calculates the Jacobian determinant of the
  ! mapping. This is on contrast to elem_calcRealCoords, which not only 
  ! calculates this determinant but also maps the point. So this routine
  ! can be used to speed up the code if the coordinates of the
  ! mapped point already exist.
  !
  ! Before this routine can be called, the auxiliary factors DJF
  ! have to be calculated with elem_calcJacPrepare for the considered element.
  
  !</description>  

  !<input>
  
  ! Coordinates of the fout points forming the element.
  ! Dcoord(1,.) = x-coordinates,
  ! Dcoord(2,.) = y-coordinates.
  REAL(DP), DIMENSION(NDIM2D,TRAFO_MAXNVE), INTENT(IN) :: Dcoord (2,4)
  
  ! Auxiliary constants for the considered element with
  ! coordinates in Dcoord; have to be computed previously
  ! by elem_calcJacPrepare.
  REAL(DP), DIMENSION(TRAFO_NAUXJAC), INTENT(IN) :: DjacPrep
  
  ! Coordinates of a point on the reference element
  REAL(DP), INTENT(IN) :: dparx,dpary
  
  !</input>
  
  !<output>
  
  ! The Jacobian matrix of the mapping from the reference to the 
  ! real element.
  !  Djac(1) = J(1,1)
  !  Djac(2) = J(2,1)
  !  Djac(3) = J(1,2)
  !  Djac(4) = J(2,2)
  REAL(DP), DIMENSION(TRAFO_NJACENTRIES), INTENT(OUT) :: Djac
  
  ! The determinant of the mapping.
  REAL(DP), INTENT(OUT) :: ddetj
  
  !</output>
  
!</subroutine>

  ! Jacobian matrix of the mapping:

  Djac(1) = 0.5E0_DP * (Dcoord(1,2)-Dcoord(1,1)+DjacPrep(2)) &
            + 0.5E0_DP * DjacPrep(2)*dpary
  Djac(2) = 0.5E0_DP * DjacPrep(4) - 0.5E0_DP * DjacPrep(3)*dpary
  Djac(3) = 0.5E0_DP * DjacPrep(1) + 0.5E0_DP * DjacPrep(2)*dparx
  Djac(4) = 0.5E0_DP * (Dcoord(2,3)-Dcoord(2,1)-DjacPrep(4)) &
            - 0.5E0_DP * DjacPrep(3)*dparx

  ! Determinant of the mapping

  ddetj = Djac(1)*Djac(4) - Djac(3)*Djac(2)
  
  END SUBROUTINE

  !************************************************************************

  PURE SUBROUTINE trafo_calcRealCoords2 (Dcoord,DjacPrep,dparx,dpary,dxreal,dyreal)

  !<description>

  ! This subroutine computes the real coordinates of a point which 
  ! is given by parameter values (dparx,dpary) on the reference element. 
  ! It is assumed that the array DjacPrep was initialised before using
  ! with elem_calcJacPrepare.
  
  !</description>  

  !<input>
  
  ! Coordinates of the fout points forming the element.
  ! Dcoord(1,.) = x-coordinates,
  ! Dcoord(2,.) = y-coordinates.
  REAL(DP), DIMENSION(NDIM2D,TRAFO_MAXNVE), INTENT(IN) :: Dcoord
  
  ! Auxiliary constants for the considered element with
  ! coordinates in Dcoord; have to be computed previously
  ! by elem_calcJacPrepare.
  REAL(DP), DIMENSION(TRAFO_NAUXJAC), INTENT(IN) :: DjacPrep
  
  ! Coordinates of a point on the reference element
  REAL(DP), INTENT(IN) :: dparx,dpary
  
  !</input>
  
  !<output>
  
  ! X-Coordinates of a point on the real element
  REAL(DP), INTENT(OUT) :: dxreal
  
  ! X-Coordinates of a point on the real element
  REAL(DP), INTENT(OUT) :: dyreal

  !</output>
  
!</subroutine>

  ! Map the point to the real element
  
  dxreal = 0.5E0_DP*((Dcoord(1,1)+Dcoord(1,2)+DjacPrep(1)) + &
                     (Dcoord(1,2)-Dcoord(1,1)+DjacPrep(2))*dparx + &
                     DjacPrep(1)*dpary + &
                     DjacPrep(2)*dparx*dpary )
  dyreal = 0.5E0_DP*((Dcoord(2,1)+Dcoord(2,3)+DjacPrep(3)) + &
                     DjacPrep(4)*dparx + &
                     (Dcoord(2,3)-Dcoord(2,1)-DjacPrep(4))*dpary - &
                     DjacPrep(3)*dparx*dpary )
  
  END SUBROUTINE

END MODULE 
