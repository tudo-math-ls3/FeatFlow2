**********************************************************************
* This file contains functions for volume- and boundary-integration.
* These routines implement a general way to integrate about
* boundries/volumes. The caller always have to specify some
* callback-routines which collect the information.
**********************************************************************

**********************************************************************
* General volume integration.
*
* Performs a loop over all elements to perform a volume integration.
* Calls user provided callback routines to calculate the integral
* value.
*
* In:
*  DU1, DU2, DP - solution vector with velocity and pressure
*  TRIA         - array [1..SZTRIA] of integer
*                 Triangulation structure of the underlying mesg
*  KVERT, KMID  - Arrays with information about the triangulation;
*                 must correspond to TRIA!
*  DCORVG       - Coordinates of vertices
*  ELE          - Element that's used for calculating the 
*                 continuous velocity
*  BNONPR       - Type of element: 
*                  false = parametric element (E0xx)
*                  true  = nonparametric element (EMxx)
*  ICUB         - Cubature formula to use for vol. integration
*
*  CINIT        - Initialisation routine for calculation
*  CVDER        - Calculate values/derivatives in cubature point
*  CINSUM       - Integrate and sum up calculated values
*  CDONE        - Postprocessing routine for calculation
*
*  IINFO        - Integer array with information for calculation
*                 routines; structure is administrated by the caller
*  NIINFO       - Length of IINFO
*  DINFO        - Double array with information for calculation
*                 routines; structure is administrated by the caller
*  NDINFO       - Length of DINFO
*  BUDEFC       - User defined calculation only.
*                 =false: normal calculation
*                 =true: omit standard calculation of information
*                  about velocity and pressure; *only* perform
*                  user defined calculation by Cxxxx-routines.
*                  In this case DU1, DU2, DP is allowed to be
*                  undefined (more exactly: can point to dummy
*                  variables) if not necessary for the calculation.
*  IGEOM  - array [1..*] of integer 
*  DGEOM  - array [1..*] of double 
*           Integer- and double-precision parameter blocks with
*           geometry information. Passed to boundary
*           routines. Not used in this routine.
*  
* Out:
*  IINFO,
*  DINFO        - As calculated by the calculation routines
*  
* Functions the caller has to provide:
*
* SUBROUTINE CINIT (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,
*                   TRIA,ELE,BNONPR)
* --------------------------------------------------------------------
* IMPLICIT NONE
* INTEGER NIINFO,IINFO(NIINFO),NDINFO,TRIA(*)
* LOGICAL BNONPR
* DOUBLE PRECISION DINFO(NDINFO),DU1(*),DU2(*),DP(*)
* EXTERNAL ELE
*  -> Initialises the calculation routines, e.g. by setting all values
*     in IINFO/DINFO to 0 and initialising necessary structures there
*  -> The cubature formula is already initalised before the call to
*     CINIT.
*  -> If a conforming element is used, the element is already
*     initialised, but no values are calculated before.
* 
* SUBROUTINE CDONE (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,
*                   TRIA,ELE,BNONPR)
* --------------------------------------------------------------------
* IMPLICIT NONE
* INTEGER NIINFO,IINFO(NIINFO),NDINFO,TRIA(*)
* LOGICAL BNONPR
* DOUBLE PRECISION DINFO(NDINFO),DU1(*),DU2(*),DP(*)
* EXTERNAL ELE
*  -> Can perform some postprocessing with the IINFO/DINFO values
*  -> Has to perform cleanup of everything that was prepared in CINIT
*
* The following two subroutines are called for every quadrature point
* to calculate the values of interest and to sum these values up:
*
* SUBROUTINE CVDERE (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,
*                    IELEM,TRIA,ELE,BNONPR,IDFL,KDFG,KDFL,IGEOM,DGEOM)
* --------------------------------------------------------------------
* INTEGER NIINFO,IINFO(NIINFO),NDINFO,IELEM,TRIA(*),IGEOM(*)
* LOGICAL BNONPR
* DOUBLE PRECISION DINFO(NDINFO),DU1(*),DU2(*),DP(*),DGEOM(*)
* INTEGER IDFL,KDFG(*),KDFL(*)
* EXTERNAL ELE
*  -> Is called once for each element before the loop through the
*     cubature points
*  -> Can make any preparations for the current element IELEM
*  -> The element ELE is already prepared for obtaining the values
*     in the cubature points. The local and global degrees of
*     freedom are already initialised.
*  -> IELEM = number of current element
*  -> IDFL,KDFG,KDFL contain information about number of DOF's per
*     element as well as local and global DOF's in element IELEM
*  -> The function has to store any values of interest in IINFO/DINFO
*     for later (user defined) use.
*
* SUBROUTINE CVDERP (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,
*                    XX,YY,XX1,YY1,ICUBP,IELEM,TRIA,ELE,BNONPR,
*                    IDFL,KDFG,KDFL)
* --------------------------------------------------------------------
* IMPLICIT NONE
* INTEGER NIINFO,IINFO(NIINFO),NDINFO,IELEM,ICUBP,TRIA(*)
* LOGICAL BNONPR
* DOUBLE PRECISION DINFO(NDINFO),DU1(*),DU2(*),DP(*)
* DOUBLE PRECISION XX,YY,XX1,YY1
* INTEGER IDFL,KDFG(*),KDFL(*)
* EXTERNAL ELE
*  -> Is called once for each cubature point
*  -> Calculates user specified values and derivatives in the current
*     cubature point (XX,YY) on the real or (XX1,YY1) on the
*     reference element, respectively. 
*  -> The points (XX,YY),(XX1,YY1) and the element ELE are only
*     provided for completeness. The variables DBAS(KDFL(I),J)
*     are already initialised with the values of the local basis
*     functions in this cubature point.
*  -> ICUBP = number of current cubature point on 
*  -> IELEM = number of current element
*  -> IDFL,KDFG,KDFL contain information about number of DOF's per
*     element as well as local and global DOF's in element IELEM
*  -> The function has to store these values in IINFO/DINFO
*     for later use in CINSUM.
* Remark: The function values + 1st derivatives of the velocity and
*  the value of the pressure on the current element are
*  automatically calculated. This routine is for the calculation of 
*  user specified values.
*
* SUBROUTINE CINSUM (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,
*                    DETJ,ICUBP,IELEM,OMEGA,
*                    DU1V,DU1X,DU1Y,DU2V,DU2X,DU2Y,DPV)
* --------------------------------------------------------------------
* IMPLICIT NONE
* INTEGER NIINFO,IINFO(NIINFO),NDINFO,ICUBP,IELEM
* DOUBLE PRECISION DINFO(NDINFO),DU1(*),DU2(*),DP(*)
* DOUBLE PRECISION DETJ,OMEGA
* DOUBLE PRECISION DU1V,DU1X,DU1Y,DU2V,DU2X,DU2Y,DPV
* EXTERNAL ELE
*  -> Forms the integrand and sums up the calculated values
*     in the current cubature point. This is the main calculation
*     routine for the values to be calculated. It has to store 
*     the values of interest in the used-defined arrays IINFO/DINFO.
*  -> Is called directly after CVDER, with the values of velocity/
*     pressure and the cubature formula weight as a parameter.
*  -> DU1,DU2,DP are the vectors representing velocity and pressure.
*  -> DETJ = value of the Jacobian determinant of the mapping
*     to the reference element.
*  -> IELEM = number of current element
*  -> ICUB = number of current cubature point on 
*  -> OMEGA = Weight of current cubature point in cubature formula
*  -> DU1V,DU1X,DU1Y,DU2V,DU2X,DU2Y,DPV
*     These represent function value, value of X- and value of
*     Y-derivative for X- and Y-comp. of the velocity in the
*     current cubature point. DPV is the value of the pressure there.
*     Other necessary values have to be taken from IINFO/DINFO which
*     had been filled by CVDER before.
*     In case of a pure user defined calculation (BUDEFC=TRUE)
*     these variables are all set to 0D0.
*
* All parameters are calculated in linint and passed to these
* functions. No function has an output parameter to pass back to
* VOLINT. The functions can use the parameters for their need but are
* not allowed to modify anything except for the user defined arrays
* DINFO/IINFO!
* After the calculation the caller of VOLINT can find the values of
* interest in the user defined arrays DINFO/IINFO as calculated
* in the provided function CINSUM.
*
* The routine provides information about the current element in the
* variable /ELEM/.IEL.
**********************************************************************

      SUBROUTINE VOLINT (DU1,DU2,DP,TRIA,KVERT,KMID,
     *    DCORVG,ELE,BNONPR,ICUB,
     *    IINFO, NIINFO, DINFO, NDINFO,
     *    CINIT,CVDERE,CVDERP,CINSUM,CDONE,BUDEFC,
     *    IGEOM,DGEOM)

      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasictria.inc'

      INCLUDE 'ccub.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      
      INCLUDE 'stria.inc'
      
C parameters

      DOUBLE PRECISION DU1(*),DU2(*),DP(*),DCORVG(2,*),DGEOM(*)
      INTEGER KVERT(NNVE,*),KMID(NNVE,*),IGEOM(*)
      INTEGER ICUB
      LOGICAL BNONPR
      
      INTEGER NIINFO,IINFO(NIINFO),NDINFO,TRIA(SZTRIA)
      DOUBLE PRECISION DINFO(NDINFO)
      EXTERNAL CINIT,CVDERE,CVDERP,CINSUM,CDONE
      LOGICAL BUDEFC

C externals

      EXTERNAL ELE
      INTEGER NDFL
      EXTERNAL NDFL
      
C local variables

      INTEGER IELTYP,I,IVE,JP,IG
      DOUBLE PRECISION DU1V,DU1X,DU1Y,DU2V,DU2X,DU2Y,DPV
      DOUBLE PRECISION DJF(2,2),DCOORD(2,4),OM
      DOUBLE PRECISION XX,YY,XI1,XI2
      INTEGER IDFL,KDFG(NNCUBP),KDFL(NNCUBP)

c     *** preparation - evaluation of parameters
      IER=0

c     *** which derivatives of basis functions are needed?
      DO  I = 1,NNDER
        BDER(I)=.FALSE.
      ENDDO
      
      BDER(1)=.TRUE.
      BDER(2)=.TRUE.
      BDER(3)=.TRUE.
      
c     *** dummy call of ele sets number of element

      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)

      IDFL=NDFL(IELTYP)

C initialise cubature formula

      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
      
C Prepare parametric element for cubature. Save the current cubature
C formula into the COMMON block variable ICUBP - either for now or
C for later use. Perform a dummy call to the element in the conforming
C case.
C This is done for saving arithmetic operations in later calls.
C
C In the nonconforming case the element has to be initialised
C separately on every element, so there's no use in initialising it
C globally.

      ICUBP=ICUB
      
      IF (.NOT.BNONPR) THEN
        CALL ELE(0D0,0D0,-2)
      END IF

C User defined initialisation

      CALL CINIT (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,TRIA,ELE,BNONPR)

C loop over all elements
C
C The variable IEL is stored in the COMMON block /ELEM/ in case that
C any procedure has to know in which element we are...

      DO IEL=1,TRIA(ONEL)
      
C Get the degrees of freedom:

        CALL NDFGLX(TRIA,IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
        IF (IER.LT.0) GOTO 99999

C Store the coordinates of the corners of the current element in 
C the COMMON block variables DX/DY (necessary for the element) 
C as well as in the DCOORD array for the coordinate transformation.

        DO IVE = 1, TRIA(ONVE)
          JP=KVERT(IVE,IEL)
          KVE(IVE)=JP
          DX(IVE)=DCORVG(1,JP)
          DY(IVE)=DCORVG(2,JP)
          DCOORD (1,IVE) = DCORVG(1,JP)
          DCOORD (2,IVE) = DCORVG(2,JP)
        END DO

C Prepare nonparametric elements for cubature on current element.
C The COMMON block variable ICUBP was initialised earlier...
C
C Because the later loop will change ICUBP as it is used for looping 
C through the cubature points on the element, so we have to reset it 
C here.

        ICUBP = ICUB
        
        IF (BNONPR) THEN
          CALL ELE(0D0,0D0,-2)
        END IF

C Initialise auxiliary Jacobian factors DJF for transformation

        CALL QINIJF (DCOORD,DJF)
       
C Perform user defined element-calculations:
       
        CALL CVDERE (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,
     *               IEL,TRIA,ELE,BNONPR,IDFL,KDFG,KDFL,IGEOM,DGEOM)

c Loop over all cubature points

        DO ICUBP = 1, NCUBP
          
c Cubature point on the reference element

          XI1=DXI(ICUBP,1)
          XI2=DXI(ICUBP,2)
          
C Calculate Jacobian matrix, determinant and mapping of the cubature
C point (XI1,YI1) on the reference element into (XX,YY) on the "real"
C element. The Jacobian matrix is stored in the COMMON block variable 
C DJAC - necessary for the element.

          CALL  QTRAF (DCOORD,DJF,DJAC,DETJ,XI1,XI2,XX,YY)

C Calculate the weighting factor for the current cubature point
C with the help of the Jacobian determinant

          OM = DOMEGA(ICUBP)*DETJ
          
C Evaluate the basis functions in the cubature point
C for the velocities

          IF(BNONPR) THEN
            CALL ELE(XX,YY,-3)
          ELSE
            CALL ELE(XI1,XI2,-3)
          ENDIF
          IF (IER.LT.0) GOTO 99999

C Perform standard calculations

          IF (.NOT.BUDEFC) THEN

C Evaluate the solution values and derivatives in the cubature point:

C X-Component:

            DU1V=0D0 ! value
            DU1X=0D0 ! x dreiv.
            DU1Y=0D0 ! y deriv
            DO I=1,IDFL
              IG=KDFG(I)
              DU1V=DU1V+DU1(IG)*DBAS(KDFL(I),1)
              DU1X=DU1X+DU1(IG)*DBAS(KDFL(I),2)
              DU1Y=DU1Y+DU1(IG)*DBAS(KDFL(I),3)
            END DO
          
C Y-Component:
          
            DU2V=0D0 ! value
            DU2X=0D0 ! x dreiv.
            DU2Y=0D0 ! y deriv
            DO I=1,IDFL
              IG=KDFG(I)
              DU2V=DU2V+DU2(IG)*DBAS(KDFL(I),1)
              DU2X=DU2X+DU2(IG)*DBAS(KDFL(I),2)
              DU2Y=DU2Y+DU2(IG)*DBAS(KDFL(I),3)
            END DO
          
C Pressure:
          
            DPV=DP(IEL)
            
          ELSE
          
C User defined calculation - don't do anything
          
            DU1V = 0D0
            DU1X = 0D0
            DU1Y = 0D0
            DU2V = 0D0
            DU2X = 0D0
            DU2Y = 0D0
            DPV  = 0D0
            
          END IF
          
C Calculate user defined values and derivatives:
          
          CALL CVDERP (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,
     *                 XX,YY,XI1,XI2,ICUBP,IEL,TRIA,ELE,BNONPR,
     *                 IDFL,KDFG,KDFL)

C Form the integrand and sum up values of interest in the
C user defined information arrays:

          CALL CINSUM (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,
     *                 DETJ,ICUB,IEL,OM,
     *                 DU1V,DU1X,DU1Y,DU2V,DU2X,DU2Y,DPV)
        
        END DO
        
      END DO

C Do some postprocessing with the calculated values and
C the cleanup if necessary

      CALL CDONE (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,TRIA,ELE,BNONPR)

99999 END
      
      
**********************************************************************
* General line integration.
*
* Performs a loop over all elements to perform a line integration
* (i.e. on the reconstructed interface of fictitious boundary
* components).
* Calls user provided callback routines to calculate the integral
* value.
*
* In:
*  DU1, DU2, DP - solution vector with velocity and pressure
*  TRIA         - array [1..SZTRIA] of integer
*                 Triangulation information about the underlying mesh.
*  KVERT, KMID, 
*  DCORVG,
*  DCORMG       - Usual geometry information; must correspont to TRIA!
*  ELE          - Element that's used for calculating the 
*                 continuous velocity
*  BNONPR       - Type of element: 
*                  false = parametric element (E0xx)
*                  true  = nonparametric element (EMxx)
*  ICUB         - Cubature formula to use for vol. integration,
*                 see CB2QL.
*
*  CINIT        - Initialisation routine for calculation
*  CVDER        - Calculate values/derivatives in cubature point
*  CINSUM       - Integrate and sum up calculated values
*  CDONE        - Postprocessing routine for calculation
*  CGETLN       - Get information about the interface where the
*                 integration takes part
*
*  IINFO        - Integer array with information for calculation
*                 routines; structure is administrated by the caller
*  NIINFO       - Length of IINFO
*  DINFO        - Double array with information for calculation
*                 routines; structure is administrated by the caller
*  NDINFO       - Length of DINFO
*  BUDEFC       - User defined calculation only.
*                 =false: normal calculation
*                 =true: omit standard calculation of information
*                  about velocity and pressure; *only* perform
*                  user defined calculation by Cxxxx-routines.
*                  In this case DU1, DU2, DP is allowed to be
*                  undefined (more exactly: can point to dummy
*                  variables) if not necessary for the calculation.
*  IGEOM  - array [1..*] of integer 
*  DGEOM  - array [1..*] of double 
*           Integer- and double-precision parameter blocks with
*           geometry information. Passed to boundary
*           routines. Not used in this routine.
*  IINFO  - array [1..*] if integer
*  DINFO  - array [1..*] of double
*           Integer- and double precision data block. Passed to
*           boundary routines. Not used in this routine.
*  
* Out:
*  IINFO,
*  DINFO        - As calculated by the calculation routines
*  
* Functions the caller has to provide:
*
* SUBROUTINE CINIT (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,
*                   TRIA,ELE,BNONPR)
* --------------------------------------------------------------------
* see above
* 
* SUBROUTINE CDONE (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,
*                   TRIA,ELE,BNONPR)
* --------------------------------------------------------------------
* see above
*
* SUBROUTINE CGETLN (IINFO,NIINFO,DINFO,NDINFO,IELEM,TRIA,
**                   DCORVG,KVERT,DCORMG,KMID,                 
**                   X1,Y1,X2,Y2,BDOINT,IOGEOM,DGEOM)
* --------------------------------------------------------------------
* IMPLICIT NONE
* INTEGER NIINFO,IINFO(NIINFO),NDINFO,IELEM,TRIA(*),IGEOM(*)
* DOUBLE PRECISION DINFO(NDINFO),X1,Y1,X2,Y2
* DOUBLE PRECISION DCORVG(2,*),DCORMG(2,*),DGEOM(*)
* INTEGER KVERT(NNVE,*),KMID(NNVE,*),NEL,NVE
* LOGICAL BDOINT
* -> Is called for every element.
* -> If there's nothing to integrate on element IELT, BDOINT must be
*    set to .FALSE.
* -> If the integration on element IELEM should be performed, BDOINT
*    must be set to .TRUE. In this case the variables (X1,Y1) and
*    (X2,Y2) must be initialised with the coordinates of starting- and
*    ending point of the line where to integrate on element IELT.
*    The coordinates have to be given in "real" coordinates, not on the
*    reference element.
*    The points do not necessarily have to be positioned on the
*    boundary edges of the quadrilateral.
* -> If the any summation routine below needs information about the
*    reconstructed interface, this routine should save the starting-
*    and ending-point of the line in the DINFO array!
*
* SUBROUTINE CVDERE (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,
*                    IELEM,TRIA,ELE,BNONPR,IDFL,KDFG,KDFL)
* --------------------------------------------------------------------
* see above
*
* SUBROUTINE CVDERP (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,
*                    XX,YY,XX1,YY1,ICUBP,IELEM,TRIA,ELE,BNONPR,
*                    IDFL,KDFG,KDFL)
* --------------------------------------------------------------------
* see above. 
*
* SUBROUTINE CINSUM (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,
*                    DETJ,ICUBP,IELEM,OMEGA,
*                    DU1V,DU1X,DU1Y,DU2V,DU2X,DU2Y,DPV)
* --------------------------------------------------------------------
* see above
*
* SUBROUTINE CBPROJ (IINFO,NIINFO,DINFO,NDINFO,
**                   TRIA,DCORVG,KVERT,DCORMG,KMID,   
**                   DCOORD,IELEM,X1,Y1,X2,Y2,DCUBP,DJAC)
* --------------------------------------------------------------------
* IMPLICIT NONE
* INTEGER NIINFO,IINFO(NIINFO),NDINFO,IELEM,TRIA(*)
* INTEGER KVERT(NNVE,*),KMID(NNVE,*),NEL,NVE
* DOUBLE PRECISION DINFO(NDINFO),DCORVG(2,*),DCORMG(2,*)
* DOUBLE PRECISION X1,Y1,X2,Y2,DJAC
* DOUBLE PRECISION DCOORD(2,4),DCUBP(2,NNCUBP)
*
* -> Cubature projection subroutine. Must project the parameter values
*    on the reference interval [-1,1] into coordinates on the element 
*    IEL whose corner vertices are given by the coordinates DCOORD.
* -> IEL may be 0 if there's no element belonging to the coordinates
*    in DCOORD.
* -> The caller should provide a CBPROJ routine which calls the
*    fictitious boundary subroutine FBCBPR with the correct
*    number of the fictitious boundary component. The number
*    of this boundary component can be obtained during the previously
*    executed call to CGETLN.
* -> The routine CGETLN is always called before this routine, i.e.
*    one call of CGETLN can follow multiple calls to CBPROJ
*    regarding the same line segment. This way the caller can
*    save some information about the current line (perhaps number
*    of the boundary component it belongs to) and then pass this
*    to the projection method.
* -> For standard distribution of the cubature points along the
*    line [(X1,Y1),(X2,Y2)] XCBLPR can be used here. This wrapper
*    routine calls CBLPRQ directly, which performs a distribution
*    of the cubature points along a line.
*    For the exact treatment of fictitious boundary interfaces
*    CBPROJ should call the user specified fictitious boundary
*    handling routine FBCBPR.
* -> For a detailed description of the parameters see CBLPRQ.
*
* All parameters are calculated in LININT and passed to these
* functions. No function has an output parameter to pass back to
* LININT. The functions can use the parameters for their need but are
* not allowed to modify anything except for the user defined arrays
* DINFO/IINFO!
* After the calculation the caller of LININT can find the values of
* interest in the user defined arrays DINFO/IINFO as calculated
* in the provided function CINSUM.
*
* The routine provides information about the current element in the
* variable /ELEM/.IEL.
*
* The routine only supports finite elements that do not depend on
* the quadrature formula ICUB, since the cubature points on the 
* line can be arbitrary in the element!
**********************************************************************

      SUBROUTINE LININT (DU1,DU2,DP,TRIA,DCORVG,KVERT,DCORMG,KMID,
     *    ELE,BNONPR,ICUB,IINFO, NIINFO, DINFO, NDINFO,
     *    CINIT,CVDERE,CVDERP,CINSUM,CDONE,CGETLN,CBPROJ,BUDEFC,
     *    IGEOM,DGEOM)

      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasictria.inc'

      INCLUDE 'ccub.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      
      INCLUDE 'stria.inc'
      
C parameters

      DOUBLE PRECISION DU1(*),DU2(*),DP(*),DCORVG(2,*),DCORMG(2,*)
      INTEGER KVERT(NNVE,*),KMID(NNVE,*)
      INTEGER ICUB
      LOGICAL BNONPR
      
      DOUBLE PRECISION DGEOM(*)
      INTEGER IGEOM(*)
      
      INTEGER NIINFO,IINFO(NIINFO),NDINFO,TRIA(SZTRIA)
      DOUBLE PRECISION DINFO(NDINFO)
      EXTERNAL CINIT,CVDERE,CVDERP,CINSUM,CDONE,CGETLN,CBPROJ
      LOGICAL BUDEFC

C externals

      EXTERNAL ELE
      INTEGER NDFL
      EXTERNAL NDFL
      
C local variables

      INTEGER IELTYP,I,IVE,JP,IG
      DOUBLE PRECISION DU1V,DU1X,DU1Y,DU2V,DU2X,DU2Y,DPV
      DOUBLE PRECISION DCOORD(2,4),OM
      DOUBLE PRECISION DCUBP(2,NNCUBP)
      DOUBLE PRECISION XX,YY,XI1,XI2
      DOUBLE PRECISION X1,Y1,X2,Y2
      LOGICAL BDOINT
      INTEGER IDFL,KDFG(NNCUBP),KDFL(NNCUBP)
      
c     *** preparation - evaluation of parameters
      IER=0

c     *** which derivatives of basis functions are needed?
      DO  I = 1,NNDER
        BDER(I)=.FALSE.
      ENDDO
      
      BDER(1)=.TRUE.
      BDER(2)=.TRUE.
      BDER(3)=.TRUE.
      
c     *** dummy call of ele sets number of element

      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)

      IDFL=NDFL(IELTYP)

C Prepare conforming element for cubature. Perform a dummy call 
C to the element in the conforming case.
C This is done for saving arithmetic operations in later calls.
C
C In the nonconforming case the element has to be initialised
C separately on every element, so there's no use in initialising it
C globally.
C
C Normally we have to set ICUBP before calling the element.
C But because our later line cubature formula uses quadrature points
C that are more or less arbitrary placed in the element,
C we only support elements here not depending on the type
C of the cubature formula!
C Therefore we don't set ICUBP here before we call the element.
C Elements that depend on this information will either give
C false values or halt the program!

      ICUBP=0
      
      IF (.NOT.BNONPR) THEN
        CALL ELE(0D0,0D0,-2)
      END IF

C User defined initialisation

      CALL CINIT (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,TRIA,ELE,BNONPR)

C Loop over all elements.
C
C The variable IEL is stored in the COMMON block /ELEM/ in case that
C any procedure has to know in which element we are...

      DO IEL=1,TRIA(ONEL)
      
C Get the degrees of freedom:

        CALL NDFGLX(TRIA,IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
        IF (IER.LT.0) GOTO 99999

C Store the coordinates of the corners of the current element in 
C the COMMON block variables DX/DY (necessary for the element) 
C as well as in the DCOORD array for the coordinate transformation.

        DO IVE = 1, TRIA(ONVE)
          JP=KVERT(IVE,IEL)
          KVE(IVE)=JP
          DX(IVE)=DCORVG(1,JP)
          DY(IVE)=DCORVG(2,JP)
          DCOORD (1,IVE) = DCORVG(1,JP)
          DCOORD (2,IVE) = DCORVG(2,JP)
        END DO
        
C Is there at all a line on the current element?

        BDOINT = .FALSE.
        CALL CGETLN (IINFO,NIINFO,DINFO,NDINFO,IEL,TRIA,
     *               DCORVG,KVERT,DCORMG,KMID,
     *               X1,Y1,X2,Y2,BDOINT)

C If the two points are identical, we have the case that exactly the
C corner of one element is hit. Of course there is nothing to integrate
C then! Btw.: There must not be done anything as the length of the
C line is 0, which would produce NAN when dividing by that!

        IF ((X1.EQ.X2).AND.(Y1.EQ.Y2)) BDOINT=.FALSE.

        IF (BDOINT) THEN

C Initialise the line-cubature formula für the current element.
C This will set up the
C quadrature points on the real element as well as on the reference
C element. The cubature points on the reference element are
C stored in DXI/DYI, which can be found in the COMMON block.
C The coordinates on the cubature points on the real element
C will be stored in DCUBP.
C
C DETJ will receive the length of the line and must be used
C as a weighting factor in addition to the weighting factors
C of the quadrature points.

C Unfortunately we can't use CB2LQ here like
C          CALL CB2LQ (ICUB,DCOORD,IEL,X1,Y1,X2,Y2,DCUBP,DETJ,CBPROJ)
C CBPROJ is not parameter compatible!
C as this routine does not call our callback routines. 
C But as the implementation is short, it can directly be done here:

C The CB1 routine will initialise the weighting factors, the NCUBP
C information and the distribution of the cubature points on the
C reference interval [-1,1]:

          CALL CB1 (ICUB)
      
C Transform the parametrisation into real coordinates.
C This is done by a user specified reconstruction routine. A standard
C implementation of such a routine is the routine XCBLPR below!

          CALL CBPROJ (IINFO,NIINFO,DINFO,NDINFO,
     *             TRIA,DCORVG,KVERT,DCORMG,KMID,
     *             DCOORD,IEL,X1,Y1,X2,Y2,DCUBP,DETJ)

C Transform the real coordinates back into coordinates on 
C the reference element - to support conforming
C FE-approaches. This requires the inverse transformation...
C This of course is not necessary in the case of non-parametric elements,
C but we do this in all cases here for sure. Maybe that another
C user-provided routine depend on this.

          DO I=1,NCUBP
            CALL QBTRAF (DCOORD, DXI(I,1), DXI(I,2),
     *                   DCUBP(1,I), DCUBP(2,I))
          END DO

C Prepare nonparametric elements for cubature on current element.
C The COMMON block variable ICUBP was initialised earlier...
C
C We set ICUBP to 0 for the same reason as explained above.
C The later loop will change ICUBP as it is used for looping through
C the cubature points on the element, so we have to reset it here.

          ICUBP=0
          IF (BNONPR) THEN
            CALL ELE(0D0,0D0,-2)
          END IF

C Perform user defined element-calculations:
       
          CALL CVDERE (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,
     *               IEL,TRIA,ELE,BNONPR,IDFL,KDFG,KDFL)

c Loop over all cubature points

          DO ICUBP = 1, NCUBP
          
c Cubature point on the reference element

            XI1 = DXI(ICUBP,1)
            XI2 = DXI(ICUBP,2)
          
C Cubature point on the real element:

            XX = DCUBP(1,ICUBP)
            YY = DCUBP(2,ICUBP)
            
C Calculate the weighting factor for the current cubature point
C with the help of the Jacobian determinant = (length of the line)/2

            OM = DOMEGA(ICUBP)*DETJ
          
C Evaluate the basis functions in the cubature point
C for the velocities

            IF(BNONPR) THEN
              CALL ELE(XX,YY,-3)
            ELSE
              CALL ELE(XI1,XI2,-3)
            ENDIF
            IF (IER.LT.0) GOTO 99999

C Perform standard calculations

            IF (.NOT.BUDEFC) THEN

C Evaluate the solution values and derivatives in the cubature point:

C X-Component:

              DU1V=0D0 ! value
              DU1X=0D0 ! x dreiv.
              DU1Y=0D0 ! y deriv
              DO I=1,IDFL
                IG=KDFG(I)
                DU1V=DU1V+DU1(IG)*DBAS(KDFL(I),1)
                DU1X=DU1X+DU1(IG)*DBAS(KDFL(I),2)
                DU1Y=DU1Y+DU1(IG)*DBAS(KDFL(I),3)
              END DO
          
C Y-Component:
          
              DU2V=0D0 ! value
              DU2X=0D0 ! x dreiv.
              DU2Y=0D0 ! y deriv
              DO I=1,IDFL
                IG=KDFG(I)
                DU2V=DU2V+DU2(IG)*DBAS(KDFL(I),1)
                DU2X=DU2X+DU2(IG)*DBAS(KDFL(I),2)
                DU2Y=DU2Y+DU2(IG)*DBAS(KDFL(I),3)
              END DO
          
C Pressure:
          
              DPV=DP(IEL)
              
            ELSE
          
C User defined calculation - don't do anything
          
              DU1V = 0D0
              DU1X = 0D0
              DU1Y = 0D0
              DU2V = 0D0
              DU2X = 0D0
              DU2Y = 0D0
              DPV  = 0D0
              
            END IF
          
C Calculate user defined values and derivatives:
          
            CALL CVDERP (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,
     *                 XX,YY,XI1,XI2,ICUBP,IEL,TRIA,ELE,BNONPR,
     *                 IDFL,KDFG,KDFL)

C Form the integrand and sum up values of interest in the
C user defined information arrays:

            CALL CINSUM (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,
     *                 DETJ,ICUB,IEL,OM,
     *                 DU1V,DU1X,DU1Y,DU2V,DU2X,DU2Y,DPV)
        
          END DO
          
        END IF
        
      END DO

C Do some postprocessing with the calculated values and
C the cleanup if necessary

      CALL CDONE (IINFO,NIINFO,DINFO,NDINFO,DU1,DU2,DP,TRIA,ELE,BNONPR)

99999 END
            
            
**********************************************************************
* Line cubature on quadrilateral element
*
* This routine initialises a set of quadrature points for integration
* on a line of a quadrilateral element. DCOORD must provide
* the coordinates of the corners of the element. DCUBP will receive
* the coordinates of the quadrature points on the real element.
* DJAC will receive the Jacobian determinant of the mapping between
* the reference interval [-1,1] to the line given by the start- and
* end-coordinates (X1,Y1) / (X2,Y2). 
*
* Furthermore the cubature point coordinates on the real element 
* are transformed back into coordinates on the reference element,
* stored in
*               (DXI(1..NCUBP,1),DXI(1..NCUBP,2))
*
* In:
*  ICUB   - Cubature formula according to CB1
*           (1 = 1-point Gauss, 2=Trapezoidal rule, 3=2 point Gauss,
*            ...)
*  DCOORD - array [1..2,1..4] of double
*           Coordinates of the four corners of the real quadrilateral.
*           DCOORD(1,.) saves the X-, DCOORD(2,.) the Y-coordinates.
*  IEL    - The number of current element or 0, if there is no current
*           element belonging to the coordinates in DCOORD.
*  (X1,Y1),
*  (X2,Y2) - start- and endpoint if the line in "real" coordinates
*            on the "real" element
*  CBPROJ  - Cubature projection subroutine; see below
* 
* Out:
*  DJAC    - Jacobi determinant of the mapping between [-1,1] and the
*            line in real coordinates.
*  DCUBP   - array [1..2,1..NNCUBP] of double
*            Coordinates of cubature points in real coordinates
*            on the quadrilateral given by DCOORD.
*
* Functions the caller has to provide:
*
* SUBROUTINE CBPROJ (DCOORD,IEL,X1,Y1,X2,Y2,DCUBP,DJAC)
* --------------------------------------------------------------------
* INTEGER I,IEL
* DOUBLE PRECISION X1,Y1,X2,Y2,DJAC
* DOUBLE PRECISION DCOORD(2,4),DCUBP(2,NNCUBP)
*
* -> Cubature projection subroutine. Must project the parameter values
*    on the reference interval [-1,1] into coordinates on the element 
*    IEL whose corner vertices are given by the coordinates DCOORD.
* -> IDUMMY is an integer dummy parameter, but is usually used pro
*    provide a number of a boundary component to the projection
*    routine.
* -> IEL may be 0 if there's no element belonging to the coordinates
*    in DCOORD.
* -> For standard distribution of the cubature points along the
*    line [(X1,Y1),(X2,Y2)] the routine CBLPRQ can be passed here.
* -> For a detailed description of the parameters see CBLPRQ!
**********************************************************************
            
      SUBROUTINE CB2LQ (ICUB,DCOORD,IEL,X1,Y1,X2,Y2,DCUBP,DJAC,CBPROJ)
      
      IMPLICIT NONE
      
      INCLUDE 'ccub.inc'
      
      INTEGER ICUB,IEL,I
      DOUBLE PRECISION X1,Y1,X2,Y2,DJAC
      DOUBLE PRECISION DCOORD(2,4),DCUBP(2,NNCUBP)
      EXTERNAL CBPROJ
      
C The CB1 routine will initialise the weighting factors, the NCUBP
C information and the distribution of the cubature points on the
C reference interval [-1,1]:

      CALL CB1 (ICUB)
      
C Transform the parametrisation into real coordinates.
C This is done by a user specified reconstruction routine. A standard
C implementation of such a routine is the routine CBLPRQ below!

      CALL CBPROJ (DCOORD,IEL,X1,Y1,X2,Y2,DCUBP,DJAC)

C Transform the real coordinates back into coordinates on 
C the reference element - to support conforming
C FE-approaches. This requires the inverse transformation...

      DO I=1,NCUBP
        CALL QBTRAF (DCOORD, DXI(I,1), DXI(I,2), DCUBP(1,I), DCUBP(2,I))
      END DO

      END
      
      
**********************************************************************
* Wrapper routine for treatment of linear reconstruction
* of boundary interfaces.
*
* Calls the routine CBLPRQ which performs a linear approximation of
* the line where the integration should be performed, given by the
* points (X1,Y1) and (X2,Y2). CBLPRQ will then distribute 
* cubature points along this line.
*
* This routine can be used as standard implementation of the
* CBPROJ parameter in the line integration routine.
**********************************************************************

      SUBROUTINE XCBLPR (IINFO,NIINFO,DINFO,NDINFO,
     *                   TRIA,DCORVG,KVERT,DCORMG,KMID,    
     *                   DCOORD,IELEM,X1,Y1,X2,Y2,DCUBP,DJAC)
            
      IMPLICIT NONE

      INCLUDE 'cbasictria.inc'

      INTEGER NIINFO,IINFO(NIINFO),NDINFO,IELEM,TRIA(*)
      INTEGER KVERT(NNVE,*),KMID(NNVE,*)
      DOUBLE PRECISION DINFO(NDINFO),DCORVG(2,*),DCORMG(2,*)
      DOUBLE PRECISION X1,Y1,X2,Y2,DJAC
      DOUBLE PRECISION DCOORD(2,4),DCUBP(2,*)
      
C Call the linear distribution routine directly.

      CALL CBLPRQ (DCOORD,IELEM,X1,Y1,X2,Y2,DCUBP,DJAC)
      
      END
      
      
**********************************************************************
* Cubature point projection method, linear version
*
* This subroutine projects cubature points on the reference interval
* [-1,1] to real cubature point coordinates on a given quadrilateral.
* The cubature points on the reference interval are given
* in the COMMON-block variable DXI(1..NCUBP,1), with /CUB/.NCUBP the
* number of cubature points. The points (X1,Y1) and (X2,Y2) are
* the endpoints of a line that is parametrised by the DXI(.,.) 
* parameter values.
*
* The routine now translates the parameter values in DXI(.,.) into
* real coordinates on the line [(X1,Y1),(X2,Y2)]. The X/Y-coordinates
* of these points are stored in (DCUBP(1,1..NCUBP),DCUBP(2,1..NCUBP)).
*
* DJAC will receive the length of the line / 2, which is the Jacobian 
* determinant of  the mapping between the reference interval [-1,1] to 
* the line given by the start- and end-coordinates (X1,Y1) / (X2,Y2).
* This is because the weighting factors in the quadrature rule are 
* normalised to this interval.
*
* In:
*  DCOORD - array [1..2,1..4] of double
*           Coordinates of the four corners of the real quadrilateral.
*           DCOORD(1,.) saves the X-, DCOORD(2,.) the Y-coordinates.
*  IEL    - The number of current element or 0, if there is no current
*           element belonging to the coordinates in DCOORD.
*  (X1,Y1),
*  (X2,Y2) - start- and endpoint if the line in "real" coordinates
*            on the "real" element
* 
* Out:
*  DJAC    - Jacobi determinant of the mapping between [-1,1] and the
*            line in real coordinates.
*  DCUBP   - array [1..2,1..NCUBP] of double
*            Coordinates of cubature points in real coordinates
*            on the quadrilateral given by DCOORD.
**********************************************************************
            
      SUBROUTINE CBLPRQ (DCOORD,IEL,X1,Y1,X2,Y2,DCUBP,DJAC)
      
      IMPLICIT NONE
      
      INCLUDE 'ccub.inc'
      
C parameters
      
      INTEGER I,IEL
      DOUBLE PRECISION X1,Y1,X2,Y2,DJAC
      DOUBLE PRECISION DCOORD(2,4),DCUBP(2,NNCUBP)
      DOUBLE PRECISION A1,A2,B1,B2
      
C The variable IEL is actually not used here. It was added for
C compatibility with other - perhaps user defined - routines
C for reconstructing the quadrature points, which might later
C need this information.
C The varaible IDUMMY is also not used here but declared so that
C the parameter list is conform to other routines of the same kind.
      
C Length of the line / 2 = Jacobi determinant of the mapping:

      DJAC = 0.5D0 * DSQRT((X2-X1)**2+(Y2-Y1)**2)

C Now we have to map the coordinates of the cubature points on the
C reference interval onto our real element.
C For t in [-1,1] we use the transformation 
C    s(t) = (a1*t + b1)
C           (a2*t + b2)
C that fulfills:  
C    s(-1)=(x1,y1), s(1)=(x2,y2)
C
C The linear system gives:
C    a1=(x2-x1)/2, b1=(x2+x1)/2, a2=(y2-y1)/2, b2=(y2+y1)/2

      A1=0.5D0*(X2-X1)
      B1=0.5D0*(X2+X1)
      
      A2=0.5D0*(Y2-Y1)
      B2=0.5D0*(Y2+Y1)

      DO I=1,NCUBP

C In DXI(I,1) we can find the position of the quadrature point.
C Calculate the position in "real" coordinates:

        DCUBP(1,I) = A1*DXI(I,1) + B1
        DCUBP(2,I) = A2*DXI(I,1) + B2

      END DO

      END
            