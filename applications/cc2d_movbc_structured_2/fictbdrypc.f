**********************************************************************
* This file contains additional routines for fictitious boundary
* handling that use precalculated information. The routines are
* user-defineable and an extension to the routines in FICTBDRY.F.
*
* A central routine FBDPRC is called by the framework in order to
* precalculate information. The precalculatied information is then
* passed to new, extended, mesh-dependent fictitious boundary
* routines, which can either use them or fall back the standard
* implementation defined in FICTBDRY.F. 
*
* The node-dependent routines for fictitious boundary handling
* provided by this file normally receive a node-variable INODE
* with a node number to handle (determine distance to this node,...).
* The specification of the value INODE is mesh-dependent. The
* standard treatment of this value is:
*   INODE = 1..NVT                  - vertices
*         = NVT+1..NVT+NMT          - edge midpoints
*         = NVT+NMT+1..NVT+NMT+NEL  - element midpoints
* If the geometry is more complex (e.g. more than one node per
* edge or element), the user has to adapt the code according to that!
* The routines here have therefore to be synchronized to the node
* numbering defined/administrated in structria.f.
*
* Remark: The framework normally tries to work with precalculated
*  information to save computational time. If precalculation is not
*  possible (due to complexity e.g.), the routines in FICTBDRY.F are
*  used directly. 
*  The routines are very much linked and dependent to the routines 
*  in FICTBDRY.F; if new routines/calculation methods are introduced 
*  there, the appropriate changes must be done in this file, too!
*  This way the framework can provide optimal computational costs:
*  Either by using precalculated values or by using the standard
*  routines directly.
**********************************************************************

**********************************************************************
* Fictitious boundary mesh precalculation method
*
* This routine is called by the framework in order to precalculate
* information for a given mesh. It accepts the mesh infomation in
* a TRIA structure array, described in STRIA.INC. From this 
* infomation the routine can pre-calculate any information. The
* routine should return two handles to arrays containing the
* precalculated information -- one handle points to an integer-array
* and one handle to a double array. These handles are passed to
* various routines below to obtain information directly according
* to a mesh (distance-functions,...). These functions can then
* decide on their own whether to use the precalculated values
* or to recalculate the values with the direct routines above.
* If the arrays are not used anymore (e.g. because the mesh is
* deleted), the framework will automatically release these arrays.
* 
* In:
*  TRIA   - array [1..SZTRIA] of integer
*           Triangulation structure STRIA.
*  IGEOM  - array [1..*] of integer
*  DGEOM  - array [1..*] of double precision
*           Parameter blocks describing the geometry
*
* Out:
*  LIARR  - integer
*           Handle to integer array with precalculated information.
*           Is set to 0 if not used.
*  LDARR  - integer
*           Handle to double array with precalculated information.
*           Is set to 0 if not used.
*
* The handles LIARR and LDARR are later directly passed to all FBMxxx-
* routines which can use the information stored there. If
* LIARR and LDARR are set to 0 by this routine, this indicates that
* there is no precalculated information. The FBMxxx-routines then have
* to switch back to the standard implementation FBDxxx defined in the
* file FICTBDRY.F!
**********************************************************************

      SUBROUTINE FBDPRC (TRIA,LIARR,LDARR,IGEOM,DGEOM)

      IMPLICIT NONE
      INCLUDE 'cmem.inc'
      INCLUDE 'stria.inc'
      INCLUDE 'sinigeometry.inc'
      
      INTEGER TRIA(SZTRIA),IGEOM(*)
      INTEGER LIARR,LDARR
      DOUBLE PRECISION DGEOM(*)
      
      EXTERNAL GPCWRP, GIPWRP, GNPWRP, GDPWRP, GPPWRP
      INTEGER GIPWRP, GDPWRP, GPPWRP
      
      INTEGER LTMP,NIMEM,NDMEM,I,INDI,INDD
      INTEGER IINFO(2)
      DOUBLE PRECISION DINFO
      
      INTEGER L1,L2,L0
      INTEGER ICIRTP,IGEOCN
      
C     Get rotation, radius,... from the parameter blocks

      ICIRTP = IGEOM(OICIRTP)
      IGEOCN = IGEOM(OIGEOCN)
      
C     standard implementation: no precalculated information

      LIARR = 0
      LDARR = 0
      
C     More complex geometries are precalculated:

      IF (ICIRTP.GE.50) THEN

C       First we have to count how much memory the geometries need.
C       Create an array for saving all the memory usage:

        CALL ZNEW (2*IGEOCN,3,LTMP,'KTMP  ')
        
C       Loop throuh the objects, count the necessary memory and save
C       it in the array.
        
        NIMEM = 0
        NDMEM = 0
        DO I=1,IGEOCN

C         Get the handle of the integer- and double-prec. parameter
C         block from IGEOMS(2,I) and IGEOMS(3,I), resp.:
C
C          L1 = IGEOMS(2,I)
C          L2 = IGEOMS(3,I)

          L0 = IGEOM(OIGEOMS+3*(I-1))
          L1 = IGEOM(OIGEOMS+3*(I-1)+1)
          L2 = IGEOM(OIGEOMS+3*(I-1)+2)
          
C         Check the type of the object; is it an object with 
C         precalculated information?

          IF (L0.EQ.12) THEN
          
C           Yes, there is precomputed information in that object.
C           Ask how much...
          
            CALL GPCCMP (KWORK(L(L1)), 
     *                   DWORK(L(L2)),
     *                   GPCWRP, GIPWRP, GNPWRP, GDPWRP, GPPWRP,
     *                   TRIA, 0, IINFO, DINFO)
     
            NIMEM = NIMEM + IINFO(1)
            NDMEM = NDMEM + IINFO(2)
            KWORK(L(LTMP)+(I-1)) = IINFO(1)
            KWORK(L(LTMP)+IGEOCN+(I-1)) = IINFO(2)
            
          ELSE

C           No precalculated information in here. Set the memory
C           consumption of that object to 0.

            KWORK(L(LTMP)+(I-1)) = 0
            KWORK(L(LTMP)+IGEOCN+(I-1)) = 0
          
          END IF
          
        END DO
        
C       Now we know how much memory we need. Allocate the memory.
C       At least allocate 1 byte to prevent allocating all the memory
C       if the returned memory usage is 0.

        CALL ZNEW (MAX(1,2*IGEOCN+NIMEM),3,LIARR,'KIARR ')
        CALL ZNEW (MAX(1,NDMEM),1,LDARR,'DDARR ')
        
C In the first 2*IGEOCN entries of KIARR we save the starting addresses
C of the precalculated information for every object. The first object
C has:
C   Starting address in int-array    = 2*IGEOCN
C   Starting address in double-array = 0
C The second object:
C   Starting address in int-array    = 2*IGEOCN+mem. usage of obj. 1
C   Starting address in double-array = mem. usage of obj. 1
C and so on.
C The starting addresses/memory usages are organized as:
C  KIARR(1..IGEOCN) = starting addresses of precalculated integers
C  KIARR(IGEOCN+1,2*IGEOCN) = st. adr. of prec. doubles
C The KTMP-array is prepared for this structure above...
C
C Use the information about the memory usage in KTMP to build this
C information array consecutively:

        KWORK(L(LIARR)) = 2*IGEOCN
        KWORK(L(LIARR)+IGEOCN) = 0

        DO I=1,IGEOCN-1
          KWORK(L(LIARR)+I) = KWORK(L(LIARR)+I-1) + KWORK(L(LTMP)+I-1)
          KWORK(L(LIARR)+IGEOCN+I) = KWORK(L(LIARR)+IGEOCN+I-1)
     *                               + KWORK(L(LTMP)+IGEOCN+(I-1))
        END DO
        
C       Now the Temp-array can be released.
        
        CALL ZDISP (0,LTMP,'KTMP  ')

C       We make a second loop over the objects. This time we
C       calculate the information and store it in the arrays:

        DO I=1,IGEOCN
        
C         Calculate the starting addresses of that object in
C         the work-array

          INDI = L(LIARR) + KWORK(L(LIARR)+I-1)
          INDD = L(LDARR) + KWORK(L(LIARR)+IGEOCN+I-1)
          
C         ...and call the calculation routine.
          
C         Get the handle of the integer- and double-prec. parameter
C         block from IGEOMS(2,I) and IGEOMS(3,I), resp.;
C         L0 gets the type of the object:
C
C          L0 = IGEOMS(1,I)
C          L1 = IGEOMS(2,I)
C          L2 = IGEOMS(3,I)

          L0 = IGEOM(OIGEOMS+3*(I-1))
          L1 = IGEOM(OIGEOMS+3*(I-1)+1)
          L2 = IGEOM(OIGEOMS+3*(I-1)+2)
          
C         When precalculated information is used...

          IF (L0.EQ.12) THEN
          
C           Calculate the information to the arrays identified by L1 
C           and L2

            CALL GPCCMP (KWORK(L(L1)), 
     *                   DWORK(L(L2)),
     *                   GPCWRP, GIPWRP, GNPWRP, GDPWRP, GPPWRP,
     *                   TRIA, 1, KWORK(INDI), DWORK(INDD))
     
          END IF
     
        END DO

      END IF
      
      END

************************************************************************
* ISFBDM - Test for being in a fictitious boundary object
*          Precalculated version
*
* This is the node/mesh-dependent version of ISFBDY.
*
* It has to check whether the given node INODE is in an object 
* surrounded by a fictitious boundary.
* The return value must be = 0 if the point is a normal point and
* != 0 if the node is in an object with a fictitious boundary.
*
* In:
*  INODE  - Node number of the node 
*  TRIA   - array [1..SZTRIA] of integer
*           Triangulation structure STRIA.
*  LIARR  - integer; might be 0 depending on FBDPRC
*           Handle to integer array with precalculated information.
*  LDARR  - integer; might be 0 depending on FBDPRC
*           Handle to double array with precalculated information.
*  IFCB   - =0 or NBCT+1..NBCT+NFBDYC
*           =0: test, if (X,Y) is in any fictitious boundary component;
*               the return value is either 0 or the number of any
*               fictitious boundary component containing this point
*           >0: test, if (X,Y) is contained especially in boundary
*               component IFBC. The return value is either 0 or
*               +-IFBC, depending on whether IFBC contains (X,Y) or not.
*           Fictitious boundary components are numbered
*           in the range NBCT+1..NBCT+NFBDYC with numbers behind
*           the normal boundary components.
*   IGEOM - array [1..*] of integer 
*   DGEOM - array [1..*] of double 
*           Integer- and double-precision parameter blocks with
*           geometry information. 
*
* result:
*   0 - If (x,y) is not in a fictitious boundary object
*  >0 - If (x,y) is in a rigid fictitious boundary object.
*  <0 - If (x,y) if in a virtual fictitious boundary object
*
* For a detailed documentation about the meaning of the different
* quantities here, see ISFBDY!
************************************************************************

      INTEGER FUNCTION ISFBDM (INODE, TRIA, LIARR, LDARR, IFBC, 
     *                         IGEOM, DGEOM)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'stria.inc'
      INCLUDE 'sinigeometry.inc'
      INCLUDE 'sgeometries.inc'
      
C     parameters

      INTEGER INODE, LIARR, LDARR, IFBC, IGEOM(*)
      INTEGER TRIA(SZTRIA)
      DOUBLE PRECISION DGEOM(*)
      
C     externals

      INTEGER ISFBDY,GIPCMP,GIPWRP,TNBC,GINWRP
      EXTERNAL ISFBDY,GIPCMP,GIPWRP,TNBC,GINWRP

C     local variables

      DOUBLE PRECISION X,Y
      INTEGER INDI,INDD,I,J,L1,L2,L0
      DOUBLE PRECISION CORSYS(SZCORS)

      INTEGER ICIRTP,IGEOCN
      
C     Get rotation, radius,... from the parameter blocks

      ICIRTP = IGEOM(OICIRTP)
      IGEOCN = IGEOM(OIGEOCN)
      
C     standard return value:
      
      ISFBDM = 0

C     Standard implementation: Determine coordinates and call ISFBDY!

      IF ((LIARR.EQ.0).OR.(LDARR.EQ.0).OR.(ICIRTP.LT.50)) THEN
      
C       We have exactly 1 simple-type geometry particle - circle,
C       square or something like that.
        
        CALL NDE2XY (INODE,TRIA,X,Y)
        ISFBDM = ISFBDY (X,Y,IFBC,IGEOM,DGEOM)

      ELSE 

C       Composed geometry or multiple simple-type geometries,
C       some probably with precomputed information.
C       Here we have to take into account whether to test for a special
C       geometry or generally for all geometries!

        IF (IFBC.NE.0) THEN

C         Calculate the starting addresses of that object in 
C         the work-array

          INDI = L(LIARR) + KWORK(L(LIARR)+IFBC-TNBC()-1)
          INDD = L(LDARR) + KWORK(L(LIARR)+IGEOCN+IFBC-TNBC()-1)
          
C         Get the handle of the integer- and double-prec. parameter
C         block from IGEOMS(2,IFBC-TNBC()) and
C         IGEOMS(3,IFBC-TNBC()), resp.;
C
C          L0 = IGEOMS(1,IFBC-TNBC())
C          L1 = IGEOMS(2,IFBC-TNBC())
C          L2 = IGEOMS(3,IFBC-TNBC())

          L0 = IGEOM(OIGEOMS+3*(IFBC-TNBC()-1))
          L1 = IGEOM(OIGEOMS+3*(IFBC-TNBC()-1)+1)
          L2 = IGEOM(OIGEOMS+3*(IFBC-TNBC()-1)+2)
          
C         Do we have precalculated information?

          IF (L0.EQ.12) THEN
          
C           Yes! use the precalc-routines to determine the return 
C           value.
          
            IF (GIPCMP (KWORK(L(L1)), 
     *                DWORK(L(L2)), 
     *                GIPWRP, TRIA,
     *                KWORK(INDI),DWORK(INDD),
     *                INODE).GT.0) THEN
              ISFBDM = IFBC
            END IF
          
          ELSE
          
C           No precalculated information. Do we have a simple-type
C           geometry?

            IF ((L0.EQ.0).OR.(L0.EQ.1)) THEN
            
C             Create an identity coordinate system for the current
C             object

              CALL GCRSYS (CORSYS, 0D0, 0D0, 0D0, 1D0)

C             The geometry data is stored in IGEODL at the starting
C             address given by L1!

              L2 = IGEOM(OIGEODL)
              IF (L2.EQ.0) THEN
                WRITE (*,*) 'Geometry data array IGEODL undefined!'
                STOP
              END IF
              L2 = L(L2)+L1-1
              
C             Call the standard wrapper routine to get the desired
C             information!

              CALL NDE2XY (INODE,TRIA,X,Y)
              ISFBDM = GINWRP (DWORK(L2),CORSYS, X, Y)
            
            ELSE
            
C             The geometry is somehow more complex.
C             Fall back to the standard calculation routine ISFBDY.
          
              CALL NDE2XY (INODE,TRIA,X,Y)
              ISFBDM = ISFBDY (X,Y,IFBC,IGEOM,DGEOM)
              
            END IF ! L0=0 or L0=1
            
          END IF ! L0=12
          
        ELSE
        
C         General test. Find out the object that contains
C         our point - if there is one.
C         We make a linear search over all non-group objects
C         here.
          
          DO I=1,IGEOCN

C           Calculate the starting addresses of that object 
C           in the work-array
              
            INDI = L(LIARR) + KWORK(L(LIARR)+I-1)
            INDD = L(LDARR) + KWORK(L(LIARR)+IGEOCN+I-1)
            
C           Get the handle of the integer- and double-prec. parameter
C           block from IGEOMS(2,I) and IGEOMS(3,I), resp.; L0 gets
C           the type of the object.
C
C            L0 = IGEOM(1,I)
C            L1 = IGEOM(2,I)
C            L2 = IGEOM(3,I)

            L0 = IGEOM(OIGEOMS+3*(I-1))
            L1 = IGEOM(OIGEOMS+3*(I-1)+1)
            L2 = IGEOM(OIGEOMS+3*(I-1)+2)
            
C           Is it a composed geometry with precalculated information?

            IF (L0.EQ.12) THEN
            
C             Use the precalculated information to calculate the
C             information

              IF (GIPCMP (KWORK(L(L1)), 
     *                    DWORK(L(L2)), 
     *                    GIPWRP, TRIA,
     *                    KWORK(INDI),DWORK(INDD),
     *                    INODE).GT.0) THEN
     
C               Stop at the first found geometry containing the point

                ISFBDM = TNBC()+I
                GOTO 10

              END IF

            ELSE
            
C             No precalculated information. Do we have a simple-type
C             geometry (L0=0)?
C             Simple type geometries of type L0=1 which belong to a
C             group of objects are not handled here! For these objects
C             there is a "master"-object that takes care of the sub-
C             objects! That master-object is handles in ISFBDY.

              IF (L0.EQ.0) THEN

C               Create an identity coordinate system for the current
C               object

                CALL GCRSYS (CORSYS, 0D0, 0D0, 0D0, 1D0)

C               The geometry data is stored in IGEODL at the starting
C               address given by L1!

                L2 = IGEOM(OIGEODL)
                IF (L2.EQ.0) THEN
                 WRITE (*,*) 'Geometry data array IGEODL undefined!'
                 STOP
                END IF
                L2 = L(L2)+L1-1
                
C               Call the standard wrapper routine to get the desired
C               information!

                CALL NDE2XY (INODE,TRIA,X,Y)
                IF (GINWRP (DWORK(L2),CORSYS, X, Y).GT.0) THEN
                
C                 Stop at the first found geometry containing the point

                  ISFBDM = TNBC()+I
                  GOTO 10
                  
                END IF
              
              ELSE

C               The geometry is somehow more complex.
C               Fall back to the standard calculation routine ISFBDY.
              
                CALL NDE2XY (INODE,TRIA,X,Y)
                IF (ISFBDY (X,Y,I,IGEOM,DGEOM).GT.0) THEN

C                 Stop at the first found geometry containing the point

                  ISFBDM = TNBC()+I
                  GOTO 10

                END IF
              
              END IF ! L0=0 
              
            END IF ! L0=12

          END DO ! I
        
10        CONTINUE        
        
        END IF

      END IF

      END
      
************************************************************************
* Prescribe Dirichlet boundary conditions on fictitious boundary
*
* This is the node-dependent version of FBINDT.
*
* This routine is called by the framework when implementing Dirichlet 
* boundary conditions in a point of the fictitious boundary domain.
* It can be used to prescribe e.g. tangential velocity initiated
* by rotation.
* ITYP describes the information that has to be returned in the
* return value.
*
* In:
*  ITYP   - 0=X-velocity, 1=Y-velocity
*  INODE  - Node number of the node
*  TRIA   - array [1..SZTRIA] of integer
*           Triangulation structure STRIA.
*  TIMENS - Current time in instationary Navier-Stokes calculation
*           =0 for stationary simulation.
*  RE     - Reynolds number; from the DAT file
*  LIARR  - integer; might be 0 depending on FBDPRC
*           Handle to integer array with precalculated information.
*  LDARR  - integer; might be 0 depending on FBDPRC
*           Handle to double array with precalculated information.
*  IPARAM - array [1..*] of integer 
*  DPARAM - array [1..*] of integer 
*           TIntAssembly/TDoubleAssembly assembly structures; gives
*           additional information about the discretization.
*           This is passed to user defined callback routines so that 
*           they can access more detailed information about the
*           problem. Not used in this routine.
*  IGEOM  - array [1..*] of integer 
*  DGEOM  - array [1..*] of double 
*           Integer- and double-precision parameter blocks with
*           geometry information. Passed to boundary
*           routines. Not used in this routine.
* 
* Out:
*  Return value = the desired information
*
* For a detailed documentation about the meaning of the different
* quantities here, see FBINDT!
************************************************************************

      DOUBLE PRECISION FUNCTION FBMIND (ITYP, INODE, TIMENS,RE, TRIA, 
     *                                  LIARR, LDARR,
     *                                  IPARAM,DPARAM,IGEOM,DGEOM)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'stria.inc'
      
C     parameters

      INTEGER INODE, LIARR, LDARR, ITYP
      INTEGER TRIA(SZTRIA),IPARAM(*),IGEOM(*)
      DOUBLE PRECISION TIMENS,RE,DPARAM(*),DGEOM(*)
      
C     externals

      DOUBLE PRECISION FBINDT
      EXTERNAL FBINDT

C     local variables

      DOUBLE PRECISION X,Y

C     Standard implementation: 
C     Determine coordinates and call FBINDT!

      IF ((LIARR.EQ.0).OR.(LDARR.EQ.0)) THEN
        CALL NDE2XY (INODE,TRIA,X,Y)
        FBMIND = FBINDT (ITYP,X,Y,RE,INODE,TRIA,IPARAM,DPARAM,
     *                   IGEOM,DGEOM)
      END IF

99999 END
      
      
************************************************************************
* FBMNML - Get normal vector of fictitious boundary
*
* This is the node-dependent version of FBDNML.
*
* This routine is called by the framework to compute analytical normal
* vectors of the fictitious boundary. To a given node INODE -
* which may be inside the fictitious boundary or not, depending
* on the approximation - this routine has to calculate a normalised
* vector (XN,YN) perpendicular to the fictitious boundary domain IFBC.
*
* In:
*  INODE  - Node number of the node 
*  TRIA   - array [1..SZTRIA] of integer
*           Triangulation structure STRIA.
*  LIARR  - integer; might be 0 depending on FBDPRC
*           Handle to integer array with precalculated information.
*  LDARR  - integer; might be 0 depending on FBDPRC
*           Handle to double array with precalculated information.
*  IFBC   - Number of fictitious boundary component to calculate
*           the volume for.
*           0 = automatically determine a fictitious boundary
*               component which boundary is taken for the calculation
*               of the normal vector (i.e. by analytically analysing
*               the distance)
*           NBCT+1..NBCT+NFBDYC = Number of fictitious boundary
*               component which is used for the calculation of the
*               normal vector.
*  IINFO  - array [1..*] of integer
*  DINFO  - array [1..*] of double precision
*           Parameter blocks describing the geometry
*
* Out:
*  XN,YN  - normal vector of the fictitious boundary component IFBC
*           in the point (X,Y); must be normalised to length=1D0 !
* 
* If calculation of (XN,YN) is possible, this routine should set
* (XN,YN)=(0D0,0D0) to indicate this to the framework.
*
* For a detailed documentation about the meaning of the different
* quantities here, see FBDNML!
************************************************************************

      SUBROUTINE FBMNML (INODE, TRIA, LIARR, LDARR, IFBC, XN, YN,
     *                   IINFO,DINFO)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'stria.inc'
      
C     parameters

      INTEGER INODE, LIARR, LDARR, IFBC,DINFO(*)
      INTEGER TRIA(SZTRIA),XN,YN,IINFO(*)
      
C     local variables

      DOUBLE PRECISION X,Y

C     Standard implementation:
C     Determine coordinates and call FBINML!

      IF ((LIARR.EQ.0).OR.(LDARR.EQ.0)) THEN
        CALL NDE2XY (INODE,TRIA,X,Y)
        CALL FBDNML (X,Y,IFBC,XN,YN,IINFO,DINFO)
      END IF

      END
      
************************************************************************
* FBMDST - Get distance of point to fictitious boundary interface  
*
* This is the node-dependent version of FBDDST.
*
* This routine is called by the framework to compute the minimum 
* distance of a node INODE to the interface of a fictitious
* boundary component. 
*
* In:
*  INODE  - Node number of the node 
*  TRIA   - array [1..SZTRIA] of integer
*           Triangulation structure STRIA.
*  LIARR  - integer; might be 0 depending on FBDPRC
*           Handle to integer array with precalculated information.
*  LDARR  - integer; might be 0 depending on FBDPRC
*           Handle to double array with precalculated information.
*  IFBC   - Number of fictitious boundary component to calculate
*           the distance to.
*           0 = Compute the minimum distance to all fictitious
*               boundaty interfaces
*           NBCT+1..NBCT+NFBDYC = Number of fictitious boundary
*               component where to calculate the distance to.
*  IINFO  - array [1..*] of integer
*  DINFO  - array [1..*] of double precision
*           Parameter blocks describing the geometry
* Out:
*  DIST   - Distance of (X,Y) to the interface of a/the fictitious
*           boundary
* 
* Return:
*  >  0   - If the distance was successfully computed
*  <= 0   - If the computation is not possible
* 
* For a detailed documentation about the meaning of the different
* quantities here, see FBDDST!
************************************************************************

      INTEGER FUNCTION FBMDST (INODE, TRIA, LIARR, LDARR, IFBC, DIST,
     *                         IINFO,DINFO)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'stria.inc'
      INCLUDE 'sinigeometry.inc'
      INCLUDE 'sgeometries.inc'
      
C parameters

      INTEGER INODE, LIARR, LDARR, IFBC
      INTEGER TRIA(SZTRIA),IINFO(*)
      DOUBLE PRECISION DIST,DINFO(*)
      
C local variables

      DOUBLE PRECISION X,Y,TMP
      INTEGER INDI, INDD, I
      
C externals

      INTEGER FBDDST,GDPCMP,GDPWRP,GDSWRP
      EXTERNAL FBDDST,GDPCMP,GDPWRP,GDSWRP

      INTEGER ICIRTP,IGEOCN,L1,L2,L0
      DOUBLE PRECISION CORSYS(SZCORS)
      
C     Get rotation, radius,... from the parameter blocks

      ICIRTP = IINFO(OICIRTP)
      IGEOCN = IINFO(OIGEOCN)
      
C Standard implementation: Determine coordinates and call FBDDST!

      IF ((LIARR.EQ.0).OR.(LDARR.EQ.0).OR.(ICIRTP.LT.50)) THEN
      
        CALL NDE2XY (INODE,TRIA,X,Y)
        FBMDST = FBDDST (X,Y,IFBC, DIST, IINFO,DINFO)
        
      ELSE
      
C Composed geometries. Loop through all geometries and
C determine the minimum distance.

        TMP = 1D99

        DO I=1,IGEOCN

C Calculate the starting addresses of that object in the work-array

          INDI = L(LIARR) + KWORK(L(LIARR)+I-1)
          INDD = L(LDARR) + KWORK(L(LIARR)+IGEOCN+I-1)

C There's precalculated information available - use it!

C         Get the handle of the integer- and double-prec. parameter
C         block from IGEOMS(1,I) and IGEOMS(2,I), resp.; L0 gets
C         the type of the object.
C
C          L0 = IGEOMS(1,I)
C          L1 = IGEOMS(2,I)
C          L2 = IGEOMS(3,I)

          L0 = IINFO(OIGEOMS+3*(I-1))
          L1 = IINFO(OIGEOMS+3*(I-1)+1)
          L2 = IINFO(OIGEOMS+3*(I-1)+2)

C         Is precalculated information available?

          IF (L0.EQ.12) THEN

C           Use precalculated information to get the distance for
C           that object.

            IF (GDPCMP (KWORK(L(L1)), 
     *                  DWORK(L(L2)),
     *                  GDPWRP, TRIA, 
     *                  KWORK(INDI), DWORK(INDD), INODE, DIST).GT.0) 
     *      THEN
              IF (ABS(DIST).LT.ABS(TMP)) THEN
                TMP = DIST
              END IF
            END IF
            
          ELSE
          
C           No precalculated information. Do we have a simple-type
C           geometry?
C           Don't handle sub-objects of object-groups! (L0=1)

            IF (L0.EQ.0) THEN
            
C             Create an identity coordinate system for the current
C             object

              CALL GCRSYS (CORSYS, 0D0, 0D0, 0D0, 1D0)

C             The geometry data is stored in IGEODL at the starting
C             address given by L1!

              L2 = IINFO(OIGEODL)
              IF (L2.EQ.0) THEN
                WRITE (*,*) 'Geometry data array IGEODL undefined!'
                STOP
              END IF
              L2 = L(L2)+L1-1
              
C             Call the standard wrapper routine to get the desired
C             information!

              CALL NDE2XY (INODE,TRIA,X,Y)
              IF (GDSWRP (DWORK(L2),CORSYS, X, Y,DIST).GE.0) THEN
                IF (ABS(DIST).LT.ABS(TMP)) THEN
                  TMP = DIST
                END IF
              END IF
            
            ELSE
            
C             The geometry is somehow more complex.
C             Fall back to the standard calculation routine ISFBDY.
          
              CALL NDE2XY (INODE,TRIA,X,Y)
              DIST = FBDDST (X,Y,IFBC, DIST, IINFO,DINFO)
              
              IF (ABS(DIST).LT.ABS(TMP)) THEN
                TMP = DIST
              END IF

            END IF ! L0=0 
            
          END IF ! L0=12

        END DO

        FBMDST = TMP

      END IF

      END

************************************************************************
* FBMMON - Get monitor function value of fictitious boundary domain
*
* This is the node-dependent version of FBDMON.
*
* This routine is called by the framework to compute a value for the
* monitor function in the point (X,Y) when performing grid adation. 
* It is called once with IFBC=0 to compute a general value for the
* monitor function and once for every fictitious boundary domain
* IFBC=NBCT+1..NBCT+NFBDYC. The monitor function itself is computed
* by taking the minimum value of all returned values in that point.
*
* In:
*  INODE  - Node number of the node 
*  TRIA   - array [1..SZTRIA] of integer
*           Triangulation structure STRIA.
*  LIARR  - integer; might be 0 depending on FBDPRC
*           Handle to integer array with precalculated information.
*  LDARR  - integer; might be 0 depending on FBDPRC
*           Handle to double array with precalculated information.
*  IFBC   - Number of fictitious boundary component to calculate
*           the volume for.
*           0 = General computation of the monitor function value;
*               the routine can decide how to calculate
*           NBCT+1..NBCT+NFBDYC = Number of fictitious boundary
*               component whose weight should be calculated
*           If one case is not handled by this routine, 1D0 should
*           be returned.
*  EPS0   - minimum monitor function values; values below that
*           will be truncated.
*  IGEOM  - array [1..*] of integer
*  DGEOM  - array [1..*] of double precision
*           Parameter blocks describing the geometry
*
* Return:
*  Value of monitor function, either for all boundary components (if
*  IFBC=0) or for the fictitious boundary component IFBC.
* 
* For a detailed documentation about the meaning of the different
* quantities here, see FBDMON!
************************************************************************

      DOUBLE PRECISION FUNCTION FBMMON (INODE, TRIA, LIARR, LDARR,
     *                                  IFBC,EPS0,IGEOM,DGEOM)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'stria.inc'
      INCLUDE 'sinigeometry.inc'
      INCLUDE 'sgeometries.inc'
      
C     parameters

      INTEGER INODE, LIARR, LDARR, IFBC
      INTEGER TRIA(SZTRIA),IGEOM(*)
      DOUBLE PRECISION DGEOM(*),EPS0
      
C     local variables

      DOUBLE PRECISION X,Y,TMP,DIST
      INTEGER INDI, INDD, I

C     externals

      INTEGER FBDMON,GDPCMP,GDPWRP,GDSWRP
      EXTERNAL FBDMON,GDPCMP,GDPWRP,GDSWRP

      INTEGER ICIRTP,IGEOCN,L1,L2,L0
      DOUBLE PRECISION DMNISF,DMNOSF
      DOUBLE PRECISION CORSYS(SZCORS)
      
C     Get rotation, radius,... from the parameter blocks

      ICIRTP = IGEOM(OICIRTP)
      IGEOCN = IGEOM(OIGEOCN) 
      DMNISF = DGEOM(ODMNISF)
      DMNOSF = DGEOM(ODMNOSF)
      
C     Standard implementation: Determine coordinates and call FBDDST!

      IF ((LIARR.EQ.0).OR.(LDARR.EQ.0)) THEN
      
        CALL NDE2XY (INODE,TRIA,X,Y)
        FBMMON = FBDMON (X,Y,IFBC, EPS0, IGEOM,DGEOM)
        
      ELSE
        
C       Multiple geometries. Loop through all geometries and
C       determine the minimum monitor function value.

        TMP = 1D0

        DO I=1,IGEOCN

C         Calculate the starting addresses of that object in the
C         work-array

          INDI = L(LIARR) + KWORK(L(LIARR)+I-1)
          INDD = L(LDARR) + KWORK(L(LIARR)+IGEOCN+I-1)

C         Get the handle of the integer- and double-prec. parameter
C         block from IGEOMS(1,I) and IGEOMS(2,I), resp.; L0 gets
C         the type of the object.
C
C          L0 = IGEOMS(1,I)
C          L1 = IGEOMS(2,I)
C          L2 = IGEOMS(3,I)

          L0 = IGEOM(OIGEOMS+3*(I-1))
          L1 = IGEOM(OIGEOMS+3*(I-1)+1)
          L2 = IGEOM(OIGEOMS+3*(I-1)+2)

C         Is precalculated information available?

          IF (L0.EQ.12) THEN

C           There's precalculated information available - use it!

            IF (GDPCMP (KWORK(L(L1)), 
     *                  DWORK(L(L2)),
     *                  GDPWRP, TRIA, 
     *                  KWORK(INDI), DWORK(INDD), INODE, DIST).GT.0) 
     *        THEN
     
C             Correct distance depending on the sign with the help of the
C             multiplication factors. This enables us to absorb improper
C             geometry scaling manually.

              IF (DIST.LT.0D0) THEN
                DIST=DIST*DMNISF
              ELSE
                DIST=DIST*DMNOSF
              END IF
              
C             Save the distance in TMP
              
              IF (I.EQ.1) THEN
                TMP = ABS(DIST)
              ELSE
                IF (IGEOM(OIOVMNT).EQ.0) THEN
                  TMP = MIN(TMP,ABS(DIST))
                ELSE
                  TMP = TMP*ABS(DIST)
                END IF
              END IF
            END IF

          ELSE
          
C           No precalculated information. Do we have a simple-type
C           geometry?
C           Don't handle sub-objects of object-groups! (L0=1)

            IF (L0.EQ.0) THEN
            
C             Create an identity coordinate system for the current
C             object

              CALL GCRSYS (CORSYS, 0D0, 0D0, 0D0, 1D0)

C             The geometry data is stored in IGEODL at the starting
C             address given by L1!

              L2 = IGEOM(OIGEODL)
              IF (L2.EQ.0) THEN
                WRITE (*,*) 'Geometry data array IGEODL undefined!'
                STOP
              END IF
              L2 = L(L2)+L1-1
              
C             Call the standard wrapper routine to get the desired
C             information!

              CALL NDE2XY (INODE,TRIA,X,Y)
              IF (GDSWRP(DWORK(L2),CORSYS, X,Y, DIST).GE.0) THEN

C               Correct distance depending on the sign with the help of the
C               multiplication factors. This enables us to absorb improper
C               geometry scaling manually.

                IF (DIST.LT.0D0) THEN
                  DIST=DIST*DMNISF
                ELSE
                  DIST=DIST*DMNOSF
                END IF
                
C               Save the distance in TMP
              
                DIST = MAX(EPS0,ABS(DIST))
                IF (I.EQ.1) THEN
                  TMP = ABS(DIST)
                ELSE
                  IF (IGEOM(OIOVMNT).EQ.0) THEN
                    TMP = MIN(TMP,ABS(DIST))
                  ELSE
                    TMP = TMP*ABS(DIST)
                  END IF
                END IF
              END IF
              
            END IF ! L0 = 0

          END IF ! L0 = 12
          
        END DO

        FBMMON = MAX(EPS0,TMP)

      END IF
      
C      FBMMON = MIN(FBMMON,8D0*EPS0*(1-X/1.1)+1*(X/1.1))      

C     Monitor function is between 0 and 1. Scale it to be
C     between EPS0 and 1

      IF (IGEOM(OISCMNT).NE.0) THEN
        FBMMON = EPS0 + (1D0-EPS0)*FBMMON
      END IF

      END
