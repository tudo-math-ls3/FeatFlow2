************************************************************************
* This file defines a parameter structure for maintaining geometries
* inside of the solver.
*
* The basic function is as follows:
*  - Using RDGEOM, geometry data is loaded from a DAT file into the
*    COMMON bloock variables defined in the cinigeometry.inc file.
*  - Using INIGEO, the geometry data is transferred from the
*    COMMON block variables to the parameter structures defined here
*  - When geometry data is changed (e.g. position,...), the routine
*    UPTGEO should be called to update geometry data.
*    This is especially important if complex geometries are used,
*    as UPDGEO updates the coordinate systems of all sub-objects.
*  - The structure is passed in the program to the subroutines, so
*    the COMMON-block variables need not to be touched.
*    For updating the structure, UPDGEO can be used.
*  - When the main program is finished, DONGEO releases dynamically
*    allocated information from the structure.
*
* The realization of the structures as arrays is made e.g. by:
*
*   IGEODT = array [1..SZGEOI] of integer
*          = TIGeometryData
*   DGEODT = array [1..SZGEOD] of double
*          = TDGeometryData
*
* The Oxxxxx-constants below can be used to access the different
* variables inside of the structure-array, e.g.:
*
*   ICIRTP = IPARAM(OICIRTP)
*
* The two structure array have the following content:
*
* TIGeometryData = record
*
*   // For simple geometries, the following type identifier
*   // defines the obstacle directly
*
*   ICIRTP - Type of obstacle (circle, square, ellipse, ...)
*
*   // The following variables are used to modify the way, the
*   // triangulation is adapted to the geometry:
*
*   IFBMNT - Defines the type of the monitor function for adapting
*            the triangulation to the geometry.
*            This is passed to the FBDMON function for the
*            fictitious boundary definition and more or less user
*            defined.
*   ISCMNT - Whether or not to rescale the monitor function for
*            adapting the triangulation to the geometry.
*            This is passed to the FBDMON function for the
*            fictitious boundary definition and more or less user
*            defined.
*   IOVMNT - Method of how to overlay monitor functions of multiple
*            geometries.
*            =0: Standard method, take minimum distance
*            =1: Digital filter technique, multiply all distances
*
*   // A user defined block where arbitrary information can be stored
*
*   IGEOUD - array [1..16] of integer
*            A user defined block where arbitrary information can
*            be stored, eng. by fictitious boundary routines.
*
*   // More complex objects can be composed by simple ones using
*   // composed geometries. ICIRTP then contains an identifier
*   // that tells the framework to use the following two variables:
*
*   IGEOCN - Saves the number of present geometries, i.e.
*            the number of handles allocated in IGEOMS.
*   IGEOIL - Handle to an array [1..*] of integer on the heap
*            to store additional information, or =0 if not used
*   IGEODL - Handle to an array [1..*] of double on the heap
*            to store additional information
*   IGEOMS - array [1..3,1..100] of integer
*            Saves handles/information for up to 100 (composed) geometry
*            objects. The first entry in each 3-tuple decides on
*            how the framework handles the object:
*
*            IGEOMS(1,.) = 0
*              This is a simple-type object without precomputed
*              information. IGEOMS(2,.) gives the starting address of
*              the corresponding TGEOxxx-structure in the IGEODL
*              array.
*            IGEOMS(1,.) = 1
*              Simple-type object that belongs to a group of objects.
*              No precalculated information.
*              IGEOMS(2,.) gives the starting address of the
*              corresponding TGEOxxx-structure in the IGEODL
*              array. IGEOMS(3,.) gives the number of the guiding
*              object that defines the group.
*
*            IGEOMS(1,.) = 3..9 - reserved
*            IGEOMS(1,.) = 10
*              Simple composed geometry, no precalculated information.
*              IGEOMS(2,.) is the starting address in IGEODL of a
*              TGEOSIM structure that defines the coordinate system of
*              the whole group. IGEOMS(3,.) is the starting address of
*              a user defined block in IGEOIL that gives additional
*              information about the group (e.g. how many and which
*              objects are contained in the group,...)
*            IGEOMS(1,.) = 11
*              Complex-type composed geometry object.
*              The object noes not use precalculated information.
*              IGEOMS(2,.) = handle to an array in KWORK which
*                represents TGEOCOMIDX
*              IGEOMS(3,.) = handle to an array in DWORK which
*                represents the corresponding TGEOCOM.
*            IGEOMS(1,.) = 12
*              Complex-type composed geometry object.
*              The object uses precalculated information.
*              IGEOMS(2,.) = handle to an array in KWORK which
*                represents TGEOCOMIDX
*              IGEOMS(3,.) = handle to an array in DWORK which
*                represents the corresponding TGEOCOM.
*
*   // A flag tells the routines whether we are calculating stationary
*   // or nonstationary
*
*   INSGEO - =0, if the simulation is stationary (standard)
*            =1, if the simulation is nonstationary
*
* end
*
* TDGeometryData = record
*
*   // The following parameters define the basic position / rotation
*   // of a geometry object:
*
*   DCXPOS - X-Position of the obstacle midpoint
*   DCYPOS - Y-Position of the obstacle midpoint
*   DCRAD  - Radius of a circle / half length of an edge on a square
*   DCRDX  - X-Radius of an ellipse
*   DCRDY  - Y-Radius of an ellipse
*   DCROT  - Roration of the obstacle, 0..360°.
*
*   DCRLX,
*   DCRLY  - DCRLX and DCRLY define a "radius relaxation" - a small
*            value that is added to the "radius" of simple-type
*            objects like circle or square when checking which points
*            are inside an object and which not. This gives usually
*            better results as it prevents some rounding difficulties.
*            Normally, these values are set to 0.0, but values like
*            0.00001 without negative influence of the computation.
*
*   // For saving computational time the sin/cos of the rotation is
*   // precalculated. Every time DCROT is changed, these values
*   // must be updated!
*
*   DCRSIN - SIN(DCROT*PI/180); must always be updated if DCROT is changed
*   DCRCOS - COS(DCROT*PI/180); must always be updated if DCROT is changed
*
*   // The following variables are used to modify the way, the
*   // triangulation is adapted to the geometry:
*
*   DMNISF,
*   DMNOSF - Steepness factors for monitor function inside/outside of
*            objects; is multiplied with the monitor function value
*            to correct improper large values through distance
*            calculation when large geometries are used.
*
*   // The following variables are used during instationary simulations
*   // as status variables
*
*   GEOTIM - Geometry time.
*            During nonstationary simulations, this variable defines
*            the current time stamp of the geometry. This is typically
*            set to the simulational time of the whole simulation.
*            In stationary simulations, GEOTIM is set to 0.0.
*
*   // A user defined block where arbitrary information can be stored
*
*   DGEOUD - array [1..16] of integer
*            A user defined block where arbitrary information can
*            be stored, eng. by fictitious boundary routines.
* end
*
************************************************************************

C     Maximum number of composed geometries.
C     Each composed geometry is composed of multiple single objects.

      INTEGER MAXGEO
      PARAMETER (MAXGEO  = 100)

C     TIGeometryData:
C     ---------------

      INTEGER OICIRTP,OIGEOCN,OIGEOMS,OIFBMNT,OISCMNT,OINSGEO,OIGEOUI
      INTEGER OIGEOIL,OIGEODL,OIOVMNT
      
      PARAMETER (OICIRTP = 1)
      PARAMETER (OIFBMNT = 2)
      PARAMETER (OISCMNT = 3)
      PARAMETER (OIGEOCN = 4)
      PARAMETER (OIGEOUI = 5)
      PARAMETER (OIGEOIL = 6)
      PARAMETER (OIGEODL = 7)
      PARAMETER (OIGEOMS = 8)
      PARAMETER (OINSGEO = 9)
      PARAMETER (OIOVMNT = 10)
      
      INTEGER SZGEOI
      PARAMETER (SZGEOI  = 16+16+3*MAXGEO)

C     TDGeometryData:
C     ---------------

      INTEGER ODCXPOS,ODCYPOS,ODCROT
      INTEGER ODCRAD ,ODCRDX,ODCRDY,ODCRSIN,ODCRCOS
      INTEGER ODCRLX ,ODCRLY
      INTEGER ODMNISF,ODMNOSF,OGEOTIM,ODGEOUD

      PARAMETER (ODCXPOS = 1)
      PARAMETER (ODCYPOS = 2)
      PARAMETER (ODCROT  = 3)

      PARAMETER (ODCRAD  = 4)
      PARAMETER (ODCRDX  = 5)
      PARAMETER (ODCRDY  = 6)
      PARAMETER (ODCRLX  = 7)
      PARAMETER (ODCRLY  = 8)
      PARAMETER (ODCRSIN = 9)
      PARAMETER (ODCRCOS = 10)
      PARAMETER (ODMNISF = 11)
      PARAMETER (ODMNOSF = 12)
      
      PARAMETER (OGEOTIM = 13)
      
      PARAMETER (ODGEOUD = 14)

      INTEGER SZGEOD
      PARAMETER (SZGEOD  = 16+16)
