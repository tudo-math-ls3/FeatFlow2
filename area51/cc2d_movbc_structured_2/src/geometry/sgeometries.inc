************************************************************************
* This file defines constants for the collection of geometries.
* These constants represent index-numbers inside of the structure-array
* and allow the access of the different variables of the structure
* by a name.
************************************************************************

************************************************************************
* Basic geometry structure
*
* The basic structure of all simple geometry is defined as follows:
*
* TGEOSIM = record
*   TYPE    : double (only takes integer-values)
*   TAG     : double
*   INVER   : double (only takes integer values)
*   REFX    : double
*   REFY    : double
*   ROT     : double
*   ROTSIN  : double
*   ROTCOS  : double
*   SCALE   : double
* end
*
* Parameters:
*   TYPE    = Type-ID of the geometry
*   TAG     = a user-definable tag
*   INVER   = whether the object is "inversed"
*             0=normal object, 1=inversed object
*   REFX    = Reference X-position
*   REFY    = Reference Y-position
*   ROT     = Angle of rotation (0..360)
*   ROTSIN  = SIND(ROT) = SIN(ROT*PI/180D0)
*   ROTCOS  = COSD(ROT) = COS(ROT*PI/180D0)
*   SCALE   = scaling factor; 1=standard=no scaling
*             e.g. 0.5=object is half of its orig. size.
*
* Remarks:
* 1.) The type of the geometry structure is identified with the
*     TYPE entry. The number stored here is integer-valued and
*     unique for every object-type
* 2.) The TAG is a user-defined value. The user might save e.g. the
*     element number or an internal number of the object here.
* 3.) (REFX,REFY) are the coordinates of a "reference point" of that
*     object (e.g. midpoint of a circle). The position of the reference
*     point in the object (if it's a midpoint, a corner,...) is 
*     object-specific.
* 4.) ROT defines the rotation of the object around the reference
*     point. For saving computational time, ROTSIN/ROTCOS save the
*     SIN- and COS-value of ROT. Whenever the user changes ROT,
*     ROTSIN/ROTCOS have to be updated by the user, too!
*     For user convenience, the rotation is used as degrees 0..360
*     instead of arc length (0..2*PI). Therefore ROTSIN/ROTCOS
*     has to be calculated by SIND/COSD instead of SIN/COS
*     (or use SIN/COS (x*PI/180D0), respectively, as SIND/COSD is not
*     available on all compilers)!
* 5.) For initialisation all elements in the structure can be set to 0
*     except for SCALE, which must be set to 1.
* 6.) SCALE gives the scaling of the coordinate system.
*     SCALE=1 means no scaling, SCALE=0.5 means halfing the size of
*     the object,...
* 7.) INVER allowes the object to be specified as "inversed" - i.e.
*     the inside is treated as outside and vice versa. This is realized
*     by inverting the sign of the distance, the normal vector,
*     switching 1<->0 when testing of being inside/outside of the
*     object,...
* 8.) The values starting at REFX share the coordinate system structure
*     defined below. For convenience, a constant OCSYS is defined
*     that point to REFX inside of the structure.
*
* This is realised as a double-array:
*
*   GEOSIM = array [1..SZGSIM] of double precision
*
* with GEOSIM[OTYPE] = TYPE, ...
* 
* The actual geometry structure that describes more complex objects
* is build as an extension to this structure, occupying the
* array-positions SZGSIM+1,SZGSIM+2,... depending on the object.
*  
* Remark: The SZGxxx-constants are used for defining the size of the
*  structures for better integration into Fortran programs.
************************************************************************

      INTEGER OTYPE, OTAG, OREFX, OREFY, OROT, OROTSIN, OROTCOS, OSCALE
      INTEGER OINVER, OCSYS
      INTEGER SZGSIM
      PARAMETER (OTYPE   =  1)
      PARAMETER (OTAG    =  2)
      PARAMETER (OINVER  =  3)
      PARAMETER (OREFX   =  4)
      PARAMETER (OREFY   =  5)
      PARAMETER (OROT    =  6)
      PARAMETER (OROTSIN =  7)
      PARAMETER (OROTCOS =  8)
      PARAMETER (OSCALE  =  9)

C Coordinate system starts at REFX:

      PARAMETER (OCSYS   = OREFX)

C We reserve some more space for future developments:

      PARAMETER (SZGSIM  = 10)

************************************************************************
* Coordinate sytstem structure
*
* The following structure is used for defining a coordinate system.
* It specifies the coordinates of the origin as well as the rotation
* of the coordinate system.
*
* TCORSYS = record
*   ORIGX    : double
*   ORIGY    : double
*   ORIGRT   : double
*   ORISIN   : double
*   ORICOS   : double
*   ORISCL   : double
* end
*
* with:
*  ORIGX  - X-coordinate of origin of coordinate system
*  ORIGY  - Y-coordinate of origin of coordinate system
*  ORIGRT - rotation of coordinate system; 0..360°
*  ORISIN - saves SIND(ORIGRT) = SIN(ORIGRT*PI/180D0)
*  ORICOS - saves COSD(ORIGRT) = COS(ORIGRT*PI/180D0)
*  ORSCL  - scaling factor; 1=standard=no scaling
*           e.g. 0.5=object is half of its orig. size.
*
* It's realized as an array:
*
*    CORSYS = array [1..SZCORS] of double
*
* The structure is part of the TGEOSIM-structure, starting at position
* GEOSIM(OREFX). The standard coordinate system is described by
* (ORIGX=0,ORIGY=0,ORIGRT=0).
* Whenever ORIGRT is changed, ORISIN/ORICOS has to be updated, too!
* This allowes quicker computation without calculating SIN/COS too
* often.
* For user convenience, the rotation is used as degrees 0..360
* instead of arc length (0..2*PI). Therefore ROTSIN/ROTCOS
* has to be calculated by SIND/COSD instead of SIN/COS
* (or by use of SIN/COS (x*PI/180D0), respectively, as SIND/COSD is not
* available on all compilers)!
*
* The routine GCRSYS offeres a standard initialization for such a
* coordinate system.
************************************************************************

      INTEGER ORIGX, ORIGY, ORIGRT, ORISIN, ORICOS, ORISCL
      INTEGER SZCORS
      
      PARAMETER (ORIGX   =  1)
      PARAMETER (ORIGY   =  2)
      PARAMETER (ORIGRT  =  3)
      PARAMETER (ORISIN  =  4)
      PARAMETER (ORICOS  =  5)
      PARAMETER (ORISCL  =  6)

C We reserve some more space for future developments:

      PARAMETER (SZCORS  =  6)

************************************************************************
* TYPE-identifiers for simple geometries.
* Identifiers 0..100 are reserved for simple-type geometries
* (circle, ellipse, ...)
************************************************************************

      INTEGER SGTNON,SGTCIR,SGTELL,SGTSQR,SGTREC,SGTSWV,SGTPSL
      
C Dummy-type, no geometry:
      
      PARAMETER (SGTNON =  0)
      
C Circle:

      PARAMETER (SGTCIR =  1)
      
C Ellipse
      
      PARAMETER (SGTELL =  2)
      
C Square

      PARAMETER (SGTSQR =  3)
      
C Rectangle
      
      PARAMETER (SGTREC =  4)
      
C SIN-Wave

      PARAMETER (SGTSWV = 5)

C PIE Slice

      PARAMETER (SGTPSL = 6)

************************************************************************
* Simple Geometry: Circle
*
* TGEOCIRCLE = record
*   BASIC   : TGEOSIM
*   RAD     : double
*             Radius of the circle
* end
************************************************************************

      INTEGER OCRAD, SZGCIR
      PARAMETER (OCRAD   = SZGSIM+1)

      PARAMETER (SZGCIR  = SZGSIM+1)
      
************************************************************************
* Simple geometry: Ellipse
*
* TGEOELLIPSE = record
*   BASIC   : TGEOSIM
*   RADX    : double
*             X-Radius of the ellipse
*   RADY    : double
*             Y-Radius of the ellipse
* end
************************************************************************
      
      INTEGER OERADX, OERADY, SZGELL
      PARAMETER (OERADX  = SZGSIM+1)
      PARAMETER (OERADY  = SZGSIM+2)

      PARAMETER (SZGELL  = SZGSIM+2)

************************************************************************
* Simple Geometry: Square
*
* TGEOSQUARE= record
*   BASIC   : TGEOSIM
*   LEN     : double
*             Length of each edge
* end
************************************************************************

      INTEGER OSLEN, SZGSQR
      PARAMETER (OSLEN   = SZGSIM+1)

      PARAMETER (SZGSQR  = SZGSIM+1)
      
************************************************************************
* Simple geometry: Rectangle
*
* TGEORECTANGLE = record
*   BASIC   : TGEOSIM
*   LENX    : double
*             Length of the X-edges
*   LENY    : double
*             Length of the Y-edges
* end
************************************************************************
      
      INTEGER ORLENX, ORLENY, SZGREC
      PARAMETER (ORLENX  = SZGSIM+1)
      PARAMETER (ORLENY  = SZGSIM+2)

      PARAMETER (SZGREC  = SZGSIM+2)

************************************************************************
* Simple geometry: SIN-Wave
*
* Two SIN-Waves above each other, connected to a closed object.
*
* TGEOSINWAVE = record
*   BASIC   : TGEOSIM
*   WIDTH   : double
*             Width of the SIN-Wave
*   POS1    : double
*             Y-position of the upper SIN-Wave; should be > 0
*   HEIGH1  : double
*             Height of the upper SIN-Wave
*   PHAS1   : double
*             Phase-angle shift of the upper SIN-Wave; 0..360°
*   FREQ1   : double
*             Frequency of the upper SIN-Wave; =number of waves
*   POS2    : double
*             Y-position of the lower SIN-Wave; should be < 0
*   HEIGH2  : double
*             Height of the bottom SIN-Wave
*   PHAS2   : double
*             Phase-angle shift of the bottom SIN-Wave; 0..360°
*   FREQ2   : double
*             Frequency of the lower SIN-Wave; =number of waves
* end
************************************************************************

      INTEGER OWWID, OWHEI1, OWPHAS1, OWFREQ1, OWHEI2, OWPHAS2, OWFREQ2
      INTEGER OWPOS1, OWPOS2, SZGSWV
      PARAMETER (OWWID   = SZGSIM+1)
      PARAMETER (OWPOS1  = SZGSIM+2)
      PARAMETER (OWHEI1  = SZGSIM+3)
      PARAMETER (OWPHAS1 = SZGSIM+4)
      PARAMETER (OWFREQ1 = SZGSIM+5)
      PARAMETER (OWPOS2  = SZGSIM+6)
      PARAMETER (OWHEI2  = SZGSIM+7)
      PARAMETER (OWPHAS2 = SZGSIM+8)
      PARAMETER (OWFREQ2 = SZGSIM+9)

      PARAMETER (SZGSWV  = SZGSIM+9)

************************************************************************
* Simple geometry: PIE SLICE
*
* A user defined geometry made of line-segments that cover a part of
* a circle. The line segments form the boundary. The start- and end-
* point is connected to a reference point (e.g. (0,0)) to form
* a pie slice.
*
* TGEOPIESLICE = record
*   BASIC   : TGEOSIM
*   COUNT   : integer
*             Number of points forming the boundary interface
*   ICLOS   : =0: The pie slice is "open"
*             =1: The pie slice is "closed"
*             This is merely a flag for GMV-output to indicate
*             whether a line should be drawn between the starting- and
*             the ending-point
*   XBAS,
*   YBAS    : Coordinates of the "basis" point.
*             The "basis" point specifies any point that is "inside" of
*             the geometry. For circle-like objects this is 
*             normally (0,0).
*   ILMP    : Number of the leftmost point / point with the smallest
*             angle
*   IRMP    : Number of the rightmost point / point with the 
*             largest angle
*   POINTS  : array [1..2,1..COUNT] of double
*             X- and Y-Coordinates of the points that form the
*             interface. This array defines a list of lines forming
*             a polygon:
*               [(x1,y1),(x2,y2),...,(xn,yn)]
*             with lines between the points: 
*               (x1,y1)->(x2,y2), (x2,y2)->(x3,y3),...
*             If the polygon is closed (ICLOS=1), the last point
*             must be set to the same coordinates as the first point:
*               (xn,yn)=(x1,y1)
*             (This is not automatically assumed!)
*             The points must be stored counterclockwise!
*             All points must be different!
* end
*
* A pie slice of a circle has (in rotation coordinates) a point
* with a smallest and a point with a largest angle. Typically
* this is the starting and ending point of the polygon. The numbers
* of these both "boundary" points are stored in ILMP/IRMP.
*
* As this structure is somehow "dynamic" (depending on the number of
* points forming the boundary), the SZGPSL-constant only gives the
* size of the structure without the points. To calculate the actual
* size of the object, 2*(number of points) has to be added to SZGPSL.
************************************************************************

      INTEGER OPCNT, OPCLOS, OPXBAS, OPYBAS, OPPOIN, OPILMP, OPIRMP 
      INTEGER SZGPSL
      PARAMETER (OPCNT   = SZGSIM+1)
      PARAMETER (OPCLOS  = SZGSIM+2)
      PARAMETER (OPXBAS  = SZGSIM+3)
      PARAMETER (OPYBAS  = SZGSIM+4)
      PARAMETER (OPIRMP  = SZGSIM+5)
      PARAMETER (OPILMP  = SZGSIM+6)
      
      PARAMETER (OPPOIN  = SZGSIM+7)

      PARAMETER (SZGPSL  = SZGSIM+6)
