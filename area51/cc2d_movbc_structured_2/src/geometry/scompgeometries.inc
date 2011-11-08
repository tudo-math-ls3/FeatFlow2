************************************************************************
* This file defines constants for the collection of composed geometries.
* By composing a set of simple geometries together, a composed
* geometry arises. The outines here allow to directly combine
* such a set of simple geometries with the help of wrapper functions.
*
* This include-file is an extension to the SGEOMETRIES.INC-file, so
* SGEOMETRIES.INC must be included prior to this file!
************************************************************************

************************************************************************
* Composed geometry structure
*
* A composed geometry is based on two arrays, which must be passed
* to all routines:
*
* TGEOCOMIDX = record
*   COUNT   : integer
*   INDEX   : array [1..COUNT] of integer
* end
*
* realized as:
*   GEOCOMIDX = array [1..COUNT+1] of integer
*
* TGEOCOM = record
*   MAIN    : TGEOSIM
*   OBJS    : array [1..COUNT] of TGEOxxx
* end
*
* realized as:
*   GETCOM = array [1..SZGSIM+*] of double
*
* The TGEOCOM-structure is the main structure defining the composed
* object and its components. It consists of two parts:
* 1.) MAIN - This is the description of the composed object in the style
*      of a simple-type object. The TYPE-identifier identifies the
*      object as composed-object. All other parameters (coordinate
*      of the reference point, rotation,...) affect the whole composed
*      object.
* 2.) OBJS - This is a list of all sub-objects contained in the composed
*      object. All objects are realised by their appropriate structure.
*      The coordinates of the reference point and the rotation of the
*      sub-elements are all relative to thereference point of the
*      composed object!
*      (e.g. an object reference-point of (0,0) places the sub-object
*      at the same position as the complete composed object)
* Obviously the sizes of the structures in the OBJS-array differs
* depending on the object. To provide direct access to each sub-object,
* the TGEOCOMIDX-structure is used. COUNT describes the number of
* sub-objects in the composed object and INDEX describes the starting
* index of each sub-object in the GETCOM-array.
*
* Remark:
* 1.) The routines provided here can basically be divided into the same
*  groups as described in the simple-object case. Nevertheless because
*  of difficulties in determining the appropriate values (like
*  normal vector,...) the routines don't always work reliable for all
*  type of combination of simple objects.
* 2.) For easy access to the simple objects, wrapper-routines are
*  necessary to call the appropriate routine depending on the object
*  type. A standard implementation of these wrapper-routines
*  is provided by the GxxWRP-routines. If the user wants to use more
*  complex "simple" geometries, wither these routines must be changed
*  or different routines must be provided to perform this wrapping.
*  For this purpose every composed-geometry-routine allowes to
*  specify such a wrapper routine.
* 3.) Composed geometries allow the subtraction of one object from
*  another. Every sub-object that should be geometrically subtracted
*  from the group must have the "INVER" flag set. Every subtracted
*  object only subtracts one layer from the group, i.e.
*  if there are 3 objects overlaying each other and only one has
*  the INVERT flag set, there remains one object!
************************************************************************

************************************************************************
* Offsets of the entries in the TGEOCOMIDX-structure-array
************************************************************************

      INTEGER OGCCNT,OGCIDX
      PARAMETER (OGCCNT =  1)
      PARAMETER (OGCIDX =  2)

************************************************************************
* Offsets of the entries in the TGEOCOM-structure-array
************************************************************************

      INTEGER OGCMAIN,OGCOBJS
      PARAMETER (OGCMAIN =  1)
      PARAMETER (OGCOBJS =  SZGSIM+1)

      