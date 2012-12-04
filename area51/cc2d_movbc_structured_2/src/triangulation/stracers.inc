************************************************************************
* This file describes the structure that is maintained by the tracer
* routines.
*
* "Tracers" are points inside of the domain with no mass that move
* with the fluid. All routines dealing with tracers use the following
* tracer structure to maintain the tracers:
*
*   TDTracers = record
*     TRCNT  : Number of tracers in the structure
*     TRPEL  : sqrt(Number of tracers per element) if tracers are
*              regularly distributed in all elements during the
*              initialisation phase;
*              =0, if the tracers are generally placed in the
*                  domain without any regular structure
*     TRDEF  : array [1..TRCNT] of TDTracer = record
*                TRXPOS : current X-position of the tracer
*                TRYPOS : current Y-position of the tracer
*                TRXVEL : X-velocity of the tracer
*                TRYVEL : Y-velocity of the tracer
*                TRELEM : Number of the element that currently contains
*                         the tracer;
*                         = 0 : if element is unknown
*                         =-1 : if tracer is outside of the domain and must
*                               be reintroduced anywhere
*                         <=-2: User defined tracer status. The tracer is
*                               ignored by the standard routines.
*                TRTAG  : A user defined tag assigned to the tracer.
*                         Using GMVTRS, the value in TRTAG is writen out
*                         as function value assigned to the tracer.
*              end
*   end
*
* The structure is realised as DOUBLE-array:
*
*   TDTRACERS : array [1..*] of double
*
* and each element in the structure can be accessed using the offset
* constants declared in this file.
************************************************************************

      INTEGER OTRXPOS,OTRYPOS,OTRXVEL,OTRYVEL,OTRELEM,OTRTAG
      INTEGER SZTRCR
      
C     Offset positions in the "substructure" TDTracer
      
      PARAMETER (OTRXPOS = 1)
      PARAMETER (OTRYPOS = 2)
      PARAMETER (OTRXVEL = 3)
      PARAMETER (OTRYVEL = 4)
      PARAMETER (OTRELEM = 5)
      PARAMETER (OTRTAG  = 6)
      
      PARAMETER (SZTRCR  = 6)
      
C     Offset positions in the main structure:

      INTEGER OTRCNT,OTRDEF,OTRPEL
      INTEGER SZTRCS
      
      PARAMETER (OTRCNT = 1)
      PARAMETER (OTRPEL = 2)
      PARAMETER (OTRDEF = 3)
      
C     The actual "size" of the tracers structure is 
C         SZTRCS + SZTRCR*TRCNT
      
      PARAMETER (SZTRCS = 2)
      
