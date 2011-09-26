************************************************************************
* This file describes the structure that contains the handles for the
* tracking of solution values.
*
* Depending on the parameters in the INDATxD.F, solutions in different
* vertices are tracked: either they are saved at the end of a stationary
* simulation, or they are written out over time in the nonstationary
* simulation. The structure in this file defines the handles where to
* write the tracked values to.
*
* The parameter block is defined as a structure and realized as 
* integer/double precision arrays:
*   IPARAM = array [1..SZTRKI] of integer 
*          = TNSDefPostprocessingIParams
*
* Elements of the structure can be accessed with the help of the offset
* constants in this file, e.g.
*
*   HBFB1 = IPARAM(OHBFB1)
*
* The structures are defined as follows:
*
* TTrackingIParams = record
*
*   HBFB1  : File handle for output of body forces on real boundary.
*            Horizontal force.
*            =0, if tracked body forces should not be written into
*                a tracking file.
*   HBFB2  : File handle for output of body forces on real boundary.
*            Vertical force.
*            =0, if tracked body forces should not be written into
*                a tracking file.
*   HBFV1  : File handle for output of body forces on fict. boundary.
*            Horizontal force.
*            =0, if tracked body forces should not be written into
*                a tracking file.
*   HBFV2  : File handle for output of body forces on fict. boundary.
*            Vertical force.
*            =0, if tracked body forces should not be written into
*                a tracking file.
*   HBFL1  : File handle for output of body forces on fict. boundary
*            by integration on reconstructed interface.
*            Horizontal force.
*            =0, if tracked body forces should not be written into
*                a tracking file.
*   HBFL2  : File handle for output of body forces on fict. boundary
*            by integration on reconstructed interface.
*            Vertical force.
*            =0, if tracked body forces should not be written into
*                a tracking file.
*   HVEL1  : File handle for output of velocity values in a couple
*            of user-defined points.
*            X-Velocity.
*            =0, if values should not be written into a tracking file.
*   HVEL2  : File handle for output of velocity values in a couple
*            of user-defined points.
*            Y-Velocity.
*            =0, if values should not be written into a tracking file.
*   HPRES  : File handle for output of pressure values in a couple
*            of user-defined points.
*            =0, if values should not be written into a tracking file.
*   HSTRM  : File handle for output of streamfunction values in a couple
*            of user-defined points.
*            =0, if values should not be written into a tracking file.
*   HIPRS  : File handle for output of mean integral pressure on user
*            defined segments of real boundary.
*            =0, if values should not be written into a tracking file.
*
* end
* 
************************************************************************

      INTEGER OHBFB1, OHBFB2, OHBFV1,OHBFV2,OHBFL1,OHBFL2,OHVEL1,OHVEL2
      INTEGER OHPRES, OHSTRM, OHIPRS

      PARAMETER (OHBFB1  =  1)
      PARAMETER (OHBFB2  =  2)
      PARAMETER (OHBFV1  =  3)
      PARAMETER (OHBFV2  =  4)
      PARAMETER (OHBFL1  =  5)
      PARAMETER (OHBFL2  =  6)
      PARAMETER (OHVEL1  =  7)
      PARAMETER (OHVEL2  =  8)
      PARAMETER (OHPRES  =  9)
      PARAMETER (OHSTRM  = 10)
      PARAMETER (OHIPRS  = 11)

      INTEGER SZTRKI
      PARAMETER (SZTRKI  = 16)
      
