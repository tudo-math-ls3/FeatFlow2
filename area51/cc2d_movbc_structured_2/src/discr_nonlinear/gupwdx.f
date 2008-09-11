************************************************************************
* This file contains upwinding-based routines to set up the
* nonlinear convective part of the system matrix in the Navier
* Stokes equations.
************************************************************************

************************************************************************
* Standard upwinding
*
* Purpose: -  Adds the upwind-part on matrix block A after
*             it was initialized by the Stokes matrix
*          -  The input vector (U1 ,U2 ) is the old velocity field
*          -  The input vectors (UiL1,UiL2) in linear combination  
*             Ail are the transport direction
*
* In:
*   A,
*   KCOLA,
*   LLDA   - array [1..*} of double/integer
*            Structure arrays of system matrix, 
*            maybe initialised with linear parts.  
*   KVERT,
*   KMID,
*   DCORVG - usual geometry information
*
*   U1L1,
*   U1L2   - array [1..NU] of double
*            main velocity field used for assembling of 
*            the nonlinearity. Can be undefined if A1L=0.
*   U2L1,
*   U2L2   - array [1..NU] of double
*            secondary velocity field, used for the assembling
*            of the nonlinearity. Can be undefined if A2L=0.
*   A1L    - double; weighting factor for U1L1/U1L2
*   A2L    - double; weighting factor for U2L1/U2L2
*
*   IDEF   - Controls the behaviour of this routine.
*            =0: modify system matrix, add nonlinearity.
*                Defect vectors D1,D2 and velocity vectors U1,U2
*                can be undefined.
*            =1: modify both, system matrix and defect vector
*            =2: modify defect vectors, include nonlinearity.
*                A can be undefined (but not KCOLA,KLDA!)
*
*   U1,
*   U2     - array [1..NU] of double
*            Solution vector for modifying the defect vector.
*            Only if IDEF=1,2, otherwise it can be undefined.
*   D1,
*   D2     - array [1..NU] of double
*            Defect vector, modified by U1,U2.
*            Only if IDEF=1,2, otherwise it can be undefined.
*   
*   UPSAM  - control parameter.
*            -1=simple upwind
*            >=0: Samarskji upwind
*   RE     - 1/nu = viscosity
*   THSTEP - Current Step-size of the Theta-scheme.
*
* Out:
*   A      - system matrix;
*            the nonlinearity is added to that matrix
*   D1,
*   D2     - Modified defect vector; only if IDEF=1,2.
*
* Remarks:
*  
* 1.) In a typical call of the upwinding, the caller can use:
*     A1L = 1, U1L1/U1L2 = velocity field
*     A2L = 0, U2L1/U2L2 = undefined
*   So the upwinding scheme only uses one velocity field.
*   Such a call e.g. adds the integral
*                ( U1Lx*grad(.) , v )_Omega
*   to the system matrix.
*
*  2.) In case that there are two velocity fields representing
*   the solution (may happen in a nonstationary simulation where
*   U1L1/U1L2 represents the solution in the current and U2L1/U2L2
*   that of the previous time step), A1L/A2L defines how these both
*   velocity vectors should be weighted to compute the actual
*   velocity field for the assembling:
*                U_act = A1L*U1Lx + A2L*U2Lx
*   This is e.g. used for the linear extrapolation technique to
*   reconstruct a velocity from two previous time steps...
*
*  3.) In the nonlinear iteration, as a right hand side there arises
*   a defect vector D, which linear part can easily being assembled.
*   However, there is a nonlinearity to be included into that vector,
*   too. By setting IDEF=1,2, this routine incorporates the nonlinearity
*   into that vector, using the formula
*
*             D = D - THSTEP * UUx * grad (Ux)
*   
************************************************************************
      SUBROUTINE GUPWD (U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,
     *                  A,KCOLA,KLDA,KVERT,KMID,DCORVG,NEL,NVT,IDEF,
     *                  UPSAM,RE,THSTEP)

      IMPLICIT NONE
      
C standard COMMON blocks

      INCLUDE 'cout.inc'
      
      INCLUDE 'cbasictria.inc'
      
C parameters
      
      DOUBLE PRECISION U1L1(*),U1L2(*),U2L1(*),U2L2(*),U1(*),U2(*)
      DOUBLE PRECISION D1(*),D2(*)
      DOUBLE PRECISION A(*)
      INTEGER KCOLA(*),KLDA(*),IDEF,NEL,NVT
      DOUBLE PRECISION A1L,A2L
      
      DOUBLE PRECISION UPSAM,RE,THSTEP
      
C *** Usual data for mesh management
      INTEGER KVERT(NNVE,*),KMID(NNVE,*)
      DOUBLE PRECISION DCORVG(2,*)

C local variables

C *** Local arrays for informations about one element

      INTEGER IMID(4),ISTORE(4,4)
      DOUBLE PRECISION FLUX(4),UU1(4),UU2(4),XV(4),YV(4)
      DOUBLE PRECISION ELMA(4,4)

      DOUBLE PRECISION XE, YE, XN, YN, G1, G2, UPSRE, ELMH
      DOUBLE PRECISION DL0, DL2, H00, H22, FL0, FL2
      INTEGER I,II, IEL, IV
      INTEGER IM1, IA1, IA2, J, JJ, IA, IM0, IM2
      
C implicit function definitions

*********************************************************************
*   Weighted Samarski upwind
*
*   This implementation follows the documentation of
*   [F. Schieweck, Parallele Lösung der stationären inkompressiblen
*    Navier-Stokes Gleichungen, Habilitation, Fakultät für
*    Mathematik, Otto-von-Guericke-Universität Magdeburg]
*
*********************************************************************

      DOUBLE PRECISION X, PHIP, PHIM

      PHIP(X)=(0.5D0+X)/(1D0+X)
      PHIM(X)=    0.5D0/(1D0-X)
      
C *******************************************************************
C     What we want to discretize here is the term
C
C         n(z,u,v) = ( (z*grad (.))u , v )_Omega
C
C     Let's assume we have two elements next to each other:
C
C       X---------------X---------------X
C       |            /  |               |
C       |          /    |               |
C       |        /      |               |
C       |       X       Gl              I
C       |        \      |               |
C       |         Glk   |               |
C       |            \  |               |
C       X------Gk-------X---------------X
C
C     The edges Gl and Gk enclose a diagonal edge Gklof the above
C     triangle. The operator can now be rewritten by decomposition
C     onto the elements as
C
C       n(z,u,v) ~= sum_l sum_k int_Gkl (z*n_lk) (u-u(Bl)) v(Bl) dGamma
C
C     with Bl and Bk being the midpoints of the edges and n_lk being the
C     outer normal vector of the edge Glk. 
C
C       X---------------X              X---------------X
C       |            /  |              |            /  |
C       |          /    |              |          /    |
C       |        /      |              |        /      |
C       |       X       X Bl           |       X       u_l
C       |        \      |              |        \      |
C       |          \    |              |        u_upw  |
C       |            \  |              |            \  |
C       X-------X-------X              X------u_k------X
C               Bk
C
C     The integral at the end of this term is replaced by 1x-Gauss 
C     rule, thus u can be replaced by an approximation u_upw on the
C     edge Glk - which is calculated with the help of the velocities 
C     in the neighborhood u_l and u_k.
C
C     The main task in Upwinding is thus to calc u_upw - a re-
C     constructed velocity, based on u_l and u_k !
C     (Remember: With nonconforming, rotated bilinear elements,
C      we have the velocity in the midpoints/edges - but not
C      in the midpoints of the diagonals of those triangles.
C      But there we need it because of an integration along this
C      edge with simple Gauss rule.)
C
C     What's the approach here? As said, u_upw is reconstructed
C     from u_1 and u_2 by a simple mean formula:
C
C          u_upw = Lambda u_l  + (1-Lambda) u_k
C
C     What is Lambda? 0<Lambda<1 is chosen depending on the flow
C     crossing the diagonal Glk. More precisely, it's chosen
C     depending on the flow:
C       Flow direction       lambda        u_upw
C          Bk -> Bl            ~0          ~u_k
C          Bl -> Bk            ~1          ~u_l
C          equal              ~0.5    ~mean between u_k and u_l
C
C     The "flow" is described by z. The "flow through Glk" or 
C     "flux" is described by the line integral
C
C             t = 1/nu int_Glk (z*nlk) dGamma
C
C     (again with nlk being the normal vector of the edge Glk).
C
C     The parameter lambda is now chosen as:
C
C        lambda =    1 - 1/(2+2theta*t))  ,  t >= 0
C               =    1/(2+2theta*t))      ,  t < 0
C
C     with theta being the UPSAM-parameter from the DAT-file.
C     So UPSAM controls the weighting of the two neighboring
C     velocities when reconstructing the velocity in the
C     midpoint of the diagonal.
C
C     For theta=UPSAM = 0, we have central difference.
C     For theta=UPSAM -> infinity, we obtain the simple upwind
C     (lambda=0 for t<0, lambda=1 for t>=0).
C
C *******************************************************************

C     Loop over all elements in the current grid. The data is
C     elementwise collected and added to the matrix.

      DO IEL=1,NEL

C       XE,YE will be coordinates of the center of the element

        XE=0.D0
        YE=0.D0
      
C       Loop over all 4 U-nodes.
C       Calculate IMID,XV,YV,XE,YE,UU1,UU2

        DO II=1,4
        
          I=KMID(II,IEL)-NVT
          
C         Store the number of the edge in IMID:
          
          IMID(II)=I
          
C         Store the coordinates of the corner vertices of that
C         element in XV/YV:
          
          IV=KVERT(II,IEL)
          
          XV(II)=DCORVG(1,IV)
          YV(II)=DCORVG(2,IV)

C         Sum up the coordinates if the element - will later result
C         in the element midpoint:
          
          XE=XE+XV(II)
          YE=YE+YV(II)
          
C         Compute the actual velocity in the edge II (following the
C         node II) by a weighted mean of the both velocity vectors
C         U1Lx and U2Lx. This allowes e.g. to reconstruct a velocity
C         vector by linear extrapolation of two previous time steps.
          
          UU1(II)=A1L*U1L1(I)+A2L*U2L1(I)
          UU2(II)=A1L*U1L2(I)+A2L*U2L2(I)
          
        END DO

C       Divide XE/YE by 4 - so XE/YE receives the coordinate
C       of the element midpoint

        XE = 0.25D0*XE
        YE = 0.25D0*YE
        
C       After this procedure we have the following variable setting:
C
C   (XV(4),YV(4))               (XV(3),YV(3))
C               X----IMID(3)----X
C               |               | 
C               |               | 
C               |               | 
C         IMID(4)       X       IMID(2)
C               |    (XE,YE)    | 
C               |               | 
C               |               |  
C               X----IMID(1)----X
C   (XV(1),YV(1))               (XV(2),YV(2))
C
C       UU1/UU2 contains the velocity along the edges following
C       the four corner points.
C        
C       Loop over all 4 U-nodes.
C       Calculate FLUX(.), ISTORE(.,.), ELMA(.,.)

        DO II=1,4

C *** Setting II-1 modulo 4 on IM1
C
C         II1 corresponds to the current node in question.
C         IM1 receives the predecessor of II1 in an anticlockwise
C         sense:

          IM1=II-1
          IF (IM1.LT.1) IM1=4
      
C *** Calculation of the flux  FLUX(II)
C
C                    /    |                        |
C                  /      |                        |
C              XE/YE      |             XE/YE      u_l
C                  \      |                 \      |
C                   /\    |                  GP    |
C                  /   \  |                     \  |
C         IM------/-------II       IM-----u_k------II
C                / n
C               v

C         From the mitpoint XE/YE and the current corner II,
C         calculate the outer normal vector n of the edge Glk 
C         that connects II with the midpoint:

          XN=-YV(II)+YE
          YN= XV(II)-XE
          
C         Calculate the (scaled) flux
C
C           t = int_Glk (z*nlk) dGamma
C             ~= nlk * u(GP)
C              = 1/2 * nlk * (u_l+u_k)
C
C         approximating the integral with 1-point Gauss in the
C         Gauss-point (=midpoint) of the edge. Save t in FLUX(II).
C         So FLUX(II) saves the flux along the edge that is pointing
C         from II to the element midpoint.
          
          G1=0.5D0*(UU1(IM1)+UU1(II))
          G2=0.5D0*(UU2(IM1)+UU2(II))
          FLUX(II)=XN*G1+YN*G2
      
C *** Determine the indices ISTORE(II,JJ) to store the element matrix
C     entry ELMA(II,JJ) on array A

C         Go into line I of the matrix, corresponding to the current
C         edge IMID(II). IA1 and IA2 saves the start- and end-index
C         of this row of the matrix.

          I=IMID(II)
          IA1=KLDA(I)
          IA2=KLDA(I+1)-1
          
C         Loop over the edges of the element
          
          DO JJ=1,4

C           In the current row, search for column J=IMID(JJ).
C           We know from FE-Theory, that our current node I interacts
C           only with all the nodes on the edges of the current element,
C           so only there can be an additional value to be included
C           into the system matrix.
C
C             |---J---|
C             |       |
C             J       J
C             |       |
C             |---I---|
C             |       |
C
C           So search in the current row for the matrix entry that
C           corresponds to the common support of DOF I and J.

            J=IMID(JJ)

            DO IA=IA1,IA2
              IF (KCOLA(IA).EQ.J)  GOTO 121
            END DO

C *** Error case

            WRITE(MTERM,*) 'ERROR in GUPWD: entry index IA not found'
            RETURN

121         CONTINUE

C           Save the matrix index in ISTORE(II,JJ) so we can find it
C           later.

            ISTORE(II,JJ)=IA      

C           Initialize ELMA(.,.) to 0. ELMA will assemble the "local"
C           matrix, i.e. the values that are later incorporated into
C           the system matrix at the positions stored in ISTORE.

            ELMA(II,JJ)=0D0

          END DO ! JJ
          
        END DO ! II
        
C       What have we calculated up to here? Let's collect...
C
C       FLUX   - The flux along the edges of the triangles
C       ISTORE - The indices in the matrix A of the entries that
C                are affected by integration on the current element IEL.
C       ELMA   - Local element matrix - initialized with 0 up to now.
C
C       So the next step is to assemble the local element matrix,
C       which is to be incorporated into the system matrix later.
C
C       Loop over the nodes to calculate ELMA:

        DO II=1,4
      
C         Set IM1=predecessor of II, IM2=successor of II,
C         in counterclockwise sense.

          IM0=II
          IM1=II-1
          IF (IM1.LT.1) IM1=4
          IM2=II+1
          IF (IM2.GT.4) IM2=1

C         We interpret II here as the local number of the edge
C         (rather than the corner), following the corner II.
C         The current triangle-edge GAMMA is the edge between
C         corner II and the midpoint.
C
C          +------IM2------IM2
C          |             / |
C          |          FL2  |
C          |         /     |
C          |       X       II=IM0
C          |         \     |
C          |          FL0  |
C          |             \ |
C        IM1------IM1------II=IM0

C *** Calculate the part corresponding to GAMMA(IM0) and GAMMA(IM2)

C         Get the flux of the edges IM1->center, IM2->center and
C         save it to FL0,FL2:

          FL0=FLUX(IM0)
          FL2=FLUX(IM2)

C         Now choose the Upwind scheme depending on UPSAM:
C         UPSAM>0: Samarskji-Upwind
C         UPSAM<0: Simple upwind

          IF (UPSAM.GE.0) THEN

C           The user wants Samarskji-Upwind.
C           Weight the UPSAM-parameter by 1/nu.
C           Remember: In the previous calculation of the line-integral
C           ro calculate t, we didn't incorporate 1/nu - this is 
C           repaired here:

            UPSRE=UPSAM*RE

C           Analyze the two fluxes on the edges of the triangle.
C           Calculate the lambda-value by Samarskji-upwind,
C           depending on whether the flux goes "in" ou "out" of the
C           triangle. DL0 saves the lambda of the triangle-edge
C           corresponding to II, DL2 that of the triangle edge
C           corresponding to  IM2.

            IF (FL0.GE.0D0) THEN
              DL0=PHIP(UPSRE*FL0)
            ELSE
              DL0=PHIM(UPSRE*FL0)
            ENDIF

            IF (FL2.GE.0D0) THEN
              DL2=PHIP(UPSRE*FL2)
            ELSE
              DL2=PHIM(UPSRE*FL2)
            ENDIF

          ELSE

C           Simple Upwinding scheme.
C           The corresponding lambda (wighting factor of the
C           "adjacent velocities) is either 0 or 1, depending on
C           whether the flux goes "into" or "out of" the triangle
C           on that edge.

            DL0=0D0
            DL2=0D0
            IF (FL0.GE.0.D0) DL0=1D0
            IF (FL2.GE.0.D0) DL2=1D0

          ENDIF

C         Calculate the local element matrix with the integrals
C         being evaluated by 1-point Gauss rule in the midpoint
C         of the edges of the triangle. 
C
C         In fact, this calculates the local contribution of 
C            UUx * grad ( . )

          H00=DL0*FL0
          H22=(1D0-DL2)*FL2
          ELMA(IM0,IM0) = H00-H22
          ELMA(IM0,IM2) =     H22
          ELMA(IM0,IM1) =-H00
          
        END DO ! II
      
C       We are nearly done. Now we only have to incorporate
C       the local element matrix ELMA into the global matrix.
C       This is simple: Grab the index to be modified from ISTORE
C       and add the corresponding ELMA-value to that matrix entry.
C
C       Another matter of concern is the update of a defect vector,
C       which arises during a nonlinear iteration. If IDEF=1,2,
C       we modify such a defect vector D in the form
C
C           D = D - THSTEP * UUx * grad (Ux)
C
C       i.e. using a solution vector U (independent of the velocity
C       fields), we subtract the nonlinearity from D. 
      
        DO II=1,4
        
          DO JJ=1,4
          
C           Weight the local matrix by THSTEP. Remember, in the 
C           nonlinear iteration we have to incorporate the term
C
C               THETA*K*u*grad(u)
C
C           which is realized by a multiplication of the u*grad(.)-
C           matrix by THETA*K = THSTEP here!
          
            ELMH=THSTEP*ELMA(II,JJ)

C           If IDEF=0,1, modify the system matrix:

            IF (IDEF.LT.2) THEN
              IA   =ISTORE(II,JJ)
              A(IA)=A(IA)+ELMH
            ENDIF

C           If IDEF=1,2, modify the defect vectors.
C
C           Multiply the velocity Ux by the local matrix, resulting in
C           the term UUx*grad(Ux) - and subtract that from the defect
C           vector.

            IF (IDEF.GT.0) THEN 
              D1(IMID(II))=D1(IMID(II))-ELMH*U1(IMID(JJ))
              D2(IMID(II))=D2(IMID(II))-ELMH*U2(IMID(JJ))
            ENDIF 

          END DO ! JJ
          
        END DO ! II
      
C-----------------------------------------------------------------------

      END DO ! IEL

      END
      
************************************************************************
*  Extended upwinding with 1st order ALE method 
*
*  Purpose: -  Adds the upwind-part on matrix block A after
*              it was initialized by the Stokes matrix
*           -  The input vector (U1 ,U2 ) is the old velocity field
*           -  The input vectors (UiL1,UiL2) in linear combination  
*              Ail are the transport direction
*           -  Adds a mesh velocity field to the nonlinearity if desired
*
* In:
*   A,
*   KCOLA,
*   LLDA   - array [1..*] of double/integer
*            Structure arrays of system matrix, 
*            maybe initialised with linear parts.  
*   TRIA   - array [1..SZTRIA] of integer
*            Triangulation structure about the current mesh
*   KVERT,
*   KMID,
*   DCORVG - Usual geometry information. These must correspond
*            to the information in the TRIA-structure, but are
*            passed separately for faster access.
*   U1L1,
*   U1L2   - array [1..NU] of double
*            main velocity field used for assembling of 
*            the nonlinearity. Can be undefined if A1L=0.
*   U2L1,
*   U2L2   - array [1..NU] of double
*            secondary velocity field, used for the assembling
*            of the nonlinearity. Can be undefined if A2L=0.
*   A1L    - double; weighting factor for U1L1/U1L2
*   A2L    - double; weighting factor for U2L1/U2L2
*
*   IDEF   - Controls the behaviour of this routine.
*            =0: modify system matrix, add nonlinearity.
*                Defect vectors D1,D2 and velocity vectors U1,U2
*                can be undefined.
*            =1: modify both, system matrix and defect vector
*            =2: modify defect vectors, include nonlinearity.
*                A can be undefined (but not KCOLA,KLDA!)
*
*   U1,
*   U2     - array [1..NU] of double
*            Solution vector for modifying the defect vector.
*            Only if IDEF=1,2, otherwise it can be undefined.
*   D1,
*   D2     - array [1..NU] of double
*            Defect vector, modified by U1,U2.
*            Only if IDEF=1,2, otherwise it can be undefined.
*   
*   IALE   - Controls whether or not the ALE method should be
*            used.
*            =0: no ALE, standard UPWIND. DMVALE can be undefined.
*            =1: use ALE. DMVALE must be defined.
*   DMVALE - array [1..2,1..NVT] of double
*            Mesh velocity field to include when setting up the
*            nonlinearity. (DMVALE(1,.),DMVALE(2,.)) represents
*            the X- and Y-derivative, respectively, of the speed
*            of the 1..NVT vertices in the mest.
*
*   UPSAM  - control parameter.
*            -1=simple upwind
*            >=0: Samarskji upwind
*   RE     - 1/nu = viscosity
*   THSTEP - Current Step-size of the Theta-scheme.
*
* Out:
*   A      - system matrix;
*            the nonlinearity is added to that matrix
*   D1,
*   D2     - Modified defect vector; only if IDEF=1,2.
*
* Remarks:
*  
* 1.) In a typical call of the upwinding, the caller can use:
*     A1L = 1, U1L1/U1L2 = velocity field
*     A2L = 0, U2L1/U2L2 = undefined
*   So the upwinding scheme only uses one velocity field.
*   Such a call e.g. adds the integral
*                ( U1Lx*grad(.) , v )_Omega
*   to the system matrix.
*
*  2.) In case that there are two velocity fields representing
*   the solution (may happen in a nonstationary simulation where
*   U1L1/U1L2 represents the solution in the current and U2L1/U2L2
*   that of the previous time step), A1L/A2L defines how these both
*   velocity vectors should be weighted to compute the actual
*   velocity field for the assembling:
*                U_act = A1L*U1Lx + A2L*U2Lx
*   This is e.g. used for the linear extrapolation technique to
*   reconstruct a velocity from two previous time steps...
*
*  3.) In the nonlinear iteration, as a right hand side there arises
*   a defect vector D, which linear part can easily being assembled.
*   However, there is a nonlinearity to be included into that vector,
*   too. By setting IDEF=1,2, this routine incorporates the nonlinearity
*   into that vector, using the formula
*
*             D = D - THSTEP * UUx * grad (Ux)
*   
*  4.) This routine makes no use of COMMON-blocks about the
*   triangulation; all information must be passed in the TRIA structure.
*   The triangulation information KVERT,KMID,DCORVG must correspond
*   to the arrays identified by the handles in TRIA, but are passed
*   separately for faster access.
*
*  5.) If IALE=1, a mesh velocity field is added to the nonlineareity
*   according to the formula  "U * grad (U-DMVALE)".
*   For IALE=0, the simple nonlinearity "U * grad (U)" is used.
************************************************************************
      SUBROUTINE GUPWDX (U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,
     *                   A,KCOLA,KLDA,TRIA,KVERT,KMID,DCORVG,IDEF,
     *                   UPSAM,RE,THSTEP, IALE, DMVALE)

      IMPLICIT NONE
      
C standard COMMON blocks

      INCLUDE 'cout.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
C parameters
      
      DOUBLE PRECISION U1L1(*),U1L2(*),U2L1(*),U2L2(*),U1(*),U2(*)
      DOUBLE PRECISION D1(*),D2(*)
      DOUBLE PRECISION A(*)
      INTEGER KCOLA(*),KLDA(*),IDEF
      DOUBLE PRECISION A1L,A2L
      
      DOUBLE PRECISION UPSAM,RE,THSTEP
      
      INTEGER IALE
      DOUBLE PRECISION DMVALE(2,*)

C *** Usual data for mesh management
      INTEGER TRIA(SZTRIA)
      INTEGER KVERT(NNVE,*),KMID(NNVE,*)
      DOUBLE PRECISION DCORVG(2,*)

C local variables

C *** Local arrays for informations about one element

      INTEGER IMID(4),ISTORE(4,4)
      DOUBLE PRECISION FLUX(4),UU1(4),UU2(4),XV(4),YV(4)
      DOUBLE PRECISION ELMA(4,4)

      DOUBLE PRECISION XE, YE, XN, YN, G1, G2, UPSRE, ELMH
      DOUBLE PRECISION DL0, DL2, H00, H22, FL0, FL2
      INTEGER I,II, IEL, IV
      INTEGER IM1, IA1, IA2, J, JJ, IA, IM0, IM2
      
      DOUBLE PRECISION UALE1,UALE2
      INTEGER IVT1,IVT2
      
C implicit function definitions

*********************************************************************
*   Weighted Samarski upwind
*
*   This implementation follows the documentation of
*   [F. Schieweck, Parallele Lösung der stationären inkompressiblen
*    Navier-Stokes Gleichungen, Habilitation, Fakultät für
*    Mathematik, Otto-von-Guericke-Universität Magdeburg]
*
*********************************************************************

      DOUBLE PRECISION X, PHIP, PHIM

      PHIP(X)=(0.5D0+X)/(1D0+X)
      PHIM(X)=    0.5D0/(1D0-X)
      
C *******************************************************************
C     What we want to discretize here is the term
C
C         n(z,u,v) = ( (z*grad (.))u , v )_Omega
C
C     Let's assume we have two elements next to each other:
C
C       X---------------X---------------X
C       |            /  |               |
C       |          /    |               |
C       |        /      |               |
C       |       X       Gl              I
C       |        \      |               |
C       |         Glk   |               |
C       |            \  |               |
C       X------Gk-------X---------------X
C
C     The edges Gl and Gk enclose a diagonal edge Gklof the above
C     triangle. The operator can now be rewritten by decomposition
C     onto the elements as
C
C       n(z,u,v) ~= sum_l sum_k int_Gkl (z*n_lk) (u-u(Bl)) v(Bl) dGamma
C
C     with Bl and Bk being the midpoints of the edges and n_lk being the
C     outer normal vector of the edge Glk. 
C
C       X---------------X              X---------------X
C       |            /  |              |            /  |
C       |          /    |              |          /    |
C       |        /      |              |        /      |
C       |       X       X Bl           |       X       u_l
C       |        \      |              |        \      |
C       |          \    |              |        u_upw  |
C       |            \  |              |            \  |
C       X-------X-------X              X------u_k------X
C               Bk
C
C     The integral at the end of this term is replaced by 1x-Gauss 
C     rule, thus u can be replaced by an approximation u_upw on the
C     edge Glk - which is calculated with the help of the velocities 
C     in the neighborhood u_l and u_k.
C
C     The main task in Upwinding is thus to calc u_upw - a re-
C     constructed velocity, based on u_l and u_k !
C     (Remember: With nonconforming, rotated bilinear elements,
C      we have the velocity in the midpoints/edges - but not
C      in the midpoints of the diagonals of those triangles.
C      But there we need it because of an integration along this
C      edge with simple Gauss rule.)
C
C     What's the approach here? As said, u_upw is reconstructed
C     from u_1 and u_2 by a simple mean formula:
C
C          u_upw = Lambda u_l  + (1-Lambda) u_k
C
C     What is Lambda? 0<Lambda<1 is chosen depending on the flow
C     crossing the diagonal Glk. More precisely, it's chosen
C     depending on the flow:
C       Flow direction       lambda        u_upw
C          Bk -> Bl            ~0          ~u_k
C          Bl -> Bk            ~1          ~u_l
C          equal              ~0.5    ~mean between u_k and u_l
C
C     The "flow" is described by z. The "flow through Glk" or 
C     "flux" is described by the line integral
C
C             t = 1/nu int_Glk (z*nlk) dGamma
C
C     (again with nlk being the normal vector of the edge Glk).
C
C     The parameter lambda is now chosen as:
C
C        lambda =    1 - 1/(2+2theta*t))  ,  t >= 0
C               =    1/(2+2theta*t))      ,  t < 0
C
C     with theta being the UPSAM-parameter from the DAT-file.
C     So UPSAM controls the weighting of the two neighboring
C     velocities when reconstructing the velocity in the
C     midpoint of the diagonal.
C
C     For theta=UPSAM = 0, we have central difference.
C     For theta=UPSAM -> infinity, we obtain the simple upwind
C     (lambda=0 for t<0, lambda=1 for t>=0).
C
C *******************************************************************

C     Loop over all elements in the current grid. The data is
C     elementwise collected and added to the matrix.

      DO IEL=1,TRIA(ONEL)

C       XE,YE will be coordinates of the center of the element

        XE=0.D0
        YE=0.D0
      
C       Loop over all 4 U-nodes.
C       Calculate IMID,XV,YV,XE,YE,UU1,UU2

        DO II=1,4
        
          I=KMID(II,IEL)-TRIA(ONVT)
          
C         Store the number of the edge in IMID:
          
          IMID(II)=I
          
C         Store the coordinates of the corner vertices of that
C         element in XV/YV:
          
          IV=KVERT(II,IEL)
          
          XV(II)=DCORVG(1,IV)
          YV(II)=DCORVG(2,IV)

C         Sum up the coordinates if the element - will later result
C         in the element midpoint:
          
          XE=XE+XV(II)
          YE=YE+YV(II)
          
C         Now we want to compute the velocity on the edge II
C         (following the node II). Here's the point where we include
C         our ALE-stuff:

          UALE1 = 0D0
          UALE2 = 0D0

          IF (IALE.GT.0) THEN
          
C           In the ALE-equation, the nonlinear part is slightly modified:
C
C                     u * grad(u)    --->    (u-v) * grad(u)
C
C           with v being the "mesh velocity field", i.e. an approximation
C           to the velocity of the vertices in the mesh. Our mesh
C           velocity field DMVALE is given in the nodes. To compute
C           it in the midpoints of the edges, we take the mean and
C           save the result in UALE1/2 for the X- and Y-coordinate.
C
C           At first, compute the start/endpoint of the edge II:

            IVT1 = KVERT(II,IEL)
            IVT2 = KVERT(MOD(II,4)+1,IEL)

C           And then calculate the mean velocity field:

            UALE1 = 0.5D0 * (DMVALE (1,IVT1) + DMVALE(1,IVT2) )
            UALE2 = 0.5D0 * (DMVALE (2,IVT1) + DMVALE(2,IVT2) )
          
          END IF
          
C         Compute the actual velocity in the edge II (following the
C         node II) by a weighted mean of the both velocity vectors
C         U1Lx and U2Lx. This allowes e.g. to reconstruct a velocity
C         vector by linear extrapolation of two previous time steps.
C
C         Subtract the mesh-velocity on the current edge as described
C         above if we are using an ALE-formulation.
          
          UU1(II) = A1L*U1L1(I)+A2L*U2L1(I) - UALE1
          UU2(II) = A1L*U1L2(I)+A2L*U2L2(I) - UALE2
          
        END DO

C       Divide XE/YE by 4 - so XE/YE receives the coordinate
C       of the element midpoint

        XE = 0.25D0*XE
        YE = 0.25D0*YE
        
C       After this procedure we have the following variable setting:
C
C   (XV(4),YV(4))               (XV(3),YV(3))
C               X----IMID(3)----X
C               |               | 
C               |               | 
C               |               | 
C         IMID(4)       X       IMID(2)
C               |    (XE,YE)    | 
C               |               | 
C               |               |  
C               X----IMID(1)----X
C   (XV(1),YV(1))               (XV(2),YV(2))
C
C       UU1/UU2 contains the velocity along the edges following
C       the four corner points.
C        
C       Loop over all 4 U-nodes.
C       Calculate FLUX(.), ISTORE(.,.), ELMA(.,.)

        DO II=1,4

C *** Setting II-1 modulo 4 on IM1
C
C         II1 corresponds to the current node in question.
C         IM1 receives the predecessor of II1 in an anticlockwise
C         sense:

          IM1=II-1
          IF (IM1.LT.1) IM1=4
      
C *** Calculation of the flux  FLUX(II)
C
C                    /    |                        |
C                  /      |                        |
C              XE/YE      |             XE/YE      u_l
C                  \      |                 \      |
C                   /\    |                  GP    |
C                  /   \  |                     \  |
C         IM------/-------II       IM-----u_k------II
C                / n
C               v

C         From the mitpoint XE/YE and the current corner II,
C         calculate the outer normal vector n of the edge Glk 
C         that connects II with the midpoint:

          XN=-YV(II)+YE
          YN= XV(II)-XE
          
C         Calculate the (scaled) flux
C
C           t = int_Glk (z*nlk) dGamma
C             ~= nlk * u(GP)
C              = 1/2 * nlk * (u_l+u_k)
C
C         approximating the integral with 1-point Gauss in the
C         Gauss-point (=midpoint) of the edge. Save t in FLUX(II).
C         So FLUX(II) saves the flux along the edge that is pointing
C         from II to the element midpoint.
          
          G1=0.5D0*(UU1(IM1)+UU1(II))
          G2=0.5D0*(UU2(IM1)+UU2(II))
          FLUX(II)=XN*G1+YN*G2
      
C *** Determine the indices ISTORE(II,JJ) to store the element matrix
C     entry ELMA(II,JJ) on array A

C         Go into line I of the matrix, corresponding to the current
C         edge IMID(II). IA1 and IA2 saves the start- and end-index
C         of this row of the matrix.

          I=IMID(II)
          IA1=KLDA(I)
          IA2=KLDA(I+1)-1
          
C         Loop over the edges of the element
          
          DO JJ=1,4

C           In the current row, search for column J=IMID(JJ).
C           We know from FE-Theory, that our current node I interacts
C           only with all the nodes on the edges of the current element,
C           so only there can be an additional value to be included
C           into the system matrix.
C
C             |---J---|
C             |       |
C             J       J
C             |       |
C             |---I---|
C             |       |
C
C           So search in the current row for the matrix entry that
C           corresponds to the common support of DOF I and J.

            J=IMID(JJ)

            DO IA=IA1,IA2
              IF (KCOLA(IA).EQ.J)  GOTO 121
            END DO

C *** Error case

            WRITE(MTERM,*) 'ERROR in GUPWD: entry index IA not found'
            RETURN

121         CONTINUE

C           Save the matrix index in ISTORE(II,JJ) so we can find it
C           later.

            ISTORE(II,JJ)=IA      

C           Initialize ELMA(.,.) to 0. ELMA will assemble the "local"
C           matrix, i.e. the values that are later incorporated into
C           the system matrix at the positions stored in ISTORE.

            ELMA(II,JJ)=0D0

          END DO ! JJ
          
        END DO ! II
        
C       What have we calculated up to here? Let's collect...
C
C       FLUX   - The flux along the edges of the triangles
C       ISTORE - The indices in the matrix A of the entries that
C                are affected by integration on the current element IEL.
C       ELMA   - Local element matrix - initialized with 0 up to now.
C
C       So the next step is to assemble the local element matrix,
C       which is to be incorporated into the system matrix later.
C
C       Loop over the nodes to calculate ELMA:

        DO II=1,4
      
C         Set IM1=predecessor of II, IM2=successor of II,
C         in counterclockwise sense.

          IM0=II
          IM1=II-1
          IF (IM1.LT.1) IM1=4
          IM2=II+1
          IF (IM2.GT.4) IM2=1

C         We interpret II here as the local number of the edge
C         (rather than the corner), following the corner II.
C         The current triangle-edge GAMMA is the edge between
C         corner II and the midpoint.
C
C          +------IM2------IM2
C          |             / |
C          |          FL2  |
C          |         /     |
C          |       X       II=IM0
C          |         \     |
C          |          FL0  |
C          |             \ |
C        IM1------IM1------II=IM0

C *** Calculate the part corresponding to GAMMA(IM0) and GAMMA(IM2)

C         Get the flux of the edges IM1->center, IM2->center and
C         save it to FL0,FL2:

          FL0=FLUX(IM0)
          FL2=FLUX(IM2)

C         Now choose the Upwind scheme depending on UPSAM:
C         UPSAM>0: Samarskji-Upwind
C         UPSAM<0: Simple upwind

          IF (UPSAM.GE.0) THEN

C           The user wants Samarskji-Upwind.
C           Weight the UPSAM-parameter by 1/nu.
C           Remember: In the previous calculation of the line-integral
C           ro calculate t, we didn't incorporate 1/nu - this is 
C           repaired here:

            UPSRE=UPSAM*RE

C           Analyze the two fluxes on the edges of the triangle.
C           Calculate the lambda-value by Samarskji-upwind,
C           depending on whether the flux goes "in" ou "out" of the
C           triangle. DL0 saves the lambda of the triangle-edge
C           corresponding to II, DL2 that of the triangle edge
C           corresponding to  IM2.

            IF (FL0.GE.0D0) THEN
              DL0=PHIP(UPSRE*FL0)
            ELSE
              DL0=PHIM(UPSRE*FL0)
            ENDIF

            IF (FL2.GE.0D0) THEN
              DL2=PHIP(UPSRE*FL2)
            ELSE
              DL2=PHIM(UPSRE*FL2)
            ENDIF

          ELSE

C           Simple Upwinding scheme.
C           The corresponding lambda (wighting factor of the
C           "adjacent velocities) is either 0 or 1, depending on
C           whether the flux goes "into" or "out of" the triangle
C           on that edge.

            DL0=0D0
            DL2=0D0
            IF (FL0.GE.0.D0) DL0=1D0
            IF (FL2.GE.0.D0) DL2=1D0

          ENDIF

C         Calculate the local element matrix with the integrals
C         being evaluated by 1-point Gauss rule in the midpoint
C         of the edges of the triangle. 
C
C         In fact, this calculates the local contribution of 
C            UUx * grad ( . )

          H00=DL0*FL0
          H22=(1D0-DL2)*FL2
          ELMA(IM0,IM0) = H00-H22
          ELMA(IM0,IM2) =     H22
          ELMA(IM0,IM1) =-H00
          
        END DO ! II
      
C       We are nearly done. Now we only have to incorporate
C       the local element matrix ELMA into the global matrix.
C       This is simple: Grab the index to be modified from ISTORE
C       and add the corresponding ELMA-value to that matrix entry.
C
C       Another matter of concern is the update of a defect vector,
C       which arises during a nonlinear iteration. If IDEF=1,2,
C       we modify such a defect vector D in the form
C
C           D = D - THSTEP * UUx * grad (Ux)
C
C       i.e. using a solution vector U (independent of the velocity
C       fields), we subtract the nonlinearity from D. 
      
        DO II=1,4
        
          DO JJ=1,4
          
C           Weight the local matrix by THSTEP. Remember, in the 
C           nonlinear iteration we have to incorporate the term
C
C               THETA*K*u*grad(u)
C
C           which is realized by a multiplication of the u*grad(.)-
C           matrix by THETA*K = THSTEP here!
          
            ELMH=THSTEP*ELMA(II,JJ)

C           If IDEF=0,1, modify the system matrix:

            IF (IDEF.LT.2) THEN
              IA   =ISTORE(II,JJ)
              A(IA)=A(IA)+ELMH
            ENDIF

C           If IDEF=1,2, modify the defect vectors.
C
C           Multiply the velocity Ux by the local matrix, resulting in
C           the term UUx*grad(Ux) - and subtract that from the defect
C           vector.

            IF (IDEF.GT.0) THEN 
              D1(IMID(II))=D1(IMID(II))-ELMH*U1(IMID(JJ))
              D2(IMID(II))=D2(IMID(II))-ELMH*U2(IMID(JJ))
            ENDIF 

          END DO ! JJ
          
        END DO ! II
      
C-----------------------------------------------------------------------

      END DO ! IEL

      END
      
