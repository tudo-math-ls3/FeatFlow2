************************************************************************
* Prolongation/Restriction for the P2 element
************************************************************************

************************************************************************
* P2 prolongation
*
* Transfers a P2-solution from a coarser level to a finer level.
*
* In:
*   DUC    : array [1..NVTC+NMTC] of double
*            Coarse grid solution vector.
*   NVTC   : Number of vertices on the coarse grid
*   NMTC   : Number of edge midpoints on the coarse grid
*   NELC   : Number of elements on the coarse grid
*   DCORVC : array [1..2,1..NVTC] of double
*            Coordinates of all vertices on the coarse grid
*   KVERTC : array [1..NNVE,1..NVTC] of double
*            Vertices on each element on the coarse grid
*   KMIDC  : array [1..NNVE,1..NVTC] of double
*            Midpoints on each element on the coarse grid
*   KADJC  : array [1..NNVE,1..NVTC]
*            Numbers of adjacent elements to each element on the
*            coarse grid.
*   KVERTF : array [1..NNVE,1..NVTF] of double
*            Vertices on each element on the fine grid
*   KMIDF  : array [1..NNVE,1..NVTF] of double
*            Midpoints on each element on the fine grid
*
* Out:
*   DUF    : array [1..NVTF+NMTF] of double
*            Fine grid solution vector
************************************************************************

      SUBROUTINE MP002 (DUC,DUF,NVTC,NMTC,NELC,
     *                  KVERTC,KMIDC,KADJC,
     *                  KMIDF,KADJF)
      
      IMPLICIT NONE

      INTEGER NNVE      
      PARAMETER (NNVE=4)
      
C     parameters
      
      DOUBLE PRECISION DUC(*),DUF(*)
      INTEGER KVERTC(NNVE,*),KMIDC(NNVE,*),KADJC(NNVE,*),NVTC,NMTC,NELC
      INTEGER KMIDF(NNVE,*),KADJF(NNVE,*)
      
C     local variables

      INTEGER IEL,IEL1,IEL2
      DOUBLE PRECISION DUC1,DUC2,DUC3,DUM1,DUM2,DUM3
      
      DOUBLE PRECISION Q2,Q4,Q8
      PARAMETER (Q2=0.5D0)
      PARAMETER (Q4=0.25D0)
      PARAMETER (Q8=0.125D0)
      
C     First we remember the refinement scheme to clarify the
C     local numbering of the vertices in the coarse- and fine-
C     grid triangles.
C     Let a coarse grid triangle be locally numbered as:
C
C       2 
C       |  \
C       |    \
C       | IEL  \
C       |        \
C       3----------1
C     
C     Then the refinement process assigns the following numbers:
C
C       2
C       |  \
C       |     \
C       |        \
C       2*-------- 1* 
C       | \   IEL  |   \
C       |    \     |      \
C       |       \  |         \
C       3----------3*----------1
C
C     i.e. the element number "IEL" is put into the middle with
C     the corner vertices numbered according to the local edge
C     numbers of the coarse grid element.
C
C     To access information on the edges of the fine grid element,
C     we have to work with adjacencies!
C      
C     First copy DUC to DUF. This will transfer the corner values
C     from the coarse grid to the fine grid. More precisely, this
C     will map:
C       Coarse grid vertices 1..NVTC  -> Fine grid vertices 1..NVTC
C       Coarse grid midpoints 1..NMTC -> Fine grid vertices NVTC+1..NVTC+NMTC = NVTF
C     Afterwards, we only have to create the missing midpoint values!

      CALL LCP1(DUC,DUF,NVTC+NMTC)
      
C     loop over the elements

      DO IEL = 1,NELC
      
C       We fetch the function values of the coarse grid element
C       into variables following the following scheme:
C
C          DUC2
C           |   \
C           |      \
C           |         \
C          DUM2         DUM1
C           |                \
C           |                   \
C           |                      \
C          DUC3 --------DUM3-------  DUC1

        DUC1 = DUC(KVERTC(1,IEL))
        DUC2 = DUC(KVERTC(2,IEL))
        DUC3 = DUC(KVERTC(3,IEL))
        DUM1 = DUC(KMIDC(1,IEL))
        DUM2 = DUC(KMIDC(2,IEL))
        DUM3 = DUC(KMIDC(3,IEL))

C       We have to calculate the function values in the new DOF's on the
C       fine grid, which are located here:
C
C          DUC2
C           |   \
C           X     X
C           |        \
C          DUM2 --X-- DUM1
C           |   \      |   \
C           X     X    X     X
C           |       \  |        \
C          DUC3---X---DUM3---X---DUC1
C
C       On this trangle, the function is a quadratic polynomial:
C
C         P(X,Y) = c1 + c2*x + c3*y + c4*x^2 + c5*x*y + c6*y^2
C
C       Solving for the coefficients such that the polynomial takes
C       the DUxy-values in the corners/midpoints, we obtain the 
C       polynomial as:
C
C         P(X,Y) = DUC3 + 
C                  (-3*DUC3-DUC1+4*DUM3)*x + 
C                  (-3*DUC3-DUC2+4*DUM2)*y + 
C                  (4*DUC3-4*DUM3+4*DUM1-4*DUM2)*x*y + 
C                  (2*DUC3+2*DUC1-4*DUM3)*x^2 + 
C                  (2*DUC3+2*DUC2-4*DUM2)*y^2
C
C       This has to be evaluated in the new points, marked as "X"
C       in the above sketch.
C
C       Remember, that the coarse grid element IEL is moved
C       to the inner fine-grid element. The corners of the coars
C       grid element are always the first vertices of the triangles
C       on the fine grid elements.
C
C       First calculate the edge mitpoint values of that one!
C
C         |        \
C        DUM2---X---DUM1
C         | \   IEL  |   \
C              X     X        
C                 \  |           
C                 --DUM3--        

C        DUF(KMIDF(1,IEL)) = P(1/4,1/2)
C                          = -1/8*DUC3-1/8*DUC1+1/4*DUM3+1/2*DUM2+1/2*DUM1
C        DUF(KMIDF(2,IEL)) = P(1/4,1/2)
C                          = -1/8*DUC1+1/2*DUM3-1/8*DUC2+1/2*DUM2+1/4*DUM1
C        DUF(KMIDF(3,IEL)) = P(1/2,1/4)
C                          = -1/8*DUC3+1/2*DUM3-1/8*DUC2+1/4*DUM2+1/2*DUM1

        DUF(KMIDF(1,IEL)) = -Q8*DUC3-Q8*DUC1+Q4*DUM3+Q2*DUM2+Q2*DUM1
        DUF(KMIDF(2,IEL)) = -Q8*DUC1+Q2*DUM3-Q8*DUC2+Q2*DUM2+Q4*DUM1
        DUF(KMIDF(3,IEL)) = -Q8*DUC3+Q2*DUM3-Q8*DUC2+Q4*DUM2+Q2*DUM1

C       Now to the values on the other edges.
C       Only calculate information on edges that are adjacent to
C       an element with lower element number. This prevents
C       information from being computed twice.
C
C       Check the first edge:

        IF (KADJC(1,IEL).LT.IEL) THEN
        
C         Neighbor element has a smaller number than our element.
C         Calculate the information on the current edge.
C         This is the edge corresponding to the edge midpoint DUM1!
C         We have to calculate the values in the new edge midpoints.
C
C          DUC2
C           |   \
C           |     X
C           |  IEL1  \
C          DUM2 ----- DUM1
C               \ IEL  |   \
C                 \    |     X
C                   \  |  IEL2  \
C                     DUM3-------DUC1
C
C         Use adjacencies to get the fine-grid elements IEL1 and IEL2.

          IEL1 = KADJF(1,IEL)
          IEL2 = KADJF(3,IEL)
          
C         Calculate the new edge midpoints:

C          DUF(KMIDF(3,IEL1)) = P(1/4,3/4)
C                             = -1/8*DUC1+3/8*DUC2+3/4*DUM1
C          DUF(KMIDF(1,IEL2)) = P(3/4,1/4)
C                             = 3/8*DUC1-1/8*DUC2+3/4*DUM1
          
          DUF(KMIDF(3,IEL1)) = -Q8*DUC1+3D0*Q8*DUC2+3D0*Q4*DUM1
          DUF(KMIDF(1,IEL2)) = 3D0*Q8*DUC1-Q8*DUC2+3D0*Q4*DUM1
        
        END IF
      
C       Check the next edge:

        IF (KADJC(2,IEL).LT.IEL) THEN
        
C          DUC2
C           |   \
C           X     \
C           |   IEL1 \
C          DUM2 ----- DUM1
C           |   \      |     
C           X     \IEL |       
C           |  IEL2 \  |          
C          DUC3-------DUM3               
C
C         Use adjacencies to get the fine-grid elements IEL1 and IEL2.

          IEL1 = KADJF(1,IEL)
          IEL2 = KADJF(2,IEL)
          
C         Calculate the new edge midpoints:

C          DUF(KMIDF(1,IEL1)) = P(0,3/4)
C                             = -1/8*DUC3+3/8*DUC2+3/4*DUM2
C          DUF(KMIDF(3,IEL2)) = P(0,1/4)
C                             = 3/8*DUC3-1/8*DUC2+3/4*DUM2
          
          DUF(KMIDF(1,IEL1)) = -Q8*DUC3+3D0*Q8*DUC2+3D0*Q4*DUM2
          DUF(KMIDF(3,IEL2)) = 3D0*Q8*DUC3-1D0*Q8*DUC2+3D0*Q4*DUM2
        
        END IF

C       Check the last edge

        IF (KADJC(3,IEL).LT.IEL) THEN
        
C          DUM2 ----- DUM1
C           |   \  IEL |   \
C           | IEL2\    | IEL1\
C           |       \  |        \
C          DUC3---X---DUM3---X---DUC1
C
C         Use adjacencies to get the fine-grid elements IEL1 and IEL2.

          IEL1 = KADJF(3,IEL)
          IEL2 = KADJF(2,IEL)
          
C         Calculate the new edge midpoints:

C          DUF(KMIDF(3,IEL1)) = P(3/4,0)
C                             = -1/8*DUC3+3/8*DUC1+3/4*DUM3
C          DUF(KMIDF(1,IEL2)) = P(1/4,0)
C                             = 3/8*DUC3-1/8*DUC1+3/4*DUM3
          
          DUF(KMIDF(3,IEL1)) = -Q8*DUC3+3D0*Q8*DUC1+3D0*Q4*DUM3
          DUF(KMIDF(1,IEL2)) = 3D0*Q8*DUC3-Q8*DUC1+3D0*Q4*DUM3
        
        END IF
        
C       That's it - next element.
        
      END DO ! IEL
      
      END
      
************************************************************************
* P2 restriction
*
* Restricts a P2-RHS vector from a finer level to a coarser level.
*
* In:
*   DUF    : array [1..NVTF+NMTF] of double
*            Fine grid RHS vector
*   NVTC   : Number of vertices on the coarse grid
*   NMTC   : Number of edge midpoints on the coarse grid
*   NELC   : Number of elements on the coarse grid
*   DCORVC : array [1..2,1..NVTC] of double
*            Coordinates of all vertices on the coarse grid
*   KVERTC : array [1..NNVE,1..NVTC] of double
*            Vertices on each element on the coarse grid
*   KMIDC  : array [1..NNVE,1..NVTC] of double
*            Midpoints on each element on the coarse grid
*   KADJC  : array [1..NNVE,1..NVTC]
*            Numbers of adjacent elements to each element on the
*            coarse grid.
*   KVERTF : array [1..NNVE,1..NVTF] of double
*            Vertices on each element on the fine grid
*   KMIDF  : array [1..NNVE,1..NVTF] of double
*            Midpoints on each element on the fine grid
*
* Out:
*   DUC    : array [1..NVTC+NMTC] of double
*            Coarse grid RHS vector.
************************************************************************

      SUBROUTINE MR002 (DUC,DUF,NVTC,NMTC,NELC,
     *                  KVERTC,KMIDC,KADJC,
     *                  KMIDF,KADJF)
      
      IMPLICIT NONE
      
      INTEGER NNVE
      PARAMETER (NNVE=4)
      
C     parameters
      
      DOUBLE PRECISION DUC(*),DUF(*)
      INTEGER KVERTC(NNVE,*),KMIDC(NNVE,*),KADJC(NNVE,*),NVTC,NMTC,NELC
      INTEGER KMIDF(NNVE,*),KADJF(NNVE,*)
      
C     local variables

      INTEGER IEL,IEL2,IEL3,IEL4
      DOUBLE PRECISION DUF11,DUF12,DUF13,DUF21,DUF23,DUF33,DUF31
      DOUBLE PRECISION DUF41,DUF43
      DOUBLE PRECISION N1,N2,N3
      
      DOUBLE PRECISION Q2,Q4,Q8
      PARAMETER (Q2=0.5D0)
      PARAMETER (Q4=0.25D0)
      PARAMETER (Q8=0.125D0)
      
C     First we remember the refinement scheme to clarify the
C     local numbering of the vertices in the coarse- and fine-
C     grid triangles.
C     Let a coarse grid triangle be locally numbered as:
C
C       2 
C       |  \
C       |    \
C       | IEL  \
C       |        \
C       3----------1
C     
C     Then the refinement process assigns the following numbers:
C
C       2
C       |  \
C       |     \
C       |        \
C       2*-------- 1* 
C       | \   IEL  |   \
C       |    \     |      \
C       |       \  |         \
C       3----------3*----------1
C
C     i.e. the element number "IEL" is put into the middle with
C     the corner vertices numbered according to the local edge
C     numbers of the coarse grid element.
C
C     To access information on the edges of the fine grid element,
C     we have to work with adjacencies!
C
C     Copy the first NVTC+NMTC values from DUF to DUC. This will
C     transfer the contribution of the values from the 
C     fine-grid vertices that are coarse-grid vertices or midpoints
C     as well. More precisely, this will transfer:
C
C       Fine grid vertices 1..NVTC -> Coarse grid vertices 1..NVTC  
C       Fine grid vertices NVTC+1..NVTC+NMTC = NVTF -> Coarse grid midpoints 1..NMTC
C
C     Afterwards, we have to add only the contribution of the fine grid
C     edge midpoints to the coarse grid vertices/midpoints.

      CALL LCP1(DUF,DUC,NVTC+NMTC)
      
C     loop over the elements

      DO IEL = 1,NELC
      
C       The prolongation created from DUC1-3 and DUM1-3 in the 
C       following sketch the fine grid values DUFxy:
C
C          DUC2
C           |  \
C           |    \
C         DUF21  DUF23
C           |        \
C           |          \
C          DUM2--DUF11--DUM1
C           |  \         |  \
C           |    \       |    \
C         DUF33  DUF12 DUF13  DUF41
C           |        \   |        \
C           |          \ |          \
C          DUC3--DUF31--DUM3--DUF43--DUC1
C
C       This was done by a weighted averaging. The restriction now
C       builds the DUC1-3 and DUM1-3 values as a weighted sum
C       of themselves (DUM1-3, DUC1-3) and the new midpoints
C       DUFxy using the same weights.
C
C       We had 
C         DUCx(fine grid) := 1*DUCx(coarse grid)
C         DUMy(fine grid) := 1*DUCy(coarse grid)
C       Therefore:
C         DUCx(coarse grid) := 1*DUCx(fine grid) + ...
C         DUCy(coarse grid) := 1*DUMy(fine grid) + ...
C
C       This part of the sum is already written to DUC with the
C       above LCP1 command! now comes the rest of the sum.
C
C       The prolongation used the following formulas to calculate
C       the fine grid vertices:
C
C        DUF11 = P(1/4,1/2)
C              = -1/8*DUC3-1/8*DUC1+1/4*DUM3+1/2*DUM2+1/2*DUM1
C        DUF12 = P(1/4,1/4)
C              = -1/8*DUC1+1/2*DUM3-1/8*DUC2+1/2*DUM2+1/4*DUM1
C        DUF13 = P(1/2,1/4)
C              = -1/8*DUC3+1/2*DUM3-1/8*DUC2+1/4*DUM2+1/2*DUM1
C
C        DUF23 = P(1/4,3/4)
C              = -1/8*DUC1+3/8*DUC2+3/4*DUM1
C        DUF41 = P(3/4,1/4)
C              = 3/8*DUC1-1/8*DUC2+3/4*DUM1
C
C        DUF21 = P(0,3/4)
C              = -1/8*DUC3+3/8*DUC2+3/4*DUM2
C        DUF33 = P(0,1/4)
C              = 3/8*DUC3-1/8*DUC2+3/4*DUM2
C
C        DUF43 = P(3/4,0)
C              = -1/8*DUC3+3/8*DUC1+3/4*DUM3
C        DUF31 = P(1/4,0)
C              = 3/8*DUC3-1/8*DUC1+3/4*DUM3
C
C       This is equivalent to the system
C
C        DUF11    [-1/8,    0, -1/8, 1/2,  1/2,  1/4]   DUC1
C        DUF12    [-1/8, -1/8 ,   0, 1/4,  1/2,  1/2]   DUC2
C        DUF13    [   0, -1/8, -1/8, 1/2,  1/4,  1/2]   DUC3
C        DUF21    [   0,  3/8, -1/8,   0,  3/4,    0]   DUM1
C        DUF23 =  [-1/8,  3/8,    0, 3/4,    0,    0] * DUM2
C        DUF31    [-1/8,    0,  3/8,   0,    0,  3/4]   DUM3
C        DUF33    [   0, -1/8,  3/8,   0,  3/4,    0]
C        DUF41    [ 3/8, -1/8 ,   0, 3/4,    0,    0]
C        DUF43    [ 3/8,    0, -1/8,   0,    0,  3/4]
C
C       Transposing it says what is left to add to DUC1-3/DUM1-3:
C
C         DUC1    [-1/8, -1/8,    0,    0, -1/8, -1/8,    0,  3/8,  3/8]   DUF11
C         DUC2    [   0, -1/8, -1/8,  3/8,  3/8,    0, -1/8, -1/8,    0]   DUF12
C         DUC3 += [-1/8,    0, -1/8, -1/8,    0,  3/8,  3/8,    0, -1/8] * DUF13
C         DUM1    [ 1/2,  1/4,  1/2,    0,  3/4,    0,    0,  3/4,    0]   DUF21
C         DUM2    [ 1/2,  1/2,  1/4, 3/4,     0,    0,  3/4,    0,    0]   DUF23
C         DUM3    [ 1/4,  1/2,  1/2,    0,    0,  3/4,    0,    0,  3/4]   DUF31
C                                                                          DUF33
C                                                                          DUF41
C                                                                          DUF43
C
C       Fetch the fine grid values into local variables:

        IEL2 = KADJF(1,IEL)
        IEL3 = KADJF(2,IEL)
        IEL4 = KADJF(3,IEL)
        DUF11 = DUF(KMIDF(1,IEL))
        DUF12 = DUF(KMIDF(2,IEL))
        DUF13 = DUF(KMIDF(3,IEL))
        DUF21 = DUF(KMIDF(1,IEL2))
        DUF23 = DUF(KMIDF(3,IEL2))
        DUF31 = DUF(KMIDF(1,IEL3))
        DUF33 = DUF(KMIDF(3,IEL3))
        DUF41 = DUF(KMIDF(1,IEL4))
        DUF43 = DUF(KMIDF(3,IEL4))

C       When we add the information to DUC1-3/DUM1-3 we have to take
C       into account whether there is a neighbor element or not!
C       If there is a neighbor, we only add half of the information
C       of the edge with the neighbor to DUC1-3/DUM1-3. 
C       In the loop over the elements here, we will
C       later reach the neighbor and add another time half of the
C       information, which that way completes that edge.
C
C       Set Nx=0.5 if there is a neighbor on edge x on the
C       coarse grid or =1.0 if there is no neighbor

        N1 = 0.5D0
        N2 = 0.5D0
        N3 = 0.5D0
        
        IF (KADJC(1,IEL).EQ.0) N1 = 1D0
        IF (KADJC(2,IEL).EQ.0) N2 = 1D0
        IF (KADJC(3,IEL).EQ.0) N3 = 1D0
        
C       Now sum the restriction together.
C
C       DUC1:

        DUC(KVERTC(1,IEL)) = DUC(KVERTC(1,IEL)) 
     *         -Q8*DUF11 -Q8*DUF12 
     *    +N1*(-Q8*DUF23 +3D0*Q8*DUF41)
     *    +N3*(-Q8*DUF31 +3D0*Q8*DUF43)
          
C       DUC2:

        DUC(KVERTC(2,IEL)) = DUC(KVERTC(2,IEL)) 
     *         -Q8*DUF12 -Q8*DUF13
     *    +N2*(3D0*Q8*DUF21-Q8*DUF33)
     *    +N1*(3D0*Q8*DUF23-Q8*DUF41)
     
C       DUC3:

        DUC(KVERTC(3,IEL)) = DUC(KVERTC(3,IEL)) 
     *         -Q8*DUF11 -Q8*DUF13
     *    +N2*(-Q8*DUF21 +3D0*Q8*DUF33)
     *    +N3*(-Q8*DUF43 +3D0*Q8*DUF31)
     
C       DUM1:

        DUC(KMIDC(1,IEL)) = DUC(KMIDC(1,IEL)) 
     *   +     Q2*DUF11 +Q4*DUF12 +Q2*DUF13
     *   + N1*(3D0*Q4*DUF23 +3*Q4*DUF41)

C       DUM2:

        DUC(KMIDC(2,IEL)) = DUC(KMIDC(2,IEL)) 
     *   +     Q2*DUF11 +Q2*DUF12 +Q4*DUF13
     *   + N2*(3D0*Q4*DUF21+3D0*Q4*DUF33)

C       DUM3:

        DUC(KMIDC(3,IEL)) = DUC(KMIDC(3,IEL)) 
     *   +     Q4*DUF11 +Q2*DUF12 +Q2*DUF13
     *   +N3*(3D0*Q4*DUF43 +3D0*Q4*DUF31)

C       That's it - next element.
        
      END DO ! IEL
      
      END
            
************************************************************************
* P2 prolongation, linear interpolation
*
* Transfers a P2-solution from a coarser level to a finer level.
* The prolongation is performed "linearly" on a once refined grid.
*
* In:
*   DUC    : array [1..NVTC+NMTC] of double
*            Coarse grid solution vector.
*   NVTC   : Number of vertices on the coarse grid
*   NMTC   : Number of edge midpoints on the coarse grid
*   NELC   : Number of elements on the coarse grid
*   DCORVC : array [1..2,1..NVTC] of double
*            Coordinates of all vertices on the coarse grid
*   KVERTC : array [1..NNVE,1..NVTC] of double
*            Vertices on each element on the coarse grid
*   KMIDC  : array [1..NNVE,1..NVTC] of double
*            Midpoints on each element on the coarse grid
*   KADJC  : array [1..NNVE,1..NVTC]
*            Numbers of adjacent elements to each element on the
*            coarse grid.
*   KVERTF : array [1..NNVE,1..NVTF] of double
*            Vertices on each element on the fine grid
*   KMIDF  : array [1..NNVE,1..NVTF] of double
*            Midpoints on each element on the fine grid
*
* Out:
*   DUF    : array [1..NVTF+NMTF] of double
*            Fine grid solution vector
************************************************************************

      SUBROUTINE MP002L (DUC,DUF,NVTC,NMTC,NELC,
     *                   KVERTC,KMIDC,KADJC,
     *                   KMIDF,KADJF)
      
      IMPLICIT NONE
      
      INTEGER NNVE
      PARAMETER (NNVE=4)
      
C     parameters
      
      DOUBLE PRECISION DUC(*),DUF(*)
      INTEGER KVERTC(NNVE,*),KMIDC(NNVE,*),KADJC(NNVE,*),NVTC,NMTC,NELC
      INTEGER KMIDF(NNVE,*),KADJF(NNVE,*)
      
C     local variables

      INTEGER IEL,IEL1,IEL2
      DOUBLE PRECISION DUC1,DUC2,DUC3,DUM1,DUM2,DUM3
      
      DOUBLE PRECISION Q2,Q4,Q8
      PARAMETER (Q2=0.5D0)
      PARAMETER (Q4=0.25D0)
      PARAMETER (Q8=0.125D0)
      
C     First we remember the refinement scheme to clarify the
C     local numbering of the vertices in the coarse- and fine-
C     grid triangles.
C     Let a coarse grid triangle be locally numbered as:
C
C       2 
C       |  \
C       |    \
C       | IEL  \
C       |        \
C       3----------1
C     
C     Then the refinement process assigns the following numbers:
C
C       2
C       |  \
C       |     \
C       |        \
C       2*-------- 1* 
C       | \   IEL  |   \
C       |    \     |      \
C       |       \  |         \
C       3----------3*----------1
C
C     i.e. the element number "IEL" is put into the middle with
C     the corner vertices numbered according to the local edge
C     numbers of the coarse grid element.
C
C     To access information on the edges of the fine grid element,
C     we have to work with adjacencies!
C      
C     First copy DUC to DUF. This will transfer the corner values
C     from the coarse grid to the fine grid. More precisely, this
C     will map:
C       Coarse grid vertices 1..NVTC  -> Fine grid vertices 1..NVTC
C       Coarse grid midpoints 1..NMTC -> Fine grid vertices NVTC+1..NVTC+NMTC = NVTF
C     Afterwards, we only have to create the missing midpoint values!
C
C     As we use a simple "linear" mapping, this is all we have to do
C     for the corners/midpoints of the first grid. The values
C     in the new midpoints on the fine grid are calculated by
C     linear interpolation later.

      CALL LCP1(DUC,DUF,NVTC+NMTC)
      
C     loop over the elements

      DO IEL = 1,NELC
      
C       We fetch the function values of the coarse grid element
C       into variables following the following scheme:
C
C          DUC2
C           |   \
C           |      \
C           |         \
C          DUM2         DUM1
C           |                \
C           |                   \
C           |                      \
C          DUC3 --------DUM3-------  DUC1

        DUC1 = DUC(KVERTC(1,IEL))
        DUC2 = DUC(KVERTC(2,IEL))
        DUC3 = DUC(KVERTC(3,IEL))
        DUM1 = DUC(KMIDC(1,IEL))
        DUM2 = DUC(KMIDC(2,IEL))
        DUM3 = DUC(KMIDC(3,IEL))

C       We have to calculate the function values in the new DOF's on the
C       fine grid, which are located here:
C
C          DUC2
C           |   \
C           X     X
C           |        \
C          DUM2 --X-- DUM1
C           |   \      |   \
C           X     X    X     X
C           |       \  |        \
C          DUC3---X---DUM3---X---DUC1
C
C       We have to calculate all new midpoints, marked by X.
C       These are calculated by simple linear interpolation.
C
C       First calculate the edge mitpoint values of the inner
C       fine-grid element.
C
C         |        \
C        DUM2---X---DUM1
C         | \   IEL  |   \
C              X     X        
C                 \  |           
C                 --DUM3--        

        DUF(KMIDF(1,IEL)) = Q2*(DUM1+DUM2)
        DUF(KMIDF(2,IEL)) = Q2*(DUM2+DUM3)
        DUF(KMIDF(3,IEL)) = Q2*(DUM3+DUM1)

C       Now to the values on the other edges.
C       Only calculate information on edges that are adjacent to
C       an element with lower element number. This prevents
C       information from being computed twice.
C
C       Check the first edge:

        IF (KADJC(1,IEL).LT.IEL) THEN
        
C         Neighbor element has a smaller number than our element.
C         Calculate the information on the current edge.
C         This is the edge corresponding to the edge midpoint DUM1!
C         We have to calculate the values in the new edge midpoints.
C
C          DUC2
C           |   \
C           |     X
C           |  IEL1  \
C          DUM2 ----- DUM1
C               \ IEL  |   \
C                 \    |     X
C                   \  |  IEL2  \
C                     DUM3-------DUC1
C
C         Use adjacencies to get the fine-grid elements IEL1 and IEL2.

          IEL1 = KADJF(1,IEL)
          IEL2 = KADJF(3,IEL)
          
C         Calculate the new edge midpoints:

          DUF(KMIDF(3,IEL1)) = Q2*(DUM1+DUC2)
          DUF(KMIDF(1,IEL2)) = Q2*(DUC1+DUM1)
        
        END IF
      
C       Check the next edge:

        IF (KADJC(2,IEL).LT.IEL) THEN
        
C          DUC2
C           |   \
C           X     \
C           |   IEL1 \
C          DUM2 ----- DUM1
C           |   \      |     
C           X     \IEL |       
C           |  IEL2 \  |          
C          DUC3-------DUM3               
C
C         Use adjacencies to get the fine-grid elements IEL1 and IEL2.

          IEL1 = KADJF(1,IEL)
          IEL2 = KADJF(2,IEL)
          
C         Calculate the new edge midpoints:

          DUF(KMIDF(1,IEL1)) = Q2*(DUC2+DUM2)
          DUF(KMIDF(3,IEL2)) = Q2*(DUM2+DUC3)
        
        END IF

C       Check the last edge

        IF (KADJC(3,IEL).LT.IEL) THEN
        
C          DUM2 ----- DUM1
C           |   \  IEL |   \
C           | IEL2\    | IEL1\
C           |       \  |        \
C          DUC3---X---DUM3---X---DUC1
C
C         Use adjacencies to get the fine-grid elements IEL1 and IEL2.

          IEL1 = KADJF(3,IEL)
          IEL2 = KADJF(2,IEL)
          
C         Calculate the new edge midpoints:

          DUF(KMIDF(3,IEL1)) = Q2*(DUM3+DUC1)
          DUF(KMIDF(1,IEL2)) = Q2*(DUC3+DUM3)
        
        END IF
        
C       That's it - next element.
        
      END DO ! IEL
      
      END
      
************************************************************************
* P2 restriction, linear interpolation
*
* Restricts a P2-RHS vector from a finer level to a coarser level.
* The restriction is done "linearly" on a once refined grid.
*
* In:
*   DUF    : array [1..NVTF+NMTF] of double
*            Fine grid RHS vector
*   NVTC   : Number of vertices on the coarse grid
*   NMTC   : Number of edge midpoints on the coarse grid
*   NELC   : Number of elements on the coarse grid
*   DCORVC : array [1..2,1..NVTC] of double
*            Coordinates of all vertices on the coarse grid
*   KVERTC : array [1..NNVE,1..NVTC] of double
*            Vertices on each element on the coarse grid
*   KMIDC  : array [1..NNVE,1..NVTC] of double
*            Midpoints on each element on the coarse grid
*   KADJC  : array [1..NNVE,1..NVTC]
*            Numbers of adjacent elements to each element on the
*            coarse grid.
*   KVERTF : array [1..NNVE,1..NVTF] of double
*            Vertices on each element on the fine grid
*   KMIDF  : array [1..NNVE,1..NVTF] of double
*            Midpoints on each element on the fine grid
*
* Out:
*   DUC    : array [1..NVTC+NMTC] of double
*            Coarse grid RHS vector.
************************************************************************

      SUBROUTINE MR002L (DUC,DUF,NVTC,NMTC,NELC,
     *                   KVERTC,KMIDC,KADJC,
     *                   KMIDF,KADJF)
      
      IMPLICIT NONE
      
      INTEGER NNVE
      PARAMETER (NNVE=4)
      
C     parameters
      
      DOUBLE PRECISION DUC(*),DUF(*)
      INTEGER KVERTC(NNVE,*),KMIDC(NNVE,*),KADJC(NNVE,*),NVTC,NMTC,NELC
      INTEGER KMIDF(NNVE,*),KADJF(NNVE,*)
      
C     local variables

      INTEGER IEL,IEL2,IEL3,IEL4
      DOUBLE PRECISION DUF11,DUF12,DUF13,DUF21,DUF23,DUF33,DUF31
      DOUBLE PRECISION DUF41,DUF43
      DOUBLE PRECISION N1,N2,N3
      
      DOUBLE PRECISION Q2,Q4,Q8
      PARAMETER (Q2=0.5D0)
      PARAMETER (Q4=0.25D0)
      PARAMETER (Q8=0.125D0)
      
C     First we remember the refinement scheme to clarify the
C     local numbering of the vertices in the coarse- and fine-
C     grid triangles.
C     Let a coarse grid triangle be locally numbered as:
C
C       2 
C       |  \
C       |    \
C       | IEL  \
C       |        \
C       3----------1
C     
C     Then the refinement process assigns the following numbers:
C
C       2
C       |  \
C       |     \
C       |        \
C       2*-------- 1* 
C       | \   IEL  |   \
C       |    \     |      \
C       |       \  |         \
C       3----------3*----------1
C
C     i.e. the element number "IEL" is put into the middle with
C     the corner vertices numbered according to the local edge
C     numbers of the coarse grid element.
C
C     To access information on the edges of the fine grid element,
C     we have to work with adjacencies!
C
C     Copy the first NVTC+NMTC values from DUF to DUC. This will
C     transfer the contribution of the values from the 
C     fine-grid vertices that are coarse-grid vertices or midpoints
C     as well. More precisely, this will transfer:
C
C       Fine grid vertices 1..NVTC -> Coarse grid vertices 1..NVTC  
C       Fine grid vertices NVTC+1..NVTC+NMTC = NVTF -> Coarse grid midpoints 1..NMTC
C
C     Afterwards, we have to add only the contribution of the fine grid
C     edge midpoints to the coarse grid vertices/midpoints.

      CALL LCP1(DUF,DUC,NVTC+NMTC)
      
C     loop over the elements

      DO IEL = 1,NELC
      
C       The prolongation created from DUC1-3 and DUM1-3 in the 
C       following sketch the fine grid values DUFxy:
C
C          DUC2
C           |  \
C           |    \
C         DUF21  DUF23
C           |        \
C           |          \
C          DUM2--DUF11--DUM1
C           |  \         |  \
C           |    \       |    \
C         DUF33  DUF12 DUF13  DUF41
C           |        \   |        \
C           |          \ |          \
C          DUC3--DUF31--DUM3--DUF43--DUC1
C
C       This was done by a weighted averaging. The restriction now
C       builds the DUC1-3 and DUM1-3 values as a weighted sum
C       of themselves (DUM1-3, DUC1-3) and the new midpoints
C       DUFxy using the same weights.
C
C       We had 
C         DUCx(fine grid) := 1*DUCx(coarse grid)
C         DUMy(fine grid) := 1*DUCy(coarse grid)
C       Therefore:
C         DUCx(coarse grid) := 1*DUCx(fine grid) + ...
C         DUCy(coarse grid) := 1*DUMy(fine grid) + ...
C
C       This part of the sum is already written to DUC with the
C       above LCP1 command! now comes the rest of the sum.
C
C       The prolongation used the following formulas to calculate
C       the fine grid vertices:
C
C        DUF11 = (DUM1+DUM2)/2
C        DUF12 = (DUM2+DUM3)/2
C        DUF13 = (DUM3+DUM1)/2
C
C        DUF23 = (DUM1+DUC2)/2
C        DUF41 = (DUC1+DUM1)/2
C
C        DUF21 = (DUC2+DUM2)/2
C        DUF33 = (DUC1+DUM2)/2
C
C        DUF43 = (DUC1+DUM3)/2
C        DUF31 = (DUC3+DUM3)/2
C
C       This is equivalent to the system
C
C        DUF11    [    ,    ,    , 1/2, 1/2,    ]   DUC1
C        DUF12    [    ,    ,    ,    , 1/2, 1/2]   DUC2
C        DUF13    [    ,    ,    , 1/2,    , 1/2]   DUC3
C        DUF21    [    , 1/2,    ,    , 1/2,    ]   DUM1
C        DUF23 =  [    , 1/2,    , 1/2,    ,    ] * DUM2
C        DUF31    [    ,    , 1/2,    ,    , 1/2]   DUM3
C        DUF33    [    ,    , 1/2,    , 1/2,    ]
C        DUF41    [ 1/2,    ,    , 1/2,    ,    ]
C        DUF43    [ 1/2,    ,    ,    ,    , 1/2]
C
C       Transposing it says what is left to add to DUC1-3/DUM1-3:
C
C         DUC1    [   0,   0,   0,   0,   0,   0,   0, 1/2, 1/2]   DUF11
C         DUC2    [   0,   0,   0, 1/2, 1/2,   0,   0,   0,   0]   DUF12
C         DUC3 += [   0,   0,   0,   0,   0, 1/2, 1/2,   0,   0] * DUF13
C         DUM1    [ 1/2,   0, 1/2,   0, 1/2,   0,   0, 1/2,   0]   DUF21
C         DUM2    [ 1/2, 1/2,   0, 1/2,   0,   0, 1/2,   0,   0]   DUF23
C         DUM3    [   0, 1/2, 1/2,   0,   0, 1/2,   0,   0, 1/2]   DUF31
C                                                                  DUF33
C                                                                  DUF41
C                                                                  DUF43
C
C       Fetch the fine grid values into local variables:

        IEL2 = KADJF(1,IEL)
        IEL3 = KADJF(2,IEL)
        IEL4 = KADJF(3,IEL)
        DUF11 = DUF(KMIDF(1,IEL))
        DUF12 = DUF(KMIDF(2,IEL))
        DUF13 = DUF(KMIDF(3,IEL))
        DUF21 = DUF(KMIDF(1,IEL2))
        DUF23 = DUF(KMIDF(3,IEL2))
        DUF31 = DUF(KMIDF(1,IEL3))
        DUF33 = DUF(KMIDF(3,IEL3))
        DUF41 = DUF(KMIDF(1,IEL4))
        DUF43 = DUF(KMIDF(3,IEL4))

C       When we add the information to DUC1-3/DUM1-3 we have to take
C       into account whether there is a neighbor element or not!
C       If there is a neighbor, we only add half of the information
C       of the edge with the neighbor to DUC1-3/DUM1-3. 
C       In the loop over the elements here, we will
C       later reach the neighbor and add another time half of the
C       information, which that way completes that edge.
C
C       Set Nx=0.5 if there is a neighbor on edge x on the
C       coarse grid or =1.0 if there is no neighbor

        N1 = 0.5D0
        N2 = 0.5D0
        N3 = 0.5D0
        
        IF (KADJC(1,IEL).EQ.0) N1 = 1D0
        IF (KADJC(2,IEL).EQ.0) N2 = 1D0
        IF (KADJC(3,IEL).EQ.0) N3 = 1D0
        
C       Now sum the restriction together.
C
C       DUC1:

        DUC(KVERTC(1,IEL)) = DUC(KVERTC(1,IEL)) 
     *    + N1*Q2*DUF41 + N3*Q2*DUF43
          
C       DUC2:

        DUC(KVERTC(2,IEL)) = DUC(KVERTC(2,IEL)) 
     *    + N1*Q2*DUF23 + N2*Q2*DUF21
     
C       DUC3:

        DUC(KVERTC(3,IEL)) = DUC(KVERTC(3,IEL)) 
     *     + N2*Q2*DUF33 + N3*Q2*DUF31   
     
C       DUM1:

        DUC(KMIDC(1,IEL)) = DUC(KMIDC(1,IEL)) 
     *     +    Q2*DUF13 + Q2*DUF11 
     *     + N1*(Q2*DUF41 + Q2*DUF23)

C       DUM2:

        DUC(KMIDC(2,IEL)) = DUC(KMIDC(2,IEL)) 
     *     +    Q2*DUF11 + Q2*DUF12
     *     + N2*(Q2*DUF21 + Q2*DUF33)

C       DUM3:

        DUC(KMIDC(3,IEL)) = DUC(KMIDC(3,IEL)) 
     *     +    Q2*DUF12 + Q2*DUF13
     *     + N3*(Q2*DUF31 + Q2*DUF43)

C       That's it - next element.
        
      END DO ! IEL
      
      END
            