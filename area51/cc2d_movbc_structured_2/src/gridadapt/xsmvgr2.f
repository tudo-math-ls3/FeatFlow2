************************************************************************
* Move Gridpoint, Method 1 = Explicit Euler
*
* This routine is used in the static grid adaption method. It moves
* one grid-vertex IVT inside of the grid DCRVSR. The resulting 
* coordinates of the point are written into the destination array 
* DCRVDS. The routine is designed only to move (corner-)vertices and
* boundary-vertices, not midpoints!
* The point will be moved along the gradient of the solution vector
* DSOL, which is given by DGRADX, DGRADY.
*
* This routine uses the explicit Euler method with NSTPS steps to move
* the point to its destination coordinates.
*
* In:
*  TRIA    - array [1..SZTRIA] of integer
*            Triangulation structure for "source", undeformed 
*            triangulation
*  DCRVSR  - Array DCORVG with source-coordinates of all points in the
*            grid. Must correspond to LCORVG in the TRIA structure.
*  DCRVDS  - array [1..2,1..NVT] of double
*            Array receiving the destination coordinates of the moved
*            gridpoint.
*  IVT     - Number of the node that should be moved; range 1..NVT
*  IBND    - If IVT is a boundary node, IBND is the index of the node
*            inside of the KVBD-array. If IVT is not a boundary node,
*            IBND=0 must be set.
*  IBCT    - If IVT is a boundary node, IBCT must be the number of
*            the boundary component that contains IVT. Otherwise,
*            IBCT=0 must be set.
*  LSOL    - Handle to array [1..NEQ] of double precision
*            Solution of the grid adaption system
*  LGRADX  - Handle to array [1..NEQ] of double precision
*            X-gradient of the grid adaption system
*  LGRADY  - Handle to array [1..NEQ] of double precision
*            Y-gradient of the grid adaption system
*  LMON    - Handle to array [1..NEQ] of double precision
*            Monitor function
*  LMNOLD  - Handle to array [1..NEQ] of double precision
*            Old reference function determining the area distribution
*            on the original grid
*
* Method specific input-parameters:
*  NSTPS   - number of steps
*
* Out:
*  DCRVDS [1..2,IVT] - receives the destination coordinates of the
*                      node IVT
*  DPARM             - if IVT is a boundary node, DPARM receives
*                      the new parameter value. Otherwise DPARM is
*                      unchanged.
*  Return value = 0, if successful
*               = 1, if a point was moved but could not be found
*                    in the grid after moving.
************************************************************************

      INTEGER FUNCTION MVGPT1 (TRIA, DCRVSR, DCRVDS, 
     *                   IVT, IBND, IBCT, DPARM,
     *                   LSOL, LGRADX, LGRADY,
     *                   LMON, LMNOLD,
     *                   NSTPS)
     
      IMPLICIT NONE
      
      INCLUDE 'stria.inc'
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
 
C parameters
      
      INTEGER TRIA(SZTRIA), IVT, NSTPS, IBND, IBCT
      DOUBLE PRECISION DCRVSR(2,*), DCRVDS(2,*), DPARM
      INTEGER LSOL, LGRADX, LGRADY,KVEL
      INTEGER LMON, LMNOLD
      
C     local variables
      
      DOUBLE PRECISION DX, DY, DPRM, DXOLD, DYOLD, DT, DDT
      DOUBLE PRECISION DGRX,DGRY,DF, DG, DTMP1
      INTEGER IEL, IELOLD, I
      
C     externals

      EXTERNAL E011
      INTEGER TNDCNT,BSRCH4
      EXTERNAL TNDCNT,BSRCH4
      DOUBLE PRECISION TMAX
      EXTERNAL TMAX
      
C     Some arrays with constants so that we don't have to set them
C     on runtime. This array is used for the multiple scalar evaluation
C     function. It says that for all handles that we pass to that 
C     routine we only want to calculate function values.

      INTEGER FCLC(4)
      DATA FCLC /1,1,1,1/
      SAVE FCLC
      
C     Furthermore a small array for the multiple scalar evaluation
C     function to pass handles as parameters and obtain return values:

      INTEGER LSOLAR (4)
      DOUBLE PRECISION DVAL(4), DVGRX(4), DVGRY(4)
      
C     Get the initial position of the point and the original element
C     that contains the point (IEL = KVEL(1,NVEL))

      DX = DCRVSR(1,IVT)
      DY = DCRVSR(2,IVT)
      
      KVEL = L(TRIA(OLVEL))
      
      IEL = KWORK(KVEL+TRIA(ONVEL)*(IVT-1))
      
C     Calculate the time-step for the explicit Euler

      DDT = 1D0/DBLE(NSTPS)
      
C     Loop through all time steps. Our first time step is 0, our last
C     is NSTPS-1.

      DO I=0,NSTPS-1
      
C     What's the current time step?

        DT = DBLE(I)*DDT
        
C       We want to solve the following ODE:
C
C           d/dt phi(x,t) = eta ( phi(x,t), t),     t=0..1
C
C       for the right hand side eta:(Omega x [0,1])->Omega,
C
C           eta (y,s) = - v(y) / ( s/f(y) + (1-s)/g(y) )
C
C       and with the initial condition:  phi (x,0) = x
C       Here f corresponds to 1/DMON and g to 1/DMNOLD!
        
C Use the interpolation routine to calculate the gradient of LSOL
C in the current point, based on the reconstructed gradient in
C LGRADX/LGRADY.

C        CALL SCEVQR (TNDCNT(TRIA),DSOL,DGRADX,DGRADY,E011,.FALSE.,
C     *               DX,DY,IEL,
C     *               TRIA,KWORK(L(TRIA(OLVERT))),KWORK(L(TRIA(OLMID))),
C     *               DVAL, DGRX, DGRY)

C Furthermore we need for the denominator the both values 1/DMON
C and 1/DMNOLD in the point (DX,DY). Both functiond DMON and
C DMONOLD are given as function values in the vertices of the old,
C undeformed grid. We interpret the values in the vertices as the
C function values of a piecewise linear FE-function. Therefore to
C obtain a function value in (DX,DY), use the scalar evaluation
C subroutine and the E011 element. The gradient this function returns
C can be ignored:

C        CALL SCEVLQ (TNDCNT(TRIA),DWORK(L(LMON)),E011,.FALSE.,
C     *               DX,DY,IEL,
C     *               TRIA,KWORK(L(TRIA(OLVERT))),KWORK(L(TRIA(OLMID))),
C     *               DF, DTMP1, DTMP1)

C        CALL SCEVLQ (TNDCNT(TRIA),DMNOLD,E011,.FALSE.,
C     *               DX,DY,IEL,
C     *               TRIA,KWORK(L(TRIA(OLVERT))),KWORK(L(TRIA(OLMID))),
C     *               DG, DTMP1, DTMP2)

C       For the real evaluation of these values, we use the "multiple 
C       evaluation function" MSCVQR. this allowes us to calculate multiple
C       values at once, thus increasing the speed of the overall adaption
C       by about 100% !

C       For the multiple solution evaluation function we have to
C       set up some parameters. We calculate:
C         - the X-gradient of DSOL with the recovered gradient vector
C         - the Y-gradient of DSOL with the recovered gradient vector
C         - the function value of DMON giving the value DF=f(x)
C         - the function value of DMNOLF giving the value DG=g(x)
C
C       So set up the handle arrays to pass to MSCVQR:

        LSOLAR(1) = LGRADX
        LSOLAR(2) = LGRADY
        LSOLAR(3) = LMON
        LSOLAR(4) = LMNOLD
        
C       *What* to calculate is prescribed by the constant array FCLC.
C       This says to evaluate only function values. As the first
C       two solution vectors in LSOLAR point to the recovered gradient
C       of LSOL, the calculated values in DVAL(1) and DVAL(2) are
C       in fact the recovered gradient vector of LSOL!

        CALL MSCVQR (4,TNDCNT(TRIA),LSOLAR,FCLC,E011,.FALSE.,
     *               DX,DY,IEL,
     *               TRIA,DCRVSR,
     *               KWORK(L(TRIA(OLVERT))),KWORK(L(TRIA(OLMID))),
     *               DVAL, DVGRX, DVGRY)

C       Transfer the gradient into GRDX/DGRDY and the other DVAL-
C       values to DF/DG:

        DGRX = DVAL(1)
        DGRY = DVAL(2)
        DF   = DVAL(3)
        DG   = DVAL(4)

C       Before moving the point, save the current position

        DXOLD = DX
        DYOLD = DY
        IELOLD = IEL

C       We have the gradient in DGRX/DGRY, the solution in DVAL and the
C       values DF, DG. Build the denominator of the right hand side
C
C                 - 1 / ( s/DF + (1-s)/DG )
C
C       and use this together with the vector field v(x)=(DGRDX,DGRDY)
C       to move the point (DXOLD,DYOLD) to (DX,DY).
      
        DTMP1 = 1D0 / ( DT/DF + (1D0-DT)/DG )
      
        DX = DXOLD + DDT * DGRX*DTMP1
        DY = DYOLD + DDT * DGRY*DTMP1

C       Search in which element the point has been moved to. This must
C       also be done in the last step for boundary nodes to ensure 
C       proper projection with RXBNPT later!
          
        IF (BSRCH4(DX,DY,DCRVSR,
     *             KWORK(L(TRIA(OLVERT))),KWORK(L(TRIA(OLADJ))),
     *             IELOLD,IEL).NE.0) THEN
     
C         Element not found. Set it to 0. May frequently happen on boundary
C         nodes (then we have to set IEL to 0 for RXBNPT), but for inner
C         nodes this is an error!

          IF (IBND.GT.0) THEN
          
            IEL = 0

          ELSE
            IF (MT.GE.2) THEN
              WRITE (MTERM,*) 
     *               'MVGPT1-Warning: Inner node moved to nowhere!'
              WRITE (MTERM,*) 'IVT=',IVT
              WRITE (MTERM,*) 
     *               ' (X,Y)=(',DCRVSR(1,IVT),',',DCRVSR(2,IVT),')'
              WRITE (MTERM,*) 
     *               ' -> (',DX,',',DY,')'
            END IF

            MVGPT1 = 1
            RETURN
            
          END IF
        END IF
     
        IF (IBND.GT.0) THEN
            
C         If we have a boundary node, we have to project the coordinates
C         back onto the boundary. Also determine the new element number if
C         the coordinate felt out of the domain, since there is IEL=0 
C         in that case!

          CALL RXBNPT(TRIA,IBND,IBCT,TMAX(IBCT),IEL,DX,DY,DPRM,IEL)
            
        END IF

      END DO
      
C     The final position is in DX/DY. Save it in the grid.

      DCRVDS(1,IVT) = DX
      DCRVDS(2,IVT) = DY

C     And return the new parameter value in DPARM...
      
      IF (IBND.GT.0) THEN
        DPARM = DPRM
      END IF
      
      MVGPT1 = 0
      
      END
            
************************************************************************
* Move Gridpoint, Method 2 = Explicit Euler on multiple grids
*
* This routine is used in the static grid adaption method. It moves
* one grid-vertex IVT inside of the grid DCRVSR. The resulting 
* coordinates of the point are written into the destination array 
* DCRVDS. The routine is designed only to move (corner-)vertices and
* boundary-vertices, not midpoints!
* The point will be moved along the gradient of the solution vector
* DSOL, which is given by DGRADX, DGRADY.
*
* This routine uses the explicit Euler method with NSTPS steps to move
* the point to its destination coordinates.
*
* The method is the same as in MVGPT1. In contrast to MVGPT1,
* there are slightly different parameters in the call and the routine
* is capable of using different refinement levels for the solution
* of the FE-system and the grid which is deformed!
*
* In:
*  TRIA    - array [1..SZTRIA] of integer
*            Triangulation structure for "source", undeformed 
*            triangulation
*  DCRVSR  - Array DCORVG with source-coordinates of all points in the
*            grid. Must correspond to LCORVG in the TRIA structure.
*  DCRVDS  - array [1..2,1..NVT] of double
*            Array receiving the destination coordinates of the moved
*            gridpoint.
*  IVT     - Number of the node that should be moved; range 1..NVT
*  IBND    - If IVT is a boundary node, IBND is the index of the node
*            inside of the KVBD-array. If IVT is not a boundary node,
*            IBND=0 must be set.
*  IBCT    - If IVT is a boundary node, IBCT must be the number of
*            the boundary component that contains IVT. Otherwise,
*            IBCT=0 must be set.
*  TRIAS   - array [1..SZTRIA,1..*] of integer
*            Triangulation structures of all underlying triangulations
*            that were used to calculate the solution of the linear
*            system, monitor function,...
*  ILVSOL  - Level in TRIAS that corresponds to the data in LSOL,...
*  ILEV    - Level in TRIAS that corresponds to TRIA.
*            TRIA=TRIAS(.,ILEV) is not necessary, but the connectivity
*            of TRIA and TRIAS(.,ILEV) must be the same! (I.e.
*            coordinates of points may differ, but the logical structure
*            must be the same.)
*  LSOL    - Handle to array [1..NEQ] of double precision
*            Solution of the grid adaption system
*  LGRADX  - Handle to array [1..NEQ] of double precision
*            X-gradient of the grid adaption system
*  LGRADY  - Handle to array [1..NEQ] of double precision
*            Y-gradient of the grid adaption system
*  LMON    - Handle to array [1..NEQ] of double precision
*            Monitor function
*  LMNOLD  - Handle to array [1..NEQ] of double precision
*            Old reference function determining the area distribution
*            on the original grid
*
* Method specific input-parameters:
*  NSTPS   - number of steps
*
* Out:
*  DCRVDS [1..2,IVT] - receives the destination coordinates of the
*                      node IVT
*  DPARM             - if IVT is a boundary node, DPARM receives
*                      the new parameter value. Otherwise DPARM is
*                      unchanged.
*  Return value = 0, if successful
*               = 1, if a point was moved but could not be found
*                    in the grid after moving.
************************************************************************

      INTEGER FUNCTION MVGPT2 (TRIA, DCRVSR, DCRVDS, 
     *                   IVT, IBND, IBCT, DPARM,
     *                   TRIAS,ILEV,ILVSOL,
     *                   LSOL, LGRADX, LGRADY,
     *                   LMON, LMNOLD,
     *                   NSTPS)
     
      IMPLICIT NONE
      
      INCLUDE 'stria.inc'
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
 
C parameters
      
      INTEGER TRIA(SZTRIA), IVT, NSTPS, IBND, IBCT
      INTEGER TRIAS(SZTRIA,*), ILVSOL, ILEV
      DOUBLE PRECISION DCRVSR(2,*), DCRVDS(2,*), DPARM
      INTEGER LSOL, LGRADX, LGRADY
      INTEGER LMON, LMNOLD
      
C     local variables
      
      DOUBLE PRECISION DX, DY, DPRM, DXOLD, DYOLD, DT, DDT
      DOUBLE PRECISION DGRX,DGRY,DF, DG, DTMP1
      INTEGER IEL, IELOLD, I, J, KVEL, IELSOL
      
C     externals

      EXTERNAL E011
      INTEGER TNDCNT,BSRCH4
      EXTERNAL TNDCNT,BSRCH4
      DOUBLE PRECISION TMAX
      EXTERNAL TMAX
      
C     Some arrays with constants so that we don't have to set them
C     on runtime. This array is used for the multiple scalar evaluation
C     function. It says that for all handles that we pass to that 
C     routine we only want to calculate function values.

      INTEGER FCLC(4)
      DATA FCLC /1,1,1,1/
      SAVE FCLC
      
C     Furthermore a small array for the multiple scalar evaluation
C     function to pass handles as parameters and obtain return values:

      INTEGER LSOLAR (4)
      DOUBLE PRECISION DVAL(4), DVGRX(4), DVGRY(4)
      
C     Get the initial position of the point and the original element
C     that contains the point (IEL = KVEL(1,NVEL))

      DX = DCRVSR(1,IVT)
      DY = DCRVSR(2,IVT)
      
      KVEL = L(TRIA(OLVEL))
      
      IEL = KWORK(KVEL+TRIA(ONVEL)*(IVT-1))
      
C     Calculate the time-step for the explicit Euler

      DDT = 1D0/DBLE(NSTPS)
      
C     Loop through all time steps. Our first time step is 0, our last
C     is NSTPS-1.

      DO I=0,NSTPS-1
      
C     What's the current time step?

        DT = DBLE(I)*DDT
        
C       We want to solve the following ODE:
C
C           d/dt phi(x,t) = eta ( phi(x,t), t),     t=0..1
C
C       for the right hand side eta:(Omega x [0,1])->Omega,
C
C           eta (y,s) = - v(y) / ( s/f(y) + (1-s)/g(y) )
C
C       and with the initial condition:  phi (x,0) = x
C       Here f corresponds to 1/DMON and g to 1/DMNOLD!
        
C Use the interpolation routine to calculate the gradient of LSOL
C in the current point, based on the reconstructed gradient in
C LGRADX/LGRADY.

C        CALL SCEVQR (TNDCNT(TRIA),DSOL,DGRADX,DGRADY,E011,.FALSE.,
C     *               DX,DY,IEL,
C     *               TRIA,KWORK(L(TRIA(OLVERT))),KWORK(L(TRIA(OLMID))),
C     *               DVAL, DGRX, DGRY)

C Furthermore we need for the denominator the both values 1/DMON
C and 1/DMNOLD in the point (DX,DY). Both functiond DMON and
C DMONOLD are given as function values in the vertices of the old,
C undeformed grid. We interpret the values in the vertices as the
C function values of a piecewise linear FE-function. Therefore to
C obtain a function value in (DX,DY), use the scalar evaluation
C subroutine and the E011 element. The gradient this function returns
C can be ignored:

C        CALL SCEVLQ (TNDCNT(TRIA),DWORK(L(LMON)),E011,.FALSE.,
C     *               DX,DY,IEL,
C     *               TRIA,KWORK(L(TRIA(OLVERT))),KWORK(L(TRIA(OLMID))),
C     *               DF, DTMP1, DTMP1)

C        CALL SCEVLQ (TNDCNT(TRIA),DMNOLD,E011,.FALSE.,
C     *               DX,DY,IEL,
C     *               TRIA,KWORK(L(TRIA(OLVERT))),KWORK(L(TRIA(OLMID))),
C     *               DG, DTMP1, DTMP2)

C       For the real evaluation of these values, we use the "multiple 
C       evaluation function" MSCVQR. this allowes us to calculate multiple
C       values at once, thus increasing the speed of the overall adaption
C       by about 100% !

C       For the multiple solution evaluation function we have to
C       set up some parameters. We calculate:
C         - the X-gradient of DSOL with the recovered gradient vector
C         - the Y-gradient of DSOL with the recovered gradient vector
C         - the function value of DMON giving the value DF=f(x)
C         - the function value of DMNOLF giving the value DG=g(x)
C
C       So set up the handle arrays to pass to MSCVQR:

        LSOLAR(1) = LGRADX
        LSOLAR(2) = LGRADY
        LSOLAR(3) = LMON
        LSOLAR(4) = LMNOLD
        
C       Find the element on the level of the solution of the linear
C       system. If the level of the solution is smaller, we can obtain
C       the element number with the help of the 2-level ordering;
C       otherwise we find it by hierarchical search!
        
        IF (ILVSOL.LT.ILEV) THEN
        
          IELSOL = IEL
          DO J=ILEV-1,ILVSOL,-1
            IF (IELSOL.GT.TRIAS(ONEL,J)) THEN
              IELSOL = (IELSOL-1-TRIAS(ONEL,J)) / 3 + 1
            END IF
          END DO
          
C         Unfortunately, this may not be the element that contains the
C         vertex. Imagine:
C
C         |----|\
C         |    | X  <- X is not in the coarse grid element anymore
C         |----|/
C         
C         Therefore, we start with a search on the coarser
C         level to finally find the element that contains X!
C         We use the hierarchical search on one level only, since
C         that will fall back to line search in case that the element
C         could not be found...

          CALL PSRCH5(TRIAS,ILVSOL,ILVSOL,DX,DY,IELSOL,IELSOL)
        
        ELSE IF (ILVSOL.GT.ILEV) THEN
        
          CALL PSRCH5(TRIAS,ILEV,ILVSOL,DX,DY,IEL,IELSOL)
        
        ELSE 
        
          IELSOL = IEL
          
        END IF
        
C       *What* to calculate is prescribed by the constant array FCLC.
C       This says to evaluate only function values. As the first
C       two solution vectors in LSOLAR point to the recovered gradient
C       of LSOL, the calculated values in DVAL(1) and DVAL(2) are
C       in fact the recovered gradient vector of LSOL!

        CALL MSCVQR (4,TNDCNT(TRIAS(1,ILVSOL)),LSOLAR,FCLC,E011,.FALSE.,
     *               DX,DY,IELSOL,
     *               TRIAS(1,ILVSOL),DWORK(L(TRIAS(OLCORVG,ILVSOL))),
     *               KWORK(L(TRIAS(OLVERT,ILVSOL))),
     *               KWORK(L(TRIAS(OLMID,ILVSOL))),
     *               DVAL, DVGRX, DVGRY)

C       Transfer the gradient into GRDX/DGRDY and the other DVAL-
C       values to DF/DG:

        DGRX = DVAL(1)
        DGRY = DVAL(2)
        DF   = DVAL(3)
        DG   = DVAL(4)

C       Before moving the point, save the current position

        DXOLD = DX
        DYOLD = DY
        IELOLD = IEL

C       We have the gradient in DGRX/DGRY, the solution in DVAL and the
C       values DF, DG. Build the denominator of the right hand side
C
C                 - 1 / ( s/DF + (1-s)/DG )
C
C       and use this together with the vector field v(x)=(DGRDX,DGRDY)
C       to move the point (DXOLD,DYOLD) to (DX,DY).
      
        DTMP1 = 1D0 / ( DT/DF + (1D0-DT)/DG )
      
        DX = DXOLD + DDT * DGRX*DTMP1
        DY = DYOLD + DDT * DGRY*DTMP1

C       Search in which element the point has been moved to. This must
C       also be done in the last step for boundary nodes to ensure 
C       proper projection with RXBNPT later!
          
        IF (BSRCH4(DX,DY,DCRVSR,
     *             KWORK(L(TRIA(OLVERT))),KWORK(L(TRIA(OLADJ))),
     *             IELOLD,IEL).NE.0) THEN
     
C         Element not found. Set it to 0. May frequently happen on boundary
C         nodes (then we have to set IEL to 0 for RXBNPT), but for inner
C         nodes this is an error!

          IF (IBND.GT.0) THEN
          
            IEL = 0

          ELSE
            IF (MT.GE.2) THEN
              WRITE (MTERM,*) 
     *               'MVGPT2-Warning: Inner node moved to nowhere!'
              WRITE (MTERM,*) 'IVT=',IVT
              WRITE (MTERM,*) 
     *               ' (X,Y)=(',DCRVSR(1,IVT),',',DCRVSR(2,IVT),')'
              WRITE (MTERM,*) 
     *               ' -> (',DX,',',DY,')'
            END IF

            MVGPT2 = 1
            RETURN
            
          END IF
        END IF
     
        IF (IBND.GT.0) THEN
            
C         If we have a boundary node, we have to project the coordinates
C         back onto the boundary. Also determine the new element number if
C         the coordinate felt out of the domain, since there is IEL=0 
C         in that case!

          CALL RXBNPT(TRIA,IBND,IBCT,TMAX(IBCT),IEL,DX,DY,DPRM,IEL)
            
        END IF

      END DO
      
C     The final position is in DX/DY. Save it in the grid.

      DCRVDS(1,IVT) = DX
      DCRVDS(2,IVT) = DY

C     And return the new parameter value in DPARM...
      
      IF (IBND.GT.0) THEN
        DPARM = DPRM
      END IF
      
      MVGPT2 = 0
      
      END
            
************************************************************************
* Move Gridpoint, Method 3 = Runge-Kutta 3
*
* This routine is used in the static grid adaption method. It moves
* one grid-vertex IVT inside of the grid DCRVSR. The resulting 
* coordinates of the point are written into the destination array 
* DCRVDS. The routine is designed only to move (corner-)vertices and
* boundary-vertices, not midpoints!
* The point will be moved along the gradient of the solution vector
* DSOL, which is given by DGRADX, DGRADY.
*
* This routine uses the explicit RK3 method with NSTPS steps to move
* the point to its destination coordinates.
*
* In:
*  TRIA    - array [1..SZTRIA] of integer
*            Triangulation structure for "source", undeformed 
*            triangulation
*  DCRVSR  - Array DCORVG with source-coordinates of all points in the
*            grid. Must correspond to LCORVG in the TRIA structure.
*  DCRVDS  - array [1..2,1..NVT] of double
*            Array receiving the destination coordinates of the moved
*            gridpoint.
*  IVT     - Number of the node that should be moved; range 1..NVT
*  IBND    - If IVT is a boundary node, IBND is the index of the node
*            inside of the KVBD-array. If IVT is not a boundary node,
*            IBND=0 must be set.
*  IBCT    - If IVT is a boundary node, IBCT must be the number of
*            the boundary component that contains IVT. Otherwise,
*            IBCT=0 must be set.
*  LSOL    - Handle to array [1..NEQ] of double precision
*            Solution of the grid adaption system
*  LGRADX  - Handle to array [1..NEQ] of double precision
*            X-gradient of the grid adaption system
*  LGRADY  - Handle to array [1..NEQ] of double precision
*            Y-gradient of the grid adaption system
*  LMON    - Handle to array [1..NEQ] of double precision
*            Monitor function
*  LMNOLD  - Handle to array [1..NEQ] of double precision
*            Old reference function determining the area distribution
*            on the original grid
*
* Method specific input-parameters:
*  NSTPS   - number of steps
*
* Out:
*  DCRVDS [1..2,IVT] - receives the destination coordinates of the
*                      node IVT
*  DPARM             - if IVT is a boundary node, DPARM receives
*                      the new parameter value. Otherwise DPARM is
*                      unchanged.
*  Return value = 0, if successful
*               = 1, if a point was moved but could not be found
*                    in the grid after moving.
************************************************************************

      INTEGER FUNCTION MVGPT3 (TRIA, DCRVSR, DCRVDS, 
     *                   IVT, IBND, IBCT, DPARM,
     *                   LSOL, LGRADX, LGRADY,
     *                   LMON, LMNOLD,
     *                   NSTPS)
     
      IMPLICIT NONE
      
      INCLUDE 'stria.inc'
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
 
C     parameters
      
      INTEGER TRIA(SZTRIA), IVT, NSTPS, IBND, IBCT
      DOUBLE PRECISION DCRVSR(2,*), DCRVDS(2,*), DPARM
      INTEGER LSOL, LGRADX, LGRADY, KVEL
      INTEGER LMON, LMNOLD
      
C     local variables
      
      DOUBLE PRECISION DX, DY, DPRM, DXOLD, DYOLD, DT, DDT
      DOUBLE PRECISION DGRX,DGRY,DF   , DG   , DTMP1, DTMP2, K(2,3)
      DOUBLE PRECISION DKX(2), DKY(2)
      INTEGER IEL, IELOLD, I, IELK(2), IBC
      
C     externals

      EXTERNAL E011
      INTEGER TNDCNT,BSRCH4
      EXTERNAL TNDCNT,BSRCH4
      DOUBLE PRECISION TMAX
      EXTERNAL TMAX
      
C     Some arrays with constants so that we don't have to set them
C     on runtime. This array is used for the multiple scalar evaluation
C     function. It says that for all handles that we pass to that
C     routine we only want to calculate function values.

      INTEGER FCLC(8)
      DATA FCLC /1,1,1,1, 1,1, 1,1/
      SAVE FCLC
      
C     Furthermore a small array for the multiple scalar evaluation
C     function to pass handles as parameters and obtain return values:

      INTEGER LSOLAR (8)
      DOUBLE PRECISION DVAL(8), DVGRX(8), DVGRY(8)
      
C     Get the initial position of the point and the original element
C     that contains the point (IEL = KVEL(1,NVEL))

      DX = DCRVSR(1,IVT)
      DY = DCRVSR(2,IVT)
      
      KVEL = L(TRIA(OLVEL))
      
      IEL = KWORK(KVEL+TRIA(ONVEL)*(IVT-1))
      
C     Calculate the time-step for the explicit Euler

      DDT = 1D0/DBLE(NSTPS)
      
C     Loop through all time steps. Because we are looking "ahead", our 
C     first time step is 0 and the last time step is NSTPS-1.

      DO I=0,NSTPS-1
      
C       What's the current time step?

        DT = DBLE(I)*DDT
        
C       We want to solve the following ODE:
C
C           d/dt phi(x,t) = eta ( phi(x,t), t),     t=0..1
C
C       for the right hand side eta:(Omega x [0,1])->Omega,
C
C           eta (y,s) = - v(y) / ( s/f(y) + (1-s)/g(y) )
C
C       and with the initial condition:  phi (x,0) = x
C       Here f corresponds to 1/DMON and g to 1/DMNOLD!
              
C       Use the interpolation routine to calculate the gradient of LSOL
C       in the current point, based on the reconstructed gradient in
C       LGRADX/LGRADY.

        LSOLAR(1) = LGRADX
        LSOLAR(2) = LGRADY

C       Furthermore we need for the denominator the both values 1/DMON
C       and 1/DMNOLD in the point (DX,DY). Both functiond DMON and
C       DMONOLD are given as function values in the vertices of the old,
C       undeformed grid. We interpret the values in the vertices as the
C       function values of a piecewise linear FE-function. Therefore to
C       obtain a function value in (DX,DY), use the scalar evaluation
C       subroutine and the E011 element. 

        LSOLAR(3) = LMON
        LSOLAR(4) = LMNOLD

        LSOLAR(5) = LMON
        LSOLAR(6) = LMNOLD

        LSOLAR(7) = LMON
        LSOLAR(8) = LMNOLD
        
C For the real evaluation of these values, we use the "multiple 
C evaluation function" MSCVQR. this allowes us to calculate multiple
C values at once, thus increasing the speed of the overall adaption
C by about 100% !

C       For the multiple solution evaluation function we have 
C       set up some parameters above. We prepared to calculate:
C         - the X-gradient of DSOL with the recovered gradient vector
C         - the Y-gradient of DSOL with the recovered gradient vector
C         - the function value of DMON giving the value DF=f(x)
C         - the function value of DMNOLF giving the value DG=g(x)
C
C       So set up the handle arrays to pass to MSCVQR:

C       *What* to calculate is prescribed by the constant array FCLC.
C       This says to evaluate only function values. As the first
C       two solution vectors in LSOLAR point to the recovered gradient
C       of LSOL, the calculated values in DVAL(1) and DVAL(2) are
C       in fact the recovered gradient vector of LSOL!
C       We call MSCVQR to evaluate DGRADX,DGRADY,DMON and DMNOLD
C       in the current point (DX,DY):

        CALL MSCVQR (4,TNDCNT(TRIA),LSOLAR,FCLC,E011,.FALSE.,
     *               DX,DY,IEL,
     *               TRIA,DCRVSR,
     *               KWORK(L(TRIA(OLVERT))),KWORK(L(TRIA(OLMID))),
     *               DVAL, DVGRX, DVGRY)

C       Transfer the gradient into GRDX/DGRDY and the other DVAL-
C       values to DF/DG:

        DGRX = DVAL(1)
        DGRY = DVAL(2)
        DF   = DVAL(3)
        DG   = DVAL(4)
        
C       Remember: for Runge-Kutta-3 we need more!
C       The RK3-Formula was (for y_n being the X- and/or the 
C       Y-coordinate of our point (DX,DY) ):
C
C       y_n = y_(n-1) + 1/6 * h * (k1 + 4*k2 + k3)
C
C       With:
C             k1 = f( t_(n-1)   , y(n-1) )
C             k2 = f( t_(n-1/2) , y_(n-1) + 1/2*h*k1 )
C             k3 = f( t_n       , y_(n-1) - h*k1 + 2*h*k2 )
C
C       With the values calculated in DVAL(1..4) we are able to
C       construct k1 for the X- and Y-coordinate:
C
C       We have the gradient in DGRX/DGRY and the
C       values DF, DG. Build the denominator of the right hand side
C
C           - 1 / ( s/DF + (1-s)/DG )
C
C       and use this together with the vector field v(x)=(DGRDX,DGRDY)
C       to build that k1. K(1,1) is the k1 of the 1D-ODE for the X-
C       coordinate, K(2,1) the k1 of the 1D-ODE for the Y-coordinate -
C       since our 2D-system decomposes into two 1D ODE's:
      
        DTMP1 = 1D0 / ( DT/DF + (1D0-DT)/DG )
      
        K(1,1) = DGRX*DTMP1
        K(2,1) = DGRY*DTMP1
        
C       Now we have to calculate k2, i.e. we again have to build the 
C       denominator and build K(1,2) and K(2,2). Changing the time
C       by 1/2*DDT is not the problem, but changing the point is!
C       Because by changing the point we might left our current element!
C
C       Therefore before doing anything else we calculate in advance
C       the elements that contain the spatial position that is necessary
C       for k2. We have the following coordinates:

        DKX(1) = DX + 0.5D0 * DDT * K(1,1)
        DKY(1) = DY + 0.5D0 * DDT * K(2,1)
        
C       Use the search routine to determine the element of that point:

        IF (BSRCH4(DKX(1),DKY(1),DCRVSR,
     *                KWORK(L(TRIA(OLVERT))),KWORK(L(TRIA(OLADJ))),
     *                IEL,IELK(1)).NE.0) THEN

C         Element not found. Set it to 0. May happen on boundary nodes 
C         (then we have to set IEL to 0 for RXBNPT), but should not 
C         happen for inner nodes!
C         If this happens for a boundary node, use the boundary 
C         projection routine to recalculate the current element.

          IF (IBND.GT.0) THEN
            IELK(1) = 0
            CALL RXBNPT(TRIA,IBND,IBCT,TMAX(IBCT),IELK(1),DKX(1),DKY(1),
     *                  DPRM,IELK(1))
          ELSE
            IF (MT.GE.2) THEN
              WRITE (MTERM,*) 'MVGPT3-Warning: RK3 substep 1: '//
     *        ' Inner node moved to nowhere!'
              WRITE (MTERM,*) 'IVT=',IVT
              WRITE (MTERM,*) 
     *               ' (X,Y)=(',DCRVSR(1,IVT),',',DCRVSR(2,IVT),')'
              WRITE (MTERM,*) ' -> (',DKX(1),',',DKY(1),')'
            END IF
              
C            MVGPT3 = 1
C            RETURN
 
C           It's only a warning, we don't halt the program here.
C           This case should normally not occur. But since we computed
C           only an intermediate point, we don't want to treat this as 
C           an error here. Instead, if the point left the domain,
C           we project it back and continue working with the
C           projected point on the boundary. this is not full order
C           anymore, but perhaps it helps to repair a nearly
C           corrupted mesh...

            IELK(1) = 0
            IBC = 0
            CALL RXBNPT(TRIA,IBND,IBC,TMAX(IBC),IELK(1),DKX(1),DKY(1),
     *                  DPRM,IELK(1))
            
          END IF
        END IF

C       Again use the multiple-evaluation routine to calculate the
C       values of DF and DG - this time in the new point. This time
C       save the results in DVAL (5..6). DVGRX/DVGRY are
C       passed as dummy-parameters.
     
        CALL MSCVQR (2,TNDCNT(TRIA),LSOLAR(5),FCLC(5),E011,.FALSE.,
     *               DKX(1),DKY(1),IELK(1),
     *               TRIA,DCRVSR,
     *               KWORK(L(TRIA(OLVERT))),KWORK(L(TRIA(OLMID))),
     *               DVAL(5), DVGRX(5), DVGRY(5))
     
C       Calculate the denominator of the RHS and with that K2.
C       Remember that here we evaluate at the time step DT_(n-1/3)!

        DF   = DVAL(5)
        DG   = DVAL(6)

        DTMP2 = DT + 0.5D0 * DDT
        DTMP1 = 1D0 / ( DTMP2/DF + (1D0-DTMP2)/DG )
      
        K(1,2) = DGRX*DTMP1
        K(2,2) = DGRY*DTMP1
        
C       Now use k1 and k2 to calculate the spacial coordinates of the
C       point to evaluate for k3:
        
        DKX(2) = DX - DDT * K(1,1) + 2D0*DDT*K(1,2)
        DKY(2) = DY - DDT * K(2,1) + 2D0*DDT*K(2,2)

C       Again use the search routine to determine the element of the
C       new spacial point. Again start searching from element IEL:

        IF (BSRCH4(DKX(2),DKY(2),DCRVSR,
     *                KWORK(L(TRIA(OLVERT))),KWORK(L(TRIA(OLADJ))),
     *                IEL,IELK(2)).NE.0) THEN

C         Element not found. Set it to 0. May happen on boundary nodes 
C         (then we have to set IEL to 0 for RXBNPT), but should not 
C         happen for inner nodes!
C         If this happens for a boundary node, use the boundary 
C         projection routine to recalculate the current element.

          IF (IBND.GT.0) THEN
            IELK(2) = 0
            CALL RXBNPT(TRIA,IBND,IBCT,TMAX(IBCT),IELK(2),DKX(2),DKY(2),
     *                  DPRM,IELK(2))
          ELSE
            IF (MT.GE.2) THEN
              WRITE (MTERM,*) 'MVGPT3-Warning: RK3 substep 2: '//
     *        ' Inner node moved to nowhere!'
              WRITE (MTERM,*) 'IVT=',IVT
              WRITE (MTERM,*) 
     *               ' (X,Y)=(',DCRVSR(1,IVT),',',DCRVSR(2,IVT),')'
              WRITE (MTERM,*) ' -> (',DKY(2),',',DKY(2),')'
            END IF
              
C            MVGPT3 = 1
C            RETURN

C           It's only a warning, we don't halt the program here.
C           This case should normally not occur. But since we computed
C           only an intermediate point, we don't want to treat this as 
C           an error here. Instead, if the point left the domain,
C           we project it back and continue working with the
C           projected point on the boundary. this is not full order
C           anymore, but perhaps it helps to repair a nearly
C           corrupted mesh...

            IELK(2) = 0
            IBC = 0
            CALL RXBNPT(TRIA,IBND,IBC,TMAX(IBC),IELK(2),DKX(2),DKY(2),
     *                  DPRM,IELK(2))
            
          END IF
        END IF

C       Again use the multiple-evaluation routine to calculate the
C       values of DF and DG - this time in the new point. This time
C       save the results in DVAL(7..8).
     
        CALL MSCVQR (2,TNDCNT(TRIA),LSOLAR(7),FCLC(7),E011,.FALSE.,
     *               DKX(2),DKY(2),IELK(2),
     *               TRIA,DCRVSR,
     *               KWORK(L(TRIA(OLVERT))),KWORK(L(TRIA(OLMID))),
     *               DVAL(7), DVGRX(7), DVGRY(7))
     
C       Calculate the denominator of the RHS and with that K3.
C       Remember that here we evaluate at the time step DT_n!

        DF   = DVAL(7)
        DG   = DVAL(8)

        DTMP2 = DT + DDT
        DTMP1 = 1D0 / ( DTMP2/DF + (1D0-DTMP2)/DG )
      
        K(1,3) = DGRX*DTMP1
        K(2,3) = DGRY*DTMP1

C       Before moving the point, save the current position

        DXOLD = DX
        DYOLD = DY
        IELOLD = IEL

C       Finally move the point according to the RK-3 formula:

        DX = DXOLD + (1D0/6D0) * DDT * (K(1,1) + 4D0*K(1,2) + K(1,3))
        DY = DYOLD + (1D0/6D0) * DDT * (K(2,1) + 4D0*K(2,2) + K(2,3))

C       Search in which element the point has been moved to. This must also
C       be done in the last step for boundary nodes to ensure proper 
C       projection with RXBNPT later!
          
        IF (BSRCH4(DX,DY,DCRVSR,
     *             KWORK(L(TRIA(OLVERT))),KWORK(L(TRIA(OLADJ))),
     *             IELOLD,IEL).NE.0) THEN
     
C         Element not found. Set it to 0. May frequently happen on boundary
C         nodes (then we have to set IEL to 0 for RXBNPT), but for inner
C         nodes this is an error!

          IF (IBND.GT.0) THEN

            IEL = 0

C           If we have a boundary node, we have to project the coordinates
C           back onto the boundary. Also determine the new element number of
C           the node in (X,Y), since there is IEL=0 at the moment!

            CALL RXBNPT(TRIA,IBND,IBCT,TMAX(IBCT),IEL,DX,DY,DPRM,IEL)
            
          ELSE
            IF (MT.GE.2) THEN
              WRITE (MTERM,*) 
     *               'MVGPT3-Warning: Inner node moved to nowhere!'
              WRITE (MTERM,*) 'IVT=',IVT
              WRITE (MTERM,*) 
     *               ' (X,Y)=(',DCRVSR(1,IVT),',',DCRVSR(2,IVT),')'
              WRITE (MTERM,*) ' -> (',DX,',',DY,')'
            END IF
            
            MVGPT3 = 1
            RETURN
            
          END IF
        END IF
     
        IF (IBND.GT.0) THEN

C         If we have a boundary node, we have to project the coordinates
C         back onto the boundary. Also determine the new element number if
C         the coordinate felt out of the domain, since there is IEL=0 
C         in that case!

          CALL RXBNPT(TRIA,IBND,IBCT,TMAX(IBCT),IEL,DX,DY,DPRM,IEL)
            
        END IF
     
      END DO
      
C     The final position is in DX/DY. Save it in the grid.

      DCRVDS(1,IVT) = DX
      DCRVDS(2,IVT) = DY
      
C     And return the new parameter value in DPARM...
      
      IF (IBND.GT.0) THEN
        DPARM = DPRM
      END IF
      
      MVGPT3 = 0
      
      END
            
************************************************************************
* Move Gridpoint, Method 4 = Runge-Kutta 3 on multiple levels
*
* This routine is used in the static grid adaption method. It moves
* one grid-vertex IVT inside of the grid DCRVSR. The resulting 
* coordinates of the point are written into the destination array 
* DCRVDS. The routine is designed only to move (corner-)vertices and
* boundary-vertices, not midpoints!
* The point will be moved along the gradient of the solution vector
* DSOL, which is given by DGRADX, DGRADY.
*
* This routine uses the explicit RK3 method with NSTPS steps to move
* the point to its destination coordinates.
*
* The method is the same as in MVGPT3. In contrast to MVGPT3,
* there are slightly different parameters in the call and the routine
* is capable of using different refinement levels for the solution
* of the FE-system and the grid which is deformed!
*
* In:
*  TRIA    - array [1..SZTRIA] of integer
*            Triangulation structure for "source", undeformed 
*            triangulation
*  DCRVSR  - Array DCORVG with source-coordinates of all points in the
*            grid. Must correspond to LCORVG in the TRIA structure.
*  DCRVDS  - array [1..2,1..NVT] of double
*            Array receiving the destination coordinates of the moved
*            gridpoint.
*  IVT     - Number of the node that should be moved; range 1..NVT
*  IBND    - If IVT is a boundary node, IBND is the index of the node
*            inside of the KVBD-array. If IVT is not a boundary node,
*            IBND=0 must be set.
*  IBCT    - If IVT is a boundary node, IBCT must be the number of
*            the boundary component that contains IVT. Otherwise,
*            IBCT=0 must be set.
*  TRIAS   - array [1..SZTRIA,1..*] of integer
*            Triangulation structures of all underlying triangulations
*            that were used to calculate the solution of the linear
*            system, monitor function,...
*  ILVSOL  - Level in TRIAS that corresponds to the data in LSOL,...
*  ILEV    - Level in TRIAS that corresponds to TRIA.
*            TRIA=TRIAS(.,ILEV) is not necessary, but the connectivity
*            of TRIA and TRIAS(.,ILEV) must be the same! (I.e.
*            coordinates of points may differ, but the logical structure
*            must be the same.)
*  LSOL    - Handle to array [1..NEQ] of double precision
*            Solution of the grid adaption system
*  LGRADX  - Handle to array [1..NEQ] of double precision
*            X-gradient of the grid adaption system
*  LGRADY  - Handle to array [1..NEQ] of double precision
*            Y-gradient of the grid adaption system
*  LMON    - Handle to array [1..NEQ] of double precision
*            Monitor function
*  LMNOLD  - Handle to array [1..NEQ] of double precision
*            Old reference function determining the area distribution
*            on the original grid
*
* Method specific input-parameters:
*  NSTPS   - number of steps
*
* Out:
*  DCRVDS [1..2,IVT] - receives the destination coordinates of the
*                      node IVT
*  DPARM             - if IVT is a boundary node, DPARM receives
*                      the new parameter value. Otherwise DPARM is
*                      unchanged.
*  Return value = 0, if successful
*               = 1, if a point was moved but could not be found
*                    in the grid after moving.
************************************************************************

      INTEGER FUNCTION MVGPT4 (TRIA, DCRVSR, DCRVDS, 
     *                   IVT, IBND, IBCT, DPARM,
     *                   TRIAS,ILEV,ILVSOL,
     *                   LSOL, LGRADX, LGRADY,
     *                   LMON, LMNOLD,
     *                   NSTPS)
     
      IMPLICIT NONE
      
      INCLUDE 'stria.inc'
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
 
C     parameters
      
      INTEGER TRIA(SZTRIA), IVT, NSTPS, IBND, IBCT
      INTEGER TRIAS(SZTRIA,*),ILEV,ILVSOL
      DOUBLE PRECISION DCRVSR(2,*), DCRVDS(2,*), DPARM
      INTEGER LSOL, LGRADX, LGRADY
      INTEGER LMON, LMNOLD
      
C     local variables
      
      DOUBLE PRECISION DX, DY, DPRM, DXOLD, DYOLD, DT, DDT
      DOUBLE PRECISION DGRX,DGRY,DF   , DG   , DTMP1, DTMP2, K(2,3)
      DOUBLE PRECISION DKX(2), DKY(2)
      INTEGER IEL, IELOLD, I, J, IELK(2), IBC, IELSOL, KVEL
      
C     externals

      EXTERNAL E011
      INTEGER TNDCNT,BSRCH4
      EXTERNAL TNDCNT,BSRCH4
      DOUBLE PRECISION TMAX
      EXTERNAL TMAX
      
C     Some arrays with constants so that we don't have to set them
C     on runtime. This array is used for the multiple scalar evaluation
C     function. It says that for all handles that we pass to that
C     routine we only want to calculate function values.

      INTEGER FCLC(8)
      DATA FCLC /1,1,1,1, 1,1, 1,1/
      SAVE FCLC
      
C     Furthermore a small array for the multiple scalar evaluation
C     function to pass handles as parameters and obtain return values:

      INTEGER LSOLAR (8)
      DOUBLE PRECISION DVAL(8), DVGRX(8), DVGRY(8)
      
C     Get the initial position of the point and the original element
C     that contains the point (IEL = KVEL(1,NVEL))

      DX = DCRVSR(1,IVT)
      DY = DCRVSR(2,IVT)
      
      KVEL = L(TRIA(OLVEL))
      
      IEL = KWORK(KVEL+TRIA(ONVEL)*(IVT-1))
      
C     Calculate the time-step for the explicit Euler

      DDT = 1D0/DBLE(NSTPS)
      
C     Loop through all time steps. Because we are looking "ahead", our 
C     first time step is 0 and the last time step is NSTPS-1.

      DO I=0,NSTPS-1
      
C       What's the current time step?

        DT = DBLE(I)*DDT
        
C       We want to solve the following ODE:
C
C           d/dt phi(x,t) = eta ( phi(x,t), t),     t=0..1
C
C       for the right hand side eta:(Omega x [0,1])->Omega,
C
C           eta (y,s) = - v(y) / ( s/f(y) + (1-s)/g(y) )
C
C       and with the initial condition:  phi (x,0) = x
C       Here f corresponds to 1/DMON and g to 1/DMNOLD!
              
C       Use the interpolation routine to calculate the gradient of LSOL
C       in the current point, based on the reconstructed gradient in
C       LGRADX/LGRADY.

        LSOLAR(1) = LGRADX
        LSOLAR(2) = LGRADY

C       Furthermore we need for the denominator the both values 1/DMON
C       and 1/DMNOLD in the point (DX,DY). Both functiond DMON and
C       DMONOLD are given as function values in the vertices of the old,
C       undeformed grid. We interpret the values in the vertices as the
C       function values of a piecewise linear FE-function. Therefore to
C       obtain a function value in (DX,DY), use the scalar evaluation
C       subroutine and the E011 element. 

        LSOLAR(3) = LMON
        LSOLAR(4) = LMNOLD

        LSOLAR(5) = LMON
        LSOLAR(6) = LMNOLD

        LSOLAR(7) = LMON
        LSOLAR(8) = LMNOLD
        
C       Find the element on the level of the solution of the linear
C       system. If the level of the solution is smaller, we can obtain
C       the element number with the help of the 2-level ordering;
C       otherwise we find it by hierarchical search!
        
        IF (ILVSOL.LT.ILEV) THEN
        
          IELSOL = IEL
          DO J=ILEV-1,ILVSOL,-1
            IF (IELSOL.GT.TRIAS(ONEL,J)) THEN
              IELSOL = (IELSOL-1-TRIAS(ONEL,J)) / 3 + 1
            END IF
          END DO
          
C         Unfortunately, this may not be the element that contains the
C         vertex. Imagine:
C
C         |----|\
C         |    | X  <- X is not in the coarse grid element anymore
C         |----|/
C         
C         Therefore, we start with a search on the coarser
C         level to finally find the element that contains X!
C         We use the hierarchical search on one level only, since
C         that will fall back to line search in case that the element
C         could not be found...

          CALL PSRCH5(TRIAS,ILVSOL,ILVSOL,DX,DY,IELSOL,IELSOL)
        
        ELSE IF (ILVSOL.GT.ILEV) THEN
        
          CALL PSRCH5(TRIAS,ILEV,ILVSOL,DX,DY,IEL,IELSOL)
        
        ELSE 
        
          IELSOL = IEL
          
        END IF
        
C For the real evaluation of these values, we use the "multiple 
C evaluation function" MSCVQR. this allowes us to calculate multiple
C values at once, thus increasing the speed of the overall adaption
C by about 100% !

C       For the multiple solution evaluation function we have 
C       set up some parameters above. We prepared to calculate:
C         - the X-gradient of DSOL with the recovered gradient vector
C         - the Y-gradient of DSOL with the recovered gradient vector
C         - the function value of DMON giving the value DF=f(x)
C         - the function value of DMNOLF giving the value DG=g(x)
C
C       So set up the handle arrays to pass to MSCVQR:

C       *What* to calculate is prescribed by the constant array FCLC.
C       This says to evaluate only function values. As the first
C       two solution vectors in LSOLAR point to the recovered gradient
C       of LSOL, the calculated values in DVAL(1) and DVAL(2) are
C       in fact the recovered gradient vector of LSOL!
C       We call MSCVQR to evaluate DGRADX,DGRADY,DMON and DMNOLD
C       in the current point (DX,DY):

        CALL MSCVQR (4,TNDCNT(TRIAS(1,ILVSOL)),LSOLAR,FCLC,E011,.FALSE.,
     *               DX,DY,IELSOL,
     *               TRIAS(1,ILVSOL),DWORK(L(TRIAS(OLCORVG,ILVSOL))),
     *               KWORK(L(TRIAS(OLVERT,ILVSOL))),
     *               KWORK(L(TRIAS(OLMID,ILVSOL))),
     *               DVAL, DVGRX, DVGRY)

C       Transfer the gradient into GRDX/DGRDY and the other DVAL-
C       values to DF/DG:

        DGRX = DVAL(1)
        DGRY = DVAL(2)
        DF   = DVAL(3)
        DG   = DVAL(4)
        
C       Remember: for Runge-Kutta-3 we need more!
C       The RK3-Formula was (for y_n being the X- and/or the 
C       Y-coordinate of our point (DX,DY) ):
C
C       y_n = y_(n-1) + 1/6 * h * (k1 + 4*k2 + k3)
C
C       With:
C             k1 = f( t_(n-1)   , y(n-1) )
C             k2 = f( t_(n-1/2) , y_(n-1) + 1/2*h*k1 )
C             k3 = f( t_n       , y_(n-1) - h*k1 + 2*h*k2 )
C
C       With the values calculated in DVAL(1..4) we are able to
C       construct k1 for the X- and Y-coordinate:
C
C       We have the gradient in DGRX/DGRY and the
C       values DF, DG. Build the denominator of the right hand side
C
C           - 1 / ( s/DF + (1-s)/DG )
C
C       and use this together with the vector field v(x)=(DGRDX,DGRDY)
C       to build that k1. K(1,1) is the k1 of the 1D-ODE for the X-
C       coordinate, K(2,1) the k1 of the 1D-ODE for the Y-coordinate -
C       since our 2D-system decomposes into two 1D ODE's:
      
        DTMP1 = 1D0 / ( DT/DF + (1D0-DT)/DG )
      
        K(1,1) = DGRX*DTMP1
        K(2,1) = DGRY*DTMP1
        
C       Now we have to calculate k2, i.e. we again have to build the 
C       denominator and build K(1,2) and K(2,2). Changing the time
C       by 1/2*DDT is not the problem, but changing the point is!
C       Because by changing the point we might left our current element!
C
C       Therefore before doing anything else we calculate in advance
C       the elements that contain the spatial position that is necessary
C       for k2. We have the following coordinates:

        DKX(1) = DX + 0.5D0 * DDT * K(1,1)
        DKY(1) = DY + 0.5D0 * DDT * K(2,1)
        
C       Use the search routine to determine the element of that point:

        IF (BSRCH4(DKX(1),DKY(1),DCRVSR,
     *                KWORK(L(TRIA(OLVERT))),KWORK(L(TRIA(OLADJ))),
     *                IEL,IELK(1)).NE.0) THEN

C         Element not found. Set it to 0. May happen on boundary nodes 
C         (then we have to set IEL to 0 for RXBNPT), but should not 
C         happen for inner nodes!
C         If this happens for a boundary node, use the boundary 
C         projection routine to recalculate the current element.

          IF (IBND.GT.0) THEN
            IELK(1) = 0
            CALL RXBNPT(TRIA,IBND,IBCT,TMAX(IBCT),IELK(1),DKX(1),DKY(1),
     *                  DPRM,IELK(1))
          ELSE
            IF (MT.GE.2) THEN
              WRITE (MTERM,*) 'MVGPT3-Warning: RK3 substep 1: '//
     *        ' Inner node moved to nowhere!'
              WRITE (MTERM,*) 'IVT=',IVT
              WRITE (MTERM,*) 
     *               ' (X,Y)=(',DCRVSR(1,IVT),',',DCRVSR(2,IVT),')'
              WRITE (MTERM,*) ' -> (',DKX(1),',',DKY(1),')'
            END IF
              
C            MVGPT3 = 1
C            RETURN
 
C           It's only a warning, we don't halt the program here.
C           This case should normally not occur. But since we computed
C           only an intermediate point, we don't want to treat this as 
C           an error here. Instead, if the point left the domain,
C           we project it back and continue working with the
C           projected point on the boundary. this is not full order
C           anymore, but perhaps it helps to repair a nearly
C           corrupted mesh...

            IELK(1) = 0
            IBC = 0
            CALL RXBNPT(TRIA,IBND,IBC,TMAX(IBC),IELK(1),DKX(1),DKY(1),
     *                  DPRM,IELK(1))
            
          END IF
        END IF

C       Find the element on the level of the solution of the linear
C       system. If the level of the solution is smaller, we can obtain
C       the element number with the help of the 2-level ordering;
C       otherwise we find it by hierarchical search!
        
        IF (ILVSOL.LT.ILEV) THEN
        
          IELSOL = IELK(1)
          DO J=ILEV-1,ILVSOL,-1
            IF (IELSOL.GT.TRIAS(ONEL,J)) THEN
              IELSOL = (IELSOL-1-TRIAS(ONEL,J)) / 3 + 1
            END IF
          END DO
          
C         Unfortunately, this may not be the element that contains the
C         vertex. Imagine:
C
C         |----|\
C         |    | X  <- X is not in the coarse grid element anymore
C         |----|/
C         
C         Therefore, we start with a search on the coarser
C         level to finally find the element that contains X!
C         We use the hierarchical search on one level only, since
C         that will fall back to line search in case that the element
C         could not be found...

          CALL PSRCH5(TRIAS,ILVSOL,ILVSOL,DKX(1),DKY(1),IELSOL,IELSOL)
        
        ELSE IF (ILVSOL.GT.ILEV) THEN
        
          CALL PSRCH5(TRIAS,ILEV,ILVSOL,DKX(1),DKY(1),IELK(1),IELSOL)
        
        ELSE 
        
          IELSOL = IELK(1)
          
        END IF
        
C       Again use the multiple-evaluation routine to calculate the
C       values of DF and DG - this time in the new point. This time
C       save the results in DVAL (5..6). DVGRX/DVGRY are
C       passed as dummy-parameters.
     
        CALL MSCVQR (2,TNDCNT(TRIAS(1,ILVSOL)),LSOLAR(5),FCLC(5),
     *               E011,.FALSE.,
     *               DKX(1),DKY(1),IELSOL,
     *               TRIAS(1,ILVSOL),DWORK(L(TRIAS(OLCORVG,ILVSOL))),
     *               KWORK(L(TRIAS(OLVERT,ILVSOL))),
     *               KWORK(L(TRIAS(OLMID,ILVSOL))),
     *               DVAL(5), DVGRX(5), DVGRY(5))
     
C       Calculate the denominator of the RHS and with that K2.
C       Remember that here we evaluate at the time step DT_(n-1/3)!

        DF   = DVAL(5)
        DG   = DVAL(6)

        DTMP2 = DT + 0.5D0 * DDT
        DTMP1 = 1D0 / ( DTMP2/DF + (1D0-DTMP2)/DG )
      
        K(1,2) = DGRX*DTMP1
        K(2,2) = DGRY*DTMP1
        
C       Now use k1 and k2 to calculate the spacial coordinates of the
C       point to evaluate for k3:
        
        DKX(2) = DX - DDT * K(1,1) + 2D0*DDT*K(1,2)
        DKY(2) = DY - DDT * K(2,1) + 2D0*DDT*K(2,2)

C       Again use the search routine to determine the element of the
C       new spacial point. Again start searching from element IEL:

        IF (BSRCH4(DKX(2),DKY(2),DCRVSR,
     *                KWORK(L(TRIA(OLVERT))),KWORK(L(TRIA(OLADJ))),
     *                IEL,IELK(2)).NE.0) THEN

C         Element not found. Set it to 0. May happen on boundary nodes 
C         (then we have to set IEL to 0 for RXBNPT), but should not 
C         happen for inner nodes!
C         If this happens for a boundary node, use the boundary 
C         projection routine to recalculate the current element.

          IF (IBND.GT.0) THEN
            IELK(2) = 0
            CALL RXBNPT(TRIA,IBND,IBCT,TMAX(IBCT),IELK(2),DKX(2),DKY(2),
     *                  DPRM,IELK(2))
          ELSE
            IF (MT.GE.2) THEN
              WRITE (MTERM,*) 'MVGPT3-Warning: RK3 substep 2: '//
     *        ' Inner node moved to nowhere!'
              WRITE (MTERM,*) 'IVT=',IVT
              WRITE (MTERM,*) 
     *               ' (X,Y)=(',DCRVSR(1,IVT),',',DCRVSR(2,IVT),')'
              WRITE (MTERM,*) ' -> (',DKY(2),',',DKY(2),')'
            END IF
              
C            MVGPT3 = 1
C            RETURN

C           It's only a warning, we don't halt the program here.
C           This case should normally not occur. But since we computed
C           only an intermediate point, we don't want to treat this as 
C           an error here. Instead, if the point left the domain,
C           we project it back and continue working with the
C           projected point on the boundary. this is not full order
C           anymore, but perhaps it helps to repair a nearly
C           corrupted mesh...

            IELK(2) = 0
            IBC = 0
            CALL RXBNPT(TRIA,IBND,IBC,TMAX(IBC),IELK(2),DKX(2),DKY(2),
     *                  DPRM,IELK(2))
            
          END IF
        END IF

C       Find the element on the level of the solution of the linear
C       system. If the level of the solution is smaller, we can obtain
C       the element number with the help of the 2-level ordering;
C       otherwise we find it by hierarchical search!
        
        IF (ILVSOL.LT.ILEV) THEN
        
          IELSOL = IELK(2)
          DO J=ILEV-1,ILVSOL,-1
            IF (IELSOL.GT.TRIAS(ONEL,J)) THEN
              IELSOL = (IELSOL-1-TRIAS(ONEL,J)) / 3 + 1
            END IF
          END DO
          
C         Unfortunately, this may not be the element that contains the
C         vertex. Imagine:
C
C         |----|\
C         |    | X  <- X is not in the coarse grid element anymore
C         |----|/
C         
C         Therefore, we start with a search on the coarser
C         level to finally find the element that contains X!
C         We use the hierarchical search on one level only, since
C         that will fall back to line search in case that the element
C         could not be found...

          CALL PSRCH5(TRIAS,ILVSOL,ILVSOL,DKX(2),DKY(2),IELSOL,IELSOL)
        
        ELSE IF (ILVSOL.GT.ILEV) THEN
        
          CALL PSRCH5(TRIAS,ILEV,ILVSOL,DKX(2),DKY(2),IELK(2),IELSOL)
        
        ELSE 
        
          IELSOL = IELK(2)
          
        END IF
        
C       Again use the multiple-evaluation routine to calculate the
C       values of DF and DG - this time in the new point. This time
C       save the results in DVAL(7..8).
     
        CALL MSCVQR (2,TNDCNT(TRIAS(1,ILVSOL)),LSOLAR(7),FCLC(7),
     *               E011,.FALSE.,
     *               DKX(2),DKY(2),IELSOL,
     *               TRIAS(1,ILVSOL),DCRVSR,
     *               KWORK(L(TRIAS(OLVERT,ILVSOL))),
     *               KWORK(L(TRIAS(OLMID,ILVSOL))),
     *               DVAL(7), DVGRX(7), DVGRY(7))
     
C       Calculate the denominator of the RHS and with that K3.
C       Remember that here we evaluate at the time step DT_n!

        DF   = DVAL(7)
        DG   = DVAL(8)

        DTMP2 = DT + DDT
        DTMP1 = 1D0 / ( DTMP2/DF + (1D0-DTMP2)/DG )
      
        K(1,3) = DGRX*DTMP1
        K(2,3) = DGRY*DTMP1

C       Before moving the point, save the current position

        DXOLD = DX
        DYOLD = DY
        IELOLD = IEL

C       Finally move the point according to the RK-3 formula:

        DX = DXOLD + (1D0/6D0) * DDT * (K(1,1) + 4D0*K(1,2) + K(1,3))
        DY = DYOLD + (1D0/6D0) * DDT * (K(2,1) + 4D0*K(2,2) + K(2,3))

C       Search in which element the point has been moved to. This must also
C       be done in the last step for boundary nodes to ensure proper 
C       projection with RXBNPT later!
          
        IF (BSRCH4(DX,DY,DCRVSR,
     *             KWORK(L(TRIA(OLVERT))),KWORK(L(TRIA(OLADJ))),
     *             IELOLD,IEL).NE.0) THEN
     
C         Element not found. Set it to 0. May frequently happen on boundary
C         nodes (then we have to set IEL to 0 for RXBNPT), but for inner
C         nodes this is an error!

          IF (IBND.GT.0) THEN

            IEL = 0

C           If we have a boundary node, we have to project the coordinates
C           back onto the boundary. Also determine the new element number of
C           the node in (X,Y), since there is IEL=0 at the moment!

            CALL RXBNPT(TRIA,IBND,IBCT,TMAX(IBCT),IEL,DX,DY,DPRM,IEL)
            
          ELSE
            IF (MT.GE.2) THEN
              WRITE (MTERM,*) 
     *               'MVGPT4-Warning: Inner node moved to nowhere!'
              WRITE (MTERM,*) 'IVT=',IVT
              WRITE (MTERM,*) 
     *               ' (X,Y)=(',DCRVSR(1,IVT),',',DCRVSR(2,IVT),')'
              WRITE (MTERM,*) ' -> (',DX,',',DY,')'
            END IF
            
            MVGPT4 = 1
            RETURN
            
          END IF
        END IF
     
        IF (IBND.GT.0) THEN

C         If we have a boundary node, we have to project the coordinates
C         back onto the boundary. Also determine the new element number if
C         the coordinate felt out of the domain, since there is IEL=0 
C         in that case!

          CALL RXBNPT(TRIA,IBND,IBCT,TMAX(IBCT),IEL,DX,DY,DPRM,IEL)
            
        END IF
     
      END DO
      
C     The final position is in DX/DY. Save it in the grid.

      DCRVDS(1,IVT) = DX
      DCRVDS(2,IVT) = DY
      
C     And return the new parameter value in DPARM...
      
      IF (IBND.GT.0) THEN
        DPARM = DPRM
      END IF
      
      MVGPT4 = 0
      
      END
                        
************************************************************************
* Grid adaption with static method.
*
* This routine performs dynamic grid adaption for one computational
* grid with the static method.
* Every point is moved separately, thus allowing adaptive time-stepping
* techniques. The parameters determine the behaviour of the routine.
*
* In:
*   NLMIN  - Minimum level in TRIAS
*   NLMAX  - Maximum level in TRIAS
*   TRIAS  - array [1..NLMAX] of STRIA
*            = array [1..SZTRIA,1..NLMAX] of integer
*            Triangulation structures for every level determining
*            the grids. This is a grid hierarchy that allowes the 
*            use of a multigrid solver.
*            Level NLMIN..NLMAX must be filled with data.
*   ILEV   - Level inside of the grid hierarchy that should be
*            adapted. NLMIN <= ILEV <= NLMAX
*   ILVSLV - Level inside of the grid hierarchy that should be
*            used for computing the solution of the linear systems.
*            NLMIN <= ILVSLV <= NLMAX
*   IPARAM - array [1..SZGACI+SZGLSI+SZGASI] of integer
*            Integer parameter-array determining the behaviour of the
*            routine.
*            This parameter array must contain the following
*            three blocks, right after another:
*              Index [ 1 .. SZGACI ]
*                    TIGridAdaptionCommon-structure
*              Index [ SZGACI+1 .. SZGACI+SZGLSI ]
*                    TIGridAdaptionSolver-structure
*              Index [ SZGACI+SZGLSI+1 .. SZGACI+SZGLSI+SZGASI ]
*                    TIGridAdaptionStaticParams-structure
*            
*   DPARAM - array [1..*] of double precision
*            Double precision parameter-array determining the
*            behaviour of the routine.
*            This parameter array must contain the following
*            three blocks, right after another:
*              Index [ 1 .. SZGACD ]
*                    TDGridAdaptionCommon-structure
*              Index [ SZGACD+1 .. SZGACD+SZGLSD ]
*                    TDGridAdaptionSolver-structure
*              Index [ SZGACD+SZGLSD+1 .. SZGACD+SZGLSD+SZGASD ]
*                    TDGridAdaptionStaticParams-structure
*
*  IGEOM  - array [1..*] of integer 
*  DGEOM  - array [1..*] of double 
*           Integer- and double-precision parameter blocks with
*           geometry information. Passed to boundary
*           routines. Not used in this routine.
*  IDATA  - array [1..*] of integer 
*  DDATA  - array [1..*] of double 
*           User defined integer- and double-precision parameter 
*           blocks. Passed to callback routines like monitor
*           function. Not used in this routine.
*
* Out:
*   TRIA   - modified triangulation
*   The output-variables in IPARAM/DPARAM are set according to the
*   success of the routine.
*
* For a description about the structure of IPARAM and DPARAM, see
* the documentation in the file "sgridadaptstatic.inc".
************************************************************************

      SUBROUTINE XGASTA (NLMIN,NLMAX,TRIAS,ILEV,ILVSLV,
     *                   IPARAM, DPARAM,IGEOM,DGEOM,IDATA,DDATA)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'stria.inc'
      
      INCLUDE 'sgridadaptgeneralparams.inc'
      INCLUDE 'sgridadaptsolverparams.inc'
      INCLUDE 'sgridadaptstaticparams.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'ssolvers.inc'
      INCLUDE 'm020.inc'
      INCLUDE 'sgridadaptstaticsolver.inc'
      
      INTEGER NLMIN,NLMAX,TRIAS(SZTRIA,*),IGEOM(*),IDATA(*)
      INTEGER IPARAM(SZGACI+SZGLSI+SZGASI)
      DOUBLE PRECISION DPARAM(SZGACD+SZGLSD+SZGASD),DGEOM(*),DDATA(*)
      INTEGER ILEV,ILVSLV
      
C     Structure-array for multigrid solver. These structures are filled
C     with data by the preparation routines. They use these as a basis 
C     of a multigrid solver with different sub-components. All solver-
C     specific data is stored here and released by the destruction 
C     routines when the grid deformation is finished.
C
C     These arrays define instances of the TM020IExtParams and
C     TM020DExtParams structure blocks.
      
      INTEGER MGPARI(SZGSSI)
      DOUBLE PRECISION MGPARD(SZGSSD)
      
C     externals

      INTEGER MVGPT1,MVGPT2,MVGPT3,MVGPT4
      EXTERNAL MVGPT1,MVGPT2,MVGPT3,MVGPT4
      EXTERNAL E011
      INTEGER  NDFGX
      EXTERNAL NDFGX
      EXTERNAL YAXG2,YPG2Q1,YRG2Q1,YSMG2
      EXTERNAL YEXG2,YDBCG2,YSG2Q1,YNMFL2
      
C     local variables

      INTEGER LU, LGRADX, LGRADY, LTMP, LSOLMG, LRHSMG
      INTEGER LMON,LMNOLD,LSIZE,LSIZEM,LCOUNT,LRATIO,LCOLOR
      
      INTEGER CTRIA(SZTRIA),TRIFL
      INTEGER NVT,NEL,NEQ
      INTEGER IGSOUT,NGASTP, ISMTH,NGAPCS,ILVMAX
      INTEGER ITER,I,J,KAUX,KXNPR,KVBDP,IBND, IBCT,IODE,MXODIT
      INTEGER KCORVG
      DOUBLE PRECISION DMSR,DPARM,D
      
C     Block for saving timing-information

      DOUBLE PRECISION TIMIN(SZGACD)
      
C     Clear the timing-statistic parameters we set in this routine:

      DPARAM(OTGTOT)  = 0D0
      DPARAM(OTGGAPR) = 0D0

      DPARAM(OTGLSPR)= 0D0
      DPARAM(OTGLS)  = 0D0
      DPARAM(OTGODE) = 0D0
      DPARAM(OTGGSM) = 0D0
      DPARAM(OTGGRC) = 0D0
      DPARAM(OTGMSM) = 0D0
      
      IPARAM(ONLSITE)= 0
      IPARAM(ONLSYS) = 0
      IPARAM(ONCODE) = 0
      IPARAM(ONMFEVL)= 0
      
C     Save the start time - in our local timing-block!

      CALL GTMAUX (TIMIN,DPARAM,OTGTOT,0)

C     Let's go...

      IF (MT.GE.2) THEN
        WRITE (MTERM,'(A)') '======================================='
        WRITE (MTERM,'(A)') 'Starting grid deformation'
      END IF

      NVT = TRIAS(ONVT,ILEV)
      NEL = TRIAS(ONEL,ILEV)
      
      IF (NVT.LE.0) THEN
        WRITE (MTERM,'(A)') 'XGASTA: Error: NVT <= 0!'
        RETURN
      END IF
      IF (NEL.LE.0) THEN
        WRITE (MTERM,'(A)') 'XGASTA: Error: NEL <= 0!'
        RETURN
      END IF
      
      IGSOUT = IPARAM (OIGSOUT)
      NGASTP = IPARAM (ONGASTP)
      NGAPCS = IPARAM (SZGACI+SZGLSI+ONGAPCS)
      
C     Convert the TRIA-structure on the level that we should deform
C     into an extended triangulation structure if necessary; 
C     allowes quicker computation due to extended information.

      TRIFL = TRIAS(OTRIFL,ILEV)
      IF (IAND(TRIFL,1).EQ.0) CALL GENETR (TRIAS(1,ILEV))
      
C     General approach:
C
C     CTRIA saves the current deformed triangulation, while
C     TRIAS(1,ILEV) (and the fine-grid midpoints on the higher levels)
C     are only changed after a deformation.
C     A number of macrosteps is performed to consecutively enhance
C     the grid. In every macrostep one Laplacian problem is solved
C     to generate a vector field. This vector field is used to
C     move all points in the grid to a new position.
C
C     At first we copy the current grid to CTRIA:

      CALL LCP3(TRIAS(1,ILEV),CTRIA,SZTRIA)

C     This only copied the handles, not the data. But during the grid
C     deformation we need the original as well as the deformed grid.
C     Starting from the original grid DCORVG and DVBDP in CTRIA will 
C     be modified. We therefore copy these both arrays so that they can
C     be modified without danger:

      CALL ZNEW (2*NVT,-1,CTRIA(OLCORVG),'DCORV2')
      CALL ZNEW (TRIAS(ONVBD,ILEV),-1,CTRIA(OLVBDP),'LVBDP2')
      CALL LCP1 (DWORK(L( TRIAS(OLCORVG,ILEV) )),
     *           DWORK(L( CTRIA(OLCORVG) )),2*NVT)
      CALL LCP1 (DWORK(L( TRIAS(OLVBDP,ILEV) )),
     *           DWORK(L( CTRIA(OLVBDP) )),CTRIA(ONVBD))
      
C     Reserve some memory.
C
C     Calculate the maximum level that arises in our grid hierarchy.
C     We need this for reserving memory:
C     Solution vectors, gradient vectors,... are computed only on
C     the "solution" level ILVSLV. Actually the monitor function also,
C     but not always! For GMV output, we need to calculate some stuff
C     like monitor function etc. also on level ILEV!
C     Therefore, some vectors have to be allocated on the maximum
C     occurring level, so they are capable of holding information
C     of both levels!

      ILVMAX = MAX(ILEV,ILVSLV)

C     At first we need memory for a solution vector,
C     and two vectors representing the derivative of the solution.
C     Therefore, call NDFGL to determine the number of degrees of
C     freedom for our element E011 on the level where we compute
C     the solution of the corresponding LGS:

      I = -1
      CALL E011 (0D0,0D0,I)
      NEQ = NDFGX(I,TRIAS(1,ILVSLV))

C     Allocate the memory...

      CALL ZNEW(NEQ,1,LU,'DU    ')
      CALL ZNEW(NEQ,1,LGRADX,'GRADX ')
      CALL ZNEW(NEQ,1,LGRADY,'GRADY ')

      CALL ZNEW(NEQ,1,LTMP,'DTMP  ')
      
C     Allocate a handle for the current monitor function LMON,
C     as well as a handle for the cell-size distribution in the
C     adaption itself.
C     For output purposes, LMON must be large enough to hold 
C     the monitor function on the maximum occurring level.
C     LMNOLD must only hold the monitor function on the level
C     where the LGS is solved.
      
      CALL ZNEW(TRIAS(ONVT,ILVMAX),1,LMON  ,'DMON  ')
      CALL ZNEW(TRIAS(ONVT,ILVSLV),1,LMNOLD,'DMNOLD')
      
C     Average size of the elements around each vertex:

      CALL ZNEW(TRIAS(ONVT,ILVMAX),1,LSIZE ,'SIZE  ')
      
C     Number of elements meeting in a vertex:

      CALL ZNEW(TRIAS(ONVT,ILVMAX),3,LCOUNT,'COUNT ')
      
C     Area of each element:

      CALL ZNEW(TRIAS(ONEL,ILVMAX),1,LSIZEM,'SIZEM ')
      
C     Ratio of area of each element in comparison to all 
C     neighbour elements:
      
      CALL ZNEW(TRIAS(ONEL,ILVMAX),1,LRATIO,'RATIO ')
      
C     Colorize the macros using the coarse grid - for GMV-output.
C     As the color array is used for GMV output on the level of
C     the solver as well as on the level of the main grid, we
C     have to calculate it for the maximum occurring level!

      CALL ZNEW(TRIAS(ONEL,ILVMAX),3,LCOLOR,'COLOR ')
      
      IF (IPARAM(OICOLMC).NE.0) THEN
        CALL COLMCR(TRIAS(ONEL,NLMIN),KWORK(L(TRIAS(OLVERT,NLMIN))),
     *            DWORK(L(TRIAS(OLCORVG,NLMIN))),TRIAS(ONEL,ILVMAX),
     *            KWORK(L(TRIAS(OLVERT,ILVMAX))),KWORK(L(LCOLOR)))
      ELSE
        DO I=0,TRIAS(ONEL,ILVMAX)-1
          KWORK(L(LCOLOR)+I) = 1
        END DO
      END IF
      
C     Write out the initial grid

      IF ( IAND(IGSOUT,65).NE.0 ) THEN
      
C       Initialise element ststistics for writing the GMV file:

        CALL MEASR2(DWORK(L(TRIAS(OLCORVG,ILEV))),
     *              KWORK(L(TRIAS(OLVERT,ILEV))),
     *              KWORK(L(TRIAS(OLADJ,ILEV))),NEL,
     *              NVT,KWORK(L(LCOUNT)),DWORK(L(LSIZE)),
     *              DWORK(L(LSIZEM)),DWORK(L(LRATIO)))
     
C       Furthermore determine the monitor function f(x) on the
C       calculation level. Don't weight this by AREA, i.e.
C       call the routine with IMNREL=1!
C       This is only a temporary monitor function for the output
C       of the initial grid. The actual monitor function, which
C       is used for the deformation, is recreated later on.

        CALL XMON0(TRIAS(1,ILVSLV),LSIZE,LMON,1,DPARAM(OEPS0GS),
     *             IGEOM,DGEOM,IDATA,DDATA)

C       Write the GMV file

        CALL GAGMVS (CTRIA,-1,ILEV,0,0,
     *               TRIAS(ONEL,NLMIN),LCOLOR,LSIZE,LSIZEM,LRATIO,
     *               LMON,0,0,0,0,
     *               'gmv/INITIAL')    
     
      END IF

C     In each macro-step we have to solve a linear system. We can already
C     initialize our solver for that system in advance. The definition
C     how the solver should be initialized can be found in the
C     TxGridAdaptionSolver substructure in IPARAM/DPARAM, which starts
C     at index SZGACx+1 of IPARAM/DPARAM.
C     This also returns handles for solution and RHS-vector, so
C     we know where to write the RHS to and where to read the
C     calculated solution from. Furthermore it returns the number of
C     equations NEQ - which is the same as computed earlier.

      CALL GTMAUX (TIMIN,DPARAM,OTGLSPR,0)
      CALL INMGG2 (IPARAM(SZGACI+1), DPARAM(SZGACD+1), ILVSLV, TRIAS, 
     *             MGPARI, MGPARD, LSOLMG, LRHSMG, NEQ)
      CALL GTMAUX (TIMIN,DPARAM,OTGLSPR,1)
     
C     Start with the macro-steps

      DO ITER=1,NGASTP+NGAPCS
      
        IF (MT.GE.2) THEN
          WRITE (MTERM,'(A)') '======================================='
          IF ( ITER.LE.NGASTP ) THEN
            WRITE (MTERM,'(A,I5)') 'Macrostep ',ITER
          ELSE
            WRITE (MTERM,'(A,I5,A)') 'Macrostep ',ITER,
     *            ' (Postcorrection)'
          END IF
          WRITE (MTERM,'(A)') '======================================='
        END IF
      
        IF ((IPARAM(ONPRGSM).GT.0).AND.(MT.GE.2)) THEN
          WRITE (MTERM,'(A)') 'Grid-Presmoothing'
        END IF 
        
C       Before we really perform the grid deformation, we will first
C       smooth the grid a little bit using Laplacian grid smoothing.
C       This can get rid of numerical instzalilities and cruel grid
C       configurations, which might destroy the grid when performing
C       multiple macro-steps.

C       Perform NPRGSM smoothing steps, as indicated in the parameter
C       block. The smoothing is performed on the original grid. The
C       resulting grid serves as a new source-grid for the solution of
C       the linear system below as well as for the grid deformation.
C       Performing the smoothing only on CTRIA is not enough, as we need
C       a better "starting" mesh also for our linear system!

        CALL GTMAUX (TIMIN,DPARAM,OTGGSM,0)

        DO ISMTH=1,IPARAM(ONPRGSM)
        
          IF ((MT.GE.2).AND.(MOD(ISMTH,20).EQ.1)) THEN
            WRITE (MTERM,'(A,I5)') 'Grid-smoothing iteration ',ISMTH
          END IF
          
C         Smooth the grid. This affects the coordinates of all points 
C         as well as the parameter values of the boundary points.

          IF (IPARAM(OISMMTH).EQ.1) THEN 
            CALL GSMTH6(CTRIA)
          ELSE
            CALL GSMTH5(CTRIA)
          END IF

C         Write the output of the laplacian grid smoothing iterations
C         to files - every iterations / 10 iterations / 100 
C         iterations...

          IF ( (IAND(IGSOUT,4).NE.0) .OR.
     *     ( (IAND(IGSOUT,8).NE.0)  .AND. (MOD(ISMTH,10).EQ.0)  ) .OR.
     *     ( (IAND(IGSOUT,16).NE.0) .AND. (MOD(ISMTH,100).EQ.0) ) .OR.
     *     ( (IAND(IGSOUT,32).NE.0) .AND. 
     *       (ISMTH.EQ.IPARAM(ONPRGSM)) )) THEN

            CALL GTMAUX (TIMIN,DPARAM,OTGGSM,1)
     
C           Recalculate element statistics for output.
C           If no GMV output is to be written, the element statistics
C           are only calculated after the smoothing, which is *much*
C           faster :-) ...
     
            CALL MEASR2(DWORK(L(TRIAS(OLCORVG,ILEV))),
     *              KWORK(L(TRIAS(OLVERT,ILEV))),
     *              KWORK(L(TRIAS(OLADJ,ILEV))),NEL,
     *              NVT,KWORK(L(LCOUNT)),DWORK(L(LSIZE)),
     *              DWORK(L(LSIZEM)),DWORK(L(LRATIO)))

            CALL GAGMVS (CTRIA,-1,ILEV,ITER,ISMTH,
     *                   TRIAS(ONEL,NLMIN),LCOLOR,LSIZE,LSIZEM,LRATIO,
     *                   0,0,0,0,0,
     *                   'gmv/SMPRE')    

            CALL GTMAUX (TIMIN,DPARAM,OTGGSM,0)
     
          END IF

        END DO

        CALL GTMAUX (TIMIN,DPARAM,OTGGSM,1)

        IF ((MT.GE.2).AND.(IPARAM(ONPRGSM).NE.0).AND.
     *      (MOD(IPARAM(ONPRGSM),20).NE.1)) THEN
          WRITE (MTERM,'(A,I5)') 'Grid-smoothing iteration ',ISMTH
        END IF
        
C       Finally use the resulting grid as new base grid for the linear
C       system, e.g. for the whole grid adaption process. Copy CTRIA
C       to TRIA so they are the same again.
        
        IF (IPARAM(ONPRGSM).GT.0) THEN
        
          CALL LCP1 (DWORK(L( CTRIA(OLCORVG) )),
     *               DWORK(L( TRIAS(OLCORVG,ILEV) )),2*NVT)
          CALL LCP1 (DWORK(L( CTRIA(OLVBDP) )),
     *               DWORK(L( TRIAS(OLVBDP,ILEV) )),TRIAS(ONVBD,ILEV))
     
C         Correct midpoint information on all finer levels
C         if the current level is less than the solver level.
C         If the current level is larger than the solver level,
C         the sharing of the coordinates will automatically ensure
C         that the coarse grid coordinates are corrected as well.

          DO I = ILEV,ILVSLV-1
            CALL PRC2FG (DWORK(L( TRIAS(OLCORVG,I) )),
     *                   KWORK(L( TRIAS(OLVERT,I) )),
     *                   TRIAS(ONEL,I),
     *                   DWORK(L( TRIAS(OLCORVG,I+1) )),
     *                   KWORK(L( TRIAS(OLVERT,I+1) )),
     *                   KWORK(L( TRIAS(OLADJ,I+1) )),
     *                   TRIAS(ONEL,I+1))
          END DO
     
          IF (MT.GE.2) THEN
            WRITE (MTERM,9000) 
          END IF
        END IF
      
C       The grid deformation in each macro-step is now performed
C       in three steps. For out domain Omega we want to find a
C       bijective mapping Phi:Omega->Omega with:
C
C          g(x) |J Phi(x)| = f(Phi(x))    for all x
C
C       Here g is the initial cell size distribution and f is a given
C       monitor function. If this mapping is found, we can move the
C       points of our grid by y=Phi(x) and obtain a grid with a cell
C       distribution determined by f(.): With g(.) being the initial
C       cell size and |J Phi(x)| being the relative growth of a cell
C       around x, g(x) |J Phi(x)| is the absolute size of the cell
C       after it has been moved. f(Phi(x)) is now a given monitor
C       function at the point Phi(x) where a node x will move to.
C       So if Phi(x) is found fulfilling the above equation,
C       it is designed that way that determines (up to a constant
C       factor) the absolute cell size at the point Phi(x) where a 
C       given node x is moved to. 
C
C       This mapping Phi(.) can't be determined explicitly, it's
C       only implicitly given by a couple of formulas. We will not
C       determine the actual mapping Phi(.), but instead we will
C       directly determine the position y=Phi(x) of all nodes Phi(.).
C       The approach is the following:
C       1.) Determine a scalar function v(x) being a potential
C           of a vector field of "grid velocity vectors"
C       2.) Determine the vector field of the "grid velocity":
C           w(x) := grad v(x)
C       3.) Move all grid nodes x along this velocity field w(x)
C           using an ODE-solver
C
C       At first calculate the current distribution of the cells 
C       to determine the function g(x); the handle LSIZE points to
C       this function:

        CALL GTMAUX (TIMIN,DPARAM,OTGGAPR,0)

        CALL MEASR2(DWORK(L(TRIAS(OLCORVG,ILVSLV))),
     *              KWORK(L(TRIAS(OLVERT,ILVSLV))),
     *              KWORK(L(TRIAS(OLADJ,ILVSLV))),TRIAS(ONEL,ILVSLV),
     *              TRIAS(ONVT,ILVSLV),KWORK(L(LCOUNT)),DWORK(L(LSIZE)),
     *              DWORK(L(LSIZEM)),DWORK(L(LRATIO)))

C       Furthermore determine the monitor function f(x) on the
C       calculation level. Don't weight this by AREA, i.e.
C       call the routine with IMNREL=1!

        CALL XMON0(TRIAS(1,ILVSLV),LSIZE,LMON,1,DPARAM(OEPS0GS),
     *             IGEOM,DGEOM,IDATA,DDATA)
        
C       Smooth that monitor function a little bit if desired to
C       avoid jumps

        IF ((IPARAM(ONFSMTH).GT.0).AND.(MT.GE.2)) THEN
          WRITE (MTERM,'(A)') 'Monitor function smoothing'
        END IF 

        CALL GTMAUX (TIMIN,DPARAM,OTGMSM,0)

        DO ISMTH=1,IPARAM(ONFSMTH)
        
          IF ((MT.GE.2).AND.(MOD(ISMTH,20).EQ.1)) THEN
            WRITE (MTERM,'(A,I5)') 'Monitor function smoothing '//
     *                             'iteration ',ISMTH
          END IF
          
C         Smooth the grid. This affects the coordinates of all points 
C         as well as the parameter values of the boundary points.

          CALL GSMTH9(TRIAS(1,ILVSLV),TRIAS(ONVT,ILVSLV),DWORK(L(LMON)))

        END DO

        CALL GTMAUX (TIMIN,DPARAM,OTGMSM,0)

        IF (MT.GE.2) THEN
          IF ((IPARAM(ONFSMTH).NE.0).AND.
     *      (MOD(IPARAM(ONFSMTH),20).NE.1)) THEN
            WRITE (MTERM,'(A,I5)') 'Monitor function smoothing '//
     *                             'iteration ',ISMTH
          END IF
          WRITE (MTERM,9000) 
        END IF
        
C       Furthermore we need the cell distribution as comparison
C       function g(x). We copy DSIZE to DMNOLD and use DMNOLD as g(x):

        CALL XLCP1(LSIZE,LMNOLD,NEQ)

C       Implementation of a relative monitor function. We scale the
C       monitor function f(.) by g(.), which is the original area.
C       As the previous f(.) represents the absolute cell size 
C       distribution, scaling f(.) by area(.) will result in a
C       relative cell size distribution.

        IF (IPARAM(OIMNREL).EQ.1) THEN
          DO I=0,NEQ-1
            DWORK(L(LMON)+I) = DWORK(L(LMON)+I)*DWORK(L(LMNOLD)+I)
          END DO
        END IF
        
C       Calculate the area of the domain

        CALL MEASURE(TRIAS(ONEL,ILVSLV),DWORK(L(LSIZEM)),DMSR)
        
C       Calculate the quality measure. This might differ from
C       the calculated quality measure at the end of the previous
C       macrostep because of the grid smoothing!

        CALL GAQMSR (TRIAS(ONVT,ILVSLV),DWORK(L(LMON)),
     *               DWORK(L(LSIZE)),DMSR,D)
        
        IF (MT.GE.2) THEN
          WRITE (MTERM,'(A,F10.5)') 'Initial quality measure: ',D
        END IF

C       In the first NGASTP steps, we perform some regularization to f(x) 
C       to force slower and saver convergence. In this case we make a 
C       blending:
C                        f(x) := s*f(x) + (1-s)*g(x)
C       with s a parameter in [0,1]. We choose this parameter according
C       to the current macrostep number.
C       The second NGAPCS mactosteps are then performed without blending,
C       so performing a fine post-correction to the grid.

        IF ( ITER.LT.NGASTP ) THEN
        
C         We need some normalization of f(x) and g(x) in order for the
C         blending to work properly. In detail we now scale f(x) and
C         g(x) such that 
C               Int(f) = Int(g) = |Omega|
C         The new f(x) will then also fulfill the equation
C         Int(f)=|Omega|.

          IF (IPARAM(OICUNGS).EQ.0) THEN
            CALL NORMIF(DWORK(L(LSIZE)),NEQ,
     *                  DMSR,DWORK(L(LMNOLD)))
            CALL NORMIF(DWORK(L(LSIZE)),NEQ,
     *                  DMSR,DWORK(L(LMON)))
          ELSE
            CALL NRMM2(TRIAS(1,ILVSLV),DWORK(L(LMNOLD)),
     *                 DMSR,NEQ,E011,IPARAM(OICUNGS))
            CALL NRMM2(TRIAS(1,ILVSLV),DWORK(L(LMON)),
     *                 DMSR,NEQ,E011,IPARAM(OICUNGS))
          END IF
        
C         Calculate blending parameter, final grid should be archived
C         in step NGASTP. D is a value between 0..1. 

          IF ((DPARAM(ODMNSFC).LT.0D0).OR.(DPARAM(ODMNSFC).GT.1D0)) THEN

C           Formula 1: Unique distribution.
C           Take 2x a square to regularize the adaption speed (because the 
C           adaption is the faster, the more D is around 1!).
        
            D = SQRT(SQRT (DBLE(ITER)/DBLE(NGASTP)))

          ELSE

C           Formula 2: Take 1-DMNSFC, 1-DMNSFC^2, 1-DMNSFC^3,...
C           to respect the fast adaption in the last steps
      
            IF (ITER.EQ.NGASTP) THEN
              D = 1D0
            ELSE
              D = 1D0 - DPARAM(ODMNSFC)**ITER
            END IF
            
          END IF
          
C         Calculate the blended monitor function

          DO I=0,NEQ-1
            DWORK(L(LMON)+I) =       D*DWORK(L(LMON)+I) + 
     *                         (1D0-D)*DWORK(L(LMNOLD)+I)
          END DO
        
        END IF
        
C       Scale both the monitor function and the cell-distribution
C       function to have
C               integral (1/f) = integral (1/g) = |Omega|
C       The last normalization "= |Omega|" is normally not necessary 
C       for the method, integral (1/f) = integral (1/g) would be
C       enough. Nevertheless normalizing both functions is easier
C       to handle.
C
C       Use quick-and-dirty weighting with mean area around vertices
C       (corresponds to a 1x1-Gauss-Formula), or use real
C       quadrature.

        IF (IPARAM(OICUNGS).EQ.0) THEN
          CALL XNORM(LSIZE,NEQ,DMSR,LMNOLD)
          CALL XNORM(LSIZE,NEQ,DMSR,LMON)
        ELSE
          CALL NRMM2R(TRIAS(1,ILVSLV),DWORK(L(LMNOLD)),
     *                DMSR,NEQ,E011,IPARAM(OICUNGS))
          CALL NRMM2R(TRIAS(1,ILVSLV),DWORK(L(LMON)),  
     *                DMSR,NEQ,E011,IPARAM(OICUNGS))
        END IF
      
C       Now we have to solve the system:
C
C              div v(x) = 1/f - 1/g
C
C       This gives our potential v(x) with the derivative being the
C       grid-velocity vector field.
C
C       The first step in setting up this system is to build the
C       RHS vector.
     
        CALL GTMAUX (TIMIN,DPARAM,OTGGAPR,1)
      
C       If the cubature formula for building the RHS has a number > 0, 
C       call PRRLX2 to calculate the "exact" RHS. If the cubature 
C       formula is < 0, call PRRLG2 to calculate an approximate RHS, 
C       which is a little bit faster.

        IF (MT.GE.2) THEN
          WRITE (MTERM,9000) 
          WRITE (MTERM,'(A)') 'Preparing linear system'
        END IF 
        
        CALL GTMAUX (TIMIN,DPARAM,OTGLSPR,0)

        IF (IPARAM(SZGACI+OICURGS).LT.0) THEN

C         We set up - pointwise - the right hand side 1/f-1/g and
C         save this (finite difference) vector in the right-hand-side
C         vector of the solver structure. The resulting vector will
C         be saved au DTMP:

          KAUX = L( LTMP )
          DO I=0,NEQ-1
            DWORK(KAUX+I) = 1D0/DWORK(L(LMON)+I) - 
     *                      1D0/DWORK(L(LMNOLD)+I)
          END DO
      
C         Now interpret this as a piecewise linear FE-function in the
C         nodes of the grid, which values are determined by the
C         piecewise linear element E011. Call the routine PRRLG2 to
C         create a FE-RHS from this function. Save this RHS in the
C         DU-vector temporarily.
      
          CALL PRRLG2 (IPARAM(SZGACI+1), DPARAM(SZGACD+1), 
     *                 TRIAS(1,ILVSLV),LTMP,LU,NEQ)
        ELSE
     
C         We use PRRLX2 to calculate the "exact" RHS-vector by 
C         integration. PRRLX2 accepts both DMON and DMNOLD and sets 
C         up the FE-RHS-vector, which was defined as (1/DMON+1/DMNOLD,phi).
C         Save the RHS vector in the DU-vector temporarily.
        
          CALL PRRLX2 (IPARAM(SZGACI+1), DPARAM(SZGACD+1),
     *                 TRIAS(1,ILVSLV),
     *                 DWORK(L(LMON)),DWORK(L(LMNOLD)),DWORK(L(LU)),
     *                 NEQ,E011)
     
        END IF
        
C       Call the preparation routine of the solver to calculate
C       system matrix according to the current triangulation
C       structure TRIA of the underlying grid. It will also
C       set up all index/coefficient vectors and the smoothing
C       count arrays:

        CALL PRMGG2 (IPARAM(SZGACI+1), DPARAM(SZGACD+1), 
     *               TRIAS, MGPARI, MGPARD)
     
C       Don't forget to set up the right hand side - that is
C       still temporarily saved in LU!
C       Also set the starting vector to 0...

        CALL LCP1 (DWORK(L(LU)),DWORK(L(LRHSMG)),NEQ)
        CALL LCL1 (DWORK(L(LSOLMG)),NEQ)

C       Normalize the right hand side to have integral-mean-value=0.
C       Necessary because the solution will also posess this
C       property; otherwise the residuum will never get 0, it
C       will converge to a constant <> 0!
C       Normally this is not necessary - the RHS should have integral
C       mean value=0. But we do this for sure, because integral
C       formulas of lower order (like 1x1 Gauss) might introduce
C       an error large enough to disturb this property!
C
C       Remark: This trick to make 1x1-Gauss work properly works well
C       in FEAT - but does not work with FEAST because the RHS-vector
C       in FEAST can be split between multiple macros, and therefore
C       integral-mean-value=0 on all macros does not mean (because of
C       duplicate values on the boundary) integral-mean-value=0 totally!
C       In FEAST this command can simply be dropped, as long as a 
C       higher-order integration formula is used.

        CALL IMVZER (DWORK(L(LRHSMG)),NEQ)
        
        CALL GTMAUX (TIMIN,DPARAM,OTGLSPR,1)
        
C       So far so good, the system has been set up. Now it's time
C       to solve that system to calculate our potential:

        IF (MT.GE.2) THEN
          WRITE (MTERM,9000) 
          WRITE (MTERM,'(A)') 'Solving linear system'
        END IF 

C       Set the output level of the solver to the general 
C       output level - 2

        MGPARI(OMSGTRM) = MT-2

C       Set the output level of the coarse grid solver to the 
C       general output level - 3

        MGPARI(OMCGTRM) = MT-3

C       Set the output level of the smoother to the 
C       general output level - 4

        MGPARI(OMCGTRM) = MT-4

        CALL GTMAUX (TIMIN,DPARAM,OTGLS,0)

C       IDATA/DDATA in the call to M020 is not used and 
C       simply set to MGPARI,MGPARD.

        CALL M020 (MGPARI,MGPARD, MGPARI,MGPARD, 
     *             DWORK(1),DWORK(1),DWORK(1),
     *             YAXG2,YPG2Q1,YRG2Q1,YSMG2,YSMG2,
     *             YEXG2,YEXG2,YDBCG2,YSG2Q1,YNMFL2)

        CALL GTMAUX (TIMIN,DPARAM,OTGLS,1)
        
C       Collect statistical information from MG-solver:

        IPARAM(ONLSITE) = IPARAM(ONLSITE) + MGPARI(OITE)
        IPARAM(ONLSYS)  = IPARAM(ONLSYS) + 1
     
C       Copy the solution to LU, away from the multigrid solver.
C       Let's hope that NEQ=NVT for now ;)

        CALL LCP1 (DWORK(L(LSOLMG)),DWORK(L(LU)),NEQ)
        
C       Now we have it in LU. In the next step build the
C       reconstructed gradient, which is the velocity vector
C       field for our grid movement:

        CALL GTMAUX (TIMIN,DPARAM,OTGGRC,0)

        CALL XLCL1 (LGRADX,NEQ)
        CALL XLCL1 (LGRADY,NEQ)
        
C       Use ZZ-technique or standard technique for reconstructing 
C       gradients:

        IF (IPARAM(SZGACI+OIRCGRD).EQ.0) THEN
        
C         standard interpolation

          CALL XGRDR3(TRIAS(1,ILVSLV),DWORK(L(LU)),DWORK(L(LGRADX)),
     *                DWORK(L(LGRADY)),NEQ,E011,.FALSE.)

        ELSE

C         ZZ-technique

          CALL XGRDR2(TRIAS(1,ILVSLV),DWORK(L(LU)),DWORK(L(LGRADX)),
     *                DWORK(L(LGRADY)),NEQ,E011)
     
        END IF

        CALL GTMAUX (TIMIN,DPARAM,OTGGRC,1)
        
C       Write the solution into a GMV-file

        IF ((IAND(IGSOUT,1).NE.0).OR.(IAND(IGSOUT,2).NE.0)) THEN
        
C         We still use ILEV and not ILVSLV as "identifier" for the
C         level in the following call to name all files of the whole
C         method in the same way!
C         The parameter is only used for the filename, not for the
C         content of the file.

          CALL GAGMVS (TRIAS(1,ILVSLV),-1,ILEV,ITER,0,
     *                 TRIAS(ONEL,NLMIN),LCOLOR,LSIZE,LSIZEM,LRATIO,
     *                 LMON,LRHSMG,LU,
     *                 LGRADX,LGRADY,'gmv/SOL')    
     
        END IF
        
C       Start the actual moving of the grid - vertex by vertex...

        IF (MT.GE.2) THEN
          WRITE (MTERM,9000) 
          WRITE (MTERM,'(A)') 'Moving grid'
        END IF 
        
C       The velocity vector field can be found in LGRADX/LGRADY;
C       Now it's time to move our grid! We have to move
C       each node separately. This means, for every node we have
C       to solve an ODE how the gridpoint is moved. 

        CALL GTMAUX (TIMIN,DPARAM,OTGODE,0)

C       Let's assume we have a gridpoint x. The goal is now to solve
C       the ODE:
C
C           d/dt phi(x,t) = eta ( phi(x,t), t),     t=0..1
C
C       for the right hand side eta:(Omega x [0,1])->Omega,
C
C           eta (y,s) = - v(y) / ( s/f(y) + (1-s)/g(y) )
C
C       and with the initial condition:  phi (x,0) = x
C
C       We don't want to find phi(.,.) explicitly, as this pseudo-time 
C       dependent function is only an auxiliary construction. What we
C       are interested in is the point:
C
C           Phi (x) = phi (x,1)
C
C       which is the destinatioun point, where the source point x is 
C       moved to. This point can now be calculated via different ODE 
C       sovers, e.g. explicit Euler, Runge Kutta or others.

        KXNPR = L(CTRIA(OLXNPR))
        KVBDP = L(CTRIA(OLVBDP))
        KCORVG = L(CTRIA(OLCORVG))
        IODE = IPARAM(SZGACI+SZGLSI+OODESLV)
        MXODIT = IPARAM(SZGACI+SZGLSI+ONITDT)

C       Make a loop over all vertices on our grid that we want 
C       to deform:

        DO I=1,NVT
 
C         Now we have node I. First decide whether this node is an
C         inner node or a boundary node. This determines how much
C         information we have to collect for the call to our ODE-
C         solver subroutine.
C         We use the extended nodal property array to decide whether
C         we have a boundary point or not:

          IF ( IAND(KWORK(KXNPR+2*(I-1)),2**8).NE.0 ) THEN
C           Boundary node. Search for the boundary node index:
            CALL SRBDNI (CTRIA, I, IBND, IBCT)
          ELSE
C           Inner node. There's not much to collect...
            IBND = 0
            IBCT = 0
          END IF
          
C         If we have a boundary node with an integer-valued parameter,
C         don't move that node. These are our "guiding" nodes, the
C         corners of the geometry!

          IF ( (IAND(KWORK(KXNPR+2*(I-1)),2**8).EQ.0) .OR.
     *        (DWORK(KVBDP+IBND-1).NE.AINT(DWORK(KVBDP+IBND-1)) ) ) THEN
     
C           Depending on the parameter we now choose an ODE-solver
C           to find the "target" coordinates of our point.

            J = 0
            IF (IODE.EQ.-1) THEN
              
C             UNDOCUMENTED PARAMETER:
C             "Old" Euler method, which works only on
C             one level!
              
              J = MVGPT1 (TRIAS(1,ILEV), 
     *                  DWORK(L( TRIAS(OLCORVG,ILEV) )), 
     *                  DWORK(L( CTRIA(OLCORVG) )), 
     *                  I, IBND, IBCT, DPARM,
     *                  LU,LGRADX,LGRADY,
     *                  LMON,LMNOLD,
     *                  MXODIT)
     
C             Increase NMFEVL by MXODIT since the ode needed
C             MXODIT evaluations to be solved
     
              IPARAM(ONMFEVL) = IPARAM(ONMFEVL) + MXODIT
              
            ELSE IF (IODE.EQ.1) THEN
            
C             Explicit Euler, multiple-grid version
            
              J = MVGPT2 (TRIAS(1,ILEV), 
     *                    DWORK(L( TRIAS(OLCORVG,ILEV) )), 
     *                    DWORK(L( CTRIA(OLCORVG) )), 
     *                    I, IBND, IBCT, DPARM,
     *                    TRIAS, ILEV,ILVSLV,
     *                    LU,LGRADX,LGRADY,
     *                    LMON,LMNOLD,
     *                    MXODIT)

C             Increase NMFEVL by MXODIT since the ode needed
C             MXODIT evaluations to be solved

              IPARAM(ONMFEVL) = IPARAM(ONMFEVL) + MXODIT

            ELSE IF (IODE.EQ.-3) THEN
            
C             UNDOCUMENTED PARAMETER:
C             "Old" Runge-Kutta-3 method, which works only on
C             one level!
            
              J = MVGPT3 (TRIAS(1,ILEV), 
     *                  DWORK(L( TRIAS(OLCORVG,ILEV) )), 
     *                  DWORK(L( CTRIA(OLCORVG) )), 
     *                  I, IBND, IBCT, DPARM,
     *                  LU,LGRADX,LGRADY,
     *                  LMON,LMNOLD,
     *                  MXODIT)
     
C             Increase NMFEVL by 3xMXODIT since the ode needed
C             3xMXODIT evaluations to be solved, 3 per step
     
              IPARAM(ONMFEVL) = IPARAM(ONMFEVL) + 3*MXODIT

            ELSE IF (IODE.EQ.3) THEN
            
              J = MVGPT4 (TRIAS(1,ILEV), 
     *                  DWORK(L( TRIAS(OLCORVG,ILEV) )), 
     *                  DWORK(L( CTRIA(OLCORVG) )), 
     *                  I, IBND, IBCT, DPARM,
     *                  TRIAS, ILEV,ILVSLV,
     *                  LU,LGRADX,LGRADY,
     *                  LMON,LMNOLD,
     *                  MXODIT)
     
C             Increase NMFEVL by 3xMXODIT since the ode needed
C             3xMXODIT evaluations to be solved, 3 per step
     
              IPARAM(ONMFEVL) = IPARAM(ONMFEVL) + 3*MXODIT

            END IF
          
C           The routine will result a success-flag in J and
C           (in case of a boundary node) the parameter value in DPARM.
C           It will write the target coordinates of the point into the
C           DCORVG-array of CTRIA.

            IF (J.NE.0) THEN
              WRITE (MTERM,'(A,I5)') 'Fatal: Error moving node ',I
              WRITE (MTERM,'(A,I5)') 
     *        'Grid adaption aborted in macrostep ',ITER
              WRITE (MTERM,*) 
              WRITE (MTERM,'(A)') 
     *        'Try to increase NITDT (=better solution of ODE''s).'
              WRITE (MTERM,'(A)') 
     *        'Try to increase NGASTP (=performs more macrosteps) and'
              WRITE (MTERM,'(A)') 
     *        'put DMNSFC closer to 1 (=decrease mnesh adaption speed)'
              WRITE (MTERM,'(A)') 
     *        'Perhaps perform some presmoothing with NPRGSM>0.'
              GOTO 99990
            END IF

C           In case of a parameter node we have to save the new
C           parameter value. As KVBDP points to the parameter values
C           of the new grid, the original KVBDP-array of TRIA remains
C           unchanged for now.

            IF ( IAND(KWORK(KXNPR+2*(I-1)),2**8).NE.0 ) THEN
              DWORK(KVBDP+IBND-1) = DPARM
            END IF
            
          END IF

        END DO
        
C       Increase NCODE by NVT since we solved NVT ODE's:
        
        IPARAM(ONCODE)  = IPARAM(ONCODE) + NVT
        
        CALL GTMAUX (TIMIN,DPARAM,OTGODE,1)
        
C       Finally perform NPOGSM postsmoothing steps, as indicated in the 
C       parameter block. The smoothing is performed on the original grid. 

        IF ((IPARAM(ONPOGSM).GT.0).AND.(MT.GE.2)) THEN
          WRITE (MTERM,9000) 
          WRITE (MTERM,'(A)') 'Grid-Postsmoothing'
        END IF 
        
        CALL GTMAUX (TIMIN,DPARAM,OTGGSM,0)

        DO ISMTH=1,IPARAM(ONPOGSM)
        
          IF ((MT.GE.2).AND.(MOD(ISMTH,20).EQ.1)) THEN
            WRITE (MTERM,'(A,I5)') 'Grid-smoothing iteration ',ISMTH
          END IF
          
C         Smooth the grid. This affects the coordinates of all points 
C         as well as the parameter values of the boundary points.

          IF (IPARAM(OISMMTH).EQ.1) THEN 
            CALL GSMTH6(CTRIA)
          ELSE
            CALL GSMTH5(CTRIA)
          END IF

C         Write the output of the laplacian grid smoothing iterations
C         to files - every iterations / 10 iterations / 100 
C         iterations...

          IF ( (IAND(IGSOUT,4).NE.0) .OR.
     *     ( (IAND(IGSOUT,8).NE.0)  .AND. (MOD(ISMTH,10).EQ.0)  ) .OR.
     *     ( (IAND(IGSOUT,16).NE.0) .AND. (MOD(ISMTH,100).EQ.0) ) .OR.
     *     ( (IAND(IGSOUT,32).NE.0) .AND. 
     *       (ISMTH.EQ.IPARAM(ONPRGSM)) )) THEN

            CALL GTMAUX (TIMIN,DPARAM,OTGGSM,1)
     
C           Recalculate element statistics for output.
C           The element statistics are recalculated later on anyway
C           and no more used, so it's save to overwrite them here
C           with the current values.
     
            CALL MEASR2(DWORK(L(TRIAS(OLCORVG,ILEV))),
     *              KWORK(L(TRIAS(OLVERT,ILEV))),
     *              KWORK(L(TRIAS(OLADJ,ILEV))),NEL,
     *              NVT,KWORK(L(LCOUNT)),DWORK(L(LSIZE)),
     *              DWORK(L(LSIZEM)),DWORK(L(LRATIO)))

            CALL GAGMVS (CTRIA,-1,ILEV,ITER,ISMTH,
     *                   TRIAS(ONEL,NLMIN),LCOLOR,LSIZE,LSIZEM,LRATIO,
     *                   0,0,0,0,0,
     *                   'gmv/SMPOS')    

            CALL GTMAUX (TIMIN,DPARAM,OTGGSM,0)
     
           END IF

        END DO

        IF ((MT.GE.2).AND.(IPARAM(ONPOGSM).NE.0).AND.
     *      (MOD(IPARAM(ONPRGSM),20).EQ.1)) THEN
          WRITE (MTERM,'(A,I5)') 'Grid-smoothing iteration ',ISMTH
        END IF
        
        CALL GTMAUX (TIMIN,DPARAM,OTGGSM,1)
        
C       Ok, at this point the whole grid was deformed. 
C       Before we can go to the next macrostep, we have to do some
C       postprocessing.
C
C       For the time being, we have the grid in CTRIA. Because the
C       grid adaption with that grid is completed now, we can
C       throw away our old grid and use the new one:

        CALL LCP1 (DWORK(L( CTRIA(OLCORVG) )),
     *             DWORK(L( TRIAS(OLCORVG,ILEV) )),2*NVT)
        CALL LCP1 (DWORK(L( CTRIA(OLVBDP) )),
     *             DWORK(L( TRIAS(OLVBDP,ILEV) )),TRIAS(ONVBD,ILEV))

C       Correct midpoint information on all finer levels
C       if the current level is less than the solver level.
C       If the current level is larger than the solver level,
C       the sharing of the coordinates will automatically ensure
C       that the coarse grid coordinates are corrected.

        DO I = ILEV,ILVSLV-1
          CALL PRC2FG (DWORK(L( TRIAS(OLCORVG,I) )),
     *                 KWORK(L( TRIAS(OLVERT,I) )),
     *                 TRIAS(ONEL,I),
     *                 DWORK(L( TRIAS(OLCORVG,I+1) )),
     *                 KWORK(L( TRIAS(OLVERT,I+1) )),
     *                 KWORK(L( TRIAS(OLADJ,I+1) )),
     *                 TRIAS(ONEL,I+1))
        END DO

C       At the end of the macrostep we want to print out our 
C       quality measure. Therefore we have to recalculate some
C       statistical information...

        CALL GTMAUX (TIMIN,DPARAM,OTGGAPR,0)

C       Write out the final grid of this macrostep into a GMV-file
C       if desired.

        IF ( IAND(IGSOUT,65).NE.0 ) THEN
    
C         Calculate area distribution to update statistical data

          CALL MEASR2(DWORK(L(TRIAS(OLCORVG,ILEV))),
     *                KWORK(L(TRIAS(OLVERT,ILEV))),
     *                KWORK(L(TRIAS(OLADJ,ILEV))),NEL,
     *                NVT,KWORK(L(LCOUNT)),DWORK(L(LSIZE)),
     *                DWORK(L(LSIZEM)),DWORK(L(LRATIO)))

C         Calculate the monitor function on the current grid.
C         DMON is large enough to hold the information not only
C         on the "solution" level, but also on the level which is
C         deformed now.
         
          CALL XMON0(TRIAS(1,ILEV),LSIZE,LMON,1,DPARAM(OEPS0GS),
     *               IGEOM,DGEOM,IDATA,DDATA)
          
          CALL GAGMVS (TRIAS(1,ILEV),-1,ILEV,ITER,0,
     *                 TRIAS(ONEL,NLMIN),LCOLOR,LSIZE,LSIZEM,LRATIO,
     *                 LMON,0,0,0,0,
     *                 'gmv/MACRO')    
     
        END IF

C       Calculate area distribution on solver level

        CALL MEASR2(DWORK(L(TRIAS(OLCORVG,ILVSLV))),
     *              KWORK(L(TRIAS(OLVERT,ILVSLV))),
     *              KWORK(L(TRIAS(OLADJ,ILVSLV))),TRIAS(ONEL,ILVSLV),
     *              TRIAS(ONVT,ILVSLV),KWORK(L(LCOUNT)),
     *              DWORK(L(LSIZE)),
     *              DWORK(L(LSIZEM)),DWORK(L(LRATIO)))
     
C       Monitor function on the solver level

        CALL XMON0(TRIAS(1,ILVSLV),LSIZE,LMON,1,DPARAM(OEPS0GS),
     *             IGEOM,DGEOM,IDATA,DDATA)
          
C       Domain measure

        CALL MEASURE(NEL,DWORK(L(LSIZEM)),DMSR)
          
C       And finally: Calculate the quality measure:

        CALL GAQMSR (NVT,DWORK(L(LMON)),DWORK(L(LSIZE)),DMSR,D)
        
        CALL GTMAUX (TIMIN,DPARAM,OTGGAPR,1)
        
        IF (MT.GE.2) THEN
          WRITE (MTERM,9000) 
          WRITE (MTERM,'(A,F10.5)') 'Current quality measure: ',D
          WRITE (MTERM,9000) 
        END IF
        
        IF ((IPARAM(OIARCTR).NE.0).AND.
     *      (DPARAM(ODMNTOL).GT.0D0).AND.(D.LE.DPARAM(ODMNTOL))) THEN

C         Stop the grid adaption if the quality measure drops below
C         the desired value.

          GOTO 110

        END IF

      END DO
      
110   CONTINUE
      
C     For further output purposes reset ITER to maximum
      
      IF (ITER.GT.NGASTP) ITER = NGASTP+NGAPCS
      
C     Grid adaption finished. We need some postprocessing...
      
      IF (IAND(IGSOUT,64).NE.0) THEN

C       Calculate area distribution to update statistical data

        CALL MEASR2(DWORK(L(TRIAS(OLCORVG,ILEV))),
     *              KWORK(L(TRIAS(OLVERT,ILEV))),
     *              KWORK(L(TRIAS(OLADJ,ILEV))),NEL,
     *              NVT,KWORK(L(LCOUNT)),DWORK(L(LSIZE)),
     *              DWORK(L(LSIZEM)),DWORK(L(LRATIO)))

        CALL GAGMVS (CTRIA,-1,ILEV,ITER,0,
     *               TRIAS(ONEL,NLMIN),LCOLOR,LSIZE,LSIZEM,LRATIO,LMON,
     *               0,0,0,0,
     *               'gmv/FINAL')    

      END IF
      
C     The following jump mark will be used in case of an error in the
C     grid adaption:
      
99990 CONTINUE      
      
C     Release memory that was used for the solving process of the 
C     linear systems: Handles in MGPARI

      CALL DISGS2(IPARAM(SZGACI+1), DPARAM(SZGACD+1), 
     *            MGPARI, MGPARD)
      
C     Release memory of solution/gradient,monitor function,...

      CALL ZDISP (0,LCOLOR,'KCOLOR')
      CALL ZDISP (0,LRATIO,'DRATIO')
      CALL ZDISP (0,LSIZEM,'DSIZEM')
      CALL ZDISP (0,LCOUNT,'KCOUNT')
      CALL ZDISP (0,LSIZE,'DSIZE ')
      CALL ZDISP (0,LMNOLD,'DMNOLD')
      CALL ZDISP (0,LMON,'DMON  ')
      CALL ZDISP (0,LTMP,'DTMP  ')
      CALL ZDISP (0,LGRADY,'GRADY ')
      CALL ZDISP (0,LGRADX,'GRADX ')
      CALL ZDISP (0,LU,'DU    ')
      
C     Release extended triangulation information from CTRIA.
C     Don't release anything if the structure was converted on entry
C     of this routine, as we then would delete information
C     generated by the caller!

      IF (IAND(TRIFL,1).EQ.0) CALL DISETR (TRIAS(1,ILEV),0)
      
C     Release memory of DCORVG/DVBDP in CTRIA

      CALL ZDISP (0,CTRIA(OLVBDP),'LVBDP2')
      CALL ZDISP (0,CTRIA(OLCORVG),'DCORV2')
            
      IF (MT.GE.2) THEN
        WRITE (MTERM,'(A)') '======================================='
        WRITE (MTERM,'(A)') 'Grid deformation finished.'
        WRITE (MTERM,'(A)') '======================================='
        WRITE (MTERM,'(A)')
      END IF

C     Calculate total time

      CALL GTMAUX (TIMIN,DPARAM,OTGTOT,1)

9000  FORMAT(30('-'))

      END 

************************************************************************
* Static grid adaption with correction of multiple levels
*
* This subroutine performs grid-adaption according to the parameter
* set. It accepts a set of levels, which are expected to be formed
* by regular refinement of the coarsest grid. The grid adaption
* takes place on level ILEV. After the grid adaption the coarser levels 
* are corrected in a post-correction step to fit again to the data on 
* the finest grid. If the grid adaption doesn't take place on the
* finest grid, the finer grids are obtained from the coarser ones
* by prolongation.
*
* In:
*   TRIAS  - array [1..NNLEV] of STRIA
*            = array [1..SZTRIA,1..NNLEV] of integer
*            Set of triangulations for all possible levels
*   NLMIN  - Minimum existing level in TRIAS
*   NLMAX  - Maximum existing level in TRIAS
*   ILEV   - Actual level of the grid adaption.
*            NLMIN <= ILEV <= NLMAX
*   ILVSLV - Level where the linear system should be solved.
*            NLMIN <= ILEV <= NLMAX
*   IPARAM - array [1..SZGACI+SZGLSI+SZGASI] of integer
*            Integer parameter-array determining the behaviour of the
*            routine.
*            This parameter array must contain the following
*            three blocks, right after another:
*              Index [ 1 .. SZGACI ]
*                    TIGridAdaptionCommon-structure
*              Index [ SZGACI+1 .. SZGACI+SZGLSI ]
*                    TIGridAdaptionSolver-structure
*              Index [ SZGACI+SZGLSI+1 .. SZGACI+SZGLSI+SZGASI ]
*                    TIGridAdaptionStaticParams-structure
*            
*   DPARAM - array [1..*] of double precision
*            Double precision parameter-array determining the
*            behaviour of the routine.
*            This parameter array must contain the following
*            three blocks, right after another:
*              Index [ 1 .. SZGACD ]
*                    TDGridAdaptionCommon-structure
*              Index [ SZGACD+1 .. SZGACD+SZGLSD ]
*                    TDGridAdaptionSolver-structure
*              Index [ SZGACD+SZGLSD+1 .. SZGACD+SZGLSD+SZGASD ]
*                    TDGridAdaptionStaticParams-structure
* 
*  IGEOM  - array [1..*] of integer 
*  DGEOM  - array [1..*] of double 
*           Integer- and double-precision parameter blocks with
*           geometry information. Passed to boundary
*           routines. Not used in this routine.
*  IDATA  - array [1..*] of integer 
*  DDATA  - array [1..*] of double 
*           User defined integer- and double-precision parameter 
*           blocks. Passed to callback routines like monitor
*           function. Not used in this routine.
*
* Out:
*   TRIAS  - adapted grids on all levels
*
* The parameter blocks must be initialized properly, e.g. with INGAST.
************************************************************************
      
      SUBROUTINE XSMMGS (TRIAS,NLMIN,NLMAX,ILEV,ILVSLV,IPARAM,DPARAM,
     *                   IGEOM,DGEOM,IDATA,DDATA)
      
      IMPLICIT NONE

      
      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'stria.inc'
      
      INCLUDE 'sgridadaptgeneralparams.inc'
      INCLUDE 'sgridadaptsolverparams.inc'
      INCLUDE 'sgridadaptstaticparams.inc'
      
      INCLUDE 'cbasicmg.inc'
      
      INTEGER TRIAS(SZTRIA,*),IGEOM(*),IDATA(*)
      INTEGER IPARAM(SZGACI+SZGLSI+SZGASI)
      DOUBLE PRECISION DPARAM(SZGACD+SZGLSD+SZGASD),DGEOM(*),DDATA(*)
      DOUBLE PRECISION TIMIN(SZGACD)
      INTEGER NLMIN,NLMAX,ILEV,ILVSLV

C     local variables

      INTEGER LV,OLMILV,I,KVBD,STARTI,LVSLV
      INTEGER KVERT ,KCORVG ,KADJ ,KMID ,KCORMG ,KBCT ,KVBDP ,NEL ,NMT
      INTEGER KMBDP ,NVT ,NBCT
      INTEGER KVERT2,KCORV2 ,KADJ2,KMID2,KCORM2 ,KBCT2,KVBDP2,NEL2,NMT2
      INTEGER KMBDP2,NVT2,NBCT2

C     Clear the timing-statistic parameters we set in this routine:

      DPARAM(OTGTOT) = 0D0
      DPARAM(OTGGPR) = 0D0

C     Save the start time - in our local timing-block!

      CALL GTMAUX (TIMIN,DPARAM,OTGTOT,0)

C     Depending on the parameters in the parameter list set the actual
C     level LV where the adaption should be performed:

      LV = ILEV
      IF (LV.GT.NLMAX) LV = NLMAX
      IF (LV.LT.NLMIN) LV = NLMIN

      LVSLV = ILVSLV
      IF (LVSLV.GT.NLMAX) LVSLV = NLMAX
      IF (LVSLV.LT.NLMIN) LVSLV = NLMIN

C     The minimum level is to be found in the solver-substructure
C     of the parameter-array, starting at position SZGACI.

      STARTI = SZGACI+1-1
      
C     Transfer minimum level to structure - because the value in the
C     structure must be between NLMIN and NLMAX in order for the
C     solver to work! The old minimum level is restored later.

      OLMILV = IPARAM(STARTI+OICGSLV)

      IF (IPARAM(STARTI+OICGSLV).LT.0)  
     *  IPARAM(STARTI+OICGSLV) = LVSLV+IPARAM(STARTI+OICGSLV)
      IF (IPARAM(STARTI+OICGSLV).LT.NLMIN)
     *  IPARAM(STARTI+OICGSLV) = NLMIN
      IF (IPARAM(STARTI+OICGSLV).GT.LVSLV)   
     *  IPARAM(STARTI+OICGSLV) = LV

C     Call the static grid adaption for grid on level LV
      
      CALL XGASTA (NLMIN,NLMAX,TRIAS,LV,LVSLV,
     *             IPARAM, DPARAM,IGEOM,DGEOM,IDATA,DDATA)

C     Restore the old minimum level in the parameter block immediately

      IPARAM(STARTI+OICGSLV) = OLMILV
        
C     Stop time of grid prol/rest

      CALL GTMAUX (TIMIN,DPARAM,OTGGPR,0)
      
C     Recalculate the information about the higher levels if necessary

      DO I = LV,NLMAX-1
          
C       Transfer variables for easier access:

        NVT    = TRIAS(ONVT    ,I)
        NMT    = TRIAS(ONMT    ,I)
        NEL    = TRIAS(ONEL    ,I)
        NBCT   = L(TRIAS(ONBCT   ,I))
        KCORVG = L(TRIAS(OLCORVG ,I))
        KVERT  = L(TRIAS(OLVERT  ,I))
        KCORMG = L(TRIAS(OLCORMG ,I))
        KMID   = L(TRIAS(OLMID   ,I))
        KADJ   = L(TRIAS(OLADJ   ,I))
        KBCT   = L(TRIAS(OLBCT   ,I))
        KVBDP  = L(TRIAS(OLVBDP  ,I))
        KMBDP  = L(TRIAS(OLMBDP  ,I))

        NVT2    = TRIAS(ONVT    ,I+1)
        NMT2    = TRIAS(ONMT    ,I+1)
        NEL2    = TRIAS(ONEL    ,I+1)
        NBCT2   = L(TRIAS(ONBCT   ,I+1))
        KCORV2  = L(TRIAS(OLCORVG ,I+1))
        KVERT2  = L(TRIAS(OLVERT  ,I+1))
        KCORM2  = L(TRIAS(OLCORMG ,I+1))
        KMID2   = L(TRIAS(OLMID   ,I+1))
        KADJ2   = L(TRIAS(OLADJ   ,I+1))
        KBCT2   = L(TRIAS(OLBCT   ,I+1))
        KVBDP2  = L(TRIAS(OLVBDP  ,I+1))
        KMBDP2  = L(TRIAS(OLMBDP  ,I+1))
          
C       Prolongate corner vertices:

        CALL PRC2FG (DWORK(KCORVG),KWORK(KVERT),
     *               NEL,DWORK(KCORV2),KWORK(KVERT2),
     *               KWORK(KADJ2),NEL2)          
     
C       Recalculate edge midpoint data from corners

        CALL RXMPTS (DWORK(KCORV2),KWORK(KVERT2),
     *               KWORK(KMID2),NEL2,NVT2,NMT2,
     *               DWORK(KCORM2))
     
C       Prolongate boundary points

        CALL PRC2FP (KWORK(KBCT),DWORK(KVBDP),NBCT,
     *               KWORK(KBCT2),DWORK(KVBDP2))
     
C       Recalculate boundary midpoints

        CALL RXBMPR (KWORK(KBCT2),DWORK(KVBDP2),
     *               DWORK(KMBDP2),NBCT)

      END DO
        
C     Ok, so far so good, the grid adaption on the higher levels is 
C     finished. But we are not finished yet! We adapted all meshes, 
C     but the grid adaption on finer grids affected the grids on lower 
C     levels, too. 
C     Therefore the parameter values of the boundary nodes of lower levels 
C     generally have changed! The information on the finer grids is 
C     appropriate, but that on coarser grids no more :( But it's necessary 
C     for the computation!
C     So we have no other chance to again loop through all levels and
C     recalculate the boundary information :(

      DO I=NLMIN,NLMAX
      
C       Recalculate midpoint data - important for nonconforming 
C       elements!!! Especially for all drag-/lift postprocessing 
C       routines, which rely on this...
C
C       Transfer variables for easier access:

        NVT    = TRIAS(ONVT    ,I)
        NMT    = TRIAS(ONMT    ,I)
        NEL    = TRIAS(ONEL    ,I)
        NBCT   = L(TRIAS(ONBCT   ,I))
        KCORVG = L(TRIAS(OLCORVG ,I))
        KVERT  = L(TRIAS(OLVERT  ,I))
        KMID   = L(TRIAS(OLMID   ,I))
        KADJ   = L(TRIAS(OLADJ   ,I))
        KBCT   = L(TRIAS(OLBCT   ,I))
        KVBD   = L(TRIAS(OLVBD   ,I))
        KVBDP  = L(TRIAS(OLVBDP  ,I))

        IF (TRIAS(OLCORMG ,I).NE.0) THEN
        
C         Recalculate edge midpoint data from corners

          KCORMG = L(TRIAS(OLCORMG ,I))
          CALL RXMPTS (DWORK(KCORVG),KWORK(KVERT),
     *                 KWORK(KMID),NEL,NVT,NMT,
     *                 DWORK(KCORMG))
        END IF
          
C       Correct the boundary parameter values:

        CALL RXBNPC(DWORK(KCORVG),KWORK(KVBD),
     *              KWORK(KBCT),DWORK(KVBDP))

C       Recalculate boundary midpoints

        IF (TRIAS(OLMBDP  ,I).NE.0) THEN
          KMBDP  = L(TRIAS(OLMBDP  ,I))
          CALL RXBMPR (KWORK(KBCT),DWORK(KVBDP),
     *                 DWORK(KMBDP),NBCT)
        END IF

C       Recalculate area-array of the cells

        IF (TRIAS(OLAREA  ,I).NE.0) THEN
          CALL GENET6(TRIAS(1,I))
        END IF
          
      END DO

C     Time of grid prol/rest

      CALL GTMAUX (TIMIN,DPARAM,OTGGPR,1)
      
C     Calculate total time, overwrite time calculated by subroutine

      DPARAM(OTGTOT) = 0D0
      CALL GTMAUX (TIMIN,DPARAM,OTGTOT,1)

      END 
