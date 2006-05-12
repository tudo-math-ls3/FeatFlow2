************************************************************************
* This file contains various optimization algorithms to optimize
* a value of interest.
************************************************************************

************************************************************************
* Optimization algorithm 00: Brute-force search
*
* This routine implements brute force search
*
* It's assumed that there is a calculation routine DCALC which 
* represents a mapping
*     DCALC: R^IDIMS -> R
* i.e. that maps the IDIMS dimensional design variable space R^IDIMS
* to a scalar value. The optimization algorithm tries to find a
* parameter setting x such that
*     DCALC(x) -> min!
* For that purpose DCALC is frequently called and the results are
* then used to find a better parameter setting x.
*
* DSTART specified the starting point where to start the search
* in R^IDIMS. DSTEP is the stepsize in each direction.
* For every dimension, DSTEP contains the stepsize parameter. 
* DEND specifies for each direction the maximum allowed parameter
* value. 
*
* Brute-force search loops through all parameter values from DSTART
* to <= DEND in the step sizes according to DSTEP. The minimum
* calculated fittness-value and the corresponding parameter
* setting is returned.
*
* In:
*   IDIMS  - integer
*            Number of dimensions / variables to optimize a value of
*            interest for.
*   DSTART - array [1..IDIMS] of double
*            Start vector where to begin the search.
*   DSTEP  - array [1..IDIMS] of double
*            For every dimension: initial resolution / stepsize.
*            If the step-length parameter of a dimension is <=0,
*            this dimension is not tested.
*   DEND   - array [1..IDIMS] of double
*            For each dimension: maximum allowedparameter value
*   IINFO  - array [1..*] of integer
*            Integer information block for the probing routine;
*            is directly passed to DPROBE on call without change.
*   DINFO  - array [1..*] of double
*            Integer information block for the probing routine;
*            is directly passed to DPROBE on call without change.
*   NRES   - integer
*            Size of double-precision array that will be used as
*            temporary space for saving the results of each
*            probing; see DPROBE.
*   DPROBE - SUBROUTINE (IDIMS,LCOORDS,NRES,LRES,IINFO,DINFO)
*            User defined probing function to calculate the value 
*            of interest.
*            In:
*              IDIMS   - Number of dimensions of the design parameters
*              LCOORD  - Handle to array [1..IDIMS] of double
*                        Design parameter setting where to evaluate
*              NRES    - integer
*                        Size of result array
*              LRES    - integer
*                        Handle to array [1..NRES] of double,
*                        where to write the results to
*              PITER   - Number of current probing iterate
*              IINFO,
*              DINFO   - Integer / double precision parameter block
*                        from above; passed to DPROBE unchanged.
*            Out:
*              RES     - array [1..NRES] of double
*                        RES(1) is the value of the point DCOORDS,
*                        which is used by the optimizer to
*                        optimize for. RES(2..NRES) can be set by
*                        DPROBE to arbitrary values that are later
*                        passed to DNXSTP as result values.
*            Remark: The RES-array and the COORD-array are identified by 
*             handles, not by starting addresses. The reason is that
*             the DPROBE-routine is assumed to reassign many handles
*             on the heap, thus destroying a possible starting address
*             of a dynamically allocated array for the results on the
*             heap. By passing the handle to that result array instead
*             of the result array itself, DPROBE will be able to write
*             the results to the correct address.
*     
*            It may happen that the parameter setting is outside
*            of a valid domain. In this case DPROBE must return
*            1D99 to indicate this.
*
*            DPROBE and DNXSTP always belong together. It's acceptable
*            do do memory allocation + calculation in DPROBE, while
*            doing postprocessing and releasing of memory in DNXSTP!
*
*   DNXSTP - SUBROUTINE (IDIMS,DCOORDS,DCOOPT,IINFO,DINFO,
*                        PITER,STEP,IACPT,NRES,DRES,DRSOPT)
*            This subroutine is called after each probing to inform the
*            caller whether a new iterate was accepted or not
*            In:
*              DCOORDS - array [1..IDIMS] of double precision
*                        Parameter set of probed iterate.
*              DCOOPT  - array [1..IDIMS] of double precision
*                        Parameter set of current optimal iterate
*              IINFO,
*              DINFO   - Integer / double precision parameter
*                        block from above
*              PITER   - Number of current probing iterate
*              STEP    - Number of current accepted iterate
*              IACPT   - whether DCOORDS will be accepted as the new
*                        iterate.
*                        =0: DCOORDS will not be accepted
*                        =1: DCOORDS will be accepted
*                        Can be modified by DNXSTP to prevent the
*                        point from being accepted.
*              NRES    - integer
*                        Size of result array RES
*              DRES    - Result array that was calculated 
*                        in DPROBE for the parameter setting DCOORDS.
*              DRSOPT  - Result array that was calculated 
*                        in DPROBE for the parameter setting DCOOPT.
*
* Out:
*   DFINAL - array [1..IDIMS] of double
*            best configuration
*   STEPS  - integer
*            Number of steps the algorithm used.
************************************************************************

      SUBROUTINE OPTA00 (IDIMS,DSTART,DSTEP,DEND,
     *                   IINFO,DINFO,NRES,DPROBE,DNXSTP,
     *                   DFINAL,STEPS)
     
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
C     parameters

      INTEGER IDIMS,STEPS,IINFO(*),NRES
      DOUBLE PRECISION DSTART(IDIMS),DSTEP(IDIMS),DEND(IDIMS)
      DOUBLE PRECISION DINFO(*),DFINAL(IDIMS)
      
      EXTERNAL DPROBE,DNXSTP
      
C     local variables

      INTEGER STEP,CNT
      INTEGER LCOORDS,LRES,LRESTM,J,I
      DOUBLE PRECISION NEWVAL,VAL
      LOGICAL BSTOP
      
C     Reserve some memory for the current parameter setting:

      CALL ZNEW (IDIMS,-1,LCOORDS,'DCOORDS')
      
C     We work with parameters in DCOORDS; copy start point
C     from DSTART:

      CALL LCP1 (DSTART,DWORK(L(LCOORDS)),IDIMS)
      
C     And for storing the results. 

      CALL ZNEW (NRES,-1,LRES,'DRES  ')
      CALL ZNEW (NRES,-1,LRESTM,'DRESTM')
      
C     Let's go...
      
      STEP = 1
      CNT = 0
      VAL = 1D99
      
C     Check if the initial parameter values are out-of-range - 
C     in this case do not even start.
C     (We slightly enlarge the upper bound because some Fortran
C     compilers have problems with comparing doubles)

      BSTOP = .TRUE.
      DO I=1,IDIMS 
        BSTOP = BSTOP.AND.(DSTART(I).GT.(DEND(I)+1D-10))
      END DO
      
      DO WHILE (.NOT.BSTOP)
      
C       Increase probing counter

        CNT = CNT + 1
       
C       Calculate the value of interest:

        CALL DPROBE (IDIMS,LCOORDS,NRES,LRESTM,CNT,
     *               IINFO,DINFO)
     
C       Result in first variable of 2nd part of RES-array:
     
        NEWVAL = DWORK(L(LRESTM))
       
C       Is it better? -> Change the iterate...

        I = 0
        IF (NEWVAL.LT.VAL) I = 1
       
C       Indicate the framework whether we will accept the iterate
C       or not

        CALL DNXSTP (IDIMS,DWORK(L(LCOORDS)),DFINAL,
     *               IINFO,DINFO,CNT,STEP,I,NRES,
     *               DWORK(L(LRESTM)),DWORK(L(LRES)))
     
C       Are we allowed to accept it? If yes, do that!
      
        IF (I.NE.0) THEN

C         Accept the new iterate
        
          CALL LCP1(DWORK(L(LCOORDS)),DFINAL,IDIMS)
          CALL LCP1(DWORK(L(LRESTM)),DWORK(L(LRES)),NRES)
          VAL = NEWVAL
          STEP = STEP + 1
        
        END IF

C       Increase the parameter setting according to DSTEP.
C       Search the first parameter that can be increased by 
C       DSTEP(.).
C       (We slightly enlarge the upper bound because some Fortran
C       compilers have problems with comparing doubles)
      
        I=1
        DO WHILE 
     *    ( ( ((DWORK(L(LCOORDS)+I-1)+DSTEP(I)).GT.(DEND(I)+1D-10)) .OR.
     *         (DSTEP(I).LE.0D0) ) .AND.
     *      (I.LE.IDIMS) ) 
          I=I+1
        END DO
        
C       Check if all configurations appeared:
        
        BSTOP = I.GT.IDIMS
        
        IF (.NOT.BSTOP) THEN
        
C         Increase the next possible parameter value and reset
C         the others.
  
          DWORK(L(LCOORDS)+I-1) = DWORK(L(LCOORDS)+I-1) + DSTEP(I)
          
          DO J=1,I-1
            DWORK(L(LCOORDS)+J-1) = DSTART(J)
          END DO
          
C         Go to the next step...

        END IF
      
      END DO

C Finally release memory:

      CALL ZDISP (0,LRESTM,'DRESTM')
      CALL ZDISP (0,LRES,'DRES  ')
      CALL ZDISP (0,LCOORDS,'DCOORDS')

C Set statistical data, finish

      STEPS = STEP

      END
      
************************************************************************
* Optimization algorithm 01: Compass-Search
*
* This routine implements basic compass search.
*
* It's assumed that there is a calculation routine DCALC which 
* represents a mapping
*     DCALC: R^IDIMS -> R
* i.e. that maps the IDIMS dimensional design variable space R^IDIMS
* to a scalar value. The optimization algorithm tries to find a
* parameter setting x such that
*     DCALC(x) -> min!
* For that purpose DCALC is frequently called and the results are
* then used to find a better parameter setting x.
*
* DSTART specified the starting point where to start the search.
* HSTART is an initial resolution / stepsize parameter array, i.e. 
* specifies the size of the domain around DSTART where to search for 
* better points. For every dimension, HSTART contains the initial 
* stepsize parameter. EPS specifies a stopping criterion.
*
* Compass search stops the search if the maximum of the stepsizes
* in all dimensions drop belop EPS, or if > MAXSTP optimization steps 
* have been performed.
*
* For a description of the algorithm, see:
* [Kolda, T.G.; Lewis, R.M.; Torczon, V.; Optimization by Direct
*  Search: New Perspectives on Some Classical and Modern Methods;
*  SIAM REVIEW, Vol. 45, No. 3, 385-482]
*
* In:
*   IDIMS  - integer
*            Number of dimensions / variables to optimize a value of
*            interest for.
*   DSTART - array [1..IDIMS] of double
*            Start vector where to begin the search.
*   HSTART - array [1..IDIMS] of double
*            For every dimension: initial resolution / stepsize
*   EPS    - double
*            Stopping criterion. Algorithm stops if all stepsizes
*            (which arise from modifying HSTART) drop below EPS.
*   MAXSTP - Maximum number of optimization steps
*            =-1: not specified.
*   IINFO  - array [1..*] of integer
*            Integer information block for the probing routine;
*            is directly passed to DPROBE on call without change.
*   DINFO  - array [1..*] of double
*            Integer information block for the probing routine;
*            is directly passed to DPROBE on call without change.
*   NRES   - integer
*            Size of double-precision array that will be used as
*            temporary space for saving the results of each
*            probing; see DPROBE.
*   DPROBE - SUBROUTINE (IDIMS,LCOORDS,NRES,LRES,IINFO,DINFO)
*            User defined probing function to calculate the value 
*            of interest.
*            In:
*              IDIMS   - Number of dimensions of the design parameters
*              LCOORD  - Handle to array [1..IDIMS] of double
*                        Design parameter setting where to evaluate
*              NRES    - integer
*                        Size of result array
*              LRES    - integer
*                        Handle to array [1..NRES] of double,
*                        where to write the results to
*              PITER   - Number of current probing iterate
*              IINFO,
*              DINFO   - Integer / double precision parameter block
*                        from above; passed to DPROBE unchanged.
*            Out:
*              RES     - array [1..NRES] of double
*                        RES(1) is the value of the point DCOORDS,
*                        which is used by the optimizer to
*                        optimize for. RES(2..NRES) can be set by
*                        DPROBE to arbitrary values that are later
*                        passed to DNXSTP as result values.
*            Remark: The RES-array and the COORD-array are identified by 
*             handles, not by starting addresses. The reason is that
*             the DPROBE-routine is assumed to reassign many handles
*             on the heap, thus destroying a possible starting address
*             of a dynamically allocated array for the results on the
*             heap. By passing the handle to that result array instead
*             of the result array itself, DPROBE will be able to write
*             the results to the correct address.
*     
*            It may happen that the parameter setting is outside
*            of a valid domain. In this case DPROBE must return
*            1D99 to indicate this.
*
*            DPROBE and DNXSTP always belong together. It's acceptable
*            do do memory allocation + calculation in DPROBE, while
*            doing postprocessing and releasing of memory in DNXSTP!
*
*   DNXSTP - SUBROUTINE (IDIMS,DCOORDS,DCOOPT,IINFO,DINFO,
*                        PITER,STEP,IACPT,NRES,DRES,DRSOPT)
*            This subroutine is called after each probing to inform the
*            caller whether a new iterate was accepted or not
*            In:
*              DCOORDS - array [1..IDIMS] of double precision
*                        Parameter set of probed iterate.
*              DCOOPT  - array [1..IDIMS] of double precision
*                        Parameter set of current optimal iterate
*              IINFO,
*              DINFO   - Integer / double precision parameter
*                        block from above
*              PITER   - Number of current probing iterate
*              STEP    - Number of current accepted iterate
*              IACPT   - whether DCOORDS will be accepted as the new
*                        iterate.
*                        =0: DCOORDS will not be accepted
*                        =1: DCOORDS will be accepted
*                        Can be modified by DNXSTP to prevent the
*                        point from being accepted.
*              NRES    - integer
*                        Size of result array RES
*              DRES    - array [1..NRES] of double precision
*                        Result array that was calculated 
*                        in DPROBE for the parameter setting DCOORDS.
*              DRSOPT  - array [1..NRES] of double precision
*                        Result array that was calculated 
*                        in DPROBE for the parameter setting DCOOPT.
*
* Out:
*   DFINAL - array [1..IDIMS] of double
*            Final parameter setting where the search was stopped.
*   FINTOL - array [1..IDIMS] of double
*            Final step size when the algorithm was stopped
*   STEPS  - integer
*            Number of steps the algorithm used.
************************************************************************

      SUBROUTINE OPTA01 (IDIMS,DSTART,HSTART,EPS,MAXSTP,
     *                   IINFO,DINFO,NRES,DPROBE,DNXSTP,
     *                   DFINAL,FINTOL,STEPS)
     
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
C parameters

      INTEGER IDIMS,MAXSTP,STEPS,IINFO(*),NRES
      DOUBLE PRECISION DSTART(IDIMS),HSTART(IDIMS),DFINAL(IDIMS),EPS
      DOUBLE PRECISION DINFO(*)
      DOUBLE PRECISION FINTOL(IDIMS)
      
      EXTERNAL DPROBE,DNXSTP
      
C local variables

      INTEGER STEP,CNT
      INTEGER LCOORDS,LRES,LRESTM,CDIM,J,I
      DOUBLE PRECISION NEWVAL,VAL,H
      
C We work with parameters in DFINAL; copy start point from DSTART:

      CALL LCP1(DSTART,DFINAL,IDIMS)
      
C Current stepsize in every dimension

      CALL LCP1(HSTART,FINTOL,IDIMS)
      
C Reserve some memory for the current parameter setting:

      CALL ZNEW (IDIMS,-1,LCOORDS,'DCOORDS')
      
C And for storing the results. 

      CALL ZNEW (NRES,-1,LRES,'DRES  ')
      CALL ZNEW (NRES,-1,LRESTM,'DRESTM')
      
C Let's go...
      
      STEP = 1
      CNT = 1
      
C Calculate the value of the function in the starting point:

      CALL LCP1(DSTART,DWORK(L(LCOORDS)),IDIMS)
      CALL DPROBE (IDIMS,LCOORDS,NRES,LRES,CNT,IINFO,DINFO)

C     Result in first variable of RES-array
      
      VAL = DWORK(L(LRES))
      
      I = 1
      CALL DNXSTP (IDIMS,DWORK(L(LCOORDS)),DWORK(L(LCOORDS)),
     *             IINFO,DINFO,CNT,STEP,I,NRES,
     *             DWORK(L(LRES)),DWORK(L(LRES)))
      
C Crippled DO-WHILE loop to realize a maximum of MAXSTP steps:
      
10    CONTINUE
      
C       Loop over all possible dimensions:

        DO CDIM = 0,IDIMS-1
        
C         In every dimension: 1x forward, 1x backward
        
          DO J=0,1
          
C           Increase probing counter

            CNT = CNT + 1
        
C           Copy our current point to the probe coordinate array

            CALL LCP1(DFINAL,DWORK(L(LCOORDS)),IDIMS)
          
C           Modify the CDIM'th entry according to the stepsize
C           to create a neighbour point of our current midpoint:

            H = FINTOL(CDIM+1)
            
            DWORK(L(LCOORDS)+CDIM) = DWORK(L(LCOORDS)+CDIM) 
     *                               + H - DBLE(J)*2D0*H

            IF (H.NE.0D0) THEN
           
C             Calculate the value of interest:

              CALL DPROBE (IDIMS,LCOORDS,NRES,LRESTM,CNT,
     *                     IINFO,DINFO)
     
C             Result in first variable of 2nd part of RES-array:
     
              NEWVAL = DWORK(L(LRESTM))
          
C             Is it better? -> Change the iterate...

              I = 0
              IF (NEWVAL.LT.VAL) I = 1
          
C             Indicate the framework whether we will accept the iterate
C             or not

              CALL DNXSTP (IDIMS,DWORK(L(LCOORDS)),DFINAL,
     *                     IINFO,DINFO,CNT,STEP,I,NRES,
     *                     DWORK(L(LRESTM)),DWORK(L(LRES)))
     
            ELSE
            
C             H=0 is not allowed, we don't compute anything in this
C             step!            
            
              I = 0
            
            END IF
     
C           Are we allowed to accept it? If yes, do that!
      
            IF (I.NE.0) THEN

C             Accept the new iterate: Remember parameter setting and 
C             value
            
              CALL LCP1(DWORK(L(LCOORDS)),DFINAL,IDIMS)
              CALL LCP1(DWORK(L(LRESTM)),DWORK(L(LRES)),NRES)
              VAL = NEWVAL
              STEP = STEP + 1
            
              GOTO 80000

            END IF
        
C         Otherwise check next dimension...

          END DO
        
        END DO
      
C       All dimensions testet, no one was better. Halfen the step-
C       length and calculate the maximum step length

        H = 0.5D0 * FINTOL(1)
        DO J=1,IDIMS
          FINTOL(J) = 0.5D0 * FINTOL(J)
          IF (FINTOL(J).GT.H) H = FINTOL(J)
        END DO
        
C       Stop the algorithm if our current maximum step-length drops
C       below the EPS-value...

        IF (H.LT.EPS) GOTO 90000
      
80000   CONTINUE
      
C Perform next step if we haven't reached the maximum
      
      IF ((MAXSTP.EQ.-1).OR.(STEP.LT.MAXSTP)) GOTO 10

90000 CONTINUE

C Finally release memory:

      CALL ZDISP (0,LRESTM,'DRESTM')
      CALL ZDISP (0,LRES,'DRES  ')
      CALL ZDISP (0,LCOORDS,'DCOORDS')

C Set statistical data, finish

      STEPS = STEP

      END
      
************************************************************************
* Optimization algorithm 02: Nelder-Mead
*
* This routine implements the Nelder-Mead search algorithm
*
* It's assumed that there is a calculation routine DCALC which 
* represents a mapping
*     DCALC: R^IDIMS -> R
* i.e. that maps the IDIMS dimensional design variable space R^IDIMS
* to a scalar value. The optimization algorithm tries to find a
* parameter setting x such that
*     DCALC(x) -> min!
* For that purpose DCALC is frequently called and the results are
* then used to find a better parameter setting x.
*
* DSTART specified the starting point where to start the search.
* HSTART is an initial resolution / stepsize parameter, i.e. specifies
* the size of the domain around DSTART where to search for better 
* points. EPS specifies a stopping criterion.
*
* The routine is configured by a "reflection coefficient" ALPHA>0,
* an "expansion coefficient" GAMMA>1 and a "positive contraction
* coefficient" 0<BETA<1.
*
* Nelder Mead stops if > MAXSTP optimization steps have been performed
* or if the last IDIMS found iterates are in a ball with radius EPS.
*
* For a description of the algorithm, see:
* [Bazeraa, M. S.; Sherali, H. D.; Shetty, C. M.; Nonlinear Programming-
*  Theory and Algorithms; 2nd ed., J. Wiley & Sons, New York, 1993;
*  S. 353]
*
* In:
*   IDIMS  - integer
*            Number of dimensions / variables to optimize a value of
*            interest for.
*   DSTART - array [1..IDIMS] of double
*            Start vector where to begin the search.
*            The distance of the starting point to the boundary of the
*            possible parameter settings must be at least HSTART in 
*            order for the algorithm to work correctly!
*   HSTART - array [1..IDIMS] of double
*            For every dimension: initial resolution / stepsize
*   EPS    - double
*            Stopping criterion
*   MAXSTP - Maximum number of optimization steps
*            =-1: not specified.
*   ALPHA  - reflection coefficient; > 0; standard=0.5
*   BETA   - contraction coefficient; > 0, < 1; standard=0.5
*   GAMMA  - expansion coefficient; > 1; standard=2.0
*   IINFO  - array [1..*] of integer
*            Integer information block for the probing routine;
*            is directly passed to DPROBE on call without change.
*   DINFO  - array [1..*] of double
*            Integer information block for the probing routine;
*            is directly passed to DPROBE on call without change.
*   NRES   - integer
*            Size of double-precision array that will be used as
*            temporary space for saving the results of each
*            probing; see DPROBE.
*   DPROBE - SUBROUTINE (IDIMS,LCOORDS,NRES,LRES,IINFO,DINFO)
*            User defined probing function to calculate the value 
*            of interest.
*            In:
*              IDIMS   - Number of dimensions of the design parameters
*              LCOORD  - Handle to array [1..IDIMS] of double
*                        Design parameter setting where to evaluate
*              NRES    - integer
*                        Size of result array
*              LRES    - integer
*                        Handle to array [1..NRES] of double,
*                        where to write the results to
*              PITER   - Number of current probing iterate
*              IINFO,
*              DINFO   - Integer / double precision parameter block
*                        from above; passed to DPROBE unchanged.
*            Out:
*              RES     - array [1..NRES] of double
*                        RES(1) is the value of the point DCOORDS,
*                        which is used by the optimizer to
*                        optimize for. RES(2..NRES) can be set by
*                        DPROBE to arbitrary values that are later
*                        passed to DNXSTP as result values.
*            Remark: The RES-array and the COORD-array are identified by 
*             handles, not by starting addresses. The reason is that
*             the DPROBE-routine is assumed to reassign many handles
*             on the heap, thus destroying a possible starting address
*             of a dynamically allocated array for the results on the
*             heap. By passing the handle to that result array instead
*             of the result array itself, DPROBE will be able to write
*             the results to the correct address.
*     
*            It may happen that the parameter setting is outside
*            of a valid domain. In this case DPROBE must return
*            1D99 to indicate this.
*
*            DPROBE and DNXSTP always belong together. It's acceptable
*            do do memory allocation + calculation in DPROBE, while
*            doing postprocessing and releasing of memory in DNXSTP!
*
*   DNXSTP - SUBROUTINE (IDIMS,DCOORDS,DCOOPT,IINFO,DINFO,
*                        PITER,STEP,IACPT,NRES,LRES,LRSOPT)
*            This subroutine is called after each probing to inform the
*            caller whether a new iterate was accepted or not
*            In:
*              DCOORDS - array [1..IDIMS] of double precision
*                        Parameter set of probed iterate.
*              DCOOPT  - array [1..IDIMS] of double precision
*                        Parameter set of current optimal iterate
*              IINFO,
*              DINFO   - Integer / double precision parameter
*                        block from above
*              PITER   - Number of current probing iterate
*              STEP    - Number of current accepted iterate
*              IACPT   - whether DCOORDS will be accepted as the new
*                        iterate.
*                        =0: DCOORDS will not be accepted
*                        =1: DCOORDS will be accepted
*                        Can be modified by DNXSTP to prevent the
*                        point from being accepted.
*              NRES    - integer
*                        Size of result array RES
*              DRES    - array [1..NRES] of double precision
*                        Result array that was calculated 
*                        in DPROBE for the parameter setting DCOORDS.
*              DRSOPT  - array [1..NRES] of double precision
*                        Result array that was calculated 
*                        in DPROBE for the parameter setting DCOOPT.
*
* Out:
*   DFINAL - array [1..IDIMS] of double
*            Final parameter setting where the search was stopped.
*   FINTOL - double
*            Final step size when the algorithm was stopped
*   STEPS  - integer
*            Number of successful (i.e. accepted because of better 
*            function value) calculations in the algorithm.
************************************************************************

      SUBROUTINE OPTA02 (IDIMS,DSTART,HSTART,EPS,MAXSTP,
     *                   ALPHA, BETA, GAMMA,
     *                   IINFO,DINFO,NRES,DPROBE,DNXSTP,
     *                   DFINAL,FINTOL,STEPS)
     
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
C parameters

      INTEGER IDIMS,MAXSTP,STEPS,IINFO(*),NRES
      DOUBLE PRECISION DSTART(IDIMS),HSTART(IDIMS),DFINAL(IDIMS),EPS
      DOUBLE PRECISION DINFO(*)
      DOUBLE PRECISION FINTOL, ALPHA, BETA, GAMMA
      
      EXTERNAL DPROBE,DNXSTP
      
C     local variables

      INTEGER STEP,CNT,IDX,IDX2,CRSTP
      INTEGER LCOORDS,LRES,CDIM,J,I,IRECLC,ACTDIM
      DOUBLE PRECISION VALMI, VALMX, VAL, DISTMX
      INTEGER IDXMI, IDXMX
      
      INTEGER ARST
      PARAMETER (ARST=4)
      
C     We work with parameters in DFINAL; copy start point from DSTART:

      CALL LCP1 (DSTART,DFINAL,IDIMS)
      
C     Reserve some memory for the current parameter setting as well as
C     for the collection of the points that form our simplex with
C     IDIMS+1 corners.
C     DCOORDS(.,1) is a temporary array for the calculation,
C     DCOORDS(.,2..ARST) are temporary arrays for stored points
C     DCOORDS(.,ARST+1..IDIMS+ARST+1) saves the coordinates of the corners of 
C     the simplex.

      CALL ZNEW (IDIMS*(IDIMS+ARST+1),-1,LCOORDS,'DCOORDS')
      
C     and a little bit of memory for storing the results. 

      CALL ZNEW (NRES*(IDIMS+ARST+1),-1,LRES,'DRES  ')
      
C     Let's go...
      
      STEP = 0
      CNT = 1
      DISTMX = 0D0
      
C     At first we have to build a simplex of IDIMS+1 points which
C     is used for searching for the optimum. This is done by
C     evaluating in DSTART and DSTART+HSTART*e_i, i=1..IDIMS. As DSTART
C     is required to have distance to the boundary of the parameter
C     set of at least HSTART, this gives a nondegenrated, valid
C     simplex.

C     At first calculate the value of the function at DSTART:

      CALL LCP1 (DSTART,DWORK(L(LCOORDS)),IDIMS)
      CALL DPROBE (IDIMS,LCOORDS,NRES,LRES,CNT,IINFO,DINFO)

      I = 1
      CALL DNXSTP (IDIMS,DWORK(L(LCOORDS)),DWORK(L(LCOORDS)),
     *             IINFO,DINFO,CNT,STEP,I,NRES,
     *             DWORK(L(LRES)),DWORK(L(LRES)))
      VALMI = DWORK(L(LRES))
      IDXMI = 0
      STEP = STEP+1
     
C     Save parameter- and result-set as first point

      CALL LCP1 (DWORK(L(LCOORDS)),DWORK(L(LCOORDS)+ARST*IDIMS),IDIMS)
      CALL LCP1 (DWORK(L(LRES)),DWORK(L(LRES)+ARST*NRES),NRES)
      
C     Loop over all possible dimensions:

      ACTDIM = IDIMS

      DO CDIM = 1,IDIMS

        IDX = L(LCOORDS)+(ARST+CDIM)*IDIMS
        
C       Copy DSTART to the probe coordinate array

        CALL LCP1(DSTART,DWORK(IDX),IDIMS)
          
C       Modify the CDIM'th entry according to the stepsize
C       to create a neighbour point of our current midpoint:

        DWORK(IDX+CDIM-1) = DWORK(IDX+CDIM-1) + HSTART (CDIM)
           
C       Calculate the actual number of dimensions.
C       If the initial step length of a dimension is =0, it does
C       not count!

        IF (HSTART(CDIM).EQ.0D0) ACTDIM=ACTDIM-1
           
      END DO

C     Now we have prepared the following simplex:
C
C        DSTART+H*e_2
C            |       \_
C            |         \
C         DSTART--------DSTART+H*e_1
C
C     In the first step of the optimization loop that follows now, we have
C     to calculate the function values of all these points. This is done
C     by setting IRECLC<>-1. More precisely we set IRECLC to 0 which
C     not recalculates all values, but prevents coordinate 0 from being
C     recalculated, too - since it has been calculated, and we don't have
C     to repeat that...

      IRECLC = 0

C     We can now start our optimization loop.
C
C     Crippled Do-While-Loop

      CRSTP = 0
100   CONTINUE

        CRSTP = CRSTP+1
      
C       If IRECLC <> -1, this indicates that the values of all points 
C       of the simplex except for IRECLC have to be recalculated:
      
        IF (IRECLC.NE.-1) THEN
      
C         For the recalculation loop over all points in the simplex:
      
          DO CDIM = 0,IDIMS
        
C           Prevent point IRECLC from being recalculated

            IF (CDIM.NE.IRECLC) THEN
        
C             Increase probing counter

              CNT = CNT + 1
        
C             Copy the coordinate of the point to the array that is used
C             in the call to DPROBE - i.e. the beginning of DCOORDS:
              
              IDX = L(LCOORDS)+(ARST+CDIM)*IDIMS
              CALL LCP1(DWORK(IDX),DWORK(L(LCOORDS)),IDIMS)
          
C             If the step length is =0 in any direction, we have to ignore
C             that dimension. Set the value to 1D99 there and do nothing.

              IF ((CDIM.NE.0).AND.(HSTART (CDIM).EQ.0D0)) THEN
              
                DWORK(L(LRES)) = 1D99
              
              ELSE
          
C               Calculate the value of interest:

                CALL DPROBE (IDIMS,LCOORDS,NRES,LRES,CNT,IINFO,DINFO)
     
C               Try to accept the point as minimum:
     
                IF (DWORK(L(LRES)).LT.VALMI) THEN     
                  I = 1
                ELSE
                  I = 0
                END IF
        
                IDX = L(LCOORDS)+(ARST+IDXMI)*IDIMS
                IDX2 = L(LRES)+(ARST+IDXMI)*NRES
                CALL DNXSTP (IDIMS,DWORK(L(LCOORDS)),DWORK(IDX),
     *                       IINFO,DINFO,CNT,STEP,I,NRES,
     *                       DWORK(L(LRES)),DWORK(IDX2))
     
C               Remember our current minimum
     
                IF (I.EQ.1) THEN
                  IDXMI = CDIM
                  VALMI = DWORK(L(LRES))
                  STEP = STEP + 1
                END IF
       
C               Anyway save the results in the result-array

                CALL LCP1 (DWORK(L(LRES)),
     *                     DWORK(L(LRES)+(ARST+CDIM)*NRES),NRES)
     
              END IF ! (dim<>0 and hstart(dim)=0)
          
            END IF ! CDIM<>IRECLC
            
          END DO

C         Set IRECLC to -1 to mark the simplex as calculated:

          IRECLC = -1
        
        END IF

C       At first we check if all points are in a ball with radius EPS
C       around the midpoint of the simplex. For that purpose we have to
C       build the midpoint of the current simplex:

        IDX = L(LCOORDS) + IDIMS
        CALL LCL1(DWORK(IDX),IDIMS)
        DO J = 0,IDIMS
          IDX2 = L(LCOORDS)+(ARST+J)*IDIMS
          CALL LLC1(DWORK(IDX2),DWORK(IDX),IDIMS,1D0,1D0)
        END DO
        CALL LSC1 (DWORK(IDX),IDIMS,1D0/DBLE(IDIMS+1))

C       Check all the distances:

        DISTMX = 0D0
        DO J=0,IDIMS
          VAL = 0D0
          IDX2 = L(LCOORDS)+(ARST+J)*IDIMS
          DO I=1,IDIMS
            VAL = VAL + (DWORK(IDX+I-1)-DWORK(IDX2+I-1))**2
          END DO
          VAL = SQRT(VAL)
          IF (VAL.GT.DISTMX) DISTMX = VAL
        END DO
        
C       Cancal the loop if the simplex is small enough.

        IF (DISTMX.LE.EPS) GOTO 80000

C       At least one point is more far away from the midpoint than EPS.
C       So we can proceed with the real Nelder-Mead step.

C       Calculate maximum function value from the calculated points. 
C       The result can always be found in the first entry of each
C       RES-array. Ignore those dimensions which have initially a 
C       step size parameter =0:

        IDXMX = 0
        VALMX = VALMI
        
        DO I=0,IDIMS
          IF ( ((I.EQ.0).OR.((I.NE.0).AND.(HSTART(I).NE.0D0)) ).AND.
     *         (DWORK(L(LRES)+(ARST+I)*NRES).GT.VALMX)) THEN
            VALMX = DWORK(L(LRES)+(ARST+I)*NRES)
            IDXMX = I
          END IF
        END DO

C       Create the midpoint of the simplex that arises when not respecting
C       the point with the maximum value. Save the coordinates in 
C       DCOORDS(.,2).

        IDX = L(LCOORDS) + IDIMS
        CALL LCL1(DWORK(IDX),IDIMS)
        DO J = 0,IDIMS
          IF ((J.NE.IDXMX).AND.
     *        ((J.EQ.0).OR.((J.NE.0).AND.(HSTART(J).NE.0D0)))) THEN
            IDX2 = L(LCOORDS)+(ARST+J)*IDIMS
            CALL LLC1(DWORK(IDX2),DWORK(IDX),IDIMS,1D0,1D0)
          END IF
        END DO
        CALL LSC1 (DWORK(IDX),IDIMS,1D0/DBLE(ACTDIM))

C       Create x^ = x_mid + ALPHA*(xmid-xmax) and probe that point:
C
C                O
C              / |
C            /   |
C          /     |
C        Xmax---Xmid---> x^
C          \     |
C            \   |
C              \ |
C               Xmin

        IDX = L(LCOORDS)+IDIMS
        IDX2 = L(LCOORDS)+(ARST+IDXMX)*IDIMS
        CALL LCP1(DWORK(IDX2),DWORK(L(LCOORDS)),IDIMS)
        CALL LLC1(DWORK(IDX),DWORK(L(LCOORDS)),IDIMS,(1D0+ALPHA),-ALPHA)

        CNT = CNT + 1

        CALL DPROBE (IDIMS,LCOORDS,NRES,LRES,CNT,IINFO,DINFO)
     
C       Probably accept it as new minimum

        IF (DWORK(L(LRES)).LT.VALMI) THEN     
          I = 1
        ELSE
          I = 0
        END IF
        
        IDX = L(LCOORDS)+(ARST+IDXMI)*IDIMS
        IDX2 = L(LRES)+(ARST+IDXMI)*NRES
        CALL DNXSTP (IDIMS,DWORK(L(LCOORDS)),DWORK(IDX),
     *               IINFO,DINFO,CNT,STEP,I,NRES,
     *               DWORK(L(LRES)),DWORK(IDX2))

C       Don't update IDXMI; for the moment if the new point is a 
C       minimum, it's a temporary minimum!

        IF (I.EQ.1) THEN
        
C         The new point is a better minimum, thus replace our
C         maximum by that:

          CALL LCP1 (DWORK(L(LCOORDS)),
     *               DWORK(L(LCOORDS)+(ARST+IDXMX)*IDIMS),IDIMS)
          CALL LCP1 (DWORK(L(LRES)),
     *               DWORK(L(LRES)+(ARST+IDXMX)*NRES),NRES)

          IDXMI = IDXMX
          VALMI = DWORK(L(LRES))
          STEP = STEP + 1
          
C         But we don't give up here! Try to create an even better
C         point! Build a second probing point and overwrite our old:
C
C           x_e = xmid + GAMMA(x^ - xmid)

          IDX = L(LCOORDS)+IDIMS
          IDX2 = L(LCOORDS)
          CALL LLC1(DWORK(IDX),DWORK(IDX2),IDIMS,1D0-GAMMA,GAMMA)
          
C         Probe it
          
          CNT = CNT + 1

          CALL DPROBE (IDIMS,LCOORDS,NRES,LRES,CNT,IINFO,DINFO)
     
C         Probably accept it as new minimum

          IF (DWORK(L(LRES)).LT.VALMI) THEN     
            I = 1
          ELSE
            I = 0
          END IF
        
          IDX = L(LCOORDS)+(ARST+IDXMI)*IDIMS
          IDX2 = L(LRES)+(ARST+IDXMI)*NRES
          CALL DNXSTP (IDIMS,DWORK(L(LCOORDS)),DWORK(IDX),
     *                 IINFO,DINFO,CNT,STEP,I,NRES,
     *                 DWORK(L(LRES)),DWORK(IDX2))
     
          IF (I.EQ.1) THEN
          
C           Wow, we got even better; take that point to replace our
C           former maximum!
          
            CALL LCP1 (DWORK(L(LCOORDS)),
     *                 DWORK(L(LCOORDS)+(ARST+IDXMX)*IDIMS),IDIMS)
            CALL LCP1 (DWORK(L(LRES)),
     *                 DWORK(L(LRES)+(ARST+IDXMX)*NRES),NRES)
          
            IDXMI = IDXMX
            VALMI = DWORK(L(LRES))
            STEP = STEP + 1
          
          END IF
        
        ELSE
        
C         Ok, we haven't got better directly. Nelder Mead now replaces
C         the maximum if the value is less than the second but largest 
C         value! So at first determine the second but largest value.
C         For dimensions with initial stepsize =0, the value is 1D99,
C         so they are for sure not respected here.

          VAL = VALMI
          J = IDXMI
          DO I=0,IDIMS
            IF ((I.NE.IDXMX).AND.
     *          (DWORK(L(LRES)+(ARST+I)*NRES).GT.VAL)) THEN
              J = I
              VAL = DWORK(L(LRES)+(ARST+I)*NRES)
            END IF
          END DO
          
          IF (DWORK(L(LRES)).LT.VAL) THEN
          
C           The new probed point was better than the second largest
C           value. In this case replace our maximum by the new point.
C           The minimum has not changed.
          
            CALL LCP1 (DWORK(L(LCOORDS)),
     *                 DWORK(L(LCOORDS)+(ARST+IDXMX)*IDIMS),IDIMS)
            CALL LCP1 (DWORK(L(LRES)),
     *                 DWORK(L(LRES)+(ARST+IDXMX)*NRES),NRES)
          
          ELSE
          
C           In this case the probed point is either larger or smaller
C           than our current maximum, but always larger than the second
C           largest value:
C
C             f(xi) < ... < f(x2ndmax) < ( f(x^2) <> f(xmax) )
C
C           Let's find out which one is smaller, and copy that to 
C           DCOORDS(.,3). We will denote that by x':

            IF (DWORK(L(LRES)).LT.VALMX) THEN
              CALL LCP1 (DWORK(L(LCOORDS)),
     *                   DWORK(L(LCOORDS)+2*IDIMS),IDIMS)
              CALL LCP1 (DWORK(L(LRES)),
     *                   DWORK(L(LRES)+2*NRES),NRES)
            ELSE
              CALL LCP1 (DWORK(L(LCOORDS)+(ARST+IDXMX)*IDIMS),
     *                   DWORK(L(LCOORDS)+2*IDIMS),IDIMS)
              CALL LCP1 (DWORK(L(LRES)+(ARST+IDXMX)*NRES),
     *                   DWORK(L(LRES)+2*NRES),NRES)
            END IF

C           Define a new probing point:
C
C             x'' = xmid + BETA(x'-xmid):

            CALL LCP1(DWORK(L(LCOORDS)+IDIMS),
     *                DWORK(L(LCOORDS)),IDIMS)
            CALL LLC1(DWORK(L(LCOORDS)+2*IDIMS),
     *                DWORK(L(LCOORDS)),
     *                IDIMS,BETA,1D0-BETA)
     
C           Probe that point:

            CNT = CNT + 1

            CALL DPROBE (IDIMS,LCOORDS,NRES,LRES,CNT,IINFO,DINFO)
     
C           Probably accept it as new minimum

            IF (DWORK(L(LRES)).LT.VALMI) THEN     
              I = 1
            ELSE
              I = 0
            END IF

            IDX = L(LCOORDS)+(ARST+IDXMI)*IDIMS
            IDX2 = L(LRES)+(ARST+IDXMI)*NRES
            CALL DNXSTP (IDIMS,DWORK(L(LCOORDS)),DWORK(IDX),
     *                   IINFO,DINFO,CNT,STEP,I,NRES,
     *                   DWORK(L(LRES)),DWORK(IDX2))
     
C           Check whether the new point has better or worse
C           value than x':

            IF (DWORK(L(LRES)).LT.VAL) THEN
            
C             Replace our maximum value by that point:
            
              CALL LCP1 (DWORK(L(LCOORDS)),
     *                   DWORK(L(LCOORDS)+(ARST+IDXMX)*IDIMS),IDIMS)
              CALL LCP1 (DWORK(L(LRES)),
     *                   DWORK(L(LRES)+(ARST+IDXMX)*NRES),NRES)
            
C             Have we even to accept it as new minimum?

              IF (I.EQ.1) THEN
                VALMI = DWORK(L(LRES))
                IDXMI = IDXMX
                STEP = STEP + 1
              END IF
            
            ELSE 
            
C             x'' has worse value than x', so surely worse than our
C             minimum. We have no chance probing around here, so we
C             move our whole simplex to a new part of the parameter
C             space. We move all coordinates:
C
C                 xj = xj + 1/2 * (xmin - xj)
C
C             This will not move xmin of course. Afterwards we must
C             recalculate everything except for xmin.

              DO I=0,IDIMS
                IF (I.NE.IDXMI) THEN
                  CALL LLC1 (DWORK(L(LCOORDS)+(ARST+IDXMI)*IDIMS),
     *                       DWORK(L(LCOORDS)+(ARST+I)*IDIMS),
     *                       IDIMS,0.5D0,0.5D0)
                END IF
              END DO
            
C             Recalculate everything in the next sweep except for 
C             the current minimum - which coordinates have not changed:

              IRECLC = IDXMI
            
            END IF
          
          END IF
        
        END IF
      
C     Perform next sweep

      IF ((MAXSTP.LE.-1).OR.(CRSTP.LT.MAXSTP)) GOTO 100

C     -----
      
80000 CONTINUE

C     At this point Nelder Mead is finished. Either we have reached 
C     the maximum number of steps of the stopping criterion.
C
C     In IDXMI we have the index of the best point so far. Copy this
C     to the result array:

      CALL LCP1 (DWORK(L(LCOORDS)+(ARST+IDXMI)*IDIMS),
     *           DFINAL,IDIMS)
     
90000 CONTINUE

C Finally release memory:

      CALL ZDISP (0,LRES,'DRES  ')
      CALL ZDISP (0,LCOORDS,'DCOORDS')

C Set statistical data, finish.
C The maximum distance was calculated before every loop:
     
      FINTOL = DISTMX
      STEPS = STEP-1

      END
      
************************************************************************
* Optimization algorithm 03: Gradient based line-search
*
* This routine implements basic line search with gradient
* information, reconstructed with a simple finite-difference approach,
* Armijo step length principle.
*
* It's assumed that there is a calculation routine DCALC which 
* represents a mapping
*     DCALC: R^IDIMS -> R
* i.e. that maps the IDIMS dimensional design variable space R^IDIMS
* to a scalar value. The optimization algorithm tries to find a
* parameter setting x such that
*     DCALC(x) -> min!
* For that purpose DCALC is frequently called and the results are
* then used to find a better parameter setting x.
*
* DSTART specified the starting point where to start the search.
* HSTART is an initial resolution parameter array, which specifies
* the initial steplength in each direction that is used for the
* (re)construction of the gradient information.
* EPS specifies a stopping criterion.
*
* Line-search stops the search if the stepsize of the line-search
* drops belop EPS, or if > MAXSTP optimization steps 
* have been performed.
*
* For a description of the algorithm, see:
* [Großmann, Ch.; Terno, J.; Numerik der Optimierung; Teubner 
*  Studienführer Mathematik, pp. 66ff]
*
* In:
*   IDIMS  - integer
*            Number of dimensions / variables to optimize a value of
*            interest for.
*   DSTART - array [1..IDIMS] of double
*            Start vector where to begin the search.
*   HSTART - array [1..IDIMS] of double
*            For every dimension: initial resolution / stepsize
*   DELTA  - double
*            Control parameter of Armijo step length control.
*            0<DELTA<1; standard)=0.5D0
*   EPS    - double
*            Stopping criterion. Algorithm stops if the stepsize
*            of the line-length drops below EPS.
*            (which arise from modifying HSTART) drop below EPS.
*   MAXSTP - Maximum number of optimization steps
*            =-1: not specified.
*   IINFO  - array [1..*] of integer
*            Integer information block for the probing routine;
*            is directly passed to DPROBE on call without change.
*   DINFO  - array [1..*] of double
*            Integer information block for the probing routine;
*            is directly passed to DPROBE on call without change.
*   NRES   - integer
*            Size of double-precision array that will be used as
*            temporary space for saving the results of each
*            probing; see DPROBE.
*   DPROBE - SUBROUTINE (IDIMS,LCOORDS,NRES,LRES,IINFO,DINFO)
*            User defined probing function to calculate the value 
*            of interest.
*            In:
*              IDIMS   - Number of dimensions of the design parameters
*              LCOORD  - Handle to array [1..IDIMS] of double
*                        Design parameter setting where to evaluate
*              NRES    - integer
*                        Size of result array
*              LRES    - integer
*                        Handle to array [1..NRES] of double,
*                        where to write the results to
*              PITER   - Number of current probing iterate
*              IINFO,
*              DINFO   - Integer / double precision parameter block
*                        from above; passed to DPROBE unchanged.
*            Out:
*              RES     - array [1..NRES] of double
*                        RES(1) is the value of the point DCOORDS,
*                        which is used by the optimizer to
*                        optimize for. RES(2..NRES) can be set by
*                        DPROBE to arbitrary values that are later
*                        passed to DNXSTP as result values.
*            Remark: The RES-array and the COORD-array are identified by 
*             handles, not by starting addresses. The reason is that
*             the DPROBE-routine is assumed to reassign many handles
*             on the heap, thus destroying a possible starting address
*             of a dynamically allocated array for the results on the
*             heap. By passing the handle to that result array instead
*             of the result array itself, DPROBE will be able to write
*             the results to the correct address.
*     
*            It may happen that the parameter setting is outside
*            of a valid domain. In this case DPROBE must return
*            1D99 to indicate this.
*
*            DPROBE and DNXSTP always belong together. It's acceptable
*            do do memory allocation + calculation in DPROBE, while
*            doing postprocessing and releasing of memory in DNXSTP!
*
*   DNXSTP - SUBROUTINE (IDIMS,DCOORDS,DCOOPT,IINFO,DINFO,
*                        PITER,STEP,IACPT,NRES,DRES,DRSOPT)
*            This subroutine is called after each probing to inform the
*            caller whether a new iterate was accepted or not
*            In:
*              DCOORDS - array [1..IDIMS] of double precision
*                        Parameter set of probed iterate.
*              DCOOPT  - array [1..IDIMS] of double precision
*                        Parameter set of current optimal iterate
*              IINFO,
*              DINFO   - Integer / double precision parameter
*                        block from above
*              PITER   - Number of current probing iterate
*              STEP    - Number of current accepted iterate
*              IACPT   - whether DCOORDS will be accepted as the new
*                        iterate.
*                        =0: DCOORDS will not be accepted
*                        =1: DCOORDS will be accepted
*                        Can be modified by DNXSTP to prevent the
*                        point from being accepted.
*              NRES    - integer
*                        Size of result array RES
*              DRES    - array [1..NRES] of double precision
*                        Result array that was calculated 
*                        in DPROBE for the parameter setting DCOORDS.
*              DRSOPT  - array [1..NRES] of double precision
*                        Result array that was calculated 
*                        in DPROBE for the parameter setting DCOOPT.
*
* Out:
*   DFINAL - array [1..IDIMS] of double
*            Final parameter setting where the search was stopped.
*   FINTOL - array [1..IDIMS] of double
*            Final step size in each direction, used for the gradient
*            calculation when the algorithm was stopped
*   STEPS  - integer
*            Number of steps the algorithm used.
************************************************************************

      SUBROUTINE OPTA03 (IDIMS,DSTART,HSTART,DELTA,EPS,MAXSTP,
     *                   IINFO,DINFO,NRES,DPROBE,DNXSTP,
     *                   DFINAL,FINTOL,STEPS)
     
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
C parameters

      INTEGER IDIMS,MAXSTP,STEPS,IINFO(*),NRES
      DOUBLE PRECISION DSTART(IDIMS),HSTART(IDIMS),DFINAL(IDIMS),EPS
      DOUBLE PRECISION DINFO(*),DELTA
      DOUBLE PRECISION FINTOL(IDIMS)
      
      EXTERNAL DPROBE,DNXSTP
      
C     local variables

      INTEGER STEP,CNT
      INTEGER LCOORDS,LCOORFI,LRES,LRESTM,CDIM,J,I,K
      DOUBLE PRECISION NEWVAL,VAL,H,VAL1,VALM,VAL2,GR,ALPHA,MINSL,FT
      
C     During the algorithm:
C
C      DFINAL  - saves the current best parameter setting
C      LCOORDS - temporary array that stores the current parameter
C                setting for call to DPROBE
C      LRES    - Handle to:
C                array [1..NRES] of double
C                --> temporary array for saving the result of the 
C                    call to DPROBE
C                array [NRES+1..2*NRES] of doubole
C                --> saves the best result of any call to DPROBE,
C                    corresponding to the parameter set DFINAL.
C      LCOORFI - Handle to:
C                array [1..IDIMS,1] of double
C                --> parameter set of the midpoint of the FD-stencil
C                array [1..IDIMS,2..2*IDIMS+1] of double
C                --> parameter set of the corners of the FD-stencil
C                array [1..IDIMS,2*IDIMS+2] of double
C                --> gradient information
C
C     We work with parameters in DFINAL; copy start point from DSTART:

      CALL LCP1 (DSTART,DFINAL,IDIMS)
      
C     Current stepsize in every dimension

      CALL LCP1(HSTART,FINTOL,IDIMS)
      
C     Reserve memory for current parameter setting:

      CALL ZNEW (IDIMS,-1,LCOORDS,'DCOORDS')
      
C     Reserve some memory for the parameter settings in all directions.
C     We make a typical 2nd order finite difference approach, which results
C     in the typical "stencil", e.g. 4-point stencil in 2D, 9-point-stencil
C     in 3D,...
C        -1
C     -1  4 -1
C        -1
C     1 midpoint + 2*IDIMS corners + gradient information:

      CALL ZNEW ((2+2*IDIMS)*IDIMS,-1,LCOORFI,'DCOORDS')
      
C     Memory for storing the results:

      CALL ZNEW (2*NRES,-1,LRES,'DRES  ')
      CALL ZNEW ((1+2*IDIMS)*NRES,-1,LRESTM,'DRESTM')

C     Copy the start point to the array containing our current
C     midpoint of the stencil:
      
      CALL LCP1(DFINAL,DWORK(L(LCOORFI)),IDIMS)
      
C     Calculate the value of the functional in the midpoint of
C     the stencil:

      CNT = 1
      STEP = 0

      CALL LCP1(DWORK(L(LCOORFI)),DWORK(L(LCOORDS)),IDIMS)
      CALL DPROBE (IDIMS,LCOORDS,NRES,LRES,CNT,IINFO,DINFO)
      CALL LCP1(DWORK(L(LRES)),DWORK(L(LRESTM)),NRES)

C     Value of interest in first variable of RES-array:
    
      VAL = DWORK(L(LRES))
    
C     The very first calculation is always accepted as 
C     first iterate:

      I = 1
      CALL DNXSTP (IDIMS,DWORK(L(LCOORDS)),DFINAL,
     *                 IINFO,DINFO,CNT,STEP,I,NRES,
     *                 DWORK(L(LRES)),DWORK(L(LRES)))
    
C     Remember parameter setting and 
C     value in DFINAL and the second part of DRES:
          
      CALL LCP1(DWORK(L(LCOORDS)),DFINAL,IDIMS)
      CALL LCP1(DWORK(L(LRES)),DWORK(L(LRES)+NRES),NRES)

C     Let's go... 
      
      STEP = STEP + 1
      
C     Crippled DO-WHILE loop to realize a maximum of MAXSTP steps:
      
10    CONTINUE

C       Loop over all possible dimensions:

        DO CDIM = 0,IDIMS-1
        
C         In every dimension: 1x forward, 1x backward
        
          DO J=0,1
          
C           Increase probing counter

            CNT = CNT + 1
        
C           Copy our current point to the probe coordinate array

            CALL LCP1(DWORK(L(LCOORFI)),DWORK(L(LCOORDS)),IDIMS)
          
C           Modify the CDIM'th entry according to the stepsize
C           to create a neighbour point of our current midpoint:

            H = FINTOL(CDIM+1)
            DWORK(L(LCOORDS)+CDIM) = DWORK(L(LCOORDS)+CDIM) 
     *                                 + H - DBLE(J)*2D0*H
           
C           Calculate the value of interest.
C           If H=0, save 0.0 in this dimension as value.

            IF (H.NE.0D0) THEN
              CALL DPROBE (IDIMS,LCOORDS,NRES,LRES,CNT,
     *                     IINFO,DINFO)
            ELSE
              DWORK(L(LRES)) = 0D0
            END IF
     
C           Save the result for later gradient calculation.
C           The first NRES entries belong to the midpoint, the
C           other entries pairwise for the "left" and "right"
C           neighbouring corner of the stencil.

            I = NRES + ( 2*CDIM + J ) * NRES
            CALL LCP1 (DWORK(L(LRES)),DWORK(L(LRESTM)+I),NRES)
     
C           Also save the coordinate of that point:

            I = IDIMS + ( 2*CDIM + J ) * IDIMS
            CALL LCP1(DWORK(L(LCOORDS)),DWORK(L(LCOORFI)+I),IDIMS)
          
C           Don't do anything if H=0, since the solver was not
C           called in this case.

            IF (H.NE.0) THEN          
            
C             Result in first variable of 2nd part of RES-array:
     
              NEWVAL = DWORK(L(LRES))
          
C             Is it better? -> Change the iterate...

              I = 0
              IF (NEWVAL.LT.VAL) I = 1
          
C             Indicate the framework whether we will accept the iterate
C             or not. 

              CALL DNXSTP (IDIMS,DWORK(L(LCOORDS)),DFINAL,
     *                     IINFO,DINFO,CNT,STEP,I,NRES,
     *                     DWORK(L(LRES)),DWORK(L(LRES)+NRES))
            ELSE
              I = 0
            END IF
     
C           Are we allowed to accept it? If yes, do that!
      
            IF (I.NE.0) THEN

C             Accept the new iterate: Remember parameter setting and 
C             value
            
              CALL LCP1(DWORK(L(LCOORDS)),DFINAL,IDIMS)
              CALL LCP1(DWORK(L(LRES)),DWORK(L(LRES)+NRES),NRES)
              VAL = NEWVAL
              STEP = STEP + 1
            
            END IF
        
C         Otherwise check next dimension...

          END DO
        
        END DO
      
C       All dimensions probed. Now we try to calculate an approximate
C       gradient from the gained information.
      
        DO CDIM=0,IDIMS-1
        
C         If both information tags are available,
C         use 2nd order approximation. 

          I = NRES + 2*CDIM * NRES
          VALM = DWORK(L(LRESTM))
          VAL1 = DWORK(L(LRESTM)+I)
          VAL2 = DWORK(L(LRESTM)+I+NRES)
          
          FT = FINTOL(CDIM+1)

C         A step length of 0.0 in a dimension means VAL1=VAL2=0.
C         Set the denominator of the finite-difference quotient
C         to 0 to get a valid FD-fraction 0 in the later calculation.

          IF (FT.EQ.0D0) FT = 1D0
          
          IF ((VAL1.NE.1D99).AND.(VAL2.NE.1D99)) THEN
          
C           Calculate the gradient and save it in DCOORFI
          
            GR = (VAL1-VAL2) / (2D0*FT)

          ELSE IF (VAL1.NE.1D99) THEN

C           If only one information tag is available, use 1st order
C           approximation.

            GR = (VAL1-VALM) / FT

          ELSE IF (VAL2.NE.1D99) THEN

            GR = (VALM-VAL2) / FT
          
          ELSE
 
C           If nothing is available... ehm...
C           stop calculation, halfen the step length and try again!

            GR = 1D99
 
          END IF

          IF (GR.NE.1D99) THEN
          
C           Save the *negative* gradient!
          
            I = IDIMS + 2*IDIMS*IDIMS
            DWORK(L(LCOORFI)+ I + CDIM) = -GR
            
          ELSE
          
C           Halfen all step lengths
          
            H = 0.5D0 * FINTOL(1)
            DO J=1,IDIMS
              FINTOL(J) = 0.5D0 * FINTOL(J)
              IF (FINTOL(J).GT.H) H = FINTOL(J)
            END DO
            
C           And repeat the calculation of the corners of the stencil,
C           if the step length is not too small.
C           H=0 is only possible if it was =0 initially, which means
C           that this dimension should be ignored!
            
            IF ((H.NE.0D0).AND.(H.GE.EPS)) THEN
              GOTO 10
            ELSE
              GOTO 90000
            END IF
            
          END IF
        
        END DO
      
C       If we are here, the gradient information was successfully
C       calculated. 
C       The calculated gradient (or the negative of it, respectively)
C       describes our current line where to search for a minimum.
C       We now use the Armijo-principle to determine a step-length,
C       i.e. we search for a step-length ALPHA_k such that:
C
C       ALPHA = MAX (a : f(x_i + a*d_i) <=
C                        a*DELTA*grad(f(x_i))*d_i)
C
C       with d_k = search direction, a = 2^0, 2^-1, 2^-2, 2^-3,...
C
C       So at first calculate that DELTA*grad(f(x_i))*d_i =: MINSL

        MINSL = 0D0
        I = IDIMS + 2*IDIMS*IDIMS
        DO J=0,IDIMS-1
          MINSL = -DELTA*DWORK(L(LCOORFI)+I+J)**2
        END DO

C       Crippled DO-WHILE-Loop 

        ALPHA = 1D0
        
100     CONTINUE

C         Calculate the current step-length ( ||ALPHA*d_i|| ) to check,
C         if it drops below EPS. If yes, stop the algorithm.

          I = IDIMS + 2*IDIMS*IDIMS
          CALL LL21 (DWORK(L(LCOORFI)+I),IDIMS,H)
C          print *,'*** Current step-length: ',ALPHA*H
          IF (ALPHA*H.LT.EPS) GOTO 90000

C         Calculate coordinates of new point: x_i+ALPHA*d_i
C         with: d_i = -gradf(x_i)

          CALL LCP1(DWORK(L(LCOORFI)),DWORK(L(LCOORDS)),IDIMS)
          I = IDIMS + 2*IDIMS*IDIMS
          DO K=0,IDIMS-1
            DWORK(L(LCOORDS)+K) = DWORK(L(LCOORDS)+K) + 
     *                            ALPHA*DWORK(L(LCOORFI)+I+K)
          END DO

C         Increase probing counter

          CNT = CNT + 1

C         Probe the function at that point:

          CALL DPROBE (IDIMS,LCOORDS,NRES,LRES,CNT,
     *                 IINFO,DINFO)
     
C         Result in first variable of 2nd part of RES-array:
     
          NEWVAL = DWORK(L(LRES))
          
C         Is it better? -> Change the iterate...

          I = 0
          IF (NEWVAL.LT.VAL) I = 1
          
C         Indicate the framework whether we will accept the iterate
C         or not

          CALL DNXSTP (IDIMS,DWORK(L(LCOORDS)),DFINAL,
     *                   IINFO,DINFO,CNT,STEP,I,NRES,
     *                   DWORK(L(LRES)),DWORK(L(LRES)+NRES))
     
C         Are we allowed to accept it? If yes, do that!
      
          IF (I.NE.0) THEN

C           Accept the new iterate: Remember parameter setting and 
C           value
            
            CALL LCP1(DWORK(L(LCOORDS)),DFINAL,IDIMS)
            CALL LCP1(DWORK(L(LRES)),DWORK(L(LRES)+NRES),NRES)
            VAL = NEWVAL
            STEP = STEP + 1
            
          END IF
          
C         Check if the new function value is 
C         <= ALPHA*DELTA*grad(f(x_i))*d_i). If yes, we found
C         our ALPHA and the point that lets us decrease:

          IF (NEWVAL.LE. DWORK(L(LRESTM))+ALPHA*MINSL ) THEN

C           Juhey, we can continue searching at that point.
C           Take that iterate as new center:

            CALL LCP1(DWORK(L(LCOORDS)),DWORK(L(LCOORFI)),IDIMS)
            CALL LCP1(DWORK(L(LRES)),DWORK(L(LRESTM)),NRES)
        
C           Decrease the step-length for the gradient calculation,
C           if our search radius is too large, set the appropriate
C           coordinate step length to the half ALPHA*H.

C            I = IDIMS + 2*IDIMS*IDIMS
C            DO J=1,IDIMS
C              H=DWORK(L(LCOORFI)+I+J)
C              IF (FINTOL(J).GE.0.5D0*ALPHA*H) FINTOL(J) = 0.5D0*ALPHA*H
C            END DO
        
C           Perform next step if we haven't reached the maximum
      
            IF ((MAXSTP.EQ.-1).OR.(STEP.LT.MAXSTP)) GOTO 10

          ELSE

C           No that's not good enough. Halfen out step length ALPHA

            ALPHA = 0.5D0*ALPHA

            GOTO 100
      
          END IF
          
C     If we reach this point, ALPHA dropped below EPS or we reached
C     the maximum number of iterations - we can stop here. The best
C     point is to be found in DFINAL.

90000 CONTINUE

C Finally release memory:

      CALL ZDISP (0,LRESTM,'DRESTM')
      CALL ZDISP (0,LRES,'DRES  ')
      CALL ZDISP (0,LCOORFI,'DCOORDS')
      CALL ZDISP (0,LCOORDS,'DCOORDS')

C Set statistical data, finish

      STEPS = STEP
      
C      PRINT *,(FINTOL(I),I=1,IDIMS)

      END
      
