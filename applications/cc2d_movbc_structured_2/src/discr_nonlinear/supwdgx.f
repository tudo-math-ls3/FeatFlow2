************************************************************************
* This file implements the assembling of the nonlinear part of the
* system matrix with Streamline Diffusion stabilization technique
* in extended calling convention.
************************************************************************

************************************************************************
* Streamline diffusion
*
* Adds the SUPG-part on matrix block A after it was initialized by the
* linear part or builds up the complete nonlinear system matrix.
* The input vector Ui is the old velocity field.
* The input vectors UjLi are the transport directions.
*
* Parametric version
*
* Extended calling convention
*
* In:
*   A, NA,
*   KCOLA,
*   LLDA   - array [1..*] of double/integer
*            Structure arrays of system matrix, maybe initialised
*            with linear parts of the nonlinear system matrix.
*   KVERT,
*   KMID,
*   DCORVG - usual geometry information; must correspond to TRIA
*   
*   TRIA   - array [1..SZTRIA] of integer
*            Triangulation structure of the underlying mesh
*
*   NEQ    - length of solution/defect vectors
*
*   U1L1,
*   U1L2   - array [1..NEQ] of double
*            main velocity field used for assembling of 
*            the nonlinearity. Can be undefined if A1L=0.
*   U2L1,
*   U2L2   - array [1..NEQ] of double
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
*                A can be undefined (but not KCOLA,KLDA!).
*                Can be used to modify the defect vector matrix
*                free. Note that the rules how to set up the
*                system matrix (NY, DCMASS,...) still hold,
*                although the matrix itself is not modified!
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
*   BFULMT - Whether or not to generate the full nonlinear system
*            matrix.
*            =false: generate only the nonlinear part 
*                         THWEIG * u grad(.)
*                    of the system matrix and add it to the existing
*                    matrix in A
*            =true : Generate the full system matrix
*               [  M*I  +  THWEIG * (-nu * Laplace(.) * u * grad(.)) ] 
*                    and add it to A
*            
*   THSTEP - weighting factor of the nonlinearity N(u) which is to be
*            added to the matrix
*            (is used as step-size of a Theta-scheme)
*
*   DCMASS - This defines how to modify the system matrix - if the
*            matrix is modified at all (see IDEF!).
*             =0: subtract the mass matrix from the system matrix;
*                 The nonlinearity is not build.
*            <>0: Add the nonlinear term including the stabilization
*                 to the system matrix.
*                 If the full matrix with not-lumped mass matrix 
*                 is to be build (IPRECA=4 and IMASS=1),
*                 add DCMASS*(Mass matrix) to the system matrix.
*
*   UPSAM  - control parameter of the streamline diffusion
*   RE     - 1/nu = viscosity
*   ISTOK  - =1 if the Stokes equation is to be discretized
*            rather than Navier-Stokes.
*   IMASS  - if BFULMT=TRUE:
*            =1: add the real mass matrix to the system matrix
*            =0: don't add the mass matrix to the system matrix
*
*   ELE    - used element for the discretisation
*   ICUBN  - Cubature formula to use for the discretization
*
* Out:
*   A      - system matrix; only if IDEF=0,1.
*            The nonlinearity is added to that matrix or
*            the linearized nonlinear system matrix is build
*            completely, respectively.
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
*   If BFULMT=true, IMASS=1, DCMASS<>0, D is updated according to
*
*             D = D - DCMASS*M*UUx - THSTEP * UUx * grad (Ux)
*   
*   If BFULMT=true, IMASS=1, DCMASS=0, D is updated according to
*
*             D = D + M*UUx - THSTEP * UUx * grad (Ux)
************************************************************************

      SUBROUTINE SUPGPX (NEQ,U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,
     *                  A,NA,KCOLA,KLDA,KVERT,KMID,DCORVG,TRIA,
     *                  ELE,ICUBN,
     *                  BFULMT, IMASS, ISTOK, UPSAM, RE,  
     *                  IDEF,DCMASS,THSTEP)
                   
      IMPLICIT NONE
      
C include necessary COMMON blocks

      INCLUDE 'cout.inc' 
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      INCLUDE 'ccub.inc'
      
C parameters
      
      INTEGER NEQ
      DOUBLE PRECISION A(*)
      DOUBLE PRECISION U1L1(NEQ),U1L2(NEQ),U2L1(NEQ),U2L2(NEQ)
      DOUBLE PRECISION U1(NEQ),U2(NEQ),D1(NEQ),D2(NEQ)
      DOUBLE PRECISION DCORVG(2,*)
      DOUBLE PRECISION DENTRY(NNBAS,NNBAS)
      DOUBLE PRECISION A1L,A2L,DCMASS
      INTEGER NA
      INTEGER KCOLA(*),KLDA(*)
      INTEGER KVERT(NNVE,*),KMID(NNVE,*)
      INTEGER KENTRY(NNBAS,NNBAS)
      INTEGER IDEF
      INTEGER KDFG(NNBAS), KDFL(NNBAS), IDFL
      
      INTEGER TRIA(SZTRIA),ICUBN
      
      LOGICAL BFULMT
      
      INTEGER ISTOK, IMASS
      DOUBLE PRECISION UPSAM, RE
      
      DOUBLE PRECISION THSTEP
      
      EXTERNAL ELE

C externals

      INTEGER NDFL
      EXTERNAL NDFL

C local variables

      INTEGER ICUB,I,IELTYP,IEQ,JDOFE,ILD,JCOL,JCOL0,IDOFE,IDFG
      INTEGER IVE,JP,JDFL,JDFG,IDOFEH,JDOFEH,IA
      DOUBLE PRECISION DNY,CT0,DUMAX,DU1,DU2,DUNORM,DELTA,DUMAXR
      DOUBLE PRECISION DJ1,DJ2,DJ3,DJ4
      DOUBLE PRECISION XI1,XI2
      DOUBLE PRECISION OM,HBAS,HBASJ1,HBASJ2,HBASJ3,HSUMJ,AH
      DOUBLE PRECISION HBASI1,HBASI2,HBASI3,HSUMI,DENTH

C     In case that IPRECA=4, create the Laplace-part of the system
C     and include it into the system matrix each time this routine
C     is called. Otherwise don't create the Laplace part, just
C     add the nonlinearity.
C     This handling is simply realized by setting the factor NU in
C     front of the Laplace part to 0.

      IF (BFULMT) THEN
        DNY=1D0/RE
      ELSE
        DNY=0D0
      ENDIF

C     A similar handling holds for the case that the (full!) mass
C     matrix is to be included into the system matrix while the
C     matrix is to be rebuild completely. 
C     DCMASS represents the factor in front of the mass matrix.
C     Note that we cannot use DCMASS directly! In the Theta scheme,
C     we typically calculate 
C         [ DCMASS*M + THSTEP*(Laplace) + ... ]
C     so we have to weight everything except for the mass matrix!
C     We make a little trick here to realize that. We always weight
C     everything by THSTEP including the mass matrix - but divide
C     DCMASS by THSTEP before to compensate that!
C         THSTEP * [ CT0*M + (Laplace) + ... ]
C     with CT0 = DCMASS/THSTEP.
C
C     If only the nonlinear part is to be build (which is the normal
C     case, as rebuilding the comlete matrix including the full mass
C     matrix is normally not necessary), simply set CT0=0 which 
C     prevents calculating the mass matrix.

      IF (BFULMT.AND.(IMASS.EQ.1)) THEN
        CT0=DCMASS/THSTEP
      ELSE
        CT0=0D0
      ENDIF

C     Initialize BDER for our element. We want the element to calculate
C     function values as well as first X/Y derivatives:

      DO I = 1,NNDER
        BDER(I)=.FALSE.
      END DO

      DO I=1,3
        BDER(I)=.TRUE.
      END DO

C     Ask the element about its type:

      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      
C     Get the local number of degrees of freedom:
      
      IDFL=NDFL(IELTYP)
      
C     Initialize the cubature formula identifier in the COMMON block
C     with our cubature formula we have to use here:

      ICUB=ICUBN
      CALL CB2Q(ICUB)
      
      IF (IER.NE.0) GOTO 99999

************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************

C     Initialize the values in the cubature points. Allowes the element
C     to precalculate some information for faster access later

      CALL ELE(0D0,0D0,-2)

C     Calculate the maximum norm of the actual velocity field
C     U = A1*U1 + A2*U2 into DUMAX. 
C     Round up the norm to 1D-8 if it's too small...

      DUMAX=0D0
      IF (A2L.EQ.0D0) THEN
        DO IEQ=1,NEQ
          DU1=U1L1(IEQ)
          DU2=U1L2(IEQ)
          DUNORM=SQRT(DU1**2+DU2**2)
          DUMAX=MAX(DUMAX,DUNORM)
        END DO
      ELSE       
        DO IEQ=1,NEQ
          DU1=A1L*U1L1(IEQ)+A2L*U2L1(IEQ)
          DU2=A1L*U1L2(IEQ)+A2L*U2L2(IEQ)
          DUNORM=SQRT(DU1**2+DU2**2)
          DUMAX=MAX(DUMAX,DUNORM)
        END DO
      ENDIF       

      IF (DUMAX.LT.1D-8) DUMAX=1D-8
      DUMAXR = 1D0/DUMAX

C *** Loop over all elements

      DO IEL=1,TRIA(ONEL)
      
C       Calculate the local and global degrees of freedom on the
C       current element IEL:
      
        CALL NDFGLX(TRIA,IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
        IF (IER.LT.0) GOTO 99999

C       Determine local DELTA for streamline-diffusion
C       (cf. p. 121 in Turek's CFD book).
C
C       For Stokes flow, we have the equation
C
C                -nu*div(u) + grad(p) = f
C       
C       not containing any convective part. Therefore we don't need
C       any stabilization technique. So in this case, switch of the
C       stabilization by setting DELTA to 0:

        IF (ISTOK.NE.1) THEN
          CALL DELSDX (U1L1,U1L2,U2L1,U2L2,A1L,A2L,IEL,DUMAXR,DELTA,
     *                 KVERT,KMID,DCORVG,KDFG,IDFL,UPSAM,RE)
        ELSE
          DELTA=0D0
        END IF

C       Determine entry positions in matrix.
C
C       Here we build a small matrix DENTRY/KENTRY which
C       corresponds to the current element. We will assemble the
C       contributions of our current element into this matrix and
C       will add the result into the main matrix later.
C
C       To successfully integrate the contributions into the main
C       matrix, we compute in advance the positions in the main
C       matrix where we have to add the contribution to.
C       KENTRY(X,Y) corresponds to the index in the array A
C       of the entry in line Y, row KCOL(line start + X).

        DO JDOFE=1,IDFL
          ILD=KLDA(KDFG(JDOFE))
          KENTRY(JDOFE,JDOFE)=ILD
          DENTRY(JDOFE,JDOFE)=0D0
          JCOL0=ILD
          DO IDOFE=1,IDFL
            IF (IDOFE.NE.JDOFE) THEN
              IDFG=KDFG(IDOFE)
              DO JCOL=JCOL0,NA
                IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
              END DO
113           JCOL0=JCOL+1
              KENTRY(JDOFE,IDOFE)=JCOL
              DENTRY(JDOFE,IDOFE)=0D0
            END IF
          END DO
        END DO

C       Calculate auxiliary Jacobian factors for the transformation 
C       onto the reference element. See QTRAFO.F for details...

        DO IVE = 1, TRIA(ONVE)
          JP=KVERT(IVE,IEL)
          KVE(IVE)=JP
          DX(IVE)=DCORVG(1,JP)
          DY(IVE)=DCORVG(2,JP)
        END DO

        DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
        DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
        DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
        DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))

C       Loop over the cubature points in our element to calculate
C       its contribution to each of the degrees of freedom on our
C       element

        DO ICUBP = 1, NCUBP
        
C         Get the coordinates of the cubature point on the reference
C         element
        
          XI1=DXI(ICUBP,1)
          XI2=DXI(ICUBP,2)

C         Calculate the Jacobian of the bilinear mapping onto the
C         reference element and the weight OM of the cubature formula:

          DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
          DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
          DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
          DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
          DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)
          OM=DOMEGA(ICUBP)*DETJ

C         Call the element to calculate the values in the current
C         cubature point on the reference element:

          CALL ELE(XI1,XI2,-3)
          IF (IER.LT.0) GOTO 99999

C         Now we have to assemble the "local" matrix DENTRY/KENTRY.
C         This assembling decomposes now into different parts,
C         depending on what has to me assembled.
C
C         We want to set up the nonlinear part of the matrix
C
C           n~_h (u_h, u_h, v_h) 
C
C         = n_h (u_h, u_h, v_h) + sum_T ( delta_T ( u_h*grad u_h, u_h*grad v_h)_T )
C           ^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C          standard nonlin. part                  stabilization
C
C         More precisely, as we want to assemble the matrix which is 
C         later multiplied with coefficient vectors, we have to insert
C         basis functions in the above terms instead of u_h and v_h.
C         Assuming the representation u_h=sum_j(u_j*Phi_j) and 
C         v_h=sum_i(u_i,Phi_i), the above term is evaluated in the
C         DOF's as:
C         
C           n_h (u_h, Phi_j, Phi_i) 
C         + sum_T ( delta_T ( u_h*grad Phi_j, u_h*grad Phi_i )_T )
C
C         In nonstationary simulations, the system matrix typically
C         contains a mass matrix to respect the time derivative.
C         The matrix has the form
C
C         [  dcmass*M*I  +  THWEIG * (-nu * Laplace(.))  ] + THWEIG * u grad(.)
C
C         so if DCMASS<>0, incorporate the (real, not lumped!) mass matrix
C         into the local matrix.

          IF (DCMASS.NE.0D0) THEN

C           Calculate the actual velocity in the current cubature point
C           into (DU1,DU2). If we only have a primary velocity field
C           (A2L=0), we can calculate that only by summing up the
C           velocities in U1Lx, otherwise we have to sum up
C           A1*U1Lx + A2*U2Lx.

            DU1=0D0
            DU2=0D0
            IF (A2L.EQ.0D0) THEN
              DO JDFL=1,IDFL
                HBAS=DBAS(KDFL(JDFL),1)
                IF (ABS(HBAS).GE.1D-8) THEN
                  JDFG=KDFG(JDFL)
                  DU1=DU1+U1L1(JDFG)*HBAS
                  DU2=DU2+U1L2(JDFG)*HBAS
                ENDIF
              END DO
            ELSE
              DO JDFL=1,IDFL
                HBAS=DBAS(KDFL(JDFL),1)
                IF (ABS(HBAS).GE.1D-8) THEN
                  JDFG=KDFG(JDFL)
                  DU1=DU1+(A1L*U1L1(JDFG)+A2L*U2L1(JDFG))*HBAS
                  DU2=DU2+(A1L*U1L2(JDFG)+A2L*U2L2(JDFG))*HBAS
                ENDIF
              END DO
            ENDIF

C           We take a more detailed look onto the last scalar product
C           of n~_h (u_h, u_h, v_h) what we want to calculate here.
C
C           The vector u_h=(DU1,DU2) contains both velocity components,
C           for the X as well as for the Y velocity. On the other hand
C           the system matrix we want to build here will be designed for 
C           one velocity component only! Therefore, Phi_i and Phi_j
C           are scalar functions, so grad(Phi_i), grad(Phi_j) are vectors
C           with two components. Therefore, the last scalar product is more 
C           in detail:
C
C               ( u_h*grad Phi_j, u_h*grad Phi_i )_T
C
C           =   ( < (DU1) , (grad(Phi_j)_1) > , < (DU1) , (grad(Phi_i)_1) > )_T
C                   (DU2) , (grad(Phi_j)_2)       (DU2) , (grad(Phi_i)_2)  
C
C           =   < (DU1) , (grad(Phi_j)_1) >  *  < (DU1) , (grad(Phi_j)_1) >
C                 (DU2) , (grad(Phi_j)_2)         (DU2) , (grad(Phi_j)_2)
C
C           =   HSUMJ * HSUMI
C
C           i.e. a product of two scalar values!
C
C           Summing up over all pairs of multiindices.
C
C           Outer loop over the DOF's j=1..IDFL on our current element, 
C           which corresponds to the basis functions Phi_j:

            DO JDOFE=1,IDFL
            
C             Fetch the contributions of the basis functions Phi_j for
C             function value and first derivatives for the current
C             DOF into HBASxy:
            
              JDOFEH=KDFL(JDOFE)
              HBASJ1=DBAS(JDOFEH,1)
              HBASJ2=DBAS(JDOFEH,2)
              HBASJ3=DBAS(JDOFEH,3)

C             Calculate 
C
C                 U * grad(Phi_j)  =  < grad(Phi_j), U >
C     
C               = ( grad(Phi_j)_1 , (DU1) )
C                 ( grad(Phi_j)_2   (DU2) )
              
              HSUMJ = HBASJ2*DU1 + HBASJ3*DU2

C             Inner loop over the DOF's i=1..IDFL, which corresponds to
C             the basis function Phi_i:

              DO IDOFE=1,IDFL
              
                IF (IDOFE.EQ.JDOFE) THEN
                
C                 Short version of the evaluation of the matrix
C                 contribution - see below for a more detailed
C                 description what is added together here!
                
                  AH = HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *               + DNY*(HBASJ2**2+HBASJ3**2)
     *               + CT0*HBASJ1**2
     
                ELSE
                
C                 Fetch the contributions of the basis function Phi_i for
C                 function value and first derivatives for the current
C                 DOF into HBASIy:
                
                  IDOFEH=KDFL(IDOFE)
                  HBASI1=DBAS(IDOFEH,1)
                  HBASI2=DBAS(IDOFEH,2)
                  HBASI3=DBAS(IDOFEH,3)

C                 Calculate 
C
C                     U * grad(Phi_i)  =  < grad(Phi_i), U >
C     
C                   = ( grad(Phi_i)_1 , (DU1) )
C                     ( grad(Phi_i)_2   (DU2) )

                  HSUMI=HBASI2*DU1+HBASI3*DU2
     
C                 Finally calculate the contribution to the system
C                 matrix. Depending on the configuration of DNY,
C                 IPRECA, DCMASS,... this decomposes into three
C                 different parts:
C
C                 AH = n~_h(u_h,phi_j,phi_i)        | nonlinear part
C                    + DNY*(grad(phi_j,grad(phi_i)) | -nu*Laplace(u)
C                    + CT0*(phi_j*phi_i)            | Mass matrix
C
C                 The last two parts are probably not added to the
C                 matrix by setting DNY or CT0 to 0, respectively.
C
C                 For saving some numerical operations, we write:
C     
C                     HSUMI * (Delta * HSUMJ + HBASJ1)
C
C                 =   Delta * HSUMI * HSUMJ  
C                   + HSUMI * HBASJ1
C     
C                 =   Delta * ( U*grad(Phi_j), U*grad(Phi_i) )
C                   + (grad(Phi_i)*U,Phi_j)
C
C               <->   sum_T ( delta_T ( u_h*grad Phi_j, u_h*grad Phi_i) )_T
C                   + n_h (u_h, Phi_j, Phi_i)
                  
                  AH = HSUMI*(DELTA*HSUMJ+HBASJ1)
     *               + DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *               + CT0*HBASJ1*HBASI1
     
                ENDIF ! (IDOFE.EQ.JDOFE)

C               Weighten the calculated value AH by the cubature
C               weight OM and add it to the local matrix. After the
C               loop over all DOF's is finished, each entry contains
C               the calculated integral.

                DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+OM*AH
                
              END DO ! IDOFE
              
            END DO ! JDOFE

          ELSE

C           Coefficient in front of the mass matrix is 0.
C           Subtract the mass matrix from the system matrix.
C
C           Outer loop over the DOF's j=1..IDFL corresponding
C           to Phi_j:

            DO JDOFE=1,IDFL
            
C             Fetch the function value in the current DOF:
            
              HBASJ1=DBAS(KDFL(JDOFE),1)

C             Inner loop over the DOF's i=1..IDFL corresponding
C             to Phi_i:

              DO IDOFE=1,IDFL
                
C               Fetch the function value of the other DOF of the
C               current element
                
                HBASI1=DBAS(KDFL(IDOFE),1)
                
C               Calculate the contribution for the entry. The factor
C               THSTEP is compensated later when the local matrix
C               is included into the global matrix and/or in the
C               modification of the RHS vector.
                
                AH=-1D0/THSTEP*HBASJ1*HBASI1
                
C               ...and incorporate it into the local matrix.
                
                DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+OM*AH
                
              END DO ! IDOFE
              
            END DO ! JDOFE

          END IF ! (DCMASS.NE.0D0)

        END DO ! ICUBP
      
C       Now we have set up a "local" system matrix. We can either
C       include it into the real matrix or we can use it to simply
C       modify the RHS vector to create a defect vector (throwing
C       away the information about the matrix afterwards, which would
C       result in a matrix free modification of the RHS vector).
      
        DO JDOFE=1,IDFL
        
          DO IDOFE=1,IDFL
          
C           Get the entry from the local matrix and weight it according
C           to the current THETA in the Theta scheme, given by the
C           parameter. 
C           (Remark: For stationary simulations, THSTEP is typically
C            1D0 which includes the local matrix into the global one
C            directly).
          
            DENTH=THSTEP*DENTRY(JDOFE,IDOFE)

C           For IDEF=0/1, incorporate our "local" system matrix into 
C           the global matrix. The position of each entry DENTRY(X,Y) 
C           in the global matrix array A was saved in element KENTRY(X,Y)
C           before.

            IF (IDEF.LT.2) THEN
              IA   =KENTRY(JDOFE,IDOFE)
              A(IA)=A(IA)+DENTH
            ENDIF

C           For IDEF=1,2, build the defect vector
C               D = RHS - A*U
C           This is done matrix free, only with the help of the local 
C           matrix.
C           In this case, D=(D1,D2) is expected to be the RHS on
C           entry and will be updated to be the defect vector when
C           this routine is left.

            IF (IDEF.GT.0) THEN 
              IDFG=KDFG(IDOFE)
              JDFG=KDFG(JDOFE)
              D1(JDFG)= D1(JDFG)-DENTH*U1(IDFG)
              D2(JDFG)= D2(JDFG)-DENTH*U2(IDFG)
            ENDIF 

          END DO ! IDOFE
          
        END DO ! JDOFE

C       Matrix/defect vector updated for that element; proceed to the
C       next one...

      END DO ! IEL

99999 END

************************************************************************
* Streamline diffusion
*
* Adds the SUPG-part on matrix block A after it was initialized by the
* linear part or builds up the complete nonlinear system matrix.
* The input vector Ui is the old velocity field.
* The input vectors UjLi are the transport directions.
*
* Nonparametric version
*
* Extended calling convention
*
* In:
*   A, NA,
*   KCOLA,
*   LLDA   - array [1..*] of double/integer
*            Structure arrays of system matrix, maybe initialised
*            with linear parts of the nonlinear system matrix.
*   KVERT,
*   KMID,
*   DCORVG - usual geometry information; must correspond to TRIA
*   
*   TRIA   - array [1..SZTRIA] of integer
*            Triangulation structure of the underlying mesh
*
*   NEQ    - length of solution/defect vectors
*
*   U1L1,
*   U1L2   - array [1..NEQ] of double
*            main velocity field used for assembling of 
*            the nonlinearity. Can be undefined if A1L=0.
*   U2L1,
*   U2L2   - array [1..NEQ] of double
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
*                A can be undefined (but not KCOLA,KLDA!).
*                Can be used to modify the defect vector matrix
*                free. Note that the rules how to set up the
*                system matrix (NY, DCMASS,...) still hold,
*                although the matrix itself is not modified!
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
*   BFULMT - Whether or not to generate the full nonlinear system
*            matrix.
*            =false: generate only the nonlinear part 
*                         THWEIG * u grad(.)
*                    of the system matrix and add it to the existing
*                    matrix in A
*            =true : Generate the full system matrix
*               [  M*I  +  THWEIG * (-nu * Laplace(.) * u * grad(.)) ] 
*                    and add it to A
*            
*   THSTEP - weighting factor of the nonlinearity N(u) which is to be
*            added to the matrix
*            (is used as step-size of a Theta-scheme)
*
*   DCMASS - This defines how to modify the system matrix - if the
*            matrix is modified at all (see IDEF!).
*             =0: subtract the mass matrix from the system matrix;
*                 The nonlinearity is not build.
*            <>0: Add the nonlinear term including the stabilization
*                 to the system matrix.
*                 If the full matrix with not-lumped mass matrix 
*                 is to be build (IPRECA=4 and IMASS=1),
*                 add DCMASS*(Mass matrix) to the system matrix.
*
*   UPSAM  - control parameter of the streamline diffusion
*   RE     - 1/nu = viscosity
*   ISTOK  - =1 if the Stokes equation is to be discretized
*            rather than Navier-Stokes.
*   IMASS  - if BFULMT=TRUE:
*            =1: add the real mass matrix to the system matrix
*            =0: don't add the mass matrix to the system matrix
*
*   ELE    - used element for the discretisation
*   ICUBN  - Cubature formula to use for the discretization
*
* Out:
*   A      - system matrix; only if IDEF=0,1.
*            The nonlinearity is added to that matrix or
*            the linearized nonlinear system matrix is build
*            completely, respectively.
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
*   If BFULMT=true, IMASS=1, DCMASS<>0, D is updated according to
*
*             D = D - DCMASS*M*UUx - THSTEP * UUx * grad (Ux)
*   
*   If BFULMT=true, IMASS=1, DCMASS=0, D is updated according to
*
*             D = D + M*UUx - THSTEP * UUx * grad (Ux)
************************************************************************

      SUBROUTINE SUPGNX (NEQ,U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,
     *                  A,NA,KCOLA,KLDA,KVERT,KMID,DCORVG,TRIA,
     *                  ELE,ICUBN,
     *                  BFULMT, IMASS, ISTOK, UPSAM, RE,  
     *                  IDEF,DCMASS,THSTEP)

      IMPLICIT NONE
      
C include necessary COMMON blocks

      INCLUDE 'cout.inc' 
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      INCLUDE 'ccub.inc'
      
C parameters
      
      INTEGER NEQ
      DOUBLE PRECISION A(*)
      DOUBLE PRECISION U1L1(NEQ),U1L2(NEQ),U2L1(NEQ),U2L2(NEQ)
      DOUBLE PRECISION U1(NEQ),U2(NEQ),D1(NEQ),D2(NEQ)
      DOUBLE PRECISION DCORVG(2,*)
      DOUBLE PRECISION DENTRY(NNBAS,NNBAS)
      DOUBLE PRECISION A1L,A2L,DCMASS
      INTEGER NA
      INTEGER KCOLA(*),KLDA(*)
      INTEGER KVERT(NNVE,*),KMID(NNVE,*)
      INTEGER KENTRY(NNBAS,NNBAS)
      INTEGER IDEF
      INTEGER KDFG(NNBAS), KDFL(NNBAS), IDFL
      
      INTEGER TRIA(SZTRIA),ICUBN
      
      LOGICAL BFULMT
      
      INTEGER ISTOK, IMASS
      DOUBLE PRECISION UPSAM, RE
      
      DOUBLE PRECISION THSTEP
      
      EXTERNAL ELE

C externals

      INTEGER NDFL
      EXTERNAL NDFL

C local variables

      INTEGER ICUB,I,IELTYP,IEQ,JDOFE,ILD,JCOL,JCOL0,IDOFE,IDFG
      INTEGER IVE,JP,JDFL,JDFG,IDOFEH,JDOFEH,IA
      DOUBLE PRECISION DNY,CT0,DUMAX,DU1,DU2,DUNORM,DELTA,DUMAXR
      DOUBLE PRECISION DJ1,DJ2,DJ3,DJ4
      DOUBLE PRECISION XI1,XI2,XX,YY
      DOUBLE PRECISION OM,HBAS,HBASJ1,HBASJ2,HBASJ3,HSUMJ,AH
      DOUBLE PRECISION HBASI1,HBASI2,HBASI3,HSUMI,DENTH

C     In case that IPRECA=4, create the Laplace-part of the system
C     and include it into the system matrix each time this routine
C     is called. Otherwise don't create the Laplace part, just
C     add the nonlinearity.
C     This handling is simply realized by setting the factor NU in
C     front of the Laplace part to 0.
      
      IF (BFULMT) THEN
        DNY=1D0/RE
      ELSE
        DNY=0D0
      ENDIF

C     A similar handling holds for the case that the (full!) mass
C     matrix is to be included into the system matrix while the
C     matrix is to be rebuild completely. 
C     DCMASS represents the factor in front of the mass matrix.
C     Note that we cannot use DCMASS directly! In the Theta scheme,
C     we typically calculate 
C         [ DCMASS*M + THSTEP*(Laplace) + ... ]
C     so we have to weight everything except for the mass matrix!
C     We make a little trick here to realize that. We always weight
C     everything by THSTEP including the mass matrix - but divide
C     DCMASS by THSTEP before to compensate that!
C         THSTEP * [ CT0*M + (Laplace) + ... ]
C     with CT0 = DCMASS/THSTEP.
C
C     If only the nonlinear part is to be build (which is the normal
C     case, as rebuilding the comlete matrix including the full mass
C     matrix is normally not necessary), simply set CT0=0 which 
C     prevents calculating the mass matrix.

      IF (BFULMT.AND.(IMASS.EQ.1)) THEN
        CT0=DCMASS/THSTEP
      ELSE
        CT0=0D0
      ENDIF

C     Initialize BDER for our element. We want the element to calculate
C     function values as well as first X/Y derivatives:

      DO I = 1,NNDER
        BDER(I)=.FALSE.
      END DO

      DO I=1,3
        BDER(I)=.TRUE.
      END DO

C     Ask the element about its type:

      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      
C     Get the local number of degrees of freedom:
      
      IDFL=NDFL(IELTYP)
      
C     Initialize the cubature formula identifier in the COMMON block
C     with our cubature formula we have to use here:

      ICUB=ICUBN
      CALL CB2Q(ICUB)
      
      IF (IER.NE.0) GOTO 99999

************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************

C     Calculate the maximum norm of the actual velocity field
C     U = A1*U1 + A2*U2 into DUMAX. 
C     Round up the norm to 1D-8 if it's too small...

      DUMAX=0D0
      IF (A2L.EQ.0D0) THEN
        DO IEQ=1,NEQ
          DU1=U1L1(IEQ)
          DU2=U1L2(IEQ)
          DUNORM=SQRT(DU1**2+DU2**2)
          DUMAX=MAX(DUMAX,DUNORM)
        END DO
      ELSE       
        DO IEQ=1,NEQ
          DU1=A1L*U1L1(IEQ)+A2L*U2L1(IEQ)
          DU2=A1L*U1L2(IEQ)+A2L*U2L2(IEQ)
          DUNORM=SQRT(DU1**2+DU2**2)
          DUMAX=MAX(DUMAX,DUNORM)
        END DO
      ENDIF       

      IF (DUMAX.LT.1D-8) DUMAX=1D-8
      DUMAXR = 1D0/DUMAX

C *** Loop over all elements

      DO IEL=1,TRIA(ONEL)
      
C       Calculate the local and global degrees of freedom on the
C       current element IEL:
      
        CALL NDFGLX(TRIA,IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
        IF (IER.LT.0) GOTO 99999

C       Determine local DELTA for streamline-diffusion
C       (cf. p. 121 in Turek's CFD book).
C
C       For Stokes flow, we have the equation
C
C                -nu*div(u) + grad(p) = f
C       
C       not containing any convective part. Therefore we don't need
C       any stabilization technique. So in this case, switch of the
C       stabilization by setting DELTA to 0:

        IF (ISTOK.NE.1) THEN
          CALL DELSDX (U1L1,U1L2,U2L1,U2L2,A1L,A2L,IEL,DUMAXR,DELTA,
     *                 KVERT,KMID,DCORVG,KDFG,IDFL,UPSAM,RE)
        ELSE
          DELTA=0D0
        END IF

C       Determine entry positions in matrix.
C
C       Here we build a small matrix DENTRY/KENTRY which
C       corresponds to the current element. We will assemble the
C       contributions of our current element into this matrix and
C       will add the result into the main matrix later.
C
C       To successfully integrate the contributions into the main
C       matrix, we compute in advance the positions in the main
C       matrix where we have to add the contribution to.
C       KENTRY(X,Y) corresponds to the index in the array A
C       of the entry in line Y, row KCOL(line start + X).

        DO JDOFE=1,IDFL
          ILD=KLDA(KDFG(JDOFE))
          KENTRY(JDOFE,JDOFE)=ILD
          DENTRY(JDOFE,JDOFE)=0D0
          JCOL0=ILD
          DO IDOFE=1,IDFL
            IF (IDOFE.NE.JDOFE) THEN
              IDFG=KDFG(IDOFE)
              DO JCOL=JCOL0,NA
                IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
              END DO
113           JCOL0=JCOL+1
              KENTRY(JDOFE,IDOFE)=JCOL
              DENTRY(JDOFE,IDOFE)=0D0
            END IF
          END DO
        END DO

C       Calculate auxiliary Jacobian factors for the transformation 
C       onto the reference element. See QTRAFO.F for details...

        DO IVE = 1, TRIA(ONVE)
          JP=KVERT(IVE,IEL)
          KVE(IVE)=JP
          DX(IVE)=DCORVG(1,JP)
          DY(IVE)=DCORVG(2,JP)
        END DO

        DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
        DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
        DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
        DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))

C       Initialize the values in the cubature points. Allowes the element
C       to precalculate some information for faster access later

        CALL ELE(0D0,0D0,-2)
        IF (IER.NE.0) GOTO 99999

C       Loop over the cubature points in our element to calculate
C       its contribution to each of the degrees of freedom on our
C       element

        DO ICUBP = 1, NCUBP
        
C         Get the coordinates of the cubature point on the reference
C         element
        
          XI1=DXI(ICUBP,1)
          XI2=DXI(ICUBP,2)

C         Calculate the Jacobian of the bilinear mapping onto the
C         reference element and the weight OM of the cubature formula:

          DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
          DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
          DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
          DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
          DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)
          OM=DOMEGA(ICUBP)*DETJ
          
C         Calculate the real coordinates of the cubature point
C         with the bilinear transformation

          XX=0.5D0*(DX(1)+DX(2)+DJ1)+0.5D0*(DX(2)-DX(1)+DJ2)*XI1
     *      +0.5D0*DJ1*XI2+0.5D0*DJ2*XI1*XI2
          YY=0.5D0*(DY(1)+DY(3)+DJ3)+0.5D0*DJ4*XI1
     *      +0.5D0*(DY(3)-DY(1)-DJ4)*XI2-0.5D0*DJ3*XI1*XI2

C         Call the element to calculate the values in the current
C         cubature point on the reference element:

          CALL ELE(XX,YY,-3)
          IF (IER.LT.0) GOTO 99999

C         Now we have to assemble the "local" matrix DENTRY/KENTRY.
C         This assembling decomposes now into different parts,
C         depending on what has to me assembled.
C
C         We want to set up the nonlinear part of the matrix
C
C           n~_h (u_h, u_h, v_h) 
C
C         = n_h (u_h, u_h, v_h) + sum_T ( delta_T ( u_h*grad u_h, u_h*grad v_h)_T )
C           ^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C          standard nonlin. part                  stabilization
C
C         More precisely, as we want to assemble the matrix which is 
C         later multiplied with coefficient vectors, we have to insert
C         basis functions in the above terms instead of u_h and v_h.
C         Assuming the representation u_h=sum_j(u_j*Phi_j) and 
C         v_h=sum_i(u_i,Phi_i), the above term is evaluated in the
C         DOF's as:
C         
C           n_h (u_h, Phi_j, Phi_i) 
C         + sum_T ( delta_T ( u_h*grad Phi_j, u_h*grad Phi_i )_T )
C
C         In nonstationary simulations, the system matrix typically
C         contains a mass matrix to respect the time derivative.
C         The matrix has the form
C
C         [  dcmass*M*I  +  THWEIG * (-nu * Laplace(.))  ] + THWEIG * u grad(.)
C
C         so if DCMASS<>0, incorporate the (real, not lumped!) mass matrix
C         into the local matrix.

          IF (DCMASS.NE.0D0) THEN

C           Calculate the actual velocity in the current cubature point
C           into (DU1,DU2). If we only have a primary velocity field
C           (A2L=0), we can calculate that only by summing up the
C           velocities in U1Lx, otherwise we have to sum up
C           A1*U1Lx + A2*U2Lx.

            DU1=0D0
            DU2=0D0
            IF (A2L.EQ.0D0) THEN
              DO JDFL=1,IDFL
                HBAS=DBAS(KDFL(JDFL),1)
                IF (ABS(HBAS).GE.1D-8) THEN
                  JDFG=KDFG(JDFL)
                  DU1=DU1+U1L1(JDFG)*HBAS
                  DU2=DU2+U1L2(JDFG)*HBAS
                ENDIF
              END DO
            ELSE
              DO JDFL=1,IDFL
                HBAS=DBAS(KDFL(JDFL),1)
                IF (ABS(HBAS).GE.1D-8) THEN
                  JDFG=KDFG(JDFL)
                  DU1=DU1+(A1L*U1L1(JDFG)+A2L*U2L1(JDFG))*HBAS
                  DU2=DU2+(A1L*U1L2(JDFG)+A2L*U2L2(JDFG))*HBAS
                ENDIF
              END DO
            ENDIF

C           We take a more detailed look onto the last scalar product
C           of n~_h (u_h, u_h, v_h) what we want to calculate here.
C
C           The vector u_h=(DU1,DU2) contains both velocity components,
C           for the X as well as for the Y velocity. On the other hand
C           the system matrix we want to build here will be designed for 
C           one velocity component only! Therefore, Phi_i and Phi_j
C           are scalar functions, so grad(Phi_i), grad(Phi_j) are vectors
C           with two components. Therefore, the last scalar product is more 
C           in detail:
C
C               ( u_h*grad Phi_j, u_h*grad Phi_i )_T
C
C           =   ( < (DU1) , (grad(Phi_j)_1) > , < (DU1) , (grad(Phi_i)_1) > )_T
C                   (DU2) , (grad(Phi_j)_2)       (DU2) , (grad(Phi_i)_2)  
C
C           =   < (DU1) , (grad(Phi_j)_1) >  *  < (DU1) , (grad(Phi_j)_1) >
C                 (DU2) , (grad(Phi_j)_2)         (DU2) , (grad(Phi_j)_2)
C
C           =   HSUMJ * HSUMI
C
C           i.e. a product of two scalar values!
C
C           Summing up over all pairs of multiindices.
C
C           Outer loop over the DOF's j=1..IDFL on our current element, 
C           which corresponds to the basis functions Phi_j:

            DO JDOFE=1,IDFL
            
C             Fetch the contributions of the basis functions Phi_j for
C             function value and first derivatives for the current
C             DOF into HBASxy:
            
              JDOFEH=KDFL(JDOFE)
              HBASJ1=DBAS(JDOFEH,1)
              HBASJ2=DBAS(JDOFEH,2)
              HBASJ3=DBAS(JDOFEH,3)

C             Calculate 
C
C                 U * grad(Phi_j)  =  < grad(Phi_j), U >
C     
C               = ( grad(Phi_j)_1 , (DU1) )
C                 ( grad(Phi_j)_2   (DU2) )
              
              HSUMJ = HBASJ2*DU1 + HBASJ3*DU2

C             Inner loop over the DOF's i=1..IDFL, which corresponds to
C             the basis function Phi_i:

              DO IDOFE=1,IDFL
              
                IF (IDOFE.EQ.JDOFE) THEN
                
C                 Short version of the evaluation of the matrix
C                 contribution - see below for a more detailed
C                 description what is added together here!
                
                  AH = HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *               + DNY*(HBASJ2**2+HBASJ3**2)
     *               + CT0*HBASJ1**2
     
                ELSE
                
C                 Fetch the contributions of the basis function Phi_i for
C                 function value and first derivatives for the current
C                 DOF into HBASIy:
                
                  IDOFEH=KDFL(IDOFE)
                  HBASI1=DBAS(IDOFEH,1)
                  HBASI2=DBAS(IDOFEH,2)
                  HBASI3=DBAS(IDOFEH,3)

C                 Calculate 
C
C                     U * grad(Phi_i)  =  < grad(Phi_i), U >
C     
C                   = ( grad(Phi_i)_1 , (DU1) )
C                     ( grad(Phi_i)_2   (DU2) )

                  HSUMI=HBASI2*DU1+HBASI3*DU2
     
C                 Finally calculate the contribution to the system
C                 matrix. Depending on the configuration of DNY,
C                 IPRECA, DCMASS,... this decomposes into three
C                 different parts:
C
C                 AH = n~_h(u_h,phi_j,phi_i)        | nonlinear part
C                    + DNY*(grad(phi_j,grad(phi_i)) | -nu*Laplace(u)
C                    + CT0*(phi_j*phi_i)            | Mass matrix
C
C                 The last two parts are probably not added to the
C                 matrix by setting DNY or CT0 to 0, respectively.
C
C                 For saving some numerical operations, we write:
C     
C                     HSUMI * (Delta * HSUMJ + HBASJ1)
C
C                 =   Delta * HSUMI * HSUMJ  
C                   + HSUMI * HBASJ1
C     
C                 =   Delta * ( U*grad(Phi_j), U*grad(Phi_i) )
C                   + (grad(Phi_i)*U,Phi_j)
C
C               <->   sum_T ( delta_T ( u_h*grad Phi_j, u_h*grad Phi_i) )_T
C                   + n_h (u_h, Phi_j, Phi_i)
                  
                  AH = HSUMI*(DELTA*HSUMJ+HBASJ1)
     *               + DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *               + CT0*HBASJ1*HBASI1
     
                ENDIF ! (IDOFE.EQ.JDOFE)

C               Weighten the calculated value AH by the cubature
C               weight OM and add it to the local matrix. After the
C               loop over all DOF's is finished, each entry contains
C               the calculated integral.

                DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+OM*AH
                
              END DO ! IDOFE
              
            END DO ! JDOFE

          ELSE

C           Coefficient in front of the mass matrix is 0.
C           Subtract the mass matrix from the system matrix.
C
C           Outer loop over the DOF's j=1..IDFL corresponding
C           to Phi_j:

            DO JDOFE=1,IDFL
            
C             Fetch the function value in the current DOF:
            
              HBASJ1=DBAS(KDFL(JDOFE),1)

C             Inner loop over the DOF's i=1..IDFL corresponding
C             to Phi_i:

              DO IDOFE=1,IDFL
                
C               Fetch the function value of the other DOF of the
C               current element
                
                HBASI1=DBAS(KDFL(IDOFE),1)
                
C               Calculate the contribution for the entry. The factor
C               THSTEP is compensated later when the local matrix
C               is included into the global matrix and/or in the
C               modification of the RHS vector.
                
                AH=-1D0/THSTEP*HBASJ1*HBASI1
                
C               ...and incorporate it into the local matrix.
                
                DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+OM*AH
                
              END DO ! IDOFE
              
            END DO ! JDOFE

          END IF ! (DCMASS.NE.0D0)

        END DO ! ICUBP
      
C       Now we have set up a "local" system matrix. We can either
C       include it into the real matrix or we can use it to simply
C       modify the RHS vector to create a defect vector (throwing
C       away the information about the matrix afterwards, which would
C       result in a matrix free modification of the RHS vector).
      
        DO JDOFE=1,IDFL
        
          DO IDOFE=1,IDFL
          
C           Get the entry from the local matrix and weight it according
C           to the current THETA in the Theta scheme, given by the
C           parameter. 
C           (Remark: For stationary simulations, THSTEP is typically
C            1D0 which includes the local matrix into the global one
C            directly).
          
            DENTH=THSTEP*DENTRY(JDOFE,IDOFE)

C           For IDEF=0/1, incorporate our "local" system matrix into 
C           the global matrix. The position of each entry DENTRY(X,Y) 
C           in the global matrix array A was saved in element KENTRY(X,Y)
C           before.

            IF (IDEF.LT.2) THEN
              IA   =KENTRY(JDOFE,IDOFE)
              A(IA)=A(IA)+DENTH
            ENDIF

C           For IDEF=1,2, build the defect vector
C               D = RHS - A*U
C           This is done matrix free, only with the help of the local 
C           matrix.
C           In this case, D=(D1,D2) is expected to be the RHS on
C           entry and will be updated to be the defect vector when
C           this routine is left.

            IF (IDEF.GT.0) THEN 
              IDFG=KDFG(IDOFE)
              JDFG=KDFG(JDOFE)
              D1(JDFG)= D1(JDFG)-DENTH*U1(IDFG)
              D2(JDFG)= D2(JDFG)-DENTH*U2(IDFG)
            ENDIF 

          END DO ! IDOFE
          
        END DO ! JDOFE

C       Matrix/defect vector updated for that element; proceed to the
C       next one...

      END DO ! IEL

99999 END

************************************************************************
* Calculation of a local DELTA
*
* Extended calling convention
*
* This routine calculates a local DELTA=DELTA_T for a finite element
* T=IEL. This can be used by the streamline diffusion stabilization
* technique as a multiplier of the (local) bilinear form.
*
* In:
*   IEL    - Element where the Delta should be calculated
*
*   KVERT,
*   KMID,
*   DCORVG - Usual geometry information
*
*   IDFL   - Number of degrees of freedom on element IEL
*   KDFG   - array [1..IDFL] of integer
*            Array with global degrees of freedom, corresponding to
*            local degrees of freedom 1..IDFL on element IEL.
*   U1L1,
*   U1L2   - array [1..NU] of double
*            Main velocity field. 
*   U2L1,
*   U2L2   - array [1..NU] of double
*            Secondary velocity field. 
*   A1L    - double; weighting factor for U1L1/U1L2
*   A2L    - double; weighting factor for U2L1/U2L2
*   UPSAM  - double; user defined parameter for configuring the 
*            streamline diffusion.
*            < 0: Simple calculation of Delta, using
*                 DELTA = |UPSAM| * h_T
*            > 0: usually UPSAM = 0.1 .. 2;
*                 Samarskji-like calculation of DELTA using:
*                 DELTA = UPSAM * h_t/||u||_T * 2*Re_T/(1+Re_T)
*   NUREC  - Reciprocal value 1/NU of coefficient NU in front of the
*            Laplacian term of the Navier-Stokes equation
*               NU * Laplace(u) + u*grad(u) + ...
*   DUMAXR - Reciprocal of the maximum norm of velocity in the domain:
*            1/DUMAXR = 1/||u||_Omega
*            
*
* Out:
*   DELTA  - local Delta
*
* Remarks:
*
*   The effective velocity that is used for calculating the DELTA
*   is combined by a weighted mean of the two velocity fields U1,U2
*   by:
*                     Ux = A1*U1Lx + A2*U2Lx
*   The coefficients A1,A2 allow the caller to take influence on which
*   velocity field to weight more.
*
************************************************************************

      SUBROUTINE  DELSDX (U1L1,U1L2,U2L1,U2L2,A1L,A2L,IEL,DUMAXR,DELTA,
     *                    KVERT,KMID,DCORVG,KDFG,IDFL,UPSAM,NUREC)

      IMPLICIT NONE
      
C include necessary COMMON blocks
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'cbasictria.inc'
      
C parameters
      
      DOUBLE PRECISION U1L1(*),U1L2(*),U2L1(*),U2L2(*)
      DOUBLE PRECISION A1L,A2L,DUMAXR,DELTA,NUREC,UPSAM
      INTEGER IEL,IDFL,KDFG(*)
      INTEGER KVERT(NNVE,*),KMID(NNVE,*)
      DOUBLE PRECISION DCORVG(2,*)

C local variables
      
      DOUBLE PRECISION HLOCAL,DU1,DU2,UNORM,RELOC
      INTEGER IDOF

C     Loop through the local degrees of freedom on element IEL.
C     Sum up the velocities on these DOF's. This will result
C     in the vector (DU1,DU2) representing the (mean) X/Y-velocity
C     through element IEL.
C
C     For Ex30/Ex31 element, U1/U2 represent the mean velocity
C     along an egde/on the midpoint of each edge, so U1/U2 is
C     clearly an approximation to the velocity in element T.

      DU1=0D0
      DU2=0D0
      DO IDOF=1,IDFL
        DU1=DU1+(A1L*U1L1(KDFG(IDOF))+A2L*U2L1(KDFG(IDOF)))
        DU2=DU2+(A1L*U1L2(KDFG(IDOF))+A2L*U2L2(KDFG(IDOF)))
      END DO

C     Calculate the norm of that local velocity:

      UNORM = SQRT(DU1**2+DU2**2) / DBLE(IDFL)
      
C     Now we have:   UNORM = ||u||_T
C     and:           u_T = a1*u1_T + a2*u2_T
C
C     If the norm of the velocity is small, we choose DELTA = 0,
C     which results in central difference in the streamline diffusion
C     matrix assembling:

      IF (UNORM.LE.1D-8) THEN
      
        DELTA=0D0

      ELSE

C       u_T defines the "slope" of the velocity through
C       the element T. At next, calculate the local mesh width
C       HLOCAL = h = h_T on our element T=IEL:

        CALL HLOCLX (HLOCAL,UNORM, DU1, DU2, IEL, KVERT,KMID,DCORVG)

C       Calculate DELTA... (cf. p. 121 in Turek's CFD book)

        IF (UPSAM.LT.0D0) THEN

C         For UPSAM<0, we use simple calculation of Delta:        
        
          DELTA = ABS(UPSAM)*HLOCAL
          
        ELSE
        
C         For Delta >= 0, we use standard Samarskji-like calculation
C         of Delta. At first calculate the local Reynolds number
C         RELOC = Re_T = ||u||_T * h_T / NU
          
          RELOC=UNORM*HLOCAL*NUREC
          
C         and then the DELTA = UPSAM * h_t/||u|| * 2*Re_T/(1+Re_T)
          
          DELTA = UPSAM * HLOCAL*DUMAXR * 2D0*(RELOC/(1D0+RELOC))
          
        ENDIF ! (UPSAM.LT.0D0)
        
      ENDIF ! (UNORM.LE.1D-8)

      END

************************************************************************
* Determine the local mesh width for an element JEL of a triangulation.
* 
* In:
*   JEL    - Element where the local h should be calculated
*   UNORM  - norm ||u||_T = mean velocity through element T=JEL
*   XBETA1,
*   XBETA2 - mean velocity u_T = (xbeta1,xbeta2) through element T=JEL
*
*   KVERT,
*   KMID,
*   DCORVG - Usual geometry information
*
* Out:
*   HLOCAL - local mesh width
************************************************************************

      SUBROUTINE HLOCLX (HLOCAL, UNORM,  XBETA1, 
     *                   XBETA2, JEL,KVERT,KMID,DCORVG)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'

C parameters

      INTEGER JEL, KVERT(NNVE,*),KMID(NNVE,*)
      DOUBLE PRECISION DCORVG(2,*), HLOCAL, UNORM, XBETA1, XBETA2
      
C local variables
      
      DOUBLE PRECISION LAMBDA
      INTEGER NECK1,NECK2,NECK3,NECK4
      DOUBLE PRECISION X1,Y1,X2,Y2,X3,Y3,X4,Y4
      DOUBLE PRECISION ALPMAX, ALPHA

C     Fetch the numbers of the four corners of element JEL

      neck1=KVERT(1,JEL)
      neck2=KVERT(2,JEL)
      neck3=KVERT(3,JEL)
      neck4=KVERT(4,JEL)

C     Fetch the coordinates of these corners

      x1=DCORVG(1,neck1)
      y1=DCORVG(2,neck1)
      x2=DCORVG(1,neck2)
      y2=DCORVG(2,neck2)
      x3=DCORVG(1,neck3)
      y3=DCORVG(2,neck3)
      x4=DCORVG(1,neck4)
      y4=DCORVG(2,neck4)

C     Scale: (deactivated)
C
C      skal=max(xbeta1,xbeta2)

C      xbeta1=xbeta1
C      xbeta2=xbeta2

      ALPMAX=0D0
      
C     Loop through the four corners of element JEL and check
C     of a line with slope BETA=(xbeta1,xbeta2) starting in this
C     corner really intersects with one of the edges of the element.
C     Remark that we only have to check the two opposite edges
C     to the current corner!
C
C     -----------------------------------------------------------------
C     Check the first corner:

      CALL ISCT2L(X1,Y1,ALPHA,XBETA1,XBETA2,
     *            X3,Y3,LAMBDA,X2,Y2)
      ALPMAX=MAX(ALPHA,ALPMAX)

      CALL ISCT2L(X1,Y1,ALPHA,XBETA1,XBETA2,
     *            X3,Y3,LAMBDA,X4,Y4)
      ALPMAX=MAX(ALPHA,ALPMAX)
      
C     -----------------------------------------------------------------
C     The second one...
      
      CALL ISCT2L(X2,Y2,ALPHA,XBETA1,XBETA2,
     *            X4,Y4,LAMBDA,X1,Y1)
      ALPMAX=MAX(ALPHA,ALPMAX)

      CALL ISCT2L(X2,Y2,ALPHA,XBETA1,XBETA2,
     *            X4,Y4,LAMBDA,X3,Y3)
      ALPMAX=MAX(ALPHA,ALPMAX)
      
C     -----------------------------------------------------------------
C     The third one...
      
      CALL ISCT2L(X3,Y3,ALPHA,XBETA1,XBETA2,
     *            X1,Y1,LAMBDA,X2,Y2)
      ALPMAX=MAX(ALPHA,ALPMAX)

      CALL ISCT2L(X3,Y3,ALPHA,XBETA1,XBETA2,
     *            X1,Y1,LAMBDA,X4,Y4)
      ALPMAX=MAX(ALPHA,ALPMAX)
      
C     -----------------------------------------------------------------
C     And the fourth=last one...
      
      CALL ISCT2L(X4,Y4,ALPHA,XBETA1,XBETA2,
     *            X2,Y2,LAMBDA,X1,Y1)
      ALPMAX=MAX(ALPHA,ALPMAX)

      CALL ISCT2L(X4,Y4,ALPHA,XBETA1,XBETA2,
     *            X2,Y2,LAMBDA,X3,Y3)
      ALPMAX=MAX(ALPHA,ALPMAX)

C     -----------------------------------------------------------------
C     finally determine the local h=h_T

      HLOCAL=ALPMAX*4D0*UNORM

      END

*******************************************************************
* Intersect two lines in R^2
*******************************************************************

      SUBROUTINE ISCT2L (XO,YO,ALPHA,BETA1,BETA2,
     *                   XA,YA,LAMBDA,XB,YB)

      IMPLICIT NONE
      
C parameters

      DOUBLE PRECISION XO, YO, BETA1, BETA2, ALPHA, XA, YA, XB, YB
      DOUBLE PRECISION LAMBDA
      
C local variables

      DOUBLE PRECISION SKAL
      LOGICAL BFLAG

      SKAL=BETA2*(XB-XA)-BETA1*(YB-YA)
      
      if (skal.eq.0D0) then
  
C        beta and the vector are parallel

         ALPHA=0D0
C         WRITE(*,*) 'EINS'
         BFLAG=.FALSE.
      ELSE  
           LAMBDA=(BETA1*(YA-YO)-BETA2*(XA-XO))/SKAL
           BFLAG=.TRUE.      
      ENDIF

C     is the intersection point insode of the element?

      IF (BFLAG) THEN
        IF ((LAMBDA.GE.-1E-1).AND.(LAMBDA.LE.1.11E0)) THEN
          IF (BETA1.NE.0D0) THEN
              ALPHA=((XA-XO)+LAMBDA*(XB-XA))/BETA1
          ELSE
              IF (BETA2.NE.0D0) THEN
                ALPHA=((YA-YO)+LAMBDA*(YB-YA))/BETA2
              ELSE
                ALPHA=0D0
            ENDIF
          ENDIF
        ELSE
C           WRITE(*,*) 'DREI'
          ALPHA=0D0
        ENDIF
      ENDIF
   
      END

************************************************************************
* Streamline diffusion, ALE method
*
* Adds the SUPG-part on matrix block A after it was initialized by the
* linear part or builds up the complete nonlinear system matrix.
* The input vector Ui is the old velocity field.
* The input vectors UjLi are the transport directions.
*
* Parametric version
*
* In contrast to SUPGNX, this routine supports the 1st order ALE
* method. The caller has to specify a grid velocity vector field,
* which gives the velocity of the grid points in each time step.
* This vector field is interpreted piecewise linear with the element
* E011 and subtracted from the main velocity field to calculate
* the nonlinearity (cf. [Duarte, Formaz, Natesan; "Arbitrary Lagrangian-
* Euler Method for Navier-Stokes equations with moving boundaries";
* Comput. Methods Appl. Mech. Engrg. 193 (2004), 4819-4836).
*
* Extended calling convention
*
* In:
*   A, NA,
*   KCOLA,
*   LLDA   - array [1..*] of double/integer
*            Structure arrays of system matrix, maybe initialised
*            with linear parts of the nonlinear system matrix.
*   KVERT,
*   KMID,
*   DCORVG - usual geometry information; must correspond to TRIA
*   
*   TRIA   - array [1..SZTRIA] of integer
*            Triangulation structure of the underlying mesh
*
*   NEQ    - length of solution/defect vectors
*
*   U1L1,
*   U1L2   - array [1..NEQ] of double
*            main velocity field used for assembling of 
*            the nonlinearity. Can be undefined if A1L=0.
*   U2L1,
*   U2L2   - array [1..NEQ] of double
*            secondary velocity field, used for the assembling
*            of the nonlinearity. Can be undefined if A2L=0.
*   A1L    - double; weighting factor for U1L1/U1L2
*   A2L    - double; weighting factor for U2L1/U2L2
*
*   IALE   - =0: don't use ALE, U1MVEL/U2MVEl can be undefined
*            =1: use ALE method; U1MVEL/U2MVEl must be undefined
*   UMVEL  - array [1..2,1..NVT] of double
*            Grid velocity field. UMVEL(1,.) is the X-velocity,
*            UMVEL(2,.) the Y-velocity in each corner grid point.
*
*   IDEF   - Controls the behaviour of this routine.
*            =0: modify system matrix, add nonlinearity.
*                Defect vectors D1,D2 and velocity vectors U1,U2
*                can be undefined.
*            =1: modify both, system matrix and defect vector
*            =2: modify defect vectors, include nonlinearity.
*                A can be undefined (but not KCOLA,KLDA!).
*                Can be used to modify the defect vector matrix
*                free. Note that the rules how to set up the
*                system matrix (NY, DCMASS,...) still hold,
*                although the matrix itself is not modified!
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
*   BFULMT - Whether or not to generate the full nonlinear system
*            matrix.
*            =false: generate only the nonlinear part 
*                         THWEIG * u grad(.)
*                    of the system matrix and add it to the existing
*                    matrix in A
*            =true : Generate the full system matrix
*               [  M*I  +  THWEIG * (-nu * Laplace(.) * u * grad(.)) ] 
*                    and add it to A
*            
*   THSTEP - weighting factor of the nonlinearity N(u) which is to be
*            added to the matrix
*            (is used as step-size of a Theta-scheme)
*
*   DCMASS - This defines how to modify the system matrix - if the
*            matrix is modified at all (see IDEF!).
*             =0: subtract the mass matrix from the system matrix;
*                 The nonlinearity is not build.
*            <>0: Add the nonlinear term including the stabilization
*                 to the system matrix.
*                 If the full matrix with not-lumped mass matrix 
*                 is to be build (IPRECA=4 and IMASS=1),
*                 add DCMASS*(Mass matrix) to the system matrix.
*
*   UPSAM  - control parameter of the streamline diffusion
*   RE     - 1/nu = viscosity
*   ISTOK  - =1 if the Stokes equation is to be discretized
*            rather than Navier-Stokes.
*   IMASS  - if BFULMT=TRUE:
*            =1: add the real mass matrix to the system matrix
*            =0: don't add the mass matrix to the system matrix
*
*   ELE    - used element for the discretisation
*   ICUBN  - Cubature formula to use for the discretization
*
* Out:
*   A      - system matrix; only if IDEF=0,1.
*            The nonlinearity is added to that matrix or
*            the linearized nonlinear system matrix is build
*            completely, respectively.
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
*   If BFULMT=true, IMASS=1, DCMASS<>0, D is updated according to
*
*             D = D - DCMASS*M*UUx - THSTEP * UUx * grad (Ux)
*   
*   If BFULMT=true, IMASS=1, DCMASS=0, D is updated according to
*
*             D = D + M*UUx - THSTEP * UUx * grad (Ux)
************************************************************************

      SUBROUTINE SUPAPX (NEQ,U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,
     *                   A,NA,KCOLA,KLDA,KVERT,KMID,DCORVG,TRIA,
     *                   ELE,ICUBN,
     *                   BFULMT, IMASS, ISTOK, UPSAM, RE,  
     *                   IDEF,DCMASS,THSTEP,
     *                   IALE,UMVEL)
                   
      IMPLICIT NONE
      
C include necessary COMMON blocks

      INCLUDE 'cout.inc' 
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      INCLUDE 'ccub.inc'
      
C parameters
      
      INTEGER NEQ
      DOUBLE PRECISION A(*)
      DOUBLE PRECISION U1L1(NEQ),U1L2(NEQ),U2L1(NEQ),U2L2(NEQ)
      DOUBLE PRECISION U1(NEQ),U2(NEQ),D1(NEQ),D2(NEQ)
      DOUBLE PRECISION DCORVG(2,*)
      DOUBLE PRECISION DENTRY(NNBAS,NNBAS)
      DOUBLE PRECISION A1L,A2L,DCMASS
      INTEGER NA
      INTEGER KCOLA(*),KLDA(*)
      INTEGER KVERT(NNVE,*),KMID(NNVE,*)
      INTEGER KENTRY(NNBAS,NNBAS)
      INTEGER IDEF
      INTEGER KDFG(NNBAS), KDFL(NNBAS), IDFL
      INTEGER KDFGLI(NNBAS), KDFLLI(NNBAS), IDFLLI, IELTLI
      
      INTEGER TRIA(SZTRIA),ICUBN
      
      INTEGER IALE
      DOUBLE PRECISION UMVEL(2,*)
      
      LOGICAL BFULMT
      
      INTEGER ISTOK, IMASS
      DOUBLE PRECISION UPSAM, RE
      
      DOUBLE PRECISION THSTEP
      
      EXTERNAL ELE

C externals

      INTEGER NDFL
      EXTERNAL NDFL

C local variables

      INTEGER ICUB,I,IELTYP,IEQ,JDOFE,ILD,JCOL,JCOL0,IDOFE,IDFG
      INTEGER IVE,JP,JDFL,JDFG,IDOFEH,JDOFEH,IA
      DOUBLE PRECISION DNY,CT0,DUMAX,DU1,DU2,DUNORM,DELTA,DUMAXR
      DOUBLE PRECISION DJ1,DJ2,DJ3,DJ4
      DOUBLE PRECISION XI1,XI2
      DOUBLE PRECISION OM,HBAS,HBASJ1,HBASJ2,HBASJ3,HSUMJ,AH
      DOUBLE PRECISION HBASI1,HBASI2,HBASI3,HSUMI,DENTH
      DOUBLE PRECISION DU1MV,DU2MV

C     In case that IPRECA=4, create the Laplace-part of the system
C     and include it into the system matrix each time this routine
C     is called. Otherwise don't create the Laplace part, just
C     add the nonlinearity.
C     This handling is simply realized by setting the factor NU in
C     front of the Laplace part to 0.

      IF (BFULMT) THEN
        DNY=1D0/RE
      ELSE
        DNY=0D0
      ENDIF

C     A similar handling holds for the case that the (full!) mass
C     matrix is to be included into the system matrix while the
C     matrix is to be rebuild completely. 
C     DCMASS represents the factor in front of the mass matrix.
C     Note that we cannot use DCMASS directly! In the Theta scheme,
C     we typically calculate 
C         [ DCMASS*M + THSTEP*(Laplace) + ... ]
C     so we have to weight everything except for the mass matrix!
C     We make a little trick here to realize that. We always weight
C     everything by THSTEP including the mass matrix - but divide
C     DCMASS by THSTEP before to compensate that!
C         THSTEP * [ CT0*M + (Laplace) + ... ]
C     with CT0 = DCMASS/THSTEP.
C
C     If only the nonlinear part is to be build (which is the normal
C     case, as rebuilding the comlete matrix including the full mass
C     matrix is normally not necessary), simply set CT0=0 which 
C     prevents calculating the mass matrix.

      IF (BFULMT.AND.(IMASS.EQ.1)) THEN
        CT0=DCMASS/THSTEP
      ELSE
        CT0=0D0
      ENDIF

C     Initialize BDER for our element. We want the element to calculate
C     function values as well as first X/Y derivatives:

      DO I = 1,NNDER
        BDER(I)=.FALSE.
      END DO

      DO I=1,3
        BDER(I)=.TRUE.
      END DO

C     Ask the element E011 about its type:

      IELTLI=-1
      CALL E011(0D0,0D0,IELTLI)
      
C     Get the local number of degrees of freedom of E011:
      
      IDFLLI=NDFL(IELTLI)
      
C     Ask the element about its type:

      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      
C     Get the local number of degrees of freedom:
      
      IDFL=NDFL(IELTYP)
      
C     Initialize the cubature formula identifier in the COMMON block
C     with our cubature formula we have to use here:

      ICUB=ICUBN
      CALL CB2Q(ICUB)
      
      IF (IER.NE.0) GOTO 99999

************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************

C     Initialize the values in the cubature points. Allowes the element
C     to precalculate some information for faster access later

      CALL ELE(0D0,0D0,-2)

C     Calculate the maximum norm of the actual velocity field
C     U = A1*U1 + A2*U2 into DUMAX. 
C     Round up the norm to 1D-8 if it's too small...

      DUMAX=0D0
      IF (A2L.EQ.0D0) THEN
        DO IEQ=1,NEQ
          DU1=U1L1(IEQ)
          DU2=U1L2(IEQ)
          DUNORM=SQRT(DU1**2+DU2**2)
          DUMAX=MAX(DUMAX,DUNORM)
        END DO
      ELSE       
        DO IEQ=1,NEQ
          DU1=A1L*U1L1(IEQ)+A2L*U2L1(IEQ)
          DU2=A1L*U1L2(IEQ)+A2L*U2L2(IEQ)
          DUNORM=SQRT(DU1**2+DU2**2)
          DUMAX=MAX(DUMAX,DUNORM)
        END DO
      ENDIF       

      IF (DUMAX.LT.1D-8) DUMAX=1D-8
      DUMAXR = 1D0/DUMAX

C *** Loop over all elements

      DO IEL=1,TRIA(ONEL)
      
C       Calculate the local and global degrees of freedom on the
C       current element IEL:
      
        CALL NDFGLX(TRIA,IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
        IF (IER.LT.0) GOTO 99999

C       If ALE is active, calculate the local and global DOF's
C       of the current element IEL with respect to E011:
        
        IF (IALE.NE.0) THEN
          CALL NDFGLX(TRIA,IEL,1,IELTLI,KVERT,KMID,KDFGLI,KDFLLI)
          IF (IER.LT.0) GOTO 99999
        END IF

C       Determine local DELTA for streamline-diffusion
C       (cf. p. 121 in Turek's CFD book).
C
C       For Stokes flow, we have the equation
C
C                -nu*div(u) + grad(p) = f
C       
C       not containing any convective part. Therefore we don't need
C       any stabilization technique. So in this case, switch of the
C       stabilization by setting DELTA to 0:

        IF (ISTOK.NE.1) THEN
          CALL DELSDX (U1L1,U1L2,U2L1,U2L2,A1L,A2L,IEL,DUMAXR,DELTA,
     *                 KVERT,KMID,DCORVG,KDFG,IDFL,UPSAM,RE)
        ELSE
          DELTA=0D0
        END IF

C       Determine entry positions in matrix.
C
C       Here we build a small matrix DENTRY/KENTRY which
C       corresponds to the current element. We will assemble the
C       contributions of our current element into this matrix and
C       will add the result into the main matrix later.
C
C       To successfully integrate the contributions into the main
C       matrix, we compute in advance the positions in the main
C       matrix where we have to add the contribution to.
C       KENTRY(X,Y) corresponds to the index in the array A
C       of the entry in line Y, row KCOL(line start + X).

        DO JDOFE=1,IDFL
          ILD=KLDA(KDFG(JDOFE))
          KENTRY(JDOFE,JDOFE)=ILD
          DENTRY(JDOFE,JDOFE)=0D0
          JCOL0=ILD
          DO IDOFE=1,IDFL
            IF (IDOFE.NE.JDOFE) THEN
              IDFG=KDFG(IDOFE)
              DO JCOL=JCOL0,NA
                IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
              END DO
113           JCOL0=JCOL+1
              KENTRY(JDOFE,IDOFE)=JCOL
              DENTRY(JDOFE,IDOFE)=0D0
            END IF
          END DO
        END DO

C       Calculate auxiliary Jacobian factors for the transformation 
C       onto the reference element. See QTRAFO.F for details...

        DO IVE = 1, TRIA(ONVE)
          JP=KVERT(IVE,IEL)
          KVE(IVE)=JP
          DX(IVE)=DCORVG(1,JP)
          DY(IVE)=DCORVG(2,JP)
        END DO

        DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
        DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
        DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
        DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))

C       Loop over the cubature points in our element to calculate
C       its contribution to each of the degrees of freedom on our
C       element

        DO ICUBP = 1, NCUBP
        
C         Get the coordinates of the cubature point on the reference
C         element
        
          XI1=DXI(ICUBP,1)
          XI2=DXI(ICUBP,2)

C         Calculate the Jacobian of the bilinear mapping onto the
C         reference element and the weight OM of the cubature formula:

          DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
          DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
          DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
          DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
          DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)
          OM=DOMEGA(ICUBP)*DETJ

C         If ALE is active, calculate the mesh velocity field in 
C         the cubature point.
C         The mesh velocity field is interpreted piecewise linear
C         using the element E011.
C         If ALE is not active, we simply set DU1MV and DU2MV to 0.
C
C         This must be done prior to the evaluation of the actual
C         element because otherwise it would overwrite the
C         DBAS entries!

          DU1MV = 0D0
          DU2MV = 0D0

          IF ((DCMASS.NE.0D0).AND.(IALE.NE.0)) THEN
          
            CALL E011(XI1,XI2,0)
            IF (IER.NE.0) GOTO 99999

            DO JDFL=1,IDFLLI
              HBAS=DBAS(KDFLLI(JDFL),1)
              JDFG=KDFGLI(JDFL)
              DU1MV=DU1MV+UMVEL(1,JDFG)*HBAS
              DU2MV=DU2MV+UMVEL(2,JDFG)*HBAS
            END DO
            
          END IF

C         Call the element to calculate the values in the current
C         cubature point on the reference element:

          CALL ELE(XI1,XI2,-3)
          IF (IER.LT.0) GOTO 99999

C         Now we have to assemble the "local" matrix DENTRY/KENTRY.
C         This assembling decomposes now into different parts,
C         depending on what has to me assembled.
C
C         We want to set up the nonlinear part of the matrix
C
C           n~_h (u_h, u_h, v_h) 
C
C         = n_h (u_h, u_h, v_h) + sum_T ( delta_T ( u_h*grad u_h, u_h*grad v_h)_T )
C           ^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C          standard nonlin. part                  stabilization
C
C         More precisely, as we want to assemble the matrix which is 
C         later multiplied with coefficient vectors, we have to insert
C         basis functions in the above terms instead of u_h and v_h.
C         Assuming the representation u_h=sum_j(u_j*Phi_j) and 
C         v_h=sum_i(u_i,Phi_i), the above term is evaluated in the
C         DOF's as:
C         
C           n_h (u_h, Phi_j, Phi_i) 
C         + sum_T ( delta_T ( u_h*grad Phi_j, u_h*grad Phi_i )_T )
C
C         In nonstationary simulations, the system matrix typically
C         contains a mass matrix to respect the time derivative.
C         The matrix has the form
C
C         [  dcmass*M*I  +  THWEIG * (-nu * Laplace(.))  ] + THWEIG * u grad(.)
C
C         so if DCMASS<>0, incorporate the (real, not lumped!) mass matrix
C         into the local matrix.

          IF (DCMASS.NE.0D0) THEN

C           Calculate the actual velocity in the current cubature point
C           into (DU1,DU2). If we only have a primary velocity field
C           (A2L=0), we can calculate that only by summing up the
C           velocities in U1Lx, otherwise we have to sum up
C           A1*U1Lx + A2*U2Lx.

            DU1=0D0
            DU2=0D0
            IF (A2L.EQ.0D0) THEN
              DO JDFL=1,IDFL
                HBAS=DBAS(KDFL(JDFL),1)
                IF (ABS(HBAS).GE.1D-8) THEN
                  JDFG=KDFG(JDFL)
                  DU1=DU1+U1L1(JDFG)*HBAS
                  DU2=DU2+U1L2(JDFG)*HBAS
                ENDIF
              END DO
            ELSE
              DO JDFL=1,IDFL
                HBAS=DBAS(KDFL(JDFL),1)
                IF (ABS(HBAS).GE.1D-8) THEN
                  JDFG=KDFG(JDFL)
                  DU1=DU1+(A1L*U1L1(JDFG)+A2L*U2L1(JDFG))*HBAS
                  DU2=DU2+(A1L*U1L2(JDFG)+A2L*U2L2(JDFG))*HBAS
                ENDIF
              END DO
            ENDIF

C           We take a more detailed look onto the last scalar product
C           of n~_h (u_h, u_h, v_h) what we want to calculate here.
C
C           The vector u_h=(DU1,DU2) contains both velocity components,
C           for the X as well as for the Y velocity. On the other hand
C           the system matrix we want to build here will be designed for 
C           one velocity component only! Therefore, Phi_i and Phi_j
C           are scalar functions, so grad(Phi_i), grad(Phi_j) are vectors
C           with two components. Therefore, the last scalar product is more 
C           in detail:
C
C               ( u_h*grad Phi_j, u_h*grad Phi_i )_T
C
C           =   ( < (DU1) , (grad(Phi_j)_1) > , < (DU1) , (grad(Phi_i)_1) > )_T
C                   (DU2) , (grad(Phi_j)_2)       (DU2) , (grad(Phi_i)_2)  
C
C           =   < (DU1) , (grad(Phi_j)_1) >  *  < (DU1) , (grad(Phi_j)_1) >
C                 (DU2) , (grad(Phi_j)_2)         (DU2) , (grad(Phi_j)_2)
C
C           =   HSUMJ * HSUMI
C
C           i.e. a product of two scalar values!
C
C           Summing up over all pairs of multiindices.
C
C           Outer loop over the DOF's j=1..IDFL on our current element, 
C           which corresponds to the basis functions Phi_j:

            DO JDOFE=1,IDFL
            
C             Fetch the contributions of the basis functions Phi_j for
C             function value and first derivatives for the current
C             DOF into HBASxy:
            
              JDOFEH=KDFL(JDOFE)
              HBASJ1=DBAS(JDOFEH,1)
              HBASJ2=DBAS(JDOFEH,2)
              HBASJ3=DBAS(JDOFEH,3)

C             If ALE is not active, calculate 
C
C                 U * grad(Phi_j)  =  < grad(Phi_j), U >
C     
C               = ( grad(Phi_j)_1 , (DU1) )
C                 ( grad(Phi_j)_2   (DU2) )
C
C             Remember: DU1MV=DU2MV=0 in this case.
C             If ALE is active, use v=mesh velocity and calculate 
C
C                   (U-v) * grad(Phi_j)  =  < grad(Phi_j), U-v >
C     
C                 = ( grad(Phi_j)_1 , (DU1-v) )
C                   ( grad(Phi_j)_2   (DU2-v) )
              
              HSUMJ = HBASJ2*(DU1-DU1MV) + HBASJ3*(DU2-DU2MV)

C             Inner loop over the DOF's i=1..IDFL, which corresponds to
C             the basis function Phi_i:

              DO IDOFE=1,IDFL
              
                IF (IDOFE.EQ.JDOFE) THEN
                
C                 Short version of the evaluation of the matrix
C                 contribution - see below for a more detailed
C                 description what is added together here!
                
                  AH = HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *               + DNY*(HBASJ2**2+HBASJ3**2)
     *               + CT0*HBASJ1**2
     
                ELSE
                
C                 Fetch the contributions of the basis function Phi_i for
C                 function value and first derivatives for the current
C                 DOF into HBASIy:
                
                  IDOFEH=KDFL(IDOFE)
                  HBASI1=DBAS(IDOFEH,1)
                  HBASI2=DBAS(IDOFEH,2)
                  HBASI3=DBAS(IDOFEH,3)

C                 Calculate 
C
C                     U * grad(Phi_i)  =  < grad(Phi_i), U >
C     
C                   = ( grad(Phi_i)_1 , (DU1) )
C                     ( grad(Phi_i)_2   (DU2) )
C
C                 Remember: DU1MV=DU2MV=0 in this case.
C
C                 If ALE is active, use v=mesh velocity and calculate 
C
C                     (U-v) * grad(Phi_i)  =  < grad(Phi_i), U-v >
C     
C                   = ( grad(Phi_i)_1 , (DU1-DU1MV) )
C                     ( grad(Phi_i)_2   (DU2-DU2MV) )

                  HSUMI=HBASI2*(DU1-DU1MV)+HBASI3*(DU2-DU2MV)
     
C                 Finally calculate the contribution to the system
C                 matrix. Depending on the configuration of DNY,
C                 IPRECA, DCMASS,... this decomposes into three
C                 different parts:
C
C                 AH = n~_h(u_h,phi_j,phi_i)        | nonlinear part
C                    + DNY*(grad(phi_j,grad(phi_i)) | -nu*Laplace(u)
C                    + CT0*(phi_j*phi_i)            | Mass matrix
C
C                 The last two parts are probably not added to the
C                 matrix by setting DNY or CT0 to 0, respectively.
C
C                 For saving some numerical operations, we write:
C     
C                     HSUMI * (Delta * HSUMJ + HBASJ1)
C
C                 =   Delta * HSUMI * HSUMJ  
C                   + HSUMI * HBASJ1
C     
C                 =   Delta * ( U*grad(Phi_j), U*grad(Phi_i) )
C                   + (grad(Phi_i)*U,Phi_j)
C
C               <->   sum_T ( delta_T ( u_h*grad Phi_j, u_h*grad Phi_i) )_T
C                   + n_h (u_h, Phi_j, Phi_i)
                  
                  AH = HSUMI*(DELTA*HSUMJ+HBASJ1)
     *               + DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *               + CT0*HBASJ1*HBASI1
     
                ENDIF ! (IDOFE.EQ.JDOFE)

C               Weighten the calculated value AH by the cubature
C               weight OM and add it to the local matrix. After the
C               loop over all DOF's is finished, each entry contains
C               the calculated integral.

                DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+OM*AH
                
              END DO ! IDOFE
              
            END DO ! JDOFE

          ELSE

C           Coefficient in front of the mass matrix is 0.
C           Subtract the mass matrix from the system matrix.
C
C           Outer loop over the DOF's j=1..IDFL corresponding
C           to Phi_j:

            DO JDOFE=1,IDFL
            
C             Fetch the function value in the current DOF:
            
              HBASJ1=DBAS(KDFL(JDOFE),1)

C             Inner loop over the DOF's i=1..IDFL corresponding
C             to Phi_i:

              DO IDOFE=1,IDFL
                
C               Fetch the function value of the other DOF of the
C               current element
                
                HBASI1=DBAS(KDFL(IDOFE),1)
                
C               Calculate the contribution for the entry. The factor
C               THSTEP is compensated later when the local matrix
C               is included into the global matrix and/or in the
C               modification of the RHS vector.
                
                AH=-1D0/THSTEP*HBASJ1*HBASI1
                
C               ...and incorporate it into the local matrix.
                
                DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+OM*AH
                
              END DO ! IDOFE
              
            END DO ! JDOFE

          END IF ! (DCMASS.NE.0D0)

        END DO ! ICUBP
      
C       Now we have set up a "local" system matrix. We can either
C       include it into the real matrix or we can use it to simply
C       modify the RHS vector to create a defect vector (throwing
C       away the information about the matrix afterwards, which would
C       result in a matrix free modification of the RHS vector).
      
        DO JDOFE=1,IDFL
        
          DO IDOFE=1,IDFL
          
C           Get the entry from the local matrix and weight it according
C           to the current THETA in the Theta scheme, given by the
C           parameter. 
C           (Remark: For stationary simulations, THSTEP is typically
C            1D0 which includes the local matrix into the global one
C            directly).
          
            DENTH=THSTEP*DENTRY(JDOFE,IDOFE)

C           For IDEF=0/1, incorporate our "local" system matrix into 
C           the global matrix. The position of each entry DENTRY(X,Y) 
C           in the global matrix array A was saved in element KENTRY(X,Y)
C           before.

            IF (IDEF.LT.2) THEN
              IA   =KENTRY(JDOFE,IDOFE)
              A(IA)=A(IA)+DENTH
            ENDIF

C           For IDEF=1,2, build the defect vector
C               D = RHS - A*U
C           This is done matrix free, only with the help of the local 
C           matrix.
C           In this case, D=(D1,D2) is expected to be the RHS on
C           entry and will be updated to be the defect vector when
C           this routine is left.

            IF (IDEF.GT.0) THEN 
              IDFG=KDFG(IDOFE)
              JDFG=KDFG(JDOFE)
              D1(JDFG)= D1(JDFG)-DENTH*U1(IDFG)
              D2(JDFG)= D2(JDFG)-DENTH*U2(IDFG)
            ENDIF 

          END DO ! IDOFE
          
        END DO ! JDOFE

C       Matrix/defect vector updated for that element; proceed to the
C       next one...

      END DO ! IEL

99999 END

************************************************************************
* Streamline diffusion, ALE-method
*
* Adds the SUPG-part on matrix block A after it was initialized by the
* linear part or builds up the complete nonlinear system matrix.
* The input vector Ui is the old velocity field.
* The input vectors UjLi are the transport directions.
*
* Nonparametric version
*
* In contrast to SUPGNX, this routine supports the 1st order ALE
* method. The caller has to specify a grid velocity vector field,
* which gives the velocity of the grid points in each time step.
* This vector field is interpreted piecewise linear with the element
* E011 and subtracted from the main velocity field to calculate
* the nonlinearity (cf. [Duarte, Formaz, Natesan; "Arbitrary Lagrangian-
* Euler Method for Navier-Stokes equations with moving boundaries";
* Comput. Methods Appl. Mech. Engrg. 193 (2004), 4819-4836).
*
* Extended calling convention
*
* In:
*   A, NA,
*   KCOLA,
*   LLDA   - array [1..*] of double/integer
*            Structure arrays of system matrix, maybe initialised
*            with linear parts of the nonlinear system matrix.
*   KVERT,
*   KMID,
*   DCORVG - usual geometry information; must correspond to TRIA
*   
*   TRIA   - array [1..SZTRIA] of integer
*            Triangulation structure of the underlying mesh
*
*   NEQ    - length of solution/defect vectors
*
*   U1L1,
*   U1L2   - array [1..NEQ] of double
*            main velocity field used for assembling of 
*            the nonlinearity. Can be undefined if A1L=0.
*   U2L1,
*   U2L2   - array [1..NEQ] of double
*            secondary velocity field, used for the assembling
*            of the nonlinearity. Can be undefined if A2L=0.
*   A1L    - double; weighting factor for U1L1/U1L2
*   A2L    - double; weighting factor for U2L1/U2L2
*
*   IALE   - =0: don't use ALE, U1MVEL/U2MVEl can be undefined
*            =1: use ALE method; U1MVEL/U2MVEl must be undefined
*   UMVEL  - array [1..2,1..NVT] of double
*            Grid velocity field. UMVEL(1,.) is the X-velocity,
*            UMVEL(2,.) the Y-velocity in each corner grid point.
*
*   IDEF   - Controls the behaviour of this routine.
*            =0: modify system matrix, add nonlinearity.
*                Defect vectors D1,D2 and velocity vectors U1,U2
*                can be undefined.
*            =1: modify both, system matrix and defect vector
*            =2: modify defect vectors, include nonlinearity.
*                A can be undefined (but not KCOLA,KLDA!).
*                Can be used to modify the defect vector matrix
*                free. Note that the rules how to set up the
*                system matrix (NY, DCMASS,...) still hold,
*                although the matrix itself is not modified!
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
*   BFULMT - Whether or not to generate the full nonlinear system
*            matrix.
*            =false: generate only the nonlinear part 
*                         THWEIG * u grad(.)
*                    of the system matrix and add it to the existing
*                    matrix in A
*            =true : Generate the full system matrix
*               [  M*I  +  THWEIG * (-nu * Laplace(.) * u * grad(.)) ] 
*                    and add it to A
*            
*   THSTEP - weighting factor of the nonlinearity N(u) which is to be
*            added to the matrix
*            (is used as step-size of a Theta-scheme)
*
*   DCMASS - This defines how to modify the system matrix - if the
*            matrix is modified at all (see IDEF!).
*             =0: subtract the mass matrix from the system matrix;
*                 The nonlinearity is not build.
*            <>0: Add the nonlinear term including the stabilization
*                 to the system matrix.
*                 If the full matrix with not-lumped mass matrix 
*                 is to be build (IPRECA=4 and IMASS=1),
*                 add DCMASS*(Mass matrix) to the system matrix.
*
*   UPSAM  - control parameter of the streamline diffusion
*   RE     - 1/nu = viscosity
*   ISTOK  - =1 if the Stokes equation is to be discretized
*            rather than Navier-Stokes.
*   IMASS  - if BFULMT=TRUE:
*            =1: add the real mass matrix to the system matrix
*            =0: don't add the mass matrix to the system matrix
*
*   ELE    - used element for the discretisation
*   ICUBN  - Cubature formula to use for the discretization
*
* Out:
*   A      - system matrix; only if IDEF=0,1.
*            The nonlinearity is added to that matrix or
*            the linearized nonlinear system matrix is build
*            completely, respectively.
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
*   If BFULMT=true, IMASS=1, DCMASS<>0, D is updated according to
*
*             D = D - DCMASS*M*UUx - THSTEP * UUx * grad (Ux)
*   
*   If BFULMT=true, IMASS=1, DCMASS=0, D is updated according to
*
*             D = D + M*UUx - THSTEP * UUx * grad (Ux)
************************************************************************

      SUBROUTINE SUPANX (NEQ,U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,
     *                   A,NA,KCOLA,KLDA,KVERT,KMID,DCORVG,TRIA,
     *                   ELE,ICUBN,
     *                   BFULMT, IMASS, ISTOK, UPSAM, RE,  
     *                   IDEF,DCMASS,THSTEP,
     *                   IALE,UMVEL)

      IMPLICIT NONE
      
C include necessary COMMON blocks

      INCLUDE 'cout.inc' 
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      INCLUDE 'ccub.inc'
      
C parameters
      
      INTEGER NEQ
      DOUBLE PRECISION A(*)
      DOUBLE PRECISION U1L1(NEQ),U1L2(NEQ),U2L1(NEQ),U2L2(NEQ)
      DOUBLE PRECISION U1(NEQ),U2(NEQ),D1(NEQ),D2(NEQ)
      DOUBLE PRECISION DCORVG(2,*)
      DOUBLE PRECISION DENTRY(NNBAS,NNBAS)
      DOUBLE PRECISION A1L,A2L,DCMASS
      INTEGER NA
      INTEGER KCOLA(*),KLDA(*)
      INTEGER KVERT(NNVE,*),KMID(NNVE,*)
      INTEGER KENTRY(NNBAS,NNBAS)
      INTEGER IDEF
      INTEGER KDFG(NNBAS), KDFL(NNBAS), IDFL
      INTEGER KDFGLI(NNBAS), KDFLLI(NNBAS), IDFLLI, IELTLI

      INTEGER IALE
      DOUBLE PRECISION UMVEL(2,*)
      
      INTEGER TRIA(SZTRIA),ICUBN
      
      LOGICAL BFULMT
      
      INTEGER ISTOK, IMASS
      DOUBLE PRECISION UPSAM, RE
      
      DOUBLE PRECISION THSTEP
      
      EXTERNAL ELE

C externals

      INTEGER NDFL
      EXTERNAL NDFL
      EXTERNAL E011

C local variables

      INTEGER ICUB,I,IELTYP,IEQ,JDOFE,ILD,JCOL,JCOL0,IDOFE,IDFG
      INTEGER IVE,JP,JDFL,JDFG,IDOFEH,JDOFEH,IA
      DOUBLE PRECISION DNY,CT0,DUMAX,DU1,DU2,DUNORM,DELTA,DUMAXR
      DOUBLE PRECISION DJ1,DJ2,DJ3,DJ4
      DOUBLE PRECISION XI1,XI2,XX,YY
      DOUBLE PRECISION OM,HBAS,HBASJ1,HBASJ2,HBASJ3,HSUMJ,AH
      DOUBLE PRECISION HBASI1,HBASI2,HBASI3,HSUMI,DENTH
      DOUBLE PRECISION DU1MV,DU2MV

C     In case that IPRECA=4, create the Laplace-part of the system
C     and include it into the system matrix each time this routine
C     is called. Otherwise don't create the Laplace part, just
C     add the nonlinearity.
C     This handling is simply realized by setting the factor NU in
C     front of the Laplace part to 0.
      
      IF (BFULMT) THEN
        DNY=1D0/RE
      ELSE
        DNY=0D0
      ENDIF

C     A similar handling holds for the case that the (full!) mass
C     matrix is to be included into the system matrix while the
C     matrix is to be rebuild completely. 
C     DCMASS represents the factor in front of the mass matrix.
C     Note that we cannot use DCMASS directly! In the Theta scheme,
C     we typically calculate 
C         [ DCMASS*M + THSTEP*(Laplace) + ... ]
C     so we have to weight everything except for the mass matrix!
C     We make a little trick here to realize that. We always weight
C     everything by THSTEP including the mass matrix - but divide
C     DCMASS by THSTEP before to compensate that!
C         THSTEP * [ CT0*M + (Laplace) + ... ]
C     with CT0 = DCMASS/THSTEP.
C
C     If only the nonlinear part is to be build (which is the normal
C     case, as rebuilding the comlete matrix including the full mass
C     matrix is normally not necessary), simply set CT0=0 which 
C     prevents calculating the mass matrix.

      IF (BFULMT.AND.(IMASS.EQ.1)) THEN
        CT0=DCMASS/THSTEP
      ELSE
        CT0=0D0
      ENDIF

C     Initialize BDER for our element. We want the element to calculate
C     function values as well as first X/Y derivatives:

      DO I = 1,NNDER
        BDER(I)=.FALSE.
      END DO

      DO I=1,3
        BDER(I)=.TRUE.
      END DO

C     Ask the element E011 about its type:

      IELTLI=-1
      CALL E011(0D0,0D0,IELTLI)
      
C     Get the local number of degrees of freedom of E011:
      
      IDFLLI=NDFL(IELTLI)
      
C     Ask the element about its type:

      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      
C     Get the local number of degrees of freedom:
      
      IDFL=NDFL(IELTYP)
      
C     Initialize the cubature formula identifier in the COMMON block
C     with our cubature formula we have to use here:

      ICUB=ICUBN
      CALL CB2Q(ICUB)
      
      IF (IER.NE.0) GOTO 99999

************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************

C     Calculate the maximum norm of the actual velocity field
C     U = A1*U1 + A2*U2 into DUMAX. 
C     Round up the norm to 1D-8 if it's too small...

      DUMAX=0D0
      IF (A2L.EQ.0D0) THEN
        DO IEQ=1,NEQ
          DU1=U1L1(IEQ)
          DU2=U1L2(IEQ)
          DUNORM=SQRT(DU1**2+DU2**2)
          DUMAX=MAX(DUMAX,DUNORM)
        END DO
      ELSE       
        DO IEQ=1,NEQ
          DU1=A1L*U1L1(IEQ)+A2L*U2L1(IEQ)
          DU2=A1L*U1L2(IEQ)+A2L*U2L2(IEQ)
          DUNORM=SQRT(DU1**2+DU2**2)
          DUMAX=MAX(DUMAX,DUNORM)
        END DO
      ENDIF       

      IF (DUMAX.LT.1D-8) DUMAX=1D-8
      DUMAXR = 1D0/DUMAX

C *** Loop over all elements

      DO IEL=1,TRIA(ONEL)
      
C       Calculate the local and global degrees of freedom on the
C       current element IEL:
      
        CALL NDFGLX(TRIA,IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
        IF (IER.LT.0) GOTO 99999

C       If ALE is active, calculate the local and global DOF's
C       of the current element IEL with respect to E011:
        
        IF (IALE.NE.0) THEN
          CALL NDFGLX(TRIA,IEL,1,IELTLI,KVERT,KMID,KDFGLI,KDFLLI)
          IF (IER.LT.0) GOTO 99999
        END IF

C       Determine local DELTA for streamline-diffusion
C       (cf. p. 121 in Turek's CFD book).
C
C       For Stokes flow, we have the equation
C
C                -nu*div(u) + grad(p) = f
C       
C       not containing any convective part. Therefore we don't need
C       any stabilization technique. So in this case, switch of the
C       stabilization by setting DELTA to 0:

        IF (ISTOK.NE.1) THEN
          CALL DELSDX (U1L1,U1L2,U2L1,U2L2,A1L,A2L,IEL,DUMAXR,DELTA,
     *                 KVERT,KMID,DCORVG,KDFG,IDFL,UPSAM,RE)
        ELSE
          DELTA=0D0
        END IF

C       Determine entry positions in matrix.
C
C       Here we build a small matrix DENTRY/KENTRY which
C       corresponds to the current element. We will assemble the
C       contributions of our current element into this matrix and
C       will add the result into the main matrix later.
C
C       To successfully integrate the contributions into the main
C       matrix, we compute in advance the positions in the main
C       matrix where we have to add the contribution to.
C       KENTRY(X,Y) corresponds to the index in the array A
C       of the entry in line Y, row KCOL(line start + X).

        DO JDOFE=1,IDFL
          ILD=KLDA(KDFG(JDOFE))
          KENTRY(JDOFE,JDOFE)=ILD
          DENTRY(JDOFE,JDOFE)=0D0
          JCOL0=ILD
          DO IDOFE=1,IDFL
            IF (IDOFE.NE.JDOFE) THEN
              IDFG=KDFG(IDOFE)
              DO JCOL=JCOL0,NA
                IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
              END DO
113           JCOL0=JCOL+1
              KENTRY(JDOFE,IDOFE)=JCOL
              DENTRY(JDOFE,IDOFE)=0D0
            END IF
          END DO
        END DO

C       Calculate auxiliary Jacobian factors for the transformation 
C       onto the reference element. See QTRAFO.F for details...

        DO IVE = 1, TRIA(ONVE)
          JP=KVERT(IVE,IEL)
          KVE(IVE)=JP
          DX(IVE)=DCORVG(1,JP)
          DY(IVE)=DCORVG(2,JP)
        END DO

        DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
        DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
        DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
        DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))

C       Initialize the values in the cubature points. Allowes the element
C       to precalculate some information for faster access later

        CALL ELE(0D0,0D0,-2)
        IF (IER.NE.0) GOTO 99999

C       Loop over the cubature points in our element to calculate
C       its contribution to each of the degrees of freedom on our
C       element

        DO ICUBP = 1, NCUBP
        
C         Get the coordinates of the cubature point on the reference
C         element
        
          XI1=DXI(ICUBP,1)
          XI2=DXI(ICUBP,2)

C         Calculate the Jacobian of the bilinear mapping onto the
C         reference element and the weight OM of the cubature formula:

          DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
          DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
          DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
          DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
          DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)
          OM=DOMEGA(ICUBP)*DETJ
          
C         Calculate the real coordinates of the cubature point
C         with the bilinear transformation

          XX=0.5D0*(DX(1)+DX(2)+DJ1)+0.5D0*(DX(2)-DX(1)+DJ2)*XI1
     *      +0.5D0*DJ1*XI2+0.5D0*DJ2*XI1*XI2
          YY=0.5D0*(DY(1)+DY(3)+DJ3)+0.5D0*DJ4*XI1
     *      +0.5D0*(DY(3)-DY(1)-DJ4)*XI2-0.5D0*DJ3*XI1*XI2

C         If ALE is active, calculate the mesh velocity field in 
C         the cubature point.
C         The mesh velocity field is interpreted piecewise linear
C         using the element E011.
C         If ALE is not active, we simply set DU1MV and DU2MV to 0.
C
C         This must be done prior to the evaluation of the actual
C         element because otherwise it would overwrite the
C         DBAS entries!

          DU1MV = 0D0
          DU2MV = 0D0

          IF ((DCMASS.NE.0D0).AND.(IALE.NE.0)) THEN
          
            CALL E011(XI1,XI2,0)
            IF (IER.NE.0) GOTO 99999

            DO JDFL=1,IDFLLI
              HBAS=DBAS(KDFLLI(JDFL),1)
              JDFG=KDFGLI(JDFL)
              DU1MV=DU1MV+UMVEL(1,JDFG)*HBAS
              DU2MV=DU2MV+UMVEL(2,JDFG)*HBAS
            END DO
            
          END IF

C         Call the element to calculate the values in the current
C         cubature point on the reference element:

          CALL ELE(XX,YY,-3)
          IF (IER.LT.0) GOTO 99999

C         Now we have to assemble the "local" matrix DENTRY/KENTRY.
C         This assembling decomposes now into different parts,
C         depending on what has to me assembled.
C
C         We want to set up the nonlinear part of the matrix
C
C           n~_h (u_h, u_h, v_h) 
C
C         = n_h (u_h, u_h, v_h) + sum_T ( delta_T ( u_h*grad u_h, u_h*grad v_h)_T )
C           ^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C          standard nonlin. part                  stabilization
C
C         More precisely, as we want to assemble the matrix which is 
C         later multiplied with coefficient vectors, we have to insert
C         basis functions in the above terms instead of u_h and v_h.
C         Assuming the representation u_h=sum_j(u_j*Phi_j) and 
C         v_h=sum_i(u_i,Phi_i), the above term is evaluated in the
C         DOF's as:
C         
C           n_h (u_h, Phi_j, Phi_i) 
C         + sum_T ( delta_T ( u_h*grad Phi_j, u_h*grad Phi_i )_T )
C
C         In nonstationary simulations, the system matrix typically
C         contains a mass matrix to respect the time derivative.
C         The matrix has the form
C
C         [  dcmass*M*I  +  THWEIG * (-nu * Laplace(.))  ] + THWEIG * u grad(.)
C
C         so if DCMASS<>0, incorporate the (real, not lumped!) mass matrix
C         into the local matrix.

          IF (DCMASS.NE.0D0) THEN

C           Calculate the actual velocity in the current cubature point
C           into (DU1,DU2). If we only have a primary velocity field
C           (A2L=0), we can calculate that only by summing up the
C           velocities in U1Lx, otherwise we have to sum up
C           A1*U1Lx + A2*U2Lx.

            DU1=0D0
            DU2=0D0
            IF (A2L.EQ.0D0) THEN
              DO JDFL=1,IDFL
                HBAS=DBAS(KDFL(JDFL),1)
                IF (ABS(HBAS).GE.1D-8) THEN
                  JDFG=KDFG(JDFL)
                  DU1=DU1+U1L1(JDFG)*HBAS
                  DU2=DU2+U1L2(JDFG)*HBAS
                ENDIF
              END DO
            ELSE
              DO JDFL=1,IDFL
                HBAS=DBAS(KDFL(JDFL),1)
                IF (ABS(HBAS).GE.1D-8) THEN
                  JDFG=KDFG(JDFL)
                  DU1=DU1+(A1L*U1L1(JDFG)+A2L*U2L1(JDFG))*HBAS
                  DU2=DU2+(A1L*U1L2(JDFG)+A2L*U2L2(JDFG))*HBAS
                ENDIF
              END DO
            ENDIF
            
C           We take a more detailed look onto the last scalar product
C           of n~_h (u_h, u_h, v_h) what we want to calculate here.
C
C           The vector u_h=(DU1,DU2) contains both velocity components,
C           for the X as well as for the Y velocity. On the other hand
C           the system matrix we want to build here will be designed for 
C           one velocity component only! Therefore, Phi_i and Phi_j
C           are scalar functions, so grad(Phi_i), grad(Phi_j) are vectors
C           with two components. Therefore, the last scalar product is more 
C           in detail:
C
C               ( u_h*grad Phi_j, u_h*grad Phi_i )_T
C
C           =   ( < (DU1) , (grad(Phi_j)_1) > , < (DU1) , (grad(Phi_i)_1) > )_T
C                   (DU2) , (grad(Phi_j)_2)       (DU2) , (grad(Phi_i)_2)  
C
C           =   < (DU1) , (grad(Phi_j)_1) >  *  < (DU1) , (grad(Phi_j)_1) >
C                 (DU2) , (grad(Phi_j)_2)         (DU2) , (grad(Phi_j)_2)
C
C           =   HSUMJ * HSUMI
C
C           i.e. a product of two scalar values!
C
C           Summing up over all pairs of multiindices.
C
C           Outer loop over the DOF's j=1..IDFL on our current element, 
C           which corresponds to the basis functions Phi_j:

            DO JDOFE=1,IDFL
            
C             Fetch the contributions of the basis functions Phi_j for
C             function value and first derivatives for the current
C             DOF into HBASxy:
            
              JDOFEH=KDFL(JDOFE)
              HBASJ1=DBAS(JDOFEH,1)
              HBASJ2=DBAS(JDOFEH,2)
              HBASJ3=DBAS(JDOFEH,3)
              
C             If ALE is not active, calculate 
C
C                 U * grad(Phi_j)  =  < grad(Phi_j), U >
C     
C               = ( grad(Phi_j)_1 , (DU1) )
C                 ( grad(Phi_j)_2   (DU2) )
C
C             Remember: DU1MV=DU2MV=0 in this case.
C             If ALE is active, use v=mesh velocity and calculate 
C
C                   (U-v) * grad(Phi_j)  =  < grad(Phi_j), U-v >
C     
C                 = ( grad(Phi_j)_1 , (DU1-v) )
C                   ( grad(Phi_j)_2   (DU2-v) )
              
              HSUMJ = HBASJ2*(DU1-DU1MV) + HBASJ3*(DU2-DU2MV)

C             Inner loop over the DOF's i=1..IDFL, which corresponds to
C             the basis function Phi_i:

              DO IDOFE=1,IDFL
              
                IF (IDOFE.EQ.JDOFE) THEN
                
C                 Short version of the evaluation of the matrix
C                 contribution - see below for a more detailed
C                 description what is added together here!
                
                  AH = HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *               + DNY*(HBASJ2**2+HBASJ3**2)
     *               + CT0*HBASJ1**2
     
                ELSE
                
C                 Fetch the contributions of the basis function Phi_i for
C                 function value and first derivatives for the current
C                 DOF into HBASIy:
                
                  IDOFEH=KDFL(IDOFE)
                  HBASI1=DBAS(IDOFEH,1)
                  HBASI2=DBAS(IDOFEH,2)
                  HBASI3=DBAS(IDOFEH,3)

C                 Calculate 
C
C                     U * grad(Phi_i)  =  < grad(Phi_i), U >
C     
C                   = ( grad(Phi_i)_1 , (DU1) )
C                     ( grad(Phi_i)_2   (DU2) )
C
C                 Remember: DU1MV=DU2MV=0 in this case.
C
C                 If ALE is active, use v=mesh velocity and calculate 
C
C                     (U-v) * grad(Phi_i)  =  < grad(Phi_i), U-v >
C     
C                   = ( grad(Phi_i)_1 , (DU1-DU1MV) )
C                     ( grad(Phi_i)_2   (DU2-DU2MV) )

                  HSUMI=HBASI2*(DU1-DU1MV)+HBASI3*(DU2-DU2MV)
     
C                 Finally calculate the contribution to the system
C                 matrix. Depending on the configuration of DNY,
C                 IPRECA, DCMASS,... this decomposes into three
C                 different parts:
C
C                 AH = n~_h(u_h,phi_j,phi_i)        | nonlinear part
C                    + DNY*(grad(phi_j,grad(phi_i)) | -nu*Laplace(u)
C                    + CT0*(phi_j*phi_i)            | Mass matrix
C
C                 The last two parts are probably not added to the
C                 matrix by setting DNY or CT0 to 0, respectively.
C
C                 For saving some numerical operations, we write:
C     
C                     HSUMI * (Delta * HSUMJ + HBASJ1)
C
C                 =   Delta * HSUMI * HSUMJ  
C                   + HSUMI * HBASJ1
C     
C                 =   Delta * ( U*grad(Phi_j), U*grad(Phi_i) )
C                   + (grad(Phi_i)*U,Phi_j)
C
C               <->   sum_T ( delta_T ( u_h*grad Phi_j, u_h*grad Phi_i) )_T
C                   + n_h (u_h, Phi_j, Phi_i)
                  
                  AH = HSUMI*(DELTA*HSUMJ+HBASJ1)
     *               + DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *               + CT0*HBASJ1*HBASI1
     
                ENDIF ! (IDOFE.EQ.JDOFE)

C               Weighten the calculated value AH by the cubature
C               weight OM and add it to the local matrix. After the
C               loop over all DOF's is finished, each entry contains
C               the calculated integral.

                DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+OM*AH
                
              END DO ! IDOFE
              
            END DO ! JDOFE

          ELSE

C           Coefficient in front of the mass matrix is 0.
C           Subtract the mass matrix from the system matrix.
C
C           Outer loop over the DOF's j=1..IDFL corresponding
C           to Phi_j:

            DO JDOFE=1,IDFL
            
C             Fetch the function value in the current DOF:
            
              HBASJ1=DBAS(KDFL(JDOFE),1)

C             Inner loop over the DOF's i=1..IDFL corresponding
C             to Phi_i:

              DO IDOFE=1,IDFL
                
C               Fetch the function value of the other DOF of the
C               current element
                
                HBASI1=DBAS(KDFL(IDOFE),1)
                
C               Calculate the contribution for the entry. The factor
C               THSTEP is compensated later when the local matrix
C               is included into the global matrix and/or in the
C               modification of the RHS vector.
                
                AH=-1D0/THSTEP*HBASJ1*HBASI1
                
C               ...and incorporate it into the local matrix.
                
                DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+OM*AH
                
              END DO ! IDOFE
              
            END DO ! JDOFE

          END IF ! (DCMASS.NE.0D0)

        END DO ! ICUBP
      
C       Now we have set up a "local" system matrix. We can either
C       include it into the real matrix or we can use it to simply
C       modify the RHS vector to create a defect vector (throwing
C       away the information about the matrix afterwards, which would
C       result in a matrix free modification of the RHS vector).
      
        DO JDOFE=1,IDFL
        
          DO IDOFE=1,IDFL
          
C           Get the entry from the local matrix and weight it according
C           to the current THETA in the Theta scheme, given by the
C           parameter. 
C           (Remark: For stationary simulations, THSTEP is typically
C            1D0 which includes the local matrix into the global one
C            directly).
          
            DENTH=THSTEP*DENTRY(JDOFE,IDOFE)

C           For IDEF=0/1, incorporate our "local" system matrix into 
C           the global matrix. The position of each entry DENTRY(X,Y) 
C           in the global matrix array A was saved in element KENTRY(X,Y)
C           before.

            IF (IDEF.LT.2) THEN
              IA   =KENTRY(JDOFE,IDOFE)
              A(IA)=A(IA)+DENTH
            ENDIF

C           For IDEF=1,2, build the defect vector
C               D = RHS - A*U
C           This is done matrix free, only with the help of the local 
C           matrix.
C           In this case, D=(D1,D2) is expected to be the RHS on
C           entry and will be updated to be the defect vector when
C           this routine is left.

            IF (IDEF.GT.0) THEN 
              IDFG=KDFG(IDOFE)
              JDFG=KDFG(JDOFE)
              D1(JDFG)= D1(JDFG)-DENTH*U1(IDFG)
              D2(JDFG)= D2(JDFG)-DENTH*U2(IDFG)
            ENDIF 

          END DO ! IDOFE
          
        END DO ! JDOFE

C       Matrix/defect vector updated for that element; proceed to the
C       next one...

      END DO ! IEL

99999 END

