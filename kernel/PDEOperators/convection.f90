!##############################################################################
!# ****************************************************************************
!# <name> convection </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains discretisation and application routines for basic 
!# convection: Upwind and streamline diffusion as used in CCxD / PPxD.
!#
!# The following routines can be found in this module:
!#
!# 1.) conv_upwind2d
!#     -> Apply upwind to a vector, a matrix or both.
!#
!# 2.) conv_streamlinediffusion2d
!#     -> Apply streamline diffusion to a vector, a matrix or both.
!#
!# 3.) conv_JumpStabilisation2d
!#     -> Apply jump stabilisation to a vector, a matrix or both.
!# 
!# Some auxiliary routines:
!#
!# 1.) conv_ueoJumpStabil_double_uni
!#     -> Adds the Unified Edge Oriented stabilisation to a scalar matrix
!#        which is discretised based on a uniform discretisation.
!#
!# </purpose>
!##############################################################################

MODULE convection

  USE fsystem
  USE linearsystemscalar
  USE linearsystemblock
  USE cubature
  USE vectorio
  USE matrixio
  USE domainintegration
  USE bilinearformevaluation
  
  IMPLICIT NONE

!<constants>

!<constantblock description="Constants for cdef parameter of convection routine">

  ! Modify the matrix
  INTEGER, PARAMETER :: CONV_MODMATRIX = 2**0
  
  ! Set up defect vector
  INTEGER, PARAMETER :: CONV_MODDEFECT = 2**1
  
  ! Set up both, matrix and defect vector
  INTEGER, PARAMETER :: CONV_MODBOTH   = CONV_MODMATRIX+CONV_MODDEFECT

!</constantblock>

!<constantblock description="Constants that define the upwind type">

  ! Standard Samarskji or simple upwind.
  INTEGER, PARAMETER :: CONV_UPW_SAMARSKJI = 0
  
!</constantblock>

!<constantblock description="Constants that define the jump stabilisation type">

  ! No jump stabilisation
  INTEGER, PARAMETER :: CONV_JUMP_NONE        = 0
  
  ! Unified edge jump stabilisation
  INTEGER, PARAMETER :: CONV_JUMP_UNIFIEDEDGE = 1
  
!</constantblock>

!</constants>

!<types>

!<typeblock>
  
  ! Configuration block for standard UPWIND scheme.
  TYPE t_convUpwind
  
    ! Type of upwind. One of the CONV_UPW_xxxx constants.
    ! CONV_UPW_SAMARSKJI: Standard Samarskji or simple-type upwind
    INTEGER :: cupwType = CONV_UPW_SAMARSKJI
    
    ! Stabilisation parameter.
    ! If cupwType=CONV_UPW_SAMARSKJI:
    !  -1.0 = simple Upwind,
    !  >= 0:  Samarskji upwind with parameter $\theta=$ dupsam.
    !         Standard value = 0.1.
    REAL(DP) :: dupsam = 0.1_DP
    
    ! Whether the viscosity is constant.
    LOGICAL :: bconstViscosity = .TRUE.
    
    ! Viscosity parameter $\nu = 1/Re$ if viscosity is constant.
    ! We set this to infinity, what quickly leads to a program crash if the
    ! application does not initialise that properly!
    REAL(DP) :: dnu = SYS_INFINITY
    
    ! Weighting factor of the convective operator: $\theta * u*grad(u)$. 
    ! For time-dependent problems, this can be set to the step size
    ! in the $\Theta$-scheme.
    REAL(DP) :: dtheta = 1.0_DP
    
    ! Whether to use the ALE method for computing the convective operator.
    LOGICAL :: bALE = .FALSE.
    
  END TYPE
  
!</typeblock>

!<typeblock>
  
  ! Configuration block for Streamline Diffusion discretisation scheme
  TYPE t_convStreamlineDiffusion
  
    ! Stabilisation parameter.
    ! Standard value = 1.0_DP
    ! Note: A value of 0.0_DP sets up the convection part without
    ! any stabilisation (central-difference like discretisation).
    REAL(DP) :: dupsam = 1.0_DP
    
    ! Whether the viscosity is constant.
    LOGICAL :: bconstViscosity = .TRUE.
    
    ! Viscosity parameter $\nu = 1/Re$ if viscosity is constant
    REAL(DP) :: dnu = 1.0_DP
    
    ! Weighting factor for the mass matrix.
    ! =0.0: don't incorporate mass matrix into the operator.
    REAL(DP) :: dalpha = 0.0_DP
    
    ! Weighting factor for the Stokes matrix. (Stokes matrix = dnu*Laplace)
    ! =0.0: don't incorporate Stokes matrix into the operator.
    REAL(DP) :: dbeta = 0.0_DP
    
    ! Weighting factor for the convective part.
    ! =0.0: don't incorporate convective part (=Stokes problem)
    ! =1.0: incorporate full convective part (=Navier Stokes problem)
    REAL(DP) :: ddelta = 1.0_DP
    
    ! Weighting factor of the complete operator.
    ! For time-dependent problems, this can be set to the step size
    ! in the $\Theta$-scheme. For stationary problems, 1.0_DP must
    ! be used to assembly the not-weighted operator matrix.
    REAL(DP) :: dtheta = 1.0_DP
    
    ! Whether to use the ALE method for computing the convective operator.
    LOGICAL :: bALE = .FALSE.
    
  END TYPE
  
!</typeblock>

!<typeblock>

  ! Parameter block for Jump stabilisation
  TYPE t_jumpStabilisation
  
    ! Whether the viscosity is constant.
    LOGICAL :: bconstViscosity = .TRUE.
    
    ! Viscosity parameter $\nu = 1/Re$ if viscosity is constant
    REAL(DP) :: dnu            = 1.0_DP
    
    ! Line integral cubature formula for discretising the Jump.
    ! One of the CUB_xxxx_1D-constants of the cubature.f90 module.
    ! Standard is Gauss 2-point formula.
    INTEGER               :: ccubType = CUB_G2_1D
  
    ! Type of Jump stabilisation.
    ! One of the CONV_JUMP_xxxx-constants. Standard is unified edge
    ! stabilisation.
    INTEGER               :: cjump = CONV_JUMP_UNIFIEDEDGE
  
    ! 1st Relaxation parameter for the Jump stabilisation.
    ! =0.0: No Jump stabilisation
    ! >0.0: Add Jump stabilisation with relaxation parameter dgamma=djump.
    REAL(DP)              :: dgamma = 0.01_DP

    ! 2nd Relaxation parameter for the Jump stabilisation.
    ! =0.0: No Jump stabilisation
    ! >0.0: Add Jump stabilisation with relaxation parameter dgammastar
    REAL(DP)              :: dgammastar = 0.01_DP
  
    ! Weighting factor of the complete operator.
    ! For time-dependent problems, this can be set to the step size
    ! in the $\Theta$-scheme. For stationary problems, 1.0_DP must
    ! be used to assembly the not-weighted operator matrix.
    REAL(DP)              :: dtheta = 1.0_DP
    
  END TYPE
!</typeblock>

!</types>

CONTAINS

  ! ***************************************************************************

!<subroutine>

 SUBROUTINE conv_upwind2d (rvecPrimary, rvecSecondary, dprimWeight, dsecWeight,&
                           rconfig, cdef, &
                           rmatrix, rsolution, rdefect, DmeshVelocity, &
                           IvelocityComp)

!<description>
  ! Standard 1st order upwinding method to set up the convection operator
  !  
  !            $$ u_1 * grad(u_2) $$
  !
  ! in a matrix or to build a defect vector.
  ! 2D-version (X- and Y-velocity).
  !
  ! rvecPrimary, rvecSecondary are two velocity field vectors for the X-
  ! and Y-veclocity; IvelocityComp defines which components of these
  ! vectors contains the X- and which contains the Y-velocity.
  ! The final velocity vector field is then computed as a weighted average
  ! of these two:
  !
  !  $$ u_1  =  dprimWeight * rvecPrimary  +  dsecWeight * rvecSecondary $$
  !
  ! $u_2 = rsolution(.)$ defines the second velocity field inside of
  ! the grad-term.
  !
  ! The switch cdef decides on whether the routine sets up the nonlinear
  ! defect, the nonlinear matrix or both.
  !
  ! rmeshVelocity is an optional mesh velocity field that must be present
  ! if the ALE method should be used. 
  !
  ! The configuration how the routine should react is to be configured
  ! in the configuration block rconfig.
!</description>

!<input>

  ! Primary velocity field for the computation of $u_1$
  TYPE(t_vectorBlock), INTENT(IN), TARGET :: rvecPrimary
  
  ! Secondary velocity field for the computation of $u_1$
  TYPE(t_vectorBlock), INTENT(IN), TARGET :: rvecSecondary
  
  ! Weighting factor for rvecPrimary.
  REAL(DP), INTENT(IN) :: dprimWeight
  
  ! Weighting factor for rvecSecondary.
  REAL(DP), INTENT(IN) :: dsecWeight
  
  ! Configuration block for the upwind scheme
  TYPE(t_convUpwind), INTENT(IN) :: rconfig
  
  ! Computation/defect correction method. One of the CONV_MODxxxx constants:
  ! CONV_MODMATRIX: Set up the nonlinear matrix. rmatrix must be present, the 
  !                 nonlinear part is added to the matrix.
  ! CONV_MODDEFECT: Set up the nonlinear defect. rdefect and rsolution must be 
  !                 present.
  ! CONV_MODBOTH  : Set up the nonlinear matrix as well as the nonlinear defect.
  !                 rmatrix, rdefect and rsolution must all be present.
  INTEGER, INTENT(IN) :: cdef

  ! OPTIONAL: Solution vector u_2.
  ! Must be present if cdef=CONV_MODDEFECT or =CONV_MODBOTH.
  TYPE(t_vectorBlock), INTENT(IN), TARGET, OPTIONAL :: rsolution
  
  ! OPTIONAL: Mesh velocity field.
  ! DmeshVelocity(1,ivt) gives the X-velocity of the mesh, i.e. the X-velocity
  !   of the corner vertex ivt.
  ! DmeshVelocity(2,ivt) gives the Y-velocity of the mesh, i.e. the Y-velocity
  !   of the corner vertex ivt.
  ! The parameter must be present if ALE is activated in the
  ! configuration parameter block.
  REAL(DP), DIMENSION(:,:), INTENT(IN), OPTIONAL :: DmeshVelocity

  ! OPTIONAL: 
  ! Index block that specifies which component in rvecPrimary / rvecSecondary /
  ! rsolution / rdefect is the X-velocity and which one is the Y-velocity.
  !  IvelocityComp(1) gives the number of the X-velocity (usually = 1),
  !  IvelocityComp(2) gives the number of the Y-velocity (usually = 2).
  ! If not present, IvelocityComp=(/1,2/) is assumed, thus the X-velocity
  ! must be in rsolution\%RvectorBlock(1) and the Y-velocity in
  ! rsolution\%RvectorBlock(2).
  INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL :: IvelocityComp
!</input>

!<inputoutput>
  ! System matrix.
  ! The content of the matrix must be present if cdef=CONV_MODMATRIX or 
  ! =CONV_MODBOTH, otherwise only the structure is used.
  ! The nonlinear operator is added to the matrix.
  TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
  
  ! OPTIONAL: Defect vector.
  ! Must have the same structure as rsolution/rvecPrimary/rvecSecondary.
  ! Must be present if cdef=CONV_MODDEFECT or =CONV_MODBOTH.
  ! The nonlinear part is subtracted from this vector: 
  ! $r = r - \theta * u_1*grad(u_2)$
  TYPE(t_vectorBlock), INTENT(INOUT), OPTIONAL, TARGET :: rdefect
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: i
    INTEGER, DIMENSION(2) :: Icomp
    TYPE(t_vectorScalar), POINTER :: p_rvelX1,p_rvelX2,p_rvelY1,p_rvelY2
    TYPE(t_vectorScalar), POINTER :: p_rsolX,p_rsolY,p_rdefectX,p_rdefectY
    REAL(DP), DIMENSION(:), POINTER :: p_DvelX1,p_DvelX2,p_DvelY1,p_DvelY2
    REAL(DP), DIMENSION(:), POINTER :: p_DsolX,p_DsolY,p_DdefectX,p_DdefectY
    
    ! At first check the input parameters that everything is present what
    ! we need:
    IF ((cdef .EQ. CONV_MODDEFECT) .OR. (cdef .EQ. CONV_MODBOTH)) THEN
      IF ((.NOT. PRESENT(rsolution)) .OR. (.NOT. PRESENT(rdefect))) THEN
        PRINT *,'UPWIND: Solution/defect vector not present!'
        STOP
      END IF
    END IF
    
    IF (rconfig%bALE) THEN
      IF (.NOT. PRESENT(DmeshVelocity)) THEN
        PRINT *,'UPWIND: Mesh velocity vector not present!'
        STOP
      END IF
    END IF
    
    ! Get the actual subvectors from the velocity vectors that define
    ! the X- and Y-velocity.
    IF (PRESENT(IvelocityComp)) THEN
      Icomp = IvelocityComp
    ELSE
      Icomp = (/1,2/)
    END IF
    
    p_rvelX1 => rvecPrimary%RvectorBlock(Icomp(1))
    p_rvelY1 => rvecPrimary%RvectorBlock(Icomp(2))
    p_rvelX2 => rvecSecondary%RvectorBlock(Icomp(1))
    p_rvelY2 => rvecSecondary%RvectorBlock(Icomp(2))
    
    IF (PRESENT(rsolution)) THEN
      p_rsolX => rsolution%RvectorBlock(Icomp(1))
      p_rsolY => rsolution%RvectorBlock(Icomp(2))
    ELSE
      NULLIFY(p_rsolX)
      NULLIFY(p_rsolY)
    END IF
    
    IF (PRESENT(rsolution)) THEN
      p_rdefectX => rdefect%RvectorBlock(Icomp(1))
      p_rdefectY => rdefect%RvectorBlock(Icomp(2))
    ELSE
      NULLIFY(p_rdefectX)
      NULLIFY(p_rdefectY)
    END IF
      
    ! At the moment, we only support a rather limited set of configurations:
    ! Matrix and vectors must all be double precision, matrix must be format 
    ! 7 or 9, discretisation must be Q1~, constant viscosity.
    IF ((rmatrix%cmatrixFormat .NE. LSYSSC_MATRIX9) .AND. &
        (rmatrix%cmatrixFormat .NE. LSYSSC_MATRIX7)) THEN
      PRINT *,'UPWIND: Unsupported matrix format'
      STOP
    END IF

    i = rmatrix%p_rspatialDiscretisation%RelementDistribution(1)%itrialElement
    IF ((rmatrix%p_rspatialDiscretisation%ccomplexity .NE. SPDISC_UNIFORM) .OR. &
        ((i .NE. EL_E030) .AND. (i .NE. EL_E031) .AND. &
         (i .NE. EL_EM30) .AND. (i .NE. EL_EM31))) THEN
      PRINT *,'UPWIND: Unsupported discretisation.'
      STOP
    END IF

    IF ((rvecPrimary%cdataType .NE. ST_DOUBLE) .OR. &
        (rvecSecondary%cdataType .NE. ST_DOUBLE)) THEN
      PRINT *,'UPWIND: Unsupported vector data type in velocity.'
      STOP
    END IF
    
    IF (PRESENT(rdefect)) THEN
      IF ((rsolution%cdataType .NE. ST_DOUBLE) .OR. &
          (rdefect%cdataType .NE. ST_DOUBLE)) THEN
        PRINT *,'UPWIND: Unsupported vector data type in solution/defect'
        STOP
      END IF
    END IF
    
    IF (.NOT. rconfig%bconstViscosity) THEN
      PRINT *,'UPWIND: Only constant viscosity supported at the moment!'
      STOP
    END IF
    
    IF (rconfig%dnu .EQ. SYS_INFINITY) THEN
      PRINT *,'UPWIND: Viscosity parameter nu not initialised!'
      STOP
    END IF
    
    ! Call the actual calculation routine.
    ! Hide the p_rsol...-parameters to prevent passing the NULL()-pointer
    ! if rsolution is not present -- some compilers don't like that ^^

    CALL lsyssc_getbase_double (p_rvelX1,p_DvelX1)
    CALL lsyssc_getbase_double (p_rvelY1,p_DvelY1)
    CALL lsyssc_getbase_double (p_rvelX2,p_DvelX2)
    CALL lsyssc_getbase_double (p_rvelY2,p_DvelY2)
    
    !!! DEBUG:
    !WHERE (abs(p_DvelX1) .LT. 1E-12_DP) p_DvelX1 = 0.0_DP
    !WHERE (abs(p_DvelY1) .LT. 1E-12_DP) p_DvelY1 = 0.0_DP
    !CALL vecio_writeArray_Dble (p_DvelX1, 'vecx1', &
    !                               0, 'vectorx1.txt', '(D10.3)')
    !CALL vecio_writeArray_Dble (p_DvelY1, 'vecx2', &
    !                               0, 'vectorx2.txt', '(D10.3)')
    
    IF (PRESENT(rdefect)) THEN
      CALL lsyssc_getbase_double (p_rsolX   ,p_DsolX   )
      CALL lsyssc_getbase_double (p_rsolY   ,p_DsolY   )
      CALL lsyssc_getbase_double (p_rdefectX,p_DdefectX)
      CALL lsyssc_getbase_double (p_rdefectY,p_DdefectY)
      
      CALL conv_upwind2dALE_Q1Tdouble ( &
                    p_DvelX1,p_DvelY1,p_DvelX2,p_DvelY2,dprimWeight,dsecWeight, &
                    rmatrix,rmatrix%p_rspatialDiscretisation%p_rtriangulation, &
                    cdef, rconfig%dupsam, rconfig%dnu, rconfig%dtheta, &
                    rconfig%bALE, &
                    p_DsolX,p_DsolY,p_DdefectX,p_DdefectY, DmeshVelocity)
                    
    ELSE
    
      CALL conv_upwind2dALE_Q1Tdouble ( &
                    p_DvelX1,p_DvelY1,p_DvelX2,p_DvelY2,dprimWeight,dsecWeight, &
                    rmatrix,rmatrix%p_rspatialDiscretisation%p_rtriangulation, &
                    cdef, rconfig%dupsam, rconfig%dnu, rconfig%dtheta, &
                    rconfig%bALE, DmeshVelocity=DmeshVelocity)

      !!! DEBUG:
      !CALL matio_writeMatrixHR (rmatrix, 'matrix',&
      !                          .TRUE., 0, 'matrixL.txt', '(D10.3)')
                    
    END IF

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE conv_upwind2dALE_Q1Tdouble ( &
                  u1Xvel,u1Yvel,u2Xvel,u2Yvel,dweight1,dweight2,&
                  rmatrix,rtriangulation,cdef, &
                  dupsam,dnu,dtheta, bALE, &
                  Du1,Du2,Ddef1,Ddef2, DmeshVelocity)
!<description>
  ! Standard 1st order upwinding method to set up the convection operator
  !  
  !            $$ u_1 * grad(u_2) $$
  !
  ! in a matrix or to build a defect vector.
  ! 2D-version (X- and Y-velocity), uniform $\tilde Q_1$ discretisation,
  ! double precision vectors/matrix.
  !
  ! u1Xvel,u1Yvel, u2Xvel,u2Yvel are two velocity field vectors, 
  ! (u1Xvel,u1Yvel) a primary and (u2Xvel,u2Yvel) a secondary velocity field.
  ! The final velocity vector field is then computed as a weighted average
  ! of these two:
  !
  !  $$ u_1  =  dweight1 * u1vel  +  dweight2 * u2vel $$
  !
  ! $u_2 = rsolution(.)$ defines the second velocity field inside of
  ! the grad-term.
  !
  ! The switch cdef decides on whether the routine sets up the nonlinear
  ! defect, the nonlinear matrix or both.
  !
  ! DmeshVelocity is an optional mesh velocity field that must be present
  ! if the ALE method should be used (bALE=true). In this case, the nonlinear
  ! term is modified to include the mesh velocity.\\
  !
  ! For a reference about the ALE method, see
  ! [Duarte, Formaz, Natesan; "Arbitrary Lagrangian-Euler Method 
  ! for Navier-Stokes equations with moving boundaries";
  ! Comput. Methods Appl. Mech. Engrg. 193 (2004), 4819-4836]
  !
  ! Remarks:\\
  !  
  ! 1.) In a typical call of the upwinding, the caller can use:
  !     dweight1 = 1, u1Xvel/u1Yvel = velocity field
  !     dweight2 = 0, u2Xvel/u2Yvel = undefined
  !   So the upwinding scheme only uses one velocity field.
  !   Such a call e.g. adds the integral
  !                $$ ( u_1 * grad(.) , v )_{\Omega} $$
  !   to the system matrix.\\
  !
  !  2.) In case that there are two velocity fields representing
  !   the solution (may happen in a nonstationary simulation where
  !   u1Xvel/u1Yvel represents the solution in the current and u2Xvel/u2Yvel
  !   that of the previous time step), dweight1/dweight2 defines how these both
  !   velocity vectors should be weighted to compute the actual
  !   velocity field for the assembling:
  !               $$ U_act = dweight1*u1vel + dweight2*u2vel $$
  !   This is e.g. used for the linear extrapolation technique to
  !   reconstruct a velocity from two previous time steps...\\
  !
  !  3.) In the nonlinear iteration, as a right hand side there arises
  !   a defect vector D, which linear part can easily being assembled.
  !   However, there is a nonlinearity to be included into that vector,
  !   too. By setting cdef=1,2, this routine incorporates the nonlinearity
  !   into that vector, using the formula
  !
  !            $$ D = D - dtheta * UUx * grad (Ux) $$
  !   
  !  4.) If bALE=true, a mesh velocity field is added to the nonlineareity
  !   according to the formula  "U * grad (U-DmeshVelocity)".
  !   For bALE=false, the simple nonlinearity "U * grad (U)" is used.
  
!</description>

!<input>

  ! Primary X-velocity of $u_1$
  REAL(DP), DIMENSION(:), INTENT(IN) :: u1Xvel
  
  ! Primary Y-velocity of $u_1$
  REAL(DP), DIMENSION(:), INTENT(IN) :: u1Yvel
  
  ! Secondary X-velocity of $u_1$
  REAL(DP), DIMENSION(:), INTENT(IN) :: u2Xvel
  
  ! Secondary Y-velocity of $u_1$
  REAL(DP), DIMENSION(:), INTENT(IN) :: u2Yvel
  
  ! Computation/defect correction method. One of the CONV_MODxxxx constants:
  ! CONV_MODMATRIX: Set up the nonlinear matrix. rmatrix must be present, the 
  !                 nonlinear part is added to the matrix.
  ! CONV_MODDEFECT: Set up the nonlinear defect. rdefect and rsolution must be 
  !                 present.
  ! CONV_MODBOTH  : Set up the nonlinear matrix as well as the nonlinear defect.
  !                 rmatrix, rdefect and rsolution must all be present.
  INTEGER, INTENT(IN) :: cdef

  ! Weighting factor for u1Xvel/u1Yvel.
  REAL(DP), INTENT(IN) :: dweight1
  
  ! Weighting factor for u2Xvel/u2Yvel.
  REAL(DP), INTENT(IN) :: dweight2
  
  ! dupsam  - control parameter.
  !          -1: simple upwind,
  !          =0: Samarskji upwind
  REAL(DP), INTENT(IN) :: dupsam
  
  ! Viscosity parameter $\nu = 1/Re$ if viscosity is constant
  REAL(DP), INTENT(IN) :: dnu 

  ! Weighting factor of the convective operator: $\theta * u*grad(u)$. 
  ! For time-dependent problems, this can be set to the step size
  ! in the $\Theta$-scheme.
  REAL(DP), INTENT(IN) :: dtheta 
      
  ! Whether or not to use the ALE method
  LOGICAL, INTENT(IN) :: bALE
      
  ! OPTIONAL: Mesh velocity field. Must be present if bALE=TRUE.
  ! DmeshVelocity(1,:) gives the X-velocity of all the corner points of the mesh,
  ! DmeshVelocity(2,:) gives the Y-velocity.
  REAL(DP), DIMENSION(:,:), INTENT(IN), OPTIONAL :: DmeshVelocity(:,:)
  
  ! Triangulation structure specifying the underlying mesh.
  TYPE(t_triangulation), INTENT(IN) :: rtriangulation

  ! OPTIONAL: X-velocity of $u_2$. Must be present if cdef=CONV_MODDEFECT
  ! or cdef=CONV_MODBOTH.
  REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: Du1
  
  ! Y-velocity of $u_2$. Must be present if cdef=CONV_MODDEFECT
  ! or cdef=CONV_MODBOTH.
  REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: Du2
  
!</input>

!<inputoutput>
  ! The system matrix. Must be format 7 or 9.
  TYPE(t_matrixScalar), INTENT(INOUT), TARGET :: rmatrix
  
  ! OPTIONAL: X-defect vector. Must be present if cdef=CONV_MODDEFECT
  ! or =CONV_MODBOTH.
  REAL(DP), DIMENSION(:), INTENT(INOUT), OPTIONAL :: Ddef1(:)
  
  ! OPTIONAL: Y-defect vector. Must be present if cdef=CONV_MODDEFECT
  ! or =CONV_MODBOTH.
  REAL(DP), DIMENSION(:), INTENT(INOUT), OPTIONAL :: Ddef2(:)
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER(PREC_EDGEIDX), DIMENSION(4) :: Iedge
    INTEGER(PREC_MATIDX) :: IlocalMatrix(4,4)
    REAL(DP), DIMENSION(4) :: Dflux(4),Duu1(4),Duu2(4),XV(4),YV(4)
    REAL(DP), DIMENSION(4,4) :: DlocalMatrix(4,4)

    REAL(DP) :: dcenterX, dcenterY, XN, YN, G1, G2, dupsre, ELMH
    REAL(DP) :: DL0, DL2, H00, H22, dflux0, dflux2
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER(PREC_POINTIDX) :: iv, ivt1,ivt2
    INTEGER(PREC_EDGEIDX) :: im1, im0, im2
    INTEGER I,II
    INTEGER ia1, ia2, J, JJ, ia
    REAL(DP), DIMENSION(4) :: DuALE1,DuALE2

    REAL(DP), DIMENSION(:), POINTER :: p_Da
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld
    
    INTEGER(PREC_POINTIDX) :: NVT
    INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_Kvert
    INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_Kmid
    REAL(DP), DIMENSION(:,:), POINTER :: p_Dcorvg
    
    ! There is no additional type/completeness check in the parameters here. 
    ! This routine should never be called from outside, otherwise the main 
    ! application may have problems if the parameters are not specified correctly!
    !
    ! Get a pointer to the matrix entries and the matrix structure
    IF (rmatrix%h_DA .NE. ST_NOHANDLE) THEN
      ! The entries may be undefined - allowed if cdef=CONV_MODDEFECT!
      CALL storage_getbase_double (rmatrix%h_DA,p_Da)
    END IF
    CALL storage_getbase_int (rmatrix%h_Kcol,p_Kcol)
    CALL storage_getbase_int (rmatrix%h_Kld,p_Kld)
    
    ! Get a pointer to triangulation arrays.
    CALL storage_getbase_int2D (rtriangulation%h_IverticesAtElement,p_Kvert)
    CALL storage_getbase_int2D (rtriangulation%h_IedgesAtElement,p_Kmid)
    CALL storage_getbase_double2D (rtriangulation%h_DcornerCoordinates,p_Dcorvg)
    NVT = rtriangulation%NVT

    !********************************************************************
    !    Weighted Samarski upwind
    ! 
    !    This implementation follows the documentation of
    !    [F. Schieweck, Parallele L�sung der station�ren inkompressiblen
    !     Navier-Stokes Gleichungen, Habilitation, Fakult�t f�r
    !     Mathematik, Otto-von-Guericke-Universit�t Magdeburg]
    ! 
    !********************************************************************
    !
    !********************************************************************
    !     What we want to discretize here is the term
    !
    !         n(z,u,v) = ( (z*grad (.))u , v )_Omega
    !
    !     Let's assume we have two elements next to each other:
    !
    !       X---------------X---------------X
    !       |            /  |               |
    !       |          /    |               |
    !       |        /      |               |
    !       |       X       Gl              I
    !       |        \      |               |
    !       |         Glk   |               |
    !       |            \  |               |
    !       X------Gk-------X---------------X
    !
    !     The edges Gl and Gk enclose a diagonal edge Gklof the above
    !     triangle. The operator can now be rewritten by decomposition
    !     onto the elements as
    !
    !       n(z,u,v) ~= sum_l sum_k int_Gkl (z*n_lk) (u-u(Bl)) v(Bl) dGamma
    !
    !     with Bl and Bk being the midpoints of the edges and n_lk being the
    !     outer normal vector of the edge Glk. 
    !
    !       X---------------X              X---------------X
    !       |            /  |              |            /  |
    !       |          /    |              |          /    |
    !       |        /      |              |        /      |
    !       |       X       X Bl           |       X       u_l
    !       |        \      |              |        \      |
    !       |          \    |              |        u_upw  |
    !       |            \  |              |            \  |
    !       X-------X-------X              X------u_k------X
    !               Bk
    !
    !     The integral at the end of this term is replaced by 1x-Gauss 
    !     rule, thus u can be replaced by an approximation u_upw on the
    !     edge Glk - which is calculated with the help of the velocities 
    !     in the neighborhood u_l and u_k.
    !
    !     The main task in Upwinding is thus to calc u_upw -
    !     a reconstructed velocity, based on u_l and u_k !
    !     (Remember: With nonconforming, rotated bilinear elements,
    !      we have the velocity in the midpoints/edges - but not
    !      in the midpoints of the diagonals of those triangles.
    !      But there we need it because of an integration along this
    !      edge with simple Gauss rule.)
    !
    !     What's the approach here? As said, u_upw is reconstructed
    !     from u_1 and u_2 by a simple mean formula:
    !
    !          u_upw = Lambda u_l  + (1-Lambda) u_k
    !
    !     What is Lambda? 0 < Lambda < 1 is chosen depending on the flow
    !     crossing the diagonal Glk. More precisely, it's chosen
    !     depending on the flow:
    !       Flow direction       lambda        u_upw
    !          Bk -> Bl            ~0          ~u_k
    !          Bl -> Bk            ~1          ~u_l
    !          equal              ~0.5    ~mean between u_k and u_l
    !
    !     The "flow" is described by z. The "flow through Glk" or 
    !     "flux" is described by the line integral
    !
    !             t = 1/nu int_Glk (z*nlk) dGamma
    !
    !     (again with nlk being the normal vector of the edge Glk).
    !
    !     The parameter lambda is now chosen as:
    !
    !        lambda =    1 - 1/(2+2theta*t))  ,  t >= 0
    !               =    1/(2+2theta*t))      ,  t < 0
    !
    !     with theta being the dupsam-parameter from the DAT-file.
    !     So dupsam controls the weighting of the two neighboring
    !     velocities when reconstructing the velocity in the
    !     midpoint of the diagonal.
    !
    !     For theta=dupsam = 0, we have central difference.
    !     For theta=dupsam -> infinity, we obtain the simple upwind
    !     (lambda=0 for t<0, lambda=1 for t>=0).
    !
    !********************************************************************
    
    ! If the user wants Samarskji-Upwind, calculate an auxiliary variable
    ! dupsre: Weight the dupsam-parameter by 1/nu.
    ! This is needed later...

    IF (dupsam .GE. 0.0_DP) dupsre = dupsam / dnu
    
    ! Set dUale to 0, which is the standard contribution of ALE to the
    ! discretisation. If ALE is activated, this is changed on each cell later
    ! to the actual ALE contribution.
    DuALE1 = 0.0_DP
    DuALE2 = 0.0_DP
    
    ! We have a uniform grid, so we can simply loop over all elements
    ! and assemble the data element by element....
    !
    ! Loop over all elements in the current grid. The data is
    ! elementwise collected and added to the matrix.

    DO iel=1,rtriangulation%NEL

      ! dcenterX,dcenterY will be coordinates of the center of the element

      dcenterX=0.0_DP
      dcenterY=0.0_DP
      
      ! ALE handling on the current element.
      ! Loop over all 4 nodes to calculate the contribution of ALE:
      IF (bALE) THEN
        DO II=1,4
          ! In the ALE-equation, the nonlinear part is slightly modified:
          !
          !           u * grad(u)    --->    (u-v) * grad(u)
          !
          ! with v=dUale being the "mesh velocity field", i.e. an approximation
          ! to the velocity of the vertices in the mesh. Our mesh
          ! velocity field DmeshVelocity is given in the nodes. To compute
          ! it in the midpoints of the edges, we take the mean and
          ! save the result in DuALE1/2 for the X- and Y-coordinate.
          !
          ! At first, compute the start/endpoint of the edge II:

          ivt1 = p_Kvert(II,iel)
          ivt2 = p_Kvert(MOD(II,4)+1,iel)

          ! And then calculate the mean velocity field:

          DuALE1(II) = 0.5_DP * ( DmeshVelocity (1,ivt1) + DmeshVelocity(1,ivt2) )
          DuALE2(II) = 0.5_DP * ( DmeshVelocity (2,ivt1) + DmeshVelocity(2,ivt2) )
        END DO      
      END IF
      
      ! Loop over all 4 U-nodes.
      ! Calculate Iedge,XV,YV,dcenterX,dcenterY,Duu1,Duu2

      DO II=1,4
        
        ! Get the number of the II'th edge - which is at the same time the
        ! DOF in the velocity vector(s).
        I = p_Kmid(II,iel) - NVT
          
        ! Store the number of the edge/DOF in Iedge:
        Iedge(II)=I
          
        ! Store the coordinates of the corner vertices of that
        ! element in XV/YV:
          
        iv=p_Kvert(II,iel)
          
        XV(II)=p_Dcorvg(1,iv)
        YV(II)=p_Dcorvg(2,iv)

        ! Sum up the coordinates if the element - will later result
        ! in the element midpoint:
          
        dcenterX=dcenterX+XV(II)
        dcenterY=dcenterY+YV(II)
          
        ! Now we want to compute the velocity on the edge II
        ! (following the node II). 
        !
        ! Compute the actual velocity in the edge II (following the
        ! node II) by a weighted mean of the both velocity vectors
        ! U1Lx and U2Lx. This allowes e.g. to reconstruct a velocity
        ! vector by linear extrapolation of two previous time steps.
        !
        ! Subtract the mesh-velocity on the current edge as described
        ! above; if ALE is not used, DuALE is = 0, so nothing happens.
          
        Duu1(II) = dweight1*u1Xvel(I)+dweight2*u2Xvel(I) - DuALE1(II)
        Duu2(II) = dweight1*u1Yvel(I)+dweight2*u2Yvel(I) - DuALE2(II)
        
      END DO

      ! Divide dcenterX/dcenterY by 4 - so dcenterX/dcenterY receives the coordinate
      ! of the element midpoint

      dcenterX = 0.25_DP*dcenterX
      dcenterY = 0.25_DP*dcenterY
        
      !       After this procedure we have the following variable setting:
      !
      !   (XV(4),YV(4))               (XV(3),YV(3))
      !               X----Iedge(3)----X
      !               |                | 
      !               |                | 
      !               |                | 
      !         Iedge(4)       X      Iedge(2)
      !               |    (dcenter)   |
      !               |                | 
      !               |                |  
      !               X----Iedge(1)-- --X
      !   (XV(1),YV(1))               (XV(2),YV(2))
      !
      ! Duu1/Duu2 contains the velocity along the edges following
      ! the four corner points.
      !  
      ! Initialize DlocalMatrix(.,.) to 0. DlocalMatrix will assemble the "local"
      ! matrix, i.e. the values that are later incorporated into
      ! the system matrix at the positions stored in IlocalMatrix.
      !
      ! Loop over all 4 U-nodes.
      ! Calculate Dflux(.), IlocalMatrix(.,.), DlocalMatrix(.,.)
      
      DlocalMatrix = 0.0_DP

      DO II=1,4

        ! II1 corresponds to the current node in question.
        ! im1=II-1 MOD 4 receives the predecessor of II1 in an anticlockwise
        ! sense:
        im1=II-1
        IF (im1.LT.1) im1=4
      
        ! Calculation of the flux Dflux(II)
        !
        !                    /    |                        |
        !                  /      |                        |
        !              center     |            center      u_l
        !                  \      |                 \      |
        !                   /\    |                  GP    |
        !                  /   \  |                     \  |
        !         IM------/-------II       IM-----u_k------II
        !                / n
        !               v
        !
        ! From the mitpoint dcenterX/dcenterY and the current corner II,
        ! calculate the outer normal vector n of the edge Glk 
        ! that connects II with the midpoint:

        XN=-YV(II)+dcenterY
        YN= XV(II)-dcenterX
          
        ! Calculate the (scaled) flux
        !
        !   t = int_Glk (z*nlk) dGamma
        !     ~= nlk * u(GP)
        !      = 1/2 * nlk * (u_l+u_k)
        !
        ! approximating the integral with 1-point Gauss in the
        ! Gauss-point (=midpoint) of the edge. Save t in Dflux(II).
        ! So Dflux(II) saves the flux along the edge that is pointing
        ! from II to the element midpoint.
          
        G1=0.5_DP*(Duu1(im1)+Duu1(II))
        G2=0.5_DP*(Duu2(im1)+Duu2(II))
        Dflux(II)=XN*G1+YN*G2
      
        ! Determine the indices IlocalMatrix(II,JJ) to store the element matrix
        ! entry DlocalMatrix(II,JJ) on array A
        !
        ! Go into line I of the matrix, corresponding to the current
        ! edge Iedge(II). ia1 and ia2 saves the start- and end-index
        ! of this row of the matrix.

        I=Iedge(II)
        ia1=p_Kld(I)
        ia2=p_Kld(I+1)-1
          
        ! Loop over the edges of the element
          
        dofsearch: DO JJ=1,4

          ! In the current row, search for column J=Iedge(JJ).
          ! We know from FE-Theory, that our current node I interacts
          ! only with all the nodes on the edges of the current element,
          ! so only there can be an additional value to be included
          ! into the system matrix.
          !
          !   |---J---|
          !   |       |
          !   J       J
          !   |       |
          !   |---I---|
          !   |       |
          !
          ! So search in the current row for the matrix entry that
          ! corresponds to the common support of DOF I and J.

          J=Iedge(JJ)

          DO ia=ia1,ia2
            IF (p_Kcol(ia) .EQ. J) THEN
              ! Save the matrix index in IlocalMatrix(II,JJ) so we can find 
              ! the matrix entry later without searching for it.

              IlocalMatrix(II,JJ)=ia      
              
              ! Next DOF, no error
              CYCLE dofsearch

            END IF
          END DO

          ! Error case

          PRINT *,'ERROR in UPWIND: entry index ia not found'
          RETURN

        END DO dofsearch ! JJ
        
      END DO ! II
        
      ! What have we calculated up to here? Let's collect...
      !
      ! Dflux        - The flux along the edges of the triangles
      ! IlocalMatrix - The indices in the matrix A of the entries that
      !                are affected by integration on the current element iel.
      ! DlocalMatrix - Local element matrix - initialized with 0 up to now.
      !
      ! So the next step is to assemble the local element matrix,
      ! which is to be incorporated into the system matrix later.
      !
      ! Loop over the nodes to calculate DlocalMatrix:

      DO II=1,4
      
        ! Set im1=predecessor of II, im2=successor of II,
        ! in counterclockwise sense.

        im0=II
        im1=II-1
        IF (im1.LT.1) im1=4
        im2=II+1
        IF (im2.GT.4) im2=1

        ! We interpret II here as the local number of the edge
        ! (rather than the corner), following the corner II.
        ! The current triangle-edge GAMMA is the edge between
        ! corner II and the midpoint.
        !
        !    +------im2------im2
        !    |             / |
        !    |       dflux2  |
        !    |         /     |
        !    |       X       II=im0
        !    |         \     |
        !    |       dflux0  |
        !    |             \ |
        !  im1------im1------II=im0
        !
        ! Calculate the part corresponding to GAMMA(im0) and GAMMA(im2)
        !
        ! Get the flux of the edges im1->center, im2->center and
        ! save it to dflux0,dflux2:

        dflux0=Dflux(im0)
        dflux2=Dflux(im2)

        ! Now choose the Upwind scheme depending on dupsam:
        ! dupsam>0: Samarskji-Upwind
        ! dupsam<0: Simple upwind

        IF (dupsam .GE. 0.0_DP) THEN

          ! The user wants Samarskji-Upwind.
          ! Take dupsre, the dupsam-parameter, weighted by 1/nu.
          ! Remember: In the previous calculation of the line-integral
          ! to calculate t, we didn't incorporate 1/nu - this is 
          ! repaired here:
          !
          ! Analyze the two fluxes on the edges of the triangle.
          ! Calculate the lambda-value by Samarskji-upwind,
          ! depending on whether the flux goes "in" ou "out" of the
          ! triangle. DL0 saves the lambda of the triangle-edge
          ! corresponding to II, DL2 that of the triangle edge
          ! corresponding to  im2.

          IF (dflux0 .GE. 0.0_DP) THEN
            DL0=PHIP(dupsre*dflux0)
          ELSE
            DL0=PHIM(dupsre*dflux0)
          ENDIF

          IF (dflux2 .GE. 0.0_DP) THEN
            DL2=PHIP(dupsre*dflux2)
          ELSE
            DL2=PHIM(dupsre*dflux2)
          ENDIF

        ELSE

          ! Simple Upwinding scheme.
          ! The corresponding lambda (wighting factor of the
          ! "adjacent velocities) is either 0 or 1, depending on
          ! whether the flux goes "into" or "out of" the triangle
          ! on that edge.

          DL0=0.0_DP
          DL2=0.0_DP
          IF (dflux0 .GE. 0.0_DP) DL0=1.0_DP
          IF (dflux2 .GE. 0.0_DP) DL2=1.0_DP

        ENDIF

        ! Calculate the local element matrix with the integrals
        ! being evaluated by 1-point Gauss rule in the midpoint
        ! of the edges of the triangle. 
        !
        ! In fact, this calculates the local contribution of 
        !    UUx * grad ( . )

        H00=DL0*dflux0
        H22=(1.0_DP - DL2)*dflux2
        DlocalMatrix(im0,im0) =  H00-H22
        DlocalMatrix(im0,im2) =      H22
        DlocalMatrix(im0,im1) = -H00
        
      END DO ! II
      
      ! We are nearly done. Now we only have to incorporate
      ! the local element matrix DlocalMatrix into the global matrix.
      ! This is simple: Grab the index to be modified from IlocalMatrix
      ! and add the corresponding DlocalMatrix-value to that matrix entry.
      !
      ! Another matter of concern is the update of a defect vector,
      ! which arises during a nonlinear iteration. If cdef=1,2,
      ! we modify such a defect vector D in the form
      !
      !     D = D - dtheta * UUx * grad (Ux)
      !
      ! i.e. using a solution vector U (independent of the velocity
      ! fields), we subtract the nonlinearity from D. 
      !
      ! Weight the local matrix by dtheta. Remember, in the 
      ! nonlinear iteration we have to incorporate the term
      !
      !     THETA*K*u*grad(u)
      !
      ! which is realized by a multiplication of the u*grad(.)-
      ! matrix by THETA*K = dtheta here!
      
      DlocalMatrix = dtheta*DlocalMatrix
      
      ! Then incorporate the local matrix into the global matrix
      ! and/or the global defect.
      
      IF (IAND(cdef,CONV_MODMATRIX) .NE. 0) THEN

        DO JJ=1,4
          DO II=1,4
            ia       = IlocalMatrix(II,JJ)
            p_Da(ia) = p_Da(ia) + DlocalMatrix(II,JJ)
          END DO ! II
        END DO ! JJ

      END IF
      
      IF (IAND(cdef,CONV_MODDEFECT) .NE. 0) THEN
      
        DO II=1,4
          DO JJ=1,4
            ! Accessing the matrix via (II,1..4) is slower than accessing it
            ! via (1..4,JJ) would be, but we have no choice. Otherwise, we would 
            ! 'jump' through the memory in the defect creation process below 
            ! (Iedge(1..4)), which is even slower...
            
            ELMH = DlocalMatrix(II,JJ)

            ! Multiply the velocity Ux by the local matrix, resulting in
            ! the term UUx*grad(Ux) - and subtract that from the defect
            ! vector.

            Ddef1(Iedge(II)) = Ddef1(Iedge(II)) - ELMH*Du1(Iedge(JJ))
            Ddef2(Iedge(II)) = Ddef2(Iedge(II)) - ELMH*Du2(Iedge(JJ))
          END DO ! JJ
        END DO ! II
        
      END IF
      
    END DO ! IEL

  CONTAINS

    ! Auxiliary function 1
    ELEMENTAL REAL(DP) FUNCTION PHIP(x)
    REAL(DP), INTENT(IN) :: x
      PHIP = (0.5_DP+x)/(1.0_DP+x)
    END FUNCTION
    
    ! Auxiliary function 2
    ELEMENTAL REAL(DP) FUNCTION PHIM(x)
    REAL(DP), INTENT(IN) :: x
      PHIM = 0.5_DP/(1.0_DP-x)
    END FUNCTION

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE conv_streamlineDiffusion2d ( &
                           rvecPrimary, rvecSecondary, dprimWeight, dsecWeight,&
                           rconfig, cdef, &
                           rmatrix, rsolution, rdefect, DmeshVelocity, &
                           IvelocityComp)

!<description>
  ! Standard streamline diffusion method to set up the operator
  !  
  ! $$ dtheta  *  ( dalpha * MASS  +  dbeta * STOKES  +  ddelta * u_1 * grad(u_2) ) $$
  !
  ! in a matrix or to build a defect vector.
  ! 2D-version (X- and Y-velocity).
  !
  ! rvecPrimary, rvecSecondary are two velocity field vectors for the X-
  ! and Y-veclocity; IvelocityComp defines which components of these
  ! vectors contains the X- and which contains the Y-velocity.
  ! The final velocity vector field is then computed as a weighted average
  ! of these two:
  !
  !  $$ u_1  =  dprimWeight * rvecPrimary  +  dsecWeight * rvecSecondary $$
  !
  ! $u_2 = rsolution(.)$ defines the second velocity field inside of
  ! the grad-term.
  !
  ! The switch cdef decides on whether the routine sets up the nonlinear
  ! defect, the nonlinear matrix or both.
  !
  ! rmeshVelocity is an optional mesh velocity field that must be present
  ! if the ALE method should be used. 
  !
  ! The configuration how the routine should react is to be configured
  ! in the configuration block rconfig.
!</description>

!<input>

  ! Primary velocity field for the computation of $u_1$
  TYPE(t_vectorBlock), INTENT(IN), TARGET :: rvecPrimary
  
  ! Secondary velocity field for the computation of $u_1$
  TYPE(t_vectorBlock), INTENT(IN), TARGET :: rvecSecondary
  
  ! Weighting factor for rvecPrimary.
  REAL(DP), INTENT(IN) :: dprimWeight
  
  ! Weighting factor for rvecSecondary.
  REAL(DP), INTENT(IN) :: dsecWeight
  
  ! Configuration block for the streamline diffusion scheme
  TYPE(t_convStreamlineDiffusion), INTENT(IN) :: rconfig
  
  ! Computation/defect correction method. One of the CONV_MODxxxx constants:
  ! CONV_MODMATRIX: Set up the nonlinear matrix. rmatrix must be present, the 
  !                 nonlinear part is added to the matrix.
  ! CONV_MODDEFECT: Set up the nonlinear defect. rdefect and rsolution must be 
  !                 present.
  ! CONV_MODBOTH  : Set up the nonlinear matrix as well as the nonlinear defect.
  !                 rmatrix, rdefect and rsolution must all be present.
  INTEGER, INTENT(IN) :: cdef

  ! OPTIONAL: Solution vector u_2.
  ! Must be present if cdef=CONV_MODDEFECT or =CONV_MODBOTH.
  TYPE(t_vectorBlock), INTENT(IN), TARGET, OPTIONAL :: rsolution
  
  ! OPTIONAL: Mesh velocity field.
  ! DmeshVelocity(1,ivt) gives the X-velocity of the mesh, i.e. the X-velocity
  !   of the corner vertex ivt.
  ! DmeshVelocity(2,ivt) gives the Y-velocity of the mesh, i.e. the Y-velocity
  !   of the corner vertex ivt.
  ! The parameter must be present if ALE is activated in the
  ! configuration parameter block by bALE=true.
  REAL(DP), DIMENSION(:,:), INTENT(IN), OPTIONAL :: DmeshVelocity

  ! OPTIONAL: 
  ! Index block that specifies which component in rvecPrimary / rvecSecondary /
  ! rsolution / rdefect is the X-velocity and which one is the Y-velocity.
  !  IvelocityComp(1) gives the number of the X-velocity (usually = 1),
  !  IvelocityComp(2) gives the number of the Y-velocity (usually = 2).
  ! If not present, IvelocityComp=(/1,2/) is assumed, thus the X-velocity
  ! must be in rsolution\%RvectorBlock(1) and the Y-velocity in
  ! rsolution\%RvectorBlock(2).
  INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL :: IvelocityComp
!</input>

!<inputoutput>
  ! System matrix.
  ! The content of the matrix must be present if cdef=CONV_MODMATRIX or 
  ! =CONV_MODBOTH, otherwise only the structure is used.
  ! The nonlinear operator is added to the matrix.
  TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
  
  ! OPTIONAL: Defect vector.
  ! Must have the same structure as rsolution/rvecPrimary/rvecSecondary.
  ! Must be present if cdef=CONV_MODDEFECT or =CONV_MODBOTH.
  ! The nonlinear part is subtracted from this vector: 
  ! $r = r - \theta * u_1*grad(u_2)$
  TYPE(t_vectorBlock), INTENT(INOUT), OPTIONAL, TARGET :: rdefect
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: i
    INTEGER, DIMENSION(2) :: Icomp
    TYPE(t_vectorScalar), POINTER :: p_rvelX1,p_rvelX2,p_rvelY1,p_rvelY2
    TYPE(t_vectorScalar), POINTER :: p_rsolX,p_rsolY,p_rdefectX,p_rdefectY
    REAL(DP), DIMENSION(:), POINTER :: p_DvelX1,p_DvelX2,p_DvelY1,p_DvelY2
    REAL(DP), DIMENSION(:), POINTER :: p_DsolX,p_DsolY,p_DdefectX,p_DdefectY
    
    ! At first check the input parameters that everything is present what
    ! we need:
    IF ((cdef .EQ. CONV_MODDEFECT) .OR. (cdef .EQ. CONV_MODBOTH)) THEN
      IF ((.NOT. PRESENT(rsolution)) .OR. (.NOT. PRESENT(rdefect))) THEN
        PRINT *,'SD: Solution/defect vector not present!'
        STOP
      END IF
    END IF
    
    IF (rconfig%bALE) THEN
      IF (.NOT. PRESENT(DmeshVelocity)) THEN
        PRINT *,'SD: Mesh velocity vector not present!'
        STOP
      END IF
    END IF
    
    ! Get the actual subvectors from the velocity vectors that define
    ! the X- and Y-velocity.
    IF (PRESENT(IvelocityComp)) THEN
      Icomp = IvelocityComp
    ELSE
      Icomp = (/1,2/)
    END IF
    
    p_rvelX1 => rvecPrimary%RvectorBlock(Icomp(1))
    p_rvelY1 => rvecPrimary%RvectorBlock(Icomp(2))
    p_rvelX2 => rvecSecondary%RvectorBlock(Icomp(1))
    p_rvelY2 => rvecSecondary%RvectorBlock(Icomp(2))
    
    IF (PRESENT(rsolution)) THEN
      p_rsolX => rsolution%RvectorBlock(Icomp(1))
      p_rsolY => rsolution%RvectorBlock(Icomp(2))
    ELSE
      NULLIFY(p_rsolX)
      NULLIFY(p_rsolY)
    END IF
    
    IF (PRESENT(rsolution)) THEN
      p_rdefectX => rdefect%RvectorBlock(Icomp(1))
      p_rdefectY => rdefect%RvectorBlock(Icomp(2))
    ELSE
      NULLIFY(p_rdefectX)
      NULLIFY(p_rdefectY)
    END IF
      
    ! At the moment, we only support a rather limited set of configurations:
    ! Matrix and vectors must all be double precision, matrix must be format 
    ! 7 or 9, discretisation must be Q1~, constant viscosity.
    IF ((rmatrix%cmatrixFormat .NE. LSYSSC_MATRIX9) .AND. &
        (rmatrix%cmatrixFormat .NE. LSYSSC_MATRIX7)) THEN
      PRINT *,'SD: Unsupported matrix format'
      STOP
    END IF

    i = rmatrix%p_rspatialDiscretisation%RelementDistribution(1)%itrialElement
    IF (rmatrix%p_rspatialDiscretisation%ccomplexity .NE. SPDISC_UNIFORM) THEN
      PRINT *,'SD: Unsupported discretisation.'
      STOP
    END IF

    IF ((rvecPrimary%cdataType .NE. ST_DOUBLE) .OR. &
        (rvecSecondary%cdataType .NE. ST_DOUBLE)) THEN
      PRINT *,'SD: Unsupported vector data type in velocity.'
      STOP
    END IF
    
    IF (PRESENT(rdefect)) THEN
      IF ((rsolution%cdataType .NE. ST_DOUBLE) .OR. &
          (rdefect%cdataType .NE. ST_DOUBLE)) THEN
        PRINT *,'SD: Unsupported vector data type in solution/defect'
        STOP
      END IF
    END IF
    
    IF (.NOT. rconfig%bconstViscosity) THEN
      PRINT *,'SD: Only constant viscosity supported at the moment!'
      STOP
    END IF
    
    IF (rconfig%dnu .EQ. SYS_INFINITY) THEN
      PRINT *,'SD: Viscosity parameter nu not initialised!'
      STOP
    END IF
    
    ! Hide the p_rsol...-parameters to prevent passing the NULL()-pointer
    ! if rsolution is not present -- some compilers don't like that ^^

    CALL lsyssc_getbase_double (p_rvelX1,p_DvelX1)
    CALL lsyssc_getbase_double (p_rvelY1,p_DvelY1)
    CALL lsyssc_getbase_double (p_rvelX2,p_DvelX2)
    CALL lsyssc_getbase_double (p_rvelY2,p_DvelY2)
    
    !!! DEBUG:
    !WHERE (abs(p_DvelX1) .LT. 1E-12_DP) p_DvelX1 = 0.0_DP
    !WHERE (abs(p_DvelY1) .LT. 1E-12_DP) p_DvelY1 = 0.0_DP
    !CALL vecio_writeArray_Dble (p_DvelX1, 'vecx1', &
    !                               0, 'vectorx1.txt', '(D10.3)')
    !CALL vecio_writeArray_Dble (p_DvelY1, 'vecx2', &
    !                               0, 'vectorx2.txt', '(D10.3)')
    
    IF (PRESENT(rdefect)) THEN
      CALL lsyssc_getbase_double (p_rsolX   ,p_DsolX   )
      CALL lsyssc_getbase_double (p_rsolY   ,p_DsolY   )
      CALL lsyssc_getbase_double (p_rdefectX,p_DdefectX)
      CALL lsyssc_getbase_double (p_rdefectY,p_DdefectY)
      
      CALL conv_strdiff2dALE_double ( &
                    p_DvelX1,p_DvelY1,p_DvelX2,p_DvelY2,dprimWeight,dsecWeight, &
                    rmatrix,cdef, rconfig%dupsam, rconfig%dnu, &
                    rconfig%dalpha, rconfig%dbeta, rconfig%dtheta, rconfig%ddelta, &
                    rconfig%bALE, &
                    p_DsolX,p_DsolY,p_DdefectX,p_DdefectY, DmeshVelocity)
                    
    ELSE
    
      CALL conv_strdiff2dALE_double ( &
                    p_DvelX1,p_DvelY1,p_DvelX2,p_DvelY2,dprimWeight,dsecWeight, &
                    rmatrix, cdef, rconfig%dupsam, rconfig%dnu, &
                    rconfig%dalpha, rconfig%dbeta, rconfig%dtheta, rconfig%ddelta, &
                    rconfig%bALE, DmeshVelocity=DmeshVelocity)

      !!! DEBUG:
      !CALL matio_writeMatrixHR (rmatrix, 'matrix',&
      !                          .TRUE., 0, 'matrixL.txt', '(D10.3)')
                    
    END IF

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  SUBROUTINE conv_strdiff2dALE_double ( &
                  u1Xvel,u1Yvel,u2Xvel,u2Yvel,dweight1,dweight2,&
                  rmatrix,cdef, &
                  dupsam,dnu,dalpha,dbeta,dtheta, ddelta, bALE, &
                  Du1,Du2,Ddef1,Ddef2, DmeshVelocity)
!<description>
  ! Standard streamline diffusion method to set up the operator
  !  
  ! $$ dtheta  *  ( dalpha * MASS  +  dbeta * STOKES  +  ddelta * u_1 * grad(u_2) ) $$
  !
  ! in a matrix or to build a defect vector with that.
  ! 2D-version (X- and Y-velocity), uniform $\tilde Q_1$ discretisation,
  ! double precision vectors/matrix.
  !
  ! u1Xvel,u1Yvel, u2Xvel,u2Yvel are two velocity field vectors, 
  ! (u1Xvel,u1Yvel) a primary and (u2Xvel,u2Yvel) a secondary velocity field.
  ! The final velocity vector field is then computed as a weighted average
  ! of these two:
  !
  !  $$ u_1  =  dweight1 * u1vel  +  dweight2 * u2vel $$
  !
  ! $u_2 = rsolution(.)$ defines the second velocity field inside of
  ! the grad-term.
  !
  ! The switch cdef decides on whether the routine sets up the nonlinear
  ! defect, the nonlinear matrix or both.
  !
  ! DmeshVelocity is an optional mesh velocity field that must be present
  ! if the ALE method should be used (bALE=true). In this case, the nonlinear
  ! term is modified to include the mesh velocity.\\
  !
  ! For a reference about the ALE method, see
  ! [Duarte, Formaz, Natesan; "Arbitrary Lagrangian-Euler Method 
  ! for Navier-Stokes equations with moving boundaries";
  ! Comput. Methods Appl. Mech. Engrg. 193 (2004), 4819-4836]
  !
  ! Remarks:\\
  !  
  ! 1.) In a typical call of the upwinding, the caller can use:
  !     dweight1 = 1, u1Xvel/u1Yvel = velocity field
  !     dweight2 = 0, u2Xvel/u2Yvel = undefined
  !   So the upwinding scheme only uses one velocity field.
  !   Such a call e.g. adds the integral
  !                $$ ( u_1 * grad(.) , v )_{\Omega} $$
  !   to the system matrix.\\
  !
  !  2.) In case that there are two velocity fields representing
  !   the solution (may happen in a nonstationary simulation where
  !   u1Xvel/u1Yvel represents the solution in the current and u2Xvel/u2Yvel
  !   that of the previous time step), dweight1/dweight2 defines how these both
  !   velocity vectors should be weighted to compute the actual
  !   velocity field for the assembling:
  !               $$ U_act = dweight1*u1vel + dweight2*u2vel $$
  !   This is e.g. used for the linear extrapolation technique to
  !   reconstruct a velocity from two previous time steps...\\
  !
  !  3.) In the nonlinear iteration, as a right hand side there arises
  !   a defect vector D, which linear part can easily being assembled.
  !   However, there is a nonlinearity to be included into that vector,
  !   too. By setting cdef=1,2, this routine incorporates the nonlinearity
  !   into that vector, using the formula
  !
  !            $$ D = D - dtheta * UUx * grad (Ux) $$
  !   
  !  4.) If bALE=true, a mesh velocity field is added to the nonlineareity
  !   according to the formula  "U * grad (U-DmeshVelocity)".
  !   For bALE=false, the simple nonlinearity "U * grad (U)" is used.
  
!</description>

!<input>

  ! Primary X-velocity of $u_1$
  REAL(DP), DIMENSION(:), INTENT(IN) :: u1Xvel
  
  ! Primary Y-velocity of $u_1$
  REAL(DP), DIMENSION(:), INTENT(IN) :: u1Yvel
  
  ! Secondary X-velocity of $u_1$
  REAL(DP), DIMENSION(:), INTENT(IN) :: u2Xvel
  
  ! Secondary Y-velocity of $u_1$
  REAL(DP), DIMENSION(:), INTENT(IN) :: u2Yvel
  
  ! Computation/defect correction method. One of the CONV_MODxxxx constants:
  ! CONV_MODMATRIX: Set up the nonlinear matrix. rmatrix must be present, the 
  !                 nonlinear part is added to the matrix.
  ! CONV_MODDEFECT: Set up the nonlinear defect. rdefect and rsolution must be 
  !                 present.
  ! CONV_MODBOTH  : Set up the nonlinear matrix as well as the nonlinear defect.
  !                 rmatrix, rdefect and rsolution must all be present.
  INTEGER, INTENT(IN) :: cdef

  ! Weighting factor for u1Xvel/u1Yvel.
  REAL(DP), INTENT(IN) :: dweight1
  
  ! Weighting factor for u2Xvel/u2Yvel.
  REAL(DP), INTENT(IN) :: dweight2
  
  ! dupsam  - control parameter.
  !          -1: simple upwind,
  !          =0: Samarskji upwind
  REAL(DP), INTENT(IN) :: dupsam
  
  ! Viscosity parameter $\nu = 1/Re$ if viscosity is constant
  REAL(DP), INTENT(IN) :: dnu 
  
  ! Weighting factor for the mass matrix.
  REAL(DP), INTENT(IN) :: dalpha

  ! Weighting factor for the Stokes matrix. (Stokes matrix = 1/Re * Laplace)
  REAL(DP), INTENT(IN) :: dbeta

  ! Weighting factor of the convective operator: $\theta * u*grad(u)$. 
  ! For time-dependent problems, this can be set to the step size
  ! in the $\Theta$-scheme.
  REAL(DP), INTENT(IN) :: dtheta 
  
  ! Weighting factor for the nonlinear term
  REAL(DP), INTENT(IN) :: ddelta
      
  ! Whether or not to use the ALE method
  LOGICAL, INTENT(IN) :: bALE
      
  ! OPTIONAL: Mesh velocity field. Must be present if bALE=TRUE.
  ! DmeshVelocity(1,:) gives the X-velocity of all the corner points of the mesh,
  ! DmeshVelocity(2,:) gives the Y-velocity.
  REAL(DP), DIMENSION(:,:), INTENT(IN), OPTIONAL :: DmeshVelocity(:,:)
  
  ! OPTIONAL: X-velocity of $u_2$. Must be present if cdef=CONV_MODDEFECT
  ! or cdef=CONV_MODBOTH.
  REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: Du1
  
  ! Y-velocity of $u_2$. Must be present if cdef=CONV_MODDEFECT
  ! or cdef=CONV_MODBOTH.
  REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: Du2
  
!</input>

!<inputoutput>
  ! The system matrix. Must be format 7 or 9.
  TYPE(t_matrixScalar), INTENT(INOUT), TARGET :: rmatrix
  
  ! OPTIONAL: X-defect vector. Must be present if cdef=CONV_MODDEFECT
  ! or =CONV_MODBOTH.
  REAL(DP), DIMENSION(:), INTENT(INOUT), OPTIONAL :: Ddef1(:)
  
  ! OPTIONAL: Y-defect vector. Must be present if cdef=CONV_MODDEFECT
  ! or =CONV_MODBOTH.
  REAL(DP), DIMENSION(:), INTENT(INOUT), OPTIONAL :: Ddef2(:)
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: indof,indofALE,csysTrial,IEQ,I,J,K,IDOFE,JDOFE,icubp
  INTEGER(PREC_DOFIDX) :: JCOL0,IDFG,JDFG,JCOL
  INTEGER(PREC_ELEMENTIDX) :: IEL,IELset,IELmax
  LOGICAL, DIMENSION(EL_MAXNDER) :: Bder,BderALE
  LOGICAL :: bnonpar
  REAL(DP) :: dumax,dumaxr, du1loc, du2loc, dunorm,db,OM,AH,denth,dre,dny
  REAL(DP) :: HBASI1,HBASI2,HBASI3,HBASJ1,HBASJ2,HBASJ3,HSUMI,HSUMJ
  INTEGER :: NVE
  
  ! Matrix structure arrays
  INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol
  INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld
  REAL(DP), DIMENSION(:), POINTER :: p_Da
  
  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPtsRef

  ! An array receiving the coordinates of cubature points on
  ! the real element for all elements in a set.
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPtsReal

  ! The discretisation - for easier access
  TYPE(t_spatialDiscretisation), POINTER :: p_rdiscretisation
  
  ! Triangulation
  TYPE(t_triangulation), POINTER :: p_rtriangulation
  REAL(DP), DIMENSION(:,:), POINTER :: p_DcornerCoordinates
  INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement,p_IverticesAtElement

  ! Number of elements in a block. Normally =BILF_NELEMSIM,
  ! except if there are less elements in the discretisation.
  INTEGER :: nelementsPerBlock

  ! One and only element distribution
  TYPE(t_elementDistribution), POINTER :: p_relementDistribution

  ! Cubature point coordinates on the reference element
  REAL(DP), DIMENSION(CUB_MAXCUBP, NDIM3D) :: Dxi

  ! For every cubature point on the reference element,
  ! the corresponding cubature weight
  REAL(DP), DIMENSION(CUB_MAXCUBP) :: Domega
  
  ! number of cubature points on the reference element
  INTEGER :: ncubp

  ! Pointer to the point coordinates to pass to the element function.
  ! Point either to p_DcubPtsRef or to p_DcubPtsReal, depending on whether
  ! the trial/test element is parametric or not.
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPts

  ! A t_domainIntSubset structure for integrating over the domain.
  TYPE(t_domainIntSubset) :: rintSubset

  ! Array with coordinates of the corners that form the real element.
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_Dcoords
  
  ! Arrays for saving Jacobian determinants and matrices
  REAL(DP), DIMENSION(:,:), POINTER :: p_Ddetj
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_Djac
  
  ! An allocateable array accepting the DOF's of a set of elements.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE, TARGET :: Idofs, IdofsALE
  
  ! Allocateable arrays for the values of the basis functions - 
  ! for test and trial spaces.
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: Dbas,DbasALE

  ! Local matrices, used during the assembly.
  ! Values and positions of values in the global matrix.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:,:), ALLOCATABLE :: Kentry
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Dentry

  ! A pointer to an element-number list
  INTEGER(I32), DIMENSION(:), POINTER :: p_IelementList

  ! Pointer to the velocity field in the cubature points.
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Dvelocity
  
  ! An array with local DELTA's, each DELTA for one element
  REAL(DP), DIMENSION(:), ALLOCATABLE :: DlocalDelta
  
    ! Initialise the derivative flags
    Bder = .FALSE.
    Bder(DER_FUNC) = .TRUE.
    Bder(DER_DERIV_X) = .TRUE.
    Bder(DER_DERIV_Y) = .TRUE.

    ! For ALE we don't even need so much
    BderALE = .FALSE.
    BderALE(DER_FUNC) = .TRUE.
    
    ! Shortcut to the spatial discretisation
    p_rdiscretisation => rmatrix%p_rspatialDiscretisation
    
    ! Get the element distribution. Here, we can find information about
    ! the cubature formula etc...
    p_relementDistribution => p_rdiscretisation%RelementDistribution(1)
    
    ! Get some information about the triangulation
    p_rtriangulation => p_rdiscretisation%p_rtriangulation
    CALL storage_getbase_double2d (p_rtriangulation%h_DcornerCoordinates,&
                                   p_DcornerCoordinates)
    CALL storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
                                p_IverticesAtElement)
    CALL storage_getbase_int2d (p_rtriangulation%h_IedgesAtElement,&
                                p_IedgesAtElement)
    
    ! Get the number of local DOF's for trial/test functions.
    ! We assume trial and test functions to be the same.
    indof = elem_igetNDofLoc(p_relementDistribution%itrialElement)

    ! Get the number of local DOF's Q1 -- we need them for ALE.
    indofALE = elem_igetNDofLoc(p_relementDistribution%itrialElement)
    
    ! Number of local DOF's
    NVE = elem_igetNVE(p_relementDistribution%itrialElement)
    
    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it's at most
    ! BILF_NELEMSIM. This is only used for allocaing some arrays.
    nelementsPerBlock = MIN(BILF_NELEMSIM,p_rtriangulation%NEL)
    
    ! For cdef containing CONV_MODDEFECT, we build the defect vector                     
    !     D = RHS - A*U                                         
    ! In this case, the defect(rhs vectors must be present
    
    IF (IAND(cdef,CONV_MODDEFECT) .NE. 0) THEN
      IF (.NOT. (PRESENT(Ddef1) .AND. PRESENT(Ddef2) .AND. &
                 PRESENT(Du1) .AND. PRESENT(Du2))) THEN
        PRINT *,'conv_strdiff2dALE_double: Necessary arguments missing!'
        STOP
      END IF
    END IF
    
    IF (IAND(cdef,CONV_MODMATRIX) .NE. 0) THEN
      ! Get matrix arrays
      CALL storage_getbase_double (rmatrix%h_Da,p_Da)
    END IF
    CALL storage_getbase_int (rmatrix%h_Kcol,p_Kcol)
    CALL storage_getbase_int (rmatrix%h_Kld,p_Kld)
    
    ! Initialise the cubature formula,
    ! get cubature weights and point coordinates on the reference element
    CALL cub_getCubPoints(p_relementDistribution%ccubTypeBilForm, ncubp, Dxi, Domega)
    
    ! Get from the trial element space the type of coordinate system
    ! that is used there:
    csysTrial = elem_igetCoordSystem(p_relementDistribution%itrialElement)

    ! Allocate memory and get local references to it.
    CALL domint_initIntegration (rintSubset,nelementsPerBlock,ncubp,csysTrial,NDIM2D)
    p_DcubPtsRef =>  rintSubset%p_DcubPtsRef
    p_DcubPtsReal => rintSubset%p_DcubPtsReal
    p_Djac =>        rintSubset%p_Djac
    p_Ddetj =>       rintSubset%p_Ddetj
    p_Dcoords =>     rintSubset%p_DCoords
    
    ! Put the cubature point coordinates in the right format to the
    ! cubature-point array.
    ! Initialise all entries in p_DcubPtsRef with the same coordinates -
    ! as the cubature point coordinates are identical on all elements
    DO j=1,SIZE(p_DcubPtsRef,3)
      DO i=1,ncubp
        DO k=1,SIZE(p_DcubPtsRef,1)
          p_DcubPtsRef(k,i,j) = Dxi(i,k)
        END DO
      END DO
    END DO

    ! Allocate an array saving the coordinates of corner vertices of elements
    
    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  ALLOCATE(Dbas(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    ! would lead to nonused memory blocks in these arrays during the assembly, 
    ! which reduces the speed by 50%!
    ALLOCATE(Dbas(indof,elem_getMaxDerivative(p_relementDistribution%itrialElement), &
             ncubp,nelementsPerBlock))

    ! Allocate memory for the DOF's of all the elements.
    ALLOCATE(Idofs(indof,nelementsPerBlock))
    
    ! The same for the ALE-space
    ALLOCATE(DbasALE(indofALE,elem_getMaxDerivative(EL_Q1), &
             ncubp,nelementsPerBlock))

    ! Allocate memory for the DOF's of all the elements.
    ALLOCATE(IdofsALE(indofALE,nelementsPerBlock))
    
    ! Allocate memory for array with local DELTA's
    ALLOCATE(DlocalDelta(nelementsPerBlock))

    ! Check if one of the trial/test elements is nonparametric
    bnonpar = elem_isNonparametric(p_relementDistribution%itrialElement)
    
    ! Let p_DcubPtspoint either to p_DcubPtsReal or
    ! p_DcubPtsRef - depending on whether the space is parametric or not.
    IF (bnonpar) THEN
      p_DcubPts => p_DcubPtsReal
    ELSE
      p_DcubPts => p_DcubPtsRef
    END IF
    
    ! Allocate an array saving the local matrices for all elements
    ! in an element set.
    ! We could also allocate EL_MAXNBAS*EL_MAXNBAS*BILF_NELEMSIM integers
    ! for this local matrix, but this would normally not fit to the cache
    ! anymore! indofTrial*indofTest*BILF_NELEMSIM is normally much smaller!
    ALLOCATE(Kentry(indof,indof,nelementsPerBlock))
    ALLOCATE(Dentry(indof,indof))
    
    ! Allocate memory for the velocity in the cubature points.
    ALLOCATE(Dvelocity(NDIM2D,ncubp,nelementsPerBlock))
    
    ! What is the reciprocal of nu? We need it later.
    IF (dnu .NE. 0.0_DP) THEN
      dre = 1.0_DP/dnu
      
      ! dny gets the actual multiplier for the Laplace matrix.
      ! Remember: dbeta*Stokes = dbeta*dnu*Laplace = dny*Laplace.
      ! This may be =0.0 if the Stokes operator should not be included into
      ! the matrix.
      dny = dbeta*dnu
    ELSE
      PRINT *,'SD: NU=0 not allowed! Set dbeta=0 to prevent Stokes operator'// &
              ' from being build!'
      STOP
    END IF

    ! If ddelta=0, we have to neglect the nonlinearity. In both cases,
    ! set DlocalDelta=0 which disables the nonlinear term in the assembly.
    ! If dupsam=0, we neglect the stabilisation term (central difference like
    ! discretisation), so we set DlocalDelta=0 as well.
    IF ((ddelta .EQ. 0.0_DP) .OR. (dupsam .EQ. 0.0_DP)) THEN
      CALL lalg_clearVectorDble (DlocalDelta)
    END IF
    
    ! Calculate the maximum norm of the actual velocity field
    ! U = A1*U1 + A2*U2 into DUMAX. 
    ! Round up the norm to 1D-8 if it's too small...

    dumax=0.0_DP
    IF (dweight2 .EQ. 0.0_DP) THEN
      DO IEQ=1,SIZE(u1Xvel)
        du1loc = u1Xvel(IEQ)
        du2loc = u1Yvel(IEQ)
        dunorm = SQRT(du1loc**2+du2loc**2)
        dumax = MAX(DUMAX,DUNORM)
      END DO
    ELSE       
      DO ieq=1,SIZE(u1Xvel)
        du1loc = dweight1*u1Xvel(IEQ)+dweight2*u2Xvel(IEQ)
        du2loc = dweight1*u1Yvel(IEQ)+dweight2*u2Yvel(IEQ)
        dunorm = SQRT(du1loc**2+du2loc**2)
        dumax = MAX(dumax,dunorm)
      END DO
    ENDIF       

    IF (dumax.LT.1E-8_DP) dumax=1E-8_DP
    dumaxr = 1.0_DP/dumax

    ! p_IelementList must point to our set of elements in the discretisation
    ! with that combination of trial/test functions
    CALL storage_getbase_int (p_relementDistribution%h_IelementList, &
                              p_IelementList)
    
    ! Loop over the elements - blockwise.
    DO IELset = 1, p_rtriangulation%NEL, BILF_NELEMSIM
    
      ! We always handle BILF_NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most BILF_NELEMSIM
      ! elements simultaneously.
      
      IELmax = MIN(p_rtriangulation%NEL,IELset-1+BILF_NELEMSIM)
    
      ! The outstanding feature with finite elements is: A basis
      ! function for a DOF on one element has common support only
      ! with the DOF's on the same element! E.g. for Q1:
      !
      !        #. . .#. . .#. . .#
      !        .     .     .     .
      !        .  *  .  *  .  *  .
      !        #-----O-----O. . .#
      !        |     |     |     .
      !        |     | IEL |  *  .
      !        #-----X-----O. . .#
      !        |     |     |     .
      !        |     |     |  *  .
      !        #-----#-----#. . .#
      !
      ! --> On element IEL, the basis function at "X" only interacts
      !     with the basis functions in "O". Elements in the 
      !     neighbourhood ("*") have no support, therefore we only have
      !     to collect all "O" DOF's.
      !
      ! Calculate the global DOF's into IdofsTrial / IdofsTest.
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF's of our BILF_NELEMSIM elements simultaneously.
      CALL dof_locGlobMapping_mult(p_rdiscretisation, p_IelementList(IELset:IELmax), &
                                  .TRUE.,Idofs)
                                  
      ! In case ALE is used, do this also for the ALE stuff.
      IF (bALE) THEN
        CALL dof_locGlobMapping_mult(p_rdiscretisation, &
                                    p_IelementList(IELset:IELmax), &
                                    .TRUE.,IdofsALE)
      END IF
      
      ! Calculate local DELTA's for streamline diffusion method.
      ! (cf. p. 121 in Turek's CFD book).
      ! For every element, we need a local DELTA.
      ! Every local delta is weighted by the global "ddelta".
      ! If ddelta=0, we don't do anything as this disables the
      ! nonlinear term.
      ! If UPSAM=0.0, we have a central-difference like discretisation, which
      ! is one can see as the local stabilisation weight Delta is also = 0.0.
      ! In this case, we even switch of the calculation of the local Delta,
      ! as it is always =0.0, so we save a little bit time.
      IF ((ddelta .NE. 0.0_DP) .AND. (dupsam .NE. 0.0_DP))THEN
        DO IEL=1,IELmax-IELset+1
          CALL getLocalDelta (u1Xvel,u1Yvel,u2Xvel,u2Yvel,dweight1,dweight2, &
                      IEL,DUMAXR,DlocalDelta(IEL), &
                      p_IverticesAtElement,p_IedgesAtElement,&
                      p_DcornerCoordinates,Idofs(:,IEL),indof, &
                      dupsam,dre)
          DlocalDelta(IEL) = ddelta*DlocalDelta(IEL)
        END DO ! IEL
      END IF
                                   
      ! For the assembly of the global matrix, we use a "local"
      ! approach. At first we build a "local" system matrix according
      ! to the current element. This contains all additive
      ! contributions of element IEL, which are later added at the
      ! right positions to the elements in the global system matrix.
      !
      ! We have indofTrial trial DOF's per element and
      ! indofTest test DOF's per element. Therefore there are
      ! indofTrial*indofTest tupel of basis-/testfunctions (phi_i,psi_j) 
      ! "active" (i.e. have common support) on our current element, each 
      ! giving an additive contribution to the system matrix.
      !
      ! We build a quadratic indofTrial*indofTest local matrix:
      ! Kentry(1..indofTrial,1..indofTest) receives the position 
      !   in the global system matrix, where the corresponding value 
      !   has to be added to.
      ! (The corresponding contrbutions can be saved separately, 
      !  but we directly add them to the global matrix in this 
      !  approach.)
      !
      ! We build local matrices for all our elements 
      ! in the set simultaneously.
      ! Loop through elements in the set and for each element,
      ! loop through the local matrices to initialise them:
      DO IEL=1,IELmax-IELset+1
      
        ! For building the local matrices, we have first to
        ! loop through the test functions (the "O"'s), as these
        ! define the rows in the matrix.
        DO IDOFE=1,indof
        
          ! Row IDOFE of the local matrix corresponds 
          ! to row=global DOF KDFG(IDOFE) in the global matrix.
          ! This is one of the the "O"'s in the above picture.
          ! Get the starting position of the corresponding row
          ! to JCOL0:

          JCOL0=p_KLD(Idofs(IDOFE,IEL))
          
          ! Now we loop through the other DOF's on the current element
          ! (the "O"'s).
          ! All these have common support with our current basis function
          ! and will therefore give an additive value to the global
          ! matrix.
          
          DO JDOFE=1,indof
            
            ! Get the global DOF of the "X" which interacts with 
            ! our "O".
            
            JDFG=Idofs(JDOFE,IEL)
            
            ! Starting in JCOL0 (which points to the beginning of
            ! the line initially), loop through the elements in
            ! the row to find the position of column IDFG.
            ! Jump out of the DO loop if we find the column.
            
            DO JCOL=JCOL0,rmatrix%NA
              IF (p_KCOL(JCOL) .EQ. JDFG) EXIT
            END DO

            ! Because columns in the global matrix are sorted 
            ! ascendingly (except for the diagonal element),
            ! the next search can start after the column we just found.
            
            ! JCOL0=JCOL+1
            
            ! Save the position of the matrix entry into the local
            ! matrix.
            ! Note that a column in Kentry corresponds to a row in
            ! the real matrix. We aligned Kentry/DENTRY this way to get
            ! higher speed of the assembly routine, since this leads
            ! to better data locality.
            
            Kentry(JDOFE,IDOFE,IEL)=JCOL
            
          END DO ! IDOFE
          
        END DO ! JDOFE
        
      END DO ! IEL
      
      ! Ok, we found the positions of the local matrix entries
      ! that we have to change.
      ! To calculate the matrix contributions, we have to evaluate
      ! the elements to give us the values of the basis functions
      ! in all the DOF's in all the elements in our set.
      !
      ! We have the coordinates of the cubature points saved in the
      ! coordinate array from above. Unfortunately for nonparametric
      ! elements, we need the real coordinate.
      ! Furthermore, we anyway need the coordinates of the element
      ! corners and the Jacobian determinants corresponding to
      ! all the points.
      !
      ! At first, get the coordinates of the corners of all the
      ! elements in the current set. 
      
!      DO IEL=1,IELmax-IELset+1
!        p_Dcoords(:,:,IEL) = p_DcornerCoordinates(:, &
!                            p_IverticesAtElement(:,p_IelementList(IELset+IEL-1)))
!      END DO
      DO IEL=1,IELmax-IELset+1
        DO J = 1,NVE
          DO I = 1,NDIM2D
            p_Dcoords(I,J,IEL) = p_DcornerCoordinates(I, &
                               p_IverticesAtElement(J,p_IelementList(IELset+IEL-1)))
          END DO
        END DO
      END DO
      
      ! Depending on the type of transformation, we must now choose
      ! the mapping between the reference and the real element.
      ! In case we use a nonparametric element or a nonconstant coefficient function,
      ! we need the coordinates of the points on the real element, too.
      IF (bnonpar) THEN
      
        CALL trafo_calctrafo_sim (p_relementDistribution%ctrafoType,&
             IELmax-IELset+1,ncubp,p_Dcoords,&
             p_DcubPtsRef,p_Djac(:,:,1:IELmax-IELset+1),p_Ddetj(:,1:IELmax-IELset+1), &
             p_DcubPtsReal)
      
      ELSE
      
        CALL trafo_calctrafo_sim (p_relementDistribution%ctrafoType,&
             IELmax-IELset+1,ncubp,p_Dcoords,&
             p_DcubPtsRef,p_Djac(:,:,1:IELmax-IELset+1),p_Ddetj(:,1:IELmax-IELset+1))
             
      END IF
      
      ! Calculate the values of the basis functions.
      ! Pass p_DcubPts as point coordinates, which point either to the
      ! coordinates on the reference element (the same for all elements)
      ! or on the real element - depending on whether this is a 
      ! parametric or nonparametric element.
      CALL elem_generic_sim (p_relementDistribution%itrialElement, p_Dcoords, &
            p_Djac(:,:,1:IELmax-IELset+1), p_Ddetj(:,1:IELmax-IELset+1), &
            Bder, Dbas, ncubp, IELmax-IELset+1, p_DcubPts)
            
      ! We want to set up the nonlinear part of the matrix
      !
      !   n~_h (u_h, u_h, v_h) 
      !
      ! = n_h (u_h, u_h, v_h) + sum_T ( delta_T ( u_h*grad u_h, u_h*grad v_h)_T )
      !   ^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !  standard nonlin. part                  stabilization
      !
      ! More precisely, as we want to assemble the matrix which is 
      ! later multiplied with coefficient vectors, we have to insert
      ! basis functions in the above terms instead of u_h and v_h.
      ! Assuming the representation u_h=sum_j(u_j*Phi_j) and 
      ! v_h=sum_i(u_i,Phi_i), the above term is evaluated in the
      ! DOF's as:
      ! 
      !   n_h (u_h, Phi_j, Phi_i) 
      ! + sum_T ( delta_T ( u_h*grad Phi_j, u_h*grad Phi_i )_T )
      !
      ! In nonstationary simulations, the system matrix typically
      ! contains a mass matrix to respect the time derivative.
      ! The matrix has the form
      !
      ! [  dcmass*M*I  +  THWEIG * (-nu * Laplace(.))  ] + THWEIG * u grad(.)
      !
      ! In a first step, we calculate the velocity field in all
      ! cubature points on all elements of the current block.
      ! If we only have a primary velocity field
      ! (dweight2=0), we can calculate that only by summing up the
      ! velocities in U1Lx, otherwise we have to sum up
      ! dweight1*u1vel + dweight2*u2vel
      
      IF (dweight2 .EQ. 0) THEN
      
        ! Loop over all elements in the current set
        DO IEL=1,IELmax-IELset+1
        
          ! Loop over all cubature points on the current element
          DO ICUBP = 1, ncubp
          
            du1loc = 0.0_DP
            du2loc = 0.0_DP
          
            ! Perform a loop through the trial DOF's.
            DO JDOFE=1,indof

              ! Get the value of the (test) basis function 
              ! phi_i (our "O") in the cubature point:
              db = Dbas(JDOFE,1,ICUBP,IEL)
              
              ! Sum up to the value in the cubature point
              JDFG = Idofs(JDOFE,IEL)
              du1loc = du1loc + u1Xvel(JDFG)*db
              du2loc = du2loc + u1Yvel(JDFG)*db

            END DO ! JDOFE
            
            ! Save the computed velocity
            Dvelocity(1,ICUBP,IEL) = du1loc
            Dvelocity(2,ICUBP,IEL) = du2loc
          
          END DO ! ICUBP
          
        END DO ! IEL
          
      ELSE

        DO IEL=1,IELmax-IELset+1
        
          ! Loop over all cubature points on the current element
          DO ICUBP = 1, ncubp
          
            du1loc = 0.0_DP
            du2loc = 0.0_DP
          
            ! Perform a loop through the trial DOF's.
            DO JDOFE=1,indof

              ! Get the value of the (trial) basis function 
              ! phi_i in the cubature point:
              db = Dbas(JDOFE,1,ICUBP,IEL)

              ! Sum up to the value in the cubature point
              JDFG = Idofs(JDOFE,IEL)
              du1loc = du1loc + (dweight1*u1Xvel(JDFG) + dweight2*u2Xvel(JDFG))*db
              du2loc = du2loc + (dweight1*u1Yvel(JDFG) + dweight2*u2Yvel(JDFG))*db

            END DO ! JDOFE
            
            ! Save the computed velocity
            Dvelocity(1,ICUBP,IEL) = du1loc
            Dvelocity(2,ICUBP,IEL) = du2loc
          
          END DO ! ICUBP
          
        END DO ! IEL
      
      END IF
      
      ! If ALE is not active, calculate 
      !
      !     U * grad(Phi_j)  =  < grad(Phi_j), U >
      !
      !   = ( grad(Phi_j)_1 , (DU1) )
      !     ( grad(Phi_j)_2   (DU2) )
      !
      ! If ALE is active, use v=mesh velocity and calculate 
      !
      !       (U-v) * grad(Phi_j)  =  < grad(Phi_j), U-v >
      !
      !     = ( grad(Phi_j)_1 , (DU1-v) )
      !       ( grad(Phi_j)_2   (DU2-v) )
      !
      ! That means, we have to modify Dvelocity in that way that
      ! we have to substract the mesh velocity field in the cubature
      ! points.
      
      IF (bALE) THEN
        
        ! Calculate the values of the basis functions in all the points
        ! on all the elements
        CALL elem_generic_sim (EL_Q1, p_Dcoords, &
              p_Djac(:,:,1:IELmax-IELset+1), p_Ddetj(:,1:IELmax-IELset+1), &
              Bder, DbasALE, ncubp, IELmax-IELset+1, p_DcubPtsRef)
        
        ! Loop over all elements in the current set
        DO IEL=1,IELmax-IELset+1
        
          ! Loop over all cubature points on the current element
          DO ICUBP = 1, ncubp
          
            du1loc = 0.0_DP
            du2loc = 0.0_DP
          
            ! Perform a loop through the trial DOF's.
            DO JDOFE=1,indof

              ! Get the value of the (trial) basis function 
              db= Dbas(JDOFE,1,ICUBP,IEL)
              
              ! Sum up to the value in the cubature point
              JDFG = IdofsALE(IDOFE,IEL)
              du1loc = du1loc + DmeshVelocity(1,JDFG)*db
              du2loc = du2loc + DmeshVelocity(2,JDFG)*db

            END DO ! JDOFE
              
            ! Save the computed velocity
            Dvelocity(1,ICUBP,IEL) = Dvelocity(1,ICUBP,IEL) - du1loc
            Dvelocity(2,ICUBP,IEL) = Dvelocity(2,ICUBP,IEL) - du2loc
          
          END DO ! ICUBP
          
        END DO ! IEL
        
      END IF
      
      ! Ok, we now use Dvelocity as coefficient array in the assembly
      ! of a biliinear form!
      !
      ! Loop over the elements in the current set.

      DO IEL=1,IELmax-IELset+1
        
        ! Clear the local matrix
        Dentry = 0.0_DP
        
        ! Loop over all cubature points on the current element
        DO ICUBP = 1, ncubp

          ! Calculate the current weighting factor in the cubature formula
          ! in that cubature point.

          OM = Domega(ICUBP)*p_Ddetj(ICUBP,IEL)

          ! Current velocity in this cubature point:
          du1loc = Dvelocity (1,ICUBP,IEL)
          du2loc = Dvelocity (2,ICUBP,IEL)
          
          ! We take a more detailed look onto the last scalar product
          ! of n~_h (u_h, u_h, v_h) what we want to calculate here.
          !
          ! The vector u_h=(DU1,DU2) contains both velocity components,
          ! for the X as well as for the Y velocity. On the other hand
          ! the system matrix we want to build here will be designed for 
          ! one velocity component only! Therefore, Phi_i and Phi_j
          ! are scalar functions, so grad(Phi_i), grad(Phi_j) are vectors
          ! with two components. Therefore, the last scalar product is more 
          ! in detail:
          !
          !     ( u_h*grad Phi_j, u_h*grad Phi_i )_T
          !
          ! =   ( < (DU1) , (grad(Phi_j)_1) > , < (DU1) , (grad(Phi_i)_1) > )_T
          !         (DU2) , (grad(Phi_j)_2)       (DU2) , (grad(Phi_i)_2)  
          !
          ! =   < (DU1) , (grad(Phi_j)_1) >  *  < (DU1) , (grad(Phi_j)_1) >
          !       (DU2) , (grad(Phi_j)_2)         (DU2) , (grad(Phi_j)_2)
          !
          ! =   HSUMJ * HSUMI
          !
          ! i.e. a product of two scalar values!
          !
          ! Summing up over all pairs of multiindices.
          !
          ! Outer loop over the DOF's i=1..indof on our current element, 
          ! which corresponds to the basis functions Phi_i:

          DO IDOFE=1,indof
          
            ! Fetch the contributions of the (test) basis functions Phi_i
            ! (our "O")  for function value and first derivatives for the 
            ! current DOF into HBASIy:
          
            HBASI1 = Dbas(IDOFE,1,ICUBP,IEL)
            HBASI2 = Dbas(IDOFE,2,ICUBP,IEL)
            HBASI3 = Dbas(IDOFE,3,ICUBP,IEL)
           
            ! Calculate 
            !
            !     U * grad(Phi_i)  =  < grad(Phi_i), U >
            !
            !   = ( grad(Phi_i)_1 , (DU1) )
            !     ( grad(Phi_i)_2   (DU2) )
            !
            ! Remember: DU1MV=DU2MV=0 in this case.
            !
            ! If ALE is active, use v=mesh velocity and calculate 
            !
            !     (U-v) * grad(Phi_i)  =  < grad(Phi_i), U-v >
            !
            !   = ( grad(Phi_i)_1 , (DU1-DU1MV) )
            !     ( grad(Phi_i)_2   (DU2-DU2MV) )

            HSUMI = HBASI2*du1loc + HBASI3*du2loc

            ! Inner loop over the DOF's j=1..indof, which corresponds to
            ! the basis function Phi_j:

            DO JDOFE=1,indof
              
              IF (IDOFE.EQ.JDOFE) THEN
              
                ! Short version of the evaluation of the matrix
                ! contribution - see below for a more detailed
                ! description what is added together here!
              
                AH = HSUMI*(DlocalDelta(IEL)*HSUMI+HBASI1) &
                    + dny*(HBASI2**2+HBASI3**2) &
                    + dalpha*HBASI1**2
    
              ELSE
              
                ! Fetch the contributions of the (trial) basis function Phi_j
                ! (out "X") for function value and first derivatives for the 
                ! current DOF into HBASJy:
              
                HBASJ1 = Dbas(JDOFE,1,ICUBP,IEL)
                HBASJ2 = Dbas(JDOFE,2,ICUBP,IEL)
                HBASJ3 = Dbas(JDOFE,3,ICUBP,IEL)

                ! Calculate 
                !
                !     U * grad(Phi_j)  =  < grad(Phi_j), U >
                !
                !   = ( grad(Phi_j)_1 , (DU1) )
                !     ( grad(Phi_j)_2   (DU2) )
                !
                ! Remember: DU1MV=DU2MV=0 in this case.
                !
                ! If ALE is active, use v=mesh velocity and calculate 
                !
                !     (U-v) * grad(Phi_j)  =  < grad(Phi_j), U-v >
                !
                !   = ( grad(Phi_j)_1 , (DU1-DU1MV) )
                !     ( grad(Phi_j)_2   (DU2-DU2MV) )
                !
                ! But as v is already incorporated into DVelocity,
                ! we don't have to worry about that.

                HSUMJ = HBASJ2*du1loc+HBASJ3*du2loc
    
                ! Finally calculate the contribution to the system
                ! matrix. Depending on the configuration of DNU,
                ! dalpha,ddelta,... this decomposes into three
                ! different parts:
                !
                ! AH = n~_h(u_h,phi_j,phi_i)        | nonlinear part
                !    + dny*(grad(phi_j,grad(phi_i)) | -dny*Laplace(u) = -dbeta*Stokes
                !    + dalpha*(phi_j*phi_i)         | Mass matrix
                !
                ! The last two parts are probably not added to the
                ! matrix by setting DNY or CT0 to 0, respectively.
                !
                ! For saving some numerical operations, we write:
                !
                !     HSUMJ * (Delta * HSUMI + HBASI1)
                !
                ! =   Delta * HSUMJ * HSUMI
                !   + HSUMJ * HBASI1
                !
                ! =   Delta * ( U*grad(Phi_j), U*grad(Phi_i) )
                !   + (U*grad(Phi_j),Phi_i)
                !
                ! <->   sum_T ( delta_T ( u_h*grad Phi_j, u_h*grad Phi_i) )_T
                !     + n_h (u_h, Phi_j, Phi_i)
                !
                ! plus the terms for the Stokes and Mass matrix,
                ! if their coefficient is <> 0.
                
                AH = HSUMJ*(DlocalDelta(IEL)*HSUMI+HBASI1) &
                    + dny*(HBASI2*HBASJ2+HBASI3*HBASJ3) &
                    + dalpha*HBASI1*HBASJ1
    
              ENDIF ! (IDOFE.EQ.JDOFE)

              ! Weighten the calculated value AH by the cubature
              ! weight OM and add it to the local matrix. After the
              ! loop over all DOF's is finished, each entry contains
              ! the calculated integral.

              Dentry(JDOFE,IDOFE) = Dentry(JDOFE,IDOFE)+OM*AH
              
            END DO ! IDOFE
            
          END DO ! JDOFE

        END DO ! ICUBP 
        
        ! Now we have set up a "local" system matrix. We can either    
        ! include it into the real matrix or we can use it to simply   
        ! modify the RHS vector to create a defect vector (throwing    
        ! away the information about the matrix afterwards, which would
        ! result in a matrix free modification of the RHS vector).     
        !
        ! For cdef= containing CONV_MODMATRIX, incorporate our "local" system matrix
        ! into the global matrix. The position of each entry DENTRY(X,Y)    
        ! in the global matrix array A was saved in element Kentry(X,Y)
        ! before.                                                      
        ! Kentry gives the position of the additive contributions in Dentry.
        ! The entry is weighted by the current dtheta, which is usually
        ! the weighting parameter of the corresponding THETA-scheme of a
        ! nonstationary simulation. For stationary simulations, dtheta is typically
        ! 1.0 which includes the local matrix into the global one directly.)
        
        IF (IAND(cdef,CONV_MODMATRIX) .NE. 0) THEN
          DO IDOFE=1,indof
            DO JDOFE=1,indof
              p_DA(Kentry(JDOFE,IDOFE,IEL)) = p_DA(Kentry(JDOFE,IDOFE,IEL)) + &
                dtheta * Dentry(JDOFE,IDOFE)
            END DO
          END DO
        END IF
        
        ! For cdef containing CONV_MODDEFECT, build the defect vector                     
        !     D = RHS - A*U                                         
        ! This is done matrix free, only with the help of the local 
        ! matrix.                                                   
        ! In this case, D=(D1,D2) is expected to be the RHS on      
        ! entry and will be updated to be the defect vector when    
        ! this routine is left.                                     
        
        IF (IAND(cdef,CONV_MODDEFECT) .NE. 0) THEN
          DO IDOFE=1,indof

            IDFG=Idofs(IDOFE,IEL)

            DO JDOFE=1,indof

              denth = dtheta*Dentry(JDOFE,IDOFE)         
    
              JDFG=Idofs(JDOFE,IEL)
              Ddef1(IDFG)= Ddef1(IDFG) - denth*Du1(JDFG)
              Ddef2(IDFG)= Ddef2(IDFG) - denth*Du2(JDFG)

            END DO
          END DO
        END IF
        
      END DO ! IEL

    END DO ! IELset
    
    ! Release memory
    CALL domint_doneIntegration(rintSubset)

    DEALLOCATE(DlocalDelta)
    DEALLOCATE(Dvelocity)
    DEALLOCATE(Dentry)
    DEALLOCATE(Kentry)
    DEALLOCATE(IdofsALE)
    DEALLOCATE(Idofs)
    DEALLOCATE(DbasALE)
    DEALLOCATE(Dbas)

  CONTAINS
  
    SUBROUTINE getLocalDelta (U1L1,U1L2,U2L1,U2L2,A1L,A2L,IEL,duMaxR,ddelta, &
                        Kvert,Kmid,Dcorvg,KDFG,IDFL,UPSAM,NUREC)

    ! This routine calculates a local ddelta=DELTA_T for a finite element
    ! T=IEL. This can be used by the streamline diffusion stabilization
    ! technique as a multiplier of the (local) bilinear form.
    !
    ! The effective velocity that is used for calculating the ddelta
    ! is combined by a weighted mean of the two velocity fields U1,U2
    ! by:
    !                   Ux = A1*U1Lx + A2*U2Lx
    ! The coefficients A1,A2 allow the caller to take influence on which
    ! velocity field to weight more.
    
    ! Main velocity field.
    REAL(DP), DIMENSION(*), INTENT(IN) :: U1L1,U1L2
    
    ! Secondary velocity field. 
    REAL(DP), DIMENSION(*), INTENT(IN) :: U2L1,U2L2
    
    ! weighting factor for U1L1/U1L2
    REAL(DP), INTENT(IN) :: A1L
    
    ! weighting factor for U2L1/U2L2
    REAL(DP), INTENT(IN) :: A2L
    
    ! Reciprocal of the maximum norm of velocity in the domain:
    ! 1/duMaxR = 1/||u||_Omega
    REAL(DP), INTENT(IN) :: duMaxR
    
    ! Reciprocal value 1/NU of coefficient NU in front of the
    ! Laplacian term of the Navier-Stokes equation
    !   NU * Laplace(u) + u*grad(u) + ...
    REAL(DP), INTENT(IN) :: NUREC
    
    ! user defined parameter for configuring the streamline diffusion.
    ! < 0: Simple calculation of ddelta, using 
    !      ddelta = |UPSAM| * h_T.
    ! > 0: usually UPSAM = 0.1 .. 2; Samarskji-like calculation of ddelta using:
    !      ddelta = UPSAM * h_t/||u||_T * 2*Re_T/(1+Re_T)
    REAL(DP), INTENT(IN) :: UPSAM
    
    ! Element where the ddelta should be calculated
    INTEGER(PREC_ELEMENTIDX) :: IEL
    
    ! Number of degrees of freedom on element IEL
    INTEGER(PREC_DOFIDX) :: IDFL
    
    ! Array with global degrees of freedom, corresponding to
    ! local degrees of freedom 1..IDFL on element IEL.
    INTEGER(PREC_DOFIDX), DIMENSION(*) :: KDFG
    
    INTEGER(PREC_POINTIDX), DIMENSION(TRIA_MAXNVE2D,*) :: Kvert
    INTEGER(PREC_POINTIDX), DIMENSION(TRIA_MAXNME2D,*) :: Kmid
    REAL(DP), DIMENSION(NDIM2D,*), INTENT(IN) :: Dcorvg

    ! local ddelta
    REAL(DP), INTENT(OUT) :: ddelta

    ! local variables
    REAL(DP) :: dlocalH,DU1,DU2,dunorm,RELOC
    INTEGER(PREC_DOFIDX) :: idof

      ! Loop through the local degrees of freedom on element IEL.
      ! Sum up the velocities on these DOF's. This will result
      ! in the vector (DU1,DU2) representing the (mean) X/Y-velocity
      ! through element IEL.
  
      ! For elements whose DOF's represent directly the velocity, U1/U2 
      ! represent the mean velocity
      ! along an egde/on the midpoint of each edge, so U1/U2 is
      ! clearly an approximation to the velocity in element T.

      DU1=0.0_DP
      DU2=0.0_DP
      DO idof=1,IDFL
        DU1=DU1+(A1L*U1L1(KDFG(idof))+A2L*U2L1(KDFG(idof)))
        DU2=DU2+(A1L*U1L2(KDFG(idof))+A2L*U2L2(KDFG(idof)))
      END DO

      ! Calculate the norm of that local velocity:

      dunorm = SQRT(DU1**2+DU2**2) / DBLE(IDFL)
      
      ! Now we have:   dunorm = ||u||_T
      ! and:           u_T = a1*u1_T + a2*u2_T
  
      ! If the norm of the velocity is small, we choose ddelta = 0,
      ! which results in central difference in the streamline diffusion
      ! matrix assembling:

      IF (dunorm.LE.1D-8) THEN
      
        ddelta = 0.0_DP

      ELSE

        ! u_T defines the "slope" of the velocity through
        ! the element T. At next, calculate the local mesh width
        ! dlocalH = h = h_T on our element T=IEL:

        CALL getLocalMeshWidth (dlocalH,dunorm, DU1, DU2, IEL, Kvert,Kmid,Dcorvg)

        ! Calculate ddelta... (cf. p. 121 in Turek's CFD book)

        IF (UPSAM.LT.0.0_DP) THEN

          ! For UPSAM<0, we use simple calculation of ddelta:        
        
          ddelta = ABS(UPSAM)*dlocalH
          
        ELSE
        
          ! For ddelta >= 0, we use standard Samarskji-like calculation
          ! of ddelta. At first calculate the local Reynolds number
          ! RELOC = Re_T = ||u||_T * h_T / NU
          
          RELOC = dunorm*dlocalH*NUREC
          
          ! and then the ddelta = UPSAM * h_t/||u|| * 2*Re_T/(1+Re_T)
          
          ddelta = UPSAM * dlocalH*duMaxR * 2.0_DP*(RELOC/(1.0_DP+RELOC))
          
        ENDIF ! (UPSAM.LT.0.0)
        
      END IF ! (dunorm.LE.1D-8)

    END SUBROUTINE

    ! ----------------------------------------------------------------------

    SUBROUTINE getLocalMeshWidth (dlocalH, dunorm,  XBETA1, &
                       XBETA2, JEL,Kvert,Kmid,Dcorvg)
    
    ! Determine the local mesh width for an element JEL of a 
    ! triangulation.
    
    ! Element where the local h should be calculated
    INTEGER(PREC_ELEMENTIDX), INTENT(IN)               :: JEL
    
    INTEGER(PREC_POINTIDX), DIMENSION(TRIA_MAXNVE2D,*) :: Kvert
    INTEGER(PREC_EDGEIDX), DIMENSION(TRIA_MAXNVE2D,*)  :: Kmid
    REAL(DP), DIMENSION(NDIM2D,*), INTENT(IN)          :: Dcorvg
    
    ! norm ||u||_T = mean velocity through element T=JEL
    REAL(DP), INTENT(IN)  :: dunorm
    
    ! mean velocity u_T = (xbeta1,xbeta2) through element T=JEL
    REAL(DP), INTENT(IN)  :: XBETA1, XBETA2
    
    ! local mesh width
    REAL(DP), INTENT(OUT) :: dlocalH
    
    ! local variables
    REAL(DP) :: dlambda
    INTEGER(PREC_POINTIDX) :: NECK1,NECK2,NECK3,NECK4
    REAL(DP) :: X1,Y1,X2,Y2,X3,Y3,X4,Y4
    REAL(DP) :: dalphaMax, dalpha

      ! Fetch the numbers of the four corners of element JEL

      neck1=Kvert(1,JEL)
      neck2=Kvert(2,JEL)
      neck3=Kvert(3,JEL)
      neck4=Kvert(4,JEL)

      ! Fetch the coordinates of these corners

      x1=Dcorvg(1,neck1)
      y1=Dcorvg(2,neck1)
      x2=Dcorvg(1,neck2)
      y2=Dcorvg(2,neck2)
      x3=Dcorvg(1,neck3)
      y3=Dcorvg(2,neck3)
      x4=Dcorvg(1,neck4)
      y4=Dcorvg(2,neck4)

      ! Scale: (deactivated)

      !  dsp=max(xbeta1,xbeta2)

      !  xbeta1=xbeta1
      !  xbeta2=xbeta2

      dalphaMax=0.0_DP
      
      ! Loop through the four corners of element JEL and check
      ! of a line with slope BETA=(xbeta1,xbeta2) starting in this
      ! corner really intersects with one of the edges of the element.
      ! Remark that we only have to check the two opposite edges
      ! to the current corner!

      ! -----------------------------------------------------------------
      ! Check the first corner:

      CALL intersectLines2D(X1,Y1,dalpha,XBETA1,XBETA2, &
                  X3,Y3,dlambda,X2,Y2)
      dalphaMax=MAX(dalpha,dalphaMax)

      CALL intersectLines2D(X1,Y1,dalpha,XBETA1,XBETA2, &
                  X3,Y3,dlambda,X4,Y4)
      dalphaMax=MAX(dalpha,dalphaMax)
      
      ! -----------------------------------------------------------------
      ! The second one...
      
      CALL intersectLines2D(X2,Y2,dalpha,XBETA1,XBETA2, &
                  X4,Y4,dlambda,X1,Y1)
      dalphaMax=MAX(dalpha,dalphaMax)

      CALL intersectLines2D(X2,Y2,dalpha,XBETA1,XBETA2, &
                  X4,Y4,dlambda,X3,Y3)
      dalphaMax=MAX(dalpha,dalphaMax)
      
      ! -----------------------------------------------------------------
      ! The third one...
      
      CALL intersectLines2D(X3,Y3,dalpha,XBETA1,XBETA2, &
                  X1,Y1,dlambda,X2,Y2)
      dalphaMax=MAX(dalpha,dalphaMax)

      CALL intersectLines2D(X3,Y3,dalpha,XBETA1,XBETA2, &
                  X1,Y1,dlambda,X4,Y4)
      dalphaMax=MAX(dalpha,dalphaMax)
      
      ! -----------------------------------------------------------------
      ! And the fourth=last one...
      
      CALL intersectLines2D(X4,Y4,dalpha,XBETA1,XBETA2, &
                  X2,Y2,dlambda,X1,Y1)
      dalphaMax=MAX(dalpha,dalphaMax)

      CALL intersectLines2D(X4,Y4,dalpha,XBETA1,XBETA2, &
                  X2,Y2,dlambda,X3,Y3)
      dalphaMax=MAX(dalpha,dalphaMax)

      ! -----------------------------------------------------------------
      ! finally determine the local h=h_T

      dlocalH=dalphaMax*4.0_DP*dunorm

    END SUBROUTINE
    
    ! ----------------------------------------------------------------------

    SUBROUTINE intersectLines2D (XO,YO,dalpha,BETA1,BETA2, &
                       XA,YA,dlambda,XB,YB)

    ! Intersect two lines in R^2

    ! Origin of line 1
    REAL(DP), INTENT(IN) :: XO,YO
    
    ! Direction of line 1
    REAL(DP), INTENT(IN) :: BETA1,BETA2
    
    ! One point on the second line
    REAL(DP), INTENT(IN) :: XA,YA
    
    ! Another point on the second line
    REAL(DP), INTENT(IN) :: XB,YB
    
    ! Parameter value of the intersection point on line 1.
    ! =0.0, if there is no intersection point
    REAL(DP), INTENT(OUT) :: dalpha
    
    REAL(DP), INTENT(OUT) :: dlambda
    
    ! local variables
    DOUBLE PRECISION :: dsp

      dsp=BETA2*(XB-XA)-BETA1*(YB-YA)
      
      IF (dsp.eq.0.0_DP) THEN
      
        ! beta and the vector are parallel
        dalpha=0.0_DP
        
      ELSE  

        dlambda=(BETA1*(YA-YO)-BETA2*(XA-XO))/dsp

        ! is the intersection point inside of the element?
        IF ((dlambda.GE.-1E-1_DP).AND.(dlambda.LE.1.11E0_DP)) THEN
          IF (BETA1 .NE. 0.0_DP) THEN
            dalpha=((XA-XO)+dlambda*(XB-XA))/BETA1
          ELSE
            IF (BETA2 .NE. 0.0_DP) THEN
              dalpha=((YA-YO)+dlambda*(YB-YA))/BETA2
            ELSE
              dalpha=0.0_DP
            ENDIF
          ENDIF
        ELSE
          dalpha=0.0_DP
        ENDIF
        
      ENDIF

    END SUBROUTINE

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE conv_ueoJumpStabil_double_uni (rmatrixScalar,rconfig)
  
!<description>
  ! Unified edge oriented jump stabilisation.
  !
  ! Adds the unified edge oriented jump stabilisation to the matrix rmatrix:
  ! $$< Ju,v > = \sum_E \max(\gamma^{*}*\nu*h_E, \gamma h_E^2) 
  !            \int_E [grad u] [grad v] ds$$
  !
  ! Uniform discretisation, double precision structure-7 and 9 matrix.
  !
  ! For a rerefence about the stabilisation, see
  ! [Ouazzi, A.; Finite Element Simulation of Nonlinear Fluids, Application
  ! to Granular Material and Powder; Shaker Verlag, ISBN 3-8322-5201-0, p. 55ff]
  !
  ! WARNING: For edge oriented stabilisation, the underlying matrix rmatrix
  !   must have an extended stencil! The matrix structure must be set up with
  !   the BILF_MATC_EDGEBASED switch!!!
!</description>

!<input>
  ! Configuration block for the streamline diffusion scheme
  TYPE(t_jumpStabilisation), INTENT(IN) :: rconfig
!</input>

!<inputoutput>
  ! The matrix which is to be modified. The jump stabilisation is added
  ! to that matrix. 
  TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrixScalar
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER(PREC_DOFIDX) :: irow, jcol, idof
  INTEGER(PREC_EDGEIDX) :: IMT
  INTEGER(PREC_POINTIDX) :: ivt,ivt1,ivt2,NVT
  INTEGER(PREC_ELEMENTIDX) :: IEL
  LOGICAL :: bIdenticalTrialAndTest
  INTEGER :: IELcount,IDOFE, JDOFE, i, NVE, iedge
  REAL(DP) :: dval,dedgelength,dedgeweight,dphidx,dphidy,dpsidx,dpsidy,dcoeff
  
  ! Pointer to KLD, KCOL, DA
  INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld
  INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol
  REAL(DP), DIMENSION(:), POINTER :: p_Da
  
  ! An allocateable array accepting the DOF's of a set of elements.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsTestTempl
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsTrialTempl
  INTEGER(PREC_DOFIDX), DIMENSION(EL_MAXNBAS*2), TARGET :: IdofsTest
  INTEGER(PREC_DOFIDX), DIMENSION(EL_MAXNBAS*2), TARGET :: IdofsTrial
  INTEGER(PREC_DOFIDX), DIMENSION(:), POINTER :: p_IdofsTrial
  
  ! Arrays saving the local DOF numbers belonging to the global
  ! DOF numbers in IdofsTest/IdofsTrial  
  INTEGER, DIMENSION(EL_MAXNBAS*2),TARGET :: IlocalDofsTest,IlocalDofsTrial
  INTEGER(PREC_DOFIDX), DIMENSION(:), POINTER :: p_IlocalDofsTrial
  
  ! Renumbering strategy for local DOF's
  INTEGER, DIMENSION(EL_MAXNBAS), TARGET :: IlocalDofTestRenum,IlocalDofTrialRenum
  INTEGER, DIMENSION(:), POINTER :: p_IlocalDofTrialRenum
  
  ! Number of local DOF's on the patch
  INTEGER :: ndofTest,ndofTrial
  
  ! Number of local degees of freedom for trial and test functions
  INTEGER :: indofTrialPerElement, indofTestPerElement
  
  ! The triangulation structure - to shorten some things...
  TYPE(t_triangulation), POINTER :: p_rtriangulation
  
  ! Some triangulation arrays we need frequently
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IneighboursAtElement
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IelementsAtEdge
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
  INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
  INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtEdge
  REAL(DP), DIMENSION(:,:), POINTER :: p_DcornerCoordinates
  
  ! Current element distribution
  TYPE(t_elementDistribution), POINTER :: p_elementDistribution
  
  ! Underlying discretisation structure
  TYPE(t_spatialDiscretisation), POINTER :: p_rdiscretisation

  ! Arrays for saving Jacobian determinants and matrices
  REAL(DP), DIMENSION(:,:), POINTER :: p_Ddetj
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_Djac

  ! Allocateable arrays for the values of the basis functions - 
  ! for test and trial spaces.
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: DbasTest,DbasTrial
  REAL(DP), DIMENSION(:,:,:,:), POINTER :: p_DbasTrial

  ! Local matrices, used during the assembly.
  ! Values and positions of values in the global matrix.
  INTEGER(PREC_DOFIDX), DIMENSION(:), ALLOCATABLE :: Kentry
  REAL(DP), DIMENSION(:), ALLOCATABLE :: Dentry

  ! A t_domainIntSubset structure that is used for storing information
  ! and passing it to callback routines.
  TYPE(t_domainIntSubset) :: rintSubset
  
  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: DcubPtsRefOnAllEdges

  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPtsRef

  ! An array receiving the coordinates of cubature points on
  ! the real element for all elements in a set.
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPtsReal
  
  ! Pointer to the cubature points for trial and test space.
  ! These point either to p_DcubPtsRef or p_DcubPtsReal, depending on
  ! whether the element is parametric or not.
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPtsTrial
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPtsTest
  
  ! Array with coordinates of the corners that form the real element.
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_Dcoords
  
  ! Cubature point weights
  REAL(DP), DIMENSION(CUB_MAXCUBP) :: Domega
  
  ! Cubature point coordinates on the reference element
  REAL(DP), DIMENSION(CUB_MAXCUBP, NDIM3D) :: Dxi1D,Dxi2D
  
  ! number of cubature points on the reference element
  INTEGER :: ncubp,icubp
  
  ! Test/trial functions nonparameric?
  LOGICAL :: bnonparTrial,bnonparTest
  
  ! Derivative specifiers
  LOGICAL, DIMENSION(EL_MAXNDER) :: BderTrial, BderTest

    ! Get a pointer to the triangulation and discretisation.
    p_rtriangulation => rmatrixScalar%p_rspatialDiscretisation%p_rtriangulation
    p_rdiscretisation => rmatrixScalar%p_rspatialDiscretisation
    
    ! Get Kvert, Kadj, Kmid, Kmel,...
    CALL storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
                                p_IverticesAtElement)
    CALL storage_getbase_int2d (p_rtriangulation%h_IneighboursAtElement,&
                                p_IneighboursAtElement)
    CALL storage_getbase_int2d (p_rtriangulation%h_IedgesAtElement,&
                                p_IedgesAtElement)
    CALL storage_getbase_int2d (p_rtriangulation%h_IelementsAtEdge,&
                                p_IelementsAtEdge)
    CALL storage_getbase_int2d (p_rtriangulation%h_IverticesAtEdge,&
                                p_IverticesAtEdge)
    CALL storage_getbase_double2d (p_rtriangulation%h_DcornerCoordinates,&
                                  p_DcornerCoordinates)
                                
    NVT = p_rtriangulation%NVT
    
    ! Get KLD, KCol...
    CALL storage_getbase_int (rmatrixScalar%h_KLD,p_KLD)
    CALL storage_getbase_int (rmatrixScalar%h_Kcol,p_Kcol)
    CALL storage_getbase_double (rmatrixScalar%h_Da,p_Da)
    
    ! Activate the one and only element distribution
    p_elementDistribution => p_rdiscretisation%RelementDistribution(1)

    ! Get the number of local DOF's for trial and test functions
    indofTrialPerElement = elem_igetNDofLoc(p_elementDistribution%itrialElement)
    indofTestPerElement = elem_igetNDofLoc(p_elementDistribution%itestElement)
    
    ! Triangle elements? Quad elements?
    NVE = elem_igetNVE(p_elementDistribution%itrialElement)
    
    ! Get the number of corner vertices of the element
    IF (NVE .NE. elem_igetNVE(p_elementDistribution%itestElement)) THEN
      PRINT *,'conv_ueoJumpStabil2d_double_uni: element spaces incompatible!'
      STOP
    END IF
    
    ! Initialise the cubature formula,
    ! Get cubature weights and point coordinates on the reference line [-1,1]
    CALL cub_getCubPoints(rconfig%ccubType, ncubp, Dxi1D, Domega)

    ! Allocate arrays accepting cubature point coordinates.
    ! We have ncubp cubature points on every egde.
    ! DcubPtsRef saves all possible combinations, which edge of
    ! one element might interact with one edge of one neighbour
    ! element.
    ALLOCATE(DcubPtsRefOnAllEdges(NDIM2D,ncubp,NVE))
    
    ! Put the 1D cubature points from to all of the edges on the
    ! reference element, so we can access them quickly.
    DO iedge = 1,NVE
      CALL trafo_mapCubPts1Dto2DRefQuad(iedge, ncubp, Dxi1D, Dxi2D)
      DO i=1,ncubp
        DcubPtsRefOnAllEdges(1,i,iedge) = Dxi2D(i,1)
        DcubPtsRefOnAllEdges(2,i,iedge) = Dxi2D(i,2)
      END DO
    END DO
    
    ! Get from the trial element space the type of coordinate system
    ! that is used there:
    i = elem_igetCoordSystem(p_elementDistribution%itrialElement)
    
    ! Allocate memory and get local references to it.
    CALL domint_initIntegration (rintSubset,2,ncubp,i,NDIM2D)
    p_DcubPtsRef =>  rintSubset%p_DcubPtsRef
    p_DcubPtsReal => rintSubset%p_DcubPtsReal
    p_Djac =>        rintSubset%p_Djac
    p_Ddetj =>       rintSubset%p_Ddetj
    p_Dcoords =>     rintSubset%p_DCoords
    
    ! Check if one of the trial/test elements is nonparametric
    bnonparTrial = elem_isNonparametric(p_elementDistribution%itrialElement)
    bnonparTest  = elem_isNonparametric(p_elementDistribution%itestElement)
                    
    ! Let p_DcubPtsTrial / p_DcubPtsTest point either to DcubPtsRealEval or
    ! DcubPtsRefEval - depending on whether the space is parametric or not.
    IF (bnonparTrial) THEN
      p_DcubPtsTrial => p_DcubPtsReal
    ELSE
      p_DcubPtsTrial => p_DcubPtsRef
    END IF
    
    IF (bnonparTest) THEN
      p_DcubPtsTest => p_DcubPtsReal
    ELSE
      p_DcubPtsTest => p_DcubPtsRef
    END IF
    
    ! Allocate arrays saving the local matrices for all elements
    ! in an element set. We allocate the arrays large enough...
    ALLOCATE(Kentry(indofTestPerElement*2*indofTrialPerElement*2))
    ALLOCATE(Dentry(indofTestPerElement*2*indofTrialPerElement*2))
    
    ! Allocate memory for obtaining DOF's:
    ALLOCATE(IdofsTrialTempl(indofTrialPerElement,2))
    ALLOCATE(IdofsTestTempl(indofTestPerElement,2))
    
    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  ALLOCATE(DbasTest(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock+1))
    !  ALLOCATE(DbasTrial(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock+1))
    ! would lead to nonused memory blocks in these arrays during the assembly, 
    ! which reduces the speed by 50%!
    !
    ! We allocate space for 3 instead of 2 elements. The reason is that
    ! we later permute the values of the basis functions to get
    ! a local numbering on the patch. That's also the reason, we allocate
    ! not indofTestPerElement elements, but even indofTestPerElement,
    ! which is more than enough space to hold the values of the DOF's of
    ! a whole element patch.
    
    ALLOCATE(DbasTest(indofTestPerElement*2, &
            elem_getMaxDerivative(p_elementDistribution%itestElement),&
            ncubp,3))
    ALLOCATE(DbasTrial(indofTrialPerElement*2,&
            elem_getMaxDerivative(p_elementDistribution%itrialElement), &
            ncubp,3))
             
    ! Test if trial/test functions are identical.
    ! We don't rely on bidenticalTrialAndTest purely, as this does not
    ! indicate whether there are identical trial and test functions
    ! in one block!
    bIdenticalTrialAndTest = &
      p_elementDistribution%itrialElement .EQ. p_elementDistribution%itestElement
      
    ! Let p_IdofsTrial point either to IdofsTrial or to the DOF's of the test
    ! space IdofTest (if both spaces are identical). 
    ! We create a pointer for the trial space and not for the test space to
    ! prevent pointer-arithmetic in the innerst loop below!
    ! The same for p_IlocalDofsTrial and the permutation array p_IlocalDofTrialRenum.
    ! Let p_DbasTrial point either to DbasTrial or DbasTest, depending on
    ! whether the spaces are identical.
    IF (bIdenticalTrialAndTest) THEN
      p_IdofsTrial => IdofsTest
      p_IlocalDofsTrial => IlocalDofsTest
      p_IlocalDofTrialRenum => IlocalDofTestRenum
      p_DbasTrial => DbasTest
    ELSE
      p_IdofsTrial => IdofsTrial
      p_IlocalDofsTrial => IlocalDofsTrial
      p_IlocalDofTrialRenum => IlocalDofTrialRenum
      p_DbasTrial => DbasTrial
    END IF
    
    ! Set up which derivatives to compute in the basis functions: X/Y-derivative
    BderTrial = .FALSE.
    BderTrial(DER_DERIV_X) = .TRUE.
    BderTrial(DER_DERIV_Y) = .TRUE.

    BderTest = .FALSE.
    BderTest(DER_DERIV_X) = .TRUE.
    BderTest(DER_DERIV_Y) = .TRUE.

    ! We have always 2 elements in each patch...
    IELcount = 2
      
    ! We loop through all edges
    DO IMT = 1,p_rtriangulation%NMT
    
      ! Check if we have 1 or 2 elements on the current edge
      IF (p_IelementsAtEdge (2,IMT) .EQ. 0) THEN
      
        ! This only happens if we have a boundary edge.
        ! The boundary handling is... doing nothing! So we can skip
        ! the handling of that edge completely!
        CYCLE
      
      END IF
      
      ! Fill the basis function arrays with 0. Essential, as only parts
      ! of them are overwritten later.
      DbasTest = 0.0_DP
      DbasTrial = 0.0_DP
      
      ! On an example, we now show the relationship between all the DOF's
      ! we have to consider in the current situation. We have two elements,
      ! let's say IEL1 and IEL2, with their local and global DOF's
      ! (example for Q1~):
      !
      !   local DOF's on the           corresponding global
      !   element and on the patch     DOF's
      !
      !    +----4----+----7----+       +----20---+----50---+
      !    |    4    |    4    |       |         |         |
      !    |         |         |       |         |         |
      !    1 1     3 3 1     3 6       10        60        70
      !    |         |         |       |         |         |
      !    |    2    |    2    |       |         |         |
      !    +----2----+----3----+       +----40---+----30---+
      !        IEL1      IEL2              IEL1      IEL2
      !
      !
      ! On every element, we have 4 local DOF's (1..4). On the other hand,
      ! we have "local DOF's on the patch" (1..7), ehich we call "patch DOF's"
      ! from now on. To every local DOF, there belongs a global DOF
      ! (10,20,30,... or whatever), which gives the coefficient of the basis
      ! function.
      !
      ! Numbering that way, the local DOF's of IEL1 obviously coincide 
      ! with the first couple of local DOF's of the element patch! Only
      ! the local DOF's of element IEL2 make trouble, as they have another
      ! numbering.
      
      ! Get the global DOF's of the 1 or two elements
      CALL dof_locGlobMapping_mult(p_rdiscretisation, &
                                  p_IelementsAtEdge (1:IELcount,IMT), &
                                  .TRUE., IdofsTestTempl)
                                   
      ! Some of the DOF's on element 2 may coincide with DOF's on element 1.
      ! More precisely, some must coincide! Therefore, we now have to collect the
      ! DOF's uniquely and to figure out, which local DOF's of element 2
      ! must renumbered to the appropriate local patch-DOF (like in the
      ! example, where local DOF 1 of element IEL2 must be renumbered
      ! to local patch DOF 3!
      !
      ! As the first couple of local DOF's of IEL1 coincide with the local
      ! DOF's of the patch, we can simply copy them:
      
      ndofTest = indofTestPerElement
      IdofsTest(1:ndofTest) = IdofsTestTempl(1:ndofTest,1)
      
      ! Furthermore, we read IdofsTestTempl and store the DOF's there in IdofsTest,
      ! skipping all DOF's we already have and setting up the renumbering strategy
      ! of local DOF's on element IEL2 to patch DOF's.
      
      skiploop: DO IDOFE = 1,indofTestPerElement
        
        ! Do we have the DOF? 
        idof = IdofsTestTempl(IDOFE,IELcount)
        DO JDOFE=1,ndofTest
          IF (IdofsTest(JDOFE) .EQ. idof) THEN
            ! Yes, we have it.
            ! That means, the local DOF idof has to be mapped to the
            ! local patch dof...
            IlocalDofTestRenum (IDOFE) = JDOFE
            
            ! Proceed with the next one.
            CYCLE skiploop
          END IF
        END DO
        
        ! We don't have that DOF! Append it to IdofsTest.
        ndofTest = ndofTest+1
        IdofsTest(ndofTest) = idof
        IlocalDofTestRenum (IDOFE) = ndofTest
        
        ! Save also the number of the local DOF.
        ! Note that the global DOF's in IdofsTestTempl(1..indofTestPerElement)
        ! belong to the local DOF's 1..indofTestPerElement -- in that order!
        IlocalDofsTest(ndofTest) = IDOFE
        
      END DO skiploop
        
      ! If the trial functions are different, get those DOF's
      IF (.NOT. bIdenticalTrialAndTest) THEN
        CALL dof_locGlobMapping_mult(p_rdiscretisation, &
                                    p_IelementsAtEdge (1:IELcount,IMT), &
                                    .FALSE., IdofsTrialTempl)

        ! Collect the DOF's uniquely in IdofsTrial.
        ! The DOF's of the first element are taken as they are.
        
        ndofTrial = indofTrialPerElement
        IdofsTrial(1:ndofTrial) = IdofsTrialTempl(1:ndofTrial,1)
        
        ! From the 2nd element on, we have to check which
        ! DOF's already appear in the first element.
        
        skiploop2: DO IDOFE = 1,indofTrialPerElement
          
          ! Do we have the DOF? 
          idof = IdofsTrialTempl(IDOFE,IELcount)
          DO JDOFE=1,ndofTrial
            IF (IdofsTrial(JDOFE) .EQ. idof) THEN
              ! Yes, we have it.
              ! That means, the local DOF idof has to be mapped to the
              ! local patch dof...
              IlocalDofTrialRenum (IDOFE) = JDOFE
              
              ! Proceed with the next one.
              CYCLE skiploop2
            END IF
          END DO
          
          ! We don't have that DOF! Append it to IdofsTrial.
          ndofTrial = ndofTrial+1
          IdofsTrial(ndofTrial) = idof
          IlocalDofTrialRenum (IDOFE) = ndofTrial
          
          ! Save also the number of the local DOF.
          ! Note that the global DOF's in IdofsTrialTempl(1..indofTrialPerElement)
          ! belong to the local DOF's 1..indofTrialPerElement -- in that order!
          IlocalDofsTrial(ndofTrial) = IDOFE
          
        END DO skiploop2
          
      ELSE
       
        ! Number of DOF's in test and trial space the same, as DOF's are the same.
        ndofTrial = ndofTest
        
      END IF
      
      ! Now we know: Our 'local' matrix (consisting of only these DOF's we just
      ! calculated) is a ndofsTest*ndofsTrial matrix.
      !
      ! Now extract the corresponding entries from the matrix.
      ! Kentry is an index where the entries of the local matrix Dentry
      ! (which we build up later) can be found in the global matrix.
      !
      ! The test space gives us the rows; loop through them
      
      DO IDOFE = 0,ndofTest-1
      
        ! The 'local' DOF IDOFE corresponds in the global matrix to line...
        irow = IdofsTest (1+IDOFE)
        
        ! Loop through that line to find the columns, indexed by the local
        ! DOF's in the trial space.
        trialspaceloop: DO JDOFE = 1,ndofTrial
          
          DO jcol = p_Kld(irow),p_Kld(irow+1)-1
            
            IF (p_Kcol(jcol) .EQ. p_IdofsTrial(JDOFE)) THEN
            
              ! Found! Put as destination pointer to Kentry
              Kentry (IDOFE*ndofTrial+JDOFE) = jcol
              
              ! Clear local matrix
              Dentry (IDOFE*ndofTrial+JDOFE) = 0.0_DP
              
              ! Jump out of the loop, proceed with next column
              CYCLE trialspaceloop
            
            END IF
            
          END DO

          PRINT *,'Matrix invalid! Trial-DOF not found!'
          STOP
            
        END DO trialspaceloop
      
      END DO ! JDOFE
      
      ! Now we can set up the local matrix in Dentry. Later, we'll plug it into
      ! the global matrix using the positions in Kentry.
      !
      ! The next step is to evaluate the basis functions in the cubature
      ! points on the edge. To compute the jump, this has to be done
      ! for both elements on the edge. 
      !
      ! Figure out which edge on the current element is IMT.
      ! We need this local numbering later to determine on which edge we
      ! have to place the cubature points.
      DO i = 1,IELcount
        
        IEL = p_IelementsAtEdge(i,IMT)
        DO iedge = 1,NVE
          IF (p_IedgesAtElement (iedge,IEL) .EQ. IMT+NVT) EXIT
        END DO
        
        ! Copy the coordinates of the corresponding cubature points
        ! to DcubPtsEval. We calculated the cubature points on the
        ! reference element in advance, so we can simply take them.
        
        p_DcubPtsRef(:,:,i) = DcubPtsRefOnAllEdges (:,:,iedge)
        
      END DO

      ! For the transformation from the reference to the real element, get the 
      ! coordinates of the corners.
      
      DO IEL=1,IELcount
        DO ivt = 1,NVE
          DO i = 1,NDIM2D
            ! ... = DCORVG(i,KVERT(ivt,KMEL(iel,imt)))
            p_Dcoords(i,ivt,IEL) = p_DcornerCoordinates(i, &
                                p_IverticesAtElement(ivt,p_IelementsAtEdge(IEL,IMT)))
          END DO
        END DO
      END DO

      ! Depending on the type of transformation, we must now choose
      ! the mapping between the reference and the real element.
      ! In case we use a nonparametric element or a nonconstant coefficient function,
      ! we need the coordinates of the points on the real element, too.
      IF (bnonparTrial .OR. bnonparTest .OR. (.NOT. rconfig%bconstViscosity)) THEN
      
        CALL trafo_calctrafo_sim (&
              p_elementDistribution%ctrafoType,&
              IELcount,ncubp,p_Dcoords,&
              p_DcubPtsRef,p_Djac(:,:,1:IELcount),p_Ddetj(:,1:IELcount),p_DcubPtsReal)
      
      ELSE
      
        CALL trafo_calctrafo_sim (p_elementDistribution%ctrafoType,&
              IELcount,ncubp,p_Dcoords,&
              p_DcubPtsRef,p_Djac(:,:,1:IELcount),p_Ddetj(:,1:IELcount))
              
      END IF

      ! Cubature prepared. Now we call the element subroutine to evaluate in the
      ! cubature points.
      !
      ! Pass p_DcubPts as point coordinates, which point either to the
      ! coordinates on the reference element (the same for all elements)
      ! or on the real element - depending on whether this is a 
      ! parametric or nonparametric element.
      CALL elem_generic_sim (p_elementDistribution%itestElement, p_Dcoords, &
            p_Djac(:,:,1:IELcount), p_Ddetj(:,1:IELcount), &
            BderTest, DbasTest, ncubp, IELcount, p_DcubPtsTest)
            
      ! Apply the permutation of the local DOF's on the test functions
      ! on element 2. The numbers of the local DOF's on element 1
      ! coincides with the numbers of the local DOF's on the patch.
      ! Those on element 2, we have to renumber according to the permutation
      ! so that they are in the correct order according to the DOF's on the patch.
      !
      ! We copy the values of the basis functions to the space in
      ! p_DcubPtsTest which is reserved for the 3rd element!
      ! That way, DbasTest has:
      !
      ! local DOF on element   :       1       2       3       4
      !
      ! Global DOF on element 1:      10      40      60      20
      ! DbasTest(1..4,*,*,1):    d1psi10 d1psi40 d1psi60 d1psi20
      !
      ! Global DOF on elemenr 2:      60      30      70      50
      ! DbasTest(1..4,*,*,2):    d2psi60 d2psi30 d2psi70 d2psi50
      ! 
      ! and will be transformed to:
      !
      ! local patch-DOF:            1       2       3       4       5      6        7
      ! global DOF:                10      40      60      20      30     70       50
      ! DbasTest(1..7,*,*,1): d1psi10 d1psi40 d1psi60 d1psi20     0.0    0.0      0.0
      ! DbasTest(1..7,*,*,2): --------------------------unused-----------------------
      ! DbasTest(1..7,*,*,3):     0.0     0.0 d2psi60     0.0 d2psi30 d2psi70 d2psi50
      !
      ! ("d1psi10" = grad(psi_10) on element 1, 
      !  "psi_10" = test function on global DOF #10)
      !
      ! That space DbasTest(:,:,:,3) is unused up to now, initialised by 0.0 and 
      ! now overwritten.
      DbasTest (IlocalDofTestRenum(1:indofTestPerElement),:,:,3) &
        = DbasTest (1:indofTestPerElement,:,:,2)
      
      ! Omit the calculation of the trial function values if they
      ! are identical to the test function values.
      IF (.NOT. bidenticalTrialAndTest) THEN
        CALL elem_generic_sim (p_elementDistribution%itrialElement, p_Dcoords, &
            p_Djac(:,:,1:IELcount), p_Ddetj(:,1:IELcount), &
            BderTrial, DbasTrial, ncubp, IELcount, p_DcubPtsTrial)

        ! Apply the renumbering for the 2nd element, store the result in the
        ! space of the 3rd element.          
        DbasTrial (IlocalDofTrialRenum(1:indofTrialPerElement),:,:,3) &
          = DbasTrial (1:indofTrialPerElement,:,:,2)
            
      END IF
      
      ! What do we have now? When looking at the example, we have:
      !
      ! ndofTest = 7
      ! IdofsTest(1..7)             = global DOF's in local numbering 1..7
      ! DbasTest(1..7,*,1..ncubp,1) = values of basis functions on element 1
      !                               in the cubature points, filled by 0.0 in the
      !                               DOF's only appearing at element 2
      ! DbasTest(1..7,*,1..ncubp,3) = values of basis functions on element 2
      !                               in the cubature points, filled by 0.0 in the
      !                               DOF's only appearing at element 1
      !
      ! Now we can start to integrate using this.
      
      ! Get the length of the current edge. It serves as a "determinant"
      ! in the cubature, so we have to divide it by 2 as an edge on the unit inverval
      ! [-1,1] has length 2.
      ivt1 = p_IverticesAtEdge (1,IMT)
      ivt2 = p_IverticesAtEdge (2,IMT)
      dedgelength = &
        SQRT ((p_DcornerCoordinates(1,ivt2)-p_DcornerCoordinates(1,ivt1))**2 &
            + (p_DcornerCoordinates(2,ivt2)-p_DcornerCoordinates(2,ivt1))**2 )
      dedgeweight = 0.5_DP * dedgelength
      
      ! Compute the coefficient in front of the integral:
      ! < Ju,v > = sum_E max(gammastar*nu*h_E, gamma*h_E^2) int_E [grad u] [grad v] ds
      dcoeff = MAX(rconfig%dgammastar * rconfig%dnu * dedgelength, &
                  rconfig%dgamma * dedgelength**2)
      
      ! Now we have the values of the basis functions in all the cubature 
      ! points.
      !
      ! Integrate the jump over the edges. This calculates the local matrix.
      !
      ! Loop through the test basis functions
      DO IDOFE = 1,ndofTest
      
        ! Loop through the trial basis functions
        DO JDOFE = 1,ndofTrial
    
          dval = 0.0_DP
          
          ! Loop through the cubature points to calculate the integral
          ! of the jump. Note that for computing the jump, we have to
          ! look in the inverse order to the cubature points of the neighbour
          ! element!
          ! Remember that the values of the basis functions on the first element
          ! are in p_DbasXXXX (.,.,.,1), while those of the 2nd element are
          ! in p_DbasXXXX (.,.,.,3) (rather than in p_DbasXXXX (.,.,.,2))
          ! by the above construction!
          DO icubp = 1,ncubp
          
            ! [ grad phi ]   ( jump in the derivative of trial basis function)
            ! = [ (d/dx) phi  ,  (d/dy) phi ]
            dphidx = p_DbasTrial (JDOFE,DER_DERIV_X,icubp,1) &
                  - p_DbasTrial (JDOFE,DER_DERIV_X,ncubp-icubp+1,3)

            dphidy = p_DbasTrial (JDOFE,DER_DERIV_Y,icubp,1) &
                  - p_DbasTrial (JDOFE,DER_DERIV_Y,ncubp-icubp+1,3)

            ! [ grad psi ]   ( jump in the derivative of test basis function)
            ! = [ (d/dx) phi  ,  (d/dy) phi ]
            dpsidx = DbasTest (IDOFE,DER_DERIV_X,icubp,1) &
                  - DbasTest (IDOFE,DER_DERIV_X,ncubp-icubp+1,3)

            dpsidy = DbasTest (IDOFE,DER_DERIV_Y,icubp,1) &
                  - DbasTest (IDOFE,DER_DERIV_Y,ncubp-icubp+1,3)
            
          
            ! Compute int_edge ( [grad phi_i] [grad phi_j] )
            dval = dval + Domega(icubp) * dedgeweight * &
                          (dphidx*dpsidx + dphidy*dpsidy)
          
          END DO

          ! Add the contribution to the local matrix -- weighted by the
          ! Omega from the cubature formula and the length of the edge.
          Dentry ((IDOFE-1)*ndofTrial+JDOFE) = &
            Dentry ((IDOFE-1)*ndofTrial+JDOFE) + dcoeff*dval

        END DO
      
      END DO
      
      ! Add the local matrix to the global matrix
      DO IDOFE = 0,ndofTest-1
        DO JDOFE = 1,ndofTrial
        
          p_Da(Kentry(IDOFE*ndofTrial+JDOFE)) = &
            p_Da(Kentry(IDOFE*ndofTrial+JDOFE)) + &
            rconfig%dtheta*Dentry (IDOFE*ndofTrial+JDOFE)
        
        END DO
      END DO
      
      ! Proceed with next edge
    
    END DO ! IMT

    ! Clean up allocated arrays and memory.
    DEALLOCATE(DcubPtsRefOnAllEdges)
    CALL domint_doneIntegration(rintSubset)
    
    DEALLOCATE(DbasTest)
    DEALLOCATE(DbasTrial)

    DEALLOCATE(Kentry) 
    DEALLOCATE(Dentry)

    DEALLOCATE(IdofsTestTempl) 
    DEALLOCATE(IdofsTrialTempl)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE conv_JumpStabilisation2d ( &
                           rvecPrimary, rvecSecondary, dprimWeight, dsecWeight,&
                           rconfig, cdef, &
                           rmatrix, rsolution, rdefect, &
                           IvelocityComp)

!<description>
  ! Edge oriented stabilisation technique. Wrapper routine.
  ! 2D-version (X- and Y-velocity).
  !
  ! rvecPrimary, rvecSecondary are two velocity field vectors for the X-
  ! and Y-veclocity; IvelocityComp defines which components of these
  ! vectors contains the X- and which contains the Y-velocity.
  ! The final velocity vector field is then computed as a weighted average
  ! of these two:
  !
  !  $$ u_1  =  dprimWeight * rvecPrimary  +  dsecWeight * rvecSecondary $$
  !
  ! $u_2 = rsolution(.)$ defines the second velocity field inside of
  ! the grad-term.
  !
  ! The switch cdef decides on whether the routine sets up the nonlinear
  ! defect, the nonlinear matrix or both.
  !
  ! The configuration how the routine should react is to be configured
  ! in the configuration block rconfig.
  !
  ! Remark: Actually, rvecPrimary/rvecSecondary are not used, but might be
  !  needed in future versions...
!</description>

!<input>

  ! Primary velocity field for the computation of $u_1$
  TYPE(t_vectorBlock), INTENT(IN), TARGET :: rvecPrimary
  
  ! Secondary velocity field for the computation of $u_1$
  TYPE(t_vectorBlock), INTENT(IN), TARGET :: rvecSecondary
  
  ! Weighting factor for rvecPrimary.
  REAL(DP), INTENT(IN) :: dprimWeight
  
  ! Weighting factor for rvecSecondary.
  REAL(DP), INTENT(IN) :: dsecWeight
  
  ! Configuration block for the streamline diffusion scheme
  TYPE(t_jumpStabilisation), INTENT(IN) :: rconfig
  
  ! Computation/defect correction method. One of the CONV_MODxxxx constants:
  ! CONV_MODMATRIX: Set up the nonlinear matrix. rmatrix must be present, the 
  !                 nonlinear part is added to the matrix.
  ! CONV_MODDEFECT: Set up the nonlinear defect. rdefect and rsolution must be 
  !                 present.
  ! CONV_MODBOTH  : Set up the nonlinear matrix as well as the nonlinear defect.
  !                 rmatrix, rdefect and rsolution must all be present.
  INTEGER, INTENT(IN) :: cdef

  ! OPTIONAL: Solution vector u_2.
  ! Must be present if cdef=CONV_MODDEFECT or =CONV_MODBOTH.
  TYPE(t_vectorBlock), INTENT(IN), TARGET, OPTIONAL :: rsolution
  
  ! OPTIONAL: 
  ! Index block that specifies which component in rvecPrimary / rvecSecondary /
  ! rsolution / rdefect is the X-velocity and which one is the Y-velocity.
  !  IvelocityComp(1) gives the number of the X-velocity (usually = 1),
  !  IvelocityComp(2) gives the number of the Y-velocity (usually = 2).
  ! If not present, IvelocityComp=(/1,2/) is assumed, thus the X-velocity
  ! must be in rsolution\%RvectorBlock(1) and the Y-velocity in
  ! rsolution\%RvectorBlock(2).
  INTEGER, DIMENSION(2), INTENT(IN), OPTIONAL :: IvelocityComp
!</input>

!<inputoutput>
  ! System matrix.
  ! The content of the matrix must be present if cdef=CONV_MODMATRIX or 
  ! =CONV_MODBOTH, otherwise only the structure is used.
  ! The nonlinear operator is added to the matrix.
  TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrix
  
  ! OPTIONAL: Defect vector.
  ! Must have the same structure as rsolution/rvecPrimary/rvecSecondary.
  ! Must be present if cdef=CONV_MODDEFECT or =CONV_MODBOTH.
  ! The nonlinear part is subtracted from this vector: 
  ! $r = r - \theta * u_1*grad(u_2)$
  TYPE(t_vectorBlock), INTENT(INOUT), OPTIONAL, TARGET :: rdefect
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER, DIMENSION(2) :: Icomp
    TYPE(t_vectorScalar), POINTER :: p_rsolX,p_rsolY,p_rdefectX,p_rdefectY
    REAL(DP), DIMENSION(:), POINTER :: p_DsolX,p_DsolY,p_DdefectX,p_DdefectY
    
    ! At first check the input parameters that everything is present what
    ! we need:
    IF ((cdef .EQ. CONV_MODDEFECT) .OR. (cdef .EQ. CONV_MODBOTH)) THEN
      IF ((.NOT. PRESENT(rsolution)) .OR. (.NOT. PRESENT(rdefect))) THEN
        PRINT *,'EOS: Solution/defect vector not present!'
        STOP
      END IF
    END IF
    
    ! Get the actual subvectors from the velocity vectors that define
    ! the X- and Y-velocity.
    IF (PRESENT(IvelocityComp)) THEN
      Icomp = IvelocityComp
    ELSE
      Icomp = (/1,2/)
    END IF
    
    IF (PRESENT(rsolution)) THEN
      p_rsolX => rsolution%RvectorBlock(Icomp(1))
      p_rsolY => rsolution%RvectorBlock(Icomp(2))
    ELSE
      NULLIFY(p_rsolX)
      NULLIFY(p_rsolY)
    END IF
    
    IF (PRESENT(rsolution)) THEN
      p_rdefectX => rdefect%RvectorBlock(Icomp(1))
      p_rdefectY => rdefect%RvectorBlock(Icomp(2))
    ELSE
      NULLIFY(p_rdefectX)
      NULLIFY(p_rdefectY)
    END IF
      
    ! At the moment, we only support a rather limited set of configurations:
    ! Matrix and vectors must all be double precision, matrix must be format 
    ! 7 or 9, discretisation must be Q1~, constant viscosity.
    IF ((rmatrix%cmatrixFormat .NE. LSYSSC_MATRIX9) .AND. &
        (rmatrix%cmatrixFormat .NE. LSYSSC_MATRIX7)) THEN
      PRINT *,'EOS: Unsupported matrix format'
      STOP
    END IF

    IF (rmatrix%p_rspatialDiscretisation%ccomplexity .NE. SPDISC_UNIFORM) THEN
      PRINT *,'EOS: Unsupported discretisation.'
      STOP
    END IF

    !IF ((rvecPrimary%cdataType .NE. ST_DOUBLE) .OR. &
    !    (rvecSecondary%cdataType .NE. ST_DOUBLE)) THEN
    !  PRINT *,'EOS: Unsupported vector data type in velocity.'
    !  STOP
    !END IF
    ! 
    !IF (PRESENT(rdefect)) THEN
    !  IF ((rsolution%cdataType .NE. ST_DOUBLE) .OR. &
    !      (rdefect%cdataType .NE. ST_DOUBLE)) THEN
    !    PRINT *,'EOS: Unsupported vector data type in solution/defect'
    !    STOP
    !  END IF
    !END IF
    
    IF (.NOT. rconfig%bconstViscosity) THEN
      PRINT *,'EOS: Only constant viscosity supported at the moment!'
      STOP
    END IF
    
    IF (rconfig%dnu .EQ. SYS_INFINITY) THEN
      PRINT *,'EOS: Viscosity parameter nu not initialised!'
      STOP
    END IF
    
    IF (rconfig%cjump .EQ. CONV_JUMP_UNIFIEDEDGE) THEN
      IF (PRESENT(rdefect)) THEN
        CALL lsyssc_getbase_double (p_rsolX   ,p_DsolX   )
        CALL lsyssc_getbase_double (p_rsolY   ,p_DsolY   )
        CALL lsyssc_getbase_double (p_rdefectX,p_DdefectX)
        CALL lsyssc_getbase_double (p_rdefectY,p_DdefectY)
        
        CALL conv_ueoJumpStabil2d_mv_unidble ( &
                    rmatrix,cdef,rconfig, &
                    p_DsolX,p_DsolY,p_DdefectX,p_DdefectY)
                      
      ELSE
      
        CALL conv_ueoJumpStabil_double_uni (rmatrix,rconfig)

      END IF
    ELSE
      PRINT *,'EOS: Unknown jump stabilisation!'
      STOP
    END IF

  END SUBROUTINE
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE conv_ueoJumpStabil2d_mv_unidble ( &
                  rmatrixScalar,cdef,rconfig, &
                  Du1,Du2,Ddef1,Ddef2)
!<description>
  ! Unified edge oriented jump stabilisation.
  !
  ! Adds the unified edge oriented jump stabilisation to the matrix rmatrix:
  ! $$< Ju,v > = \sum_E \max(\gamma^{*}*\nu*h_E, \gamma h_E^2) 
  !            \int_E [grad u] [grad v] ds$$
  ! or updates a defect vector.
  !
  ! Uniform discretisation, double precision structure-7 and 9 matrix.
  !
  ! For a rerefence about the stabilisation, see
  ! [Ouazzi, A.; Finite Element Simulation of Nonlinear Fluids, Application
  ! to Granular Material and Powder; Shaker Verlag, ISBN 3-8322-5201-0, p. 55ff]
  !
  ! Remarks:\\
  !  
  !  1.) In the nonlinear iteration, as a right hand side there arises
  !   a defect vector D, which linear part can easily being assembled.
  !   However, there is a nonlinearity to be included into that vector,
  !   too. By setting cdef=1,2, this routine incorporates the nonlinearity
  !   into that vector, using the formula
  !
  !            $$ D = D - dtheta * UUx * grad (Ux) $$
  !   
  ! 2.) WARNING: For edge oriented stabilisation, the underlying matrix rmatrix
  !   must have an extended stencil! The matrix structure must be set up with
  !   the BILF_MATC_EDGEBASED switch!!!
  !
!</description>
  
!<input>

  ! Computation/defect correction method. One of the CONV_MODxxxx constants:
  ! CONV_MODMATRIX: Set up the nonlinear matrix. rmatrix must be present, the 
  !                 nonlinear part is added to the matrix.
  ! CONV_MODDEFECT: Set up the nonlinear defect. rdefect and rsolution must be 
  !                 present.
  ! CONV_MODBOTH  : Set up the nonlinear matrix as well as the nonlinear defect.
  !                 rmatrix, rdefect and rsolution must all be present.
  INTEGER, INTENT(IN) :: cdef
  
  ! Configuration block for the jump stabilisation
  TYPE(t_jumpStabilisation), INTENT(IN) :: rconfig

  ! OPTIONAL: X-velocity of $u_2$. Must be present if cdef=CONV_MODDEFECT
  ! or cdef=CONV_MODBOTH.
  REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: Du1
  
  ! Y-velocity of $u_2$. Must be present if cdef=CONV_MODDEFECT
  ! or cdef=CONV_MODBOTH.
  REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: Du2
  
!</input>

!<inputoutput>
  ! The system matrix. Must be format 7 or 9.
  TYPE(t_matrixScalar), INTENT(INOUT), TARGET :: rmatrixScalar
  
  ! OPTIONAL: X-defect vector. Must be present if cdef=CONV_MODDEFECT
  ! or =CONV_MODBOTH.
  REAL(DP), DIMENSION(:), INTENT(INOUT), OPTIONAL :: Ddef1(:)
  
  ! OPTIONAL: Y-defect vector. Must be present if cdef=CONV_MODDEFECT
  ! or =CONV_MODBOTH.
  REAL(DP), DIMENSION(:), INTENT(INOUT), OPTIONAL :: Ddef2(:)
!</inputoutput>

!</subroutine>


  ! local variables
  INTEGER(PREC_DOFIDX) :: irow, jcol, idof
  INTEGER(PREC_EDGEIDX) :: IMT
  INTEGER(PREC_POINTIDX) :: ivt,ivt1,ivt2,NVT
  INTEGER(PREC_ELEMENTIDX) :: IEL
  LOGICAL :: bIdenticalTrialAndTest
  INTEGER :: IELcount,IDOFE, JDOFE, i, NVE, iedge
  REAL(DP) :: dval,dedgelength,dedgeweight,dphidx,dphidy,dpsidx,dpsidy,dcoeff
  
  ! Pointer to KLD, KCOL, DA
  INTEGER(PREC_MATIDX), DIMENSION(:), POINTER :: p_Kld
  INTEGER(PREC_VECIDX), DIMENSION(:), POINTER :: p_Kcol
  REAL(DP), DIMENSION(:), POINTER :: p_Da
  
  ! An allocateable array accepting the DOF's of a set of elements.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsTestTempl
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsTrialTempl
  INTEGER(PREC_DOFIDX), DIMENSION(EL_MAXNBAS*2), TARGET :: IdofsTest
  INTEGER(PREC_DOFIDX), DIMENSION(EL_MAXNBAS*2), TARGET :: IdofsTrial
  INTEGER(PREC_DOFIDX), DIMENSION(:), POINTER :: p_IdofsTrial
  
  ! Arrays saving the local DOF numbers belonging to the global
  ! DOF numbers in IdofsTest/IdofsTrial  
  INTEGER, DIMENSION(EL_MAXNBAS*2),TARGET :: IlocalDofsTest,IlocalDofsTrial
  INTEGER(PREC_DOFIDX), DIMENSION(:), POINTER :: p_IlocalDofsTrial
  
  ! Renumbering strategy for local DOF's
  INTEGER, DIMENSION(EL_MAXNBAS), TARGET :: IlocalDofTestRenum,IlocalDofTrialRenum
  INTEGER, DIMENSION(:), POINTER :: p_IlocalDofTrialRenum
  
  ! Number of local DOF's on the patch
  INTEGER :: ndofTest,ndofTrial
  
  ! Number of local degees of freedom for trial and test functions
  INTEGER :: indofTrialPerElement, indofTestPerElement
  
  ! The triangulation structure - to shorten some things...
  TYPE(t_triangulation), POINTER :: p_rtriangulation
  
  ! Some triangulation arrays we need frequently
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IneighboursAtElement
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IelementsAtEdge
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement
  INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
  INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtEdge
  REAL(DP), DIMENSION(:,:), POINTER :: p_DcornerCoordinates
  
  ! Current element distribution
  TYPE(t_elementDistribution), POINTER :: p_elementDistribution
  
  ! Underlying discretisation structure
  TYPE(t_spatialDiscretisation), POINTER :: p_rdiscretisation

  ! Arrays for saving Jacobian determinants and matrices
  REAL(DP), DIMENSION(:,:), POINTER :: p_Ddetj
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_Djac

  ! Allocateable arrays for the values of the basis functions - 
  ! for test and trial spaces.
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: DbasTest,DbasTrial
  REAL(DP), DIMENSION(:,:,:,:), POINTER :: p_DbasTrial

  ! Local matrices, used during the assembly.
  ! Values and positions of values in the global matrix.
  INTEGER(PREC_DOFIDX), DIMENSION(:), ALLOCATABLE :: Kentry
  REAL(DP), DIMENSION(:), ALLOCATABLE :: Dentry

  ! A t_domainIntSubset structure that is used for storing information
  ! and passing it to callback routines.
  TYPE(t_domainIntSubset) :: rintSubset
  
  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: DcubPtsRefOnAllEdges

  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPtsRef

  ! An array receiving the coordinates of cubature points on
  ! the real element for all elements in a set.
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPtsReal
  
  ! Pointer to the cubature points for trial and test space.
  ! These point either to p_DcubPtsRef or p_DcubPtsReal, depending on
  ! whether the element is parametric or not.
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPtsTrial
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPtsTest
  
  ! Array with coordinates of the corners that form the real element.
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_Dcoords
  
  ! Cubature point weights
  REAL(DP), DIMENSION(CUB_MAXCUBP) :: Domega
  
  ! Cubature point coordinates on the reference element
  REAL(DP), DIMENSION(CUB_MAXCUBP, NDIM3D) :: Dxi1D,Dxi2D
  
  ! number of cubature points on the reference element
  INTEGER :: ncubp,icubp
  
  ! Test/trial functions nonparameric?
  LOGICAL :: bnonparTrial,bnonparTest
  
  ! Derivative specifiers
  LOGICAL, DIMENSION(EL_MAXNDER) :: BderTrial, BderTest

    ! Get a pointer to the triangulation and discretisation.
    p_rtriangulation => rmatrixScalar%p_rspatialDiscretisation%p_rtriangulation
    p_rdiscretisation => rmatrixScalar%p_rspatialDiscretisation
    
    ! Get Kvert, Kadj, Kmid, Kmel,...
    CALL storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
                                p_IverticesAtElement)
    CALL storage_getbase_int2d (p_rtriangulation%h_IneighboursAtElement,&
                                p_IneighboursAtElement)
    CALL storage_getbase_int2d (p_rtriangulation%h_IedgesAtElement,&
                                p_IedgesAtElement)
    CALL storage_getbase_int2d (p_rtriangulation%h_IelementsAtEdge,&
                                p_IelementsAtEdge)
    CALL storage_getbase_int2d (p_rtriangulation%h_IverticesAtEdge,&
                                p_IverticesAtEdge)
    CALL storage_getbase_double2d (p_rtriangulation%h_DcornerCoordinates,&
                                  p_DcornerCoordinates)
                                
    NVT = p_rtriangulation%NVT
    
    ! Get KLD, KCol...
    CALL storage_getbase_int (rmatrixScalar%h_KLD,p_KLD)
    CALL storage_getbase_int (rmatrixScalar%h_Kcol,p_Kcol)
    CALL storage_getbase_double (rmatrixScalar%h_Da,p_Da)
    
    ! Activate the one and only element distribution
    p_elementDistribution => p_rdiscretisation%RelementDistribution(1)

    ! Get the number of local DOF's for trial and test functions
    indofTrialPerElement = elem_igetNDofLoc(p_elementDistribution%itrialElement)
    indofTestPerElement = elem_igetNDofLoc(p_elementDistribution%itestElement)
    
    ! Triangle elements? Quad elements?
    NVE = elem_igetNVE(p_elementDistribution%itrialElement)
    
    ! Get the number of corner vertices of the element
    IF (NVE .NE. elem_igetNVE(p_elementDistribution%itestElement)) THEN
      PRINT *,'conv_ueoJumpStabil2d_double_uni: element spaces incompatible!'
      STOP
    END IF
    
    ! Initialise the cubature formula,
    ! Get cubature weights and point coordinates on the reference line [-1,1]
    CALL cub_getCubPoints(rconfig%ccubType, ncubp, Dxi1D, Domega)

    ! Allocate arrays accepting cubature point coordinates.
    ! We have ncubp cubature points on every egde.
    ! DcubPtsRef saves all possible combinations, which edge of
    ! one element might interact with one edge of one neighbour
    ! element.
    ALLOCATE(DcubPtsRefOnAllEdges(NDIM2D,ncubp,NVE))
    
    ! Put the 1D cubature points from to all of the edges on the
    ! reference element, so we can access them quickly.
    DO iedge = 1,NVE
      CALL trafo_mapCubPts1Dto2DRefQuad(iedge, ncubp, Dxi1D, Dxi2D)
      DO i=1,ncubp
        DcubPtsRefOnAllEdges(1,i,iedge) = Dxi2D(i,1)
        DcubPtsRefOnAllEdges(2,i,iedge) = Dxi2D(i,2)
      END DO
    END DO
    
    ! Get from the trial element space the type of coordinate system
    ! that is used there:
    i = elem_igetCoordSystem(p_elementDistribution%itrialElement)
    
    ! Allocate memory and get local references to it.
    CALL domint_initIntegration (rintSubset,2,ncubp,i,NDIM2D)
    p_DcubPtsRef =>  rintSubset%p_DcubPtsRef
    p_DcubPtsReal => rintSubset%p_DcubPtsReal
    p_Djac =>        rintSubset%p_Djac
    p_Ddetj =>       rintSubset%p_Ddetj
    p_Dcoords =>     rintSubset%p_DCoords
    
    ! Check if one of the trial/test elements is nonparametric
    bnonparTrial = elem_isNonparametric(p_elementDistribution%itrialElement)
    bnonparTest  = elem_isNonparametric(p_elementDistribution%itestElement)
                    
    ! Let p_DcubPtsTrial / p_DcubPtsTest point either to DcubPtsRealEval or
    ! DcubPtsRefEval - depending on whether the space is parametric or not.
    IF (bnonparTrial) THEN
      p_DcubPtsTrial => p_DcubPtsReal
    ELSE
      p_DcubPtsTrial => p_DcubPtsRef
    END IF
    
    IF (bnonparTest) THEN
      p_DcubPtsTest => p_DcubPtsReal
    ELSE
      p_DcubPtsTest => p_DcubPtsRef
    END IF
    
    ! Allocate arrays saving the local matrices for all elements
    ! in an element set. We allocate the arrays large enough...
    ALLOCATE(Kentry(indofTestPerElement*2*indofTrialPerElement*2))
    ALLOCATE(Dentry(indofTestPerElement*2*indofTrialPerElement*2))
    
    ! Allocate memory for obtaining DOF's:
    ALLOCATE(IdofsTrialTempl(indofTrialPerElement,2))
    ALLOCATE(IdofsTestTempl(indofTestPerElement,2))
    
    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  ALLOCATE(DbasTest(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock+1))
    !  ALLOCATE(DbasTrial(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock+1))
    ! would lead to nonused memory blocks in these arrays during the assembly, 
    ! which reduces the speed by 50%!
    !
    ! We allocate space for 3 instead of 2 elements. The reason is that
    ! we later permute the values of the basis functions to get
    ! a local numbering on the patch. That's also the reason, we allocate
    ! not indofTestPerElement elements, but even indofTestPerElement,
    ! which is more than enough space to hold the values of the DOF's of
    ! a whole element patch.
    
    ALLOCATE(DbasTest(indofTestPerElement*2, &
            elem_getMaxDerivative(p_elementDistribution%itestElement),&
            ncubp,3))
    ALLOCATE(DbasTrial(indofTrialPerElement*2,&
            elem_getMaxDerivative(p_elementDistribution%itrialElement), &
            ncubp,3))
             
    ! Test if trial/test functions are identical.
    ! We don't rely on bidenticalTrialAndTest purely, as this does not
    ! indicate whether there are identical trial and test functions
    ! in one block!
    bIdenticalTrialAndTest = &
      p_elementDistribution%itrialElement .EQ. p_elementDistribution%itestElement
      
    ! Let p_IdofsTrial point either to IdofsTrial or to the DOF's of the test
    ! space IdofTest (if both spaces are identical). 
    ! We create a pointer for the trial space and not for the test space to
    ! prevent pointer-arithmetic in the innerst loop below!
    ! The same for p_IlocalDofsTrial and the permutation array p_IlocalDofTrialRenum.
    ! Let p_DbasTrial point either to DbasTrial or DbasTest, depending on
    ! whether the spaces are identical.
    IF (bIdenticalTrialAndTest) THEN
      p_IdofsTrial => IdofsTest
      p_IlocalDofsTrial => IlocalDofsTest
      p_IlocalDofTrialRenum => IlocalDofTestRenum
      p_DbasTrial => DbasTest
    ELSE
      p_IdofsTrial => IdofsTrial
      p_IlocalDofsTrial => IlocalDofsTrial
      p_IlocalDofTrialRenum => IlocalDofTrialRenum
      p_DbasTrial => DbasTrial
    END IF
    
    ! Set up which derivatives to compute in the basis functions: X/Y-derivative
    BderTrial = .FALSE.
    BderTrial(DER_DERIV_X) = .TRUE.
    BderTrial(DER_DERIV_Y) = .TRUE.

    BderTest = .FALSE.
    BderTest(DER_DERIV_X) = .TRUE.
    BderTest(DER_DERIV_Y) = .TRUE.

    ! We have always 2 elements in each patch...
    IELcount = 2
      
    ! We loop through all edges
    DO IMT = 1,p_rtriangulation%NMT
    
      ! Check if we have 1 or 2 elements on the current edge
      IF (p_IelementsAtEdge (2,IMT) .EQ. 0) THEN
      
        ! This only happens if we have a boundary edge.
        ! The boundary handling is... doing nothing! So we can skip
        ! the handling of that edge completely!
        CYCLE
      
      END IF
      
      ! Fill the basis function arrays with 0. Essential, as only parts
      ! of them are overwritten later.
      DbasTest = 0.0_DP
      DbasTrial = 0.0_DP
      
      ! On an example, we now show the relationship between all the DOF's
      ! we have to consider in the current situation. We have two elements,
      ! let's say IEL1 and IEL2, with their local and global DOF's
      ! (example for Q1~):
      !
      !   local DOF's on the           corresponding global
      !   element and on the patch     DOF's
      !
      !    +----4----+----7----+       +----20---+----50---+
      !    |    4    |    4    |       |         |         |
      !    |         |         |       |         |         |
      !    1 1     3 3 1     3 6       10        60        70
      !    |         |         |       |         |         |
      !    |    2    |    2    |       |         |         |
      !    +----2----+----3----+       +----40---+----30---+
      !        IEL1      IEL2              IEL1      IEL2
      !
      !
      ! On every element, we have 4 local DOF's (1..4). On the other hand,
      ! we have "local DOF's on the patch" (1..7), ehich we call "patch DOF's"
      ! from now on. To every local DOF, there belongs a global DOF
      ! (10,20,30,... or whatever), which gives the coefficient of the basis
      ! function.
      !
      ! Numbering that way, the local DOF's of IEL1 obviously coincide 
      ! with the first couple of local DOF's of the element patch! Only
      ! the local DOF's of element IEL2 make trouble, as they have another
      ! numbering.
      
      ! Get the global DOF's of the 1 or two elements
      CALL dof_locGlobMapping_mult(p_rdiscretisation, &
                                  p_IelementsAtEdge (1:IELcount,IMT), &
                                  .TRUE., IdofsTestTempl)
                                   
      ! Some of the DOF's on element 2 may coincide with DOF's on element 1.
      ! More precisely, some must coincide! Therefore, we now have to collect the
      ! DOF's uniquely and to figure out, which local DOF's of element 2
      ! must renumbered to the appropriate local patch-DOF (like in the
      ! example, where local DOF 1 of element IEL2 must be renumbered
      ! to local patch DOF 3!
      !
      ! As the first couple of local DOF's of IEL1 coincide with the local
      ! DOF's of the patch, we can simply copy them:
      
      ndofTest = indofTestPerElement
      IdofsTest(1:ndofTest) = IdofsTestTempl(1:ndofTest,1)
      
      ! Furthermore, we read IdofsTestTempl and store the DOF's there in IdofsTest,
      ! skipping all DOF's we already have and setting up the renumbering strategy
      ! of local DOF's on element IEL2 to patch DOF's.
      
      skiploop: DO IDOFE = 1,indofTestPerElement
        
        ! Do we have the DOF? 
        idof = IdofsTestTempl(IDOFE,IELcount)
        DO JDOFE=1,ndofTest
          IF (IdofsTest(JDOFE) .EQ. idof) THEN
            ! Yes, we have it.
            ! That means, the local DOF idof has to be mapped to the
            ! local patch dof...
            IlocalDofTestRenum (IDOFE) = JDOFE
            
            ! Proceed with the next one.
            CYCLE skiploop
          END IF
        END DO
        
        ! We don't have that DOF! Append it to IdofsTest.
        ndofTest = ndofTest+1
        IdofsTest(ndofTest) = idof
        IlocalDofTestRenum (IDOFE) = ndofTest
        
        ! Save also the number of the local DOF.
        ! Note that the global DOF's in IdofsTestTempl(1..indofTestPerElement)
        ! belong to the local DOF's 1..indofTestPerElement -- in that order!
        IlocalDofsTest(ndofTest) = IDOFE
        
      END DO skiploop
        
      ! If the trial functions are different, get those DOF's
      IF (.NOT. bIdenticalTrialAndTest) THEN
        CALL dof_locGlobMapping_mult(p_rdiscretisation, &
                                    p_IelementsAtEdge (1:IELcount,IMT), &
                                    .FALSE., IdofsTrialTempl)

        ! Collect the DOF's uniquely in IdofsTrial.
        ! The DOF's of the first element are taken as they are.
        
        ndofTrial = indofTrialPerElement
        IdofsTrial(1:ndofTrial) = IdofsTrialTempl(1:ndofTrial,1)
        
        ! From the 2nd element on, we have to check which
        ! DOF's already appear in the first element.
        
        skiploop2: DO IDOFE = 1,indofTrialPerElement
          
          ! Do we have the DOF? 
          idof = IdofsTrialTempl(IDOFE,IELcount)
          DO JDOFE=1,ndofTrial
            IF (IdofsTrial(JDOFE) .EQ. idof) THEN
              ! Yes, we have it.
              ! That means, the local DOF idof has to be mapped to the
              ! local patch dof...
              IlocalDofTrialRenum (IDOFE) = JDOFE
              
              ! Proceed with the next one.
              CYCLE skiploop2
            END IF
          END DO
          
          ! We don't have that DOF! Append it to IdofsTrial.
          ndofTrial = ndofTrial+1
          IdofsTrial(ndofTrial) = idof
          IlocalDofTrialRenum (IDOFE) = ndofTrial
          
          ! Save also the number of the local DOF.
          ! Note that the global DOF's in IdofsTrialTempl(1..indofTrialPerElement)
          ! belong to the local DOF's 1..indofTrialPerElement -- in that order!
          IlocalDofsTrial(ndofTrial) = IDOFE
          
        END DO skiploop2
          
      ELSE
       
        ! Number of DOF's in test and trial space the same, as DOF's are the same.
        ndofTrial = ndofTest
        
      END IF
      
      ! Now we know: Our 'local' matrix (consisting of only these DOF's we just
      ! calculated) is a ndofsTest*ndofsTrial matrix.
      !
      ! Now extract the corresponding entries from the matrix.
      ! Kentry is an index where the entries of the local matrix Dentry
      ! (which we build up later) can be found in the global matrix.
      !
      ! The test space gives us the rows; loop through them
      
      DO IDOFE = 0,ndofTest-1
      
        ! The 'local' DOF IDOFE corresponds in the global matrix to line...
        irow = IdofsTest (1+IDOFE)
        
        ! Loop through that line to find the columns, indexed by the local
        ! DOF's in the trial space.
        trialspaceloop: DO JDOFE = 1,ndofTrial
          
          DO jcol = p_Kld(irow),p_Kld(irow+1)-1
            
            IF (p_Kcol(jcol) .EQ. p_IdofsTrial(JDOFE)) THEN
            
              ! Found! Put as destination pointer to Kentry
              Kentry (IDOFE*ndofTrial+JDOFE) = jcol
              
              ! Clear local matrix
              Dentry (IDOFE*ndofTrial+JDOFE) = 0.0_DP
              
              ! Jump out of the loop, proceed with next column
              CYCLE trialspaceloop
            
            END IF
            
          END DO

          PRINT *,'Matrix invalid! Trial-DOF not found!'
          STOP
            
        END DO trialspaceloop
      
      END DO ! JDOFE
      
      ! Now we can set up the local matrix in Dentry. Later, we'll plug it into
      ! the global matrix using the positions in Kentry.
      !
      ! The next step is to evaluate the basis functions in the cubature
      ! points on the edge. To compute the jump, this has to be done
      ! for both elements on the edge. 
      !
      ! Figure out which edge on the current element is IMT.
      ! We need this local numbering later to determine on which edge we
      ! have to place the cubature points.
      DO i = 1,IELcount
        
        IEL = p_IelementsAtEdge(i,IMT)
        DO iedge = 1,NVE
          IF (p_IedgesAtElement (iedge,IEL) .EQ. IMT+NVT) EXIT
        END DO
        
        ! Copy the coordinates of the corresponding cubature points
        ! to DcubPtsEval. We calculated the cubature points on the
        ! reference element in advance, so we can simply take them.
        
        p_DcubPtsRef(:,:,i) = DcubPtsRefOnAllEdges (:,:,iedge)
        
      END DO

      ! For the transformation from the reference to the real element, get the 
      ! coordinates of the corners.
      
      DO IEL=1,IELcount
        DO ivt = 1,NVE
          DO i = 1,NDIM2D
            ! ... = DCORVG(i,KVERT(ivt,KMEL(iel,imt)))
            p_Dcoords(i,ivt,IEL) = p_DcornerCoordinates(i, &
                                p_IverticesAtElement(ivt,p_IelementsAtEdge(IEL,IMT)))
          END DO
        END DO
      END DO

      ! Depending on the type of transformation, we must now choose
      ! the mapping between the reference and the real element.
      ! In case we use a nonparametric element or a nonconstant coefficient function,
      ! we need the coordinates of the points on the real element, too.
      IF (bnonparTrial .OR. bnonparTest .OR. (.NOT. rconfig%bconstViscosity)) THEN
      
        CALL trafo_calctrafo_sim (&
              p_elementDistribution%ctrafoType,&
              IELcount,ncubp,p_Dcoords,&
              p_DcubPtsRef,p_Djac(:,:,1:IELcount),p_Ddetj(:,1:IELcount),p_DcubPtsReal)
      
      ELSE
      
        CALL trafo_calctrafo_sim (p_elementDistribution%ctrafoType,&
              IELcount,ncubp,p_Dcoords,&
              p_DcubPtsRef,p_Djac(:,:,1:IELcount),p_Ddetj(:,1:IELcount))
              
      END IF

      ! Cubature prepared. Now we call the element subroutine to evaluate in the
      ! cubature points.
      !
      ! Pass p_DcubPts as point coordinates, which point either to the
      ! coordinates on the reference element (the same for all elements)
      ! or on the real element - depending on whether this is a 
      ! parametric or nonparametric element.
      CALL elem_generic_sim (p_elementDistribution%itestElement, p_Dcoords, &
            p_Djac(:,:,1:IELcount), p_Ddetj(:,1:IELcount), &
            BderTest, DbasTest, ncubp, IELcount, p_DcubPtsTest)
            
      ! Apply the permutation of the local DOF's on the test functions
      ! on element 2. The numbers of the local DOF's on element 1
      ! coincides with the numbers of the local DOF's on the patch.
      ! Those on element 2, we have to renumber according to the permutation
      ! so that they are in the correct order according to the DOF's on the patch.
      !
      ! We copy the values of the basis functions to the space in
      ! p_DcubPtsTest which is reserved for the 3rd element!
      ! That way, DbasTest has:
      !
      ! local DOF on element   :       1       2       3       4
      !
      ! Global DOF on element 1:      10      40      60      20
      ! DbasTest(1..4,*,*,1):    d1psi10 d1psi40 d1psi60 d1psi20
      !
      ! Global DOF on elemenr 2:      60      30      70      50
      ! DbasTest(1..4,*,*,2):    d2psi60 d2psi30 d2psi70 d2psi50
      ! 
      ! and will be transformed to:
      !
      ! local patch-DOF:            1       2       3       4       5      6        7
      ! global DOF:                10      40      60      20      30     70       50
      ! DbasTest(1..7,*,*,1): d1psi10 d1psi40 d1psi60 d1psi20     0.0    0.0      0.0
      ! DbasTest(1..7,*,*,2): --------------------------unused-----------------------
      ! DbasTest(1..7,*,*,3):     0.0     0.0 d2psi60     0.0 d2psi30 d2psi70 d2psi50
      !
      ! ("d1psi10" = grad(psi_10) on element 1, 
      !  "psi_10" = test function on global DOF #10)
      !
      ! That space DbasTest(:,:,:,3) is unused up to now, initialised by 0.0 and 
      ! now overwritten.
      DbasTest (IlocalDofTestRenum(1:indofTestPerElement),:,:,3) &
        = DbasTest (1:indofTestPerElement,:,:,2)
      
      ! Omit the calculation of the trial function values if they
      ! are identical to the test function values.
      IF (.NOT. bidenticalTrialAndTest) THEN
        CALL elem_generic_sim (p_elementDistribution%itrialElement, p_Dcoords, &
            p_Djac(:,:,1:IELcount), p_Ddetj(:,1:IELcount), &
            BderTrial, DbasTrial, ncubp, IELcount, p_DcubPtsTrial)

        ! Apply the renumbering for the 2nd element, store the result in the
        ! space of the 3rd element.          
        DbasTrial (IlocalDofTrialRenum(1:indofTrialPerElement),:,:,3) &
          = DbasTrial (1:indofTrialPerElement,:,:,2)
            
      END IF
      
      ! What do we have now? When looking at the example, we have:
      !
      ! ndofTest = 7
      ! IdofsTest(1..7)             = global DOF's in local numbering 1..7
      ! DbasTest(1..7,*,1..ncubp,1) = values of basis functions on element 1
      !                               in the cubature points, filled by 0.0 in the
      !                               DOF's only appearing at element 2
      ! DbasTest(1..7,*,1..ncubp,3) = values of basis functions on element 2
      !                               in the cubature points, filled by 0.0 in the
      !                               DOF's only appearing at element 1
      !
      ! Now we can start to integrate using this.
      
      ! Get the length of the current edge. It serves as a "determinant"
      ! in the cubature, so we have to divide it by 2 as an edge on the unit inverval
      ! [-1,1] has length 2.
      ivt1 = p_IverticesAtEdge (1,IMT)
      ivt2 = p_IverticesAtEdge (2,IMT)
      dedgelength = &
        SQRT ((p_DcornerCoordinates(1,ivt2)-p_DcornerCoordinates(1,ivt1))**2 &
            + (p_DcornerCoordinates(2,ivt2)-p_DcornerCoordinates(2,ivt1))**2 )
      dedgeweight = 0.5_DP * dedgelength
      
      ! Compute the coefficient in front of the integral:
      ! < Ju,v > = sum_E max(gammastar*nu*h_E, gamma*h_E^2) int_E [grad u] [grad v] ds
      dcoeff = MAX(rconfig%dgammastar * rconfig%dnu * dedgelength, &
                  rconfig%dgamma * dedgelength**2)
      
      ! Now we have the values of the basis functions in all the cubature 
      ! points.
      !
      ! Integrate the jump over the edges. This calculates the local matrix.
      !
      ! Loop through the test basis functions
      DO IDOFE = 1,ndofTest
      
        ! Loop through the trial basis functions
        DO JDOFE = 1,ndofTrial
    
          dval = 0.0_DP
          
          ! Loop through the cubature points to calculate the integral
          ! of the jump. Note that for computing the jump, we have to
          ! look in the inverse order to the cubature points of the neighbour
          ! element!
          ! Remember that the values of the basis functions on the first element
          ! are in p_DbasXXXX (.,.,.,1), while those of the 2nd element are
          ! in p_DbasXXXX (.,.,.,3) (rather than in p_DbasXXXX (.,.,.,2))
          ! by the above construction!
          DO icubp = 1,ncubp
          
            ! [ grad phi ]   ( jump in the derivative of trial basis function)
            ! = [ (d/dx) phi  ,  (d/dy) phi ]
            dphidx = p_DbasTrial (JDOFE,DER_DERIV_X,icubp,1) &
                  - p_DbasTrial (JDOFE,DER_DERIV_X,ncubp-icubp+1,3)

            dphidy = p_DbasTrial (JDOFE,DER_DERIV_Y,icubp,1) &
                  - p_DbasTrial (JDOFE,DER_DERIV_Y,ncubp-icubp+1,3)

            ! [ grad psi ]   ( jump in the derivative of test basis function)
            ! = [ (d/dx) phi  ,  (d/dy) phi ]
            dpsidx = DbasTest (IDOFE,DER_DERIV_X,icubp,1) &
                  - DbasTest (IDOFE,DER_DERIV_X,ncubp-icubp+1,3)

            dpsidy = DbasTest (IDOFE,DER_DERIV_Y,icubp,1) &
                  - DbasTest (IDOFE,DER_DERIV_Y,ncubp-icubp+1,3)
            
          
            ! Compute int_edge ( [grad phi_i] [grad phi_j] )
            dval = dval + Domega(icubp) * dedgeweight * &
                          (dphidx*dpsidx + dphidy*dpsidy)
          
          END DO

          ! Add the contribution to the local matrix -- weighted by the
          ! Omega from the cubature formula and the length of the edge.
          Dentry ((IDOFE-1)*ndofTrial+JDOFE) = &
            Dentry ((IDOFE-1)*ndofTrial+JDOFE) + dcoeff*dval

        END DO
      
      END DO
      
      ! For cdef= containing CONV_MODMATRIX, incorporate our "local" system matrix
      ! into the global matrix. The position of each entry DENTRY(X,Y)    
      ! in the global matrix array A was saved in element Kentry(X,Y)
      ! before.                                                      
      ! Kentry gives the position of the additive contributions in Dentry.
      ! The entry is weighted by the current dtheta, which is usually
      ! the weighting parameter of the corresponding THETA-scheme of a
      ! nonstationary simulation. For stationary simulations, dtheta is typically
      ! 1.0 which includes the local matrix into the global one directly.)
      
      IF (IAND(cdef,CONV_MODMATRIX) .NE. 0) THEN
        DO IDOFE = 0,ndofTest-1
          DO JDOFE = 1,ndofTrial
          
            p_Da(Kentry(IDOFE*ndofTrial+JDOFE)) = &
              p_Da(Kentry(IDOFE*ndofTrial+JDOFE)) + &
              rconfig%dtheta * Dentry (IDOFE*ndofTrial+JDOFE)
          
          END DO
        END DO
      END IF
      
      ! For cdef containing CONV_MODDEFECT, build the defect vector                     
      !     D = RHS - A*U                                         
      ! This is done matrix free, only with the help of the local 
      ! matrix.                                                   
      ! In this case, D=(D1,D2) is expected to be the RHS on      
      ! entry and will be updated to be the defect vector when    
      ! this routine is left.                                     
      
      IF (IAND(cdef,CONV_MODDEFECT) .NE. 0) THEN
        DO IDOFE=0,ndofTest-1

          irow=IdofsTest(1+IDOFE)

          DO JDOFE=1,ndofTrial

            dval = rconfig%dtheta * Dentry(IDOFE*ndofTrial+JDOFE)         
  
            jcol=p_IdofsTrial(JDOFE)
            Ddef1(irow)= Ddef1(irow) - dval*Du1(jcol)
            Ddef2(irow)= Ddef2(irow) - dval*Du2(jcol)

          END DO
        END DO
      END IF
      
      ! Proceed with next edge
    
    END DO ! IMT

    ! Clean up allocated arrays and memory.
    DEALLOCATE(DcubPtsRefOnAllEdges)
    CALL domint_doneIntegration(rintSubset)
    
    DEALLOCATE(DbasTest)
    DEALLOCATE(DbasTrial)

    DEALLOCATE(Kentry) 
    DEALLOCATE(Dentry)

    DEALLOCATE(IdofsTestTempl) 
    DEALLOCATE(IdofsTrialTempl)

  END SUBROUTINE

END MODULE
