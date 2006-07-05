!##############################################################################
!# ****************************************************************************
!# <name> convection </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains discretisation and application routines for basic 
!# convection: Upwind and streamline diffusion as used in CCxD / PPxD.
!#
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
    ! Standard value = 0.1.
    REAL(DP) :: dupsam = 0.1_DP
    
    ! Whether the viscosity is constant.
    LOGICAL :: bconstViscosity = .TRUE.
    
    ! Viscosity parameter $\nu = 1/Re$ if viscosity is constant
    REAL(DP) :: dnu = 1.0_DP
    
    ! Type of cubature formula to use. One of the CUB_xxxs constants from
    ! the cubature.f90 module. Standard is Gauss-2x2 for 2D equations.
    INTEGER :: ccubType = CUB_G2X2
    
    ! Weighting factor of the convective operator: $\theta * u*grad(u)$. 
    ! For time-dependent problems, this can be set to the step size
    ! in the $\Theta$-scheme.
    REAL(DP) :: dtheta = 1.0_DP
    
    ! Whether to use the ALE method for computing the convective operator.
    LOGICAL :: bALE = .FALSE.
    
  END TYPE
  
!</typeblock>

!</types>


!</constants>


CONTAINS

  ! ***************************************************************************

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
  ! if the ALE method should be used(bALE=true). In this case, the nonlinear
  ! term is modified to include the mesh velocity.\\
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
  !   too. By setting IDEF=1,2, this routine incorporates the nonlinearity
  !   into that vector, using the formula
  !
  !            $$ D = D - dtheta * UUx * grad (Ux) $$
  !   
  !  4.) This routine makes no use of COMMON-blocks about the
  !   triangulation; all information must be passed in the TRIA structure.
  !   The triangulation information KVERT,KMID,DCORVG must correspond
  !   to the arrays identified by the handles in TRIA, but are passed
  !   separately for faster access.\\
  !
  !  5.) If bALE=true, a mesh velocity field is added to the nonlineareity
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
    !    [F. Schieweck, Parallele Lösung der stationären inkompressiblen
    !     Navier-Stokes Gleichungen, Habilitation, Fakultät für
    !     Mathematik, Otto-von-Guericke-Universität Magdeburg]
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
    !     What is Lambda? 0<Lambda<1 is chosen depending on the flow
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
      ! which arises during a nonlinear iteration. If IDEF=1,2,
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

END MODULE
