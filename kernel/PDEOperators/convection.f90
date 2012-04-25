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
!# 2.) conv_streamlineDiffusion2d
!#     -> Apply streamline diffusion to a vector, a matrix or both.
!#
!# 3.) conv_streamlineDiffusion3d
!#     -> Apply streamline diffusion to a vector, a matrix or both.
!#
!# 4.) conv_JumpStabilisation1d,
!#     conv_JumpStabilisation2d,
!#     conv_JumpStabilisation3d
!#     -> Apply jump stabilisation to a vector, a matrix or both.
!#
!# 5.) conv_streamlineDiffusionBlk2d
!#     -> Apply streamline diffusion to a block vector, a block matrix or both.
!#        Extended method compared to conv_streamlineDiffusion2d.
!#        Allows to assemble the Newton matrix directly.
!#
!# 6.) conv_streamlineDiffusionBlk3d
!#     -> Apply streamline diffusion to a block vector, a block matrix or both.
!#        Extended method compared to conv_streamlineDiffusion3d.
!#        Allows to assemble the Newton matrix directly.
!#
!# 7.) conv_streamDiff2Blk2dMat
!#     -> Impose the SD operator into a matrix. 2D variant.
!#        New implementation, element independent.
!#
!# 8.) conv_streamDiff2Blk2dDef
!#     -> Impose the SD operator into a defect vector. 2D variant.
!#        New implementation, element independent.
!# </purpose>
!##############################################################################

module convection

!$use omp_lib
  use fsystem
  use genoutput
  use basicgeometry
  use bilinearformevaluation
  use collection
  use cubature
  use derivatives
  use dofmapping
  use domainintegration
  use element
  use elementpreprocessing
  use extstdassemblyinfo
  use feevaluation
  use geometryaux
  use jumpstabilisation
  use linearalgebra
  use linearsystemblock
  use linearsystemscalar
  use matrixio
  use perfconfig
  use spatialdiscretisation
  use statistics
  use storage
  use transformation
  use triangulation
  use vectorio
  use extstdassemblyinfo

  implicit none

  private

!<constants>

!<constantblock description="Constants for cdef parameter of convection routine">

  ! Modify the matrix
  integer, parameter, public :: CONV_MODMATRIX = 2**0

  ! Set up defect vector
  integer, parameter, public :: CONV_MODDEFECT = 2**1

  ! Set up both, matrix and defect vector
  integer, parameter, public :: CONV_MODBOTH   = CONV_MODMATRIX+CONV_MODDEFECT

!</constantblock>

!<constantblock description="Constants that define the upwind type">

  ! Standard Samarskji or simple upwind.
  integer, parameter, public :: CONV_UPW_SAMARSKJI = 0

!</constantblock>

!<constantblock description="Constants that define the jump stabilisation type">

  ! No jump stabilisation
  integer, parameter, public :: CONV_JUMP_NONE        = 0

  ! Unified edge jump stabilisation
  integer, parameter, public :: CONV_JUMP_UNIFIEDEDGE = 1

  ! Reactive jump stabilisation
  integer, parameter, public :: CONV_JUMP_REACTIVE    = 2

!</constantblock>

!</constants>

!<types>

!<typeblock>

  ! Configuration block for standard UPWIND scheme.
  type t_convUpwind

    ! Type of upwind. One of the CONV_UPW_xxxx constants.
    ! CONV_UPW_SAMARSKJI: Standard Samarskji or simple-type upwind
    integer :: cupwType = CONV_UPW_SAMARSKJI

    ! Stabilisation parameter.
    ! If cupwType=CONV_UPW_SAMARSKJI:
    !  -1.0 = simple Upwind,
    !  >= 0:  Samarskji upwind with parameter <tex>$ \theta $</tex> dupsam.
    !         Standard value = 0.1.
    real(DP) :: dupsam = 0.1_DP

    ! Whether the viscosity is constant.
    logical :: bconstViscosity = .true.

    ! Viscosity parameter <tex>$ \nu = 1/Re $</tex> if viscosity is constant.
    ! We set this to infinity, what quickly leads to a program crash if the
    ! application does not initialise that properly!
    real(DP) :: dnu = SYS_INFINITY_DP

    ! Weighting factor of the convective operator: <tex>$ \theta * u*grad(u) $</tex>.
    ! For time-dependent problems, this can be set to the step size
    ! in the <tex>$ \Theta $</tex>-scheme.
    real(DP) :: dtheta = 1.0_DP

    ! Whether to use the ALE method for computing the convective operator.
    logical :: bALE = .false.

  end type

  public :: t_convUpwind

!</typeblock>

!<typeblock>

  ! Configuration block for Streamline Diffusion discretisation scheme
  type t_convStreamlineDiffusion

    ! Stabilisation parameter.
    ! Standard value = 1.0_DP
    ! Note: A value of 0.0_DP sets up the convection part without
    ! any stabilisation (central-difference like discretisation).
    real(DP) :: dupsam = 1.0_DP

    ! Whether the viscosity is constant.
    logical :: bconstViscosity = .true.

    ! Viscosity parameter <tex>$ \nu = 1/Re $</tex> if viscosity is constant
    real(DP) :: dnu = 1.0_DP

    ! Weighting factor for the mass matrix.
    ! =0.0: do not incorporate mass matrix into the operator.
    real(DP) :: dalpha = 0.0_DP

    ! Weighting factor for the Stokes matrix. (Stokes matrix = dnu*Laplace)
    ! =0.0: do not incorporate Stokes matrix into the operator.
    real(DP) :: dbeta = 0.0_DP

    ! Weighting factor for the convective part.
    ! =0.0: do not incorporate convective part (=Stokes problem)
    ! =1.0: incorporate full convective part (=Navier Stokes problem)
    real(DP) :: ddelta = 1.0_DP

    ! Weighting factor for the transposed convective part.
    ! =0.0: do not incorporate transposed convective part
    ! =1.0: incorporate full transposed convective part
    real(DP) :: ddeltaTransposed = 0.0_DP

    ! Weighting factor of the complete operator.
    ! For time-dependent problems, this can be set to the step size
    ! in the <tex>$ \Theta $</tex>-scheme. For stationary problems, 1.0_DP must
    ! be used to assembly the not-weighted operator matrix.
    real(DP) :: dtheta = 1.0_DP

    ! Weighting factor for the Newton matrix (Frechet derivative <tex>$ (\Nabla u \cdot,\cdot) $</tex> of
    ! the convective operator <tex>$ u\Nabla u $</tex>, used for preconditioning).
    ! A value of 0.0 deactivates the Newton matrix.
    real(DP) :: dnewton = 0.0_DP

    ! Weighting factor for the transposed Newton matrix ((\Nabla u)^t\cdot)
    ! A value of 0.0 deactivates the transposed Newton matrix.
    real(DP) :: dnewtonTransposed = 0.0_DP

    ! Calculation of local H.
    ! =0: 2D: Use the root of the area of the hexahedron as local H
    !     3D: Use the cube-root of the volume of the hexahedron as local H
    ! =1: Use the length of the way that a particle travels through
    !     the hexahedron in direction of the flow
    integer :: clocalH = 1

    ! Whether to use the ALE method for computing the convective operator.
    logical :: bALE = .false.

  end type

  public :: t_convStreamlineDiffusion

!</typeblock>

!<typeblock>

  ! Configuration block for Streamline Diffusion discretisation scheme
  type t_convStreamDiff2

    ! Type of SD method to apply.
    ! = 0: Use simple SD stabilisation
    ! = 1: Use Samarskji SD stabilisation
    integer :: cstabilType = 1

    ! Weighting factor of the complete operator.
    ! For time-dependent problems, this can be set to the step size
    ! in the <tex>$ \Theta $</tex>-scheme. For stationary problems, 1.0_DP must
    ! be used to assembly the not-weighted operator matrix.
    real(DP) :: dtheta = 1.0_DP

    ! Stabilisation parameter.
    ! Standard value = 1.0_DP
    ! Note: A value of 0.0_DP sets up the convection part without
    ! any stabilisation (central-difference like discretisation).
    real(DP) :: dupsam = 0.0_DP

    ! Whether the viscosity dnu in front of the Laplace operator
    ! is constant.
    ! A value of .true. will use the coefficient dalpha below while
    ! a value .false. will lead to a computation of alpha using the
    ! callback routine.
    logical :: bconstNu = .true.

    ! Viscosity parameter <tex>$ \nu = 1/Re $</tex> if viscosity is constant.
    real(DP) :: dnu = 1.0_DP

    ! Whether the mass matrix coefficient dalpha is constant.
    ! A value of .true. will use the coefficient dalpha below while
    ! a value .false. will lead to a computation of alpha using the
    ! callback routine.
    logical :: bconstAlpha = .true.

    ! Weighting factor for the mass matrix M.
    ! =0.0: do not incorporate mass matrix into the operator.
    ! This factor is multiplied to a probably nonconstant alpha
    ! specified by the callback routine!
    real(DP) :: dalpha = 0.0_DP

    ! Weighting factor for the Stokes matrix nu*div(grad(.)) = nu*Laplace.
    ! =0.0: do not incorporate Stokes matrix into the operator.
    ! This factor is multiplied to a probably nonconstant dnu
    ! specified by the callback routine!
    real(DP) :: dbeta = 0.0_DP

    ! Weighting factor for the transposed Stokes matrix
    ! nu*div(grad(.)^T) which is typically used in the deformatino tensor.
    ! =0.0: do not incorporate.
    ! This factor is multiplied to a probably nonconstant dnu
    ! specified by the callback routine!
    real(DP) :: dbetaT = 0.0_DP

    ! Weighting factor for the convective part (u*grad(.).,.)
    ! =0.0: do not incorporate convective part (=Stokes problem)
    ! =1.0: incorporate full convective part (=Navier Stokes problem)
    real(DP) :: ddelta = 0.0_DP

    ! Weighting factor for the adjoint convective part (.,u*gred(.).).
    ! This is the adjoint operator to the convection operator above,
    ! i.e., the convection is applied to the test space.
    ! =0.0: do not incorporate convective part (=Stokes problem)
    ! =1.0: incorporate full convective part (=Navier Stokes problem)
    real(DP) :: ddeltaAdj = 0.0_DP

    ! Weighting factor for the transposed convective part ((u^T)*grad(.).,.)
    ! =0.0: do not incorporate transposed convective part
    ! =1.0: incorporate full transposed convective part
    real(DP) :: ddeltaT = 0.0_DP

    ! Weighting factor for the adjoint transposed convective part (.,(u^T)*grad(.).)
    ! =0.0: do not incorporate transposed convective part
    ! =1.0: incorporate full transposed convective part
    real(DP) :: ddeltaTadj = 0.0_DP

    ! Weighting factor for the Newton matrix ((grad(u))*(.),.)
    ! (-> Frechet derivative <tex>$ \cdot\Nabla u $</tex> of
    ! the convective operator <tex>$ u\Nabla u $</tex>, used for preconditioning).
    ! A value of 0.0 deactivates the Newton matrix.
    real(DP) :: dnewton = 0.0_DP

    ! Weighting factor for the adjoint Newton matrix 
    ! (Frechet derivative <tex>$ (\Nabla u \cdot,\cdot) $</tex> of
    ! the convective operator <tex>$ u\Nabla u $</tex> with test and trial
    ! functions exchanged, used for preconditioning).
    ! A value of 0.0 deactivates the Newton matrix.
    real(DP) :: dnewtonAdj = 0.0_DP

    ! Weighting factor for the transposed Newton matrix ((grad(u)^T)*(.),.)
    ! ((\Nabla u)^t\cdot)
    ! A value of 0.0 deactivates the transposed Newton matrix.
    real(DP) :: dnewtonT = 0.0_DP

    ! Calculation of local H.
    ! =0: 2D: Use the root of the area of the hexahedron as local H
    !     3D: Use the cube-root of the volume of the hexahedron as local H
    ! =1: Use the length of the way that a particle travels through
    !     the hexahedron in direction of the flow
    integer :: clocalH = 1

  end type

  public :: t_convStreamDiff2

!</typeblock>

!<typeblock>

  ! Parameter block for Jump stabilisation
  type t_jumpStabilisation

    ! Whether the viscosity is constant.
    logical :: bconstViscosity = .true.

    ! Viscosity parameter <tex>$ \nu = 1/Re $</tex> if viscosity is constant
    real(DP) :: dnu            = 1.0_DP

    ! Line integral cubature formula for discretising the Jump.
    ! One of the CUB_xxxx_1D-constants of the cubature.f90 module.
    ! Standard is Gauss 2-point formula.
    integer(I32) :: ccubType = CUB_G2_1D

    ! Type of Jump stabilisation.
    ! One of the CONV_JUMP_xxxx-constants. Standard is unified edge
    ! stabilisation.
    ! CONV_JUMP_UNIFIEDEDGE: <tex>$ \sum_E \gamma h_E^2 \int_E [grad u] [grad v] ds $</tex>
    ! CONV_JUMP_REACTIVE:    <tex>$ \sum_E gamma nu 1/|E| \int_E [u] [v] ds $</tex>
    integer :: cjump = CONV_JUMP_UNIFIEDEDGE

    ! 1st Relaxation parameter for the Jump stabilisation.
    ! =0.0: No Jump stabilisation
    ! >0.0: Add Jump stabilisation with relaxation parameter dgamma=djump.
    real(DP) :: dgamma = 0.01_DP

    ! 2nd Relaxation parameter for the Jump stabilisation.
    ! =0.0: No Jump stabilisation (Standard)
    ! >0.0: Add Jump stabilisation with relaxation parameter dgammastar
    real(DP) :: dgammastar = 0.0_DP

    ! Exponent for edge length weight in the jump stabilisation.
    ! A value of 2 corresponds to a weight h_E^2, but this can be changed here.
    real(dp) :: deojEdgeExp = 2.0_DP

    ! Weighting factor of the complete operator.
    ! For time-dependent problems, this can be set to the step size
    ! in the <tex>$ \Theta $</tex>-scheme. For stationary problems, 1.0_DP must
    ! be used to assembly the not-weighted operator matrix.
    real(DP) :: dtheta = 1.0_DP

  end type

  public :: t_jumpStabilisation

!</typeblock>

!</types>

  !************************************************************************

  ! global performance configuration
  type(t_perfconfig), target, save :: conv_perfconfig

  !************************************************************************

  public :: conv_upwind2d
  public :: conv_streamlineDiffusion2d
  public :: conv_streamlineDiffusion3d
  public :: conv_jumpstabilisation1d
  public :: conv_JumpStabilisation2d
  public :: conv_jumpstabilisation3d
  public :: conv_streamlineDiffusionBlk2d
  public :: conv_streamlineDiffusionBlk3d
  public :: conv_streamDiff2Blk2dMat
  public :: conv_streamDiff2Blk2dDef

contains

  ! ***************************************************************************

!<subroutine>

 subroutine conv_upwind2d (rvecPrimary, rvecSecondary, dprimWeight, dsecWeight,&
                           rconfig, cdef, &
                           rmatrix, rsolution, rdefect, DmeshVelocity, &
                           IvelocityComp, rperfconfig)

!<description>
  ! Standard 1st order upwinding method to set up the convection operator
  !
  !            <tex> $$ u_1 * grad(u_2) $$ </tex>
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
  !  <tex> $$ u_1  =  dprimWeight * rvecPrimary  +  dsecWeight * rvecSecondary $$ </tex>
  !
  ! <tex>$ u_2 = rsolution(.) $</tex> defines the second velocity field inside of
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

  ! Primary velocity field for the computation of <tex>$ u_1 $</tex>
  type(t_vectorBlock), intent(in), target :: rvecPrimary

  ! Secondary velocity field for the computation of <tex>$ u_1 $</tex>
  type(t_vectorBlock), intent(in), target :: rvecSecondary

  ! Weighting factor for rvecPrimary.
  real(DP), intent(in) :: dprimWeight

  ! Weighting factor for rvecSecondary.
  real(DP), intent(in) :: dsecWeight

  ! Configuration block for the upwind scheme
  type(t_convUpwind), intent(in) :: rconfig

  ! Computation/defect correction method. One of the CONV_MODxxxx constants:
  ! CONV_MODMATRIX: Set up the nonlinear matrix. rmatrix must be present, the
  !                 nonlinear part is added to the matrix.
  ! CONV_MODDEFECT: Set up the nonlinear defect. rdefect and rsolution must be
  !                 present.
  ! CONV_MODBOTH  : Set up the nonlinear matrix as well as the nonlinear defect.
  !                 rmatrix, rdefect and rsolution must all be present.
  integer, intent(in) :: cdef

  ! OPTIONAL: Solution vector u_2.
  ! Must be present if cdef=CONV_MODDEFECT or =CONV_MODBOTH.
  type(t_vectorBlock), intent(in), target, optional :: rsolution

  ! OPTIONAL: Mesh velocity field.
  ! DmeshVelocity(1,ivt) gives the X-velocity of the mesh, i.e. the X-velocity
  !   of the corner vertex ivt.
  ! DmeshVelocity(2,ivt) gives the Y-velocity of the mesh, i.e. the Y-velocity
  !   of the corner vertex ivt.
  ! The parameter must be present if ALE is activated in the
  ! configuration parameter block.
  real(DP), dimension(:,:), intent(in), optional :: DmeshVelocity

  ! OPTIONAL:
  ! Index block that specifies which component in rvecPrimary / rvecSecondary /
  ! rsolution / rdefect is the X-velocity and which one is the Y-velocity.
  !  IvelocityComp(1) gives the number of the X-velocity (usually = 1),
  !  IvelocityComp(2) gives the number of the Y-velocity (usually = 2).
  ! If not present, IvelocityComp=(/1,2/) is assumed, thus the X-velocity
  ! must be in rsolution\%RvectorBlock(1) and the Y-velocity in
  ! rsolution\%RvectorBlock(2).
  integer, dimension(2), intent(in), optional :: IvelocityComp

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>

!<inputoutput>
  ! System matrix.
  ! The content of the matrix must be present if cdef=CONV_MODMATRIX or
  ! =CONV_MODBOTH, otherwise only the structure is used.
  ! The nonlinear operator is added to the matrix.
  type(t_matrixScalar), intent(inout) :: rmatrix

  ! OPTIONAL: Defect vector.
  ! Must have the same structure as rsolution/rvecPrimary/rvecSecondary.
  ! Must be present if cdef=CONV_MODDEFECT or =CONV_MODBOTH.
  ! The nonlinear part is subtracted from this vector:
  ! <tex>$ r = r - \theta * u_1*grad(u_2) $</tex>
  type(t_vectorBlock), intent(inout), optional, target :: rdefect
!</inputoutput>

!</subroutine>

    ! local variables
    integer(I32) :: celement
    integer, dimension(2) :: Icomp
    type(t_vectorScalar), pointer :: p_rvelX1,p_rvelX2,p_rvelY1,p_rvelY2
    type(t_vectorScalar), pointer :: p_rsolX,p_rsolY,p_rdefectX,p_rdefectY
    real(DP), dimension(:), pointer :: p_DvelX1,p_DvelX2,p_DvelY1,p_DvelY2
    real(DP), dimension(:), pointer :: p_DsolX,p_DsolY,p_DdefectX,p_DdefectY

    ! At first check the input parameters that everything is present what
    ! we need:
    if ((cdef .eq. CONV_MODDEFECT) .or. (cdef .eq. CONV_MODBOTH)) then
      if ((.not. present(rsolution)) .or. (.not. present(rdefect))) then
        call output_line ("Solution/defect vector not present!", &
            OU_CLASS_ERROR,OU_MODE_STD,"conv_upwind2d")
        call sys_halt()
      end if
    end if

    if (rconfig%bALE) then
      if (.not. present(DmeshVelocity)) then
        call output_line ("Mesh velocity vector not present!", &
            OU_CLASS_ERROR,OU_MODE_STD,"conv_upwind2d")
        call sys_halt()
      end if
    end if

    ! Get the actual subvectors from the velocity vectors that define
    ! the X- and Y-velocity.
    if (present(IvelocityComp)) then
      Icomp = IvelocityComp
    else
      Icomp = (/1,2/)
    end if

    p_rvelX1 => rvecPrimary%RvectorBlock(Icomp(1))
    p_rvelY1 => rvecPrimary%RvectorBlock(Icomp(2))
    p_rvelX2 => rvecSecondary%RvectorBlock(Icomp(1))
    p_rvelY2 => rvecSecondary%RvectorBlock(Icomp(2))

    if (present(rsolution)) then
      p_rsolX => rsolution%RvectorBlock(Icomp(1))
      p_rsolY => rsolution%RvectorBlock(Icomp(2))
    else
      nullify(p_rsolX)
      nullify(p_rsolY)
    end if

    if (present(rdefect)) then
      p_rdefectX => rdefect%RvectorBlock(Icomp(1))
      p_rdefectY => rdefect%RvectorBlock(Icomp(2))
    else
      nullify(p_rdefectX)
      nullify(p_rdefectY)
    end if

    ! At the moment, we only support a rather limited set of configurations:
    ! Matrix and vectors must all be double precision, matrix must be format
    ! 7 or 9, discretisation must be Q1~, constant viscosity.
    if ((rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
        (rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
      call output_line ("Unsupported matrix format!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_upwind2d")
      call sys_halt()
    end if

    celement = rmatrix%p_rspatialDiscrTest%RelementDistr(1)%celement
    if ((rmatrix%p_rspatialDiscrTest%ccomplexity .ne. SPDISC_UNIFORM) .or. &
        (elem_getPrimaryElement(celement) .ne. EL_Q1T)) then
      call output_line ("Unsupported discretisation!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_upwind2d")
      call sys_halt()
    end if

    if ((rvecPrimary%cdataType .ne. ST_DOUBLE) .or. &
        (rvecSecondary%cdataType .ne. ST_DOUBLE)) then
      call output_line ("Unsupported vector data type in velocity!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_upwind2d")
      call sys_halt()
    end if

    if (present(rdefect)) then
      if ((rsolution%cdataType .ne. ST_DOUBLE) .or. &
          (rdefect%cdataType .ne. ST_DOUBLE)) then
        call output_line ("Unsupported vector data type in solution/defect!", &
            OU_CLASS_ERROR,OU_MODE_STD,"conv_upwind2d")
        call sys_halt()
      end if
    end if

    if (.not. rconfig%bconstViscosity) then
      call output_line ("Only constant viscosity supported at the moment!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_upwind2d")
      call sys_halt()
    end if

    if (rconfig%dnu .eq. SYS_INFINITY_DP) then
      call output_line ("Viscosity parameter nu not initialised!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_upwind2d")
      call sys_halt()
    end if

    ! Call the actual calculation routine.
    ! Hide the p_rsol...-parameters to prevent passing the NULL()-pointer
    ! if rsolution is not present -- some compilers do not like that ^^

    call lsyssc_getbase_double (p_rvelX1,p_DvelX1)
    call lsyssc_getbase_double (p_rvelY1,p_DvelY1)
    call lsyssc_getbase_double (p_rvelX2,p_DvelX2)
    call lsyssc_getbase_double (p_rvelY2,p_DvelY2)

    !!! DEBUG:
    !WHERE (abs(p_DvelX1) .LT. 1E-12_DP) p_DvelX1 = 0.0_DP
    !WHERE (abs(p_DvelY1) .LT. 1E-12_DP) p_DvelY1 = 0.0_DP
    !call vecio_writeArray_Dble (p_DvelX1, 'vecx1', &
    !                               0, 'vectorx1.txt', '(D10.3)')
    !call vecio_writeArray_Dble (p_DvelY1, 'vecx2', &
    !                               0, 'vectorx2.txt', '(D10.3)')

    if (present(rdefect)) then
      call lsyssc_getbase_double (p_rsolX   ,p_DsolX   )
      call lsyssc_getbase_double (p_rsolY   ,p_DsolY   )
      call lsyssc_getbase_double (p_rdefectX,p_DdefectX)
      call lsyssc_getbase_double (p_rdefectY,p_DdefectY)

      call conv_upwind2dALE_Q1Tdouble ( &
                    p_DvelX1,p_DvelY1,p_DvelX2,p_DvelY2,dprimWeight,dsecWeight, &
                    rmatrix,rmatrix%p_rspatialDiscrTest%p_rtriangulation, &
                    cdef, rconfig%dupsam, rconfig%dnu, rconfig%dtheta, &
                    rconfig%bALE, &
                    p_DsolX,p_DsolY,p_DdefectX,p_DdefectY, DmeshVelocity, rperfconfig)

    else

      call conv_upwind2dALE_Q1Tdouble ( &
                    p_DvelX1,p_DvelY1,p_DvelX2,p_DvelY2,dprimWeight,dsecWeight, &
                    rmatrix,rmatrix%p_rspatialDiscrTest%p_rtriangulation, &
                    cdef, rconfig%dupsam, rconfig%dnu, rconfig%dtheta, &
                    rconfig%bALE, DmeshVelocity=DmeshVelocity, rperfconfig=rperfconfig)

      !!! DEBUG:
      !call matio_writeMatrixHR (rmatrix, 'matrix',&
      !                          .TRUE., 0, 'matrixL.txt', '(D10.3)')

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine conv_upwind2dALE_Q1Tdouble ( &
                  u1Xvel,u1Yvel,u2Xvel,u2Yvel,dweight1,dweight2,&
                  rmatrix,rtriangulation,cdef, &
                  dupsam,dnu,dtheta,bALE, &
                  Du1,Du2,Ddef1,Ddef2,DmeshVelocity,rperfconfig)
!<description>
  ! Standard 1st order upwinding method to set up the convection operator
  !
  !            <tex> $$ u_1 * grad(u_2) $$ </tex>
  !
  ! in a matrix or to build a defect vector.
  ! 2D-version (X- and Y-velocity), uniform <tex>$\tilde Q_1$</tex> discretisation,
  ! double precision vectors/matrix.
  !
  ! u1Xvel,u1Yvel, u2Xvel,u2Yvel are two velocity field vectors,
  ! (u1Xvel,u1Yvel) a primary and (u2Xvel,u2Yvel) a secondary velocity field.
  ! The final velocity vector field is then computed as a weighted average
  ! of these two:
  !
  !  <tex> $$ u_1  =  dweight1 * u1vel  +  dweight2 * u2vel $$ <tex>
  !
  ! <tex>$ u_2 = rsolution(.) $</tex> defines the second velocity field inside of
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
  ! [Duarte, Formaz, Natesan; `Arbitrary Lagrangian-Euler Method
  ! for Navier-Stokes equations with moving boundaries`;
  ! Comput. Methods Appl. Mech. Engrg. 193 (2004), 4819-4836]
  !
  ! Remarks:\\
  !
  ! 1.) In a typical call of the upwinding, the caller can use:
  !     dweight1 = 1, u1Xvel/u1Yvel = velocity field
  !     dweight2 = 0, u2Xvel/u2Yvel = undefined
  !   So the upwinding scheme only uses one velocity field.
  !   Such a call e.g. adds the integral
  !                <tex> $$ ( u_1 * grad(.) , v )_{\Omega} $$ </tex>
  !   to the system matrix.\\
  !
  !  2.) In case that there are two velocity fields representing
  !   the solution (may happen in a nonstationary simulation where
  !   u1Xvel/u1Yvel represents the solution in the current and u2Xvel/u2Yvel
  !   that of the previous time step), dweight1/dweight2 defines how these both
  !   velocity vectors should be weighted to compute the actual
  !   velocity field for the assembling:
  !               <tex> $$ U_act = dweight1*u1vel + dweight2*u2vel $$ </tex>
  !   This is e.g. used for the linear extrapolation technique to
  !   reconstruct a velocity from two previous time steps...\\
  !
  !  3.) In the nonlinear iteration, as a right hand side there arises
  !   a defect vector D, which linear part can easily being assembled.
  !   However, there is a nonlinearity to be included into that vector,
  !   too. By setting cdef=1,2, this routine incorporates the nonlinearity
  !   into that vector, using the formula
  !
  !            <tex> $$ D = D - dtheta * UUx * grad (Ux) $$ </tex>
  !
  !  4.) If bALE=true, a mesh velocity field is added to the nonlineareity
  !   according to the formula  "U * grad (U-DmeshVelocity)".
  !   For bALE=false, the simple nonlinearity "U * grad (U)" is used.

!</description>

!<input>

  ! Primary X-velocity of <tex>$ u_1 $</tex>
  real(DP), dimension(:), intent(in) :: u1Xvel

  ! Primary Y-velocity of <tex>$ u_1 $</tex>
  real(DP), dimension(:), intent(in) :: u1Yvel

  ! Secondary X-velocity of <tex>$ u_1 $</tex>
  real(DP), dimension(:), intent(in) :: u2Xvel

  ! Secondary Y-velocity of <tex>$ u_1 $</tex>
  real(DP), dimension(:), intent(in) :: u2Yvel

  ! Computation/defect correction method. One of the CONV_MODxxxx constants:
  ! CONV_MODMATRIX: Set up the nonlinear matrix. rmatrix must be present, the
  !                 nonlinear part is added to the matrix.
  ! CONV_MODDEFECT: Set up the nonlinear defect. rdefect and rsolution must be
  !                 present.
  ! CONV_MODBOTH  : Set up the nonlinear matrix as well as the nonlinear defect.
  !                 rmatrix, rdefect and rsolution must all be present.
  integer, intent(in) :: cdef

  ! Weighting factor for u1Xvel/u1Yvel.
  real(DP), intent(in) :: dweight1

  ! Weighting factor for u2Xvel/u2Yvel.
  real(DP), intent(in) :: dweight2

  ! dupsam  - control parameter.
  !          -1: simple upwind,
  !          =0: Samarskji upwind
  real(DP), intent(in) :: dupsam

  ! Viscosity parameter <tex>$ \nu = 1/Re $</tex> if viscosity is constant
  real(DP), intent(in) :: dnu

  ! Weighting factor of the convective operator: <tex>$ \theta * u*grad(u) $</tex>.
  ! For time-dependent problems, this can be set to the step size
  ! in the <tex>$ \Theta $</tex>-scheme.
  real(DP), intent(in) :: dtheta

  ! Whether or not to use the ALE method
  logical, intent(in) :: bALE

  ! OPTIONAL: Mesh velocity field. Must be present if bALE=TRUE.
  ! DmeshVelocity(1,:) gives the X-velocity of all the corner points of the mesh,
  ! DmeshVelocity(2,:) gives the Y-velocity.
  real(DP), dimension(:,:), intent(in), optional :: DmeshVelocity(:,:)

  ! Triangulation structure specifying the underlying mesh.
  type(t_triangulation), intent(in) :: rtriangulation

  ! OPTIONAL: X-velocity of <tex>$ u_2 $</tex>. Must be present if cdef=CONV_MODDEFECT
  ! or cdef=CONV_MODBOTH.
  real(DP), dimension(:), intent(in), optional :: Du1

  ! Y-velocity of <tex>$ u_2 $</tex>. Must be present if cdef=CONV_MODDEFECT
  ! or cdef=CONV_MODBOTH.
  real(DP), dimension(:), intent(in), optional :: Du2

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>

!<inputoutput>
  ! The system matrix. Must be format 7 or 9.
  type(t_matrixScalar), intent(inout), target :: rmatrix

  ! OPTIONAL: X-defect vector. Must be present if cdef=CONV_MODDEFECT
  ! or =CONV_MODBOTH.
  real(DP), dimension(:), intent(inout), optional :: Ddef1

  ! OPTIONAL: Y-defect vector. Must be present if cdef=CONV_MODDEFECT
  ! or =CONV_MODBOTH.
  real(DP), dimension(:), intent(inout), optional :: Ddef2
!</inputoutput>

!</subroutine>

    ! local variables
    integer, dimension(4) :: Iedge
    integer :: IlocalMatrix(4,4)
    real(DP), dimension(4) :: Dflux(4),Duu1(4),Duu2(4),XV(4),YV(4)
    real(DP), dimension(4,4) :: DlocalMatrix(4,4)

    real(DP) :: dcenterX, dcenterY, XN, YN, G1, G2, dupsre, ELMH
    real(DP) :: DL0, DL2, H00, H22, dflux0, dflux2
    integer :: iel
    integer :: iv, ivt1,ivt2
    integer :: im1, im0, im2
    integer :: I,II
    integer ia1, ia2, J, JJ, ia
    real(DP), dimension(4) :: DuALE1,DuALE2

    real(DP), dimension(:), pointer :: p_Da
    integer, dimension(:), pointer :: p_Kcol
    integer, dimension(:), pointer :: p_Kld

    integer :: NVT
    integer, dimension(:,:), pointer :: p_Kvert
    integer, dimension(:,:), pointer :: p_Kmid
    real(DP), dimension(:,:), pointer :: p_Dcorvg

    ! There is no additional type/completeness check in the parameters here.
    ! This routine should never be called from outside, otherwise the main
    ! application may have problems if the parameters are not specified correctly!
    !
    ! Get a pointer to the matrix entries and the matrix structure
    if (rmatrix%h_DA .ne. ST_NOHANDLE) then
      ! The entries may be undefined - allowed if cdef=CONV_MODDEFECT!
      call lsyssc_getbase_double (rmatrix,p_Da)
    end if
    call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
    call lsyssc_getbase_Kld (rmatrix,p_Kld)

    ! Get a pointer to triangulation arrays.
    call storage_getbase_int2D (rtriangulation%h_IverticesAtElement,p_Kvert)
    call storage_getbase_int2D (rtriangulation%h_IedgesAtElement,p_Kmid)
    call storage_getbase_double2D (rtriangulation%h_DvertexCoords,p_Dcorvg)
    NVT = rtriangulation%NVT

    !********************************************************************
    !    Weighted Samarski upwind
    !
    !    This implementation follows the documentation of
    !    [F. Schieweck, Parallele Loesung der stationaeren inkompressiblen
    !     Navier-Stokes Gleichungen, Habilitation, Fakultaet fuer
    !     Mathematik, Otto-von-Guericke-Universitaet Magdeburg]
    !
    !********************************************************************
    !
    !********************************************************************
    !     What we want to discretise here is the term
    !
    !         n(z,u,v) = ( (z*grad (.))u , v )_Omega
    !
    !     Let us assume we have two elements next to each other:
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
    !     The edges Gl and Gk enclose a diagonal edge Glk of the above
    !     triangle. The operator can now be rewritten by decomposition
    !     into the elements as
    !
    !       n(z,u,v) ~= sum_l sum_k int_Glk (z*n_lk) (u-u(Bl)) v(Bl) dGamma
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
    !     What is the approach here? As said, u_upw is reconstructed
    !     from u_1 and u_2 by a simple mean formula:
    !
    !          u_upw = Lambda u_l  + (1-Lambda) u_k
    !
    !     What is Lambda? 0 < Lambda < 1 is chosen depending on the flow
    !     crossing the diagonal Glk. More precisely, it is chosen
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

    if (dupsam .ge. 0.0_DP) dupsre = dupsam / dnu

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

    do iel=1,rtriangulation%NEL

      ! dcenterX,dcenterY will be coordinates of the center of the element

      dcenterX=0.0_DP
      dcenterY=0.0_DP

      ! ALE handling on the current element.
      ! Loop over all 4 nodes to calculate the contribution of ALE:
      if (bALE) then
        do II=1,4
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
          ivt2 = p_Kvert(mod(II,4)+1,iel)

          ! And then calculate the mean velocity field:

          DuALE1(II) = 0.5_DP * ( DmeshVelocity (1,ivt1) + DmeshVelocity(1,ivt2) )
          DuALE2(II) = 0.5_DP * ( DmeshVelocity (2,ivt1) + DmeshVelocity(2,ivt2) )
        end do
      end if

      ! Loop over all 4 U-nodes.
      ! Calculate Iedge,XV,YV,dcenterX,dcenterY,Duu1,Duu2

      do II=1,4

        ! Get the number of the II-th edge - which is at the same time the
        ! DOF in the velocity vector(s).
        I = p_Kmid(II,iel)

        ! Store the number of the edge/DOF in Iedge:
        Iedge(II)=I

        ! Store the coordinates of the corner vertices of that
        ! element in XV/YV:

        iv=p_Kvert(II,iel)

        XV(II)=p_Dcorvg(1,iv)
        YV(II)=p_Dcorvg(2,iv)

        ! Sum up the coordinates of the element - will later result
        ! in the element midpoint:

        dcenterX=dcenterX+XV(II)
        dcenterY=dcenterY+YV(II)

        ! Now we want to compute the velocity on the edge II
        ! (following the node II).
        !
        ! Compute the actual velocity in the edge II (following the
        ! node II) by a weighted mean of the both velocity vectors
        ! U1Lx and U2Lx. This allows e.g. to reconstruct a velocity
        ! vector by linear extrapolation of two previous time steps.
        !
        ! Subtract the mesh-velocity on the current edge as described
        ! above; if ALE is not used, DuALE is = 0, so nothing happens.

        Duu1(II) = dweight1*u1Xvel(I)+dweight2*u2Xvel(I) - DuALE1(II)
        Duu2(II) = dweight1*u1Yvel(I)+dweight2*u2Yvel(I) - DuALE2(II)

      end do

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
      ! Initialise DlocalMatrix(.,.) to 0. DlocalMatrix will assemble the "local"
      ! matrix, i.e. the values that are later incorporated into
      ! the system matrix at the positions stored in IlocalMatrix.
      !
      ! Loop over all 4 U-nodes.
      ! Calculate Dflux(.), IlocalMatrix(.,.), DlocalMatrix(.,.)

      DlocalMatrix = 0.0_DP

      do II=1,4

        ! II1 corresponds to the current node in question.
        ! im1=II-1 MOD 4 receives the predecessor of II1 in an anticlockwise
        ! sense:
        im1=II-1
        if (im1.lt.1) im1=4

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

        dofsearch: do JJ=1,4

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

          do ia=ia1,ia2
            if (p_Kcol(ia) .eq. J) then
              ! Save the matrix index in IlocalMatrix(II,JJ) so we can find
              ! the matrix entry later without searching for it.

              IlocalMatrix(II,JJ)=ia

              ! Next DOF, no error
              cycle dofsearch

            end if
          end do

          ! Error case

          call output_line ("Entry index ia not found!", &
              OU_CLASS_ERROR,OU_MODE_STD,"conv_upwind2dALE_Q1Tdouble")
          return

        end do dofsearch ! JJ

      end do ! II

      ! What have we calculated up to here? Let us collect...
      !
      ! Dflux        - The flux along the edges of the triangles
      ! IlocalMatrix - The indices in the matrix A of the entries that
      !                are affected by integration on the current element iel.
      ! DlocalMatrix - Local element matrix - initialised with 0 up to now.
      !
      ! So the next step is to assemble the local element matrix,
      ! which is to be incorporated into the system matrix later.
      !
      ! Loop over the nodes to calculate DlocalMatrix:

      do II=1,4

        ! Set im1=predecessor of II, im2=successor of II,
        ! in counterclockwise sense.

        im0=II
        im1=II-1
        if (im1.lt.1) im1=4
        im2=II+1
        if (im2.gt.4) im2=1

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

        if (dupsam .ge. 0.0_DP) then

          ! The user wants Samarskji-Upwind.
          ! Take dupsre, the dupsam-parameter, weighted by 1/nu.
          ! Remember: In the previous calculation of the line-integral
          ! to calculate t, we did not incorporate 1/nu - this is
          ! repaired here:
          !
          ! Analyze the two fluxes on the edges of the triangle.
          ! Calculate the lambda-value by Samarskji-upwind,
          ! depending on whether the flux goes "in" or "out" of the
          ! triangle. DL0 saves the lambda of the triangle-edge
          ! corresponding to II, DL2 that of the triangle edge
          ! corresponding to  im2.

          if (dflux0 .ge. 0.0_DP) then
            DL0=PHIP(dupsre*dflux0)
          else
            DL0=PHIM(dupsre*dflux0)
          end if

          if (dflux2 .ge. 0.0_DP) then
            DL2=PHIP(dupsre*dflux2)
          else
            DL2=PHIM(dupsre*dflux2)
          end if

        else

          ! Simple Upwinding scheme.
          ! The corresponding lambda (wighting factor of the
          ! "adjacent velocities") is either 0 or 1, depending on
          ! whether the flux goes "into" or "out of" the triangle
          ! on that edge.

          DL0=0.0_DP
          DL2=0.0_DP
          if (dflux0 .ge. 0.0_DP) DL0=1.0_DP
          if (dflux2 .ge. 0.0_DP) DL2=1.0_DP

        end if

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

      end do ! II

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
      ! which is realised by a multiplication of the u*grad(.)-
      ! matrix by THETA*K = dtheta here!

      DlocalMatrix = dtheta*DlocalMatrix

      ! Then incorporate the local matrix into the global matrix
      ! and/or the global defect.

      if (iand(cdef,CONV_MODMATRIX) .ne. 0) then

        do JJ=1,4
          do II=1,4
            ia       = IlocalMatrix(II,JJ)
            p_Da(ia) = p_Da(ia) + DlocalMatrix(II,JJ)
          end do ! II
        end do ! JJ

      end if

      if (iand(cdef,CONV_MODDEFECT) .ne. 0) then

        do II=1,4
          do JJ=1,4
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
          end do ! JJ
        end do ! II

      end if

    end do ! IEL

  contains

    ! Auxiliary function 1
    elemental real(DP) function PHIP(x)
    real(DP), intent(in) :: x
      PHIP = (0.5_DP+x)/(1.0_DP+x)
    end function

    ! Auxiliary function 2
    elemental real(DP) function PHIM(x)
    real(DP), intent(in) :: x
      PHIM = 0.5_DP/(1.0_DP-x)
    end function

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine conv_streamlineDiffusion2d ( &
      rvecPrimary, rvecSecondary, dprimWeight, dsecWeight,&
      rconfig, cdef, &
      rmatrix, rsolution, rdefect, DmeshVelocity, rcubatureInfo, rperfconfig)

!<description>
  ! Standard streamline diffusion method to set up the operator
  ! <tex>
  !   $$ dtheta  *  ( dalpha * MASS  +  dbeta * STOKES  +  ddelta * u_1 * grad(u_2) ) $$
  ! </tex>
  ! in a matrix or to build a defect vector.
  ! 2D-version (X- and Y-velocity).
  !
  ! rvecPrimary, rvecSecondary are two velocity field vectors for the X-
  ! and Y-veclocity; IvelocityComp defines which components of these
  ! vectors contains the X- and which contains the Y-velocity.
  ! The final velocity vector field is then computed as a weighted average
  ! of these two:
  ! <tex>
  !  $$ u_1  =  dprimWeight * rvecPrimary  +  dsecWeight * rvecSecondary $$
  ! </tex>
  ! <tex>$ u_2 = rsolution(.) $</tex> defines the second velocity field inside of
  ! the grad-term.
  !
  ! The switch cdef decides on whether the routine sets up the nonlinear
  ! defect, the nonlinear matrix or both.
  !
  ! If the operator is to be applied to a defect vector, it is applied
  ! to all components in the vector!
  ! (Note: In a Navier-Stokes case e.g., the rdefect/rsolution vectors
  ! must therefore only contain the velocity components, not the pressure!
  ! If necessary, derive a subvector containing only the velocity blocks.)
  !
  ! rmeshVelocity is an optional mesh velocity field that must be present
  ! if the ALE method should be used.
  !
  ! The configuration how the routine should react is to be configured
  ! in the configuration block rconfig.
!</description>

!<input>

  ! Primary velocity field for the computation of <tex>$ u_1 $</tex>
  type(t_vectorBlock), intent(in), target :: rvecPrimary

  ! Secondary velocity field for the computation of <tex>$ u_1 $</tex>
  type(t_vectorBlock), intent(in), target :: rvecSecondary

  ! Weighting factor for rvecPrimary.
  real(DP), intent(in) :: dprimWeight

  ! Weighting factor for rvecSecondary.
  real(DP), intent(in) :: dsecWeight

  ! Configuration block for the streamline diffusion scheme
  type(t_convStreamlineDiffusion), intent(in) :: rconfig

  ! Computation/defect correction method. One of the CONV_MODxxxx constants:
  ! CONV_MODMATRIX: Set up the nonlinear matrix. rmatrix must be present, the
  !                 nonlinear part is added to the matrix.
  ! CONV_MODDEFECT: Set up the nonlinear defect. rdefect and rsolution must be
  !                 present.
  ! CONV_MODBOTH  : Set up the nonlinear matrix as well as the nonlinear defect.
  !                 rmatrix, rdefect and rsolution must all be present.
  integer, intent(in) :: cdef

  ! OPTIONAL: Solution vector u_2.
  ! Must be present if cdef=CONV_MODDEFECT or =CONV_MODBOTH.
  type(t_vectorBlock), intent(in), target, optional :: rsolution

  ! OPTIONAL: Mesh velocity field.
  ! DmeshVelocity(1,ivt) gives the X-velocity of the mesh, i.e. the X-velocity
  !   of the corner vertex ivt.
  ! DmeshVelocity(2,ivt) gives the Y-velocity of the mesh, i.e. the Y-velocity
  !   of the corner vertex ivt.
  ! The parameter must be present if ALE is activated in the
  ! configuration parameter block by bALE=true.
  real(DP), dimension(:,:), intent(in), optional :: DmeshVelocity

  ! OPTIONAL: A scalar cubature information structure that specifies the cubature
  ! formula(s) to use. If not specified, default settings are used.
  type(t_scalarCubatureInfo), intent(in), optional, target :: rcubatureInfo

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! System matrix.
  ! The content of the matrix must be present if cdef=CONV_MODMATRIX or
  ! =CONV_MODBOTH, otherwise only the structure is used.
  ! The nonlinear operator is added to the matrix.
  type(t_matrixScalar), intent(inout) :: rmatrix

  ! OPTIONAL: Defect vector.
  ! Must have the same structure as rsolution/rvecPrimary/rvecSecondary.
  ! Must be present if cdef=CONV_MODDEFECT or =CONV_MODBOTH.
  ! The nonlinear part is subtracted from this vector:
  ! <tex>$ r = r - \theta * u_1*grad(u_2) $</tex>
  type(t_vectorBlock), intent(inout), optional, target :: rdefect
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: icomponent
    integer(I32) :: celement
    real(DP), dimension(:), pointer :: p_DvelX1,p_DvelX2,p_DvelY1,p_DvelY2
    real(DP), dimension(:), pointer :: p_DsolX,p_DsolY,p_DdefectX,p_DdefectY

    ! At first check the input parameters that everything is present what
    ! we need:
    if ((cdef .eq. CONV_MODDEFECT) .or. (cdef .eq. CONV_MODBOTH)) then
      if ((.not. present(rsolution)) .or. (.not. present(rdefect))) then
        call output_line ("Solution/defect vector not present!", &
            OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusion2d")
        call sys_halt()
      end if
    end if

    if (rconfig%bALE) then
      if (.not. present(DmeshVelocity)) then
        call output_line ("Mesh velocity vector not present!", &
            OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusion2d")
        call sys_halt()
      end if
    end if


    ! At the moment, we only support a rather limited set of configurations:
    ! Matrix and vectors must all be double precision, matrix must be format
    ! 7 or 9, discretisation must be Q1~, constant viscosity.
    if ((rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
        (rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
      call output_line ("Unsupported matrix format!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusion2d")
      call sys_halt()
    end if

    celement = rmatrix%p_rspatialDiscrTest%RelementDistr(1)%celement
    if (rmatrix%p_rspatialDiscrTest%ccomplexity .ne. SPDISC_UNIFORM) then
      call output_line ("Unsupported discretisation!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusion2d")
      call sys_halt()
    end if

    if ((rvecPrimary%cdataType .ne. ST_DOUBLE) .or. &
        (rvecSecondary%cdataType .ne. ST_DOUBLE)) then
      call output_line ("Unsupported vector data type in velocity!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusion2d")
      call sys_halt()
    end if

    if (present(rdefect)) then
      if ((rsolution%cdataType .ne. ST_DOUBLE) .or. &
          (rdefect%cdataType .ne. ST_DOUBLE)) then
        call output_line ("Unsupported vector data type in solution/defect!", &
            OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusion2d")
        call sys_halt()
      end if
    end if

    if (.not. rconfig%bconstViscosity) then
      call output_line ("Only constant viscosity supported at the moment!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusion2d")
      call sys_halt()
    end if

    if (rconfig%dnu .eq. SYS_INFINITY_DP) then
      call output_line ("Viscosity parameter nu not initialised!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusion2d")
      call sys_halt()
    end if

    ! Hide the p_rsol...-parameters to prevent passing the NULL()-pointer
    ! if rsolution is not present -- some compilers do not like that ^^

    call lsyssc_getbase_double (rvecPrimary%RvectorBlock(1),p_DvelX1)
    call lsyssc_getbase_double (rvecPrimary%RvectorBlock(2),p_DvelY1)
    call lsyssc_getbase_double (rvecSecondary%RvectorBlock(1),p_DvelX2)
    call lsyssc_getbase_double (rvecSecondary%RvectorBlock(2),p_DvelY2)

    !!! DEBUG:
    !WHERE (abs(p_DvelX1) .LT. 1E-12_DP) p_DvelX1 = 0.0_DP
    !WHERE (abs(p_DvelY1) .LT. 1E-12_DP) p_DvelY1 = 0.0_DP
    !call vecio_writeArray_Dble (p_DvelX1, 'vecx1', &
    !                               0, 'vectorx1.txt', '(D10.3)')
    !call vecio_writeArray_Dble (p_DvelY1, 'vecx2', &
    !                               0, 'vectorx2.txt', '(D10.3)')

    if (present(rdefect)) then

      if (rdefect%nblocks .eq. 2) then

        ! Special 2D variant.
        call lsyssc_getbase_double (rsolution%RvectorBlock(1),p_DsolX)
        call lsyssc_getbase_double (rsolution%RvectorBlock(2),p_DsolY)
        call lsyssc_getbase_double (rdefect%RvectorBlock(1),p_DdefectX)
        call lsyssc_getbase_double (rdefect%RvectorBlock(2),p_DdefectY)

        call conv_strdiff2dALE_double ( &
                      p_DvelX1,p_DvelY1,p_DvelX2,p_DvelY2,dprimWeight,dsecWeight, &
                      rmatrix,cdef, rconfig%dupsam, rconfig%dnu, &
                      rconfig%dalpha, rconfig%dbeta, rconfig%dtheta, rconfig%ddelta, &
                      rconfig%bALE, rconfig%clocalh,&
                      p_DsolX,p_DsolY,p_DdefectX,p_DdefectY, DmeshVelocity,&
                      rcubatureInfo,rperfconfig)

      else

        ! Apply the operator to all blocks in the vector, block by block.
        !
        ! This sets up:
        !
        !   b1      b1        A           x1
        !   b2   =  b2    -      A        x2
        !   ...     ...             ...   ...
        !
        do icomponent = 1,rdefect%nblocks

          call lsyssc_getbase_double (rsolution%RvectorBlock(icomponent),p_DsolX)
          call lsyssc_getbase_double (rdefect%RvectorBlock(icomponent),p_DdefectX)

          if (rmatrix%NEQ .ne. rdefect%RvectorBlock(icomponent)%NEQ) then
            call output_line ("Velocity matrix not compatible to component "//&
                trim(sys_siL(icomponent,10))//" of the defect/solution vector!", &
                OU_CLASS_ERROR,OU_MODE_STD,'conv_streamlineDiffusion2d')
            call sys_halt()
          end if

          call conv_strdiff2dALEsingle_double ( &
                        p_DvelX1,p_DvelY1,p_DvelX2,p_DvelY2,dprimWeight,dsecWeight, &
                        rmatrix,cdef, rconfig%dupsam, rconfig%dnu, &
                        rconfig%dalpha, rconfig%dbeta, rconfig%dtheta, rconfig%ddelta, &
                        rconfig%bALE, rconfig%clocalh,&
                        p_DsolX,p_DdefectX, DmeshVelocity,rcubatureInfo,rperfconfig)
        end do

      end if

    else

      ! Calculate only the matrix. "The" matrix is used by the caller
      ! on all diagonal blocks of a block matrix!
      call conv_strdiff2dALEsingle_double ( &
                    p_DvelX1,p_DvelY1,p_DvelX2,p_DvelY2,dprimWeight,dsecWeight, &
                    rmatrix, cdef, rconfig%dupsam, rconfig%dnu, &
                    rconfig%dalpha, rconfig%dbeta, rconfig%dtheta, rconfig%ddelta, &
                    rconfig%bALE, rconfig%clocalh,DmeshVelocity=DmeshVelocity,&
                    rcubatureInfo=rcubatureInfo,rperfconfig=rperfconfig)

      !!! DEBUG:
      !call matio_writeMatrixHR (rmatrix, 'matrix',&
      !                          .TRUE., 0, 'matrixL.txt', '(D10.3)')

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine conv_strdiff2dALEsingle_double ( &
                  u1Xvel,u1Yvel,u2Xvel,u2Yvel,dweight1,dweight2,&
                  rmatrix,cdef, &
                  dupsam,dnu,dalpha,dbeta,dtheta, ddelta, bALE, &
                  clocalh, Du1,Ddef1, DmeshVelocity, rcubatureInfo, rperfconfig)
!<description>
  ! Standard streamline diffusion method to set up the operator
  ! <tex>
  !   $$ dtheta  *  ( dalpha * MASS  +  dbeta * STOKES  +  ddelta * u_1 * grad(u_2) ) $$
  ! </tex>
  ! in a matrix or to build a defect vector with that.
  ! 2D-version (X- and Y-velocity), uniform <tex>$\tilde Q_1$</tex> discretisation,
  ! double precision vectors/matrix.
  !
  ! u1Xvel,u1Yvel, u2Xvel,u2Yvel are two velocity field vectors,
  ! (u1Xvel,u1Yvel) a primary and (u2Xvel,u2Yvel) a secondary velocity field.
  ! The final velocity vector field is then computed as a weighted average
  ! of these two:
  ! <tex>
  !  $$ u_1  =  dweight1 * u1vel  +  dweight2 * u2vel $$
  ! </tex>
  ! <tex>$ u_2 = rsolution(.) $</tex> defines the second velocity field inside of
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
  ! [Duarte, Formaz, Natesan; `Arbitrary Lagrangian-Euler Method
  ! for Navier-Stokes equations with moving boundaries`;
  ! Comput. Methods Appl. Mech. Engrg. 193 (2004), 4819-4836]
  !
  ! Remarks:\\
  !
  ! 1.) In a typical call of the upwinding, the caller can use:
  !     dweight1 = 1, u1Xvel/u1Yvel = velocity field
  !     dweight2 = 0, u2Xvel/u2Yvel = undefined
  !   So the upwinding scheme only uses one velocity field.
  !   Such a call e.g. adds the integral
  !                <tex> $$ ( u_1 * grad(.) , v )_{\Omega} $$ </tex>
  !   to the system matrix.\\
  !
  !  2.) In case that there are two velocity fields representing
  !   the solution (may happen in a nonstationary simulation where
  !   u1Xvel/u1Yvel represents the solution in the current and u2Xvel/u2Yvel
  !   that of the previous time step), dweight1/dweight2 defines how these both
  !   velocity vectors should be weighted to compute the actual
  !   velocity field for the assembling:
  !               <tex> $$ U_act = dweight1*u1vel + dweight2*u2vel $$ </tex>
  !   This is e.g. used for the linear extrapolation technique to
  !   reconstruct a velocity from two previous time steps...\\
  !
  !  3.) In the nonlinear iteration, as a right hand side there arises
  !   a defect vector D whose linear part can easily be assembled.
  !   However, there is a nonlinearity to be included into that vector,
  !   too. By setting cdef=1,2, this routine incorporates the nonlinearity
  !   into that vector, using the formula
  !
  !            <tex> $$ D = D - dtheta * UUx * grad (Ux) $$ </tex>
  !
  !  4.) If bALE=true, a mesh velocity field is added to the nonlinearity
  !   according to the formula  "U * grad (U-DmeshVelocity)".
  !   For bALE=false, the simple nonlinearity "U * grad (U)" is used.
  !
  ! The target in this routine is either a single matrix block (for setting up
  ! a matrix) or a single block in a block matrix. A variant of this
  ! which acts simultaneously on two blocks of a velocity field can be found
  ! below!

!</description>

!<input>

  ! Primary X-velocity of <tex>$ u_1 $</tex>
  real(DP), dimension(:), intent(in) :: u1Xvel

  ! Primary Y-velocity of <tex>$ u_1 $</tex>
  real(DP), dimension(:), intent(in) :: u1Yvel

  ! Secondary X-velocity of <tex>$ u_1 $</tex>
  real(DP), dimension(:), intent(in) :: u2Xvel

  ! Secondary Y-velocity of <tex>$ u_1 $</tex>
  real(DP), dimension(:), intent(in) :: u2Yvel

  ! Computation/defect correction method. One of the CONV_MODxxxx constants:
  ! CONV_MODMATRIX: Set up the nonlinear matrix. rmatrix must be present, the
  !                 nonlinear part is added to the matrix.
  ! CONV_MODDEFECT: Set up the nonlinear defect. rdefect and rsolution must be
  !                 present.
  ! CONV_MODBOTH  : Set up the nonlinear matrix as well as the nonlinear defect.
  !                 rmatrix, rdefect and rsolution must all be present.
  integer, intent(in) :: cdef

  ! Weighting factor for u1Xvel/u1Yvel.
  real(DP), intent(in) :: dweight1

  ! Weighting factor for u2Xvel/u2Yvel.
  real(DP), intent(in) :: dweight2

  ! dupsam  - control parameter.
  !          -1: simple upwind,
  !          =0: Samarskji upwind
  real(DP), intent(in) :: dupsam

  ! Viscosity parameter <tex>$ \nu = 1/Re $</tex> if viscosity is constant
  real(DP), intent(in) :: dnu

  ! Weighting factor for the mass matrix.
  real(DP), intent(in) :: dalpha

  ! Weighting factor for the Stokes matrix. (Stokes matrix = 1/Re * Laplace)
  real(DP), intent(in) :: dbeta

  ! Weighting factor of the convective operator: <tex>$ \theta * u*grad(u) $</tex>.
  ! For time-dependent problems, this can be set to the step size
  ! in the <tex>$ \Theta $</tex>-scheme.
  real(DP), intent(in) :: dtheta

  ! Weighting factor for the nonlinear term
  real(DP), intent(in) :: ddelta

  ! Whether or not to use the ALE method
  logical, intent(in) :: bALE

  ! Method how to compute the local h
  integer, intent(in) :: clocalh

  ! OPTIONAL: Mesh velocity field. Must be present if bALE=TRUE.
  ! DmeshVelocity(1,:) gives the X-velocity of all the corner points of the mesh,
  ! DmeshVelocity(2,:) gives the Y-velocity.
  real(DP), dimension(:,:), intent(in), optional :: DmeshVelocity(:,:)

  ! OPTIONAL: velocity vector <tex>$ u_2 $</tex>. Must be present if cdef=CONV_MODDEFECT
  ! or cdef=CONV_MODBOTH.
  real(DP), dimension(:), intent(in), optional :: Du1

  ! OPTIONAL: A scalar cubature information structure that specifies the cubature
  ! formula(s) to use. If not specified, default settings are used.
  type(t_scalarCubatureInfo), intent(in), optional, target :: rcubatureInfo

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! The system matrix. Must be format 7 or 9.
  type(t_matrixScalar), intent(inout), target :: rmatrix

  ! OPTIONAL: defect vector. Must be present if cdef=CONV_MODDEFECT
  ! or =CONV_MODBOTH.
  real(DP), dimension(:), intent(inout), optional :: Ddef1
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: indof,indofALE,IEQ,IDOFE,JDOFE,icubp
  integer :: JCOL0,IDFG,JDFG,JCOL
  integer :: IEL,IELset,IELmax
  logical, dimension(EL_MAXNDER) :: Bder,BderALE
  real(DP) :: dumax,dumaxr, du1loc, du2loc, dunorm,db,OM,AH,denth,dre,dny
  real(DP) :: HBASI1,HBASI2,HBASI3,HBASJ1,HBASJ2,HBASJ3,HSUMI,HSUMJ
  integer :: NVE

  ! Matrix structure arrays
  integer, dimension(:), pointer :: p_Kcol
  integer, dimension(:), pointer :: p_Kld
  real(DP), dimension(:), pointer :: p_Da

  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  real(DP), dimension(:,:), allocatable :: p_DcubPtsRef

  ! The discretisation - for easier access
  type(t_spatialDiscretisation), pointer :: p_rdiscretisation

  ! Triangulation
  type(t_triangulation), pointer :: p_rtriangulation
  real(DP), dimension(:,:), pointer :: p_DvertexCoords
  integer, dimension(:,:), pointer :: p_IedgesAtElement,p_IverticesAtElement

  ! Number of elements in a block. Normally =NELEMSIM,
  ! except if there are less elements in the discretisation.
  integer :: nelementsPerBlock

  ! One and only element distribution
  type(t_elementDistribution), pointer :: p_relementDistribution

  ! For every cubature point on the reference element,
  ! the corresponding cubature weight
  real(DP), dimension(:), allocatable :: Domega

  ! number of cubature points on the reference element
  integer :: ncubp

  ! An element evaluation set for evaluating elements.
  type(t_evalElementSet) :: revalElementSet

  ! Arrays for saving Jacobian determinants
  real(DP), dimension(:,:), pointer :: p_Ddetj

  ! An allocateable array accepting the DOF`s of a set of elements.
  integer, dimension(:,:), allocatable, target :: Idofs, IdofsALE

  ! Allocateable arrays for the values of the basis functions -
  ! for test and trial spaces.
  real(DP), dimension(:,:,:,:), allocatable, target :: Dbas,DbasALE

  ! Local matrices, used during the assembly.
  ! Values and positions of values in the global matrix.
  integer, dimension(:,:,:), allocatable :: Kentry
  real(DP), dimension(:,:), allocatable :: Dentry

  ! A pointer to an element-number list
  integer, dimension(:), pointer :: p_IelementList

  ! Pointer to the velocity field in the cubature points.
  real(DP), dimension(:,:,:), allocatable :: Dvelocity

  ! An array with local DELTA`s, each DELTA for one element
  real(DP), dimension(:), allocatable :: DlocalDelta

  ! Type of transformation from the reference to the real element
  integer(I32) :: ctrafoType

  ! Element evaluation tag; collects some information necessary for evaluating
  ! the elements.
  integer(I32) :: cevaluationTag

  ! Cubature formula
  integer(I32) :: ccubature
  type(t_scalarCubatureInfo), target :: rtempCubatureInfo
  type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo

  ! Pointer to the performance configuration
  type(t_perfconfig), pointer :: p_rperfconfig

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => conv_perfconfig
    end if

    ! Initialise the derivative flags
    Bder = .false.
    Bder(DER_FUNC) = .true.
    Bder(DER_DERIV_X) = .true.
    Bder(DER_DERIV_Y) = .true.

    ! For ALE we do not even need so much
    BderALE = .false.
    BderALE(DER_FUNC) = .true.

    ! Shortcut to the spatial discretisation
    p_rdiscretisation => rmatrix%p_rspatialDiscrTest

    ! Get the element distribution. Here, we can find information about
    ! the cubature formula etc...
    p_relementDistribution => p_rdiscretisation%RelementDistr(1)

    ! Get some information about the triangulation
    p_rtriangulation => p_rdiscretisation%p_rtriangulation
    call storage_getbase_double2d (p_rtriangulation%h_DvertexCoords,&
                                   p_DvertexCoords)
    call storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
                                p_IverticesAtElement)
    call storage_getbase_int2d (p_rtriangulation%h_IedgesAtElement,&
                                p_IedgesAtElement)

    ! Get the number of local DOF`s for trial/test functions.
    ! We assume trial and test functions to be the same.
    indof = elem_igetNDofLoc(p_relementDistribution%celement)

    ! Get the number of local DOF`s Q1 -- we need them for ALE.
    indofALE = elem_igetNDofLoc(p_relementDistribution%celement)

    ! Number of local DOF`s
    NVE = elem_igetNVE(p_relementDistribution%celement)

    ! Do we have an assembly structure?
    ! If we do not have it, create a cubature info structure that
    ! defines how to do the assembly.
    if (.not. present(rcubatureInfo)) then
      call spdiscr_createDefCubStructure(p_rdiscretisation,&
          rtempCubatureInfo,CUB_GEN_DEPR_EVAL)
      p_rcubatureInfo => rtempCubatureInfo
    else
      p_rcubatureInfo => rcubatureInfo
    end if

    ! Cubature formula. Only one cubature formula supported.
    ccubature = p_rcubatureInfo%p_RinfoBlocks(1)%ccubature

    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it is at most
    ! NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = min(p_rperfconfig%NELEMSIM,p_rtriangulation%NEL)

    ! For cdef containing CONV_MODDEFECT, we build the defect vector
    !     D = RHS - A*U
    ! In this case, the defect(rhs vectors must be present

    if (iand(cdef,CONV_MODDEFECT) .ne. 0) then
      if (.not. (present(Ddef1) .and. present(Du1) )) then
        call output_line ("Necessary arguments missing!", &
            OU_CLASS_ERROR,OU_MODE_STD,"conv_strdiff2dALEsingle_double")
        call sys_halt()
      end if
    end if

    if (iand(cdef,CONV_MODMATRIX) .ne. 0) then
      ! Get matrix arrays
      call lsyssc_getbase_double (rmatrix,p_Da)
    end if
    call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
    call lsyssc_getbase_Kld (rmatrix,p_Kld)

    ! Get from the trial element space the type of coordinate system
    ! that is used there:
    ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)

    ! Get the number of cubature points for the cubature formula
    ncubp = cub_igetNumPts(ccubature)

    ! Allocate two arrays for the points and the weights
    allocate(Domega(ncubp))
    allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),ncubp))

    ! Get the cubature formula
    call cub_getCubature(ccubature,p_DcubPtsRef, Domega)

    ! Allocate an array saving the coordinates of corner vertices of elements

    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  allocate(Dbas(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    ! would lead to nonused memory blocks in these arrays during the assembly,
    ! which reduces the speed by 50%!
    allocate(Dbas(indof,elem_getMaxDerivative(p_relementDistribution%celement), &
             ncubp,nelementsPerBlock))

    ! Allocate memory for the DOF`s of all the elements.
    allocate(Idofs(indof,nelementsPerBlock))

    ! The same for the ALE-space
    allocate(DbasALE(indofALE,elem_getMaxDerivative(EL_Q1), &
             ncubp,nelementsPerBlock))

    ! Allocate memory for the DOF`s of all the elements.
    allocate(IdofsALE(indofALE,nelementsPerBlock))

    ! Allocate memory for array with local DELTA`s
    allocate(DlocalDelta(nelementsPerBlock))

    ! Allocate an array saving the local matrices for all elements
    ! in an element set.
    ! We could also allocate EL_MAXNBAS*EL_MAXNBAS*NELEMSIM integers
    ! for this local matrix, but this would normally not fit to the cache
    ! anymore! indofTrial*indofTest*NELEMSIM is normally much smaller!
    allocate(Kentry(indof,indof,nelementsPerBlock))
    allocate(Dentry(indof,indof))

    ! Allocate memory for the velocity in the cubature points.
    allocate(Dvelocity(NDIM2D,ncubp,nelementsPerBlock))

    ! Initialisation of the element set.
    call elprep_init(revalElementSet)

    ! What is the reciprocal of nu? We need it later.
    if (dnu .ne. 0.0_DP) then
      dre = 1.0_DP/dnu

      ! dny gets the actual multiplier for the Laplace matrix.
      ! Remember: dbeta*Stokes = dbeta*dnu*Laplace = dny*Laplace.
      ! This may be =0.0 if the Stokes operator should not be included into
      ! the matrix.
      dny = dbeta*dnu
    else
      call output_line ("NU=0 not allowed! Set dbeta=0 to prevent Stokes operator "//&
          "from being build!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_strdiff2dALEsingle_double")
      call sys_halt()
    end if

    ! If ddelta=0, we have to neglect the nonlinearity. In both cases,
    ! set DlocalDelta=0 which disables the nonlinear term in the assembly.
    ! If dupsam=0, we neglect the stabilisation term (central difference like
    ! discretisation), so we set DlocalDelta=0 as well.
    if ((ddelta .eq. 0.0_DP) .or. (dupsam .eq. 0.0_DP)) then
      call lalg_clearVectorDble (DlocalDelta)
    end if

    ! Calculate the maximum norm of the actual velocity field
    ! U = A1*U1 + A2*U2 into DUMAX.
    ! Round up the norm to 1D-8 if it is too small...

    dumax=0.0_DP
    if (dweight2 .eq. 0.0_DP) then
      do IEQ=1,size(u1Xvel)
        du1loc = dweight1*u1Xvel(IEQ)
        du2loc = dweight1*u1Yvel(IEQ)
        dunorm = sqrt(du1loc**2+du2loc**2)
        dumax = max(DUMAX,DUNORM)
      end do
    else
      do ieq=1,size(u1Xvel)
        du1loc = dweight1*u1Xvel(IEQ)+dweight2*u2Xvel(IEQ)
        du2loc = dweight1*u1Yvel(IEQ)+dweight2*u2Yvel(IEQ)
        dunorm = sqrt(du1loc**2+du2loc**2)
        dumax = max(dumax,dunorm)
      end do
    end if

    if (dumax.lt.1E-8_DP) dumax=1E-8_DP
    dumaxr = 1.0_DP/dumax

    ! p_IelementList must point to our set of elements in the discretisation
    ! with that combination of trial/test functions
    call storage_getbase_int (p_relementDistribution%h_IelementList, &
                              p_IelementList)

    ! Loop over the elements - blockwise.
    do IELset = 1, size(p_IelementList), nelementsPerBlock

      ! We always handle NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most NELEMSIM
      ! elements simultaneously.

      IELmax = min(size(p_IelementList),IELset-1+p_rperfconfig%NELEMSIM)

      ! The outstanding feature with finite elements is: A basis
      ! function for a DOF on one element has common support only
      ! with the DOF`s on the same element! E.g. for Q1:
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
      !     to collect all "O" DOF`s.
      !
      ! Calculate the global DOF`s into IdofsTrial / IdofsTest.
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF`s of our NELEMSIM elements simultaneously.
      call dof_locGlobMapping_mult(p_rdiscretisation, p_IelementList(IELset:IELmax), &
                                  Idofs)

      ! In case ALE is used, do this also for the ALE stuff.
      if (bALE) then
        call dof_locGlobMapping_mult(p_rdiscretisation, &
                                    p_IelementList(IELset:IELmax), &
                                    IdofsALE)
      end if

      ! Calculate local DELTA`s for streamline diffusion method.
      ! (cf. p. 121 in Turek`s CFD book).
      ! For every element, we need a local DELTA.
      ! Every local delta is weighted by the global "ddelta".
      ! If ddelta=0, we do not do anything as this disables the
      ! nonlinear term.
      ! If UPSAM=0.0, we have a central-difference like discretisation, which
      ! is one can see as the local stabilisation weight Delta is also = 0.0.
      ! In this case, we even switch of the calculation of the local Delta,
      ! as it is always =0.0, so we save a little bit time.
      if ((ddelta .ne. 0.0_DP) .and. (dupsam .ne. 0.0_DP))then
        call getLocalDeltaQuad (clocalh,&
                      u1Xvel,u1Yvel,u2Xvel,u2Yvel,dweight1,dweight2,&
                      p_IelementList(IELset:IELmax),&
                      duMaxR,DlocalDelta,p_rtriangulation,Idofs,dupsam,dre)
      end if

      ! For the assembly of the global matrix, we use a "local"
      ! approach. At first we build a "local" system matrix according
      ! to the current element. This contains all additive
      ! contributions of element IEL, which are later added at the
      ! right positions to the elements in the global system matrix.
      !
      ! We have indofTrial trial DOF`s per element and
      ! indofTest test DOF`s per element. Therefore there are
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
      do IEL=1,IELmax-IELset+1

        ! For building the local matrices, we have first to
        ! loop through the test functions (the "O"`s), as these
        ! define the rows in the matrix.
        do IDOFE=1,indof

          ! Row IDOFE of the local matrix corresponds
          ! to row=global DOF KDFG(IDOFE) in the global matrix.
          ! This is one of the the "O"`s in the above picture.
          ! Get the starting position of the corresponding row
          ! to JCOL0:

          JCOL0=p_KLD(Idofs(IDOFE,IEL))

          ! Now we loop through the other DOF`s on the current element
          ! (the "O"`s).
          ! All these have common support with our current basis function
          ! and will therefore give an additive value to the global
          ! matrix.

          do JDOFE=1,indof

            ! Get the global DOF of the "X" which interacts with
            ! our "O".

            JDFG=Idofs(JDOFE,IEL)

            ! Starting in JCOL0 (which points to the beginning of
            ! the line initially), loop through the elements in
            ! the row to find the position of column IDFG.
            ! Jump out of the do loop if we find the column.

            do JCOL=JCOL0,rmatrix%NA
              if (p_KCOL(JCOL) .eq. JDFG) exit
            end do

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

          end do ! IDOFE

        end do ! JDOFE

      end do ! IEL

      ! Ok, we found the positions of the local matrix entries
      ! that we have to change.
      ! To calculate the matrix contributions, we have to evaluate
      ! the elements to give us the values of the basis functions
      ! in all the DOF`s in all the elements in our set.
      !
      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag.
      cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)

      ! In the first loop, calculate the coordinates on the reference element.
      ! In all later loops, use the precalculated information.
      if (IELset .eq. 1) then
        cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REFPOINTS)
      else
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
      end if

      ! Calculate all information that is necessary to evaluate the finite element
      ! on all cells of our subset. This includes the coordinates of the points
      ! on the cells.
      call elprep_prepareSetForEvaluation (revalElementSet,&
          cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
          ctrafoType, p_DcubPtsRef(:,1:ncubp), rperfconfig=rperfconfig)
      p_Ddetj => revalElementSet%p_Ddetj

      ! Calculate the values of the basis functions.
      call elem_generic_sim2 (p_relementDistribution%celement, &
          revalElementSet, Bder, Dbas)

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
      ! DOF`s as:
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

      if (dweight2 .eq. 0.0_DP) then

        ! Loop over all elements in the current set
        do IEL=1,IELmax-IELset+1

          ! Loop over all cubature points on the current element
          do ICUBP = 1, ncubp

            du1loc = 0.0_DP
            du2loc = 0.0_DP

            ! Perform a loop through the trial DOF`s.
            do JDOFE=1,indof

              ! Get the value of the (test) basis function
              ! phi_i (our "O") in the cubature point:
              db = Dbas(JDOFE,1,ICUBP,IEL)

              ! Sum up to the value in the cubature point
              JDFG = Idofs(JDOFE,IEL)
              du1loc = du1loc + u1Xvel(JDFG)*db
              du2loc = du2loc + u1Yvel(JDFG)*db

            end do ! JDOFE

            ! Save the computed velocity
            Dvelocity(1,ICUBP,IEL) = dweight1*du1loc
            Dvelocity(2,ICUBP,IEL) = dweight1*du2loc

          end do ! ICUBP

        end do ! IEL

      else

        do IEL=1,IELmax-IELset+1

          ! Loop over all cubature points on the current element
          do ICUBP = 1, ncubp

            du1loc = 0.0_DP
            du2loc = 0.0_DP

            ! Perform a loop through the trial DOF`s.
            do JDOFE=1,indof

              ! Get the value of the (trial) basis function
              ! phi_i in the cubature point:
              db = Dbas(JDOFE,1,ICUBP,IEL)

              ! Sum up to the value in the cubature point
              JDFG = Idofs(JDOFE,IEL)
              du1loc = du1loc + (dweight1*u1Xvel(JDFG) + dweight2*u2Xvel(JDFG))*db
              du2loc = du2loc + (dweight1*u1Yvel(JDFG) + dweight2*u2Yvel(JDFG))*db

            end do ! JDOFE

            ! Save the computed velocity
            Dvelocity(1,ICUBP,IEL) = du1loc
            Dvelocity(2,ICUBP,IEL) = du2loc

          end do ! ICUBP

        end do ! IEL

      end if

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

      if (bALE) then

        ! Calculate the values of the basis functions in all the points
        ! on all the elements
        call elem_generic_sim2 (EL_Q1, &
            revalElementSet, Bder, DbasALE)

        ! Loop over all elements in the current set
        do IEL=1,IELmax-IELset+1

          ! Loop over all cubature points on the current element
          do ICUBP = 1, ncubp

            du1loc = 0.0_DP
            du2loc = 0.0_DP

            ! Perform a loop through the trial DOF`s.
            do JDOFE=1,indof

              ! Get the value of the (trial) basis function
              db= Dbas(JDOFE,1,ICUBP,IEL)

              ! Sum up to the value in the cubature point
              JDFG = IdofsALE(IDOFE,IEL)
              du1loc = du1loc + DmeshVelocity(1,JDFG)*db
              du2loc = du2loc + DmeshVelocity(2,JDFG)*db

            end do ! JDOFE

            ! Save the computed velocity
            Dvelocity(1,ICUBP,IEL) = Dvelocity(1,ICUBP,IEL) - du1loc
            Dvelocity(2,ICUBP,IEL) = Dvelocity(2,ICUBP,IEL) - du2loc

          end do ! ICUBP

        end do ! IEL

      end if

      ! Ok, we now use Dvelocity as coefficient array in the assembly
      ! of a biliinear form!
      !
      ! Loop over the elements in the current set.

      do IEL=1,IELmax-IELset+1

        ! Clear the local matrix
        Dentry = 0.0_DP

        ! Loop over all cubature points on the current element
        do ICUBP = 1, ncubp

          ! Calculate the current weighting factor in the cubature formula
          ! in that cubature point.
          !
          ! Normally, we have to take the absolut value of the determinant
          ! of the mapping here!
          ! In 2D, the determinant is always positive, whereas in 3D,
          ! the determinant might be negative -- that is normal!
          ! But because this routine only works in 2D, we can skip
          ! the ABS here!

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
          ! Outer loop over the DOF`s i=1..indof on our current element,
          ! which corresponds to the basis functions Phi_i:

          do IDOFE=1,indof

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

            ! Inner loop over the DOF`s j=1..indof, which corresponds to
            ! the basis function Phi_j:

            do JDOFE=1,indof

              if (IDOFE.eq.JDOFE) then

                ! Short version of the evaluation of the matrix
                ! contribution - see below for a more detailed
                ! description what is added together here!

                AH = ddelta*HSUMI*(DlocalDelta(IEL)*HSUMI+HBASI1) &
                    + dny*(HBASI2**2+HBASI3**2) &
                    + dalpha*HBASI1**2

              else

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
                ! we do not have to worry about that.

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

                AH = ddelta*HSUMJ*(DlocalDelta(IEL)*HSUMI+HBASI1) &
                    + dny*(HBASI2*HBASJ2+HBASI3*HBASJ3) &
                    + dalpha*HBASI1*HBASJ1

              end if ! (IDOFE.EQ.JDOFE)

              ! Weighten the calculated value AH by the cubature
              ! weight OM and add it to the local matrix. After the
              ! loop over all DOF`s is finished, each entry contains
              ! the calculated integral.

              Dentry(JDOFE,IDOFE) = Dentry(JDOFE,IDOFE)+OM*AH

            end do ! IDOFE

          end do ! JDOFE

        end do ! ICUBP

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

        if (iand(cdef,CONV_MODMATRIX) .ne. 0) then
          do IDOFE=1,indof
            do JDOFE=1,indof
              p_DA(Kentry(JDOFE,IDOFE,IEL)) = p_DA(Kentry(JDOFE,IDOFE,IEL)) + &
                dtheta * Dentry(JDOFE,IDOFE)
            end do
          end do
        end if

        ! For cdef containing CONV_MODDEFECT, build the defect vector
        !     D = RHS - A*U
        ! This is done matrix free, only with the help of the local
        ! matrix.
        ! In this case, D=(D1,D2) is expected to be the RHS on
        ! entry and will be updated to be the defect vector when
        ! this routine is left.

        if (iand(cdef,CONV_MODDEFECT) .ne. 0) then
          do IDOFE=1,indof

            IDFG=Idofs(IDOFE,IEL)

            do JDOFE=1,indof

              denth = dtheta*Dentry(JDOFE,IDOFE)

              JDFG=Idofs(JDOFE,IEL)
              Ddef1(IDFG)= Ddef1(IDFG) - denth*Du1(JDFG)

            end do
          end do
        end if

      end do ! IEL

    end do ! IELset

    ! Release memory
    call elprep_releaseElementSet(revalElementSet)

    ! Release the assembly structure if necessary.
    if (.not. present(rcubatureInfo)) then
      call spdiscr_releaseCubStructure(rtempCubatureInfo)
    end if

    deallocate(p_DcubPtsRef)
    deallocate(Domega)
    deallocate(DlocalDelta)
    deallocate(Dvelocity)
    deallocate(Dentry)
    deallocate(Kentry)
    deallocate(IdofsALE)
    deallocate(Idofs)
    deallocate(DbasALE)
    deallocate(Dbas)

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine conv_strdiff2dALE_double ( &
                  u1Xvel,u1Yvel,u2Xvel,u2Yvel,dweight1,dweight2,&
                  rmatrix,cdef, &
                  dupsam,dnu,dalpha,dbeta,dtheta, ddelta, bALE, &
                  clocalh,Du1,Du2,Ddef1,Ddef2,DmeshVelocity,&
                  rcubatureInfo,rperfconfig)
!<description>
  ! Standard streamline diffusion method to set up the operator
  ! <tex>
  !   $$ dtheta  *  ( dalpha * MASS  +  dbeta * STOKES  +  ddelta * u_1 * grad(u_2) ) $$
  ! </tex>
  ! in a matrix or to build a defect vector with that.
  ! 2D-version (X- and Y-velocity), uniform <tex>$\tilde Q_1$</tex> discretisation,
  ! double precision vectors/matrix.
  !
  ! u1Xvel,u1Yvel, u2Xvel,u2Yvel are two velocity field vectors,
  ! (u1Xvel,u1Yvel) a primary and (u2Xvel,u2Yvel) a secondary velocity field.
  ! The final velocity vector field is then computed as a weighted average
  ! of these two:
  ! <tex>
  !   $$ u_1  =  dweight1 * u1vel  +  dweight2 * u2vel $$
  ! </tex>
  ! <tex>$ u_2 = rsolution(.) $</tex> defines the second velocity field inside of
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
  ! [Duarte, Formaz, Natesan; `Arbitrary Lagrangian-Euler Method
  ! for Navier-Stokes equations with moving boundaries`;
  ! Comput. Methods Appl. Mech. Engrg. 193 (2004), 4819-4836]
  !
  ! Remarks:\\
  !
  ! 1.) In a typical call of the upwinding, the caller can use:
  !     dweight1 = 1, u1Xvel/u1Yvel = velocity field
  !     dweight2 = 0, u2Xvel/u2Yvel = undefined
  !   So the upwinding scheme only uses one velocity field.
  !   Such a call e.g. adds the integral
  !                <tex> $$ ( u_1 * grad(.) , v )_{\Omega} $$ </tex>
  !   to the system matrix.\\
  !
  !  2.) In case that there are two velocity fields representing
  !   the solution (may happen in a nonstationary simulation where
  !   u1Xvel/u1Yvel represents the solution in the current and u2Xvel/u2Yvel
  !   that of the previous time step), dweight1/dweight2 defines how these both
  !   velocity vectors should be weighted to compute the actual
  !   velocity field for the assembling:
  !               <tex> $$ U_act = dweight1*u1vel + dweight2*u2vel $$ </tex>
  !   This is e.g. used for the linear extrapolation technique to
  !   reconstruct a velocity from two previous time steps...\\
  !
  !  3.) In the nonlinear iteration, as a right hand side there arises
  !   a defect vector D whose linear part can easily be assembled.
  !   However, there is a nonlinearity to be included into that vector,
  !   too. By setting cdef=1,2, this routine incorporates the nonlinearity
  !   into that vector, using the formula
  !
  !            <tex> $$ D = D - dtheta * UUx * grad (Ux) $$ </tex>
  !
  !  4.) If bALE=true, a mesh velocity field is added to the nonlinearity
  !   according to the formula  "U * grad (U-DmeshVelocity)".
  !   For bALE=false, the simple nonlinearity "U * grad (U)" is used.

!</description>

!<input>

  ! Primary X-velocity of <tex>$ u_1 $</tex>
  real(DP), dimension(:), intent(in) :: u1Xvel

  ! Primary Y-velocity of <tex>$ u_1 $</tex>
  real(DP), dimension(:), intent(in) :: u1Yvel

  ! Secondary X-velocity of <tex>$ u_1 $</tex>
  real(DP), dimension(:), intent(in) :: u2Xvel

  ! Secondary Y-velocity of <tex>$ u_1 $</tex>
  real(DP), dimension(:), intent(in) :: u2Yvel

  ! Computation/defect correction method. One of the CONV_MODxxxx constants:
  ! CONV_MODMATRIX: Set up the nonlinear matrix. rmatrix must be present, the
  !                 nonlinear part is added to the matrix.
  ! CONV_MODDEFECT: Set up the nonlinear defect. rdefect and rsolution must be
  !                 present.
  ! CONV_MODBOTH  : Set up the nonlinear matrix as well as the nonlinear defect.
  !                 rmatrix, rdefect and rsolution must all be present.
  integer, intent(in) :: cdef

  ! Weighting factor for u1Xvel/u1Yvel.
  real(DP), intent(in) :: dweight1

  ! Weighting factor for u2Xvel/u2Yvel.
  real(DP), intent(in) :: dweight2

  ! dupsam  - control parameter.
  !          -1: simple upwind,
  !          =0: Samarskji upwind
  real(DP), intent(in) :: dupsam

  ! Viscosity parameter <tex>$ \nu = 1/Re $</tex> if viscosity is constant
  real(DP), intent(in) :: dnu

  ! Weighting factor for the mass matrix.
  real(DP), intent(in) :: dalpha

  ! Weighting factor for the Stokes matrix. (Stokes matrix = 1/Re * Laplace)
  real(DP), intent(in) :: dbeta

  ! Weighting factor of the convective operator: <tex>$ \theta * u*grad(u) $</tex>.
  ! For time-dependent problems, this can be set to the step size
  ! in the <tex>$ \Theta $</tex>-scheme.
  real(DP), intent(in) :: dtheta

  ! Weighting factor for the nonlinear term
  real(DP), intent(in) :: ddelta

  ! Whether or not to use the ALE method
  logical, intent(in) :: bALE

  ! Method how to compute the local h
  integer, intent(in) :: clocalh

  ! OPTIONAL: Mesh velocity field. Must be present if bALE=TRUE.
  ! DmeshVelocity(1,:) gives the X-velocity of all the corner points of the mesh,
  ! DmeshVelocity(2,:) gives the Y-velocity.
  real(DP), dimension(:,:), intent(in), optional :: DmeshVelocity(:,:)

  ! OPTIONAL: X-velocity of <tex>$ u_2 $</tex>. Must be present if cdef=CONV_MODDEFECT
  ! or cdef=CONV_MODBOTH.
  real(DP), dimension(:), intent(in), optional :: Du1

  ! Y-velocity of <tex>$ u_2 $</tex>. Must be present if cdef=CONV_MODDEFECT
  ! or cdef=CONV_MODBOTH.
  real(DP), dimension(:), intent(in), optional :: Du2

  ! OPTIONAL: A scalar cubature information structure that specifies the cubature
  ! formula(s) to use. If not specified, default settings are used.
  type(t_scalarCubatureInfo), intent(in), optional, target :: rcubatureInfo

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! The system matrix. Must be format 7 or 9.
  type(t_matrixScalar), intent(inout), target :: rmatrix

  ! OPTIONAL: X-defect vector. Must be present if cdef=CONV_MODDEFECT
  ! or =CONV_MODBOTH.
  real(DP), dimension(:), intent(inout), optional :: Ddef1

  ! OPTIONAL: Y-defect vector. Must be present if cdef=CONV_MODDEFECT
  ! or =CONV_MODBOTH.
  real(DP), dimension(:), intent(inout), optional :: Ddef2
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: indof,indofALE,IEQ,IDOFE,JDOFE,icubp
  integer :: JCOL0,IDFG,JDFG,JCOL
  integer :: IEL,IELset,IELmax
  logical, dimension(EL_MAXNDER) :: Bder,BderALE
  real(DP) :: dumax,dumaxr, du1loc, du2loc, dunorm,db,OM,AH,denth,dre,dny
  real(DP) :: HBASI1,HBASI2,HBASI3,HBASJ1,HBASJ2,HBASJ3,HSUMI,HSUMJ
  integer :: NVE

  ! Matrix structure arrays
  integer, dimension(:), pointer :: p_Kcol
  integer, dimension(:), pointer :: p_Kld
  real(DP), dimension(:), pointer :: p_Da

  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  real(DP), dimension(:,:), allocatable :: p_DcubPtsRef

  ! The discretisation - for easier access
  type(t_spatialDiscretisation), pointer :: p_rdiscretisation

  ! Triangulation
  type(t_triangulation), pointer :: p_rtriangulation
  real(DP), dimension(:,:), pointer :: p_DvertexCoords
  integer, dimension(:,:), pointer :: p_IedgesAtElement,p_IverticesAtElement

  ! Number of elements in a block. Normally =NELEMSIM,
  ! except if there are less elements in the discretisation.
  integer :: nelementsPerBlock

  ! One and only element distribution
  type(t_elementDistribution), pointer :: p_relementDistribution

  ! For every cubature point on the reference element,
  ! the corresponding cubature weight
  real(DP), dimension(:), allocatable :: Domega

  ! number of cubature points on the reference element
  integer :: ncubp

  ! An element evaluation set for evaluating elements.
  type(t_evalElementSet) :: revalElementSet

  ! Arrays for saving Jacobian determinants
  real(DP), dimension(:,:), pointer :: p_Ddetj

  ! An allocateable array accepting the DOF`s of a set of elements.
  integer, dimension(:,:), allocatable, target :: Idofs, IdofsALE

  ! Allocateable arrays for the values of the basis functions -
  ! for test and trial spaces.
  real(DP), dimension(:,:,:,:), allocatable, target :: Dbas,DbasALE

  ! Local matrices, used during the assembly.
  ! Values and positions of values in the global matrix.
  integer, dimension(:,:,:), allocatable :: Kentry
  real(DP), dimension(:,:), allocatable :: Dentry

  ! A pointer to an element-number list
  integer, dimension(:), pointer :: p_IelementList

  ! Pointer to the velocity field in the cubature points.
  real(DP), dimension(:,:,:), allocatable :: Dvelocity

  ! An array with local DELTA`s, each DELTA for one element
  real(DP), dimension(:), allocatable :: DlocalDelta

  ! Type of transformation from the reference to the real element
  integer(I32) :: ctrafoType

  ! Element evaluation tag; collects some information necessary for evaluating
  ! the elements.
  integer(I32) :: cevaluationTag

  ! Pointer to the performance configuration
  type(t_perfconfig), pointer :: p_rperfconfig

  ! Cubature information structure
  type(t_scalarCubatureInfo), target :: rtempCubatureInfo
  type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => conv_perfconfig
    end if

    ! If we do not have it, create a cubature info structure that
    ! defines how to do the assembly.
    if (.not. present(rcubatureInfo)) then
      call spdiscr_createDefCubStructure(rmatrix%p_rspatialDiscrTrial,&
          rtempCubatureInfo,CUB_GEN_DEPR_BILFORM)
      p_rcubatureInfo => rtempCubatureInfo
    else
      p_rcubatureInfo => rcubatureInfo
    end if

    ! Initialise the derivative flags
    Bder = .false.
    Bder(DER_FUNC) = .true.
    Bder(DER_DERIV_X) = .true.
    Bder(DER_DERIV_Y) = .true.

    ! For ALE we do not even need so much
    BderALE = .false.
    BderALE(DER_FUNC) = .true.

    ! Shortcut to the spatial discretisation
    p_rdiscretisation => rmatrix%p_rspatialDiscrTest

    ! Get the element distribution. Here, we can find information about
    ! the cubature formula etc...
    p_relementDistribution => p_rdiscretisation%RelementDistr(1)

    ! Get some information about the triangulation
    p_rtriangulation => p_rdiscretisation%p_rtriangulation
    call storage_getbase_double2d (p_rtriangulation%h_DvertexCoords,&
                                   p_DvertexCoords)
    call storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
                                p_IverticesAtElement)
    call storage_getbase_int2d (p_rtriangulation%h_IedgesAtElement,&
                                p_IedgesAtElement)

    ! Get the number of local DOF`s for trial/test functions.
    ! We assume trial and test functions to be the same.
    indof = elem_igetNDofLoc(p_relementDistribution%celement)

    ! Get the number of local DOF`s Q1 -- we need them for ALE.
    indofALE = elem_igetNDofLoc(p_relementDistribution%celement)

    ! Number of local DOF`s
    NVE = elem_igetNVE(p_relementDistribution%celement)

    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it is at most
    ! NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = min(p_rperfconfig%NELEMSIM,p_rtriangulation%NEL)

    ! For cdef containing CONV_MODDEFECT, we build the defect vector
    !     D = RHS - A*U
    ! In this case, the defect(rhs vectors must be present

    if (iand(cdef,CONV_MODDEFECT) .ne. 0) then
      if (.not. (present(Ddef1) .and. present(Ddef2) .and. &
                 present(Du1) .and. present(Du2))) then
        call output_line ("Necessary arguments missing!", &
            OU_CLASS_ERROR,OU_MODE_STD,"conv_strdiff2dALE_double")
        call sys_halt()
      end if
    end if

    if (iand(cdef,CONV_MODMATRIX) .ne. 0) then
      ! Get matrix arrays
      call lsyssc_getbase_double (rmatrix,p_Da)
    end if
    call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
    call lsyssc_getbase_Kld (rmatrix,p_Kld)

    ! Get from the trial element space the type of coordinate system
    ! that is used there:
    ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)

    ! Get the number of cubature points for the cubature formula
    ncubp = cub_igetNumPts(p_rcubatureInfo%p_RinfoBlocks(1)%ccubature)

    ! Allocate two arrays for the points and the weights
    allocate(Domega(ncubp))
    allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),ncubp))

    ! Get the cubature formula
    call cub_getCubature(p_rcubatureInfo%p_RinfoBlocks(1)%ccubature,p_DcubPtsRef, Domega)

    ! Allocate an array saving the coordinates of corner vertices of elements

    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  allocate(Dbas(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    ! would lead to nonused memory blocks in these arrays during the assembly,
    ! which reduces the speed by 50%!
    allocate(Dbas(indof,elem_getMaxDerivative(p_relementDistribution%celement), &
             ncubp,nelementsPerBlock))

    ! Allocate memory for the DOF`s of all the elements.
    allocate(Idofs(indof,nelementsPerBlock))

    ! The same for the ALE-space
    allocate(DbasALE(indofALE,elem_getMaxDerivative(EL_Q1), &
             ncubp,nelementsPerBlock))

    ! Allocate memory for the DOF`s of all the elements.
    allocate(IdofsALE(indofALE,nelementsPerBlock))

    ! Allocate memory for array with local DELTA`s
    allocate(DlocalDelta(nelementsPerBlock))

    ! Allocate an array saving the local matrices for all elements
    ! in an element set.
    ! We could also allocate EL_MAXNBAS*EL_MAXNBAS*NELEMSIM integers
    ! for this local matrix, but this would normally not fit to the cache
    ! anymore! indofTrial*indofTest*NELEMSIM is normally much smaller!
    allocate(Kentry(indof,indof,nelementsPerBlock))
    allocate(Dentry(indof,indof))

    ! Allocate memory for the velocity in the cubature points.
    allocate(Dvelocity(NDIM2D,ncubp,nelementsPerBlock))

    ! Initialisation of the element set.
    call elprep_init(revalElementSet)

    ! What is the reciprocal of nu? We need it later.
    if (dnu .ne. 0.0_DP) then
      dre = 1.0_DP/dnu

      ! dny gets the actual multiplier for the Laplace matrix.
      ! Remember: dbeta*Stokes = dbeta*dnu*Laplace = dny*Laplace.
      ! This may be =0.0 if the Stokes operator should not be included into
      ! the matrix.
      dny = dbeta*dnu
    else
      call output_line ("NU=0 not allowed! Set dbeta=0 to prevent Stokes operator "//&
          "from being build!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_strdiff2dALE_double")
      call sys_halt()
    end if

    ! If ddelta=0, we have to neglect the nonlinearity. In both cases,
    ! set DlocalDelta=0 which disables the nonlinear term in the assembly.
    ! If dupsam=0, we neglect the stabilisation term (central difference like
    ! discretisation), so we set DlocalDelta=0 as well.
    if ((ddelta .eq. 0.0_DP) .or. (dupsam .eq. 0.0_DP)) then
      call lalg_clearVectorDble (DlocalDelta)
    end if

    ! Calculate the maximum norm of the actual velocity field
    ! U = A1*U1 + A2*U2 into DUMAX.
    ! Round up the norm to 1D-8 if it is too small...

    dumax=0.0_DP
    if (dweight2 .eq. 0.0_DP) then
      do IEQ=1,size(u1Xvel)
        du1loc = dweight1*u1Xvel(IEQ)
        du2loc = dweight1*u1Yvel(IEQ)
        dunorm = sqrt(du1loc**2+du2loc**2)
        dumax = max(DUMAX,DUNORM)
      end do
    else
      do ieq=1,size(u1Xvel)
        du1loc = dweight1*u1Xvel(IEQ)+dweight2*u2Xvel(IEQ)
        du2loc = dweight1*u1Yvel(IEQ)+dweight2*u2Yvel(IEQ)
        dunorm = sqrt(du1loc**2+du2loc**2)
        dumax = max(dumax,dunorm)
      end do
    end if

    if (dumax.lt.1E-8_DP) dumax=1E-8_DP
    dumaxr = 1.0_DP/dumax

    ! p_IelementList must point to our set of elements in the discretisation
    ! with that combination of trial/test functions
    call storage_getbase_int (p_relementDistribution%h_IelementList, &
                              p_IelementList)

    ! Loop over the elements - blockwise.
    do IELset = 1, size(p_IelementList), nelementsPerBlock

      ! We always handle NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most NELEMSIM
      ! elements simultaneously.

      IELmax = min(size(p_IelementList),IELset-1+p_rperfconfig%NELEMSIM)

      ! The outstanding feature with finite elements is: A basis
      ! function for a DOF on one element has common support only
      ! with the DOF`s on the same element! E.g. for Q1:
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
      !     to collect all "O" DOF`s.
      !
      ! Calculate the global DOF`s into IdofsTrial / IdofsTest.
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF`s of our NELEMSIM elements simultaneously.
      call dof_locGlobMapping_mult(p_rdiscretisation, p_IelementList(IELset:IELmax), &
                                  Idofs)

      ! In case ALE is used, do this also for the ALE stuff.
      if (bALE) then
        call dof_locGlobMapping_mult(p_rdiscretisation, &
                                    p_IelementList(IELset:IELmax), &
                                    IdofsALE)
      end if

      ! Calculate local DELTA`s for streamline diffusion method.
      ! (cf. p. 121 in Turek`s CFD book).
      ! For every element, we need a local DELTA.
      ! Every local delta is weighted by the global "ddelta".
      ! If ddelta=0, we do not do anything as this disables the
      ! nonlinear term.
      ! If UPSAM=0.0, we have a central-difference like discretisation, which
      ! is one can see as the local stabilisation weight Delta is also = 0.0.
      ! In this case, we even switch of the calculation of the local Delta,
      ! as it is always =0.0, so we save a little bit time.
      if ((ddelta .ne. 0.0_DP) .and. (dupsam .ne. 0.0_DP))then
        call getLocalDeltaQuad (clocalh,&
                      u1Xvel,u1Yvel,u2Xvel,u2Yvel,dweight1,dweight2,&
                      p_IelementList(IELset:IELmax),&
                      duMaxR,DlocalDelta,p_rtriangulation,Idofs,dupsam,dre)
      end if

      ! For the assembly of the global matrix, we use a "local"
      ! approach. At first we build a "local" system matrix according
      ! to the current element. This contains all additive
      ! contributions of element IEL, which are later added at the
      ! right positions to the elements in the global system matrix.
      !
      ! We have indofTrial trial DOF`s per element and
      ! indofTest test DOF`s per element. Therefore there are
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
      do IEL=1,IELmax-IELset+1

        ! For building the local matrices, we have first to
        ! loop through the test functions (the "O"`s), as these
        ! define the rows in the matrix.
        do IDOFE=1,indof

          ! Row IDOFE of the local matrix corresponds
          ! to row=global DOF KDFG(IDOFE) in the global matrix.
          ! This is one of the the "O"`s in the above picture.
          ! Get the starting position of the corresponding row
          ! to JCOL0:

          JCOL0=p_KLD(Idofs(IDOFE,IEL))

          ! Now we loop through the other DOF`s on the current element
          ! (the "O"`s).
          ! All these have common support with our current basis function
          ! and will therefore give an additive value to the global
          ! matrix.

          do JDOFE=1,indof

            ! Get the global DOF of the "X" which interacts with
            ! our "O".

            JDFG=Idofs(JDOFE,IEL)

            ! Starting in JCOL0 (which points to the beginning of
            ! the line initially), loop through the elements in
            ! the row to find the position of column IDFG.
            ! Jump out of the do loop if we find the column.

            do JCOL=JCOL0,rmatrix%NA
              if (p_KCOL(JCOL) .eq. JDFG) exit
            end do

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

          end do ! IDOFE

        end do ! JDOFE

      end do ! IEL

      ! Ok, we found the positions of the local matrix entries
      ! that we have to change.
      ! To calculate the matrix contributions, we have to evaluate
      ! the elements to give us the values of the basis functions
      ! in all the DOF`s in all the elements in our set.
      !
      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag.
      cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)

      ! In the first loop, calculate the coordinates on the reference element.
      ! In all later loops, use the precalculated information.
      if (IELset .eq. 1) then
        cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REFPOINTS)
      else
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
      end if

      ! Calculate all information that is necessary to evaluate the finite element
      ! on all cells of our subset. This includes the coordinates of the points
      ! on the cells.
      call elprep_prepareSetForEvaluation (revalElementSet,&
          cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
          ctrafoType, p_DcubPtsRef(:,1:ncubp), rperfconfig=rperfconfig)
      p_Ddetj => revalElementSet%p_Ddetj

      ! Calculate the values of the basis functions.
      call elem_generic_sim2 (p_relementDistribution%celement, &
          revalElementSet, Bder, Dbas)

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
      ! DOF`s as:
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

      if (dweight2 .eq. 0.0_DP) then

        ! Loop over all elements in the current set
        do IEL=1,IELmax-IELset+1

          ! Loop over all cubature points on the current element
          do ICUBP = 1, ncubp

            du1loc = 0.0_DP
            du2loc = 0.0_DP

            ! Perform a loop through the trial DOF`s.
            do JDOFE=1,indof

              ! Get the value of the (test) basis function
              ! phi_i (our "O") in the cubature point:
              db = Dbas(JDOFE,1,ICUBP,IEL)

              ! Sum up to the value in the cubature point
              JDFG = Idofs(JDOFE,IEL)
              du1loc = du1loc + u1Xvel(JDFG)*db
              du2loc = du2loc + u1Yvel(JDFG)*db

            end do ! JDOFE

            ! Save the computed velocity
            Dvelocity(1,ICUBP,IEL) = dweight1*du1loc
            Dvelocity(2,ICUBP,IEL) = dweight1*du2loc

          end do ! ICUBP

        end do ! IEL

      else

        do IEL=1,IELmax-IELset+1

          ! Loop over all cubature points on the current element
          do ICUBP = 1, ncubp

            du1loc = 0.0_DP
            du2loc = 0.0_DP

            ! Perform a loop through the trial DOF`s.
            do JDOFE=1,indof

              ! Get the value of the (trial) basis function
              ! phi_i in the cubature point:
              db = Dbas(JDOFE,1,ICUBP,IEL)

              ! Sum up to the value in the cubature point
              JDFG = Idofs(JDOFE,IEL)
              du1loc = du1loc + (dweight1*u1Xvel(JDFG) + dweight2*u2Xvel(JDFG))*db
              du2loc = du2loc + (dweight1*u1Yvel(JDFG) + dweight2*u2Yvel(JDFG))*db

            end do ! JDOFE

            ! Save the computed velocity
            Dvelocity(1,ICUBP,IEL) = du1loc
            Dvelocity(2,ICUBP,IEL) = du2loc

          end do ! ICUBP

        end do ! IEL

      end if

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

      if (bALE) then

        ! Calculate the values of the basis functions in all the points
        ! on all the elements
        call elem_generic_sim2 (EL_Q1, &
            revalElementSet, Bder, DbasALE)

        ! Loop over all elements in the current set
        do IEL=1,IELmax-IELset+1

          ! Loop over all cubature points on the current element
          do ICUBP = 1, ncubp

            du1loc = 0.0_DP
            du2loc = 0.0_DP

            ! Perform a loop through the trial DOF`s.
            do JDOFE=1,indof

              ! Get the value of the (trial) basis function
              db= Dbas(JDOFE,1,ICUBP,IEL)

              ! Sum up to the value in the cubature point
              JDFG = IdofsALE(IDOFE,IEL)
              du1loc = du1loc + DmeshVelocity(1,JDFG)*db
              du2loc = du2loc + DmeshVelocity(2,JDFG)*db

            end do ! JDOFE

            ! Save the computed velocity
            Dvelocity(1,ICUBP,IEL) = Dvelocity(1,ICUBP,IEL) - du1loc
            Dvelocity(2,ICUBP,IEL) = Dvelocity(2,ICUBP,IEL) - du2loc

          end do ! ICUBP

        end do ! IEL

      end if

      ! Ok, we now use Dvelocity as coefficient array in the assembly
      ! of a biliinear form!
      !
      ! Loop over the elements in the current set.

      do IEL=1,IELmax-IELset+1

        ! Clear the local matrix
        Dentry = 0.0_DP

        ! Loop over all cubature points on the current element
        do ICUBP = 1, ncubp

          ! Calculate the current weighting factor in the cubature formula
          ! in that cubature point.
          !
          ! Normally, we have to take the absolut value of the determinant
          ! of the mapping here!
          ! In 2D, the determinant is always positive, whereas in 3D,
          ! the determinant might be negative -- that is normal!
          ! But because this routine only works in 2D, we can skip
          ! the ABS here!

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
          ! Outer loop over the DOF`s i=1..indof on our current element,
          ! which corresponds to the basis functions Phi_i:

          do IDOFE=1,indof

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

            ! Inner loop over the DOF`s j=1..indof, which corresponds to
            ! the basis function Phi_j:

            do JDOFE=1,indof

              if (IDOFE.eq.JDOFE) then

                ! Short version of the evaluation of the matrix
                ! contribution - see below for a more detailed
                ! description what is added together here!

                AH = ddelta*HSUMI*(DlocalDelta(IEL)*HSUMI+HBASI1) &
                    + dny*(HBASI2**2+HBASI3**2) &
                    + dalpha*HBASI1**2

              else

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
                ! we do not have to worry about that.

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

                AH = ddelta*HSUMJ*(DlocalDelta(IEL)*HSUMI+HBASI1) &
                    + dny*(HBASI2*HBASJ2+HBASI3*HBASJ3) &
                    + dalpha*HBASI1*HBASJ1

              end if ! (IDOFE.EQ.JDOFE)

              ! Weighten the calculated value AH by the cubature
              ! weight OM and add it to the local matrix. After the
              ! loop over all DOF`s is finished, each entry contains
              ! the calculated integral.

              Dentry(JDOFE,IDOFE) = Dentry(JDOFE,IDOFE)+OM*AH

            end do ! IDOFE

          end do ! JDOFE

        end do ! ICUBP

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

        if (iand(cdef,CONV_MODMATRIX) .ne. 0) then
          do IDOFE=1,indof
            do JDOFE=1,indof
              p_DA(Kentry(JDOFE,IDOFE,IEL)) = p_DA(Kentry(JDOFE,IDOFE,IEL)) + &
                dtheta * Dentry(JDOFE,IDOFE)
            end do
          end do
        end if

        ! For cdef containing CONV_MODDEFECT, build the defect vector
        !     D = RHS - A*U
        ! This is done matrix free, only with the help of the local
        ! matrix.
        ! In this case, D=(D1,D2) is expected to be the RHS on
        ! entry and will be updated to be the defect vector when
        ! this routine is left.

        if (iand(cdef,CONV_MODDEFECT) .ne. 0) then
          do IDOFE=1,indof

            IDFG=Idofs(IDOFE,IEL)

            do JDOFE=1,indof

              denth = dtheta*Dentry(JDOFE,IDOFE)

              JDFG=Idofs(JDOFE,IEL)
              Ddef1(IDFG)= Ddef1(IDFG) - denth*Du1(JDFG)
              Ddef2(IDFG)= Ddef2(IDFG) - denth*Du2(JDFG)

            end do
          end do
        end if

      end do ! IEL

    end do ! IELset

    ! Release memory
    call elprep_releaseElementSet(revalElementSet)

    deallocate(p_DcubPtsRef)
    deallocate(Domega)
    deallocate(DlocalDelta)
    deallocate(Dvelocity)
    deallocate(Dentry)
    deallocate(Kentry)
    deallocate(IdofsALE)
    deallocate(Idofs)
    deallocate(DbasALE)
    deallocate(Dbas)

    ! Release the assembly structure if necessary.
    if (.not. present(rcubatureInfo)) then
      call spdiscr_releaseCubStructure(rtempCubatureInfo)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine conv_streamlineDiffusionBlk2d ( &
                           rvecPrimary, rvecSecondary, dprimWeight, dsecWeight,&
                           rconfig, cdef, &
                           rmatrix, rsolution, rdefect, DmeshVelocity, &
                           rcubatureInfo,rperfconfig)

!<description>
  ! Standard streamline diffusion method to set up the operator
  ! <tex>
  ! $$ dtheta  *  (                dalpha * MASS
  !                  +              dbeta * STOKES
  !                  +             ddelta * u_1 * grad(.)
  !                  +            dnewton * (.) * grad(u_1)
  !                  +   ddeltaTransposed * grad(.)^T * u_1
  !                  +  dnewtonTransposed * grad(u_1)^T * (.) ) $$
  ! </tex>
  ! in a matrix or to build a defect vector.
  ! 2D-version (X- and Y-velocity).
  !
  ! rvecPrimary, rvecSecondary are two velocity field vectors for the X-
  ! and Y-veclocity; IvelocityComp defines which components of these
  ! vectors contains the X- and which contains the Y-velocity.
  ! The final velocity vector field is then computed as a weighted average
  ! of these two:
  ! <tex>
  !  $$ u_1  =  dprimWeight * rvecPrimary  +  dsecWeight * rvecSecondary $$
  ! </tex>
  ! <tex>$ u_2 = rsolution(.) $</tex> defines the second velocity field inside of
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
  !
  ! The routine works with a block matrix rmatrix and allows to include
  ! extended operators like the Newton operator (Frechet-derivative of the
  ! convective part).
!</description>

!<input>
  ! Primary velocity field for the computation of <tex>$ u_1 $</tex>
  type(t_vectorBlock), intent(in), target :: rvecPrimary

  ! Secondary velocity field for the computation of <tex>$ u_1 $</tex>
  type(t_vectorBlock), intent(in), target :: rvecSecondary

  ! Weighting factor for rvecPrimary.
  real(DP), intent(in) :: dprimWeight

  ! Weighting factor for rvecSecondary.
  real(DP), intent(in) :: dsecWeight

  ! Configuration block for the streamline diffusion scheme
  type(t_convStreamlineDiffusion), intent(in) :: rconfig

  ! Computation/defect correction method. One of the CONV_MODxxxx constants:
  ! CONV_MODMATRIX: Set up the nonlinear matrix. rmatrix must be present, the
  !                 nonlinear part is added to the matrix.
  ! CONV_MODDEFECT: Set up the nonlinear defect. rdefect and rsolution must be
  !                 present.
  ! CONV_MODBOTH  : Set up the nonlinear matrix as well as the nonlinear defect.
  !                 rmatrix, rdefect and rsolution must all be present.
  integer, intent(in) :: cdef

  ! OPTIONAL: Solution vector u_2.
  ! Must be present if cdef=CONV_MODDEFECT or =CONV_MODBOTH.
  type(t_vectorBlock), intent(in), target, optional :: rsolution

  ! OPTIONAL: Mesh velocity field.
  ! DmeshVelocity(1,ivt) gives the X-velocity of the mesh, i.e. the X-velocity
  !   of the corner vertex ivt.
  ! DmeshVelocity(2,ivt) gives the Y-velocity of the mesh, i.e. the Y-velocity
  !   of the corner vertex ivt.
  ! The parameter must be present if ALE is activated in the
  ! configuration parameter block by bALE=true.
  real(DP), dimension(:,:), intent(in), optional :: DmeshVelocity

  ! OPTIONAL: A scalar cubature information structure that specifies the cubature
  ! formula(s) to use. If not specified, default settings are used.
  type(t_scalarCubatureInfo), intent(in), optional, target :: rcubatureInfo

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! System block matrix.
  ! The content of the matrix must be present if cdef=CONV_MODMATRIX or
  ! =CONV_MODBOTH, otherwise only the structure is used.
  ! The nonlinear operator is added to the matrix.
  ! The blocks A11,A12,A21 and A22 of this matrix are tackled by streamline
  ! diffusion.
  type(t_matrixBlock), intent(inout) :: rmatrix

  ! OPTIONAL: Defect vector.
  ! Must have the same structure as rsolution/rvecPrimary/rvecSecondary.
  ! Must be present if cdef=CONV_MODDEFECT or =CONV_MODBOTH.
  ! The nonlinear part is subtracted from this vector:
  ! <tex>$ r = r - \theta * u_1*grad(u_2) $</tex>
  type(t_vectorBlock), intent(inout), optional, target :: rdefect
!</inputoutput>

!</subroutine>

    ! local variables
    integer(I32) :: celement
    type(t_vectorScalar), pointer :: p_rvelX1,p_rvelX2,p_rvelY1,p_rvelY2
    type(t_vectorScalar), pointer :: p_rsolX,p_rsolY,p_rdefectX,p_rdefectY
    real(DP), dimension(:), pointer :: p_DvelX1,p_DvelX2,p_DvelY1,p_DvelY2
    real(DP), dimension(:), pointer :: p_DsolX,p_DsolY,p_DdefectX,p_DdefectY

    ! At first check the input parameters that everything is present what
    ! we need:
    if ((cdef .eq. CONV_MODDEFECT) .or. (cdef .eq. CONV_MODBOTH)) then
      if ((.not. present(rsolution)) .or. (.not. present(rdefect))) then
        call output_line ("Solution/defect vector not present!", &
            OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusionBlk2d")
        call sys_halt()
      end if
    end if

    if (rconfig%bALE) then
      if (.not. present(DmeshVelocity)) then
        call output_line ("Mesh velocity vector not present!", &
            OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusionBlk2d")
        call sys_halt()
      end if
    end if

    ! Get the actual subvectors from the velocity vectors that define
    ! the X- and Y-velocity.
    p_rvelX1 => rvecPrimary%RvectorBlock(1)
    p_rvelY1 => rvecPrimary%RvectorBlock(2)
    p_rvelX2 => rvecSecondary%RvectorBlock(1)
    p_rvelY2 => rvecSecondary%RvectorBlock(2)

    if (present(rsolution)) then
      p_rsolX => rsolution%RvectorBlock(1)
      p_rsolY => rsolution%RvectorBlock(2)
    else
      nullify(p_rsolX)
      nullify(p_rsolY)
    end if

    if (present(rdefect)) then
      p_rdefectX => rdefect%RvectorBlock(1)
      p_rdefectY => rdefect%RvectorBlock(2)
    else
      nullify(p_rdefectX)
      nullify(p_rdefectY)
    end if

    ! At the moment, we only support a rather limited set of configurations:
    ! Matrix and vectors must all be double precision, matrix must be format
    ! 7 or 9, discretisation must be Q1~, constant viscosity.
    if ((rmatrix%RmatrixBlock(1,1)%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
        (rmatrix%RmatrixBlock(1,1)%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
      call output_line ("Unsupported matrix format!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusionBlk2d")
      call sys_halt()
    end if

    if ((rmatrix%nblocksPerRow .ne. 1) .or. (rmatrix%nblocksPerCol .ne. 1)) then

      if ((rmatrix%RmatrixBlock(2,2)%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
          (rmatrix%RmatrixBlock(2,2)%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
        call output_line ("Unsupported matrix format!", &
            OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusionBlk2d")
        call sys_halt()
      end if

      if (lsysbl_isSubmatrixPresent(rmatrix,1,2) .and. &
          (rmatrix%RmatrixBlock(1,2)%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
          (rmatrix%RmatrixBlock(1,2)%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
        call output_line ("Unsupported matrix format!", &
            OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusionBlk2d")
        call sys_halt()
      end if

      if (lsysbl_isSubmatrixPresent(rmatrix,2,1) .and. &
          (rmatrix%RmatrixBlock(2,1)%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
          (rmatrix%RmatrixBlock(2,1)%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
        call output_line ("Unsupported matrix format!", &
            OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusionBlk2d")
        call sys_halt()
      end if

    end if

    ! If Newton must be calculated, make sure A12 and A21 exists and that
    ! all A11, A12, A21 and A22 are independent of each other!
    if ((rconfig%dnewton .ne. 0.0_DP) .or. (rconfig%dnewtonTransposed .ne. 0.0_DP)) then
      if (.not. lsysbl_isSubmatrixPresent(rmatrix,1,2) .or. &
          .not. lsysbl_isSubmatrixPresent(rmatrix,2,1)) then
        call output_line ("For the Newton matrix, A12 and A21 must be defined!", &
            OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusionBlk2d")
        call sys_halt()
      end if
      if (lsyssc_isMatrixContentShared ( &
              rmatrix%RmatrixBlock(1,1),rmatrix%RmatrixBlock(2,2)) .or. &
          lsyssc_isMatrixContentShared ( &
              rmatrix%RmatrixBlock(1,2),rmatrix%RmatrixBlock(2,1)) ) then
        call output_line ("For the Newton matrix, the matrix blocks must be indepentent!", &
            OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusionBlk2d")
        call sys_halt()
      end if
    end if

    celement = rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest% &
                RelementDistr(1)%celement
    if (rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%ccomplexity &
        .ne. SPDISC_UNIFORM) then
      call output_line ("Unsupported discretisation!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusionBlk2d")
      call sys_halt()
    end if

    if ((rvecPrimary%cdataType .ne. ST_DOUBLE) .or. &
        (rvecSecondary%cdataType .ne. ST_DOUBLE)) then
      call output_line ("Unsupported vector data type in velocity!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusionBlk2d")
      call sys_halt()
    end if

    if (present(rdefect)) then
      if ((rsolution%cdataType .ne. ST_DOUBLE) .or. &
          (rdefect%cdataType .ne. ST_DOUBLE)) then
        call output_line ("Unsupported vector data type in solution/defect!", &
            OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusionBlk2d")
        call sys_halt()
      end if
    end if

    if (.not. rconfig%bconstViscosity) then
      call output_line ("Only constant viscosity supported at the moment!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusionBlk2d")
      call sys_halt()
    end if

    if (rconfig%dnu .eq. SYS_INFINITY_DP) then
      call output_line ("Viscosity parameter nu not initialised!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusionBlk2d")
      call sys_halt()
    end if

    ! Hide the p_rsol...-parameters to prevent passing the NULL()-pointer
    ! if rsolution is not present -- some compilers do not like that ^^

    call lsyssc_getbase_double (p_rvelX1,p_DvelX1)
    call lsyssc_getbase_double (p_rvelY1,p_DvelY1)
    call lsyssc_getbase_double (p_rvelX2,p_DvelX2)
    call lsyssc_getbase_double (p_rvelY2,p_DvelY2)

    !!! DEBUG:
    !WHERE (abs(p_DvelX1) .LT. 1E-12_DP) p_DvelX1 = 0.0_DP
    !WHERE (abs(p_DvelY1) .LT. 1E-12_DP) p_DvelY1 = 0.0_DP
    !call vecio_writeArray_Dble (p_DvelX1, 'vecx1', &
    !                               0, 'vectorx1.txt', '(D10.3)')
    !call vecio_writeArray_Dble (p_DvelY1, 'vecx2', &
    !                               0, 'vectorx2.txt', '(D10.3)')

    if (present(rdefect)) then
      call lsyssc_getbase_double (p_rsolX   ,p_DsolX   )
      call lsyssc_getbase_double (p_rsolY   ,p_DsolY   )
      call lsyssc_getbase_double (p_rdefectX,p_DdefectX)
      call lsyssc_getbase_double (p_rdefectY,p_DdefectY)

      call conv_strdiff2dALEblk_double ( &
                    p_DvelX1,p_DvelY1,p_DvelX2,p_DvelY2,dprimWeight,dsecWeight, &
                    rmatrix,cdef, rconfig%dupsam, rconfig%dnu, &
                    rconfig%dalpha, rconfig%dbeta, rconfig%dtheta, rconfig%ddelta, &
                    rconfig%dnewton, rconfig%ddeltaTransposed, &
                    rconfig%dnewtonTransposed, rconfig%bALE, rconfig%clocalh,&
                    p_DsolX,p_DsolY,p_DdefectX,p_DdefectY, DmeshVelocity, &
                    rcubatureInfo,rperfconfig)

    else

      call conv_strdiff2dALEblk_double ( &
                    p_DvelX1,p_DvelY1,p_DvelX2,p_DvelY2,dprimWeight,dsecWeight, &
                    rmatrix, cdef, rconfig%dupsam, rconfig%dnu, &
                    rconfig%dalpha, rconfig%dbeta, rconfig%dtheta, rconfig%ddelta, &
                    rconfig%dnewton, rconfig%ddeltaTransposed, rconfig%dnewtonTransposed, &
                    rconfig%bALE, rconfig%clocalh,DmeshVelocity=DmeshVelocity,&
                    rcubatureInfo=rcubatureInfo,rperfconfig=rperfconfig)

      !!! DEBUG:
      !call matio_writeMatrixHR (rmatrix, 'matrix',&
      !                          .TRUE., 0, 'matrixL.txt', '(D10.3)')

    end if

  end subroutine

  ! ***************************************************************************

!                denth = dtheta*Dentry(JDOFE,IDOFE,IEL)
!
!                JDFG=Idofs(JDOFE,IEL)
!                Ddef1(IDFG)= Ddef1(IDFG) - denth*Du1(JDFG)
!                Ddef2(IDFG)= Ddef2(IDFG) - denth*Du2(JDFG)

!<subroutine>
  subroutine conv_strdiff2dALEblk_double ( &
                  u1Xvel,u1Yvel,u2Xvel,u2Yvel,dweight1,dweight2,&
                  rmatrix,cdef, &
                  dupsam,dnu,dalpha,dbeta,dtheta, ddelta, dnewton, &
                  ddeltaTransposed,dnewtonTransposed, bALE, clocalh, &
                  Du1,Du2,Ddef1,Ddef2, DmeshVelocity, &
                  rcubatureInfo, rperfconfig)
!<description>
  ! Standard streamline diffusion method to set up the operator
  ! <tex>
  ! $$ dtheta  *  (                dalpha * MASS
  !                  +              dbeta * STOKES
  !                  +             ddelta * u_1 * grad(.)
  !                  +            dnewton * (.) * grad(u_1)
  !                  +   ddeltaTransposed * grad(.)^T * u_1
  !                  +  dnewtonTransposed * grad(u_1)^T * (.) ) $$
  ! </tex>
  ! in a matrix or to build a defect vector with that.
  ! 2D-version (X- and Y-velocity), uniform <tex>$\tilde Q_1$</tex> discretisation,
  ! double precision vectors/matrix.
  !
  ! The routine supports fully coupled matrices, and the generation of the Newton
  ! matrix.
  !
  ! u1Xvel,u1Yvel, u2Xvel,u2Yvel are two velocity field vectors,
  ! (u1Xvel,u1Yvel) a primary and (u2Xvel,u2Yvel) a secondary velocity field.
  ! The final velocity vector field is then computed as a weighted average
  ! of these two:
  ! <tex>
  !  $$ u_1  =  dweight1 * u1vel  +  dweight2 * u2vel $$
  ! </tex>
  ! <tex>$ u_2 = rsolution(.) $</tex> defines the second velocity field inside of
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
  ! [Duarte, Formaz, Natesan; `Arbitrary Lagrangian-Euler Method
  ! for Navier-Stokes equations with moving boundaries`;
  ! Comput. Methods Appl. Mech. Engrg. 193 (2004), 4819-4836]
  !
  ! Remarks:\\
  !
  ! 1.) In a typical call of the upwinding, the caller can use:
  !     dweight1 = 1, u1Xvel/u1Yvel = velocity field
  !     dweight2 = 0, u2Xvel/u2Yvel = undefined
  !   So the upwinding scheme only uses one velocity field.
  !   Such a call e.g. adds the integral
  !                <tex> $$ ( u_1 * grad(.) , v )_{\Omega} $$ </tex>
  !   to the system matrix.\\
  !
  !  2.) In case that there are two velocity fields representing
  !   the solution (may happen in a nonstationary simulation where
  !   u1Xvel/u1Yvel represents the solution in the current and u2Xvel/u2Yvel
  !   that of the previous time step), dweight1/dweight2 defines how these both
  !   velocity vectors should be weighted to compute the actual
  !   velocity field for the assembling:
  !               <tex> $$ U_act = dweight1*u1vel + dweight2*u2vel $$ </tex>
  !   This is e.g. used for the linear extrapolation technique to
  !   reconstruct a velocity from two previous time steps...\\
  !
  !  3.) In the nonlinear iteration, as a right hand side there arises
  !   a defect vector D, which linear part can easily being assembled.
  !   However, there is a nonlinearity to be included into that vector,
  !   too. By setting cdef=1,2, this routine incorporates the nonlinearity
  !   into that vector, using the formula
  !
  !            <tex> $$ D = D - dtheta * UUx * grad (Ux) $$ </tex>
  !
  !  4.) If bALE=true, a mesh velocity field is added to the nonlineareity
  !   according to the formula  "U * grad (U-DmeshVelocity)".
  !   For bALE=false, the simple nonlinearity "U * grad (U)" is used.

!</description>

!<input>

  ! Primary X-velocity of <tex>$ u_1 $</tex>
  real(DP), dimension(:), intent(in) :: u1Xvel

  ! Primary Y-velocity of <tex>$ u_1 $</tex>
  real(DP), dimension(:), intent(in) :: u1Yvel

  ! Secondary X-velocity of <tex>$ u_1 $</tex>
  real(DP), dimension(:), intent(in) :: u2Xvel

  ! Secondary Y-velocity of <tex>$ u_1 $</tex>
  real(DP), dimension(:), intent(in) :: u2Yvel

  ! Computation/defect correction method. One of the CONV_MODxxxx constants:
  ! CONV_MODMATRIX: Set up the nonlinear matrix. rmatrix must be present, the
  !                 nonlinear part is added to the matrix.
  ! CONV_MODDEFECT: Set up the nonlinear defect. rdefect and rsolution must be
  !                 present.
  ! CONV_MODBOTH  : Set up the nonlinear matrix as well as the nonlinear defect.
  !                 rmatrix, rdefect and rsolution must all be present.
  integer, intent(in) :: cdef

  ! Weighting factor for u1Xvel/u1Yvel.
  real(DP), intent(in) :: dweight1

  ! Weighting factor for u2Xvel/u2Yvel.
  real(DP), intent(in) :: dweight2

  ! dupsam  - control parameter.
  !          -1: simple upwind,
  !          =0: Samarskji upwind
  real(DP), intent(in) :: dupsam

  ! Viscosity parameter <tex>$ \nu = 1/Re $</tex> if viscosity is constant
  real(DP), intent(in) :: dnu

  ! Weighting factor for the mass matrix.
  real(DP), intent(in) :: dalpha

  ! Weighting factor for the Stokes matrix. (Stokes matrix = 1/Re * Laplace)
  real(DP), intent(in) :: dbeta

  ! Weighting factor of the convective operator: <tex>$ \theta * u*grad(u) $</tex>.
  ! For time-dependent problems, this can be set to the step size
  ! in the <tex>$ \Theta $</tex>-scheme.
  real(DP), intent(in) :: dtheta

  ! Weighting factor for the nonlinear term
  real(DP), intent(in) :: ddelta

  ! Weighting factor of the Newton matrix. A value of 0.0 deactivates the
  ! Newton part. A value != 0.0 activates Newton; in this case the submatrices
  ! A12 and A21 must be present in rmatrix.
  real(DP), intent(in) :: dnewton

  ! Weighting factor of the transposed convection matrix. A value of 0.0 deactivates
  ! this operator.
  real(DP), intent(in) :: ddeltaTransposed

  ! Weighting factor of the transposed Newton matrix. A value of 0.0 deactivates
  ! this operator.
  real(DP), intent(in) :: dnewtonTransposed

  ! Whether or not to use the ALE method
  logical, intent(in) :: bALE

  ! Method how to compute the local h
  integer, intent(in) :: clocalh

  ! OPTIONAL: Mesh velocity field. Must be present if bALE=TRUE.
  ! DmeshVelocity(1,:) gives the X-velocity of all the corner points of the mesh,
  ! DmeshVelocity(2,:) gives the Y-velocity.
  real(DP), dimension(:,:), intent(in), optional :: DmeshVelocity(:,:)

  ! OPTIONAL: X-velocity of <tex>$ u_2 $</tex>. Must be present if cdef=CONV_MODDEFECT
  ! or cdef=CONV_MODBOTH.
  real(DP), dimension(:), intent(in), optional :: Du1

  ! Y-velocity of <tex>$ u_2 $</tex>. Must be present if cdef=CONV_MODDEFECT
  ! or cdef=CONV_MODBOTH.
  real(DP), dimension(:), intent(in), optional :: Du2

  ! OPTIONAL: A scalar cubature information structure that specifies the cubature
  ! formula(s) to use. If not specified, default settings are used.
  type(t_scalarCubatureInfo), intent(in), optional, target :: rcubatureInfo

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! The system matrix. The submatrices for the velocity must be in block
  ! A11, A12, A21 and A22 and must be in matrix format 7 or 9.
  ! A11 and A22 must have the same structure. A12 and A21 must have
  ! the same structure.
  type(t_matrixBlock), intent(inout), target :: rmatrix

  ! OPTIONAL: X-defect vector. Must be present if cdef=CONV_MODDEFECT
  ! or =CONV_MODBOTH.
  real(DP), dimension(:), intent(inout), optional :: Ddef1

  ! OPTIONAL: Y-defect vector. Must be present if cdef=CONV_MODDEFECT
  ! or =CONV_MODBOTH.
  real(DP), dimension(:), intent(inout), optional :: Ddef2
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: indof,indofALE,IEQ,IDOFE,JDOFE,icubp
  integer :: JCOL0,IDFG,JDFG,JCOL
  integer :: IEL,IELset,IELmax
  logical, dimension(EL_MAXNDER) :: Bder,BderALE
  real(DP) :: dumax,dumaxr, du1loc, du2loc, dunorm,db,OM,AH,denth,dre,dny
  real(DP) :: du1locx,du2locx,du1locy,du2locy,dbx,dby
  real(DP) :: AH11,AH12,AH21,AH22
  real(DP) :: HBASI1,HBASI2,HBASI3,HBASJ1,HBASJ2,HBASJ3,HSUMI,HSUMJ
  integer :: NVE

  ! Matrix structure arrays
  integer, dimension(:), pointer :: p_Kcol
  integer, dimension(:), pointer :: p_Kld
  real(DP), dimension(:), pointer :: p_Da11,p_Da22

  integer, dimension(:), pointer :: p_Kcol12
  integer, dimension(:), pointer :: p_Kld12
  real(DP), dimension(:), pointer :: p_Da12,p_Da21

  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  real(DP), dimension(:,:), pointer :: p_DcubPtsRef

  ! The discretisation - for easier access
  type(t_spatialDiscretisation), pointer :: p_rdiscretisation

  ! Triangulation
  type(t_triangulation), pointer :: p_rtriangulation
  real(DP), dimension(:,:), pointer :: p_DvertexCoords
  integer, dimension(:,:), pointer :: p_IedgesAtElement,p_IverticesAtElement

  ! Number of elements in a block. Normally =NELEMSIM,
  ! except if there are less elements in the discretisation.
  integer :: nelementsPerBlock

  ! One and only element distribution
  type(t_elementDistribution), pointer :: p_relementDistribution

  ! For every cubature point on the reference element,
  ! the corresponding cubature weight
  real(DP), dimension(:), allocatable :: Domega

  ! number of cubature points on the reference element
  integer :: ncubp

  ! An element evaluation set for evaluating elements.
  type(t_evalElementSet) :: revalElementSet
  logical :: bcubPtsInitialised

  ! Arrays for saving Jacobian determinants and matrices
  real(DP), dimension(:,:), pointer :: p_Ddetj

  ! An allocateable array accepting the DOF`s of a set of elements.
  integer, dimension(:,:), allocatable, target :: Idofs, IdofsALE

  ! Allocateable arrays for the values of the basis functions -
  ! for test and trial spaces.
  real(DP), dimension(:,:,:,:), allocatable, target :: Dbas,DbasALE

  ! Local matrices, used during the assembly.
  ! Values and positions of values in the global matrix.
  integer, dimension(:,:,:), allocatable :: Kentry
  real(DP), dimension(:,:,:), allocatable :: Dentry

  ! Additional contributions for the submatrices A11, A12, A21, A22 stemming from Newton.
  integer, dimension(:,:,:), allocatable :: Kentry12
  real(DP), dimension(:,:,:), allocatable :: DentryA11
  real(DP), dimension(:,:,:), allocatable :: DentryA12
  real(DP), dimension(:,:,:), allocatable :: DentryA21
  real(DP), dimension(:,:,:), allocatable :: DentryA22

  ! A pointer to an element-number list
  integer, dimension(:), pointer :: p_IelementList

  ! Pointer to the velocity field in the cubature points.
  real(DP), dimension(:,:,:), allocatable :: Dvelocity

  ! Pointer to the velocity X- and Y-derivative in the cubature points
  real(DP), dimension(:,:,:), allocatable :: DvelocityUderiv
  real(DP), dimension(:,:,:), allocatable :: DvelocityVderiv

  ! An array with local DELTA`s, each DELTA for one element
  real(DP), dimension(:), allocatable :: DlocalDelta

  ! Type of transformation from the reference to the real element
  integer(I32) :: ctrafoType

  ! Element evaluation tag; collects some information necessary for evaluating
  ! the elements.
  integer(I32) :: cevaluationTag

  ! Cubature formula and cubature info structure
  type(t_scalarCubatureInfo), target :: rtempCubatureInfo
  type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo
  integer(I32) :: ccubature

  ! Pointer to the performance configuration
  type(t_perfconfig), pointer :: p_rperfconfig

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => conv_perfconfig
    end if

    ! Initialise the derivative flags
    Bder = .false.
    Bder(DER_FUNC) = .true.
    Bder(DER_DERIV_X) = .true.
    Bder(DER_DERIV_Y) = .true.

    ! For ALE we do not even need so much
    BderALE = .false.
    BderALE(DER_FUNC) = .true.
    BderALE(DER_DERIV_X) = .true.
    BderALE(DER_DERIV_X) = .true.

    ! Shortcut to the spatial discretisation.
    ! We assume the same for all, A11, A12, A21 and A22.
    p_rdiscretisation => rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest

    ! Get the element distribution. Here, we can find information about
    ! the cubature formula etc...
    p_relementDistribution => p_rdiscretisation%RelementDistr(1)

    ! Get some information about the triangulation
    p_rtriangulation => p_rdiscretisation%p_rtriangulation
    call storage_getbase_double2d (p_rtriangulation%h_DvertexCoords,&
                                   p_DvertexCoords)
    call storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
                                p_IverticesAtElement)
    call storage_getbase_int2d (p_rtriangulation%h_IedgesAtElement,&
                                p_IedgesAtElement)

    ! Get the number of local DOF`s for trial/test functions.
    ! We assume trial and test functions to be the same.
    indof = elem_igetNDofLoc(p_relementDistribution%celement)

    ! Get the number of local DOF`s Q1 -- we need them for ALE.
    indofALE = elem_igetNDofLoc(EL_Q1)

    ! Number of local DOF`s
    NVE = elem_igetNVE(p_relementDistribution%celement)

    ! Do we have an assembly structure?
    ! If we do not have it, create a cubature info structure that
    ! defines how to do the assembly.
    if (.not. present(rcubatureInfo)) then
      call spdiscr_createDefCubStructure(p_rdiscretisation,&
          rtempCubatureInfo,CUB_GEN_DEPR_EVAL)
      p_rcubatureInfo => rtempCubatureInfo
    else
      p_rcubatureInfo => rcubatureInfo
    end if

    ! Cubature
    ccubature = p_rcubatureInfo%p_RinfoBlocks(1)%ccubature

    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it is at most
    ! NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = min(p_rperfconfig%NELEMSIM,p_rtriangulation%NEL)

    ! For cdef containing CONV_MODDEFECT, we build the defect vector
    !     D = RHS - A*U
    ! In this case, the defect(rhs vectors must be present

    if (iand(cdef,CONV_MODDEFECT) .ne. 0) then
      if (.not. (present(Ddef1) .and. present(Ddef2) .and. &
                 present(Du1) .and. present(Du2))) then
        call output_line ("Necessary arguments missing!", &
            OU_CLASS_ERROR,OU_MODE_STD,"conv_strdiff2dALEblk_double")
        call sys_halt()
      end if
    end if

    ! Get pointers to the matrix content (if necessary)
    if (iand(cdef,CONV_MODMATRIX) .ne. 0) then
      ! Get matrix arrays
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,1),p_Da11)
      if (lsysbl_isSubmatrixPresent(rmatrix,2,2)) then
        call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,2),p_Da22)
      else
        ! Take A22=A11. In case there is only one block,
        ! a later check ensures that A22 is not modified if
        ! p_Da22 == p_Da11, so that's ok here.
        p_Da22 => p_Da11
      end if

      if (lsysbl_isSubmatrixPresent(rmatrix,1,2)) then
        call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,2),p_Da12)
        call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,1),p_Da21)
      else
        nullify(p_Da12,p_Da21)
      end if
    end if

    ! Get pointers to the matrix structure(s).
    call lsyssc_getbase_Kcol (rmatrix%RmatrixBlock(1,1),p_Kcol)
    call lsyssc_getbase_Kld (rmatrix%RmatrixBlock(1,1),p_Kld)

    if (lsysbl_isSubmatrixPresent(rmatrix,1,2)) then
      call lsyssc_getbase_Kcol (rmatrix%RmatrixBlock(1,2),p_Kcol12)
      call lsyssc_getbase_Kld (rmatrix%RmatrixBlock(1,2),p_Kld12)
    else
      nullify(p_Kcol12,p_Kld12)
    end if

    ! Get from the trial element space the type of coordinate system
    ! that is used there:
    ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)

    ! Get the number of cubature points for the cubature formula
    ncubp = cub_igetNumPts(ccubature)

    ! Allocate two arrays for the points and the weights
    allocate(Domega(ncubp))
    allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),ncubp))

    ! Get the cubature formula
    call cub_getCubature(ccubature,p_DcubPtsRef, Domega)

    ! OpenMP-Extension: Open threads here.
    ! "csysTrial" is declared as private; shared gave errors with the Intel compiler
    ! in Windows!?!
    ! Each thread will allocate its own local memory...

    !%omp parallel private( &
    !%omp p_Ddetj, i,k,Dbas,Idofs,DbasALE, &
    !%omp IdofsALE,DlocalDelta,Kentry,Kentry12,Dentry, &
    !%omp DentryA11,DentryA12,DentryA21,DentryA22,Dvelocity, &
    !%omp DvelocityUderiv,DvelocityVderiv,dre,IEL,db,icubp,&
    !%omp IDOFE,JCOL0,JDOFE,JDFG,jcol,du1loc,du2loc,dbx,dby, &
    !%omp du1locx,du1locy,du2locx,du2locy,OM,AH,HBASI1,HBASI2,&
    !%omp HBASI3,HBASJ1,HBASJ2,HBASJ3,HSUMI,HSUMJ,AH11,AH12,AH21, &
    !%omp AH22,IELmax,revalElementSet,dny)

    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  allocate(Dbas(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    ! would lead to nonused memory blocks in these arrays during the assembly,
    ! which reduces the speed by 50%!

    allocate(Dbas(indof,elem_getMaxDerivative(p_relementDistribution%celement), &
             ncubp,nelementsPerBlock))

    ! Allocate memory for the DOF`s of all the elements.
    allocate(Idofs(indof,nelementsPerBlock))

    ! The same for the ALE-space
    allocate(DbasALE(indofALE,elem_getMaxDerivative(EL_Q1), &
             ncubp,nelementsPerBlock))

    ! Allocate memory for the DOF`s of all the elements.
    allocate(IdofsALE(indofALE,nelementsPerBlock))

    ! Allocate memory for array with local DELTA`s
    allocate(DlocalDelta(nelementsPerBlock))

    ! Allocate an array saving the local matrices for all elements
    ! in an element set.
    ! We could also allocate EL_MAXNBAS*EL_MAXNBAS*NELEMSIM integers
    ! for this local matrix, but this would normally not fit to the cache
    ! anymore! indofTrial*indofTest*NELEMSIM is normally much smaller!
    !
    ! Kentry (:,:,:) defines the positions of the local matrices
    ! in the submatrices A11 and A22.
    ! KentryA12 (:,:,:) defines the positions of the local matrices
    ! in the submatrices A12 and A21.
    allocate(Kentry(indof,indof,nelementsPerBlock))

    if ((dnewton .ne. 0.0_DP) .or. (dnewtonTransposed .ne. 0.0_DP) &
        .or. (ddeltaTransposed .ne. 0.0_DP)) then
      allocate(Kentry12(indof,indof,nelementsPerBlock))
    end if

    ! Dentry (:,:,:) fetches the 'main' matrix entries (Laplace, Mass,
    ! Convection).
    ! DentryA11, DentryA12, DentryA21 and DentryA22 fetches additional entries in
    ! A11, A12, A21 and A22 of the Newton matrix, which is not always calculated
    ! and therefore not always used!
    allocate(Dentry(indof,indof,nelementsPerBlock))

    if ((dnewton .ne. 0.0_DP) .or. (dnewtonTransposed .ne. 0.0_DP) &
        .or. (ddeltaTransposed .ne. 0.0_DP)) then
      allocate(DentryA11(indof,indof,nelementsPerBlock))
      allocate(DentryA12(indof,indof,nelementsPerBlock))
      allocate(DentryA21(indof,indof,nelementsPerBlock))
      allocate(DentryA22(indof,indof,nelementsPerBlock))
    end if

    ! Allocate memory for the velocity in the cubature points.
    allocate(Dvelocity(NDIM2D,ncubp,nelementsPerBlock))

    ! Initialisation of the element set.
    call elprep_init(revalElementSet)

    ! Indicate that cubature points must still be initialised in the element set.
    bcubPtsInitialised = .false.

    if ((dnewton .ne. 0.0_DP) .or. (dnewtonTransposed .ne. 0.0_DP)) then
      allocate(DvelocityUderiv(NDIM2D,ncubp,nelementsPerBlock))
      allocate(DvelocityVderiv(NDIM2D,ncubp,nelementsPerBlock))
    end if

    ! What is the reciprocal of nu? We need it later.
    if (dnu .ne. 0.0_DP) then
      dre = 1.0_DP/dnu

      ! dny gets the actual multiplier for the Laplace matrix.
      ! Remember: dbeta*Stokes = dbeta*dnu*Laplace = dny*Laplace.
      ! This may be =0.0 if the Stokes operator should not be included into
      ! the matrix.
      dny = dbeta*dnu
    else
      call output_line ("NU=0 not allowed! Set dbeta=0 to prevent Stokes operator "//&
          "from being build!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_strdiff2dALEblk_double")
      call sys_halt()
    end if

    ! If ddelta=0, we have to neglect the nonlinearity. In both cases,
    ! set DlocalDelta=0 which disables the nonlinear term in the assembly.
    ! If dupsam=0, we neglect the stabilisation term (central difference like
    ! discretisation), so we set DlocalDelta=0 as well.
    if ((ddelta .eq. 0.0_DP) .or. (dupsam .eq. 0.0_DP)) then
      call lalg_clearVectorDble (DlocalDelta)
    end if

    ! Calculate the maximum norm of the actual velocity field
    ! U = A1*U1 + A2*U2 into DUMAX.
    ! Round up the norm to 1D-8 if it is too small...
    !%omp single
    dumax=0.0_DP
    if (dweight2 .eq. 0.0_DP) then


      do IEQ=1,size(u1Xvel)
        du1loc = dweight1*u1Xvel(IEQ)
        du2loc = dweight1*u1Yvel(IEQ)
        dunorm = sqrt(du1loc**2+du2loc**2)
        dumax = max(DUMAX,DUNORM)
      end do

    else

      do ieq=1,size(u1Xvel)
        du1loc = dweight1*u1Xvel(IEQ)+dweight2*u2Xvel(IEQ)
        du2loc = dweight1*u1Yvel(IEQ)+dweight2*u2Yvel(IEQ)
        dunorm = sqrt(du1loc**2+du2loc**2)
        dumax = max(dumax,dunorm)
      end do

    end if

    !print *,"dumax: ",dumax
    if (dumax.lt.1E-8_DP) dumax=1E-8_DP
    dumaxr = 1.0_DP/dumax
    !%omp end single

    ! p_IelementList must point to our set of elements in the discretisation
    ! with that combination of trial/test functions
    call storage_getbase_int (p_relementDistribution%h_IelementList, &
                              p_IelementList)


    ! Loop over the elements - blockwise.
    !
    ! OpenMP-Extension: Each loop cycle is executed in a different thread,
    ! so NELEMSIM local matrices are simultaneously calculated in the
    ! inner loop(s).
    ! The blocks have all the same size, so we can use static scheduling.
    !%omp do schedule(dynamic,1)
    do IELset = 1, size(p_IelementList), nelementsPerBlock

      ! We always handle NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most NELEMSIM
      ! elements simultaneously.

      IELmax = min(size(p_IelementList),IELset-1+p_rperfconfig%NELEMSIM)

      ! The outstanding feature with finite elements is: A basis
      ! function for a DOF on one element has common support only
      ! with the DOF`s on the same element! E.g. for Q1:
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
      !     to collect all "O" DOF`s.
      !
      ! Calculate the global DOF`s into IdofsTrial / IdofsTest.
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF`s of our NELEMSIM elements simultaneously.
      call dof_locGlobMapping_mult(p_rdiscretisation, p_IelementList(IELset:IELmax), &
                                  Idofs)

      ! In case ALE is used, do this also for the ALE stuff.
      if (bALE) then
        call dof_locGlobMapping_mult(p_rdiscretisation, &
                                    p_IelementList(IELset:IELmax), &
                                    IdofsALE)
      end if

      ! Calculate local DELTA`s for streamline diffusion method.
      ! (cf. p. 121 in Turek`s CFD book).
      ! For every element, we need a local DELTA.
      ! Every local delta is weighted by the global "ddelta".
      ! If ddelta=0, we do not do anything as this disables the
      ! nonlinear term.
      ! If UPSAM=0.0, we have a central-difference like discretisation, which
      ! is one can see as the local stabilisation weight Delta is also = 0.0.
      ! In this case, we even switch of the calculation of the local Delta,
      ! as it is always =0.0, so we save a little bit time.
      if ((ddelta .ne. 0.0_DP) .and. (dupsam .ne. 0.0_DP))then
        call getLocalDeltaQuad (clocalh,&
                      u1Xvel,u1Yvel,u2Xvel,u2Yvel,dweight1,dweight2,&
                      p_IelementList(IELset:IELmax),&
                      duMaxR,DlocalDelta,p_rtriangulation,Idofs,dupsam,dre)
      end if

      ! For the assembly of the global matrix, we use a "local"
      ! approach. At first we build a "local" system matrix according
      ! to the current element. This contains all additive
      ! contributions of element IEL, which are later added at the
      ! right positions to the elements in the global system matrix.
      !
      ! We have indofTrial trial DOF`s per element and
      ! indofTest test DOF`s per element. Therefore there are
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
      do IEL=1,IELmax-IELset+1

        ! For building the local matrices, we have first to
        ! loop through the test functions (the "O"`s), as these
        ! define the rows in the matrix.
        do IDOFE=1,indof

          ! Row IDOFE of the local matrix corresponds
          ! to row=global DOF KDFG(IDOFE) in the global matrix.
          ! This is one of the the "O"`s in the above picture.
          ! Get the starting position of the corresponding row
          ! to JCOL0:

          JCOL0=p_KLD(Idofs(IDOFE,IEL))

          ! Now we loop through the other DOF`s on the current element
          ! (the "O"`s).
          ! All these have common support with our current basis function
          ! and will therefore give an additive value to the global
          ! matrix.

          do JDOFE=1,indof

            ! Get the global DOF of the "X" which interacts with
            ! our "O".

            JDFG=Idofs(JDOFE,IEL)

            ! Starting in JCOL0 (which points to the beginning of
            ! the line initially), loop through the elements in
            ! the row to find the position of column IDFG.
            ! Jump out of the do loop if we find the column.

            do JCOL=JCOL0,rmatrix%RmatrixBlock(1,1)%NA
              if (p_KCOL(JCOL) .eq. JDFG) exit
            end do

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

          end do ! IDOFE

        end do ! JDOFE

      end do ! IEL

      ! If the Newton part is to be calculated, we also need the matrix positions
      ! in A12 and A21. We can skip this part if the column structure is
      ! exactly the same!
      if ((dnewton .ne. 0.0_DP) .or. (dnewtonTransposed .ne. 0.0_DP) &
          .or. (ddeltaTransposed .ne. 0.0_DP)) then
        if (associated(p_Kcol,p_Kcol12)) then

          Kentry12(:,:,:) = Kentry(:,:,:)

        else

          do IEL=1,IELmax-IELset+1

            ! For building the local matrices, we have first to
            ! loop through the test functions (the "O"`s), as these
            ! define the rows in the matrix.
            do IDOFE=1,indof

              ! Row IDOFE of the local matrix corresponds
              ! to row=global DOF KDFG(IDOFE) in the global matrix.
              ! This is one of the the "O"`s in the above picture.
              ! Get the starting position of the corresponding row
              ! to JCOL0:

              JCOL0=p_KLD12(Idofs(IDOFE,IEL))

              ! Now we loop through the other DOF`s on the current element
              ! (the "O"`s).
              ! All these have common support with our current basis function
              ! and will therefore give an additive value to the global
              ! matrix.

              do JDOFE=1,indof

                ! Get the global DOF of the "X" which interacts with
                ! our "O".

                JDFG=Idofs(JDOFE,IEL)

                ! Starting in JCOL0 (which points to the beginning of
                ! the line initially), loop through the elements in
                ! the row to find the position of column IDFG.
                ! Jump out of the do loop if we find the column.

                do JCOL=JCOL0,rmatrix%RmatrixBlock(1,2)%NA
                  if (p_KCOL12(JCOL) .eq. JDFG) exit
                end do

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

                Kentry12(JDOFE,IDOFE,IEL)=JCOL

              end do ! IDOFE

            end do ! JDOFE

          end do ! IEL

        end if

      end if ! dnewton != 0


      ! Ok, we found the positions of the local matrix entries
      ! that we have to change.
      ! To calculate the matrix contributions, we have to evaluate
      ! the elements to give us the values of the basis functions
      ! in all the DOF`s in all the elements in our set.

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag.
      cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)
      cevaluationTag = ior(cevaluationTag,elem_getEvaluationTag(EL_Q1))

      ! In the first loop, calculate the coordinates on the reference element.
      ! In all later loops, use the precalculated information.
      !
      ! Note: Why not using
      !   if (IELset .EQ. 1) then
      ! here, but this strange concept with the boolean variable?
      ! Because the if-command does not work with OpenMP! bcubPtsInitialised
      ! is a local variable and will therefore ensure that every thread
      ! is initialising its local set of cubature points!
      if (.not. bcubPtsInitialised) then
        bcubPtsInitialised = .true.
        cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REFPOINTS)
      else
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
      end if

      ! Calculate all information that is necessary to evaluate the finite element
      ! on all cells of our subset. This includes the coordinates of the points
      ! on the cells.
      call elprep_prepareSetForEvaluation (revalElementSet,&
          cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
          ctrafoType, p_DcubPtsRef(:,1:ncubp), rperfconfig=rperfconfig)
      p_Ddetj => revalElementSet%p_Ddetj

      ! Calculate the values of the basis functions.
      ! Pass p_DcubPts as point coordinates, which point either to the
      ! coordinates on the reference element (the same for all elements)
      ! or on the real element - depending on whether this is a
      ! parametric or nonparametric element.
      call elem_generic_sim2 (p_relementDistribution%celement, &
          revalElementSet, Bder, Dbas)

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
      ! DOF`s as:
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

      ! only primary velocity field
      if (dweight2 .eq. 0.0_DP) then
!      print *,"dweight2 .EQ. 0.0"

        ! Loop over all elements in the current set
        do IEL=1,IELmax-IELset+1

          ! Loop over all cubature points on the current element
          do ICUBP = 1, ncubp

            du1loc = 0.0_DP
            du2loc = 0.0_DP

            ! Perform a loop through the trial DOF`s.
            do JDOFE=1,indof

              ! Get the value of the (test) basis function
              ! phi_i (our "O") in the cubature point:

              db = Dbas(JDOFE,1,ICUBP,IEL)

              ! Sum up to the value in the cubature point

              JDFG = Idofs(JDOFE,IEL)
              du1loc = du1loc +u1Xvel(JDFG)*db
              du2loc = du2loc +u1Yvel(JDFG)*db

            end do ! JDOFE

            ! Save the computed velocity
            Dvelocity(1,ICUBP,IEL) = dweight1*du1loc
            Dvelocity(2,ICUBP,IEL) = dweight1*du2loc

          end do ! ICUBP

        end do ! IEL

        ! Compute X- and Y-derivative of the velocity?
        if ((dnewton .ne. 0.0_DP) .or. (dnewtonTransposed .ne. 0.0_DP)) then

          do IEL=1,IELmax-IELset+1

            ! Loop over all cubature points on the current element
            do ICUBP = 1, ncubp

              du1locx = 0.0_DP
              du1locy = 0.0_DP
              du2locx = 0.0_DP
              du2locy = 0.0_DP

              ! Perform a loop through the trial DOF`s.
              do JDOFE=1,indof

                ! Get the value of the (trial) basis function
                ! phi_i in the cubature point:
                dbx = Dbas(JDOFE,DER_DERIV_X,ICUBP,IEL)
                dby = Dbas(JDOFE,DER_DERIV_Y,ICUBP,IEL)

                ! Sum up to the value in the cubature point
                JDFG = Idofs(JDOFE,IEL)
                du1locx = du1locx + u1Xvel(JDFG)*dbx
                du1locy = du1locy + u1Xvel(JDFG)*dby
                du2locx = du2locx + u1Yvel(JDFG)*dbx
                du2locy = du2locy + u1Yvel(JDFG)*dby

              end do ! JDOFE


              ! Save the computed velocity derivative
              DvelocityUderiv(1,ICUBP,IEL) = dweight1*du1locx
              DvelocityUderiv(2,ICUBP,IEL) = dweight1*du1locy
              DvelocityVderiv(1,ICUBP,IEL) = dweight1*du2locx
              DvelocityVderiv(2,ICUBP,IEL) = dweight1*du2locy

            end do ! ICUBP

          end do ! IEL

        end if ! dnewton != 0

      else
!        print *,"dweight2 .ne. 0"
        do IEL=1,IELmax-IELset+1

          ! Loop over all cubature points on the current element
          do ICUBP = 1, ncubp

            du1loc = 0.0_DP
            du2loc = 0.0_DP

            ! Perform a loop through the trial DOF`s.
            do JDOFE=1,indof

              ! Get the value of the (trial) basis function
              ! phi_i in the cubature point:
              db = Dbas(JDOFE,1,ICUBP,IEL)

              ! Sum up to the value in the cubature point
              JDFG = Idofs(JDOFE,IEL)
              du1loc = du1loc + (dweight1*u1Xvel(JDFG) + dweight2*u2Xvel(JDFG))*db
              du2loc = du2loc + (dweight1*u1Yvel(JDFG) + dweight2*u2Yvel(JDFG))*db

            end do ! JDOFE

            ! Save the computed velocity
            Dvelocity(1,ICUBP,IEL) = du1loc
            Dvelocity(2,ICUBP,IEL) = du2loc

          end do ! ICUBP

        end do ! IEL

        ! Compute X- and Y-derivative of the velocity?
        if ((dnewton .ne. 0.0_DP) .or. (dnewtonTransposed .ne. 0.0_DP)) then

          do IEL=1,IELmax-IELset+1

            ! Loop over all cubature points on the current element
            do ICUBP = 1, ncubp

              du1locx = 0.0_DP
              du1locy = 0.0_DP
              du2locx = 0.0_DP
              du2locy = 0.0_DP

              ! Perform a loop through the trial DOF`s.
              do JDOFE=1,indof

                ! Get the value of the (trial) basis function
                ! phi_i in the cubature point:
                dbx = Dbas(JDOFE,DER_DERIV_X,ICUBP,IEL)
                dby = Dbas(JDOFE,DER_DERIV_Y,ICUBP,IEL)

                ! Sum up to the value in the cubature point
                JDFG = Idofs(JDOFE,IEL)
                du1locx = du1locx + (dweight1*u1Xvel(JDFG) + dweight2*u2Xvel(JDFG))*dbx
                du1locy = du1locy + (dweight1*u1Xvel(JDFG) + dweight2*u2Xvel(JDFG))*dby
                du2locx = du2locx + (dweight1*u1Yvel(JDFG) + dweight2*u2Yvel(JDFG))*dbx
                du2locy = du2locy + (dweight1*u1Yvel(JDFG) + dweight2*u2Yvel(JDFG))*dby

              end do ! JDOFE

              ! Save the computed velocity derivative
              DvelocityUderiv(1,ICUBP,IEL) = du1locx
              DvelocityUderiv(2,ICUBP,IEL) = du1locy
              DvelocityVderiv(1,ICUBP,IEL) = du2locx
              DvelocityVderiv(2,ICUBP,IEL) = du2locy

            end do ! ICUBP

          end do ! IEL

        end if ! dnewton != 0

      end if

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

      if (bALE) then

        ! Calculate the values of the basis functions in all the points
        ! on all the elements
        call elem_generic_sim2 (EL_Q1, &
            revalElementSet, Bder, DbasALE)

        ! Loop over all elements in the current set
        do IEL=1,IELmax-IELset+1

          ! Loop over all cubature points on the current element
          do ICUBP = 1, ncubp

            du1loc = 0.0_DP
            du2loc = 0.0_DP

            ! Perform a loop through the trial DOF`s.
            do JDOFE=1,indof

              ! Get the value of the (trial) basis function
              db= Dbas(JDOFE,1,ICUBP,IEL)

              ! Sum up to the value in the cubature point
              JDFG = IdofsALE(IDOFE,IEL)
              du1loc = du1loc + DmeshVelocity(1,JDFG)*db
              du2loc = du2loc + DmeshVelocity(2,JDFG)*db

            end do ! JDOFE

            ! Save the computed velocity
            Dvelocity(1,ICUBP,IEL) = Dvelocity(1,ICUBP,IEL) - du1loc
            Dvelocity(2,ICUBP,IEL) = Dvelocity(2,ICUBP,IEL) - du2loc

          end do ! ICUBP

        end do ! IEL

        ! Subtract the X- and Y-derivative of the mesh velocity to the
        ! velocity derivative field if Newton is active.
        if ((dnewton .ne. 0.0_DP) .or. (dnewtonTransposed .ne. 0.0_DP)) then

          do IEL=1,IELmax-IELset+1

            ! Loop over all cubature points on the current element
            do ICUBP = 1, ncubp

              du1locx = 0.0_DP
              du1locy = 0.0_DP
              du2locx = 0.0_DP
              du2locy = 0.0_DP

              ! Perform a loop through the trial DOF`s.
              do JDOFE=1,indof

                ! Get the value of the (trial) basis function
                ! phi_i in the cubature point:
                dbx = Dbas(JDOFE,DER_DERIV_X,ICUBP,IEL)
                dby = Dbas(JDOFE,DER_DERIV_Y,ICUBP,IEL)

                ! Sum up to the value in the cubature point
                JDFG = Idofs(JDOFE,IEL)
                du1locx = du1locx + DmeshVelocity(1,JDFG)*dbx
                du1locy = du1locy + DmeshVelocity(1,JDFG)*dby
                du2locx = du2locx + DmeshVelocity(2,JDFG)*dbx
                du2locy = du2locy + DmeshVelocity(2,JDFG)*dby

              end do ! JDOFE

              ! Subtract the velocity derivative to the previously calculated one.
              DvelocityUderiv(1,ICUBP,IEL) = DvelocityUderiv(1,ICUBP,IEL)-du1locx
              DvelocityUderiv(2,ICUBP,IEL) = DvelocityUderiv(2,ICUBP,IEL)-du1locy
              DvelocityVderiv(1,ICUBP,IEL) = DvelocityVderiv(1,ICUBP,IEL)-du2locx
              DvelocityVderiv(2,ICUBP,IEL) = DvelocityVderiv(2,ICUBP,IEL)-du2locy

            end do ! ICUBP

          end do ! IEL

        end if ! dnewton != 0

      end if

      ! Ok, we now use Dvelocity as coefficient array in the assembly
      ! of a bilinear form!
      !
      ! Clear the local matrices. If the Newton part is to be calculated,
      ! we must clear everything, otherwise only Dentry.
      if ((dnewton .ne. 0.0_DP) .or. (dnewtonTransposed .ne. 0.0_DP) &
          .or. (ddeltaTransposed .ne. 0.0_DP)) then
        Dentry = 0.0_DP
        DentryA11 = 0.0_DP
        DentryA12 = 0.0_DP
        DentryA21 = 0.0_DP
        DentryA22 = 0.0_DP
      else
        Dentry = 0.0_DP
      end if

      ! If ddelta != 0, set up the nonlinearity U*grad(u), probably with
      ! streamline diffusion stabilisation.
      if (ddelta .ne. 0.0_DP) then

        ! Loop over the elements in the current set.
        do IEL=1,IELmax-IELset+1

          ! Nonlinearity:
          !    ddelta * u_1 * grad(.)
          !  = ddelta * [ (DU1) (dx)            ]
          !             [            (DU2) (dy) ]
          !
          ! Loop over all cubature points on the current element
          do ICUBP = 1, ncubp

            ! Calculate the current weighting factor in the cubature formula
            ! in that cubature point.
            !
            ! Normally, we have to take the absolut value of the determinant
            ! of the mapping here!
            ! In 2D, the determinant is always positive, whereas in 3D,
            ! the determinant might be negative -- that is normal!
            ! But because this routine only works in 2D, we can skip
            ! the ABS here!

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
            ! Outer loop over the DOF`s i=1..indof on our current element,
            ! which corresponds to the basis functions Phi_i:

            do IDOFE=1,indof

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

              ! Inner loop over the DOF`s j=1..indof, which corresponds to
              ! the basis function Phi_j:

              do JDOFE=1,indof

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
                ! we do not have to worry about that.

                HSUMJ = HBASJ2*du1loc+HBASJ3*du2loc

                ! Finally calculate the contribution to the system
                ! matrix. Depending on the configuration of ddelta,...
                ! this is:
                !
                ! AH = n~_h(u_h,phi_j,phi_i)        | nonlinear part
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

                AH = ddelta * HSUMJ*(DlocalDelta(IEL)*HSUMI+HBASI1)

                ! Weighten the calculated value AH by the cubature
                ! weight OM and add it to the local matrix. After the
                ! loop over all DOF`s is finished, each entry contains
                ! the calculated integral.

                Dentry(JDOFE,IDOFE,IEL) = Dentry(JDOFE,IDOFE,IEL)+OM*AH

              end do ! IDOFE

            end do ! JDOFE

          end do ! ICUBP

        end do ! IEL

      end if



      ! If dny != 0 or dalpha != 0, add the Laplace/Mass matrix to the
      ! local matrices.
      if ((dalpha .ne. 0.0_DP) .or. (dny .ne. 0.0_DP)) then

        ! Loop over the elements in the current set.
        do IEL=1,IELmax-IELset+1

          ! Loop over all cubature points on the current element
          do ICUBP = 1, ncubp

            ! Calculate the current weighting factor in the cubature formula
            ! in that cubature point.
            !
            ! Normally, we have to take the absolut value of the determinant
            ! of the mapping here!
            ! In 2D, the determinant is always positive, whereas in 3D,
            ! the determinant might be negative -- that is normal!
            ! But because this routine only works in 2D, we can skip
            ! the ABS here!

            OM = Domega(ICUBP)*p_Ddetj(ICUBP,IEL)

            ! Current velocity in this cubature point:
            du1loc = Dvelocity (1,ICUBP,IEL)
            du2loc = Dvelocity (2,ICUBP,IEL)

            ! Outer loop over the DOF`s i=1..indof on our current element,
            ! which corresponds to the basis functions Phi_i:

            do IDOFE=1,indof

              ! Fetch the contributions of the (test) basis functions Phi_i
              ! (our "O")  for function value and first derivatives for the
              ! current DOF into HBASIy:

              HBASI1 = Dbas(IDOFE,1,ICUBP,IEL)
              HBASI2 = Dbas(IDOFE,2,ICUBP,IEL)
              HBASI3 = Dbas(IDOFE,3,ICUBP,IEL)

              ! Inner loop over the DOF`s j=1..indof, which corresponds to
              ! the basis function Phi_j:

              do JDOFE=1,indof

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! (out "X") for function value and first derivatives for the
                ! current DOF into HBASJy:

                HBASJ1 = Dbas(JDOFE,1,ICUBP,IEL)
                HBASJ2 = Dbas(JDOFE,2,ICUBP,IEL)
                HBASJ3 = Dbas(JDOFE,3,ICUBP,IEL)

                ! Finally calculate the contribution to the system
                ! matrix. Depending on the configuration of DNU,
                ! dalpha,... this decomposes into:
                !
                ! AH = dny*(grad(phi_j,grad(phi_i)) | -dny*Laplace(u) = -dbeta*Stokes
                !    + dalpha*(phi_j*phi_i)         | Mass matrix

                AH = dny*(HBASI2*HBASJ2+HBASI3*HBASJ3) &
                    + dalpha*HBASI1*HBASJ1

                ! Weighten the calculated value AH by the cubature
                ! weight OM and add it to the local matrix. After the
                ! loop over all DOF`s is finished, each entry contains
                ! the calculated integral.

                Dentry(JDOFE,IDOFE,IEL) = Dentry(JDOFE,IDOFE,IEL)+OM*AH

              end do ! IDOFE

            end do ! JDOFE

          end do ! ICUBP

        end do ! IEL

      end if

      ! Should we assemble the Newton matrices?
      if (dnewton .ne. 0.0_DP) then

        ! Newton operator
        !
        !    dnewton * (.) * grad(u_1)
        !  = dnewton * [ (dx DU1) (dy DU1) ]
        !              [ (dx DU2) (dy DU2) ]
        !
        ! Loop over the elements in the current set.
        do IEL=1,IELmax-IELset+1

          ! Loop over all cubature points on the current element
          do ICUBP = 1, ncubp

            ! Calculate the current weighting factor in the cubature formula
            ! in that cubature point.
            !
            ! Normally, we have to take the absolut value of the determinant
            ! of the mapping here!
            ! In 2D, the determinant is always positive, whereas in 3D,
            ! the determinant might be negative -- that is normal!
            ! But because this routine only works in 2D, we can skip
            ! the ABS here!

            OM = Domega(ICUBP)*p_Ddetj(ICUBP,IEL)

            ! Current velocity in this cubature point:
            du1locx = DvelocityUderiv (1,ICUBP,IEL)
            du1locy = DvelocityUderiv (2,ICUBP,IEL)
            du2locx = DvelocityVderiv (1,ICUBP,IEL)
            du2locy = DvelocityVderiv (2,ICUBP,IEL)

            ! Outer loop over the DOF`s i=1..indof on our current element,
            ! which corresponds to the basis functions Phi_i:

            do IDOFE=1,indof

              ! Fetch the contributions of the (test) basis functions Phi_i
              ! (our "O")  for function value and first derivatives for the
              ! current DOF into HBASIy:

              HBASI1 = Dbas(IDOFE,1,ICUBP,IEL)

              ! Inner loop over the DOF`s j=1..indof, which corresponds to
              ! the basis function Phi_j:

              do JDOFE=1,indof

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! (out "X") for function value and first derivatives for the
                ! current DOF into HBASJy:

                HBASJ1 = Dbas(JDOFE,1,ICUBP,IEL)

                ! Finally calculate the contribution to the system
                ! matrices A11, A12, A21 and A22.
                !
                ! The Newton part is calculated as follows:
                !
                ! U * grad(V)  =  ( U * grad(.) ) V
                !
                !              =  ( U * grad(V1) )
                !                 ( U * grad(V2) )
                !
                !              =  ( (U1) * (V1x) )
                !                 ( (U2)   (V1y) )
                !                 (              )
                !                 ( (U1) * (V2x) )
                !                 ( (U2)   (V2y) )
                !
                !              =  ( U1 * V1x  + U2 * V1y )
                !                 ( U1 * V2x  + U2 * V2y )
                !
                !              =  ( V1_x  V1_y ) ( U1 )
                !                 ( V2_x  V2_y ) ( U2 )
                !
                !              -> ( A11   A12  )
                !                 ( A21   A22  )
                !
                ! With the velocity V=(u,v), we have to assemble:
                ! grad(V)*U, which is realised in each cubature point as:
                !   du/dx * phi_j*phi_i -> A11
                !   du/dy * phi_j*phi_i -> A12
                !   dv/dx * phi_j*phi_i -> A21
                !   dv/dy * phi_j*phi_i -> A22

                AH11 = dnewton * du1locx * HBASJ1*HBASI1
                AH12 = dnewton * du1locy * HBASJ1*HBASI1
                AH21 = dnewton * du2locx * HBASJ1*HBASI1
                AH22 = dnewton * du2locy * HBASJ1*HBASI1

                ! Weighten the calculated value AHxy by the cubature
                ! weight OM and add it to the local matrices. After the
                ! loop over all DOF`s is finished, each entry contains
                ! the calculated integral.

                DentryA11(JDOFE,IDOFE,IEL) = DentryA11(JDOFE,IDOFE,IEL)+OM*AH11
                DentryA12(JDOFE,IDOFE,IEL) = DentryA12(JDOFE,IDOFE,IEL)+OM*AH12
                DentryA21(JDOFE,IDOFE,IEL) = DentryA21(JDOFE,IDOFE,IEL)+OM*AH21
                DentryA22(JDOFE,IDOFE,IEL) = DentryA22(JDOFE,IDOFE,IEL)+OM*AH22

              end do ! IDOFE

            end do ! JDOFE

          end do ! ICUBP

        end do ! IEL

      end if

      ! Transposed convection operator
      !
      !  ddeltaTransposed * grad(.)^T * u_1
      !  = ddeltaTransposed * [ (DU1+DU2) dx               ]
      !                       [               (DU1+DU2) dy ]
      !
      ! Should we assemble the transposed convection matrices?
      if (ddeltaTransposed .ne. 0.0_DP) then

        ! Loop over the elements in the current set.
        do IEL=1,IELmax-IELset+1

          ! Loop over all cubature points on the current element
          do ICUBP = 1, ncubp

            ! Calculate the current weighting factor in the cubature formula
            ! in that cubature point.
            !
            ! Normally, we have to take the absolut value of the determinant
            ! of the mapping here!
            ! In 2D, the determinant is always positive, whereas in 3D,
            ! the determinant might be negative -- that is normal!
            ! But because this routine only works in 2D, we can skip
            ! the ABS here!

            OM = Domega(ICUBP)*p_Ddetj(ICUBP,IEL)

            ! Current velocity in this cubature point:
            du1loc = Dvelocity (1,ICUBP,IEL)
            du2loc = Dvelocity (2,ICUBP,IEL)

            ! Outer loop over the DOF`s i=1..indof on our current element,
            ! which corresponds to the basis functions Phi_i:

            do IDOFE=1,indof

              ! Fetch the contributions of the (test) basis functions Phi_i
              ! (our "O")  for function value and first derivatives for the
              ! current DOF into HBASIy:

              HBASI1 = Dbas(IDOFE,1,ICUBP,IEL)
              HBASI2 = Dbas(IDOFE,2,ICUBP,IEL)
              HBASI3 = Dbas(IDOFE,3,ICUBP,IEL)

              ! Inner loop over the DOF`s j=1..indof, which corresponds to
              ! the basis function Phi_j:

              do JDOFE=1,indof

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! (out "X") for function value and first derivatives for the
                ! current DOF into HBASJy:

                HBASJ1 = Dbas(JDOFE,1,ICUBP,IEL)
                HBASJ2 = Dbas(JDOFE,2,ICUBP,IEL)
                HBASJ3 = Dbas(JDOFE,3,ICUBP,IEL)

                ! Finally calculate the contribution to the system
                ! matrices A11, A12, A21 and A22.

                AH11 = ddeltaTransposed * du1loc * HBASJ2*HBASI1
                AH12 = ddeltaTransposed * du2loc * HBASJ2*HBASI1
                AH21 = ddeltaTransposed * du1loc * HBASJ3*HBASI1
                AH22 = ddeltaTransposed * du2loc * HBASJ3*HBASI1

                ! Weighten the calculated value AHxy by the cubature
                ! weight OM and add it to the local matrices. After the
                ! loop over all DOF`s is finished, each entry contains
                ! the calculated integral.

                DentryA11(JDOFE,IDOFE,IEL) = DentryA11(JDOFE,IDOFE,IEL)+OM*AH11
                DentryA12(JDOFE,IDOFE,IEL) = DentryA12(JDOFE,IDOFE,IEL)+OM*AH12
                DentryA21(JDOFE,IDOFE,IEL) = DentryA21(JDOFE,IDOFE,IEL)+OM*AH21
                DentryA22(JDOFE,IDOFE,IEL) = DentryA22(JDOFE,IDOFE,IEL)+OM*AH22

              end do ! IDOFE

            end do ! JDOFE

          end do ! ICUBP

        end do ! IEL

      end if

      ! Transposed Newton operator
      !
      ! <tex> $$ dnewtonTransposed * grad(u_1)^T * (.) ) $$ </tex>
      !        = dnewtonTransposed * [ (dx DU1) (dx DU2) ]
      !                              [ (dy DU1) (dy DU2) ]
      !
      ! Should we assemble the transposed Newton matrices?
      if (dnewtonTransposed .ne. 0.0_DP) then

        ! Loop over the elements in the current set.
        do IEL=1,IELmax-IELset+1

          ! Loop over all cubature points on the current element
          do ICUBP = 1, ncubp

            ! Calculate the current weighting factor in the cubature formula
            ! in that cubature point.
            !
            ! Normally, we have to take the absolut value of the determinant
            ! of the mapping here!
            ! In 2D, the determinant is always positive, whereas in 3D,
            ! the determinant might be negative -- that is normal!
            ! But because this routine only works in 2D, we can skip
            ! the ABS here!

            OM = Domega(ICUBP)*p_Ddetj(ICUBP,IEL)

            ! Current velocity in this cubature point:
            du1locx = DvelocityUderiv (1,ICUBP,IEL)
            du1locy = DvelocityUderiv (2,ICUBP,IEL)
            du2locx = DvelocityVderiv (1,ICUBP,IEL)
            du2locy = DvelocityVderiv (2,ICUBP,IEL)

            ! Outer loop over the DOF`s i=1..indof on our current element,
            ! which corresponds to the basis functions Phi_i:

            do IDOFE=1,indof

              ! Fetch the contributions of the (test) basis functions Phi_i
              ! (our "O")  for function value and first derivatives for the
              ! current DOF into HBASIy:

              HBASI1 = Dbas(IDOFE,1,ICUBP,IEL)

              ! Inner loop over the DOF`s j=1..indof, which corresponds to
              ! the basis function Phi_j:

              do JDOFE=1,indof

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! (out "X") for function value and first derivatives for the
                ! current DOF into HBASJy:

                HBASJ1 = Dbas(JDOFE,1,ICUBP,IEL)

                ! Finally calculate the contribution to the system
                ! matrices A11, A12, A21 and A22.
                !
                ! With the velocity V=(u,v), we have to assemble:
                ! grad(V)^t*U, which is realised in each cubature point as:
                !   du/dx * phi_j*phi_i -> A11
                !   dv/dx * phi_j*phi_i -> A12
                !   du/dy * phi_j*phi_i -> A21
                !   dv/dy * phi_j*phi_i -> A22

                AH11 = dnewtonTransposed * du1locx * HBASJ1*HBASI1
                AH12 = dnewtonTransposed * du2locx * HBASJ1*HBASI1
                AH21 = dnewtonTransposed * du1locy * HBASJ1*HBASI1
                AH22 = dnewtonTransposed * du2locy * HBASJ1*HBASI1

                ! Weighten the calculated value AHxy by the cubature
                ! weight OM and add it to the local matrices. After the
                ! loop over all DOF`s is finished, each entry contains
                ! the calculated integral.

                DentryA11(JDOFE,IDOFE,IEL) = DentryA11(JDOFE,IDOFE,IEL)+OM*AH11
                DentryA12(JDOFE,IDOFE,IEL) = DentryA12(JDOFE,IDOFE,IEL)+OM*AH12
                DentryA21(JDOFE,IDOFE,IEL) = DentryA21(JDOFE,IDOFE,IEL)+OM*AH21
                DentryA22(JDOFE,IDOFE,IEL) = DentryA22(JDOFE,IDOFE,IEL)+OM*AH22

              end do ! IDOFE

            end do ! JDOFE

          end do ! ICUBP

        end do ! IEL

      end if

      ! Now we have set up "local" system matrices. We can either
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

      if (iand(cdef,CONV_MODMATRIX) .ne. 0) then

        ! With or without Newton?
        if ((dnewton .eq. 0.0_DP) .and. (dnewtonTransposed .eq. 0.0_DP)) then

          ! Include the local matrices into the global system matrix,
          ! subblock A11 and (if different from A11) also into A22.
          !%omp critical
          do IEL=1,IELmax-IELset+1
            do IDOFE=1,indof
              do JDOFE=1,indof
                p_Da11(Kentry(JDOFE,IDOFE,IEL)) = p_Da11(Kentry(JDOFE,IDOFE,IEL)) + &
                    dtheta * Dentry(JDOFE,IDOFE,IEL)
              end do
            end do
          end do
          !%omp end critical

          if (.not. associated(p_Da11,p_Da22)) then
            !%omp critical
            do IEL=1,IELmax-IELset+1
              do IDOFE=1,indof
                do JDOFE=1,indof
                  p_Da22(Kentry(JDOFE,IDOFE,IEL)) = &
                      p_Da22(Kentry(JDOFE,IDOFE,IEL)) + &
                      dtheta * Dentry(JDOFE,IDOFE,IEL)
                end do
              end do
            end do
            !%omp end critical

          end if

        else

          ! Include the local matrices into the global system matrix,
          ! subblock A11 and A22 (both must exist and be independent from
          ! each other).
          !%omp critical
          do IEL=1,IELmax-IELset+1
            do IDOFE=1,indof
              do JDOFE=1,indof
                ! Kentry (:,:,:) -> positions of local matrix in A11 and A22.
                !
                ! DentryA11 (:,:,:) -> Newton part of A11
                p_Da11(Kentry(JDOFE,IDOFE,IEL)) = p_Da11(Kentry(JDOFE,IDOFE,IEL)) + &
                    dtheta * ( Dentry(JDOFE,IDOFE,IEL) + &
                               DentryA11(JDOFE,IDOFE,IEL) )

                ! DentryA22 (:,:,:) -> Newton part of A22
                p_Da22(Kentry(JDOFE,IDOFE,IEL)) = p_Da22(Kentry(JDOFE,IDOFE,IEL)) + &
                    dtheta * ( Dentry(JDOFE,IDOFE,IEL) + &
                               DentryA22(JDOFE,IDOFE,IEL) )
              end do
            end do
          end do
          !%omp end critical

          !%omp critical
          ! Include the local Newton matrix parts into A12 and A21.
          do IEL=1,IELmax-IELset+1
            do IDOFE=1,indof
              do JDOFE=1,indof
                ! Kentry12 (:,:,:) -> positions of local matrix in A12 and A21.
                !
                ! Dentry (:,:,:) -> Newton part of A12
                p_Da12(Kentry12(JDOFE,IDOFE,IEL)) = p_Da12(Kentry12(JDOFE,IDOFE,IEL)) + &
                    dtheta * DentryA12(JDOFE,IDOFE,IEL)

                ! Dentry (:,:,:) -> Newton part of A21
                p_Da21(Kentry12(JDOFE,IDOFE,IEL)) = p_Da21(Kentry12(JDOFE,IDOFE,IEL)) + &
                    dtheta * DentryA21(JDOFE,IDOFE,IEL)
              end do
            end do
          end do
          !%omp end critical

        end if

      end if

      ! For cdef containing CONV_MODDEFECT, build the defect vector
      !     D = RHS - A*U
      ! This is done matrix free, only with the help of the local
      ! matrix.
      ! In this case, D=(D1,D2) is expected to be the RHS on
      ! entry and will be updated to be the defect vector when
      ! this routine is left.

      if (iand(cdef,CONV_MODDEFECT) .ne. 0) then

        ! With or without Newton?
        if ((dnewton .eq. 0.0_DP) .and. (dnewtonTransposed .eq. 0.0_DP)) then
          !%omp critical
          do IEL=1,IELmax-IELset+1
            do IDOFE=1,indof

              IDFG=Idofs(IDOFE,IEL)

              do JDOFE=1,indof

                denth = dtheta*Dentry(JDOFE,IDOFE,IEL)

                JDFG=Idofs(JDOFE,IEL)
                Ddef1(IDFG)= Ddef1(IDFG) - denth*Du1(JDFG)
                Ddef2(IDFG)= Ddef2(IDFG) - denth*Du2(JDFG)

              end do
            end do
          end do
          !%omp end critical
        else
          !%omp critical
          do IEL=1,IELmax-IELset+1
            do IDOFE=1,indof

              IDFG=Idofs(IDOFE,IEL)

              do JDOFE=1,indof

                denth = dtheta*Dentry(JDOFE,IDOFE,IEL)

                JDFG=Idofs(JDOFE,IEL)
                Ddef1(IDFG)= Ddef1(IDFG) - denth*Du1(JDFG)
                Ddef2(IDFG)= Ddef2(IDFG) - denth*Du2(JDFG)

                ! Newton part
                Ddef1(IDFG)= Ddef1(IDFG) &
                           - dtheta*DentryA11(JDOFE,IDOFE,IEL)*Du1(JDFG) &
                           - dtheta*DentryA12(JDOFE,IDOFE,IEL)*Du2(JDFG)
                Ddef2(IDFG)= Ddef2(IDFG) &
                           - dtheta*DentryA21(JDOFE,IDOFE,IEL)*Du1(JDFG) &
                           - dtheta*DentryA22(JDOFE,IDOFE,IEL)*Du2(JDFG)

              end do
            end do
          end do
          !%omp end critical
        end if

      end if


    end do ! IELset
    !%omp end do

    ! Release memory
    call elprep_releaseElementSet(revalElementSet)

    ! Release the assembly structure if necessary.
    if (.not. present(rcubatureInfo)) then
      call spdiscr_releaseCubStructure(rtempCubatureInfo)
    end if

    deallocate(DlocalDelta)
    if ((dnewton .ne. 0.0_DP) .or. (dnewtonTransposed .ne. 0.0_DP)) then
      deallocate(DentryA22)
      deallocate(DentryA21)
      deallocate(DentryA12)
      deallocate(DentryA11)
      deallocate(Kentry12)
      deallocate(DvelocityUderiv)
      deallocate(DvelocityVderiv)
    end if
    deallocate(Dvelocity)
    deallocate(Dentry)
    deallocate(Kentry)
    deallocate(IdofsALE)
    deallocate(Idofs)
    deallocate(DbasALE)
    deallocate(Dbas)
    !%omp end parallel
    deallocate(Domega)
    deallocate(p_DcubPtsRef)

  end subroutine

  ! ----------------------------------------------------------------------

  subroutine getLocalDeltaQuad (clocalh,&
                      Du1x,Du1y,Du2x,Du2y,da1,da2,Ielements,&
                      duMaxR,Ddelta,rtriangulation,Idofs,dupsam,dnurec)

  ! This routine calculates a local ddelta=DELTA_T for a set of finite
  ! elements Ielements. This can be used by the streamline diffusion
  ! stabilisation technique as a multiplier of the (local) bilinear form.
  !
  ! The effective velocity that is used for calculating the ddelta
  ! is combined by a weighted mean of the two velocity fields U1,U2
  ! by:
  !                   du = da1*Du1 + da2*Du2
  ! The coefficients A1,A2 allow the caller to take influence on which
  ! velocity field to weight more.

  ! Method how to compute the local h.
  ! =0: Use the root of the area of the element as local H
  ! =1: Use the length of the way that a particle travels through
  !     the element in direction of the flow
  integer, intent(in) :: clocalH

  ! Main velocity field.
  real(DP), dimension(*), intent(in) :: du1x,du1y

  ! Secondary velocity field.
  real(DP), dimension(*), intent(in) :: du2x,du2y

  ! weighting factor for Du1
  real(DP), intent(in) :: da1

  ! weighting factor for Du2
  real(DP), intent(in) :: da2

  ! Reciprocal of the maximum norm of velocity in the domain:
  ! 1/duMaxR = 1/||u||_Omega
  real(DP), intent(in) :: duMaxR

  ! Reciprocal value 1/NU of coefficient NU in front of the
  ! Laplacian term of the Navier-Stokes equation
  !   NU * Laplace(u) + u*grad(u) + ...
  real(DP), intent(in) :: dnuRec

  ! user defined parameter for configuring the streamline diffusion.
  ! < 0: Simple calculation of ddelta, using
  !      ddelta = |UPSAM| * h_T.
  ! > 0: usually UPSAM = 0.1 .. 2; Samarskji-like calculation of ddelta using:
  !      ddelta = UPSAM * h_t/||u||_T * 2*Re_T/(1+Re_T)
  real(DP), intent(in) :: dupsam

  ! List of elements where the Ddelta should be calculated
  integer, dimension(:), intent(in) :: Ielements

  ! Array with global degrees of freedom on the elements
  integer, dimension(:,:), intent(in) :: Idofs

  ! Triangulation that defines the mesh.
  type(t_triangulation), intent(in) :: rtriangulation

  ! Out: local Ddelta on all elements
  real(DP), dimension(:), intent(out) :: ddelta

  ! local variables
  real(DP) :: dlocalH,du1,du2,dunorm,dreLoc
  integer :: iel,ielidx
  integer :: idof
  integer, dimension(:,:), pointer :: p_IverticesAtElement
  real(DP), dimension(:,:), pointer :: p_DvertexCoords
  real(DP), dimension(:), pointer :: p_DelementVolume

    ! Get some crucial data
    if (clocalh .eq. 0) then
      call storage_getbase_double (rtriangulation%h_DelementVolume,p_DelementVolume)

      ! Loop through all elements
      do ielidx = 1,size(Ielements)

        iel = Ielements(ielidx)

        ! Loop through the local degrees of freedom on element IEL.
        ! Sum up the velocities on these DOF`s. This will result
        ! in the vector (DU1,DU2) representing the (mean) X/Y-velocity
        ! through element IEL.

        ! For elements whose DOF`s represent directly the velocity, U1/U2
        ! represent the mean velocity
        ! along an egde/on the midpoint of each edge, so U1/U2 is
        ! clearly an approximation to the velocity in element T.

        du1=0.0_DP
        du2=0.0_DP
        do idof=1,ubound(Idofs,1)
          du1=du1+(da1*du1x(Idofs(idof,ielidx))+da2*du2x(Idofs(idof,ielidx)))
          du2=du2+(da1*du1y(Idofs(idof,ielidx))+da2*du2y(Idofs(idof,ielidx)))
        end do

        ! Calculate the norm of that local velocity:

        dunorm = sqrt(du1**2+du2**2) / real(ubound(Idofs,1),DP)

        ! Now we have:   dunorm = ||u||_T
        ! and:           u_T = a1*u1_T + a2*u2_T

        ! If the norm of the velocity is small, we choose ddelta = 0,
        ! which results in central difference in the streamline diffusion
        ! matrix assembling:

        if (dunorm .le. 1E-8_DP) then

          Ddelta(ielidx) = 0.0_DP

        else

          ! Calculate the local h from the area of the element
          dlocalH = sqrt(p_DelementVolume(iel))

          ! Calculate ddelta... (cf. p. 121 in Turek`s CFD book)

          if (dupsam .lt. 0.0_DP) then

            ! For UPSAM<0, we use simple calculation of ddelta:

            Ddelta(ielidx) = abs(dupsam)*dlocalH

          else

            ! For UPSAM >= 0, we use standard Samarskji-like calculation
            ! of ddelta. At first calculate the local Reynolds number
            ! RELOC = Re_T = ||u||_T * h_T / NU

            dreLoc = dunorm*dlocalH*dnuRec

            ! and then the ddelta = UPSAM * h_t/||u|| * 2*Re_T/(1+Re_T)

            Ddelta(ielidx) = dupsam * dlocalH*duMaxR * 2.0_DP*(dreLoc/(1.0_DP+dreLoc))

          end if ! (UPSAM.LT.0.0)

        end if ! (dunorm.LE.1D-8)

      end do

    else

      call storage_getbase_double2d (rtriangulation%h_DvertexCoords,p_DvertexCoords)
      call storage_getbase_int2d (rtriangulation%h_IverticesAtElement,p_IverticesAtElement)

      ! Loop through all elements
      do ielidx = 1,size(Ielements)

        iel = Ielements(ielidx)

        ! Loop through the local degrees of freedom on element IEL.
        ! Sum up the velocities on these DOF`s. This will result
        ! in the vector (DU1,DU2) representing the (mean) X/Y-velocity
        ! through element IEL.

        ! For elements whose DOF`s represent directly the velocity, U1/U2
        ! represent the mean velocity
        ! along an egde/on the midpoint of each edge, so U1/U2 is
        ! clearly an approximation to the velocity in element T.

        du1=0.0_DP
        du2=0.0_DP
        do idof=1,ubound(Idofs,1)
          du1=du1+(da1*du1x(Idofs(idof,ielidx))+da2*du2x(Idofs(idof,ielidx)))
          du2=du2+(da1*du1y(Idofs(idof,ielidx))+da2*du2y(Idofs(idof,ielidx)))
        end do

        ! Calculate the norm of the local velocity and the velocity itself.
        du1 = du1 / real(ubound(Idofs,1),DP)
        du2 = du2 / real(ubound(Idofs,1),DP)
        dunorm = sqrt(du1**2+du2**2)

        ! Now we have:   dunorm = ||u||_T
        ! and:           u_T = a1*u1_T + a2*u2_T

        ! If the norm of the velocity is small, we choose ddelta = 0,
        ! which results in central difference in the streamline diffusion
        ! matrix assembling:

        if (dunorm .le. 1E-8_DP) then

          Ddelta(ielidx) = 0.0_DP

        else

          ! u_T defines the "slope" of the velocity through
          ! the element T. At next, calculate the local mesh width
          ! dlocalH = h = h_T on our element T=IEL:

          call getLocalMeshWidthQuad (dlocalH,dunorm, du1, du2, iel, &
              p_IverticesAtElement,p_DvertexCoords)

          ! Calculate ddelta... (cf. p. 121 in Turek`s CFD book)

          if (dupsam .lt. 0.0_DP) then

            ! For UPSAM<0, we use simple calculation of ddelta:

            Ddelta(ielidx) = abs(dupsam)*dlocalH

          else

            ! For UPSAM >= 0, we use standard Samarskji-like calculation
            ! of ddelta. At first calculate the local Reynolds number
            ! RELOC = Re_T = ||u||_T * h_T / NU

            dreLoc = dunorm*dlocalH*dnuRec

            ! and then the ddelta = UPSAM * h_t/||u|| * 2*Re_T/(1+Re_T)

            Ddelta(ielidx) = dupsam * dlocalH*duMaxR * 2.0_DP*(dreLoc/(1.0_DP+dreLoc))

          end if ! (UPSAM.LT.0.0)

        end if ! (dunorm.LE.1D-8)

      end do

    end if

  end subroutine


  ! ----------------------------------------------------------------------


!<subroutine>
  pure subroutine getLocalMeshWidthQuad (dlocalH, dunorm,  XBETA1, &
                      XBETA2, JEL,Kvert,Dcorvg)
  !<description>
    ! This routine determines the length of the maximum path a particle can travel in
    ! element JEL of a triangulation in the direction given by vector beta.
    ! This length is one possible definition of a local mesh width.
  !</description>

  !<input>
    ! Element where the local h should be calculated
    integer, intent(in) :: JEL

    integer, dimension(TRIA_MAXNVE2D,*), intent(in) :: Kvert
    real(DP), dimension(NDIM2D,*), intent(in) :: Dcorvg

    ! norm ||u||_T = mean velocity through element T=JEL
    real(DP), intent(in) :: dunorm

    ! mean velocity u_T = (xbeta1,xbeta2) through element T=JEL
    real(DP), intent(in) :: XBETA1, XBETA2
  !</input>

  !<output>
    ! local mesh width
    real(DP), intent(out) :: dlocalH
  !</output>
!</subroutine>

    ! local variables
    real(DP) :: dlambda
    integer :: NECK1,NECK2,NECK3,NECK4
    real(DP) :: X1,Y1,X2,Y2,X3,Y3,X4,Y4
    real(DP) :: dalphaMax, dalpha

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

    ! In the next step, we calculate the `maximum possible mesh width
    ! in direction of the flow`; this is the maximum possible length
    ! that a particle can cross in the current element.
    ! The picture in mind is the following:
    !
    !          G3
    !   +-------------X-------+
    !   |            /        |
    !   |           /         |
    !   |          /          |
    !   |         /           |
    !   |        /            |
    ! G4|       /             | G2
    !   |      ^ (beta1,beta2)|
    !   |     /               |
    !   |    /                |
    !   |   /                 |
    !   |  /                  |
    !   | /                   |
    !   |/                    |
    !   O---------------------+
    !            G1
    !
    ! The vector (beta1,beta2) indicates the direction of the flow.
    ! A particle starting in point O moves at most up to point X.
    ! The length of the line (O,X) is defined to be the local mesh width h.
    !
    ! Loop through the four corners of element JEL and check for a line
    ! with slope BETA=(xbeta1,xbeta2) starting in this corner whether it
    ! really intersects with one of the edges of the element. Remark
    ! that we only have to check the two opposite edges to the current
    ! corner!

    ! -----------------------------------------------------------------
    ! Check the first corner:

    call intersectLines2D(X1,Y1,dalpha,XBETA1,XBETA2, &
         X3,Y3,dlambda,X2,Y2)
    dalphaMax=max(abs(dalpha),dalphaMax)

    call intersectLines2D(X1,Y1,dalpha,XBETA1,XBETA2, &
         X3,Y3,dlambda,X4,Y4)
    dalphaMax=max(abs(dalpha),dalphaMax)

    ! -----------------------------------------------------------------
    ! The second one...

    call intersectLines2D(X2,Y2,dalpha,XBETA1,XBETA2, &
         X4,Y4,dlambda,X1,Y1)
    dalphaMax=max(abs(dalpha),dalphaMax)

    call intersectLines2D(X2,Y2,dalpha,XBETA1,XBETA2, &
         X4,Y4,dlambda,X3,Y3)
    dalphaMax=max(abs(dalpha),dalphaMax)

    ! -----------------------------------------------------------------
    ! The third one...

    call intersectLines2D(X3,Y3,dalpha,XBETA1,XBETA2, &
         X1,Y1,dlambda,X2,Y2)
    dalphaMax=max(abs(dalpha),dalphaMax)

    call intersectLines2D(X3,Y3,dalpha,XBETA1,XBETA2, &
         X1,Y1,dlambda,X4,Y4)
    dalphaMax=max(abs(dalpha),dalphaMax)

    ! -----------------------------------------------------------------
    ! And the fourth=last one...

    call intersectLines2D(X4,Y4,dalpha,XBETA1,XBETA2, &
         X2,Y2,dlambda,X1,Y1)
    dalphaMax=max(abs(dalpha),dalphaMax)

    call intersectLines2D(X4,Y4,dalpha,XBETA1,XBETA2, &
         X2,Y2,dlambda,X3,Y3)
    dalphaMax=max(abs(dalpha),dalphaMax)

    ! -----------------------------------------------------------------
    ! finally determine the local h=h_T
    !
    ! dalphaMax is the stretching/shrink factor for the vector
    ! (dbeta1,dbeta2) to cover the longest distance in the current quad.
    ! We multiply with dunorm=|(dbeta1,dbeta2)| to get the actual length
    ! of the vector which can be placed inside of the element.
    dlocalH = dalphaMax * dunorm

  end subroutine getLocalMeshWidthQuad


  ! ----------------------------------------------------------------------


!<subroutine>
  pure subroutine intersectLines2D (xo, yo, dalpha, beta1, beta2, &
       xa, ya, dlambda, xb, yb)
  !<description>
    ! This routine determines whether two lines intersect in <tex>$R^2$</tex>.
    ! If they do, the parameters along the two lines for the intersection point
    ! are returned. Zero is returned, if the two lines are parallel or coincide.
  !</description>

  !<input>
    ! Origin of line 1
    real(DP), intent(in) :: xo, yo

    ! Direction of line 1
    real(DP), intent(in) :: beta1, beta2

    ! Start point of segment defining the second line
    real(DP), intent(in) :: xa,ya

    ! End point of segment defining the second line
    real(DP), intent(in) :: xb,yb
  !</input>

  !<output>
    ! Parameter value of the intersection point on line 1.
    ! = 0.0, if there is no intersection point
    real(DP), intent(out) :: dalpha

    ! Parameter value of the intersection point on line 2.
    ! Has to be between 0 and 1, otherwise there is no intersection point within the given
    ! segment
    real(DP), intent(out) :: dlambda
  !</output>
!</subroutine>


    ! local variables
    ! Determinant of equation system, denominator when applying Cramer`s rule to
    ! solve equation system to determine intersection point between two lines
    real(DP) :: denominator

    ! tolerance to apply when checking whether 0 <= lambda <= 1
    real(DP), parameter :: TOL = 0.1_DP


    ! Be n1 the clockwise normal of (beta1,beta2). Then holds: n1 = (beta2, -beta1).
    ! Calculate the scalar product of the line (xa,ya)->(xb,yb) with n1:
    denominator = beta2 * (xb-xa) - beta1 * (yb-ya)

    if (denominator .eq. 0.0_DP) then

      ! beta and the line (xa,ya)->(xb,yb) are either parallel disjoint or do coincide
      dalpha = 0.0_DP

    else

      ! Now, solve the following 2D equation system to determine the intersection point of
      ! the two given lines:
      !                  o + \alpha \beta      = a + \lambda (b - a)
      ! <=>   \alpha \beta + \lambda (a - b  ) = a - 0
      ! <=> \alpha \beta_1 + \lambda (xa - xb) = (xa - xo)   (*)
      !     \alpha \beta_2 + \lambda (ya - yb) = (ya - yo)   (**)
      ! Apply Cramer`s rule to determine \alpha and \lambda:
      !              / xa - xo     xa - xb \
      !           det\ ya - yo     ya - yb /
      !  \alpha = --------------------------
      !              / \beta_1     xa - xb \
      !           det\ \beta_2     ya - yb /
      !
      !               / \beta_1     xa - xo \
      !            det\ \beta_2     ya - yo /
      !  \lambda = --------------------------
      !               / \beta_1     xa - xb \
      !            det\ \beta_2     ya - yb /
      dlambda = (beta1 * (ya - yo) + beta2 * (xo - xa)) / denominator

      ! Is the intersection point inside along teh segment given by start- and endpoint?
      if ((dlambda .ge. 0.0_DP - TOL) .and. (dlambda .le. 1.0_DP + TOL)) then
        ! variant 1: determine \alpha also by applying Cramer`s rule
        dalpha = ((xa - xo) * (ya - yb) + (xb - xa) * (ya - yo)) / denominator

        ! variant 2: use \lambda to solve the equation system (*)-(**)
        ! (probably slower because it involves if conditions)
!        if (beta1 .ne. 0.0_DP) then
!          dalpha=((xa-xo)+dlambda*(xb-xa))/beta1
!        else
!          if (beta2 .ne. 0.0_DP) then
!            dalpha=((ya-yo)+dlambda*(yb-ya))/beta2
!          else
!            ! This case cannot happen. It would mean
!            ! BETA1 = BETA2 = 0 and that would have meant denominator = 0,
!            ! a case already catched before.
!            dalpha=0.0_DP
!          end if
!        end if
      else
        dalpha = 0.0_DP
      end if

    end if

  end subroutine intersectLines2D


  ! ----------------------------------------------------------------------


!<subroutine>
  subroutine conv_streamlineDiffusion3d ( &
                           rvecPrimary, rvecSecondary, dprimWeight, dsecWeight,&
                           rconfig, cdef, &
                           rmatrix, rsolution, rdefect, DmeshVelocity, &
                           IvelocityComp, rcubatureInfo, rperfconfig)

!<description>
  ! Standard streamline diffusion method to set up the operator
  ! <tex>
  ! $$ dtheta  *  ( dalpha * MASS  +  dbeta * STOKES  +  ddelta * u_1 * grad(u_2) ) $$
  ! </tex>
  ! in a matrix or to build a defect vector.
  ! 3D-version (X-, Y- and Z-velocity).
  !
  ! rvecPrimary, rvecSecondary are two velocity field vectors for the X-, Y-
  ! and Z-veclocity; IvelocityComp defines which components of these
  ! vectors contains the X-, Y- and which contains the Z-velocity.
  ! The final velocity vector field is then computed as a weighted average
  ! of these two:
  ! <tex>
  !  $$ u_1  =  dprimWeight * rvecPrimary  +  dsecWeight * rvecSecondary $$
  ! </tex>
  ! <tex>$ u_2 = rsolution(.) $</tex> defines the second velocity field inside of
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

  ! Primary velocity field for the computation of <tex>$ u_1 $</tex>
  type(t_vectorBlock), intent(in), target :: rvecPrimary

  ! Secondary velocity field for the computation of <tex>$ u_1 $</tex>
  type(t_vectorBlock), intent(in), target :: rvecSecondary

  ! Weighting factor for rvecPrimary.
  real(DP), intent(in) :: dprimWeight

  ! Weighting factor for rvecSecondary.
  real(DP), intent(in) :: dsecWeight

  ! Configuration block for the streamline diffusion scheme
  type(t_convStreamlineDiffusion), intent(in) :: rconfig

  ! Computation/defect correction method. One of the CONV_MODxxxx constants:
  ! CONV_MODMATRIX: Set up the nonlinear matrix. rmatrix must be present, the
  !                 nonlinear part is added to the matrix.
  ! CONV_MODDEFECT: Set up the nonlinear defect. rdefect and rsolution must be
  !                 present.
  ! CONV_MODBOTH  : Set up the nonlinear matrix as well as the nonlinear defect.
  !                 rmatrix, rdefect and rsolution must all be present.
  integer, intent(in) :: cdef

  ! OPTIONAL: Solution vector u_2.
  ! Must be present if cdef=CONV_MODDEFECT or =CONV_MODBOTH.
  type(t_vectorBlock), intent(in), target, optional :: rsolution

  ! OPTIONAL: Mesh velocity field.
  ! DmeshVelocity(1,ivt) gives the X-velocity of the mesh, i.e. the X-velocity
  !   of the corner vertex ivt.
  ! DmeshVelocity(2,ivt) gives the Y-velocity of the mesh, i.e. the Y-velocity
  !   of the corner vertex ivt.
  ! DmeshVelocity(3,ivt) gives the Z-velocity of the mesh, i.e. the Z-velocity
  !   of the corner vertex ivt.
  ! The parameter must be present if ALE is activated in the
  ! configuration parameter block by bALE=true.
  real(DP), dimension(:,:), intent(in), optional :: DmeshVelocity

  ! OPTIONAL:
  ! Index block that specifies which component in rvecPrimary / rvecSecondary /
  ! rsolution / rdefect is the X-, Y- and Z-velocity.
  !  IvelocityComp(1) gives the number of the X-velocity (usually = 1),
  !  IvelocityComp(2) gives the number of the Y-velocity (usually = 2).
  !  IvelocityComp(3) gives the number of the Z-velocity (usually = 3).
  ! If not present, IvelocityComp=(/1,2,3/) is assumed.
  integer, dimension(3), intent(in), optional :: IvelocityComp

  ! OPTIONAL: A scalar cubature information structure that specifies the cubature
  ! formula(s) to use. If not specified, default settings are used.
  type(t_scalarCubatureInfo), intent(in), optional, target :: rcubatureInfo

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! System matrix.
  ! The content of the matrix must be present if cdef=CONV_MODMATRIX or
  ! =CONV_MODBOTH, otherwise only the structure is used.
  ! The nonlinear operator is added to the matrix.
  type(t_matrixScalar), intent(inout) :: rmatrix

  ! OPTIONAL: Defect vector.
  ! Must have the same structure as rsolution/rvecPrimary/rvecSecondary.
  ! Must be present if cdef=CONV_MODDEFECT or =CONV_MODBOTH.
  ! The nonlinear part is subtracted from this vector:
  ! <tex>$ r = r - \theta * u_1*grad(u_2) $</tex>
  type(t_vectorBlock), intent(inout), optional, target :: rdefect
!</inputoutput>

!</subroutine>

    ! local variables
    integer(I32) :: celement
    integer, dimension(3) :: Icomp
    type(t_vectorScalar), pointer :: p_rvelX1,p_rvelX2,p_rvelY1,p_rvelY2,&
                                     p_rvelZ1,p_rvelZ2
    type(t_vectorScalar), pointer :: p_rsolX,p_rsolY,p_rdefectX,p_rdefectY,&
                                     p_rsolZ,p_rdefectZ
    real(DP), dimension(:), pointer :: p_DvelX1,p_DvelX2,p_DvelY1,p_DvelY2,&
                                       p_DvelZ1,p_DvelZ2
    real(DP), dimension(:), pointer :: p_DsolX,p_DsolY,p_DdefectX,p_DdefectY,&
                                       p_DsolZ,p_DdefectZ

    ! At first check the input parameters that everything is present what
    ! we need:
    if ((cdef .eq. CONV_MODDEFECT) .or. (cdef .eq. CONV_MODBOTH)) then
      if ((.not. present(rsolution)) .or. (.not. present(rdefect))) then
        call output_line ("Solution/defect vector not present!", &
            OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusion3d")
        call sys_halt()
      end if
    end if

    if (rconfig%bALE) then
      if (.not. present(DmeshVelocity)) then
        call output_line ("Mesh velocity vector not present!", &
            OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusion3d")
        call sys_halt()
      end if
    end if

    ! Get the actual subvectors from the velocity vectors that define
    ! the X- and Y-velocity.
    if (present(IvelocityComp)) then
      Icomp = IvelocityComp
    else
      Icomp = (/1,2,3/)
    end if

    p_rvelX1 => rvecPrimary%RvectorBlock(Icomp(1))
    p_rvelY1 => rvecPrimary%RvectorBlock(Icomp(2))
    p_rvelZ1 => rvecPrimary%RvectorBlock(Icomp(3))
    p_rvelX2 => rvecSecondary%RvectorBlock(Icomp(1))
    p_rvelY2 => rvecSecondary%RvectorBlock(Icomp(2))
    p_rvelZ2 => rvecSecondary%RvectorBlock(Icomp(3))

    if (present(rsolution)) then
      p_rsolX => rsolution%RvectorBlock(Icomp(1))
      p_rsolY => rsolution%RvectorBlock(Icomp(2))
      p_rsolZ => rsolution%RvectorBlock(Icomp(3))
    else
      nullify(p_rsolX)
      nullify(p_rsolY)
      nullify(p_rsolZ)
    end if

    if (present(rdefect)) then
      p_rdefectX => rdefect%RvectorBlock(Icomp(1))
      p_rdefectY => rdefect%RvectorBlock(Icomp(2))
      p_rdefectZ => rdefect%RvectorBlock(Icomp(3))
    else
      nullify(p_rdefectX)
      nullify(p_rdefectY)
      nullify(p_rdefectZ)
    end if

    ! At the moment, we only support a rather limited set of configurations:
    ! Matrix and vectors must all be double precision, matrix must be format
    ! 7 or 9, discretisation must be Q1~, constant viscosity.
    if ((rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
        (rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
      call output_line ("Unsupported matrix format!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusion3d")
      call sys_halt()
    end if

    celement = rmatrix%p_rspatialDiscrTest%RelementDistr(1)%celement
    if (rmatrix%p_rspatialDiscrTest%ccomplexity .ne. SPDISC_UNIFORM) then
      call output_line ("Unsupported discretisation!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusion3d")
      call sys_halt()
    end if

    if ((rvecPrimary%cdataType .ne. ST_DOUBLE) .or. &
        (rvecSecondary%cdataType .ne. ST_DOUBLE)) then
      call output_line ("Unsupported vector data type in velocity!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusion3d")
      call sys_halt()
    end if

    if (present(rdefect)) then
      if ((rsolution%cdataType .ne. ST_DOUBLE) .or. &
          (rdefect%cdataType .ne. ST_DOUBLE)) then
        call output_line ("Unsupported vector data type in solution/defect!", &
            OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusion3d")
        call sys_halt()
      end if
    end if

    if (.not. rconfig%bconstViscosity) then
      call output_line ("Only constant viscosity supported at the moment!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusion3d")
      call sys_halt()
    end if

    if (rconfig%dnu .eq. SYS_INFINITY_DP) then
      call output_line ("Viscosity parameter nu not initialised!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusion3d")
      call sys_halt()
    end if

    ! Hide the p_rsol...-parameters to prevent passing the NULL()-pointer
    ! if rsolution is not present -- some compilers do not like that ^^

    call lsyssc_getbase_double (p_rvelX1,p_DvelX1)
    call lsyssc_getbase_double (p_rvelY1,p_DvelY1)
    call lsyssc_getbase_double (p_rvelZ1,p_DvelZ1)
    call lsyssc_getbase_double (p_rvelX2,p_DvelX2)
    call lsyssc_getbase_double (p_rvelY2,p_DvelY2)
    call lsyssc_getbase_double (p_rvelZ2,p_DvelZ2)

    if (present(rsolution) .and. present(rdefect)) then
      call lsyssc_getbase_double (p_rsolX   ,p_DsolX   )
      call lsyssc_getbase_double (p_rsolY   ,p_DsolY   )
      call lsyssc_getbase_double (p_rsolZ   ,p_DsolZ   )
      call lsyssc_getbase_double (p_rdefectX,p_DdefectX)
      call lsyssc_getbase_double (p_rdefectY,p_DdefectY)
      call lsyssc_getbase_double (p_rdefectZ,p_DdefectZ)

      call conv_strdiff3dALE_double ( &
                    p_DvelX1,p_DvelY1,p_DvelZ1,p_DvelX2,p_DvelY2,p_DvelZ2,&
                    dprimWeight,dsecWeight,rmatrix,cdef,rconfig%dupsam, &
                    rconfig%dnu,rconfig%dalpha,rconfig%dbeta,rconfig%dtheta,&
                    rconfig%ddelta,rconfig%clocalH,rconfig%bALE,p_DsolX,p_DsolY,&
                    p_DsolZ,p_DdefectX,p_DdefectY,p_DdefectZ,DmeshVelocity,&
                    rcubatureInfo,rperfconfig)

    else

      call conv_strdiff3dALE_double ( &
                    p_DvelX1,p_DvelY1,p_DvelZ1,p_DvelX2,p_DvelY2,p_DvelZ2,&
                    dprimWeight,dsecWeight,rmatrix,cdef,rconfig%dupsam, &
                    rconfig%dnu,rconfig%dalpha,rconfig%dbeta,rconfig%dtheta,&
                    rconfig%ddelta,rconfig%clocalH,rconfig%bALE,&
                    DmeshVelocity=DmeshVelocity,rcubatureInfo=rcubatureInfo,&
                    rperfconfig=rperfconfig)

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine conv_strdiff3dALE_double ( &
                  u1Xvel,u1Yvel,u1Zvel,u2Xvel,u2Yvel,u2Zvel,&
                  dweight1,dweight2,rmatrix,cdef, &
                  dupsam,dnu,dalpha,dbeta,dtheta,ddelta,clocalH,bALE, &
                  Du1,Du2,Du3,Ddef1,Ddef2,Ddef3,DmeshVelocity,rcubatureInfo,&
                  rperfconfig)
!<description>
  ! Standard streamline diffusion method to set up the operator
  ! <tex>
  ! $$ dtheta  *  ( dalpha * MASS  +  dbeta * STOKES  +  ddelta * u_1 * grad(u_2) ) $$
  ! </tex>
  ! in a matrix or to build a defect vector with that.
  ! 3D-version (X-, Y- and Z-velocity), uniform <tex>$\tilde Q_1$</tex> discretisation,
  ! double precision vectors/matrix.
  !
  ! u1Xvel,u1Yvel,u1Zvel,u2Xvel,u2Yvel,u2Zvel are two velocity field vectors,
  ! (u1Xvel,u1Yvel,u1Zvel) a primary and (u2Xvel,u2Yvel,u2Zvel) a secondary
  ! velocity field.
  ! The final velocity vector field is then computed as a weighted average
  ! of these two:
  ! <tex>
  !  $$ u_1  =  dweight1 * u1vel  +  dweight2 * u2vel $$
  ! </tex>
  ! <tex>$ u_2 = rsolution(.) $</tex> defines the second velocity field inside of
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
  ! [Duarte, Formaz, Natesan; `Arbitrary Lagrangian-Euler Method
  ! for Navier-Stokes equations with moving boundaries`;
  ! Comput. Methods Appl. Mech. Engrg. 193 (2004), 4819-4836]
  !
  ! Remarks:\\
  !
  ! 1.) In a typical call of the upwinding, the caller can use:
  !     dweight1 = 1, u1Xvel/u1Yvel = velocity field
  !     dweight2 = 0, u2Xvel/u2Yvel = undefined
  !   So the upwinding scheme only uses one velocity field.
  !   Such a call e.g. adds the integral
  !                <tex> $$ ( u_1 * grad(.) , v )_{\Omega} $$ </tex>
  !   to the system matrix.\\
  !
  !  2.) In case that there are two velocity fields representing
  !   the solution (may happen in a nonstationary simulation where
  !   u1Xvel/u1Yvel represents the solution in the current and u2Xvel/u2Yvel
  !   that of the previous time step), dweight1/dweight2 defines how these both
  !   velocity vectors should be weighted to compute the actual
  !   velocity field for the assembling:
  !               <tex> $$ U_act = dweight1*u1vel + dweight2*u2vel $$ </tex>
  !   This is e.g. used for the linear extrapolation technique to
  !   reconstruct a velocity from two previous time steps...\\
  !
  !  3.) In the nonlinear iteration, as a right hand side there arises
  !   a defect vector D whose linear part can easily be assembled.
  !   However, there is a nonlinearity to be included into that vector,
  !   too. By setting cdef=1,2, this routine incorporates the nonlinearity
  !   into that vector, using the formula
  !
  !            <tex> $$ D = D - dtheta * UUx * grad (Ux) $$ </tex>
  !
  !  4.) If bALE=true, a mesh velocity field is added to the nonlinearity
  !   according to the formula  "U * grad (U-DmeshVelocity)".
  !   For bALE=false, the simple nonlinearity "U * grad (U)" is used.

!</description>

!<input>

  ! Primary X-velocity of <tex>$ u_1 $</tex>
  real(DP), dimension(:), intent(in) :: u1Xvel

  ! Primary Y-velocity of <tex>$ u_1 $</tex>
  real(DP), dimension(:), intent(in) :: u1Yvel

  ! Primary Z-velocity of <tex>$ u_1 $</tex>
  real(DP), dimension(:), intent(in) :: u1Zvel

  ! Secondary X-velocity of <tex>$ u_1 $</tex>
  real(DP), dimension(:), intent(in) :: u2Xvel

  ! Secondary Y-velocity of <tex>$ u_1 $</tex>
  real(DP), dimension(:), intent(in) :: u2Yvel

  ! Secondary Z-velocity of <tex>$ u_1 $</tex>
  real(DP), dimension(:), intent(in) :: u2Zvel

  ! Computation/defect correction method. One of the CONV_MODxxxx constants:
  ! CONV_MODMATRIX: Set up the nonlinear matrix. rmatrix must be present, the
  !                 nonlinear part is added to the matrix.
  ! CONV_MODDEFECT: Set up the nonlinear defect. rdefect and rsolution must be
  !                 present.
  ! CONV_MODBOTH  : Set up the nonlinear matrix as well as the nonlinear defect.
  !                 rmatrix, rdefect and rsolution must all be present.
  integer, intent(in) :: cdef

  ! Weighting factor for u1Xvel/u1Yvel.
  real(DP), intent(in) :: dweight1

  ! Weighting factor for u2Xvel/u2Yvel.
  real(DP), intent(in) :: dweight2

  ! dupsam  - control parameter.
  !          -1: simple upwind,
  !          =0: Samarskji upwind
  real(DP), intent(in) :: dupsam

  ! Viscosity parameter <tex>$ \nu = 1/Re $</tex> if viscosity is constant
  real(DP), intent(in) :: dnu

  ! Weighting factor for the mass matrix.
  real(DP), intent(in) :: dalpha

  ! Weighting factor for the Stokes matrix. (Stokes matrix = 1/Re * Laplace)
  real(DP), intent(in) :: dbeta

  ! Weighting factor of the convective operator: <tex>$ \theta * u*grad(u) $</tex>.
  ! For time-dependent problems, this can be set to the step size
  ! in the <tex>$ \Theta $</tex>-scheme.
  real(DP), intent(in) :: dtheta

  ! Weighting factor for the nonlinear term
  real(DP), intent(in) :: ddelta

  ! How to calculate local H?
  integer, intent(in) :: clocalH

  ! Whether or not to use the ALE method
  logical, intent(in) :: bALE

  ! OPTIONAL: Mesh velocity field. Must be present if bALE=TRUE.
  ! DmeshVelocity(1,:) gives the X-velocity of all the corner points of the mesh,
  ! DmeshVelocity(2,:) gives the Y-velocity.
  real(DP), dimension(:,:), intent(in), optional :: DmeshVelocity(:,:)

  ! OPTIONAL: X-velocity of <tex>$ u_2 $</tex>. Must be present if cdef=CONV_MODDEFECT
  ! or cdef=CONV_MODBOTH.
  real(DP), dimension(:), intent(in), optional :: Du1

  ! OPTIONAL: Y-velocity of <tex>$ u_2 $</tex>. Must be present if cdef=CONV_MODDEFECT
  ! or cdef=CONV_MODBOTH.
  real(DP), dimension(:), intent(in), optional :: Du2

  ! OPTIONAL: Z-velocity of <tex>$ u_2 $</tex>. Must be present if cdef=CONV_MODDEFECT
  ! or cdef=CONV_MODBOTH.
  real(DP), dimension(:), intent(in), optional :: Du3

  ! OPTIONAL: A scalar cubature information structure that specifies the cubature
  ! formula(s) to use. If not specified, default settings are used.
  type(t_scalarCubatureInfo), intent(in), optional, target :: rcubatureInfo

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! The system matrix. Must be format 7 or 9.
  type(t_matrixScalar), intent(inout), target :: rmatrix

  ! OPTIONAL: X-defect vector. Must be present if cdef=CONV_MODDEFECT
  ! or =CONV_MODBOTH.
  real(DP), dimension(:), intent(inout), optional :: Ddef1

  ! OPTIONAL: Y-defect vector. Must be present if cdef=CONV_MODDEFECT
  ! or =CONV_MODBOTH.
  real(DP), dimension(:), intent(inout), optional :: Ddef2

  ! OPTIONAL: Z-defect vector. Must be present if cdef=CONV_MODDEFECT
  ! or =CONV_MODBOTH.
  real(DP), dimension(:), intent(inout), optional :: Ddef3
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: indof,indofALE,IEQ,IDOFE,JDOFE,icubp
  integer :: JCOL0,IDFG,JDFG,JCOL
  integer :: IEL,IELset,IELmax
  logical, dimension(EL_MAXNDER) :: Bder,BderALE
  real(DP) :: dumax,dumaxr, du1loc, du2loc, du3loc, dunorm,db,OM,AH,denth,dre,dny
  real(DP) :: HBASI1,HBASI2,HBASI3,HBASI4,HBASJ1,HBASJ2,HBASJ3,HBASJ4,HSUMI,HSUMJ
  integer :: NVE

  ! Matrix structure arrays
  integer, dimension(:), pointer :: p_Kcol
  integer, dimension(:), pointer :: p_Kld
  real(DP), dimension(:), pointer :: p_Da

  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  real(DP), dimension(:,:), allocatable :: p_DcubPtsRef

  ! The discretisation - for easier access
  type(t_spatialDiscretisation), pointer :: p_rdiscretisation

  ! Triangulation
  type(t_triangulation), pointer :: p_rtriangulation
  real(DP), dimension(:,:), pointer :: p_DvertexCoords
  integer, dimension(:,:), pointer :: p_IverticesAtElement

  ! Number of elements in a block. Normally =NELEMSIM,
  ! except if there are less elements in the discretisation.
  integer :: nelementsPerBlock

  ! One and only element distribution
  type(t_elementDistribution), pointer :: p_relementDistribution

  ! For every cubature point on the reference element,
  ! the corresponding cubature weight
  real(DP), dimension(:), allocatable :: Domega

  ! number of cubature points on the reference element
  integer :: ncubp

  ! An element evaluation set for evaluating elements.
  type(t_evalElementSet) :: revalElementSet

  ! Arrays for saving Jacobian determinants
  real(DP), dimension(:,:), pointer :: p_Ddetj

  ! An allocateable array accepting the DOF`s of a set of elements.
  integer, dimension(:,:), allocatable, target :: Idofs, IdofsALE

  ! Allocateable arrays for the values of the basis functions -
  ! for test and trial spaces.
  real(DP), dimension(:,:,:,:), allocatable, target :: Dbas,DbasALE

  ! Local matrices, used during the assembly.
  ! Values and positions of values in the global matrix.
  integer, dimension(:,:,:), allocatable :: Kentry
  real(DP), dimension(:,:), allocatable :: Dentry

  ! A pointer to an element-number list
  integer, dimension(:), pointer :: p_IelementList

  ! Pointer to the velocity field in the cubature points.
  real(DP), dimension(:,:,:), allocatable :: Dvelocity

  ! An array with local DELTA`s, each DELTA for one element
  real(DP), dimension(:), allocatable :: DlocalDelta

  ! Type of transformation from the reference to the real element
  integer(I32) :: ctrafoType

  ! Element evaluation tag; collects some information necessary for evaluating
  ! the elements.
  integer(I32) :: cevaluationTag

  ! Pointer to the performance configuration
  type(t_perfconfig), pointer :: p_rperfconfig
  
  ! Cubature information structure
  type(t_scalarCubatureInfo), target :: rtempCubatureInfo
  type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => conv_perfconfig
    end if

    ! If we do not have it, create a cubature info structure that
    ! defines how to do the assembly.
    if (.not. present(rcubatureInfo)) then
      call spdiscr_createDefCubStructure(rmatrix%p_rspatialDiscrTrial,&
          rtempCubatureInfo,CUB_GEN_DEPR_BILFORM)
      p_rcubatureInfo => rtempCubatureInfo
    else
      p_rcubatureInfo => rcubatureInfo
    end if

    ! Initialise the derivative flags
    Bder = .false.
    Bder(DER_FUNC3D) = .true.
    Bder(DER_DERIV3D_X) = .true.
    Bder(DER_DERIV3D_Y) = .true.
    Bder(DER_DERIV3D_Z) = .true.

    ! For ALE we do not even need so much
    BderALE = .false.
    BderALE(DER_FUNC3D) = .true.

    ! Shortcut to the spatial discretisation
    p_rdiscretisation => rmatrix%p_rspatialDiscrTest

    ! Get the element distribution. Here, we can find information about
    ! the cubature formula etc...
    p_relementDistribution => p_rdiscretisation%RelementDistr(1)

    ! Get some information about the triangulation
    p_rtriangulation => p_rdiscretisation%p_rtriangulation
    call storage_getbase_double2d (p_rtriangulation%h_DvertexCoords,&
                                   p_DvertexCoords)
    call storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
                                p_IverticesAtElement)
    !call storage_getbase_int2d (p_rtriangulation%h_IedgesAtElement,&
    !                            p_IedgesAtElement)

    ! Get the number of local DOF`s for trial/test functions.
    ! We assume trial and test functions to be the same.
    indof = elem_igetNDofLoc(p_relementDistribution%celement)

    ! Get the number of local DOF`s Q1 -- we need them for ALE.
    indofALE = elem_igetNDofLoc(EL_Q1_3D)

    ! Number of local DOF`s
    NVE = elem_igetNVE(p_relementDistribution%celement)

    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it is at most
    ! NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = min(p_rperfconfig%NELEMSIM,p_rtriangulation%NEL)

    ! For cdef containing CONV_MODDEFECT, we build the defect vector
    !     D = RHS - A*U
    ! In this case, the defect(rhs vectors must be present

    if (iand(cdef,CONV_MODDEFECT) .ne. 0) then
      if (.not. (present(Ddef1) .and. present(Ddef2) .and. &
                 present(Ddef3) .and. present(Du1) .and. &
                 present(Du2) .and. present(Du3))) then
        call output_line ("Necessary arguments missing!", &
            OU_CLASS_ERROR,OU_MODE_STD,"conv_strdiff3dALE_double")
        call sys_halt()
      end if
    end if

    if (iand(cdef,CONV_MODMATRIX) .ne. 0) then
      ! Get matrix arrays
      call lsyssc_getbase_double (rmatrix,p_Da)
    end if
    call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
    call lsyssc_getbase_Kld (rmatrix,p_Kld)

    ! Get from the trial element space the type of coordinate system
    ! that is used there:
    ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)

    ! Get the number of cubature points for the cubature formula
    ncubp = cub_igetNumPts(p_rcubatureInfo%p_RinfoBlocks(1)%ccubature)

    ! Allocate two arrays for the points and the weights
    allocate(Domega(ncubp))
    allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),ncubp))

    ! Get the cubature formula
    call cub_getCubature(p_rcubatureInfo%p_RinfoBlocks(1)%ccubature,p_DcubPtsRef, Domega)

    ! Allocate an array saving the coordinates of corner vertices of elements

    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  allocate(Dbas(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    ! would lead to nonused memory blocks in these arrays during the assembly,
    ! which reduces the speed by 50%!
    allocate(Dbas(indof,elem_getMaxDerivative(p_relementDistribution%celement), &
             ncubp,nelementsPerBlock))

    ! Allocate memory for the DOF`s of all the elements.
    allocate(Idofs(indof,nelementsPerBlock))

    ! The same for the ALE-space
    allocate(DbasALE(indofALE,elem_getMaxDerivative(EL_Q1_3D), &
             ncubp,nelementsPerBlock))

    ! Allocate memory for the DOF`s of all the elements.
    allocate(IdofsALE(indofALE,nelementsPerBlock))

    ! Allocate memory for array with local DELTA`s
    allocate(DlocalDelta(nelementsPerBlock))

    ! Allocate an array saving the local matrices for all elements
    ! in an element set.
    ! We could also allocate EL_MAXNBAS*EL_MAXNBAS*NELEMSIM integers
    ! for this local matrix, but this would normally not fit to the cache
    ! anymore! indofTrial*indofTest*NELEMSIM is normally much smaller!
    allocate(Kentry(indof,indof,nelementsPerBlock))
    allocate(Dentry(indof,indof))

    ! Allocate memory for the velocity in the cubature points.
    allocate(Dvelocity(NDIM3D,ncubp,nelementsPerBlock))

    ! Initialisation of the element set.
    call elprep_init(revalElementSet)

    ! What is the reciprocal of nu? We need it later.
    if (dnu .ne. 0.0_DP) then
      dre = 1.0_DP/dnu

      ! dny gets the actual multiplier for the Laplace matrix.
      ! Remember: dbeta*Stokes = dbeta*dnu*Laplace = dny*Laplace.
      ! This may be =0.0 if the Stokes operator should not be included into
      ! the matrix.
      dny = dbeta*dnu
    else
      call output_line ("NU=0 not allowed! Set dbeta=0 to prevent Stokes operator "//&
          "from being build!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_strdiff3dALE_double")
      call sys_halt()
    end if

    ! If ddelta=0, we have to neglect the nonlinearity. In both cases,
    ! set DlocalDelta=0 which disables the nonlinear term in the assembly.
    ! If dupsam=0, we neglect the stabilisation term (central difference like
    ! discretisation), so we set DlocalDelta=0 as well.
    if ((ddelta .eq. 0.0_DP) .or. (dupsam .eq. 0.0_DP)) then
      call lalg_clearVectorDble (DlocalDelta)
    end if

    ! Calculate the maximum norm of the actual velocity field
    ! U = A1*U1 + A2*U2 into DUMAX.
    ! Round up the norm to 1D-8 if it is too small...

    dumax=0.0_DP
    if (dweight2 .eq. 0.0_DP) then
      do IEQ=1,size(u1Xvel)
        du1loc = dweight1*u1Xvel(IEQ)
        du2loc = dweight1*u1Yvel(IEQ)
        du3loc = dweight1*u1Zvel(IEQ)
        dunorm = sqrt(du1loc**2+du2loc**2+du3loc**2)
        dumax = max(DUMAX,DUNORM)
      end do
    else
      do ieq=1,size(u1Xvel)
        du1loc = dweight1*u1Xvel(IEQ)+dweight2*u2Xvel(IEQ)
        du2loc = dweight1*u1Yvel(IEQ)+dweight2*u2Yvel(IEQ)
        du3loc = dweight1*u1Zvel(IEQ)+dweight2*u2Zvel(IEQ)
        dunorm = sqrt(du1loc**2+du2loc**2+du3loc**2)
        dumax = max(dumax,dunorm)
      end do
    end if

    if (dumax.lt.1E-8_DP) dumax=1E-8_DP
    dumaxr = 1.0_DP/dumax

    ! p_IelementList must point to our set of elements in the discretisation
    ! with that combination of trial/test functions
    call storage_getbase_int (p_relementDistribution%h_IelementList, &
                              p_IelementList)

    ! Loop over the elements - blockwise.
    do IELset = 1, size(p_IelementList), nelementsPerBlock

      ! We always handle NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most NELEMSIM
      ! elements simultaneously.

      IELmax = min(size(p_IelementList),IELset-1+p_rperfconfig%NELEMSIM)

      ! The outstanding feature with finite elements is: A basis
      ! function for a DOF on one element has common support only
      ! with the DOF`s on the same element! E.g. for Q1:
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
      !     to collect all "O" DOF`s.
      !
      ! Calculate the global DOF`s into IdofsTrial / IdofsTest.
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF`s of our NELEMSIM elements simultaneously.
      call dof_locGlobMapping_mult(p_rdiscretisation, p_IelementList(IELset:IELmax), &
                                  Idofs)

      ! In case ALE is used, do this also for the ALE stuff.
      if (bALE) then
        call dof_locGlobMapping_mult(p_rdiscretisation, &
                                    p_IelementList(IELset:IELmax), &
                                    IdofsALE)
      end if

      ! Calculate local DELTA`s for streamline diffusion method.
      ! (cf. p. 121 in Turek`s CFD book).
      ! For every element, we need a local DELTA.
      ! Every local delta is weighted by the global "ddelta".
      ! If ddelta=0, we do not do anything as this disables the
      ! nonlinear term.
      ! If UPSAM=0.0, we have a central-difference like discretisation, which
      ! is one can see as the local stabilisation weight Delta is also = 0.0.
      ! In this case, we even switch of the calculation of the local Delta,
      ! as it is always =0.0, so we save a little bit time.
      if ((ddelta .ne. 0.0_DP) .and. (dupsam .ne. 0.0_DP))then
        ! How do we calculate the local H?
        if (clocalH .eq. 1) then
          ! Use length of the way a particle travels
          do IEL=1,IELmax-IELset+1
            call getLocalDeltaHexaRay (u1Xvel,u1Yvel,u1Zvel,u2Xvel,u2Yvel,u2Zvel,&
                        dweight1,dweight2,IEL+IELset-1,&
                        DUMAXR,DlocalDelta(IEL),p_IverticesAtElement,&
                        p_DvertexCoords,Idofs(:,IEL),indof,dupsam,dre)
          end do ! IEL
        else
          ! Use volume of the cell
          do IEL=1,IELmax-IELset+1
            call getLocalDeltaHexaVol (u1Xvel,u1Yvel,u1Zvel,u2Xvel,u2Yvel,u2Zvel,&
                        dweight1,dweight2,IEL+IELset-1,&
                        DUMAXR,DlocalDelta(IEL),p_IverticesAtElement,&
                        p_DvertexCoords,Idofs(:,IEL),indof,dupsam,dre)
          end do ! IEL
        end if
      end if

      ! For the assembly of the global matrix, we use a "local"
      ! approach. At first we build a "local" system matrix according
      ! to the current element. This contains all additive
      ! contributions of element IEL, which are later added at the
      ! right positions to the elements in the global system matrix.
      !
      ! We have indofTrial trial DOF`s per element and
      ! indofTest test DOF`s per element. Therefore there are
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
      do IEL=1,IELmax-IELset+1

        ! For building the local matrices, we have first to
        ! loop through the test functions (the "O"`s), as these
        ! define the rows in the matrix.
        do IDOFE=1,indof

          ! Row IDOFE of the local matrix corresponds
          ! to row=global DOF KDFG(IDOFE) in the global matrix.
          ! This is one of the the "O"`s in the above picture.
          ! Get the starting position of the corresponding row
          ! to JCOL0:

          JCOL0=p_KLD(Idofs(IDOFE,IEL))

          ! Now we loop through the other DOF`s on the current element
          ! (the "O"`s).
          ! All these have common support with our current basis function
          ! and will therefore give an additive value to the global
          ! matrix.

          do JDOFE=1,indof

            ! Get the global DOF of the "X" which interacts with
            ! our "O".

            JDFG=Idofs(JDOFE,IEL)

            ! Starting in JCOL0 (which points to the beginning of
            ! the line initially), loop through the elements in
            ! the row to find the position of column IDFG.
            ! Jump out of the do loop if we find the column.

            do JCOL=JCOL0,rmatrix%NA
              if (p_KCOL(JCOL) .eq. JDFG) exit
            end do

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

          end do ! IDOFE

        end do ! JDOFE

      end do ! IEL

      ! Ok, we found the positions of the local matrix entries
      ! that we have to change.
      ! To calculate the matrix contributions, we have to evaluate
      ! the elements to give us the values of the basis functions
      ! in all the DOF`s in all the elements in our set.
      !
      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag.
      cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)
      cevaluationTag = ior(cevaluationTag,elem_getEvaluationTag(EL_Q1_3D))

      ! In the first loop, calculate the coordinates on the reference element.
      ! In all later loops, use the precalculated information.
      if (IELset .eq. 1) then
        cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REFPOINTS)
      else
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
      end if

      ! Calculate all information that is necessary to evaluate the finite element
      ! on all cells of our subset. This includes the coordinates of the points
      ! on the cells.
      call elprep_prepareSetForEvaluation (revalElementSet,&
          cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
          ctrafoType, p_DcubPtsRef(:,1:ncubp), rperfconfig=rperfconfig)
      p_Ddetj => revalElementSet%p_Ddetj

      ! Calculate the values of the basis functions.
      call elem_generic_sim2 (p_relementDistribution%celement, &
          revalElementSet, Bder, Dbas)

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
      ! DOF`s as:
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

      if (dweight2 .eq. 0.0_DP) then

        ! Loop over all elements in the current set
        do IEL=1,IELmax-IELset+1

          ! Loop over all cubature points on the current element
          do ICUBP = 1, ncubp

            du1loc = 0.0_DP
            du2loc = 0.0_DP
            du3loc = 0.0_DP

            ! Perform a loop through the trial DOF`s.
            do JDOFE=1,indof

              ! Get the value of the (test) basis function
              ! phi_i (our "O") in the cubature point:
              db = Dbas(JDOFE,1,ICUBP,IEL)

              ! Sum up to the value in the cubature point
              JDFG = Idofs(JDOFE,IEL)
              du1loc = du1loc + u1Xvel(JDFG)*db
              du2loc = du2loc + u1Yvel(JDFG)*db
              du3loc = du3loc + u1Zvel(JDFG)*db

            end do ! JDOFE

            ! Save the computed velocity
            Dvelocity(1,ICUBP,IEL) = dweight1*du1loc
            Dvelocity(2,ICUBP,IEL) = dweight1*du2loc
            Dvelocity(3,ICUBP,IEL) = dweight1*du3loc

          end do ! ICUBP

        end do ! IEL

      else

        do IEL=1,IELmax-IELset+1

          ! Loop over all cubature points on the current element
          do ICUBP = 1, ncubp

            du1loc = 0.0_DP
            du2loc = 0.0_DP
            du3loc = 0.0_DP

            ! Perform a loop through the trial DOF`s.
            do JDOFE=1,indof

              ! Get the value of the (trial) basis function
              ! phi_i in the cubature point:
              db = Dbas(JDOFE,1,ICUBP,IEL)

              ! Sum up to the value in the cubature point
              JDFG = Idofs(JDOFE,IEL)
              du1loc = du1loc + (dweight1*u1Xvel(JDFG) + dweight2*u2Xvel(JDFG))*db
              du2loc = du2loc + (dweight1*u1Yvel(JDFG) + dweight2*u2Yvel(JDFG))*db
              du3loc = du3loc + (dweight1*u1Zvel(JDFG) + dweight2*u2Zvel(JDFG))*db

            end do ! JDOFE

            ! Save the computed velocity
            Dvelocity(1,ICUBP,IEL) = du1loc
            Dvelocity(2,ICUBP,IEL) = du2loc
            Dvelocity(3,ICUBP,IEL) = du3loc

          end do ! ICUBP

        end do ! IEL

      end if

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

      if (bALE) then

        ! Calculate the values of the basis functions in all the points
        ! on all the elements
        call elem_generic_sim2 (EL_Q1_3D, &
            revalElementSet, Bder, DbasALE)

        ! Loop over all elements in the current set
        do IEL=1,IELmax-IELset+1

          ! Loop over all cubature points on the current element
          do ICUBP = 1, ncubp

            du1loc = 0.0_DP
            du2loc = 0.0_DP
            du3loc = 0.0_DP

            ! Perform a loop through the trial DOF`s.
            do JDOFE=1,indof

              ! Get the value of the (trial) basis function
              db= Dbas(JDOFE,1,ICUBP,IEL)

              ! Sum up to the value in the cubature point
              JDFG = IdofsALE(IDOFE,IEL)
              du1loc = du1loc + DmeshVelocity(1,JDFG)*db
              du2loc = du2loc + DmeshVelocity(2,JDFG)*db
              du3loc = du3loc + DmeshVelocity(3,JDFG)*db

            end do ! JDOFE

            ! Save the computed velocity
            Dvelocity(1,ICUBP,IEL) = Dvelocity(1,ICUBP,IEL) - du1loc
            Dvelocity(2,ICUBP,IEL) = Dvelocity(2,ICUBP,IEL) - du2loc
            Dvelocity(3,ICUBP,IEL) = Dvelocity(3,ICUBP,IEL) - du3loc

          end do ! ICUBP

        end do ! IEL

      end if

      ! Ok, we now use Dvelocity as coefficient array in the assembly
      ! of a biliinear form!
      !
      ! Loop over the elements in the current set.

      do IEL=1,IELmax-IELset+1

        ! Clear the local matrix
        Dentry = 0.0_DP

        ! Loop over all cubature points on the current element
        do ICUBP = 1, ncubp

          ! Calculate the current weighting factor in the cubature formula
          ! in that cubature point.
          !
          ! Normally, we have to take the absolut value of the determinant
          ! of the mapping here!
          ! In 2D, the determinant is always positive, whereas in 3D,
          ! the determinant might be negative...
          OM = Domega(ICUBP)* abs(p_Ddetj(ICUBP,IEL))

          ! Current velocity in this cubature point:
          du1loc = Dvelocity (1,ICUBP,IEL)
          du2loc = Dvelocity (2,ICUBP,IEL)
          du3loc = Dvelocity (3,ICUBP,IEL)

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
          ! Outer loop over the DOF`s i=1..indof on our current element,
          ! which corresponds to the basis functions Phi_i:

          do IDOFE=1,indof

            ! Fetch the contributions of the (test) basis functions Phi_i
            ! (our "O")  for function value and first derivatives for the
            ! current DOF into HBASIy:

            HBASI1 = Dbas(IDOFE,1,ICUBP,IEL)
            HBASI2 = Dbas(IDOFE,2,ICUBP,IEL)
            HBASI3 = Dbas(IDOFE,3,ICUBP,IEL)
            HBASI4 = Dbas(IDOFE,4,ICUBP,IEL)

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

            HSUMI = HBASI2*du1loc + HBASI3*du2loc + HBASI4*du3loc

            ! Inner loop over the DOF`s j=1..indof, which corresponds to
            ! the basis function Phi_j:

            do JDOFE=1,indof

              !if (IDOFE.EQ.JDOFE) then

                ! Short version of the evaluation of the matrix
                ! contribution - see below for a more detailed
                ! description what is added together here!

              !  AH = ddelta*HSUMI*(DlocalDelta(IEL)*HSUMI+HBASI1) &
              !      + dny*(HBASI2**2+HBASI3**2+HBASI4**2) &
              !      + dalpha*HBASI1**2

              !else

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! (out "X") for function value and first derivatives for the
                ! current DOF into HBASJy:

                HBASJ1 = Dbas(JDOFE,1,ICUBP,IEL)
                HBASJ2 = Dbas(JDOFE,2,ICUBP,IEL)
                HBASJ3 = Dbas(JDOFE,3,ICUBP,IEL)
                HBASJ4 = Dbas(JDOFE,4,ICUBP,IEL)

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
                ! we do not have to worry about that.

                HSUMJ = HBASJ2*du1loc+HBASJ3*du2loc+HBASJ4*du3loc

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

                AH = ddelta*HSUMJ*(DlocalDelta(IEL)*HSUMI+HBASI1) &
                    + dny*(HBASI2*HBASJ2+HBASI3*HBASJ3+HBASI4*HBASJ4) &
                    + dalpha*HBASI1*HBASJ1

              !end if ! (IDOFE.EQ.JDOFE)

              ! Weighten the calculated value AH by the cubature
              ! weight OM and add it to the local matrix. After the
              ! loop over all DOF`s is finished, each entry contains
              ! the calculated integral.

              Dentry(JDOFE,IDOFE) = Dentry(JDOFE,IDOFE)+OM*AH

            end do ! IDOFE

          end do ! JDOFE

        end do ! ICUBP

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

        if (iand(cdef,CONV_MODMATRIX) .ne. 0) then
          do IDOFE=1,indof
            do JDOFE=1,indof
              p_DA(Kentry(JDOFE,IDOFE,IEL)) = p_DA(Kentry(JDOFE,IDOFE,IEL)) + &
                dtheta * Dentry(JDOFE,IDOFE)
            end do
          end do
        end if

        ! For cdef containing CONV_MODDEFECT, build the defect vector
        !     D = RHS - A*U
        ! This is done matrix free, only with the help of the local
        ! matrix.
        ! In this case, D=(D1,D2) is expected to be the RHS on
        ! entry and will be updated to be the defect vector when
        ! this routine is left.

        if (iand(cdef,CONV_MODDEFECT) .ne. 0) then
          do IDOFE=1,indof

            IDFG=Idofs(IDOFE,IEL)

            do JDOFE=1,indof

              denth = dtheta*Dentry(JDOFE,IDOFE)

              JDFG=Idofs(JDOFE,IEL)
              Ddef1(IDFG)= Ddef1(IDFG) - denth*Du1(JDFG)
              Ddef2(IDFG)= Ddef2(IDFG) - denth*Du2(JDFG)
              Ddef3(IDFG)= Ddef3(IDFG) - denth*Du3(JDFG)

            end do
          end do
        end if

      end do ! IEL

    end do ! IELset

    ! Release memory
    call elprep_releaseElementSet(revalElementSet)

    deallocate(p_DcubPtsRef)
    deallocate(Domega)
    deallocate(DlocalDelta)
    deallocate(Dvelocity)
    deallocate(Dentry)
    deallocate(Kentry)
    deallocate(IdofsALE)
    deallocate(Idofs)
    deallocate(DbasALE)
    deallocate(Dbas)

    ! Release the assembly structure if necessary.
    if (.not. present(rcubatureInfo)) then
      call spdiscr_releaseCubStructure(rtempCubatureInfo)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine conv_streamlineDiffusionBlk3d (rvecPrimary, rvecSecondary,&
                                   dprimWeight, dsecWeight, rconfig, cdef, &
                                   rmatrix, rsolution, rdefect, DmeshVelocity,&
                                   rperfconfig)

!<description>
  ! Standard streamline diffusion method to set up the operator
  ! <tex>
  ! $$ dtheta  *  ( dalpha * MASS  +  dbeta * STOKES  +  ddelta * u_1 * grad(u_2) ) $$
  ! </tex>
  ! in a matrix or to build a defect vector.
  ! 3D-version (X-, Y- and Z-velocity).
  !
  ! rvecPrimary, rvecSecondary are two velocity field vectors for the X-, Y-
  ! and Z-veclocity; IvelocityComp defines which components of these
  ! vectors contains the X-, Y- and which contains the Z-velocity.
  ! The final velocity vector field is then computed as a weighted average
  ! of these two:
  !
  !  <tex> $$ u_1  =  dprimWeight * rvecPrimary  +  dsecWeight * rvecSecondary $$ </tex>
  !
  ! <tex>$ u_2 = rsolution(.) $</tex> defines the second velocity field inside of
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
  !
  ! The routine works with a block matrix rmatrix and allows to include
  ! extended operators like the Newton operator (Frechet-derivative of the
  ! convective part).
!</description>

!<input>

  ! Primary velocity field for the computation of <tex>$ u_1 $</tex>
  type(t_vectorBlock), intent(in), target :: rvecPrimary

  ! Secondary velocity field for the computation of <tex>$ u_1 $</tex>
  type(t_vectorBlock), intent(in), target :: rvecSecondary

  ! Weighting factor for rvecPrimary.
  real(DP), intent(in) :: dprimWeight

  ! Weighting factor for rvecSecondary.
  real(DP), intent(in) :: dsecWeight

  ! Configuration block for the streamline diffusion scheme
  type(t_convStreamlineDiffusion), intent(in) :: rconfig

  ! Computation/defect correction method. One of the CONV_MODxxxx constants:
  ! CONV_MODMATRIX: Set up the nonlinear matrix. rmatrix must be present, the
  !                 nonlinear part is added to the matrix.
  ! CONV_MODDEFECT: Set up the nonlinear defect. rdefect and rsolution must be
  !                 present.
  ! CONV_MODBOTH  : Set up the nonlinear matrix as well as the nonlinear defect.
  !                 rmatrix, rdefect and rsolution must all be present.
  integer, intent(in) :: cdef

  ! OPTIONAL: Solution vector u_2.
  ! Must be present if cdef=CONV_MODDEFECT or =CONV_MODBOTH.
  type(t_vectorBlock), intent(in), target, optional :: rsolution

  ! OPTIONAL: Mesh velocity field.
  ! DmeshVelocity(1,ivt) gives the X-velocity of the mesh, i.e. the X-velocity
  !   of the corner vertex ivt.
  ! DmeshVelocity(2,ivt) gives the Y-velocity of the mesh, i.e. the Y-velocity
  !   of the corner vertex ivt.
  ! DmeshVelocity(3,ivt) gives the Z-velocity of the mesh, i.e. the Z-velocity
  !   of the corner vertex ivt.
  ! The parameter must be present if ALE is activated in the
  ! configuration parameter block by bALE=true.
  real(DP), dimension(:,:), intent(in), optional :: DmeshVelocity

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! System block matrix.
  ! The content of the matrix must be present if cdef=CONV_MODMATRIX or
  ! =CONV_MODBOTH, otherwise only the structure is used.
  ! The nonlinear operator is added to the matrix.
  ! The blocks A11,A12,A21 and A22 of this matrix are tackled by streamline
  ! diffusion.
  type(t_matrixBlock), intent(inout) :: rmatrix

  ! OPTIONAL: Defect vector.
  ! Must have the same structure as rsolution/rvecPrimary/rvecSecondary.
  ! Must be present if cdef=CONV_MODDEFECT or =CONV_MODBOTH.
  ! The nonlinear part is subtracted from this vector:
  ! <tex>$ r = r - \theta * u_1*grad(u_2) $</tex>
  type(t_vectorBlock), intent(inout), optional, target :: rdefect
!</inputoutput>

!</subroutine>

    ! local variables
    integer(I32) :: celement
    type(t_vectorScalar), pointer :: p_rvelX1,p_rvelX2,p_rvelY1,p_rvelY2,&
                                     p_rvelZ1,p_rvelZ2
    type(t_vectorScalar), pointer :: p_rsolX,p_rsolY,p_rsolZ,&
                                     p_rdefectX,p_rdefectY,p_rdefectZ
    real(DP), dimension(:), pointer :: p_DvelX1,p_DvelX2,p_DvelY1,p_DvelY2,&
                                       p_DvelZ1,p_DvelZ2
    real(DP), dimension(:), pointer :: p_DsolX,p_DsolY,p_DsolZ,&
                                       p_DdefectX,p_DdefectY,p_DdefectZ

    ! At first check the input parameters that everything is present what
    ! we need:
    if ((cdef .eq. CONV_MODDEFECT) .or. (cdef .eq. CONV_MODBOTH)) then
      if ((.not. present(rsolution)) .or. (.not. present(rdefect))) then
        call output_line ("Solution/defect vector not present!", &
            OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusionBlk3d")
        call sys_halt()
      end if
    end if

    if (rconfig%bALE) then
      if (.not. present(DmeshVelocity)) then
        call output_line ("Mesh velocity vector not present!", &
            OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusionBlk3d")
        call sys_halt()
      end if
    end if

    ! Get the actual subvectors from the velocity vectors that define
    ! the X-, Y- and Z-velocity.
    p_rvelX1 => rvecPrimary%RvectorBlock(1)
    p_rvelY1 => rvecPrimary%RvectorBlock(2)
    p_rvelZ1 => rvecPrimary%RvectorBlock(3)
    p_rvelX2 => rvecSecondary%RvectorBlock(1)
    p_rvelY2 => rvecSecondary%RvectorBlock(2)
    p_rvelZ2 => rvecSecondary%RvectorBlock(3)

    if (present(rsolution)) then
      p_rsolX => rsolution%RvectorBlock(1)
      p_rsolY => rsolution%RvectorBlock(2)
      p_rsolZ => rsolution%RvectorBlock(3)
    else
      nullify(p_rsolX)
      nullify(p_rsolY)
      nullify(p_rsolZ)
    end if

    if (present(rdefect)) then
      p_rdefectX => rdefect%RvectorBlock(1)
      p_rdefectY => rdefect%RvectorBlock(2)
      p_rdefectZ => rdefect%RvectorBlock(3)
    else
      nullify(p_rdefectX)
      nullify(p_rdefectY)
      nullify(p_rdefectZ)
    end if

    ! At the moment, we only support a rather limited set of configurations:
    ! Matrix and vectors must all be double precision, matrix must be format
    ! 7 or 9, discretisation must be Q1~, constant viscosity.
    if ((rmatrix%RmatrixBlock(1,1)%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
        (rmatrix%RmatrixBlock(1,1)%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
      call output_line ("Unsupported matrix format!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusionBlk3d")
      call sys_halt()
    end if

    if ((rmatrix%RmatrixBlock(2,2)%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
        (rmatrix%RmatrixBlock(2,2)%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
      call output_line ("Unsupported matrix format!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusionBlk3d")
      call sys_halt()
    end if

    if ((rmatrix%RmatrixBlock(3,3)%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
        (rmatrix%RmatrixBlock(3,3)%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
      call output_line ("Unsupported matrix format!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusionBlk3d")
      call sys_halt()
    end if

    if (lsysbl_isSubmatrixPresent(rmatrix,1,2) .and. &
        (rmatrix%RmatrixBlock(1,2)%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
        (rmatrix%RmatrixBlock(1,2)%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
      call output_line ("Unsupported matrix format!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusionBlk3d")
      call sys_halt()
    end if

    if (lsysbl_isSubmatrixPresent(rmatrix,2,1) .and. &
        (rmatrix%RmatrixBlock(2,1)%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
        (rmatrix%RmatrixBlock(2,1)%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
      call output_line ("Unsupported matrix format!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusionBlk3d")
      call sys_halt()
    end if

    if (lsysbl_isSubmatrixPresent(rmatrix,1,3) .and. &
        (rmatrix%RmatrixBlock(1,3)%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
        (rmatrix%RmatrixBlock(1,3)%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
      call output_line ("Unsupported matrix format!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusionBlk3d")
      call sys_halt()
    end if

    if (lsysbl_isSubmatrixPresent(rmatrix,3,1) .and. &
        (rmatrix%RmatrixBlock(3,1)%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
        (rmatrix%RmatrixBlock(3,1)%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
      call output_line ("Unsupported matrix format!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusionBlk3d")
      call sys_halt()
    end if

    if (lsysbl_isSubmatrixPresent(rmatrix,2,3) .and. &
        (rmatrix%RmatrixBlock(2,3)%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
        (rmatrix%RmatrixBlock(2,3)%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
      call output_line ("Unsupported matrix format!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusionBlk3d")
      call sys_halt()
    end if

    if (lsysbl_isSubmatrixPresent(rmatrix,3,2) .and. &
        (rmatrix%RmatrixBlock(3,2)%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
        (rmatrix%RmatrixBlock(3,2)%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
      call output_line ("Unsupported matrix format!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusionBlk3d")
      call sys_halt()
    end if

    ! If Newton must be calculated, make sure A12,A21,A13,A31,A23,A32 exists
    ! and that all Aij are independent of each other!
    if (rconfig%dnewton .ne. 0.0_DP) then
      if ((.not. lsysbl_isSubmatrixPresent(rmatrix,1,2)) .or. &
          (.not. lsysbl_isSubmatrixPresent(rmatrix,2,1)) .or. &
          (.not. lsysbl_isSubmatrixPresent(rmatrix,1,3)) .or. &
          (.not. lsysbl_isSubmatrixPresent(rmatrix,3,1)) .or. &
          (.not. lsysbl_isSubmatrixPresent(rmatrix,2,3)) .or. &
          (.not. lsysbl_isSubmatrixPresent(rmatrix,3,2))) then
        call output_line ("For the Newton matrix, A12 and A21 must be defined!", &
            OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusionBlk3d")
        call sys_halt()
      end if
      if (lsyssc_isMatrixContentShared (rmatrix%RmatrixBlock(1,1),&
                                        rmatrix%RmatrixBlock(2,2)) .or. &
          lsyssc_isMatrixContentShared (rmatrix%RmatrixBlock(1,1),&
                                        rmatrix%RmatrixBlock(3,3)) .or. &
          lsyssc_isMatrixContentShared (rmatrix%RmatrixBlock(2,2),&
                                        rmatrix%RmatrixBlock(3,3)) .or. &
          lsyssc_isMatrixContentShared (rmatrix%RmatrixBlock(1,2),&
                                        rmatrix%RmatrixBlock(2,1)) .or. &
          lsyssc_isMatrixContentShared (rmatrix%RmatrixBlock(1,3),&
                                        rmatrix%RmatrixBlock(3,1)) .or. &
          lsyssc_isMatrixContentShared (rmatrix%RmatrixBlock(2,3),&
                                        rmatrix%RmatrixBlock(3,2))) then
        call output_line ("For the Newton matrix, the matrix blocks must be indepentent!", &
            OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusionBlk3d")
        call sys_halt()
      end if
    end if

    celement = rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest% &
                RelementDistr(1)%celement
    if (rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest%ccomplexity &
        .ne. SPDISC_UNIFORM) then
      call output_line ("Unsupported discretisation!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusionBlk3d")
      call sys_halt()
    end if

    if ((rvecPrimary%cdataType .ne. ST_DOUBLE) .or. &
        (rvecSecondary%cdataType .ne. ST_DOUBLE)) then
      call output_line ("Unsupported vector data type in velocity!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusionBlk3d")
      call sys_halt()
    end if

    if (present(rdefect)) then
      if ((rsolution%cdataType .ne. ST_DOUBLE) .or. &
          (rdefect%cdataType .ne. ST_DOUBLE)) then
        call output_line ("Unsupported vector data type in solution/defect!", &
            OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusionBlk3d")
        call sys_halt()
      end if
    end if

    if (.not. rconfig%bconstViscosity) then
      call output_line ("Only constant viscosity supported at the moment!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusionBlk3d")
      call sys_halt()
    end if

    if (rconfig%dnu .eq. SYS_INFINITY_DP) then
      call output_line ("Viscosity parameter nu not initialised!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_streamlineDiffusionBlk3d")
      call sys_halt()
    end if

    ! Hide the p_rsol...-parameters to prevent passing the NULL()-pointer
    ! if rsolution is not present -- some compilers do not like that ^^

    call lsyssc_getbase_double (p_rvelX1,p_DvelX1)
    call lsyssc_getbase_double (p_rvelY1,p_DvelY1)
    call lsyssc_getbase_double (p_rvelZ1,p_DvelZ1)
    call lsyssc_getbase_double (p_rvelX2,p_DvelX2)
    call lsyssc_getbase_double (p_rvelY2,p_DvelY2)
    call lsyssc_getbase_double (p_rvelZ2,p_DvelZ2)

    if (present(rdefect)) then
      call lsyssc_getbase_double (p_rsolX   ,p_DsolX   )
      call lsyssc_getbase_double (p_rsolY   ,p_DsolY   )
      call lsyssc_getbase_double (p_rsolZ   ,p_DsolZ   )
      call lsyssc_getbase_double (p_rdefectX,p_DdefectX)
      call lsyssc_getbase_double (p_rdefectY,p_DdefectY)
      call lsyssc_getbase_double (p_rdefectZ,p_DdefectZ)

      call conv_strdiff3dALEblk_double (&
              p_DvelX1,p_DvelY1,p_DvelZ1,p_DvelX2,p_DvelY2,p_DvelZ2, &
              dprimWeight, dsecWeight, rmatrix,cdef, rconfig%dupsam, &
              rconfig%dnu, rconfig%dalpha, rconfig%dbeta, rconfig%dtheta,&
              rconfig%ddelta, rconfig%dnewton, rconfig%clocalH,rconfig%bALE, &
              p_DsolX,p_DsolY,p_DsolZ,p_DdefectX,p_DdefectY,p_DdefectZ,&
              DmeshVelocity,rperfconfig)

    else

      call conv_strdiff3dALEblk_double ( &
                    p_DvelX1,p_DvelY1,p_DvelZ1,p_DvelX2,p_DvelY2,p_DvelZ2,&
                    dprimWeight,dsecWeight, rmatrix, cdef, rconfig%dupsam, &
                    rconfig%dnu,rconfig%dalpha, rconfig%dbeta, rconfig%dtheta,&
                    rconfig%ddelta, rconfig%dnewton, rconfig%clocalH,rconfig%bALE,&
                    DmeshVelocity=DmeshVelocity,rperfconfig=rperfconfig)

    end if

  end subroutine

  ! ***************************************************************************

!                denth = dtheta*Dentry(JDOFE,IDOFE,IEL)
!
!                JDFG=Idofs(JDOFE,IEL)
!                Ddef1(IDFG)= Ddef1(IDFG) - denth*Du1(JDFG)
!                Ddef2(IDFG)= Ddef2(IDFG) - denth*Du2(JDFG)

!<subroutine>
  subroutine conv_strdiff3dALEblk_double (&
                  u1Xvel,u1Yvel,u1Zvel,u2Xvel,u2Yvel,u2Zvel,dweight1,dweight2,&
                  rmatrix,cdef,dupsam,dnu,dalpha,dbeta,dtheta, ddelta, dnewton, &
                  clocalH,bALE, Du1,Du2,Du3,Ddef1,Ddef2,Ddef3, DmeshVelocity,&
                  rperfconfig)
!<description>
  ! Standard streamline diffusion method to set up the operator
  ! <tex>
  ! $$ dtheta  *  ( dalpha * MASS  +  dbeta * STOKES  +  ddelta * u_1 * grad(u_2) ) $$
  ! </tex>
  ! in a matrix or to build a defect vector with that.
  ! 3D-version (X-, Y- and Z-velocity), uniform <tex>$\tilde Q_1$</tex> discretisation,
  ! double precision vectors/matrix.
  !
  ! The routine supports fully coupled matrices, and the generation of the Newton
  ! matrix.
  !
  ! u1Xvel,u1Yvel, u2Xvel,u2Yvel are two velocity field vectors,
  ! (u1Xvel,u1Yvel) a primary and (u2Xvel,u2Yvel) a secondary velocity field.
  ! The final velocity vector field is then computed as a weighted average
  ! of these two:
  ! <tex>
  !  $$ u_1  =  dweight1 * u1vel  +  dweight2 * u2vel $$
  ! </tex>
  ! <tex>$ u_2 = rsolution(.) $</tex> defines the second velocity field inside of
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
  ! [Duarte, Formaz, Natesan; `Arbitrary Lagrangian-Euler Method
  ! for Navier-Stokes equations with moving boundaries`;
  ! Comput. Methods Appl. Mech. Engrg. 193 (2004), 4819-4836]
  !
  ! Remarks:\\
  !
  ! 1.) In a typical call of the upwinding, the caller can use:
  !     dweight1 = 1, u1Xvel/u1Yvel = velocity field
  !     dweight2 = 0, u2Xvel/u2Yvel = undefined
  !   So the upwinding scheme only uses one velocity field.
  !   Such a call e.g. adds the integral
  !                <tex> $$ ( u_1 * grad(.) , v )_{\Omega} $$ </tex>
  !   to the system matrix.\\
  !
  !  2.) In case that there are two velocity fields representing
  !   the solution (may happen in a nonstationary simulation where
  !   u1Xvel/u1Yvel represents the solution in the current and u2Xvel/u2Yvel
  !   that of the previous time step), dweight1/dweight2 defines how these both
  !   velocity vectors should be weighted to compute the actual
  !   velocity field for the assembling:
  !               <tex> $$ U_act = dweight1*u1vel + dweight2*u2vel $$ </tex>
  !   This is e.g. used for the linear extrapolation technique to
  !   reconstruct a velocity from two previous time steps...\\
  !
  !  3.) In the nonlinear iteration, as a right hand side there arises
  !   a defect vector D, which linear part can easily being assembled.
  !   However, there is a nonlinearity to be included into that vector,
  !   too. By setting cdef=1,2, this routine incorporates the nonlinearity
  !   into that vector, using the formula
  !
  !            <tex> $$ D = D - dtheta * UUx * grad (Ux) $$ </tex>
  !
  !  4.) If bALE=true, a mesh velocity field is added to the nonlineareity
  !   according to the formula  "U * grad (U-DmeshVelocity)".
  !   For bALE=false, the simple nonlinearity "U * grad (U)" is used.

!</description>

!<input>

  ! Primary X-velocity of <tex>$ u_1 $</tex>
  real(DP), dimension(:), intent(in) :: u1Xvel

  ! Primary Y-velocity of <tex>$ u_1 $</tex>
  real(DP), dimension(:), intent(in) :: u1Yvel

  ! Primary Z-velocity of <tex>$ u_1 $</tex>
  real(DP), dimension(:), intent(in) :: u1Zvel

  ! Secondary X-velocity of <tex>$ u_1 $</tex>
  real(DP), dimension(:), intent(in) :: u2Xvel

  ! Secondary Y-velocity of <tex>$ u_1 $</tex>
  real(DP), dimension(:), intent(in) :: u2Yvel

  ! Secondary Z-velocity of <tex>$ u_1 $</tex>
  real(DP), dimension(:), intent(in) :: u2Zvel

  ! Computation/defect correction method. One of the CONV_MODxxxx constants:
  ! CONV_MODMATRIX: Set up the nonlinear matrix. rmatrix must be present, the
  !                 nonlinear part is added to the matrix.
  ! CONV_MODDEFECT: Set up the nonlinear defect. rdefect and rsolution must be
  !                 present.
  ! CONV_MODBOTH  : Set up the nonlinear matrix as well as the nonlinear defect.
  !                 rmatrix, rdefect and rsolution must all be present.
  integer, intent(in) :: cdef

  ! Weighting factor for u1Xvel/u1Yvel.
  real(DP), intent(in) :: dweight1

  ! Weighting factor for u2Xvel/u2Yvel.
  real(DP), intent(in) :: dweight2

  ! dupsam  - control parameter.
  !          -1: simple upwind,
  !          =0: Samarskji upwind
  real(DP), intent(in) :: dupsam

  ! Viscosity parameter <tex>$ \nu = 1/Re $</tex> if viscosity is constant
  real(DP), intent(in) :: dnu

  ! Weighting factor for the mass matrix.
  real(DP), intent(in) :: dalpha

  ! Weighting factor for the Stokes matrix. (Stokes matrix = 1/Re * Laplace)
  real(DP), intent(in) :: dbeta

  ! Weighting factor of the convective operator: <tex>$ \theta * u*grad(u) $</tex>.
  ! For time-dependent problems, this can be set to the step size
  ! in the <tex>$ \Theta $</tex>-scheme.
  real(DP), intent(in) :: dtheta

  ! Weighting factor for the nonlinear term
  real(DP), intent(in) :: ddelta

  ! Weighting factor of the Newton matrix. A value of 0.0 deactivates the
  ! Newton part. A value != 0.0 activates Newton; in this case the submatrices
  ! A12 and A21 must be present in rmatrix.
  real(DP), intent(in) :: dnewton

  ! How to calculate the local H?
  integer, intent(in) :: clocalH

  ! Whether or not to use the ALE method
  logical, intent(in) :: bALE

  ! OPTIONAL: Mesh velocity field. Must be present if bALE=TRUE.
  ! DmeshVelocity(1,:) gives the X-velocity of all the corner points of the mesh,
  ! DmeshVelocity(2,:) gives the Y-velocity,
  ! DmeshVelocity(3,:) gives the Z-velocity.
  real(DP), dimension(:,:), intent(in), optional :: DmeshVelocity(:,:)

  ! OPTIONAL: X-velocity of <tex>$ u_2 $</tex>. Must be present if cdef=CONV_MODDEFECT
  ! or cdef=CONV_MODBOTH.
  real(DP), dimension(:), intent(in), optional :: Du1

  ! OPTIONAL: Y-velocity of <tex>$ u_2 $</tex>. Must be present if cdef=CONV_MODDEFECT
  ! or cdef=CONV_MODBOTH.
  real(DP), dimension(:), intent(in), optional :: Du2

  ! OPTIONAL: Z-velocity of <tex>$ u_2 $</tex>. Must be present if cdef=CONV_MODDEFECT
  ! or cdef=CONV_MODBOTH.
  real(DP), dimension(:), intent(in), optional :: Du3

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! The system matrix. The submatrices for the velocity must be in block
  ! A11, A12, A21 and A22 and must be in matrix format 7 or 9.
  ! A11 and A22 must have the same structure. A12 and A21 must have
  ! the same structure.
  type(t_matrixBlock), intent(inout), target :: rmatrix

  ! OPTIONAL: X-defect vector. Must be present if cdef=CONV_MODDEFECT
  ! or =CONV_MODBOTH.
  real(DP), dimension(:), intent(inout), optional :: Ddef1

  ! OPTIONAL: Y-defect vector. Must be present if cdef=CONV_MODDEFECT
  ! or =CONV_MODBOTH.
  real(DP), dimension(:), intent(inout), optional :: Ddef2

  ! OPTIONAL: Z-defect vector. Must be present if cdef=CONV_MODDEFECT
  ! or =CONV_MODBOTH.
  real(DP), dimension(:), intent(inout), optional :: Ddef3
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: indof,indofALE,IEQ,IDOFE,JDOFE,icubp
  integer :: JCOL0,IDFG,JDFG,JCOL
  integer :: IEL,IELset,IELmax
  logical, dimension(EL_MAXNDER) :: Bder,BderALE
  real(DP) :: dumax,dumaxr, du1loc, du2loc, du3loc, dunorm,db,OM,AH,denth,dre,dny
  real(DP) :: du1locx,du2locx,du3locx,du1locy,du2locy,du3locy,&
              du1locz,du2locz,du3locz,dbx,dby,dbz
  real(DP) :: AH11,AH12,AH13,AH21,AH22,AH23,AH31,AH32,AH33
  real(DP) :: HBASI1,HBASI2,HBASI3,HBASI4,HBASJ1,HBASJ2,HBASJ3,HBASJ4,HSUMI,HSUMJ
  integer :: NVE

  ! Matrix structure arrays
  integer, dimension(:), pointer :: p_Kcol
  integer, dimension(:), pointer :: p_Kld
  real(DP), dimension(:), pointer :: p_Da11,p_Da12,p_Da13,&
                 p_Da21,p_Da22,p_Da23,p_Da31,p_Da32,p_Da33

  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  real(DP), dimension(:,:), pointer :: p_DcubPtsRef

  ! The discretisation - for easier access
  type(t_spatialDiscretisation), pointer :: p_rdiscretisation

  ! Triangulation
  type(t_triangulation), pointer :: p_rtriangulation
  real(DP), dimension(:,:), pointer :: p_DvertexCoords
  integer, dimension(:,:), pointer :: p_IfacesAtElement,&
                                    p_IverticesAtElement

  ! Number of elements in a block. Normally =NELEMSIM,
  ! except if there are less elements in the discretisation.
  integer :: nelementsPerBlock

  ! One and only element distribution
  type(t_elementDistribution), pointer :: p_relementDistribution

  ! For every cubature point on the reference element,
  ! the corresponding cubature weight
  real(DP), dimension(:), allocatable :: Domega

  ! number of cubature points on the reference element
  integer :: ncubp

  ! An element evaluation set for evaluating elements.
  type(t_evalElementSet) :: revalElementSet
  logical :: bcubPtsInitialised

  ! Arrays for saving Jacobian determinants and matrices
  real(DP), dimension(:,:), pointer :: p_Ddetj

  ! An allocateable array accepting the DOF`s of a set of elements.
  integer, dimension(:,:), allocatable, target :: Idofs, IdofsALE

  ! Allocateable arrays for the values of the basis functions -
  ! for test and trial spaces.
  real(DP), dimension(:,:,:,:), allocatable, target :: Dbas,DbasALE

  ! Local matrices, used during the assembly.
  ! Values and positions of values in the global matrix.
  integer, dimension(:,:,:), allocatable :: Kentry
  real(DP), dimension(:,:,:), allocatable :: Dentry

  ! Additional contributions for the submatrices Aij stemming from Newton.
  real(DP), dimension(:,:,:), allocatable :: DentryA11
  real(DP), dimension(:,:,:), allocatable :: DentryA12
  real(DP), dimension(:,:,:), allocatable :: DentryA13
  real(DP), dimension(:,:,:), allocatable :: DentryA21
  real(DP), dimension(:,:,:), allocatable :: DentryA22
  real(DP), dimension(:,:,:), allocatable :: DentryA23
  real(DP), dimension(:,:,:), allocatable :: DentryA31
  real(DP), dimension(:,:,:), allocatable :: DentryA32
  real(DP), dimension(:,:,:), allocatable :: DentryA33

  ! A pointer to an element-number list
  integer, dimension(:), pointer :: p_IelementList

  ! Pointer to the velocity field in the cubature points.
  real(DP), dimension(:,:,:), allocatable :: Dvelocity

  ! Pointer to the velocity X-, Y- and Z-derivative in the cubature points
  real(DP), dimension(:,:,:), allocatable :: DvelocityUderiv
  real(DP), dimension(:,:,:), allocatable :: DvelocityVderiv
  real(DP), dimension(:,:,:), allocatable :: DvelocityWderiv

  ! An array with local DELTA`s, each DELTA for one element
  real(DP), dimension(:), allocatable :: DlocalDelta

  ! Type of transformation from the reference to the real element
  integer(I32) :: ctrafoType

  ! Element evaluation tag; collects some information necessary for evaluating
  ! the elements.
  integer(I32) :: cevaluationTag

  ! Pointer to the performance configuration
  type(t_perfconfig), pointer :: p_rperfconfig

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => conv_perfconfig
    end if

    ! Initialise the derivative flags
    Bder = .false.
    Bder(DER_FUNC3D) = .true.
    Bder(DER_DERIV3D_X) = .true.
    Bder(DER_DERIV3D_Y) = .true.
    Bder(DER_DERIV3D_Z) = .true.

    ! For ALE we do not even need so much
    BderALE = .false.
    BderALE(DER_FUNC3D) = .true.
    BderALE(DER_DERIV3D_X) = .true.
    !BderALE(DER_DERIV3D_X) = .TRUE.

    ! Shortcut to the spatial discretisation.
    ! We assume the same for all Aij.
    p_rdiscretisation => rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest

    ! Get the element distribution. Here, we can find information about
    ! the cubature formula etc...
    p_relementDistribution => p_rdiscretisation%RelementDistr(1)

    ! Get some information about the triangulation
    p_rtriangulation => p_rdiscretisation%p_rtriangulation
    call storage_getbase_double2d (p_rtriangulation%h_DvertexCoords,&
                                   p_DvertexCoords)
    call storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
                                p_IverticesAtElement)
    call storage_getbase_int2d (p_rtriangulation%h_IfacesAtElement,&
                                p_IfacesAtElement)

    ! Get the number of local DOF`s for trial/test functions.
    ! We assume trial and test functions to be the same.
    indof = elem_igetNDofLoc(p_relementDistribution%celement)

    ! Get the number of local DOF`s Q1 -- we need them for ALE.
    indofALE = elem_igetNDofLoc(EL_Q1_3D)

    ! Number of local DOF`s
    NVE = elem_igetNVE(p_relementDistribution%celement)

    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it is at most
    ! NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = min(p_rperfconfig%NELEMSIM,p_rtriangulation%NEL)

    ! For cdef containing CONV_MODDEFECT, we build the defect vector
    !     D = RHS - A*U
    ! In this case, the defect(rhs vectors must be present

    if (iand(cdef,CONV_MODDEFECT) .ne. 0) then
      if (.not. (present(Ddef1) .and. present(Ddef2) .and. present(Ddef3) .and. &
                 present(Du1) .and. present(Du2) .and. present(Du3))) then
        call output_line ("Necessary arguments missing!", &
            OU_CLASS_ERROR,OU_MODE_STD,"conv_strdiff3dALEblk_double")
        call sys_halt()
      end if
    end if

    ! Get pointers to the matrix content (if necessary)
    if (iand(cdef,CONV_MODMATRIX) .ne. 0) then
      ! Get matrix arrays
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,1),p_Da11)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,2),p_Da22)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(3,3),p_Da33)

      if (dnewton .ne. 0.0_DP) then
        call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,2),p_Da12)
        call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,1),p_Da21)
        call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,3),p_Da13)
        call lsyssc_getbase_double (rmatrix%RmatrixBlock(3,1),p_Da31)
        call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,3),p_Da23)
        call lsyssc_getbase_double (rmatrix%RmatrixBlock(3,2),p_Da32)
      else
        nullify(p_Da12,p_Da21,p_Da13,p_Da31,p_Da23,p_Da32)
      end if
    end if

    ! Get pointers to the matrix structure(s).
    call lsyssc_getbase_Kcol (rmatrix%RmatrixBlock(1,1),p_Kcol)
    call lsyssc_getbase_Kld (rmatrix%RmatrixBlock(1,1),p_Kld)

    ! Get from the trial element space the type of coordinate system
    ! that is used there:
    ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)

    ! Get the number of cubature points for the cubature formula
    ncubp = cub_igetNumPts(p_relementDistribution%ccubTypeBilForm)

    ! Allocate two arrays for the points and the weights
    allocate(Domega(ncubp))
    allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),ncubp))

    ! Get the cubature formula
    call cub_getCubature(p_relementDistribution%ccubTypeBilForm,p_DcubPtsRef, Domega)

    ! OpenMP-Extension: Open threads here.
    ! Each thread will allocate its own local memory...

    !%omp parallel private(csysTrial, p_DcubPtsReal, &
    !%omp p_Ddetj, j,i,k,Dbas,Idofs,DbasALE, &
    !%omp IdofsALE,DlocalDelta,bnonpar,Kentry,Kentry12,Dentry, &
    !%omp DentryA11,DentryA12,DentryA21,DentryA22,Dvelocity, &
    !%omp DvelocityUderiv,DvelocityVderiv,dre,IEL,db,icubp,&
    !%omp IDOFE,JCOL0,JDOFE,JDFG,jcol,du1loc,du2loc,dbx,dby, &
    !%omp du1locx,du1locy,du2locx,du2locy,OM,AH,HBASI1,HBASI2,&
    !%omp HBASI3,HBASJ1,HBASJ2,HBASJ3,HSUMI,HSUMJ,AH11,AH12,AH21, &
    !%omp AH22,IELmax,revalElementSet,dny,p_DcubPts)

    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  allocate(Dbas(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    ! would lead to nonused memory blocks in these arrays during the assembly,
    ! which reduces the speed by 50%!
    allocate(Dbas(indof,elem_getMaxDerivative(p_relementDistribution%celement), &
             ncubp,nelementsPerBlock))

    ! Allocate memory for the DOF`s of all the elements.
    allocate(Idofs(indof,nelementsPerBlock))

    ! The same for the ALE-space
    allocate(DbasALE(indofALE,elem_getMaxDerivative(EL_Q1_3D), &
             ncubp,nelementsPerBlock))

    ! Allocate memory for the DOF`s of all the elements.
    allocate(IdofsALE(indofALE,nelementsPerBlock))

    ! Allocate memory for array with local DELTA`s
    allocate(DlocalDelta(nelementsPerBlock))

    ! Allocate an array saving the local matrices for all elements
    ! in an element set.
    ! We could also allocate EL_MAXNBAS*EL_MAXNBAS*NELEMSIM integers
    ! for this local matrix, but this would normally not fit to the cache
    ! anymore! indofTrial*indofTest*NELEMSIM is normally much smaller!
    !
    ! Kentry (:,:,:) defines the positions of the local matrices
    ! in the submatrices Aij.
    allocate(Kentry(indof,indof,nelementsPerBlock))

    ! Dentry (:,:,:) fetches the 'main' matrix entries (Laplace, Mass,
    ! Convection).
    ! DentryA11, DentryA12, DentryA21 and DentryA22 fetches additional entries in
    ! A11, A12, A21 and A22 of the Newton matrix, which is not always calculated
    ! and therefore not always used!
    allocate(Dentry(indof,indof,nelementsPerBlock))

    if (dnewton .ne. 0.0_DP) then
      allocate(DentryA11(indof,indof,nelementsPerBlock))
      allocate(DentryA12(indof,indof,nelementsPerBlock))
      allocate(DentryA13(indof,indof,nelementsPerBlock))
      allocate(DentryA21(indof,indof,nelementsPerBlock))
      allocate(DentryA22(indof,indof,nelementsPerBlock))
      allocate(DentryA23(indof,indof,nelementsPerBlock))
      allocate(DentryA31(indof,indof,nelementsPerBlock))
      allocate(DentryA32(indof,indof,nelementsPerBlock))
      allocate(DentryA33(indof,indof,nelementsPerBlock))
    end if

    ! Allocate memory for the velocity in the cubature points.
    allocate(Dvelocity(NDIM3D,ncubp,nelementsPerBlock))

    if (dnewton .ne. 0.0_DP) then
      allocate(DvelocityUderiv(NDIM3D,ncubp,nelementsPerBlock))
      allocate(DvelocityVderiv(NDIM3D,ncubp,nelementsPerBlock))
      allocate(DvelocityWderiv(NDIM3D,ncubp,nelementsPerBlock))
    end if

    ! Initialisation of the element set.
    call elprep_init(revalElementSet)

    ! Indicate that cubature points must still be initialised in the element set.
    bcubPtsInitialised = .false.

    ! What is the reciprocal of nu? We need it later.
    if (dnu .ne. 0.0_DP) then
      dre = 1.0_DP/dnu

      ! dny gets the actual multiplier for the Laplace matrix.
      ! Remember: dbeta*Stokes = dbeta*dnu*Laplace = dny*Laplace.
      ! This may be =0.0 if the Stokes operator should not be included into
      ! the matrix.
      dny = dbeta*dnu
    else
      call output_line ("NU=0 not allowed! Set dbeta=0 to prevent Stokes operator "//&
          "from being build!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_strdiff3dALEblk_double")
      call sys_halt()
    end if

    ! If ddelta=0, we have to neglect the nonlinearity. In both cases,
    ! set DlocalDelta=0 which disables the nonlinear term in the assembly.
    ! If dupsam=0, we neglect the stabilisation term (central difference like
    ! discretisation), so we set DlocalDelta=0 as well.
    if ((ddelta .eq. 0.0_DP) .or. (dupsam .eq. 0.0_DP)) then
      call lalg_clearVectorDble (DlocalDelta)
    end if

    ! Calculate the maximum norm of the actual velocity field
    ! U = A1*U1 + A2*U2 into DUMAX.
    ! Round up the norm to 1D-8 if it is too small...
    !%omp single
    dumax=0.0_DP
    if (dweight2 .eq. 0.0_DP) then


      do IEQ=1,size(u1Xvel)
        du1loc = dweight1*u1Xvel(IEQ)
        du2loc = dweight1*u1Yvel(IEQ)
        du3loc = dweight1*u1Zvel(IEQ)
        dunorm = sqrt(du1loc**2 + du2loc**2 + du3loc**2)
        dumax = max(DUMAX,DUNORM)
      end do

    else

      do ieq=1,size(u1Xvel)
        du1loc = dweight1*u1Xvel(IEQ)+dweight2*u2Xvel(IEQ)
        du2loc = dweight1*u1Yvel(IEQ)+dweight2*u2Yvel(IEQ)
        du3loc = dweight1*u1Zvel(IEQ)+dweight2*u2Zvel(IEQ)
        dunorm = sqrt(du1loc**2 + du2loc**2 + du3loc**2)
        dumax = max(dumax,dunorm)
      end do

    end if

    !print *,"dumax: ",dumax
    if (dumax.lt.1E-8_DP) dumax=1E-8_DP
    dumaxr = 1.0_DP/dumax
    !%omp end single

    ! p_IelementList must point to our set of elements in the discretisation
    ! with that combination of trial/test functions
    call storage_getbase_int (p_relementDistribution%h_IelementList, &
                              p_IelementList)


    ! Loop over the elements - blockwise.
    !
    ! OpenMP-Extension: Each loop cycle is executed in a different thread,
    ! so NELEMSIM local matrices are simultaneously calculated in the
    ! inner loop(s).
    ! The blocks have all the same size, so we can use static scheduling.
    !%omp do schedule(dynamic,1)
    do IELset = 1, size(p_IelementList), nelementsPerBlock

      ! We always handle NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most NELEMSIM
      ! elements simultaneously.

      IELmax = min(size(p_IelementList),IELset-1+p_rperfconfig%NELEMSIM)

      ! The outstanding feature with finite elements is: A basis
      ! function for a DOF on one element has common support only
      ! with the DOF`s on the same element! E.g. for Q1:
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
      !     to collect all "O" DOF`s.
      !
      ! Calculate the global DOF`s into IdofsTrial / IdofsTest.
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF`s of our NELEMSIM elements simultaneously.
      call dof_locGlobMapping_mult(p_rdiscretisation, p_IelementList(IELset:IELmax), &
                                  Idofs)

      ! In case ALE is used, do this also for the ALE stuff.
      if (bALE) then
        call dof_locGlobMapping_mult(p_rdiscretisation, &
                                    p_IelementList(IELset:IELmax), &
                                    IdofsALE)
      end if

      ! Calculate local DELTA`s for streamline diffusion method.
      ! (cf. p. 121 in Turek`s CFD book).
      ! For every element, we need a local DELTA.
      ! Every local delta is weighted by the global "ddelta".
      ! If ddelta=0, we do not do anything as this disables the
      ! nonlinear term.
      ! If UPSAM=0.0, we have a central-difference like discretisation, which
      ! is one can see as the local stabilisation weight Delta is also = 0.0.
      ! In this case, we even switch of the calculation of the local Delta,
      ! as it is always =0.0, so we save a little bit time.
      if ((ddelta .ne. 0.0_DP) .and. (dupsam .ne. 0.0_DP))then
        if (clocalH .eq. 1) then
          do IEL=1,IELmax-IELset+1
            call getLocalDeltaHexaRay (u1Xvel,u1Yvel,u1Zvel,u2Xvel,u2Yvel,u2Zvel,&
                        dweight1,dweight2, IEL+IELset-1,&
                        DUMAXR,DlocalDelta(IEL),p_IverticesAtElement, &
                        p_DvertexCoords,Idofs(:,IEL),indof, dupsam,dre)
          end do ! IEL
        else
          do IEL=1,IELmax-IELset+1
            call getLocalDeltaHexaVol (u1Xvel,u1Yvel,u1Zvel,u2Xvel,u2Yvel,u2Zvel,&
                        dweight1,dweight2, IEL+IELset-1,&
                        DUMAXR,DlocalDelta(IEL),p_IverticesAtElement, &
                        p_DvertexCoords,Idofs(:,IEL),indof, dupsam,dre)
          end do ! IEL
        end if
      end if

      ! For the assembly of the global matrix, we use a "local"
      ! approach. At first we build a "local" system matrix according
      ! to the current element. This contains all additive
      ! contributions of element IEL, which are later added at the
      ! right positions to the elements in the global system matrix.
      !
      ! We have indofTrial trial DOF`s per element and
      ! indofTest test DOF`s per element. Therefore there are
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
      do IEL=1,IELmax-IELset+1

        ! For building the local matrices, we have first to
        ! loop through the test functions (the "O"`s), as these
        ! define the rows in the matrix.
        do IDOFE=1,indof

          ! Row IDOFE of the local matrix corresponds
          ! to row=global DOF KDFG(IDOFE) in the global matrix.
          ! This is one of the the "O"`s in the above picture.
          ! Get the starting position of the corresponding row
          ! to JCOL0:

          JCOL0=p_KLD(Idofs(IDOFE,IEL))

          ! Now we loop through the other DOF`s on the current element
          ! (the "O"`s).
          ! All these have common support with our current basis function
          ! and will therefore give an additive value to the global
          ! matrix.

          do JDOFE=1,indof

            ! Get the global DOF of the "X" which interacts with
            ! our "O".

            JDFG=Idofs(JDOFE,IEL)

            ! Starting in JCOL0 (which points to the beginning of
            ! the line initially), loop through the elements in
            ! the row to find the position of column IDFG.
            ! Jump out of the do loop if we find the column.

            do JCOL=JCOL0,rmatrix%RmatrixBlock(1,1)%NA
              if (p_KCOL(JCOL) .eq. JDFG) exit
            end do

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

          end do ! IDOFE

        end do ! JDOFE

      end do ! IEL

      ! Ok, we found the positions of the local matrix entries
      ! that we have to change.
      ! To calculate the matrix contributions, we have to evaluate
      ! the elements to give us the values of the basis functions
      ! in all the DOF`s in all the elements in our set.

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag.
      cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)
      cevaluationTag = ior(cevaluationTag,&
                      elem_getEvaluationTag(EL_Q1_3D))

      ! In the first loop, calculate the coordinates on the reference element.
      ! In all later loops, use the precalculated information.
      !
      ! Note: Why not using
      !   if (IELset .EQ. 1) then
      ! here, but this strange concept with the boolean variable?
      ! Because the if-command does not work with OpenMP! bcubPtsInitialised
      ! is a local variable and will therefore ensure that every thread
      ! is initialising its local set of cubature points!
      if (.not. bcubPtsInitialised) then
        bcubPtsInitialised = .true.
        cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REFPOINTS)
      else
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
      end if

      ! Calculate all information that is necessary to evaluate the finite element
      ! on all cells of our subset. This includes the coordinates of the points
      ! on the cells.
      call elprep_prepareSetForEvaluation (revalElementSet,&
          cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
          ctrafoType, p_DcubPtsRef(:,1:ncubp), rperfconfig=rperfconfig)
      p_Ddetj => revalElementSet%p_Ddetj

      ! Calculate the values of the basis functions.
      ! Pass p_DcubPts as point coordinates, which point either to the
      ! coordinates on the reference element (the same for all elements)
      ! or on the real element - depending on whether this is a
      ! parametric or nonparametric element.
      call elem_generic_sim2 (p_relementDistribution%celement, &
          revalElementSet, Bder, Dbas)

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
      ! DOF`s as:
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

      ! only primary velocity field
      if (dweight2 .eq. 0.0_DP) then
!      print *,"dweight2 .EQ. 0.0"

        ! Loop over all elements in the current set
        do IEL=1,IELmax-IELset+1

          ! Loop over all cubature points on the current element
          do ICUBP = 1, ncubp

            du1loc = 0.0_DP
            du2loc = 0.0_DP
            du3loc = 0.0_DP

            ! Perform a loop through the trial DOF`s.
            do JDOFE=1,indof

              ! Get the value of the (test) basis function
              ! phi_i (our "O") in the cubature point:

              db = Dbas(JDOFE,1,ICUBP,IEL)

              ! Sum up to the value in the cubature point

              JDFG = Idofs(JDOFE,IEL)
              du1loc = du1loc +u1Xvel(JDFG)*db
              du2loc = du2loc +u1Yvel(JDFG)*db
              du3loc = du3loc +u1Zvel(JDFG)*db

            end do ! JDOFE

            ! Save the computed velocity
            Dvelocity(1,ICUBP,IEL) = dweight1*du1loc
            Dvelocity(2,ICUBP,IEL) = dweight1*du2loc
            Dvelocity(3,ICUBP,IEL) = dweight1*du3loc

          end do ! ICUBP

        end do ! IEL

        ! Compute X-, Y- and Z-derivative of the velocity?
        if (dnewton .ne. 0.0_DP) then

          do IEL=1,IELmax-IELset+1

            ! Loop over all cubature points on the current element
            do ICUBP = 1, ncubp

              du1locx = 0.0_DP
              du1locy = 0.0_DP
              du1locz = 0.0_DP
              du2locx = 0.0_DP
              du2locy = 0.0_DP
              du2locz = 0.0_DP
              du3locx = 0.0_DP
              du3locy = 0.0_DP
              du3locz = 0.0_DP

              ! Perform a loop through the trial DOF`s.
              do JDOFE=1,indof

                ! Get the value of the (trial) basis function
                ! phi_i in the cubature point:
                dbx = Dbas(JDOFE,DER_DERIV3D_X,ICUBP,IEL)
                dby = Dbas(JDOFE,DER_DERIV3D_Y,ICUBP,IEL)
                dbz = Dbas(JDOFE,DER_DERIV3D_Z,ICUBP,IEL)

                ! Sum up to the value in the cubature point
                JDFG = Idofs(JDOFE,IEL)
                du1locx = du1locx + u1Xvel(JDFG)*dbx
                du1locy = du1locy + u1Xvel(JDFG)*dby
                du1locz = du1locz + u1Xvel(JDFG)*dbz
                du2locx = du2locx + u1Yvel(JDFG)*dbx
                du2locy = du2locy + u1Yvel(JDFG)*dby
                du2locz = du2locz + u1Yvel(JDFG)*dbz
                du3locx = du3locx + u1Zvel(JDFG)*dbx
                du3locy = du3locy + u1Zvel(JDFG)*dby
                du3locz = du3locz + u1Zvel(JDFG)*dbz

              end do ! JDOFE

              ! Save the computed velocity derivative
              DvelocityUderiv(1,ICUBP,IEL) = dweight1*du1locx
              DvelocityUderiv(2,ICUBP,IEL) = dweight1*du1locy
              DvelocityUderiv(3,ICUBP,IEL) = dweight1*du1locz
              DvelocityVderiv(1,ICUBP,IEL) = dweight1*du2locx
              DvelocityVderiv(2,ICUBP,IEL) = dweight1*du2locy
              DvelocityVderiv(3,ICUBP,IEL) = dweight1*du2locz
              DvelocityWderiv(1,ICUBP,IEL) = dweight1*du3locx
              DvelocityWderiv(2,ICUBP,IEL) = dweight1*du3locy
              DvelocityWderiv(3,ICUBP,IEL) = dweight1*du3locz

            end do ! ICUBP

          end do ! IEL

        end if ! dnewton != 0

      else

        do IEL=1,IELmax-IELset+1

          ! Loop over all cubature points on the current element
          do ICUBP = 1, ncubp

            du1loc = 0.0_DP
            du2loc = 0.0_DP
            du3loc = 0.0_DP

            ! Perform a loop through the trial DOF`s.
            do JDOFE=1,indof

              ! Get the value of the (trial) basis function
              ! phi_i in the cubature point:
              db = Dbas(JDOFE,1,ICUBP,IEL)

              ! Sum up to the value in the cubature point
              JDFG = Idofs(JDOFE,IEL)
              du1loc = du1loc + (dweight1*u1Xvel(JDFG) + dweight2*u2Xvel(JDFG))*db
              du2loc = du2loc + (dweight1*u1Yvel(JDFG) + dweight2*u2Yvel(JDFG))*db
              du3loc = du3loc + (dweight1*u1Zvel(JDFG) + dweight2*u2Zvel(JDFG))*db

            end do ! JDOFE

            ! Save the computed velocity
            Dvelocity(1,ICUBP,IEL) = du1loc
            Dvelocity(2,ICUBP,IEL) = du2loc
            Dvelocity(3,ICUBP,IEL) = du3loc

          end do ! ICUBP

        end do ! IEL

        ! Compute X-, Y- and Z-derivative of the velocity?
        if (dnewton .ne. 0.0_DP) then

          do IEL=1,IELmax-IELset+1

            ! Loop over all cubature points on the current element
            do ICUBP = 1, ncubp

              du1locx = 0.0_DP
              du1locy = 0.0_DP
              du1locz = 0.0_DP
              du2locx = 0.0_DP
              du2locy = 0.0_DP
              du2locz = 0.0_DP
              du3locx = 0.0_DP
              du3locy = 0.0_DP
              du3locz = 0.0_DP

              ! Perform a loop through the trial DOF`s.
              do JDOFE=1,indof

                ! Get the value of the (trial) basis function
                ! phi_i in the cubature point:
                dbx = Dbas(JDOFE,DER_DERIV3D_X,ICUBP,IEL)
                dby = Dbas(JDOFE,DER_DERIV3D_Y,ICUBP,IEL)
                dbz = Dbas(JDOFE,DER_DERIV3D_Z,ICUBP,IEL)

                ! Sum up to the value in the cubature point
                JDFG = Idofs(JDOFE,IEL)
                du1locx = du1locx + (dweight1*u1Xvel(JDFG) + dweight2*u2Xvel(JDFG))*dbx
                du1locy = du1locy + (dweight1*u1Xvel(JDFG) + dweight2*u2Xvel(JDFG))*dby
                du1locz = du1locz + (dweight1*u1Xvel(JDFG) + dweight2*u2Xvel(JDFG))*dbz
                du2locx = du2locx + (dweight1*u1Yvel(JDFG) + dweight2*u2Yvel(JDFG))*dbx
                du2locy = du2locy + (dweight1*u1Yvel(JDFG) + dweight2*u2Yvel(JDFG))*dby
                du2locz = du2locz + (dweight1*u1Yvel(JDFG) + dweight2*u2Yvel(JDFG))*dbz
                du3locx = du3locx + (dweight1*u1Zvel(JDFG) + dweight2*u2Zvel(JDFG))*dbx
                du3locy = du3locy + (dweight1*u1Zvel(JDFG) + dweight2*u2Zvel(JDFG))*dby
                du3locz = du3locz + (dweight1*u1Zvel(JDFG) + dweight2*u2Zvel(JDFG))*dbz

              end do ! JDOFE

              ! Save the computed velocity derivative
              DvelocityUderiv(1,ICUBP,IEL) = du1locx
              DvelocityUderiv(2,ICUBP,IEL) = du1locy
              DvelocityUderiv(3,ICUBP,IEL) = du1locz
              DvelocityVderiv(1,ICUBP,IEL) = du2locx
              DvelocityVderiv(2,ICUBP,IEL) = du2locy
              DvelocityVderiv(3,ICUBP,IEL) = du2locz
              DvelocityWderiv(1,ICUBP,IEL) = du3locx
              DvelocityWderiv(2,ICUBP,IEL) = du3locy
              DvelocityWderiv(3,ICUBP,IEL) = du3locz

            end do ! ICUBP

          end do ! IEL

        end if ! dnewton != 0

      end if

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

      if (bALE) then

        ! Calculate the values of the basis functions in all the points
        ! on all the elements
        call elem_generic_sim2 (EL_Q1_3D, &
            revalElementSet, Bder, DbasALE)

        ! Loop over all elements in the current set
        do IEL=1,IELmax-IELset+1

          ! Loop over all cubature points on the current element
          do ICUBP = 1, ncubp

            du1loc = 0.0_DP
            du2loc = 0.0_DP
            du3loc = 0.0_DP

            ! Perform a loop through the trial DOF`s.
            do JDOFE=1,indof

              ! Get the value of the (trial) basis function
              db= Dbas(JDOFE,1,ICUBP,IEL)

              ! Sum up to the value in the cubature point
              JDFG = IdofsALE(IDOFE,IEL)
              du1loc = du1loc + DmeshVelocity(1,JDFG)*db
              du2loc = du2loc + DmeshVelocity(2,JDFG)*db
              du3loc = du3loc + DmeshVelocity(3,JDFG)*db

            end do ! JDOFE

            ! Save the computed velocity
            Dvelocity(1,ICUBP,IEL) = Dvelocity(1,ICUBP,IEL) - du1loc
            Dvelocity(2,ICUBP,IEL) = Dvelocity(2,ICUBP,IEL) - du2loc
            Dvelocity(3,ICUBP,IEL) = Dvelocity(3,ICUBP,IEL) - du3loc

          end do ! ICUBP

        end do ! IEL

        ! Subtract the X-, Y- and Z-derivative of the mesh velocity to the
        ! velocity derivative field if Newton is active.
        if (dnewton .ne. 0.0_DP) then

          do IEL=1,IELmax-IELset+1

            ! Loop over all cubature points on the current element
            do ICUBP = 1, ncubp

              du1locx = 0.0_DP
              du1locy = 0.0_DP
              du1locz = 0.0_DP
              du2locx = 0.0_DP
              du2locy = 0.0_DP
              du2locz = 0.0_DP
              du3locx = 0.0_DP
              du3locy = 0.0_DP
              du3locz = 0.0_DP

              ! Perform a loop through the trial DOF`s.
              do JDOFE=1,indof

                ! Get the value of the (trial) basis function
                ! phi_i in the cubature point:
                dbx = Dbas(JDOFE,DER_DERIV3D_X,ICUBP,IEL)
                dby = Dbas(JDOFE,DER_DERIV3D_Y,ICUBP,IEL)
                dbz = Dbas(JDOFE,DER_DERIV3D_Z,ICUBP,IEL)

                ! Sum up to the value in the cubature point
                JDFG = Idofs(JDOFE,IEL)
                du1locx = du1locx + DmeshVelocity(1,JDFG)*dbx
                du1locy = du1locy + DmeshVelocity(1,JDFG)*dby
                du1locz = du1locz + DmeshVelocity(1,JDFG)*dbz
                du2locx = du2locx + DmeshVelocity(2,JDFG)*dbx
                du2locy = du2locy + DmeshVelocity(2,JDFG)*dby
                du2locz = du2locz + DmeshVelocity(2,JDFG)*dbz
                du3locx = du3locx + DmeshVelocity(3,JDFG)*dbx
                du3locy = du3locy + DmeshVelocity(3,JDFG)*dby
                du3locz = du3locz + DmeshVelocity(3,JDFG)*dbz

              end do ! JDOFE

              ! Subtract the velocity derivative to the previously calculated one.
              DvelocityUderiv(1,ICUBP,IEL) = DvelocityUderiv(1,ICUBP,IEL)-du1locx
              DvelocityUderiv(2,ICUBP,IEL) = DvelocityUderiv(2,ICUBP,IEL)-du1locy
              DvelocityUderiv(3,ICUBP,IEL) = DvelocityUderiv(3,ICUBP,IEL)-du1locz
              DvelocityVderiv(1,ICUBP,IEL) = DvelocityVderiv(1,ICUBP,IEL)-du2locx
              DvelocityVderiv(2,ICUBP,IEL) = DvelocityVderiv(2,ICUBP,IEL)-du2locy
              DvelocityVderiv(3,ICUBP,IEL) = DvelocityVderiv(3,ICUBP,IEL)-du2locz
              DvelocityWderiv(1,ICUBP,IEL) = DvelocityWderiv(1,ICUBP,IEL)-du3locx
              DvelocityWderiv(2,ICUBP,IEL) = DvelocityWderiv(2,ICUBP,IEL)-du3locy
              DvelocityWderiv(3,ICUBP,IEL) = DvelocityWderiv(3,ICUBP,IEL)-du3locz

            end do ! ICUBP

          end do ! IEL

        end if ! dnewton != 0

      end if

      ! Ok, we now use Dvelocity as coefficient array in the assembly
      ! of a bilinear form!
      !
      ! Clear the local matrices. If the Newton part is to be calculated,
      ! we must clear everything, otherwise only Dentry.
      Dentry = 0.0_DP
      if (dnewton .ne. 0) then
        DentryA11 = 0.0_DP
        DentryA12 = 0.0_DP
        DentryA13 = 0.0_DP
        DentryA21 = 0.0_DP
        DentryA22 = 0.0_DP
        DentryA23 = 0.0_DP
        DentryA31 = 0.0_DP
        DentryA32 = 0.0_DP
        DentryA33 = 0.0_DP
      end if

      ! If ddelta != 0, set up the nonlinearity U*grad(u), probably with
      ! streamline diffusion stabilisation.
      if (ddelta .ne. 0.0_DP) then

        ! Loop over the elements in the current set.
        do IEL=1,IELmax-IELset+1

          ! Loop over all cubature points on the current element
          do ICUBP = 1, ncubp

            ! Calculate the current weighting factor in the cubature formula
            ! in that cubature point.
            !
            ! Normally, we have to take the absolut value of the determinant
            ! of the mapping here!
            ! In 2D, the determinant is always positive, whereas in 3D,
            ! the determinant might be negative -- that is normal!
            ! But because this routine only works in 2D, we can skip
            ! the ABS here!

            OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))

            ! Current velocity in this cubature point:
            du1loc = Dvelocity (1,ICUBP,IEL)
            du2loc = Dvelocity (2,ICUBP,IEL)
            du3loc = Dvelocity (3,ICUBP,IEL)

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
            ! Outer loop over the DOF`s i=1..indof on our current element,
            ! which corresponds to the basis functions Phi_i:

            do IDOFE=1,indof

              ! Fetch the contributions of the (test) basis functions Phi_i
              ! (our "O")  for function value and first derivatives for the
              ! current DOF into HBASIy:

              HBASI1 = Dbas(IDOFE,1,ICUBP,IEL)
              HBASI2 = Dbas(IDOFE,2,ICUBP,IEL)
              HBASI3 = Dbas(IDOFE,3,ICUBP,IEL)
              HBASI4 = Dbas(IDOFE,4,ICUBP,IEL)

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

              HSUMI = HBASI2*du1loc + HBASI3*du2loc + HBASI4*du3loc

              ! Inner loop over the DOF`s j=1..indof, which corresponds to
              ! the basis function Phi_j:

              do JDOFE=1,indof

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! (out "X") for function value and first derivatives for the
                ! current DOF into HBASJy:

                HBASJ1 = Dbas(JDOFE,1,ICUBP,IEL)
                HBASJ2 = Dbas(JDOFE,2,ICUBP,IEL)
                HBASJ3 = Dbas(JDOFE,3,ICUBP,IEL)
                HBASJ4 = Dbas(JDOFE,4,ICUBP,IEL)

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
                ! we do not have to worry about that.

                HSUMJ = HBASJ2*du1loc + HBASJ3*du2loc + HBASJ4*du3loc

                ! Finally calculate the contribution to the system
                ! matrix. Depending on the configuration of ddelta,...
                ! this is:
                !
                ! AH = n~_h(u_h,phi_j,phi_i)        | nonlinear part
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

                AH = ddelta * HSUMJ*(DlocalDelta(IEL)*HSUMI+HBASI1)

                ! Weighten the calculated value AH by the cubature
                ! weight OM and add it to the local matrix. After the
                ! loop over all DOF`s is finished, each entry contains
                ! the calculated integral.

                Dentry(JDOFE,IDOFE,IEL) = Dentry(JDOFE,IDOFE,IEL) + OM*AH

              end do ! IDOFE

            end do ! JDOFE

          end do ! ICUBP

        end do ! IEL

      end if

      ! If dny != 0 or dalpha != 0, add the Laplace/Mass matrix to the
      ! local matrices.
      if ((dalpha .ne. 0.0_DP) .or. (dny .ne. 0.0_DP)) then

        ! Loop over the elements in the current set.
        do IEL=1,IELmax-IELset+1

          ! Loop over all cubature points on the current element
          do ICUBP = 1, ncubp

            ! Calculate the current weighting factor in the cubature formula
            ! in that cubature point.
            !
            ! Normally, we have to take the absolut value of the determinant
            ! of the mapping here!
            ! In 2D, the determinant is always positive, whereas in 3D,
            ! the determinant might be negative -- that is normal!
            ! But because this routine only works in 2D, we can skip
            ! the ABS here!

            OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))

            ! Current velocity in this cubature point:
            du1loc = Dvelocity (1,ICUBP,IEL)
            du2loc = Dvelocity (2,ICUBP,IEL)
            du3loc = Dvelocity (3,ICUBP,IEL)

            ! Outer loop over the DOF`s i=1..indof on our current element,
            ! which corresponds to the basis functions Phi_i:

            do IDOFE=1,indof

              ! Fetch the contributions of the (test) basis functions Phi_i
              ! (our "O")  for function value and first derivatives for the
              ! current DOF into HBASIy:

              HBASI1 = Dbas(IDOFE,1,ICUBP,IEL)
              HBASI2 = Dbas(IDOFE,2,ICUBP,IEL)
              HBASI3 = Dbas(IDOFE,3,ICUBP,IEL)
              HBASI4 = Dbas(IDOFE,4,ICUBP,IEL)

              ! Inner loop over the DOF`s j=1..indof, which corresponds to
              ! the basis function Phi_j:

              do JDOFE=1,indof

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! (out "X") for function value and first derivatives for the
                ! current DOF into HBASJy:

                HBASJ1 = Dbas(JDOFE,1,ICUBP,IEL)
                HBASJ2 = Dbas(JDOFE,2,ICUBP,IEL)
                HBASJ3 = Dbas(JDOFE,3,ICUBP,IEL)
                HBASJ4 = Dbas(JDOFE,4,ICUBP,IEL)

                ! Finally calculate the contribution to the system
                ! matrix. Depending on the configuration of DNU,
                ! dalpha,... this decomposes into:
                !
                ! AH = dny*(grad(phi_j,grad(phi_i)) | -dny*Laplace(u) = -dbeta*Stokes
                !    + dalpha*(phi_j*phi_i)         | Mass matrix

                AH = dny*(HBASI2*HBASJ2 + HBASI3*HBASJ3 + HBASI4*HBASJ4) &
                    + dalpha*HBASI1*HBASJ1

                ! Weighten the calculated value AH by the cubature
                ! weight OM and add it to the local matrix. After the
                ! loop over all DOF`s is finished, each entry contains
                ! the calculated integral.

                Dentry(JDOFE,IDOFE,IEL) = Dentry(JDOFE,IDOFE,IEL) + OM*AH

              end do ! IDOFE

            end do ! JDOFE

          end do ! ICUBP

        end do ! IEL

      end if

      ! Should we assemble the Newton matrices?
      if (dnewton .ne. 0.0_DP) then

        ! Loop over the elements in the current set.
        do IEL=1,IELmax-IELset+1

          ! Loop over all cubature points on the current element
          do ICUBP = 1, ncubp

            ! Calculate the current weighting factor in the cubature formula
            ! in that cubature point.
            !
            ! Normally, we have to take the absolut value of the determinant
            ! of the mapping here!
            OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))

            ! Current velocity in this cubature point:
            du1locx = DvelocityUderiv (1,ICUBP,IEL)
            du1locy = DvelocityUderiv (2,ICUBP,IEL)
            du1locz = DvelocityUderiv (3,ICUBP,IEL)
            du2locx = DvelocityVderiv (1,ICUBP,IEL)
            du2locy = DvelocityVderiv (2,ICUBP,IEL)
            du2locz = DvelocityVderiv (3,ICUBP,IEL)
            du3locx = DvelocityWderiv (1,ICUBP,IEL)
            du3locy = DvelocityWderiv (2,ICUBP,IEL)
            du3locz = DvelocityWderiv (3,ICUBP,IEL)

            ! Outer loop over the DOF`s i=1..indof on our current element,
            ! which corresponds to the basis functions Phi_i:

            do IDOFE=1,indof

              ! Fetch the contributions of the (test) basis functions Phi_i
              ! (our "O")  for function value and first derivatives for the
              ! current DOF into HBASIy:
              HBASI1 = Dbas(IDOFE,1,ICUBP,IEL)

              ! Inner loop over the DOF`s j=1..indof, which corresponds to
              ! the basis function Phi_j:
              do JDOFE=1,indof

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! (out "X") for function value and first derivatives for the
                ! current DOF into HBASJy:
                HBASJ1 = Dbas(JDOFE,1,ICUBP,IEL)

                ! Finally calculate the contribution to the system
                ! matrices A11, A12, A21 and A22.
                !
                ! The Newton part is calculated as follows:
                !
                ! U * grad(V)  =  ( U * grad(.) ) V
                !
                !              =  ( U * grad(V1) )
                !                 ( U * grad(V2) )
                !                 ( U * grad(V3) )
                !
                !              =  ( (U1)   (V1x) )
                !                 ( (U2) * (V1y) )
                !                 ( (U3)   (V1z) )
                !                 (              )
                !                 ( (U1)   (V2x) )
                !                 ( (U2) * (V2y) )
                !                 ( (U3)   (V2z) )
                !                 (              )
                !                 ( (U1)   (V3x) )
                !                 ( (U2) * (V3y) )
                !                 ( (U3)   (V3z) )
                !
                !              =  ( U1 * V1x  + U2 * V1y + U3 * V1z)
                !                 ( U1 * V2x  + U2 * V2y + U3 * V2z)
                !                 ( U1 * V3x  + U2 * V3y + U3 * V3z)
                !
                !              =  ( V1x  V1y  V1z ) ( U1 )
                !                 ( V2x  V2y  V2z ) ( U2 )
                !                 ( V3x  V3y  V3z ) ( U3 )
                !
                !              -> ( A11  A12  A13 )
                !                 ( A21  A22  A23 )
                !                 ( A31  A32  A33 )
                !
                ! With the velocity V=(u,v,w), we have to assemble:
                ! grad(V)*U, which is realised in each cubature point as:
                !   du/dx * phi_j*phi_i -> A11
                !   du/dy * phi_j*phi_i -> A12
                !   du/dz * phi_j*phi_i -> A13
                !   dv/dx * phi_j*phi_i -> A21
                !   dv/dy * phi_j*phi_i -> A22
                !   dv/dz * phi_j*phi_i -> A23
                !   dw/dx * phi_j*phi_i -> A31
                !   dw/dy * phi_j*phi_i -> A32
                !   dw/dz * phi_j*phi_i -> A33

                AH11 = du1locx * HBASJ1*HBASI1
                AH12 = du1locy * HBASJ1*HBASI1
                AH13 = du1locz * HBASJ1*HBASI1
                AH21 = du2locx * HBASJ1*HBASI1
                AH22 = du2locy * HBASJ1*HBASI1
                AH23 = du2locz * HBASJ1*HBASI1
                AH31 = du3locx * HBASJ1*HBASI1
                AH32 = du3locy * HBASJ1*HBASI1
                AH33 = du3locz * HBASJ1*HBASI1

                ! Weighten the calculated value AHxy by the cubature
                ! weight OM and add it to the local matrices. After the
                ! loop over all DOF`s is finished, each entry contains
                ! the calculated integral.

                DentryA11(JDOFE,IDOFE,IEL) = DentryA11(JDOFE,IDOFE,IEL)+OM*AH11
                DentryA12(JDOFE,IDOFE,IEL) = DentryA12(JDOFE,IDOFE,IEL)+OM*AH12
                DentryA13(JDOFE,IDOFE,IEL) = DentryA13(JDOFE,IDOFE,IEL)+OM*AH13
                DentryA21(JDOFE,IDOFE,IEL) = DentryA21(JDOFE,IDOFE,IEL)+OM*AH21
                DentryA22(JDOFE,IDOFE,IEL) = DentryA22(JDOFE,IDOFE,IEL)+OM*AH22
                DentryA23(JDOFE,IDOFE,IEL) = DentryA23(JDOFE,IDOFE,IEL)+OM*AH23
                DentryA31(JDOFE,IDOFE,IEL) = DentryA31(JDOFE,IDOFE,IEL)+OM*AH31
                DentryA32(JDOFE,IDOFE,IEL) = DentryA32(JDOFE,IDOFE,IEL)+OM*AH32
                DentryA33(JDOFE,IDOFE,IEL) = DentryA33(JDOFE,IDOFE,IEL)+OM*AH33

              end do ! IDOFE

            end do ! JDOFE

          end do ! ICUBP

        end do ! IEL

      end if



      ! Now we have set up "local" system matrices. We can either
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

      if (iand(cdef,CONV_MODMATRIX) .ne. 0) then

        ! With or without Newton?
        if (dnewton .eq. 0.0_DP) then

          ! Include the local matrices into the global system matrix,
          ! subblock A11 and (if different from A11) also into A22 and A33.
          !%omp critical
          do IEL=1,IELmax-IELset+1
            do IDOFE=1,indof
              do JDOFE=1,indof
                p_Da11(Kentry(JDOFE,IDOFE,IEL)) = p_Da11(Kentry(JDOFE,IDOFE,IEL)) + &
                    dtheta * Dentry(JDOFE,IDOFE,IEL)
              end do
            end do
          end do
          if (.not. associated(p_Da22,p_Da11)) then
            do IEL=1,IELmax-IELset+1
              do IDOFE=1,indof
                do JDOFE=1,indof
                  p_Da22(Kentry(JDOFE,IDOFE,IEL)) = p_Da22(Kentry(JDOFE,IDOFE,IEL)) + &
                      dtheta * Dentry(JDOFE,IDOFE,IEL)
                end do
              end do
            end do
          end if
          if ((.not. associated(p_Da33,p_Da11)) .and. &
              (.not. associated(p_Da33,p_Da22))) then
            do IEL=1,IELmax-IELset+1
              do IDOFE=1,indof
                do JDOFE=1,indof
                  p_Da33(Kentry(JDOFE,IDOFE,IEL)) = p_Da33(Kentry(JDOFE,IDOFE,IEL)) + &
                      dtheta * Dentry(JDOFE,IDOFE,IEL)
                end do
              end do
            end do
          end if
          !%omp end critical

        else

          ! Include the local matrices into the global system matrix,
          ! subblock A11 and A22 (both must exist and be independent from
          ! each other).
          !%omp critical
          do IEL=1,IELmax-IELset+1
            do IDOFE=1,indof
              do JDOFE=1,indof
                ! Kentry (:,:,:) -> positions of local matrix in A11 and A22.
                !
                ! DentryA11 (:,:,:) -> Newton part of A11
                p_Da11(Kentry(JDOFE,IDOFE,IEL)) = p_Da11(Kentry(JDOFE,IDOFE,IEL)) + &
                    dtheta * ( Dentry(JDOFE,IDOFE,IEL) + &
                               dnewton*DentryA11(JDOFE,IDOFE,IEL) )

                ! DentryA22 (:,:,:) -> Newton part of A22
                p_Da22(Kentry(JDOFE,IDOFE,IEL)) = p_Da22(Kentry(JDOFE,IDOFE,IEL)) + &
                    dtheta * ( Dentry(JDOFE,IDOFE,IEL) + &
                               dnewton*DentryA22(JDOFE,IDOFE,IEL) )

                ! Dentry12 (:,:,:) -> Newton part of A12
                p_Da12(Kentry(JDOFE,IDOFE,IEL)) = p_Da12(Kentry(JDOFE,IDOFE,IEL)) + &
                    dtheta * dnewton * DentryA12(JDOFE,IDOFE,IEL)

                ! Dentry21 (:,:,:) -> Newton part of A21
                p_Da21(Kentry(JDOFE,IDOFE,IEL)) = p_Da21(Kentry(JDOFE,IDOFE,IEL)) + &
                    dtheta * dnewton * DentryA21(JDOFE,IDOFE,IEL)

                ! Dentry13 (:,:,:) -> Newton part of A13
                p_Da13(Kentry(JDOFE,IDOFE,IEL)) = p_Da13(Kentry(JDOFE,IDOFE,IEL)) + &
                    dtheta * dnewton * DentryA13(JDOFE,IDOFE,IEL)

                ! Dentry31 (:,:,:) -> Newton part of A31
                p_Da31(Kentry(JDOFE,IDOFE,IEL)) = p_Da31(Kentry(JDOFE,IDOFE,IEL)) + &
                    dtheta * dnewton * DentryA31(JDOFE,IDOFE,IEL)

                ! Dentry23 (:,:,:) -> Newton part of A23
                p_Da23(Kentry(JDOFE,IDOFE,IEL)) = p_Da23(Kentry(JDOFE,IDOFE,IEL)) + &
                    dtheta * dnewton * DentryA23(JDOFE,IDOFE,IEL)

                ! Dentry32 (:,:,:) -> Newton part of A32
                p_Da32(Kentry(JDOFE,IDOFE,IEL)) = p_Da32(Kentry(JDOFE,IDOFE,IEL)) + &
                    dtheta * dnewton * DentryA32(JDOFE,IDOFE,IEL)
              end do
            end do
          end do
          !%omp end critical

        end if

      end if

      ! For cdef containing CONV_MODDEFECT, build the defect vector
      !     D = RHS - A*U
      ! This is done matrix free, only with the help of the local
      ! matrix.
      ! In this case, D=(D1,D2) is expected to be the RHS on
      ! entry and will be updated to be the defect vector when
      ! this routine is left.

      if (iand(cdef,CONV_MODDEFECT) .ne. 0) then

        ! With or without Newton?
        if (dnewton .eq. 0.0_DP) then
          !%omp critical
          do IEL=1,IELmax-IELset+1
            do IDOFE=1,indof

              IDFG=Idofs(IDOFE,IEL)

              do JDOFE=1,indof

                denth = dtheta*Dentry(JDOFE,IDOFE,IEL)

                JDFG=Idofs(JDOFE,IEL)
                Ddef1(IDFG)= Ddef1(IDFG) - denth*Du1(JDFG)
                Ddef2(IDFG)= Ddef2(IDFG) - denth*Du2(JDFG)
                Ddef3(IDFG)= Ddef3(IDFG) - denth*Du3(JDFG)

              end do
            end do
          end do
          !%omp end critical
        else
          !%omp critical
          do IEL=1,IELmax-IELset+1
            do IDOFE=1,indof

              IDFG=Idofs(IDOFE,IEL)

              do JDOFE=1,indof

                denth = dtheta*Dentry(JDOFE,IDOFE,IEL)

                JDFG=Idofs(JDOFE,IEL)
                Ddef1(IDFG)= Ddef1(IDFG) - denth*Du1(JDFG)
                Ddef2(IDFG)= Ddef2(IDFG) - denth*Du2(JDFG)
                Ddef3(IDFG)= Ddef3(IDFG) - denth*Du3(JDFG)

                ! Newton part
                Ddef1(IDFG)= Ddef1(IDFG) &
                           - dtheta*dnewton*DentryA11(JDOFE,IDOFE,IEL)*Du1(JDFG) &
                           - dtheta*dnewton*DentryA12(JDOFE,IDOFE,IEL)*Du2(JDFG) &
                           - dtheta*dnewton*DentryA13(JDOFE,IDOFE,IEL)*Du3(JDFG)
                Ddef2(IDFG)= Ddef2(IDFG) &
                           - dtheta*dnewton*DentryA21(JDOFE,IDOFE,IEL)*Du1(JDFG) &
                           - dtheta*dnewton*DentryA22(JDOFE,IDOFE,IEL)*Du2(JDFG) &
                           - dtheta*dnewton*DentryA23(JDOFE,IDOFE,IEL)*Du3(JDFG)
                Ddef3(IDFG)= Ddef3(IDFG) &
                           - dtheta*dnewton*DentryA31(JDOFE,IDOFE,IEL)*Du1(JDFG) &
                           - dtheta*dnewton*DentryA32(JDOFE,IDOFE,IEL)*Du2(JDFG) &
                           - dtheta*dnewton*DentryA33(JDOFE,IDOFE,IEL)*Du3(JDFG)

              end do
            end do
          end do
          !%omp end critical
        end if

      end if


    end do ! IELset
    !%omp end do

    ! Release memory
    call elprep_releaseElementSet(revalElementSet)

    deallocate(DlocalDelta)
    if (dnewton .ne. 0.0_DP) then
      deallocate(DentryA33)
      deallocate(DentryA32)
      deallocate(DentryA31)
      deallocate(DentryA23)
      deallocate(DentryA22)
      deallocate(DentryA21)
      deallocate(DentryA13)
      deallocate(DentryA12)
      deallocate(DentryA11)
      deallocate(DvelocityWderiv)
      deallocate(DvelocityVderiv)
      deallocate(DvelocityUderiv)
    end if
    deallocate(Dvelocity)
    deallocate(Dentry)
    deallocate(Kentry)
    deallocate(IdofsALE)
    deallocate(Idofs)
    deallocate(DbasALE)
    deallocate(Dbas)
    !%omp end parallel
    deallocate(Domega)
    deallocate(p_DcubPtsRef)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine getLocalDeltaHexaRay (U1L1,U1L2,U1L3,U2L1,U2L2,U2L3,&
           A1L,A2L,IEL,duMaxR,ddelta,Kvert,Dcorvg,KDFG,IDFL,UPSAM,NUREC)
!<description>
  ! This routine calculates a local ddelta=DELTA_T for a hexahedral finite
  ! element T=IEL. This can be used by the streamline diffusion stabilisation
  ! technique as a multiplier of the (local) bilinear form.
  !
  ! The effective velocity that is used for calculating the ddelta
  ! is combined by a weighted mean of the two velocity fields U1,U2
  ! by:
  !                   Ux = A1*U1Lx + A2*U2Lx
  ! The coefficients A1,A2 allow the caller to take influence on which
  ! velocity field to weight more.
!</description>

!<input>
  ! Main velocity field.
  real(DP), dimension(:), intent(in) :: U1L1,U1L2,U1L3

  ! Secondary velocity field.
  real(DP), dimension(:), intent(in) :: U2L1,U2L2,U2L3

  ! weighting factor for U1L1/U1L2
  real(DP), intent(in) :: A1L

  ! weighting factor for U2L1/U2L2
  real(DP), intent(in) :: A2L

  ! Reciprocal of the maximum norm of velocity in the domain:
  ! 1/duMaxR = 1/||u||_Omega
  real(DP), intent(in) :: duMaxR

  ! Reciprocal value 1/NU of coefficient NU in front of the
  ! Laplacian term of the Navier-Stokes equation
  !   NU * Laplace(u) + u*grad(u) + ...
  real(DP), intent(in) :: NUREC

  ! user defined parameter for configuring the streamline diffusion.
  ! < 0: Simple calculation of ddelta, using
  !      ddelta = |UPSAM| * h_T.
  ! > 0: usually UPSAM = 0.1 .. 2; Samarskji-like calculation of ddelta using:
  !      ddelta = UPSAM * h_t/||u||_T * 2*Re_T/(1+Re_T)
  real(DP), intent(in) :: UPSAM

  ! Element where the ddelta should be calculated
  integer, intent(in) :: IEL

  ! Number of degrees of freedom on element IEL
  integer, intent(in) :: IDFL

  ! Array with global degrees of freedom, corresponding to
  ! local degrees of freedom 1..IDFL on element IEL.
  integer, dimension(:), intent(in) :: KDFG

  ! The IverticesAtElement array from the triangulation
  integer, dimension(:,:), intent(in) :: Kvert

  ! The DvertexCoords array from the triangulation
  real(DP), dimension(:,:), intent(in) :: Dcorvg
!</input>

!<output>
  ! The local delta for this quadrilateral
  real(DP), intent(out) :: ddelta

!</output>
!</subroutine>

  ! local variables
  real(DP) :: dlocalH,dunorm,RELOC
  real(DP), dimension(3) :: Du
  integer :: idof

    ! Loop through the local degrees of freedom on element IEL.
    ! Sum up the velocities on these DOF`s. This will result
    ! in the vector (DU1,DU2) representing the (mean) X/Y-velocity
    ! through element IEL.

    ! For elements whose DOF`s represent directly the velocity, U1/U2
    ! represent the mean velocity
    ! along an egde/on the midpoint of each edge, so U1/U2 is
    ! clearly an approximation to the velocity in element T.
    Du = 0.0_DP
    do idof=1, IDFL
      Du(1) = Du(1) + (A1L*U1L1(KDFG(idof)) + A2L*U2L1(KDFG(idof)))
      Du(2) = Du(2) + (A1L*U1L2(KDFG(idof)) + A2L*U2L2(KDFG(idof)))
      Du(3) = Du(3) + (A1L*U1L3(KDFG(idof)) + A2L*U2L3(KDFG(idof)))
    end do

    ! Calculate the norm of that local velocity:
    dunorm = sqrt(Du(1)**2 + Du(2)**2 + Du(3)**2) / dble(IDFL)

    ! Now we have:   dunorm = ||u||_T
    ! and:           u_T = a1*u1_T + a2*u2_T

    ! If the norm of the velocity is small, we choose ddelta = 0,
    ! which results in central difference in the streamline diffusion
    ! matrix assembling:

    if (dunorm .le. 1.0E-8_DP) then

      ddelta = 0.0_DP

    else

      ! u_T defines the "slope" of the velocity through
      ! the element T. At next, calculate the local mesh width
      ! dlocalH = h = h_T on our element T=IEL:
      call getLocalMeshWidthHexa (dlocalH,Du,dunorm,IEL,Kvert,Dcorvg)

      ! Calculate ddelta... (cf. p. 121 in Turek`s CFD book)

      if (UPSAM.lt.0.0_DP) then

        ! For UPSAM<0, we use simple calculation of ddelta:

        ddelta = abs(UPSAM)*dlocalH

      else

        ! For UPSAM >= 0, we use standard Samarskji-like calculation
        ! of ddelta. At first calculate the local Reynolds number
        ! RELOC = Re_T = ||u||_T * h_T / NU

        RELOC = dunorm*dlocalH*NUREC

        ! and then the ddelta = UPSAM * h_t/||u|| * 2*Re_T/(1+Re_T)

        ddelta = UPSAM * dlocalH*duMaxR * 2.0_DP*(RELOC/(1.0_DP+RELOC))

      end if ! (UPSAM.LT.0.0)

    end if ! (dunorm.LE.1D-8)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine getLocalDeltaHexaVol (U1L1,U1L2,U1L3,U2L1,U2L2,U2L3,&
           A1L,A2L,IEL,duMaxR,ddelta,Kvert,Dcorvg,KDFG,IDFL,UPSAM,NUREC)
!<description>
  ! This routine calculates a local ddelta=DELTA_T for a hexahedral finite
  ! element T=IEL. This can be used by the streamline diffusion stabilisation
  ! technique as a multiplier of the (local) bilinear form.
  !
  ! The effective velocity that is used for calculating the ddelta
  ! is combined by a weighted mean of the two velocity fields U1,U2
  ! by:
  !                   Ux = A1*U1Lx + A2*U2Lx
  ! The coefficients A1,A2 allow the caller to take influence on which
  ! velocity field to weight more.
!</description>

!<input>
  ! Main velocity field.
  real(DP), dimension(:), intent(in) :: U1L1,U1L2,U1L3

  ! Secondary velocity field.
  real(DP), dimension(:), intent(in) :: U2L1,U2L2,U2L3

  ! weighting factor for U1L1/U1L2
  real(DP), intent(in) :: A1L

  ! weighting factor for U2L1/U2L2
  real(DP), intent(in) :: A2L

  ! Reciprocal of the maximum norm of velocity in the domain:
  ! 1/duMaxR = 1/||u||_Omega
  real(DP), intent(in) :: duMaxR

  ! Reciprocal value 1/NU of coefficient NU in front of the
  ! Laplacian term of the Navier-Stokes equation
  !   NU * Laplace(u) + u*grad(u) + ...
  real(DP), intent(in) :: NUREC

  ! user defined parameter for configuring the streamline diffusion.
  ! < 0: Simple calculation of ddelta, using
  !      ddelta = |UPSAM| * h_T.
  ! > 0: usually UPSAM = 0.1 .. 2; Samarskji-like calculation of ddelta using:
  !      ddelta = UPSAM * h_t/||u||_T * 2*Re_T/(1+Re_T)
  real(DP), intent(in) :: UPSAM

  ! Element where the ddelta should be calculated
  integer, intent(in) :: IEL

  ! Number of degrees of freedom on element IEL
  integer, intent(in) :: IDFL

  ! Array with global degrees of freedom, corresponding to
  ! local degrees of freedom 1..IDFL on element IEL.
  integer, dimension(:), intent(in) :: KDFG

  ! The IverticesAtElement array from the triangulation
  integer, dimension(:,:), intent(in) :: Kvert

  ! The DvertexCoords array from the triangulation
  real(DP), dimension(:,:), intent(in) :: Dcorvg
!</input>

!<output>
  ! The local delta for this quadrilateral
  real(DP), intent(out) :: ddelta

!</output>
!</subroutine>

  ! local variables
  real(DP) :: dlocalH,dunorm,RELOC
  real(DP), dimension(3) :: Du
  integer :: idof

    ! Calculate the local mesh width dlocalH = h = h_T on our element T=IEL:
    call getHexaVolume(dlocalH,IEL,Kvert,Dcorvg)
    dlocalH = dlocalH**(1.0_DP / 3.0_DP)

    ! Calculate ddelta... (cf. p. 121 in Turek`s CFD book)
    if (UPSAM.lt.0.0_DP) then

      ! For UPSAM<0, we use simple calculation of ddelta:
      ddelta = abs(UPSAM)*dlocalH

    else

      ! Loop through the local degrees of freedom on element IEL.
      ! Sum up the velocities on these DOF`s. This will result
      ! in the vector (DU1,DU2) representing the (mean) X/Y-velocity
      ! through element IEL.

      ! For elements whose DOF`s represent directly the velocity, U1/U2
      ! represent the mean velocity
      ! along an egde/on the midpoint of each edge, so U1/U2 is
      ! clearly an approximation to the velocity in element T.
      Du = 0.0_DP
      do idof=1, IDFL
        Du(1) = Du(1) + (A1L*U1L1(KDFG(idof)) + A2L*U2L1(KDFG(idof)))
        Du(2) = Du(2) + (A1L*U1L2(KDFG(idof)) + A2L*U2L2(KDFG(idof)))
        Du(3) = Du(3) + (A1L*U1L3(KDFG(idof)) + A2L*U2L3(KDFG(idof)))
      end do

      ! Calculate the norm of that local velocity:
      dunorm = sqrt(Du(1)**2 + Du(2)**2 + Du(3)**2) / dble(IDFL)

      ! For UPSAM >= 0, we use standard Samarskji-like calculation
      ! of ddelta. At first calculate the local Reynolds number
      ! RELOC = Re_T = ||u||_T * h_T / NU
      RELOC = dunorm*dlocalH*NUREC

      ! and then the ddelta = UPSAM * h_t/||u|| * 2*Re_T/(1+Re_T)
      ddelta = UPSAM * dlocalH*duMaxR * 2.0_DP*(RELOC/(1.0_DP+RELOC))

    end if ! (UPSAM.LT.0.0)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine getLocalMeshWidthHexa (dlocalH,Du,dunorm,iel,&
                                         IverticesAtElement,DvertexCoords)

!<description>
  ! Determine the local mesh width for a hexahedral element iel of a
  ! triangulation.
!</description>

!<input>
  ! mean velocity u_T through element T=iel
  real(DP), dimension(3), intent(in) :: Du

  ! norm ||u||_T = mean velocity through element T=JEL
  real(DP), intent(in) :: dunorm

  ! Element where the local h should be calculated
  integer, intent(in) :: iel

  ! The IverticesAtElement array from the triangulation
  integer, dimension(:,:), intent(in) :: IverticesAtElement

  ! The DvertexCoords array from the triangulation
  real(DP), dimension(:,:), intent(in) :: DvertexCoords
!</input>

!<output>
  ! The local mesh width for the hexahedron
  real(DP), intent(out) :: dlocalH
!</output>
!</subroutine>

  ! local variables
  integer :: i,isrc
  real(DP) :: ds,dt,dmin,dalpha
  real(DP), dimension(3,8) :: Dv
  real(DP), dimension(3,6) :: Dtan1,Dtan2,Dnormal,Dmid
  real(DP), dimension(3) :: Dray

  ! Description of the local Mesh width calculation
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! Basically, this method corresponds to the mesh width calculation
  ! in the 2D quadrilateral case, however, as the faces of a hexhedron
  ! are not necessarily planes, this gets a bit trickier.
  ! To avoid calculating intersections with manifolds we will approximate
  ! the "real" hexahedron by a hexahedron whose faces are planes, and
  ! calculate the intersections with that approximating hexahedron.
  !
  ! For each face of the "real" hexahedron, we will calculate 4 vectors:
  ! 1. the face midpoint
  ! 2. two vectors approximating the tangential space of the face midpoint
  ! 3. the outer normal vector of the tangential space
  !
  ! Imagine that this is a face of our "real" hexahedron:
  !
  !                        4------c------3
  !                        |      ^      |
  !                        |      |      |
  !                        d--t1--+----->b
  !                        |      |      |
  !                        |      t2     |
  !                        1------a------2
  !
  ! The face midpoint (+) will be calculated as the average of the face`s
  ! four corner vertices.
  ! The tangential vectors will be approximated by the vector between
  ! the midpoints of two opposite edges, i.e. t1 := b-d, t2 := c-a.
  ! The outer normal vector will be calculated using the 3D cross product,
  ! so in the image above it is pointing right between your eyes... ^_^

    ! Normalise the mean velocity to get a ray vector
    dt = sqrt(Du(1)**2 + Du(2)**2 + Du(3)**2)
    Dray = Du / dt

    ! Get the coordinates of our eight corner vertices
    Dv(1:3,1) = DvertexCoords(1:3,IverticesAtElement(1,iel))
    Dv(1:3,2) = DvertexCoords(1:3,IverticesAtElement(2,iel))
    Dv(1:3,3) = DvertexCoords(1:3,IverticesAtElement(3,iel))
    Dv(1:3,4) = DvertexCoords(1:3,IverticesAtElement(4,iel))
    Dv(1:3,5) = DvertexCoords(1:3,IverticesAtElement(5,iel))
    Dv(1:3,6) = DvertexCoords(1:3,IverticesAtElement(6,iel))
    Dv(1:3,7) = DvertexCoords(1:3,IverticesAtElement(7,iel))
    Dv(1:3,8) = DvertexCoords(1:3,IverticesAtElement(8,iel))

    ! Calculate the midpoints of the six faces
    Dmid(1:3,1) = 0.25_DP * (Dv(1:3,4)+Dv(1:3,3)+Dv(1:3,2)+Dv(1:3,1))
    Dmid(1:3,2) = 0.25_DP * (Dv(1:3,1)+Dv(1:3,2)+Dv(1:3,5)+Dv(1:3,6))
    Dmid(1:3,3) = 0.25_DP * (Dv(1:3,2)+Dv(1:3,3)+Dv(1:3,7)+Dv(1:3,6))
    Dmid(1:3,4) = 0.25_DP * (Dv(1:3,3)+Dv(1:3,4)+Dv(1:3,8)+Dv(1:3,7))
    Dmid(1:3,5) = 0.25_DP * (Dv(1:3,4)+Dv(1:3,1)+Dv(1:3,5)+Dv(1:3,8))
    Dmid(1:3,6) = 0.25_DP * (Dv(1:3,5)+Dv(1:3,6)+Dv(1:3,7)+Dv(1:3,8))

    ! Calculate the tangential vectors of the six faces
    Dtan1(1:3,1) = (Dv(1:3,1)+Dv(1:3,4)) - (Dv(1:3,2)+Dv(1:3,3))
    Dtan1(1:3,2) = (Dv(1:3,2)+Dv(1:3,6)) - (Dv(1:3,1)+Dv(1:3,5))
    Dtan1(1:3,3) = (Dv(1:3,3)+Dv(1:3,7)) - (Dv(1:3,2)+Dv(1:3,6))
    Dtan1(1:3,4) = (Dv(1:3,4)+Dv(1:3,8)) - (Dv(1:3,3)+Dv(1:3,7))
    Dtan1(1:3,5) = (Dv(1:3,1)+Dv(1:3,5)) - (Dv(1:3,4)+Dv(1:3,8))
    Dtan1(1:3,6) = (Dv(1:3,6)+Dv(1:3,7)) - (Dv(1:3,5)+Dv(1:3,8))
    Dtan2(1:3,1) = (Dv(1:3,3)+Dv(1:3,4)) - (Dv(1:3,1)+Dv(1:3,2))
    Dtan2(1:3,2) = (Dv(1:3,5)+Dv(1:3,6)) - (Dv(1:3,1)+Dv(1:3,2))
    Dtan2(1:3,3) = (Dv(1:3,6)+Dv(1:3,7)) - (Dv(1:3,2)+Dv(1:3,3))
    Dtan2(1:3,4) = (Dv(1:3,7)+Dv(1:3,8)) - (Dv(1:3,3)+Dv(1:3,4))
    Dtan2(1:3,5) = (Dv(1:3,8)+Dv(1:3,5)) - (Dv(1:3,4)+Dv(1:3,1))
    Dtan2(1:3,6) = (Dv(1:3,7)+Dv(1:3,8)) - (Dv(1:3,5)+Dv(1:3,6))

    ! Go through all faces and calculate the normals
    do i = 1,6
      ! Calculate outer normal vector by 3D cross product
      Dnormal(1,i) = Dtan1(2,i)*Dtan2(3,i) - Dtan1(3,i)*Dtan2(2,i)
      Dnormal(2,i) = Dtan1(3,i)*Dtan2(1,i) - Dtan1(1,i)*Dtan2(3,i)
      Dnormal(3,i) = Dtan1(1,i)*Dtan2(2,i) - Dtan1(2,i)*Dtan2(1,i)

      ! And normalise it
      dt = 1.0_DP / sqrt(Dnormal(1,i)**2 + Dnormal(2,i)**2 + Dnormal(3,i)**2)
      Dnormal(:,i) = Dnormal(:,i) * dt
    end do

    ! Now we have a little problem:
    ! Since the refinement algorithm for hexahedron cells produces cells
    ! which are flipped, it may (and will) happen that the normal vectors
    ! we have calculated are inner normal vectors instead of outer ones.
    ! So we need to check whether the hexahedron is flipped.
    if (gaux_isFlipped_hexa3D(Dv)) Dnormal = -Dnormal

    ! Now we have calculated all the vectors we need - we have the face
    ! midpoints, normal vectors and a ray vector.
    ! The next step is to choose a "source face", i.e. the face at whose
    ! midpoint our ray should start. As both the ray vector and the normal
    ! vectors are normalised, the scalar product of them is in range [-1,1].
    !
    ! So, let r be the ray vector and n_1,...,n_6 be the six normal vectors.
    ! Our "source face" will be the face whose normal vector points
    ! at most into the opposite direction of the ray vector, i.e.:
    !
    !          ( r, n_isrc ) <= ( r, n_i )     for all i=1,...,6
    isrc = 0
    dmin = 10.0_DP
    do i = 1,6
      ! Calculate scalar product of ray vector and outer normal vector
      dt = Dray(1)*Dnormal(1,i) + Dray(2)*Dnormal(2,i) + Dray(3)*Dnormal(3,i)

      ! Check if this is the "source face"
      if (dt .lt. dmin) then
        isrc = i
        dmin = dt
      end if
    end do

    ! Now we have chosen a source face, so we now need to calculate the
    ! intersection of the ray vector starting at the source face with
    ! all the other five faces. Our local mesh width will be defined
    ! as the minimum of all s_i, where s_i is the ray scaling factor
    ! such that:
    !
    !                    m + s_i*r \in f_i   for i=1,...,6
    !
    ! where:
    ! m is the face midpoint of the "source face"
    ! r is the normalised ray vector
    ! f_i is the i-th face of the approximating hexahedron
    !
    ! The ray scaling factor for a given face i is calculated as:
    !
    !         - ( n_i, r )
    ! s_i := ----------------
    !        ( n_i, m_i - m )
    !
    ! where:
    ! r is the normalised ray vector
    ! m is the face midpoint of the "source face"
    ! m_i is the midpoint of the i-th face
    ! n_i is the normal vector of the i-th face

    ! Go through all faces again
    dalpha = 1.0E99_DP
    do i = 1, 6

      ! If this is the "source face" then skip it
      if (i .eq. isrc) cycle

      ! Now compute the scalar product of the ray and the outer normal
      ! of this face.
      dt = Dray(1)*Dnormal(1,i) + Dray(2)*Dnormal(2,i) + Dray(3)*Dnormal(3,i)

      ! If dt = 0, then the ray is parallel to the face, and if dt < 0,
      ! then the ray points away from the face - in either case, we can
      ! skip this face.
      ! Note: As we need to divide by dt later, we will check it against
      !       machine exactness instead of 0.
      if (dt .le. SYS_EPSREAL_DP) cycle

      ! Now calculate the scalar product of the face normal and the vector
      ! between the face midpoint and the "source face" midpoint:
      ds = Dnormal(1,i)*(Dmid(1,isrc) - Dmid(1,i)) &
         + Dnormal(2,i)*(Dmid(2,isrc) - Dmid(2,i)) &
         + Dnormal(3,i)*(Dmid(3,isrc) - Dmid(3,i))

      ! Now divide ds by -dt to get the ray scaling factor
      ds = -ds / dt

      ! Now ds is always positive
      dalpha = min(dalpha,ds)

    end do

    ! Now if alpha is in a valid range, then use it as the local H
    if ((dalpha .gt. 0.0_DP) .and. (dalpha .lt. 1.0E10_DP)) then
      dlocalH = dalpha
    else
      dlocalH = 0.0_DP
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine getHexaVolume (dlocalH,iel,IverticesAtElement,DvertexCoords)

!<description>
  ! Calculates the volume of a hexahedron.
!</description>

!<input>
  ! Element where the local h should be calculated
  integer, intent(in) :: iel

  ! The IverticesAtElement array from the triangulation
  integer, dimension(:,:), intent(in) :: IverticesAtElement

  ! The DvertexCoords array from the triangulation
  real(DP), dimension(:,:), intent(in) :: DvertexCoords
!</input>

!<output>
  ! The local mesh width for the hexahedron
  real(DP), intent(out) :: dlocalH
!</output>
!</subroutine>

  ! local variables
  real(DP), dimension(3,8) :: Dv

    ! Get the coordinates of our eight corner vertices
    Dv(1:3,1) = DvertexCoords(1:3,IverticesAtElement(1,iel))
    Dv(1:3,2) = DvertexCoords(1:3,IverticesAtElement(2,iel))
    Dv(1:3,3) = DvertexCoords(1:3,IverticesAtElement(3,iel))
    Dv(1:3,4) = DvertexCoords(1:3,IverticesAtElement(4,iel))
    Dv(1:3,5) = DvertexCoords(1:3,IverticesAtElement(5,iel))
    Dv(1:3,6) = DvertexCoords(1:3,IverticesAtElement(6,iel))
    Dv(1:3,7) = DvertexCoords(1:3,IverticesAtElement(7,iel))
    Dv(1:3,8) = DvertexCoords(1:3,IverticesAtElement(8,iel))

    ! Reset the volume
    dlocalH = gaux_getVolume_hexa3D(Dv)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine conv_JumpStabilisation1d ( &
      rconfig, cdef, rmatrix, rsolution, rdefect, rdiscretisation, rperfconfig)

!<description>
  ! Edge oriented stabilisation technique. Wrapper routine.
  ! 1D-version.
  !
  ! rvecPrimary, rvecSecondary are two velocity field vectors for the X-
  ! and Y-veclocity; IvelocityComp defines which components of these
  ! vectors contains the X- and which contains the Y-velocity.
  ! The final velocity vector field is then computed as a weighted average
  ! of these two:
  ! <tex>
  !  $$ u_1  =  dprimWeight * rvecPrimary  +  dsecWeight * rvecSecondary $$
  ! </tex>
  ! <tex>$ u_2 = rsolution(.) $</tex> defines the second velocity field inside of
  ! the grad-term.
  !
  ! The switch cdef decides on whether the routine sets up the nonlinear
  ! defect, the nonlinear matrix or both.
  !
  ! The configuration how the routine should react is to be configured
  ! in the configuration block rconfig.
!</description>

!<input>
  ! Configuration block for the streamline diffusion scheme
  type(t_jumpStabilisation), intent(in) :: rconfig

  ! Computation/defect correction method. One of the CONV_MODxxxx constants:
  ! CONV_MODMATRIX: Set up the nonlinear matrix. rmatrix must be present, the
  !                 nonlinear part is added to the matrix.
  ! CONV_MODDEFECT: Set up the nonlinear defect. rdefect and rsolution must be
  !                 present.
  ! CONV_MODBOTH  : Set up the nonlinear matrix as well as the nonlinear defect.
  !                 rmatrix, rdefect and rsolution must all be present.
  integer, intent(in) :: cdef

  ! OPTIONAL: Solution vector u_2.
  ! Must be present if cdef=CONV_MODDEFECT or =CONV_MODBOTH.
  type(t_vectorBlock), intent(in), target, optional :: rsolution

  ! OPTIONAL: Alternative discretisation structure to use for setting up
  ! the jump stabilisaton. This allows to use a different FE pair for
  ! setting up the stabilisation than the matrix itself.
  type(t_spatialDiscretisation), intent(in), optional :: rdiscretisation

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>

!<inputoutput>
  ! System matrix.
  ! The content of the matrix must be present if cdef=CONV_MODMATRIX or
  ! =CONV_MODBOTH, otherwise only the structure is used.
  ! The nonlinear operator is added to the matrix.
  type(t_matrixScalar), intent(inout) :: rmatrix

  ! OPTIONAL: Defect vector.
  ! Must have the same structure as rsolution/rvecPrimary/rvecSecondary.
  ! Must be present if cdef=CONV_MODDEFECT or =CONV_MODBOTH.
  ! The nonlinear part is subtracted from this vector:
  ! <tex>$ r = r - \theta * u_1*grad(u_2) $</tex>
  type(t_vectorBlock), intent(inout), optional, target :: rdefect
!</inputoutput>

!</subroutine>

    ! At first check the input parameters that everything is present what
    ! we need:
    if ((cdef .eq. CONV_MODDEFECT) .or. (cdef .eq. CONV_MODBOTH)) then
!      if ((.not. present(rsolution)) .or. (.not. present(rdefect))) then
!        call output_line ('Solution/defect vector not present', &
!            OU_CLASS_ERROR,OU_MODE_STD,'conv_JumpStabilisation3d')
!        call sys_halt()
!      end if
      call output_line ('Defect modification currently not supported', &
          OU_CLASS_ERROR,OU_MODE_STD,'conv_JumpStabilisation1d')
    end if

    ! At the moment, we only support a rather limited set of configurations:
    ! Matrix and vectors must all be double precision, matrix must be format
    ! 7 or 9, discretisation must be Q1~, constant viscosity.
    if ((rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
        (rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
      call output_line ('Unsupported matrix format', &
          OU_CLASS_ERROR,OU_MODE_STD,'conv_JumpStabilisation1d')
      call sys_halt()
    end if

    if (rmatrix%p_rspatialDiscrTest%ccomplexity .ne. SPDISC_UNIFORM) then
      call output_line ('Unsupported discretisation', &
          OU_CLASS_ERROR,OU_MODE_STD,'conv_JumpStabilisation1d')
      call sys_halt()
    end if

    !if ((rvecPrimary%cdataType .NE. ST_DOUBLE) .OR. &
    !    (rvecSecondary%cdataType .NE. ST_DOUBLE)) then
    !  PRINT *,'EOS: Unsupported vector data type in velocity.'
    !  call sys_halt()
    !end if
    !
    !if (PRESENT(rdefect)) then
    !  if ((rsolution%cdataType .NE. ST_DOUBLE) .OR. &
    !      (rdefect%cdataType .NE. ST_DOUBLE)) then
    !    PRINT *,'EOS: Unsupported vector data type in solution/defect'
    !    call sys_halt()
    !  end if
    !end if

    if (.not. rconfig%bconstViscosity) then
      call output_line ('Only constant viscosity supported at the moment', &
          OU_CLASS_ERROR,OU_MODE_STD,'conv_JumpStabilisation1d')
      call sys_halt()
    end if

    if (rconfig%dnu .eq. SYS_INFINITY_DP) then
      call output_line ('Viscosity parameter nu not initialised', &
          OU_CLASS_ERROR,OU_MODE_STD,'conv_JumpStabilisation1d')
      call sys_halt()
    end if

    if (rconfig%cjump .eq. CONV_JUMP_UNIFIEDEDGE) then
!      if (present(rdefect)) then
!
!        ! Modify the defect?
!        if (iand(cdef,CONV_MODDEFECT) .ne. 0) then
!          call jstab_matvecUEOJumpStabilBlk3d ( &
!              rconfig%dgamma,rconfig%dgammastar,rconfig%ccubType,rconfig%dnu,&
!              rmatrix,rsolution,rdefect,-rconfig%dtheta,1.0_DP,&
!              rdiscretisation)
!        end if
!
!      end if

      ! Modify the matrix?
      if (iand(cdef,CONV_MODMATRIX) .ne. 0) then
        call jstab_calcUEOJumpStabilisation (&
          rmatrix,rconfig%dgamma,rconfig%dgammastar,rconfig%deojEdgeExp,&
          rconfig%dtheta,rconfig%ccubType,rconfig%dnu,rdiscretisation,&
          rperfconfig=rperfconfig)
      end if

!    else if (rconfig%cjump .eq. CONV_JUMP_REACTIVE) then
!      if (present(rdefect)) then
!
!        ! Modify the defect?
!        if (iand(cdef,CONV_MODDEFECT) .ne. 0) then
!          call jstab_matvecReacJumpStabilBlk3d ( &
!              rconfig%dgamma,rconfig%ccubType,rconfig%dnu,&
!              rmatrix,rsolution,rdefect,-rconfig%dtheta,1.0_DP,&
!              rdiscretisation)
!        end if
!
!      end if
!
!      ! Modify the matrix?
!      if (iand(cdef,CONV_MODMATRIX) .ne. 0) then
!        call jstab_calcReacJumpStabilisation (&
!          rmatrix,rconfig%dgamma,rconfig%dtheta,&
!          rconfig%ccubType,rconfig%dnu,rdiscretisation)
!      end if

    else
      call output_line ('Unknown jump stabilisation', &
          OU_CLASS_ERROR,OU_MODE_STD,'conv_JumpStabilisation1d')
      call sys_halt()
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine conv_JumpStabilisation2d ( &
      rconfig, cdef, rmatrix, rsolution, rdefect, rdiscretisation,&
      InodeList, rperfconfig)

!<description>
  ! Edge oriented stabilisation technique. Wrapper routine.
  ! 2D-version (X- and Y-velocity).
  !
  ! rvecPrimary, rvecSecondary are two velocity field vectors for the X-
  ! and Y-veclocity; IvelocityComp defines which components of these
  ! vectors contains the X- and which contains the Y-velocity.
  ! The final velocity vector field is then computed as a weighted average
  ! of these two:
  ! <tex>
  !  $$ u_1  =  dprimWeight * rvecPrimary  +  dsecWeight * rvecSecondary $$
  ! </tex>
  ! <tex>$ u_2 = rsolution(.) $</tex> defines the second velocity field inside of
  ! the grad-term.
  !
  ! The switch cdef decides on whether the routine sets up the nonlinear
  ! defect, the nonlinear matrix or both.
  !
  ! The configuration how the routine should react is to be configured
  ! in the configuration block rconfig.
!</description>

!<input>
  ! Configuration block for the streamline diffusion scheme
  type(t_jumpStabilisation), intent(in) :: rconfig

  ! Computation/defect correction method. One of the CONV_MODxxxx constants:
  ! CONV_MODMATRIX: Set up the nonlinear matrix. rmatrix must be present, the
  !                 nonlinear part is added to the matrix.
  ! CONV_MODDEFECT: Set up the nonlinear defect. rdefect and rsolution must be
  !                 present.
  ! CONV_MODBOTH  : Set up the nonlinear matrix as well as the nonlinear defect.
  !                 rmatrix, rdefect and rsolution must all be present.
  integer, intent(in) :: cdef

  ! OPTIONAL: Solution vector u_2.
  ! Must be present if cdef=CONV_MODDEFECT or =CONV_MODBOTH.
  type(t_vectorBlock), intent(in), target, optional :: rsolution

  ! OPTIONAL: Alternative discretisation structure to use for setting up
  ! the jump stabilisaton. This allows to use a different FE pair for
  ! setting up the stabilisation than the matrix itself.
  type(t_spatialDiscretisation), intent(in), optional :: rdiscretisation

  ! OPTIONAL: List of edges/faces where the operator should be computed.
  ! If not present, the operator will be computed on all edges/faces.
  integer, dimension(:), intent(in), optional :: InodeList

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>

!<inputoutput>
  ! System matrix.
  ! The content of the matrix must be present if cdef=CONV_MODMATRIX or
  ! =CONV_MODBOTH, otherwise only the structure is used.
  ! The nonlinear operator is added to the matrix.
  type(t_matrixScalar), intent(inout) :: rmatrix

  ! OPTIONAL: Defect vector.
  ! Must have the same structure as rsolution/rvecPrimary/rvecSecondary.
  ! Must be present if cdef=CONV_MODDEFECT or =CONV_MODBOTH.
  ! The nonlinear part is subtracted from this vector:
  ! <tex>$ r = r - \theta * u_1*grad(u_2) $</tex>
  type(t_vectorBlock), intent(inout), optional, target :: rdefect
!</inputoutput>

!</subroutine>

    ! At first check the input parameters that everything is present what
    ! we need:
    if ((cdef .eq. CONV_MODDEFECT) .or. (cdef .eq. CONV_MODBOTH)) then
      if ((.not. present(rsolution)) .or. (.not. present(rdefect))) then
        call output_line ("Solution/defect vector not present!", &
            OU_CLASS_ERROR,OU_MODE_STD,"conv_JumpStabilisation2d")
        call sys_halt()
      end if
    end if

    ! At the moment, we only support a rather limited set of configurations:
    ! Matrix and vectors must all be double precision, matrix must be format
    ! 7 or 9, discretisation must be Q1~, constant viscosity.
    if ((rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
        (rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
      call output_line ("Unsupported matrix format!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_JumpStabilisation2d")
      call sys_halt()
    end if

    if (rmatrix%p_rspatialDiscrTest%ccomplexity .ne. SPDISC_UNIFORM) then
      call output_line ("Unsupported discretisation!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_JumpStabilisation2d")
      call sys_halt()
    end if

    !if ((rvecPrimary%cdataType .NE. ST_DOUBLE) .OR. &
    !    (rvecSecondary%cdataType .NE. ST_DOUBLE)) then
    !  PRINT *,'EOS: Unsupported vector data type in velocity.'
    !  call sys_halt()
    !end if
    !
    !if (PRESENT(rdefect)) then
    !  if ((rsolution%cdataType .NE. ST_DOUBLE) .OR. &
    !      (rdefect%cdataType .NE. ST_DOUBLE)) then
    !    PRINT *,'EOS: Unsupported vector data type in solution/defect'
    !    call sys_halt()
    !  end if
    !end if

    if (.not. rconfig%bconstViscosity) then
      call output_line ("Only constant viscosity supported at the moment!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_JumpStabilisation2d")
      call sys_halt()
    end if

    if (rconfig%dnu .eq. SYS_INFINITY_DP) then
      call output_line ("Viscosity parameter nu not initialised!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_JumpStabilisation2d")
      call sys_halt()
    end if

    if (rconfig%cjump .eq. CONV_JUMP_UNIFIEDEDGE) then
      if (present(rdefect)) then

        ! Modify the defect?
        if (iand(cdef,CONV_MODDEFECT) .ne. 0) then
          call jstab_matvecUEOJumpStabilBlk2d ( &
              rconfig%dgamma,rconfig%dgammastar,rconfig%deojEdgeExp,rconfig%ccubType,&
              rconfig%dnu,rmatrix,rsolution,rdefect,-rconfig%dtheta,1.0_DP,&
              rdiscretisation,InodeList,rperfconfig)
        end if

      end if

      ! Modify the matrix?
      if (iand(cdef,CONV_MODMATRIX) .ne. 0) then
        call jstab_calcUEOJumpStabilisation (&
          rmatrix,rconfig%dgamma,rconfig%dgammastar,rconfig%deojEdgeExp,rconfig%dtheta,&
          rconfig%ccubType,rconfig%dnu,rdiscretisation,InodeList,rperfconfig)
      end if

    else if (rconfig%cjump .eq. CONV_JUMP_REACTIVE) then
      if (present(rdefect)) then

        ! Modify the defect?
        if (iand(cdef,CONV_MODDEFECT) .ne. 0) then
          call jstab_matvecReacJumpStabilBlk2d ( &
              rconfig%dgamma,rconfig%ccubType,rconfig%dnu,&
              rmatrix,rsolution,rdefect,-rconfig%dtheta,1.0_DP,&
              rdiscretisation,rperfconfig)
        end if

      end if

      ! Modify the matrix?
      if (iand(cdef,CONV_MODMATRIX) .ne. 0) then
        call jstab_calcReacJumpStabilisation (&
          rmatrix,rconfig%dgamma,rconfig%dtheta,&
          rconfig%ccubType,rconfig%dnu,rdiscretisation, rperfconfig)
      end if

    else
      call output_line ("Unknown jump stabilisation!", &
          OU_CLASS_ERROR,OU_MODE_STD,"conv_JumpStabilisation2d")
      call sys_halt()
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine conv_JumpStabilisation3d ( &
      rconfig, cdef, rmatrix, rsolution, rdefect, rdiscretisation, rperfconfig)

!<description>
  ! Edge oriented stabilisation technique. Wrapper routine.
  ! 3D-version (X-, Y- and Z-velocity).
  !
  ! rvecPrimary, rvecSecondary are two velocity field vectors for the X-
  ! and Y-veclocity; IvelocityComp defines which components of these
  ! vectors contains the X- and which contains the Y-velocity.
  ! The final velocity vector field is then computed as a weighted average
  ! of these two:
  ! <tex>
  !  $$ u_1  =  dprimWeight * rvecPrimary  +  dsecWeight * rvecSecondary $$
  ! </tex>
  ! <tex>$ u_2 = rsolution(.) $</tex> defines the second velocity field inside of
  ! the grad-term.
  !
  ! The switch cdef decides on whether the routine sets up the nonlinear
  ! defect, the nonlinear matrix or both.
  !
  ! The configuration how the routine should react is to be configured
  ! in the configuration block rconfig.
!</description>

!<input>
  ! Configuration block for the streamline diffusion scheme
  type(t_jumpStabilisation), intent(in) :: rconfig

  ! Computation/defect correction method. One of the CONV_MODxxxx constants:
  ! CONV_MODMATRIX: Set up the nonlinear matrix. rmatrix must be present, the
  !                 nonlinear part is added to the matrix.
  ! CONV_MODDEFECT: Set up the nonlinear defect. rdefect and rsolution must be
  !                 present.
  ! CONV_MODBOTH  : Set up the nonlinear matrix as well as the nonlinear defect.
  !                 rmatrix, rdefect and rsolution must all be present.
  integer, intent(in) :: cdef

  ! OPTIONAL: Solution vector u_2.
  ! Must be present if cdef=CONV_MODDEFECT or =CONV_MODBOTH.
  type(t_vectorBlock), intent(in), target, optional :: rsolution

  ! OPTIONAL: Alternative discretisation structure to use for setting up
  ! the jump stabilisaton. This allows to use a different FE pair for
  ! setting up the stabilisation than the matrix itself.
  type(t_spatialDiscretisation), intent(in), optional :: rdiscretisation

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>

!<inputoutput>
  ! System matrix.
  ! The content of the matrix must be present if cdef=CONV_MODMATRIX or
  ! =CONV_MODBOTH, otherwise only the structure is used.
  ! The nonlinear operator is added to the matrix.
  type(t_matrixScalar), intent(inout) :: rmatrix

  ! OPTIONAL: Defect vector.
  ! Must have the same structure as rsolution/rvecPrimary/rvecSecondary.
  ! Must be present if cdef=CONV_MODDEFECT or =CONV_MODBOTH.
  ! The nonlinear part is subtracted from this vector:
  ! <tex>$ r = r - \theta * u_1*grad(u_2) $</tex>
  type(t_vectorBlock), intent(inout), optional, target :: rdefect
!</inputoutput>

!</subroutine>

    ! At first check the input parameters that everything is present what
    ! we need:
    if ((cdef .eq. CONV_MODDEFECT) .or. (cdef .eq. CONV_MODBOTH)) then
!      if ((.not. present(rsolution)) .or. (.not. present(rdefect))) then
!        call output_line ('Solution/defect vector not present', &
!            OU_CLASS_ERROR,OU_MODE_STD,'conv_JumpStabilisation3d')
!        call sys_halt()
!      end if
      call output_line ('Defect modification currently not supported', &
          OU_CLASS_ERROR,OU_MODE_STD,'conv_JumpStabilisation3d')
    end if

    ! At the moment, we only support a rather limited set of configurations:
    ! Matrix and vectors must all be double precision, matrix must be format
    ! 7 or 9, discretisation must be Q1~, constant viscosity.
    if ((rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX9) .and. &
        (rmatrix%cmatrixFormat .ne. LSYSSC_MATRIX7)) then
      call output_line ('Unsupported matrix format', &
          OU_CLASS_ERROR,OU_MODE_STD,'conv_JumpStabilisation3d')
      call sys_halt()
    end if

    if (rmatrix%p_rspatialDiscrTest%ccomplexity .ne. SPDISC_UNIFORM) then
      call output_line ('Unsupported discretisation', &
          OU_CLASS_ERROR,OU_MODE_STD,'conv_JumpStabilisation3d')
      call sys_halt()
    end if

    !if ((rvecPrimary%cdataType .NE. ST_DOUBLE) .OR. &
    !    (rvecSecondary%cdataType .NE. ST_DOUBLE)) then
    !  PRINT *,'EOS: Unsupported vector data type in velocity.'
    !  call sys_halt()
    !end if
    !
    !if (PRESENT(rdefect)) then
    !  if ((rsolution%cdataType .NE. ST_DOUBLE) .OR. &
    !      (rdefect%cdataType .NE. ST_DOUBLE)) then
    !    PRINT *,'EOS: Unsupported vector data type in solution/defect'
    !    call sys_halt()
    !  end if
    !end if

    if (.not. rconfig%bconstViscosity) then
      call output_line ('Only constant viscosity supported at the moment', &
          OU_CLASS_ERROR,OU_MODE_STD,'conv_JumpStabilisation3d')
      call sys_halt()
    end if

    if (rconfig%dnu .eq. SYS_INFINITY_DP) then
      call output_line ('Viscosity parameter nu not initialised', &
          OU_CLASS_ERROR,OU_MODE_STD,'conv_JumpStabilisation3d')
      call sys_halt()
    end if

    if (rconfig%cjump .eq. CONV_JUMP_UNIFIEDEDGE) then
!      if (present(rdefect)) then
!
!        ! Modify the defect?
!        if (iand(cdef,CONV_MODDEFECT) .ne. 0) then
!          call jstab_matvecUEOJumpStabilBlk3d ( &
!              rconfig%dgamma,rconfig%dgammastar,rconfig%ccubType,rconfig%dnu,&
!              rmatrix,rsolution,rdefect,-rconfig%dtheta,1.0_DP,&
!              rdiscretisation)
!        end if
!
!      end if

      ! Modify the matrix?
      if (iand(cdef,CONV_MODMATRIX) .ne. 0) then
        call jstab_calcUEOJumpStabilisation (&
          rmatrix,rconfig%dgamma,rconfig%dgammastar,rconfig%deojEdgeExp,&
          rconfig%dtheta,rconfig%ccubType,rconfig%dnu,rdiscretisation,&
          rperfconfig=rperfconfig)
      end if

    else if (rconfig%cjump .eq. CONV_JUMP_REACTIVE) then
!      if (present(rdefect)) then
!
!        ! Modify the defect?
!        if (iand(cdef,CONV_MODDEFECT) .ne. 0) then
!          call jstab_matvecReacJumpStabilBlk3d ( &
!              rconfig%dgamma,rconfig%ccubType,rconfig%dnu,&
!              rmatrix,rsolution,rdefect,-rconfig%dtheta,1.0_DP,&
!              rdiscretisation)
!        end if
!
!      end if
!
      ! Modify the matrix?
      if (iand(cdef,CONV_MODMATRIX) .ne. 0) then
        call jstab_calcReacJumpStabilisation (&
          rmatrix,rconfig%dgamma,rconfig%dtheta,&
          rconfig%ccubType,rconfig%dnu,rdiscretisation,&
        rperfconfig=rperfconfig)
      end if

    else
      call output_line ('Unknown jump stabilisation', &
          OU_CLASS_ERROR,OU_MODE_STD,'conv_JumpStabilisation3d')
      call sys_halt()
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine conv_streamDiff2Blk2dMat (rconfig,rmatrix,rvelocity,&
      ffunctionCoefficient,rcollection,rcubatureInfo,rperfconfig)

!<description>
  ! Standard streamline diffusion method to set up the operator
  ! <tex>
  ! $$ dtheta  *  (                dalpha * MASS
  !                  +              dbeta * dnu * LAPLACE
  !                  +             ddelta * u * grad(.)
  !                  +            dnewton * (.) * grad(u)
  !                  +   ddeltaTransposed * grad(.)^T * u
  !                  +  dnewtonTransposed * grad(u)^T * (.) ) $$
  ! </tex>
  ! into a matrix rmatrix. 2D-version (X- and Y-velocity). The optional
  ! parameter u=rvelocity defines the evaluation point of the nonlinearity
  ! if there is one.
  !
  ! The configuration how the routine should react is to be configured
  ! in the configuration block rconfig.
  !
  ! The routine works with a block matrix rmatrix and allows to include
  ! extended operators like the Newton operator (Frechet-derivative of the
  ! convective part).
!</description>

!<input>
  ! Configuration block for the streamline diffusion scheme
  type(t_convStreamDiff2), intent(in) :: rconfig

  ! OPTIONAL: Velocity field where to evaluate the nonlinearity.
  ! Can be omitted if there is no nonlinearity to be assembled.
  type(t_vectorBlock), intent(in), target, optional :: rvelocity

  ! OPTIONAL: A callback routine for a nonconstant <tex>$ \nu $</tex>.
  ! Must be present if dnu is set to be nonconstant by rconfig%bconstNu=.false.
  ! or rconfig%bconstAlpha=.false..
  include '../DOFMaintenance/intf_sdcoefficient.inc'
  optional :: ffunctionCoefficient

  ! OPTIONAL: A scalar cubature information structure that specifies the cubature
  ! formula(s) to use. If not specified, default settings are used.
  type(t_scalarCubatureInfo), intent(in), optional, target :: rcubatureInfo

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! OPTIONAL: A collection structure. This structure is given to the
  ! callback function for nonconstant coefficients to provide additional
  ! information.
  type(t_collection), intent(inout), target, optional :: rcollection

  ! System block matrix.
  ! The nonlinear operator is added to the matrix.
  type(t_matrixBlock), intent(inout), target :: rmatrix
!</inputoutput>

!</subroutine>

    type(t_collection) :: rincorporateCollection
    logical :: bsimpleAij

    ! If there is only one block, specify the matrix as being 'simple'.
    ! Otherwise, the matrix is simple if the two diagonal blocks are the same.
    bsimpleAij = (rmatrix%nblocksPerCol .eq. 1) .and. (rmatrix%nblocksPerRow .eq. 1)
    if (.not. bsimpleAij) &
      bsimpleAij = lsyssc_isMatrixContentShared(&
            rmatrix%RmatrixBlock(1,1),rmatrix%RmatrixBlock(2,2))

    ! Calculate the operator with conv_streamDiff2Blk2dCalc.
    ! Use the callback routine conv_sdIncorpToMatrix2D to incorporate
    ! the data into the given block matrix.

    rincorporateCollection%p_rmatrixQuickAccess1 => rmatrix

    call conv_streamDiff2Blk2dCalc (rconfig,rmatrix%p_rblockDiscrTrial,&
        bsimpleAij,&
        rvelocity,&
        conv_sdIncorpToMatrix2D,rincorporateCollection,&
        ffunctionCoefficient,rcollection,rcubatureInfo,rperfconfig)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine conv_streamDiff2Blk2dDef (rconfig,rmatrix,rx,rd,rvelocity,&
      ffunctionCoefficient,rcollection,rcubatureInfo,rperfconfig)

!<description>
  ! Standard streamline diffusion method to set up the operator
  ! <tex>
  ! $$ dtheta  *  (                dalpha * MASS
  !                  +              dbeta * dnu * LAPLACE
  !                  +             ddelta * u * grad(.)
  !                  +            dnewton * (.) * grad(u)
  !                  +   ddeltaTransposed * grad(.)^T * u
  !                  +  dnewtonTransposed * grad(u)^T * (.) ) $$
  ! </tex>
  ! into a defect vector rd:
  !  rd = rd - operator(rvelocity)*rx. 2D-version (X- and Y-velocity). The optional
  ! parameter u=rvelocity defines the evaluation point of the nonlinearity
  ! if there is one.
  !
  ! The configuration how the routine should react is to be configured
  ! in the configuration block rconfig.
  !
  ! The routine works with a block matrix rmatrix and allows to include
  ! extended operators like the Newton operator (Frechet-derivative of the
  ! convective part).
!</description>

!<input>
  ! Configuration block for the streamline diffusion scheme
  type(t_convStreamDiff2), intent(in) :: rconfig

  ! Block matrix specifying the structure of the velocoity submatrices.
  type(t_matrixBlock), intent(in) :: rmatrix

  ! Solution vector for creating the defect.
  type(t_vectorBlock), intent(in), target :: rx

  ! OPTIONAL: Velocity field where to evaluate the nonlinearity.
  ! Can be omitted if there is no nonlinearity to be assembled.
  type(t_vectorBlock), intent(in), target, optional :: rvelocity

  ! OPTIONAL: A callback routine for a nonconstant <tex>$ \nu $</tex>.
  ! Must be present if dnu is set to be nonconstant by rconfig%bconstNu=.false.
  ! or rconfig%bconstAlpha=.false..
  include '../DOFMaintenance/intf_sdcoefficient.inc'
  optional :: ffunctionCoefficient

  ! OPTIONAL: A scalar cubature information structure that specifies the cubature
  ! formula(s) to use. If not specified, default settings are used.
  type(t_scalarCubatureInfo), intent(in), optional, target :: rcubatureInfo

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! OPTIONAL: A collection structure. This structure is given to the
  ! callback function for nonconstant coefficients to provide additional
  ! information.
  type(t_collection), intent(inout), target, optional :: rcollection

  ! Defect vector where to incorporate the defect
  type(t_vectorBlock), intent(inout), target :: rd
!</inputoutput>

!</subroutine>

    type(t_collection) :: rincorporateCollection
    logical :: bsimpleAij

    ! Calculate the operator with conv_streamDiff2Blk2dCalc.
    ! Use the callback routine conv_sdIncorpToMatrix2D to incorporate
    ! the data into the given block matrix.

    rincorporateCollection%p_rvectorQuickAccess1 => rx
    rincorporateCollection%p_rvectorQuickAccess2 => rd

    ! If there is only one block, specify the matrix as being 'simple'.
    ! Otherwise, the matrix is simple if the two diagonal blocks are the same.
    bsimpleAij = (rmatrix%nblocksPerCol .eq. 1) .and. (rmatrix%nblocksPerRow .eq. 1)
    if (.not. bsimpleAij) &
      bsimpleAij = lsyssc_isMatrixContentShared(&
            rmatrix%RmatrixBlock(1,1),rmatrix%RmatrixBlock(2,2))

    call conv_streamDiff2Blk2dCalc (rconfig,rmatrix%p_rblockDiscrTrial,&
        bsimpleAij,&
        rvelocity,&
        conv_sdIncorpToDefect2D,rincorporateCollection,&
        ffunctionCoefficient,rcollection,rcubatureInfo,rperfconfig)

    ! Note: When calculating the defect, it does not matter if
    ! one uses a call to lsyssc_isMatrixContentShared or .TRUE. in the above
    ! call, but the latter case probably needs slightly more time.

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine conv_streamDiff2Blk2dCalc (rconfig,rvelocityDiscr,bsimpleAij,rvelocity,&
      fincorporate,rincorporateCollection,ffunctionCoefficient,rcollection,&
      rcubatureInfo,rperfconfig)

!<description>
  ! Standard streamline diffusion method to set up the operator
  ! <tex>
  ! $$ dtheta  *  (                dalpha * MASS
  !                  +              dbeta * dnu * LAPLACE
  !                  +             ddelta * u * grad(.)
  !                  +            dnewton * (.) * grad(u)
  !                  +   ddeltaTransposed * grad(.)^T * u
  !                  +  dnewtonTransposed * grad(u)^T * (.) ) $$
  ! </tex>
  ! 2D-version (X- and Y-velocity). The optional
  ! parameter u=rvelocity defines the evaluation point of the nonlinearity
  ! if there is one.
  !
  ! The configuration how the routine should react is to be configured
  ! in the configuration block rconfig.
  !
  ! The routine works with a block matrix rmatrix and allows to include
  ! extended operators like the Newton operator (Frechet-derivative of the
  ! convective part).
!</description>

!<input>
  ! Configuration block for the streamline diffusion scheme
  type(t_convStreamDiff2), intent(in) :: rconfig

  ! Callback routine that incorporates local element matrices in
  ! a defect vector or a matrix.
  interface

    subroutine fincorporate (inonlinComplexity,nelements,indof,dtheta,Idofs,&
      DentryA11,DentryA12,DentryA21,DentryA22,rcollection,KentryA11,KentryA12)

      use fsystem
      use collection, only: t_collection

      ! Complexity of the matrix
      integer, intent(in) :: inonlinComplexity

      ! Number of elements and vertices on each element
      integer, intent(in) :: nelements,indof

      ! Weight for the local matrices
      real(DP), intent(in) :: dtheta

      ! The DOF`s on all elements the routine should work on
      integer, dimension(:,:), intent(in) :: Idofs

      ! Temporary arrays for positions of the local matrices in
      ! the global matrix for A11/A22 and A12/A21. Can be undefined.
      integer, dimension(:,:,:), intent(inout) :: KentryA11,KentryA12

      ! Values of the local matrices.
      real(DP), dimension(:,:,:), intent(in) :: DentryA11,DentryA12,DentryA21,DentryA22

      ! Collection structure. p_rmatrixQuickAccess1 points to the matrix
      ! where to incorporate the data.
      type(t_collection), intent(in) :: rcollection

    end subroutine

  end interface

  ! Discretisation stucture for the velocity components
  type(t_blockDiscretisation), intent(in) :: rvelocityDiscr

  ! Can be set to TRUE to specify that the entries in the velocity submatrices
  ! are saved 'compressed' -- which is the case if A11=A22 share its data.
  logical, intent(in) :: bsimpleAij

  ! Collection structure which is passed to fincorporate.
  type(t_collection), intent(in) :: rincorporatecollection

  ! OPTIONAL: Velocity field where to evaluate the nonlinearity.
  ! Can be omitted if there is no nonlinearity to be assembled.
  type(t_vectorBlock), intent(in), target, optional :: rvelocity

  ! OPTIONAL: A callback routine for a nonconstant <tex>$ \nu $</tex>.
  ! Must be present if dnu is set to be nonconstant by rconfig%bconstNu=.false.
  ! or rconfig%bconstAlpha=.false..
  include '../DOFMaintenance/intf_sdcoefficient.inc'
  optional :: ffunctionCoefficient

  ! OPTIONAL: A scalar cubature information structure that specifies the cubature
  ! formula(s) to use. If not specified, default settings are used.
  type(t_scalarCubatureInfo), intent(in), optional, target :: rcubatureInfo

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! OPTIONAL: A collection structure. This structure is given to the
  ! callback function for nonconstant coefficients to provide additional
  ! information.
  type(t_collection), intent(inout), target, optional :: rcollection
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iel,idofe,jdofe,jdfg,NVE
    real(DP) :: du1loc,du2loc,du1locx,du1locy,du2locx,du2locy,db,dbx,dby,OM
    real(DP) :: HBASI1,HBASI2,HBASI3,HBASJ1,HBASJ2,HBASJ3,HSUMI,HSUMJ
    real(DP) :: AH11,AH22,AH12,AH21

    ! Underlying triangulation
    type(t_triangulation), pointer :: p_rtriangulation
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    integer, dimension(:,:), pointer :: p_IedgesAtElement,p_IverticesAtElement

    ! An array receiving the coordinates of cubature points on
    ! the reference element for all elements in a set.
    real(DP), dimension(:,:), pointer :: p_DcubPtsRef

    ! Arrays for saving Jacobian determinants and matrices
    real(DP), dimension(:,:), pointer :: p_Ddetj

    ! For every cubature point on the reference element,
    ! the corresponding cubature weight
    real(DP), dimension(:), pointer :: p_Domega

    ! number of cubature points on the reference element
    integer :: ncubp,icubp

    ! Derivative qualifiers for evaluating the finite elements.
    logical, dimension(EL_MAXNDER) :: Bder

    ! Element evaluation tag; collects some information necessary for evaluating
    ! the elements.
    integer(I32) :: cevaluationTag

    ! An element evaluation set for evaluating elements.
    type(t_evalElementSet) :: revalElementSet
    type(t_domainIntSubset) :: rintSubset
    logical :: bcubPtsInitialised

    ! Number of elements in the current block
    integer :: nelementsPerBlock

    ! A pointer to an element-number list
    integer, dimension(:), pointer :: p_IelementList

    ! An allocateable array accepting the DOF`s of a set of elements.
    integer, dimension(:,:), pointer :: p_Idofs

    ! Allocateable arrays for the values of the basis functions -
    ! for test and trial spaces.
    real(DP), dimension(:,:,:,:), pointer :: p_Dbas

    ! Index arrays that define the positions of local matrices in the global matrix.
    integer, dimension(:,:,:), allocatable :: Kentry11, Kentry12

    ! Contributions for the submatrices A11, A12, A21, A22.
    real(DP), dimension(:,:,:), allocatable :: DentryA11
    real(DP), dimension(:,:,:), allocatable :: DentryA12
    real(DP), dimension(:,:,:), allocatable :: DentryA21
    real(DP), dimension(:,:,:), allocatable :: DentryA22

    ! An array with local DELTA`s, each DELTA for one element
    real(DP), dimension(:), allocatable :: DlocalDelta

    ! An array for the viscosity coefficients in all cubature points on
    ! all elements of the current element set
    real(DP), dimension(:,:), allocatable :: Dnu,Dalpha

    ! Local values in a cubature point
    real(DP) :: dnuloc,dalphaloc

    ! The discretisation - for easier access
    type(t_spatialDiscretisation), pointer :: p_rdiscr

    ! Variables specifying the current element set
    integer :: IELset,IELmax

    ! Number of DOF`s in the current element distribution
    integer :: indof

    ! Maximum norm of the vector field and its reciprocal
    real(DP) :: dumax,dumaxR

    ! Some specifies which type of nonlinearity we have.
    ! =0: linear problem
    ! =1: nonlinear problem with A11=A22 and A12=A21=0
    ! =2: nonlinear problem with A11 and A22 different, A12=A21=0
    ! =3: nonlinear problem with all Aij being present, but grad(u) not needed
    ! =4: nonlinear problem with all Aij being present, grad(u) needed
    integer :: inonlinComplexity

    ! Pointer to the velocity field in the cubature points.
    real(DP), dimension(:,:,:), allocatable :: Dvelocity

    ! Pointer to the velocity X- and Y-derivative in the cubature points
    real(DP), dimension(:,:,:), allocatable :: DvelocityUderiv
    real(DP), dimension(:,:,:), allocatable :: DvelocityVderiv

    ! Pointers to the FE vector data of the velocity field
    real(DP), dimension(:), pointer :: p_Du1,p_Du2

    ! A temporary block vector containing only the velocity
    type(t_vectorBlock) :: rvelocitytemp

    ! Current assembly block, cubature formula, element type,...
    integer :: ielementDistr,icubatureBlock,NEL
    integer(I32) :: cevalTag, celement, ccubature, ctrafoType
    type(t_scalarCubatureInfo), target :: rtempCubatureInfo
    type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo
    type(t_stdCubatureData) :: rcubatureData
    type(t_stdFEBasisEvalData) :: rfeBasisEvalData

    ! Pointer to the performance configuration
    type(t_perfconfig), pointer :: p_rperfconfig

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => conv_perfconfig
    end if

    ! Initialise the derivative flags. We need function values and
    ! 1st order derivatives.
    Bder = .false.
    Bder(DER_FUNC2D) = .true.
    Bder(DER_DERIV2D_X) = .true.
    Bder(DER_DERIV2D_Y) = .true.

    ! Shortcut to the spatial discretisation.
    ! We assume the same for all, A11, A12, A21 and A22.
    p_rdiscr => rvelocityDiscr%RspatialDiscr(1)

    ! Get some information about the triangulation
    p_rtriangulation => p_rdiscr%p_rtriangulation
    call storage_getbase_double2d (p_rtriangulation%h_DvertexCoords,&
                                   p_DvertexCoords)
    call storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
                                p_IverticesAtElement)
    call storage_getbase_int2d (p_rtriangulation%h_IedgesAtElement,&
                                p_IedgesAtElement)

    ! How complex is our nonlinearity.
    inonlinComplexity = 0

    if (rconfig%dupsam .ne. 0.0_DP) inonlinComplexity = 1
    if (rconfig%ddelta .ne. 0.0_DP) inonlinComplexity = 1

    if ((inonlinComplexity .eq. 1) .and. .not. bsimpleAij) &
      inonlinComplexity = 2

    ! inonlinComplexity = 3 does not happen up to now.

    if ((rconfig%ddeltaT .ne. 0.0_DP) .or. (rconfig%dnewtonT .ne. 0.0_DP) .or. &
        (rconfig%ddeltaTadj .ne. 0.0_DP) .or. &
        (rconfig%dnewton .ne. 0.0_DP) .or. (rconfig%dbetaT .ne. 0.0_DP)) &
      inonlinComplexity = 4

    ! Get the maximum velocity -- if we have a nonlinearity
    dumax = 0.0_DP
    dumaxR = 0.0_DP

    if (inonlinComplexity .gt. 0) then

      if (.not. present(rvelocity)) then
        call output_line ('No velocity field present!', &
            OU_CLASS_ERROR,OU_MODE_STD,'conv_streamDiff2Blk2dMat')
        call sys_halt()
      end if

      ! Create a temp vector with the velocity components to get the
      ! maximum velocity magnitude.
      call lsysbl_deriveSubvector(rvelocity,rvelocitytemp,1,NDIM2D,.true.)

      call fevl_getVectorMagnitude(rvelocitytemp,dumax)
      dumax = max(dumax,1.0E-8_DP)
      dumaxr = 1.0_DP/dumax

      call lsysbl_releaseVector (rvelocitytemp)

      ! Get pointers to the velocity field
      call lsyssc_getbase_double (rvelocity%RvectorBlock(1),p_Du1)
      call lsyssc_getbase_double (rvelocity%RvectorBlock(2),p_Du2)

    end if

    ! Do we have an assembly structure?
    ! If we do not have it, create a cubature info structure that
    ! defines how to do the assembly.
    if (.not. present(rcubatureInfo)) then
      call spdiscr_createDefCubStructure(p_rdiscr,rtempCubatureInfo,CUB_GEN_DEPR_EVAL)
      p_rcubatureInfo => rtempCubatureInfo
    else
      p_rcubatureInfo => rcubatureInfo
    end if

    ! Loop over the element blocks. Each defines a separate cubature formula.
    do icubatureBlock = 1,p_rcubatureInfo%ninfoBlockCount

      ! Get typical information: Number of elements, element list,...
      call spdiscr_getStdDiscrInfo (icubatureBlock,p_rcubatureInfo,p_rdiscr,&
          ielementDistr,celement,ccubature,NEL,p_IelementList)

      ! Cancel if this element list is empty.
      if (NEL .le. 0) cycle

      ! For saving some memory in smaller discretisations, we calculate
      ! the number of elements per block. For smaller triangulations,
      ! this is NEL. If there are too many elements, it is at most
      ! NELEMSIM. This is only used for allocating some arrays.
      nelementsPerBlock = min(p_rperfconfig%NELEMSIM, NEL)

      ! Initialise cubature and evaluation of the FE basis
      call easminfo_initStdCubature(ccubature,rcubatureData)

      ! For faster access...
      ncubp = cub_igetNumPts(ccubature)
      p_DcubPtsRef => rcubatureData%p_DcubPtsRef
      p_Domega => rcubatureData%p_Domega

      ! Initialise the evaluation structure for the FE basis
      call easminfo_initStdFEBasisEval(celement,&
          elem_getMaxDerivative(celement),rcubatureData%ncubp,nelementsPerBlock,rfeBasisEvalData)

      ! Get the basic element shape in the element distribution;
      ! 3=triangles, 4=quads
      NVE = rfeBasisEvalData%NVE

      ! For faster access...
      p_Dbas => rfeBasisEvalData%p_Dbas
      p_Idofs => rfeBasisEvalData%p_Idofs
      ctrafoType = rfeBasisEvalData%ctrafoType

      ! Get the number of local DOF`s for trial/test functions.
      ! We assume trial and test functions to be the same.
      indof = rfeBasisEvalData%ndofLocal

      ! Allocate memory for array with local DELTA`s
      allocate(DlocalDelta(nelementsPerBlock))
      call lalg_clearVectorDble (DlocalDelta,nelementsPerBlock)

      ! Allocate memory for the coefficients.
      ! If the coefficient is constant, fill it with the value.
      ! Otherwise, it is filled with the values later.
      allocate(Dnu(ncubp,nelementsPerBlock))
      if (rconfig%bconstNu) then
        Dnu(:,:) = rconfig%dnu
      end if

      ! The same for the ALPHA-coefficients
      allocate(Dalpha(ncubp,nelementsPerBlock))
      if (rconfig%bconstAlpha) then
        Dalpha(:,:) = 1.0_DP
      end if

      ! Allocate memory for the indices where to find the local matrices
      ! in the global matrix. These are just dummy arrays for the callback routine.
      allocate(Kentry11(indof,indof,nelementsPerBlock))
      allocate(Kentry12(indof,indof,nelementsPerBlock))

      ! Allocate memory for the local matrices
      allocate(DentryA11(indof,indof,nelementsPerBlock))
      allocate(DentryA12(indof,indof,nelementsPerBlock))
      allocate(DentryA21(indof,indof,nelementsPerBlock))
      allocate(DentryA22(indof,indof,nelementsPerBlock))

      ! Allocate memory for the velocity in the cubature points.
      allocate(Dvelocity(NDIM2D,ncubp,nelementsPerBlock))
      allocate(DvelocityUderiv(NDIM2D,ncubp,nelementsPerBlock))
      allocate(DvelocityVderiv(NDIM2D,ncubp,nelementsPerBlock))

      ! Initialisation of the element set.
      call elprep_init(revalElementSet)

      ! Indicate that cubature points must still be initialised in the element set.
      bcubPtsInitialised = .false.

      ! Loop over the elements - blockwise.
      do IELset = 1, size(p_IelementList), nelementsPerBlock

        ! Number of the last element in the current set
        IELmax = min(size(p_IelementList),IELset-1+nelementsPerBlock)

        ! Calculate the global DOF`s into IdofsTrial / IdofsTest.
        !
        ! More exactly, we call dof_locGlobMapping_mult to calculate all the
        ! global DOF`s of our NELEMSIM elements simultaneously.
        call dof_locGlobMapping_mult(p_rdiscr, p_IelementList(IELset:IELmax), &
                                     p_Idofs)

        ! Ok, we found the positions of the local matrix entries
        ! that we have to change.
        ! To calculate the matrix contributions, we have to evaluate
        ! the elements to give us the values of the basis functions
        ! in all the DOF`s in all the elements in our set.

        ! Get the element evaluation tag of all FE spaces. We need it to evaluate
        ! the elements later. All of them can be combined with OR, what will give
        ! a combined evaluation tag.
        cevaluationTag = elem_getEvaluationTag(celement)
        cevaluationTag = ior(cevaluationTag,elem_getEvaluationTag(EL_Q1))

        ! In the first loop, calculate the coordinates on the reference element.
        ! In all later loops, use the precalculated information.
        !
        ! Note: Why not using
        !   if (IELset .EQ. 1) then
        ! here, but this strange concept with the boolean variable?
        ! Because the if-command does not work with OpenMP! bcubPtsInitialised
        ! is a local variable and will therefore ensure that every thread
        ! is initialising its local set of cubature points!
        if (.not. bcubPtsInitialised) then
          bcubPtsInitialised = .true.
          cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REFPOINTS)
        else
          cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
        end if

        ! Calculate all information that is necessary to evaluate the finite element
        ! on all cells of our subset. This includes the coordinates of the points
        ! on the cells.
        call elprep_prepareSetForEvaluation (revalElementSet,&
            cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
            ctrafoType, p_DcubPtsRef(:,1:ncubp), rperfconfig=rperfconfig)
        p_Ddetj => revalElementSet%p_Ddetj

        ! Calculate the values of the basis functions.
        ! Pass p_DcubPts as point coordinates, which point either to the
        ! coordinates on the reference element (the same for all elements)
        ! or on the real element - depending on whether this is a
        ! parametric or nonparametric element.
        call elem_generic_sim2 (celement,revalElementSet, Bder, p_Dbas)

        ! Probably evaluate the velocity field in the cubature points
        if ((inonlinComplexity .ge. 1) .and. (inonlinComplexity .le. 3)) then

          if (present(rvelocity)) then

            ! We need only the u.
            ! Loop over all elements in the current set
            do IEL=1,IELmax-IELset+1

              ! Loop over all cubature points on the current element
              do ICUBP = 1, ncubp

                du1loc = 0.0_DP
                du2loc = 0.0_DP

                ! Perform a loop through the trial DOF`s.
                do JDOFE=1,indof

                  ! Get the value of the (test) basis function
                  ! phi_i (our "O") in the cubature point:

                  db = p_Dbas(JDOFE,1,ICUBP,IEL)

                  ! Sum up to the value in the cubature point

                  JDFG = p_Idofs(JDOFE,IEL)
                  du1loc = du1loc + p_Du1(JDFG)*db
                  du2loc = du2loc + p_Du2(JDFG)*db

                end do ! JDOFE

                ! Save the computed velocity
                Dvelocity(1,ICUBP,IEL) = du1loc
                Dvelocity(2,ICUBP,IEL) = du2loc

              end do ! ICUBP

            end do ! IEL

          end if

        else if (inonlinComplexity .ge. 4) then

          if (present(rvelocity)) then

            ! We need the velocity an the gradients.
            ! Loop over all elements in the current set
            do IEL=1,IELmax-IELset+1

              ! Loop over all cubature points on the current element
              do ICUBP = 1, ncubp

                du1loc = 0.0_DP
                du2loc = 0.0_DP
                du1locx = 0.0_DP
                du1locy = 0.0_DP
                du2locx = 0.0_DP
                du2locy = 0.0_DP

                ! Perform a loop through the trial DOF`s.
                do JDOFE=1,indof

                  ! Get the value of the (test) basis function
                  db = p_Dbas(JDOFE,1,ICUBP,IEL)
                  dbx = p_Dbas(JDOFE,DER_DERIV_X,ICUBP,IEL)
                  dby = p_Dbas(JDOFE,DER_DERIV_Y,ICUBP,IEL)

                  ! Sum up to the value in the cubature point
                  JDFG = p_Idofs(JDOFE,IEL)
                  du1loc = du1loc + p_Du1(JDFG)*db
                  du2loc = du2loc + p_Du2(JDFG)*db
                  du1locx = du1locx + p_Du1(JDFG)*dbx
                  du1locy = du1locy + p_Du1(JDFG)*dby
                  du2locx = du2locx + p_Du2(JDFG)*dbx
                  du2locy = du2locy + p_Du2(JDFG)*dby

                end do ! JDOFE

                ! Save the computed velocity
                Dvelocity(1,ICUBP,IEL) = du1loc
                Dvelocity(2,ICUBP,IEL) = du2loc
                DvelocityUderiv(1,ICUBP,IEL) = du1locx
                DvelocityUderiv(2,ICUBP,IEL) = du1locy
                DvelocityVderiv(1,ICUBP,IEL) = du2locx
                DvelocityVderiv(2,ICUBP,IEL) = du2locy

              end do ! ICUBP

            end do ! IEL

          end if

        end if

        ! In case we have nonconstant coefficients, calculate
        ! the viscosity coefficients in the cubature points
        ! using our callback routine.
        if (.not. (rconfig%bconstNu .and. rconfig%bconstAlpha)) then
          ! Prepare the call to the evaluation routine of the analytic function.
          call domint_initIntegrationByEvalSet (revalElementSet,rintSubset)
          !rintSubset%ielementDistribution = ielementDistr
          rintSubset%ielementStartIdx = IELset
          rintSubset%p_Ielements => p_IelementList(IELset:IELmax)
          rintSubset%p_IdofsTrial => p_Idofs
          rintSubset%celement = celement

          if (.not. rconfig%bconstNu) then
            call ffunctionCoefficient (0,p_rdiscr, &
                      int(IELmax-IELset+1),ncubp,&
                      revalElementSet%p_DpointsReal,&
                      p_Idofs,rintSubset,&
                      Dnu(:,1:IELmax-IELset+1_I32),rcollection)
          end if

          if (.not. rconfig%bconstAlpha) then
            call ffunctionCoefficient (1,p_rdiscr, &
                      int(IELmax-IELset+1),ncubp,&
                      revalElementSet%p_DpointsReal,&
                      p_Idofs,rintSubset,&
                      Dalpha(:,1:IELmax-IELset+1_I32),rcollection)
          end if

          ! Release the temporary domain integration structure again
          call domint_doneIntegration (rintSubset)
        end if

        ! Calculate local DELTA`s for streamline diffusion method.
        ! (cf. p. 121 in Turek`s CFD book).
        ! For every element, we need a local DELTA.
        ! Every local delta is weighted by the global "ddelta".
        ! If ddelta=0, we do not do anything as this disables the
        ! nonlinear term.
        ! If UPSAM=0.0, we have a central-difference like discretisation, which
        ! is one can see as the local stabilisation weight Delta is also = 0.0.
        ! In this case, we even switch of the calculation of the local Delta,
        ! as it is always =0.0, so we save a little bit time.
        if ((rconfig%dupsam .ne. 0.0_DP) .and. present(rvelocity)) then
          if (NVE .eq. 3) then
            call getLocalDeltaTriSim (rconfig%clocalh,&
                          Dvelocity,Dnu,p_Domega,p_Ddetj,duMaxR,&
                          rconfig%cstabiltype,rconfig%dupsam,&
                          p_IelementList(IELset:IELmax),p_rtriangulation,DlocalDelta)
          else
            call getLocalDeltaQuadSim (rconfig%clocalh,&
                          Dvelocity,Dnu,p_Domega,p_Ddetj,duMaxR,&
                          rconfig%cstabiltype,rconfig%dupsam,&
                          p_IelementList(IELset:IELmax),p_rtriangulation,DlocalDelta)
          end if
        end if

        ! Ok, we now use Dvelocity as coefficient array in the assembly
        ! of a bilinear form!
        !
        ! Clear the local matrices. If the Newton part is to be calculated,
        ! we must clear everything, otherwise only Dentry.
        select case (inonlinComplexity)
        case (0,1)
          ! There is only data for A11
          DentryA11 = 0.0_DP
        case (2)
          DentryA11 = 0.0_DP
          DentryA22 = 0.0_DP
        case default
          ! Clear all local Aij
          DentryA11 = 0.0_DP
          DentryA12 = 0.0_DP
          DentryA21 = 0.0_DP
          DentryA22 = 0.0_DP
        end select

        ! If ddelta != 0, set up the nonlinearity U*grad(u), probably with
        ! streamline diffusion stabilisation.
        if (rconfig%ddelta .ne. 0.0_DP) then

          if (inonlinComplexity .le. 1) then

            ! Build only A11, there is A11=A22.

            ! Loop over the elements in the current set.
            do IEL=1,IELmax-IELset+1

              ! Nonlinearity:
              !    ddelta * u_1 * grad(.)
              !  = ddelta * [ (DU1 dx + DU2 dy)                   ]
              !             [                   (DU1 dx + DU2 dy) ]
              !
              ! Loop over all cubature points on the current element
              do ICUBP = 1, ncubp

                ! Calculate the current weighting factor in the cubature formula
                ! in that cubature point.
                !
                ! Normally, we have to take the absolut value of the determinant
                ! of the mapping here!
                ! In 2D, the determinant is always positive, whereas in 3D,
                ! the determinant might be negative -- that is normal!
                ! But because this routine only works in 2D, we can skip
                ! the ABS here!

                OM = p_Domega(ICUBP)*p_Ddetj(ICUBP,IEL)

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
                ! Outer loop over the DOF`s i=1..indof on our current element,
                ! which corresponds to the basis functions Phi_i:

                do IDOFE=1,indof

                  ! Fetch the contributions of the (test) basis functions Phi_i
                  ! (our "O")  for function value and first derivatives for the
                  ! current DOF into HBASIy:

                  HBASI1 = p_Dbas(IDOFE,1,ICUBP,IEL)
                  HBASI2 = p_Dbas(IDOFE,2,ICUBP,IEL)
                  HBASI3 = p_Dbas(IDOFE,3,ICUBP,IEL)

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

                  ! Inner loop over the DOF`s j=1..indof, which corresponds to
                  ! the basis function Phi_j:

                  do JDOFE=1,indof

                    ! Fetch the contributions of the (trial) basis function Phi_j
                    ! (out "X") for function value and first derivatives for the
                    ! current DOF into HBASJy:

                    HBASJ1 = p_Dbas(JDOFE,1,ICUBP,IEL)
                    HBASJ2 = p_Dbas(JDOFE,2,ICUBP,IEL)
                    HBASJ3 = p_Dbas(JDOFE,3,ICUBP,IEL)

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
                    ! we do not have to worry about that.

                    HSUMJ = HBASJ2*du1loc+HBASJ3*du2loc

                    ! Finally calculate the contribution to the system
                    ! matrix. Depending on the configuration of ddelta,...
                    ! this is:
                    !
                    ! AH = n~_h(u_h,phi_j,phi_i)        | nonlinear part
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

                    AH11 = rconfig%ddelta * HSUMJ*(DlocalDelta(IEL)*HSUMI+HBASI1)

                    ! Weighten the calculated value AH by the cubature
                    ! weight OM and add it to the local matrix. After the
                    ! loop over all DOF`s is finished, each entry contains
                    ! the calculated integral.

                    DentryA11(JDOFE,IDOFE,IEL) = DentryA11(JDOFE,IDOFE,IEL)+OM*AH11

                  end do ! IDOFE

                end do ! JDOFE

              end do ! ICUBP

            end do ! IEL

          else

            ! Same loop, but affect both, A11 and A22.

            ! Loop over the elements in the current set.
            do IEL=1,IELmax-IELset+1

              do ICUBP = 1, ncubp

                OM = p_Domega(ICUBP)*p_Ddetj(ICUBP,IEL)

                du1loc = Dvelocity (1,ICUBP,IEL)
                du2loc = Dvelocity (2,ICUBP,IEL)

                do IDOFE=1,indof

                  HBASI1 = p_Dbas(IDOFE,1,ICUBP,IEL)
                  HBASI2 = p_Dbas(IDOFE,2,ICUBP,IEL)
                  HBASI3 = p_Dbas(IDOFE,3,ICUBP,IEL)

                  HSUMI = HBASI2*du1loc + HBASI3*du2loc

                  do JDOFE=1,indof

                    HBASJ1 = p_Dbas(JDOFE,1,ICUBP,IEL)
                    HBASJ2 = p_Dbas(JDOFE,2,ICUBP,IEL)
                    HBASJ3 = p_Dbas(JDOFE,3,ICUBP,IEL)

                    HSUMJ = HBASJ2*du1loc+HBASJ3*du2loc

                    AH11 = rconfig%ddelta * HSUMJ*(DlocalDelta(IEL)*HSUMI+HBASI1)

                    DentryA11(JDOFE,IDOFE,IEL) = DentryA11(JDOFE,IDOFE,IEL)+OM*AH11
                    DentryA22(JDOFE,IDOFE,IEL) = DentryA22(JDOFE,IDOFE,IEL)+OM*AH11

                  end do ! IDOFE

                end do ! JDOFE

              end do ! ICUBP

            end do ! IEL

          end if

        end if


        ! If ddeltaAdj != 0, set up the adjoint nonlinearity (.,U*grad(u)),
        ! i.e., the convection is applied to the test space.
        ! Basically, this exchanges trial and test space.

        if (rconfig%ddeltaAdj .ne. 0.0_DP) then

          if (inonlinComplexity .le. 1) then

            ! Build only A11, there is A11=A22.

            ! Loop over the elements in the current set.
            do IEL=1,IELmax-IELset+1

              ! The adjoint convection is defined by
              !
              !    ddeltaAdj * ( (trial), grad(u_1) * (test) )
              !
              ! which is obtained by a partial integration of the (negative) convection
              ! exploiting the solenoidality of the vector field,
              !
              !     ( - grad(u_1) * (trial), (test) )
              !         = ( (trial), grad(u_1) * (test) )
              !
              ! Loop over all cubature points on the current element
              do ICUBP = 1, ncubp

                ! Calculate the current weighting factor in the cubature formula
                ! in that cubature point.
                !
                ! Normally, we have to take the absolut value of the determinant
                ! of the mapping here!
                ! In 2D, the determinant is always positive, whereas in 3D,
                ! the determinant might be negative -- that is normal!
                ! But because this routine only works in 2D, we can skip
                ! the ABS here!

                OM = p_Domega(ICUBP)*p_Ddetj(ICUBP,IEL)

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
                ! Outer loop over the DOF`s i=1..indof on our current element,
                ! which corresponds to the basis functions Phi_i:

                do IDOFE=1,indof

                  ! Fetch the contributions of the (test) basis functions Phi_i
                  ! (our "O")  for function value and first derivatives for the
                  ! current DOF into HBASIy:

                  HBASI1 = p_Dbas(IDOFE,1,ICUBP,IEL)
                  HBASI2 = p_Dbas(IDOFE,2,ICUBP,IEL)
                  HBASI3 = p_Dbas(IDOFE,3,ICUBP,IEL)

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

                  ! Inner loop over the DOF`s j=1..indof, which corresponds to
                  ! the basis function Phi_j:

                  do JDOFE=1,indof

                    ! Fetch the contributions of the (trial) basis function Phi_j
                    ! (out "X") for function value and first derivatives for the
                    ! current DOF into HBASJy:

                    HBASJ1 = p_Dbas(JDOFE,1,ICUBP,IEL)
                    HBASJ2 = p_Dbas(JDOFE,2,ICUBP,IEL)
                    HBASJ3 = p_Dbas(JDOFE,3,ICUBP,IEL)

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
                    ! we do not have to worry about that.

                    HSUMJ = HBASJ2*du1loc+HBASJ3*du2loc

                    ! Finally calculate the contribution to the system
                    ! matrix. Depending on the configuration of ddelta,...
                    ! this is:
                    !
                    ! AH = n~_h(u_h,phi_j,phi_i)        | nonlinear part
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
                    ! The adjoint operator exchanges trial and test space,
                    ! so one obtains the formula
                    !
                    !     HSUMI * (Delta * HSUMJ + HBASJ1)
                    !
                    ! plus the terms for the Stokes and Mass matrix,
                    ! if their coefficient is <> 0.

                    AH11 = rconfig%ddeltaAdj * HSUMI*(DlocalDelta(IEL)*HSUMJ+HBASJ1)

                    ! Weighten the calculated value AH by the cubature
                    ! weight OM and add it to the local matrix. After the
                    ! loop over all DOF`s is finished, each entry contains
                    ! the calculated integral.

                    DentryA11(JDOFE,IDOFE,IEL) = DentryA11(JDOFE,IDOFE,IEL)+OM*AH11

                  end do ! IDOFE

                end do ! JDOFE

              end do ! ICUBP

            end do ! IEL

          else

            ! Same loop, but affect both, A11 and A22.

            ! Loop over the elements in the current set.
            do IEL=1,IELmax-IELset+1

              do ICUBP = 1, ncubp

                OM = p_Domega(ICUBP)*p_Ddetj(ICUBP,IEL)

                du1loc = Dvelocity (1,ICUBP,IEL)
                du2loc = Dvelocity (2,ICUBP,IEL)

                do IDOFE=1,indof

                  HBASI1 = p_Dbas(IDOFE,1,ICUBP,IEL)
                  HBASI2 = p_Dbas(IDOFE,2,ICUBP,IEL)
                  HBASI3 = p_Dbas(IDOFE,3,ICUBP,IEL)

                  HSUMI = HBASI2*du1loc + HBASI3*du2loc

                  do JDOFE=1,indof

                    HBASJ1 = p_Dbas(JDOFE,1,ICUBP,IEL)
                    HBASJ2 = p_Dbas(JDOFE,2,ICUBP,IEL)
                    HBASJ3 = p_Dbas(JDOFE,3,ICUBP,IEL)

                    HSUMJ = HBASJ2*du1loc+HBASJ3*du2loc

                    AH11 = rconfig%ddeltaAdj * HSUMI*(DlocalDelta(IEL)*HSUMJ+HBASJ1)

                    DentryA11(JDOFE,IDOFE,IEL) = DentryA11(JDOFE,IDOFE,IEL)+OM*AH11
                    DentryA22(JDOFE,IDOFE,IEL) = DentryA22(JDOFE,IDOFE,IEL)+OM*AH11

                  end do ! IDOFE

                end do ! JDOFE

              end do ! ICUBP

            end do ! IEL

          end if

        end if


        ! If dbeta != 0 or dalpha != 0, add the Laplace/Mass matrix to the
        ! local matrices.
        if ((.not. rconfig%bconstAlpha) .or. (rconfig%dalpha .ne. 0.0_DP) .or. &
            (rconfig%dbeta .ne. 0.0_DP) ) then

          if (inonlinComplexity .le. 1) then

            ! Build only A11, there is A11=A22.

            ! Loop over the elements in the current set.
            do IEL=1,IELmax-IELset+1

              ! Loop over all cubature points on the current element
              do ICUBP = 1, ncubp

                ! Current local nu / alpha
                dnuloc = rconfig%dbeta * Dnu(icubp,iel)
                dalphaloc = rconfig%dalpha * Dalpha(icubp,iel)

                ! Calculate the current weighting factor in the cubature formula
                ! in that cubature point.
                !
                ! Normally, we have to take the absolut value of the determinant
                ! of the mapping here!
                ! In 2D, the determinant is always positive, whereas in 3D,
                ! the determinant might be negative -- that is normal!
                ! But because this routine only works in 2D, we can skip
                ! the ABS here!

                OM = p_Domega(ICUBP)*p_Ddetj(ICUBP,IEL)

                ! Outer loop over the DOF`s i=1..indof on our current element,
                ! which corresponds to the basis functions Phi_i:

                do IDOFE=1,indof

                  ! Fetch the contributions of the (test) basis functions Phi_i
                  ! (our "O")  for function value and first derivatives for the
                  ! current DOF into HBASIy:

                  HBASI1 = p_Dbas(IDOFE,1,ICUBP,IEL)
                  HBASI2 = p_Dbas(IDOFE,2,ICUBP,IEL)
                  HBASI3 = p_Dbas(IDOFE,3,ICUBP,IEL)

                  ! Inner loop over the DOF`s j=1..indof, which corresponds to
                  ! the basis function Phi_j:

                  do JDOFE=1,indof

                    ! Fetch the contributions of the (trial) basis function Phi_j
                    ! (out "X") for function value and first derivatives for the
                    ! current DOF into HBASJy:

                    HBASJ1 = p_Dbas(JDOFE,1,ICUBP,IEL)
                    HBASJ2 = p_Dbas(JDOFE,2,ICUBP,IEL)
                    HBASJ3 = p_Dbas(JDOFE,3,ICUBP,IEL)

                    ! Finally calculate the contribution to the system
                    ! matrix. Depending on the configuration of DNU,
                    ! dalpha,... this decomposes into:
                    !
                    ! AH = dny*(grad(phi_j,grad(phi_i)) | -dny*Laplace(u) = -dbeta*Stokes
                    !    + dalpha*(phi_j*phi_i)         | Mass matrix

                    AH11 = dnuloc*(HBASI2*HBASJ2+HBASI3*HBASJ3) &
                          + dalphaloc*HBASI1*HBASJ1

                    ! Weighten the calculated value AH by the cubature
                    ! weight OM and add it to the local matrix. After the
                    ! loop over all DOF`s is finished, each entry contains
                    ! the calculated integral.

                    DentryA11(JDOFE,IDOFE,IEL) = DentryA11(JDOFE,IDOFE,IEL)+OM*AH11

                  end do ! IDOFE

                end do ! JDOFE

              end do ! ICUBP

            end do ! IEL

          else

            ! Same loop, but affect both, A11 and A22.

            ! Loop over the elements in the current set.
            do IEL=1,IELmax-IELset+1

              ! Loop over all cubature points on the current element
              do ICUBP = 1, ncubp

                dnuloc = rconfig%dbeta * Dnu(icubp,iel)
                dalphaloc = rconfig%dalpha * Dalpha(icubp,iel)

                OM = p_Domega(ICUBP)*p_Ddetj(ICUBP,IEL)

                do IDOFE=1,indof

                  HBASI1 = p_Dbas(IDOFE,1,ICUBP,IEL)
                  HBASI2 = p_Dbas(IDOFE,2,ICUBP,IEL)
                  HBASI3 = p_Dbas(IDOFE,3,ICUBP,IEL)

                  do JDOFE=1,indof

                    HBASJ1 = p_Dbas(JDOFE,1,ICUBP,IEL)
                    HBASJ2 = p_Dbas(JDOFE,2,ICUBP,IEL)
                    HBASJ3 = p_Dbas(JDOFE,3,ICUBP,IEL)

                    AH11 = dnuloc*(HBASI2*HBASJ2+HBASI3*HBASJ3) &
                          + dalphaloc*HBASI1*HBASJ1

                    DentryA11(JDOFE,IDOFE,IEL) = DentryA11(JDOFE,IDOFE,IEL)+OM*AH11
                    DentryA22(JDOFE,IDOFE,IEL) = DentryA22(JDOFE,IDOFE,IEL)+OM*AH11

                  end do ! IDOFE

                end do ! JDOFE

              end do ! ICUBP

            end do ! IEL

          end if

        end if

        ! Should we assemble the 'transposed Stokes' operator
        if (rconfig%dbetaT .ne. 0.0_DP) then

          ! Build only A11, there is A11=A22.

          ! Loop over the elements in the current set.
          do IEL=1,IELmax-IELset+1

            ! Loop over all cubature points on the current element
            do ICUBP = 1, ncubp

              ! Current local nu
              dnuloc = rconfig%dbetaT * Dnu(icubp,iel)

              ! Calculate the current weighting factor in the cubature formula
              ! in that cubature point.
              !
              ! Normally, we have to take the absolut value of the determinant
              ! of the mapping here!
              ! In 2D, the determinant is always positive, whereas in 3D,
              ! the determinant might be negative -- that is normal!
              ! But because this routine only works in 2D, we can skip
              ! the ABS here!

              OM = p_Domega(ICUBP)*p_Ddetj(ICUBP,IEL)

              ! Outer loop over the DOF`s i=1..indof on our current element,
              ! which corresponds to the basis functions Phi_i:

              do IDOFE=1,indof

                ! Fetch the contributions of the (test) basis functions Phi_i
                ! (our "O")  for function value and first derivatives for the
                ! current DOF into HBASIy:

                HBASI2 = p_Dbas(IDOFE,2,ICUBP,IEL)
                HBASI3 = p_Dbas(IDOFE,3,ICUBP,IEL)

                ! Inner loop over the DOF`s j=1..indof, which corresponds to
                ! the basis function Phi_j:

                do JDOFE=1,indof

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! (out "X") for function value and first derivatives for the
                  ! current DOF into HBASJy:

                  HBASJ2 = p_Dbas(JDOFE,2,ICUBP,IEL)
                  HBASJ3 = p_Dbas(JDOFE,3,ICUBP,IEL)

                  ! Finally calculate the contribution to the system
                  ! matrix. Depending on the configuration of DNU,
                  ! dalpha,... this decomposes into:
                  !
                  ! AH = dny*(grad(phi_j,grad(phi_i)) | -dny*Laplace(u) = -dbeta*Stokes
                  !    + dalpha*(phi_j*phi_i)         | Mass matrix

                  AH11 = dnuloc*(HBASI2*HBASJ2)
                  AH12 = dnuloc*(HBASI3*HBASJ2)
                  AH21 = dnuloc*(HBASI2*HBASJ3)
                  ! AH22 = AH11, so do not compute.

                  ! Weighten the calculated value AH by the cubature
                  ! weight OM and add it to the local matrix. After the
                  ! loop over all DOF`s is finished, each entry contains
                  ! the calculated integral.

                  DentryA11(JDOFE,IDOFE,IEL) = DentryA11(JDOFE,IDOFE,IEL)+OM*AH11
                  DentryA12(JDOFE,IDOFE,IEL) = DentryA12(JDOFE,IDOFE,IEL)+OM*AH12
                  DentryA21(JDOFE,IDOFE,IEL) = DentryA21(JDOFE,IDOFE,IEL)+OM*AH21
                  DentryA22(JDOFE,IDOFE,IEL) = DentryA22(JDOFE,IDOFE,IEL)+OM*AH11

                end do ! IDOFE

              end do ! JDOFE

            end do ! ICUBP

          end do ! IEL

        end if

        ! Should we assemble the Newton matrices?
        if (rconfig%dnewton .ne. 0.0_DP) then

          ! Newton operator
          !
          !    dnewton * ( grad(u_1) * (trial), (test) )
          !  = dnewton * [ (dx DU1) (dy DU1) ]
          !              [ (dx DU2) (dy DU2) ]
          !
          ! Loop over the elements in the current set.
          do IEL=1,IELmax-IELset+1

            ! Loop over all cubature points on the current element
            do ICUBP = 1, ncubp

              ! Calculate the current weighting factor in the cubature formula
              ! in that cubature point.
              !
              ! Normally, we have to take the absolut value of the determinant
              ! of the mapping here!
              ! In 2D, the determinant is always positive, whereas in 3D,
              ! the determinant might be negative -- that is normal!
              ! But because this routine only works in 2D, we can skip
              ! the ABS here!

              OM = p_Domega(ICUBP)*p_Ddetj(ICUBP,IEL)

              ! Current velocity in this cubature point:
              du1locx = DvelocityUderiv (1,ICUBP,IEL)
              du1locy = DvelocityUderiv (2,ICUBP,IEL)
              du2locx = DvelocityVderiv (1,ICUBP,IEL)
              du2locy = DvelocityVderiv (2,ICUBP,IEL)

              ! Outer loop over the DOF`s i=1..indof on our current element,
              ! which corresponds to the basis functions Phi_i:

              do IDOFE=1,indof

                ! Fetch the contributions of the (test) basis functions Phi_i
                ! (our "O")  for function value and first derivatives for the
                ! current DOF into HBASIy:

                HBASI1 = p_Dbas(IDOFE,1,ICUBP,IEL)

                ! Inner loop over the DOF`s j=1..indof, which corresponds to
                ! the basis function Phi_j:

                do JDOFE=1,indof

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! (out "X") for function value and first derivatives for the
                  ! current DOF into HBASJy:

                  HBASJ1 = p_Dbas(JDOFE,1,ICUBP,IEL)

                  ! Finally calculate the contribution to the system
                  ! matrices A11, A12, A21 and A22.
                  !
                  ! The Newton part is calculated as follows:
                  !
                  ! U * grad(V)  =  ( U * grad(.) ) V
                  !
                  !              =  ( U * grad(V1) )
                  !                 ( U * grad(V2) )
                  !
                  !              =  ( (U1) * (V1x) )
                  !                 ( (U2)   (V1y) )
                  !                 (              )
                  !                 ( (U1) * (V2x) )
                  !                 ( (U2)   (V2y) )
                  !
                  !              =  ( U1 * V1x  + U2 * V1y )
                  !                 ( U1 * V2x  + U2 * V2y )
                  !
                  !              =  ( V1_x  V1_y ) ( U1 )
                  !                 ( V2_x  V2_y ) ( U2 )
                  !
                  !              -> ( A11   A12  )
                  !                 ( A21   A22  )
                  !
                  ! With the velocity V=(u,v), we have to assemble:
                  ! grad(V)*U, which is realised in each cubature point as:
                  !   du/dx * phi_j*phi_i -> A11
                  !   du/dy * phi_j*phi_i -> A12
                  !   dv/dx * phi_j*phi_i -> A21
                  !   dv/dy * phi_j*phi_i -> A22

                  AH11 = rconfig%dnewton * du1locx * HBASJ1*HBASI1
                  AH12 = rconfig%dnewton * du1locy * HBASJ1*HBASI1
                  AH21 = rconfig%dnewton * du2locx * HBASJ1*HBASI1
                  AH22 = rconfig%dnewton * du2locy * HBASJ1*HBASI1

                  ! Weighten the calculated value AHxy by the cubature
                  ! weight OM and add it to the local matrices. After the
                  ! loop over all DOF`s is finished, each entry contains
                  ! the calculated integral.

                  DentryA11(JDOFE,IDOFE,IEL) = DentryA11(JDOFE,IDOFE,IEL)+OM*AH11
                  DentryA12(JDOFE,IDOFE,IEL) = DentryA12(JDOFE,IDOFE,IEL)+OM*AH12
                  DentryA21(JDOFE,IDOFE,IEL) = DentryA21(JDOFE,IDOFE,IEL)+OM*AH21
                  DentryA22(JDOFE,IDOFE,IEL) = DentryA22(JDOFE,IDOFE,IEL)+OM*AH22

                end do ! IDOFE

              end do ! JDOFE

            end do ! ICUBP

          end do ! IEL

        end if

        ! Should we assemble the ajoint Newton matrices?
        if (rconfig%dnewtonAdj .ne. 0.0_DP) then

          ! The Newton operator was
          !
          !    dnewton * ( grad(u_1) * (trial), (test) )
          !
          ! The adjoint Newton operator is a bit more complicated.
          ! It comes from the following calculation which invokes
          ! a partial integration, exploiting solenoidality,
          !
          !    ( - grad(u_1) * (trial), (test) )
          !        =    ( - ((trial)*grad) (u_1), (test) )
          !        =    ( u_1 , ((trial)*grad) (test) )
          !
          ! i.e., test and trial functions are both on the right-hand
          ! side of the scalar product. The "adjoint newton" operator
          ! is defined as
          !
          !    dnewtonAdj * ( u_1 , ((trial)*grad) (test) )
          !
          !
          ! Loop over the elements in the current set.
          do IEL=1,IELmax-IELset+1

            ! Loop over all cubature points on the current element
            do ICUBP = 1, ncubp

              ! Calculate the current weighting factor in the cubature formula
              ! in that cubature point.
              !
              ! Normally, we have to take the absolut value of the determinant
              ! of the mapping here!
              ! In 2D, the determinant is always positive, whereas in 3D,
              ! the determinant might be negative -- that is normal!
              ! But because this routine only works in 2D, we can skip
              ! the ABS here!

              OM = p_Domega(ICUBP)*p_Ddetj(ICUBP,IEL)

              ! Current velocity in this cubature point:
              du1loc = Dvelocity (1,ICUBP,IEL)
              du2loc = Dvelocity (2,ICUBP,IEL)

              ! Outer loop over the DOF`s i=1..indof on our current element,
              ! which corresponds to the basis functions Phi_i:

              do IDOFE=1,indof

                ! Fetch the contributions of the (test) basis functions Phi_i
                ! (our "O")  for function value and first derivatives for the
                ! current DOF into HBASIy:

                HBASI1 = p_Dbas(IDOFE,1,ICUBP,IEL)
                HBASI2 = p_Dbas(IDOFE,2,ICUBP,IEL)
                HBASI3 = p_Dbas(IDOFE,3,ICUBP,IEL)

                ! Inner loop over the DOF`s j=1..indof, which corresponds to
                ! the basis function Phi_j:

                do JDOFE=1,indof

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! (out "X") for function value and first derivatives for the
                  ! current DOF into HBASJy:

                  HBASJ1 = p_Dbas(JDOFE,1,ICUBP,IEL)
                  HBASJ2 = p_Dbas(JDOFE,2,ICUBP,IEL)
                  HBASJ3 = p_Dbas(JDOFE,3,ICUBP,IEL)

                  ! Finally calculate the contribution to the system
                  ! matrices A11, A12, A21 and A22.

                  AH11 = rconfig%dnewtonAdj * du1loc * HBASJ1*HBASI2
                  AH12 = rconfig%dnewtonAdj * du1loc * HBASJ1*HBASI3
                  AH21 = rconfig%dnewtonAdj * du2loc * HBASJ1*HBASI2
                  AH22 = rconfig%dnewtonAdj * du2loc * HBASJ1*HBASI3

                  ! Weighten the calculated value AHxy by the cubature
                  ! weight OM and add it to the local matrices. After the
                  ! loop over all DOF`s is finished, each entry contains
                  ! the calculated integral.

                  DentryA11(JDOFE,IDOFE,IEL) = DentryA11(JDOFE,IDOFE,IEL)+OM*AH11
                  DentryA12(JDOFE,IDOFE,IEL) = DentryA12(JDOFE,IDOFE,IEL)+OM*AH12
                  DentryA21(JDOFE,IDOFE,IEL) = DentryA21(JDOFE,IDOFE,IEL)+OM*AH21
                  DentryA22(JDOFE,IDOFE,IEL) = DentryA22(JDOFE,IDOFE,IEL)+OM*AH22

                end do ! IDOFE

              end do ! JDOFE

            end do ! ICUBP

          end do ! IEL

        end if

        ! Transposed convection operator
        !
        !  ddeltaTransposed * ( grad(trial)^T * u_1 , (test) )
        !  = ddeltaTransposed * [ (DU1 dx) (DU2 dx) ]
        !                       [ (DU1 dy) (DU2 dy) ]
        !
        ! Should we assemble the transposed convection matrices?
        if (rconfig%ddeltaT .ne. 0.0_DP) then

          ! Loop over the elements in the current set.
          do IEL=1,IELmax-IELset+1

            ! Loop over all cubature points on the current element
            do ICUBP = 1, ncubp

              ! Calculate the current weighting factor in the cubature formula
              ! in that cubature point.
              !
              ! Normally, we have to take the absolut value of the determinant
              ! of the mapping here!
              ! In 2D, the determinant is always positive, whereas in 3D,
              ! the determinant might be negative -- that is normal!
              ! But because this routine only works in 2D, we can skip
              ! the ABS here!

              OM = p_Domega(ICUBP)*p_Ddetj(ICUBP,IEL)

              ! Current velocity in this cubature point:
              du1loc = Dvelocity (1,ICUBP,IEL)
              du2loc = Dvelocity (2,ICUBP,IEL)

              ! Outer loop over the DOF`s i=1..indof on our current element,
              ! which corresponds to the basis functions Phi_i:

              do IDOFE=1,indof

                ! Fetch the contributions of the (test) basis functions Phi_i
                ! (our "O")  for function value and first derivatives for the
                ! current DOF into HBASIy:

                HBASI1 = p_Dbas(IDOFE,1,ICUBP,IEL)
                HBASI2 = p_Dbas(IDOFE,2,ICUBP,IEL)
                HBASI3 = p_Dbas(IDOFE,3,ICUBP,IEL)

                ! Inner loop over the DOF`s j=1..indof, which corresponds to
                ! the basis function Phi_j:

                do JDOFE=1,indof

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! (out "X") for function value and first derivatives for the
                  ! current DOF into HBASJy:

                  HBASJ1 = p_Dbas(JDOFE,1,ICUBP,IEL)
                  HBASJ2 = p_Dbas(JDOFE,2,ICUBP,IEL)
                  HBASJ3 = p_Dbas(JDOFE,3,ICUBP,IEL)

                  ! Finally calculate the contribution to the system
                  ! matrices A11, A12, A21 and A22.

                  AH11 = rconfig%ddeltaT * du1loc * HBASJ2*HBASI1
                  AH12 = rconfig%ddeltaT * du2loc * HBASJ2*HBASI1
                  AH21 = rconfig%ddeltaT * du1loc * HBASJ3*HBASI1
                  AH22 = rconfig%ddeltaT * du2loc * HBASJ3*HBASI1

                  ! Weighten the calculated value AHxy by the cubature
                  ! weight OM and add it to the local matrices. After the
                  ! loop over all DOF`s is finished, each entry contains
                  ! the calculated integral.

                  DentryA11(JDOFE,IDOFE,IEL) = DentryA11(JDOFE,IDOFE,IEL)+OM*AH11
                  DentryA12(JDOFE,IDOFE,IEL) = DentryA12(JDOFE,IDOFE,IEL)+OM*AH12
                  DentryA21(JDOFE,IDOFE,IEL) = DentryA21(JDOFE,IDOFE,IEL)+OM*AH21
                  DentryA22(JDOFE,IDOFE,IEL) = DentryA22(JDOFE,IDOFE,IEL)+OM*AH22

                end do ! IDOFE

              end do ! JDOFE

            end do ! ICUBP

          end do ! IEL

        end if

        ! Adjoint transposed convection operator.
        !
        ! The standard convection operator is defined by
        !
        !  ( grad(trial)^T * u_1 , (test) )
        !
        ! By exchanging trial and test functions, one obtains the
        ! adjoint transposed convection operator,
        !
        !  ddeltaTransposedAdj * ( (trial) , grad(test)^T * u_1 )
        !
        if (rconfig%ddeltaTadj .ne. 0.0_DP) then

          ! Loop over the elements in the current set.
          do IEL=1,IELmax-IELset+1

            ! Loop over all cubature points on the current element
            do ICUBP = 1, ncubp

              ! Calculate the current weighting factor in the cubature formula
              ! in that cubature point.
              !
              ! Normally, we have to take the absolut value of the determinant
              ! of the mapping here!
              ! In 2D, the determinant is always positive, whereas in 3D,
              ! the determinant might be negative -- that is normal!
              ! But because this routine only works in 2D, we can skip
              ! the ABS here!

              OM = p_Domega(ICUBP)*p_Ddetj(ICUBP,IEL)

              ! Current velocity in this cubature point:
              du1loc = Dvelocity (1,ICUBP,IEL)
              du2loc = Dvelocity (2,ICUBP,IEL)

              ! Outer loop over the DOF`s i=1..indof on our current element,
              ! which corresponds to the basis functions Phi_i:

              do IDOFE=1,indof

                ! Fetch the contributions of the (test) basis functions Phi_i
                ! (our "O")  for function value and first derivatives for the
                ! current DOF into HBASIy:

                HBASI1 = p_Dbas(IDOFE,1,ICUBP,IEL)
                HBASI2 = p_Dbas(IDOFE,2,ICUBP,IEL)
                HBASI3 = p_Dbas(IDOFE,3,ICUBP,IEL)

                ! Inner loop over the DOF`s j=1..indof, which corresponds to
                ! the basis function Phi_j:

                do JDOFE=1,indof

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! (out "X") for function value and first derivatives for the
                  ! current DOF into HBASJy:

                  HBASJ1 = p_Dbas(JDOFE,1,ICUBP,IEL)
                  HBASJ2 = p_Dbas(JDOFE,2,ICUBP,IEL)
                  HBASJ3 = p_Dbas(JDOFE,3,ICUBP,IEL)

                  ! Finally calculate the contribution to the system
                  ! matrices A11, A12, A21 and A22.

                  AH11 = rconfig%ddeltaTadj * du1loc * HBASI2*HBASJ1
                  AH12 = rconfig%ddeltaTadj * du2loc * HBASI2*HBASJ1
                  AH21 = rconfig%ddeltaTadj * du1loc * HBASI3*HBASJ1
                  AH22 = rconfig%ddeltaTadj * du2loc * HBASI3*HBASJ1

                  ! Weighten the calculated value AHxy by the cubature
                  ! weight OM and add it to the local matrices. After the
                  ! loop over all DOF`s is finished, each entry contains
                  ! the calculated integral.

                  DentryA11(JDOFE,IDOFE,IEL) = DentryA11(JDOFE,IDOFE,IEL)+OM*AH11
                  DentryA12(JDOFE,IDOFE,IEL) = DentryA12(JDOFE,IDOFE,IEL)+OM*AH12
                  DentryA21(JDOFE,IDOFE,IEL) = DentryA21(JDOFE,IDOFE,IEL)+OM*AH21
                  DentryA22(JDOFE,IDOFE,IEL) = DentryA22(JDOFE,IDOFE,IEL)+OM*AH22

                end do ! IDOFE

              end do ! JDOFE

            end do ! ICUBP

          end do ! IEL

        end if

        ! Transposed Newton operator
        !
        !  <tex> $$ dnewtonTransposed * grad(u_1)^T * (.) ) $$ </tex>
        !         = dnewtonTransposed * [ (dx DU1) (dx DU2) ]
        !                               [ (dy DU1) (dy DU2) ]
        !
        ! Should we assemble the transposed Newton matrices?
        if (rconfig%dnewtonT .ne. 0.0_DP) then

          ! Loop over the elements in the current set.
          do IEL=1,IELmax-IELset+1

            ! Loop over all cubature points on the current element
            do ICUBP = 1, ncubp

              ! Calculate the current weighting factor in the cubature formula
              ! in that cubature point.
              !
              ! Normally, we have to take the absolut value of the determinant
              ! of the mapping here!
              ! In 2D, the determinant is always positive, whereas in 3D,
              ! the determinant might be negative -- that is normal!
              ! But because this routine only works in 2D, we can skip
              ! the ABS here!

              OM = p_Domega(ICUBP)*p_Ddetj(ICUBP,IEL)

              ! Current velocity in this cubature point:
              du1locx = DvelocityUderiv (1,ICUBP,IEL)
              du1locy = DvelocityUderiv (2,ICUBP,IEL)
              du2locx = DvelocityVderiv (1,ICUBP,IEL)
              du2locy = DvelocityVderiv (2,ICUBP,IEL)

              ! Outer loop over the DOF`s i=1..indof on our current element,
              ! which corresponds to the basis functions Phi_i:

              do IDOFE=1,indof

                ! Fetch the contributions of the (test) basis functions Phi_i
                ! (our "O")  for function value and first derivatives for the
                ! current DOF into HBASIy:

                HBASI1 = p_Dbas(IDOFE,1,ICUBP,IEL)

                ! Inner loop over the DOF`s j=1..indof, which corresponds to
                ! the basis function Phi_j:

                do JDOFE=1,indof

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! (out "X") for function value and first derivatives for the
                  ! current DOF into HBASJy:

                  HBASJ1 = p_Dbas(JDOFE,1,ICUBP,IEL)

                  ! Finally calculate the contribution to the system
                  ! matrices A11, A12, A21 and A22.
                  !
                  ! With the velocity V=(u,v), we have to assemble:
                  ! grad(V)^t*U, which is realised in each cubature point as:
                  !   du/dx * phi_j*phi_i -> A11
                  !   dv/dx * phi_j*phi_i -> A12
                  !   du/dy * phi_j*phi_i -> A21
                  !   dv/dy * phi_j*phi_i -> A22

                  AH11 = rconfig%dnewtonT * du1locx * HBASJ1*HBASI1
                  AH12 = rconfig%dnewtonT * du2locx * HBASJ1*HBASI1
                  AH21 = rconfig%dnewtonT * du1locy * HBASJ1*HBASI1
                  AH22 = rconfig%dnewtonT * du2locy * HBASJ1*HBASI1

                  ! Weighten the calculated value AHxy by the cubature
                  ! weight OM and add it to the local matrices. After the
                  ! loop over all DOF`s is finished, each entry contains
                  ! the calculated integral.

                  DentryA11(JDOFE,IDOFE,IEL) = DentryA11(JDOFE,IDOFE,IEL)+OM*AH11
                  DentryA12(JDOFE,IDOFE,IEL) = DentryA12(JDOFE,IDOFE,IEL)+OM*AH12
                  DentryA21(JDOFE,IDOFE,IEL) = DentryA21(JDOFE,IDOFE,IEL)+OM*AH21
                  DentryA22(JDOFE,IDOFE,IEL) = DentryA22(JDOFE,IDOFE,IEL)+OM*AH22

                end do ! IDOFE

              end do ! JDOFE

            end do ! ICUBP

          end do ! IEL

        end if

        ! Now do something with the calculated matrices...
        call fincorporate (inonlinComplexity,IELmax-IELset+1,indof,rconfig%dtheta,p_Idofs,&
          DentryA11,DentryA12,DentryA21,DentryA22,rincorporatecollection,Kentry11,Kentry12)

      end do ! IELset

      ! Release the element set
      call elprep_releaseElementSet(revalElementSet)

      ! Release the memory
      deallocate(Dvelocity)
      deallocate(DvelocityUderiv)
      deallocate(DvelocityVderiv)
      deallocate(Dnu)
      deallocate(Dalpha)
      deallocate(DlocalDelta)
      deallocate(Kentry11)
      deallocate(Kentry12)
      deallocate(DentryA11)
      deallocate(DentryA12)
      deallocate(DentryA21)
      deallocate(DentryA22)

      ! Release FE evaluation and cubature information
      call easminfo_doneStdFEBasisEval(rfeBasisEvalData)
      call easminfo_doneStdCubature(rcubatureData)

    end do ! icubatureBlock

    ! Release the assembly structure if necessary.
    if (.not. present(rcubatureInfo)) then
      call spdiscr_releaseCubStructure(rtempCubatureInfo)
    end if

  end subroutine

  ! ----------------------------------------------------------------------

  subroutine getLocalDeltaTriSim (clocalh,&
                      Dvelocity,Dnu,DcubWeights,Ddetj,duMaxR,&
                      cstabiltype,dupsam,Ielements,rtriangulation,Ddelta)

  ! This routine calculates a local ddelta=DELTA_T for a set of finite
  ! elements. This can be used by the streamline diffusion
  ! stabilisation technique as a multiplier of the (local) bilinear form.
  !
  ! Triangular version
  !
  ! Method how to compute the local h.
  ! =0: Use the root of the 2*area of the element as local H
  integer, intent(in) :: clocalH

  ! Array with the values of the velocity field in all cubature points
  ! on all the elements.
  ! dimension(ndim2d, #cubature points per element, #elements)
  real(DP), dimension(:,:,:), intent(in) :: Dvelocity

  ! Reciprocal of the maximum norm of velocity in the domain:
  ! 1/duMaxR = 1/||u||_Omega
  real(DP), intent(in) :: duMaxR

  ! Viscosity coefficient in all cubature points on all selement
  real(DP), dimension(:,:), intent(in) :: Dnu

  ! Cubature weights in the cubature points.
  real(DP), dimension(:), intent(in) :: DcubWeights

  ! Determinant of the mapping from the reference to the real element
  ! in every cubature point.
  real(DP), dimension(:,:), intent(in) :: Ddetj

  ! Type of SD method to apply.
  ! = 0: Use simple SD stabilisation:
  !      ddelta = dupsam * h_T.
  ! = 1: Use Samarskji SD stabilisation; usually dupsam = 0.1 .. 2.
  !      ddelta = dupsam * h_t/||u||_T * 2*Re_T/(1+Re_T)
  integer :: cstabilType

  ! user defined parameter for configuring the streamline diffusion.
  real(DP), intent(in) :: dupsam

  ! List of elements where to calculate the local delta.
  integer, dimension(:), intent(in) :: Ielements

  ! Triangulation that defines the mesh.
  type(t_triangulation), intent(in) :: rtriangulation

  ! Out: local Ddelta on all elements
  real(DP), dimension(:), intent(out) :: Ddelta

  ! local variables
  real(DP) :: dlocalH,du1,du2,dunorm,dreLoc,dnuRec
  integer :: iel,ielidx,icubp
  real(DP), dimension(:), pointer :: p_DelementVolume
  real(DP) :: domega

    ! Currently, we only support clocalh=0!

    call storage_getbase_double (rtriangulation%h_DelementVolume,p_DelementVolume)

    ! Loop through all elements
    do ielidx = 1,size(Ielements)

      ! Get the element number
      iel = Ielements(ielidx)

      ! Loop through the cubature points on the current element
      ! and calculate the mean velocity there.
      du1 = 0.0_DP
      du2 = 0.0_DP
      dnuRec = 0.0_DP
      do icubp = 1,ubound(Dvelocity,2)
        domega = DcubWeights(icubp)*Ddetj(icubp,ielidx)
        du1 = du1 + domega*Dvelocity(1,icubp,ielidx)
        du2 = du2 + domega*Dvelocity(2,icubp,ielidx)
        dnuRec = dnuRec + domega*Dnu(icubp,ielidx)
      end do

      ! Calculate the norm of the mean local velocity
      ! as well as the mean local velocity
      du1 = du1 / p_DelementVolume(iel)
      du2 = du2 / p_DelementVolume(iel)
      dunorm = sqrt(du1**2+du2**2)

      ! Calculate the mean viscosity coefficient -- or more precisely,
      ! its reciprocal.
      dnuRec = p_DelementVolume(iel) / dnuRec

      ! Now we have:   dunorm = ||u||_T
      ! and:           u_T = a1*u1_T + a2*u2_T

      ! If the norm of the velocity is small, we choose ddelta = 0,
      ! which results in central difference in the streamline diffusion
      ! matrix assembling:

      if (dunorm .le. 1E-8_DP) then

        Ddelta(ielidx) = 0.0_DP

      else

        ! Calculate the local h from the area of the element.
        ! As the element is a tri, multiply the volume by 2.
        dlocalH = sqrt(2.0_DP*p_DelementVolume(iel))

        ! Calculate ddelta... (cf. p. 121 in Turek`s CFD book)

        if (cstabiltype .eq. 0) then

          ! For UPSAM<0, we use simple calculation of ddelta:

          Ddelta(ielidx) = abs(dupsam)*dlocalH

        else

          ! For UPSAM >= 0, we use standard Samarskji-like calculation
          ! of ddelta. At first calculate the local Reynolds number
          ! RELOC = Re_T = ||u||_T * h_T / NU

          dreLoc = dunorm*dlocalH*dnuRec

          ! and then the ddelta = UPSAM * h_t/||u|| * 2*Re_T/(1+Re_T)

          Ddelta(ielidx) = dupsam * dlocalH*duMaxR * 2.0_DP*(dreLoc/(1.0_DP+dreLoc))

        end if ! (UPSAM.LT.0.0)

      end if ! (dunorm.LE.1D-8)

    end do

  end subroutine

  ! ----------------------------------------------------------------------

  subroutine getLocalDeltaQuadSim (clocalh,&
                      Dvelocity,Dnu,DcubWeights,Ddetj,duMaxR,&
                      cstabiltype,dupsam,Ielements,rtriangulation,Ddelta)

  ! This routine calculates a local ddelta=DELTA_T for a set of finite
  ! elements. This can be used by the streamline diffusion
  ! stabilisation technique as a multiplier of the (local) bilinear form.
  !
  ! Method how to compute the local h.
  ! =0: Use the root of the area of the element as local H
  ! =1: Use the length of the way that a particle travels through
  !     the element in direction of the flow
  integer, intent(in) :: clocalH

  ! Array with the values of the velocity field in all cubature points
  ! on all the elements.
  ! dimension(ndim2d, #cubature points per element, #elements)
  real(DP), dimension(:,:,:), intent(in) :: Dvelocity

  ! Reciprocal of the maximum norm of velocity in the domain:
  ! 1/duMaxR = 1/||u||_Omega
  real(DP), intent(in) :: duMaxR

  ! Viscosity coefficient in all cubature points on all selement
  real(DP), dimension(:,:), intent(in) :: Dnu

  ! Cubature weights in the cubature points on the reference element
  real(DP), dimension(:), intent(in) :: DcubWeights

  ! Determinant of the mapping from the reference to the real element
  ! in every cubature point.
  real(DP), dimension(:,:), intent(in) :: Ddetj

  ! Type of SD method to apply.
  ! = 0: Use simple SD stabilisation:
  !      ddelta = dupsam * h_T.
  ! = 1: Use Samarskji SD stabilisation; usually dupsam = 0.1 .. 2.
  !      ddelta = dupsam * h_t/||u||_T * 2*Re_T/(1+Re_T)
  integer :: cstabilType

  ! user defined parameter for configuring the streamline diffusion.
  real(DP), intent(in) :: dupsam

  ! List of elements where to calculate the local delta.
  integer, dimension(:), intent(in) :: Ielements

  ! Triangulation that defines the mesh.
  type(t_triangulation), intent(in) :: rtriangulation

  ! Out: local Ddelta on all elements
  real(DP), dimension(:), intent(out) :: Ddelta

  ! local variables
  real(DP) :: dlocalH,du1,du2,dunorm,dreLoc,dnuRec
  integer :: iel,ielidx,icubp
  integer, dimension(:,:), pointer :: p_IverticesAtElement
  real(DP), dimension(:,:), pointer :: p_DvertexCoords
  real(DP), dimension(:), pointer :: p_DelementVolume
  real(DP) :: domega

    call storage_getbase_double (rtriangulation%h_DelementVolume,p_DelementVolume)

    ! Get some crucial data
    if (clocalh .eq. 0) then

      ! Loop through all elements
      do ielidx = 1,size(Ielements)

        ! Get the element number
        iel = Ielements(ielidx)

        ! Loop through the cubature points on the current element
        ! and calculate the mean velocity there.
        du1 = 0.0_DP
        du2 = 0.0_DP
        dnuRec = 0.0_DP
        do icubp = 1,ubound(Dvelocity,2)
          domega = DcubWeights(icubp)*Ddetj(icubp,ielidx)
          du1 = du1 + domega*Dvelocity(1,icubp,ielidx)
          du2 = du2 + domega*Dvelocity(2,icubp,ielidx)
          dnuRec = dnuRec + domega*Dnu(icubp,ielidx)
        end do

        ! Calculate the norm of the mean local velocity
        ! as well as the mean local velocity
        du1 = du1 / p_DelementVolume(iel)
        du2 = du2 / p_DelementVolume(iel)
        dunorm = sqrt(du1**2+du2**2)

        ! Calculate the mean viscosity coefficient -- or more precisely,
        ! its reciprocal.
        dnuRec = p_DelementVolume(iel) / dnuRec

        ! Now we have:   dunorm = ||u||_T
        ! and:           u_T = a1*u1_T + a2*u2_T

        ! If the norm of the velocity is small, we choose ddelta = 0,
        ! which results in central difference in the streamline diffusion
        ! matrix assembling:

        if (dunorm .le. 1E-8_DP) then

          Ddelta(ielidx) = 0.0_DP

        else

          ! Calculate the local h from the area of the element
          dlocalH = sqrt(p_DelementVolume(iel))

          ! Calculate ddelta... (cf. p. 121 in Turek`s CFD book)

          if (cstabiltype .eq. 0) then

            ! For UPSAM<0, we use simple calculation of ddelta:

            Ddelta(ielidx) = abs(dupsam)*dlocalH

          else

            ! For UPSAM >= 0, we use standard Samarskji-like calculation
            ! of ddelta. At first calculate the local Reynolds number
            ! RELOC = Re_T = ||u||_T * h_T / NU

            dreLoc = dunorm*dlocalH*dnuRec

            ! and then the ddelta = UPSAM * h_t/||u|| * 2*Re_T/(1+Re_T)

            Ddelta(ielidx) = dupsam * dlocalH*duMaxR * 2.0_DP*(dreLoc/(1.0_DP+dreLoc))

          end if ! (UPSAM.LT.0.0)

        end if ! (dunorm.LE.1D-8)

      end do

    else

      call storage_getbase_double2d (rtriangulation%h_DvertexCoords,p_DvertexCoords)
      call storage_getbase_int2d (rtriangulation%h_IverticesAtElement,p_IverticesAtElement)

      ! Loop through all elements
      do ielidx = 1,size(Ielements)

        ! Get the element number
        iel = Ielements(ielidx)

        ! Loop through the cubature points on the current element
        ! and calculate the mean velocity there.
        du1 = 0.0_DP
        du2 = 0.0_DP
        dnuRec = 0.0_DP
        do icubp = 1,ubound(Dvelocity,2)
          domega = DcubWeights(icubp)*Ddetj(icubp,ielidx)
          du1 = du1 + domega*Dvelocity(1,icubp,ielidx)
          du2 = du2 + domega*Dvelocity(2,icubp,ielidx)
          dnuRec = dnuRec + domega*Dnu(icubp,ielidx)
        end do

        ! Calculate the norm of the mean local velocity
        ! as well as the mean local velocity
        du1 = du1 / p_DelementVolume(iel)
        du2 = du2 / p_DelementVolume(iel)
        dunorm = sqrt(du1**2+du2**2)

        ! Calculate the mean viscosity coefficient -- or more precisely,
        ! its reciprocal.
        dnuRec = p_DelementVolume(iel) / dnuRec

        ! Now we have:   dunorm = ||u||_T
        ! and:           u_T = a1*u1_T + a2*u2_T

        ! If the norm of the velocity is small, we choose ddelta = 0,
        ! which results in central difference in the streamline diffusion
        ! matrix assembling:

        if (dunorm .le. 1E-8_DP) then

          Ddelta(ielidx) = 0.0_DP

        else

          ! u_T defines the "slope" of the velocity through
          ! the element T. At next, calculate the local mesh width
          ! dlocalH = h = h_T on our element T=IEL:

          call getLocalMeshWidthQuad (dlocalH,dunorm, du1, du2, iel, &
              p_IverticesAtElement,p_DvertexCoords)

          ! Calculate ddelta... (cf. p. 121 in Turek`s CFD book)

          if (cstabiltype .eq. 0) then

            ! For UPSAM<0, we use simple calculation of ddelta:

            Ddelta(ielidx) = abs(dupsam)*dlocalH

          else

            ! For UPSAM >= 0, we use standard Samarskji-like calculation
            ! of ddelta. At first calculate the local Reynolds number
            ! RELOC = Re_T = ||u||_T * h_T / NU

            dreLoc = dunorm*dlocalH*dnuRec

            ! and then the ddelta = UPSAM * h_t/||u|| * 2*Re_T/(1+Re_T)

            Ddelta(ielidx) = dupsam * dlocalH*duMaxR * 2.0_DP*(dreLoc/(1.0_DP+dreLoc))

          end if ! (UPSAM.LT.0.0)

        end if ! (dunorm.LE.1D-8)

      end do

    end if

  end subroutine

  ! ---------------------------------------------------------------------------

  subroutine conv_sdIncorpToMatrix2D (inonlinComplexity,nelements,indof,dtheta,Idofs,&
      DentryA11,DentryA12,DentryA21,DentryA22,rcollection,KentryA11,KentryA12)

  ! Callback routine. Incorporates local matrices into a global matrix.
  ! The global matrix must be referred by rcollection%p_rmatrixQuickAccess1.

  ! Complexity of the matrix
  integer, intent(in) :: inonlinComplexity

  ! Number of elements and vertices on each element
  integer, intent(in) :: nelements,indof

  ! Weight for the local matrices
  real(DP), intent(in) :: dtheta

  ! Dummy array. Receives positions of the local matrices in
  ! the global matrix for A11/A22 and A12/A21.
  ! Can be undefined upon entering the routine.
  integer, dimension(:,:,:), intent(inout) :: KentryA11,KentryA12

  ! Values of the local matrices.
  real(DP), dimension(:,:,:), intent(in) :: DentryA11,DentryA12,DentryA21,DentryA22

  ! The DOF`s on all elements the routine should work on
  integer, dimension(:,:), intent(in) :: Idofs

  ! Collection structure. p_rmatrixQuickAccess1 points to the matrix
  ! where to incorporate the data.
  type(t_collection), intent(in) :: rcollection

    ! local variables
    integer :: iel,idofe,jdofe

    ! Matrix structure arrays.
    ! The matrix structure of A11 and A22 must be the same.
    ! The matrix structure of A12 and A21 must be the same.
    ! However, A11 and A12 may have different structure.
    real(DP), dimension(:), pointer :: p_Da11,p_Da12,p_Da21,p_Da22

    ! We assemble local matrices and incorporate them later into the
    ! global matrix.
    ! Calculate the positions of the local matrix in the global matrix.
    call bilf_getLocalMatrixIndices (rcollection%p_rmatrixQuickAccess1%RmatrixBlock(1,1),&
        Idofs,Idofs,KentryA11,indof,indof,nelements)

    ! Now incorporate the local system matrices into the global one
    select case (inonlinComplexity)
    case (0,1)
      ! Get the matrix data
      call lsyssc_getbase_double (rcollection%p_rmatrixQuickAccess1%RmatrixBlock(1,1),&
          p_Da11)

      ! There is only data for A11
      do IEL=1,nelements
        do IDOFE=1,indof
          do JDOFE=1,indof
            p_Da11(KentryA11(JDOFE,IDOFE,IEL)) = p_Da11(KentryA11(JDOFE,IDOFE,IEL)) + &
                dtheta * DentryA11(JDOFE,IDOFE,IEL)
          end do
        end do
      end do

    case (2)
      ! Get the matrix data
      call lsyssc_getbase_double (rcollection%p_rmatrixQuickAccess1%RmatrixBlock(1,1),&
          p_Da11)
      call lsyssc_getbase_double (rcollection%p_rmatrixQuickAccess1%RmatrixBlock(2,2),&
          p_Da22)

      ! Include all A11, A12

      do IEL=1,nelements
        do IDOFE=1,indof
          do JDOFE=1,indof
            ! Kentry (:,:,:) -> positions of local matrix in A11 and A22.
            !
            ! DentryA11 (:,:,:) -> part of A11
            p_Da11(KentryA11(JDOFE,IDOFE,IEL)) = p_Da11(KentryA11(JDOFE,IDOFE,IEL)) + &
                dtheta * DentryA11(JDOFE,IDOFE,IEL)

            ! DentryA22 (:,:,:) -> part of A22
            p_Da22(KentryA11(JDOFE,IDOFE,IEL)) = p_Da22(KentryA11(JDOFE,IDOFE,IEL)) + &
                dtheta * DentryA22(JDOFE,IDOFE,IEL)
          end do
        end do
      end do

    case default
      ! A12/A21 is there
      call bilf_getLocalMatrixIndices (rcollection%p_rmatrixQuickAccess1%RmatrixBlock(1,2),&
          Idofs,Idofs,KentryA12,indof,indof,nelements)

      ! Get the matrix data
      call lsyssc_getbase_double (rcollection%p_rmatrixQuickAccess1%RmatrixBlock(1,1),&
          p_Da11)
      call lsyssc_getbase_double (rcollection%p_rmatrixQuickAccess1%RmatrixBlock(1,2),&
          p_Da12)
      call lsyssc_getbase_double (rcollection%p_rmatrixQuickAccess1%RmatrixBlock(2,1),&
          p_Da21)
      call lsyssc_getbase_double (rcollection%p_rmatrixQuickAccess1%RmatrixBlock(2,2),&
          p_Da22)

      ! Include all Aij

      do IEL=1,nelements
        do IDOFE=1,indof
          do JDOFE=1,indof
            ! Kentry (:,:,:) -> positions of local matrix in A11 and A22.
            !
            ! DentryA11 (:,:,:) -> Newton part of A11
            p_Da11(KentryA11(JDOFE,IDOFE,IEL)) = p_Da11(KentryA11(JDOFE,IDOFE,IEL)) + &
                dtheta * DentryA11(JDOFE,IDOFE,IEL)

            ! DentryA22 (:,:,:) -> Newton part of A22
            p_Da22(KentryA11(JDOFE,IDOFE,IEL)) = p_Da22(KentryA11(JDOFE,IDOFE,IEL)) + &
                dtheta * DentryA22(JDOFE,IDOFE,IEL)

            ! DentryA12 (:,:,:) -> Newton part of A12
            p_Da12(KentryA12(JDOFE,IDOFE,IEL)) = p_Da12(KentryA12(JDOFE,IDOFE,IEL)) + &
                dtheta * DentryA12(JDOFE,IDOFE,IEL)

            ! Dentry21 (:,:,:) -> Newton part of A21
            p_Da21(KentryA12(JDOFE,IDOFE,IEL)) = p_Da21(KentryA12(JDOFE,IDOFE,IEL)) + &
                dtheta * DentryA21(JDOFE,IDOFE,IEL)
          end do
        end do
      end do

    end select

  end subroutine

  ! ---------------------------------------------------------------------------

  subroutine conv_sdIncorpToDefect2D (inonlinComplexity,nelements,indof,dtheta,Idofs,&
      DentryA11,DentryA12,DentryA21,DentryA22,rcollection,KentryA11,KentryA12)

  ! Callback routine. Incorporates local matrices into a defect vector.
  ! The solution vector must be referred to by rcollection%p_rvectorQuickAccess1.
  ! The RHS and destination defect vector must be referred to by
  ! rcollection%p_rvectorQuickAccess2.

  ! Complexity of the matrix
  integer, intent(in) :: inonlinComplexity

  ! Number of elements and vertices on each element
  integer, intent(in) :: nelements,indof

  ! Weight for the local matrices
  real(DP), intent(in) :: dtheta

  ! Positions of the local matrices in the global matrix for A11/A22 and A12/A21.
  integer, dimension(:,:,:), intent(inout) :: KentryA11,KentryA12

  ! Values of the local matrices.
  real(DP), dimension(:,:,:), intent(in) :: DentryA11,DentryA12,DentryA21,DentryA22

  ! The DOF`s on all elements the routine should work on
  integer, dimension(:,:), intent(in) :: Idofs

  ! Collection structure. p_rmatrixQuickAccess1 points to the matrix
  ! where to incorporate the data.
  type(t_collection), intent(in) :: rcollection

    ! local variables
    integer :: iel,idofe,jdofe,jdfg,idfg

    ! Vector arrays.
    real(DP), dimension(:), pointer :: p_Dd1,p_Dd2,p_Dx1,p_Dx2

    ! Get the vector data
    call lsyssc_getbase_double (rcollection%p_rvectorQuickAccess1%RvectorBlock(1),&
        p_Dx1)
    call lsyssc_getbase_double (rcollection%p_rvectorQuickAccess1%RvectorBlock(2),&
        p_Dx2)
    call lsyssc_getbase_double (rcollection%p_rvectorQuickAccess2%RvectorBlock(1),&
        p_Dd1)
    call lsyssc_getbase_double (rcollection%p_rvectorQuickAccess2%RvectorBlock(2),&
        p_Dd2)

    ! Now incorporate the local system matrices into the global one
    select case (inonlinComplexity)
    case (0,1)

      ! There is only data for A11
      do IEL=1,nelements
        do IDOFE=1,indof
          idfg=Idofs(IDOFE,IEL)
          do JDOFE=1,indof
            jdfg =Idofs(JDOFE,IEL)
            p_Dd1(idfg) = p_Dd1(idfg) - dtheta * DentryA11(JDOFE,IDOFE,IEL) * p_Dx1(jdfg)
            p_Dd2(idfg) = p_Dd2(idfg) - dtheta * DentryA11(JDOFE,IDOFE,IEL) * p_Dx2(jdfg)
          end do
        end do
      end do

    case (2)
      ! Include all A11, A12

      do IEL=1,nelements
        do IDOFE=1,indof
          idfg=Idofs(IDOFE,IEL)
          do JDOFE=1,indof
            jdfg =Idofs(JDOFE,IEL)
            p_Dd1(idfg) = p_Dd1(idfg) - dtheta * DentryA11(JDOFE,IDOFE,IEL) * p_Dx1(jdfg)
            p_Dd2(idfg) = p_Dd2(idfg) - dtheta * DentryA22(JDOFE,IDOFE,IEL) * p_Dx2(jdfg)
          end do
        end do
      end do

    case default
      ! Include all Aij

      do IEL=1,nelements
        do IDOFE=1,indof
          idfg=Idofs(IDOFE,IEL)
          do JDOFE=1,indof
            jdfg =Idofs(JDOFE,IEL)
            p_Dd1(idfg) = p_Dd1(idfg) - &
              dtheta * (DentryA11(JDOFE,IDOFE,IEL) * p_Dx1(jdfg) + &
                        DentryA12(JDOFE,IDOFE,IEL) * p_Dx2(jdfg) )
            p_Dd2(idfg) = p_Dd2(idfg) - &
              dtheta * (DentryA21(JDOFE,IDOFE,IEL) * p_Dx1(jdfg) + &
                        DentryA22(JDOFE,IDOFE,IEL) * p_Dx2(jdfg) )
          end do
        end do
      end do

    end select

  end subroutine

end module

