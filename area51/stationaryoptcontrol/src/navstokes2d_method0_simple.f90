!##############################################################################
!# ****************************************************************************
!# <name> navstokes2d_method0_simple </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a demonstration program how to solve a simple Navier-Stokes
!# optimal control problem with constant coefficients on a simple domain.
!# </purpose>
!#
!# Problem and system definition
!# -----------------------------
!# The system to solve is the following:
!#
!#                                 -Laplace(y) + y grad y + grad(p) = f + u
!#                                                           -div y = 0
!#  -Laplace(lambda) + y grad lambda + (grad y)^t lambda + grad(xi) = y - z
!#                                                      -div lambda = 0
!#                                                                u = P(-1/alpha l)
!#
!# with:
!#   u/p = primal solution
!#   l/xi = dual solution
!#   z = target velocity
!#
!# The boundary conditions for the dual solution are =0 on the
!# Dirichlet boundary. 
!# The projection operator
!#
!#   P(g) = P(g)g = P[a,b](g)g = max(a,min(b,g))
!#
!# with
!#
!#   P(h)g = a     , for all x with h(x) < a
!#         = g(x)  , for all x with a <= h(x) <= b
!#         = b     , for all x with h(x) > b
!#
!# for a function g=g(x) defines the projection to the interval [a,b].
!#
!##############################################################################

module navstokes2d_method0_simple

  use fsystem
  use genoutput
  use storage
  use boundary
  use cubature
  use linearalgebra
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use ucd
  use pprocerror
  use genoutput
  use stdoperators
  use spdiscprojection
  use convection
  use analyticprojection
  use optcontrolconvection
  use matrixio
  use collection
  use filtersupport
  use scalarpde
  use bilinearformevaluation
  use linearformevaluation
  use multilevelprojection
  use linearsolver
  use matrixmodification
  use feevaluation
  use vectorio
    
  use navstokes2d_callback
  
  implicit none


  ! Encapsules a target flow

  type t_targetFlow
  
    ! Type of the target flow.
    ! =0: Poiseuille-flow
    ! =1: Precalculated solution
    integer :: itype
    
    ! Triangulation of the target flow. Underlying domain is the 
    ! same as in the main problem.
    type(t_triangulation) :: rtriangulation
    
    ! Underlying discretisation
    type(t_blockDiscretisation) :: rdiscretisation
  
    ! Block vector encapsuling the flow
    type(t_vectorBlock) :: rvector
  
  end type

 ! Information on one level
  type t_level
  
    ! Triangulation
    type(t_triangulation) :: rtriangulation
    
    ! Discretisation
    type(t_blockDiscretisation) :: rdiscretisation
  
    ! Laplace matrix of that level
    type(t_matrixScalar) :: rmatrixLaplace
    
    ! Mass matrix
    type(t_matrixScalar) :: rmatrixMass
    
    ! B1-matrix
    type(t_matrixScalar) :: rmatrixB1
    
    ! B2-matrix
    type(t_matrixScalar) :: rmatrixB2
    
    ! Diagonal pressure matrix
    type(t_matrixScalar) :: rmatrixDiagP
    
    ! Boundary conditions
    type(t_discreteBC) :: rdiscreteBC
    
    ! Projection structure for MG solvers
    type(t_interlevelProjectionBlock) :: rprojection
    
    type(t_vectorBlock) :: rtempVector

  end type
  
  ! Configuration of a matrix
  type t_matrixConfig

    ! Diffusion parameter
    real(DP) :: dnu
    
    ! Emulate one timestep of the nonstationary Navier-Stokes.
    logical :: bemulateTimestep
    
    ! Length of the timestep
    real(dp) :: dt
    
    ! Stabilisation parameter
    real(dp) :: dupsamPrimal

    ! Stabilisation parameter
    real(dp) :: dupsamDual
    
    ! Pure Dirichlet problem
    logical :: bpureDirichlet

    ! Couple the dual equation to the primal one
    logical :: bdualcoupledtoprimal
    
    ! Calculate Navier-Stokes
    logical :: bnavierStokes
    
    ! Regularisation parameter for the control
    real(dp) :: dalpha

    ! Calculate the Newton matrix
    logical :: bnewton
    
    ! Activate the effect of the control
    logical :: bcontrolactive
    
    ! Activate constraints in the control
    logical :: bboundsActive
    
    ! Use exact derivatives for implementing constraints
    logical :: bexactderiv
    
    ! Min/Max bound for X-velocity
    real(dp) :: dmin1,dmax1
    
    ! Min/Max bound for Y-velocity
    real(dp) :: dmin2,dmax2

  end type

contains

  ! ***************************************************************************

!<subroutine>

  subroutine calcStandardMatrices (rlevel)
  
!<description>
  ! Calculate the basic matrices on that level.
!</description>

!<inputoutput>
  ! Information for that level
  type(t_level), intent(inout) :: rlevel
!</inputoutput>

!</subroutine>

    ! Now as the discretisation is set up, we can start to generate
    ! the structure of the system matrix which is to solve.
    ! We create a scalar matrix, based on the discretisation structure
    ! for our one and only solution component.
    call bilf_createMatrixStructure (rlevel%rdiscretisation%RspatialDiscr(1),&
                                     LSYSSC_MATRIX9,rlevel%rmatrixLaplace)
                                     
    ! Assemble the Laplace matrix.
    call stdop_assembleLaplaceMatrix (rlevel%rmatrixLaplace)
    
    ! And a mass matrix; we need it for the coupling.
    ! Share the structure with the Laplace matrix.
    call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rlevel%rmatrixMass,&
        LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
    call stdop_assembleSimpleMatrix (rlevel%rmatrixMass,DER_FUNC,DER_FUNC)

    ! Create B-matrices
    call bilf_createMatrixStructure (rlevel%rdiscretisation%RspatialDiscr(3),&
        LSYSSC_MATRIX9,rlevel%rmatrixB1,rlevel%rdiscretisation%RspatialDiscr(1))
    call lsyssc_duplicateMatrix (rlevel%rmatrixB1,rlevel%rmatrixB2,&
        LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)

    ! Diagonal pressure matrix
    call bilf_createMatrixStructure (rlevel%rdiscretisation%RspatialDiscr(3),&
        LSYSSC_MATRIX9,rlevel%rmatrixDiagP)
        
    call stdop_assembleSimpleMatrix (rlevel%rmatrixB1,DER_FUNC,DER_DERIV_X,-1.0_DP,.true.)
    call stdop_assembleSimpleMatrix (rlevel%rmatrixB2,DER_FUNC,DER_DERIV_Y,-1.0_DP,.true.)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getDefect (rdefectBlock,rvectorBlock,rrhsBlock,rparams,rlevel)
    
!<description>
  ! Calculate the defect.
!</description>
  
!<inputoutput>
  ! The defect vector to create
  type(t_vectorBlock), intent(inout) :: rdefectBlock
!</inputoutput>

!<input>
  ! Solution vector
  type(t_vectorBlock), intent(inout) :: rvectorBlock

  ! RHS vector
  type(t_vectorBlock), intent(inout) :: rrhsBlock

  ! Parameter that configure the matrix
  type(t_matrixConfig), intent(in) :: rparams

  ! Information for that level
  type(t_level), intent(in) :: rlevel
!</input>

!</subroutine>

    ! local variables
    type(t_vectorScalar) :: rvectorTmp
    type(t_matrixBlock) :: rtempMatrix
    type(t_vectorBlock) :: rtempVectorSol,rtempVectorB,rtempVectorX
    type(t_convStreamlineDiffusion) :: rstreamlineDiffPrimalDef,rstreamlineDiffDualDef
    real(dp), dimension(:), pointer :: p_DdataX,p_DdataY,p_DdataP
    type(t_optcoperator) :: roptcoperator
    
    logical, parameter :: bnewmethod = .true.
    
    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Dx,p_Dy
    call lsysbl_getbase_double (rvectorBlock,p_Dx)
    call lsysbl_getbase_double (rrhsBlock,p_Dy)
    
    if (bnewmethod) then

      ! Start with the RHS...
      call lsysbl_copyVector (rrhsBlock,rdefectBlock)

      ! Matrix-vector with the B-matrices
      call lsyssc_scalarMatVec (rlevel%rmatrixB1,&
          rvectorBlock%RvectorBlock(3),rdefectBlock%RvectorBlock(1),-1.0_DP,1.0_DP)
      call lsyssc_scalarMatVec (rlevel%rmatrixB2,&
          rvectorBlock%RvectorBlock(3),rdefectBlock%RvectorBlock(2),-1.0_DP,1.0_DP)

      call lsyssc_scalarMatVec (rlevel%rmatrixB1,&
          rvectorBlock%RvectorBlock(1),rdefectBlock%RvectorBlock(3),-1.0_DP,1.0_DP,&
          .true.)
      call lsyssc_scalarMatVec (rlevel%rmatrixB2,&
          rvectorBlock%RvectorBlock(2),rdefectBlock%RvectorBlock(3),-1.0_DP,1.0_DP,&
          .true.)
      
      call lsyssc_scalarMatVec (rlevel%rmatrixB1,&
          rvectorBlock%RvectorBlock(6),rdefectBlock%RvectorBlock(4),-1.0_DP,1.0_DP)
      call lsyssc_scalarMatVec (rlevel%rmatrixB2,&
          rvectorBlock%RvectorBlock(6),rdefectBlock%RvectorBlock(5),-1.0_DP,1.0_DP)

      call lsyssc_scalarMatVec (rlevel%rmatrixB1,&
          rvectorBlock%RvectorBlock(4),rdefectBlock%RvectorBlock(6),-1.0_DP,1.0_DP,&
          .true.)
      call lsyssc_scalarMatVec (rlevel%rmatrixB2,&
          rvectorBlock%RvectorBlock(5),rdefectBlock%RvectorBlock(6),-1.0_DP,1.0_DP,&
          .true.)

      ! Initialise the operator structure for what we need.
      roptcoperator%dupsamPrimal = rparams%dupsamPrimal
      roptcoperator%dupsamDual = rparams%dupsamDual
      
      ! Timestep-weights
      if (rparams%bemulateTimestep) then
        roptcoperator%dprimalAlpha = 1.0_DP/rparams%dt
        roptcoperator%ddualAlpha   = 1.0_DP/rparams%dt
      end if

      ! Stokes operator
      roptcoperator%dnu = rparams%dnu
      roptcoperator%dprimalBeta = 1.0_DP
      roptcoperator%ddualBeta   = 1.0_DP
      
      ! Nonlinearity
      if (rparams%bnavierStokes) then
        roptcoperator%dprimalDelta =  1.0_DP
        roptcoperator%ddualDelta   = -1.0_DP
        roptcoperator%ddualNewtonTrans = 1.0_DP
        
        ! Whether or not Newton is active has no influence to the
        ! defect, so the following lines are commented out.
        ! if (rparams%bnewton) then
        !   roptcoperator%dprimalNewton    = 1.0_DP
        !   roptcoperator%ddualRDeltaTrans = 1.0_DP
        !   roptcoperator%ddualRNewton     = -1.0_DP
        ! end if
        
      end if
      
      ! Coupling matrices
      if (rparams%bdualcoupledtoprimal) then
        roptcoperator%ddualRAlpha = -1.0_DP
      end if

      if (rparams%bcontrolactive) then
        roptcoperator%dcontrolWeight = -1.0_DP
        roptcoperator%dcontrolMultiplier = -1.0_DP/rparams%dalpha
      end if
      
      if (rparams%bboundsactive) then
        if (.not. rparams%bexactderiv) then
          roptcoperator%ccontrolProjection = 2
        else
          roptcoperator%ccontrolProjection = 1
        end if
        roptcoperator%dmin1 = rparams%dmin1
        roptcoperator%dmax1 = rparams%dmax1
        roptcoperator%dmin2 = rparams%dmin2
        roptcoperator%dmax2 = rparams%dmax2
      end if
      
      ! Calculate the velocity-dependent part of the system matrix.
      call conv_strdiffOptC2dgetDefect (rlevel%rmatrixMass,roptcoperator,&
          rvectorBlock,rvectorBlock,&
          1.0_DP,rvectorBlock,rdefectBlock)
      
      !call vecio_writeBlockVectorHR (rdefectblock, 'vector', .false.,&
      !                               0, 'vector.txt', '(E12.5)')
      
    else    
      ! DEBUG!!!
      call lsyssc_getbase_double (rdefectBlock%RvectorBlock(1),p_DdataX)
      call lsyssc_getbase_double (rdefectBlock%RvectorBlock(2),p_DdataY)
      call lsyssc_getbase_double (rdefectBlock%RvectorBlock(3),p_DdataP)
      
      ! Prepare the streamline diffusion structure for assembling the
      ! nonlinearities.
      !
      ! Primal Defect
      rstreamlineDiffPrimalDef%dupsam = rparams%dupsamPrimal
      rstreamlineDiffPrimalDef%dnu = rparams%dnu
      rstreamlineDiffPrimalDef%ddelta = 1.0_DP
      rstreamlineDiffPrimalDef%dnewton = 0.0_DP

      ! Dual Defect
      rstreamlineDiffDualDef%dupsam = rparams%dupsamDual
      rstreamlineDiffDualDef%dnu = rparams%dnu
      rstreamlineDiffDualDef%ddelta = -1.0_DP
      rstreamlineDiffDualDef%dnewtonTransposed = 1.0_DP

      ! *************************************************************
      ! Create the nonlinear defect
      !  (d1) =  (  f ) - (  A    -P(u)(-1/alpha .) ) ( y ) 
      !  (d2)    ( -z )   ( -M    A                 ) ( p ) 
      ! We do this manually...
      
      call lsysbl_copyVector (rvectorBlock,rdefectBlock)
      call lsysbl_copyVector (rrhsBlock,rdefectBlock)
      
      call lsyssc_scalarMatVec (rlevel%rmatrixLaplace,&
          rvectorBlock%RvectorBlock(1),rdefectBlock%RvectorBlock(1),-rparams%dnu,1.0_DP)
      call lsyssc_scalarMatVec (rlevel%rmatrixLaplace,&
          rvectorBlock%RvectorBlock(2),rdefectBlock%RvectorBlock(2),-rparams%dnu,1.0_DP)
    
      if (rparams%bemulateTimestep) then
        ! One timestep from zero solution to t=1. Emulated by adding a mass matrix
        ! to the Laplace.
        call lsyssc_scalarMatVec (rlevel%rmatrixMass,&
            rvectorBlock%RvectorBlock(1),rdefectBlock%RvectorBlock(1),-1.0_DP/rparams%dt,1.0_DP)
        call lsyssc_scalarMatVec (rlevel%rmatrixMass,&
            rvectorBlock%RvectorBlock(2),rdefectBlock%RvectorBlock(2),-1.0_DP/rparams%dt,1.0_DP)
      end if

      call lsyssc_scalarMatVec (rlevel%rmatrixB1,&
          rvectorBlock%RvectorBlock(3),rdefectBlock%RvectorBlock(1),-1.0_DP,1.0_DP)
      call lsyssc_scalarMatVec (rlevel%rmatrixB2,&
          rvectorBlock%RvectorBlock(3),rdefectBlock%RvectorBlock(2),-1.0_DP,1.0_DP)

      call lsyssc_scalarMatVec (rlevel%rmatrixB1,&
          rvectorBlock%RvectorBlock(1),rdefectBlock%RvectorBlock(3),-1.0_DP,1.0_DP,&
          .true.)
      call lsyssc_scalarMatVec (rlevel%rmatrixB2,&
          rvectorBlock%RvectorBlock(2),rdefectBlock%RvectorBlock(3),-1.0_DP,1.0_DP,&
          .true.)
      
      call lsyssc_scalarMatVec (rlevel%rmatrixLaplace,&
          rvectorBlock%RvectorBlock(4),rdefectBlock%RvectorBlock(4),-rparams%dnu,1.0_DP)
      call lsyssc_scalarMatVec (rlevel%rmatrixLaplace,&
          rvectorBlock%RvectorBlock(5),rdefectBlock%RvectorBlock(5),-rparams%dnu,1.0_DP)

      if (rparams%bemulateTimestep) then
        ! One timestep from zero solution to t=1. Emulated by adding a mass matrix
        ! to the Laplace.
        call lsyssc_scalarMatVec (rlevel%rmatrixMass,&
            rvectorBlock%RvectorBlock(4),rdefectBlock%RvectorBlock(4),-1.0_DP/rparams%dt,1.0_DP)
        call lsyssc_scalarMatVec (rlevel%rmatrixMass,&
            rvectorBlock%RvectorBlock(5),rdefectBlock%RvectorBlock(5),-1.0_DP/rparams%dt,1.0_DP)
      end if

      call lsyssc_scalarMatVec (rlevel%rmatrixB1,&
          rvectorBlock%RvectorBlock(6),rdefectBlock%RvectorBlock(4),-1.0_DP,1.0_DP)
      call lsyssc_scalarMatVec (rlevel%rmatrixB2,&
          rvectorBlock%RvectorBlock(6),rdefectBlock%RvectorBlock(5),-1.0_DP,1.0_DP)

      call lsyssc_scalarMatVec (rlevel%rmatrixB1,&
          rvectorBlock%RvectorBlock(4),rdefectBlock%RvectorBlock(6),-1.0_DP,1.0_DP,&
          .true.)
      call lsyssc_scalarMatVec (rlevel%rmatrixB2,&
          rvectorBlock%RvectorBlock(5),rdefectBlock%RvectorBlock(6),-1.0_DP,1.0_DP,&
          .true.)
      
      if (rparams%bdualcoupledtoprimal) then
        call lsyssc_scalarMatVec (rlevel%rmatrixMass,&
            rvectorBlock%RvectorBlock(1),rdefectBlock%RvectorBlock(4),1.0_DP,1.0_DP)
        call lsyssc_scalarMatVec (rlevel%rmatrixMass,&
            rvectorBlock%RvectorBlock(2),rdefectBlock%RvectorBlock(5),1.0_DP,1.0_DP)
      end if
      
      if (rparams%bcontrolactive) then
        ! No constraints: -P(u_n) = 1/alpha p_n -> multiply p_n by weighted Mass matrix.
        ! Constraints active: Create the projection and multiply by the weighted mass matrix.
        call lsyssc_copyVector (rvectorBlock%RvectorBlock(4),rvectorTmp)
        call lsyssc_scaleVector (rvectorTmp,-1.0_DP/rparams%dalpha)
        if (rparams%bboundsActive) then
          call projectControlTimestep (rvectorTmp,rparams%dmin1,rparams%dmax1)
        end if
        call lsyssc_scalarMatVec (rlevel%rmatrixMass,&
            rvectorTmp,rdefectBlock%RvectorBlock(1),1.0_DP,1.0_DP)

        call lsyssc_copyVector (rvectorBlock%RvectorBlock(5),rvectorTmp)
        call lsyssc_scaleVector (rvectorTmp,-1.0_DP/rparams%dalpha)
        if (rparams%bboundsActive) then
          call projectControlTimestep (rvectorTmp,rparams%dmin2,rparams%dmax2)
        end if
        call lsyssc_scalarMatVec (rlevel%rmatrixMass,&
            rvectorTmp,rdefectBlock%RvectorBlock(2),1.0_DP,1.0_DP)
        call lsyssc_releaseVector (rvectorTmp)
            
      end if

      ! Now the actual nonlinearity
      if (rparams%bnavierStokes) then

        ! Create a dummy matrix that specifies the matrix structure.
        ! Needed for the assembly.
        call lsysbl_createEmptyMatrix (rtempmatrix,2,2)
        call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rtempmatrix%RmatrixBlock(1,1),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)
        call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rtempmatrix%RmatrixBlock(1,2),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)
        call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rtempmatrix%RmatrixBlock(2,1),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)
        call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rtempmatrix%RmatrixBlock(2,2),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_IGNORE)
        call lsysbl_updateMatStrucInfo(rtempmatrix)
      
        ! Primal equation.
        ! Assemble b = b - y grad(y)
        call lsysbl_deriveSubvector(rvectorBlock,rtempVectorSol, 1,2,.true.)
        call lsysbl_deriveSubvector(rdefectBlock,rtempVectorB, 1,2,.true.)
        call conv_streamlineDiffusionBlk2d ( &
                            rtempVectorSol, rtempVectorSol, 1.0_DP, 0.0_DP,&
                            rstreamlineDiffPrimalDef, CONV_MODDEFECT, &
                            rtempmatrix, rtempVectorSol, rtempVectorB)        

        ! Dual equation.
        ! Assemble b = b + y grad lambda - (grad lambda)^t y
        call lsysbl_deriveSubvector(rvectorBlock,rtempVectorSol, 1,2,.true.)
        call lsysbl_deriveSubvector(rdefectBlock,rtempVectorB, 4,5,.true.)
        call lsysbl_deriveSubvector(rvectorBlock,rtempVectorX, 4,5,.true.)
        call conv_streamlineDiffusionBlk2d ( &
                            rtempVectorSol, rtempVectorSol, 1.0_DP, 0.0_DP,&
                            rstreamlineDiffDualDef, CONV_MODDEFECT, &
                            rtempmatrix, rtempVectorX, rtempVectorB)        
      
        call lsysbl_releaseVector (rtempVectorSol)
        call lsysbl_releaseVector (rtempVectorB)
        call lsysbl_releaseVector (rtempVectorX)
        call lsysbl_releaseMatrix (rtempmatrix)
      
      end if
      
      !call vecio_writeBlockVectorHR (rdefectblock, 'vector', .false.,&
      !                               0, 'vector1.txt', '(E12.5)')
      
    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getRHS (rlevel, rtargetFlow, rrhsBlock)

!<description>
  ! Calculates the RHS vector.
  ! Note: BC are not implemented.
!</desctiption>
  
!<input>
  ! Target flow structure.
  type(t_targetFlow), intent(inout), target :: rtargetFlow

  ! Information about the discretisation
  type(t_level), intent(in) :: rlevel
!</input>

!<output>
  ! RHS vector, to be created
  type(t_vectorBlock), intent(inout) :: rrhsBlock
!</output>

!</subroutine>
  
    ! A bilinear and linear form describing the analytic problem to solve
    type(t_linearForm) :: rlinform
    type(t_collection) :: rcollection
    
    ! Create the RHS vector f of the primal equation
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    call linf_buildVectorScalar (rlevel%rdiscretisation%RspatialDiscr(1),&
                                  rlinform,.true.,rrhsBlock%RvectorBlock(1),&
                                  coeff_RHS_primal_2D)
    call linf_buildVectorScalar (rlevel%rdiscretisation%RspatialDiscr(1),&
                                  rlinform,.true.,rrhsBlock%RvectorBlock(2),&
                                  coeff_RHS_primal_2D)
    call lsyssc_clearVector (rrhsBlock%RvectorBlock(3))
                                  
    ! Create the RHS vector -z of the dual equation.
    ! For that purpose, transfer the target flow parameters via the collection
    ! to the callback routine.
    rcollection%IquickAccess(1) = rtargetFlow%itype
    rcollection%p_rvectorQuickAccess1 => rtargetFlow%rvector
    call linf_buildVectorScalar (rlevel%rdiscretisation%RspatialDiscr(1),&
                                  rlinform,.true.,rrhsBlock%RvectorBlock(4),&
                                  coeff_RHSx_dual_2D,rcollection)
    call linf_buildVectorScalar (rlevel%rdiscretisation%RspatialDiscr(1),&
                                  rlinform,.true.,rrhsBlock%RvectorBlock(5),&
                                  coeff_RHSy_dual_2D,rcollection)
    call lsyssc_clearVector (rrhsBlock%RvectorBlock(6))
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getSystemMatrix (rmatrixBlock,rvectorBlock,rparams,rlevel,bsimple,bdirectSolver)

!<description>
  ! Calculate the system matrix
!</description>
  
!<inputoutput>
  ! The system matrix to create
  type(t_matrixBlock), intent(inout) :: rmatrixBlock
!</inputoutput>

!<input>
  ! Parameter that configure the matrix
  type(t_matrixConfig), intent(in) :: rparams

  ! Information for that level
  type(t_level), intent(in) :: rlevel

  ! Current solution vector
  type(t_vectorBlock), intent(in), target :: rvectorBlock

  ! Simple matrix, the B^T matrices are virtually transposed.
  logical, intent(in) :: bsimple
  
  ! Create matrix for direct solver. Needs a different handling in the
  ! pure Dirichlet case.
  logical, intent(in) :: bdirectSolver
!</input>

!</subroutine>

    ! local variables
    type(t_collection) :: rcollection
    type(t_convStreamlineDiffusion) :: rstreamlineDiffPrimal,rstreamlineDiffDual
    type(t_convStreamlineDiffusion) :: rstreamlineDiffDualR
    type(t_matrixBlock) :: rtempMatrix
    type(t_vectorBlock) :: rtempVectorSol
    type(t_bilinearForm) :: rbilinearForm
    type(t_optcoperator) :: roptcoperator
    
    integer, parameter :: inewmethod = 2
    
    if ((inewmethod .eq. 2) .and. rparams%bnewton) then
      if (.not. lsysbl_isSubmatrixPresent(rmatrixBlock,1,1)) then
        ! First call, initialise empty submatrices.
        call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rmatrixBlock%RmatrixBlock(1,1),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rmatrixBlock%RmatrixBlock(2,2),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
            
        if (rparams%bnewton) then
          call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rmatrixBlock%RmatrixBlock(1,2),&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rmatrixBlock%RmatrixBlock(2,1),&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        end if
        
        call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rmatrixBlock%RmatrixBlock(4,4),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rmatrixBlock%RmatrixBlock(5,5),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        
        if (rparams%bnavierStokes) then
          call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rmatrixBlock%RmatrixBlock(4,5),&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rmatrixBlock%RmatrixBlock(5,4),&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        end if

        if (rparams%bnavierStokes .and. rparams%bnewton) then
          call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rmatrixBlock%RmatrixBlock(4,2),&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rmatrixBlock%RmatrixBlock(5,1),&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        end if

        !if (rparams%bdualcoupledtoprimal) then
          call lsyssc_duplicateMatrix (rlevel%rmatrixMass,rmatrixBlock%RmatrixBlock(4,1),&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          call lsyssc_duplicateMatrix (rlevel%rmatrixMass,rmatrixBlock%RmatrixBlock(5,2),&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        !end if

        !if (rparams%bcontrolactive) then
          call lsyssc_duplicateMatrix (rlevel%rmatrixMass,rmatrixBlock%RmatrixBlock(1,4),&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          call lsyssc_duplicateMatrix (rlevel%rmatrixMass,rmatrixBlock%RmatrixBlock(2,5),&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        !end if
      
      else
        
        ! Remove the B-matrices, attach them later.
        call lsyssc_releaseMatrix(rmatrixBlock%RmatrixBlock(3,1))
        call lsyssc_releaseMatrix(rmatrixBlock%RmatrixBlock(3,2))
        call lsyssc_releaseMatrix(rmatrixBlock%RmatrixBlock(6,4))
        call lsyssc_releaseMatrix(rmatrixBlock%RmatrixBlock(6,5))
        
      end if
      
      ! Clear the matrix
      call lsysbl_clearMatrix (rmatrixBlock)
    
      ! Initialise the operator structure for what we need.
      roptcoperator%dupsamPrimal = rparams%dupsamPrimal
      roptcoperator%dupsamDual = rparams%dupsamDual
      
      ! Timestep-weights
      if (rparams%bemulateTimestep) then
        roptcoperator%dprimalAlpha = 1.0_DP/rparams%dt
        roptcoperator%ddualAlpha   = 1.0_DP/rparams%dt
      end if

      ! Stokes operator
      roptcoperator%dnu = rparams%dnu
      roptcoperator%dprimalBeta = 1.0_DP
      roptcoperator%ddualBeta   = 1.0_DP
      
      ! Nonlinearity
      if (rparams%bnavierStokes) then
        roptcoperator%dprimalDelta =  1.0_DP
        roptcoperator%ddualDelta   = -1.0_DP
        roptcoperator%ddualNewtonTrans = 1.0_DP
        
        ! Whether or not Newton is active has no influence to the
        ! defect, so the following lines are commented out.
        ! if (rparams%bnewton) then
        !   roptcoperator%dprimalNewton    = 1.0_DP
        !   roptcoperator%ddualRDeltaTrans = 1.0_DP
        !   roptcoperator%ddualRNewton     = -1.0_DP
        ! end if
        
      end if
      
      ! Coupling matrices
      if (rparams%bdualcoupledtoprimal) then
        roptcoperator%ddualRAlpha = -1.0_DP
      end if

      if (rparams%bcontrolactive) then
        roptcoperator%dcontrolWeight = -1.0_DP
        roptcoperator%dcontrolMultiplier = -1.0_DP/rparams%dalpha
      end if
      
      if (rparams%bboundsactive) then
        if (.not. rparams%bexactderiv) then
          roptcoperator%ccontrolProjection = 2
        else
          roptcoperator%ccontrolProjection = 1
        end if
        roptcoperator%dmin1 = rparams%dmin1
        roptcoperator%dmax1 = rparams%dmax1
        roptcoperator%dmin2 = rparams%dmin2
        roptcoperator%dmax2 = rparams%dmax2
      end if
      
      ! Calculate the velocity-dependent part of the system matrix.
      call conv_strdiffOptC2dgetDerMatrix (rmatrixBlock,roptcoperator,1.0_DP,&
          rvectorBlock,rvectorBlock,0.001_DP)
      
      ! Finally, attach the B-matrices.
      call lsyssc_duplicateMatrix (rlevel%rmatrixB1,rmatrixBlock%RmatrixBlock(1,3),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      call lsyssc_duplicateMatrix (rlevel%rmatrixB2,rmatrixBlock%RmatrixBlock(2,3),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)

      if (.not. bsimple) then
        call lsyssc_transposeMatrix (rlevel%rmatrixB1,rmatrixBlock%RmatrixBlock(3,1))
        call lsyssc_transposeMatrix (rlevel%rmatrixB2,rmatrixBlock%RmatrixBlock(3,2))
      else
        call lsyssc_transposeMatrix (rlevel%rmatrixB1,rmatrixBlock%RmatrixBlock(3,1),LSYSSC_TR_VIRTUAL)
        call lsyssc_transposeMatrix (rlevel%rmatrixB2,rmatrixBlock%RmatrixBlock(3,2),LSYSSC_TR_VIRTUAL)
      end if
      
      call lsyssc_duplicateMatrix (rlevel%rmatrixB1,rmatrixBlock%RmatrixBlock(4,6),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      call lsyssc_duplicateMatrix (rlevel%rmatrixB2,rmatrixBlock%RmatrixBlock(5,6),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)

      if (.not. bsimple) then
        call lsyssc_transposeMatrix (rlevel%rmatrixB1,rmatrixBlock%RmatrixBlock(6,4))
        call lsyssc_transposeMatrix (rlevel%rmatrixB2,rmatrixBlock%RmatrixBlock(6,5))
      else
        call lsyssc_transposeMatrix (rlevel%rmatrixB1,rmatrixBlock%RmatrixBlock(6,4),LSYSSC_TR_VIRTUAL)
        call lsyssc_transposeMatrix (rlevel%rmatrixB2,rmatrixBlock%RmatrixBlock(6,5),LSYSSC_TR_VIRTUAL)
        call lsyssc_duplicateMatrix (rmatrixBlock%RmatrixBlock(1,3),rmatrixBlock%RmatrixBlock(4,6),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        call lsyssc_duplicateMatrix (rmatrixBlock%RmatrixBlock(2,3),rmatrixBlock%RmatrixBlock(5,6),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      end if
      
      if (rparams%bpureDirichlet .and. bdirectSolver) then
        if (bsimple) then
          ! We cannot use virtually transposed matrices here as we have to overwrite
          ! some data...
          call lsyssc_releaseMatrix(rmatrixBlock%RmatrixBlock(3,1))
          call lsyssc_releaseMatrix(rmatrixBlock%RmatrixBlock(3,2))
          call lsyssc_releaseMatrix(rmatrixBlock%RmatrixBlock(6,4))
          call lsyssc_releaseMatrix(rmatrixBlock%RmatrixBlock(6,5))
          call lsyssc_transposeMatrix (rlevel%rmatrixB1,rmatrixBlock%RmatrixBlock(3,1))
          call lsyssc_transposeMatrix (rlevel%rmatrixB2,rmatrixBlock%RmatrixBlock(3,2))
          call lsyssc_transposeMatrix (rlevel%rmatrixB1,rmatrixBlock%RmatrixBlock(6,4))
          call lsyssc_transposeMatrix (rlevel%rmatrixB2,rmatrixBlock%RmatrixBlock(6,5))
        end if

        ! Fixed pressure; create zero diag matrix and fix first DOF.
        call lsyssc_releaseMatrix (rmatrixBlock%RmatrixBlock(3,3))
        call lsyssc_duplicateMatrix (rlevel%rmatrixDiagP,rmatrixBlock%RmatrixBlock(3,3),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        call lsyssc_clearMatrix(rmatrixBlock%RmatrixBlock(3,3))
        call mmod_replaceLinesByUnitBlk (rmatrixBlock,3,(/1/))

        call lsyssc_releaseMatrix (rmatrixBlock%RmatrixBlock(6,6))
        call lsyssc_duplicateMatrix (rlevel%rmatrixDiagP,rmatrixBlock%RmatrixBlock(6,6),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        call lsyssc_clearMatrix(rmatrixBlock%RmatrixBlock(6,6))
        call mmod_replaceLinesByUnitBlk (rmatrixBlock,6,(/1/))
      else
        rmatrixBlock%RmatrixBlock(6,6)%dscaleFactor = 0.0_DP
        rmatrixBlock%RmatrixBlock(3,3)%dscaleFactor = 0.0_DP
      end if

      !call matio_writeBlockMatrixHR (rmatrixBlock, 'matrix',&
      !    .true., 0, 'matrix2.txt', '(E12.5)', 1E-10_DP)
      
    else if (inewmethod .gt. 0) then
    
      if (.not. lsysbl_isSubmatrixPresent(rmatrixBlock,1,1)) then
        ! First call, initialise empty submatrices.
        call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rmatrixBlock%RmatrixBlock(1,1),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rmatrixBlock%RmatrixBlock(2,2),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
            
        if (rparams%bnewton) then
          call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rmatrixBlock%RmatrixBlock(1,2),&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rmatrixBlock%RmatrixBlock(2,1),&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        end if
        
        call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rmatrixBlock%RmatrixBlock(4,4),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rmatrixBlock%RmatrixBlock(5,5),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        
        if (rparams%bnavierStokes) then
          call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rmatrixBlock%RmatrixBlock(4,5),&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rmatrixBlock%RmatrixBlock(5,4),&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        end if

        if (rparams%bnavierStokes .and. rparams%bnewton) then
          call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rmatrixBlock%RmatrixBlock(4,2),&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rmatrixBlock%RmatrixBlock(5,1),&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        end if

        !if (rparams%bdualcoupledtoprimal) then
          call lsyssc_duplicateMatrix (rlevel%rmatrixMass,rmatrixBlock%RmatrixBlock(4,1),&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          call lsyssc_duplicateMatrix (rlevel%rmatrixMass,rmatrixBlock%RmatrixBlock(5,2),&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        !end if

        !if (rparams%bcontrolactive) then
          call lsyssc_duplicateMatrix (rlevel%rmatrixMass,rmatrixBlock%RmatrixBlock(1,4),&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          call lsyssc_duplicateMatrix (rlevel%rmatrixMass,rmatrixBlock%RmatrixBlock(2,5),&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        !end if
      
      else
        
        ! Remove the B-matrices, attach them later.
        call lsyssc_releaseMatrix(rmatrixBlock%RmatrixBlock(3,1))
        call lsyssc_releaseMatrix(rmatrixBlock%RmatrixBlock(3,2))
        call lsyssc_releaseMatrix(rmatrixBlock%RmatrixBlock(6,4))
        call lsyssc_releaseMatrix(rmatrixBlock%RmatrixBlock(6,5))
        
      end if
      
      ! Clear the matrix
      call lsysbl_clearMatrix (rmatrixBlock)
    
      ! Initialise the operator structure for what we need.
      roptcoperator%dupsamPrimal = rparams%dupsamPrimal
      roptcoperator%dupsamDual = rparams%dupsamDual
      
      ! Timestep-weights
      if (rparams%bemulateTimestep) then
        roptcoperator%dprimalAlpha = 1.0_DP/rparams%dt
        roptcoperator%ddualAlpha   = 1.0_DP/rparams%dt
      end if

      ! Stokes operator
      roptcoperator%dnu = rparams%dnu
      roptcoperator%dprimalBeta = 1.0_DP
      roptcoperator%ddualBeta   = 1.0_DP
      
      ! Nonlinearity
      if (rparams%bnavierStokes) then
        roptcoperator%dprimalDelta =  1.0_DP
        roptcoperator%ddualDelta   = -1.0_DP
        roptcoperator%ddualNewtonTrans = 1.0_DP
        
        ! If Newton is active, we need some more terms
        if (rparams%bnewton) then
          roptcoperator%dprimalNewton    = 1.0_DP
          roptcoperator%ddualRDeltaTrans = 1.0_DP
          roptcoperator%ddualRNewton     = -1.0_DP
        end if
        
      end if
      
      ! Coupling matrices
      if (rparams%bdualcoupledtoprimal) then
        roptcoperator%ddualRAlpha = -1.0_DP
      end if

      if (rparams%bcontrolactive) then
        roptcoperator%dcontrolWeight = -1.0_DP
        roptcoperator%dcontrolMultiplier = -1.0_DP/rparams%dalpha
      end if
      
      if (rparams%bboundsactive) then
        if (.not. rparams%bexactderiv) then
          roptcoperator%ccontrolProjection = 1
        else
          roptcoperator%ccontrolProjection = 2
        end if
        roptcoperator%dmin1 = rparams%dmin1
        roptcoperator%dmax1 = rparams%dmax1
        roptcoperator%dmin2 = rparams%dmin2
        roptcoperator%dmax2 = rparams%dmax2
      end if
      
      ! Calculate the velocity-dependent part of the system matrix.
      call conv_strdiffOptC2dgetMatrix (rmatrixBlock,roptcoperator,1.0_DP,rvectorBlock,rvectorBlock)
      
      ! Finally, attach the B-matrices.
      call lsyssc_duplicateMatrix (rlevel%rmatrixB1,rmatrixBlock%RmatrixBlock(1,3),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      call lsyssc_duplicateMatrix (rlevel%rmatrixB2,rmatrixBlock%RmatrixBlock(2,3),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)

      if (.not. bsimple) then
        call lsyssc_transposeMatrix (rlevel%rmatrixB1,rmatrixBlock%RmatrixBlock(3,1))
        call lsyssc_transposeMatrix (rlevel%rmatrixB2,rmatrixBlock%RmatrixBlock(3,2))
      else
        call lsyssc_transposeMatrix (rlevel%rmatrixB1,rmatrixBlock%RmatrixBlock(3,1),LSYSSC_TR_VIRTUAL)
        call lsyssc_transposeMatrix (rlevel%rmatrixB2,rmatrixBlock%RmatrixBlock(3,2),LSYSSC_TR_VIRTUAL)
      end if
      
      call lsyssc_duplicateMatrix (rlevel%rmatrixB1,rmatrixBlock%RmatrixBlock(4,6),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      call lsyssc_duplicateMatrix (rlevel%rmatrixB2,rmatrixBlock%RmatrixBlock(5,6),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)

      if (.not. bsimple) then
        call lsyssc_transposeMatrix (rlevel%rmatrixB1,rmatrixBlock%RmatrixBlock(6,4))
        call lsyssc_transposeMatrix (rlevel%rmatrixB2,rmatrixBlock%RmatrixBlock(6,5))
      else
        call lsyssc_transposeMatrix (rlevel%rmatrixB1,rmatrixBlock%RmatrixBlock(6,4),LSYSSC_TR_VIRTUAL)
        call lsyssc_transposeMatrix (rlevel%rmatrixB2,rmatrixBlock%RmatrixBlock(6,5),LSYSSC_TR_VIRTUAL)
        call lsyssc_duplicateMatrix (rmatrixBlock%RmatrixBlock(1,3),rmatrixBlock%RmatrixBlock(4,6),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        call lsyssc_duplicateMatrix (rmatrixBlock%RmatrixBlock(2,3),rmatrixBlock%RmatrixBlock(5,6),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      end if
      
      if (rparams%bpureDirichlet .and. bdirectSolver) then
        if (bsimple) then
          ! We cannot use virtually transposed matrices here as we have to overwrite
          ! some data...
          call lsyssc_releaseMatrix(rmatrixBlock%RmatrixBlock(3,1))
          call lsyssc_releaseMatrix(rmatrixBlock%RmatrixBlock(3,2))
          call lsyssc_releaseMatrix(rmatrixBlock%RmatrixBlock(6,4))
          call lsyssc_releaseMatrix(rmatrixBlock%RmatrixBlock(6,5))
          call lsyssc_transposeMatrix (rlevel%rmatrixB1,rmatrixBlock%RmatrixBlock(3,1))
          call lsyssc_transposeMatrix (rlevel%rmatrixB2,rmatrixBlock%RmatrixBlock(3,2))
          call lsyssc_transposeMatrix (rlevel%rmatrixB1,rmatrixBlock%RmatrixBlock(6,4))
          call lsyssc_transposeMatrix (rlevel%rmatrixB2,rmatrixBlock%RmatrixBlock(6,5))
        end if

        ! Fixed pressure; create zero diag matrix and fix first DOF.
        call lsyssc_releaseMatrix (rmatrixBlock%RmatrixBlock(3,3))
        call lsyssc_duplicateMatrix (rlevel%rmatrixDiagP,rmatrixBlock%RmatrixBlock(3,3),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        call lsyssc_clearMatrix(rmatrixBlock%RmatrixBlock(3,3))
        call mmod_replaceLinesByUnitBlk (rmatrixBlock,3,(/1/))

        call lsyssc_releaseMatrix (rmatrixBlock%RmatrixBlock(6,6))
        call lsyssc_duplicateMatrix (rlevel%rmatrixDiagP,rmatrixBlock%RmatrixBlock(6,6),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        call lsyssc_clearMatrix(rmatrixBlock%RmatrixBlock(6,6))
        call mmod_replaceLinesByUnitBlk (rmatrixBlock,6,(/1/))
      else
        rmatrixBlock%RmatrixBlock(6,6)%dscaleFactor = 0.0_DP
        rmatrixBlock%RmatrixBlock(3,3)%dscaleFactor = 0.0_DP
      end if

      !call matio_writeBlockMatrixHR (rmatrixBlock, 'matrix',&
      !    .true., 0, 'matrix.txt', '(E12.5)', 1E-10_DP)
    
    else

      ! Primal preconditioner
      rstreamlineDiffPrimal%dupsam = rparams%dupsamPrimal
      rstreamlineDiffPrimal%dnu = rparams%dnu
      rstreamlineDiffPrimal%ddelta = 1.0_DP
      rstreamlineDiffPrimal%dnewton = 0.0_DP
      if (rparams%bnewton) rstreamlineDiffPrimal%dnewton = 1.0_DP

      ! Dual preconditioner, velocity block
      rstreamlineDiffDual%dupsam = rparams%dupsamDual
      rstreamlineDiffDual%dnu = rparams%dnu
      rstreamlineDiffDual%ddelta = -1.0_DP
      !if (rparams%bnewton) 
      rstreamlineDiffDual%dnewtonTransposed = 1.0_DP

      ! Dual preconditioner, reactive mass matrix block
      rstreamlineDiffDualR%dupsam = rparams%dupsamDual
      rstreamlineDiffDualR%dnu = rparams%dnu
      rstreamlineDiffDualR%ddelta = 0.0_DP
      if (rparams%bnewton) then
        rstreamlineDiffDualR%ddeltaTransposed = 1.0_DP
        rstreamlineDiffDualR%dnewton = -1.0_DP
      end if
      
      ! *************************************************************
      ! Prepare the preconditioner matrix (Newton).
      call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rmatrixBlock%RmatrixBlock(1,1),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rmatrixBlock%RmatrixBlock(2,2),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      call lsyssc_scaleMatrix (rmatrixBlock%RmatrixBlock(1,1),rparams%dnu)
      call lsyssc_scaleMatrix (rmatrixBlock%RmatrixBlock(2,2),rparams%dnu)
      
      ! Newton implies A12/A21
      if (rparams%bnewton) then
        call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rmatrixBlock%RmatrixBlock(1,2),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rmatrixBlock%RmatrixBlock(2,1),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        call lsyssc_clearMatrix (rmatrixBlock%RmatrixBlock(1,2))
        call lsyssc_clearMatrix (rmatrixBlock%RmatrixBlock(2,1))
      end if

      if (rparams%bemulateTimestep) then
        ! One timestep from zero solution to t=1. Emulated by adding a mass matrix
        ! to the Laplace.
        call lsyssc_matrixLinearComb (rlevel%rmatrixMass,1.0_DP/rparams%dt,&
            rmatrixBlock%RmatrixBlock(1,1),1.0_DP,rmatrixBlock%RmatrixBlock(1,1),&
            .false.,.false.,.true.,.true.)
        call lsyssc_matrixLinearComb (rlevel%rmatrixMass,1.0_DP/rparams%dt,&
            rmatrixBlock%RmatrixBlock(2,2),1.0_DP,rmatrixBlock%RmatrixBlock(2,2),&
            .false.,.false.,.true.,.true.)
      end if      

      call lsyssc_duplicateMatrix (rlevel%rmatrixB1,rmatrixBlock%RmatrixBlock(1,3),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      call lsyssc_duplicateMatrix (rlevel%rmatrixB2,rmatrixBlock%RmatrixBlock(2,3),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)

      call lsyssc_releaseMatrix(rmatrixBlock%RmatrixBlock(3,1))
      call lsyssc_releaseMatrix(rmatrixBlock%RmatrixBlock(3,2))
      if (.not. bsimple) then
        call lsyssc_transposeMatrix (rlevel%rmatrixB1,rmatrixBlock%RmatrixBlock(3,1))
        call lsyssc_transposeMatrix (rlevel%rmatrixB2,rmatrixBlock%RmatrixBlock(3,2))
      else
        call lsyssc_transposeMatrix (rlevel%rmatrixB1,rmatrixBlock%RmatrixBlock(3,1),LSYSSC_TR_VIRTUAL)
        call lsyssc_transposeMatrix (rlevel%rmatrixB2,rmatrixBlock%RmatrixBlock(3,2),LSYSSC_TR_VIRTUAL)
      end if

      if (rparams%bpureDirichlet .and. bdirectSolver) then
        if (bsimple) then
          ! Virtually transposed matrices don't allow to replace rows by zero rows...
          call lsyssc_releaseMatrix (rmatrixBlock%RmatrixBlock(3,1))
          call lsyssc_releaseMatrix (rmatrixBlock%RmatrixBlock(3,2))
          call lsyssc_transposeMatrix (rlevel%rmatrixB1,rmatrixBlock%RmatrixBlock(3,1))
          call lsyssc_transposeMatrix (rlevel%rmatrixB2,rmatrixBlock%RmatrixBlock(3,2))
        end if
        ! Fixed pressure; create zero diag matrix and fix first DOF.
        call lsyssc_releaseMatrix (rmatrixBlock%RmatrixBlock(3,3))
        call lsyssc_duplicateMatrix (rlevel%rmatrixDiagP,rmatrixBlock%RmatrixBlock(3,3),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        call lsyssc_clearMatrix(rmatrixBlock%RmatrixBlock(3,3))
        call mmod_replaceLinesByUnitBlk (rmatrixBlock,3,(/1/))
      else
        rmatrixBlock%RmatrixBlock(3,3)%dscaleFactor = 0.0_DP
      end if

      ! ---
      ! Dual equation
      call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rmatrixBlock%RmatrixBlock(4,4),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rmatrixBlock%RmatrixBlock(5,5),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      call lsyssc_scaleMatrix (rmatrixBlock%RmatrixBlock(4,4),rparams%dnu)
      call lsyssc_scaleMatrix (rmatrixBlock%RmatrixBlock(5,5),rparams%dnu)

      if (rparams%bemulateTimestep) then
        ! One timestep from zero solution to t=1. Emulated by adding a mass matrix
        ! to the Laplace.
        call lsyssc_matrixLinearComb (rlevel%rmatrixMass,1.0_DP/rparams%dt,&
            rmatrixBlock%RmatrixBlock(4,4),1.0_DP,rmatrixBlock%RmatrixBlock(4,4),&
            .false.,.false.,.true.,.true.)
        call lsyssc_matrixLinearComb (rlevel%rmatrixMass,1.0_DP/rparams%dt,&
            rmatrixBlock%RmatrixBlock(5,5),1.0_DP,rmatrixBlock%RmatrixBlock(5,5),&
            .false.,.false.,.true.,.true.)
      end if

      call lsyssc_duplicateMatrix (rlevel%rmatrixB1,rmatrixBlock%RmatrixBlock(4,6),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
      call lsyssc_duplicateMatrix (rlevel%rmatrixB2,rmatrixBlock%RmatrixBlock(5,6),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)

      call lsyssc_releaseMatrix(rmatrixBlock%RmatrixBlock(6,4))
      call lsyssc_releaseMatrix(rmatrixBlock%RmatrixBlock(6,5))
      if (.not. bsimple) then
        call lsyssc_transposeMatrix (rlevel%rmatrixB1,rmatrixBlock%RmatrixBlock(6,4))
        call lsyssc_transposeMatrix (rlevel%rmatrixB2,rmatrixBlock%RmatrixBlock(6,5))
      else
        call lsyssc_transposeMatrix (rlevel%rmatrixB1,rmatrixBlock%RmatrixBlock(6,4),LSYSSC_TR_VIRTUAL)
        call lsyssc_transposeMatrix (rlevel%rmatrixB2,rmatrixBlock%RmatrixBlock(6,5),LSYSSC_TR_VIRTUAL)
        call lsyssc_duplicateMatrix (rmatrixBlock%RmatrixBlock(1,3),rmatrixBlock%RmatrixBlock(4,6),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        call lsyssc_duplicateMatrix (rmatrixBlock%RmatrixBlock(2,3),rmatrixBlock%RmatrixBlock(5,6),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      end if

      if (rparams%bpureDirichlet .and. bdirectSolver) then
        if (bsimple) then
          ! Virtually transposed matrices don't allow to replace rows by zero rows...
          call lsyssc_releaseMatrix (rmatrixBlock%RmatrixBlock(6,4))
          call lsyssc_releaseMatrix (rmatrixBlock%RmatrixBlock(6,5))
          call lsyssc_transposeMatrix (rlevel%rmatrixB1,rmatrixBlock%RmatrixBlock(6,4))
          call lsyssc_transposeMatrix (rlevel%rmatrixB2,rmatrixBlock%RmatrixBlock(6,5))
        end if
        ! Fixed pressure; create zero diag matrix and fix first DOF.
        call lsyssc_releaseMatrix (rmatrixBlock%RmatrixBlock(6,6))
        call lsyssc_duplicateMatrix (rlevel%rmatrixDiagP,rmatrixBlock%RmatrixBlock(6,6),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        call lsyssc_clearMatrix(rmatrixBlock%RmatrixBlock(6,6))
        call mmod_replaceLinesByUnitBlk (rmatrixBlock,6,(/1/))
      else
        rmatrixBlock%RmatrixBlock(6,6)%dscaleFactor = 0.0_DP
      end if

      if (rparams%bdualcoupledtoprimal) then
        call lsyssc_duplicateMatrix (rlevel%rmatrixMass,rmatrixBlock%RmatrixBlock(4,1),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
        call lsyssc_duplicateMatrix (rlevel%rmatrixMass,rmatrixBlock%RmatrixBlock(5,2),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
        call lsyssc_scaleMatrix (rmatrixBlock%RmatrixBlock(4,1),-1.0_DP)
        call lsyssc_scaleMatrix (rmatrixBlock%RmatrixBlock(5,2),-1.0_DP)
      else
        call lsyssc_clearMatrix (rmatrixBlock%RmatrixBlock(4,1))
        call lsyssc_clearMatrix (rmatrixBlock%RmatrixBlock(5,2))
      end if

      ! Navier-Stokes implies A45/A54
      if (rparams%bnavierStokes) then
        call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rmatrixBlock%RmatrixBlock(4,5),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rmatrixBlock%RmatrixBlock(5,4),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        call lsyssc_clearMatrix (rmatrixBlock%RmatrixBlock(4,5))
        call lsyssc_clearMatrix (rmatrixBlock%RmatrixBlock(5,4))
      end if

      ! Navier-Stokes+Newton implies A42/A51
      if (rparams%bnavierStokes .and. rparams%bnewton) then
        call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rmatrixBlock%RmatrixBlock(4,2),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        call lsyssc_duplicateMatrix (rlevel%rmatrixLaplace,rmatrixBlock%RmatrixBlock(5,1),&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
        call lsyssc_clearMatrix (rmatrixBlock%RmatrixBlock(4,2))
        call lsyssc_clearMatrix (rmatrixBlock%RmatrixBlock(5,1))
      end if

      if (rparams%bcontrolactive) then
        if (.not. rparams%bexactderiv) then
          ! Copy also the mass matrix for the 4th block. In the case of no
          ! control constraints, that's it!
          call lsyssc_duplicateMatrix (rlevel%rmatrixMass,rmatrixBlock%RmatrixBlock(1,4),&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
          call lsyssc_duplicateMatrix (rlevel%rmatrixMass,rmatrixBlock%RmatrixBlock(2,5),&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_COPY)
          rmatrixBlock%RmatrixBlock(1,4)%dscaleFactor = -(-1.0_DP/rparams%dalpha)
          rmatrixBlock%RmatrixBlock(2,5)%dscaleFactor = -(-1.0_DP/rparams%dalpha)

          if (rparams%bboundsActive .and. rparams%bnewton) then
            ! Filter the mass matrix. Set those rows to zero where the DOF's are out
            ! of bounds. The result is the derivative of the projection operator...
            call massmatfilter (rmatrixBlock%RmatrixBlock(1,4), &
                rvectorBlock%RvectorBlock(4), rparams%dalpha, rparams%dmin1, rparams%dmax1)
            call massmatfilter (rmatrixBlock%RmatrixBlock(2,5), &
                rvectorBlock%RvectorBlock(5), rparams%dalpha, rparams%dmin2, rparams%dmax2)
          end if
        
        else

          call lsyssc_duplicateMatrix (rlevel%rmatrixMass,rmatrixBlock%RmatrixBlock(1,4),&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
          call lsyssc_duplicateMatrix (rlevel%rmatrixMass,rmatrixBlock%RmatrixBlock(2,5),&
              LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
    
          rbilinearForm%itermCount = 1
          rbilinearForm%Idescriptors(1,1) = DER_FUNC
          rbilinearForm%Idescriptors(2,1) = DER_FUNC
          rbilinearForm%ballCoeffConstant = .false.
          rbilinearForm%BconstantCoeff(1) = .false.
          rcollection%p_rvectorQuickAccess1 => rvectorBlock
          rcollection%DquickAccess(3) = rparams%dalpha
          rcollection%IquickAccess(1) = 1
          if (rparams%bboundsActive .and. rparams%bnewton) then
            rcollection%DquickAccess(1) = rparams%dmin1
            rcollection%DquickAccess(2) = rparams%dmax1
          else
            rcollection%DquickAccess(1) = -SYS_MAXREAL_DP
            rcollection%DquickAccess(2) = SYS_MAXREAL_DP
          end if
          call bilf_buildMatrixScalar (rbilinearForm,.true.,rmatrixBlock%RmatrixBlock(1,4),&
                                      coeff_ProjMass,rcollection)

          rcollection%IquickAccess(1) = 2
          if (rparams%bboundsActive .and. rparams%bnewton) then
            rcollection%DquickAccess(1) = rparams%dmin2
            rcollection%DquickAccess(2) = rparams%dmax2
          else
            rcollection%DquickAccess(1) = -SYS_MAXREAL_DP
            rcollection%DquickAccess(2) = SYS_MAXREAL_DP
          end if
          call bilf_buildMatrixScalar (rbilinearForm,.true.,rmatrixBlock%RmatrixBlock(2,5),&
                                      coeff_ProjMass,rcollection)
                                      
        end if
      else
        call lsyssc_clearMatrix (rmatrixBlock%RmatrixBlock(1,4))
        call lsyssc_clearMatrix (rmatrixBlock%RmatrixBlock(2,5))
      end if
    
      ! Now the actual nonlinearity
      if (rparams%bnavierStokes) then

        ! Primal equation.
        ! Assemble y grad(.) ( or   y grad(.) + (.) grad (y)   in case of Newton)
        call lsysbl_deriveSubmatrix (rmatrixBlock,rtempmatrix,&
                                      LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE,1,2)
        call lsysbl_deriveSubvector(rvectorBlock,rtempVectorSol, 1,2,.true.)
        call conv_streamlineDiffusionBlk2d ( &
                            rtempVectorSol, rtempVectorSol, 1.0_DP, 0.0_DP,&
                            rstreamlineDiffPrimal, CONV_MODMATRIX,rtempmatrix)

        ! Dual equation.
        ! Assemble - y grad (.) + (grad .)^t y
        call lsysbl_deriveSubmatrix (rmatrixBlock,rtempmatrix,&
                                      LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE,4,5)
        call lsysbl_deriveSubvector(rvectorBlock,rtempVectorSol, 1,2,.true.)
        call conv_streamlineDiffusionBlk2d ( &
                            rtempVectorSol, rtempVectorSol, 1.0_DP, 0.0_DP,&
                            rstreamlineDiffDual, CONV_MODMATRIX,rtempmatrix)  

        ! Newton?
        if (rparams%bnewton .and. rparams%bdualcoupledtoprimal) then

          ! Assemble - lambda grad(.) + grad(.)^t lambda
          call lsysbl_deriveSubmatrix (rmatrixBlock,rtempmatrix,&
                                      LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE,4,5,1,2)
          call lsysbl_deriveSubvector(rvectorBlock,rtempVectorSol, 4,5,.true.)
          call conv_streamlineDiffusionBlk2d ( &
                            rtempVectorSol, rtempVectorSol, 1.0_DP, 0.0_DP,&
                            rstreamlineDiffDualR, CONV_MODMATRIX,rtempmatrix)  
                            
        end if
        
        call lsysbl_releaseVector (rtempVectorSol)
        call lsysbl_releaseMatrix (rtempmatrix)
      
      end if    
      
      !call matio_writeBlockMatrixHR (rmatrixBlock, 'matrix',&
      !    .true., 0, 'matrix1.txt', '(E12.5)', 1E-10_DP)
    
    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine initTargetFlow (rtargetFlow,itype,rboundary,smesh,ilevel,sfilename)
  
!<description>
  ! Initialises a target flow.
!</description>
  
!<output>
  ! Target flow structure.
  type(t_targetFlow), intent(out) :: rtargetFlow
!</output>
  
!<input>
  ! Type of the target flow.
  ! =0: Poiseuille-flow
  ! =1: Precalculated solution in file sfilename.
  ! =2: Zero target flow
  integer, intent(in) :: itype

  ! Underlying domain.
  type(t_boundary), intent(in), target :: rboundary
  
  ! Filename of the mesh
  character(len=*), intent(in) :: smesh
  
  ! Level of the solution.
  integer, intent(in) :: ilevel
  
  ! Filename of the solution file.
  character(len=*), intent(in) :: sfilename
!</input>

!</subroutine>

    character(len=SYS_STRLEN) :: sarray

    rtargetFlow%itype = itype
    
    ! Read from file!
    if (itype .eq. 1) then
      
      ! Get the mesh
      call tria_readTriFile2D (rtargetFlow%rtriangulation, smesh, rboundary)
       
      ! Refine it.
      call tria_quickRefine2LevelOrdering (ilevel-1,rtargetFlow%rtriangulation,rboundary)
      
      ! And create information about adjacencies and everything one needs from
      ! a triangulation.
      call tria_initStandardMeshFromRaw (rtargetFlow%rtriangulation,rboundary)
      
      ! Set up a block discretisation structure for two components:
      ! Primal and dual solution vector.
      call spdiscr_initBlockDiscr (rtargetFlow%rdiscretisation,3,&
                                     rtargetFlow%rtriangulation, rboundary)

      call spdiscr_initDiscr_simple (rtargetFlow%rdiscretisation%RspatialDiscr(1), &
          EL_EM30,CUB_G3X3,rtargetFlow%rtriangulation, rboundary)
      call spdiscr_duplicateDiscrSc (rtargetFlow%rdiscretisation%RspatialDiscr(1), &
          rtargetFlow%rdiscretisation%RspatialDiscr(2), .true.)
      call spdiscr_deriveSimpleDiscrSc (rtargetFlow%rdiscretisation%RspatialDiscr(1), &
          EL_Q0,CUB_G2X2,rtargetFlow%rdiscretisation%RspatialDiscr(3))
      
      ! Read the solution vector.
      call lsysbl_createVecBlockByDiscr(rtargetFlow%rdiscretisation,&
          rtargetFlow%rvector)
      call vecio_readBlockVectorHR (rtargetFlow%rvector, sarray, .true.,&
                                    0, sfilename, .true.)
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine doneTargetFlow (rtargetFlow)
  
!<description>
  ! Cleans up a target flow.
!</description>
  
!<inputoutput>
  ! Target flow structure.
  type(t_targetFlow), intent(inout) :: rtargetFlow
!</inputoutput>
  
!</subroutine>

    if (rtargetFlow%itype .eq. 1) then
      call lsysbl_releaseVector (rtargetFlow%rvector)
      call spdiscr_releaseBlockDiscr (rtargetFlow%rdiscretisation)
      call tria_done (rtargetFlow%rtriangulation)
    end if
    
    rtargetFlow%itype = 0

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine getBC (rdiscretisation,ibcType,rdiscreteBC,bprimaldual,bpureDirichlet)
  
!<description>
  ! Assemble the boundary conditions for a given discretisation.
!</description>
  
!<input>
  ! Underlying discretisation
  type(t_blockDiscretisation), intent(in) :: rdiscretisation
  
  ! Type of boundary conditions.
  ! =0: Channel without flow.
  ! =1: Channel with inflow.
  ! =2: Driven Cavity
  integer, intent(in) :: ibcType
  
  ! TRUE: Both, primal and dual. FALSE: Only primal
  logical, intent(in) :: bprimaldual
!</input>

!<output>
  ! Receives the discrete boundary conditions
  type(t_discreteBC), intent(out) :: rdiscreteBC

  ! Set to TRUE if there are no Neumann boundary components in the
  ! flow configuration of the target flow.
  logical, intent(out), optional :: bpureDirichlet
!</output>

!</subroutine>

    ! local variables
  
    ! A set of variables describing the discrete boundary conditions.    
    type(t_boundaryRegion) :: rboundaryRegion
    logical :: binflowBC

    binflowBC = ibcType .eq. 1
    if (present(bpureDirichlet)) then
      bpureDirichlet = ibcType .eq. 2
    end if

    call bcasm_initDiscreteBC(rdiscreteBC)

    if ((ibcType .eq. 0) .or. (ibcType .eq. 1)) then
      ! -----
      ! Boundary conditions
      !
      ! In our example, we have pure Dirichlet-0-BC on all boundary components.
      
      ! Edge 1 of boundary component 1 the domain.
      call boundary_createRegion(rdiscretisation%p_rboundary,1,1,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesZero_2D)
                               
      ! Now to the edge 2 of boundary component 1 the domain.
      call boundary_createRegion(rdiscretisation%p_rboundary,1,3,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesZero_2D)

      if (binflowBC) then
        ! Now to the edge 2 of boundary component 1 the domain.
        call boundary_createRegion(rdiscretisation%p_rboundary,1,4,rboundaryRegion)
        rboundaryRegion%iproperties = 0
        call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                          rboundaryRegion,rdiscreteBC,&
                                          getBoundaryValuesX_2D)
      end if
                               
      ! Edge 3 of boundary component 1.
      call boundary_createRegion(rdiscretisation%p_rboundary,1,1,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesZero_2D)
      
      ! Edge 4 of boundary component 1. That's it.
      call boundary_createRegion(rdiscretisation%p_rboundary,1,3,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesZero_2D)

      if (binflowBC) then
        ! Edge 4 of boundary component 1. That's it.
        call boundary_createRegion(rdiscretisation%p_rboundary,1,4,rboundaryRegion)
        call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                          rboundaryRegion,rdiscreteBC,&
                                          getBoundaryValuesY_2D)
      end if

      if (.not. bprimaldual) return

      ! ---
      ! The same for the dual variable.
      call boundary_createRegion(rdiscretisation%p_rboundary,1,1,rboundaryRegion)
      rboundaryRegion%iproperties = 0
      call bcasm_newDirichletBConRealBD (rdiscretisation,4,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesZero_2D)
                               
      ! Now to the edge 2 of boundary component 1 the domain.
      call boundary_createRegion(rdiscretisation%p_rboundary,1,3,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,4,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesZero_2D)

      if (binflowBC) then
        ! Now to the edge 2 of boundary component 1 the domain.
        call boundary_createRegion(rdiscretisation%p_rboundary,1,4,rboundaryRegion)
        call bcasm_newDirichletBConRealBD (rdiscretisation,4,&
                                          rboundaryRegion,rdiscreteBC,&
                                          getBoundaryValuesZero_2D)
      end if
                               
      ! Edge 3 of boundary component 1.
      call boundary_createRegion(rdiscretisation%p_rboundary,1,1,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,5,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesZero_2D)
      
      ! Edge 4 of boundary component 1. That's it.
      call boundary_createRegion(rdiscretisation%p_rboundary,1,3,rboundaryRegion)
      call bcasm_newDirichletBConRealBD (rdiscretisation,5,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesZero_2D)

      if (binflowBC) then
        ! Edge 4 of boundary component 1. That's it.
        call boundary_createRegion(rdiscretisation%p_rboundary,1,4,rboundaryRegion)
        call bcasm_newDirichletBConRealBD (rdiscretisation,5,&
                                          rboundaryRegion,rdiscreteBC,&
                                          getBoundaryValuesZero_2D)
      end if
    else if (ibcType .eq. 2) then
      ! -----
      ! Boundary conditions
      !
      ! Dirichlet-0-BC on all boundary components except the upper edge
      
      ! Edge 1 of boundary component 1 the domain.
      call boundary_createRegion(rdiscretisation%p_rboundary,1,1,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesZero_2D)

      ! Edge 2 of boundary component 1 the domain.
      call boundary_createRegion(rdiscretisation%p_rboundary,1,2,rboundaryRegion)
      rboundaryRegion%iproperties = 0
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesZero_2D)
                               
      ! Now to the edge 2 of boundary component 1 the domain.
      call boundary_createRegion(rdiscretisation%p_rboundary,1,3,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValues1_2D)

      ! Now to the edge 2 of boundary component 1 the domain.
      call boundary_createRegion(rdiscretisation%p_rboundary,1,4,rboundaryRegion)
      rboundaryRegion%iproperties = 0
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesZero_2D)
                               
      ! Edge 3 of boundary component 1.
      call boundary_createRegion(rdiscretisation%p_rboundary,1,1,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesZero_2D)
      
      ! Edge 3 of boundary component 1.
      call boundary_createRegion(rdiscretisation%p_rboundary,1,2,rboundaryRegion)
      rboundaryRegion%iproperties = 0
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesZero_2D)

      ! Edge 4 of boundary component 1. That's it.
      call boundary_createRegion(rdiscretisation%p_rboundary,1,3,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesZero_2D)

      ! Edge 4 of boundary component 1. That's it.
      call boundary_createRegion(rdiscretisation%p_rboundary,1,4,rboundaryRegion)
      rboundaryRegion%iproperties = 0
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesZero_2D)
      if (.not. bprimaldual) return

      ! ---
      ! The same for the dual variable.
      call boundary_createRegion(rdiscretisation%p_rboundary,1,1,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rdiscretisation,4,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesZero_2D)

      call boundary_createRegion(rdiscretisation%p_rboundary,1,2,rboundaryRegion)
      rboundaryRegion%iproperties = 0
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rdiscretisation,4,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesZero_2D)
                               
      ! Now to the edge 2 of boundary component 1 the domain.
      call boundary_createRegion(rdiscretisation%p_rboundary,1,3,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rdiscretisation,4,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesZero_2D)

      ! Now to the edge 2 of boundary component 1 the domain.
      call boundary_createRegion(rdiscretisation%p_rboundary,1,4,rboundaryRegion)
      rboundaryRegion%iproperties = 0
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rdiscretisation,4,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesZero_2D)

      ! Edge 3 of boundary component 1.
      call boundary_createRegion(rdiscretisation%p_rboundary,1,1,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rdiscretisation,5,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesZero_2D)

      ! Edge 3 of boundary component 1.
      call boundary_createRegion(rdiscretisation%p_rboundary,1,2,rboundaryRegion)
      rboundaryRegion%iproperties = 0
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rdiscretisation,5,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesZero_2D)
      
      ! Edge 4 of boundary component 1. That's it.
      call boundary_createRegion(rdiscretisation%p_rboundary,1,3,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rdiscretisation,5,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesZero_2D)

      ! Edge 4 of boundary component 1. That's it.
      call boundary_createRegion(rdiscretisation%p_rboundary,1,4,rboundaryRegion)
      rboundaryRegion%iproperties = 0
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rdiscretisation,5,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesZero_2D)
    else if (ibcType .eq. 3) then
      ! -----
      ! Boundary conditions
      !
      ! Dirichlet-0-BC the bottom, Dirichlet-1-BC at the top
      
      ! Edge 1 of boundary component 1 the domain.
      call boundary_createRegion(rdiscretisation%p_rboundary,1,1,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesZero_2D)

      ! Now to the edge 2 of boundary component 1 the domain.
      call boundary_createRegion(rdiscretisation%p_rboundary,1,3,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rdiscretisation,1,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValues1_2D)

      ! Edge 3 of boundary component 1.
      call boundary_createRegion(rdiscretisation%p_rboundary,1,1,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesZero_2D)
      
      ! Edge 4 of boundary component 1. That's it.
      call boundary_createRegion(rdiscretisation%p_rboundary,1,3,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rdiscretisation,2,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesZero_2D)

      if (.not. bprimaldual) return

      ! ---
      ! The same for the dual variable.
      call boundary_createRegion(rdiscretisation%p_rboundary,1,1,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rdiscretisation,4,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesZero_2D)

      ! Now to the edge 2 of boundary component 1 the domain.
      call boundary_createRegion(rdiscretisation%p_rboundary,1,3,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rdiscretisation,4,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesZero_2D)

      ! Edge 3 of boundary component 1.
      call boundary_createRegion(rdiscretisation%p_rboundary,1,1,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rdiscretisation,5,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesZero_2D)

      ! Edge 4 of boundary component 1. That's it.
      call boundary_createRegion(rdiscretisation%p_rboundary,1,3,rboundaryRegion)
      rboundaryRegion%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
      call bcasm_newDirichletBConRealBD (rdiscretisation,5,&
                                        rboundaryRegion,rdiscreteBC,&
                                        getBoundaryValuesZero_2D)

    end if

  end subroutine

  ! ***************************************************************************

  subroutine massmatfilter (rmatrix, rvector, dalphaC, dmin, dmax)
      
  ! Filters a mass matrix. The lines in the matrix rmatrix corresponding
  ! to all entries in the (control-)vector violating the constraints
  ! of the problem.
  
  ! Matrix to be filtered
  type(t_matrixScalar), intent(inout) :: rmatrix

  ! Vector containing a dial solution lambda. Whereever -1/alpha*lambda
  ! violates the control constraints given by rmatrixComponents, the corresponding
  ! lines are set to 0.
  type(t_vectorScalar), intent(in) :: rvector
  
  ! ALPHA regularisation parameter from the space-time matrix
  real(dp), intent(in) :: dalphaC
  
  ! minimum bound for the control
  real(dp), intent(in) :: dmin

  ! maximum bound for the control
  real(dp), intent(in) :: dmax
  
    ! local variables
    real(dp), dimension(:), pointer :: p_Ddata
    integer, dimension(:), allocatable :: p_Idofs
    integer :: i,nviolate
    real(dp) :: du
    
    ! Get the vector data
    call lsyssc_getbase_double (rvector,p_Ddata)
    
    ! Figure out the DOF's violating the constraints
    allocate(p_Idofs(rvector%NEQ))
    
    nviolate = 0
    do i=1,rvector%NEQ
      du = -p_Ddata(i)/dalphaC
      if ((du .le. dmin) .or. (du .ge. dmax)) then
        nviolate = nviolate + 1
        p_Idofs(nviolate) = i
      end if
    end do
    
    if (nviolate .gt. 0) then
      ! Filter the matrix
      call mmod_replaceLinesByZero (rmatrix,p_Idofs(1:nviolate))
    end if
    
    deallocate(p_Idofs)

  end subroutine    

  ! ***************************************************************************
  
!<subroutine>

  subroutine projectControlTimestep (rdualSolution,dumin,dumax)

!<description>
  ! Projects a dual solution vector u in such a way, that
  ! dumin <= u <= dumax holds.
!</description>

!<input>
  ! Minimum value for u
  real(DP), intent(in) :: dumin

  ! Maximum value for u
  real(DP), intent(in) :: dumax
!</input>

!<inputoutput>
  ! Vector to be restricted
  type(t_vectorScalar), intent(inout) :: rdualSolution
!</inputoutput>

!</subroutine>
 
    ! local variables
    real(DP), dimension(:), pointer :: p_Ddata
    integer :: i
    
    ! Get the vector array
    call lsyssc_getbase_double (rdualSolution,p_Ddata)
    
    ! Restrict the vector
    do i=1,rdualSolution%NEQ
      p_Ddata(i) = min(max(p_Ddata(i),dumin),dumax)
    end do

  end subroutine 
  
  ! ***************************************************************************

!<subroutine>

  subroutine navstokes2d_0_simple
  
!<description>
  ! This is an all-in-one poisson solver for directly solving a Poisson
  ! problem without making use of special features like collections
  ! and so on. The routine performs the following tasks:
  !
  ! 1.) Read in parametrisation
  ! 2.) Read in triangulation
  ! 3.) Set up RHS
  ! 4.) Set up matrix
  ! 5.) Create solver structure
  ! 6.) Solve the problem
  ! 7.) Write solution to GMV file
  ! 8.) Release all variables, finish
!</description>

!</subroutine>

    ! Definitions of variables.
    !
    ! We need a couple of variables for this problem. Let's see...
    !
    ! An object for saving the domain:
    type(t_boundary) :: rboundary
    
    ! An object specifying the discretisation.
    ! This contains also information about trial/test functions,...
    type(t_blockDiscretisation) :: rdiscrOutBlock
    type(t_spatialDiscretisation) :: rdiscrOutput
    
    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    type(t_matrixBlock), dimension(:), allocatable :: RmatrixBlock
    type(t_vectorBlock), target :: rvectorBlock
    type(t_vectorBlock) :: rrhsBlock,rvecout
    type(t_vectorScalar), target :: rvectorOutputX,rvectorOutputY,rvectorTmp
    type(t_vectorBlock) :: rdefectBlock

    ! A solver node that accepts parameters for the linear solver    
    type(t_linsolNode), pointer :: p_rsolverNode, p_rprecSolver
    type(t_filterChain), dimension(3), target :: RfilterChain
    type(t_filterChain), dimension(:), pointer :: p_RfilterChain
    type(t_linsolMG2LevelInfo), pointer     :: p_rlevelInfo

    ! NLMAX receives the level where we want to solve.
    integer :: NLMIN,NLMAX
    
    ! Error indicator during initialisation of the solver
    integer :: ierror, ilevel, ilinearsolver
    integer :: ibctype
    integer :: itargetFlow, ioutputlevel
    type(t_discreteBC) :: rdiscreteBCout
    
    ! Data for the Newton iteration
    integer :: ite, nmaxiterations, idiscretisation, nmg
    real(dp) :: dinitRes, dcurrentRes
    
    type(t_targetFlow), target :: rtargetFlow
    
    real(DP), dimension(10) :: Derror
    
    type(t_collection) :: rcollection
    
    ! Output block for UCD output to GMV file
    type(t_ucdExport) :: rexport
    real(DP), dimension(:), pointer :: p_Ddata,p_DdataX,p_DdataY
    !real(DP), dimension(:), pointer :: p_Dtmp1,p_Dtmp2,p_Dtmp3,p_Dtmp4

    ! All information for the discretsiation
    type(t_level), dimension(:), allocatable :: Rlevel
    
    ! Information for setting up the system matrix
    type(t_matrixConfig) :: rparams
    
    integer :: icurrentmaxre,icurrentre,icurrentalpha,icurrentupsam
    integer, dimension(4), parameter :: imaxre = (/1000,500,250,100/) 
    integer, dimension(9), parameter :: ire = (/1000,500,250,100,50,25,10,5,1/) 
    real(dp), dimension(3), parameter :: Dalpha = (/0.01_DP,0.1_DP,0.001_DP/)
    real(dp), dimension(6), parameter :: Dupsam = (/0.0_DP,0.2_DP,0.4_DP,0.6_DP,0.8_DP,1.0_DP/)
    !integer, dimension(1), parameter :: ire = (/100/)
    !real(dp), dimension(1), parameter :: Dalpha = (/0.01_DP/)

    ! Ok, let's start. 
    !
    ! We want to solve our Poisson problem on level... 
    NLMAX = 5
    
    ! Minimum level in the MG solver
    NLMIN = 2
    
    ! Newton iteration counter
    nmaxiterations = 20 !200
    
    ! Relaxation parameter
    !rparams%dalpha = 0.1_DP
    
    ! Control active or not.
    rparams%bcontrolactive = .true.
    rparams%bdualcoupledtoprimal = .true.
    
    ! Bounds on the control
    rparams%bboundsActive = .false.
    rparams%dmin1 = -0.05
    rparams%dmax1 = 0.05
    rparams%dmin2 = -0.05
    rparams%dmax2 = 0.05
    
    ! Viscosity constant
    !rparams%dnu = 1.0_DP/1.0_DP !/200.0_DP
    
    ! Discretisation. 1=EM30, 2=Q2
    idiscretisation = 1
    
    ! TRUE emulates a timestep with implicit Euler, timestep length dt, 
    ! starting from a zero solution
    rparams%bemulateTimestep = .false.
    
    ! Timestep length
    rparams%dt = 1.0_DP
    
    ! Stabilisation parameter
    rparams%dupsamPrimal = 0.7_DP
    rparams%dupsamDual = 0.7_DP
    
    ! TRUE: Use Newton iteration
    rparams%bnewton = .true.
    
    ! TRUE: Use Navier-Stokes. FALSE: Use Stokes
    rparams%bnavierStokes = .true.
    
    ! TRUE: Use exact derivative of the semismooth operator
    ! (mass matrix set up with an appropriate coefficient) instead
    ! of zero rows in the mass matrix.
    rparams%bexactderiv = .false.
    
    ! Activate inflow BC's
    ! =0: no inflow
    ! =1: Poiseuille
    ! =2: Driven Cavity
    ! =3: ?
    ibctype = 2
    
    ! Type of the target flow.
    ! =0: Poiseuille
    ! =1: Read from file
    ! =2: Zero
    itargetFlow = 1
    
    ! Type of linear solver.
    ! =0: UMFPACK
    ! =1: MG
    ilinearsolver = 1
    
    ! Output level
    ioutputlevel = 2
    
    ! Allocate level information
    allocate(Rlevel(NLMAX))
    allocate(RmatrixBlock(NLMAX))
    
    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    call boundary_read_prm(rboundary, './pre/QUAD.prm')

    ! Create the target flow.
    do icurrentmaxre = 1,size(imaxre)
      call initTargetFlow (rtargetFlow,itargetFlow,rboundary,&
          './pre/QUAD.tri',NLMAX,'./ns/navstdc'//trim(sys_siL(NLMAX,10))//&
          're'//sys_siL(imaxre(icurrentmaxre),10))
          
      ! Now read in the basic triangulation.
      call tria_readTriFile2D (Rlevel(1)%rtriangulation, './pre/QUAD.tri', rboundary)
      call tria_initStandardMeshFromRaw (Rlevel(1)%rtriangulation,rboundary)
       
      call tria_infoStatistics(Rlevel(1)%rtriangulation,.true.)
      
      do ilevel = 2,NLMAX
        ! Refine it.
        call tria_refine2LevelOrdering (Rlevel(ilevel-1)%rtriangulation,&
            Rlevel(ilevel)%rtriangulation,rboundary)
      
        ! And create information about adjacencies and everything one needs from
        ! a triangulation.
        call tria_initStandardMeshFromRaw (Rlevel(ilevel)%rtriangulation,rboundary)

        ! Print mesh statistics
        call tria_infoStatistics(Rlevel(ilevel)%rtriangulation,.false.)
        
        ! Release the old mesh if it's not used
        if (ilevel-1 .lt. NLMIN) then
          call tria_done (Rlevel(ilevel-1)%rtriangulation)
        end if
      end do
      
      ! Remove redundant mesh data
      do ilevel=NLMAX-1,NLMIN,-1
        call tria_compress2LevelOrdHierarchy (Rlevel(ilevel+1)%rtriangulation,&
            Rlevel(ilevel)%rtriangulation)
      end do
      
      call output_lbrk()
      
      do ilevel = NLMIN,NLMAX
        ! Set up a block discretisation structure for two components:
        ! Primal and dual solution vector.
        call spdiscr_initBlockDiscr (Rlevel(ilevel)%rdiscretisation,6,&
                                      Rlevel(ilevel)%rtriangulation, rboundary)
        
        ! Set up the blocks. Both are discretised with the same finite element.
        if (idiscretisation .eq. 1) then
        
          call spdiscr_initDiscr_simple (Rlevel(ilevel)%rdiscretisation%RspatialDiscr(1), &
              EL_EM30,CUB_G3X3,Rlevel(ilevel)%rtriangulation, rboundary)
          call spdiscr_duplicateDiscrSc (Rlevel(ilevel)%rdiscretisation%RspatialDiscr(1), &
              Rlevel(ilevel)%rdiscretisation%RspatialDiscr(2), .true.)
          call spdiscr_duplicateDiscrSc (Rlevel(ilevel)%rdiscretisation%RspatialDiscr(1), &
              Rlevel(ilevel)%rdiscretisation%RspatialDiscr(4), .true.)
          call spdiscr_duplicateDiscrSc (Rlevel(ilevel)%rdiscretisation%RspatialDiscr(1), &
              Rlevel(ilevel)%rdiscretisation%RspatialDiscr(5), .true.)
              
          call spdiscr_deriveSimpleDiscrSc (Rlevel(ilevel)%rdiscretisation%RspatialDiscr(1), &
              EL_Q0,CUB_G2X2,Rlevel(ilevel)%rdiscretisation%RspatialDiscr(3))
          call spdiscr_duplicateDiscrSc (Rlevel(ilevel)%rdiscretisation%RspatialDiscr(3), &
              Rlevel(ilevel)%rdiscretisation%RspatialDiscr(6), .true.)
              
        else if (idiscretisation .eq. 2) then
        
          call spdiscr_initDiscr_simple (Rlevel(ilevel)%rdiscretisation%RspatialDiscr(1), &
              EL_Q2,CUB_G4X4,Rlevel(ilevel)%rtriangulation, rboundary)
          call spdiscr_duplicateDiscrSc (Rlevel(ilevel)%rdiscretisation%RspatialDiscr(1), &
              Rlevel(ilevel)%rdiscretisation%RspatialDiscr(2), .true.)
          call spdiscr_duplicateDiscrSc (Rlevel(ilevel)%rdiscretisation%RspatialDiscr(1), &
              Rlevel(ilevel)%rdiscretisation%RspatialDiscr(4), .true.)
          call spdiscr_duplicateDiscrSc (Rlevel(ilevel)%rdiscretisation%RspatialDiscr(1), &
              Rlevel(ilevel)%rdiscretisation%RspatialDiscr(5), .true.)
              
          call spdiscr_deriveSimpleDiscrSc (Rlevel(ilevel)%rdiscretisation%RspatialDiscr(1), &
              EL_QP1,CUB_G2X2,Rlevel(ilevel)%rdiscretisation%RspatialDiscr(3))
          call spdiscr_duplicateDiscrSc (Rlevel(ilevel)%rdiscretisation%RspatialDiscr(3), &
              Rlevel(ilevel)%rdiscretisation%RspatialDiscr(6), .true.)
              
        end if
                
        ! Create not-changing matrices
        call calcStandardMatrices (Rlevel(ilevel))
        
        ! Create a temp vector
        call lsysbl_createVecBlockByDiscr (Rlevel(ilevel)%rdiscretisation,Rlevel(ilevel)%rtempVector)

        ! Assemble the boundary conditions
        call getBC (Rlevel(ilevel)%rdiscretisation,ibctype,&
            Rlevel(ilevel)%rdiscreteBC,.true.,rparams%bpureDirichlet)

      end do

      ! -----
      ! Linear system: Basic structure.
      !
      ! Create a 2x2 matrix and 2-block vectors from the discretisation.
      do ilevel=NLMIN,NLMAX
        call lsysbl_createMatBlockByDiscr (Rlevel(ilevel)%rdiscretisation,RmatrixBlock(ilevel))
      end do
      
      call lsysbl_createVecBlockByDiscr (Rlevel(NLMAX)%rdiscretisation,rrhsBlock)
      call lsysbl_createVecBlockByDiscr (Rlevel(NLMAX)%rdiscretisation,rvectorBlock)
      call lsysbl_createVecBlockByDiscr (Rlevel(NLMAX)%rdiscretisation,rdefectBlock)
      call lsyssc_createVecByDiscr (Rlevel(NLMAX)%rdiscretisation%RspatialDiscr(1),rvectorTmp)
      
      call getRHS (Rlevel(NLMAX), rtargetFlow, rrhsBlock)

      ! Implement BC's
      call vecfil_discreteBCrhs (rrhsBlock,Rlevel(NLMAX)%rdiscreteBC)
      
      do icurrentupsam = 1,size(Dupsam)
      
        rparams%dupsamPrimal = Dupsam(icurrentupsam)
        rparams%dupsamDual = 0.0_DP
      
        do icurrentalpha = 1,size(Dalpha)

          rparams%dalpha = Dalpha(icurrentalpha)
        
          do icurrentre = 1,size(Ire)

            ! Don't calculate higher RE as the target flow
            if (Ire(icurrentre) .gt. imaxre(icurrentmaxre)) cycle

            if (ioutputLevel .ge. 1) then
            
              call output_line ('---')
              call output_line ('maxre  = '//sys_siL(Imaxre(icurrentmaxre),10))
              call output_line ('re     = '//sys_siL(Ire(icurrentre),10))
              call output_line ('alpha  = '//sys_sdL(rparams%dalpha,10))
              call output_line ('dupsam = '//sys_sdL(rparams%dupsamPrimal,10))
             
            end if

            ! Current RE number
            rparams%dnu = 1.0_DP/real(Ire(icurrentre))

            ! -----
            ! Nonlinear iteration
            !
            ! Initialise the solution vector.
            call lsysbl_clearVector (rvectorBlock)
            
            call vecfil_discreteBCsol (rvectorBlock,Rlevel(NLMAX)%rdiscreteBC)

            ! Filter chain for boundary conditions in the linear solver
            RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL
            if (rparams%bpureDirichlet) then
              RfilterChain(2)%ifilterType = FILTER_TOL20
              RfilterChain(2)%itoL20component = 3
              RfilterChain(3)%ifilterType = FILTER_TOL20
              RfilterChain(3)%itoL20component = 6
            end if
            p_RfilterChain => RfilterChain
            
            ! Prepare an UMFPACK solver for the lienar systems
            select case(ilinearsolver)
            case (0)
              call linsol_initUMFPACK4(p_rsolverNode)
              
            case (1)
              call linsol_initMultigrid2(p_rsolverNode,NLMAX-NLMIN+1,p_RfilterChain)
              
              ! Add the coarse grid solver
              call linsol_getMultigrid2Level (p_rsolverNode,1,p_rlevelInfo)
              call linsol_initVANKA (p_rprecSolver,0.7_DP,LINSOL_VANKA_GENERAL)
              call linsol_initBiCGStab (p_rlevelInfo%p_rcoarseGridSolver,p_rprecSolver,p_RfilterChain)
              p_rlevelInfo%p_rcoarseGridSolver%nmaxIterations = 1000
              !call linsol_initUMFPACK4(p_rlevelInfo%p_rcoarseGridSolver)
              
              do ilevel = NLMIN+1,NLMAX
                call linsol_getMultigrid2Level (p_rsolverNode,ilevel-NLMIN+1,p_rlevelInfo)
              
                ! Add the smoother
                call linsol_initVANKA (p_rprecSolver,1.0_DP,LINSOL_VANKA_GENERAL)
                call linsol_initBiCGStab (p_rlevelInfo%p_rpresmoother,p_rprecSolver,p_RfilterChain)
                call linsol_convertToSmoother (p_rlevelInfo%p_rpresmoother,4,1.0_DP)
                p_rlevelInfo%p_rpostsmoother => p_rlevelInfo%p_rpresmoother
                
                ! Add the interlevel projection structure
                call mlprj_initProjectionDiscr (Rlevel(ilevel)%rprojection,&
                    Rlevel(ilevel)%rdiscretisation)
                    
                ! Set up the interlevel projection structure on all levels
                call linsol_initProjMultigrid2Level(p_rlevelInfo,&
                    Rlevel(ilevel)%rprojection)

              end do
              
            end select
            !p_rsolverNode%p_rsubnodeUMFPACK4%imatrixDebugOutput = 1
          
            ! Set the output level of the solver to 2 for some output
            p_rsolverNode%ioutputLevel = ioutputlevel-2
            
            !p_rsolverNode%p_rsubnodeMultigrid2%rcoarseGridCorrection%ccorrectionType = CGCOR_SCALARENERGYMIN
            
            ! Count the number of steps in the linear solver
            nmg = 0
            
            do ite = 1,nmaxIterations+1
            
        !      if (ite .eq. 2) then
        !        bcontrolactive = .false.
        !        rmatrixBlock%RmatrixBlock(1,4)%dscaleFactor = 0.0_DP
        !        rmatrixBlock%RmatrixBlock(2,5)%dscaleFactor = 0.0_DP
        !      end if
            
              ! Calculate the system matrix on all levels
              call lsysbl_copyVector (rvectorBlock,Rlevel(NLMAX)%rtempVector)
              do ilevel = NLMAX,NLMIN,-1
                
                ! Project down the solution vector; this is the evaluation point for
                ! the matrix.
                if (ilevel .lt. NLMAX) then
                  ! Abuse the defect vector as temp vector
                  call mlprj_performInterpolation (Rlevel(ilevel+1)%rprojection,&
                      Rlevel(ilevel)%rtempVector, Rlevel(ilevel+1)%rtempVector,rdefectBlock%RvectorBlock(1))
                end if
              
                call getSystemMatrix (RmatrixBlock(ilevel),Rlevel(ilevel)%rtempVector,rparams,Rlevel(ilevel),&
                    .false.,ilinearsolver .eq. 0)
                call matfil_discreteBC (RmatrixBlock(ilevel),Rlevel(ilevel)%rdiscreteBC)
              end do

              ! Calculate the nonlinear defect
              call getDefect (rdefectBlock,rvectorBlock,rrhsBlock,rparams,Rlevel(NLMAX))

              ! Include boundary conditions
              call vecfil_discreteBCdef (rdefectBlock,Rlevel(NLMAX)%rdiscreteBC)

              ! Check for convergence
              if (ite .eq. 1) then
                dinitRes = lsysbl_vectorNorm(rdefectBlock,LINALG_NORML2)
                if (ioutputlevel .ge. 2) then
                  call output_line ('Iteration: '//TRIM(sys_siL(ite-1,10))//&
                      ', Initial defect: '//sys_sdEL(dinitRes,10))
                end if
              else
                dcurrentRes = lsysbl_vectorNorm(rdefectBlock,LINALG_NORML2)
                if (ioutputlevel .ge. 2) then
                  call output_line ('Iteration: '//TRIM(sys_siL(ite-1,10))//&
                      ', Current defect: '//sys_sdEL(dcurrentRes,10))
                end if
                if (dcurrentRes .lt. 1E-13_DP) exit
                if (dcurrentRes .lt. dinitRes * 1E-13_DP) exit
                if (ite .gt. nmaxIterations) exit
              end if
              
              ! Attach the system matrix to the solver.
              ! First create an array with the matrix data (on all levels, but we
              ! only have one level here), then call the initialisation 
              ! routine to attach all these matrices.
              ! Remark: Don't make a call like
              !    CALL linsol_setMatrices(p_RsolverNode,(/p_rmatrix/))
              ! This doesn't work on all compilers, since the compiler would have
              ! to create a temp array on the stack - which does not always work!
              call linsol_setMatrices(p_RsolverNode,RmatrixBlock(NLMIN:NLMAX))
              
              ! Initialise structure/data of the solver. This allows the
              ! solver to allocate memory / perform some precalculation
              ! to the problem.
              if (ite .eq. 1) then
                ! In the first iteration, initialise the linear solver
                call linsol_initStructure (p_rsolverNode, ierror)
                if (ierror .ne. LINSOL_ERR_NOERROR) stop
              end if
              
              call linsol_initData (p_rsolverNode, ierror)
              if (ierror .ne. LINSOL_ERR_NOERROR) stop
              
              ! Preconditioning of the defect by the inverse of the Newton matrix
              call linsol_precondDefect (p_rsolverNode,rdefectBlock)
              if (p_rsolverNode%iresult .ne. 0) then
                call output_line("Iteration canceled, solver broke down.")
                nmg = 0
                exit
              end if
              
              nmg = nmg + p_rsolverNode%iiterations
              
              ! Release solver data and structure
              call linsol_doneData (p_rsolverNode)
              
              ! Sum up the correction to the current solution.
              call lsysbl_vectorLinearComb (rdefectBlock,rvectorBlock,1.0_DP,1.0_DP)
              
              ! Plug in the BC's to the solution
              call vecfil_discreteBCsol (rvectorBlock,Rlevel(NLMAX)%rdiscreteBC)
              
            end do

            ! Clean up the linear solver    
            if (ite .ge. 1) then
              call linsol_doneStructure (p_rsolverNode)
            end if
            
            ! Release the solver node and all subnodes attached to it (if at all):
            call linsol_releaseSolver (p_rsolverNode)
            
            select case(ilinearsolver)
            case (1)
              do ilevel = NLMIN+1,NLMAX
                ! Delete the interlevel projection structure
                call mlprj_doneProjection (Rlevel(ilevel)%rprojection)
              end do
            
            end select
                        
            if (ioutputlevel .ge. 2) then
              call output_lbrk()
            end if
            
            if (ioutputlevel .ge. 1) then
              call output_line("#mg    = "//sys_siL(nmg,10))
              call output_line("#nl    = "//sys_siL(ite-1,10))
            end if
            
            if (ioutputlevel .ge. 2) then
              call output_lbrk()
            end if
            
            if (nmg .eq. 0) then
              call output_line("|y-z|  = -")
              call output_line("|u|    = -")
              call output_line("J(y,u) = -")
              cycle
            end if
            
            ! Calculate the error in y:
            rcollection%IquickAccess(1) = rtargetFlow%itype
            rcollection%p_rvectorQuickAccess1 => rtargetFlow%rvector
            call pperr_scalar (rvectorBlock%RvectorBlock(1),PPERR_L2ERROR,Derror(1),&
                              getReferenceFunctionX_2D,rcollection)
            call pperr_scalar (rvectorBlock%RvectorBlock(2),PPERR_L2ERROR,Derror(2),&
                              getReferenceFunctionY_2D,rcollection)
            Derror(3) = sqrt(0.5_DP*(Derror(1)**2+Derror(2)**2))
            if (ioutputlevel .ge. 1) then
              call output_line("|y-z|  = "//sys_sdEL(Derror(3),15))
            end if

            call pperr_scalar (rvectorBlock%RvectorBlock(4),PPERR_L2ERROR,Derror(4))
            call pperr_scalar (rvectorBlock%RvectorBlock(5),PPERR_L2ERROR,Derror(5))
            Derror(6) =  sqrt(0.5_DP*(Derror(4)**2+Derror(5)**2))
            
            if (ioutputlevel .ge. 1) then
              call output_line("|u|    = "//sys_sdEL(Derror(6)/rparams%dalpha,15))
            
              call output_line("J(y,u) = "//&
                sys_sdEL(0.5_DP*Derror(3)**2 + 0.5_DP*rparams%dalpha*Derror(6)**2,15))
            end if

            ! That's it, rvectorBlock now contains our solution. We can now
            ! start the postprocessing. 
            ! Start UCD export to GMV file:
            call ucd_startGMV (rexport,UCD_FLAG_STANDARD,Rlevel(NLMAX)%rtriangulation,&
                              'gmv/unst2d_'//trim(sys_siL(imaxre(icurrentmaxre),10))//&
                              '_'//trim(sys_siL(ire(icurrentre),10))//'.gmv')
                               
            call spdp_stdProjectionToP1Q1Scalar (rvectorBlock%RvectorBlock(1),&
              rvectorOutputX,rdiscrOutput)
            call spdp_stdProjectionToP1Q1Scalar (rvectorBlock%RvectorBlock(2),&
              rvectorOutputY,rdiscrOutput)
              
            ! Create a discretisation for Q1.
            call spdiscr_initBlockDiscr (rdiscrOutBlock,2,&
                                          Rlevel(NLMAX)%rtriangulation, rboundary)
            call spdiscr_duplicateDiscrSc (rdiscrOutput,&
                rdiscrOutBlock%RspatialDiscr(1), .true.)
            call spdiscr_duplicateDiscrSc (rdiscrOutput,&
                rdiscrOutBlock%RspatialDiscr(2), .true.)
            call lsysbl_createVecBlockByDiscr(rdiscrOutBlock,rvecout)
            call lsyssc_duplicateVector(rvectorOutputX,rvecout%RvectorBlock(1),&
                LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
            call lsyssc_duplicateVector(rvectorOutputY,rvecout%RvectorBlock(2),&
                LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

            ! Assemble the BC for Q1 and implement them
            call getBC (rdiscrOutBlock,ibcType,rdiscreteBCout,.false.)

            call vecfil_discreteBCsol (rvecout,rdiscreteBCout)
            
            ! Release temporary data
            call lsysbl_releaseVector (rvecout)
            call spdiscr_releaseBlockDiscr(rdiscrOutBlock)
            call bcasm_releaseDiscreteBC (rdiscreteBCout)
            
            call lsyssc_getbase_double (rvectorOutputX,p_DdataX)
            call lsyssc_getbase_double (rvectorOutputY,p_DdataY)
            
            call ucd_addVarVertBasedVec (rexport, 'primalvel', p_DdataX, p_DdataY)

            call lsyssc_getbase_double (rvectorBlock%RvectorBlock(3),p_Ddata)
            call ucd_addVariableElementBased (rexport,'primalP',UCD_VAR_STANDARD, p_Ddata)
            
            call spdp_stdProjectionToP1Q1Scalar (rvectorBlock%RvectorBlock(4),&
                rvectorOutputX,rdiscrOutput)
            call spdp_stdProjectionToP1Q1Scalar (rvectorBlock%RvectorBlock(5),&
                rvectorOutputY,rdiscrOutput)
            call ucd_addVariableVertexBased (rexport,'dualvelX',UCD_VAR_STANDARD, p_DdataX)
            call ucd_addVariableVertexBased (rexport,'dualvelY',UCD_VAR_STANDARD, p_DdataY)

            call lsyssc_getbase_double (rvectorBlock%RvectorBlock(6),p_Ddata)
            call ucd_addVariableElementBased (rexport,'dualP',UCD_VAR_STANDARD, p_Ddata)
            
            call lsyssc_copyVector (rvectorBlock%RvectorBlock(4),rvectorTmp)
            call lsyssc_scaleVector (rvectorTmp,-1.0_DP/rparams%dalpha)
            if (rparams%bboundsActive) then
              call projectControlTimestep (rvectorTmp,rparams%dmin1,rparams%dmax1)
            end if
            call spdp_stdProjectionToP1Q1Scalar (rvectorTmp,rvectorOutputX,rdiscrOutput)
            
            call lsyssc_copyVector (rvectorBlock%RvectorBlock(5),rvectorTmp)
            call lsyssc_scaleVector (rvectorTmp,-1.0_DP/rparams%dalpha)
            if (rparams%bboundsActive) then
              call projectControlTimestep (rvectorTmp,rparams%dmin2,rparams%dmax2)
            end if
            call spdp_stdProjectionToP1Q1Scalar (rvectorTmp,rvectorOutputY,rdiscrOutput)
            
            call ucd_addVariableVertexBased (rexport,'controlX',UCD_VAR_STANDARD, p_DdataX)
            call ucd_addVariableVertexBased (rexport,'controlY',UCD_VAR_STANDARD, p_DdataY)
            
            ! residual
            call spdp_stdProjectionToP1Q1Scalar (rdefectBlock%RvectorBlock(1),&
                rvectorOutputX,rdiscrOutput)
            call spdp_stdProjectionToP1Q1Scalar (rdefectBlock%RvectorBlock(2),&
                rvectorOutputY,rdiscrOutput)
            call ucd_addVariableVertexBased (rexport,'respx',UCD_VAR_STANDARD, p_DdataX)
            call ucd_addVariableVertexBased (rexport,'respy',UCD_VAR_STANDARD, p_DdataY)

            call spdp_stdProjectionToP1Q1Scalar (rdefectBlock%RvectorBlock(4),&
                rvectorOutputX,rdiscrOutput)
            call spdp_stdProjectionToP1Q1Scalar (rdefectBlock%RvectorBlock(5),&
                rvectorOutputY,rdiscrOutput)
            call ucd_addVariableVertexBased (rexport,'resdx',UCD_VAR_STANDARD, p_DdataX)
            call ucd_addVariableVertexBased (rexport,'resdy',UCD_VAR_STANDARD, p_DdataY)
            
            ! rhs
            call spdp_stdProjectionToP1Q1Scalar (rrhsBlock%RvectorBlock(1),&
                rvectorOutputX,rdiscrOutput)
            call spdp_stdProjectionToP1Q1Scalar (rrhsBlock%RvectorBlock(2),&
                rvectorOutputY,rdiscrOutput)
            call ucd_addVariableVertexBased (rexport,'rhspx',UCD_VAR_STANDARD, p_DdataX)
            call ucd_addVariableVertexBased (rexport,'rhspy',UCD_VAR_STANDARD, p_DdataY)

            call spdp_stdProjectionToP1Q1Scalar (rrhsBlock%RvectorBlock(4),&
                rvectorOutputX,rdiscrOutput)
            call spdp_stdProjectionToP1Q1Scalar (rrhsBlock%RvectorBlock(5),&
                rvectorOutputY,rdiscrOutput)
            call ucd_addVariableVertexBased (rexport,'rhsdx',UCD_VAR_STANDARD, p_DdataX)
            call ucd_addVariableVertexBased (rexport,'rhsdy',UCD_VAR_STANDARD, p_DdataY)

            ! Error
            !call lsyssc_clearVector (rvectorOutputX)
            !call lsyssc_clearVector (rvectorOutputY)
            !rcollection%p_rvectorQuickAccess1 => rvectorBlock
            !call anprj_discrDirect (rvectorOutputX,getErrorX_2D,rcollection,2)
            !call anprj_discrDirect (rvectorOutputY,getErrorY_2D,rcollection,2)
            !call ucd_addVariableVertexBased (rexport,'errorx',UCD_VAR_STANDARD, p_DdataX)
            !call ucd_addVariableVertexBased (rexport,'errory',UCD_VAR_STANDARD, p_DdataY)
            call lsysbl_clearVector (rdefectBlock)
            rcollection%p_rvectorQuickAccess1 => rtargetFlow%rvector
            rcollection%p_rvectorQuickAccess2 => rvectorBlock
            call anprj_discrDirect (rdefectBlock%RvectorBlock(1),getErrorX_2D,rcollection,2)
            call anprj_discrDirect (rdefectBlock%RvectorBlock(2),getErrorY_2D,rcollection,2)
            call spdp_stdProjectionToP1Q1Scalar (rdefectBlock%RvectorBlock(1),&
                rvectorOutputX,rdiscrOutput)
            call spdp_stdProjectionToP1Q1Scalar (rdefectBlock%RvectorBlock(2),&
                rvectorOutputY,rdiscrOutput)
            call ucd_addVariableVertexBased (rexport,'errorx',UCD_VAR_STANDARD, p_DdataX)
            call ucd_addVariableVertexBased (rexport,'errory',UCD_VAR_STANDARD, p_DdataY)
            
            ! Write the file to disc, that's it.
            call ucd_write (rexport)
            call ucd_release (rexport)

            call lsyssc_releaseVector (rvectorOutputY)
            call lsyssc_releaseVector (rvectorOutputX)
            
          end do
          
        end do
        
      end do
      
      ! Release our discrete version of the boundary conditions
      do ilevel=NLMIN,NLMAX
        call bcasm_releaseDiscreteBC (Rlevel(ilevel)%rdiscreteBC)
      end do
      
      ! We are finished - but not completely!
      ! Now, clean up so that all the memory is available again.
      
      ! Release the block matrix/vectors
      call lsyssc_releaseVector (rvectorTmp)
      call lsysbl_releaseVector (rdefectBlock)
      call lsysbl_releaseVector (rvectorBlock)
      call lsysbl_releaseVector (rrhsBlock)

      if (ilinearsolver .ne. 1) then
        do ilevel=NLMIN+1,NLMAX
          call mlprj_doneProjection (Rlevel(ilevel)%rprojection)
        end do
      end if

      do ilevel=NLMIN,NLMAX
      
        call lsysbl_releaseVector (Rlevel(ilevel)%rtempVector)
        call lsysbl_releaseMatrix (RmatrixBlock(ilevel))
        call lsyssc_releaseMatrix (Rlevel(ilevel)%rmatrixB2)
        call lsyssc_releaseMatrix (Rlevel(ilevel)%rmatrixB1)
        call lsyssc_releaseMatrix (Rlevel(ilevel)%rmatrixMass)
        call lsyssc_releaseMatrix (Rlevel(ilevel)%rmatrixLaplace)
        call lsyssc_releaseMatrix (Rlevel(ilevel)%rmatrixDiagP)
        call spdiscr_releaseBlockDiscr(Rlevel(ilevel)%rdiscretisation)
      
        ! Release the triangulation. 
        call tria_done (Rlevel(ilevel)%rtriangulation)
      
      end do
      
      ! Release the discretisation structure and all spatial discretisation
      ! structures in it.
      call spdiscr_releaseDiscr(rdiscrOutput)
      
      ! Clean up the target flow
      call doneTargetFlow (rtargetFlow)
      
    end do
    
    ! Finally release the domain, that's it.
    call boundary_release (rboundary)
    
  end subroutine

end module

! Tests: Den L2-Fehler zur echten Lösung bei verschiedenen Konfigurationen.
! |u| berechnen, J() berechnen.
! Konvergenzverlauf aufzeichnen für DefCorr/Newton bei versch. Konfigurationen.
! -> Mit/Ohne Beschränkung der Kontrolle.

