!##############################################################################
!# ****************************************************************************
!# <name> spacematvecassembly </name>
!# ****************************************************************************
!#
!# This module realises the basic martix-vector multiplication which is
!# used during the forward/backward iteration as well as the assembly
!# of the linearised matrices. It can be described as the heart of the
!# discretisation as it has to respect all the different settings in the
!# discretisation.
!#
!# </purpose>
!##############################################################################

module spacematvecassembly
  
  use fsystem
  use genoutput
  
  use spatialdiscretisation
  use timediscretisation
  use linearsystemscalar
  use linearsystemblock
  
  use scalarpde
  use linearformevaluation
  use bilinearformevaluation
  use feevaluation2
  use blockmatassemblybase
  use blockmatassembly
  use collection
  
  use spacetimevectors
  use analyticsolution
  
  use structuresdiscretisation
  use structuresoptcontrol
  use structuresgeneral
  use assemblytemplates
  
  use kktsystemspaces
  
  use user_callback
    
  implicit none
  
  private
  
!<constants>

!<constantblock description = "A list of operator tyes for the assembly">

  ! Navier-Stokes operator, primal equation
  integer, parameter :: OPTP_PRIMAL = 0

  ! Simple linearised Navier-Stokes, primal equation
  integer, parameter :: OPTP_PRIMALLIN_SIMPLE = 1
  
  ! Full linearised Navier-Stokes (Newton), primal equation
  integer, parameter :: OPTP_PRIMALLIN = 2
  
  ! Navier-Stokes operator, dual equation
  integer, parameter :: OPTP_DUAL = 3

  ! Navier-Stokes operator, linearised dual equation
  integer, parameter :: OPTP_DUALLIN_SIMPLE = 4

  ! Navier-Stokes operator, linearised dual equation (Newton)
  integer, parameter :: OPTP_DUALLIN = 5
  
  ! Full linearised Navier-Stokes (Newton), dual equation,
  ! nonlinear part in the RHS.
  integer, parameter :: OPTP_DUALLIN_RHS = 6

!</constantblock>

!</constants>

!<types>

!<typeblock>

  ! This type encapsules structures necessary for the
  ! assembly of space-time operators.
  type t_spacetimeOperatorAsm
  
    ! Space discretisation.
    type(t_blockDiscretisation), pointer :: p_rspaceDiscr => null()
    
    ! Time discretisation. Defines the time stepping scheme
    type(t_timeDiscretisation), pointer :: p_rtimeDiscr => null()
    
    ! Physics of the problem
    type(t_settings_physics), pointer :: p_rphysics => null()
    
    ! Parameters for the discretisation in space
    type(t_settings_discr), pointer :: p_rsettingsDiscr => null()

    ! Optimal-control parameters
    type(t_settings_optcontrol), pointer :: p_rsettingsOptControl => null()
    
    ! Assembly templates corresponding to the above space discretisation
    type(t_staticSpaceAsmTemplates), pointer :: p_rasmTemplates => null()

    ! Reference to the analytic solution defining the RHS of the primal equation
    type(t_anSolution), pointer :: p_rrhsPrimal => null()

    ! Reference to the analytic solution defining the RHS of the dual equation
    type(t_anSolution), pointer :: p_rrhsDual => null()

    ! Reference to the analytic solution defining the target
    type(t_anSolution), pointer :: p_rtargetFlow => null()

    ! Reference to global data
    type(t_globalData), pointer :: p_rglobalData => null()
    
    ! Reference to the debug flags
    type(t_optcDebugFlags), pointer :: p_rdebugFlags => null()
  end type

!</typeblock>

  public :: t_spacetimeOperatorAsm

!<typeblock>

  ! Pointer to t_spacetimeOperator for being passed to callback routines.
  type p_t_spacetimeOperatorAsm
    type(t_spacetimeOperatorAsm), pointer :: p_rspatialOperatorAsm
  end type

!</typeblock>

!<typeblock>

  ! Temporary assembly data
  type t_assemblyTempDataSpace
  
    private
  
    ! A temporary block matrix which can be used for
    ! intermediate calculations.
    type(t_matrixBlock) :: rmatrix

    ! A set of temporary block vectors which can be used for
    ! intermediate calculations.
    type(t_vectorBlock), dimension(1) :: rvector

  end type

!</typeblock>

public :: t_assemblyTempDataSpace

!</types>

  ! Allocate temporary memory
  public :: smva_allocTempData
  
  ! Release temprary memory
  public :: smva_releaseTempData
  
  ! Calculate the defect in the primal equation
  public :: smva_getDef_primal
  
  ! Calculate the defect in the dual equation
  public :: smva_getDef_dual
  
  ! Calculate the defect in the linearised primal equation
  public :: smva_getDef_primalLin
  
  ! Calculate the defect in the linearised dual equation
  public :: smva_getDef_dualLin
  
  ! Assemble the matrices of the operator of the linearised primal equation
  public :: smva_assembleMatrix_primal
  
  ! Assemble the matrices of the operator of the linearised dual equation
  public :: smva_assembleMatrix_dual
  
contains

  ! ***************************************************************************
  
!<subroutine>

  subroutine smva_allocTempData (rspatialOperatorAsm,rtempData)

!<description>
  ! Allocates temporary data for the assembly of an operator
!</description>

!<input>
  ! Definition of the operator. 
  type(t_spacetimeOperatorAsm) :: rspatialOperatorAsm
!</input>

!<output>
  ! Structure containing temporary data.
  type(t_assemblyTempDataSpace), intent(out) :: rtempData
!</output>

!</subroutine>

    integer :: i
    
    ! Create temp vectors
    do i=1,ubound(rtempData%Rvector,1)
      call lsysbl_createVectorBlock(rspatialOperatorAsm%p_rspaceDiscr,rtempData%Rvector(i),.false.)
    end do
  
    ! Create a full matrix
    call lsysbl_createMatBlockByDiscr (rspatialOperatorAsm%p_rspaceDiscr,rtempData%rmatrix)
    
    ! Fill it with data from the assembly template structure.
    ! The structure depends on the current equation...
    select case (rspatialOperatorAsm%p_rphysics%cequation)

    ! *************************************************************
    ! Stokes/Navier Stokes.
    ! *************************************************************
    case (0,1)
    
      ! ---------------------------------------------------
      ! 2x2 block for the velocity
      call lsyssc_duplicateMatrix (&
          rspatialOperatorAsm%p_rasmTemplates%rmatrixTemplateFEM,&
          rtempData%rmatrix%RmatrixBlock(1,1),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (&
          rspatialOperatorAsm%p_rasmTemplates%rmatrixTemplateFEMoffdiag,&
          rtempData%rmatrix%RmatrixBlock(2,1),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (&
          rspatialOperatorAsm%p_rasmTemplates%rmatrixTemplateFEMoffdiag,&
          rtempData%rmatrix%RmatrixBlock(1,2),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (&
          rspatialOperatorAsm%p_rasmTemplates%rmatrixTemplateFEM,&
          rtempData%rmatrix%RmatrixBlock(2,2),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
          
      ! ---------------------------------------------------
      ! B-matrices
      call lsyssc_duplicateMatrix (&
          rspatialOperatorAsm%p_rasmTemplates%rmatrixTemplateGradient,&
          rtempData%rmatrix%RmatrixBlock(1,3),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (&
          rspatialOperatorAsm%p_rasmTemplates%rmatrixTemplateGradient,&
          rtempData%rmatrix%RmatrixBlock(2,3),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
      
      ! ---------------------------------------------------
      ! D-matrices
      call lsyssc_duplicateMatrix (&
          rspatialOperatorAsm%p_rasmTemplates%rmatrixTemplateDivergence,&
          rtempData%rmatrix%RmatrixBlock(3,1),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)

      call lsyssc_duplicateMatrix (&
          rspatialOperatorAsm%p_rasmTemplates%rmatrixTemplateDivergence,&
          rtempData%rmatrix%RmatrixBlock(3,2),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
      
      ! ---------------------------------------------------
      ! Pressure block
      call lsyssc_duplicateMatrix (&
          rspatialOperatorAsm%p_rasmTemplates%rmatrixTemplateFEMPressure,&
          rtempData%rmatrix%RmatrixBlock(3,3),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
         
    case default
      
      call output_line("Unknown equation",&
          OU_CLASS_ERROR,OU_MODE_STD,"smva_allocTempData")
      call sys_halt()
    
    end select
    
    ! Update the matrix structure
    call lsysbl_updateMatStrucInfo(rtempData%rmatrix)
        
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine smva_releaseTempData (rtempData)

!<description>
  ! Releases allocated memory.
!</description>

!<inputoutput>
  ! Structure containing temporary data to be released
  type(t_assemblyTempDataSpace), intent(inout) :: rtempData
!</inputoutput>

!</subroutine>

    integer :: i
    
    ! Release everything
    call lsysbl_releaseMatrix (rtempData%rmatrix)

    do i=1,ubound(rtempData%Rvector,1)
      call lsysbl_releaseVector (rtempData%Rvector(i))
    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine smva_apply_primal (rspatialOperatorAsm,rprimalSol,idofTime,rdest,rtempData)
  
!<description>
  ! Applies the operator of the primal equation at DOF idofTime in time
  ! to the primal solurion rprimalSol.
!</description>

!<input>
  ! Definition of the operator
  type(t_spacetimeOperatorAsm) :: rspatialOperatorAsm

  ! Structure that defines the current primal solution.
  type(t_primalSpace), intent(inout) :: rprimalSol

  ! Number of the DOF in time which should be calculated into rdest.
  integer, intent(in) :: idofTime
!</input>

!<inputoutput>
  ! Structure with temporary assembly data.
  type(t_assemblyTempDataSpace), intent(inout) :: rtempData
    
  ! Destination vector
  type(t_vectorBlock), intent(inout) :: rdest
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_rvector
    real(DP) :: dtheta, dtstep, dtimeend, dtimestart
    
    ! Clear the output vector
    call lsysbl_clearVector (rdest)

    ! Timestepping technique?
    select case (rspatialOperatorAsm%p_rtimeDiscr%ctype)
    
    ! ***********************************************************
    ! Standard Theta one-step scheme.
    ! ***********************************************************
    case (TDISCR_ONESTEPTHETA)
    
      ! Theta-scheme identifier
      dtheta = rspatialOperatorAsm%p_rtimeDiscr%dtheta

      ! Characteristics of the current timestep.
      call tdiscr_getTimestep(rspatialOperatorAsm%p_rtimeDiscr,idofTime-1,&
          dtimeend,dtstep,dtimestart)

      ! itag=0: old 1-step scheme.
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      select case (rspatialOperatorAsm%p_rtimeDiscr%itag)
      
      ! ***********************************************************
      ! itag=0: old/standard 1-step scheme.
      ! ***********************************************************
      case (0)

        call output_line("Old 1-step-scheme not implemented",&
            OU_CLASS_ERROR,OU_MODE_STD,"smva_getRhs_Primal")
        call sys_halt()

      ! ***********************************************************
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      ! ***********************************************************
      case (1)

        ! Which equation do we have?    
        select case (rspatialOperatorAsm%p_rphysics%cequation)

        ! *************************************************************
        ! Stokes/Navier Stokes.
        ! *************************************************************
        case (0,1)
      
          ! ***********************************************
          ! PREVIOUS TIMESTEP
          ! ***********************************************
          if (idofTime .gt. 1) then

            ! ===============================================
            ! Prepare the linear parts of the matrix.
            ! ===============================================

            ! Clear the temporary matrix
            call lsysbl_clearMatrix (rtempData%rmatrix)
        
            ! -----------------------------------------
            ! Mass matrix for timestepping
            call smva_getMassMatrix (rspatialOperatorAsm,rtempData%rmatrix,1.0_DP-dtstep)
            
            ! -----------------------------------------
            ! Laplace -- if the viscosity is constant
            if (rspatialOperatorAsm%p_rphysics%cviscoModel .eq. 0) then
              call smva_getLaplaceMatrix (&
                  rspatialOperatorAsm,rtempData%rmatrix,rspatialOperatorAsm%p_rphysics%dnuConst)
            else
              call output_line("Nonconstant viscosity not supported.",&
                  OU_CLASS_ERROR,OU_MODE_STD,"smva_apply_primal")
              call sys_halt()
            end if

            ! -----------------------------------------
            ! Realise the defect
            call sptivec_getVectorFromPool (rprimalSol%p_rvectorAccess,idofTime-1,p_rvector)
            call lsysbl_blockMatVec (rtempData%rmatrix, p_rvector, rdest, -1.0_DP, 1.0_DP)
          end if

          ! ***********************************************
          ! CURRENT TIMESTEP
          ! ***********************************************

          ! ===============================================
          ! Prepare the linear parts of the matrix.
          ! ===============================================
          
          ! Clear the temporary matrix
          call lsysbl_clearMatrix (rtempData%rmatrix)
          
          ! -----------------------------------------
          ! Mass matrix for timestepping
          call smva_getMassMatrix (rspatialOperatorAsm,rtempData%rmatrix,dtstep)
          
          ! -----------------------------------------
          ! Laplace -- if the viscosity is constant
          if (rspatialOperatorAsm%p_rphysics%cviscoModel .eq. 0) then
            call smva_getLaplaceMatrix (&
                rspatialOperatorAsm,rtempData%rmatrix,rspatialOperatorAsm%p_rphysics%dnuConst)
          else
            call output_line("Nonconstant viscosity not supported.",&
                OU_CLASS_ERROR,OU_MODE_STD,"smva_apply_primal")
            call sys_halt()
          end if
          
          ! -----------------------------------------
          ! B-matrices
          call smva_getBMatrix (rspatialOperatorAsm,rtempData%rmatrix,1.0_DP)
          
          ! -----------------------------------------
          ! D-matrices
          call smva_getDMatrix (rspatialOperatorAsm,rtempData%rmatrix,1.0_DP)

          ! -----------------------------------------
          ! EOJ-stabilisation
          !if (rspatialOperatorAsm%p_rsettingsDiscr%rstabilConvecPrimal%cupwind .eq. 4) then
          !  call smva_getEOJMatrix (rspatialOperatorAsm,rtempData%rmatrix,1.0_DP)
          !end if
            
          ! ===============================================
          ! Prepare the semilinear parts of the matrix.
          ! ===============================================
          
          ! The semilinear parts of the matrix can be set up with
          ! the block matrix assembly routines. Invoke them using
          ! a separate subroutine
          call smva_getSemilinMat_primal (rspatialOperatorAsm,idofTime,&
              rprimalSol,1.0_DP,rtempData%rmatrix)

          ! -----------------------------------------
          ! Realise the defect with the linear and
          ! semilinear parts of the operator
          ! -----------------------------------------

          call sptivec_getVectorFromPool (rprimalSol%p_rvectorAccess,idofTime,p_rvector)
          call lsysbl_blockMatVec (rtempData%rmatrix, p_rvector, rdest, 1.0_DP, 1.0_DP)

        end select ! Equation
      
      end select ! Timestep sub-scheme
      
    end select ! Timestep scheme

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine smva_getDef_primal (rspatialOperatorAsm,rprimalSol,rcontrol,idofTime,rdest,rtempData)
  
!<description>
  ! Calculates the defect in timestep idofTime of the nonlinear primaö equation
  ! in one DOF in time.
!</description>

!<input>
  ! Definition of the operator
  type(t_spacetimeOperatorAsm) :: rspatialOperatorAsm

  ! Structure that defines the current primal solution.
  type(t_primalSpace), intent(inout) :: rprimalSol

  ! Structure that defines the current control.
  type(t_controlSpace), intent(inout) :: rcontrol

  ! Number of the DOF in time which should be calculated into rdest.
  integer, intent(in) :: idofTime
!</input>

!<inputoutput>
  ! Structure with temporary assembly data.
  type(t_assemblyTempDataSpace), intent(inout) :: rtempData
    
  ! Destination vector
  type(t_vectorBlock), intent(inout) :: rdest
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_rvector
    real(DP) :: dtheta, dtstep, dtimeend, dtimestart
    
    ! Clear the output vector
    call lsysbl_clearVector (rdest)

    ! Timestepping technique?
    select case (rspatialOperatorAsm%p_rtimeDiscr%ctype)
    
    ! ***********************************************************
    ! Standard Theta one-step scheme.
    ! ***********************************************************
    case (TDISCR_ONESTEPTHETA)
    
      ! Theta-scheme identifier
      dtheta = rspatialOperatorAsm%p_rtimeDiscr%dtheta

      ! Characteristics of the current timestep.
      call tdiscr_getTimestep(rspatialOperatorAsm%p_rtimeDiscr,idofTime-1,&
          dtimeend,dtstep,dtimestart)

      ! itag=0: old 1-step scheme.
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      select case (rspatialOperatorAsm%p_rtimeDiscr%itag)
      
      ! ***********************************************************
      ! itag=0: old/standard 1-step scheme.
      ! ***********************************************************
      case (0)

        call output_line("Old 1-step-scheme not implemented",&
            OU_CLASS_ERROR,OU_MODE_STD,"smva_getRhs_Primal")
        call sys_halt()

      ! ***********************************************************
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      ! ***********************************************************
      case (1)

        ! Which equation do we have?    
        select case (rspatialOperatorAsm%p_rphysics%cequation)

        ! *************************************************************
        ! Stokes/Navier Stokes.
        ! *************************************************************
        case (0,1)
      
          ! ***********************************************
          ! PREVIOUS TIMESTEP
          ! ***********************************************
          if (idofTime .gt. 1) then

            ! ===============================================
            ! RHS assembly
            ! ===============================================
            
            ! nothing to do here

            ! ===============================================
            ! Prepare the linear parts of the matrix.
            ! ===============================================

            ! Clear the temporary matrix
            call lsysbl_clearMatrix (rtempData%rmatrix)
          
            ! -----------------------------------------
            ! Mass matrix for timestepping
            call smva_getMassMatrix (rspatialOperatorAsm,rtempData%rmatrix,1.0_DP-dtstep)
            
            ! -----------------------------------------
            ! Laplace -- if the viscosity is constant
            if (rspatialOperatorAsm%p_rphysics%cviscoModel .eq. 0) then
              call smva_getLaplaceMatrix (&
                  rspatialOperatorAsm,rtempData%rmatrix,rspatialOperatorAsm%p_rphysics%dnuConst)
            else
              call output_line("Nonconstant viscosity not supported.",&
                  OU_CLASS_ERROR,OU_MODE_STD,"smva_getDef_primal")
              call sys_halt()
            end if

            ! -----------------------------------------
            ! Realise the defect
            call sptivec_getVectorFromPool (rprimalSol%p_rvectorAccess,idofTime-1,p_rvector)
            call lsysbl_blockMatVec (rtempData%rmatrix, p_rvector, rdest, -1.0_DP, 1.0_DP)
            
          end if

          ! ***********************************************
          ! CURRENT TIMESTEP
          ! ***********************************************

          ! ===============================================
          ! RHS assembly
          ! ===============================================
          
          ! Get the RHS vector.
          call smva_getRhs_Primal (rspatialOperatorAsm,idofTime,rcontrol,1.0_DP,rdest)
          
          ! ===============================================
          ! Prepare the linear parts of the matrix.
          ! ===============================================
          
          ! Clear the temporary matrix
          call lsysbl_clearMatrix (rtempData%rmatrix)
          
          ! -----------------------------------------
          ! Mass matrix for timestepping
          call smva_getMassMatrix (rspatialOperatorAsm,rtempData%rmatrix,dtstep)
          
          ! -----------------------------------------
          ! Laplace -- if the viscosity is constant
          if (rspatialOperatorAsm%p_rphysics%cviscoModel .eq. 0) then
            call smva_getLaplaceMatrix (&
                rspatialOperatorAsm,rtempData%rmatrix,rspatialOperatorAsm%p_rphysics%dnuConst)
          else
            call output_line("Nonconstant viscosity not supported.",&
                OU_CLASS_ERROR,OU_MODE_STD,"smva_getDef_primal")
            call sys_halt()
          end if
          
          ! -----------------------------------------
          ! B-matrices
          call smva_getBMatrix (rspatialOperatorAsm,rtempData%rmatrix,1.0_DP)
          
          ! -----------------------------------------
          ! D-matrices
          call smva_getDMatrix (rspatialOperatorAsm,rtempData%rmatrix,1.0_DP)

          ! -----------------------------------------
          ! EOJ-stabilisation
          !if (rspatialOperatorAsm%p_rsettingsDiscr%rstabilConvecPrimal%cupwind .eq. 4) then
          !  call smva_getEOJMatrix (rspatialOperatorAsm,rtempData%rmatrix,1.0_DP)
          !end if
            
          ! ===============================================
          ! Prepare the semilinear parts of the matrix.
          ! ===============================================
          
          ! The semilinear parts of the matrix can be set up with
          ! the block matrix assembly routines. Invoke them using
          ! a separate subroutine
          call smva_getSemilinMat_primal (rspatialOperatorAsm,idofTime,&
              rprimalSol,1.0_DP,rtempData%rmatrix)

          ! -----------------------------------------
          ! Realise the defect with the linear and
          ! semilinear parts of the operator
          ! -----------------------------------------

          call sptivec_getVectorFromPool (rprimalSol%p_rvectorAccess,idofTime,p_rvector)
          call lsysbl_blockMatVec (rtempData%rmatrix, p_rvector, rdest, -1.0_DP, 1.0_DP)

        end select ! Equation
      
      end select ! Timestep sub-scheme
      
    end select ! Timestep scheme

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine smva_getDef_dual (rspatialOperatorAsm,rprimalSol,rdualSol,idofTime,&
      rdest,rtempData)
  
!<description>
  ! Calculates the defect in timestep idofTime of the nonlinear dual equation
  ! in one DOF in time.
!</description>

!<input>
  ! Definition of the operator
  type(t_spacetimeOperatorAsm) :: rspatialOperatorAsm

  ! Space-time vector that describes the current primal solution.
  type(t_primalSpace), intent(inout) :: rprimalSol

  ! Structure that defines the vector rdualSol.
  type(t_dualSpace), intent(inout) :: rdualSol
  
  ! Number of the DOF in time which should be calculated into rdest.
  integer, intent(in) :: idofTime
!</input>

!<inputoutput>
  ! Structure with temporary assembly data.
  type(t_assemblyTempDataSpace), intent(inout) :: rtempData
    
  ! Destination vector
  type(t_vectorBlock), intent(inout) :: rdest
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_rvector
    real(DP) :: dtheta, dtstep, dtimeend, dtimestart
    
    ! Clear the output vector
    call lsysbl_clearVector (rdest)

    ! Timestepping technique?
    select case (rspatialOperatorAsm%p_rtimeDiscr%ctype)
    
    ! ***********************************************************
    ! Standard Theta one-step scheme.
    ! ***********************************************************
    case (TDISCR_ONESTEPTHETA)
    
      ! Theta-scheme identifier
      dtheta = rspatialOperatorAsm%p_rtimeDiscr%dtheta

      ! Characteristics of the current timestep.
      call tdiscr_getTimestep(rspatialOperatorAsm%p_rtimeDiscr,idofTime-1,&
          dtimeend,dtstep,dtimestart)

      ! itag=0: old 1-step scheme.
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      select case (rspatialOperatorAsm%p_rtimeDiscr%itag)
      
      ! ***********************************************************
      ! itag=0: old/standard 1-step scheme.
      ! ***********************************************************
      case (0)

        call output_line("Old 1-step-scheme not implemented",&
            OU_CLASS_ERROR,OU_MODE_STD,"smva_getDef_dual")
        call sys_halt()

      ! ***********************************************************
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      ! ***********************************************************
      case (1)

        ! Which equation do we have?    
        select case (rspatialOperatorAsm%p_rphysics%cequation)

        ! *************************************************************
        ! Stokes/Navier Stokes.
        ! *************************************************************
        case (0,1)

          ! ***********************************************
          ! NEXT TIMESTEP
          ! ***********************************************
          if (idofTime .lt. rprimalSol%p_rvector%NEQtime) then

            ! ===============================================
            ! RHS assembly
            ! ===============================================
          
            ! nothing to do here

            ! ===============================================
            ! Prepare the linear parts of the matrix.
            ! ===============================================

            ! Clear the temporary matrix
            call lsysbl_clearMatrix (rtempData%rmatrix)

            ! -----------------------------------------
            ! Mass matrix for timestepping
            call smva_getMassMatrix (rspatialOperatorAsm,rtempData%rmatrix,1.0_DP-dtstep)
            
            ! -----------------------------------------
            ! Laplace -- if the viscosity is constant
            if (rspatialOperatorAsm%p_rphysics%cviscoModel .eq. 0) then
              call smva_getLaplaceMatrix (&
                  rspatialOperatorAsm,rtempData%rmatrix,rspatialOperatorAsm%p_rphysics%dnuConst)
            else
              call output_line("Nonconstant viscosity not supported.",&
                  OU_CLASS_ERROR,OU_MODE_STD,"smva_getDef_dual")
              call sys_halt()
            end if

            ! -----------------------------------------
            ! Realise the defect
            call sptivec_getVectorFromPool (rdualSol%p_rvectorAccess,idofTime-1,p_rvector)
            call lsysbl_blockMatVec (rtempData%rmatrix, p_rvector, rdest, -1.0_DP, 1.0_DP)
          end if
            
          ! ***********************************************
          ! CURRENT TIMESTEP
          ! ***********************************************

          ! ===============================================
          ! RHS assembly
          ! ===============================================
          
          ! Get the RHS vector
          call smva_getRhs_Dual (rspatialOperatorAsm,idofTime,rprimalSol,1.0_DP,rdest)

          ! ===============================================
          ! Prepare the linear parts of the matrix.
          ! ===============================================

          ! Clear the temporary matrix
          call lsysbl_clearMatrix (rtempData%rmatrix)

          ! -----------------------------------------
          ! Mass matrix for timestepping
          call smva_getMassMatrix (rspatialOperatorAsm,rtempData%rmatrix,dtstep)
          
          ! -----------------------------------------
          ! Laplace -- if the viscosity is constant
          if (rspatialOperatorAsm%p_rphysics%cviscoModel .eq. 0) then
            call smva_getLaplaceMatrix (rspatialOperatorAsm,&
                rtempData%rmatrix,rspatialOperatorAsm%p_rphysics%dnuConst)
          else
            call output_line("Nonconstant viscosity not supported.",&
                OU_CLASS_ERROR,OU_MODE_STD,"smva_getDef_dual")
            call sys_halt()
          end if
          
          ! -----------------------------------------
          ! B-matrices
          call smva_getBMatrix (rspatialOperatorAsm,rtempData%rmatrix,1.0_DP)
          
          ! -----------------------------------------
          ! D-matrices
          call smva_getDMatrix (rspatialOperatorAsm,rtempData%rmatrix,1.0_DP)

          ! -----------------------------------------
          ! EOJ-stabilisation
          !if (rspatialOperatorAsm%p_rsettingsDiscr%rstabilConvecPrimal%cupwind .eq. 4) then
          !  call smva_getEOJMatrix (rspatialOperatorAsm,rtempData%rmatrix,1.0_DP)
          !end if
            
          ! ===============================================
          ! Prepare the semilinear parts of the matrix.
          ! ===============================================
          
          ! The semilinear parts of the matrix can be set up with
          ! the block matrix assembly routines. Invoke them using
          ! a separate subroutine
          call smva_getSemilinMat_Dual (rspatialOperatorAsm,idofTime,&
              rprimalSol,1.0_DP,rtempData%rmatrix)

          ! -----------------------------------------
          ! Realise the defect with the linear and
          ! semilinear parts of the operator
          ! -----------------------------------------

          call sptivec_getVectorFromPool (rdualSol%p_rvectorAccess,idofTime,p_rvector)
          call lsysbl_blockMatVec (rtempData%rmatrix, p_rvector, rdest, -1.0_DP, 1.0_DP)
          
        end select ! Equation
      
      end select ! Timestep sub-scheme
      
    end select ! Timestep scheme

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine smva_getDef_primalLin (rspatialOperatorAsm,&
      rprimalSol,rcontrol,rprimalLinSol,rcontrolLin,idofTime,rdest,rtempData)
  
!<description>
  ! Calculates the defect in timestep idofTime of the linearised primal equation
  ! in one DOF in time.
!</description>

!<input>
  ! Definition of the operator
  type(t_spacetimeOperatorAsm) :: rspatialOperatorAsm

  ! Structure that describes the solution of the forward equation.
  type(t_primalSpace), intent(inout) :: rprimalSol

  ! Structure that describes the solution of the control equation.
  type(t_controlSpace), intent(inout) :: rcontrol

  ! Structure that describes the solution of the linearised primal equation.
  type(t_primalSpace), intent(inout) :: rprimalLinSol

  ! Structure that describes the solution of the linearised control equation.
  type(t_controlSpace), intent(inout) :: rcontrolLin

  ! Number of the DOF in time which should be calculated into rdest.
  integer, intent(in) :: idofTime
!</input>

!<inputoutput>
  ! Structure with temporary assembly data.
  type(t_assemblyTempDataSpace), intent(inout) :: rtempData
    
  ! Destination vector
  type(t_vectorBlock), intent(inout) :: rdest
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_vectorBlock), pointer :: p_rvector
    real(DP) :: dtheta, dtstep, dtimeend, dtimestart
    
    ! Clear the output vector
    call lsysbl_clearVector (rdest)

    ! Timestepping technique?
    select case (rspatialOperatorAsm%p_rtimeDiscr%ctype)
    
    ! ***********************************************************
    ! Standard Theta one-step scheme.
    ! ***********************************************************
    case (TDISCR_ONESTEPTHETA)
    
      ! Theta-scheme identifier
      dtheta = rspatialOperatorAsm%p_rtimeDiscr%dtheta

      ! Characteristics of the current timestep.
      call tdiscr_getTimestep(rspatialOperatorAsm%p_rtimeDiscr,idofTime-1,&
          dtimeend,dtstep,dtimestart)

      ! itag=0: old 1-step scheme.
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      select case (rspatialOperatorAsm%p_rtimeDiscr%itag)
      
      ! ***********************************************************
      ! itag=0: old/standard 1-step scheme.
      ! ***********************************************************
      case (0)

        call output_line("Old 1-step-scheme not implemented",&
            OU_CLASS_ERROR,OU_MODE_STD,"smva_getDef_primalLin")
        call sys_halt()

      ! ***********************************************************
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      ! ***********************************************************
      case (1)

        ! Which equation do we have?    
        select case (rspatialOperatorAsm%p_rphysics%cequation)

        ! *************************************************************
        ! Stokes/Navier Stokes.
        ! *************************************************************
        case (0,1)

          ! ***********************************************
          ! PREVIOUS TIMESTEP
          ! ***********************************************
          if (idofTime .gt. 1) then

            ! ===============================================
            ! Prepare the linear parts of the matrix.
            ! ===============================================

            ! Clear the temporary matrix
            call lsysbl_clearMatrix (rtempData%rmatrix)
          
            ! -----------------------------------------
            ! Mass matrix for timestepping
            call smva_getMassMatrix (rspatialOperatorAsm,rtempData%rmatrix,1.0_DP-dtstep)
            
            ! -----------------------------------------
            ! Laplace -- if the viscosity is constant
            if (rspatialOperatorAsm%p_rphysics%cviscoModel .eq. 0) then
              call smva_getLaplaceMatrix (&
                  rspatialOperatorAsm,rtempData%rmatrix,rspatialOperatorAsm%p_rphysics%dnuConst)
            else
              call output_line("Nonconstant viscosity not supported.",&
                  OU_CLASS_ERROR,OU_MODE_STD,"smva_getDef_primalLin")
              call sys_halt()
            end if

            ! -----------------------------------------
            ! Realise the defect
            call sptivec_getVectorFromPool (rprimalLinSol%p_rvectorAccess,idofTime-1,p_rvector)
            call lsysbl_blockMatVec (rtempData%rmatrix, p_rvector, rdest, -1.0_DP, 1.0_DP)
            
          end if

          ! ***********************************************
          ! CURRENT TIMESTEP
          ! ***********************************************

          ! ===============================================
          ! RHS assembly
          ! ===============================================
          call smva_getRhs_primalLin (rspatialOperatorAsm,idofTime,&
              rcontrol,rcontrolLin,1.0_DP,rdest)
          
          ! ===============================================
          ! Prepare the linear parts of the matrix.
          ! ===============================================
          
          ! Clear the temporary matrix
          call lsysbl_clearMatrix (rtempData%rmatrix)
          
          ! -----------------------------------------
          ! Mass matrix for timestepping
          call smva_getMassMatrix (rspatialOperatorAsm,rtempData%rmatrix,dtstep)
          
          ! -----------------------------------------
          ! Laplace -- if the viscosity is constant
          if (rspatialOperatorAsm%p_rphysics%cviscoModel .eq. 0) then
            call smva_getLaplaceMatrix (&
                rspatialOperatorAsm,rtempData%rmatrix,rspatialOperatorAsm%p_rphysics%dnuConst)
          else
            call output_line("Nonconstant viscosity not supported.",&
                OU_CLASS_ERROR,OU_MODE_STD,"smva_getDef_primalLin")
            call sys_halt()
          end if
          
          ! -----------------------------------------
          ! B-matrices
          call smva_getBMatrix (rspatialOperatorAsm,rtempData%rmatrix,1.0_DP)
          
          ! -----------------------------------------
          ! D-matrices
          call smva_getDMatrix (rspatialOperatorAsm,rtempData%rmatrix,1.0_DP)

          ! -----------------------------------------
          ! EOJ-stabilisation
          !if (rspatialOperatorAsm%p_rsettingsDiscr%rstabilConvecPrimal%cupwind .eq. 4) then
          !  call smva_getEOJMatrix (rspatialOperatorAsm,rtempData%rmatrix,1.0_DP)
          !end if
            
          ! ===============================================
          ! Prepare the semilinear parts of the matrix.
          ! ===============================================
          
          ! The semilinear parts of the matrix can be set up with
          ! the block matrix assembly routines. Invoke them using
          ! a separate subroutine
          call smva_getSemilinMat_primalLin (rspatialOperatorAsm,idofTime,&
              rprimalSol,1.0_DP,.true.,rtempData%rmatrix)

          ! -----------------------------------------
          ! Realise the defect with the linear and
          ! semilinear parts of the operator
          ! -----------------------------------------

          call sptivec_getVectorFromPool (rprimalLinSol%p_rvectorAccess,idofTime,p_rvector)
          call lsysbl_blockMatVec (rtempData%rmatrix, p_rvector, rdest, -1.0_DP, 1.0_DP)

        end select ! Equation
      
      end select ! Timestep sub-scheme
      
    end select ! Timestep scheme

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine smva_getDef_dualLin (rspatialOperatorAsm,&
      rprimalSol,rdualSol,rprimalLinSol,rdualLinSol,idofTime,bfull,rdest,rtempData)
  
!<description>
  ! Calculates the defect in timestep idofTime of the linearised dual equation
  ! in one DOF in time.
!</description>

!<input>
  ! Definition of the operator
  type(t_spacetimeOperatorAsm) :: rspatialOperatorAsm

  ! Space-time vector that describes the solution of the forward equation.
  type(t_primalSpace), intent(inout) :: rprimalSol

  ! Space-time vector that describes the solution of the backward equation.
  type(t_dualSpace), intent(inout) :: rdualSol

  ! Space-time vector that describes the solution of the linearised forward equation.
  type(t_dualSpace), intent(inout) :: rprimalLinSol

  ! Structure that defines the vector rdualLinSol.
  type(t_dualSpace), intent(inout) :: rdualLinSol

  ! Number of the DOF in time which should be calculated into rdest.
  integer, intent(in) :: idofTime

  ! TRUE activates the full linearised operator (Newton).
  ! FALSE activates a partially linearised operator without the Newton part.
  logical, intent(in) :: bfull
!</input>

!<inputoutput>
  ! Structure with temporary assembly data.
  type(t_assemblyTempDataSpace), intent(inout) :: rtempData
    
  ! Destination vector
  type(t_vectorBlock), intent(inout) :: rdest
!</inputoutput>

!</subroutine>

    ! local variables
    real(DP) :: dtheta, dtstep, dtimeend, dtimestart
    type(t_vectorBlock), pointer :: p_rvector
    
    ! Clear the output vector
    call lsysbl_clearVector (rdest)

    ! Timestepping technique?
    select case (rspatialOperatorAsm%p_rtimeDiscr%ctype)
    
    ! ***********************************************************
    ! Standard Theta one-step scheme.
    ! ***********************************************************
    case (TDISCR_ONESTEPTHETA)
    
      ! Theta-scheme identifier
      dtheta = rspatialOperatorAsm%p_rtimeDiscr%dtheta

      ! Characteristics of the current timestep.
      call tdiscr_getTimestep(rspatialOperatorAsm%p_rtimeDiscr,idofTime-1,&
          dtimeend,dtstep,dtimestart)

      ! itag=0: old 1-step scheme.
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      select case (rspatialOperatorAsm%p_rtimeDiscr%itag)
      
      ! ***********************************************************
      ! itag=0: old/standard 1-step scheme.
      ! ***********************************************************
      case (0)

        call output_line("Old 1-step-scheme not implemented",&
            OU_CLASS_ERROR,OU_MODE_STD,"smva_getDef_dualLin")
        call sys_halt()

      ! ***********************************************************
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      ! ***********************************************************
      case (1)

        ! Which equation do we have?    
        select case (rspatialOperatorAsm%p_rphysics%cequation)

        ! *************************************************************
        ! Stokes/Navier Stokes.
        ! *************************************************************
        case (0,1)

          ! ***********************************************
          ! NEXT TIMESTEP
          ! ***********************************************
          if (idofTime .lt. rprimalSol%p_rvector%NEQtime) then

            ! ===============================================
            ! Prepare the linear parts of the matrix.
            ! ===============================================

            ! Clear the temporary matrix
            call lsysbl_clearMatrix (rtempData%rmatrix)
            
            ! -----------------------------------------
            ! Mass matrix for timestepping
            call smva_getMassMatrix (rspatialOperatorAsm,rtempData%rmatrix,1.0_DP-dtstep)
            
            ! -----------------------------------------
            ! Laplace -- if the viscosity is constant
            if (rspatialOperatorAsm%p_rphysics%cviscoModel .eq. 0) then
              call smva_getLaplaceMatrix (&
                  rspatialOperatorAsm,rtempData%rmatrix,rspatialOperatorAsm%p_rphysics%dnuConst)
            else
              call output_line("Nonconstant viscosity not supported.",&
                  OU_CLASS_ERROR,OU_MODE_STD,"smva_getDef_dualLin")
              call sys_halt()
            end if

            ! -----------------------------------------
            ! Realise the defect
            call sptivec_getVectorFromPool (rdualLinSol%p_rvectorAccess,idofTime-1,p_rvector)
            call lsysbl_blockMatVec (rtempData%rmatrix, p_rvector, rdest, -1.0_DP, 1.0_DP)
          
          end if
            
          ! ***********************************************
          ! CURRENT TIMESTEP
          ! ***********************************************

          ! ===============================================
          ! RHS assembly
          ! ===============================================
               
          call smva_getRhs_dualLin (rspatialOperatorAsm,idofTime,rdualSol,rprimalLinSol,&
              bfull,1.0_DP,rdest)
            
          ! ===============================================
          ! Prepare the linear parts of the matrix.
          ! ===============================================

          ! Clear the temporary matrix
          call lsysbl_clearMatrix (rtempData%rmatrix)
          
          ! -----------------------------------------
          ! Mass matrix for timestepping
          call smva_getMassMatrix (rspatialOperatorAsm,rtempData%rmatrix,dtstep)
          
          ! -----------------------------------------
          ! Laplace -- if the viscosity is constant
          if (rspatialOperatorAsm%p_rphysics%cviscoModel .eq. 0) then
            call smva_getLaplaceMatrix (&
                rspatialOperatorAsm,rtempData%rmatrix,rspatialOperatorAsm%p_rphysics%dnuConst)
          else
            call output_line("Nonconstant viscosity not supported.",&
                OU_CLASS_ERROR,OU_MODE_STD,"smva_getDef_dualLin")
            call sys_halt()
          end if
          
          ! -----------------------------------------
          ! B-matrices
          call smva_getBMatrix (rspatialOperatorAsm,rtempData%rmatrix,1.0_DP)
          
          ! -----------------------------------------
          ! D-matrices
          call smva_getDMatrix (rspatialOperatorAsm,rtempData%rmatrix,1.0_DP)

          ! -----------------------------------------
          ! EOJ-stabilisation
          !if (rspatialOperatorAsm%p_rsettingsDiscr%rstabilConvecPrimal%cupwind .eq. 4) then
          !  call smva_getEOJMatrix (rspatialOperatorAsm,rtempData%rmatrix,1.0_DP)
          !end if
            
          ! ===============================================
          ! Prepare the semilinear parts of the matrix.
          ! ===============================================
          
          ! The semilinear parts of the matrix can be set up with
          ! the block matrix assembly routines. Invoke them using
          ! a separate subroutine
          call smva_getSemilinMat_Dual (rspatialOperatorAsm,idofTime,&
              rprimalSol,1.0_DP,rtempData%rmatrix)

          ! -----------------------------------------
          ! Realise the defect with the linear and
          ! semilinear parts of the operator
          ! -----------------------------------------

          call sptivec_getVectorFromPool (rdualLinSol%p_rvectorAccess,idofTime,p_rvector)
          call lsysbl_blockMatVec (rtempData%rmatrix, p_rvector, rdest, -1.0_DP, 1.0_DP)
          
          ! ***********************************************
          ! ADDITIONAL RHS TERMS
          ! ***********************************************
          
          ! For the linearised backward equation, there are some additional
          ! terms to be subtracted, which stem from the full Frechet derivative
          ! of the dual equation. There is:
          !
          !    (dual operator) = rhs - (rdualSol, (rprimalLinSol grad) (phi))
          !                          - (rdualSol,  grad(rprimalLinSol) phi)
          !
          ! For the assembly, the RHS assembly has to be invoked.
          call smva_getSemilinRhs_dualLin (rspatialOperatorAsm,idofTime,&
              rdualSol,rprimalLinSol,bfull,-1.0_DP,rdest)
          
        end select ! Equation
      
      end select ! Timestep sub-scheme
      
    end select ! Timestep scheme

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>

  subroutine smva_getMassMatrix (rspatialOperatorAsm,rmatrix,dweight)

!<description>
  ! Implements the mass matrix into rmatrix, weighted by dweight
!</description>

!<input>
  ! Definition of the operator. Must be a linearised operator.
  type(t_spacetimeOperatorAsm) :: rspatialOperatorAsm
  
  ! Weight for the operator
  real(DP), intent(in) :: dweight
!</input>

!<inputoutput>
  ! Block matrix receiving the result.
  type(t_matrixBlock), intent(inout) :: rmatrix
!</inputoutput>

!</subroutine>
    
    ! Which equation do we have?    
    select case (rspatialOperatorAsm%p_rphysics%cequation)

    ! ***********************************************************
    ! Stokes/Navier Stokes.
    ! ***********************************************************
    case (0,1)
    
      call lsyssc_matrixLinearComb (&
          rspatialOperatorAsm%p_rasmTemplates%rmatrixMassVelocity,dweight,&
          rmatrix%RmatrixBlock(1,1),1.0_DP,rmatrix%RmatrixBlock(1,1),&
          .false.,.false.,.true.,.true.)

      call lsyssc_matrixLinearComb (&
          rspatialOperatorAsm%p_rasmTemplates%rmatrixMassVelocity,dweight,&
          rmatrix%RmatrixBlock(2,2),1.0_DP,rmatrix%RmatrixBlock(2,2),&
          .false.,.false.,.true.,.true.)
    end select
    
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine smva_getLaplaceMatrix (rspatialOperatorAsm,rmatrix,dweight)

!<description>
  ! Implements the Laplace matrix into rmatrix, weighted by dweight
!</description>

!<input>
  ! Definition of the operator. Must be a linearised operator.
  type(t_spacetimeOperatorAsm) :: rspatialOperatorAsm
  
  ! Weight for the operator
  real(DP), intent(in) :: dweight
!</input>

!<inputoutput>
  ! Block matrix receiving the result.
  type(t_matrixBlock), intent(inout) :: rmatrix
!</inputoutput>

!</subroutine>
    
    ! Which equation do we have?    
    select case (rspatialOperatorAsm%p_rphysics%cequation)

    ! ***********************************************************
    ! Stokes/Navier Stokes.
    ! ***********************************************************
    case (0,1)
    
      call lsyssc_matrixLinearComb (&
          rspatialOperatorAsm%p_rasmTemplates%rmatrixLaplace,dweight,&
          rmatrix%RmatrixBlock(1,1),1.0_DP,rmatrix%RmatrixBlock(1,1),&
          .false.,.false.,.true.,.true.)

      call lsyssc_matrixLinearComb (&
          rspatialOperatorAsm%p_rasmTemplates%rmatrixLaplace,dweight,&
          rmatrix%RmatrixBlock(2,2),1.0_DP,rmatrix%RmatrixBlock(2,2),&
          .false.,.false.,.true.,.true.)
    end select
    
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine smva_getBMatrix (rspatialOperatorAsm,rmatrix,dweight)

!<description>
  ! Implements the B-matrix into rmatrix, weighted by dweight
!</description>

!<input>
  ! Definition of the operator. Must be a linearised operator.
  type(t_spacetimeOperatorAsm) :: rspatialOperatorAsm
  
  ! Weight for the operator
  real(DP), intent(in) :: dweight
!</input>

!<inputoutput>
  ! Block matrix receiving the result.
  type(t_matrixBlock), intent(inout) :: rmatrix
!</inputoutput>

!</subroutine>
    
    ! Which equation do we have?    
    select case (rspatialOperatorAsm%p_rphysics%cequation)

    ! ***********************************************************
    ! Stokes/Navier Stokes.
    ! ***********************************************************
    case (0,1)
      call lsyssc_matrixLinearComb (&
          rspatialOperatorAsm%p_rasmTemplates%rmatrixB1,dweight,&
          rmatrix%RmatrixBlock(1,3),1.0_DP,rmatrix%RmatrixBlock(1,3),&
          .false.,.false.,.true.,.true.)

      call lsyssc_matrixLinearComb (&
          rspatialOperatorAsm%p_rasmTemplates%rmatrixB2,dweight,&
          rmatrix%RmatrixBlock(2,3),1.0_DP,rmatrix%RmatrixBlock(2,3),&
          .false.,.false.,.true.,.true.)
          
    end select
    
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine smva_getDMatrix (rspatialOperatorAsm,rmatrix,dweight)

!<description>
  ! Implements the D-matrix into rmatrix, weighted by dweight
!</description>

!<input>
  ! Definition of the operator. Must be a linearised operator.
  type(t_spacetimeOperatorAsm) :: rspatialOperatorAsm
  
  ! Weight for the operator
  real(DP), intent(in) :: dweight
!</input>

!<inputoutput>
  ! Block matrix receiving the result.
  type(t_matrixBlock), intent(inout) :: rmatrix
!</inputoutput>

!</subroutine>
    
    ! Which equation do we have?    
    select case (rspatialOperatorAsm%p_rphysics%cequation)
    
    ! ***********************************************************
    ! Stokes/Navier Stokes.
    ! ***********************************************************
    case (0,1)
      call lsyssc_matrixLinearComb (&
          rspatialOperatorAsm%p_rasmTemplates%rmatrixD1,dweight,&
          rmatrix%RmatrixBlock(3,1),1.0_DP,rmatrix%RmatrixBlock(3,1),&
          .false.,.false.,.true.,.true.)

      call lsyssc_matrixLinearComb (&
          rspatialOperatorAsm%p_rasmTemplates%rmatrixD2,dweight,&
          rmatrix%RmatrixBlock(3,2),1.0_DP,rmatrix%RmatrixBlock(3,2),&
          .false.,.false.,.true.,.true.)
    end select
    
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine smva_getEOJMatrix (rspatialOperatorAsm,rmatrix,dweight)

!<description>
  ! Implements the D-matrix into rmatrix, weighted by dweight
!</description>

!<input>
  ! Definition of the operator. Must be a linearised operator.
  type(t_spacetimeOperatorAsm) :: rspatialOperatorAsm
  
  ! Weight for the operator
  real(DP), intent(in) :: dweight
!</input>

!<inputoutput>
  ! Block matrix receiving the result.
  type(t_matrixBlock), intent(inout) :: rmatrix
!</inputoutput>

!</subroutine>

    ! not yet implemented    
    
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine smva_getSemilinMat_primal (rspatialOperatorAsm,idofTime,&
      rprimalSol,dweight,rmatrix)

!<description>
  ! Implements the semilinear part of the operator described by
  ! rspatialOperatorAsm into rmatrix, weighted by dweight.
  !
  ! Nonlinear backward equation.
!</description>

!<input>
  ! Definition of the operator. Must be a linearised operator.
  type(t_spacetimeOperatorAsm), target :: rspatialOperatorAsm

  ! Number of the DOF in time which should be calculated into rmatrix.
  integer, intent(in) :: idofTime
  
  ! Weight for the operator
  real(DP), intent(in) :: dweight

  ! Space-time vector which describes the nonlinearity in the operator.
  ! This is the primal solution of the last nonlinear iteration.
  type(t_primalSpace), intent(inout) :: rprimalSol
!</input>

!<inputoutput>
  ! Block matrix receiving the result.
  type(t_matrixBlock), intent(inout) :: rmatrix
!</inputoutput>

!</subroutine>

    ! local variables    
    type(p_t_spacetimeOperatorAsm), pointer :: p_p_rspatialOperatorAsm
    type(t_collection) :: rcollection
    type(t_fev2Vectors) :: rvectorEval
    type(t_vectorBlock), pointer :: p_rvector

    ! Cancel if nothing to do
    if (dweight .eq. 0.0_DP) return

    ! Prepare a pointer to our operator for callback routines.
    ! It is passed via rcollection.
    p_p_rspatialOperatorAsm%p_rspatialOperatorAsm => rspatialOperatorAsm
    rcollection%IquickAccess(8:) = transfer(p_p_rspatialOperatorAsm,rcollection%IquickAccess(8:))
    rcollection%DquickAccess(1) = dweight

    ! Which equation do we have?
    select case (rspatialOperatorAsm%p_rphysics%cequation)

    ! ***********************************************************
    ! Navier Stokes
    ! ***********************************************************
    case (1)
    
      ! Notify the callback routine what to assemble.
      rcollection%IquickAccess(1) = OPTP_PRIMAL
      
      ! Prepare the evaluation
      !
      ! Vector 1+2 = primal velocity (for the nonlinearity)
      call sptivec_getVectorFromPool (rprimalSol%p_rvectorAccess,idofTime,p_rvector)
      call fev2_addVectorToEvalList(rvectorEval,p_rvector%RvectorBlock(1),0)
      call fev2_addVectorToEvalList(rvectorEval,p_rvector%RvectorBlock(2),0)

      ! Build the matrix      
      call bma_buildMatrix (rmatrix,BMA_CALC_STANDARD,&
          smva_fcalc_semilinearMat, rcollection, revalVectors=rvectorEval,&
          rcubatureInfo=rspatialOperatorAsm%p_rasmTemplates%rcubatureInfoMassVelocity)

      ! Cleanup
      call fev2_releaseVectorList(rvectorEval)

    end select
    
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine smva_getSemilinMat_primalLin (rspatialOperatorAsm,idofTime,&
      rprimalSol,dweight,bfull,rmatrix)

!<description>
  ! Implements the semilinear part of the operator described by
  ! rspatialOperatorAsm into rmatrix, weighted by dweight.
  !
  ! Linearised forward equation.
!</description>

!<input>
  ! Definition of the operator. Must be a linearised operator.
  type(t_spacetimeOperatorAsm), target :: rspatialOperatorAsm
  
  ! Number of the DOF in time which should be calculated into rmatrix.
  integer, intent(in) :: idofTime
  
  ! Weight for the operator
  real(DP), intent(in) :: dweight

  ! Space-time vector which describes the nonlinearity in the operator.
  ! This is the primal solution of the last nonlinear iteration.
  type(t_primalSpace), intent(inout) :: rprimalSol

  ! TRUE activates the full linearised operator (Newton).
  ! FALSE activates a partially linearised operator without the Newton part.
  logical, intent(in) :: bfull
!</input>

!<inputoutput>
  ! Block matrix receiving the result.
  type(t_matrixBlock), intent(inout) :: rmatrix
!</inputoutput>

!</subroutine>

    ! local variables    
    type(p_t_spacetimeOperatorAsm), pointer :: p_p_rspatialOperatorAsm
    type(t_collection) :: rcollection
    type(t_fev2Vectors) :: rvectorEval
    type(t_vectorBlock), pointer :: p_rvector

    ! Cancel if nothing to do
    if (dweight .eq. 0.0_DP) return

    ! Prepare a pointer to our operator for callback routines.
    ! It is passed via rcollection.
    p_p_rspatialOperatorAsm%p_rspatialOperatorAsm => rspatialOperatorAsm
    rcollection%IquickAccess(8:) = transfer(p_p_rspatialOperatorAsm,rcollection%IquickAccess(8:))
    rcollection%DquickAccess(1) = dweight

    ! Which equation do we have?
    select case (rspatialOperatorAsm%p_rphysics%cequation)

    ! ***********************************************************
    ! Navier Stokes
    ! ***********************************************************
    case (1)
    
      ! Notify the callback routine what to assemble.
      if (bfull) then
        rcollection%IquickAccess(1) = OPTP_PRIMALLIN
      else
        rcollection%IquickAccess(1) = OPTP_PRIMALLIN_SIMPLE
      end if
      
      ! Prepare the evaluation
      !
      ! Vector 1+2 = primal velocity.
      call sptivec_getVectorFromPool (rprimalSol%p_rvectorAccess,idofTime,p_rvector)
      call fev2_addVectorToEvalList(rvectorEval,p_rvector%RvectorBlock(1),0)
      call fev2_addVectorToEvalList(rvectorEval,p_rvector%RvectorBlock(2),0)
      
      ! Build the matrix
      call bma_buildMatrix (rmatrix,BMA_CALC_STANDARD,&
          smva_fcalc_semilinearMat, rcollection, revalVectors=rvectorEval,&
          rcubatureInfo=rspatialOperatorAsm%p_rasmTemplates%rcubatureInfoMassVelocity)
          
      ! Cleanup
      call fev2_releaseVectorList(rvectorEval)
          
    end select
    
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine smva_getSemilinMat_Dual (rspatialOperatorAsm,idofTime,rprimalSol,dweight,rmatrix)

!<description>
  ! Implements the semilinear part of the operator described by
  ! rspatialOperatorAsm into rmatrix, weighted by dweight.
  !
  ! Backward equation, nonlinear as well as linearised version.
!</description>

!<input>
  ! Definition of the operator. Must be a linearised operator.
  type(t_spacetimeOperatorAsm), target :: rspatialOperatorAsm

  ! Number of the DOF in time which should be calculated into rmatrix.
  integer, intent(in) :: idofTime
  
  ! Weight for the operator
  real(DP), intent(in) :: dweight

  ! Space-time vector which describes the nonlinearity in the operator.
  ! This is the primal solution of the last nonlinear iteration.
  type(t_primalSpace), intent(inout) :: rprimalSol
!</input>

!<inputoutput>
  ! Block matrix receiving the result.
  type(t_matrixBlock), intent(inout) :: rmatrix
!</inputoutput>

!</subroutine>

    ! local variables    
    type(p_t_spacetimeOperatorAsm), pointer :: p_p_rspatialOperatorAsm
    type(t_collection) :: rcollection
    type(t_fev2Vectors) :: rvectorEval
    type(t_vectorBlock), pointer :: p_rvector

    ! Cancel if nothing to do
    if (dweight .eq. 0.0_DP) return

    ! Prepare a pointer to our operator for callback routines.
    ! It is passed via rcollection.
    p_p_rspatialOperatorAsm%p_rspatialOperatorAsm => rspatialOperatorAsm
    rcollection%IquickAccess(8:) = transfer(p_p_rspatialOperatorAsm,rcollection%IquickAccess(8:))
    rcollection%DquickAccess(1) = dweight

    ! Notify the callback routine what to assemble.
    rcollection%IquickAccess(1) = OPTP_DUAL
          
    ! Which equation do we have?
    select case (rspatialOperatorAsm%p_rphysics%cequation)

    ! ***********************************************************
    ! Navier Stokes
    ! ***********************************************************
    case (1)
    
      ! Prepare the evaluation.
      !
      ! Vector 1+2 = primal velocity, including 1st derivative.
      call sptivec_getVectorFromPool (rprimalSol%p_rvectorAccess,idofTime,p_rvector)
      call fev2_addVectorToEvalList(rvectorEval,p_rvector%RvectorBlock(1),1)
      call fev2_addVectorToEvalList(rvectorEval,p_rvector%RvectorBlock(2),1)
      
      ! Build the matrix
      call bma_buildMatrix (rmatrix,BMA_CALC_STANDARD,&
          smva_fcalc_semilinearMat, rcollection, revalVectors=rvectorEval,&
          rcubatureInfo=rspatialOperatorAsm%p_rasmTemplates%rcubatureInfoMassVelocity)
          
      ! Cleanup
      call fev2_releaseVectorList(rvectorEval)
          
    end select
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine smva_fcalc_semilinearMat(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the semilinear part of the RHS.
    ! Calculates the local vectors Dentry in all rvectorData structures.
!</description>

!<inputoutput>
    ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
    ! have to be filled with data.
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</input>

!<subroutine>

    ! Local variables
    real(DP) :: dweight
    real(DP) :: dbasI, dbasJ, dbasIx, dbasJx, dbasIy, dbasJy
    real(DP) :: du1, du2, du1x, du1y, du2x, du2y
    integer :: iel, icubp, idofe, jdofe, i
    integer :: coptype
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix11
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix12
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix21
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix22
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial,p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaMatrixData), pointer :: p_rmatrixData11,p_rmatrixData22
    type(t_bmaMatrixData), pointer :: p_rmatrixData12,p_rmatrixData21
    real(DP), dimension(:,:,:), pointer :: p_DnonlinearityY1 => null()
    real(DP), dimension(:,:,:), pointer :: p_DnonlinearityY2 => null()

    ! local variables    
    type(p_t_spacetimeOperatorAsm), pointer :: p_p_rspatialOperatorAsm
    type(t_spacetimeOperatorAsm), pointer :: p_rspatialOperatorAsm

    ! From the collection, fetch our operator structure, nonlinearity,...
    p_p_rspatialOperatorAsm = transfer(rcollection%IquickAccess(8:),p_p_rspatialOperatorAsm)
    p_rspatialOperatorAsm => p_p_rspatialOperatorAsm%p_rspatialOperatorAsm
    
    ! Type of the operator to compute.
    copType = rcollection%IquickAccess(1)

    ! Get the nonlinearity Y -- X-velocity and Y-velocity.
    p_DnonlinearityY1 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,:)
    p_DnonlinearityY2 => revalVectors%p_RvectorData(2)%p_Ddata(:,:,:)

    ! Weight of the operator.
    dweight = rcollection%DquickAccess(1)

    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Get local data of the velocity submatrices.
    ! They are all discretised with the same FEM space.
    p_rmatrixData11 => RmatrixData(1,1)
    p_rmatrixData21 => RmatrixData(2,1)
    p_rmatrixData12 => RmatrixData(1,2)
    p_rmatrixData22 => RmatrixData(2,2)
    p_DbasTrial => RmatrixData(i,i)%p_DbasTrial
    p_DbasTest => RmatrixData(i,i)%p_DbasTest

    ! Get the matrix data      
    p_DlocalMatrix11 => RmatrixData(1,1)%p_Dentry
    p_DlocalMatrix12 => RmatrixData(1,2)%p_Dentry
    p_DlocalMatrix21 => RmatrixData(2,1)%p_Dentry
    p_DlocalMatrix22 => RmatrixData(2,2)%p_Dentry

    ! Which equation do we have?
    select case (p_rspatialOperatorAsm%p_rphysics%cequation)
    case (1)
    
      ! -------------------------------------------------------------
      ! Navier Stokes.
      ! -------------------------------------------------------------
      
      ! Primal or dual equation?
      select case (copType)

      ! ***********************************************************
      ! Navier-Stokes. Forward equation.
      ! ***********************************************************
      case (OPTP_PRIMAL,OPTP_PRIMALLIN_SIMPLE)
      
        ! In case of the Stokes/Navier-Stokes equations, we can
        ! save some arithmetic operations in some situations.
        ! Only for the forward equation, A11 may coincide with
        ! A22, so we only have to calculate this matrix once:
        ! Either in the A11 block or in A11 and A22, depending on
        ! whether A11 and A22 share their data.
        !
        ! In all other situationsm A11 does not coincide
        ! with A22, so we do not have to respect that possibility.
        
        if (p_rmatrixData22%bsharedMatrixData) then
        
          ! ---------------------------------------------------------
          ! Assemble the nonlinearity "((y grad) (phi_j) , phi_i)"
          do iel = 1,nelements
            do icubp = 1,npointsPerElement
            
              ! Get the X-velocity and the Y-velocity in that point
              ! du1 = y_1, du2 = y_2
              du1 = p_DnonlinearityY1(icubp,iel,DER_FUNC)
              du2 = p_DnonlinearityY2(icubp,iel,DER_FUNC)
            
              do idofe=1,p_rmatrixData11%ndofTest

                dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

                do jdofe=1,p_rmatrixData11%ndofTrial

                  dbasJx = p_DbasTrial(jdofe,DER_DERIV2D_X,icubp,iel)
                  dbasJy = p_DbasTrial(jdofe,DER_DERIV2D_Y,icubp,iel)

                  p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                      dweight * p_DcubWeight(icubp,iel) * (dbasJx*du1+dbasJy*du2)*dbasI

                end do ! idofe
              end do ! jdofe
            end do ! icubp
          end do ! iel
          
        else

          ! ---------------------------------------------------------
          ! Assemble the nonlinearity "((y grad) (phi_j) , phi_i)"
          do iel = 1,nelements
            do icubp = 1,npointsPerElement

              ! Get the X-velocity and the Y-velocity in that point:
              ! du1 = y_1, du2 = y_2
              du1 = p_DnonlinearityY1(icubp,iel,DER_FUNC)
              du2 = p_DnonlinearityY2(icubp,iel,DER_FUNC)

              do idofe=1,p_rmatrixData11%ndofTest

                dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

                do jdofe=1,p_rmatrixData11%ndofTrial

                  dbasJx = p_DbasTrial(jdofe,DER_DERIV2D_X,icubp,iel)
                  dbasJy = p_DbasTrial(jdofe,DER_DERIV2D_Y,icubp,iel)

                  p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                      dweight * p_DcubWeight(icubp,iel) * (dbasJx*du1+dbasJy*du2)*dbasI

                  p_DlocalMatrix22(jdofe,idofe,iel) = p_DlocalMatrix22(jdofe,idofe,iel) + &
                      dweight * p_DcubWeight(icubp,iel) * (dbasJx*du1+dbasJy*du2)*dbasI

                end do ! idofe
              end do ! jdofe
            end do ! icubp
          end do ! iel
        
        end if
        
      ! ***********************************************************
      ! Linearised Navier-Stokes. Forward equation.
      ! ***********************************************************
      case (OPTP_PRIMALLIN)

        ! ---------------------------------------------------------
        ! The nonlinearity consists of two terms:
        ! a) "((y grad) (phi_j) , phi_i)"  
        ! b) "(grad(y) phi_j , phi_i)"
        do iel = 1,nelements
          do icubp = 1,npointsPerElement

            ! Get the X-velocity and the Y-velocity in that point:
            ! du1 = y_1, du2 = y_2
            du1 = p_DnonlinearityY1(icubp,iel,DER_FUNC)
            du2 = p_DnonlinearityY2(icubp,iel,DER_FUNC)

            ! as well as their derivatives
            du1x = p_DnonlinearityY1(icubp,iel,DER_DERIV2D_X)
            du1y = p_DnonlinearityY1(icubp,iel,DER_DERIV2D_Y)
            du2x = p_DnonlinearityY2(icubp,iel,DER_DERIV2D_X)
            du2y = p_DnonlinearityY2(icubp,iel,DER_DERIV2D_Y)

            do idofe=1,p_rmatrixData11%ndofTest

              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              do jdofe=1,p_rmatrixData11%ndofTrial

                dbasJ  = p_DbasTrial(jdofe,DER_FUNC,icubp,iel)
                dbasJx = p_DbasTrial(jdofe,DER_DERIV2D_X,icubp,iel)
                dbasJy = p_DbasTrial(jdofe,DER_DERIV2D_Y,icubp,iel)

                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dweight * p_DcubWeight(icubp,iel) * &
                    ( (dbasJx*du1+dbasJy*du2)*dbasI + &    ! "((y grad) (phi_j) , phi_i)"
                      du1x * dbasJ * dbasI )               !  "( grad(y) phi_j , phi_i)"

                p_DlocalMatrix12(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dweight * p_DcubWeight(icubp,iel) * &
                    ( du1y * dbasJ * dbasI )               !  "( grad(y) phi_j , phi_i)"

                p_DlocalMatrix21(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dweight * p_DcubWeight(icubp,iel) * &
                    ( du2x * dbasJ * dbasI )               !  "( grad(y) phi_j , phi_i)"

                p_DlocalMatrix22(jdofe,idofe,iel) = p_DlocalMatrix22(jdofe,idofe,iel) + &
                    dweight * p_DcubWeight(icubp,iel) * &
                    ( ( dbasJx*du1+dbasJy*du2)*dbasI + &   ! "((y grad) (phi_j) , phi_i)"
                      du2y * dbasJ * dbasI )               !  "( grad(y) phi_j , phi_i)"

              end do ! idofe
            end do ! jdofe
          end do ! icubp
        end do ! iel
        
      ! ***********************************************************
      ! Navier-Stokes. Backward equation.
      ! ***********************************************************
      case (OPTP_DUAL)
        
        if (p_rmatrixData22%bsharedMatrixData) then
          call output_line("Shared matrix data not supported.",&
              OU_CLASS_ERROR,OU_MODE_STD,"smva_fcalc_semilinear")
          call sys_halt()
        end if
        
        ! ---------------------------------------------------------
        ! Assemble the nonlinearity. There are two terms:
        ! a)  "( phi_j , (y grad) phi_i)"
        ! b)  "( phi_j ,  grad(y) phi_i)"
        ! The dual equation is always linear, so there is de-facto
        ! only one type of equation. However, whether the dual equation
        ! is used for the Newton iteration or not, the right-hand side
        ! changes.
        do iel = 1,nelements
          do icubp = 1,npointsPerElement

            ! Get the X-velocity and the Y-velocity in that point:
            ! du1 = y_1, du2 = y_2
            du1 = p_DnonlinearityY1(icubp,iel,DER_FUNC)
            du2 = p_DnonlinearityY2(icubp,iel,DER_FUNC)

            ! as well as their derivatives
            du1x = p_DnonlinearityY1(icubp,iel,DER_DERIV2D_X)
            du1y = p_DnonlinearityY1(icubp,iel,DER_DERIV2D_Y)
            du2x = p_DnonlinearityY2(icubp,iel,DER_DERIV2D_X)
            du2y = p_DnonlinearityY2(icubp,iel,DER_DERIV2D_Y)

            do idofe=1,p_rmatrixData11%ndofTest

              dbasIx = p_DbasTest(idofe,DER_DERIV2D_X,icubp,iel)
              dbasIy = p_DbasTest(idofe,DER_DERIV2D_Y,icubp,iel)

              do jdofe=1,p_rmatrixData11%ndofTrial

                dbasJ = p_DbasTrial(jdofe,DER_FUNC,icubp,iel)

                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dweight * p_DcubWeight(icubp,iel) * &
                      ( dbasJ * (du1*dbasIx+du2*dbasIy) + &  !  "( phi_j , (y grad) phi_i)"
                        dbasJ * du1x * dbasI )               !  "( phi_j ,  grad(y) phi_i)"

                p_DlocalMatrix12(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dweight * p_DcubWeight(icubp,iel) * &
                      ( dbasJ * du1y * dbasI )               !  "( phi_j ,  grad(y) phi_i)"

                p_DlocalMatrix12(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dweight * p_DcubWeight(icubp,iel) * &
                      ( dbasJ * du2x * dbasI )               !  "( phi_j ,  grad(y) phi_i)"

                p_DlocalMatrix22(jdofe,idofe,iel) = p_DlocalMatrix22(jdofe,idofe,iel) + &
                    dweight * p_DcubWeight(icubp,iel) * &
                      ( dbasJ * (du1*dbasIx+du2*dbasIy) + &  !  "( phi_j , (y grad) phi_i)"
                        dbasJ * du2y * dbasI )               !  "( phi_j ,  grad(y) phi_i)"

              end do ! idofe
            end do ! jdofe
          end do ! icubp
        end do ! iel
        
      end select ! Operator type
      
    end select ! Equation

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine smva_getRhs_Primal (rspatialOperatorAsm,idofTime,rcontrol,dweight,rrhs)

!<description>
  ! Calculates the RHS of the primal equation, based on a 'current'
  ! dual solution.
!</description>

!<input>
  ! Definition of the operator. Must be a linearised operator.
  type(t_spacetimeOperatorAsm), target :: rspatialOperatorAsm

  ! Number of the DOF in time which should be calculated into rrhs.
  integer, intent(in) :: idofTime
  
  ! Weight for the operator.
  real(DP), intent(in) :: dweight

  ! Space-time vector which contains the control.
  type(t_controlSpace), intent(inout) :: rcontrol
!</input>

!<inputoutput>
  ! Block vector receiving the result.
  type(t_vectorBlock), intent(inout) :: rrhs
!</inputoutput>

!</subroutine>

    ! local variables    
    type(p_t_spacetimeOperatorAsm), pointer :: p_p_rspatialOperatorAsm
    type(t_collection) :: rcollection
    type(t_collection), target :: ruserCollection
    type(t_fev2Vectors) :: rvectorEval
    type(t_vectorBlock), pointer :: p_rvector1
    real(DP) :: dtheta, dtstep, dtime, dtimeend, dtimestart

    ! Cancel if nothing to do
    if (dweight .eq. 0.0_DP) return

    ! Prepare a local and a user-defined collection
    call collct_init (rcollection)
    call collct_init (ruserCollection)
    rcollection%p_rnextCollection => ruserCollection
    
    ! Prepare a pointer to our operator for callback routines.
    ! It is passed via rcollection.
    p_p_rspatialOperatorAsm%p_rspatialOperatorAsm => rspatialOperatorAsm
    rcollection%IquickAccess(8:) = transfer(p_p_rspatialOperatorAsm,rcollection%IquickAccess(8:))
    rcollection%DquickAccess(1) = dweight
    
    ! Notify the callback routine what to assemble.
    rcollection%IquickAccess(1) = OPTP_PRIMAL
    
    ! Timestepping technique?
    select case (rspatialOperatorAsm%p_rtimeDiscr%ctype)
    
    ! ***********************************************************
    ! Standard Theta one-step scheme.
    ! ***********************************************************
    case (TDISCR_ONESTEPTHETA)
    
      ! Theta-scheme identifier
      dtheta = rspatialOperatorAsm%p_rtimeDiscr%dtheta
      
      ! Characteristics of the current timestep.
      call tdiscr_getTimestep(rspatialOperatorAsm%p_rtimeDiscr,idofTime-1,&
          dtimeend,dtstep,dtimestart)

      ! itag=0: old 1-step scheme.
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      select case (rspatialOperatorAsm%p_rtimeDiscr%itag)
      
      ! ***********************************************************
      ! itag=0: old/standard 1-step scheme.
      ! ***********************************************************
      case (0)

        call output_line("Old 1-step-scheme not implemented",&
            OU_CLASS_ERROR,OU_MODE_STD,"smva_getRhs_Primal")
        call sys_halt()

      ! ***********************************************************
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      ! ***********************************************************
      case (1)

        ! Which equation do we have?
        select case (rspatialOperatorAsm%p_rphysics%cequation)

        ! ***********************************************************
        ! Stokes/Navier Stokes
        ! ***********************************************************
        case (0,1)
      
          ! ---------------------------------------------------------
          ! The one-and-only RHS
          ! ---------------------------------------------------------

          ! Evaluation point of the RHS in time
          dtime = dtimestart + (1.0_DP-dtheta) * dtstep
          
          ! Prepare the user-defined collection for the assembly
          call user_initCollectForAssembly (rspatialOperatorAsm%p_rglobalData,dtime,rusercollection)
          
          ! Prepare the evaluation of the primal RHS.
          call ansol_prepareEval (rspatialOperatorAsm%p_rrhsPrimal,rcollection,"RHS",dtime)

          ! Prepare the evaluation.
          !
          ! Vector 1+2 = Temp-vectors for the RHS.
          call fev2_addDummyVectorToEvalList(rvectorEval)
          call fev2_addDummyVectorToEvalList(rvectorEval)
          
          ! Vector 3+4 = dual velocity -- for the calculation of distributed control
          ! in the velocity space.
          ! Only add this if we have distributed control. Otherwise add dummy
          ! vectors which take no time in being computed.
          if (rspatialOperatorAsm%p_rsettingsOptControl%dalphaC .ge. 0.0_DP) then
            call sptivec_getVectorFromPool (rcontrol%p_rvectorAccess,idofTime,p_rvector1)
            call fev2_addVectorToEvalList(rvectorEval,p_rvector1%RvectorBlock(1),0)
            call fev2_addVectorToEvalList(rvectorEval,p_rvector1%RvectorBlock(2),0)
          else
            call fev2_addDummyVectorToEvalList(rvectorEval)
            call fev2_addDummyVectorToEvalList(rvectorEval)
          end if
          
          ! Build the vector
          call bma_buildVector (rrhs,BMA_CALC_STANDARD,&
              smva_fcalc_rhs, rcollection, revalVectors=rvectorEval,&
              rcubatureInfo=rspatialOperatorAsm%p_rasmTemplates%rcubatureInfoRHScontinuity)
          
          ! Cleanup
          call fev2_releaseVectorList(rvectorEval)
          call ansol_doneEvalCollection (rcollection,"RHS")
          call user_doneCollectForAssembly (rspatialOperatorAsm%p_rglobalData,rusercollection)

        end select ! Equation

      end select ! Timestep sub-scheme
            
    end select ! Timestep scheme
    
    ! Release the collections
    call collct_done (ruserCollection)
    call collct_init (rcollection)
    
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine smva_getRhs_Dual (rspatialOperatorAsm,idofTime,rprimalSol,dweight,rrhs)

!<description>
  ! Calculates the RHS of the dual equation, based on a 'current'
  ! primal solution.
!</description>

!<input>
  ! Definition of the operator. Must be a linearised operator.
  type(t_spacetimeOperatorAsm), target :: rspatialOperatorAsm

  ! Number of the DOF in time which should be calculated into rrhs.
  integer, intent(in) :: idofTime
  
  ! Weight for the operator
  real(DP), intent(in) :: dweight

  ! Space-time vector which contains the solution of the dual equation.
  type(t_primalSpace), intent(inout) :: rprimalSol
!</input>

!<inputoutput>
  ! Block vector receiving the result.
  type(t_vectorBlock), intent(inout) :: rrhs
!</inputoutput>

!</subroutine>

    ! local variables    
    type(p_t_spacetimeOperatorAsm), pointer :: p_p_rspatialOperatorAsm
    type(t_fev2Vectors) :: rvectorEval
    type(t_vectorBlock), pointer :: p_rvector1
    type(t_collection) :: rcollection
    type(t_collection), target :: ruserCollection
    real(DP) :: dtheta, dtstep, dtime, dtimeend, dtimestart

    ! Cancel if nothing to do
    if (dweight .eq. 0.0_DP) return

    ! Prepare a local and a user-defined collection
    call collct_init (rcollection)
    call collct_init (ruserCollection)
    rcollection%p_rnextCollection => ruserCollection
    
    ! Prepare a pointer to our operator for callback routines.
    ! It is passed via rcollection.
    p_p_rspatialOperatorAsm%p_rspatialOperatorAsm => rspatialOperatorAsm
    rcollection%IquickAccess(8:) = transfer(p_p_rspatialOperatorAsm,rcollection%IquickAccess(8:))
    rcollection%DquickAccess(1) = dweight
    
    ! Notify the callback routine what to assemble.
    rcollection%IquickAccess(1) = OPTP_DUAL
    
    ! Timestepping technique?
    select case (rspatialOperatorAsm%p_rtimeDiscr%ctype)
    
    ! ***********************************************************
    ! Standard Theta one-step scheme.
    ! ***********************************************************
    case (TDISCR_ONESTEPTHETA)
    
      ! Theta-scheme identifier
      dtheta = rspatialOperatorAsm%p_rtimeDiscr%dtheta
      
      ! Characteristics of the current timestep.
      call tdiscr_getTimestep(rspatialOperatorAsm%p_rtimeDiscr,idofTime-1,&
          dtimeend,dtstep,dtimestart)

      ! itag=0: old 1-step scheme.
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      select case (rspatialOperatorAsm%p_rtimeDiscr%itag)
      
      ! ***********************************************************
      ! itag=0: old/standard 1-step scheme.
      ! ***********************************************************
      case (0)

        call output_line("Old 1-step-scheme not implemented",&
            OU_CLASS_ERROR,OU_MODE_STD,"smva_getRhs_Dual")
        call sys_halt()

      ! ***********************************************************
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      ! ***********************************************************
      case (1)

        ! Which equation do we have?
        select case (rspatialOperatorAsm%p_rphysics%cequation)

        ! ***********************************************************
        ! Stokes/Navier Stokes
        ! ***********************************************************
        case (0,1)
        
          ! Evaluation point of the RHS in time
          dtime = dtimestart
          
          ! Prepare the user-defined collection for the assembly
          call user_initCollectForAssembly (rspatialOperatorAsm%p_rglobalData,dtime,rusercollection)

          ! Prepare the evaluation of the primal RHS.
          call ansol_prepareEval (rspatialOperatorAsm%p_rrhsDual,rcollection,"RHS",dtime)
          call ansol_prepareEval (rspatialOperatorAsm%p_rtargetFlow,rcollection,"TARGET",dtime)
          
          ! Prepare the evaluation.
          !
          ! Vector 1+2 = Temp-vectors for the RHS.
          call fev2_addDummyVectorToEvalList(rvectorEval)
          call fev2_addDummyVectorToEvalList(rvectorEval)

          ! Vector 3+4 = Primal velocity
          call sptivec_getVectorFromPool (rprimalSol%p_rvectorAccess,idofTime,p_rvector1)
          call fev2_addVectorToEvalList(rvectorEval,p_rvector1%RvectorBlock(3),0)
          call fev2_addVectorToEvalList(rvectorEval,p_rvector1%RvectorBlock(4),0)

          ! Build the vector
          call bma_buildVector (rrhs,BMA_CALC_STANDARD,&
              smva_fcalc_rhs, rcollection, revalVectors=rvectorEval,&
              rcubatureInfo=rspatialOperatorAsm%p_rasmTemplates%rcubatureInfoRHScontinuity)
          
          ! Cleanup
          call fev2_releaseVectorList(rvectorEval)
          call ansol_doneEvalCollection (rcollection,"TARGER")
          call ansol_doneEvalCollection (rcollection,"RHS")
          call user_doneCollectForAssembly (rspatialOperatorAsm%p_rglobalData,rusercollection)
            
        end select ! Equation
        
      end select ! Timestep sub-scheme
      
    end select ! Timestep scheme
    
    ! Release the collections
    call collct_done (ruserCollection)
    call collct_init (rcollection)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine smva_getRhs_primalLin (rspatialOperatorAsm,idofTime,&
      rcontrol,rcontrolLin,dweight,rrhs)

!<description>
  ! Calculates the RHS of the linearised primal equation, based on a 'current'
  ! dual solution and dual linearised solution.
!</description>

!<input>
  ! Definition of the operator. Must be a linearised operator.
  type(t_spacetimeOperatorAsm), target :: rspatialOperatorAsm
  
  ! Number of the DOF in time which should be calculated into rrhs.
  integer, intent(in) :: idofTime

  ! Weight for the operator
  real(DP), intent(in) :: dweight

  ! Current control, solution of the control equation.
  type(t_controlSpace), intent(inout) :: rcontrol

  ! Solution of the linearised control equation.
  type(t_controlSpace), intent(inout) :: rcontrolLin
!</input>

!<inputoutput>
  ! Block vector receiving the result.
  type(t_vectorBlock), intent(inout) :: rrhs
!</inputoutput>

!</subroutine>

    ! local variables    
    type(p_t_spacetimeOperatorAsm), pointer :: p_p_rspatialOperatorAsm
    type(t_collection) :: rcollection
    type(t_collection), target :: ruserCollection
    type(t_fev2Vectors) :: rvectorEval
    type(t_vectorBlock), pointer :: p_rvector1,p_rvector2
    real(DP) :: dtheta, dtstep, dtime, dtimeend, dtimestart

    ! Cancel if nothing to do
    if (dweight .eq. 0.0_DP) return

    ! Prepare a local and a user-defined collection
    call collct_init (rcollection)
    call collct_init (ruserCollection)
    rcollection%p_rnextCollection => ruserCollection
    
    ! Prepare a pointer to our operator for callback routines.
    ! It is passed via rcollection.
    p_p_rspatialOperatorAsm%p_rspatialOperatorAsm => rspatialOperatorAsm
    rcollection%IquickAccess(8:) = transfer(p_p_rspatialOperatorAsm,rcollection%IquickAccess(8:))
    rcollection%DquickAccess(1) = dweight
    
    ! Notify the callback routine what to assemble.
    rcollection%IquickAccess(1) = OPTP_PRIMALLIN
    
    ! Timestepping technique?
    select case (rspatialOperatorAsm%p_rtimeDiscr%ctype)
    
    ! ***********************************************************
    ! Standard Theta one-step scheme.
    ! ***********************************************************
    case (TDISCR_ONESTEPTHETA)
    
      ! Theta-scheme identifier
      dtheta = rspatialOperatorAsm%p_rtimeDiscr%dtheta
      
      ! Characteristics of the current timestep.
      call tdiscr_getTimestep(rspatialOperatorAsm%p_rtimeDiscr,idofTime-1,&
          dtimeend,dtstep,dtimestart)

      ! itag=0: old 1-step scheme.
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      select case (rspatialOperatorAsm%p_rtimeDiscr%itag)
      
      ! ***********************************************************
      ! itag=0: old/standard 1-step scheme.
      ! ***********************************************************
      case (0)

        call output_line("Old 1-step-scheme not implemented",&
            OU_CLASS_ERROR,OU_MODE_STD,"smva_getRhs_primalLin")
        call sys_halt()

      ! ***********************************************************
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      ! ***********************************************************
      case (1)

        ! Which equation do we have?
        select case (rspatialOperatorAsm%p_rphysics%cequation)

        ! ***********************************************************
        ! Stokes/Navier Stokes
        ! ***********************************************************
        case (0,1)
        
          ! Evaluation point of the RHS in time
          dtime = dtimestart + (1.0_DP-dtheta) * dtstep
          
          ! Prepare the user-defined collection for the assembly
          call user_initCollectForAssembly (rspatialOperatorAsm%p_rglobalData,dtime,rusercollection)

          ! Prepare the evaluation.
          !
          ! Add the dual velocity if we have distributed control. Otherwise add dummy
          ! vectors which take no time in being computed.
          if (rspatialOperatorAsm%p_rsettingsOptControl%dalphaC .ge. 0.0_DP) then
            ! Position 1+2 = control
            call sptivec_getVectorFromPool (rcontrol%p_rvectorAccess,idofTime,p_rvector1)
            call fev2_addVectorToEvalList(rvectorEval,p_rvector1%RvectorBlock(1),0)
            call fev2_addVectorToEvalList(rvectorEval,p_rvector1%RvectorBlock(2),0)

            ! Position 3+4 = update for the control
            call sptivec_getVectorFromPool (rcontrolLin%p_rvectorAccess,idofTime,p_rvector2)
            call fev2_addVectorToEvalList(rvectorEval,p_rvector2%RvectorBlock(1),0)
            call fev2_addVectorToEvalList(rvectorEval,p_rvector2%RvectorBlock(2),0)

          else
            call fev2_addDummyVectorToEvalList(rvectorEval)
            call fev2_addDummyVectorToEvalList(rvectorEval)

            call fev2_addDummyVectorToEvalList(rvectorEval)
            call fev2_addDummyVectorToEvalList(rvectorEval)
          end if

          ! Build the vector
          call bma_buildVector (rrhs,BMA_CALC_STANDARD,&
              smva_fcalc_rhs, rcollection, revalVectors=rvectorEval,&
              rcubatureInfo=rspatialOperatorAsm%p_rasmTemplates%rcubatureInfoRHScontinuity)
          
          ! Cleanup
          call fev2_releaseVectorList(rvectorEval)
          call user_doneCollectForAssembly (rspatialOperatorAsm%p_rglobalData,rusercollection)
          
        end select ! Equation
      
      end select ! Timestep sub-scheme
      
    end select ! Timestep scheme
    
    ! Release the collections
    call collct_done (ruserCollection)
    call collct_init (rcollection)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine smva_getRhs_dualLin (rspatialOperatorAsm,idofTime,rdualSol,rprimalLinSol,&
      bfull,dweight,rrhs)

!<description>
  ! Implements semilinear parts of the linearised dual operator
  ! into the right-hand side vector rrhs, weighted by dweight.
!</description>

!<input>
  ! Definition of the operator. Must be a linearised operator.
  type(t_spacetimeOperatorAsm), target :: rspatialOperatorAsm
  
  ! Number of the DOF in time which should be calculated into rrhs.
  integer, intent(in) :: idofTime

  ! TRUE activates the full linearised operator (Newton).
  ! FALSE activates a partially linearised operator without the Newton part.
  logical, intent(in) :: bfull

  ! Weight for the operator
  real(DP), intent(in) :: dweight

  ! Space-time vector which contains the solution of the dual equation.
  type(t_dualSpace), intent(inout) :: rdualSol

  ! Space-time vector which contains the solution of the linearised primal equation.
  type(t_dualSpace), intent(inout) :: rprimalLinSol
!</input>

!<inputoutput>
  ! Block vector receiving the result.
  type(t_vectorBlock), intent(inout) :: rrhs
!</inputoutput>

!</subroutine>

    ! local variables    
    type(p_t_spacetimeOperatorAsm), pointer :: p_p_rspatialOperatorAsm
    type(t_collection) :: rcollection
    type(t_collection), target :: ruserCollection
    type(t_fev2Vectors) :: rvectorEval
    type(t_vectorBlock), pointer :: p_rvector1,p_rvector2
    real(DP) :: dtheta, dtstep, dtime, dtimeend, dtimestart

    ! Cancel if nothing to do
    if (dweight .eq. 0.0_DP) return

    ! Prepare a local and a user-defined collection
    call collct_init (rcollection)
    call collct_init (ruserCollection)
    rcollection%p_rnextCollection => ruserCollection
    
    ! Prepare a pointer to our operator for callback routines.
    ! It is passed via rcollection.
    p_p_rspatialOperatorAsm%p_rspatialOperatorAsm => rspatialOperatorAsm
    rcollection%IquickAccess(8:) = transfer(p_p_rspatialOperatorAsm,rcollection%IquickAccess(8:))
    rcollection%DquickAccess(1) = dweight
    
    ! Notify the callback routine what to assemble.
    if (bfull) then
      rcollection%IquickAccess(1) = OPTP_DUALLIN
    else
      rcollection%IquickAccess(1) = OPTP_DUALLIN_SIMPLE
    end if
    
    ! Timestepping technique?
    select case (rspatialOperatorAsm%p_rtimeDiscr%ctype)
    
    ! ***********************************************************
    ! Standard Theta one-step scheme.
    ! ***********************************************************
    case (TDISCR_ONESTEPTHETA)
    
      ! Theta-scheme identifier
      dtheta = rspatialOperatorAsm%p_rtimeDiscr%dtheta
      
      ! Characteristics of the current timestep.
      call tdiscr_getTimestep(rspatialOperatorAsm%p_rtimeDiscr,idofTime-1,&
          dtimeend,dtstep,dtimestart)

      ! itag=0: old 1-step scheme.
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      select case (rspatialOperatorAsm%p_rtimeDiscr%itag)
      
      ! ***********************************************************
      ! itag=0: old/standard 1-step scheme.
      ! ***********************************************************
      case (0)

        call output_line("Old 1-step-scheme not implemented",&
            OU_CLASS_ERROR,OU_MODE_STD,"smva_getRhs_dualLin")
        call sys_halt()

      ! ***********************************************************
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      ! ***********************************************************
      case (1)

        ! Which equation do we have?
        select case (rspatialOperatorAsm%p_rphysics%cequation)

        ! ***********************************************************
        ! Stokes/Navier Stokes
        ! ***********************************************************
        case (0,1)
        
          ! Evaluation point of the RHS in time
          dtime = dtimestart
          
          ! Prepare the user-defined collection for the assembly
          call user_initCollectForAssembly (rspatialOperatorAsm%p_rglobalData,dtime,rusercollection)

          ! Prepare the evaluation.
          !
          ! Vector 1+2 = dual velocity.
          call sptivec_getVectorFromPool (rdualSol%p_rvectorAccess,idofTime,p_rvector1)
          call fev2_addVectorToEvalList(rvectorEval,p_rvector1%RvectorBlock(1),0)
          call fev2_addVectorToEvalList(rvectorEval,p_rvector1%RvectorBlock(2),0)

          ! Vector 3+4 = linearised primal velocity. We need the 1st
          ! derivative as well.
          call sptivec_getVectorFromPool (rprimalLinSol%p_rvectorAccess,idofTime,p_rvector2)
          call fev2_addVectorToEvalList(rvectorEval,p_rvector2%RvectorBlock(1),1)
          call fev2_addVectorToEvalList(rvectorEval,p_rvector2%RvectorBlock(2),1)
          
          ! Build the vector
          call bma_buildVector (rrhs,BMA_CALC_STANDARD,&
              smva_fcalc_semilinRhs, rcollection, revalVectors=rvectorEval,&
              rcubatureInfo=rspatialOperatorAsm%p_rasmTemplates%rcubatureInfoRHScontinuity)
          
          ! Cleanup
          call fev2_releaseVectorList(rvectorEval)
          call user_doneCollectForAssembly (rspatialOperatorAsm%p_rglobalData,rusercollection)
          
        end select ! Equation
        
      end select ! Timestep sub-scheme
      
    end select ! Timestep scheme
    
    ! Release the collections
    call collct_done (ruserCollection)
    call collct_init (rcollection)
    
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine smva_getSemilinRhs_dualLin (rspatialOperatorAsm,idofTime,rdualSol,rprimalLinSol,&
      bfull,dweight,rrhs)

!<description>
  ! Implements semilinear parts of the linearised dual operator
  ! into the right-hand side vector rrhs, weighted by dweight.
!</description>

!<input>
  ! Definition of the operator. Must be a linearised operator.
  type(t_spacetimeOperatorAsm), target :: rspatialOperatorAsm
  
  ! Number of the DOF in time which should be calculated into rrhs.
  integer, intent(in) :: idofTime

  ! TRUE activates the full linearised operator (Newton).
  ! FALSE activates a partially linearised operator without the Newton part.
  logical, intent(in) :: bfull

  ! Weight for the operator
  real(DP), intent(in) :: dweight

  ! Space-time vector which contains the solution of the dual equation.
  type(t_dualSpace), intent(inout) :: rdualSol

  ! Space-time vector which contains the solution of the linearised primal equation.
  type(t_dualSpace), intent(inout) :: rprimalLinSol
!</input>

!<inputoutput>
  ! Block vector receiving the result.
  type(t_vectorBlock), intent(inout) :: rrhs
!</inputoutput>

!</subroutine>

    ! local variables    
    type(p_t_spacetimeOperatorAsm), pointer :: p_p_rspatialOperatorAsm
    type(t_collection) :: rcollection
    type(t_collection), target :: ruserCollection
    type(t_fev2Vectors) :: rvectorEval
    type(t_vectorBlock), pointer :: p_rvector1,p_rvector2
    real(DP) :: dtheta, dtstep, dtimeend, dtimestart

    ! Cancel if nothing to do
    if (dweight .eq. 0.0_DP) return

    ! Prepare a local and a user-defined collection
    call collct_init (rcollection)
    call collct_init (ruserCollection)
    rcollection%p_rnextCollection => ruserCollection
    
    ! Prepare a pointer to our operator for callback routines.
    ! It is passed via rcollection.
    p_p_rspatialOperatorAsm%p_rspatialOperatorAsm => rspatialOperatorAsm
    rcollection%IquickAccess(8:) = transfer(p_p_rspatialOperatorAsm,rcollection%IquickAccess(8:))
    rcollection%DquickAccess(1) = dweight
    
    ! Notify the callback routine what to assemble.
    if (bfull) then
      rcollection%IquickAccess(1) = OPTP_DUALLIN
    else
      rcollection%IquickAccess(1) = OPTP_DUALLIN_SIMPLE
    end if
    
    ! Timestepping technique?
    select case (rspatialOperatorAsm%p_rtimeDiscr%ctype)
    
    ! ***********************************************************
    ! Standard Theta one-step scheme.
    ! ***********************************************************
    case (TDISCR_ONESTEPTHETA)
    
      ! Theta-scheme identifier
      dtheta = rspatialOperatorAsm%p_rtimeDiscr%dtheta
      
      ! Characteristics of the current timestep.
      call tdiscr_getTimestep(rspatialOperatorAsm%p_rtimeDiscr,idofTime-1,&
          dtimeend,dtstep,dtimestart)

      ! itag=0: old 1-step scheme.
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      select case (rspatialOperatorAsm%p_rtimeDiscr%itag)
      
      ! ***********************************************************
      ! itag=0: old/standard 1-step scheme.
      ! ***********************************************************
      case (0)

        call output_line("Old 1-step-scheme not implemented",&
            OU_CLASS_ERROR,OU_MODE_STD,"smva_getSemilinRhs_dualLin")
        call sys_halt()

      ! ***********************************************************
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      ! ***********************************************************
      case (1)

        ! Which equation do we have?
        select case (rspatialOperatorAsm%p_rphysics%cequation)

        ! ***********************************************************
        ! Navier Stokes
        ! ***********************************************************
        case (1)
        
          ! Prepare the evaluation.
          !
          ! Vector 1+2 = dual velocity.
          call sptivec_getVectorFromPool (rdualSol%p_rvectorAccess,idofTime,p_rvector1)
          call fev2_addVectorToEvalList(rvectorEval,p_rvector1%RvectorBlock(1),0)
          call fev2_addVectorToEvalList(rvectorEval,p_rvector1%RvectorBlock(2),0)

          ! Vector 3+4 = linearised primal velocity.
          call sptivec_getVectorFromPool (rprimalLinSol%p_rvectorAccess,idofTime,p_rvector2)
          call fev2_addVectorToEvalList(rvectorEval,p_rvector2%RvectorBlock(1),1)
          call fev2_addVectorToEvalList(rvectorEval,p_rvector2%RvectorBlock(2),1)
          
          ! Build the vector
          call bma_buildVector (rrhs,BMA_CALC_STANDARD,&
              smva_fcalc_semilinRhs, rcollection, revalVectors=rvectorEval,&
              rcubatureInfo=rspatialOperatorAsm%p_rasmTemplates%rcubatureInfoRHScontinuity)
          
          ! Cleanup
          call fev2_releaseVectorList(rvectorEval)
          
        end select ! Equation

      end select ! Timestep sub-scheme
            
    end select ! Timestep scheme
    
    ! Release the collections
    call collct_done (ruserCollection)
    call collct_init (rcollection)
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine smva_fcalc_semilinRhs(rvectorData,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the semilinear part of the system matrix.
!</description>

!<inputoutput>
    ! Vector data of all subvectors. The arrays p_Dentry of all subvectors
    ! have to be filled with data.
    type(t_bmaVectorData), dimension(:), intent(inout), target :: rvectorData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaVectorAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaVectorAssembly), intent(in) :: rvectorAssembly
    
    ! Number of points per element
    integer, intent(in) :: npointsPerElement
    
    ! Number of elements
    integer, intent(in) :: nelements
    
    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</input>

!<subroutine>

    ! Local variables
    real(DP) :: dbasI, dbasIx, dbasIy
    real(DP) :: dlambda1, dlambda2, dylin1, dylin2, dylin1x, dylin1y, dylin2x, dylin2y
    real(DP) :: dweight
    integer :: iel, icubp, idofe
    integer :: copType
    real(DP), dimension(:,:), pointer :: p_DlocalVector1,p_DlocalVector2
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaVectorData), pointer :: p_rvectorData1,p_rvectorData2
    real(DP), dimension(:,:,:), pointer :: p_Dlambda1,p_Dlambda2
    real(DP), dimension(:,:,:), pointer :: p_Dylin1,p_Dylin2

    type(p_t_spacetimeOperatorAsm), pointer :: p_p_rspatialOperatorAsm
    type(t_spacetimeOperatorAsm), pointer :: p_rspatialOperatorAsm

    ! From the collection, fetch our operator structure, nonlinearity,...
    p_p_rspatialOperatorAsm = transfer(rcollection%IquickAccess(8:),p_p_rspatialOperatorAsm)
    p_rspatialOperatorAsm => p_p_rspatialOperatorAsm%p_rspatialOperatorAsm

    ! Type of the operator to compute.
    copType = rcollection%IquickAccess(1)

    ! Weight of the operator.
    dweight = rcollection%DquickAccess(1)
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    
    ! Which equation do we have?
    select case (p_rspatialOperatorAsm%p_rphysics%cequation)
    case (0,1)
    
      ! -------------------------------------------------------------
      ! Stokes/Navier Stokes.
      ! -------------------------------------------------------------
      
      ! Primal or dual equation?
      select case (copType)

      ! ***********************************************************
      ! Navier-Stokes. Linearised backward equation.
      ! ***********************************************************
      case (OPTP_DUALLIN)
      
        ! The additional operator assembled here and added to the
        ! right-hand side vector reads as follows:
        !
        !   dweight * ( lambda , (ylin grad) (phi) + grad(y) phi )
        !
        ! with lambda being the solution of the dual equation
        ! and ylin being the solution of the linearised primal equation.
        ! Writing this operator in a component-wise way is a bit
        ! tideous, see below...
        !
        ! This calculation shall not be included in the RHS calculation
        ! as it is used with different timestep weight than the RHS.

        ! Get the nonlinearity
        p_Dlambda1 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,:)
        p_Dlambda2 => revalVectors%p_RvectorData(2)%p_Ddata(:,:,:)
        p_Dylin1 => revalVectors%p_RvectorData(3)%p_Ddata(:,:,:)
        p_Dylin2 => revalVectors%p_RvectorData(4)%p_Ddata(:,:,:)
      
        ! Get the data arrays of the subvector
        p_rvectorData1 => RvectorData(1)
        p_rvectorData2 => RvectorData(2)
        
        p_DlocalVector1 => RvectorData(1)%p_Dentry
        p_DlocalVector2 => RvectorData(2)%p_Dentry
        
        p_DbasTest => RvectorData(1)%p_DbasTest
        
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Phi_i:
            do idofe=1,p_rvectorData1%ndofTest
            
              ! Fetch the contributions of the (test) basis functions Phi_i
              ! into dbasI
              dbasI  = p_DbasTest(idofe,DER_FUNC,icubp,iel)
              dbasIx = p_DbasTest(idofe,DER_DERIV2D_X,icubp,iel)
              dbasIy = p_DbasTest(idofe,DER_DERIV2D_Y,icubp,iel)
              
              ! Get the values of lambda and ylin.
              dlambda1 = p_Dlambda1(icubp,iel,DER_FUNC)
              dlambda2 = p_Dlambda2(icubp,iel,DER_FUNC)
              
              dylin1  = p_Dylin1(icubp,iel,DER_FUNC)
              dylin1x = p_Dylin1(icubp,iel,DER_DERIV_X)
              dylin1y = p_Dylin1(icubp,iel,DER_DERIV_Y)
              
              dylin2  = p_Dylin2(icubp,iel,DER_FUNC)
              dylin2x = p_Dylin2(icubp,iel,DER_DERIV_X)
              dylin2y = p_Dylin2(icubp,iel,DER_DERIV_Y)
              
              ! Multiply the values of the basis functions
              ! (1st derivatives) by the cubature weight and sum up
              ! into the local vectors.
              p_DlocalVector1(idofe,iel) = p_DlocalVector1(idofe,iel) + &
                  dweight * p_DcubWeight(icubp,iel) * &
                  ( dlambda1 * ( dylin1*dbasIx + dylin2*dbasIy + dylin1x*dbasI ) + &
                    dlambda2 * dylin2x * dbasI )

              p_DlocalVector2(idofe,iel) = p_DlocalVector2(idofe,iel) + &
                  dweight * p_DcubWeight(icubp,iel) * &
                  ( dlambda1 * dylin1y * dbasI + &
                    dlambda2 * ( dylin1*dbasIx + dylin2*dbasIy + dylin2y*dbasI ) )
                  
            end do ! jdofe

          end do ! icubp
        
        end do ! iel

      end select
      
    end select    
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine smva_fcalc_rhs(rvectorData,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the general right-hand side vector of any equation.
    ! The RHS to calculate is specified via the collection.
!</description>

!<inputoutput>
    ! Vector data of all subvectors. The arrays p_Dentry of all subvectors
    ! have to be filled with data.
    type(t_bmaVectorData), dimension(:), intent(inout), target :: rvectorData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaVectorAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaVectorAssembly), intent(in) :: rvectorAssembly
    
    ! Number of points per element
    integer, intent(in) :: npointsPerElement
    
    ! Number of elements
    integer, intent(in) :: nelements
    
    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), optional :: rcollection
!</input>

!<subroutine>

    ! Local variables
    real(DP) :: dbasI
    real(DP) :: dylin1, dylin2
    real(DP) :: dy1,dy2,dz1,dz2,du1,du2,dulin1,dulin2,drhs1,drhs2
    real(DP) :: dx,dy,dweight
    real(DP) :: dumin1,dumax1,dumin2,dumax2
    real(DP) :: dalpha
    integer :: iel, icubp, idofe, ierror, iid
    integer :: copType
    real(DP), dimension(:,:), pointer :: p_DlocalVector1,p_DlocalVector2
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaVectorData), pointer :: p_rvectorData1,p_rvectorData2
    real(DP), dimension(:,:,:), pointer :: p_Du1,p_Du2
    real(DP), dimension(:,:,:), pointer :: p_DuLin1,p_DuLin2
    real(DP), dimension(:,:,:), pointer :: p_Dylin1,p_Dylin2
    real(DP), dimension(:,:,:), pointer :: p_Dy1,p_Dy2
    real(DP), dimension(:,:,:), pointer :: p_Dz1,p_Dz2
    real(DP), dimension(:,:,:), pointer :: p_Drhs1,p_Drhs2
    real(DP), dimension(:), pointer :: p_DobservationArea => null()

    type(p_t_spacetimeOperatorAsm), pointer :: p_p_rspatialOperatorAsm
    type(t_spacetimeOperatorAsm), pointer :: p_rspatialOperatorAsm

    ! From the collection, fetch our operator structure, nonlinearity,...
    p_p_rspatialOperatorAsm = transfer(rcollection%IquickAccess(8:),p_p_rspatialOperatorAsm)
    p_rspatialOperatorAsm => p_p_rspatialOperatorAsm%p_rspatialOperatorAsm

    ! Type of the operator to compute.
    copType = rcollection%IquickAccess(1)

    ! Weight of the operator.
    dweight = rcollection%DquickAccess(1)
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    
    ! Which equation do we have?
    select case (p_rspatialOperatorAsm%p_rphysics%cequation)
    case (0,1)
    
      ! -------------------------------------------------------------
      ! Stokes/Navier Stokes.
      ! -------------------------------------------------------------
      
      ! Primal or dual equation?
      select case (copType)

      ! ***********************************************************
      ! Navier-Stokes. Forward equation.
      ! ***********************************************************
      case (OPTP_PRIMAL)
      
        ! Get the data arrays of the subvector
        p_rvectorData1 => RvectorData(1)
        p_rvectorData2 => RvectorData(2)
        
        p_DlocalVector1 => RvectorData(1)%p_Dentry
        p_DlocalVector2 => RvectorData(2)%p_Dentry
        
        p_DbasTest => RvectorData(1)%p_DbasTest

        ! ------------------------------------------------
        ! Calculate the user-defined RHS in the
        ! cubature points
        ! ------------------------------------------------
      
        ! Get memory for the user-defined right-hand side f.
        p_Drhs1 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,:)
        p_Drhs2 => revalVectors%p_RvectorData(2)%p_Ddata(:,:,:)

        ! The right-hand side is given as an analytical function
        ! or as a discrete function on an arbitrary mesh, discretised
        ! with an arbitrary finite element. We therefore cannot
        ! automatically precalculate the values in the cubature
        ! points, but have to do it right here.
        !
        ! Evaluate the analytic function in the cubature points. 1st component
        call ansol_evaluate (rcollection,"RHS",1,p_Drhs1(:,:,DER_FUNC),&
            npointsPerElement,nelements,rassemblyData%revalElementSet%p_DpointsReal,&
            rassemblyData%p_IelementList,ierror,iid)

        if (ierror .eq. 1) then
        
          ! This is an error, something went wrong.
          call output_line("Cannot evaluate function",&
              OU_CLASS_ERROR,OU_MODE_STD,"smva_fcalc_semilinRhs")
          call sys_halt()
        
        else  if (ierror .eq. -1) then
        
          ! Oops, target function is a user-defined one.
          ! Call the user-defined evaluation routine, providing iid
          ! as function id.
          
          call output_line("User defined function not implemented.",&
              OU_CLASS_ERROR,OU_MODE_STD,"smva_fcalc_semilinRhs")
          call sys_halt()
        
        end if

        ! Evaluate the analytic function in the cubature points. 2nd component.
        call ansol_evaluate (rcollection,"RHS",2,p_Drhs2(:,:,DER_FUNC),&
            npointsPerElement,nelements,rassemblyData%revalElementSet%p_DpointsReal,&
            rassemblyData%p_IelementList,ierror,iid)

        if (ierror .eq. 1) then
        
          ! This is an error, something went wrong.
          call output_line("Cannot target function",&
              OU_CLASS_ERROR,OU_MODE_STD,"smva_fcalc_semilinRhs")
          call sys_halt()
        
        else  if (ierror .eq. -1) then
        
          ! Oops, target function is a user-defined one.
          ! Call the user-defined evaluation routine, providing iid
          ! as function id.
          call output_line("User defined function not implemented.",&
              OU_CLASS_ERROR,OU_MODE_STD,"smva_fcalc_semilinRhs")
          call sys_halt()
        
        end if

        ! Ok, now we have the right-and side.

        ! ------------------------------------------------
        ! Calculate the user-defined RHS 
        ! ------------------------------------------------

        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Phi_i:
            do idofe=1,p_rvectorData1%ndofTest
            
              ! Fetch the contributions of the (test) basis functions Phi_i
              ! into dbasI
              dbasI  = p_DbasTest(idofe,DER_FUNC,icubp,iel)
              
              ! Get the RHS
              drhs1 = p_Drhs1(icubp,iel,DER_FUNC)
              drhs2 = p_Drhs2(icubp,iel,DER_FUNC)

              ! Calculate the entries in the RHS.
              p_DlocalVector1(idofe,iel) = p_DlocalVector1(idofe,iel) + &
                  dweight * p_DcubWeight(icubp,iel) * &
                  drhs1 * dbasI 

              p_DlocalVector2(idofe,iel) = p_DlocalVector2(idofe,iel) + &
                  dweight * p_DcubWeight(icubp,iel) * &
                  drhs2 * dbasI 
                  
            end do ! jdofe

          end do ! icubp
        
        end do ! iel

        ! ------------------------------------------------
        ! Calculate the problem-specific RHS
        ! ------------------------------------------------

        ! For the primal equation, the right-hand side reads
        !
        !    rhs = (user defined) + u
        !
        ! with u being the control. The control is already evaluated
        ! in the cubature points, we can just take the values.
        
        ! Check the regularisation parameter ALPHA. Do we have
        ! distributed control?
        dalpha = p_rspatialOperatorAsm%p_rsettingsOptControl%dalphaC
        
        if (dalpha .ge. 0.0_DP) then

          ! Get the control.
          p_Du1 => revalVectors%p_RvectorData(3)%p_Ddata(:,:,:)
          p_Du2 => revalVectors%p_RvectorData(4)%p_Ddata(:,:,:)

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Phi_i:
              do idofe=1,p_rvectorData1%ndofTest
              
                ! Fetch the contributions of the (test) basis functions Phi_i
                ! into dbasI
                dbasI  = p_DbasTest(idofe,DER_FUNC,icubp,iel)
                
                ! Get the values of the control.
                du1 = p_Du1(icubp,iel,DER_FUNC)
                du2 = p_Du2(icubp,iel,DER_FUNC)
                
                ! Calculate the entries in the RHS.
                p_DlocalVector1(idofe,iel) = p_DlocalVector1(idofe,iel) + &
                    dweight * p_DcubWeight(icubp,iel) * &
                    du1 * dbasI 

                p_DlocalVector2(idofe,iel) = p_DlocalVector2(idofe,iel) + &
                    dweight * p_DcubWeight(icubp,iel) * &
                    du2 * dbasI 
                    
              end do ! jdofe

            end do ! icubp
          
          end do ! iel
          
        end if
            
      ! ***********************************************************
      ! Navier-Stokes. Linearised forward equation.
      ! ***********************************************************
      case (OPTP_PRIMALLIN)
      
        ! For the linearised primal equation, the right-hand side reads
        !
        !    rhs = (user defined) + u'
        !
        ! with u being the control. If we have constraints, the 
        ! right-hand side reads
        !
        !    rhs = (user defined) + DP(u) u'
        !
        ! DP(.) is the Newton derivative of the projection P(.). 
        ! In case of Box constraints, there is
        !
        !   DP(u) = identity  where u is in the bounds
        !         = 0         where u violates the bounds
        !
        ! Check the regularisation parameter ALPHA. Do we have
        ! distributed control?
        dalpha = p_rspatialOperatorAsm%p_rsettingsOptControl%dalphaC
        
        if (dalpha .gt. 0.0_DP) then

          ! Get the nonlinearity lambda and its current linearisation LambdaLin=lambda'
          p_Du1 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,:)
          p_Du2 => revalVectors%p_RvectorData(2)%p_Ddata(:,:,:)

          ! Get the nonlinearity lambda and its current linearisation LambdaLin=lambda'
          p_DuLin1 => revalVectors%p_RvectorData(3)%p_Ddata(:,:,:)
          p_DuLin2 => revalVectors%p_RvectorData(4)%p_Ddata(:,:,:)

          ! Get the data arrays of the subvector
          p_rvectorData1 => RvectorData(1)
          p_rvectorData2 => RvectorData(2)
          
          p_DlocalVector1 => RvectorData(1)%p_Dentry
          p_DlocalVector2 => RvectorData(2)%p_Dentry
          
          p_DbasTest => RvectorData(1)%p_DbasTest

          select case (p_rspatialOperatorAsm%p_rsettingsOptControl%rconstraints%cdistVelConstraints)
            
          ! ---------------------------------------------------------
          ! No constraints
          ! ---------------------------------------------------------
          case (0)
            
            ! Loop over the elements in the current set.
            do iel = 1,nelements

              ! Loop over all cubature points on the current element
              do icubp = 1,npointsPerElement

                ! Outer loop over the DOF's i=1..ndof on our current element,
                ! which corresponds to the (test) basis functions Phi_i:
                do idofe=1,p_rvectorData1%ndofTest
                
                  ! Fetch the contributions of the (test) basis functions Phi_i
                  ! into dbasI
                  dbasI  = p_DbasTest(idofe,DER_FUNC,icubp,iel)
                  
                  ! Get the values of lambda'
                  duLin1 = p_Du1(icubp,iel,DER_FUNC)
                  duLin2 = p_Du2(icubp,iel,DER_FUNC)

                  ! Calculate the entries in the RHS.
                  p_DlocalVector1(idofe,iel) = p_DlocalVector1(idofe,iel) + &
                      dweight * p_DcubWeight(icubp,iel) * &
                      dulin1 * dbasI 

                  p_DlocalVector2(idofe,iel) = p_DlocalVector2(idofe,iel) + &
                      dweight * p_DcubWeight(icubp,iel) * &
                      dulin2 * dbasI 
                      
                end do ! jdofe

              end do ! icubp
            
            end do ! iel
            
          ! ---------------------------------------------------------
          ! Constant box constraints
          ! ---------------------------------------------------------
          case (1)
            
            ! Get the box constraints
            dumin1 = p_rspatialOperatorAsm%p_rsettingsOptControl%rconstraints%ddistVelUmin1
            dumin2 = p_rspatialOperatorAsm%p_rsettingsOptControl%rconstraints%ddistVelUmin2
            dumax1 = p_rspatialOperatorAsm%p_rsettingsOptControl%rconstraints%ddistVelUmax1
            dumax2 = p_rspatialOperatorAsm%p_rsettingsOptControl%rconstraints%ddistVelUmax2

            ! Loop over the elements in the current set.
            do iel = 1,nelements

              ! Loop over all cubature points on the current element
              do icubp = 1,npointsPerElement

                ! Get the control
                du1 = p_Du1(icubp,iel,DER_FUNC)
                du2 = p_Du2(icubp,iel,DER_FUNC)
                
                ! If the control violates the bounds, the Newton derivative
                ! is zero, so only calculate something if the bounds
                ! are not violated.

                if ((du1 .gt. dumin1) .and. (du1 .lt. dumax1)) then
                  
                  ! Outer loop over the DOF's i=1..ndof on our current element,
                  ! which corresponds to the (test) basis functions Phi_i:
                  do idofe=1,p_rvectorData1%ndofTest
                  
                    ! Fetch the contributions of the (test) basis functions Phi_i
                    ! into dbasI
                    dbasI  = p_DbasTest(idofe,DER_FUNC,icubp,iel)
                  
                    ! Calculate the linearised control
                    duLin1 = p_DuLin1(icubp,iel,DER_FUNC)

                    ! Calculate the entries in the RHS.
                    p_DlocalVector1(idofe,iel) = p_DlocalVector1(idofe,iel) + &
                        dweight * p_DcubWeight(icubp,iel) * &
                        duLin1 * dbasI

                  end do ! idofe
                  
                end if
                  
                if ((du2 .gt. dumin2) .and. (du2 .lt. dumax2)) then

                  ! Outer loop over the DOF's i=1..ndof on our current element,
                  ! which corresponds to the (test) basis functions Phi_i:
                  do idofe=1,p_rvectorData1%ndofTest

                    ! Calculate the linearised control
                    duLin2 = p_DuLin2(icubp,iel,DER_FUNC)
      
                    ! Calculate the entries in the RHS.
                    p_DlocalVector2(idofe,iel) = p_DlocalVector2(idofe,iel) + &
                        dweight * p_DcubWeight(icubp,iel) * &
                        duLin2 * dbasI 
  
                  end do ! jdofe

                end if

              end do ! icubp
            
            end do ! iel

          end select ! constraints
          
        end if

      ! ***********************************************************
      ! Navier-Stokes. Backward equation.
      ! ***********************************************************
      case (OPTP_DUAL)
      
        ! Get the data arrays of the subvector
        p_rvectorData1 => RvectorData(1)
        p_rvectorData2 => RvectorData(2)
        
        p_DlocalVector1 => RvectorData(1)%p_Dentry
        p_DlocalVector2 => RvectorData(2)%p_Dentry
        
        p_DbasTest => RvectorData(1)%p_DbasTest
      
        ! ------------------------------------------------
        ! Calculate the user-defined RHS in the
        ! cubature points
        ! ------------------------------------------------

        ! Get memory for the user-defined right-hand side f.
        p_Drhs1 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,:)
        p_Drhs2 => revalVectors%p_RvectorData(2)%p_Ddata(:,:,:)

        ! The right-hand side is given as an analytical function
        ! or as a discrete function on an arbitrary mesh, discretised
        ! with an arbitrary finite element. We therefore cannot
        ! automatically precalculate the values in the cubature
        ! points, but have to do it right here.
        !
        ! Evaluate the analytic function in the cubature points. 1st component
        call ansol_evaluate (rcollection,"RHS",1,p_Drhs1(:,:,DER_FUNC),&
            npointsPerElement,nelements,rassemblyData%revalElementSet%p_DpointsReal,&
            rassemblyData%p_IelementList,ierror,iid)

        if (ierror .eq. 1) then
        
          ! This is an error, something went wrong.
          call output_line("Cannot evaluate function",&
              OU_CLASS_ERROR,OU_MODE_STD,"smva_fcalc_semilinRhs")
          call sys_halt()
        
        else  if (ierror .eq. -1) then
        
          ! Oops, target function is a user-defined one.
          ! Call the user-defined evaluation routine, providing iid
          ! as function id.
          
          call output_line("User defined function not implemented.",&
              OU_CLASS_ERROR,OU_MODE_STD,"smva_fcalc_semilinRhs")
          call sys_halt()
        
        end if

        ! Evaluate the analytic function in the cubature points. 2nd component.
        call ansol_evaluate (rcollection,"RHS",2,p_Drhs2(:,:,DER_FUNC),&
            npointsPerElement,nelements,rassemblyData%revalElementSet%p_DpointsReal,&
            rassemblyData%p_IelementList,ierror,iid)

        if (ierror .eq. 1) then
        
          ! This is an error, something went wrong.
          call output_line("Cannot target function",&
              OU_CLASS_ERROR,OU_MODE_STD,"smva_fcalc_semilinRhs")
          call sys_halt()
        
        else  if (ierror .eq. -1) then
        
          ! Oops, target function is a user-defined one.
          ! Call the user-defined evaluation routine, providing iid
          ! as function id.
          
          call output_line("User defined function not implemented.",&
              OU_CLASS_ERROR,OU_MODE_STD,"smva_fcalc_semilinRhs")
          call sys_halt()
        
        end if

        ! Ok, we have the right-and side now.

        ! ------------------------------------------------
        ! Calculate the user-defined RHS 
        ! ------------------------------------------------

        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Phi_i:
            do idofe=1,p_rvectorData1%ndofTest
            
              ! Fetch the contributions of the (test) basis functions Phi_i
              ! into dbasI
              dbasI  = p_DbasTest(idofe,DER_FUNC,icubp,iel)
              
              ! Get the RHS
              drhs1 = p_Drhs1(icubp,iel,DER_FUNC)
              drhs2 = p_Drhs2(icubp,iel,DER_FUNC)

              ! Calculate the entries in the RHS.
              p_DlocalVector1(idofe,iel) = p_DlocalVector1(idofe,iel) + &
                  dweight * p_DcubWeight(icubp,iel) * &
                  drhs1 * dbasI 

              p_DlocalVector2(idofe,iel) = p_DlocalVector2(idofe,iel) + &
                  dweight * p_DcubWeight(icubp,iel) * &
                  drhs2 * dbasI 
                  
            end do ! jdofe

          end do ! icubp
        
        end do ! iel

        ! ------------------------------------------------
        ! Calculate the target flow in the cubature points
        ! ------------------------------------------------

        ! For the dual equation, the right-hand side reads
        !
        !    rhs = (user defined) + (y - z)
        !
        ! so we have to calculate (y-z) and add/subtract it
        ! from the RHS. In case the observation area is restricted,
        ! the right-hand side reads
        !
        !    rhs = (user defined) +- Chi_Omega~ (y - z)
        !
        ! so it has to be forced to zero outside of the observation
        ! area Omega~. Thus, the RHS is basically nonlinear.

        ! Get memory for the target field z.
        ! We can recycle the temp memory vectors from above (position 1+2)
        p_Dz1 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,:)
        p_Dz2 => revalVectors%p_RvectorData(2)%p_Ddata(:,:,:)
        
        ! Ok, we have the right-and side now.

        ! The target field is given as an analytical function
        ! or as a discrete function on an arbitrary mesh, discretised
        ! with an arbitrary finite element. We therefore cannot
        ! automatically precalculate the values in the cubature
        ! points, but have to do it right here.
        !
        ! Evaluate the analytic function in the cubature points. 1st component.
        call ansol_evaluate (rcollection,"TARGET",1,p_Dz1(:,:,DER_FUNC),&
            npointsPerElement,nelements,rassemblyData%revalElementSet%p_DpointsReal,&
            rassemblyData%p_IelementList,ierror,iid)

        if (ierror .eq. 1) then
        
          ! This is an error, something went wrong.
          call output_line("Cannot evaluate function",&
              OU_CLASS_ERROR,OU_MODE_STD,"smva_fcalc_semilinRhs")
          call sys_halt()
        
        else  if (ierror .eq. -1) then
        
          ! Oops, target function is a user-defined one.
          ! Call the user-defined evaluation routine, providing iid
          ! as function id.
          
          call output_line("User defined function not implemented.",&
              OU_CLASS_ERROR,OU_MODE_STD,"smva_fcalc_semilinRhs")
          call sys_halt()
        
        end if

        ! Evaluate the analytic function in the cubature points. 2nd component.
        call ansol_evaluate (rcollection,"TARGET",2,p_Dz2(:,:,DER_FUNC),&
            npointsPerElement,nelements,rassemblyData%revalElementSet%p_DpointsReal,&
            rassemblyData%p_IelementList,ierror,iid)

        if (ierror .eq. 1) then
        
          ! This is an error, something went wrong.
          call output_line("Cannot evaluate function",&
              OU_CLASS_ERROR,OU_MODE_STD,"smva_fcalc_semilinRhs")
          call sys_halt()
        
        else  if (ierror .eq. -1) then
        
          ! Oops, target function is a user-defined one.
          ! Call the user-defined evaluation routine, providing iid
          ! as function id.
          
          call output_line("User defined function not implemented.",&
              OU_CLASS_ERROR,OU_MODE_STD,"smva_fcalc_semilinRhs")
          call sys_halt()
        
        end if

        ! Ok, we have the target now.

        ! ------------------------------------------------
        ! Calculate the problem-specific RHS
        ! ------------------------------------------------
      
        ! Get the nonlinearity y.
        p_Dy1 => revalVectors%p_RvectorData(3)%p_Ddata(:,:,:)
        p_Dy2 => revalVectors%p_RvectorData(4)%p_Ddata(:,:,:)

        if (.not. associated(p_rspatialOperatorAsm%p_rsettingsOptControl%p_DobservationArea)) then
          
          ! ---------------------------------------------------------
          ! Observation area is the complete domain
          ! ---------------------------------------------------------
          
          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Phi_i:
              do idofe=1,p_rvectorData1%ndofTest
              
                ! Fetch the contributions of the (test) basis functions Phi_i
                ! into dbasI
                dbasI  = p_DbasTest(idofe,DER_FUNC,icubp,iel)
                
                ! Get the values of the RHS.
                drhs1 = p_Drhs1(icubp,iel,DER_FUNC)
                drhs2 = p_Drhs2(icubp,iel,DER_FUNC)
                
                ! Get the values of lambda and ylin.
                dy1 = p_Dy1(icubp,iel,DER_FUNC)
                dy2 = p_Dy2(icubp,iel,DER_FUNC)

                ! Get the target flow
                dz1 = p_Dz1(icubp,iel,DER_FUNC)
                dz2 = p_Dz2(icubp,iel,DER_FUNC)
                
                ! Calculate the entries in the RHS.
                p_DlocalVector1(idofe,iel) = p_DlocalVector1(idofe,iel) + &
                    dweight * p_DcubWeight(icubp,iel) * &
                    (dy1 - dz1) * dbasI 

                p_DlocalVector2(idofe,iel) = p_DlocalVector2(idofe,iel) + &
                    dweight * p_DcubWeight(icubp,iel) * &
                    (dy2 - dz2) * dbasI 
                    
              end do ! jdofe

            end do ! icubp
          
          end do ! iel
          
        else
        
          ! ---------------------------------------------------------
          ! Observation area is a rectangle
          ! ---------------------------------------------------------
          p_DobservationArea => p_rspatialOperatorAsm%p_rsettingsOptControl%p_DobservationArea
        
          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement
            
              ! Coordinates of the cubature point
              dx = rassemblyData%revalElementSet%p_DpointsReal(1,icubp,iel)
              dy = rassemblyData%revalElementSet%p_DpointsReal(1,icubp,iel)

              ! Calculate the entries in the RHS inside of the observation area
              if ((dx .ge. p_DobservationArea(1)) .and. (dy .ge. p_DobservationArea(2)) .and. &
                  (dx .le. p_DobservationArea(3)) .and. (dy .le. p_DobservationArea(4))) then
                  
                ! Outer loop over the DOF's i=1..ndof on our current element,
                ! which corresponds to the (test) basis functions Phi_i:
                do idofe=1,p_rvectorData1%ndofTest
                
                  ! Fetch the contributions of the (test) basis functions Phi_i
                  ! into dbasI
                  dbasI  = p_DbasTest(idofe,DER_FUNC,icubp,iel)
                  
                  ! Get the values of lambda and ylin.
                  dy1 = p_Dy1(icubp,iel,DER_FUNC)
                  dy2 = p_Dy2(icubp,iel,DER_FUNC)

                  ! Get the target flow
                  dz1 = p_Dz1(icubp,iel,DER_FUNC)
                  dz2 = p_Dz2(icubp,iel,DER_FUNC)
                  
                  p_DlocalVector1(idofe,iel) = p_DlocalVector1(idofe,iel) + &
                      dweight * p_DcubWeight(icubp,iel) * &
                      ( dy1 - dz1 ) * dbasI 

                  p_DlocalVector2(idofe,iel) = p_DlocalVector2(idofe,iel) + &
                      dweight * p_DcubWeight(icubp,iel) * &
                      ( dy2 - dz2 ) * dbasI 
                      
                end do ! jdofe
                
              end if ! Observation area

            end do ! icubp
          
          end do ! iel
        
        end if

      ! ***********************************************************
      ! Navier-Stokes. Linearised backward equation.
      ! ***********************************************************
      case (OPTP_DUALLIN)
      
        ! From the right-hand side 
        !
        !    rhs = y-z
        !
        ! there follows the addutional term
        !
        !    rhs = y'
        !
        ! in case the observation area is the complete domain.
        ! If the observation area is only a restricted domain Omega~,
        ! the RHS of the dual equation was
        !
        !    rhs = Chi_Omega~ (y-z)
        !
        ! so the rhs of the linearised dual equation is
        !
        !    rhs = Chi_Omega~ y'

        ! Get the nonlinearity
        p_Dylin1 => revalVectors%p_RvectorData(3)%p_Ddata(:,:,:)
        p_Dylin2 => revalVectors%p_RvectorData(4)%p_Ddata(:,:,:)
      
        ! Get the data arrays of the subvector
        p_rvectorData1 => RvectorData(1)
        p_rvectorData2 => RvectorData(2)
        
        p_DlocalVector1 => RvectorData(1)%p_Dentry
        p_DlocalVector2 => RvectorData(2)%p_Dentry
        
        p_DbasTest => RvectorData(1)%p_DbasTest
        
        if (.not. associated(p_rspatialOperatorAsm%p_rsettingsOptControl%p_DobservationArea)) then
          
          ! ---------------------------------------------------------
          ! Observation area is the complete domain
          ! ---------------------------------------------------------

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Phi_i:
              do idofe=1,p_rvectorData1%ndofTest
              
                ! Fetch the contributions of the (test) basis functions Phi_i
                ! into dbasI
                dbasI  = p_DbasTest(idofe,DER_FUNC,icubp,iel)
                
                dylin1  = p_Dylin1(icubp,iel,DER_FUNC)
                dylin2  = p_Dylin2(icubp,iel,DER_FUNC)
                
                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local vectors.
                p_DlocalVector1(idofe,iel) = p_DlocalVector1(idofe,iel) + &
                    dweight * p_DcubWeight(icubp,iel) * &
                    ( dylin1 * dbasI ) ! (y',phi)

                p_DlocalVector2(idofe,iel) = p_DlocalVector2(idofe,iel) + &
                    dweight * p_DcubWeight(icubp,iel) * &
                    ( dylin2 * dbasI ) ! (y',phi)
                    
              end do ! jdofe

            end do ! icubp
          
          end do ! iel
       
        else
        
          ! ---------------------------------------------------------
          ! Observation area is a rectangle
          ! ---------------------------------------------------------
          p_DobservationArea => p_rspatialOperatorAsm%p_rsettingsOptControl%p_DobservationArea

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Coordinates of the cubature point
              dx = rassemblyData%revalElementSet%p_DpointsReal(1,icubp,iel)
              dy = rassemblyData%revalElementSet%p_DpointsReal(1,icubp,iel)

              ! Calculate the entries in the RHS inside of the observation area
              if ((dx .ge. p_DobservationArea(1)) .and. (dy .ge. p_DobservationArea(2)) .and. &
                  (dx .le. p_DobservationArea(3)) .and. (dy .le. p_DobservationArea(4))) then

                ! Outer loop over the DOF's i=1..ndof on our current element,
                ! which corresponds to the (test) basis functions Phi_i:
                do idofe=1,p_rvectorData1%ndofTest
                
                  ! Fetch the contributions of the (test) basis functions Phi_i
                  ! into dbasI
                  dbasI  = p_DbasTest(idofe,DER_FUNC,icubp,iel)
                  
                  dylin1  = p_Dylin1(icubp,iel,DER_FUNC)
                  dylin2  = p_Dylin2(icubp,iel,DER_FUNC)
                  
                  ! Multiply the values of the basis functions
                  ! (1st derivatives) by the cubature weight and sum up
                  ! into the local vectors.
                  p_DlocalVector1(idofe,iel) = p_DlocalVector1(idofe,iel) + &
                      dweight * p_DcubWeight(icubp,iel) * &
                      ( dylin1 * dbasI ) ! (y',phi)

                  p_DlocalVector2(idofe,iel) = p_DlocalVector2(idofe,iel) + &
                      dweight * p_DcubWeight(icubp,iel) * &
                      ( dylin2 * dbasI ) ! (y',phi)
                      
                end do ! jdofe

              end if ! observation area

            end do ! icubp
          
          end do ! iel
        
        end if ! observation area type

      end select
      
    end select    
    
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>

  subroutine smva_assembleMatrix_primal (rmatrix,&
      rspatialOperatorAsm,rprimalSol,idofTime,bfull)

!<description>
  ! Assembles a linearised operator A'(.) which can be used for linear
  ! solvers for the primal equation.
!</description>

!<input>
  ! Definition of the operator. Must be a linearised operator.
  type(t_spacetimeOperatorAsm) :: rspatialOperatorAsm

  ! Structure that defines the nonlinearity in the
  ! primal equation.
  type(t_primalSpace), intent(inout) :: rprimalSol

  ! Number of the DOF in time which should be calculated into rmatrix.
  integer, intent(in) :: idofTime

  ! TRUE activates the full linearised operator (Newton).
  ! FALSE activates a partially linearised operator without the Newton part.
  logical, intent(in) :: bfull
!</input>

!<inputoutput>
  ! Block matrix receiving the result.
  type(t_matrixBlock) :: rmatrix
!</inputoutput>

!</subroutine>
    
    ! local variables
    real(DP) :: dtheta, dtstep, dtimeend, dtimestart
    
    ! Clear the output matrix
    call lsysbl_clearMatrix (rmatrix)

    ! Timestepping technique?
    select case (rspatialOperatorAsm%p_rtimeDiscr%ctype)
    
    ! ***********************************************************
    ! Standard Theta one-step scheme.
    ! ***********************************************************
    case (TDISCR_ONESTEPTHETA)
    
      ! Theta-scheme identifier
      dtheta = rspatialOperatorAsm%p_rtimeDiscr%dtheta

      ! Characteristics of the current timestep.
      call tdiscr_getTimestep(rspatialOperatorAsm%p_rtimeDiscr,idofTime-1,&
          dtimeend,dtstep,dtimestart)

      ! itag=0: old 1-step scheme.
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      select case (rspatialOperatorAsm%p_rtimeDiscr%itag)
      
      ! ***********************************************************
      ! itag=0: old/standard 1-step scheme.
      ! ***********************************************************
      case (0)

        call output_line("Old 1-step-scheme not implemented",&
            OU_CLASS_ERROR,OU_MODE_STD,"smva_assembleMatrix_primal")
        call sys_halt()

      ! ***********************************************************
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      ! ***********************************************************
      case (1)

        ! Which equation do we have?    
        select case (rspatialOperatorAsm%p_rphysics%cequation)

        ! *************************************************************
        ! Stokes/Navier Stokes.
        ! *************************************************************
        case (0,1)

          ! ===============================================
          ! Prepare the linear parts of the matrix.
          ! ===============================================
          
          ! -----------------------------------------
          ! Mass matrix for timestepping
          call smva_getMassMatrix (rspatialOperatorAsm,rmatrix,dtstep)
          
          ! -----------------------------------------
          ! Laplace -- if the viscosity is constant
          if (rspatialOperatorAsm%p_rphysics%cviscoModel .eq. 0) then
            call smva_getLaplaceMatrix (&
                rspatialOperatorAsm,rmatrix,rspatialOperatorAsm%p_rphysics%dnuConst)
          else
            call output_line("Nonconstant viscosity not supported.",&
                OU_CLASS_ERROR,OU_MODE_STD,"smva_assembleMatrix_primal")
            call sys_halt()
          end if
          
          ! -----------------------------------------
          ! B-matrices
          call smva_getBMatrix (rspatialOperatorAsm,rmatrix,1.0_DP)
          
          ! -----------------------------------------
          ! D-matrices
          call smva_getDMatrix (rspatialOperatorAsm,rmatrix,1.0_DP)

          ! -----------------------------------------
          ! EOJ-stabilisation
          !if (rspatialOperatorAsm%p_rsettingsDiscr%rstabilConvecPrimal%cupwind .eq. 4) then
          !  call smva_getEOJMatrix (rspatialOperatorAsm,rtempData%rmatrix,1.0_DP)
          !end if
            
          ! ===============================================
          ! Prepare the semilinear parts of the matrix.
          ! ===============================================
          
          ! The semilinear parts of the matrix can be set up with
          ! the block matrix assembly routines. Invoke them using
          ! a separate subroutine
          call smva_getSemilinMat_primalLin (rspatialOperatorAsm,idofTime,&
              rprimalSol,1.0_DP,bfull,rmatrix)

        end select ! Equation
      
      end select ! Timestep sub-scheme
      
    end select ! Timestep scheme
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine smva_assembleMatrix_dual (rmatrix,&
      rspatialOperatorAsm,rprimalSol,idofTime)

!<description>
  ! Assembles a linearised operator A'(.) which can be used for linear
  ! solvers for the dual equation.
!</description>

!<input>
  ! Definition of the operator. Must be a linearised operator.
  type(t_spacetimeOperatorAsm) :: rspatialOperatorAsm

  ! Structure that defines the nonlinearity in the
  ! primal equation.
  type(t_primalSpace), intent(inout) :: rprimalSol

  ! Number of the DOF in time which should be calculated into rmatrix.
  integer, intent(in) :: idofTime
!</input>

!<inputoutput>
  ! Block matrix receiving the result.
  type(t_matrixBlock) :: rmatrix
!</inputoutput>

!</subroutine>
    
    ! local variables
    real(DP) :: dtheta, dtstep, dtimeend, dtimestart
    
    ! Clear the output matrix
    call lsysbl_clearMatrix (rmatrix)

    ! Timestepping technique?
    select case (rspatialOperatorAsm%p_rtimeDiscr%ctype)
    
    ! ***********************************************************
    ! Standard Theta one-step scheme.
    ! ***********************************************************
    case (TDISCR_ONESTEPTHETA)
    
      ! Theta-scheme identifier
      dtheta = rspatialOperatorAsm%p_rtimeDiscr%dtheta

      ! Characteristics of the current timestep.
      call tdiscr_getTimestep(rspatialOperatorAsm%p_rtimeDiscr,idofTime-1,&
          dtimeend,dtstep,dtimestart)

      ! itag=0: old 1-step scheme.
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      select case (rspatialOperatorAsm%p_rtimeDiscr%itag)
      
      ! ***********************************************************
      ! itag=0: old/standard 1-step scheme.
      ! ***********************************************************
      case (0)

        call output_line("Old 1-step-scheme not implemented",&
            OU_CLASS_ERROR,OU_MODE_STD,"smva_getDef_dualLin")
        call sys_halt()

      ! ***********************************************************
      ! itag=1: new 1-step scheme, dual solutions inbetween primal solutions.
      ! ***********************************************************
      case (1)

        ! Which equation do we have?    
        select case (rspatialOperatorAsm%p_rphysics%cequation)

        ! *************************************************************
        ! Stokes/Navier Stokes.
        ! *************************************************************
        case (0,1)

          ! ===============================================
          ! Prepare the linear parts of the matrix.
          ! ===============================================

          ! Clear the temporary matrix
          call lsysbl_clearMatrix (rmatrix)
          
          ! -----------------------------------------
          ! Mass matrix for timestepping
          call smva_getMassMatrix (rspatialOperatorAsm,rmatrix,dtstep)
          
          ! -----------------------------------------
          ! Laplace -- if the viscosity is constant
          if (rspatialOperatorAsm%p_rphysics%cviscoModel .eq. 0) then
            call smva_getLaplaceMatrix (&
                rspatialOperatorAsm,rmatrix,rspatialOperatorAsm%p_rphysics%dnuConst)
          else
            call output_line("Nonconstant viscosity not supported.",&
                OU_CLASS_ERROR,OU_MODE_STD,"smva_getDef_dualLin")
            call sys_halt()
          end if
          
          ! -----------------------------------------
          ! B-matrices
          call smva_getBMatrix (rspatialOperatorAsm,rmatrix,1.0_DP)
          
          ! -----------------------------------------
          ! D-matrices
          call smva_getDMatrix (rspatialOperatorAsm,rmatrix,1.0_DP)

          ! -----------------------------------------
          ! EOJ-stabilisation
          !if (rspatialOperatorAsm%p_rsettingsDiscr%rstabilConvecPrimal%cupwind .eq. 4) then
          !  call smva_getEOJMatrix (rspatialOperatorAsm,rtempData%rmatrix,1.0_DP)
          !end if
            
          ! ===============================================
          ! Prepare the semilinear parts of the matrix.
          ! ===============================================
          
          ! The semilinear parts of the matrix can be set up with
          ! the block matrix assembly routines. Invoke them using
          ! a separate subroutine
          call smva_getSemilinMat_Dual (rspatialOperatorAsm,idofTime,&
              rprimalSol,1.0_DP,rmatrix)

        end select ! Equation
      
      end select ! Timestep sub-scheme
      
    end select ! Timestep scheme
  
  end subroutine

end module
