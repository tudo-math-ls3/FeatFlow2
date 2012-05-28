!##############################################################################
!# ****************************************************************************
!# <name> structurespostproc </name>
!# ****************************************************************************
!#
!# <purpose>
!# Underlying structures of the space-time optimal control solver.
!# </purpose>
!##############################################################################

module structurespostproc

  use fsystem
  use storage
  use basicgeometry
  use boundary
  use triangulation
  use cubature
  use paramlist
  use discretebc
  use discretefbc
  use fparser
  use linearsystemscalar
  use multilevelprojection
  
  use collection

  use spatialdiscretisation
  use timediscretisation
  use timescalehierarchy

  use meshhierarchy
  use fespacehierarchybase
  use fespacehierarchy
  use spacetimehierarchy

  use spacetimevectors
  use analyticsolution
  use spacetimeinterlevelprojection
  
  use assemblytemplates

  use constantsdiscretisation
  use structuresdiscretisation
  use structuresboundaryconditions
  use structuresdiscretisation
  
  implicit none
  
  private
  
!<types>

!<typeblock>

  ! Postprocessing structure. Whenever the postprocessing calculation routines are
  ! called, they fills this structure using the current solution vector (and other
  ! information if necessary). The information in this structure can then be used
  ! for GMV output e.g.
  type t_optcPostprocessing

    ! <!-- Input parameters -->
    
    ! Physics of the problem
    type(t_settings_physics), pointer :: p_rphysics
  
    ! Type of output file to generate from solutions.
    ! 0=disabled
    ! 1=GMV
    ! 2=AVS
    ! 3=Paraview (VTK)
    ! 4=Matlab
    integer :: ioutputUCD = 0
    
    ! Filename for UCD output.
    ! A timestep index '.0001','.0002',... is appended to this.
    character(len=SYS_STRLEN) :: sfilenameUCD = "./gmv/u"
    
    ! Output format for the final solution.
    ! =0: don't write
    ! =1: write out, use formatted output (default).
    ! =2: write out, use unformatted output.
    integer :: cwriteFinalSolution = 1

    ! Filename of a file sequence where the final solution is saved to.
    ! ="": Disable.
    character(len=SYS_STRLEN) :: sfinalSolutionFileName = ""

    ! Output format for the final control.
    ! =0: don't write
    ! =1: write out, use formatted output (default).
    ! =2: write out, use unformatted output.
    integer :: cwriteFinalControl = 1

    ! Filename of a file sequence where the final control is saved to.
    ! ="": Disable.
    character(len=SYS_STRLEN) :: sfinalControlFileName = ""
    
    ! Whether to calculate the values of the optimisation functional
    ! J(.) as well as ||y-y0|| etc. during the postprocessing of space-time vectors.
    integer :: icalcFunctionalValues = 0

    ! Whether to calculate the error to the analytic reference function
    ! ranalyticRefFunction during the postprocessing of space-time vectors.
    integer :: icalcError = 0

    ! Analytic reference function to be used during error calculation
    ! if icalcError > 0.
    type(t_anSolution) :: ranalyticRefFunction
    
    ! Whether to calculate drag/lift forces.
    integer :: icalcForces = 0
    
    ! Boundary component where to calculate body forces
    integer :: ibodyForcesBdComponent = 0
    
    ! 1st coefficient in the boundary integral of the drag coefficient.
    ! If this is commented out, 1/RE is assumed.
    real(DP) :: dbdForcesCoeff1 = 0.0_DP
    
    ! 2nd coefficient in the boundary integral of the drag coefficient.
    ! If this is commented out, 0.004 is assumed (corresonds to flow
    ! around cylinder with RE=1000: Umean=0.2, len=0.1
    ! -> coeff = Umean^2*len = 0.04*0.1 = 0.004 )
    real(DP) :: dbdForcesCoeff2 = 0.0_DP

    ! Whether to write drag/lift forces to hard disc.
    integer :: iwriteBodyForces = 0

    ! Filename for the body forces
    character(len=SYS_STRLEN) :: sfilenameBodyForces = ""

    ! Whether to calculate flux
    integer :: icalcFlux = 0

    ! Start/End coordinates of the lines along which to calculate the flux.
    real(DP), dimension(4) :: Dfluxline = (/0.0_DP,0.0_DP,0.0_DP,0.0_DP/)

    ! Whether to write the flux to a file
    integer :: iwriteFlux = 0

    ! Filename for the flux
    character(len=SYS_STRLEN) :: sfilenameFlux = ""

    ! Whether to calculate flux
    integer :: icalcKineticEnergy = 0

    ! Whether to write the flux to a file
    integer :: iwriteKineticEnergy = 0

    ! Filename for the flux
    character(len=SYS_STRLEN) :: sfilenameKineticEnergy = ""
    ! <!-- the following parameters are automatically maintained during a simulation -->
    
    ! Space that is available in rsolution. One of the CCSPACE_xxxx constants.
    integer :: cspace = CCSPACE_PRIMAL
    
    ! Underlying space discretisation of the primal space
    type(t_blockDiscretisation), pointer :: p_rspaceDiscrPrimal

    ! Underlying space discretisation
    type(t_blockDiscretisation), pointer :: p_rspaceDiscr

    ! Underlying time discretisation
    type(t_timeDiscretisation), pointer :: p_rtimeDiscr
    
    ! Space discretisation based on P1/Q1/P0/Q0 for visualisation output.
    type(t_blockDiscretisation) :: rspaceDiscrLinear

    ! Boundary conditions to use.
    type(t_optcBDC), pointer :: p_rboundaryConditions => null()
    
    ! Points coordinates where to evaluate point values. =NULL: Do not evaluate
    ! point values.
    real(DP), dimension(:,:), pointer :: p_DcoordsPointEval => null()
    
    ! Type of the point value to evaluate. Every entry corresponds to one
    ! point coordinate in p_DcoordsPointEval. The tuples are formed by
    ! (type,der) with
    !   type=1: primal x-velocity, =2: primal y-velocity, =3: primal pressure,
    !       =4: dual x-velocity, =5: dual y-velocity, =6: dual pressure,
    !   der =0: function value, =1: x-derivative, =2: y-derivative
    integer, dimension(:,:), pointer :: p_ItypePointEval => null()
    
    ! Whether or not to write the point values to a file.
    integer :: iwritePointValues = 0
    
    ! Filename for the point values if iwritePointValues <> 0.
    character(len=SYS_STRLEN) :: sfilenamePointValues = ""
        
  end type

!</typeblock>

  public :: t_optcPostprocessing

!</types>

contains

  ! ***************************************************************************
  
!<subroutine>

  subroutine struc_initPostprocParams (rparlist,ssection,rphysics,rpostproc)

!<description>
  ! Reads postprocessing parameters from a parameter list.
  ! Remark: Substructures are not initialised!
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section where the parameters can be found.
  character(len=*), intent(in) :: ssection

  ! Information about the physics.
  type(t_settings_physics), intent(in) :: rphysics
!</input>

!<output>
  ! Postprocessing structure, to set up with data.
  type(t_optcPostprocessing), intent(inout) :: rpostproc
!</output>

!</subroutine>

    ! local variables
    character(len=SYS_STRLEN) :: sstr, sparam
    integer :: npoints, i

    ! Read remaining parameters from the DAT file.
    !
    ! UCD export
    call parlst_getvalue_int (rparlist,ssection,&
        "ioutputUCD",rpostproc%ioutputUCD,0)

    call parlst_getvalue_string (rparlist,ssection,&
        "sfilenameUCD",rpostproc%sfilenameUCD,"",bdequote=.true.)

    ! Export of solution and control.
    call parlst_getvalue_int (rparlist,ssection,&
        "cwriteFinalSolution",rpostproc%cwriteFinalSolution,1)

    call parlst_getvalue_string (rparlist,ssection,&
        "sfinalSolutionFileName",rpostproc%sfinalSolutionFileName,&
        "",bdequote=.true.)

    call parlst_getvalue_int (rparlist,ssection,&
        "cwriteFinalControl",rpostproc%cwriteFinalControl,1)

    call parlst_getvalue_string (rparlist,ssection,&
        "sfinalControlFileName",rpostproc%sfinalControlFileName,&
        "",bdequote=.true.)

    ! function value calculation

    call parlst_getvalue_int (rparlist,ssection,&
        "icalcFunctionalValues",rpostproc%icalcFunctionalValues,0)

    ! Body forces

    call parlst_getvalue_int (rparlist,ssection,&
        "icalcForces",rpostproc%icalcForces,0)

    call parlst_getvalue_int (rparlist,ssection,&
        "ibodyForcesBdComponent",rpostproc%ibodyForcesBdComponent,2)

    call parlst_getvalue_double (rparlist,ssection,&
        "dbdForcesCoeff1",rpostproc%dbdForcesCoeff1,rphysics%dnuConst)

    call parlst_getvalue_double (rparlist,ssection,&
        "dbdForcesCoeff2",rpostproc%dbdForcesCoeff2,0.1_DP * 0.2_DP**2)

    call parlst_getvalue_int (rparlist,ssection,&
        "iwriteBodyForces",rpostproc%iwriteBodyForces,0)

    call parlst_getvalue_string (rparlist,ssection,&
        "sfilenameBodyForces",rpostproc%sfilenameBodyForces,&
        "",bdequote=.true.)

    ! Flux

    call parlst_getvalue_int (rparlist,ssection,&
        "icalcFlux",rpostproc%icalcFlux,0)

    call parlst_getvalue_int (rparlist,ssection,&
        "iwriteFlux",rpostproc%iwriteFlux,0)

    call parlst_getvalue_string (rparlist,ssection,&
        'sfilenameFlux',rpostproc%sfilenameFlux,"",bdequote=.true.)

    call parlst_getvalue_string (rparlist,ssection,&
        'dfluxline',sstr,"",bdequote=.true.)
        
    call parlst_getvalue_string (rparlist,ssection,&
        'dfluxline',sstr,"",bdequote=.true.)
    if (sstr .ne. "") then
      ! Read the start/end coordinates
      read(sstr,*) rpostproc%Dfluxline(1),rpostproc%Dfluxline(2),&
          rpostproc%Dfluxline(3),rpostproc%Dfluxline(4)
    end if

    ! internal Energy

    call parlst_getvalue_int (rparlist,ssection,&
        "icalcKineticEnergy",rpostproc%icalcKineticEnergy,1)

    call parlst_getvalue_int (rparlist,ssection,&
        "iwriteKineticEnergy",rpostproc%iwriteKineticEnergy,0)

    call parlst_getvalue_string (rparlist,ssection,&
        'sfilenameKineticEnergy',rpostproc%sfilenameKineticEnergy,"",bdequote=.true.)

    ! Error calculation

    call parlst_getvalue_int (rparlist,ssection,&
        "icalcError",rpostproc%icalcError,0)

    call parlst_getvalue_string (rparlist,ssection,&
        "ssectionReferenceFunction",sstr,"",bdequote=.true.)
    
    if (sstr .eq. "") &
        rpostproc%icalcError = 0
    
    ! Init the points to evaluate
    npoints = parlst_querysubstrings (rparlist, ssection, &
        "CEVALUATEPOINTVALUES")
 
    if (npoints .gt. 0) then
      allocate (rpostproc%p_DcoordsPointEval(NDIM2D,npoints))
      allocate (rpostproc%p_ItypePointEval(NDIM2D,npoints))
    
      ! Read the points
      do i=1,npoints
        call parlst_getvalue_string (rparlist, ssection, &
            "CEVALUATEPOINTVALUES", sparam, "", i)
        read (sparam,*) rpostproc%p_DcoordsPointEval(1,i),rpostproc%p_DcoordsPointEval(2,i),&
            rpostproc%p_ItypePointEval(1,i),rpostproc%p_ItypePointEval(2,i)
      end do
    end if

    call parlst_getvalue_int (rparlist,ssection,&
        "iwritePointValues",rpostproc%iwritePointValues,0)

    call parlst_getvalue_string (rparlist,ssection,&
        'sfilenamePointValues',rpostproc%sfilenamePointValues,"",bdequote=.true.)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine struc_donePostprocParams (rpostproc)
  
!<description>
  ! Clean up postprocessing parameters.
!</description>
  
!<inputoutput>
  ! Postprocessing parameters to be cleaned up.
  type(t_optcPostprocessing), intent(inout) :: rpostproc
!</inputoutput>
  
!</subroutine>

    ! Release allocated parameter arrays.
    if (associated(rpostproc%p_DcoordsPointEval)) deallocate(rpostproc%p_DcoordsPointEval)
    if (associated(rpostproc%p_ItypePointEval)) deallocate(rpostproc%p_ItypePointEval)

  end subroutine
  
end module
