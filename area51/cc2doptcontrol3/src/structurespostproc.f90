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
  
  implicit none
  
  private
  
  public :: t_optcPostprocessing
  
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

!</types>

end module
