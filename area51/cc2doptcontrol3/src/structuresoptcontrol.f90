!##############################################################################
!# ****************************************************************************
!# <name> structuresoptcontrol </name>
!# ****************************************************************************
!#
!# <purpose>
!# Underlying structures of the space-time optimal control solver.
!#
!# Routines in this module:
!#
!# 1.) soptc_initParOptControl
!#     -> Reads the parameters for the optimal control from the parameter list
!#        and stores them in the structure.
!#
!# 2.) soptc_doneParOptControl
!#     -> Releases data allocated in soptc_initParOptControl.
!#
!# 3.) soptc_initStepPathFollower
!#     -> Adapts the parameters in the optimal control structure
!#        according to path following information.
!# </purpose>
!##############################################################################

module structuresoptcontrol

  use fsystem
  use storage
  use genoutput
  use boundary
  use triangulation
  use paramlist
  use discretebc
  use discretefbc
  use fparser
  
  use collection
  use paramlist
  use linearsystemblock

  use analyticsolution
  
  use timediscretisation
  
  use structuresdiscretisation
  
  implicit none
  
  private
  
!<types>

!<typeblock>

!  ! Type block encapsuling constraints, only in space.
!  type t_optcconstraintsSpace
!  
!    ! Type of constraints to apply to the control u.
!    ! =0: No constraints.
!    ! =1: Constant constraints on u active:
!    !     ddistVelUmin1 <= u_1 <= ddistVelUmax1, ddistVelUmin2 <= u_2 <= ddistVelUmax2
!    !     Implementation by DOF.
!    ! =2: Constant constraints on u active:
!    !     ddistVelUmin1 <= u_1 <= ddistVelUmax1, ddistVelUmin2 <= u_2 <= ddistVelUmax2
!    !     Implemented by cubature point.
!    ! =3: Constant constraints on u active: 
!    !     ddistVelUmin1 <= u_1 <= ddistVelUmax1, ddistVelUmin2 <= u_2 <= ddistVelUmax2.
!    !     Implemented by DOF. Newton matrix is approximative.
!    ! =4: Constant constraints on u active: 
!    !     ddistVelUmin1 <= u_1 <= ddistVelUmax1, ddistVelUmin2 <= u_2 <= ddistVelUmax2.
!    !     Implemented by cubature point (more exact). Adaptive integration.
!    integer :: ccontrolConstraints = 0
!
!    ! Type of definition of the constraints if ccontrolConstraints <> 0.
!    ! =0: constants specified in ddistVelUmin1/2, ddistVelUmax1/2.
!    ! =1: analytical functions defined in p_rumin1/2, p_rumax1/2.
!    integer :: ccontrolConstraintsType = 0
!
!    ! Constraints on u_1
!    real(DP) :: ddistVelUmin1 = -1.0E10
!    real(DP) :: ddistVelUmax1 = 1.0E10
!    
!    ! Constraints in u_2
!    real(DP) :: ddistVelUmin2 = -1.0E10
!    real(DP) :: ddistVelUmax2 = 1.0E10
!
!    ! Analytical constraints for u_1
!    type(t_anSolution), pointer :: p_rumin1 => null()
!    type(t_anSolution), pointer :: p_rumax1 => null()
!    
!    ! Analytical constraints for u_2
!    type(t_anSolution), pointer :: p_rumin2 => null()
!    type(t_anSolution), pointer :: p_rumax2 => null()
!    
!    ! Discrete constraints for u_1, u_2
!    type(t_vectorBlock), pointer :: p_rvectorumin => null()
!    type(t_vectorBlock), pointer :: p_rvectorumax => null()
!
!
!    ! Type of constraints to apply to the state y.
!    ! =0: No constraints.
!    ! =1: Constant constraints on u active:
!    !     dymin1 <= y_1 <= dymax1, dymin2 <= y_2 <= dymax2
!    !     Implementation by DOF.
!    ! =2: Constant constraints on y active:
!    !     dymin1 <= y_1 <= dymax1, dymin2 <= y_2 <= dymax2
!    !     Implemented by cubature point.
!    integer :: cstateConstraints = 0
!
!    ! Type of definition of the constraints if ccontrolConstraints <> 0.
!    ! =0: constants specified in ddistVelUmin1/2, ddistVelUmax1/2.
!    ! =1: analytical functions defined in p_rumin1/2, p_rumax1/2.
!    integer :: cstateConstraintsType = 0
!
!    ! Regularisation parameter for the Moreau-Yosida regularisation
!    real(DP) :: dstateConstrReg = 1.0_DP
!
!    ! Constraints on u_1
!    real(DP) :: dymin1 = -1.0E10
!    real(DP) :: dymax1 = 1.0E10
!    
!    ! Constraints in u_2
!    real(DP) :: dymin2 = -1.0E10
!    real(DP) :: dymax2 = 1.0E10
!
!    ! Analytical constraints for u_1
!    type(t_anSolution), pointer :: p_rymin1 => null()
!    type(t_anSolution), pointer :: p_rymax1 => null()
!    
!    ! Analytical constraints for u_2
!    type(t_anSolution), pointer :: p_rymin2 => null()
!    type(t_anSolution), pointer :: p_rymax2 => null()
!    
!    ! Discrete constraints for u_1, u_2
!    type(t_vectorBlock), pointer :: p_rvectorymin => null()
!    type(t_vectorBlock), pointer :: p_rvectorymax => null()
!    
!    
!    ! Time used for the evaluation of the analytical control constraints
!    ! p_rumin1/2, p_rumax1/2.
!    real(DP) :: dstateConstrTime = 0.0_DP
!
!    ! Time used for the evaluation of the analytical state constraints
!    ! p_rymin1/2, p_rymax1/2.
!    real(DP) :: dcontrolConstrTime = 0.0_DP
!    
!  end type
!  
!!</types>
!
!  public :: t_optcconstraintsSpace

!<types>

  ! Type block encapsuling constraints, in space/time.
  type t_optcconstraintsSpaceTime
  
    ! Type of constraints to apply to the distributed velocity control.
    ! =0: No constraints.
    ! =1: Constant constraints on u active:
    !     ddistVelUmin1 <= u_1 <= ddistVelUmax1, ddistVelUmin2 <= u_2 <= ddistVelUmax2
    !     Implementation by DOF.
    ! =2: Constant constraints on u active:
    !     ddistVelUmin1 <= u_1 <= ddistVelUmax1, ddistVelUmin2 <= u_2 <= ddistVelUmax2
    !     Implemented by cubature point.
    integer :: cdistVelConstraints = 0

    ! Type of definition of the constraints if ccontrolConstraints <> 0.
    ! =0: constants specified in ddistVelUmin1/2, ddistVelUmax1/2.
    ! =1: analytical functions defined in p_rumin1/2, p_rumax1/2.
    integer :: cdistVelConstType = 0

    ! Constraints on u_1
    real(DP) :: ddistVelUmin1 = -1.0E10
    real(DP) :: ddistVelUmax1 = 1.0E10

    ! Analytical constraints for u_1.
    type(t_anSolution), pointer :: p_rdistVelUmin1 => null()
    type(t_anSolution), pointer :: p_rdistVelUmax1 => null()

    ! Constraints in u_2
    real(DP) :: ddistVelUmin2 = -1.0E10
    real(DP) :: ddistVelUmax2 = 1.0E10
  
    ! Analytical constraints for u_2
    type(t_anSolution), pointer :: p_rdistVelUmin2 => null()
    type(t_anSolution), pointer :: p_rdistVelUmax2 => null()
  
    ! Type of constraints to apply to the state y.
    ! =0: No constraints.
    ! =1: Constant constraints on u active:
    !     dymin1 <= y_1 <= dymax1, dymin2 <= y_2 <= dymax2
    !     Implementation by DOF.
    ! =2: Constant constraints on y active:
    !     dymin1 <= y_1 <= dymax1, dymin2 <= y_2 <= dymax2
    !     Implemented by cubature point.
    integer :: cstateConstraints = 0

    ! Type of definition of the constraints if ccontrolConstraints <> 0.
    ! =0: constants specified in ddistVelUmin1/2, ddistVelUmax1/2.
    ! =1: analytical functions defined in p_rumin1/2, p_rumax1/2.
    integer :: cstateConstraintsType = 0

    ! Regularisation parameter for the Moreau-Yosida regularisation
    real(DP) :: dstateConstrReg = 1.0_DP

    ! Constraints on u_1
    real(DP) :: dymin1 = -1.0E10
    real(DP) :: dymax1 = 1.0E10
    
    ! Constraints in u_2
    real(DP) :: dymin2 = -1.0E10
    real(DP) :: dymax2 = 1.0E10

    ! Analytical constraints for u_1
    type(t_anSolution), pointer :: p_rymin1 => null()
    type(t_anSolution), pointer :: p_rymax1 => null()
    
    ! Analytical constraints for u_2
    type(t_anSolution), pointer :: p_rymin2 => null()
    type(t_anSolution), pointer :: p_rymax2 => null()
    
  end type

!</typeblock>

  public :: t_optcconstraintsSpaceTime

!<typeblock>

  ! This structure collects a set of path following data.
  ! During the nonlinear iteration, this data is automatically
  ! adapted from nonlinear to nonlinear step.
  type t_optcpathfollowingdata
    
    ! List of $\alpha$ parameter of the optimal control functional parameters. 
    ! If this is not =>NULL, a path following algorithm will be applied. 
    ! Every nonlinear step, the next value is taken.
    real(DP), dimension(:), pointer :: p_DalphaC => null()
    
    ! List of $\beta$ parameter of the optimal control functional parameters. 
    ! If this is not =>NULL, a path following algorithm will be applied. 
    ! Every nonlinear step, the next value is taken.
    real(DP), dimension(:), pointer :: p_DbetaC => null()
    
    ! List of $\gamma$ parameter of the optimal control functional parameters. 
    ! If this is not =>NULL, a path following algorithm will be applied. 
    ! Every nonlinear step, the next value is taken.
    real(DP), dimension(:), pointer :: p_DgammaC => null()
  
  end type

!</typeblock>

  public :: t_optcpathfollowingdata

!<typeblock>

  ! This block saves parameters of the optimal control problem
  type t_settings_optcontrol
  
    ! $\alpha$ parameter of the optimal control functional
    real(DP) :: dalphaC = 1.0_DP
    
    ! $\beta$ parameter of the optimal control functional via boundary control
    real(DP) :: dbetaC = 1.0_DP
    
    ! Penalty parameter for the dirichlet boundary control
    real(DP) :: ddirichletBCPenalty = 100.0_DP
    
    ! $\gamma$ parameter of the nonstationary optimal control functional
    real(DP) :: dgammaC = 0.0_DP
  
    ! Formulation of the Space-time problem.
    ! =0: usual formulation as specified in the DFG applicance
    ! =1: Formulation for the generation of reference results from papers
    ! The two formulations differ in a "-"-sign in front of the dual velocity.
    integer :: ispaceTimeFormulation = 0
  
    ! Whether to treat the convection explicitly or implicitly.
    ! =0: Treat the convection implicitely.
    ! =1: Treat the convection explicitely. (This is the formulation of the paper
    !     of Baerwolff and Hinze!)
    integer :: iconvectionExplicit = 0

    !<!-- --------------- -->
    !<!-- TARGET FUNCTION -->
    !<!-- --------------- -->

    ! Analytic solution defining the target function
    type(t_anSolution) :: rtargetFunction
    
    !<!-- ----------- -->
    !<!-- CONSTRAINTS -->
    !<!-- ----------- -->
  
    ! Possible constraints in the problem.
    type(t_optcconstraintsSpaceTime) :: rconstraints

    !<!-- ---------------- -->
    !<!-- OBSERVATION AREA -->
    !<!-- ---------------- -->
    
    ! Observation area. If this points to NULL(), the whole
    ! domain is observed. Otherwise, this is an array with four entries
    ! in the form
    !   p_DobservationArea = (x1 y1 x2 y2)
    ! specifying the area to be observed.
    real(DP), dimension(:), pointer :: p_DobservationArea => null()

    !<!-- ---------------- -->
    !<!-- OBSERVATION AREA -->
    !<!-- ---------------- -->

    ! Contains data for path following strategies
    type(t_optcpathfollowingdata) :: rpathfollowingdata
  end type

!</typeblock>

  public :: t_settings_optcontrol

!</types>

  public :: soptc_initParOptControl
  public :: soptc_doneParOptControl
  public :: soptc_initStepPathFollower

contains

  ! ***************************************************************************

!<subroutine>

  subroutine sys_getParameter (sstring,stoken,istart)
  
!<description>
  ! Gets a parameter out of a list of parameters.
!</description>

!<input>
  ! A string with parameters, e.g. "0.1 0.2 0.3",...
  character(len=*), intent(in) :: sstring
  
  ! String that receives the parameter.
  character(len=*), intent(out) :: stoken
!</input>

!<inputoutput>
  ! On input: Must be set to =1 for the first call.
  !   Position from where to search for the next parameter
  ! On output:
  !   =0, if this was the last parameter
  !   Otherwise: position in sstring where the next parameter starts.
  integer, intent(inout) :: istart
!</inputoutput>

!</subroutine>

    integer :: i,slen
    
    ! If istart=0, do not do anything. There is nothing in the string
    stoken = ""
    if (istart .le. 0) return
    
    ! If the string length is =0, finish.
    slen = len_trim(sstring)
    if (slen .eq. 0) return
    
    ! Skip whitespaces.
    do while (istart .le. slen)
      if (sstring(istart:istart) .eq. " ") then
        istart = istart + 1
      else
        exit
      end if
    end do
    
    ! End of the string?
    i = istart
    do while (i .le. slen)
      if (sstring(i:i) .ne. " ") then
        i = i + 1
      else
        exit
      end if
    end do
    
    ! Copy
    if (istart .le. slen) then
      stoken = sstring(istart:i)
    end if
    
    ! Put istart behind the parameter or set it to 0.
    if (i .ge. slen) then
      ! No more parameters
      istart = 0
    else
      ! Put istart behind the parameter -- on the whitespace which
      ! is skipped in the next loop.
      istart = i+1
    end if
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sys_countParameters (sstring,ntokens)
  
!<description>
  ! Returns the number of parameters in sstring
!</description>

!<input>
  ! A string with parameters, e.g. "0.1 0.2 0.3",...
  character(len=*), intent(in) :: sstring
!</input>

!<output>
  ! Number of parameters in sstring
  integer, intent(out) :: ntokens
!</output>

!</subroutine>

    integer :: slen, istart
    
    ! If the string length is =0, finish.
    ntokens = 0
    
    slen = len_trim(sstring)
    if (slen .eq. 0) return
    
    ! Find all substrings  
    istart = 1
    
    do
      ! Cancel if nothing left.
      if (istart .gt. slen) exit

      ! Skip whitespaces.
      do while (istart .le. slen)
        if (sstring(istart:istart) .eq. " ") then
          istart = istart + 1
        else
          exit
        end if
      end do
      
      ! Cancel if we reached the string end
      if (istart .gt. slen) exit
      
      ! End of the string?
      do while (istart .le. slen)
        if (sstring(istart:istart) .ne. " ") then
          istart = istart + 1
        else
          exit
        end if
      end do
      
      ! One more
      ntokens = ntokens + 1
      
    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine soptc_initParOptControl (rparlist,ssectionOptC,roptcontrol)
  
!<description>
  ! Reads the parameters for the optimal control from the parameter list
  ! and stores them in the structure.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section where the parameters of the optimal control can be found.
  character(len=*), intent(in) :: ssectionOptC
!</input>

!<inputoutput>
  ! Structure defining the parameters for the optimal control.
  type(t_settings_optcontrol), intent(inout) :: roptcontrol
!</inputoutput>

!</subroutine>
    
    character(len=SYS_STRLEN) :: sstring,stoken
    integer :: ntokens,istart,itoken

    ! Alpha/Gamma parameters.
    !
    ! These parameters may have two types: either it is a single
    ! parameter or a list of parameters. In the second case,
    ! this defines a path following strategy, i.e., the parameters
    ! change according to the list in every nonlinear iteration.
    call parlst_getvalue_string (rparlist,ssectionOptC,&
        "dalphaC",sstring,"-1.0",bdequote=.true.)

    ! Is this a list?
    call sys_countParameters (sstring,ntokens)
    
    if (ntokens .eq. 1) then
    
      read(sstring,*) roptcontrol%dalphaC
      
    else if (ntokens .gt. 1) then
      ! Read the list
      roptcontrol%dalphaC = 0.0_DP
      allocate(roptcontrol%rpathfollowingdata%p_DalphaC(ntokens))
      
      ! Loop through all tokens and save them
      istart = 1
      do itoken = 1,ntokens
        call sys_getParameter (sstring,stoken,istart)
        read(stoken,*) roptcontrol%rpathfollowingdata%p_DalphaC(itoken)
      end do
    end if
        
    call parlst_getvalue_string (rparlist,ssectionOptC,&
        "dbetaC",sstring,"-1.0",bdequote=.true.)

    ! Is this a list?
    call sys_countParameters (sstring,ntokens)

    if (ntokens .eq. 1) then
    
      read(sstring,*) roptcontrol%dbetaC
      
    else if (ntokens .gt. 1) then
      ! Read the list
      roptcontrol%dbetaC = 0.0_DP
      allocate(roptcontrol%rpathfollowingdata%p_DbetaC(ntokens))
      
      ! Loop through all tokens and save them
      istart = 1
      do itoken = 1,ntokens
        call sys_getParameter (sstring,stoken,istart)
        read(stoken,*) roptcontrol%rpathfollowingdata%p_DbetaC(itoken)
      end do
    end if

    call parlst_getvalue_string (rparlist,ssectionOptC,&
        "dgammaC",sstring,"0.0",bdequote=.true.)

    ! Is this a list?
    call sys_countParameters (sstring,ntokens)

    if (ntokens .eq. 1) then
    
      read(sstring,*) roptcontrol%dgammaC

    else if (ntokens .gt. 1) then
      ! Read the list
      roptcontrol%dgammaC = 0.0_DP
      allocate(roptcontrol%rpathfollowingdata%p_DgammaC(ntokens))
      
      ! Loop through all tokens and save them
      istart = 1
      do itoken = 1,ntokens
        call sys_getParameter (sstring,stoken,istart)
        read(stoken,*) roptcontrol%rpathfollowingdata%p_DgammaC(itoken)
      end do
    end if

            
    call parlst_getvalue_double (rparlist,ssectionOptC,&
        "ddirichletBCPenalty",roptcontrol%ddirichletBCPenalty,100.0_DP)
        
    ! Type of the formulation
    call parlst_getvalue_int (rparlist,ssectionOptC,&
        "ispaceTimeFormulation",roptcontrol%ispaceTimeFormulation,0)
    
    call parlst_getvalue_int (rparlist,ssectionOptC,&
        "iconvectionExplicit",roptcontrol%iconvectionExplicit,0)

    ! Parameters defining the control constraints
    call parlst_getvalue_int (rparlist,ssectionOptC,&
        "ccontrolConstraints",roptcontrol%rconstraints%ccontrolConstraints,0)

    roptcontrol%rconstraints%cdistVelConstType = 0

    call parlst_getvalue_double (rparlist,ssectionOptC,&
        "ddistVelUmin1",roptcontrol%rconstraints%ddistVelUmin1,-1.0E10_DP)

    call parlst_getvalue_double (rparlist,ssectionOptC,&
        "ddistVelUmax1",roptcontrol%rconstraints%ddistVelUmax1,1.0E10_DP)

    call parlst_getvalue_double (rparlist,ssectionOptC,&
        "ddistVelUmin2",roptcontrol%rconstraints%ddistVelUmin2,-1.0E10_DP)

    call parlst_getvalue_double (rparlist,ssectionOptC,&
        "ddistVelUmax2",roptcontrol%rconstraints%ddistVelUmax2,1.0E10_DP)

    ! Parameters defining the state constraints
    call parlst_getvalue_int (rparlist,ssectionOptC,&
        "cstateConstraints",roptcontrol%rconstraints%cstateConstraints,0)

    roptcontrol%rconstraints%cstateConstraintsType = 0

    call parlst_getvalue_double (rparlist,ssectionOptC,&
        "dymin1",roptcontrol%rconstraints%dymin1,-1.0E10_DP)

    call parlst_getvalue_double (rparlist,ssectionOptC,&
        "dymax1",roptcontrol%rconstraints%dymax1,1.0E10_DP)

    call parlst_getvalue_double (rparlist,ssectionOptC,&
        "dymin2",roptcontrol%rconstraints%dymin2,-1.0E10_DP)

    call parlst_getvalue_double (rparlist,ssectionOptC,&
        "dymax2",roptcontrol%rconstraints%dymax2,1.0E10_DP)

    call parlst_getvalue_double (rparlist,ssectionOptC,&
        "dstateConstrReg",roptcontrol%rconstraints%dstateConstrReg,1.0_DP)
        
    ! Observation area
    call parlst_getvalue_string (rparlist,ssectionOptC,&
        "DobservationArea",sstring,"")
        
    if (sstring .ne. "") then
      ! Read the observation area. This is a box, the parameters
      ! have the format "x1, y1, x2, y2".
      allocate (roptcontrol%p_DobservationArea(4))
      read (sstring,*) &
          roptcontrol%p_DobservationArea(1), &
          roptcontrol%p_DobservationArea(2), &
          roptcontrol%p_DobservationArea(3), &
          roptcontrol%p_DobservationArea(4)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine soptc_doneParOptControl (roptcontrol)
  
!<description>
  ! Cleans up information in the optimal control structure.
!</description>

!<inputoutput>
  ! Structure defining the parameters for the optimal control.
  type(t_settings_optcontrol), intent(inout) :: roptcontrol
!</inputoutput>

!</subroutine>

    ! Release the observation area
    if (associated(roptcontrol%p_DobservationArea)) then
      deallocate(roptcontrol%p_DobservationArea)
    end if
    
    ! Release path following data
    if (associated(roptcontrol%rpathfollowingdata%p_DalphaC)) then
      deallocate(roptcontrol%rpathfollowingdata%p_DalphaC)
    end if

    if (associated(roptcontrol%rpathfollowingdata%p_DbetaC)) then
      deallocate(roptcontrol%rpathfollowingdata%p_DbetaC)
    end if

    if (associated(roptcontrol%rpathfollowingdata%p_DgammaC)) then
      deallocate(roptcontrol%rpathfollowingdata%p_DgammaC)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine soptc_initStepPathFollower (roptcontrol,istep,bprint)
  
!<description>
  ! Adapts the parameters in the optimal control structure
  ! according to path following information.
  ! The step istep defines the current setting of the parameters
  ! which are adapted via path following methods.
!</description>

!<input>
  ! Current step in the path following strategy
  integer, intent(in) :: istep
  
  ! Whether or not to print information about the path following parameters.
  logical, intent(in) :: bprint
!</input>

!<inputoutput>
  ! Structure defining the parameters for the optimal control.
  ! The parameters in this structure which are influenced
  ! via the path following data are modified according to
  ! step istep.
  type(t_settings_optcontrol), intent(inout) :: roptcontrol
!</inputoutput>

!</subroutine>

    ! Change the data in roptcontrol.
    ! istep defines the array position of the element to take.
    if (associated(roptcontrol%rpathfollowingdata%p_DalphaC)) then
      roptControl%dalphaC = roptControl%rpathfollowingdata%p_DalphaC(&
              min(size(roptcontrol%rpathfollowingdata%p_DalphaC),istep))
      if (bprint) then
        call output_line ("Path-Follower: Parameter ALPHA = "//&
            trim(sys_sdEL(roptcontrol%dalphaC,10)))
      end if
    end if

    if (associated(roptcontrol%rpathfollowingdata%p_DbetaC)) then
      roptcontrol%dbetaC = roptControl%rpathfollowingdata%p_DbetaC(&
              min(size(roptcontrol%rpathfollowingdata%p_DbetaC),istep))
      if (bprint) then
        call output_line ("Path-Follower: Parameter BETA  = "//&
            trim(sys_sdEL(roptcontrol%dbetaC,10)))
      end if
    end if

    if (associated(roptcontrol%rpathfollowingdata%p_DgammaC)) then
      roptcontrol%dgammaC = roptControl%rpathfollowingdata%p_DgammaC(&
              min(size(roptcontrol%rpathfollowingdata%p_DgammaC),istep))
      if (bprint) then
        call output_line ("Path-Follower: Parameter GAMMA = "//&
            trim(sys_sdEL(roptcontrol%dgammaC,10)))
      end if
    end if

  end subroutine

end module
