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
!    !     rconstraintsDistCtrl%dmin1 <= u_1 <= rconstraintsDistCtrl%dmax1, rconstraintsDistCtrl%dmin2 <= u_2 <= rconstraintsDistCtrl%dmax2
!    !     Implementation by DOF.
!    ! =2: Constant constraints on u active:
!    !     rconstraintsDistCtrl%dmin1 <= u_1 <= rconstraintsDistCtrl%dmax1, rconstraintsDistCtrl%dmin2 <= u_2 <= rconstraintsDistCtrl%dmax2
!    !     Implemented by cubature point.
!    ! =3: Constant constraints on u active: 
!    !     rconstraintsDistCtrl%dmin1 <= u_1 <= rconstraintsDistCtrl%dmax1, rconstraintsDistCtrl%dmin2 <= u_2 <= rconstraintsDistCtrl%dmax2.
!    !     Implemented by DOF. Newton matrix is approximative.
!    ! =4: Constant constraints on u active: 
!    !     rconstraintsDistCtrl%dmin1 <= u_1 <= rconstraintsDistCtrl%dmax1, rconstraintsDistCtrl%dmin2 <= u_2 <= rconstraintsDistCtrl%dmax2.
!    !     Implemented by cubature point (more exact). Adaptive integration.
!    integer :: ccontrolConstraints = 0
!
!    ! Type of definition of the constraints if ccontrolConstraints <> 0.
!    ! =0: constants specified in rconstraintsDistCtrl%dmin1/2, rconstraintsDistCtrl%dmax1/2.
!    ! =1: analytical functions defined in p_rumin1/2, p_rumax1/2.
!    integer :: ccontrolConstraintsType = 0
!
!    ! Constraints on u_1
!    real(DP) :: rconstraintsDistCtrl%dmin1 = -1.0E10
!    real(DP) :: rconstraintsDistCtrl%dmax1 = 1.0E10
!    
!    ! Constraints in u_2
!    real(DP) :: rconstraintsDistCtrl%dmin2 = -1.0E10
!    real(DP) :: rconstraintsDistCtrl%dmax2 = 1.0E10
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
!    !     rconstraintsState%dmin1 <= y_1 <= rconstraintsState%dmax1, rconstraintsState%dmin2 <= y_2 <= rconstraintsState%dmax2
!    !     Implementation by DOF.
!    ! =2: Constant constraints on y active:
!    !     rconstraintsState%dmin1 <= y_1 <= rconstraintsState%dmax1, rconstraintsState%dmin2 <= y_2 <= rconstraintsState%dmax2
!    !     Implemented by cubature point.
!    integer :: rconstraintsState%cconstrType%cconstraints = 0
!
!    ! Type of definition of the constraints if ccontrolConstraints <> 0.
!    ! =0: constants specified in rconstraintsDistCtrl%dmin1/2, rconstraintsDistCtrl%dmax1/2.
!    ! =1: analytical functions defined in p_rumin1/2, p_rumax1/2.
!    integer :: rconstraintsState%cconstrType%cconstraintsType = 0
!
!    ! Regularisation parameter for the Moreau-Yosida regularisation
!    real(DP) :: rconstraintsState%dregularisation = 1.0_DP
!
!    ! Constraints on u_1
!    real(DP) :: rconstraintsState%dmin1 = -1.0E10
!    real(DP) :: rconstraintsState%dmax1 = 1.0E10
!    
!    ! Constraints in u_2
!    real(DP) :: rconstraintsState%dmin2 = -1.0E10
!    real(DP) :: rconstraintsState%dmax2 = 1.0E10
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

!<typeblock>

  ! Defines a set of contraints to any function
  type t_optcConstraintsDef
  
    ! Type of constraints to apply.
    ! =0: No constraints.
    ! =1: Constant constraints active:
    !     dmin1 <= u_1 <= dmax1, dmin2 <= u_2 <= dmax2
    !     Implementation by DOF.
    ! =2: Constant constraints active:
    !     dmin1 <= u_1 <= dmax1, dmin2 <= u_2 <= dmax2
    !     Implemented by cubature point.
    integer :: cconstraints = 0

    ! Type of definition of the constraints if ccontrolConstraints <> 0.
    ! =0: constants specified in dmin1/2, dmax1/2.
    ! =1: analytical functions defined in p_rmin1/2, p_rmax1/2.
    integer :: cconstrType = 0

    ! Constraints on u_1
    real(DP) :: dmin1 = -1.0E10_DP
    real(DP) :: dmax1 = 1.0E10_DP
    
    ! Section defining the minimum/maximum for the 1st component
    character(LEN=SYS_STRLEN) :: ssectionMin1 = ""
    character(LEN=SYS_STRLEN) :: ssectionMax1 = ""

    ! Analytical constraints for u_1.
    type(t_anSolution), pointer :: p_rmin1 => null()
    type(t_anSolution), pointer :: p_rmax1 => null()

    ! Constraints in u_2
    real(DP) :: dmin2 = -1.0E10_DP
    real(DP) :: dmax2 = 1.0E10_DP

    ! Section defining the minimum/maximum for the 1st component
    character(LEN=SYS_STRLEN) :: ssectionMin2 = ""
    character(LEN=SYS_STRLEN) :: ssectionMax2 = ""
  
    ! Analytical constraints for u_2
    type(t_anSolution), pointer :: p_rmin2 => null()
    type(t_anSolution), pointer :: p_rmax2 => null()
    
    ! Regularisation parameter. Not used for all types of constraints
    real(DP) :: dregularisation = 0.0_DP

  end type

!</typeblock>

  public :: t_optcConstraintsDef

!<typeblock>

  ! Type block encapsuling constraints, in space/time.
  type t_optcconstraintsSpaceTime
  
    ! Constraints in the distributed control
    type(t_optcConstraintsDef) :: rconstraintsDistCtrl
    
    ! Constraints in the L2 boundary control
    type(t_optcConstraintsDef) :: rconstraintsL2BdC

    ! Constraints in the $H^{1/2}$ boundary control
    type(t_optcConstraintsDef) :: rconstraintsH12BdC
    
    ! Constraints in the state
    type(t_optcConstraintsDef) :: rconstraintsState
    
  end type

!</typeblock>

  public :: t_optcconstraintsSpaceTime

!<typeblock>

  ! This structure collects a set of path following data.
  ! During the nonlinear iteration, this data is automatically
  ! adapted from nonlinear to nonlinear step.
  type t_optcpathfollowingdata
    
    ! List of $\alpha$ parameter of the optimal control functional parameters.
    ! Distributed control.
    ! If this is not =>NULL, a path following algorithm will be applied. 
    ! Every nonlinear step, the next value is taken.
    real(DP), dimension(:), pointer :: p_dalphaDistC => null()
    
    ! List of $\alpha$ parameter of the optimal control functional parameters. 
    ! $L_2$ boundary control.
    ! If this is not =>NULL, a path following algorithm will be applied. 
    ! Every nonlinear step, the next value is taken.
    real(DP), dimension(:), pointer :: p_dalphaL2BdC => null()

    ! List of $\alpha$ parameter of the optimal control functional parameters. 
    ! $H^{1/2}$ boundary control.
    ! If this is not =>NULL, a path following algorithm will be applied. 
    ! Every nonlinear step, the next value is taken.
    real(DP), dimension(:), pointer :: p_dalphaH12BdC => null()
    
    ! List of $\gamma$ parameter of the optimal control functional parameters. 
    ! If this is not =>NULL, a path following algorithm will be applied. 
    ! Every nonlinear step, the next value is taken.
    real(DP), dimension(:), pointer :: p_dalphaEndTimeC => null()
  
  end type

!</typeblock>

  public :: t_optcpathfollowingdata

!<typeblock>

  ! This block saves parameters of the optimal control problem
  type t_settings_optcontrol
  
    ! $\alpha$ parameter of the optimal control functional
    real(DP) :: dalphaDistC = 1.0_DP
    
    ! $\alpha$ parameter of the optimal control functional via L2 boundary control
    real(DP) :: dalphaL2BdC = 1.0_DP

    ! $\alpha$ parameter of the optimal control functional via $H^{1/2}$ boundary control
    real(DP) :: dalphaH12BdC = 1.0_DP
    
    ! Penalty parameter for the dirichlet boundary control
    real(DP) :: ddirichletBCPenalty = 100.0_DP
    
    ! Penalty parameter for the end time condition in the dual equation.
    real(DP) :: dendTimeCondDualPenalty = 0.0_DP

    ! $\gamma$ parameter of the nonstationary optimal control functional
    real(DP) :: dalphaEndTimeC = 0.0_DP

    ! $\delta$ parameter of the terminal condition
    real(DP) :: ddeltaC = 0.0_DP
    
    ! $\sigma$-parameter that controls the blending in
    ! $F(u) = u - P( u-\sigma(\alpha u + \lambda) )$. Defaults to 1.0.
    ! A value of sigma=1/alpha would give
    ! $F(u) = u - P( -1/alpha \lambda) )$ which seems to be more instable.
    real(DP) :: dsigmaC = 1.0_DP
  
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
  !
  ! NOTE: Analytical constraints are not read from disc!
  ! The caller has to do that!
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

    ! Sigma-parameter that controls the formulation.
    call parlst_getvalue_double (rparlist,ssectionOptC,&
        "dsigmaC",roptcontrol%dsigmaC,roptcontrol%dsigmaC)
        
    ! Alpha/Gamma parameters.
    !
    ! These parameters may have two types: either it is a single
    ! parameter or a list of parameters. In the second case,
    ! this defines a path following strategy, i.e., the parameters
    ! change according to the list in every nonlinear iteration.
    call parlst_getvalue_string (rparlist,ssectionOptC,&
        "dalphaDistC",sstring,"-1.0",bdequote=.true.)

    ! Is this a list?
    call sys_countParameters (sstring,ntokens)
    
    if (ntokens .eq. 1) then
    
      read(sstring,*) roptcontrol%dalphaDistC
      
    else if (ntokens .gt. 1) then
      ! Read the list
      roptcontrol%dalphaDistC = 0.0_DP
      allocate(roptcontrol%rpathfollowingdata%p_dalphaDistC(ntokens))
      
      ! Loop through all tokens and save them
      istart = 1
      do itoken = 1,ntokens
        call sys_getParameter (sstring,stoken,istart)
        read(stoken,*) roptcontrol%rpathfollowingdata%p_dalphaDistC(itoken)
      end do
    end if
        
    call parlst_getvalue_string (rparlist,ssectionOptC,&
        "dalphaL2BdC",sstring,"-1.0",bdequote=.true.)

    ! Is this a list?
    call sys_countParameters (sstring,ntokens)

    if (ntokens .eq. 1) then
    
      read(sstring,*) roptcontrol%dalphaL2BdC
      
    else if (ntokens .gt. 1) then
      ! Read the list
      roptcontrol%dalphaL2BdC = 0.0_DP
      allocate(roptcontrol%rpathfollowingdata%p_dalphaL2BdC(ntokens))
      
      ! Loop through all tokens and save them
      istart = 1
      do itoken = 1,ntokens
        call sys_getParameter (sstring,stoken,istart)
        read(stoken,*) roptcontrol%rpathfollowingdata%p_dalphaL2BdC(itoken)
      end do
    end if

    call parlst_getvalue_string (rparlist,ssectionOptC,&
        "dalphaH12BdC",sstring,"-1.0",bdequote=.true.)

    ! Is this a list?
    call sys_countParameters (sstring,ntokens)

    if (ntokens .eq. 1) then
    
      read(sstring,*) roptcontrol%dalphaH12BdC
      
    else if (ntokens .gt. 1) then
      ! Read the list
      roptcontrol%dalphaH12BdC = 0.0_DP
      allocate(roptcontrol%rpathfollowingdata%p_dalphaH12BdC(ntokens))
      
      ! Loop through all tokens and save them
      istart = 1
      do itoken = 1,ntokens
        call sys_getParameter (sstring,stoken,istart)
        read(stoken,*) roptcontrol%rpathfollowingdata%p_dalphaH12BdC(itoken)
      end do
    end if

    call parlst_getvalue_string (rparlist,ssectionOptC,&
        "dalphaEndTimeC",sstring,"0.0",bdequote=.true.)

    ! Is this a list?
    call sys_countParameters (sstring,ntokens)

    if (ntokens .eq. 1) then
    
      read(sstring,*) roptcontrol%dalphaEndTimeC

    else if (ntokens .gt. 1) then
      ! Read the list
      roptcontrol%dalphaEndTimeC = 0.0_DP
      allocate(roptcontrol%rpathfollowingdata%p_dalphaEndTimeC(ntokens))
      
      ! Loop through all tokens and save them
      istart = 1
      do itoken = 1,ntokens
        call sys_getParameter (sstring,stoken,istart)
        read(stoken,*) roptcontrol%rpathfollowingdata%p_dalphaEndTimeC(itoken)
      end do
    end if

    ! Penalty parameter for the dirichlet boundary control            
    call parlst_getvalue_double (rparlist,ssectionOptC,&
        "ddirichletBCPenalty",roptcontrol%ddirichletBCPenalty,100.0_DP)
     
    # Penalty parameter for the end time condition in the dual equation.
    call parlst_getvalue_double (rparlist,ssectionOptC,&
        "dendTimeCondDualPenalty",roptcontrol%dendTimeCondDualPenalty,0.0_DP)
        
    ! Type of the formulation
    call parlst_getvalue_int (rparlist,ssectionOptC,&
        "ispaceTimeFormulation",roptcontrol%ispaceTimeFormulation,0)
    
    call parlst_getvalue_int (rparlist,ssectionOptC,&
        "iconvectionExplicit",roptcontrol%iconvectionExplicit,0)

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

    ! -------------------------------------------
    ! Read the constraints
    ! -------------------------------------------

    ! Distributed control
    call parlst_getvalue_string (rparlist,ssectionOptC,&
        "ssectionConstrDistCtrl",sstring,"",bdequote=.true.)
        
    call soptc_getConstraints (rparlist,sstring,&
        roptcontrol%rconstraints%rconstraintsDistCtrl)

    ! L2 boundary control
    call parlst_getvalue_string (rparlist,ssectionOptC,&
        "ssectionConstrL2BdC",sstring,"",bdequote=.true.)
        
    call soptc_getConstraints (rparlist,sstring,&
        roptcontrol%rconstraints%rconstraintsL2BdC)

    ! H^^1/2 boundary control
    call parlst_getvalue_string (rparlist,ssectionOptC,&
        "ssectionConstrH12BdC",sstring,"",bdequote=.true.)
        
    call soptc_getConstraints (rparlist,sstring,&
        roptcontrol%rconstraints%rconstraintsH12BdC)

    ! State
    call parlst_getvalue_string (rparlist,ssectionOptC,&
        "ssectionConstrState",sstring,"",bdequote=.true.)
        
    call soptc_getConstraints (rparlist,sstring,&
        roptcontrol%rconstraints%rconstraintsState)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine soptc_getConstraints (rparlist,ssection,rconstraints)
  
!<description>
  ! Reads the parameters that defines a set of constraints.
  !
  ! NOTE: Analytical constraints are not read from disc!
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section where the parameters of the constraints can be found.
  character(len=*), intent(in) :: ssection
!</input>

!<inputoutput>
  ! Constraints-structure, to be set up.
  type(t_optcConstraintsDef), intent(out) :: rconstraints
!</inputoutput>

!</subroutine>
    
    call parlst_getvalue_int (rparlist,ssection,&
        "cconstraints",rconstraints%cconstraints,0)

    call parlst_getvalue_int (rparlist,ssection,&
        "cconstrType",rconstraints%cconstrType,0)

    call parlst_getvalue_double (rparlist,ssection,&
        "dmin1",rconstraints%dmin1,-1.0E10_DP)

    call parlst_getvalue_double (rparlist,ssection,&
        "dmax1",rconstraints%dmax1,1.0E10_DP)

    call parlst_getvalue_double (rparlist,ssection,&
        "dmin2",rconstraints%dmin2,-1.0E10_DP)

    call parlst_getvalue_double (rparlist,ssection,&
        "dmax2",rconstraints%dmax2,1.0E10_DP)

    call parlst_getvalue_double (rparlist,ssection,&
        "dregularisation",rconstraints%dregularisation,1.0_DP)
        
    call parlst_getvalue_string (rparlist,ssection,&
        "ssectionMin1",rconstraints%ssectionMin1,"",bdequote=.true.)

    call parlst_getvalue_string (rparlist,ssection,&
        "ssectionMax1",rconstraints%ssectionMax1,"",bdequote=.true.)

    call parlst_getvalue_string (rparlist,ssection,&
        "ssectionMin2",rconstraints%ssectionMin2,"",bdequote=.true.)

    call parlst_getvalue_string (rparlist,ssection,&
        "ssectionMax2",rconstraints%ssectionMax2,"",bdequote=.true.)
        
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
    if (associated(roptcontrol%rpathfollowingdata%p_dalphaDistC)) then
      deallocate(roptcontrol%rpathfollowingdata%p_dalphaDistC)
    end if

    if (associated(roptcontrol%rpathfollowingdata%p_dalphaH12BdC)) then
      deallocate(roptcontrol%rpathfollowingdata%p_dalphaH12BdC)
    end if

    if (associated(roptcontrol%rpathfollowingdata%p_dalphaL2BdC)) then
      deallocate(roptcontrol%rpathfollowingdata%p_dalphaL2BdC)
    end if

    if (associated(roptcontrol%rpathfollowingdata%p_dalphaEndTimeC)) then
      deallocate(roptcontrol%rpathfollowingdata%p_dalphaEndTimeC)
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
    if (associated(roptcontrol%rpathfollowingdata%p_dalphaDistC)) then
      roptControl%dalphaDistC = roptControl%rpathfollowingdata%p_dalphaDistC(&
              min(size(roptcontrol%rpathfollowingdata%p_dalphaDistC),istep))
      if (bprint) then
        call output_line ("Path-Follower: Parameter ALPHA_DIST = "//&
            trim(sys_sdEL(roptcontrol%dalphaDistC,10)))
      end if
    end if

    if (associated(roptcontrol%rpathfollowingdata%p_dalphaH12BdC)) then
      roptcontrol%dalphaL2BdC = roptControl%rpathfollowingdata%p_dalphaH12BdC(&
              min(size(roptcontrol%rpathfollowingdata%p_dalphaH12BdC),istep))
      if (bprint) then
        call output_line ("Path-Follower: Parameter ALPHA_L2BDC  = "//&
            trim(sys_sdEL(roptcontrol%dalphaL2BdC,10)))
      end if
    end if

    if (associated(roptcontrol%rpathfollowingdata%p_dalphaH12BdC)) then
      roptcontrol%dalphaH12BdC = roptControl%rpathfollowingdata%p_dalphaH12BdC(&
              min(size(roptcontrol%rpathfollowingdata%p_dalphaH12BdC),istep))
      if (bprint) then
        call output_line ("Path-Follower: Parameter ALPHA_H12BDC  = "//&
            trim(sys_sdEL(roptcontrol%dalphaH12BdC,10)))
      end if
    end if

    if (associated(roptcontrol%rpathfollowingdata%p_dalphaEndTimeC)) then
      roptcontrol%dalphaEndTimeC = roptControl%rpathfollowingdata%p_dalphaEndTimeC(&
              min(size(roptcontrol%rpathfollowingdata%p_dalphaEndTimeC),istep))
      if (bprint) then
        call output_line ("Path-Follower: Parameter GAMMA = "//&
            trim(sys_sdEL(roptcontrol%dalphaEndTimeC,10)))
      end if
    end if

  end subroutine

end module
