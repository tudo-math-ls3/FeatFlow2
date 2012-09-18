!##############################################################################
!# ****************************************************************************
!# <name> structuresdiscretisation </name>
!# ****************************************************************************
!#
!# <purpose>
!# Discretisation related structures.
!#
!# 1.) struc_getSpaceDiscrSettings
!#     -> Initialise the parameters defining the discretisation
!#        (element, cubature, stabilisation)
!#
!# 2.) struc_initPhysics
!#     -> Initialise parameters concerning the physics of the problem
!#        (type of equation,...)
!# </purpose>
!##############################################################################

module structuresdiscretisation

  use fsystem
  use cubature
  use paramlist
  use spatialdiscretisation
  
  use constantsdiscretisation
  
  implicit none
  
  private
  
!<types>

!<typeblock>

  ! This type block encapsules the settings for the stabilisation of the convection.
  type t_settings_stabil
  
    ! Type of stabilization of the convective term.
    ! 3=Streamline Diffusion (new implementation)
    integer :: cupwind = 3
    
    ! Relaxation parameter for upwind/Streamline Diffusion/Jump stabilisation.
    ! Standard values: Streamline diffusion=1.0, Upwind=0.1, Jump stabil=0.01.
    ! dUpsam1 holds for the primal equation, dUpsam2 for the dual one.
    real(DP) :: dupsam = 0.0_DP
    
    ! Defines whether or not the convection operator of the dual equation
    ! is included into the solution.
    ! =0: no, =1: yes
    integer :: cconvectionOnBoundaryDefect = 1

    ! Defines whether or not the convection operator of the dual equation
    ! is set up for the preconditioners on the boundary.
    ! =0: no, =1: yes
    integer :: cconvectionOnBoundaryMatrix = 1

  end type

!</typeblock>

  public :: t_settings_stabil

!<typeblock>
  
  ! Structure encapsuling the main discretisation
  type t_settings_spacediscr
  
    ! Type of element pair to use for the discretisation.
    ! 2D Stokes/Navier-Stokes:
    !  0 = Q1~(E031) / Q1~(E031) / Q0
    !  1 = Q1~(E030) / Q1~(E030) / Q0
    !  2 = Q1~(EM31) / Q1~(EM31) / Q0
    !  3 = Q1~(EM30) / Q1~(EM30) / Q0 = standard
    !  4 = Q2 (E013) / Q2 (E013) / QP1
    !  5 = Q1~(EM30) / Q1~(EM30) / Q0 unpivoted (much faster than 3 but less stable)
    !  6 = Q1~(EM30) / Q1~(EM30) / Q0 unscaled (slightly faster than 3 but less stable)
    ! (EM30 = nonparametric, nonconformal Rannacher-Turek element)
    ! (QP1  = Quadrilateral discontinuous P1 element)
    integer :: ielementType = 3
    
    ! cubature formula for Mass matrix
    integer(i32) :: icubMass = CUB_GEN_AUTO

    ! cubature formula for Stokes/Laplacian matrix
    integer(i32) :: icubStokes = CUB_GEN_AUTO

    ! cubature formula for Pressure matrices B
    integer(i32) :: icubB = CUB_GEN_AUTO

    ! cubature formula for RHS F
    integer(i32) :: icubF = CUB_GEN_AUTO

    ! Support for integral mean value constraints.
    ! If set =1, integral mean value constraints are supported. The matrix
    ! structure is enlarged appropriately. This slows down the computation.
    integer :: csupportIntMeanConstr = 0

    ! Stabilisation parameters for the convection in the primal system.
    type(t_settings_stabil) :: rstabilConvecPrimal

    ! Stabilisation parameters for the convection in the dual system.
    type(t_settings_stabil) :: rstabilConvecDual

  end type

!</typeblock>

  public :: t_settings_spacediscr

!<typeblock>

  ! This type block encapsules all physical constants and configuration
  ! parameters for the primal equation. This includes e.g. the type of the equation,
  ! viscosity parameter etc.
  type t_settings_physics
  
    ! Type of problem.
    ! =0: Navier-Stokes 2D.
    ! =1: Stokes 2D.
    integer :: cequation = CCEQ_NAVIERSTOKES2D
    
    ! Type of subproblem of the main problem. Depending on iequationType.
    ! If iequationType=0 or =1:
    ! =0: (Navier-)Stokes with gradient tensor
    ! =1: (Navier-)Stokes with deformation tensor
    integer :: isubEquation = 0
    
    ! Model for the viscosity.
    ! =0: Constant viscosity.
    ! =1: Power law: nu = nu_0 * z^(dviscoexponent/2 - 1), nu_0 = 1/RE, z=||D(u)||^2+dviscoEps
    ! =2: Bingham fluid: nu = nu_0 + dviscoyield / sqrt(|D(u)||^2+dviscoEps^2), nu_0 = 1/RE
    integer :: cviscoModel = 0
        
    ! Exponent parameter for the viscosity model
    real(DP) :: dviscoexponent = 2.0_DP

    ! Epsilon regularisation for the viscosity model
    real(DP) :: dviscoEps = 0.01_DP
    
    ! Yield stress for Bingham fluid
    real(DP) :: dviscoYield = 1.0_DP

    ! Viscosity parameter nu = 1/Re if cviscoModel=0 (constant viscosity).
    ! Otherwise not used.
    real(DP) :: dnuConst = 0.001_DP
      
  end type

!</typeblock>

  public :: t_settings_physics

!</types>

  ! Extracts main discretisation settings in space from the DAT file.
  public :: struc_getSpaceDiscrSettings
  
  ! Reads information about physics from the DAT file and saves them to rphysics.
  public :: struc_initPhysics

contains

  ! ***************************************************************************

!<subroutine>

  subroutine struc_initStabil (rparlist,rstabilPrimal,rstabilDual,ssection)
  
!<description>
  ! Reads information about stabilisation from the DAT file and saves them
  ! to rstabilPrimal/rstabilDual.
!</description>
  
!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section where the parameters can be found.
  character(len=*), intent(in) :: ssection
!</input>
  
!<output>
  ! Stabilisation parameters for the primal equation
  type(t_settings_stabil), intent(out) :: rstabilPrimal

  ! Stabilisation parameters for the dual equation
  type(t_settings_stabil), intent(out) :: rstabilDual
!</output>

!</subroutine>

    ! Get the parameters
    call parlst_getvalue_int (rparlist, ssection, &
        "IUPWIND1", rstabilPrimal%cupwind)
    call parlst_getvalue_int (rparlist, ssection, &
        "IUPWIND2", rstabilDual%cupwind)
    call parlst_getvalue_double (rparlist, ssection, &
        "DUPSAM1", rstabilPrimal%dupsam)
    call parlst_getvalue_double (rparlist, ssection, &
        "DUPSAM2", rstabilDual%dupsam)

    ! Get the flags that decides if in the dual equation, the convection is
    ! set up on the boundary. For the primal equation, this is always set up,
    ! so pass the values only to the stabilisation structure of the dual equation.
    call parlst_getvalue_int (rparlist, ssection, &
        "CCONVECTIONONBOUNDARYDEFECT", rstabilDual%cconvectionOnBoundaryDefect,1)

    call parlst_getvalue_int (rparlist, ssection, &
        "CCONVECTIONONBOUNDARYMATRIX", rstabilDual%cconvectionOnBoundaryMatrix,1)

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>

  subroutine struc_getSpaceDiscrSettings (rparlist,ssection,rsettingsSpaceDiscr)

!<description>
  ! Extracts main discretisation settings in space from the DAT file.
!</description>

!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section containing the parameters
  character(len=*), intent(in) :: ssection
!</input>

!<output>
   ! Structure receiving the main discretisation settings
   type(t_settings_spacediscr), intent(out) :: rsettingsSpaceDiscr
!</output>

!</subroutine>

    character(len=SYS_STRLEN) :: sstr
    integer :: icub

    call parlst_getvalue_int (rparlist,ssection,&
        "ielementType",rsettingsSpaceDiscr%ielementType,3)
                              
    call parlst_getvalue_string (rparlist,ssection,"scubStokes",sstr,"")
    if (sstr .eq. "") then
      call parlst_getvalue_int (rparlist,ssection,&
          "icubStokes",icub,int(SPDISC_CUB_AUTOMATIC))
    else
      icub = cub_igetID(sstr)
    end if
    rsettingsSpaceDiscr%icubStokes = icub

    call parlst_getvalue_string (rparlist,ssection,"scubB",sstr,"")
    if (sstr .eq. "") then
      call parlst_getvalue_int (rparlist,ssection,&
          "icubB",icub,int(SPDISC_CUB_AUTOMATIC))
    else
      icub = cub_igetID(sstr)
    end if
    rsettingsSpaceDiscr%icubB = icub

    call parlst_getvalue_string (rparlist,ssection,"scubF",sstr,"")
    if (sstr .eq. "") then
      call parlst_getvalue_int (rparlist,ssection,&
          "icubF",icub,int(SPDISC_CUB_AUTOMATIC))
    else
      icub = cub_igetID(sstr)
    end if
    rsettingsSpaceDiscr%icubF = icub
    
    ! Which cubature rule to use for mass matrices?
    call parlst_getvalue_string (rparlist,ssection,"scubMass",sstr,"")
    if (sstr .eq. "") then
      call parlst_getvalue_int (rparlist,ssection,&
          "icubM",icub,int(SPDISC_CUB_AUTOMATIC))
    else
      icub = cub_igetID(sstr)
    end if
    rsettingsSpaceDiscr%icubMass = icub

    ! Initialise stabilisation parameters
    call struc_initStabil (rparlist,&
        rsettingsSpaceDiscr%rstabilConvecPrimal,&
        rsettingsSpaceDiscr%rstabilConvecDual,ssection)

    ! Support for integral mean value constraints.
    call parlst_getvalue_int (rparlist,ssection,&
        "csupportIntMeanConstr",rsettingsSpaceDiscr%csupportIntMeanConstr,0)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine struc_initPhysics (rparlist,rphysics,ssection)
  
!<description>
  ! Reads information about physics from the DAT file and saves them to rphysics.
!</description>
  
!<input>
  ! Parameter list
  type(t_parlist), intent(in) :: rparlist
  
  ! Section where the parameters can be found.
  character(len=*), intent(in) :: ssection
!</input>
  
!<output>
  ! A structure receiving the physics of the equation
  type(t_settings_physics), intent(out) :: rphysics
!</output>

!</subroutine>

    ! Which type of problem to discretise? (Stokes, Navier-Stokes,...)
    call parlst_getvalue_int (rparlist,ssection,&
        "cequation",rphysics%cequation,CCEQ_NAVIERSTOKES2D)

    ! Type of subproblem (gradient tensor, deformation tensor,...)
    call parlst_getvalue_int (rparlist,ssection,&
        "csubEquation",rphysics%isubEquation,0)

    ! Get the viscosity model
    ! Standard = 0 = constant viscosity
    call parlst_getvalue_int (rparlist,ssection,&
        "cviscoModel",rphysics%cviscoModel,0)
                                 
    call parlst_getvalue_double (rparlist,ssection,&
        "dviscoexponent",rphysics%dviscoexponent,2.0_DP)
                                 
    call parlst_getvalue_double (rparlist,ssection,&
        "dviscoEps",rphysics%dviscoEps,0.01_DP)

    call parlst_getvalue_double (rparlist,ssection,&
        "dviscoYield",rphysics%dviscoYield,1.0_DP)

    ! Get the viscosity parameter.
    ! Note that the parameter in the DAT file is 1/nu !
    call parlst_getvalue_double (rparlist,ssection,&
        "RE",rphysics%dnuConst,1000.0_DP)
    rphysics%dnuConst = 1E0_DP/rphysics%dnuConst

  end subroutine
  
end module
