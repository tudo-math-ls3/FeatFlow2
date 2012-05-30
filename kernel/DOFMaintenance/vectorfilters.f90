!##############################################################################
!# ****************************************************************************
!# <name> vectorfilters </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains a set of filters that can be applied to a scalar
!# or block vector. These can be used e.g. during a solution process to impose
!# different quantities into solution and/or defect vectors.
!#
!# The discrete boundary conditions realised in the module 'bcassembly' are
!# one type of filter. While being created in the module 'bcassembly',
!# this module now provides the functionality to impose discrete boundary
!# conditions into a vector.
!# Also other filters can be found here, e.g. normalisation ov vectors, etc.
!#
!# 'Linear' and 'nonlinear' vector filters \\
!# --------------------------------------- \\
!# There are two elemental types of vector filters realised here:
!# 'Linear' filters and 'nonlinear' filters:
!#
!# a) 'Linear' filters
!#
!#   These filters realise simple vector filters like
!#   - implementation of Dirichlet boundary conditions
!#   - filter a vector to be in <tex>$L^2_0$</tex>
!#   or similar.
!#   The name comes from the ability to be able to be used as a filter
!#   when solving a *linear system*. Some iterative solvers like BiCGStab
!#   or Multigrid allow a vector to be filtered during the iteration.
!#   All filters marked as 'linear filters' can be collected to a
!#   'filter chain' that is applied to a vector in such an iteration.
!#   For this purpose, there exists a higher-level module 'filtersupport',
!#   which realises such a filter chain.
!#   Of course, they can be applied to a vector anytime manually, e.g. for
!#   implementing Dirichlet boundary conditions 'by hand'.
!#
!# b) 'Nonlinear' filters
!#
!#   All filters that can not be collected in a filter chain to be
!#   applied during the solution of a linear system are called
!#   'nonlinear filters'. These filters are usually called inside of
!#   a nonlinear iteration to perform a special filtering to a vector
!#   which is probably not possible to formulate in the calling convention
!#   of a linear filter - e.g. if additional auxiliary vectors are used
!#   or if a filter consists of a predictor-corrector step or similar.
!#   Example:
!#   - implementation of Slip boundary conditions.
!#   - implementation of Pressure Drop boundary conditions.
!#
!# Note that filters are allowed consist a linear and a nonlinear part!
!# In such a case, the 'nonlinear' filter part is usually called during
!# the nonlinear iteration, while the 'linear' part is used during the
!# solution of a linear system. An example for this may be the
!# predictor-corrector implementation of Slip boundary conditions (not
!# realised here), where the nonlinear iteration treats the BC while
!# in the solution process of the linear system, all respective nodes
!# are handled as Dirichlet.
!#
!# The following routines can be found here:
!#
!#  1.) vecfil_normaliseToL20Sca
!#      -> Linear filter
!#      -> Normalise a scalar vector to be in the space <tex>$L^2_0$</tex>.
!#
!#  2.) vecfil_discreteBCsol
!#      -> Linear filter
!#      -> Apply the 'discrete boundary conditions for solution vectors' filter
!#         onto a given (block) solution vector.
!#
!#  3.) vecfil_discreteBCrhs
!#      -> Linear filter
!#      -> Apply the 'discrete boundary conditions for RHS vectors' filter
!#         onto a given (block) vector.
!#
!#  4.) vecfil_discreteBCdef
!#      -> Linear filter
!#      -> Apply the 'discrete boundary conditions for defect vectors' filter
!#         onto a given (block) vector.
!#
!#  5.) vecfil_discreteFBCsol
!#      -> Linear filter
!#      -> Apply the 'discrete fictitious boundary conditions for solution vectors'
!#         filter onto a given (block) solution vector.
!#
!#  6.) vecfil_discreteFBCrhs
!#      -> Linear filter
!#      -> Apply the 'discrete fictitious boundary conditions for RHS vectors'
!#         filter onto a given (block) vector.
!#
!#  7.) vecfil_discreteFBCdef
!#      -> Linear filter
!#      -> Apply the 'discrete fictitious boundary conditions for defect vectors'
!#         filter onto a given (block) vector.
!#
!#  8.) vecfil_subvectorToL20
!#      -> Linear filter
!#      -> Normalise a subvector of a block vector to be in the space <tex>$L^2_0$</tex>.
!#
!#  9.) vecfil_discreteNLSlipBCdef
!#      -> Nonlinear filter
!#      -> Implements discrete nonlinear slip BC`s into a scalar defect vector.
!#
!# 10.) vecfil_discreteNLPDropBCrhs
!#      -> Nonlinear filter
!#      -> Implements discrete pressure drop BC`s into a RHS vector.
!#
!# 11.) vecfil_subvectorSmallL1To0
!#      -> Linear filter
!#      -> Normalise a subvector to have an l1 vector sum = 0.
!#
!# Auxiliary routines, usually not called by the main program:
!#
!#  1.) vecfil_imposeDirichletBC
!#      -> Implements discrete Dirichlet BC`s into a scalar vector.
!#
!#  2.) vecfil_imposeDirichletDefectBC
!#      -> Implements discrete Dirichlet BC`s into a scalar defect vector.
!#
!#  3.)  vecfil_imposeDirichletFBC (rx,icomponent,rdbcStructure)
!#      -> Implements discrete Dirichlet fictitious boundary conditions into a
!#      -> scalar vector.
!#
!#  4.) vecfil_normaliseToL20Sca (rx)
!#      -> Normalises a scalar vector to bring it into the space <tex>$L^2_0$</tex>.
!#
!#  5.) vecfil_imposePressureDropBC (rx,dtimeweight,rpdbcStructure)
!#      -> Implements discrete pressure drop BC`s into a block vector.
!#
!#  6.) vecfil_imposeNLSlipDefectBC
!#      -> Implements discrete nonlinear slip BC`s into a scalar defect vector
!#         as configured in the slip BC structure.
!#
!#  7.) vecfil_normaliseSmallL1To0Sca (rx)
!#      -> Normalise a scalar vector to have an l1 vector sum = 0.
!#
!# </purpose>
!##############################################################################

module vectorfilters

!$use omp_lib
  use fsystem
  use genoutput
  use storage
  use basicgeometry
  use linearsystemscalar
  use linearsystemblock
  use element
  use elementbase
  use elementpreprocessing
  use element_quad2d
  use discretebc
  use discretefbc
  use dofmapping
  use spatialdiscretisation
  use triangulation
  use derivatives
  use transformation

  implicit none

  private

  public :: vecfil_normaliseToL20Sca
  public :: vecfil_discreteBCsol
  public :: vecfil_discreteBCrhs
  public :: vecfil_discreteBCdef
  public :: vecfil_discreteFBCsol
  public :: vecfil_discreteFBCrhs
  public :: vecfil_discreteFBCdef
  public :: vecfil_subvectorToL20
  public :: vecfil_discreteNLSlipBCdef
  public :: vecfil_discreteNLPDropBCrhs
  public :: vecfil_subvectorSmallL1To0
  public :: vecfil_imposeDirichletBC
  public :: vecfil_imposeDirichletDefectBC
  public :: vecfil_imposeDirichletFBC
  public :: vecfil_imposePressureDropBC
  public :: vecfil_imposeNLSlipDefectBC
  public :: vecfil_normaliseSmallL1To0Sca

contains

! *****************************************************************************
! Scalar vector filters
! *****************************************************************************

!<subroutine>

  subroutine vecfil_imposeDirichletBC (rx,rdbcStructure)

!<description>
  ! Implements discrete Dirichlet BC`s into a scalar vector.
!</description>

!<input>
  ! The t_discreteBCDirichlet that describes the discrete Dirichlet BC`s
  type(t_discreteBCDirichlet), intent(in), target  :: rdbcStructure
!</input>

!<inputoutput>

  ! The scalar vector where the boundary conditions should be imposed.
  type(t_vectorScalar), intent(inout), target :: rx

!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i
    real(DP), dimension(:), pointer    :: p_vec
    integer, dimension(:), pointer :: p_idx
    real(DP), dimension(:), pointer    :: p_val
    integer, dimension(:), pointer :: p_Iperm

    ! If nDOF=0, there are no DOF`s the current boundary condition segment,
    ! so we do not have to do anything. Maybe the case if the user selected
    ! a boundary region that is just too small.
    if (rdbcStructure%nDOF .eq. 0) return

    ! Get pointers to the structures. For the vector, get the pointer from
    ! the storage management.

    if (rdbcStructure%h_DdirichletValues .eq. ST_NOHANDLE) then
      call output_line('Dirichlet BC''s not correctly configured!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeDirichletBC')
      call output_line('Are the BC''s only configured for defect values?!?',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeDirichletBC')
      call sys_halt()
    end if

    if (rdbcStructure%h_IdirichletDOFs .eq. ST_NOHANDLE) then
      call output_line('DBC not configured',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeDirichletBC')
      call sys_halt()
    end if

    call storage_getbase_int(rdbcStructure%h_IdirichletDOFs,p_idx)
    call storage_getbase_double(rdbcStructure%h_DdirichletValues,p_val)

    ! Impose the DOF value directly into the vector - more precisely, into the
    ! components of the subvector that is indexed by icomponent.

    if ((.not.associated(p_idx)).or.(.not.associated(p_val))) then
      call output_line('DBC not configured',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeDirichletBC')
      call sys_halt()
    end if

    call lsyssc_getbase_double (rx, p_vec)

    if (.not.associated(p_vec)) then
      call output_line('No vector!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeDirichletBC')
      call sys_halt()
    end if

    ! Only handle nDOF DOF`s, not the complete array!
    ! Probably, the array is longer (e.g. has the length of the vector), but
    ! contains only some entries...

    ! Is the vector sorted?
    if (rx%isortStrategy .le. 0) then
      ! No. Implement directly.
      do i=1,rdbcStructure%nDOF
        p_vec(p_idx(i)) = p_val(i)
      end do
    else
      ! Ups, vector sorted. At first get the permutation how its sorted
      ! - or more precisely, the back-permutation, as we need this one for
      ! the loop below.
      call storage_getbase_int (rx%h_IsortPermutation,p_Iperm)
      p_Iperm => p_Iperm(rx%NEQ+1:)

      ! And 'filter' each DOF during the boundary value implementation!
      do i=1,rdbcStructure%nDOF
        p_vec(p_Iperm(p_idx(i))) = p_val(i)
      end do
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vecfil_imposeDirichletDefectBC (rx,rdbcStructure)

!<description>
  ! Implements discrete Dirichlet BC`s into a scalar defect vector.
!</description>

!<input>
  ! The t_discreteBCDirichlet that describes the discrete Dirichlet BC`s
  type(t_discreteBCDirichlet), intent(in),target  :: rdbcStructure
!</input>

!<inputoutput>
  ! The scalar vector where the boundary conditions should be imposed.
  type(t_vectorScalar), intent(inout), target :: rx
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i
    real(DP), dimension(:), pointer :: p_vec
    integer, dimension(:), pointer :: p_idx
    integer, dimension(:), pointer :: p_Iperm

    ! If nDOF=0, there are no DOF`s the current boundary condition segment,
    ! so we do not have to do anything. Maybe the case if the user selected
    ! a boundary region that is just too small.
    if (rdbcStructure%nDOF .eq. 0) return

    ! Get pointers to the structures. For the vector, get the pointer from
    ! the storage management.

    if (rdbcStructure%h_IdirichletDOFs .eq. ST_NOHANDLE) then
      call output_line('DBC not configured',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeDirichletDefectBC')
      call sys_halt()
    end if

    call storage_getbase_int(rdbcStructure%h_IdirichletDOFs,p_idx)

    if (.not.associated(p_idx)) then
      call output_line('DBC not configured',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeDirichletDefectBC')
      call sys_halt()
    end if

    call lsyssc_getbase_double (rx, p_vec)

    if (.not.associated(p_vec)) then
      call output_line('No vector!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeDirichletDefectBC')
      call sys_halt()
    end if

    ! Impose the BC-DOF`s directly - more precisely, into the
    ! components of the subvector that is indexed by icomponent.
    !
    ! Only handle nDOF DOF`s, not the complete array!
    ! Probably, the array is longer (e.g. has the length of the vector), but
    ! contains only some entries...

    ! Is the vector sorted?
    if (rx%isortStrategy .le. 0) then
      ! No. Implement directly.
      do i=1,rdbcStructure%nDOF
        p_vec(p_idx(i)) = 0.0_DP
      end do
    else
      ! Ups, vector sorted. At first get the permutation how its sorted -
      ! or more precisely, the back-permutation, as we need this one for
      ! the loop below.
      call storage_getbase_int (rx%h_IsortPermutation,p_Iperm)
      p_Iperm => p_Iperm(rx%NEQ+1:)

      ! And 'filter' each DOF during the boundary value implementation!
      do i=1,rdbcStructure%nDOF
        p_vec(p_Iperm(p_idx(i))) = 0.0_DP
      end do
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vecfil_imposeDirichletFBC (rx,icomponent,rdbcStructure)

!<description>
  ! Implements discrete Dirichlet fictitious boundary conditions into a
  ! scalar vector.
!</description>

!<input>
  ! The t_discreteFBCDirichlet that describes the discrete Dirichlet BC`s
  type(t_discreteFBCDirichlet), intent(in), target  :: rdbcStructure

  ! Index of the solution component in rdbcStructure\%Icomponent
  integer, intent(in) :: icomponent
!</input>

!<inputoutput>

  ! The scalar vector where the boundary conditions should be imposed.
  type(t_vectorScalar), intent(inout), target :: rx

!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i
    real(DP), dimension(:), pointer    :: p_vec
    integer, dimension(:), pointer :: p_idx
    real(DP), dimension(:,:), pointer    :: p_val
    integer, dimension(:), pointer :: p_Iperm

    ! If nDOF=0, there are no DOF`s the current boundary condition segment,
    ! so we do not have to do anything. Maybe the case if the user selected
    ! a boundary region that is just too small.
    if (rdbcStructure%nDOF .eq. 0) return

    ! Get pointers to the structures. For the vector, get the pointer from
    ! the storage management.

    if ((rdbcStructure%h_IdirichletDOFs .eq. ST_NOHANDLE) .or. &
        (rdbcStructure%h_DdirichletValues .eq. ST_NOHANDLE)) then
      call output_line('Dirichlet BC''s not correctly configured!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeDirichletFBC')
      call output_line('Are the BC''s only configured for defect values?!?',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeDirichletFBC')
      call sys_halt()
    end if

    call storage_getbase_int(rdbcStructure%h_IdirichletDOFs,p_idx)
    call storage_getbase_double2d(rdbcStructure%h_DdirichletValues,p_val)

    if ((.not.associated(p_idx)).or.(.not.associated(p_val))) then
      call output_line('DBC not configured',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeDirichletFBC')
      call sys_halt()
    end if

    ! Impose the DOF value directly into the vector - more precisely, into the
    ! components of the subvector that is indexed by icomponent.

    call lsyssc_getbase_double (rx, p_vec)

    if (.not.associated(p_vec)) then
      call output_line('Error: No vector',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeDirichletFBC')
      call sys_halt()
    end if

    ! Only handle nDOF DOF`s, not the complete array!
    ! Probably, the array is longer (e.g. has the length of the vector), but
    ! contains only some entries...

    ! Is the vector sorted?
    if (rx%isortStrategy .le. 0) then
      ! No. Implement directly.
      do i=1,rdbcStructure%nDOF
        p_vec(p_idx(i)) = p_val(icomponent,i)
      end do
    else
      ! Ups, vector sorted. At first get the permutation how its sorted
      ! - or more precisely, the back-permutation, as we need this one for
      ! the loop below.
      call storage_getbase_int (rx%h_IsortPermutation,p_Iperm)
      p_Iperm => p_Iperm(rx%NEQ+1:)

      ! And 'filter' each DOF during the boundary value implementation!
      do i=1,rdbcStructure%nDOF
        p_vec(p_Iperm(p_idx(i))) = p_val(icomponent,i)
      end do
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vecfil_imposeDirichletDefectFBC (rx,rdbcStructure)

!<description>
  ! Implements discrete Dirichlet fictitious boundary conditions into
  ! a scalar defect vector.
!</description>

!<input>
  ! The t_discreteFBCDirichlet that describes the discrete Dirichlet BC`s
  type(t_discreteFBCDirichlet), intent(in),target  :: rdbcStructure
!</input>

!<inputoutput>
  ! The scalar vector where the boundary conditions should be imposed.
  type(t_vectorScalar), intent(inout), target :: rx
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i
    real(DP), dimension(:), pointer :: p_vec
    integer, dimension(:), pointer :: p_idx
    integer, dimension(:), pointer :: p_Iperm

    ! If nDOF=0, there are no DOF`s the current boundary condition segment,
    ! so we do not have to do anything. Maybe the case if the user selected
    ! a boundary region that is just too small.
    if (rdbcStructure%nDOF .eq. 0) return

    ! Get pointers to the structures. For the vector, get the pointer from
    ! the storage management.

    if (rdbcStructure%h_IdirichletDOFs .eq. ST_NOHANDLE) then
      call output_line('DBC not configured',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeDirichletDefectFBC')
      call sys_halt()
    end if

    call storage_getbase_int(rdbcStructure%h_IdirichletDOFs,p_idx)

    if (.not.associated(p_idx)) then
      call output_line('DBC not configured',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeDirichletDefectFBC')
      call sys_halt()
    end if

    call lsyssc_getbase_double (rx, p_vec)

    if (.not.associated(p_vec)) then
      call output_line('No vector!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeDirichletDefectFBC')
      call sys_halt()
    end if

    ! Impose the BC-DOF`s directly - more precisely, into the
    ! components of the subvector that is indexed by icomponent.
    !
    ! Only handle nDOF DOF`s, not the complete array!
    ! Probably, the array is longer (e.g. has the length of the vector), but
    ! contains only some entries...

    ! Is the vector sorted?
    if (rx%isortStrategy .le. 0) then
      ! No. Implement directly.
      do i=1,rdbcStructure%nDOF
        p_vec(p_idx(i)) = 0.0_DP
      end do
    else
      ! Ups, vector sorted. At first get the permutation how its sorted -
      ! or more precisely, the back-permutation, as we need this one for
      ! the loop below.
      call storage_getbase_int (rx%h_IsortPermutation,p_Iperm)
      p_Iperm => p_Iperm(rx%NEQ+1:)

      ! And 'filter' each DOF during the boundary value implementation!
      do i=1,rdbcStructure%nDOF
        p_vec(p_Iperm(p_idx(i))) = 0.0_DP
      end do
    end if

  end subroutine

  ! ***************************************************************************
  ! Other scalar filters
  ! ***************************************************************************

!<subroutine>

  subroutine vecfil_normaliseToL20Sca (rx)

!<description>
  ! This routine normalises a scalar vector to bring it into the space <tex>$L^2_0$</tex>.
!</description>

!<inputoutput>
  ! The vector which is to be normalised.
  type(t_vectorScalar), intent(inout),target :: rx
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_spatialDiscretisation), pointer :: p_rdiscretisation
    real(DP), dimension(:), pointer        :: p_DelementArea,p_Ddata
    real(DP) :: dpintegral,c
    integer :: iel,nel,itrialspace

    ! Get the discretisation...
    if (.not. associated(rx%p_rspatialDiscr)) return

    p_rdiscretisation => rx%p_rspatialDiscr

    ! ... and check it. If we have a uniform discretisation with P_0/Q_0,
    ! we have the easy case, that the integral of the function rx is
    ! representing is area*rx. Otherwise we have to calculate the integral
    ! which is somehow more costly...


    if (p_rdiscretisation%ccomplexity .eq. SPDISC_UNIFORM) then

      itrialspace = elem_getPrimaryElement(&
          p_rdiscretisation%RelementDistr(1)%celement)

      select case (itrialspace)
      case (EL_P0, EL_Q0)

        ! Ok, easy case. Get from the triangulation the AREA-array for calculating
        ! a simple integral of rx:

        call storage_getbase_double (p_rdiscretisation%p_rtriangulation%h_DelementVolume, &
                                     p_DelementArea)

        ! Get the vector data of rx
        call lsyssc_getbase_double (rx,p_Ddata)

        nel = size(p_DelementArea)-1

        ! Build the integral
        !   int_Omega p dx
        ! This is approximated by
        !   dpintegral = SUM_Elements P(Element)*Volume(Element)

        dpintegral=0.0_DP
        do iel=1,nel
          dpintegral = dpintegral + p_Ddata(iel)*p_DelementArea(iel)
        end do

        ! Divide dpintegral by the volume of the domain; this gives the integral
        ! mean value of the pressure:

        C = dpintegral / p_DelementArea(nel+1)

        ! Subtract the integral mean value C of the pressure from all
        ! pressure components. Afterwards, P has integral mean value = 0.

        do iel=1,nel
          p_Ddata(iel) = p_Ddata(iel) - C
        end do

      case (EL_QP1)

        ! Ok, quadrilateral P1 element. Get from the triangulation the AREA-array for
        ! calculating a simple integral of rx:

        call storage_getbase_double (p_rdiscretisation%p_rtriangulation%h_DelementVolume, &
                                    p_DelementArea)

        ! Get the vector data of rx
        call lsyssc_getbase_double (rx,p_Ddata)

        nel = size(p_DelementArea)-1

        ! Build the integral
        !   int_Omega p dx
        ! This is approximated by
        !   dpintegral = SUM_Elements P(Element)*Volume(Element)
        ! Because taking the value in the element midpoint approximates the
        ! integral exactly by means of the 1x1 Gauss formula, the implementation
        ! is the same as in the Q0 and P0 case, respectively.

        dpintegral=0.0_DP
        do iel=1,nel
          dpintegral = dpintegral + p_Ddata(iel)*p_DelementArea(iel)
        end do

        ! Divide dpintegral by the volume of the domain; this gives the integral
        ! mean value of the pressure:

        C = dpintegral / p_DelementArea(nel+1)

        ! Subtract the integral mean value C of the pressure from all
        ! pressure components. Afterwards, P has integral mean value = 0.

        do iel=1,nel
          p_Ddata(iel) = p_Ddata(iel) - C
        end do

      case (EL_QPW4DCP1_2D)
        ! piecewise discontinuous P1 on quads; including a fancy constraint
        ! we'll call a separate subroutine for the dirty work here
        call vecfil_zim_QPW4DCP1_2D(rx)

      case (EL_DCQP1_2D)
        ! discontinuous non-parametric P1 on quads
        call vecfil_zim_DCQP1_2D(rx)

      case (EL_DCQP2_2D)
        ! discontinuous non-parametric P2 on quads
        call vecfil_zim_DCQP2_2D(rx)

      case default

        call output_line('Unsupported discretisation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'vecfil_normaliseToL20Sca')
        call sys_halt()

      end select

    else

      call output_line('Unsupported discretisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_normaliseToL20Sca')
      call sys_halt()

    end if

  end subroutine

! *****************************************************************************

!<subroutine>

  subroutine vecfil_zim_QPW4DCP1_2D (rx)

!<description>
  ! AUXILIARY ROUTINE:
  ! This routine imposes the zero-integral-mean constraint onto a QPW4DCP1_2D FE vector.
!</description>

!<inputoutput>
  ! The vector which is to be normalised.
  type(t_vectorScalar), intent(inout), target :: rx
!</inputoutput>

!</subroutine>

  ! local variables
  real(DP) :: dmean, dvolume
  type(t_triangulation), pointer :: p_rtria
  integer(I32) :: ctrafo
  integer :: iel, j, idx
  real(DP), dimension(:), pointer :: p_Dx
  integer, dimension(:,:), pointer :: p_IvertAtElem
  real(DP), dimension(:,:), pointer :: p_Dvtx
  real(DP), dimension(2,4) :: Dpoints
  real(DP), dimension(2,5) :: Dcoords
  real(DP), dimension(4,4) :: Djac
  real(DP), dimension(4) :: Ddetj

  real(DP), parameter :: d13 = 1.0_DP / 3.0_DP
  real(DP), parameter :: d23 = 2.0_DP / 3.0_DP

    dmean = 0.0_DP
    dvolume = 0.0_DP

    ! fetch the vector's data array
    call lsyssc_getbase_double(rx, p_Dx)

    ! fetch the triangulation from the vector
    p_rtria => rx%p_rspatialDiscr%p_rtriangulation

    ! fetch the vertex coordinates and vertices-at-element array
    call storage_getbase_double2D(p_rtria%h_DvertexCoords, p_Dvtx)
    call storage_getbase_int2D(p_rtria%h_IverticesAtElement, p_IvertAtElem)

    ! fetch the trafo type
    ctrafo = elem_igetTrafoType(EL_QPW4DCP1_2D)

    ! set up the points array; these correspond to the barycenters of the four sub-triangles
    Dpoints(1,1) = 0.0_DP
    Dpoints(2,1) = -d23
    Dpoints(1,2) = d23
    Dpoints(2,2) = 0.0_DP
    Dpoints(1,3) = 0.0_DP
    Dpoints(2,3) = d23
    Dpoints(1,4) = -d23
    Dpoints(2,4) = 0.0_DP

    ! loop over all elements of the triangulation
    do iel = 1, p_rtria%NEL

      ! fetch the vertex corner coords; someone came up with the fantastic idea that the
      ! piecewise affine trafo has to have 5 vertices, so we'll need to calculate the
      ! midpoint as the fifth vertice, too...
      Dcoords(:, 5) = 0.0_DP
      do j = 1, 4
        Dcoords(1:2, j) = p_Dvtx(1:2, p_IvertAtElem(j, iel))
        Dcoords(:,5) = Dcoords(:,5) + Dcoords(:,j)
      end do
      Dcoords(:,5) = 0.25_DP * Dcoords(:,5)

      ! call the trafo to calculate the jacobian determinants
      call trafo_calctrafoabs_mult (ctrafo, 4, Dcoords, Dpoints, Djac, Ddetj)

      ! update domain volume
      dvolume = dvolume + Ddetj(1) + Ddetj(2) + Ddetj(3) + Ddetj(4)

      ! update the integral mean
      ! The following lines correspond to integration using a piecewise 1-Gauss cubature rule.
      ! The cubature points are the barycenters of the four sub-triangles, therefore by definition
      ! of the QPW4DCP1 basis functions, each basis function has a value of +/- 1/3 or 0 in each
      ! of the cubtaute points.
      j = 11*(iel-1)
      dmean = dmean + Ddetj(1)*d13*(p_Dx(j+1) + p_Dx(j+2) + p_Dx(j+9) + p_Dx(j+10) + p_Dx(j+11))
      dmean = dmean + Ddetj(2)*d13*(p_Dx(j+3) + p_Dx(j+4) + p_Dx(j+9) + p_Dx(j+10) - p_Dx(j+11))
      dmean = dmean + Ddetj(3)*d13*(p_Dx(j+5) + p_Dx(j+6) + p_Dx(j+9) - p_Dx(j+10) - p_Dx(j+11))
      dmean = dmean + Ddetj(4)*d13*(p_Dx(j+7) + p_Dx(j+8) + p_Dx(j+9) - p_Dx(j+10) + p_Dx(j+11))

    end do

    ! divide mean by domain volume
    dmean = dmean / dvolume

    ! loop over all elements
    do iel = 1, p_rtria%NEL

      ! and correct the first 9 basis functions per element
      ! The basis functions 10 and 11 have a zero coefficient for the constant one function,
      ! therefore we do not correct the coefficients corresponding to those basis functions.
      do j = 1, 9
        idx = 11*(iel-1) + j
        p_Dx(idx) = p_Dx(idx) - dmean
      end do

    end do

  end subroutine

! *****************************************************************************

!<subroutine>

  subroutine vecfil_zim_DCQP1_2D (rx)

!<description>
  ! AUXILIARY ROUTINE:
  ! This routine imposes the zero-integral-mean constraint onto a DCQP1_2D FE vector.
!</description>

!<inputoutput>
  ! The vector which is to be normalised.
  type(t_vectorScalar), intent(inout), target :: rx
!</inputoutput>

!</subroutine>

  ! local variables
  real(DP) :: dmean, dvolume, dalpha
  type(t_triangulation), pointer :: p_rtria
  integer(I32) :: ctrafo, cevalTag
  integer :: iel, i,j, idx
  real(DP), dimension(:), pointer :: p_Dx
  real(DP), dimension(2,4) :: Dpoints
  integer, dimension(1) :: IelList
  type(t_evalElementSet) :: reval
  logical, dimension(DER_MAXNDER) :: Bder
  real(DP), dimension(3,1,4,1) :: Dbas
  real(DP) :: G2

    G2 = sqrt(1.0_DP / 3.0_DP)

    dmean = 0.0_DP
    dvolume = 0.0_DP

    ! evaluate function values only
    Bder = .false.
    Bder(DER_FUNC2D) = .true.

    ! fetch the vector's data array
    call lsyssc_getbase_double(rx, p_Dx)

    ! fetch the triangulation from the vector
    p_rtria => rx%p_rspatialDiscr%p_rtriangulation

    ! fetch the trafo type
    ctrafo = elem_igetTrafoType(EL_DCQP1_2D)

    ! fetch the evaluation tag
    cevalTag = ior(elem_getEvaluationTag(EL_DCQP1_2D), &
                   EL_EVLTAG_DETJ + EL_EVLTAG_REFPOINTS)

    ! set up the points array; these correspond to the 2x2 Gauss formula
    Dpoints(1,1) = -G2
    Dpoints(2,1) = -G2
    Dpoints(1,2) =  G2
    Dpoints(2,2) = -G2
    Dpoints(1,3) = -G2
    Dpoints(2,3) =  G2
    Dpoints(1,4) =  G2
    Dpoints(2,4) =  G2

    ! loop over all elements of the triangulation
    do iel = 1, p_rtria%NEL

      IelList(1) = iel

      ! prepare the element evaluation set
      call elprep_prepareSetForEvaluation (reval, cevalTag, &
          p_rtria, IelList, ctrafo, Dpoints)

      ! update domain volume
      dvolume = dvolume + 0.25_DP * (reval%p_Ddetj(1,1) + reval%p_Ddetj(2,1) &
                                   + reval%p_Ddetj(3,1) + reval%p_Ddetj(4,1))

      ! Evaluate the element
      call elem_eval_DCQP1_2D(EL_DCQP1_2D, reval, Bder, Dbas)

      ! update the integral mean
      idx = 3*(iel-1)
      do i = 1, 4
        dalpha = 0.0_DP
        do j = 1, 3
          dalpha = dalpha + p_Dx(idx+j)*Dbas(j,DER_FUNC,i,1)
        end do
        dmean = dmean + 0.25_DP * reval%p_Ddetj(i,1) * dalpha
      end do

    end do

    ! release element set
    call elprep_releaseElementSet(reval)

    ! divide mean by domain volume
    dmean = dmean / dvolume

    ! loop over all elements
    do iel = 1, p_rtria%NEL

      ! correct the first basis function on each element
      p_Dx(3*iel-2) = p_Dx(3*iel-2) - dmean

    end do

  end subroutine

! *****************************************************************************

!<subroutine>

  subroutine vecfil_zim_DCQP2_2D (rx)

!<description>
  ! AUXILIARY ROUTINE:
  ! This routine imposes the zero-integral-mean constraint onto a DCQP2_2D FE vector.
!</description>

!<inputoutput>
  ! The vector which is to be normalised.
  type(t_vectorScalar), intent(inout), target :: rx
!</inputoutput>

!</subroutine>

  ! local variables
  real(DP) :: dmean, dvolume, dalpha
  type(t_triangulation), pointer :: p_rtria
  integer(I32) :: ctrafo, cevalTag
  integer :: iel, i,j, idx
  real(DP), dimension(:), pointer :: p_Dx
  real(DP), dimension(2,4) :: Dpoints
  integer, dimension(1) :: IelList
  type(t_evalElementSet) :: reval
  logical, dimension(DER_MAXNDER) :: Bder
  real(DP), dimension(6,1,4,1) :: Dbas
  real(DP) :: G2

    G2 = sqrt(1.0_DP / 3.0_DP)

    dmean = 0.0_DP
    dvolume = 0.0_DP

    ! evaluate function values only
    Bder = .false.
    Bder(DER_FUNC2D) = .true.

    ! fetch the vector's data array
    call lsyssc_getbase_double(rx, p_Dx)

    ! fetch the triangulation from the vector
    p_rtria => rx%p_rspatialDiscr%p_rtriangulation

    ! fetch the trafo type
    ctrafo = elem_igetTrafoType(EL_DCQP2_2D)

    ! fetch the evaluation tag
    cevalTag = ior(elem_getEvaluationTag(EL_DCQP2_2D), &
                   EL_EVLTAG_DETJ + EL_EVLTAG_REFPOINTS)

    ! set up the points array; these correspond to the 2x2 Gauss formula
    Dpoints(1,1) = -G2
    Dpoints(2,1) = -G2
    Dpoints(1,2) =  G2
    Dpoints(2,2) = -G2
    Dpoints(1,3) = -G2
    Dpoints(2,3) =  G2
    Dpoints(1,4) =  G2
    Dpoints(2,4) =  G2

    ! loop over all elements of the triangulation
    do iel = 1, p_rtria%NEL

      IelList(1) = iel

      ! prepare the element evaluation set
      call elprep_prepareSetForEvaluation (reval, cevalTag, &
          p_rtria, IelList, ctrafo, Dpoints)

      ! update domain volume
      dvolume = dvolume + 0.25_DP * (reval%p_Ddetj(1,1) + reval%p_Ddetj(2,1) &
                                   + reval%p_Ddetj(3,1) + reval%p_Ddetj(4,1))

      ! Evaluate the element
      call elem_eval_DCQP2_2D(EL_DCQP2_2D, reval, Bder, Dbas)

      ! update the integral mean
      idx = 6*(iel-1)
      do i = 1, 4
        dalpha = 0.0_DP
        do j = 1, 6
          dalpha = dalpha + p_Dx(idx+j)*Dbas(j,DER_FUNC,i,1)
        end do
        dmean = dmean + 0.25_DP * reval%p_Ddetj(i,1) * dalpha
      end do

    end do

    ! release element set
    call elprep_releaseElementSet(reval)

    ! divide mean by domain volume
    dmean = dmean / dvolume

    ! loop over all elements
    do iel = 1, p_rtria%NEL

      ! correct the first basis function on each element
      p_Dx(6*iel-5) = p_Dx(6*iel-5) - dmean

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vecfil_normaliseSmallL1To0Sca (rx)

!<description>
  ! This routine realises the 'vector sum to 0' filter.
  ! The vector rxis normalised to bring it to the vector sum = 0
  ! (which corresponds to an l1-norm = 0).
!</description>

!<inputoutput>
  ! The vector which is to be normalised.
  type(t_vectorScalar), intent(inout),target :: rx
!</inputoutput>

!</subroutine>

    real(DP), dimension(:), pointer :: p_Ddata
    integer :: ieq
    real(DP) :: dsum

    ! Sum up all components of the vector
    call lsyssc_getbase_double (rx,p_Ddata)
    dsum = 0
    do ieq = 1,rx%NEQ
      dsum = dsum + p_Ddata(ieq)
    end do

    ! Divide by NEQ to get the mean value
    dsum = dsum / real(rx%NEQ,DP)

    ! Substract this from all components; this finishes the filter.
    do ieq = 1,rx%NEQ
      p_Ddata(ieq) = p_Ddata(ieq) - dsum
    end do

  end subroutine

! ***************************************************************************
! Block vector filters
! ***************************************************************************

!<subroutine>

  subroutine vecfil_imposePressureDropBC (rx,dtimeweight,rpdbcStructure)

!<description>
  ! Implements discrete pressure drop BC`s into a block vector.
!</description>

!<input>
  ! The t_discreteBCpressureDrop that describes the discrete pressure
  ! drop BC`s
  type(t_discreteBCpressureDrop), intent(in), target  :: rpdbcStructure

  ! Weighting factor for time-dependent problems.
  ! =1.0 for stationary simulation.
  real(DP), intent(in) :: dtimeweight
!</input>

!<inputoutput>

  ! The block vector where the boundary conditions should be imposed.
  type(t_vectorblock), intent(inout), target :: rx

!</inputoutput>

!</subroutine>

    ! local variables
    integer :: j,icp
    integer :: i
    real(DP), dimension(:), pointer    :: p_vec
    integer, dimension(:), pointer :: p_idx
    real(DP), dimension(:,:), pointer    :: p_val
    integer, dimension(:), pointer :: p_Iperm

    ! If nDOF=0, there are no DOF`s the current boundary condition segment,
    ! so we do not have to do anything. Maybe the case if the user selected
    ! a boundary region that is just too small.
    if (rpdbcStructure%nDOF .eq. 0) return

    ! Get pointers to the structures. For the vector, get the pointer from
    ! the storage management.

    call storage_getbase_int(rpdbcStructure%h_IpressureDropDOFs,p_idx)
    call storage_getbase_double2d(rpdbcStructure%h_Dmodifier,p_val)

    ! Impose the DOF value directly into the vector - more precisely, into the
    ! components of the subvectors that is indexed by Icomponent(1..NDIM2D).

    if ((.not.associated(p_idx)).or.(.not.associated(p_val))) then
      call output_line('Pressure drop BC not configured!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposePressureDropBC')
      call sys_halt()
    end if

    ! Currently, only 2D is supported.
    ! First handle the X-velocity (j=1), then the Y-velocity (j=2)
    do j=1,NDIM2D

      ! Get the subvector
      icp = rpdbcStructure%Icomponents(j)
      call lsyssc_getbase_double (rx%RvectorBlock(icp), p_vec)

      if (.not.associated(p_vec)) then
        call output_line('No vector!',&
            OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposePressureDropBC')
        call sys_halt()
      end if

      ! Only handle nDOF DOF`s, not the complete array!
      ! Probably, the array is longer (e.g. has the length of the vector), but
      ! contains only some entries...

      ! Is the vector sorted?
      if (rx%RvectorBlock(j)%isortStrategy .le. 0) then
        ! No. Implement directly.
        do i=1,rpdbcStructure%nDOF
          p_vec(p_idx(i)) = p_vec(p_idx(i)) + p_val(j,i)*dtimeweight
        end do
      else
        ! Oops, vector sorted. At first get the permutation how its sorted
        ! - or more precisely, the back-permutation, as we need this one for
        ! the loop below.
        call storage_getbase_int (rx%RvectorBlock(j)%h_IsortPermutation,p_Iperm)
        p_Iperm => p_Iperm(rx%RvectorBlock(j)%NEQ+1:)

        ! And 'filter' each DOF during the boundary value implementation!
        do i=1,rpdbcStructure%nDOF
          p_vec(p_Iperm(p_idx(i))) = p_vec(p_Iperm(p_idx(i))) + &
                                     p_val(j,i) * dtimeweight
        end do
      end if

    end do

  end subroutine

! *****************************************************************************

!<subroutine>

  subroutine vecfil_imposeNLSlipDefectBC (rx,rslipBCStructure)

!<description>
  ! Implements discrete nonlinear slip BC`s into a scalar defect vector
  ! as configured in the slip BC structure.
  ! This routine performs a special filtering to the defect vector
  ! of the type $r_m := r_m - (n*r_m)*n$ as described in
  ! [Kuzmin, Turek, Haario: Finite element simulation of turbulent
  ! bubble flows in gas-liquid reactors. Technical Report 298,
  ! September 2005, Chair of mathematics III, University of Dortmund]
!</description>

!<input>
  ! The t_discreteBCSlip that describes the discrete Dirichlet BC`s
  type(t_discreteBCSlip), intent(in), target  :: rslipBCStructure
!</input>

!<inputoutput>

  ! The block vector where the boundary conditions should be imposed.
  type(t_vectorBlock), intent(inout), target :: rx

!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,idof
    real(DP), dimension(:), pointer :: p_vecX,p_vecY
    integer, dimension(:), pointer :: p_idx
    integer, dimension(:), pointer :: p_Iperm
    real(DP), dimension(:,:), pointer :: p_Dnormals
    real(DP) :: d

    ! If nDOF=0, there are no DOF`s the current boundary condition segment,
    ! so we do not have to do anything. Maybe the case if the user selected
    ! a boundary region that is just too small.
    if (rslipBCStructure%nDOF .eq. 0) return

    ! Only 2D supported at the moment
    if (rslipBCStructure%ncomponents .ne. NDIM2D) then
      call output_line('Only 2D supported!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeNLSlipDefectBC')
      call sys_halt()
    end if

    ! Only double precision vectors supported.
    if (rx%cdataType .ne. ST_DOUBLE) then
      call output_line('Only double precision supported!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeNLSlipDefectBC')
      call sys_halt()
    end if

    ! Get pointers to the structures. For the vector, get the pointer from
    ! the storage management.

    if (rslipBCStructure%h_IslipDOFs .eq. ST_NOHANDLE) then
      call output_line('Slip-BC not configured!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeNLSlipDefectBC')
      call sys_halt()
    end if

    call storage_getbase_int(rslipBCStructure%h_IslipDOFs,p_idx)

    if (.not.associated(p_idx)) then
      call output_line('Slip-BC not configured!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeNLSlipDefectBC')
      call sys_halt()
    end if

    if (rx%RvectorBlock(rslipBCStructure%Icomponents(1))%isortStrategy .ne.&
        rx%RvectorBlock(rslipBCStructure%Icomponents(2))%isortStrategy) then
      call output_line('Subectors differently sorted!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeNLSlipDefectBC')
      call sys_halt()
    end if

    call lsyssc_getbase_double ( &
           rx%RvectorBlock(rslipBCStructure%Icomponents(1)), p_vecX)
    call lsyssc_getbase_double ( &
           rx%RvectorBlock(rslipBCStructure%Icomponents(2)), p_vecY)

    if ( (.not.associated(p_vecX)) .or. (.not.associated(p_vecX)) )then
      call output_line('No vector!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeNLSlipDefectBC')
      call sys_halt()
    end if

    call storage_getbase_double2d(rslipBCStructure%h_DnormalVectors,p_Dnormals)

    ! Impose the BC-DOF`s directly - more precisely, into the
    ! components of all the subvectors.
    !
    ! Only handle nDOF DOF`s, not the complete array!
    ! Probably, the array is longer (e.g. has the length of the vector), but
    ! contains only some entries...

    ! Is the vector sorted? If yes, all vectors are sorted the same way!
    if (rx%RvectorBlock(1)%isortStrategy .le. 0) then
      ! No. Implement directly.
      do i=1,rslipBCStructure%nDOF
        ! Get the DOF:
        idof = p_idx(i)

        ! Build n*r
        d = p_Dnormals(1,i)*p_vecX(idof) + p_Dnormals(2,i)*p_vecY(idof)

        ! Compute: r := r - (n*r)*n
        p_vecX(idof) = p_vecX(idof) - d*p_Dnormals(1,i)
        p_vecY(idof) = p_vecY(idof) - d*p_Dnormals(2,i)
      end do
    else
      ! Ups, vector sorted. At first get the permutation how its sorted -
      ! or more precisely, the back-permutation, as we need this one for
      ! the loop below.
      call storage_getbase_int (rx%RvectorBlock(1)%h_IsortPermutation,p_Iperm)
      p_Iperm => p_Iperm(rx%NEQ+1:)

      ! And 'filter' each DOF during the boundary value implementation!
      do i=1,rslipBCStructure%nDOF
        ! Get the DOF:
        idof = p_Iperm(p_idx(i))

        ! Build n*r
        d = p_Dnormals(1,i)*p_vecX(idof) + p_Dnormals(2,i)*p_vecY(idof)

        ! Compute: r := r - (n*r)*n
        p_vecX(idof) = p_vecX(idof) - d*p_Dnormals(1,idof)
        p_vecY(idof) = p_vecY(idof) - d*p_Dnormals(2,idof)
      end do
    end if

  end subroutine

! *****************************************************************************

!<subroutine>

  subroutine vecfil_imposeFeastMirrorBC (rx,rfmbcStructure)

!<description>
  ! Implements discrete Feast Mirror BC`s into a scalar vector.
!</description>

!<input>
  ! The t_discreteBCFeastMirror that describes the discrete Feast Mirror BC`s
  type(t_discreteBCFeastMirror), intent(in), target  :: rfmbcStructure
!</input>

!<inputoutput>
  ! The scalar vector where the boundary conditions should be imposed.
  type(t_vectorBlock), intent(inout), target :: rx
!</inputoutput>

!</subroutine>

  ! local variables
  integer, dimension(:), pointer :: p_ImirrorDOFs
  integer :: i
  real(DP), dimension(:), pointer    :: p_Dvec
  integer, dimension(:), pointer :: p_Iperm
  real(DP) :: dmirrorWeight

  ! Impose the DOF value directly into the vector - more precisely, into the
  ! components of the subvector that is indexed by icomponent.

  if (rx%cdataType .ne. ST_DOUBLE) then
    call output_line('Matrix must be double precision!',&
        OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeFeastMirrorBC')
    call sys_halt()
  end if

  if (rfmbcStructure%icomponent .eq. 0) then
    call output_line('FMBC not configured!',&
        OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeFeastMirrorBC')
    call sys_halt()
  end if

  if (rfmbcStructure%h_ImirrorDOFs .eq. ST_NOHANDLE) then
    ! No data inside of this structure.
    ! May happen if the region is not large enough to cover at least one DOF.
    return
  end if

  ! Get the weight of the entries.
  ! =2 on finest level, =1.5 on level NLMAX-1,...
  !dmirrorWeight = 1.0_DP+REAL(4**rfmbcStructure%icoarseningLevel,DP)
  dmirrorWeight = 1.0_DP+1.0_DP*real(2**rfmbcStructure%icoarseningLevel,DP)

  ! Get the vector data
  call lsyssc_getbase_double (rx%RvectorBlock(rfmbcStructure%icomponent),p_Dvec)

  ! Get pointers to the list of DOF`s that belong to that region and have
  ! to be tackled.
  call storage_getbase_int(rfmbcStructure%h_ImirrorDOFs,p_ImirrorDOFs)

  if ((rfmbcStructure%isubtype .eq. 1) .or. (rfmbcStructure%isubtype .eq. 3)) then
    ! The DOF`s should be treated as Dirichlet-DOF`s.
    ! Only the matrix is modified according to the FEAST mirror boundary conditions!
    !
    ! For the implementation, we just set dmirrorWeight to 0.0.
    ! This clears all DOF entries and thus treats the DOF`s like Dirichlet DOF`s.
    dmirrorWeight = 0.0_DP
  end if

  ! The vector entry corresponds to the DOF. For every DOF decide on
  ! whether it is on the FEAST mirror boundary component or not.
  ! If yes, double the entry entry.

  ! Is the vector sorted?
  if (rx%RvectorBlock(rfmbcStructure%icomponent)%isortStrategy .le. 0) then

    ! Loop through the DOF`s. Each DOF gives us the number of an entry
    ! which is to be doubled.
    do i=1,size(p_ImirrorDOFs)
      p_Dvec(p_ImirrorDOFs(i)) = dmirrorWeight * p_Dvec(p_ImirrorDOFs(i))
    end do

  else

    ! Ok, vector is sorted, so we have to filter all the DOF`s through the
    ! permutation before using them for implementing boundary conditions.
    !
    ! Get the permutation (or more precisely, the inverse permutation)
    ! from the vector to renumber the columns into
    ! the actual DOF numbers.
    call storage_getbase_int (&
        rx%RvectorBlock(rfmbcStructure%icomponent)%h_IsortPermutation,p_Iperm)
    p_Iperm => p_Iperm(rx%RvectorBlock(rfmbcStructure%icomponent)%NEQ+1:)

    ! Loop through the DOF`s. Each DOF gives us the number of an entry
    ! which is to be doubled.
    do i=1,size(p_ImirrorDOFs)
      p_Dvec(p_Iperm(p_ImirrorDOFs(i))) = dmirrorWeight * p_Dvec(p_Iperm(p_ImirrorDOFs(i)))
    end do

  end if

  end subroutine

! *****************************************************************************

!<subroutine>

  subroutine vecfil_imposeFeastMirrorDefBC (rx,rfmbcStructure)

!<description>
  ! Implements discrete Feast Mirror BC`s into a scalar defect vector.
!</description>

!<input>
  ! The t_discreteBCFeastMirror that describes the discrete Feast Mirror BC`s
  type(t_discreteBCFeastMirror), intent(in), target  :: rfmbcStructure
!</input>

!<inputoutput>
  ! The scalar vector where the boundary conditions should be imposed.
  type(t_vectorBlock), intent(inout), target :: rx
!</inputoutput>

!</subroutine>

  ! local variables
  integer, dimension(:), pointer :: p_ImirrorDOFs
  integer :: i
  real(DP), dimension(:), pointer    :: p_Dvec
  integer, dimension(:), pointer :: p_Iperm

  ! Impose the DOF value directly into the vector - more precisely, into the
  ! components of the subvector that is indexed by icomponent.

  if (rx%cdataType .ne. ST_DOUBLE) then
    call output_line('Matrix must be double precision!',&
        OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeFeastMirrorDefBC')
    call sys_halt()
  end if

  if (rfmbcStructure%icomponent .eq. 0) then
    call output_line('FMBC not configured!',&
        OU_CLASS_ERROR,OU_MODE_STD,'vecfil_imposeFeastMirrorDefBC')
    call sys_halt()
  end if

  if (rfmbcStructure%h_ImirrorDOFs .eq. ST_NOHANDLE) then
    ! No data inside of this structure.
    ! May happen if the region is not large enough to cover at least one DOF.
    return
  end if

  ! Get the vector data
  call lsyssc_getbase_double (rx%RvectorBlock(rfmbcStructure%icomponent),p_Dvec)

  ! Get pointers to the list of DOF`s that belong to that region and have
  ! to be tackled.
  call storage_getbase_int(rfmbcStructure%h_ImirrorDOFs,p_ImirrorDOFs)

  if ((rfmbcStructure%isubtype .eq. 1) .or. (rfmbcStructure%isubtype .eq. 3)) then
    ! The DOF`s should be treated as Dirichlet-DOF`s.
    !
    ! The vector entry corresponds to the DOF. For every DOF decide on
    ! whether it is on the FEAST mirror boundary component or not.
    ! If yes, put the entry to zero.

    ! Is the vector sorted?
    if (rx%RvectorBlock(rfmbcStructure%icomponent)%isortStrategy .le. 0) then

      ! Loop through the DOF`s. Each DOF gives us the number of an entry
      ! which is to be doubled.
      do i=1,size(p_ImirrorDOFs)
        p_Dvec(p_ImirrorDOFs(i)) = 0.0_DP
      end do

    else

      ! Ok, vector is sorted, so we have to filter all the DOF`s through the
      ! permutation before using them for implementing boundary conditions.
      !
      ! Get the permutation (or more precisely, the inverse permutation)
      ! from the vector to renumber the columns into
      ! the actual DOF numbers.
      call storage_getbase_int (&
          rx%RvectorBlock(rfmbcStructure%icomponent)%h_IsortPermutation,p_Iperm)
      p_Iperm => p_Iperm(rx%RvectorBlock(rfmbcStructure%icomponent)%NEQ+1:)

      ! Loop through the DOF`s. Each DOF gives us the number of an entry
      ! which is to be doubled.
      do i=1,size(p_ImirrorDOFs)
        p_Dvec(p_Iperm(p_ImirrorDOFs(i))) = 0.0_DP
      end do

    end if

  end if

  end subroutine

  ! ***************************************************************************
  ! Implementation of discrete boundary conditions into block solution vectors
  ! ***************************************************************************

!<subroutine>

  subroutine vecfil_discreteBCsol (rx,rdiscreteBC)

!<description>
  ! This routine realises the 'impose discrete boundary conditions to solution'
  ! filter. This filter imposes the discrete boundary conditions rdiscreteBC
  ! (if specified) or (if rdiscreteBC is not specified) the boundary conditions
  ! which are  associated to the vector rx (with rx%p_discreteBC) to this
  ! (block) vector.
!</description>

!<input>
  ! OPTIONAL: boundary conditions to impose into the vector.
  ! If not specified, the default boundary conditions associated to the
  ! vector rx are imposed to the matrix.
  type(t_discreteBC), optional, intent(in), target :: rdiscreteBC
!</input>

!<inputoutput>
  ! The block vector where the boundary conditions should be imposed.
  type(t_vectorBlock), intent(inout),target :: rx
!</inputoutput>

!</subroutine>

    integer :: iblock,i
    type(t_discreteBC), pointer :: p_rdiscreteBC

    if (.not. present(rdiscreteBC)) then
      ! Grab the boundary condition entry list from the vector. This
      ! is a list of all discretised boundary conditions in the system.
      if (associated(rx%p_rdiscreteBC)) then
        p_rdiscreteBC => rx%p_rdiscreteBC
      else
        ! No BC
        nullify(p_RdiscreteBC)
      end if
    else
      p_RdiscreteBC => rdiscreteBC
    end if

    if (.not. associated(p_RdiscreteBC)) return

    ! Now loop through all entries in this list:
    !DO i=1,SIZE(p_RdiscreteBC)
    do i=1, p_rdiscreteBC%inumEntriesUsed

      ! What for BC`s do we have here?
      select case (p_rdiscreteBC%p_RdiscBCList(i)%itype)
      case (DISCBC_TPUNDEFINED)
        ! Do-nothing

      case (DISCBC_TPDIRICHLET)
        ! Dirichlet boundary conditions. Not time-dependent.
        ! On which component are they defined? The component specifies
        ! the row of the block matrix that is to be altered.
        iblock = p_rdiscreteBC%p_RdiscBCList(i)%rdirichletBCs%icomponent

        ! Implement the Dirichlet boundary conditions into that component
        ! of the vector.
        call vecfil_imposeDirichletBC (rx%RvectorBlock(iblock),&
            p_rdiscreteBC%p_RdiscBCList(i)%rdirichletBCs)

      case (DISCBC_TPPRESSUREDROP)
        ! Nothing to do; pressure drop BC`s are implemented only into the RHS.

      case (DISCBC_TPSLIP)
        ! Nothing to do

      case (DISCBC_TPFEASTMIRROR)
        ! Nothing to do

      case default
        call output_line(&
            'Unknown boundary condition:'//&
            sys_siL(p_rdiscreteBC%p_RdiscBCList(i)%itype,5),&
            OU_CLASS_ERROR,OU_MODE_STD,'vecfil_discreteBCsol')

        call sys_halt()

      end select
    end do

  end subroutine

  ! ***************************************************************************
  ! Implementation of discrete boundary conditions into block RHS vectors
  ! ***************************************************************************

!<subroutine>

  subroutine vecfil_discreteBCrhs (rx,rdiscreteBC)

!<description>
  ! This routine realises the 'impose discrete boundary conditions to RHS'
  ! filter. This filter imposes the discrete boundary conditions rdiscreteBC
  ! (if specified) or (if rdiscreteBC is not specified) the boundary conditions
  ! which are  associated to the vector rx (with rx%p_discreteBC) to this
  ! (block) vector.
!</description>

!<input>
  ! OPTIONAL: boundary conditions to impose into the vector.
  ! If not specified, the default boundary conditions associated to the
  ! vector rx are imposed to the matrix.
  type(t_discreteBC), optional, intent(in), target :: rdiscreteBC
!</input>

!<inputoutput>
  ! The block vector where the boundary conditions should be imposed.
  type(t_vectorBlock), intent(inout),target :: rx
!</inputoutput>

!</subroutine>

    integer :: iblock,i
    type(t_discreteBC), pointer :: p_rdiscreteBC

    if (.not. present(rdiscreteBC)) then
      ! Grab the boundary condition entry list from the vector. This
      ! is a list of all discretised boundary conditions in the system.
      p_RdiscreteBC => rx%p_rdiscreteBC
    else
      p_RdiscreteBC => rdiscreteBC
    end if

    if (.not. associated(p_RdiscreteBC)) return

    ! Now loop through all entries in this list:
    !DO i=1,SIZE(p_RdiscreteBC)
    do i=1, p_rdiscreteBC%inumEntriesUsed

      ! What for BC`s do we have here?
      select case (p_rdiscreteBC%p_RdiscBCList(i)%itype)
      case (DISCBC_TPUNDEFINED)
        ! Do-nothing

      case (DISCBC_TPDIRICHLET)
        ! Dirichlet boundary conditions. Not time-dependent.
        ! On which component are they defined? The component specifies
        ! the row of the block matrix that is to be altered.
        iblock = p_rdiscreteBC%p_RdiscBCList(i)%rdirichletBCs%icomponent

        ! Implement the Dirichlet boundary conditions into that component
        ! of the vector.
        call vecfil_imposeDirichletBC (rx%RvectorBlock(iblock),&
                                      p_rdiscreteBC%p_RdiscBCList(i)%rdirichletBCs)

      case (DISCBC_TPPRESSUREDROP)
        ! Nothing to do.

      case (DISCBC_TPSLIP)
        ! Nothing to do.

      case (DISCBC_TPFEASTMIRROR)
        call vecfil_imposeFeastMirrorBC (rx,&
            p_rdiscreteBC%p_RdiscBCList(i)%rfeastMirrorBCs)

      case default
        call output_line(&
            'Unknown boundary condition:'//&
            sys_siL(p_rdiscreteBC%p_RdiscBCList(i)%itype,5),&
            OU_CLASS_ERROR,OU_MODE_STD,'vecfil_discreteBCrhs')

        call sys_halt()

      end select
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vecfil_discreteNLPDropBCrhs (rx,dtimeWeight,rdiscreteBC)

!<description>
  ! Implements discrete pressure drop BC`s into a block RHS vector.
  ! To be called inside of a nonlinear or time-stepping loop.
  ! This routine performs a special filtering to the RHS vector
  ! of the type $r_m := r_m - \sum_j P_j \int_{S_j} \phi n ds$ as described
  ! on page 257 (235) Turek`s book.
  !
  ! The filtering is applied to the boundary components configured
  ! for pressure drop when setting up the BC`s.
!</description>

!<input>
  ! OPTIONAL: Time-step weight. This weight is multiplied to time-dependent
  ! boundary conditions before these are added to the vector rx.
  ! The parameter can be omitted in stationary simulations.
  real(DP), intent(in), optional :: dtimeWeight

  ! OPTIONAL: boundary conditions to impose into the vector.
  ! If not specified, the default boundary conditions associated to the
  ! vector rx are imposed to the vector.
  type(t_discreteBC), optional, intent(in), target :: rdiscreteBC
!</input>

!<inputoutput>
  ! The block vector where the boundary conditions should be imposed.
  type(t_vectorBlock), intent(inout),target :: rx
!</inputoutput>

!</subroutine>

    integer :: i
    real(DP) :: dtweight
    type(t_discreteBC), pointer :: p_rdiscreteBC

    if (.not. present(rdiscreteBC)) then
      ! Grab the boundary condition entry list from the vector. This
      ! is a list of all discretised boundary conditions in the system.
      if (associated(rx%p_rdiscreteBC)) then
        p_RdiscreteBC => rx%p_rdiscreteBC
      else
        ! No BC
        nullify(p_RdiscreteBC)
      end if
    else
      p_RdiscreteBC => rdiscreteBC
    end if

    if (.not. associated(p_RdiscreteBC)) return

    ! If the time-weight is not specified, 1.0 is assumed.
    if (present(dtimeWeight)) then
      dtweight = dtimeWeight
    else
      dtweight = 1.0_DP
    end if
    ! Note: Time-step weight not used by any filter up to now!
    ! Perhaps in a later implementation it is needed anywhere...

    ! Now loop through all entries in this list:
    !DO i=1,SIZE(p_RdiscreteBC)
    do i=1, p_rdiscreteBC%inumEntriesUsed

      ! Only implement discrete pressure drop BC`s.
      if (p_rdiscreteBC%p_RdiscBCList(i)%itype .eq. DISCBC_TPPRESSUREDROP) then
        call vecfil_imposePressureDropBC (rx,dtweight,&
            p_rdiscreteBC%p_RdiscBCList(i)%rpressureDropBCs)
      end if

    end do

  end subroutine

  ! ***************************************************************************
  ! Implementation of discrete boundary conditions into block defect vectors
  ! ***************************************************************************

!<subroutine>

  subroutine vecfil_discreteBCdef (rx,rdiscreteBC)

!<description>
  ! This routine realises the 'impose discrete boundary conditions to defect'
  ! filter. This filter imposes the discrete boundary conditions rdiscreteBC
  ! (if specified) or (if rdiscreteBC is not specified) the boundary conditions
  ! which are  associated to the defect vector rx (with rx%p_discreteBC) to
  ! this (block) vector.
!</description>

!<input>
  ! OPTIONAL: boundary conditions to impose into the vector.
  ! If not specified, the default boundary conditions associated to the
  ! vector rx are imposed to the matrix.
  type(t_discreteBC), optional, intent(in), target :: rdiscreteBC
!</input>

!<inputoutput>
  ! The block vector where the boundary conditions should be imposed.
  type(t_vectorBlock), intent(inout),target :: rx
!</inputoutput>

!</subroutine>

    integer :: iblock,i !,icp
    type(t_discreteBC), pointer :: p_rdiscreteBC

    if (.not. present(rdiscreteBC)) then
      ! Grab the boundary condition entry list from the vector. This
      ! is a list of all discretised boundary conditions in the system.
      if (associated(rx%p_rdiscreteBC)) then
        p_rdiscreteBC => rx%p_rdiscreteBC
      else
        ! No BC
        nullify(p_rdiscreteBC)
      end if
    else
      p_rdiscreteBC => rdiscreteBC
    end if

    if (.not. associated(p_rdiscreteBC)) return

    ! Now loop through all entries in this list:
    !DO i=1,SIZE(p_RdiscreteBC)
    do i=1, p_rdiscreteBC%inumEntriesUsed

      ! What for BC`s do we have here?
      select case (p_rdiscreteBC%p_RdiscBCList(i)%itype)
      case (DISCBC_TPUNDEFINED)
        ! Do-nothing

      case (DISCBC_TPDIRICHLET)
        ! Dirichlet boundary conditions. Not time-dependent.
        ! On which component are they defined? The component specifies
        ! the row of the block matrix that is to be altered.
        iblock = p_rdiscreteBC%p_RdiscBCList(i)%rdirichletBCs%icomponent

        ! Implement the Dirichlet boundary conditions into that component
        ! of the vector.
        call vecfil_imposeDirichletDefectBC (rx%RvectorBlock(iblock),&
            p_rdiscreteBC%p_RdiscBCList(i)%rdirichletBCs)

      case (DISCBC_TPPRESSUREDROP)
        ! Nothing to do; pressure drop BC`s are implemented only into the RHS.

      case (DISCBC_TPSLIP)
        ! Slip boundary conditions in the linear case are implemented
        ! in a nonlinear loop - so there is nothing to do here.

      case (DISCBC_TPFEASTMIRROR)
        ! Routine is on purpose not commented in! Not used for now!
        ! CALL vecfil_imposeFeastMirrorDefBC (rx,p_RdiscreteBC(i)%rfeastMirrorBCs)

      case default
        call output_line(&
            'Unknown boundary condition:'//&
            sys_siL(p_rdiscreteBC%p_RdiscBCList(i)%itype,5),&
            OU_CLASS_ERROR,OU_MODE_STD,'vecfil_discreteBCdef')
        call sys_halt()

      end select
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vecfil_discreteNLSlipBCdef (rx,rdiscreteBC)

!<description>
  ! Implements discrete nonlinear slip BC`s into a defect vector.
  ! Nonlinear filter, to be called inside of a nonlinear loop.
  ! This routine performs a special filtering to the defect vector
  ! of the type $r_m := r_m - (n*r_m)*n$ as described in
  ! [Kuzmin, Turek, Haario: Finite element simulation of turbulent
  ! bubble flows in gas-liquid reactors. Technical Report 298,
  ! September 2005, Chair of mathematics III, University of Dortmund]
  !
  ! The filtering is applied to all boundary components configured
  ! as slip in rx or rdiscreteBC (if given), respectively.
!</description>

!<input>
  ! OPTIONAL: boundary conditions to impose into the vector.
  ! If not specified, the default boundary conditions associated to the
  ! vector rx are imposed to the vector.
  type(t_discreteBC), optional, intent(in), target :: rdiscreteBC
!</input>

!<inputoutput>
  ! The block vector where the boundary conditions should be imposed.
  type(t_vectorBlock), intent(inout),target :: rx
!</inputoutput>

!</subroutine>

    integer :: i
    type(t_discreteBC), pointer :: p_rdiscreteBC

    if (.not. present(rdiscreteBC)) then
      ! Grab the boundary condition entry list from the vector. This
      ! is a list of all discretised boundary conditions in the system.
      if (associated(rx%p_rdiscreteBC)) then
        p_rdiscreteBC => rx%p_rdiscreteBC
      else
        ! No BC
        nullify(p_rdiscreteBC)
      end if
    else
      p_rdiscreteBC => rdiscreteBC
    end if

    if (.not. associated(p_rdiscreteBC)) return

    ! Now loop through all entries in this list:
    !DO i=1,SIZE(p_RdiscreteBC)
    do i=1, p_rdiscreteBC%inumEntriesUsed

      ! Only implement discrete slip BC`s.
      if (p_rdiscreteBC%p_RdiscBCList(i)%itype .eq. DISCBC_TPSLIP) then
        call vecfil_imposeNLSlipDefectBC (rx,p_rdiscreteBC%p_RdiscBCList(i)%rslipBCs)
      end if

    end do

  end subroutine

  ! ***************************************************************************
  ! Implementation of discrete fictitious boundary conditions into
  ! block solution vectors
  ! ***************************************************************************

!<subroutine>

  subroutine vecfil_discreteFBCsol (rx,rdiscreteFBC)

!<description>
  ! This routine realises the `impose discrete fictitious boundary conditions
  ! to solution` filter.
  ! This filter imposes the discrete fictitious boundary conditions rdiscreteFBC
  ! (if specified) or (if rdiscreteFBC is not specified) the boundary conditions
  ! which are  associated to the vector rx (with rx%p_discreteBC) to this
  ! (block) vector.
!</description>

!<input>
  ! OPTIONAL: boundary conditions to impose into the vector.
  ! If not specified, the default fictitious boundary conditions associated
  ! to the vector rx are imposed to the matrix.
  type(t_discreteFBC), optional, intent(in), target :: rdiscreteFBC
!</input>

!<inputoutput>
  ! The block vector where the boundary conditions should be imposed.
  type(t_vectorBlock), intent(inout),target :: rx
!</inputoutput>

!</subroutine>

    integer :: iblock,i,j
    type(t_discreteFBCEntry), dimension(:), pointer :: p_RdiscreteFBC

    if (.not. present(rdiscreteFBC)) then
      ! Grab the boundary condition entry list from the vector. This
      ! is a list of all discretised boundary conditions in the system.
      if (associated(rx%p_rdiscreteBCfict)) then
        p_RdiscreteFBC => rx%p_rdiscreteBCfict%p_RdiscFBCList
      else
        ! No BC
        nullify(p_RdiscreteFBC)
      end if
    else
      p_RdiscreteFBC => rdiscreteFBC%p_RdiscFBCList
    end if

    if (.not. associated(p_RdiscreteFBC)) return

    ! Now loop through all entries in this list:
    do i=1,size(p_RdiscreteFBC)

      ! What for BC`s do we have here?
      select case (p_RdiscreteFBC(i)%itype)
      case (DISCFBC_TPUNDEFINED)
        ! Do-nothing

      case (DISCFBC_TPDIRICHLET)
        ! Dirichlet boundary conditions. Not time-dependent.
        ! Loop through all components, these boundary conditions should apply to.

        do j=1,p_RdiscreteFBC(i)%rdirichletFBCs%ncomponents
          iblock = p_RdiscreteFBC(i)%rdirichletFBCs%Icomponents(j)
          ! Implement the Dirichlet boundary conditions into that component
          ! of the vector.
          call vecfil_imposeDirichletFBC (rx%RvectorBlock(iblock),j,&
                                          p_RdiscreteFBC(i)%rdirichletFBCs)
        end do

      case (DISCBC_TPSLIP)
        ! Nothing to do.

      case (DISCBC_TPFEASTMIRROR)
        ! Nothing to do

      case default
        call output_line(&
            'Unknown boundary condition:'//sys_siL(p_RdiscreteFBC(i)%itype,5),&
            OU_CLASS_ERROR,OU_MODE_STD,'vecfil_discreteBCsol')
        call sys_halt()

      end select
    end do

  end subroutine

  ! ***************************************************************************
  ! Implementation of discrete fictitious boundary conditions into
  ! block RHS vectors
  ! ***************************************************************************

!<subroutine>

  subroutine vecfil_discreteFBCrhs (rx,rdiscreteFBC)

!<description>
  ! This routine realises the `impose discrete fictitious boundary conditions
  ! to RHS` filter.
  ! This filter imposes the discrete fictitious boundary conditions rdiscreteFBC
  ! (if specified) or (if rdiscreteFBC is not specified) the boundary conditions
  ! which are  associated to the vector rx (with rx%p_discreteFBC) to this
  ! (block) vector.
!</description>

!<input>
  ! OPTIONAL: boundary conditions to impose into the vector.
  ! If not specified, the default boundary conditions associated to the
  ! vector rx are imposed to the matrix.
  type(t_discreteFBC), optional, intent(in), target :: rdiscreteFBC
!</input>

!<inputoutput>
  ! The block vector where the boundary conditions should be imposed.
  type(t_vectorBlock), intent(inout),target :: rx
!</inputoutput>

!</subroutine>

    integer :: iblock,i,j
    type(t_discreteFBCEntry), dimension(:), pointer :: p_RdiscreteFBC

    if (.not. present(rdiscreteFBC)) then
      ! Grab the boundary condition entry list from the vector. This
      ! is a list of all discretised boundary conditions in the system.
      if (associated(rx%p_rdiscreteBCfict)) then
        p_RdiscreteFBC => rx%p_rdiscreteBCfict%p_RdiscFBCList
      else
        ! No BC
        nullify(p_RdiscreteFBC)
      end if
    else
      p_RdiscreteFBC => rdiscreteFBC%p_RdiscFBCList
    end if

    if (.not. associated(p_RdiscreteFBC)) return

    ! Now loop through all entries in this list:
    do i=1,size(p_RdiscreteFBC)

      ! What for BC`s do we have here?
      select case (p_RdiscreteFBC(i)%itype)
      case (DISCFBC_TPUNDEFINED)
        ! Do-nothing

      case (DISCFBC_TPDIRICHLET)
        ! Loop through all components, these boundary conditions should apply to.

        do j=1,p_RdiscreteFBC(i)%rdirichletFBCs%ncomponents
          iblock = p_RdiscreteFBC(i)%rdirichletFBCs%Icomponents(j)

          ! Implement the Dirichlet boundary conditions into that component
          ! of the vector.
          call vecfil_imposeDirichletFBC (rx%RvectorBlock(iblock),j,&
                                          p_RdiscreteFBC(i)%rdirichletFBCs)
        end do

      case (DISCBC_TPSLIP)
        ! Nothing to do.

      case (DISCBC_TPFEASTMIRROR)
        ! Nothing to do

      case default
        call output_line(&
            'Unknown boundary condition:'//sys_siL(p_RdiscreteFBC(i)%itype,5),&
            OU_CLASS_ERROR,OU_MODE_STD,'vecfil_discreteFBCrhs')
        call sys_halt()

      end select
    end do

  end subroutine

  ! ***************************************************************************
  ! Implementation of discrete fictitious boundary conditions into
  ! block defect vectors
  ! ***************************************************************************

!<subroutine>

  subroutine vecfil_discreteFBCdef (rx,rdiscreteFBC)

!<description>
  ! This routine realises the `impose discrete fictitious boundary conditions
  ! to defect` filter.
  ! This filter imposes the discrete fictitious boundary conditions rdiscretFeBC
  ! (if specified) or (if rdiscreteFBC is not specified) the boundary conditions
  ! which are  associated to the defect vector rx (with rx%p_discreteFBC) to
  ! this (block) vector.
!</description>

!<input>
  ! OPTIONAL: boundary conditions to impose into the vector.
  ! If not specified, the default boundary conditions associated to the
  ! vector rx are imposed to the matrix.
  type(t_discreteFBC), optional, intent(in), target :: rdiscreteFBC
!</input>

!<inputoutput>
  ! The block vector where the boundary conditions should be imposed.
  type(t_vectorBlock), intent(inout),target :: rx
!</inputoutput>

!</subroutine>

    integer :: iblock,i,j
    type(t_discreteFBCEntry), dimension(:), pointer :: p_RdiscreteFBC

    if (.not. present(rdiscreteFBC)) then
      ! Grab the boundary condition entry list from the vector. This
      ! is a list of all discretised boundary conditions in the system.
      if (associated(rx%p_rdiscreteBCfict)) then
        p_RdiscreteFBC => rx%p_rdiscreteBCfict%p_RdiscFBCList
      else
        ! No BC.
        nullify(p_RdiscreteFBC)
      end if
    else
      p_RdiscreteFBC => rdiscreteFBC%p_RdiscFBCList
    end if

    if (.not. associated(p_RdiscreteFBC)) return

    ! Now loop through all entries in this list:
    do i=1,size(p_RdiscreteFBC)

      ! What for BC`s do we have here?
      select case (p_RdiscreteFBC(i)%itype)
      case (DISCFBC_TPUNDEFINED)
        ! Do-nothing

      case (DISCFBC_TPDIRICHLET)
        ! Dirichlet boundary conditions. Not time-dependent.

        ! Loop over all blocks where to implement these FBC`s.
        do j=1,p_RdiscreteFBC(i)%rdirichletFBCs%ncomponents
          iblock = p_RdiscreteFBC(i)%rdirichletFBCs%Icomponents(j)

          ! Implement the Dirichlet boundary conditions into that component
          ! of the vector.
          call vecfil_imposeDirichletDefectFBC (rx%RvectorBlock(iblock),&
                                                p_RdiscreteFBC(i)%rdirichletFBCs)
        end do

      case (DISCBC_TPSLIP)
        ! Nothing to do.

      case (DISCBC_TPFEASTMIRROR)
        ! Nothing to do

      case default
        call output_line(&
            'Unknown boundary condition:'//sys_siL(p_RdiscreteFBC(i)%itype,5),&
            OU_CLASS_ERROR,OU_MODE_STD,'vecfil_discreteFBCdef')
        call sys_halt()

      end select
    end do

  end subroutine

  ! ***************************************************************************
  ! Other block filters
  ! ***************************************************************************

!<subroutine>

  subroutine vecfil_subvectorToL20 (rx,isubvector)

!<description>
  ! This routine realises the 'subvector to <tex>$L^2_0$</tex>' filter.
  ! The subvector isubvector of the block vector rx is normalised
  ! with vecfil_normaliseScalarToL20 to bring it to the space <tex>$L^2_0$</tex>.
!</description>

!<inputoutput>

  ! The block vector which is partially to be normalised.
  type(t_vectorBlock), intent(inout),target :: rx

  ! The number of the subvector in rx which is to be normalised.
  integer, intent(in) :: isubvector

!</inputoutput>

!</subroutine>

    if ((isubvector .le. 0) .or. (isubvector .gt. size(rx%RvectorBlock))) then
      call output_line('isubvector out of allowed range!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_subvectorToL20')
      call sys_halt()
    end if

    ! Do not use one single IF here as this may lead to errors depending
    ! on the compiler (subvector=0 and access to vector(isubvector)).
    if ((rx%RvectorBlock(isubvector)%h_Ddata .eq. ST_NOHANDLE) .or. &
        (rx%RvectorBlock(isubvector)%NEQ .le. 0)) &
      return

    ! Normalise the subvector isubvector
    call vecfil_normaliseToL20Sca (rx%RvectorBlock(isubvector))

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vecfil_subvectorSmallL1To0 (rx,isubvector)

!<description>
  ! This routine realises the 'vector sum to 0' filter.
  ! The subvector isubvector of the block vector rx is normalised
  ! with vecfil_normaliseSmallL1To0Sca to bring it to the vector sum = 0
  ! (which corresponds to an l1-norm = 0).
!</description>

!<inputoutput>

  ! The block vector which is partially to be normalised.
  type(t_vectorBlock), intent(inout),target :: rx

  ! The number of the subvector in rx which is to be normalised.
  integer, intent(in) :: isubvector

!</inputoutput>

!</subroutine>

    if ((isubvector .le. 0) .or. (isubvector .gt. size(rx%RvectorBlock))) then
      call output_line('isubvector out of allowed range!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vecfil_normaliseSmallL1To0')
      call sys_halt()
    end if

    ! Do not use one single IF here as this may lead to errors depending
    ! on the compiler (subvector=0 and access to vector(isubvector)).
    if ((rx%RvectorBlock(isubvector)%h_Ddata .eq. ST_NOHANDLE) .or. &
        (rx%RvectorBlock(isubvector)%NEQ .le. 0)) &
      return

    ! Normalise the subvector isubvector
    call vecfil_normaliseSmallL1To0Sca (rx%RvectorBlock(isubvector))

  end subroutine

end module
