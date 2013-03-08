!##############################################################################
!# ****************************************************************************
!# <name> vanka_navst2d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the vanka driver for 2D Navier-Stokes systems.
!# The most general version of a 2D Navier-System looks as follows:
!#
!# <verb>
!#                     / A11 A12 B1 \
!#                     | A21 A22 B2 |
!#                     \ D1  D2  C  /
!# </verb>
!#
!# The following matrices are optional: A12, A21, C
!# It is silently assumed that the following conditions hold:
!#
!# 1.) A11 and A22 have the same matrix structure.
!#
!# 2.) Either both A12 and A21 exist or none of them.
!#     If they exist, then they have the same matrix structure.
!#
!# 3.) B1 and B2 have the same matrix structure.
!#
!# 4.) D1 and D2 have the same matrix structure.
!#
!# Please note that in contrast to the "old" Navier-Stokes Vanka methods,
!# the D-matrices must NOT be virtually transposed!
!#
!#
!# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\\
!# How do I add Vanka support for a new discretisation? \\
!# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!#
!# Unless you have some extraordinary wishes, this job is pretty easy
!#
!# At the bottom of this file you will find two templates: one for the "full"
!# and one for the "diagonal" Vanka variant. Search this file for $TEMPLATE$
!# to get to the template section. There is also a small description what
!# needs to be done with that template.
!#
!# After you have converted the template to an implementation (or wrote your
!# own implementation without using the template) you will need to modify
!# the vanka_init_NavSt2D and vanka_solve_NavSt2D routines to use your new
!# Vanka implementation (search the routines for $TODO$).
!#
!# Once you are finished with that, Vanka should be able to work with your
!# new implementation.
!#
!# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\\
!# Solving the local Navier-Stokes System \\
!# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!#
!# Our local Navier-Stokes system looks as follows:
!#
!# <verb>
!#                   / A11 A12 B1 \   / u1 \   / f_u1 \
!#                   | A21 A22 B2 | * | u2 | = | f_u2 |
!#                   \ D1  D2  C  /   \ p  /   \ f_p  /
!# </verb>
!#
!# First, we will rewrite our system:
!#
!# <verb>
!# A := / A11 A12 \     B := / B1 \     D := ( D1 D2 )
!#      \ A21 A22 /          \ B2 /
!#
!# u := / u1 \     f_u := / f_u1 \
!#      \ u1 /            \ f_u2 /
!# </verb>
!#
!# Now our system looks as follows:
!#
!# <verb>
!#                   / A B \ * / u \ = / f_u \
!#                   \ D C /   \ p /   \ f_p /
!# </verb>
!#
!#
!# 1.) Calculate the Schur-Complement of A:
!#
!#     <tex> $$          S := C - D * A^{-1} * B   $$ </tex>
!#
!# 2.) Calculate pressure:
!#
!#     <tex> $$          p := S^{-1} * (f_p - D * A^{-1} * f_u)  $$ </tex>
!#
!# 3.) Calculate velocity:
!#
!#     <tex> $$          u := A^{-1} * (f_u - B * p)  $$ </tex>
!#
!# </purpose>
!##############################################################################

module vanka_navst2d

!$use omp_lib
  use fsystem
  use storage
  use genoutput
  use mprimitives
  use element
  use triangulation
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use dofmapping

  implicit none

  private

  public :: t_vanka_NavSt2D
  public :: vanka_init_NavSt2D
  public :: vanka_solve_NavSt2D

!<constants>

!<constantblock description="Vanka type identifiers for the 2D Navier-Stokes Class">

  ! Diagonal-type VANKA
  integer, parameter, public :: VANKATP_NAVST2D_DIAG = 0

  ! "Full" VANKA
  integer, parameter, public :: VANKATP_NAVST2D_FULL = 1

!</constantblock>
!</constants>


!<types>

!<typeblock>

  ! A structure that saves matrix pointers for the 2D Navier-Stokes Vanka driver.
  type t_vanka_NavSt2D

    private

    ! The Vanka subtype that is to be used
    integer :: csubtype = VANKATP_NAVST2D_DIAG

    ! Pointer to the column structure of the velocity matrix A11/A22
    integer, dimension(:), pointer :: p_KcolA => null()

    ! Pointer to the row structure of the velocity matrix A11/A22
    integer, dimension(:), pointer :: p_KldA => null()

    ! Pointer to diagonal entries in the velocity matrix A11/A22
    integer, dimension(:), pointer :: p_KdiagonalA => null()

    ! Pointer to the matrix entries of the velocity matrix A11
    real(DP), dimension(:), pointer :: p_DA11 => null()

    ! Pointer to the matrix entries of the velocity matrix A22
    real(DP), dimension(:), pointer :: p_DA22 => null()

    ! Pointer to the column structure of the matrix A12/A21
    integer, dimension(:), pointer :: p_KcolA12 => null()

    ! Pointer to the row structure of the matrix A12/A21
    integer, dimension(:), pointer :: p_KldA12 => null()

    ! Pointer to the matrix entries of the velocity matrix A12
    real(DP), dimension(:), pointer :: p_DA12 => null()

    ! Pointer to the matrix entries of the velocity matrix A21
    real(DP), dimension(:), pointer :: p_DA21 => null()

    ! Pointer to the column structure of the B-matrices.
    integer, dimension(:), pointer :: p_KcolB => null()

    ! Pointer to the row structure of the B-matrices
    integer, dimension(:), pointer :: p_KldB => null()

    ! Pointer to the entries of the B1-matrix
    real(DP), dimension(:), pointer :: p_DB1 => null()

    ! Pointer to the entries of the B2-matrix
    real(DP), dimension(:), pointer :: p_DB2 => null()

    ! Pointer to the column structure of the D-matrices.
    integer, dimension(:), pointer :: p_KcolD => null()

    ! Pointer to the row structure of the D-matrices
    integer, dimension(:), pointer :: p_KldD => null()

    ! Pointer to the entries of the D1-matrix
    real(DP), dimension(:), pointer :: p_DD1 => null()

    ! Pointer to the entries of the D2-matrix
    real(DP), dimension(:), pointer :: p_DD2 => null()

    ! Pointer to the column structure of the C-matrix.
    integer, dimension(:), pointer :: p_KcolC => null()

    ! Pointer to the row structure of the C-matrix.
    integer, dimension(:), pointer :: p_KldC => null()

    ! Pointer to diagonal entries in the C-matrix
    integer, dimension(:), pointer :: p_KdiagonalC => null()

    ! Pointer to the entries of the C-matrix
    real(DP), dimension(:), pointer :: p_DC => null()

    ! Spatial discretisation structure for velocity
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscrV => null()

    ! Spatial discretisation structure for pressure
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscrP => null()

    ! Multiplication factors for the submatrices; taken from the system matrix.
    real(DP), dimension(3,3) :: Dmultipliers

  end type

!</typeblock>

!</types>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_init_NavSt2D (rvanka, rmatrix, csubtype)

!<description>
  ! Initialises the Vanka variant for 2D Navier-Stokes problems.
!</description>

!<input>
  ! The system matrix of the linear system.
  type(t_matrixBlock), intent(in), target :: rmatrix

  ! Desired subtype. One of the VANKATP_NAVST2D_XXXX constants.
  integer, intent(in) :: csubtype
!</input>

!<output>
  ! The data structure that saves algorithm-specific parameters.
  type(t_vanka_NavSt2D), intent(out) :: rvanka
!</output>

!</subroutine>

  integer :: i
  integer(I32) :: celemV, celemP
  type(t_blockDiscretisation), pointer :: p_rblockDiscr
  type(t_elementDistribution), pointer :: p_relementDistrV, p_relementDistrP

    ! Matrix must be 3x3.
    if ((rmatrix%nblocksPerCol .ne. 3) .or. (rmatrix%nblocksPerRow .ne. 3)) then
      call output_line ("System matrix is not 3x3.",&
          OU_CLASS_ERROR,OU_MODE_STD,"vanka_init_NavSt2D")
      call sys_halt()
    end if

    ! Todo: Do more checks

    ! The structure of A(1,3) must be identical to A(3,1) and
    ! that of A(2,3) must be identical to A(3,2).
    if ((rmatrix%RmatrixBlock(1,3)%NA .ne. rmatrix%RmatrixBlock(2,3)%NA) .or. &
        (rmatrix%RmatrixBlock(1,3)%NEQ .ne. rmatrix%RmatrixBlock(2,3)%NEQ)) then
      call output_line ("Structure of B1 and B2 different!",&
          OU_CLASS_ERROR,OU_MODE_STD,"vanka_init_NavSt2D")
      call sys_halt()
    end if

    if ((rmatrix%RmatrixBlock(3,1)%NA .ne. rmatrix%RmatrixBlock(3,2)%NA) .or. &
        (rmatrix%RmatrixBlock(3,1)%NEQ .ne. rmatrix%RmatrixBlock(3,2)%NEQ)) then
      call output_line ("Structure of D1 and D2 different!",&
          OU_CLASS_ERROR,OU_MODE_STD,"vanka_init_NavSt2D")
      call sys_halt()
    end if

    ! Virtually transposed matrices not supported
    if ((iand(rmatrix%RmatrixBlock(3,1)%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0) .or. &
        (iand(rmatrix%RmatrixBlock(3,1)%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0)) then
      call output_line ("Virtually transposed matrices not suppored!",&
          OU_CLASS_ERROR,OU_MODE_STD,"vanka_init_NavSt2D")
      call sys_halt()
    end if

    ! Sorted matrices not supported.
    if ((rmatrix%RmatrixBlock(1,1)%bcolumnsSorted .or. rmatrix%RmatrixBlock(1,1)%browsSorted) .or. &
        (rmatrix%RmatrixBlock(2,2)%bcolumnsSorted .or. rmatrix%RmatrixBlock(2,2)%browsSorted)) then
      call output_line ("Sorted matrices not supported!",&
          OU_CLASS_ERROR,OU_MODE_STD,"vanka_init_NavSt2D")
      call sys_halt()
    end if

    ! Store the subtype
    rvanka%csubtype = csubtype

    ! Fill the output structure with data of the matrices.
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(1,1),&
        rvanka%p_DA11)
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(2,2),&
        rvanka%p_DA22)
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(1,3),&
        rvanka%p_DB1)
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(2,3),&
        rvanka%p_DB2)
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(3,1),&
        rvanka%p_DD1)
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(3,2),&
        rvanka%p_DD2)
    call lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,1),&
        rvanka%p_KcolA)
    call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,1), &
        rvanka%p_KldA )
    call lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,3),&
        rvanka%p_KcolB)
    call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,3), &
        rvanka%p_KldB)
    call lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(3,1),&
        rvanka%p_KcolD)
    call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(3,1), &
        rvanka%p_KldD)
    if (rmatrix%RmatrixBlock(1,1)%cmatrixFormat .eq. LSYSSC_MATRIX9) then
      call lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(1,1), &
                              rvanka%p_KdiagonalA)
    else
      rvanka%p_KdiagonalA => rvanka%p_KldA
    end if

    ! Are the A12/A21 matrices present?
    if (lsysbl_isSubmatrixPresent(rmatrix,1,2)) then

      call lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,2),&
          rvanka%p_KcolA12)
      call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,2), &
          rvanka%p_KldA12)

      call lsyssc_getbase_double(rmatrix%RmatrixBlock(1,2),&
          rvanka%p_DA12)
      call lsyssc_getbase_double(rmatrix%RmatrixBlock(2,1),&
          rvanka%p_DA21)

    end if

    ! Is the C-Matrix present?
    if (lsysbl_isSubmatrixPresent(rmatrix,3,3)) then

      call lsyssc_getbase_double(rmatrix%RmatrixBlock(3,3),&
          rvanka%p_DC)

      call lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(3,3),&
          rvanka%p_KcolC)
      call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(3,3), &
          rvanka%p_KldC)

      if (rmatrix%RmatrixBlock(3,3)%cmatrixFormat .eq. LSYSSC_MATRIX9) then
        call lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(3,3), &
                                rvanka%p_KdiagonalC)
      else
        rvanka%p_KdiagonalC => rvanka%p_KldC
      end if

    end if

    ! Get the multiplication factors of the submatrices.
    rvanka%Dmultipliers(1:3,1:3) = &
        rmatrix%RmatrixBlock(1:3,1:3)%dscaleFactor

    ! Get the block discretisation structure from the matrix.
    p_rblockDiscr => rmatrix%p_rblockDiscrTest

    if (.not. associated(p_rblockDiscr)) then
      call output_line ("No block discretisation assigned to matrix!",&
          OU_CLASS_ERROR,OU_MODE_STD,"vanka_init_NavSt2D")
      call sys_halt()
    end if

    ! Get the discretisation structure of V and P from the block
    ! discretisation structure.
    rvanka%p_rspatialDiscrV => p_rblockDiscr%RspatialDiscr(1)
    rvanka%p_rspatialDiscrP => p_rblockDiscr%RspatialDiscr(3)

    if (p_rblockDiscr%RspatialDiscr(1)%inumFESpaces .ne. &
        p_rblockDiscr%RspatialDiscr(2)%inumFESpaces) then
      call output_line (&
          "Discretisation structures of X- and Y-velocity incompatible!",&
          OU_CLASS_ERROR,OU_MODE_STD,"vanka_init_NavSt2D")
      call sys_halt()
    end if

    if ((rvanka%p_rspatialDiscrP%inumFESpaces .ne. 1) .and. &
        (rvanka%p_rspatialDiscrP%inumFESpaces .ne. &
          rvanka%p_rspatialDiscrV%inumFESpaces)) then
      ! Either there must be only one element type for the pressure, or there one
      ! pressure element distribution for every velocity element distribution!
      ! If this is not the case, we cannot determine (at least not in reasonable time)
      ! which element type the pressure represents on a cell!
      call output_line (&
          "Discretisation structures of velocity and pressure incompatible!",&
          OU_CLASS_ERROR,OU_MODE_STD,"vanka_init_NavSt2D")
      call sys_halt()
    end if


    ! Loop through all FE spaces
    do i = 1, rvanka%p_rspatialDiscrV%inumFESpaces

      ! Get the corresponding element distributions of V.
      p_relementDistrV => &
          rvanka%p_rspatialDiscrV%RelementDistr(i)

      ! Either the same element for P everywhere, or there must be given one
      ! element distribution in the pressure for every velocity element distribution.
      if (rvanka%p_rspatialDiscrP%inumFESpaces .gt. 1) then
        p_relementDistrP => &
            rvanka%p_rspatialDiscrP%RelementDistr(i)
      else
        p_relementDistrP => &
            rvanka%p_rspatialDiscrP%RelementDistr(1)
      end if

      ! Get the elements
      celemV = elem_getPrimaryElement(p_relementDistrV%celement)
      celemP = elem_getPrimaryElement(p_relementDistrP%celement)

      ! Check whether we support the discretisation.
      if((celemV .eq. EL_Q1T) .and. (celemP .eq. EL_Q0)) cycle
      if((celemV .eq. EL_Q2) .and. (celemP .eq. EL_QP1)) cycle

      ! $TODO$: Check whether the element combination matches your new Vanka
      !         implementation and, if so, cycle the loop to avoid that the
      !         error message below is printed and the program is aborted.

      ! If we come out here, the discretisation is not supported
      call output_line ("Discretisation not supported!",&
          OU_CLASS_ERROR,OU_MODE_STD,"vanka_init_NavSt2D")
      call sys_halt()

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_solve_NavSt2D (rvanka, rsol, rrhs, niterations, domega)

!<description>
  ! Performs a desired number of iterations of the Vanka solver for
  ! 2D Navier-Stokes problems.
!</description>

!<input>
  ! The Vanka structure that saves algorithm-specific parameters.
  type(t_vanka_NavSt2D), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! The number of iterations that are to be performed
  integer, intent(in) :: niterations

  ! The relaxation parameter in range (0,2)
  real(DP), intent(in) :: domega
!</input>

!<inputoutput>
  ! The iteration vector that is to be updated.
  type(t_vectorBlock), intent(inout) :: rsol
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: ielementdist
  integer(I32) :: celemV, celemP
  integer, dimension(:), pointer :: p_IelementList
  type(t_elementDistribution), pointer :: p_relementDistrV
  type(t_elementDistribution), pointer :: p_relementDistrP

    ! Nothing to do?
    if(niterations .le. 0) return

    ! Loop through the element distributions of the velocity.
    do ielementdist = 1,rvanka%p_rspatialDiscrV%inumFESpaces

      ! Get the corresponding element distributions of u and p.
      p_relementDistrV => &
          rvanka%p_rspatialDiscrV%RelementDistr(ielementdist)

      ! Either the same element for P everywhere, or there must be given one
      ! element distribution in the pressure for every velocity element distribution.
      if (rvanka%p_rspatialDiscrP%inumFESpaces .gt. 1) then
        p_relementDistrP => &
            rvanka%p_rspatialDiscrP%RelementDistr(ielementdist)
      else
        p_relementDistrP => &
            rvanka%p_rspatialDiscrP%RelementDistr(1)
      end if

      ! Get the list of the elements to process.
      ! We take the element list of the X-velocity as "primary" element list
      ! and assume that it coincides to that of the Y-velocity (and to that
      ! of the pressure).
      call storage_getbase_int (p_relementDistrV%h_IelementList,p_IelementList)

      ! Which element combination do we have now?
      celemV = elem_getPrimaryElement(p_relementDistrV%celement)
      celemP = elem_getPrimaryElement(p_relementDistrP%celement)
      if ((celemV .eq. EL_Q1T) .and. (celemP .eq. EL_Q0)) then

        ! Q1~/Q0 discretisation

        ! Which VANKA subtype do we have? The diagonal VANKA or the full VANKA?
        select case (rvanka%csubtype)
        case (VANKATP_NAVST2D_DIAG)
          ! Call the jacobi-style vanka
          call vanka_NS2D_Q1TQ0_js(rvanka, rsol, rrhs, niterations, &
                                   domega, p_IelementList)

        case (VANKATP_NAVST2D_FULL)
          if (.not. associated(rvanka%p_DA12)) then
            ! Call the block-diagonal vanka
            call vanka_NS2D_Q1TQ0_bd(rvanka, rsol, rrhs, niterations, &
                                     domega, p_IelementList)
          else
            ! Call the fully coupled vanka
            call vanka_NS2D_Q1TQ0_fc(rvanka, rsol, rrhs, niterations, &
                                     domega, p_IelementList)
          end if

        case default
          call output_line ("Unknown Vanka subtype!",&
              OU_CLASS_ERROR,OU_MODE_STD,"vanka_NavierStokes2D")
          call sys_halt()

        end select

      else if ((celemV .eq. EL_Q2) .and. (celemP .eq. EL_QP1)) then

        ! Q2/QP1 discretisation

        ! Which VANKA subtype do we have? The diagonal VANKA or the full VANKA?
        select case (rvanka%csubtype)
        case (VANKATP_NAVST2D_DIAG)
          ! Call the jacobi-style vanka
          call vanka_NS2D_Q2QP1_js(rvanka, rsol, rrhs, niterations, &
                                   domega, p_IelementList)
        
        case (VANKATP_NAVST2D_FULL)
          ! Call the fully coupled vanka
          call vanka_NS2D_Q2QP1_fc(rvanka, rsol, rrhs, niterations, &
                                    domega, p_IelementList)
        
        case default
          call output_line ("Unknown Vanka subtype!",&
              OU_CLASS_ERROR,OU_MODE_STD,"vanka_NavierStokes2D")
          call sys_halt()

        end select

      ! $TODO$: Add an else-if case for your new Vanka implementation, and
      !         call the corresponding subroutine to perform the dirty work.

      end if

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_NS2D_Q1TQ0_js(rvanka, rsol, rrhs, niterations, domega, &
                                 IelementList)

!<description>
  ! Performs a desired number of iterations of the Vanka solver for
  ! 2D Navier-Stokes problem, "diagonal" variant for Q1~/Q0 discretisations.
!</description>

!<input>
  ! t_vanka_NavSt2D structure that saves algorithm-specific parameters.
  type(t_vanka_NavSt2D), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! The number of iterations that are to be performed
  integer, intent(in) :: niterations

  ! Relaxation parameter.
  real(DP), intent(in) :: domega

  ! A list of element numbers where Vanka should be applied to.
  integer, dimension(:), intent(in) :: IelementList
!</input>

!<inputoutput>
  ! The iteration vector that is to be updated.
  type(t_vectorBlock), intent(inout) :: rsol
!</inputoutput>

!</subroutine>

  ! How many local DOFs do we have?
  integer, parameter :: ndofV = 4   ! Dofs per velocity
  integer, parameter :: ndofP = 1   ! Dofs per pressure
  integer, parameter :: ndof = 2*ndofV+ndofP

  ! Triangulation information
  type(t_triangulation), pointer :: p_rtria
  integer, dimension(:,:), pointer :: p_IedgesAtElement

  ! DOF-mapping arrays
  integer, dimension(ndofV) :: IdofV
  integer, dimension(ndofP) :: IdofP

  ! Variables for the local system
  real(DP), dimension(ndofV) :: Da1, Da2, Du1, Du2, Df1, Df2
  real(DP), dimension(ndofV,ndofP) :: Db1, Db2
  real(DP), dimension(ndofP,ndofV) :: Dd1, Dd2
  real(DP), dimension(ndofP) :: Ds, Dc, Dup, Dfp

  ! Multiplication factors
  real(DP), dimension(3,3) :: Dmult

  ! Quick access for the matrix arrays
  integer, dimension(:), pointer :: p_KldA,p_KldA12,p_KldB,p_KldC,p_KldD,&
      p_KcolA,p_KcolA12,p_KcolB,p_KcolC,p_KcolD,p_KdiagA,p_KdiagC
  real(DP), dimension(:), pointer :: p_DA11,p_DA12,p_DA21,p_DA22,p_DB1,p_DB2,&
      p_DD1,p_DD2,p_DC

  ! Data arrays of vectors
  real(DP), dimension(:), pointer :: p_DrhsU,p_DrhsV,p_DrhsP,&
                                     p_DvecU,p_DvecV,p_DvecP

  ! local variables
  integer :: i,j,iel,ielidx,i1,i2,k,l,iter
  real(DP) :: daux,daux1,daux2
  logical :: bHaveA12, bHaveC

    ! Get the arrays from the triangulation
    p_rtria => rvanka%p_rspatialDiscrV%p_rtriangulation
    call storage_getbase_int2d (p_rtria%h_IedgesAtElement, p_IedgesAtElement)

    ! Get the pointers to the vector data
    call lsyssc_getbase_double(rsol%RvectorBlock(1), p_DvecU)
    call lsyssc_getbase_double(rsol%RvectorBlock(2), p_DvecV)
    call lsyssc_getbase_double(rsol%RvectorBlock(3), p_DvecP)
    call lsyssc_getbase_double(rrhs%RvectorBlock(1), p_DrhsU)
    call lsyssc_getbase_double(rrhs%RvectorBlock(2), p_DrhsV)
    call lsyssc_getbase_double(rrhs%RvectorBlock(3), p_DrhsP)

    ! Let us assume we do not have the optional matrices
    bHaveA12 = .false.
    bHaveC = .false.

    ! Get the pointers from the vanka structure
    p_KldA => rvanka%p_KldA
    p_KcolA => rvanka%p_KcolA
    p_KdiagA => rvanka%p_KdiagonalA
    p_DA11 => rvanka%p_DA11
    p_DA22 => rvanka%p_DA22
    p_KldB => rvanka%p_KldB
    p_KcolB => rvanka%p_KcolB
    p_DB1 => rvanka%p_DB1
    p_DB2 => rvanka%p_DB2
    p_KldD => rvanka%p_KldD
    p_KcolD => rvanka%p_KcolD
    p_DD1 => rvanka%p_DD1
    p_DD2 => rvanka%p_DD2

    if(associated(rvanka%p_DA12)) then
      bHaveA12 = .true.
      p_KldA12 => rvanka%p_KldA12
      p_KcolA12 => rvanka%p_KcolA12
      p_DA12 => rvanka%p_DA12
      p_DA21 => rvanka%p_DA21
    end if

    if(associated(rvanka%p_DC)) then
      bHaveC = .true.
      p_KldC => rvanka%p_KldC
      p_KcolC => rvanka%p_KcolC
      p_KdiagC => rvanka%p_KdiagonalC
      p_DC => rvanka%p_DC
    end if

    ! Get the multiplication factors
    Dmult = rvanka%Dmultipliers

    ! Take care of the "soft-deactivation" of the sub-matrices
    bHaveA12 = bHaveA12 .and. ((Dmult(1,2) .ne. 0.0_DP) .or. &
                               (Dmult(2,1) .ne. 0.0_DP))
    bHaveC = bHaveC .and. (Dmult(3,3) .ne. 0.0_DP)

    ! Clear the optional matrices
    Dc = 0.0_DP

    ! Now which of the optional matrices are present?
    if((.not. bHaveA12) .and. (.not. bHaveC)) then

      do iter = 1, niterations
        ! No optional matrices
        do ielidx = 1, size(IelementList)

          ! Get the element number which is to be processed.
          iel = IelementList(ielidx)

          ! Get all DOFs for this element
          IdofV(:) = p_IedgesAtElement(:,iel)
          IdofP(1) = iel

          ! First of all, fetch the local RHS
          do i = 1, ndofV
            Df1(i) = p_DrhsU(IdofV(i))   ! f_u
            Df2(i) = p_DrhsV(IdofV(i))   ! f_v
          end do
          do i = 1, ndofP
            Dfp(i) = p_DrhsP(IdofP(i))    ! f_p
          end do

          ! Let us update the local RHS vector by subtracting A*u from it:
          ! f_u := f_u - A11*u
          ! f_v := f_v - A22*v
          do k = 1, ndofV
            i1 = p_KldA(IdofV(k))
            i2 = p_KldA(IdofV(k)+1)-1
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = i1, i2
              j = p_KcolA(i)
              daux1 = daux1 + p_DA11(i)*p_DvecU(j)
              daux2 = daux2 + p_DA22(i)*p_DvecV(j)
            end do
            Df1(k) = Df1(k) - Dmult(1,1)*daux1
            Df2(k) = Df2(k) - Dmult(2,2)*daux2
            ! Get the main diagonal entries
            j = p_KdiagA(IdofV(k))
            Da1(k) = Dmult(1,1)*p_DA11(j)
            Da2(k) = Dmult(2,2)*p_DA22(j)
          end do

          ! Now we also need to subtract B*p from our RHS, and by the same time,
          ! we will build the local B matrices.
          ! f_u := f_u - B1*p
          ! f_v := f_v - B2*p
          do k = 1, ndofV
            i1 = p_KldB(IdofV(k))
            i2 = p_KldB(IdofV(k)+1)-1
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = i1, i2
              j = p_KcolB(i)
              daux = p_DvecP(j)
              daux1 = daux1 + p_DB1(i)*daux
              daux2 = daux2 + p_DB2(i)*daux
              do l = 1, ndofP
                if(j .eq. IdofP(l)) then
                  Db1(k,l) = Dmult(1,3)*p_DB1(i)
                  Db2(k,l) = Dmult(2,3)*p_DB2(i)
                  exit
                end if
              end do
            end do
            Df1(k) = Df1(k) - Dmult(1,3)*daux1
            Df2(k) = Df2(k) - Dmult(2,3)*daux2
          end do

          ! Now we also need to subtract D*u from our RHS, and by the same time,
          ! we will build the local D matrices.
          ! f_p := f_p - D1*u - D2*v
          do k = 1, ndofP
            i1 = p_KldD(IdofP(k))
            i2 = p_KldD(IdofP(k)+1)-1
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = i1, i2
              j = p_KcolD(i)
              daux1 = daux1 + p_DD1(i)*p_DvecU(j)
              daux2 = daux2 + p_DD2(i)*p_DvecV(j)
              do l = 1, ndofV
                if(j .eq. IdofV(l)) then
                  Dd1(k,l) = Dmult(3,1)*p_DD1(i)
                  Dd2(k,l) = Dmult(3,2)*p_DD2(i)
                  exit
                end if
              end do
            end do
            Dfp(k) = Dfp(k) - Dmult(3,1)*daux1 - Dmult(3,2)*daux2
          end do

          ! Invert A1 and A2
          do i = 1, ndofV
            Da1(i) = 1.0_DP / Da1(i)
            Da2(i) = 1.0_DP / Da2(i)
          end do

          ! Precalculate D * A^-1
          do i = 1, ndofV
            do j = 1, ndofP
              Dd1(j,i) = Dd1(j,i)*Da1(i)
              Dd2(j,i) = Dd2(j,i)*Da2(i)
            end do
          end do

          ! Calculate Schur-Complement of A
          ! S := -C + D * A^-1 * B
          do j = 1, ndofP
            Ds(j) = -Dc(j)
            do i = 1, ndofV
              Ds(j) = Ds(j) + Dd1(j,i)*Db1(i,j) &
                            + Dd2(j,i)*Db2(i,j)
            end do
          end do

          ! Calculate pressure
          ! p := S^-1 * (D * A^-1 * f_u - f_p)
          do j = 1, ndofP
            daux = -Dfp(j)
            do i = 1, ndofV
              daux = daux + Dd1(j,i)*Df1(i) &
                          + Dd2(j,i)*Df2(i)
            end do
            Dup(j) = daux / Ds(j)
          end do

          ! Calculate X- and Y-velocity
          do i = 1, ndofV
            do j = 1, ndofP
              Df1(i) = Df1(i) - Db1(i,j)*Dup(j)
              Df2(i) = Df2(i) - Db2(i,j)*Dup(j)
            end do
            Du1(i) = Da1(i)*Df1(i)
            Du2(i) = Da2(i)*Df2(i)
          end do

          ! Incorporate our local solution into the global one.
          do i = 1, ndofV
            j = IdofV(i)
            p_DvecU(j) = p_DvecU(j) + domega * Du1(i)
            p_DvecV(j) = p_DvecV(j) + domega * Du2(i)
          end do
          do i = 1, ndofP
            j = IdofP(i)
            p_DvecP(j) = p_DvecP(j) + domega * Dup(i)
          end do

        end do ! ielidx

      end do ! iter

    else

      do iter = 1, niterations

        ! General case
        do ielidx=1, size(IelementList)

          ! Get the element number which is to be processed.
          iel = IelementList(ielidx)

          ! Get all DOFs for this element
          IdofV(:) = p_IedgesAtElement(:,iel)
          IdofP(1) = iel

          ! First of all, fetch the local RHS
          do i = 1, ndofV
            Df1(i) = p_DrhsU(IdofV(i))   ! f_u
            Df2(i) = p_DrhsV(IdofV(i))   ! f_v
          end do
          do i = 1, ndofP
            Dfp(i) = p_DrhsP(IdofP(i))    ! f_p
          end do

          ! Let us update the local RHS vector by subtracting A*u from it:
          ! f_u := f_u - A11*u
          ! f_v := f_v - A22*v
          do k = 1, ndofV
            i1 = p_KldA(IdofV(k))
            i2 = p_KldA(IdofV(k)+1)-1
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = i1, i2
              j = p_KcolA(i)
              daux1 = daux1 + p_DA11(i)*p_DvecU(j)
              daux2 = daux2 + p_DA22(i)*p_DvecV(j)
            end do
            Df1(k) = Df1(k) - Dmult(1,1)*daux1
            Df2(k) = Df2(k) - Dmult(2,2)*daux2
            ! Get the main diagonal entries
            j = p_KdiagA(IdofV(k))
            Da1(k) = Dmult(1,1)*p_DA11(j)
            Da2(k) = Dmult(2,2)*p_DA22(j)
          end do

          ! What about A12/A21?
          if(bHaveA12) then
            ! f_u := f_u - A12*v
            ! f_v := f_v - A21*u
            do k = 1, ndofV
              i1 = p_KldA12(IdofV(k))
              i2 = p_KldA12(IdofV(k)+1)-1
              daux1 = 0.0_DP
              daux2 = 0.0_DP
              do i = i1, i2
                j = p_KcolA12(i)
                daux1 = daux1 + p_DA12(i)*p_DvecV(j)
                daux2 = daux2 + p_DA21(i)*p_DvecU(j)
              end do
              Df1(k) = Df1(k) - Dmult(1,2)*daux1
              Df2(k) = Df2(k) - Dmult(2,1)*daux2
            end do
          end if

          ! Now we also need to subtract B*p from our RHS, and by the same time,
          ! we will build the local B matrices.
          ! f_u := f_u - B1*p
          ! f_v := f_v - B2*p
          do k = 1, ndofV
            i1 = p_KldB(IdofV(k))
            i2 = p_KldB(IdofV(k)+1)-1
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = i1, i2
              j = p_KcolB(i)
              daux = p_DvecP(j)
              daux1 = daux1 + p_DB1(i)*daux
              daux2 = daux2 + p_DB2(i)*daux
              do l = 1, ndofP
                if(j .eq. IdofP(l)) then
                  Db1(k,l) = Dmult(1,3)*p_DB1(i)
                  Db2(k,l) = Dmult(2,3)*p_DB2(i)
                  exit
                end if
              end do
            end do
            Df1(k) = Df1(k) - Dmult(1,3)*daux1
            Df2(k) = Df2(k) - Dmult(2,3)*daux2
          end do

          ! Now we also need to subtract D*u from our RHS, and by the same time,
          ! we will build the local D matrices.
          ! f_p := f_p - D1*u - D2*v
          do k = 1, ndofP
            i1 = p_KldD(IdofP(k))
            i2 = p_KldD(IdofP(k)+1)-1
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = i1, i2
              j = p_KcolD(i)
              daux1 = daux1 + p_DD1(i)*p_DvecU(j)
              daux2 = daux2 + p_DD2(i)*p_DvecV(j)
              do l = 1, ndofV
                if(j .eq. IdofV(l)) then
                  Dd1(k,l) = Dmult(3,1)*p_DD1(i)
                  Dd2(k,l) = Dmult(3,2)*p_DD2(i)
                  exit
                end if
              end do
            end do
            Dfp(k) = Dfp(k) - Dmult(3,1)*daux1 - Dmult(3,2)*daux2
          end do

          ! Do we have a C-matrix?
          if(bHaveC) then
            ! Yes, so update the local RHS
            ! f_p := f_p - C*p
            do k = 1, ndofP
              i1 = p_KldC(IdofP(k))
              i2 = p_KldC(IdofP(k)+1)-1
              daux1 = 0.0_DP
              do i = i1, i2
                daux1 = daux1 + p_DC(i)*p_DvecP(p_KcolC(i))
              end do
              Dfp(k) = Dfp(k) - Dmult(3,3)*daux1
              ! Get the main diagonal entry of C
              Dc(k) = Dmult(3,3)*p_DC(p_KdiagC(IdofP(k)))
            end do
          end if

          ! Invert A1 and A2
          do i = 1, ndofV
            Da1(i) = 1.0_DP / Da1(i)
            Da2(i) = 1.0_DP / Da2(i)
          end do

          ! Precalculate D * A^-1
          do i = 1, ndofV
            do j = 1, ndofP
              Dd1(j,i) = Dd1(j,i)*Da1(i)
              Dd2(j,i) = Dd2(j,i)*Da2(i)
            end do
          end do

          ! Calculate Schur-Complement of A
          ! S := -C + D * A^-1 * B
          do j = 1, ndofP
            Ds(j) = -Dc(j)
            do i = 1, ndofV
              Ds(j) = Ds(j) + Dd1(j,i)*Db1(i,j) &
                            + Dd2(j,i)*Db2(i,j)
            end do
          end do

          ! Calculate pressure
          ! p := S^-1 * (D * A^-1 * f_u - f_p)
          do j = 1, ndofP
            daux = -Dfp(j)
            do i = 1, ndofV
              daux = daux + Dd1(j,i)*Df1(i) &
                          + Dd2(j,i)*Df2(i)
            end do
            Dup(j) = daux / Ds(j)
          end do

          ! Calculate X- and Y-velocity
          do i = 1, ndofV
            do j = 1, ndofP
              Df1(i) = Df1(i) - Db1(i,j)*Dup(j)
              Df2(i) = Df2(i) - Db2(i,j)*Dup(j)
            end do
            Du1(i) = Da1(i)*Df1(i)
            Du2(i) = Da2(i)*Df2(i)
          end do

          ! Incorporate our local solution into the global one.
          do i = 1, ndofV
            j = IdofV(i)
            p_DvecU(j) = p_DvecU(j) + domega * Du1(i)
            p_DvecV(j) = p_DvecV(j) + domega * Du2(i)
          end do
          do i = 1, ndofP
            j = IdofP(i)
            p_DvecP(j) = p_DvecP(j) + domega * Dup(i)
          end do

        end do ! ielidx

      end do ! iter

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_NS2D_Q2QP1_js(rvanka, rsol, rrhs, niterations, domega, &
                                 IelementList)

!<description>
  ! Performs a desired number of iterations of the Vanka solver for
  ! 2D Navier-Stokes problem, "diagonal" variant for Q2/QP1 discretisations.
!</description>

!<input>
  ! t_vanka_NavSt2D structure that saves algorithm-specific parameters.
  type(t_vanka_NavSt2D), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! The number of iterations that are to be performed
  integer, intent(in) :: niterations

  ! Relaxation parameter.
  real(DP), intent(in) :: domega

  ! A list of element numbers where Vanka should be applied to.
  integer, dimension(:), intent(in) :: IelementList
!</input>

!<inputoutput>
  ! The iteration vector that is to be updated.
  type(t_vectorBlock), intent(inout) :: rsol
!</inputoutput>

!</subroutine>

  ! How many local DOFs do we have?
  integer, parameter :: ndofV = 9   ! Dofs per velocity
  integer, parameter :: ndofP = 3   ! Dofs per pressure
  integer, parameter :: ndof = 2*ndofV+ndofP

  ! Triangulation information
  type(t_triangulation), pointer :: p_rtria
  integer, dimension(:,:), pointer :: p_IedgesAtElement
  integer, dimension(:,:), pointer :: p_IverticesAtElement

  ! DOF-mapping arrays
  integer, dimension(ndofV) :: IdofV
  integer, dimension(ndofP) :: IdofP

  ! Variables for the local system
  real(DP), dimension(ndofV) :: Da1, Da2, Du1, Du2, Df1, Df2
  real(DP), dimension(ndofV,ndofP) :: Db1, Db2
  real(DP), dimension(ndofP,ndofV) :: Dd1, Dd2
  real(DP), dimension(ndofP) :: Ds, Dc, Dup, Dfp

  ! Multiplication factors
  real(DP), dimension(3,3) :: Dmult

  ! Quick access for the matrix arrays
  integer, dimension(:), pointer :: p_KldA,p_KldA12,p_KldB,p_KldC,p_KldD,&
      p_KcolA,p_KcolA12,p_KcolB,p_KcolC,p_KcolD,p_KdiagA,p_KdiagC
  real(DP), dimension(:), pointer :: p_DA11,p_DA12,p_DA21,p_DA22,p_DB1,p_DB2,&
      p_DD1,p_DD2,p_DC

  ! Data arrays of vectors
  real(DP), dimension(:), pointer :: p_DrhsU,p_DrhsV,p_DrhsP,&
                                     p_DvecU,p_DvecV,p_DvecP

  ! local variables
  integer :: i,j,iel,ielidx,i1,i2,k,l,iter,nvt,nmt,nel
  real(DP) :: daux,daux1,daux2
  logical :: bHaveA12, bHaveC

    ! Get the arrays from the triangulation
    p_rtria => rvanka%p_rspatialDiscrV%p_rtriangulation
    call storage_getbase_int2d (p_rtria%h_IedgesAtElement, p_IedgesAtElement)
    call storage_getbase_int2d (p_rtria%h_IverticesAtElement, p_IverticesAtElement)
    nvt = p_rtria%nvt
    nmt = p_rtria%nmt
    nel = p_rtria%nel

    ! Get the pointers to the vector data
    call lsyssc_getbase_double(rsol%RvectorBlock(1), p_DvecU)
    call lsyssc_getbase_double(rsol%RvectorBlock(2), p_DvecV)
    call lsyssc_getbase_double(rsol%RvectorBlock(3), p_DvecP)
    call lsyssc_getbase_double(rrhs%RvectorBlock(1), p_DrhsU)
    call lsyssc_getbase_double(rrhs%RvectorBlock(2), p_DrhsV)
    call lsyssc_getbase_double(rrhs%RvectorBlock(3), p_DrhsP)

    ! Let us assume we do not have the optional matrices
    bHaveA12 = .false.
    bHaveC = .false.

    ! Get the pointers from the vanka structure
    p_KldA => rvanka%p_KldA
    p_KcolA => rvanka%p_KcolA
    p_KdiagA => rvanka%p_KdiagonalA
    p_DA11 => rvanka%p_DA11
    p_DA22 => rvanka%p_DA22
    p_KldB => rvanka%p_KldB
    p_KcolB => rvanka%p_KcolB
    p_DB1 => rvanka%p_DB1
    p_DB2 => rvanka%p_DB2
    p_KldD => rvanka%p_KldD
    p_KcolD => rvanka%p_KcolD
    p_DD1 => rvanka%p_DD1
    p_DD2 => rvanka%p_DD2

    if(associated(rvanka%p_DA12)) then
      bHaveA12 = .true.
      p_KldA12 => rvanka%p_KldA12
      p_KcolA12 => rvanka%p_KcolA12
      p_DA12 => rvanka%p_DA12
      p_DA21 => rvanka%p_DA21
    end if

    if(associated(rvanka%p_DC)) then
      bHaveC = .true.
      p_KldC => rvanka%p_KldC
      p_KcolC => rvanka%p_KcolC
      p_KdiagC => rvanka%p_KdiagonalC
      p_DC => rvanka%p_DC
    end if

    ! Get the multiplication factors
    Dmult = rvanka%Dmultipliers

    ! Take care of the "soft-deactivation" of the sub-matrices
    bHaveA12 = bHaveA12 .and. ((Dmult(1,2) .ne. 0.0_DP) .or. &
                               (Dmult(2,1) .ne. 0.0_DP))
    bHaveC = bHaveC .and. (Dmult(3,3) .ne. 0.0_DP)

    ! Clear the optional matrices
    Dc = 0.0_DP

    ! Now which of the optional matrices are present?
    if((.not. bHaveA12) .and. (.not. bHaveC)) then

      do iter = 1, niterations
        ! No optional matrices
        do ielidx = 1, size(IelementList)

          ! Get the element number which is to be processed.
          iel = IelementList(ielidx)

          ! Get all DOFs for this element
          IdofV(1:4) = p_IverticesAtElement(1:4,iel)
          IdofV(5:8) = p_IedgesAtElement(1:4,iel)+nvt
          IdofV(9) = iel+nvt+nmt
          IdofP(1) = iel
          IdofP(2) = iel+nel
          IdofP(3) = iel+2*nel

          ! First of all, fetch the local RHS
          do i = 1, ndofV
            Df1(i) = p_DrhsU(IdofV(i))   ! f_u
            Df2(i) = p_DrhsV(IdofV(i))   ! f_v
          end do
          do i = 1, ndofP
            Dfp(i) = p_DrhsP(IdofP(i))    ! f_p
          end do

          ! Let us update the local RHS vector by subtracting A*u from it:
          ! f_u := f_u - A11*u
          ! f_v := f_v - A22*v
          do k = 1, ndofV
            i1 = p_KldA(IdofV(k))
            i2 = p_KldA(IdofV(k)+1)-1
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = i1, i2
              j = p_KcolA(i)
              daux1 = daux1 + p_DA11(i)*p_DvecU(j)
              daux2 = daux2 + p_DA22(i)*p_DvecV(j)
            end do
            Df1(k) = Df1(k) - Dmult(1,1)*daux1
            Df2(k) = Df2(k) - Dmult(2,2)*daux2
            ! Get the main diagonal entries
            j = p_KdiagA(IdofV(k))
            Da1(k) = Dmult(1,1)*p_DA11(j)
            Da2(k) = Dmult(2,2)*p_DA22(j)
          end do

          ! Now we also need to subtract B*p from our RHS, and by the same time,
          ! we will build the local B matrices.
          ! f_u := f_u - B1*p
          ! f_v := f_v - B2*p
          do k = 1, ndofV
            i1 = p_KldB(IdofV(k))
            i2 = p_KldB(IdofV(k)+1)-1
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = i1, i2
              j = p_KcolB(i)
              daux = p_DvecP(j)
              daux1 = daux1 + p_DB1(i)*daux
              daux2 = daux2 + p_DB2(i)*daux
              do l = 1, ndofP
                if(j .eq. IdofP(l)) then
                  Db1(k,l) = Dmult(1,3)*p_DB1(i)
                  Db2(k,l) = Dmult(2,3)*p_DB2(i)
                  exit
                end if
              end do
            end do
            Df1(k) = Df1(k) - Dmult(1,3)*daux1
            Df2(k) = Df2(k) - Dmult(2,3)*daux2
          end do

          ! Now we also need to subtract D*u from our RHS, and by the same time,
          ! we will build the local D matrices.
          ! f_p := f_p - D1*u - D2*v
          do k = 1, ndofP
            i1 = p_KldD(IdofP(k))
            i2 = p_KldD(IdofP(k)+1)-1
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = i1, i2
              j = p_KcolD(i)
              daux1 = daux1 + p_DD1(i)*p_DvecU(j)
              daux2 = daux2 + p_DD2(i)*p_DvecV(j)
              do l = 1, ndofV
                if(j .eq. IdofV(l)) then
                  Dd1(k,l) = Dmult(3,1)*p_DD1(i)
                  Dd2(k,l) = Dmult(3,2)*p_DD2(i)
                  exit
                end if
              end do
            end do
            Dfp(k) = Dfp(k) - Dmult(3,1)*daux1 - Dmult(3,2)*daux2
          end do

          ! Invert A1 and A2
          do i = 1, ndofV
            Da1(i) = 1.0_DP / Da1(i)
            Da2(i) = 1.0_DP / Da2(i)
          end do

          ! Precalculate D * A^-1
          do i = 1, ndofV
            do j = 1, ndofP
              Dd1(j,i) = Dd1(j,i)*Da1(i)
              Dd2(j,i) = Dd2(j,i)*Da2(i)
            end do
          end do

          ! Calculate Schur-Complement of A
          ! S := -C + D * A^-1 * B
          do j = 1, ndofP
            Ds(j) = -Dc(j)
            do i = 1, ndofV
              Ds(j) = Ds(j) + Dd1(j,i)*Db1(i,j) &
                            + Dd2(j,i)*Db2(i,j)
            end do
          end do

          ! Calculate pressure
          ! p := S^-1 * (D * A^-1 * f_u - f_p)
          do j = 1, ndofP
            daux = -Dfp(j)
            do i = 1, ndofV
              daux = daux + Dd1(j,i)*Df1(i) &
                          + Dd2(j,i)*Df2(i)
            end do
            Dup(j) = daux / Ds(j)
          end do

          ! Calculate X- and Y-velocity
          do i = 1, ndofV
            do j = 1, ndofP
              Df1(i) = Df1(i) - Db1(i,j)*Dup(j)
              Df2(i) = Df2(i) - Db2(i,j)*Dup(j)
            end do
            Du1(i) = Da1(i)*Df1(i)
            Du2(i) = Da2(i)*Df2(i)
          end do

          ! Incorporate our local solution into the global one.
          do i = 1, ndofV
            j = IdofV(i)
            p_DvecU(j) = p_DvecU(j) + domega * Du1(i)
            p_DvecV(j) = p_DvecV(j) + domega * Du2(i)
          end do
          do i = 1, ndofP
            j = IdofP(i)
            p_DvecP(j) = p_DvecP(j) + domega * Dup(i)
          end do

        end do ! ielidx

      end do ! iter

    else

      do iter = 1, niterations

        ! General case
        do ielidx=1, size(IelementList)

          ! Get the element number which is to be processed.
          iel = IelementList(ielidx)

          ! Get all DOFs for this element
          IdofV(1:4) = p_IverticesAtElement(1:4,iel)
          IdofV(5:8) = p_IedgesAtElement(1:4,iel)+nvt
          IdofV(9) = iel+nvt+nmt
          IdofP(1) = iel
          IdofP(2) = iel+nel
          IdofP(3) = iel+2*nel

          ! First of all, fetch the local RHS
          do i = 1, ndofV
            Df1(i) = p_DrhsU(IdofV(i))   ! f_u
            Df2(i) = p_DrhsV(IdofV(i))   ! f_v
          end do
          do i = 1, ndofP
            Dfp(i) = p_DrhsP(IdofP(i))    ! f_p
          end do

          ! Let us update the local RHS vector by subtracting A*u from it:
          ! f_u := f_u - A11*u
          ! f_v := f_v - A22*v
          do k = 1, ndofV
            i1 = p_KldA(IdofV(k))
            i2 = p_KldA(IdofV(k)+1)-1
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = i1, i2
              j = p_KcolA(i)
              daux1 = daux1 + p_DA11(i)*p_DvecU(j)
              daux2 = daux2 + p_DA22(i)*p_DvecV(j)
            end do
            Df1(k) = Df1(k) - Dmult(1,1)*daux1
            Df2(k) = Df2(k) - Dmult(2,2)*daux2
            ! Get the main diagonal entries
            j = p_KdiagA(IdofV(k))
            Da1(k) = Dmult(1,1)*p_DA11(j)
            Da2(k) = Dmult(2,2)*p_DA22(j)
          end do

          ! What about A12/A21?
          if(bHaveA12) then
            ! f_u := f_u - A12*v
            ! f_v := f_v - A21*u
            do k = 1, ndofV
              i1 = p_KldA12(IdofV(k))
              i2 = p_KldA12(IdofV(k)+1)-1
              daux1 = 0.0_DP
              daux2 = 0.0_DP
              do i = i1, i2
                j = p_KcolA12(i)
                daux1 = daux1 + p_DA12(i)*p_DvecV(j)
                daux2 = daux2 + p_DA21(i)*p_DvecU(j)
              end do
              Df1(k) = Df1(k) - Dmult(1,2)*daux1
              Df2(k) = Df2(k) - Dmult(2,1)*daux2
            end do
          end if

          ! Now we also need to subtract B*p from our RHS, and by the same time,
          ! we will build the local B matrices.
          ! f_u := f_u - B1*p
          ! f_v := f_v - B2*p
          do k = 1, ndofV
            i1 = p_KldB(IdofV(k))
            i2 = p_KldB(IdofV(k)+1)-1
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = i1, i2
              j = p_KcolB(i)
              daux = p_DvecP(j)
              daux1 = daux1 + p_DB1(i)*daux
              daux2 = daux2 + p_DB2(i)*daux
              do l = 1, ndofP
                if(j .eq. IdofP(l)) then
                  Db1(k,l) = Dmult(1,3)*p_DB1(i)
                  Db2(k,l) = Dmult(2,3)*p_DB2(i)
                  exit
                end if
              end do
            end do
            Df1(k) = Df1(k) - Dmult(1,3)*daux1
            Df2(k) = Df2(k) - Dmult(2,3)*daux2
          end do

          ! Now we also need to subtract D*u from our RHS, and by the same time,
          ! we will build the local D matrices.
          ! f_p := f_p - D1*u - D2*v
          do k = 1, ndofP
            i1 = p_KldD(IdofP(k))
            i2 = p_KldD(IdofP(k)+1)-1
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = i1, i2
              j = p_KcolD(i)
              daux1 = daux1 + p_DD1(i)*p_DvecU(j)
              daux2 = daux2 + p_DD2(i)*p_DvecV(j)
              do l = 1, ndofV
                if(j .eq. IdofV(l)) then
                  Dd1(k,l) = Dmult(3,1)*p_DD1(i)
                  Dd2(k,l) = Dmult(3,2)*p_DD2(i)
                  exit
                end if
              end do
            end do
            Dfp(k) = Dfp(k) - Dmult(3,1)*daux1 - Dmult(3,2)*daux2
          end do

          ! Do we have a C-matrix?
          if(bHaveC) then
            ! Yes, so update the local RHS
            ! f_p := f_p - C*p
            do k = 1, ndofP
              i1 = p_KldC(IdofP(k))
              i2 = p_KldC(IdofP(k)+1)-1
              daux1 = 0.0_DP
              do i = i1, i2
                daux1 = daux1 + p_DC(i)*p_DvecP(p_KcolC(i))
              end do
              Dfp(k) = Dfp(k) - Dmult(3,3)*daux1
              ! Get the main diagonal entry of C
              Dc(k) = Dmult(3,3)*p_DC(p_KdiagC(IdofP(k)))
            end do
          end if

          ! Invert A1 and A2
          do i = 1, ndofV
            Da1(i) = 1.0_DP / Da1(i)
            Da2(i) = 1.0_DP / Da2(i)
          end do

          ! Precalculate D * A^-1
          do i = 1, ndofV
            do j = 1, ndofP
              Dd1(j,i) = Dd1(j,i)*Da1(i)
              Dd2(j,i) = Dd2(j,i)*Da2(i)
            end do
          end do

          ! Calculate Schur-Complement of A
          ! S := -C + D * A^-1 * B
          do j = 1, ndofP
            Ds(j) = -Dc(j)
            do i = 1, ndofV
              Ds(j) = Ds(j) + Dd1(j,i)*Db1(i,j) &
                            + Dd2(j,i)*Db2(i,j)
            end do
          end do

          ! Calculate pressure
          ! p := S^-1 * (D * A^-1 * f_u - f_p)
          do j = 1, ndofP
            daux = -Dfp(j)
            do i = 1, ndofV
              daux = daux + Dd1(j,i)*Df1(i) &
                          + Dd2(j,i)*Df2(i)
            end do
            Dup(j) = daux / Ds(j)
          end do

          ! Calculate X- and Y-velocity
          do i = 1, ndofV
            do j = 1, ndofP
              Df1(i) = Df1(i) - Db1(i,j)*Dup(j)
              Df2(i) = Df2(i) - Db2(i,j)*Dup(j)
            end do
            Du1(i) = Da1(i)*Df1(i)
            Du2(i) = Da2(i)*Df2(i)
          end do

          ! Incorporate our local solution into the global one.
          do i = 1, ndofV
            j = IdofV(i)
            p_DvecU(j) = p_DvecU(j) + domega * Du1(i)
            p_DvecV(j) = p_DvecV(j) + domega * Du2(i)
          end do
          do i = 1, ndofP
            j = IdofP(i)
            p_DvecP(j) = p_DvecP(j) + domega * Dup(i)
          end do

        end do ! ielidx

      end do ! iter

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_NS2D_Q2QP1_fc(rvanka, rsol, rrhs, niterations, domega, &
                                 IelementList)

!<description>
  ! Performs a desired number of iterations of the Vanka solver for
  ! 2D Navier-Stokes problem, "full" variant for Q2/Qp1 discretisations.
!</description>

!<input>
  ! t_vanka_NavSt2D structure that saves algorithm-specific parameters.
  type(t_vanka_NavSt2D), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! The number of iterations that are to be performed
  integer, intent(in) :: niterations

  ! Relaxation parameter.
  real(DP), intent(in) :: domega

  ! A list of element numbers where Vanka should be applied to.
  integer, dimension(:), intent(in) :: IelementList
!</input>

!<inputoutput>
  ! The iteration vector that is to be updated.
  type(t_vectorBlock), intent(inout) :: rsol
!</inputoutput>

!</subroutine>

  ! How many local DOFs do we have?
  integer, parameter :: ndofV = 9   ! Dofs per velocity
  integer, parameter :: ndofP = 3   ! Dofs per pressure
  integer, parameter :: ndof = 2*ndofV+ndofP

  ! Triangulation information
  type(t_triangulation), pointer :: p_rtria
  integer, dimension(:,:), pointer :: p_IedgesAtElement
  integer, dimension(:,:), pointer :: p_IverticesAtElement
  integer :: nvt,nmt,nel

  ! Data arrays of vectors
  real(DP), dimension(:), pointer :: p_DrhsU,p_DrhsV,p_DrhsP,&
                                     p_DvecU,p_DvecV,p_DvecP

  ! DOF-mapping arrays
  integer, dimension(ndofV) :: IdofV
  integer, dimension(ndofP) :: IdofP

  ! Variables for the local system
  real(DP), dimension(ndof) :: Df
  real(DP), dimension(ndof,ndof) :: Da

  ! Multiplication factors
  real(DP), dimension(3,3) :: Dmult

  ! Quick access for the matrix arrays
  integer, dimension(:), pointer :: p_KldA,p_KldA12,p_KldB,p_KldC,p_KldD,&
      p_KcolA,p_KcolA12,p_KcolB,p_KcolC,p_KcolD
  real(DP), dimension(:), pointer :: p_DA11,p_DA12,p_DA21,p_DA22,p_DB1,p_DB2,&
      p_DD1,p_DD2,p_DC

  ! local variables
  integer :: i,j,iel,ielidx,i1,i2,k,l,o,p,iter
  real(DP) :: daux,daux1,daux2
  logical :: bHaveA12, bHaveC

  ! variables for LAPACK`s DGESV routine
  integer, dimension(ndof) :: Ipivot
  integer :: info


    ! Get the arrays from the triangulation
    p_rtria => rvanka%p_rspatialDiscrV%p_rtriangulation
    call storage_getbase_int2d (p_rtria%h_IedgesAtElement, p_IedgesAtElement)
    call storage_getbase_int2d (p_rtria%h_IverticesAtElement, p_IverticesAtElement)
    nvt = p_rtria%nvt
    nmt = p_rtria%nmt
    nel = p_rtria%nel


    ! Get the pointers to the vector data
    call lsyssc_getbase_double(rsol%RvectorBlock(1), p_DvecU)
    call lsyssc_getbase_double(rsol%RvectorBlock(2), p_DvecV)
    call lsyssc_getbase_double(rsol%RvectorBlock(3), p_DvecP)
    call lsyssc_getbase_double(rrhs%RvectorBlock(1), p_DrhsU)
    call lsyssc_getbase_double(rrhs%RvectorBlock(2), p_DrhsV)
    call lsyssc_getbase_double(rrhs%RvectorBlock(3), p_DrhsP)

    ! Let us assume we do not have the optional matrices
    bHaveA12 = .false.
    bHaveC = .false.

    ! Get the pointers from the vanka structure
    p_KldA => rvanka%p_KldA
    p_KcolA => rvanka%p_KcolA
    p_DA11 => rvanka%p_DA11
    p_DA22 => rvanka%p_DA22
    p_KldB => rvanka%p_KldB
    p_KcolB => rvanka%p_KcolB
    p_DB1 => rvanka%p_DB1
    p_DB2 => rvanka%p_DB2
    p_KldD => rvanka%p_KldD
    p_KcolD => rvanka%p_KcolD
    p_DD1 => rvanka%p_DD1
    p_DD2 => rvanka%p_DD2

    if(associated(rvanka%p_DA12)) then
      bHaveA12 = .true.
      p_KldA12 => rvanka%p_KldA12
      p_KcolA12 => rvanka%p_KcolA12
      p_DA12 => rvanka%p_DA12
      p_DA21 => rvanka%p_DA21
    end if

    if(associated(rvanka%p_DC)) then
      bHaveC = .true.
      p_KldC => rvanka%p_KldC
      p_KcolC => rvanka%p_KcolC
      p_DC => rvanka%p_DC
    end if

    ! Get the multiplication factors
    Dmult = rvanka%Dmultipliers

    ! Take care of the "soft-deactivation" of the sub-matrices
    bHaveC = bHaveC .and. (Dmult(3,3) .ne. 0.0_DP)
    bHaveA12 = bHaveA12 .and. &
      ((Dmult(1,2) .ne. 0.0_DP) .or. (Dmult(2,1) .ne. 0.0_DP))

    ! Perform the desired number of iterations
    do iter = 1, niterations

      ! Loop over all elements in the list
      do ielidx = 1, size(IelementList)

        ! Get the element number which is to be processed.
        iel = IelementList(ielidx)

        ! Get all DOFs for this element
        IdofV(1:4) = p_IverticesAtElement(1:4,iel)
        IdofV(5:8) = p_IedgesAtElement(1:4,iel)+nvt
        IdofV(9) = iel+nvt+nmt
        IdofP(1) = iel
        IdofP(2) = iel+nel
        IdofP(3) = iel+2*nel

        ! First of all, fetch the local RHS
        do i = 1, ndofV
          Df(i)       = p_DrhsU(IdofV(i))   ! f_u
          Df(ndofV+i) = p_DrhsV(IdofV(i))   ! f_v
        end do
        o = 2*ndofV
        do i = 1, ndofP
          Df(o+i) = p_DrhsP(IdofP(i))       ! f_p
        end do

        ! Clear the local system matrix
        Da = 0.0_DP

        ! Let us update the local RHS vector by subtracting A*u from it:
        ! f_u := f_u - A11*u
        ! f_v := f_v - A22*v
        p = ndofV
        do k = 1, ndofV
          i1 = p_KldA(IdofV(k))
          i2 = p_KldA(IdofV(k)+1)-1
          daux1 = 0.0_DP
          daux2 = 0.0_DP
          do i = i1, i2
            j = p_KcolA(i)
            daux1 = daux1 + p_DA11(i)*p_DvecU(j)
            daux2 = daux2 + p_DA22(i)*p_DvecV(j)
            do l = 1, ndofV
              if(j .eq. IdofV(l)) then
                Da(  k,  l) = Dmult(1,1)*p_DA11(i)
                Da(p+k,p+l) = Dmult(2,2)*p_DA22(i)
                exit
              end if
            end do
          end do
          Df(      k) = Df(      k) - Dmult(1,1)*daux1
          Df(ndofV+k) = Df(ndofV+k) - Dmult(2,2)*daux2
        end do

        if(bHaveA12) then
          ! f_u := f_u - A12*v
          ! f_v := f_v - A21*u
          do k = 1, ndofV
            i1 = p_KldA12(IdofV(k))
            i2 = p_KldA12(IdofV(k)+1)-1
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = i1, i2
              j = p_KcolA12(i)
              daux1 = daux1 + p_DA12(i)*p_DvecV(j)
              daux2 = daux2 + p_DA21(i)*p_DvecU(j)
              do l = 1, ndofV
                if(j .eq. IdofV(l)) then
                  Da(  k,p+l) = Dmult(1,2)*p_DA12(i)
                  Da(p+k,  l) = Dmult(2,1)*p_DA21(i)
                  exit
                end if
              end do
            end do
            Df(      k) = Df(      k) - Dmult(1,2)*daux1
            Df(ndofV+k) = Df(ndofV+k) - Dmult(2,1)*daux2
          end do
        end if ! bHaveA12

        ! Now we also need to subtract B*p from our RHS, and by the same time,
        ! we will build the local B matrices.
        ! f_u := f_u - B1*p
        ! f_v := f_v - B2*p
        o = ndofV
        p = 2*ndofV
        do k = 1, ndofV
          i1 = p_KldB(IdofV(k))
          i2 = p_KldB(IdofV(k)+1)-1
          daux1 = 0.0_DP
          daux2 = 0.0_DP
          do i = i1, i2
            j = p_KcolB(i)
            daux = p_DvecP(j)
            daux1 = daux1 + p_DB1(i)*daux
            daux2 = daux2 + p_DB2(i)*daux
            do l = 1, ndofP
              if(j .eq. IdofP(l)) then
                Da(  k,p+l) = Dmult(1,3)*p_DB1(i)
                Da(o+k,p+l) = Dmult(2,3)*p_DB2(i)
                exit
              end if
            end do
          end do
          Df(      k) = Df(      k) - Dmult(1,3)*daux1
          Df(ndofV+k) = Df(ndofV+k) - Dmult(2,3)*daux2
        end do

        ! Now we also need to subtract D*u from our RHS, and by the same time,
        ! we will build the local D matrices.
        ! f_p := f_p - D1*u - D2*v
        o = ndofV
        p = 2*ndofV
        do k = 1, ndofP
          i1 = p_KldD(IdofP(k))
          i2 = p_KldD(IdofP(k)+1)-1
          daux1 = 0.0_DP
          daux2 = 0.0_DP
          do i = i1, i2
            j = p_KcolD(i)
            daux1 = daux1 + p_DD1(i)*p_DvecU(j)
            daux2 = daux2 + p_DD2(i)*p_DvecV(j)
            do l = 1, ndofV
              if(j .eq. IdofV(l)) then
                Da(p+k,  l) = Dmult(3,1)*p_DD1(i)
                Da(p+k,o+l) = Dmult(3,2)*p_DD2(i)
                exit
              end if
            end do
          end do
          Df(p+k) = Df(p+k) - Dmult(3,1)*daux1 - Dmult(3,2)*daux2
        end do

        if(bHaveC) then
          ! f_p := f_p - C*p
          o = 2*ndofV
          do k = 1, ndofP
            i1 = p_KldC(IdofP(k))
            i2 = p_KldC(IdofP(k)+1)-1
            do i = i1, i2
              j = p_KcolC(i)
              Df(o+k) = Df(o+k) - Dmult(3,3)*p_DC(i)*p_DvecP(j)
              do l = 1, ndofP
                if(j .eq. IdofP(l)) then
                  Da(o+k,o+l) = Dmult(3,3)*p_DC(i)
                  exit
                end if
              end do
            end do
          end do
        end if ! bHaveC

        ! Solve the local system
        call DGESV(ndof,1,Da,ndof,Ipivot,Df,ndof,info)

        ! Did DGESV fail?
        if(info .ne. 0) cycle

        ! Incorporate our local solution into the global one.
        do i = 1, ndofV
          j = IdofV(i)
          p_DvecU(j) = p_DvecU(j) + domega * Df(i)
          p_DvecV(j) = p_DvecV(j) + domega * Df(ndofV+i)
        end do
        o = 2*ndofV
        do i = 1, ndofP
          j = IdofP(i)
          p_DvecP(j) = p_DvecP(j) + domega * Df(o+i)
        end do

      end do ! ielidx

    end do ! iter

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_NS2D_Q1TQ0_bd(rvanka, rsol, rrhs, niterations, domega, &
                                 IelementList)

!<description>
  ! Performs a desired number of iterations of the Vanka solver for
  ! 2D Navier-Stokes problem, "full" variant for Q1~/Q0 discretisations,
  ! no off-diagonal A matrices (A12,A21).
!</description>

!<input>
  ! t_vanka_NavSt2D structure that saves algorithm-specific parameters.
  type(t_vanka_NavSt2D), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! The number of iterations that are to be performed
  integer, intent(in) :: niterations

  ! Relaxation parameter.
  real(DP), intent(in) :: domega

  ! A list of element numbers where Vanka should be applied to.
  integer, dimension(:), intent(in) :: IelementList
!</input>

!<inputoutput>
  ! The iteration vector that is to be updated.
  type(t_vectorBlock), intent(inout) :: rsol
!</inputoutput>

!</subroutine>

  ! How many local DOFs do we have?
  integer, parameter :: ndofV = 4   ! Dofs per velocity
  integer, parameter :: ndofP = 1   ! Dofs per pressure
  integer, parameter :: ndof = 2*ndofV+ndofP

  ! Triangulation information
  type(t_triangulation), pointer :: p_rtria
  integer, dimension(:,:), pointer :: p_IedgesAtElement

  ! Data arrays of vectors
  real(DP), dimension(:), pointer :: p_DrhsU,p_DrhsV,p_DrhsP,&
                                     p_DvecU,p_DvecV,p_DvecP

  ! DOF-mapping arrays
  integer, dimension(ndofV) :: IdofV
  integer, dimension(ndofP) :: IdofP

  ! Variables for the local system
  real(DP), dimension(ndofV) :: Du1, Du2, Df1, Df2
  real(DP), dimension(ndofP) :: Dup, Dfp
  real(DP), dimension(ndofV,ndofV) :: Da1, Da2
  real(DP), dimension(ndofV,ndofP) :: Db1, Db2
  real(DP), dimension(ndofP,ndofV) :: Dd1, Dd2
  real(DP), dimension(ndofP,ndofP) :: Dc

  ! Local variables
  real(DP), dimension(ndofV,ndofV) :: Di1,Di2
  !real(DP), dimension(ndofP,ndofV) :: Dt1, Dt2
  real(DP), dimension(ndofV) :: Dt1, Dt2
  real(DP), dimension(ndofP,ndofP) :: Ds

  ! Multiplication factors
  real(DP), dimension(3,3) :: Dmult

  ! Quick access for the matrix arrays
  integer, dimension(:), pointer :: p_KldA,p_KldB,p_KldC,p_KldD,&
      p_KcolA,p_KcolB,p_KcolC,p_KcolD,p_KdiagA,p_KdiagC
  real(DP), dimension(:), pointer :: p_DA11,p_DA22,p_DB1,p_DB2,&
      p_DD1,p_DD2,p_DC

  ! local variables
  integer :: i,j,iel,ielidx,i1,i2,k,l,iter
  real(DP) :: daux,daux1,daux2
  logical :: bHaveC
  logical :: bsuccess1,bsuccess2

    ! Get the arrays from the triangulation
    p_rtria => rvanka%p_rspatialDiscrV%p_rtriangulation
    call storage_getbase_int2d (p_rtria%h_IedgesAtElement, p_IedgesAtElement)

    ! Get the pointers to the vector data
    call lsyssc_getbase_double(rsol%RvectorBlock(1), p_DvecU)
    call lsyssc_getbase_double(rsol%RvectorBlock(2), p_DvecV)
    call lsyssc_getbase_double(rsol%RvectorBlock(3), p_DvecP)
    call lsyssc_getbase_double(rrhs%RvectorBlock(1), p_DrhsU)
    call lsyssc_getbase_double(rrhs%RvectorBlock(2), p_DrhsV)
    call lsyssc_getbase_double(rrhs%RvectorBlock(3), p_DrhsP)

    ! Let us assume we do not have the optional matrices
    bHaveC = .false.

    ! Get the pointers from the vanka structure
    p_KldA => rvanka%p_KldA
    p_KcolA => rvanka%p_KcolA
    p_KdiagA => rvanka%p_KdiagonalA
    p_DA11 => rvanka%p_DA11
    p_DA22 => rvanka%p_DA22
    p_KldB => rvanka%p_KldB
    p_KcolB => rvanka%p_KcolB
    p_DB1 => rvanka%p_DB1
    p_DB2 => rvanka%p_DB2
    p_KldD => rvanka%p_KldD
    p_KcolD => rvanka%p_KcolD
    p_DD1 => rvanka%p_DD1
    p_DD2 => rvanka%p_DD2

    if(associated(rvanka%p_DC)) then
      bHaveC = .true.
      p_KldC => rvanka%p_KldC
      p_KcolC => rvanka%p_KcolC
      p_KdiagC => rvanka%p_KdiagonalC
      p_DC => rvanka%p_DC
    end if

    ! Get the multiplication factors
    Dmult = rvanka%Dmultipliers

    ! Take care of the "soft-deactivation" of the sub-matrices
    bHaveC = bHaveC .and. (Dmult(3,3) .ne. 0.0_DP)

    ! Clear the optional matrices
    Dc = 0.0_DP

    if(bHaveC) then

      do iter = 1, niterations

        ! C exists
        do ielidx=1, size(IelementList)

          ! Get the element number which is to be processed.
          iel = IelementList(ielidx)

          ! Get all DOFs for this element
          IdofV(:) = p_IedgesAtElement(:,iel)
          IdofP(1) = iel

          ! First of all, fetch the local RHS
          do i = 1, ndofV
            Df1(i) = p_DrhsU(IdofV(i))   ! f_u
            Df2(i) = p_DrhsV(IdofV(i))   ! f_v
          end do
          do i = 1, ndofP
            Dfp(i) = p_DrhsP(IdofP(i))   ! f_p
          end do

          ! Let us update the local RHS vector by subtracting A*u from it:
          ! f_u := f_u - A11*u
          ! f_v := f_v - A22*v
          do k = 1, ndofV
            i1 = p_KldA(IdofV(k))
            i2 = p_KldA(IdofV(k)+1)-1
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = i1, i2
              j = p_KcolA(i)
              daux1 = daux1 + p_DA11(i)*p_DvecU(j)
              daux2 = daux2 + p_DA22(i)*p_DvecV(j)
              do l = 1, ndofV
                if(j .eq. IdofV(l)) then
                  Da1(k,l) = Dmult(1,1)*p_DA11(i)
                  Da2(k,l) = Dmult(2,2)*p_DA22(i)
                  exit
                end if
              end do
            end do
            Df1(k) = Df1(k) - Dmult(1,1)*daux1
            Df2(k) = Df2(k) - Dmult(2,2)*daux2
          end do

          ! Now we also need to subtract B*p from our RHS, and by the same time,
          ! we will build the local B matrices.
          ! f_u := f_u - B1*p
          ! f_v := f_v - B2*p
          do k = 1, ndofV
            i1 = p_KldB(IdofV(k))
            i2 = p_KldB(IdofV(k)+1)-1
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = i1, i2
              j = p_KcolB(i)
              daux = p_DvecP(j)
              daux1 = daux1 + p_DB1(i)*daux
              daux2 = daux2 + p_DB2(i)*daux
              do l = 1, ndofP
                if(j .eq. IdofP(l)) then
                  Db1(k,l) = Dmult(1,3)*p_DB1(i)
                  Db2(k,l) = Dmult(2,3)*p_DB2(i)
                  exit
                end if
              end do
            end do
            Df1(k) = Df1(k) - Dmult(1,3)*daux1
            Df2(k) = Df2(k) - Dmult(2,3)*daux2
          end do

          ! Now we also need to subtract D*u from our RHS, and by the same time,
          ! we will build the local D matrices.
          ! f_p := f_p - D1*u - D2*v
          do k = 1, ndofP
            i1 = p_KldD(IdofP(k))
            i2 = p_KldD(IdofP(k)+1)-1
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = i1, i2
              j = p_KcolD(i)
              daux = p_DvecP(j)
              daux1 = daux1 + p_DD1(i)*p_DvecU(j)
              daux2 = daux2 + p_DD2(i)*p_DvecV(j)
              do l = 1, ndofV
                if(j .eq. IdofV(l)) then
                  Dd1(k,l) = Dmult(3,1)*p_DD1(i)
                  Dd2(k,l) = Dmult(3,2)*p_DD2(i)
                  exit
                end if
              end do
            end do
            Dfp(k) = Dfp(k) - Dmult(3,1)*daux1 - Dmult(3,2)*daux2
          end do

          ! Let us update the local RHS vector by subtracting C*p from it:
          ! f_p := f_p - C*p
          do k = 1, ndofP
            i1 = p_KldC(IdofP(k))
            i2 = p_KldC(IdofP(k)+1)-1
            daux1 = 0.0_DP
            do i = i1, i2
              j = p_KcolC(i)
              daux1 = daux1 + p_DC(i)*p_DvecP(j)
              do l = 1, ndofP
                if(j .eq. IdofP(l)) then
                  Dc(k,l) = Dmult(3,3)*p_DC(i)
                  exit
                end if
              end do
            end do
            Dfp(k) = Dfp(k) - Dmult(3,3)*daux1
          end do

          ! Invert A1 and A2
          call mprim_invert4x4MatrixDirect(Da1, Di1,bsuccess1)
          call mprim_invert4x4MatrixDirect(Da2, Di2,bsuccess2)

          if (bsuccess1 .and. bsuccess2) then

            ! Precalculate D * A^-1
            Dt1(1) = Dd1(1,1)*Di1(1,1)+Dd1(1,2)*Di1(2,1)+Dd1(1,3)*Di1(3,1)+Dd1(1,4)*Di1(4,1)
            Dt1(2) = Dd1(1,1)*Di1(1,2)+Dd1(1,2)*Di1(2,2)+Dd1(1,3)*Di1(3,2)+Dd1(1,4)*Di1(4,2)
            Dt1(3) = Dd1(1,1)*Di1(1,3)+Dd1(1,2)*Di1(2,3)+Dd1(1,3)*Di1(3,3)+Dd1(1,4)*Di1(4,3)
            Dt1(4) = Dd1(1,1)*Di1(1,4)+Dd1(1,2)*Di1(2,4)+Dd1(1,3)*Di1(3,4)+Dd1(1,4)*Di1(4,4)
            Dt2(1) = Dd2(1,1)*Di2(1,1)+Dd2(1,2)*Di2(2,1)+Dd2(1,3)*Di2(3,1)+Dd2(1,4)*Di2(4,1)
            Dt2(2) = Dd2(1,1)*Di2(1,2)+Dd2(1,2)*Di2(2,2)+Dd2(1,3)*Di2(3,2)+Dd2(1,4)*Di2(4,2)
            Dt2(3) = Dd2(1,1)*Di2(1,3)+Dd2(1,2)*Di2(2,3)+Dd2(1,3)*Di2(3,3)+Dd2(1,4)*Di2(4,3)
            Dt2(4) = Dd2(1,1)*Di2(1,4)+Dd2(1,2)*Di2(2,4)+Dd2(1,3)*Di2(3,4)+Dd2(1,4)*Di2(4,4)

            ! Calculate Schur-Complement of A
            ! S := -C + D * A^-1 * B
            Ds(1,1) = -Dc(1,1) &
                    + Dt1(1)*Db1(1,1)+Dt1(2)*Db1(2,1)+Dt1(3)*Db1(3,1)+Dt1(4)*Db1(4,1)&
                    + Dt2(1)*Db2(1,1)+Dt2(2)*Db2(2,1)+Dt2(3)*Db2(3,1)+Dt2(4)*Db2(4,1)

            ! Calculate pressure
            ! p := S^-1 * (D * A^-1 * f_u - f_p)
            Dup(1) = (-Dfp(1) &
                  + Dt1(1)*Df1(1)+Dt1(2)*Df1(2)+Dt1(3)*Df1(3)+Dt1(4)*Df1(4) &
                  + Dt2(1)*Df2(1)+Dt2(2)*Df2(2)+Dt2(3)*Df2(3)+Dt2(4)*Df2(4)) / Ds(1,1)

            ! Update RHS
            ! f_u := f_u - B * p
            Df1(1) = Df1(1) - Db1(1,1)*Dup(1)
            Df1(2) = Df1(2) - Db1(2,1)*Dup(1)
            Df1(3) = Df1(3) - Db1(3,1)*Dup(1)
            Df1(4) = Df1(4) - Db1(4,1)*Dup(1)
            Df2(1) = Df2(1) - Db2(1,1)*Dup(1)
            Df2(2) = Df2(2) - Db2(2,1)*Dup(1)
            Df2(3) = Df2(3) - Db2(3,1)*Dup(1)
            Df2(4) = Df2(4) - Db2(4,1)*Dup(1)

            ! Calculate X- and Y-velocity
            ! u := A^-1 * f_u
            Du1(1) = Di1(1,1)*Df1(1)+Di1(1,2)*Df1(2)+Di1(1,3)*Df1(3)+Di1(1,4)*Df1(4)
            Du1(2) = Di1(2,1)*Df1(1)+Di1(2,2)*Df1(2)+Di1(2,3)*Df1(3)+Di1(2,4)*Df1(4)
            Du1(3) = Di1(3,1)*Df1(1)+Di1(3,2)*Df1(2)+Di1(3,3)*Df1(3)+Di1(3,4)*Df1(4)
            Du1(4) = Di1(4,1)*Df1(1)+Di1(4,2)*Df1(2)+Di1(4,3)*Df1(3)+Di1(4,4)*Df1(4)
            Du2(1) = Di2(1,1)*Df2(1)+Di2(1,2)*Df2(2)+Di2(1,3)*Df2(3)+Di2(1,4)*Df2(4)
            Du2(2) = Di2(2,1)*Df2(1)+Di2(2,2)*Df2(2)+Di2(2,3)*Df2(3)+Di2(2,4)*Df2(4)
            Du2(3) = Di2(3,1)*Df2(1)+Di2(3,2)*Df2(2)+Di2(3,3)*Df2(3)+Di2(3,4)*Df2(4)
            Du2(4) = Di2(4,1)*Df2(1)+Di2(4,2)*Df2(2)+Di2(4,3)*Df2(3)+Di2(4,4)*Df2(4)

            ! Incorporate our local solution into the global one.
            do i = 1, ndofV
              j = IdofV(i)
              p_DvecU(j) = p_DvecU(j) + domega * Du1(i)
              p_DvecV(j) = p_DvecV(j) + domega * Du2(i)
            end do
            do i = 1, ndofP
              j = IdofP(i)
              p_DvecP(j) = p_DvecP(j) + domega * Dup(i)
            end do

          end if

        end do ! ielidx

      end do ! iter

    else

      do iter = 1, niterations

        ! C does not exist
        do ielidx=1, size(IelementList)

          ! Get the element number which is to be processed.
          iel = IelementList(ielidx)

          ! Get all DOFs for this element
          IdofV(:) = p_IedgesAtElement(:,iel)
          IdofP(1) = iel

          ! First of all, fetch the local RHS
          do i = 1, ndofV
            Df1(i) = p_DrhsU(IdofV(i))   ! f_u
            Df2(i) = p_DrhsV(IdofV(i))   ! f_v
          end do
          do i = 1, ndofP
            Dfp(i) = p_DrhsP(IdofP(i))   ! f_p
          end do

          ! Let us update the local RHS vector by subtracting A*u from it:
          ! f_u := f_u - A11*u
          ! f_v := f_v - A22*v
          do k = 1, ndofV
            i1 = p_KldA(IdofV(k))
            i2 = p_KldA(IdofV(k)+1)-1
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = i1, i2
              j = p_KcolA(i)
              daux1 = daux1 + p_DA11(i)*p_DvecU(j)
              daux2 = daux2 + p_DA22(i)*p_DvecV(j)
              do l = 1, ndofV
                if(j .eq. IdofV(l)) then
                  Da1(k,l) = Dmult(1,1)*p_DA11(i)
                  Da2(k,l) = Dmult(2,2)*p_DA22(i)
                  exit
                end if
              end do
            end do
            Df1(k) = Df1(k) - Dmult(1,1)*daux1
            Df2(k) = Df2(k) - Dmult(2,2)*daux2
          end do

          ! Now we also need to subtract B*p from our RHS, and by the same time,
          ! we will build the local B matrices.
          ! f_u := f_u - B1*p
          ! f_v := f_v - B2*p
          do k = 1, ndofV
            i1 = p_KldB(IdofV(k))
            i2 = p_KldB(IdofV(k)+1)-1
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = i1, i2
              j = p_KcolB(i)
              daux = p_DvecP(j)
              daux1 = daux1 + p_DB1(i)*daux
              daux2 = daux2 + p_DB2(i)*daux
              do l = 1, ndofP
                if(j .eq. IdofP(l)) then
                  Db1(k,l) = Dmult(1,3)*p_DB1(i)
                  Db2(k,l) = Dmult(2,3)*p_DB2(i)
                  exit
                end if
              end do
            end do
            Df1(k) = Df1(k) - Dmult(1,3)*daux1
            Df2(k) = Df2(k) - Dmult(2,3)*daux2
          end do

          ! Now we also need to subtract D*u from our RHS, and by the same time,
          ! we will build the local D matrices.
          ! f_p := f_p - D1*u - D2*v
          do k = 1, ndofP
            i1 = p_KldD(IdofP(k))
            i2 = p_KldD(IdofP(k)+1)-1
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = i1, i2
              j = p_KcolD(i)
              daux1 = daux1 + p_DD1(i)*p_DvecU(j)
              daux2 = daux2 + p_DD2(i)*p_DvecV(j)
              do l = 1, ndofV
                if(j .eq. IdofV(l)) then
                  Dd1(k,l) = Dmult(3,1)*p_DD1(i)
                  Dd2(k,l) = Dmult(3,2)*p_DD2(i)
                  exit
                end if
              end do
            end do
            Dfp(k) = Dfp(k) - Dmult(3,1)*daux1 - Dmult(3,2)*daux2
          end do

          ! Invert A1 and A2
          call mprim_invert4x4MatrixDirect(Da1, Di1,bsuccess1)
          call mprim_invert4x4MatrixDirect(Da2, Di2,bsuccess2)

          if (bsuccess1 .and. bsuccess2) then

            ! Precalculate D * A^-1
            Dt1(1) = Dd1(1,1)*Di1(1,1)+Dd1(1,2)*Di1(2,1)+Dd1(1,3)*Di1(3,1)+Dd1(1,4)*Di1(4,1)
            Dt1(2) = Dd1(1,1)*Di1(1,2)+Dd1(1,2)*Di1(2,2)+Dd1(1,3)*Di1(3,2)+Dd1(1,4)*Di1(4,2)
            Dt1(3) = Dd1(1,1)*Di1(1,3)+Dd1(1,2)*Di1(2,3)+Dd1(1,3)*Di1(3,3)+Dd1(1,4)*Di1(4,3)
            Dt1(4) = Dd1(1,1)*Di1(1,4)+Dd1(1,2)*Di1(2,4)+Dd1(1,3)*Di1(3,4)+Dd1(1,4)*Di1(4,4)
            Dt2(1) = Dd2(1,1)*Di2(1,1)+Dd2(1,2)*Di2(2,1)+Dd2(1,3)*Di2(3,1)+Dd2(1,4)*Di2(4,1)
            Dt2(2) = Dd2(1,1)*Di2(1,2)+Dd2(1,2)*Di2(2,2)+Dd2(1,3)*Di2(3,2)+Dd2(1,4)*Di2(4,2)
            Dt2(3) = Dd2(1,1)*Di2(1,3)+Dd2(1,2)*Di2(2,3)+Dd2(1,3)*Di2(3,3)+Dd2(1,4)*Di2(4,3)
            Dt2(4) = Dd2(1,1)*Di2(1,4)+Dd2(1,2)*Di2(2,4)+Dd2(1,3)*Di2(3,4)+Dd2(1,4)*Di2(4,4)

            ! Calculate Schur-Complement of A
            ! S := D * A^-1 * B
            Ds(1,1) = Dt1(1)*Db1(1,1)+Dt1(2)*Db1(2,1)+Dt1(3)*Db1(3,1)+Dt1(4)*Db1(4,1)&
                    + Dt2(1)*Db2(1,1)+Dt2(2)*Db2(2,1)+Dt2(3)*Db2(3,1)+Dt2(4)*Db2(4,1)

            ! Calculate pressure
            ! p := S^-1 * (D * A^-1 * f_u - f_p)
            Dup(1) = (-Dfp(1) &
                  + Dt1(1)*Df1(1)+Dt1(2)*Df1(2)+Dt1(3)*Df1(3)+Dt1(4)*Df1(4) &
                  + Dt2(1)*Df2(1)+Dt2(2)*Df2(2)+Dt2(3)*Df2(3)+Dt2(4)*Df2(4)) / Ds(1,1)

            ! Update RHS
            ! f_u := f_u - B * p
            Df1(1) = Df1(1) - Db1(1,1)*Dup(1)
            Df1(2) = Df1(2) - Db1(2,1)*Dup(1)
            Df1(3) = Df1(3) - Db1(3,1)*Dup(1)
            Df1(4) = Df1(4) - Db1(4,1)*Dup(1)
            Df2(1) = Df2(1) - Db2(1,1)*Dup(1)
            Df2(2) = Df2(2) - Db2(2,1)*Dup(1)
            Df2(3) = Df2(3) - Db2(3,1)*Dup(1)
            Df2(4) = Df2(4) - Db2(4,1)*Dup(1)

            ! Calculate X- and Y-velocity
            ! u := A^-1 * f_u
            Du1(1) = Di1(1,1)*Df1(1)+Di1(1,2)*Df1(2)+Di1(1,3)*Df1(3)+Di1(1,4)*Df1(4)
            Du1(2) = Di1(2,1)*Df1(1)+Di1(2,2)*Df1(2)+Di1(2,3)*Df1(3)+Di1(2,4)*Df1(4)
            Du1(3) = Di1(3,1)*Df1(1)+Di1(3,2)*Df1(2)+Di1(3,3)*Df1(3)+Di1(3,4)*Df1(4)
            Du1(4) = Di1(4,1)*Df1(1)+Di1(4,2)*Df1(2)+Di1(4,3)*Df1(3)+Di1(4,4)*Df1(4)
            Du2(1) = Di2(1,1)*Df2(1)+Di2(1,2)*Df2(2)+Di2(1,3)*Df2(3)+Di2(1,4)*Df2(4)
            Du2(2) = Di2(2,1)*Df2(1)+Di2(2,2)*Df2(2)+Di2(2,3)*Df2(3)+Di2(2,4)*Df2(4)
            Du2(3) = Di2(3,1)*Df2(1)+Di2(3,2)*Df2(2)+Di2(3,3)*Df2(3)+Di2(3,4)*Df2(4)
            Du2(4) = Di2(4,1)*Df2(1)+Di2(4,2)*Df2(2)+Di2(4,3)*Df2(3)+Di2(4,4)*Df2(4)

            ! Incorporate our local solution into the global one.
            do i = 1, ndofV
              j = IdofV(i)
              p_DvecU(j) = p_DvecU(j) + domega * Du1(i)
              p_DvecV(j) = p_DvecV(j) + domega * Du2(i)
            end do
            do i = 1, ndofP
              j = IdofP(i)
              p_DvecP(j) = p_DvecP(j) + domega * Dup(i)
            end do

          end if

        end do ! ielidx

      end do ! iter

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_NS2D_Q1TQ0_fc(rvanka, rsol, rrhs, niterations, domega, &
                                 IelementList)

!<description>
  ! Performs a desired number of iterations of the Vanka solver for
  ! 2D Navier-Stokes problem, "full" variant for Q1~/Q0 discretisations,
  ! off-diagonal A-matrices exist (A12,A21).
!</description>

!<input>
  ! t_vanka_NavSt2D structure that saves algorithm-specific parameters.
  type(t_vanka_NavSt2D), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs

  ! The number of iterations that are to be performed
  integer, intent(in) :: niterations

  ! Relaxation parameter.
  real(DP), intent(in) :: domega

  ! A list of element numbers where VANKA should be applied to.
  integer, dimension(:), intent(in) :: IelementList
!</input>

!<inputoutput>
  ! The iteration vector that is to be updated.
  type(t_vectorBlock), intent(inout) :: rsol
!</inputoutput>

!</subroutine>

  ! How many local DOFs do we have?
  integer, parameter :: ndofV = 4   ! Dofs per velocity
  integer, parameter :: ndofP = 1   ! Dofs per pressure
  integer, parameter :: ndof = 2*ndofV+ndofP

  ! Triangulation information
  type(t_triangulation), pointer :: p_rtria
  integer, dimension(:,:), pointer :: p_IedgesAtElement

  ! Data arrays of vectors
  real(DP), dimension(:), pointer :: p_DrhsU,p_DrhsV,p_DrhsP,&
                                     p_DvecU,p_DvecV,p_DvecP

  ! DOF-mapping arrays
  integer, dimension(ndofV) :: IdofV
  integer, dimension(ndofP) :: IdofP

  ! Variables for the local system
  real(DP), dimension(ndof) :: Df
  real(DP), dimension(ndof,ndof) :: Da

  ! Multiplication factors
  real(DP), dimension(3,3) :: Dmult

  ! Quick access for the matrix arrays
  integer, dimension(:), pointer :: p_KldA,p_KldA12,p_KldB,p_KldC,p_KldD,&
      p_KcolA,p_KcolA12,p_KcolB,p_KcolC,p_KcolD
  real(DP), dimension(:), pointer :: p_DA11,p_DA12,p_DA21,p_DA22,p_DB1,p_DB2,&
      p_DD1,p_DD2,p_DC

  ! local variables
  integer :: i,j,iel,ielidx,i1,i2,k,l,o,p,iter
  real(DP) :: daux,daux1,daux2
  logical :: bHaveC

  ! variables for LAPACK`s DGESV routine
  integer, dimension(ndof) :: Ipivot
  integer :: info

    ! Get the arrays from the triangulation
    p_rtria => rvanka%p_rspatialDiscrV%p_rtriangulation
    call storage_getbase_int2d (p_rtria%h_IedgesAtElement, p_IedgesAtElement)

    ! Get the pointers to the vector data
    call lsyssc_getbase_double(rsol%RvectorBlock(1), p_DvecU)
    call lsyssc_getbase_double(rsol%RvectorBlock(2), p_DvecV)
    call lsyssc_getbase_double(rsol%RvectorBlock(3), p_DvecP)
    call lsyssc_getbase_double(rrhs%RvectorBlock(1), p_DrhsU)
    call lsyssc_getbase_double(rrhs%RvectorBlock(2), p_DrhsV)
    call lsyssc_getbase_double(rrhs%RvectorBlock(3), p_DrhsP)

    ! Let us assume we do not have the optional matrices
    bHaveC = .false.

    ! Get the pointers from the vanka structure
    p_KldA => rvanka%p_KldA
    p_KcolA => rvanka%p_KcolA
    p_DA11 => rvanka%p_DA11
    p_DA22 => rvanka%p_DA22
    p_KldA12 => rvanka%p_KldA12
    p_KcolA12 => rvanka%p_KcolA12
    p_DA12 => rvanka%p_DA12
    p_DA21 => rvanka%p_DA21
    p_KldB => rvanka%p_KldB
    p_KcolB => rvanka%p_KcolB
    p_DB1 => rvanka%p_DB1
    p_DB2 => rvanka%p_DB2
    p_KldD => rvanka%p_KldD
    p_KcolD => rvanka%p_KcolD
    p_DD1 => rvanka%p_DD1
    p_DD2 => rvanka%p_DD2

    if(associated(rvanka%p_DC)) then
      bHaveC = .true.
      p_KldC => rvanka%p_KldC
      p_KcolC => rvanka%p_KcolC
      p_DC => rvanka%p_DC
    end if

    ! Get the multiplication factors
    Dmult = rvanka%Dmultipliers

    ! Take care of the "soft-deactivation" of the sub-matrices
    bHaveC = bHaveC .and. (Dmult(3,3) .ne. 0.0_DP)

    if(bHaveC) then

      do iter = 1, niterations

        ! C exists
        do ielidx = 1, size(IelementList)

          ! Get the element number which is to be processed.
          iel = IelementList(ielidx)

          ! Get all doFs for this element
          IdofV(:) = p_IedgesAtElement(:,iel)
          IdofP(1) = iel

          ! First of all, fetch the local RHS
          do i = 1, ndofV
            Df(i)       = p_DrhsU(IdofV(i))   ! f_u
            Df(ndofV+i) = p_DrhsV(IdofV(i))   ! f_v
          end do
          o = 2*ndofV
          do i = 1, ndofP
            Df(o+i) = p_DrhsP(IdofP(i))       ! f_p
          end do

          ! Clear the local system matrix
          Da = 0.0_DP

          ! Let us update the local RHS vector by subtracting A*u from it:
          ! f_u := f_u - A11*u
          ! f_v := f_v - A22*v
          p = ndofV
          do k = 1, ndofV
            i1 = p_KldA(IdofV(k))
            i2 = p_KldA(IdofV(k)+1)-1
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = i1, i2
              j = p_KcolA(i)
              daux1 = daux1 + p_DA11(i)*p_DvecU(j)
              daux2 = daux2 + p_DA22(i)*p_DvecV(j)
              do l = 1, ndofV
                if(j .eq. IdofV(l)) then
                  Da(  k,  l) = Dmult(1,1)*p_DA11(i)
                  Da(p+k,p+l) = Dmult(2,2)*p_DA22(i)
                  exit
                end if
              end do
            end do
            Df(      k) = Df(      k) - Dmult(1,1)*daux1
            Df(ndofV+k) = Df(ndofV+k) - Dmult(2,2)*daux2
          end do

          ! f_u := f_u - A12*v
          ! f_v := f_v - A21*u
          do k = 1, ndofV
            i1 = p_KldA12(IdofV(k))
            i2 = p_KldA12(IdofV(k)+1)-1
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = i1, i2
              j = p_KcolA12(i)
              daux1 = daux1 + p_DA12(i)*p_DvecV(j)
              daux2 = daux2 + p_DA21(i)*p_DvecU(j)
              do l = 1, ndofV
                if(j .eq. IdofV(l)) then
                  Da(  k,p+l) = Dmult(1,2)*p_DA12(i)
                  Da(p+k,  l) = Dmult(2,1)*p_DA21(i)
                  exit
                end if
              end do
            end do
            Df(      k) = Df(      k) - Dmult(1,2)*daux1
            Df(ndofV+k) = Df(ndofV+k) - Dmult(2,1)*daux2
          end do

          ! Now we also need to subtract B*p from our RHS, and by the same time,
          ! we will build the local B matrices.
          ! f_u := f_u - B1*p
          ! f_v := f_v - B2*p
          o = ndofV
          p = 2*ndofV
          do k = 1, ndofV
            i1 = p_KldB(IdofV(k))
            i2 = p_KldB(IdofV(k)+1)-1
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = i1, i2
              j = p_KcolB(i)
              daux = p_DvecP(j)
              daux1 = daux1 + p_DB1(i)*daux
              daux2 = daux2 + p_DB2(i)*daux
              do l = 1, ndofP
                if(j .eq. IdofP(l)) then
                  Da(  k,p+l) = Dmult(1,3)*p_DB1(i)
                  Da(o+k,p+l) = Dmult(2,3)*p_DB2(i)
                  exit
                end if
              end do
            end do
            Df(      k) = Df(      k) - Dmult(1,3)*daux1
            Df(ndofV+k) = Df(ndofV+k) - Dmult(2,3)*daux2
          end do

          ! Now we also need to subtract D*u from our RHS, and by the same time,
          ! we will build the local D matrices.
          ! f_p := f_p - D1*u - D2*v
          o = ndofV
          p = 2*ndofV
          do k = 1, ndofP
            i1 = p_KldD(IdofP(k))
            i2 = p_KldD(IdofP(k)+1)-1
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = i1, i2
              j = p_KcolD(i)
              daux1 = daux1 + p_DD1(i)*p_DvecU(j)
              daux2 = daux2 + p_DD2(i)*p_DvecV(j)
              do l = 1, ndofV
                if(j .eq. IdofV(l)) then
                  Da(p+k,  l) = Dmult(3,1)*p_DD1(i)
                  Da(p+k,o+l) = Dmult(3,2)*p_DD2(i)
                  exit
                end if
              end do
            end do
            Df(p+k) = Df(p+k) - Dmult(3,1)*daux1 - Dmult(3,2)*daux2
          end do

          ! f_p := f_p - C*p
          o = 2*ndofV
          do k = 1, ndofP
            i1 = p_KldC(IdofP(k))
            i2 = p_KldC(IdofP(k)+1)-1
            do i = i1, i2
              j = p_KcolC(i)
              Df(o+k) = Df(o+k) - Dmult(3,3)*p_DC(i)*p_DvecP(j)
              do l = 1, ndofP
                if(j .eq. IdofP(l)) then
                  Da(o+k,o+l) = Dmult(3,3)*p_DC(i)
                  exit
                end if
              end do
            end do
          end do

          ! Solve the local system
          call DGESV(ndof,1,Da,ndof,Ipivot,Df,ndof,info)

          ! Did DGESV fail?
          if(info .ne. 0) cycle

          ! Incorporate our local solution into the global one.
          do i = 1, ndofV
            j = IdofV(i)
            p_DvecU(j) = p_DvecU(j) + domega * Df(i)
            p_DvecV(j) = p_DvecV(j) + domega * Df(ndofV+i)
          end do
          o = 2*ndofV
          do i = 1, ndofP
            j = IdofP(i)
            p_DvecP(j) = p_DvecP(j) + domega * Df(o+i)
          end do

        end do ! ielidx

      end do ! iter

    else

      do iter = 1, niterations

        ! C exists
        do ielidx=1, size(IelementList)

          ! Get the element number which is to be processed.
          iel = IelementList(ielidx)

          ! Get all DOFs for this element
          IdofV(:) = p_IedgesAtElement(:,iel)
          IdofP(1) = iel

          ! First of all, fetch the local RHS
          do i = 1, ndofV
            Df(i)       = p_DrhsU(IdofV(i))   ! f_u
            Df(ndofV+i) = p_DrhsV(IdofV(i))   ! f_v
          end do
          o = 2*ndofV
          do i = 1, ndofP
            Df(o+i) = p_DrhsP(IdofP(i))       ! f_p
          end do

          ! Clear the local system matrix
          Da = 0.0_DP

          ! Let us update the local RHS vector by subtracting A*u from it:
          ! f_u := f_u - A11*u
          ! f_v := f_v - A22*v
          p = ndofV
          do k = 1, ndofV
            i1 = p_KldA(IdofV(k))
            i2 = p_KldA(IdofV(k)+1)-1
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = i1, i2
              j = p_KcolA(i)
              daux1 = daux1 + p_DA11(i)*p_DvecU(j)
              daux2 = daux2 + p_DA22(i)*p_DvecV(j)
              do l = 1, ndofV
                if(j .eq. IdofV(l)) then
                  Da(  k,  l) = Dmult(1,1)*p_DA11(i)
                  Da(p+k,p+l) = Dmult(2,2)*p_DA22(i)
                  exit
                end if
              end do
            end do
            Df(      k) = Df(      k) - Dmult(1,1)*daux1
            Df(ndofV+k) = Df(ndofV+k) - Dmult(2,2)*daux2
          end do

          ! f_u := f_u - A12*v
          ! f_v := f_v - A21*u
          do k = 1, ndofV
            i1 = p_KldA12(IdofV(k))
            i2 = p_KldA12(IdofV(k)+1)-1
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = i1, i2
              j = p_KcolA12(i)
              daux1 = daux1 + p_DA12(i)*p_DvecV(j)
              daux2 = daux2 + p_DA21(i)*p_DvecU(j)
              do l = 1, ndofV
                if(j .eq. IdofV(l)) then
                  Da(  k,p+l) = Dmult(1,2)*p_DA12(i)
                  Da(p+k,  l) = Dmult(2,1)*p_DA21(i)
                  exit
                end if
              end do
            end do
            Df(      k) = Df(      k) - Dmult(1,2)*daux1
            Df(ndofV+k) = Df(ndofV+k) - Dmult(2,1)*daux2
          end do

          ! Now we also need to subtract B*p from our RHS, and by the same time,
          ! we will build the local B matrices.
          ! f_u := f_u - B1*p
          ! f_v := f_v - B2*p
          o = ndofV
          p = 2*ndofV
          do k = 1, ndofV
            i1 = p_KldB(IdofV(k))
            i2 = p_KldB(IdofV(k)+1)-1
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = i1, i2
              j = p_KcolB(i)
              daux = p_DvecP(j)
              daux1 = daux1 + p_DB1(i)*daux
              daux2 = daux2 + p_DB2(i)*daux2
              do l = 1, ndofP
                if(j .eq. IdofP(l)) then
                  Da(  k,p+l) = Dmult(1,3)*p_DB1(i)
                  Da(o+k,p+l) = Dmult(2,3)*p_DB2(i)
                  exit
                end if
              end do
            end do
            Df(      k) = Df(      k) - Dmult(1,3)*daux1
            Df(ndofV+k) = Df(ndofV+k) - Dmult(2,3)*daux2
          end do

          ! Now we also need to subtract D*u from our RHS, and by the same time,
          ! we will build the local D matrices.
          ! f_p := f_p - D1*u - D2*v
          o = ndofV
          p = 2*ndofV
          do k = 1, ndofP
            i1 = p_KldD(IdofP(k))
            i2 = p_KldD(IdofP(k)+1)-1
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = i1, i2
              j = p_KcolD(i)
              daux1 = daux1 + p_DD1(i)*p_DvecU(j)
              daux2 = daux2 + p_DD2(i)*p_DvecV(j)
              do l = 1, ndofV
                if(j .eq. IdofV(l)) then
                  Da(p+k,  l) = Dmult(3,1)*p_DD1(i)
                  Da(p+k,o+l) = Dmult(3,2)*p_DD2(i)
                  exit
                end if
              end do
            end do
            Df(p+k) = Df(p+k) - Dmult(3,1)*daux1 - Dmult(3,2)*daux2
          end do

          ! Solve the local system
          call DGESV(ndof,1,Da,ndof,Ipivot,Df,ndof,info)

          ! Did DGESV fail?
          if(info .ne. 0) cycle

          ! Incorporate our local solution into the global one.
          do i = 1, ndofV
            j = IdofV(i)
            p_DvecU(j) = p_DvecU(j) + domega * Df(i)
            p_DvecV(j) = p_DvecV(j) + domega * Df(ndofV+i)
          end do
          o = 2*ndofV
          do i = 1, ndofP
            j = IdofP(i)
            p_DvecP(j) = p_DvecP(j) + domega * Df(o+i)
          end do

        end do ! ielidx

      end do ! iter

    end if

  end subroutine

  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !
  !      T E M P L A T E S   F O R   N E W   I M P L E M E N T A T I O N S
  !
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ! $TEMPLATE$
  ! To add a new 2D Navier-Stokes Vanka implementation for a specific element
  ! combination you can use the templates below.
  !
  ! To use a template, do the following steps:
  ! 1. Copy-&-Paste the template
  ! 2. Uncomment the template
  ! 3. Replace all $TODO$ tags of the template by the necessary information
  !    that fits your new implementation
  ! 4. Optional: Customise the implementation
  ! 5. Add a check for your discretisation in vanka_init_NavSt2D
  !    (see $TODO$ tag in that routine)
  ! 6. Add your discretisation to vanka_solve_NavSt2D and call your the
  !    routine you created from the template (see $TODO$ tag in that routine)

  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  ! $TEMPLATE$ for "diagonal" Vanka

!<!-- // hide from automatic documentation parser
!<subroutine>
!
!  subroutine vanka_NS2D_$TODO$(rvanka, rsol, rrhs, niterations, domega, &
!                               IelementList)
!
!!<description>
!  ! Performs a desired number of iterations of the Vanka solver for
!  ! 2D Navier-Stokes problem, "diagonal" variant for $TODO$ discretisations.
!!</description>
!
!!<input>
!  ! t_vanka_NavSt2D structure that saves algorithm-specific parameters.
!  type(t_vanka_NavSt2D), intent(in) :: rvanka
!
!  ! The right-hand-side vector of the system
!  type(t_vectorBlock), intent(in) :: rrhs
!
!  ! The number of iterations that are to be performed
!  integer, intent(in) :: niterations
!
!  ! Relaxation parameter.
!  real(DP), intent(in) :: domega
!
!  ! A list of element numbers where Vanka should be applied to.
!  integer, dimension(:), intent(in) :: IelementList
!!</input>
!
!!<inputoutput>
!  ! The iteration vector that is to be updated.
!  type(t_vectorBlock), intent(inout) :: rsol
!!</inputoutput>
!
!!</subroutine>
!
!  ! How many local DOFs do we have?
!  integer, parameter :: ndofV = $TODO$   ! Dofs per velocity
!  integer, parameter :: ndofP = $TODO$   ! Dofs per pressure
!  integer, parameter :: ndof = 2*ndofV+ndofP
!
!  ! Triangulation information
!  type(t_triangulation), pointer :: p_rtria
!  $TODO$: Declare the necessary variables that need to be accessed from
!          the triangulation to perform the DOF-mapping.
!
!  ! Data arrays of vectors
!  real(DP), dimension(:), pointer :: p_DrhsU,p_DrhsV,p_DrhsP,&
!                                     p_DvecU,p_DvecV,p_DvecP
!
!  ! DOF-mapping arrays
!  integer, dimension(ndofV) :: IdofV
!  integer, dimension(ndofP) :: IdofP
!
!  ! Variables for the local system
!  real(DP), dimension(ndofV) :: Da1, Da2, Du1, Du2, Df1, Df2
!  real(DP), dimension(ndofV,ndofP) :: Db1, Db2
!  real(DP), dimension(ndofP,ndofV) :: Dd1, Dd2
!  real(DP), dimension(ndofP) :: Ds, Dc, Dup, Dfp
!
!  ! Multiplication factors
!  real(DP), dimension(3,3) :: Dmult
!
!  ! Quick access for the matrix arrays
!  integer, dimension(:), pointer :: p_KldA,p_KldA12,p_KldB,p_KldC,p_KldD,&
!      p_KcolA,p_KcolA12,p_KcolB,p_KcolC,p_KcolD,p_KdiagA,p_KdiagC
!  real(DP), dimension(:), pointer :: p_DA11,p_DA12,p_DA21,p_DA22,p_DB1,p_DB2,&
!      p_DD1,p_DD2,p_DC
!
!  ! local variables
!  integer :: i,j,iel,ielidx,i1,i2,k,l,iter
!  real(DP) :: daux,daux1,daux2
!  logical :: bHaveA12, bHaveC
!
!    ! Get the arrays from the triangulation
!    p_rtria => rvanka%p_rspatialDiscrV%p_rtriangulation
!    $TODO$ Get the arrays from the triangulation that are necessary to
!           perform the DOF-mapping
!
!    ! Get the pointers to the vector data
!    call lsyssc_getbase_double(rsol%RvectorBlock(1), p_DvecU)
!    call lsyssc_getbase_double(rsol%RvectorBlock(2), p_DvecV)
!    call lsyssc_getbase_double(rsol%RvectorBlock(3), p_DvecP)
!    call lsyssc_getbase_double(rrhs%RvectorBlock(1), p_DrhsU)
!    call lsyssc_getbase_double(rrhs%RvectorBlock(2), p_DrhsV)
!    call lsyssc_getbase_double(rrhs%RvectorBlock(3), p_DrhsP)
!
!    ! Let us assume we do not have the optional matrices
!    bHaveA12 = .FALSE.
!    bHaveC = .FALSE.
!
!    ! Get the pointers from the vanka structure
!    p_KldA => rvanka%p_KldA
!    p_KcolA => rvanka%p_KcolA
!    p_KdiagA => rvanka%p_KdiagonalA
!    p_DA11 => rvanka%p_DA11
!    p_DA22 => rvanka%p_DA22
!    p_KldB => rvanka%p_KldB
!    p_KcolB => rvanka%p_KcolB
!    p_DB1 => rvanka%p_DB1
!    p_DB2 => rvanka%p_DB2
!    p_KldD => rvanka%p_KldD
!    p_KcolD => rvanka%p_KcolD
!    p_DD1 => rvanka%p_DD1
!    p_DD2 => rvanka%p_DD2
!
!    if(associated(rvanka%p_DA12)) then
!      bHaveA12 = .TRUE.
!      p_KldA12 => rvanka%p_KldA12
!      p_KcolA12 => rvanka%p_KcolA12
!      p_DA12 => rvanka%p_DA12
!      p_DA21 => rvanka%p_DA21
!    end if
!
!    if(associated(rvanka%p_DC)) then
!      bHaveC = .TRUE.
!      p_KldC => rvanka%p_KldC
!      p_KcolC => rvanka%p_KcolC
!      p_KdiagC => rvanka%p_KdiagonalC
!      p_DC => rvanka%p_DC
!    end if
!
!    ! Get the multiplication factors
!    Dmult = rvanka%Dmultipliers
!
!    ! Take care of the "soft-deactivation" of the sub-matrices
!    bHaveA12 = bHaveA12 .and. ((Dmult(1,2) .ne. 0.0_DP) .or. &
!                               (Dmult(2,1) .ne. 0.0_DP))
!    bHaveC = bHaveC .and. (Dmult(3,3) .ne. 0.0_DP)
!
!    ! Clear the optional matrices
!    Dc = 0.0_DP
!
!    ! Perform the desired number of iterations
!    do iter = 1, niterations
!
!      ! Loop over all elements in the list
!      do ielidx = 1, size(IelementList)
!
!        ! Get the element number which is to be processed.
!        iel = IelementList(ielidx)
!
!        ! Get all DOFs for this element
!        IdofV(:) = $TODO$
!        IdofP(:) = $TODO$
!
!        ! First of all, fetch the local RHS
!        do i = 1, ndofV
!          Df1(i) = p_DrhsU(IdofV(i))   ! f_u
!          Df2(i) = p_DrhsV(IdofV(i))   ! f_v
!        end do
!        do i = 1, ndofP
!          Dfp(i) = p_DrhsP(IdofP(i))    ! f_p
!        end do
!
!        ! Let us update the local RHS vector by subtracting A*u from it:
!        ! f_u := f_u - A11*u
!        ! f_v := f_v - A22*v
!        do k = 1, ndofV
!          i1 = p_KldA(IdofV(k))
!          i2 = p_KldA(IdofV(k)+1)-1
!          daux1 = 0.0_DP
!          daux2 = 0.0_DP
!          do i = i1, i2
!            j = p_KcolA(i)
!            daux1 = daux1 + p_DA11(i)*p_DvecU(j)
!            daux2 = daux2 + p_DA22(i)*p_DvecV(j)
!          end do
!          Df1(k) = Df1(k) - Dmult(1,1)*daux1
!          Df2(k) = Df2(k) - Dmult(2,2)*daux2
!          ! Get the main diagonal entries
!          j = p_KdiagA(IdofV(k))
!          Da1(k) = Dmult(1,1)*p_DA11(j)
!          Da2(k) = Dmult(2,2)*p_DA22(j)
!        end do
!
!        if(bHaveA12) then
!          ! f_u := f_u - A12*v
!          ! f_v := f_v - A21*u
!          do k = 1, ndofV
!            i1 = p_KldA12(IdofV(k))
!            i2 = p_KldA12(IdofV(k)+1)-1
!            daux1 = 0.0_DP
!            daux2 = 0.0_DP
!            do i = i1, i2
!              j = p_KcolA12(i)
!              daux1 = daux1 + p_DA12(i)*p_DvecV(j)
!              daux2 = daux2 + p_DA21(i)*p_DvecU(j)
!            end do
!            Df1(k) = Df1(k) - Dmult(1,2)*daux1
!            Df2(k) = Df2(k) - Dmult(2,1)*daux2
!          end do
!        end if
!
!        ! Now we also need to subtract B*p from our RHS, and by the same time,
!        ! we will build the local B matrices.
!        ! f_u := f_u - B1*p
!        ! f_v := f_v - B2*p
!        do k = 1, ndofV
!          i1 = p_KldB(IdofV(k))
!          i2 = p_KldB(IdofV(k)+1)-1
!          daux1 = 0.0_DP
!          daux2 = 0.0_DP
!          do i = i1, i2
!            j = p_KcolB(i)
!            daux = p_DvecP(j)
!            daux1 = daux1 + p_DB1(i)*daux
!            daux2 = daux2 + p_DB2(i)*daux
!            do l = 1, ndofP
!              if(j .eq. IdofP(l)) then
!                Db1(k,l) = Dmult(1,3)*p_DB1(i)
!                Db2(k,l) = Dmult(2,3)*p_DB2(i)
!                exit
!              end if
!            end do
!          end do
!          Df1(k) = Df1(k) - Dmult(1,3)*daux1
!          Df2(k) = Df2(k) - Dmult(2,3)*daux2
!        end do
!
!        ! Now we also need to subtract D*u from our RHS, and by the same time,
!        ! we will build the local D matrices.
!        ! f_p := f_p - D1*u - D2*v
!        do k = 1, ndofP
!          i1 = p_KldD(IdofP(k))
!          i2 = p_KldD(IdofP(k)+1)-1
!          daux1 = 0.0_DP
!          daux2 = 0.0_DP
!          do i = i1, i2
!            j = p_KcolD(i)
!            daux1 = daux1 + p_DD1(i)*p_DvecU(j)
!            daux2 = daux2 + p_DD2(i)*p_DvecV(j)
!            do l = 1, ndofV
!              if(j .eq. IdofV(l)) then
!                Dd1(k,l) = Dmult(3,1)*p_DD1(i)
!                Dd2(k,l) = Dmult(3,2)*p_DD2(i)
!                exit
!              end if
!            end do
!          end do
!          Dfp(k) = Dfp(k) - Dmult(3,1)*daux1 - Dmult(3,2)*daux2
!        end do
!
!        if(bHaveC) then
!          ! f_p := f_p - C*p
!          do k = 1, ndofP
!            i1 = p_KldC(IdofP(k))
!            i2 = p_KldC(IdofP(k)+1)-1
!            daux1 = 0.0_DP
!            do i = i1, i2
!              daux1 = daux1 + p_DC(i)*p_DvecP(p_KcolC(i))
!            end do
!            Dfp(k) = Dfp(k) - Dmult(3,3)*daux1
!            ! Get the main diagonal entry of C
!            Dc(k) = Dmult(3,3)*p_DC(p_KdiagC(IdofP(k)))
!          end do
!        end if
!
!        ! Invert A1 and A2
!        do i = 1, ndofV
!          Da1(i) = 1.0_DP / Da1(i)
!          Da2(i) = 1.0_DP / Da2(i)
!        end do
!
!        ! Precalculate D * A^-1
!        do i = 1, ndofV
!          do j = 1, ndofP
!            Dd1(j,i) = Dd1(j,i)*Da1(i)
!            Dd2(j,i) = Dd2(j,i)*Da2(i)
!          end do
!        end do
!
!        ! Calculate Schur-Complement of A
!        ! S := -C + D * A^-1 * B
!        do j = 1, ndofP
!          Ds(j) = -Dc(j)
!          do i = 1, ndofV
!            Ds(j) = Ds(j) + Dd1(j,i)*Db1(i,j) &
!                          + Dd2(j,i)*Db2(i,j)
!          end do
!        end do
!
!        ! Calculate pressure
!        ! p := S^-1 * (D * A^-1 * f_u - f_p)
!        do j = 1, ndofP
!          daux = -Dfp(j)
!          do i = 1, ndofV
!            daux = daux + Dd1(j,i)*Df1(i) &
!                        + Dd2(j,i)*Df2(i)
!          end do
!          Dup(j) = daux / Ds(j)
!        end do
!
!        ! Calculate X- and Y-velocity
!        do i = 1, ndofV
!          do j = 1, ndofP
!            Df1(i) = Df1(i) - Db1(i,j)*Dup(j)
!            Df2(i) = Df2(i) - Db2(i,j)*Dup(j)
!          end do
!          Du1(i) = Da1(i)*Df1(i)
!          Du2(i) = Da2(i)*Df2(i)
!        end do
!
!        ! Incorporate our local solution into the global one.
!        do i = 1, ndofV
!          j = IdofV(i)
!          p_DvecU(j) = p_DvecU(j) + domega * Du1(i)
!          p_DvecV(j) = p_DvecV(j) + domega * Du2(i)
!        end do
!        do i = 1, ndofP
!          j = IdofP(i)
!          p_DvecP(j) = p_DvecP(j) + domega * Dup(i)
!        end do
!
!      end do ! ielidx
!
!    end do ! iter
!
!  end subroutine

  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ! $TEMPLATE$ for "full" Vanka

!<subroutine>
!
!  subroutine vanka_NS2D_$TODO$(rvanka, rsol, rrhs, niterations, domega, &
!                               IelementList)
!
!!<description>
!  ! Performs a desired number of iterations of the Vanka solver for
!  ! 2D Navier-Stokes problem, "full" variant for $TODO$ discretisations.
!!</description>
!
!!<input>
!  ! t_vanka_NavSt2D structure that saves algorithm-specific parameters.
!  type(t_vanka_NavSt2D), intent(in) :: rvanka
!
!  ! The right-hand-side vector of the system
!  type(t_vectorBlock), intent(in) :: rrhs
!
!  ! The number of iterations that are to be performed
!  integer, intent(in) :: niterations
!
!  ! Relaxation parameter.
!  real(DP), intent(in) :: domega
!
!  ! A list of element numbers where Vanka should be applied to.
!  integer, dimension(:), intent(in) :: IelementList
!!</input>
!
!!<inputoutput>
!  ! The iteration vector that is to be updated.
!  type(t_vectorBlock), intent(inout) :: rsol
!!</inputoutput>
!
!!</subroutine>
!
!  ! How many local DOFs do we have?
!  integer, parameter :: ndofV = $TODO$   ! Dofs per velocity
!  integer, parameter :: ndofP = $TODO$   ! Dofs per pressure
!  integer, parameter :: ndof = 2*ndofV+ndofP
!
!  ! Triangulation information
!  ! $TODO$: Declare the necessary variables that need to be accessed from
!  !         the triangulation to perform the DOF-mapping.
!  type(t_triangulation), pointer :: p_rtria
!
!
!  ! Data arrays of vectors
!  real(DP), dimension(:), pointer :: p_DrhsU,p_DrhsV,p_DrhsP,&
!                                     p_DvecU,p_DvecV,p_DvecP
!
!  ! DOF-mapping arrays
!  integer, dimension(ndofV) :: IdofV
!  integer, dimension(ndofP) :: IdofP
!
!  ! Variables for the local system
!  real(DP), dimension(ndof) :: Df
!  real(DP), dimension(ndof,ndof) :: Da
!
!  ! Multiplication factors
!  real(DP), dimension(3,3) :: Dmult
!
!  ! Quick access for the matrix arrays
!  integer, dimension(:), pointer :: p_KldA,p_KldA12,p_KldB,p_KldC,p_KldD,&
!      p_KcolA,p_KcolA12,p_KcolB,p_KcolC,p_KcolD
!  real(DP), dimension(:), pointer :: p_DA11,p_DA12,p_DA21,p_DA22,p_DB1,p_DB2,&
!      p_DD1,p_DD2,p_DC
!
!  ! local variables
!  integer :: i,j,iel,ielidx,i1,i2,k,l,o,p,iter
!  real(DP) :: daux,daux1,daux2
!  logical :: bHaveA12, bHaveC
!
!  ! variables for LAPACK`s DGESV routine
!  integer, dimension(ndof) :: Ipivot
!  integer :: info
!
!
!    ! Get the arrays from the triangulation
!    $TODO$ Get the arrays from the triangulation that are necessary to
!           perform the DOF-mapping
!    p_rtria => rvanka%p_rspatialDiscrV%p_rtriangulation
!
!
!    ! Get the pointers to the vector data
!    call lsyssc_getbase_double(rsol%RvectorBlock(1), p_DvecU)
!    call lsyssc_getbase_double(rsol%RvectorBlock(2), p_DvecV)
!    call lsyssc_getbase_double(rsol%RvectorBlock(3), p_DvecP)
!    call lsyssc_getbase_double(rrhs%RvectorBlock(1), p_DrhsU)
!    call lsyssc_getbase_double(rrhs%RvectorBlock(2), p_DrhsV)
!    call lsyssc_getbase_double(rrhs%RvectorBlock(3), p_DrhsP)
!
!    ! Let us assume we do not have the optional matrices
!    bHaveA12 = .false.
!    bHaveC = .false.
!
!    ! Get the pointers from the vanka structure
!    p_KldA => rvanka%p_KldA
!    p_KcolA => rvanka%p_KcolA
!    p_DA11 => rvanka%p_DA11
!    p_DA22 => rvanka%p_DA22
!    p_KldB => rvanka%p_KldB
!    p_KcolB => rvanka%p_KcolB
!    p_DB1 => rvanka%p_DB1
!    p_DB2 => rvanka%p_DB2
!    p_KldD => rvanka%p_KldD
!    p_KcolD => rvanka%p_KcolD
!    p_DD1 => rvanka%p_DD1
!    p_DD2 => rvanka%p_DD2
!
!    if(associated(rvanka%p_DA12)) then
!      bHaveA12 = .true.
!      p_KldA12 => rvanka%p_KldA12
!      p_KcolA12 => rvanka%p_KcolA12
!      p_DA12 => rvanka%p_DA12
!      p_DA21 => rvanka%p_DA21
!    end if
!
!    if(associated(rvanka%p_DC)) then
!      bHaveC = .true.
!      p_KldC => rvanka%p_KldC
!      p_KcolC => rvanka%p_KcolC
!      p_DC => rvanka%p_DC
!    end if
!
!    ! Get the multiplication factors
!    Dmult = rvanka%Dmultipliers
!
!    ! Take care of the "soft-deactivation" of the sub-matrices
!    bHaveC = bHaveC .and. (Dmult(3,3) .ne. 0.0_DP)
!    bHaveA12 = bHaveA12 .and. &
!      ((Dmult(1,2) .ne. 0.0_DP) .or. (Dmult(2,1) .ne. 0.0_DP))
!
!    ! Perform the desired number of iterations
!    do iter = 1, niterations
!
!      ! Loop over all elements in the list
!      do ielidx = 1, size(IelementList)
!
!        ! Get the element number which is to be processed.
!        iel = IelementList(ielidx)
!
!        ! Get all DOFs for this element
!        IdofV(:) = $TODO$
!        IdofP(:) = $TODO$
!
!        ! First of all, fetch the local RHS
!        do i = 1, ndofV
!          Df(i)       = p_DrhsU(IdofV(i))   ! f_u
!          Df(ndofV+i) = p_DrhsV(IdofV(i))   ! f_v
!        end do
!        o = 2*ndofV
!        do i = 1, ndofP
!          Df(o+i) = p_DrhsP(IdofP(i))       ! f_p
!        end do
!
!        ! Clear the local system matrix
!        Da = 0.0_DP
!
!        ! Let us update the local RHS vector by subtracting A*u from it:
!        ! f_u := f_u - A11*u
!        ! f_v := f_v - A22*v
!        p = ndofV
!        do k = 1, ndofV
!          i1 = p_KldA(IdofV(k))
!          i2 = p_KldA(IdofV(k)+1)-1
!          daux1 = 0.0_DP
!          daux2 = 0.0_DP
!          do i = i1, i2
!            j = p_KcolA(i)
!            daux1 = daux1 + p_DA11(i)*p_DvecU(j)
!            daux2 = daux2 + p_DA22(i)*p_DvecV(j)
!            do l = 1, ndofV
!              if(j .eq. IdofV(l)) then
!                Da(  k,  l) = Dmult(1,1)*p_DA11(i)
!                Da(p+k,p+l) = Dmult(2,2)*p_DA22(i)
!                exit
!              end if
!            end do
!          end do
!          Df(      k) = Df(      k) - Dmult(1,1)*daux1
!          Df(ndofV+k) = Df(ndofV+k) - Dmult(2,2)*daux2
!        end do
!
!        if(bHaveA12) then
!          ! f_u := f_u - A12*v
!          ! f_v := f_v - A21*u
!          do k = 1, ndofV
!            i1 = p_KldA12(IdofV(k))
!            i2 = p_KldA12(IdofV(k)+1)-1
!            daux1 = 0.0_DP
!            daux2 = 0.0_DP
!            do i = i1, i2
!              j = p_KcolA12(i)
!              daux1 = daux1 + p_DA12(i)*p_DvecV(j)
!              daux2 = daux2 + p_DA21(i)*p_DvecU(j)
!              do l = 1, ndofV
!                if(j .eq. IdofV(l)) then
!                  Da(  k,p+l) = Dmult(1,2)*p_DA12(i)
!                  Da(p+k,  l) = Dmult(2,1)*p_DA21(i)
!                  exit
!                end if
!              end do
!            end do
!            Df(      k) = Df(      k) - Dmult(1,2)*daux1
!            Df(ndofV+k) = Df(ndofV+k) - Dmult(2,1)*daux2
!          end do
!        end if ! bHaveA12
!
!        ! Now we also need to subtract B*p from our RHS, and by the same time,
!        ! we will build the local B matrices.
!        ! f_u := f_u - B1*p
!        ! f_v := f_v - B2*p
!        o = ndofV
!        p = 2*ndofV
!        do k = 1, ndofV
!          i1 = p_KldB(IdofV(k))
!          i2 = p_KldB(IdofV(k)+1)-1
!          daux1 = 0.0_DP
!          daux2 = 0.0_DP
!          do i = i1, i2
!            j = p_KcolB(i)
!            daux = p_DvecP(j)
!            daux1 = daux1 + p_DB1(i)*daux
!            daux2 = daux2 + p_DB2(i)*daux
!            do l = 1, ndofP
!              if(j .eq. IdofP(l)) then
!                Da(  k,p+l) = Dmult(1,3)*p_DB1(i)
!                Da(o+k,p+l) = Dmult(2,3)*p_DB2(i)
!                exit
!              end if
!            end do
!          end do
!          Df(      k) = Df(      k) - Dmult(1,3)*daux1
!          Df(ndofV+k) = Df(ndofV+k) - Dmult(2,3)*daux2
!        end do
!
!        ! Now we also need to subtract D*u from our RHS, and by the same time,
!        ! we will build the local D matrices.
!        ! f_p := f_p - D1*u - D2*v
!        o = ndofV
!        p = 2*ndofV
!        do k = 1, ndofP
!          i1 = p_KldD(IdofP(k))
!          i2 = p_KldD(IdofP(k)+1)-1
!          daux1 = 0.0_DP
!          daux2 = 0.0_DP
!          do i = i1, i2
!            j = p_KcolD(i)
!            daux1 = daux1 + p_DD1(i)*p_DvecU(j)
!            daux2 = daux2 + p_DD2(i)*p_DvecV(j)
!            do l = 1, ndofV
!              if(j .eq. IdofV(l)) then
!                Da(p+k,  l) = Dmult(3,1)*p_DD1(i)
!                Da(p+k,o+l) = Dmult(3,2)*p_DD2(i)
!                exit
!              end if
!            end do
!          end do
!          Df(p+k) = Df(p+k) - Dmult(3,1)*daux1 - Dmult(3,2)*daux2
!        end do
!
!        if(bHaveC) then
!          ! f_p := f_p - C*p
!          o = 2*ndofV
!          do k = 1, ndofP
!            i1 = p_KldC(IdofP(k))
!            i2 = p_KldC(IdofP(k)+1)-1
!            do i = i1, i2
!              j = p_KcolC(i)
!              Df(o+k) = Df(o+k) - Dmult(3,3)*p_DC(i)*p_DvecP(j)
!              do l = 1, ndofP
!                if(j .eq. IdofP(l)) then
!                  Da(o+k,o+l) = Dmult(3,3)*p_DC(i)
!                  exit
!                end if
!              end do
!            end do
!          end do
!        end if ! bHaveC
!
!        ! Solve the local system
!        call DGESV(ndof,1,Da,ndof,Ipivot,Df,ndof,info)
!
!        ! Did DGESV fail?
!        if(info .ne. 0) cycle
!
!        ! Incorporate our local solution into the global one.
!        do i = 1, ndofV
!          j = IdofV(i)
!          p_DvecU(j) = p_DvecU(j) + domega * Df(i)
!          p_DvecV(j) = p_DvecV(j) + domega * Df(ndofV+i)
!        end do
!        o = 2*ndofV
!        do i = 1, ndofP
!          j = IdofP(i)
!          p_DvecP(j) = p_DvecP(j) + domega * Df(o+i)
!        end do
!
!      end do ! ielidx
!
!    end do ! iter
!
!  end subroutine
! // unhide from automatic documentation parser -->

end module
