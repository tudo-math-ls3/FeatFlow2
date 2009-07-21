!##############################################################################
!# ****************************************************************************
!# <name> vanka_optcontrol </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the Vanka driver for optimal control problems.
!# The routines have the usual Vanka interface and can be used as usual.
!#
!# Routines in this module:
!#
!# 1.) vanka_init_NavStOptC2D / vanka_solve_NavStOptC2D / vanka_done_NavStOptC2D
!#     -> Initialise / solve / clean up the Vanka variant for 2D
!#        optimal control problems of the Navier-Stokes equation.
!#        Applies for all FEM spaces with discontinuous pressure.
!# </purpose>
!##############################################################################

module vanka_optcontrol

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
  
  public :: t_vanka_NavStOptC2D
  public :: vanka_init_NavStOptC2D
  public :: vanka_solve_NavStOptC2D
  public :: vanka_done_NavStOptC2D

!<constants>

!<constantblock description="Vanka type identifiers for the 2D Navier-Stokes \
! Optimal Control Class">

  ! Diagonal-type VANKA
  integer, parameter, public :: VANKATP_NAVSTOPTC2D_DIAG = 0

!</constantblock>

!<types>

!<typeblock>

  ! Element distribution information for the 2D Navier-Stokes 
  ! optimal control Vanka driver. COntains some preallocated
  ! data arrays. Private type.
  type t_vanka_NavStOptC2D_eldist

    integer :: ndofu, ndofp
    real(DP), dimension(:,:), pointer :: p_DS1
    real(DP), dimension(:,:), pointer :: p_DS2
    integer, dimension(:), pointer :: p_Ipiv
    integer, dimension(:,:), pointer :: p_IdofsP
    integer, dimension(:,:), pointer :: p_KentryLocalB
    real(DP), dimension(:,:), pointer :: p_DaInv
    integer, dimension(:), pointer :: p_IdofsU
    real(DP), dimension(:,:), pointer :: p_DdefectU
    real(DP), dimension(:,:), pointer :: p_DdefectP
    
  end type

!</typeblock>

!<typeblock>
  
  ! A structure that saves matrix pointers for the 2D Navier-Stokes 
  ! optimal control Vanka driver.
  type t_vanka_NavStOptC2D
  
    private
    
    ! The Vanka subtype that is to be used
    integer :: csubtype = VANKATP_NAVSTOPTC2D_DIAG
    
    ! Pointer to the column structure of the velocity matrix A11/A22/A44/A55
    integer, dimension(:), pointer :: p_KcolA11 => null()
    
    ! Pointer to the row structure of the velocity matrix A11/A22/A44/A55
    integer, dimension(:), pointer :: p_KldA11 => null()

    ! Pointer to the diagonal structure of the velocity matrix A11/A22/A44/A55
    integer, dimension(:), pointer :: p_KdiagA11 => null()
    
    ! Pointer to the matrix entries of the velocity matrix A11
    real(DP), dimension(:), pointer :: p_DA11 => null()

    ! Pointer to the matrix entries of the velocity matrix A22
    real(DP), dimension(:), pointer :: p_DA22 => null()

    ! Pointer to the matrix entries of the velocity matrix A44
    real(DP), dimension(:), pointer :: p_DA44 => null()

    ! Pointer to the matrix entries of the velocity matrix A55
    real(DP), dimension(:), pointer :: p_DA55 => null()

    ! Pointer to the matrix entries of the velocity matrix A14
    real(DP), dimension(:), pointer :: p_DA14 => null()

    ! Pointer to the matrix entries of the velocity matrix A25
    real(DP), dimension(:), pointer :: p_DA25 => null()

    ! Pointer to the matrix entries of the velocity matrix A15
    real(DP), dimension(:), pointer :: p_DA15 => null()

    ! Pointer to the matrix entries of the velocity matrix A24
    real(DP), dimension(:), pointer :: p_DA24 => null()

    ! Pointer to the matrix entries of the velocity matrix A41
    real(DP), dimension(:), pointer :: p_DA41 => null()

    ! Pointer to the matrix entries of the velocity matrix A52
    real(DP), dimension(:), pointer :: p_DA52 => null()

    ! Pointer to the matrix entries of the velocity matrix A42
    real(DP), dimension(:), pointer :: p_DA42 => null()

    ! Pointer to the matrix entries of the velocity matrix A51
    real(DP), dimension(:), pointer :: p_DA51 => null()

    ! Pointer to the column structure of the matrix A12/A21/A45/A54/A24/A15/A42/A51
    integer, dimension(:), pointer :: p_KcolA12 => null()
    
    ! Pointer to the row structure of the matrix A12/A21/A45/A54/A24/A15/A42/A51
    integer, dimension(:), pointer :: p_KldA12 => null()

    ! Pointer to the matrix entries of the velocity matrix A12
    real(DP), dimension(:), pointer :: p_DA12 => null()

    ! Pointer to the matrix entries of the velocity matrix A21
    real(DP), dimension(:), pointer :: p_DA21 => null()

    ! Pointer to the matrix entries of the velocity matrix A45
    real(DP), dimension(:), pointer :: p_DA45 => null()

    ! Pointer to the matrix entries of the velocity matrix A54
    real(DP), dimension(:), pointer :: p_DA54 => null()

    ! Pointer to the column structure of the B-matrices.
    integer, dimension(:), pointer :: p_KcolB => null()
    
    ! Pointer to the row structure of the B-matrices
    integer, dimension(:), pointer :: p_KldB => null()
    
    ! Pointer to the entries of the B1-matrix
    real(DP), dimension(:), pointer :: p_DB1 => null()

    ! Pointer to the entries of the B2-matrix
    real(DP), dimension(:), pointer :: p_DB2 => null()

    ! Pointer to the entries of the B4-matrix
    real(DP), dimension(:), pointer :: p_DB4 => null()

    ! Pointer to the entries of the B5-matrix
    real(DP), dimension(:), pointer :: p_DB5 => null()

    ! Pointer to the column structure of the D-matrices.
    integer, dimension(:), pointer :: p_KcolD => null()
    
    ! Pointer to the row structure of the D-matrices
    integer, dimension(:), pointer :: p_KldD => null()

    ! Pointer to the entries of the D1-matrix
    real(DP), dimension(:), pointer :: p_DD1 => null()

    ! Pointer to the entries of the D2-matrix
    real(DP), dimension(:), pointer :: p_DD2 => null()

    ! Pointer to the entries of the D4-matrix
    real(DP), dimension(:), pointer :: p_DD4 => null()

    ! Pointer to the entries of the D5-matrix
    real(DP), dimension(:), pointer :: p_DD5 => null()

    ! Pointer to the column structure of the C-matrix.
    integer, dimension(:), pointer :: p_KcolC => null()
    
    ! Pointer to the row structure of the C-matrix.
    integer, dimension(:), pointer :: p_KldC => null()
    
    ! Pointer to the entries of the C1-matrix
    real(DP), dimension(:), pointer :: p_DC1 => null()

    ! Pointer to the entries of the C2-matrix
    real(DP), dimension(:), pointer :: p_DC2 => null()

    ! Multiplication factors for the submatrices; taken from the system matrix.
    real(DP), dimension(6,6) :: Dmultipliers = 0.0_DP
    
    ! Whether extra matrices exist
    logical :: bhaveA12 = .false.
    logical :: bhaveA45 = .false.
    logical :: bhaveA14 = .false.
    logical :: bhaveA41 = .false.
    logical :: bhaveA15 = .false.
    logical :: bhaveA51 = .false.
    logical :: bhaveC = .false.
    
    ! Spatial discretisation structure for velocity
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscrV => null()
    
    ! Spatial discretisation structure for pressure
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscrP => null()
    
    ! Information about the FE spaces. Pointers to preallocated data arrays.
    type(t_vanka_NavStOptC2D_eldist), dimension(:), pointer :: p_rfemdata
    
  end type
  
!</typeblock>

!</types>

contains

  ! ***************************************************************************

!<subroutine>
  
  subroutine vanka_init_NavStOptC2D (rvanka, rmatrix, csubtype)
  
!<description>
  ! Initialises the Vanka variant for 2D Navier-Stokes problems.
  ! Applies for all FEM spaces with discontinuous pressure.
!</description>

!<input>
  ! The system matrix of the linear system.
  type(t_matrixBlock), intent(in), target :: rmatrix

  ! Desired subtype. One of the VANKATP_NAVST2D_XXXX constants.
  integer, intent(in) :: csubtype  
!</input>

!<output>
  ! The data structure that saves algorithm-specific parameters.
  type(t_vanka_NavStOptC2D), intent(out) :: rvanka
!</output>

!</subroutine>

  integer :: ndofu,ndofp,ielementdist
  integer(I32) :: celemV, celemP
  type(t_blockDiscretisation), pointer :: p_rblockDiscr
  type(t_elementDistribution), pointer :: p_relementDistrV, p_relementDistrP

    ! Matrix must be 3x3.
    if ((rmatrix%nblocksPerCol .ne. 6) .or. (rmatrix%nblocksPerRow .ne. 6)) then
      call output_line ('System matrix is not 6x6.',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init_NavStOptC2D')
      call sys_halt()
    end if
    
    ! Todo: Do more checks
    
    ! The structure of A(1,3) must be identical to A(3,1) and
    ! that of A(2,3) must be identical to A(3,2).
    if ((rmatrix%RmatrixBlock(1,3)%NA .ne. rmatrix%RmatrixBlock(2,3)%NA) .OR. &
        (rmatrix%RmatrixBlock(1,3)%NEQ .ne. rmatrix%RmatrixBlock(2,3)%NEQ)) then
      call output_line ('Structure of B1 and B2 different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init_NavStOptC2D')
      call sys_halt()
    end if

    if ((rmatrix%RmatrixBlock(3,1)%NA .ne. rmatrix%RmatrixBlock(3,2)%NA) .OR. &
        (rmatrix%RmatrixBlock(3,1)%NEQ .ne. rmatrix%RmatrixBlock(3,2)%NEQ)) then
      call output_line ('Structure of D1 and D2 different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init_NavStOptC2D')
      call sys_halt()
    end if
    
    ! Get the block discretisation structure from the matrix.
    p_rblockDiscr => rmatrix%p_rblockDiscrTest
    
    if (.NOT. associated(p_rblockDiscr)) then
      call output_line ('No block discretisation assigned to matrix!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init_NavStOptC2D')
      call sys_halt()
    end if
    
    ! Get the discretisation structure of V and P from the block
    ! discretisation structure.
    rvanka%p_rspatialDiscrV => p_rblockDiscr%RspatialDiscr(1)
    rvanka%p_rspatialDiscrP => p_rblockDiscr%RspatialDiscr(3)
    
    if ((p_rblockDiscr%RspatialDiscr(1)%inumFESpaces .ne. &
         p_rblockDiscr%RspatialDiscr(2)%inumFESpaces) .or. &
        (p_rblockDiscr%RspatialDiscr(1)%inumFESpaces .ne. &
         p_rblockDiscr%RspatialDiscr(4)%inumFESpaces) .or. &
        (p_rblockDiscr%RspatialDiscr(1)%inumFESpaces .ne. &
         p_rblockDiscr%RspatialDiscr(5)%inumFESpaces)) then
      call output_line (&
          'Discretisation structures incompatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init_NavStOptC2D')
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
          'Discretisation structures of velocity and pressure incompatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_init_NavStOptC2D')
      call sys_halt()
    end if

    ! -----
    ! Store the subtype
    rvanka%csubtype = csubtype
  
    ! Get the pointers for the vanka structure
    call lsyssc_getbase_Kld (rmatrix%RmatrixBlock(1,1),rvanka%p_KldA11)
    call lsyssc_getbase_Kcol (rmatrix%RmatrixBlock(1,1),rvanka%p_KcolA11)
    call lsyssc_getbase_Kdiagonal (rmatrix%RmatrixBlock(1,1),rvanka%p_KdiagA11)

    call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,1),rvanka%p_Da11)
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,2),rvanka%p_Da22)

    rvanka%bhaveA12 = lsysbl_isSubmatrixPresent(rmatrix,1,2)
    if (rvanka%bhaveA12) then
      call lsyssc_getbase_Kld (rmatrix%RmatrixBlock(1,2),rvanka%p_KldA12)
      call lsyssc_getbase_Kcol (rmatrix%RmatrixBlock(1,2),rvanka%p_KcolA12)

      call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,2),rvanka%p_Da12)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,1),rvanka%p_Da21)
    end if

    call lsyssc_getbase_Kld (rmatrix%RmatrixBlock(1,3),rvanka%p_KldB)
    call lsyssc_getbase_Kcol (rmatrix%RmatrixBlock(1,3),rvanka%p_KcolB)

    call lsyssc_getbase_Kld (rmatrix%RmatrixBlock(3,1),rvanka%p_KldD)
    call lsyssc_getbase_Kcol (rmatrix%RmatrixBlock(3,1),rvanka%p_KcolD)

    call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,3),rvanka%p_Db1)
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,3),rvanka%p_Db2)

    call lsyssc_getbase_double (rmatrix%RmatrixBlock(3,1),rvanka%p_Dd1)
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(3,2),rvanka%p_Dd2)
    
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(4,4),rvanka%p_Da44)
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(5,5),rvanka%p_Da55)

    rvanka%bhaveA45 = lsysbl_isSubmatrixPresent(rmatrix,4,5)
    if (rvanka%bhaveA45) then
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(4,5),rvanka%p_Da45)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(5,4),rvanka%p_Da54)
    end if
    
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(4,6),rvanka%p_Db4)
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(5,6),rvanka%p_Db5)

    call lsyssc_getbase_double (rmatrix%RmatrixBlock(6,4),rvanka%p_Dd4)
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(6,5),rvanka%p_Dd5)

    rvanka%bhaveC = lsysbl_isSubmatrixPresent(rmatrix,3,3)
    if (rvanka%bhaveC) then
      call lsyssc_getbase_Kld (rmatrix%RmatrixBlock(3,3),rvanka%p_KldC)
      call lsyssc_getbase_Kcol (rmatrix%RmatrixBlock(3,3),rvanka%p_KcolC)

      call lsyssc_getbase_double (rmatrix%RmatrixBlock(3,3),rvanka%p_DC1)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(6,6),rvanka%p_DC2)
    end if
    
    rvanka%bhaveA14 = lsysbl_isSubmatrixPresent(rmatrix,1,4)
    if (rvanka%bhaveA14) then
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,4),rvanka%p_Da14)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,5),rvanka%p_Da25)
    end if

    rvanka%bhaveA15 = lsysbl_isSubmatrixPresent(rmatrix,1,5)
    if (rvanka%bhaveA15) then
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,5),rvanka%p_Da15)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,4),rvanka%p_Da24)
    end if

    rvanka%bhaveA41 = lsysbl_isSubmatrixPresent(rmatrix,4,1)
    if (rvanka%bhaveA41) then
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(4,1),rvanka%p_Da41)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(5,2),rvanka%p_Da52)
    end if

    rvanka%bhaveA51 = lsysbl_isSubmatrixPresent(rmatrix,5,1)
    if (rvanka%bhaveA51) then
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(5,1),rvanka%p_Da51)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(4,2),rvanka%p_Da42)
    end if

    ! Get the multiplication factors of the submatrices.
    rvanka%Dmultipliers(1:6,1:6) = &
        rmatrix%RmatrixBlock(1:6,1:6)%dscaleFactor

    ! Take care of the 'soft-deactivation' of the sub-matrices
    rvanka%bHaveA12 = rvanka%bHaveA12 .and. &
        ((rvanka%Dmultipliers(1,2) .ne. 0.0_DP) .or. &
         (rvanka%Dmultipliers(2,1) .ne. 0.0_DP))
    rvanka%bHaveA45 = rvanka%bHaveA45 .and. &
        ((rvanka%Dmultipliers(4,5) .ne. 0.0_DP) .or. &
         (rvanka%Dmultipliers(5,4) .ne. 0.0_DP))
    rvanka%bHaveC = rvanka%bHaveC .and. &
        (rvanka%Dmultipliers(3,3) .ne. 0.0_DP)
    
    ! Preallocate memory for the FEM data.
    allocate(rvanka%p_rfemdata(p_rblockDiscr%RspatialDiscr(1)%inumFESpaces))
    
    do ielementdist = 1,size(rvanka%p_rfemdata)
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
      
      ! Which element combination do we have now?
      celemV = elem_getPrimaryElement(p_relementDistrV%celement)
      celemP = elem_getPrimaryElement(p_relementDistrP%celement)

      ! #dofs in the FEM-spaces?
      ndofu = elem_igetNDofLoc(celemV)
      ndofp = elem_igetNDofLoc(celemP)
      
      rvanka%p_rfemdata(ielementdist)%ndofu = ndofu
      rvanka%p_rfemdata(ielementdist)%ndofp = ndofp

      ! Allocate an array for the pressure DOF's.
      ! Allocate memory for the correction on each element.
      ! Note that all P-dofs on one element are connected to all the V-dofs on that
      ! element. Therefore, each row in D corresponding to an arbitrary P-dof on an
      ! element will return all the V-dof's on that element!
      allocate (rvanka%p_rfemdata(ielementdist)%p_IdofsP(ndofp,1))
      allocate (rvanka%p_rfemdata(ielementdist)%p_DdefectU(ndofu,4))
      allocate (rvanka%p_rfemdata(ielementdist)%p_DdefectP(ndofp,2))
      allocate (rvanka%p_rfemdata(ielementdist)%p_Ds1(ndofp,ndofp))
      allocate (rvanka%p_rfemdata(ielementdist)%p_Ds2(ndofp,ndofp))
      allocate (rvanka%p_rfemdata(ielementdist)%p_Ipiv(ndofp))
      allocate (rvanka%p_rfemdata(ielementdist)%p_DaInv(4,ndofu))
      allocate (rvanka%p_rfemdata(ielementdist)%p_IdofsU(ndofu))
      allocate (rvanka%p_rfemdata(ielementdist)%p_KentryLocalB(ndofu,ndofp))
      
    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine vanka_done_NavStOptC2D (rvanka)
  
!<description>
  ! Releases memory.
!</description>

!<inputoutput>
  ! The data structure to clean up
  type(t_vanka_NavStOptC2D), intent(inout) :: rvanka
!</inputoutput>

!</subroutine>

  integer :: i
    
    do i = 1,size(rvanka%p_rfemdata)
      ! Deallocate the memory
      deallocate (rvanka%p_rfemdata(i)%p_IdofsP)
      deallocate (rvanka%p_rfemdata(i)%p_DdefectU)
      deallocate (rvanka%p_rfemdata(i)%p_DdefectP)
      deallocate (rvanka%p_rfemdata(i)%p_Ds1)
      deallocate (rvanka%p_rfemdata(i)%p_Ds2)
      deallocate (rvanka%p_rfemdata(i)%p_Ipiv)
      deallocate (rvanka%p_rfemdata(i)%p_DaInv)
      deallocate (rvanka%p_rfemdata(i)%p_IdofsU)
      deallocate (rvanka%p_rfemdata(i)%p_KentryLocalB)
    end do
    
    deallocate (rvanka%p_rfemdata)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine vanka_solve_NavStOptC2D (rvanka, rsol, rrhs, niterations, domega)
  
!<description>
  ! Performs a desired number of iterations of the Vanka solver for
  ! 2D Navier-Stokes problems.
!</description>

!<input>
  ! The Vanka structure that saves algorithm-specific parameters.
  type(t_vanka_NavStOptC2D), intent(in) :: rvanka

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
      ! We take the element list of the X-velocity as 'primary' element list
      ! and assume that it coincides to that of the Y-velocity (and to that
      ! of the pressure).
      call storage_getbase_int (p_relementDistrV%h_IelementList,p_IelementList)
      
      ! Which element combination do we have now?
      celemV = elem_getPrimaryElement(p_relementDistrV%celement)
      celemP = elem_getPrimaryElement(p_relementDistrP%celement)
        
      ! Which VANKA subtype do we have? The diagonal VANKA of the full VANKA?
      select case (rvanka%csubtype)
      case (VANKATP_NAVSTOPTC2D_DIAG)
        ! Call the jacobi-style vanka
        call vanka_NavStOptC2D(rvanka, rsol, rrhs, niterations, &
                               domega, p_IelementList,rvanka%p_rfemdata(ielementdist))

      case default
        call output_line ('Unknown Vanka subtype!',&
            OU_CLASS_ERROR,OU_MODE_STD,'vanka_solve_NavStOptC2D')
        call sys_halt()  
      
      end select
      
    end do
      
  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine vanka_NavStOptC2D (rvanka,rvector,rrhs,niterations,domega,&
      IelementList,rfemdata)
  
!<description>
  ! Applies the Navier-Stokes optimal control diagonal VANKA variant
  ! to a vector for a list of elements and one pair of velocity/pressure
  ! spaces.
!</description>

!<input>
  ! The Vanka structure that saves algorithm-specific parameters.
  type(t_vanka_NavStOptC2D), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in)                    :: domega
  
  ! The number of iterations that are to be performed
  integer, intent(in) :: niterations
  
  ! A list of element numbers where VANKA should be applied to.
  integer, dimension(:), intent(in) :: IelementList

  ! FemData-Structure of the current FEM space.
  type(t_vanka_NavStOptC2D_eldist), intent(in) :: rfemData
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(inout)         :: rvector
!</inputoutput>

!</subroutine>

  ! Multiplication factors
  real(DP), dimension(6,6) :: Dmult
  
  ! Array for the pressure DOF's on the element
  integer, dimension(1) :: IelIdx2
  
  ! Quick access for the matrix arrays
  integer, dimension(:), pointer :: p_KldA11,p_KldA12,p_KldB,p_KldC,p_KldD,&
      p_KcolA11,p_KcolA12,p_KcolB,p_KcolC,p_KcolD,p_KdiagA11
  real(DP), dimension(:), pointer :: p_DA11,p_DA12,p_DA21,p_DA22
  real(DP), dimension(:), pointer :: p_DA14,p_DA15,p_DA24,p_DA25
  real(DP), dimension(:), pointer :: p_DA41,p_DA42,p_DA51,p_DA52
  real(DP), dimension(:), pointer :: p_DA44,p_DA45,p_DA54,p_DA55
  real(DP), dimension(:), pointer :: p_DB1,p_DB2,p_DD1,p_DD2
  real(DP), dimension(:), pointer :: p_DB4,p_DB5,p_DD4,p_DD5
  real(DP), dimension(:), pointer :: p_DC1, p_DC2

  real(DP), dimension(:,:), pointer :: DS1, DS2
  integer, dimension(:), pointer :: Ipiv
  integer, dimension(:,:), pointer :: IdofsP,KentryLocalB
  real(DP), dimension(:,:), pointer :: DaInv
  integer, dimension(:), pointer :: IdofsU
  real(DP), dimension(:,:), pointer :: DdefectU
  real(DP), dimension(:,:), pointer :: DdefectP

  ! Quick access for the vector arrays
  real(DP), dimension(:), pointer :: p_DrhsU1,p_DrhsV1,p_DrhsP1,&
                                     p_DvecU1,p_DvecV1,p_DvecP1
  real(DP), dimension(:), pointer :: p_DrhsU2,p_DrhsV2,p_DrhsP2,&
                                     p_DvecU2,p_DvecV2,p_DvecP2
  
  ! local variables
  logical :: bHaveA12, bHaveA45, bhaveA14,bhaveA15,bhaveA41,bhaveA51,bHaveC
  integer :: idxu,idxp,idxp2,idofp,idofu,i,j,id1,id2,ndofu,ndofp,info,ielidx,iter
  real(DP) :: daux1,daux2,daux4,daux5
  real(DP) :: dp1,dp2
  
    ! Get the pointers to the vector data
    call lsyssc_getbase_double(rvector%RvectorBlock(1), p_DvecU1)
    call lsyssc_getbase_double(rvector%RvectorBlock(2), p_DvecV1)
    call lsyssc_getbase_double(rvector%RvectorBlock(3), p_DvecP1)
    call lsyssc_getbase_double(rrhs%RvectorBlock(1), p_DrhsU1)
    call lsyssc_getbase_double(rrhs%RvectorBlock(2), p_DrhsV1)
    call lsyssc_getbase_double(rrhs%RvectorBlock(3), p_DrhsP1)

    call lsyssc_getbase_double(rvector%RvectorBlock(4), p_DvecU2)
    call lsyssc_getbase_double(rvector%RvectorBlock(5), p_DvecV2)
    call lsyssc_getbase_double(rvector%RvectorBlock(6), p_DvecP2)
    call lsyssc_getbase_double(rrhs%RvectorBlock(4), p_DrhsU2)
    call lsyssc_getbase_double(rrhs%RvectorBlock(5), p_DrhsV2)
    call lsyssc_getbase_double(rrhs%RvectorBlock(6), p_DrhsP2)
    
    ! Get wqhich optional matrices we have.
    bhaveA12 = rvanka%bhaveA12
    bhaveA45 = rvanka%bhaveA45
    bhaveC = rvanka%bhaveC
    bhaveA14 = rvanka%bhaveA14
    bhaveA15 = rvanka%bhaveA15
    bhaveA41 = rvanka%bhaveA41
    bhaveA51 = rvanka%bhaveA51
    
    ! Get the pointers from the Vanka structure for faster access
    p_KldA11   => rvanka%p_KldA11   
    p_KcolA11  => rvanka%p_KcolA11 
    p_KdiagA11 => rvanka%p_KdiagA11
    p_Da11     => rvanka%p_Da11    
    p_Da22     => rvanka%p_Da22    
    p_KldA12   => rvanka%p_KldA12  
    p_KcolA12  => rvanka%p_KcolA12 
    p_Da12     => rvanka%p_Da12    
    p_Da21     => rvanka%p_Da21    
    p_KldB     => rvanka%p_KldB    
    p_KcolB    => rvanka%p_KcolB   
    p_KldD     => rvanka%p_KldD    
    p_KcolD    => rvanka%p_KcolD   
    p_Db1      => rvanka%p_Db1     
    p_Db2      => rvanka%p_Db2     
    p_Dd1      => rvanka%p_Dd1     
    p_Dd2      => rvanka%p_Dd2     
    p_Da44     => rvanka%p_Da44    
    p_Da55     => rvanka%p_Da55    
    p_Da45     => rvanka%p_Da45    
    p_Da54     => rvanka%p_Da54    
    p_Db4      => rvanka%p_Db4     
    p_Db5      => rvanka%p_Db5     
    p_Dd4      => rvanka%p_Dd4     
    p_Dd5      => rvanka%p_Dd5     
    p_KldC     => rvanka%p_KldC    
    p_KcolC    => rvanka%p_KcolC   
    p_DC1      => rvanka%p_DC1     
    p_DC2      => rvanka%p_DC2     
    p_Da14     => rvanka%p_Da14    
    p_Da25     => rvanka%p_Da25    
    p_Da15     => rvanka%p_Da15    
    p_Da24     => rvanka%p_Da24    
    p_Da41     => rvanka%p_Da41    
    p_Da52     => rvanka%p_Da52    
    p_Da51     => rvanka%p_Da51    
    p_Da42     => rvanka%p_Da42

    ! Get the multiplication factors
    Dmult(:,:) = rvanka%Dmultipliers(:,:)

    ! Get pointers and data from the FEM data structure.
    IdofsP => rfemData%p_IdofsP
    DdefectU => rfemData%p_DdefectU
    DdefectP => rfemData%p_DdefectP
    Ds1 => rfemData%p_Ds1
    Ds2 => rfemData%p_Ds2
    Ipiv => rfemData%p_Ipiv
    DaInv => rfemData%p_DaInv
    IdofsU => rfemData%p_IdofsU
    KentryLocalB => rfemData%p_KentryLocalB
    ndofu = rfemData%ndofu
    ndofp = rfemData%ndofp

    ! Perform niterations iterations
    do iter = 1,niterations

      ! Loop through all elements 
      do ielidx = 1,size(IelementList)
      
        ! On the element, get the local DOF's in the pressure space
        IelIdx2(1) = ielidx
        call dof_locGlobMapping_mult(rvanka%p_rspatialDiscrP, IelIdx2, IdofsP)
            
        ! Get A^-1, which is a diagonal matrix in our case.
        ! We can fetch it by going through the the first line of D corresponding to our
        ! element.
        ! Simultaneously get the DOF's in U.
        idofp = IdofsP(1,1)
        do id1 = p_KldD(idofp), p_KldD(idofp+1)-1
          idofu = p_KcolD(id1)
          idxu = id1 - p_KldD(idofp) + 1
          
          ! Save the DOF for future use.
          IdofsU(idxu) = idofu

          ! Get the main diagonal entries of the A-matrices
          i = p_KdiagA11(idofu)
          DaInv(1,idxu) = 1.0_DP/(Dmult(1,1)*p_DA11(i))
          DaInv(2,idxu) = 1.0_DP/(Dmult(2,2)*p_DA22(i))
          DaInv(3,idxu) = 1.0_DP/(Dmult(4,4)*p_DA44(i))
          DaInv(4,idxu) = 1.0_DP/(Dmult(5,5)*p_DA55(i))

          ! Loop through the B matrix and determine the positions in the
          ! matrix of those P-dofs that belong to our current element.
          do id2 = p_KldB(idofu), p_KldB(idofu+1)-1
            do j=1,ndofp
              if (p_KcolB(id2) .eq. IdofsP(j,1)) then
                KentryLocalB(idxu,j) = id2
                exit
              end if
            end do
          end do

        end do
        
        ! Clear the local defect, fetch the local RHS.
        idofP = IdofsP(1,1)
        do id1 = p_KldD(idofp), p_KldD(idofp+1)-1
          idxu = id1-p_KldD(idofp)+1
          DdefectU(idxu,1) = p_DrhsU1(IdofsU(idxu))
          DdefectU(idxu,2) = p_DrhsV1(IdofsU(idxu))
          DdefectU(idxu,3) = p_DrhsU2(IdofsU(idxu))
          DdefectU(idxu,4) = p_DrhsV2(IdofsU(idxu))
        end do

        do idxp = 1, ndofp
          idofp = IdofsP(idxp,1)
          DdefectP(idxp,1) = p_DrhsP1(idofp)
          DdefectP(idxp,2) = p_DrhsP2(idofp)
        end do

        ! Does the C matrix exist? If yes, then update the local RHS:
        ! f_p := f_p - C*p
        if(bHaveC) then

          ! So let's loop all pressure DOFs on the current element
          do idxp = 1, ndofp
          
            idofP = IdofsP(idxp,1)
          
            ! Get the corresponding RHS entry in pressure space
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = p_KldC(idofp), p_KldC(idofp+1)-1
              daux1 = daux1 + p_DC1(i)*p_DvecP1(p_KcolC(i))
              daux2 = daux2 + p_DC2(i)*p_DvecP2(p_KcolC(i))
            end do
            DdefectP(idxp,1) = DdefectP(idxp,1) - Dmult(3,3)*daux1
            DdefectP(idxp,2) = DdefectP(idxp,2) - Dmult(6,6)*daux2
            
          end do
          
        end if
          
        ! Create: f_u = f_u - A u - B p
        do idxu = 1,ndofu
        
          ! The column index gives us the index of a velocity DOF which is
          ! adjacent to the current pressure dof - so get its index.
          idofu = IdofsU(idxu)
          
          ! Subtract A*u from the local RHS
          ! f_u := f_u - A11*u
          ! f_v := f_v - A22*v
          daux1 = 0.0_DP
          daux2 = 0.0_DP
          daux4 = 0.0_DP
          daux5 = 0.0_DP
          do i = p_KldA11(idofu), p_KldA11(idofu+1)-1
            j = p_KcolA11(i)
            daux1 = daux1 + p_DA11(i)*p_DvecU1(j)
            daux2 = daux2 + p_DA22(i)*p_DvecV1(j)
            daux4 = daux4 + p_DA44(i)*p_DvecU2(j)
            daux5 = daux5 + p_DA55(i)*p_DvecV2(j)
          end do
          DdefectU(idxu,1) = DdefectU(idxu,1) - Dmult(1,1)*daux1
          DdefectU(idxu,2) = DdefectU(idxu,2) - Dmult(2,2)*daux2
          DdefectU(idxu,3) = DdefectU(idxu,3) - Dmult(4,4)*daux4
          DdefectU(idxu,4) = DdefectU(idxu,4) - Dmult(5,5)*daux5

          ! Subtract the pressure stuff.
          ! f_u := f_u - B1*p
          ! f_v := f_v - B2*p
          daux1 = 0.0_DP
          daux2 = 0.0_DP
          daux4 = 0.0_DP
          daux5 = 0.0_DP
          do i = p_KldB(idofu), p_KldB(idofu+1)-1
            dp1 = p_DvecP1(p_KcolB(i))
            dp2 = p_DvecP2(p_KcolB(i))
            daux1 = daux1 + p_DB1(i)*dp1
            daux2 = daux2 + p_DB2(i)*dp1
            daux4 = daux4 + p_DB4(i)*dp2
            daux5 = daux5 + p_DB5(i)*dp2
          end do
          DdefectU(idxu,1) = DdefectU(idxu,1) - Dmult(1,3)*daux1
          DdefectU(idxu,2) = DdefectU(idxu,2) - Dmult(2,3)*daux2
          DdefectU(idxu,3) = DdefectU(idxu,3) - Dmult(4,6)*daux4
          DdefectU(idxu,4) = DdefectU(idxu,4) - Dmult(5,6)*daux5

        end do
        
        ! Create the defect in the divergence space.
        do idxp = 1, ndofp
          idofp = IdofsP(idxp,1)
          daux1 = 0.0_DP
          daux2 = 0.0_DP
          daux4 = 0.0_DP
          daux5 = 0.0_DP
          do id1 = p_KldD(idofp), p_KldD(idofp+1)-1
            j = p_KcolD(id1)
            daux1 = daux1 + p_Dd1(id1)*p_DvecU1(j)
            daux2 = daux2 + p_Dd2(id1)*p_DvecV1(j)
            daux4 = daux4 + p_Dd4(id1)*p_DvecU2(j)
            daux5 = daux5 + p_Dd5(id1)*p_DvecV2(j) 
          end do
          DdefectP(idxp,1) = DdefectP(idxp,1) - Dmult(3,1)*daux1 - Dmult(3,2)*daux2
          DdefectP(idxp,2) = DdefectP(idxp,2) - Dmult(6,4)*daux4 - Dmult(6,5)*daux5
        end do
            
        ! Do the A12/A21 matrices exist? If yes, then we will also need to
        ! update the local defect by these matrices.
        if(bHaveA12) then
        
          ! Create the defect
          do idxu = 1,ndofu
            idofu = IdofsU(idxu)
            
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = p_KldA12(idofu), p_KldA12(idofu+1)-1
              j = p_KcolA12(i)
              daux1 = daux1 + p_DA12(i)*p_DvecV1(j)
              daux2 = daux2 + p_DA21(i)*p_DvecU1(j)
            end do
            DdefectU(idxu,1) = DdefectU(idxu,1) - Dmult(1,2)*daux1
            DdefectU(idxu,2) = DdefectU(idxu,2) - Dmult(2,1)*daux2
          end do
        
        end if

        ! Do the A45/A54 matrices exist? If yes, then we will also need to
        ! update the local defect by these matrices.
        if(bHaveA45) then
        
          ! Create the defect
          do idxu = 1,ndofu
            idofu = IdofsU(idxu)
            
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = p_KldA12(idofu), p_KldA12(idofu+1)-1
              j = p_KcolA12(i)
              daux1 = daux1 + p_DA45(i)*p_DvecV2(j)
              daux2 = daux2 + p_DA54(i)*p_DvecU2(j)
            end do
            DdefectU(idxu,3) = DdefectU(idxu,3) - Dmult(4,5)*daux1
            DdefectU(idxu,4) = DdefectU(idxu,4) - Dmult(5,4)*daux2
          end do

        end if

        ! Do the A14/A25 matrices exist? If yes, then we will also need to
        ! update the local defect by these matrices.
        if(bHaveA14) then

          ! Create the defect
          do idxu = 1,ndofu
            idofu = IdofsU(idxu)
            
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = p_KldA11(idofu), p_KldA11(idofu+1)-1
              j = p_KcolA11(i)
              daux1 = daux1 + p_DA14(i)*p_DvecU2(j)
              daux2 = daux2 + p_DA25(i)*p_DvecV2(j)
            end do
            DdefectU(idxu,1) = DdefectU(idxu,1) - Dmult(1,4)*daux1
            DdefectU(idxu,2) = DdefectU(idxu,2) - Dmult(2,5)*daux2
          end do

        end if

        ! Do the A15/A24 matrices exist? If yes, then we will also need to
        ! update the local defect by these matrices.
        if(bHaveA15) then

          ! Create the defect
          do idxu = 1,ndofu
            idofu = IdofsU(idxu)
            
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = p_KldA12(idofu), p_KldA12(idofu+1)-1
              j = p_KcolA12(i)
              daux1 = daux1 + p_DA15(i)*p_DvecV2(j)
              daux2 = daux2 + p_DA24(i)*p_DvecU2(j)
            end do
            DdefectU(idxu,1) = DdefectU(idxu,1) - Dmult(1,5)*daux1
            DdefectU(idxu,2) = DdefectU(idxu,2) - Dmult(2,4)*daux2
          end do

        end if

        ! Do the A41/A52 matrices exist? If yes, then we will also need to
        ! update the local defect by these matrices.
        if(bHaveA41) then

          ! Create the defect
          do idxu = 1,ndofu
            idofu = IdofsU(idxu)
            
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = p_KldA11(idofu), p_KldA11(idofu+1)-1
              j = p_KcolA11(i)
              daux1 = daux1 + p_DA41(i)*p_DvecU1(j)
              daux2 = daux2 + p_DA52(i)*p_DvecV1(j)
            end do
            DdefectU(idxu,3) = DdefectU(idxu,3) - Dmult(4,1)*daux1
            DdefectU(idxu,4) = DdefectU(idxu,4) - Dmult(5,2)*daux2
          end do

        end if
            

        ! Do the A42/A51 matrices exist? If yes, then we will also need to
        ! update the local defect by these matrices.
        if(bHaveA51) then

          ! Create the defect
          do idxu = 1,ndofu
            idofu = IdofsU(idxu)
            
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = p_KldA12(idofu), p_KldA12(idofu+1)-1
              j = p_KcolA12(i)
              daux1 = daux1 + p_DA42(i)*p_DvecV1(j)
              daux2 = daux2 + p_DA51(i)*p_DvecU1(j)
            end do
            DdefectU(idxu,3) = DdefectU(idxu,3) - Dmult(4,2)*daux1
            DdefectU(idxu,4) = DdefectU(idxu,4) - Dmult(5,1)*daux2
          end do

        end if

        ! Now we have the defect "d = f-Bp-Au". 
        !
        ! In the next step, we apply a local preconditioner P^-1 to get an element update:
        !
        !   x  =  x + omega * P^-1 d  
        !      =  x + omega * P^-1 (f_u-Bp-Au , f_p - Du - Cp)^T
        !
        ! For the preconditioner, we choose
        !   P = ( diag(A) B )
        !       ( D       C )
        !
        ! from the local system matrix
        !
        !       ( A B )
        !       ( D C )
        !
        ! The local matrices A,B,C,D are rather small. We apply a Schur complement
        ! decomposition to get an update. For the full local systm matrix, this
        ! has the form:
        !
        !  P^-1 d  = ( A B ) ^-1  ( d1 )
        !            ( D C )      ( d2 )
        !
        !          = ( A^-1 ( d1 - B ( S^-1 ( d2 - DA^-1 d1 ) ) )
        !            (                 S^-1 ( d2 - DA^-1 d1 )   )
        !
        ! where  S = C - D A^-1 B.
        !
        ! In a first step, we set upo the vector v=(d2 - DA^-1 d1)
        ! which is later multiplied to S^-1.
        ! The matrix A here is in our case actually the diagonal
        ! of the original matrix.

        do idxp = 1, ndofp
        
          idofp = IdofsP(idxp,1)
          
          do id1 = p_KldD(idofp), p_KldD(idofp+1)-1
   
            idofu = p_KcolD(id1)
            idxu = id1-p_KldD(idofp)+1

            ! Write v into d2.
            DdefectP(idxp,1) = DdefectP(idxp,1) - p_Dd1(id1)*DaInv(1,idxu)*DdefectU(idxu,1) &
                                                - p_Dd2(id1)*DaInv(2,idxu)*DdefectU(idxu,2)
            DdefectP(idxp,2) = DdefectP(idxp,2) - p_Dd4(id1)*DaInv(3,idxu)*DdefectU(idxu,3) &
                                                - p_Dd5(id1)*DaInv(4,idxu)*DdefectU(idxu,4)
          end do
        end do
        
        ! Now, we have to apply S^-1. That means in a first step, we have to
        ! set up S = C - D A^-1 B. We ignore the C at first as we yet don't know
        ! if it exists.
        do idxp2 = 1,ndofp   ! Columns to compute
          do idxp = 1, ndofp     ! Rows to compute
            ! Compute an entry
            idofp = IdofsP (idxp,1)
            Ds1(idxp,idxp2) = 0.0_DP
            Ds2(idxp,idxp2) = 0.0_DP
            do idxu = 1,ndofu
              id1 = p_KldD(idofp)+idxu-1
              id2 = KentryLocalB(idxu,idxp2)
              Ds1(idxp,idxp2) = Ds1(idxp,idxp2) &
                                -p_Dd1(id1)*p_Db1(id2)*DaInv(1,idxu)&
                                -p_Dd2(id1)*p_Db2(id2)*DaInv(2,idxu) 
              Ds2(idxp,idxp2) = Ds2(idxp,idxp2) &
                                -p_Dd4(id1)*p_Db4(id2)*DaInv(3,idxu)&
                                -p_Dd5(id1)*p_Db5(id2)*DaInv(4,idxu)
            end do
          end do
        end do
        
        ! If we have C, sum it up to S.
        if(bHaveC) then

          ! So let's loop all pressure DOFs on the current element
          do idxp = 1, ndofp
          
            idofp = IdofsP(idxp,1)
          
            do id2 = p_KldC(idofp), p_KldC(idofp+1)-1
              idxp2 = id2-p_KldC(idofp)+1
              Ds1(idxp,idxp2) = Ds1(idxp,idxp2) + p_DC1(id2)
              Ds2(idxp,idxp2) = Ds2(idxp,idxp2) + p_DC2(id2)
            end do
            
          end do
          
        end if
        
        ! Apply S^-1 to d2 to get the update dp=S^-1 ( d2 - DA^-1 d1 )  for p.
        call DGESV(ndofp,1,Ds1,ndofp,Ipiv,DdefectP(:,1),ndofp,info)
        
        ! Did DGESV fail?
        if(info .ne. 0) cycle

        ! Apply S^-1 to d2 to get the update dp=S^-1 ( d2 - DA^-1 d1 )  for p.
        call DGESV(ndofp,1,Ds2,ndofp,Ipiv,DdefectP(:,2),ndofp,info)
        
        ! Did DGESV fail?
        if(info .ne. 0) cycle
        
        ! Get the update for u:
        ! du = A^-1 ( d1 - B dp )
        
        do idxu = 1,ndofu
          idofu = IdofsU (idxu)
          do idxp = 1,ndofp
            id1 = KentryLocalB(idxu,idxp)
            DdefectU(idxu,1) = DdefectU(idxu,1) - p_Db1(id1)*DdefectP(idxp,1)
            DdefectU(idxu,2) = DdefectU(idxu,2) - p_Db2(id1)*DdefectP(idxp,1)
            DdefectU(idxu,3) = DdefectU(idxu,3) - p_Db4(id1)*DdefectP(idxp,2)
            DdefectU(idxu,4) = DdefectU(idxu,4) - p_Db5(id1)*DdefectP(idxp,2)
          end do
          DdefectU(idxu,1) = DaInv(1,idxu)*DdefectU(idxu,1)
          DdefectU(idxu,2) = DaInv(2,idxu)*DdefectU(idxu,2)
          DdefectU(idxu,3) = DaInv(3,idxu)*DdefectU(idxu,3)
          DdefectU(idxu,4) = DaInv(4,idxu)*DdefectU(idxu,4)
        end do
        
        ! Do the update: x_n+1 = x_n + omega*(du,dp)
        do i=1,ndofu
          idofu = IdofsU(i)
          p_DvecU1(idofu) = p_DvecU1(idofu) + domega*DdefectU(i,1)
          p_DvecV1(idofu) = p_DvecV1(idofu) + domega*DdefectU(i,2)
          p_DvecU2(idofu) = p_DvecU2(idofu) + domega*DdefectU(i,3)
          p_DvecV2(idofu) = p_DvecV2(idofu) + domega*DdefectU(i,4)
        end do
              
        do i=1,ndofp
          idofp = IdofsP(i,1)
          p_DvecP1(idofp) = p_DvecP1(idofp) + domega*DdefectP(i,1)
          p_DvecP2(idofp) = p_DvecP2(idofp) + domega*DdefectP(i,2)
        end do

      end do ! ielidx

    end do ! iter

    ! That's it

  end subroutine

end module
