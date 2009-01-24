!##############################################################################
!# ****************************************************************************
!# <name> vanka_navst2d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the vanka driver for 2D Navier-Stokes systems.
!# The most general version of a 2D Navier-System looks as follows:
!#
!#                     / A11 A12 B1 \
!#                     | A21 A22 B2 |
!#                     \ D1  D2  C  /
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
!# Please note that in contrast to the 'old' Navier-Stokes Vanka methods,
!# the D-matrices must NOT be virtually transposed!
!#
!# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!# Solving the local Navier-Stokes System
!# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!# Our local Navier-Stokes system looks as follows:
!# 
!#                   / A11 A12 B1 \   / u1 \   / f_u1 \
!#                   | A21 A22 B2 | * | u2 | = | f_u2 |
!#                   \ D1  D2  C  /   \ p  /   \ f_p  /
!#
!# First, we will rewrite our system:
!#
!# A := / A11 A12 \     B := / B1 \     D := ( D1 D2 )
!#      \ A21 A22 /          \ B2 /
!#
!# u := / u1 \     f_u := / f_u1 \
!#      \ u1 /            \ f_u2 /
!#
!# Now our system looks as follows:
!#
!#                   / A B \ * / u \ = / f_u \
!#                   \ D C /   \ p /   \ f_p /
!#
!# 
!# 1.) Calculate the Schur-Complement of A:
!#
!#                       S := -C + D * A^-1 * B
!#
!# 2.) Calculate pressure:
!#
!#                  p := S^-1 * (D * A^-1 * f_u - f_p)
!#
!# 3.) Calculate velocity:
!#
!#                      u := A^-1 * (f_u - B * p)
!#
!# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!# Pressure-DOF based Vanka information
!# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!# This module offers a so called 'pressure-DOF based Vanka' variant, which,
!# in contrast to the other Vanka variants, is a 'black-box' algorithm which
!# does not need any additional information about the underlying discretisation
!# and therefore *should* work for all 2D Navier-Stokes systems, indepentent
!# of the element combination that is used for the velocity and pressure.
!#
!# Instead of looping over the elements of a mesh and performing the
!# DOF-mapping for each element to determine the set of DOFs which is used
!# for a local system, this Vanka variant loops over the DOFs in pressure
!# space.
!# 
!# </purpose>
!##############################################################################

module vanka_navst2d

use fsystem
use genoutput
use mprimitives
use spatialdiscretisation
use linearsystemscalar
use linearsystemblock
use dofmapping

implicit none

  private
  
  public :: t_vankaPointerNavSt2D
  public :: vanka_initNavierStokes2D
  public :: vanka_doneNavierStokes2D
  public :: vanka_NavierStokes2D
  
  public :: vanka_NS2D_precSPSOR
  public :: vanka_NS2D_precSPSSOR

!<constants>

!<constantblock description="Vanka type identifiers for the 2D Navier-Stokes Class">

  ! Diagonal-type VANKA
  integer, parameter, public :: VANKATP_NAVST2D_DIAG      = 0

  ! 'Full' VANKA
  integer, parameter, public :: VANKATP_NAVST2D_FULL      = 1
  
  ! Pressure-DOF based VANKA
  integer, parameter, public :: VANKATP_NAVST2D_PDOF      = 2
  
  ! SP-SOR
  integer, parameter, public :: VANKATP_NAVST2D_SPSOR     = 3
  
  ! SP-SSOR, symmetric variant of SP-SOR
  ! Warning!
  ! This variant can only be called via vanka_NS2D_precSPSSOR - it can not
  ! be called via the default wrapper routine vanka_NavierStokes2D, as it
  ! needs to work on a defect vector rather than rhs and solution vectors!
  integer, parameter, public :: VANKATP_NAVST2D_SPSSOR    = 4

!</constantblock>

!<constantblock description="Constants for the universal 2D Navier-Stokes Vanka">

  ! Number of elements that are handled simultaneously
  integer, parameter :: VANKA_NAVST2D_NELEMSIM = 10000

!</constantblock>

!<types>

!<typeblock>

  ! A structure holding information for the universal 2D Navier-Stokes Vanka driver.
  type t_vankaNavSt2D_uni_diag
    
    ! Arrays for the doF-mapping
    integer, dimension(:,:), pointer :: IdofV => null()
    integer, dimension(:,:), pointer :: IdofP => null()
    
    ! Local RHS/Solution vectors
    real(DP), dimension(:), pointer :: Du1 => null()
    real(DP), dimension(:), pointer :: Du2 => null()
    real(DP), dimension(:), pointer :: Dup => null()
    real(DP), dimension(:), pointer :: Df1 => null()
    real(DP), dimension(:), pointer :: Df2 => null()
    real(DP), dimension(:), pointer :: Dg => null()
    
    ! Main diagonal entries of A
    real(DP), dimension(:), pointer :: Da1 => null()
    real(DP), dimension(:), pointer :: Da2 => null()
    
    ! B/D matrices
    real(DP), dimension(:,:), pointer :: Db1 => null()
    real(DP), dimension(:,:), pointer :: Db2 => null()
    real(DP), dimension(:,:), pointer :: Dd1 => null()
    real(DP), dimension(:,:), pointer :: Dd2 => null()
    
    ! C matrix entries
    real(DP), dimension(:), pointer :: Dc => null()
    
    ! Schur-complement entries
    real(DP), dimension(:), pointer :: Ds => null()
    
    ! Temporary vectors
    real(DP), dimension(:,:), pointer :: Dt1 => null()
    real(DP), dimension(:,:), pointer :: Dt2 => null()
    
  end type

!</typeblock>

!<typeblock>

  ! A structure holding information for the universal 2D Navier-Stokes Vanka driver.
  type t_vankaNavSt2D_uni_full
    
    ! Arrays for the doF-mapping
    integer, dimension(:,:), pointer :: IdofV => null()
    integer, dimension(:,:), pointer :: IdofP => null()
    
    ! Local RHS/Solution vector
    real(DP), dimension(:), pointer :: Du => null()
    real(DP), dimension(:), pointer :: Df => null()
    
    ! Local matrix
    real(DP), dimension(:,:), pointer :: Da => null()
    
    ! Pivot array for LAPACK
    integer, dimension(:), pointer :: Ipiv => null()
    
  end type

!</typeblock>

!<typeblock>
  
  ! A structure that saves matrix pointers for the 2D Navier-Stokes Vanka driver.
  type t_vankaPointerNavSt2D
  
    private
    
    ! Pointer to the column structure of the velocity matrix A11/A22
    integer, dimension(:), pointer              :: p_KcolA => null()
    
    ! Pointer to the row structure of the velocity matrix A11/A22
    integer, dimension(:), pointer              :: p_KldA => null()
    
    ! Pointer to diagonal entries in the velocity matrix A11/A22
    integer, dimension(:), pointer              :: p_KdiagonalA => null()

    ! Pointer to the matrix entries of the velocity matrix A11
    real(DP), dimension(:), pointer             :: p_DA11 => null()

    ! Pointer to the matrix entries of the velocity matrix A22
    real(DP), dimension(:), pointer             :: p_DA22 => null()

    ! Pointer to the column structure of the matrix A12/A21
    integer, dimension(:), pointer              :: p_KcolA12 => null()
    
    ! Pointer to the row structure of the matrix A12/A21
    integer, dimension(:), pointer              :: p_KldA12 => null()

    ! Pointer to the matrix entries of the velocity matrix A12
    real(DP), dimension(:), pointer             :: p_DA12 => null()

    ! Pointer to the matrix entries of the velocity matrix A21
    real(DP), dimension(:), pointer             :: p_DA21 => null()

    ! Pointer to the column structure of the B-matrices.
    integer, dimension(:), pointer              :: p_KcolB => null()
    
    ! Pointer to the row structure of the B-matrices
    integer, dimension(:), pointer              :: p_KldB => null()
    
    ! Pointer to the entries of the B1-matrix
    real(DP), dimension(:), pointer             :: p_DB1 => null()

    ! Pointer to the entries of the B2-matrix
    real(DP), dimension(:), pointer             :: p_DB2 => null()

    ! Pointer to the column structure of the D-matrices.
    integer, dimension(:), pointer              :: p_KcolD => null()
    
    ! Pointer to the row structure of the D-matrices
    integer, dimension(:), pointer              :: p_KldD => null()

    ! Pointer to the entries of the D1-matrix
    real(DP), dimension(:), pointer             :: p_DD1 => null()

    ! Pointer to the entries of the D2-matrix
    real(DP), dimension(:), pointer             :: p_DD2 => null()

    ! Pointer to the column structure of the C-matrix.
    integer, dimension(:), pointer              :: p_KcolC => null()
    
    ! Pointer to the row structure of the C-matrix.
    integer, dimension(:), pointer              :: p_KldC => null()
    
    ! Pointer to diagonal entries in the C-matrix
    integer, dimension(:), pointer              :: p_KdiagonalC => null()
    
    ! Pointer to the entries of the C-matrix
    real(DP), dimension(:), pointer             :: p_DC => null()

    ! Spatial discretisation structure for velocity
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscrV => null()
    
    ! Spatial discretisation structure for pressure
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscrP => null()
    
    ! Multiplication factors for the submatrices; taken from the system matrix.
    real(DP), dimension(3,3) :: Dmultipliers
    
    ! Structure for the universal diagonal 2D Navier-Stokes Vanka
    type(t_vankaNavSt2D_uni_diag) :: runi_diag
    
    ! Structure for the universal full 2D Navier-Stokes Vanka
    type(t_vankaNavSt2D_uni_full) :: runi_full
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Support for Pressure-DOF based Vanka
    
    ! Total number of DOFs in velocity space component
    integer :: ndofVelo = 0
    
    ! Total number of DOFs in pressure space
    integer :: ndofPres = 0
    
    ! Pointer to an array storing the Schur-Complements
    real(DP), dimension(:), pointer :: p_Dschur => null()

  end type
  
!</typeblock>

!</types>

contains

!<subroutine>
  
  subroutine vanka_initNavierStokes2D (rmatrix,rvanka,csubtype)
  
!<description>
  ! Initialises the VANKA variant for 2D Navier-Stokes problems 
  ! for conformal discretisations.
  ! Checks if the "2D-Navier-Stokes" VANKA variant 
  ! for conformal discretisations can be applied to the system given by rmatrix.
  ! If not, the program is stopped.
!</description>

!<input>
  ! The system matrix of the linear system.
  type(t_matrixBlock), intent(IN), target :: rmatrix

  ! Desired subtype
  integer, intent(IN) :: csubtype  
!</input>

!<inputoutput>
  ! t_vankaPointerNavSt2D structure that saves algorithm-specific parameters.
  type(t_vankaPointerNavSt2D), intent(INOUT) :: rvanka
!</inputoutput>

!</subroutine>

  integer :: i,nmaxdofV,nmaxdofP,ndofV,ndofP,nmaxdof,nmaxel
  integer(I32) :: elV, elP
  type(t_blockDiscretisation), pointer :: p_rblockDiscr
  type(t_elementDistribution), pointer :: p_relementDistrV, p_relementDistrP
  
    ! Matrix must be 3x3.
    if ((rmatrix%nblocksPerCol .ne. 3) .or. (rmatrix%nblocksPerRow .ne. 3)) then
      call output_line ('System matrix is not 3x3.',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_initNavierStokes2D')
      call sys_halt()
    end if
    
    ! TOdo: Do more checks
    
    ! The structure of A(1,3) must be identical to A(3,1) and
    ! that of A(2,3) must be identical to A(3,2).
    if ((rmatrix%RmatrixBlock(1,3)%NA .ne. rmatrix%RmatrixBlock(2,3)%NA) .OR. &
        (rmatrix%RmatrixBlock(1,3)%NEQ .ne. rmatrix%RmatrixBlock(2,3)%NEQ)) then
      call output_line ('Structure of B1 and B2 different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_initNavierStokes2D')
      call sys_halt()
    end if

    if ((rmatrix%RmatrixBlock(3,1)%NA .ne. rmatrix%RmatrixBlock(3,2)%NA) .OR. &
        (rmatrix%RmatrixBlock(3,1)%NEQ .ne. rmatrix%RmatrixBlock(3,2)%NEQ)) then
      call output_line ('Structure of D1 and D2 different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_initNavierStokes2D')
      call sys_halt()
    end if      
  
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
      
      call lsyssc_getbase_double(rmatrix%RmatrixBlock(1,2),&
          rvanka%p_DA12 )
      
      call lsyssc_getbase_double(rmatrix%RmatrixBlock(2,1),&
          rvanka%p_DA21 )
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

    ! Do we use pressure-DOF based Vanka or any of the SP-SOR algorithms?
    if((csubtype .eq. VANKATP_NAVST2D_PDOF) .or. &
       (csubtype .eq. VANKATP_NAVST2D_SPSOR) .or. &
       (csubtype .eq. VANKATP_NAVST2D_SPSSOR)) then
    
      ! Yes, we do. In this case we need to allocate an array for the local
      ! Schur-Complement matrices.
      
      ! Determine the total number of DOFs in pressure space - this is equal
      ! to the number of rows (equations) of the matrix located in block (3,1)
      rvanka%ndofPres = rmatrix%RmatrixBlock(3,1)%NEQ
      
      ! And determine the total number of DOFs in one velocity component -
      ! this is needed by the 'fast' variant.
      rvanka%ndofVelo = rmatrix%RmatrixBlock(1,1)%NEQ
      
      ! Allocate the array.
      allocate(rvanka%p_Dschur(rvanka%ndofPres))

      ! Perform the data-depenent initialisation
      call vanka_initDataNS2D_pdof(rvanka)
      
      ! And we can immediately return here as the rest of the code in this
      ! routine is only used by the other Vanka variants.
      return
    
    end if

    ! Get the block discretisation structure from the matrix.
    p_rblockDiscr => rmatrix%p_rblockDiscrTest
    
    if (.NOT. associated(p_rblockDiscr)) then
      call output_line ('No discretisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_initNavierStokes2D')
      call sys_halt()
    end if
    
    ! Get the discretisation structure of V and P from the block
    ! discretisation structure.
    rvanka%p_rspatialDiscrV => p_rblockDiscr%RspatialDiscr(1)
    rvanka%p_rspatialDiscrP => p_rblockDiscr%RspatialDiscr(3)
    
    if (p_rblockDiscr%RspatialDiscr(1)%inumFESpaces .ne. &
        p_rblockDiscr%RspatialDiscr(2)%inumFESpaces) then
      call output_line (&
          'Discretisation structures of X- and Y-velocity incompatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_initNavierStokes2D')
      call sys_halt()
    end if

    if ((rvanka%p_rspatialDiscrP%inumFESpaces .ne. 1) .AND. &
        (rvanka%p_rspatialDiscrP%inumFESpaces .ne. &
          rvanka%p_rspatialDiscrV%inumFESpaces)) then
      ! Either there must be only one element type for the pressure, or there one
      ! pressure element distribution for every velocity element distribution!
      ! If this is not the case, we cannot determine (at least not in reasonable time)
      ! which element type the pressure represents on a cell!
      call output_line (&
          'Discretisation structures of velocity and pressure incompatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_initNavierStokes2D')
      call sys_halt()
    end if
    
    
    ! Loop through all discretisations
    nmaxdofV = 0
    nmaxdofP = 0
    nmaxel = 0
    do i = 1, rvanka%p_rspatialDiscrV%inumFESpaces
    
      ! Get the corresponding element distributions of V.
      p_relementDistrV => &
          rvanka%p_rspatialDiscrV%RelementDistr(i)
      
      ! Either the same element for P everywhere, or there must be given one
      ! element distribution in the pressure for every velocity element distribution.
      if (rvanka%p_rspatialDiscrP%inumFESpaces .GT. 1) then
        p_relementDistrP => &
            rvanka%p_rspatialDiscrP%RelementDistr(i)
      else
        p_relementDistrP => &
            rvanka%p_rspatialDiscrP%RelementDistr(1)
      end if
      
      ! Get the number of elements
      nmaxel = max(nmaxel,p_relementDistrV%NEL)

      ! Get the elements
      elV = elem_getPrimaryElement(p_relementDistrV%celement)
      elP = elem_getPrimaryElement(p_relementDistrP%celement)
      
      ! Check if we have a specialised variant for this combination - and if yes,
      ! then continue with the next element distribution.
      if((elV .eq. EL_Q1T) .and. (elP .eq. EL_Q0)) cycle
      
      ! Otherwise get the number of local doFs
      ndofV = elem_igetNDofLoc(elV)
      ndofP = elem_igetNDofLoc(elP)
      
      nmaxdofV = MAX(nmaxdofV,ndofV)
      nmaxdofP = MAX(nmaxdofP,ndofP)
      
    end do
    
    ! Don't we use the universal Vanka?
    if((nmaxdofV .eq. 0) .and. (nmaxdofP .eq. 0)) return
    
    nmaxel = min(nmaxel,VANKA_NAVST2D_NELEMSIM)
    
    if(csubtype .eq. VANKATP_NAVST2D_DIAG) then

      ! Okay, allocate the necessary arrays
      allocate(rvanka%runi_diag%IdofV(nmaxdofV,nmaxel))
      allocate(rvanka%runi_diag%IdofP(nmaxdofP,nmaxel))
      allocate(rvanka%runi_diag%Du1(nmaxdofV))
      allocate(rvanka%runi_diag%Du2(nmaxdofV))
      allocate(rvanka%runi_diag%Dup(nmaxdofP))
      allocate(rvanka%runi_diag%Df1(nmaxdofV))
      allocate(rvanka%runi_diag%Df2(nmaxdofV))
      allocate(rvanka%runi_diag%Dg(nmaxdofP))
      allocate(rvanka%runi_diag%Da1(nmaxdofV))
      allocate(rvanka%runi_diag%Da2(nmaxdofV))
      allocate(rvanka%runi_diag%Db1(nmaxdofV,nmaxdofP))
      allocate(rvanka%runi_diag%Db2(nmaxdofV,nmaxdofP))
      allocate(rvanka%runi_diag%Dd1(nmaxdofP,nmaxdofV))
      allocate(rvanka%runi_diag%Dd2(nmaxdofP,nmaxdofV))
      allocate(rvanka%runi_diag%Dc(nmaxdofV))
      allocate(rvanka%runi_diag%Ds(nmaxdofP))
      allocate(rvanka%runi_diag%Dt1(nmaxdofP,nmaxdofV))
      allocate(rvanka%runi_diag%Dt2(nmaxdofP,nmaxdofV))
    
    else
    
      nmaxdof = 2*nmaxdofV + nmaxdofP
    
      ! Allocate the necessary arrays
      allocate(rvanka%runi_full%IdofV(nmaxdofV,nmaxel))
      allocate(rvanka%runi_full%IdofP(nmaxdofP,nmaxel))
      allocate(rvanka%runi_full%Du(nmaxdof))
      allocate(rvanka%runi_full%Df(nmaxdof))
      allocate(rvanka%runi_full%Da(nmaxdof,nmaxdof))
      allocate(rvanka%runi_full%Ipiv(nmaxdof))
    
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine vanka_doneNavierStokes2D (rvanka)
  
!<description>
  ! Releases the VANKA variant for 2D Navier-Stokes problems 
  ! for conformal discretisations.
!</description>

!<inputoutput>
  ! t_vankaPointerNavSt2D structure that saves algorithm-specific parameters.
  type(t_vankaPointerNavSt2D), intent(INOUT) :: rvanka
!</inputoutput>

!</subroutine>

    ! Okay, start releasing everything
    if(associated(rvanka%runi_diag%IdofV)) deallocate(rvanka%runi_diag%IdofV)
    if(associated(rvanka%runi_diag%IdofP)) deallocate(rvanka%runi_diag%IdofP)
    if(associated(rvanka%runi_diag%Du1)) deallocate(rvanka%runi_diag%Du1)
    if(associated(rvanka%runi_diag%Du2)) deallocate(rvanka%runi_diag%Du2)
    if(associated(rvanka%runi_diag%Dup)) deallocate(rvanka%runi_diag%Dup)
    if(associated(rvanka%runi_diag%Df1)) deallocate(rvanka%runi_diag%Df1)
    if(associated(rvanka%runi_diag%Df2)) deallocate(rvanka%runi_diag%Df2)
    if(associated(rvanka%runi_diag%Dg)) deallocate(rvanka%runi_diag%Dg)
    if(associated(rvanka%runi_diag%Da1)) deallocate(rvanka%runi_diag%Da1)
    if(associated(rvanka%runi_diag%Da2)) deallocate(rvanka%runi_diag%Da2)
    if(associated(rvanka%runi_diag%Db1)) deallocate(rvanka%runi_diag%Db1)
    if(associated(rvanka%runi_diag%Db2)) deallocate(rvanka%runi_diag%Db2)
    if(associated(rvanka%runi_diag%Dd1)) deallocate(rvanka%runi_diag%Dd1)
    if(associated(rvanka%runi_diag%Dd2)) deallocate(rvanka%runi_diag%Dd2)
    if(associated(rvanka%runi_diag%Dc)) deallocate(rvanka%runi_diag%Dc)
    if(associated(rvanka%runi_diag%Ds)) deallocate(rvanka%runi_diag%Ds)
    if(associated(rvanka%runi_diag%Dt1)) deallocate(rvanka%runi_diag%Dt1)
    if(associated(rvanka%runi_diag%Dt2)) deallocate(rvanka%runi_diag%Dt2)
    
    ! Okay, start releasing everything
    if(associated(rvanka%runi_full%IdofV)) deallocate(rvanka%runi_full%IdofV)
    if(associated(rvanka%runi_full%IdofP)) deallocate(rvanka%runi_full%IdofP)
    if(associated(rvanka%runi_full%Du)) deallocate(rvanka%runi_full%Du)
    if(associated(rvanka%runi_full%Df)) deallocate(rvanka%runi_full%Df)
    if(associated(rvanka%runi_full%Da)) deallocate(rvanka%runi_full%Da)
    if(associated(rvanka%runi_full%Ipiv)) deallocate(rvanka%runi_full%Ipiv)
    
    if(associated(rvanka%p_Dschur)) deallocate(rvanka%p_Dschur)

  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine vanka_NavierStokes2D (rvanka,rvector,rrhs,domega,csubtype)
  
!<description>
  ! This routine applies the VANKA variant for 2D Navier-Stokes problems
  ! to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
  !
  ! The routine supports arbitrary conformal discretisations, but works
  ! 'specialised'! That means, there must exist a specialised implementation
  ! for every type of problem, which is called by this routine.
  ! So, if the user has a special problem, this routine must tackled to
  ! call the corresponding specialised VANKA variant for that problem!
  ! If a matrix/problem structure is not supported, the routine will stop
  ! the program!
!</description>

!<input>
  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(IN)                    :: domega

  ! The subtype of VANKA that should handle the above problem class.
  ! One of the VANKATP_BOUSS2D_xxxx constants, e.g. VANKATP_BOUSS2D_DIAG.
  integer, intent(IN)                     :: csubtype
  
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(INOUT)         :: rvector

  ! t_vanka structure that saves algorithm-specific parameters.
  type(t_vankaPointerNavSt2D), intent(INOUT) :: rvanka
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ielementdist, ndofV, ndofP
    integer(I32) :: elV, elP
    integer, dimension(:), pointer :: p_IelementList
    type(t_elementDistribution), pointer :: p_relementDistrV
    type(t_elementDistribution), pointer :: p_relementDistrP
    
    ! Do we use the pressure-DOF based Vanka?
    select case(csubtype)
    case (VANKATP_NAVST2D_PDOF)
      
      ! Yes, we do. So call the corresponding routine here.
      ! The pressure-DOF based Vanka is a 'black-box' algorithm which does
      ! not need any information about the underlying discretisations.
      call vanka_NS2D_pdof(rrhs, domega, rvector, rvanka)
      
      ! And return here, as the rest of the code in this routine is only
      ! used by the other Vanka variants.
      return
    
    case(VANKATP_NAVST2D_SPSOR)
      
      ! Yes, we do, but we use the SP-SOR variant.
      call vanka_NS2D_spsor(rrhs, domega, rvector, rvanka)
      
      ! And return here.
      return
    
    case(VANKATP_NAVST2D_SPSSOR)
    
      ! This is bad: The SP-SSOR algorithm can only be applied onto a defect
      ! vector, so a call via this wrapper routine is not allowed!
      call output_line ('SP-SSOR cannot be called via vanka_NavierStokes2D!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'vanka_NavierStokes2D')
      
      ! Stop here
      call sys_halt()
      
    end select
    
    ! Loop through the element distributions of the velocity.
    do ielementdist = 1,rvanka%p_rspatialDiscrV%inumFESpaces
    
      ! Get the corresponding element distributions of U, V and P.
      p_relementDistrV => &
          rvanka%p_rspatialDiscrV%RelementDistr(ielementdist)
      
      ! Either the same element for P everywhere, or there must be given one
      ! element distribution in the pressure for every velocity element distribution.
      if (rvanka%p_rspatialDiscrP%inumFESpaces .GT. 1) then
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
      elV = elem_getPrimaryElement(p_relementDistrV%celement)
      elP = elem_getPrimaryElement(p_relementDistrP%celement)
      if ((elV .eq. EL_Q1T) .AND. (elP .eq. EL_Q0)) then
        ! Q1~/Q1~/Q0 discretisation
        
        ! Which VANKA subtype do we have? The diagonal VANKA of the full VANKA?
        select case (csubtype)
        case (VANKATP_NAVST2D_DIAG)
          ! Call the jacobi-style vanka
          call vanka_NS2D_Q1TQ0_js(rvanka, &
                  rvector, rrhs, domega, p_IelementList)

        case (VANKATP_NAVST2D_FULL)
          if (.NOT. associated(rvanka%p_DA12)) then
            ! Call the block-diagonal vanka
            call vanka_NS2D_Q1TQ0_bd(rvanka, &
                  rvector, rrhs, domega, p_IelementList)
          else
            ! Call the fully coupled vanka
            call vanka_NS2D_Q1TQ0_fc(rvanka, &
                  rvector, rrhs, domega, p_IelementList)
          end if
        
        case default
          call output_line ('Unknown VANKA subtype!',&
              OU_CLASS_ERROR,OU_MODE_STD,'vanka_NavierStokes2D')
          call sys_halt()  
        
        end select
    
      else
        ! Get the number of local doFs
        ndofV = elem_igetNDofLoc(elV)
        ndofP = elem_igetNDofLoc(elP)
        
        ! Choose the universal 2D Navier-Stokes Vanka
        select case(csubtype)
        case (VANKATP_NAVST2D_DIAG)
          call vanka_NS2D_universal_diag(rvanka, rvector, rrhs, domega, &
              p_IelementList, ndofV, ndofP)
        
        case (VANKATP_NAVST2D_FULL)
          call vanka_NS2D_universal_full(rvanka, rvector, rrhs, domega, &
              p_IelementList, ndofV, ndofP)
        end select
        
      end if
      
    end do
      
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_initDataNS2D_pdof(rvanka)

!<description>
  ! This routine performs a matrix-data-dependent initialisation of the
  ! pressure-DOF based Vanka.
!</description>

!<inputoutput>

  ! t_vanka structure that saves algorithm-specific parameters.
  type(t_vankaPointerNavSt2D), intent(INOUT) :: rvanka
  
!</inputoutput>

!</subroutine>

  ! Multiplication factors
  real(DP), dimension(3,3) :: Dmult
  
  ! Quick access for the matrix arrays
  integer, dimension(:), pointer :: p_KldA,p_KldA12,p_KldB,p_KldC,p_KldD,&
      p_KcolA,p_KcolA12,p_KcolB,p_KcolC,p_KcolD,p_KdiagA,p_KdiagC
  real(DP), dimension(:), pointer :: p_DA11,p_DA12,p_DA21,p_DA22,p_DB1,p_DB2,&
      p_DD1,p_DD2,p_DC, p_DS

  ! Temporary local matrices
  !real(DP), dimension(:,:), allocatable :: DA
  !real(DP), dimension(:), allocatable :: DB
  real(DP), dimension(:,:), pointer :: DA
  real(DP), dimension(:), pointer :: DB
  
  ! pivot array for LAPACK
  !integer, dimension(:), allocatable :: Ipivot
  integer, dimension(:), pointer :: Ipivot
  
  ! local variables
  logical :: bHaveA12, bHaveC
  real(DP) :: dc,daux1,daux2
  integer :: idofp,idofu,i,j1,j2,k,id1,id2,nmaxdofV,ndofV,info

    ! Let's assume we do not have the optional matrices
    bHaveA12 = .FALSE.
    bHaveC = .FALSE.
    
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
      bHaveA12 = .TRUE.
      p_KldA12 => rvanka%p_KldA12
      p_KcolA12 => rvanka%p_KcolA12
      p_DA12 => rvanka%p_DA12
      p_DA21 => rvanka%p_DA21
    end if
    
    if(associated(rvanka%p_DC)) then
      bHaveC = .TRUE.
      p_KldC => rvanka%p_KldC
      p_KcolC => rvanka%p_KcolC
      p_KdiagC => rvanka%p_KdiagonalC
      p_DC => rvanka%p_DC
    end if
    
    ! Get the Schur-complement array
    p_DS => rvanka%p_Dschur
    
    ! Get the multiplication factors
    Dmult = rvanka%Dmultipliers
    
    ! Take care of the 'soft-deactivation' of the sub-matrices
    bHaveA12 = bHaveA12 .and. ((Dmult(1,2) .ne. 0.0_DP) .or. &
                               (Dmult(2,1) .ne. 0.0_DP))
    bHaveC = bHaveC .and. (Dmult(3,3) .ne. 0.0_DP)
    
    ! First of all determine the maximum number of velocity DOFs that one
    ! pressure DOF may be adjacent to.
    nmaxdofV = 0
    do idofp = 1, rvanka%ndofPres
      nmaxdofV = max(nmaxdofV, p_KldD(idofp+1)-p_KldD(idofp))
    end do ! idofp
    nmaxdofV = 2*nmaxdofV
    
    ! Okay, let's allocate the temporary data
    allocate(DA(nmaxdofV,nmaxdofV))
    allocate(DB(nmaxdofV))
    allocate(Ipivot(nmaxdofV))
    
    ! Let's loop over the pressure DOFs
    do idofp = 1, rvanka%ndofPres
    
      ! Format the local matrices
      DA = 0.0_DP
      DB = 0.0_DP
      dc = 0.0_DP
      
      ! Format the Schur-complement entry for this pressure DOF
      p_DS(idofP) = 0.0_DP
      
      ! Fetch the number of velocity DOFs adjacent to this pressure DOF
      ndofV = p_KldD(idofp+1) - p_KldD(idofp)
      
      ! If the C matrix exists, grab it's main diagonal entry
      if(bHaveC) then
        dc = Dmult(3,3)*p_DC(p_KdiagC(idofp))
      end if
      
      ! Let's loop over all velocity DOFs which are adjacent to the current
      ! pressure DOF.
      do id1 = p_KldD(idofp), p_KldD(idofp+1)-1
      
        ! Get the index of the velocity DOF
        idofu = p_KcolD(id1)
        
        ! Let's fetch the local A11/A22 matrices
        do i = p_KldA(idofu), p_KldA(idofu+1)-1
        
          ! Get the column index
          k = p_KcolA(i)
          
          ! Let's see if this corresponds to one of our local velocity DOFs
          do id2 = p_KldD(idofp), p_KldD(idofp+1)-1
            if(k .eq. p_KcolD(id2)) then
              ! Okay, incorporate the entries into the local matrix
              j1 = id1 - p_KldD(idofp) + 1
              j2 = id2 - p_KldD(idofp) + 1
              DA(      j1,      j2) = Dmult(1,1)*p_DA11(i)
              DA(ndofV+j1,ndofV+j2) = Dmult(2,2)*p_DA22(i)
              exit
            end if
          end do ! id2
        end do ! i

        ! Do the A12/A21 matrices exist? If yes, then we also need to grab
        ! their local sub-matrices.
        if(bHaveA12) then
          ! Let's fetch the local A12/A21 matrices
          do i = p_KldA12(idofu), p_KldA12(idofu+1)-1
          
            ! Get the column index
            k = p_KcolA12(i)
            
            ! Let's see if this corresponds to one of our local velocity DOFs
            do id2 = p_KldD(idofp), p_KldD(idofp+1)-1
              if(k .eq. p_KcolD(id2)) then
                ! Okay, incorporate the entries into the local matrix
                j1 = id1 - p_KldD(idofp) + 1
                j2 = id2 - p_KldD(idofp) + 1
                DA(      j1,ndofV+j2) = Dmult(1,2)*p_DA12(i)
                DA(ndofV+j1,      j2) = Dmult(2,1)*p_DA21(i)
                exit
              end if
            end do ! id2
          end do ! i
        end if

        ! Let's fetch the local B matrices
        do i = p_KldB(idofu), p_KldB(idofu+1)
          if(p_KcolB(i) .eq. idofP) then
            j1 = id1 - p_KldD(idofp) + 1
            DB(      j1) = Dmult(1,3)*p_DB1(i)
            DB(ndofV+j1) = Dmult(2,3)*p_DB2(i)
            exit
          end if
        end do ! i
      
      end do ! id1
      
      ! Okay, try to solve the local system A*X=B
      call DGESV(2*ndofV,1,DA,nmaxdofV,Ipivot,DB,nmaxdofV,info)
      
      ! Did LAPACK fail? If yes, simply continue with the next pressure DOF.
      if(info .ne. 0) cycle
      
      ! Okay, let's calculate the Schur-complement dc = C - D * A^-1 * B
      daux1 = 0.0_DP
      daux2 = 0.0_DP
      k = 1
      do id1 = p_KldD(idofp), p_KldD(idofp+1)-1
        daux1 = daux1 + p_DD1(id1)*DB(k)
        daux2 = daux2 + p_DD2(id1)*DB(ndofV+k)
        k = k+1
      end do
      dc = dc - Dmult(3,1)*daux1 - Dmult(3,2)*daux2
      
      ! Now if the Schur-complement matrix is regular, we'll store the inverse
      ! in the corresponding entry in p_DS.
      if(dabs(dc) .gt. SYS_EPSREAL) p_DS(idofP) = 1.0_DP / dc
      
    end do ! idofp
    
    ! And release the temporary memory
    deallocate(Ipivot)
    deallocate(DB)
    deallocate(DA)
    
    ! That's it
  
  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine vanka_NS2D_pdof (rrhs,domega,rvector,rvanka)
  
!<description>
  ! This routine applies the pressure-DOF based Vanka algorithm onto a 2D
  ! Navier-Stokes system Ax = b.
  ! In contrast to the other Vanka variants, this algorithm is a 'black-box'
  ! method which does not need any information about the discretisations.
!</description>

!<input>
  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(IN)                    :: domega
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(INOUT)         :: rvector

  ! t_vanka structure that saves algorithm-specific parameters.
  type(t_vankaPointerNavSt2D), intent(INOUT) :: rvanka
!</inputoutput>

!</subroutine>

  ! Multiplication factors
  real(DP), dimension(3,3) :: Dmult
  
  ! Quick access for the matrix arrays
  integer, dimension(:), pointer :: p_KldA,p_KldA12,p_KldB,p_KldC,p_KldD,&
      p_KcolA,p_KcolA12,p_KcolB,p_KcolC,p_KcolD,p_KdiagA,p_KdiagC
  real(DP), dimension(:), pointer :: p_DA11,p_DA12,p_DA21,p_DA22,p_DB1,p_DB2,&
      p_DD1,p_DD2,p_DC, p_DS

  ! Quick access for the vector arrays
  real(DP), dimension(:), pointer :: p_DrhsU,p_DrhsV,p_DrhsP,&
                                     p_DvecU,p_DvecV,p_DvecP
  
  ! local variables
  logical :: bHaveA12, bHaveC
  integer :: idofp,idofu,i,j,id1
  real(DP) :: dt,daux1,daux2,dfu,dfv,dfp,dd1u,dd2v

    ! Get the pointers to the vector data
    call lsyssc_getbase_double(rvector%RvectorBlock(1), p_DvecU)
    call lsyssc_getbase_double(rvector%RvectorBlock(2), p_DvecV)
    call lsyssc_getbase_double(rvector%RvectorBlock(3), p_DvecP)
    call lsyssc_getbase_double(rrhs%RvectorBlock(1), p_DrhsU)
    call lsyssc_getbase_double(rrhs%RvectorBlock(2), p_DrhsV)
    call lsyssc_getbase_double(rrhs%RvectorBlock(3), p_DrhsP)
    
    ! Let's assume we do not have the optional matrices
    bHaveA12 = .FALSE.
    bHaveC = .FALSE.
    
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
      bHaveA12 = .TRUE.
      p_KldA12 => rvanka%p_KldA12
      p_KcolA12 => rvanka%p_KcolA12
      p_DA12 => rvanka%p_DA12
      p_DA21 => rvanka%p_DA21
    end if
    
    if(associated(rvanka%p_DC)) then
      bHaveC = .TRUE.
      p_KldC => rvanka%p_KldC
      p_KcolC => rvanka%p_KcolC
      p_KdiagC => rvanka%p_KdiagonalC
      p_DC => rvanka%p_DC
    end if
    
    ! Get the Schur-complement array
    p_DS => rvanka%p_Dschur
    
    ! Get the multiplication factors
    Dmult = rvanka%Dmultipliers

    ! Take care of the 'soft-deactivation' of the sub-matrices
    bHaveA12 = bHaveA12 .and. ((Dmult(1,2) .ne. 0.0_DP) .or. &
                               (Dmult(2,1) .ne. 0.0_DP))
    bHaveC = bHaveC .and. (Dmult(3,3) .ne. 0.0_DP)

    ! What case do we have here?
    if((.not. bHaveA12) .and. (.not. bHaveC)) then
      ! No optional matrices
      ! So let's loop over all pressure DOFs
      do idofp = 1, rvanka%ndofPres
      
        ! Get the corresponding RHS entry in pressure space
        dfp = p_DrhsP(idofp)
        
        ! Reset auxiliary variables - these will recieve (D * A^-1 * f_u) in
        ! the next loop.
        dd1u = 0.0_DP
        dd2v = 0.0_DP
        
        ! Now let's loop over the entries of row idofp in the D-matrices
        do id1 = p_KldD(idofp), p_KldD(idofp+1)-1
        
          ! The column index gives us the index of a velocity DOF which is
          ! adjacent to the current pressure dof - so get its index.
          idofu = p_KcolD(id1)
          
          ! Fetch local RHS entries for this velocity DOF
          dfu = p_DrhsU(idofu)
          dfv = p_DrhsV(idofu)
          
          ! The first thing we want to do is to perform:
          ! f_u := f_u - B1*p
          ! f_v := f_v - B2*p
          daux1 = 0.0_DP
          daux2 = 0.0_DP
          do i = p_KldB(idofu), p_KldB(idofu+1)-1
            dt = p_DvecP(p_KcolB(i))
            daux1 = daux1 + p_DB1(i)*dt
            daux2 = daux2 + p_DB2(i)*dt
          end do
          dfu = dfu - Dmult(1,3)*daux1
          dfv = dfv - Dmult(2,3)*daux2
          
          ! Now we'll also subtract A*u from the local RHS
          ! f_u := f_u - A11*u
          ! f_v := f_v - A22*v
          daux1 = 0.0_DP
          daux2 = 0.0_DP
          do i = p_KldA(idofu), p_KldA(idofu+1)-1
            j = p_KcolA(i)
            daux1 = daux1 + p_DA11(i)*p_DvecU(j)
            daux2 = daux2 + p_DA22(i)*p_DvecV(j)
          end do
          dfu = dfu - Dmult(1,1)*daux1
          dfv = dfv - Dmult(2,2)*daux2

          ! Now we have built up the local defect for the velocity space, so
          ! we can now solve the local system and by the same time update the
          ! velocity vector components.
          
          ! Get the main diagonal entries of the A-matrices
          i = p_KdiagA(idofu)
          
          ! And update the velocity DOFs:
          ! u := u + omega * (A11)^-1 * f_u
          ! v := v + omega * (A22)^-1 * f_v
          p_DvecU(idofu) = p_DvecU(idofu) + domega*dfu / (Dmult(1,1)*p_DA11(i))
          p_DvecV(idofu) = p_DvecV(idofu) + domega*dfv / (Dmult(2,2)*p_DA22(i))
          
          ! Finally, update the auxiliary variables for the local defect in
          ! pressure space using the new velocity.
          dd1u = dd1u + p_DD1(id1)*p_DvecU(idofu)
          dd2v = dd2v + p_DD2(id1)*p_DvecV(idofu)
        
        end do ! id1
        
        ! Okay, let's calculate the correction entry for the pressure DOF
        ! p := p + omega * S^-1 * (f_p - D * A^-1 * f_u)
        p_DvecP(idofp) = p_DvecP(idofp) + domega * p_DS(idofp) * &
                        (dfp - Dmult(3,1)*dd1u - Dmult(3,2)*dd2v)
      
      end do ! idofp

    else
      ! General case
      ! So let's loop over all pressure DOFs
      do idofp = 1, rvanka%ndofPres
      
        ! Get the corresponding RHS entry in pressure space
        dfp = p_DrhsP(idofp)
        
        ! Does the C matrix exist? If yes, then update the local RHS:
        ! f_p := f_p - C*p
        if(bHaveC) then
          daux1 = 0.0_DP
          do i = p_KldC(idofp), p_KldC(idofp+1)-1
            daux1 = daux1 + p_DC(i)*p_DvecP(p_KcolC(i))
          end do
          dfp = dfp - Dmult(3,3)*daux1
        end if
        
        ! Reset auxiliary variables - these will recieve (D * A^-1 * f_u) in
        ! the next loop.
        dd1u = 0.0_DP
        dd2v = 0.0_DP
        
        ! Now let's loop over the entries of row idofp in the D-matrices
        do id1 = p_KldD(idofp), p_KldD(idofp+1)-1
        
          ! The column index gives us the index of a velocity DOF which is
          ! adjacent to the current pressure dof - so get its index.
          idofu = p_KcolD(id1)
          
          ! Fetch local RHS entries for this velocity DOF
          dfu = p_DrhsU(idofu)
          dfv = p_DrhsV(idofu)
          
          ! The first thing we want to do is to perform:
          ! f_u := f_u - B1*p
          ! f_v := f_v - B2*p
          daux1 = 0.0_DP
          daux2 = 0.0_DP
          do i = p_KldB(idofu), p_KldB(idofu+1)-1
            dt = p_DvecP(p_KcolB(i))
            daux1 = daux1 + p_DB1(i)*dt
            daux2 = daux2 + p_DB2(i)*dt
          end do
          dfu = dfu - Dmult(1,3)*daux1
          dfv = dfv - Dmult(2,3)*daux2
          
          ! Now we'll also subtract A*u from the local RHS
          ! f_u := f_u - A11*u
          ! f_v := f_v - A22*v
          daux1 = 0.0_DP
          daux2 = 0.0_DP
          do i = p_KldA(idofu), p_KldA(idofu+1)-1
            j = p_KcolA(i)
            daux1 = daux1 + p_DA11(i)*p_DvecU(j)
            daux2 = daux2 + p_DA22(i)*p_DvecV(j)
          end do
          dfu = dfu - Dmult(1,1)*daux1
          dfv = dfv - Dmult(2,2)*daux2
          
          ! Do the A12/A21 matrices exist? If yes, then we will also need to
          ! update the local defect by these matrices.
          if(bHaveA12) then
            ! f_u := f_u - A12*v
            ! f_v := f_v - A21*u
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = p_KldA12(idofu), p_KldA12(idofu+1)-1
              j = p_KcolA12(i)
              daux1 = daux1 + p_DA12(i)*p_DvecV(j)
              daux2 = daux2 + p_DA21(i)*p_DvecU(j)
            end do
            dfu = dfu - Dmult(1,2)*daux1
            dfv = dfv - Dmult(2,1)*daux2
          end if

          ! Now we have built up the local defect for the velocity space, so
          ! we can now solve the local system and by the same time update the
          ! velocity vector components.
          
          ! Get the main diagonal entries of the A-matrices
          i = p_KdiagA(idofu)
          
          ! And update the velocity DOFs:
          ! u := u + omega * (A11)^-1 * f_u
          ! v := v + omega * (A22)^-1 * f_v
          p_DvecU(idofu) = p_DvecU(idofu) + domega*dfu / (Dmult(1,1)*p_DA11(i))
          p_DvecV(idofu) = p_DvecV(idofu) + domega*dfv / (Dmult(2,2)*p_DA22(i))
          
          ! Finally, update the auxiliary variables for the local defect in
          ! pressure space using the new velocity.
          dd1u = dd1u + p_DD1(id1)*p_DvecU(idofu)
          dd2v = dd2v + p_DD2(id1)*p_DvecV(idofu)
        
        end do ! id1
        
        ! Okay, let's calculate the correction entry for the pressure DOF
        ! p := p + omega * S^-1 * (f_p - D * A^-1 * f_u)
        p_DvecP(idofp) = p_DvecP(idofp) + domega * p_DS(idofp) * &
                        (dfp - Dmult(3,1)*dd1u - Dmult(3,2)*dd2v)
      
      end do ! idofp

    end if
    
    ! That's it

  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine vanka_NS2D_spsor (rrhs,domega,rvector,rvanka)
  
!<description>
  ! This routine applies the SP-SOR algorithm onto a 2D Navier-Stokes system
  ! Ax = b.
!</description>

!<input>
  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(IN)                    :: domega
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(INOUT)         :: rvector

  ! t_vanka structure that saves algorithm-specific parameters.
  type(t_vankaPointerNavSt2D), intent(INOUT) :: rvanka
!</inputoutput>

!</subroutine>

  ! Multiplication factors
  real(DP), dimension(3,3) :: Dmult
  
  ! Quick access for the matrix arrays
  integer, dimension(:), pointer :: p_KldA,p_KldA12,p_KldB,p_KldC,p_KldD,&
      p_KcolA,p_KcolA12,p_KcolB,p_KcolC,p_KcolD,p_KdiagA,p_KdiagC
  real(DP), dimension(:), pointer :: p_DA11,p_DA12,p_DA21,p_DA22,p_DB1,p_DB2,&
      p_DD1,p_DD2,p_DC, p_DS

  ! Quick access for the vector arrays
  real(DP), dimension(:), pointer :: p_DrhsU,p_DrhsV,p_DrhsP,&
                                     p_DvecU,p_DvecV,p_DvecP
  
  ! local variables
  logical :: bHaveA12, bHaveC
  integer :: idofp,idofu,i,j
  real(DP) :: dt,daux1,daux2,dfu,dfv,dfp

    ! Get the pointers to the vector data
    call lsyssc_getbase_double(rvector%RvectorBlock(1), p_DvecU)
    call lsyssc_getbase_double(rvector%RvectorBlock(2), p_DvecV)
    call lsyssc_getbase_double(rvector%RvectorBlock(3), p_DvecP)
    call lsyssc_getbase_double(rrhs%RvectorBlock(1), p_DrhsU)
    call lsyssc_getbase_double(rrhs%RvectorBlock(2), p_DrhsV)
    call lsyssc_getbase_double(rrhs%RvectorBlock(3), p_DrhsP)
    
    ! Let's assume we do not have the optional matrices
    bHaveA12 = .FALSE.
    bHaveC = .FALSE.
    
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
      bHaveA12 = .TRUE.
      p_KldA12 => rvanka%p_KldA12
      p_KcolA12 => rvanka%p_KcolA12
      p_DA12 => rvanka%p_DA12
      p_DA21 => rvanka%p_DA21
    end if
    
    if(associated(rvanka%p_DC)) then
      bHaveC = .TRUE.
      p_KldC => rvanka%p_KldC
      p_KcolC => rvanka%p_KcolC
      p_KdiagC => rvanka%p_KdiagonalC
      p_DC => rvanka%p_DC
    end if
    
    ! Get the Schur-complement array
    p_DS => rvanka%p_Dschur
    
    ! Get the multiplication factors
    Dmult = rvanka%Dmultipliers

    ! Take care of the 'soft-deactivation' of the sub-matrices
    bHaveA12 = bHaveA12 .and. ((Dmult(1,2) .ne. 0.0_DP) .or. &
                               (Dmult(2,1) .ne. 0.0_DP))
    bHaveC = bHaveC .and. (Dmult(3,3) .ne. 0.0_DP)
   
    ! Step 1: Update velocity
    if(.not. bHaveA12) then
    
      ! Let's loop over the velocity DOFs
      do idofu = 1, rvanka%ndofVelo
      
        ! Get the corresponding RHS entries in velocity space
        dfu = p_DrhsU(idofu)
        dfv = p_DrhsV(idofu)
        
        ! Calculate A(i,.) * u(.)
        daux1 = 0.0_DP
        daux2 = 0.0_DP
        do i = p_KldA(idofu), p_KldA(idofu+1)-1
          j = p_KcolA(i)
          daux1 = daux1 + p_DA11(i)*p_DvecU(j)
          daux2 = daux2 + p_DA22(i)*p_DvecV(j)
        end do ! i
        dfu = dfu - Dmult(1,1)*daux1
        dfv = dfv - Dmult(2,2)*daux2
        
        ! Calculate B(i,.) * p(.)
        daux1 = 0.0_DP
        daux2 = 0.0_DP
        do i = p_KldB(idofu), p_KldB(idofu+1)-1
          dt = p_DvecP(p_KcolB(i))
          daux1 = daux1 + p_DB1(i)*dt
          daux2 = daux2 + p_DB2(i)*dt
        end do ! i
        dfu = dfu - Dmult(1,3)*daux1
        dfv = dfv - Dmult(2,3)*daux2
        
        ! Divide by A(i,i) and update velocity
        i = p_KdiagA(idofu)
        p_DvecU(idofu) = p_DvecU(idofu) + domega*dfu / (Dmult(1,1)*p_DA11(i))
        p_DvecV(idofu) = p_DvecV(idofu) + domega*dfv / (Dmult(2,2)*p_DA22(i))
      
      end do ! idofu
    
    else
      
      ! Let's loop over the velocity DOFs - we'll process A11/A12 now
      do idofu = 1, rvanka%ndofVelo
      
        ! Get the corresponding RHS entries in velocity space
        dfu = p_DrhsU(idofu)
        
        ! Calculate A11(i,.) * u(.)
        daux1 = 0.0_DP
        do i = p_KldA(idofu), p_KldA(idofu+1)-1
          daux1 = daux1 + p_DA11(i)*p_DvecU(p_KcolA(i))
        end do ! i
        dfu = dfu - Dmult(1,1)*daux1
        
        ! Calculate A12(i,.) * v(.)
        daux1 = 0.0_DP
        do i = p_KldA12(idofu), p_KldA12(idofu+1)-1
          daux1 = daux1 + p_DA12(i)*p_DvecV(p_KcolA12(i))
        end do ! i
        dfu = dfu - Dmult(1,2)*daux2

        ! Calculate B(i,.) * p(.)
        daux1 = 0.0_DP
        do i = p_KldB(idofu), p_KldB(idofu+1)-1
          daux1 = daux1 + p_DB1(i)*p_DvecP(p_KcolB(i))
        end do ! i
        dfu = dfu - Dmult(1,3)*daux1
        
        ! Divide by A(i,i) and update velocity
        p_DvecU(idofu) = p_DvecU(idofu) + domega*dfu &
                       / (Dmult(1,1)*p_DA11(p_KdiagA(idofu)))
      
      end do ! idofu

      ! And now let's go for A21/A22
      do idofu = 1, rvanka%ndofVelo
      
        ! Get the corresponding RHS entries in velocity space
        dfv = p_DrhsV(idofu)
        
        ! Calculate A21(i,.) * u(.)
        daux1 = 0.0_DP
        do i = p_KldA12(idofu), p_KldA12(idofu+1)-1
          daux1 = daux1 + p_DA21(i)*p_DvecU(p_KcolA12(i))
        end do ! i
        dfv = dfv - Dmult(2,1)*daux1
        
        ! Calculate A22(i,.) * v(.)
        daux1 = 0.0_DP
        do i = p_KldA(idofu), p_KldA(idofu+1)-1
          daux1 = daux1 + p_DA22(i)*p_DvecV(p_KcolA(i))
        end do ! i
        dfv = dfv - Dmult(2,2)*daux1

        ! Calculate B(i,.) * p(.)
        daux1 = 0.0_DP
        do i = p_KldB(idofu), p_KldB(idofu+1)-1
          daux1 = daux1 + p_DB2(i)*p_DvecP(p_KcolB(i))
        end do ! i
        dfv = dfv - Dmult(2,3)*daux1
        
        ! Divide by A(i,i) and update velocity
        p_DvecV(idofu) = p_DvecV(idofu) + domega*dfv &
                       / (Dmult(2,2)*p_DA22(p_KdiagA(idofu)))
      
      end do ! idofu

    end if
    
    ! Step 2: Update pressure
    if(.not. bHaveC) then
    
      ! Let's loop over all pressure DOFs
      do idofp = 1, rvanka%ndofPres
      
        ! Skip this pressure DOF if the Schur complement is zero.
        if(p_DS(idofp) .eq. 0.0_DP) cycle
      
        ! Get the corresponding RHS entry in pressure space
        dfp = p_DrhsP(idofp)
        
        ! Loop through the D-matrices to calculate the local defect
        daux1 = 0.0_DP
        daux2 = 0.0_DP
        do i = p_KldD(idofp), p_KldD(idofp+1)-1
          j = p_KcolD(i)
          daux1 = daux1 + p_DD1(i)*p_DvecU(j)
          daux2 = daux2 + p_DD2(i)*p_DvecV(j)
        end do ! id1
        dfp = dfp - Dmult(3,1)*daux1 - Dmult(3,2)*daux2
        
        ! Update pressure DOF
        p_DvecP(idofp) = p_DvecP(idofp) + domega*dfp*p_DS(idofp)
      
      end do ! idofp
    
    else

      ! Let's loop over all pressure DOFs
      do idofp = 1, rvanka%ndofPres
      
        ! Skip this pressure DOF if the Schur complement is zero.
        if(p_DS(idofp) .eq. 0.0_DP) cycle
      
        ! Get the corresponding RHS entry in pressure space
        dfp = p_DrhsP(idofp)
        
        ! Calculate f_p := f_p - C*p
        daux1 = 0.0_DP
        do i = p_KldC(idofp), p_KldC(idofp+1)-1
          daux1 = daux1 + p_DC(i)*p_DvecP(p_KcolC(i))
        end do
        dfp = dfp - Dmult(3,3)*daux1
        
        ! Loop through the D-matrices to calculate:
        ! f_p := f_p - D1*u - D2*v
        daux1 = 0.0_DP
        daux2 = 0.0_DP
        do i = p_KldD(idofp), p_KldD(idofp+1)-1
          j = p_KcolD(i)
          daux1 = daux1 + p_DD1(i)*p_DvecU(j)
          daux2 = daux2 + p_DD2(i)*p_DvecV(j)
        end do ! id1
        dfp = dfp - Dmult(3,1)*daux1 - Dmult(3,2)*daux2
        
        ! Update pressure DOF
        p_DvecP(idofp) = p_DvecP(idofp) + domega*dfp*p_DS(idofp)
      
      end do ! idofp

    end if
    
    ! That's it

  end subroutine


  ! ***************************************************************************

!<subroutine>
  
  subroutine vanka_NS2D_precSPSOR (rvanka,rdef,domega)
  
!<description>
  ! This routine applies the SP-SOR preconditioner onto a 2D Navier-Stokes
  ! system Ax = b.
  ! 
!</description>

!<input>
  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(IN)                    :: domega
!</input>

!<inputoutput>
  ! On entry, the defect vector that is to be preconditioned.
  ! On exit, the preconditioned defect.
  type(t_vectorBlock), intent(INOUT)         :: rdef

  ! t_vanka structure that saves algorithm-specific parameters.
  type(t_vankaPointerNavSt2D), intent(INOUT) :: rvanka
!</inputoutput>

!</subroutine>

  ! Multiplication factors
  real(DP), dimension(3,3) :: Dmult
  
  ! Quick access for the matrix arrays
  integer, dimension(:), pointer :: p_KldA,p_KldA12,p_KldB,p_KldD,&
      p_KcolA,p_KcolA12,p_KcolB,p_KcolD,p_KdiagA
  real(DP), dimension(:), pointer :: p_DA11,p_DA12,p_DA21,p_DA22,p_DB1,p_DB2,&
      p_DD1,p_DD2,p_DS

  ! Quick access for the vector arrays
  real(DP), dimension(:), pointer :: p_DdefU,p_DdefV,p_DdefP
  
  ! local variables
  logical :: bHaveA21
  integer :: idofp,idofu,i,j,idiag
  real(DP) :: daux1,daux2,dfu,dfv,dfp,dsa1,dsa2

    ! Get the pointers to the vector data
    call lsyssc_getbase_double(rdef%RvectorBlock(1), p_DdefU)
    call lsyssc_getbase_double(rdef%RvectorBlock(2), p_DdefV)
    call lsyssc_getbase_double(rdef%RvectorBlock(3), p_DdefP)
    
    ! Let's assume we do not have the optional matrices
    bHaveA21 = .false.
    
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
      bHaveA21 = .true.
      p_KldA12 => rvanka%p_KldA12
      p_KcolA12 => rvanka%p_KcolA12
      p_DA12 => rvanka%p_DA12
      p_DA21 => rvanka%p_DA21
    end if
    
    ! Get the Schur-complement array
    p_DS => rvanka%p_Dschur
    
    ! Get the multiplication factors
    Dmult = rvanka%Dmultipliers
    
    ! Take care of the 'soft-deactivation' of the sub-matrices
    bHaveA21 = bHaveA21 .and. (Dmult(2,1) .ne. 0.0_DP)
   
    ! Step 1: Update velocity
    if(.not. bHaveA21) then
    
      ! Pre-calculate scaling factors
      dsa1 = 1.0_DP / Dmult(1,1)
      dsa2 = 1.0_DP / Dmult(2,2)
    
      ! Let's loop over the velocity DOFs
      do idofu = 1, rvanka%ndofVelo
      
        ! Get pre-scaled old defect vector entries
        dfu = dsa1*p_DdefU(idofu)
        dfv = dsa2*p_DdefV(idofu)
        
        ! Get index of main diagonal entry
        idiag = p_KdiagA(idofu)

        ! Calculate A(i,.) * u(.) for the lower-triangular part of A11/A22
        do i = p_KldA(idofu), idiag-1
          j = p_KcolA(i)
          dfu = dfu - p_DA11(i)*p_DdefU(j)
          dfv = dfv - p_DA22(i)*p_DdefV(j)
        end do ! i
        
        ! Update velocity
        p_DdefU(idofu) = domega*dfu / p_DA11(idiag)
        p_DdefV(idofu) = domega*dfv / p_DA22(idiag)
      
      end do ! idofu
    
    else
      
      ! Let's loop over the velocity DOFs - we'll process A11 now
      dsa1 = 1.0_DP / Dmult(1,1)
      do idofu = 1, rvanka%ndofVelo
      
        ! Get pre-scaled old defect vector entries
        dfu = dsa1*p_DdefU(idofu)
        
        ! Get index of main diagonal entry
        idiag = p_KdiagA(idofu)

        ! Calculate A(i,.) * u(.) for the lower-triangular part of A11
        do i = p_KldA(idofu), idiag-1
          dfu = dfu - p_DA11(i)*p_DdefU(p_KcolA(i))
        end do ! i
        
        ! Update velocity
        p_DdefU(idofu) = domega*dfu / p_DA11(idiag)
      
      end do ! idofu

      ! And take care of A21/A22 now
      dsa1 = 1.0_DP / Dmult(2,2)
      dsa2 = Dmult(1,2) / Dmult(2,2)
      do idofu = 1, rvanka%ndofVelo
      
        ! Get pre-scaled old defect vector entries
        dfv = dsa1*p_DdefV(idofu)
        
        ! Get index of main diagonal entry
        idiag = p_KdiagA(idofu)

        ! Subtract A21*u from defect
        daux1 = 0.0_DP
        do i = p_KldA12(idofu), p_KldA12(idofu+1)-1
          daux1 = daux1 + p_DA21(i)*p_DdefU(p_KcolA12(i))
        end do ! i
        dfv = dfv - dsa2*daux1

        ! Calculate A(i,.) * u(.) for the lower-triangular part of A22
        do i = p_KldA(idofu), idiag-1
          dfv = dfv - p_DA22(i)*p_DdefV(p_KcolA(i))
        end do ! i
        
        ! Update velocity
        p_DdefV(idofu) = domega*dfv / p_DA22(idiag)
      
      end do ! idofu
      
    end if
    
    ! Step 2: Update pressure
    do idofp = 1, rvanka%ndofPres
    
      ! Skip this pressure DOF if the Schur complement is zero.
      if(p_DS(idofp) .eq. 0.0_DP) cycle
    
      ! Loop through the D-matrices to calculate the local defect
      daux1 = 0.0_DP
      daux2 = 0.0_DP
      do i = p_KldD(idofp), p_KldD(idofp+1)-1
        j = p_KcolD(i)
        daux1 = daux1 + p_DD1(i)*p_DdefU(j)
        daux2 = daux2 + p_DD2(i)*p_DdefV(j)
      end do ! id1
      dfp = p_DdefP(idofp) - Dmult(3,1)*daux1 - Dmult(3,2)*daux2
      
      ! Update pressure DOF
      p_DdefP(idofp) = domega * dfp * p_DS(idofp)
    
    end do ! idofp

    ! That's it

  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine vanka_NS2D_precSPSSOR (rvanka,rdef,domega)
  
!<description>
  ! This routine applies the SP-SSOR preconditioner onto a 2D Navier-Stokes
  ! system Ax = b.
  ! 
!</description>

!<input>
  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(IN)                    :: domega
!</input>

!<inputoutput>
  ! On entry, the defect vector that is to be preconditioned.
  ! On exit, the preconditioned defect.
  type(t_vectorBlock), intent(INOUT)         :: rdef

  ! t_vanka structure that saves algorithm-specific parameters.
  type(t_vankaPointerNavSt2D), intent(INOUT) :: rvanka
!</inputoutput>

!</subroutine>

  ! Multiplication factors
  real(DP), dimension(3,3) :: Dmult
  
  ! Quick access for the matrix arrays
  integer, dimension(:), pointer :: p_KldA,p_KldA12,p_KldB,p_KldD,&
      p_KcolA,p_KcolA12,p_KcolB,p_KcolD,p_KdiagA
  real(DP), dimension(:), pointer :: p_DA11,p_DA12,p_DA21,p_DA22,p_DB1,p_DB2,&
      p_DD1,p_DD2,p_DS

  ! Quick access for the vector arrays
  real(DP), dimension(:), pointer :: p_DdefU,p_DdefV,p_DdefP
  
  ! local variables
  logical :: bHaveA12,bHaveA21
  integer :: idofp,idofu,i,j,idiag
  real(DP) :: dt,daux1,daux2,dfu,dfv,dfp,dsa1,dsa2

    ! Get the pointers to the vector data
    call lsyssc_getbase_double(rdef%RvectorBlock(1), p_DdefU)
    call lsyssc_getbase_double(rdef%RvectorBlock(2), p_DdefV)
    call lsyssc_getbase_double(rdef%RvectorBlock(3), p_DdefP)
    
    ! Let's assume we do not have the optional matrices
    bHaveA12 = .false.
    bHaveA21 = .false.
    
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
      bHaveA21 = .true.
      p_KldA12 => rvanka%p_KldA12
      p_KcolA12 => rvanka%p_KcolA12
      p_DA12 => rvanka%p_DA12
      p_DA21 => rvanka%p_DA21
    end if
    
    ! Get the Schur-complement array
    p_DS => rvanka%p_Dschur
    
    ! Get the multiplication factors
    Dmult = rvanka%Dmultipliers
    
    ! Take care of the 'soft-deactivation' of the sub-matrices
    bHaveA12 = bHaveA12 .and. (Dmult(1,2) .ne. 0.0_DP)
    bHaveA21 = bHaveA21 .and. (Dmult(2,1) .ne. 0.0_DP)
   
    ! Step 1: Update velocity
    if(.not. bHaveA21) then
    
      ! Pre-calculate scaling factors
      dsa1 = 1.0_DP / Dmult(1,1)
      dsa2 = 1.0_DP / Dmult(2,2)
    
      ! Let's loop over the velocity DOFs
      do idofu = 1, rvanka%ndofVelo
      
        ! Get pre-scaled old defect vector entries
        dfu = dsa1*p_DdefU(idofu)
        dfv = dsa2*p_DdefV(idofu)
        
        ! Get index of main diagonal entry
        idiag = p_KdiagA(idofu)

        ! Calculate A(i,.) * u(.) for the lower-triangular part of A11/A22
        do i = p_KldA(idofu), idiag-1
          j = p_KcolA(i)
          dfu = dfu - p_DA11(i)*p_DdefU(j)
          dfv = dfv - p_DA22(i)*p_DdefV(j)
        end do ! i
        
        ! Update velocity
        p_DdefU(idofu) = domega*dfu / p_DA11(idiag)
        p_DdefV(idofu) = domega*dfv / p_DA22(idiag)
      
      end do ! idofu
    
    else
      
      ! Let's loop over the velocity DOFs - we'll process A11 now
      dsa1 = 1.0_DP / Dmult(1,1)
      do idofu = 1, rvanka%ndofVelo
      
        ! Get pre-scaled old defect vector entries
        dfu = dsa1*p_DdefU(idofu)
        
        ! Get index of main diagonal entry
        idiag = p_KdiagA(idofu)

        ! Calculate A(i,.) * u(.) for the lower-triangular part of A11
        do i = p_KldA(idofu), idiag-1
          dfu = dfu - p_DA11(i)*p_DdefU(p_KcolA(i))
        end do ! i
        
        ! Update velocity
        p_DdefU(idofu) = domega*dfu / p_DA11(idiag)
      
      end do ! idofu

      ! And take care of A21/A22 now
      dsa1 = 1.0_DP / Dmult(2,2)
      dsa2 = Dmult(1,2) / Dmult(2,2)
      do idofu = 1, rvanka%ndofVelo
      
        ! Get pre-scaled old defect vector entries
        dfv = dsa1*p_DdefV(idofu)
        
        ! Get index of main diagonal entry
        idiag = p_KdiagA(idofu)

        ! Subtract A21*u from defect
        daux1 = 0.0_DP
        do i = p_KldA12(idofu), p_KldA12(idofu+1)-1
          daux1 = daux1 + p_DA21(i)*p_DdefU(p_KcolA12(i))
        end do ! i
        dfv = dfv - dsa2*daux1

        ! Calculate A(i,.) * u(.) for the lower-triangular part of A22
        do i = p_KldA(idofu), idiag-1
          dfv = dfv - p_DA22(i)*p_DdefV(p_KcolA(i))
        end do ! i
        
        ! Update velocity
        p_DdefV(idofu) = domega*dfv / p_DA22(idiag)
      
      end do ! idofu
            
    end if
    
    ! Step 2: Update pressure
    do idofp = 1, rvanka%ndofPres
    
      ! Skip this pressure DOF if the Schur complement is zero.
      if(p_DS(idofp) .eq. 0.0_DP) cycle
    
      ! Loop through the D-matrices to calculate the local defect
      daux1 = 0.0_DP
      daux2 = 0.0_DP
      do i = p_KldD(idofp), p_KldD(idofp+1)-1
        j = p_KcolD(i)
        daux1 = daux1 + p_DD1(i)*p_DdefU(j)
        daux2 = daux2 + p_DD2(i)*p_DdefV(j)
      end do ! id1
      dfp = p_DdefP(idofp) - Dmult(3,1)*daux1 - Dmult(3,2)*daux2
      
      ! Update pressure DOF
      p_DdefP(idofp) = domega * dfp * p_DS(idofp)
    
    end do ! idofp
    
    ! Step 3: Update velocity
    if(.not. bHaveA12) then
    
      ! Pre-calculate scaling factors
      dsa1 = Dmult(1,3) / Dmult(1,1)
      dsa2 = Dmult(2,3) / Dmult(2,2)
    
      ! Let's loop over the velocity DOFs
      do idofu = rvanka%ndofVelo, 1, -1
      
        ! Calculate B * p
        dfu = 0.0_DP
        dfv = 0.0_DP
        do i = p_KldB(idofu), p_KldB(idofu+1)-1
          dt = p_DdefP(p_KcolB(i))
          dfu = dfu - p_DB1(i)*dt
          dfv = dfv - p_DB2(i)*dt
        end do ! i
        dfu = dsa1*dfu
        dfv = dsa2*dfv
        
        ! Get index of main diagonal entry
        idiag = p_KdiagA(idofu)

        ! Calculate A(i,.) * u(.) for the upper-triangular part of A11/A22
        do i = p_KldA(idofu+1)-1, idiag+1, -1
          j = p_KcolA(i)
          dfu = dfu - p_DA11(i)*p_DdefU(j)
          dfv = dfv - p_DA22(i)*p_DdefV(j)
        end do ! i
        
        ! Update velocity
        p_DdefU(idofu) = p_DdefU(idofu) + domega*dfu / p_DA11(idiag)
        p_DdefV(idofu) = p_DdefV(idofu) + domega*dfv / p_DA22(idiag)
      
      end do ! idofu
    
    else
      
      ! Let's loop over the velocity DOFs - we'll process A22 now
      dsa1 = Dmult(2,3) / Dmult(2,2)
      do idofu = rvanka%ndofVelo, 1, -1
      
        ! Calculate B * p
        dfv = 0.0_DP
        do i = p_KldB(idofu), p_KldB(idofu+1)-1
          dfv = dfv - p_DB2(i)*p_DdefP(p_KcolB(i))
        end do ! i
        dfv = dsa1*dfv
        
        ! Get index of main diagonal entry
        idiag = p_KdiagA(idofu)

        ! Calculate A(i,.) * u(.) for the upper-triangular part of A11/A22
        do i = p_KldA(idofu+1)-1, idiag+1, -1
          dfv = dfv - p_DA22(i)*p_DdefV(p_KcolA(i))
        end do ! i
        
        ! Update velocity
        p_DdefV(idofu) = p_DdefV(idofu) + domega*dfv / p_DA22(idiag)
      
      end do ! idofu
      
      ! Take care of A11/A12 now
      dsa1 = Dmult(1,3) / Dmult(1,1)
      dsa2 = Dmult(1,2) / Dmult(1,1)
      do idofu = rvanka%ndofVelo, 1, -1
      
        ! Calculate B * p
        dfu = 0.0_DP
        do i = p_KldB(idofu), p_KldB(idofu+1)-1
          dfu = dfu - p_DB1(i)*p_DdefP(p_KcolB(i))
        end do ! i
        dfu = dsa1*dfu
        
        ! Calculate A12*v
        daux1 = 0.0_DP
        do i = p_KldA12(idofu), p_KldA12(idofu+1)-1
          daux1 = daux1 + p_DA12(i)*p_DdefV(p_KcolA12(i))
        end do
        dfu = dfu - dsa2*daux1
        
        ! Get index of main diagonal entry
        idiag = p_KdiagA(idofu)

        ! Calculate A(i,.) * u(.) for the upper-triangular part of A11/A22
        do i = p_KldA(idofu+1)-1, idiag+1, -1
          dfu = dfu - p_DA11(i)*p_DdefU(p_KcolA(i))
        end do ! i
        
        ! Update velocity
        p_DdefU(idofu) = p_DdefU(idofu) + domega*dfu / p_DA11(idiag)
      
      end do ! idofu

    end if
    
    ! That's it

  end subroutine

  ! ***************************************************************************

  subroutine vanka_NS2D_Q1TQ0_js(rvanka, rvector, rrhs, domega, IelementList)
  
!<description>
  ! This routine applies the specialised diagonal VANKA algorithm for
  ! 2D Navier-Stokes problems with Q1~/Q0 discretisation to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
!</description>

!<input>
  ! t_vankaPointerNavSt2D structure that saves algorithm-specific parameters.
  type(t_vankaPointerNavSt2D), intent(IN) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(IN)                    :: domega
  
  ! A list of element numbers where VANKA should be applied to.
  integer, dimension(:), intent(IN)       :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(IN)         :: rvector
!</inputoutput>

!</subroutine>

  ! How many local doFs do we have?
  integer, parameter :: ndofV = 4   ! Dofs per velocity
  integer, parameter :: ndofP = 1   ! Dofs per pressure
  integer, parameter :: ndof = 2*ndofV+ndofP
  
  ! doFs
  integer, dimension(ndofV) :: IdofV
  integer, dimension(ndofP) :: IdofP
  
  ! Variables for the local system
  real(DP), dimension(ndofV) :: Da1, Da2, Du1, Du2, Df1, Df2
  real(DP), dimension(ndofV,ndofP) :: Db1, Db2
  real(DP), dimension(ndofP,ndofV) :: Dd1, Dd2
  real(DP), dimension(ndofP) :: Ds, Dc, Dup, Dfp
  
  ! temporary vectors
  real(DP), dimension(ndofP,ndofV) :: Dt1, Dt2

  ! Multiplication factors
  real(DP), dimension(3,3) :: Dmult
  
  ! Quick access for the matrix arrays
  integer, dimension(:), pointer :: p_KldA,p_KldA12,p_KldB,p_KldC,p_KldD,&
      p_KcolA,p_KcolA12,p_KcolB,p_KcolC,p_KcolD,p_KdiagA,p_KdiagC
  real(DP), dimension(:), pointer :: p_DA11,p_DA12,p_DA21,p_DA22,p_DB1,p_DB2,&
      p_DD1,p_DD2,p_DC

  ! Quick acces for the vector arrays
  real(DP), dimension(:), pointer :: p_DrhsU,p_DrhsV,p_DrhsP,&
                                     p_DvecU,p_DvecV,p_DvecP
    
  ! Triangulation information
  integer, dimension(:,:), pointer :: p_IedgesAtElement

  ! local variables
  integer :: i,j,iel,ielidx,i1,i2,k,l
  real(DP) :: daux,daux1,daux2
  logical :: bHaveA12, bHaveC

    ! Get pointers to the  triangulation information
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
        p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    
    ! Get the pointers to the vector data
    call lsyssc_getbase_double(rvector%RvectorBlock(1), p_DvecU)
    call lsyssc_getbase_double(rvector%RvectorBlock(2), p_DvecV)
    call lsyssc_getbase_double(rvector%RvectorBlock(3), p_DvecP)
    call lsyssc_getbase_double(rrhs%RvectorBlock(1), p_DrhsU)
    call lsyssc_getbase_double(rrhs%RvectorBlock(2), p_DrhsV)
    call lsyssc_getbase_double(rrhs%RvectorBlock(3), p_DrhsP)
    
    ! Let's assume we do not have the optional matrices
    bHaveA12 = .FALSE.
    bHaveC = .FALSE.
    
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
      bHaveA12 = .TRUE.
      p_KldA12 => rvanka%p_KldA12
      p_KcolA12 => rvanka%p_KcolA12
      p_DA12 => rvanka%p_DA12
      p_DA21 => rvanka%p_DA21
    end if
    
    if(associated(rvanka%p_DC)) then
      bHaveC = .TRUE.
      p_KldC => rvanka%p_KldC
      p_KcolC => rvanka%p_KcolC
      p_KdiagC => rvanka%p_KdiagonalC
      p_DC => rvanka%p_DC
    end if
    
    ! Get the multiplication factors
    Dmult = rvanka%Dmultipliers

    ! Take care of the 'soft-deactivation' of the sub-matrices
    bHaveA12 = bHaveA12 .and. ((Dmult(1,2) .ne. 0.0_DP) .or. &
                               (Dmult(2,1) .ne. 0.0_DP))
    bHaveC = bHaveC .and. (Dmult(3,3) .ne. 0.0_DP)
    
    ! Clear the optional matrices
    Dc = 0.0_DP
    
    ! Now which of the optional matrices are present?
    if((.not. bHaveA12) .and. (.not. bHaveC)) then
      ! No optional matrices
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
        
        ! Let's update the local RHS vector by subtracting A*u from it:
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
            Dt1(j,i) = Dd1(j,i)*Da1(i)
            Dt2(j,i) = Dd2(j,i)*Da2(i)
          end do
        end do

        ! Calculate Schur-Complement of A
        ! S := -C + D * A^-1 * B 
        do j = 1, ndofP
          Ds(j) = -Dc(j)
          do i = 1, ndofV
            Ds(j) = Ds(j) + Dt1(j,i)*Db1(i,j) &
                          + Dt2(j,i)*Db2(i,j)
          end do
        end do
        
        ! Calculate pressure
        ! p := S^-1 * (D * A^-1 * f_u - f_p)
        do j = 1, ndofP
          daux = -Dfp(j)
          do i = 1, ndofV
            daux = daux + Dt1(j,i)*Df1(i) &
                        + Dt2(j,i)*Df2(i)
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
      
    else
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
        
        ! Let's update the local RHS vector by subtracting A*u from it:
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
            Dt1(j,i) = Dd1(j,i)*Da1(i)
            Dt2(j,i) = Dd2(j,i)*Da2(i)
          end do
        end do

        ! Calculate Schur-Complement of A
        ! S := -C + D * A^-1 * B 
        do j = 1, ndofP
          Ds(j) = -Dc(j)
          do i = 1, ndofV
            Ds(j) = Ds(j) + Dt1(j,i)*Db1(i,j) &
                          + Dt2(j,i)*Db2(i,j)
          end do
        end do
        
        ! Calculate pressure
        ! p := S^-1 * (D * A^-1 * f_u - f_p)
        do j = 1, ndofP
          daux = -Dfp(j)
          do i = 1, ndofV
            daux = daux + Dt1(j,i)*Df1(i) &
                        + Dt2(j,i)*Df2(i)
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
    end if

  end subroutine

  ! ***************************************************************************
  
  subroutine vanka_NS2D_Q1TQ0_bd(rvanka, rvector, rrhs, domega, IelementList)
  
!<description>
  ! This routine applies the specialised diagonal VANKA algorithm for
  ! 2D Navier-Stokes problems with Q1~/Q0 discretisation to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
!</description>

!<input>
  ! t_vankaPointerNavSt2D structure that saves algorithm-specific parameters.
  type(t_vankaPointerNavSt2D), intent(IN) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(IN)                    :: domega
  
  ! A list of element numbers where VANKA should be applied to.
  integer, dimension(:), intent(IN)       :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(IN)         :: rvector
!</inputoutput>

!</subroutine>

  ! How many local doFs do we have?
  integer, parameter :: ndofV = 4   ! Dofs per velocity
  integer, parameter :: ndofP = 1   ! Dofs per pressure
  integer, parameter :: ndof = 2*ndofV+ndofP
  
  ! Triangulation information
  integer, dimension(:,:), pointer :: p_IedgesAtElement
  real(DP), dimension(:), pointer :: p_DrhsU,p_DrhsV,p_DrhsP,&
                                     p_DvecU,p_DvecV,p_DvecP
  
  ! DOFs
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
  integer :: i,j,iel,ielidx,i1,i2,k,l
  real(DP) :: daux,daux1,daux2
  logical :: bHaveC

    ! Get pointers to the  triangulation information
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
        p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    
    ! Get the pointers to the vector data
    call lsyssc_getbase_double(rvector%RvectorBlock(1), p_DvecU)
    call lsyssc_getbase_double(rvector%RvectorBlock(2), p_DvecV)
    call lsyssc_getbase_double(rvector%RvectorBlock(3), p_DvecP)
    call lsyssc_getbase_double(rrhs%RvectorBlock(1), p_DrhsU)
    call lsyssc_getbase_double(rrhs%RvectorBlock(2), p_DrhsV)
    call lsyssc_getbase_double(rrhs%RvectorBlock(3), p_DrhsP)
    
    ! Let's assume we do not have the optional matrices
    bHaveC = .FALSE.
    
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
      bHaveC = .TRUE.
      p_KldC => rvanka%p_KldC
      p_KcolC => rvanka%p_KcolC
      p_KdiagC => rvanka%p_KdiagonalC
      p_DC => rvanka%p_DC
    end if

    ! Get the multiplication factors
    Dmult = rvanka%Dmultipliers

    ! Take care of the 'soft-deactivation' of the sub-matrices
    bHaveC = bHaveC .and. (Dmult(3,3) .ne. 0.0_DP)
    
    ! Clear the optional matrices
    Dc = 0.0_DP
    
    if(bHaveC) then
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
        
        ! Let's update the local RHS vector by subtracting A*u from it:
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

        ! Let's update the local RHS vector by subtracting C*p from it:
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
        call mprim_invert4x4MatrixDirectDble(Da1, Di1)
        call mprim_invert4x4MatrixDirectDble(Da2, Di2)
        
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
      
      end do ! ielidx

    else
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
        
        ! Let's update the local RHS vector by subtracting A*u from it:
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
        call mprim_invert4x4MatrixDirectDble(Da1, Di1)
        call mprim_invert4x4MatrixDirectDble(Da2, Di2)
        
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
      
      end do ! ielidx

    end if

  end subroutine

  ! ***************************************************************************
  
  subroutine vanka_NS2D_Q1TQ0_fc(rvanka, rvector, rrhs, domega, IelementList)
  
!<description>
  ! This routine applies the specialised diagonal VANKA algorithm for
  ! 2D Navier-Stokes problems with Q1~/Q0 discretisation to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
!</description>

!<input>
  ! t_vankaPointerNavSt2D structure that saves algorithm-specific parameters.
  type(t_vankaPointerNavSt2D), intent(IN) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(IN)                    :: domega
  
  ! A list of element numbers where VANKA should be applied to.
  integer, dimension(:), intent(IN)       :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(IN)         :: rvector
!</inputoutput>

!</subroutine>

  ! How many local DOFs do we have?
  integer, parameter :: ndofV = 4   ! Dofs per velocity
  integer, parameter :: ndofP = 1   ! Dofs per pressure
  integer, parameter :: ndof = 2*ndofV+ndofP
  
  ! Triangulation information
  integer, dimension(:,:), pointer :: p_IedgesAtElement
  real(DP), dimension(:), pointer :: p_DrhsU,p_DrhsV,p_DrhsP,&
                                     p_DvecU,p_DvecV,p_DvecP
  
  ! DOFs
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
  integer :: i,j,iel,ielidx,i1,i2,k,l,o,p
  real(DP) :: daux,daux1,daux2
  logical :: bHaveC
  
  ! variables for LAPACK's DGESV routine
  integer, dimension(ndof) :: Ipivot
  integer :: info

    ! Get pointers to the  triangulation information
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
        p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    
    ! Get the pointers to the vector data
    call lsyssc_getbase_double(rvector%RvectorBlock(1), p_DvecU)
    call lsyssc_getbase_double(rvector%RvectorBlock(2), p_DvecV)
    call lsyssc_getbase_double(rvector%RvectorBlock(3), p_DvecP)
    call lsyssc_getbase_double(rrhs%RvectorBlock(1), p_DrhsU)
    call lsyssc_getbase_double(rrhs%RvectorBlock(2), p_DrhsV)
    call lsyssc_getbase_double(rrhs%RvectorBlock(3), p_DrhsP)
    
    ! Let's assume we do not have the optional matrices
    bHaveC = .FALSE.
    
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
      bHaveC = .TRUE.
      p_KldC => rvanka%p_KldC
      p_KcolC => rvanka%p_KcolC
      p_DC => rvanka%p_DC
    end if
    
    ! Get the multiplication factors
    Dmult = rvanka%Dmultipliers

    ! Take care of the 'soft-deactivation' of the sub-matrices
    bHaveC = bHaveC .and. (Dmult(3,3) .ne. 0.0_DP)
    
    if(bHaveC) then
      ! C exists
      do ielidx=1, size(IelementList)
      
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
        
        ! Let's update the local RHS vector by subtracting A*u from it:
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
    
    else
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
        
        ! Let's update the local RHS vector by subtracting A*u from it:
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
    
    end if
    
  end subroutine

  ! ***************************************************************************

  subroutine vanka_NS2D_universal_diag(rvanka, rvector, rrhs, domega, &
                                       IelementList, ndofV, ndofP)
  
!<description>
  ! This routine applies the specialised diagonal VANKA algorithm for
  ! 2D Navier-Stokes problems with arbitrary discretisation to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
!</description>

!<input>
  ! t_vankaPointerNavSt2D structure that saves algorithm-specific parameters.
  type(t_vankaPointerNavSt2D), target, intent(IN) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(IN)                    :: domega
  
  ! A list of element numbers where VANKA should be applied to.
  integer, dimension(:), intent(IN)       :: IelementList
  
  ! The number of local doFs.
  integer, intent(IN)                     :: ndofV, ndofP
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(IN)         :: rvector
!</inputoutput>

!</subroutine>

  type(t_vankaNavSt2D_uni_diag), pointer :: p_runi

  real(DP), dimension(:), pointer :: p_DrhsU,p_DrhsV,p_DrhsP,&
                                     p_DvecU,p_DvecV,p_DvecP
  
  ! Multiplication factors
  real(DP), dimension(3,3) :: Dmult
  
  ! Quick access for the matrix arrays
  integer, dimension(:), pointer :: p_KldA,p_KldA12,p_KldB,p_KldC,p_KldD,&
      p_KcolA,p_KcolA12,p_KcolB,p_KcolC,p_KcolD,p_KdiagA,p_KdiagC
  real(DP), dimension(:), pointer :: p_DA11,p_DA12,p_DA21,p_DA22,p_DB1,p_DB2,&
      p_DD1,p_DD2,p_DC
    
  ! local variables
  integer :: i,j,iel,i1,i2,k,l
  integer :: NEL,NELdone,NELtodo
  real(DP) :: daux, daux1,daux2
  logical :: bHaveA12, bHaveC
  
    ! Get the number of elements to process
    NEL = rvanka%p_rspatialDiscrV%p_rtriangulation%NEL
    
    ! Get the diagonal vanka structure
    p_runi => rvanka%runi_diag
    
    ! Get the pointers to the vector data
    call lsyssc_getbase_double(rvector%RvectorBlock(1), p_DvecU)
    call lsyssc_getbase_double(rvector%RvectorBlock(2), p_DvecV)
    call lsyssc_getbase_double(rvector%RvectorBlock(3), p_DvecP)
    call lsyssc_getbase_double(rrhs%RvectorBlock(1), p_DrhsU)
    call lsyssc_getbase_double(rrhs%RvectorBlock(2), p_DrhsV)
    call lsyssc_getbase_double(rrhs%RvectorBlock(3), p_DrhsP)
    
    ! Let's assume we do not have the optional matrices
    bHaveA12 = .FALSE.
    bHaveC = .FALSE.
    
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
      bHaveA12 = .TRUE.
      p_KldA12 => rvanka%p_KldA12
      p_KcolA12 => rvanka%p_KcolA12
      p_DA12 => rvanka%p_DA12
      p_DA21 => rvanka%p_DA21
    end if
    
    if(associated(rvanka%p_DC)) then
      bHaveC = .TRUE.
      p_KldC => rvanka%p_KldC
      p_KcolC => rvanka%p_KcolC
      p_KdiagC => rvanka%p_KdiagonalC
      p_DC => rvanka%p_DC
    end if
    
    ! Get the multiplication factors
    Dmult = rvanka%Dmultipliers

    ! Take care of the 'soft-deactivation' of the sub-matrices
    bHaveA12 = bHaveA12 .and. ((Dmult(1,2) .ne. 0.0_DP) .or. &
                               (Dmult(2,1) .ne. 0.0_DP))
    bHaveC = bHaveC .and. (Dmult(3,3) .ne. 0.0_DP)
    
    ! Clear the optional matrices
    p_runi%Dc = 0.0_DP
    
    ! Okay, loop through the element sets
    NELdone = 0
    do while(NELdone .lt. NEL)
    
      ! How many elements do we process this time?
      NELtodo = MIN(NEL-NELdone,VANKA_NAVST2D_NELEMSIM)
      
      ! Perform the DOF-mapping for the pressure and
      call dof_locGlobMapping_mult(rvanka%p_rspatialDiscrV,&
          IelementList(NELdone+1:NELdone+NELtodo), p_runi%IdofV)
      call dof_locGlobMapping_mult(rvanka%p_rspatialDiscrP,&
          IelementList(NELdone+1:NELdone+NELtodo), p_runi%IdofP)
    
      ! General case
      do iel=1, NELtodo
      
        ! First of all, fetch the local RHS
        do i = 1, ndofV
          p_runi%Df1(i) = p_DrhsU(p_runi%IdofV(i,iel))   ! f_u
          p_runi%Df2(i) = p_DrhsV(p_runi%IdofV(i,iel))   ! f_v
        end do
        do i = 1, ndofP
          p_runi%Dg(i) = p_DrhsP(p_runi%IdofP(i,iel))    ! f_p
        end do
        
        ! Let's update the local RHS vector by subtracting A*u from it:
        ! f_u := f_u - A11*u
        ! f_v := f_v - A22*v
        do k = 1, ndofV
          i1 = p_KldA(p_runi%IdofV(k,iel))
          i2 = p_KldA(p_runi%IdofV(k,iel)+1)-1
          daux1 = 0.0_DP
          daux2 = 0.0_DP
          do i = i1, i2
            j = p_KcolA(i)
            daux1 = daux1 + p_DA11(i)*p_DvecU(j)
            daux2 = daux2 + p_DA22(i)*p_DvecV(j)
          end do
          p_runi%Df1(k) = p_runi%Df1(k) - Dmult(1,1)*daux1
          p_runi%Df2(k) = p_runi%Df2(k) - Dmult(2,2)*daux2
          ! Get the main diagonal entries
          j = p_KdiagA(p_runi%IdofV(k,iel))
          p_runi%Da1(k) = Dmult(1,1)*p_DA11(j)
          p_runi%Da2(k) = Dmult(2,2)*p_DA22(j)
        end do
        
        ! What about A12/A21?
        if(bHaveA12) then
          ! f_u := f_u - A12*v
          ! f_v := f_v - A21*u
          do k = 1, ndofV
            i1 = p_KldA12(p_runi%IdofV(k,iel))
            i2 = p_KldA12(p_runi%IdofV(k,iel)+1)-1
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = i1, i2
              j = p_KcolA12(i)
              daux1 = daux1 + p_DA12(i)*p_DvecV(j)
              daux2 = daux2 + p_DA21(i)*p_DvecU(j)
            end do
            p_runi%Df1(k) = p_runi%Df1(k) - Dmult(1,2)*daux1
            p_runi%Df2(k) = p_runi%Df2(k) - Dmult(2,1)*daux2
          end do
        end if
        
        ! Now we also need to subtract B*p from our RHS, and by the same time,
        ! we will build the local B matrices.
        ! f_u := f_u - B1*p
        ! f_v := f_v - B2*p
        do k = 1, ndofV
          i1 = p_KldB(p_runi%IdofV(k,iel))
          i2 = p_KldB(p_runi%IdofV(k,iel)+1)-1
          daux1 = 0.0_DP
          daux2 = 0.0_DP
          do i = i1, i2
            j = p_KcolB(i)
            daux = p_DvecP(j)
            daux1 = daux1 + p_DB1(i)*daux
            daux2 = daux2 + p_DB2(i)*daux
            do l = 1, ndofP
              if(j .eq. p_runi%IdofP(l,iel)) then
                p_runi%Db1(k,l) = Dmult(1,3)*p_DB1(i)
                p_runi%Db2(k,l) = Dmult(2,3)*p_DB2(i)
                exit
              end if
            end do
          end do
          p_runi%Df1(k) = p_runi%Df1(k) - Dmult(1,3)*daux1
          p_runi%Df2(k) = p_runi%Df2(k) - Dmult(2,3)*daux2
        end do
        
        ! Now we also need to subtract D*u from our RHS, and by the same time,
        ! we will build the local D matrices.
        ! f_p := f_p - D1*u - D2*v
        do k = 1, ndofP
          i1 = p_KldD(p_runi%IdofP(k,iel))
          i2 = p_KldD(p_runi%IdofP(k,iel)+1)-1
          daux1 = 0.0_DP
          daux2 = 0.0_DP
          do i = i1, i2
            j = p_KcolD(i)
            daux1 = daux1 + p_DD1(i)*p_DvecU(j)
            daux2 = daux2 + p_DD2(i)*p_DvecV(j)
            do l = 1, ndofV
              if(j .eq. p_runi%IdofV(l,iel)) then
                p_runi%Dd1(k,l) = Dmult(3,1)*p_DD1(i)
                p_runi%Dd2(k,l) = Dmult(3,2)*p_DD2(i)
                exit
              end if
            end do
          end do
          p_runi%Dg(k) = p_runi%Dg(k) - Dmult(3,1)*daux1 - Dmult(3,2)*daux2
        end do
        
        ! Do we have a C-matrix?
        if(bHaveC) then
          ! Yes, so update the local RHS
          ! f_p := f_p - C*p
          do k = 1, ndofP
            i1 = p_KldC(p_runi%IdofP(k,iel))
            i2 = p_KldC(p_runi%IdofP(k,iel)+1)-1
            daux1 = 0.0_DP
            do i = i1, i2
              daux1 = daux1 + p_DC(i)*p_DvecP(p_KcolC(i))
            end do
            p_runi%Dg(k) = p_runi%Dg(k) - Dmult(3,3)*daux1
            ! Get the main diagonal entry of C
            p_runi%Dc(k) = Dmult(3,3)*p_DC(p_KdiagC(p_runi%IdofP(k,iel)))
          end do
        end if
        
        ! Invert A1 and A2
        do i = 1, ndofV
          p_runi%Da1(i) = 1.0_DP / p_runi%Da1(i)
          p_runi%Da2(i) = 1.0_DP / p_runi%Da2(i)
        end do
        
        ! Precalculate D * A^-1
        do i = 1, ndofV
          do j = 1, ndofP
            p_runi%Dt1(j,i) = p_runi%Dd1(j,i)*p_runi%Da1(i)
            p_runi%Dt2(j,i) = p_runi%Dd2(j,i)*p_runi%Da2(i)
          end do
        end do

        ! Calculate Schur-Complement of A
        ! S := -C + D * A^-1 * B 
        do j = 1, ndofP
          p_runi%Ds(j) = -p_runi%Dc(j)
          do i = 1, ndofV
            p_runi%Ds(j) = p_runi%Ds(j) + p_runi%Dt1(j,i)*p_runi%Db1(i,j) &
                                        + p_runi%Dt2(j,i)*p_runi%Db2(i,j)
          end do
        end do
        
        ! Calculate pressure
        ! p := S^-1 * (D * A^-1 * f_u - f_p)
        do j = 1, ndofP
          daux = -p_runi%Dg(j)
          do i = 1, ndofV
            daux = daux + p_runi%Dt1(j,i)*p_runi%Df1(i) &
                        + p_runi%Dt2(j,i)*p_runi%Df2(i)
          end do
          p_runi%Dup(j) = daux / p_runi%Ds(j)
        end do
        
        ! Calculate X- and Y-velocity
        do i = 1, ndofV
          do j = 1, ndofP
            p_runi%Df1(i) = p_runi%Df1(i) - p_runi%Db1(i,j)*p_runi%Dup(j)
            p_runi%Df2(i) = p_runi%Df2(i) - p_runi%Db2(i,j)*p_runi%Dup(j)
          end do
          p_runi%Du1(i) = p_runi%Da1(i)*p_runi%Df1(i)
          p_runi%Du2(i) = p_runi%Da2(i)*p_runi%Df2(i)
        end do

        ! Incorporate our local solution into the global one.
        do i = 1, ndofV
          j = p_runi%IdofV(i,iel)
          p_DvecU(j) = p_DvecU(j) + domega * p_runi%Du1(i)
          p_DvecV(j) = p_DvecV(j) + domega * p_runi%Du2(i)
        end do
        do i = 1, ndofP
          j = p_runi%IdofP(i,iel)
          p_DvecP(j) = p_DvecP(j) + domega * p_runi%Dup(i)
        end do
      
      end do ! iel
      
      NELdone = NELdone + NELtodo
    end do

  end subroutine

  ! ***************************************************************************
  
  subroutine vanka_NS2D_universal_full(rvanka, rvector, rrhs, domega, &
                                       IelementList, ndofV, ndofP)
  
!<description>
  ! This routine applies the specialised diagonal VANKA algorithm for
  ! 2D Navier-Stokes problems with arbitrary discretisation to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
!</description>

!<input>
  ! t_vankaPointerNavSt2D structure that saves algorithm-specific parameters.
  type(t_vankaPointerNavSt2D), target, intent(IN) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(IN)                    :: domega
  
  ! A list of element numbers where VANKA should be applied to.
  integer, dimension(:), intent(IN)       :: IelementList
  
  ! The number of local doFs.
  integer, intent(IN)                     :: ndofV, ndofP
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(IN)         :: rvector
!</inputoutput>

!</subroutine>
  
  type(t_vankaNavSt2D_uni_full), pointer :: p_runi

  real(DP), dimension(:), pointer :: p_DrhsU,p_DrhsV,p_DrhsP,&
                                     p_DvecU,p_DvecV,p_DvecP
  
  ! Multiplication factors
  real(DP), dimension(3,3) :: Dmult
  
  ! Quick access for the matrix arrays
  integer, dimension(:), pointer :: p_KldA,p_KldA12,p_KldB,p_KldC,p_KldD,&
      p_KcolA,p_KcolA12,p_KcolB,p_KcolC,p_KcolD
  real(DP), dimension(:), pointer :: p_DA11,p_DA12,p_DA21,p_DA22,p_DB1,p_DB2,&
      p_DD1,p_DD2,p_DC
  
  ! local variables
  integer :: i,j,iel,i1,i2,k,l,o,p,ndof
  integer :: NEL,NELdone,NELtodo
  real(DP) :: daux,daux1,daux2
  logical :: bHaveA12,bHaveC
  
  ! variables for LAPACK's DGESV routine
  integer :: info

    ! Get the number of elements to process
    NEL = rvanka%p_rspatialDiscrV%p_rtriangulation%NEL
    
    ! Calculate the total number of DOFs
    ndof = 2*ndofV + ndofP
    
    ! Get the diagonal vanka structure
    p_runi => rvanka%runi_full
    
    ! Get the pointers to the vector data
    call lsyssc_getbase_double(rvector%RvectorBlock(1), p_DvecU)
    call lsyssc_getbase_double(rvector%RvectorBlock(2), p_DvecV)
    call lsyssc_getbase_double(rvector%RvectorBlock(3), p_DvecP)
    call lsyssc_getbase_double(rrhs%RvectorBlock(1), p_DrhsU)
    call lsyssc_getbase_double(rrhs%RvectorBlock(2), p_DrhsV)
    call lsyssc_getbase_double(rrhs%RvectorBlock(3), p_DrhsP)
    
    ! Let's assume we do not have the optional matrices
    bHaveA12 = .FALSE.
    bHaveC = .FALSE.
    
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
      bHaveC = .TRUE.
      p_KldC => rvanka%p_KldC
      p_KcolC => rvanka%p_KcolC
      p_DC => rvanka%p_DC
    end if
    
    ! Get the multiplication factors
    Dmult = rvanka%Dmultipliers

    ! Take care of the 'soft-deactivation' of the sub-matrices
    bHaveA12 = bHaveA12 .and. ((Dmult(1,2) .ne. 0.0_DP) .or. &
                               (Dmult(2,1) .ne. 0.0_DP))
    bHaveC = bHaveC .and. (Dmult(3,3) .ne. 0.0_DP)
    
    ! Okay, loop through the element sets
    NELdone = 0
    do while(NELdone .lt. NEL)
    
      ! How many elements do we process this time?
      NELtodo = MIN(NEL-NELdone,VANKA_NAVST2D_NELEMSIM)
      
      ! Perform the DOF-mapping for the pressure and
      call dof_locGlobMapping_mult(rvanka%p_rspatialDiscrV,&
          IelementList(NELdone+1:NELdone+NELtodo), p_runi%IdofV)
      call dof_locGlobMapping_mult(rvanka%p_rspatialDiscrP,&
          IelementList(NELdone+1:NELdone+NELtodo), p_runi%IdofP)

      ! General case
      do iel=1, NELtodo
      
        ! First of all, fetch the local RHS
        do i = 1, ndofV
          p_runi%Df(i)       = p_DrhsU(p_runi%IdofV(i,iel))   ! f_u
          p_runi%Df(ndofV+i) = p_DrhsV(p_runi%IdofV(i,iel))   ! f_v
        end do
        o = 2*ndofV
        do i = 1, ndofP
          p_runi%Df(o+i) = p_DrhsP(p_runi%IdofP(i,iel))       ! f_p
        end do
        
        ! Clear the local system matrix
        p_runi%Da = 0.0_DP
        
        ! Let's update the local RHS vector by subtracting A*u from it:
        ! f_u := f_u - A11*u
        ! f_v := f_v - A22*v
        p = ndofV
        do k = 1, ndofV
          i1 = p_KldA(p_runi%IdofV(k,iel))
          i2 = p_KldA(p_runi%IdofV(k,iel)+1)-1
          daux1 = 0.0_DP
          daux2 = 0.0_DP
          do i = i1, i2
            j = p_KcolA(i)
            daux1 = daux1 + p_DA11(i)*p_DvecU(j)
            daux2 = daux2 + p_DA22(i)*p_DvecV(j)
            do l = 1, ndofV
              if(j .eq. p_runi%IdofV(l,iel)) then
                p_runi%Da(  k,  l) = Dmult(1,1)*p_DA11(i)
                p_runi%Da(p+k,p+l) = Dmult(2,2)*p_DA22(i)
                exit
              end if
            end do
          end do
          p_runi%Df(      k) = p_runi%Df(      k) - Dmult(1,1)*daux1
          p_runi%Df(ndofV+k) = p_runi%Df(ndofV+k) - Dmult(2,2)*daux2
        end do

        if(bHaveA12) then
          ! f_u := f_u - A12*v
          ! f_v := f_v - A21*u
          do k = 1, ndofV
            i1 = p_KldA12(p_runi%IdofV(k,iel))
            i2 = p_KldA12(p_runi%IdofV(k,iel)+1)-1
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = i1, i2
              j = p_KcolA12(i)
              daux1 = daux1 + p_DA12(i)*p_DvecV(j)
              daux2 = daux2 + p_DA21(i)*p_DvecU(j)
              do l = 1, ndofV
                if(j .eq. p_runi%IdofV(l,iel)) then
                  p_runi%Da(  k,p+l) = Dmult(1,2)*p_DA12(i)
                  p_runi%Da(p+k,  l) = Dmult(2,1)*p_DA21(i)
                  exit
                end if
              end do
            end do
            p_runi%Df(      k) = p_runi%Df(      k) - Dmult(1,2)*daux1
            p_runi%Df(ndofV+k) = p_runi%Df(ndofV+k) - Dmult(2,1)*daux2
          end do
        end if
        
        ! Now we also need to subtract B*p from our RHS, and by the same time,
        ! we will build the local B matrices.
        ! f_u := f_u - B1*p
        ! f_v := f_v - B2*p
        o = ndofV
        p = 2*ndofV
        do k = 1, ndofV
          i1 = p_KldB(p_runi%IdofV(k,iel))
          i2 = p_KldB(p_runi%IdofV(k,iel)+1)-1
          daux1 = 0.0_DP
          daux2 = 0.0_DP
          do i = i1, i2
            j = p_KcolB(i)
            daux = p_DvecP(j)
            daux1 = daux1 + p_DB1(i)*daux
            daux2 = daux2 + p_DB2(i)*daux2
            do l = 1, ndofP
              if(j .eq. p_runi%IdofP(l,iel)) then
                p_runi%Da(  k,p+l) = Dmult(1,3)*p_DB1(i)
                p_runi%Da(o+k,p+l) = Dmult(2,3)*p_DB2(i)
                exit
              end if
            end do
          end do
          p_runi%Df(      k) = p_runi%Df(      k) - Dmult(1,3)*daux1
          p_runi%Df(ndofV+k) = p_runi%Df(ndofV+k) - Dmult(2,3)*daux2
        end do
        
        ! Now we also need to subtract D*u from our RHS, and by the same time,
        ! we will build the local D matrices.
        ! f_p := f_p - D1*u - D2*v
        o = ndofV
        p = 2*ndofV
        do k = 1, ndofP
          i1 = p_KldD(p_runi%IdofP(k,iel))
          i2 = p_KldD(p_runi%IdofP(k,iel)+1)-1
          daux1 = 0.0_DP
          daux2 = 0.0_DP
          do i = i1, i2
            j = p_KcolD(i)
            daux1 = daux1 + p_DD1(i)*p_DvecU(j)
            daux2 = daux2 + p_DD2(i)*p_DvecV(j)
            do l = 1, ndofV
              if(j .eq. p_runi%IdofV(l,iel)) then
                p_runi%Da(p+k,  l) = Dmult(3,1)*p_DD1(i)
                p_runi%Da(p+k,o+l) = Dmult(3,2)*p_DD2(i)
                exit
              end if
            end do
          end do
          p_runi%Df(p+k) = p_runi%Df(p+k) - Dmult(3,1)*daux1 - Dmult(3,2)*daux2
        end do

        if(bHaveC) then
          ! f_p := f_p - C*p
          o = 2*ndofV
          do k = 1, ndofP
            i1 = p_KldC(p_runi%IdofP(k,iel))
            i2 = p_KldC(p_runi%IdofP(k,iel)+1)-1
            daux1 = 0.0_DP
            do i = i1, i2
              j = p_KcolC(i)
              daux1 = daux1 + p_DC(i)*p_DvecP(j)
              do l = 1, ndofP
                if(j .eq. p_runi%IdofP(l,iel)) then
                  p_runi%Da(o+k,o+l) = Dmult(3,3)*p_DC(i)
                  exit
                end if
              end do
            end do
            p_runi%Df(o+k) = p_runi%Df(o+k) - Dmult(3,3)*daux1
          end do
        end if

        ! Solve the local system
        call DGESV(ndof,1,p_runi%Da,ndof,p_runi%Ipiv,p_runi%Df,ndof,info)
        
        ! Did DGESV fail?
        if(info .ne. 0) cycle

        ! Incorporate our local solution into the global one.
        do i = 1, ndofV
          j = p_runi%IdofV(i,iel)
          p_DvecU(j) = p_DvecU(j) + domega * p_runi%Df(i)
          p_DvecV(j) = p_DvecV(j) + domega * p_runi%Df(ndofV+i)
        end do
        o = 2*ndofV
        do i = 1, ndofP
          j = p_runi%IdofP(i,iel)
          p_DvecP(j) = p_DvecP(j) + domega * p_runi%Df(o+i)
        end do
      
      end do ! ielidx
      
      NELdone = NELdone + NELtodo

    end do
    
  end subroutine

end module