!##############################################################################
!# ****************************************************************************
!# <name> vanka_bouss2d </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the vanka driver for 2D Boussinesq systems.
!# The most general version of a 2D Boussinesq-System looks as follows:
!#
!# <verb>
!#                     / A11 A12 B1 M1 \
!#                     | A21 A22 B2 M2 |
!#                     | D1  D2  C  K  |
!#                     \ 0   0   0  N  /
!# </verb>
!#
!# The following matrices are optional: A12, A21, C, M1, M2, K
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
!# 5.) Either both M1 and M2 exist or none of them.
!#     If they exist, then they have the same matrix structure.
!#
!# Please note that in contrast to the 'old' Navier-Stokes Vanka methods,
!# the D-matrices must NOT be virtually transposed!
!#
!# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- \\
!# Solving the local Boussinesq System \\
!# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- \\
!# Our local Boussinesq system looks as follows:
!#
!# <verb>
!#                   / A11 A12 B1 M1 \   / u1 \   / f_u1 \
!#                   | A21 A22 B2 M2 | * | u2 | = | f_u2 |
!#                   | D1  D2  C  K  |   | p  |   | f_p  |
!#                   \ 0   0   0  N  /   \ t  /   \ f_t  /
!# </verb>
!#
!# First, we will rewrite our system:
!#
!# <verb>
!# A := / A11 A12 \     B := / B1 \     M := / M1 \    D := ( D1 D2 )
!#      \ A21 A22 /          \ B2 /          \ M2 /
!#
!# u := / u1 \     f_u := / f_u1 \
!#      \ u1 /            \ f_u2 /
!# </verb>
!#
!# Now our system looks as follows:
!#
!# <verb>
!#                   / A B M \   / u \   / f_u \
!#                   | D C K | * | p | = | f_p |
!#                   \ 0 0 N /   \ t /   \ f_t /
!# </verb>
!#
!#
!# 1.) Solve the last equation to get the temperature t:
!#
!#     <tex> $$                  t := N^{-1} * f_t  $$ </tex>
!#
!# 2.) Calculate the Schur-Complement of A:
!#
!#     <tex> $$                  S := -C + D * A^{-1} * B  $$ </tex>
!#
!# 3.) Calculate pressure:
!#
!#     <tex> $$ p := S^{-1} * (D * A^{-1} * (f_u - M * t) + K * t - f_p) $$ </tex>
!#
!# 4.) Calculate velocity:
!#
!#     <tex> $$ u := A^{-1} * (f_u - M * t - B * p) $$ </tex>
!#
!# </purpose>
!##############################################################################

module vanka_bouss2d

!$use omp_lib
  use fsystem
  use genoutput
  use storage
  use mprimitives
  use element
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock

  implicit none

  private
  
  public :: t_vankaPointerBouss2D
  public :: vanka_initBoussinesq2D
  public :: vanka_Boussinesq2D

!<constants>

!<constantblock description="Vanka type identifiers for the 2D Boussinesq Class">

  ! Diagonal-type VANKA
  integer, parameter, public :: VANKATP_BOUSS2D_DIAG       = 0

  ! 'Full' VANKA
  integer, parameter, public :: VANKATP_BOUSS2D_FULL       = 1
  
  ! Pressure-DOF based VANKA
  integer, parameter, public :: VANKATP_BOUSS2D_PDOF       = 2
  
  ! Pressure-DOF based VANKA, fast variant ("SPSOR")
  integer, parameter, public :: VANKATP_BOUSS2D_PDOF_FAST  = 3

!</constantblock>

!</constants>


!<types>

!<typeblock>
  
  ! A structure that saves matrix pointers for the 2D-Boussinesq Vanka driver.
  type t_vankaPointerBouss2D
  
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

    ! Pointer to the column structure of the B/D-matrices.
    integer, dimension(:), pointer              :: p_KcolB => null()
    
    ! Pointer to the row structure of the B/D-matrices
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

    ! Pointer to the column structure of the M1/M2-matrices.
    integer, dimension(:), pointer              :: p_KcolM => null()
    
    ! Pointer to the row structure of the M1/M2-matrices
    integer, dimension(:), pointer              :: p_KldM => null()
    
    ! Pointer to the entries of the M1-matrix
    real(DP), dimension(:), pointer             :: p_DM1 => null()

    ! Pointer to the entries of the M2-matrix
    real(DP), dimension(:), pointer             :: p_DM2 => null()

    ! Pointer to the column structure of the K-matrix.
    integer, dimension(:), pointer              :: p_KcolK => null()
    
    ! Pointer to the row structure of the K-matrix.
    integer, dimension(:), pointer              :: p_KldK => null()
    
    ! Pointer to the entries of the K-matrix
    real(DP), dimension(:), pointer             :: p_DK => null()

    ! Pointer to the column structure of the N-matrix.
    integer, dimension(:), pointer              :: p_KcolN => null()
    
    ! Pointer to the row structure of the N-matrix
    integer, dimension(:), pointer              :: p_KldN => null()
    
    ! Pointer to diagonal entries in the N-matrix
    integer, dimension(:), pointer              :: p_KdiagonalN => null()
    
    ! Pointer to the entries of the N-matrix
    real(DP), dimension(:), pointer             :: p_DN => null()

    ! Spatial discretisation structure for X-velocity
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscrU => null()
    
    ! Spatial discretisation structure for Y-velocity
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscrV => null()
    
    ! Spatial discretisation structure for pressure
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscrP => null()
    
    ! Spatial discretisation structure for temperature
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscrT => null()
    
    ! Multiplication factors for the submatrices; taken from the system matrix.
    real(DP), dimension(4,4) :: Dmultipliers
    
    ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
    ! Support for Pressure-DOF based Vanka
    
    ! Total number of DOFs in velocity space component
    integer :: ndofVelo = 0
    
    ! Total number of DOFs in pressure space
    integer :: ndofPres = 0
    
    ! Total number of DOFs in temperature space
    integer :: ndofTemp = 0
    
    ! Pointer to an array storing the Schur-Complements
    real(DP), dimension(:), pointer :: p_Dschur => null()

  end type
  
!</typeblock>

!</types>

contains

!<subroutine>
  
  subroutine vanka_initBoussinesq2D (rmatrix,rvanka,csubtype)
  
!<description>
  ! Initialises the VANKA variant for 2D Boussinesq problems
  ! for conformal discretisations.
  ! Checks if the "2D-Boussinesq" VANKA variant
  ! for conformal discretisations can be applied to the system given by rmatrix.
  ! If not, the program is stopped.
  !
  ! The substructure rvanka%rvanka2DBouss is intitialised according
  ! to the information provided in rmatrix.
!</description>

!<input>
  ! The system matrix of the linear system.
  type(t_matrixBlock), intent(in), target :: rmatrix

  ! Desired subtype
  integer, intent(in) :: csubtype
!</input>

!<inputoutput>
  ! t_vankaPointer2DSPNavSt structure that saves algorithm-specific parameters.
  type(t_vankaPointerBouss2D), intent(inout) :: rvanka
!</inputoutput>

!</subroutine>

    type(t_blockDiscretisation), pointer :: p_rblockDiscr
    
    ! Matrix must be 4x4.
    if ((rmatrix%nblocksPerCol .ne. 4) .or. (rmatrix%nblocksPerRow .ne. 4)) then
      call output_line ('System matrix is not 4x4.',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_initBoussinesq2D')
      call sys_halt()
    end if
    
    ! TODO: Do more checks
    
    ! The structure of A(1,3) must be identical to A(3,1) and
    ! that of A(2,3) must be identical to A(3,2).
    if ((rmatrix%RmatrixBlock(1,3)%NA .ne. rmatrix%RmatrixBlock(3,1)%NA) .or. &
        (rmatrix%RmatrixBlock(1,3)%NEQ .ne. rmatrix%RmatrixBlock(3,1)%NCOLS)) then
      call output_line ('Structure of B1 and B1^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_initBoussinesq2D')
      call sys_halt()
    end if

    if ((rmatrix%RmatrixBlock(2,3)%NA .ne. rmatrix%RmatrixBlock(3,2)%NA) .or. &
        (rmatrix%RmatrixBlock(2,3)%NEQ .ne. rmatrix%RmatrixBlock(3,2)%NCOLS)) then
      call output_line ('Structure of B2 and B2^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_initBoussinesq2D')
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
    call lsyssc_getbase_double(rmatrix%RmatrixBlock(4,4),&
        rvanka%p_DN)
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
    call lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(4,4),&
        rvanka%p_KcolN)
    call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(4,4), &
        rvanka%p_KldN)
    if (rmatrix%RmatrixBlock(1,1)%cmatrixFormat .eq. LSYSSC_MATRIX9) then
      call lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(1,1), &
                              rvanka%p_KdiagonalA)
    else
      rvanka%p_KdiagonalA => rvanka%p_KldA
    end if
    if (rmatrix%RmatrixBlock(4,4)%cmatrixFormat .eq. LSYSSC_MATRIX9) then
      call lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(4,4), &
                              rvanka%p_KdiagonalN)
    else
      rvanka%p_KdiagonalN => rvanka%p_KldN
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
    
    ! Are the M1/M2 matrices present?
    if(lsysbl_isSubmatrixPresent(rmatrix,1,4)) then
    
      call lsyssc_getbase_double(rmatrix%RmatrixBlock(1,4),&
          rvanka%p_DM1)
      call lsyssc_getbase_double(rmatrix%RmatrixBlock(2,4),&
          rvanka%p_DM2)

      call lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,4),&
          rvanka%p_KcolM)
      call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,4), &
          rvanka%p_KldM)
          
    end if
    
    ! Is the K-Matrix present?
    if(lsysbl_isSubmatrixPresent(rmatrix,3,4)) then

      call lsyssc_getbase_double(rmatrix%RmatrixBlock(3,4),&
          rvanka%p_DK)
          
      call lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(3,4),&
          rvanka%p_KcolK)
      call lsyssc_getbase_Kld(rmatrix%RmatrixBlock(3,4), &
          rvanka%p_KldK)
    
    end if

    ! Get the multiplication factors of the submatrices.
    rvanka%Dmultipliers(1:4,1:4) = &
        rmatrix%RmatrixBlock(1:4,1:4)%dscaleFactor
        
    ! Do we use the pressure-DOF based Vanka?
    if((csubtype .eq. VANKATP_BOUSS2D_PDOF_FAST)) then
    
      ! Yes, we do. In this case we need to allocate an array for the local
      ! Schur-Complement matrices.
      
      ! Determine the total number of DOFs in pressure space - this is equal
      ! to the number of rows (equations) of the matrix located in block (3,1)
      rvanka%ndofPres = rmatrix%RmatrixBlock(3,1)%NEQ
      
      ! And determine the total number of DOFs in one velocity component -
      ! this is needed by the 'fast' variant.
      rvanka%ndofVelo = rmatrix%RmatrixBlock(1,1)%NEQ
      
      ! And one more time for the temperature...
      rvanka%ndofTemp = rmatrix%RmatrixBlock(4,4)%NEQ
      
      ! Allocate the array.
      allocate(rvanka%p_Dschur(rvanka%ndofPres))

      ! Perform the data-depenent initialisation
      call vanka_initDataBS2D_pdof(rvanka)
      
      ! And we can immediately return here as the rest of the code in this
      ! routine is only used by the other Vanka variants.
      return

    end if

    ! Get the block discretisation structure from the matrix.
    p_rblockDiscr => rmatrix%p_rblockDiscrTest
    
    if (.not. associated(p_rblockDiscr)) then
      call output_line ('No discretisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_initBoussinesq2D')
      call sys_halt()
    end if
    
    ! Get the discretisation structure of U,V,P and T from the block
    ! discretisation structure.
    rvanka%p_rspatialDiscrU => p_rblockDiscr%RspatialDiscr(1)
    rvanka%p_rspatialDiscrV => p_rblockDiscr%RspatialDiscr(2)
    rvanka%p_rspatialDiscrP => p_rblockDiscr%RspatialDiscr(3)
    rvanka%p_rspatialDiscrT => p_rblockDiscr%RspatialDiscr(4)
    
    if (rvanka%p_rspatialDiscrU%inumFESpaces .ne. &
        rvanka%p_rspatialDiscrV%inumFESpaces) then
      call output_line (&
          'Discretisation structures of X- and Y-velocity incompatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_initBoussinesq2D')
      call sys_halt()
    end if

    if ((rvanka%p_rspatialDiscrP%inumFESpaces .ne. 1) .and. &
        (rvanka%p_rspatialDiscrP%inumFESpaces .ne. &
          rvanka%p_rspatialDiscrU%inumFESpaces)) then
      ! Either there must be only one element type for the pressure, or there one
      ! pressure element distribution for every velocity element distribution!
      ! If this is not the case, we cannot determine (at least not in reasonable time)
      ! which element type the pressure represents on a cell!
      call output_line (&
          'Discretisation structures of velocity and pressure incompatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_initBoussinesq2D')
      call sys_halt()
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine vanka_doneBoussinesq2D (rvanka)
  
!<description>
  ! Releases the VANKA variant for 2D Boussinesq problems
  ! for conformal discretisations.
!</description>

!<inputoutput>
  ! t_vankaPointerBouss2D structure that saves algorithm-specific parameters.
  type(t_vankaPointerBouss2D), intent(inout) :: rvanka
!</inputoutput>

!</subroutine>

    ! Release Schur-complement array
    if(associated(rvanka%p_Dschur)) deallocate(rvanka%p_Dschur)

  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine vanka_Boussinesq2D (rvanka,rvector,rrhs,domega,csubtype)
  
!<description>
  ! This routine applies the VANKA variant for 2D Boussinesq problems
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
  type(t_vectorBlock), intent(in)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in)                    :: domega

  ! The subtype of VANKA that should handle the above problem class.
  ! One of the VANKATP_BOUSS2D_xxxx constants, e.g. VANKATP_BOUSS2D_DIAG.
  integer :: csubtype
  
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(inout)         :: rvector

  ! t_vanka structure that saves algorithm-specific parameters.
  type(t_vankaPointerBouss2D), intent(inout) :: rvanka
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: ielementdist
  integer, dimension(:), pointer :: p_IelementList
  type(t_elementDistribution), pointer :: p_relementDistrU
  type(t_elementDistribution), pointer :: p_relementDistrV
  type(t_elementDistribution), pointer :: p_relementDistrP
  type(t_elementDistribution), pointer :: p_relementDistrT
    
    ! Do we use the pressure-DOF based Vanka?
    if(csubtype .eq. VANKATP_BOUSS2D_PDOF) then
      
      ! Yes, we do. So call the corresponding routine here.
      ! The pressure-DOF based Vanka is a 'black-box' algorithm which does
      ! not need any information about the underlying discretisations.
      !call vanka_BS2D_pdof(rrhs, domega, rvector, rvanka)
      print *, 'ERROR: vanka_Boussinesq2D'
      stop
      
      ! And return here, as the rest of the code in this routine is only
      ! used by the other Vanka variants.
      return
    
    else if(csubtype .eq. VANKATP_BOUSS2D_PDOF_FAST) then
    
      ! Yes, we do, but we use the 'fast variant'.
      !call vanka_BS2D_pdof_fast(rrhs, domega, rvector, rvanka)
      print *, 'ERROR: vanka_Boussinesq2D'
      stop
      
      ! And return here.
      return
    
    end if

    ! Loop through the element distributions of the velocity.
    do ielementdist = 1,rvanka%p_rspatialDiscrU%inumFESpaces
    
      ! Get the corresponding element distributions of U, V and P.
      p_relementDistrU => &
          rvanka%p_rspatialDiscrU%RelementDistr(ielementdist)
      p_relementDistrV => &
          rvanka%p_rspatialDiscrV%RelementDistr(ielementdist)
      p_relementDistrT => &
          rvanka%p_rspatialDiscrT%RelementDistr(ielementdist)
      
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
      call storage_getbase_int (p_relementDistrU%h_IelementList,p_IelementList)
      
      ! Which element combination do we have now?
      if ((elem_getPrimaryElement(p_relementDistrU%celement) .eq. EL_Q1T) .and. &
          (elem_getPrimaryElement(p_relementDistrV%celement) .eq. EL_Q1T) .and. &
          (elem_getPrimaryElement(p_relementDistrP%celement) .eq. EL_Q0) .and. &
          (elem_getPrimaryElement(p_relementDistrT%celement) .eq. EL_Q1)) then
        ! Q1~/Q1~/Q0/Q1 discretisation
        
        ! Which VANKA subtype do we have? The diagonal VANKA of the full VANKA?
        select case (csubtype)
        case (VANKATP_BOUSS2D_DIAG)
          ! Call the jacobi-style vanka
          call vanka_BS2D_Q1TQ0Q1_js(rvanka, &
                  rvector, rrhs, domega, p_IelementList)

        case (VANKATP_BOUSS2D_FULL)
          if (.not. associated(rvanka%p_DA12)) then
            ! Call the block-diagonal vanka
            call vanka_BS2D_Q1TQ0Q1_bd(rvanka, &
                  rvector, rrhs, domega, p_IelementList)
          else
            ! Call the fully coupled vanka
            call vanka_BS2D_Q1TQ0Q1_fc(rvanka, &
                  rvector, rrhs, domega, p_IelementList)
          end if
        
        case default
          call output_line ('Unknown VANKA subtype!',&
              OU_CLASS_ERROR,OU_MODE_STD,'vanka_Boussinesq2D')
          call sys_halt()
        
        end select
    
      else
        call output_line ('Unsupported discretisation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'vanka_Boussinesq2D')
        call sys_halt()
        
      end if
      
    end do
      
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_initDataBS2D_pdof(rvanka)

!<description>
  ! This routine performs a matrix-data-dependent initialisation of the
  ! pressure-DOF based Vanka.
!</description>

!<inputoutput>

  ! t_vankaPointerBouss2D structure that saves algorithm-specific parameters.
  type(t_vankaPointerBouss2D), intent(inout) :: rvanka
  
!</inputoutput>

!</subroutine>

  ! Multiplication factors
  real(DP), dimension(4,4) :: Dmult
  
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
    
    ! Get the Schur-complement array
    p_DS => rvanka%p_Dschur
    
    ! Get the multiplication factors
    Dmult = rvanka%Dmultipliers
    
    ! First of all determine the maximum number of velocity DOFs that one
    ! pressure DOF may be adjacent to.
    nmaxdofV = 0
    do idofp = 1, rvanka%ndofPres
      nmaxdofV = max(nmaxdofV, p_KldD(idofp+1)-p_KldD(idofp))
    end do ! idofp
    nmaxdofV = 2*nmaxdofV
    
    ! Okay, let us allocate the temporary data
    allocate(DA(nmaxdofV,nmaxdofV))
    allocate(DB(nmaxdofV))
    allocate(Ipivot(nmaxdofV))
    
    ! Let us loop over the pressure DOFs
    do idofp = 1, rvanka%ndofPres
    
      ! Format the local matrices
      DA = 0.0_DP
      DB = 0.0_DP
      dc = 0.0_DP
      
      ! Format the Schur-complement entry for this pressure DOF
      p_DS(idofP) = 0.0_DP
      
      ! Fetch the number of velocity DOFs adjacent to this pressure DOF
      ndofV = p_KldD(idofp+1) - p_KldD(idofp)
      
      ! If the C matrix exists, grab it is main diagonal entry
      if(bHaveC) then
        dc = Dmult(3,3)*p_DC(p_KdiagC(idofp))
      end if
      
      ! Let us loop over all velocity DOFs which are adjacent to the current
      ! pressure DOF.
      do id1 = p_KldD(idofp), p_KldD(idofp+1)-1
      
        ! Get the index of the velocity DOF
        idofu = p_KcolD(id1)
        
        ! Let us fetch the local A11/A22 matrices
        do i = p_KldA(idofu), p_KldA(idofu+1)-1
        
          ! Get the column index
          k = p_KcolA(i)
          
          ! Let us see if this corresponds to one of our local velocity DOFs
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
          ! Let us fetch the local A12/A21 matrices
          do i = p_KldA12(idofu), p_KldA12(idofu+1)-1
          
            ! Get the column index
            k = p_KcolA12(i)
            
            ! Let us see if this corresponds to one of our local velocity DOFs
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

        ! Let us fetch the local B matrices
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
      
      ! Okay, let us calculate the Schur-complement dc = C - D * A^-1 * B
      daux1 = 0.0_DP
      daux2 = 0.0_DP
      k = 1
      do id1 = p_KldD(idofp), p_KldD(idofp+1)-1
        daux1 = daux1 + p_DD1(id1)*DB(k)
        daux2 = daux2 + p_DD2(id1)*DB(ndofV+k)
        k = k+1
      end do
      dc = dc - Dmult(3,1)*daux1 - Dmult(3,2)*daux2
      
      ! Now if the Schur-complement matrix is regular, we will store the inverse
      ! in the corresponding entry in p_DS.
      
      ! Remark:
      ! For some unknown reason the Intel Fortran Compiler 10.1 somehow screws
      ! up with the following IF-statement when compiling in Release-Mode under
      ! Windows:
      !
      !   if(dabs(dc) .gt. SYS_EPSREAL_DP) p_DS(idofP) = 1.0_DP / dc
      !
      ! The problem has been 'solved' by using the following equivalent
      ! alternative:
      if((dc .gt. SYS_EPSREAL_DP) .or. (dc .lt. -SYS_EPSREAL_DP)) then
        p_DS(idofP) = 1.0_DP / dc
      end if
      
    end do ! idofp
    
    ! And release the temporary memory
    deallocate(Ipivot)
    deallocate(DB)
    deallocate(DA)
    
    ! That is it
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine vanka_BS2D_Q1TQ0Q1_js(rvanka, rvector, rrhs, domega, IelementList)
  
!<description>
  ! This routine applies the specialised diagonal VANKA algorithm for
  ! 2D Boussinesq problems with Q1~/Q0/Q1 discretisation to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
!</description>

!<input>
  ! t_vankaPointerBouss2D structure that saves algorithm-specific parameters.
  type(t_vankaPointerBouss2D), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in)                    :: domega
  
  ! A list of element numbers where VANKA should be applied to.
  integer, dimension(:), intent(in)       :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(in)         :: rvector
!</inputoutput>

!</subroutine>

  ! How many local DOFs do we have?
  integer, parameter :: ndofV = 4   ! Dofs per velocity
  integer, parameter :: ndofP = 1   ! Dofs per pressure
  integer, parameter :: ndofT = 4   ! Dofs per temperature
  integer, parameter :: ndof = 2*ndofV+ndofP+ndofT
  
  ! Triangulation information
  integer, dimension(:,:), pointer :: p_IverticesAtElement
  integer, dimension(:,:), pointer :: p_IedgesAtElement
  real(DP), dimension(:), pointer :: p_DrhsU,p_DrhsV,p_DrhsP,p_DrhsT,&
                                     p_DvecU,p_DvecV,p_DvecP,p_DvecT
  
  ! DOFs
  integer, dimension(ndofV) :: IdofV
  integer, dimension(ndofP) :: IdofP
  integer, dimension(ndofT) :: IdofT
  
  ! Variables for the local system
  real(DP), dimension(ndofV) :: Da1, Da2, Du1, Du2, Df1, Df2
  real(DP), dimension(ndofV,ndofP) :: Db1,Db2
  real(DP), dimension(ndofP,ndofV) :: Dd1, Dd2
  real(DP), dimension(ndofP) :: Dc, Ds, Dup, Dfp
  real(DP), dimension(ndofV,ndofT) :: Dm1, Dm2
  real(DP), dimension(ndofP,ndofT) :: Dk
  real(DP), dimension(ndofT) :: Dn, Dut, Dft
  
  ! temporary vectors
  real(DP), dimension(ndofP,ndofV) :: Dt1, Dt2
  
  ! Multiplication factors
  real(DP), dimension(4,4) :: Dmult
  
  ! Quick access for the matrix arrays
  integer, dimension(:), pointer :: p_KldA,p_KldA12,p_KldB,p_KldD,p_KldC,&
      p_KldM,p_KldK,p_KldN,p_KcolA,p_KcolA12,p_KcolB,p_KcolD,p_KcolC,&
      p_KcolM,p_KcolK,p_KcolN,p_KdiagA,p_KdiagC,p_KdiagN
  real(DP), dimension(:), pointer :: p_DA11,p_DA12,p_DA21,p_DA22,p_DB1,p_DB2,&
      p_DD1,p_DD2,p_DC,p_DM1,p_DM2,p_DK,p_DN
    
    ! local variables
    integer :: i,j,iel,ielidx,i1,i2,k,l
    real(DP) :: daux,daux1, daux2
    logical :: bHaveA12, bHaveC, bHaveK, bHaveM

    ! Get pointers to the  triangulation information
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
        p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
        p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    
    ! Get the pointers to the vector data
    call lsyssc_getbase_double(rvector%RvectorBlock(1), p_DvecU)
    call lsyssc_getbase_double(rvector%RvectorBlock(2), p_DvecV)
    call lsyssc_getbase_double(rvector%RvectorBlock(3), p_DvecP)
    call lsyssc_getbase_double(rvector%RvectorBlock(4), p_DvecT)
    call lsyssc_getbase_double(rrhs%RvectorBlock(1), p_DrhsU)
    call lsyssc_getbase_double(rrhs%RvectorBlock(2), p_DrhsV)
    call lsyssc_getbase_double(rrhs%RvectorBlock(3), p_DrhsP)
    call lsyssc_getbase_double(rrhs%RvectorBlock(4), p_DrhsT)
    
    ! Let us assume we do not have the optional matrices
    bHaveA12 = .false.
    bHaveC = .false.
    bHaveM = .false.
    bHaveK = .false.
    
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
    p_KldN => rvanka%p_KldN
    p_KcolN => rvanka%p_KcolN
    p_KdiagN => rvanka%p_KdiagonalN
    p_DN => rvanka%p_DN
    
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
    
    if(associated(rvanka%p_DM1)) then
      bHaveM = .true.
      p_KldM => rvanka%p_KldM
      p_KcolM => rvanka%p_KcolM
      p_DM1 => rvanka%p_DM1
      p_DM2 => rvanka%p_DM2
    end if
    
    if(associated(rvanka%p_DK)) then
      bHaveK = .true.
      p_Kldk => rvanka%p_KldK
      p_KcolK => rvanka%p_KcolK
      p_DK => rvanka%p_DK
    end if
    
    ! Get the multiplication factors
    Dmult = rvanka%Dmultipliers
    
    ! Clear the optional matrices
    Dc  = 0.0_DP
    Dm1 = 0.0_DP
    Dm2 = 0.0_DP
    Dk  = 0.0_DP
    
    ! Now which of the optional matrices are present?
    ! In the following if-tree we will take care of the most common combinations,
    ! and if none of them fits, there is an else-case which can take care of any
    ! combination - but it is a bit slower.
!    if((.not. bHaveA12) .and. (.not. bHaveC) .and. bHaveM .and. (.not. bHaveK)) then
!    else if((.not. bHaveA12) .and. (.not. bHaveC) .and. (.not. bHaveM) .and. bHaveK) then
!    else if((.not. bHaveA12) .and. bHaveC .and. (.not. bHaveM) .and. bHaveK) then
!    else
      ! General case
      do ielidx=1, size(IelementList)
      
        ! Get the element number which is to be processed.
        iel = IelementList(ielidx)
        
        ! Get all DOFs for this element
        IdofV(:) = p_IedgesAtElement(:,iel)
        IdofP(1) = iel
        IdofT(:) = p_IverticesAtElement(:,iel)

        ! First of all, fetch the local RHS
        do i = 1, ndofV
          Df1(i) = p_DrhsU(IdofV(i))   ! f_u
          Df2(i) = p_DrhsV(IdofV(i))   ! f_v
        end do
        do i = 1, ndofP
          Dfp(i) = p_DrhsP(IdofP(i))   ! f_p
        end do
        do i = 1, ndofT
          Dft(i) = p_DrhsT(IdofT(i))   ! f_t
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
        
        ! Do the M1/M2 matrices exist?
        if(bHaveM) then
          ! Now subtract M*t from our local RHS
          ! f_u := f_u - M1*t
          ! f_v := f_v - M2*t
          do k = 1, ndofV
            i1 = p_KldM(IdofV(k))
            i2 = p_KldM(IdofV(k)+1)-1
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = i1, i2
              j = p_KcolM(i)
              daux = p_DvecT(j)
              daux1 = daux1 + p_DM1(i)*daux
              daux2 = daux2 + p_DM2(i)*daux
              do l = 1, ndofT
                if(j .eq. IdofT(l)) then
                  Dm1(k,l) = Dmult(1,4)*p_DM1(i)
                  Dm2(k,l) = Dmult(2,4)*p_DM2(i)
                  exit
                end if
              end do
            end do
            Df1(k) = Df1(k) - Dmult(1,4)*daux1
            Df2(k) = Df2(k) - Dmult(2,4)*daux2
          end do
        end if
        
        ! What about K? Does it exist?
        if(bHaveK) then
          ! Yes, so update the local RHS
          ! f_p := f_p - K*t
          do k = 1, ndofP
            i1 = p_KldK(IdofP(k))
            i2 = p_KldK(IdofP(k)+1)-1
            daux1 = 0.0_DP
            do i = i1, i2
              j = p_KcolK(i)
              daux1 = daux1 + p_DK(i)*p_DvecT(j)
              do l = 1, ndofT
                if(j .eq. IdofT(l)) then
                  Dk(k,l) = Dmult(3,4)*p_DK(i)
                  exit
                end if
              end do
            end do
            Dfp(k) = Dfp(k) - Dmult(3,4)*daux1
          end do
        end if
        
        ! Okay, one task is still left - the N-matrix...
        ! f_t := f_t - N*t
        do k = 1, ndofT
          i1 = p_KldN(IdofT(k))
          i2 = p_KldN(IdofT(k)+1)-1
          daux1 = 0.0_DP
          do i = i1, i2
            daux1 = daux1 + p_DN(i)*p_DvecT(p_KcolN(i))
          end do
          Dft(k) = Dft(k) - Dmult(4,4)*daux1
          ! Fetch the main diagonal entry of N
          Dn(k) = Dmult(4,4)*p_DN(p_KdiagN(IdofT(k)))
        end do

        ! Invert A1 and A2
        do i = 1, ndofV
          Da1(i) = 1.0_DP / Da1(i)
          Da2(i) = 1.0_DP / Da2(i)
        end do

        ! Calculate temperature
        ! t := N^-1 * f_t
        do i = 1, ndofT
          Dut(i) = Dft(i) / Dn(i)
        end do
        
        ! Calculate new RHS
        ! f_u := f_u - M1*t
        ! f_v := f_v - M2*t
        do i = 1, ndofV
          do j = 1, ndofT
            Df1(i) = Df1(i) - Dm1(i,j)*Dut(j)
            Df2(i) = Df2(i) - Dm2(i,j)*Dut(j)
          end do
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
        ! p := S^-1 * (D * A^-1 * f_u + K * t - f_p)
        do j = 1, ndofP
          daux = -Dfp(j)
          do i = 1, ndofV
            daux = daux + Dt1(j,i)*Df1(i) &
                        + Dt2(j,i)*Df2(i)
          end do
          do i = 1, ndofT
            daux = daux + Dk(j,i)*Dut(i)
          end do
          Dup(j) = daux / Ds(j)
        end do
        
        ! Calculate X- and Y-velocity
        ! u := A^-1 * (g_u - B * p)
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
        do i = 1, ndofT
          j = IdofT(i)
          p_DvecT(j) = p_DvecT(j) + domega * Dut(i)
        end do
      
      end do ! ielidx
    !end if

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine vanka_BS2D_Q1TQ0Q1_bd(rvanka, rvector, rrhs, domega, IelementList)
  
!<description>
  ! This routine applies the specialised diagonal VANKA algorithm for
  ! 2D Boussinesq problems with Q1~/Q0/Q1 discretisation to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
!</description>

!<input>
  ! t_vankaPointerBouss2D structure that saves algorithm-specific parameters.
  type(t_vankaPointerBouss2D), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in)                    :: domega
  
  ! A list of element numbers where VANKA should be applied to.
  integer, dimension(:), intent(in)       :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(in)         :: rvector
!</inputoutput>

!</subroutine>

  ! How many local DOFs do we have?
  integer, parameter :: ndofV = 4   ! Dofs per velocity
  integer, parameter :: ndofP = 1   ! Dofs per pressure
  integer, parameter :: ndofT = 4   ! Dofs per temperature
  integer, parameter :: ndof = 2*ndofV+ndofP+ndofT
  
  ! Triangulation information
  integer, dimension(:,:), pointer :: p_IverticesAtElement
  integer, dimension(:,:), pointer :: p_IedgesAtElement
  real(DP), dimension(:), pointer :: p_DrhsU,p_DrhsV,p_DrhsP,p_DrhsT,&
                                     p_DvecU,p_DvecV,p_DvecP,p_DvecT
  
  ! DOFs
  integer, dimension(ndofV) :: IdofV
  integer, dimension(ndofP) :: IdofP
  integer, dimension(ndofT) :: IdofT
  
  ! Variables for the local system
  real(DP), dimension(ndofV) :: Du1, Du2, Df1, Df2
  real(DP), dimension(ndofP) :: Dup, Dfp
  real(DP), dimension(ndofT) :: Dut, Dft
  real(DP), dimension(ndofV,ndofV) :: Da1, Da2
  real(DP), dimension(ndofV,ndofP) :: Db1,Db2
  real(DP), dimension(ndofP,ndofV) :: Dd1, Dd2
  real(DP), dimension(ndofP,ndofP) :: Dc
  real(DP), dimension(ndofV,ndofT) :: Dm1, Dm2
  real(DP), dimension(ndofP,ndofT) :: Dk
  real(DP), dimension(ndofT,ndofT) :: Dn
  
  ! Local variables
  real(DP), dimension(ndofV,ndofV) :: Di1,Di2
  real(DP), dimension(ndofT,ndofT) :: Dj
  real(DP), dimension(ndofP,ndofV) :: Dt1, Dt2
  real(DP), dimension(ndofP,ndofP) :: Ds, Dsi
  
  ! Multiplication factors
  real(DP), dimension(4,4) :: Dmult
  
  ! Quick access for the matrix arrays
  integer, dimension(:), pointer :: p_KldA,p_KldB,p_KldC,p_KldD,p_KldM,&
      p_KldK,p_KldN,p_KcolA,p_KcolB,p_KcolC,p_KcolD,p_KcolM,p_KcolK,&
      p_KcolN,p_KdiagA,p_KdiagC,p_KdiagN
  real(DP), dimension(:), pointer :: p_DA11,p_DA22,p_DB1,p_DB2,&
      p_DD1,p_DD2,p_DC,p_DM1,p_DM2,p_DK,p_DN
  
  ! local variables
  integer :: i,j,iel,ielidx,i1,i2,k,l
  real(DP) :: daux, daux1, daux2
  logical :: bHaveC, bHaveM, bHaveK
  logical :: bsuccess1, bsuccess2, bsuccess3

    ! Get pointers to the  triangulation information
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
        p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
        p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    
    ! Get the pointers to the vector data
    call lsyssc_getbase_double(rvector%RvectorBlock(1), p_DvecU)
    call lsyssc_getbase_double(rvector%RvectorBlock(2), p_DvecV)
    call lsyssc_getbase_double(rvector%RvectorBlock(3), p_DvecP)
    call lsyssc_getbase_double(rvector%RvectorBlock(4), p_DvecT)
    call lsyssc_getbase_double(rrhs%RvectorBlock(1), p_DrhsU)
    call lsyssc_getbase_double(rrhs%RvectorBlock(2), p_DrhsV)
    call lsyssc_getbase_double(rrhs%RvectorBlock(3), p_DrhsP)
    call lsyssc_getbase_double(rrhs%RvectorBlock(4), p_DrhsT)
    
    ! Let us assume we do not have the optional matrices
    bHaveC = .false.
    bHaveM = .false.
    bHaveK = .false.
    
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
    p_KldN => rvanka%p_KldN
    p_KcolN => rvanka%p_KcolN
    p_KdiagN => rvanka%p_KdiagonalN
    p_DN => rvanka%p_DN
    
    if(associated(rvanka%p_DC)) then
      bHaveC = .true.
      p_KldC => rvanka%p_KldC
      p_KcolC => rvanka%p_KcolC
      p_KdiagC => rvanka%p_KdiagonalC
      p_DC => rvanka%p_DC
    end if
    
    if(associated(rvanka%p_DM1)) then
      bHaveM = .true.
      p_KldM => rvanka%p_KldM
      p_KcolM => rvanka%p_KcolM
      p_DM1 => rvanka%p_DM1
      p_DM2 => rvanka%p_DM2
    end if
    
    if(associated(rvanka%p_DK)) then
      bHaveK = .true.
      p_KldK => rvanka%p_KldK
      p_KcolK => rvanka%p_KcolK
      p_DK => rvanka%p_DK
    end if
    
    ! Get the multiplication factors
    Dmult = rvanka%Dmultipliers
    
    ! Clear the optional matrices
    Dc  = 0.0_DP
    Dm1 = 0.0_DP
    Dm2 = 0.0_DP
    Dk  = 0.0_DP
    
    !if((.not. bHaveC) .and. bHaveM .and. (.not. bHaveM3)) then
    !else if(bHaveC .and. (.not. bHaveM) .and. bHaveM3) then
    !else if((.not. bHaveC) .and. (.not. bHaveM) .and. bHaveM3) then
    !else
      ! General case
      do ielidx=1, size(IelementList)
      
        ! Get the element number which is to be processed.
        iel = IelementList(ielidx)
        
        ! Get all DOFs for this element
        IdofV(:) = p_IedgesAtElement(:,iel)
        IdofP(1) = iel
        IdofT(:) = p_IverticesAtElement(:,iel)

        ! First of all, fetch the local RHS
        do i = 1, ndofV
          Df1(i) = p_DrhsU(IdofV(i))   ! f_u
          Df2(i) = p_DrhsV(IdofV(i))   ! f_v
        end do
        do i = 1, ndofP
          Dfp(i) = p_DrhsP(IdofP(i))   ! f_p
        end do
        do i = 1, ndofT
          Dft(i) = p_DrhsT(IdofT(i))   ! f_t
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
        
        ! Do we have a C-matrix?
        if(bHaveC) then
          ! Yes, so update the local RHS
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
        end if
        
        if(bHaveM) then
          ! Now subtract M*t from our local RHS
          ! f_u := f_u - M1*t
          ! f_v := f_v - M2*t
          do k = 1, ndofV
            i1 = p_KldM(IdofV(k))
            i2 = p_KldM(IdofV(k)+1)-1
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = i1, i2
              j = p_KcolM(i)
              daux = p_DvecT(j)
              daux1 = daux1 + p_DM1(i)*daux
              daux2 = daux2 + p_DM2(i)*daux
              do l = 1, ndofT
                if(j .eq. IdofT(l)) then
                  Dm1(k,l) = Dmult(1,4)*p_DM1(i)
                  Dm2(k,l) = Dmult(2,4)*p_DM2(i)
                  exit
                end if
              end do
            end do
            Df1(k) = Df1(k) - Dmult(1,4)*daux1
            Df2(k) = Df2(k) - Dmult(2,4)*daux2
          end do
        end if
        
        ! What about K? Does it exist?
        if(bHaveK) then
          ! Yes, so update the local RHS
          ! f_p := f_p - K*t
          do k = 1, ndofP
            i1 = p_KldK(IdofP(k))
            i2 = p_KldK(IdofP(k)+1)-1
            daux1 = 0.0_DP
            do i = i1, i2
              j = p_KcolK(i)
              daux1 = daux1 + p_DK(i)*p_DvecT(j)
              do l = 1, ndofT
                if(j .eq. IdofT(l)) then
                  Dk(k,l) = Dmult(3,4)*p_DK(i)
                  exit
                end if
              end do
            end do
            Dfp(k) = Dfp(k) - Dmult(3,4)*daux1
          end do
        end if
        
        ! Okay, one task is still left - the N-matrix...
        ! f_t := f_t - N*t
        do k = 1, ndofT
          i1 = p_KldN(IdofT(k))
          i2 = p_KldN(IdofT(k)+1)-1
          daux1 = 0.0_DP
          do i = i1, i2
            j = p_KcolN(i)
            daux1 = daux1 + p_DN(i)*p_DvecT(j)
            do l = 1, ndofT
              if(j .eq. IdofT(l)) then
                Dn(k,l) = Dmult(4,4)*p_DN(i)
                exit
              end if
            end do
          end do
          Dft(k) = Dft(k) - Dmult(4,4)*daux1
        end do

        ! Invert A1, A2 and N
        call mprim_invert4x4MatrixDirectDble(Da1, Di1, bsuccess1)
        call mprim_invert4x4MatrixDirectDble(Da2, Di2, bsuccess2)
        call mprim_invert4x4MatrixDirectDble(Dn, Dj, bsuccess3)

        if (bsuccess1 .and. bsuccess2 .and. bsuccess3) then

          ! Calculate temperature
          ! t := N^-1 * f_t
          do i = 1, ndofT
            Dut(i) = 0.0_DP
            do j = 1, ndofT
              Dut(i) = Dut(i) + Dj(i,j)*Dft(j)
            end do
          end do
          
          ! Calculate new RHS
          ! f_u := f_u - M1*t
          ! f_v := f_v - M2*t
          do i = 1, ndofV
            do j = 1, ndofT
              Df1(i) = Df1(i) - Dm1(i,j)*Dut(j)
              Df2(i) = Df2(i) - Dm2(i,j)*Dut(j)
            end do
          end do
          
          ! Precalculate D * A^-1
          do i = 1, ndofV
            do j = 1, ndofP
              Dt1(j,i) = 0.0_DP
              Dt2(j,i) = 0.0_DP
              do k = 1, ndofV
                Dt1(j,i) = Dt1(j,i) + Dd1(j,k)*Di1(k,i)
                Dt2(j,i) = Dt2(j,i) + Dd2(j,k)*Di2(k,i)
              end do
            end do
          end do
          
          ! Calculate Schur-Complement of A
          ! S := -C + D * A^-1 * B
          do j = 1, ndofP
            do k = 1, ndofP
              Ds(j,k) = -Dc(j,k)
              do i = 1, ndofV
                Ds(j,k) = Ds(j,k) + Dt1(j,i)*Db1(i,k) &
                                  + Dt2(j,i)*Db2(i,k)
              end do
            end do
          end do
          
          ! Invert S
          Dsi(1,1) = 1.0_DP / Ds(1,1)
          
          ! Update pressure RHS
          do j = 1, ndofP
            daux = -Dfp(j)
            do i = 1, ndofV
              daux = daux + Dt1(j,i)*Df1(i) &
                          + Dt2(j,i)*Df2(i)
            end do
            do i = 1, ndofT
              daux = daux + Dk(j,i)*Dut(i)
            end do
            Dfp(j) = daux
          end do
            
          ! Calculate pressure
          ! p := S^-1 * f_p
          do j = 1, ndofP
            Dup(j) = 0.0_DP
            do i = 1, ndofP
              Dup(j) = Dup(j) + Dsi(j,i)*Dfp(i)
            end do
          end do

          ! Update velocity RHS
          ! f_u := f_u - B * p
          do i = 1, ndofV
            do j = 1, ndofP
              Df1(i) = Df1(i) - Db1(i,j)*Dup(j)
              Df2(i) = Df2(i) - Db2(i,j)*Dup(j)
            end do
          end do
          
          ! Calculate X- and Y-velocity
          ! u := A^-1 * f_u
          do i = 1, ndofV
            Du1(i) = 0.0_DP
            Du2(i) = 0.0_DP
            do j = 1, ndofV
              Du1(i) = Du1(i) + Di1(i,j)*Df1(j)
              Du2(i) = Du2(i) + Di2(i,j)*Df2(j)
            end do
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
          do i = 1, ndofT
            j = IdofT(i)
            p_DvecT(j) = p_DvecT(j) + domega * Dut(i)
          end do
          
        end if
      
      end do ! ielidx
    !end if

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine vanka_BS2D_Q1TQ0Q1_fc(rvanka, rvector, rrhs, domega, IelementList)
  
!<description>
  ! This routine applies the specialised diagonal VANKA algorithm for
  ! 2D Boussinesq problems with Q1~/Q0/Q1 discretisation to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
!</description>

!<input>
  ! t_vankaPointerBouss2D structure that saves algorithm-specific parameters.
  type(t_vankaPointerBouss2D), intent(in) :: rvanka

  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in)                    :: domega
  
  ! A list of element numbers where VANKA should be applied to.
  integer, dimension(:), intent(in)       :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(in)         :: rvector
!</inputoutput>

!</subroutine>

  ! How many local DOFs do we have?
  integer, parameter :: ndofV = 4   ! Dofs per velocity
  integer, parameter :: ndofP = 1   ! Dofs per pressure
  integer, parameter :: ndofT = 4   ! Dofs per temperature
  integer, parameter :: ndof = 2*ndofV+ndofP+ndofT
  
  ! Triangulation information
  integer, dimension(:,:), pointer :: p_IverticesAtElement
  integer, dimension(:,:), pointer :: p_IedgesAtElement
  real(DP), dimension(:), pointer :: p_DrhsU,p_DrhsV,p_DrhsP,p_DrhsT,&
                                     p_DvecU,p_DvecV,p_DvecP,p_DvecT
  
  ! DOFs
  integer, dimension(ndofV) :: IdofV
  integer, dimension(ndofP) :: IdofP
  integer, dimension(ndofT) :: IdofT
  
  ! Variables for the local system
  real(DP), dimension(ndof) :: Df
  real(DP), dimension(ndof,ndof) :: Da
  
  ! Multiplication factors
  real(DP), dimension(4,4) :: Dmult
  
  ! Quick access for the matrix arrays
  integer, dimension(:), pointer :: p_KldA,p_KldA12,p_KldB,p_KldD,p_KldC,&
      p_KldM,p_KldK,p_KldN,p_KcolA,p_KcolA12,p_KcolB,p_KcolD,p_KcolC,&
      p_KcolM,p_KcolK,p_KcolN
  real(DP), dimension(:), pointer :: p_DA11,p_DA12,p_DA21,p_DA22,p_DB1,p_DB2,&
      p_DD1,p_DD2,p_DC,p_DM1,p_DM2,p_DK,p_DN
  
  ! local variables
  integer :: i,j,iel,ielidx,i1,i2,k,l,o,p
  real(DP) :: daux, daux1, daux2
  logical :: bHaveC, bHaveM, bHaveK
  
  ! variables for LAPACK`s DGESV routine
  integer, dimension(ndof) :: Ipivot
  integer :: info

    ! Get pointers to the  triangulation information
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
        p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
        p_rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    
    ! Get the pointers to the vector data
    call lsyssc_getbase_double(rvector%RvectorBlock(1), p_DvecU)
    call lsyssc_getbase_double(rvector%RvectorBlock(2), p_DvecV)
    call lsyssc_getbase_double(rvector%RvectorBlock(3), p_DvecP)
    call lsyssc_getbase_double(rvector%RvectorBlock(4), p_DvecT)
    call lsyssc_getbase_double(rrhs%RvectorBlock(1), p_DrhsU)
    call lsyssc_getbase_double(rrhs%RvectorBlock(2), p_DrhsV)
    call lsyssc_getbase_double(rrhs%RvectorBlock(3), p_DrhsP)
    call lsyssc_getbase_double(rrhs%RvectorBlock(4), p_DrhsT)
    
    ! Let us assume we do not have the optional matrices
    bHaveC = .false.
    bHaveM = .false.
    bHaveK = .false.
    
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
    p_KldN => rvanka%p_KldN
    p_KcolN => rvanka%p_KcolN
    p_DN => rvanka%p_DN
    
    if(associated(rvanka%p_DC)) then
      bHaveC = .true.
      p_KldC => rvanka%p_KldC
      p_KcolC => rvanka%p_KcolC
      p_DC => rvanka%p_DC
    end if
    
    if(associated(rvanka%p_DM1)) then
      bHaveM = .true.
      p_KldM => rvanka%p_KldM
      p_KcolM => rvanka%p_KcolM
      p_DM1 => rvanka%p_DM1
      p_DM2 => rvanka%p_DM2
    end if
    
    if(associated(rvanka%p_DK)) then
      bHaveK = .true.
      p_KldK => rvanka%p_KldK
      p_KcolK => rvanka%p_KcolK
      p_DK => rvanka%p_DK
    end if
    
    ! Get the multiplication factors
    Dmult = rvanka%Dmultipliers
    
    ! So we start with a loop over all elements in the list
    do ielidx=1, size(IelementList)
    
      ! Get the element number which is to be processed.
      iel = IelementList(ielidx)
      
      ! Get all DOFs for this element
      IdofV(:) = p_IedgesAtElement(:,iel)
      IdofP(1) = iel
      IdofT(:) = p_IverticesAtElement(:,iel)

      ! First of all, fetch the local RHS
      do i = 1, ndofV
        Df(i)       = p_DrhsU(IdofV(i))   ! f_u
        Df(ndofV+i) = p_DrhsV(IdofV(i))   ! f_v
      end do
      o = 2*ndofV
      do i = 1, ndofP
        Df(o+i) = p_DrhsP(IdofP(i))       ! f_p
      end do
      o = 2*ndofV + ndofP
      do i = 1, ndofT
        Df(o+i) = p_DrhsT(IdofT(i))       ! f_t
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
      
      ! Do we have a C-matrix?
      if(bHaveC) then
        ! Yes, so update the local RHS
        ! f_p := f_p - C*p
        o = 2*ndofV
        do k = 1, ndofP
          i1 = p_KldC(IdofP(k))
          i2 = p_KldC(IdofP(k)+1)-1
          daux1 = 0.0_DP
          do i = i1, i2
            j = p_KcolC(i)
            daux1 = daux1 + p_DC(i)*p_DvecP(j)
            do l = 1, ndofP
              if(j .eq. IdofP(l)) then
                Da(o+k,o+l) = Dmult(3,3)*p_DC(i)
                exit
              end if
            end do
          end do
          Df(o+k) = Df(o+k) - Dmult(3,3)*daux1
        end do
      end if
      
      ! Do the M1/M2 matrices exist?
      if(bHaveM) then
        ! Now subtract M*t from our local RHS
        ! f_u := f_u - M1*t
        ! f_v := f_v - M2*t
        o = ndofV
        p = 2*ndofV+ndofP
        do k = 1, ndofV
          i1 = p_KldM(IdofV(k))
          i2 = p_KldM(IdofV(k)+1)-1
          daux1 = 0.0_DP
          daux2 = 0.0_DP
          do i = i1, i2
            j = p_KcolM(i)
            daux = p_DvecT(j)
            daux1 = daux1 + p_DM1(i)*daux
            daux2 = daux2 + p_DM2(i)*daux
            do l = 1, ndofT
              if(j .eq. IdofT(l)) then
                Da(  k,p+l) = Dmult(1,4)*p_DM1(i)
                Da(o+k,p+l) = Dmult(2,4)*p_DM2(i)
                exit
              end if
            end do
          end do
          Df(      k) = Df(      k) - Dmult(1,4)*daux1
          Df(ndofV+k) = Df(ndofV+k) - Dmult(2,4)*daux2
        end do
      end if
      
      ! What about K? Does it exist?
      if(bHaveK) then
        ! Yes, so update the local RHS
        ! f_p := f_p - K*t
        o = 2*ndofV
        p = 2*ndofV+ndofP
        do k = 1, ndofP
          i1 = p_KldK(IdofP(k))
          i2 = p_KldK(IdofP(k)+1)-1
          daux1 = 0.0_DP
          do i = i1, i2
            j = p_KcolK(i)
            daux1 = daux1 + p_DK(i)*p_DvecT(j)
            do l = 1, ndofT
              if(j .eq. IdofT(l)) then
                Da(o+k,p+l) = Dmult(3,4)*p_DK(i)
                exit
              end if
            end do
          end do
          Df(o+k) = Df(o+k) - Dmult(3,4)*daux1
        end do
      end if
      
      ! Okay, one task is still left - the N-matrix...
      ! f_t := f_t - N*t
      o = 2*ndofV+ndofP
      do k = 1, ndofT
        i1 = p_KldN(IdofT(k))
        i2 = p_KldN(IdofT(k)+1)-1
        daux1 = 0.0_DP
        do i = i1, i2
          j = p_KcolN(i)
          daux1 = daux1 + p_DN(i)*p_DvecT(j)
          do l = 1, ndofT
            if(j .eq. IdofT(l)) then
              Da(o+k,o+l) = Dmult(4,4)*p_DN(i)
              exit
            end if
          end do
        end do
        Df(o+k) = Df(o+k) - Dmult(4,4)*daux1
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
      o = 2*ndofV + ndofP
      do i = 1, ndofT
        j = IdofT(i)
        p_DvecT(j) = p_DvecT(j) + domega * Df(o+i)
      end do
    
    end do ! ielidx
  
  end subroutine

end module
