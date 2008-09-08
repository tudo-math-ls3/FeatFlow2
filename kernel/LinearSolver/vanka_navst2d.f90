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
!# 3.) B1, B2, D1 and D2 have the same matrix structure and both
!#     D1 and D2 are 'virtually transposed'.
!#
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
!# </purpose>
!##############################################################################

module vanka_navst2d

use fsystem
use genoutput
use mprimitives
use spatialdiscretisation
use linearsystemblock

implicit none

!<constants>

!<constantblock description="Vanka type identifiers for the 2D Navier-Stokes Class">

!</constantblock>

  ! Diagonal-type VANKA
  integer, parameter :: VANKATP_NAVST2D_DIAG  = 0

  ! 'Full' VANKA
  integer, parameter :: VANKATP_NAVST2D_FULL  = 1

!<constants>


!<types>

!<typeblock>
  
  ! A structure that saves matrix pointers for the 2D Navier-Stokes Vanka driver.
  TYPE t_vankaPointerNavSt2D
    
    ! Pointer to the column structure of the velocity matrix A11/A22
    INTEGER, DIMENSION(:), POINTER              :: p_KcolA => NULL()
    
    ! Pointer to the row structure of the velocity matrix A11/A22
    INTEGER, DIMENSION(:), POINTER              :: p_KldA => NULL()
    
    ! Pointer to diagonal entries in the velocity matrix A11/A22
    INTEGER, DIMENSION(:), POINTER              :: p_KdiagonalA => NULL()

    ! Pointer to the matrix entries of the velocity matrix A11
    REAL(DP), DIMENSION(:), POINTER             :: p_DA11 => NULL()

    ! Pointer to the matrix entries of the velocity matrix A22
    REAL(DP), DIMENSION(:), POINTER             :: p_DA22 => NULL()

    ! Pointer to the column structure of the matrix A12/A21
    INTEGER, DIMENSION(:), POINTER              :: p_KcolA12 => NULL()
    
    ! Pointer to the row structure of the matrix A12/A21
    INTEGER, DIMENSION(:), POINTER              :: p_KldA12 => NULL()

    ! Pointer to the matrix entries of the velocity matrix A12
    REAL(DP), DIMENSION(:), POINTER             :: p_DA12 => NULL()

    ! Pointer to the matrix entries of the velocity matrix A21
    REAL(DP), DIMENSION(:), POINTER             :: p_DA21 => NULL()

    ! Pointer to the column structure of the B/D-matrices.
    INTEGER, DIMENSION(:), POINTER              :: p_KcolB => NULL()
    
    ! Pointer to the row structure of the B/D-matrices
    INTEGER, DIMENSION(:), POINTER              :: p_KldB => NULL()
    
    ! Pointer to the entries of the B1-matrix
    REAL(DP), DIMENSION(:), POINTER             :: p_DB1 => NULL()

    ! Pointer to the entries of the B2-matrix
    REAL(DP), DIMENSION(:), POINTER             :: p_DB2 => NULL()
    
    ! Pointer to the entries of the D1-matrix
    REAL(DP), DIMENSION(:), POINTER             :: p_DD1 => NULL()

    ! Pointer to the entries of the D2-matrix
    REAL(DP), DIMENSION(:), POINTER             :: p_DD2 => NULL()

    ! Pointer to the column structure of the C-matrix.
    INTEGER, DIMENSION(:), POINTER              :: p_KcolC => NULL()
    
    ! Pointer to the row structure of the C-matrix.
    INTEGER, DIMENSION(:), POINTER              :: p_KldC => NULL()
    
    ! Pointer to diagonal entries in the C-matrix
    INTEGER, DIMENSION(:), POINTER              :: p_KdiagonalC => NULL()
    
    ! Pointer to the entries of the C-matrix
    REAL(DP), DIMENSION(:), POINTER             :: p_DC => NULL()

    ! Spatial discretisation structure for X-velocity
    TYPE(t_spatialDiscretisation), POINTER :: p_rspatialDiscrU => NULL()
    
    ! Spatial discretisation structure for Y-velocity
    TYPE(t_spatialDiscretisation), POINTER :: p_rspatialDiscrV => NULL()
    
    ! Spatial discretisation structure for pressure
    TYPE(t_spatialDiscretisation), POINTER :: p_rspatialDiscrP => NULL()
    
    ! Multiplication factors for the submatrices; taken from the system matrix.
    REAL(DP), DIMENSION(3,3) :: Dmultipliers
    
  END TYPE
  
!</typeblock>

!</types>

contains

!<subroutine>
  
  SUBROUTINE vanka_initNavierStokes2D (rmatrix,rvanka,csubtype)
  
!<description>
  ! Initialises the VANKA variant for 2D Navier-Stokes problems 
  ! for conformal discretisations.
  ! Checks if the "2D-Navier-Stokes" VANKA variant 
  ! for conformal discretisations can be applied to the system given by rmatrix.
  ! If not, the program is stopped.
!</description>

!<input>
  ! The system matrix of the linear system.
  TYPE(t_matrixBlock), INTENT(IN), TARGET :: rmatrix

  ! Desired subtype
  INTEGER, INTENT(IN) :: csubtype  
!</input>

!<inputoutput>
  ! t_vankaPointerNavSt2D structure that saves algorithm-specific parameters.
  TYPE(t_vankaPointerNavSt2D), INTENT(INOUT) :: rvanka
!</inputoutput>

!</subroutine>

    INTEGER :: i,j
    TYPE(t_blockDiscretisation), POINTER :: p_rblockDiscr
    
    ! Matrix must be 3x3.
    IF (rmatrix%ndiagBlocks .NE. 3) THEN
      CALL output_line ('System matrix is not 3x3.',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_initNavierStokes2D')
      CALL sys_halt()
    END IF
    
    ! TODO: Do more checks
    
    ! The structure of A(1,3) must be identical to A(3,1) and
    ! that of A(2,3) must be identical to A(3,2).
    IF ((rmatrix%RmatrixBlock(1,3)%NA .NE. rmatrix%RmatrixBlock(3,1)%NA) .OR. &
        (rmatrix%RmatrixBlock(1,3)%NEQ .NE. rmatrix%RmatrixBlock(3,1)%NCOLS)) THEN
      CALL output_line ('Structure of B1 and B1^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_initNavierStokes2D')
      CALL sys_halt()
    END IF

    IF ((rmatrix%RmatrixBlock(2,3)%NA .NE. rmatrix%RmatrixBlock(3,2)%NA) .OR. &
        (rmatrix%RmatrixBlock(2,3)%NEQ .NE. rmatrix%RmatrixBlock(3,2)%NCOLS)) THEN
      CALL output_line ('Structure of B2 and B2^T different!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_initNavierStokes2D')
      CALL sys_halt()
    END IF      
  
    ! Fill the output structure with data of the matrices.
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(1,1),&
        rvanka%p_DA11)
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(2,2),&
        rvanka%p_DA22)
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(1,3),&
        rvanka%p_DB1)
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(2,3),&
        rvanka%p_DB2)
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(3,1),&
        rvanka%p_DD1)
    CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(3,2),&
        rvanka%p_DD2)
    CALL lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,1),&
        rvanka%p_KcolA)
    CALL lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,1), &
        rvanka%p_KldA )
    CALL lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(1,3),&
        rvanka%p_KcolB)
    CALL lsyssc_getbase_Kld(rmatrix%RmatrixBlock(1,3), &
        rvanka%p_KldB)
    IF (rmatrix%RmatrixBlock(1,1)%cmatrixFormat .EQ. LSYSSC_MATRIX9) THEN
      CALL lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(1,1), &
                              rvanka%p_KdiagonalA)
    ELSE
      rvanka%p_KdiagonalA => rvanka%p_KldA
    END IF
    
    ! Are the A12/A21 matrices present?
    IF (lsysbl_isSubmatrixPresent(rmatrix,1,2)) THEN
      
      CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(1,2),&
          rvanka%p_DA12 )
      
      CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(2,1),&
          rvanka%p_DA21 )
    END IF
    
    ! Is the C-Matrix present?
    IF (lsysbl_isSubmatrixPresent(rmatrix,3,3)) THEN
    
      CALL lsyssc_getbase_double(rmatrix%RmatrixBlock(3,3),&
          rvanka%p_DC)

      CALL lsyssc_getbase_Kcol(rmatrix%RmatrixBlock(3,3),&
          rvanka%p_KcolC)
      CALL lsyssc_getbase_Kld(rmatrix%RmatrixBlock(3,3), &
          rvanka%p_KldC)

      IF (rmatrix%RmatrixBlock(3,3)%cmatrixFormat .EQ. LSYSSC_MATRIX9) THEN
        CALL lsyssc_getbase_Kdiagonal(rmatrix%RmatrixBlock(3,3), &
                                rvanka%p_KdiagonalC)
      ELSE
        rvanka%p_KdiagonalC => rvanka%p_KldC
      END IF
    
    END IF
    
    ! Get the multiplication factors of the submatrices.
    rvanka%Dmultipliers(1:3,1:3) = &
        rmatrix%RmatrixBlock(1:3,1:3)%dscaleFactor

    ! Get the block discretisation structure from the matrix.
    p_rblockDiscr => rmatrix%p_rblockDiscrTest
    
    IF (.NOT. ASSOCIATED(p_rblockDiscr)) THEN
      CALL output_line ('No discretisation!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_initNavierStokes2D')
      CALL sys_halt()
    END IF
    
    ! Get the discretisation structure of U,V and P from the block
    ! discretisation structure.
    rvanka%p_rspatialDiscrU => p_rblockDiscr%RspatialDiscr(1)
    rvanka%p_rspatialDiscrV => p_rblockDiscr%RspatialDiscr(2)
    rvanka%p_rspatialDiscrP => p_rblockDiscr%RspatialDiscr(3)
    
    IF (rvanka%p_rspatialDiscrU%inumFESpaces .NE. &
        rvanka%p_rspatialDiscrV%inumFESpaces) THEN
      CALL output_line (&
          'Discretisation structures of X- and Y-velocity incompatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_initNavierStokes2D')
      CALL sys_halt()
    END IF

    IF ((rvanka%p_rspatialDiscrP%inumFESpaces .NE. 1) .AND. &
        (rvanka%p_rspatialDiscrP%inumFESpaces .NE. &
          rvanka%p_rspatialDiscrU%inumFESpaces)) THEN
      ! Either there must be only one element type for the pressure, or there one
      ! pressure element distribution for every velocity element distribution!
      ! If this is not the case, we cannot determine (at least not in reasonable time)
      ! which element type the pressure represents on a cell!
      CALL output_line (&
          'Discretisation structures of velocity and pressure incompatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'vanka_initNavierStokes2D')
      CALL sys_halt()
    END IF

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE vanka_NavierStokes2D (rvanka,rvector,rrhs,domega,csubtype)
  
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
  TYPE(t_vectorBlock), INTENT(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  REAL(DP), INTENT(IN)                    :: domega

  ! The subtype of VANKA that should handle the above problem class.
  ! One of the VANKATP_BOUSS2D_xxxx constants, e.g. VANKATP_BOUSS2D_DIAG.
  INTEGER :: csubtype
  
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  TYPE(t_vectorBlock), INTENT(INOUT)         :: rvector

  ! t_vanka structure that saves algorithm-specific parameters.
  TYPE(t_vankaPointerNavSt2D), INTENT(INOUT) :: rvanka
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: ielementdist
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementList
    TYPE(t_elementDistribution), POINTER :: p_relementDistrU
    TYPE(t_elementDistribution), POINTER :: p_relementDistrV
    TYPE(t_elementDistribution), POINTER :: p_relementDistrP
    
    ! Loop through the element distributions of the velocity.
    DO ielementdist = 1,rvanka%p_rspatialDiscrU%inumFESpaces
    
      ! Get the corresponding element distributions of U, V and P.
      p_relementDistrU => &
          rvanka%p_rspatialDiscrU%RelementDistr(ielementdist)
      p_relementDistrV => &
          rvanka%p_rspatialDiscrV%RelementDistr(ielementdist)
      
      ! Either the same element for P everywhere, or there must be given one
      ! element distribution in the pressure for every velocity element distribution.
      IF (rvanka%p_rspatialDiscrP%inumFESpaces .GT. 1) THEN
        p_relementDistrP => &
            rvanka%p_rspatialDiscrP%RelementDistr(ielementdist)
      ELSE
        p_relementDistrP => &
            rvanka%p_rspatialDiscrP%RelementDistr(1)
      END IF
      
      ! Get the list of the elements to process.
      ! We take the element list of the X-velocity as 'primary' element list
      ! and assume that it coincides to that of the Y-velocity (and to that
      ! of the pressure).
      CALL storage_getbase_int (p_relementDistrU%h_IelementList,p_IelementList)
      
      ! Which element combination do we have now?
      IF ((elem_getPrimaryElement(p_relementDistrU%celement) .EQ. EL_Q1T) .AND. &
          (elem_getPrimaryElement(p_relementDistrV%celement) .EQ. EL_Q1T) .AND. &
          (elem_getPrimaryElement(p_relementDistrP%celement) .EQ. EL_Q0)) THEN
        ! Q1~/Q1~/Q0/Q1 discretisation
        
        ! Which VANKA subtype do we have? The diagonal VANKA of the full VANKA?
        SELECT CASE (csubtype)
        CASE (VANKATP_NAVST2D_DIAG)
          ! Call the jacobi-style vanka
          CALL vanka_NS2D_Q1TQ0_js(rvanka, &
                  rvector, rrhs, domega, p_IelementList)

        CASE (VANKATP_NAVST2D_FULL)
          IF (.NOT. ASSOCIATED(rvanka%p_DA12)) THEN
            ! Call the block-diagonal vanka
            CALL vanka_NS2D_Q1TQ0_bd(rvanka, &
                  rvector, rrhs, domega, p_IelementList)
          ELSE
            ! Call the fully coupled vanka
            CALL vanka_NS2D_Q1TQ0_fc(rvanka, &
                  rvector, rrhs, domega, p_IelementList)
          END IF
        
        CASE DEFAULT
          CALL output_line ('Unknown VANKA subtype!',&
              OU_CLASS_ERROR,OU_MODE_STD,'vanka_NavierStokes2D')
          CALL sys_halt()  
        
        END SELECT
    
      ELSE
        CALL output_line ('Unsupported discretisation!',&
            OU_CLASS_ERROR,OU_MODE_STD,'vanka_NavierStokes2D')
        CALL sys_halt()
        
      END IF
      
    END DO
      
  END SUBROUTINE

  ! ***************************************************************************

  SUBROUTINE vanka_NS2D_Q1TQ0_js(rvanka, rvector, rrhs, domega, IelementList)
  
!<description>
  ! This routine applies the specialised diagonal VANKA algorithm for
  ! 2D Navier-Stokes problems with Q1~/Q01 discretisation to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
!</description>

!<input>
  ! t_vankaPointerNavSt2D structure that saves algorithm-specific parameters.
  TYPE(t_vankaPointerNavSt2D), INTENT(IN) :: rvanka

  ! The right-hand-side vector of the system
  TYPE(t_vectorBlock), INTENT(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  REAL(DP), INTENT(IN)                    :: domega
  
  ! A list of element numbers where VANKA should be applied to.
  INTEGER, DIMENSION(:), INTENT(IN)       :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  TYPE(t_vectorBlock), INTENT(IN)         :: rvector
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
  real(DP), dimension(ndof) :: Du, Df
  real(DP), dimension(ndofV) :: Da1, Da2
  real(DP), dimension(ndofV,ndofP) :: Db1,Db2
  real(DP), dimension(ndofP,ndofV) :: Dd1, Dd2
  real(DP), dimension(ndofP) :: Dc
  
  ! temporary vectors
  real(DP), dimension(ndofV) :: Dt1, Dt2
  real(DP) :: dS
  
  ! Multiplication factors
  real(DP), dimension(3,3) :: Dmult
  
  ! Quick access for the matrix arrays
  integer, dimension(:), pointer :: p_KldA,p_KldA12,p_KldB,p_KldC,&
      p_KcolA,p_KcolA12,p_KcolB,p_KcolC,p_KdiagA,p_KdiagC
  real(DP), dimension(:), pointer :: p_DA11,p_DA12,p_DA21,p_DA22,p_DB1,p_DB2,&
      p_DD1,p_DD2,p_DC
    
  ! local variables
  integer :: i,j,iel,ielidx,i1,i2,k,l,o
  real(DP) :: daux
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
          Df(i)       = p_DrhsU(IdofV(i))   ! f_u
          Df(ndofV+i) = p_DrhsV(IdofV(i))   ! f_v
        end do
        o = 2*ndofV
        do i = 1, ndofP
          Df(o+i) = p_DrhsP(IdofP(i))       ! f_p
        end do
        
        ! Now fetch the main diagonal of A
        do i = 1, ndofV
          j = p_KdiagA(IdofV(i))
          Da1(i) = Dmult(1,1)*p_DA11(j)
          Da2(i) = Dmult(2,2)*p_DA22(j)
        end do
        
        ! Let's update the local RHS vector by subtracting A*u from it:
        ! f_u := f_u - A11*u
        ! f_v := f_v - A22*v
        do k = 1, ndofV
          i1 = p_KldA(IdofV(k))
          i2 = p_KldA(IdofV(k)+1)-1
          do i = i1, i2
            j = p_KcolA(i)
            Df(k)       = Df(k)       - Dmult(1,1)*p_DA11(i)*p_DvecU(j)
            Df(ndofV+k) = Df(ndofV+k) - Dmult(2,2)*p_DA22(i)*p_DvecV(j)
          end do
        end do
        
        ! Now we also need to subtract B*p from our RHS, and by the same time,
        ! we also have to subtract D*u - and we will build the local B/D matrices.
        ! f_u := f_u - B1*p
        ! f_v := f_v - B2*p
        ! f_p := f_p - D1*u - D2*v
        do k = 1, ndofV
          i1 = p_KldB(IdofV(k))
          i2 = p_KldB(IdofV(k)+1)-1
          do i = i1, i2
            j = p_KcolB(i)
            daux = p_DvecP(j)
            Df(k)       = Df(k)       - Dmult(1,3)*p_DB1(i)*daux
            Df(ndofV+k) = Df(ndofV+k) - Dmult(2,3)*p_DB2(i)*daux
            do l = 1, ndofP
              if(j .ne. IdofP(l)) cycle
              
              ! Get the local entries of the B/D matrices
              Db1(k,l) = Dmult(1,3)*p_DB1(i)
              Db2(k,l) = Dmult(2,3)*p_DB2(i)
              Dd1(l,k) = Dmult(3,1)*p_DD1(i)
              Dd2(l,k) = Dmult(3,2)*p_DD2(i)
              
              ! Update f_p
              Df(2*ndofV+l) = Df(2*ndofV+l) - Dd1(l,k)*p_DvecU(IdofV(k)) &
                                            - Dd2(l,k)*p_DvecV(IdofV(k))
            end do
          end do
        end do
        
        ! Invert A1 and A2
        Da1(1) = 1.0_DP / Da1(1)
        Da1(2) = 1.0_DP / Da1(2)
        Da1(3) = 1.0_DP / Da1(3)
        Da1(4) = 1.0_DP / Da1(4)
        Da2(1) = 1.0_DP / Da2(1)
        Da2(2) = 1.0_DP / Da2(2)
        Da2(3) = 1.0_DP / Da2(3)
        Da2(4) = 1.0_DP / Da2(4)
        
        ! Precalculate D * A^-1
        Dt1(1) = Dd1(1,1)*Da1(1)
        Dt1(2) = Dd1(1,2)*Da1(2)
        Dt1(3) = Dd1(1,3)*Da1(3)
        Dt1(4) = Dd1(1,4)*Da1(4)
        Dt2(1) = Dd2(1,1)*Da2(1)
        Dt2(2) = Dd2(1,2)*Da2(2)
        Dt2(3) = Dd2(1,3)*Da2(3)
        Dt2(4) = Dd2(1,4)*Da2(4)

        ! Calculate Schur-Complement of A
        ! S := D * A^-1 * B 
        dS = Dt1(1)*Db1(1,1)+Dt1(2)*Db1(2,1)+Dt1(3)*Db1(3,1)+Dt1(4)*Db1(4,1) &
           + Dt2(1)*Db2(1,1)+Dt2(2)*Db2(2,1)+Dt2(3)*Db2(3,1)+Dt2(4)*Db2(4,1)
        
        ! Calculate pressure
        ! p := D * A^-1 * f_u - f_p
        Du(9) = (-Df(9) &
           + Dt1(1)*Df(1)+Dt1(2)*Df(2)+Dt1(3)*Df(3)+Dt1(4)*Df(4) &
           + Dt2(1)*Df(5)+Dt2(2)*Df(6)+Dt2(3)*Df(7)+Dt2(4)*Df(7)) / dS
        
        ! Calculate X- and Y-velocity
        ! u := A^-1 * (f_u - B * p)
        Du(1) = Da1(1)*(Df(1) - Db1(1,1)*Du(9))
        Du(2) = Da1(2)*(Df(2) - Db1(2,1)*Du(9))
        Du(3) = Da1(3)*(Df(3) - Db1(3,1)*Du(9))
        Du(4) = Da1(4)*(Df(4) - Db1(4,1)*Du(9))
        Du(5) = Da2(1)*(Df(5) - Db2(1,1)*Du(9))
        Du(6) = Da2(2)*(Df(6) - Db2(2,1)*Du(9))
        Du(7) = Da2(3)*(Df(7) - Db2(3,1)*Du(9))
        Du(8) = Da2(4)*(Df(8) - Db2(4,1)*Du(9))

        ! Incorporate our local solution into the global one.
        do i = 1, ndofV
          j = IdofV(i)
          p_DvecU(j) = p_DvecU(j) + domega * Du(i)
          p_DvecV(j) = p_DvecV(j) + domega * Du(ndofV+i)
        end do
        o = 2*ndofV
        do i = 1, ndofP
          j = IdofP(i)
          p_DvecP(j) = p_DvecP(j) + domega * Du(o+i)
        end do
      
      end do ! ielidx
      
    else if(bHaveA12 .and. (.not. bHaveC)) then
      ! Only A12/A21 present
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
        
        ! Now fetch the main diagonal of A
        do i = 1, ndofV
          j = p_KdiagA(IdofV(i))
          Da1(i) = Dmult(1,1)*p_DA11(j)
          Da2(i) = Dmult(2,2)*p_DA22(j)
        end do
        
        ! Let's update the local RHS vector by subtracting A*u from it:
        ! f_u := f_u - A11*u
        ! f_v := f_v - A22*v
        do k = 1, ndofV
          i1 = p_KldA(IdofV(k))
          i2 = p_KldA(IdofV(k)+1)-1
          do i = i1, i2
            j = p_KcolA(i)
            Df(k)       = Df(k)       - Dmult(1,1)*p_DA11(i)*p_DvecU(j)
            Df(ndofV+k) = Df(ndofV+k) - Dmult(2,2)*p_DA22(i)*p_DvecV(j)
          end do
        end do
        
        ! f_u := f_u - A12*v
        ! f_v := f_v - A21*u
        do k = 1, ndofV
          i1 = p_KldA12(IdofV(k))
          i2 = p_KldA12(IdofV(k)+1)-1
          do i = i1, i2
            j = p_KcolA12(i)
            Df(k)       = Df(k)       - Dmult(1,2)*p_DA12(i)*p_DvecV(j)
            Df(ndofV+k) = Df(ndofV+k) - Dmult(2,1)*p_DA21(i)*p_DvecU(j)
          end do
        end do
        
        ! Now we also need to subtract B*p from our RHS, and by the same time,
        ! we also have to subtract D*u - and we will build the local B/D matrices.
        ! f_u := f_u - B1*p
        ! f_v := f_v - B2*p
        ! f_p := f_p - D1*u - D2*v
        do k = 1, ndofV
          i1 = p_KldB(IdofV(k))
          i2 = p_KldB(IdofV(k)+1)-1
          do i = i1, i2
            j = p_KcolB(i)
            daux = p_DvecP(j)
            Df(k)       = Df(k)       - Dmult(1,3)*p_DB1(i)*daux
            Df(ndofV+k) = Df(ndofV+k) - Dmult(2,3)*p_DB2(i)*daux
            do l = 1, ndofP
              if(j .ne. IdofP(l)) cycle
              
              ! Get the local entries of the B/D matrices
              Db1(k,l) = Dmult(1,3)*p_DB1(i)
              Db2(k,l) = Dmult(2,3)*p_DB2(i)
              Dd1(l,k) = Dmult(3,1)*p_DD1(i)
              Dd2(l,k) = Dmult(3,2)*p_DD2(i)
              
              ! Update f_p
              Df(2*ndofV+l) = Df(2*ndofV+l) - Dd1(l,k)*p_DvecU(IdofV(k)) &
                                            - Dd2(l,k)*p_DvecV(IdofV(k))
            end do
          end do
        end do
        
        ! Invert A1 and A2
        Da1(1) = 1.0_DP / Da1(1)
        Da1(2) = 1.0_DP / Da1(2)
        Da1(3) = 1.0_DP / Da1(3)
        Da1(4) = 1.0_DP / Da1(4)
        Da2(1) = 1.0_DP / Da2(1)
        Da2(2) = 1.0_DP / Da2(2)
        Da2(3) = 1.0_DP / Da2(3)
        Da2(4) = 1.0_DP / Da2(4)
        
        ! Precalculate D * A^-1
        Dt1(1) = Dd1(1,1)*Da1(1)
        Dt1(2) = Dd1(1,2)*Da1(2)
        Dt1(3) = Dd1(1,3)*Da1(3)
        Dt1(4) = Dd1(1,4)*Da1(4)
        Dt2(1) = Dd2(1,1)*Da2(1)
        Dt2(2) = Dd2(1,2)*Da2(2)
        Dt2(3) = Dd2(1,3)*Da2(3)
        Dt2(4) = Dd2(1,4)*Da2(4)

        ! Calculate Schur-Complement of A
        ! S := D * A^-1 * B 
        dS = Dt1(1)*Db1(1,1)+Dt1(2)*Db1(2,1)+Dt1(3)*Db1(3,1)+Dt1(4)*Db1(4,1) &
           + Dt2(1)*Db2(1,1)+Dt2(2)*Db2(2,1)+Dt2(3)*Db2(3,1)+Dt2(4)*Db2(4,1)
        
        ! Calculate pressure
        ! p := D * A^-1 * f_u - f_p
        Du(9) = (-Df(9) &
           + Dt1(1)*Df(1)+Dt1(2)*Df(2)+Dt1(3)*Df(3)+Dt1(4)*Df(4) &
           + Dt2(1)*Df(5)+Dt2(2)*Df(6)+Dt2(3)*Df(7)+Dt2(4)*Df(7)) / dS
        
        ! Calculate X- and Y-velocity
        ! u := A^-1 * (f_u - B * p)
        Du(1) = Da1(1)*(Df(1) - Db1(1,1)*Du(9))
        Du(2) = Da1(2)*(Df(2) - Db1(2,1)*Du(9))
        Du(3) = Da1(3)*(Df(3) - Db1(3,1)*Du(9))
        Du(4) = Da1(4)*(Df(4) - Db1(4,1)*Du(9))
        Du(5) = Da2(1)*(Df(5) - Db2(1,1)*Du(9))
        Du(6) = Da2(2)*(Df(6) - Db2(2,1)*Du(9))
        Du(7) = Da2(3)*(Df(7) - Db2(3,1)*Du(9))
        Du(8) = Da2(4)*(Df(8) - Db2(4,1)*Du(9))

        ! Incorporate our local solution into the global one.
        do i = 1, ndofV
          j = IdofV(i)
          p_DvecU(j) = p_DvecU(j) + domega * Du(i)
          p_DvecV(j) = p_DvecV(j) + domega * Du(ndofV+i)
        end do
        o = 2*ndofV
        do i = 1, ndofP
          j = IdofP(i)
          p_DvecP(j) = p_DvecP(j) + domega * Du(o+i)
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
          Df(i)       = p_DrhsU(IdofV(i))   ! f_u
          Df(ndofV+i) = p_DrhsV(IdofV(i))   ! f_v
        end do
        o = 2*ndofV
        do i = 1, ndofP
          Df(o+i) = p_DrhsP(IdofP(i))       ! f_p
        end do
        
        ! Now fetch the main diagonal of A
        do i = 1, ndofV
          j = p_KdiagA(IdofV(i))
          Da1(i) = Dmult(1,1)*p_DA11(j)
          Da2(i) = Dmult(2,2)*p_DA22(j)
        end do
        
        ! Let's update the local RHS vector by subtracting A*u from it:
        ! f_u := f_u - A11*u
        ! f_v := f_v - A22*v
        do k = 1, ndofV
          i1 = p_KldA(IdofV(k))
          i2 = p_KldA(IdofV(k)+1)-1
          do i = i1, i2
            j = p_KcolA(i)
            Df(k)       = Df(k)       - Dmult(1,1)*p_DA11(i)*p_DvecU(j)
            Df(ndofV+k) = Df(ndofV+k) - Dmult(2,2)*p_DA22(i)*p_DvecV(j)
          end do
        end do
        
        ! What about A12/A21?
        if(bHaveA12) then
          ! f_u := f_u - A12*v
          ! f_v := f_v - A21*u
          do k = 1, ndofV
            i1 = p_KldA12(IdofV(k))
            i2 = p_KldA12(IdofV(k)+1)-1
            do i = i1, i2
              j = p_KcolA12(i)
              Df(k)       = Df(k)       - Dmult(1,2)*p_DA12(i)*p_DvecV(j)
              Df(ndofV+k) = Df(ndofV+k) - Dmult(2,1)*p_DA21(i)*p_DvecU(j)
            end do
          end do
        end if
        
        ! Now we also need to subtract B*p from our RHS, and by the same time,
        ! we also have to subtract D*u - and we will build the local B/D matrices.
        ! f_u := f_u - B1*p
        ! f_v := f_v - B2*p
        ! f_p := f_p - D1*u - D2*v
        do k = 1, ndofV
          i1 = p_KldB(IdofV(k))
          i2 = p_KldB(IdofV(k)+1)-1
          do i = i1, i2
            j = p_KcolB(i)
            daux = p_DvecP(j)
            Df(k)       = Df(k)       - Dmult(1,3)*p_DB1(i)*daux
            Df(ndofV+k) = Df(ndofV+k) - Dmult(2,3)*p_DB2(i)*daux
            do l = 1, ndofP
              if(j .ne. IdofP(l)) cycle
              
              ! Get the local entries of the B/D matrices
              Db1(k,l) = Dmult(1,3)*p_DB1(i)
              Db2(k,l) = Dmult(2,3)*p_DB2(i)
              Dd1(l,k) = Dmult(3,1)*p_DD1(i)
              Dd2(l,k) = Dmult(3,2)*p_DD2(i)
              
              ! Update f_p
              Df(2*ndofV+l) = Df(2*ndofV+l) - Dd1(l,k)*p_DvecU(IdofV(k)) &
                                            - Dd2(l,k)*p_DvecV(IdofV(k))
            end do
          end do
        end do
        
        ! Do we have a C-matrix?
        if(bHaveC) then
          ! Yes, so update the local RHS
          ! f_p := f_p - C*p
          o = 2*ndofV
          do k = 1, ndofP
            i1 = p_KldC(IdofP(k))
            i2 = p_KldC(IdofP(k)+1)-1
            do i = i1, i2
              Df(o+k) = Df(o+k) - Dmult(3,3)*p_DC(i)*p_DvecP(p_KcolC(i))
            end do
            ! Get the main diagonal entry of C
            Dc(k) = Dmult(3,3)*p_DC(p_KdiagC(IdofP(k)))
          end do
        end if
        
        ! Invert A1 and A2
        Da1(1) = 1.0_DP / Da1(1)
        Da1(2) = 1.0_DP / Da1(2)
        Da1(3) = 1.0_DP / Da1(3)
        Da1(4) = 1.0_DP / Da1(4)
        Da2(1) = 1.0_DP / Da2(1)
        Da2(2) = 1.0_DP / Da2(2)
        Da2(3) = 1.0_DP / Da2(3)
        Da2(4) = 1.0_DP / Da2(4)
        
        ! Precalculate D * A^-1
        Dt1(1) = Dd1(1,1)*Da1(1)
        Dt1(2) = Dd1(1,2)*Da1(2)
        Dt1(3) = Dd1(1,3)*Da1(3)
        Dt1(4) = Dd1(1,4)*Da1(4)
        Dt2(1) = Dd2(1,1)*Da2(1)
        Dt2(2) = Dd2(1,2)*Da2(2)
        Dt2(3) = Dd2(1,3)*Da2(3)
        Dt2(4) = Dd2(1,4)*Da2(4)

        ! Calculate Schur-Complement of A
        ! S := -C + D * A^-1 * B 
        dS = -Dc(1) &
           + Dt1(1)*Db1(1,1)+Dt1(2)*Db1(2,1)+Dt1(3)*Db1(3,1)+Dt1(4)*Db1(4,1) &
           + Dt2(1)*Db2(1,1)+Dt2(2)*Db2(2,1)+Dt2(3)*Db2(3,1)+Dt2(4)*Db2(4,1)
        
        ! Calculate pressure
        ! p := D * A^-1 * f_u - f_p
        Du(9) = (-Df(9) &
           + Dt1(1)*Df(1)+Dt1(2)*Df(2)+Dt1(3)*Df(3)+Dt1(4)*Df(4) &
           + Dt2(1)*Df(5)+Dt2(2)*Df(6)+Dt2(3)*Df(7)+Dt2(4)*Df(7)) / dS
        
        ! Calculate X- and Y-velocity
        ! u := A^-1 * (f_u - B * p)
        Du(1) = Da1(1)*(Df(1) - Db1(1,1)*Du(9))
        Du(2) = Da1(2)*(Df(2) - Db1(2,1)*Du(9))
        Du(3) = Da1(3)*(Df(3) - Db1(3,1)*Du(9))
        Du(4) = Da1(4)*(Df(4) - Db1(4,1)*Du(9))
        Du(5) = Da2(1)*(Df(5) - Db2(1,1)*Du(9))
        Du(6) = Da2(2)*(Df(6) - Db2(2,1)*Du(9))
        Du(7) = Da2(3)*(Df(7) - Db2(3,1)*Du(9))
        Du(8) = Da2(4)*(Df(8) - Db2(4,1)*Du(9))

        ! Incorporate our local solution into the global one.
        do i = 1, ndofV
          j = IdofV(i)
          p_DvecU(j) = p_DvecU(j) + domega * Du(i)
          p_DvecV(j) = p_DvecV(j) + domega * Du(ndofV+i)
        end do
        o = 2*ndofV
        do i = 1, ndofP
          j = IdofP(i)
          p_DvecP(j) = p_DvecP(j) + domega * Du(o+i)
        end do
      
      end do ! ielidx
    end if

  end subroutine

  ! ***************************************************************************
  
  SUBROUTINE vanka_NS2D_Q1TQ0_bd(rvanka, rvector, rrhs, domega, IelementList)
  
!<description>
  ! This routine applies the specialised diagonal VANKA algorithm for
  ! 2D Navier-Stokes problems with Q1~/Q0 discretisation to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
!</description>

!<input>
  ! t_vankaPointerNavSt2D structure that saves algorithm-specific parameters.
  TYPE(t_vankaPointerNavSt2D), INTENT(IN) :: rvanka

  ! The right-hand-side vector of the system
  TYPE(t_vectorBlock), INTENT(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  REAL(DP), INTENT(IN)                    :: domega
  
  ! A list of element numbers where VANKA should be applied to.
  INTEGER, DIMENSION(:), INTENT(IN)       :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  TYPE(t_vectorBlock), INTENT(IN)         :: rvector
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
  real(DP), dimension(ndof) :: Du, Df
  real(DP), dimension(ndofV,ndofV) :: Da1, Da2
  real(DP), dimension(ndofV,ndofP) :: Db1,Db2
  real(DP), dimension(ndofP,ndofV) :: Dd1, Dd2
  real(DP), dimension(ndofP,ndofP) :: Dc
  
  ! Local variables
  real(DP), dimension(ndofV,ndofV) :: Di1,Di2
  real(DP), dimension(ndofV) :: Dg1, Dg2, Dt1, Dt2
  real(DP) :: Ds
  
  ! Multiplication factors
  real(DP), dimension(3,3) :: Dmult
  
  ! Quick access for the matrix arrays
  integer, dimension(:), pointer :: p_KldA,p_KldB,p_KldC,&
      p_KcolA,p_KcolB,p_KcolC,p_KdiagA,p_KdiagC
  real(DP), dimension(:), pointer :: p_DA11,p_DA22,p_DB1,p_DB2,&
      p_DD1,p_DD2,p_DC
  
  ! local variables
  integer :: i,j,iel,ielidx,i1,i2,k,l,o
  real(DP) :: daux
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
    
    ! Clear the optional matrices
    Dc = 0.0_DP
    
    if(bHaveC) then
      ! General case
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
        
        ! Let's update the local RHS vector by subtracting A*u from it:
        ! f_u := f_u - A11*u
        ! f_v := f_v - A22*v
        do k = 1, ndofV
          i1 = p_KldA(IdofV(k))
          i2 = p_KldA(IdofV(k)+1)-1
          do i = i1, i2
            j = p_KcolA(i)
            Df(k)       = Df(k)       - Dmult(1,1)*p_DA11(i)*p_DvecU(j)
            Df(ndofV+k) = Df(ndofV+k) - Dmult(2,2)*p_DA22(i)*p_DvecV(j)
            do l = 1, ndofV
              if(j .ne. IdofV(l)) cycle
              
              ! Get the local entries of the A1/A2 matrices
              Da1(k,l) = Dmult(1,1)*p_DA11(i)
              Da2(k,l) = Dmult(2,2)*p_DA22(i)
            end do
          end do
        end do
        
        ! Now we also need to subtract B*p from our RHS, and by the same time,
        ! we also have to subtract D*u - and we will build the local B/D matrices.
        ! f_u := f_u - B1*p
        ! f_v := f_v - B2*p
        ! f_p := f_p - D1*u - D2*v
        do k = 1, ndofV
          i1 = p_KldB(IdofV(k))
          i2 = p_KldB(IdofV(k)+1)-1
          do i = i1, i2
            j = p_KcolB(i)
            daux = p_DvecP(j)
            Df(k)       = Df(k)       - Dmult(1,3)*p_DB1(i)*daux
            Df(ndofV+k) = Df(ndofV+k) - Dmult(2,3)*p_DB2(i)*daux
            do l = 1, ndofP
              if(j .ne. IdofP(l)) cycle
              
              ! Get the local entries of the B/D matrices
              Db1(k,l) = Dmult(1,3)*p_DB1(i)
              Db2(k,l) = Dmult(2,3)*p_DB2(i)
              Dd1(l,k) = Dmult(3,1)*p_DD1(i)
              Dd2(l,k) = Dmult(3,2)*p_DD2(i)
              
              ! Update f_p
              Df(2*ndofV+l) = Df(2*ndofV+l) - Dd1(l,k)*p_DvecU(IdofV(k)) &
                                            - Dd2(l,k)*p_DvecV(IdofV(k))
            end do
          end do
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
        dS = Dt1(1)*Db1(1,1)+Dt1(2)*Db1(2,1)+Dt1(3)*Db1(3,1)+Dt1(4)*Db1(4,1)&
           + Dt2(1)*Db2(1,1)+Dt2(2)*Db2(2,1)+Dt2(3)*Db2(3,1)+Dt2(4)*Db2(4,1)
        
        ! Calculate pressure
        ! p := D * A^-1 * f_u - f_p
        Du(9) = (-Df(9) &
              + Dt1(1)*Df(1)+Dt1(2)*Df(2)+Dt1(3)*Df(3)+Dt1(4)*Df(4) &
              + Dt2(1)*Df(5)+Dt2(2)*Df(6)+Dt2(3)*Df(7)+Dt2(4)*Df(8)) / dS

        ! Update RHS
        ! g_u := f_u - B * p
        Dg1(1) = Df(1) - Db1(1,1)*Du(9)
        Dg1(2) = Df(2) - Db1(2,1)*Du(9)
        Dg1(3) = Df(3) - Db1(3,1)*Du(9)
        Dg1(4) = Df(4) - Db1(4,1)*Du(9)
        Dg2(1) = Df(5) - Db2(1,1)*Du(9)
        Dg2(2) = Df(6) - Db2(2,1)*Du(9)
        Dg2(3) = Df(7) - Db2(3,1)*Du(9)
        Dg2(4) = Df(8) - Db2(4,1)*Du(9)
        
        ! Calculate X- and Y-velocity
        ! u := A^-1 * g_u
        Du(1) = Di1(1,1)*Dg1(1)+Di1(1,2)*Dg1(2)+Di1(1,3)*Dg1(3)+Di1(1,4)*Dg1(4)
        Du(2) = Di1(2,1)*Dg1(1)+Di1(2,2)*Dg1(2)+Di1(2,3)*Dg1(3)+Di1(2,4)*Dg1(4)
        Du(3) = Di1(3,1)*Dg1(1)+Di1(3,2)*Dg1(2)+Di1(3,3)*Dg1(3)+Di1(3,4)*Dg1(4)
        Du(4) = Di1(4,1)*Dg1(1)+Di1(4,2)*Dg1(2)+Di1(4,3)*Dg1(3)+Di1(4,4)*Dg1(4)
        Du(5) = Di2(1,1)*Dg2(1)+Di2(1,2)*Dg2(2)+Di2(1,3)*Dg2(3)+Di2(1,4)*Dg2(4)
        Du(6) = Di2(2,1)*Dg2(1)+Di2(2,2)*Dg2(2)+Di2(2,3)*Dg2(3)+Di2(2,4)*Dg2(4)
        Du(7) = Di2(3,1)*Dg2(1)+Di2(3,2)*Dg2(2)+Di2(3,3)*Dg2(3)+Di2(3,4)*Dg2(4)
        Du(8) = Di2(4,1)*Dg2(1)+Di2(4,2)*Dg2(2)+Di2(4,3)*Dg2(3)+Di2(4,4)*Dg2(4)
        
        ! Incorporate our local solution into the global one.
        do i = 1, ndofV
          j = IdofV(i)
          p_DvecU(j) = p_DvecU(j) + domega * Du(i)
          p_DvecV(j) = p_DvecV(j) + domega * Du(ndofV+i)
        end do
        o = 2*ndofV
        do i = 1, ndofP
          j = IdofP(i)
          p_DvecP(j) = p_DvecP(j) + domega * Du(o+i)
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
          Df(i)       = p_DrhsU(IdofV(i))   ! f_u
          Df(ndofV+i) = p_DrhsV(IdofV(i))   ! f_v
        end do
        o = 2*ndofV
        do i = 1, ndofP
          Df(o+i) = p_DrhsP(IdofP(i))       ! f_p
        end do
        
        ! Let's update the local RHS vector by subtracting A*u from it:
        ! f_u := f_u - A11*u
        ! f_v := f_v - A22*v
        do k = 1, ndofV
          i1 = p_KldA(IdofV(k))
          i2 = p_KldA(IdofV(k)+1)-1
          do i = i1, i2
            j = p_KcolA(i)
            Df(k)       = Df(k)       - Dmult(1,1)*p_DA11(i)*p_DvecU(j)
            Df(ndofV+k) = Df(ndofV+k) - Dmult(2,2)*p_DA22(i)*p_DvecV(j)
            do l = 1, ndofV
              if(j .ne. IdofV(l)) cycle
              
              ! Get the local entries of the A1/A2 matrices
              Da1(k,l) = Dmult(1,1)*p_DA11(i)
              Da2(k,l) = Dmult(2,2)*p_DA22(i)
            end do
          end do
        end do
        
        ! Now we also need to subtract B*p from our RHS, and by the same time,
        ! we also have to subtract D*u - and we will build the local B/D matrices.
        ! f_u := f_u - B1*p
        ! f_v := f_v - B2*p
        ! f_p := f_p - D1*u - D2*v
        do k = 1, ndofV
          i1 = p_KldB(IdofV(k))
          i2 = p_KldB(IdofV(k)+1)-1
          do i = i1, i2
            j = p_KcolB(i)
            daux = p_DvecP(j)
            Df(k)       = Df(k)       - Dmult(1,3)*p_DB1(i)*daux
            Df(ndofV+k) = Df(ndofV+k) - Dmult(2,3)*p_DB2(i)*daux
            do l = 1, ndofP
              if(j .ne. IdofP(l)) cycle
              
              ! Get the local entries of the B/D matrices
              Db1(k,l) = Dmult(1,3)*p_DB1(i)
              Db2(k,l) = Dmult(2,3)*p_DB2(i)
              Dd1(l,k) = Dmult(3,1)*p_DD1(i)
              Dd2(l,k) = Dmult(3,2)*p_DD2(i)
              
              ! Update f_p
              Df(2*ndofV+l) = Df(2*ndofV+l) - Dd1(l,k)*p_DvecU(IdofV(k)) &
                                            - Dd2(l,k)*p_DvecV(IdofV(k))
            end do
          end do
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
              if(j .ne. IdofP(l)) cycle
              
              ! Get the local entries of the C matrix
              Dc(k,l) = Dmult(3,3)*p_DC(i)
            end do
          end do
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
        dS = -Dc(1,1)+Dt1(1)*Db1(1,1)+Dt1(2)*Db1(2,1)+Dt1(3)*Db1(3,1)+Dt1(4)*Db1(4,1)&
                     +Dt2(1)*Db2(1,1)+Dt2(2)*Db2(2,1)+Dt2(3)*Db2(3,1)+Dt2(4)*Db2(4,1)
        
        ! Calculate pressure
        ! p := D * A^-1 * f_u - f_p
        Du(9) = (-Df(9) &
              + Dt1(1)*Df(1)+Dt1(2)*Df(2)+Dt1(3)*Df(3)+Dt1(4)*Df(4) &
              + Dt2(1)*Df(5)+Dt2(2)*Df(6)+Dt2(3)*Df(7)+Dt2(4)*Df(8)) / dS

        ! Update RHS
        ! g_u := f_u - B * p
        Dg1(1) = Df(1) - Db1(1,1)*Du(9)
        Dg1(2) = Df(2) - Db1(2,1)*Du(9)
        Dg1(3) = Df(3) - Db1(3,1)*Du(9)
        Dg1(4) = Df(4) - Db1(4,1)*Du(9)
        Dg2(1) = Df(5) - Db2(1,1)*Du(9)
        Dg2(2) = Df(6) - Db2(2,1)*Du(9)
        Dg2(3) = Df(7) - Db2(3,1)*Du(9)
        Dg2(4) = Df(8) - Db2(4,1)*Du(9)
        
        ! Calculate X- and Y-velocity
        ! u := A^-1 * g_u
        Du(1) = Di1(1,1)*Dg1(1)+Di1(1,2)*Dg1(2)+Di1(1,3)*Dg1(3)+Di1(1,4)*Dg1(4)
        Du(2) = Di1(2,1)*Dg1(1)+Di1(2,2)*Dg1(2)+Di1(2,3)*Dg1(3)+Di1(2,4)*Dg1(4)
        Du(3) = Di1(3,1)*Dg1(1)+Di1(3,2)*Dg1(2)+Di1(3,3)*Dg1(3)+Di1(3,4)*Dg1(4)
        Du(4) = Di1(4,1)*Dg1(1)+Di1(4,2)*Dg1(2)+Di1(4,3)*Dg1(3)+Di1(4,4)*Dg1(4)
        Du(5) = Di2(1,1)*Dg2(1)+Di2(1,2)*Dg2(2)+Di2(1,3)*Dg2(3)+Di2(1,4)*Dg2(4)
        Du(6) = Di2(2,1)*Dg2(1)+Di2(2,2)*Dg2(2)+Di2(2,3)*Dg2(3)+Di2(2,4)*Dg2(4)
        Du(7) = Di2(3,1)*Dg2(1)+Di2(3,2)*Dg2(2)+Di2(3,3)*Dg2(3)+Di2(3,4)*Dg2(4)
        Du(8) = Di2(4,1)*Dg2(1)+Di2(4,2)*Dg2(2)+Di2(4,3)*Dg2(3)+Di2(4,4)*Dg2(4)
        
        ! Incorporate our local solution into the global one.
        do i = 1, ndofV
          j = IdofV(i)
          p_DvecU(j) = p_DvecU(j) + domega * Du(i)
          p_DvecV(j) = p_DvecV(j) + domega * Du(ndofV+i)
        end do
        o = 2*ndofV
        do i = 1, ndofP
          j = IdofP(i)
          p_DvecP(j) = p_DvecP(j) + domega * Du(o+i)
        end do
      
      end do ! ielidx
    end if

  END SUBROUTINE

  ! ***************************************************************************
  
  SUBROUTINE vanka_NS2D_Q1TQ0_fc(rvanka, rvector, rrhs, domega, IelementList)
  
!<description>
  ! This routine applies the specialised diagonal VANKA algorithm for
  ! 2D Navier-Stokes problems with Q1~/Q0 discretisation to the system $Ax=b$.
  ! x=rvector is the initial solution vector and b=rrhs the right-hand-side
  ! vector. The rvanka structure has to be initialised before calling
  ! this routine, as this holds a reference to the system matrix.
!</description>

!<input>
  ! t_vankaPointerNavSt2D structure that saves algorithm-specific parameters.
  TYPE(t_vankaPointerNavSt2D), INTENT(IN) :: rvanka

  ! The right-hand-side vector of the system
  TYPE(t_vectorBlock), INTENT(IN)         :: rrhs
  
  ! Relaxation parameter. Standard=1.0_DP.
  REAL(DP), INTENT(IN)                    :: domega
  
  ! A list of element numbers where VANKA should be applied to.
  INTEGER, DIMENSION(:), INTENT(IN)       :: IelementList
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  TYPE(t_vectorBlock), INTENT(IN)         :: rvector
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
  integer, dimension(:), pointer :: p_KldA,p_KldA12,p_KldB,p_KldC,&
      p_KcolA,p_KcolA12,p_KcolB,p_KcolC
  real(DP), dimension(:), pointer :: p_DA11,p_DA12,p_DA21,p_DA22,p_DB1,p_DB2,&
      p_DD1,p_DD2,p_DC
  
  ! local variables
  integer :: i,j,iel,ielidx,i1,i2,k,l,o,p
  real(DP) :: daux
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
          do i = i1, i2
            j = p_KcolA(i)
            Df(k)       = Df(k)       - Dmult(1,1)*p_DA11(i)*p_DvecU(j)
            Df(ndofV+k) = Df(ndofV+k) - Dmult(2,2)*p_DA22(i)*p_DvecV(j)
            do l = 1, ndofV
              if(j .ne. IdofV(l)) cycle
              
              ! Get the local entries of the A matrices
              Da(  k,  l) = Dmult(1,1)*p_DA11(i)
              Da(p+k,p+l) = Dmult(2,2)*p_DA22(i)
            end do
          end do
        end do

        ! f_u := f_u - A12*v
        ! f_v := f_v - A21*u
        do k = 1, ndofV
          i1 = p_KldA12(IdofV(k))
          i2 = p_KldA12(IdofV(k)+1)-1
          do i = i1, i2
            j = p_KcolA12(i)
            Df(k)       = Df(k)       - Dmult(1,2)*p_DA12(i)*p_DvecV(j)
            Df(ndofV+k) = Df(ndofV+k) - Dmult(2,1)*p_DA21(i)*p_DvecU(j)
            do l = 1, ndofV
              if(j .ne. IdofV(l)) cycle
              
              ! Get the local entries of the A matrices
              Da(  k,p+l) = Dmult(1,2)*p_DA12(i)
              Da(p+k,  l) = Dmult(2,1)*p_DA21(i)
            end do
          end do
        end do
        
        ! Now we also need to subtract B*p from our RHS, and by the same time,
        ! we also have to subtract D*u - and we will build the local B/D matrices.
        ! f_u := f_u - B1*p
        ! f_v := f_v - B2*p
        ! f_p := f_p - D1*u - D2*v
        o = ndofV
        p = 2*ndofV
        do k = 1, ndofV
          i1 = p_KldB(IdofV(k))
          i2 = p_KldB(IdofV(k)+1)-1
          do i = i1, i2
            j = p_KcolB(i)
            daux = p_DvecP(j)
            Df(k)       = Df(k)       - Dmult(1,3)*p_DB1(i)*daux
            Df(ndofV+k) = Df(ndofV+k) - Dmult(2,3)*p_DB2(i)*daux
            do l = 1, ndofP
              if(j .ne. IdofP(l)) cycle
              
              ! Get the local entries of the B/D matrices
              Da(  k,p+l) = Dmult(1,3)*p_DB1(i)
              Da(o+k,p+l) = Dmult(2,3)*p_DB2(i)
              Da(p+l,  k) = Dmult(3,1)*p_DD1(i)
              Da(p+l,o+k) = Dmult(3,2)*p_DD2(i)
              
              ! Update f_p
              Df(2*ndofV+l) = Df(2*ndofV+l) - Da(p+l,k)*p_DvecU(IdofV(k)) &
                                            - Da(p+l,k)*p_DvecV(IdofV(k))
            end do
          end do
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
              if(j .ne. IdofP(l)) cycle
              
              ! Get the local entries of the C matrix
              Da(o+k,o+l) = Dmult(3,3)*p_DC(i)
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
    
      ! C does not exist
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
          do i = i1, i2
            j = p_KcolA(i)
            Df(k)       = Df(k)       - Dmult(1,1)*p_DA11(i)*p_DvecU(j)
            Df(ndofV+k) = Df(ndofV+k) - Dmult(2,2)*p_DA22(i)*p_DvecV(j)
            do l = 1, ndofV
              if(j .ne. IdofV(l)) cycle
              
              ! Get the local entries of the A matrices
              Da(  k,  l) = Dmult(1,1)*p_DA11(i)
              Da(p+k,p+l) = Dmult(2,2)*p_DA22(i)
            end do
          end do
        end do

        ! f_u := f_u - A12*v
        ! f_v := f_v - A21*u
        do k = 1, ndofV
          i1 = p_KldA12(IdofV(k))
          i2 = p_KldA12(IdofV(k)+1)-1
          do i = i1, i2
            j = p_KcolA12(i)
            Df(k)       = Df(k)       - Dmult(1,2)*p_DA12(i)*p_DvecV(j)
            Df(ndofV+k) = Df(ndofV+k) - Dmult(2,1)*p_DA21(i)*p_DvecU(j)
            do l = 1, ndofV
              if(j .ne. IdofV(l)) cycle
              
              ! Get the local entries of the A matrices
              Da(  k,p+l) = Dmult(1,2)*p_DA12(i)
              Da(p+k,  l) = Dmult(2,1)*p_DA21(i)
            end do
          end do
        end do
        
        ! Now we also need to subtract B*p from our RHS, and by the same time,
        ! we also have to subtract D*u - and we will build the local B/D matrices.
        ! f_u := f_u - B1*p
        ! f_v := f_v - B2*p
        ! f_p := f_p - D1*u - D2*v
        o = ndofV
        p = 2*ndofV
        do k = 1, ndofV
          i1 = p_KldB(IdofV(k))
          i2 = p_KldB(IdofV(k)+1)-1
          do i = i1, i2
            j = p_KcolB(i)
            daux = p_DvecP(j)
            Df(k)       = Df(k)       - Dmult(1,3)*p_DB1(i)*daux
            Df(ndofV+k) = Df(ndofV+k) - Dmult(2,3)*p_DB2(i)*daux
            do l = 1, ndofP
              if(j .ne. IdofP(l)) cycle
              
              ! Get the local entries of the B/D matrices
              Da(  k,p+l) = Dmult(1,3)*p_DB1(i)
              Da(o+k,p+l) = Dmult(2,3)*p_DB2(i)
              Da(p+l,  k) = Dmult(3,1)*p_DD1(i)
              Da(p+l,o+k) = Dmult(3,2)*p_DD2(i)
              
              ! Update f_p
              Df(2*ndofV+l) = Df(2*ndofV+l) - Da(p+l,k)*p_DvecU(IdofV(k)) &
                                            - Da(p+l,k)*p_DvecV(IdofV(k))
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

    end if
    
  end subroutine

end module