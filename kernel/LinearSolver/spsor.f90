!##############################################################################
!# ****************************************************************************
!# <name> spsor </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# The following routines can be found in this module:
!#
!# .) spsor_initNavSt2D
!#    -> Initialises a SP-SOR data structure for 2D Navier-Stokes systems.
!#
!# .) spsor_initData
!#    -> Performs data-dependent initialisation of the SP-SOR solver.
!#
!# .) spsor_solve
!#    -> Performs one iteration of the SP-SOR solver.
!#
!# .) spsor_done
!#    -> Releases a SP-SOR data structure.
!#
!# </purpose>
!##############################################################################

module spsor

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

public :: t_spsor
public :: spsor_initNavSt2D
public :: spsor_initData
public :: spsor_done
public :: spsor_solve


!<constants>

!<constantblock>

  ! Use symmetric preconditioner (SP-SSOR)
  integer(I32), parameter, public :: SPSOR_FLAG_SYM  = 1_I32
  
  ! Use diagonal version
  integer(I32), parameter, public :: SPSOR_FLAG_DIAG = 2_I32

!</constantblock>

!<constantblock>
  
  ! SP-SOR solver not initialised
  integer, parameter :: SPSOR_SYSTEM_NONE          = 0
  
  ! SP-SOR solver for 2D Navier-Stokes systems
  integer, parameter :: SPSOR_SYSTEM_NAVST2D       = 1
  
!</constantblock>

!<constantblock>

  ! uninitialised velocity
  integer, parameter :: SPSOR_VELO_NONE  = 0
  
  ! arbitrary velocity
  integer, parameter :: SPSOR_VELO_ARBIT = 1
  
  ! P1~ velocity
  integer, parameter :: SPSOR_VELO_P1T   = 2
  
  ! Q1~ velocity
  integer, parameter :: SPSOR_VELO_Q1T   = 3
  
  ! Q1~ with bubble velocity
  integer, parameter :: SPSOR_VELO_Q1TB  = 4
  
!</constantblock>

!<constantblock>

  ! uninitialised pressure
  integer, parameter :: SPSOR_PRES_NONE  = 0
  
  ! arbitrary pressure
  integer, parameter :: SPSOR_PRES_ARBIT = 1
  
  ! P0/Q0 pressure
  integer, parameter :: SPSOR_PRES_P0    = 2

!</constantblock>

!</constants>


!<types>

!<typeblock>

  ! SP-SOR data structure
  type t_spsor
  
    private
  
    ! System of SP-SOR solver
    integer :: csystem = SPSOR_SYSTEM_NONE

    ! Flags for the SP-SOR solver
    integer(I32) :: cflags = 0_I32
    
    ! Velocity discretisation
    integer :: cvelocity = SPSOR_VELO_NONE
    
    ! Pressure discretisation
    integer :: cpressure = SPSOR_PRES_NONE

    ! Matrix pointers to A11/A22 sub-matrices
    type(t_matrixScalar), pointer :: p_rA11 => null()
    type(t_matrixScalar), pointer :: p_rA22 => null()
    
    ! Matrix pointers to A12/A21 sub-matrices
    type(t_matrixScalar), pointer :: p_rA12 => null()
    type(t_matrixScalar), pointer :: p_rA21 => null()
    
    ! Matrix pointers to B1/B2 sub-matrices
    type(t_matrixScalar), pointer :: p_rB1 => null()
    type(t_matrixScalar), pointer :: p_rB2 => null()
  
    ! Matrix pointers to D1/D2 sub-matrices
    type(t_matrixScalar), pointer :: p_rD1 => null()
    type(t_matrixScalar), pointer :: p_rD2 => null()
    
    ! Matrix pointer to C sub-matrix
    type(t_matrixScalar), pointer :: p_rC => null()
    
    ! Schur-complement matrix
    type(t_matrixScalar) :: rS
    
    ! Tempoary scalar vector.
    ! Is only allocated if cpressure .ne. SPSOR_PRES_Q0.
    type(t_vectorScalar) :: rw
  
  end type

!</typeblock>

!</types>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine spsor_initNavSt2D(rspsor, rmatrix, cflags)

!<description>
  ! Initialises a SP-SOR data structure for 2D Navier-Stokes systems.
!</description>

!<input>
  ! The system matrix of the linear system.
  type(t_matrixBlock), target, intent(IN) :: rmatrix
  
  ! Optional: Flags
  integer(I32), optional, intent(IN) :: cflags
!</input>

!<output>
  ! The SP-SOR data structure that is to be initialised.
  type(t_spsor), intent(OUT) :: rspsor
!</output>

!</subroutine>

  ! Some local variables
  type(t_spatialDiscretisation), pointer :: p_rdiscrV, p_rdiscrP
  integer(I32) :: celemP, celemV
  integer :: ndofPres, ndofVelo, ieldist
  logical :: bdiag

    ! Set the flags
    if(present(cflags)) rspsor%cflags = cflags
    
    ! Set the type
    rspsor%csystem = SPSOR_SYSTEM_NAVST2D

    ! Matrix must be 3x3.
    if ((rmatrix%nblocksPerCol .ne. 3) .or. (rmatrix%nblocksPerRow .ne. 3)) then
      call output_line ('System matrix is not 3x3.',&
          OU_CLASS_ERROR,OU_MODE_STD,'spsor_initNavSt2D')
      call sys_halt()
    end if
    
    ! Hang in the matrix pointers
    rspsor%p_rA11 => rmatrix%RmatrixBlock(1,1)
    rspsor%p_rA22 => rmatrix%RmatrixBlock(2,2)
    rspsor%p_rB1 => rmatrix%RmatrixBlock(1,3)
    rspsor%p_rB2 => rmatrix%RmatrixBlock(2,3)
    rspsor%p_rD1 => rmatrix%RmatrixBlock(3,1)
    rspsor%p_rD2 => rmatrix%RmatrixBlock(3,2)
    
    ! Are the A12/A21 sub-matrices present?
    if(rmatrix%RmatrixBlock(1,2)%NEQ .gt. 0) then
      rspsor%p_rA12 => rmatrix%RmatrixBlock(1,2)
      rspsor%p_rA21 => rmatrix%RmatrixBlock(2,1)
    end if
    
    ! Is the C sub-matrix present?
    if(rmatrix%RmatrixBlock(3,3)%NEQ .gt. 0) &
      rspsor%p_rC => rmatrix%RmatrixBlock(3,3)
    
    ! Todo: Perform some checks
    
    ! Get the total number of dofs
    ndofVelo = rspsor%p_rD1%NCOLS
    ndofPres = rspsor%p_rD1%NEQ
    
    ! Get the discretisations of velocity and pressure spaces
    p_rdiscrV => rspsor%p_rB1%p_rspatialDiscrTest
    p_rdiscrP => rspsor%p_rB1%p_rspatialDiscrTrial
    
    ! Loop over all element distributions
    do ieldist = 1, p_rdiscrP%inumFESpaces
    
      ! Get the element of the spaces
      celemV = elem_getPrimaryElement(p_rdiscrV%RelementDistr(ieldist)%celement)
      celemP = elem_getPrimaryElement(p_rdiscrP%RelementDistr(ieldist)%celement)

      ! What's the velocity element?
      select case(celemV)
      case(EL_P1T_2D)
      
        select case(rspsor%cvelocity)
        case(SPSOR_VELO_NONE, SPSOR_VELO_P1T)
          rspsor%cvelocity = SPSOR_VELO_P1T
        
        case default
          rspsor%cvelocity = SPSOR_VELO_ARBIT
        end select
        
      case(EL_Q1T_2D)
      
        select case(rspsor%cvelocity)
        case(SPSOR_VELO_NONE, SPSOR_VELO_Q1T)
          rspsor%cvelocity = SPSOR_VELO_Q1T
        
        case default
          rspsor%cvelocity = SPSOR_VELO_ARBIT
        end select
        
      case(EL_Q1TB_2D)
      
        select case(rspsor%cvelocity)
        case(SPSOR_VELO_NONE, SPSOR_VELO_Q1TB)
          rspsor%cvelocity = SPSOR_VELO_Q1TB
        
        case default
          rspsor%cvelocity = SPSOR_VELO_ARBIT
        end select
      
      case default
        rspsor%cvelocity = SPSOR_VELO_ARBIT
        
      end select
      
      ! What's the pressure element?
      select case(celemP)
      case(EL_P0_2D, EL_Q0_2D)

        select case(rspsor%cpressure)
        case(SPSOR_PRES_NONE, SPSOR_PRES_P0)
          rspsor%cpressure = SPSOR_PRES_P0
        
        case default
          rspsor%cpressure = SPSOR_PRES_ARBIT
        end select
          
      case default
        rspsor%cpressure = SPSOR_PRES_ARBIT
        
      end select
      
    end do ! ieldist
    
    ! Now we need to set up the structure of the S matrix.
    
    ! Will S be a diagonal matrix?
    bdiag = (rspsor%cpressure .eq. SPSOR_PRES_P0) .or. &
            (iand(rspsor%cflags,SPSOR_FLAG_DIAG) .ne. 0_I32)
    
    if(bdiag) then
    
      ! Create a diagonal matrix
      call lsyssc_createDiagMatrixStruc(rspsor%rS, ndofPres, &
                                        LSYSSC_MATRIX9)

      ! DEBUG !
      call lsyssc_createVector (rspsor%rw, ndofPres, .false.)
    
    else
    
      ! Call the auxiliary routine to assemble the structure of S.
      call spsor_asmMatStructSchur(rspsor%rS, rspsor%p_rD1, rspsor%p_rB1)
      
      ! Allocate a temporary vector
      call lsyssc_createVector (rspsor%rw, ndofPres, .false.)
    
    end if
    
    ! Allocate an empty matrix
    call lsyssc_allocEmptyMatrix (rspsor%rS, LSYSSC_SETM_UNDEFINED)

    ! That's it

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spsor_done(rspsor)

!<description>
  ! Releases a SP-SOR data structure and frees all allocated memory.
!</description>

!<inputoutput>
  ! The SP-SOR data structure that is to be destroyed.
  type(t_spsor), target, intent(INOUT) :: rspsor
!</inputoutput>

!</subroutine>

    ! Release Schur-complement matrix
    if(lsyssc_hasMatrixStructure(rspsor%rS)) &
      call lsyssc_releaseMatrix(rspsor%rS)
    
    ! Release temporary vector
    if(rspsor%rw%NEQ .gt. 0) &
      call lsyssc_releaseVector(rspsor%rw)
      
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spsor_asmMatStructSchur(rS, rD, rB)

!<description>
  ! PRIVATE AUXILIARY ROUTINE:
  ! Assembles the matrix structure of the Schur-Complement matrix.
!</description>

!<input>
  ! The matrix D
  type(t_matrixScalar), intent(IN) :: rD

  ! The matrix B
  type(t_matrixScalar), intent(IN) :: rB
!</input>

!<output>
  ! The matrix S
  type(t_matrixScalar), intent(OUT) :: rS
!</output>

!</subroutine>

  ! access for the matrix arrays
  integer, dimension(:), pointer :: p_KldB,p_KldD,p_KldS,&
                                    p_KcolB,p_KcolD,p_KcolS,p_KdiagS
  
  ! auxiliary arrays
  integer, dimension(:,:), allocatable :: Iadj, Idofs
  integer, dimension(:), allocatable :: Imask
  
  ! some scalars
  integer :: idofp,idofv,ndofp,ndofv,ndegD,ndegB
  integer :: i,j,k,id1,ib1,ind,idx

    ! Fetch the arrays of the input matrices
    call lsyssc_getbase_Kld(rB, p_KldB)
    call lsyssc_getbase_Kcol(rB, p_KcolB)
    call lsyssc_getbase_Kld(rD, p_KldD)
    call lsyssc_getbase_Kcol(rD, p_KcolD)
    
    ndofp = rD%NEQ
    ndofv = rD%NCOLS
    
    ! Prepare the matrix S and allocate the row-pointer array
    rS%cmatrixFormat = LSYSSC_MATRIX9
    rS%NEQ = ndofp
    rS%NCOLS = ndofp
    call storage_new ('spsor_asmMatStructSchur', 'Kld', ndofp+1, ST_INT,&
                      rS%h_Kld, ST_NEWBLOCK_NOINIT)
    call storage_getbase_int(rS%h_Kld, p_KldS)
    
    ! First of all, calculate the degree of B and D
    ndegD = 0
    do i = 1, ndofp
      ndegD = max(ndegD, p_KldD(i+1) - p_KldD(i))
    end do
    ndegB = 0
    do i = 1, ndofv
      ndegB = max(ndegB, p_KldB(i+1) - p_KldB(i))
    end do
    
    ! allocate auxiliary arrays
    allocate(Iadj(ndegB,ndegD))
    allocate(Imask(ndegB))
    allocate(Idofs(ndegB,ndofp))
    
    ! Calculate the number of non-zero entries of S
    p_KldS(1) = 1
    rS%NA = 0
    do idofp = 1, ndofp
    
      ! Format arrays
      Iadj(:,:) = 0
      Imask(:) = 1
      Idofs(:,idofp) = 0
    
      ! Loop over all non-zeroes of D
      id1 = p_KldD(idofp)
      do i = id1, p_KldD(idofp+1)-1
      
        ! Get the velocity DOF index
        idofv = p_KcolD(i)
        
        ! Now loop over the row idofv of B
        ib1 = p_KldB(idofv)
        do j = ib1, p_KldB(idofv+1)-1
        
          ! Store the index of the pressure dof
          Iadj(j-ib1+1, i-id1+1) = p_KcolB(j)
          
        end do ! j
      
      end do ! i
    
      ! Calculate the number of entries in the current row of D
      ind = p_KldD(idofp+1) - id1
      
      ! Now Iadj contains the indices of all pressure DOFs to which the
      ! current pressure DOF idofp is adjacent via one of the velocity DOFs.
      ! We now need to find out which pressure DOFs are on the same element
      ! of the mesh.
      ! If a pressure DOF i is on the same element as our current pressure
      ! DOF idofp, then i exists in each row of Iadj.
      do i = 2, ind
      
        ! Loop over the row i of Iadj
        do j = 1, ndegB
          
          ! Get the index of the adjacent DOF
          idx = Iadj(j,i)
          if(idx .eq. 0) exit
          
          ! Try to find the DOF in the first row of Iadj
          do k = 1, ndegB
            if(Iadj(k,1) .eq. idx) then
              Imask(k) = Imask(k) + 1
              exit
            end if
          end do
          
        end do ! j
      
      end do ! j
      
      ! Now each adjacent DOF in the first row of Iadj whose entry in Imask
      ! is equal to ind is on the same element as idofp.
      k = 0
      do i = 1, ind
        if(Imask(i) .eq. ind) then
          k = k+1
          Idofs(k,idofp) = Iadj(i,1)
        end if
      end do
      
      ! Set up row pointer for next row
      rS%NA = rS%NA + k
      p_KldS(idofp+1) = p_KldS(idofp) + k
    
    end do ! idofp
    
    ! Now we know the number of non-zeroes in S, so let's allocate column index
    ! array.
    call storage_new ('spsor_asmMatStructSchur', 'Kcol', rS%NA, ST_INT,&
                      rS%h_Kcol, ST_NEWBLOCK_NOINIT)
    call storage_getbase_int(rS%h_Kcol, p_KcolS)
    
    ! Now loop over all pressure DOFs
    do idofp = 1, ndofp
    
      ! Get the number of non-zeroes
      j = p_KldS(idofp)
      k = p_KldS(idofp+1) - j
    
      ! Set up the column indices
      do i = 1, k
        p_KcolS(j+i-1) = Idofs(i,idofp)
      end do
      
    end do
    
    ! Release auxiliary arrays
    deallocate(Idofs)
    deallocate(Imask)
    deallocate(Iadj)
    
    ! Finally, set up the diagonal pointer array
    call storage_new ('spsor_asmMatStructSchur', 'Kdiagonal', ndofp, ST_INT,&
                      rS%h_Kdiagonal, ST_NEWBLOCK_NOINIT)
    call storage_getbase_int(rS%h_Kdiagonal, p_KdiagS)
    
    ! Build the diagonal pointer array
    do i = 1, ndofp
      do j = p_KldS(i), p_KldS(i+1)-1
        if(p_KcolS(j) .ge. i) then
          p_KdiagS(i) = j
          exit
        end if
      end do ! j
    end do
    
    ! That's it

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine spsor_initData(rdata)

!<description>
  ! This subroutine performs data-dependent initialisation of the SP-SOR
  ! solver. 
  ! The data structure must have already been set using one of the
  ! spsor_initXXXX routines (e.g. spsor_initNavSt2D).
!</description>

!<inputoutput>
  ! The SP-SOR data structure
  type(t_spsor), target, intent(INOUT) :: rdata
!</inputoutput>

!</subroutine>

    ! Okay, which SP-SOR solver do we have here?
    select case(rdata%csystem)
    case (SPSOR_SYSTEM_NAVST2D)
      
      ! SP-SOR for 2D Navier-Stokes systems
      
      ! Do we use the diagonal or full version?
      if(iand(rdata%cflags,SPSOR_FLAG_DIAG) .ne. 0_I32) then

        ! Call the diagonal version
        call spsor_initNavSt2D_diag(rdata)
      
      else
      
        ! Do we use ony P0/Q0 for pressure?
        ! We have some fast implementations for this case...
        select case(rdata%cpressure)
        case (SPSOR_PRES_P0)
        
          ! Do we even have one of our special velocity elements?
          ! This would speed up the whole process even more...
          select case(rdata%cvelocity)
!          case (SPSOR_VELO_P1T)
!            todo
!          case (SPSOR_VELO_Q1T)
!            todo
!          case (SPSOR_VELO_Q1TB)
!            todo
        
          case default
            ! P0/Q0 pressure, arbitrary velocity
            call spsor_initNavSt2D_P0_full(rdata)
          
          end select
      
        case default
          ! arbitrary pressure, arbitrary velocity
          call spsor_initNavSt2D_full(rdata)
          
        end select

      end if

    case default
      call output_line ('SP-SOR data structure not initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'spsor_initData')
      call sys_halt()

    end select
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spsor_initNavSt2D_diag(rdata)

!<description>
  ! Internal subroutine:
  ! Peforms data-dependent initialisation of the diagonal SP-SOR algorithm for
  ! 2D Navier-Stokes systems: arbitrary pressure, arbitrary velocity
!</description>

!<inputoutput>
  ! The SP-SOR data structure
  type(t_spsor), target, intent(INOUT) :: rdata
!</inputoutput>

!</subroutine>

  ! Scaling factors
  real(DP) :: dsf1, dsf2, dsfC
  
  ! Quick access for the matrix arrays
  integer, dimension(:), pointer :: p_KldA,p_KldB,p_KldC,p_KldD,&
      p_KcolA,p_KcolB,p_KcolC,p_KcolD,p_KdiagA,p_KdiagC
  real(DP), dimension(:), pointer :: p_DA11,p_DA22,p_DB1,p_DB2,&
      p_DD1,p_DD2,p_DC,p_DS
  
  ! local variables
  logical :: bHaveC
  real(DP) :: dc,daux1,daux2
  integer :: idofp,idofu,i,j,k,ndofV

    ! Let's assume we do not have the optional matrices
    bHaveC = .FALSE.
    
    ! Fetch the sub-matrix arrays
    call lsyssc_getbase_Kld(rdata%p_rA11, p_KldA)
    call lsyssc_getbase_Kcol(rdata%p_rA11, p_KcolA)
    call lsyssc_getbase_Kdiagonal(rdata%p_rA11, p_KdiagA)
    call lsyssc_getbase_double(rdata%p_rA11, p_DA11)
    call lsyssc_getbase_double(rdata%p_rA22, p_DA22)

    call lsyssc_getbase_Kld(rdata%p_rB1, p_KldB)
    call lsyssc_getbase_Kcol(rdata%p_rB1, p_KcolB)
    call lsyssc_getbase_double(rdata%p_rB1, p_DB1)
    call lsyssc_getbase_double(rdata%p_rB2, p_DB2)

    call lsyssc_getbase_Kld(rdata%p_rD1, p_KldD)
    call lsyssc_getbase_Kcol(rdata%p_rD1, p_KcolD)
    call lsyssc_getbase_double(rdata%p_rD1, p_DD1)
    call lsyssc_getbase_double(rdata%p_rD2, p_DD2)
    
    if(associated(rdata%p_rC)) then
      bHaveC = .true.
      call lsyssc_getbase_Kld(rdata%p_rC, p_KldC)
      call lsyssc_getbase_Kcol(rdata%p_rC, p_KcolC)
      call lsyssc_getbase_Kdiagonal(rdata%p_rC, p_KdiagC)
      call lsyssc_getbase_double(rdata%p_rC, p_DC)
    end if

    call lsyssc_getbase_double(rdata%rS, p_DS)
    
    ! Calculate the scaling factors
    dsf1 = (rdata%p_rB1%dscaleFactor * rdata%p_rD1%dscaleFactor) & 
         / rdata%p_rA11%dscaleFactor
    dsf2 = (rdata%p_rB2%dscaleFactor * rdata%p_rD2%dscaleFactor) & 
         / rdata%p_rA22%dscaleFactor
    
    if(bHaveC) then
      dsfC = rdata%p_rC%dscaleFactor
      if(dsfC .eq. 0.0_DP) bHaveC = .false.
    end if

    ! Let's loop over all rows of S
    do idofp = 1, rdata%rS%NEQ
    
      ! Format the local matrices
      dc = 0.0_DP
      
      ! Format the Schur-complement entry for this pressure DOF
      p_DS(idofp) = 0.0_DP

      ! Fetch the number of velocity DOFs adjacent to this pressure DOF
      ndofV = p_KldD(idofp+1) - p_KldD(idofp)
      
      ! If the C matrix exists, grab it's main diagonal entry
      if(bHaveC) then
        dc = dsfC*p_DC(p_KdiagC(idofp))
      end if
      
      ! Let's loop over all velocity DOFs which are adjacent to the current
      ! pressure DOF.
      daux1 = 0.0_DP
      daux2 = 0.0_DP
      do j = p_KldD(idofp), p_KldD(idofp+1)-1
      
        ! Get the index of the velocity DOF
        idofu = p_KcolD(j)
        
        ! Fetch the main diagonal index of A
        k = p_KdiagA(idofu)

        ! Try to find the corrent entry in B1/B2
        do i = p_KldB(idofu), p_KldB(idofu+1)
          if(p_KcolB(i) .eq. idofp) then
            daux1 = daux1 + (p_DB1(i)*p_DD1(j)) / p_DA11(k)
            daux2 = daux2 + (p_DB2(i)*p_DD2(j)) / p_DA22(k)
            exit
          end if
        end do ! i
      
      end do ! id1
      
      ! Okay, let's calculate the Schur-complement dc = C - D * A^-1 * B
      dc = dc - dsf1*daux1 - dsf2*daux2
      
      ! Now if the Schur-complement matrix is regular, we'll store the inverse
      ! in the corresponding entry in p_DS.
      if(abs(dc) .gt. SYS_EPSREAL) p_DS(idofp) = 1.0_DP / dc
      
    end do ! idofp
    
    ! That's it
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spsor_initNavSt2D_P0_full(rdata)

!<description>
  ! Internal subroutine:
  ! Peforms data-dependent initialisation of the SP-SOR algorithm for
  ! 2D Navier-Stokes systems: P0/Q0 pressure, arbitrary velocity
!</description>

!<inputoutput>
  ! The SP-SOR data structure
  type(t_spsor), target, intent(INOUT) :: rdata
!</inputoutput>

!</subroutine>

  ! Multiplication factors
  real(DP) :: dsfA11, dsfA12, dsfA21, dsfA22, dsfB1, dsfB2, dsfD1, dsfD2, dsfC
  
  ! Quick access for the matrix arrays
  integer, dimension(:), pointer :: p_KldA,p_KldA12,p_KldB,p_KldC,p_KldD,&
      p_KcolA,p_KcolA12,p_KcolB,p_KcolC,p_KcolD,p_KdiagC
  real(DP), dimension(:), pointer :: p_DA11,p_DA12,p_DA21,p_DA22,p_DB1,p_DB2,&
      p_DD1,p_DD2,p_DC,p_DS

  ! Temporary local matrices
  real(DP), dimension(:,:), allocatable :: DA
  real(DP), dimension(:), allocatable :: DB
  
  ! pivot array for LAPACK
  integer, dimension(:), allocatable :: Ipivot
  
  ! local variables
  logical :: bHaveA12, bHaveC
  real(DP) :: dc,daux1,daux2
  integer :: idofp,idofu,i,j1,j2,k,id1,id2,nmaxdofV,ndofV,info

    ! Let's assume we do not have the optional matrices
    bHaveA12 = .FALSE.
    bHaveC = .FALSE.
    
    ! Fetch the sub-matrix arrays
    call lsyssc_getbase_Kld(rdata%p_rA11, p_KldA)
    call lsyssc_getbase_Kcol(rdata%p_rA11, p_KcolA)
    call lsyssc_getbase_double(rdata%p_rA11, p_DA11)
    call lsyssc_getbase_double(rdata%p_rA22, p_DA22)

    call lsyssc_getbase_Kld(rdata%p_rB1, p_KldB)
    call lsyssc_getbase_Kcol(rdata%p_rB1, p_KcolB)
    call lsyssc_getbase_double(rdata%p_rB1, p_DB1)
    call lsyssc_getbase_double(rdata%p_rB2, p_DB2)

    call lsyssc_getbase_Kld(rdata%p_rD1, p_KldD)
    call lsyssc_getbase_Kcol(rdata%p_rD1, p_KcolD)
    call lsyssc_getbase_double(rdata%p_rD1, p_DD1)
    call lsyssc_getbase_double(rdata%p_rD2, p_DD2)
    
    if(associated(rdata%p_rC)) then
      bHaveC = .true.
      call lsyssc_getbase_Kld(rdata%p_rC, p_KldC)
      call lsyssc_getbase_Kcol(rdata%p_rC, p_KcolC)
      call lsyssc_getbase_Kdiagonal(rdata%p_rC, p_KdiagC)
      call lsyssc_getbase_double(rdata%p_rC, p_DC)
    end if
    if(associated(rdata%p_rA12)) then
      bHaveA12 = .true.
      call lsyssc_getbase_Kld(rdata%p_rA12, p_KldA12)
      call lsyssc_getbase_Kcol(rdata%p_rA12, p_KcolA12)
      call lsyssc_getbase_double(rdata%p_rA12, p_DA12)
      call lsyssc_getbase_double(rdata%p_rA21, p_DA21)
    end if

    call lsyssc_getbase_double(rdata%rS, p_DS)
    
    ! Fetch the multiplication factors
    dsfA11 = rdata%p_rA11%dscaleFactor
    dsfA22 = rdata%p_rA22%dscaleFactor
    dsfB1 = rdata%p_rB1%dscaleFactor
    dsfB2 = rdata%p_rB2%dscaleFactor
    dsfD1 = rdata%p_rD1%dscaleFactor
    dsfD2 = rdata%p_rD2%dscaleFactor
    
    if(bHaveC) then
      dsfC = rdata%p_rC%dscaleFactor
      if(dsfC .eq. 0.0_DP) bHaveC = .false.
    end if
    
    if(bHaveA12) then
      dsfA12 = rdata%p_rA12%dscaleFactor
      dsfA21 = rdata%p_rA21%dscaleFactor
      if((dsfA12 .eq. 0.0_DP) .and. (dsfA21 .eq. 0.0_DP)) &
        bHaveA12 = .false.
    end if
    
    ! Determine the maximum number of velocity DOFs that one
    ! pressure DOF may be adjacent to.
    nmaxdofV = 0
    do idofp = 1, rdata%rS%NEQ
      nmaxdofV = max(nmaxdofV, p_KldD(idofp+1)-p_KldD(idofp))
    end do ! idofp
    nmaxdofV = 2*nmaxdofV
    
    ! Okay, let's allocate the temporary data
    allocate(DA(nmaxdofV,nmaxdofV))
    allocate(DB(nmaxdofV))
    allocate(Ipivot(nmaxdofV))
    
    ! Let's loop over all rows of S
    do idofp = 1, rdata%rS%NEQ
    
      ! Format the local matrices
      DA = 0.0_DP
      DB = 0.0_DP
      dc = 0.0_DP
      
      ! Format the Schur-complement entry for this pressure DOF
      p_DS(idofp) = 0.0_DP

      ! Fetch the number of velocity DOFs adjacent to this pressure DOF
      ndofV = p_KldD(idofp+1) - p_KldD(idofp)
      
      ! If the C matrix exists, grab it's main diagonal entry
      if(bHaveC) then
        dc = dsfC*p_DC(p_KdiagC(idofp))
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
              DA(      j1,      j2) = dsfA11*p_DA11(i)
              DA(ndofV+j1,ndofV+j2) = dsfA22*p_DA22(i)
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
                DA(      j1,ndofV+j2) = dsfA12*p_DA12(i)
                DA(ndofV+j1,      j2) = dsfA21*p_DA21(i)
                exit
              end if
            end do ! id2
          end do ! i
        end if

        ! Let's fetch the local B matrices
        do i = p_KldB(idofu), p_KldB(idofu+1)
          if(p_KcolB(i) .eq. idofP) then
            j1 = id1 - p_KldD(idofp) + 1
            DB(      j1) = dsfB1*p_DB1(i)
            DB(ndofV+j1) = dsfB2*p_DB2(i)
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
      dc = dc - dsfD1*daux1 - dsfD2*daux2
      
      ! Now if the Schur-complement matrix is regular, we'll store the inverse
      ! in the corresponding entry in p_DS.
      if(abs(dc) .gt. SYS_EPSREAL) p_DS(idofp) = 1.0_DP / dc
      
    end do ! idofp
    
    ! And release the temporary memory
    deallocate(Ipivot)
    deallocate(DB)
    deallocate(DA)
    
    ! That's it
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spsor_initNavSt2D_full(rdata)

!<description>
  ! Internal subroutine:
  ! Peforms data-dependent initialisation of the SP-SOR algorithm for
  ! 2D Navier-Stokes systems: discontinous pressure, arbitrary velocity
!</description>

!<inputoutput>
  ! The SP-SOR data structure
  type(t_spsor), target, intent(INOUT) :: rdata
!</inputoutput>

!</subroutine>

  ! Multiplication factors
  real(DP) :: dsfA11, dsfA12, dsfA21, dsfA22, dsfB1, dsfB2, dsfD1, dsfD2, dsfC
  
  ! Quick access for the matrix arrays
  integer, dimension(:), pointer :: p_KldA,p_KldA12,p_KldB,p_KldC,p_KldD,p_KldS,&
      p_KcolA,p_KcolA12,p_KcolB,p_KcolC,p_KcolD,p_KcolS,p_KdiagC
  real(DP), dimension(:), pointer :: p_DA11,p_DA12,p_DA21,p_DA22,p_DB1,p_DB2,&
      p_DD1,p_DD2,p_DC,p_DS

  ! Temporary local matrices
  real(DP), dimension(:,:), allocatable :: DA,DB,DC,DS
  
  integer, dimension(:), allocatable :: IdofsV, IdofsP
  
  ! pivot array for LAPACK
  integer, dimension(:), allocatable :: Ipivot
  
  ! local variables
  logical :: bHaveA12, bHaveC
  real(DP) :: daux1,daux2
  integer :: irow,idofp,idofu,i,j,k,idx,nmaxdofV,nmaxdofP,ndofV,ndofP,info
  integer :: id1,id2,is1
  
    ! Let's assume we do not have the optional matrices
    bHaveA12 = .FALSE.
    bHaveC = .FALSE.
    
    ! Fetch the sub-matrix arrays
    call lsyssc_getbase_Kld(rdata%p_rA11, p_KldA)
    call lsyssc_getbase_Kcol(rdata%p_rA11, p_KcolA)
    call lsyssc_getbase_double(rdata%p_rA11, p_DA11)
    call lsyssc_getbase_double(rdata%p_rA22, p_DA22)

    call lsyssc_getbase_Kld(rdata%p_rB1, p_KldB)
    call lsyssc_getbase_Kcol(rdata%p_rB1, p_KcolB)
    call lsyssc_getbase_double(rdata%p_rB1, p_DB1)
    call lsyssc_getbase_double(rdata%p_rB2, p_DB2)

    call lsyssc_getbase_Kld(rdata%p_rD1, p_KldD)
    call lsyssc_getbase_Kcol(rdata%p_rD1, p_KcolD)
    call lsyssc_getbase_double(rdata%p_rD1, p_DD1)
    call lsyssc_getbase_double(rdata%p_rD2, p_DD2)
    
    if(associated(rdata%p_rC)) then
      bHaveC = .true.
      call lsyssc_getbase_Kld(rdata%p_rC, p_KldC)
      call lsyssc_getbase_Kcol(rdata%p_rC, p_KcolC)
      call lsyssc_getbase_Kdiagonal(rdata%p_rC, p_KdiagC)
      call lsyssc_getbase_double(rdata%p_rC, p_DC)
    end if
    if(associated(rdata%p_rA12)) then
      bHaveA12 = .true.
      call lsyssc_getbase_Kld(rdata%p_rA12, p_KldA12)
      call lsyssc_getbase_Kcol(rdata%p_rA12, p_KcolA12)
      call lsyssc_getbase_double(rdata%p_rA12, p_DA12)
      call lsyssc_getbase_double(rdata%p_rA21, p_DA21)
    end if

    call lsyssc_getbase_Kld(rdata%rS, p_KldS)
    call lsyssc_getbase_Kcol(rdata%rS, p_KcolS)
    call lsyssc_getbase_double(rdata%rS, p_DS)
    
    ! Fetch the multiplication factors
    dsfA11 = rdata%p_rA11%dscaleFactor
    dsfA22 = rdata%p_rA22%dscaleFactor
    dsfB1 = rdata%p_rB1%dscaleFactor
    dsfB2 = rdata%p_rB2%dscaleFactor
    dsfD1 = rdata%p_rD1%dscaleFactor
    dsfD2 = rdata%p_rD2%dscaleFactor
    
    if(bHaveC) then
      dsfC = rdata%p_rC%dscaleFactor
      if(dsfC .eq. 0.0_DP) bHaveC = .false.
    end if
    
    if(bHaveA12) then
      dsfA12 = rdata%p_rA12%dscaleFactor
      dsfA21 = rdata%p_rA21%dscaleFactor
      if((dsfA12 .eq. 0.0_DP) .and. (dsfA21 .eq. 0.0_DP)) &
        bHaveA12 = .false.
    end if
    
    ! Determine the maximum number of velocity DOFs that one
    ! pressure DOF may be adjacent to.
    nmaxdofV = 0
    nmaxdofP = 0
    do idofp = 1, rdata%rS%NEQ
      nmaxdofV = max(nmaxdofV, p_KldD(idofp+1)-p_KldD(idofp))
      nmaxdofP = max(nmaxdofP, p_KldS(idofp+1)-p_KldS(idofp))
    end do ! idofp
    nmaxdofV = 2*nmaxdofV
    
    ! Okay, let's allocate the auxiliary arrays
    allocate(DA(nmaxdofV,nmaxdofV))
    allocate(DB(nmaxdofV,nmaxdofP))
    allocate(DC(nmaxdofP,nmaxdofP))
    allocate(DS(nmaxdofP,nmaxdofP))
    allocate(Ipivot(max(nmaxdofV,nmaxdofP)))
    allocate(IdofsV(nmaxdofV))
    allocate(IdofsP(nmaxdofP))
    
    ! Let's loop over all rows of S
    do irow = 1, rdata%rS%NEQ
    
      ! Check the first non-zero entry of S. If its column index is less
      ! than irow, then we can to skip it.
      if(p_KcolS(p_KldS(irow)) .lt. irow) cycle
    
      ! Format the local matrices
      DA = 0.0_DP
      DB = 0.0_DP
      DC = 0.0_DP
      
      ! Fetch the number of velocity and pressure DOFs
      ndofV = p_KldD(irow+1) - p_KldD(irow)
      ndofP = p_KldS(irow+1) - p_KldS(irow)
      
      ! Fetch the indices of the velocity DOFs.
      k = p_KldD(irow)
      do i = k, p_KldD(irow+1)-1
        IdofsV(i-k+1) = p_KcolD(i)
      end do ! i

      ! Fetch the indices of the pressure DOFs
      k = p_KldS(irow)
      do i = k, p_KldS(irow+1)-1
        IdofsP(i-k+1) = p_KcolS(i)
      end do ! i
      
      ! If the C matrix exists, grab it's main diagonal entry
      if(bHaveC) then
      
        ! Loop over all local pressure DOFs
        do i = 1, ndofP

          ! Get the index of the pressure DOF
          idofp = IdofsP(i)
          
          ! Loop over the row of C
          do k = p_KldC(idofp), p_KldC(idofp+1)-1
          
            ! Get the index of the pressure DOF
            idx = p_KcolC(k)
            
            ! Let's see if this corresponds to one of our local pressure DOFs
            do j = 1, ndofP
              if(idx .eq. IdofsP(j)) then
                DC(i,j) = dsfC*p_DC(k)
                exit
              end if
            end do ! j
          end do ! k
        end do ! i

      end if
      
      ! Let's loop over all velocity DOFs which are adjacent to the current
      ! pressure DOF.
      do i = 1, ndofV
      
        ! Get the index of the velocity DOF
        idofu = IdofsV(i)
        
        ! Let's fetch the local A11/A22 matrices
        do k = p_KldA(idofu), p_KldA(idofu+1)-1
        
          ! Get the column index
          idx = p_KcolA(k)
          
          ! Let's see if this corresponds to one of our local velocity DOFs
          do j = 1, ndofV
            if(idx .eq. IdofsV(j)) then
              DA(      i,      j) = dsfA11*p_DA11(k)
              DA(ndofV+i,ndofV+j) = dsfA22*p_DA22(k)
              exit
            end if
          end do ! j
        end do ! k

        ! Do the A12/A21 matrices exist? If yes, then we also need to grab
        ! their local sub-matrices.
        if(bHaveA12) then
          ! Let's fetch the local A12/A21 matrices
          do k = p_KldA12(idofu), p_KldA12(idofu+1)-1
          
            ! Get the column index
            idx = p_KcolA12(k)
            
            ! Let's see if this corresponds to one of our local velocity DOFs
            do j = 1, ndofV
              if(idx .eq. IdofsV(j)) then
                DA(      i,ndofV+j) = dsfA12*p_DA12(k)
                DA(ndofV+i,      j) = dsfA21*p_DA21(k)
                exit
              end if
            end do ! j
          end do ! k
        end if

        ! Let's fetch the local B matrices
        do k = p_KldB(idofu), p_KldB(idofu+1)
        
          ! Get the column index
          idx = p_KcolB(k)
          
          ! Let's see if this corresponds to one of our local pressure DOFs
          do j = 1, ndofP
            if(idx .eq. IdofsP(j)) then
              DB(      i,j) = dsfB1*p_DB1(k)
              DB(ndofV+i,j) = dsfB2*p_DB2(k)
              exit
            end if
          end do ! j
        end do ! k
      
      end do ! i
      
      ! Okay, try to solve the local system A*X = B
      call DGESV(2*ndofV,1,DA,nmaxdofV,Ipivot,DB,nmaxdofV,info)
      
      ! Did LAPACK fail? If yes, simply continue with the next pressure DOF.
      if(info .ne. 0) cycle
      
      ! Okay, let's calculate the Schur-complement C := C - D * A^-1 * B
      ! Loop over all local pressure DOFs
      do i = 1, ndofP
      
        ! Get the index of the pressure DOF
        idofp = IdofsP(i)

        id1 = p_KldD(idofp)
        id2 = p_KldD(idofp+1)-1
        
        ! Another loop over all local pressure DOFs
        do j = 1, ndofP
        
          daux1 = 0.0_DP
          daux2 = 0.0_DP
          do k = id1, id2
            daux1 = daux1 + p_DD1(k)*DB(      k-id1+1,j)
            daux2 = daux2 + p_DD2(k)*DB(ndofP+k-id1+1,j)
          end do ! k
          DC(i,j) = DC(i,j) - dsfD1*daux1 - dsfD2*daux2
        
        end do ! j
      
      end do ! i

      ! Okay, try to solve the local system C*S = id
      DS = 0.0_DP
      do i = 1, ndofP
        DS(i,i) = 1.0_DP
      end do
      call DGESV(ndofP,1,DC,nmaxdofP,Ipivot,DS,nmaxdofP,info)
      
      ! Now if the Schur-complement matrix is regular, we'll store the inverse
      ! in the corresponding entry in p_DS.
      if(info .ne. 0) cycle
      
      ! Loop over all local pressure DOFs
      do i = 1, ndofP

        ! Get the index of the pressure DOF
        idofp = IdofsP(i)
        
        ! Loop over the row of S
        is1 = p_KldS(idofp)
        do k = is1, p_KldS(idofp+1)-1
          p_DS(k) = DS(i, k-is1+1)
        end do ! k
      end do ! i
      
    end do ! idofp
    
    ! And release the temporary memory
    deallocate(Ipivot)
    deallocate(DB)
    deallocate(DA)
    
    ! That's it
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spsor_solve(rdata, rx, rf, drelax)

!<description>
  ! This routine performs one SP-SOR iteration.
!</description>

!<input>
  ! The SP-SOR data structure
  type(t_spsor), target, intent(INOUT) :: rdata
  
  ! The right-hand-side vector
  type(t_vectorBlock), intent(IN) :: rf
  
  ! The relaxation parameter. Default = 1.0
  real(DP), intent(IN) :: drelax
!</input>

!<inputoutput>
  ! The iteration vector that is to be updated.
  type(t_vectorblock), intent(INOUT) :: rx
!</inputoutput>

!</subroutine>

    ! Make sure we do not apply the symmetric variant ('SP-SSOR')
    if(iand(rdata%cflags, SPSOR_FLAG_SYM) .ne. 0_I32) then
      call output_line ('SP-SSOR cannot be called via spsor_solve',&
          OU_CLASS_ERROR,OU_MODE_STD,'spsor_solve')
      call sys_halt()
    end if

    ! Okay, which SP-SOR solver do we have here?
    select case(rdata%csystem)
    case (SPSOR_SYSTEM_NAVST2D)
    
      ! Call the routine for 2D Navier-Stokes systems
      call spsor_solve_NavSt2D(rdata, rx, rf, drelax)
      
    case default
      call output_line ('SP-SOR data structure not initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'spsor_solve')
      call sys_halt()
    
    end select
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spsor_solve_NavSt2D(rdata, rx, rf, drelax)

!<description>
  ! INTERNAL ROUTINE:
  ! Perform one SP-SOR iteration for 2D Navier-Stokes systems.
!</description>

!<input>
  ! The SP-SOR data structure
  type(t_spsor), target, intent(INOUT) :: rdata
  
  ! The right-hand-side vector
  type(t_vectorBlock), intent(IN) :: rf
  
  ! The relaxation parameter. Default = 1.0
  real(DP), intent(IN) :: drelax
!</input>

!<inputoutput>
  ! The iteration vector that is to be updated.
  type(t_vectorblock), intent(INOUT) :: rx
!</inputoutput>

!</subroutine>

  ! Multiplication factors
  real(DP) :: dsfA11, dsfA12, dsfA21, dsfA22, dsfB1, dsfB2, dsfD1, dsfD2
  
  ! access for the matrix arrays
  integer, dimension(:), pointer :: p_KldA,p_KldA12,p_KldB,p_KldD,p_KldS,&
      p_KcolA,p_KcolA12,p_KcolB,p_KcolD,p_KcolS,p_KdiagA,p_KdiagS
  real(DP), dimension(:), pointer :: p_DA11,p_DA12,p_DA21,p_DA22,p_DB1,p_DB2,&
      p_DD1,p_DD2,p_DS
  
  ! vector arrays
  real(DP), dimension(:), pointer :: p_Du1, p_Du2, p_Dup, p_Df1, p_Df2, p_Dfp,&
      p_Dw

  integer :: m,n
  logical :: bHaveA12, bdiag
  
    ! Let's assume we do not have the optional matrices
    bHaveA12 = .FALSE.
    
    ! Fetch the sub-matrix arrays
    call lsyssc_getbase_Kld(rdata%p_rA11, p_KldA)
    call lsyssc_getbase_Kcol(rdata%p_rA11, p_KcolA)
    call lsyssc_getbase_Kdiagonal(rdata%p_rA11, p_KdiagA)
    call lsyssc_getbase_double(rdata%p_rA11, p_DA11)
    call lsyssc_getbase_double(rdata%p_rA22, p_DA22)
    
    if(associated(rdata%p_rA12)) then
      bHaveA12 = .true.
      call lsyssc_getbase_Kld(rdata%p_rA12, p_KldA12)
      call lsyssc_getbase_Kcol(rdata%p_rA12, p_KcolA12)
      call lsyssc_getbase_double(rdata%p_rA12, p_DA12)
      call lsyssc_getbase_double(rdata%p_rA21, p_DA21)
    end if

    call lsyssc_getbase_Kld(rdata%p_rB1, p_KldB)
    call lsyssc_getbase_Kcol(rdata%p_rB1, p_KcolB)
    call lsyssc_getbase_double(rdata%p_rB1, p_DB1)
    call lsyssc_getbase_double(rdata%p_rB2, p_DB2)

    call lsyssc_getbase_Kld(rdata%p_rD1, p_KldD)
    call lsyssc_getbase_Kcol(rdata%p_rD1, p_KcolD)
    call lsyssc_getbase_double(rdata%p_rD1, p_DD1)
    call lsyssc_getbase_double(rdata%p_rD2, p_DD2)

    call lsyssc_getbase_Kld(rdata%rS, p_KldS)
    call lsyssc_getbase_Kcol(rdata%rS, p_KcolS)
    call lsyssc_getbase_Kdiagonal(rdata%rS, p_KdiagS)
    call lsyssc_getbase_double(rdata%rS, p_DS)
    
    ! Fetch the multiplication factors
    dsfA11 = rdata%p_rA11%dscaleFactor
    dsfA22 = rdata%p_rA22%dscaleFactor
    dsfB1 = rdata%p_rB1%dscaleFactor
    dsfB2 = rdata%p_rB2%dscaleFactor
    dsfD1 = rdata%p_rD1%dscaleFactor
    dsfD2 = rdata%p_rD2%dscaleFactor
    
    if(bHaveA12) then
      dsfA12 = rdata%p_rA12%dscaleFactor
      dsfA21 = rdata%p_rA21%dscaleFactor
      if((dsfA12 .eq. 0.0_DP) .and. (dsfA21 .eq. 0.0_DP)) &
        bHaveA12 = .false.
    end if
    
    ! Fetch the vector arrays
    call lsyssc_getbase_double(rx%RvectorBlock(1), p_Du1)
    call lsyssc_getbase_double(rx%RvectorBlock(2), p_Du2)
    call lsyssc_getbase_double(rx%RvectorBlock(3), p_Dup)
    call lsyssc_getbase_double(rf%RvectorBlock(1), p_Df1)
    call lsyssc_getbase_double(rf%RvectorBlock(2), p_Df2)
    call lsyssc_getbase_double(rf%RvectorBlock(3), p_Dfp)
    if(rdata%cpressure .ne. SPSOR_PRES_P0) then
      call lsyssc_getbase_double(rdata%rw, p_Dw)
    end if
    
    ! Fetch the number of velocity and pressure DOFs
    m = rdata%p_rA11%NEQ
    n = rdata%rS%NEQ
    
    
    ! Step One
    ! --------
    ! Solve the following equation using SOR:
    !
    !                        A * u = f - B*p
    
    ! Do we have the A12/A21 matrices?
    if(bHaveA12) then
    
      ! Solve A11 * u1 = f1 - A12*u2 - B1*p
      call spsor_aux_solveAAB(m,p_Du1,p_Du2,p_Dup,p_Df1,&
          p_KldA,p_KcolA,p_KdiagA,p_DA11,dsfA11,&
          p_KldA12,p_KcolA12,p_DA12,dsfA12,&
          p_KldB,p_KcolB,p_DB1,dsfB1,drelax)

      ! Solve A22 * u2 = f2 - A21*u1 - B2*p
      call spsor_aux_solveAAB(m,p_Du2,p_Du1,p_Dup,p_Df2,&
          p_KldA,p_KcolA,p_KdiagA,p_DA22,dsfA22,&
          p_KldA12,p_KcolA12,p_DA21,dsfA21,&
          p_KldB,p_KcolB,p_DB2,dsfB2,drelax)

    else
    
      ! A12/A21 are not present
      call spsor_aux_solveAB_2D(m,p_Du1,p_Du2,p_Dup,p_Df1,p_Df2,&
          p_KldA,p_KcolA,p_KdiagA,p_DA11,p_DA22,dsfA11,dsfA22,&
          p_KldB,p_KcolB,p_DB1,p_DB2,dsfB1,dsfB2,drelax)
      
    end if
    
    
    ! Step Two
    ! --------
    ! Calculate:
    !
    !  p := p + relax * S * (g - D*u)
    
    ! Is S a diagonal matrix?
    bdiag = (rdata%cpressure .eq. SPSOR_PRES_P0) .or. &
            (iand(rdata%cflags,SPSOR_FLAG_DIAG) .ne. 0_I32)
    if(bdiag) then

      ! Calculate p := p - relax * S * (g - D*u)
      call spsor_aux_multDDS(n,p_Du1,p_Du2,p_Dup,p_Dfp,&
          p_KldD,p_KcolD,p_DD1,p_DD2,dsfD1,dsfD2,p_DS,drelax)

    else
    
      ! Calculate w := g - D*u
      call spsor_aux_multDD(n,p_Du1,p_Du2,p_Dw,p_Dfp,&
          p_KldD,p_KcolD,p_DD1,p_DD2,dsfD1,dsfD2)

      ! Calculate p := p + relax * S * w
      call lsyssc_scalarMatVec(rdata%rS, rdata%rw, &
          rx%rvectorBlock(3), drelax, 1.0_DP)
    
    end if
    
    ! That's it
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine spsor_aux_solveAB_2D(n,Du1,Du2,Dup,Df1,Df2,KldA,KcolA,KdiagA,&
                                       DA1,DA2,dsfA1,dsfA2,KldB,KcolB,&
                                       DB1,DB2,dsfB1,dsfB2,drelax)

!<description>
  ! INTERNAL AUXILIARY ROUTINE:
  ! Solves:
  !              ( A1  0  ) * ( u ) = ( f ) - ( B1 ) * ( p )
  !              ( 0   A2 )   ( v )   ( g )   ( B2 )
!</description>

!<inputoutput>
  real(DP), dimension(*), intent(INOUT) :: Du1, Du2
!</inputoutput>

!<input>
  integer, intent(IN) :: n
  
  real(DP), dimension(*), intent(IN) :: Dup, Df1, Df2
  
  integer, dimension(*), intent(IN) :: KldA, KcolA, KdiagA
  real(DP), dimension(*), intent(IN) :: DA1, DA2
  real(DP), intent(IN) :: dsfA1, dsfA2

  integer, dimension(*), intent(IN) :: KldB, KcolB
  real(DP), dimension(*), intent(IN) :: DB1, DB2
  real(DP), intent(IN) :: dsfB1, dsfB2
  
  real(DP), intent(IN) :: drelax
!</input>

!</subroutine>

  integer :: i,j,k
  real(DP) :: dt,df,dg,daux1, daux2
  real(DP) :: dsf1,dsf2,dsb1,dsb2
  
    ! Pre-calculate scaling factors
    dsf1 = 1.0_DP / dsfA1
    dsf2 = 1.0_DP / dsfA2
    dsb1 = dsfB1 / dsfA1
    dsb2 = dsfB2 / dsfA2
  
    ! Loop over all rows
    do i = 1, n
    
      ! Fetch the rhs entries
      df = Df1(i)*dsf1
      dg = Df2(i)*dsf2
    
      ! Subtract B*p from rhs
      daux1 = 0.0_DP
      daux2 = 0.0_DP
      do j = KldB(i), KldB(i+1)-1
        dt = Dup(KcolB(j))
        daux1 = daux1 + DB1(j)*dt
        daux2 = daux2 + DB2(j)*dt
      end do
      df = df - dsb1*daux1
      dg = dg - dsb2*daux2
      
      ! Subtract A*u from rhs
      daux1 = 0.0_DP
      daux2 = 0.0_DP
      do j = KldA(i), KldA(i+1)-1
        k = KcolA(j)
        daux1 = daux1 + DA1(j)*Du1(k)
        daux2 = daux2 + DA2(j)*Du2(k)
      end do
      df = df - daux1
      dg = dg - daux2
      
      ! Solve
      k = KdiagA(i)
      Du1(i) = Du1(i) + drelax*df / DA1(k)
      Du2(i) = Du2(i) + drelax*dg / DA2(k)
    
    end do ! i
    
  end subroutine ! spsor_aux_solveAB_2D

!  ! ***************************************************************************
!
!!<subroutine>
!
!  pure subroutine spsor_aux_solveAB_3D(n,Du1,Du2,Du3,Dup,Df1,Df2,Df3,&
!                                       KldA,KcolA,KdiagA,DA1,DA2,DA3,&
!                                       dsfA1,dsfA2,dsfA3,KldB,KcolB,&
!                                       DB1,DB2,DB3,dsfB1,dsfB2,dsfB3,drelax)
!
!!<description>
!  ! INTERNAL AUXILIARY ROUTINE:
!  ! Solves:
!  !              ( A1  0  0  ) * ( u1 ) = ( f1 ) - ( B1 ) * ( p )
!  !              ( 0   A2 0  )   ( u2 )   ( f2 )   ( B2 )
!  !              ( 0   0  A3 )   ( u3 )   ( f3 )   ( B3 )
!!</description>
!
!!<inputoutput>
!  real(DP), dimension(*), intent(INOUT) :: Du1, Du2, Du3
!!</inputoutput>
!
!!<input>
!  integer, intent(IN) :: n
!  
!  real(DP), dimension(*), intent(IN) :: Dup, Df1, Df2, Df3
!  
!  integer, dimension(*), intent(IN) :: KldA, KcolA, KdiagA
!  real(DP), dimension(*), intent(IN) :: DA1, DA2, DA3
!  real(DP), intent(IN) :: dsfA1, dsfA2, dsfA3
!
!  integer, dimension(*), intent(IN) :: KldB, KcolB
!  real(DP), dimension(*), intent(IN) :: DB1, DB2, DB3
!  real(DP), intent(IN) :: dsfB1, dsfB2, dsfB3
!  
!  real(DP), intent(IN) :: drelax
!!</input>
!
!!</subroutine>
!
!  integer :: i,j,k
!  real(DP) :: dt,dfu1,dfu2,dfu3,daux1,daux2,daux3
!  real(DP) :: dsf1,dsf2,dsf3,dsb1,dsb2,dsb3
!  
!    ! Pre-calculate scaling factors
!    dsf1 = 1.0_DP / dsfA1
!    dsf2 = 1.0_DP / dsfA2
!    dsf3 = 1.0_DP / dsfA3
!    dsb1 = dsfB1 / dsfA1
!    dsb2 = dsfB2 / dsfA2
!    dsb3 = dsfB3 / dsfA3
!  
!    ! Loop over all rows
!    do i = 1, n
!    
!      ! Fetch the rhs entries
!      dfu1 = Df1(i)*dsf1
!      dfu2 = Df2(i)*dsf2
!      dfu3 = Df3(i)*dsf3
!
!      ! Subtract B*p from rhs
!      daux1 = 0.0_DP
!      daux2 = 0.0_DP
!      daux3 = 0.0_DP
!      do j = KldB(i), KldB(i+1)-1
!        dt = Dup(KcolB(j))
!        daux1 = daux1 + DB1(j)*dt
!        daux2 = daux2 + DB2(j)*dt
!        daux3 = daux3 + DB3(j)*dt
!      end do
!      dfu1 = dfu1 - dsb1*daux1
!      dfu2 = dfu2 - dsb2*daux2
!      dfu3 = dfu3 - dsb3*daux3
!      
!      ! Subtract A*u from rhs
!      daux1 = 0.0_DP
!      daux2 = 0.0_DP
!      daux2 = 0.0_DP
!      do j = KldA(i), KldA(i+1)-1
!        k = KcolA(j)
!        daux1 = daux1 + DA1(j)*Du1(k)
!        daux2 = daux2 + DA2(j)*Du2(k)
!        daux3 = daux3 + DA3(j)*Du3(k)
!      end do
!      dfu1 = dfu1 - daux1
!      dfu2 = dfu2 - daux2
!      dfu3 = dfu3 - daux3
!      
!      ! Solve
!      k = KdiagA(i)
!      Du1(i) = Du1(i) + drelax*dfu1 / DA1(k)
!      Du2(i) = Du2(i) + drelax*dfu2 / DA2(k)
!      Du3(i) = Du3(i) + drelax*dfu3 / DA3(k)
!    
!    end do ! i
!    
!  end subroutine ! spsor_aux_solveAB_3D

  ! ***************************************************************************

!<subroutine>

  pure subroutine spsor_aux_solveAAB(n,Du1,Du2,Dup,Df,KldA,KcolA,KdiagA,&
                                     DA,dsfA,KldA2,KcolA2,DA2,dsfA2,&
                                     KldB,KcolB,DB,dsfB,drelax)

!<description>
  ! INTERNAL AUXILIARY ROUTINE:
  ! Solves:
  !                A * u1 = f - A2 * u2 - B*p
!</description>

!<inputoutput>
  real(DP), dimension(*), intent(INOUT) :: Du1
!</inputoutput>

!<input>
  integer, intent(IN) :: n
  
  real(DP), dimension(*), intent(IN) :: Du2, Dup, Df
  
  integer, dimension(*), intent(IN) :: KldA, KcolA, KdiagA
  real(DP), dimension(*), intent(IN) :: DA
  real(DP), intent(IN) :: dsfA
  
  integer, dimension(*), intent(IN) :: KldA2, KcolA2
  real(DP), dimension(*), intent(IN) :: DA2
  real(DP), intent(IN) :: dsfA2

  integer, dimension(*), intent(IN) :: KldB, KcolB
  real(DP), dimension(*), intent(IN) :: DB
  real(DP), intent(IN) :: dsfB
  
  real(DP), intent(IN) :: drelax
!</input>

!</subroutine>

  integer :: i,j
  real(DP) :: df1,daux
  real(DP) :: dsa2,dsf,dsb
  
    ! Pre-calculate scaling factors
    dsf = 1.0_DP / dsfA
    dsa2 = dsfA2 / dsfA
    dsb = dsfB / dsfA
  
    ! Loop over all rows
    do i = 1, n
    
      ! Fetch the rhs entries
      df1 = Df(i)*dsf
    
      ! Subtract B*p from rhs
      daux = 0.0_DP
      do j = KldB(i), KldB(i+1)-1
        daux = daux + DB(j)*Dup(KcolB(j))
      end do
      df1 = df1 - dsb*daux
      
      ! Subtract A2*u2 from rhs
      daux = 0.0_DP
      do j = KldA2(i), KldA2(i+1)-1
        daux = daux + DA2(j)*Du2(KcolA2(j))
      end do
      df1 = df1 - dsa2*daux
      
      ! Subtract A*u1 from rhs
      daux = 0.0_DP
      do j = KldA(i), KldA(i+1)-1
        daux = daux + DA(j)*Du1(KcolA(j))
      end do
      df1 = df1 - daux
      
      ! Solve
      Du1(i) = Du1(i) + drelax*df1 / DA(KdiagA(i))
    
    end do ! i
    
  end subroutine ! spsor_aux_solveAAB

!  ! ***************************************************************************
!
!!<subroutine>
!
!  pure subroutine spsor_aux_solveAAAB(n,Du1,Du2,Du3,Dup,Df,KldA,KcolA,KdiagA,&
!                                     DA,dsfA,KldA2,KcolA2,DA2,DA3,dsfA2,dsfA3,&
!                                     KldB,KcolB,DB,dsfB,drelax)
!
!!<description>
!  ! INTERNAL AUXILIARY ROUTINE:
!  ! Solves:
!  !                A * u1 = f - A2 * u2 - A3 * u3 - B*p
!!</description>
!
!!<inputoutput>
!  real(DP), dimension(*), intent(INOUT) :: Du1
!!</inputoutput>
!
!!<input>
!  integer, intent(IN) :: n
!  
!  real(DP), dimension(*), intent(IN) :: Du2, Du3, Dup, Df
!  
!  integer, dimension(*), intent(IN) :: KldA, KcolA, KdiagA
!  real(DP), dimension(*), intent(IN) :: DA
!  real(DP), intent(IN) :: dsfA
!  
!  integer, dimension(*), intent(IN) :: KldA2, KcolA2
!  real(DP), dimension(*), intent(IN) :: DA2, DA3
!  real(DP), intent(IN) :: dsfA2, dsfA3
!
!  integer, dimension(*), intent(IN) :: KldB, KcolB
!  real(DP), dimension(*), intent(IN) :: DB
!  real(DP), intent(IN) :: dsfB
!  
!  real(DP), intent(IN) :: drelax
!!</input>
!
!!</subroutine>
!
!  integer :: i,j,k
!  real(DP) :: df1,daux,daux2
!  real(DP) :: dsa2,dsa3,dsf,dsb
!  
!    ! Pre-calculate scaling factors
!    dsf = 1.0_DP / dsfA
!    dsa2 = dsfA2 / dsfA
!    dsa3 = dsfA3 / dsfA
!    dsb = dsfB / dsfA
!  
!    ! Loop over all rows
!    do i = 1, n
!    
!      ! Fetch the rhs entries
!      df1 = Df(i)*dsf
!    
!      ! Subtract B*p from rhs
!      daux = 0.0_DP
!      do j = KldB(i), KldB(i+1)-1
!        daux = daux + DB(j)*Dup(KcolB(j))
!      end do
!      df1 = df1 - dsb*daux
!      
!      ! Subtract A2*u2 and A3*u3 from rhs
!      daux = 0.0_DP
!      daux2 = 0.0_DP
!      do j = KldA2(i), KldA2(i+1)-1
!        k = KcolA2(j)
!        daux  = daux  + DA2(j)*Du2(k)
!        daux2 = daux2 + DA3(j)*Du3(k)
!      end do
!      df1 = df1 - dsa2*daux - dsa3*daux2
!      
!      ! Subtract A*u1 from rhs
!      daux = 0.0_DP
!      do j = KldA(i), KldA(i+1)-1
!        daux = daux + DA(j)*Du1(KcolA(j))
!      end do
!      df1 = df1 - daux
!      
!      ! Solve
!      Du1(i) = Du1(i) + drelax*df1 / DA(KdiagA(i))
!    
!    end do ! i
!    
!  end subroutine ! spsor_aux_solveAAAB

  ! ***************************************************************************

!<subroutine>

  pure subroutine spsor_aux_multDDS(n,Du1,Du2,Dup,Df,KldD,KcolD,DD1,DD2,&
                                    dsfD1,dsfD2,DS,drelax)

!<description>
  ! INTERNAL AUXILIARY ROUTINE:
  ! Calculates:
  !
  !                p := p + relax * S * (f - D1 * u1 - D2 * u2)
  !
  ! Where S is a diagonal matrix
!</description>

!<inputoutput>
  real(DP), dimension(*), intent(INOUT) :: Dup
!</inputoutput>

!<input>
  integer, intent(IN) :: n
  
  real(DP), dimension(*), intent(IN) :: Du1, Du2, Df
  
  integer, dimension(*), intent(IN) :: KldD, KcolD
  real(DP), dimension(*), intent(IN) :: DD1, DD2
  real(DP), intent(IN) :: dsfD1, dsfD2

  real(DP), dimension(*), intent(IN) :: DS
  
  real(DP), intent(IN) :: drelax
!</input>

!</subroutine>

  integer :: i,j,k
  real(DP) :: daux1, daux2
  
    ! Loop over all rows
    !$omp parallel do private(j,k,daux1,daux2) if(n .gt. 1000)
    do i = 1, n
    
      ! Calculate D1*u1 and D2*u2
      daux1 = 0.0_DP
      daux2 = 0.0_DP
      do j = KldD(i), KldD(i+1)-1
        k = KcolD(j)
        daux1 = daux1 + DD1(j)*Du1(k)
        daux2 = daux2 + DD2(j)*Du2(k)
      end do
      
      ! Solve
      Dup(i) = Dup(i) + drelax*DS(i)*(Df(i) - dsfD1*daux1 - dsfD2*daux2)
    
    end do ! i
    !$omp end parallel do
    
  end subroutine ! spsor_aux_multDDS

!  ! ***************************************************************************
!
!!<subroutine>
!
!  pure subroutine spsor_aux_multDDDS(n,Du1,Du2,Du3,Dup,Df,KldD,KcolD,&
!                                     DD1,DD2,DD3,dsfD1,dsfD2,dsfD3,DS,drelax)
!
!!<description>
!  ! INTERNAL AUXILIARY ROUTINE:
!  ! Calculates:
!  !
!  !                p := p + relax * S * (f - D1 * u1 - D2 * u2 - D3 * u3)
!  !
!  ! Where S is a diagonal matrix
!!</description>
!
!!<inputoutput>
!  real(DP), dimension(*), intent(INOUT) :: Dup
!!</inputoutput>
!
!!<input>
!  integer, intent(IN) :: n
!  
!  real(DP), dimension(*), intent(IN) :: Du1, Du2, Du3, Df
!  
!  integer, dimension(*), intent(IN) :: KldD, KcolD
!  real(DP), dimension(*), intent(IN) :: DD1, DD2, DD3
!  real(DP), intent(IN) :: dsfD1, dsfD2, dsfD3
!
!  real(DP), dimension(*), intent(IN) :: DS
!  
!  real(DP), intent(IN) :: drelax
!!</input>
!
!!</subroutine>
!
!  integer :: i,j,k
!  real(DP) :: daux1, daux2, daux3
!  
!    ! Loop over all rows
!    !$omp parallel do private(j,k,daux1,daux2,daux3) if(n .gt. 1000)
!    do i = 1, n
!    
!      ! Calculate D1*u and D2*v
!      daux1 = 0.0_DP
!      daux2 = 0.0_DP
!      daux3 = 0.0_DP
!      do j = KldD(i), KldD(i+1)-1
!        k = KcolD(j)
!        daux1 = daux1 + DD1(j)*Du1(k)
!        daux2 = daux2 + DD2(j)*Du2(k)
!        daux3 = daux3 + DD3(j)*Du3(k)
!      end do
!      
!      ! Solve
!      Dup(i) = Dup(i) + drelax*DS(i)*(Df(i) - dsfD1*daux1 - dsfD2*daux2 - dsfD3*daux3)
!    
!    end do ! i
!    !$omp end parallel do
!    
!  end subroutine ! spsor_aux_multDDDS

  ! ***************************************************************************

!<subroutine>

  pure subroutine spsor_aux_multDD(n,Du1,Du2,Dw,Df,KldD,KcolD,DD1,DD2,&
                                   dsfD1,dsfD2)

!<description>
  ! INTERNAL AUXILIARY ROUTINE:
  ! Calculates:
  !
  !                w = f - D1 * u1 - D2 * u2
!</description>

!<inputoutput>
  real(DP), dimension(*), intent(INOUT) :: Dw
!</inputoutput>

!<input>
  integer, intent(IN) :: n
  
  real(DP), dimension(*), intent(IN) :: Du1, Du2, Df
  
  integer, dimension(*), intent(IN) :: KldD, KcolD
  real(DP), dimension(*), intent(IN) :: DD1, DD2
  real(DP), intent(IN) :: dsfD1, dsfD2
!</input>

!</subroutine>

  integer :: i,j,k
  real(DP) :: daux1, daux2
  
    ! Loop over all rows
    !$omp parallel do private(j,k,daux1,daux2) if(n .gt. 1000)
    do i = 1, n
    
      daux1 = 0.0_DP
      daux2 = 0.0_DP
      do j = KldD(i), KldD(i+1)-1
        k = KcolD(j)
        daux1 = daux1 + DD1(j)*Du1(k)
        daux2 = daux2 + DD2(j)*Du2(k)
      end do
      Dw(i) = Df(i) - dsfD1*daux1 - dsfD2*daux2
    
    end do ! i
    !$omp end parallel do
    
  end subroutine ! spsor_aux_multDD

!  ! ***************************************************************************
!
!!<subroutine>
!
!  pure subroutine spsor_aux_multDDD(n,Du1,Du2,Du3,Dw,Df,KldD,KcolD,&
!                                    DD1,DD2,DD3,dsfD1,dsfD2,dsfD3)
!
!!<description>
!  ! INTERNAL AUXILIARY ROUTINE:
!  ! Calculates:
!  !
!  !                w = f - D1 * u1 - D2 * u2 - D3 * u3
!!</description>
!
!!<inputoutput>
!  real(DP), dimension(*), intent(INOUT) :: Dw
!!</inputoutput>
!
!!<input>
!  integer, intent(IN) :: n
!  
!  real(DP), dimension(*), intent(IN) :: Du1, Du2, Du3, Df
!  
!  integer, dimension(*), intent(IN) :: KldD, KcolD
!  real(DP), dimension(*), intent(IN) :: DD1, DD2, DD3
!  real(DP), intent(IN) :: dsfD1, dsfD2, dsfD3
!!</input>
!
!!</subroutine>
!
!  integer :: i,j,k
!  real(DP) :: daux1, daux2, daux3
!  
!    ! Loop over all rows
!    !$omp parallel do private(j,k,daux1,daux2,daux3) if(n .gt. 1000)
!    do i = 1, n
!    
!      daux1 = 0.0_DP
!      daux2 = 0.0_DP
!      daux3 = 0.0_DP
!      do j = KldD(i), KldD(i+1)-1
!        k = KcolD(j)
!        daux1 = daux1 + DD1(j)*Du1(k)
!        daux2 = daux2 + DD2(j)*Du2(k)
!        daux3 = daux3 + DD3(j)*Du3(k)
!      end do
!      Dw(i) = Df(i) - dsfD1*daux1 - dsfD2*daux2 - dsfD3*daux3
!    
!    end do ! i
!    !$omp end parallel do
!    
!  end subroutine ! spsor_aux_multDDD

end module
