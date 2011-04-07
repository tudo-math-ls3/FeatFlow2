!##############################################################################
!# ****************************************************************************
!# <name> spsor </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module implements the saddle-point SOR ('SP-SOR') solver.
!#
!# The SP-SOR algorithm is a special solver for saddle-point systems,
!# i.e. systems of the following structure:
!#
!# <verb>
!#                        ( A  B ) * ( u ) = ( f )
!#                        ( D  0 )   ( p )   ( g )
!# </verb>
!#
!#
!# WARNING !
!# ---------
!# The SP-SOR solver works only for discontinous pressure spaces!
!# Although the 'diagonal' version of SP-SOR can also be defined for continous
!# pressure spaces, one should not expect that the solver will necessarily
!# converge.
!#
!#
!# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=\\
!# SP-SOR Algorithm Description\\
!# ----------------------------\\
!#
!# Basically, one SP-SOR iteration consists of two steps:
!#
!# Step 1: Update u:\\
!# - - - - - - - - -
!#
!# <tex> $$     u_{k+1} := u_k + P^{-1} * ( f - A * u_k - B * p_k )  $$ </tex>
!#
!# where P is the SOR-preconditioner matrix:
!#
!# <tex> $$     P := ( 1/omega * diagonal(A) + lower_triangular(A) )  $$ </tex>
!#
!#
!# Step 2: Update p:\\
!# - - - - - - - - -
!#
!# <tex> $$     p_{k+1} := p_k + omega * S * ( g - D * u_{k+1} )  $$ </tex>
!#
!# where S is an approximation of the inverse of the Schur-Complement of A:
!#
!# <tex> $$      S \approx ( -D * A^{-1} * B )^{-1}   $$ </tex>
!#
!#
!# TODO: document the matrix S
!#
!#
!#
!# The following routines can be found in this module:
!#
!# .) spsor_init
!#    -> Initialises a SP-SOR data structure.
!#
!# .) spsor_initNavSt2D
!#    -> Initialises a SP-SOR data structure for 2D Navier-Stokes systems.
!#
!# .) spsor_initData
!#    -> Performs data-dependent initialisation of the SP-SOR solver.
!#
!# .) spsor_solve
!#    -> Performs a number of iterations of the SP-SOR solver.
!#
!# .) spsor_precond
!#    -> Applies the SP-SOR / SP-SSOR preconditioner onto a defect vector.
!#
!# .) spsor_done
!#    -> Releases the SP-SOR data structure.
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
public :: spsor_init
public :: spsor_initNavSt2D
public :: spsor_initData
public :: spsor_done
public :: spsor_solve
public :: spsor_precond


!<constants>

!<constantblock description="Flags for the SP-SOR solver">

  ! Use symmetric preconditioner (SP-SSOR)
  integer(I32), parameter, public :: SPSOR_FLAG_SYM  = 1_I32
  
  ! Use diagonal version
  integer(I32), parameter, public :: SPSOR_FLAG_DIAG = 2_I32

!</constantblock>

!<constantblock description="Internal system identifiers">
  
  ! SP-SOR solver not initialised
  integer, parameter, public :: SPSOR_SYSTEM_NONE          = 0
  
  ! SP-SOR solver for 2D Navier-Stokes systems
  integer, parameter, public :: SPSOR_SYSTEM_NAVST2D       = 1
  
  ! SP-SOR solver for 3D Navier-Stokes systems
  !integer, parameter, public :: SPSOR_SYSTEM_NAVST3D       = 2
  
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
    
    ! Specifies whether S is a diagonal matrix.
    logical :: bdiagS = .false.
    
    ! Tempoary scalar vector of the size of the pressure space.
    ! Is only allocated if bdiagS is FALSE.
    type(t_vectorScalar) :: rw
  
  end type

!</typeblock>

!</types>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine spsor_init(rdata, rmatrix, csystem, cflags)

!<description>
  ! Initialises a SP-SOR data structure. This is just a wrapper for the
  ! system-specific initialisation routines.
!</description>

!<input>
  ! The system matrix of the linear system.
  type(t_matrixBlock), target, intent(in) :: rmatrix
  
  ! The system for which the SP-SOR solver is to be created.
  ! One of the SPSOR_SYSTEM_XXXX constants.
  integer, intent(in) :: csystem
  
  ! Optional: Flags for the solver.
  ! An OR-ed combination of SPSOR_FLAG_XXXX constants.
  integer(I32), optional, intent(in) :: cflags
!</input>

!<output>
  ! The SP-SOR data structure that is to be initialised.
  type(t_spsor), intent(out) :: rdata
!</output>

!</subroutine>

  integer(I32) :: cflag
  
    cflag = 0_I32
    if(present(cflags)) cflag = cflags
    
    ! Which system?
    select case(csystem)
    case(SPSOR_SYSTEM_NAVST2D)
      ! 2D Navier-Stokes system
      call spsor_initNavSt2D(rdata, rmatrix, cflag)
      
    case default
      call output_line ('Invalid system identifier',&
          OU_CLASS_ERROR,OU_MODE_STD,'spsor_init')
      call sys_halt()
      
    end select

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine spsor_initData(rdata)

!<description>
  ! This subroutine performs data-dependent initialisation of the SP-SOR
  ! solver. The data structure must have already been set up by spsor_init().
!</description>

!<inputoutput>
  ! The SP-SOR data structure
  type(t_spsor), intent(inout) :: rdata
!</inputoutput>

!</subroutine>

    ! Okay, which SP-SOR solver do we have here?
    select case(rdata%csystem)
    case (SPSOR_SYSTEM_NAVST2D)
      
      ! SP-SOR for 2D Navier-Stokes systems
      call spsor_initData_NavSt2D(rdata)

    case default
      call output_line ('SP-SOR data structure not initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'spsor_initData')
      call sys_halt()

    end select
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spsor_solve(rdata, rsol, rrhs, niterations, domega)

!<description>
  ! This routine performs a specified number of SP-SOR iterations.
!</description>

!<input>
  ! The SP-SOR data structure
  type(t_spsor), intent(in) :: rdata
  
  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs
  
  ! The number of iterations that are to be performed
  integer, intent(in) :: niterations
  
  ! Optional: The relaxation parameter in range (0,2).
  ! If not given, 1.0 is used.
  real(DP), optional, intent(in) :: domega
!</input>

!<inputoutput>
  ! The iteration vector that is to be updated.
  type(t_vectorblock), intent(inout) :: rsol
!</inputoutput>

!</subroutine>

  real(DP) :: dom
  
    dom = 1.0_DP
    if(present(domega)) dom = domega

    ! Make sure we do not apply the symmetric variant ('SP-SSOR')
    if(iand(rdata%cflags, SPSOR_FLAG_SYM) .ne. 0_I32) then
      call output_line ('SP-SSOR cannot be called via spsor_solve',&
          OU_CLASS_ERROR,OU_MODE_STD,'spsor_solve')
      call sys_halt()
    end if
    
    ! Nothing to do?
    if(niterations .le. 0) return

    ! Okay, which SP-SOR solver do we have here?
    select case(rdata%csystem)
    case (SPSOR_SYSTEM_NAVST2D)
    
      ! Call the routine for 2D Navier-Stokes systems
      call spsor_solve_NavSt2D(rdata, rsol, rrhs, niterations, dom)
      
    case default
      call output_line ('SP-SOR data structure not initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'spsor_solve')
      call sys_halt()
    
    end select
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spsor_precond(rdata, rdef, domega)

!<description>
  ! This routine performs applies the SP-SOR / SP-SSOR preconditioner onto
  ! a defect vector.
!</description>

!<input>
  ! The SP-SOR data structure
  type(t_spsor), intent(in) :: rdata
  
  ! Optional: The relaxation parameter in range (0,2).
  ! If not given, 1.0 is used.
  real(DP), optional, intent(in) :: domega
!</input>

!<inputoutput>
  ! The defect vector that is to be preconditioned
  type(t_vectorblock), intent(inout) :: rdef
!</inputoutput>

!</subroutine>

  real(DP) :: dom
  
    dom = 1.0_DP
    if(present(domega)) dom = domega

    ! Okay, which SP-SOR solver do we have here?
    select case(rdata%csystem)
    case (SPSOR_SYSTEM_NAVST2D)
    
      ! Do we apply SP-SOR or SP-SSOR?
      if(iand(rdata%cflags, SPSOR_FLAG_SYM) .ne. 0_I32) then
        ! Apply SP-SSOR for 2D Navier-Stokes systems
        call spsor_prec_sym_NavSt2D(rdata, rdef, dom)
      else
        ! Apply SP-SOR for 2D Navier-Stokes systems
        call spsor_prec_NavSt2D(rdata, rdef, dom)
      end if
      
    case default
      call output_line ('SP-SOR data structure not initialised',&
          OU_CLASS_ERROR,OU_MODE_STD,'spsor_precond')
      call sys_halt()
    
    end select
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spsor_done(rdata)

!<description>
  ! Releases a SP-SOR data structure and frees all allocated memory.
!</description>

!<inputoutput>
  ! The SP-SOR data structure that is to be destroyed.
  type(t_spsor), intent(inout) :: rdata
!</inputoutput>

!</subroutine>

    ! Release Schur-complement matrix
    if(lsyssc_hasMatrixStructure(rdata%rS)) &
      call lsyssc_releaseMatrix(rdata%rS)
    
    ! Release temporary vector
    if(rdata%rw%NEQ .gt. 0) &
      call lsyssc_releaseVector(rdata%rw)
      
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spsor_asmMatStructS(rS, rD, rB)

!<description>
  ! INTERNAL AUXILIARY ROUTINE:
  ! Assembles the matrix structure of the S matrix.
!</description>

!<input>
  ! The matrix D
  type(t_matrixScalar), intent(in) :: rD

  ! The matrix B
  type(t_matrixScalar), intent(in) :: rB
!</input>

!<output>
  ! The matrix S
  type(t_matrixScalar), intent(out) :: rS
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
    call storage_new ('spsor_asmMatStructS', 'Kld', ndofp+1, ST_INT,&
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
      do i = 1, ndegB
        if(Imask(i) .eq. ind) then
          k = k+1
          Idofs(k,idofp) = Iadj(i,1)
        end if
      end do
      
      ! Set up row pointer for next row
      rS%NA = rS%NA + k
      p_KldS(idofp+1) = p_KldS(idofp) + k
    
    end do ! idofp
    
    ! Now we know the number of non-zeroes in S, so let us allocate column index
    ! array.
    call storage_new ('spsor_asmMatStructS', 'Kcol', rS%NA, ST_INT,&
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
    call storage_new ('spsor_asmMatStructS', 'Kdiagonal', ndofp, ST_INT,&
                      rS%h_Kdiagonal, ST_NEWBLOCK_NOINIT)
    call storage_getbase_int(rS%h_Kdiagonal, p_KdiagS)
    
    ! Build the diagonal pointer array
    call lsyssc_rebuildKdiagonal(p_KcolS,p_KldS,p_KdiagS,ndofp)
    
    ! That is it

  end subroutine
  

  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !
  !  R O U T I N E S   F O R   2 D   N A V I E R - S T O K E S   S Y S T E M S
  !
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


!<subroutine>

  subroutine spsor_initNavSt2D(rdata, rmatrix, cflags)

!<description>
  ! Initialises a SP-SOR data structure for 2D Navier-Stokes systems.
!</description>

!<input>
  ! The system matrix of the linear system.
  type(t_matrixBlock), target, intent(in) :: rmatrix
  
  ! Flags for the SP-SOR solver.
  integer(I32), intent(in) :: cflags
!</input>

!<output>
  ! The SP-SOR data structure that is to be initialised.
  type(t_spsor), intent(out) :: rdata
!</output>

!</subroutine>

  ! Some local variables
  type(t_spatialDiscretisation), pointer :: p_rdiscrV, p_rdiscrP
  integer(I32) :: celemP, celemV
  integer :: ndofPres, ndofVelo, ieldist

    ! Set the system and flags
    rdata%csystem = SPSOR_SYSTEM_NAVST2D
    rdata%cflags = cflags
    
    ! Check if the matrix is valid
    call spsor_checkMatrixNavSt2D(rmatrix)
    
    ! Hang in the matrix pointers
    rdata%p_rA11 => rmatrix%RmatrixBlock(1,1)
    rdata%p_rA22 => rmatrix%RmatrixBlock(2,2)
    rdata%p_rB1 => rmatrix%RmatrixBlock(1,3)
    rdata%p_rB2 => rmatrix%RmatrixBlock(2,3)
    rdata%p_rD1 => rmatrix%RmatrixBlock(3,1)
    rdata%p_rD2 => rmatrix%RmatrixBlock(3,2)
    
    ! Are the A12/A21 sub-matrices present?
    if(rmatrix%RmatrixBlock(1,2)%NEQ .gt. 0) then
      rdata%p_rA12 => rmatrix%RmatrixBlock(1,2)
      rdata%p_rA21 => rmatrix%RmatrixBlock(2,1)
    end if
    
    ! Is the C sub-matrix present?
    if(rmatrix%RmatrixBlock(3,3)%NEQ .gt. 0) &
      rdata%p_rC => rmatrix%RmatrixBlock(3,3)
    
    ! Get the total number of dofs
    ndofVelo = rdata%p_rD1%NCOLS
    ndofPres = rdata%p_rD1%NEQ
    
    ! Get the discretisation of the velocity spaces
    p_rdiscrV => rdata%p_rB1%p_rspatialDiscrTest
    
    if(associated(p_rdiscrV)) then
    
      rdata%cvelocity = SPSOR_VELO_NONE
    
      ! Loop over all element distributions
      do ieldist = 1, p_rdiscrV%inumFESpaces
      
        ! Get the element of the spaces
        celemV = elem_getPrimaryElement(&
                 p_rdiscrV%RelementDistr(ieldist)%celement)

        ! What is the velocity element?
        select case(celemV)
        case(EL_P1T_2D)
        
          select case(rdata%cvelocity)
          case(SPSOR_VELO_NONE, SPSOR_VELO_P1T)
            rdata%cvelocity = SPSOR_VELO_P1T
          
          case default
            rdata%cvelocity = SPSOR_VELO_ARBIT
          end select
          
        case(EL_Q1T_2D)
        
          select case(rdata%cvelocity)
          case(SPSOR_VELO_NONE, SPSOR_VELO_Q1T)
            rdata%cvelocity = SPSOR_VELO_Q1T
          
          case default
            rdata%cvelocity = SPSOR_VELO_ARBIT
          end select
          
        case(EL_Q1TB_2D)
        
          select case(rdata%cvelocity)
          case(SPSOR_VELO_NONE, SPSOR_VELO_Q1TB)
            rdata%cvelocity = SPSOR_VELO_Q1TB
          
          case default
            rdata%cvelocity = SPSOR_VELO_ARBIT
            
          end select
        
        case default
          rdata%cvelocity = SPSOR_VELO_ARBIT
          
        end select
        
      end do ! ieldist
    
    else
    
      ! No velocity discretisation available
      rdata%cvelocity = SPSOR_VELO_ARBIT
      
    end if
    
    ! Get the discretisation of the pressure space
    p_rdiscrP => rdata%p_rB1%p_rspatialDiscrTrial
    
    if(associated(p_rdiscrP)) then
    
      rdata%cpressure = SPSOR_PRES_NONE
      
      ! Loop over all element distributions
      do ieldist = 1, p_rdiscrP%inumFESpaces
      
        ! Get the element of the spaces
        celemP = elem_getPrimaryElement(&
                 p_rdiscrP%RelementDistr(ieldist)%celement)
        
        ! What is the pressure element?
        select case(celemP)
        case(EL_P0_2D, EL_Q0_2D)

          select case(rdata%cpressure)
          case(SPSOR_PRES_NONE, SPSOR_PRES_P0)
            rdata%cpressure = SPSOR_PRES_P0
          
          case default
            rdata%cpressure = SPSOR_PRES_ARBIT
            
          end select
            
        case default
          rdata%cpressure = SPSOR_PRES_ARBIT
          
        end select
        
      end do ! ieldist
    
    else
      
      ! No pressure discretisation available
      rdata%cpressure = SPSOR_PRES_ARBIT
      
    end if
    
    ! Now we need to set up the structure of the S matrix.
    
    ! Will S be a diagonal matrix?
    rdata%bdiagS = (rdata%cpressure .eq. SPSOR_PRES_P0) .or. &
                   (iand(rdata%cflags,SPSOR_FLAG_DIAG) .ne. 0_I32)
    
    if(rdata%bdiagS) then
    
      ! Create a diagonal matrix
      call lsyssc_createDiagMatrixStruc(rdata%rS, ndofPres, &
                                        LSYSSC_MATRIX9)
    
    else
    
      ! Call the auxiliary routine to assemble the structure of S.
      call spsor_asmMatStructS(rdata%rS, rdata%p_rD1, rdata%p_rB1)
      
      ! Allocate a temporary vector
      call lsyssc_createVector (rdata%rw, ndofPres, .false.)
    
    end if
    
    ! Allocate an empty matrix
    call lsyssc_allocEmptyMatrix (rdata%rS, LSYSSC_SETM_UNDEFINED)

    ! That is it

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spsor_checkMatrixNavSt2D(rmatrix, bvalid)

!<description>
  ! Checks whether a given block matrix is valid for the 2D Navier-Stokes
  ! SP-SOR solver.
!</description>

!<input>
  ! The system matrix that is to be checked.
  type(t_matrixBlock), intent(in) :: rmatrix
!</input>

!<output>
  ! OPTIONAL: Recieves whether the matrix is valid.
  ! If the given matrix cannot be handled, then:
  !   if given, bvalid is set to FALSE.
  !   if not given, this routine prints an error and aborts program execution.
  logical, optional, intent(out) :: bvalid
!</output>

!</subroutine>

  ! local variables
  integer :: i,j
  logical :: bHaveA12, bHaveAij
  
    bHaveA12 = .false.
  
    ! Assume that the matrix is invalid
    if(present(bvalid)) bvalid = .false.
  
    ! First of all, the matrix must be 3x3
    if((rmatrix%nblocksPerRow .ne. 3) .or. (rmatrix%nblocksPerCol .ne. 3)) then
      if(present(bvalid)) then
        return
      else
        call output_line ('Matrix must be 3x3',&
            OU_CLASS_ERROR,OU_MODE_STD,'spsor_checkMatrixNavSt2D')
        call sys_halt()
      end if
    end if
    
    ! Check all matrices
    do i = 1, 3
      do j = 1, 3
      
        ! Skip the block if it is empty
        if(rmatrix%RmatrixBlock(i,j)%NEQ .le. 0) cycle

        ! Make sure it is a type-9 matrix
        if(rmatrix%RmatrixBlock(i,j)%cmatrixFormat .ne. LSYSSC_MATRIX9) then
          if(present(bvalid)) then
            return
          else
            call output_line ('Matrices must be of format LSYSSC_MATRIX9',&
                OU_CLASS_ERROR,OU_MODE_STD,'spsor_checkMatrixNavSt2D')
            call sys_halt()
          end if
        end if
        
        ! Make sure it is not virtually transposed
        if(iand(rmatrix%RmatrixBlock(i,j)%imatrixSpec,&
                LSYSSC_MSPEC_TRANSPOSED) .ne. 0) then
          if(present(bvalid)) then
            return
          else
            call output_line ('Matrices format must not be transposed',&
                OU_CLASS_ERROR,OU_MODE_STD,'spsor_checkMatrixNavSt2D')
            call sys_halt()
          end if
        end if

      end do ! j
    end do ! i
    
    ! Now loop over the A matrices
    do i = 1, 2
      do j = 1, 2
        if(i .eq. j) then
          
          ! Main diagonal block
          
          ! Make sure the matrix exists
          if((.not. lsyssc_hasMatrixStructure(rmatrix%RmatrixBlock(i,j))) .or. &
             (abs(rmatrix%RmatrixBlock(i,j)%dscaleFactor) .le. SYS_EPSREAL_DP)) then
            if(present(bvalid)) then
              return
            else
              call output_line ('Sub-Matrix A(' // trim(sys_sil(i,2)) // ',' // &
                  trim(sys_sil(j,2)) // ') does not exist',&
                  OU_CLASS_ERROR,OU_MODE_STD,'spsor_checkMatrixNavSt2D')
              call sys_halt()
            end if
          end if
          
          if(i .gt. 1) then
            ! Make all main diagonal blocks have the same structure
            if(.not. lsyssc_isMatrixStructureShared(&
                rmatrix%RmatrixBlock(1,1), rmatrix%RmatrixBlock(i,i))) then
              if(present(bvalid)) then
                return
              else
                call output_line ('All A(i,i) sub-matrices must share the same structure',&
                    OU_CLASS_ERROR,OU_MODE_STD,'spsor_checkMatrixNavSt2D')
                call sys_halt()
              end if
            end if
          end if
          
        else
        
          ! Off-diagonal block
          
          if((i .eq. 1) .and. (j .eq. 2)) then
            
            ! Block A12 - does it exist?
            bHaveA12 = lsyssc_hasMatrixStructure(rmatrix%RmatrixBlock(i,j)) .and. &
                      (abs(rmatrix%RmatrixBlock(i,j)%dscaleFactor) .gt. SYS_EPSREAL_DP)
            
          else

            ! Off-diagonal block other than A12 - does it exist?
            bHaveAij = lsyssc_hasMatrixStructure(rmatrix%RmatrixBlock(i,j)) .and. &
                      (abs(rmatrix%RmatrixBlock(i,j)%dscaleFactor) .gt. SYS_EPSREAL_DP)

            !if(bHaveAij .xor. bHaveA12) then
            if((bHaveAij .and. (.not. bHaveA12)) .or. ((.not. bHaveAij) .and. bHaveA12)) then
              ! Either A12 does not exist, but Aij does - or the other way
              if(present(bvalid)) then
                return
              else
                call output_line ('Only either all or none off-diagonal A(i,j) may exist',&
                    OU_CLASS_ERROR,OU_MODE_STD,'spsor_checkMatrixNavSt2D')
                call sys_halt()
              end if
            end if
          end if
        end if
      end do ! j
    end do ! i
    
    ! Loop over the B matrices
    do i = 1, 2
    
      ! Make sure the matrix exists
      if((.not. lsyssc_hasMatrixStructure(rmatrix%RmatrixBlock(i,3))) .or. &
         (abs(rmatrix%RmatrixBlock(i,3)%dscaleFactor) .le. SYS_EPSREAL_DP)) then
        if(present(bvalid)) then
          return
        else
          call output_line ('Sub-Matrix B(' // trim(sys_sil(i,2)) // ') does not exist',&
              OU_CLASS_ERROR,OU_MODE_STD,'spsor_checkMatrixNavSt2D')
          call sys_halt()
        end if
      end if
      
!      if(i .gt. 1) then
!        ! Make all B matrices have the same structure
!        if(.not. lsyssc_isMatrixStructureShared(&
!            rmatrix%RmatrixBlock(1,3), rmatrix%RmatrixBlock(i,3))) then
!          if(present(bvalid)) then
!            return
!          else
!            call output_line ('All B(i) sub-matrices must share the same structure',&
!                OU_CLASS_ERROR,OU_MODE_STD,'spsor_checkMatrixNavSt2D')
!            call sys_halt()
!          end if
!        end if
!      end if
    end do ! i

    ! Loop over the D matrices
    do j = 1, 2
    
      ! Make sure the matrix exists
      if((.not. lsyssc_hasMatrixStructure(rmatrix%RmatrixBlock(3,j))) .or. &
         (abs(rmatrix%RmatrixBlock(3,j)%dscaleFactor) .le. SYS_EPSREAL_DP)) then
        if(present(bvalid)) then
          return
        else
          call output_line ('Sub-Matrix D(' // trim(sys_sil(j,2)) // ') does not exist',&
              OU_CLASS_ERROR,OU_MODE_STD,'spsor_checkMatrixNavSt2D')
          call sys_halt()
        end if
      end if
      
!      if(j .gt. 1) then
!        ! Make all D matrices have the same structure
!        if(.not. lsyssc_isMatrixStructureShared(&
!            rmatrix%RmatrixBlock(3,1), rmatrix%RmatrixBlock(3,j))) then
!          if(present(bvalid)) then
!            return
!          else
!            call output_line ('All D(j) sub-matrices must share the same structure',&
!                OU_CLASS_ERROR,OU_MODE_STD,'spsor_checkMatrixNavSt2D')
!            call sys_halt()
!          end if
!        end if
!      end if
    end do ! j
    
    ! All checks passed, so the matrix is valid
    if(present(bvalid)) bvalid = .true.
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine spsor_initData_NavSt2D(rdata)

!<description>
  ! This subroutine performs data-dependent initialisation of the SP-SOR
  ! solver for 2D Navier-Stokes systems. 
!</description>

!<inputoutput>
  ! The SP-SOR data structure
  type(t_spsor), intent(inout) :: rdata
!</inputoutput>

!</subroutine>

  ! local variables
  logical :: bdone, bHaveA12
  

    ! Do we use the 'diagonal' or the 'full' SP-SOR variant?
    if(iand(rdata%cflags,SPSOR_FLAG_DIAG) .ne. 0_I32) then

      ! Call the diagonal version
      call spsor_initData_NavSt2D_diag(rdata)
    
    else
    
      ! Do we use only P0/Q0 for pressure?
      ! We have some fast implementations for this case...
      select case(rdata%cpressure)
      case (SPSOR_PRES_P0)
      
        bdone = .false.
        
        ! Do we have the off-diagonal A-matrices?
        bHaveA12 = associated(rdata%p_rA12)
      
        ! Do we even have one of our special velocity elements?
        ! This would speed up the whole process even more...
        select case(rdata%cvelocity)
!        case (SPSOR_VELO_P1T)
!          todo

        case (SPSOR_VELO_Q1T)
          
          ! Q1~/Q0 discretisation
          if(.not. bHaveA12) then
          
            ! Okay, we have a special implementation for this case...
            call spsor_initData_NavSt2D_Q1T_Q0(rdata)
            bdone = .true.
          
          end if
          
!        case (SPSOR_VELO_Q1TB)
!          todo
        end select
        
        ! If none of the above special cases applied, we will take care of
        ! the initialisation now...
        if(.not. bdone) then
        
          ! P0/Q0 pressure, arbitrary velocity
          call spsor_initData_NavSt2D_P0(rdata)

        end if
      
      case default
    
        ! arbitrary pressure, arbitrary velocity
        call spsor_initData_NavSt2D_full(rdata)

      end select
      
    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spsor_initData_NavSt2D_diag(rdata)

!<description>
  ! Internal subroutine:
  ! Peforms data-dependent initialisation of the 'diagonal' SP-SOR algorithm
  ! for 2D Navier-Stokes systems: arbitrary pressure, arbitrary velocity
!</description>

!<inputoutput>
  ! The SP-SOR data structure
  type(t_spsor), intent(inout) :: rdata
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

    ! Let us assume we do not have the optional matrices
    bHaveC = .false.
    
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

    ! Let us loop over all rows of S
    do idofp = 1, rdata%rS%NEQ
    
      ! Format the local matrices
      dc = 0.0_DP
      
      ! Format the Schur-complement entry for this pressure DOF
      p_DS(idofp) = 0.0_DP

      ! Fetch the number of velocity DOFs adjacent to this pressure DOF
      ndofV = p_KldD(idofp+1) - p_KldD(idofp)
      
      ! If the C matrix exists, grab it is main diagonal entry
      if(bHaveC) then
        dc = dsfC*p_DC(p_KdiagC(idofp))
      end if
      
      ! Let us loop over all velocity DOFs which are adjacent to the current
      ! pressure DOF.
      daux1 = 0.0_DP
      daux2 = 0.0_DP
      do j = p_KldD(idofp), p_KldD(idofp+1)-1
      
        ! Get the index of the velocity DOF
        idofu = p_KcolD(j)
        
        ! Fetch the main diagonal index of A
        k = p_KdiagA(idofu)

        ! Try to find the corrent entry in B1/B2
        do i = p_KldB(idofu), p_KldB(idofu+1)-1
          if(p_KcolB(i) .eq. idofp) then
            daux1 = daux1 + (p_DB1(i)*p_DD1(j)) / p_DA11(k)
            daux2 = daux2 + (p_DB2(i)*p_DD2(j)) / p_DA22(k)
            exit
          end if
        end do ! i
      
      end do ! id1
      
      ! Okay, let us calculate the Schur-complement dc = C - D * A^-1 * B
      dc = dc - dsf1*daux1 - dsf2*daux2
      
      ! Now if the Schur-complement matrix is regular, we will store the inverse
      ! in the corresponding entry in p_DS.
      if(abs(dc) .gt. SYS_EPSREAL_DP) p_DS(idofp) = 1.0_DP / dc
      
    end do ! idofp
    
    ! That is it
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spsor_initData_NavSt2D_Q1T_Q0(rdata)

!<description>
  ! Internal subroutine:
  ! Peforms data-dependent initialisation of the 'full' SP-SOR algorithm
  ! for 2D Navier-Stokes systems: Q0 pressure, Q1~ velocity
!</description>

!<inputoutput>
  ! The SP-SOR data structure
  type(t_spsor), intent(inout) :: rdata
!</inputoutput>

!</subroutine>

  ! Multiplication factors
  real(DP) :: dsfA11, dsfA22, dsfB1, dsfB2, dsfD1, dsfD2, dsfC
  
  ! Quick access for the matrix arrays
  integer, dimension(:), pointer :: p_KldA,p_KldB,p_KldC,p_KldD,&
      p_KcolA,p_KcolB,p_KcolC,p_KcolD,p_KdiagC
  real(DP), dimension(:), pointer :: p_DA11,p_DA22,p_DB1,p_DB2,&
      p_DD1,p_DD2,p_DC,p_DS

  ! Temporary local matrices
  real(DP), dimension(4,4) :: DA1, DA2, DI1, DI2
  real(DP), dimension(4) :: DB1, DB2, DIB1, DIB2
  
  ! local variables
  logical :: bHaveC
  real(DP) :: dc,daux1,daux2
  integer :: idofp,idofu,i,j1,j2,k,id1,id2
  logical :: bsuccess1,bsuccess2

    ! Let us assume we do not have the optional matrices
    bHaveC = .false.
    
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

    ! Let us loop over all rows of S
    do idofp = 1, rdata%rS%NEQ
    
      ! Format the local matrices
      DA1 = 0.0_DP
      DA2 = 0.0_DP
      DB1 = 0.0_DP
      DB2 = 0.0_DP
      dc = 0.0_DP
      
      ! Format the Schur-complement entry for this pressure DOF
      p_DS(idofp) = 0.0_DP

      ! If the C matrix exists, grab it is main diagonal entry
      if(bHaveC) then
        dc = dsfC*p_DC(p_KdiagC(idofp))
      end if
      
      ! Let us loop over all velocity DOFs which are adjacent to the current
      ! pressure DOF.
      do id1 = p_KldD(idofp), p_KldD(idofp+1)-1
      
        ! Get the index of the velocity DOF
        idofu = p_KcolD(id1)

        ! Calculate the local index of the velocity DOF
        j1 = id1 - p_KldD(idofp) + 1
        
        ! Let us fetch the local A11/A22 matrices
        do i = p_KldA(idofu), p_KldA(idofu+1)-1
        
          ! Get the column index
          k = p_KcolA(i)
          
          ! Let us see if this corresponds to one of our local velocity DOFs
          do id2 = p_KldD(idofp), p_KldD(idofp+1)-1
            if(k .eq. p_KcolD(id2)) then
              ! Okay, incorporate the entries into the local matrix
              j2 = id2 - p_KldD(idofp) + 1
              DA1(j1,j2) = dsfA11*p_DA11(i)
              DA2(j1,j2) = dsfA22*p_DA22(i)
              exit
            end if
          end do ! id2
        end do ! i

        ! Let us fetch the local B matrices
        do i = p_KldB(idofu), p_KldB(idofu+1)-1
          if(p_KcolB(i) .eq. idofP) then
            DB1(j1) = dsfB1*p_DB1(i)
            DB2(j1) = dsfB2*p_DB2(i)
            exit
          end if
        end do ! i
      
      end do ! id1
      
      ! Invert the A matrices
      call mprim_invert4x4MatrixDirectDble(DA1, DI1,bsuccess1)
      call mprim_invert4x4MatrixDirectDble(DA2, DI2,bsuccess2)
      
      if (bsuccess1 .and. bsuccess2) then
      
        ! Calculate A^{-1} * B
        DIB1 = matmul(DI1, DB1)
        DIB2 = matmul(DI2, DB2)
        
        ! Okay, let us calculate the Schur-complement dc = C - D * A^-1 * B
        daux1 = 0.0_DP
        daux2 = 0.0_DP
        k = 1
        do id1 = p_KldD(idofp), p_KldD(idofp+1)-1
          daux1 = daux1 + p_DD1(id1)*DIB1(k)
          daux2 = daux2 + p_DD2(id1)*DIB2(k)
          k = k+1
        end do
        dc = dc - dsfD1*daux1 - dsfD2*daux2
        
        ! Now if the Schur-complement matrix is regular, we will store the inverse
        ! in the corresponding entry in p_DS.
        if(abs(dc) .gt. SYS_EPSREAL_DP) p_DS(idofp) = 1.0_DP / dc
        
      end if
      
    end do ! idofp
    
    ! That is it
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spsor_initData_NavSt2D_P0(rdata)

!<description>
  ! Internal subroutine:
  ! Peforms data-dependent initialisation of the 'full' SP-SOR algorithm
  ! for 2D Navier-Stokes systems: P0/Q0 pressure, arbitrary velocity
!</description>

!<inputoutput>
  ! The SP-SOR data structure
  type(t_spsor), intent(inout) :: rdata
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

    ! Let us assume we do not have the optional matrices
    bHaveA12 = .false.
    bHaveC = .false.
    
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
    
    ! Okay, let us allocate the temporary data
    allocate(DA(nmaxdofV,nmaxdofV))
    allocate(DB(nmaxdofV))
    allocate(Ipivot(nmaxdofV))
    
    ! Let us loop over all rows of S
    do idofp = 1, rdata%rS%NEQ
    
      ! Format the local matrices
      DA = 0.0_DP
      DB = 0.0_DP
      dc = 0.0_DP
      
      ! Format the Schur-complement entry for this pressure DOF
      p_DS(idofp) = 0.0_DP

      ! Fetch the number of velocity DOFs adjacent to this pressure DOF
      ndofV = p_KldD(idofp+1) - p_KldD(idofp)
      
      ! If the C matrix exists, grab it is main diagonal entry
      if(bHaveC) then
        dc = dsfC*p_DC(p_KdiagC(idofp))
      end if
      
      ! Let us loop over all velocity DOFs which are adjacent to the current
      ! pressure DOF.
      do id1 = p_KldD(idofp), p_KldD(idofp+1)-1
      
        ! Get the index of the velocity DOF
        idofu = p_KcolD(id1)

        ! Calculate the local index of the velocity DOF
        j1 = id1 - p_KldD(idofp) + 1
        
        ! Let us fetch the local A11/A22 matrices
        do i = p_KldA(idofu), p_KldA(idofu+1)-1
        
          ! Get the column index
          k = p_KcolA(i)
          
          ! Let us see if this corresponds to one of our local velocity DOFs
          do id2 = p_KldD(idofp), p_KldD(idofp+1)-1
            if(k .eq. p_KcolD(id2)) then
              ! Okay, incorporate the entries into the local matrix
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
          ! Let us fetch the local A12/A21 matrices
          do i = p_KldA12(idofu), p_KldA12(idofu+1)-1
          
            ! Get the column index
            k = p_KcolA12(i)
            
            ! Let us see if this corresponds to one of our local velocity DOFs
            do id2 = p_KldD(idofp), p_KldD(idofp+1)-1
              if(k .eq. p_KcolD(id2)) then
                ! Okay, incorporate the entries into the local matrix
                j2 = id2 - p_KldD(idofp) + 1
                DA(      j1,ndofV+j2) = dsfA12*p_DA12(i)
                DA(ndofV+j1,      j2) = dsfA21*p_DA21(i)
                exit
              end if
            end do ! id2
          end do ! i
        end if

        ! Let us fetch the local B matrices
        do i = p_KldB(idofu), p_KldB(idofu+1)-1
          if(p_KcolB(i) .eq. idofP) then
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
      
      ! Okay, let us calculate the Schur-complement dc = C - D * A^-1 * B
      daux1 = 0.0_DP
      daux2 = 0.0_DP
      k = 1
      do id1 = p_KldD(idofp), p_KldD(idofp+1)-1
        daux1 = daux1 + p_DD1(id1)*DB(k)
        daux2 = daux2 + p_DD2(id1)*DB(ndofV+k)
        k = k+1
      end do
      dc = dc - dsfD1*daux1 - dsfD2*daux2
      
      ! Now if the Schur-complement matrix is regular, we will store the inverse
      ! in the corresponding entry in p_DS.
      if(abs(dc) .gt. SYS_EPSREAL_DP) p_DS(idofp) = 1.0_DP / dc
      
    end do ! idofp
    
    ! And release the temporary memory
    deallocate(Ipivot)
    deallocate(DB)
    deallocate(DA)
    
    ! That is it
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spsor_initData_NavSt2D_full(rdata)

!<description>
  ! Internal subroutine:
  ! Peforms data-dependent initialisation of the 'full' SP-SOR algorithm
  ! for 2D Navier-Stokes systems: arbitrary pressure, arbitrary velocity
!</description>

!<inputoutput>
  ! The SP-SOR data structure
  type(t_spsor), intent(inout) :: rdata
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
  
    ! Let us assume we do not have the optional matrices
    bHaveA12 = .false.
    bHaveC = .false.
    
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
    call lsyssc_clearMatrix(rdata%rS)
    
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
    
    ! Okay, let us allocate the auxiliary arrays
    allocate(DA(nmaxdofV,nmaxdofV))
    allocate(DB(nmaxdofV,nmaxdofP))
    allocate(DC(nmaxdofP,nmaxdofP))
    allocate(DS(nmaxdofP,nmaxdofP))
    allocate(Ipivot(max(nmaxdofV,nmaxdofP)))
    allocate(IdofsV(nmaxdofV))
    allocate(IdofsP(nmaxdofP))
    
    ! Let us loop over all rows of S
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
      
      ! If the C matrix exists, grab it is main diagonal entry
      if(bHaveC) then
      
        ! Loop over all local pressure DOFs
        do i = 1, ndofP

          ! Get the index of the pressure DOF
          idofp = IdofsP(i)
          
          ! Loop over the row of C
          do k = p_KldC(idofp), p_KldC(idofp+1)-1
          
            ! Get the index of the pressure DOF
            idx = p_KcolC(k)
            
            ! Let us see if this corresponds to one of our local pressure DOFs
            do j = 1, ndofP
              if(idx .eq. IdofsP(j)) then
                DC(i,j) = dsfC*p_DC(k)
                exit
              end if
            end do ! j
          end do ! k
        end do ! i

      end if
      
      ! Let us loop over all velocity DOFs which are adjacent to the current
      ! pressure DOF.
      do i = 1, ndofV
      
        ! Get the index of the velocity DOF
        idofu = IdofsV(i)
        
        ! Let us fetch the local A11/A22 matrices
        do k = p_KldA(idofu), p_KldA(idofu+1)-1
        
          ! Get the column index
          idx = p_KcolA(k)
          
          ! Let us see if this corresponds to one of our local velocity DOFs
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
          ! Let us fetch the local A12/A21 matrices
          do k = p_KldA12(idofu), p_KldA12(idofu+1)-1
          
            ! Get the column index
            idx = p_KcolA12(k)
            
            ! Let us see if this corresponds to one of our local velocity DOFs
            do j = 1, ndofV
              if(idx .eq. IdofsV(j)) then
                DA(      i,ndofV+j) = dsfA12*p_DA12(k)
                DA(ndofV+i,      j) = dsfA21*p_DA21(k)
                exit
              end if
            end do ! j
          end do ! k
        end if

        ! Let us fetch the local B matrices
        do k = p_KldB(idofu), p_KldB(idofu+1)-1
        
          ! Get the column index
          idx = p_KcolB(k)
          
          ! Let us see if this corresponds to one of our local pressure DOFs
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
      call DGESV(2*ndofV,ndofP,DA,nmaxdofV,Ipivot,DB,nmaxdofV,info)
      
      ! Did LAPACK fail? If yes, simply continue with the next pressure DOF.
      if(info .ne. 0) cycle
      
      ! Okay, let us calculate the Schur-complement C := C - D * A^-1 * B
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
          do k = 1, ndofV
            daux1 = daux1 + p_DD1(id1+k-1)*DB(      k,j)
            daux2 = daux2 + p_DD2(id1+k-1)*DB(ndofV+k,j)
          end do ! k
          DC(i,j) = DC(i,j) - dsfD1*daux1 - dsfD2*daux2
        
        end do ! j
      
      end do ! i

      ! Okay, try to solve the local system C*S = id
      DS = 0.0_DP
      do i = 1, ndofP
        DS(i,i) = 1.0_DP
      end do
      call DGESV(ndofP,ndofP,DC,nmaxdofP,Ipivot,DS,nmaxdofP,info)
      
      ! Now if the Schur-complement matrix is regular, we will store the inverse
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
    
    ! That is it
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spsor_solve_NavSt2D(rdata, rsol, rrhs, niterations, domega)

!<description>
  ! INTERNAL ROUTINE:
  ! Performs a specifies number of SP-SOR iterations for 2D Navier-Stokes systems.
!</description>

!<input>
  ! The SP-SOR data structure
  type(t_spsor), target, intent(in) :: rdata
  
  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in) :: rrhs
  
  ! The number of iterations that are to be performed
  integer, intent(in) :: niterations
  
  ! The relaxation parameter in range (0,2)
  real(DP), intent(in) :: domega
!</input>

!<inputoutput>
  ! The iteration vector that is to be updated.
  type(t_vectorblock), intent(inout) :: rsol
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

  integer :: m,n, iter
  logical :: bHaveA12
  real(DP) :: domegaA, domegaS
    
    ! Let us assume we do not have the optional matrices
    bHaveA12 = .false.
    
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
    call lsyssc_getbase_double(rsol%RvectorBlock(1), p_Du1)
    call lsyssc_getbase_double(rsol%RvectorBlock(2), p_Du2)
    call lsyssc_getbase_double(rsol%RvectorBlock(3), p_Dup)
    call lsyssc_getbase_double(rrhs%RvectorBlock(1), p_Df1)
    call lsyssc_getbase_double(rrhs%RvectorBlock(2), p_Df2)
    call lsyssc_getbase_double(rrhs%RvectorBlock(3), p_Dfp)
    if(.not. rdata%bdiagS) call lsyssc_getbase_double(rdata%rw, p_Dw)
        
    ! Fetch the number of velocity and pressure DOFs
    m = rdata%p_rA11%NEQ
    n = rdata%rS%NEQ

    ! Set up relaxation parameters
    ! TODO: Find out which combination is useful!!!
    domegaA = domega
    !domegaA = 1.0_DP
    domegaS = domega
    !domegaS = 1.0_DP
    
    ! perform the desired number of iterations
    do iter = 1, niterations
    
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
            p_KldB,p_KcolB,p_DB1,dsfB1,domegaA)

        ! Solve A22 * u2 = f2 - A21*u1 - B2*p
        call spsor_aux_solveAAB(m,p_Du2,p_Du1,p_Dup,p_Df2,&
            p_KldA,p_KcolA,p_KdiagA,p_DA22,dsfA22,&
            p_KldA12,p_KcolA12,p_DA21,dsfA21,&
            p_KldB,p_KcolB,p_DB2,dsfB2,domegaA)

      else
      
        ! A12/A21 are not present
        call spsor_aux_solveAB_2D(m,p_Du1,p_Du2,p_Dup,p_Df1,p_Df2,&
            p_KldA,p_KcolA,p_KdiagA,p_DA11,p_DA22,dsfA11,dsfA22,&
            p_KldB,p_KcolB,p_DB1,p_DB2,dsfB1,dsfB2,domegaA)
        
      end if
      
      
      ! Step Two
      ! --------
      ! Calculate:
      !
      !  p := p + omega * S * (g - D*u)
      
      if(rdata%bdiagS) then

        ! Calculate p := p - omega * S * (g - D*u)
        call spsor_aux_solveDDS(n,p_Du1,p_Du2,p_Dup,p_Dfp,&
            p_KldD,p_KcolD,p_DD1,p_DD2,dsfD1,dsfD2,p_DS,domegaS)

      else
      
        ! Calculate w := g - D*u
        call spsor_aux_multDD(n,p_Du1,p_Du2,p_Dw,p_Dfp,&
            p_KldD,p_KcolD,p_DD1,p_DD2,dsfD1,dsfD2)

        ! Calculate p := p + omega * S * w
        call lsyssc_scalarMatVec(rdata%rS, rdata%rw, &
            rsol%rvectorBlock(3), domegaS, 1.0_DP)
      
      end if
    
    end do ! iter
    
    ! That is it
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spsor_prec_NavSt2D(rdata, rdef, domega)

!<description>
  ! INTERNAL ROUTINE:
  ! Applies the SP-SOR preconditioner for 2D Navier-Stokes systems.
!</description>

!<input>
  ! The SP-SOR data structure
  type(t_spsor), target, intent(in) :: rdata
 
  ! The relaxation parameter. Default = 1.0
  real(DP), intent(in) :: domega
!</input>

!<inputoutput>
  ! The defect vector that is to be preconditioned.
  type(t_vectorblock), intent(inout) :: rdef
!</inputoutput>

!</subroutine>

  ! Multiplication factors
  real(DP) :: dsfA11, dsfA21, dsfA22, dsfD1, dsfD2
  
  ! access for the matrix arrays
  integer, dimension(:), pointer :: p_KldA,p_KldA12,p_KldD,p_KldS,&
      p_KcolA,p_KcolA12,p_KcolD,p_KcolS,p_KdiagA,p_KdiagS
  real(DP), dimension(:), pointer :: p_DA11,p_DA21,p_DA22,&
      p_DD1,p_DD2,p_DS
  
  ! vector arrays
  real(DP), dimension(:), pointer :: p_Du1, p_Du2, p_Dup, p_Dw

  integer :: m, n
  logical :: bHaveA12
  real(DP) :: domegaA, domegaS
  
    ! Let us assume we do not have the optional matrices
    bHaveA12 = .false.
    
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
      call lsyssc_getbase_double(rdata%p_rA21, p_DA21)
    end if

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
    dsfD1 = rdata%p_rD1%dscaleFactor
    dsfD2 = rdata%p_rD2%dscaleFactor
    
    if(bHaveA12) then
      dsfA21 = rdata%p_rA21%dscaleFactor
      if(dsfA21 .eq. 0.0_DP) bHaveA12 = .false.
    end if
    
    ! Fetch the vector arrays
    call lsyssc_getbase_double(rdef%RvectorBlock(1), p_Du1)
    call lsyssc_getbase_double(rdef%RvectorBlock(2), p_Du2)
    call lsyssc_getbase_double(rdef%RvectorBlock(3), p_Dup)
    if(.not. rdata%bdiagS) call lsyssc_getbase_double(rdata%rw, p_Dw)
        
    ! Fetch the number of velocity and pressure DOFs
    m = rdata%p_rA11%NEQ
    n = rdata%rS%NEQ
    
    ! Set up relaxation parameters
    ! TODO: Find out which combination is useful!!!
    domegaA = domega
    !domegaA = 1.0_DP
    domegaS = domega
    !domegaS = 1.0_DP    
    
    ! Step One
    ! --------
    ! Solve the following equation:
    !
    !                        L * u = u
    !
    ! where L := 1/omega * diag(A) + Ltri(A)
    
    ! Do we have the A12/A21 matrices?
    if(bHaveA12) then
    
      ! Solve for A11 * u1 = u1
      call spsor_aux_precLA(m,p_Du1,p_KldA,p_KcolA,p_KdiagA,&
          p_DA11,dsfA11,domegaA)

      ! Solve for A22 * u2 = u2 - A21 * u1
      call spsor_aux_precLAA(m,p_Du1,p_Du2,p_KldA,p_KcolA,p_KdiagA,&
          p_DA22,dsfA22,p_KldA12,p_KcolA12,p_DA21,dsfA21,domegaA)

    else
    
      ! A12/A21 are not present
      call spsor_aux_precLA_2D(m,p_Du1,p_Du2,p_KldA,p_KcolA,p_KdiagA,&
          p_DA11,p_DA22,dsfA11,dsfA22,domegaA)

    end if
    
    
    ! Step Two
    ! --------
    ! Calculate:
    !
    !  p := omega * S * (p - D*u)
    !
    
    if(rdata%bdiagS) then

      ! Calculate p := omega * S * (p - D*u)
      call spsor_aux_precDDS(n,p_Du1,p_Du2,p_Dup,p_KldD,p_KcolD,&
          p_DD1,p_DD2,dsfD1,dsfD2,p_DS,domegaS)

    else
    
      ! Calculate w := p - D*u
      call spsor_aux_multDD(n,p_Du1,p_Du2,p_Dw,p_Dup,&
          p_KldD,p_KcolD,p_DD1,p_DD2,dsfD1,dsfD2)

      ! Calculate p := omega * S * w
      call lsyssc_scalarMatVec(rdata%rS, rdata%rw, &
          rdef%rvectorBlock(3), domegaS, 0.0_DP)
    
    end if
    
    ! That is it
    
  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine spsor_prec_sym_NavSt2D(rdata, rdef, domega)

!<description>
  ! INTERNAL ROUTINE:
  ! Applies the SP-SSOR preconditioner for 2D Navier-Stokes systems.
!</description>

!<input>
  ! The SP-SOR data structure
  type(t_spsor), target, intent(in) :: rdata
 
  ! The relaxation parameter. Default = 1.0
  real(DP), intent(in) :: domega
!</input>

!<inputoutput>
  ! The defect vector that is to be preconditioned.
  type(t_vectorblock), intent(inout) :: rdef
!</inputoutput>

!</subroutine>

  ! Multiplication factors
  real(DP) :: dsfA11, dsfA12, dsfA21, dsfA22, dsfB1, dsfB2, dsfD1, dsfD2
  
  ! access for the matrix arrays
  integer, dimension(:), pointer :: p_KldA,p_KldA12,p_KldB,p_KldD,p_KldS,&
      p_KcolA,p_KcolA12,p_KcolB,p_KcolD,p_KcolS,p_KdiagA,p_KdiagS
  real(DP), dimension(:), pointer :: p_DA11,p_DA12,p_DA21,p_DA22,&
      p_DB1,p_DB2,p_DD1,p_DD2,p_DS
  
  ! vector arrays
  real(DP), dimension(:), pointer :: p_Du1, p_Du2, p_Dup, p_Dw

  integer :: m, n
  logical :: bHaveA12
  real(DP) :: domegaA, domegaS
  
    ! Let us assume we do not have the optional matrices
    bHaveA12 = .false.
    
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
      if((dsfA12 .eq. 0.0_DP) .or. (dsfA21 .eq. 0.0_DP)) &
        bHaveA12 = .false.
    end if
    
    ! Fetch the vector arrays
    call lsyssc_getbase_double(rdef%RvectorBlock(1), p_Du1)
    call lsyssc_getbase_double(rdef%RvectorBlock(2), p_Du2)
    call lsyssc_getbase_double(rdef%RvectorBlock(3), p_Dup)
    if(.not. rdata%bdiagS) call lsyssc_getbase_double(rdata%rw, p_Dw)
        
    ! Fetch the number of velocity and pressure DOFs
    m = rdata%p_rA11%NEQ
    n = rdata%rS%NEQ
    
    ! Set up relaxation parameters
    ! TODO: Find out which combination is useful!!!
    domegaA = domega
    !domegaA = 1.0_DP
    domegaS = domega
    !domegaS = 1.0_DP
    
    ! Step One
    ! --------
    ! Solve the following equation:
    !
    !                        L * u = u
    !
    ! where L := 1/omega * diag(A) + Ltri(A)
    
    ! Do we have the A12/A21 matrices?
    if(bHaveA12) then
    
      ! Solve for A11 * u1 = u1
      call spsor_aux_precLA(m,p_Du1,p_KldA,p_KcolA,p_KdiagA,&
          p_DA11,dsfA11,domegaA)

      ! Solve for A22 * u2 = u2 - A21 * u1
      call spsor_aux_precLAA(m,p_Du1,p_Du2,p_KldA,p_KcolA,p_KdiagA,&
          p_DA22,dsfA22,p_KldA12,p_KcolA12,p_DA21,dsfA21,domegaA)

    else
    
      ! A12/A21 are not present
      call spsor_aux_precLA_2D(m,p_Du1,p_Du2,p_KldA,p_KcolA,p_KdiagA,&
          p_DA11,p_DA22,dsfA11,dsfA22,domegaA)

    end if
    
    
    ! Step Two
    ! --------
    ! Calculate:
    !
    !  p := omega * S * (p - D*u)
    !
    
    ! Is S a diagonal matrix?
    if(rdata%bdiagS) then

      ! Calculate p := omega * S * (p - D*u)
      call spsor_aux_precDDS(n,p_Du1,p_Du2,p_Dup,p_KldD,p_KcolD,&
          p_DD1,p_DD2,dsfD1,dsfD2,p_DS,domegaS)

    else
    
      ! Calculate w := p - D*u
      call spsor_aux_multDD(n,p_Du1,p_Du2,p_Dw,p_Dup,&
          p_KldD,p_KcolD,p_DD1,p_DD2,dsfD1,dsfD2)

      ! Calculate p := omega * S * w
      call lsyssc_scalarMatVec(rdata%rS, rdata%rw, &
          rdef%rvectorBlock(3), domegaS, 0.0_DP)
    
    end if
    
    ! Step Three
    ! ----------
    ! Solve the following equation:
    !
    !                        U * u = u - B*p
    !
    ! where L := id + ( 1/omega * diag(A) )^{-1} * Utri(A)
    
    ! Do we have the A12/A21 matrices?
    if(bHaveA12) then
    
      ! Solve for A22 * u2 = u2 - B2*p
      call spsor_aux_precUAB(m,p_Du2,p_Dup,p_KldA,p_KcolA,p_KdiagA,&
          p_DA22,dsfA22,p_KldB,p_KcolB,p_DB2,dsfB2,domegaA)

      ! Solve for A11 * u1 = u2 - A12 * u2 - B1*p
      call spsor_aux_precUAAB(m,p_Du1,p_Du2,p_Dup,p_KldA,p_KcolA,p_KdiagA,&
          p_DA11,dsfA11,p_KldA12,p_KcolA12,p_DA12,dsfA12,p_KldB,p_KcolB,&
          p_DB1,dsfB1,domegaA)

    else
    
      ! A12/A21 are not present
      call spsor_aux_precUAB_2D(m,p_Du1,p_Du2,p_Dup,p_KldA,p_KcolA,p_KdiagA,&
          p_DA11,p_DA22,dsfA11,dsfA22,p_KldB,p_KcolB,p_DB1,p_DB2,&
          dsfB1,dsfB2,domegaA)

    end if

    ! That is it
    
  end subroutine


  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !
  !       A U X I L I A R Y   R O U T I N E S   F O R   S O L V I N G
  !
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


!<subroutine>

  pure subroutine spsor_aux_solveAB_2D(n,Du1,Du2,Dup,Df1,Df2,KldA,KcolA,KdiagA,&
                                       DA1,DA2,dsfA1,dsfA2,KldB,KcolB,&
                                       DB1,DB2,dsfB1,dsfB2,domega)

!<description>
  ! INTERNAL AUXILIARY ROUTINE:
  ! Performs an SOR iteration on the following equation:
  !
  ! <verb>
  !              ( A1  0  ) * ( u1 ) = ( f1 ) - ( B1 ) * ( p )
  !              ( 0   A2 )   ( u2 )   ( f2 )   ( B2 )
  ! </verb>
  !
!</description>

!<inputoutput>
  real(DP), dimension(*), intent(inout) :: Du1, Du2
!</inputoutput>

!<input>
  integer, intent(in) :: n
  
  real(DP), dimension(*), intent(in) :: Dup, Df1, Df2
  
  integer, dimension(*), intent(in) :: KldA, KcolA, KdiagA
  real(DP), dimension(*), intent(in) :: DA1, DA2
  real(DP), intent(in) :: dsfA1, dsfA2

  integer, dimension(*), intent(in) :: KldB, KcolB
  real(DP), dimension(*), intent(in) :: DB1, DB2
  real(DP), intent(in) :: dsfB1, dsfB2
  
  real(DP), intent(in) :: domega
!</input>

!</subroutine>

  integer :: i,j,k
  real(DP) :: dt,df,dg,daux1,daux2
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
      
      ! Subtract A*u from rhs
      daux1 = 0.0_DP
      daux2 = 0.0_DP
      do j = KldA(i), KldA(i+1)-1
        k = KcolA(j)
        daux1 = daux1 + DA1(j)*Du1(k)
        daux2 = daux2 + DA2(j)*Du2(k)
      end do
      df = domega*(df - daux1)
      dg = domega*(dg - daux2)

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
      
      ! Solve
      k = KdiagA(i)
      Du1(i) = Du1(i) + df / DA1(k)
      Du2(i) = Du2(i) + dg / DA2(k)
    
    end do ! i
    
  end subroutine ! spsor_aux_solveAB_2D

  ! ***************************************************************************

!<subroutine>

  pure subroutine spsor_aux_solveAB_3D(n,Du1,Du2,Du3,Dup,Df1,Df2,Df3,&
                                       KldA,KcolA,KdiagA,DA1,DA2,DA3,&
                                       dsfA1,dsfA2,dsfA3,KldB,KcolB,&
                                       DB1,DB2,DB3,dsfB1,dsfB2,dsfB3,domega)

!<description>
  ! INTERNAL AUXILIARY ROUTINE:
  ! Performs an SOR iteration on the following equation:
  !
  ! <verb>
  !              ( A1  0  0  ) * ( u1 ) = ( f1 ) - ( B1 ) * ( p )
  !              ( 0   A2 0  )   ( u2 )   ( f2 )   ( B2 )
  !              ( 0   0  A3 )   ( u3 )   ( f3 )   ( B3 )
  ! </verb>
  !
!</description>

!<inputoutput>
  real(DP), dimension(*), intent(inout) :: Du1, Du2, Du3
!</inputoutput>

!<input>
  integer, intent(in) :: n
  
  real(DP), dimension(*), intent(in) :: Dup, Df1, Df2, Df3
  
  integer, dimension(*), intent(in) :: KldA, KcolA, KdiagA
  real(DP), dimension(*), intent(in) :: DA1, DA2, DA3
  real(DP), intent(in) :: dsfA1, dsfA2, dsfA3

  integer, dimension(*), intent(in) :: KldB, KcolB
  real(DP), dimension(*), intent(in) :: DB1, DB2, DB3
  real(DP), intent(in) :: dsfB1, dsfB2, dsfB3
  
  real(DP), intent(in) :: domega
!</input>

!</subroutine>

  integer :: i,j,k
  real(DP) :: dt,dfu1,dfu2,dfu3,daux1,daux2,daux3
  real(DP) :: dsf1,dsf2,dsf3,dsb1,dsb2,dsb3
  
    ! Pre-calculate scaling factors
    dsf1 = 1.0_DP / dsfA1
    dsf2 = 1.0_DP / dsfA2
    dsf3 = 1.0_DP / dsfA3
    dsb1 = dsfB1 / dsfA1
    dsb2 = dsfB2 / dsfA2
    dsb3 = dsfB3 / dsfA3
  
    ! Loop over all rows
    do i = 1, n
    
      ! Fetch the rhs entries
      dfu1 = Df1(i)*dsf1
      dfu2 = Df2(i)*dsf2
      dfu3 = Df3(i)*dsf3
      
      ! Subtract A*u from rhs
      daux1 = 0.0_DP
      daux2 = 0.0_DP
      daux2 = 0.0_DP
      do j = KldA(i), KldA(i+1)-1
        k = KcolA(j)
        daux1 = daux1 + DA1(j)*Du1(k)
        daux2 = daux2 + DA2(j)*Du2(k)
        daux3 = daux3 + DA3(j)*Du3(k)
      end do
      dfu1 = domega*(dfu1 - daux1)
      dfu2 = domega*(dfu2 - daux2)
      dfu3 = domega*(dfu3 - daux3)

      ! Subtract B*p from rhs
      daux1 = 0.0_DP
      daux2 = 0.0_DP
      daux3 = 0.0_DP
      do j = KldB(i), KldB(i+1)-1
        dt = Dup(KcolB(j))
        daux1 = daux1 + DB1(j)*dt
        daux2 = daux2 + DB2(j)*dt
        daux3 = daux3 + DB3(j)*dt
      end do
      dfu1 = dfu1 - dsb1*daux1
      dfu2 = dfu2 - dsb2*daux2
      dfu3 = dfu3 - dsb3*daux3
      
      ! Solve
      k = KdiagA(i)
      Du1(i) = Du1(i) + dfu1 / DA1(k)
      Du2(i) = Du2(i) + dfu2 / DA2(k)
      Du3(i) = Du3(i) + dfu3 / DA3(k)
    
    end do ! i
    
  end subroutine ! spsor_aux_solveAB_3D

  ! ***************************************************************************

!<subroutine>

  pure subroutine spsor_aux_solveAAB(n,Du1,Du2,Dup,Df,KldA,KcolA,KdiagA,&
                                     DA,dsfA,KldA2,KcolA2,DA2,dsfA2,&
                                     KldB,KcolB,DB,dsfB,domega)

!<description>
  ! INTERNAL AUXILIARY ROUTINE:
  ! Performs an SOR iteration on the following equation:
  !
  ! <verb>
  !                A * u1 = f - A2 * u2 - B*p
  ! </verb>
  !
!</description>

!<inputoutput>
  real(DP), dimension(*), intent(inout) :: Du1
!</inputoutput>

!<input>
  integer, intent(in) :: n
  
  real(DP), dimension(*), intent(in) :: Du2, Dup, Df
  
  integer, dimension(*), intent(in) :: KldA, KcolA, KdiagA
  real(DP), dimension(*), intent(in) :: DA
  real(DP), intent(in) :: dsfA
  
  integer, dimension(*), intent(in) :: KldA2, KcolA2
  real(DP), dimension(*), intent(in) :: DA2
  real(DP), intent(in) :: dsfA2

  integer, dimension(*), intent(in) :: KldB, KcolB
  real(DP), dimension(*), intent(in) :: DB
  real(DP), intent(in) :: dsfB
  
  real(DP), intent(in) :: domega
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
      
      ! Subtract A*u1 from rhs
      daux = 0.0_DP
      do j = KldA(i), KldA(i+1)-1
        daux = daux + DA(j)*Du1(KcolA(j))
      end do
      df1 = domega*(df1 - daux)

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
      
      ! Solve
      Du1(i) = Du1(i) + df1 / DA(KdiagA(i))
    
    end do ! i
    
  end subroutine ! spsor_aux_solveAAB

  ! ***************************************************************************

!<subroutine>

  pure subroutine spsor_aux_solveAAAB(n,Du1,Du2,Du3,Dup,Df,KldA,KcolA,KdiagA,&
                                     DA,dsfA,KldA2,KcolA2,DA2,DA3,dsfA2,dsfA3,&
                                     KldB,KcolB,DB,dsfB,domega)

!<description>
  ! INTERNAL AUXILIARY ROUTINE:
  ! Performs an SOR iteration on the following equation:
  !
  ! <verb>
  !                A * u1 = f - A2 * u2 - A3 * u3 - B*p
  ! </verb>
  !
  ! where A2 and A3 share the same structure
!</description>

!<inputoutput>
  real(DP), dimension(*), intent(inout) :: Du1
!</inputoutput>

!<input>
  integer, intent(in) :: n
  
  real(DP), dimension(*), intent(in) :: Du2, Du3, Dup, Df
  
  integer, dimension(*), intent(in) :: KldA, KcolA, KdiagA
  real(DP), dimension(*), intent(in) :: DA
  real(DP), intent(in) :: dsfA
  
  integer, dimension(*), intent(in) :: KldA2, KcolA2
  real(DP), dimension(*), intent(in) :: DA2, DA3
  real(DP), intent(in) :: dsfA2, dsfA3

  integer, dimension(*), intent(in) :: KldB, KcolB
  real(DP), dimension(*), intent(in) :: DB
  real(DP), intent(in) :: dsfB
  
  real(DP), intent(in) :: domega
!</input>

!</subroutine>

  integer :: i,j,k
  real(DP) :: df1,daux,daux2
  real(DP) :: dsa2,dsa3,dsf,dsb
  
    ! Pre-calculate scaling factors
    dsf = 1.0_DP / dsfA
    dsa2 = dsfA2 / dsfA
    dsa3 = dsfA3 / dsfA
    dsb = dsfB / dsfA
  
    ! Loop over all rows
    do i = 1, n
    
      ! Fetch the rhs entries
      df1 = Df(i)*dsf
      
      ! Subtract A*u1 from rhs
      daux = 0.0_DP
      do j = KldA(i), KldA(i+1)-1
        daux = daux + DA(j)*Du1(KcolA(j))
      end do
      df1 = domega*(df1 - daux)

      ! Subtract B*p from rhs
      daux = 0.0_DP
      do j = KldB(i), KldB(i+1)-1
        daux = daux + DB(j)*Dup(KcolB(j))
      end do
      df1 = df1 - dsb*daux
      
      ! Subtract A2*u2 and A3*u3 from rhs
      daux = 0.0_DP
      daux2 = 0.0_DP
      do j = KldA2(i), KldA2(i+1)-1
        k = KcolA2(j)
        daux  = daux  + DA2(j)*Du2(k)
        daux2 = daux2 + DA3(j)*Du3(k)
      end do
      df1 = df1 - dsa2*daux - dsa3*daux2
      
      ! Solve
      Du1(i) = Du1(i) + df1 / DA(KdiagA(i))
    
    end do ! i
    
  end subroutine ! spsor_aux_solveAAAB

  ! ***************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine spsor_aux_solveDDS(n,Du1,Du2,Dup,Df,KldD,KcolD,DD1,DD2,&
                                dsfD1,dsfD2,DS,domega)

!<description>
  ! INTERNAL AUXILIARY ROUTINE:
  ! Calculates:
  !
  ! <verb>
  !                p := p + relax * S * (f - D1 * u1 - D2 * u2)
  ! </verb>
  !
  ! where S is a diagonal matrix, and D1 and D2 share the same structure
!</description>

!<inputoutput>
  real(DP), dimension(*), intent(inout) :: Dup
!</inputoutput>

!<input>
  integer, intent(in) :: n
  
  real(DP), dimension(*), intent(in) :: Du1, Du2, Df
  
  integer, dimension(*), intent(in) :: KldD, KcolD
  real(DP), dimension(*), intent(in) :: DD1, DD2
  real(DP), intent(in) :: dsfD1, dsfD2

  real(DP), dimension(*), intent(in) :: DS
  
  real(DP), intent(in) :: domega
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
      Dup(i) = Dup(i) + domega*DS(i)*(Df(i) - dsfD1*daux1 - dsfD2*daux2)
    
    end do ! i
    !$omp end parallel do
    
  end subroutine ! spsor_aux_solveDDS

  ! ***************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine spsor_aux_solveDDDS(n,Du1,Du2,Du3,Dup,Df,KldD,KcolD,&
                                 DD1,DD2,DD3,dsfD1,dsfD2,dsfD3,DS,domega)

!<description>
  ! INTERNAL AUXILIARY ROUTINE:
  ! Calculates:
  !
  ! <verb>
  !                p := p + relax * S * (f - D1 * u1 - D2 * u2 - D3 * u3)
  ! </verb>
  !
  ! where S is a diagonal matrix, and D1, D2 and D3 share the same structure
!</description>

!<inputoutput>
  real(DP), dimension(*), intent(inout) :: Dup
!</inputoutput>

!<input>
  integer, intent(in) :: n
  
  real(DP), dimension(*), intent(in) :: Du1, Du2, Du3, Df
  
  integer, dimension(*), intent(in) :: KldD, KcolD
  real(DP), dimension(*), intent(in) :: DD1, DD2, DD3
  real(DP), intent(in) :: dsfD1, dsfD2, dsfD3

  real(DP), dimension(*), intent(in) :: DS
  
  real(DP), intent(in) :: domega
!</input>

!</subroutine>

  integer :: i,j,k
  real(DP) :: daux1, daux2, daux3
  
    ! Loop over all rows
    !$omp parallel do private(j,k,daux1,daux2,daux3) if(n .gt. 1000)
    do i = 1, n
    
      ! Calculate D1*u and D2*v
      daux1 = 0.0_DP
      daux2 = 0.0_DP
      daux3 = 0.0_DP
      do j = KldD(i), KldD(i+1)-1
        k = KcolD(j)
        daux1 = daux1 + DD1(j)*Du1(k)
        daux2 = daux2 + DD2(j)*Du2(k)
        daux3 = daux3 + DD3(j)*Du3(k)
      end do
      
      ! Solve
      Dup(i) = Dup(i) + domega*DS(i)*(Df(i) - dsfD1*daux1 - dsfD2*daux2 - dsfD3*daux3)
    
    end do ! i
    !$omp end parallel do
    
  end subroutine ! spsor_aux_solveDDDS

  ! ***************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine spsor_aux_multDD(n,Du1,Du2,Dw,Df,KldD,KcolD,DD1,DD2,&
                              dsfD1,dsfD2)

!<description>
  ! INTERNAL AUXILIARY ROUTINE:
  ! Calculates:
  !
  ! <verb>
  !                w := f - D1 * u1 - D2 * u2
  ! </verb>
  !
  ! where D1 and D2 share the same structure
!</description>

!<inputoutput>
  real(DP), dimension(*), intent(inout) :: Dw
!</inputoutput>

!<input>
  integer, intent(in) :: n
  
  real(DP), dimension(*), intent(in) :: Du1, Du2, Df
  
  integer, dimension(*), intent(in) :: KldD, KcolD
  real(DP), dimension(*), intent(in) :: DD1, DD2
  real(DP), intent(in) :: dsfD1, dsfD2
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

  ! ***************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine spsor_aux_multDDD(n,Du1,Du2,Du3,Dw,Df,KldD,KcolD,&
                               DD1,DD2,DD3,dsfD1,dsfD2,dsfD3)

!<description>
  ! INTERNAL AUXILIARY ROUTINE:
  ! Calculates:
  !
  ! <verb>
  !                w := f - D1 * u1 - D2 * u2 - D3 * u3
  ! </verb>
  !
  ! where D1, D2 and D3 share the same structure
!</description>

!<inputoutput>
  real(DP), dimension(*), intent(inout) :: Dw
!</inputoutput>

!<input>
  integer, intent(in) :: n
  
  real(DP), dimension(*), intent(in) :: Du1, Du2, Du3, Df
  
  integer, dimension(*), intent(in) :: KldD, KcolD
  real(DP), dimension(*), intent(in) :: DD1, DD2, DD3
  real(DP), intent(in) :: dsfD1, dsfD2, dsfD3
!</input>

!</subroutine>

  integer :: i,j,k
  real(DP) :: daux1, daux2, daux3
  
    ! Loop over all rows
    !$omp parallel do private(j,k,daux1,daux2,daux3) if(n .gt. 1000)
    do i = 1, n
    
      daux1 = 0.0_DP
      daux2 = 0.0_DP
      daux3 = 0.0_DP
      do j = KldD(i), KldD(i+1)-1
        k = KcolD(j)
        daux1 = daux1 + DD1(j)*Du1(k)
        daux2 = daux2 + DD2(j)*Du2(k)
        daux3 = daux3 + DD3(j)*Du3(k)
      end do
      Dw(i) = Df(i) - dsfD1*daux1 - dsfD2*daux2 - dsfD3*daux3
    
    end do ! i
    !$omp end parallel do
    
  end subroutine ! spsor_aux_multDDD


  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  !
  ! A U X I L I A R Y   R O U T I N E S   F O R   P R E C O N D I T I O N I N G
  !
  ! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


!<subroutine>

  pure subroutine spsor_aux_precLA(n,Du,KldA,KcolA,KdiagA,DA,dsfA,domega)

!<description>
  ! INTERNAL AUXILIARY ROUTINE:
  !
  ! Solves the following system:
  !
  ! <verb>
  !                L * u = f
  ! </verb>
  !
  ! where L := 1/omega * diag(A) + ltri(A)
  !
!</description>

!<inputoutput>
  ! On entry, the right-hand-side f of the system
  ! On exit, the solution u of the system
  real(DP), dimension(*), intent(inout) :: Du
!</inputoutput>

!<input>
  ! The dimension of A
  integer, intent(in) :: n
  
  ! The matrix structure of A
  integer, dimension(*), intent(in) :: KldA, KcolA, KdiagA
  
  ! The matrix data of A
  real(DP), dimension(*), intent(in) :: DA

  ! The scaling factor of A
  real(DP), intent(in) :: dsfA
  
  ! The relaxation parameter in range (0,2)
  real(DP), intent(in) :: domega
!</input>

!</subroutine>

  integer :: i,j,k
  real(DP) :: daux, dsa
  
    ! Pre-calculate scaling factor
    dsa = 1.0_DP / dsfA
  
    ! Loop over all rows
    do i = 1, n
    
      k = KdiagA(i)
    
      daux = 0.0_DP
      do j = KldA(i), k-1
        daux = daux + DA(j)*Du(KcolA(j))
      end do
      
      Du(i) = domega*(dsa*Du(i) - daux) / DA(k)
    
    end do ! i
    
  end subroutine ! spsor_aux_precLA

  ! ***************************************************************************

!<subroutine>

  pure subroutine spsor_aux_precLA_2D(n,Du1,Du2,KldA,KcolA,KdiagA,DA1,DA2,&
                                      dsfA1,dsfA2,domega)

!<description>
  ! INTERNAL AUXILIARY ROUTINE:
  !
  ! Solves the following system:
  !
  ! <verb>
  !              ( L1  0  ) * ( u1 ) = ( f1 )
  !              ( 0   L2 )   ( u2 )   ( f2 )
  ! </verb>
  !
  ! where Li := 1/omega * diag(Ai) + ltri(Ai)
  !
!</description>

!<inputoutput>
  ! On entry, the right-hand-side (f1,f2) of the system
  ! On exit, the solution (u1,u2) of the system
  real(DP), dimension(*), intent(inout) :: Du1, Du2
!</inputoutput>

!<input>
  ! The dimension of Ai
  integer, intent(in) :: n
  
  ! The matrix structure of Ai
  integer, dimension(*), intent(in) :: KldA, KcolA, KdiagA
  
  ! The matrix data of Ai
  real(DP), dimension(*), intent(in) :: DA1, DA2

  ! The scaling factors of Ai
  real(DP), intent(in) :: dsfA1, dsfA2
  
  ! The relaxation parameter in range (0,2)
  real(DP), intent(in) :: domega
!</input>

!</subroutine>

  integer :: i,j,k,m
  real(DP) :: daux1,daux2,dsa1,dsa2
  
    ! pre-calculate scaling factors
    dsa1 = 1.0_DP / dsfA1
    dsa2 = 1.0_DP / dsfA2
  
    ! Loop over all rows
    do i = 1, n
    
      k = KdiagA(i)
    
      daux1 = 0.0_DP
      daux2 = 0.0_DP
      do j = KldA(i), k-1
        m = KcolA(j)
        daux1 = daux1 + DA1(j)*Du1(m)
        daux2 = daux2 + DA2(j)*Du2(m)
      end do
      
      Du1(i) = domega*(dsa1*Du1(i) - daux1) / DA1(k)
      Du2(i) = domega*(dsa2*Du2(i) - daux2) / DA2(k)
    
    end do ! i
    
  end subroutine ! spsor_aux_precLA_2D

  ! ***************************************************************************

!<subroutine>

  pure subroutine spsor_aux_precLA_3D(n,Du1,Du2,Du3,KldA,KcolA,KdiagA,&
                                      DA1,DA2,DA3,dsfA1,dsfA2,dsfA3,domega)

!<description>
  ! INTERNAL AUXILIARY ROUTINE:
  !
  ! Solves the following system:
  !
  ! <verb>
  !              ( L1  0  0  ) * ( u1 ) = ( f1 )
  !              ( 0   L2 0  )   ( u2 )   ( f2 )
  !              ( 0   0  L3 )   ( u3 )   ( f3 )
  ! </verb>
  !
  ! where Li := 1/omega * diag(Ai) + ltri(Ai)
  !
!</description>

!<inputoutput>
  ! On entry, the right-hand-side (f1,f2) of the system
  ! On exit, the solution (u1,u2) of the system
  real(DP), dimension(*), intent(inout) :: Du1, Du2, Du3
!</inputoutput>

!<input>
  ! The dimension of Ai
  integer, intent(in) :: n
  
  ! The matrix structure of Ai
  integer, dimension(*), intent(in) :: KldA, KcolA, KdiagA
  
  ! The matrix data of Ai
  real(DP), dimension(*), intent(in) :: DA1, DA2, DA3

  ! The scaling factors of Ai
  real(DP), intent(in) :: dsfA1, dsfA2, dsfA3
  
  ! The relaxation parameter in range (0,2)
  real(DP), intent(in) :: domega
!</input>

!</subroutine>

  integer :: i,j,k,m
  real(DP) :: daux1,daux2,daux3,dsa1,dsa2,dsa3
  
    ! pre-calculate scaling factors
    dsa1 = 1.0_DP / dsfA1
    dsa2 = 1.0_DP / dsfA2
    dsa3 = 1.0_DP / dsfA3
  
    ! Loop over all rows
    do i = 1, n
    
      k = KdiagA(i)
    
      daux1 = 0.0_DP
      daux2 = 0.0_DP
      daux3 = 0.0_DP
      do j = KldA(i), k-1
        m = KcolA(j)
        daux1 = daux1 + DA1(j)*Du1(m)
        daux2 = daux2 + DA2(j)*Du2(m)
        daux3 = daux3 + DA3(j)*Du3(m)
      end do
      
      Du1(i) = domega*(dsa1*Du1(i) - daux1) / DA1(k)
      Du2(i) = domega*(dsa2*Du2(i) - daux2) / DA2(k)
      Du3(i) = domega*(dsa3*Du3(i) - daux3) / DA3(k)
    
    end do ! i
    
  end subroutine ! spsor_aux_precLA_3D

  ! ***************************************************************************

!<subroutine>

  pure subroutine spsor_aux_precLAA(n,Du,Du2,KldA,KcolA,KdiagA,DA,dsfA,&
                                    KldA2,KcolA2,DA2,dsfA2,domega)

!<description>
  ! INTERNAL AUXILIARY ROUTINE:
  !
  ! Solves the following system:
  !
  ! <verb>
  !              L * u = f - A2 * u2
  ! </verb>
  !
  ! where L := 1/omega * diag(A) + ltri(A)
  !
!</description>

!<inputoutput>
  ! On entry, the right-hand-side f of the system
  ! On exit, the solution u of the system
  real(DP), dimension(*), intent(inout) :: Du
!</inputoutput>

!<input>
  ! The dimension of Ai
  integer, intent(in) :: n

  ! The input vector u2
  real(DP), dimension(*), intent(in) :: Du2
  
  ! The matrix structure of A
  integer, dimension(*), intent(in) :: KldA, KcolA, KdiagA
  
  ! The matrix data of A
  real(DP), dimension(*), intent(in) :: DA

  ! The scaling factors of A
  real(DP), intent(in) :: dsfA
  
  ! The matrix structure of A2
  integer, dimension(*), intent(in) :: KldA2, KcolA2
  
  ! The matrix data of A2
  real(DP), dimension(*), intent(in) :: DA2

  ! The scaling factors of A2
  real(DP), intent(in) :: dsfA2
  
  ! The relaxation parameter in range (0,2)
  real(DP), intent(in) :: domega
!</input>

!</subroutine>

  integer :: i,j,k
  real(DP) :: dom,daux,daux2
  
    ! pre-calculate weight factor
    dom = domega / dsfA
  
    ! Loop over all rows
    do i = 1, n
    
      daux2 = 0.0_DP
      do j = KldA2(i), KldA2(i+1)-1
        daux2 = daux2 + DA2(j)*Du2(KcolA2(j))
      end do

      k = KdiagA(i)
    
      daux = 0.0_DP
      do j = KldA(i), k-1
        daux = daux + DA(j)*Du(KcolA(j))
      end do
      
      Du(i) = dom*(Du(i) - dsfA*daux - dsfA2*daux2) / DA(k)
    
    end do ! i
    
  end subroutine ! spsor_aux_precLAA

  ! ***************************************************************************

!<subroutine>

  pure subroutine spsor_aux_precLAAA(n,Du,Du2,Du3,KldA,KcolA,KdiagA,DA,dsfA,&
                                     KldA2,KcolA2,DA2,DA3,dsfA2,dsfA3,domega)

!<description>
  ! INTERNAL AUXILIARY ROUTINE:
  !
  ! Solves the following system:
  !
  ! <verb>
  !              L * u = f - A2 * u2 - A3 * u3
  ! </verb>
  !
  ! where L := 1/omega * diag(A) + ltri(A)
  !
!</description>

!<inputoutput>
  ! On entry, the right-hand-side f of the system
  ! On exit, the solution u of the system
  real(DP), dimension(*), intent(inout) :: Du
!</inputoutput>

!<input>
  ! The dimension of Ai
  integer, intent(in) :: n

  ! The input vectors u2, u3
  real(DP), dimension(*), intent(in) :: Du2, Du3
  
  ! The matrix structure of A
  integer, dimension(*), intent(in) :: KldA, KcolA, KdiagA
  
  ! The matrix data of A
  real(DP), dimension(*), intent(in) :: DA

  ! The scaling factor of A
  real(DP), intent(in) :: dsfA

  ! The matrix structure of A2, A3
  integer, dimension(*), intent(in) :: KldA2, KcolA2
  
  ! The matrix data of A2, A3
  real(DP), dimension(*), intent(in) :: DA2, DA3

  ! The scaling factors of A2, A2
  real(DP), intent(in) :: dsfA2, dsfA3
  
  ! The relaxation parameter in range (0,2)
  real(DP), intent(in) :: domega
!</input>

!</subroutine>

  integer :: i,j,k,m
  real(DP) :: dom,daux,daux2,daux3
  
    ! pre-calculate weight factor
    dom = domega / dsfA
  
    ! Loop over all rows
    do i = 1, n
    
      daux2 = 0.0_DP
      daux3 = 0.0_DP
      do j = KldA2(i), KldA2(i+1)-1
        m = KcolA2(j)
        daux2 = daux2 + DA2(j)*Du2(m)
        daux3 = daux3 + DA3(j)*Du3(m)
      end do

      k = KdiagA(i)
    
      daux = 0.0_DP
      do j = KldA(i), k-1
        daux = daux + DA(j)*Du(KcolA(j))
      end do
      
      Du(i) = dom*(Du(i) - dsfA*daux - dsfA2*daux2 - dsfA3*daux3) / DA(k)
    
    end do ! i
    
  end subroutine ! spsor_aux_precLAAA


  ! ***************************************************************************

!<subroutine>

  pure subroutine spsor_aux_precUAB(n,Du,Dup,KldA,KcolA,KdiagA,DA,dsfA,&
                                    KldB,KcolB,DB,dsfB,domega)

!<description>
  ! INTERNAL AUXILIARY ROUTINE:
  !
  ! Solves the following system:
  !
  ! <verb>
  !                U * u = f - omega*B*p
  ! </verb>
  !
  ! where <tex>$ U := id + ( 1/omega * diag(A) )^{-1} * ltri(A) $</tex>
  !
!</description>

!<inputoutput>
  ! On entry, the right-hand-side f of the system
  ! On exit, the solution u of the system
  real(DP), dimension(*), intent(inout) :: Du
!</inputoutput>

!<input>
  ! The dimension of A
  integer, intent(in) :: n

  ! The input vector p
  real(DP), dimension(*), intent(inout) :: Dup
  
  ! The matrix structure of A
  integer, dimension(*), intent(in) :: KldA, KcolA, KdiagA
  
  ! The matrix data of A
  real(DP), dimension(*), intent(in) :: DA

  ! The scaling factor of A
  real(DP), intent(in) :: dsfA
  
  ! The matrix structure of B
  integer, dimension(*), intent(in) :: KldB, KcolB
  
  ! The matrix data of B
  real(DP), dimension(*), intent(in) :: DB

  ! The scaling factor of B
  real(DP), intent(in) :: dsfB
  
  ! The relaxation parameter in range (0,2)
  real(DP), intent(in) :: domega
!</input>

!</subroutine>

  integer :: i,j,k
  real(DP) :: daux,dauxb,dsb
  
    ! Pre-calculate scaling factor
    dsb = (domega*dsfB) / dsfA
  
    ! Loop over all rows
    do i = n, 1, -1
    
      dauxb = 0.0_DP
      do j = KldB(i), KldB(i+1)-1
        dauxb = dauxb + DB(j)*Dup(KcolB(j))
      end do
    
      k = KdiagA(i)
    
      daux = 0.0_DP
      do j = KldA(i+1)-1, k+1, -1
        daux = daux + DA(j)*Du(KcolA(j))
      end do
      
      Du(i) = Du(i) - (domega*daux + dsb*dauxb) / DA(k)
    
    end do ! i
    
  end subroutine ! spsor_aux_precUAB

  ! ***************************************************************************

!<subroutine>

  pure subroutine spsor_aux_precUAB_2D(n,Du1,Du2,Dup,KldA,KcolA,KdiagA,&
                                       DA1,DA2,dsfA1,dsfA2,KldB,KcolB,&
                                       DB1,DB2,dsfB1,dsfB2,domega)

!<description>
  ! INTERNAL AUXILIARY ROUTINE:
  !
  ! Solves the following system:
  !
  ! <verb>
  !               ( U1  0  ) * ( u1 ) = ( f1 ) - omega * ( B1 ) * ( p )
  !               ( 0   U2 )   ( u2 )   ( f2 )           ( B2 )
  ! </verb>
  !
  ! where <tex>$ Ui := id + ( 1/omega * diag(Ai) )^{-1} * ltri(Ai) $</tex>
  !
!</description>

!<inputoutput>
  ! On entry, the right-hand-side (f1,f2) of the system
  ! On exit, the solution (u1,u2) of the system
  real(DP), dimension(*), intent(inout) :: Du1, Du2
!</inputoutput>

!<input>
  ! The dimension of A
  integer, intent(in) :: n

  ! The input vector p
  real(DP), dimension(*), intent(inout) :: Dup
  
  ! The matrix structure of Ai
  integer, dimension(*), intent(in) :: KldA, KcolA, KdiagA
  
  ! The matrix data of Ai
  real(DP), dimension(*), intent(in) :: DA1, DA2

  ! The scaling factor of Ai
  real(DP), intent(in) :: dsfA1, dsfA2
  
  ! The matrix structure of Bi
  integer, dimension(*), intent(in) :: KldB, KcolB
  
  ! The matrix data of Bi
  real(DP), dimension(*), intent(in) :: DB1, DB2

  ! The scaling factor of Bi
  real(DP), intent(in) :: dsfB1, dsfB2
  
  ! The relaxation parameter in range (0,2)
  real(DP), intent(in) :: domega
!</input>

!</subroutine>

  integer :: i,j,k,m
  real(DP) :: dt,daux1,daux2,dauxb1,dauxb2,dsb1,dsb2
  
    ! Pre-calculate scaling factor
    dsb1 = (domega*dsfB1) / dsfA1
    dsb2 = (domega*dsfB2) / dsfA2
  
    ! Loop over all rows
    do i = n, 1, -1
    
      dauxb1 = 0.0_DP
      dauxb2 = 0.0_DP
      do j = KldB(i), KldB(i+1)-1
        dt = Dup(KcolB(j))
        dauxb1 = dauxb1 + DB1(j)*dt
        dauxb2 = dauxb2 + DB2(j)*dt
      end do
    
      k = KdiagA(i)
    
      daux1 = 0.0_DP
      daux2 = 0.0_DP
      do j = KldA(i+1)-1, k+1, -1
        m = KcolA(j)
        daux1 = daux1 + DA1(j)*Du1(m)
        daux2 = daux2 + DA2(j)*Du2(m)
      end do
      
      Du1(i) = Du1(i) - (domega*daux1 + dsb1*dauxb1) / DA1(k)
      Du2(i) = Du2(i) - (domega*daux2 + dsb2*dauxb2) / DA2(k)
    
    end do ! i
    
  end subroutine ! spsor_aux_precUAB_2D


  ! ***************************************************************************

!<subroutine>

  pure subroutine spsor_aux_precUAB_3D(n,Du1,Du2,Du3,Dup,KldA,KcolA,KdiagA,&
                                       DA1,DA2,DA3,dsfA1,dsfA2,dsfA3,KldB,KcolB,&
                                       DB1,DB2,DB3,dsfB1,dsfB2,dsfB3,domega)

!<description>
  ! INTERNAL AUXILIARY ROUTINE:
  !
  ! Solves the following system:
  !
  ! <verb>
  !               ( U1  0  0  ) * ( u1 ) = ( f1 ) - omega * ( B1 ) * ( p )
  !               ( 0   U2 0  )   ( u2 )   ( f2 )           ( B2 )
  !               ( 0   0  U3 )   ( u3 )   ( f3 )           ( B3 )
  ! </verb>
  !
  ! where <tex>$ Ui := id + ( 1/omega * diag(Ai) )^{-1} * ltri(Ai) $</tex>
  !
!</description>

!<inputoutput>
  ! On entry, the right-hand-side (f1,f2,f3) of the system
  ! On exit, the solution (u1,u2,u3) of the system
  real(DP), dimension(*), intent(inout) :: Du1, Du2, Du3
!</inputoutput>

!<input>
  ! The dimension of A
  integer, intent(in) :: n

  ! The input vector p
  real(DP), dimension(*), intent(inout) :: Dup
  
  ! The matrix structure of Ai
  integer, dimension(*), intent(in) :: KldA, KcolA, KdiagA
  
  ! The matrix data of Ai
  real(DP), dimension(*), intent(in) :: DA1, DA2, DA3

  ! The scaling factor of Ai
  real(DP), intent(in) :: dsfA1, dsfA2, dsfA3
  
  ! The matrix structure of Bi
  integer, dimension(*), intent(in) :: KldB, KcolB
  
  ! The matrix data of Bi
  real(DP), dimension(*), intent(in) :: DB1, DB2, DB3

  ! The scaling factor of Bi
  real(DP), intent(in) :: dsfB1, dsfB2, dsfB3
  
  ! The relaxation parameter in range (0,2)
  real(DP), intent(in) :: domega
!</input>

!</subroutine>

  integer :: i,j,k,m
  real(DP) :: dt,daux1,daux2,daux3,dauxb1,dauxb2,dauxb3,dsb1,dsb2,dsb3
  
    ! Pre-calculate scaling factor
    dsb1 = (domega*dsfB1) / dsfA1
    dsb2 = (domega*dsfB2) / dsfA2
    dsb3 = (domega*dsfB3) / dsfA3
  
    ! Loop over all rows
    do i = n, 1, -1
    
      dauxb1 = 0.0_DP
      dauxb2 = 0.0_DP
      dauxb3 = 0.0_DP
      do j = KldB(i), KldB(i+1)-1
        dt = Dup(KcolB(j))
        dauxb1 = dauxb1 + DB1(j)*dt
        dauxb2 = dauxb2 + DB2(j)*dt
        dauxb3 = dauxb3 + DB3(j)*dt
      end do
    
      k = KdiagA(i)
    
      daux1 = 0.0_DP
      daux2 = 0.0_DP
      daux3 = 0.0_DP
      do j = KldA(i+1)-1, k+1, -1
        m = KcolA(j)
        daux1 = daux1 + DA1(j)*Du1(m)
        daux2 = daux2 + DA2(j)*Du2(m)
        daux3 = daux3 + DA3(j)*Du3(m)
      end do
      
      Du1(i) = Du1(i) - (domega*daux1 + dsb1*dauxb1) / DA1(k)
      Du2(i) = Du2(i) - (domega*daux2 + dsb2*dauxb2) / DA2(k)
      Du3(i) = Du3(i) - (domega*daux3 + dsb3*dauxb3) / DA3(k)
    
    end do ! i
    
  end subroutine ! spsor_aux_precUAB_3D

  ! ***************************************************************************

!<subroutine>

  pure subroutine spsor_aux_precUAAB(n,Du,Du2,Dup,KldA,KcolA,KdiagA,DA,dsfA,&
                                     KldA2,KcolA2,DA2,dsfA2,KldB,KcolB,DB,&
                                     dsfB,domega)

!<description>
  ! INTERNAL AUXILIARY ROUTINE:
  !
  ! Solves the following system:
  !
  ! <verb>
  !              U * u = f - A2 * u2 - B * p
  ! </verb>
  !
  ! where <tex>$ U := id + ( 1/omega * diag(A) )^{-1} * utri(A) $</tex>
  !
!</description>

!<inputoutput>
  ! On entry, the right-hand-side f of the system
  ! On exit, the solution u of the system
  real(DP), dimension(*), intent(inout) :: Du
!</inputoutput>

!<input>
  ! The dimension of Ai
  integer, intent(in) :: n

  ! The input vector u2
  real(DP), dimension(*), intent(in) :: Du2, Dup
  
  ! The matrix structure of A
  integer, dimension(*), intent(in) :: KldA, KcolA, KdiagA
  
  ! The matrix data of A
  real(DP), dimension(*), intent(in) :: DA

  ! The scaling factors of A
  real(DP), intent(in) :: dsfA
  
  ! The matrix structure of A2
  integer, dimension(*), intent(in) :: KldA2, KcolA2
  
  ! The matrix data of A2
  real(DP), dimension(*), intent(in) :: DA2

  ! The scaling factors of A2
  real(DP), intent(in) :: dsfA2
  
  ! The matrix structure of B
  integer, dimension(*), intent(in) :: KldB, KcolB
  
  ! The matrix data of B
  real(DP), dimension(*), intent(in) :: DB

  ! The scaling factor of B
  real(DP), intent(in) :: dsfB

  ! The relaxation parameter in range (0,2)
  real(DP), intent(in) :: domega
!</input>

!</subroutine>

  integer :: i,j,k
  real(DP) :: dom,daux,daux2,dauxb
  
    ! pre-calculate weight factor
    dom = domega / dsfA
  
    ! Loop over all rows
    do i = 1, n
    
      dauxb = 0.0_DP
      do j = KldB(i), KldB(i+1)-1
        dauxb = dauxb + DB(j)*Dup(KcolB(j))
      end do
    
      daux2 = 0.0_DP
      do j = KldA2(i), KldA2(i+1)-1
        daux2 = daux2 + DA2(j)*Du2(KcolA2(j))
      end do

      k = KdiagA(i)
    
      daux = 0.0_DP
      do j = KldA(i), k-1
        daux = daux + DA(j)*Du(KcolA(j))
      end do
      
      Du(i) = dom*(Du(i) - dsfA*daux - dsfA2*daux2 - dsfB*dauxb) / DA(k)
    
    end do ! i
    
  end subroutine ! spsor_aux_precUAAB

  ! ***************************************************************************

!<subroutine>

  pure subroutine spsor_aux_precUAAAB(n,Du,Du2,Du3,Dup,KldA,KcolA,KdiagA,DA,&
                                      dsfA,KldA2,KcolA2,DA2,DA3,dsfA2,dsfA3,&
                                      KldB,KcolB,DB,dsfB,domega)

!<description>
  ! INTERNAL AUXILIARY ROUTINE:
  !
  ! Solves the following system:
  !
  ! <verb>
  !              U * u = f - A2 * u2 - A3 * u3 - B * p
  ! </verb>
  !
  ! where <tex>$ U := id + ( 1/omega * diag(A) )^{-1} * utri(A) $</tex>
  !
!</description>

!<inputoutput>
  ! On entry, the right-hand-side f of the system
  ! On exit, the solution u of the system
  real(DP), dimension(*), intent(inout) :: Du
!</inputoutput>

!<input>
  ! The dimension of Ai
  integer, intent(in) :: n

  ! The input vector u2
  real(DP), dimension(*), intent(in) :: Du2, Du3, Dup
  
  ! The matrix structure of A
  integer, dimension(*), intent(in) :: KldA, KcolA, KdiagA
  
  ! The matrix data of A
  real(DP), dimension(*), intent(in) :: DA

  ! The scaling factors of A
  real(DP), intent(in) :: dsfA
  
  ! The matrix structure of A2
  integer, dimension(*), intent(in) :: KldA2, KcolA2
  
  ! The matrix data of A2
  real(DP), dimension(*), intent(in) :: DA2, DA3

  ! The scaling factors of A2
  real(DP), intent(in) :: dsfA2, dsfA3
  
  ! The matrix structure of B
  integer, dimension(*), intent(in) :: KldB, KcolB
  
  ! The matrix data of B
  real(DP), dimension(*), intent(in) :: DB

  ! The scaling factor of B
  real(DP), intent(in) :: dsfB

  ! The relaxation parameter in range (0,2)
  real(DP), intent(in) :: domega
!</input>

!</subroutine>

  integer :: i,j,k,m
  real(DP) :: dom,daux,daux2,daux3,dauxb
  
    ! pre-calculate weight factor
    dom = domega / dsfA
  
    ! Loop over all rows
    do i = 1, n
    
      dauxb = 0.0_DP
      do j = KldB(i), KldB(i+1)-1
        dauxb = dauxb + DB(j)*Dup(KcolB(j))
      end do
    
      daux2 = 0.0_DP
      daux3 = 0.0_DP
      do j = KldA2(i), KldA2(i+1)-1
        m = KcolA2(j)
        daux2 = daux2 + DA2(j)*Du2(m)
        daux3 = daux3 + DA3(j)*Du3(m)
      end do

      k = KdiagA(i)
    
      daux = 0.0_DP
      do j = KldA(i), k-1
        daux = daux + DA(j)*Du(KcolA(j))
      end do
      
      Du(i) = dom*(Du(i) - dsfA*daux - dsfA2*daux2 - dsfA3*daux3 &
                         - dsfB*dauxb) / DA(k)
    
    end do ! i
    
  end subroutine ! spsor_aux_precUAAAB

  ! ***************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine spsor_aux_precDDS(n,Du1,Du2,Dup,KldD,KcolD,DD1,DD2,&
                               dsfD1,dsfD2,DS,domega)

!<description>
  ! INTERNAL AUXILIARY ROUTINE:
  ! Calculates:
  !
  ! <verb>
  !                p := relax * S * (p - D1 * u1 - D2 * u2)
  ! </verb>
  !
  ! where S is a diagonal matrix, and D1 and D2 share the same structure
!</description>

!<inputoutput>
  real(DP), dimension(*), intent(inout) :: Dup
!</inputoutput>

!<input>
  integer, intent(in) :: n
  
  real(DP), dimension(*), intent(in) :: Du1, Du2
  
  integer, dimension(*), intent(in) :: KldD, KcolD
  real(DP), dimension(*), intent(in) :: DD1, DD2
  real(DP), intent(in) :: dsfD1, dsfD2

  real(DP), dimension(*), intent(in) :: DS
  
  real(DP), intent(in) :: domega
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
      Dup(i) = domega*DS(i)*(Dup(i) - dsfD1*daux1 - dsfD2*daux2)
    
    end do ! i
    !$omp end parallel do
    
  end subroutine ! spsor_aux_precDDS

  ! ***************************************************************************

!<subroutine>

#ifndef USE_OPENMP
  pure &
#endif

  subroutine spsor_aux_precDDDS(n,Du1,Du2,Du3,Dup,KldD,KcolD,&
                                DD1,DD2,DD3,dsfD1,dsfD2,dsfD3,DS,domega)

!<description>
  ! INTERNAL AUXILIARY ROUTINE:
  ! Calculates:
  !
  ! <verb>
  !                p := relax * S * (p - D1 * u1 - D2 * u2 - D3 * u3)
  ! </verb>
  !
  ! where S is a diagonal matrix, and D1, D2 and D3 share the same structure
!</description>

!<inputoutput>
  real(DP), dimension(*), intent(inout) :: Dup
!</inputoutput>

!<input>
  integer, intent(in) :: n
  
  real(DP), dimension(*), intent(in) :: Du1, Du2, Du3
  
  integer, dimension(*), intent(in) :: KldD, KcolD
  real(DP), dimension(*), intent(in) :: DD1, DD2, DD3
  real(DP), intent(in) :: dsfD1, dsfD2, dsfD3

  real(DP), dimension(*), intent(in) :: DS
  
  real(DP), intent(in) :: domega
!</input>

!</subroutine>

  integer :: i,j,k
  real(DP) :: daux1, daux2, daux3
  
    ! Loop over all rows
    !$omp parallel do private(j,k,daux1,daux2,daux3) if(n .gt. 1000)
    do i = 1, n
    
      ! Calculate D1*u and D2*v
      daux1 = 0.0_DP
      daux2 = 0.0_DP
      daux3 = 0.0_DP
      do j = KldD(i), KldD(i+1)-1
        k = KcolD(j)
        daux1 = daux1 + DD1(j)*Du1(k)
        daux2 = daux2 + DD2(j)*Du2(k)
        daux3 = daux3 + DD3(j)*Du3(k)
      end do
      
      ! Solve
      Dup(i) = domega*DS(i)*(Dup(i) - dsfD1*daux1 - dsfD2*daux2 - dsfD3*daux3)
    
    end do ! i
    !$omp end parallel do
    
  end subroutine ! spsor_aux_precDDDS

end module
