!##############################################################################
!# ****************************************************************************
!# <name> schurapprox </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains a set of routines which calculate an approximation
!# of a Schur-complement matrix S of a saddle-point matrix of the style
!#
!# <verb>
!#                             ( A  B )
!#                             ( D  0 )
!# </verb>
!#
!# where the exact Schur-complement matrix S is given as
!#
!#             <tex>$$     S := -D * A^{-1} * B   $$</tex>
!#
!# See the comment block below the routine list for more detailed information
!# on how exactly the approximation is assembled, as there are multiple
!# different variants available.
!#
!# The following routines can be found in this module:
!#
!# 1.) schur_assembleApprox2D
!#     -> Assembles a Schur-complement approximation matrix for a
!#        2D saddle-point system.
!#
!# 2.) schur_initPerfConfig
!#     -> Initialises the global performance configuration
!#
!#
!# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\\
!# Detailed information on the approximation assembly variants\\
!# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\\
!# Basically, all currently implemented algorithms calculate the approximation
!# of S by replacing <tex>$A^{-1}$</tex> by an approximation.
!#
!# All the currently implemented algorithms also support matrices S which have
!# an extended stencil, however, it is currently not clear under which
!# circumstances the approximation using an extended stencil is more accurate
!# than an approximation with the 'normal' stencil.
!#
!#
!#  SCHUR_TYPE_MAIN_DIAGONAL  \\
!# -------------------------- \\
!# This variant approximates <tex>$A^{-1}$</tex> by the inverse of the main diagonal of A.
!# This is the most simple and cheapest variant, which should do a good job
!# for Stokes systems.
!#
!#
!#  SCHUR_TYPE_LUMP_DIAGONAL  \\
!# -------------------------- \\
!# This variant approximates <tex>$A^{-1}$</tex> by the inverse of the lumped diagonal
!# matrix of A, i.e. a diagonal matrix L with
!#
!#                                  
!#          <tex>$$       L_{ii} :=  \sum_{j=1}^n  a_ij      $$</tex>
!#                                 
!#
!# This variant calculates the L_ii 'on the fly' without modifying the original
!# matrix A. This variant should only be used if A is (mostly) a mass matrix.
!#
!# Remark:
!# Even if the off-diagonal block sub-matrices of A (e.g. A12) exist, they are
!# ignored by this 'lumping' variant, i.e. the lumping is only performed on the
!# main diagonal blocks A11, A22 and (in 3D) A33 of A.
!#
!#
!#  SCHUR_TYPE_LOCAL_INVERSE  \\
!# -------------------------- \\
!# THIS ONE IS EXPERIMENTAL, USE AT YOUR OWN RISK!
!# This variant is a little difficult to explain, as it approximates <tex>$A^{-1}$</tex>
!# by local inverse matrices. The algorithm assembles the entries of S by
!# proceeding row-wise, and for each row i of S, there is a different local
!# matrix used, which we shall denote by T_i here.
!#
!# Assume that the algorithm is currently processing row i of S, then the
!# entry S(i,j) is defined as:
!#
!#            <tex>$$ S(i,j) := D(i,:) * T_i^{-1} * B(:,j)  $$</tex>
!#
!# where T_i is the local m by m sub-matrix of A, with m being the number of
!# non-zero entries in row i of D, where T_i contains all non-zero entries A_kl
!# of A, such that k and l are non-zero entries in row i of D, i.e. D(i,k) and
!# D(i,l) are non-zero entries in the sparsity pattern of D.
!#
!# This variant is (much) more expensive than the diagonal versions, however,
!# it *might* provide a better approximation of S.
!#
!# Remark:
!# Comparing the 'main diagonal' and 'local inverse' variants here is similar
!# to the comparision of 'diagonal' and 'full' Vanka or SP-SOR preconditioners;
!# see the vanka.f90 module for more information on this.
!#
!# </purpose>
!##############################################################################

module schurapprox

  use fsystem
  use genoutput
  use linearsystemscalar
  use perfconfig
  
  implicit none
  
  private
  
  public :: t_schurApprox
  public :: schur_initPerfConfig
  public :: schur_assembleApprox2D

!<constants>

!<constantblock>

  ! Approximate <tex>$A^{-1}$</tex> by the inverse of the main diagonal of A.
  integer, parameter, public :: SCHUR_TYPE_MAIN_DIAGONAL = 1

  ! Approximate <tex>$A^{-1}$</tex> by the inverse of the lumped diagonal of A.
  integer, parameter, public :: SCHUR_TYPE_LUMP_DIAGONAL = 2

  ! EXPERIMENTAL: Approximate <tex>$A^{-1}$</tex> by local inverse matrices.
  integer, parameter, public :: SCHUR_TYPE_LOCAL_INVERSE = 10

!</constantblock>

!</constants>

!<types>

!<typeblock>

  type t_schurApprox

    ! A pointer to a scalar matrix containing A11.
    type(t_matrixScalar), pointer :: p_rmatrixA11 => null()

    ! A pointer to a scalar matrix containing A22.
    type(t_matrixScalar), pointer :: p_rmatrixA22 => null()

    ! A pointer to a scalar matrix containing A33 (3D only).
    type(t_matrixScalar), pointer :: p_rmatrixA33 => null()

    ! A pointer to a scalar matrix containing B1.
    type(t_matrixScalar), pointer :: p_rmatrixB1 => null()

    ! A pointer to a scalar matrix containing B2.
    type(t_matrixScalar), pointer :: p_rmatrixB2 => null()

    ! A pointer to a scalar matrix containing B3 (3D only).
    type(t_matrixScalar), pointer :: p_rmatrixB3 => null()

    ! A pointer to a scalar matrix containing D1.
    type(t_matrixScalar), pointer :: p_rmatrixD1 => null()

    ! A pointer to a scalar matrix containing D2.
    type(t_matrixScalar), pointer :: p_rmatrixD2 => null()

    ! A pointer to a scalar matrix containing D3 (3D only).
    type(t_matrixScalar), pointer :: p_rmatrixD3 => null()

    ! A pointer to a scalar matrix containing S.
    type(t_matrixScalar), pointer :: p_rmatrixS => null()

    ! A factor that the approximation of S is to be multiplied by when adding
    ! it onto the content of S.
    real(DP) :: dtheta = 1.0_DP

    ! One of the SCHUR_TYPE_XXXX constants which specifies how the approximation
    ! is to be calculated.
    integer :: ctype = SCHUR_TYPE_MAIN_DIAGONAL

    ! Specifies whether S is to be cleared before assembling the approximation
    ! of the Schur-complement matrix.
    logical :: bclearS = .true.

  end type

!</typeblock>

!</types>

  !*****************************************************************************
  
  ! global performance configuration
  type(t_perfconfig), target, save :: schur_perfconfig

  !*****************************************************************************

contains

  ! ****************************************************************************

!<subroutine>

  subroutine schur_initPerfConfig(rperfconfig)

!<description>
  ! This routine initialises the global performance configuration
!</description>

!<input>
  ! OPTIONAL: performance configuration that should be used to initialise
  ! the global performance configuration. If not present, the values of
  ! the legacy constants is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>
!</subroutine>

    if (present(rperfconfig)) then
      schur_perfconfig = rperfconfig
    else
      call pcfg_initPerfConfig(schur_perfconfig)
    end if
  
  end subroutine schur_initPerfConfig
  
  ! ***************************************************************************

!<subroutine>

  subroutine schur_assembleApprox2D(rschur, rperfconfig)

!<description>
  ! TODO
!</description>

!<input>
  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! A t_schurApprox structure that contains the information of how the
  ! Schur-complement matrix approximation is to be assembled.
  type(t_schurApprox), intent(inout) :: rschur
!</inputoutput>

!</subroutine>

  ! The scaling factor for the approximation.
  real(DP) :: dtheta

  ! The scaling factors for the corresponding parts of the assembly.
  real(DP) :: dsf1, dsf2

  ! Arrays for the matrix structures
  integer, dimension(:), pointer :: p_IrowA, p_IcolA, &!p_IrowA12, p_IcolA12,&
      p_IdiagA, p_IrowB, p_IcolB, p_IrowD, p_IcolD, p_IrowS, p_IcolS

  ! Arrays for the matrix entries
  real(DP), dimension(:), pointer :: p_DA11, p_DA22, &!p_DA12, p_DA21, &
      p_DB1, p_DB2, p_DD1, p_DD2, p_DS

  ! The assembly type
  integer :: ctype

  ! The degree of D
  integer :: ndegree, neq

  

  

    ! First of all, let us make sure that all necessary matrices are present.
    if(.not. associated(rschur%p_rmatrixA11)) then
      call output_line ('Submatrix A11 is missing!', &
                        OU_CLASS_ERROR,OU_MODE_STD, 'schur_assembleApprox2D')
      call sys_halt()
    end if
    if(.not. associated(rschur%p_rmatrixA22)) then
      call output_line ('Submatrix A22 is missing!', &
                        OU_CLASS_ERROR,OU_MODE_STD, 'schur_assembleApprox2D')
      call sys_halt()
    end if
    if(.not. associated(rschur%p_rmatrixB1)) then
      call output_line ('Submatrix B1 is missing!', &
                        OU_CLASS_ERROR,OU_MODE_STD, 'schur_assembleApprox2D')
      call sys_halt()
    end if
    if(.not. associated(rschur%p_rmatrixB2)) then
      call output_line ('Submatrix B2 is missing!', &
                        OU_CLASS_ERROR,OU_MODE_STD, 'schur_assembleApprox2D')
      call sys_halt()
    end if
    if(.not. associated(rschur%p_rmatrixD1)) then
      call output_line ('Submatrix D1 is missing!', &
                        OU_CLASS_ERROR,OU_MODE_STD, 'schur_assembleApprox2D')
      call sys_halt()
    end if
    if(.not. associated(rschur%p_rmatrixD2)) then
      call output_line ('Submatrix D2 is missing!', &
                        OU_CLASS_ERROR,OU_MODE_STD, 'schur_assembleApprox2D')
      call sys_halt()
    end if
    if(.not. associated(rschur%p_rmatrixS)) then
      call output_line ('Matrix S is missing!', &
                        OU_CLASS_ERROR,OU_MODE_STD, 'schur_assembleApprox2D')
      call sys_halt()
    end if

    ! Make sure A11 and A22 are compatible
    if(.not. lsyssc_isMatrixStructureShared(rschur%p_rmatrixA11, &
                                            rschur%p_rmatrixA22)) then
      call output_line ('Submatrices A11 and A22 must have the same structure!', &
                        OU_CLASS_ERROR,OU_MODE_STD, 'schur_assembleApprox2D')
      call sys_halt()
    end if
    if((iand(rschur%p_rmatrixA11%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0) .or. &
       (iand(rschur%p_rmatrixA22%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0)) then
      call output_line ('A11 and A22 must not be transposed!', &
                        OU_CLASS_ERROR,OU_MODE_STD, 'schur_assembleApprox2D')
      call sys_halt()
    end if

    ! Make sure B1 and B2 are compatible
    if(.not. lsyssc_isMatrixStructureShared(rschur%p_rmatrixB1, &
                                            rschur%p_rmatrixB2)) then
      call output_line ('Submatrices B1 and B2 must have the same structure!', &
                        OU_CLASS_ERROR,OU_MODE_STD, 'schur_assembleApprox2D')
      call sys_halt()
    end if
    if((iand(rschur%p_rmatrixB1%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0) .or. &
       (iand(rschur%p_rmatrixB2%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0)) then
      call output_line ('B1 and B2 must not be transposed!', &
                        OU_CLASS_ERROR,OU_MODE_STD, 'schur_assembleApprox2D')
      call sys_halt()
    end if

    ! Make sure D1 and D2 are compatible
    if(.not. lsyssc_isMatrixStructureShared(rschur%p_rmatrixD1, &
                                            rschur%p_rmatrixD2)) then
      call output_line ('Submatrices D1 and D2 must have the same structure!', &
                        OU_CLASS_ERROR,OU_MODE_STD, 'schur_assembleApprox2D')
      call sys_halt()
    end if
    if((iand(rschur%p_rmatrixD1%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0) .or. &
       (iand(rschur%p_rmatrixD2%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0)) then
      call output_line ('D1 and D2 must not be transposed!', &
                        OU_CLASS_ERROR,OU_MODE_STD, 'schur_assembleApprox2D')
      call sys_halt()
    end if

    ! Make sure A, B and D are compatible
    if(rschur%p_rmatrixA11%NEQ .ne. rschur%p_rmatrixB1%NEQ) then
      call output_line ('B1/B2 submatrices are incompatible to A11/A22!', &
                        OU_CLASS_ERROR,OU_MODE_STD, 'schur_assembleApprox2D')
      call sys_halt()
    end if
    if(rschur%p_rmatrixA11%NCOLS .ne. rschur%p_rmatrixD1%NCOLS) then
      call output_line ('D1/D2 submatrices are incompatible to A11/A22!', &
                        OU_CLASS_ERROR,OU_MODE_STD, 'schur_assembleApprox2D')
      call sys_halt()
    end if
    if(rschur%p_rmatrixD1%NEQ .ne. rschur%p_rmatrixB1%NCOLS) then
      call output_line ('D1/D2 submatrices are incompatible to B1/B2!', &
                        OU_CLASS_ERROR,OU_MODE_STD, 'schur_assembleApprox2D')
      call sys_halt()
    end if

    ! Make sure S is valid
    if(iand(rschur%p_rmatrixS%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0) then
      call output_line ('S must not be transposed!', &
                        OU_CLASS_ERROR,OU_MODE_STD, 'schur_assembleApprox2D')
      call sys_halt()
    end if

    ! Make sure S is compatible to A, B and D
    if((rschur%p_rmatrixS%NEQ   .ne. rschur%p_rmatrixD1%NEQ) .or. &
       (rschur%p_rmatrixS%NCOLS .ne. rschur%p_rmatrixB1%NCOLS)) then
      call output_line ('S is incompatible to saddle point system!', &
                        OU_CLASS_ERROR,OU_MODE_STD, 'schur_assembleApprox2D')
      call sys_halt()
    end if

    ! Okay, the matrices are given and seem to be valid.

    ! Do we need to clear S?
    if(rschur%bclearS) call lsyssc_clearMatrix(rschur%p_rmatrixS)

    ! If theta is zero, then we have nothing more to do...
    dtheta = rschur%dtheta
    if(dtheta .eq. 0.0_DP) return

    ! Okay, let us calculate the scaling factors for the diagonal assembly
    ! types.
    dsf1 = dtheta * (rschur%p_rmatrixD1%dscaleFactor * &
                     rschur%p_rmatrixB1%dscaleFactor) / &
                    (rschur%p_rmatrixA11%dscaleFactor * &
                     rschur%p_rmatrixS%dscaleFactor)
    dsf2 = dtheta * (rschur%p_rmatrixD2%dscaleFactor * &
                     rschur%p_rmatrixB2%dscaleFactor) / &
                    (rschur%p_rmatrixA22%dscaleFactor * &
                     rschur%p_rmatrixS%dscaleFactor)

    ! Fetch the arrays of the necessary matrices
    call lsyssc_getbase_Kld (rschur%p_rmatrixA11, p_IrowA)
    call lsyssc_getbase_Kcol(rschur%p_rmatrixA11, p_IcolA)
    call lsyssc_getbase_Kdiagonal(rschur%p_rmatrixA11, p_IdiagA)
    call lsyssc_getbase_double(rschur%p_rmatrixA11, p_DA11)
    call lsyssc_getbase_double(rschur%p_rmatrixA22, p_DA22)

    call lsyssc_getbase_Kld (rschur%p_rmatrixB1, p_IrowB)
    call lsyssc_getbase_Kcol(rschur%p_rmatrixB1, p_IcolB)
    call lsyssc_getbase_double(rschur%p_rmatrixB1, p_DB1)
    call lsyssc_getbase_double(rschur%p_rmatrixB2, p_DB2)

    call lsyssc_getbase_Kld (rschur%p_rmatrixD1, p_IrowD)
    call lsyssc_getbase_Kcol(rschur%p_rmatrixD1, p_IcolD)
    call lsyssc_getbase_double(rschur%p_rmatrixD1, p_DD1)
    call lsyssc_getbase_double(rschur%p_rmatrixD2, p_DD2)

    call lsyssc_getbase_Kld (rschur%p_rmatrixS, p_IrowS)
    call lsyssc_getbase_Kcol(rschur%p_rmatrixS, p_IcolS)
    call lsyssc_getbase_double(rschur%p_rmatrixS, p_DS)

    ! Fetch the number of equations in S
    neq = rschur%p_rmatrixS%NEQ

    ! Fetch the assembly type
    ctype = rschur%ctype

    ! Okay, now what are we going to assemble?
    select case(ctype)
    case(SCHUR_TYPE_MAIN_DIAGONAL)

      ! This one is easy.
      call schur_mainDiagonal2D(neq, p_IdiagA, p_DA11, p_DA22, p_IrowB, &
          p_IcolB, p_DB1, p_DB2, p_IrowD, p_IcolD, p_DD1, p_DD2, p_IrowS, &
          p_IcolS, p_DS, dsf1, dsf2, rperfconfig)

    case(SCHUR_TYPE_LUMP_DIAGONAL)

      ! This one is also easy.
      call schur_lumpDiagonal2D(neq, p_IrowA, p_DA11, p_DA22, p_IrowB, &
          p_IcolB, p_DB1, p_DB2, p_IrowD, p_IcolD, p_DD1, p_DD2, p_IrowS, &
          p_IcolS, p_DS, dsf1, dsf2, rperfconfig)

    case(SCHUR_TYPE_LOCAL_INVERSE)

      ! Okay, this one is a little more interesting.

!      celemVelo = 0_I32
!      celemPres = 0_I32
!
!      ! Fetch the discretisations of B (if they exist)
!      p_rdiscrVelo => p_rB%RmatrixBlock(1,1)%p_rspatialDiscrTest
!      p_rdiscrPres => p_rB%RmatrixBlock(1,1)%p_rspatialDiscrTrial
!      if(associated(p_rdiscrVelo) .and. associated(p_rdiscrPres)) then
!
!        ! Check whether both discretisations are uniform
!        if((p_rdiscrVelo%inumFESpaces .eq. 1) .and. &
!           (p_rdiscrPres%inumFESpaces .eq. 1)) then
!
!          ! Yes, they are uniform, so get the element identifiers.
!          celemVelo = p_rdiscrVelo%RelementDistr(1)%celement
!          celemPres = p_rdiscrPres%RelementDistr(1)%celement
!
!        end if
!
!      end if
!
!      ! Now let us see if there is an interesting combination of velocity and
!      ! pressure spaces that we have a special implementation for. If not, then
!      ! there still is a 'generic' routine which we can call which performs
!      ! exactly the same job. However, the 'generic' routine is much slower
!      ! than the specialised versions...
!
!      if(.false.) then
!
!        ! ...
!
!      else

        ! Okay, we do not have a specialised version. So we will call the 'generic'
        ! routine, but first we need to calculate the degree of D - this is the
        ! maximum number of velocity DOFs which are adjacent to one pressure DOF.
        ! This information is needed by the 'generic' routine.
        ndegree = schur_aux_calcDegree(p_IrowD)

        ! And call the routine
        call schur_localInverse2D(neq, ndegree, p_IrowA, p_IcolA, &
            p_DA11, p_DA22, p_IrowB, p_IcolB, p_DB1, p_DB2, p_IrowD, p_IcolD, &
            p_DD1, p_DD2, p_IrowS, p_IcolS, p_DS, dsf1, dsf2)

!      end if

    case default

      ! stupid caller...
      call output_line ('Invalid assembly type!', &
                        OU_CLASS_ERROR,OU_MODE_STD, 'schur_assembleApprox2D')
      call sys_halt()

    end select

    ! That is it

  end subroutine

  ! ***************************************************************************

!<function>

  integer function schur_aux_calcDegree(p_IrowPtr) result(ndegree)

!<description>
  ! AUXILIARY FUNCTION:
  ! This function calculates the degree of the adjacency graph of a CSR matrix.
!</description>

!<input>
  ! A scalar matrix from which the degree is to be calculated.
  integer, dimension(:), pointer :: p_IrowPtr
!</input>

!<result>
  ! The degree of the adjacency graph of the matrix.
!</result>

!</function>

  integer :: i

    ! Get the row-pointer array of the matrix
    ndegree = 0
    do i = 1, ubound(p_IrowPtr,1)-1
      ndegree = max(ndegree, p_IrowPtr(i+1)-p_IrowPtr(i))
    end do

  end function

  ! ***************************************************************************

!<subroutine>

  subroutine schur_mainDiagonal2D(n, p_IdiagA, p_DA11, p_DA22, p_IrowB, p_IcolB, &
             p_DB1, p_DB2, p_IrowD, p_IcolD, p_DD1, p_DD2, p_IrowS, p_IcolS, &
             p_DS, dsf1, dsf2, rperfconfig)

!<description>
  ! INTERNAL ROUTINE:
  ! Assembles a Schur-complement approximation for a 2D saddle-point system
  ! by approximating <tex>$A^{-1}$</tex> by the inverse of the main diagonal of A.
!</description>

!<input>
  ! The number of rows of S (= #DOFs in pressure space)
  integer, intent(in) :: n

  ! The diagonal pointer array of A
  integer, dimension(*), intent(in) :: p_IdiagA

  ! The data arrays of A11 and A22
  real(DP), dimension(*), intent(in) :: p_DA11, p_DA22

  ! The matrix structure of B
  integer, dimension(*), intent(in) :: p_IrowB, p_IcolB

  ! The data arrays of B1 and B2
  real(DP), dimension(*), intent(in) :: p_DB1, p_DB2

  ! The matrix structure of D
  integer, dimension(*), intent(in) :: p_IrowD, p_IcolD

  ! The data arrays of D1 and D2
  real(DP), dimension(*), intent(in) :: p_DD1, p_DD2

  ! The matrix structure of S
  integer, dimension(*), intent(in) :: p_IrowS, p_IcolS

  ! The scaling factor for D1 * A11^{-1} * B1
  real(DP), intent(in) :: dsf1

  ! The scaling factor for D2 * A22^{-1} * B2
  real(DP), intent(in) :: dsf2

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! The data array of S
  real(DP), dimension(*), intent(inout) :: p_DS
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i,j,k,l,idxB,idxD,idxS,inextS
  real(DP) :: dDA1, dDA2

  ! Pointer to the performance configuration
  type(t_perfconfig), pointer :: p_rperfconfig
  
  if (present(rperfconfig)) then
    p_rperfconfig => rperfconfig
  else
    p_rperfconfig => schur_perfconfig
  end if

    ! Loop over all rows of S/D
    !$omp parallel do if(n > p_rperfconfig%NEQMIN_OMP)&
    !$omp private(j,k,l,idxB,idxD,idxS,inextS,dDA1,dDA2)
    do i = 1, n

      ! Loop over all non-zeroes of row i of D
      do j = p_IrowD(i), p_IrowD(i+1)-1

        ! Get the column index of D_ij
        idxD = p_IcolD(j)

        ! Precalculate D1_ij / A11_jj and D2 / A22_jj, scaled by the
        ! corresponding scaling factors.
        k = p_IdiagA(idxD)
        dDA1 = (dsf1 * p_DD1(j)) / p_DA11(k)
        dDA2 = (dsf2 * p_DD2(j)) / p_DA22(k)

        ! Reset next non-zero pointer of row i of S. See the comment inside the
        ! loop below for closer information.
        inextS = p_IrowS(i)

        ! Loop over all non-zeroes of row j of B
        do k = p_IrowB(idxD), p_IrowB(idxD+1)-1

          ! Get the column index of B_jk
          idxB = p_IcolB(k)

          ! We are current processing entry B_jk, so we are going to calculate
          ! D_ij * A_jj * B_jk. This is only interesting if S_ik exists in the
          ! sparsity pattern of S. So let us loop over all non-zeroes of row i
          ! of S now to check whether S_ik exists and, if it exists, update its
          ! entry.

          ! We will exploit a little trick here: Since the non-zero entries of a
          ! row of B and S are sorted in ascending order in respect to their
          ! column indices, we do not need to loop through all non-zero entries
          ! of row i of S for each non-zero entry in row j of B.
          do l = inextS, p_IrowS(i+1)-1

            ! Get the column index
            idxS = p_IcolS(l)

            ! Is this entry S_ik?
            if(idxS .eq. idxB) then

              ! Yes, S_ik exists, so update it:
              ! S_ik := S_ik - D_ij * A_jj^-1 * B_jk
              p_DS(l) = p_DS(l) - dDA1*p_DB1(k) - dDA2*p_DB2(k)

              ! Also, update the next non-zero pointer of S.
              inextS = l+1
              exit

            else if(idxS .gt. idxB) then

              ! If we come out here, then we are at entry S_il with l > k, so
              ! S_ik does not exist, so let us update the next non-zero pointer
              ! of S and exit the inner-most loop here.
              inextS = l
              exit

            end if

          end do ! l

        end do ! k

      end do ! j

    end do ! i
    !$omp end parallel do

    ! That is it

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine schur_lumpDiagonal2D(n, p_IrowA, p_DA11, p_DA22, p_IrowB, p_IcolB, &
             p_DB1, p_DB2, p_IrowD, p_IcolD, p_DD1, p_DD2, p_IrowS, p_IcolS, &
             p_DS, dsf1, dsf2, rperfconfig)

!<description>
  ! INTERNAL ROUTINE:
  ! Assembles a Schur-complement approximation for a 2D saddle-point system
  ! by approximating <tex>$A^{-1}$</tex> by the inverse of the lumped diagonal matrix of A.
!</description>

!<input>
  ! The number of rows of S (= #DOFs in pressure space)
  integer, intent(in) :: n

  ! The row-pointer array of A
  integer, dimension(*), intent(in) :: p_IrowA

  ! The data arrays of A11 and A22
  real(DP), dimension(*), intent(in) :: p_DA11, p_DA22

  ! The matrix structure of B
  integer, dimension(*), intent(in) :: p_IrowB, p_IcolB

  ! The data arrays of B1 and B2
  real(DP), dimension(*), intent(in) :: p_DB1, p_DB2

  ! The matrix structure of D
  integer, dimension(*), intent(in) :: p_IrowD, p_IcolD

  ! The data arrays of D1 and D2
  real(DP), dimension(*), intent(in) :: p_DD1, p_DD2

  ! The matrix structure of S
  integer, dimension(*), intent(in) :: p_IrowS, p_IcolS

  ! The scaling factor for D1 * A11^{-1} * B1
  real(DP), intent(in) :: dsf1

  ! The scaling factor for D2 * A22^{-1} * B2
  real(DP), intent(in) :: dsf2

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! The data array of S
  real(DP), dimension(*), intent(inout) :: p_DS
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i,j,k,l,idxB,idxD,idxS,inextS
  real(DP) :: dDA1, dDA2

  ! Pointer to the performance configuration
  type(t_perfconfig), pointer :: p_rperfconfig
  
  if (present(rperfconfig)) then
    p_rperfconfig => rperfconfig
  else
    p_rperfconfig => schur_perfconfig
  end if

    ! Loop over all rows of S/D
    !$omp parallel do if(n > p_rperfconfig%NEQMIN_OMP)&
    !$omp private(j,k,l,idxB,idxD,idxS,inextS,dDA1,dDA2)
    do i = 1, n

      ! Loop over all non-zeros of row i of D
      do j = p_IrowD(i), p_IrowD(i+1)-1

        ! Get the column index of D_ij
        idxD = p_IcolD(j)

        ! First of all, lump row j of A11/A22
        dDA1 = 0.0_DP
        dDA2 = 0.0_DP
        do k = p_IrowA(idxD), p_IrowA(idxD+1)-1
          dDA1 = dDA1 + p_DA11(k)
          dDA2 = dDA2 + p_DA22(k)
        end do

        ! Precalculate D1_ij / A11_jj and D2 / A22_jj, scaled by the
        ! corresponding scaling factors.
        dDA1 = (dsf1 * p_DD1(j)) / dDA1
        dDA2 = (dsf2 * p_DD2(j)) / dDA2

        ! Reset next non-zero pointer of row i of S.
        inextS = p_IrowS(i)

        ! Loop over all non-zeroes of row j of B
        do k = p_IrowB(idxD), p_IrowB(idxD+1)-1

          ! Get the column index of B_jk
          idxB = p_IcolB(k)

          do l = inextS, p_IrowS(i+1)-1

            ! Get the column index
            idxS = p_IcolS(l)

            ! Is this entry S_ik?
            if(idxS .eq. idxB) then

              ! Yes, S_ik exists, so update it:
              ! S_ik := S_ik - D_ij * A_jj^-1 * B_jk
              p_DS(l) = p_DS(l) - dDA1*p_DB1(k) - dDA2*p_DB2(k)

              ! Also, update the next non-zero pointer of S.
              inextS = l+1
              exit

            else if(idxS .gt. idxB) then

              inextS = l
              exit

            end if

          end do ! l

        end do ! k

      end do ! j

    end do ! i
    !$omp end parallel do

    ! That is it

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine schur_localInverse2D(n, ndegree, p_IrowA, p_IcolA, p_DA11, p_DA22, &
             p_IrowB, p_IcolB, p_DB1, p_DB2, p_IrowD, p_IcolD, p_DD1, p_DD2, &
             p_IrowS, p_IcolS, p_DS, dsf1, dsf2)

!<description>
  ! INTERNAL ROUTINE:
  ! Assembles a Schur-complement approximation for a 2D saddle-point system
  ! by approximating <tex>$A^{-1}$</tex> by local inverse sub-matrices.
!</description>

!<input>
  ! The number of rows of S (= #DOFs in pressure space)
  integer, intent(in) :: n

  ! The degree of the matrix D.
  integer, intent(in) :: ndegree

  ! The matrix structure of A
  integer, dimension(*), intent(in) :: p_IrowA, p_IcolA

  ! The data arrays of A11 and A22
  real(DP), dimension(*), intent(in) :: p_DA11, p_DA22

  ! The matrix structure of B
  integer, dimension(*), intent(in) :: p_IrowB, p_IcolB

  ! The data arrays of B1 and B2
  real(DP), dimension(*), intent(in) :: p_DB1, p_DB2

  ! The matrix structure of D
  integer, dimension(*), intent(in) :: p_IrowD, p_IcolD

  ! The data arrays of D1 and D2
  real(DP), dimension(*), intent(in) :: p_DD1, p_DD2

  ! The matrix structure of S
  integer, dimension(*), intent(in) :: p_IrowS, p_IcolS

  ! The scaling factor for D1 * A11^{-1} * B1
  real(DP), intent(in) :: dsf1

  ! The scaling factor for D2 * A22^{-1} * B2
  real(DP), intent(in) :: dsf2
!</input>

!<inputoutput>
  ! The data array of S
  real(DP), dimension(*), intent(inout) :: p_DS
!</inputoutput>

!</subroutine>

  integer :: i,j,k,l,m,idxA,idxB,idxD,idxS,inextD,inextS,irowD1,irowD2

  ! Local matrices
  real(DP), dimension(ndegree) :: DDA1, DDA2
  real(DP), dimension(ndegree,ndegree) :: DA11, DA22

  ! Column index array for current row of D
  integer, dimension(ndegree) :: IidxD

  ! A pivot array for LAPACK
  integer, dimension(ndegree) :: Ipivot
  integer :: info1, info2

    ! Remark:
    ! For an explaination of what 'inextD' and 'inextS' are used for here,
    ! see the comments in the routine 'schur_mainDiagonal2D'.

    ! Loop over all rows of S/D
    do i = 1, n

      ! Clear the local matrices
      !DDA1 = 0.0_DP
      !DDA2 = 0.0_DP
      DA11 = 0.0_DP
      DA22 = 0.0_DP

      ! Fetch the current row of D
      irowD1 = p_IrowD(i)
      irowD2 = p_IrowD(i+1)-1
      m = irowD2 - irowD1 + 1
      DDA1(1:m) = dsf1 * p_DD1(irowD1:irowD2)
      DDA2(1:m) = dsf2 * p_DD2(irowD1:irowD2)
      IidxD(1:m) = p_IcolD(irowD1:irowD2)

      ! Loop over all non-zeroes of row i of D
      do j = 1, m

        ! Get the column index of D_ij
        idxD = IidxD(j)

        ! Reset next non-zero pointer of D
        inextD = 1

        ! Loop over all non-zeroes of row j of A
        do k = p_IrowA(idxD), p_IrowA(idxD+1)-1

          ! Get the column index of A_jk
          idxA = p_IcolA(k)

          ! Check whether A_jk exists in our local A matrix and, if it exists,
          ! store the entry.
          do l = inextD, m

            if(IidxD(l) .eq. idxA) then

              ! Yes, so store A_jk to our local matrix.
              ! Remark: We store the local matrix in transposed format, as this
              ! will be needed by LAPACK's DGESV routine below!
              DA11(l,j) = p_DA11(k)
              DA22(l,j) = p_DA22(k)

              inextD = l+1
              exit

            else if(IidxD(l) .gt. idxA) then

              inextD = l
              exit

            end if

          end do ! l

        end do ! k

      end do ! j

      ! Okay, we now have row i of D and the local matrix of A stored in
      ! DDA1/DA2 and DA11/DA22. Now we are going to calculate D * A^-1.
      ! Since D is a row vector, we stored A in transposed format in the code
      ! above, so we will solve now (A^T)^-1 * D^T.
      call DGESV(m,1,DA11,ndegree,Ipivot,DDA1,ndegree,info1)
      call DGESV(m,1,DA22,ndegree,Ipivot,DDA2,ndegree,info2)

      ! Make sure LAPACK did not fail
      if((info1 .ne. 0) .or. (info2 .ne. 0)) then
        call output_line ('Local matrix singular!', &
                          OU_CLASS_ERROR,OU_MODE_STD, 'schur_localInverse2D')
        call sys_halt()
      end if

      ! Okay, loop over all non-zeroes of row D again.
      do j = 1, m

        ! Get the column index of D_ij
        idxD = IidxD(j)

        ! Reset next non-zero pointer of S
        inextS = p_IrowS(i)

        ! Loop over all non-zeroes of row j of B
        do k = p_IrowB(idxD), p_IrowB(idxD+1)-1

          ! Get the column index of B_jk
          idxB = p_IcolB(k)

          ! Try to find S_ik
          do l = inextS, p_IrowS(i+1)-1

            ! Get column index of S
            idxS = p_IcolS(l)

            ! Is this entry S_ik?
            if(idxS .eq. idxB) then

              ! Yes, so update it.
              p_DS(l) = p_DS(l) - DDA1(j)*p_DB1(k) - DDA2(j)*p_DB2(k)

              inextS = l+1
              exit

            else if(idxS .gt. idxB) then

              inextS = l
              exit

            end if

          end do ! l

        end do ! k

      end do ! j

    end do ! i

    ! That is it

  end subroutine

end module
