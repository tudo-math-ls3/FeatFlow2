!##############################################################################
!# ****************************************************************************
!# <name> spacepreconditioner </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the routines to solve the stationary core equation
!# of the problem with a nonlinear solver. The core equation aims to solve
!# the coupled KKT system of the minimisation problem
!#
!#   $$J(y)  :=  1/2 ||y-z||_{L_2}  +  \gamma/2 ||y(T)-z(T)||_{L_2}  +  \alpha_C||u||^2  ->  min! $$
!#
!# with $z$ being a given 'desired' flow field in the domain and
!# u being an unknown control.
!#
!# The discretised core equation reads at the moment:
!#
!#  $$        A_1 y   +  \eta_1 B p   +  \mu_1 M \lambda         = f_1 $$
!#  $$ \tau_1 B^T y   +  \kappa_1 I p                            = f_2 $$
!#
!#  $$        R_2 y   +  A_2 \lambda  + \eta_2   B \xi           = f_3 $$
!#  $$            \tau_2 B^T \lambda  + \kappa_2 I \xi           = f_4 $$
!#
!# with
!#
!#   $$ A_1 = \iota_1 I  +  \alpha_1 M  +  \theta_1 L  +  \gamma_1 N(y) + dnewton_1 N*(y)$$
!#   $$ A_2 = \iota_2 I  +  \alpha_2 M  +  \theta_2 L  +  \gamma_2 N(y) + dnewton_2 N*(y)$$
!#   $$ R_2 =               \mu_2    M  +                                 dr_2      N*(\lambda) $$
!#  
!# and
!#
!#   $I$     = identity matrix,
!#   $M$     = mass matrix,
!#   $L$     = Stokes matrix ($\nu$*Laplace),
!#   $N(y)$  = Nonlinearity including stabilisation, depending on the 
!#             primal velocity, i.e.
!#                   $$ (y\Delta)\cdot $$
!#   $N*(y)$ = Newton matrix, depending on the primal velocity, i.e.
!#                  $$ (\Delta y)\cdot $$
!#   $N*(\lambda)$ = Newton matrix, depending on the dual velocity, i.e.
!#                  $$ (\Delta \lambda)\cdot $$
!#   
!#   $\iota_i$  = 0/1     - switches the identity matrix on/off,
!#   $\alpha_i$ = 0/1     - switches the mass matrix on/off;
!#                          =0 for stationary problem,
!#   $\theta_í$           - weight for the Laplace matrix,
!#   $\gamma_i$ = 0/1     - Switches the nonlinearity on/off;
!#                          =0 for Stokes system,
!#   $dnewton_i \in R$    - Switches the Newton matrix on/off.
!#   $\eta_i$   = 0/1     - Switches the 'B'-term on/off,
!#   $\tau_i$   = 0/1     - Switches the 'B^T'-term on/off,
!#   $\mu_i$              - Weight for the 'coupling' mass matrix.
!#   $\kappa_i$ = 0/1     - Switches of the identity matrix I for the pressure
!#                          in the continuity equation
!#   $\dr_i \in R$        - Switches the 'reactive coupling mass matrix' on/off
!#                    
!# (y,p) is the velocity/pressure solution pair.
!# (lambda,xi) is the dual velocity/pressure solution.
!#
!# The core equation is abstractly written a nonlinear system of the form
!#
!#  $$ A(x)x = b $$
!#
!# and solved with the nonlinear solver from the kernel, using the defect
!# correction approach
!#
!#  $$  x_{n+1}  =  x_n  +  \omega_n C^{-1} ( b - A(x_n) x_n )  $$
!#
!# where $C^{-1}$ means to apply a suitable preconditioner (inverse mass
!# matrix, apply the linearised $A(x_n)^-1$ with multigrid, apply Newton or 
!# do something similar). 
!#
!# The following routines can be found here:
!#
!# 1.) cc_createSpacePreconditioner / cc_releasePreconditioner
!#     -> Creates/Releases a basic preconditioner structure for a spatial 
!#        preconditioner.
!#
!# 2.) cc_precondSpaceDefect
!#     -> Executes spatial preconditioning on a given defect vector
!#
!# </purpose>
!##############################################################################

module spacepreconditioner

  use fsystem
  use storage
  use boundary
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use paramlist
  use filtersupport
  use bilinearformevaluation
  use linearformevaluation
  use multilevelprojection
  use matrixmodification
  use linearsolver
  use linearsolverautoinitialise
  use matrixrestriction
  use trilinearformevaluation
  
  use collection
  use convection
    
  use cc2dmedium_callback
  
  use matrixio
  use vectorio
  
  use optcanalysis
  use spacematvecassembly
  
  implicit none
  
!<constants>

!<constantblock description="Preconditioner identifiers for the defect in the nonlinear iteration">

  ! No preconditioning
  integer, parameter :: CCPREC_NONE         = -1

  ! Preconditioning with inverse mass matrix (not yet implemented)
  integer, parameter :: CCPREC_INVERSEMASS   = 0

  ! Preconditioning by linear solver, solving the linearised system
  integer, parameter :: CCPREC_LINEARSOLVER  = 1

  ! Preconditioning by Newton-Iteration
  integer, parameter :: CCPREC_NEWTON        = 2
  
  ! Preconditioning by inexact/adaptive Newton iteration
  integer, parameter :: CCPREC_INEXACTNEWTON = 3

!</constantblock>

!</constants>

  
!<types>
  ! This type is used to save some situation specific assembly information
  ! during the setup phase of the nonlinear solver. Here it's noted, if
  ! and whose matrices exist and/or must be assmebled transposed to be
  ! compatible with the preconditioner and more. It's more or less
  ! a collection if different flags.
  type t_ccPreconditionerSpecials
  
    ! Whether to use 'adaptive matrices', i.e. set up coarse grid matrices
    ! with the help of fine grid matrices. This is used for very special
    ! discretisations only (e.g. Q1~/Q0). =0: deactivate
    integer :: iadaptiveMatrices    = 0
    
    ! A configuration parameter for adaptive matrices.
    real(DP) :: dadMatThreshold     = 0.0_DP

    ! If the preconditioner is a linear solver:
    ! Type of solver.
    ! =0: Gauss elimination (UMFPACK)
    ! =1: Multigrid solver
    integer :: isolverType = 0
    
    ! If the preconditioner is the linear multigrid solver:
    ! Type of smoother.
    ! =0: general VANKA (slow, but independent of the discretisation and of the problem)
    ! =1: general VANKA; 'direct' method, bypassing the defect correction approach.
    !     (-> specialised variant of 0, but slightly faster)
    ! =2: Simple Jacobi-like VANKA, 2D Navier Stokes problem, general discretisation
    !     (i.e. automatically chooses the best suitable VANKA variant).
    ! =3: Simple Jacobi-like VANKA, 2D Navier Stokes problem, general discretisation
    !     (i.e. automatically chooses the best suitable VANKA variant).
    !     'direct' method, bypassing the defect correction approach.
    !     (-> specialised variant of 8, but faster)
    ! =4: Full VANKA, 2D Navier Stokes problem, general discretisation
    !     (i.e. automatically chooses the best suitable VANKA variant).
    ! =5: Full VANKA, 2D Navier Stokes problem, general discretisation
    !     (i.e. automatically chooses the best suitable VANKA variant).
    !     'direct' method, bypassing the defect correction approach.
    !     (-> specialised variant of 10, but faster)
    integer :: ismootherType = 3
    
    ! If the preconditioner is the linear multigrid solver:
    ! Type of coarse grid solver.    
    ! =0: Gauss elimination (UMFPACK)
    ! =1: Defect correction with diagonal VANKA preconditioning.
    ! =2: BiCGStab with diagonal VANKA preconditioning
    integer :: icoarseGridSolverType = 1
        
    ! This flag is set to .TRUE. if there are no Neumann boundary
    ! components. In that case, the pressure matrices of direct
    ! solvers must be changed.
    logical :: bneedPressureDiagonalBlock = .false.
    
    ! Set to TRUE if the preconditioner needs virtually transposed B matrices
    ! as D matrices on all levels except for the coarse mesh.
    logical :: bneedVirtTransposedD = .false.
    
    ! Set to TRUE if the preconditioner needs virtually transposed B matrices
    ! as D matrices on the coarse mesh.
    logical :: bneedVirtTransposedDonCoarse = .false.
    
  end type

!</typeblock>

!<typeblock>

  ! Represents the core equation on one level of the discretisation.
  ! Collects all information that are necessary to assemble the 
  ! (linearised) system matrix and RHS vector.
  type t_cccoreEquationOneLevel
  
    ! The (linearised) system matrix for that specific level. 
    type(t_matrixBlock), pointer :: p_rmatrix => null()

    ! Reference to the static matrices on this level (Stokes, B,...)
    type(t_staticLevelInfo), pointer :: p_rstaticInfo

    ! Temporary vectors for the interpolation of a solution to a lower level.
    ! Exists only on levels NLMIN..NLMAX-1 !
    type(t_vectorBlock), pointer :: p_rtempVector1 => null()
    type(t_vectorBlock), pointer :: p_rtempVector2 => null()
    type(t_vectorBlock), pointer :: p_rtempVector3 => null()

    ! Block matrix, which is used in the defect correction / Newton
    ! algorithm as preconditioner matrix of the correspnding underlying
    ! linear sytem. Is usually the (linearise) system matrix or
    ! a Newton matrix. This matrix is changed during the
    ! nonlinear iteration and used e.g. if a linear solver (Multigrid) is
    ! used for preconditioning.
    type(t_matrixBlock), pointer :: p_rmatrixPreconditioner => null()
    
  end type

!</typeblock>

!<typeblock>

!<typeblock>

  ! Preconditioner structure for CCxD. This structure saves the configuration of the
  ! spatial preconditioner.
  
  type t_ccspatialPreconditioner
  
    ! Type of preconditioner.
    ! This is one of the CCPREC_xxxx flags as defined above (CCPREC_INVERSEMASS for
    ! preconditioning with inverse mass matrix, CCPREC_LINEARSOLVER for solving a linear
    ! system, CCPREC_NEWTON / CCPREC_INEXACTNEWTON for a Newton iteration,...)
    integer :: ctypePreconditioning = CCPREC_NONE
    
    ! Name of the section in the DAT file configuring this preconditioner.
    character(LEN=SYS_STRLEN) :: spreconditionerSection = ''
    
    ! Minimum discretisation level
    integer :: nlmin = 0
    
    ! Maximum discretisation level
    integer :: nlmax = 0
    
    ! A t_ccPreconditionerSpecials structure that saves information about
    ! special 'tweaks' in matrices such that everything works.
    type(t_ccPreconditionerSpecials) :: rprecSpecials
    
    ! An array of t_cccoreEquationOneLevel structures for all levels
    ! of the discretisation.
    type(t_cccoreEquationOneLevel), dimension(:), pointer :: RcoreEquation => null()

    ! Pointer to linear solver node if a linear solver is the preconditioner.
    ! (Thus, this applies for the defect correction and the Newton preconditioner).
    type(t_linsolNode), pointer :: p_rsolverNode
    
    ! An interlevel projection structure for changing levels
    type(t_interlevelProjectionBlock), dimension(:), pointer :: p_Rprojection
    
    ! Temporary scalar vector; used for calculating the nonlinear matrix
    ! on lower levels / projecting the solution from higher to lower levels.
    type(t_vectorScalar), pointer :: p_rtempVectorSc

    ! A filter chain that is used for implementing boundary conditions or other
    ! things when invoking the linear solver.
    type(t_filterChain), dimension(:), pointer :: p_RfilterChain
    
  end type

!</typeblock>

!</types>

contains
  
  ! ***************************************************************************

!<subroutine>
  
  subroutine vanka_NSOptC2D (rrhs,domega,rvector,rmatrix,IelementList,celementU,celementP)
  
!<description>
!</description>

!<input>
  ! The right-hand-side vector of the system
  type(t_vectorBlock), intent(in)         :: rrhs
  
  ! Matrix
  type(t_matrixBlock), intent(in) :: rmatrix
  
  ! Relaxation parameter. Standard=1.0_DP.
  real(DP), intent(in)                    :: domega
!</input>

!<inputoutput>
  ! The initial solution vector. Is replaced by a new iterate.
  type(t_vectorBlock), intent(inout)         :: rvector

  ! A list of element numbers where VANKA should be applied to.
  integer, dimension(:)     :: IelementList

  ! Element type for the FE spaces
  integer(I32), intent(in) :: celementU,celementP
!</inputoutput>

!</subroutine>

  ! Multiplication factors
  real(DP), dimension(6,6) :: Dmult
  
  ! Array for the pressure DOF's on the element
  integer, dimension(1) :: IelIdx2
  integer, dimension(:,:), allocatable :: IdofsP
  real(DP), dimension(:,:), allocatable :: DaInv
  integer, dimension(:), allocatable :: IdofsU
  
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
  real(DP), dimension(:,:), allocatable :: DS1, DS2
  integer, dimension(:), allocatable :: Ipiv

  ! Quick access for the vector arrays
  real(DP), dimension(:), pointer :: p_DrhsU1,p_DrhsV1,p_DrhsP1,&
                                     p_DvecU1,p_DvecV1,p_DvecP1
  real(DP), dimension(:), pointer :: p_DrhsU2,p_DrhsV2,p_DrhsP2,&
                                     p_DvecU2,p_DvecV2,p_DvecP2
  
  ! local variables
  logical :: bHaveA12, bHaveA45, bhaveA14,bhaveA15,bhaveA41,bhaveA51,bHaveC
  integer :: idxu,idxp,idxp2,idofp,idofu,i,j,id1,id2,ndofu,ndofp,info,ielidx
  real(DP) :: daux1,daux2,daux4,daux5
  real(DP) :: dp1,dp2
  real(DP), dimension(:,:), allocatable :: DdefectU
  real(DP), dimension(:,:), allocatable :: DdefectP

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
    
    ! Let's assume we do not have the optional matrices
    bHaveA12 = .FALSE.
    bHaveC = .FALSE.
    
    ! Get the pointers from the vanka structure
    call lsyssc_getbase_Kld (rmatrix%RmatrixBlock(1,1),p_KcolA11)
    call lsyssc_getbase_Kcol (rmatrix%RmatrixBlock(1,1),p_KldA11)
    call lsyssc_getbase_Kdiagonal (rmatrix%RmatrixBlock(1,1),p_KdiagA11)

    call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,1),p_Da11)
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,2),p_Da22)

    bhaveA12 = lsysbl_isSubmatrixPresent(rmatrix,1,2)
    if (bhaveA12) then
      call lsyssc_getbase_Kld (rmatrix%RmatrixBlock(1,2),p_KldA12)
      call lsyssc_getbase_Kcol (rmatrix%RmatrixBlock(1,2),p_KcolA12)

      call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,2),p_Da12)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,1),p_Da21)
    end if

    call lsyssc_getbase_Kld (rmatrix%RmatrixBlock(1,3),p_KldB)
    call lsyssc_getbase_Kcol (rmatrix%RmatrixBlock(1,3),p_KcolB)

    call lsyssc_getbase_Kld (rmatrix%RmatrixBlock(3,1),p_KldD)
    call lsyssc_getbase_Kcol (rmatrix%RmatrixBlock(3,2),p_KcolD)

    call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,3),p_Db1)
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,3),p_Db2)

    call lsyssc_getbase_double (rmatrix%RmatrixBlock(3,1),p_Dd1)
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(3,2),p_Dd2)
    
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(4,4),p_Da11)
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(5,5),p_Da22)

    bhaveA45 = lsysbl_isSubmatrixPresent(rmatrix,4,5)
    if (bhaveA45) then
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(4,5),p_Da12)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(5,4),p_Da21)
    end if
    
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(4,6),p_Db1)
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(5,6),p_Db2)

    call lsyssc_getbase_double (rmatrix%RmatrixBlock(6,4),p_Dd1)
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(6,5),p_Dd2)

    bhaveC = lsysbl_isSubmatrixPresent(rmatrix,3,3)
    if (bhaveC) then
      call lsyssc_getbase_Kld (rmatrix%RmatrixBlock(3,3),p_KldC)
      call lsyssc_getbase_Kcol (rmatrix%RmatrixBlock(3,3),p_KcolC)

      call lsyssc_getbase_double (rmatrix%RmatrixBlock(3,3),p_DC1)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(6,6),p_DC2)
    end if
    
    bhaveA14 = lsysbl_isSubmatrixPresent(rmatrix,1,4)
    if (bhaveA14) then
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,4),p_Da14)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,5),p_Da25)
    end if

    bhaveA15 = lsysbl_isSubmatrixPresent(rmatrix,1,5)
    if (bhaveA14) then
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,5),p_Da15)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,4),p_Da24)
    end if

    bhaveA41 = lsysbl_isSubmatrixPresent(rmatrix,4,1)
    if (bhaveA41) then
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(4,1),p_Da41)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(5,2),p_Da52)
    end if

    bhaveA51 = lsysbl_isSubmatrixPresent(rmatrix,5,1)
    if (bhaveA51) then
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(5,1),p_Da51)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(4,2),p_Da42)
    end if

    ! Get the multiplication factors
    Dmult = rmatrix%RmatrixBlock(:,:)%dscaleFactor

    ! Take care of the 'soft-deactivation' of the sub-matrices
    bHaveA12 = bHaveA12 .and. ((Dmult(1,2) .ne. 0.0_DP) .or. &
                               (Dmult(2,1) .ne. 0.0_DP))
    bHaveA45 = bHaveA45 .and. ((Dmult(4,5) .ne. 0.0_DP) .or. &
                               (Dmult(5,4) .ne. 0.0_DP))
    bHaveC = bHaveC .and. (Dmult(3,3) .ne. 0.0_DP)

    ! Allocate an array for the pressure DOF's.
    ndofp = elem_igetNDofLoc(celementP)
    allocate(IdofsP(ndofp,1))

    ! Allocate memory for the correction on each element.
    ! Note that all P-dofs on one element are connected to all the V-dofs on that
    ! element. Therefore, each row in D corresponding to an arbitrary P-dof on an
    ! element will return all the V-dof's on that element!
    ndofu = elem_igetNDofLoc(celementU)
    allocate (DdefectU(ndofu,4))
    allocate (DdefectP(ndofp,2))
    allocate (Ds1(ndofp,ndofp))
    allocate (Ds2(ndofp,ndofp))
    allocate (Ipiv(ndofp))
    allocate (DaInv(4,ndofu))
    allocate (IdofsU(ndofu))
    
    ! Loop through all elements 
    do ielidx = 1,size(IelementList)
    
      ! On the element, get the local DOF's in the pressure space
      IelIdx2(1) = ielidx
      call dof_locGlobMapping_mult(rmatrix%p_rblockDiscrTest%RspatialDiscr(3), &
          IelIdx2, IdofsP)
      
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
        DaInv(3,idxu) = 1.0_DP/(Dmult(4,4)*p_DA11(i))
        DaInv(4,idxu) = 1.0_DP/(Dmult(5,5)*p_DA22(i))
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
        
      ! Create: f_u = f_u - B p - A u
      ! with A being the diagonal velocity matrix.  
      do idxp = 1, ndofp
              
        idofP = IdofsP(idxp,1)
              
        ! Now let's loop over the entries of row idofp in the D-matrices
        do id1 = p_KldD(idofp), p_KldD(idofp+1)-1
        
          ! The column index gives us the index of a velocity DOF which is
          ! adjacent to the current pressure dof - so get its index.
          idofu = p_KcolD(id1)
          idxu = id1-p_KldD(idofp)+1
          
          ! The first thing we want to do is to perform:
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
          
          ! Now we'll also subtract A*u from the local RHS
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
          
        end do
        
      end do
          
      ! Do the A12/A21 matrices exist? If yes, then we will also need to
      ! update the local defect by these matrices.
      if(bHaveA12) then
      
        do idxp = 1, ndofp
        
          idofp = IdofsP(idxp,1)
                
          ! Now let's loop over the entries of row idofp in the D-matrices
          do id1 = p_KldD(idofp), p_KldD(idofp+1)-1

            idofu = p_KcolD(id1)
            idxu = id1-p_KldD(idofp)+1
          
            ! f_u := f_u - A12*v
            ! f_v := f_v - A21*u
            daux1 = 0.0_DP
            daux2 = 0.0_DP
            do i = p_KldA12(idofu), p_KldA12(idofu+1)-1
              j = p_KcolA12(i)
              daux1 = daux1 + p_DA12(i)*p_DvecV1(j)
              daux2 = daux2 + p_DA21(i)*p_DvecU1(j)
            end do
            DdefectU(idofu,1) = DdefectU(idofu,1) - Dmult(1,2)*daux1
            DdefectU(idofu,2) = DdefectU(idofu,2) - Dmult(2,1)*daux2
            
          end do
        end do
              
      end if

      ! Do the A45/A54 matrices exist? If yes, then we will also need to
      ! update the local defect by these matrices.
      if(bHaveA45) then
      
        do idxp = 1, ndofp
        
          idofp = IdofsP(idxp,1)
                
          ! Now let's loop over the entries of row idofp in the D-matrices
          do id1 = p_KldD(idofp), p_KldD(idofp+1)-1

            idofu = p_KcolD(id1)
            idxu = id1-p_KldD(idofp)+1
          
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
        
        end do
        
      end if

      ! Do the A14/A25 matrices exist? If yes, then we will also need to
      ! update the local defect by these matrices.
      if(bHaveA14) then

        do idxp = 1, ndofp
        
          idofp = IdofsP(idxp,1)
                
          ! Now let's loop over the entries of row idofp in the D-matrices
          do id1 = p_KldD(idofp), p_KldD(idofp+1)-1

            idofu = p_KcolD(id1)
            idxu = id1-p_KldD(idofp)+1
          
            daux4 = 0.0_DP
            daux5 = 0.0_DP
            do i = p_KldA11(idofu), p_KldA11(idofu+1)-1
              j = p_KcolA11(i)
              daux4 = daux4 + p_DA14(i)*p_DvecU2(j)
              daux5 = daux5 + p_DA25(i)*p_DvecV2(j)
            end do
            
            DdefectU(idxu,1) = DdefectU(idxu,1) - Dmult(4,5)*daux4
            DdefectU(idxu,2) = DdefectU(idxu,2) - Dmult(5,4)*daux5
            
          end do
          
        end do
      
      end if

      ! Do the A15/A24 matrices exist? If yes, then we will also need to
      ! update the local defect by these matrices.
      if(bHaveA15) then

        do idxp = 1, ndofp
        
          idofp = IdofsP(idxp,1)
                
          ! Now let's loop over the entries of row idofp in the D-matrices
          do id1 = p_KldD(idofp), p_KldD(idofp+1)-1

            idofu = p_KcolD(id1)
            idxu = id1-p_KldD(idofp)+1

            daux4 = 0.0_DP
            daux5 = 0.0_DP
            do i = p_KldA12(idofu), p_KldA12(idofu+1)-1
              j = p_KcolA12(i)
              daux4 = daux4 + p_DA15(i)*p_DvecV2(j)
              daux5 = daux5 + p_DA24(i)*p_DvecU2(j)
            end do
            
            DdefectU(idxu,1) = DdefectU(idxu,1) - Dmult(1,5)*daux4
            DdefectU(idxu,2) = DdefectU(idxu,2) - Dmult(2,4)*daux5

          end do
          
        end do
        
      end if

      ! Do the A41/A52 matrices exist? If yes, then we will also need to
      ! update the local defect by these matrices.
      if(bHaveA41) then

        do idxp = 1, ndofp
        
          idofp = IdofsP(idxp,1)
                
          ! Now let's loop over the entries of row idofp in the D-matrices
          do id1 = p_KldD(idofp), p_KldD(idofp+1)-1

            idofu = p_KcolD(id1)
            idxu = id1-p_KldD(idofp)+1

            daux4 = 0.0_DP
            daux5 = 0.0_DP
            do i = p_KldA11(idofu), p_KldA11(idofu+1)-1
              j = p_KcolA11(i)
              daux4 = daux4 + p_DA41(i)*p_DvecU1(j)
              daux5 = daux5 + p_DA52(i)*p_DvecV1(j)
            end do
            
            DdefectU(idxu,3) = DdefectU(idxu,3) - Dmult(4,1)*daux4
            DdefectU(idxu,4) = DdefectU(idxu,4) - Dmult(5,2)*daux5

          end do
          
        end do
        
      end if
          

      ! Do the A42/A51 matrices exist? If yes, then we will also need to
      ! update the local defect by these matrices.
      if(bHaveA51) then

        do idxp = 1, ndofp
        
          idofp = IdofsP(idxp,1)
                
          ! Now let's loop over the entries of row idofp in the D-matrices
          do id1 = p_KldD(idofp), p_KldD(idofp+1)-1

            idofu = p_KcolD(id1)
            idxu = id1-p_KldD(idofp)+1

            daux4 = 0.0_DP
            daux5 = 0.0_DP
            do i = p_KldA12(idofu), p_KldA12(idofu+1)-1
              j = p_KcolA12(i)
              daux4 = daux4 + p_DA42(i)*p_DvecV1(j)
              daux5 = daux5 + p_DA51(i)*p_DvecU1(j)
            end do
            
            DdefectU(idxu,3) = DdefectU(idxu,3) - Dmult(4,2)*daux4
            DdefectU(idxu,4) = DdefectU(idxu,4) - Dmult(5,1)*daux5

          end do
          
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
      do idxp = 1, ndofp
        idofp = IdofsP(idxp,1)
        do id1 = p_KldD(idofp), p_KldD(idofp+1)-1
          idofu = p_KcolD(id1)
          
          do id2 = p_KldB(idofu), p_KldB(idofu+1)-1
            idxp2 = id2-p_KldB(idofu)+1
            Ds1(idxp,idxp2) = -p_Dd1(id1)*p_Db1(id2)*DaInv(1,idxu)&
                              -p_Dd2(id1)*p_Db2(id2)*DaInv(2,idxu)
            Ds2(idxp,idxp2) = -p_Dd4(id1)*p_Db5(id2)*DaInv(1,idxu)&
                              -p_Dd5(id1)*p_Db5(id2)*DaInv(2,idxu)
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
      
      ! Get the update for u:
      ! du = A^-1 ( d1 - B dp )
      
      do i=1,ndofu
        idofu = IdofsU (i)
        do id1 = p_KldB(idofu),p_KldB(idofu+1)-1
          idxp = id1 - p_KldB(idofu) + 1
          DdefectU(i,1) = DdefectU(i,1) - p_Db1(id1)*DdefectP(idxp,1)
          DdefectU(i,2) = DdefectU(i,2) - p_Db2(id1)*DdefectP(idxp,2)
          DdefectU(i,3) = DdefectU(i,3) - p_Db4(id1)*DdefectP(idxp,3)
          DdefectU(i,4) = DdefectU(i,4) - p_Db5(id1)*DdefectP(idxp,4)
        end do
        DdefectU(i,1) = DaInv(1,i)*DdefectU(i,1)
        DdefectU(i,2) = DaInv(2,i)*DdefectU(i,2)
        DdefectU(i,3) = DaInv(3,i)*DdefectU(i,3)
        DdefectU(i,4) = DaInv(4,i)*DdefectU(i,4)
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

    ! That's it

  end subroutine

  ! ***************************************************************************

  !<subroutine>
  
    subroutine cc_initNonlinMatrix (rnonlinearSpatialMatrix,rproblem,&
        rdiscretisation,rstaticLevelInfo)
  
  !<description>
    ! Initialises the rnonlinearCCMatrix structure with parameters and pointers
    ! from the main problem and precalculated information.
  !</description>

  !<input>
    ! Global problem structure.
    type(t_problem), intent(inout) :: rproblem
    
    ! Discretisation of the level where the matrix is to be assembled.
    type(t_blockDiscretisation), intent(in), target :: rdiscretisation
    
    ! Core equation structure of one level.
    type(t_staticLevelInfo), intent(in), target :: rstaticLevelInfo
  !</input>
  
  !<inputoutput>
    ! Nonlinear matrix structure.
    ! Basic parameters in this structure are filled with data.
    type(t_nonlinearSpatialMatrix), intent(inout) :: rnonlinearSpatialMatrix
  !</inputoutput>
               
  !</subroutine>
      
      ! Initialise the matrix assembly structure rnonlinearCCMatrix 
      ! with basic global information.
      !
      ! 1.) Model, stabilisation
      rnonlinearSpatialMatrix%dnu = collct_getvalue_real (rproblem%rcollection,'NU')

      ! Get stabilisation parameters
      call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                                'iUpwind1',rnonlinearSpatialMatrix%iupwind1,0)
      call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                                'iUpwind2',rnonlinearSpatialMatrix%iupwind2,0)
      
      call parlst_getvalue_double (rproblem%rparamList,'CC-DISCRETISATION',&
                                  'dUpsam1',rnonlinearSpatialMatrix%dupsam1,0.0_DP)
      call parlst_getvalue_double (rproblem%rparamList,'CC-DISCRETISATION',&
                                  'dUpsam2',rnonlinearSpatialMatrix%dupsam2,0.0_DP)

      ! Change the sign of dupsam2 for a consistent stabilisation.
      ! Reason: The stablisation is added to the dual operator by the SD/
      ! EOJ stabilisation in the following way:
      !
      !    ... - (u grad lamda + dupsam2*stabilisation) + ... = rhs
      !
      ! We want to *add* the stabilisation, so we have to introduce a "-" sign
      ! in dupsam2 to get
      !
      !    ... - (u grad lamda) - (-dupsam2*stabilisation) + ... = rhs
      ! <=>
      !    ... - (u grad lamda) + dupsam2*stabilisation + ... = rhs
      
      rnonlinearSpatialMatrix%dupsam2 = -rnonlinearSpatialMatrix%dupsam2
      
      ! 2.) Pointers to global precalculated matrices.
      rnonlinearSpatialMatrix%p_rdiscretisation => rdiscretisation
      rnonlinearSpatialMatrix%p_rstaticInfo => rstaticLevelInfo
      
    end subroutine

  ! ***************************************************************************

  !<subroutine>
  
    subroutine cc_preparePrecondMatrixAssembly (rnonlinearSpatialMatrix,&
        ilev,nlmin,nlmax,rprecSpecials)
  
  !<description>
    ! Prepares a rnonlinearCCMatrix structure for the assembly according
    ! to a preconditioner. rprecSpecials specifies a couple of preconditioner
    ! flags that configure the shape of the system matrix that the preconditioner
    ! needs.
    !
    ! cc_initNonlinMatrix must have been called prior to this routine to
    ! initialise the basic matrix. cc_preparePrecondMatrixAssembly will then
    ! add assembly-specific parameters of the preconditioner.
  !</description>

  !<input>
    ! Current assembly level.
    integer, intent(in) :: ilev
    
    ! Minimum assembly level.
    integer, intent(in) :: nlmin
    
    ! Maximum assembly level.
    integer, intent(in) :: nlmax
  
    ! Structure with assembly-specific parameters of the preconditioner.
    type(t_ccPreconditionerSpecials), intent(in) :: rprecSpecials
  !</input>
  
  !<inputoutput>
    ! Nonlinear matrix structure.
    ! Basic parameters in this structure are filled with data.
    type(t_nonlinearSpatialMatrix), intent(inout) :: rnonlinearSpatialMatrix
  !</inputoutput>
               
  !</subroutine>
      
      ! Parameters for adaptive matrices for Q1~ with anisotropic elements
      rnonlinearSpatialMatrix%iadaptiveMatrices = rprecSpecials%iadaptiveMatrices
      rnonlinearSpatialMatrix%dadmatthreshold = rprecSpecials%dadmatthreshold
      
      ! Depending on the level, we have to set up information about
      ! transposing B-matrices.
      if (ilev .eq. nlmin) then
        rnonlinearSpatialMatrix%bvirtualTransposedD = rprecSpecials%bneedVirtTransposedDonCoarse
      else
        rnonlinearSpatialMatrix%bvirtualTransposedD = rprecSpecials%bneedVirtTransposedD
      end if
      
    end subroutine

  ! ***************************************************************************
  ! Routines to create a nonlinear iteration structure, to save it
  ! to a collection, to rebuild it from there and to clean it up.
  ! ***************************************************************************

!<subroutine>

  subroutine cc_createSpacePreconditioner (rpreconditioner,NLMIN,NLMAX)
  
!<description>
  ! This routine creates a spational preconditioner structure. The structure is
  ! initialised to handle NLMAX-NLMIN+1 discretisation levels.
!</description>

!<input>
  ! Minimum discretisation level to be maintained
  integer, intent(IN) :: NLMIN
  
  ! Maximum discretisation level to be maintained. The maximum level coincides
  ! with the level where to solve the system.
  integer, intent(IN) :: NLMAX
!</input>

!<output>
  ! A spatial preconditioner structure to be initialised.
  type(t_ccspatialPreconditioner), intent(OUT) :: rpreconditioner
!</output>

!</subroutine>

    rpreconditioner%NLMIN = NLMIN
    rpreconditioner%NLMAX = NLMAX

    ! Initialise the matrix pointers on all levels that we have to maintain.
    allocate(rpreconditioner%RcoreEquation(NLMIN:NLMAX))

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_releaseSpacePreconditioner (rpreconditioner)
  
!<description>
  ! Releases allocated memory in the spatial preconditioner structure.
!</description>

!<inputoutput>
  ! The spatial preconditioner structure that should be cleaned up.
  type(t_ccspatialPreconditioner), intent(INOUT) :: rpreconditioner
!</inputoutput>

!</subroutine>
    
    if (associated(rpreconditioner%RcoreEquation)) &
      deallocate(rpreconditioner%RcoreEquation)

    rpreconditioner%NLMIN = 0
    rpreconditioner%NLMAX = 0

  end subroutine

  ! ***************************************************************************

  !<subroutine>

    subroutine cc_precondSpaceDefect (rpreconditioner,rnonlinearSpatialMatrix,&
        rd,rx1,rx2,rx3,bsuccess,rcollection)
  
    use linearsystemblock
    use collection
    
  !<description>
    ! Defect preconditioning routine. Based on the current iteration 
    ! vector rx, this routine has to perform preconditioning on the defect 
    ! vector rd.
  !</description>

  !<input>
    ! Configuration of the core equation on the maximum level.
    type(t_nonlinearSpatialMatrix), intent(IN)      :: rnonlinearSpatialMatrix

  !</input>

  !<inputoutput>
    ! Spatial preconditioner structure that defines all parameters how to perform
    ! preconditioning.
    type(t_ccspatialPreconditioner), intent(INOUT) :: rpreconditioner

    ! Defect vector b-A(rx)x. This must be replaced by J^{-1} rd by a preconditioner.
    type(t_vectorBlock), intent(INOUT)            :: rd

    ! Ccollection structure of the application.
    type(t_collection)                            :: rcollection
    
    ! If the preconditioning was a success. Is normally automatically set to
    ! TRUE. If there is an error in the preconditioner, this flag can be
    ! set to FALSE. In this case, the nonlinear solver breaks down with
    ! the error flag set to 'preconditioner broke down'.
    logical, intent(INOUT)                        :: bsuccess
  !</inputoutput>
  
  !<input>
    ! Iteration vector of the 'previous' timestep. Can be undefined if there
    ! is no previous timestep.
    type(t_vectorBlock), intent(IN), target       :: rx1

    ! Iteration vector of the 'current' timestep.
    type(t_vectorBlock), intent(IN), target       :: rx2

    ! Iteration vector of the 'next' timestep. Can be undefined if there
    ! is no next timestep.
    type(t_vectorBlock), intent(IN), target       :: rx3
  !</input>
  
  !</subroutine>
  
    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    integer :: ierror
    type(t_linsolNode), pointer :: p_rsolverNode
    
    integer :: i
    logical :: bassembleNewton
    type(t_matrixBlock), dimension(:), pointer :: Rmatrices
    type(t_filterChain), dimension(:), pointer :: p_RfilterChain
    

    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Ddata
    type(t_vectorBlock) :: rxtemp
    integer, dimension(:), pointer :: p_IelementList

    ! DEBUG!!!
    real(dp), dimension(:), pointer :: p_vec,p_def,p_da
    call lsysbl_getbase_double (rd,p_def)
    call lsysbl_getbase_double (rx2,p_vec)

      select case (rpreconditioner%ctypePreconditioning)
      case (CCPREC_NONE)
        ! No preconditioning. Do nothing.
      case (CCPREC_LINEARSOLVER,CCPREC_NEWTON,CCPREC_INEXACTNEWTON)
        ! Preconditioning with a linear solver.
        !
        ! At first, assemble the preconditioner matrices on all levels
        ! and incorporate all boundary conditions.
        !
        ! Should we assemble the Newton part?
        
        bassembleNewton = .false.
        
        if ((rpreconditioner%ctypePreconditioning .eq. CCPREC_NEWTON) .or. &
            (rpreconditioner%ctypePreconditioning .eq. CCPREC_INEXACTNEWTON)) then
            
          ! Use Newton in any case.
          bassembleNewton = .true.
          
        end if
        
        ! Assemble the preconditioner matrices in rpreconditioner
        ! on all levels that the solver uses.
        call assembleLinsolMatrices (rpreconditioner,rnonlinearSpatialMatrix,&
            rcollection,bassembleNewton,rx1,rx2,rx3)
          
        ! Our 'parent' (the caller of the nonlinear solver) has prepared
        ! a preconditioner node for us (a linear solver with symbolically
        ! factorised matrices). Get this from the collection.
      
        p_rsolverNode => rpreconditioner%p_rsolverNode

        ! Re-attach the system matrices to the solver.
        ! Note that no pointers and no handles are changed, so we can savely do
        ! that without calling linsol_doneStructure/linsol_doneStructure.
        ! This simply informs the solver about possible new scaling factors
        ! in the matrices in case they have changed...
        allocate(Rmatrices(rpreconditioner%NLMIN:rpreconditioner%NLMAX))
        
        ! Attach the system matrix
        do i=rpreconditioner%NLMIN,rpreconditioner%NLMAX
          call lsysbl_duplicateMatrix ( &
            rpreconditioner%RcoreEquation(i)%p_rmatrixPreconditioner, &
            Rmatrices(i), LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
        end do
          
        ! DEBUG!!!
        !CALL matio_writeBlockMatrixHR (Rmatrices(rpreconditioner%NLMIN), 'matrix',&
        !                               .TRUE., 0, 'matrixstat.txt','(E10.2)')
        
        call linsol_setMatrices(rpreconditioner%p_rsolverNode,Rmatrices(:))
        
        ! DEBUG!!!
        !DO i=rnonlinearIteration%NLMIN,rnonlinearIteration%NLMAX
        !  CALL storage_getbase_double (Rmatrices(i)% &
        !      RmatrixBlock(4,1)%h_Da,p_Ddata)
        !END DO
        
        ! Initialise data of the solver. This in fact performs a numeric
        ! factorisation of the matrices in UMFPACK-like solvers.
        call linsol_updateStructure (rpreconditioner%p_rsolverNode,ierror)
        call linsol_initData (p_rsolverNode, ierror)
        if (ierror .ne. LINSOL_ERR_NOERROR) then
          print *,'linsol_initData failed!'
          call sys_halt()
        end if
        
        ! The solver got the matrices; clean up Rmatrices, it was only of temporary
        ! nature...
        do i=rpreconditioner%NLMIN,rpreconditioner%NLMAX
          call lsysbl_releaseMatrix (Rmatrices(i))
        end do
        deallocate(Rmatrices)

        ! Solve the system. As we want to solve Ax=b with
        ! b being the real RHS and x being the real solution vector,
        ! we use linsol_solveAdaptively. If b is a defect
        ! RHS and x a defect update to be added to a solution vector,
        ! we would have to use linsol_precondDefect instead.
        call linsol_precondDefect (p_rsolverNode,rd)
!        call lsysbl_createVecBlockIndirect(rd,rxtemp,.true.)
!        call storage_getbase_int(rd%p_rblockDiscr%RspatialDiscr(1)%RelementDistr(1)%h_IelementList,&
!          p_IelementList)
!        do i=1,1
!          call vanka_NSOptC2D (rd,1.0_DP,rxtemp,Rmatrices(rpreconditioner%NLMAX),&
!              p_IelementList,EL_EM30,EL_Q0)
!        end do
!        call lsysbl_copyVector (rxtemp,rd)
!        call lsysbl_releaseVector (rxtemp)

        ! Release the numeric factorisation of the matrix.
        ! We don't release the symbolic factorisation, as we can use them
        ! for the next iteration.
        call linsol_doneData (p_rsolverNode)
        
        ! Did the preconditioner work?
        bsuccess = p_rsolverNode%iresult .eq. 0
        
        if (bsuccess) then
          ! Filter the final defect
          p_RfilterChain => rpreconditioner%p_RfilterChain
          call filter_applyFilterChainVec (rd, p_RfilterChain)
        end if
        
        if (p_rsolverNode%dfinalDefect .gt. p_rsolverNode%dinitialDefect*0.99_DP) then
          ! Ignore the correction, it cannot be good enough!
          call output_line (&
            'Space-Time-Preconditioner: Warning. Solution ignored for missing accuracy.')
            
          call lsysbl_clearVector (rd)
        end if
        
      end select
      
    contains
      
      subroutine assembleLinsolMatrices (rpreconditioner,rnonlinearSpatialMatrix,rcollection,&
          bassembleNewton,rx1,rx2,rx3)

      use linearsystemblock
      use collection

      ! Assembles on every level a matrix for the linear-solver/Newton preconditioner.
      ! bnewton allows to specify whether the Newton matrix or only the standard
      ! system matrix is evaluated. The output is written to the p_rpreconditioner 
      ! matrices specified in the rnonlinearIteration structure.

      ! Spatial preconditioner structure that defines all parameters how to perform
      ! preconditioning.
      type(t_ccspatialPreconditioner), intent(IN)    :: rpreconditioner

      ! Level independent configuration of the core equation
      type(t_nonlinearSpatialMatrix), intent(IN)      :: rnonlinearSpatialMatrix

      ! Reference to a collection structure that contains all parameters of the
      ! discretisation (for nonlinearity, etc.).
      type(t_collection), intent(INOUT)                :: rcollection

      ! TRUE  = Assemble the Newton preconditioner.
      ! FALSE = Assemble the standard defect correction preconditioner
      !         (i.e. the linearised system matrix).
      logical, intent(IN) :: bassembleNewton
      
      ! Current iteration vector of the 'previous' timestep. May be undefined
      ! if there is no previous timestep. 
      type(t_vectorBlock), intent(IN), target          :: rx1

      ! Current iteration vector. 
      type(t_vectorBlock), intent(IN), target          :: rx2

      ! Current iteration vector of the 'next' timestep. May be undefined
      ! if there is no previous timestep. 
      type(t_vectorBlock), intent(IN), target          :: rx3

      ! local variables
      real(DP) :: dnewton
      integer :: ilev
      type(t_matrixBlock), pointer :: p_rmatrix,p_rmatrixFine
      type(t_vectorScalar), pointer :: p_rvectorTemp
      type(t_vectorBlock), pointer :: p_rvectorFine1,p_rvectorFine2,p_rvectorFine3
      type(t_vectorBlock), pointer :: p_rvectorCoarse1,p_rvectorCoarse2,p_rvectorCoarse3
      type(t_interlevelProjectionBlock), pointer :: p_rprojection
      type(t_nonlinearSpatialMatrix) :: rlocalNonlSpatialMatrix
      integer, dimension(1), parameter :: Irows = (/1/)

      ! A filter chain for the linear solver
      type(t_filterChain), dimension(:), pointer :: p_RfilterChain
      
      ! DEBUG!!!
      type(t_cccoreEquationOneLevel), pointer :: p_rcore
    !    real(dp), dimension(:), pointer :: p_vec,p_def,p_da
    !    call lsysbl_getbase_double (rd,p_def)
    !    call lsysbl_getbase_double (rx,p_vec)

        ! Get the interlevel projection structure and the temporary vector
        ! from the collection.
        ! Our 'parent' prepared there how to interpolate the solution on the
        ! fine grid to coarser grids.
        p_rvectorTemp => rpreconditioner%p_rtempVectorSc

        ! Get the filter chain. We need that later to filter the matrices.        
        p_RfilterChain => rpreconditioner%p_RfilterChain

        ! Initialise the matrix assembly structure rlocalNonlSpatialMatrix to describe the
        ! matrix we want to have. We have to initialise the adaptivity constants,
        ! which are not part of the standard initialisation.
        rlocalNonlSpatialMatrix = rnonlinearSpatialMatrix
        rlocalNonlSpatialMatrix%iadaptiveMatrices = &
            rpreconditioner%rprecSpecials%iadaptiveMatrices
        rlocalNonlSpatialMatrix%dadmatthreshold = &
            rpreconditioner%rprecSpecials%dadmatthreshold

        ! On all levels, we have to set up the nonlinear system matrix,
        ! so that the linear solver can be applied to it.
        
        nullify(p_rmatrix)

        do ilev=rpreconditioner%NLMAX,rpreconditioner%NLMIN,-1
        
          ! Get the matrix on the current level.
          ! Shift the previous matrix to the pointer of the fine grid matrix.
          p_rmatrixFine => p_rmatrix
          p_rmatrix => rpreconditioner%RcoreEquation(ilev)%p_rmatrixPreconditioner
        
          ! On the highest level, we use rx as solution to build the nonlinear
          ! matrix. On lower levels, we have to create a solution
          ! on that level from a fine-grid solution before we can use
          ! it to build the matrix!
          if (ilev .eq. rpreconditioner%NLMAX) then
          
            p_rvectorCoarse1 => rx1
            p_rvectorCoarse2 => rx2
            p_rvectorCoarse3 => rx3
            
          else
            ! We have to discretise a level hierarchy and are on a level < NLMAX.

            ! Get the mujltilevel projection structure that describes the
            ! projection from the finer level to the current one.
            p_rprojection => rpreconditioner%p_Rprojection(ilev+1)
            
            ! Get the temporary vector on level i. Will receive the solution
            ! vector on that level. 
            p_rvectorCoarse1 => rpreconditioner%RcoreEquation(ilev)%p_rtempVector1
            p_rvectorCoarse2 => rpreconditioner%RcoreEquation(ilev)%p_rtempVector2
            p_rvectorCoarse3 => rpreconditioner%RcoreEquation(ilev)%p_rtempVector3
            
            ! Get the solution vector on level i+1. This is either the temporary
            ! vector on that level, or the solution vector on the maximum level.
            if (ilev .lt. rpreconditioner%NLMAX-1) then
              p_rvectorFine1 => rpreconditioner%RcoreEquation(ilev+1)%p_rtempVector1
              p_rvectorFine2 => rpreconditioner%RcoreEquation(ilev+1)%p_rtempVector2
              p_rvectorFine3 => rpreconditioner%RcoreEquation(ilev+1)%p_rtempVector3
            else
              p_rvectorFine1 => rx1
              p_rvectorFine2 => rx2
              p_rvectorFine3 => rx3
            end if

            ! Interpolate the solution from the finer grid to the coarser grid.
            ! The interpolation is configured in the interlevel projection
            ! structure we got from the collection.
            call mlprj_performInterpolation (p_rprojection,p_rvectorCoarse1, &
                                             p_rvectorFine1,p_rvectorTemp)
            call mlprj_performInterpolation (p_rprojection,p_rvectorCoarse2, &
                                             p_rvectorFine2,p_rvectorTemp)
            call mlprj_performInterpolation (p_rprojection,p_rvectorCoarse3, &
                                             p_rvectorFine3,p_rvectorTemp)

            ! Apply the filter chain to the temp vector.
            ! This implements the boundary conditions that are attached to it.
            ! NOTE: Deactivated for standard CC2D compatibility -- and because
            ! it has to be checked whether the correct boundary conditions
            ! are attached to that vector!
            ! CALL filter_applyFilterChainVec (p_rvectorCoarse, p_RfilterChain)

          end if

          ! Set the pointers in the rlocalNonlSpatialMatrix structure according
          ! to the current level.
          p_rcore => rpreconditioner%RcoreEquation(ilev)

          rlocalNonlSpatialMatrix%p_rdiscretisation         => &
              p_rmatrix%p_rblockDiscrTrial

          rlocalNonlSpatialMatrix%p_rstaticInfo         => &
              rpreconditioner%RcoreEquation(ilev)%p_rstaticInfo

          ! Assemble the matrix.
          ! If we are on a lower level, we can specify a 'fine-grid' matrix.
          if (ilev .eq. rpreconditioner%NLMAX) then
            call cc_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
                p_rmatrix,rlocalNonlSpatialMatrix,&
                p_rvectorCoarse1,p_rvectorCoarse2,p_rvectorCoarse3)
          else
            call cc_assembleMatrix (CCMASM_COMPUTE,CCMASM_MTP_AUTOMATIC,&
                p_rmatrix,rlocalNonlSpatialMatrix,&
                p_rvectorCoarse1,p_rvectorCoarse2,p_rvectorCoarse3,&
                p_rmatrixFine)
          end if

          ! Boundary conditions
          ! ---------------------------------------------------

          if (associated(p_RfilterChain)) then
            ! Apply the filter chain to the matrix.
            ! As the filter consists only of an implementation filter for
            ! boundary conditions, this implements the boundary conditions
            ! into the system matrix.
            call filter_applyFilterChainMat (p_rmatrix, p_RfilterChain)
          else
            ! Call the matrix filter for the boundary conditions to include the BC's
            ! into the matrix.
            call matfil_discreteBC (p_rmatrix)
            call matfil_discreteFBC (p_rmatrix)
          end if
            
          ! 'Nonlinear' boundary conditions like slip boundary conditions
          ! are not implemented with a filter chain into a matrix.
          ! Call the appropriate matrix filter of 'nonlinear' boundary
          ! conditions manually:
          call matfil_discreteNLSlipBC (p_rmatrix,.true.)
            
          ! DEBUG!!!
          !CALL matio_writeBlockMatrixHR (p_rmatrix, 'matrix',&
          !                              .TRUE., 0, 'matrix.txt','(E20.10)')

        end do
        
        if (rpreconditioner%rprecSpecials%bneedPressureDiagonalBlock) then
          
          ! The 3,3-matrix must exist! This is ensured by the initialisation routine.
          !
          ! We have a pure Dirichlet problem. This may give us some difficulties
          ! in the case, the preconditioner uses a direct solver (UMFPACK).
          ! In this case, we have to include a unit vector to the pressure
          ! matrix to make the problem definite!
          if (rpreconditioner%rprecSpecials%isolverType .eq. 0) then
            p_rmatrix => rpreconditioner%RcoreEquation(rpreconditioner%NLMAX)%&
                p_rmatrixPreconditioner
            
            ! Include a unit vector to the matrix part of the pressure in
            ! the primal equation -- as long as there is not a full identity
            ! matrix in the pressure matrix (what would be the case for 
            ! the initial condition).
            if (rlocalNonlSpatialMatrix%Dkappa(1,1) .eq. 0.0_DP) then
              ! Switch the pressure matrix on and clear it; we don't know what is inside.
              p_rmatrix%RmatrixBlock(3,3)%dscaleFactor = 1.0_DP
              call lsyssc_clearMatrix (p_rmatrix%RmatrixBlock(3,3))
              call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,1),Irows)
              call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,2),Irows)
              call mmod_replaceLinesByUnit(p_rmatrix%RmatrixBlock(3,3),Irows)
              call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,4),Irows)
              call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,5),Irows)
              call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,6),Irows)
            end if

            ! Also in the dual equation, as the BC type coincides
            if (rlocalNonlSpatialMatrix%Dkappa(2,2) .eq. 0.0_DP) then
              ! Switch the pressure matrix on and clear it; we don't know what is inside.
              p_rmatrix%RmatrixBlock(6,6)%dscaleFactor = 1.0_DP
              call lsyssc_clearMatrix (p_rmatrix%RmatrixBlock(6,6))
              call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(6,1),Irows)
              call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(6,2),Irows)
              call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(6,3),Irows)
              call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(6,4),Irows)
              call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(6,5),Irows)
              call mmod_replaceLinesByUnit(p_rmatrix%RmatrixBlock(6,6),Irows)
            end if
            
          end if
          
          if (rpreconditioner%rprecSpecials%isolverType .eq. 1) then
          
            ! If we have a MG solver, We also check the coarse grid solver for 
            ! the same thing!
            ! What we don't check is the smoother, thus we assume that smoothers
            ! are always solvers that allow the applicance of a filter chain.
            if (rpreconditioner%rprecSpecials%icoarseGridSolverType .eq. 0) then
              p_rmatrix => rpreconditioner%RcoreEquation(rpreconditioner%NLMIN)%&
                  p_rmatrixPreconditioner
              
              ! Include a unit vector to the matrix part of the pressure in
              ! the primal equation -- as long as there is not a full identity
              ! matrix in the pressure matrix (what would be the case for 
              ! the initial condition).
              if (rlocalNonlSpatialMatrix%Dkappa(1,1) .eq. 0.0_DP) then
                ! Switch the pressure matrix on and clear it; we don't know what is inside.
                p_rmatrix%RmatrixBlock(3,3)%dscaleFactor = 1.0_DP
                call lsyssc_clearMatrix (p_rmatrix%RmatrixBlock(3,3))
                call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,1),Irows)
                call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,2),Irows)
                call mmod_replaceLinesByUnit(p_rmatrix%RmatrixBlock(3,3),Irows)
                call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,4),Irows)
                call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,5),Irows)
                call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(3,6),Irows)
              end if

              ! Also in the dual equation, as the BC type coincides
              if (rlocalNonlSpatialMatrix%Dkappa(2,2) .eq. 0.0_DP) then
                ! Switch the pressure matrix on and clear it; we don't know what is inside.
                p_rmatrix%RmatrixBlock(6,6)%dscaleFactor = 1.0_DP
                call lsyssc_clearMatrix (p_rmatrix%RmatrixBlock(6,6))
                call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(6,1),Irows)
                call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(6,2),Irows)
                call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(6,3),Irows)
                call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(6,4),Irows)
                call mmod_replaceLinesByZero(p_rmatrix%RmatrixBlock(6,5),Irows)
                call mmod_replaceLinesByUnit(p_rmatrix%RmatrixBlock(6,6),Irows)
              end if
              
            end if
            
          end if
            
        end if        
        
      end subroutine
      
    end subroutine

end module
