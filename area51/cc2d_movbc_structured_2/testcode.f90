module testcode

  use fsystem
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  
  implicit none


contains

  subroutine coeff_Laplace (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial,IdofsTest,rdomainIntSubset, p_rcollection,&
                  Dcoefficients)
    
    use BasicGeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the matrix assembly. It has to compute
    ! the coefficients in front of the terms of the bilinear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
    ! According to the terms in the bilinear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the bilinear form
    ! the corresponding coefficients in front of the terms.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
    
    ! The bilinear form which is currently being evaluated:
    type(t_bilinearForm), intent(IN)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(NDIM2D,npointsPerElement,nelements), intent(IN)  :: Dpoints
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in trial space,Number of elements)
    integer(PREC_DOFIDX), dimension(:,:), intent(IN) :: IdofsTrial
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,Number of elements)
    integer(PREC_DOFIDX), dimension(:,:), intent(IN) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A pointer to a collection structure to provide additional
    ! information to the coefficient routine. May point to NULL() if not defined.
    type(t_collection), pointer                      :: p_rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

  Dcoefficients = 1.0_DP

  end subroutine

!***********************************

  !<subroutine>

    subroutine coeff_RHS (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,p_rcollection, &
                  Dcoefficients)
    
    use BasicGeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in in real coordinates.
    ! According to the terms in the linear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the linear form
    ! the corresponding coefficients in front of the terms.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(IN)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    real(DP), dimension(NDIM2D,npointsPerElement,nelements), intent(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,Number of elements)
    integer(PREC_DOFIDX), dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! A pointer to a collection structure to provide additional
    ! information to the coefficient routine. May point to NULL() if not defined.
    type(t_collection), pointer                      :: p_rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    Dcoefficients = 1.0_DP
  
    end subroutine

  !<subroutine>

    subroutine getBoundaryValues (rdiscretisation,rbcRegion,ielement, &
                                   cinfoNeeded,iwhere,dwhere, p_rcollection, Dvalues)
    
    use collection
    use spatialdiscretisation
    use discretebc
    
  !<description>
    ! This subroutine is called during the discretisation of boundary
    ! conditions. It calculates a special quantity on the boundary, which is
    ! then used by the discretisation routines to generate a discrete
    ! 'snapshot' of the (actually analytic) boundary conditions.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN)                   :: rdiscretisation
    
    ! Boundary condition region that is currently being processed.
    ! (This e.g. defines the type of boundary conditions that are
    !  currently being calculated, as well as information about the current
    !  boundary segment 'where we are at the moment'.)
    type(t_bcRegion), intent(IN)                                :: rbcRegion
    
    
    ! The element number on the boundary which is currently being processed
    integer(I32), intent(IN)                                    :: ielement
    
    ! The type of information, the routine should calculate. One of the
    ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
    ! to return one or multiple information value in the result array.
    integer, intent(IN)                                         :: cinfoNeeded
    
    ! A reference to a geometric object where information should be computed.
    ! cinfoNeeded=DISCBC_NEEDFUNC :
    !   iwhere = number of the point in the triangulation or
    !          = 0, if only the parameter value of the point is known; this
    !               can be found in dwhere,
    ! cinfoNeeded=DISCBC_NEEDDERIV :
    !   iwhere = number of the point in the triangulation or
    !          = 0, if only the parameter value of the point is known; this
    !               can be found in dwhere,
    ! cinfoNeeded=DISCBC_NEEDINTMEAN :
    !   iwhere = number of the edge where the value integral mean value
    !            should be computed
    integer, intent(IN)                                         :: iwhere

    ! A reference to a geometric object where information should be computed.
    ! cinfoNeeded=DISCBC_NEEDFUNC :
    !   dwhere = parameter value of the point where the value should be computed,
    ! cinfoNeeded=DISCBC_NEEDDERIV :
    !   dwhere = parameter value of the point where the value should be computed,
    ! cinfoNeeded=DISCBC_NEEDINTMEAN :
    !   dwhere = 0 (not used)
    real(DP), intent(IN)                                        :: dwhere
     
    ! A pointer to a collection structure to provide additional
    ! information to the coefficient routine. May point to NULL() if not defined.
    type(t_collection), pointer                  :: p_rcollection

  !</input>
  
  !<output>
    ! This array receives the calculated information. If the caller
    ! only needs one value, the computed quantity is put into Dvalues(1).
    ! If multiple values are needed, they are collected here (e.g. for
    ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
    real(DP), dimension(:), intent(OUT)                         :: Dvalues
  !</output>
    
  !</subroutine>
  
    Dvalues(1) = 0.0_DP
  
    end subroutine

!***************************************

  subroutine test_MatrixQ1
  
    use collection
    use triangulation
    use spatialdiscretisation

    implicit none

    !INCLUDE 'stria.inc'
    include 'cout.inc'
    include 'cerr.inc'
    include 'cmem.inc'
    include 'cparametrization.inc'

    ! Variables
    
    type(t_boundary), pointer :: p_rboundary
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_spatialDiscretisation), pointer :: p_rdiscretisation
    real(DP) :: x,y,t
    integer, dimension(SZTRIA,20) :: TRIAS
    character(LEN=60) :: CFILE
    real(DP), dimension(:), pointer :: p_Ddata
    
    type(t_matrixScalar) :: rmatrix
    type(t_vectorScalar) :: rvector
    
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinform
    
    integer :: lv, LCOL,LLD,NA,NEQ,LA,LB
    external E011,COEAGS,EM30,FEATRHS,NDFGX
    double precision FEATRHS
    integer :: NDFGX
    integer, dimension(2,2) :: KAB
    real(DP) :: TTT0,TTT1,TTT2

    call storage_init(999, 100)
    call boundary_read_prm(p_rboundary, 'pre/QUAD.prm')
    
    lv = 10
    M = 0
    ICHECK = 0
    IMESH = 1

    CFILE = 'pre/QUAD.prm'
    call GENPAR (.true.,IMESH,CFILE)

    CFILE = 'pre/QUAD.tri'
    call INMTRI (2,TRIAS,lv,lv,0,CFILE)
    call tria_wrp_tria2Structure(TRIAS(:,lv),p_rtriangulation)
    
    call spdiscr_initDiscr_simple (p_rdiscretisation,EL_E011,CUB_TRZ,&
                                   p_rtriangulation, p_rboundary, NULL())
    call ZTIME(TTT0)
    call XAP7X(LCOL,LLD,NA,NEQ,TRIAS(:,lv),E011,0)
    
    call ZTIME(TTT1)
    call bilf_createMatrixStructure (p_rdiscretisation,LSYSSC_MATRIX9,rmatrix)
    
    call ZTIME(TTT2)
    print *,'Time1 structure: ',TTT1-TTT0
    print *,'Time2 structure: ',TTT2-TTT1
    
    call ZTIME(TTT0)
    KAB = reshape((/2,2,3,3/),(/2,2/))
    LA = 0
    call XAB7X(LA,LCOL,LLD,NA,NEQ,1,1,TRIAS(:,lv),E011,.false., &
               COEAGS,.true.,KAB,2,2,0,'MATRIX',0,0.0_DP)
    call ZTIME(TTT1)
    
    !CFILE = 'MATRIX1.TXT'
    !CALL OWM17(DWORK(L(LA)),KWORK(L(LCOL)),KWORK(L(LLD)),&
    !           NEQ,NEQ,.TRUE.,0,'MAT1  ',CFILE,'(D20.10)')
    
    ! set up bilinear form
    rform%itermCount = 2
    rform%Idescriptors(1,1) = DER_DERIV_X
    rform%Idescriptors(2,1) = DER_DERIV_X
    rform%Idescriptors(1,2) = DER_DERIV_Y
    rform%Idescriptors(2,2) = DER_DERIV_Y
    rform%Dcoefficients(1)  = 1.0
    rform%Dcoefficients(2)  = 1.0
    rform%ballCoeffConstant = .true.
    rform%BconstantCoeff = .true.
    !BILF_NELEMSIM = 100
    !BILF_NELEMSIM = 300000
    call bilf_buildMatrixScalar (rform,.true.,rmatrix,coeff_Laplace)
    call ZTIME(TTT2)

    print *,'Time1 assembly: ',TTT1-TTT0
    print *,'Time2 assembly: ',TTT2-TTT1
    
    print *,IWORK,IWMAX,NNWORK
    
    !CALL lsyssc_releaseVector (rvector)
    call lsyssc_releaseMatrix (rmatrix)
    call boundary_release (p_rboundary)
    
  end subroutine

!***************************************

  subroutine test_MatrixEM30
  
    use collection
    use triangulation
    use spatialdiscretisation

    implicit none

    !INCLUDE 'stria.inc'
    include 'cout.inc'
    include 'cerr.inc'
    include 'cmem.inc'
    include 'cparametrization.inc'

    ! Variables
    
    type(t_boundary), pointer :: p_rboundary
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_spatialDiscretisation),pointer :: p_rdiscretisation
    real(DP) :: x,y,t
    integer, dimension(SZTRIA,20) :: TRIAS
    character(LEN=60) :: CFILE
    real(DP), dimension(:), pointer :: p_Ddata
    
    type(t_matrixScalar) :: rmatrix
    type(t_vectorScalar) :: rvector
    
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinform
    
    integer :: lv, LCOL,LLD,NA,NEQ,LA,LB
    external E011,COEAGS,EM30,FEATRHS,NDFGX
    double precision FEATRHS
    integer :: NDFGX
    integer, dimension(2,2) :: KAB
    real(DP) :: TTT0,TTT1,TTT2

    call storage_init(999, 100)
    call boundary_read_prm(p_rboundary, 'pre/QUAD.prm')
    
    lv = 10
    M = 0
    ICHECK = 0
    IMESH = 1

    CFILE = 'pre/QUAD.prm'
    call GENPAR (.true.,IMESH,CFILE)

    CFILE = 'pre/QUAD.tri'
    call INMTRI (2,TRIAS,lv,lv,0,CFILE)
    call tria_wrp_tria2Structure(TRIAS(:,lv),p_rtriangulation)
    
    call spdiscr_initDiscr_simple (p_rdiscretisation,EL_EM30,CUB_TRZ,&
                                   p_rtriangulation, p_rboundary, NULL())
    call ZTIME(TTT0)
    call XAP7X(LCOL,LLD,NA,NEQ,TRIAS(:,lv),EM30,0)
    
    call ZTIME(TTT1)
    call bilf_createMatrixStructure (p_rdiscretisation,LSYSSC_MATRIX9,rmatrix)
    
    call ZTIME(TTT2)
    print *,'Time1 structure: ',TTT1-TTT0
    print *,'Time2 structure: ',TTT2-TTT1
    
    call ZTIME(TTT0)
    KAB = reshape((/2,2,3,3/),(/2,2/))
    LA = 0
    call XAB7X(LA,LCOL,LLD,NA,NEQ,1,1,TRIAS(:,lv),EM30,.true., &
               COEAGS,.true.,KAB,2,2,0,'MATRIX',0,0.0_DP)
    call ZTIME(TTT1)
    
    !CFILE = 'MATRIX1.TXT'
    !CALL OWM17(DWORK(L(LA)),KWORK(L(LCOL)),KWORK(L(LLD)),&
    !           NEQ,NEQ,.TRUE.,0,'MAT1  ',CFILE,'(D20.10)')
    
    ! set up bilinear form
    rform%itermCount = 2
    rform%Idescriptors(1,1) = DER_DERIV_X
    rform%Idescriptors(2,1) = DER_DERIV_X
    rform%Idescriptors(1,2) = DER_DERIV_Y
    rform%Idescriptors(2,2) = DER_DERIV_Y
    rform%Dcoefficients(1)  = 1.0
    rform%Dcoefficients(2)  = 1.0
    rform%ballCoeffConstant = .true.
    rform%BconstantCoeff = .true.
    !BILF_NELEMSIM = 100
    !BILF_NELEMSIM = 300000
    call bilf_buildMatrixScalar (rform,.true.,rmatrix,coeff_Laplace)
    call ZTIME(TTT2)

    print *,'Time1 assembly: ',TTT1-TTT0
    print *,'Time2 assembly: ',TTT2-TTT1
    
    print *,IWORK,IWMAX,NNWORK
    
    !CALL lsyssc_releaseVector (rvector)
    call lsyssc_releaseMatrix (rmatrix)
    call boundary_release (p_rboundary)
    
  end subroutine

!***************************************

  subroutine test_RHSQ1
  
    use collection
    use triangulation
    use spatialdiscretisation

    implicit none

    !INCLUDE 'stria.inc'
    include 'cout.inc'
    include 'cerr.inc'
    include 'cmem.inc'
    include 'cparametrization.inc'

    ! Variables
    
    type(t_boundary), pointer :: p_rboundary
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_spatialDiscretisation), pointer :: p_rdiscretisation
    real(DP) :: x,y,t
    integer, dimension(SZTRIA,20) :: TRIAS
    character(LEN=60) :: CFILE
    real(DP), dimension(:), pointer :: p_Ddata
    
    type(t_matrixScalar) :: rmatrix
    type(t_vectorScalar) :: rvector
    
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinform
    
    integer :: lv, LCOL,LLD,NA,NEQ,LA,LB
    external E011,COEAGS,EM30,FEATRHS,NDFGX
    double precision FEATRHS
    integer :: NDFGX
    integer, dimension(2,2) :: KAB
    real(DP) :: TTT0,TTT1,TTT2

    call storage_init(999, 100)
    call boundary_read_prm(p_rboundary, 'pre/QUAD.prm')
    
    lv = 10
    M = 0
    ICHECK = 0
    IMESH = 1

    CFILE = 'pre/QUAD.prm'
    call GENPAR (.true.,IMESH,CFILE)

    CFILE = 'pre/QUAD.tri'
    call INMTRI (2,TRIAS,lv,lv,0,CFILE)
    call tria_wrp_tria2Structure(TRIAS(:,lv),p_rtriangulation)
    
    call spdiscr_initDiscr_simple (p_rdiscretisation,EL_E011,CUB_TRZ,&
                                   p_rtriangulation, p_rboundary, NULL())
    
    call ZTIME(TTT0)
    
    LB = 0
    NEQ = NDFGX(30,TRIAS(:,lv))
    call XVB0X(TRIAS(:,lv), 0, 0.0,0,0.0,&
                       1,LB,NEQ,1,E011,.false.,&
                       FEATRHS,1,1,.true.,2,'RHS   ')
    call ZTIME(TTT1)
    
    ! set up linear form
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    call linf_buildVectorScalar (p_rdiscretisation,rlinform,.true.,rvector,coeff_RHS)
    
    call ZTIME(TTT2)
    print *,'Time1 RHS: ',TTT1-TTT0
    print *,'Time2 RHS: ',TTT2-TTT1

    !CALL storage_getbase_double (rvector%h_Ddata,p_Ddata)
    !p_Ddata = p_Ddata - DWORK(L(LB):L(LB)+NEQ-1)
    
    print *,IWORK,IWMAX,NNWORK
    
    call storage_getbase_double (rmatrix%h_DA,p_Ddata)
    
    call lsyssc_releaseVector (rvector)
    !CALL lsyssc_releaseMatrix (rmatrix)
    call boundary_release (p_rboundary)
    
  end subroutine

!***************************************

  subroutine test_RHSEM30
  
    use collection
    use triangulation
    use spatialdiscretisation

    implicit none

    !INCLUDE 'stria.inc'
    include 'cout.inc'
    include 'cerr.inc'
    include 'cmem.inc'
    include 'cparametrization.inc'

    ! Variables
    
    type(t_boundary), pointer :: p_rboundary
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_spatialDiscretisation), pointer :: p_rdiscretisation
    real(DP) :: x,y,t
    integer, dimension(SZTRIA,20) :: TRIAS
    character(LEN=60) :: CFILE
    real(DP), dimension(:), pointer :: p_Ddata
    
    type(t_matrixScalar) :: rmatrix
    type(t_vectorScalar) :: rvector
    
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinform
    
    integer :: lv, LCOL,LLD,NA,NEQ,LA,LB
    external E011,COEAGS,EM30,FEATRHS,NDFGX
    double precision FEATRHS
    integer :: NDFGX
    integer, dimension(2,2) :: KAB
    real(DP) :: TTT0,TTT1,TTT2

    call storage_init(999, 100)
    call boundary_read_prm(p_rboundary, 'pre/QUAD.prm')
    
    lv = 10
    M = 0
    ICHECK = 0
    IMESH = 1

    CFILE = 'pre/QUAD.prm'
    call GENPAR (.true.,IMESH,CFILE)

    CFILE = 'pre/QUAD.tri'
    call INMTRI (2,TRIAS,lv,lv,0,CFILE)
    call tria_wrp_tria2Structure(TRIAS(:,lv),p_rtriangulation)
    
    call spdiscr_initDiscr_simple (p_rdiscretisation,EL_EM30,CUB_TRZ,&
                                   p_rtriangulation, p_rboundary, NULL())
    
    call ZTIME(TTT0)
    
    LB = 0
    NEQ = NDFGX(30,TRIAS(:,lv))
    call XVB0X(TRIAS(:,lv), 0, 0.0,0,0.0,&
                       1,LB,NEQ,1,EM30,.true.,&
                       FEATRHS,1,1,.true.,2,'RHS   ')
    call ZTIME(TTT1)
    
    ! set up linear form
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    call linf_buildVectorScalar (p_rdiscretisation,rlinform,.true.,rvector,coeff_RHS)
    
    call ZTIME(TTT2)
    
    print *,'Time1 RHS: ',TTT1-TTT0
    print *,'Time2 RHS: ',TTT2-TTT1

    print *,IWORK,IWMAX,NNWORK

    !CALL storage_getbase_double (rvector%h_Ddata,p_Ddata)
    !p_Ddata = p_Ddata - DWORK(L(LB):L(LB)+NEQ-1)
    
    call storage_getbase_double (rmatrix%h_DA,p_Ddata)
    
    call lsyssc_releaseVector (rvector)
    !CALL lsyssc_releaseMatrix (rmatrix)
    call boundary_release (p_rboundary)
    
  end subroutine

!***************************************

  subroutine test_LGS_Q1
  
    use collection
    use triangulation
    use spatialdiscretisation

    implicit none

    !INCLUDE 'stria.inc'
    include 'cout.inc'
    include 'cerr.inc'
    include 'cmem.inc'
    include 'cparametrization.inc'

    ! Variables
    
    type(t_boundary), pointer :: p_rboundary
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_spatialDiscretisation),pointer :: p_rdiscretisation
    real(DP) :: x,y,t
    integer, dimension(SZTRIA,20) :: TRIAS
    character(LEN=60) :: CFILE
    real(DP), dimension(:), pointer :: p_Ddata
    integer(I32), dimension(:), pointer :: p_Kcol,p_Kld
    
    type(t_matrixScalar) :: rmatrix
    type(t_vectorScalar) :: rvector
    
    type(t_bilinearForm) :: rform
    type(t_linearForm) :: rlinform
    
    integer :: lv, LCOL,LLD,NA,NEQ,LA,LB,ierror
    external E011,COEAGS,EM30,FEATRHS,NDFGX
    double precision FEATRHS
    integer :: NDFGX
    integer, dimension(2,2) :: KAB
    real(DP) :: TTT0,TTT1,TTT2
    
    type(t_boundaryConditions), pointer :: p_rboundaryConditions
    type(t_boundaryRegion) :: rboundaryRegion
    type(t_bcRegion), pointer :: p_rbcRegion
    type(t_discreteBC), pointer :: p_rdiscreteBC
    
    type(t_matrixBlock) :: rmatrixBlock
    type(t_vectorBlock) :: rvectorBlock,rrhsBlock,rtempBlock
    
    type(t_linsolNode), pointer :: p_rsolverNode
    
    type(t_filterChain), dimension(1), target :: RfilterChain
    type(t_filterChain), dimension(:), pointer :: p_RfilterChain
    
    integer NCELLS,NVERTS

    call storage_init(999, 100)
    call boundary_read_prm(p_rboundary, 'pre/QUAD.prm')
        
    lv = 7
    M = 0
    ICHECK = 0
    IMESH = 1

    CFILE = 'pre/QUAD.prm'
    call GENPAR (.true.,IMESH,CFILE)

    CFILE = 'pre/QUAD.tri'
    call INMTRI (2,TRIAS,lv,lv,0,CFILE)
    call tria_wrp_tria2Structure(TRIAS(:,lv),p_rtriangulation)
    
    call spdiscr_initDiscr_simple (p_rdiscretisation,EL_E011,CUB_TRZ,&
                                   p_rtriangulation, p_rboundary, NULL())
    call ZTIME(TTT0)
    call XAP7X(LCOL,LLD,NA,NEQ,TRIAS(:,lv),E011,0)
    
    call ZTIME(TTT1)
    call bilf_createMatrixStructure (p_rdiscretisation,LSYSSC_MATRIX9,rmatrix)
    
    call ZTIME(TTT2)
    print *,'Time1 structure: ',TTT1-TTT0
    print *,'Time2 structure: ',TTT2-TTT1
    
    call ZTIME(TTT0)
    KAB = reshape((/2,2,3,3/),(/2,2/))
    LA = 0
    call XAB7X(LA,LCOL,LLD,NA,NEQ,1,1,TRIAS(:,lv),E011,.false., &
               COEAGS,.true.,KAB,2,2,0,'MATRIX',0,0.0_DP)
    call ZTIME(TTT1)
    
    !CFILE = 'MATRIX1.TXT'
    !CALL OWM17(DWORK(L(LA)),KWORK(L(LCOL)),KWORK(L(LLD)),&
    !           NEQ,NEQ,.TRUE.,0,'MAT1  ',CFILE,'(D20.10)')
    
    ! set up bilinear form
    rform%itermCount = 2
    rform%Idescriptors(1,1) = DER_DERIV_X
    rform%Idescriptors(2,1) = DER_DERIV_X
    rform%Idescriptors(1,2) = DER_DERIV_Y
    rform%Idescriptors(2,2) = DER_DERIV_Y
    rform%Dcoefficients(1)  = 1.0
    rform%Dcoefficients(2)  = 1.0
    rform%ballCoeffConstant = .true.
    rform%BconstantCoeff = .true.
    !BILF_NELEMSIM = 100
    !BILF_NELEMSIM = 300000
    call bilf_buildMatrixScalar (rform,.true.,rmatrix,coeff_Laplace)
    call ZTIME(TTT2)

    print *,'Time1 assembly: ',TTT1-TTT0
    print *,'Time2 assembly: ',TTT2-TTT1
    
    print *,IWORK,IWMAX,NNWORK
    
    LB = 0
    NEQ = NDFGX(30,TRIAS(:,lv))
    call XVB0X(TRIAS(:,lv), 0, 0.0,0,0.0,&
                       1,LB,NEQ,1,E011,.false.,&
                       FEATRHS,1,1,.true.,2,'RHS   ')
    call ZTIME(TTT1)
    
    ! set up linear form
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    call linf_buildVectorScalar (p_rdiscretisation,rlinform,.true.,rvector,coeff_RHS)
    
    call ZTIME(TTT2)
    print *,'Time1 RHS: ',TTT1-TTT0
    print *,'Time2 RHS: ',TTT2-TTT1

    call storage_getbase_double (rvector%h_Ddata,p_Ddata)
    call storage_getbase_double (rmatrix%h_DA,p_Ddata)
    
    ! Structure for boundary conditions
    call scbc_initScalarBC (p_rboundaryConditions,p_rboundary)
    
    ! Add all four edges as Dirichlet
    call boundary_createRegion(p_rboundary,1,1,rboundaryRegion)
    call scbc_newBConRealBD (BC_DIRICHLET,BC_RTYPE_REAL,p_rboundaryConditions, &
                             rboundaryRegion,p_rbcRegion)
                             
    call boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
    call scbc_newBConRealBD (BC_DIRICHLET,BC_RTYPE_REAL,p_rboundaryConditions, &
                             rboundaryRegion,p_rbcRegion)
                             
    call boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
    call scbc_newBConRealBD (BC_DIRICHLET,BC_RTYPE_REAL,p_rboundaryConditions, &
                             rboundaryRegion,p_rbcRegion)
                             
    call boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
    call scbc_newBConRealBD (BC_DIRICHLET,BC_RTYPE_REAL,p_rboundaryConditions, &
                             rboundaryRegion,p_rbcRegion)
                             
    p_rdiscretisation%p_rboundaryConditions => p_rboundaryConditions

    ! Discretise the boundary conditions; this gives p_RdiscreteBC
    nullify(p_RdiscreteBC)
    call bcasm_discretiseBC (p_rdiscretisation,p_RdiscreteBC,.false., &
                             getBoundaryValues,NULL())
                             
    ! Hang them in into the vector and matrix
    rmatrix%p_RdiscreteBC => p_RdiscreteBC
    rvector%p_RdiscreteBC => p_RdiscreteBC
                             
    ! Create block matrices/vectors for the solver
    call lsysbl_createMatFromScalar (rmatrix,rmatrixBlock)
    call lsysbl_createVecFromScalar (rvector,rrhsBlock)
    call lsysbl_createVecBlockIndirect (rrhsBlock, rvectorBlock, .true.)
    call lsysbl_createVecBlockIndirect (rrhsBlock, rtempBlock, .false.)
    
    ! Initialise the first filter of the filter chain as boundary
    ! implementation filter:
    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL
    
    ! Apply the filter chain to the matrix and the vectors.
    ! This implement them into the vectors and matrices
    call filter_applyFilterChainVec (rrhsBlock, RfilterChain)
    call filter_applyFilterChainVec (rvectorBlock, RfilterChain)
    call filter_applyFilterChainMat (rmatrixBlock, RfilterChain)
    
    ! Change the filter to work with defect vectors during the
    ! solution process
    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

    ! Create a BiCGStab-solver with the above filter chain
    p_RfilterChain => RfilterChain
    call linsol_initBiCGStab (p_rsolverNode,NULL(),p_RfilterChain)
    
    ! Initialise the solver with the system matrix
    call linsol_setMatrices(p_RsolverNode,(/rmatrixBlock/))
    
    ! Initialise structure/data of the solver
    call linsol_initStructure (p_rsolverNode,ierror)
    call linsol_initData (p_rsolverNode,ierror)
    
    ! Solve the system
    call linsol_solveAdaptively (p_rsolverNode,&
                                 rvectorBlock,rrhsBlock,rtempBlock)
    
    ! Release unnecessary data
    call linsol_doneData (p_rsolverNode)
    call linsol_doneStructure (p_rsolverNode)
    call linsol_releaseSolver (p_rsolverNode)

    call bcasm_releaseDiscreteBC (p_RdiscreteBC)
    call lsyssc_releaseVector (rvector)
    call lsyssc_releaseMatrix (rmatrix)
    call boundary_release (p_rboundary)
    
    call GMVOF0 (69,-2,'u.gmv')
    call GMVHEA (69)
    call GMVTRI (69,p_rtriangulation%Itria,0,NCELLS,NVERTS)
    
    call storage_getbase_double (rvectorBlock%RvectorBlock(1)%h_Ddata,p_Ddata)
    call GMVSCA (69,p_rtriangulation%Itria,1,NVERTS,&
                 rvectorBlock%RvectorBlock(1)%NEQ,p_Ddata,'sol')
    
    call GMVFOT (69)
    close(69)
    
  end subroutine

end module

    double precision function FEATRHS (X,Y,IA,IBLOC,BFIRST,&
                                    TRIA,IPARAM,DPARAM,IGEOM,DGEOM)
    
    implicit none
    
    double precision X,Y
    integer IA, IB, IBLOC
    logical BFIRST
    
    integer IPARAM(*),TRIA(*),IGEOM(*)
    double precision DPARAM(*),DGEOM(*)

    double precision TIMENS,RE
    
    FEATRHS=1.0D0

    end function


program testf90
  use testcode
  implicit none
  include 'cmem.inc'
  call ZINIT(NNWORK,'feat.msg','data/cc2d.err','data/cc2d.prt',&
             'data/cc2d.sys','data/cc2d.trc')
  call test_MatrixEM30
end program
