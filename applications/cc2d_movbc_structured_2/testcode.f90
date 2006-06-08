MODULE testcode

  USE fsystem
  USE storage
  USE linearsolver
  USE boundary
  USE bilinearformevaluation
  USE linearformevaluation
  USE cubature
  USE matrixfilters
  USE vectorfilters
  USE bcassembly
  
  IMPLICIT NONE


CONTAINS

  SUBROUTINE coeff_Laplace (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTrial,IdofsTest,rdomainIntSubset, p_rcollection,&
                  Dcoefficients)
    
    USE BasicGeometry
    USE triangulation
    USE collection
    USE scalarpde
    USE domainintegration
    
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
    TYPE(t_spatialDiscretisation), INTENT(IN)                   :: rdiscretisation
    
    ! The bilinear form which is currently being evaluated:
    TYPE(t_bilinearForm), INTENT(IN)                            :: rform
    
    ! Number of elements, where the coefficients must be computed.
    INTEGER, INTENT(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    INTEGER, INTENT(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    REAL(DP), DIMENSION(NDIM2D,npointsPerElement,nelements), INTENT(IN)  :: Dpoints
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in trial space,Number of elements)
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(IN) :: IdofsTrial
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,Number of elements)
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(IN) :: IdofsTest
    
    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    TYPE(t_domainIntSubset), INTENT(IN)              :: rdomainIntSubset

    ! A pointer to a collection structure to provide additional 
    ! information to the coefficient routine. May point to NULL() if not defined.
    TYPE(t_collection), POINTER                      :: p_rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    REAL(DP), DIMENSION(:,:,:), INTENT(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

  Dcoefficients = 1.0_DP

  END SUBROUTINE

!***********************************

  !<subroutine>

    SUBROUTINE coeff_RHS (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,p_rcollection, &
                  Dcoefficients)
    
    USE BasicGeometry
    USE triangulation
    USE collection
    USE scalarpde
    USE domainintegration
    
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
    TYPE(t_spatialDiscretisation), INTENT(IN)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    TYPE(t_linearForm), INTENT(IN)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    INTEGER, INTENT(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    INTEGER, INTENT(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    REAL(DP), DIMENSION(NDIM2D,npointsPerElement,nelements), INTENT(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,Number of elements)
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    TYPE(t_domainIntSubset), INTENT(IN)              :: rdomainIntSubset

    ! A pointer to a collection structure to provide additional 
    ! information to the coefficient routine. May point to NULL() if not defined.
    TYPE(t_collection), POINTER                      :: p_rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    REAL(DP), DIMENSION(:,:,:), INTENT(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    Dcoefficients = 1.0_DP
  
    END SUBROUTINE

  !<subroutine>

    SUBROUTINE getBoundaryValues (rdiscretisation,rbcRegion,ielement, &
                                   cinfoNeeded,iwhere,dwhere, p_rcollection, Dvalues)
    
    USE collection
    USE spatialdiscretisation
    USE discretebc
    
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
    TYPE(t_spatialDiscretisation), INTENT(IN)                   :: rdiscretisation
    
    ! Boundary condition region that is currently being processed.
    ! (This e.g. defines the type of boundary conditions that are
    !  currently being calculated, as well as information about the current
    !  boundary segment 'where we are at the moment'.)
    TYPE(t_bcRegion), INTENT(IN)                                :: rbcRegion
    
    
    ! The element number on the boundary which is currently being processed
    INTEGER(I32), INTENT(IN)                                    :: ielement
    
    ! The type of information, the routine should calculate. One of the
    ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
    ! to return one or multiple information value in the result array.
    INTEGER, INTENT(IN)                                         :: cinfoNeeded
    
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
    INTEGER, INTENT(IN)                                         :: iwhere

    ! A reference to a geometric object where information should be computed.
    ! cinfoNeeded=DISCBC_NEEDFUNC : 
    !   dwhere = parameter value of the point where the value should be computed,
    ! cinfoNeeded=DISCBC_NEEDDERIV : 
    !   dwhere = parameter value of the point where the value should be computed,
    ! cinfoNeeded=DISCBC_NEEDINTMEAN : 
    !   dwhere = 0 (not used)
    REAL(DP), INTENT(IN)                                        :: dwhere
     
    ! A pointer to a collection structure to provide additional 
    ! information to the coefficient routine. May point to NULL() if not defined.
    TYPE(t_collection), POINTER                  :: p_rcollection

  !</input>
  
  !<output>
    ! This array receives the calculated information. If the caller
    ! only needs one value, the computed quantity is put into Dvalues(1). 
    ! If multiple values are needed, they are collected here (e.g. for 
    ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
    REAL(DP), DIMENSION(:), INTENT(OUT)                         :: Dvalues
  !</output>
    
  !</subroutine>
  
    Dvalues(1) = 0.0_DP
  
    END SUBROUTINE

!***************************************

  SUBROUTINE test_MatrixQ1
  
    use collection
    use triangulation
    use spatialdiscretisation

    implicit none

    !INCLUDE 'stria.inc'
    INCLUDE 'cout.inc'
    INCLUDE 'cerr.inc'
    INCLUDE 'cmem.inc'
    INCLUDE 'cparametrization.inc'

    ! Variables
    
    TYPE(t_boundary), POINTER :: p_rboundary
    TYPE(t_triangulation), POINTER :: p_rtriangulation
    TYPE(t_spatialDiscretisation), POINTER :: p_rdiscretisation
    REAL(DP) :: x,y,t
    INTEGER, DIMENSION(SZTRIA,20) :: TRIAS
    CHARACTER(LEN=60) :: CFILE
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    
    TYPE(t_matrixScalar) :: rmatrix
    TYPE(t_vectorScalar) :: rvector
    
    TYPE(t_bilinearForm) :: rform
    TYPE(t_linearForm) :: rlinform
    
    INTEGER :: lv, LCOL,LLD,NA,NEQ,LA,LB
    EXTERNAL E011,COEAGS,EM30,FEATRHS,NDFGX
    DOUBLE PRECISION FEATRHS 
    INTEGER :: NDFGX
    INTEGER, DIMENSION(2,2) :: KAB
    REAL(DP) :: TTT0,TTT1,TTT2

    CALL storage_init(999, 100)
    CALL boundary_read_prm(p_rboundary, 'pre/QUAD.prm')
    
    lv = 10
    M = 0
    ICHECK = 0
    IMESH = 1

    CFILE = 'pre/QUAD.prm'
    CALL GENPAR (.TRUE.,IMESH,CFILE)

    CFILE = 'pre/QUAD.tri'
    CALL INMTRI (2,TRIAS,lv,lv,0,CFILE)
    CALL tria_wrp_tria2Structure(TRIAS(:,lv),p_rtriangulation)
    
    CALL spdiscr_initDiscr_simple (p_rdiscretisation,EL_E011,CUB_TRZ,&
                                   p_rtriangulation, p_rboundary, NULL())
    CALL ZTIME(TTT0)
    CALL XAP7X(LCOL,LLD,NA,NEQ,TRIAS(:,lv),E011,0)
    
    CALL ZTIME(TTT1)
    CALL bilf_createMatrixStructure (p_rdiscretisation,LSYSSC_MATRIX9,rmatrix)
    
    CALL ZTIME(TTT2)
    PRINT *,'Time1 structure: ',TTT1-TTT0
    PRINT *,'Time2 structure: ',TTT2-TTT1
    
    CALL ZTIME(TTT0)
    KAB = RESHAPE((/2,2,3,3/),(/2,2/))
    LA = 0
    CALL XAB7X(LA,LCOL,LLD,NA,NEQ,1,1,TRIAS(:,lv),E011,.FALSE., &
               COEAGS,.TRUE.,KAB,2,2,0,'MATRIX',0,0.0_DP)
    CALL ZTIME(TTT1)
    
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
    rform%ballCoeffConstant = .TRUE.
    rform%BconstantCoeff = .TRUE.
    !BILF_NELEMSIM = 100
    !BILF_NELEMSIM = 300000
    CALL bilf_buildMatrixScalar (p_rdiscretisation,rform,.TRUE.,rmatrix,coeff_Laplace)
    CALL ZTIME(TTT2)

    PRINT *,'Time1 assembly: ',TTT1-TTT0
    PRINT *,'Time2 assembly: ',TTT2-TTT1
    
    PRINT *,IWORK,IWMAX,NNWORK
    
    !CALL lsyssc_releaseVector (rvector)
    CALL lsyssc_releaseMatrix (rmatrix)
    CALL boundary_release (p_rboundary)
    
  END SUBROUTINE

!***************************************

  SUBROUTINE test_MatrixEM30
  
    use collection
    use triangulation
    use spatialdiscretisation

    implicit none

    !INCLUDE 'stria.inc'
    INCLUDE 'cout.inc'
    INCLUDE 'cerr.inc'
    INCLUDE 'cmem.inc'
    INCLUDE 'cparametrization.inc'

    ! Variables
    
    TYPE(t_boundary), POINTER :: p_rboundary
    TYPE(t_triangulation), POINTER :: p_rtriangulation
    TYPE(t_spatialDiscretisation),POINTER :: p_rdiscretisation
    REAL(DP) :: x,y,t
    INTEGER, DIMENSION(SZTRIA,20) :: TRIAS
    CHARACTER(LEN=60) :: CFILE
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    
    TYPE(t_matrixScalar) :: rmatrix
    TYPE(t_vectorScalar) :: rvector
    
    TYPE(t_bilinearForm) :: rform
    TYPE(t_linearForm) :: rlinform
    
    INTEGER :: lv, LCOL,LLD,NA,NEQ,LA,LB
    EXTERNAL E011,COEAGS,EM30,FEATRHS,NDFGX
    DOUBLE PRECISION FEATRHS 
    INTEGER :: NDFGX
    INTEGER, DIMENSION(2,2) :: KAB
    REAL(DP) :: TTT0,TTT1,TTT2

    CALL storage_init(999, 100)
    CALL boundary_read_prm(p_rboundary, 'pre/QUAD.prm')
    
    lv = 10
    M = 0
    ICHECK = 0
    IMESH = 1

    CFILE = 'pre/QUAD.prm'
    CALL GENPAR (.TRUE.,IMESH,CFILE)

    CFILE = 'pre/QUAD.tri'
    CALL INMTRI (2,TRIAS,lv,lv,0,CFILE)
    CALL tria_wrp_tria2Structure(TRIAS(:,lv),p_rtriangulation)
    
    CALL spdiscr_initDiscr_simple (p_rdiscretisation,EL_EM30,CUB_TRZ,&
                                   p_rtriangulation, p_rboundary, NULL())
    CALL ZTIME(TTT0)
    CALL XAP7X(LCOL,LLD,NA,NEQ,TRIAS(:,lv),EM30,0)
    
    CALL ZTIME(TTT1)
    CALL bilf_createMatrixStructure (p_rdiscretisation,LSYSSC_MATRIX9,rmatrix)
    
    CALL ZTIME(TTT2)
    PRINT *,'Time1 structure: ',TTT1-TTT0
    PRINT *,'Time2 structure: ',TTT2-TTT1
    
    CALL ZTIME(TTT0)
    KAB = RESHAPE((/2,2,3,3/),(/2,2/))
    LA = 0
    CALL XAB7X(LA,LCOL,LLD,NA,NEQ,1,1,TRIAS(:,lv),EM30,.TRUE., &
               COEAGS,.TRUE.,KAB,2,2,0,'MATRIX',0,0.0_DP)
    CALL ZTIME(TTT1)
    
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
    rform%ballCoeffConstant = .TRUE.
    rform%BconstantCoeff = .TRUE.
    !BILF_NELEMSIM = 100
    !BILF_NELEMSIM = 300000
    CALL bilf_buildMatrixScalar (p_rdiscretisation,rform,.TRUE.,rmatrix,coeff_Laplace)
    CALL ZTIME(TTT2)

    PRINT *,'Time1 assembly: ',TTT1-TTT0
    PRINT *,'Time2 assembly: ',TTT2-TTT1
    
    PRINT *,IWORK,IWMAX,NNWORK
    
    !CALL lsyssc_releaseVector (rvector)
    CALL lsyssc_releaseMatrix (rmatrix)
    CALL boundary_release (p_rboundary)
    
  END SUBROUTINE

!***************************************

  SUBROUTINE test_RHSQ1
  
    use collection
    use triangulation
    use spatialdiscretisation

    implicit none

    !INCLUDE 'stria.inc'
    INCLUDE 'cout.inc'
    INCLUDE 'cerr.inc'
    INCLUDE 'cmem.inc'
    INCLUDE 'cparametrization.inc'

    ! Variables
    
    TYPE(t_boundary), POINTER :: p_rboundary
    TYPE(t_triangulation), POINTER :: p_rtriangulation
    TYPE(t_spatialDiscretisation), POINTER :: p_rdiscretisation
    REAL(DP) :: x,y,t
    INTEGER, DIMENSION(SZTRIA,20) :: TRIAS
    CHARACTER(LEN=60) :: CFILE
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    
    TYPE(t_matrixScalar) :: rmatrix
    TYPE(t_vectorScalar) :: rvector
    
    TYPE(t_bilinearForm) :: rform
    TYPE(t_linearForm) :: rlinform
    
    INTEGER :: lv, LCOL,LLD,NA,NEQ,LA,LB
    EXTERNAL E011,COEAGS,EM30,FEATRHS,NDFGX
    DOUBLE PRECISION FEATRHS 
    INTEGER :: NDFGX
    INTEGER, DIMENSION(2,2) :: KAB
    REAL(DP) :: TTT0,TTT1,TTT2

    CALL storage_init(999, 100)
    CALL boundary_read_prm(p_rboundary, 'pre/QUAD.prm')
    
    lv = 10
    M = 0
    ICHECK = 0
    IMESH = 1

    CFILE = 'pre/QUAD.prm'
    CALL GENPAR (.TRUE.,IMESH,CFILE)

    CFILE = 'pre/QUAD.tri'
    CALL INMTRI (2,TRIAS,lv,lv,0,CFILE)
    CALL tria_wrp_tria2Structure(TRIAS(:,lv),p_rtriangulation)
    
    CALL spdiscr_initDiscr_simple (p_rdiscretisation,EL_E011,CUB_TRZ,&
                                   p_rtriangulation, p_rboundary, NULL())
    
    CALL ZTIME(TTT0)
    
    LB = 0
    NEQ = NDFGX(30,TRIAS(:,lv))
    CALL XVB0X(TRIAS(:,lv), 0, 0.0,0,0.0,&
                       1,LB,NEQ,1,E011,.FALSE.,&
                       FEATRHS,1,1,.TRUE.,2,'RHS   ') 
    CALL ZTIME(TTT1)
    
    ! set up linear form
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    CALL linf_buildVectorScalar (p_rdiscretisation,rlinform,.TRUE.,rvector,coeff_RHS)
    
    CALL ZTIME(TTT2)
    PRINT *,'Time1 RHS: ',TTT1-TTT0
    PRINT *,'Time2 RHS: ',TTT2-TTT1

    !CALL storage_getbase_double (rvector%h_Ddata,p_Ddata)
    !p_Ddata = p_Ddata - DWORK(L(LB):L(LB)+NEQ-1)
    
    PRINT *,IWORK,IWMAX,NNWORK
    
    CALL storage_getbase_double (rmatrix%h_DA,p_Ddata)
    
    CALL lsyssc_releaseVector (rvector)
    !CALL lsyssc_releaseMatrix (rmatrix)
    CALL boundary_release (p_rboundary)
    
  END SUBROUTINE

!***************************************

  SUBROUTINE test_RHSEM30
  
    use collection
    use triangulation
    use spatialdiscretisation

    implicit none

    !INCLUDE 'stria.inc'
    INCLUDE 'cout.inc'
    INCLUDE 'cerr.inc'
    INCLUDE 'cmem.inc'
    INCLUDE 'cparametrization.inc'

    ! Variables
    
    TYPE(t_boundary), POINTER :: p_rboundary
    TYPE(t_triangulation), POINTER :: p_rtriangulation
    TYPE(t_spatialDiscretisation), POINTER :: p_rdiscretisation
    REAL(DP) :: x,y,t
    INTEGER, DIMENSION(SZTRIA,20) :: TRIAS
    CHARACTER(LEN=60) :: CFILE
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    
    TYPE(t_matrixScalar) :: rmatrix
    TYPE(t_vectorScalar) :: rvector
    
    TYPE(t_bilinearForm) :: rform
    TYPE(t_linearForm) :: rlinform
    
    INTEGER :: lv, LCOL,LLD,NA,NEQ,LA,LB
    EXTERNAL E011,COEAGS,EM30,FEATRHS,NDFGX
    DOUBLE PRECISION FEATRHS 
    INTEGER :: NDFGX
    INTEGER, DIMENSION(2,2) :: KAB
    REAL(DP) :: TTT0,TTT1,TTT2

    CALL storage_init(999, 100)
    CALL boundary_read_prm(p_rboundary, 'pre/QUAD.prm')
    
    lv = 10
    M = 0
    ICHECK = 0
    IMESH = 1

    CFILE = 'pre/QUAD.prm'
    CALL GENPAR (.TRUE.,IMESH,CFILE)

    CFILE = 'pre/QUAD.tri'
    CALL INMTRI (2,TRIAS,lv,lv,0,CFILE)
    CALL tria_wrp_tria2Structure(TRIAS(:,lv),p_rtriangulation)
    
    CALL spdiscr_initDiscr_simple (p_rdiscretisation,EL_EM30,CUB_TRZ,&
                                   p_rtriangulation, p_rboundary, NULL())
    
    CALL ZTIME(TTT0)
    
    LB = 0
    NEQ = NDFGX(30,TRIAS(:,lv))
    CALL XVB0X(TRIAS(:,lv), 0, 0.0,0,0.0,&
                       1,LB,NEQ,1,EM30,.TRUE.,&
                       FEATRHS,1,1,.TRUE.,2,'RHS   ') 
    CALL ZTIME(TTT1)
    
    ! set up linear form
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    CALL linf_buildVectorScalar (p_rdiscretisation,rlinform,.TRUE.,rvector,coeff_RHS)
    
    CALL ZTIME(TTT2)
    
    PRINT *,'Time1 RHS: ',TTT1-TTT0
    PRINT *,'Time2 RHS: ',TTT2-TTT1

    PRINT *,IWORK,IWMAX,NNWORK

    !CALL storage_getbase_double (rvector%h_Ddata,p_Ddata)
    !p_Ddata = p_Ddata - DWORK(L(LB):L(LB)+NEQ-1)
    
    CALL storage_getbase_double (rmatrix%h_DA,p_Ddata)
    
    CALL lsyssc_releaseVector (rvector)
    !CALL lsyssc_releaseMatrix (rmatrix)
    CALL boundary_release (p_rboundary)
    
  END SUBROUTINE

!***************************************

  SUBROUTINE test_LGS_Q1
  
    use collection
    use triangulation
    use spatialdiscretisation

    implicit none

    !INCLUDE 'stria.inc'
    INCLUDE 'cout.inc'
    INCLUDE 'cerr.inc'
    INCLUDE 'cmem.inc'
    INCLUDE 'cparametrization.inc'

    ! Variables
    
    TYPE(t_boundary), POINTER :: p_rboundary
    TYPE(t_triangulation), POINTER :: p_rtriangulation
    TYPE(t_spatialDiscretisation),POINTER :: p_rdiscretisation
    REAL(DP) :: x,y,t
    INTEGER, DIMENSION(SZTRIA,20) :: TRIAS
    CHARACTER(LEN=60) :: CFILE
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    INTEGER(I32), DIMENSION(:), POINTER :: p_Kcol,p_Kld
    
    TYPE(t_matrixScalar) :: rmatrix
    TYPE(t_vectorScalar) :: rvector
    
    TYPE(t_bilinearForm) :: rform
    TYPE(t_linearForm) :: rlinform
    
    INTEGER :: lv, LCOL,LLD,NA,NEQ,LA,LB,ierror
    EXTERNAL E011,COEAGS,EM30,FEATRHS,NDFGX
    DOUBLE PRECISION FEATRHS 
    INTEGER :: NDFGX
    INTEGER, DIMENSION(2,2) :: KAB
    REAL(DP) :: TTT0,TTT1,TTT2
    
    TYPE(t_boundaryConditions), POINTER :: p_rboundaryConditions
    TYPE(t_boundaryRegion) :: rboundaryRegion
    TYPE(t_bcRegion), POINTER :: p_rbcRegion
    TYPE(t_discreteBC), POINTER :: p_rdiscreteBC
    
    TYPE(t_matrixBlock) :: rmatrixBlock
    TYPE(t_vectorBlock) :: rvectorBlock,rrhsBlock,rtempBlock
    
    TYPE(t_linsolNode), POINTER :: p_rsolverNode
    
    TYPE(t_filterChain), DIMENSION(1), TARGET :: RfilterChain
    TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
    
    INTEGER NCELLS,NVERTS

    CALL storage_init(999, 100)
    CALL boundary_read_prm(p_rboundary, 'pre/QUAD.prm')
        
    lv = 7
    M = 0
    ICHECK = 0
    IMESH = 1

    CFILE = 'pre/QUAD.prm'
    CALL GENPAR (.TRUE.,IMESH,CFILE)

    CFILE = 'pre/QUAD.tri'
    CALL INMTRI (2,TRIAS,lv,lv,0,CFILE)
    CALL tria_wrp_tria2Structure(TRIAS(:,lv),p_rtriangulation)
    
    CALL spdiscr_initDiscr_simple (p_rdiscretisation,EL_E011,CUB_TRZ,&
                                   p_rtriangulation, p_rboundary, NULL())
    CALL ZTIME(TTT0)
    CALL XAP7X(LCOL,LLD,NA,NEQ,TRIAS(:,lv),E011,0)
    
    CALL ZTIME(TTT1)
    CALL bilf_createMatrixStructure (p_rdiscretisation,LSYSSC_MATRIX9,rmatrix)
    
    CALL ZTIME(TTT2)
    PRINT *,'Time1 structure: ',TTT1-TTT0
    PRINT *,'Time2 structure: ',TTT2-TTT1
    
    CALL ZTIME(TTT0)
    KAB = RESHAPE((/2,2,3,3/),(/2,2/))
    LA = 0
    CALL XAB7X(LA,LCOL,LLD,NA,NEQ,1,1,TRIAS(:,lv),E011,.FALSE., &
               COEAGS,.TRUE.,KAB,2,2,0,'MATRIX',0,0.0_DP)
    CALL ZTIME(TTT1)
    
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
    rform%ballCoeffConstant = .TRUE.
    rform%BconstantCoeff = .TRUE.
    !BILF_NELEMSIM = 100
    !BILF_NELEMSIM = 300000
    CALL bilf_buildMatrixScalar (p_rdiscretisation,rform,.TRUE.,rmatrix,coeff_Laplace)
    CALL ZTIME(TTT2)

    PRINT *,'Time1 assembly: ',TTT1-TTT0
    PRINT *,'Time2 assembly: ',TTT2-TTT1
    
    PRINT *,IWORK,IWMAX,NNWORK
    
    LB = 0
    NEQ = NDFGX(30,TRIAS(:,lv))
    CALL XVB0X(TRIAS(:,lv), 0, 0.0,0,0.0,&
                       1,LB,NEQ,1,E011,.FALSE.,&
                       FEATRHS,1,1,.TRUE.,2,'RHS   ') 
    CALL ZTIME(TTT1)
    
    ! set up linear form
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    CALL linf_buildVectorScalar (p_rdiscretisation,rlinform,.TRUE.,rvector,coeff_RHS)
    
    CALL ZTIME(TTT2)
    PRINT *,'Time1 RHS: ',TTT1-TTT0
    PRINT *,'Time2 RHS: ',TTT2-TTT1

    CALL storage_getbase_double (rvector%h_Ddata,p_Ddata)
    CALL storage_getbase_double (rmatrix%h_DA,p_Ddata)
    
    ! Structure for boundary conditions
    CALL scbc_initScalarBC (p_rboundaryConditions,p_rboundary)
    
    ! Add all four edges as Dirichlet
    CALL boundary_createRegion(p_rboundary,1,1,rboundaryRegion)
    CALL scbc_newBConRealBD (BC_DIRICHLET,BC_RTYPE_REAL,p_rboundaryConditions, &
                             rboundaryRegion,p_rbcRegion)
                             
    CALL boundary_createRegion(p_rboundary,1,2,rboundaryRegion)
    CALL scbc_newBConRealBD (BC_DIRICHLET,BC_RTYPE_REAL,p_rboundaryConditions, &
                             rboundaryRegion,p_rbcRegion)
                             
    CALL boundary_createRegion(p_rboundary,1,3,rboundaryRegion)
    CALL scbc_newBConRealBD (BC_DIRICHLET,BC_RTYPE_REAL,p_rboundaryConditions, &
                             rboundaryRegion,p_rbcRegion)
                             
    CALL boundary_createRegion(p_rboundary,1,4,rboundaryRegion)
    CALL scbc_newBConRealBD (BC_DIRICHLET,BC_RTYPE_REAL,p_rboundaryConditions, &
                             rboundaryRegion,p_rbcRegion)
                             
    p_rdiscretisation%p_rboundaryConditions => p_rboundaryConditions

    ! Discretise the boundary conditions; this gives p_RdiscreteBC
    NULLIFY(p_RdiscreteBC)
    CALL bcasm_discretiseBC (p_rdiscretisation,p_RdiscreteBC,.FALSE., &
                             getBoundaryValues,NULL())
                             
    ! Hang them in into the vector and matrix
    rmatrix%p_RdiscreteBC => p_RdiscreteBC
    rvector%p_RdiscreteBC => p_RdiscreteBC
                             
    ! Create block matrices/vectors for the solver
    CALL lsysbl_createMatFromScalar (rmatrix,rmatrixBlock)
    CALL lsysbl_createVecFromScalar (rvector,rrhsBlock)
    CALL lsysbl_createVecBlockIndirect (rrhsBlock, rvectorBlock, .TRUE.)
    CALL lsysbl_createVecBlockIndirect (rrhsBlock, rtempBlock, .FALSE.)
    
    ! Initialise the first filter of the filter chain as boundary
    ! implementation filter:
    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL
    
    ! Apply the filter chain to the matrix and the vectors.
    ! This implement them into the vectors and matrices
    CALL filter_applyFilterChainVec (rrhsBlock, RfilterChain)
    CALL filter_applyFilterChainVec (rvectorBlock, RfilterChain)
    CALL filter_applyFilterChainMat (rmatrixBlock, RfilterChain)
    
    ! Change the filter to work with defect vectors during the 
    ! solution process
    RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

    ! Create a BiCGStab-solver with the above filter chain
    p_RfilterChain => RfilterChain
    CALL linsol_initBiCGStab (p_rsolverNode,NULL(),p_RfilterChain)
    
    ! Initialise the solver with the system matrix
    CALL linsol_setMatrices(p_RsolverNode,(/rmatrixBlock/))
    
    ! Initialise structure/data of the solver
    CALL linsol_initStructure (p_rsolverNode,ierror)
    CALL linsol_initData (p_rsolverNode,ierror)
    
    ! Solve the system
    CALL linsol_solveAdaptively (p_rsolverNode,&
                                 rvectorBlock,rrhsBlock,rtempBlock)
    
    ! Release unnecessary data
    CALL linsol_doneData (p_rsolverNode)
    CALL linsol_doneStructure (p_rsolverNode)
    CALL linsol_releaseSolver (p_rsolverNode)

    CALL bcasm_releaseDiscreteBC (p_RdiscreteBC)
    CALL lsyssc_releaseVector (rvector)
    CALL lsyssc_releaseMatrix (rmatrix)
    CALL boundary_release (p_rboundary)
    
    CALL GMVOF0 (69,-2,'u.gmv')
    CALL GMVHEA (69)
    CALL GMVTRI (69,p_rtriangulation%Itria,0,NCELLS,NVERTS)
    
    CALL storage_getbase_double (rvectorBlock%RvectorBlock(1)%h_Ddata,p_Ddata)
    CALL GMVSCA (69,p_rtriangulation%Itria,1,NVERTS,&
                 rvectorBlock%RvectorBlock(1)%NEQ,p_Ddata,'sol')
    
    CALL GMVFOT (69)
    CLOSE(69)
    
  END SUBROUTINE

END MODULE

    DOUBLE PRECISION FUNCTION FEATRHS (X,Y,IA,IBLOC,BFIRST,&
                                    TRIA,IPARAM,DPARAM,IGEOM,DGEOM)
    
    IMPLICIT NONE
    
    DOUBLE PRECISION X,Y
    INTEGER IA, IB, IBLOC
    LOGICAL BFIRST
    
    INTEGER IPARAM(*),TRIA(*),IGEOM(*)
    DOUBLE PRECISION DPARAM(*),DGEOM(*)

    DOUBLE PRECISION TIMENS,RE
    
    FEATRHS=1.0D0

    END FUNCTION


PROGRAM testf90
  USE testcode
  IMPLICIT NONE
  INCLUDE 'cmem.inc'
  CALL ZINIT(NNWORK,'feat.msg','data/cc2d.err','data/cc2d.prt',&
             'data/cc2d.sys','data/cc2d.trc') 
  CALL test_MatrixEM30
END PROGRAM
