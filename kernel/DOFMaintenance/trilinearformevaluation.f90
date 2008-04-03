!##############################################################################
!# ****************************************************************************
!# <name> trilinearformevaluation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines for the discretisation of trilinear forms,
!# i.e. the creation of matrices and matrix structures. It contains the
!# following set of routines:
!#
!# 1.) trilf_buildMatrixScalar
!#     -> Assembles the entries of a matrix, which structure was build
!#        with bilf_createMatrixStructure before.
!#
!# A 'trilinear form' is in our context a bilinear form whose coefficient
!# function may depend on a Finite Element function:
!#
!#   $$ a(u,phi_i,psi_j)  =  \int c(x,y) f(u) g(\psi_j) h(\phi_i) $$
!#
!# with f,g,h derivative operators. Therefore, there's no
!# special 'create structure' routine; the structure is generated based on the
!# underlying bilinear form, which steps from the idea of treating u as
!# a variable coefficient:
!#
!#   $$ a_u(phi_i,psi_j)  =  \int c_u(x,y) g(\psi_j) h(\phi_i), $$
!#
!#   $$ c_u = c*f(u) $$
!# 
!# </purpose>
!##############################################################################

MODULE trilinearformevaluation

  USE fsystem
  USE linearsystemscalar
  USE spatialdiscretisation
  USE scalarpde
  USE derivatives
  USE cubature
  USE collection
  USE domainintegration
  USE element
  USE elementpreprocessing
  
  IMPLICIT NONE

!<constants>

!<constantblock description="Constants defining the blocking of the assembly">

  ! Number of elements to handle simultaneously when building matrices
  INTEGER, PARAMETER :: TRILF_NELEMSIM   = 1000
  
!</constantblock>
!</constants>

CONTAINS

  !****************************************************************************

!<subroutine>

  SUBROUTINE trilf_buildMatrixScalar (rform,bclear,rmatrixScalar,rvector,&
                                      fcoeff_buildTrilMatrixSc_sim,rcollection)
  
!<description>
  ! This routine calculates the entries of a finite element matrix using
  ! a trilinear form
  !     $$ a(u,phi_i,psi_j)  =  \int c(x,y) f(u) g(\psi_j) h(\phi_i) $$
  ! The function $u$ is specified in rvector. The derivative quantifier
  ! rform(Idescriptors(1,.) specifies f(.), i.e. whether to take function 
  ! values, derivatives or what else to get the function value f(u).
  !
  ! The matrix structure must be prepared with bilf_createMatrixStructure
  ! in advance.
  ! In case the array for the matrix entries does not exist, the routine
  ! allocates memory in size of the matrix of the heap for the matrix entries.
  !
  ! For setting up the entries, the discretisation structure attached to
  ! the matrix is used (rmatrixScalar%p_rdiscretisation). This is
  ! normally attached to the matrix by bilf_createMatrixStructure.
  !
  ! For the routine to work properly, it's important that the discretisation
  ! structure of rvector is 'compatible' with the discretisation structure
  ! of the matrix! I.e., the element distributions must be the same 
  ! (in number and ordering of the elements) except for the element type!
  !
  ! The matrix must be unsorted when this routine is called, 
  ! otherwise an error is thrown.
!</description>

!<input>

  ! The trilinear form specifying the underlying PDE of the discretisation.
  TYPE(t_trilinearForm), INTENT(IN) :: rform
  
  ! Whether to clear the matrix before calculating the entries.
  ! If .FALSE., the new matrix entries are added to the existing entries.
  LOGICAL, INTENT(IN) :: bclear
  
  ! A finite element function $u$ to be used as multiplier in the trilinear form.
  ! Must be a double precision vector.
  TYPE(t_vectorScalar), INTENT(IN) :: rvector

  ! OPTIONAL: A collection structure. This structure is given to the
  ! callback function for nonconstant coefficients to provide additional
  ! information. 
  TYPE(t_collection), INTENT(INOUT), TARGET, OPTIONAL :: rcollection
  
  ! OPTIONAL: A callback routine for nonconstant coefficient matrices.
  ! Must be present if the matrix has nonconstant coefficients!
  INCLUDE 'intf_coefficientTrilMatrixSc.inc'
  OPTIONAL :: fcoeff_buildTrilMatrixSc_sim
!</input>

!<inputoutput>
  ! The FE matrix. Calculated matrix entries are imposed to this matrix.
  TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrixScalar
!</inputoutput>

!</subroutine>

  ! The matrix must be unsorted, otherwise we can't set up the matrix.
  ! Note that we cannot switch off the sorting as easy as in the case
  ! of a vector, since there's a structure behind the matrix! So the caller
  ! has to make sure, the matrix is unsorted when this routine is called.
  IF (rmatrixScalar%isortStrategy .GT. 0) THEN
    PRINT *,'trilf_buildMatrixScalar: Matrix-structure must be unsorted!'
    CALL sys_halt()
  END IF

  IF (.NOT. ASSOCIATED(rmatrixScalar%p_rspatialDiscretisation)) THEN
    PRINT *,'trilf_buildMatrixScalar: No discretisation associated!'
    CALL sys_halt()
  END IF

  ! Do we have a uniform triangulation? Would simplify a lot...
  SELECT CASE (rmatrixScalar%p_rspatialDiscretisation%ccomplexity)
  CASE (SPDISC_UNIFORM) 
    ! Uniform discretisation; only one type of elements, e.g. P1 or Q1
    SELECT CASE (rmatrixScalar%cdataType)
    CASE (ST_DOUBLE) 
      ! Which matrix structure do we have?
      SELECT CASE (rmatrixScalar%cmatrixFormat) 
      CASE (LSYSSC_MATRIX9)
        !IF (PRESENT(fcoeff_buildMatrixSc_sim)) THEN
          CALL trilf_buildMatrix9d_conf2 (rform,bclear,rmatrixScalar,rvector,&  
                                          fcoeff_buildTrilMatrixSc_sim,rcollection)
      CASE (LSYSSC_MATRIX7)
        ! Convert structure 7 to structure 9.
        CALL lsyssc_convertMatrix (rmatrixScalar,LSYSSC_MATRIX9)
        
        ! Create the matrix in structure 9
        CALL trilf_buildMatrix9d_conf2 (rform,bclear,rmatrixScalar,rvector,&  
                                        fcoeff_buildTrilMatrixSc_sim,rcollection)
                                       
        ! Convert back to structure 7
        CALL lsyssc_convertMatrix (rmatrixScalar,LSYSSC_MATRIX7)
                                       
      CASE DEFAULT
        PRINT *,'trilf_buildMatrixScalar: Not supported matrix structure!'
        CALL sys_halt()
      END SELECT
    CASE DEFAULT
      PRINT *,'trilf_buildMatrixScalar: Single precision matrices currently not supported!'
      CALL sys_halt()
    END SELECT
    
  CASE (SPDISC_CONFORMAL) 
    
    ! Conformal discretisation; may have mixed P1/Q1 elements e.g.
    SELECT CASE (rmatrixScalar%cdataType)
    CASE (ST_DOUBLE) 
      ! Which matrix structure do we have?
      SELECT CASE (rmatrixScalar%cmatrixFormat) 
      CASE (LSYSSC_MATRIX9)
        !IF (PRESENT(fcoeff_buildMatrixSc_sim)) THEN
          CALL trilf_buildMatrix9d_conf2 (rform,bclear,rmatrixScalar,rvector,&  
                                          fcoeff_buildTrilMatrixSc_sim,rcollection)
        
      CASE (LSYSSC_MATRIX7)
        ! Convert structure 7 to structure 9
        CALL lsyssc_convertMatrix (rmatrixScalar,LSYSSC_MATRIX9)
        
        ! Create the matrix in structure 9
        CALL trilf_buildMatrix9d_conf2 (rform,bclear,rmatrixScalar,rvector,&  
                                        fcoeff_buildTrilMatrixSc_sim,rcollection)
                                       
        ! Convert back to structure 7
        CALL lsyssc_convertMatrix (rmatrixScalar,LSYSSC_MATRIX7)

      CASE DEFAULT
        PRINT *,'bilf_buildMatrixScalar: Not supported matrix structure!'
        CALL sys_halt()
      END SELECT
    CASE DEFAULT
      PRINT *,'bilf_buildMatrixScalar: Single precision matrices &
              &currently not supported!'
      CALL sys_halt()
    END SELECT
  CASE DEFAULT
    PRINT *,'bilf_buildMatrixScalar: General discretisation not implemented!'
    CALL sys_halt()
  END SELECT

  END SUBROUTINE
  
  !****************************************************************************

!<subroutine>

  SUBROUTINE trilf_buildMatrix9d_conf2 (rform,bclear,rmatrixScalar,rvector,&
                                        fcoeff_buildTrilMatrixSc_sim,rcollection)
  
!<description>
  ! This routine calculates the entries of a finite element matrix.
  ! The matrix structure must be prepared with bilf_createMatrixStructure
  ! in advance. The discretisation is assumed to be conformal, i.e. the DOF's
  ! of all finite elements must 'match'. Trial and test functions may be
  ! different.
  ! In case the array for the matrix entries does not exist, the routine
  ! allocates memory in size of the matrix of the heap for the matrix entries.
  !
  ! For setting up the entries, the discretisation structure attached to
  ! the matrix is used (rmatrixScalar%p_rdiscretisation). This is
  ! normally attached to the matrix by bilf_createMatrixStructure.
  !
  ! For the routine to work properly, it's important that the discretisation
  ! structure of rvector is 'compatible' with the discretisation structure
  ! of the matrix! I.e., the element distributions must be the same 
  ! (in number and ordering of the elements) except for the element type!
  !
  ! Double-precision version.
!</description>

!<input>
  ! The trilinear form specifying the underlying PDE of the discretisation.
  TYPE(t_trilinearForm), INTENT(IN) :: rform
  
  ! Whether to clear the matrix before calculating the entries.
  ! If .FALSE., the new matrix entries are added to the existing entries.
  LOGICAL, INTENT(IN) :: bclear
  
  ! A finite element function $u$ to be used as multiplier in the trilinear form.
  ! Must be a double precision vector.
  TYPE(t_vectorScalar), INTENT(IN) :: rvector

  ! OPTIONAL: A pointer to a collection structure. This structure is given to the
  ! callback function for nonconstant coefficients to provide additional
  ! information. 
  TYPE(t_collection), INTENT(INOUT), TARGET, OPTIONAL :: rcollection
  
  ! OPTIONAL: A callback routine for nonconstant coefficient matrices.
  ! Must be present if the matrix has nonconstant coefficients!
  INCLUDE 'intf_coefficientTrilMatrixSc.inc'
  OPTIONAL :: fcoeff_buildTrilMatrixSc_sim
!</input>

!<inputoutput>
  ! The FE matrix. Calculated matrix entries are imposed to this matrix.
  TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrixScalar
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i,i1,k,icurrentElementDistr,JDFG, ICUBP, IALBET, IA, IB, ifunc
  LOGICAL :: bIdenticalTrialAndTest, bIdenticalFuncAndTrial, bIdenticalFuncAndTest
  INTEGER(I32) :: IEL, IELmax, IELset, IDOFE, JDOFE
  INTEGER(PREC_DOFIDX) :: JCOL0,JCOL,idertype
  REAL(DP) :: OM,AUX, DB
  
  ! Array to tell the element which derivatives to calculate
  LOGICAL, DIMENSION(EL_MAXNDER) :: BderTrialTempl, BderTestTempl, BderFuncTempl
  LOGICAL, DIMENSION(EL_MAXNDER) :: BderTrial, BderTest, BderFunc
  
  ! Cubature point coordinates on the reference element
  REAL(DP), DIMENSION(CUB_MAXCUBP, NDIM3D) :: Dxi

  ! For every cubature point on the reference element,
  ! the corresponding cubature weight
  REAL(DP), DIMENSION(CUB_MAXCUBP) :: Domega
  
  ! number of cubature points on the reference element
  INTEGER :: ncubp
  
  ! Pointer to KLD, KCOL, DA
  INTEGER(I32), DIMENSION(:), POINTER :: p_KLD, p_KCOL
  REAL(DP), DIMENSION(:), POINTER :: p_DA
  
  ! An allocateable array accepting the DOF's of a set of elements.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsTest, IdofsTrial
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsFunc
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), POINTER :: p_IdofsTrial,p_IdofsFunc
  !INTEGER(PREC_DOFIDX), DIMENSION(EL_MAXNBAS,BILF_NELEMSIM), TARGET :: IdofsTest, IdofsTrial
  !INTEGER(PREC_DOFIDX), DIMENSION(:,:), POINTER :: p_IdofsTrial
  
  ! Allocateable arrays for the values of the basis functions - 
  ! for test and trial spaces.
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: DbasFunc,DbasTest,DbasTrial
  REAL(DP), DIMENSION(:,:,:,:), POINTER :: p_DbasFunc,p_DbasTrial
  
  ! Number of entries in the matrix - for quicker access
  INTEGER(PREC_DOFIDX) :: NA,NVE
  INTEGER(I32) :: NEQ
  
  ! Type of transformation from the reference to the real element 
  INTEGER :: ctrafoType
  
  ! Element evaluation tag; collects some information necessary for evaluating
  ! the elements.
  INTEGER(I32) :: cevaluationTag
  
  ! Number of local degees of freedom for trial and test functions
  INTEGER :: indofFunc, indofTrial, indofTest
  
  ! The triangulation structure - to shorten some things...
  TYPE(t_triangulation), POINTER :: p_rtriangulation
  
  ! A pointer to an element-number list
  INTEGER(I32), DIMENSION(:), POINTER :: p_IelementList
  
  ! Local matrices, used during the assembly.
  ! Values and positions of values in the global matrix.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:,:), ALLOCATABLE :: Kentry
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Dentry
  
  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: p_DcubPtsRef

  ! Pointer to the jacobian determinants
  REAL(DP), DIMENSION(:,:), POINTER :: p_Ddetj
  
  ! Current element distribution for discretisation and function $u$.
  TYPE(t_elementDistribution), POINTER :: p_elementDistribution
  TYPE(t_elementDistribution), POINTER :: p_elementDistributionFunc
  
  ! Number of elements in the current element distribution
  INTEGER(PREC_ELEMENTIDX) :: NEL
  
  ! DOF-Data of the vector
  REAL(DP), DIMENSION(:), POINTER :: p_Ddata
  
  ! Number of elements in a block. Normally =BILF_NELEMSIM,
  ! except if there are less elements in the discretisation.
  INTEGER :: nelementsPerBlock
  
  ! Some variables to support nonconstant coefficients in the matrix.
  
  ! Pointer to the coefficients that are computed by the callback routine.
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Dcoefficients
  
  ! A t_domainIntSubset structure that is used for storing information
  ! and passing it to callback routines.
  TYPE(t_domainIntSubset) :: rintSubset
  
  ! The discretisation - for easier access
  TYPE(t_spatialDiscretisation), POINTER :: p_rdiscretisation
  TYPE(t_spatialDiscretisation), POINTER :: p_rdiscretisationFunc
  
  IF (.NOT. ASSOCIATED(rmatrixScalar%p_rspatialDiscretisation)) THEN
    PRINT *,'trilf_buildMatrix9d_conf2: No discretisation associated!'
    CALL sys_halt()
  END IF

  ! Which derivatives of basis functions are needed?
  ! Check the descriptors of the bilinear form and set BDERxxxx
  ! according to these.

  BderTrialTempl = .FALSE.
  BderTestTempl = .FALSE.
  BderFuncTempl = .FALSE.
  
  ! Loop through the additive terms
  DO i=1,rform%itermCount
    ! The desriptor Idescriptors gives directly the derivative
    ! which is to be computed! Build template's for BDER.
    ! We don't compute the actual BDER here, as there might be some special
    ! processing if trial/test functions are identical!
    !
    ! At first build the descriptors for the trial functions
    I1=rform%Idescriptors(1,I)
    
    IF ((I1 .LT.0) .OR. (I1 .GT. DER_MAXNDER)) THEN
      PRINT *,'trilf_buildMatrix9d_conf2: Invalid descriptor'
      CALL sys_halt()
    ENDIF
    
    IF (I1 .NE. 0) BderFuncTempl(I1)=.TRUE.

    I1=rform%Idescriptors(2,I)
    
    IF ((I1 .LE.0) .OR. (I1 .GT. DER_MAXNDER)) THEN
      PRINT *,'trilf_buildMatrix9d_conf2: Invalid descriptor'
      CALL sys_halt()
    ENDIF
    
    BderTrialTempl(I1)=.TRUE.

    ! Then those of the test functions
    I1=rform%Idescriptors(3,I)
    
    IF ((I1 .LE.0) .OR. (I1 .GT. DER_MAXNDER)) THEN
      PRINT *,'trilf_buildMatrix9d_conf2: Invalid descriptor'
      CALL sys_halt()
    ENDIF
    
    BderTestTempl(I1)=.TRUE.
  END DO
  
  ! Get information about the matrix:
  NA = rmatrixScalar%NA
  NEQ = rmatrixScalar%NEQ
  
  ! We need KCOL/KLD of our matric
  CALL lsyssc_getbase_Kcol (rmatrixScalar,p_KCOL)
  CALL lsyssc_getbase_Kld (rmatrixScalar,p_KLD)
  
  ! Check if the matrix entries exist. If not, allocate the matrix.
  IF (rmatrixScalar%h_DA .EQ. ST_NOHANDLE) THEN

    ! Clear the entries in the matrix - we need to start with zero
    ! when assembling a new matrix!
    CALL storage_new1D ('trilf_buildMatrix9d_conf2', 'DA', &
                        NA, ST_DOUBLE, rmatrixScalar%h_DA, &
                        ST_NEWBLOCK_ZERO)
    CALL lsyssc_getbase_double (rmatrixScalar,p_DA)

  ELSE
  
    CALL lsyssc_getbase_double (rmatrixScalar,p_DA)

    ! If desired, clear the matrix before assembling.
    IF (bclear) THEN
      CALL lalg_clearVectorDble (p_DA)
    END IF
    
  END IF
  
  ! Get the discretisation(s)
  p_rdiscretisation => rmatrixScalar%p_rspatialDiscretisation
  p_rdiscretisationFunc => rvector%p_rspatialDiscretisation
  
  IF (.NOT. ASSOCIATED(p_rdiscretisation)) THEN
    PRINT *,'trilf_buildMatrix9d_conf2 error: No discretisation attached to the matrix!'
    CALL sys_halt()
  END IF
  
  IF (.NOT. ASSOCIATED(p_rdiscretisationFunc)) THEN
    PRINT *,'trilf_buildMatrix9d_conf2 error: No discretisation attached to the vector!'
    CALL sys_halt()
  END IF
  
  IF (p_rdiscretisation%inumFESpaces .NE. p_rdiscretisationFunc%inumFESpaces) THEN
    PRINT *,'trilf_buildMatrix9d_conf2 error: Discretisations not compatible!'
    CALL sys_halt()
  END IF
  
  ! Get a pointer to the triangulation - for easier access.
  p_rtriangulation => p_rdiscretisation%p_rtriangulation
  
  ! For saving some memory in smaller discretisations, we calculate
  ! the number of elements per block. For smaller triangulations,
  ! this is NEL. If there are too many elements, it's at most
  ! BILF_NELEMSIM. This is only used for allocating some arrays.
  nelementsPerBlock = MIN(TRILF_NELEMSIM,p_rtriangulation%NEL)
  
  ! Now loop over the different element distributions (=combinations
  ! of trial and test functions) in the discretisation.

  DO icurrentElementDistr = 1,p_rdiscretisation%inumFESpaces
  
    ! Activate the current element distribution(s)
    p_elementDistribution => &
      p_rdiscretisation%RelementDistribution(icurrentElementDistr)

    p_elementDistributionFunc => &
      p_rdiscretisationFunc%RelementDistribution(icurrentElementDistr)
  
    ! Cancel if this element distribution is empty.
    IF (p_elementDistribution%NEL .EQ. 0) CYCLE

    ! Get the number of local DOF's for trial and test functions
    indofFunc = elem_igetNDofLoc(p_elementDistributionFunc%itrialElement)
    indofTrial = elem_igetNDofLoc(p_elementDistribution%itrialElement)
    indofTest = elem_igetNDofLoc(p_elementDistribution%itestElement)
    
    ! Get the number of corner vertices of the element
    NVE = elem_igetNVE(p_elementDistribution%itrialElement)
    IF (NVE .NE. elem_igetNVE(p_elementDistribution%itestElement)) THEN
      PRINT *,'trilf_buildMatrix9d_conf2: element spaces incompatible!'
      CALL sys_halt()
    END IF
    IF (NVE .NE. elem_igetNVE(p_elementDistributionFunc%itrialElement)) THEN
      PRINT *,'trilf_buildMatrix9d_conf2: element spaces incompatible!'
      CALL sys_halt()
    END IF
    
    ! Get from the trial element space the type of coordinate system
    ! that is used there:
    ctrafoType = elem_igetTrafoType(p_elementDistribution%itrialElement)
    
    ! Allocate some memory to hold the cubature points on the reference element
    ALLOCATE(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),CUB_MAXCUBP))

    ! Initialise the cubature formula,
    ! Get cubature weights and point coordinates on the reference element
    CALL cub_getCubPoints(p_elementDistribution%ccubTypeBilForm, ncubp, Dxi, Domega)
    
    ! Reformat the cubature points; they are in the wrong shape!
    DO i=1,ncubp
      DO k=1,UBOUND(p_DcubPtsRef,1)
        p_DcubPtsRef(k,i) = Dxi(i,k)
      END DO
    END DO
    
    ! Open-MP-Extension: Open threads here.
    ! Each thread will allocate its own local memory...
    !
    !$OMP PARALLEL PRIVATE(rintSubset, p_DcubPtsRef,p_DcubPtsReal, &
    !$OMP   p_Ddetj,kentry, dentry,DbasTest,DbasTrial, &
    !$OMP   IdofsTest,IdofsTrial,p_DcubPtsTrial,cevaluationTag, &
    !$OMP   p_DcubPtsTest,Dcoefficients, bIdenticalTrialandTest, p_IdofsTrial, &
    !$OMP   p_DbasTrial, BderTrial, BderTest, ielmax,IEL, idofe,jdofe,jcol0, &
    !$OMP   jcol,JDFG,ICUBP, IALBET,OM,ia,ib,ifunc,aux, db,IdofsFunc,DbasFunc, &
    !$OMP   p_DcubPtsFunc,bIdenticalFuncAndTrial,p_IdofsFunc, &
    !$OMP   p_DbasFunc,iderType,p_IelementList,p_Ddata)    
    
    ! Quickly check if one of the specified derivatives is out of the allowed range:
    !$OMP SINGLE
    DO IALBET = 1,rform%itermcount
      ifunc = rform%Idescriptors(1,IALBET)
      IA = rform%Idescriptors(2,IALBET)
      IB = rform%Idescriptors(3,IALBET)      

      IF ((ifunc.LT.0) .OR. &
          (ifunc .GT. elem_getMaxDerivative(p_elementDistribution%itrialElement))) THEN
        PRINT *,'trilf_buildMatrix9d_conf2: Specified function-derivative',ifunc,&
                ' not available'
        CALL sys_halt()
      END IF
      
      IF ((IA.LE.0) .OR. &
          (IA .GT. elem_getMaxDerivative(p_elementDistribution%itrialElement))) THEN
        PRINT *,'trilf_buildMatrix9d_conf2: Specified trial-derivative',IA,&
                ' not available'
        CALL sys_halt()
      END IF

      IF ((IB.LE.0) .OR. &
          (IB .GT. elem_getMaxDerivative(p_elementDistribution%itestElement))) THEN
        PRINT *,'trilf_buildMatrix9d_conf2: Specified test-derivative',IB,&
                ' not available'
        CALL sys_halt()
      END IF
    END DO
    !$OMP END SINGLE
    
    ! Allocate an array saving the coordinates of corner vertices of elements
    
    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  ALLOCATE(DbasTest(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    !  ALLOCATE(DbasTrial(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    ! would lead to nonused memory blocks in these arrays during the assembly, 
    ! which reduces the speed by 50%!
    
    ALLOCATE(DbasTest(indofTest,elem_getMaxDerivative(p_elementDistribution%itestElement),&
             ncubp,nelementsPerBlock))
    ALLOCATE(DbasTrial(indofTrial,elem_getMaxDerivative(p_elementDistribution%itrialElement), &
             ncubp,nelementsPerBlock))
    ALLOCATE(DbasFunc(indofTrial,&
             elem_getMaxDerivative(p_elementDistributionFunc%itrialElement), &
             ncubp,nelementsPerBlock))

    ! Allocate memory for the DOF's of all the elements.
    ALLOCATE(IdofsTest(indofTest,nelementsPerBlock))
    ALLOCATE(IdofsTrial(indofTrial,nelementsPerBlock))
    ALLOCATE(IdofsFunc(indofFunc,nelementsPerBlock))

    ! Allocate an array saving the local matrices for all elements
    ! in an element set.
    ! We could also allocate EL_MAXNBAS*EL_MAXNBAS*BILF_NELEMSIM integers
    ! for this local matrix, but this would normally not fit to the cache
    ! anymore! indofTrial*indofTest*BILF_NELEMSIM is normally much smaller!
    ALLOCATE(Kentry(indofTrial,indofTest,nelementsPerBlock))
    ALLOCATE(Dentry(indofTrial,indofTest))
    
    ! Allocate an array for the coefficients c(x,y) f(u(x,y)). Because
    ! of the finite element function u included there, we have definitely
    ! no constant coefficients.
    ALLOCATE(Dcoefficients(rform%itermCount,ncubp,nelementsPerBlock))
                    
    ! p_IdofsTest points either to the just created array or to the
    ! array with the DOF's of the trial functions - when trial and
    ! test functions are identical.
    ! We don't rely on bidenticalTrialAndTest purely, as this does not
    ! indicate whether there are identical trial and test functions
    ! in one block!
    bIdenticalTrialAndTest = p_elementDistribution%itrialElement .EQ. &
                             p_elementDistribution%itestElement
                             
    ! Let p_IdofsTrial point either to IdofsTrial or to the DOF's of the test
    ! space IdofTest (if both spaces are identical). 
    ! We create a pointer for the trial space and not for the test space to
    ! prevent pointer-arithmetic in the innerst loop below!
    IF (bIdenticalTrialAndTest) THEN
      p_IdofsTrial => IdofsTest
      p_DbasTrial  => DbasTest
      ! Build the actual combination of what the element should calculate.
      ! As we evaluate only once, what the element must calculate is an
      ! OR combination of the BDER from trial and test functions.
      BderTrial = BderTrialTempl .OR. BderTestTempl
      BderTest = BderTrial
    ELSE
      p_IdofsTrial => IdofsTrial
      p_DbasTrial  => DbasTrial
      
      ! Build the actual combination of what the element should calculate.
      ! Copy BDERxxxx to BDERxxxxAct
      BderTrial = BderTrialTempl
      BderTest = BderTestTempl
    END IF
    
    ! Test whether the u trial functions are identical to the trial functions
    ! of the discretisation.
    bIdenticalFuncAndTest = p_elementDistribution%itestElement .EQ. &
                            p_elementDistributionFunc%itrialElement
    bIdenticalFuncAndTrial = p_elementDistribution%itrialElement .EQ. &
                             p_elementDistributionFunc%itrialElement
                             
    ! If yes, we can use the data calculated for the trial functions.
    IF (bIdenticalFuncAndTest) THEN
      p_IdofsFunc => IdofsTest
      p_DbasFunc  => DbasTest
      ! Build the actual combination of what the element should calculate.
      ! As we evaluate only once, what the element must calculate is an
      ! OR combination of the BDER from trial and test functions.
      BderTest = BderTest .OR. BderFuncTempl
    ELSE IF (bIdenticalFuncAndTrial) THEN
      p_IdofsFunc => p_IdofsTrial
      p_DbasFunc  => p_DbasTrial
      BderTrial = BderTrial .OR. BderFuncTempl
    ELSE
      p_IdofsFunc => IdofsFunc
      p_DbasFunc  => DbasFunc
      BderFunc = BderFuncTempl
    END IF

    ! p_IelementList must point to our set of elements in the discretisation
    ! with that combination of trial/test functions
    CALL storage_getbase_int (p_elementDistribution%h_IelementList, &
                              p_IelementList)
                         
    ! Get the data array from the vector
    CALL lsyssc_getbase_double(rvector,p_Ddata)
                              
    ! Get the number of elements there.
    NEL = p_elementDistribution%NEL

    ! Loop over the elements - blockwise.
    !
    ! Open-MP-Extension: Each loop cycle is executed in a different thread,
    ! so TRILF_NELEMSIM local matrices are simultaneously calculated in the
    ! inner loop(s).
    ! The blocks have all the same size, so we can use static scheduling.
    !
    !$OMP do schedule(static,1)
    DO IELset = 1, NEL, TRILF_NELEMSIM
    
      ! We always handle BILF_NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most BILF_NELEMSIM
      ! elements simultaneously.
      
      IELmax = MIN(NEL,IELset-1+TRILF_NELEMSIM)
    
      ! --------------------- DOF SEARCH PHASE ------------------------
    
      ! The outstanding feature with finite elements is: A basis
      ! function for a DOF on one element has common support only
      ! with the DOF's on the same element! E.g. for Q1:
      !
      !        #. . .#. . .#. . .#
      !        .     .     .     .
      !        .  *  .  *  .  *  .
      !        #-----O-----O. . .#
      !        |     |     |     .
      !        |     | IEL |  *  .
      !        #-----X-----O. . .#
      !        |     |     |     .
      !        |     |     |  *  .
      !        #-----#-----#. . .#
      !
      ! --> On element IEL, the basis function at "X" only interacts
      !     with the basis functions in "O". Elements in the 
      !     neighbourhood ("*") have no support, therefore we only have
      !     to collect all "O" DOF's.
      !
      ! Calculate the global DOF's into IdofsTrial / IdofsTest.
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF's of our BILF_NELEMSIM elements simultaneously.
      CALL dof_locGlobMapping_mult(p_rdiscretisation, p_IelementList(IELset:IELmax), &
                                   .TRUE.,IdofsTest)
                                   
      ! If the DOF's for the test functions are different, calculate them, too.
      IF (.NOT.bIdenticalTrialAndTest) THEN
        CALL dof_locGlobMapping_mult(p_rdiscretisation, p_IelementList(IELset:IELmax), &
                                     .FALSE.,IdofsTrial)
      END IF

      ! If the DOF's for the coefficient function values are different, 
      ! calculate them, too.
      IF ((.NOT. bIdenticalFuncAndTest) .AND. (.NOT. bIdenticalFuncAndTrial)) THEN
        CALL dof_locGlobMapping_mult(p_rdiscretisationFunc, p_IelementList(IELset:IELmax), &
                                     .FALSE.,IdofsFunc)
      END IF
      
      ! ------------------- LOCAL MATRIX SETUP PHASE -----------------------
      
      ! For the assembly of the global matrix, we use a "local"
      ! approach. At first we build a "local" system matrix according
      ! to the current element. This contains all additive
      ! contributions of element IEL, which are later added at the
      ! right positions to the elements in the global system matrix.
      !
      ! We have indofTrial trial DOF's per element and
      ! indofTest test DOF's per element. Therefore there are
      ! indofTrial*indofTest tupel of basis-/testfunctions (phi_i,psi_j) 
      ! "active" (i.e. have common support) on our current element, each 
      ! giving an additive contribution to the system matrix.
      !
      ! We build a quadratic indofTrial*indofTest local matrix:
      ! Kentry(1..indofTrial,1..indofTest) receives the position 
      !   in the global system matrix, where the corresponding value 
      !   has to be added to.
      ! (The corresponding contrbutions can be saved separately, 
      !  but we directly add them to the global matrix in this 
      !  approach.)
      !
      ! We build local matrices for all our elements 
      ! in the set simultaneously.
      ! Loop through elements in the set and for each element,
      ! loop through the local matrices to initialise them:
      DO IEL=1,IELmax-IELset+1
      
        ! For building the local matrices, we have first to
        ! loop through the test functions (the "O"'s), as these
        ! define the rows in the matrix.
        DO IDOFE=1,indofTest
        
          ! Row IDOFE of the local matrix corresponds 
          ! to row=global DOF KDFG(IDOFE) in the global matrix.
          ! This is one of the the "O"'s in the above picture.
          ! Get the starting position of the corresponding row
          ! to JCOL0:

          JCOL0=p_KLD(IdofsTest(IDOFE,IEL))
          
          ! Now we loop through the other DOF's on the current element
          ! (the "O"'s).
          ! All these have common support with our current basis function
          ! and will therefore give an additive value to the global
          ! matrix.
          
          DO JDOFE=1,indofTrial
            
            ! Get the global DOF of the "X" which interacts with 
            ! our "O".
            
            JDFG=p_IdofsTrial(JDOFE,IEL)
            
            ! Starting in JCOL0 (which points to the beginning of
            ! the line initially), loop through the elements in
            ! the row to find the position of column IDFG.
            ! Jump out of the DO loop if we find the column.
            
            DO JCOL=JCOL0,NA
              IF (p_KCOL(JCOL) .EQ. JDFG) EXIT
            END DO

            ! Because columns in the global matrix are sorted 
            ! ascendingly (except for the diagonal element),
            ! the next search can start after the column we just found.
            
            ! JCOL0=JCOL+1
            
            ! Save the position of the matrix entry into the local
            ! matrix.
            ! Note that a column in Kentry corresponds to a row in
            ! the real matrix. We aligned Kentry/DENTRY this way to get
            ! higher speed of the assembly routine, since this leads
            ! to better data locality.
            
            Kentry(JDOFE,IDOFE,IEL)=JCOL
            
          END DO ! IDOFE
          
        END DO ! JDOFE
        
      END DO ! IEL
      
      ! -------------------- ELEMENT EVALUATION PHASE ----------------------
      
      ! Ok, we found the positions of the local matrix entries
      ! that we have to change.
      ! To calculate the matrix contributions, we have to evaluate
      ! the elements to give us the values of the basis functions
      ! in all the DOF's in all the elements in our set.
      !
      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag. 
      cevaluationTag = elem_getEvaluationTag(p_elementDistribution%itrialElement)
      cevaluationTag = IOR(cevaluationTag,&
                      elem_getEvaluationTag(p_elementDistribution%itestElement))
                      
      IF (.NOT. rform%ballCoeffConstant) THEN
        ! Evaluate real coordinates if not necessary.
        cevaluationTag = IOR(cevaluationTag,EL_EVLTAG_REALCOORDS)
      END IF
      
      ! In the first loop, calculate the coordinates on the reference element.
      ! In all later loops, use the precalculated information.
      IF (IELset .EQ. 1) THEN
        cevaluationTag = IOR(cevaluationTag,EL_EVLTAG_REFCOORDS)
      ELSE
        cevaluationTag = IAND(cevaluationTag,NOT(EL_EVLTAG_REFCOORDS))
      END IF

      ! Calculate all information that is necessary to evaluate the finite element
      ! on all cells of our subset. This includes the coordinates of the points
      ! on the cells.
      CALL elprep_prepareSetForEvaluation (rintSubset%revalElementSet,&
          cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
          ctrafoType, p_DcubPtsRef(:,1:ncubp))
      p_Ddetj => rintSubset%revalElementSet%p_Ddetj

      ! If the matrix has nonconstant coefficients c(x,y) , calculate the 
      ! coefficients now.
      IF (.NOT. rform%ballCoeffConstant) THEN
        rintSubset%ielementDistribution = icurrentElementDistr
        rintSubset%ielementStartIdx = IELset
        rintSubset%p_Ielements => p_IelementList(IELset:IELmax)
        CALL fcoeff_buildMatrixSc_sim (p_rdiscretisation,rform, &
                  IELmax-IELset+1_I32,ncubp,&
                  rintSubset%revalElementSet%p_DpointsReal,&
                  p_IdofsTrial,IdofsTest,rintSubset, Dcoefficients,rcollection)
      END IF
      
      ! Calculate the values of the basis functions.
      CALL elem_generic_sim2 (p_elementDistribution%itestElement, &
          rintSubset%revalElementSet, BderTest, DbasTest)
      
      ! Omit the calculation of the trial function values if they
      ! are identical to the test function values.
      IF (.NOT. bidenticalTrialAndTest) THEN
        CALL elem_generic_sim2 (p_elementDistribution%itrialElement, &
            rintSubset%revalElementSet, BderTrial, DbasTrial)
      END IF
      
      ! Omit the calculation of the coefficient function values if they
      ! are identical to the trial function values.
      IF ((.NOT. bIdenticalFuncAndTest) .AND. (.NOT. bIdenticalFuncAndTrial)) THEN
        CALL elem_generic_sim2 (p_elementDistributionFunc%itrialElement, &
            rintSubset%revalElementSet, BderFunc, DbasFunc)
      END IF
      
      ! ----------------- COEFFICIENT EVALUATION PHASE ---------------------
      
      ! Now its time to form the actual coefficients for the integral:
      ! Dcoefficients = c(x,y) f(u(x,y)). So we must evaluate f(u(x,y))
      ! and multiply it with the coefficients in Dcoefficient to get the
      ! correct coefficients for later use.
      IF (rform%ballCoeffConstant) THEN
        ! Constant coefficients. Take the coefficients from the bilinear form
        ! and multiply with the values of f(u).
        DO IALBET = 1,rform%itermcount
          iderType = rform%Idescriptors(1,IALBET)
          IF (iderType .NE. 0) THEN
            DO iel=1,IELmax-IELset+1
              DO ICUBP = 1,ncubp
                ! Calculate the value in the point
                DB = 0.0_DP
                DO IDOFE = 1,indofTrial
                  DB = DB + &
                    p_Ddata(p_IdofsFunc(IDOFE,iel)) * p_DbasFunc(IDOFE,iderType,ICUBP,iel)
                END DO
                ! Save the value in the point, multiplied with the coefficient
                Dcoefficients(IALBET,ICUBP,iel) = rform%Dcoefficients(IALBET) * DB
              END DO
            END DO
          ELSE
            Dcoefficients(IALBET,1:ncubp,1:IELmax-IELset+1) = rform%Dcoefficients(IALBET)
          END IF
        END DO
      ELSE
        ! Nonconstant coefficients. Take the calculated coefficients in Dcoefficients
        ! and multiply with the values of f(u).      
        DO IALBET = 1,rform%itermcount
          iderType = rform%Idescriptors(1,IALBET)
          IF (iderType .NE. 0) THEN
            DO iel=1,IELmax-IELset+1
              DO ICUBP = 1,ncubp
                ! Calculate the value in the point
                DB = 0.0_DP
                DO IDOFE = 1,indofTrial
                  DB = DB + &
                    p_Ddata(p_IdofsFunc(IDOFE,iel)) * p_DbasFunc(IDOFE,iderType,ICUBP,iel)
                END DO
                ! Save the value in the point, multiplied with the existing coefficient
                Dcoefficients(IALBET,ICUBP,iel) = Dcoefficients(IALBET,ICUBP,iel) * DB
              END DO
            END DO
          END IF
        END DO
      END IF
      
      ! --------------------- DOF COMBINATION PHASE ------------------------
      
      ! Values of all basis functions calculated. Now we can start 
      ! to integrate! As we have never the case of constant coefficients
      ! (because of the FE function u involved), we have only the 'complex'
      ! loop here in comparison to the standard bilinear form.
      !
      ! Loop over the elements in the current set.
      DO IEL=1,IELmax-IELset+1
        
        ! Clear the local matrix
        Dentry = 0.0_DP
        ! Loop over all cubature points on the current element
        DO ICUBP = 1, ncubp

          ! calculate the current weighting factor in the cubature formula
          ! in that cubature point.
          !
          ! Take the absolut value of the determinant of the mapping.
          ! In 2D, the determinant is always positive, whereas in 3D,
          ! the determinant might be negative -- that's normal!

          OM = Domega(ICUBP)*ABS(p_Ddetj(ICUBP,IEL))

          ! Loop over the additive factors in the bilinear form.
          DO IALBET = 1,rform%itermcount
          
            ! Get from Idescriptors the type of the derivatives for the 
            ! test and trial functions. The summand we calculate
            ! here will be added to the matrix entry:
            !
            ! a_ij  =  int_... ( psi_j )_IA  *  ( phi_i )_IB
            !
            ! -> Ix=0: function value, 
            !      =1: first derivative, ...
            !    as defined in the module 'derivative'.
            
            IA = rform%Idescriptors(2,IALBET)
            IB = rform%Idescriptors(3,IALBET)
            
            ! Multiply OM with the coefficient of the form.
            ! This gives the actual value to multiply the
            ! function value with before summing up to the integral.
            ! Get the precalculated coefficient from the coefficient array.
            AUX = OM * Dcoefficients(IALBET,ICUBP,IEL)
            
            ! Now loop through all possible combinations of DOF's
            ! in the current cubature point. The outer loop
            ! loops through the "O" in the above picture,
            ! the test functions:

            DO IDOFE=1,indofTest
              
              ! Get the value of the (test) basis function 
              ! phi_i (our "O") in the cubature point:
              DB = DbasTest(IDOFE,IB,ICUBP,IEL)
              
              ! Perform an inner loop through the other DOF's
              ! (the "X"). 

              DO JDOFE=1,indofTrial
            
                ! Get the value of the basis function 
                ! psi_j (our "X") in the cubature point. 
                ! Them multiply:
                !    DB * DBAS(..) * AUX
                ! ~= phi_i * psi_j * coefficient * cub.weight
                ! Summing this up gives the integral, so the contribution
                ! to the global matrix. 
                !
                ! Simply summing up DB * DBAS(..) * AUX would give
                ! the coefficient of the local matrix. We save this
                ! contribution in the local matrix of element IEL.

                !JCOLB = Kentry(JDOFE,IDOFE,IEL)
                !p_DA(JCOLB) = p_DA(JCOLB) + DB*p_DbasTrial(JDOFE,IA,ICUBP,IEL)*AUX
                Dentry(JDOFE,IDOFE) = Dentry(JDOFE,IDOFE)+DB*p_DbasTrial(JDOFE,IA,ICUBP,IEL)*AUX
              
              END DO
            
            END DO ! JDOFE
            
          END DO ! IALBET

        END DO ! ICUBP 
        
        ! Incorporate the local matrices into the global one.
        ! Kentry gives the position of the additive contributions in Dentry.
        !$OMP CRITICAL
        DO IDOFE=1,indofTest
          DO JDOFE=1,indofTrial
            p_DA(Kentry(JDOFE,IDOFE,IEL)) = p_DA(Kentry(JDOFE,IDOFE,IEL)) + Dentry(JDOFE,IDOFE)
          END DO
        END DO
        !$OMP END CRITICAL

      END DO ! IEL

    END DO ! IELset
    !$OMP END DO
    
    ! Release memory
    CALL elprep_releaseElementSet(rintSubset%revalElementSet)

    DEALLOCATE(p_DcubPtsRef)
    DEALLOCATE(Dcoefficients)
    DEALLOCATE(IdofsTrial)
    DEALLOCATE(IdofsTest)
    DEALLOCATE(DbasTrial)
    DEALLOCATE(DbasTest)
    DEALLOCATE(Kentry)
    DEALLOCATE(Dentry)
    
    !$OMP END PARALLEL

  END DO ! icurrentElementDistr

  ! Finish
  
  END SUBROUTINE
  
END MODULE
