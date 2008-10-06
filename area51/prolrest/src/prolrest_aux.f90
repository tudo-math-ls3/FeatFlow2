!##############################################################################
!# ****************************************************************************
!# <name> prolrest_aux </name>
!# ****************************************************************************
!# <purpose>
!# This module contains a set of auxiliary routines which are used by the
!# prol/rest examples as well as some callback functions for the assembly
!# of the RHS vectors.
!# </purpose>
!##############################################################################

MODULE prolrest_aux

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

  ! ***************************************************************************

!<subroutine>

  PURE SUBROUTINE vec_filterByEps(Dvector, deps)
  
!<description>
  ! This routine filters all vector entries which are below a specified
  ! threshold, i.e. sets vector entries to zero.
!</description>

!<input>
  ! OPTIONAL: A threshold parameter. All vector entries whose absolute value
  ! drops below the threshold will be set to zero. If not specified,
  ! 100*SYS_EPSREAL is used.
  REAL(DP), OPTIONAL, INTENT(IN) :: deps
!</input>

!<inputoutput>
  ! The vector that is to be filtered.
  REAL(DP), DIMENSION(:), INTENT(INOUT) :: Dvector
!</inputoutput>

!<subroutine>

  ! local variables
  INTEGER :: i
  REAL(DP) :: dthreshold
  
    ! Is eps specified?
    IF (PRESENT(deps)) THEN
      dthreshold = deps
    ELSE
      dthreshold = 100.0_DP * SYS_EPSREAL
    END IF
    
    ! Do through the vector
    DO i = LBOUND(Dvector,1), UBOUND(Dvector,1)
      IF(ABS(Dvector(i)) .LE. dthreshold) Dvector(i) = 0.0_DP
    END DO
    
    ! That's it
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  
  SUBROUTINE mat_densify(Da, rmatrix)

!<description>
  ! This routine densifies a type 7/9 sparse matrix, i.e. copies the matrix
  ! entries into a dense matrix.
!</description>

!<input>
  ! The scalar matrix which is to be densified.
  type(t_matrixScalar), intent(IN) :: rmatrix
!</input>

!<output>
  ! The dense ouput matrix. The matrix is assumed to have correct dimensions.
  real(DP), dimension(:,:), intent(OUT) :: Da
!</output>

!</subroutine>

  ! local variables
  real(DP), dimension(:), pointer :: p_DA
  integer, dimension(:), pointer :: p_Kld, p_Kcol
  integer :: i,j
  
    ! Get the pointers from the storage
    call storage_getbase_double(rmatrix%h_DA,p_DA)
    call storage_getbase_int(rmatrix%h_Kld, p_Kld)
    call storage_getbase_int(rmatrix%h_Kcol, p_Kcol)
    
    ! Format output matrix
    Da = 0.0_DP
    
    ! densify
    do i = 1, rmatrix%NEQ
      do j = p_Kld(i), p_Kld(i+1)-1
        Da(i,p_Kcol(j)) = p_DA(j)
      end do ! j
    end do ! i
  
  END subroutine

  ! ***************************************************************************

!<subroutine>
  
  pure subroutine mat_identity(Da)
  
!<description>
  ! This routine formats a dense matrix to identity.
!</description>

!<inputoutput>
  ! The matrix which should recieve the identity matrix.
  real(DP), dimension(:,:), intent(INOUT) :: Da
!</inputoutput>

!<subroutine>
  
  integer :: i, n
  
    Da = 0.0_DP
    n = min(ubound(Da,1),ubound(Da,2))
    do i = 1, n
      Da(i,i) = 1.0_DP
    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  PURE SUBROUTINE mat_filterByEps(Da, deps)
  
!<description>
  ! This routine filters all matrix entries which are below a specified
  ! threshold, i.e. sets matrix entries to zero.
!</description>

!<input>
  ! OPTIONAL: A threshold parameter. All matrix entries whose absolute value
  ! drops below the threshold will be set to zero. If not specified,
  ! 100*SYS_EPSREAL is used.
  REAL(DP), OPTIONAL, INTENT(IN) :: deps
!</input>

!<inputoutput>
  ! The matrix that is to be filtered.
  REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: Da
!</inputoutput>

!<subroutine>

  ! local variables
  INTEGER :: i,j
  REAL(DP) :: dthreshold
  
    ! Is eps specified?
    IF (PRESENT(deps)) THEN
      dthreshold = deps
    ELSE
      dthreshold = 100.0_DP * SYS_EPSREAL
    END IF
    
    ! Do through the vector
    DO i = LBOUND(Da,1), UBOUND(Da,1)
      DO j = LBOUND(Da,2), UBOUND(Da,2)
        IF(ABS(Da(i,j)) .LE. dthreshold) Da(i,j) = 0.0_DP
      END DO
    END DO
    
    ! That's it
  
  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  subroutine prolrest_buildProlMatrix(rdiscrC, rdiscrF, Dprol)
  
!<description>
  ! Builds the local prolongation matrix based on an inter-level projection
  ! structure defined in the multilevelprojection module.
!</description>

!<input>
  ! The coarse mesh discretisation.
  type(t_blockDiscretisation), intent(IN) :: rdiscrC

  ! The fine mesh discretisation.
  type(t_blockDiscretisation), intent(IN) :: rdiscrF
!</input>

!<output>
  ! The local prolongation matrix.
  real(DP), dimension(:,:), intent(OUT) :: Dprol
!</output>

!</subroutine>

  ! Local variables
  type(t_vectorBlock) :: rvecC, rvecF
  type(t_vectorScalar) :: rtmp
  TYPE(t_interlevelProjectionBlock) :: rproj
  real(DP), dimension(:), pointer :: p_DvecC, p_DvecF
  integer :: i,j,m,n
  
    ! Create the projection structure
    CALL mlprj_initProjectionDiscr(rproj, rdiscrF)
    
    ! Create two block vectors
    CALL lsysbl_createVecBlockByDiscr(rdiscrC, rvecC)
    CALL lsysbl_createVecBlockByDiscr(rdiscrF, rvecF)
    
    ! Get the data arrays
    call storage_getbase_double(rvecC%h_Ddata, p_DvecC)
    call storage_getbase_double(rvecF%h_Ddata, p_DvecF)
    
    ! Get the dimensions
    m = rvecC%NEQ
    n = rvecF%NEQ
    
    ! Create temporary vector
    i = mlprj_getTempMemoryVec (rproj,rvecC,rvecF)
    CALL lsyssc_createVector (rtmp,i,.FALSE.)
    
    ! Calculate prolongation matrix
    Dprol = 0.0_DP
    
    ! Go through all coarse mesh DOFs
    do i = 1, m
    
      ! Set coarse vector
      p_DvecC = 0.0_DP
      p_DvecC(i) = 1.0_DP
      
      ! Prolongate vector
      CALL mlprj_performProlongation(rproj, rvecC, rvecF, rtmp)
      
      ! Store result
      do j = 1, n
        Dprol(j,i) = p_DvecF(j)
      end do
    end do

    ! Release temporary vector
    IF(rtmp%NEQ .gt. 0) then
      CALL lsyssc_releaseVector(rtmp)
    end if
    
    ! Release the vectors
    CALL lsysbl_releaseVector(rvecF)
    CALL lsysbl_releaseVector(rvecC)
    
    ! Release the projection structure
    CALL mlprj_doneProjection(rproj)
    
    ! That's it

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine prolrest_buildRestMatrix(rdiscrC, rdiscrF, Drest)
  
!<description>
  ! Builds the local restriction matrix based on an inter-level projection
  ! structure defined in the multilevelprojection module.
!</description>

!<input>
  ! The coarse mesh discretisation.
  type(t_blockDiscretisation), intent(IN) :: rdiscrC

  ! The fine mesh discretisation.
  type(t_blockDiscretisation), intent(IN) :: rdiscrF
!</input>

!<output>
  ! The local restriction matrix.
  real(DP), dimension(:,:), intent(OUT) :: Drest
!</output>

!</subroutine>

  ! Local variables
  type(t_vectorBlock) :: rvecC, rvecF
  type(t_vectorScalar) :: rtmp
  TYPE(t_interlevelProjectionBlock) :: rproj
  real(DP), dimension(:), pointer :: p_DvecC, p_DvecF
  integer :: i,j,m,n
  
    ! Create the projection structure
    CALL mlprj_initProjectionDiscr(rproj, rdiscrF)
    
    ! Create two block vectors
    CALL lsysbl_createVecBlockByDiscr(rdiscrC, rvecC)
    CALL lsysbl_createVecBlockByDiscr(rdiscrF, rvecF)
    
    ! Get the data arrays
    call storage_getbase_double(rvecC%h_Ddata, p_DvecC)
    call storage_getbase_double(rvecF%h_Ddata, p_DvecF)
    
    ! Get the dimensions
    m = rvecC%NEQ
    n = rvecF%NEQ
    
    ! Create temporary vector
    i = mlprj_getTempMemoryVec (rproj,rvecC,rvecF)
    CALL lsyssc_createVector (rtmp,i,.FALSE.)

    ! Calculate restriction matrix
    Drest = 0.0_DP
    
    ! Go through all fine mesh DOFs
    do i = 1, n
    
      ! Set fine vector
      p_DvecF = 0.0_DP
      p_DvecF(i) = 1.0_DP
      
      ! restrict vector
      CALL mlprj_performRestriction(rproj, rvecC, rvecF, rtmp)
      
      ! Store result
      do j = 1, m
        Drest(j,i) = p_DvecC(j)
      end do
    end do
      
    ! Release temporary vector
    IF(rtmp%NEQ .gt. 0) then
      CALL lsyssc_releaseVector(rtmp)
    end if
    
    ! Release the vectors
    CALL lsysbl_releaseVector(rvecF)
    CALL lsysbl_releaseVector(rvecC)
    
    ! Release the projection structure
    CALL mlprj_doneProjection(rproj)
    
    ! That's it

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine prolrest_buildInterpMatrix(rdiscrC, rdiscrF, Dinterp)
  
!<description>
  ! Builds the local interpolation matrix based on an inter-level projection
  ! structure defined in the multilevelprojection module.
!</description>

!<input>
  ! The coarse mesh discretisation.
  type(t_blockDiscretisation), intent(IN) :: rdiscrC

  ! The fine mesh discretisation.
  type(t_blockDiscretisation), intent(IN) :: rdiscrF
!</input>

!<output>
  ! The local interpolation matrix.
  real(DP), dimension(:,:), intent(OUT) :: Dinterp
!</output>

!</subroutine>

  ! Local variables
  type(t_vectorBlock) :: rvecC, rvecF
  type(t_vectorScalar) :: rtmp
  TYPE(t_interlevelProjectionBlock) :: rproj
  real(DP), dimension(:), pointer :: p_DvecC, p_DvecF
  integer :: i,j,m,n
  
    ! Create the projection structure
    CALL mlprj_initProjectionDiscr(rproj, rdiscrF)
    
    ! Create two block vectors
    CALL lsysbl_createVecBlockByDiscr(rdiscrC, rvecC)
    CALL lsysbl_createVecBlockByDiscr(rdiscrF, rvecF)
    
    ! Get the data arrays
    call storage_getbase_double(rvecC%h_Ddata, p_DvecC)
    call storage_getbase_double(rvecF%h_Ddata, p_DvecF)
    
    ! Get the dimensions
    m = rvecC%NEQ
    n = rvecF%NEQ
    
    ! Create temporary vector
    i = mlprj_getTempMemoryVec (rproj,rvecC,rvecF)
    CALL lsyssc_createVector (rtmp,i,.FALSE.)

    ! Calculate interpolation matrix
    Dinterp = 0.0_DP
    
    ! Go through all fine mesh DOFs
    do i = 1, n
    
      ! Set fine vector
      p_DvecF = 0.0_DP
      p_DvecF(i) = 1.0_DP
      
      ! interpolate vector
      CALL mlprj_performInterpolation(rproj, rvecC, rvecF, rtmp)
      
      ! Store result
      do j = 1, m
        Dinterp(j,i) = p_DvecC(j)
      end do
    end do

    ! Release temporary vector
    IF(rtmp%NEQ .gt. 0) then
      CALL lsyssc_releaseVector(rtmp)
    end if
    
    ! Release the vectors
    CALL lsysbl_releaseVector(rvecF)
    CALL lsysbl_releaseVector(rvecC)
    
    ! Release the projection structure
    CALL mlprj_doneProjection(rproj)
    
    ! That's it

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  SUBROUTINE procRHS_One (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
    USE basicgeometry
    USE triangulation
    USE collection
    USE scalarpde
    USE domainintegration
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
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
    INTEGER(PREC_ELEMENTIDX), INTENT(IN)                        :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    INTEGER, INTENT(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    REAL(DP), DIMENSION(:,:,:), INTENT(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    TYPE(t_domainIntSubset), INTENT(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    TYPE(t_collection), INTENT(INOUT), OPTIONAL      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    REAL(DP), DIMENSION(:,:,:), INTENT(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    !    u(x,y) = 1
    Dcoefficients (1,:,:) = 1.0_DP

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE procRHS_Q2_2D (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
    USE basicgeometry
    USE triangulation
    USE collection
    USE scalarpde
    USE domainintegration
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
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
    INTEGER(PREC_ELEMENTIDX), INTENT(IN)                        :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    INTEGER, INTENT(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    REAL(DP), DIMENSION(:,:,:), INTENT(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    TYPE(t_domainIntSubset), INTENT(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    TYPE(t_collection), INTENT(INOUT), OPTIONAL      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    REAL(DP), DIMENSION(:,:,:), INTENT(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    !    u(x,y) = 16*x*(1-x)*y*(1-y)
    Dcoefficients (1,:,:) = 16.0_DP * Dpoints(1,:,:)*(1.0_DP-Dpoints(1,:,:)) * &
                              Dpoints(2,:,:)*(1.0_DP-Dpoints(2,:,:))

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE procRHS_Sin_2D (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
    USE basicgeometry
    USE triangulation
    USE collection
    USE scalarpde
    USE domainintegration
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points) in real coordinates.
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
    INTEGER(PREC_ELEMENTIDX), INTENT(IN)                        :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    INTEGER, INTENT(IN)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    REAL(DP), DIMENSION(:,:,:), INTENT(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    TYPE(t_domainIntSubset), INTENT(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    TYPE(t_collection), INTENT(INOUT), OPTIONAL      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    REAL(DP), DIMENSION(:,:,:), INTENT(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>

    !    u(x,y) = sin(pi*x)*(sin(pi*y)
    Dcoefficients (1,:,:) = SIN(SYS_PI * Dpoints(1,:,:))*SIN(SYS_PI * Dpoints(2,:,:))

  END SUBROUTINE

END MODULE
