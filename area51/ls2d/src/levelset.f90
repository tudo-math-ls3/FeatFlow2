MODULE levelset

  USE fsystem
  USE linearsystemscalar
  USE linearsystemblock
  USE spatialdiscretisation
  USE scalarpde
  USE derivatives
  USE cubature
  USE collection
  USE domainintegration
  USE element
  USE elementpreprocessing
  USE multilevelprojection
  use filtersupport
  
  USE ccbasic
  
  IMPLICIT NONE

 
  
  CONTAINS
  
! ***************************************************************************
 SUBROUTINE initLS_analytical (rproblem,rvector)
 
   type(t_problem), intent(INOUT) :: rproblem
   type(t_vectorBlock), intent(INOUT) :: rvector
   
   ! local variables
   type(t_vectorBlock) :: rvectortmp,rvector2
   integer :: ierror
   type(t_linearForm) :: rlinform
   type(t_blockDiscretisation), pointer :: p_rdiscretisation
   type(t_linsolNode), pointer :: p_rsolverNode,p_rpreconditioner
   type(t_matrixBlock), dimension(1) :: Rmatrices
   TYPE(t_filterChain), DIMENSION(1), TARGET :: RfilterChain
   TYPE(t_filterChain), DIMENSION(:), POINTER :: p_RfilterChain
      
      p_rdiscretisation => rproblem%RlevelInfo(rproblem%NLMAX)%rdiscretisationLS
      call lsysbl_createVecBlockByDiscr (p_rdiscretisation,rvectortmp,.true.)
   
      rlinform%itermCount = 1
      rlinform%Idescriptors(1) = DER_FUNC
   
      call linf_buildVectorScalar (&
                p_rdiscretisation%RspatialDiscr(1),rlinform,.true.,&
                rvectortmp%RvectorBlock(1),coeff_AnalyticSolution_LS)
                
      !call linsol_initJacobi (p_rpreconditioner, 0.8_DP)
      !call linsol_initDefCorr (p_rsolverNode, p_rpreconditioner)
      
      !p_rsolverNode%depsRel = SYS_EPSREAL_DP * 100.0_DP
      !p_rsolverNode%nmaxiterations = 1000
      !p_rsolverNode%ioutputLevel = 1
     
      
      
      RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL
      p_RfilterChain => RfilterChain
      NULLIFY(p_rpreconditioner)
      CALL linsol_initUMFPACK4 (p_rsolverNode)
      p_rsolverNode%ioutputLevel = 2
      call lsysbl_createMatFromScalar(rproblem%rlevelinfo(rproblem%NLMAX)%rmatrixLS,Rmatrices(1))
      
      call linsol_setMatrices (p_rsolverNode,Rmatrices)
      call linsol_initStructure (p_rsolverNode,ierror)
      call linsol_initData (p_rsolverNode,ierror)
      
      call lsysbl_duplicateVector (rvector,rvector2,&
          LSYSSC_DUP_COPY,LSYSSC_DUP_EMPTY)
       
      call linsol_solveAdaptively (p_rsolverNode,rvector,rvectortmp,rvector2)
      
 
      call lsysbl_releaseVector (rvector2)
      call lsysbl_releaseVector (rvectortmp)
      call lsysbl_releaseMatrix (Rmatrices(1))
      
      call linsol_doneData (p_rsolverNode,ierror)
      call linsol_doneStructure (p_rsolverNode,ierror)
      call linsol_releaseSolver (p_rsolverNode)
                

 end subroutine
! *************************************************************************
  subroutine coeff_AnalyticSolution_ls (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! This routine is called upon program start if ctypeInitialSolution=3.
    ! It returns analytical values for the X-velocity in the
    ! initial solution vector.
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
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN)  :: Dpoints

    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(\#local DOF's in test space,nelements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It's usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
   integer :: i, j
   real(DP) :: dr
   
   
   do i = 1, nelements
      do j = 1, npointsPerElement
    
         dr = SQRT((Dpoints(1,j,i)-0.1_DP)**2+(Dpoints(2,j,i)-0.205)**2)
         Dcoefficients(:,j,i) = 0.05_DP-dr
         
     end do
   end do

  end subroutine



! ***************************************************************************
  SUBROUTINE prjF2C(rvector,rdiscretisation)
 
    TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
    !TYPE(t_problem_lvl), INTENT(IN),TARGET :: rlevelInfo
    type(t_blockDiscretisation),intent(in):: rdiscretisation
    
    !local
    TYPE(t_interlevelProjectionBlock) :: rprojection
    TYPE(t_vectorBlock):: rtmpvector
    TYPE(t_vectorScalar) :: rvectorTemp
    INTEGER:: NEQ

    
    CALL lsysbl_createVecBlockByDiscr (&
            rdiscretisation,rtmpvector,.false.)
        
          
    CALL mlprj_initProjectionVec (rprojection,rtmpVector)
    
    NEQ = mlprj_getTempMemoryVec (rprojection,rtmpVector,rvector)
   
    IF (NEQ .NE. 0) CALL lsyssc_createVector (rvectorTemp,NEQ,.FALSE.)
   
    CALL mlprj_performInterpolation (rprojection,rtmpvector,rvector, &
                                       rvectorTemp)
       
    IF (NEQ .NE. 0) CALL lsyssc_releaseVector (rvectorTemp)
    CALL lsysbl_swapVectors (rvector,rtmpvector)
    CALL lsysbl_releaseVector (rtmpvector)
    CALL mlprj_doneProjection (rprojection)
  
  END SUBROUTINE
  
   

END MODULE
