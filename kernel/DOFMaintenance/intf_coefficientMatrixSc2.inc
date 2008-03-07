  ! This file describes the interface for a coefficient function that is called
  ! during the matrix assembly of scalar matrices.
  
  ! WARNING: THIS IS THE OLD STYLE. WILL BE SWITCHED OFF IN NEAR FUTURE!
  
  INTERFACE
  
  !<subroutine>

    SUBROUTINE fcoeff_buildMatrixSc_sim (rdiscretisation,ielementDistribution, rform, &
                  ielementStartIdx,nelements,npointsPerElement,Ielements,Dcoords, &
                  DcubPtsRef,DcubPtsReal,IdofsTrial,IdofsTest,Djac,Ddetj, &
                  Dcoefficients,rcollection)
    
    USE basicgeometry
    USE triangulation
    USE collection
    USE scalarpde
    
  !<description>
    ! This subroutine is called during the matrix assembly. It has to compute
    ! the coefficients in front of the terms of the bilinear form.
    !
    ! The routine accepts a set of elements and a set of points on these
    ! elements (cubature points), in reference as well as in real coordinates.
    ! According to the terms in the bilinear form, the routine has to compute
    ! simultaneously for all these points and all the terms in the bilinear form
    ! the corresponding coefficients in front of the terms.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    TYPE(t_spatialDiscretisation), INTENT(IN)                   :: rdiscretisation
    
    ! The currently active element distribution in the discretisation.
    ! Allows the routine to get the currently active element type for
    ! trial and test functions.
    INTEGER, INTENT(IN)                                         :: ielementDistribution
    
    ! The bilinear form which is currently to be evaluated:
    TYPE(t_bilinearForm), INTENT(IN)                            :: rform
    
    ! Start index of the current element block Ielements in the current element 
    ! distribution ielementDistribution. If this is =1, the routine is called the 
    ! first time for the current element distribution.
    INTEGER(I32), INTENT(IN)                                    :: ielementStartIdx

    ! Number of elements, where the coefficients must be computed.
    ! This is always a part of the element distribution.
    INTEGER, INTENT(IN)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    INTEGER, INTENT(IN)                                         :: npointsPerElement
    
    ! A list of elements of length nelements where coefficients must
    ! be computed by this routine.
    INTEGER(I32), DIMENSION(nelements), INTENT(IN)              :: Ielements
    
    ! A list of the corner vertices of all elements in Ielements
    REAL(DP), DIMENSION(NDIM2D,TRIA_MAXNVE2D,nelements), INTENT(IN)  :: Dcoords
    
    ! A list of points in coordinates on the reference element.
    ! Each set of points corresponds to the corresponding element
    ! in Ielements
    REAL(DP), DIMENSION(NDIM2D,npointsPerElement,nelements), INTENT(IN) :: DcubPtsRef

    ! A list of points, corresponding to DcubPtsRef, in real coordinates.
    ! Each set of points corresponds to the corresponding element
    ! in Ielements.
    REAL(DP), DIMENSION(NDIM2D,npointsPerElement,nelements), INTENT(IN) :: DcubPtsReal
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in trial space,nelements)
    INTEGER(PREC_DOFIDX), DIMENSION(:,:) :: IdofsTrial
    
    ! An array accepting the DOF's on all elements trial in the trial space.
    ! DIMENSION(#local DOF's in test space,nelements)
    INTEGER(PREC_DOFIDX), DIMENSION(:,:) :: IdofsTest
    
    ! The Jacobian matrix of the mapping between the reference and each
    ! real element, for all points on all elements.
    ! DIMENSION(#entries in the matrix,npointsPerEl,nelements)
    REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: Djac
    
    ! The Jacobian determinant of the mapping of each point from the
    ! reference element to each real element
    REAL(DP), DIMENSION(npointsPerElement,nelements), INTENT(IN) :: Ddetj
    
    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    TYPE(t_collection), INTENT(INOUT), OPTIONAL      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    REAL(DP), DIMENSION(:,:,:), INTENT(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    END SUBROUTINE
    
  END INTERFACE
  