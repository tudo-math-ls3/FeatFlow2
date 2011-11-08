  ! -*- mode: f90; -*-
  
  ! This file describes the interface for a coefficient function that is called
  ! during the matrix assembly of scalar matrices.
  
  ! WARNING: THIS IS THE OLD STYLE. WILL BE SWITCHED OFF IN NEAR FUTURE!
  
  interface
  
  !<subroutine>

    subroutine fcoeff_buildMatrixSc_sim (rdiscretisation,ielementDistribution, rform, &
                  ielementStartIdx,nelements,npointsPerElement,Ielements,Dcoords, &
                  DcubPtsRef,DcubPtsReal,IdofsTrial,IdofsTest,Djac,Ddetj, &
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    
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
    type(t_spatialDiscretisation), intent(in) :: rdiscretisation
    
    ! The currently active element distribution in the discretisation.
    ! Allows the routine to get the currently active element type for
    ! trial and test functions.
    integer, intent(in) :: ielementDistribution
    
    ! The bilinear form which is currently to be evaluated:
    type(t_bilinearForm), intent(in) :: rform
    
    ! Start index of the current element block Ielements in the current element
    ! distribution ielementDistribution. If this is =1, the routine is called the
    ! first time for the current element distribution.
    integer, intent(in) :: ielementStartIdx

    ! Number of elements, where the coefficients must be computed.
    ! This is always a part of the element distribution.
    integer, intent(in) :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in) :: npointsPerElement
    
    ! A list of elements of length nelements where coefficients must
    ! be computed by this routine.
    integer, dimension(nelements), intent(in) :: Ielements
    
    ! A list of the corner vertices of all elements in Ielements
    real(DP), dimension(NDIM2D,TRIA_MAXNVE2D,nelements), intent(in) :: Dcoords
    
    ! A list of points in coordinates on the reference element.
    ! Each set of points corresponds to the corresponding element
    ! in Ielements
    real(DP), dimension(NDIM2D,npointsPerElement,nelements), intent(in) :: DcubPtsRef

    ! A list of points, corresponding to DcubPtsRef, in real coordinates.
    ! Each set of points corresponds to the corresponding element
    ! in Ielements.
    real(DP), dimension(NDIM2D,npointsPerElement,nelements), intent(in) :: DcubPtsReal
    
    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in trial space,nelements)
    integer, dimension(:,:) :: IdofsTrial
    
    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(#local DOF`s in test space,nelements)
    integer, dimension(:,:) :: IdofsTest
    
    ! The Jacobian matrix of the mapping between the reference and each
    ! real element, for all points on all elements.
    ! DIMENSION(#entries in the matrix,npointsPerEl,nelements)
    real(DP), dimension(:,:,:), intent(in) :: Djac
    
    ! The Jacobian determinant of the mapping of each point from the
    ! reference element to each real element
    real(DP), dimension(npointsPerElement,nelements), intent(in) :: Ddetj
  !</input>

  !<inputoutput>
    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional :: rcollection
  !</inputoutput>
  
  !<output>
    ! A list of all coefficients in front of all terms in the bilinear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the bilinear form.
    real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    end subroutine
    
  end interface
  
