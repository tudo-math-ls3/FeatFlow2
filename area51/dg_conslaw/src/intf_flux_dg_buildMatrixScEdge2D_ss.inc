  ! -*- mode: f90; -*-
  
  ! This file describes the interface for a coefficient function that is called
  ! during the matrix assembly of scalar vectors.
  
  interface
  
  !<subroutine>

    subroutine flux_dg_buildMatrixScEdge2D_sim (&
!              Dcoefficients,&
!              DsolVals,&
			  DfluxValues,&
			  rvectorSol,&
			  IelementList,&
			  Dside,&
              normal,&
!              DpointsReal,&
              rintSubSet,&
              rcollection )
    
    use fsystem
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    
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
!  real(DP), dimension(:,:,:), intent(inout) :: DsolVals
  real(DP), dimension(:,:,:), intent(out) :: DfluxValues
  real(DP), dimension(:,:), intent(in) :: normal
!  real(DP), dimension(:,:,:), intent(in) :: DpointsReal
  type(t_domainIntSubset), dimension(2), intent(in) :: rintSubset
  type(t_collection), intent(inout), target, optional :: rcollection
  type(t_vectorScalar), intent(in) :: rvectorSol
  integer, dimension(:), intent(in) :: IelementList
  real(DP), dimension(:,:), intent(out) :: Dside
    
  !</input>
  
  !<output>
!  real(DP), dimension(:,:), intent(out) :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    end subroutine
    
  end interface
  
  
  
  
  
  
  
  
  
  
  
  
  
