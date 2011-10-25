  ! -*- mode: f90; -*-
  
  ! This file describes the interface for a coefficient function that is called
  ! during the matrix assembly of scalar vectors.
  
  interface
  
  !<subroutine>

    subroutine flux_dg_buildVectorScEdge2D_sim (&
              Dcoefficients,&
              IelementList,&
              normal,&
              rintSubSet,&
              rcollection )
    
    use fsystem
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    use linearsystemscalar
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    ! These are the flux*normal.
  !</description>
    
  !<input>
  real(DP), dimension(:,:), intent(in) :: normal
  type(t_domainIntSubset), dimension(2), intent(in) :: rintSubset
  integer, dimension(:) , intent(in) :: IelementList
  !</input>
  
  !<inputoutput>
  type(t_collection), intent(inout), target, optional :: rcollection
  !</inputoutput>
  
  !<output>
  ! Dimension: nterms,ncubpoints,nedges
  real(DP), dimension(:,:,:), intent(out) :: Dcoefficients
  !</output>
    
  !</subroutine>
  
  end subroutine
    
  end interface
  
  
  
  
  
  
  
  
  
  
  
  
  
