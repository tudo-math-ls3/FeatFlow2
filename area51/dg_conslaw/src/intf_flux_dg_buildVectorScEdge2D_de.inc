  ! -*- mode: f90; -*-
  
  ! This file describes the interface for a coefficient function that is called
  ! during the matrix assembly of scalar vectors.
  
  interface
  
  !<subroutine>

    subroutine flux_dg_buildVectorScEdge2D_sim (&
              Dcoefficients,&
              IelementList,&
              Dnormal,&
              DedgeLength,&
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
    ! These are, for \int f * v ds (with rhs f and testfct v)
    ! \int dcoeff(1,...) f^i * v^i + \int dcoeff(2,...) f^a * v^a ds
    ! Where i=inside first element, a=outside first element.
    ! The normal is the one for the first element, -normal for second.
  !</description>
    
  !<input>
  real(DP), dimension(:,:), intent(in) :: Dnormal
  real(DP), dimension(:), intent(in) :: DedgeLength
  type(t_domainIntSubset), dimension(2), intent(in) :: rintSubset
  integer, dimension(:) , intent(in) :: IelementList
  !</input>
  
  !<inputoutput>
  type(t_collection), intent(inout), target, optional :: rcollection
  !</inputoutput>
  
  !<output>
  ! Dimension: 2 (testfcts on each side of the edge),nterms,ncubpoints,nedges
  real(DP), dimension(:,:,:,:), intent(out) :: Dcoefficients
  !</output>
    
  !</subroutine>
  
  end subroutine
    
  end interface
  
  
  
  
  
  
  
  
  
  
  
  
  
