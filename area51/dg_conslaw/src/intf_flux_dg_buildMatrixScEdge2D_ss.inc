  ! -*- mode: f90; -*-
  
  ! This file describes the interface for a coefficient function that is called
  ! during the matrix assembly of scalar vectors.
  
  interface
  
  !<subroutine>

    subroutine flux_dg_buildMatrixScEdge2D_sim (&
			  DfluxValues,&
			  rvectorSol,&
			  IelementList,&
			  Dside,&
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
    ! This subroutine is called during the matrix assembly. It has to compute
    ! the values of flux*normal on the edges of the elements.
  !</description>
    
  !<input>
  real(DP), dimension(:,:), intent(in) :: normal
  type(t_domainIntSubset), dimension(2), intent(in) :: rintSubset
  type(t_vectorScalar), intent(in) :: rvectorSol
  integer, dimension(:), intent(in) :: IelementList
  !</input>
  
  !<inputoutput>
  type(t_collection), intent(inout), target, optional :: rcollection
  !</inputoutput>
  
  !<output>
  ! The flux vector * normal vector
  ! Dimension: nterms,ncubp,nedges
  ! DfluxValues(iterm, icubp, iedge)
  real(DP), dimension(:,:,:), intent(out) :: DfluxValues
  ! The coefficients for the DOFs from each side of the edge
  ! Dimension: 2 sides,nterms
  real(DP), dimension(:,:), intent(out) :: Dside
  !</output>
    
  !</subroutine>
  
    end subroutine
    
  end interface