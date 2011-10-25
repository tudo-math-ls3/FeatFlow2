  ! -*- mode: f90; -*-
  
  ! This file describes the interface for a coefficient function that is called
  ! during the matrix assembly of scalar vectors.
  
  interface
  
  !<subroutine>

    subroutine flux_dg_buildVectorBlEdge2D_sim (&
			  DfluxValues,&
			  rvectorSolBlock,&
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
    use linearsystemblock
    
  !<description>
    ! This subroutine is called during the vector assembly. It has to compute
    ! the coefficients in front of the terms of the linear form.
    ! These are the flux vectors (for each block) * normal vector.
  !</description>
    
  !<input>
  real(DP), dimension(:,:), intent(in) :: normal
  type(t_domainIntSubset), dimension(2), intent(in) :: rintSubset
  type(t_collection), intent(inout), target, optional :: rcollection
  type(t_vectorBlock), intent(in) :: rvectorSolBlock
  integer, dimension(:), intent(in) :: IelementList
    
  !</input>
  
  !<output>
  ! Dimension: nblocks,nterms,ncubpoints,nedges
  real(DP), dimension(:,:,:,:), intent(out) :: DfluxValues
  !</output>
    
  !</subroutine>
  
    end subroutine
    
  end interface
  
  
  
  
  
  
  
  
  
  
  
  
  
