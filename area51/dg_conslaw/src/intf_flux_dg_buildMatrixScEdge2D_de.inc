  ! -*- mode: f90; -*-
  
  ! This file describes the interface for a coefficient function that is called
  ! during the matrix assembly of scalar vectors.
  
  interface
  
  !<subroutine>

    subroutine flux_dg_buildMatrixScEdge2D_sim (&
			  DCoefficients,&
			  IelementList,&
        normal,&
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
    ! This subroutine is called during the matrix assembly. It has to compute
    ! the coeffitients for all combinations of trial and testfunctions in (i/1)
    ! and outside (a/2) of the first element (so in the second element).
    ! If we are at the boundary, only the i,i (1,1) combination is written,
    ! as the other combination of trial/test functions do not exist.
    ! For an edge on the boundary, it is Ielementlist(iedge) = 0.
    !
    ! DCoefficients
    ! (2 sides for the trial functions, 2 sides for the test functions,
    !  iterm, icubp, iedge)
  !</description>
    
  !<input>
  real(DP), dimension(:,:), intent(in) :: normal
  real(dp), dimension(:), intent(in) :: DedgeLength
  type(t_domainIntSubset), dimension(2), intent(in) :: rintSubset
  integer, dimension(:), intent(in) :: IelementList
  !</input>
  
  !<inputoutput>
  type(t_collection), intent(inout), target, optional :: rcollection
  !</inputoutput>
  
  !<output>
  ! The coeffitient
  ! Dimension: 2,2,nterms,ncubp,nedges
  ! DCoefficients
  ! (2 sides for the trial functions, 2 sides for the test functions,
  !  iterm, icubp, iedge)
  real(DP), dimension(:,:,:,:,:), intent(out) :: DCoefficients
  !</output>
    
  !</subroutine>
  
    end subroutine
    
  end interface