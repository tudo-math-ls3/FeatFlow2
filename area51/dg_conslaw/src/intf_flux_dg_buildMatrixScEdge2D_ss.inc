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
    ! Let us take a look at the bilinearform
    ! $\int_{all Edges} w (v*n) u ds$
    ! The testfunktion w on the two sides of each edge will be computed automatically
    ! in this version of the matrixbuildroutine - so we don't have to take care about that.
    ! But we have to take care about the flux (v*n) u.
    ! The integral is split into two terms
    ! $\int_{all Edges} w (v*n) Dside(1,1) u^i   +    w (v*n) Dside(2,1) u^a ds$
    ! where u^i and u^a are the inner and outer values of the solution u on the edge
    ! with respect to the (outward) normal vector n in normals
    !
    !      u^i | u^a
    !          |
    !          |---> normal
    !          |
    !          |
    !          Edge
    !
    ! So if we wanted a centered flux, we would choose
    ! Dside(1,1) = 0.5_dp
    ! Dside(2,1) = 0.5_dp
    ! and then looping over all edges and quadrature points set
    ! DfluxValues(1,ipoint,iedge) = beta1*normal(1,iedge) + beta2*normal(2,iedge)
    ! Then we would call it with
    ! rform%itermCount = 1
    ! rform%Idescriptors(1,1) = DER_FUNC
    ! rform%Idescriptors(2,1) = DER_FUNC
    ! rform%ballCoeffConstant = .false.
    ! rform%BconstantCoeff = .false.
    ! call bilf_dg_buildMatrixScEdge2D_ss (rform, CUB_G5_1D, .false., rmatrix,&
    !      rsol, raddTriaData, flux_dg_buildMatrixScEdge2D_sim,rcollection)
    !
    !
    ! For an upwind flux, we would split the sum in two terms, one taking
    ! only the inner solution value into account and one taking only the
    ! outer value into account
    ! Dside(1,1) = 1.0_dp
    ! Dside(2,1) = 0.0_dp
    ! Dside(1,2) = 0.0_dp
    ! Dside(2,2) = 1.0_dp
    ! Then wo would loop over all edges and integration points, for example:
    ! do iedge = 1, size(DfluxValues,3)
    !   do ipoint = 1, size(DfluxValues,2)
    !     ! Get x and y value of domain
    !     dx = rintSubset(1)%p_DcubPtsReal(1,ipoint,iedge)
    !     dy = rintSubset(1)%p_DcubPtsReal(2,ipoint,iedge)
    !     
    !     ! Get velocity 
    !     velx = dx
    !     vely = dy
    !
    !     ! Calculate v*n
    !     vn = velx*normal(1,iedge) + vely*normal(2,iedge)
    !      
    !     ! Check if velocity points in or out of first element
    !     if (DfluxValues(1,ipoint,iedge).ge.0.0_dp) then
    !        ! Points out so u^i is the upwind value
    !        DfluxValues(1,ipoint,iedge) = vn
    !        DfluxValues(2,ipoint,iedge) = 0.0_dp
    !     else
    !        ! Points in so u^a is the upwind value
    !        DfluxValues(1,ipoint,iedge) = 0.0_dp
    !        DfluxValues(2,ipoint,iedge) = vn
    !     end if
    !      
    !     ! Take care of boundary
    !     if (ielementList(iedge)==0) then
    !       ! We are at the boundary
    !       
    !       ! This value is not written anyways because u^a does not exist
    !       ! on the boundary
    !       DfluxValues(2,ipoint,iedge) = 0.0_dp
    !
    !       ! If we are at the inflow part of the boundary
    !       if (vn<0.0_dp) then
    !         ! Delete the term because it is usually put into the RHS
    !         DfluxValues(1,ipoint,iedge) = 0.0_dp
    !       end if
    !     end if
    !
    !   end do
    ! end do
    !
    ! Then call it with
    ! rform%itermCount = 2
    ! rform%Idescriptors(1,1) = DER_FUNC
    ! rform%Idescriptors(2,1) = DER_FUNC
    ! rform%Idescriptors(1,2) = DER_FUNC
    ! rform%Idescriptors(2,2) = DER_FUNC
    ! rform%ballCoeffConstant = .false.
    ! rform%BconstantCoeff = .false.
    ! call bilf_dg_buildMatrixScEdge2D_ss (rform, CUB_G5_1D, .false., rmatrix,&
    !        rsol, raddTriaData, flux_dg_buildMatrixScEdge2D_sim,rcollection)
    
    
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