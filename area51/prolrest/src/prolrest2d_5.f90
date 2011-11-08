!##############################################################################
!# ****************************************************************************
!# <name> prolrest2d_test5 </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module compares the local prolongation, restriction and
!# interpolation matrices for an element using L2-projection and using the
!# inter-level projection structure defined in multilevelprojection.f90 and
!# prints the matrix entries and the errors to the output stream.
!#
!# This module can be helpful if one has to implemented hard-coded projection
!# routines in the multilevelprojection module and wants to compare the
!# local matrices with the ones derived using L2-projection.
!# </purpose>
!##############################################################################

module prolrest2d_test5

  use fsystem
  use genoutput
  use storage
  use boundary
  use cubature
  use triangulation
  use element
  use spatialdiscretisation
  use genoutput
  use multileveloperators
  use multilevelprojection
  use stdoperators
    
  use prolrest_aux
  
  implicit none

contains

  ! ***************************************************************************

!<subroutine>

  subroutine prolrest2d_5
  
!<description>
!</description>

!</subroutine>

    ! Definitions of variables.
    !
    ! We need a couple of variables for this problem. Let's see...
    !
    ! An object for saving the domain:
    type(t_boundary) :: rboundary
    
    ! An object for saving the triangulation on the domain
    type(t_triangulation) :: rtriaC, rtriaF

    ! An object specifying the discretisation.
    ! This contains also information about trial/test functions,...
    type(t_blockDiscretisation) :: rdiscrC, rdiscrF
    
    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    type(t_matrixScalar) :: rmatrixC, rmatrixF, rmatrix2Lvl
    
    ! NLMAX receives the level where we want to solve.
    integer :: i,j
    
    ! Local dense matrices
    real(DP), dimension(:,:), allocatable :: Dmc,Dmf,Dn,DnT,Dimc,Dimf,&
        DAp,DAr,DAi,DBp,DBr,DBi,DEp,DEr,DEi
    integer, dimension(:), allocatable :: Ipivot
    integer :: nc, nf, info
    integer(I32) :: ielem
    logical :: bProl, bRest, bInterp, bFilterMat
    
    
    ! Specify the element to be used
    ielem = EL_Q2TB
    
    ! Compare Prolongation?
    bProl = .true.
    
    ! Compare Restriction?
    bRest = .true.
    
    ! Compare Interpolation?
    bInterp = .true.
    
    ! Filter local matrices?
    bFilterMat = .true.
    

    call output_lbrk()
    call output_separator(OU_SEP_STAR)
    call output_line('Test 5: Local Projection Matrix Comparison ')
    call output_separator(OU_SEP_STAR)
    
    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    call boundary_read_prm(rboundary, './pre/REFQUAD.prm')
        
    ! Now read in the basic coarse triangulation.
    call tria_readTriFile2D (rtriaC, './pre/QUAD.tri', rboundary)
     
    ! Create the coarse mesh
    call tria_initStandardMeshFromRaw (rtriaC,rboundary)
    
    ! Create the fine mesh
    call tria_refine2LevelOrdering(rtriaC,rtriaF,rboundary)
    call tria_initStandardMeshFromRaw (rtriaF,rboundary)

    ! Create 2 discretisations for both meshes
    call spdiscr_initBlockDiscr (rdiscrC,1,rtriaC, rboundary)
    call spdiscr_initBlockDiscr (rdiscrF,1,rtriaF, rboundary)

    call spdiscr_initDiscr_simple (rdiscrC%RspatialDiscr(1), &
                                   ielem,CUB_G5X5,rtriaC, rboundary)
    call spdiscr_initDiscr_simple (rdiscrF%RspatialDiscr(1), &
                                   ielem,CUB_G5X5,rtriaF, rboundary)
    
    ! And create the matrix structures
    call bilf_createMatrixStructure (rdiscrC%RspatialDiscr(1),&
                              LSYSSC_MATRIX9,rmatrixC)
    call bilf_createMatrixStructure (rdiscrF%RspatialDiscr(1),&
                              LSYSSC_MATRIX9,rmatrixF)
    
    ! Assemble the two mass matrices
    call stdop_assembleSimpleMatrix(rmatrixC,DER_FUNC,DER_FUNC)
    call stdop_assembleSimpleMatrix(rmatrixF,DER_FUNC,DER_FUNC)
    
    ! Create the 2-Level matrix structure
    call mlop_create2LvlMatrixStruct(rdiscrC%RspatialDiscr(1),&
            rdiscrF%RspatialDiscr(1),LSYSSC_MATRIX9,rmatrix2Lvl)
    
    ! Build the 2-Level-Mass matrix
    call mlop_build2LvlMassMatrix(rdiscrC%RspatialDiscr(1),&
                 rdiscrF%RspatialDiscr(1),.true.,rmatrix2Lvl)
    
    
    
    ! Get the number of DOFs for the coarse and fine mesh
    nc = rmatrixC%NEQ
    nf = rmatrixF%NEQ
    
    ! Allocate the local dense matrices
    allocate(Dmc(nc,nc))    ! coarse mass matrix
    allocate(Dimc(nc,nc))   ! coarse mass matrix inverse
    allocate(Dmf(nf,nf))    ! fine mass matrix
    allocate(Dimf(nf,nf))   ! fine mass matrix inverse
    allocate(Dn(nf,nc))     ! 2-level mass matrix
    allocate(DnT(nc,nf))    ! 2-level mass matrix transposed
    allocate(DAp(nf,nc))    ! prolongation matrix (L2-projection)
    allocate(DAr(nc,nf))    ! restriction matrix (L2-projection)
    allocate(DAi(nc,nf))    ! interpolation matrix (L2-projection)
    allocate(DBp(nf,nc))    ! prolongation matrix (hard-coded)
    allocate(DBr(nc,nf))    ! restriction matrix (hard-coded)
    allocate(DBi(nc,nf))    ! interpolation matrix (hard-coded)
    allocate(DEp(nf,nc))    ! prolongation erro matrix
    allocate(DEr(nc,nf))    ! restriction error matrix
    allocate(DEi(nc,nf))    ! interpolation error matrix
    
    ! Allocate pivot array
    allocate(Ipivot(max(nc,nf)))
    
    ! Densify the matrices
    call mat_densify(Dmc, rmatrixC)
    call mat_densify(Dmf, rmatrixF)
    call mat_densify(Dn, rmatrix2Lvl)
    
    ! Transpose 2-Level mass
    DnT = transpose(Dn)
    
    ! Invert the mass matrices
    call mat_identity(Dimc)
    call DGESV(nc, nc, Dmc, nc, Ipivot, Dimc, nc, info)
    call mat_identity(Dimf)
    call DGESV(nf, nf, Dmf, nf, Ipivot, Dimf, nf, info)
    
    if(bProl) then
    
      ! Calculate prolongation matrices
      DAp = matmul(Dimf, Dn)
      call prolrest_buildProlMatrix(rdiscrC, rdiscrF, DBp)
      
      ! Filter matrices if desired
      if(bFilterMat) then
        call mat_filterByEps(DAp)
        call mat_filterByEps(DBp)
      end if
      
      ! Calculate error matrix
      DEp = DAp - DBp
      call mat_filterByEps(DEp)
      
      ! Print Prolongation matrix
      call output_separator(OU_SEP_MINUS)
      call output_line('Prolongation Matrix Entries')
      call output_line('Entry      L2-projection           Hard-Coded              Error')
      do i = 1, nf
        do j = 1, nc
          call output_line('(' // trim(sys_si(i,2)) // ',' // trim(sys_si(j,2)) &
            // ') = ' // trim(sys_sdEP(DAp(i,j),20,13)) &
            // '    ' // trim(sys_sdEP(DBp(i,j),20,13)) &
            // '    ' // trim(sys_sdEP(DEp(i,j),20,13)))
        end do
      end do
      
    end if

    if(bRest) then
    
      ! Calculate restriction matrices
      DAr = matmul(DnT, Dimf)
      call prolrest_buildRestMatrix(rdiscrC, rdiscrF, DBr)
      
      ! Filter matrices if desired
      if(bFilterMat) then
        call mat_filterByEps(DAr)
        call mat_filterByEps(DBr)
      end if

      ! Calculate error matrix
      DEr = DAr - DBr
      call mat_filterByEps(DEr)
      
      ! Print Restriction matrix
      call output_separator(OU_SEP_MINUS)
      call output_line('Restriction Matrix Entries')
      call output_line('Entry      L2-projection           Hard-Coded              Error')
      do i = 1, nc
        do j = 1, nf
          call output_line('(' // trim(sys_si(i,2)) // ',' // trim(sys_si(j,2)) &
            // ') = ' // trim(sys_sdEP(DAr(i,j),20,13)) &
            // '    ' // trim(sys_sdEP(DBr(i,j),20,13)) &
            // '    ' // trim(sys_sdEP(DEr(i,j),20,13)))
        end do
      end do
    
    end if
    
    if(bInterp) then

      ! Calculate interpolation matrices
      DAi = matmul(Dimc, DnT)
      call prolrest_buildInterpMatrix(rdiscrC, rdiscrF, DBi)

      ! Filter matrices if desired
      if(bFilterMat) then
        call mat_filterByEps(DAi)
        call mat_filterByEps(DBi)
      end if

      ! Calculate error matrix
      DEi = DAi - DBi
      call mat_filterByEps(DEi)

      ! Print Interpolation matrix
      call output_separator(OU_SEP_MINUS)
      call output_line('Interpolation Matrix Entries')
      call output_line('Entry      L2-projection           Hard-Coded              Error')
      do i = 1, nc
        do j = 1, nf
          call output_line('(' // trim(sys_si(i,2)) // ',' // trim(sys_si(j,2)) &
            // ') = ' // trim(sys_sdEP(DAi(i,j),20,13)) &
            // '    ' // trim(sys_sdEP(DBi(i,j),20,13)) &
            // '    ' // trim(sys_sdEP(DEi(i,j),20,13)))
        end do
      end do
    
    end if
    
    ! Release pivot array
    deallocate(Ipivot)
    
    ! Release dense matrices
    deallocate(DEi)
    deallocate(DEr)
    deallocate(DEp)
    deallocate(DBi)
    deallocate(DBr)
    deallocate(DBp)
    deallocate(DAi)
    deallocate(DAr)
    deallocate(DAp)
    deallocate(DnT)
    deallocate(Dn)
    deallocate(Dimf)
    deallocate(Dmf)
    deallocate(Dimc)
    deallocate(Dmc)
    
    ! Release the matrices
    call lsyssc_releaseMatrix (rmatrix2Lvl)
    call lsyssc_releaseMatrix (rmatrixF)
    call lsyssc_releaseMatrix (rmatrixC)

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    call spdiscr_releaseBlockDiscr(rdiscrF)
    call spdiscr_releaseBlockDiscr(rdiscrC)
    
    ! Release the triangulation.
    call tria_done (rtriaF)
    call tria_done (rtriaC)
    
    ! Finally release the domain, that's it.
    call boundary_release (rboundary)
    
  end subroutine

end module
