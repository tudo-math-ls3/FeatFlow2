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

MODULE prolrest2d_test5

  USE fsystem
  USE genoutput
  USE storage
  USE boundary
  USE cubature
  USE triangulation
  USE spatialdiscretisation
  USE genoutput
  USE multileveloperators
  USE multilevelprojection
  USE stdoperators
    
  USE prolrest_aux
  
  IMPLICIT NONE

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE prolrest2d_5
  
!<description>
!</description>

!</subroutine>

    ! Definitions of variables.
    !
    ! We need a couple of variables for this problem. Let's see...
    !
    ! An object for saving the domain:
    TYPE(t_boundary) :: rboundary
    
    ! An object for saving the triangulation on the domain
    TYPE(t_triangulation) :: rtriaC, rtriaF

    ! An object specifying the discretisation.
    ! This contains also information about trial/test functions,...
    TYPE(t_blockDiscretisation) :: rdiscrC, rdiscrF
    
    ! A block matrix and a couple of block vectors. These will be filled
    ! with data for the linear solver.
    TYPE(t_matrixScalar) :: rmatrixC, rmatrixF, rmatrix2Lvl
    
    ! NLMAX receives the level where we want to solve.
    INTEGER :: i,j
    
    ! Local dense matrices
    real(DP), dimension(:,:), allocatable :: Dmc,Dmf,Dn,DnT,Dimc,Dimf,&
        DAp,DAr,DAi,DBp,DBr,DBi,DEp,DEr,DEi
    integer, dimension(:), allocatable :: Ipivot
    integer :: nc, nf, info
    integer(I32) :: ielem
    logical :: bProl, bRest, bInterp
    
    
    ! Specify the element to be used
    ielem = EL_E037
    
    ! Compare Prolongation?
    bProl = .TRUE.
    
    ! Compare Restriction?
    bRest = .TRUE.
    
    ! Compare Interpolation?
    bInterp = .TRUE.
    

    CALL output_lbrk()
    CALL output_separator(OU_SEP_STAR)
    CALL output_line('Test 5: Local Projection Matrix Comparison ')
    CALL output_separator(OU_SEP_STAR)
    
    ! At first, read in the parametrisation of the boundary and save
    ! it to rboundary.
    CALL boundary_read_prm(rboundary, './pre/REFQUAD.prm')
        
    ! Now read in the basic coarse triangulation.
    CALL tria_readTriFile2D (rtriaC, './pre/QUAD.tri', rboundary)
     
    ! Create the coarse mesh
    CALL tria_initStandardMeshFromRaw (rtriaC,rboundary)
    
    ! Create the fine mesh
    CALL tria_refine2LevelOrdering(rtriaC,rtriaF,rboundary)
    CALL tria_initStandardMeshFromRaw (rtriaF,rboundary)

    ! Create 2 discretisations for both meshes
    CALL spdiscr_initBlockDiscr2D (rdiscrC,1,rtriaC, rboundary)
    CALL spdiscr_initBlockDiscr2D (rdiscrF,1,rtriaF, rboundary)

    CALL spdiscr_initDiscr_simple (rdiscrC%RspatialDiscr(1), &
                                   ielem,CUB_G5X5,rtriaC, rboundary)
    CALL spdiscr_initDiscr_simple (rdiscrF%RspatialDiscr(1), &
                                   ielem,CUB_G5X5,rtriaF, rboundary)
    
    ! And create the matrix structures
    CALL bilf_createMatrixStructure (rdiscrC%RspatialDiscr(1),&
                              LSYSSC_MATRIX9,rmatrixC)
    CALL bilf_createMatrixStructure (rdiscrF%RspatialDiscr(1),&
                              LSYSSC_MATRIX9,rmatrixF)
    
    ! Assemble the two mass matrices
    CALL stdop_assembleSimpleMatrix(rmatrixC,DER_FUNC,DER_FUNC)
    CALL stdop_assembleSimpleMatrix(rmatrixF,DER_FUNC,DER_FUNC)
    
    ! Create the 2-Level matrix structure
    CALL mlop_create2LvlMatrixStruct(rdiscrC%RspatialDiscr(1),&
            rdiscrF%RspatialDiscr(1),LSYSSC_MATRIX9,rmatrix2Lvl)
    
    ! Build the 2-Level-Mass matrix
    CALL mlop_build2LvlMassMatrix(rdiscrC%RspatialDiscr(1),&
                 rdiscrF%RspatialDiscr(1),.TRUE.,rmatrix2Lvl)
    
    
    
    ! Get the number of DOFs for the coarse and fine mesh
    nc = rmatrixC%NEQ
    nf = rmatrixF%NEQ
    
    ! Allocate the local dense matrices
    ALLOCATE(Dmc(nc,nc))    ! coarse mass matrix
    ALLOCATE(Dimc(nc,nc))   ! coarse mass matrix inverse
    ALLOCATE(Dmf(nf,nf))    ! fine mass matrix
    ALLOCATE(Dimf(nf,nf))   ! fine mass matrix inverse
    ALLOCATE(Dn(nf,nc))     ! 2-level mass matrix
    ALLOCATE(DnT(nc,nf))    ! 2-level mass matrix transposed
    ALLOCATE(DAp(nf,nc))    ! prolongation matrix (L2-projection)
    ALLOCATE(DAr(nc,nf))    ! restriction matrix (L2-projection)
    ALLOCATE(DAi(nc,nf))    ! interpolation matrix (L2-projection)
    ALLOCATE(DBp(nf,nc))    ! prolongation matrix (hard-coded)
    ALLOCATE(DBr(nc,nf))    ! restriction matrix (hard-coded)
    ALLOCATE(DBi(nc,nf))    ! interpolation matrix (hard-coded)
    ALLOCATE(DEp(nf,nc))    ! prolongation erro matrix
    ALLOCATE(DEr(nc,nf))    ! restriction error matrix
    ALLOCATE(DEi(nc,nf))    ! interpolation error matrix
    
    ! Allocate pivot array
    ALLOCATE(Ipivot(MAX(nc,nf)))
    
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
      DAp = MATMUL(Dimf, Dn)
      call prolrest_buildProlMatrix(rdiscrC, rdiscrF, DBp)
      
      ! Calculate error matrix
      DEp = DAp - DBp
      call mat_filterByEps(DEp)
      
      ! Print Prolongation matrix
      CALL output_separator(OU_SEP_MINUS)
      CALL output_line('Prolongation Matrix Entries')
      CALL output_line('Entry      L2-projection           Hard-Coded              Error')
      DO i = 1, nf
        do j = 1, nc
          CALL output_line('(' // TRIM(sys_si(i,2)) // ',' // TRIM(sys_si(j,2)) &
            // ') = ' // TRIM(sys_sdEP(DAp(i,j),20,13)) &
            // '    ' // TRIM(sys_sdEP(DBp(i,j),20,13)) &
            // '    ' // TRIM(sys_sdEP(DEp(i,j),20,13)))
        end do
      END DO
      
    end if

    if(bRest) then
    
      ! Calculate restriction matrices
      DAr = MATMUL(DnT, Dimf)
      call prolrest_buildRestMatrix(rdiscrC, rdiscrF, DBr)
      
      ! Calculate error matrix
      DEr = DAr - DBr
      call mat_filterByEps(DEr)
      
      ! Print Restriction matrix
      CALL output_separator(OU_SEP_MINUS)
      CALL output_line('Restriction Matrix Entries')
      CALL output_line('Entry      L2-projection           Hard-Coded              Error')
      DO i = 1, nc
        do j = 1, nf
          CALL output_line('(' // TRIM(sys_si(i,2)) // ',' // TRIM(sys_si(j,2)) &
            // ') = ' // TRIM(sys_sdEP(DAr(i,j),20,13)) &
            // '    ' // TRIM(sys_sdEP(DBr(i,j),20,13)) &
            // '    ' // TRIM(sys_sdEP(DEr(i,j),20,13)))
        end do
      END DO
    
    end if
    
    if(bInterp) then

      ! Calculate interpolation matrices
      DAi = MATMUL(Dimc, DnT)
      call prolrest_buildInterpMatrix(rdiscrC, rdiscrF, DBi)

      ! Calculate error matrix
      DEi = DAi - DBi
      call mat_filterByEps(DEi)

      ! Print Interpolation matrix
      CALL output_separator(OU_SEP_MINUS)
      CALL output_line('Interpolation Matrix Entries')
      CALL output_line('Entry      L2-projection           Hard-Coded              Error')
      DO i = 1, nc
        do j = 1, nf
          CALL output_line('(' // TRIM(sys_si(i,2)) // ',' // TRIM(sys_si(j,2)) &
            // ') = ' // TRIM(sys_sdEP(DAi(i,j),20,13)) &
            // '    ' // TRIM(sys_sdEP(DBi(i,j),20,13)) &
            // '    ' // TRIM(sys_sdEP(DEi(i,j),20,13)))
        end do
      END DO
    
    end if
    
    ! Release pivot array
    DEALLOCATE(Ipivot)
    
    ! Release dense matrices
    DEALLOCATE(DEi)
    DEALLOCATE(DEr)
    DEALLOCATE(DEp)
    DEALLOCATE(DBi)
    DEALLOCATE(DBr)
    DEALLOCATE(DBp)
    DEALLOCATE(DAi)
    DEALLOCATE(DAr)
    DEALLOCATE(DAp)
    DEALLOCATE(DnT)
    DEALLOCATE(Dn)
    DEALLOCATE(Dimf)
    DEALLOCATE(Dmf)
    DEALLOCATE(Dimc)
    DEALLOCATE(Dmc)
    
    ! Release the matrices
    CALL lsyssc_releaseMatrix (rmatrix2Lvl)
    CALL lsyssc_releaseMatrix (rmatrixF)
    CALL lsyssc_releaseMatrix (rmatrixC)

    ! Release the discretisation structure and all spatial discretisation
    ! structures in it.
    CALL spdiscr_releaseBlockDiscr(rdiscrF)
    CALL spdiscr_releaseBlockDiscr(rdiscrC)
    
    ! Release the triangulation. 
    CALL tria_done (rtriaF)
    CALL tria_done (rtriaC)
    
    ! Finally release the domain, that's it.
    CALL boundary_release (rboundary)
    
  END SUBROUTINE

END MODULE
