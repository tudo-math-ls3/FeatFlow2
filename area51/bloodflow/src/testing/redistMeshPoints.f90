
  !*****************************************************************************

!<subroutine>

  subroutine bloodflow_redistMeshPoints(rbloodflow)

!<description>

    ! This subroutine redistributes the mesh points based on the
    ! equilibration of a two-dimensional mass spring system

!</description>

!<inputoutput>

    ! Bloodflow structure
    type(t_bloodflow), intent(inout) :: rbloodflow

!</inputoutput>
!</subroutine>

    ! local variables
    type(t_blockDiscretisation) :: rdiscretisation
    type(t_matrixBlock), dimension(1) :: Rmatrix
    type(t_matrixScalar) :: rmatrixScalar
    type(t_vectorBlock) :: rparticles, rforces, rincrement, rtemp, raux1, raux2
    type(t_linsolNode), pointer :: p_rsolverNode
    real(DP), dimension(:,:), pointer :: p_DobjectCoords, p_DvertexCoords
    real(DP), dimension(:), pointer :: p_Dmatrix, p_Dparticles, p_Dforces, p_Dindicator, p_Daux1, p_Daux2
    real(DP), dimension(1) :: dnorm1, dnorm2
    real(DP) :: dforcing, dared, dpred, daux, dscale, g0, g1, dg0
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:), pointer :: p_Kld, p_Kcol, p_Kdiagonal, p_Ksep
    integer :: h_Ksep, ierror, ite, ibacktrack


    ! Create matrix structure for standard P1 discretization
    call spdiscr_initBlockDiscr(rdiscretisation, 1,&
        rbloodflow%rtriangulation, rbloodflow%rboundary)
    call spdiscr_initDiscr_simple(rdiscretisation%RspatialDiscr(1), &
        EL_E001, SPDISC_CUB_AUTOMATIC, rbloodflow%rtriangulation, rbloodflow%rboundary)
    call bilf_createMatrixStructure (rdiscretisation%RspatialDiscr(1),&
        LSYSSC_MATRIX9, rmatrixScalar)
    
    ! Create scalar matrix with local 2x2 blocks
    rmatrixScalar%NVAR = 2
    rmatrixScalar%cmatrixFormat = LSYSSC_MATRIX9INTL
    rmatrixScalar%cinterleavematrixFormat = LSYSSC_MATRIX1
    call lsyssc_allocEmptyMatrix(rmatrixScalar, LSYSSC_SETM_UNDEFINED)

    ! Create 1-block system matrix and set points
    call lsysbl_createMatFromScalar(rmatrixScalar, Rmatrix(1), rdiscretisation)
    call lsyssc_getbase_double(rmatrixScalar, p_Dmatrix)
    call lsyssc_getbase_Kld(rmatrixScalar, p_Kld)
    call lsyssc_getbase_Kcol(rmatrixScalar, p_Kcol)
    call lsyssc_getbase_Kdiagonal(rmatrixScalar, p_Kdiagonal)
    
    ! Create diagonal separator
    h_Ksep = ST_NOHANDLE
    call storage_copy(rmatrixScalar%h_Kld, h_Ksep)
    call storage_getbase_int(h_Ksep, p_Ksep, rmatrixScalar%NEQ+1)
    
    ! Create vectors
    call lsysbl_createVecBlockIndMat(Rmatrix(1), rparticles, .true.)
    call lsysbl_createVecBlockIndMat(Rmatrix(1), rforces, .true.)
    call lsysbl_createVecBlockIndMat(Rmatrix(1), rincrement, .false.)
    call lsysbl_createVecBlockIndMat(Rmatrix(1), rtemp, .false.)
    call lsysbl_createVecBlockIndMat(Rmatrix(1), raux1, .false.)
    call lsysbl_createVecBlockIndMat(Rmatrix(1), raux2, .false.)
    call lsysbl_getbase_double(rparticles, p_Dparticles)
    call lsysbl_getbase_double(rforces, p_Dforces)
    call lsysbl_getbase_double(raux1, p_Daux1)
    call lsysbl_getbase_double(raux2, p_Daux2)
    
    ! Get indicator
    call lsyssc_getbase_double(rbloodflow%rindicator, p_Dindicator)

    ! Set pointers to coordinates vectors and triangulation data
    call storage_getbase_double2d(&
        rbloodflow%rtriangulation%h_DvertexCoords, p_DvertexCoords)
    call storage_getbase_int2d(&
        rbloodflow%rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    call storage_getbase_double2d(rbloodflow%h_DobjectCoords, p_DobjectCoords)

    ! Create a linear BiCGSTAB solver
    call linsol_initBiCGStab(p_rsolverNode)
    p_rsolverNode%ioutputLevel       = 0
    p_rsolverNode%istoppingCriterion = LINSOL_STOP_ONEOF
    p_rsolverNode%iresNorm           = LINALG_NORMEUCLID
    p_rsolverNode%depsRel            = 1e-2
    p_rsolverNode%depsAbs            = 1e-12
        
    ! Attach system matrix and initialize structures and data
    call linsol_setMatrices(p_RsolverNode, Rmatrix)
    call linsol_initStructure (p_rsolverNode, ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) stop
    call linsol_initData (p_rsolverNode, ierror)
    if (ierror .ne. LINSOL_ERR_NOERROR) stop

    ! Initialize particle positions by point coordinates
    call setParticles(p_DvertexCoords, rbloodflow%rtriangulation%NVT, p_Dparticles)

    ! Move point by hand
!    p_Dparticles(2*1013+1) = p_Dparticles(2*1013+1) + 0.15
!    p_Dparticles(2*1013+2) = p_Dparticles(2*1013+2) - 0.15

    ! Clear matrix and force vector
    call lalg_clearVector(p_Dmatrix)
    call lalg_clearVector(p_Dforces)
    
    ! Assemble matrix and force vector
    call assembleForces(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, p_IverticesAtElement,&
        rmatrixScalar%NA, rmatrixScalar%NEQ, rbloodflow%rtriangulation%NEL, &
        p_Dmatrix, p_Dparticles, p_Dforces, p_DvertexCoords, p_DobjectCoords)
    
    ! Set fixed points
    call setFixedPoints(p_Kld, p_Kcol, p_Kdiagonal, rmatrixScalar%NA,&
        rmatrixScalar%NEQ, p_Dmatrix, p_Dforces, p_DvertexCoords)
    
    ! Perform nonlinear iterations
    newton: do ite = 1,1000
      
      ! Solve linear problem
      call lsysbl_clearVector(rincrement)
      call linsol_solveAdaptively (p_rsolverNode, rincrement, rforces, rtemp)
      
      call lsysbl_copyVector(rparticles, rtemp)
      call lsysbl_copyVector(rforces, raux1)
      call lsysbl_blockMatVec(Rmatrix(1), rincrement, raux2, 1.0_DP, 0.0_DP)

      g0  =  lalg_scalarProductDble(p_Daux1, p_Daux1)
      dg0 = -lalg_scalarProductDble(p_Daux1, p_Daux2)
      
      dscale = 1.0_DP
      do ibacktrack = 1, 5
        
        ! Update particle positions
        call lsysbl_vectorLinearComb(rincrement, rtemp, dscale, 1.0_DP, rparticles)
        
        ! Compute relative changes
        if (ibacktrack .eq. 1) then
          call lsysbl_vectorNormBlock(rincrement, (/LINALG_NORMEUCLID/), Dnorm1)
          call lsysbl_vectorNormBlock(rparticles, (/LINALG_NORMEUCLID/), Dnorm2)
        end if

        ! Clear matrix and force vector
        call lalg_clearVector(p_Dmatrix)
        call lalg_clearVector(p_Dforces)
        call lalg_copyVector(p_Kld, p_Ksep)
        
        ! Assemble matrix and force vector
        call assembleForces(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, p_IverticesAtElement,&
            rmatrixScalar%NA, rmatrixScalar%NEQ, rbloodflow%rtriangulation%NEL, &
            p_Dmatrix, p_Dparticles, p_Dforces, p_DvertexCoords, p_DobjectCoords)
        
        ! Set fixed points
        call setFixedPoints(p_Kld, p_Kcol, p_Kdiagonal, rmatrixScalar%NA,&
            rmatrixScalar%NEQ, p_Dmatrix, p_Dforces, p_DvertexCoords)

        ! Compute actual and preducted reduction
        dared = p_rsolverNode%dinitialDefect - lalg_norm(p_Dforces, LINALG_NORMEUCLID)
        dpred = p_rsolverNode%dinitialDefect - sum((-p_Daux1 + dscale*p_Daux2)**2)
        daux  = (1-p_rsolverNode%depsRel) * p_rsolverNode%dinitialDefect

        ! Check sufficient decrease condition
        if (dpred .ge. daux .and. dared .ge. 0.01*daux) exit

        ! Perform line search
        g1   = lalg_scalarProductDble(p_Dforces, p_Dforces)
        daux = -dg0 / (g1-g0-2*dg0)

        ! Scale solution and reduction criterion
        call lsysbl_scaleVector(rincrement, daux)
      end do

      print *, "Relative changes", Dnorm1/Dnorm2, p_rsolverNode%dinitialDefect
      
      ! Check relative changes and exit Newton algorithm
      if (ite .gt. 1 .and. Dnorm1(1)/Dnorm2(1) .le. 1e-14) exit newton
      
      ! Update the new forcing term
      dforcing = abs(lalg_norm(p_Dforces, p_rsolverNode%iresNorm) -&
                     p_rsolverNode%dfinalDefect) / p_rsolverNode%dinitialDefect
      p_rsolverNode%depsRel = min(0.9_DP, max(dforcing, p_rsolverNode%depsRel**(0.5+sqrt(5.0)/2.0)))
                                       
    end do newton


    call lalg_copyVector(p_Kld, p_Ksep)
    
    ! Assemble matrix and force vector
    call assembleForces(p_Kld, p_Kcol, p_Kdiagonal, p_Ksep, p_IverticesAtElement,&
        rmatrixScalar%NA, rmatrixScalar%NEQ, rbloodflow%rtriangulation%NEL, &
        p_Dmatrix, p_Dparticles, p_Dforces, p_DvertexCoords, p_DobjectCoords, p_Dindicator)
    
    ! Set fixed points
    call setFixedPoints(p_Kld, p_Kcol, p_Kdiagonal, rmatrixScalar%NA,&
        rmatrixScalar%NEQ, p_Dmatrix, p_Dforces, p_DvertexCoords)
    

    print *, lalg_norm(p_Dforces, LINALG_NORML2)


    ! Update point coordinates
    call setVertexCoords(rbloodflow%rtriangulation%NVT, p_Dparticles, p_DvertexCoords)

    ! Release temporal memory
    call storage_free(h_Ksep)
    call lsysbl_releaseVector(rparticles)
    call lsysbl_releaseVector(rforces)
    call lsysbl_releaseVector(rincrement)
    call lsysbl_releaseVector(rtemp)
    call lsysbl_releaseVector(raux1)
    call lsysbl_releaseVector(raux2)
    call lsysbl_releaseMatrix(Rmatrix(1))
    call lsyssc_releaseMatrix(rmatrixScalar)
    call spdiscr_releaseBlockDiscr(rdiscretisation, .true.)
    call linsol_doneData(p_rsolverNode)
    call linsol_doneStructure(p_rsolverNode)
    call linsol_releaseSolver(p_rsolverNode)

  contains

    !***************************************************************************

    subroutine setParticles(DvertexCoords, n, Dparticles)

      ! This subroutine initializes the particle positions by the
      ! physical coordinates of the grid points

      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      integer, intent(in) :: n
      real(DP), dimension(NDIM2D, n), intent(out) :: Dparticles

      
      ! Copy data
      call lalg_copyVector(DvertexCoords, Dparticles)
    end subroutine setParticles

    !***************************************************************************

    subroutine setVertexCoords(n, Dparticles, DvertexCoords)

      ! This subroutine updates the physical coordinates of the grid
      ! points by the particle positions

      integer, intent(in) :: n
      real(DP), dimension(NDIM2D, n), intent(in) :: Dparticles
      real(DP), dimension(:,:), intent(inout) :: DvertexCoords
      
      ! Copy data
      call lalg_copyVector(Dparticles, DvertexCoords)
      
    end subroutine setVertexCoords
    
    !***************************************************************************
    
    subroutine assembleForces(Kld, Kcol, Kdiagonal, Ksep, IverticesAtElement,&
        na, neq, nel, Dmatrix, Dparticles, Dforces, DvertexCoords, DobjectCoords, Dindicator)

      ! This subroutine assembles the Jacobian matrix of the
      ! spring-mass system and the force vector which contains the
      ! forces for the linear edge springs satisfying Hooke's law.
      ! The force vector contains both the contributions from the edge
      ! springs and from the non-linear projection springs which are
      ! required to prevent collapsing of elements.
      
      real(DP), dimension(NDIM2D,neq), intent(in) :: Dparticles
      real(DP), dimension(:,:), intent(in) :: DvertexCoords, DobjectCoords
      real(DP), dimension(:), intent(in), optional :: Dindicator
      integer, dimension(:,:), intent(in) :: IverticesAtElement
      integer, dimension(:), intent(in) :: Kld, Kcol, Kdiagonal
      integer, intent(in) :: na, neq, nel
      
      real(DP), dimension(NDIM2D,NDIM2D,na), intent(inout) :: Dmatrix
      real(DP), dimension(NDIM2D,neq), intent(inout) :: Dforces
      integer, dimension(:), intent(inout) :: Ksep

      ! local variables
      real(DP), dimension(NDIM2D,NDIM2D) :: J_ij, J_il, J_ik, J_lj, J_lk, J_ls
      real(DP), dimension(NDIM2D) :: f_ij, f_il, d_ij, d_kj, d_il, p_l
      real(DP) :: dlength, dlength0, dprj, d2_ij, d3_ij, d2_il, d2_kj
      integer :: i,j,k,ii,ij,ji,jj,iel,ive


      !---------------------------------------------------------------------------
      ! (1) Assemble the stiffness spring forces and Jacobian matrix
      !     from their derivatives.
      !---------------------------------------------------------------------------
      
      ! Loop over all rows
      rows: do i = 1, neq
        
        ! Get position of diagonal entry
        ii = Kdiagonal(i)

        ! Loop over all off-diagonal matrix entries such that i<j
        cols: do ij = Kdiagonal(i)+1, Kld(i+1)-1

          ! Get node number j, the corresponding matrix positions ji
          j = Kcol(ij); jj = Kdiagonal(j); ji = Ksep(j); Ksep(j) = Ksep(j)+1
          
          ! Compute spring length and rest length
          dlength  = sqrt((Dparticles(1,i)-Dparticles(1,j))**2 +&
                          (Dparticles(2,i)-Dparticles(2,j))**2)
          dlength0 = sqrt((DvertexCoords(1,i)-DvertexCoords(1,j))**2 +&
                          (DvertexCoords(2,i)-DvertexCoords(2,j))**2)
          
          ! Compute coefficients and direction
          d_ij  = (Dparticles(:,j)-Dparticles(:,i))
          d2_ij = d_ij(1)*d_ij(1) + d_ij(2)*d_ij(2)
          d3_ij = d2_ij*sqrt(d2_ij)

          ! Compute local Jacobian matrix
          J_ij(1,1) = d2_ij - d_ij(1)*d_ij(1)
          J_ij(2,1) =       - d_ij(2)*d_ij(1)
          J_ij(1,2) =       - d_ij(1)*d_ij(2)
          J_ij(2,2) = d2_ij - d_ij(2)*d_ij(2)

          ! Scale local Jacobian matrix
          J_ij = J_ij / d3_ij

          ! Update local Jacobian matrix
          J_ij(1,1) = SPRING_STIFFNESS - SPRING_STIFFNESS * dlength0 * J_ij(1,1)
          J_ij(1,2) =                  - SPRING_STIFFNESS * dlength0 * J_ij(1,2)
          J_ij(2,1) =                  - SPRING_STIFFNESS * dlength0 * J_ij(2,1)
          J_ij(2,2) = SPRING_STIFFNESS - SPRING_STIFFNESS * dlength0 * J_ij(2,2)

          ! Apply off-diagonal matrix entries
          Dmatrix(:,:,ij) = J_ij
          Dmatrix(:,:,ji) = J_ij

          ! Update diagonal matrix entris
          Dmatrix(:,:,ii) = Dmatrix(:,:,ii) - J_ij
          Dmatrix(:,:,jj) = Dmatrix(:,:,jj) - J_ij
          
          ! Compute force vector by Hooke's law for edge (i,j)
          f_ij =  SPRING_STIFFNESS * (dlength-dlength0) * d_ij / dlength
    
          ! Apply force vector to nodes i and j
          Dforces(:,i) = Dforces(:,i) - f_ij
          Dforces(:,j) = Dforces(:,j) + f_ij
          
        end do cols
      end do rows

      
      !---------------------------------------------------------------------------
      ! (2) Assemble the repulsion spring forces and Jacobian
      !     matrix from their derivatives. The convention is as
      !     follows: For each element, we loop over all three vertices
      !     (i,j,k).  Vertex $i$ is the one which is projected onto
      !     the line through vertices $j$ and $k$. The projection
      !     point is referred to as $p_l$. The Jacobian matrix is
      !     computed making use of the chain rule. In short, a 2x2
      !     Jacobian $J_{il}$ is assembled which corresponds to the
      !     derivative of the the force vector with respect to the
      !     projection point $p_l$. A second 2x2 Jacobian $J_{ls}$
      !     stands for the derivative of the projection point $p_l$
      !     with respect to the position of the grid point $m_s$.
      !---------------------------------------------------------------------------

      ! Loop over all elements
      elems: do iel = 1, nel

        ! Loop over all corners
        do ive = 1, 3
          
          ! Get global vertex numbers
          i = p_IverticesAtElement(mod(ive-1, 3)+1, iel)
          j = p_IverticesAtElement(mod(ive,   3)+1, iel)
          k = p_IverticesAtElement(mod(ive+1, 3)+1, iel)
          
          ! Compute auxiliary coefficients
          d_ij  = DvertexCoords(:,j) - DvertexCoords(:,i)
          d2_ij = d_ij(1)*d_ij(1) + d_ij(2)*d_ij(2)

          d_kj  = DvertexCoords(:,j) - DvertexCoords(:,k)
          d2_kj = d_kj(1)*d_kj(1) + d_kj(2)*d_kj(2)

          ! Compute position of initial projection point
          p_l = DvertexCoords(:,j) - (d_ij(1)*d_kj(1) + d_ij(2)*d_kj(2)) * d_kj / d2_kj
     
          ! Compute square of rest length of the projection spring
          dlength0 = (DvertexCoords(1,i)-p_l(1))**2 +&
                     (DvertexCoords(2,i)-p_l(2))**2

          ! Compute auxiliary coefficients
          d_ij  = Dparticles(:,j) - Dparticles(:,i)
          d2_ij = d_ij(1)*d_ij(1) + d_ij(2)*d_ij(2)

          d_kj  = Dparticles(:,j) - Dparticles(:,k)
          d2_kj = d_kj(1)*d_kj(1) + d_kj(2)*d_kj(2)

          ! Compute position of corrent projection point
          p_l = Dparticles(:,j) - (d_ij(1)*d_kj(1) + d_ij(2)*d_kj(2)) * d_kj / d2_kj
          
          
          ! Compute nonlinear spring force
          d_il  = p_l - Dparticles(:,i)
          d2_il = d_il(1)*d_il(1) + d_il(2)*d_il(2)
          f_il  = SPRING_STIFFNESS * (d2_il - dlength0) * d_il / d2_il

          ! Apply force vector to node i
          Dforces(:,i) = Dforces(:,i) - f_il
          
          ! Compute local Jacobian matrix $dF_i/Dp_l$
          J_il(1,1) = d2_il - 2 * d_il(1)*d_il(1)
          J_il(2,1) =       - 2 * d_il(2)*d_il(1)
          J_il(1,2) =       - 2 * d_il(1)*d_il(2)
          J_il(2,2) = d2_il - 2 * d_il(2)*d_il(2)

          ! Scale local Jacobian matrix $dF_i/Dp_l$
          J_il = -SPRING_STIFFNESS * dlength0 * J_il / (d2_il*d2_il)

          ! Update local Jacobian matrix $dF_i/Dp_l$
          J_il(1,1) = SPRING_STIFFNESS + J_il(1,1)
          J_il(2,2) = SPRING_STIFFNESS + J_il(2,2)
          
          
          ! Compute local Jacobian matrix $dp_l/Dm_s$ for $s=i$
          J_ls(1,1) = d_kj(1)*d_kj(1) / d2_kj
          J_ls(2,1) = d_kj(2)*d_kj(1) / d2_kj
          J_ls(1,2) = d_kj(1)*d_kj(2) / d2_kj
          J_ls(2,2) = d_kj(2)*d_kj(2) / d2_kj

          ! Compute the matrix-matrix product $J_{il} \times K_{ls}$
          J_ij = matmul(J_il, J_ls)

          ! Apply Jacobian matrix to the diagonal
          ii = Kdiagonal(i)
          Dmatrix(:,:,ii) = Dmatrix(:,:,ii) + J_ij

          
          ! Compute local Jacobian matrix $dp_l/Dm_s$ for $s=j,k$
          J_ls(1,1) = d2_kj - 2 * d_kj(1)*d_kj(1)
          J_ls(2,1) =       - 2 * d_kj(2)*d_kj(1)
          J_ls(1,2) =       - 2 * d_kj(1)*d_kj(2)
          J_ls(2,2) = d2_kj - 2 * d_kj(2)*d_kj(2)
          
          ! Scale local Jacobian matrix $dp_l/Dm_s$ for $s=j,k$
          J_ls = (d_ij(1)*d_kj(1)+d_ij(2)*d_kj(2)) * J_ls / (d2_kj*d2_kj)

          ! Compute local Jacobian matrix $dp_l/dm_j$
          J_lj(1,1) = 1 - J_ls(1,1) - d_kj(1)*(d_ij(1)+d_kj(1)) / d2_kj
          J_lj(2,1) =   - J_ls(2,1) - d_kj(2)*(d_ij(1)+d_kj(1)) / d2_kj
          J_lj(1,2) =   - J_ls(1,2) - d_kj(1)*(d_ij(2)+d_kj(2)) / d2_kj
          J_lj(2,2) = 1 - J_ls(2,2) - d_kj(2)*(d_ij(2)+d_kj(2)) / d2_kj

          
          ! Compute local Jacobian matrix $dp_l/dm_k$
          J_lk(1,1) = J_ls(1,1) + d_kj(1)*d_ij(1) / d2_kj
          J_lk(2,1) = J_ls(2,1) + d_kj(2)*d_ij(1) / d2_kj
          J_lk(1,2) = J_ls(1,2) + d_kj(1)*d_ij(2) / d2_kj
          J_lk(2,2) = J_ls(2,2) + d_kj(2)*d_ij(2) / d2_kj


          ! Compute the matrix-matrix products $J_{il} \times K_{ls}$ for $s=j,k$
          J_ij = matmul(J_il, J_lj)
          J_ik = matmul(J_il, J_lk)


          ! Loop over i-th row and search for column $j$ and $k$
          do ij = Kld(i), Kld(i+1)-1
            if (Kcol(ij) .eq. j) Dmatrix(:,:,ij) = Dmatrix(:,:,ij) + J_ij
            if (Kcol(ij) .eq. k) Dmatrix(:,:,ij) = Dmatrix(:,:,ij) + J_ik
          end do

        end do
      end do elems
      
      !---------------------------------------------------------------------------
      ! (3) Assemble the attraction spring forces and the Jacobian
      !     matrix from their derivatives.
      !---------------------------------------------------------------------------

      i = 410

      do j = 83,88
        
        d_ij = DobjectCoords(:,j) - Dparticles(:,i)
        d2_ij = d_ij(1)*d_ij(1) + d_ij(2)*d_ij(2)
        d3_ij = sqrt(d2_ij)**3
        
        f_ij = SPRING_ATTRACTION * d_ij
        Dforces(:,i) = Dforces(:,i) - f_ij

        J_ij(1,1) = -SPRING_ATTRACTION
        J_ij(2,1) = 0
        J_ij(1,2) = 0
        J_ij(2,2) = -SPRING_ATTRACTION

        ii = Kdiagonal(i)
        Dmatrix(:,:,ii) = Dmatrix(:,:,ii) + J_ij
      end do
      
    end subroutine assembleForces

    !***************************************************************************
    
    subroutine setFixedPoints(Kld, Kcol, Kdiagonal, na, neq, Dmatrix, Dforces, DvertexCoords)

      ! This subroutine nullifies the force vector for fixed points
      ! and replaces the corresponding row of the matrix by that of
      ! the identity matrix.

      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      integer, dimension(:), intent(in) :: Kld, Kcol, Kdiagonal
      integer, intent(in) :: na, neq
      
      real(DP), dimension(NDIM2D,NDIM2D,na), intent(inout) :: Dmatrix
      real(DP), dimension(NDIM2D,neq), intent(inout) :: Dforces

      ! local variables
      integer :: i,ii,ij


      ! Loop over all rows
      rows: do i = 1, neq

        if (DvertexCoords(2,i) .ge. 1.5_DP .or. DvertexCoords(2,i) .le. -1.5_DP .or.&
            DvertexCoords(1,i) .ge. 2_DP-SYS_EPSREAL_DP .or. DvertexCoords(1,i) .le.  SYS_EPSREAL_DP) then

          ! Nullify forces
          Dforces(:,i) = 0.0_DP

          ! Nullify row in matrix
          do ij = Kld(i), Kld(i+1)-1
            Dmatrix(:,:,ij) = 0.0_DP
          end do

          ! Set identity matrix at diagonal block
          ii = Kdiagonal(i)
          Dmatrix(1,1,ii) = 1.0_DP
          Dmatrix(2,2,ii) = 1.0_DP

        end if
      end do rows
      
    end subroutine setFixedPoints

  end subroutine bloodflow_redistMeshPoints

