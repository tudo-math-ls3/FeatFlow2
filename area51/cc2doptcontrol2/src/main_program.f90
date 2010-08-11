!##############################################################################
!# ****************************************************************************
!# <name> main_program </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module solves an optimal control problem for the stationary and
!# nonstationary Navier-Stokes optimal control problem 
!#
!#  $$ min J(y,u) = 1/2||y-z||_{L^2} + \gamma/2||y(T)-z(T)||_{L^2} + \alphga/2||u||^2 $$
!#
!#  $$- \nu Delta(y) + y*\Nabla(y) + \Nabla p = f $$
!#  $$ \Nabla \cdot y = 0$$
!#  $$- \nu Delta(\lambda) - y*\Nabla(\lambda) + \lambda\Nabla y + \Nabla \xi = y-z $$
!#  $$ \Nabla \cdot \lambda = 0$$
!#              
!#
!# on a 2D domain for a 2D function $y=(y_1,y_2)$, a pressure $p$,
!# a dual velocity $\lambda$ and a dual pressure $\xi$. $u$ is the control
!# and $z$ a desired target function field.
!#
!# The routine splits up the tasks of reading the domain, creating 
!# triangulations, discretisation, solving, postprocessing and creanup into
!# different subroutines. The communication between these subroutines
!# is done using an application-specific structure saving problem data
!# as well as a collection structure for the communication with callback
!# routines.
!#
!# For the nonlinearity, the nonlinear solver is invoked. The
!# defect that is setted up there is preconditioned by a linear Multigrid
!# solver with a simple-VANKA smoother/preconditioner for
!# 2D saddle point problems, Jacobi-Type. As coarse grid solver,
!# UMFPACK is used.
!# </purpose>
!##############################################################################

module main_program

  use fsystem
  use genoutput
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use paramlist
  use fparser
  use statistics
  
  use collection
  use convection
    
!  use basicstructures
!  use paramtriainit
!  use spatialbc
!  use spacediscretisation

  use postprocessing
  
  use externalstorage
  use paramlist
  
  use spacetimevectors
  use analyticsolution
  
  use structuresoptflow
  use initsolver
  use structuresoptcspacetimenlsol
  use nonlinearoneshotspacetimesolver
  
  use structuresmain
  
  implicit none

!<globals>

  ! Directory containing the data files.
  character(LEN=SYS_STRLEN), save :: DIR_DATA = "./data";

!</globals>
  
contains

!  ! ***************************************************************************
!
!!<subroutine>
!  
!  subroutine vanka (rvanka, rvector, rrhs, domega, IelementList,celementP)
!  
!!<description>
!!</description>
!
!!<input>
!  ! t_vankaPointer2DNavSt structure that saves algorithm-specific parameters.
!  type(t_vankaPointer2DNavStOptC), intent(in) :: rvanka
!
!  ! The right-hand-side vector of the system
!  type(t_vectorBlock), intent(in)         :: rrhs
!  
!  ! Relaxation parameter. Standard=1.0_DP.
!  real(DP), intent(in)                    :: domega
!
!  ! A list of element numbers where VANKA should be applied to.
!  integer, dimension(:)     :: IelementList
!  
!  ! Element type for the pressure spaces
!  integer(I32), intent(in) :: celementP
!!</input>
!
!!<inputoutput>
!  ! The initial solution vector. Is replaced by a new iterate.
!  type(t_vectorBlock), intent(in)         :: rvector
!!</inputoutput>
!
!!</subroutine>
!
!    ! local variables
!    integer :: ielidx
!    integer, dimension(1) :: IelIdx2
!    integer, dimension(:,:), allocatable :: IdofP
!    
!    ! B- and D-matrices
!    integer, dimension(:), pointer :: p_Db1, p_Db2, p_Dd1, p_Dd2
!    integer, dimension(:), pointer :: p_KcolB, p_KcolD
!    
!    ! Allocate an array for the pressure DOF's.
!    allocate(IdofP(elem_igetNDofLoc(celementP),1))
!    
!    ! Loop through all elements in the list
!    do ielidx = 1,size(IelementList)
!    
!      ! On the element, get the local DOF's in the pressure space
!      IelIdx2(1) = ielidx
!      call dof_locGlobMapping_mult(rdiscretisation, IelIdx2, IdofP)
!    
!      ! The D-matrix allows us to access the velocity DOF's on the element,
!      ! as all P-DOFs are connected exactly to all velocity DOF's.
!    
!    
!    end do ! ielidx
!
!
!
!
!
!
!
!
!
!
!
!
!
!    ! local vairables
!    integer :: iel,ielidx,nmaxelements,nelementsinset,ilocalel
!    integer :: inode,idof
!    
!    integer, dimension(:), pointer :: p_KcolA11
!    integer, dimension(:), pointer :: p_KldA11
!    real(DP), dimension(:), pointer :: p_DA11,p_DA12,p_DA21,p_DA22
!    integer, dimension(:), pointer :: p_KcolA45
!    integer, dimension(:), pointer :: p_KldA45
!    real(DP), dimension(:), pointer :: p_DA44,p_DA45,p_DA54,p_DA55
!    real(DP), dimension(:), pointer :: p_DR41,p_DR52,p_DR51,p_DR42
!    integer, dimension(:), pointer :: p_KcolA12
!    integer, dimension(:), pointer :: p_KldA12
!    integer, dimension(:), pointer :: p_KcolM
!    integer, dimension(:), pointer :: p_KldM
!    real(DP), dimension(:), pointer             :: p_DM
!    integer, dimension(:), pointer :: p_KcolB
!    integer, dimension(:), pointer :: p_KldB
!    real(DP), dimension(:), pointer :: p_DB1
!    real(DP), dimension(:), pointer :: p_DB2
!    real(DP), dimension(:), pointer :: p_DD1
!    real(DP), dimension(:), pointer :: p_DD2
!    real(DP), dimension(:), pointer :: p_Da33,p_Da66
!    integer, dimension(:), pointer :: p_KdiagonalA33,p_KdiagonalA66
!    
!    ! Triangulation information
!    integer :: NEL
!    integer, dimension(:,:), pointer :: p_IedgesAtElement
!    real(DP), dimension(:), pointer :: p_Drhs,p_Dvector
!    
!    ! Local arrays for informations about one element
!    integer, parameter :: nnvel = 9      ! Q2 = 9 DOF's per velocity
!    integer, parameter :: nnpressure = 3 ! QP1 = 3 DOF's per pressure
!    integer, parameter :: nndualvel = 9      ! Q2 = 9 DOF's per dual velocity
!    integer, parameter :: nndualpressure = 3 ! QP1 = 3 DOF's per dual pressure
!    integer, parameter :: nnprimal = 2*nnvel+nnpressure ! Q2/Q2/QP10 = 4+4+1 = 9 DOF's per element
!    integer, parameter :: nnld = 2*nnprimal
!    integer, dimension(:,:),allocatable :: IdofGlobal
!    real(DP), dimension(:,:,:),allocatable :: AA
!    real(DP), dimension(nnld) :: FF
!    
!    ! Offsets of the 'local' solution parts in the 'local' solution vector
!    integer, parameter :: lofsu = 0
!    integer, parameter :: lofsv = nnvel
!    integer, parameter :: lofsp = 2*nnvel
!    integer, parameter :: lofsl1 = 2*nnvel+1
!    integer, parameter :: lofsl2 = 2*nnvel+1+nnvel
!    integer, parameter :: lofsxi = 2*nnvel+1+2*nnvel
!    
!    ! LAPACK temporary space
!    integer :: Ipiv(nnld),ilapackInfo
!    
!    ! Offset information in arrays.
!    ! Primal variables
!    integer     :: ioffsetu,ioffsetv,ioffsetp,j
!    
!    ! Dual variables
!    integer     :: ioffsetl1,ioffsetl2,ioffsetxi
!    
!    integer :: ia1,ia2,ib1,ib2,ia,ib,k
!    real(DP) :: daux,daux2
!    
!    ! Get pointers to the system matrix, so we don't have to write
!    ! so much - and it's probably faster.
!    
!    ! Structure of A11 is assumed to be the same as A22
!    p_KcolA11 => rvanka%p_KcolA11
!    p_KldA11 => rvanka%p_KldA11
!    p_DA11 => rvanka%p_DA11
!    p_DA22 => rvanka%p_DA22
!
!    ! Structure of A12 is assumed to be the same as A21.
!    ! Get A12 and A21 -- except for if the multipliers are =0, then
!    ! we switch them off by nullifying the pointers.
!    if (rvanka%Dmultipliers(1,2) .ne. 0.0_DP) then
!      p_KcolA12 => rvanka%p_KcolA12
!      p_KldA12 => rvanka%p_KldA12
!      p_DA12 => rvanka%p_DA12
!      p_DA21 => rvanka%p_DA21
!    else
!      nullify(p_KcolA12)
!      nullify(p_KldA12) 
!      nullify(p_DA12 )
!      nullify(p_DA21 )
!    end if
!    
!    p_KcolB => rvanka%p_KcolB
!    p_KldB => rvanka%p_KldB
!    p_DB1 => rvanka%p_DB1
!    p_DB2 => rvanka%p_DB2
!    p_DD1 => rvanka%p_DD1
!    p_DD2 => rvanka%p_DD2
!    
!    ! Structure of A44 is assumed to be the same as A55, A11 and A22
!    p_DA44 => rvanka%p_DA44
!    p_DA55 => rvanka%p_DA55
!
!    ! Structure of A45 is assumed to be the same as A54
!    p_KcolA45 => rvanka%p_KcolA45
!    p_KldA45 => rvanka%p_KldA45
!    p_DA45 => rvanka%p_DA45
!    p_DA54 => rvanka%p_DA54
!    
!    ! Mass matrix - if it's given, otherwise the pointers will be set to NULL
!    ! because of the initialisation of the structure!
!    p_KcolM => rvanka%p_KcolM
!    p_KldM => rvanka%p_KldM
!    p_DM => rvanka%p_DM14
!    
!    ! Coupling matrix in the dual equation at position (4:5,1:2). For a standard
!    ! system, there is A(4,1) = A(5,2) = M and A(5,1) = A(4,2) = 0.
!    ! For a Newton system, this block is completely decoupled!
!    p_DR41 => rvanka%p_DR41
!    p_DR42 => rvanka%p_DR42
!    p_DR51 => rvanka%p_DR51
!    p_DR52 => rvanka%p_DR52
!    
!    ! Diagonal submatrices A33 and A66 (if they exist)
!    if (rvanka%Dmultipliers(3,3) .ne. 0.0_DP) then
!      p_Da33 => rvanka%p_DA33
!      p_KdiagonalA33 => rvanka%p_KdiagonalA33
!    else
!      nullify(p_Da33)
!      nullify(p_KdiagonalA33)
!    end if
!    
!    if (rvanka%Dmultipliers(6,6) .ne. 0.0_DP) then
!      p_Da66 => rvanka%p_DA66
!      p_KdiagonalA66 => rvanka%p_KdiagonalA66
!    else
!      nullify(p_Da66)
!      nullify(p_KdiagonalA66)
!    end if
!
!    ! Get pointers to the vectors, RHS, get triangulation information
!    NEL = rvector%RvectorBlock(1)%p_rspatialDiscr%p_rtriangulation%NEL
!    call storage_getbase_int2d (rvector%RvectorBlock(1)%p_rspatialDiscr% &
!                                p_rtriangulation%h_IedgesAtElement, p_IedgesAtElement)
!    call lsysbl_getbase_double (rvector,p_Dvector)
!    call lsysbl_getbase_double (rrhs,p_Drhs)
!    
!    ! Get the relative offsets of the 2nd and 3rd solution of the component
!    ioffsetu = 0
!    ioffsetv = rvector%RvectorBlock(1)%NEQ
!    ioffsetp = ioffsetv+rvector%RvectorBlock(2)%NEQ
!
!    ! Get the offsets of lambda1, lambda2 and xi, so the offsets
!    ! of the dual solution vectors.
!    ioffsetl1 = ioffsetp+rvector%RvectorBlock(3)%NEQ
!    ioffsetl2 = ioffsetl1+rvector%RvectorBlock(4)%NEQ
!    ioffsetxi = ioffsetl2+rvector%RvectorBlock(5)%NEQ
!    
!    !=======================================================================
!    !     Block Gauss-Seidel on Schur Complement
!    !=======================================================================
!
!    ! Basic algorithm:
!    !
!    ! What are we doing here? Well, we want to perform 
!    ! *preconditioning*, i.e. we have to solve the problem
!    !
!    !   x_new  =  C^-1 (x_old)  =  C^-1 (F)  =  C^-1 (f,g)
!    !
!    ! for a "special" preconditioner C which we define in a moment.
!    ! This is equivalent to solving the system
!    !
!    !   C (x_new)  = x_old
!    !
!    ! C should be some approximation to A. Imagine our global system:
!    !
!    !     [ A   B ] (u) = (f)
!    !     [ B^t 0 ] (p)   (g)
!    !
!    ! In the Navier-Stokes equations with (u,p) being the preconditioned
!    ! vector, there should be g=0 - but this cannot be assumed
!    ! as it does not happen in general.
!    ! Now the algorithm for generating a new (u,p) vector from the old
!    ! one reads roughly as follows:
!    !
!    ! a) Restrict to a small part of the domain, in our case to one cell.
!    ! b) Fetch all the data (velocity, pressure) on that cell. On the
!    !    first cell, we have only "old" velocity entries. These values
!    !    are updated and then the calculation proceeds with the 2nd cell.
!    !
!    !           old                      new     
!    !        +---X---+                +---X---+
!    !        |       |                |       |
!    !    old X       X       -->  new X   X   X new
!    !        |   1   |                |   1   |
!    !        +---X---+                +---X---+
!    !           old                      new     
!    !
!    !    From the second cell on, there might be "old" data and "new" 
!    !    data on that cell - the old data that has not been updated and
!    !    perhaps some already updated velocity data from a neighbor cell.
!    !    
!    !           new     old                   new       new    
!    !        +---X---+---X---+             +---X---+---X---+
!    !        |     1 |     2 |             |     1 |     2 |
!    !    new X       X       X old --> new X       X       X new
!    !        |       |new    |             |       |newer  |
!    !        +---X---+---X---+             +---X---+---X---+
!    !           new     old                   new       new    
!    !
!    !    These values are updated and then the calculation proceeds
!    !    with the next cell.
!    !    As can be seen in the above picture, the "new" node in the
!    !    middle is even going to be a "newer" node when handled again
!    !    for the 2nd cell. This is meant by "Gauss-Seldel" character:
!    !    Information is updated subsequently by using "old" data and
!    !    "new" data from a previous calculation.
!    !
!    ! So we start with a loop over all elements in the list
!
!    nmaxelements = min(1000,size(IelementList))
!    allocate(AA(nnld,nnld,nmaxelements)
!    allocate(IdofGlobal(nnvel+nnpressure,nmaxelements))
!
!    ! Loop through the elements in sets.
!    do ielidx=1,size(IelementList),nmaxelements
!    
!      ! Number of elements in this set.
!      nelements = min(nmaxelements,size(IelementList)-ielidx+1)
!    
!      ! Clear the 'local system matrices'.
!      AA(:,:,:) = 0.0_DP
!      
!      ! Calculate the global DOF's on all the elements.
!      ! Corners, edge midpoints, element midpoints etc...
!      ! We assume that all velocity spaces are the same, that the
!      ! discretisation is uniform and that the pressure spaces are
!      ! the same.
!      ! DOF 1..9: Q2, DOF10..12: QP1.
!      do inode = 1,nelements
!        iel = IelementList(ielidx+inode-1)
!        IdofGlobal(1:4,inode) = p_IverticesAtElement(1:4,iel)
!        IdofGlobal(5:8,inode) = p_IedgesAtElement(1:4,iel)+NVT
!        IdofGlobal(9,inode)   = iel
!        IdofGlobal(10,inode)  = iel
!        IdofGlobal(11,inode)  = iel+NEL
!        IdofGlobal(12,inode)  = iel+2*NEL
!      end do
!      
!      ! Fetch the local A11 submatrices.
!      do ilocalel = 1,nelements
!        iel = IelementList(ielidx+ilocalel-1)
!        
!        ! Get the velocity-DOF we have to tackle:
!        idof = IdofGlobal(inode)
!        
!        
!      end do
!      
!      
!      ! Get the element number which is to be processed.
!      iel = IelementList(ielidx)
!    
!      ! We now have the element
!      !                                               
!      ! +---------+                       +----3----+
!      ! |         |                       |         |
!      ! |   IEL   |   with DOF's          4    P    2      
!      ! |         |                       |    Q0   |
!      ! +---------+                       +----1----+
!      !                                               
!      !
!      ! Fetch the pressure P on the current element into FF.
!      ! The numbers of the DOF's coincide with the definition
!      ! in dofmapping.f90!
!    
!      ! Get the primal pressure
!      FF(1+lofsp) = p_Drhs(iel+ioffsetp)
!      
!      ! Get the dual pressure
!      FF(1+lofsxi) = p_Drhs(iel+ioffsetxi)
!      
!      ! Get the velocity DOF's on the current element.
!      ! We assume: DOF 1..4 = edge.
!      ! That's the same implementation as in dofmapping.f90!
!      IdofGlobal(1:4) = p_IedgesAtElement(1:4,iel)
!
!      ! Loop over all U-nodes of that element.
!      do inode=1,nnvel
!      
!        ! Get the DOF we have to tackle:
!        idof = IdofGlobal(inode)
!        
!        ! Set FF initially to the value of the right hand
!        ! side vector that belongs to our current DOF corresponding
!        ! to inode.
!        
!        ! Primal equation
!        FF(inode+lofsu) = p_Drhs(idof+ioffsetu)
!        FF(inode+lofsv) = p_Drhs(idof+ioffsetv)
!
!        ! dual equation
!        FF(inode+lofsl1) = p_Drhs(idof+ioffsetl1)
!        FF(inode+lofsl2) = p_Drhs(idof+ioffsetl2)
!        
!        ! What do we have at this point?                           
!        ! FF     : "local" RHS vector belonging to the DOF's on the
!        !          current element                                 
!        ! AA     : Diagonal entries of A belonging to these DOF's  
!        !                                                          
!        ! And at the moment:                                       
!        ! idof      : number of current DOF on element IEL            
!        ! inode     : "local" number of DOF on element IEL, i.e.      
!        !              number of the edge         
!        !                     
!        ! Now comes the crucial point with the "update": How to         
!        ! subsequently update the vertex values, such that the whole    
!        ! thing still converges to the solution, even if a node         
!        ! is updated more than once? Here, we use a typical             
!        ! matrix-decomposition approach:                                
!        !                                                               
!        ! Again consider the problem:                                   
!        !                                                               
!        !    [ A   B  M     ] (u)  = (f  )                                        
!        !    [ B^t 0        ] (p)    (g  )                                        
!        !    [ M      A   B ] (l)    (fl )                                        
!        !    [        B^t 0 ] (xi)   (fxi)                                        
!        !                                                               
!        ! We assume, that all components in the vector (u,p) are        
!        ! given - except for the velocity and pressure unknowns 
!        ! on the current element; these 21 unknowns  
!        ! are located anywhere in the (u,p) vector. The idea is to      
!        ! shift "everything known" to the right hand side to obtain     
!        ! a system for only these unknowns!   
!        !
!        ! Extracting all the lines of the system that correspond to     
!        ! DOF's on our single element IEL results in rectangular
!        ! systems of the form                                            
!        !                                                               
!        !    [ === A^ === B~  === M^ ====   ] (| ) = (f1 )                                
!        !    [ B~^t       I1~               ] (u )   (f2 )                                
!        !    [ === M^ ===     === A^ === B~ ] (| )   (g  )                                   
!        !    [                B~^t       I2~] (p )   (fl1)
!        !                                     (| )   (fl2)
!        !                                     (l )   (flg)
!        !                                     (| )
!        !                                     (xi)
!        !                                     
!        !                                                               
!        ! B~ is a 8 x 2 matrix: As every velocity couples with at most  
!        ! 2*1 pressure elements on the adjacent cells, so we have       
!        ! 2 columns in the B-matrix.                                    
!        !                                                               
!        !        IEL                              IEL                   
!        !     |--------|             |--------|--------|                
!        !     |        |             |        |        |                
!        !     |   P    |      or     |   Q    X   P    |                
!        !     |   X    |             |        |        |                
!        !   --|--------|--           |--------|--------|                
!        !
!        !
!        ! Now, throw all summands to the RHS vector to build a local
!        ! 'defect' on our single element IEL.
!        !                                                               
!        !  (d1 ) = (f1 ) -  [ === A^ === B~  === M^ ====   ] (| )                                 
!        !  (d2 )   (f2 )    [ B~^t       I1~               ] (u )                                 
!        !  (dg )   (g  )    [ === M^ ===     === A^ === B~ ] (| )                                    
!        !  (dl1)   (fl1)    [                B~^t       I2~] (p ) 
!        !  (dl2)   (fl2)                                     (| ) 
!        !  (dlg)   (flg)                                     (l ) 
!        !                                                    (| )
!        !                                                    (xi)
!        !
!        ! Extract those entries in the A-, B- and M-matrices to our local
!        ! matrix AA, which belong to DOF's in our current solution vector.
!        !
!        ! At first build: fi = fi-Aui
!        
!        ia1 = p_KldA11(idof)
!        ia2 = p_KldA11(idof+1)-1
!        do ia = ia1,ia2
!          ! Calculate:
!          !
!          !   ( du  ) = ( du  ) - ( A11  .   .   .   .   .  ) ( u  )
!          !   ( dv  )   ( dv  )   (  .  A22  .   .   .   .  ) ( v  )
!          !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
!          !   ( dl1 )   ( dl1 )   (  .   .   .  A44  .   .  ) ( l1 )
!          !   ( dl2 )   ( dl2 )   (  .   .   .   .  A55  .  ) ( l2 )
!          !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )
!
!          J = p_KcolA11(ia)
!          
!          ! Primal equation:
!          FF(inode+lofsu) = FF(inode+lofsu) &
!                          - rvanka%Dmultipliers(1,1)*p_DA11(ia)*p_Dvector(J+ioffsetu)
!          FF(inode+lofsv) = FF(inode+lofsv) &
!                          - rvanka%Dmultipliers(2,2)*p_DA22(ia)*p_Dvector(J+ioffsetv)
!
!          ! dual equation
!          FF(inode+lofsl1) = FF(inode+lofsl1) &
!                           - rvanka%Dmultipliers(4,4)*p_DA44(ia)*p_Dvector(J+ioffsetl1)
!          FF(inode+lofsl2) = FF(inode+lofsl2) &
!                           - rvanka%Dmultipliers(5,5)*p_DA55(ia)*p_Dvector(J+ioffsetl2)
!          
!          ! Whereever we find a DOF that couples to another DOF on the 
!          ! same element, we put that to both A-blocks of our local matrix.
!          do k=1,nnvel
!            if (j .eq. IdofGlobal(k)) then
!              AA (inode+lofsu,k+lofsu) = p_DA11(ia)*rvanka%Dmultipliers(1,1)
!              AA (inode+lofsv,k+lofsv) = p_DA22(ia)*rvanka%Dmultipliers(2,2)
!              
!              AA (inode+lofsl1,k+lofsl1) = p_DA44(ia)*rvanka%Dmultipliers(4,4)
!              AA (inode+lofsl2,k+lofsl2) = p_DA55(ia)*rvanka%Dmultipliers(5,5)
!              exit
!            end if
!          end do          
!        end do
!
!        ! Handle the 'off-diagonal' matrices A12 and A21
!        
!        if (associated(p_KldA12)) then
!          ia1 = p_KldA12(idof)
!          ia2 = p_KldA12(idof+1)-1
!          do ia = ia1,ia2
!            ! Calculate:
!            !
!            !   ( du  ) = ( du  ) - (  .  A12  .   .   .   .  ) ( u  )
!            !   ( dv  )   ( dv  )   ( A21  .   .   .   .   .  ) ( v  )
!            !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
!            !   ( dl1 )   ( dl1 )   (  .   .   .   .   .   .  ) ( l1 )
!            !   ( dl2 )   ( dl2 )   (  .   .   .   .   .   .  ) ( l2 )
!            !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )
!
!            J = p_KcolA12(ia)
!            FF(inode+lofsu) = FF(inode+lofsu) &
!                            - rvanka%Dmultipliers(1,2)*p_DA12(ia)*p_Dvector(J+ioffsetv)
!            FF(inode+lofsv) = FF(inode+lofsv) &
!                            - rvanka%Dmultipliers(2,1)*p_DA21(ia)*p_Dvector(J+ioffsetu)
!            
!            ! Whereever we find a DOF that couples to another DOF on the 
!            ! same element, we put that to both A-blocks of our local matrix.
!            do k=1,nnvel
!              if (j .eq. IdofGlobal(k)) then
!                AA (inode+lofsu,k+lofsv) = p_DA12(ia)*rvanka%Dmultipliers(1,2)
!                AA (inode+lofsv,k+lofsu) = p_DA21(ia)*rvanka%Dmultipliers(2,1)
!                exit
!              end if
!            end do          
!          end do
!        end if
!        
!        ! Handle the 'off-diagonal' matrices A45 and A54 if they exist
!        if (associated(p_KldA45)) then
!          ia1 = p_KldA45(idof)
!          ia2 = p_KldA45(idof+1)-1
!          do ia = ia1,ia2
!            ! Calculate:
!            !
!            !   ( du  ) = ( du  ) - (  .   .   .   .   .   .  ) ( u  )
!            !   ( dv  )   ( dv  )   (  .   .   .   .   .   .  ) ( v  )
!            !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
!            !   ( dl1 )   ( dl1 )   (  .   .   .   .  A45  .  ) ( l1 )
!            !   ( dl2 )   ( dl2 )   (  .   .   .  A54  .   .  ) ( l2 )
!            !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )
!
!            J = p_KcolA45(ia)
!            FF(inode+lofsl1) = FF(inode+lofsl1) &
!                             - rvanka%Dmultipliers(4,5)*p_DA45(ia)*p_Dvector(J+ioffsetl2)
!            FF(inode+lofsl2) = FF(inode+lofsl2) &
!                             - rvanka%Dmultipliers(5,4)*p_DA54(ia)*p_Dvector(J+ioffsetl1)
!            
!            ! Whereever we find a DOF that couples to another DOF on the 
!            ! same element, we put that to both A-blocks of our local matrix.
!            do k=1,nnvel
!              if (j .eq. IdofGlobal(k)) then
!                AA (inode+lofsl1,k+lofsl2) = p_DA45(ia)*rvanka%Dmultipliers(4,5)
!                AA (inode+lofsl2,k+lofsl1) = p_DA54(ia)*rvanka%Dmultipliers(5,4)
!                exit
!              end if
!            end do          
!          end do
!        end if
!                
!        ! Process A33 if it exists
!                
!        if (associated(p_KdiagonalA33)) then
!
!          ! Calculate:
!          !
!          !   ( du  ) = ( du  ) - (  .   .   .   .   .   .  ) ( u  )
!          !   ( dv  )   ( dv  )   (  .   .   .   .   .   .  ) ( v  )
!          !   ( dp  )   ( dp  )   (  .   .   I1  .   .   .  ) ( p  )
!          !   ( dl1 )   ( dl1 )   (  .   .   .   .   .   .  ) ( l1 )
!          !   ( dl2 )   ( dl2 )   (  .   .   .   .   .   .  ) ( l2 )
!          !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )
!          !
!          ! IEL is the pressure DOF which we have to tackle.
!          
!          daux = rvanka%Dmultipliers(3,3)
!          FF(1+lofsp) = FF(1+lofsp) &
!                      - daux*p_DA33(p_KdiagonalA33(IEL))*p_Dvector(IEL+ioffsetp)
!          AA(1+lofsp,1+lofsp) = daux*p_DA33(p_KdiagonalA33(IEL))
!
!        end if
!                
!        ! Process A66 if it exists
!                
!        if (associated(p_KdiagonalA66)) then
!
!          ! Calculate:
!          !
!          !   ( du  ) = ( du  ) - (  .   .   .   .   .   .  ) ( u  )
!          !   ( dv  )   ( dv  )   (  .   .   .   .   .   .  ) ( v  )
!          !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
!          !   ( dl1 )   ( dl1 )   (  .   .   .   .   .   .  ) ( l1 )
!          !   ( dl2 )   ( dl2 )   (  .   .   .   .   .   .  ) ( l2 )
!          !   ( dxi )   ( dxi )   (  .   .   .   .   .   I2 ) ( xi )
!          !
!          ! IEL is the pressure DOF which we have to tackle.
!          
!          daux = rvanka%Dmultipliers(6,6)
!          FF(1+lofsxi) = FF(1+lofsxi) &
!                       - daux*p_DA66(p_KdiagonalA66(IEL))*p_Dvector(IEL+ioffsetxi)
!          AA(1+lofsxi,1+lofsxi) = daux*p_DA66(p_KdiagonalA66(IEL))
!
!        end if
!                
!        ! Then subtract B*p: f_i = (f_i-Aui) - Bi pi
!        
!        ib1=p_KldB(idof)
!        ib2=p_KldB(idof+1)-1
!        do ib = ib1,ib2
!          ! Calculate:
!          !
!          !   ( du  ) = ( du  ) - (  .   .  B1   .   .   .  ) ( u  )
!          !   ( dv  )   ( dv  )   (  .   .  B2   .   .   .  ) ( v  )
!          !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
!          !   ( dl1 )   ( dl1 )   (  .   .   .   .   .  B1  ) ( l1 )
!          !   ( dl2 )   ( dl2 )   (  .   .   .   .   .  B2  ) ( l2 )
!          !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )
!
!          J = p_KcolB(ib)
!          
!          ! primal equation
!          daux = p_Dvector(j+ioffsetp) 
!          FF(inode+lofsu) = FF(inode+lofsu)-p_DB1(ib)*daux * rvanka%Dmultipliers(1,3)
!          FF(inode+lofsv) = FF(inode+lofsv)-p_DB2(ib)*daux * rvanka%Dmultipliers(2,3)
!
!          ! dual equation
!          daux2 = p_Dvector(j+ioffsetxi)
!          FF(inode+lofsl1) = FF(inode+lofsl1)-p_DB1(ib)*daux2 * rvanka%Dmultipliers(4,6)
!          FF(inode+lofsl2) = FF(inode+lofsl2)-p_DB2(ib)*daux2 * rvanka%Dmultipliers(5,6)
!          
!          ! Don't incorporate the B-matrices into AA yet; this will come later!
!        end do
!        
!        ! The mass matrix defect.
!        if (associated(p_DM)) then
!          ! We assume: multiplier of A(1,4) = multiplier of A(2,5)
!          daux = rvanka%Dmultipliers(1,4)
!          
!          ! We assume: multiplier of A(4,1) = multiplier of A(5,2)
!          daux2 = rvanka%Dmultipliers(4,1)
!          
!          ia1 = p_KldM(idof)
!          ia2 = p_KldM(idof+1)-1
!          do ia = ia1,ia2
!
!            J = p_KcolM(ia)
!
!            ! Calculate:
!            !
!            !   ( du  ) = ( du  ) - (  .   .   .  aM   .   .  ) ( u  )
!            !   ( dv  )   ( dv  )   (  .   .   .   .  aM   .  ) ( v  )
!            !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
!            !   ( dl1 )   ( dl1 )   (  .   .   .   .   .   .  ) ( l1 )
!            !   ( dl2 )   ( dl2 )   (  .   .   .   .   .   .  ) ( l2 )
!            !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )
!
!            FF(inode+lofsu) = FF(inode+lofsu)-daux*p_DM(ia)*p_Dvector(J+ioffsetl1)
!            FF(inode+lofsv) = FF(inode+lofsv)-daux*p_DM(ia)*p_Dvector(J+ioffsetl2)
!
!            ! Whereever we find a DOF that couples to another DOF on the 
!            ! same element, we put that to both A-blocks of our local matrix.
!            do k=1,nnvel
!              if (j .eq. IdofGlobal(k)) then
!                AA (inode+lofsu,k+lofsl1) = daux*p_DM(ia)
!                AA (inode+lofsv,k+lofsl2) = daux*p_DM(ia)
!                exit
!              end if
!            end do          
!          end do
!        end if
!
!        ! The defect in the coupling matrix from the primal to the dual system
!        if (associated(p_DR41)) then
!          ! Get the multipliers
!          daux = rvanka%Dmultipliers(4,1)
!          daux2 = rvanka%Dmultipliers(5,2)
!          
!          ia1 = p_KldM(idof)
!          ia2 = p_KldM(idof+1)-1
!          do ia = ia1,ia2
!
!            J = p_KcolM(ia)
!
!            ! Calculate:
!            !
!            !   ( du  ) = ( du  ) - (  .   .   .   .   .   .  ) ( u  )
!            !   ( dv  )   ( dv  )   (  .   .   .   .   .   .  ) ( v  )
!            !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
!            !   ( dl1 )   ( dl1 )   ( bR   .   .   .   .   .  ) ( l1 )
!            !   ( dl2 )   ( dl2 )   (  .  bR   .   .   .   .  ) ( l2 )
!            !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )
!
!            FF(inode+lofsl1) = FF(inode+lofsl1)-daux*p_DR41(ia)*p_Dvector(J+ioffsetu)
!            FF(inode+lofsl2) = FF(inode+lofsl2)-daux2*p_DR52(ia)*p_Dvector(J+ioffsetv)
!            
!            ! Whereever we find a DOF that couples to another DOF on the 
!            ! same element, we put that to both A-blocks of our local matrix.
!            do k=1,nnvel
!              if (j .eq. IdofGlobal(k)) then
!                AA (inode+lofsl1,k+lofsu) = daux*p_DR41(ia)
!                AA (inode+lofsl2,k+lofsv) = daux2*p_DR52(ia)
!                exit
!              end if
!            end do          
!          end do
!        end if
!
!        if (associated(p_DR51)) then
!          ! Get the multipliers
!          daux = rvanka%Dmultipliers(5,1)
!          daux2 = rvanka%Dmultipliers(4,2)
!          
!          ia1 = p_KldM(idof)
!          ia2 = p_KldM(idof+1)-1
!          do ia = ia1,ia2
!
!            J = p_KcolM(ia)
!
!            ! Calculate:
!            !
!            !   ( du  ) = ( du  ) - (  .   .   .   .   .   .  ) ( u  )
!            !   ( dv  )   ( dv  )   (  .   .   .   .   .   .  ) ( v  )
!            !   ( dp  )   ( dp  )   (  .   .   .   .   .   .  ) ( p  )
!            !   ( dl1 )   ( dl1 )   (  .  bR   .   .   .   .  ) ( l1 )
!            !   ( dl2 )   ( dl2 )   ( bR   .   .   .   .   .  ) ( l2 )
!            !   ( dxi )   ( dxi )   (  .   .   .   .   .   .  ) ( xi )
!
!            FF(inode+lofsl1) = FF(inode+lofsl1)-daux2*p_DR42(ia)*p_Dvector(J+ioffsetv)
!            FF(inode+lofsl2) = FF(inode+lofsl2)-daux*p_DR51(ia)*p_Dvector(J+ioffsetu)
!            
!            ! Whereever we find a DOF that couples to another DOF on the 
!            ! same element, we put that to both A-blocks of our local matrix.
!            do k=1,nnvel
!              if (j .eq. IdofGlobal(k)) then
!                AA (inode+lofsl2,k+lofsu) = daux*p_DR51(ia)
!                AA (inode+lofsl1,k+lofsv) = daux2*p_DR42(ia)
!                exit
!              end if
!            end do          
!          end do
!        end if
!        
!        ! THe next loop will determine the local B1, B2, D1 and D2.
!        ! We have to find in the B-matrices the column that corresponds
!        ! to our element and pressure DOF IEL - which makes it necessary
!        ! to compare the column numbers in KcolB with IEL.
!        ! Remember: The column numbers in B correspond to the pressure-DOF's
!        ! and so to element numbers. 
!        !
!        ! Btw: Each row of B has at most two entries:
!        !
!        !      IEL                              IEL
!        !   |--------|             |--------|--------|
!        !   |        |             |        |        |
!        !   |   P1   |      or     |   P2   X   P1   |
!        !   |        |             |        |        |
!        ! --|---X----|--           |--------|--------|
!        !
!        ! Either two (if the velocity DOF is an edge with two neighbouring
!        ! elements) or one (if the velocity DOF is at an edge on the boundary
!        ! and there is no neighbour).
!        do ib = ib1,ib2
!        
!          ! Calculate:
!          !
!          !   ( du  ) = ( du  ) - (  .   .   .   .   .   .  ) ( u  )
!          !   ( dv  )   ( dv  )   (  .   .   .   .   .   .  ) ( v  )
!          !   ( dp  )   ( dp  )   ( D1  D2   .   .   .   .  ) ( p  )
!          !   ( dl1 )   ( dl1 )   (  .   .   .   .   .   .  ) ( l1 )
!          !   ( dl2 )   ( dl2 )   (  .   .   .   .   .   .  ) ( l2 )
!          !   ( dxi )   ( dxi )   (  .   .   .  D1  D2   .  ) ( xi )
!          !
!          ! In AA, we simultaneously set up (locally):
!          !
!          !   (  .   .  B1   .   .   .  ) 
!          !   (  .   .  B2   .   .   .  ) 
!          !   ( D1  D2   .   .   .   .  ) 
!          !   (  .   .   .   .   .  B1  ) 
!          !   (  .   .   .   .   .  B2  ) 
!          !   (  .   .   .  D1  D2   .  ) 
!
!          if (p_KcolB(ib) .eq. IEL) then
!          
!            J = p_KcolB(ib)
!            
!            ! Get the entries in the B-matrices.
!            ! Primal equation
!            AA(inode+lofsu,1+lofsp) = p_DB1(ib) * rvanka%Dmultipliers(1,3)
!            AA(inode+lofsv,1+lofsp) = p_DB2(ib) * rvanka%Dmultipliers(2,3)
!
!            ! The same way, get DD1 and DD2.
!            ! Note that DDi has exacty the same matrix structrure as BBi and is noted
!            ! as 'transposed matrix' only because of the transposed-flag.
!            ! So we can use "ib" as index here to access the entry of DDi:
!            AA(1+lofsp,inode+lofsu) = p_DD1(ib) * rvanka%Dmultipliers(3,1)
!            AA(1+lofsp,inode+lofsv) = p_DD2(ib) * rvanka%Dmultipliers(3,2)
!
!            ! The same for the dual equation
!            AA(inode+lofsl1,1+lofsxi) = p_DB1(ib) * rvanka%Dmultipliers(4,6)
!            AA(inode+lofsl2,1+lofsxi) = p_DB2(ib) * rvanka%Dmultipliers(5,6)
!
!            AA(1+lofsxi,inode+lofsl1) = p_DD1(ib) * rvanka%Dmultipliers(6,4)
!            AA(1+lofsxi,inode+lofsl2) = p_DD2(ib) * rvanka%Dmultipliers(6,5)
!
!            ! Build the pressure entry in the local defect vector:
!            !   f_i = (f_i-Aui) - D_i pi
!            ! or more precisely (as D is roughly B^T):
!            !   f_i = (f_i-Aui) - (B^T)_i pi
!            FF(1+lofsp) = FF(1+lofsp) &
!                        - AA(1+lofsp,inode+lofsu)*p_Dvector(idof+ioffsetu) &
!                        - AA(1+lofsp,inode+lofsv)*p_Dvector(idof+ioffsetv)
!          
!            ! The same for the dual pressure
!            FF(1+lofsxi) = FF(1+lofsxi) &
!                         - AA(1+lofsxi,inode+lofsl1)*p_Dvector(idof+ioffsetl1) &
!                         - AA(1+lofsxi,inode+lofsl2)*p_Dvector(idof+ioffsetl2)
!          
!            ! Quit the loop - the other possible entry belongs to another 
!            ! element, not to the current one
!            exit
!          end if
!        end do ! ib
!        
!      end do ! inode
!    
!      ! Now we make a defect-correction approach for this system:
!      !
!      !    x_new  =  x  +  P( \omega C^{-1} (f~ - A~ x) )
!      !                                     -----------
!      !                                        =d~
!      !
!      ! Here the 'projection' operator simply converts the small
!      ! preconditioned defect (\omega C^{-1} d~) to a 'full' defect
!      ! of the same size as x - what is easy using the number of
!      ! the DOF's on the element.
!      !
!      ! For C, we use our local AA, i.e. applying C^{-1} means to
!      ! solve the local system AA dd = FF for dd. The local defect dd is then
!      ! added back to the global solution vector.
!      
!      call DGESV (nnld, 1, AA, nnld, Ipiv, FF, nnld, ilapackInfo)
!      
!      if (ilapackInfo .eq. 0) then
!        
!        ! Ok, we got the update vector in FF. Incorporate this now into our
!        ! solution vector with the update formula
!        !
!        !  x_{n+1} = x_n + domega * y
!        
!        do inode=1,nnvel
!          ! Update of the primal velocity vectors
!          p_Dvector(idofGlobal(inode)+ioffsetu) &
!            = p_Dvector(idofGlobal(inode)+ioffsetu) + domega * FF(inode+lofsu)
!          p_Dvector(idofGlobal(inode)+ioffsetv) &
!            = p_Dvector(idofGlobal(inode)+ioffsetv) + domega * FF(inode+lofsv)
!
!          ! Update of the dual velocity vectors
!          p_Dvector(idofGlobal(inode)+ioffsetl1) &
!            = p_Dvector(idofGlobal(inode)+ioffsetl1) + domega * FF(inode+lofsl1)
!          p_Dvector(idofGlobal(inode)+ioffsetl2) &
!            = p_Dvector(idofGlobal(inode)+ioffsetl2) + domega * FF(inode+lofsl2)
!        end do
!        
!        p_Dvector(iel+ioffsetp) = p_Dvector(iel+ioffsetp) + &
!                                  domega * FF(1+lofsp)
!
!        p_Dvector(iel+ioffsetxi) = p_Dvector(iel+ioffsetxi) + &
!                                  domega * FF(1+lofsxi)
!      
!      else if (ilapackInfo .lt. 0) then
!        
!        call output_line (&
!            'LAPACK(DGESV) solver failed! Error code: '//sys_siL(ilapackInfo,10),&
!            OU_CLASS_ERROR,OU_MODE_STD,'vanka_2DNSSOCQ1TQ0fullCoupConf')
!        
!      end if
!
!      ! (ilapackInfo > 0) May happen in rare cases, e.g. if there is one element on the
!      ! coarse grid with all boundaries = Dirichlet.
!      ! In this case, nothing must be changed in the vector!
!    
!    end do ! iel
!
!  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc2doptc_getLogFiles (slogfile,serrorfile)
  
!<description>
  ! Temporarily reads the output DAT file to get the names of the output
  ! files.
!</description>

!<output>
  ! Name of the message log file.
  character(LEN=*), intent(OUT) :: slogfile
  
  ! Name of the error log file.
  character(LEN=*), intent(OUT) :: serrorfile
!</output>

!</subroutine>

    type(t_parlist) :: rparlist

    ! Init parameter list that accepts parameters for output files
    call parlst_init (rparlist)

    ! Read parameters that configure the output
    call parlst_readfromfile (rparlist, trim(DIR_DATA)//'/output.dat')
    
    ! Now the real initialisation of the output including log file stuff!
    call parlst_getvalue_string (rparlist,'GENERALOUTPUT',&
        'smsgLog',slogfile,"",bdequote=.true.)

    call parlst_getvalue_string (rparlist,'GENERALOUTPUT',&
        'serrorLog',serrorfile,"",bdequote=.true.)
    
    ! That temporary parameter list is not needed anymore.
    call parlst_done (rparlist)
    
    end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc2doptc_evalParameters ()
  
!<description>
  ! Evaluates command line parameters.
  ! In the current implementation, command line parameters are passed as
  ! text file. This routine searches in the main directory for a file
  ! "cmdline.dat". If this file is found, it's opened and evaluated.
  ! Every line may contain a command line parameter in the form of
  ! a DAT file (name=value pairs).
  !
  ! Supported command line parameters:
  !   "datdirectory = [Directory, where DAT files can be found]"
!</description>

!</subroutine>

    ! local variables
    type(t_parlist) :: rparamList
    logical :: bexists

    ! Figure out if the file exists.
    inquire(file='./cmdline.dat', exist=bexists)
    
    if (bexists) then
      ! Read the file
      call parlst_init (rparamList)
      call parlst_readfromfile (rparamList, './cmdline.dat')
      
      ! Evaluate parameters
      call parlst_getvalue_string ( &
          rparamList, "","datdirectory", DIR_DATA, DIR_DATA,bdequote=.true.)

      call parlst_done (rparamList)
      
    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine main_getDat (rparamList)
  
!<description>
  ! Reads in all DAT files into the parameter list rparlist
!</description>

!<inputoutput>
  ! The parameter list where the values of the DAT files should be stored.
  ! The structure must have been initialised, the parameters are just added
  ! to the list.
  type(t_parlist), intent(INOUT) :: rparamList
!</inputoutput>

!</subroutine>

    logical :: bexists
    character(LEN=SYS_STRLEN) :: smaster
    
    ! Check if a command line parameter specifies the master.dat file.
    call sys_getcommandLineArg(1,smaster,sdefault=trim(DIR_DATA)//'/master.dat')

    ! Read the file 'master.dat'.
    ! If that does not exist, try to manually read files with parameters from a
    ! couple of files.
    inquire(file=smaster, exist=bexists)
    
    if (bexists) then
      ! Read the master file. That either one contains all parameters or
      ! contains references to subfiles with data.
      call parlst_readfromfile (rparamList, trim(smaster),trim(DIR_DATA))
    else
      call output_line("Master file not found: "//trim(smaster),&
          OU_CLASS_WARNING,ssubroutine='main_getDat')
      call output_line("Reading standard parameters.",&
          OU_CLASS_WARNING,ssubroutine='main_getDat')
    
      ! Each 'readfromfile' command adds the parameter of the specified file 
      ! to the parameter list.
      call parlst_readfromfile (rparamList, trim(DIR_DATA)//'/main.dat')
      call parlst_readfromfile (rparamList, trim(DIR_DATA)//'/bdconditions.dat')
      call parlst_readfromfile (rparamList, trim(DIR_DATA)//'/discretisation.dat')
      call parlst_readfromfile (rparamList, trim(DIR_DATA)//'/flows.dat')
      call parlst_readfromfile (rparamList, trim(DIR_DATA)//'/linsol.dat')
      call parlst_readfromfile (rparamList, trim(DIR_DATA)//'/forwardsolver.dat')
      call parlst_readfromfile (rparamList, trim(DIR_DATA)//'/optcontrol.dat')
      call parlst_readfromfile (rparamList, trim(DIR_DATA)//'/output.dat')
      call parlst_readfromfile (rparamList, trim(DIR_DATA)//'/paramtriang.dat')
      call parlst_readfromfile (rparamList, trim(DIR_DATA)//'/postprocessing.dat')
      call parlst_readfromfile (rparamList, trim(DIR_DATA)//'/spacetimesolver.dat')
      call parlst_readfromfile (rparamList, trim(DIR_DATA)//'/timediscr.dat')
    end if
  
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine cc2doptccalculate (rsettings,rparlist)
  
!<description>
  ! This is a 'separated' Navier-Stokes solver for solving a Navier-Stokes
  ! problem. The different tasks of the problem are separated into
  ! subroutines. The problem uses a problem-specific structure for the 
  ! communication: All subroutines add their generated information to the
  ! structure, so that the other subroutines can work with them.
  ! (This is somehow a cleaner implementation than using a collection!).
  ! For the communication to callback routines of black-box subroutines
  ! (matrix-assembly), a collection is used.
  !
  ! The following tasks are performed by the subroutines:
  !
  ! 1.) Read in parametrisation
  ! 2.) Read in triangulation
  ! 3.) Set up RHS
  ! 4.) Set up matrix
  ! 5.) Create solver structure
  ! 6.) Solve the problem
  ! 7.) Write solution to GMV file
  ! 8.) Release all variables, finish
!</description>

!<input>
  ! Basic program settings
  type(t_settings_main), intent(in) :: rsettings
  
  ! Parameter list containing the DAT file parameters
  type(t_parlist), intent(in) :: rparlist
!</input>

!</subroutine>

    ! Local variables
    !
    ! The following structure configures our solver.
    ! We allocate this dynamically since the structure may be too large
    ! for the stack.
    type(t_settings_optflow), pointer :: p_rsettingsSolver
    
    ! Basic space discretisation settings
    type(t_settings_discr) :: rsettingsSpaceDiscr
    
    ! Nonlinear solve structure representing the solver.
    type(t_nlstsolver), pointer :: p_rnlstsolver
    
    ! Global solution, temp vector
    type(t_spacetimevector) :: rsolution,rtemp,rrhsDiscrete

    ! Analytic solution defining the RHS of the equation.
    type(t_anSolution) :: rrhs

    ! Postprocessing data.    
    type(t_optcPostprocessing) :: rpostproc
    
    ! Some timers
    type(t_timer) :: rtotalTime,rinitTime,rsolverTime,rtimePostProc,rstartvectime

    ! Ok, let us start. 
    !
    ! Initialise the external storage management.
    
    call exstor_init (999,100)
    !CALL exstor_attachDirectory('./ff2storage')
    
    ! Ok, parameters are read in.
    ! Print the parameters to the terminal.
    if (rsettings%routput%ioutputInit .ge. 2) then
      call parlst_info (rparlist)
      call output_separator (OU_SEP_EQUAL)
    end if
    
    ! Measure the total time.
    call stat_clearTimer (rtotalTime)
    call stat_startTimer (rtotalTime)
    
    ! Initialise the settings of the solver,
    ! allocate all template matrices etc.
    allocate(p_rsettingsSolver)
    allocate(p_rnlstsolver)
    
    call stat_clearTimer (rinitTime)
    call stat_startTimer (rinitTime)

    call output_line ("Initialising solver structures.")
    call init_initStandardSolver (rparlist,rsettings,p_rsettingsSolver,p_rnlstsolver,&
        rpostproc,rrhs,rsettings%routput%ioutputInit)

    ! Discretise the RHS according to the time stepping scheme.
    call output_lbrk()
    call output_line ("Discretising RHS.")
    call sptivec_initVector (rrhsdiscrete,&
        p_rsettingsSolver%rtimeHierarchy%p_rtimeLevels(p_rsettingsSolver%rtimeHierarchy%nlevels),&
        p_rsettingsSolver%rfeHierPrimalDual% &
        p_rfeSpaces(p_rsettingsSolver%rfeHierPrimalDual%nlevels)%p_rdiscretisation)
    call sptivec_clearVector (rrhsdiscrete)
    call init_discretiseRHS (p_rsettingsSolver,rrhs,rrhsDiscrete)
    
    ! DEBUG!!!
    !call sptivec_saveToFileSequence(rrhsdiscrete,"(""rhs"// &
    !    trim(sys_siL(p_rsettingsSolver%rtimeHierarchy%nlevels,10))// &
    !    ".txt."",I5.5)",.true.)
    
    ! Create a start vector for the solver.
    call output_lbrk()
    call output_line ("Initialising start vector.")
    call sptivec_initVector (rsolution,&
        p_rsettingsSolver%rtimeHierarchy%p_rtimeLevels(p_rsettingsSolver%rtimeHierarchy%nlevels),&
        p_rsettingsSolver%rfeHierPrimalDual% &
        p_rfeSpaces(p_rsettingsSolver%rfeHierPrimalDual%nlevels)%p_rdiscretisation)
    call init_getSpaceDiscrSettings (rparlist,rsettings%ssectionDiscrSpace,&
        rsettingsSpaceDiscr)

    call stat_clearTimer (rstartvectime)
    call stat_startTimer (rstartvectime)

    call init_initStartVector(p_rsettingsSolver,rsettingsSpaceDiscr,&
        p_rsettingsSolver%rspaceTimeHierPrimalDual%nlevels,&
        rparlist,rsettings%ssectionSpaceTimePreprocessing,&
        p_rsettingsSolver%rinitialCondition,rsolution,rrhsdiscrete,rsettings%routput%ioutputInit)
    
    call stat_stopTimer (rstartvectime)
    
    if (rsettings%routput%ioutputInit .ge. 1) then
      call output_lbrk ()
      call output_line ("Time for creation of the start vector = "//&
          sys_sdL(rstartvectime%delapsedReal,10))
    end if
    
    ! Implement the initial condition to the discrete RHS.
    call init_implementInitCondRHS (p_rsettingsSolver,rsolution,rrhsdiscrete)

    ! Create a temp vector    
    call sptivec_initVector (rtemp,&
        p_rsettingsSolver%rtimeHierarchy%p_rtimeLevels(p_rsettingsSolver%rtimeHierarchy%nlevels),&
        p_rsettingsSolver%rfeHierPrimalDual% &
        p_rfeSpaces(p_rsettingsSolver%rfeHierPrimalDual%nlevels)%p_rdiscretisation)

    call stat_stopTimer (rinitTime)

    call output_separator (OU_SEP_EQUAL)
    call output_line ("Time for initialisation            = "//&
        sys_sdL(rinitTime%delapsedReal,10))
    call output_separator (OU_SEP_EQUAL)
        
    ! Solve the system
    call stat_clearTimer (rsolverTime)
    call stat_startTimer (rsolverTime)
    call nlstslv_solve (p_rsettingsSolver,p_rnlstsolver,rpostproc,rsolution,rrhsdiscrete,rtemp)
    call stat_stopTimer (rsolverTime)
    
    call output_separator (OU_SEP_EQUAL)

    ! Pipe the solution through our postprocessing routines
    call output_line ("Postprocessing of the final solution...")
    call stat_clearTimer (rtimePostProc)
    call stat_startTimer (rtimePostProc)
    call optcpp_postprocessSpaceTimeVec (rpostproc,rsolution,rrhsDiscrete,&
        p_rsettingsSolver%rsettingsOptControl,p_rsettingsSolver)    
    call stat_stopTimer (rtimePostProc)
    
    ! Sum up the time for the postprocesing during the simulation
    call stat_addTimers (p_rnlstsolver%rtimePostprocessing,rtimePostProc)

    call output_separator (OU_SEP_EQUAL)
    
    ! Print out statistics about our solver.
    call stnlsinit_printSolverStatistics (p_rnlstsolver)
    
    ! Release all data
    call sptivec_releaseVector (rtemp)
    call sptivec_releaseVector (rsolution)
    call sptivec_releaseVector (rrhsdiscrete)
    call ansol_done(rrhs)
    call init_doneStandardSolver (p_rsettingsSolver,p_rnlstsolver,rpostproc)
    
    call stat_stopTimer (rtotalTime)

    call output_separator (OU_SEP_EQUAL)
    call output_line ("Time for initialisation            = "//&
        sys_sdL(rinitTime%delapsedReal,10))
    call output_line ("Time for creation of the start vec = "//&
        sys_sdL(rstartvectime%delapsedReal,10))
    call output_line ("Time for postprocessing            = "//&
        sys_sdL(rtimePostProc%delapsedReal,10))
    call output_line ("Time for solving                   = "//&
        sys_sdL(rsolverTime%delapsedReal,10))
    call output_line ("Total time                         = "//&
        sys_sdL(rtotalTime%delapsedReal,10))
    call output_separator (OU_SEP_EQUAL)
    
    deallocate(p_rsettingsSolver)
    deallocate(p_rnlstsolver)
    
    ! Information about external storage usage
    call output_lbrk ()
    call exstor_info (bprintHandles=.true.)
    
    ! Clean up the external storage management
    call exstor_done ()
    
  end subroutine

  ! ***************************************************************************

  subroutine cc2doptcmain
    
    ! Program parameters
    type(t_parlist) :: rparlist
    type(t_settings_main) :: rsettings
    
    ! The very first thing in every application: 
    ! Initialise system-wide settings:
    call system_init()
    
    ! Read the program parameters.
    call parlst_init (rparlist)
    call main_getDat (rparlist)
    
    ! Read the basic parameter settings from the "MAIN" section.
    call smain_initMainParams (rsettings,rparlist,"MAIN")
    
    ! Initialise log file for output.
    call output_init (rsettings%routput%smsgLog,rsettings%routput%serrorLog)
    OU_LINE_LENGTH = 132
    cdefaultDateTimeLogPolicy = OU_DTP_ADDDATETIME
    cdatetimeLogFormat = 1
    
    ! Now we can really start!
    !
    ! Initialise the storage management: 
    call storage_init(999, 100)
    
    ! Initialise the parser
    call fparser_init ()
    
    ! Call the problem to solve. 
    call output_lbrk ()
    call output_line ('Calculating cc2doptc-Problem')
    call output_separator (OU_SEP_MINUS)
    
    call cc2doptccalculate (rsettings,rparlist)

    ! Release the parser
    call fparser_done ()

    ! Print out heap statistics - just to check if everything
    ! is cleaned up.
    ! This should display 'Handles in use=0' and 'Memory in use=0'!
    call output_lbrk ()
    call storage_info(.true.)
    
    ! Clean up the storage management, parameter list, finish
    call storage_done()
    call parlst_done (rparlist)
    
  end subroutine

end module
