!#########################################################################
!# ***********************************************************************
!# <name> sortstrategy </name>
!# ***********************************************************************
!#
!# <purpose>
!# This module contains different routines to calculate resorting
!# strategies for scalar vectors.
!# A sorting strategy is simply a permutation of the numbers 1..NEQ
!# how to permute a scalar vector and its inverse permutation. 
!# Both permutations are usually assigned as one large array
!# to a (scalar) vector to indicate how it's resorted.
!#
!# Which information are necessary for the calculation of a permutation
!# is completely algorithm dependent. Therefore, this module has full
!# access to the whole (spatial) discretisation including triangulation
!# the domain, matrix data and more. It depends on the algorithm which
!# information is actually used.
!#
!# The following routines can be found in this module:
!#
!# 1.) sstrat_calcCuthillMcKee
!#     -> Calculate column numbering using the Cuthill McKee algorithm
!#
!# 2.) sstrat_calcRevCuthillMcKee
!#     -> Calculate column numbering using reverse Cuthill-McKee algorithm.
!#
!# 3.) sstrat_calcXYZsorting
!#     -> Calculates rowwise renumbering based on Cartesian coordinates
!#
!# 4.) sstrat_calcFEASTsorting
!#     -> Calculates FEAST rowwise renumbering based cell adjacencies
!#
!# 5.) sstrat_calcStochastic
!#     -> Calculates stoastic renumbering (i.e. a random permutation)
!#
!# 6.) sstrat_calcHierarchical
!#     -> Calculates a renumbering strategy based on element patches
!#        in a level hierarchy
!#
!# Auxiliary routines:
!#
!# 1.) sstrat_calcColNumberingCM7
!#     -> Calculate Cuthill-Mc-Kee resorting strategy for a structure-7
!#        matrix
!#
!# 2.) sstrat_calcColNumberingCM9
!#     -> Calculate Cuthill-Mc-Kee resorting strategy for a structure-9
!#        matrix
!#
!# 3.) sstrat_calcPermutationCM
!#     -> Auxiliary routine to calculate a permutation from the output
!#        of sstrat_calcColNumberingCM7 or sstrat_calcColNumberingCM9.
!#
!# 4.) sstrat_calcInversePermutation
!#     -> Calculates an inverse permutation
!#
!#
!# The renumbering routines in this module always calculate a
!# permutation as well as its inverse permutation. The calculated
!# permutation Ipermutation given as a parameter to the routines
!# has length 2*N for N numbers. It contains the permutation in the
!# first N and its inverse at the second N entries. Here, the
!# exact meaning of the word 'permutation' is as follows:
!#
!#  Ipermutation (position in sorted vector) = position in unsorted vector.
!#  Ipermutation (N+position in unsorted vector) = position in sorted vector.
!#
!# </purpose>
!#########################################################################

module sortstrategy

  use fsystem
  use storage
  use genoutput
  use linearalgebra
  use spatialdiscretisation
  use element
  use triangulation
  use linearsystemscalar

  implicit none
  
!<constants>

!<constantblock description="Sort strategy identifiers.">

  ! No sort strategy; this must be =0!
  integer, parameter :: SSTRAT_UNSORTED     = 0

  ! Cuthill-McKee sort strategy
  integer, parameter :: SSTRAT_CM           = 1

  ! Reverse Cuthill-McKee sort strategy
  integer, parameter :: SSTRAT_RCM          = 2
  
  ! Row-wise sorting for point coordinate. 
  ! (As calculated by sstrat_calcXYZsorting with idirection=0.)
  ! Only for special type of discretisations ($Q_1$, $\tilde Q_1$), where 
  ! the DOF's can be identified with X/Y/Z coordinates.
  ! Coincides with the sorting strategy of FEAST for simple-type domains
  ! like the unit square.
  integer, parameter :: SSTRAT_XYZCOORD     = 3
  
  ! Column-wise sorting for point coordinate. 
  ! (As calculated by sstrat_calcXYZsorting with idirection=1.)
  ! Only for special type of discretisations ($Q_1$, $\tilde Q_1$), where 
  ! the DOF's can be identified with X/Y/Z coordinates.
  integer, parameter :: SSTRAT_ZYXCOORD     = 4
  
  ! General FEAST renumbering.
  ! The DOF's are numbered rowwise, independent of the geometrical
  ! structure of the domain.
  ! Only for special type of discretisations ($Q_1$) and tensor product meshes.
  integer, parameter :: SSTRAT_FEAST        = 5
  
  ! Stochastic renumbering / Random permutation.
  ! The permutation is completely random.
  integer, parameter :: SSTRAT_STOCHASTIC   = 6

  ! Hierarchical renumbering: The permutation is calculated by a sequence of meshes,
  ! regularly refined.
  integer, parameter :: SSTRAT_HIERARCHICAL = 7
  
!</constantblock>

!</constants>

!<types>

!<typeblock description="Level hierarchy structure for te hierarchically calculated permutation">

  ! A local structure for the hierarchical sorting strategy.
  ! Represents one level in the hierarchy.
  type t_levelHirarchy
  
    ! Pointer to the refinement-patch array of the level
    integer, dimension(:), pointer :: p_IrefinementPatch
    
    ! Pointer to the refinement-patch-index array of the level
    integer, dimension(:), pointer :: p_IrefinementPatchIdx

    ! Whether the corresponding discretisation on that level is uniform or not.
    logical :: bisUniform
    
    ! Element type; only valid if the corresponding discretisation is uniform.
    integer(I32) :: celement
    
    ! Pointer to the identifier for the element distribution of an element.
    ! Only valid if the corresponding discretisation is not uniform.
    integer, dimension(:), pointer :: p_IelementDistr
    
  end type

!</typeblock>

!</types>

  private :: t_levelHirarchy
  
  interface sstrat_calcCuthillMcKee
    module procedure sstrat_calcCuthillMcKee_p
    module procedure sstrat_calcCuthillMcKee_h
  end interface

contains

  ! ***************************************************************************

!<subroutine>

  subroutine sstrat_calcStochastic (Ipermutation)
  
  !<description>
    ! Generates a random permutation.
    !
    ! The used algorithm is Knuth shuffle.
    ! (See e.g. [Knuth, 1969, 1998, The Art of Computer Programming vol. 2, 3rd ed., 
    !  145�146. ISBN 0-201-89684-2] or http://en.wikipedia.org/wiki/Knuth_shuffle])
  !</description>
    
  !<output>
    ! The permutation vector for sorting and its inverse.
    ! With NEQ=NEQ(matrix):
    !   Ipermutation(1:NEQ)       = permutation,
    !   Ipermutation(NEQ+1:2*NEQ) = inverse permutation.
    integer, dimension(:), intent(OUT) :: Ipermutation
  !</output>    

  !</subroutine>
  
    real(DP) :: d
    integer :: i,k,n
    integer :: j

    ! Fill the array with 1,2,3,...
    n = size(Ipermutation) / 2
    
    do i = 1,n
      Ipermutation(i) = i
    end do
    
    ! Loop through the array and randomly swap element i with element [i,...,n]
    do i = 1,n-1
      call random_number(d)
      k = i + int(real(n-i,DP) * d + 0.5_DP)
      
      j = Ipermutation(i)
      Ipermutation(i) = Ipermutation(k)
      Ipermutation(k) = j
    end do

    ! Calculate the inverse permutation, that's it.
    call sstrat_calcInversePermutation (Ipermutation(1:N), Ipermutation(N+1:) )

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sstrat_calcCuthillMcKee_h (rmatrix, h_Ipermutation)
  
!<description>
  ! Computes a column renumbering strategy using the algorithm
  ! of Cuthill-McKee. The algorithm acceps a scalar matrix rmatrix and
  ! uses its structure to calculate the renumbering. The result
  ! Ipermutation then receives the permutation and its inverse.
!</description>
    
!<input>
  ! Matrix which should be used to calculate the renumbering strategy
  type(t_matrixScalar), intent(IN) :: rmatrix
!</input>
    
!<output>
  ! A storage handle to the permutation array.
  integer, intent(OUT) :: h_Ipermutation
!</output>    

!</subroutine>

  ! The permutation array
  integer, dimension(:), pointer :: p_Ipermutation

    ! Allocate an array for holding the resorting strategy.
    call storage_new ('sstrat_calcCuthillMcKee_h', 'Ipermutation', &
          2*rmatrix%NEQ, ST_INT, h_Ipermutation, ST_NEWBLOCK_ZERO)
    call storage_getbase_int(h_Ipermutation,p_Ipermutation)
    
    ! Call the algorithm itself
    call sstrat_calcCuthillMcKee_p(rmatrix, p_Ipermutation)
  
  end subroutine  

  ! ***************************************************************************

!<subroutine>

  subroutine sstrat_calcRevCuthillMcKee (rmatrix, h_Ipermutation)
  
!<description>
  ! Computes a column renumbering strategy using the of reverse Cuthill-McKee
  ! algorithm. The algorithm acceps a scalar matrix rmatrix and uses its
  ! structure to calculate the renumbering.
!</description>
    
!<input>
  ! Matrix which should be used to calculate the renumbering strategy
  type(t_matrixScalar), intent(IN) :: rmatrix
!</input>
    
!<output>
  ! A storage handle to the permutation array.
  integer, intent(OUT) :: h_Ipermutation
!</output>    

!</subroutine>
  
  ! The permutation array
  integer, dimension(:), pointer :: p_Ipermutation
  integer :: i,x,n

    ! Calculate a 'forward' Cuthill-McKee permutation
    call sstrat_calcCuthillMcKee(rmatrix, h_Ipermutation)

    ! Get the calculated permutation
    call storage_getbase_int(h_Ipermutation, p_Ipermutation)
    
    ! Reverse the ordering
    n = size(p_Ipermutation) / 2
    do i = 1, n / 2
      x = p_Ipermutation(i)
      p_Ipermutation(i) = p_Ipermutation(n-i+1)
      p_Ipermutation(n-i+1) = x
    end do
    
    ! Recalculate inverse permuation
    call sstrat_calcInversePermutation(p_Ipermutation(1:n), &
                                       p_Ipermutation(n+1:2*n))
    
  end subroutine  

  ! ***************************************************************************

!<subroutine>

  subroutine sstrat_calcCuthillMcKee_p (rmatrix,Ipermutation)
  
  !<description>
    ! Computes a column renumbering strategy using the algorithm
    ! of Cuthill-McKee. The algorithm acceps a scalar matrix rmatrix and
    ! uses its structure to calculate the renumbering. The result
    ! Ipermutation then receives the permutation and its inverse.
  !</description>
    
  !<input>
    ! Matrix which should be used to calculate the renumbering strategy
    type(t_matrixScalar), intent(IN) :: rmatrix
  !</input>
    
  !<output>
    ! The permutation vector for sorting and its inverse.
    ! With NEQ=NEQ(matrix):
    !   Ipermutation(1:NEQ)       = permutation,
    !   Ipermutation(NEQ+1:2*NEQ) = inverse permutation.
    integer, dimension(2*rmatrix%neq), intent(OUT) :: Ipermutation
  !</output>    

  !</subroutine>
  
  ! local variables
  integer :: h_Ideg,h_IcolTmp
  integer, dimension(:), pointer :: p_Ideg
  integer, dimension(:), pointer :: p_Kld,p_IcolTmp, p_Kcol,p_Kdiag
  integer :: NEQ
  
  NEQ = rmatrix%NEQ
  
  ! Currently, only matrix structure 7 and 9 are supported:
  select case (rmatrix%cmatrixFormat)
  case (LSYSSC_MATRIX9)
    ! At first, duplicate KCOL and also get a temporary Ideg array
    h_IcolTmp = ST_NOHANDLE
    call storage_copy(rmatrix%h_Kcol,h_IcolTmp)
    call storage_new('sstrat_calcCuthillMcKee', 'KDEG', rmatrix%NEQ, &
                     ST_INT, h_Ideg, ST_NEWBLOCK_NOINIT)
    call storage_getbase_int(h_IcolTmp, p_IcolTmp)
    call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
    call storage_getbase_int(h_Ideg, p_Ideg)
    call lsyssc_getbase_Kld(rmatrix, p_Kld)
    call lsyssc_getbase_Kdiagonal(rmatrix, p_Kdiag)
    
    ! Calculate the strategy, calculate p_IcolTmp
    ! BETA STATUS, ROUTINE NOT TESTED!!!
    call sstrat_calcColNumberingCM9 (p_Kld, p_Kcol, p_Kdiag,&
                                     p_IcolTmp, p_Ideg, NEQ, NEQ)
    
    ! Use p_IcolTmp to calculate the actual resorting permutation
    ! and its inverse. Clear the target vector before calling the
    ! calculation routine, as we are creating a new permutation!
    call lalg_clearVectorInt(Ipermutation(1:NEQ*2))
    call sstrat_calcPermutationCM (p_Kld, p_IcolTmp, NEQ, &
                                   Ipermutation(1:NEQ), Ipermutation(NEQ+1:NEQ*2))

    ! Release temp data.
    call storage_free (h_Ideg)
    call storage_free (h_IcolTmp)    

  case (LSYSSC_MATRIX7)
    ! At first, duplicate KCOL and also get a temporary Ideg array
    h_IcolTmp = ST_NOHANDLE
    call storage_copy(rmatrix%h_Kcol,h_IcolTmp)
    call storage_new('sstrat_calcCuthillMcKee', 'KDEG', rmatrix%NEQ, &
                     ST_INT, h_Ideg, ST_NEWBLOCK_NOINIT)
    call storage_getbase_int(h_IcolTmp, p_IcolTmp)
    call lsyssc_getbase_Kcol(rmatrix, p_Kcol)
    call storage_getbase_int(h_Ideg, p_Ideg)
    call lsyssc_getbase_Kld(rmatrix, p_Kld)
    
    ! Calculate the strategy, calculate p_IcolTmp
    call sstrat_calcColNumberingCM7 (p_Kld, p_Kcol, p_IcolTmp, p_Ideg, NEQ, NEQ)
    
    ! Use p_IcolTmp to calculate the actual resorting permutation
    ! and its inverse. Clear the target vector before calling the
    ! calculation routine, as we are creating a new permutation!
    call lalg_clearVectorInt(Ipermutation(1:NEQ*2))
    call sstrat_calcPermutationCM (p_Kld, p_IcolTmp, NEQ, &
                                   Ipermutation(1:NEQ), Ipermutation(NEQ+1:NEQ*2))

    ! Release temp data.
    call storage_free (h_Ideg)
    call storage_free (h_IcolTmp)    
  
  case DEFAULT
    print *,'sstrat_calcCuthillMcKee: Unsupported matrix format'
    call sys_halt()
  end select
    

  end subroutine     

  ! ***************************************************************************

!<subroutine>
  subroutine sstrat_calcColNumberingCM7 (Ild, Icol, Icon, Ideg, neq, ndeg)
  
    !<description>
    ! Purpose: Cuthill McKee matrix renumbering
    !          Calculate column numbering
  
    ! The algorithm of Cuthill McKee interprets the system matrix as
    ! adjacense matrix. Every Row denotes a node in the corresponing graph.
    ! In the first step this function sorts the columns in every row
    ! for increasing degree of the nodes in the matrix.
    ! The matrix must have symmetric structure!
    ! (For FE matrices this is always the case...)
    !
    ! Matrix storage technique 7 version.
    !</description>
  
  !<input>
  
    ! Number of equations
    integer,intent(IN)                    :: neq

    ! Maximum number of entries != 0 in every row of the matrix
    integer, intent(IN)                   :: ndeg    
   
    ! Row description of matrix
    integer, dimension(neq+1), intent(IN) :: Ild
  
    ! Column description of matrix
    integer, dimension(:), intent(IN)     :: Icol
    
  !</input>

  !<inputoutput>
    ! Auxiliary vector; must be at least as long as the
    ! maximum number of entries != 0 in every row of the matrix
    integer, dimension(ndeg), intent(INOUT) :: Ideg

    ! Auxiliary vector; the column numbers of KCOL are assigned to this in
    ! the order of increasing degree. When calling the routine the user
    ! must copy the content of KCOL to this! These values are then
    ! resorted.
    integer, dimension(:), intent(INOUT)  :: Icon
  !</inputoutput>

!</subroutine>

    ! local variables
    integer :: ieq, idegIdx, ildIdx, idegIdx1, idegIdx2
    integer :: idegMin, iidxMin
  
    ! Clear auxiliary vector
    Ideg(1:ndeg) = 0
    !DO idegIdx=1, ndeg
    !  Ideg(idegIdx) = 0
    !END DO
    
    ! Later in the algorithm the column numbers of the non-diagonal
    ! entries are changed. The numbers of the diagonal entries are not
    ! touched. Therefore we copy the column numbers of the diagonal
    ! entries directly.
    do ieq=1, neq
      Icon(Ild(ieq)) = Icol(Ild(ieq))
    end do
    
    ! Loop about all rows in the matrix: ieq = current row.
    ! Every row corresponds to an entry in the solution vector = a node
    ! in the graph of the matrix.
    do ieq=1, neq
    
      ! Copy the column numbers of the non-diagonal entries in the current
      ! column to the auxiliary vector Ideg. Ideg contains always the
      ! following nodes of the current node ieq, which are not processed yet.
      ! The entries of the vector are set to 0 one after the other in the
      ! order of the degree of the following nodes.
      do ildIdx=Ild(ieq)+1, Ild(ieq+1)-1
        Ideg(ildIdx-Ild(ieq)) = Icol(ildIdx)
      end do
      
      ! Loop about every column in the current row. The entries in the
      ! row (=column numbers) of the matrix represent the node numbers of
      ! those nodes in the graph of the matrix, which are connected to the
      ! current node ieq.
      
      ! The algorithm now has to sort these nodes for increasing degree.
      ! For this purpose the node with the smallest degree is searched and
      ! deleted from Ideg, then the node with the 2nd smallest degree
      ! and so on.
      ! We don't have many nodes, so we use a simple O(n^2) algorithm here
      ! which simply loops over all nodes.
      ! It's only necessary to search in the non-diagonal entries,
      ! because we want to sort the adjacent nodes of the current one
      ! for increasing degree.
      do idegIdx1=1, Ild(ieq+1)-Ild(ieq)-1
      
        ! iidxMin will receive the index in the Ideg-array of the node with
        ! the smallest degree. We start with node 1 in the current column:
        iidxMin=1

        ! The variable idegMin always contains the degree of the node, which
        ! is described by the index iidxMin

        ! If Ideg(iidxMin)=0, this node has already been processed. In this
        ! case idegMin is set to neq - the largest upper bound which implies
        ! that iidxMin is later on surely replaced by a node with a smaller
        ! degree.

        ! Otherwise idegMin receives the degree of the node described by
        ! iidxMin, which is calculated by the number of elements in Icol,
        ! i.e. the difference of the indices in the Ild array.
        if (Ideg(iidxMin) .eq. 0) then
          idegMin = neq
        else
          idegMin = Ild(Ideg(iidxMin)+1)-Ild(Ideg(iidxMin))-1
        endif

        ! Compare iidxMin with every node in that line to find that with
        ! minimum degree:
        do idegIdx2=1, Ild(ieq+1)-Ild(ieq)-1

          ! If Ideg(idegIdx2)=0, idegIdx2 has already been processed. Set
          ! idegIdx=neq ; here a lower bound to prevent the current node
          ! iidxMin from being replaced by idegIdx2.
          
          ! Otherwise set idegIdx to the degree of the node with index
          ! idegIdx2
          if (Ideg(idegIdx2) .eq. 0) then
            idegIdx = neq
          else
            idegIdx = Ild(Ideg(idegIdx2)+1)-Ild(Ideg(idegIdx2))-1
          endif
                
          ! If now idegIdx=grad(iidxMin) < grad(idegIdx2)=idegMin set
          ! iidxMin to the new node with smaller degree.
          if (idegIdx .lt. idegMin) then
            iidxMin=idegIdx2
            idegMin=idegIdx
          endif

        end do

        ! At this point, iidxMin contains the index in Ideg of that node with
        ! minimal degree, which has not yet been processed.
        ! Ideg(iidxMin) contains the column number of that column of the matrix,
        ! which corresponds with the node with the next higher degree.

        ! This value is written to Icon. If this subroutine is finished,
        ! Icon therefore contains the column numbers of the Icol array -
        ! sorted for increasing degree of the nodes in the graph of the matrix.
        Icon(Ild(ieq)+idegIdx1) = Ideg(iidxMin)

        ! Set the Ideg-value of the node with minimum degree to 0 to indicate
        ! that it has been precessed. This node is then ignored on the next
        ! searching process for nodes with minimum degree.
        Ideg(iidxMin) = 0
  
      end do
      
      ! Clear auxiliary vector; only some entries were used. This is only for
      ! reasons of safetyness, as if the upper loops are processed correctly,
      ! (no nodes were forgotten), all Ideg-arrays should already be 0.
      do idegIdx=1, Ild(ieq+1)-Ild(ieq)
        Ideg(idegIdx) = 0
      end do
      
    ! Process next line:  
    end do 
  
  end subroutine 

  ! ***************************************************************************

!<subroutine>
  subroutine sstrat_calcColNumberingCM9 (Ild, Icol, Idiag, Icon, Ideg, &
                                         neq, ndeg)
  
    !<description>
    ! Purpose: Cuthill McKee matrix renumbering
    !          Calculate column numbering
  
    ! The algorithm of Cuthill McKee interprets the system matrix as
    ! adjacense matrix. Every Row denotes a node in the corresponing graph.
    ! In the first step this function sorts the columns in every row
    ! for increasing degree of the nodes in the matrix.
    ! The matrix must have symmetric structure!
    ! (For FE matrices this is always the case...)
    !
    ! Matrix storage technique 9 version.
    !</description>
  
  !<input>
  
    ! Number of equations
    integer,intent(IN)                    :: neq

    ! Maximum number of entries != 0 in every row of the matrix
    integer, intent(IN)                   :: ndeg    
   
    ! Row description of matrix
    integer, dimension(neq+1), intent(IN) :: Ild
  
    ! Column description of matrix
    integer, dimension(:), intent(IN)     :: Icol
    
    ! Incides of diagonal elements in structure 9 matrix
    integer, dimension(:), intent(IN)     :: Idiag

  !</input>

  !<inputoutput>
    ! Auxiliary vector; must be at least as long as the
    ! maximum number of entries != 0 in every row of the matrix
    integer, dimension(ndeg), intent(INOUT) :: Ideg

    ! Auxiliary vector; the column numbers of KCOL are assigned to this in
    ! the order of increasing degree. When calling the routine the user
    ! must copy the content of KCOL to this! These values are then
    ! resorted.
    integer, dimension(:), intent(INOUT)  :: Icon
  !</inputoutput>

!</subroutine>

    ! local variables
    integer :: ieq, idegIdx, ildIdx, idegIdx1, idegIdx2
    integer :: idegMin, iidxMin
  
    ! Clear auxiliary vector
    do idegIdx=1, ndeg
      Ideg(idegIdx) = 0
    end do
    
    ! Later in the algorithm the column numbers of the non-diagonal
    ! entries are changed. The numbers of the diagonal entries are not
    ! touched. Therefore we copy the column numbers of the diagonal
    ! entries directly.
    do ieq=1, neq
      Icon(Idiag(ieq)) = Icol(Idiag(ieq))
    end do
    
    ! Loop about all rows in the matrix: ieq = current row.
    ! Every row corresponds to an entry in the solution vector = a node
    ! in the graph of the matrix.
    do ieq=1, neq
    
      ! Copy the column numbers of the non-diagonal entries in the current
      ! column to the auxiliary vector Ideg. Ideg contains always the
      ! following nodes of the current node ieq, which are not processed yet.
      ! The entries of the vector are set to 0 one after the other in the
      ! order of the degree of the following nodes.
      do ildIdx=Ild(ieq), Ild(ieq+1)-1
        Ideg(ildIdx-Ild(ieq)+1) = Icol(ildIdx)
      end do

      ! Set Ideg of the diagonal entry to 0 to prevent it from being 
      ! processed later.
      Ideg(Idiag(ieq)-Ild(ieq)+1) = 0
      
      ! Loop about every column in the current row. The entries in the
      ! row (=column numbers) of the matrix represent the node numbers of
      ! those nodes in the graph of the matrix, which are connected to the
      ! current node ieq.
      
      ! The algorithm now has to sort these nodes for increasing degree.
      ! For this purpose the node with the smallest degree is searched and
      ! deleted from Ideg, then the node with the 2nd smallest degree
      ! and so on.
      ! We don't have many nodes, so we use a simple O(n^2) algorithm here
      ! which simply loops over all nodes.
      do idegIdx1=1, Ild(ieq+1)-Ild(ieq)
      
        ! It's only necessary to search in the non-diagonal entries,
        ! because we want to sort the adjacent nodes of the current one
        ! for increasing degree.
        if (Icol(Ild(ieq)+idegIdx1-1) .ne. ieq) then
        
          ! iidxMin will receive the index in the Ideg-array of the node with
          ! the smallest degree. We start with node 1 in the current column:
          iidxMin=1

          ! The variable idegMin always contains the degree of the node, which
          ! is described by the index iidxMin

          ! If Ideg(iidxMin)=0, this node has already been processed. In this
          ! case idegMin is set to neq - the largest upper bound which implies
          ! that iidxMin is later on surely replaced by a node with a smaller
          ! degree.

          ! Otherwise idegMin receives the degree of the node described by
          ! iidxMin, which is calculated by the number of elements in Icol,
          ! i.e. the difference of the indices in the Ild array.
          if (Ideg(iidxMin) .eq. 0) then
            idegMin = neq
          else
            idegMin = Ild(Ideg(iidxMin)+1)-Ild(Ideg(iidxMin))-1
          endif

          ! Compare iidxMin with every node in that line to find that with
          ! minimum degree:
          do idegIdx2=1, Ild(ieq+1)-Ild(ieq)

            ! If Ideg(idegIdx2)=0, idegIdx2 has already been processed (or it's
            ! the diagonal element which does not have to be processed). Set
            ! idegIdx=neq ; here a lower bound to prevent the current node
            ! iidxMin from being replaced by idegIdx2.
            
            ! Otherwise set idegIdx to the degree of the node with index
            ! idegIdx2
            if (Ideg(idegIdx2) .eq. 0) then
              idegIdx = neq
            else
              idegIdx = Ild(Ideg(idegIdx2)+1)-Ild(Ideg(idegIdx2))-1
            endif
                  
            ! If now idegIdx=grad(iidxMin) < grad(idegIdx2)=idegMin set
            ! iidxMin to the new node with smaller degree.
            if (idegIdx .lt. idegMin) then
              iidxMin=idegIdx2
              idegMin=idegIdx
            endif

          end do

          ! At this point, iidxMin contains the index in Ideg of that node with
          ! minimal degree, which has not yet been processed.
          ! Ideg(iidxMin) contains the column number of that column of the matrix,
          ! which corresponds with the node with the next higher degree.

          ! This value is written to Icon. If this subroutine is finished,
          ! Icon therefore contains the column numbers of the Icol array -
          ! sorted for increasing degree of the nodes in the graph of the matrix.
          Icon(Ild(ieq)+idegIdx1) = Ideg(iidxMin)

          ! Set the Ideg-value of the node with minimum degree to 0 to indicate
          ! that it has been precessed. This node is then ignored on the next
          ! searching process for nodes with minimum degree.
          Ideg(iidxMin) = 0
    
        end if
      
      end do
        
      ! Clear auxiliary vector; only some entries were used. This is only for
      ! reasons of safetyness, as if the upper loops are processed correctly,
      ! (no nodes were forgotten), all Ideg-arrays should already be 0.
      do idegIdx=1, Ild(ieq+1)-Ild(ieq)
        Ideg(idegIdx) = 0
      end do
      
    ! Process next line:  
    end do 
  
  end subroutine 

  ! ***************************************************************************

!<subroutine>
  subroutine sstrat_calcPermutationCM (Ild, Icon, neq, Itr1, Itr2)
  
    !<description>
    ! Purpose: Cuthill McKee matrix renumbering
    !          Calculate permutation

    ! Uses the Ild/Icon vectors of sstrat_calcColumnNumbering to calc-
    ! ulate permutation vectors Itr1 and Itr2; these describe how a
    ! vector must be sorted/sorted back.

    ! Itr1 must be initialized with 0 on call of this subroutine!
    ! Itr2 should be initialized with 0 on call of this subroutine,
    ! as only the entries != 0 are calculated!

    ! Matrix storage technique 7 version!
    !</description>
    
    !<input>
    
    ! Number of equations
    integer, intent(IN)                   :: neq
    
    ! Row description of matrix
    integer, dimension(neq+1), intent(IN) :: Ild
    
    ! Auxiliary vector, calculated with cmsort_calcColumnNumbering
    integer, dimension(:), intent(IN)     :: Icon
    !</input>
    
    !<output>
    
    ! The permutation vector that describes how the solution vector
    ! has to be restored
    integer, dimension(neq), intent(OUT) :: Itr1
    
    ! Describes how a resorted vector can be sorted back 
    integer, dimension(neq), intent(OUT) :: Itr2

    !</output>    

  !</subroutine>
    
    ! local variables
    integer :: ineqIdx, icount, ildIdx, icolIdx
    
!    DO icount=1, neq
!      IF (Itr1(icount) .EQ. 0) Itr1(icount) = icount
!    END DO
    ! Commented out, old FEAT-Code. The current code is using the fact that
    ! the Itr1-vector is filled with 0 on call to this routine.

    ! Itr1 should receive a permutation of the numbers 1..neq for
    ! resorting the entries, Itr2 should receive the inverse permutation.

    ! The Cuthill-McCee algorithm processes levelwise. One node is picked
    ! as the starting node and receives the first number. The nodes
    ! in the neighbourhood receive the next numbers in increasing order
    ! of the degrees of these nodes. Then the neighbourhood of these
    ! nodes is analyzed for increasing degree to give the next node
    ! numbers, and so on.
    !
    ! For this purpose the permutation vectors Itr1 and Itr2 are build
    ! subsequently, depending on which nodes have already been processed.        

    ! We fix the first node (=first line of the matrix), i.e. we set
    ! Itr1(1)=1. This of course implies that the inverse permutation
    ! does not have to process node 1, i.e. Itr2(1)=1. But this will
    ! be computed by the algorithm later.
    
!    Itr1(1)=1
!    Itr2(1)=1
    ! Old FEAT-code, commented out.
    ! The current code uses the fact, that Itr1 is initialized with 0 on
    ! call to this routine. In in the processing a 0 is found in Itr1, a
    ! new matrix block is found. This is especially the case in the
    ! beginning of the algorithm, as no matrix block has been found yet and
    ! Itr1 is completely 0. But also later on if the matrix contains
    ! independent blocks, this will work correctly in contrast to the
    ! old FEAT code.

    ! The renumbering strategy of Cuthill-McKee orders the nodes for
    ! increasing degree in the graph of the matrix. For the new numbering
    ! we use the index variable ineqIdx, which is increased for every node
    ! and specifies always the new node number for the next found node.
    ! We start with node number ineqIdx=0 and increase ineqIdx before
    ! the node number of the next found node is set.
    ineqIdx = 0

    ! Now loop over all nodes (=lines in the matrix). The following loop
    ! processes every node in a kind of quere exactly one time, so it
    ! has neq passes:
    do icount=1, neq

      ! Look at line Itr1(icount) (=0, if a new matrix block starts, like
      ! in the beginning). We look at every element in that line and save
      ! the numbers of the following nodes in Itr1. These are processed later.

      ! If Itr1(icount)=0, a new matrix block starts. This is the case if
      ! a new independent block arises in the matrix without connection
      ! to previous elements. This is normally not the case with finite
      ! elements, but may arise in special cases (e.g. compressed matrices,
      ! where in a line there is only a diagonal entry = 1x1-block).
      ! In this case we ignore the number of the first line of that block
      ! (i.e. we give the next free number directly), and continue with
      ! processing the neighbours of the corresponding node in the graph
      ! of the matrix.
      if (Itr1(icount) .eq. 0) then

        ! New block. all previous blocks are contained in line 1..icount-1.
        ! The next line number to fix is therefore icount.          
        Itr1(icount) = icount

        ! Now the block structure of the matrix (representing the connected 
        ! components in the graph of the matrix) and the structure of the
        ! loop implies:
        !    Itr1(icount) = icount = ineqIdx + 1   !!!
        ! because: If we processed icount-1 elements and start a new block,
        ! the nodes contained in the previous blocks don't contain neighbours
        ! in the current block or in later blocks. Therefore all indices
        ! are <= icount-1, and at the same time ineqIdx=icount-1 holds
        ! (is not set but implicitly the case).

        ! The inverse of the entry is calculated in the DO-loop below,
        ! which also increments ineqIdx appropriately.

      end if

      ! Now loop about the elements in the line and collect the neighbours
      ! of the corresponding node in the graph:        
      do ildIdx=Ild(Itr1(icount)), Ild(Itr1(icount)+1)-1

        ! Collect the column numbers in the current line in the order of
        ! increasing degree of the nodes. icolIdx will receive the column
        ! numbers subsequently.
        icolIdx=Icon(ildIdx)

        ! Test if we already calculated the permutation for the current node.
        ! For that purpose analyze if the Itr2-entry of the currently
        ! processed neighbour contains a value:
        if (Itr2(icolIdx) .eq. 0) then

          ! This is a new node, which follows our current node Itr1(icount).
          ! Give it the next free number.
    
          ! Obviously for the very first element there is:
          !   Itr1(1)=1, icolIdx=Icon(1)=1, Itr2(1)=0, ineqIdx=0.
          ! That way this routine correctly sets Itr2(1)=1 and correctly
          ! calculates the inverse of the permutation.
          ineqIdx       = ineqIdx+1
          Itr2(icolIdx) = ineqIdx

          ! Now insert the new node in Itr1. That way Itr1 also serves as a kind
          ! of queue of the nodes that have not yet been processed - subsequently
          ! Itr1(1), Itr1(2), Itr1(3),... is build up. Because in all cases
          ! icount <= ineqIdx holds, not calculated elements in Itr1 are never
          ! accessed.
          Itr1(ineqIdx) = icolIdx

        end if  

      end do
          
    end do

  end subroutine     

  ! ***************************************************************************

!<subroutine>
  subroutine sstrat_calcXYZsorting (rdiscretisation,Ipermutation,idirection)
  
  !<description>
    ! Computes a column renumbering strategy based on the coordinates of
    ! DOF's. The routine supports 2D and 3D triangulations; in a 2D
    ! triangulation, the Z-coordinate is ignored.
    !
    ! idirection specifies the sorting direction.
    ! The standard sorting direction is idirection=0. In this case, 
    ! the DOF's are sorted rowwise, i.e. first for the X-coordinate,
    ! then for the Y-coordinate. This sorting strategy coincides with the 
    ! FEAST sorting strategy in simple situations like a unit square.
    !
    ! This sorting strategy can only applied for special type discretisations
    ! where the DOF's of the finite elements coincides with some sort of
    ! point coordinates in the domain (like $Q_1$, edge midpoint based
    ! $\tilde Q_1$). If this is not the case, the routine will stop the 
    ! program.
    !
    ! The algorithm acceps a scalar discretisation structure rdiscretisation and
    ! uses its structure to calculate the renumbering. The result
    ! Ipermutation then receives the permutation and its inverse.
  !</description>
    
  !<input>
    ! Spatial discretisation structure that specifies the DOF's and the
    ! triangulation
    type(t_spatialDiscretisation), intent(IN) :: rdiscretisation
    
    ! OPTIONAL: Specifies the sorting direction.
    ! =0: Sort first for X-, then for Y-, then for Z-coordinate.
    !     In 2D this is rowwise sorting.
    ! =1: Sort first for Z-, then for Y-, then for X-coordinate.
    !     In 2D this is columnwise sorting.
    ! If not specified, idirection=0 is assumed.
    integer, intent(IN),optional :: idirection
  !</input>
    
  !<output>
    ! The permutation vector for sorting and its inverse.
    ! With NEQ=NEQ(matrix):
    !   Ipermutation(1:NEQ)       = permutation,
    !   Ipermutation(NEQ+1:2*NEQ) = inverse permutation.
    integer, dimension(:), intent(OUT) :: Ipermutation
  !</output>    

  !</subroutine>
  
    ! local variables
    real(DP), dimension(:,:), pointer :: p_Dcoords,p_Dcoords2
    integer, dimension(2) :: Isize
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer :: hhandle
    integer :: idir
    integer :: ivt,nvt
    integer :: imt,nmt
    integer :: iel,nel
    integer :: idim,ivtlocal
    
    if (rdiscretisation%ndimension .eq. 0) then
      call output_line ('Discretisation not initialised.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'sstrat_calcRowwise')
      call sys_halt()
    end if
    
    idir = 0
    if (present(idirection)) idir=idirection

    ! Depending on the discrisation, choose the correct implementation.
    select case (rdiscretisation%ccomplexity)
    case (SPDISC_UNIFORM)
    
      ! FE-space?
      select case (elem_getPrimaryElement(&
                       rdiscretisation%RelementDistr(1)%celement))
      case (EL_Q1,EL_P1)
      
        ! $Q_1$-element. Take the vertex coordinates as DOF's and sort for that.
        call storage_getbase_double2d (&
            rdiscretisation%p_rtriangulation%h_DvertexCoords,p_Dcoords)
            
        call sortCoords (p_Dcoords, &
          Ipermutation(1:rdiscretisation%p_rtriangulation%NVT), idir)
      
      case (EL_Q2)
      
        ! $Q_2$-element. Allocate an array for all the coordinates
        nvt = rdiscretisation%p_rtriangulation%NVT
        nmt = rdiscretisation%p_rtriangulation%NMT
        nel = rdiscretisation%p_rtriangulation%NEL
        
        Isize(1) = rdiscretisation%p_rtriangulation%ndim
        Isize(2) = nvt+nmt+nel
        call storage_new2D ('rowwiseSorting', 'Dmidpoints', Isize, ST_DOUBLE, &
                            hhandle, ST_NEWBLOCK_NOINIT)
        call storage_getbase_double2d (hhandle,p_Dcoords)
        
        ! Get triangulation information
        call storage_getbase_double2d (&
          rdiscretisation%p_rtriangulation%h_DvertexCoords,p_Dcoords2)
        call storage_getbase_int2d ( &
          rdiscretisation%p_rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
        call storage_getbase_int2d ( &
          rdiscretisation%p_rtriangulation%h_IverticesAtEdge,p_IverticesAtEdge)
        
        ! Copy the vertex coordinates
        do ivt = 1,nvt
          do idim = 1,ubound(p_Dcoords,1)
            p_Dcoords(idim,ivt) = p_Dcoords2(idim,ivt)
          end do
        end do
        
        ! Calculate edge midpoint coordinates
        do imt = 1,nmt
          p_Dcoords(:,nvt+imt) = &
              0.5_DP*p_Dcoords2(:,p_IverticesAtEdge(1,imt)) + &
              0.5_DP*p_Dcoords2(:,p_IverticesAtEdge(2,imt)) 
        end do
        
        ! Calculate element midpoint coordinates
        do iel = 1,nel
          p_Dcoords(:,nvt+nmt+iel) = 0.0_DP
          do ivtlocal = 1,ubound(p_IverticesAtElement,1)
            ivt = p_IverticesAtElement (ivtlocal,iel)
            p_Dcoords(:,nvt+nmt+iel) = &
                p_Dcoords(:,nvt+nmt+iel) &
                + 0.25_DP*p_Dcoords2(:,ivt) 
          end do
        end do
        
        ! Sort for the coordinates
        call sortCoords (p_Dcoords, &
            Ipermutation(1:rdiscretisation%p_rtriangulation%NVT), idir)

        ! Release temp memory
        call storage_free (hhandle)

      case (EL_Q1T)
      
        ! $\tilde Q_1$-element. Take the edge midpoint coordinates as DOF's
        ! and sort for that. We have to calculate the midpoints for that...
        Isize(1) = rdiscretisation%p_rtriangulation%ndim
        Isize(2) = rdiscretisation%p_rtriangulation%NMT
        call storage_new2D ('rowwiseSorting', 'Dmidpoints', Isize, ST_DOUBLE, &
                            hhandle, ST_NEWBLOCK_NOINIT)
        call storage_getbase_double2d (hhandle,p_Dcoords)
        
        ! Call tria_getPointsOnEdge with npointsPerEdge=1; this calculates
        ! the midpoint coordinates.
        call tria_getPointsOnEdge (rdiscretisation%p_rtriangulation,p_Dcoords,1)
        
        ! Sort for the midpoint coordinates
        call sortCoords (p_Dcoords, &
            Ipermutation(1:rdiscretisation%p_rtriangulation%NVT), idir)
        
        ! Release temp array, finish
        call storage_free (hhandle)
      
      case DEFAULT
        call output_line ('Element type not supported.', &
                          OU_CLASS_ERROR,OU_MODE_STD,'sstrat_calcRowwise')
        call sys_halt()
      end select
    
    case DEFAULT
      call output_line ('Discretisation too complex.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'sstrat_calcRowwise')
      call sys_halt()
    end select
  
    ! Calculate the inverse permutation, that's it.
    call sstrat_calcInversePermutation (&
        Ipermutation(1:rdiscretisation%p_rtriangulation%NVT), &
        Ipermutation(rdiscretisation%p_rtriangulation%NVT+1:) )
    
  contains
  
    subroutine sortCoords (Dcoords, Ipermutation,idirection)

    ! Calculates the rowwise sorting. Dcoords must contain a 2D or 3D
    ! array of point coordinates. In Ipermutation, the routine returns
    ! the permutation of these points.
    ! idirection specifies the sorting direction.
    ! =0: points are first ordered for X-, then for Y- and at 
    !     the end for the Z-coordinate (in 3D).
    ! =1: points are first ordered for Z- (in 3D), then for Y- and at 
    !     the end for the X-coordinate.
    
    real(DP), dimension(:,:), intent(IN) :: Dcoords
    integer, dimension(:), intent(OUT) :: Ipermutation
    integer, intent(IN) :: idirection
    
      ! local variables
      integer, dimension(2) :: Isize
      integer :: h_Dsort
      integer :: i
      real(DP), dimension(:,:), pointer :: p_Dsort
      
      ! Allocate a 2D array with (dim(Dcoords)+1,#coords) elements.
      Isize(1) = ubound(Dcoords,1)+1
      Isize(2) = ubound(Dcoords,2)
      call storage_new2D ('rowwiseSorting', 'Dsort', Isize, ST_DOUBLE, &
                          h_Dsort, ST_NEWBLOCK_NOINIT)
      call storage_getbase_double2d (h_Dsort,p_Dsort)
      
      ! In the first element of each ndim+1-tupel, store the number of
      ! the point. In the 2nd/3rd,... element, store the coordinate.
      do i=1,Isize(2)
        p_Dsort(1,i) = real(i,DP)
        p_Dsort(2:,i) = Dcoords(:,i)
      end do
      
      ! Sort the array. First for the last coordinate, then for the
      ! last but one, etc.
      ! Use a stable sorting algorithm to prevent the previous sorting
      ! from getting destroyed.
      select case (idirection)
      case (0)
        do i=1,ubound(Dcoords,1)
          call arraySort_sortByIndex_dp(p_Dsort,1+i,SORT_STABLE)
        end do
        
      case (1)
        do i=ubound(Dcoords,1),1,-1
          call arraySort_sortByIndex_dp(p_Dsort,1+i,SORT_STABLE)
        end do
      
      case DEFAULT
        call output_line ('Invalid direction!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'sstrat_calcXYZsorting')
        call sys_halt()
      end select
      
      ! The first element in each ndim+2-tupel is now the permutation.
      ! Do a type conversion to int to get it.
      do i=1,Isize(2)
        Ipermutation(i) = p_Dsort(1,i)
      end do
      
      ! Release the temp array, that's it.
      call storage_free (h_Dsort)
      
    end subroutine

  end subroutine     

  ! ***************************************************************************

!<subroutine>
  subroutine sstrat_calcFEASTsorting (rdiscretisation,Ipermutation,ifirstVertex)
  
  !<description>
    ! Computes a column renumbering strategy based on the tensor product
    ! structure of the domain. The DOF's are numbered rowwise.
    ! ifirstVertex must specify the first vertex (in the lower left corner)
    ! of the domain. From here, a geometrical search is started to find
    ! all DOF's.
    !
    ! The renumbering strategy is only applicable to special-type 
    ! discretisations (like $Q_1$) and tensor product meshes (FEAST macros).
    ! If this is not the case, the routine will stop the program.
    !
    ! The algorithm acceps a scalar discretisation structure rdiscretisation and
    ! uses its structure to calculate the renumbering. The result
    ! Ipermutation then receives the permutation and its inverse.
  !</description>
    
  !<input>
    ! Spatial discretisation structure that specifies the DOF's and the
    ! triangulation
    type(t_spatialDiscretisation), intent(IN) :: rdiscretisation
    
    ! OPTIONAL: First vertex (lower left corner) of the domain.
    ! If not specified, ifirstVertex=1 is assumed.
    integer, intent(IN), optional :: ifirstVertex
  !</input>
    
  !<output>
    ! The permutation vector for sorting and its inverse.
    ! With NEQ=NEQ(matrix):
    !   Ipermutation(1:NEQ)       = permutation,
    !   Ipermutation(NEQ+1:2*NEQ) = inverse permutation.
    integer, dimension(:), intent(OUT) :: Ipermutation
  !</output>    

  !</subroutine>
  
    ! local variables
    integer, dimension(:,:), pointer :: p_IelementsAtEdge
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer, dimension(:), pointer :: p_IelementsAtVertexIdx
    integer, dimension(:), pointer :: p_IelementsAtVertex
    integer :: iel
    integer :: ivt
    
    if (rdiscretisation%ndimension .eq. 0) then
      call output_line ('Discretisation not initialised.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'sstrat_calcFEASTsorting')
      call sys_halt()
    end if
    
    ! Depending on the discrisation, choose the correct implementation.
    select case (rdiscretisation%ccomplexity)
    case (SPDISC_UNIFORM)
    
      ! FE-space?
      select case (elem_getPrimaryElement(&
                       rdiscretisation%RelementDistr(1)%celement))
      case (EL_Q1)
      
        ! Get geometrical data
        call storage_getbase_int2d (rdiscretisation%p_rtriangulation%h_IelementsAtEdge,&
            p_IelementsAtEdge)
        call storage_getbase_int2d (&
            rdiscretisation%p_rtriangulation%h_IverticesAtElement,&
            p_IverticesAtElement)
        call storage_getbase_int2d (&
            rdiscretisation%p_rtriangulation%h_IedgesAtElement,&
            p_IedgesAtElement)
            
        ! Get the first element on vertex ifirstVertex
        ivt = 1
        if (present(ifirstVertex)) ivt = ifirstVertex
        call storage_getbase_int (rdiscretisation%p_rtriangulation%h_IelementsAtVertex,&
            p_IelementsAtVertex)
        call storage_getbase_int (&
            rdiscretisation%p_rtriangulation%h_IelementsAtVertexIdx,&
            p_IelementsAtVertexIdx)
        iel = p_IelementsAtVertex(p_IelementsAtVertexIdx(ivt))
      
        call sortForFeastQ1 (p_IelementsAtEdge, p_IverticesAtElement, p_IedgesAtElement,&
            Ipermutation(1: rdiscretisation%p_rtriangulation%NVT), ivt, iel, &
            rdiscretisation%p_rtriangulation%NVT)
      
      case DEFAULT
        call output_line ('Element type not supported.', &
                          OU_CLASS_ERROR,OU_MODE_STD,'sstrat_calcFEASTsorting')
        call sys_halt()
      end select
    
    case DEFAULT
      call output_line ('Discretisation too complex.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'sstrat_calcFEASTsorting')
      call sys_halt()
    end select
  
    ! Calculate the inverse permutation, that's it.
    call sstrat_calcInversePermutation (Ipermutation(1:rdiscretisation%p_rtriangulation%NVT), &
                                        Ipermutation(rdiscretisation%p_rtriangulation%NVT+1:) )
    
  contains
  
    subroutine sortForFeastQ1 (IelementsAtEdge,IverticesAtElement,IedgesAtElement,&
                               Ipermutation,ivt,iel,NVT)

    ! Calculates the rowwise sorting in a FEAST like style. 
    ! IelementsAtEdge, IverticesAtElement and IedgesAtElement specify the 
    ! adjacencies in the macro.
    ! In Ipermutation, the routine returns the permutation of these points.
    ! ivt specifies the lower left corner of the macro. iel specifies
    ! the element in the lower left corner that contains ivt.
    
    integer, dimension(:,:), intent(IN) :: IelementsAtEdge
    integer, dimension(:,:), intent(IN) :: IverticesAtElement
    integer, dimension(:,:), intent(IN) :: IedgesAtElement
    integer, dimension(:), intent(OUT) :: Ipermutation
    integer, intent(IN) :: ivt
    integer, intent(IN) :: iel
    integer, intent(IN) :: NVT

      ! local variables
      integer :: icornervertex,icornerelement,ipermidx
      integer :: icurrentvertex,icurrentelement
      integer :: iedgeright,iedgetop
      integer :: ilocalvertex
      integer, parameter :: NVE = 4
      
      ! Current position in the mesh
      icurrentvertex = ivt
      icurrentelement = iel
      
      ! ipermidx counts how many vertices we found.
      ipermidx = 0
      
      ! Loop through the columns until we find the macro border at the top
      do
        
        ! icornervertex remembers the current lower left corner. icornerelement
        ! the corresponding element.
        icornervertex = icurrentvertex
        icornerelement = icurrentelement
        
        ! We are here:
        !
        !    |   |
        !    +---+---
        !    |IEL|
        !  IVT---+---
        !
        ! Add the first vertex to the permutation
        ipermidx = ipermidx+1
        Ipermutation(ipermidx) = icurrentvertex
        
        ! Get the local number of the vertex in the element
        do ilocalvertex = 1,NVE
          if (IverticesAtElement(ilocalvertex,icurrentelement) .eq. icurrentvertex) exit
        end do
        
        ! Get the edges to the neighbour elements
        !
        !    |   |     
        !    +---+---
        !    |   X iedgeright
        !  IVT---+---
        
        iedgeright = IedgesAtElement(mod(ilocalvertex,NVE)+1,icurrentElement)
      
        ! Loop through the macro row until we find the right border of the macro
        do while (IelementsAtEdge(2,iedgeright) .ne. 0) 
        
          ! Step right to the next element 
          !
          !    |   |     
          !    +---+---+
          !    |   |   |           
          !    +--IVT--+
          
          icurrentvertex = IverticesAtElement(mod(ilocalvertex,NVE)+1,icurrentElement)
          
          if (IelementsAtEdge(2,iedgeright) .ne. icurrentElement) then
            icurrentElement = IelementsAtEdge(2,iedgeright)
          else
            icurrentElement = IelementsAtEdge(1,iedgeright)
          end if
          
          ! Add the vertex to the permutation
          ipermidx = ipermidx+1
          Ipermutation(ipermidx) = icurrentvertex
          
          ! Get the local number of the vertex in the element
          do ilocalvertex = 1,NVE
            if (IverticesAtElement(ilocalvertex,icurrentelement) &
                .eq. icurrentvertex) exit
          end do
        
          ! Get the edges to the neighbour elements
          !
          !    |   |   |  
          !    +---+---+--
          !    |   |   X iedgeright
          !    +--IVT--+--
          
          iedgeright = IedgesAtElement(mod(ilocalvertex,NVE)+1,icurrentElement)
          
        end do
        
        ! We have reached the end of the row
        !
        !        |   |   |
        !     ---+---+---+
        !        |   |   X iedgeright
        !     ---+--IVT--IVT2
        !
        ! Remember the last vertex IVT2 in the row
        
        icurrentvertex = IverticesAtElement(mod(ilocalvertex,NVE)+1,icurrentelement)
        ipermidx = ipermidx+1
        Ipermutation(ipermidx) = icurrentvertex
        
        ! Hop back to the first element and find the local number of the first
        ! vertex in the row again.
        !
        !    |   |
        !    +---+---
        !    |IEL|
        !  IVT---+---

        icurrentvertex = icornervertex
        icurrentelement = icornerelement
        
        do ilocalvertex = 1,NVE
          if (IverticesAtElement(ilocalvertex,icurrentelement) .eq. icurrentvertex) exit
        end do
        
        ! Get the edge that leads to the element row above us.
        !
        !    | iedgetop
        !    +-X-+---
        !    |IEL|
        !  IVT---+---
        
        iedgetop = IedgesAtElement(mod(ilocalvertex+1,NVE)+1,icurrentElement)
        
        ! Get the vertex and the element there. Note: If there is no neighbour
        ! element, the current element number gets =0 which is the terminal criterion.
        !
        !    |IEL|     
        !  IVT---+---
        !    |   |
        !    +---+---
      
        icurrentvertex = IverticesAtElement(mod(ilocalvertex+2,NVE)+1,icurrentElement)
      
        if (IelementsAtEdge(2,iedgetop) .ne. icurrentElement) then
          icurrentElement = IelementsAtEdge(2,iedgetop)
        else
          icurrentElement = IelementsAtEdge(1,iedgetop)
        end if
      
        if (icurrentelement .eq. 0) then
          exit
        else
          ! There is a neighbour element on top.
          ! Get the edge that leads to the right and continue in the new
          ! element row.
          do ilocalvertex = 1,NVE
            if (IverticesAtElement(ilocalvertex,icurrentelement) .eq. icurrentvertex) exit
          end do

          iedgeright = IedgesAtElement(mod(ilocalvertex,NVE)+1,icurrentElement)
        end if
      end do
      
      ! Ok, we have done all element rows now.
      !
      !   |   |   |   |   |
      !   +---+---+---+---+
      !   |   |   |   |   |           
      !   +---+---+---+---+
      !
      ! What's still missing is the final edge-row on top! This is another loop.
      ! Switch the element number back to the last remembered one, then we have: 
      !
      !  IVT--+---+---+---+
      !   |IEL|   |   |   |           
      !   +---+---+---+---+
      !   |   |   |   |   |
      
      icurrentelement = icornerelement
      
      ! Add the vertex to the permutation
      ipermidx = ipermidx+1
      Ipermutation(ipermidx) = icurrentvertex
      
      ! Local number and edge to the right? 
      !
      !  IVT--+--           
      !   |   X iedgeright             
      !   +---+-            
      !   |   |            

      do ilocalvertex = 1,NVE
        if (IverticesAtElement(ilocalvertex,icurrentelement) .eq. icurrentvertex) exit
      end do
      
      iedgeright = IedgesAtElement(mod(ilocalvertex+1,NVE)+1,icurrentElement)
    
      ! Loop through the macro row until we find the right border of the macro
      do while (IelementsAtEdge(2,iedgeright) .ne. 0) 
      
        ! Step right to the next element 
        !
        !    +--IVT--+
        !    |   |IEL|           
        !    +---+---+
        !    |   |     
        
        icurrentvertex = IverticesAtElement(mod(ilocalvertex+2,NVE)+1,icurrentElement)
        
        if (IelementsAtEdge(2,iedgeright) .ne. icurrentElement) then
          icurrentElement = IelementsAtEdge(2,iedgeright)
        else
          icurrentElement = IelementsAtEdge(1,iedgeright)
        end if
        
        ! Add the vertex to the permutation
        ipermidx = ipermidx+1
        Ipermutation(ipermidx) = icurrentvertex
        
        ! Get the local number of the vertex in the element
        do ilocalvertex = 1,NVE
          if (IverticesAtElement(ilocalvertex,icurrentelement) &
              .eq. icurrentvertex) exit
        end do
      
        ! Get the edges to the neighbour elements
        !
        !    +--IVT--+--
        !    |   |   X iedgeright
        !    +---+---+--
        !    |   |   |  
        
        iedgeright = IedgesAtElement(mod(ilocalvertex+1,NVE)+1,icurrentElement)
        
      end do    
      
      ! Remember the last vertex IVT2 in the row
      !
      !     --IVT--IVT2
      !        |   |           
      !     ---+---+  
      !        |   |  
      
      icurrentvertex = IverticesAtElement(mod(ilocalvertex+2,NVE)+1,icurrentelement)
      ipermidx = ipermidx+1
      Ipermutation(ipermidx) = icurrentvertex

    end subroutine

  end subroutine     

  ! ***************************************************************************

!<subroutine>

  subroutine sstrat_calcHierarchical (Rdiscretisation,Ipermutation)
  
!<description>
  ! This subroutine calculates a hierarchical renumbering strategy.
  ! Based on a set of discretisation structures defining different levels
  ! of a level hierarchy coming from a refinement process, the permutation
  ! calculated by this routine tries to group DOF's by element macros.
!</description>
    
!<input>
  
  ! Array of discretisation structures identifying the different levels
  ! of refinement. The discretisation structures must stem from a standard
  ! 2-level refinement.
  ! The last element in this array must correspond
  ! to the permutation which is to be computed.
  type(t_spatialDiscretisation), dimension(:), intent(IN) :: Rdiscretisation
  
!</input>
    
!<output>
  ! The permutation vector for sorting and its inverse.
  ! With NEQ=NEQ(matrix):
  !   Ipermutation(1:NEQ)       = permutation,
  !   Ipermutation(NEQ+1:2*NEQ) = inverse permutation.
  integer, dimension(:), intent(OUT) :: Ipermutation
!</output>    

!</subroutine>
  
    integer :: N
    
    ! Length of the permutation. Must correspond to the #DOF's
    ! on the finest level.
    N = size(Ipermutation)/2
    
    if (N .ne. dof_igetNDofGlob(Rdiscretisation(size(Rdiscretisation)))) then
      call output_line ('Permutation target vector has the wrong size!', &
          OU_CLASS_ERROR,OU_MODE_STD,'sstrat_calcHierarchical')        
      call sys_halt()
    end if
  
    select case (Rdiscretisation(1)%p_rtriangulation%ndim)
    case (NDIM2D)
      ! Regular 2-level refinement in 2D.
      call calcHierarch(Rdiscretisation,Ipermutation)
    
    case DEFAULT
      call output_line ('Invalid dimension.', &
          OU_CLASS_ERROR,OU_MODE_STD,'sstrat_calcHierarchical')        
      call sys_halt()
    end select

    ! Calculate the inverse permutation, that's it.
    call sstrat_calcInversePermutation (Ipermutation(1:N), Ipermutation(N+1:) )

  contains
  
    subroutine calcHierarch(Rdiscretisation,Ipermutation)
  
    ! Array of discretisation structures identifying the different levels
    ! of refinement. The discretisation structures must stem from a standard
    ! 2-level refinement.
    ! The last element in this array must correspond
    ! to the permutation which is to be computed.
    type(t_spatialDiscretisation), dimension(:), intent(IN), target :: Rdiscretisation
  
    ! The permutation vector for sorting and its inverse.
    ! With NEQ=NEQ(matrix):
    !   Ipermutation(1:NEQ)       = permutation,
    !   Ipermutation(NEQ+1:2*NEQ) = inverse permutation.
    integer, dimension(:), intent(OUT) :: Ipermutation

      ! local variables      
      integer, dimension(:), pointer :: p_IrefinementPatch
      integer, dimension(:), pointer :: p_IrefinementPatchIdx
      type(t_triangulation), pointer :: p_rtriaCoarse,p_rtria
      integer :: NEQ
      integer :: hmarker
      integer, dimension(:), pointer :: p_Imarker
      integer, dimension(EL_MAXNBAS) :: Idofs
      integer, dimension(size(Rdiscretisation)) :: IpatchIndex
      integer, dimension(size(Rdiscretisation)) :: ImaxIndex
      integer, dimension(size(Rdiscretisation)) :: Ielement
      type(t_levelHirarchy), dimension(size(Rdiscretisation)) :: Rhierarchy
      integer :: ilev,ndof,ieldistr,idof
      integer(I32) :: celement
      integer :: ipos
      integer :: ielcoarse
      logical :: bisUniform
      
      ! Save pointers to the element patch arrays for all levels.
      ! We'll frequently need them.
      do ilev=1,size(Rhierarchy)
        ! Refinement information
        if (ilev .gt. 1) then
          p_rtria => Rdiscretisation(ilev)%p_rtriangulation
          call storage_getbase_int (p_rtria%h_IrefinementPatch,&
              Rhierarchy(ilev)%p_IrefinementPatch)
          call storage_getbase_int (p_rtria%h_IrefinementPatchIdx,&
              Rhierarchy(ilev)%p_IrefinementPatchIdx)
        end if
            
        ! Information about the discretisation: Arrays that allow
        ! to determine the type of an element.
        bisUniform = Rdiscretisation(ilev)%ccomplexity .eq. SPDISC_UNIFORM
        Rhierarchy(ilev)%bisUniform = bisUniform
        
        if (bisUniform) then
          ! One element type for all elements
          Rhierarchy(ilev)%celement = &
              Rdiscretisation(ilev)%RelementDistr(1)%celement
        else
          ! A different element type for every element.
          ! Get a pointer to the array that defines the element distribution
          ! of the element. This allows us later to determine the element type.
          call storage_getbase_int (p_rtria%h_IrefinementPatchIdx,&
              Rhierarchy(ilev)%p_IelementDistr)
        end if
      end do

      p_rtriaCoarse => Rdiscretisation(1)%p_rtriangulation

      ! Get the number of DOF's on the finest level. This is the
      ! size of the permutation.          
      NEQ = dof_igetNDofGlob(Rdiscretisation(size(Rdiscretisation)))

      ! Set up a marker array where we remember whether we processed
      ! a DOF or not. Initialise with zero; all DOF's we already
      ! processed are marked here with a 1.
      call storage_new ('calcHierarch2Level2D', &
          'mark', NEQ, ST_INT, hmarker, ST_NEWBLOCK_ZERO)
      call storage_getbase_int (hmarker,p_Imarker)

      ipos = 0
      ilev = 1

      ! IelementPatch(i) saves an index into the IrefinementPatch
      ! array. When being on element iel on the coarser mesh i-1,
      ! IelementPatch(i) is a pointer to the current subelement
      ! in the finer mesh on level i.
      !
      ! Ielement(i) on the other hand saves the current element number
      ! on level i.
      !
      ! Loop through all elements on the coarse mesh
      do ielcoarse = 1,p_rtriaCoarse%NEL
      
        Ielement(1) = ielcoarse
      
        patchcycle: do
      
          ! From this element, figure out the DOF's of all subelements.
          ! This has to be done patchwise on the finest level.
          ! Thus, as long as we aren't on the fines level, we have to
          ! increase the current one.
          do while (ilev .lt. size(Rdiscretisation))
          
            ! Go up
            ilev = ilev + 1
            
            ! Ielement(ilev-1) is not the patch number for level ilev.
            ! Set the start index to the first element in that patch
            ! of the finer mesh. Remember the maximum index 'where the patch ends).
            IpatchIndex(ilev) = &
              Rhierarchy(ilev)%p_IrefinementPatchIdx(Ielement(ilev-1))
            ImaxIndex(ilev) = &
              Rhierarchy(ilev)%p_IrefinementPatchIdx(Ielement(ilev-1)+1)-1
              
            ! Get the element number of the first element in the patch.
            Ielement(ilev) = &
              Rhierarchy(ilev)%p_IrefinementPatch(IpatchIndex(ilev))
          
          end do
          
          ! We are now on the maximum level on element Ielement(max).
          ! Get the DOF's of that element.
          ! For that purpose, we need the element type.
          if (Rhierarchy(ilev)%bisUniform) then
            celement = Rhierarchy(ilev)%celement
          else
            ! Get the element distribution and from that the element type.
            ieldistr = Rhierarchy(ilev)%p_IelementDistr(Ielement(ilev))
            celement = Rdiscretisation(ilev)%RelementDistr(ieldistr)%celement
          end if
          
          ndof = elem_igetNDofLoc(celement)
          call dof_locGlobMapping(Rdiscretisation(ilev), Ielement(ilev),  Idofs)
          
          ! Check the DOF's. All DOF's we don't have yet, we collect into the
          ! permutation.
          do idof = 1,ndof
            if (p_Imarker(Idofs(idof)) .eq. 0) then
              ipos = ipos + 1
              Ipermutation(ipos) = Idofs(idof)
              
              ! Mark the DOF as being handled.
              p_Imarker(Idofs(idof)) = 1
            end if
          end do
        
          ! Now we have to proceed to the next element. How to do that depends
          ! on 'where we are'.
          if (ilev .gt. 1) then
          
            ! Go to the next element in the current patch.
            IpatchIndex(ilev) = IpatchIndex(ilev) + 1
          
            do while ((ilev .gt. 1) .and. &
                    (IpatchIndex(ilev) .gt. ImaxIndex(ilev)))
            
              ! All elements of the patch completed. Go down one level
              ! and proceed there to the next element patch.
              ilev = ilev - 1  
              IpatchIndex(ilev) = IpatchIndex(ilev) + 1
            
            end do
          end if
          
          ! As long as we don't reach level 1, there are elements left
          ! in the patch to proceed. So cycle the patchloop
          ! to proceed to the next element.
          if (ilev .eq. 1) then
            exit patchcycle
          else
            ! Get the new current element number
            Ielement(ilev) = &
              Rhierarchy(ilev)%p_IrefinementPatch(IpatchIndex(ilev))
          end if
          
        end do patchcycle
      
      end do
      
      ! Release temp memory.
      call storage_free (hmarker)

      ! Calculate the inverse permutation, that's it.
      call sstrat_calcInversePermutation (Ipermutation(1:NEQ), Ipermutation(NEQ+1:) )

    end subroutine

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine sstrat_calcInversePermutation (IpermutationSource, IpermutationDest)
  
  !<description>
    ! Computes the inverse of a permutation. IpermutationSource is a given
    ! permutation. IpermutationDest will receive the inverse permutation.
  !</description>
    
  !<input>
    ! A permutation.
    integer, dimension(:), intent(IN) :: IpermutationSource
  !</input>
    
  !<output>
    ! An array of the same size as IpermutationSource. Receives the inverse
    ! permutation.
    integer, dimension(:), intent(OUT) :: IpermutationDest
  !</output>    

  !</subroutine>
  
    integer :: i
    do i=1,size(IpermutationSource)
      IpermutationDest(IpermutationSource(i)) = i
    end do
  
  end subroutine

end module
