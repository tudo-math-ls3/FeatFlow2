!##############################################################################
!# ****************************************************************************
!# <name> bilinearformevaluation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines for the discretisation of bilinear forms,
!# i.e. the creation of matrices and matrix structures. It contains the
!# following set of routines:
!#
!# 1.) bilf_createMatrixStructure
!#     -> Creates a 'scalar' matrix structure in a specified matrix
!#        format according to a discretisation.
!#
!# 2.) bilf_buildMatrixScalar
!#     -> Assembles the entries of a matrix, which structure was build
!#        with bilf_createMatrixStructure before.
!#
!# Auxiliary routines
!#
!# 1.) bilf_getLocalMatrixIndices
!#     -> Calculate the positions of local matrices in a global matrix.
!#
!# </purpose>
!#
!# NOTE:
!#   In Windows if OpenMP is activated, this source file must not be processed
!#   with checking for uninitialised variables enabled! The Intel Fortran
!#   compiler usually gets messed up with that!
!##############################################################################

module bilinearformevaluation

  use fsystem
  use genoutput
  use linearsystemscalar
  use spatialdiscretisation
  use scalarpde
  use derivatives
  use cubature
  use domainintegration
  use element
  use elementpreprocessing
  use collection, only: t_collection
  
  implicit none

!<types>
  
!<typeblock>
  
  ! A private type block. Used as memory block during the creation of matrices.
  type t_matrixmem
    ! A pointer to a 1D memory block of integers - receives the
    ! column identifiers of a matrix
    integer, dimension(:), pointer :: p_Icol
    
    ! A pointer to a 1D memory block of integers - receives 
    ! indices of next entries in the list of each line in the matrix
    integer, dimension(:), pointer :: p_Iindx
  end type
  
!</typeblock>

  private :: t_matrixmem

!</types>
  
!<constants>

!<constantblock description="Method identifiers for construction of matrix structure.">

  ! Element-based matrix construction. This is the standard matrix construction method. 
  integer, parameter :: BILF_MATC_ELEMENTBASED = 0
  
  ! Edge-based matrix construction. The matrix stencil is extended in such a way,
  ! that the DOF's of one element may interact with the DOF's of all other elements
  ! that are adjacent via one of the edges.
  integer, parameter :: BILF_MATC_EDGEBASED    = 1

  ! Vertex-based matrix construction. The matrix stencil is extended in such a way,
  ! that the DOF's of one element may interact with the DOF's of all other elements
  ! that are adjacent via one of the corner vertices.
  integer, parameter :: BILF_MATC_VERTEXBASED  = 2

!</constantblock>

!<constantblock description="Constants defining the blocking of the assembly">

  ! Number of elements to handle simultaneously when building matrices
  integer, parameter :: BILF_NELEMSIM   = 1000
  
!</constantblock>
!</constants>

contains

  !****************************************************************************

!<subroutine>

  subroutine bilf_createMatrixStructure (rdiscretisationTrial,iformat,rmatrixScalar, &
                                         rdiscretisationTest,cconstrType,imemguess)
  
!<description>
  ! This routine allows to calculate the structure of a finite-element matrix
  ! on the heap. The size of the matrix is determined dynamically.
!</description>

!<input>
  ! The underlying discretisation structure which is to be used to
  ! create the matrix. Specifies the discretisation of the trial functions.
  ! If rdiscretisationTest is not specified, this also specifies
  ! the discretisation of the test functions, this test and trial
  ! functions coincide.
  type(t_spatialDiscretisation), intent(IN), target :: rdiscretisationTrial
  
  ! Format of the matrix structure to be created. One of the LSYSSC_xxxx
  ! constants.
  integer, intent(IN) :: iformat

  ! OPTIONAL: One of the BILF_MATC_xxxx constants that allow to specify
  ! the matrix construction method. If not specified,
  ! BILF_MATC_ELEMENTBASED is used.
  integer, intent(IN), optional      :: cconstrType

  ! OPTIONAL: The underlying discretisation structure for the test functions.
  ! If not specified, the trial functions coincide with the test functions.
  type(t_spatialDiscretisation), intent(IN), target, optional :: rdiscretisationTest

  ! OPTIONAL: An initial guess about how much memory the matrix needs. If set 
  ! to 0 or not given, an initial guess of 16*NEQ (but at least 10000 matrix 
  ! entries) is assumed.
  integer, intent(IN), optional :: imemGuess
!</input>

!<output>
  ! The structure of a scalar matrix, fitting to the given discretisation.
  ! Memory fo the structure is allocated dynamically on the heap.
  type(t_matrixScalar), intent(OUT) :: rmatrixScalar
!</output>

!</subroutine>

  ! local variables
  integer :: imem
  integer :: ccType
  
  imem = 0
  if (present(imemguess)) then
    imem = max(0,imemguess)
  end if
  
  ccType = BILF_MATC_ELEMENTBASED
  if (present(cconstrType)) ccType = cconstrType
  
  ! Do we have a not too complex triangulation? Would simplify a lot...
  if ( (rdiscretisationTrial%ccomplexity .eq. SPDISC_UNIFORM) .or. &
       (rdiscretisationTrial%ccomplexity .eq. SPDISC_CONFORMAL) ) then
  
    ! Which matrix structure do we have to create?
    select case (iformat) 
    
    case (LSYSSC_MATRIX9)
    
      select case (ccType)
      
      case (BILF_MATC_ELEMENTBASED)
        ! Call the creation routine for structure 9:
        call bilf_createMatStructure9_conf (rdiscretisationTrial,rmatrixScalar,&
            rdiscretisationTest,imem)
        
      case (BILF_MATC_EDGEBASED)
      
        if (present(rdiscretisationTest)) then
      
          if (rdiscretisationTest%ccomplexity .ne. SPDISC_UNIFORM) then
            call output_line ('Edge-based matrix constrution only for'//&
                    ' uniform discr., supported.', &
                    OU_CLASS_ERROR,OU_MODE_STD,'bilf_createMatrixStructure')
            call sys_halt()
          end if
        
        end if

        call bilf_createMatStructure9eb_uni (rdiscretisationTrial,rmatrixScalar,&
            rdiscretisationTest,imem)
        
      case DEFAULT
        call output_line ('Invalid matrix construction method.', &
                OU_CLASS_ERROR,OU_MODE_STD,'bilf_createMatrixStructure')
        call sys_halt()
      end select
      
    case (LSYSSC_MATRIX7)
    
      select case (ccType)
      
      case (BILF_MATC_ELEMENTBASED)
      
        ! Call the creation routine for structure 9:
        call bilf_createMatStructure9_conf (rdiscretisationTrial,rmatrixScalar,&
            rdiscretisationTest,imem)

      case (BILF_MATC_EDGEBASED)
      
        if (present(rdiscretisationTest)) then
          if (rdiscretisationTest%ccomplexity .eq. SPDISC_UNIFORM) then
            call output_line ('Edge-based matrix constrution only for'//&
                    ' uniform discr., supported.', &
                    OU_CLASS_ERROR,OU_MODE_STD,'bilf_createMatrixStructure')
            call sys_halt()
          end if
        end if
        
        call bilf_createMatStructure9eb_uni (rdiscretisationTrial,rmatrixScalar,&
          rdiscretisationTest,imem)
        
      case DEFAULT
        call output_line ('Invalid matrix construction method.', &
            OU_CLASS_ERROR,OU_MODE_STD,'bilf_createMatrixStructure')
        call sys_halt()
      end select
        
      ! Translate to matrix structure 7:
      call lsyssc_convertMatrix (rmatrixScalar,LSYSSC_MATRIX7)
      
    case DEFAULT
      call output_line ('Not supported matrix structure!', &
          OU_CLASS_ERROR,OU_MODE_STD,'bilf_createMatrixStructure')
      call sys_halt()
    end select
  
  else
    call output_line ('General discretisation not implemented!', &
        OU_CLASS_ERROR,OU_MODE_STD,'bilf_createMatrixStructure')
    call sys_halt()
  end if

  end subroutine
  
  !****************************************************************************

!<subroutine>

  subroutine bilf_buildMatrixScalar (rform,bclear,rmatrixScalar,&
                                     fcoeff_buildMatrixSc_sim,rcollection)
  
!<description>
  ! This routine calculates the entries of a finite element matrix.
  ! The matrix structure must be prepared with bilf_createMatrixStructure
  ! in advance.
  ! In case the array for the matrix entries does not exist, the routine
  ! allocates memory in size of the matrix of the heap for the matrix entries.
  !
  ! For setting up the entries, the discretisation structure attached to
  ! the matrix is used (rmatrixScalar%p_rdiscretisation). This is
  ! normally attached to the matrix by bilf_createMatrixStructure.
  !
  ! The matrix must be unsorted when this routine is called, 
  ! otherwise an error is thrown.
!</description>

!<input>
  ! The bilinear form specifying the underlying PDE of the discretisation.
  type(t_bilinearForm), intent(IN) :: rform
  
  ! Whether to clear the matrix before calculating the entries.
  ! If .FALSE., the new matrix entries are added to the existing entries.
  logical, intent(IN) :: bclear
  
  ! OPTIONAL: A collection structure. This structure is given to the
  ! callback function for nonconstant coefficients to provide additional
  ! information. 
  type(t_collection), intent(INOUT), target, optional :: rcollection
  
  ! OPTIONAL: A callback routine for nonconstant coefficient matrices.
  ! Must be present if the matrix has nonconstant coefficients!
  include 'intf_coefficientMatrixSc.inc'
  optional :: fcoeff_buildMatrixSc_sim
!</input>

!<inputoutput>
  ! The FE matrix. Calculated matrix entries are imposed to this matrix.
  type(t_matrixScalar), intent(INOUT) :: rmatrixScalar
!</inputoutput>

!</subroutine>

  ! local variables
  type(t_matrixScalar) :: rmatrixBackup
  
  ! The matrix must be unsorted, otherwise we can't set up the matrix.
  ! Note that we cannot switch off the sorting as easy as in the case
  ! of a vector, since there's a structure behind the matrix! So the caller
  ! has to make sure, the matrix is unsorted when this routine is called.
  if (rmatrixScalar%isortStrategy .gt. 0) then
    call output_line ('Matrix-structure must be unsorted!', &
        OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalar')
    call sys_halt()
  end if

  if ((.not. associated(rmatrixScalar%p_rspatialDiscrTest)) .or. &
      (.not. associated(rmatrixScalar%p_rspatialDiscrTrial))) then
    call output_line ('No discretisation associated!', &
        OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalar')
    call sys_halt()
  end if

  ! Do we have a uniform triangulation? Would simplify a lot...
  select case (rmatrixScalar%p_rspatialDiscrTest%ccomplexity)
  case (SPDISC_UNIFORM) 
    ! Uniform discretisation; only one type of elements, e.g. P1 or Q1
    select case (rmatrixScalar%cdataType)
    case (ST_DOUBLE) 
      ! Which matrix structure do we have?
      select case (rmatrixScalar%cmatrixFormat) 
      case (LSYSSC_MATRIX9)
        !IF (PRESENT(fcoeff_buildMatrixSc_sim)) THEN
          call bilf_buildMatrix9d_conf3 (rform,bclear,rmatrixScalar,&  
                                         fcoeff_buildMatrixSc_sim,rcollection)
        !ELSE
        !  CALL bilf_buildMatrix9d_conf2 (rform,bclear,rmatrixScalar)
        !END IF
      case (LSYSSC_MATRIX7)
        ! Convert structure 7 to structure 9.For that purpose, make a backup of
        ! the original matrix...
        call lsyssc_duplicateMatrix (rmatrixScalar,rmatrixBackup,&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
            
        ! Convert the matrix 
        call lsyssc_convertMatrix (rmatrixBackup,LSYSSC_MATRIX9)
        
        ! Create the matrix in structure 9
        call bilf_buildMatrix9d_conf3 (rform,bclear,rmatrixBackup,&  
                                       fcoeff_buildMatrixSc_sim,rcollection)
                                       
        ! Convert back to structure 7
        call lsyssc_convertMatrix (rmatrixBackup,LSYSSC_MATRIX7)
        
        ! Copy the entries back to the original matrix and release memory.
        call lsyssc_duplicateMatrix (rmatrixBackup,rmatrixScalar,&
            LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPYOVERWRITE)
            
        call lsyssc_releaseMatrix (rmatrixBackup)
                                       
      case DEFAULT
        call output_line ('Not supported matrix structure!', &
            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalar')
        call sys_halt()
      end select
    case DEFAULT
      call output_line ('Single precision matrices currently not supported!', &
          OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalar')
      call sys_halt()
    end select
    
  case (SPDISC_CONFORMAL) 
    
    ! Conformal discretisation; may have mixed P1/Q1 elements e.g.
    select case (rmatrixScalar%cdataType)
    case (ST_DOUBLE) 
      ! Which matrix structure do we have?
      select case (rmatrixScalar%cmatrixFormat) 
      case (LSYSSC_MATRIX9)
        !IF (PRESENT(fcoeff_buildMatrixSc_sim)) THEN
          call bilf_buildMatrix9d_conf3 (rform,bclear,rmatrixScalar,&  
                                         fcoeff_buildMatrixSc_sim,rcollection)
        !ELSE
        !  CALL bilf_buildMatrix9d_conf2 (rform,bclear,rmatrixScalar)
        !END IF
        
      case (LSYSSC_MATRIX7)
        ! Convert structure 7 to structure 9.For that purpose, make a backup of
        ! the original matrix...
        call lsyssc_duplicateMatrix (rmatrixScalar,rmatrixBackup,&
            LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
            
        ! Convert the matrix 
        call lsyssc_convertMatrix (rmatrixBackup,LSYSSC_MATRIX9)
        
        ! Create the matrix in structure 9
        call bilf_buildMatrix9d_conf3 (rform,bclear,rmatrixBackup,&  
                                       fcoeff_buildMatrixSc_sim,rcollection)
                                       
        ! Convert back to structure 7
        call lsyssc_convertMatrix (rmatrixBackup,LSYSSC_MATRIX7)
        
        ! Copy the entries back to the original matrix and release memory.
        call lsyssc_duplicateMatrix (rmatrixBackup,rmatrixScalar,&
            LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPYOVERWRITE)
            
        call lsyssc_releaseMatrix (rmatrixBackup)

      case DEFAULT
        call output_line ('Not supported matrix structure!', &
            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalar')
        call sys_halt()
      end select
    case DEFAULT
      call output_line ('Single precision matrices currently not supported!', &
          OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalar')
      call sys_halt()
    end select
  case DEFAULT
    call output_line ('General discretisation not implemented!', &
        OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalar')
    call sys_halt()
  end select

  end subroutine
  
  !****************************************************************************
  
!<subroutine>
  
  subroutine bilf_createMatStructure9_conf (rdiscretisationTrial,rmatrixScalar,&
      rdiscretisationTest,imemGuess)
  
!<description>
  ! This routine creates according to a given discretisation the matrix 
  ! structure of a structure-9 matrix. The discretisation is assumed to be
  ! conformal, i.e. the DOF's of different FE spaces in the trial space
  ! fit together. The function space for trial and test functions 
  ! may be different.
!</description>

!<input>
  ! The underlying discretisation structure which is to be used to
  ! create the matrix. Specifies the discretisation of the test functions.
  ! If rdiscretisationTrial is not specified, this also specifies
  ! the discretisation of the trial functions, this test and trial
  ! functions coincide.
  type(t_spatialDiscretisation), intent(IN), target :: rdiscretisationTrial
  
  ! OPTIONAL: The underlying discretisation structure for the trial functions.
  ! If not specified, the trial functions coincide with the test functions.
  type(t_spatialDiscretisation), intent(IN), target, optional :: rdiscretisationTest

  ! An initial guess about how much memory the matrix needs. If set to 0,
  ! an initial guess of 16*NEQ (but at least 10000 matrix entries) is assumed.
  integer, intent(IN) :: imemGuess
  
!</input>

!<output>
  ! The structure of a scalar matrix, fitting to the given discretisation.
  ! Memory fo rthe structure is allocated dynamically on the heap.
  type(t_matrixScalar), intent(OUT) :: rmatrixScalar
!</output>

!</subroutine>

  ! local variables
  integer :: NEQ, IEQ, IROW, JCOL, IPOS, istartIdx, NA, nmaxCol
  integer :: IDOFE, JDOFE, i, IHELP,NVE
  integer :: IEL, IELmax, IELset
  logical :: BSORT, bIdenticalTrialAndTest
  
  ! An allocateable list of handles for memory blocks. Size is dynamically 
  ! increased if there are too many columns in the matrix.
  integer, dimension(:), pointer :: p_Ihcol, p_Ihindx, p_IhTmp
  integer, dimension(:), pointer :: p_Isize, p_ISizeTmp
  
  ! An allocateable list of pointers to the memory blocks - corresponds
  ! to the handles in p_Ihcol/p_Ihindx
  type(t_matrixmem), dimension(:), pointer :: Rmemblock, RmemblockTmp
  
  ! Pointer to current KCOL memory block,
  ! pointer to current index memory block
  integer, dimension(:), pointer :: p_Icol, p_Iindx
  
  ! Number of currently allocated pointers in Ihmemblock
  integer :: iblocks
  
  ! Number of currently active element distribution
  integer :: icurrentElementDistr
  
  ! Currently active memory block
  integer :: icurrentblock
  
  ! Number of elements in the current element distribution
  integer :: NEL
  
  ! Size of memory blocks
  integer :: imemblkSize
  
  ! Blocksize in terms of NEQ for guessing memory.
  ! The initial guess for memory is iblkSize*iblkSize*NEQ and every time
  ! memory is needed, another iblkSize*NEQ elements are added.
  integer, parameter :: iblkSize = 4
  
  ! Number of memory blocks to allocate
  integer, parameter :: NmemBlkCount = 5

  ! Pointer to KLD, KCOL, diagonal
  integer, dimension(:), pointer :: p_KLD, p_KCOL, p_Kdiagonal
  
  ! Size of memory currently allocated
  integer :: iallocated
  
  ! An allocateable array accepting the DOF's of a set of elements.
  !INTEGER, DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsTest, IdofsTrial
  !INTEGER, DIMENSION(:,:), POINTER :: p_IdofsTrial
  integer, dimension(EL_MAXNBAS,BILF_NELEMSIM), target :: IdofsTest, IdofsTrial
  integer, dimension(:,:), pointer :: p_IdofsTrial
  
  ! Number of local degees of freedom for trial and test functions
  integer :: indofTrial, indofTest
  
  ! The triangulation structure - to shorten some things...
  type(t_triangulation), pointer :: p_rtriangulation
  
  ! A pointer to an element-number list
  integer, dimension(:), pointer :: p_IelementList
  
  ! Current element distribution
  type(t_elementDistribution), pointer :: p_relementDistrTest
  type(t_elementDistribution), pointer :: p_relementDistrTrial

  ! Number of elements in a block. Normally =BILF_NELEMSIM,
  ! except if there are less elements in the discretisation.
  integer :: nelementsPerBlock

  ! The algorithm is: Test every DOF on one element against each other
  ! DOF on the same element and save the combination into a matrix
  ! in structure 9!
  !
  ! At first, initialise the structure-9 matrix:
  
  rmatrixScalar%p_rspatialDiscrTrial => rdiscretisationTrial
  if (present(rdiscretisationTest)) then
    rmatrixScalar%p_rspatialDiscrTest => rdiscretisationTest
    rmatrixScalar%bidenticalTrialAndTest = .false.
  else
    rmatrixScalar%p_rspatialDiscrTest => rdiscretisationTrial
    rmatrixScalar%bidenticalTrialAndTest = .true.
  end if
  rmatrixScalar%cmatrixFormat = LSYSSC_MATRIX9
  
  ! Get the #DOF's of the test space - as #DOF's of the test space is
  ! the number of equations in our matrix. The #DOF's in the trial space
  ! gives the number of columns of our matrix.
  rmatrixScalar%NCOLS         = &
      dof_igetNDofGlob(rmatrixScalar%p_rspatialDiscrTrial)
  rmatrixScalar%NEQ           = &
      dof_igetNDofGlob(rmatrixScalar%p_rspatialDiscrTest)
  
  ! and get a pointer to the triangulation.
  p_rtriangulation => rdiscretisationTrial%p_rtriangulation
  
  ! Get NEQ - we need it for guessing memory...
  NEQ = rmatrixScalar%NEQ
  
  if (NEQ .eq. 0) then
    call output_line ('Empty matrix!', &
        OU_CLASS_ERROR,OU_MODE_STD,'bilf_createMatrixStructure9_uni')
    call sys_halt()
  end if
  
  ! Allocate KLD...
  call storage_new1D ('bilf_createMatStructure9_conf', 'KLD', &
                      NEQ+1, ST_INT, rmatrixScalar%h_KLD, &
                      ST_NEWBLOCK_NOINIT)
  ! This must be a storage_getbase, no lsyssc_getbase, since this is the
  ! matrix construction routine!
  call storage_getbase_int(rmatrixScalar%h_Kld,p_KLD)
  
  ! Allocate h_Kdiagonal
  call storage_new1D ('bilf_createMatStructure9_conf', 'Kdiagonal', &
                      NEQ, ST_INT, rmatrixScalar%h_Kdiagonal, &
                      ST_NEWBLOCK_NOINIT)
  ! This must be a storage_getbase, no lsyssc_getbase, since this is the
  ! matrix construction routine!
  call storage_getbase_int (rmatrixScalar%h_Kdiagonal,p_Kdiagonal)
  
  ! For saving some memory in smaller discretisations, we calculate
  ! the number of elements per block. For smaller triangulations,
  ! this is NEL. If there are too many elements, it's at most
  ! BILF_NELEMSIM. This is only used for allocating some arrays.
  nelementsPerBlock = min(BILF_NELEMSIM,p_rtriangulation%NEL)

  ! Allocate a list of handles and a list of pointers corresponding to it.
  ! Initially allocate NmemBlkCount pointers
  allocate(p_Ihcol(NmemBlkCount))
  allocate(p_Ihindx(NmemBlkCount))
  allocate(p_Isize(NmemBlkCount))
  allocate(Rmemblock(NmemBlkCount))
  
  ! Allocate the first memory block that receives a part of the
  ! temporary matrix structure.
  ! We make an initial guess of iblkSize*iblkSize*NEQ elements in the matrix,
  ! if imemguess is not given.
  if (imemguess .ne. 0) then 
    ! at least one element per line!
    iallocated = max(NEQ,imemguess) 
  else  
    iallocated = max(10000,iblkSize*iblkSize*NEQ)
  end if
  iblocks = 1
  imemblkSize = iblkSize*NEQ
  p_Isize(1) = iallocated
  
  ! imemblkSize = iallocated is necessary at the moment to simplify
  ! whether we leave a block or not.

  call storage_new1D ('bilf_createMatStructure9_conf', 'Ihicol', &
                      p_Isize(1), ST_INT, p_Ihcol(1), ST_NEWBLOCK_NOINIT)
  call storage_getbase_int (p_Ihcol(1),p_Icol)

  ! The new index array must be filled with 0 - otherwise
  ! the search routine below won't work!
  call storage_new1D ('bilf_createMatStructure9_conf', 'p_Ihindx', &
                      p_Isize(1), ST_INT, p_Ihindx(1), ST_NEWBLOCK_ZERO)
  call storage_getbase_int (p_Ihindx(1),p_Iindx)
  
  Rmemblock(1)%p_Icol => p_Icol
  Rmemblock(1)%p_Iindx => p_Iindx
  
  ! Initialise Iindx and Icol.
  ! Initially, we have only diagonal elements in our matrix.
  !
  ! The basic idea behind the building of the matrix is a linked
  ! list of column numbers in each row!
  ! We collect all upcoming columns in the whole matrix in the
  ! array Icol, i.e. each new entry is attached to that
  ! (-> the array is resorted to KCOL at the end).
  ! Iindx points for every entry in the matrix to the position
  ! inside of Icol of the next entry in the line.
  !
  ! At the beginning, we no entries in the matrix. 
  ! We initialise the "head" of this collection of linked
  ! lists with 0 to indicate this. 
  ! Iindx(IEQ) is set to 0 to indicate that there is no following
  ! element in each line, i.e. we only have diagonal elements.
  ! Later, we fill Icol(1..NEQ) with the first column number in
  ! each row. When more entries appear in a row, they are appended
  ! to KCOL at position NEQ+1..*. Iindx keeps track of the entries
  ! in each row by storing (starting from the head in Icol(IEQ))
  ! the positions of the corresponding next entry inside of KCOL1
  ! in each row - so the column numbers in each row are to be
  ! found in KCOL1 at positions IEQ,Iindx(IEQ),Iindx(Iindx(IEQ)),...
  !
  ! Example: We want to add entry (1,3) to the matrix. Then
  ! we enlarge the lists as follows:
  ! - "2" (the column number" is added to the end of Icol, i.e.
  !   Icol(NEQ+1) = 3
  ! - We add a "follower" for the diagonal element by setting
  !   Iindx(1) = NEQ+1. Furthermore we set KIND(NEQ+1)=0
  ! So we have a linked list of matrix entries for the rows:
  !   Icol:     1   2   3   ...   NEQ     3
  !   Iindx:  NEQ+1 0   0   ...    0      0
  !             |                        /:\
  !             +-------------------------|
  ! i.e. row 1 can be computed as:
  !   Icol(1)               (=1),
  !   Icol(Iindx(1))        (=3),
  !   Icol(Iindx(Iindx(1))  -> not defined, as Iindx(Iindx(1))=0, 
  !                            line finished
  
  do IEQ=1,NEQ
    p_Iindx(IEQ) = 0
    p_Icol(IEQ)  = 0
  end do
  
  ! The first NEQ elements are reserved. The matrix is assumed
  ! to have at least one element per line.
  NA = NEQ
  
  ! Now loop over the different element distributions (=combinations
  ! of trial and test functions) in the discretisation.
  
  do icurrentElementDistr = 1,rdiscretisationTrial%inumFESpaces
  
    ! Activate the current element distribution
    p_relementDistrTest => &
        rmatrixScalar%p_rspatialDiscrTest%RelementDistr(icurrentElementDistr)
    p_relementDistrTrial => &
        rmatrixScalar%p_rspatialDiscrTrial%RelementDistr(icurrentElementDistr)

    ! Cancel if this element distribution is empty.
    if (p_relementDistrTest%NEL .eq. 0) cycle

    ! Get the number of local DOF's for trial and test functions
    indofTrial = elem_igetNDofLoc(p_relementDistrTrial%celement)
    indofTest = elem_igetNDofLoc(p_relementDistrTest%celement)
    
    ! Get the number of corner vertices of the element
    NVE = elem_igetNVE(p_relementDistrTest%celement)
    if (NVE .ne. elem_igetNVE(p_relementDistrTrial%celement)) then
      call output_line ('Element spaces incompatible!', &
          OU_CLASS_ERROR,OU_MODE_STD,'bilf_createMatStructure9_conf')
      call sys_halt()
    end if
    
    ! Allocate an array saving a couple of DOF's for trial and test functions
    !ALLOCATE(IdofsTrial(indofTrial,nelementsPerBlock))
    !ALLOCATE(IdofsTest(indofTest,nelementsPerBlock))
    
    ! Test if trial/test functions are identical.
    ! We don't rely on bidenticalTrialAndTest purely, as this does not
    ! indicate whether there are identical trial and test functions
    ! in one block!
    bIdenticalTrialAndTest = &
      p_relementDistrTest%celement .eq. p_relementDistrTrial%celement
      
    ! Let p_IdofsTrial point either to IdofsTrial or to the DOF's of the test
    ! space IdofTest (if both spaces are identical). 
    ! We create a pointer for the trial space and not for the test space to
    ! prevent pointer-arithmetic in the innerst loop below!
    if (bIdenticalTrialAndTest) then
      p_IdofsTrial => IdofsTest
    else
      p_IdofsTrial => IdofsTrial
    end if
    
    ! p_IelementList must point to our set of elements in the discretisation
    ! with that combination of trial/test functions
    call storage_getbase_int (p_relementDistrTest%h_IelementList, &
                              p_IelementList)
    
    ! Get the number of elements there.
    NEL = p_relementDistrTest%NEL

    ! Set the pointers/indices to the initial position. During the
    ! search for new DOF's, these might be changed if there's not enough
    ! memory in the first block.    
    icurrentblock = 1
    istartidx = 0
    p_Icol => Rmemblock(1)%p_Icol
    p_Iindx => Rmemblock(1)%p_Iindx
    
    ! Loop over the elements. 
    do IELset = 1, NEL, BILF_NELEMSIM
    
      ! We always handle BILF_NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most BILF_NELEMSIM
      ! elements simultaneously.
      
      IELmax = min(NEL,IELset-1+BILF_NELEMSIM)
    
      ! The outstanding feature with finite elements is: A basis
      ! function for a DOF on one element has common support only
      ! with the DOF's on the same element! E.g. for Q1:
      !
      !        #. . .#. . .#. . .#
      !        .     .     .     .
      !        .  *  .  *  .  *  .
      !        #-----O-----O. . .#
      !        |     |     |     .
      !        |     | IEL |  *  .
      !        #-----X-----O. . .#
      !        |     |     |     .
      !        |     |     |  *  .
      !        #-----#-----#. . .#
      !        
      ! --> On element IEL, the basis function at "X" only interacts
      !     with the basis functions in "O". Elements in the 
      !     neighbourhood ("*") have no support, therefore we only have
      !     to collect all "O" DOF's.
      !
      ! Call dof_locGlobMapping to get the global DOF's on our current
      ! element (the "X" and all "O"'s). 
      ! We don't need the local DOF's, so by setting IPAR=0,
      ! the call will only fill KDFG.
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF's of our BILF_NELEMSIM elements simultaneously.
      ! Calculate the DOF's of the test functions:
      call dof_locGlobMapping_mult(&
          rmatrixScalar%p_rspatialDiscrTest, p_IelementList(IELset:IELmax), &
          IdofsTest)
                                   
      ! If the DOF's for the test functions are different, calculate them, too.
      if (.not.bIdenticalTrialAndTest) then
        call dof_locGlobMapping_mult(&
            rmatrixScalar%p_rspatialDiscrTrial, p_IelementList(IELset:IELmax), &
            IdofsTrial)
      end if
      
      ! Loop through all the elements in the current set
      do IEL=1,IELmax-IELset+1
        
        ! For building the local matrices, we have first to
        ! loop through the test functions (the "O"'s), as these
        ! define the rows in the matrix.
        do IDOFE=1,indofTest

          ! The DOF IDOFE is now our "O".
          ! This global DOF gives us the row we have to build.
          
          IROW = IdofsTest(IDOFE,IEL)
          
          ! Now we loop through the other DOF's on the current element
          ! (the "X"'s).
          ! All these have common support with our current basis function
          ! and will therefore give an additive value to the global
          ! matrix.

          do JDOFE=1,indofTrial
            
            ! Get the global DOF - our "X". This gives the column number
            ! in the matrix where an entry occurs in row IROW (the line of 
            ! the current global DOF "O").
              
            JCOL = p_IdofsTrial(JDOFE,IEL)

            ! This JCOL has to be inserted into line IROW.
            ! But first check, whether the element is already in that line,
            ! i.e. whether element (IROW,JCOL) is already in the matrix.
            ! This may happen because of an earlier loop in a neighbour
            ! element (imagine Q1, where there are two vertices on an edge
            ! and the edge is shared between two elements)...
            !
            ! We start walking through the linked list of row IROW to 
            ! look for column JCOL. IPOS is the position in the KCOL1 
            ! array of the column in row IROW we want to test. 
            ! KINDX(IPOS) is =0 if we reach the end of the linked list.
            !
            ! Start searching at the "head" of the list, which is the
            ! diagonal element at Icol(IROW).
            ! This is always found in the first memory block.
      
            ! Remark: This IF command gains 8% performance due to slow
            ! pointer handling!
            if (icurrentblock .ne. 1) then
              icurrentblock = 1
              istartidx = 0
              p_Icol => Rmemblock(1)%p_Icol
              p_Iindx => Rmemblock(1)%p_Iindx
            end if
            
            ! Is the list empty?

            if (p_Icol(IROW).eq.0) then

              ! Yes, row IROW is empty at the moment. Add the column as
              ! head of the list of row IROW.

              p_Icol(IROW) = JCOL
              
            else

              ! No, the list is not empty, we have a "head".
              !
              ! We start walking through the linked list of row IROW to 
              ! look for column JCOL. IPOS is the position in the KCOL1 
              ! array of the column in row IROW we want to test. 
              ! KINDX(IPOS) is =0 if we reach the end of the linked list.
              !
              ! Start searching at the "head" of the list, which is the
              ! diagonal element at p_Icol(IROW).

              IPOS=IROW
            
              ! Loop through the elements in the list until we find
              ! column JCOL - or the end of the list!
              ! IPOS must be corrected by istartidx, which is the number
              ! of elements in all blocks before the current one.
              
              searchloop: do while ( (p_Icol(IPOS-istartIdx)) .ne. JCOL)
            
                ! Did we reach the end of the list? Then we have to insert
                ! a new element...
                if (p_Iindx(IPOS-istartIdx) .eq. 0) then
                
                  ! Increase NA, which is the actual length of the matrix -
                  ! and at the same time tells us how much memory we actually
                  ! use.
                  NA = NA+1
                  
                  ! Let p_Iindx of the last element of the row point to our
                  ! new element. The new element is now the last in the row.
                
                  p_Iindx(IPOS-istartIdx) = NA
                  
                  ! Before really appending JCOL, first test
                  ! NA is now larger than the marimum amount
                  ! of storage we allocated!
                  
                  if (NA .gt. iallocated) then
                   
                    ! Hmmm, we have to allocate more memory.
                    ! Do we have enough pointers left or do we have
                    ! to enlarge our list?
                    
                    if (iblocks .ge. size(p_Ihcol)) then 
                    
                      ! Not enough blocks, we have to reallocate the pointer lists!
                      allocate (p_IhTmp(iblocks+NmemBlkCount))
                      p_IhTmp(1:iblocks) = p_Ihcol(1:iblocks)
                      deallocate(p_Ihcol)
                      p_Ihcol => p_IhTmp

                      allocate (p_IhTmp(iblocks+NmemBlkCount))
                      p_IhTmp(1:iblocks) = p_Ihindx(1:iblocks)
                      deallocate(p_Ihindx)
                      p_Ihindx => p_IhTmp
                    
                      allocate (p_IsizeTmp(iblocks+NmemBlkCount))
                      p_IsizeTmp(1:iblocks) = p_Isize(1:iblocks)
                      deallocate(p_Isize)
                      p_Isize => p_IsizeTmp

                      allocate (RmemblockTmp(iblocks+NmemBlkCount))
                      RmemblockTmp(1:iblocks) = Rmemblock(1:iblocks)
                      deallocate(Rmemblock)
                      Rmemblock => RmemblockTmp
                      
                      ! Now we have enough blocks again.
                    end if

                    ! Add a new block

                    iblocks = iblocks + 1
                    p_Isize (iblocks) = imemblkSize
                    
                    ! Move the start index behind the last completely
                    ! occupied block
                             
                    istartIdx = iallocated
                    icurrentblock = iblocks
                  
                    ! Allocate a new memory block of size imemblkSize
                    !
                    ! Use p_Icol and p_Iindx - they are not used anymore.
                    ! Allocate another imemblkSize elements for column numbers and
                    ! list pointers.

                    call storage_new1D ('bilf_createMatStructure9_conf', 'Ihicol', &
                                        p_Isize (iblocks), ST_INT, p_Ihcol(iblocks), &
                                        ST_NEWBLOCK_NOINIT)
                    call storage_getbase_int (p_Ihcol(iblocks),p_Icol)

                    ! The new index array must be filled with 0 - otherwise
                    ! the search routine below won't work!
                    call storage_new1D ('bilf_createMatStructure9_conf', 'p_Ihindx', &
                                        p_Isize (iblocks), ST_INT, p_Ihindx(iblocks), &
                                        ST_NEWBLOCK_ZERO)
                    call storage_getbase_int (p_Ihindx(iblocks),p_Iindx)
                    
                    Rmemblock(iblocks)%p_Icol => p_Icol
                    Rmemblock(iblocks)%p_Iindx => p_Iindx

                    iallocated = iallocated + p_Isize (iblocks)
                    
                  else
                    ! Be careful when leaving the current memory block
                    ! for insertion of an element at position NA!
                    !
                    ! If the new position is not in the current block,
                    ! it's in the last block... and so set the pointer
                    ! and indices appropriately!
                    
                    if ( NA .gt. (istartidx+p_Isize (icurrentblock))) then
                      istartidx = iallocated-p_Isize(iblocks)
                      icurrentblock = iblocks
                      p_Icol => Rmemblock(iblocks)%p_Icol
                      p_Iindx => Rmemblock(iblocks)%p_Iindx
                    end if
                    
                  end if
                
                  ! Append JCOL to p_Icol
                  p_Icol(NA-istartIdx) = JCOL
                  
                  ! We have to make sure that p_Indx(NA)=0 to indicate the end of
                  ! the list. Ok, this is trivial because we allocated it with 
                  ! storage_new, configured to fill the memory with 0, so it is 0.
                  !
                  ! The searchloop ends here, continue with next JDOFE
                  
                  exit 
                
                else
                
                  ! No, the list does not end here.
                  ! Take the next element in the list
                  IPOS = p_Iindx(IPOS-istartidx)
                  
                  ! Be careful when leaving the current memory block
                  do while ( IPOS .gt. (istartidx+p_Isize (icurrentblock)) )
                  
                    ! go to the next memory block and search there
                    istartidx = istartidx+p_Isize(icurrentblock)
                    icurrentblock = icurrentblock+1
                    p_Icol => Rmemblock(icurrentblock)%p_Icol
                    p_Iindx => Rmemblock(icurrentblock)%p_Iindx
                    
                  end do ! IPOS .GT. (istartidx+p_Isize (iblocks))
                
                end if ! p_Iindx(IPOS) = 0
                
              end do searchloop
            
            end if ! p_Icol(IROW) = 0
                
          end do ! JDOFE
        
        end do ! IDOFE
      
      end do ! IEL
    
    end do ! IELset

    ! Clean up the DOF's arrays    
    !DEALLOCATE(IdofsTest)
    !DEALLOCATE(IdofsTrial)
    
  end do ! icurrentElementDistr
  
  ! Ok, p_Icol is built. The hardest part is done!
  ! Now build KCOL by collecting the entries in the linear lists of 
  ! each row.
  !
  ! At first, as we now NA, we can allocate the real KCOL now!
  
  call storage_new1D ('bilf_createMatStructure9_conf', 'KCOL', &
                      NA, ST_INT, rmatrixScalar%h_KCOL, &
                      ST_NEWBLOCK_NOINIT)
  ! This must be a storage_getbase, no lsyssc_getbase, since this is the
  ! matrix construction routine!
  call storage_getbase_int (rmatrixScalar%h_Kcol,p_KCOL)
  
  ! Save NA in the matrix structure
  rmatrixScalar%NA = NA
  
  ! Set back NA to 0 at first.

  NA=0
      
  ! Loop through all of the NEQ linear lists:
  
  do IEQ=1,NEQ
  
    ! We are at the head of the list, now we have to walk
    ! through it to append the entries to KCOL.
    ! We always start in the first memory block.

    if (icurrentblock .ne. 1) then
      icurrentblock = 1
      istartidx = 0
      p_Icol => Rmemblock(1)%p_Icol
      p_Iindx => Rmemblock(1)%p_Iindx
    end if
    
    ! Add the head of the list to KCOL:
  
    NA=NA+1
    p_KCOL(NA)=p_Icol(IEQ)

    ! Set KLD appropriately:

    p_KLD(IEQ) = NA
    
    IPOS = IEQ
    
    do while (p_Iindx(IPOS-istartidx).ne.0)
    
      ! Get the position of the next entry in p_Icol:
      IPOS = p_Iindx(IPOS-istartidx)
      
      ! Be careful when leaving the current memory block
      do while ( IPOS .gt. (istartidx+p_Isize (icurrentblock)) )
      
        ! go to the next memory block and search there
        istartidx = istartidx+p_Isize(icurrentblock)
        icurrentblock = icurrentblock+1
        p_Icol => Rmemblock(icurrentblock)%p_Icol
        p_Iindx => Rmemblock(icurrentblock)%p_Iindx
        
      end do ! IPOS .GT. (istartidx+p_Isize (iblocks))
      
      ! Add the column number to the row in KCOL:
      NA=NA+1
      p_KCOL(NA)=p_Icol(IPOS-istartidx)
    
    end do ! KINDX(IPOS) <> 0

  end do ! IEQ
  
  ! Append the final entry to KLD:
  
  p_KLD(NEQ+1)=NA+1
  
  ! Sort entries on KCOL separately for each row.
  ! This is a small bubble-sort...
  !
  ! Loop through all rows:
  nmaxCol = 0

  do IEQ=1,NEQ

    ! Repeat until everything is sorted.
    
    BSORT=.false.
    do while (.not. BSORT)
    
      BSORT=.true.

      !  Loop through the line 

      do JCOL=p_KLD(IEQ),p_KLD(IEQ+1)-2
      
        ! If the next element is larger...
      
        if (p_KCOL(JCOL) .gt. p_KCOL(JCOL+1)) then
        
          ! Change position of the current and next element
        
          IHELP=p_KCOL(JCOL)
          p_KCOL(JCOL)=p_KCOL(JCOL+1)
          p_KCOL(JCOL+1)=IHELP
          
          ! And repeat the sorting of that line
          
          BSORT=.false.
          
        end if
        
      end do ! JCOL
      
    end do ! (not BSORT)      

    ! Grab the diagonal
    do JCOL=p_KLD(IEQ),p_KLD(IEQ+1)-1
      if (p_KCOL(JCOL) .ge. IEQ) then
        p_Kdiagonal(IEQ) = JCOL
        exit
      end if
    end do   
    
    ! Grab the largest column number. As the current line is sorted,
    ! we can find this using the end of the line.
    nmaxCol = max(nmaxCol,p_Kcol(p_Kld(IEQ+1)-1))

  end do ! IEQ
  
  ! HOORAY, THAT'S IT!
  ! Deallocate all temporary memory...
  
  do i=iblocks,1,-1
    call storage_free(p_Ihcol(i))
    call storage_free(p_Ihindx(i))
  end do
  
  deallocate(Rmemblock)
  deallocate(p_Isize)
  deallocate(p_Ihindx)
  deallocate(p_Ihcol)
  
  end subroutine

  !****************************************************************************
  
!<subroutine>
  
  subroutine bilf_createMatStructure9eb_uni (rdiscretisationTrial,rmatrixScalar,&
      rdiscretisationTest,imemGuess)
  
!<description>
  ! This routine creates according to a given discretisation the matrix 
  ! structure of a structure-9 matrix. The discretisation is assumed to be
  ! uniform, i.e. there is only one combination of test- and trial functions
  ! allowed. The matrix is created by an edge-based approach, which increases
  ! the matrix stencil in such a way, that the DOF's of one element may
  ! interact with the DOF's of all elements that are adjacent via the
  ! edges.
!</description>

!<input>
  ! The underlying discretisation structure which is to be used to
  ! create the matrix. Specifies the discretisation of the trial functions.
  ! If rdiscretisationTest is not specified, this also specifies
  ! the discretisation of the test functions, this test and trial
  ! functions coincide.
  type(t_spatialDiscretisation), intent(IN), target :: rdiscretisationTrial
  
  ! OPTIONAL: The underlying discretisation structure for the test functions.
  ! If not specified, the trial functions coincide with the test functions.
  type(t_spatialDiscretisation), intent(IN), target, optional :: rdiscretisationTest
  
  ! An initial guess about how much memory the matrix needs. If set to 0,
  ! an initial guess of 16*NEQ (but at least 10000 matrix entries) is assumed.
  integer, intent(IN) :: imemGuess

!</input>

!<output>
  ! The structure of a scalar matrix, fitting to the given discretisation.
  ! Memory fo rthe structure is allocated dynamically on the heap.
  type(t_matrixScalar), intent(OUT) :: rmatrixScalar
!</output>

!</subroutine>

  ! local variables
  integer :: NEQ, IEQ, IROW, JCOL, IPOS, istartIdx, NA, nmaxCol
  integer :: IDOFE, JDOFE, i, IHELP,NVE, nelemBlockCount, IELidx
  integer :: IELneighIdxJ
  integer :: IEL, IELmax, IELset
  logical :: BSORT, bIdenticalTrialAndTest
  
  ! An allocateable list of handles for memory blocks. Size is dynamically 
  ! increased if there are too many columns in the matrix.
  integer, dimension(:), pointer :: p_Ihcol, p_Ihindx, p_IhTmp
  integer, dimension(:), pointer :: p_Isize, p_ISizeTmp
  
  ! An allocateable list of pointers to the memory blocks - corresponds
  ! to the handles in p_Ihcol/p_Ihindx
  type(t_matrixmem), dimension(:), pointer :: Rmemblock, RmemblockTmp
  
  ! Pointer to current KCOL memory block,
  ! pointer to current index memory block
  integer, dimension(:), pointer :: p_Icol, p_Iindx
  
  ! Number of currently allocated pointers in Ihmemblock
  integer :: iblocks
  
  ! Currently active memory block
  integer :: icurrentblock
  
  ! Size of memory blocks
  integer :: imemblkSize
  
  ! Blocksize in terms of NEQ for guessing memory.
  ! The initial guess for memory is iblkSize*iblkSize*NEQ and every time
  ! memory is needed, another iblkSize*NEQ elements are added.
  integer, parameter :: iblkSize = 4
  
  ! Number of memory blocks to allocate
  integer, parameter :: NmemBlkCount = 5

  ! Pointer to KLD, KCOL, diagonal
  integer, dimension(:), pointer :: p_KLD, p_KCOL, p_Kdiagonal
  
  ! Size of memory currently allocated
  integer :: iallocated
  
  ! An allocateable array accepting the DOF's of a set of elements.
  integer, dimension(:,:), allocatable, target :: IdofsTest, IdofsTrial
  integer, dimension(:,:), pointer :: p_IdofsTrial
  !INTEGER, DIMENSION(EL_MAXNBAS,BILF_NELEMSIM), TARGET :: IdofsTest, IdofsTrial
  !INTEGER, DIMENSION(:,:), POINTER :: p_IdofsTrial
  
  integer, dimension(:), allocatable :: IadjPtr, IadjElem
  
  ! Number of local degees of freedom for trial and test functions
  integer :: indofTrial, indofTest
  
  ! The triangulation structure - to shorten some things...
  type(t_triangulation), pointer :: p_rtriangulation
  
  ! A pointer to an element-number list
  integer, dimension(:), pointer :: p_IelementList
  
  ! Current element distribution
  type(t_elementDistribution), pointer :: p_relementDistrTest
  type(t_elementDistribution), pointer :: p_relementDistrTrial

  ! Number of elements in a block. Normally =BILF_NELEMSIM,
  ! except if there are less elements in the discretisation.
  integer :: nelementsPerBlock
  
  ! Adjacent elements
  integer, dimension(:,:), pointer :: p_Kadj

  ! The algorithm is: Test every DOF on one element against each other
  ! DOF on the same element and save the combination into a matrix
  ! in structure 9!
  !
  ! At first, initialise the structure-9 matrix:
  
  rmatrixScalar%p_rspatialDiscrTrial => rdiscretisationTrial
  if (present(rdiscretisationTest)) then
    rmatrixScalar%p_rspatialDiscrTest => rdiscretisationTest
    rmatrixScalar%bidenticalTrialAndTest = .false.
  else
    rmatrixScalar%p_rspatialDiscrTest => rdiscretisationTrial
    rmatrixScalar%bidenticalTrialAndTest = .true.
  end if
  rmatrixScalar%cmatrixFormat = LSYSSC_MATRIX9
  
  ! Get the #DOF's of the test space - as #DOF's of the test space is
  ! the number of equations in our matrix. The #DOF's in the trial space
  ! gives the number of columns of our matrix.
  rmatrixScalar%NCOLS         = &
      dof_igetNDofGlob(rmatrixScalar%p_rspatialDiscrTrial)
  rmatrixScalar%NEQ           = &
      dof_igetNDofGlob(rmatrixScalar%p_rspatialDiscrTest)
  
  ! and get a pointer to the triangulation.
  p_rtriangulation => rdiscretisationTrial%p_rtriangulation
  
  call storage_getbase_int2d (p_rtriangulation%h_IneighboursAtElement,p_Kadj)
  
  ! Get NEQ - we need it for guessing memory...
  NEQ = rmatrixScalar%NEQ
  
  if (NEQ .eq. 0) then
    call output_line ('Empty matrix!', &
        OU_CLASS_ERROR,OU_MODE_STD,'bilf_createMatStructure9_conf')
    call sys_halt()
  end if
  
  ! Allocate KLD...
  call storage_new1D ('bilf_createMatStructure9_uni', 'KLD', &
                      NEQ+1, ST_INT, rmatrixScalar%h_KLD, &
                      ST_NEWBLOCK_NOINIT)
  ! This must be a storage_getbase, no lsyssc_getbase, since this is the
  ! matrix construction routine!
  call storage_getbase_int (rmatrixScalar%h_Kld,p_KLD)
  
  ! Allocate h_Kdiagonal
  call storage_new1D ('bilf_createMatStructure9_conf', 'Kdiagonal', &
                      NEQ, ST_INT, rmatrixScalar%h_Kdiagonal, &
                      ST_NEWBLOCK_NOINIT)
  ! This must be a storage_getbase, no lsyssc_getbase, since this is the
  ! matrix construction routine!
  call storage_getbase_int (rmatrixScalar%h_Kdiagonal,p_Kdiagonal)
  
  ! For saving some memory in smaller discretisations, we calculate
  ! the number of elements per block. For smaller triangulations,
  ! this is NEL. If there are too many elements, it's at most
  ! BILF_NELEMSIM. This is only used for allocating some arrays.
  nelementsPerBlock = min(BILF_NELEMSIM,p_rtriangulation%NEL)

  ! Allocate a list of handles and a list of pointers corresponding to it.
  ! Initially allocate NmemBlkCount pointers
  allocate(p_Ihcol(NmemBlkCount))
  allocate(p_Ihindx(NmemBlkCount))
  allocate(p_Isize(NmemBlkCount))
  allocate(Rmemblock(NmemBlkCount))
  
  ! Allocate the first memory block that receives a part of the
  ! temporary matrix structure.
  ! We make an initial guess of iblkSize*iblkSize*NEQ elements in the matrix,
  ! if imemguess is not given.
  if (imemguess .ne. 0) then 
    ! at least one element per line!
    iallocated = max(NEQ,imemguess) 
  else  
    iallocated = max(10000,iblkSize*iblkSize*NEQ)
  end if
  iblocks = 1
  imemblkSize = iblkSize*NEQ
  p_Isize(1) = iallocated
  
  ! imemblkSize = iallocated is necessary at the moment to simplify
  ! whether we leave a block or not.

  call storage_new1D ('bilf_createMatStructure9_conf', 'Ihicol', &
                      p_Isize(1), ST_INT, p_Ihcol(1), ST_NEWBLOCK_NOINIT)
  call storage_getbase_int (p_Ihcol(1),p_Icol)

  ! The new index array must be filled with 0 - otherwise
  ! the search routine below won't work!
  call storage_new1D ('bilf_createMatStructure9_conf', 'p_Ihindx', &
                      p_Isize(1), ST_INT, p_Ihindx(1), ST_NEWBLOCK_ZERO)
  call storage_getbase_int (p_Ihindx(1),p_Iindx)
  
  Rmemblock(1)%p_Icol => p_Icol
  Rmemblock(1)%p_Iindx => p_Iindx
  
  ! Initialise Iindx and Icol.
  ! Initially, we have only diagonal elements in our matrix.
  !
  ! The basic idea behind the building of the matrix is a linked
  ! list of column numbers in each row!
  ! We collect all upcoming columns in the whole matrix in the
  ! array Icol, i.e. each new entry is attached to that
  ! (-> the array is resorted to KCOL at the end).
  ! Iindx points for every entry in the matrix to the position
  ! inside of Icol of the next entry in the line.
  !
  ! At the beginning, we no entries in the matrix. 
  ! We initialise the "head" of this collection of linked
  ! lists with 0 to indicate this. 
  ! Iindx(IEQ) is set to 0 to indicate that there is no following
  ! element in each line, i.e. we only have diagonal elements.
  ! Later, we fill Icol(1..NEQ) with the first column number in
  ! each row. When more entries appear in a row, they are appended
  ! to KCOL at position NEQ+1..*. Iindx keeps track of the entries
  ! in each row by storing (starting from the head in Icol(IEQ))
  ! the positions of the corresponding next entry inside of KCOL1
  ! in each row - so the column numbers in each row are to be
  ! found in KCOL1 at positions IEQ,Iindx(IEQ),Iindx(Iindx(IEQ)),...
  !
  ! Example: We want to add entry (1,3) to the matrix. Then
  ! we enlarge the lists as follows:
  ! - "2" (the column number" is added to the end of Icol, i.e.
  !   Icol(NEQ+1) = 3
  ! - We add a "follower" for the diagonal element by setting
  !   Iindx(1) = NEQ+1. Furthermore we set KIND(NEQ+1)=0
  ! So we have a linked list of matrix entries for the rows:
  !   Icol:     1   2   3   ...   NEQ     3
  !   Iindx:  NEQ+1 0   0   ...    0      0
  !             |                        /:\
  !             +-------------------------|
  ! i.e. row 1 can be computed as:
  !   Icol(1)               (=1),
  !   Icol(Iindx(1))        (=3),
  !   Icol(Iindx(Iindx(1))  -> not defined, as Iindx(Iindx(1))=0, 
  !                            line finished
  
  do IEQ=1,NEQ
    p_Iindx(IEQ) = 0
    p_Icol(IEQ)  = 0
  end do
  
  ! The first NEQ elements are reserved. The matrix is assumed
  ! to have at least one element per line.
  NA = NEQ
  
  ! Activate the one and only element distribution
  p_relementDistrTest => rmatrixScalar%p_rspatialDiscrTest%RelementDistr(1)
  p_relementDistrTrial => rmatrixScalar%p_rspatialDiscrTrial%RelementDistr(1)

  ! Get the number of local DOF's for trial and test functions
  indofTrial = elem_igetNDofLoc(p_relementDistrTrial%celement)
  indofTest = elem_igetNDofLoc(p_relementDistrTest%celement)
  
  ! Get the number of corner vertices of the element
  NVE = elem_igetNVE(p_relementDistrTrial%celement)
  if (NVE .ne. elem_igetNVE(p_relementDistrTest%celement)) then
    call output_line ('Element spaces incompatible!', &
        OU_CLASS_ERROR,OU_MODE_STD,'bilf_createMatStructure9_conf')
    call sys_halt()
  end if
  
  ! Allocate the IadjCount array. This array counts for every element,
  ! how many elements are adjacent to that.
  allocate(IadjPtr(nelementsPerBlock+1))
  
  ! Allocate the IadjElem array. This collects for every element 
  ! the element number itself as well as all the adjacent elements. 
  ! (It's like a KLD-array...)
  ! IadjCount is a pointer into this array, so we can access directly 
  ! the numbers of the elements adjacent to one element.
  ! There is
  !   IadjElem(IadjPtr(IEL)) = IEL
  !   IadjElem(IadjPtr(IEL)+1..IadjPtr(IEL+1)-1) = adjacent elements
  ! As we know the maximum number of edges, we know the maximum size
  ! this array may need.
  allocate(IadjElem(nelementsPerBlock*(ubound(p_Kadj,1)+1)))

  ! Allocate an array saving a couple of DOF's for trial and test functions
  allocate(IdofsTrial(indofTrial,nelementsPerBlock*(ubound(p_Kadj,1)+1)))
  allocate(IdofsTest(indofTest,nelementsPerBlock*(ubound(p_Kadj,1)+1)))
  
  ! Test if trial/test functions are identical.
  ! We don't rely on bidenticalTrialAndTest purely, as this does not
  ! indicate whether there are identical trial and test functions
  ! in one block!
  bIdenticalTrialAndTest = &
    p_relementDistrTest%celement .eq. p_relementDistrTrial%celement
    
  ! Let p_IdofsTrial point either to IdofsTrial or to the DOF's of the test
  ! space IdofTest (if both spaces are identical). 
  ! We create a pointer for the trial space and not for the test space to
  ! prevent pointer-arithmetic in the innerst loop below!
  if (bIdenticalTrialAndTest) then
    p_IdofsTrial => IdofsTest
  else
    p_IdofsTrial => IdofsTrial
  end if
  
  ! p_IelementList must point to our set of elements in the discretisation
  ! with that combination of trial/test functions
  call storage_getbase_int (p_relementDistrTest%h_IelementList, &
                            p_IelementList)
  

  ! Set the pointers/indices to the initial position. During the
  ! search for new DOF's, these might be changed if there's not enough
  ! memory in the first block.    
  icurrentblock = 1
  istartidx = 0
  p_Icol => Rmemblock(1)%p_Icol
  p_Iindx => Rmemblock(1)%p_Iindx
  
  ! Loop over the elements. 
  do IELset = 1, p_rtriangulation%NEL, BILF_NELEMSIM
  
    ! We always handle BILF_NELEMSIM elements simultaneously.
    ! How many elements have we actually here?
    ! Get the maximum element number, such that we handle at most BILF_NELEMSIM
    ! elements simultaneously.
    
    IELmax = min(p_rtriangulation%NEL,IELset-1+BILF_NELEMSIM)
    
    ! --------------------- DOF SEARCH PHASE ------------------------
    
    ! In a first step, we search all the DOF's on one element
    ! and those on the elements adjacent via an edge to this.
    ! We want to get all the DOF numbers simultaneously!
    ! For this, we first set up an array that holds all the element
    ! numbers in our set as well as their neighbours.
    !
    ! Set up the IadjPtr array: Count 1 for the current element
    ! and count the number of adjacent elements to each element.
    ! Store the element number of one element to IadjElem and
    ! behind that the numbers of the adjacent elements.
    nelemBlockCount = 0
    do IELidx=1,IELmax-IELset+1
      IEL = IELidx+IELset-1      ! actual element number
      
      nelemBlockCount = nelemBlockCount+1
      IadjPtr(IELidx) = nelemBlockCount
      IadjElem(nelemBlockCount) = IEL
      
      do i=1,ubound(p_Kadj,1)
        if (p_Kadj(i,IEL) .ne. 0) then
          nelemBlockCount = nelemBlockCount+1
          IadjElem(nelemBlockCount) = p_Kadj(i,IEL)
        end if
      end do
      
    end do
    IadjPtr(IELmax-IELset+1+1) = nelemBlockCount+1
    
    ! nelemBlockCount is now the number of elements in the current
    ! block, consisting of the IELmax elements and their "edge-neighbours".
  
    ! The outstanding feature with finite elements is: A basis
    ! function for a DOF on one element has common support only
    ! with the DOF's on the same element! E.g. for Q1:
    !
    !        #. . .#. . .#. . .#
    !        .     .     .     .
    !        .  *  .  *  .  *  .
    !        #-----O-----O. . .#
    !        |     |     |     .
    !        |     | IEL |  *  .
    !        #-----X-----O. . .#
    !        |     |     |     .
    !        |     |     |  *  .
    !        #-----#-----#. . .#
    !        
    ! --> On element IEL, the basis function at "X" only interacts
    !     with the basis functions in "O". Elements in the 
    !     neighbourhood ("*") have no support, therefore we only have
    !     to collect all "O" DOF's.
    !
    ! Call dof_locGlobMapping to get the global DOF's on our current
    ! element (the "X" and all "O"'s). 
    ! We don't need the local DOF's, so by setting IPAR=0,
    ! the call will only fill KDFG.
    !
    ! More exactly, we call dof_locGlobMapping_mult to calculate all the
    ! global DOF's of our nelemBlockCount elements simultaneously.
    ! Calculate the DOF's of the test functions:
    call dof_locGlobMapping_mult(&
        rmatrixScalar%p_rspatialDiscrTest, IadjElem(1:nelemBlockCount), &
        IdofsTest)
                                 
    ! If the DOF's for the test functions are different, calculate them, too.
    if (.not.bIdenticalTrialAndTest) then
      call dof_locGlobMapping_mult(&
          rmatrixScalar%p_rspatialDiscrTrial, IadjElem(1:nelemBlockCount), &
          IdofsTrial)
    end if
  
    ! --------------------- DOF COMBINATION PHASE ------------------------
    
    ! Loop through all the elements in the current set
    do IEL=1,IELmax-IELset+1
    
      ! For building the local matrices, we have first to
      ! loop through the test functions (the "O"'s), as these
      ! define the rows in the matrix.
      do IDOFE=1,indofTest

        ! The DOF IDOFE is now our "O".
        ! This global DOF gives us the row we have to build.
        !
        ! The DOF's of element IEL start at position IadjPtr(IEL) in
        ! the IdofsTest array. 
        IROW = IdofsTest(IDOFE,IadjPtr(IEL)) 
        
        ! Loop through the "element-sets" in IadjElem, consisting
        ! of the element IEL itself as well as its neighbours - for the
        ! trial functions.
        do IELneighIdxJ = IadjPtr(IEL),IadjPtr(IEL+1)-1

          ! Now we loop through the other DOF's on the current element
          ! (the "X"'s).
          ! All these have common support with our current basis function
          ! and will therefore give an additive value to the global
          ! matrix.

          do JDOFE=1,indofTrial
            
            ! Get the global DOF - our "X". This gives the column number
            ! in the matrix where an entry occurs in row IROW (the line of 
            ! the current global DOF "O").
              
            JCOL = p_IdofsTrial(JDOFE,IELneighIdxJ)

            ! This JCOL has to be inserted into line IROW.
            ! But first check, whether the element is already in that line,
            ! i.e. whether element (IROW,JCOL) is already in the matrix.
            ! This may happen because of an earlier loop in a neighbour
            ! element (imagine Q1, where there are two vertices on an edge
            ! and the edge is shared between two elements)...
            !
            ! We start walking through the linked list of row IROW to 
            ! look for column JCOL. IPOS is the position in the KCOL1 
            ! array of the column in row IROW we want to test. 
            ! p_Iindx(IPOS) is =0 if we reach the end of the linked list.
            !
            ! Start searching at the "head" of the list, which is the
            ! diagonal element at Icol(IROW).
            ! This is always found in the first memory block.
      
            ! Remark: This IF command gains 8% performance due to slow
            ! pointer handling!
            if (icurrentblock .ne. 1) then
              icurrentblock = 1
              istartidx = 0
              p_Icol => Rmemblock(1)%p_Icol
              p_Iindx => Rmemblock(1)%p_Iindx
            end if
            
            ! Is the list empty?

            if (p_Icol(IROW).eq.0) then

              ! Yes, row IROW is empty at the moment. Add the column as
              ! head of the list of row IROW.

              p_Icol(IROW) = JCOL
              
            else

              ! No, the list is not empty, we have a "head".
              !
              ! We start walking through the linked list of row IROW to 
              ! look for column JCOL. IPOS is the position in the KCOL1 
              ! array of the column in row IROW we want to test. 
              ! KINDX(IPOS) is =0 if we reach the end of the linked list.
              !
              ! Start searching at the "head" of the list, which is the
              ! diagonal element at p_Icol(IROW).

              IPOS=IROW
            
              ! Loop through the elements in the list until we find
              ! column JCOL - or the end of the list!
              ! IPOS must be corrected by istartidx, which is the number
              ! of elements in all blocks before the current one.
              
              searchloop: do while ( (p_Icol(IPOS-istartIdx)) .ne. JCOL)
            
                ! Did we reach the end of the list? Then we have to insert
                ! a new element...
                if (p_Iindx(IPOS-istartIdx) .eq. 0) then
                
                  ! Increase NA, which is the actual length of the matrix -
                  ! and at the same time tells us how much memory we actually
                  ! use.
                  NA = NA+1
                  
                  ! Let p_Iindx of the last element of the row point to our
                  ! new element. The new element is now the last in the row.
                
                  p_Iindx(IPOS-istartIdx) = NA
                  
                  ! Before really appending JCOL, first test
                  ! NA is now larger than the marimum amount
                  ! of storage we allocated!
                  
                  if (NA .gt. iallocated) then
                   
                    ! Hmmm, we have to allocate more memory.
                    ! Do we have enough pointers left or do we have
                    ! to enlarge our list?
                    
                    if (iblocks .ge. size(p_Ihcol)) then 
                    
                      ! Not enough blocks, we have to reallocate the pointer lists!
                      allocate (p_IhTmp(iblocks+NmemBlkCount))
                      p_IhTmp(1:iblocks) = p_Ihcol(1:iblocks)
                      deallocate(p_Ihcol)
                      p_Ihcol => p_IhTmp

                      allocate (p_IhTmp(iblocks+NmemBlkCount))
                      p_IhTmp(1:iblocks) = p_Ihindx(1:iblocks)
                      deallocate(p_Ihindx)
                      p_Ihindx => p_IhTmp
                    
                      allocate (p_IsizeTmp(iblocks+NmemBlkCount))
                      p_IsizeTmp(1:iblocks) = p_Isize(1:iblocks)
                      deallocate(p_Isize)
                      p_Isize => p_IsizeTmp

                      allocate (RmemblockTmp(iblocks+NmemBlkCount))
                      RmemblockTmp(1:iblocks) = Rmemblock(1:iblocks)
                      deallocate(Rmemblock)
                      Rmemblock => RmemblockTmp
                      
                      ! Now we have enough blocks again.
                    end if

                    ! Add a new block

                    iblocks = iblocks + 1
                    p_Isize (iblocks) = imemblkSize
                    
                    ! Move the start index behind the last completely
                    ! occupied block
                             
                    istartIdx = iallocated
                    icurrentblock = iblocks
                  
                    ! Allocate a new memory block of size imemblkSize
                    !
                    ! Use p_Icol and p_Iindx - they are not used anymore.
                    ! Allocate another imemblkSize elements for column numbers and
                    ! list pointers.

                    call storage_new1D ('bilf_createMatStructure9_conf', 'Ihicol', &
                                        p_Isize (iblocks), ST_INT, p_Ihcol(iblocks), &
                                        ST_NEWBLOCK_NOINIT)
                    call storage_getbase_int (p_Ihcol(iblocks),p_Icol)

                    ! The new index array must be filled with 0 - otherwise
                    ! the search routine below won't work!
                    call storage_new1D ('bilf_createMatStructure9_conf', 'p_Ihindx', &
                                        p_Isize (iblocks), ST_INT, p_Ihindx(iblocks), &
                                        ST_NEWBLOCK_ZERO)
                    call storage_getbase_int (p_Ihindx(iblocks),p_Iindx)
                    
                    Rmemblock(iblocks)%p_Icol => p_Icol
                    Rmemblock(iblocks)%p_Iindx => p_Iindx

                    iallocated = iallocated + p_Isize (iblocks)
                    
                  else
                    ! Be careful when leaving the current memory block
                    ! for insertion of an element at position NA!
                    !
                    ! If the new position is not in the current block,
                    ! it's in the last block... and so set the pointer
                    ! and indices appropriately!
                    
                    if ( NA .gt. (istartidx+p_Isize (icurrentblock))) then
                      istartidx = iallocated-p_Isize(iblocks)
                      icurrentblock = iblocks
                      p_Icol => Rmemblock(iblocks)%p_Icol
                      p_Iindx => Rmemblock(iblocks)%p_Iindx
                    end if
                    
                  end if
                
                  ! Append JCOL to p_Icol
                  p_Icol(NA-istartIdx) = JCOL
                  
                  ! We have to make sure that p_Indx(NA)=0 to indicate the end of
                  ! the list. Ok, this is trivial because we allocated it with 
                  ! storage_new, configured to fill the memory with 0, so it is 0.
                  !
                  ! The searchloop ends here, continue with next JDOFE
                  
                  exit 
                
                else
                
                  ! No, the list does not end here.
                  ! Take the next element in the list
                  IPOS = p_Iindx(IPOS-istartidx)
                  
                  ! Be careful when leaving the current memory block
                  do while ( IPOS .gt. (istartidx+p_Isize (icurrentblock)) )
                  
                    ! go to the next memory block and search there
                    istartidx = istartidx+p_Isize(icurrentblock)
                    icurrentblock = icurrentblock+1
                    p_Icol => Rmemblock(icurrentblock)%p_Icol
                    p_Iindx => Rmemblock(icurrentblock)%p_Iindx
                    
                  end do ! IPOS .GT. (istartidx+p_Isize (iblocks))
                
                end if ! p_Iindx(IPOS) = 0
                
              end do searchloop
            
            end if ! p_Icol(IROW) = 0
                
          end do ! JDOFE
         
        end do ! IELneighIdxJ
      
      end do ! IDOFE
        
    end do ! IEL
  
  end do ! IELset

  ! Clean up the DOF's arrays    
  deallocate(IdofsTest)
  deallocate(IdofsTrial)

  ! --------------------- DOF COLLECTION PHASE ------------------------
    
  ! Ok, p_Icol is built. The hardest part is done!
  ! Now build KCOL by collecting the entries in the linear lists of 
  ! each row.
  !
  ! Save NA in the matrix structure
  rmatrixScalar%NA = NA
  
  ! As we now NA, we can allocate the real KCOL now!
  
  call storage_new1D ('bilf_createMatStructure9_conf', 'KCOL', &
                      NA, ST_INT, rmatrixScalar%h_KCOL, &
                      ST_NEWBLOCK_NOINIT)
  call storage_getbase_int (rmatrixScalar%h_Kcol,p_KCOL)
  
  ! Set back NA to 0 at first.

  NA=0
      
  ! Loop through all of the NEQ linear lists:
  
  do IEQ=1,NEQ
  
    ! We are at the head of the list, now we have to walk
    ! through it to append the entries to KCOL.
    ! We always start in the first memory block.

    if (icurrentblock .ne. 1) then
      icurrentblock = 1
      istartidx = 0
      p_Icol => Rmemblock(1)%p_Icol
      p_Iindx => Rmemblock(1)%p_Iindx
    end if
    
    ! Add the head of the list to KCOL:
  
    NA=NA+1
    p_KCOL(NA)=p_Icol(IEQ)

    ! Set KLD appropriately:

    p_KLD(IEQ) = NA
    
    IPOS = IEQ
    
    do while (p_Iindx(IPOS-istartidx).ne.0)
    
      ! Get the position of the next entry in p_Icol:
      IPOS = p_Iindx(IPOS-istartidx)
      
      ! Be careful when leaving the current memory block
      do while ( IPOS .gt. (istartidx+p_Isize (icurrentblock)) )
      
        ! go to the next memory block and search there
        istartidx = istartidx+p_Isize(icurrentblock)
        icurrentblock = icurrentblock+1
        p_Icol => Rmemblock(icurrentblock)%p_Icol
        p_Iindx => Rmemblock(icurrentblock)%p_Iindx
        
      end do ! IPOS .GT. (istartidx+p_Isize (iblocks))
      
      ! Add the column number to the row in KCOL:
      NA=NA+1
      p_KCOL(NA)=p_Icol(IPOS-istartidx)
    
    end do ! KINDX(IPOS) <> 0

  end do ! IEQ
  
  ! Append the final entry to KLD:
  
  p_KLD(NEQ+1)=NA+1
  
  ! Sort entries on KCOL separately for each row.
  ! This is a small bubble-sort...
  !
  ! Loop through all rows:
  nmaxCol = 0

  do IEQ=1,NEQ

    ! Repeat until everything is sorted.
    
    BSORT=.false.
    do while (.not. BSORT)
    
      BSORT=.true.

      !  Loop through the line 

      do JCOL=p_KLD(IEQ),p_KLD(IEQ+1)-2
      
        ! If the next element is larger...
      
        if (p_KCOL(JCOL) .gt. p_KCOL(JCOL+1)) then
        
          ! Change position of the current and next element
        
          IHELP=p_KCOL(JCOL)
          p_KCOL(JCOL)=p_KCOL(JCOL+1)
          p_KCOL(JCOL+1)=IHELP
          
          ! And repeat the sorting of that line
          
          BSORT=.false.
          
        end if
        
      end do ! JCOL
      
    end do ! (not BSORT)      

    ! Grab the diagonal
    do JCOL=p_KLD(IEQ),p_KLD(IEQ+1)-1
      if (p_KCOL(JCOL) .ge. IEQ) then
        p_Kdiagonal(IEQ) = JCOL
        exit
      end if
    end do   
    
    ! Grab the largest column number. As the current line is sorted,
    ! we can find this using the end of the line.
    nmaxCol = max(nmaxCol,p_Kcol(p_Kld(IEQ+1)-1))

  end do ! IEQ
  
  ! HOORAY, THAT'S IT!
  ! Deallocate all temporary memory...
  
  do i=iblocks,1,-1
    call storage_free(p_Ihcol(i))
    call storage_free(p_Ihindx(i))
  end do
  
  deallocate(IadjElem)
  deallocate(IadjPtr)
  deallocate(Rmemblock)
  deallocate(p_Isize)
  deallocate(p_Ihindx)
  deallocate(p_Ihcol)
  
  end subroutine

!  !****************************************************************************
!
!!<subroutine>
!
!  SUBROUTINE bilf_buildMatrix9d_conf (rdiscretisation,rform,bclear,rmatrixScalar,&
!                                      fcoeff_buildMatrixSc_sim,rcollection)
!  
!!<description>
!  ! This routine calculates the entries of a finite element matrix.
!  ! The matrix structure must be prepared with bilf_createMatrixStructure
!  ! in advance. The discretisation is assumed to be conformal, i.e. the DOF's
!  ! of all finite elements must 'match'. Trial and test functions may be
!  ! different.
!  ! In case the array for the matrix entries does not exist, the routine
!  ! allocates memory in size of the matrix of the heap for the matrix entries.
!  !
!  ! Double-precision version.
!  !
!  ! WARNING: THIS ROUTINE IS OUTDATED AND NOT USED. ONLY FOR TESTING AND
!  ! DEMONSTRATION PURPOSES! IT WEILL BE REMOVED SOONER OR LATER!
!!</description>
!
!!<input>
!  ! The underlying discretisation structure which is to be used to
!  ! create the matrix.
!  TYPE(t_spatialDiscretisation), INTENT(IN), TARGET :: rdiscretisation
!  
!  ! The bilinear form specifying the underlying PDE of the discretisation.
!  TYPE(t_bilinearForm), INTENT(IN) :: rform
!  
!  ! Whether to clear the matrix before calculating the entries.
!  ! If .FALSE., the new matrix entries are added to the existing entries.
!  LOGICAL, INTENT(IN) :: bclear
!  
!  ! OPTIONAL: A pointer to a collection structure. This structure is given to the
!  ! callback function for nonconstant coefficients to provide additional
!  ! information. 
!  TYPE(t_collection), INTENT(INOUT), TARGET, OPTIONAL :: rcollection
!  
!  ! OPTIONAL: A callback routine for nonconstant coefficient matrices.
!  ! Must be present if the matrix has nonconstant coefficients!
!  INCLUDE 'intf_coefficientMatrixSc2.inc'
!  OPTIONAL :: fcoeff_buildMatrixSc_sim
!!</input>
!
!!<inputoutput>
!  ! The FE matrix. Calculated matrix entries are imposed to this matrix.
!  TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrixScalar
!!</inputoutput>
!
!!</subroutine>
!
!  ! local variables
!  INTEGER :: i,i1,j,icurrentElementDistr,IDFG, ICUBP, IALBET, IA, IB, NVE
!  LOGICAL :: bIdenticalTrialAndTest, bnonparTest, bnonparTrial
!  INTEGER(I32) :: IEL, IELmax, IELset, IDOFE, JDOFE
!  INTEGER :: JCOL0,JCOL
!  REAL(DP) :: OM,AUX, DB
!  
!  ! Array to tell the element which derivatives to calculate
!  LOGICAL, DIMENSION(EL_MAXNDER) :: BderTrialTempl, BderTestTempl, BderTrial, BderTest
!  
!  ! Cubature point coordinates on the reference element
!  REAL(DP), DIMENSION(CUB_MAXCUBP, NDIM3D) :: Dxi
!
!  ! For every cubature point on the reference element,
!  ! the corresponding cubature weight
!  REAL(DP), DIMENSION(CUB_MAXCUBP) :: Domega
!  
!  ! number of cubature points on the reference element
!  INTEGER :: ncubp
!  
!  ! Pointer to KLD, KCOL, DA
!  INTEGER(I32), DIMENSION(:), POINTER :: p_KLD, p_KCOL
!  REAL(DP), DIMENSION(:), POINTER :: p_DA
!  
!  ! An allocateable array accepting the DOF's of a set of elements.
!  INTEGER, DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsTest, IdofsTrial
!  INTEGER, DIMENSION(:,:), POINTER :: p_IdofsTrial
!  !INTEGER, DIMENSION(EL_MAXNBAS,BILF_NELEMSIM), TARGET :: IdofsTest, IdofsTrial
!  !INTEGER, DIMENSION(:,:), POINTER :: p_IdofsTrial
!  
!  ! Allocateable arrays for the values of the basis functions - 
!  ! for test and trial spaces.
!  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: DbasTest,DbasTrial
!  REAL(DP), DIMENSION(:,:,:,:), POINTER :: p_DbasTrial
!  
!  ! Number of entries in the matrix - for quicker access
!  INTEGER :: NA
!  INTEGER(I32) :: NEQ
!  
!  ! Number of local degees of freedom for trial and test functions
!  INTEGER :: indofTrial, indofTest
!  
!  ! The triangulation structure - to shorten some things...
!  TYPE(t_triangulation), POINTER :: p_rtriangulation
!  
!  ! A pointer to an element-number list
!  INTEGER(I32), DIMENSION(:), POINTER :: p_IelementList
!  
!  ! Local matrices, used during the assembly.
!  ! Values and positions of values in the global matrix.
!  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: Kentry
!  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Dentry
!  
!  ! An array receiving the coordinates of cubature points on
!  ! the reference element for all elements in a set.
!  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: DcubPtsRef
!
!  ! An array receiving the coordinates of cubature points on
!  ! the real element for all elements in a set.
!  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: DcubPtsReal
!
!  ! Pointer to the point coordinates to pass to the element function.
!  ! Point either to DcubPtsRef or to DcubPtsReal, depending on whether
!  ! the trial/test element is parametric or not.
!  REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPtsTrial
!  REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPtsTest
!  
!  ! Array with coordinates of the corners that form the real element.
!  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Dcoords
!  
!  ! Arrays for saving Jacobian determinants and matrices
!  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Ddetj
!  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Djac
!  
!  ! Pointer to KVERT of the triangulation
!  INTEGER(I32), DIMENSION(:,:), POINTER :: p_IverticesAtElement
!  
!  ! Pointer to DCORVG of the triangulation
!  REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoords
!  
!  ! Current element distribution
!  TYPE(t_elementDistribution), POINTER :: p_relementDistribution
!  
!  ! Number of elements in a block. Normally =BILF_NELEMSIM,
!  ! except if there are less elements in the discretisation.
!  INTEGER :: nelementsPerBlock
!  
!  ! Some variables to support nonconstant coefficients in the matrix.
!  
!  ! Pointer to the coefficients that are computed by the callback routine.
!  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Dcoefficients
!  
!  !REAL(DP), DIMENSION(11) :: DT
!  
!  !CHARACTER(LEN=20) :: CFILE
!  
!  IF (.NOT. ASSOCIATED(rmatrixScalar%p_rspatialDiscretisation)) THEN
!    CALL output_line ('No discretisation associated!', &
!        OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrix9d_conf')
!    CALL sys_halt()
!  END IF
!
!  ! Which derivatives of basis functions are needed?
!  ! Check the descriptors of the bilinear form and set BDERxxxx
!  ! according to these.
!
!  !CALL ZTIME(DT(1))
!
!  BderTrialTempl = .FALSE.
!  BderTestTempl = .FALSE.
!  
!  ! Loop through the additive terms
!  DO i=1,rform%itermCount
!    ! The desriptor Idescriptors gives directly the derivative
!    ! which is to be computed! Build template's for BDER.
!    ! We don't compute the actual BDER here, as there might be some special
!    ! processing if trial/test functions are identical!
!    !
!    ! At first build the descriptors for the trial functions
!    I1=rform%Idescriptors(1,I)
!    
!    IF ((I1 .LE.0) .OR. (I1 .GT. DER_MAXNDER)) THEN
!      CALL output_line ('Invalid descriptor!', &
!          OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrix9d_conf')
!      CALL sys_halt()
!    ENDIF
!    
!    BderTrialTempl(I1)=.TRUE.
!
!    ! Then those of the test functions
!    I1=rform%Idescriptors(2,I)
!    
!    IF ((I1 .LE.0) .OR. (I1 .GT. DER_MAXNDER)) THEN
!      CALL output_line ('Invalid descriptor!', &
!          OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrix9d_conf')
!      CALL sys_halt()
!    ENDIF
!    
!    BderTestTempl(I1)=.TRUE.
!  END DO
!  
!  ! Get information about the matrix:
!  NA = rmatrixScalar%NA
!  NEQ = rmatrixScalar%NEQ
!  
!  ! We need KCOL/KLD of our matrix
!  CALL lsyssc_getbase_Kcol (rmatrixScalar,p_KCOL)
!  CALL lsyssc_getbase_Kld (rmatrixScalar,p_KLD)
!  
!  ! Check if the matrix entries exist. If not, allocate the matrix.
!  IF (rmatrixScalar%h_DA .EQ. ST_NOHANDLE) THEN
!
!    ! Clear the entries in the matrix - we need to start with zero
!    ! when assembling a new matrix!
!    CALL storage_new1D ('bilf_buildMatrix9d_conf', 'DA', &
!                        NA, ST_DOUBLE, rmatrixScalar%h_DA, &
!                        ST_NEWBLOCK_ZERO)
!    CALL lsyssc_getbase_double (rmatrixScalar,p_DA)
!
!  ELSE
!  
!    CALL lsyssc_getbase_double (rmatrixScalar,p_DA)
!
!    ! If desired, clear the matrix before assembling.
!    IF (bclear) THEN
!      CALL lalg_clearVectorDble (p_DA)
!    END IF
!    
!  END IF
!  
!  ! Get a pointer to the triangulation - for easier access.
!  p_rtriangulation => rdiscretisation%p_rtriangulation
!  
!  ! For saving some memory in smaller discretisations, we calculate
!  ! the number of elements per block. For smaller triangulations,
!  ! this is NEL. If there are too many elements, it's at most
!  ! BILF_NELEMSIM. This is only used for allocating some arrays.
!  nelementsPerBlock = MIN(BILF_NELEMSIM,p_rtriangulation%NEL)
!  
!  ! Get a pointer to the KVERT and DCORVG array
!  CALL storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement, &
!                             p_IverticesAtElement)
!  CALL storage_getbase_double2D(p_rtriangulation%h_DvertexCoords, &
!                             p_DvertexCoords)
!
!  ! Allocate memory for corner coordinates
!  ALLOCATE(DCoords(2,TRIA_MAXNVE2D,nelementsPerBlock))
!  
!  ! Now loop over the different element distributions (=combinations
!  ! of trial and test functions) in the discretisation.
!  !CALL ZTIME(DT(2))
!
!  DO icurrentElementDistr = 1,rdiscretisation%inumFESpaces
!  
!    ! Activate the current element distribution
!    p_relementDistribution => rdiscretisation%RelementDistr(icurrentElementDistr)
!  
!    ! Cancel if this element distribution is empty.
!    IF (p_relementDistribution%NEL .EQ. 0) CYCLE
!
!    ! Get the number of local DOF's for trial and test functions
!    indofTrial = elem_igetNDofLoc(p_relementDistribution%itrialElement)
!    indofTest = elem_igetNDofLoc(p_relementDistribution%itestElement)
!    
!    ! Get the number of corner vertices of the element
!    NVE = elem_igetNVE(p_relementDistribution%itrialElement)
!    IF (NVE .NE. elem_igetNVE(p_relementDistribution%itestElement)) THEN
!      CALL output_line ('Element spaces incompatible!', &
!          OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrix9d_conf')
!      CALL sys_halt()
!    END IF
!    
!    ! Initialise the cubature formula,
!    ! Get cubature weights and point coordinates on the reference element
!    CALL cub_getCubPoints(p_relementDistribution%ccubTypeBilForm, ncubp, Dxi, Domega)
!    
!    ! Allocate arrays accepting cubature point coordinates.
!    ! It's at most as large as number of elements or length
!    ! of the element set.
!    ALLOCATE(DcubPtsRef(p_rtriangulation%ndim,ncubp,nelementsPerBlock))
!    ALLOCATE(DcubPtsReal(p_rtriangulation%ndim,ncubp,nelementsPerBlock))
!    
!    ! Put the cubature point coordinates in the right format to the
!    ! cubature-point array.
!    ! Initialise all entries in DcubPtsRef with the same coordinates -
!    ! as the cubature point coordinates are identical on all elements
!    DO j=1,SIZE(DcubPtsRef,3)
!      DO i=1,ncubp
!        DcubPtsRef(1,i,j) = Dxi(i,1)
!        DcubPtsRef(2,i,j) = Dxi(i,2)
!      END DO
!    END DO
!    
!    ! Allocate an array saving the coordinates of corner vertices of elements
!    ALLOCATE(Djac(p_rtriangulation%ndim**2,ncubp,nelementsPerBlock))
!    ALLOCATE(Ddetj(ncubp,nelementsPerBlock))
!    
!    ! Allocate arrays for the values of the test- and trial functions.
!    ! This is done here in the size we need it. Allocating it in-advance
!    ! with something like
!    !  ALLOCATE(DbasTest(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
!    !  ALLOCATE(DbasTrial(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
!    ! would lead to nonused memory blocks in these arrays during the assembly, 
!    ! which reduces the speed by 50%!
!    
!    ALLOCATE(DbasTest(indofTest,elem_getMaxDerivative(p_relementDistribution%itestElement),&
!             ncubp,nelementsPerBlock))
!    ALLOCATE(DbasTrial(indofTrial,elem_getMaxDerivative(p_relementDistribution%itrialElement), &
!             ncubp,nelementsPerBlock))
!
!    ! Allocate memory for the DOF's of all the elements.
!    ALLOCATE(IdofsTest(indofTest,nelementsPerBlock))
!    ALLOCATE(IdofsTrial(indofTrial,nelementsPerBlock))
!
!    ! Check if one of the trial/test elements is nonparametric
!    bnonparTrial = elem_isNonparametric(p_relementDistribution%itrialElement)
!    bnonparTest  = elem_isNonparametric(p_relementDistribution%itestElement)
!                    
!    ! Let p_DcubPtsTrial / p_DcubPtsTest point either to DcubPtsReal or
!    ! DcubPtsRef - depending on whether the space is parametric or not.
!    IF (bnonparTrial) THEN
!      p_DcubPtsTrial => DcubPtsReal
!    ELSE
!      p_DcubPtsTrial => DcubPtsRef
!    END IF
!    
!    IF (bnonparTest) THEN
!      p_DcubPtsTest => DcubPtsReal
!    ELSE
!      p_DcubPtsTest => DcubPtsRef
!    END IF
!    
!    ! Allocate an array saving the local matrices for all elements
!    ! in an element set.
!    ! We could also allocate EL_MAXNBAS*EL_MAXNBAS*BILF_NELEMSIM integers
!    ! for this local matrix, but this would normally not fit to the cache
!    ! anymore! indofTrial*indofTest*BILF_NELEMSIM is normally much smaller!
!    ALLOCATE(Kentry(indofTest,indofTrial,nelementsPerBlock))
!    ALLOCATE(Dentry(indofTest,indofTrial))
!    
!    ! In case of nonconstant coefficients in that part of the matrix, we
!    ! need an additional array to save all the coefficients:
!    IF (.NOT. rform%BconstantCoeff(icurrentElementDistr)) THEN
!      IF (rform%ballCoeffConstant) THEN
!        CALL output_line ('Some oefficients are not constant &
!            &although thy should be!', &
!            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrix9d_conf')
!        CALL sys_halt()
!      END IF
!      IF (.NOT. PRESENT(fcoeff_buildMatrixSc_sim)) THEN
!        CALL output_line ('Coefficient function not given!',&
!            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrix9d_conf')
!        CALL sys_halt()
!      END IF
!      ALLOCATE(Dcoefficients(rform%itermCount,ncubp,nelementsPerBlock))
!    END IF
!                    
!    ! p_IdofsTest points either to the just created array or to the
!    ! array with the DOF's of the trial functions - when trial and
!    ! test functions are identical.
!    ! We don't rely on bidenticalTrialAndTest purely, as this does not
!    ! indicate whether there are identical trial and test functions
!    ! in one block!
!    bIdenticalTrialAndTest = p_relementDistribution%itrialElement .EQ. &
!                             p_relementDistribution%itestElement
!
!    ! Let p_IdofsTrial point either to IdofsTrial or to the DOF's of the test
!    ! space IdofTest (if both spaces are identical). 
!    ! We create a pointer for the trial space and not for the test space to
!    ! prevent pointer-arithmetic in the innerst loop below!
!    IF (bIdenticalTrialAndTest) THEN
!      p_IdofsTrial => IdofsTest
!      p_DbasTrial  => DbasTest
!      ! Build the actual combination of what the element should calculate.
!      ! As we evaluate only once, what the element must calculate is an
!      ! OR combination of the BDER from trial and test functions.
!      BderTrial = BderTrialTempl .OR. BderTestTempl
!      BderTest = BderTestTempl
!    ELSE
!      p_IdofsTrial => IdofsTrial
!      p_DbasTrial  => DbasTrial
!      
!      ! Build the actual combination of what the element should calculate.
!      ! Copy BDERxxxx to BDERxxxxAct
!      BderTrial = BderTrialTempl
!      BderTest = BderTestTempl
!    END IF
!    !CALL ZTIME(DT(3))
!    ! p_IelementList must point to our set of elements in the discretisation
!    ! with that combination of trial/test functions
!    CALL storage_getbase_int (p_relementDistribution%h_IelementList, &
!                              p_IelementList)
!                              
!    ! Loop over the elements - blockwise.
!    DO IELset = 1, p_rtriangulation%NEL, BILF_NELEMSIM
!    
!      ! We always handle BILF_NELEMSIM elements simultaneously.
!      ! How many elements have we actually here?
!      ! Get the maximum element number, such that we handle at most BILF_NELEMSIM
!      ! elements simultaneously.
!      
!      IELmax = MIN(p_rtriangulation%NEL,IELset-1+BILF_NELEMSIM)
!    
!      ! The outstanding feature with finite elements is: A basis
!      ! function for a DOF on one element has common support only
!      ! with the DOF's on the same element! E.g. for Q1:
!      !
!      !        #. . .#. . .#. . .#
!      !        .     .     .     .
!      !        .  *  .  *  .  *  .
!      !        #-----O-----O. . .#
!      !        |     |     |     .
!      !        |     | IEL |  *  .
!      !        #-----X-----O. . .#
!      !        |     |     |     .
!      !        |     |     |  *  .
!      !        #-----#-----#. . .#
!      !
!      ! --> On element IEL, the basis function at "X" only interacts
!      !     with the basis functions in "O". Elements in the 
!      !     neighbourhood ("*") have no support, therefore we only have
!      !     to collect all "O" DOF's.
!      !
!      ! Calculate the global DOF's into IdofsTrial / IdofsTest.
!      !
!      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
!      ! global DOF's of our BILF_NELEMSIM elements simultaneously.
!      CALL dof_locGlobMapping_mult(rdiscretisation, p_IelementList(IELset:IELmax), &
!                                  .TRUE.,IdofsTest)
!                                   
!      ! If the DOF's for the trial functions are different, calculate them, too.
!      IF (.NOT.bIdenticalTrialAndTest) THEN
!        CALL dof_locGlobMapping_mult(rdiscretisation, p_IelementList(IELset:IELmax), &
!                                    .FALSE.,IdofsTrial)
!      END IF
!      !CALL ZTIME(DT(4))
!      ! For the assembly of the global matrix, we use a "local"
!      ! approach. At first we build a "local" system matrix according
!      ! to the current element. This contains all additive
!      ! contributions of element IEL, which are later added at the
!      ! right positions to the elements in the global system matrix.
!      !
!      ! We have indofTrial trial DOF's per element and
!      ! indofTest test DOF's per element. Therefore there are
!      ! indofTrial*indofTest tupel of basis-/testfunctions (phi_i,psi_j) 
!      ! "active" (i.e. have common support) on our current element, each 
!      ! giving an additive contribution to the system matrix.
!      !
!      ! We build a quadratic indofTrial*indofTest local matrix:
!      ! Kentry(1..indofTest,1..indofTrial) receives the position 
!      !   in the global system matrix, where the corresponding value 
!      !   has to be added to.
!      ! (The corresponding contrbutions can be saved separately, 
!      !  but we directly add them to the global matrix in this 
!      !  approach.)
!      !
!      ! We build local matrices for all our elements 
!      ! in the set simultaneously.
!      ! Loop through elements in the set and for each element,
!      ! loop through the local matrices to initialise them:
!      DO IEL=1,IELmax-IELset+1
!
!        ! Loop through the trial functions
!        DO JDOFE=1,indofTrial
!        
!          ! Row JDOFE of the local matrix corresponds 
!          ! to row=global DOF KDFG(JDOFE) in the global matrix.
!          ! This is the "X" in the above picture.
!          ! Get the starting position of the corresponding row
!          ! to JCOL0:
!
!          JCOL0=p_KLD(p_IdofsTrial(JDOFE,IEL))
!          
!          ! Now we loop through the other DOF's on the current element
!          ! (the "O"'s).
!          ! All these have common support with our current basis function
!          ! and will therefore give an additive value to the global
!          ! matrix.
!          
!          DO IDOFE=1,indofTest
!            
!            ! Get the global DOF of the "O" which interacts with 
!            ! our "X".
!            
!            IDFG=IdofsTest(IDOFE,IEL)
!            
!            ! Starting in JCOL0 (which points to the beginning of
!            ! the line initially), loop through the elements in
!            ! the row to find the position of column IDFG.
!            ! Jump out of the DO loop if we find the column.
!            
!            DO JCOL=JCOL0,NA
!              IF (p_KCOL(JCOL) .EQ. IDFG) EXIT
!            END DO
!
!            ! Because columns in the global matrix are sorted 
!            ! ascendingly (except for the diagonal element),
!            ! the next search can start after the column we just found.
!            
!            ! JCOL0=JCOL+1
!            
!            ! Save the position of the matrix entry into the local
!            ! matrix.
!            
!            Kentry(IDOFE,JDOFE,IEL)=JCOL
!            
!          END DO ! IDOFE
!          
!        END DO ! JDOFE
!        
!      END DO ! IEL
!      !CALL ZTIME(DT(5))
!      ! Ok, we found the positions of the local matrix entries
!      ! that we have to change.
!      ! To calculate the matrix contributions, we have to evaluate
!      ! the elements to give us the values of the basis functions
!      ! in all the DOF's in all the elements in our set.
!      !
!      ! We have the coordinates of the cubature points saved in the
!      ! coordinate array from above. Unfortunately for nonparametric
!      ! elements, we need the real coordinate.
!      ! Furthermore, we anyway need the coordinates of the element
!      ! corners and the Jacobian determinants corresponding to
!      ! all the points.
!      !
!      ! At first, get the coordinates of the corners of all the
!      ! elements in the current set. 
!      
!!      DO IEL=1,IELmax-IELset+1
!!        DCoords(:,:,IEL) = p_DvertexCoords(:, &
!!                            p_IverticesAtElement(:,p_IelementList(IELset+IEL-1)))
!!      END DO
!      DO IEL=1,IELmax-IELset+1
!        DO J = 1,NVE
!          DO I = 1,p_rtriangulation%ndim
!            DCoords(I,J,IEL) = p_DvertexCoords(I, &
!                               p_IverticesAtElement(J,p_IelementList(IELset+IEL-1)))
!          END DO
!        END DO
!      END DO
!      !CALL ZTIME(DT(6))
!      
!      ! Depending on the type of transformation, we must now choose
!      ! the mapping between the reference and the real element.
!      ! In case we use a nonparametric element or a nonconstant coefficient function,
!      ! we need the coordinates of the points on the real element, too.
!      IF (bnonparTrial .OR. bnonparTest .OR. (.NOT. rform%ballCoeffConstant)) THEN
!      
!        CALL trafo_calctrafo_sim (&
!             rdiscretisation%RelementDistr(icurrentElementDistr)%ctrafoType,&
!             IELmax-IELset+1,ncubp,Dcoords,&
!             DcubPtsRef,Djac(:,:,1:IELmax-IELset+1),Ddetj(:,1:IELmax-IELset+1),DcubPtsReal)
!      
!      ELSE
!      
!        CALL trafo_calctrafo_sim (p_relementDistribution%ctrafoType,&
!             IELmax-IELset+1,ncubp,Dcoords,&
!             DcubPtsRef,Djac(:,:,1:IELmax-IELset+1),Ddetj(:,1:IELmax-IELset+1))
!             
!      END IF
!      
!      !CALL ZTIME(DT(7))
!      
!      ! If the matrix has nonconstant coefficients, calculate the coefficients now.
!      IF (.NOT. rform%ballCoeffConstant) THEN
!        CALL fcoeff_buildMatrixSc_sim (rdiscretisation,icurrentElementDistr, rform, &
!                  IELset,IELmax-IELset+1,ncubp,p_IelementList(IELset:IELmax),Dcoords, &
!                  DcubPtsRef,DcubPtsReal,p_IdofsTrial,IdofsTest,Djac,Ddetj, &
!                  Dcoefficients,rcollection)
!      END IF
!      
!      !CALL ZTIME(DT(8))                              
!      ! Calculate the values of the basis functions.
!      ! Pass p_DcubPts as point coordinates, which point either to the
!      ! coordinates on the reference element (the same for all elements)
!      ! or on the real element - depending on whether this is a 
!      ! parametric or nonparametric element.
!      CALL elem_generic_sim (p_relementDistribution%itestElement, Dcoords, &
!            Djac(:,:,1:IELmax-IELset+1), Ddetj(:,1:IELmax-IELset+1), &
!            BderTest, DbasTest, ncubp, IELmax-IELset+1, p_DcubPtsTest,&
!            rintSubset%p_ItwistIndex)
!            
!      ! Omit the calculation of the trial function values if they
!      ! are identical to the test function values.
!      IF (.NOT. bidenticalTrialAndTest) THEN
!        CALL elem_generic_sim (p_relementDistribution%itrialElement, Dcoords, &
!            Djac(:,:,1:IELmax-IELset+1), Ddetj(:,1:IELmax-IELset+1), &
!            BderTrial, DbasTrial, ncubp, IELmax-IELset+1, p_DcubPtsTrial,&
!            rintSubset%p_ItwistIndex)
!      END IF
!      !CALL ZTIME(DT(9))
!      ! Values of all basis functions calculated. Now we can start 
!      ! to integrate!
!      !
!      ! We have two different versions for the integration - one
!      ! with constant coefficients and one with nonconstant coefficients.
!      !
!      ! Check the bilinear form which one to use:
!      
!      IF (rform%ballCoeffConstant) THEN
!      
!        ! Constant coefficients. The coefficients are to be found in
!        ! the Dcoefficients variable of the form.
!        !
!        ! Loop over the elements in the current set.
!
!        DO IEL=1,IELmax-IELset+1
!          
!          ! Clear the local matrix
!          Dentry = 0.0_DP
!          
!          ! Loop over all cubature points on the current element
!          DO ICUBP = 1, ncubp
!
!            ! calculate the current weighting factor in the cubature formula
!            ! in that cubature point.
!            !
!            ! Take the absolut value of the determinant of the mapping.
!            ! In 2D, the determinant is always positive, whereas in 3D,
!            ! the determinant might be negative -- that's normal!
!
!            OM = Domega(ICUBP)*ABS(Ddetj(ICUBP,IEL))
!
!            ! Loop over the additive factors in the bilinear form.
!            DO IALBET = 1,rform%itermcount
!            
!              ! Get from Idescriptors the type of the derivatives for the 
!              ! test and trial functions. The summand we calculate
!              ! here will be:
!              !
!              ! int_... ( phi_i )_IA  *  ( psi_j )_IB
!              !
!              ! -> Ix=0: function value, 
!              !      =1: first derivative, 
!              !      =2: 2nd derivative,...
!              !    as defined in the module 'derivative'.
!              
!              IA = rform%Idescriptors(1,IALBET)
!              IB = rform%Idescriptors(2,IALBET)
!              
!              ! Multiply OM with the coefficient of the form.
!              ! This gives the actual value to multiply the
!              ! function value with before summing up to the integral.
!              AUX = OM * rform%Dcoefficients(IALBET)
!            
!              ! Now loop through all possible combinations of DOF's
!              ! in the current cubature point. The outer loop
!              ! loops through the "X" in the above picture,
!              ! the trial functions:
!
!              DO JDOFE=1,indofTrial
!              
!                ! Get the value of the (trial) basis function 
!                ! phi_i (our "X") in the cubature point:
!                DB = p_DbasTrial(JDOFE,IA,ICUBP,IEL)
!                
!                !Perform an inner loop through the other DOF's
!                ! (the "X"'s). 
!
!                DO IDOFE=1,indofTest
!                
!                  ! Get the value of the basis function 
!                  ! psi_j (our "O") in the cubature point. 
!                  ! Them multiply:
!                  !    DB * DBAS(..) * AUX
!                  ! ~= phi_i * psi_j * coefficient * cub.weight
!                  ! Summing this up gives the integral, so the contribution
!                  ! to the global matrix. 
!                  !
!                  ! Simply summing up DB * DBAS(..) * AUX would give
!                  ! the coefficient of the local matrix. We save this
!                  ! contriobution in the local matrix.
!
!                  !JCOLB = Kentry(IDOFE,JDOFE,IEL)
!                  !p_DA(JCOLB) = p_DA(JCOLB) + DB*DbasTest(IDOFE,IB,ICUBP,IEL)*AUX
!                  Dentry(IDOFE,JDOFE) = Dentry(IDOFE,JDOFE)+DB*DbasTest(IDOFE,IB,ICUBP,IEL)*AUX
!                
!                END DO
!              
!              END DO ! JDOFE
!              
!            END DO ! IALBET
!
!          END DO ! ICUBP 
!          
!          ! Incorporate the local matrices into the global one.
!          ! Kentry gives the position of the additive contributions in Dentry.
!          DO JDOFE=1,indofTrial
!            DO IDOFE=1,indofTest
!              p_DA(Kentry(IDOFE,JDOFE,IEL)) = p_DA(Kentry(IDOFE,JDOFE,IEL)) + Dentry(IDOFE,JDOFE)
!            END DO
!          END DO
!
!        END DO ! IEL
!        
!      ELSE
!      
!        ! Nonconstant coefficients. The coefficients are to be found in
!        ! the Dcoefficients variable as computed above.
!        !
!        ! Loop over the elements in the current set.
!
!        DO IEL=1,IELmax-IELset+1
!          
!          ! Clear the local matrix
!          Dentry = 0.0_DP
!          
!          ! Loop over all cubature points on the current element
!          DO ICUBP = 1, ncubp
!
!            ! calculate the current weighting factor in the cubature formula
!            ! in that cubature point.
!            !
!            ! Take the absolut value of the determinant of the mapping.
!            ! In 2D, the determinant is always positive, whereas in 3D,
!            ! the determinant might be negative -- that's normal!
!
!            OM = Domega(ICUBP)*ABS(Ddetj(ICUBP,IEL))
!
!            ! Loop over the additive factors in the bilinear form.
!            DO IALBET = 1,rform%itermcount
!            
!              ! Get from Idescriptors the type of the derivatives for the 
!              ! test and trial functions. The summand we calculate
!              ! here will be:
!              !
!              ! int_... ( phi_i )_IA  *  ( psi_j )_IB
!              !
!              ! -> Ix=0: function value, 
!              !      =1: first derivative, 
!              !      =2: 2nd derivative,...
!              !    as defined in the module 'derivative'.
!              
!              IA = rform%Idescriptors(1,IALBET)
!              IB = rform%Idescriptors(2,IALBET)
!              
!              ! Multiply OM with the coefficient of the form.
!              ! This gives the actual value to multiply the
!              ! function value with before summing up to the integral.
!              ! Get the precalculated coefficient from the coefficient array.
!              AUX = OM * Dcoefficients(IALBET,ICUBP,IEL)
!            
!              ! Now loop through all possible combinations of DOF's
!              ! in the current cubature point. The outer loop
!              ! loops through the "X" in the above picture,
!              ! the trial functions:
!
!              DO JDOFE=1,indofTrial
!              
!                ! Get the value of the (trial) basis function 
!                ! phi_i (our "X") in the cubature point:
!                DB = p_DbasTrial(JDOFE,IA,ICUBP,IEL)
!                
!                !Perform an inner loop through the other DOF's
!                ! (the "O"'s). 
!
!                DO IDOFE=1,indofTest
!                
!                  ! Get the value of the basis function 
!                  ! psi_j (our "O") in the cubature point. This is
!                  ! DBAS(KDFL(IDOFE),IB).
!                  ! Them multiply:
!                  !    DB * DBAS(..) * AUX
!                  ! ~= phi_i * psi_j * coefficient * cub.weight
!                  ! Summing this up gives the integral, so the contribution
!                  ! to the global matrix. 
!                  !
!                  ! Simply summing up DB * DBAS(..) * AUX would give
!                  ! the coefficient of the local matrix. We save this
!                  ! contriobution in the local matrix of element IEL.
!
!                  !JCOLB = Kentry(IDOFE,JDOFE,IEL)
!                  !p_DA(JCOLB) = p_DA(JCOLB) + DB*DbasTest(IDOFE,IB,ICUBP,IEL)*AUX
!                  Dentry(IDOFE,JDOFE) = Dentry(IDOFE,JDOFE)+DB*DbasTest(IDOFE,IB,ICUBP,IEL)*AUX
!                
!                END DO
!              
!              END DO ! JDOFE
!              
!            END DO ! IALBET
!
!          END DO ! ICUBP 
!          
!          ! Incorporate the local matrices into the global one.
!          ! Kentry gives the position of the additive contributions in Dentry.
!          DO JDOFE=1,indofTrial
!            DO IDOFE=1,indofTest
!              p_DA(Kentry(IDOFE,JDOFE,IEL)) = p_DA(Kentry(IDOFE,JDOFE,IEL)) + Dentry(IDOFE,JDOFE)
!            END DO
!          END DO
!
!        END DO ! IEL
!
!      END IF ! rform%ballCoeffConstant
!
!      !CALL ZTIME(DT(10))
!    END DO ! IELset
!    
!    IF (.NOT. rform%ballCoeffConstant) THEN
!      DEALLOCATE(Dcoefficients)
!    END IF
!    DEALLOCATE(IdofsTest)
!    DEALLOCATE(DbasTrial)
!    DEALLOCATE(DbasTest)
!    DEALLOCATE(Kentry)
!    DEALLOCATE(Dentry)
!    DEALLOCATE(Ddetj)
!    DEALLOCATE(Djac)
!    DEALLOCATE(DcubPtsReal)
!    DEALLOCATE(DcubPtsRef)
!
!  END DO ! icurrentElementDistr
!
!  ! Clean up memory, finish
!
!  DEALLOCATE(Dcoords)
!  !CALL ZTIME(DT(11))
!  
!  !DO i=2,11
!  !  PRINT *,'Time for assembly part ',i,': ',DT(i)-DT(i-1)
!  !END DO
!  
!  !CFILE = 'MATRIX2.TXT'
!  !CALL OWM17(p_DA,p_KCOL,p_KLD,&
!  !           NEQ,NEQ,.TRUE.,0,'MAT1  ',CFILE,'(D20.10)')
!
!  END SUBROUTINE

!  !****************************************************************************
!
!!<subroutine>
!
!  SUBROUTINE bilf_buildMatrix9d_conf2 (rform,bclear,rmatrixScalar,&
!                                       fcoeff_buildMatrixSc_sim,rcollection)
!  
!!<description>
!  ! This routine calculates the entries of a finite element matrix.
!  ! The matrix structure must be prepared with bilf_createMatrixStructure
!  ! in advance. The discretisation is assumed to be conformal, i.e. the DOF's
!  ! of all finite elements must 'match'. Trial and test functions may be
!  ! different.
!  ! In case the array for the matrix entries does not exist, the routine
!  ! allocates memory in size of the matrix of the heap for the matrix entries.
!  !
!  ! For setting up the entries, the discretisation structure attached to
!  ! the matrix is used (rmatrixScalar%p_rdiscretisation). This is
!  ! normally attached to the matrix by bilf_createMatrixStructure.
!  !
!  ! Double-precision version.
!!</description>
!
!!<input>
!  ! The bilinear form specifying the underlying PDE of the discretisation.
!  TYPE(t_bilinearForm), INTENT(IN) :: rform
!  
!  ! Whether to clear the matrix before calculating the entries.
!  ! If .FALSE., the new matrix entries are added to the existing entries.
!  LOGICAL, INTENT(IN) :: bclear
!  
!  ! OPTIONAL: A pointer to a collection structure. This structure is given to the
!  ! callback function for nonconstant coefficients to provide additional
!  ! information. 
!  TYPE(t_collection), INTENT(INOUT), TARGET, OPTIONAL :: rcollection
!  
!  ! OPTIONAL: A callback routine for nonconstant coefficient matrices.
!  ! Must be present if the matrix has nonconstant coefficients!
!  INCLUDE 'intf_coefficientMatrixSc.inc'
!  OPTIONAL :: fcoeff_buildMatrixSc_sim
!!</input>
!
!!<inputoutput>
!  ! The FE matrix. Calculated matrix entries are imposed to this matrix.
!  TYPE(t_matrixScalar), INTENT(INOUT) :: rmatrixScalar
!!</inputoutput>
!
!!</subroutine>
!
!  ! local variables
!  INTEGER :: i,i1,j,k,icurrentElementDistr,JDFG, ICUBP, IALBET, IA, IB
!  LOGICAL :: bIdenticalTrialAndTest, bnonparTest, bnonparTrial
!  INTEGER(I32) :: IEL, IELmax, IELset, IDOFE, JDOFE
!  INTEGER :: JCOL0,JCOL
!  REAL(DP) :: OM,AUX, DB
!  
!  ! Array to tell the element which derivatives to calculate
!  LOGICAL, DIMENSION(EL_MAXNDER) :: BderTrialTempl, BderTestTempl, BderTrial, BderTest
!  
!  ! Cubature point coordinates on the reference element
!  REAL(DP), DIMENSION(CUB_MAXCUBP, NDIM3D) :: Dxi
!
!  ! For every cubature point on the reference element,
!  ! the corresponding cubature weight
!  REAL(DP), DIMENSION(CUB_MAXCUBP) :: Domega
!  
!  ! number of cubature points on the reference element
!  INTEGER :: ncubp
!  
!  ! Pointer to KLD, KCOL, DA
!  INTEGER(I32), DIMENSION(:), POINTER :: p_KLD, p_KCOL
!  REAL(DP), DIMENSION(:), POINTER :: p_DA
!  
!  ! An allocateable array accepting the DOF's of a set of elements.
!  INTEGER, DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsTest, IdofsTrial
!  INTEGER, DIMENSION(:,:), POINTER :: p_IdofsTrial
!  !INTEGER, DIMENSION(EL_MAXNBAS,BILF_NELEMSIM), TARGET :: IdofsTest, IdofsTrial
!  !INTEGER, DIMENSION(:,:), POINTER :: p_IdofsTrial
!  
!  ! Allocateable arrays for the values of the basis functions - 
!  ! for test and trial spaces.
!  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: DbasTest,DbasTrial
!  REAL(DP), DIMENSION(:,:,:,:), POINTER :: p_DbasTrial
!  
!  ! Number of entries in the matrix - for quicker access
!  INTEGER :: NA,NVE
!  INTEGER(I32) :: NEQ
!  
!  ! Number of local degees of freedom for trial and test functions
!  INTEGER :: indofTrial, indofTest
!  
!  ! The triangulation structure - to shorten some things...
!  TYPE(t_triangulation), POINTER :: p_rtriangulation
!  
!  ! A pointer to an element-number list
!  INTEGER(I32), DIMENSION(:), POINTER :: p_IelementList
!  
!  ! Local matrices, used during the assembly.
!  ! Values and positions of values in the global matrix.
!  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: Kentry
!  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Dentry
!  
!  ! An array receiving the coordinates of cubature points on
!  ! the reference element for all elements in a set.
!  REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPtsRef
!
!  ! An array receiving the coordinates of cubature points on
!  ! the real element for all elements in a set.
!  REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPtsReal
!
!  ! Pointer to the point coordinates to pass to the element function.
!  ! Point either to p_DcubPtsRef or to p_DcubPtsReal, depending on whether
!  ! the trial/test element is parametric or not.
!  REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPtsTrial
!  REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPtsTest
!  
!  ! Array with coordinates of the corners that form the real element.
!  REAL(DP), DIMENSION(:,:,:), POINTER :: p_Dcoords
!  
!  ! Arrays for saving Jacobian determinants and matrices
!  REAL(DP), DIMENSION(:,:), POINTER :: p_Ddetj
!  REAL(DP), DIMENSION(:,:,:), POINTER :: p_Djac
!  
!  ! Pointer to KVERT of the triangulation
!  INTEGER(I32), DIMENSION(:,:), POINTER :: p_IverticesAtElement
!  
!  ! Pointer to DCORVG of the triangulation
!  REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoords
!  
!  ! Current element distribution
!  TYPE(t_elementDistribution), POINTER :: p_relementDistribution
!  
!  ! Number of elements in a block. Normally =BILF_NELEMSIM,
!  ! except if there are less elements in the discretisation.
!  INTEGER :: nelementsPerBlock
!  
!  ! Some variables to support nonconstant coefficients in the matrix.
!  
!  ! Pointer to the coefficients that are computed by the callback routine.
!  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Dcoefficients
!  
!  ! A t_domainIntSubset structure that is used for storing information
!  ! and passing it to callback routines.
!  TYPE(t_domainIntSubset) :: rintSubset
!  
!  ! The discretisation - for easier access
!  TYPE(t_spatialDiscretisation), POINTER :: p_rdiscretisation
!  
!  IF (.NOT. ASSOCIATED(rmatrixScalar%p_rspatialDiscretisation)) THEN
!    CALL output_line ('No discretisation associated!',&
!        OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrix9d_conf2')
!    CALL sys_halt()
!  END IF
!
!  ! Which derivatives of basis functions are needed?
!  ! Check the descriptors of the bilinear form and set BDERxxxx
!  ! according to these.
!
!  BderTrialTempl = .FALSE.
!  BderTestTempl = .FALSE.
!  
!  ! Loop through the additive terms
!  DO i=1,rform%itermCount
!    ! The desriptor Idescriptors gives directly the derivative
!    ! which is to be computed! Build template's for BDER.
!    ! We don't compute the actual BDER here, as there might be some special
!    ! processing if trial/test functions are identical!
!    !
!    ! At first build the descriptors for the trial functions
!    I1=rform%Idescriptors(1,I)
!    
!    IF ((I1 .LE.0) .OR. (I1 .GT. DER_MAXNDER)) THEN
!      CALL output_line ('Invalid descriptor!',&
!          OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrix9d_conf2')
!      CALL sys_halt()
!    ENDIF
!    
!    BderTrialTempl(I1)=.TRUE.
!
!    ! Then those of the test functions
!    I1=rform%Idescriptors(2,I)
!    
!    IF ((I1 .LE.0) .OR. (I1 .GT. DER_MAXNDER)) THEN
!      CALL output_line ('Invalid descriptor!',&
!          OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrix9d_conf2')
!      CALL sys_halt()
!    ENDIF
!    
!    BderTestTempl(I1)=.TRUE.
!  END DO
!  
!  ! Get information about the matrix:
!  NA = rmatrixScalar%NA
!  NEQ = rmatrixScalar%NEQ
!  
!  ! We need KCOL/KLD of our matrix
!  IF ((rmatrixScalar%h_KCOL .EQ. ST_NOHANDLE) .OR. &
!      (rmatrixScalar%h_KLD .EQ. ST_NOHANDLE)) THEN
!    CALL output_line ('No discretisation structure! Cannot assemble matrix!', &
!                      OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrix9d_conf2')
!    CALL sys_halt()
!  END IF
!  
!  CALL lsyssc_getbase_Kcol (rmatrixScalar,p_KCOL)
!  CALL lsyssc_getbase_Kld (rmatrixScalar,p_KLD)
!  
!  ! Check if the matrix entries exist. If not, allocate the matrix.
!  IF (rmatrixScalar%h_DA .EQ. ST_NOHANDLE) THEN
!
!    ! Clear the entries in the matrix - we need to start with zero
!    ! when assembling a new matrix!
!    CALL storage_new1D ('bilf_buildMatrix9d_conf', 'DA', &
!                        NA, ST_DOUBLE, rmatrixScalar%h_DA, &
!                        ST_NEWBLOCK_ZERO)
!    CALL lsyssc_getbase_double (rmatrixScalar,p_DA)
!
!  ELSE
!  
!    CALL lsyssc_getbase_double (rmatrixScalar,p_DA)
!
!    ! If desired, clear the matrix before assembling.
!    IF (bclear) THEN
!      CALL lalg_clearVectorDble (p_DA)
!    END IF
!    
!  END IF
!  
!  ! Get the discretisation
!  p_rdiscretisation => rmatrixScalar%p_rspatialDiscretisation
!  
!  IF (.NOT. ASSOCIATED(p_rdiscretisation)) THEN
!    CALL output_line ('No discretisation attached to the matrix!',&
!        OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrix9d_conf2')
!    CALL sys_halt()
!  END IF
!  
!  ! Get a pointer to the triangulation - for easier access.
!  p_rtriangulation => p_rdiscretisation%p_rtriangulation
!  
!  ! For saving some memory in smaller discretisations, we calculate
!  ! the number of elements per block. For smaller triangulations,
!  ! this is NEL. If there are too many elements, it's at most
!  ! BILF_NELEMSIM. This is only used for allocating some arrays.
!  nelementsPerBlock = MIN(BILF_NELEMSIM,p_rtriangulation%NEL)
!  
!  ! Get a pointer to the KVERT and DCORVG array
!  CALL storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement, &
!                             p_IverticesAtElement)
!  CALL storage_getbase_double2D(p_rtriangulation%h_DvertexCoords, &
!                             p_DvertexCoords)
!
!  ! Now loop over the different element distributions (=combinations
!  ! of trial and test functions) in the discretisation.
!
!  DO icurrentElementDistr = 1,p_rdiscretisation%inumFESpaces
!  
!    ! Activate the current element distribution
!    p_relementDistribution => p_rdiscretisation%RelementDistr(icurrentElementDistr)
!  
!    ! Cancel if this element distribution is empty.
!    IF (p_relementDistribution%NEL .EQ. 0) CYCLE
!
!    ! Get the number of local DOF's for trial and test functions
!    indofTrial = elem_igetNDofLoc(p_relementDistribution%itrialElement)
!    indofTest = elem_igetNDofLoc(p_relementDistribution%itestElement)
!    
!    ! Get the number of vertices of the element, specifying the transformation
!    ! form the reference to the real element.
!    NVE = elem_igetNVE(p_relementDistribution%itrialElement)
!    IF (NVE .NE. elem_igetNVE(p_relementDistribution%itestElement)) THEN
!      CALL output_line ('Element spaces incompatible!',&
!          OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrix9d_conf2')
!      CALL sys_halt()
!    END IF
!    
!    ! Initialise the cubature formula,
!    ! Get cubature weights and point coordinates on the reference element
!    CALL cub_getCubPoints(p_relementDistribution%ccubTypeBilForm, ncubp, Dxi, Domega)
!    
!    ! Open-MP-Extension: Open threads here.
!    ! "j" is declared as private; shared gave errors with the Intel compiler
!    ! in Windows!?!
!    ! Each thread will allocate its own local memory...
!    !
!    !%OMP PARALLEL PRIVATE(rintSubset, p_DcubPtsRef,p_DcubPtsReal, &
!    !%OMP   p_Djac,p_Ddetj,p_Dcoords, kentry, dentry,DbasTest,DbasTrial, &
!    !%OMP   IdofsTest,IdofsTrial,bnonparTrial,bnonparTest,p_DcubPtsTrial,&
!    !%OMP   p_DcubPtsTest,Dcoefficients, bIdenticalTrialandTest, p_IdofsTrial, &
!    !%OMP   p_DbasTrial, BderTrial, BderTest, j, ielmax,IEL, idofe,jdofe,jcol0, &
!    !%OMP   jcol,JDFG,ICUBP, IALBET,OM,ia,ib,aux, db)
!    
!    ! Get from the trial element space the type of coordinate system
!    ! that is used there:
!    j = elem_igetCoordSystem(p_relementDistribution%itrialElement)
!    
!    ! Allocate memory and get local references to it.
!    CALL domint_initIntegration (rintSubset,nelementsPerBlock,ncubp,&
!        j,p_rtriangulation%ndim,NVE,&
!        MAX(elem_getTwistIndexSize(p_relementDistribution%itrialElement),&
!            elem_getTwistIndexSize(p_relementDistribution%itestElement)))
!    p_DcubPtsRef =>  rintSubset%p_DcubPtsRef
!    p_DcubPtsReal => rintSubset%p_DcubPtsReal
!    p_Djac =>        rintSubset%p_Djac
!    p_Ddetj =>       rintSubset%p_Ddetj
!    p_Dcoords =>     rintSubset%p_DCoords
!    
!    ! Put the cubature point coordinates in the right format to the
!    ! cubature-point array.
!    ! Initialise all entries in p_DcubPtsRef with the same coordinates -
!    ! as the cubature point coordinates are identical on all elements
!    DO j=1,SIZE(p_DcubPtsRef,3)
!      DO i=1,ncubp
!        DO k=1,SIZE(p_DcubPtsRef,1)
!          p_DcubPtsRef(k,i,j) = Dxi(i,k)
!        END DO
!      END DO
!    END DO
!    
!    ! Quickly check if one of the specified derivatives is out of the allowed range:
!    DO IALBET = 1,rform%itermcount
!      IA = rform%Idescriptors(1,IALBET)
!      IB = rform%Idescriptors(2,IALBET)      
!      IF ((IA.LT.0) .OR. &
!          (IA .GT. elem_getMaxDerivative(p_relementDistribution%itrialElement))) THEN
!        CALL output_line ('Specified trial-derivative '//TRIM(sys_siL(IA,10))//&
!            ' not available!',&
!            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrix9d_conf2')
!        CALL sys_halt()
!      END IF
!
!      IF ((IB.LT.0) .OR. &
!          (IB .GT. elem_getMaxDerivative(p_relementDistribution%itestElement))) THEN
!        CALL output_line ('Specified test-derivative '//TRIM(sys_siL(IA,10))//&
!            ' not available!',&
!            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrix9d_conf2')
!        CALL sys_halt()
!      END IF
!    END DO
!    
!    ! Allocate an array saving the coordinates of corner vertices of elements
!    
!    ! Allocate arrays for the values of the test- and trial functions.
!    ! This is done here in the size we need it. Allocating it in-advance
!    ! with something like
!    !  ALLOCATE(DbasTest(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
!    !  ALLOCATE(DbasTrial(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
!    ! would lead to nonused memory blocks in these arrays during the assembly, 
!    ! which reduces the speed by 50%!
!    
!    ALLOCATE(DbasTest(indofTest,elem_getMaxDerivative(p_relementDistribution%itestElement),&
!             ncubp,nelementsPerBlock))
!    ALLOCATE(DbasTrial(indofTrial,elem_getMaxDerivative(p_relementDistribution%itrialElement), &
!             ncubp,nelementsPerBlock))
!
!    ! Allocate memory for the DOF's of all the elements.
!    ALLOCATE(IdofsTest(indofTest,nelementsPerBlock))
!    ALLOCATE(IdofsTrial(indofTrial,nelementsPerBlock))
!
!    ! Check if one of the trial/test elements is nonparametric
!    bnonparTrial = elem_isNonparametric(p_relementDistribution%itrialElement)
!    bnonparTest  = elem_isNonparametric(p_relementDistribution%itestElement)
!                    
!    ! Let p_DcubPtsTrial / p_DcubPtsTest point either to p_DcubPtsReal or
!    ! p_DcubPtsRef - depending on whether the space is parametric or not.
!    IF (bnonparTrial) THEN
!      p_DcubPtsTrial => p_DcubPtsReal
!    ELSE
!      p_DcubPtsTrial => p_DcubPtsRef
!    END IF
!    
!    IF (bnonparTest) THEN
!      p_DcubPtsTest => p_DcubPtsReal
!    ELSE
!      p_DcubPtsTest => p_DcubPtsRef
!    END IF
!    
!    ! Allocate an array saving the local matrices for all elements
!    ! in an element set.
!    ! We could also allocate EL_MAXNBAS*EL_MAXNBAS*BILF_NELEMSIM integers
!    ! for this local matrix, but this would normally not fit to the cache
!    ! anymore! indofTrial*indofTest*BILF_NELEMSIM is normally much smaller!
!    ALLOCATE(Kentry(indofTrial,indofTest,nelementsPerBlock))
!    ALLOCATE(Dentry(indofTrial,indofTest))
!    
!    ! In case of nonconstant coefficients in that part of the matrix, we
!    ! need an additional array to save all the coefficients:
!    IF (.NOT. rform%BconstantCoeff(icurrentElementDistr)) THEN
!      IF (rform%ballCoeffConstant) THEN
!        CALL output_line ('Some oefficients are not constant &
!                &although thy should be!',&
!            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrix9d_conf2')
!        CALL sys_halt()
!      END IF
!      IF (.NOT. PRESENT(fcoeff_buildMatrixSc_sim)) THEN
!        CALL output_line ('Coefficient function not given!',&
!            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrix9d_conf2')
!        CALL sys_halt()
!      END IF
!      ALLOCATE(Dcoefficients(rform%itermCount,ncubp,nelementsPerBlock))
!    END IF
!                    
!    ! p_IdofsTest points either to the just created array or to the
!    ! array with the DOF's of the trial functions - when trial and
!    ! test functions are identical.
!    ! We don't rely on bidenticalTrialAndTest purely, as this does not
!    ! indicate whether there are identical trial and test functions
!    ! in one block!
!    bIdenticalTrialAndTest = p_relementDistribution%itrialElement .EQ. &
!                             p_relementDistribution%itestElement
!
!    ! Let p_IdofsTrial point either to IdofsTrial or to the DOF's of the test
!    ! space IdofTest (if both spaces are identical). 
!    ! We create a pointer for the trial space and not for the test space to
!    ! prevent pointer-arithmetic in the innerst loop below!
!    IF (bIdenticalTrialAndTest) THEN
!      p_IdofsTrial => IdofsTest
!      p_DbasTrial  => DbasTest
!      ! Build the actual combination of what the element should calculate.
!      ! As we evaluate only once, what the element must calculate is an
!      ! OR combination of the BDER from trial and test functions.
!      BderTrial = BderTrialTempl .OR. BderTestTempl
!      BderTest = BderTrial
!    ELSE
!      p_IdofsTrial => IdofsTrial
!      p_DbasTrial  => DbasTrial
!      
!      ! Build the actual combination of what the element should calculate.
!      ! Copy BDERxxxx to BDERxxxxAct
!      BderTrial = BderTrialTempl
!      BderTest = BderTestTempl
!    END IF
!    
!    ! p_IelementList must point to our set of elements in the discretisation
!    ! with that combination of trial/test functions
!    CALL storage_getbase_int (p_relementDistribution%h_IelementList, &
!                              p_IelementList)
!                              
!    ! Loop over the elements - blockwise.
!    !
!    ! Open-MP-Extension: Each loop cycle is executed in a different thread,
!    ! so BILF_NELEMSIM local matrices are simultaneously calculated in the
!    ! inner loop(s).
!    ! The blocks have all the same size, so we can use static scheduling.
!    !
!    !%OMP do schedule(static,1)
!    DO IELset = 1, SIZE(p_IelementList), BILF_NELEMSIM
!    
!      ! We always handle BILF_NELEMSIM elements simultaneously.
!      ! How many elements have we actually here?
!      ! Get the maximum element number, such that we handle at most BILF_NELEMSIM
!      ! elements simultaneously.
!      
!      IELmax = MIN(SIZE(p_IelementList),IELset-1+BILF_NELEMSIM)
!    
!      ! --------------------- DOF SEARCH PHASE ------------------------
!    
!      ! The outstanding feature with finite elements is: A basis
!      ! function for a DOF on one element has common support only
!      ! with the DOF's on the same element! E.g. for Q1:
!      !
!      !        #. . .#. . .#. . .#
!      !        .     .     .     .
!      !        .  *  .  *  .  *  .
!      !        #-----O-----O. . .#
!      !        |     |     |     .
!      !        |     | IEL |  *  .
!      !        #-----X-----O. . .#
!      !        |     |     |     .
!      !        |     |     |  *  .
!      !        #-----#-----#. . .#
!      !
!      ! --> On element IEL, the basis function at "X" only interacts
!      !     with the basis functions in "O". Elements in the 
!      !     neighbourhood ("*") have no support, therefore we only have
!      !     to collect all "O" DOF's.
!      !
!      ! Calculate the global DOF's into IdofsTrial / IdofsTest.
!      !
!      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
!      ! global DOF's of our BILF_NELEMSIM elements simultaneously.
!      CALL dof_locGlobMapping_mult(p_rdiscretisation, p_IelementList(IELset:IELmax), &
!                                  .TRUE.,IdofsTest)
!                                   
!      ! If the DOF's for the test functions are different, calculate them, too.
!      IF (.NOT.bIdenticalTrialAndTest) THEN
!        CALL dof_locGlobMapping_mult(p_rdiscretisation, p_IelementList(IELset:IELmax), &
!                                    .FALSE.,IdofsTrial)
!      END IF
!      
!      ! ------------------- LOCAL MATRIX SETUP PHASE -----------------------
!      
!      ! For the assembly of the global matrix, we use a "local"
!      ! approach. At first we build a "local" system matrix according
!      ! to the current element. This contains all additive
!      ! contributions of element IEL, which are later added at the
!      ! right positions to the elements in the global system matrix.
!      !
!      ! We have indofTrial trial DOF's per element and
!      ! indofTest test DOF's per element. Therefore there are
!      ! indofTrial*indofTest tupel of basis-/testfunctions (phi_i,psi_j) 
!      ! "active" (i.e. have common support) on our current element, each 
!      ! giving an additive contribution to the system matrix.
!      !
!      ! We build a quadratic indofTrial*indofTest local matrix:
!      ! Kentry(1..indofTrial,1..indofTest) receives the position 
!      !   in the global system matrix, where the corresponding value 
!      !   has to be added to.
!      ! (The corresponding contrbutions can be saved separately, 
!      !  but we directly add them to the global matrix in this 
!      !  approach.)
!      !
!      ! We build local matrices for all our elements 
!      ! in the set simultaneously.
!      ! Loop through elements in the set and for each element,
!      ! loop through the local matrices to initialise them:
!      DO IEL=1,IELmax-IELset+1
!      
!        ! For building the local matrices, we have first to
!        ! loop through the test functions (the "O"'s), as these
!        ! define the rows in the matrix.
!        DO IDOFE=1,indofTest
!        
!          ! Row IDOFE of the local matrix corresponds 
!          ! to row=global DOF KDFG(IDOFE) in the global matrix.
!          ! This is one of the the "O"'s in the above picture.
!          ! Get the starting position of the corresponding row
!          ! to JCOL0:
!
!          JCOL0=p_KLD(IdofsTest(IDOFE,IEL))
!          
!          ! Now we loop through the other DOF's on the current element
!          ! (the "O"'s).
!          ! All these have common support with our current basis function
!          ! and will therefore give an additive value to the global
!          ! matrix.
!          
!          DO JDOFE=1,indofTrial
!            
!            ! Get the global DOF of the "X" which interacts with 
!            ! our "O".
!            
!            JDFG=p_IdofsTrial(JDOFE,IEL)
!            
!            ! Starting in JCOL0 (which points to the beginning of
!            ! the line initially), loop through the elements in
!            ! the row to find the position of column IDFG.
!            ! Jump out of the DO loop if we find the column.
!            
!            DO JCOL=JCOL0,NA
!              IF (p_KCOL(JCOL) .EQ. JDFG) EXIT
!            END DO
!
!            ! Because columns in the global matrix are sorted 
!            ! ascendingly (except for the diagonal element),
!            ! the next search can start after the column we just found.
!            
!            ! JCOL0=JCOL+1
!            
!            ! Save the position of the matrix entry into the local
!            ! matrix.
!            ! Note that a column in Kentry corresponds to a row in
!            ! the real matrix. We aligned Kentry/DENTRY this way to get
!            ! higher speed of the assembly routine, since this leads
!            ! to better data locality.
!            
!            Kentry(JDOFE,IDOFE,IEL)=JCOL
!            
!          END DO ! IDOFE
!          
!        END DO ! JDOFE
!        
!      END DO ! IEL
!      
!      ! -------------------- ELEMENT EVALUATION PHASE ----------------------
!      
!      ! Ok, we found the positions of the local matrix entries
!      ! that we have to change.
!      ! To calculate the matrix contributions, we have to evaluate
!      ! the elements to give us the values of the basis functions
!      ! in all the DOF's in all the elements in our set.
!      !
!
!      ! We have the coordinates of the cubature points saved in the
!      ! coordinate array from above. Unfortunately for nonparametric
!      ! elements, we need the real coordinate.
!      ! Furthermore, we anyway need the coordinates of the element
!      ! corners and the Jacobian determinants corresponding to
!      ! all the points.
!      !
!      ! At first, get the coordinates of the corners of all the
!      ! elements in the current set. 
!      
!!      DO IEL=1,IELmax-IELset+1
!!        p_Dcoords(:,:,IEL) = p_DvertexCoords(:, &
!!                            p_IverticesAtElement(:,p_IelementList(IELset+IEL-1)))
!!      END DO
!!      DO IEL=1,IELmax-IELset+1
!!        DO J = 1,NVE
!!          DO I = 1,p_rtriangulation%ndim
!!            p_Dcoords(I,J,IEL) = p_DvertexCoords(I, &
!!                               p_IverticesAtElement(J,p_IelementList(IELset+IEL-1)))
!!          END DO
!!        END DO
!!      END DO
!
!      CALL trafo_getCoords_sim (elem_igetTrafoType(p_relementDistribution%itrialElement),&
!          p_rtriangulation,p_IelementList(IELset:IELmax),p_Dcoords)
!
!      ! Depending on the type of transformation, we must now choose
!      ! the mapping between the reference and the real element.
!      ! In case we use a nonparametric element or a nonconstant coefficient function,
!      ! we need the coordinates of the points on the real element, too.
!      IF (bnonparTrial .OR. bnonparTest .OR. (.NOT. rform%ballCoeffConstant)) THEN
!      
!        CALL trafo_calctrafo_sim (&
!             p_rdiscretisation%RelementDistr(icurrentElementDistr)%ctrafoType,&
!             IELmax-IELset+1,ncubp,p_Dcoords,&
!             p_DcubPtsRef,p_Djac(:,:,1:IELmax-IELset+1),p_Ddetj(:,1:IELmax-IELset+1),&
!             p_DcubPtsReal)
!      
!      ELSE
!      
!        CALL trafo_calctrafo_sim (p_relementDistribution%ctrafoType,&
!             IELmax-IELset+1,ncubp,p_Dcoords,&
!             p_DcubPtsRef,p_Djac(:,:,1:IELmax-IELset+1),p_Ddetj(:,1:IELmax-IELset+1))
!             
!      END IF
!      
!      ! If the matrix has nonconstant coefficients, calculate the coefficients now.
!      IF (.NOT. rform%ballCoeffConstant) THEN
!        rintSubset%ielementDistribution = icurrentElementDistr
!        rintSubset%ielementStartIdx = IELset
!        rintSubset%p_Ielements => p_IelementList(IELset:IELmax)
!        CALL fcoeff_buildMatrixSc_sim (p_rdiscretisation,rform, &
!                  IELmax-IELset+1_I32,ncubp,&
!                  p_DcubPtsReal,p_IdofsTrial,IdofsTest,rintSubset, &
!                  Dcoefficients,rcollection)
!      END IF
!      
!      ! If the element needs it, calculate the twist index array.
!      IF (ASSOCIATED(rintSubset%p_ItwistIndex)) THEN
!        CALL trafo_calcTwistIndices(p_rtriangulation,&
!            p_IelementList(IELset:IELmax),rintSubset%p_ItwistIndex)
!      END IF
!      
!      ! Calculate the values of the basis functions.
!      ! Pass p_DcubPts as point coordinates, which point either to the
!      ! coordinates on the reference element (the same for all elements)
!      ! or on the real element - depending on whether this is a 
!      ! parametric or nonparametric element.
!      CALL elem_generic_sim (p_relementDistribution%itestElement, p_Dcoords, &
!            p_Djac(:,:,1:IELmax-IELset+1), p_Ddetj(:,1:IELmax-IELset+1), &
!            BderTest, DbasTest, ncubp, IELmax-IELset+1, p_DcubPtsTest,&
!            rintSubset%p_ItwistIndex)
!            
!      ! Omit the calculation of the trial function values if they
!      ! are identical to the test function values.
!      IF (.NOT. bidenticalTrialAndTest) THEN
!        CALL elem_generic_sim (p_relementDistribution%itrialElement, p_Dcoords, &
!            p_Djac(:,:,1:IELmax-IELset+1), p_Ddetj(:,1:IELmax-IELset+1), &
!            BderTrial, DbasTrial, ncubp, IELmax-IELset+1, p_DcubPtsTrial,&
!            rintSubset%p_ItwistIndex)
!      END IF
!      
!      ! --------------------- DOF COMBINATION PHASE ------------------------
!      
!      ! Values of all basis functions calculated. Now we can start 
!      ! to integrate!
!      !
!      ! We have two different versions for the integration - one
!      ! with constant coefficients and one with nonconstant coefficients.
!      !
!      ! Check the bilinear form which one to use:
!      
!      IF (rform%ballCoeffConstant) THEN
!      
!        ! Constant coefficients. The coefficients are to be found in
!        ! the Dcoefficients variable of the form.
!        !
!        ! Loop over the elements in the current set.
!
!        DO IEL=1,IELmax-IELset+1
!          
!          ! Clear the local matrix
!          Dentry = 0.0_DP
!          
!          ! Loop over all cubature points on the current element
!          DO ICUBP = 1, ncubp
!
!            ! calculate the current weighting factor in the cubature formula
!            ! in that cubature point.
!            !
!            ! Take the absolut value of the determinant of the mapping.
!            ! In 2D, the determinant is always positive, whereas in 3D,
!            ! the determinant might be negative -- that's normal!
!
!            OM = Domega(ICUBP)*ABS(p_Ddetj(ICUBP,IEL))
!
!            ! Loop over the additive factors in the bilinear form.
!            DO IALBET = 1,rform%itermcount
!            
!              ! Get from Idescriptors the type of the derivatives for the 
!              ! test and trial functions. The summand we calculate
!              ! here will be added to the matrix entry:
!              !
!              ! a_ij  =  int_... ( psi_j )_IB  *  ( phi_i )_IA
!              !
!              ! -> Ix=0: function value, 
!              !      =1: first derivative, ...
!              !    as defined in the module 'derivative'.
!              
!              IA = rform%Idescriptors(1,IALBET)
!              IB = rform%Idescriptors(2,IALBET)
!              
!              ! Multiply OM with the coefficient of the form.
!              ! This gives the actual value to multiply the
!              ! function value with before summing up to the integral.
!              AUX = OM * rform%Dcoefficients(IALBET)
!            
!              ! Now loop through all possible combinations of DOF's
!              ! in the current cubature point. The outer loop
!              ! loops through the "O"'s in the above picture,
!              ! the test functions:
!
!              DO IDOFE=1,indofTest
!              
!                ! Get the value of the (test) basis function 
!                ! phi_i (our "O") in the cubature point:
!                DB = DbasTest(IDOFE,IB,ICUBP,IEL)
!                
!                ! Perform an inner loop through the other DOF's
!                ! (the "X"). 
!
!                DO JDOFE=1,indofTrial
!                
!                  ! Get the value of the basis function 
!                  ! psi_j (our "X") in the cubature point. 
!                  ! Them multiply:
!                  !    DB * DBAS(..) * AUX
!                  ! ~= phi_i * psi_j * coefficient * cub.weight
!                  ! Summing this up gives the integral, so the contribution
!                  ! to the global matrix. 
!                  !
!                  ! Simply summing up DB * DBAS(..) * AUX would give
!                  ! the coefficient of the local matrix. We save this
!                  ! contribution in the local matrix.
!
!                  !JCOLB = Kentry(JDOFE,IDOFE,IEL)
!                  !p_DA(JCOLB) = p_DA(JCOLB) + DB*p_DbasTrial(JDOFE,IA,ICUBP,IEL)*AUX
!                  Dentry(JDOFE,IDOFE) = Dentry(JDOFE,IDOFE) + &
!                                        DB*p_DbasTrial(JDOFE,IA,ICUBP,IEL)*AUX
!                
!                END DO ! JDOFE
!              
!              END DO ! IDOFE
!              
!            END DO ! IALBET
!
!          END DO ! ICUBP 
!          
!          ! Incorporate the local matrices into the global one.
!          ! Kentry gives the position of the additive contributions in Dentry.
!          !
!          ! OpenMP-Extension: This is a critical section. Only one thread is
!          ! allowed to write to the matrix, otherwise the matrix may get
!          ! messed up.
!          ! The critical section is put around both loops as indofTest/indofTrial
!          ! are usually small and quickly to handle.
!          !
!          !%OMP CRITICAL
!          DO IDOFE=1,indofTest
!            DO JDOFE=1,indofTrial
!              p_DA(Kentry(JDOFE,IDOFE,IEL)) = p_DA(Kentry(JDOFE,IDOFE,IEL)) + &
!                                              Dentry(JDOFE,IDOFE)
!            END DO
!          END DO
!          !%OMP END CRITICAL
!
!        END DO ! IEL
!        
!      ELSE
!      
!        ! Nonconstant coefficients. The coefficients are to be found in
!        ! the Dcoefficients variable as computed above.
!        !
!        ! Loop over the elements in the current set.
!
!        DO IEL=1,IELmax-IELset+1
!          
!          ! Clear the local matrix
!          Dentry = 0.0_DP
!          
!          ! Loop over all cubature points on the current element
!          DO ICUBP = 1, ncubp
!
!            ! calculate the current weighting factor in the cubature formula
!            ! in that cubature point.
!            !
!            ! Take the absolut value of the determinant of the mapping.
!            ! In 2D, the determinant is always positive, whereas in 3D,
!            ! the determinant might be negative -- that's normal!
!
!            OM = Domega(ICUBP)*ABS(p_Ddetj(ICUBP,IEL))
!
!            ! Loop over the additive factors in the bilinear form.
!            DO IALBET = 1,rform%itermcount
!            
!              ! Get from Idescriptors the type of the derivatives for the 
!              ! test and trial functions. The summand we calculate
!              ! here will be added to the matrix entry:
!              !
!              ! a_ij  =  int_... ( psi_j )_IA  *  ( phi_i )_IB
!              !
!              ! -> Ix=0: function value, 
!              !      =1: first derivative, ...
!              !    as defined in the module 'derivative'.
!              
!              IA = rform%Idescriptors(1,IALBET)
!              IB = rform%Idescriptors(2,IALBET)
!              
!              ! Multiply OM with the coefficient of the form.
!              ! This gives the actual value to multiply the
!              ! function value with before summing up to the integral.
!              ! Get the precalculated coefficient from the coefficient array.
!              AUX = OM * Dcoefficients(IALBET,ICUBP,IEL)
!            
!              ! Now loop through all possible combinations of DOF's
!              ! in the current cubature point. The outer loop
!              ! loops through the "O" in the above picture,
!              ! the test functions:
!
!              DO IDOFE=1,indofTest
!                
!                ! Get the value of the (test) basis function 
!                ! phi_i (our "O") in the cubature point:
!                DB = DbasTest(IDOFE,IB,ICUBP,IEL)
!                
!                ! Perform an inner loop through the other DOF's
!                ! (the "X"). 
!
!                DO JDOFE=1,indofTrial
!              
!                  ! Get the value of the basis function 
!                  ! psi_j (our "X") in the cubature point. 
!                  ! Them multiply:
!                  !    DB * DBAS(..) * AUX
!                  ! ~= phi_i * psi_j * coefficient * cub.weight
!                  ! Summing this up gives the integral, so the contribution
!                  ! to the global matrix. 
!                  !
!                  ! Simply summing up DB * DBAS(..) * AUX would give
!                  ! the coefficient of the local matrix. We save this
!                  ! contribution in the local matrix of element IEL.
!
!                  !JCOLB = Kentry(JDOFE,IDOFE,IEL)
!                  !p_DA(JCOLB) = p_DA(JCOLB) + DB*p_DbasTrial(JDOFE,IA,ICUBP,IEL)*AUX
!                  Dentry(JDOFE,IDOFE) = &
!                      Dentry(JDOFE,IDOFE)+DB*p_DbasTrial(JDOFE,IA,ICUBP,IEL)*AUX
!                
!                END DO
!              
!              END DO ! JDOFE
!              
!            END DO ! IALBET
!
!          END DO ! ICUBP 
!          
!          ! Incorporate the local matrices into the global one.
!          ! Kentry gives the position of the additive contributions in Dentry.
!          !
!          ! OpenMP-Extension: This is a critical section. Only one thread is
!          ! allowed to write to the matrix, otherwise the matrix may get
!          ! messed up.
!          ! The critical section is put around both loops as indofTest/indofTrial
!          ! are usually small and quickly to handle.
!          !
!          !%OMP CRITICAL
!          DO IDOFE=1,indofTest
!            DO JDOFE=1,indofTrial
!              p_DA(Kentry(JDOFE,IDOFE,IEL)) = &
!                  p_DA(Kentry(JDOFE,IDOFE,IEL)) + Dentry(JDOFE,IDOFE)
!            END DO
!          END DO
!          !%OMP END CRITICAL
!
!        END DO ! IEL
!
!      END IF ! rform%ballCoeffConstant
!
!    END DO ! IELset
!    !%OMP END DO
!    
!    IF (.NOT. rform%ballCoeffConstant) THEN
!      DEALLOCATE(Dcoefficients)
!    END IF
!    DEALLOCATE(IdofsTrial)
!    DEALLOCATE(IdofsTest)
!    DEALLOCATE(DbasTrial)
!    DEALLOCATE(DbasTest)
!    DEALLOCATE(Kentry)
!    DEALLOCATE(Dentry)
!
!    ! Release memory
!    CALL domint_doneIntegration(rintSubset)
!
!    !%OMP END PARALLEL
!
!  END DO ! icurrentElementDistr
!
!  ! Finish
!  
!  END SUBROUTINE
!
  !****************************************************************************

!<subroutine>

  subroutine bilf_buildMatrix9d_conf3 (rform,bclear,rmatrixScalar,&
                                       fcoeff_buildMatrixSc_sim,rcollection)
  
!<description>
  ! This routine calculates the entries of a finite element matrix.
  ! The matrix structure must be prepared with bilf_createMatrixStructure
  ! in advance. The discretisation is assumed to be conformal, i.e. the DOF's
  ! of all finite elements must 'match'. Trial and test functions may be
  ! different.
  ! In case the array for the matrix entries does not exist, the routine
  ! allocates memory in size of the matrix of the heap for the matrix entries.
  !
  ! For setting up the entries, the discretisation structure attached to
  ! the matrix is used (rmatrixScalar%p_rdiscretisation). This is
  ! normally attached to the matrix by bilf_createMatrixStructure.
  !
  ! Double-precision version.
!</description>

!<input>
  ! The bilinear form specifying the underlying PDE of the discretisation.
  type(t_bilinearForm), intent(IN) :: rform
  
  ! Whether to clear the matrix before calculating the entries.
  ! If .FALSE., the new matrix entries are added to the existing entries.
  logical, intent(IN) :: bclear
  
  ! OPTIONAL: A pointer to a collection structure. This structure is given to the
  ! callback function for nonconstant coefficients to provide additional
  ! information. 
  type(t_collection), intent(INOUT), target, optional :: rcollection
  
  ! OPTIONAL: A callback routine for nonconstant coefficient matrices.
  ! Must be present if the matrix has nonconstant coefficients!
  include 'intf_coefficientMatrixSc.inc'
  optional :: fcoeff_buildMatrixSc_sim
!</input>

!<inputoutput>
  ! The FE matrix. Calculated matrix entries are imposed to this matrix.
  type(t_matrixScalar), intent(INOUT) :: rmatrixScalar
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i,i1,k,icurrentElementDistr,JDFG, ICUBP, IALBET, IA, IB
  logical :: bIdenticalTrialAndTest
  integer :: IEL, IELmax, IELset, IDOFE, JDOFE
  integer :: JCOL0,JCOL
  real(DP) :: OM,AUX, DB
  
  ! Array to tell the element which derivatives to calculate
  logical, dimension(EL_MAXNDER) :: BderTrialTempl, BderTestTempl, BderTrial, BderTest

  ! For every cubature point on the reference element,
  ! the corresponding cubature weight
  real(DP), dimension(:), allocatable :: Domega
  
  ! number of cubature points on the reference element
  integer :: ncubp
  
  ! Pointer to KLD, KCOL, DA
  integer, dimension(:), pointer :: p_KLD, p_KCOL
  real(DP), dimension(:), pointer :: p_DA
  
  ! An allocateable array accepting the DOF's of a set of elements.
  integer, dimension(:,:), allocatable, target :: IdofsTest, IdofsTrial
  integer, dimension(:,:), pointer :: p_IdofsTrial
  
  ! Allocateable arrays for the values of the basis functions - 
  ! for test and trial spaces.
  real(DP), dimension(:,:,:,:), allocatable, target :: DbasTest,DbasTrial
  real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial
  
  ! Number of entries in the matrix - for quicker access
  integer :: NA,NVE
  integer :: NEQ
  
  ! Type of transformation from the reference to the real element 
  integer(I32) :: ctrafoType
  
  ! Element evaluation tag; collects some information necessary for evaluating
  ! the elements.
  integer(I32) :: cevaluationTag
  
  ! Number of local degees of freedom for trial and test functions
  integer :: indofTrial, indofTest
  
  ! The triangulation structure - to shorten some things...
  type(t_triangulation), pointer :: p_rtriangulation
  
  ! A pointer to an element-number list
  integer, dimension(:), pointer :: p_IelementList
  
  ! Local matrices, used during the assembly.
  ! Values and positions of values in the global matrix.
  integer, dimension(:,:,:), allocatable :: Kentry
  real(DP), dimension(:,:), allocatable :: Dentry

  ! An array that takes coordinates of the cubature formula on the reference element
  real(DP), dimension(:,:), allocatable :: p_DcubPtsRef
  
  ! Pointer to the jacobian determinants
  real(DP), dimension(:,:), pointer :: p_Ddetj

  ! Current element distribution
  type(t_elementDistribution), pointer :: p_relementDistrTest
  type(t_elementDistribution), pointer :: p_relementDistrTrial
  
  ! Number of elements in a block. Normally =BILF_NELEMSIM,
  ! except if there are less elements in the discretisation.
  integer :: nelementsPerBlock
  
  ! Some variables to support nonconstant coefficients in the matrix.
  
  ! Pointer to the coefficients that are computed by the callback routine.
  real(DP), dimension(:,:,:), allocatable :: Dcoefficients
  
  ! A t_domainIntSubset structure that is used for storing information
  ! and passing it to callback routines as well as element evaluation routines.
  type(t_domainIntSubset) :: rintSubset
  type(t_evalElementSet) :: revalElementSet
  logical :: bcubPtsInitialised
  
  ! The discretisation - for easier access
  type(t_spatialDiscretisation), pointer :: p_rdiscrTest
  type(t_spatialDiscretisation), pointer :: p_rdiscrTrial
  
  if ((.not. associated(rmatrixScalar%p_rspatialDiscrTest)) .or. &
      (.not. associated(rmatrixScalar%p_rspatialDiscrTrial))) then
    call output_line ('No discretisation associated!',&
        OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrix9d_conf2')
    call sys_halt()
  end if

  ! Which derivatives of basis functions are needed?
  ! Check the descriptors of the bilinear form and set BDERxxxx
  ! according to these.

  BderTrialTempl = .false.
  BderTestTempl = .false.
  
  ! Loop through the additive terms
  do i=1,rform%itermCount
    ! The desriptor Idescriptors gives directly the derivative
    ! which is to be computed! Build template's for BDER.
    ! We don't compute the actual BDER here, as there might be some special
    ! processing if trial/test functions are identical!
    !
    ! At first build the descriptors for the trial functions
    I1=rform%Idescriptors(1,I)
    
    if ((I1 .le.0) .or. (I1 .gt. DER_MAXNDER)) then
      call output_line ('Invalid descriptor!',&
          OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrix9d_conf2')
      call sys_halt()
    endif
    
    BderTrialTempl(I1)=.true.

    ! Then those of the test functions
    I1=rform%Idescriptors(2,I)
    
    if ((I1 .le.0) .or. (I1 .gt. DER_MAXNDER)) then
      call output_line ('Invalid descriptor!',&
          OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrix9d_conf2')
      call sys_halt()
    endif
    
    BderTestTempl(I1)=.true.
  end do
  
  ! Get information about the matrix:
  NA = rmatrixScalar%NA
  NEQ = rmatrixScalar%NEQ
  
  ! We need KCOL/KLD of our matrix
  if ((rmatrixScalar%h_KCOL .eq. ST_NOHANDLE) .or. &
      (rmatrixScalar%h_KLD .eq. ST_NOHANDLE)) then
    call output_line ('No discretisation structure! Cannot assemble matrix!', &
                      OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrix9d_conf2')
    call sys_halt()
  end if
  
  call lsyssc_getbase_Kcol (rmatrixScalar,p_KCOL)
  call lsyssc_getbase_Kld (rmatrixScalar,p_KLD)
  
  ! Check if the matrix entries exist. If not, allocate the matrix.
  if (rmatrixScalar%h_DA .eq. ST_NOHANDLE) then

    ! Clear the entries in the matrix - we need to start with zero
    ! when assembling a new matrix!
    call storage_new1D ('bilf_buildMatrix9d_conf', 'DA', &
                        NA, ST_DOUBLE, rmatrixScalar%h_DA, &
                        ST_NEWBLOCK_ZERO)
    call lsyssc_getbase_double (rmatrixScalar,p_DA)

  else
  
    call lsyssc_getbase_double (rmatrixScalar,p_DA)

    ! If desired, clear the matrix before assembling.
    if (bclear) then
      call lalg_clearVectorDble (p_DA)
    end if
    
  end if
  
  ! Get the discretisation
  p_rdiscrTest => rmatrixScalar%p_rspatialDiscrTest
  p_rdiscrTrial => rmatrixScalar%p_rspatialDiscrTrial
  
  if ((.not. associated(p_rdiscrTest)) .or. &
      (.not. associated(p_rdiscrTrial))) then
    call output_line ('No discretisation attached to the matrix!',&
        OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrix9d_conf2')
    call sys_halt()
  end if
  
  ! Get a pointer to the triangulation - for easier access.
  p_rtriangulation => p_rdiscrTest%p_rtriangulation
  
  ! For saving some memory in smaller discretisations, we calculate
  ! the number of elements per block. For smaller triangulations,
  ! this is NEL. If there are too many elements, it's at most
  ! BILF_NELEMSIM. This is only used for allocating some arrays.
  nelementsPerBlock = min(BILF_NELEMSIM,p_rtriangulation%NEL)
  
  ! Now loop over the different element distributions (=combinations
  ! of trial and test functions) in the discretisation.

  do icurrentElementDistr = 1,p_rdiscrTest%inumFESpaces
  
    ! Activate the current element distribution
    p_relementDistrTest => p_rdiscrTest%RelementDistr(icurrentElementDistr)
    p_relementDistrTrial => p_rdiscrTrial%RelementDistr(icurrentElementDistr)
  
    ! Cancel if this element distribution is empty.
    if (p_relementDistrTest%NEL .eq. 0) cycle
    
    ! Get the number of local DOF's for trial and test functions
    indofTrial = elem_igetNDofLoc(p_relementDistrTrial%celement)
    indofTest = elem_igetNDofLoc(p_relementDistrTest%celement)
    
    ! Get the number of vertices of the element, specifying the transformation
    ! form the reference to the real element.
    NVE = elem_igetNVE(p_relementDistrTest%celement)
    if (NVE .ne. elem_igetNVE(p_relementDistrTrial%celement)) then
      call output_line ('Element spaces incompatible!',&
          OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrix9d_conf2')
      call sys_halt()
    end if
    
    ! Get from the element space the type of coordinate system
    ! that is used there:
    ctrafoType = elem_igetTrafoType(p_relementDistrTest%celement)
    
    ! Get the number of cubature points for the cubature formula
    ncubp = cub_igetNumPts(p_relementDistrTest%ccubTypeBilForm)
    
    ! Allocate two arrays for the points and the weights
    allocate(Domega(ncubp))
    allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),ncubp))
    
    ! Get the cubature formula
    call cub_getCubature(p_relementDistrTest%ccubTypeBilForm,p_DcubPtsRef, Domega)
    
    ! Open-MP-Extension: Open threads here.
    ! Each thread will allocate its own local memory...
    !
    !%OMP PARALLEL PRIVATE(rintSubset, revalElementSet,&
    !%OMP   cevaluationTag, kentry, dentry,DbasTest,DbasTrial, &
    !%OMP   IdofsTest,IdofsTrial,p_Ddetj,&
    !%OMP   Dcoefficients, bIdenticalTrialandTest, p_IdofsTrial, &
    !%OMP   p_DbasTrial, BderTrial, BderTest, ielmax,IEL, idofe,jdofe,jcol0, &
    !%OMP   jcol,JDFG,ICUBP, IALBET,OM,ia,ib,aux, db)
    
    ! Quickly check if one of the specified derivatives is out of the allowed range:
    do IALBET = 1,rform%itermcount
      IA = rform%Idescriptors(1,IALBET)
      IB = rform%Idescriptors(2,IALBET)      
      if ((IA.lt.0) .or. &
          (IA .gt. elem_getMaxDerivative(p_relementDistrTrial%celement))) then
        call output_line ('Specified trial-derivative '//trim(sys_siL(IA,10))//&
            ' not available!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrix9d_conf2')
        call sys_halt()
      end if

      if ((IB.lt.0) .or. &
          (IB .gt. elem_getMaxDerivative(p_relementDistrTest%celement))) then
        call output_line ('Specified test-derivative '//trim(sys_siL(IA,10))//&
            ' not available!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrix9d_conf2')
        call sys_halt()
      end if
    end do
    
    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  ALLOCATE(DbasTest(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    !  ALLOCATE(DbasTrial(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    ! would lead to nonused memory blocks in these arrays during the assembly, 
    ! which reduces the speed by 50%!
    
    allocate(DbasTest(indofTest,&
             elem_getMaxDerivative(p_relementDistrTest%celement),&
             ncubp,nelementsPerBlock))
    allocate(DbasTrial(indofTrial,&
             elem_getMaxDerivative(p_relementDistrTrial%celement), &
             ncubp,nelementsPerBlock))

    ! Allocate memory for the DOF's of all the elements.
    allocate(IdofsTest(indofTest,nelementsPerBlock))
    allocate(IdofsTrial(indofTrial,nelementsPerBlock))

    ! Allocate an array saving the local matrices for all elements
    ! in an element set.
    ! We could also allocate EL_MAXNBAS*EL_MAXNBAS*BILF_NELEMSIM integers
    ! for this local matrix, but this would normally not fit to the cache
    ! anymore! indofTrial*indofTest*BILF_NELEMSIM is normally much smaller!
    allocate(Kentry(indofTrial,indofTest,nelementsPerBlock))
    allocate(Dentry(indofTrial,indofTest))
    
    ! Initialisation of the element set.
    call elprep_init(revalElementSet)
    
    ! Indicate that cubature points must still be initialised in the element set.
    bcubPtsInitialised = .false.
    
    ! In case of nonconstant coefficients in that part of the matrix, we
    ! need an additional array to save all the coefficients:
    if (.not. rform%BconstantCoeff(icurrentElementDistr)) then
      if (rform%ballCoeffConstant) then
        call output_line ('Some oefficients are not constant &
                &although they should be!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrix9d_conf2')
        call sys_halt()
      end if
      if (.not. present(fcoeff_buildMatrixSc_sim)) then
        call output_line ('Coefficient function not given!',&
            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrix9d_conf2')
        call sys_halt()
      end if
      allocate(Dcoefficients(rform%itermCount,ncubp,nelementsPerBlock))
    end if
                    
    ! p_IdofsTest points either to the just created array or to the
    ! array with the DOF's of the trial functions - when trial and
    ! test functions are identical.
    ! We don't rely on bidenticalTrialAndTest purely, as this does not
    ! indicate whether there are identical trial and test functions
    ! in one block!
    bIdenticalTrialAndTest = p_relementDistrTrial%celement .eq. &
                             p_relementDistrTest%celement

    ! Let p_IdofsTrial point either to IdofsTrial or to the DOF's of the test
    ! space IdofTest (if both spaces are identical). 
    ! We create a pointer for the trial space and not for the test space to
    ! prevent pointer-arithmetic in the innerst loop below!
    if (bIdenticalTrialAndTest) then
      p_IdofsTrial => IdofsTest
      p_DbasTrial  => DbasTest
      ! Build the actual combination of what the element should calculate.
      ! As we evaluate only once, what the element must calculate is an
      ! OR combination of the BDER from trial and test functions.
      BderTrial = BderTrialTempl .or. BderTestTempl
      BderTest = BderTrial
    else
      p_IdofsTrial => IdofsTrial
      p_DbasTrial  => DbasTrial
      
      ! Build the actual combination of what the element should calculate.
      ! Copy BDERxxxx to BDERxxxxAct
      BderTrial = BderTrialTempl
      BderTest = BderTestTempl
    end if
    
    ! p_IelementList must point to our set of elements in the discretisation
    ! with that combination of trial/test functions
    call storage_getbase_int (p_relementDistrTest%h_IelementList, &
                              p_IelementList)
                              
    ! Loop over the elements - blockwise.
    !
    ! Open-MP-Extension: Each loop cycle is executed in a different thread,
    ! so BILF_NELEMSIM local matrices are simultaneously calculated in the
    ! inner loop(s).
    ! The blocks have all the same size, so we can use static scheduling.
    !
    !%OMP do schedule(static,1)
    do IELset = 1, size(p_IelementList), BILF_NELEMSIM
    
      ! We always handle BILF_NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most BILF_NELEMSIM
      ! elements simultaneously.
      
      IELmax = min(size(p_IelementList),IELset-1+BILF_NELEMSIM)
    
      ! --------------------- DOF SEARCH PHASE ------------------------
    
      ! The outstanding feature with finite elements is: A basis
      ! function for a DOF on one element has common support only
      ! with the DOF's on the same element! E.g. for Q1:
      !
      !        #. . .#. . .#. . .#
      !        .     .     .     .
      !        .  *  .  *  .  *  .
      !        #-----O-----O. . .#
      !        |     |     |     .
      !        |     | IEL |  *  .
      !        #-----X-----O. . .#
      !        |     |     |     .
      !        |     |     |  *  .
      !        #-----#-----#. . .#
      !
      ! --> On element IEL, the basis function at "X" only interacts
      !     with the basis functions in "O". Elements in the 
      !     neighbourhood ("*") have no support, therefore we only have
      !     to collect all "O" DOF's.
      !
      ! Calculate the global DOF's into IdofsTrial / IdofsTest.
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF's of our BILF_NELEMSIM elements simultaneously.
      call dof_locGlobMapping_mult(p_rdiscrTest, p_IelementList(IELset:IELmax), &
                                   IdofsTest)
                                   
      ! If the DOF's for the test functions are different, calculate them, too.
      if (.not.bIdenticalTrialAndTest) then
        call dof_locGlobMapping_mult(p_rdiscrTrial, p_IelementList(IELset:IELmax), &
                                     IdofsTrial)
      end if
      
      ! ------------------- LOCAL MATRIX SETUP PHASE -----------------------
      
      ! For the assembly of the global matrix, we use a "local"
      ! approach. At first we build a "local" system matrix according
      ! to the current element. This contains all additive
      ! contributions of element IEL, which are later added at the
      ! right positions to the elements in the global system matrix.
      !
      ! We have indofTrial trial DOF's per element and
      ! indofTest test DOF's per element. Therefore there are
      ! indofTrial*indofTest tupel of basis-/testfunctions (phi_i,psi_j) 
      ! "active" (i.e. have common support) on our current element, each 
      ! giving an additive contribution to the system matrix.
      !
      ! We build a quadratic indofTrial*indofTest local matrix:
      ! Kentry(1..indofTrial,1..indofTest) receives the position 
      !   in the global system matrix, where the corresponding value 
      !   has to be added to.
      ! (The corresponding contrbutions can be saved separately, 
      !  but we directly add them to the global matrix in this 
      !  approach.)
      !
      ! We build local matrices for all our elements 
      ! in the set simultaneously.
      ! Loop through elements in the set and for each element,
      ! loop through the local matrices to initialise them:
      do IEL=1,IELmax-IELset+1
      
        ! For building the local matrices, we have first to
        ! loop through the test functions (the "O"'s), as these
        ! define the rows in the matrix.
        do IDOFE=1,indofTest
        
          ! Row IDOFE of the local matrix corresponds 
          ! to row=global DOF KDFG(IDOFE) in the global matrix.
          ! This is one of the the "O"'s in the above picture.
          ! Get the starting position of the corresponding row
          ! to JCOL0:

          JCOL0=p_KLD(IdofsTest(IDOFE,IEL))
          
          ! Now we loop through the other DOF's on the current element
          ! (the "O"'s).
          ! All these have common support with our current basis function
          ! and will therefore give an additive value to the global
          ! matrix.
          
          do JDOFE=1,indofTrial
            
            ! Get the global DOF of the "X" which interacts with 
            ! our "O".
            
            JDFG=p_IdofsTrial(JDOFE,IEL)
            
            ! Starting in JCOL0 (which points to the beginning of
            ! the line initially), loop through the elements in
            ! the row to find the position of column IDFG.
            ! Jump out of the DO loop if we find the column.
            
            do JCOL=JCOL0,NA
              if (p_KCOL(JCOL) .eq. JDFG) exit
            end do

            ! Because columns in the global matrix are sorted 
            ! ascendingly (except for the diagonal element),
            ! the next search can start after the column we just found.
            
            ! JCOL0=JCOL+1
            
            ! Save the position of the matrix entry into the local
            ! matrix.
            ! Note that a column in Kentry corresponds to a row in
            ! the real matrix. We aligned Kentry/DENTRY this way to get
            ! higher speed of the assembly routine, since this leads
            ! to better data locality.
            
            Kentry(JDOFE,IDOFE,IEL)=JCOL
            
          end do ! IDOFE
          
        end do ! JDOFE
        
      end do ! IEL
      
      ! -------------------- ELEMENT EVALUATION PHASE ----------------------
      
      ! Ok, we found the positions of the local matrix entries
      ! that we have to change.
      ! To calculate the matrix contributions, we have to evaluate
      ! the elements to give us the values of the basis functions
      ! in all the DOF's in all the elements in our set.

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag. 
      cevaluationTag = elem_getEvaluationTag(p_relementDistrTrial%celement)
      cevaluationTag = ior(cevaluationTag,&
                      elem_getEvaluationTag(p_relementDistrTest%celement))
                      
      if (.not. rform%ballCoeffConstant) then
        ! Evaluate real coordinates if not necessary.
        cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REALPOINTS)
      end if
      
      ! In the first loop, calculate the coordinates on the reference element.
      ! In all later loops, use the precalculated information.
      !
      ! Note: Why not using
      !   IF (IELset .EQ. 1) THEN
      ! here, but this strange concept with the boolean variable?
      ! Because the IF-command does not work with OpenMP! bcubPtsInitialised
      ! is a local variable and will therefore ensure that every thread
      ! is initialising its local set of cubature points!
      if (.not. bcubPtsInitialised) then
        bcubPtsInitialised = .true.
        cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REFPOINTS)
      else
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
      end if

      ! Calculate all information that is necessary to evaluate the finite element
      ! on all cells of our subset. This includes the coordinates of the points
      ! on the cells.
      call elprep_prepareSetForEvaluation (revalElementSet,&
          cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
          ctrafoType, p_DcubPtsRef(:,1:ncubp))
      p_Ddetj => revalElementSet%p_Ddetj
      
      ! If the matrix has nonconstant coefficients, calculate the coefficients now.
      if (.not. rform%ballCoeffConstant) then
        call domint_initIntegrationByEvalSet (revalElementSet,rintSubset)
        rintSubset%ielementDistribution = icurrentElementDistr
        rintSubset%ielementStartIdx = IELset
        rintSubset%p_Ielements => p_IelementList(IELset:IELmax)
        rintSubset%p_IdofsTrial => p_IdofsTrial
        rintSubset%celement = p_relementDistrTrial%celement
        call fcoeff_buildMatrixSc_sim (p_rdiscrTest,p_rdiscrTrial,rform, &
                  IELmax-IELset+1,ncubp,&
                  revalElementSet%p_DpointsReal,&
                  p_IdofsTrial,IdofsTest,rintSubset, Dcoefficients,rcollection)
        call domint_doneIntegration (rintSubset)
      end if
      
      ! Calculate the values of the basis functions.
      call elem_generic_sim2 (p_relementDistrTest%celement, &
          revalElementSet, BderTest, DbasTest)
      
      ! Omit the calculation of the trial function values if they
      ! are identical to the test function values.
      if (.not. bidenticalTrialAndTest) then
        call elem_generic_sim2 (p_relementDistrTrial%celement, &
            revalElementSet, BderTrial, DbasTrial)
      end if
      
      ! --------------------- DOF COMBINATION PHASE ------------------------
      
      ! Values of all basis functions calculated. Now we can start 
      ! to integrate!
      !
      ! We have two different versions for the integration - one
      ! with constant coefficients and one with nonconstant coefficients.
      !
      ! Check the bilinear form which one to use:
      
      if (rform%ballCoeffConstant) then
      
        ! Constant coefficients. The coefficients are to be found in
        ! the Dcoefficients variable of the form.
        !
        ! Loop over the elements in the current set.

        do IEL=1,IELmax-IELset+1
          
          ! Clear the local matrix
          Dentry = 0.0_DP
          
          ! Loop over all cubature points on the current element
          do ICUBP = 1, ncubp

            ! calculate the current weighting factor in the cubature formula
            ! in that cubature point.
            !
            ! Take the absolut value of the determinant of the mapping.
            ! In 2D, the determinant is always positive, whereas in 3D,
            ! the determinant might be negative -- that's normal!

            OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))

            ! Loop over the additive factors in the bilinear form.
            do IALBET = 1,rform%itermcount
            
              ! Get from Idescriptors the type of the derivatives for the 
              ! test and trial functions. The summand we calculate
              ! here will be added to the matrix entry:
              !
              ! a_ij  =  int_... ( psi_j )_IB  *  ( phi_i )_IA
              !
              ! -> Ix=0: function value, 
              !      =1: first derivative, ...
              !    as defined in the module 'derivative'.
              
              IA = rform%Idescriptors(1,IALBET)
              IB = rform%Idescriptors(2,IALBET)
              
              ! Multiply OM with the coefficient of the form.
              ! This gives the actual value to multiply the
              ! function value with before summing up to the integral.
              AUX = OM * rform%Dcoefficients(IALBET)
            
              ! Now loop through all possible combinations of DOF's
              ! in the current cubature point. The outer loop
              ! loops through the "O"'s in the above picture,
              ! the test functions:

              do IDOFE=1,indofTest
              
                ! Get the value of the (test) basis function 
                ! phi_i (our "O") in the cubature point:
                DB = DbasTest(IDOFE,IB,ICUBP,IEL)
                
                ! Perform an inner loop through the other DOF's
                ! (the "X"). 

                do JDOFE=1,indofTrial
                
                  ! Get the value of the basis function 
                  ! psi_j (our "X") in the cubature point. 
                  ! Them multiply:
                  !    DB * DBAS(..) * AUX
                  ! ~= phi_i * psi_j * coefficient * cub.weight
                  ! Summing this up gives the integral, so the contribution
                  ! to the global matrix. 
                  !
                  ! Simply summing up DB * DBAS(..) * AUX would give
                  ! the coefficient of the local matrix. We save this
                  ! contribution in the local matrix.

                  !JCOLB = Kentry(JDOFE,IDOFE,IEL)
                  !p_DA(JCOLB) = p_DA(JCOLB) + DB*p_DbasTrial(JDOFE,IA,ICUBP,IEL)*AUX
                  Dentry(JDOFE,IDOFE) = Dentry(JDOFE,IDOFE) + &
                                        DB*p_DbasTrial(JDOFE,IA,ICUBP,IEL)*AUX
                
                end do ! JDOFE
              
              end do ! IDOFE
              
            end do ! IALBET

          end do ! ICUBP 
          
          ! Incorporate the local matrices into the global one.
          ! Kentry gives the position of the additive contributions in Dentry.
          !
          ! OpenMP-Extension: This is a critical section. Only one thread is
          ! allowed to write to the matrix, otherwise the matrix may get
          ! messed up.
          ! The critical section is put around both loops as indofTest/indofTrial
          ! are usually small and quickly to handle.
          !
          !%OMP CRITICAL
          do IDOFE=1,indofTest
            do JDOFE=1,indofTrial
              p_DA(Kentry(JDOFE,IDOFE,IEL)) = p_DA(Kentry(JDOFE,IDOFE,IEL)) + &
                                              Dentry(JDOFE,IDOFE)
            end do
          end do
          !%OMP END CRITICAL

        end do ! IEL
        
      else
      
        ! Nonconstant coefficients. The coefficients are to be found in
        ! the Dcoefficients variable as computed above.
        !
        ! Loop over the elements in the current set.

        do IEL=1,IELmax-IELset+1
          
          ! Clear the local matrix
          Dentry = 0.0_DP
          
          ! Loop over all cubature points on the current element
          do ICUBP = 1, ncubp

            ! calculate the current weighting factor in the cubature formula
            ! in that cubature point.
            !
            ! Take the absolut value of the determinant of the mapping.
            ! In 2D, the determinant is always positive, whereas in 3D,
            ! the determinant might be negative -- that's normal!

            OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,IEL))

            ! Loop over the additive factors in the bilinear form.
            do IALBET = 1,rform%itermcount
            
              ! Get from Idescriptors the type of the derivatives for the 
              ! test and trial functions. The summand we calculate
              ! here will be added to the matrix entry:
              !
              ! a_ij  =  int_... ( psi_j )_IA  *  ( phi_i )_IB
              !
              ! -> Ix=0: function value, 
              !      =1: first derivative, ...
              !    as defined in the module 'derivative'.
              
              IA = rform%Idescriptors(1,IALBET)
              IB = rform%Idescriptors(2,IALBET)
              
              ! Multiply OM with the coefficient of the form.
              ! This gives the actual value to multiply the
              ! function value with before summing up to the integral.
              ! Get the precalculated coefficient from the coefficient array.
              AUX = OM * Dcoefficients(IALBET,ICUBP,IEL)
            
              ! Now loop through all possible combinations of DOF's
              ! in the current cubature point. The outer loop
              ! loops through the "O" in the above picture,
              ! the test functions:

              do IDOFE=1,indofTest
                
                ! Get the value of the (test) basis function 
                ! phi_i (our "O") in the cubature point:
                DB = DbasTest(IDOFE,IB,ICUBP,IEL)
                
                ! Perform an inner loop through the other DOF's
                ! (the "X"). 

                do JDOFE=1,indofTrial
              
                  ! Get the value of the basis function 
                  ! psi_j (our "X") in the cubature point. 
                  ! Them multiply:
                  !    DB * DBAS(..) * AUX
                  ! ~= phi_i * psi_j * coefficient * cub.weight
                  ! Summing this up gives the integral, so the contribution
                  ! to the global matrix. 
                  !
                  ! Simply summing up DB * DBAS(..) * AUX would give
                  ! the coefficient of the local matrix. We save this
                  ! contribution in the local matrix of element IEL.

                  !JCOLB = Kentry(JDOFE,IDOFE,IEL)
                  !p_DA(JCOLB) = p_DA(JCOLB) + DB*p_DbasTrial(JDOFE,IA,ICUBP,IEL)*AUX
                  Dentry(JDOFE,IDOFE) = &
                      Dentry(JDOFE,IDOFE)+DB*p_DbasTrial(JDOFE,IA,ICUBP,IEL)*AUX
                
                end do
              
              end do ! JDOFE
              
            end do ! IALBET

          end do ! ICUBP 
          
          ! Incorporate the local matrices into the global one.
          ! Kentry gives the position of the additive contributions in Dentry.
          !
          ! OpenMP-Extension: This is a critical section. Only one thread is
          ! allowed to write to the matrix, otherwise the matrix may get
          ! messed up.
          ! The critical section is put around both loops as indofTest/indofTrial
          ! are usually small and quickly to handle.
          !
          !%OMP CRITICAL
          do IDOFE=1,indofTest
            do JDOFE=1,indofTrial
              p_DA(Kentry(JDOFE,IDOFE,IEL)) = &
                  p_DA(Kentry(JDOFE,IDOFE,IEL)) + Dentry(JDOFE,IDOFE)
            end do
          end do
          !%OMP END CRITICAL

        end do ! IEL

      end if ! rform%ballCoeffConstant

    end do ! IELset
    !%OMP END DO
    
    ! Release memory
    call elprep_releaseElementSet(revalElementSet)

    if (.not. rform%ballCoeffConstant) then
      deallocate(Dcoefficients)
    end if
    deallocate(IdofsTrial)
    deallocate(IdofsTest)
    deallocate(DbasTrial)
    deallocate(DbasTest)
    deallocate(Kentry)
    deallocate(Dentry)

    !%OMP END PARALLEL

    deallocate(p_DcubPtsRef)
    deallocate(Domega)

  end do ! icurrentElementDistr

  ! Finish
  
  end subroutine

  !****************************************************************************

  !<subroutine>

  subroutine bilf_getLocalMatrixIndices (rmatrix,Irows,Icolumns,Kentry)
  
  !<description>
  
  ! Calculates index positions of local matrices in a global matrix.
  ! For a set of elements, Icolumns and Irows define the row and column indices
  ! of local matrices which have to be accessed in a global matrix rmatrix.
  ! The routine then calculates the positions of all the corresponding matrix
  ! entries in the data array of the matrix rmatrix and saves the result
  ! to the array Kentry.
  !
  ! IMPORTANT: For performance reasons, the columns and rows in Kentry
  ! are saved *transposed* in comparison to the matrix rmatrix!
  ! That means that
  !    Kentry(j,i,:) = position of element (Irows(i,:),Icolumns(j,:))
  ! holds!
  
  !</description>
  
  !<input>
  
  ! The global matrix which has to be accessed.
  type(t_matrixScalar), intent(in) :: rmatrix
  
  ! Array identifying all rows in the global matrix which have to be
  ! accessed.
  ! DIMENSION(#rows per element, #elements).
  integer, dimension(:,:), intent(in) :: Irows

  ! Array identifying all columns in the global matrix which have to be
  ! accessed.
  ! DIMENSION(#columns per element, #elements).
  integer, dimension(:,:), intent(in) :: Icolumns
  
  !</input>
  
  !<output>
  
  ! Array receiving the positions of the local matrices in the global matrix.
  ! DIMENSION(#columns per element,#rows per element,#elements).
  ! Saved in a transposed way:
  !    Kentry(j,i,:) = position of element (Irows(i,:),Icolumns(j,:))
  integer, dimension(:,:,:), intent(out) :: Kentry
  
  !</output>
  
  !</subroutine>
  
    ! local variables
    integer, dimension(:), pointer :: p_Kcol, p_Kld
    integer :: na,iel,idofe,jdofe,indofTest,indofTrial,nelements,jcol0,jdfg,jcol

    indofTrial = ubound(Icolumns,1)
    indofTest = ubound(Irows,1)
    nelements = ubound(Irows,2)

    select case (rmatrix%cmatrixFormat)
    case (LSYSSC_MATRIX1)
    
      ! That's easy, we can directly calculate the positions
      do iel = 1,nelements
        do idofe=1,indofTest
          do jdofe=1,indofTrial
            Kentry(jdofe,idofe,iel) = &
                Irows(idofe,iel) * rmatrix%NCOLS + Icolumns(jdofe,iel)
          end do
        end do
      end do
      
    case (LSYSSC_MATRIX7,LSYSSC_MATRIX9)
    
      ! Get pointers to the row/column structure of the matrix
      call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
      call lsyssc_getbase_Kld (rmatrix,p_Kld)
      na = rmatrix%NA

      ! We build a quadratic indofTrial*indofTest local matrix:
      ! Kentry(1..indofTrial,1..indofTest) receives the position 
      !   in the global system matrix
      !
      ! Loop through elements in the set and for each element,
      ! loop through the local matrices to initialise them:
      do iel=1,nelements
      
        ! For building the local matrices, we have first to
        ! loop through the test functions (the "O"'s), as these
        ! define the rows in the matrix.
        do idofe=1,indofTest
        
          ! Row IDOFE of the local matrix corresponds 
          ! to row=global DOF KDFG(IDOFE) in the global matrix.
          ! This is one of the the "O"'s in the above picture.
          ! Get the starting position of the corresponding row
          ! to JCOL0:

          jcol0=p_Kld(Irows(idofe,iel))
          
          ! Now we loop through the other DOF's on the current element
          ! (the "O"'s).
          ! All these have common support with our current basis function
          ! and will therefore give an additive value to the global
          ! matrix.
          
          do jdofe=1,indofTrial
            
            ! Get the global DOF of the "X" which interacts with 
            ! our "O".
            
            jdfg=Icolumns(jdofe,iel)
            
            ! Starting in JCOL0 (which points to the beginning of
            ! the line initially), loop through the elements in
            ! the row to find the position of column IDFG.
            ! Jump out of the DO loop if we find the column.
            
            do jcol=jcol0,na
              if (p_Kcol(jcol) .eq. jdfg) exit
            end do

            ! Because columns in the global matrix are sorted 
            ! ascendingly (except for the diagonal element),
            ! the next search can start after the column we just found.
            
            ! JCOL0=JCOL+1
            
            ! Save the position of the matrix entry into the local
            ! matrix.
            ! Note that a column in Kentry corresponds to a row in
            ! the real matrix. We aligned Kentry this way to get
            ! higher speed of the assembly routine, since this leads
            ! to better data locality.
            
            Kentry(jdofe,idofe,iel)=jcol
            
          end do ! IDOFE
          
        end do ! JDOFE
        
      end do ! IEL
      
    end select
      
  end subroutine

end module
