!##############################################################################
!# ****************************************************************************
!# <name> multileveloperators </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines for the creation of extended multi-level
!# operators, e.g. Matrix-Creation for L2-Projection-operators.
!#
!# This module contains the following routines:
!#
!#  1.) mlop_create2LvlMatrixStruct
!#      -> Creates a scalar matrix structure in a specified matrix format
!#         for multi-level-based operations between two discretisations.
!#
!#  2.) mlop_build2LvlMassMatrix
!#      -> Assembles the entries of a 2-level-mass matrix that is needed for
!#         L2-Prolongation and restriction operators.
!#
!#  3.) mlop_build2LvlProlMatrix
!#      -> Assembles the entries of a prolongation matrix based on element-wise
!#         L2-projection.
!#
!#  4.) mlop_build2LvlInterpMatrix
!#      -> Assembles the entries of a interpolation matrix based on element-
!#         wise L2-projection.
!#
!# 5.) mlop_initPerfConfig
!#      -> Initialises the global performance configuration
!#
!# ---------------------------------------------------------------------------- \\
!# A small note on conformal discretisations \\
!# ----------------------------------------------------------------------------
!#
!# Basically, all routines implemented in this module are able to handle
!# conformal discretisations - under a few conditions:
!#
!# 1. Both the coarse and fine mesh discretisation must have the same number
!#    of FE spaces ( = number of element distributions).
!#
!# 2. If a coarse mesh element IELC, which was refined into the elements
!#    IELF_1, ..., IELF_k in the fine mesh, belongs to the element
!#    distribution i of the coarse mesh discretisation, then all its children
!#    IELF_1, ..., IELF_k must (also) belong to the element distribution i
!#    of the fine mesh discretisation.
!#
!# Please note that the above conditions are SILENTLY ASSUMED to be fulfilled.
!#
!# One of the most interesting special cases of conformal discretisations are
!# discretisations, which define "the same" FE space on differently shaped
!# elements, e.g. P1/Q1 for a mixed triangle/quadrilateral mesh. So if every
!# triangle belongs to the P1 space, and every quadrilateral belongs to the Q1
!# space on the coarse mesh, then the above two conditions are fulfilled, as
!# during refinement every triangle produces only triangles and dito for
!# quadrilaterals.
!#
!# Note:
!# The above conditions where not chosen arbitrarily to annoy application
!# programmers. If these conditions are fulfilled, then the assembly of both
!# the matrix structure and its entries can be performed patch-wise and
!# therefore relatively efficient. If (especially) the second condition was
!# violated, then the assembly would need to work element-wise and perform
!# a lot more technically non-trivial work to perform its task!
!# </purpose>
!##############################################################################

module multileveloperators

  use basicgeometry
  use collection
  use cubature
  use derivatives
  use dofmapping
  use domainintegration
  use element
  use elementpreprocessing
  use fsystem
  use genoutput
  use linearalgebra
  use linearsystemscalar
  use perfconfig
  use scalarpde
  use spatialdiscretisation
  use storage
  use transformation
  use triangulation
  
  implicit none
  
  private
  
  public :: mlop_initPerfConfig
  public :: mlop_create2LvlMatrixStruct
  public :: mlop_build2LvlMassMatrix
  public :: mlop_build2LvlProlMatrix
  public :: mlop_build2LvlInterpMatrix

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

!</types>

!<constants>

!<constantblock description="Averaging specifiers for prolongation matrix assembly">

  ! Use arithmetic averaging
  integer, parameter, public :: MLOP_AVRG_ARITHMETIC = 1
  
  ! Average by basis function L2-mass
  integer, parameter, public :: MLOP_AVRG_MASS = 2
  
  ! Average by inverse basis function L2-mass
  integer, parameter, public :: MLOP_AVRG_INV_MASS = 3

!</constantblock>

!<constantblock description="Constants defining the blocking of the assembly">

  ! *** LEGACY CONSTANT, use the more flexible performance configuration ***
  ! Number of elements to handle simultaneously when building vectors
#ifndef LINF_NELEMSIM
  integer, parameter, public :: MLOP_NELEMSIM = 100
#endif
  
!</constantblock>

!</constants>

  !************************************************************************

  ! global performance configuration
  type(t_perfconfig), target, save :: mlop_perfconfig

  !************************************************************************

contains
  
  !****************************************************************************

!<subroutine>

  subroutine mlop_initPerfConfig(rperfconfig)

!<description>
  ! This routine initialises the global performance configuration
!</description>

!<input>
  ! OPTIONAL: performance configuration that should be used to initialise
  ! the global performance configuration. If not present, the values of
  ! the legacy constants is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>
!</subroutine>

    if (present(rperfconfig)) then
      mlop_perfconfig = rperfconfig
    else
      call pcfg_initPerfConfig(mlop_perfconfig)
      mlop_perfconfig%NELEMSIM = MLOP_NELEMSIM
    end if
  
  end subroutine mlop_initPerfConfig

  !****************************************************************************

!<subroutine>

  subroutine mlop_create2LvlMatrixStruct (rdiscretisationCoarse,&
                      rdiscretisationFine, iformat, rmatrixScalar,&
                      imemguess, rperfconfig)
  
!<description>
  ! This routine allows to calculate the structure of a finite-element matrix
  ! on the heap. The size of the matrix is determined dynamically.
!</description>

!<input>
  ! The underlying discretisation structure defined on the coarse mesh which
  ! is to be used to create the matrix.
  type(t_spatialDiscretisation), intent(in), target :: rdiscretisationCoarse
  
  ! The underlying discretisation structure defined on the fine mesh which
  ! is to be used to create the matrix.
  type(t_spatialDiscretisation), intent(in), target :: rdiscretisationFine

  ! Format of the matrix structure to be created. One of the LSYSSC_xxxx
  ! constants.
  integer, intent(in) :: iformat

  ! OPTIONAL: An initial guess about how much memory the matrix needs. If set
  ! to 0 or not given, an initial guess of 16*NEQ (but at least 10000 matrix
  ! entries) is assumed.
  integer, intent(in), optional :: imemGuess
  
  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<output>
  ! The structure of a scalar matrix, fitting to the given discretisation
  ! combination. Memory fo the structure is allocated dynamically on the heap.
  type(t_matrixScalar), intent(out) :: rmatrixScalar
!</output>

!</subroutine>

  ! local variables
  integer :: imem
  
    imem = 0
    if (present(imemguess)) then
      imem = max(0,imemguess)
    end if
    
    ! Let us make sure that the discretisations are uniform or conformal.
    if (((rdiscretisationCoarse%ccomplexity .ne. SPDISC_CONFORMAL) .or. &
         (rdiscretisationFine%ccomplexity .ne. SPDISC_CONFORMAL)) .and. &
        ((rdiscretisationCoarse%ccomplexity .ne. SPDISC_UNIFORM) .or. &
         (rdiscretisationFine%ccomplexity .ne. SPDISC_UNIFORM))) then
        
      call output_line ('Discretisations must be uniform or conformal!', &
          OU_CLASS_ERROR,OU_MODE_STD,'mlop_create2LvlMatrixStruct')
      call sys_halt()
      
    end if
    
    ! In the case the discretisations are conformal, they must have the same
    ! number of FE spaces.
    if(rdiscretisationCoarse%inumFESpaces .ne. rdiscretisationFine%inumFESpaces) then

      call output_line ('Discretisations must have same number of FE spaces!', &
          OU_CLASS_ERROR,OU_MODE_STD,'mlop_create2LvlMatrixStruct')
      call sys_halt()

    end if
    
    ! Which matrix structure do we have to create?
    select case (iformat)
    
    case (LSYSSC_MATRIX9)
    
      ! Call the creation routine for structure 9:
      call mlop_create2LvlMatStruct9_conf (rdiscretisationCoarse,&
                           rdiscretisationFine,rmatrixScalar,imem,rperfconfig)
     
    case (LSYSSC_MATRIX7)
    
      ! Call the creation routine for structure 9:
      call mlop_create2LvlMatStruct9_conf (rdiscretisationCoarse,&
                           rdiscretisationFine,rmatrixScalar,imem,rperfconfig)

      ! Translate to matrix structure 7:
      call lsyssc_convertMatrix (rmatrixScalar,LSYSSC_MATRIX7)
      
    case default
      call output_line ('Not supported matrix structure!', &
          OU_CLASS_ERROR,OU_MODE_STD,'mlop_create2LvlMatrixStruct')
      call sys_halt()
      
    end select

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine mlop_build2LvlMassMatrix (rdiscretisationCoarse,&
                      rdiscretisationFine,bclear,rmatrixScalar,rperfconfig)
  
!<description>
  ! This routine calculates the entries of a 2-Level mass matrix.
  ! The matrix structure must have been initialised by the
  ! mlop_create2LvlMatrixStruct routine before calling this function.
!</description>

!<input>
  ! The underlying discretisation structure defined on the coarse mesh which
  ! is to be used to create the matrix.
  type(t_spatialDiscretisation), intent(in), target :: rdiscretisationCoarse
  
  ! The underlying discretisation structure defined on the fine mesh which
  ! is to be used to create the matrix.
  type(t_spatialDiscretisation), intent(in), target :: rdiscretisationFine
  
  ! Whether to clear the matrix before calculating the entries.
  ! If .FALSE., the new matrix entries are added to the existing entries.
  logical, intent(in) :: bclear
  
  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! The FE matrix. Calculated matrix entries are imposed to this matrix.
  type(t_matrixScalar), intent(inout) :: rmatrixScalar
!</inputoutput>

!</subroutine>

  ! local variables
  type(t_matrixScalar) :: rmatrixBackup
  
    ! The matrix must be unsorted, otherwise we can not set up the matrix.
    ! Note that we cannot switch off the sorting as easy as in the case
    ! of a vector, since there is a structure behind the matrix! So the caller
    ! has to make sure, the matrix is unsorted when this routine is called.
    if (rmatrixScalar%isortStrategy .gt. 0) then
      call output_line ('Matrix-structure must be unsorted!', &
          OU_CLASS_ERROR,OU_MODE_STD,'mlop_build2LvlMassMatrix')
      call sys_halt()
    end if

    ! Let us make sure that the discretisations are uniform or conformal.
    if (((rdiscretisationCoarse%ccomplexity .ne. SPDISC_CONFORMAL) .or. &
         (rdiscretisationFine%ccomplexity .ne. SPDISC_CONFORMAL)) .and. &
        ((rdiscretisationCoarse%ccomplexity .ne. SPDISC_UNIFORM) .or. &
         (rdiscretisationFine%ccomplexity .ne. SPDISC_UNIFORM))) then
        
      call output_line ('Discretisations must be uniform or conformal!', &
          OU_CLASS_ERROR,OU_MODE_STD,'mlop_build2LvlMassMatrix')
      call sys_halt()
      
    end if
    
    ! In the case the discretisations are conformal, they must have the same
    ! number of FE spaces.
    if(rdiscretisationCoarse%inumFESpaces .ne. rdiscretisationFine%inumFESpaces) then

      call output_line ('Discretisations must have same number of FE spaces!', &
          OU_CLASS_ERROR,OU_MODE_STD,'mlop_build2LvlMassMatrix')
      call sys_halt()

    end if

    ! Which matrix structure do we have?
    select case (rmatrixScalar%cmatrixFormat)
    case (LSYSSC_MATRIX9)
      call mlop_build2LvlMass9_conf (rdiscretisationCoarse,&
                      rdiscretisationFine,bclear,rmatrixScalar,rperfconfig)
    
    case (LSYSSC_MATRIX7)
      ! Convert structure 7 to structure 9.For that purpose, make a backup of
      ! the original matrix...
      call lsyssc_duplicateMatrix (rmatrixScalar,rmatrixBackup,&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
          
      ! Convert the matrix
      call lsyssc_convertMatrix (rmatrixBackup,LSYSSC_MATRIX9)
      
      ! Create the matrix in structure 9
      call mlop_build2LvlMass9_conf (rdiscretisationCoarse,&
                      rdiscretisationFine,bclear,rmatrixScalar,rperfconfig)
                                     
      ! Convert back to structure 7
      call lsyssc_convertMatrix (rmatrixBackup,LSYSSC_MATRIX7)
      
      ! Copy the entries back to the original matrix and release memory.
      call lsyssc_duplicateMatrix (rmatrixBackup,rmatrixScalar,&
          LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPYOVERWRITE)
          
      call lsyssc_releaseMatrix (rmatrixBackup)
                                     
    case default
      call output_line ('Not supported matrix structure!', &
          OU_CLASS_ERROR,OU_MODE_STD,'mlop_build2LvlMassMatrix')
      call sys_halt()
      
    end select


  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine mlop_build2LvlProlMatrix (rdiscretisationCoarse,&
                      rdiscretisationFine,bclear,rmatrixScalar,&
                      cavrgType,rperfconfig)
  
!<description>
  ! This routine calculates the entries of a 2-Level prolongation matrix.
  ! The matrix structure must have been initialised by the
  ! mlop_create2LvlMatrixStruct routine before calling this function.
!</description>

!<input>
  ! The underlying discretisation structure defined on the coarse mesh which
  ! is to be used to create the matrix.
  type(t_spatialDiscretisation), intent(in), target :: rdiscretisationCoarse
  
  ! The underlying discretisation structure defined on the fine mesh which
  ! is to be used to create the matrix.
  type(t_spatialDiscretisation), intent(in), target :: rdiscretisationFine
  
  ! Whether to clear the matrix before calculating the entries.
  ! If .FALSE., the new matrix entries are added to the existing entries.
  logical, intent(in) :: bclear
  
  ! OPTIONAL: Specifies which type of averaging is to be used.
  ! One of the MLOP_AVRG_XXXX constants defined above. If not given,
  ! MLOP_AVRG_ARITHMETIC is used.
  integer, optional, intent(in) :: cavrgType

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! The FE matrix. Calculated matrix entries are imposed to this matrix.
  type(t_matrixScalar), intent(inout) :: rmatrixScalar
!</inputoutput>

!</subroutine>

  ! local variables
  type(t_matrixScalar) :: rmatrixBackup
  integer :: caverage = MLOP_AVRG_ARITHMETIC
  
    if(present(cavrgType)) caverage = cavrgType

    ! The matrix must be unsorted, otherwise we can not set up the matrix.
    ! Note that we cannot switch off the sorting as easy as in the case
    ! of a vector, since there is a structure behind the matrix! So the caller
    ! has to make sure, the matrix is unsorted when this routine is called.
    if (rmatrixScalar%isortStrategy .gt. 0) then
      call output_line ('Matrix-structure must be unsorted!', &
          OU_CLASS_ERROR,OU_MODE_STD,'mlop_build2LvlProlMatrix')
      call sys_halt()
    end if

    ! Let us make sure that the discretisations are uniform or conformal.
    if (((rdiscretisationCoarse%ccomplexity .ne. SPDISC_CONFORMAL) .or. &
         (rdiscretisationFine%ccomplexity .ne. SPDISC_CONFORMAL)) .and. &
        ((rdiscretisationCoarse%ccomplexity .ne. SPDISC_UNIFORM) .or. &
         (rdiscretisationFine%ccomplexity .ne. SPDISC_UNIFORM))) then
        
      call output_line ('Discretisations must be uniform or conformal!', &
          OU_CLASS_ERROR,OU_MODE_STD,'mlop_build2LvlProlMatrix')
      call sys_halt()
      
    end if
    
    ! In the case the discretisations are conformal, they must have the same
    ! number of FE spaces.
    if(rdiscretisationCoarse%inumFESpaces .ne. rdiscretisationFine%inumFESpaces) then

      call output_line ('Discretisations must have same number of FE spaces!', &
          OU_CLASS_ERROR,OU_MODE_STD,'mlop_build2LvlProlMatrix')
      call sys_halt()

    end if

    ! Which matrix structure do we have?
    select case (rmatrixScalar%cmatrixFormat)
    case (LSYSSC_MATRIX9)
      call mlop_build2LvlProl9_conf (rdiscretisationCoarse,&
                      rdiscretisationFine,bclear,rmatrixScalar,&
                      caverage,rperfconfig)
    
    case (LSYSSC_MATRIX7)
      ! Convert structure 7 to structure 9.For that purpose, make a backup of
      ! the original matrix...
      call lsyssc_duplicateMatrix (rmatrixScalar,rmatrixBackup,&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
          
      ! Convert the matrix
      call lsyssc_convertMatrix (rmatrixBackup,LSYSSC_MATRIX9)
      
      ! Create the matrix in structure 9
      call mlop_build2LvlProl9_conf (rdiscretisationCoarse,&
                      rdiscretisationFine,bclear,rmatrixScalar,&
                      caverage,rperfconfig)
                                     
      ! Convert back to structure 7
      call lsyssc_convertMatrix (rmatrixBackup,LSYSSC_MATRIX7)
      
      ! Copy the entries back to the original matrix and release memory.
      call lsyssc_duplicateMatrix (rmatrixBackup,rmatrixScalar,&
          LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPYOVERWRITE)
          
      call lsyssc_releaseMatrix (rmatrixBackup)
                                     
    case default
      call output_line ('Not supported matrix structure!', &
          OU_CLASS_ERROR,OU_MODE_STD,'mlop_build2LvlProlMatrix')
      call sys_halt()
      
    end select
  
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine mlop_build2LvlInterpMatrix (rdiscretisationCoarse,&
                      rdiscretisationFine,bclear,rmatrixScalar,&
                      cavrgType,rperfconfig)
  
!<description>
  ! This routine calculates the entries of a 2-Level interpolation matrix.
  ! The matrix structure must have been initialised by the
  ! mlop_create2LvlMatrixStruct routine before calling this function.
!</description>

!<input>
  ! The underlying discretisation structure defined on the coarse mesh which
  ! is to be used to create the matrix.
  type(t_spatialDiscretisation), intent(in), target :: rdiscretisationCoarse
  
  ! The underlying discretisation structure defined on the fine mesh which
  ! is to be used to create the matrix.
  type(t_spatialDiscretisation), intent(in), target :: rdiscretisationFine
  
  ! Whether to clear the matrix before calculating the entries.
  ! If .FALSE., the new matrix entries are added to the existing entries.
  logical, intent(in) :: bclear
  
  ! OPTIONAL: Specifies which type of averaging is to be used.
  ! One of the MLOP_AVRG_XXXX constants defined above. If not given,
  ! MLOP_AVRG_ARITHMETIC is used.
  integer, optional, intent(in) :: cavrgType

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! The interpolation matrix that is to be assembled.
  type(t_matrixScalar), intent(inout) :: rmatrixScalar
!</inputoutput>

!</subroutine>

  ! local variables
  type(t_matrixScalar) :: rmatrixBackup
  integer :: caverage = MLOP_AVRG_ARITHMETIC
  
    if(present(cavrgType)) caverage = cavrgType

    ! The matrix must be unsorted, otherwise we can not set up the matrix.
    ! Note that we cannot switch off the sorting as easy as in the case
    ! of a vector, since there is a structure behind the matrix! So the caller
    ! has to make sure, the matrix is unsorted when this routine is called.
    if (rmatrixScalar%isortStrategy .gt. 0) then
      call output_line ('Matrix-structure must be unsorted!', &
          OU_CLASS_ERROR,OU_MODE_STD,'mlop_build2LvlInterpMatrix')
      call sys_halt()
    end if

    ! Let us make sure that the discretisations are uniform or conformal.
    if (((rdiscretisationCoarse%ccomplexity .ne. SPDISC_CONFORMAL) .or. &
         (rdiscretisationFine%ccomplexity .ne. SPDISC_CONFORMAL)) .and. &
        ((rdiscretisationCoarse%ccomplexity .ne. SPDISC_UNIFORM) .or. &
         (rdiscretisationFine%ccomplexity .ne. SPDISC_UNIFORM))) then
        
      call output_line ('Discretisations must be uniform or conformal!', &
          OU_CLASS_ERROR,OU_MODE_STD,'mlop_build2LvlInterpMatrix')
      call sys_halt()
      
    end if
    
    ! In the case the discretisations are conformal, they must have the same
    ! number of FE spaces.
    if(rdiscretisationCoarse%inumFESpaces .ne. rdiscretisationFine%inumFESpaces) then

      call output_line ('Discretisations must have same number of FE spaces!', &
          OU_CLASS_ERROR,OU_MODE_STD,'mlop_build2LvlInterpMatrix')
      call sys_halt()

    end if

    ! Which matrix structure do we have?
    select case (rmatrixScalar%cmatrixFormat)
    case (LSYSSC_MATRIX9)
      call mlop_build2LvlInterp9_conf (rdiscretisationCoarse,&
                      rdiscretisationFine,bclear,rmatrixScalar,&
                      caverage,rperfconfig)
    
    case (LSYSSC_MATRIX7)
      ! Convert structure 7 to structure 9.For that purpose, make a backup of
      ! the original matrix...
      call lsyssc_duplicateMatrix (rmatrixScalar,rmatrixBackup,&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
          
      ! Convert the matrix
      call lsyssc_convertMatrix (rmatrixBackup,LSYSSC_MATRIX9)
      
      ! Create the matrix in structure 9
      call mlop_build2LvlInterp9_conf (rdiscretisationCoarse,&
                      rdiscretisationFine,bclear,rmatrixScalar,&
                      caverage,rperfconfig)
                                     
      ! Convert back to structure 7
      call lsyssc_convertMatrix (rmatrixBackup,LSYSSC_MATRIX7)
      
      ! Copy the entries back to the original matrix and release memory.
      call lsyssc_duplicateMatrix (rmatrixBackup,rmatrixScalar,&
          LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPYOVERWRITE)
          
      call lsyssc_releaseMatrix (rmatrixBackup)
                                     
    case default
      call output_line ('Not supported matrix structure!', &
          OU_CLASS_ERROR,OU_MODE_STD,'mlop_build2LvlInterpMatrix')
      call sys_halt()
      
    end select
  
  end subroutine
    
  !****************************************************************************
  
!<subroutine>
  
  subroutine mlop_create2LvlMatStruct9_conf (rdiscrCoarse,rdiscrFine,&
                                             rmatrixScalar,imemGuess,rperfconfig)
  
!<description>
  ! This routine creates according to a given discretisation the matrix
  ! structure of a structure-9 matrix. The discretisation is assumed to be
  ! conformal, i.e. the DOF`s of different FE spaces in the trial space
  ! fit together. The function space for trial and test functions
  ! may be different.
!</description>

!<input>
  
  ! The underlying discretisation structures which are to be used to
  ! create the matrix.
  type(t_spatialDiscretisation), intent(in), target :: rdiscrCoarse
  type(t_spatialDiscretisation), intent(in), target :: rdiscrFine
  
  ! An initial guess about how much memory the matrix needs. If set to 0,
  ! an initial guess of 16*NEQ (but at least 10000 matrix entries) is assumed.
  integer, intent(in) :: imemGuess
  
  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<output>
  ! The structure of a scalar matrix, fitting to the given discretisation.
  ! Memory fo rthe structure is allocated dynamically on the heap.
  type(t_matrixScalar), intent(out) :: rmatrixScalar
!</output>

!</subroutine>

  ! local variables
  integer :: NEQ, IEQ, IROW, JCOL, IPOS, istartIdx, NA, nmaxCol
  integer :: IDOFE, JDOFE, i, IHELP, IELDIST
  integer :: IELC,IELF,IELIDX,NELREF
  logical :: BSORT
  
  ! Pointer to the performance configuration
  type(t_perfconfig), pointer :: p_rperfconfig

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
  
  ! Number of elements in the current coarse/fine mesh element distribution
  integer :: NELC, NELF
  
  ! Size of memory blocks
  integer :: imemblkSize
  
  ! Blocksize in terms of NEQ for guessing memory.
  ! The initial guess for memory is iblkSize*iblkSize*NEQ and every time
  ! memory is needed, another iblkSize*NEQ elements are added.
  integer, parameter :: iblkSize = 4
  
  ! Number of memory blocks to allocate
  integer, parameter :: NmemBlkCount = 5

  ! Pointer to KLD, KCOL
  integer, dimension(:), pointer :: p_KLD, p_KCOL
  
  ! Size of memory currently allocated
  integer :: iallocated
  
  ! An allocateable array accepting the DOF`s of a set of elements.
  integer, dimension(:,:), pointer :: p_IdofsCoarse, p_IdofsFine
  
  ! Number of local degees of freedom for trial and test functions
  integer :: indofCoarse, indofFine, inmaxdofCoarse, inmaxdofFine
  
  ! The triangulation structure - to shorten some things...
  type(t_triangulation), pointer :: p_rtriaCoarse, p_rtriaFine
  
  ! A pointer to an element-number list
  integer, dimension(:), pointer :: p_IelementList,p_IelementRef
  
  ! Current element distribution
  type(t_elementDistribution), pointer :: p_relementDistribution

  ! Number of elements that have already been processed and number of
  ! elements that are to be processed in the current run
  integer :: nelementsDone, nelementsToDo
  
  ! Number of elements that are to be processed at once
  integer :: nelementsCoarse, nelementsFine, nmaxelementsCoarse, nmaxelementsFine
    
  ! Two arrays for the refinement-patch arrays of the coarse triangulation
  integer, dimension(:), pointer :: p_IrefPatchIdx, p_IrefPatch

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => mlop_perfconfig
    end if

    ! The algorithm is: Test every DOF on one element against each other
    ! DOF on the same element and save the combination into a matrix
    ! in structure 9!
    !
    ! At first, initialise the structure-9 matrix:
    
    !rmatrixScalar%p_rspatialDiscretisation => rdiscretisation
    rmatrixScalar%cmatrixFormat = LSYSSC_MATRIX9
    
    ! Get the #DOF`s of the test space - as #DOF`s of the test space is
    ! the number of equations in our matrix. The #DOF`s in the trial space
    ! gives the number of columns of our matrix.
    rmatrixScalar%NCOLS         = dof_igetNDofGlob(rdiscrCoarse)
    rmatrixScalar%NEQ           = dof_igetNDofGlob(rdiscrFine)
    
    ! and get a pointer to the triangulation.
    p_rtriaCoarse => rdiscrCoarse%p_rtriangulation
    p_rtriaFine => rdiscrFine%p_rtriangulation
    
    ! Get the refinement patch arrays from the fine triangulation
    call storage_getbase_int(p_rtriaFine%h_IrefinementPatchIdx, p_IrefPatchIdx)
    call storage_getbase_int(p_rtriaFine%h_IrefinementPatch, p_IrefPatch)
    
    ! Get NEQ - we need it for guessing memory...
    NEQ = rmatrixScalar%NEQ
    
    if (NEQ .eq. 0) then
      call output_line ('Empty matrix!', &
          OU_CLASS_ERROR,OU_MODE_STD,'mlop_create2LvlMatStruct9_conf')
      call sys_halt()
    end if
    
    ! Allocate KLD...
    call storage_new ('mlop_create2LvlMatStruct9_conf', 'KLD', &
                        NEQ+1, ST_INT, rmatrixScalar%h_KLD, ST_NEWBLOCK_NOINIT)
                        
    ! And allocate Kdiagonal - although it is not needed, it has to allocated.
    call storage_new ('mlop_create2LvlMatStruct9_conf', 'KLD', &
                        NEQ, ST_INT, rmatrixScalar%h_Kdiagonal, ST_NEWBLOCK_NOINIT)

    ! This must be a storage_getbase, no lsyssc_getbase, since this is the
    ! matrix construction routine!
    call storage_getbase_int(rmatrixScalar%h_Kld,p_KLD)
    
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
    call storage_new ('mlop_create2LvlMatStruct9_conf', 'Ihicol', &
                        p_Isize(1), ST_INT, p_Ihcol(1), ST_NEWBLOCK_NOINIT)
    call storage_getbase_int (p_Ihcol(1),p_Icol)

    ! The new index array must be filled with 0 - otherwise
    ! the search routine below will not work!
    call storage_new ('mlop_create2LvlMatStruct9_conf', 'p_Ihindx', &
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
    ! - "2" (the column number is added to the end of Icol, i.e.
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
    
    ! Loop through all element distributions and determine maximum values
    inmaxdofCoarse = 0
    inmaxdofFine = 0
    nmaxelementsCoarse = 0
    nmaxelementsFine = 0
    do i = 1, rdiscrCoarse%inumFESpaces
    
      ! Activate the current coarse mesh element distribution
      p_relementDistribution => rdiscrCoarse%RelementDistr(i)

      ! Get the number of local DOF`s for trial and test functions
      indofCoarse = elem_igetNDofLoc(rdiscrCoarse%RelementDistr(i)%celement)
      indofFine = elem_igetNDofLoc(rdiscrFine%RelementDistr(i)%celement)
      
      ! Calculate the number of coarse mesh elements we want to process
      ! in one run.
      nelementsCoarse = min(p_rperfconfig%NELEMSIM,p_relementDistribution%NEL)
      
      ! Now calculate the number of fine mesh elements we want to process
      ! in one run.
      select case(p_rtriaFine%ndim)
      case (1)
        nelementsFine = 2*nelementsCoarse + 10
      case (2)
        nelementsFine = 4*nelementsCoarse + 20
      case (3)
        nelementsFine = 8*nelementsCoarse + 40
      end select
      
      ! Determine maximum values
      inmaxdofCoarse = max(inmaxdofCoarse, indofCoarse)
      inmaxdofFine   = max(inmaxdofFine,   indofFine)
      nmaxelementsCoarse = max(nmaxelementsCoarse, nelementsCoarse)
      nmaxelementsFine   = max(nmaxelementsFine,   nelementsFine)
    
    end do
    
    ! Allocate an array saving a couple of DOF`s for trial and test functions
    allocate(p_IdofsCoarse(inmaxdofCoarse,nmaxelementsCoarse))
    allocate(p_IdofsFine(inmaxdofFine,nmaxelementsFine))
    
    ! And allocate the refinemed element list for the test functions
    allocate(p_IelementRef(nelementsFine))
    
    ! Now let us loop over all element distributions
    do IELDIST = 1, rdiscrCoarse%inumFESpaces

      ! Activate the current coarse mesh element distribution
      p_relementDistribution => rdiscrCoarse%RelementDistr(IELDIST)
      
      if (p_relementDistribution%NEL .eq. 0) cycle
    
      ! Get the number of local DOF`s for trial and test functions
      indofCoarse = elem_igetNDofLoc(p_relementDistribution%celement)
      indofFine = elem_igetNDofLoc(rdiscrFine%RelementDistr(IELDIST)%celement)
      
      ! Calculate the number of coarse mesh elements we want to process
      ! in one run.
      nelementsCoarse = min(p_rperfconfig%NELEMSIM,p_relementDistribution%NEL)
      
      ! p_IelementList must point to our set of elements in the discretisation
      ! with that the trial functions
      call storage_getbase_int (p_relementDistribution%h_IelementList, &
                                p_IelementList)

      ! Get the number of coarse mesh elements there.
      NELC = p_relementDistribution%NEL

      ! Set the pointers/indices to the initial position. During the
      ! search for new DOF`s, these might be changed if there is not enough
      ! memory in the first block.
      icurrentblock = 1
      istartidx = 0
      p_Icol => Rmemblock(1)%p_Icol
      p_Iindx => Rmemblock(1)%p_Iindx
      
      ! Loop over the elements.
      nelementsDone = 0
      do while(nelementsDone .lt. NELC)
      
        ! We always try to handle nelementsTrial elements simultaneously.
        ! Of course, we will not handle more elements than the coarse
        ! mesh discretisation has.
        nelementsToDo = min(NELC-nelementsDone, nelementsCoarse)
        
        ! Now comes the interesting part - we have to ensure that the DOF-mapping
        ! of the fine mesh discretisation fits into our DOF-array.
        ! If, for example, a coarse mesh quad was refined into more than 4 fine
        ! mesh quads, then it might happen that we cannot handle nelementsToDo
        ! coarse mesh elements at once, but we need to decrease nelementsToDo.
        
        NELF = 0
        do IELC = 1, nelementsToDo
        
          ! Get the index of the coarse mesh element
          IELIDX = p_IelementList(nelementsDone+IELC)
          
          ! Get the number of fine mesh elements that have been refined from the
          ! currently processed coarse mesh element.
          NELREF = p_IrefPatchIdx(IELIDX+1) - p_IrefPatchIdx(IELIDX)
          
          ! Now if (NELF+NELREF) is greater than nelementsTest, then we need
          ! to decrease nelementsToDo and exit the loop...
          ! This case should never happen if the coarse mesh was refined using
          ! the 2-Level-Ordering algorithm, but might happen for more freaky
          ! refinement techniques...
          if((NELF+NELREF) .gt. nelementsFine) then
            nelementsToDo = IELC-1
            exit
          end if
          
          ! Copy the indices of the elements into the element list for the
          ! fine mesh discretisation
          do IELF = 1, NELREF
            p_IelementRef(NELF+IELF) = p_IrefPatch(p_IrefPatchIdx(IELIDX)+IELF-1)
          end do
          
          ! Add the number of refined elements to the counter
          NELF = NELF + NELREF
        
        end do
        
        ! If nelementsToDo is 0, then we have a serious problem...
        if (nelementsToDo .le. 0) then
          call output_line ('INTERNAL ERROR: nelementsToDo = 0!', &
              OU_CLASS_ERROR,OU_MODE_STD,'mlop_create2LvlMatStruct9_conf')
          call sys_halt()
        end if
        
        ! Call the DOF-mapping routine for the coarse and fine mesh
        call dof_locGlobMapping_mult(rdiscrCoarse, &
            p_IelementList(nelementsDone+1:nelementsDone+nelementsToDo), &
            p_IdofsCoarse)
        call dof_locGlobMapping_mult(rdiscrFine, p_IelementRef(1:NELF), &
            p_IdofsFine)
        
        ! Reset the counter
        NELF = 0
        
        ! Loop through all the elements of the coarse mesh set
        do IELC = 1, nelementsToDo
        
          ! Get the index of the currently processed coarse mesh element
          IELIDX = p_IelementList(nelementsDone+IELC)
          
          ! Get the number of fine mesh elements that are refined from the
          ! current coarse mesh element
          NELREF = p_IrefPatchIdx(IELIDX+1) - p_IrefPatchIdx(IELIDX)

          ! And loop through all elements of the current refinement patch
          do IELF = 1, NELREF
          
            ! For building the local matrices, we have first to
            ! loop through the test functions (the "O"`s), as these
            ! define the rows in the matrix.
            do IDOFE=1,indofFine

              ! The DOF IDOFE is now our "O".
              ! This global DOF gives us the row we have to build.
              IROW = p_IdofsFine(IDOFE,NELF+IELF)
              
              ! Now we loop through the other DOF`s on the current element
              ! (the "X"`s).
              ! All these have common support with our current basis function
              ! and will therefore give an additive value to the global
              ! matrix.

              do JDOFE=1,indofCoarse
                
                ! Get the global DOF - our "X". This gives the column number
                ! in the matrix where an entry occurs in row IROW (the line of
                ! the current global DOF "O").
                JCOL = p_IdofsCoarse(JDOFE,IELC)

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

                        call storage_new ('mlop_create2LvlMatStruct9_conf', 'Ihicol', &
                                            p_Isize (iblocks), ST_INT, p_Ihcol(iblocks), &
                                            ST_NEWBLOCK_NOINIT)
                        call storage_getbase_int (p_Ihcol(iblocks),p_Icol)

                        ! The new index array must be filled with 0 - otherwise
                        ! the search routine below will not work!
                        call storage_new ('mlop_create2LvlMatStruct9_conf', 'p_Ihindx', &
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
                        ! it is in the last block... and so set the pointer
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
          
          end do ! IELF

          ! Add the number of refined elements to the counter
          NELF = NELF + NELREF
        
        end do ! IELC
        
        ! Add the number of processed coarse mesh elements
        nelementsDone = nelementsDone + nelementsToDo
      
      end do ! WHILE(nelementsDone .LT. NEL)
    
    end do ! IELDIST
    
    ! Release the fine mesh element list
    deallocate(p_IelementRef)

    ! Clean up the DOF`s arrays
    deallocate(p_IdofsFine)
    deallocate(p_IdofsCoarse)
    
    ! Ok, p_Icol is built. The hardest part is done!
    ! Now build KCOL by collecting the entries in the linear lists of
    ! each row.
    !
    ! At first, as we now NA, we can allocate the real KCOL now!
    call storage_new ('mlop_create2LvlMatStruct9_conf', 'KCOL', &
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

      ! Grab the largest column number. As the current line is sorted,
      ! we can find this using the end of the line.
      nmaxCol = max(nmaxCol,p_Kcol(p_Kld(IEQ+1)-1))

    end do ! IEQ
    
    ! HOORAY, THAT IS IT!
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

  subroutine mlop_build2LvlMass9_conf (rdiscretisationCoarse,&
                      rdiscretisationFine,bclear,rmatrixScalar,rperfconfig)
  
!<description>
  ! This routine calculates the entries of a 2-Level mass matrix.
  ! The matrix structure must have been initialised by the
  ! mlop_create2LvlMatrixStruct routine before calling this function.
!</description>

!<input>
  ! The underlying discretisation structure defined on the coarse mesh which
  ! is to be used to create the matrix.
  type(t_spatialDiscretisation), intent(in), target :: rdiscretisationCoarse
  
  ! The underlying discretisation structure defined on the fine mesh which
  ! is to be used to create the matrix.
  type(t_spatialDiscretisation), intent(in), target :: rdiscretisationFine
  
  ! Whether to clear the matrix before calculating the entries.
  ! If .FALSE., the new matrix entries are added to the existing entries.
  logical, intent(in) :: bclear
  
  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! The FE matrix. Calculated matrix entries are imposed to this matrix.
  type(t_matrixScalar), intent(inout) :: rmatrixScalar
!</inputoutput>

!</subroutine>


  ! local variables
  integer :: i,JDFG, ICUBP, NELC,NELF, IELDIST
  integer :: IELC,IELF, IDXC, NELREF, IDOFE, JDOFE
  integer :: JCOL0,JCOL
  real(DP) :: OM, DB
  
  ! Pointer to the performance configuration
  type(t_perfconfig), pointer :: p_rperfconfig

  ! Array to tell the element which derivatives to calculate
  logical, dimension(EL_MAXNDER) :: Bder
  
  ! number of cubature points on the reference element
  integer :: ncubpFine,ncubpCoarse,nmaxcubpFine,nmaxcubpCoarse
  
  ! Pointer to KLD, KCOL, DA
  integer, dimension(:), pointer :: p_KLD, p_KCOL
  real(DP), dimension(:), pointer :: p_DA
  
  ! An allocateable array accepting the DOF`s of a set of elements.
  integer, dimension(:,:), allocatable, target :: IdofsCoarse, IdofsFine
  
  ! Allocateable arrays for the values of the basis functions -
  ! for test and trial spaces.
  real(DP), dimension(:,:,:,:), allocatable, target :: DbasCoarse,DbasFine
  
  ! Number of entries in the matrix - for quicker access
  integer :: NA
  integer :: NEQ
  
  ! Type of transformation from the reference to the real element
  integer(I32) :: ctrafoCoarse, ctrafoFine
  
  ! Element evaluation tag; collects some information necessary for evaluating
  ! the elements.
  integer(I32) :: cevalTagCoarse, cevalTagFine
  
  ! Number of local degees of freedom for trial and test functions
  integer :: indofCoarse, indofFine,inmaxdofCoarse, inmaxdofFine
  
  ! The triangulation structure - to shorten some things...
  type(t_triangulation), pointer :: p_rtriaCoarse, p_rtriaFine
  
  ! A pointer to an element-number list
  integer, dimension(:), pointer :: p_IelementList, p_IelementRef
  
  ! Local matrices, used during the assembly.
  ! Values and positions of values in the global matrix.
  integer, dimension(:,:,:), allocatable :: Kentry
  real(DP), dimension(:,:,:), allocatable :: Dentry

  ! An array that takes coordinates of the cubature formula on the reference element
  real(DP), dimension(:,:), allocatable :: p_DcubPtsRefFine, p_DcubPtsRefCoarse
  real(DP), dimension(:), allocatable :: Domega
  
  ! Pointer to the jacobian determinants
  real(DP), dimension(:,:), pointer :: p_Ddetj

  ! Current element distribution
  type(t_elementDistribution), pointer :: p_relemDistCoarse, p_relemDistFine
  
  ! Number of elements that have already been processed and number of
  ! elements that are to be processed in the current run
  integer :: nelementsDone, nelementsToDo
  
  ! Number of elements that are to be processed at once
  integer :: nelementsCoarse, nelementsFine,nmaxelementsCoarse,nmaxelementsFine
    
  ! Two arrays for the refinement-patch arrays of the coarse triangulation
  integer, dimension(:), pointer :: p_IrefPatchIdx, p_IrefPatch
  
  ! Element evaluation structures for evaluation of finite elements
  type(t_evalElementSet) :: relementSetCoarse, relementSetFine
  
  integer :: nmaxDerCoarse, nmaxDerFine
  integer :: nmaxRefDimCoarse, nmaxRefDimFine

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => mlop_perfconfig
    end if

    ! We only need the function values as we want to assemble a mass matrix.
    Bder = .false.
    Bder(DER_FUNC) = .true.
    
    ! Get information about the matrix:
    NA = rmatrixScalar%NA
    NEQ = rmatrixScalar%NEQ
    
    ! We need KCOL/KLD of our matrix
    if ((rmatrixScalar%h_KCOL .eq. ST_NOHANDLE) .or. &
        (rmatrixScalar%h_KLD .eq. ST_NOHANDLE)) then
      call output_line ('No discretisation structure! Cannot assemble matrix!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'mlop_build2LvlMass9_conf')
      call sys_halt()
    end if
    
    call lsyssc_getbase_Kcol (rmatrixScalar,p_KCOL)
    call lsyssc_getbase_Kld (rmatrixScalar,p_KLD)
    
    ! Check if the matrix entries exist. If not, allocate the matrix.
    if (rmatrixScalar%h_DA .eq. ST_NOHANDLE) then

      ! Clear the entries in the matrix - we need to start with zero
      ! when assembling a new matrix!
      call storage_new ('mlop_build2LvlMass9_conf', 'DA', &
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
    
    ! Get a pointer to the triangulation - for easier access.
    p_rtriaCoarse => rdiscretisationCoarse%p_rtriangulation
    p_rtriaFine => rdiscretisationFine%p_rtriangulation

    ! Get the refinement patch arrays from the fine triangulation
    call storage_getbase_int(p_rtriaFine%h_IrefinementPatchIdx, p_IrefPatchIdx)
    call storage_getbase_int(p_rtriaFine%h_IrefinementPatch, p_IrefPatch)
    
    ! Let us loop over all element distributions and determine the
    ! maximum values.
    inmaxdofCoarse = 0
    inmaxdofFine = 0
    nmaxelementsCoarse = 0
    nmaxelementsFine = 0
    nmaxcubpCoarse = 0
    nmaxcubpFine = 0
    nmaxDerCoarse = 0
    nmaxDerFine = 0
    nmaxRefDimCoarse = 0
    nmaxRefDimFine = 0
    do i = 1, rdiscretisationCoarse%inumFESpaces
    
      ! Activate the current coarse mesh element distribution
      p_relemDistCoarse => rdiscretisationCoarse%RelementDistr(i)
      p_relemDistFine => rdiscretisationFine%RelementDistr(i)
      
      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      ctrafoCoarse = elem_igetTrafoType(p_relemDistCoarse%celement)
      ctrafoFine = elem_igetTrafoType(p_relemDistFine%celement)
    
      ! Get the number of local DOF`s for trial and test functions
      indofCoarse = elem_igetNDofLoc(p_relemDistCoarse%celement)
      indofFine = elem_igetNDofLoc(p_relemDistFine%celement)
      
      ! Get the number of elements
      nelementsCoarse = min(p_rperfconfig%NELEMSIM, p_relemDistCoarse%NEL)
      
      ! Get the number of cubature points on the fine mesh
      ncubpFine = cub_igetNumPts(p_relemDistFine%ccubTypeBilForm)
      
      ! Get the number of cubature points on the coarse mesh
      select case(cub_igetShape(p_relemDistFine%ccubTypeBilForm))
      case (BGEOM_SHAPE_LINE)
        ncubpCoarse = 2*ncubpFine
        nelementsFine = 2*nelementsCoarse
      case (BGEOM_SHAPE_TRIA,BGEOM_SHAPE_QUAD)
        ncubpCoarse = 4*ncubpFine
        nelementsFine = 4*nelementsCoarse
      case (BGEOM_SHAPE_HEXA)
        ncubpCoarse = 8*ncubpFine
        nelementsFine = 8*nelementsCoarse
      end select
      
      ! Determine maximum values
      inmaxdofCoarse = max(inmaxdofCoarse, indofCoarse)
      inmaxdofFine   = max(inmaxdofFine,   indofFine)
      nmaxelementsCoarse = max(nmaxelementsCoarse, nelementsCoarse)
      nmaxelementsFine   = max(nmaxelementsFine,   nelementsFine)
      nmaxcubpCoarse = max(nmaxcubpCoarse, ncubpCoarse)
      nmaxcubpFine   = max(nmaxcubpFine,   ncubpFine)
      nmaxDerCoarse = max(nmaxDerCoarse, &
          elem_getMaxDerivative(p_relemDistCoarse%celement))
      nmaxDerFine = max(nmaxDerFine, &
          elem_getMaxDerivative(p_relemDistFine%celement))
      nmaxRefDimCoarse = max(nmaxRefDimCoarse, &
          trafo_igetReferenceDimension(ctrafoCoarse))
      nmaxRefDimFine = max(nmaxRefDimFine, &
          trafo_igetReferenceDimension(ctrafoFine))

    end do
    
    ! Allocate an array saving a couple of DOF`s for trial and test functions
    allocate(IdofsCoarse(inmaxdofCoarse,nmaxelementsCoarse))
    allocate(IdofsFine(inmaxdofFine,nmaxelementsFine))

    ! And allocate the refined element list for the test functions
    allocate(p_IelementRef(nmaxelementsFine))

    ! Allocate some memory to hold the cubature points on the fine mesh
    allocate(p_DcubPtsRefFine(nmaxRefDimFine,nmaxcubpFine))
    allocate(p_DcubPtsRefCoarse(nmaxRefDimCoarse,nmaxcubpCoarse))
    allocate(Domega(nmaxcubpFine))

    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  ALLOCATE(DbasTest(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    !  ALLOCATE(DbasTrial(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    ! would lead to nonused memory blocks in these arrays during the assembly,
    ! which reduces the speed by 50%!
    allocate(DbasFine(inmaxdofFine,nmaxDerFine,nmaxcubpFine,nmaxelementsFine))
    allocate(DbasCoarse(inmaxdofCoarse,nmaxDerCoarse,nmaxcubpCoarse,nmaxelementsCoarse))
    
    ! Allocate an array saving the local matrices for all elements
    ! in an element set.
    ! We could also allocate EL_MAXNBAS*EL_MAXNBAS*NELEMSIM integers
    ! for this local matrix, but this would normally not fit to the cache
    ! anymore! indofTrial*indofTest*NELEMSIM is normally much smaller!
    allocate(Kentry(inmaxdofCoarse,inmaxdofFine,nmaxelementsFine))
    allocate(Dentry(inmaxdofFine,inmaxdofCoarse,nmaxelementsFine))
    
    ! Let us run through the element distributions
    do IELDIST = 1, rdiscretisationCoarse%inumFESpaces
    
      ! Activate the current element distributions
      p_relemDistCoarse => rdiscretisationCoarse%RelementDistr(IELDIST)
      p_relemDistFine => rdiscretisationFine%RelementDistr(IELDIST)

      if (p_relemDistCoarse%NEL .eq. 0) cycle
    
      ! p_IelementList must point to our set of elements in the discretisation
      ! with that the trial functions
      call storage_getbase_int (p_relemDistCoarse%h_IelementList, p_IelementList)

      ! Get the number of local DOF`s for trial and test functions
      indofCoarse = elem_igetNDofLoc(p_relemDistCoarse%celement)
      indofFine = elem_igetNDofLoc(p_relemDistFine%celement)
        
      ! Get the number of coarse mesh elements there.
      NELC = p_relemDistCoarse%NEL

      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      ctrafoCoarse = elem_igetTrafoType(p_relemDistCoarse%celement)
      ctrafoFine = elem_igetTrafoType(p_relemDistFine%celement)
      
      ! Initialise the cubature formula, get cubature weights and point
      ! coordinates on the reference element of the fine mesh
      ncubpFine = cub_igetNumPts(p_relemDistFine%ccubTypeBilForm)
      call cub_getCubature(p_relemDistFine%ccubTypeBilForm, p_DcubPtsRefFine, Domega)

      ! Calculate the number of coarse mesh elements we want to process
      ! in one run.
      nelementsCoarse = min(p_rperfconfig%NELEMSIM, NELC)
      
      ! Now we need to transform the points from the fine mesh into the coarse mesh
      ! Please note that the following trick does only work for 2-level ordered
      ! meshes!
      select case(cub_igetShape(p_relemDistFine%ccubTypeBilForm))
      case (BGEOM_SHAPE_LINE)
        nelementsFine = 2*nelementsCoarse
        call trafo_mapCubPtsRef2LvlEdge1D(ncubpFine,p_DcubPtsRefFine,p_DcubPtsRefCoarse)
      
      case (BGEOM_SHAPE_TRIA)
        nelementsFine = 4*nelementsCoarse
        call trafo_mapCubPtsRef2LvlTri2D(ncubpFine,p_DcubPtsRefFine,p_DcubPtsRefCoarse)

      case (BGEOM_SHAPE_QUAD)
        nelementsFine = 4*nelementsCoarse
        call trafo_mapCubPtsRef2LvlQuad2D(ncubpFine,p_DcubPtsRefFine,p_DcubPtsRefCoarse)

      case (BGEOM_SHAPE_HEXA)
        nelementsFine = 8*nelementsCoarse
        call trafo_mapCubPtsRef2LvlHexa3D(ncubpFine,p_DcubPtsRefFine,p_DcubPtsRefCoarse)
        
      end select
      
      ! Loop over the elements - blockwise.
      nelementsDone = 0
      do while(nelementsDone .lt. NELC)
      
        ! We always try to handle nelementsTrial elements simultaneously.
        ! Of course, we will not handle more elements than the coarse
        ! mesh discretisation has.
        nelementsToDo = min(NELC-nelementsDone, nelementsCoarse)
        
        ! Now comes the interesting part - we have to ensure that the DOF-mapping
        ! of the fine mesh discretisation fits into our DOF-array.
        ! If, for example, a coarse mesh quad was refined into more than 4 fine
        ! mesh quads, then it might happen that we cannot handle nelementsToDo
        ! coarse mesh elements at once, but we need to decrease nelementsToDo.
        NELF = 0
        do IELC = 1, nelementsToDo
        
          ! Get the index of the coarse mesh element
          IDXC = p_IelementList(nelementsDone+IELC)
          
          ! Get the number of fine mesh elements that have been refined from the
          ! currently processed coarse mesh element.
          NELREF = p_IrefPatchIdx(IDXC+1) - p_IrefPatchIdx(IDXC)
          
          ! Now if (NELF+NELREF) is greater than nelementsFine, then we need
          ! to decrease nelementsToDo and exit the loop...
          ! This case should never happen if the coarse mesh was refined using
          ! the 2-Level-Ordering algorithm, but might happen for more freaky
          ! refinement techniques...
          if((NELF+NELREF) .gt. nelementsFine) then
            nelementsToDo = IELC-1
            exit
          end if
          
          ! Copy the indices of the elements into the element list for the
          ! fine mesh discretisation
          do IELF = 1, NELREF
            p_IelementRef(NELF+IELF) = p_IrefPatch(p_IrefPatchIdx(IDXC)+IELF-1)
          end do
          
          ! Add the number of refined elements to the counter
          NELF = NELF + NELREF
        
        end do
        
        ! If nelementsToDo is 0, then we have a serious problem...
        if (nelementsToDo .le. 0) then
          call output_line ('INTERNAL ERROR: nelementsToDo = 0!', &
              OU_CLASS_ERROR,OU_MODE_STD,'mlop_build2LvlMass9_conf')
          call sys_halt()
        end if
        
        ! Call the DOF-mapping routine for the coarse and fine mesh
        call dof_locGlobMapping_mult(rdiscretisationCoarse, &
            p_IelementList(nelementsDone+1:nelementsDone+nelementsToDo), &
            IdofsCoarse)
        call dof_locGlobMapping_mult(rdiscretisationFine, p_IelementRef(1:NELF), &
            IdofsFine)
        
        ! ------------------- LOCAL MATRIX SETUP PHASE -----------------------
        NELF = 0
        do IELC = 1, nelementsToDo
        
          ! Get the index of the currently processed coarse mesh element
          IDXC = p_IelementList(nelementsDone+IELC)
          
          ! Get the number of fine mesh elements that are refined from the
          ! current coarse mesh element
          NELREF = p_IrefPatchIdx(IDXC+1) - p_IrefPatchIdx(IDXC)

          ! And loop through all elements of the current refinement patch
          do IELF = 1, NELREF
          
            ! Get the index of the currently processed fine mesh element
            !IDXF = p_IelementRef(NELF+IELF)
        
            ! For building the local matrices, we have first to
            ! loop through the test functions (the "O"`s), as these
            ! define the rows in the matrix.
            do IDOFE=1,indofFine
            
              ! Row IDOFE of the local matrix corresponds
              ! to row=global DOF KDFG(IDOFE) in the global matrix.
              ! This is one of the the "O"`s in the above picture.
              ! Get the starting position of the corresponding row
              ! to JCOL0:
              JCOL0 = p_KLD(IdofsFine(IDOFE,NELF+IELF))
              
              ! Now we loop through the other DOF`s on the current element
              ! (the "O"`s).
              ! All these have common support with our current basis function
              ! and will therefore give an additive value to the global
              ! matrix.
              do JDOFE=1,indofCoarse
                
                ! Get the global DOF of the "X" which interacts with
                ! our "O".
                JDFG = IdofsCoarse(JDOFE,IELC)
                
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
                
                Kentry(JDOFE,IDOFE,NELF+IELF)=JCOL
                
              end do ! IDOFE
              
            end do ! JDOFE
          
          end do ! IELF
          
          NELF = NELF + NELREF
          
        end do ! IELC
        
        ! -------------------- ELEMENT EVALUATION PHASE ----------------------
        
        ! Ok, we found the positions of the local matrix entries
        ! that we have to change.
        ! To calculate the matrix contributions, we have to evaluate
        ! the elements to give us the values of the basis functions
        ! in all the DOF`s in all the elements in our set.

        ! Get the element evaluation tag of all FE spaces. We need it to evaluate
        ! the elements later. All of them can be combined with OR, what will give
        ! a combined evaluation tag.
        cevalTagCoarse = elem_getEvaluationTag(p_relemDistCoarse%celement)
        cevalTagFine = elem_getEvaluationTag(p_relemDistFine%celement)
                        
        cevalTagCoarse = ior(cevalTagCoarse,EL_EVLTAG_REFPOINTS)
        cevalTagFine   = ior(cevalTagFine,EL_EVLTAG_REFPOINTS)

        ! Calculate all information that is necessary to evaluate the finite element
        ! on all cells of our subset. This includes the coordinates of the points
        ! on the cells.
        call elprep_prepareSetForEvaluation (relementSetCoarse,&
            cevalTagCoarse, p_rtriaCoarse, &
            p_IelementList(nelementsDone+1:nelementsDone+nelementsToDo), &
            ctrafoCoarse, p_DcubPtsRefCoarse(:,1:ncubpCoarse), rperfconfig=rperfconfig)

        call elprep_prepareSetForEvaluation (relementSetFine,&
            cevalTagFine, p_rtriaFine, p_IelementRef(1:NELF), &
            ctrafoFine, p_DcubPtsRefFine(:,1:ncubpFine), rperfconfig=rperfconfig)
        p_Ddetj => relementSetFine%p_Ddetj
        
        ! Calculate the values of the basis functions.
        call elem_generic_sim2 (p_relemDistCoarse%celement, &
            relementSetCoarse, Bder, DbasCoarse)
        call elem_generic_sim2 (p_relemDistFine%celement, &
            relementSetFine, Bder, DbasFine)
        
        ! --------------------- DOF COMBINATION PHASE ------------------------
        
        ! Values of all basis functions calculated. Now we can start
        ! to integrate!
        
        ! Clear the local matrix
        Dentry = 0.0_DP

        ! Loop over the elements in the current set.
        NELF = 0
        do IELC=1,nelementsToDo
        
          ! Get the index of the currently processed coarse mesh element
          IDXC = p_IelementList(nelementsDone+IELC)
          
          ! Get the number of fine mesh elements that are refined from the
          ! current coarse mesh element
          NELREF = p_IrefPatchIdx(IDXC+1) - p_IrefPatchIdx(IDXC)

          ! And loop through all elements of the current refinement patch
          do IELF = 1, NELREF
            
            ! Loop over all cubature points on the current element
            do ICUBP = 1, ncubpFine

              ! calculate the current weighting factor in the cubature formula
              ! in that cubature point.
              !
              ! Take the absolut value of the determinant of the mapping.
              ! In 2D, the determinant is always positive, whereas in 3D,
              ! the determinant might be negative -- that is normal!
              OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,NELF+IELF))

              ! Now loop through all possible combinations of DOF`s
              ! in the current cubature point. The outer loop
              ! loops through the "O"`s in the above picture,
              ! the test functions:
              do IDOFE=1,indofFine
              
                ! Get the value of the (test) basis function
                ! phi_i (our "O") in the cubature point:
                DB = DbasFine(IDOFE,DER_FUNC,ICUBP,NELF+IELF)*OM
                
                ! Perform an inner loop through the other DOF`s
                ! (the "X").
                do JDOFE=1,indofCoarse
                
                  Dentry(JDOFE,IDOFE,NELF+IELF) = Dentry(JDOFE,IDOFE,NELF+IELF) + &
                           DB*DbasCoarse(JDOFE,DER_FUNC,&
                           ICUBP + (IELF-1)*ncubpFine,IELC)
                
                end do ! JDOFE
              
              end do ! IDOFE

            end do ! ICUBP
            
          end do ! IELF
          
          NELF = NELF + NELREF

        end do ! IELC

        ! Incorporate the local matrix into the global one.
        ! Kentry gives the position of the additive contributions in Dentry.
        do IELF = 1, NELF
          do IDOFE = 1, indofFine
            do JDOFE = 1, indofCoarse
              p_DA(Kentry(JDOFE,IDOFE,IELF)) = &
                p_DA(Kentry(JDOFE,IDOFE,IELF)) + Dentry(JDOFE,IDOFE,IELF)
            end do
          end do
        end do

        ! Release the element sets here
        call elprep_releaseElementSet(relementSetFine)
        call elprep_releaseElementSet(relementSetCoarse)
        
        ! Increase the number of done elements
        nelementsDone = nelementsDone + nelementsToDo
        
      end do ! while(nelementsDone .lt. NELC)
    
    end do ! IELDIST
    
    ! Release memory
    deallocate(Domega)
    deallocate(p_DcubPtsRefCoarse)
    deallocate(p_DcubPtsRefFine)
    deallocate(IdofsCoarse)
    deallocate(IdofsFine)
    deallocate(DbasCoarse)
    deallocate(DbasFine)
    deallocate(Kentry)
    deallocate(Dentry)
    deallocate(p_IelementRef)

    ! That is it
  
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine mlop_build2LvlProl9_conf (rdiscretisationCoarse,&
                      rdiscretisationFine,bclear,rmatrixScalar,&
                      cavrgType,rperfconfig)
  
!<description>
  ! This routine calculates the entries of a 2-Level prolongation matrix.
  ! The matrix structure must have been initialised by the
  ! mlop_create2LvlMatrixStruct routine before calling this function.
!</description>

!<input>
  ! The underlying discretisation structure defined on the coarse mesh which
  ! is to be used to create the matrix.
  type(t_spatialDiscretisation), intent(in), target :: rdiscretisationCoarse
  
  ! The underlying discretisation structure defined on the fine mesh which
  ! is to be used to create the matrix.
  type(t_spatialDiscretisation), intent(in), target :: rdiscretisationFine
  
  ! Whether to clear the matrix before calculating the entries.
  ! If .FALSE., the new matrix entries are added to the existing entries.
  logical, intent(in) :: bclear
  
  ! Specifies which type of averaging is to be used.
  ! One of the MLOP_AVRG_XXXX constants defined above.
  integer, intent(in) :: cavrgType

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! The FE matrix. Calculated matrix entries are imposed to this matrix.
  type(t_matrixScalar), intent(inout) :: rmatrixScalar
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i,k,JDFG, ICUBP, NELC,NELF, IELDIST
  integer :: IELC,IELF, IDXC, NELREF, IDOFE, JDOFE
  integer :: JCOL0,JCOL
  real(DP) :: OM, DB

  ! Pointer to the performance configuration
  type(t_perfconfig), pointer :: p_rperfconfig
  
  ! Array to tell the element which derivatives to calculate
  logical, dimension(EL_MAXNDER) :: Bder
  
  ! number of cubature points on the reference element
  integer :: ncubpFine,ncubpCoarse,nmaxcubpFine,nmaxcubpCoarse
  
  ! Pointer to KLD, KCOL, DA
  integer, dimension(:), pointer :: p_KLD, p_KCOL
  real(DP), dimension(:), pointer :: p_DA
  
  ! An allocateable array accepting the DOF`s of a set of elements.
  integer, dimension(:,:), allocatable, target :: IdofsCoarse, IdofsFine
  
  ! Allocateable arrays for the values of the basis functions -
  ! for test and trial spaces.
  real(DP), dimension(:,:,:,:), allocatable, target :: DbasCoarse,DbasFine
  
  ! Number of entries in the matrix - for quicker access
  integer :: NA
  integer :: NEQ
  
  ! Type of transformation from the reference to the real element
  integer(I32) :: ctrafoCoarse, ctrafoFine
  
  ! Element evaluation tag; collects some information necessary for evaluating
  ! the elements.
  integer(I32) :: cevalTagCoarse, cevalTagFine
  
  ! Number of local degees of freedom for trial and test functions
  integer :: indofCoarse, indofFine, inmaxdofCoarse, inmaxdofFine
  
  ! The triangulation structure - to shorten some things...
  type(t_triangulation), pointer :: p_rtriaCoarse, p_rtriaFine
  
  ! A pointer to an element-number list
  integer, dimension(:), pointer :: p_IelementList, p_IelementRef
  
  ! Local matrices, used during the assembly.
  ! Values and positions of values in the global matrix.
  integer, dimension(:,:,:), allocatable :: Kentry
  real(DP), dimension(:,:,:), allocatable :: Dentry

  ! An array that takes coordinates of the cubature formula on the reference element
  real(DP), dimension(:,:), allocatable :: p_DcubPtsRefFine, p_DcubPtsRefCoarse
  real(DP), dimension(:), allocatable :: Domega
  
  ! Pointer to the jacobian determinants
  real(DP), dimension(:,:), pointer :: p_Ddetj

  ! Current element distribution
  type(t_elementDistribution), pointer :: p_relemDistCoarse, p_relemDistFine
  
  ! Number of elements that have already been processed and number of
  ! elements that are to be processed in the current run
  integer :: nelementsDone, nelementsToDo
  
  ! Number of elements that are to be processed at once
  integer :: nelementsCoarse, nelementsFine, nmaxelementsCoarse,nmaxelementsFine
    
  ! Two arrays for the refinement-patch arrays of the coarse triangulation
  integer, dimension(:), pointer :: p_IrefPatchIdx, p_IrefPatch
  
  ! Element evaluation structures for evaluation of finite elements
  type(t_evalElementSet) :: relementSetCoarse, relementSetFine
  
  ! Two arrays for the averaging weights
  real(DP), dimension(:), allocatable :: DglobWeights
  real(DP), dimension(:,:), allocatable :: DlocWeights
  
  ! One array for the local mass matrices
  real(DP), dimension(:,:,:), allocatable :: Dmass
  
  ! One pivot array for lapack
  integer, dimension(:), allocatable :: Ipivot
  
  ! Maximum derivatives
  integer :: nmaxDerFine, nmaxDerCoarse
  
  ! Maximum reference dimension
  integer :: nmaxRefDimFine, nmaxRefDimCoarse

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => mlop_perfconfig
    end if

    ! We only need the function values as we want to assemble a mass matrix.
    Bder = .false.
    Bder(DER_FUNC) = .true.
    
    ! Get information about the matrix:
    NA = rmatrixScalar%NA
    NEQ = rmatrixScalar%NEQ
    
    ! We need KCOL/KLD of our matrix
    if ((rmatrixScalar%h_KCOL .eq. ST_NOHANDLE) .or. &
        (rmatrixScalar%h_KLD .eq. ST_NOHANDLE)) then
      call output_line ('No matrix structure! Cannot assemble matrix!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'mlop_build2LvlProl9_conf')
      call sys_halt()
    end if
    
    call lsyssc_getbase_Kcol (rmatrixScalar,p_KCOL)
    call lsyssc_getbase_Kld (rmatrixScalar,p_KLD)
    
    ! Check if the matrix entries exist. If not, allocate the matrix.
    if (rmatrixScalar%h_DA .eq. ST_NOHANDLE) then

      ! Clear the entries in the matrix - we need to start with zero
      ! when assembling a new matrix!
      call storage_new ('mlop_build2LvlProl9_conf', 'DA', &
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
    
    ! Get a pointer to the triangulation - for easier access.
    p_rtriaCoarse => rdiscretisationCoarse%p_rtriangulation
    p_rtriaFine => rdiscretisationFine%p_rtriangulation

    ! Get the refinement patch arrays from the fine triangulation
    call storage_getbase_int(p_rtriaFine%h_IrefinementPatchIdx, p_IrefPatchIdx)
    call storage_getbase_int(p_rtriaFine%h_IrefinementPatch, p_IrefPatch)
    
    ! Let us loop over all element distributions and determine the
    ! maximum values.
    inmaxdofCoarse = 0
    inmaxdofFine = 0
    nmaxelementsCoarse = 0
    nmaxelementsFine = 0
    nmaxcubpCoarse = 0
    nmaxcubpFine = 0
    nmaxDerCoarse = 0
    nmaxDerFine = 0
    nmaxRefDimCoarse = 0
    nmaxRefDimFine = 0
    do i = 1, rdiscretisationCoarse%inumFESpaces
    
      ! Activate the current coarse mesh element distribution
      p_relemDistCoarse => rdiscretisationCoarse%RelementDistr(i)
      p_relemDistFine => rdiscretisationFine%RelementDistr(i)
      
      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      ctrafoCoarse = elem_igetTrafoType(p_relemDistCoarse%celement)
      ctrafoFine = elem_igetTrafoType(p_relemDistFine%celement)
    
      ! Get the number of local DOF`s for trial and test functions
      indofCoarse = elem_igetNDofLoc(p_relemDistCoarse%celement)
      indofFine = elem_igetNDofLoc(p_relemDistFine%celement)
      
      ! Get the number of elements
      nelementsCoarse = min(p_rperfconfig%NELEMSIM, p_relemDistCoarse%NEL)
      
      ! Get the number of cubature points on the fine mesh
      ncubpFine = cub_igetNumPts(p_relemDistFine%ccubTypeBilForm)
      
      ! Get the number of cubature points on the coarse mesh
      select case(cub_igetShape(p_relemDistFine%ccubTypeBilForm))
      case (BGEOM_SHAPE_LINE)
        ncubpCoarse = 2*ncubpFine
        nelementsFine = 2*nelementsCoarse
      case (BGEOM_SHAPE_TRIA,BGEOM_SHAPE_QUAD)
        ncubpCoarse = 4*ncubpFine
        nelementsFine = 4*nelementsCoarse
      case (BGEOM_SHAPE_HEXA)
        ncubpCoarse = 8*ncubpFine
        nelementsFine = 8*nelementsCoarse
      end select
      
      ! Determine maximum values
      inmaxdofCoarse = max(inmaxdofCoarse, indofCoarse)
      inmaxdofFine   = max(inmaxdofFine,   indofFine)
      nmaxelementsCoarse = max(nmaxelementsCoarse, nelementsCoarse)
      nmaxelementsFine   = max(nmaxelementsFine,   nelementsFine)
      nmaxcubpCoarse = max(nmaxcubpCoarse, ncubpCoarse)
      nmaxcubpFine   = max(nmaxcubpFine,   ncubpFine)
      nmaxDerCoarse = max(nmaxDerCoarse, &
          elem_getMaxDerivative(p_relemDistCoarse%celement))
      nmaxDerFine = max(nmaxDerFine, &
          elem_getMaxDerivative(p_relemDistFine%celement))
      nmaxRefDimCoarse = max(nmaxRefDimCoarse, &
          trafo_igetReferenceDimension(ctrafoCoarse))
      nmaxRefDimFine = max(nmaxRefDimFine, &
          trafo_igetReferenceDimension(ctrafoFine))

    end do

    ! Allocate an array saving a couple of DOF`s for trial and test functions
    allocate(IdofsCoarse(inmaxdofCoarse,nmaxelementsCoarse))
    allocate(IdofsFine(inmaxdofFine,nmaxelementsFine))

    ! And allocate the refined element list for the test functions
    allocate(p_IelementRef(nmaxelementsFine))

    ! Allocate some memory to hold the cubature points on the fine mesh
    allocate(p_DcubPtsRefFine(nmaxRefDimFine,nmaxcubpFine))
    allocate(p_DcubPtsRefCoarse(nmaxRefDimCoarse,nmaxcubpCoarse))
    allocate(Domega(nmaxcubpFine))

    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  ALLOCATE(DbasTest(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    !  ALLOCATE(DbasTrial(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    ! would lead to nonused memory blocks in these arrays during the assembly,
    ! which reduces the speed by 50%!
    allocate(DbasFine(inmaxdofFine,nmaxDerFine,nmaxcubpFine,nmaxelementsFine))
    allocate(DbasCoarse(inmaxdofCoarse,nmaxDerCoarse,nmaxcubpCoarse,nmaxelementsCoarse))

    ! Allocate an array saving the local matrices for all elements
    ! in an element set.
    ! We could also allocate EL_MAXNBAS*EL_MAXNBAS*NELEMSIM integers
    ! for this local matrix, but this would normally not fit to the cache
    ! anymore! indofTrial*indofTest*NELEMSIM is normally much smaller!
    allocate(Kentry(inmaxdofCoarse,inmaxdofFine,nmaxelementsFine))
    allocate(Dentry(inmaxdofFine,inmaxdofCoarse,nmaxelementsFine))
    allocate(Dmass(inmaxdofFine,inmaxdofFine,nmaxelementsFine))
    
    ! Allocate a pivot array for lapack
    allocate(Ipivot(inmaxdofFine))
    
    ! And allocate two arrays for the weights
    allocate(DglobWeights(NEQ))
    allocate(DlocWeights(inmaxdofFine,nmaxelementsFine))
    
    ! And format the global weights array
    call lalg_clearVectorDble(DglobWeights,NEQ)
    
    ! Format the local weights array to 1
    call lalg_setVectorDble2D(DlocWeights,1.0_DP)

    ! Let us run through the element distributions
    do IELDIST = 1, rdiscretisationCoarse%inumFESpaces

      ! Activate the current element distributions
      p_relemDistCoarse => rdiscretisationCoarse%RelementDistr(IELDIST)
      p_relemDistFine => rdiscretisationFine%RelementDistr(IELDIST)
    
      if (p_relemDistCoarse%NEL .eq. 0) cycle
    
      ! p_IelementList must point to our set of elements in the discretisation
      ! with that the trial functions
      call storage_getbase_int (p_relemDistCoarse%h_IelementList, p_IelementList)
      
      ! Get the number of coarse mesh elements there.
      NELC = p_relemDistCoarse%NEL

      ! Get the number of local DOF`s for trial and test functions
      indofCoarse = elem_igetNDofLoc(p_relemDistCoarse%celement)
      indofFine = elem_igetNDofLoc(p_relemDistFine%celement)
        
      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      ctrafoCoarse = elem_igetTrafoType(p_relemDistCoarse%celement)
      ctrafoFine = elem_igetTrafoType(p_relemDistFine%celement)
      
      ! Initialise the cubature formula, get cubature weights and point
      ! coordinates on the reference element of the fine mesh
      ncubpFine = cub_igetNumPts(p_relemDistFine%ccubTypeBilForm)
      call cub_getCubature(p_relemDistFine%ccubTypeBilForm, p_DcubPtsRefFine, Domega)
      
      ! Calculate the number of coarse mesh elements we want to process
      ! in one run.
      nelementsCoarse = min(p_rperfconfig%NELEMSIM,p_relemDistCoarse%NEL)

      ! Now we need to transform the points from the fine mesh into the coarse mesh
      ! Please note that the following trick does only work for 2-level ordered
      ! meshes!
      select case(cub_igetShape(p_relemDistFine%ccubTypeBilForm))
      case (BGEOM_SHAPE_LINE)
        nelementsFine = 2*nelementsCoarse
        call trafo_mapCubPtsRef2LvlEdge1D(ncubpFine,p_DcubPtsRefFine,p_DcubPtsRefCoarse)
      
      case (BGEOM_SHAPE_TRIA)
        nelementsFine = 4*nelementsCoarse
        call trafo_mapCubPtsRef2LvlTri2D(ncubpFine,p_DcubPtsRefFine,p_DcubPtsRefCoarse)

      case (BGEOM_SHAPE_QUAD)
        nelementsFine = 4*nelementsCoarse
        call trafo_mapCubPtsRef2LvlQuad2D(ncubpFine,p_DcubPtsRefFine,p_DcubPtsRefCoarse)

      case (BGEOM_SHAPE_HEXA)
        nelementsFine = 8*nelementsCoarse
        call trafo_mapCubPtsRef2LvlHexa3D(ncubpFine,p_DcubPtsRefFine,p_DcubPtsRefCoarse)
        
      end select
      
      ! Loop over the elements - blockwise.
      nelementsDone = 0
      do while(nelementsDone .lt. NELC)
      
        ! We always try to handle nelementsTrial elements simultaneously.
        ! Of course, we will not handle more elements than the coarse
        ! mesh discretisation has.
        nelementsToDo = min(NELC-nelementsDone, nelementsCoarse)
        
        ! Now comes the interesting part - we have to ensure that the DOF-mapping
        ! of the fine mesh discretisation fits into our DOF-array.
        ! If, for example, a coarse mesh quad was refined into more than 4 fine
        ! mesh quads, then it might happen that we cannot handle nelementsToDo
        ! coarse mesh elements at once, but we need to decrease nelementsToDo.
        NELF = 0
        do IELC = 1, nelementsToDo
        
          ! Get the index of the coarse mesh element
          IDXC = p_IelementList(nelementsDone+IELC)
          
          ! Get the number of fine mesh elements that have been refined from the
          ! currently processed coarse mesh element.
          NELREF = p_IrefPatchIdx(IDXC+1) - p_IrefPatchIdx(IDXC)
          
          ! Now if (NELF+NELREF) is greater than nelementsFine, then we need
          ! to decrease nelementsToDo and exit the loop...
          ! This case should never happen if the coarse mesh was refined using
          ! the 2-Level-Ordering algorithm, but might happen for more freaky
          ! refinement techniques...
          if((NELF+NELREF) .gt. nelementsFine) then
            nelementsToDo = IELC-1
            exit
          end if
          
          ! Copy the indices of the elements into the element list for the
          ! fine mesh discretisation
          do IELF = 1, NELREF
            p_IelementRef(NELF+IELF) = p_IrefPatch(p_IrefPatchIdx(IDXC)+IELF-1)
          end do
          
          ! Add the number of refined elements to the counter
          NELF = NELF + NELREF
        
        end do
        
        ! If nelementsToDo is 0, then we have a serious problem...
        if (nelementsToDo .le. 0) then
          call output_line ('INTERNAL ERROR: nelementsToDo = 0!', &
              OU_CLASS_ERROR,OU_MODE_STD,'mlop_build2LvlProl9_conf')
          call sys_halt()
        end if
        
        ! Call the DOF-mapping routine for the coarse and fine mesh
        call dof_locGlobMapping_mult(rdiscretisationCoarse, &
            p_IelementList(nelementsDone+1:nelementsDone+nelementsToDo), &
            IdofsCoarse)
        call dof_locGlobMapping_mult(rdiscretisationFine, p_IelementRef(1:NELF), &
            IdofsFine)
        
        ! ------------------- LOCAL MATRIX SETUP PHASE -----------------------
        NELF = 0
        do IELC = 1, nelementsToDo
        
          ! Get the index of the currently processed coarse mesh element
          IDXC = p_IelementList(nelementsDone+IELC)
          
          ! Get the number of fine mesh elements that are refined from the
          ! current coarse mesh element
          NELREF = p_IrefPatchIdx(IDXC+1) - p_IrefPatchIdx(IDXC)

          ! And loop through all elements of the current refinement patch
          do IELF = 1, NELREF
          
            ! Get the index of the currently processed fine mesh element
            !IDXF = p_IelementRef(NELF+IELF)
        
            ! For building the local matrices, we have first to
            ! loop through the test functions (the "O"`s), as these
            ! define the rows in the matrix.
            do IDOFE=1,indofFine
            
              ! Row IDOFE of the local matrix corresponds
              ! to row=global DOF KDFG(IDOFE) in the global matrix.
              ! This is one of the the "O"`s in the above picture.
              ! Get the starting position of the corresponding row
              ! to JCOL0:
              JCOL0 = p_KLD(IdofsFine(IDOFE,NELF+IELF))
              
              ! Now we loop through the other DOF`s on the current element
              ! (the "O"`s).
              ! All these have common support with our current basis function
              ! and will therefore give an additive value to the global
              ! matrix.
              do JDOFE=1,indofCoarse
                
                ! Get the global DOF of the "X" which interacts with
                ! our "O".
                JDFG = IdofsCoarse(JDOFE,IELC)
                
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
                
                Kentry(JDOFE,IDOFE,NELF+IELF)=JCOL
                
              end do ! IDOFE
              
            end do ! JDOFE
          
          end do ! IELF
          
          NELF = NELF + NELREF
          
        end do ! IELC
        
        ! -------------------- ELEMENT EVALUATION PHASE ----------------------
        
        ! Ok, we found the positions of the local matrix entries
        ! that we have to change.
        ! To calculate the matrix contributions, we have to evaluate
        ! the elements to give us the values of the basis functions
        ! in all the DOF`s in all the elements in our set.

        ! Get the element evaluation tag of all FE spaces. We need it to evaluate
        ! the elements later. All of them can be combined with OR, what will give
        ! a combined evaluation tag.
        cevalTagCoarse = elem_getEvaluationTag(p_relemDistCoarse%celement)
        cevalTagFine = elem_getEvaluationTag(p_relemDistFine%celement)
                        
        cevalTagCoarse = ior(cevalTagCoarse,EL_EVLTAG_REFPOINTS)
        cevalTagFine   = ior(cevalTagFine,EL_EVLTAG_REFPOINTS)

        ! Calculate all information that is necessary to evaluate the finite element
        ! on all cells of our subset. This includes the coordinates of the points
        ! on the cells.
        call elprep_prepareSetForEvaluation (relementSetCoarse,&
            cevalTagCoarse, p_rtriaCoarse, &
            p_IelementList(nelementsDone+1:nelementsDone+nelementsToDo), &
            ctrafoCoarse, p_DcubPtsRefCoarse(:,1:ncubpCoarse), rperfconfig=rperfconfig)

        call elprep_prepareSetForEvaluation (relementSetFine,&
            cevalTagFine, p_rtriaFine, p_IelementRef(1:NELF), &
            ctrafoFine, p_DcubPtsRefFine(:,1:ncubpFine), rperfconfig=rperfconfig)
        p_Ddetj => relementSetFine%p_Ddetj
        
        ! Calculate the values of the basis functions.
        call elem_generic_sim2 (p_relemDistCoarse%celement, &
            relementSetCoarse, Bder, DbasCoarse)
        call elem_generic_sim2 (p_relemDistFine%celement, &
            relementSetFine, Bder, DbasFine)
        
        ! --------------------- DOF COMBINATION PHASE ------------------------
        
        ! Values of all basis functions calculated. Now we can start
        ! to integrate!
        
        ! Clear the local matrix
        Dentry = 0.0_DP
        Dmass = 0.0_DP

        ! Loop over the elements in the current set.
        NELF = 0
        do IELC=1,nelementsToDo
        
          ! Get the index of the currently processed coarse mesh element
          IDXC = p_IelementList(nelementsDone+IELC)
          
          ! Get the number of fine mesh elements that are refined from the
          ! current coarse mesh element
          NELREF = p_IrefPatchIdx(IDXC+1) - p_IrefPatchIdx(IDXC)

          ! And loop through all elements of the current refinement patch
          do IELF = 1, NELREF
            
            ! Loop over all cubature points on the current element
            do ICUBP = 1, ncubpFine

              ! calculate the current weighting factor in the cubature formula
              ! in that cubature point.
              !
              ! Take the absolut value of the determinant of the mapping.
              ! In 2D, the determinant is always positive, whereas in 3D,
              ! the determinant might be negative -- that is normal!
              OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,NELF+IELF))

              ! Now loop through all possible combinations of DOF`s
              ! in the current cubature point. The outer loop
              ! loops through the "O"`s in the above picture,
              ! the test functions:
              do IDOFE=1,indofFine
              
                ! Get the value of the (test) basis function
                ! phi_i (our "O") in the cubature point:
                DB = DbasFine(IDOFE,DER_FUNC,ICUBP,NELF+IELF)*OM
                
                ! Perform an inner loop through the other DOF`s
                ! (the "X").
                do JDOFE=1,indofCoarse
                
                  Dentry(IDOFE,JDOFE,NELF+IELF) = Dentry(IDOFE,JDOFE,NELF+IELF) + &
                           DB*DbasCoarse(JDOFE,DER_FUNC,&
                           ICUBP + (IELF-1)*ncubpFine,IELC)
                
                end do ! JDOFE
                
                do JDOFE = 1, indofFine
                
                  Dmass(JDOFE,IDOFE,NELF+IELF) = Dmass(JDOFE,IDOFE,NELF+IELF) + &
                                  DB*DbasFine(JDOFE,DER_FUNC,ICUBP,NELF+IELF)
                  
                end do ! JDOFE
              
              end do ! IDOFE

            end do ! ICUBP
            
          end do ! IELF
          
          NELF = NELF + NELREF

        end do ! IELC
        
        ! If we use any mass-based averaging, we need to backup the mass here,
        ! as in the following step, the local mass matrices will be overwritten!
        select case(cavrgType)
        case (MLOP_AVRG_MASS)
          ! mass based averaging
          do IELF = 1, NELF
            do IDOFE=1, indofFine
              DlocWeights(IDOFE,IELF) = Dmass(IDOFE,IDOFE,IELF)
            end do ! IDOFE
          end do ! IELF
        
        case (MLOP_AVRG_INV_MASS)
          ! inverse mass based averaging
          do IELF = 1, NELF
            do IDOFE=1, indofFine
              DlocWeights(IDOFE,IELF) = 1.0_DP / Dmass(IDOFE,IDOFE,IELF)
            end do ! IDOFE
          end do ! IELF
        
        end select
        
        ! Now loop through all elements and perform the local L2-projection!
        do IELF = 1, NELF
          call DGESV(indofFine,indofCoarse,Dmass(:,:,IELF),inmaxdofFine,&
                     Ipivot,Dentry(:,:,IELF),inmaxdofFine,IDOFE)
        end do

        ! Incorporate the local matrix into the global one.
        ! Kentry gives the position of the additive contributions in Dentry.
        do IELF = 1, NELF
          do IDOFE=1,indofFine
            OM = DlocWeights(IDOFE,IELF)
            do JDOFE=1,indofCoarse
              p_DA(Kentry(JDOFE,IDOFE,IELF)) = &
                p_DA(Kentry(JDOFE,IDOFE,IELF)) + OM*Dentry(IDOFE,JDOFE,IELF)
            end do
            i = IdofsFine(IDOFE,IELF)
            DglobWeights(i) = DglobWeights(i) + OM
          end do
        end do

        ! Release the element sets here
        call elprep_releaseElementSet(relementSetFine)
        call elprep_releaseElementSet(relementSetCoarse)
        
        ! Increase the number of done elements
        nelementsDone = nelementsDone + nelementsToDo
        
      end do ! while(nelementsDone .lt. NELC)
    
    end do ! IELDIST
    
    ! Okay, there is one last thing left: Scale the matrix rows by the
    ! global averaging weights!
    do i = 1, NEQ
    
      ! Get the weight
      OM = DglobWeights(i)
      
      ! Skip it if it is zero
      if(OM .eq. 0.0_DP) continue
      
      ! Invert the weight
      OM = 1.0_DP / OM
      
      ! And scale the i-th matrix row by OM
      do k = p_Kld(i), p_Kld(i+1)-1
        p_DA(k) = p_DA(k) * OM
      end do
    
    end do ! i
    
    ! Release memory
    deallocate(DlocWeights)
    deallocate(DglobWeights)
    deallocate(Ipivot)
    deallocate(Dmass)
    deallocate(Domega)
    deallocate(p_DcubPtsRefCoarse)
    deallocate(p_DcubPtsRefFine)
    deallocate(IdofsCoarse)
    deallocate(IdofsFine)
    deallocate(DbasCoarse)
    deallocate(DbasFine)
    deallocate(Kentry)
    deallocate(Dentry)
    deallocate(p_IelementRef)

    ! That is it
  
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine mlop_build2LvlInterp9_conf (rdiscretisationCoarse,&
                      rdiscretisationFine,bclear,rmatrixScalar,&
                      cavrgType,rperfconfig)
  
!<description>
  ! This routine calculates the entries of a 2-Level interpolation matrix.
  ! The matrix structure must have been initialised by the
  ! mlop_create2LvlMatrixStruct routine before calling this function.
!</description>

!<input>
  ! The underlying discretisation structure defined on the coarse mesh which
  ! is to be used to create the matrix.
  type(t_spatialDiscretisation), intent(in), target :: rdiscretisationCoarse
  
  ! The underlying discretisation structure defined on the fine mesh which
  ! is to be used to create the matrix.
  type(t_spatialDiscretisation), intent(in), target :: rdiscretisationFine
  
  ! Whether to clear the matrix before calculating the entries.
  ! If .FALSE., the new matrix entries are added to the existing entries.
  logical, intent(in) :: bclear
  
  ! Specifies which type of averaging is to be used.
  ! One of the MLOP_AVRG_XXXX constants defined above.
  integer, intent(in) :: cavrgType

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! The FE matrix. Calculated matrix entries are imposed to this matrix.
  type(t_matrixScalar), intent(inout) :: rmatrixScalar
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i,k,JDFG, ICUBP, NELC,NELF, IELDIST
  integer :: IELC,IELF, IDXC, NELREF, IDOFE, JDOFE
  integer :: JCOL0,JCOL
  real(DP) :: OM, DB
  
  ! Pointer to the performance configuration
  type(t_perfconfig), pointer :: p_rperfconfig

  ! Array to tell the element which derivatives to calculate
  logical, dimension(EL_MAXNDER) :: Bder
  
  ! number of cubature points on the reference element
  integer :: ncubpFine,ncubpCoarse,nmaxcubpFine,nmaxcubpCoarse
  
  ! Pointer to KLD, KCOL, DA
  integer, dimension(:), pointer :: p_KLD, p_KCOL
  real(DP), dimension(:), pointer :: p_DA
  
  ! An allocateable array accepting the DOF`s of a set of elements.
  integer, dimension(:,:), allocatable, target :: IdofsCoarse, IdofsFine
  
  ! Allocateable arrays for the values of the basis functions -
  ! for test and trial spaces.
  real(DP), dimension(:,:,:,:), allocatable, target :: DbasCoarse,DbasFine
  
  ! Number of entries in the matrix - for quicker access
  integer :: NA, NEQ, NCOLS
  
  ! Type of transformation from the reference to the real element
  integer(I32) :: ctrafoCoarse, ctrafoFine
  
  ! Element evaluation tag; collects some information necessary for evaluating
  ! the elements.
  integer(I32) :: cevalTagCoarse, cevalTagFine
  
  ! Number of local degees of freedom for trial and test functions
  integer :: indofCoarse, indofFine, inmaxdofCoarse, inmaxdofFine
  
  ! The triangulation structure - to shorten some things...
  type(t_triangulation), pointer :: p_rtriaCoarse, p_rtriaFine
  
  ! A pointer to an element-number list
  integer, dimension(:), pointer :: p_IelementList, p_IelementRef
  
  ! Local matrices, used during the assembly.
  ! Values and positions of values in the global matrix.
  integer, dimension(:,:,:), allocatable :: Kentry
  real(DP), dimension(:,:,:), allocatable :: Dentry

  ! An array that takes coordinates of the cubature formula on the reference element
  real(DP), dimension(:,:), allocatable :: p_DcubPtsRefFine, p_DcubPtsRefCoarse
  real(DP), dimension(:), allocatable :: Domega
  
  ! Pointer to the jacobian determinants
  real(DP), dimension(:,:), pointer :: p_Ddetj

  ! Current element distribution
  type(t_elementDistribution), pointer :: p_relemDistCoarse, p_relemDistFine
  
  ! Number of elements that have already been processed and number of
  ! elements that are to be processed in the current run
  integer :: nelementsDone, nelementsToDo
  
  ! Number of elements that are to be processed at once
  integer :: nelementsCoarse, nelementsFine, nmaxelementsCoarse,nmaxelementsFine
    
  ! Two arrays for the refinement-patch arrays of the coarse triangulation
  integer, dimension(:), pointer :: p_IrefPatchIdx, p_IrefPatch
  
  ! Element evaluation structures for evaluation of finite elements
  type(t_evalElementSet) :: relementSetCoarse, relementSetFine
  
  ! Two arrays for the averaging weights
  real(DP), dimension(:), allocatable :: DglobWeights
  real(DP), dimension(:,:), allocatable :: DlocWeights
  
  ! One array for the local mass matrices
  real(DP), dimension(:,:,:), allocatable :: Dmass
  
  ! One pivot array for lapack
  integer, dimension(:), allocatable :: Ipivot
  
  ! Maximum derivatives
  integer :: nmaxDerFine, nmaxDerCoarse
  
  ! Maximum reference dimension
  integer :: nmaxRefDimFine, nmaxRefDimCoarse

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => mlop_perfconfig
    end if

    ! We only need the function values as we want to assemble a mass matrix.
    Bder = .false.
    Bder(DER_FUNC) = .true.
    
    ! Get information about the matrix:
    NA = rmatrixScalar%NA
    NEQ = rmatrixScalar%NEQ
    NCOLS = rmatrixScalar%NCOLS
    
    ! We need KCOL/KLD of our matrix
    if ((rmatrixScalar%h_KCOL .eq. ST_NOHANDLE) .or. &
        (rmatrixScalar%h_KLD .eq. ST_NOHANDLE)) then
      call output_line ('No matrix structure! Cannot assemble matrix!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'mlop_build2LvlInterp9_conf')
      call sys_halt()
    end if
    
    call lsyssc_getbase_Kcol (rmatrixScalar,p_KCOL)
    call lsyssc_getbase_Kld (rmatrixScalar,p_KLD)
    
    ! Check if the matrix entries exist. If not, allocate the matrix.
    if (rmatrixScalar%h_DA .eq. ST_NOHANDLE) then

      ! Clear the entries in the matrix - we need to start with zero
      ! when assembling a new matrix!
      call storage_new ('mlop_build2LvlInterp9_conf', 'DA', &
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
    
    ! Get a pointer to the triangulation - for easier access.
    p_rtriaCoarse => rdiscretisationCoarse%p_rtriangulation
    p_rtriaFine => rdiscretisationFine%p_rtriangulation

    ! Get the refinement patch arrays from the fine triangulation
    call storage_getbase_int(p_rtriaFine%h_IrefinementPatchIdx, p_IrefPatchIdx)
    call storage_getbase_int(p_rtriaFine%h_IrefinementPatch, p_IrefPatch)
    
    ! Let us loop over all element distributions and determine the
    ! maximum values.
    inmaxdofCoarse = 0
    inmaxdofFine = 0
    nmaxelementsCoarse = 0
    nmaxelementsFine = 0
    nmaxcubpCoarse = 0
    nmaxcubpFine = 0
    nmaxDerCoarse = 0
    nmaxDerFine = 0
    nmaxRefDimCoarse = 0
    nmaxRefDimFine = 0
    do i = 1, rdiscretisationCoarse%inumFESpaces
    
      ! Activate the current coarse mesh element distribution
      p_relemDistCoarse => rdiscretisationCoarse%RelementDistr(i)
      p_relemDistFine => rdiscretisationFine%RelementDistr(i)
      
      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      ctrafoCoarse = elem_igetTrafoType(p_relemDistCoarse%celement)
      ctrafoFine = elem_igetTrafoType(p_relemDistFine%celement)
    
      ! Get the number of local DOF`s for trial and test functions
      indofCoarse = elem_igetNDofLoc(p_relemDistCoarse%celement)
      indofFine = elem_igetNDofLoc(p_relemDistFine%celement)
      
      ! Get the number of elements
      nelementsCoarse = min(p_rperfconfig%NELEMSIM, p_relemDistCoarse%NEL)
      
      ! Get the number of cubature points on the fine mesh
      ncubpFine = cub_igetNumPts(p_relemDistFine%ccubTypeBilForm)
      
      ! Get the number of cubature points on the coarse mesh
      select case(cub_igetShape(p_relemDistFine%ccubTypeBilForm))
      case (BGEOM_SHAPE_LINE)
        ncubpCoarse = 2*ncubpFine
        nelementsFine = 2*nelementsCoarse
      case (BGEOM_SHAPE_TRIA,BGEOM_SHAPE_QUAD)
        ncubpCoarse = 4*ncubpFine
        nelementsFine = 4*nelementsCoarse
      case (BGEOM_SHAPE_HEXA)
        ncubpCoarse = 8*ncubpFine
        nelementsFine = 8*nelementsCoarse
      end select
      
      ! Determine maximum values
      inmaxdofCoarse = max(inmaxdofCoarse, indofCoarse)
      inmaxdofFine   = max(inmaxdofFine,   indofFine)
      nmaxelementsCoarse = max(nmaxelementsCoarse, nelementsCoarse)
      nmaxelementsFine   = max(nmaxelementsFine,   nelementsFine)
      nmaxcubpCoarse = max(nmaxcubpCoarse, ncubpCoarse)
      nmaxcubpFine   = max(nmaxcubpFine,   ncubpFine)
      nmaxDerCoarse = max(nmaxDerCoarse, &
          elem_getMaxDerivative(p_relemDistCoarse%celement))
      nmaxDerFine = max(nmaxDerFine, &
          elem_getMaxDerivative(p_relemDistFine%celement))
      nmaxRefDimCoarse = max(nmaxRefDimCoarse, &
          trafo_igetReferenceDimension(ctrafoCoarse))
      nmaxRefDimFine = max(nmaxRefDimFine, &
          trafo_igetReferenceDimension(ctrafoFine))

    end do

    ! Allocate an array saving a couple of DOF`s for trial and test functions
    allocate(IdofsCoarse(inmaxdofCoarse,nmaxelementsCoarse))
    allocate(IdofsFine(inmaxdofFine,nmaxelementsFine))

    ! And allocate the refined element list for the test functions
    allocate(p_IelementRef(nmaxelementsFine))

    ! Allocate some memory to hold the cubature points on the fine mesh
    allocate(p_DcubPtsRefFine(nmaxRefDimFine,nmaxcubpFine))
    allocate(p_DcubPtsRefCoarse(nmaxRefDimCoarse,nmaxcubpCoarse))
    allocate(Domega(nmaxcubpFine))

    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  ALLOCATE(DbasTest(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    !  ALLOCATE(DbasTrial(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    ! would lead to nonused memory blocks in these arrays during the assembly,
    ! which reduces the speed by 50%!
    allocate(DbasFine(inmaxdofFine,nmaxDerFine,nmaxcubpFine,nmaxelementsFine))
    allocate(DbasCoarse(inmaxdofCoarse,nmaxDerCoarse,nmaxcubpCoarse,nmaxelementsCoarse))

    ! Allocate an array saving the local matrices for all elements
    ! in an element set.
    allocate(Kentry(inmaxdofCoarse,inmaxdofFine,nmaxelementsFine))
    allocate(Dentry(inmaxdofCoarse,inmaxdofFine,nmaxelementsFine))
    allocate(Dmass(inmaxdofCoarse,inmaxdofCoarse,nmaxelementsCoarse))
    
    ! Allocate a pivot array for lapack
    allocate(Ipivot(inmaxdofCoarse))
    
    ! And allocate two arrays for the weights
    allocate(DglobWeights(NCOLS))
    allocate(DlocWeights(inmaxdofCoarse,nmaxelementsCoarse))
    
    ! And format the global weights array
    call lalg_clearVectorDble(DglobWeights,NCOLS)
    
    ! Format the local weights array to 1
    call lalg_setVectorDble2D(DlocWeights,1.0_DP)

    ! Let us run through the element distributions
    do IELDIST = 1, rdiscretisationCoarse%inumFESpaces

      ! Activate the current element distributions
      p_relemDistCoarse => rdiscretisationCoarse%RelementDistr(IELDIST)
      p_relemDistFine => rdiscretisationFine%RelementDistr(IELDIST)
    
      if (p_relemDistCoarse%NEL .eq. 0) cycle
    
      ! p_IelementList must point to our set of elements in the discretisation
      ! with that the trial functions
      call storage_getbase_int (p_relemDistCoarse%h_IelementList, p_IelementList)
      
      ! Get the number of coarse mesh elements there.
      NELC = p_relemDistCoarse%NEL

      ! Get the number of local DOF`s for trial and test functions
      indofCoarse = elem_igetNDofLoc(p_relemDistCoarse%celement)
      indofFine = elem_igetNDofLoc(p_relemDistFine%celement)
        
      ! Get from the trial element space the type of coordinate system
      ! that is used there:
      ctrafoCoarse = elem_igetTrafoType(p_relemDistCoarse%celement)
      ctrafoFine = elem_igetTrafoType(p_relemDistFine%celement)
      
      ! Initialise the cubature formula, get cubature weights and point
      ! coordinates on the reference element of the fine mesh
      ncubpFine = cub_igetNumPts(p_relemDistFine%ccubTypeBilForm)
      call cub_getCubature(p_relemDistFine%ccubTypeBilForm, p_DcubPtsRefFine, Domega)
      
      ! Calculate the number of coarse mesh elements we want to process
      ! in one run.
      nelementsCoarse = min(p_rperfconfig%NELEMSIM,p_relemDistCoarse%NEL)

      ! Now we need to transform the points from the fine mesh into the coarse mesh
      ! Please note that the following trick does only work for 2-level ordered
      ! meshes!
      select case(cub_igetShape(p_relemDistFine%ccubTypeBilForm))
      case (BGEOM_SHAPE_LINE)
        nelementsFine = 2*nelementsCoarse
        call trafo_mapCubPtsRef2LvlEdge1D(ncubpFine,p_DcubPtsRefFine,p_DcubPtsRefCoarse)
      
      case (BGEOM_SHAPE_TRIA)
        nelementsFine = 4*nelementsCoarse
        call trafo_mapCubPtsRef2LvlTri2D(ncubpFine,p_DcubPtsRefFine,p_DcubPtsRefCoarse)

      case (BGEOM_SHAPE_QUAD)
        nelementsFine = 4*nelementsCoarse
        call trafo_mapCubPtsRef2LvlQuad2D(ncubpFine,p_DcubPtsRefFine,p_DcubPtsRefCoarse)

      case (BGEOM_SHAPE_HEXA)
        nelementsFine = 8*nelementsCoarse
        call trafo_mapCubPtsRef2LvlHexa3D(ncubpFine,p_DcubPtsRefFine,p_DcubPtsRefCoarse)
        
      end select
      
      ! Loop over the elements - blockwise.
      nelementsDone = 0
      do while(nelementsDone .lt. NELC)
      
        ! We always try to handle nelementsTrial elements simultaneously.
        ! Of course, we will not handle more elements than the coarse
        ! mesh discretisation has.
        nelementsToDo = min(NELC-nelementsDone, nelementsCoarse)
        
        ! Now comes the interesting part - we have to ensure that the DOF-mapping
        ! of the fine mesh discretisation fits into our DOF-array.
        ! If, for example, a coarse mesh quad was refined into more than 4 fine
        ! mesh quads, then it might happen that we cannot handle nelementsToDo
        ! coarse mesh elements at once, but we need to decrease nelementsToDo.
        NELF = 0
        do IELC = 1, nelementsToDo
        
          ! Get the index of the coarse mesh element
          IDXC = p_IelementList(nelementsDone+IELC)
          
          ! Get the number of fine mesh elements that have been refined from the
          ! currently processed coarse mesh element.
          NELREF = p_IrefPatchIdx(IDXC+1) - p_IrefPatchIdx(IDXC)
          
          ! Now if (NELF+NELREF) is greater than nelementsFine, then we need
          ! to decrease nelementsToDo and exit the loop...
          ! This case should never happen if the coarse mesh was refined using
          ! the 2-Level-Ordering algorithm, but might happen for more freaky
          ! refinement techniques...
          if((NELF+NELREF) .gt. nelementsFine) then
            nelementsToDo = IELC-1
            exit
          end if
          
          ! Copy the indices of the elements into the element list for the
          ! fine mesh discretisation
          do IELF = 1, NELREF
            p_IelementRef(NELF+IELF) = p_IrefPatch(p_IrefPatchIdx(IDXC)+IELF-1)
          end do
          
          ! Add the number of refined elements to the counter
          NELF = NELF + NELREF
        
        end do
        
        ! If nelementsToDo is 0, then we have a serious problem...
        if (nelementsToDo .le. 0) then
          call output_line ('INTERNAL ERROR: nelementsToDo = 0!', &
              OU_CLASS_ERROR,OU_MODE_STD,'mlop_build2LvlInterp9_conf')
          call sys_halt()
        end if
        
        ! Call the DOF-mapping routine for the coarse and fine mesh
        call dof_locGlobMapping_mult(rdiscretisationCoarse, &
            p_IelementList(nelementsDone+1:nelementsDone+nelementsToDo), &
            IdofsCoarse)
        call dof_locGlobMapping_mult(rdiscretisationFine, p_IelementRef(1:NELF), &
            IdofsFine)
        
        ! ------------------- LOCAL MATRIX SETUP PHASE -----------------------
        NELF = 0
        do IELC = 1, nelementsToDo
        
          ! Get the index of the currently processed coarse mesh element
          IDXC = p_IelementList(nelementsDone+IELC)
          
          ! Get the number of fine mesh elements that are refined from the
          ! current coarse mesh element
          NELREF = p_IrefPatchIdx(IDXC+1) - p_IrefPatchIdx(IDXC)

          ! And loop through all elements of the current refinement patch
          do IELF = 1, NELREF
          
            ! For building the local matrices, we have first to
            ! loop through the test functions (the "O"`s), as these
            ! define the rows in the matrix.
            do IDOFE=1,indofFine
            
              ! Row IDOFE of the local matrix corresponds
              ! to row=global DOF KDFG(IDOFE) in the global matrix.
              ! This is one of the the "O"`s in the above picture.
              ! Get the starting position of the corresponding row
              ! to JCOL0:
              JCOL0 = p_KLD(IdofsFine(IDOFE,NELF+IELF))
              
              ! Now we loop through the other DOF`s on the current element
              ! (the "O"`s).
              ! All these have common support with our current basis function
              ! and will therefore give an additive value to the global
              ! matrix.
              do JDOFE=1,indofCoarse
                
                ! Get the global DOF of the "X" which interacts with
                ! our "O".
                JDFG = IdofsCoarse(JDOFE,IELC)
                
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
                
                Kentry(JDOFE,IDOFE,NELF+IELF)=JCOL
                
              end do ! IDOFE
              
            end do ! JDOFE
          
          end do ! IELF
          
          NELF = NELF + NELREF
          
        end do ! IELC
        
        ! -------------------- ELEMENT EVALUATION PHASE ----------------------
        
        ! Ok, we found the positions of the local matrix entries
        ! that we have to change.
        ! To calculate the matrix contributions, we have to evaluate
        ! the elements to give us the values of the basis functions
        ! in all the DOF`s in all the elements in our set.

        ! Get the element evaluation tag of all FE spaces. We need it to evaluate
        ! the elements later. All of them can be combined with OR, what will give
        ! a combined evaluation tag.
        cevalTagCoarse = elem_getEvaluationTag(p_relemDistCoarse%celement)
        cevalTagFine = elem_getEvaluationTag(p_relemDistFine%celement)
                        
        cevalTagCoarse = ior(cevalTagCoarse,EL_EVLTAG_REFPOINTS)
        cevalTagFine   = ior(cevalTagFine,EL_EVLTAG_REFPOINTS)

        ! Calculate all information that is necessary to evaluate the finite element
        ! on all cells of our subset. This includes the coordinates of the points
        ! on the cells.
        call elprep_prepareSetForEvaluation (relementSetCoarse,&
            cevalTagCoarse, p_rtriaCoarse, &
            p_IelementList(nelementsDone+1:nelementsDone+nelementsToDo), &
            ctrafoCoarse, p_DcubPtsRefCoarse(:,1:ncubpCoarse), rperfconfig=rperfconfig)

        call elprep_prepareSetForEvaluation (relementSetFine,&
            cevalTagFine, p_rtriaFine, p_IelementRef(1:NELF), &
            ctrafoFine, p_DcubPtsRefFine(:,1:ncubpFine), rperfconfig=rperfconfig)
        p_Ddetj => relementSetFine%p_Ddetj
        
        ! Calculate the values of the basis functions.
        call elem_generic_sim2 (p_relemDistCoarse%celement, &
            relementSetCoarse, Bder, DbasCoarse)
        call elem_generic_sim2 (p_relemDistFine%celement, &
            relementSetFine, Bder, DbasFine)
        
        ! --------------------- DOF COMBINATION PHASE ------------------------
        
        ! Values of all basis functions calculated. Now we can start
        ! to integrate!
        
        ! Clear the local matrix
        Dentry = 0.0_DP
        Dmass = 0.0_DP

        ! Loop over the elements in the current set.
        NELF = 0
        do IELC=1,nelementsToDo
        
          ! Get the index of the currently processed coarse mesh element
          IDXC = p_IelementList(nelementsDone+IELC)
          
          ! Get the number of fine mesh elements that are refined from the
          ! current coarse mesh element
          NELREF = p_IrefPatchIdx(IDXC+1) - p_IrefPatchIdx(IDXC)

          ! And loop through all elements of the current refinement patch
          do IELF = 1, NELREF
            
            ! Loop over all cubature points on the current element
            do ICUBP = 1, ncubpFine

              ! calculate the current weighting factor in the cubature formula
              ! in that cubature point.
              !
              ! Take the absolut value of the determinant of the mapping.
              ! In 2D, the determinant is always positive, whereas in 3D,
              ! the determinant might be negative -- that is normal!
              OM = Domega(ICUBP)*abs(p_Ddetj(ICUBP,NELF+IELF))

              ! Build the 2-level mass matrix.
              do IDOFE=1,indofFine
              
                ! Get the value of the (test) basis function
                ! phi_i (our "O") in the cubature point:
                DB = DbasFine(IDOFE,DER_FUNC,ICUBP,NELF+IELF)*OM
                
                ! Loop over the coarse mesh DOFs.
                do JDOFE = 1, indofCoarse
                
                  Dentry(JDOFE,IDOFE,NELF+IELF) = Dentry(JDOFE,IDOFE,NELF+IELF) + &
                           DB*DbasCoarse(JDOFE,DER_FUNC,ICUBP + (IELF-1)*ncubpFine,IELC)
                
                end do ! JDOFE
              
              end do ! IDOFE
              
              ! Build the coarse mesh mass matrix.
              do IDOFE = 1, indofCoarse
              
                DB = DbasCoarse(IDOFE,DER_FUNC,ICUBP + (IELF-1)*ncubpFine,IELC)*OM
                
                do JDOFE = 1, indofCoarse
                
                  Dmass(JDOFE,IDOFE,IELC) = Dmass(JDOFE,IDOFE,IELC) + &
                            DB*DbasCoarse(JDOFE,DER_FUNC,ICUBP + (IELF-1)*ncubpFine,IELC)
                  
                end do ! JDOFE
              
              end do ! IDOFE

            end do ! ICUBP
            
          end do ! IELF
          
          NELF = NELF + NELREF

        end do ! IELC
        
        ! If we use any mass-based averaging, we need to backup the mass here,
        ! as in the following step, the local mass matrices will be overwritten!
        select case(cavrgType)
        case (MLOP_AVRG_MASS)
          ! mass based averaging
          do IELC = 1, nelementsToDo
            do IDOFE = 1, indofCoarse
              DlocWeights(IDOFE,IELC) = Dmass(IDOFE,IDOFE,IELC)
            end do ! IDOFE
          end do ! IELF
        
        case (MLOP_AVRG_INV_MASS)
          ! inverse mass based averaging
          do IELC = 1, nelementsToDo
            do IDOFE = 1, indofCoarse
              DlocWeights(IDOFE,IELC) = 1.0_DP / Dmass(IDOFE,IDOFE,IELC)
            end do ! IDOFE
          end do ! IELF
        
        end select
        
        ! Loop over all coarse mesh elements
        NELF = 0
        do IELC = 1, nelementsToDo
        
          ! Factorise the coarse mesh mass matrix
          call DGETRF(indofCoarse,indofCoarse,Dmass(:,:,IELC),&
                      inmaxdofCoarse,Ipivot,IDOFE)
          
          ! Get the index of the currently processed coarse mesh element
          IDXC = p_IelementList(nelementsDone+IELC)
          
          ! Get the number of fine mesh elements
          NELREF = p_IrefPatchIdx(IDXC+1) - p_IrefPatchIdx(IDXC)
          
          ! Loop over all 2-level mass matrices for that coarse mesh element
          do IELF = 1, NELREF
          
            ! Solve the system
            call DGETRS('N',indofCoarse,indofFine,Dmass(:,:,IELC),&
                        inmaxdofCoarse,Ipivot,Dentry(:,:,NELF+IELF),&
                        inmaxdofFine,IDOFE)
            
            ! Scale the entries by the weights
            do IDOFE = 1, indofCoarse
            
              ! Get the weight
              OM = DlocWeights(IDOFE,IELC)
              
              do JDOFE = 1, indofFine
                Dentry(IDOFE,JDOFE,NELF+IELF) = OM*Dentry(IDOFE,JDOFE,NELF+IELF)
              end do ! JDOFE
            
            end do ! IDOFE
          
          end do ! IELF
          
          ! Update counter
          NELF = NELF + NELREF
          
          ! Update the global weights
          do IDOFE = 1, indofCoarse
            i = IdofsCoarse(IDOFE,IELC)
            DglobWeights(i) = DglobWeights(i) + DlocWeights(IDOFE,IELC)
          end do ! IDOFE
        
        end do ! IELC

        ! Incorporate the local matrix into the global one.
        ! Kentry gives the position of the additive contributions in Dentry.
        do IELF = 1, NELF
          do IDOFE=1,indofFine
            do JDOFE=1,indofCoarse
              p_DA(Kentry(JDOFE,IDOFE,IELF)) = &
                p_DA(Kentry(JDOFE,IDOFE,IELF)) + Dentry(JDOFE,IDOFE,IELF)
            end do
          end do
        end do

        ! Release the element sets here
        call elprep_releaseElementSet(relementSetFine)
        call elprep_releaseElementSet(relementSetCoarse)
        
        ! Increase the number of done elements
        nelementsDone = nelementsDone + nelementsToDo
        
      end do ! while(nelementsDone .lt. NELC)
    
    end do ! IELDIST
    
    ! Okay, there is one last thing left: Scale the matrix columns by the
    ! global averaging weights!
    ! So first invert the global weights
    do i = 1, NCOLS
    
      ! Get the weight
      OM = DglobWeights(i)
      
      ! Skip it if it is zero
      if(OM .eq. 0.0_DP) continue
      
      ! Invert the weight
      DglobWeights(i) = 1.0_DP / OM
    
    end do ! i
    
    ! Low loop over the matrix` rows
    do i = 1, NEQ
    
      ! Loop over the non-zeroes
      do k = p_Kld(i), p_Kld(i+1)-1
        p_DA(k) = p_DA(k) * DglobWeights(p_Kcol(k))
      end do ! j
    
    end do ! i
    
    ! Release memory
    deallocate(DlocWeights)
    deallocate(DglobWeights)
    deallocate(Ipivot)
    deallocate(Dmass)
    deallocate(Domega)
    deallocate(p_DcubPtsRefCoarse)
    deallocate(p_DcubPtsRefFine)
    deallocate(IdofsCoarse)
    deallocate(IdofsFine)
    deallocate(DbasCoarse)
    deallocate(DbasFine)
    deallocate(Kentry)
    deallocate(Dentry)
    deallocate(p_IelementRef)

    ! That is it
  
  end subroutine

end module
