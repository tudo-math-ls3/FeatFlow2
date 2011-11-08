!##############################################################################
!# ****************************************************************************
!# <name> cfdflow </name>
!# ****************************************************************************
!#
!# <purpose>
!# </purpose>
!##############################################################################

module mlproj

  use fsystem
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use multilevelprojection
  
  implicit none
  
!<types>

!<typeblock>

  ! A type encapsuling a pointer to a discretisation structure
  type t_mldiscrPointer
  
    ! Pointer to a discretisation
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    
  end type

!</typeblock>

!<typeblock>

  ! This structure defines a hierarchy of projection structures
  ! for a hierarchy of refinement levels. It can be used to project
  ! FEM functions given as block vectors from one level to another.
  type t_interlevelProjectionHier
  
    private
  
    ! Minimum level in the hierarchy.
    integer :: nlmin = 0
    
    ! Maximum level in the hierarchy
    integer :: nlmax = 0
    
    ! An array of interlevel projection structures for all the levels.
    type(t_interlevelProjectionBlock), dimension(:), pointer :: p_Rprojection

    ! An array of discretisation structures for all the levels
    type(t_mldiscrPointer), dimension(:), pointer :: p_RdiscrPointer
  
    ! Temporary memory for the projection
    type(t_vectorScalar) :: rtempVector
  
  end type

!</typeblock>

  public :: t_interlevelProjectionHier

!</types>
  
  public :: mlprj_initPrjHierarchy
  public :: mlprj_initPrjHierarchyLevel
  public :: mlprj_commitPrjHierarchy
  public :: mlprj_releasePrjHierarchy
  public :: mlprj_performProlongationHier
  public :: mlprj_performRestrictionHier
  public :: mlprj_performInterpolationHier
  
contains

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_initPrjHierarchy (rprjHierarchy,nlmin,nlmax)

!<description>
  ! Initialises a projection hierarchy for projection between levels NLMIN
  ! and NLMAX.
!</description>

!<input>
  ! Minimum level in the hierarchy.
  integer, intent(in) :: nlmin
  
  ! Maximum level in the hierarchy.
  integer, intent(in) :: nlmax
!</input>

!<output>
  ! Projection hierarchy structure to initialise.
  type(t_interlevelProjectionHier), intent(out) :: rprjHierarchy
!</output>

!</subroutine>

    ! Basic initialisation of the structure.
    rprjHierarchy%nlmin = nlmin
    rprjHierarchy%nlmax = nlmax
    
    allocate(rprjHierarchy%p_Rprojection(nlmax-nlmin+1))
    allocate(rprjHierarchy%p_RdiscrPointer(nlmax-nlmin+1))
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_initPrjHierarchyLevel (rprjHierarchy,ilevel,rdiscretisation)

!<description>
  ! Initialises level ilevel of the projection hierarchy using the block
  ! discretisation rdiscretisation.
  ! After the initialisation of all levels, the routine
  ! mlprj_commitProjectionHierarchy must be called!
!</description>

!<input>
  ! Number of the level to initialise
  integer, intent(in) :: ilevel
  
  ! Block discretisation of this level. A pointer to this structure
  ! is written to rprjHierarchy, so the structure must persist until
  ! the hierarchy is released.
  type(t_blockDiscretisation), intent(in),target :: rdiscretisation
!</input>

!<inputoutput>
  ! Projection hierarchy structure.
  type(t_interlevelProjectionHier), intent(inout) :: rprjHierarchy
!</inputoutput>

!</subroutine>

    ! Initialise all levels except for the coarse one, this does not need
    ! initialisation.
    if (ilevel .gt. rprjHierarchy%nlmin) then
      call mlprj_initProjectionDiscr (&
          rprjHierarchy%p_Rprojection(ilevel-rprjHierarchy%nlmin+1),rdiscretisation)
          
      ! Remember the discretisation
      rprjHierarchy%p_RdiscrPointer(ilevel-rprjHierarchy%nlmin+1)%p_rdiscretisation => &
          rdiscretisation
    end if
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_commitPrjHierarchy (rprjHierarchy)

!<description>
  ! Commits a level hierarchy. This routine must be called after the
  ! initialisation of all levels using mlprj_initHierarchyLevel.
  ! It does some final changes to rprjHierarchy to make it ready to use.
!</description>

!<inputoutput>
  ! Projection hierarchy structure.
  type(t_interlevelProjectionHier), intent(inout), target :: rprjHierarchy
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,imaxmem
    type(t_blockDiscretisation), pointer :: p_rdiscrCoarse,p_rdiscrFine

    ! Calculate the amount of temp memory needed for the projection
    ! and allocate a temp vector with that size.
                 
    imaxmem = 1
    do i=rprjHierarchy%nlmin+1,rprjHierarchy%nlmax
      ! Pass the system metrices on the coarse/fine grid to
      ! mlprj_getTempMemoryDirect to specify the discretisation structures
      ! of all equations in the PDE there.
      ! This calculates the overall size of the temp vector needed for the
      ! projection between an arbitrary pair of consecutive levels.
      p_rdiscrCoarse => &
          rprjHierarchy%p_RdiscrPointer(i-rprjHierarchy%nlmin)%p_rdiscretisation
      p_rdiscrFine => &
          rprjHierarchy%p_RdiscrPointer(i-rprjHierarchy%nlmin+1)%p_rdiscretisation
      imaxmem = max(imaxmem,mlprj_getTempMemoryDirect (&
          rprjHierarchy%p_Rprojection(i),&
          p_rdiscrCoarse%RspatialDiscr(1:p_rdiscrCoarse%ncomponents),&
          p_rdiscrFine%RspatialDiscr(1:p_rdiscrFine%ncomponents)))
    end do

    ! Set up a scalar temporary vector for the projection.
    call lsyssc_createVector (rprjHierarchy%rtempVector,imaxmem,.false.)
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_releasePrjHierarchy (rprjHierarchy)

!<description>
  ! Releases a projection hierarchy.
!</description>

!<inputoutput>
  ! Projection hierarchy structure to clean up.
  type(t_interlevelProjectionHier), intent(inout) :: rprjHierarchy
!</inputoutput>

!</subroutine>

    call lsyssc_releaseVector(rprjHierarchy%rtempVector)
    deallocate(rprjHierarchy%p_Rprojection)
  
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine mlprj_performProlongationHier (rprjHierarchy,&
      ifineLevel,rcoarseVector,rfineVector)
  
!<description>
  ! Performs a prolongation for a given block vector (i.e. a projection
  ! in the primal space where the solution lives) based on a projection
  ! hierarchy.
!</description>
  
!<input>
  ! Underlying projection hierarchy structure.
  type(t_interlevelProjectionHier), intent(inout) :: rprjHierarchy

  ! Coarse grid vector
  type(t_vectorBlock), intent(inout) :: rcoarseVector
  
  ! Level of the fine grid vector.
  integer, intent(in) :: ifineLevel
!</input>
  
!<output>
  ! Fine grid vector
  type(t_vectorBlock), intent(inout) :: rfineVector
!</output>
  
!</subroutine>

    call mlprj_performProlongation (&
        rprjHierarchy%P_Rprojection(ifineLevel-rprjHierarchy%nlmin+1),&
        rcoarseVector,rfineVector,rprjHierarchy%rtempVector)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine mlprj_performRestrictionHier (rprjHierarchy,&
      ifineLevel,rcoarseVector,rfineVector)
  
!<description>
  ! Performs a restriction for a given block vector (i.e. a projection
  ! in the dual space where the RHS vector lives) based on a projection
  ! hierarchy.
!</description>
  
!<input>
  ! Underlying projection hierarchy structure.
  type(t_interlevelProjectionHier), intent(inout) :: rprjHierarchy

  ! Fine grid vector
  type(t_vectorBlock), intent(inout) :: rfineVector
  
  ! Level of the fine grid vector.
  integer, intent(in) :: ifineLevel
!</input>
  
!<output>
  ! Coarse grid vector
  type(t_vectorBlock), intent(inout) :: rcoarseVector
!</output>
  
!</subroutine>

    call mlprj_performRestriction (&
        rprjHierarchy%P_Rprojection(ifineLevel-rprjHierarchy%nlmin+1),&
        rcoarseVector,rfineVector,rprjHierarchy%rtempVector)

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine mlprj_performInterpolationHier (rprjHierarchy,&
      ifineLevel,rcoarseVector,rfineVector)
  
!<description>
  ! Performs an interpolation for a given block vector (i.e. a projection
  ! in the primal space where the solution vector lives) based on a projection
  ! hierarchy.
!</description>
  
!<input>
  ! Underlying projection hierarchy structure.
  type(t_interlevelProjectionHier), intent(inout) :: rprjHierarchy

  ! Fine grid vector
  type(t_vectorBlock), intent(inout) :: rfineVector
  
  ! Level of the fine grid vector.
  integer, intent(in) :: ifineLevel
!</input>
  
!<output>
  ! Coarse grid vector
  type(t_vectorBlock), intent(inout) :: rcoarseVector
!</output>
  
!</subroutine>

    call mlprj_performInterpolation (&
        rprjHierarchy%P_Rprojection(ifineLevel-rprjHierarchy%nlmin+1),&
        rcoarseVector,rfineVector,rprjHierarchy%rtempVector)

  end subroutine

end module
