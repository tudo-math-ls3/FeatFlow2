!##############################################################################
!# ****************************************************************************
!# <name> bcassembly </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines to discretise analytically given boundary
!# conditions. Analytically given boundary conditions are 'discretised', i.e.
!# a discrete version (realised by the structure t_discreteBCEntry in the case
!# of BC`s on the real boundary and by the structure t_discreteFBCEntry in 
!# the case of fictitious boundary) is calculated.
!# This structure is used during the solution process to impose the boundary
!# conditions into the solution vector. Therefore, this module contains
!# the bridge how to put analytic boundary conditions into a discrete vector.
!#
!# The module works in tight relationship to the module 'vectorfilters'.
!# While bcassembly provides the functionality to *create* the structures,
!# the module 'vectorfilters' contains routines to *apply* the structure
!# to a given vector.
!# (This separation is necessary to prevent circular dependencies!)
!#
!# The following routines can be found here:
!#
!# 1.) bcasm_initDiscreteBC
!#     -> Initialise a structure collecting discrete boundary conditions.
!#
!# 2.) bcasm_clearDiscreteBC
!#     -> Clear a structure with discrete BC`s. Release memory that is used
!#        by the discrete BC`s but do not destroy the structure itself.
!#
!# 3.) bcasm_releaseDiscreteBC
!#     -> Release a structure with discrete BC`s; deallocates all used memory.
!#
!# 4.) bcasm_newDirichletBC_1D
!#     -> Discretises dirichlet boundary conditions for a 1D discretisation.
!#
!# 5.) bcasm_newDirichletBConRealBd
!#     -> Discretises dirichlet boundary conditions on a 2D boundary region
!#
!# 6.) bcasm_newPdropBConRealBd
!#     -> Discretises pressure drop boundary conditions on a 2D boundary region
!#
!# 7.) bcasm_newDirichletBConMR 
!#     -> Discretises Dirichlet boundary conditions on a mesh region 
!#        (all dimensions)
!#
!# 8.) bcasm_newslipbconrealbd
!#     -> Discretises slip boundary conditions n a 2D boundary region
!#     
!# 9.) bcasm_initDiscreteFBC
!#     -> Initialise a structure collecting discrete fictitious 
!#        boundary boundary conditions.
!#
!# 10.) bcasm_clearDiscreteFBC
!#     -> Clear a structure with discrete FBC`s. Release memory that is used
!#        by the discrete FBC`s but do not destroy the structure itself.
!#
!# 11.) bcasm_releaseDiscreteFBC
!#      -> Release a structure with discrete FBC`s; deallocates all used memory.
!#
!# 12.) bcasm_newDirichletBConFBD
!#      -> Discretises dirichlet boundary conditions on a 2D fictitious boundary domain
!#
!# 14.) bcasm_getVertInBdRegion
!#      -> Determine the vertices in a boundary region
!#
!# 14.) bcasm_getEdgesInBdRegion
!#      -> Determine the edges in a boundary region
!#
!# 14.) bcasm_getElementsInBdRegion
!#      -> Determine the elements in a boundary region
!#
!# Some internal routines:
!#
!# 1.) bcasm_getElementsInBdRegion
!#     Determines elements, edges and vertices in a boundary region
!#
!# 2.) bcasm_releaseDirichlet
!#     -> Release discrete Dirichlet boundary conditions on the real bondary
!#
!# 3.) bcasm_releasePressureDrop
!#     -> Release pressure drop boundary conditions on the real bondary
!#
!# 4.) bcasm_releaseSlip
!#     -> Release slip boundary conditions on the real bondary
!#
!# </purpose>
!##############################################################################

module bcassembly

  use fsystem
  use storage
  use genoutput
  use discretebc
  use discretefbc
  use element
  use triangulation
  use dofmapping
  use mprimitives
  use basicgeometry
  use meshregion
  use spatialdiscretisation
  use sort
  use boundary
  use collection, only: t_collection
  
  implicit none
  
  private

!<constants>

!<constantblock description="Complexity of the discretised BC`s">

  ! Discretise BC`s for implementing them into a defect vector
  integer(I32), parameter, public :: BCASM_DISCFORDEF = 2**0

  ! Discretise BC`s for implementing them into a solution vector
  integer(I32), parameter, public :: BCASM_DISCFORSOL = 2**1

  ! Discretise BC`s for implementing them into a RHS vector
  integer(I32), parameter, public :: BCASM_DISCFORRHS = 2**2

  ! Discretise BC`s for implementing them into a matrix
  integer(I32), parameter, public :: BCASM_DISCFORMAT = 2**3

  ! Discretise BC`s for implementing them into matrix and defect vector
  integer(I32), parameter, public :: BCASM_DISCFORDEFMAT = BCASM_DISCFORDEF + BCASM_DISCFORMAT
  
  ! Discretise BC`s for implementing them into everything
  integer(I32), parameter, public :: BCASM_DISCFORALL = BCASM_DISCFORDEF + BCASM_DISCFORSOL + &
                                                BCASM_DISCFORRHS + BCASM_DISCFORMAT

!</constantblock>

!<constantblock description="Additional options for assembling boundary conditions">

  ! Integral mean values (e.g. for Q1) are computed exactly by the callback routine. 
  ! Without this option, integral mean values are computed by a proper cubature formula.
  integer(I32), parameter, public :: BCASM_DISCOPT_INTBYCB = 2**0

  ! Standard setting of the options during the assembly.
  integer(I32), parameter, public :: BCASM_DISCOPT_STD = 0_I32

!</constantblock>

!<constantblock description="Constants defining the blocking of the assembly">

  ! Number of entries to handle simultaneously (-> number of points or edges,
  ! depending on the situation what to discretise)
#ifndef FBCASM_MAXSIM
#ifndef ENABLE_AUTOTUNE
  integer, parameter, public :: FBCASM_MAXSIM   = 1000
#else
  integer, public            :: FBCASM_MAXSIM   = 1000
#endif
#endif
  
!</constantblock>

!</constants>


!<types>

!<typeblock>
  ! Configuration block for FEAST mirror boundary  conditions. This type of
  ! boundary condition is very special since it is level dependent and needs
  ! special parameters which are not known by the analytical definition.
  type t_configDiscreteFeastMirrorBC
  
    ! Coarsening level. Used when adding additional contributions to the 
    ! matrix/vector. Usually = 0. Must be increased for every level coarser
    ! than the maximum one in a mesh hierarchy.
    real(DP) :: icoarseningLevel = 0
    
    ! Subtype of the FEAST mirror boundary conditions.
    ! =0: apply to matrix and defect vectors.
    ! =1: apply only to matrix; apply Dirichlete BC`s to defect vectors
    integer :: isubtype
    
  end type
  
  public :: t_configDiscreteFeastMirrorBC
  
!</typeblock>

!<typeblock>

  ! Parameter block for level-dependent boundary conditions. This block collects
  ! different parameters for those boundary conditions, which cannot be
  ! represented purely analytically and are discretised by a call to
  ! bcasm_discretiseLevelDepBC.
  type t_levelDependentBC
  
    ! Parameter block for discretising the FEAST mirror boundary conditions.
    type(t_configDiscreteFeastMirrorBC) :: rconfigFeastMirrorBC
  
  end type
  
  public :: t_levelDependentBC

!</typeblock>

!</types>

  public :: bcasm_initDiscreteBC
  public :: bcasm_clearDiscreteBC
  public :: bcasm_releaseDiscreteBC
  public :: bcasm_newDirichletBC_1D
  public :: bcasm_newDirichletBConRealBd
  public :: bcasm_newPdropBConRealBd
  public :: bcasm_newDirichletBConMR 
  public :: bcasm_newSlipBConRealBd
  public :: bcasm_initDiscreteFBC
  public :: bcasm_clearDiscreteFBC
  public :: bcasm_releaseDiscreteFBC
  public :: bcasm_newDirichletBConFBD
  public :: bcasm_getVertInBdRegion
  public :: bcasm_getEdgesInBdRegion
  public :: bcasm_getElementsInBdRegion
  public :: bcasm_releaseDirichlet
  public :: bcasm_releasePressureDrop
  public :: bcasm_releaseSlip

  public :: bcasm_getVertInBCregion  ! Deprecated
  public :: bcasm_getEdgesInBCregion ! Deprecated

  ! Deprecated interface, routine was renamed!
  public :: bcasm_getElementsInBCregion
  interface bcasm_getElementsInBCregion
    module procedure bcasm_getElementsInBdRegion
  end interface
  
contains

! *****************************************************************************
! General routines for discretised boundary conditions
! *****************************************************************************

!<subroutine>

  subroutine bcasm_initDiscreteBC(rdiscreteBC)

!<description>
  ! Initialises a structure that holds information about discretised
  ! boundary conditions.
!</description>

!<output>
  ! Discrete BC structure to be initialised.
  type(t_discreteBC), intent(inout) :: rdiscreteBC
!</output>

!</subroutine>

    if (associated(rdiscreteBC%p_RdiscBCList)) then
      call output_line ('Structure is already initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'bcasm_initDiscreteBC')
      call sys_halt()
    end if

    ! Allocate memory for a first couple of entries.
    call bcasm_newBCentry(rdiscreteBC)
    
    ! Reset the number of used handles.
    rdiscreteBC%inumEntriesUsed = 0

    ! The rest of the structure is initialised by default initialisation.    

  end subroutine

! *****************************************************************************

!<subroutine>

  subroutine bcasm_clearDiscreteBC(rdiscreteBC)

!<description>
  ! Removes all information about discrete BC`s from the rdiscreteBC structure.
  ! Afterwards, the structure is ready to take a new set of discretised
  ! boundary conditions.
!</description>

!<inputoutput>
  ! Discrete BC structure to be emptied.
  type(t_discreteBC), intent(inout) :: rdiscreteBC
!</inputoutput>

!</subroutine>
    
    ! local variables
    integer :: icurrentRegion

    if (associated(rdiscreteBC%p_RdiscBCList)) then
    
      ! Destroy all the substructures in the array.
      do icurrentRegion = 1,size(rdiscreteBC%p_RdiscBCList)  
        
        ! Release all allocated information to this boundary region
        select case (rdiscreteBC%p_RdiscBCList(icurrentRegion)%itype)
        case (DISCBC_TPDIRICHLET)
          ! Discrete Dirichlet boundary conditions. Release the old structure.
          call bcasm_releaseDirichlet(&
            rdiscreteBC%p_RdiscBCList(icurrentRegion)%rdirichletBCs)

        case (DISCBC_TPPRESSUREDROP)
          ! Discrete pressure drop boundary conditions. Release the old structure.
          call bcasm_releasePressureDrop(&
              rdiscreteBC%p_RdiscBCList(icurrentRegion)%rpressureDropBCs)

        case (DISCBC_TPSLIP)
          ! Discrete Slip boundary conditions. Release the old structure.
          call bcasm_releaseSlip( &
                rdiscreteBC%p_RdiscBCList(icurrentRegion)%rslipBCs)

        case (DISCBC_TPFEASTMIRROR)
          ! Discrete FEAST mirror boundary conditions. Release the old structure.
          call bcasm_releaseFeastMirror(&
            rdiscreteBC%p_RdiscBCList(icurrentRegion)%rfeastMirrorBCs)

        end select

        ! BC released, indicate this
        rdiscreteBC%p_RdiscBCList(icurrentRegion)%itype = DISCBC_TPUNDEFINED
        
      end do  
      
      ! Reset the counters
      rdiscreteBC%inumEntriesUsed = 0
    end if
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine bcasm_releaseDiscreteBC (rdiscreteBC)
  
!<description>
  ! This routine cleans up the rdiscreteBC structure. 
  ! All allocated memory is released. 
!</description>

!<inputoutput>
  ! Discrete BC structure to be cleaned up.
  type(t_discreteBC), intent(inout) :: rdiscreteBC
!</inputoutput>

!</subroutine>

    ! Clear the structure
    call bcasm_clearDiscreteBC(rdiscreteBC)
    
    ! Release memory
    if (associated(rdiscreteBC%p_RdiscBCList)) then
      deallocate(rdiscreteBC%p_RdiscBCList)
    end if
    rdiscreteBC%inumEntriesUsed = 0
    rdiscreteBC%inumEntriesAlloc = 0

  end subroutine
      
! *****************************************************************************

!<subroutine>

  subroutine bcasm_initDiscreteFBC(rdiscreteFBC)

!<description>
  ! Initialises a structure that holds information about discretised
  ! boundary conditions.
!</description>

!<output>
  ! Discrete BC structure to be initialised.
  type(t_discreteFBC), intent(inout) :: rdiscreteFBC
!</output>

!</subroutine>

    if (associated(rdiscreteFBC%p_RdiscFBCList)) then
      call output_line ('Structure is already initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'bcasm_initDiscreteFBC')
      call sys_halt()
    end if

    ! Allocate memory for a first couple of entries.
    call bcasm_newFBCentry(rdiscreteFBC)
    
    ! Reset the number of used handles.
    rdiscreteFBC%inumEntriesUsed = 0
    
    ! The rest of the structure is initialised by default initialisation.    

  end subroutine

! *****************************************************************************

!<subroutine>

  subroutine bcasm_clearDiscreteFBC(rdiscreteFBC)

!<description>
  ! Removes all information about discrete BC`s from the rdiscreteBC structure.
  ! Afterwards, the structure is ready to take a new set of discretised
  ! boundary conditions.
!</description>

!<inputoutput>
  ! Discrete BC structure to be emptied.
  type(t_discreteFBC), intent(inout) :: rdiscreteFBC
!</inputoutput>

!</subroutine>
  
    ! local variables
    integer :: icurrentRegion
    
    ! Destroy the content of the structure completely!
    if (associated(rdiscreteFBC%p_RdiscFBCList)) then
    
      ! Destroy all the substructures in the array.
      do icurrentRegion = 1,size(rdiscreteFBC%p_RdiscFBCList)  
        
        ! Release all allocated information to this boundary region
        select case (rdiscreteFBC%p_RdiscFBCList(icurrentRegion)%itype)
        case (DISCBC_TPDIRICHLET)
          ! Discrete Dirichlet boundary conditions. Release the old structure.
          call bcasm_releaseFBCDirichlet(&
            rdiscreteFBC%p_RdiscFBCList(icurrentRegion)%rdirichletFBCs)

        end select

        ! BC released, indicate this
        rdiscreteFBC%p_RdiscFBCList(icurrentRegion)%itype = DISCBC_TPUNDEFINED
        
      end do  
      
    end if
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine bcasm_releaseDiscreteFBC (rdiscreteFBC)
  
!<description>
  ! This routine cleans up the rdiscreteBC structure. 
  ! All allocated memory is released. 
!</description>

!<inputoutput>
  ! Discrete BC structure to be cleaned up.
  type(t_discreteFBC), intent(inout) :: rdiscreteFBC
!</inputoutput>

!</subroutine>

    ! Clear the structure
    call bcasm_clearDiscreteFBC(rdiscreteFBC)
    
    ! Release memory
    deallocate(rdiscreteFBC%p_RdiscFBCList)
    rdiscreteFBC%inumEntriesUsed = 0
    rdiscreteFBC%inumEntriesAlloc = 0

  end subroutine
      
! *****************************************************************************
! Auxiliary routines
! *****************************************************************************

!<subroutine>

  subroutine bcasm_newBCentry(rdiscreteBC, iindex)

!<description>
  ! Creates a new entry for discrete BC in the rdiscreteBC. If necessary,
  ! the list is reallocated. iindex returns the index of the new entry
  ! in the rdiscreteBC%p_RdiscBCList list.
!</description>

!<inputoutput>
  ! Discrete BC structure containing the discrete BC`s.
  type(t_discreteBC), intent(inout) :: rdiscreteBC
!</inputoutput>

!<output>
  ! Optional: Index of the newly created entry for discrete BC`s.
  integer, intent(out), optional :: iindex
!</output>

!</subroutine>

    ! local variables
    integer :: inumAlloc
    type(t_discreteBCEntry), dimension(:), pointer :: p_RdbcList

    ! Allocate some new dbc entries if necessary
    if(rdiscreteBC%inumEntriesAlloc .eq. 0) then
      
      ! The list is currently empty - so allocate it
      allocate(rdiscreteBC%p_RdiscBCList(DISCBC_LISTBLOCKSIZE))
      rdiscreteBC%inumEntriesAlloc = DISCBC_LISTBLOCKSIZE
      
    else if(rdiscreteBC%inumEntriesUsed .ge. rdiscreteBC%inumEntriesAlloc) then
      
      ! We need to reallocate the list as it is full
      inumAlloc = rdiscreteBC%inumEntriesAlloc
      allocate(p_RdbcList(inumAlloc+DISCBC_LISTBLOCKSIZE))
      p_RdbcList(1:inumAlloc) = rdiscreteBC%p_RdiscBCList(1:inumAlloc)
      deallocate(rdiscreteBC%p_RdiscBCList)
      rdiscreteBC%p_RdiscBCList => p_RdbcList
      rdiscreteBC%inumEntriesAlloc = inumAlloc + DISCBC_LISTBLOCKSIZE
      
    end if  
    
    rdiscreteBC%inumEntriesUsed = rdiscreteBC%inumEntriesUsed + 1
    
    if (present(iindex)) &
      iindex = rdiscreteBC%inumEntriesUsed

  end subroutine

! *****************************************************************************

!<subroutine>

  subroutine bcasm_newFBCentry(rdiscreteFBC, iindex)

!<description>
  ! Creates a new entry for discrete BC in the rdiscreteFBC. If necessary,
  ! the list is reallocated. iindex returns the index of the new entry
  ! in the rdiscreteFBC%p_RdiscBCList list.
!</description>

!<inputoutput>
  ! Discrete FBC structure containing the discrete BC`s.
  type(t_discreteFBC), intent(inout) :: rdiscreteFBC
!</inputoutput>

!<output>
  ! Optional: Index of the newly created entry for discrete BC`s.
  integer, intent(out), optional :: iindex
!</output>

!</subroutine>

    ! local variables
    integer :: inumAlloc
    type(t_discreteFBCEntry), dimension(:), pointer :: p_RdbcList

    ! Allocate some new dbc entries if necessary
    if(rdiscreteFBC%inumEntriesAlloc .eq. 0) then
      
      ! The list is currently empty - so allocate it
      allocate(rdiscreteFBC%p_RdiscFBCList(DISCBC_LISTBLOCKSIZE))
      rdiscreteFBC%inumEntriesAlloc = DISCBC_LISTBLOCKSIZE
      
    else if(rdiscreteFBC%inumEntriesUsed .ge. rdiscreteFBC%inumEntriesAlloc) then
      
      ! We need to reallocate the list as it is full
      inumAlloc = rdiscreteFBC%inumEntriesAlloc
      allocate(p_RdbcList(inumAlloc+DISCBC_LISTBLOCKSIZE))
      p_RdbcList(1:inumAlloc) = rdiscreteFBC%p_RdiscFBCList(1:inumAlloc)
      deallocate(rdiscreteFBC%p_RdiscFBCList)
      rdiscreteFBC%p_RdiscFBCList => p_RdbcList
      rdiscreteFBC%inumEntriesAlloc = inumAlloc + DISCBC_LISTBLOCKSIZE
      
    end if  
    
    rdiscreteFBC%inumEntriesUsed = rdiscreteFBC%inumEntriesUsed + 1
    
    if (present(iindex)) &
      iindex = rdiscreteFBC%inumEntriesUsed

  end subroutine

! *****************************************************************************
! Support for boundary conditions on real boundary
! *****************************************************************************

!<subroutine>

  subroutine bcasm_newDirichletBC_1D(rblockDiscr, rdiscreteBC, dleft, dright, iequation)

!<description>
  ! Directly adds dirichlet boundary conditions for a 1D discretisation.
!</description>

!<input>
  ! The underlying block discretisation structure.
  type(t_blockDiscretisation), intent(in) :: rblockDiscr
  
  ! The dirichlet boundary values for the left and right interval ends.
  real(DP), intent(in), optional :: dleft, dright

  ! OPTIONAL: The equation for which the BC's should be applied to.
  ! If not present, this defaults to the first equation.
  integer, optional :: iequation
!</input>  

!<inputoutput>
  ! A t_discreteBC structures, representing the boundary discretised 
  ! in a discretisation-dependent way. The new BC`s are added to this structure.
  type(t_discreteBC), intent(inout) :: rdiscreteBC
!</inputoutput>
  
!</subroutine>

    ! A hand full of local variables
    integer :: idx, ibndVert, ibndElem, iDOF
    integer(I32) :: ielemType
    type(t_discreteBCEntry), pointer :: p_rdiscrBCEntry
    type(t_discreteBCDirichlet), pointer :: p_rdirichlet
    real(DP), dimension(:), pointer             :: p_DdirichletValues
    integer, dimension(:), pointer         :: p_IdirichletDOFs
    type(t_triangulation), pointer              :: p_rtria
    type(t_spatialDiscretisation), pointer      :: p_rspatialDiscr
    type(t_elementDistribution), pointer        :: p_relemDist
    integer, dimension(:), pointer :: p_IelemAtVert,p_IelemAtVertIdx, &
        p_IvertAtBnd, p_IbndCpIdx
    integer, dimension(128) :: IDOFs
    integer, dimension(2) :: IbcIndex
    
    if (rdiscreteBC%inumEntriesAlloc .eq. 0) then
      call output_line ('BC structure not initialised!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'bcasm_newDirichletBC_1D')
      call sys_halt()
    end if
    
    ! Get one or two new entries in the BC structure.
    if (.not. present(dleft) .and. .not. present(dright)) return ! nothing to do
    if (present(dleft))  call bcasm_newBCentry(rdiscreteBC, IbcIndex(1))
    if (present(dright)) call bcasm_newBCentry(rdiscreteBC, IbcIndex(2))
    
    ! Which component is to be discretised?
    p_rspatialDiscr => rblockDiscr%RspatialDiscr(1)
    p_rtria => p_rspatialDiscr%p_rtriangulation
    
    ! Make sure that there is only one FE space in the discretisation.
    if (p_rspatialDiscr%inumFESpaces .ne. 1) then
    
      ! Print an error message
      call output_line ('Spatial discretisation must have 1 FE space!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'bcasm_newDirichletBC_1D')
      
      ! And exit the program
      call sys_halt()
      
    end if
    
    ! Get the element distribution from the spatial discretisation
    p_relemDist => p_rspatialDiscr%RelementDistr(1)
    
    ! Get the element type of the trial functions
    ielemType = elem_getPrimaryElement(p_relemDist%celement)
    
    ! Get the boundary component index array
    call storage_getbase_int(p_rtria%h_IboundaryCpIdx, p_IbndCpIdx)
    
    ! Get the elements that are adjacent to the boundary
    call storage_getbase_int(p_rtria%h_IverticesAtBoundary, p_IvertAtBnd)
    call storage_getbase_int(p_rtria%h_IelementsAtVertexIdx, p_IelemAtVertIdx)
    call storage_getbase_int(p_rtria%h_IelementsAtVertex, p_IelemAtVert)

    ! Go through our both boundary conditions...
    if (present(dleft)) then

      idx = IbcIndex(1)
    
      ! Get the boundary condition entry
      p_rdiscrBCEntry => rdiscreteBC%p_RdiscBCList(idx)
      
      ! We have dirichlet conditions here
      p_rdiscrBCEntry%itype = DISCBC_TPDIRICHLET
      
      ! Get the dirichlet boundary condition
      p_rdirichlet => p_rdiscrBCEntry%rdirichletBCs
      
      ! In the current setting, there is always 1 boundary component
      ! and 1 DOF to be processed.
      p_rdirichlet%nDOF = 1

      ! Initialise the BC's for equation iequation if present.
      p_rdirichlet%icomponent = 1
      if (present(iequation)) &
        p_rdirichlet%icomponent = iequation

      ! Allocate the arrays
      call storage_new('bcasm_newDirichletBC_1D', 'h_IdirichletDOFs', &
          1, ST_INT, p_rdirichlet%h_IdirichletDOFs, ST_NEWBLOCK_NOINIT)
      call storage_new('bcasm_newDirichletBC_1D', 'h_DdirichletValues', &
          1, ST_DOUBLE, p_rdirichlet%h_DdirichletValues, ST_NEWBLOCK_NOINIT)
      
      ! Get the arrays for the dirichlet DOFs and values
      call storage_getbase_int(p_rdirichlet%h_IdirichletDOFs, p_IdirichletDOFs)
      call storage_getbase_double(p_rdirichlet%h_DdirichletValues, p_DdirichletValues)
      
      ! Set the value of the dirichlet condition
      p_DdirichletValues(1) = dleft
      
      ! Get the index of the boundary vertex
      ibndVert = p_IvertAtBnd(1)
      
      ! As we know that a boundary vertex has only one element that is
      ! adjacent to it, we can read it out directly
      ibndElem = p_IelemAtVert(p_IelemAtVertIdx(ibndVert))
      
      ! Now we need to find the DOF for this boundary condition.
      ! This is element-dependent...
      select case(ielemType)
      case (EL_P0_1D)
        ! P0_1D element
        ! This is the easy case - there is only one DOF per element and this
        ! is the one we are searching for...
        call dof_locGlobMapping(p_rspatialDiscr, ibndElem, IDOFs)
        iDOF = IDOFs(1)
      
      case (EL_P1_1D)
        ! P1_1D element
        ! In this case the boundary vertex is one of the two DOFs of the
        ! element - for the left boundary vertex it is the first DOF of
        ! the left boundary element, for the right one it is the second.
        call dof_locGlobMapping(p_rspatialDiscr, ibndElem, IDOFs)
        iDOF = IDOFs(1)
        
      case (EL_P2_1D)
        ! P2_1D element
        ! For the left boundary vertex we need to fix the first DOF, for
        ! the right boundary we need the second DOF.
        call dof_locGlobMapping(p_rspatialDiscr, ibndElem, IDOFs)
        iDOF = IDOFs(1)
      
      case (EL_S31_1D)
        ! S31_1D element
        call dof_locGlobMapping(p_rspatialDiscr, ibndElem, IDOFs)
        iDOF = IDOFs(1)
      
      case (EL_PN_1D)
        ! PN_1D element
        call dof_locGlobMapping(p_rspatialDiscr, ibndElem, IDOFs)
        iDOF = IDOFs(1)
              
      end select
          
      ! Store the DOF
      p_IdirichletDOFs(1) = iDOF
    
    end if 
    
    if (present(dright)) then
    
      idx = IbcIndex(2)
    
      ! Get the boundary condition entry
      p_rdiscrBCEntry => rdiscreteBC%p_RdiscBCList(idx)
      
      ! We have dirichlet conditions here
      p_rdiscrBCEntry%itype = DISCBC_TPDIRICHLET
      
      ! Get the dirichlet boundary condition
      p_rdirichlet => p_rdiscrBCEntry%rdirichletBCs
      
      ! In the current setting, there is always 1 boundary component
      ! and 1 DOF to be processed.
      p_rdirichlet%nDOF = 1

      ! Initialise the BC's for equation iequation if present.
      p_rdirichlet%icomponent = 1
      if (present(iequation)) &
        p_rdirichlet%icomponent = iequation

      ! Allocate the arrays
      call storage_new('bcasm_newDirichletBC_1D', 'h_IdirichletDOFs', &
          1, ST_INT, p_rdirichlet%h_IdirichletDOFs, ST_NEWBLOCK_NOINIT)
      call storage_new('bcasm_newDirichletBC_1D', 'h_DdirichletValues', &
          1, ST_DOUBLE, p_rdirichlet%h_DdirichletValues, ST_NEWBLOCK_NOINIT)
      
      ! Get the arrays for the dirichlet DOFs and values
      call storage_getbase_int(p_rdirichlet%h_IdirichletDOFs, p_IdirichletDOFs)
      call storage_getbase_double(p_rdirichlet%h_DdirichletValues, p_DdirichletValues)
      
      ! Set the value of the dirichlet condition
      p_DdirichletValues(1) = dright
      
      ! Get the index of the boundary vertex
      ibndVert = p_IvertAtBnd(2)
      
      ! As we know that a boundary vertex has only one element that is
      ! adjacent to it, we can read it out directly
      ibndElem = p_IelemAtVert(p_IelemAtVertIdx(ibndVert))
      
      ! Now we need to find the DOF for this boundary condition.
      ! This is element-dependent...
      select case(ielemType)
      case (EL_P0_1D)
        ! P0_1D element
        ! This is the easy case - there is only one DOF per element and this
        ! is the one we are searching for...
        call dof_locGlobMapping(p_rspatialDiscr, ibndElem, IDOFs)
        iDOF = IDOFs(1)
      
      case (EL_P1_1D)
        ! P1_1D element
        ! In this case the boundary vertex is one of the two DOFs of the
        ! element - for the left boundary vertex it is the first DOF of
        ! the left boundary element, for the right one it is the second.
        call dof_locGlobMapping(p_rspatialDiscr, ibndElem, IDOFs)
        iDOF = IDOFs(2)
        
      case (EL_P2_1D)
        ! P2_1D element
        ! For the left boundary vertex we need to fix the first DOF, for
        ! the right boundary we need the second DOF.
        call dof_locGlobMapping(p_rspatialDiscr, ibndElem, IDOFs)
        iDOF = IDOFs(2)

      case (EL_S31_1D)
        ! S31_1D element
        call dof_locGlobMapping(p_rspatialDiscr, ibndElem, IDOFs)
        iDOF = IDOFs(2)
      
      case (EL_PN_1D)
        ! PN_1D element
        call dof_locGlobMapping(p_rspatialDiscr, ibndElem, IDOFs)
        iDOF = IDOFs(2)

      end select
          
      ! Store the DOF
      p_IdirichletDOFs(1) = iDOF
    
    end if
    
    ! That is it
      
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine bcasm_getVertInBCregion (rtriangulation,rboundary,rregion, &
      IminIndex,ImaxIndex,icount)
  
!<description>
  ! DEPRECATED! Not used anymore, since multiply connected boundary regions
  ! must be supported. Use bcasm_getElementsInBdRegion!
  !
  ! This routine receives a boundary region rregion describing a part on the
  ! boundary. According to the triangulation, this boundary region contains
  ! some vertices and edges on the boundary. The routine now figures out,
  ! which index in the IverticesAtBoundary in the triangulation structure 
  ! of the vertices that are inside of this boundary region.
  !
  ! Each boundary region is simply connected, but they may cross the
  ! maximum parameter value. Therefore, we find wither one or two 'sets'
  ! of indices in the IverticesAtBoundary array:
  ! 1.) One simly connected index set:
  !     ..... Vertex Vertex Vertex ......
  !     Here is icount=1 and IminIndex(1) / ImaxIndex(1) gives the
  !     first/last index in IverticesAtBoundary
  ! 2.) Two sets crossing the maximum parameter value:
  !     Vertex Vertex Vertex ........... Vertex Vertex Vertex
  !     Here is icount=2, IminIndex(2) / ImaxIndex(2) gives the
  !     first/last vertex of the first set (which was found at first when
  !     starting to search in the given interval) and IminIndex(1) / ImaxIndex(1)
  !     the first/last vertex if the 2nd set.
!</description>

!<input>
  ! The triangulation structure.
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! The description of the domain boundary
  type(t_boundary), intent(in) :: rboundary
  
  ! A boundary region structure describing a part of a boundary
  type(t_boundaryRegion), intent(in) :: rregion
!</input>

!<output>
  ! The index of the first vertex in IverticesAtBoundary belonging to the
  ! bonudary region rregion. = 0, if no vertices belong to the given region.
  integer, dimension(2), intent(out) :: IminIndex

  ! The index of the last vertex in IverticesAtBoundary belonging to the
  ! bonudary region rregion. = -1, if no vertices belong to the given region.
  integer, dimension(2), intent(out) :: ImaxIndex
  
  ! Number of index sets in IverticesAtBoundary belonging
  ! to rregion. IminIndex(1)..ImaxIndex(1) is the first set,...
  ! IminIndex(icount)..ImaxIndex(icount) is the last set.
  integer, intent(out) :: icount
!</output>

!</subroutine>

  ! local variables
  logical :: binside,bfinish
  real(DP), dimension(:), pointer :: p_DvertexParameterValue
  real(DP), dimension(2) :: dbegin,dend
  real(DP) :: dmaxpar
  integer, dimension(:), pointer :: p_IboundaryCpIdx
  integer :: i, ifoundRegions
  integer, dimension(2) :: Iindex
  
  ! Get the parameter value array from the triangulation
  call storage_getbase_double (rtriangulation%h_DvertexParameterValue, &
                               p_DvertexParameterValue)
                               
  ! Get the array describing the start/end of each boundary component
  call storage_getbase_int (rtriangulation%h_IboundaryCpIdx, &
                            p_IboundaryCpIdx)
  
  ! Get the real min/max parameter value of the boundary region
  dmaxpar = boundary_dgetMaxParVal(rboundary, rregion%iboundCompIdx)
  dbegin(1) = rregion%dminParam
  dend(1) = rregion%dmaxParam
  if (dend(1)-dbegin(1) .ge. dmaxpar) then
    ! Occupy everything
    icount = 1
    dbegin(1) = 0.0_DP
    dend(1) = dmaxpar
    dbegin(2) = 0.0_DP
    dend(2) = -1.0_DP
  else
    ! Calculate the real min/max
    dbegin(1) = mod(rregion%dminParam,dmaxpar)
    dend(1) = mod(rregion%dmaxParam,dmaxpar)
    icount = 1
    if (dend(1) .lt. dbegin(1)) then
      ! probably two sets
      dbegin(2) = dend(1)
      dend(2) = dmaxpar
      dbegin(1) = 0.0_DP
      icount = 2
    end if
  end if
  
  ! All regions in [dbegin(.),dend(.)] are now in asecnding order.
  ! icount indicates the maximum number of regions we expect to have
  ! vertices inside.
  ! Remark: dbegin/dend are not used anymore here - maybe in a later
  ! version if necessary for some reasons...
  !
  ! For now, use a simple linear search. 
  ! The boundary region structure tells us which boundary component
  ! to use; loop through all vertices there.
  
  ifoundRegions = 0
  IminIndex = 0
  ImaxIndex = -1
  binside = .false.
  bfinish = .false.
  
  ! Loop through all vertices in the current boundary component
  
  do I = p_IboundaryCpIdx(rregion%iboundCompIdx), &
         p_IboundaryCpIdx(rregion%iboundCompIdx+1)-1
         
    ! Check if the vertex is inside the region.
    if (boundary_isInRegion (rregion,rregion%iboundCompIdx,p_DvertexParameterValue(I))) then
      ! We are inside
      if (.not. binside) then
        ! We are inside for the first time
        binside = .true.
        IminIndex(ifoundRegions+1) = I
      end if
    else
      ! We are outside - for the first time?
      ! If yes, quit the loop.
      if (binside) then
        binside = .false.
        ImaxIndex(ifoundRegions+1) = I-1
        
        ! We completed the current region successfully
        ifoundRegions = ifoundRegions + 1
        
        ! Decrement icount. If it is still > 0, there is another region
        ! in dbegin/dend that may contain points
        icount = icount - 1
        
        if (icount .le. 0) then
          ! Finish that, we quit the search here as no more regions
          ! are expected.
          icount = ifoundRegions
          return
        end if
      end if
    end if
         
  end do
  
  ! The loop is completed. Question: Are we still inside a region or not?
  ! If yes...
  
  if (binside) then
    ! Save the last vertex number
    ImaxIndex(ifoundRegions+1) = p_IboundaryCpIdx(rregion%iboundCompIdx+1)-1
    ! Complete the region and finish.
    icount = ifoundRegions + 1
  else 
    ! No we are not. So we were awaiting points for another part of the region
    ! that never came! Reset the last index pair and quit.
    IminIndex(ifoundRegions+1) = 0
    ImaxIndex(ifoundRegions+1) = -1
    icount = ifoundRegions
  end if
  
  if (icount .eq. 2) then
    ! We found the intervals in the 'wrong order'. When going through the boundary
    ! starting at a parameter value x, we would first find the 2nd segment
    ! x <= y..TMAX, then the first one 0..z. on the other hand set up IminIndex/
    ! ImaxIndex in the other order. So exchange IxxxIndex(1) with IxxxIndex(2)
    ! to get the right order again.
    Iindex(1) = IminIndex(2)
    IminIndex(2) = IminIndex(1)
    IminIndex(1) = Iindex(1)

    Iindex(2) = ImaxIndex(2)
    ImaxIndex(2) = ImaxIndex(1)
    ImaxIndex(1) = Iindex(2)
  end if
  
  end subroutine
    
  ! ***************************************************************************
  
!<subroutine>

  subroutine bcasm_getEdgesInBCregion (rtriangulation,rboundary,rregion, &
      IminIndex,ImaxIndex,icount)
  
!<description>
  ! DEPRECATED! Not used anymore, since multiply connected boundary regions
  ! must be supported. Use bcasm_getElementsInBdRegion!
  !
  ! This routine receives a boundary region rregion describing a part on the
  ! boundary. According to the triangulation, this boundary region contains
  ! some vertices and edges on the boundary. The routine now figures out,
  ! which  index in the IverticesAtBoundary in the triangulation structure 
  ! of the vertices that are inside of this boundary region.
  !
  ! Each boundary region is simply connected, but they may cross the
  ! maximum parameter value. Therefore, we find wither one or two 'sets'
  ! of indices in the IverticesAtBoundary array:
  ! 1.) One simly connected index set:
  !     ..... Edge Edge Edge......
  !     Here is icount=1 and IminIndex(1) / ImaxIndex(1) gives the
  !     first/last index in IverticesAtBoundary
  ! 2.) Two sets crossing the maximum parameter value:
  !     Edge Edge Edge ........... Edge Edge Edge
  !     Here is icount=2, IminIndex(2) / ImaxIndex(2) gives the
  !     first/last vertex of the first set (which was found at first when
  !     starting to search in the given interval) and IminIndex(1) / ImaxIndex(1)
  !     the first/last vertex if the 2nd set.
!</description>

!<input>
  ! The triangulation structure.
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! The description of the domain boundary
  type(t_boundary), intent(in) :: rboundary
  
  ! A boundary region structure describing a part of a boundary
  type(t_boundaryRegion), intent(in) :: rregion
!</input>

!<output>
  ! The index of the first edge in IedgesAtBoundary belonging to the
  ! bonudary region rregion. = 0, if no edges belong to the given region.
  integer, dimension(2), intent(out) :: IminIndex

  ! The index of the last edges in IedgesAtBoundary belonging to the
  ! bonudary region rregion. = -1, if no edges belong to the given region.
  integer, dimension(2), intent(out) :: ImaxIndex
  
  ! Number of index sets in IedgesAtBoundary belonging
  ! to rregion. IminIndex(1)..ImaxIndex(1) is the first set,...
  ! IminIndex(icount)..ImaxIndex(icount) is the last set.
  integer, intent(out) :: icount
!</output>

!</subroutine>

  ! local variables
  logical binside,bfinish
  real(DP), dimension(:), pointer :: p_DedgeParameterValue
  real(DP), dimension(2) :: dbegin,dend
  real(DP) :: dmaxpar
  integer, dimension(:), pointer :: p_IboundaryCpIdx
  integer :: i, ifoundRegions
  integer, dimension(2) :: Iindex
  
  ! Get the parameter value array from the triangulation
  call storage_getbase_double (rtriangulation%h_DedgeParameterValue, &
                               p_DedgeParameterValue)
                               
  ! Get the array describing the start/end of each boundary component
  call storage_getbase_int (rtriangulation%h_IboundaryCpIdx, &
                            p_IboundaryCpIdx)
  
  ! Get the real min/max parameter value of the boundary region
  dmaxpar = boundary_dgetMaxParVal(rboundary, rregion%iboundCompIdx)
  dbegin(1) = rregion%dminParam
  dend(1) = rregion%dmaxParam
  if (dend(1)-dbegin(1) .ge. dmaxpar) then
    ! Occupy everything
    icount = 1
    dbegin(1) = 0.0_DP
    dend(1) = dmaxpar
    dbegin(2) = 0.0_DP
    dend(2) = -1.0_DP
  else
    ! Calculate the real min/max
    dbegin(1) = mod(rregion%dminParam,dmaxpar)
    dend(1) = mod(rregion%dmaxParam,dmaxpar)
    icount = 1
    if (dend(1) .lt. dbegin(1)) then
      ! probably two sets
      dbegin(2) = dend(1)
      dend(2) = dmaxpar
      dbegin(1) = 0.0_DP
      icount = 2
    end if
  end if
  
  ! All regions in [dbegin(.),dend(.)] are now in asecnding order.
  ! icount indicates the maximum number of regions we expect to have
  ! vertices inside.
  ! Remark: dbegin/dend are not used anymore here - maybe in a later
  ! version if necessary for some reasons...
  !
  ! For now, use a simple linear search. 
  ! The boundary region structure tells us which boundary component
  ! to use; loop through all vertices there.
  
  ifoundRegions = 0
  IminIndex = 0
  ImaxIndex = -1
  binside = .false.
  bfinish = .false.
  
  ! Loop through all vertices in the current boundary component
  
  do I = p_IboundaryCpIdx(rregion%iboundCompIdx), &
         p_IboundaryCpIdx(rregion%iboundCompIdx+1)-1
         
    ! Check if the edge is inside the region.
    if (boundary_isInRegion (rregion,rregion%iboundCompIdx,p_DedgeParameterValue(I))) then
      ! We are inside
      if (.not. binside) then
        ! We are inside for the first time
        binside = .true.
        IminIndex(ifoundRegions+1) = I
      end if
    else
      ! We are outside - for the first time?
      ! If yes, quit the loop.
      if (binside) then
        binside = .false.
        ImaxIndex(ifoundRegions+1) = I-1
        
        ! We completed the current region successfully
        ifoundRegions = ifoundRegions + 1
        
        ! Decrement icount. If it is still > 0, there is another region
        ! in dbegin/dend that may contain points
        icount = icount - 1
        
        if (icount .le. 0) then
          ! Finish that, we quit the search here as no more regions
          ! are expected.
          icount = ifoundRegions
          return
        end if
      end if
    end if
         
  end do
  
  ! The loop is completed. Question: Are we still inside a region or not?
  ! If yes...
  
  if (binside) then
    ! Save the last edge number
    ImaxIndex(ifoundRegions+1) = p_IboundaryCpIdx(rregion%iboundCompIdx+1)-1
    ! Complete the region and finish.
    icount = ifoundRegions + 1
  else 
    ! No we are not. So we were awaiting points for another part of the region
    ! that never came! Reset the last index pair and quit.
    IminIndex(ifoundRegions+1) = 0
    ImaxIndex(ifoundRegions+1) = -1
    icount = ifoundRegions
  end if

  if (icount .eq. 2) then
    ! We found the intervals in the 'wrong order'. When going through the boundary
    ! starting at a parameter value x, we would first find the 2nd segment
    ! x <= y..TMAX, then the first one 0..z. on the other hand set up IminIndex/
    ! ImaxIndex in the other order. So exchange IxxxIndex(1) with IxxxIndex(2)
    ! to get the right order again.
    Iindex(1) = IminIndex(2)
    IminIndex(2) = IminIndex(1)
    IminIndex(1) = Iindex(1)

    Iindex(2) = ImaxIndex(2)
    ImaxIndex(2) = ImaxIndex(1)
    ImaxIndex(1) = Iindex(2)
  end if

  end subroutine
    
  ! ***************************************************************************
  
!<subroutine>

  subroutine bcasm_getVertInBdRegion(rtriangulation,rboundaryRegion, &
      ncount, Ivertices, Itemp)
  
!<description>
  ! Calculates the vertices in a boundary region.
!</description>

!<input>
  ! The triangulation structure.
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! A boundary region structure describing a part of a boundary
  type(t_boundaryRegion), intent(in) :: rboundaryRegion

  ! Temporary array. Must be large enough to hold all vertices in the region.
  ! Must be specified if Iedges is present!
  integer, dimension(:), intent(inout), optional :: Itemp
!</input>

!<output>
  ! Number of vertices in the region.
  ! If Iedges is not present: Estimated number of vertices in the boundary region;
  ! the number is larger or equal than the actual numer of vertices in the region.
  integer, intent(out) :: ncount
  
  ! A list of edges in the region. The array must be large enough to hold
  ! all vertices in the region. 
  ! If not specified, the routine only calculates ncount, the number
  ! of edges on the boundary.
  integer, dimension(:), intent(out), optional :: Ivertices 
!</output>

!</subroutine>

    ! local variables
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer :: i
    
    if (present(Ivertices) .and. .not. present (Itemp)) then
      call output_line ('Itemp not specified!', &
          OU_CLASS_ERROR,OU_MODE_STD,'bcasm_getVertInBdRegion')
      call sys_halt()
    end if
    
    ! Get ncount
    call bcasm_getElementsInBdRegion (rtriangulation,rboundaryRegion, ncount, &
        Ivertices,IvtLocal=Itemp)

    ! Calculate the edge numbers.
    if (present(Ivertices)) then
    
      call storage_getbase_int2d(rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
    
      do i=1,ncount
        Ivertices(i) = p_IverticesAtElement(Itemp(i),Ivertices(i))
      end do
      
    end if

  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine bcasm_getEdgesInBdRegion(rtriangulation,rboundaryRegion, ncount, Iedges, Itemp)
  
!<description>
  ! Calculates the edges in a boundary region.
!</description>

!<input>
  ! The triangulation structure.
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! A boundary region structure describing a part of a boundary
  type(t_boundaryRegion), intent(in) :: rboundaryRegion

  ! Temporary array. Must be large enough to hold all edges in the region.
  ! Must be specified if Iedges is present!
  integer, dimension(:), intent(inout), optional :: Itemp
!</input>

!<output>
  ! Number of edges in the region.
  ! If Iedges is not present: Estimated number of edges in the boundary region;
  ! the number is larger or equal than the actual numer of edges in the region.
  integer, intent(out) :: ncount
  
  ! A list of edges in the region. The array must be large enough to hold
  ! all edges in the region. 
  ! If not specified, the routine only calculates ncount, the number
  ! of edges on the boundary.
  integer, dimension(:), intent(out), optional :: Iedges
!</output>

!</subroutine>

    ! local variables
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer :: i
    
    if (present(Iedges) .and. .not. present (Itemp)) then
      call output_line ('Itemp not specified!', &
          OU_CLASS_ERROR,OU_MODE_STD,'bcasm_getEdgesInBdRegion')
      call sys_halt()
    end if
    
    ! Get ncount
    call bcasm_getElementsInBdRegion (rtriangulation,rboundaryRegion, ncount, &
        Iedges,IedgeLocal=Itemp)

    ! Calculate the edge numbers.
    if (present(Iedges)) then
    
      call storage_getbase_int2d(rtriangulation%h_IedgesAtElement,p_IedgesAtElement)
    
      do i=1,ncount
        Iedges(i) = p_IedgesAtElement(Itemp(i),Iedges(i))
      end do
      
    end if

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>

  subroutine bcasm_getElementsInBdRegion (rtriangulation,rregion, ncount, &
                                          IelList,IelListIdx,IvtLocal,IedgeLocal)
  
!<description>
  ! This routine receives a boundary region rregion describing a part on the
  ! boundary. According to the triangulation, this boundary region contains
  ! some vertices and edges on the boundary. 
  !
  ! The routine will search all elements that touch (with a vertex or an edge)
  ! the boundary region. All these elements are written to IelList.
  ! If a vertex of the element touches it, the corresponding entry in IvtLocal 
  ! is set to the local vertex number of that vertex, otherwise it is set to to 0.
  ! If an edge of the element touches it, the corresponding entry in IvtLocal 
  ! is set to the local edge number of that edge, otherwise it is set to to 0.
  !
  ! If neighter IvtLocal nor IedgeLocal is present, the routine calculates
  ! the maximum number of elements that the boundary region covers, including
  ! all edges and vertices.
!</description>

!<input>
  ! The triangulation structure.
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! A boundary region structure describing a part of a boundary
  type(t_boundaryRegion), intent(in) :: rregion
!</input>

!<output>
  ! Number of edges in the boundary region, including repetitions.
  integer, intent(out) :: ncount

  ! OPTIONAL: Array that receives a list of elements that touch the boundary region.
  ! Elements may be inside here more than once.
  integer, dimension(:), intent(out), optional :: IelList

  ! OPTIONAL: Index list. Receives for every element touching the BC region an index
  ! into the IelementsAtBoundary where the element can be found.
  integer, dimension(:), intent(out), optional :: IelListIdx
  
  ! OPTIONAL: For each element on the boundary, local index of a vertex touching the
  ! boundary region. If more than one vertex touches it, the element is
  ! repeated in IelList and IvertexList contains all vertices that touch.
  ! If not present, elements that touch the boundary region only with a point
  ! are ignored.
  ! If neighter IvtLocal nor IedgeLocal is present, the routine calculates
  ! the maximum number of elements that the boundary region covers, including
  ! all edges and vertices.
  integer, dimension(:), intent(out), optional :: IvtLocal

  ! OPTIONAL: For each element on the boundary, local index of an edge touching the
  ! boundary region. If more than one edge touches it, the element is
  ! repeated in IelList and IedgeList contains all edges that touch.
  ! If not present, elements that touch the boundary region only with an edge
  ! are ignored.
  ! If neighter IvtLocal nor IedgeLocal is present, the routine calculates
  ! the maximum number of elements that the boundary region covers, including
  ! all edges and vertices.
  integer, dimension(:), intent(out), optional :: IedgeLocal
  
!</output>

!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: p_DvertexParameterValue
    real(DP), dimension(:), pointer :: p_DedgeParameterValue
    integer, dimension(:), pointer :: p_IelementsAtBoundary, &
               p_IverticesAtBoundary, p_IedgesAtBoundary
    integer, dimension(:,:), pointer :: p_IverticesAtElement,p_IedgesAtElement
    integer, dimension(:), pointer :: p_IboundaryCpIdx
    integer :: i,iidx,iel,ivt,iedge
    logical :: bvertexInside, bedgeInside,bcheckall
    
    ! Get the parameter value array from the triangulation
    call storage_getbase_double (rtriangulation%h_DedgeParameterValue, &
                                p_DedgeParameterValue)
    call storage_getbase_double (rtriangulation%h_DvertexParameterValue, &
                                p_DvertexParameterValue)
    call storage_getbase_int (rtriangulation%h_IelementsAtBoundary, &
                                p_IelementsAtBoundary)
    call storage_getbase_int (rtriangulation%h_IverticesAtBoundary, &
                                p_IverticesAtBoundary)
    call storage_getbase_int (rtriangulation%h_IedgesAtBoundary, &
                                p_IedgesAtBoundary)
    call storage_getbase_int2d (rtriangulation%h_IverticesAtElement, &
                                p_IverticesAtElement)
    call storage_getbase_int2d (rtriangulation%h_IedgesAtElement, &
                                p_IedgesAtElement)
                                 
    ! Get the array describing the start/end of each boundary component
    call storage_getbase_int (rtriangulation%h_IboundaryCpIdx, p_IboundaryCpIdx)
    
    ncount = 0
    
    ! Cancel if this region does not exist.
    if ((rregion%iboundCompIdx .lt. 1) .or. &
        (rregion%iboundCompIdx .gt. rtriangulation%nbct)) then
      return
    end if

    bvertexInside = .false.
    bedgeInside = .false.
    
    ! If both destination arrays are missing, calculate the maximum number
    ! of elements in the region.
    bcheckall = .false.
    if (.not. present(IvtLocal) .and. .not. present(IedgeLocal)) &
        bcheckall = .true.
    
    ! Loop through all elements on the boundary. If we find any in the boundary
    ! region, take it!
    elementloop: do i=p_IboundaryCpIdx(rregion%iboundCompIdx),p_IboundaryCpIdx(rregion%iboundCompIdx+1)-1
    
      ! Check if the vertex or edge touches the boundary region
      if (present(IvtLocal) .or. bcheckall) &
        bvertexInside = boundary_isInRegion (rregion,rregion%iboundCompIdx,p_DvertexParameterValue(i))
        
      if (present(IedgeLocal) .or. bcheckall) &
        bedgeInside = boundary_isInRegion (rregion,rregion%iboundCompIdx,p_DedgeParameterValue(i))
      
      if (bvertexInside .or. bedgeInside) then
      
        ! Got one. Add the element
        ncount = ncount + 1
        iel = p_IelementsAtBoundary(i)
        ivt = p_IverticesAtBoundary(i)
        iedge = p_IedgesAtBoundary(i)
        if (present(IelList))    IelList(ncount) = iel
        if (present(IelListIdx)) IelListIdx(ncount) = i
        if (present(IvtLocal))   IvtLocal(ncount) = 0
        if (present(IedgeLocal)) IedgeLocal(ncount) = 0
        
        ! Remember the vertex and/or edge
        if (bvertexInside .and. present(IvtLocal)) then
        
          ! Search the local number of the vertex that we found
          do iidx = 1,ubound(p_IverticesAtElement,1)
            if (p_IverticesAtElement(iidx,iel) .eq. ivt) then
          
              IvtLocal(ncount) = iidx
              
              if (bedgeInside) then
                ! The edge is also inside -- and has the same local number.
                IedgeLocal(ncount) = iidx
              end if
              
              ! Next element
              cycle elementloop
            end if
          end do
          
          ! Oops, the triangulation is destroyed!
          call output_line ('Local vertex number not found!', &
              OU_CLASS_ERROR,OU_MODE_STD,'bcasm_getElementsInBdRegion')
          call sys_halt()
          
        end if

        ! This point is only reached, if the vertex is not inside.
        ! Check the edge.
        if (bedgeInside .and. present(IedgeLocal)) then
          ! Search the local number of the edge that we found
          do iidx = 1,ubound(p_IedgesAtElement,1)
            if (p_IedgesAtElement(iidx,iel) .eq. iedge) then
            
              IedgeLocal(ncount) = iidx
              
              ! Next element
              cycle elementloop
            end if
          end do
          
          ! Oops, the triangulation is destroyed!
          call output_line ('Local edge number not found!', &
              OU_CLASS_ERROR,OU_MODE_STD,'bcasm_getElementsInBdRegion')
          call sys_halt()
        end if
        
      end if
    end do elementloop
  
  end subroutine
    
  ! ***************************************************************************

!<subroutine>

  subroutine bcasm_newDirichletBConRealBd (rblockDiscretisation, &
      iequation, rboundaryRegion, rdiscreteBC, &
      fgetBoundaryValues,rcollection, ccomplexity, coptions)
  
!<description>
  ! Creates a discrete version of Dirichlet boundary conditions.
  ! rboundaryRegion describes the region which is to be discretised. The discretised
  ! boundary conditions are created in rdiscreteBC, which is assumed
  ! to be undefined when entering this routine.
!</description>

!<input>
  ! The discretisation structure of the underlying discretisation. The boundary
  ! conditions inside of this structure are discretised.
  type(t_blockDiscretisation), intent(in), target :: rblockDiscretisation

  ! An identifier for the equation, this boundary condition refers to.
  ! >= 1. 1=first equation (e.g. X-velocity), 2=2nd equation (e.g. 
  ! Y-velocity), etc.
  integer, intent(in) :: iequation

  ! A boundary-condition-region object, describing the position on the
  ! boundary where boundary conditions should be imposed.
  type(t_boundaryRegion), intent(in) :: rboundaryRegion

  ! A callback function that calculates values on the boundary.
  ! Is declared in the interface include file 'intf_bcassembly.inc'.
  include 'intf_bcassembly.inc'
  
  ! Optional: A collection structure to inform the callback function with
  ! additional information. 
  type(t_collection), intent(inout), optional :: rcollection

  ! Optional: A combination of BCASM_DISCFORxxx constants that specify
  ! the complexity of the discretisation that is to perform. This allows to
  ! discretise only parts of the BC`s, e.g. only setting up those
  ! information that are necessary for filtering defect vectors.
  ! If not specified, BCASM_DISCFORALL is assumed, i.e. the resulting
  ! boundary conditions can be used for everything.
  integer(I32), intent(in), optional :: ccomplexity
  
  ! Optional: A field specifying additional options for the assembly.
  ! If not specified, BCASM_DISCOPT_STD is assumed.
  integer(I32), intent(in), optional :: coptions
!</input>  

!<inputoutput>
  ! A t_discreteBC structures, representing the boundary discretised 
  ! in a discretisation-dependent way. The new BC`s are added to this structure.
  type(t_discreteBC), intent(inout) :: rdiscreteBC
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,i2,j,ilocalEdge,icount,ielidx
    integer(I32) :: celement
    integer :: ielement
    integer :: iedge,ipoint1,ipoint2,NVT
    integer, dimension(1) :: Icomponents
    type(t_discreteBCDirichlet),pointer         :: p_rdirichletBCs
    type(t_triangulation), pointer              :: p_rtriangulation
    type(t_spatialDiscretisation), pointer      :: p_rspatialDiscr
    integer, dimension(:), pointer              :: p_IelementDistr
    integer, dimension(:,:), allocatable :: Idofs
    real(DP), dimension(:,:), allocatable       :: DdofValue
    real(DP), dimension(:), pointer             :: p_DedgeParameterValue,p_DvertexParameterValue
    real(DP), dimension(:), pointer             :: p_DdirichletValues
    integer, dimension(:,:), pointer            :: p_IedgesAtElement
    integer, dimension(:,:), pointer            :: p_IverticesAtElement
    integer, dimension(:), pointer              :: p_IdirichletDOFs
    integer, dimension(:), pointer              :: p_IboundaryCpIdx
    integer, dimension(:), allocatable          :: IverticesAtBoundaryIdx
    integer, dimension(:), allocatable          :: IedgesAtBoundaryIdx
    integer, dimension(:), allocatable          :: IelementsAtBoundary
    integer, dimension(:), allocatable          :: IelementsAtBoundaryIdx
    
    real(DP), dimension(EL_MAXNDER)            :: Dvalues
    
    real(DP) :: dpar,dpar1,dpar2,dval,dval1,dval2
    integer :: nve,nnve,nvbd
    
    integer ::iidx
    type(t_discreteBCEntry), pointer :: p_rdiscreteBCentry
    
    ! Position of cubature points for 2-point Gauss formula on an edge.
    ! Used for Q2T. 
    real(DP), parameter :: Q2G1 = -0.577350269189626_DP !-SQRT(1.0_DP/3.0_DP)
    real(DP), parameter :: Q2G2 =  0.577350269189626_DP ! SQRT(1.0_DP/3.0_DP)
    
    ! List of element distributions in the discretisation structure
    type(t_elementDistribution), dimension(:), pointer :: p_RelementDistribution

    integer(I32) :: casmComplexity, cmyoptions
    
    casmComplexity = BCASM_DISCFORALL
    if (present(ccomplexity)) casmComplexity = ccomplexity
    
    cmyoptions = BCASM_DISCOPT_STD
    if (present(coptions)) cmyoptions = coptions
    
    if ((iequation .lt. 1) .or. &
        (iequation .gt. size(rblockDiscretisation%RspatialDiscr))) then
      call output_line (&
          'Specified equation number is out of range:'//trim(sys_siL(iequation,10)), &
          OU_CLASS_ERROR,OU_MODE_STD,'bcasm_newDirichletBConRealBd')
      call sys_halt()
    end if
    
    ! Which component is to be discretised?
    p_rspatialDiscr => rblockDiscretisation%RspatialDiscr(iequation)

    ! For easier access:
    p_rtriangulation => p_rspatialDiscr%p_rtriangulation
    call storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
    call storage_getbase_int2D(p_rtriangulation%h_IedgesAtElement,p_IedgesAtElement)
    call storage_getbase_int (p_rtriangulation%h_IboundaryCpIdx, p_IboundaryCpIdx)

    p_RelementDistribution => p_rspatialDiscr%RelementDistr
    
    ! The parameter value arrays may not be initialised.
    if (p_rtriangulation%h_DedgeParameterValue .ne. ST_NOHANDLE) then
      call storage_getbase_double(p_rtriangulation%h_DedgeParameterValue,&
          p_DedgeParameterValue)
    else
      nullify(p_DedgeParameterValue)
    end if

    if (p_rtriangulation%h_DvertexParameterValue .ne. ST_NOHANDLE) then
      call storage_getbase_double(p_rtriangulation%h_DvertexParameterValue,&
          p_DvertexParameterValue)
    else
      nullify(p_DvertexParameterValue)
    end if

    NVT = p_rtriangulation%NVT
    nnve = p_rtriangulation%NNVE

    if (p_rspatialDiscr%ccomplexity .ne. SPDISC_UNIFORM) then
      ! Every element can be of different type.
      call storage_getbase_int(p_rspatialDiscr%h_IelementDistr,&
          p_IelementDistr)
    else
      ! All elements are of the samne type. Get it in advance.
      celement = p_rspatialDiscr%RelementDistr(1)%celement
      nve = elem_igetNVE (celement)
    end if
    
    ! Get a new BC entry
    call bcasm_newBCentry(rdiscreteBC, iidx)
    p_rdiscreteBCentry => rdiscreteBC%p_RdiscBCList(iidx)
    
    ! We have Dirichlet boundary conditions
    p_rdiscreteBCentry%itype = DISCBC_TPDIRICHLET
    
    ! Fill the structure for discrete Dirichlet BC`s in the
    ! t_discreteBCEntry structure
    p_rdirichletBCs => p_rdiscreteBCentry%rdirichletBCs
    
    p_rdirichletBCs%icomponent = iequation
    Icomponents(1) = iequation
    
    ! We have to deal with all DOF`s on the boundary. This is highly element
    ! dependent and therefore a little bit tricky :(
    !
    ! As we are in 2D, we can use parameter values at first to figure out,
    ! which points and which edges are on the boundary.
    ! What we have is a boundary segment. Now ask the boundary-index routine
    ! to give us the vertices and edges on the boundary that belong to 
    ! this boundary segment.
    
    allocate(IverticesAtBoundaryIdx(p_rtriangulation%NVBD))
    allocate(IedgesAtBoundaryIdx(p_rtriangulation%NVBD))
    allocate(IelementsAtBoundary(p_rtriangulation%NVBD))
    allocate(IelementsAtBoundaryIdx(p_rtriangulation%NVBD))
    
    call bcasm_getElementsInBdRegion (p_rtriangulation,rboundaryRegion, &
        icount, IelementsAtBoundary, IelementsAtBoundaryIdx, &
        IverticesAtBoundaryIdx,IedgesAtBoundaryIdx)
                                   
    if (icount .eq. 0) then
      deallocate(IverticesAtBoundaryIdx)
      deallocate(IedgesAtBoundaryIdx)
      deallocate(IelementsAtBoundary)
      deallocate(IelementsAtBoundaryIdx)
      return
    end if
                                   
    ! Reserve some memory to save temporarily all DOF`s of all boundary
    ! elements.
    ! We handle all boundary elements simultaneously - let us hope that there are 
    ! never so many elements on the boundary that our memory runs out :-)
    allocate (Idofs(EL_MAXNBAS,icount))
    
    ! If the boundary conditions are only to be discretised for modifying the
    ! defect vector, we can skip calculating the values in the boundary.
    if (iand(casmComplexity,not(BCASM_DISCFORDEFMAT)) .ne. 0) then
      allocate (DdofValue(EL_MAXNBAS,icount))
    end if
    Idofs = 0

    ! Now the elements with indices iminidx..imaxidx in the ItrialElements
    ! of the triangulation are on the boundary. Some elements may appear
    ! twice (on edges e.g.) but we do not care.
    !
    ! Ask the DOF-mapping routine to get us those DOF`s belonging to elements
    ! on the boundary.
    !
    ! The 'mult' call only works on uniform discretisations. We cannot assume
    ! that and have to call dof_locGlobMapping for every element separately.
    if (p_rspatialDiscr%ccomplexity .eq. SPDISC_UNIFORM) then
      call dof_locGlobMapping_mult(p_rspatialDiscr, &
                IelementsAtBoundary(1:icount), Idofs)
    else
      do ielement = 1,icount
        call dof_locGlobMapping(p_rspatialDiscr, IelementsAtBoundary(ielement),&
            Idofs(:,ielement))
      end do
    end if
                   
    ! Loop through the elements
    do ielidx = 1,icount

      ! Get the element and information about it.
      ielement = IelementsAtBoundary (ielidx)
      
      ! Index in the boundary arrays.
      I = IelementsAtBoundaryIdx (ielidx)
      
      ! Get the element type in case we do not have a uniform triangulation.
      ! Otherwise, celement was set to the trial element type above.
      if (p_rspatialDiscr%ccomplexity .ne. SPDISC_UNIFORM) then
        celement = p_RelementDistribution(p_IelementDistr(ielement))%celement
        nve = elem_igetNVE (celement)
      end if
        
      ilocaledge = 0
      
      if (IverticesAtBoundaryIdx(ielidx) .ne. 0) then
        ! Get the local index of the edge -- it coincides with the local index
        ! of the vertex.
        ilocaledge = IverticesAtBoundaryIdx(ielidx)
        ipoint1 = p_IverticesAtElement(IverticesAtBoundaryIdx(ielidx),ielement)
        ipoint2 = p_IverticesAtElement(mod(IverticesAtBoundaryIdx(ielidx),nve)+1,ielement)
      else
        ipoint1 = 0
        ipoint2 = 0
      end if

      if (IedgesAtBoundaryIdx(ielidx) .ne. 0) then
        ! Get the local index of the edge -- it coincides with the local index
        ! of the vertex.
        ilocaledge = IedgesAtBoundaryIdx(ielidx)

        ! Get the edge
        iedge = p_IedgesAtElement(IedgesAtBoundaryIdx(ielidx),ielement)
      else
        iedge = 0
      end if
      
      dpar = -1.0_DP
    
      ! Now the element-dependent part. For each element type, we have to
      ! figure out which DOF`s are on the boundary!
      !
      ! We proceed as follows: We figure out, which DOF is on the
      ! boundary. Then, we ask our computation routine to calculate
      ! the necessary value and translate them into a DOF value.
      ! All DOF values are collected later.
      select case (elem_getPrimaryElement(celement))
      
      case (EL_P0,EL_Q0)
        ! Nice element, only one DOF :-)
        ! Either the edge or an adjacent vertex is on the boundary.
        
        ! If parameter values are available, get the parameter value.
        ! Otherwise, take the standard value from above!
        if (associated(p_DvertexParameterValue)) &
          dpar = p_DvertexParameterValue(I)
      
        ! Get the value at the corner point and accept it as
        ! Dirichlet value.
        call fgetBoundaryValues (Icomponents,p_rspatialDiscr,&
                                rboundaryRegion,ielement, DISCBC_NEEDFUNC,&
                                ipoint1,dpar, Dvalues, rcollection)
                                  
        if (iand(casmComplexity,not(BCASM_DISCFORDEFMAT)) .ne. 0) then
        
          ! Dvalues(1) gives us the function value in the point. Save it
          DdofValue(1,ielidx) = Dvalues(1)
          
        end if
        
        ! A value of SYS_INFINITY indicates a do-nothing node inside of
        ! Dirichlet boundary.
        if (Dvalues(1) .ne. SYS_INFINITY) then
          ! Set the DOF number < 0 to indicate that it is Dirichlet
          Idofs(1,ielidx) = -abs(Idofs(1,ielidx))
        end if
        
      case (EL_P1,EL_Q1)

        ! Left point inside? -> Corresponding DOF must be computed
        if ( ipoint1 .ne. 0 ) then
        
          ! If parameter values are available, get the parameter value.
          ! Otherwise, take the standard value from above!
          if (associated(p_DvertexParameterValue)) &
            dpar = p_DvertexParameterValue(I)
                
          call fgetBoundaryValues (Icomponents,p_rspatialDiscr,&
                                  rboundaryRegion,ielement, DISCBC_NEEDFUNC,&
                                  ipoint1,dpar, Dvalues, rcollection)
                                    
          if (iand(casmComplexity,not(BCASM_DISCFORDEFMAT)) .ne. 0) then  
            ! Save the computed function value
            DdofValue(ilocalEdge,ielidx) = Dvalues(1) 
          end if
          
          ! A value of SYS_INFINITY indicates a do-nothing node inside of
          ! Dirichlet boundary.
          if (Dvalues(1) .ne. SYS_INFINITY) then
            ! Set the DOF number < 0 to indicate that it is Dirichlet.
            ! ilocalEdge is the number of the local edge - and at the same
            ! time the number of the local DOF of Q1, as an edge always
            ! follows a corner vertex!
            Idofs(ilocalEdge,ielidx) = -abs(Idofs(ilocalEdge,ielidx))
          end if
        end if
        
        ! The right point does not have to be checked! It comes later
        ! with the next edge. The situation when an element crosses the
        ! maximum parameter value with its boundary is handled by the
        ! outer DO-LOOP:
        ! A boundary region with parameter value e.g. [3.0,TMAX]
        ! will produce two index sets: One index set for [0.0, 0.0]
        ! and one for [3.0, TMAX).
        
      case (EL_DG_P1_2D)

        ! Left point inside? -> Corresponding DOF must be computed
        if ( ipoint1 .ne. 0 ) then
        
          ! If parameter values are available, get the parameter value.
          ! Otherwise, take the standard value from above!
          if (associated(p_DvertexParameterValue)) &
            dpar = p_DvertexParameterValue(I)
                
          call fgetBoundaryValues (Icomponents,p_rspatialDiscr,&
                                  rboundaryRegion,ielement, DISCBC_NEEDFUNC,&
                                  ipoint1,dpar, Dvalues, rcollection)
                                    
          if (iand(casmComplexity,not(BCASM_DISCFORDEFMAT)) .ne. 0) then  
            ! Save the computed function value
            DdofValue(ilocalEdge,ielidx) = Dvalues(1) 
          end if
          
          ! A value of SYS_INFINITY indicates a do-nothing node inside of
          ! Dirichlet boundary.
          if (Dvalues(1) .ne. SYS_INFINITY) then
            ! Set the DOF number < 0 to indicate that it is Dirichlet.
            ! ilocalEdge is the number of the local edge - and at the same
            ! time the number of the local DOF of Q1, as an edge always
            ! follows a corner vertex!
            Idofs(ilocalEdge,ielidx) = -abs(Idofs(ilocalEdge,ielidx))
          end if
        end if
        
        ! Right point inside? -> Corresponding DOF must be computed
        if ( ipoint2 .ne. 0 ) then
        
          ! Index of the endpoint. If we are at the last element, the
          ! endpoint is the start point.
          if (ielidx .lt. icount) then
            i2 = IelementsAtBoundaryIdx (ielidx)
          else
            i2 = IelementsAtBoundaryIdx (1)
          end if

          ! If parameter values are available, get the parameter value.
          ! Otherwise, take the standard value from above!
          if (associated(p_DvertexParameterValue)) &
            dpar = p_DvertexParameterValue(i2)
                
          call fgetBoundaryValues (Icomponents,p_rspatialDiscr,&
                                  rboundaryRegion,ielement, DISCBC_NEEDFUNC,&
                                  ipoint1,dpar, Dvalues, rcollection)
                                    
          if (iand(casmComplexity,not(BCASM_DISCFORDEFMAT)) .ne. 0) then  
            ! Save the computed function value
            DdofValue(mod(ilocalEdge,3)+1,ielidx) = Dvalues(1) 
          end if
          
          ! A value of SYS_INFINITY indicates a do-nothing node inside of
          ! Dirichlet boundary.
          if (Dvalues(1) .ne. SYS_INFINITY) then
            ! Set the DOF number < 0 to indicate that it is Dirichlet.
            ! ilocalEdge is the number of the local edge - and at the same
            ! time the number of the local DOF of Q1, as an edge always
            ! follows a corner vertex!
            Idofs(ilocalEdge,ielidx) = -abs(Idofs(ilocalEdge,ielidx))
            Idofs(mod(ilocalEdge,3)+1,ielidx) = -abs(Idofs(mod(ilocalEdge,3)+1,ielidx))
          end if
        end if

      case (EL_P2,EL_Q2)

        ! Left point inside? -> Corresponding DOF must be computed
        if ( ipoint1 .ne. 0 ) then
        
          ! If parameter values are available, get the parameter value.
          ! Otherwise, take the standard value from above!
          if (associated(p_DvertexParameterValue)) &
            dpar = p_DvertexParameterValue(I)
                
          call fgetBoundaryValues (Icomponents,p_rspatialDiscr,&
                                  rboundaryRegion,ielement, DISCBC_NEEDFUNC,&
                                  ipoint1,dpar, Dvalues, rcollection)
                                    
          if (iand(casmComplexity,not(BCASM_DISCFORDEFMAT)) .ne. 0) then   
            ! Save the computed function value of the corner.
            DdofValue(ilocalEdge,ielidx) = Dvalues(1) 
          end if
          
          ! A value of SYS_INFINITY indicates a do-nothing node inside of
          ! Dirichlet boundary.
          if (Dvalues(1) .ne. SYS_INFINITY) then
            ! Set the DOF number < 0 to indicate that it is Dirichlet.
            ! ilocalEdge is the number of the local edge - and at the same
            ! time the number of the local DOF of Q1, as an edge always
            ! follows a corner vertex!
            Idofs(ilocalEdge,ielidx) = -abs(Idofs(ilocalEdge,ielidx))
          end if
        end if
        
        ! The right point does not have to be checked! It comes later
        ! with the next edge. The situation when an element crosses the
        ! maximum parameter value with its boundary is handled by the
        ! outer DO-LOOP:
        ! A boundary region with parameter value e.g. [3.0,TMAX]
        ! will produce two index sets: One index set for [0.0, 0.0]
        ! and one for [3.0, TMAX).
        !
        ! Edge inside? -> Calculate point value on midpoint of edge iedge
        if ( iedge .ne. 0 ) then
          
          ! If parameter values are available, get the parameter value.
          ! Otherwise, take the standard value from above!
          if (associated(p_DedgeParameterValue)) &
            dpar = p_DedgeParameterValue(I)
                
          call fgetBoundaryValues (Icomponents,p_rspatialDiscr,&
                                  rboundaryRegion,ielement, DISCBC_NEEDFUNCMID,&
                                  iedge,dpar, Dvalues, rcollection)
                                    
          if (iand(casmComplexity,not(BCASM_DISCFORDEFMAT)) .ne. 0) then
            ! Save the computed function value of the edge midpoint.
            ! This is found at position ilocalEdge+4, as the first four
            ! elements in DdofValue correspond to the corners!
            DdofValue(ilocalEdge+nve,ielidx) = Dvalues(1) 
          end if
          
          ! A value of SYS_INFINITY indicates a do-nothing node inside of
          ! Dirichlet boundary.
          if (Dvalues(1) .ne. SYS_INFINITY) then
            ! Set the DOF number < 0 to indicate that it is Dirichlet
            ! ilocalEdge is the number of the local edge, corresponding
            ! to the local DOF ilocalEdge+4, as the first four elements
            ! in this array correspond to the values in the corners.
            Idofs(ilocalEdge+nve,ielidx) = -abs(Idofs(ilocalEdge+nve,ielidx))
          end if
              
          ! The element midpoint does not have to be considered, as it cannot
          ! be on the boundary.
        end if

      case (EL_QP1)
        ! Three DOF`s: Function value in the element midpoint 
        ! and derivatives.
        ! Either the edge or an adjacent vertex is on the boundary.
        
        ! If parameter values are available, get the parameter value.
        ! Otherwise, take the standard value from above!
        if (associated(p_DvertexParameterValue)) &
          dpar = p_DvertexParameterValue(I)
              
        ! Get the value at the corner point and accept it as
        ! Dirichlet value.
        call fgetBoundaryValues (Icomponents,p_rspatialDiscr,&
                                rboundaryRegion,ielement, DISCBC_NEEDFUNC,&
                                ipoint1,dpar, Dvalues,rcollection)
                                
        if (iand(casmComplexity,not(BCASM_DISCFORDEFMAT)) .ne. 0) then

          ! Dvalues(1) gives us the function value in the point. Save it
          DdofValue(1,ielidx) = Dvalues(1)
          
          ! Save 0 as X- and Y-derivative.
          DdofValue(2,ielidx) = 0.0_DP
          DdofValue(3,ielidx) = 0.0_DP
          
        end if
        
        ! A value of SYS_INFINITY indicates a do-nothing node inside of
        ! Dirichlet boundary.
        if (Dvalues(1) .ne. SYS_INFINITY) then
          ! Set the DOF numbers < 0 to indicate that it is Dirichlet
          Idofs(1:3,ielidx) = -abs(Idofs(1:3,ielidx))
        end if
        
      case (EL_P1T)

        ! Edge midpoint based element.
        !
        ! Edge inside? -> Calculate point value on midpoint of edge iedge
        if ( iedge .ne. 0 ) then
          ! If parameter values are available, get the parameter value.
          ! Otherwise, take the standard value from above!
          if (associated(p_DedgeParameterValue)) &
              dpar = p_DedgeParameterValue(I)
          
          call fgetBoundaryValues (Icomponents,p_rspatialDiscr,&
              rboundaryRegion,ielement, DISCBC_NEEDFUNCMID,&
              iedge,dpar, Dvalues,rcollection)
          
          if (iand(casmComplexity,not(BCASM_DISCFORDEFMAT)) .ne. 0) then
            ! Save the computed function value
            DdofValue(ilocalEdge,ielidx) = Dvalues(1) 
          end if
          
          ! A value of SYS_INFINITY indicates a do-nothing node inside of
          ! Dirichlet boundary.
          if (Dvalues(1) .ne. SYS_INFINITY) then
            ! Set the DOF number < 0 to indicate that it is Dirichlet
            Idofs(ilocalEdge,ielidx) = -abs(Idofs(ilocalEdge,ielidx))
          end if
        end if

      case (EL_Q1T,EL_Q1TB)
      
        ! The Q1T-element has different variants. Check which variant we have
        ! and choose the right way to calculate boundary values.
      
        if (iand(celement,int(2**16,I32)) .ne. 0) then
        
          ! Integral mean value based element.
          
          ! Edge inside? -> Calculate integral mean value over the edge
          if ( iedge .ne. 0 ) then
          
            if (iand(cmyoptions,BCASM_DISCOPT_INTBYCB) .ne. 0) then
              ! Callback routine calculates the exact integral mean value.
              !
              ! If parameter values are available, get the parameter value.
              ! Otherwise, take the standard value from above!
              if (associated(p_DedgeParameterValue)) &
                dpar = p_DedgeParameterValue(I)

              call fgetBoundaryValues (Icomponents,p_rspatialDiscr,&
                                      rboundaryRegion,ielement, DISCBC_NEEDINTMEAN,&
                                      iedge,dpar, Dvalues,rcollection)

              if (iand(casmComplexity,not(BCASM_DISCFORDEFMAT)) .ne. 0) then
                ! Save the computed function value
                DdofValue(ilocalEdge,ielidx) = Dvalues(1) 
              end if
              
              ! A value of SYS_INFINITY indicates a do-nothing node inside of
              ! Dirichlet boundary.
              if (Dvalues(1) .ne. SYS_INFINITY) then
                ! Set the DOF number < 0 to indicate that it is Dirichlet
                Idofs(ilocalEdge,ielidx) = -abs(Idofs(ilocalEdge,ielidx))
              end if
              
            else
              ! Use a 2-point Gauss cubature formula to compute the integral
              ! mean value.
              !
              ! Ok, that's a little bit more tricky. Get the parameter values
              ! of the points at first.
              !
              ! We neet to set up two values for each edge E: On one hand the
              ! integral mean value as for Q1T (with x:[-1,1]->E):
              !             1/|E| int_[-1,1] v(x(t)) dt
              ! which is called '0th moment', on the other hand the '1st moment',
              ! which is the integral mean value:
              !             1/|E| int_[-1,1] v(x(t)) * t dt
              ! We do this by a 2-point gauss formula by asking the callback routine
              ! for function values in the Gauss points.
              !
              ! Prepare the calculation of the parameter value of the point where
              ! we evaluate the boundary.
              !
              ! 1st Gauss point
              if (associated(p_DvertexParameterValue)) then

                ! Number of vertices on the current boundary component?
                nvbd = p_IboundaryCpIdx(rboundaryRegion%iboundCompIdx+1) - &
                        p_IboundaryCpIdx(rboundaryRegion%iboundCompIdx)
                        
                ! Get the start- and end-parameter value of the vertices on the edge.
                dpar1 = p_DvertexParameterValue(I)
                if (I+1 .lt. p_IboundaryCpIdx(rboundaryRegion%iboundCompIdx+1)) then
                  dpar2 = p_DvertexParameterValue(I+1)
                else
                  dpar2 = boundary_dgetMaxParVal(p_rspatialDiscr%p_rboundary,&
                      rboundaryRegion%iboundCompIdx)
                end if

                ! Calculate the position of the first Gauss point on that edge.
                call mprim_linearRescale(Q2G1,-1.0_DP,1.0_DP,dpar1,dpar2,dpar)
              else
                ! Dummy; hopefully this is never used...
                dpar = Q2G1
              end if
              
              ! Get the value in the Gauss point
              call fgetBoundaryValues (Icomponents,p_rspatialDiscr,&
                                      rboundaryRegion,ielement, DISCBC_NEEDFUNC,&
                                      iedge,dpar, Dvalues,rcollection)
              dval1 = Dvalues(1)
              
              ! 2nd Gauss point
              if (associated(p_DvertexParameterValue)) then

                ! Calculate the position of the 2nd Gauss point on that edge.
                call mprim_linearRescale(Q2G2,-1.0_DP,1.0_DP,dpar1,dpar2,dpar)
              else
                ! Dummy; hopefully this is never used...
                dpar = Q2G2
              end if
              
              ! Get the value in the Gauss point
              call fgetBoundaryValues (Icomponents,p_rspatialDiscr,&
                                      rboundaryRegion,ielement, DISCBC_NEEDFUNC,&
                                      iedge,dpar, Dvalues,rcollection)
              dval2 = Dvalues(1)
              
              ! Calculate the integral.
              !
              ! A value of SYS_INFINITY indicates a do-nothing node inside of
              ! Dirichlet boundary.
              if ((dval1 .ne. SYS_INFINITY) .and. (dval2 .ne. SYS_INFINITY)) then
              
                ! Compute the integral mean value of the 0th moment -- that
                ! is:  1/|E| int_E v dx ~ 1/2 * (v(g1)+v(g2))
                ! This is the same value as for Q1T, but we do the integration
                ! manually by using Gauss.
                dval = ( dval1 + dval2 ) * 0.5_DP
              
                if (iand(casmComplexity,not(BCASM_DISCFORDEFMAT)) .ne. 0) then
                  ! Save the computed function value
                  DdofValue(ilocalEdge,ielidx) = dval
                end if
                
                ! A value of SYS_INFINITY indicates a do-nothing node inside of
                ! Dirichlet boundary.
                if (Dvalues(1) .ne. SYS_INFINITY) then
                  ! Set the DOF number < 0 to indicate that it is Dirichlet
                  Idofs(ilocalEdge,ielidx) = -abs(Idofs(ilocalEdge,ielidx))
                end if
              end if

            end if
            
          end if
                                      
        else
          
          ! Edge midpoint based element.
          !
          ! Edge inside? -> Calculate point value on midpoint of edge iedge
          if ( iedge .ne. 0 ) then
            ! If parameter values are available, get the parameter value.
            ! Otherwise, take the standard value from above!
            if (associated(p_DedgeParameterValue)) &
              dpar = p_DedgeParameterValue(I)

            call fgetBoundaryValues (Icomponents,p_rspatialDiscr,&
                                    rboundaryRegion,ielement, DISCBC_NEEDFUNCMID,&
                                    iedge,dpar, Dvalues,rcollection)
                                      
            if (iand(casmComplexity,not(BCASM_DISCFORDEFMAT)) .ne. 0) then
              ! Save the computed function value
              DdofValue(ilocalEdge,ielidx) = Dvalues(1) 
            end if
            
            ! A value of SYS_INFINITY indicates a do-nothing node inside of
            ! Dirichlet boundary.
            if (Dvalues(1) .ne. SYS_INFINITY) then
              ! Set the DOF number < 0 to indicate that it is Dirichlet
              Idofs(ilocalEdge,ielidx) = -abs(Idofs(ilocalEdge,ielidx))
            end if
          end if

        end if

      case (EL_Q2T,EL_Q2TB)
      
        ! The Q2T-element is only integral mean value based.
        ! On the one hand, we have integral mean values over the edges.
        ! On the other hand, we have integral mean values of function*parameter
        ! value on the edge.
        !
        ! Edge inside? 
        if ( iedge .ne. 0 ) then
          
          ! We neet to set up two values for each edge E: On one hand the
          ! integral mean value as for Q1T (with x:[-1,1]->E):
          !             1/|E| int_[-1,1] v(x(t)) dt
          ! which is called '0th moment', on the other hand the '1st moment',
          ! which is the integral mean value:
          !             1/|E| int_[-1,1] v(x(t)) * t dt
          ! We do this by a 2-point gauss formula by asking the callback routine
          ! for function values in the Gauss points.
          !
          ! Prepare the calculation of the parameter value of the point where
          ! we evaluate the boundary.
          !
          ! 1st Gauss point
          if (associated(p_DvertexParameterValue)) then

            ! Number of vertices on the current boundary component?
            nvbd = p_IboundaryCpIdx(rboundaryRegion%iboundCompIdx+1) - &
                    p_IboundaryCpIdx(rboundaryRegion%iboundCompIdx)
                    
            ! Get the start- and end-parameter value of the vertices on the edge.
            dpar1 = p_DvertexParameterValue(I)
            if (I+1 .lt. p_IboundaryCpIdx(rboundaryRegion%iboundCompIdx+1)) then
              dpar2 = p_DvertexParameterValue(I+1)
            else
              dpar2 = boundary_dgetMaxParVal(p_rspatialDiscr%p_rboundary,&
                  rboundaryRegion%iboundCompIdx)
            end if

            ! Calculate the position of the first Gauss point on that edge.
            call mprim_linearRescale(Q2G1,-1.0_DP,1.0_DP,dpar1,dpar2,dpar)
          else
            ! Dummy; hopefully this is never used...
            dpar = Q2G1
          end if
          
          ! Get the value in the Gauss point
          call fgetBoundaryValues (Icomponents,p_rspatialDiscr,&
                                  rboundaryRegion,ielement, DISCBC_NEEDFUNC,&
                                  iedge,dpar, Dvalues,rcollection)
          dval1 = Dvalues(1)
          
          ! 2nd Gauss point
          if (associated(p_DvertexParameterValue)) then

            ! Calculate the position of the 2nd Gauss point on that edge.
            call mprim_linearRescale(Q2G2,-1.0_DP,1.0_DP,dpar1,dpar2,dpar)
          else
            ! Dummy; hopefully this is never used...
            dpar = Q2G2
          end if
          
          ! Get the value in the Gauss point
          call fgetBoundaryValues (Icomponents,p_rspatialDiscr,&
                                  rboundaryRegion,ielement, DISCBC_NEEDFUNC,&
                                  iedge,dpar, Dvalues,rcollection)
          dval2 = Dvalues(1)
          
          ! A value of SYS_INFINITY indicates a do-nothing node inside of
          ! Dirichlet boundary.
          if ((dval1 .ne. SYS_INFINITY) .and. (dval2 .ne. SYS_INFINITY)) then
            
            if (iand(casmComplexity,not(BCASM_DISCFORDEFMAT)) .ne. 0) then

!                ! Convert the parameter values of the corners of the edge
!                ! into length parametrisation. We need this to calculate the length
!                ! of the edge.
!                dpar1 = boundary_convertParameter(p_rspatialDiscr%p_rboundary,&
!                    rboundaryRegion%iboundCompIdx, dpar1, &
!                    BDR_PAR_01, BDR_PAR_LENGTH)
!
!                dpar2 = boundary_convertParameter(p_rspatialDiscr%p_rboundary,&
!                    rboundaryRegion%iboundCompIdx, dpar2, &
!                    BDR_PAR_01, BDR_PAR_LENGTH)
!                    
!                IF (dpar2 .EQ. 0.0_DP) THEN
!                  ! Wrap around
!                  dpar2 = boundary_dgetMaxParVal(p_rspatialDiscr%p_rboundary,&
!                    rboundaryRegion%iboundCompIdx,BDR_PAR_LENGTH)
!                END IF

              ! Compute the integral mean value of the 0th moment -- that
              ! is:  1/|E| int_E v dx ~ 1/2 * (v(g1)+v(g2))
              ! This is the same value as for Q1T, but we do the integration
              ! manually by using Gauss.
              dval = ( dval1 + dval2 ) * 0.5_DP
              
              ! Save the computed value
              DdofValue(ilocalEdge,ielidx) = dval
              
              ! Calculate the integral mean value of the 1st moment:
              !   1/|E| int_E v(x(t))*t dx ~ 1/2 * (v(x(g1))*g1+v(x(g2))*g2)
              dval = ( (dval1*Q2G1) + (dval2*Q2G2) ) * 0.5_DP
            
              ! Save the computed function value
              DdofValue(ilocalEdge+nve,ielidx) = dval 
            end if
            
            ! Set the DOF number < 0 to indicate that it is Dirichlet
            Idofs(ilocalEdge,ielidx) = -abs(Idofs(ilocalEdge,ielidx))

            ! Set the DOF number < 0 to indicate that it is Dirichlet
            Idofs(ilocalEdge+nve,ielidx) = -abs(Idofs(ilocalEdge+nve,ielidx))
          end if
          
        end if

      case DEFAULT
      
        print *,'bcasm_newDirichletBC: Unsupported element!'
        call sys_halt()
      
      end select
      
    end do

    ! Temp arrays no more necessary    
    deallocate(IverticesAtBoundaryIdx)
    deallocate(IedgesAtBoundaryIdx)
    deallocate(IelementsAtBoundary)
    deallocate(IelementsAtBoundaryIdx)
    
    ! Now count how many values we actually have.
    icount = 0
    do J=1,size(Idofs,2)
      do I=1,size(Idofs,1)
        if (Idofs(I,J) < 0) icount = icount + 1
      end do
    end do
    
    p_rdirichletBCs%nDOF = icount
    
    if (icount .gt. 0) then
    
      ! Allocate arrays for storing these DOF`s and their values - if values are
      ! computed.
      call storage_new('bcasm_newDirichletBConRealBd', 'h_IdirichletDOFs', &
                      icount, ST_INT, p_rdirichletBCs%h_IdirichletDOFs, &
                      ST_NEWBLOCK_NOINIT)
      call storage_getbase_int(p_rdirichletBCs%h_IdirichletDOFs,p_IdirichletDOFs)
      
      if (iand(casmComplexity,int(not(BCASM_DISCFORDEFMAT),I32)) .ne. 0) then
        call storage_new('bcasm_newDirichletBConRealBd', 'h_DdirichletValues', & 
                        icount, ST_DOUBLE, p_rdirichletBCs%h_DdirichletValues, &
                        ST_NEWBLOCK_NOINIT)
        call storage_getbase_double(p_rdirichletBCs%h_DdirichletValues,p_DdirichletValues)
      end if
      
      ! Transfer the DOF`s and their values to these arrays.
      icount = 0
      do J=1,size(Idofs,2)
        do I=1,size(Idofs,1)
          if (Idofs(I,J) < 0) then
            icount = icount + 1
            p_IdirichletDOFs(icount) = abs(Idofs(I,J))
            if (iand(casmComplexity,not(BCASM_DISCFORDEFMAT)) .ne. 0) then
              p_DdirichletValues(icount) = DdofValue(I,J)
            end if
          end if
        end do
      end do
    
    else
    
      ! Let us hope there is nothing saved here :-)
      p_rdirichletBCs%h_IdirichletDOFs = ST_NOHANDLE
      p_rdirichletBCs%h_DdirichletValues = ST_NOHANDLE
    
    end if
    
    ! Remove temporary memory, finish.
    if (iand(casmComplexity,not(BCASM_DISCFORDEFMAT)) .ne. 0) then
      deallocate (DdofValue)
    end if
    deallocate (Idofs)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine bcasm_releaseDirichlet (rdiscreteBCDirichlet)
  
!<description>
  ! This routine cleans up the discrete Dirichlet boundary conditions
  ! rdiscreteBCDirichlet.
!</description>

!<inputoutput>
  ! The discrete-BC structure which is to be cleaned up
  type(t_discreteBCDirichlet), intent(inout) :: rdiscreteBCDirichlet
!</inputoutput>

!</subroutine>

  ! Release what is associated
  
  rdiscreteBCDirichlet%nDOF = 0
  if (rdiscreteBCDirichlet%h_DdirichletValues .ne. ST_NOHANDLE) &
    call storage_free(rdiscreteBCDirichlet%h_DdirichletValues)
  if (rdiscreteBCDirichlet%h_IdirichletDOFs .ne. ST_NOHANDLE) &
    call storage_free(rdiscreteBCDirichlet%h_IdirichletDOFs)

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine bcasm_discrBCFeastMirror (rblockDiscretisation, &
                                       iequation, rboundaryRegion, rdiscreteBC, &
                                       rconfigFeastMirrorBC, ccomplexity)
  
!<description>
  ! Creates a discrete version of FEAST mirror boundary conditions.
  ! rboundaryRegion describes the region which is to be discretised. The discretised
  ! boundary conditions are added to rdiscreteBC.
!</description>

!<input>
  ! The discretisation structure of the underlying discretisation. The boundary
  ! conditions inside of this structure are discretised.
  type(t_blockDiscretisation), intent(in), target :: rblockDiscretisation

  ! An identifier for the equation, this boundary condition refers to.
  ! >= 1. 1=first equation (e.g. X-velocity), 2=2nd equation (e.g. 
  ! Y-velocity), etc.
  integer, intent(in) :: iequation

  ! A boundary-condition-region object, describing the position on the
  ! boundary where boundary conditions should be imposed.
  type(t_boundaryRegion), intent(in) :: rboundaryRegion

  ! Configuration block for the FEAST mirror boundary conditions.
  type(t_configDiscreteFeastMirrorBC), intent(in) :: rconfigFeastMirrorBC

  ! Optional: A combination of BCASM_DISCFORxxx constants that specify
  ! the complexity of the discretisation that is to perform. This allows to
  ! discretise only parts of the BC`s, e.g. only setting up those
  ! information that are necessary for filtering defect vectors.
  ! If not specified, BCASM_DISCFORALL is assumed, i.e. the resulting
  ! boundary conditions can be used for everything.
  integer(I32), intent(in), optional :: ccomplexity
!</input>  

!<inputoutput>
  ! A t_discreteBC structures, representing the boundary discretised 
  ! in a discretisation-dependent way. The new BC`s are added to this structure.
  type(t_discreteBC), intent(inout) :: rdiscreteBC
!</inputoutput>

!</subroutine>

    ! local variables
    integer(I32) :: celement
    integer :: i,icount
    type(t_discreteBCFeastMirror),pointer       :: p_rfeastMirrorBCs
    type(t_triangulation), pointer              :: p_rtriangulation
    type(t_spatialDiscretisation), pointer      :: p_rspatialDiscr
    integer, dimension(:,:), pointer            :: p_IverticesAtElement
    
    integer, dimension(:), pointer              :: p_ImirrorDOFs
    
    integer, dimension(:), allocatable          :: IverticesAtBoundaryIdx
    integer, dimension(:), allocatable          :: IelementsAtBoundary
    
    type (t_boundaryregion) :: rboundaryRegionClosed

    integer ::iidx
    type(t_discreteBCEntry), pointer :: p_rdiscreteBCentry
    
    integer(I32) :: casmComplexity
    
    casmComplexity = BCASM_DISCFORALL
    if (present(ccomplexity)) casmComplexity = ccomplexity

    ! NOTE:
    ! The routine is definitively buggy!
    ! It was build to find a bug in FEAST. E.g. it works only if there are
    ! 'enough' points in the boundary region and if the discretisation is
    ! done with Q1!
    
    ! The BC`s only exist as modification of the matrix.
    ! If we should not compute them for the matrix/defect, we do not 
    ! have to do anything.
    if (iand(casmComplexity,BCASM_DISCFORMAT) .eq. 0) return

    ! Get the discretisation structures from one of the components of the solution
    ! vector that is to be modified.
    p_rspatialDiscr => rblockDiscretisation%RspatialDiscr(iequation)

    if (p_rspatialDiscr%ccomplexity .ne. SPDISC_UNIFORM) then
      print *,'Discrete FEAST mirror boundary conditions currently only supported'
      print *,'for uniform discretisations!'
      call sys_halt()
    end if
    
    ! For easier access:
    p_rtriangulation => p_rspatialDiscr%p_rtriangulation
    call storage_getbase_int2d(p_rtriangulation%h_IverticesAtElement,p_IverticesAtElement)

    if (p_rtriangulation%ndim .ne. NDIM2D) then
      print *,'FEAST mirror boundary only support 2D!'
      call sys_halt()
    end if
    
    celement = p_rspatialDiscr%RelementDistr(1)%celement
    if (elem_getPrimaryElement(celement) .ne. EL_Q1) then
      print *,'Discrete FEAST mirror boundary conditions currently only supported'
      print *,'for Q1 element!'
      call sys_halt()
    end if

    ! Get a new BC entry
    call bcasm_newBCentry(rdiscreteBC, iidx)
    p_rdiscreteBCentry => rdiscreteBC%p_RdiscBCList(iidx)

    ! We have FEAST mirror boundary conditions
    p_rdiscreteBCentry%itype = DISCBC_TPFEASTMIRROR
    
    ! Fill the structure for discrete Dirichlet BC`s in the
    ! t_discreteBCEntry structure
    p_rfeastMirrorBCs => p_rdiscreteBCentry%rfeastMirrorBCs
    
    p_rfeastMirrorBCs%icomponent = iequation
    
    ! Copy the coarsening level from the configuration block
    p_rfeastMirrorBCs%icoarseningLevel = rconfigFeastMirrorBC%icoarseningLevel
    p_rfeastMirrorBCs%isubtype = rconfigFeastMirrorBC%isubtype
    
    ! We have to deal with all DOF`s on the boundary. This is highly element
    ! dependent and therefore a little bit tricky :(
    !
    ! As we are in 2D, we can use parameter values at first to figure out,
    ! which points are on the boundary.
    allocate(IverticesAtBoundaryIdx(p_rtriangulation%NVBD))
    allocate(IelementsAtBoundary(p_rtriangulation%NVBD))
    
    call bcasm_getElementsInBdRegion (p_rtriangulation,rboundaryRegion, &
        icount, IelementsAtBoundary, IvtLocal=IverticesAtBoundaryIdx)
    
    ! Cancel if the set is empty!
    if (icount .eq. 0) then
      ! Assign ST_NOHANDLE to the h_ImirrorBCs handle to instruct the BC filter
      ! to do nothing.
      p_rfeastMirrorBCs%h_ImirrorDOFs = ST_NOHANDLE
      deallocate(IverticesAtBoundaryIdx)
      deallocate(IelementsAtBoundary)
      return
    end if
    
    ! Allocate an array for all the DOF`s
    call storage_new('bcasm_discrBCFeastMirror', 'h_ImirrorDOFs', &
                    icount, ST_INT, &
                    p_rfeastMirrorBCs%h_ImirrorDOFs, ST_NEWBLOCK_ZERO)
    call storage_getbase_int(p_rfeastMirrorBCs%h_ImirrorDOFs,p_ImirrorDOFs)
    
    ! Put all DOF`s in that array.
    do i=1,icount
      p_ImirrorDOFs(i) = p_IverticesAtElement (&
          IverticesAtBoundaryIdx(i),IelementsAtBoundary(i))
    end do
    
    ! Sort the array for quicker access.
    call sort_int(p_ImirrorDOFs)
    
    ! p_ImirrorDOFs contains all BC`s that have to be processed.
    ! But it does not contain all DOF`s that are to be doubled in the matrix.
    !
    ! Duplicate the boundary region and include start- and endpoint
    rboundaryRegionClosed = rboundaryRegion
    rboundaryRegionClosed%iproperties = BDR_PROP_WITHSTART+BDR_PROP_WITHEND
    
    ! Then collect all DOF`s inside of that like above.
    call bcasm_getElementsInBdRegion (p_rtriangulation,rboundaryRegionClosed, &
        icount, IelementsAtBoundary, IvtLocal=IverticesAtBoundaryIdx)
                                   
    ! Allocate an array for all the DOF`s
    call storage_new('bcasm_discrBCFeastMirror', 'h_ImirrorDOFsClosed', &
                    icount, ST_INT, &
                    p_rfeastMirrorBCs%h_ImirrorDOFsClosed, ST_NEWBLOCK_ZERO)
    call storage_getbase_int(p_rfeastMirrorBCs%h_ImirrorDOFsClosed,p_ImirrorDOFs)
    
    ! Put all DOF`s in that array.
    do i=1,icount
      p_ImirrorDOFs(i) = p_IverticesAtElement (&
          IverticesAtBoundaryIdx(i),IelementsAtBoundary(i))
    end do
    
    ! Sort the array for quicker access.
    call sort_int(p_ImirrorDOFs)
    
    ! Clean up, finish
    deallocate(IverticesAtBoundaryIdx)
    deallocate(IelementsAtBoundary)
                                     
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine bcasm_releaseFeastMirror (rdiscreteBCFeastMirror)
  
!<description>
  ! This routine cleans up the discrete FEAST mirror boundary conditions
  ! rdiscreteBCFeastMirror.
!</description>

!<inputoutput>
  ! The discrete-BC structure which is to be cleaned up
  type(t_discreteBCFeastMirror), intent(inout) :: rdiscreteBCFeastMirror
!</inputoutput>

!</subroutine>

  ! Release what is associated
  
  if (rdiscreteBCFeastMirror%h_ImirrorDOFsClosed .ne. ST_NOHANDLE) &
    call storage_free(rdiscreteBCFeastMirror%h_ImirrorDOFsClosed)
  if (rdiscreteBCFeastMirror%h_ImirrorDOFs .ne. ST_NOHANDLE) &
    call storage_free(rdiscreteBCFeastMirror%h_ImirrorDOFs)

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine bcasm_newPdropBConRealBd (rblockDiscretisation, &
                                     Iequations, rboundaryRegion, rdiscreteBC, &
                                     fgetBoundaryValues,rcollection, ccomplexity)
  
!<description>
  ! Creates a discrete version of pressure drop boundary conditions.
  ! rboundaryRegion describes the region which is to be discretised. The discretised
  ! boundary conditions are added to rdiscreteBC.
!</description>

!<input>
  ! The discretisation structure of the underlying discretisation. The boundary
  ! conditions inside of this structure are discretised.
  type(t_blockDiscretisation), intent(in), target :: rblockDiscretisation

  ! A list of identifiers for the velocity equations, this boundary condition
  ! modifies. Usually (1,2) for X- and Y-velocity.
  integer, dimension(:), intent(in) :: Iequations

  ! A boundary-condition-region object, describing the position on the
  ! boundary where boundary conditions should be imposed.
  type(t_boundaryRegion), intent(in) :: rboundaryRegion

  ! A callback function that calculates values on the boundary.
  ! Is declared in the interface include file 'intf_bcassembly.inc'.
  include 'intf_bcassembly.inc'
  
  ! Optional: A collection structure to inform the callback function with
  ! additional information. 
  type(t_collection), intent(inout), optional :: rcollection

  ! Optional: A combination of BCASM_DISCFORxxx constants that specify
  ! the complexity of the discretisation that is to perform. This allows to
  ! discretise only parts of the BC`s, e.g. only setting up those
  ! information that are necessary for filtering defect vectors.
  ! If not specified, BCASM_DISCFORALL is assumed, i.e. the resulting
  ! boundary conditions can be used for everything.
  integer(I32), intent(in), optional :: ccomplexity
!</input>  

!<inputoutput>
  ! A t_discreteBC structures, representing the boundary discretised 
  ! in a discretisation-dependent way. The new BC`s are added to this structure.
  type(t_discreteBC), intent(inout) :: rdiscreteBC
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,icount
    integer(I32) :: celement
    real(DP), dimension(EL_MAXNDER)            :: Dvalues
    real(DP),dimension(NDIM2D)                 :: Dtangential,Dnormal
    integer                                    :: NVT,ipoint1,ipoint2
    integer                                    :: ielement
    integer                                    :: iedge
    integer, dimension(2)                      :: ImodifierSize
    
    type(t_spatialDiscretisation), pointer      :: p_rspatialDiscr
    type(t_triangulation), pointer              :: p_rtriangulation
    integer, dimension(:,:), pointer            :: p_IedgesAtElement
    integer, dimension(:,:), pointer            :: p_IverticesAtElement
    real(DP), dimension(:), pointer             :: p_DedgeParameterValue
    real(DP), dimension(:,:), pointer           :: p_DvertexCoords
    
    type(t_discreteBCpressureDrop), pointer     :: p_rpressureDropBCs
    integer, dimension(:), pointer              :: p_IpressureDropDOFs
    real(DP), dimension(:,:), pointer           :: p_Dmodifier

    integer, dimension(:), allocatable          :: IedgesAtBoundaryIdx
    integer, dimension(:), allocatable          :: IelementsAtBoundary

    integer ::iidx
    type(t_discreteBCEntry), pointer :: p_rdiscreteBCentry

    integer(I32) :: casmComplexity
    
    casmComplexity = BCASM_DISCFORALL
    if (present(ccomplexity)) casmComplexity = ccomplexity

    ! Pressure drop BC`s only exist as modification of the RHS.
    ! If we should not compute them for the RHS, we do not have to do anything.
    if (iand(casmComplexity,BCASM_DISCFORRHS) .eq. 0) return

    ! Get a new BC entry
    call bcasm_newBCentry(rdiscreteBC, iidx)
    p_rdiscreteBCentry => rdiscreteBC%p_RdiscBCList(iidx)

    ! Fill the structure for discrete pressure drop BC`s in the
    ! t_discreteBCEntry structure
    p_rpressureDropBCs => p_rdiscreteBCentry%rpressureDropBCs
    
    ! Get the discretisation structures from one of the components of the solution
    ! vector that is to be modified.
    p_rspatialDiscr => &
      rblockDiscretisation%RspatialDiscr(Iequations(1))

    if (p_rspatialDiscr%ccomplexity .ne. SPDISC_UNIFORM) then
      print *,'Discrete pressure drop boundary conditions currently only supported'
      print *,'for uniform discretisations!'
      call sys_halt()
    end if

    celement = p_rspatialDiscr%RelementDistr(1)%celement
    if (elem_getPrimaryElement(celement) .ne. EL_Q1T) then
      print *,'Discrete pressure drop boundary conditions currently only supported'
      print *,'for Q1~ element!'
      call sys_halt()
    end if
    
    if (p_rspatialDiscr%ndimension .ne. NDIM2D) then
      print *,'Pressure drop boundary conditions only support 2D!'
      call sys_halt()
    end if
    
    ! For easier access:
    p_rtriangulation => p_rspatialDiscr%p_rtriangulation

    ! Note: All elements are of the same type celement.
    !
    ! We have pressure drop boundary conditions
    p_rdiscreteBCentry%itype = DISCBC_TPPRESSUREDROP
    
    ! Which components of the solution vector are affected by this boundary
    ! condition?
    p_rpressureDropBCs%ncomponents = NDIM2D
    allocate(p_rpressureDropBCs%Icomponents(1:NDIM2D))
    p_rpressureDropBCs%Icomponents(1:NDIM2D) = Iequations(1:NDIM2D)
    
    ! We have to deal with all DOF`s on the boundary. This is highly element
    ! dependent and therefore a little bit tricky :(
    ! But here we restrict to Q1~ only, which makes life a little bit easier.
    allocate(IedgesAtBoundaryIdx(p_rtriangulation%NVBD))
    allocate(IelementsAtBoundary(p_rtriangulation%NVBD))
    
    call bcasm_getElementsInBdRegion (p_rtriangulation,rboundaryRegion, &
        icount, IelementsAtBoundary, IedgeLocal=IedgesAtBoundaryIdx)
                                   
    ! Cancel if the set is empty!
    if (icount .eq. 0) then
      deallocate(IedgesAtBoundaryIdx)
      deallocate(IelementsAtBoundary)
      return
    end if

    ! Total number of edges?
    p_rpressureDropBCs%nDOF = icount

    ! Allocate memory to save the DOF`s as well as all modifiers.
    call storage_new('bcasm_discrBCpressureDrop', 'h_IpressureDropDOFs', &
                    icount, ST_INT, p_rpressureDropBCs%h_IpressureDropDOFs, &
                    ST_NEWBLOCK_NOINIT)
    ImodifierSize = (/NDIM2D,icount/)
    call storage_new('bcasm_discrBCpressureDrop', 'h_Dmodifier', & 
                      ImodifierSize, ST_DOUBLE, p_rpressureDropBCs%h_Dmodifier, &
                      ST_NEWBLOCK_NOINIT)
                      
    call storage_getbase_int(p_rpressureDropBCs%h_IpressureDropDOFs,p_IpressureDropDOFs)
    call storage_getbase_double2d(p_rpressureDropBCs%h_Dmodifier,p_Dmodifier)

    ! For easier access:
    call storage_getbase_double(p_rtriangulation%h_DedgeParameterValue,p_DedgeParameterValue)
    call storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
    call storage_getbase_double2D(p_rtriangulation%h_DvertexCoords,p_DvertexCoords)
    call storage_getbase_int2D(p_rtriangulation%h_IedgesAtElement,p_IedgesAtElement)
    NVT = p_rtriangulation%NVT

    ! Now calculate the pressure drop integral; cf. p. 257 (235) in Turek`s book:
    !
    ! The pressure drop boundary condition has to implement
    !
    !     - sum P_j  int_Sj  phi * n  ds
    !
    ! into the RHS vector. For each (velocity) DOF of the boundary,
    ! we save "P_j  int_Sj  phi_k * n  ds" as modifier for the DOF of
    ! the RHS!
    !
    ! Loop through all edges on the boundary belonging to our current 
    ! boundary segment.
    do i=1,icount
      ! Get information about the edge and the adjacent element
      ielement = IelementsAtBoundary(i)
      iedge = p_IedgesAtElement (IedgesAtBoundaryIdx(i),ielement)
      
      ! Get the adjacent points in counterclockwise order. The local
      ! number of the edge and the vertex coincide...
      ipoint1 = p_IverticesAtElement (IedgesAtBoundaryIdx(i),ielement)
      ipoint2 = p_IverticesAtElement (mod(IedgesAtBoundaryIdx(i),TRIA_NVEQUAD2D)+1,ielement)
    
      ! Get the coordinates of the endpoints to build the tangential
      ! vector of the edge:
      Dtangential(1:NDIM2D) = p_DvertexCoords(1:NDIM2D,ipoint2) &
                            - p_DvertexCoords(1:NDIM2D,ipoint1)
                            
      ! Get the inner normal vector. This compensates the '-' sign in front of
      ! the RHS in the formula on poage 269 in Turek`s book where the outer
      ! normal vector is used.
      Dnormal(1) = -Dtangential(2)
      Dnormal(2) =  Dtangential(1)
      
      ! Do not scale the normal vector! The scaling factor cancels out later
      ! when calculating the integral on the edge with the midpoint rule!
      !
      ! At first, we ask the boundary-value routine to give us the
      ! weight P_j in the cubature point (for Q1~ in the edge midpoint)
      ! of the current boundary segment j we are discretising here!
      !
      ! Calculate normal stress in the current midpoint of the edge.
      ! In a later implementation, we might calculate the value in a cubature
      ! point to calculate the current integral - but here we use the midpoint
      ! rule...
      call fgetBoundaryValues (p_rpressureDropBCs%Icomponents,p_rspatialDiscr,&
                                rboundaryRegion,ielement, DISCBC_NEEDNORMALSTRESS,&
                                iedge,p_DedgeParameterValue(I), Dvalues,rcollection)
      
      ! Save the current modifier (normal stress value)*n to the structure for
      ! later implementation in to the RHS vector.
      p_IpressureDropDOFs(i) = iedge
      p_Dmodifier(1:NDIM2D,i) = Dvalues(1:NDIM2D)*Dnormal(1:NDIM2D)
    
    end do ! i
    
    deallocate(IedgesAtBoundaryIdx)
    deallocate(IelementsAtBoundary)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine bcasm_releasePressureDrop (p_rdiscreteBCentryPD)
  
!<description>
  ! This routine cleans up the discrete Pressure Drop boundary conditions
  ! p_rdiscreteBCentryPD.
!</description>

!<inputoutput>
  ! The discrete-BC structure which is to be cleaned up
  type(t_discreteBCpressureDrop), intent(inout) :: p_rdiscreteBCentryPD
!</inputoutput>

!</subroutine>

    ! Release what is associated
    
    p_rdiscreteBCentryPD%nDOF = 0
    if (p_rdiscreteBCentryPD%h_IpressureDropDOFs .ne. ST_NOHANDLE) &
      call storage_free(p_rdiscreteBCentryPD%h_IpressureDropDOFs)
    if (p_rdiscreteBCentryPD%h_Dmodifier .ne. ST_NOHANDLE) &
      call storage_free(p_rdiscreteBCentryPD%h_Dmodifier)
      
    deallocate(p_rdiscreteBCentryPD%Icomponents)
    p_rdiscreteBCentryPD%ncomponents = 0

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine bcasm_newSlipBConRealBd (rblockDiscretisation, &
     Iequations, rboundaryRegion, rdiscreteBC, ccomplexity)
  
!<description>
  ! Creates a discrete version of Slip boundary conditions.
  ! rboundaryRegion describes the region which is to be discretised. The discretised
  ! boundary conditions are added to p_rdiscreteBC.
!</description>

!<input>
  ! The discretisation structure of the underlying discretisation. The boundary
  ! conditions inside of this structure are discretised.
  type(t_blockDiscretisation), intent(in), target :: rblockDiscretisation

  ! A list of identifiers for the velocity equations, this boundary condition
  ! modifies. Usually (1,2) for X- and Y-velocity.
  integer, dimension(:), intent(in) :: Iequations

  ! A boundary-condition-region object, describing the position on the
  ! boundary where boundary conditions should be imposed.
  type(t_boundaryRegion), intent(in) :: rboundaryRegion

  ! Optional: A combination of BCASM_DISCFORxxx constants that specify
  ! the complexity of the discretisation that is to perform. This allows to
  ! discretise only parts of the BC`s, e.g. only setting up those
  ! information that are necessary for filtering defect vectors.
  ! If not specified, BCASM_DISCFORALL is assumed, i.e. the resulting
  ! boundary conditions can be used for everything.
  integer(I32), intent(in), optional :: ccomplexity
!</input>  

!<inputoutput>
  ! A t_discreteBC structures, representing the boundary discretised 
  ! in a discretisation-dependent way. The new BC`s are added to this structure.
  type(t_discreteBC), intent(inout) :: rdiscreteBC
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,icount
    integer(I32) :: celement
    real(DP),dimension(NDIM2D)                  :: Dtangential,Dnormal
    integer                                     :: NVT,ipoint1,ipoint2
    integer                                     :: ielement
    integer                                     :: iedge
    integer, dimension(2)                       :: InormalsSize
    
    type(t_spatialDiscretisation), pointer      :: p_rspatialDiscr
    type(t_triangulation), pointer              :: p_rtriangulation
    integer, dimension(:,:), pointer            :: p_IedgesAtElement
    integer, dimension(:,:), pointer            :: p_IverticesAtElement
    real(DP), dimension(:), pointer             :: p_DedgeParameterValue
    real(DP), dimension(:,:), pointer           :: p_DvertexCoords
    
    type(t_discreteBCSlip), pointer             :: p_rslipBCs
    integer, dimension(:), pointer              :: p_IslipDOFs
    real(DP), dimension(:,:), pointer           :: p_Dnormals
    real(DP) :: d

    integer, dimension(:), allocatable          :: IedgesAtBoundaryIdx
    integer, dimension(:), allocatable          :: IelementsAtBoundary

    integer ::iidx
    type(t_discreteBCEntry), pointer :: p_rdiscreteBCentry

    integer(I32) :: casmComplexity
    
    casmComplexity = BCASM_DISCFORALL
    if (present(ccomplexity)) casmComplexity = ccomplexity

    ! Pressure drop BC`s only exist as modification of the defect vector
    ! and matrix.
    ! If we should not compute them for the matrix/defect, we do not 
    ! have to do anything.
    if (iand(casmComplexity,BCASM_DISCFORDEFMAT) .eq. 0) return

    ! Get a new BC entry
    call bcasm_newBCentry(rdiscreteBC, iidx)
    p_rdiscreteBCentry => rdiscreteBC%p_RdiscBCList(iidx)

    ! Fill the structure for discrete pressure drop BC`s in the
    ! t_discreteBCEntry structure
    p_rslipBCs => p_rdiscreteBCentry%rslipBCs
    
    ! Get the discretisation structures from one of the components of the solution
    ! vector that is to be modified.
    p_rspatialDiscr => &
      rblockDiscretisation%RspatialDiscr(Iequations(1))  

    if (p_rspatialDiscr%ccomplexity .ne. SPDISC_UNIFORM) then
      print *,'Discrete Slip boundary conditions currently only supported'
      print *,'for uniform discretisations!'
      call sys_halt()
    end if

    celement = p_rspatialDiscr%RelementDistr(1)%celement
    if (elem_getPrimaryElement(celement) .ne. EL_Q1T) then
      print *,'Discrete Slip boundary conditions currently only supported'
      print *,'for Q1~ element!'
      call sys_halt()
    end if
    
    if (p_rspatialDiscr%ndimension .ne. NDIM2D) then
      print *,'Pressure drop boundary conditions only support 2D!'
      call sys_halt()
    end if
    
    ! For easier access:
    p_rtriangulation => p_rspatialDiscr%p_rtriangulation

    ! Note: All elements are of the same type celement.
    !
    ! We have pressure drop boundary conditions
    p_rdiscreteBCentry%itype = DISCBC_TPSLIP
    
    ! Which components of the solution vector are affected by this boundary
    ! condition?
    p_rslipBCs%ncomponents = NDIM2D
    allocate(p_rslipBCs%Icomponents(1:NDIM2D))
    p_rslipBCs%Icomponents(1:NDIM2D) = Iequations(1:NDIM2D)
    
    ! We have to deal with all DOF`s on the boundary. This is highly element
    ! dependent and therefore a little bit tricky :(
    ! But here we restrict to Q1~ only, which makes life a little bit easier.
    ! We have to deal with all DOF`s on the boundary. This is highly element
    ! dependent and therefore a little bit tricky :(
    ! But here we restrict to Q1~ only, which makes life a little bit easier.
    allocate(IedgesAtBoundaryIdx(p_rtriangulation%NVBD))
    allocate(IelementsAtBoundary(p_rtriangulation%NVBD))
    
    call bcasm_getElementsInBdRegion (p_rtriangulation,rboundaryRegion, &
        icount, IelementsAtBoundary, IedgeLocal=IedgesAtBoundaryIdx)
                                   
    ! Cancel if the set is empty!
    if (icount .eq. 0) then
      deallocate(IedgesAtBoundaryIdx)
      deallocate(IelementsAtBoundary)
      return
    end if

    ! Total number of edges?
    p_rslipBCs%nDOF = icount

    ! Allocate memory to save the DOF`s as well as all modifiers.
    call storage_new('bcasm_discrBCSlip', 'h_IpressureDropDOFs', &
                    icount, ST_INT, p_rslipBCs%h_IslipDOFs, &
                    ST_NEWBLOCK_NOINIT)
    InormalsSize = (/NDIM2D,icount/)
    call storage_new('bcasm_discrBCSlip', 'h_Dnormals', & 
                      InormalsSize, ST_DOUBLE, p_rslipBCs%h_DnormalVectors, &
                      ST_NEWBLOCK_NOINIT)
                      
    call storage_getbase_int(p_rslipBCs%h_IslipDOFs,p_IslipDOFs)
    call storage_getbase_double2d(p_rslipBCs%h_DnormalVectors,p_Dnormals)

    ! For easier access:
    call storage_getbase_double(p_rtriangulation%h_DedgeParameterValue,p_DedgeParameterValue)
    call storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
    call storage_getbase_double2D(p_rtriangulation%h_DvertexCoords,p_DvertexCoords)
    call storage_getbase_int2D(p_rtriangulation%h_IedgesAtElement,p_IedgesAtElement)
    NVT = p_rtriangulation%NVT

    ! Now calculate the pressure drop integral; cf. p. 257 (235) in Turek`s book:
    !
    ! The pressure drop boundary condition has to implement
    !
    !     - sum P_j  int_Sj  phi * n  ds
    !
    ! into the RHS vector. For each (velocity) DOF of the boundary,
    ! we save "P_j  int_Sj  phi_k * n  ds" as modifier for the DOF of
    ! the RHS!
    !
    ! Loop through all edges on the boundary belonging to our current 
    ! boundary segment.
    do i=1,icount
    
      ! Get information about the edge and the adjacent element
      ielement = IelementsAtBoundary(i)
      iedge = p_IedgesAtElement (IedgesAtBoundaryIdx(i),ielement)
      
      ! Get the adjacent points in counterclockwise order. The local
      ! number of the edge and the vertex coincide...
      ipoint1 = p_IverticesAtElement (IedgesAtBoundaryIdx(i),ielement)
      ipoint2 = p_IverticesAtElement (mod(IedgesAtBoundaryIdx(i),TRIA_NVEQUAD2D)+1,ielement)

      ! Get the coordinates of the endpoints to build the tangential
      ! vector of the edge:
      Dtangential(1:NDIM2D) = p_DvertexCoords(1:NDIM2D,ipoint2) &
                            - p_DvertexCoords(1:NDIM2D,ipoint1)
                            
      ! Get the outer normal vector.
      Dnormal(1) =  Dtangential(2)
      Dnormal(2) = -Dtangential(1)
      
      ! Scale the vector to be of length 1.
      
      d = 1.0_DP / sqrt(Dnormal(1)**2+Dnormal(2)**2)
      Dnormal(1:2) = Dnormal(1:2) * d
      
      ! Save the DOF and the normal of the edge.            
      p_IslipDOFs(i) = iedge
      p_Dnormals(1:NDIM2D,i) = Dnormal(1:NDIM2D)
    
    end do ! i
    
    deallocate(IedgesAtBoundaryIdx)
    deallocate(IelementsAtBoundary)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine bcasm_releaseSlip (rdiscreteBCSlip)
  
!<description>
  ! This routine cleans up the discrete Pressure Drop boundary conditions
  ! rdiscreteBCPD.
!</description>

!<inputoutput>
  ! The discrete-BC structure which is to be cleaned up
  type(t_discreteBCSlip), intent(inout) :: rdiscreteBCSlip
!</inputoutput>

!</subroutine>

    ! Release what is associated
    
    rdiscreteBCSlip%nDOF = 0
    if (rdiscreteBCSlip%h_IslipDOFs .ne. ST_NOHANDLE) &
      call storage_free(rdiscreteBCSlip%h_IslipDOFs)
    if (rdiscreteBCSlip%h_DnormalVectors .ne. ST_NOHANDLE) &
      call storage_free(rdiscreteBCSlip%h_DnormalVectors)

    deallocate(rdiscreteBCSlip%Icomponents)
    rdiscreteBCSlip%ncomponents = 0

  end subroutine

! *****************************************************************************
! Support for boundary conditions on fictitious boundary components
! *****************************************************************************

!<subroutine>

  subroutine bcasm_newDirichletBConFBD (rblockDiscretisation, &
      Iequations, rdiscreteFBC, &
      fgetBoundaryValuesFBC,rcollection, ccomplexity)
  
!<description>
  ! Creates a discrete version of Dirichlet boundary conditions for a fictitious
  ! boundary component.
  ! rboundaryRegion describes the region which is to be discretised. The discretised
  ! boundary conditions are created in rdiscreteFBC, which is assumed
  ! to be undefined when entering this routine.
!</description>

!<input>
  ! The discretisation structure of the underlying discretisation. The boundary
  ! conditions inside of this structure are discretised.
  type(t_blockDiscretisation), intent(in), target :: rblockDiscretisation

  ! An array of identifiers for the equations, this boundary condition 
  ! refers to. Example: Iequations = [1 2] for X-velocity-component (1) and
  ! Y-velocity component (2).
  integer, dimension(:), intent(in) :: Iequations

  ! A callback function that calculates values in the domain.
  ! Is declared in the interface include file 'intf_fbcassembly.inc'.
  include 'intf_fbcassembly.inc'
  
  ! Optional: A collection structure to inform the callback function with
  ! additional information. 
  type(t_collection), intent(inout), optional :: rcollection

  ! Optional: A combination of BCASM_DISCFORxxx constants that specify
  ! the complexity of the discretisation that is to perform. This allows to
  ! discretise only parts of the BC`s, e.g. only setting up those
  ! information that are necessary for filtering defect vectors.
  ! If not specified, BCASM_DISCFORALL is assumed, i.e. the resulting
  ! boundary conditions can be used for everything.
  integer(I32), intent(in), optional :: ccomplexity
!</input>  

!<inputoutput>
  ! This structure receives the result of the discretisation of rbcRegion.
  ! When entering the routine, the content of this structure is undefined,
  ! all pointers are invalid. The routine fills everything with appropriate
  ! data.
  type(t_discreteFBC), intent(inout), target :: rdiscreteFBC
!</inputoutput>

!</subroutine>

    ! local variables
    integer(I32) :: casmComplexity
    type(t_triangulation), pointer :: p_rtriangulation
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscr
    type(t_discreteFBCEntry), pointer :: p_rdiscreteFBCentry
    type(t_elementDistribution), pointer :: p_relementDist
    type(t_discreteFBCDirichlet),pointer :: p_rdirichletFBCs
    integer :: nequations,ieq,nmaxInfoPerElement,icount,i,j,iel,ielidx
    integer, dimension(:), pointer :: p_Ielements
    integer, dimension(:,:), allocatable :: Idofs
    type(t_discreteFBCevaluation), dimension(:), pointer :: p_Revaluation
    integer :: h_Ddofs,h_Idofs
    integer, dimension(:), pointer :: p_Idofs
    real(DP), dimension(:,:), pointer   :: p_Ddofs
    integer :: idof,iidx,nDOFs,ieldist,nve,ndofloc
    integer :: isubsetStart,isubsetLength
    integer, dimension(2) :: IdofCount
    real(dp), dimension(:,:), pointer :: p_Dwhere
    integer, dimension(:), pointer :: p_Iwhere,p_Iinside
    logical :: bok
    
    real(dp), dimension(:,:), pointer :: p_DvertexCoords
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:,:), pointer :: p_IedgesAtElement
    integer, dimension(:,:), pointer :: p_IverticesAtEdge
    integer, dimension(:,:), pointer :: p_IverticesAtFace
    
    casmComplexity = BCASM_DISCFORALL
    if (present(ccomplexity)) casmComplexity = ccomplexity
    
    ! Get the discretisation structure of the first component.
    ! In this first rough implementation, we assume that all equations
    ! are discretised with the same discretisation!
    ! Furthermore, Dirichlet values can only be specified for all solution
    ! components; something like 'slip' where one component is left free
    ! is not supported!
    
    p_rspatialDiscr => rblockDiscretisation%RspatialDiscr(Iequations(1))

    ! For easier access:
    p_rtriangulation => p_rspatialDiscr%p_rtriangulation
    call storage_getbase_double2d (p_rtriangulation%h_DvertexCoords,p_DvertexCoords)
    call storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,p_IverticesAtElement)

    ! Get a new FBC entry
    call bcasm_newFBCentry (rdiscreteFBC,iidx)
    p_rdiscreteFBCentry => rdiscreteFBC%p_RdiscFBCList(iidx)
    
    ! We have Dirichlet boundary conditions
    p_rdiscreteFBCentry%itype = DISCFBC_TPDIRICHLET
    
    ! Fill the structure for discrete Dirichlet BC`s in the
    ! t_discreteBCEntry structure
    p_rdirichletFBCs => p_rdiscreteFBCentry%rdirichletFBCs
    
    nequations = size(Iequations)
    p_rdirichletFBCs%ncomponents = nequations
    allocate(p_rdirichletFBCs%Icomponents(1:nequations))
    p_rdirichletFBCs%Icomponents(1:nequations) = Iequations(1:nequations)
    
    ! We have to collect all DOF`s and their values that belong to our current
    ! fictitious boundary object. Depending on the element type we will loop
    ! through all vertices, edges and elements in the triangulation
    ! to collect all the important values.
    !
    ! Allocate an array as large as a solution vector would be to store the DOF
    ! values. Allocate an integer array of the same size that receives a flag
    ! whether the DOF is actually used by our current FB object.
    ! More specifically, we remember the DOF value for each of the components!
    
    nDOFs = dof_igetNDofGlob(p_rspatialDiscr)
    
    IdofCount = (/nequations,nDOFs/)
    if (iand(casmComplexity,not(BCASM_DISCFORDEFMAT)) .ne. 0) then
      call storage_new ('bcasm_discrFBCDirichlet', 'DofValues', &
                        IdofCount, ST_DOUBLE, h_Ddofs, &
                        ST_NEWBLOCK_NOINIT)
      call storage_getbase_double2d (h_Ddofs,p_Ddofs)
    end if
    
    ! Initialise p_Idofs with zero.
    call storage_new ('bcasm_discrFBCDirichlet', 'DofUsed', nDOFs, &
                      ST_INT, h_Idofs, ST_NEWBLOCK_ZERO)
    call storage_getbase_int (h_Idofs,p_Idofs)
    
    ! At most, we calculate nmaxInfoPerElement pieces of information
    ! simultaneously.
    select case (p_rtriangulation%ndim)
    case (NDIM1D)
      ! 2 corner vertices or 1 edge midpoint
      ! -> at most 2 per call
      nmaxInfoPerElement = 2
    case (NDIM2D)
      ! 4 corner vertices or 4 edges or 1 element midpoint
      ! -> at most 4 per call
      nmaxInfoPerElement = 4
    case (NDIM3D)
      ! 8 corner vertices or 12 edges or 6 faces or 1 element midpoint
      ! -> at most 12 per call
      nmaxInfoPerElement = 12
    end select
    
    ! Initialise a p_Revaluation structure array for the evaluation
    allocate(p_Revaluation(nequations))
    allocate(p_Iwhere(FBCASM_MAXSIM*nmaxInfoPerElement))
    allocate(p_Iinside(FBCASM_MAXSIM*nmaxInfoPerElement))
    allocate(p_Dwhere(p_rtriangulation%ndim,FBCASM_MAXSIM*nmaxInfoPerElement))
    do ieq=1,nequations
      allocate(p_Revaluation(ieq)%p_Dvalues(FBCASM_MAXSIM*nmaxInfoPerElement,1))
    end do
    do ieq=1,nequations
      p_Revaluation(ieq)%p_Iinside => p_Iinside
      p_Revaluation(ieq)%p_Iwhere => p_Iwhere
      p_Revaluation(ieq)%p_Dwhere => p_Dwhere
    end do
    
    ! To collect the values of all the DOF's in question, we loop in sets
    ! about the elements in every element distribution.
    !
    ! So at first, we loop over the element distribution which tells us the
    ! current element type.
    do ieldist = 1,p_rspatialDiscr%inumFESpaces
    
      ! Get the element distribution
      p_relementDist => p_rspatialDiscr%RelementDistr(ieldist)
      
      if (p_relementDist%NEL .eq. 0) cycle
      
      ! Number of corner vertices.
      nve = elem_igetNVE(p_relementDist%celement)

      ! Reserve some memory to save temporarily all DOF`s of all elements.
      ndofloc = elem_igetNDofLoc(p_relementDist%celement)
      allocate (Idofs(ndofloc,FBCASM_MAXSIM))
      Idofs(:,:) = 0
      
      ! Get the elements in that distribution
      call storage_getbase_int(p_relementDist%h_IelementList,p_Ielements)
      
      ! Loop through the elements in sets a FBCASM_MAXSIM.
      do isubsetStart = 1,p_relementDist%NEL,FBCASM_MAXSIM
      
        ! Remaining elements here...
        isubsetLength = min(p_relementDist%NEL-isubsetStart+1,FBCASM_MAXSIM)
        
        ! Get all the DOF's on all the elements in the current set.
        call dof_locGlobMapping_mult(p_rspatialDiscr, &
                  p_Ielements(isubsetStart:isubsetStart+isubsetLength-1), Idofs)
                  
        ! Now the element dependent part. What's the dimension of the current space?
        bok = .false.
        select case (p_rspatialDiscr%ndimension)
        case (NDIM1D)
        case (NDIM2D)
        
          call storage_getbase_int2d (p_rtriangulation%h_IedgesAtElement,p_IedgesAtElement)
          call storage_getbase_int2d (p_rtriangulation%h_IverticesAtEdge,p_IverticesAtEdge)
        
          ! Element type? That decides on how to handle the DOF values.
          !
          ! -----
          ! P1, Q1, P2, Q2
          if ((p_relementDist%celement .eq. EL_P1) .or. &
              (p_relementDist%celement .eq. EL_Q1) .or. &
              (p_relementDist%celement .eq. EL_P2) .or. &
              (p_relementDist%celement .eq. EL_Q2)) then
              
            bok = .true.
              
            ! Loop through the vertices of the elements and collect them.
            ! The local DOF's 1..nve correspond to the corner vertices here.
            do ielidx = 1,isubsetLength
              iel = p_Ielements(isubsetStart+ielidx-1)
              ! Index of the entry
              icount = (ielidx-1)*nve
              do j=1,nve
                ! Initialise the data
                p_Iwhere(icount+1:icount+nve) = p_IverticesAtElement(1:nve,iel)
                do i=1,p_rtriangulation%ndim
                  p_Dwhere(i,icount+j) = p_DvertexCoords(i,p_Iwhere(icount+j))
                end do
              end do
            end do
          
            ! Call the callback routine to calculate the values.
            p_Revaluation(:)%cinfoNeeded = DISCFBC_NEEDFUNC
            p_Revaluation(:)%nvalues = isubsetLength*nve
            p_Iinside(1:isubsetLength*nve) = 0
            call fgetBoundaryValuesFBC (p_rdirichletFBCs%Icomponents(1:nequations),&
                                        rblockDiscretisation,&
                                        p_Revaluation, rcollection)
                                        
            ! Save the computed values to the corresponding DOF's.
            ! Already processed DOF's are strictily positive, they have
            ! the DOF number in p_Idofs!
            do ielidx = 1,isubsetLength
              iel = p_Ielements(isubsetStart+ielidx-1)
              ! Index of the entry
              icount = (ielidx-1)*nve
              do j=1,nve
                idof = Idofs(j,ielidx)
                if ((p_Iinside(icount+j) .ne. 0) .and. (p_Idofs(idof) .eq. 0)) then
                  do i=1,nequations
                    p_Ddofs(i,idof) = p_Revaluation(i)%p_Dvalues(icount+j,1)
                  end do
                  p_Idofs(idof) = idof
                end if
              end do
            end do
          end if
          
          if ((p_relementDist%celement .eq. EL_P2) .or. &
              (p_relementDist%celement .eq. EL_Q2)) then

            ! Loop through the edges of the elements and collect them.
            ! The local DOF's nve+1..2*nve correspond to the egdes here.
            do ielidx = 1,isubsetLength
              iel = p_Ielements(isubsetStart+ielidx-1)
              ! Index of the entry
              icount = (ielidx-1)*nve
              do j=1,nve
                ! Initialise the data
                p_Iwhere(icount+1:icount+nve) = p_IedgesAtElement(1:nve,iel)
                do i=1,p_rtriangulation%ndim
                  p_Dwhere(i,icount+j) = 0.5_DP * &
                      ( p_DvertexCoords(i,p_IverticesAtEdge(1,p_Iwhere(icount+j))) &
                      + p_DvertexCoords(i,p_IverticesAtEdge(2,p_Iwhere(icount+j))) )
                end do
              end do
            end do
          
            ! Call the callback routine to calculate the values.
            p_Revaluation(:)%cinfoNeeded = DISCFBC_NEEDFUNCMID
            p_Revaluation(:)%nvalues = isubsetLength*nve
            p_Iinside(1:isubsetLength*nve) = 0
            call fgetBoundaryValuesFBC (p_rdirichletFBCs%Icomponents(1:nequations),&
                                        rblockDiscretisation,&
                                        p_Revaluation, rcollection)
                                        
            ! Save the computed values to the corresponding DOF's.
            ! Already processed DOF's are strictily positive, they have
            ! the DOF number in p_Idofs!
            do ielidx = 1,isubsetLength
              iel = p_Ielements(isubsetStart+ielidx-1)
              ! Index of the entry
              icount = (ielidx-1)*nve
              do j=1,nve
                idof = Idofs(j+nve,ielidx)
                if ((p_Iinside(icount+j) .ne. 0) .and. (p_Idofs(idof) .eq. 0)) then
                  do i=1,nequations
                    p_Ddofs(i,idof) = p_Revaluation(i)%p_Dvalues(icount+j,1)
                  end do
                  p_Idofs(idof) = idof
                end if
              end do
            end do
          end if

          if (p_relementDist%celement .eq. EL_Q2) then

            ! Getthe element midpoints, corresponding to the local DOF 2*nve+1
            do ielidx = 1,isubsetLength
              iel = p_Ielements(isubsetStart+ielidx-1)
              ! Index of the entry
              icount = ielidx

              ! Initialise the data
              p_Iwhere(icount) = iel
              do i=1,p_rtriangulation%ndim
                p_Dwhere(i,icount) = 0.25_DP * &
                                      ( p_DvertexCoords(i,p_IverticesAtElement(1,iel)) &
                                      + p_DvertexCoords(i,p_IverticesAtElement(2,iel)) &
                                      + p_DvertexCoords(i,p_IverticesAtElement(3,iel)) &
                                      + p_DvertexCoords(i,p_IverticesAtElement(4,iel)) )
              end do
            end do
          
            ! Call the callback routine to calculate the values.
            p_Revaluation(:)%cinfoNeeded = DISCFBC_NEEDFUNCELMID
            p_Revaluation(:)%nvalues = isubsetLength
            p_Iinside(1:isubsetLength*nve) = 0
            call fgetBoundaryValuesFBC (p_rdirichletFBCs%Icomponents(1:nequations),&
                                        rblockDiscretisation,&
                                        p_Revaluation, rcollection)
                                        
            ! Save the computed values to the corresponding DOF's.
            ! Already processed DOF's are strictily positive, they have
            ! the DOF number in p_Idofs!
            do ielidx = 1,isubsetLength
              iel = p_Ielements(isubsetStart+ielidx-1)
              ! Index of the entry
              icount = ielidx

              idof = Idofs(2*nve+1,ielidx)
              if ((p_Iinside(icount) .ne. 0) .and. (p_Idofs(idof) .eq. 0)) then
                do i=1,nequations
                  p_Ddofs(i,idof) = p_Revaluation(i)%p_Dvalues(icount,1)
                end do
                p_Idofs(idof) = idof
              end if
            end do
          end if
          
          ! -----
          ! Q1~, all implementations
          if (elem_getPrimaryElement(p_relementDist%celement) .eq. EL_Q1T .or.&
              elem_getPrimaryElement(p_relementDist%celement) .eq. EL_P1T) then

            bok = .true.

            ! Loop through the edges of the elements and collect them.
            ! The local DOF's nve+1..2*nve correspond to the egdes here.
            do ielidx = 1,isubsetLength
              iel = p_Ielements(isubsetStart+ielidx-1)
              ! Index of the entry
              icount = (ielidx-1)*nve
              do j=1,nve
                ! Initialise the data
                p_Iwhere(icount+1:icount+nve) = p_IedgesAtElement(1:nve,iel)
                do i=1,p_rtriangulation%ndim
                  p_Dwhere(i,icount+j) = 0.5_DP * &
                      ( p_DvertexCoords(i,p_IverticesAtEdge(1,p_Iwhere(icount+j))) &
                      + p_DvertexCoords(i,p_IverticesAtEdge(2,p_Iwhere(icount+j))) )
                end do
              end do
            end do
          
            ! Call the callback routine to calculate the values.
            p_Revaluation(:)%cinfoNeeded = DISCFBC_NEEDFUNCMID
            if (iand(p_relementDist%celement,int(2**16,I32)) .ne. 0) then
              ! Integral mean value based element
              p_Revaluation(:)%cinfoNeeded = DISCFBC_NEEDINTMEAN
            end if
            p_Revaluation(:)%nvalues = isubsetLength*nve
            p_Iinside(1:isubsetLength*nve) = 0
            call fgetBoundaryValuesFBC (p_rdirichletFBCs%Icomponents(1:nequations),&
                                        rblockDiscretisation,&
                                        p_Revaluation, rcollection)
                                        
            ! Save the computed values to the corresponding DOF's.
            ! Already processed DOF's are strictily positive, they have
            ! the DOF number in p_Idofs!
            do ielidx = 1,isubsetLength
              iel = p_Ielements(isubsetStart+ielidx-1)
              ! Index of the entry
              icount = (ielidx-1)*nve
              do j=1,nve
                idof = Idofs(j,ielidx)
                if ((p_Iinside(icount+j) .ne. 0) .and. (p_Idofs(idof) .eq. 0)) then
                  do i=1,nequations
                    p_Ddofs(i,idof) = p_Revaluation(i)%p_Dvalues(icount+j,1)
                  end do
                  p_Idofs(idof) = idof
                end if
              end do
            end do
          end if
          
        case (NDIM3D)

          call storage_getbase_int2d (p_rtriangulation%h_IverticesAtFace,p_IverticesAtFace)

          ! Element type? That decides on how to handle the DOF values.
          !
          ! -----
          ! P1, Q1, P2, Q2.
          ! Note: P2/Q2 not yet supported... too ugly to realise until now,
          ! but will have to be done similar to 2D.
          if ((p_relementDist%celement .eq. EL_P1_3D) .or. &
              (p_relementDist%celement .eq. EL_Q1_3D)) then
              
            bok = .true.
              
            ! Loop through the vertices of the elements and collect them.
            ! The local DOF's 1..nve correspond to the corner vertices here.
            do ielidx = 1,isubsetLength
              iel = p_Ielements(isubsetStart+ielidx-1)
              ! Index of the entry
              icount = (ielidx-1)*nve
              do j=1,nve
                ! Initialise the data
                p_Iwhere(icount+1:icount+nve) = p_IverticesAtElement(1:nve,iel)
                do i=1,p_rtriangulation%ndim
                  p_Dwhere(i,icount+j) = p_DvertexCoords(i,p_Iwhere(icount+j))
                end do
              end do
            end do
          
            ! Call the callback routine to calculate the values.
            p_Revaluation(:)%cinfoNeeded = DISCFBC_NEEDFUNC
            p_Revaluation(:)%nvalues = isubsetLength*nve
            p_Iinside(1:isubsetLength*nve) = 0
            call fgetBoundaryValuesFBC (p_rdirichletFBCs%Icomponents(1:nequations),&
                                        rblockDiscretisation,&
                                        p_Revaluation, rcollection)
                                        
            ! Save the computed values to the corresponding DOF's.
            ! Already processed DOF's are strictily positive, they have
            ! the DOF number in p_Idofs!
            do ielidx = 1,isubsetLength
              iel = p_Ielements(isubsetStart+ielidx-1)
              ! Index of the entry
              icount = (ielidx-1)*nve
              do j=1,nve
                idof = Idofs(j,ielidx)
                if ((p_Iinside(icount+j) .ne. 0) .and. (p_Idofs(idof) .eq. 0)) then
                  do i=1,nequations
                    p_Ddofs(i,idof) = p_Revaluation(i)%p_Dvalues(icount+j,1)
                  end do
                  p_Idofs(idof) = idof
                end if
              end do
            end do
          end if
          
          ! -----
          ! Q1~, all implementations
          if (elem_getPrimaryElement(p_relementDist%celement) .eq. EL_Q1T_3D) then

            bok = .true.

            ! Loop through the faces of the elements and collect them.
            ! The local DOF's nve+1..2*nve correspond to the egdes here.
            do ielidx = 1,isubsetLength
              iel = p_Ielements(isubsetStart+ielidx-1)
              ! Index of the entry
              icount = (ielidx-1)*nve
              do j=1,nve
                ! Initialise the data
                p_Iwhere(icount+1:icount+nve) = p_IedgesAtElement(1:nve,iel)
                do i=1,p_rtriangulation%ndim
                  p_Dwhere(i,icount+j) = 0.25_DP * &
                      ( p_DvertexCoords(i,p_IverticesAtFace(1,p_Iwhere(icount+j))) &
                      + p_DvertexCoords(i,p_IverticesAtFace(2,p_Iwhere(icount+j))) &
                      + p_DvertexCoords(i,p_IverticesAtFace(3,p_Iwhere(icount+j))) &
                      + p_DvertexCoords(i,p_IverticesAtFace(4,p_Iwhere(icount+j))) )
                end do
              end do
            end do
          
            ! Call the callback routine to calculate the values.
            p_Revaluation(:)%cinfoNeeded = DISCFBC_NEEDFUNCFACEMID
            if (iand(p_relementDist%celement,int(2**16,I32)) .ne. 0) then
              ! Integral mean value based element
              p_Revaluation(:)%cinfoNeeded = DISCFBC_NEEDFACEINTMEAN
            end if
            p_Revaluation(:)%nvalues = isubsetLength*nve
            p_Iinside(1:isubsetLength*nve) = 0
            call fgetBoundaryValuesFBC (p_rdirichletFBCs%Icomponents(1:nequations),&
                                        rblockDiscretisation,&
                                        p_Revaluation, rcollection)
                                        
            ! Save the computed values to the corresponding DOF's.
            ! Already processed DOF's are strictily positive, they have
            ! the DOF number in p_Idofs!
            do ielidx = 1,isubsetLength
              iel = p_Ielements(isubsetStart+ielidx-1)
              ! Index of the entry
              icount = (ielidx-1)*nve
              do j=1,nve
                idof = Idofs(j,ielidx)
                if ((p_Iinside(icount+j) .ne. 0) .and. (p_Idofs(idof) .eq. 0)) then
                  do i=1,nequations
                    p_Ddofs(i,idof) = p_Revaluation(i)%p_Dvalues(icount+j,1)
                  end do
                  p_Idofs(idof) = idof
                end if
              end do
            end do
          end if
                    
        end select
        
        if (.not. bok) then
          call output_line (&
              'Element space not supported!', &
              OU_CLASS_ERROR,OU_MODE_STD,'bcasm_newDirichletBConFBD')
          call sys_halt()
        end if

      end do ! isubsetStart
    
    end do

    ! Release memory
    do ieq=1,nequations
      deallocate(p_Revaluation(ieq)%p_Dvalues)
    end do
    deallocate(p_Iinside)
    deallocate(p_Iwhere)
    deallocate(p_Dwhere)
    deallocate(p_Revaluation)

    ! Now, compress the p_Idofs/p_Ddofs array on to those DOF's
    ! which must be set. These are those DOF's where p_Ddofs is positive.
    icount = 0
    do i=1,size(p_Idofs)
      if (p_Idofs(i) .gt. 0) then
        icount = icount + 1
        p_Idofs(icount) = p_Idofs(i)
        p_Ddofs(:,icount) = p_Ddofs(:,i)
      end if
    end do

    ! Cancel if we did not find any DOF.
    if (icount .gt. 0) then
      
      ! Reallocate to save memory. Store the final handles in the structure.
      if (iand(casmComplexity,not(BCASM_DISCFORDEFMAT)) .ne. 0) then
        ! In the 2D-array, the size of the 2nd dimension is changed to the
        ! number of DOF`s.
        call storage_realloc ('bcasm_discrFBCDirichlet', icount, &
                              h_Ddofs, ST_NEWBLOCK_NOINIT)
        p_rdirichletFBCs%h_DdirichletValues = h_Ddofs
      end if
      
      call storage_realloc ('bcasm_discrFBCDirichlet', icount, &
                            h_Idofs, ST_NEWBLOCK_NOINIT)
      p_rdirichletFBCs%h_IdirichletDOFs = h_Idofs
    else
      ! Nothging inside; release arrays
      if (iand(casmComplexity,not(BCASM_DISCFORDEFMAT)) .ne. 0) then
        call storage_free (h_Ddofs)
      end if
      call storage_free (h_Idofs)
    end if

    p_rdirichletFBCs%nDOF = icount

  end subroutine



!    ! local variables
!    
!    integer :: nDOFs
!    integer :: h_Ddofs, h_Idofs, i, j, iidx, ieldist
!    integer(I32) :: celement
!    integer :: nequations
!    integer, dimension(2) :: IdofCount
!    
!    integer, dimension(:), pointer :: p_Ielements
!    integer :: ndofloc, nve
!    
!    integer, dimension(:), pointer :: p_Idofs
!    real(DP), dimension(:,:), pointer   :: p_Ddofs
!    
!    integer, dimension(:), pointer :: p_IdofUsed
!    type(t_directAccessIntSet) :: rset
!
!    
!    integer :: isubsetStart, isubsetLength, icurrentDof
!    
!    type(t_discreteFBCevaluation), dimension(DISCFBC_MAXDISCBC) :: Revaluation
!    
!    
!
!    
!    if (p_rspatialDiscr%ccomplexity .ne. SPDISC_UNIFORM) then
!      call output_line (&
!          'Element space not supported!', &
!          OU_CLASS_ERROR,OU_MODE_STD,'bcasm_newDirichletBConFBD')
!      print *,'bcasm_discrFBCDirichlet: Can only handle uniform discretisation!'
!      call sys_halt()
!    end if
!
!
!    ! All elements are of the samne type. Get it in advance.
!    celement = p_rspatialDiscr%RelementDistr(1)%celement
!    
!    ! Depending on the element type, prepare the evaluation structure and call
!    ! the callback routine to calculate what we need.
!    
!    icurrentDof = 0
!    
!    if(p_rspatialDiscr%ndimension .eq. NDIM2D) then
!    
!      ! Calculate values in the vertices for Q1,Q2,P1,P2
!      if ((celement .eq. EL_P1) .or. &
!          (celement .eq. EL_Q1) .or. &
!          (celement .eq. EL_P2) .or. &
!          (celement .eq. EL_Q2)) then
!      
!        ! Let us start to collect values. This is a rather element-dependent
!        ! part. At first, loop through the vertices in case we have a
!        ! P1/P2/Q1/Q2 discretisation
!        do isubsetStart = 1,p_rtriangulation%NVT,FBCASM_MAXSIM
!        
!          isubsetLength = min(p_rtriangulation%NVT-isubsetStart+1,FBCASM_MAXSIM)
!        
!          ! Fill the subset with isubsetStart, isubsetStart+1,... to identify the
!          ! subset we evaluate.
!          call fillsubset (isubsetStart,isubsetLength,Isubset)
!          
!          ! Fill the evaluation structure with data for the callback routine
!          do i=1,nequations
!            Revaluation(i)%cinfoNeeded = DISCFBC_NEEDFUNC
!            Revaluation(i)%nvalues     = isubsetLength
!            Revaluation(i)%p_Iwhere    => Isubset
!            Revaluation(i)%p_Dvalues   => p_Dsubset(:,:,i)
!            Revaluation(i)%p_Iinside   => Iinside
!          end do
!
!          ! Clear the Iinside array
!          Iinside = 0
!          
!          ! Call the callback routine to calculate the values.
!          call fgetBoundaryValuesFBC (p_rdirichletFBCs%Icomponents(1:nequations),&
!                                      rblockDiscretisation,&
!                                      Revaluation, rcollection)
!                                      
!          ! Transfer the DOF`s that are affected
!
!          if (iand(casmComplexity,not(BCASM_DISCFORDEFMAT)) .ne. 0) then
!            
!            do i=1,isubsetLength
!              if (Iinside(i) .ne. 0) then
!                icurrentDof = icurrentDof + 1
!                do j=1,nequations
!                  p_Ddofs(j,icurrentDof) = p_Dsubset(i,1,j)
!                end do
!                p_Idofs(icurrentDof) = isubsetStart+i-1
!              end if
!            end do
!            
!          else
!          
!            do i=1,isubsetLength
!              if (Iinside(i) .ne. 0) then
!                icurrentDof = icurrentDof + 1
!                p_Idofs(icurrentDof) = isubsetStart+i-1
!              end if
!            end do
!            
!          end if
!        
!        end do
!      
!        ! In case of a P2/Q2 discretisation, also loop about the edges.
!        if ((celement .eq. EL_P2) .or. &
!            (celement .eq. EL_Q2)) then
!        
!          ! Let us start to collect values. This is a rather element-dependent
!          ! part. At first, loop through the vertices in case we have a
!          ! P1/P2/Q1/Q2 discretisation
!          do isubsetStart = 1,p_rtriangulation%NMT,FBCASM_MAXSIM
!          
!            isubsetLength = min(p_rtriangulation%NMT-isubsetStart+1,FBCASM_MAXSIM)
!          
!            ! Fill the subset with isubsetStart, isubsetStart+1,... to identify the
!            ! subset we evaluate.
!            call fillsubset (isubsetStart,isubsetLength,Isubset)
!            
!            ! Fill the evaluation structure with data for the callback routine
!            do i=1,nequations
!              Revaluation(i)%cinfoNeeded = DISCFBC_NEEDFUNCMID
!              Revaluation(i)%nvalues     = isubsetLength
!              Revaluation(i)%p_Iwhere    => Isubset
!              Revaluation(i)%p_Dvalues   => p_Dsubset(:,:,i)
!              Revaluation(i)%p_Iinside   => Iinside
!            end do
!
!            ! Clear the Iinside array
!            Iinside = 0
!            
!            ! Call the callback routine to calculate the values.
!            call fgetBoundaryValuesFBC (p_rdirichletFBCs%Icomponents(1:nequations),&
!                                        rblockDiscretisation,&
!                                        Revaluation, rcollection)
!                                        
!            ! Transfer the DOF`s that are affected
!
!            if (iand(casmComplexity,not(BCASM_DISCFORDEFMAT)) .ne. 0) then
!              
!              do i=1,isubsetLength
!                if (Iinside(i) .ne. 0) then
!                  icurrentDof = icurrentDof + 1
!                  do j=1,nequations
!                    p_Ddofs(j,icurrentDof) = p_Dsubset(i,1,j)
!                  end do
!                  ! Save the DOF ID
!                  p_Idofs(icurrentDof) = p_rtriangulation%NVT+isubsetStart+i-1
!                end if
!              end do
!              
!            else
!            
!              do i=1,isubsetLength
!                if (Iinside(i) .ne. 0) then
!                  icurrentDof = icurrentDof + 1
!                  ! Save the DOF ID
!                  p_Idofs(icurrentDof) = p_rtriangulation%NVT+isubsetStart+i-1
!                end if
!              end do
!              
!            end if
!          
!          end do
!        
!        end if 
!
!        ! In case of a Q2 discretisation, also loop about the element midpoints.
!        if ((celement .eq. EL_Q2)) then
!        
!          ! Let us start to collect values. This is a rather element-dependent
!          ! part. At first, loop through the vertices in case we have a
!          ! P1/P2/Q1/Q2 discretisation
!          do isubsetStart = 1,p_rtriangulation%NEL,FBCASM_MAXSIM
!          
!            isubsetLength = min(p_rtriangulation%NEL-isubsetStart+1,FBCASM_MAXSIM)
!          
!            ! Fill the subset with isubsetStart, isubsetStart+1,... to identify the
!            ! subset we evaluate.
!            call fillsubset (isubsetStart,isubsetLength,Isubset)
!            
!            ! Fill the evaluation structure with data for the callback routine
!            do i=1,nequations
!              Revaluation(i)%cinfoNeeded = DISCFBC_NEEDFUNCELMID
!              Revaluation(i)%nvalues     = isubsetLength
!              Revaluation(i)%p_Iwhere    => Isubset
!              Revaluation(i)%p_Dvalues   => p_Dsubset(:,:,i)
!              Revaluation(i)%p_Iinside   => Iinside
!            end do
!
!            ! Clear the Iinside array
!            Iinside = 0
!            
!            ! Call the callback routine to calculate the values.
!            call fgetBoundaryValuesFBC (p_rdirichletFBCs%Icomponents(1:nequations),&
!                                        rblockDiscretisation,&
!                                        Revaluation, rcollection)
!                                        
!            ! Transfer the DOF`s that are affected
!
!            if (iand(casmComplexity,not(BCASM_DISCFORDEFMAT)) .ne. 0) then
!              
!              do i=1,isubsetLength
!                if (Iinside(i) .ne. 0) then
!                  icurrentDof = icurrentDof + 1
!                  do j=1,nequations
!                    p_Ddofs(j,icurrentDof) = p_Dsubset(i,1,j)
!                  end do
!                  ! Save the DOF ID
!                  p_Idofs(icurrentDof) = p_rtriangulation%NVT+&
!                      p_rtriangulation%NMT+isubsetStart+i-1
!                end if
!              end do
!              
!            else
!            
!              do i=1,isubsetLength
!                if (Iinside(i) .ne. 0) then
!                  icurrentDof = icurrentDof + 1
!                  ! Save the DOF ID
!                  p_Idofs(icurrentDof) = p_rtriangulation%NVT+p_rtriangulation%NMT+&
!                      isubsetStart+i-1
!                end if
!              end do
!              
!            end if
!          
!          end do
!        
!        end if 
!      
!      end if ! end if el_typ
!    
!      ! Calculate values in the face midpoints / / integral mean values for Q1~
!      if (elem_getPrimaryElement(celement) .eq. EL_Q1T) then
!      
!        ! Let us start to collect values. This is a rather element-dependent
!        ! part. At first, loop through the vertices in case we have a
!        ! P1/Q1/Q2 discretisation
!        do isubsetStart = 1,p_rtriangulation%NMT,FBCASM_MAXSIM
!        
!          isubsetLength = min(p_rtriangulation%NMT-isubsetStart+1,FBCASM_MAXSIM)
!        
!          ! Fill the subset with isubsetStart, isubsetStart+1,... to identify the
!          ! subset we evaluate.
!          call fillsubset (isubsetStart,isubsetLength,Isubset)
!          
!          ! Fill the evaluation structure with data for the callback routine
!          do i=1,nequations
!            if ((celement .eq. EL_E031) .or. &
!                (celement .eq. EL_EM31)) then
!              Revaluation(i)%cinfoNeeded = DISCFBC_NEEDFUNCMID
!            else
!              Revaluation(i)%cinfoNeeded = DISCFBC_NEEDINTMEAN
!            end if
!            Revaluation(i)%nvalues     = isubsetLength
!            Revaluation(i)%p_Iwhere    => Isubset
!            Revaluation(i)%p_Dvalues   => p_Dsubset(:,:,i)
!            Revaluation(i)%p_Iinside   => Iinside
!          end do
!          
!          ! Clear the Iinside array
!          Iinside = 0
!          
!          ! Call the callback routine to calculate the values.
!          call fgetBoundaryValuesFBC (p_rdirichletFBCs%Icomponents(1:nequations),&
!                                      rblockDiscretisation,&
!                                      Revaluation, rcollection)
!                                      
!          ! Transfer the DOF`s that are affected
!
!          if (iand(casmComplexity,not(BCASM_DISCFORDEFMAT)) .ne. 0) then
!            
!            do i=1,isubsetLength
!              if (Iinside(i) .ne. 0) then
!                icurrentDof = icurrentDof + 1
!                do j=1,nequations
!                  p_Ddofs(j,icurrentDof) = p_Dsubset(i,1,j)
!                end do
!                p_Idofs(icurrentDof) = isubsetStart+i-1
!              end if
!            end do
!            
!          else
!          
!            do i=1,isubsetLength
!              if (Iinside(i) .ne. 0) then
!                icurrentDof = icurrentDof + 1
!                p_Idofs(icurrentDof) = isubsetStart+i-1
!              end if
!            end do
!            
!          end if
!        
!        end do
!      
!      end if ! end if q1tilde
!    
!    end if ! end dim2d
!    
!    ! q2,q1,q1t
!    if(p_rspatialDiscr%ndimension .eq. NDIM3D)then
!      if(elem_getPrimaryElement(celement) .eq. EL_Q1_3D)then      
!      
!        ! Let us start to collect values. This is a rather element-dependent
!        ! part. At first, loop through the vertices in case we have a
!        ! Q1/Q2 discretisation
!        do isubsetStart = 1,p_rtriangulation%NVT,FBCASM_MAXSIM
!        
!          isubsetLength = min(p_rtriangulation%NVT-isubsetStart+1,FBCASM_MAXSIM)
!        
!          ! Fill the subset with isubsetStart, isubsetStart+1,... to identify the
!          ! subset we evaluate.
!          call fillsubset (isubsetStart,isubsetLength,Isubset)
!          
!          ! Fill the evaluation structure with data for the callback routine
!          do i=1,nequations
!            Revaluation(i)%cinfoNeeded = DISCFBC_NEEDFUNC
!            Revaluation(i)%nvalues     = isubsetLength
!            Revaluation(i)%p_Iwhere    => Isubset
!            Revaluation(i)%p_Dvalues   => p_Dsubset(:,:,i)
!            Revaluation(i)%p_Iinside   => Iinside
!          end do
!
!          ! Clear the Iinside array
!          Iinside = 0
!          
!          ! Call the callback routine to calculate the values.
!          call fgetBoundaryValuesFBC (p_rdirichletFBCs%Icomponents(1:nequations),&
!                                      rblockDiscretisation,&
!                                      Revaluation, rcollection)
!                                      
!          ! Transfer the DOF`s that are affected
!
!          if (iand(casmComplexity,not(BCASM_DISCFORDEFMAT)) .ne. 0) then
!            
!            do i=1,isubsetLength
!              if (Iinside(i) .ne. 0) then
!                icurrentDof = icurrentDof + 1
!                do j=1,nequations
!                  p_Ddofs(j,icurrentDof) = p_Dsubset(i,1,j)
!                end do
!                p_Idofs(icurrentDof) = isubsetStart+i-1
!              end if
!            end do
!            
!          else
!          
!            do i=1,isubsetLength
!              if (Iinside(i) .ne. 0) then
!                icurrentDof = icurrentDof + 1
!                p_Idofs(icurrentDof) = isubsetStart+i-1
!              end if
!            end do
!            
!          end if
!        
!        end do
!
!      end if
!      
!      ! Calculate values in the face midpoints / / integral mean values for Q1~
!      if(elem_getPrimaryElement(celement) .eq. EL_Q1T_3D)then
!      
!        ! Let us start to collect values. This is a rather element-dependent
!        ! part. 
!        do isubsetStart = 1,p_rtriangulation%NAT,FBCASM_MAXSIM
!        
!          isubsetLength = min(p_rtriangulation%NAT-isubsetStart+1,FBCASM_MAXSIM)
!        
!          ! Fill the subset with isubsetStart, isubsetStart+1,... to identify the
!          ! subset we evaluate.
!          call fillsubset (isubsetStart,isubsetLength,Isubset)
!          
!          ! Fill the evaluation structure with data for the callback routine
!          do i=1,nequations
!            if ((celement .eq. EL_E031) .or. &
!                (celement .eq. EL_EM31)) then
!              Revaluation(i)%cinfoNeeded = DISCFBC_NEEDFUNCFACEMID
!            else
!              Revaluation(i)%cinfoNeeded = DISCFBC_NEEDFACEINTMEAN
!            end if
!            Revaluation(i)%nvalues     = isubsetLength
!            Revaluation(i)%p_Iwhere    => Isubset
!            Revaluation(i)%p_Dvalues   => p_Dsubset(:,:,i)
!            Revaluation(i)%p_Iinside   => Iinside
!          end do
!          
!          ! Clear the Iinside array
!          Iinside = 0
!          
!          ! Call the callback routine to calculate the values.
!          call fgetBoundaryValuesFBC (p_rdirichletFBCs%Icomponents(1:nequations),&
!                                      rblockDiscretisation,&
!                                      Revaluation, rcollection)
!                                      
!          ! Transfer the DOF`s that are affected
!
!          if (iand(casmComplexity,not(BCASM_DISCFORDEFMAT)) .ne. 0) then
!            
!            do i=1,isubsetLength
!              if (Iinside(i) .ne. 0) then
!                icurrentDof = icurrentDof + 1
!                do j=1,nequations
!                  p_Ddofs(j,icurrentDof) = p_Dsubset(i,1,j)
!                end do
!                p_Idofs(icurrentDof) = isubsetStart+i-1
!              end if
!            end do
!            
!          else
!          
!            do i=1,isubsetLength
!              if (Iinside(i) .ne. 0) then
!                icurrentDof = icurrentDof + 1
!                p_Idofs(icurrentDof) = isubsetStart+i-1
!              end if
!            end do
!            
!          end if
!        
!        end do
!      
!      
!      end if ! end EL_Q1T_3D
!      
!    end if ! end if NDIM3d
!    
!    
!    
!    ! Cancel if we did not find any DOF.
!    if (icurrentDof .gt. 0) then
!      
!      ! Reallocate to save memory. Store the final handles in the structure.
!      if (iand(casmComplexity,not(BCASM_DISCFORDEFMAT)) .ne. 0) then
!        ! In the 2D-array, the size of the 2nd dimension is changed to the
!        ! number of DOF`s.
!        call storage_realloc ('bcasm_discrFBCDirichlet', icurrentDof, &
!                              h_Ddofs, ST_NEWBLOCK_NOINIT)
!        p_rdirichletFBCs%h_DdirichletValues = h_Ddofs
!      end if
!      
!      call storage_realloc ('bcasm_discrFBCDirichlet', icurrentDof, &
!                            h_Idofs, ST_NEWBLOCK_NOINIT)
!      p_rdirichletFBCs%h_IdirichletDOFs = h_Idofs
!    else
!      ! Nothging inside; release arrays
!      if (iand(casmComplexity,not(BCASM_DISCFORDEFMAT)) .ne. 0) then
!        call storage_free (h_Ddofs)
!      end if
!      call storage_free (h_Idofs)
!    end if
!
!    p_rdirichletFBCs%nDOF = icurrentDof
!    
!    ! Release temporary data
!    deallocate(p_Dsubset)
!    
!  contains
!  
!    pure subroutine fillsubset (istart, ilength, Isubset)
!    integer, intent(in) :: istart, ilength
!    integer, dimension(:), intent(out) :: Isubset
!    integer :: i
!      do i=1,ilength
!        Isubset(i) = istart-1+i
!      end do
!    end subroutine
!    
!  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine bcasm_releaseFBCDirichlet (rdiscreteFBCDirichlet)
  
!<description>
  ! This routine cleans up the discrete Dirichlet conditions for fictitious 
  ! boundary objects in rdiscreteFBCDirichlet.
!</description>

!<inputoutput>
  ! The discrete-FBC structure which is to be cleaned up
  type(t_discreteFBCDirichlet), intent(inout) :: rdiscreteFBCDirichlet
!</inputoutput>

!</subroutine>

    ! Release what is associated
    
    rdiscreteFBCDirichlet%nDOF = 0
    if (rdiscreteFBCDirichlet%h_DdirichletValues .ne. ST_NOHANDLE) &
      call storage_free(rdiscreteFBCDirichlet%h_DdirichletValues)
    if (rdiscreteFBCDirichlet%h_IdirichletDOFs .ne. ST_NOHANDLE) &
      call storage_free(rdiscreteFBCDirichlet%h_IdirichletDOFs)

    deallocate(rdiscreteFBCDirichlet%Icomponents)
    rdiscreteFBCDirichlet%ncomponents = 0

  end subroutine

  ! ***************************************************************************

!<subroutine>
  
  subroutine bcasm_newDirichletBConMR (rblockDiscretisation, iequation, &
                                       rdiscreteBC, rmeshRegion, &
                                       fgetBoundaryValuesMR, rcollection)
  
!<description>
  ! Adds a Dirichlet boundary condition entry to the discrete boundary condition
  ! structure, based on a mesh region.
!</description>

!<input>
  ! The block discretisation structure of the underlying PDE.
  type(t_blockDiscretisation), intent(in) :: rblockDiscretisation
  
  ! An identifier for the equation, this boundary condition refers to.
  ! >= 1. 1=first equation (e.g. X-velocity), 2=2nd equation (e.g. 
  ! Y-velocity), etc.
  integer, intent(in)                     :: iequation

  ! The mesh region structure that is used for the boundary region.
  type(t_meshRegion), intent(in)          :: rmeshRegion
  
  ! A callback function that calculates values on the boundary.
  ! Is declared in the interface include file 'intf_discretebc.inc'.
  include 'intf_discretebc.inc'

  ! Optional: A collection structure to inform the callback function with
  ! additional information. 
  type(t_collection), intent(inout), optional :: rcollection

!</input>

!<inputoutput>
  ! A t_discreteBC structures, representing the boundary discretised 
  ! in a discretisation-dependent way. The new BC`s are added to this structure.
  type(t_discreteBC), intent(inout) :: rdiscreteBC
!</inputoutput>

!</subroutine>

  ! More than a hand full of local variables
  integer :: i,j,ivt,iel,imt,iat,idx,idofHigh,&
    inumGlobalDofs,ilenDofBitmap,idof,cinfoNeeded,ndofs
  integer, dimension(1) :: Icomponents
  real(DP), dimension(1) :: Dvalues
  integer(I32) :: ielemType, idofMask
  type(t_discreteBCEntry), pointer :: p_rdbcEntry
  type(t_discreteBCDirichlet), pointer :: p_rdirichlet
  type(t_triangulation), pointer :: p_rtria
  type(t_spatialDiscretisation), pointer :: p_rspatDisc
  integer, dimension(:), pointer :: p_IvertexIdx
  integer, dimension(:), pointer :: p_IedgeIdx
  integer, dimension(:), pointer :: p_IfaceIdx
  integer, dimension(:), pointer :: p_IelementIdx
  integer :: iregionNVT, itriaNVT
  integer :: iregionNMT, itriaNMT
  integer :: iregionNAT, itriaNAT
  integer :: iregionNEL, itriaNEL
  integer, dimension(:), pointer :: p_IdofBitmap
  integer, dimension(:), pointer :: p_IelemAtVertIdx, p_IelemAtVert,&
    p_IelemAtEdgeIdx, p_IelemAtEdge, p_IelemDist
  integer, dimension(:,:), pointer :: p_IelemAtEdge2D, p_IelemAtFace,&
    p_IedgesAtElem, p_IfacesAtElem
  real(DP), dimension(1) :: Dcoord1D
  real(DP), dimension(2) :: Dcoord2D
  real(DP), dimension(3) :: Dcoord3D
  real(DP), dimension(:,:), pointer :: p_DvertexCoords
  integer, dimension(:,:), pointer :: p_IvertsAtEdge,&
    p_IvertsAtFace, p_IvertsAtElem
  integer, dimension(4) :: Iwhere
  real(DP), dimension(3) :: Dwhere
  integer, dimension(EL_MAXNBAS) :: IdofGlob
  logical :: buniform
  integer :: iidx
  real(DP) :: daux1, daux2, daux3, daux4
  
  real(DP), parameter :: Q12 = 0.5_DP
  real(DP), parameter :: Q13 = 0.333333333333333_DP
  real(DP), parameter :: Q14 = 0.25_DP
  real(DP), parameter :: Q18 = 0.125_DP
  
  real(DP), parameter :: G2P = 0.577350269189626_DP ! = sqrt(1/3)
  real(DP), parameter :: G3P = 0.774596669241483_DP ! = sqrt(3/5)
  
  ! Value of second Legendre-Polynomial in G3P
  real(DP), parameter :: L2G3P = 0.4_DP

    ! Get the spatial discretisation and the triangulation
    p_rspatDisc => rblockDiscretisation%RspatialDiscr(iequation)
    Icomponents(1) = iequation
    cinfoNeeded = DISCBC_NEEDFUNC

    ! Is this a uniform discretisation?
    if (p_rspatDisc%inumFESpaces .gt. 1) then
  
      ! Get the element distribution array, if given
      call storage_getbase_int(p_rspatDisc%h_IelementDistr, p_IelemDist)
      buniform = .false.
  
    else
  
      ! Get the one and only element type
      ielemType = p_rspatDisc%RelementDistr(1)%celement
      buniform = .true.
  
    end if
    
    ! Get the triangulation of the spatial discretisation
    p_rtria => p_rspatDisc%p_rtriangulation

    ! Make sure that the triangulation of the spatial discretisation and the
    ! mesh region are the same!
    if(.not. associated(p_rtria, rmeshRegion%p_rtriangulation)) then

      call output_line (&
          'Different triangulations for spatial discretisation and mesh region!', &
          OU_CLASS_ERROR,OU_MODE_STD,'bcasm_newDirichletBConMR')
      call sys_halt()

    end if
    
    ! Get the vertex coordinates and vertex index arrays from the triangulation
    call storage_getbase_double2D(p_rtria%h_DvertexCoords, p_DvertexCoords)
    call storage_getbase_int2D(p_rtria%h_IverticesAtElement, p_IvertsAtElem)
    call storage_getbase_int(p_rtria%h_IelementsAtVertexIdx, p_IelemAtVertIdx)
    call storage_getbase_int(p_rtria%h_IelementsAtVertex, p_IelemAtVert)
    if (p_rtria%ndim .eq. NDIM2D) then
      call storage_getbase_int2D(p_rtria%h_IverticesAtEdge, p_IvertsAtEdge)
      call storage_getbase_int2D(p_rtria%h_IelementsAtEdge, p_IelemAtEdge2D)
      call storage_getbase_int2D(p_rtria%h_IedgesAtElement, p_IedgesAtElem)
    else if (p_rtria%ndim .eq. NDIM3D) then
      call storage_getbase_int2D(p_rtria%h_IverticesAtEdge, p_IvertsAtEdge)
      call storage_getbase_int2D(p_rtria%h_IedgesAtElement, p_IedgesAtElem)
      call storage_getbase_int2D(p_rtria%h_IfacesAtElement, p_IfacesAtElem)
      call storage_getbase_int2D(p_rtria%h_IverticesAtFace, p_IvertsAtFace)
      call storage_getbase_int(p_rtria%h_IelementsAtEdgeIdx3D, p_IelemAtEdgeIdx)
      call storage_getbase_int(p_rtria%h_IelementsAtEdge3D, p_IelemAtEdge)
      call storage_getbase_int2D(p_rtria%h_IelementsAtFace, p_IelemAtFace)
    end if
    itriaNVT = p_rtria%NVT
    itriaNMT = p_rtria%NMT
    itriaNAT = p_rtria%NAT
    itriaNEL = p_rtria%NEL

    ! Now get all the arrays and counts from the mesh region
    if(rmeshRegion%h_IvertexIdx .ne. ST_NOHANDLE) then
      call storage_getbase_int(rmeshRegion%h_IvertexIdx, p_IvertexIdx)
      iregionNVT = rmeshRegion%NVT
    else
      p_IvertexIdx => null()
      iregionNVT = 0
    end if
    if(rmeshRegion%h_IedgeIdx .ne. ST_NOHANDLE) then
      call storage_getbase_int(rmeshRegion%h_IedgeIdx, p_IedgeIdx)
      iregionNMT = rmeshRegion%NMT
    else
      p_IedgeIdx => null()
      iregionNMT = 0
    end if
    if(rmeshRegion%h_IfaceIdx .ne. ST_NOHANDLE) then
      call storage_getbase_int(rmeshRegion%h_IfaceIdx, p_IfaceIdx)
      iregionNAT = rmeshRegion%NAT
    else
      p_IfaceIdx => null()
      iregionNAT = 0
    end if
    if(rmeshRegion%h_IelementIdx .ne. ST_NOHANDLE) then
      call storage_getbase_int(rmeshRegion%h_IelementIdx, p_IelementIdx)
      iregionNEL = rmeshRegion%NEL
    else
      p_IelementIdx => null()
      iregionNEL = 0
    end if

    ! Get a new entry for the BC`s.
    call bcasm_newBCentry (rdiscreteBC,iidx)
    
    ! Get the next discrete BC entry
    p_rdbcEntry => rdiscreteBC%p_RdiscBCList(iidx)
    
    ! Set the bc type to Dirichlet
    p_rdbcEntry%itype = DISCBC_TPDIRICHLET
    ! And get a pointer to the Dirichlet BC structure
    p_rdirichlet => p_rdbcEntry%rdirichletBCs
    
    ! set the component number and the number of DOFs
    p_rdirichlet%icomponent = iequation
    
    ! Get the total number of global DOFs
    inumGlobalDofs = dof_igetNDofGlob(p_rspatDisc)
    
    ! We have to ensure that one DOF does not get added to the dirichlet BC
    ! array twice. To ensure this, we will create an array called DOF-Bitmap.
    ! This is an array of 32-bit integers, but we will interpret it as an
    ! array of bits.
    ! Assume we are currently processing the DOF with number 'idof', then we 
    ! need to calculate 2 values to access the bit of the DOF in the bitmap:
    ! 1. idofHigh: the index of the 32-bit integer inside the bitmap where 
    !              our DOF is in
    ! 2. idofMask: a 32-bit-mask for the DOF inside that 32-bit integer
    !
    ! These two values are calculated as follows:
    ! idofHigh = ISHFT(idof-1,-5) + 1
    ! idofMask = INT(ISHFT(1,IAND(idof-1,31)),I32)
    !
    ! To mark a DOF as 'processed' we need to apply an OR-operator:
    ! p_IdofBitmap(idofHigh) = IOR(p_IdofBitmap(idofHigh),idofMask)
    !
    ! To check whether a DOF has already been processed we need to apply an
    ! AND-operator:
    ! IF (IAND(p_IdofBitmap(idofHigh),idofMask) .NE. 0) THEN
    !   ! DOF has already been processed
    ! END IF
    !
    ! Remark:
    ! The DOF bitmap works for both 32- and 64-bit systems. Remember that on
    ! a 64-Bit system 'idofHigh' is 64-Bit where 'idofMask' is 32-Bit.
    !
    
    ! The length of the dof bitmap is = inumGlobalDofs/32 + 1
    ilenDofBitmap = (inumGlobalDofs / 32) + 1
    
    ! Allocate a DOF bitmap
    allocate(p_IdofBitmap(1:ilenDofBitmap))
    
    ! Format the DOF bitmap to 0
    do i=1, ilenDofBitmap
      p_IdofBitmap(i) = 0
    end do
    
    ! Calculate an initial guess of the number of DOFs that will be affected
    ! with this dirichlet condition.
    ndofs = 0
    if (p_rspatDisc%inumFESpaces .eq. 1) then
    
      ! The discretisation is uniform. In this nice case, our initial
      ! guess will be exact.
      select case(elem_getPrimaryElement(&
        p_rspatDisc%RelementDistr(1)%celement))
      
      ! P0/Q0 elements (and 2D QP1 element)
      case (EL_P0_1D,EL_P0,EL_Q0,EL_QP1,EL_P0_3D,EL_Q0_3D,EL_Y0_3D,EL_R0_3D)
        ndofs = iregionNEL
      
      ! P1/Q1 elements (and 1D S31 element)
      case (EL_P1_1D,EL_S31_1D,EL_P1,EL_Q1,EL_P1_3D,EL_Q1_3D,EL_Y1_3D,EL_R1_3D)
        ndofs = iregionNVT
      
      ! 1D P2 element
      case (EL_P2_1D)
        ndofs = iregionNVT + iregionNEL
      
      ! 2D P2 element
      case (EL_P2)
        ndofs = iregionNVT + iregionNMT
      
      ! 2D Q2 element
      case (EL_Q2)
        ndofs = iregionNVT + iregionNMT + iregionNEL
            
      ! 2D P1~/Q1~ element
      case (EL_P1T,EL_Q1T)
        ndofs = iregionNMT
      
      ! 2D Q1~ with bubble element
      case (EL_Q1TB)
        ndofs = iregionNMT + iregionNEL
      
      ! 2D Q2~ element
      case (EL_Q2T)
        ndofs = 2*iregionNMT + iregionNEL
       
      ! 2D Q2~ with bubble element
      case (EL_Q2TB)
        ndofs = 2*iregionNMT + 2*iregionNEL
      
      ! 3D Q1~ element
      case (EL_Q1T_3D)
        ndofs = iregionNAT
      
      ! 3D Q2~ element
      case (EL_Q2T_3D)
        ndofs = 2*iregionNAT + iregionNEL
      
      end select

    end if
    
    ! If the discretisation is not uniform or we did not calculate an initial guess
    ! for whatever reason, we will set it to the maximum of the array lengths.
    if (ndofs .eq. 0) ndofs = max(iregionNVT,iregionNMT,iregionNAT,iregionNEL)
    
    ! Reset Iwhere
    Iwhere = 0
    
    ! First of all, go through all vertices in the mesh region (if any at all)
    do i = 1, iregionNVT
    
      ! Get the index of the vertex
      ivt = p_IvertexIdx(i)
      
      ! Store it to Iwhere
      Iwhere(1) = ivt
      
      ! And go through all elements which are adjacent to this vertex
      do idx=p_IelemAtVertIdx(ivt), p_IelemAtVertIdx(ivt+1)-1
      
        ! Get the index of the element
        iel = p_IelemAtVert(idx)
        
        ! Store the element number to Iwhere
        Iwhere(4) = iel
        
        ! Get all global dofs on this element
        call dof_locGlobMapping(p_rspatDisc, iel, IdofGlob)
        
        ! Now get the element type for this element (if the discretisation
        ! is not uniform)
        if (.not. buniform) ielemType = &
          p_rspatDisc%RelementDistr(p_IelemDist(iel))%celement
        
        ! Now we need to find which global DOF the currently processed
        ! vertices belongs to.
        idof = 0
        select case(elem_getPrimaryElement(ielemType))
        ! 1D P1/P2/S31 element
        case (EL_P1_1D,EL_P2_1D,EL_S31_1D)
          do j = 1, 2
            if (p_IvertsAtElem(j,iel) .eq. ivt) then
              idof = IdofGlob(j)
              exit
            end if
          end do
          
        ! 2D P1/P2/P3 element
        case (EL_P1,EL_P2,EL_P3)
          do j = 1, 3
            if (p_IvertsAtElem(j,iel) .eq. ivt) then
              idof = IdofGlob(j)
              exit
            end if
          end do

        ! 2D Q1/Q2/Q3 element
        case (EL_Q1,EL_Q2,EL_Q3)
          do j = 1, 4 
            if (p_IvertsAtElem(j,iel) .eq. ivt) then
              idof = IdofGlob(j)
              exit
            end if
          end do

        ! 3D P1 element
        case (EL_P1_3D)
          do j = 1, 4
            if (p_IvertsAtElem(j,iel) .eq. ivt) then
              idof = IdofGlob(j)
              exit
            end if
          end do

        ! 3D Q1/Q2 element
        case (EL_Q1_3D,EL_Q2_3D)
          do j = 1, 8
            if (p_IvertsAtElem(j,iel) .eq. ivt) then
              idof = IdofGlob(j)
              exit
            end if
          end do

        ! 3D Y1 element
        case (EL_Y1_3D)
          do j = 1, 5
            if (p_IvertsAtElem(j,iel) .eq. ivt) then
              idof = IdofGlob(j)
              exit
            end if
          end do

        ! 3D R1 element
        case (EL_R1_3D)
          do j = 1, 6
            if (p_IvertsAtElem(j,iel) .eq. ivt) then
              idof = IdofGlob(j)
              exit
            end if
          end do

        end select
        
        ! If we come out here and idof is 0, then either the element
        ! does not have DOFs in the vertices or something unforseen has
        ! happened...
        if (idof .eq. 0) cycle
        
        ! Let us check if we have already processed this dof.
        idofHigh = ishft(idof-1,-5) + 1
        idofMask = int(ishft(1,iand(idof-1,31)),I32)
        if (iand(p_IdofBitmap(idofHigh),int(idofMask)) .ne. 0) cycle
        
        ! This is a new DOF for the list - so call the boundary values callback
        ! routine to calculate the value
        call fgetBoundaryValuesMR(Icomponents, p_rspatDisc, rmeshRegion,&
            cinfoNeeded, Iwhere, Dwhere, p_DvertexCoords(:,ivt), Dvalues,&
            rcollection)
        
        ! Okay, finally add the DOF into the list
        call addDofToDirichletEntry(p_rdirichlet, idof, Dvalues(1), ndofs)
        
        ! And mark the DOF as 'processed' in the bitmap
        p_IdofBitmap(idofHigh) = ior(p_IdofBitmap(idofHigh),int(idofMask))
        
        ! Let us go for the next element adjacent to the vertex
      
      end do ! idx
      
      ! Go for the next vertex in the mesh region
    
    end do ! ivt
    
    ! Reset Iwhere
    Iwhere = 0
    
    if (p_rtria%ndim .eq. NDIM2D) then
      
      ! Go through all edges in the mesh region
      do i=1, iregionNMT
      
        ! Get the index of the edge
        imt = p_IedgeIdx(i)
        
        ! Store the edge number
        Iwhere(2) = imt
        
        ! And go through all elements which are adjacent to this edge
        do idx = 1, 2
        
          ! Get the index of the element
          iel = p_IelemAtEdge2D(idx, imt)
          if (iel .eq. 0) cycle
          
          ! Store the element number
          Iwhere(4) = iel
        
          ! Get all global dofs on this element
          call dof_locGlobMapping(p_rspatDisc, iel, IdofGlob)
          
          ! Now get the element type for this element (if the discretisation
          ! is not uniform)
          if (.not. buniform) ielemType = &
            p_rspatDisc%RelementDistr(p_IelemDist(iel))%celement
          
          ! Now we need to find which global DOF the currently processed
          ! edge belongs to.
          idof = 0
          select case(elem_getPrimaryElement(ielemType))
          ! 2D P1~ element
          case (EL_P1T)
            do j=1,3
              if ((p_IedgesAtElem(j,iel)) .eq. imt) then
                idof = IdofGlob(j)
                exit
              end if
            end do

          ! 2D P2 element
          case (EL_P2)
            do j=1,3
              if ((p_IedgesAtElem(j,iel)) .eq. imt) then
                ! The first 3 DOFs belong to the vertices
                idof = IdofGlob(3+j)
                exit
              end if
            end do

          ! 2D Q2 element
          case (EL_Q2)
            do j=1,4
              if ((p_IedgesAtElem(j,iel)) .eq. imt) then
                ! The first 4 DOFs belong to the vertices
                idof = IdofGlob(4+j)
                exit
              end if
            end do
            
          ! 2D Q1~ element
          case (EL_Q1T,EL_Q1TB)
            do j=1,4
              if ((p_IedgesAtElem(j,iel)) .eq. imt) then
                idof = IdofGlob(j)
                exit
              end if
            end do
            
          ! 2D Q2~ element
          case (EL_Q2T,EL_Q2TB)
            do j=1,4
            
              if ((p_IedgesAtElem(j,iel)) .ne. imt) cycle

              ! Get the first DOF on this edge
              idof = IdofGlob(j)

              ! Let us check if we have already processed this dof
              idofHigh = ishft(idof-1,-5) + 1
              idofMask = int(ishft(1,iand(idof-1,31)),I32)
              if (iand(p_IdofBitmap(idofHigh),int(idofMask)) .ne. 0) exit
              
              ! Okay, the DOF has not been set yet.
              ! So let us take care of the first Gauss point.
              Dwhere(1) = -G2P
              Dcoord2D(1:2) = Q12 * (&
                (1.0_DP+G2P)*p_DvertexCoords(1:2,p_IvertsAtEdge(1,imt))+&
                (1.0_DP-G2P)*p_DvertexCoords(1:2,p_IvertsAtEdge(2,imt)))
              call fgetBoundaryValuesMR(Icomponents, p_rspatDisc, rmeshRegion,&
                  cinfoNeeded, Iwhere, Dwhere(1:1), Dcoord2D, Dvalues, rcollection)
              daux1 = Dvalues(1)

              ! And take care of the second Gauss point.
              Dwhere(1) = G2P
              Dcoord2D(1:2) = Q12 * (&
                (1.0_DP-G2P)*p_DvertexCoords(1:2,p_IvertsAtEdge(1,imt))+&
                (1.0_DP+G2P)*p_DvertexCoords(1:2,p_IvertsAtEdge(2,imt)))
              call fgetBoundaryValuesMR(Icomponents, p_rspatDisc, rmeshRegion,&
                  cinfoNeeded, Iwhere, Dwhere(1:1), Dcoord2D, Dvalues, rcollection)
              daux2 = Dvalues(1)
              
              ! Add integral-mean DOF value
              call addDofToDirichletEntry(p_rdirichlet, idof, &
                  0.5_DP*(daux1+daux2), ndofs)
              p_IdofBitmap(idofHigh) = ior(p_IdofBitmap(idofHigh),int(idofMask))
              
              ! Calculate DOF of weighted integral-mean
              idof = IdofGlob(j+4)
              idofHigh = ishft(idof-1,-5) + 1
              idofMask = int(ishft(1,iand(idof-1,31)),I32)
              
              ! Add weighted integral-mean DOF value
              call addDofToDirichletEntry(p_rdirichlet, idof, &
                  0.5_DP*G2P*(-daux1+daux2), ndofs)
              p_IdofBitmap(idofHigh) = ior(p_IdofBitmap(idofHigh),int(idofMask))
              
              ! That is it for this edge
              exit
              
            end do
            
            ! Set idof to 0 to avoid that the code below is executed
            idof = 0
          
          ! 2D P3, Q3 and others are currently not supported
          
          end select

          ! If we come out here and idof is 0, then either the element
          ! does not have DOFs in the edges or something unforseen has
          ! happened...
          if (idof .eq. 0) cycle
          
          ! Let us check if we have already processed this dof
          idofHigh = ishft(idof-1,-5) + 1
          idofMask = int(ishft(1,iand(idof-1,31)),I32)
          if (iand(p_IdofBitmap(idofHigh),int(idofMask)) .ne. 0) cycle
          
          ! Calculate the coordinates of the edge midpoint
          Dcoord2D(1:2) = Q12 * (p_DvertexCoords(1:2, p_IvertsAtEdge(1,imt)) +&
            p_DvertexCoords(1:2, p_IvertsAtEdge(2,imt)))
            
          ! Dwhere is always 0
          Dwhere = 0.0_DP

          ! This is a new DOF for the list - so call the boundary values callback
          ! routine to calculate the value
          call fgetBoundaryValuesMR(Icomponents, p_rspatDisc, rmeshRegion,&
              cinfoNeeded, Iwhere, Dwhere(1:1), Dcoord2D, Dvalues, rcollection)
          
          ! Okay, finally add the DOF into the list
          call addDofToDirichletEntry(p_rdirichlet, idof, Dvalues(1), ndofs)
          
          ! And mark the DOF as 'processed' in the bitmap
          p_IdofBitmap(idofHigh) = ior(p_IdofBitmap(idofHigh),int(idofMask))
          
          ! Let us go for the next element adjacent to the edge
                
        end do ! idx
      
        ! Go for the next edge in the mesh region
      
      end do ! imt
    
    else if (p_rtria%ndim .eq. NDIM3D) then
      
      ! Go through all edges in the mesh region
      do i=1, iregionNMT
      
        ! Get the index of the edge
        imt = p_IedgeIdx(i)
        
        ! Store the edge number
        Iwhere(2) = imt
        
        ! And go through all elements which are adjacent to this edge
        do idx = p_IelemAtEdgeIdx(imt), p_IelemAtEdgeIdx(imt+1)-1
        
          ! Get the index of the element
          iel = p_IelemAtEdge(idx)
          
          ! Store the element number
          Iwhere(4) = iel
        
          ! Get all global dofs on this element
          call dof_locGlobMapping(p_rspatDisc, iel, IdofGlob)
          
          ! Now get the element type for this element (if the discretisation
          ! is not uniform)
          if (.not. buniform) ielemType = &
            p_rspatDisc%RelementDistr(p_IelemDist(iel))%celement
          
          ! Now we need to find which global DOF the currently processed
          ! edge belongs to.
          idof = 0
          select case(elem_getPrimaryElement(ielemType))
          ! 3D Q2 element
          case (EL_Q2_3D)
            do j=1,12
              if ((p_IedgesAtElem(j,iel)) .eq. imt) then
                ! The first 8 DOFs belong to the vertices
                idof = IdofGlob(8+j)
                exit
              end if
            end do
          
          end select

          ! If we come out here and idof is 0, then either the element
          ! does not have DOFs in the edges or something unforseen has
          ! happened...
          if (idof .eq. 0) cycle
          
          ! Let us check if we have already processed this dof
          idofHigh = ishft(idof-1,-5) + 1
          idofMask = int(ishft(1,iand(idof-1,31)),I32)
          if (iand(p_IdofBitmap(idofHigh),int(idofMask)) .ne. 0) cycle
          
          ! Calculate the coordinates of the edge midpoint
          Dcoord3D(1:3) = Q12 * (p_DvertexCoords(1:3, p_IvertsAtEdge(1,imt)) +&
            p_DvertexCoords(1:3, p_IvertsAtEdge(2,imt)))
            
          ! Dwhere is always 0
          Dwhere = 0.0_DP

          ! This is a new DOF for the list - so call the boundary values callback
          ! routine to calculate the value
          call fgetBoundaryValuesMR(Icomponents, p_rspatDisc, rmeshRegion,&
              cinfoNeeded, Iwhere, Dwhere(1:1), Dcoord3D, Dvalues, rcollection)
          
          ! Okay, finally add the DOF into the list
          call addDofToDirichletEntry(p_rdirichlet, idof, Dvalues(1), ndofs)
          
          ! And mark the DOF as 'processed' in the bitmap
          p_IdofBitmap(idofHigh) = ior(p_IdofBitmap(idofHigh),int(idofMask))
          
          ! Let us go for the next element adjacent to the edge
                
        end do ! idx
      
        ! Go for the next edge in the mesh region
      
      end do ! imt
    
    end if
    
    ! Reset Iwhere
    Iwhere = 0
    
    ! Go through all faces in the mesh region
    do i=1, iregionNAT
    
      ! Get the index of the face
      iat = p_IfaceIdx(i)
      
      ! Store the face number
      Iwhere(3) = iat
      
      ! And go through all elements which are adjacent to this face
      do idx=1,2
      
        ! Get the index of the element
        iel = p_IelemAtFace(idx,iat)
        if (iel .eq. 0) cycle
        
        ! Store the element number
        Iwhere(4) = iel
      
        ! Get all global dofs on this element
        call dof_locGlobMapping(p_rspatDisc, iel, IdofGlob)
        
        ! Now get the element type for this element (if the discretisation
        ! is not uniform)
        if (.not. buniform) ielemType = &
          p_rspatDisc%RelementDistr(p_IelemDist(iel))%celement
        
        ! Now we need to find which global DOF the currently processed
        ! face belongs to.
        idof = 0
        select case(elem_getPrimaryElement(ielemType))
        ! 3D Q2 element
        case (EL_Q2_3D)
          do j=1,6
            if ((p_IfacesAtElem(j,iel)) .eq. iat) then
              idof = IdofGlob(j+20)
              exit
            end if
          end do

        ! 3D Q1~ element
        case (EL_Q1T_3D)
          do j=1,6
            if ((p_IfacesAtElem(j,iel)) .eq. iat) then
              idof = IdofGlob(j)
              exit
            end if
          end do
        
        ! 3D Q2~ element
        case (EL_Q2T_3D)
          do j = 1, 6
          
            if ((p_IfacesAtElem(j,iel)) .ne. iat) cycle

            ! Get the first DOF on this face
            idof = IdofGlob(j)

            ! Let us check if we have already processed this dof
            idofHigh = ishft(idof-1,-5) + 1
            idofMask = int(ishft(1,iand(idof-1,31)),I32)
            if (iand(p_IdofBitmap(idofHigh),int(idofMask)) .ne. 0) exit

            ! Okay, the DOF has not been set yet.
            ! So let us take care of the first Gauss point.
            Dwhere(1:2) = -G2P
            Dcoord3D(1:3) = Q14 * (&
              (1.0_DP+G2P)*(1.0_DP+G2P)*p_DvertexCoords(1:3,p_IvertsAtFace(1,iat)) +&
              (1.0_DP-G2P)*(1.0_DP+G2P)*p_DvertexCoords(1:3,p_IvertsAtFace(2,iat)) +&
              (1.0_DP-G2P)*(1.0_DP-G2P)*p_DvertexCoords(1:3,p_IvertsAtFace(3,iat)) +&
              (1.0_DP+G2P)*(1.0_DP-G2P)*p_DvertexCoords(1:3,p_IvertsAtFace(4,iat)))
            call fgetBoundaryValuesMR(Icomponents, p_rspatDisc, rmeshRegion,&
                cinfoNeeded, Iwhere, Dwhere(1:2), Dcoord3D, Dvalues, rcollection)
            daux1 = Dvalues(1)

            ! Second Gauss point.
            Dwhere(1) =  G2P
            Dwhere(2) = -G2P
            Dcoord3D(1:3) = Q14 * (&
              (1.0_DP-G2P)*(1.0_DP+G2P)*p_DvertexCoords(1:3,p_IvertsAtFace(1,iat)) +&
              (1.0_DP+G2P)*(1.0_DP+G2P)*p_DvertexCoords(1:3,p_IvertsAtFace(2,iat)) +&
              (1.0_DP+G2P)*(1.0_DP-G2P)*p_DvertexCoords(1:3,p_IvertsAtFace(3,iat)) +&
              (1.0_DP-G2P)*(1.0_DP-G2P)*p_DvertexCoords(1:3,p_IvertsAtFace(4,iat)))
            call fgetBoundaryValuesMR(Icomponents, p_rspatDisc, rmeshRegion,&
                cinfoNeeded, Iwhere, Dwhere(1:2), Dcoord3D, Dvalues, rcollection)
            daux2 = Dvalues(1)

            ! Third Gauss point.
            Dwhere(1) =  G2P
            Dwhere(2) =  G2P
            Dcoord3D(1:3) = Q14 * (&
              (1.0_DP-G2P)*(1.0_DP-G2P)*p_DvertexCoords(1:3,p_IvertsAtFace(1,iat)) +&
              (1.0_DP+G2P)*(1.0_DP-G2P)*p_DvertexCoords(1:3,p_IvertsAtFace(2,iat)) +&
              (1.0_DP+G2P)*(1.0_DP+G2P)*p_DvertexCoords(1:3,p_IvertsAtFace(3,iat)) +&
              (1.0_DP-G2P)*(1.0_DP+G2P)*p_DvertexCoords(1:3,p_IvertsAtFace(4,iat)))
            call fgetBoundaryValuesMR(Icomponents, p_rspatDisc, rmeshRegion,&
                cinfoNeeded, Iwhere, Dwhere(1:2), Dcoord3D, Dvalues, rcollection)
            daux3 = Dvalues(1)

            ! Fourth Gauss point.
            Dwhere(1) = -G2P
            Dwhere(2) =  G2P
            Dcoord3D(1:3) = Q14 * (&
              (1.0_DP+G2P)*(1.0_DP-G2P)*p_DvertexCoords(1:3,p_IvertsAtFace(1,iat)) +&
              (1.0_DP-G2P)*(1.0_DP-G2P)*p_DvertexCoords(1:3,p_IvertsAtFace(2,iat)) +&
              (1.0_DP-G2P)*(1.0_DP+G2P)*p_DvertexCoords(1:3,p_IvertsAtFace(3,iat)) +&
              (1.0_DP+G2P)*(1.0_DP+G2P)*p_DvertexCoords(1:3,p_IvertsAtFace(4,iat)))
            call fgetBoundaryValuesMR(Icomponents, p_rspatDisc, rmeshRegion,&
                cinfoNeeded, Iwhere, Dwhere(1:2), Dcoord3D, Dvalues, rcollection)
            daux4 = Dvalues(1)

            ! Add integral-mean DOF value
            call addDofToDirichletEntry(p_rdirichlet, idof, &
                0.25_DP*(daux1+daux2+daux3+daux4), ndofs)
            p_IdofBitmap(idofHigh) = ior(p_IdofBitmap(idofHigh),int(idofMask))

            ! Calculate DOF of first weighted integral-mean
            idof = IdofGlob(6+2*(j-1)+1)
            idofHigh = ishft(idof-1,-5) + 1
            idofMask = int(ishft(1,iand(idof-1,31)),I32)
              
            ! Add first weighted integral-mean DOF value
            call addDofToDirichletEntry(p_rdirichlet, idof, &
               0.25_DP*G2P*(-daux1+daux2+daux3-daux4), ndofs)
            p_IdofBitmap(idofHigh) = ior(p_IdofBitmap(idofHigh),int(idofMask))

            ! Calculate DOF of second weighted integral-mean
            idof = IdofGlob(6+2*(j-1)+2)
            idofHigh = ishft(idof-1,-5) + 1
            idofMask = int(ishft(1,iand(idof-1,31)),I32)
              
            ! Add second weighted integral-mean DOF value
            call addDofToDirichletEntry(p_rdirichlet, idof, &
               0.25_DP*G2P*(-daux1-daux2+daux3+daux4), ndofs)
            p_IdofBitmap(idofHigh) = ior(p_IdofBitmap(idofHigh),int(idofMask))
            
            ! That is it for this face
            exit
          
          end do

          ! Set idof to 0 to avoid that the code below is executed
          idof = 0
        
        end select

        ! If we come out here and idof is 0, then either the element
        ! does not have DOFs in the faces or something unforseen has
        ! happened...
        if (idof .eq. 0) cycle
        
        ! Let us check if we have already processed this dof
        idofHigh = ishft(idof-1,-5) + 1
        idofMask = int(ishft(1,iand(idof-1,31)),I32)
        if (iand(p_IdofBitmap(idofHigh),int(idofMask)) .ne. 0) cycle
        
        ! Calculate the coordinates of the face midpoint
        Dcoord3D(1:3) = Q14 * (p_DvertexCoords(1:3, p_IvertsAtFace(1,iat)) +&
          p_DvertexCoords(1:3, p_IvertsAtFace(2,iat)) +&
          p_DvertexCoords(1:3, p_IvertsAtFace(3,iat)) +&
          p_DvertexCoords(1:3, p_IvertsAtFace(4,iat)))
        
        ! Dwhere is always (0, 0)
        Dwhere = 0.0_DP

        ! This is a new DOF for the list - so call the boundary values callback
        ! routine to calculate the value
        call fgetBoundaryValuesMR(Icomponents, p_rspatDisc, rmeshRegion,&
            cinfoNeeded, Iwhere, Dwhere(1:2), Dcoord3D, Dvalues, rcollection)
        
        ! Okay, finally add the DOF into the list
        call addDofToDirichletEntry(p_rdirichlet, idof, Dvalues(1), ndofs)
        
        ! And mark the DOF as 'processed' in the bitmap
        p_IdofBitmap(idofHigh) = ior(p_IdofBitmap(idofHigh),int(idofMask))
        
        ! Let us go for the next element adjacent to the face

      end do ! idx
      
      ! Go for the next face in the mesh region
    
    end do
    
    ! Reset Iwhere
    Iwhere = 0
    
    ! Go through all elements in the mesh region
    do i=1, iregionNEL
    
      ! Get the element index
      iel = p_IelementIdx(i)
      
      ! Store the element number
      Iwhere(4) = iel

      ! Get all global dofs on this element
      call dof_locGlobMapping(p_rspatDisc, iel, IdofGlob)
      
      ! Now get the element type for this element (if the discretisation
      ! is not uniform)
      if (.not. buniform) ielemType = &
        p_rspatDisc%RelementDistr(p_IelemDist(iel))%celement
      
      ! Now we need to find which global DOF the currently processed
      ! element belongs to.
      idof = 0
      select case(elem_getPrimaryElement(ielemType))
      ! 1D P0 element
      case (EL_P0_1D)
        ! The line-midpoint DOF has number 1
        idof = IdofGlob(1)
        
        ! Calculate line-midpoint
        Dcoord1D(1) = Q12 * (p_DvertexCoords(1,p_IvertsAtElem(1,iel)) +&
          p_DvertexCoords(1,p_IvertsAtElem(2,iel)))
        
        ! Dwhere is 0
        Dwhere = 0.0_DP
        
        ! Call the boundary condition callback routine
        call fgetBoundaryValuesMR(Icomponents, p_rspatDisc, rmeshRegion,&
            cinfoNeeded, Iwhere, Dwhere(1:1), Dcoord1D, Dvalues, rcollection)

        ! Add the DOF into the list
        call addDofToDirichletEntry(p_rdirichlet, idof, Dvalues(1), ndofs)

      ! 1D P2 element
      case (EL_P2_1D)
        ! The line-midpoint DOF has number 3
        idof = IdofGlob(3)
        
        ! Calculate line-midpoint
        Dcoord1D(1) = Q12 * (p_DvertexCoords(1,p_IvertsAtElem(1,iel)) +&
          p_DvertexCoords(1,p_IvertsAtElem(2,iel)))

        ! Dwhere is 0
        Dwhere = 0.0_DP
        
        ! Call the boundary condition callback routine
        call fgetBoundaryValuesMR(Icomponents, p_rspatDisc, rmeshRegion,&
            cinfoNeeded, Iwhere, Dwhere(1:1), Dcoord1D, Dvalues, rcollection)

        ! Add the DOF into the list
        call addDofToDirichletEntry(p_rdirichlet, idof, Dvalues(1), ndofs)
      
      ! 2D P0 element
      case (EL_P0)
        ! The triangle-midpoint DOF has number 1
        idof = IdofGlob(1)
        
        ! Calculate triangle-midpoint
        Dcoord2D(1:2) = Q13 * (p_DvertexCoords(1:2,p_IvertsAtElem(1,iel)) +&
          p_DvertexCoords(1:2,p_IvertsAtElem(2,iel)) +&
          p_DvertexCoords(1:2,p_IvertsAtElem(3,iel)))
        
        ! Dwhere is (1/3, 1/3)
        Dwhere = Q13
        
        ! Call the boundary condition callback routine
        call fgetBoundaryValuesMR(Icomponents, p_rspatDisc, rmeshRegion,&
            cinfoNeeded, Iwhere, Dwhere(1:2), Dcoord2D, Dvalues, rcollection)
      
        ! Add the DOF into the list
        call addDofToDirichletEntry(p_rdirichlet, idof, Dvalues(1), ndofs)

      ! 2D P2 element
      case (EL_P2)
        ! The triangle-midpoint DOF has number 7
        idof = IdofGlob(7)
        
        ! Calculate triangle-midpoint
        Dcoord2D(1:2) = Q13 * (p_DvertexCoords(1:2,p_IvertsAtElem(1,iel)) +&
          p_DvertexCoords(1:2,p_IvertsAtElem(2,iel)) +&
          p_DvertexCoords(1:2,p_IvertsAtElem(3,iel)))
        
        ! Dwhere is (1/3, 1/3)
        Dwhere = Q13

        ! Call the boundary condition callback routine
        call fgetBoundaryValuesMR(Icomponents, p_rspatDisc, rmeshRegion,&
            cinfoNeeded, Iwhere, Dwhere(1:2),Dcoord2D, Dvalues, rcollection)
      
        ! Add the DOF into the list
        call addDofToDirichletEntry(p_rdirichlet, idof, Dvalues(1), ndofs)

      ! 2D Q0/QP1 element
      case (EL_Q0,EL_QP1)
        ! The quad-midpoint DOF has number 1
        idof = IdofGlob(1)
        
        ! Calculate quad-midpoint
        Dcoord2D(1:2) = Q14 * (p_DvertexCoords(1:2,p_IvertsAtElem(1,iel)) +&
          p_DvertexCoords(1:2,p_IvertsAtElem(2,iel)) +&
          p_DvertexCoords(1:2,p_IvertsAtElem(3,iel)) +&
          p_DvertexCoords(1:2,p_IvertsAtElem(4,iel)))
        
        ! Dwhere is (0,0)
        Dwhere = 0.0_DP
        
        ! Call the boundary condition callback routine
        call fgetBoundaryValuesMR(Icomponents, p_rspatDisc, rmeshRegion,&
            cinfoNeeded, Iwhere, Dwhere(1:2), Dcoord2D, Dvalues, rcollection)

        ! Add the DOF into the list
        call addDofToDirichletEntry(p_rdirichlet, idof, Dvalues(1), ndofs)

      ! 2D Q2 element
      case (EL_Q2)
        ! The quad-midpoint DOF has number 9
        idof = IdofGlob(9)
        
        ! Calculate quad-midpoint
        Dcoord2D(1:2) = Q14 * (p_DvertexCoords(1:2,p_IvertsAtElem(1,iel)) +&
          p_DvertexCoords(1:2,p_IvertsAtElem(2,iel)) +&
          p_DvertexCoords(1:2,p_IvertsAtElem(3,iel)) +&
          p_DvertexCoords(1:2,p_IvertsAtElem(4,iel)))
        
        ! Dwhere is (0,0)
        Dwhere = 0.0_DP
        
        ! Call the boundary condition callback routine
        call fgetBoundaryValuesMR(Icomponents, p_rspatDisc, rmeshRegion,&
            cinfoNeeded, Iwhere, Dwhere(1:2), Dcoord2D, Dvalues, rcollection)

        ! Add the DOF into the list
        call addDofToDirichletEntry(p_rdirichlet, idof, Dvalues(1), ndofs)
        
      ! 2D Q1TB element
      case (EL_Q1TB)
        ! The quad-integral-mean DOF has number 5
        idof = IdofGlob(5)
        
        ! Calculate quad-midpoint to approximate the integral mean
        Dcoord2D(1:2) = Q14 * (p_DvertexCoords(1:2,p_IvertsAtElem(1,iel)) +&
          p_DvertexCoords(1:2,p_IvertsAtElem(2,iel)) +&
          p_DvertexCoords(1:2,p_IvertsAtElem(3,iel)) +&
          p_DvertexCoords(1:2,p_IvertsAtElem(4,iel)))
        
        ! Dwhere is (0,0)
        Dwhere = 0.0_DP
        
        ! Call the boundary condition callback routine
        call fgetBoundaryValuesMR(Icomponents, p_rspatDisc, rmeshRegion,&
            cinfoNeeded, Iwhere, Dwhere(1:2), Dcoord2D, Dvalues, rcollection)

        ! Add the DOF into the list
        call addDofToDirichletEntry(p_rdirichlet, idof, Dvalues(1), ndofs)
        
      ! 2D Q2T element
      case (EL_Q2T)
        ! The quad-integral-mean DOF has number 9
        idof = IdofGlob(9)
        
        ! Gauss-Point 1
        Dwhere(1) = -G2P
        Dwhere(2) = -G2P
        Dcoord2D(1:2) = Q14 * (&
          (1.0_DP+G2P)*(1.0_DP+G2P)*p_DvertexCoords(1:2,p_IvertsAtElem(1,iel)) +&
          (1.0_DP-G2P)*(1.0_DP+G2P)*p_DvertexCoords(1:2,p_IvertsAtElem(2,iel)) +&
          (1.0_DP-G2P)*(1.0_DP-G2P)*p_DvertexCoords(1:2,p_IvertsAtElem(3,iel)) +&
          (1.0_DP+G2P)*(1.0_DP-G2P)*p_DvertexCoords(1:2,p_IvertsAtElem(4,iel)))
        call fgetBoundaryValuesMR(Icomponents, p_rspatDisc, rmeshRegion,&
            cinfoNeeded, Iwhere, Dwhere(1:2), Dcoord2D, Dvalues, rcollection)
        daux1 = Dvalues(1)

        ! Gauss-Point 2
        Dwhere(1) =  G2P
        Dwhere(2) = -G2P
        Dcoord2D(1:2) = Q14 * (&
          (1.0_DP-G2P)*(1.0_DP+G2P)*p_DvertexCoords(1:2,p_IvertsAtElem(1,iel)) +&
          (1.0_DP+G2P)*(1.0_DP+G2P)*p_DvertexCoords(1:2,p_IvertsAtElem(2,iel)) +&
          (1.0_DP+G2P)*(1.0_DP-G2P)*p_DvertexCoords(1:2,p_IvertsAtElem(3,iel)) +&
          (1.0_DP-G2P)*(1.0_DP-G2P)*p_DvertexCoords(1:2,p_IvertsAtElem(4,iel)))
        call fgetBoundaryValuesMR(Icomponents, p_rspatDisc, rmeshRegion,&
            cinfoNeeded, Iwhere, Dwhere(1:2), Dcoord2D, Dvalues, rcollection)
        daux2 = Dvalues(1)

        ! Gauss-Point 3
        Dwhere(1) =  G2P
        Dwhere(2) =  G2P
        Dcoord2D(1:2) = Q14 * (&
          (1.0_DP-G2P)*(1.0_DP-G2P)*p_DvertexCoords(1:2,p_IvertsAtElem(1,iel)) +&
          (1.0_DP+G2P)*(1.0_DP-G2P)*p_DvertexCoords(1:2,p_IvertsAtElem(2,iel)) +&
          (1.0_DP+G2P)*(1.0_DP+G2P)*p_DvertexCoords(1:2,p_IvertsAtElem(3,iel)) +&
          (1.0_DP-G2P)*(1.0_DP+G2P)*p_DvertexCoords(1:2,p_IvertsAtElem(4,iel)))
        call fgetBoundaryValuesMR(Icomponents, p_rspatDisc, rmeshRegion,&
            cinfoNeeded, Iwhere, Dwhere(1:2), Dcoord2D, Dvalues, rcollection)
        daux3 = Dvalues(1)

        ! Gauss-Point 4
        Dwhere(1) = -G2P
        Dwhere(2) =  G2P
        Dcoord2D(1:2) = Q14 * (&
          (1.0_DP+G2P)*(1.0_DP-G2P)*p_DvertexCoords(1:2,p_IvertsAtElem(1,iel)) +&
          (1.0_DP-G2P)*(1.0_DP-G2P)*p_DvertexCoords(1:2,p_IvertsAtElem(2,iel)) +&
          (1.0_DP-G2P)*(1.0_DP+G2P)*p_DvertexCoords(1:2,p_IvertsAtElem(3,iel)) +&
          (1.0_DP+G2P)*(1.0_DP+G2P)*p_DvertexCoords(1:2,p_IvertsAtElem(4,iel)))
        call fgetBoundaryValuesMR(Icomponents, p_rspatDisc, rmeshRegion,&
            cinfoNeeded, Iwhere, Dwhere(1:2), Dcoord2D, Dvalues, rcollection)
        daux4 = Dvalues(1)

        ! Add the integral-mean DOF into the list
        call addDofToDirichletEntry(p_rdirichlet, idof, &
            0.25_DP*(daux1+daux2+daux3+daux4), ndofs)
        
      ! 3D Q0 element
      case (EL_Q0_3D)
        ! The hexa-midpoint DOF has number 1
        idof = IdofGlob(1)
        
        ! Calculate hexa-midpoint
        Dcoord3D(1:2) = Q18 * (p_DvertexCoords(1:2,p_IvertsAtElem(1,iel)) +&
          p_DvertexCoords(1:2,p_IvertsAtElem(2,iel)) +&
          p_DvertexCoords(1:2,p_IvertsAtElem(3,iel)) +&
          p_DvertexCoords(1:2,p_IvertsAtElem(4,iel)) +&
          p_DvertexCoords(1:2,p_IvertsAtElem(5,iel)) +&
          p_DvertexCoords(1:2,p_IvertsAtElem(6,iel)) +&
          p_DvertexCoords(1:2,p_IvertsAtElem(7,iel)) +&
          p_DvertexCoords(1:2,p_IvertsAtElem(8,iel)))
        
        ! Dwhere is (0,0,0)
        Dwhere = 0.0_DP
        
        ! Call the boundary condition callback routine
        call fgetBoundaryValuesMR(Icomponents, p_rspatDisc, rmeshRegion,&
            cinfoNeeded, Iwhere, Dwhere, Dcoord3D, Dvalues, rcollection)

        ! Add the DOF into the list
        call addDofToDirichletEntry(p_rdirichlet, idof, Dvalues(1), ndofs)

      ! 3D Q2 element
      case (EL_Q2_3D)
        ! The quad-midpoint DOF has number 27
        idof = IdofGlob(27)
        
        ! Calculate hexa-midpoint
        Dcoord3D(1:2) = Q18 * (p_DvertexCoords(1:2,p_IvertsAtElem(1,iel)) +&
          p_DvertexCoords(1:2,p_IvertsAtElem(2,iel)) +&
          p_DvertexCoords(1:2,p_IvertsAtElem(3,iel)) +&
          p_DvertexCoords(1:2,p_IvertsAtElem(4,iel)) +&
          p_DvertexCoords(1:2,p_IvertsAtElem(5,iel)) +&
          p_DvertexCoords(1:2,p_IvertsAtElem(6,iel)) +&
          p_DvertexCoords(1:2,p_IvertsAtElem(7,iel)) +&
          p_DvertexCoords(1:2,p_IvertsAtElem(8,iel)))
        
        ! Dwhere is (0,0,0)
        Dwhere = 0.0_DP
        
        ! Call the boundary condition callback routine
        call fgetBoundaryValuesMR(Icomponents, p_rspatDisc, rmeshRegion,&
            cinfoNeeded, Iwhere, Dwhere, Dcoord3D, Dvalues, rcollection)

        ! Add the DOF into the list
        call addDofToDirichletEntry(p_rdirichlet, idof, Dvalues(1), ndofs)
      
      end select

      ! Let us go for the next element in the mesh region
      
    end do
   
    ! Deallocate the DOF map
    deallocate(p_IdofBitmap)
    
    ! That is it - unbelievable!
    
  contains
  
    ! This auxiliary routine adds a dof and a dirichlet value into a
    ! discrete dirichlet boundary condition structure.
    subroutine addDofToDirichletEntry(rdirichlet, idof, dvalue, ndofs)
    
    ! The dirichlet boundary condition structure
    type(t_discreteBCDirichlet), intent(inout) :: rdirichlet
    
    ! The number of the dof that is to be added
    integer, intent(in) :: idof
    
    ! The dirichlet value of the dof
    real(DP), intent(in) :: dvalue
    
    ! If the dirichlet entry structure is already initialized, then this
    ! corresponds to the total number of currently allocated dofs in the
    ! dirichlet entry structure.
    ! If the structure is uninitialized, then ndofs specifies the initial
    ! size which is to be used for the arrays.
    integer, intent(inout) :: ndofs
    
    ! Some local variables
    integer, dimension(:), pointer :: p_Idofs
    real(DP), dimension(:), pointer :: p_Dvalues
    
    ! Number of entries added if the list is full
    integer, parameter :: NINC = 64
    
      ! If the list is empty, then allocate it
      if (rdirichlet%nDOF .eq. 0) then
      
        ! Make sure that we do not allocate empty arrays
        if (ndofs .le. 0) ndofs = NINC
      
        ! allocate the arrays
        call storage_new('addDofToDirichletEntry', 'p_IdirichletDOFs',&
          ndofs, ST_INT, rdirichlet%h_IdirichletDOFs, ST_NEWBLOCK_ZERO)
        call storage_new('addDofToDirichletEntry', 'p_DdirichletValues',&
          ndofs, ST_DOUBLE, rdirichlet%h_DdirichletValues, ST_NEWBLOCK_ZERO)
      
      else if (rdirichlet%nDOF .ge. ndofs) then
        
        ! The list is full, so resize it
        call storage_realloc('addDofToDirichletEntry', ndofs+NINC,&
          rdirichlet%h_IdirichletDOFs, ST_NEWBLOCK_ZERO, .true.)
        call storage_realloc('addDofToDirichletEntry', ndofs+NINC,&
          rdirichlet%h_DdirichletValues, ST_NEWBLOCK_ZERO, .true.)
        
        ! And remember that we have increased the list size
        ndofs = ndofs + NINC
      
      end if

      ! Get the lists
      call storage_getbase_int(rdirichlet%h_IdirichletDOFs, p_Idofs)
      call storage_getbase_double(rdirichlet%h_DdirichletValues, p_Dvalues)
      
      ! Add the new dof to the lists
      rdirichlet%nDOF = rdirichlet%nDOF + 1
      p_Idofs(rdirichlet%nDOF) = idof
      p_Dvalues(rdirichlet%nDOF) = dvalue
          
    end subroutine ! addDofToDirichletEntry(...)

  end subroutine
  
end module
