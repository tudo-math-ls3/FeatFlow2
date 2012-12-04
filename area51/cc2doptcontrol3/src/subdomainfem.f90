!##############################################################################
!# ****************************************************************************
!# <name> subdomainfem </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module provides the functionality to work with the FEM method
!# on subdomains of a mesh. A subdomain including an appropriate
!# discretisation can be created from a discretisation structure.
!# It can be used to create vectors which represent FEM functions on
!# these subdomains. Adapted linear algebra routines provide the possibility
!# to extract the FEM values from the original mesh and write it to a 
!# subvector and vice versa.
!#
!# Routines in this module:
!#
!# 1.) sdfem_createSubdFem
!#     -> creates a subdomain-FEM structure from a list of elements or
!#        a boundary refion
!#
!# 2.) sdfem_releaseSubdFem
!#     -> releases a subdomain-FEM structure
!#
!# 3.) sdfem_extractData
!#     -> Extracts the DOFs of the subdomain from a full-domain vector
!#
!# 4.) sdfem_linearCombToGlobal / sdfem_linearCombToLocal
!#     -> linear combination of a subdomain vector into a full-domain vector
!#        or vice versa
!#
!# 5.) sdfem_createSubdBlkFem
!#     -> creates a subdomain-FEM structure for block discretisations
!#
!# 6.) sdfem_releaseSubdBlkFem
!#     -> Releases the structure created by sdfem_createSubdBlkFem
!#
!# 7.) sdfem_subdBlkFemNewBlock
!#     -> Add a discretisation to a subdomain-FEM structure for block 
!#        discretisations
!#
!# 8.) sdfem_createBlkHierarchy
!#     -> Creates a hierarchy of subdomain-block FEM structures.
!#
!# 9.) sdfem_releaseBlkHierarchy
!#     -> Releases a hierarchy created by sdfem_createBlkHierarchy.
!#
!# </purpose>
!##############################################################################

module subdomainfem

  use fsystem
  use storage
  use genoutput
  
  use boundary
  use triangulation
  use element
  use spatialdiscretisation
  use dofmapping
  use bcassemblybase
  use linearalgebra
  use linearsystemscalar
  
  use meshhierarchy
  
  implicit none
  
  private
  
!<constants>

!<constantblock description = "Constants defining shared information">

  ! Triangulation is shared.
  integer(I32), parameter, public :: SDFEM_SHAREDTRIA = 2_I32**0

!</constantblock>

!</constants>

!<types>

  !<typeblock>
  
  ! A FE-space structure that represents a FE space on a subdomain.
  type t_sdfemFeSpace
  
    ! A shared flag that specifies which information in this structure
    ! is shared with information from outside and which information
    ! belongs to this structure.
    integer(I32) :: cflags = 0

    ! Reference to the underlying domain or NULL() if no domain is attached.
    type(t_boundary), pointer :: p_rboundary => null()

    ! Reference to the underlying discretisation of the full domain.
    type(t_spatialDiscretisation), pointer :: p_rdiscrGlobal => null()
    
    ! Reference to the underlying discretisation of the subdomain.
    type(t_spatialDiscretisation), pointer :: p_rdiscrSubdomain => null()
    
    ! Reference to the underlying triangulation that represents the subdomain.
    type(t_triangulation), pointer :: p_rtriaSubdomain => null()
    
    ! Number of DOFs in this subdomain
    integer :: ndof = 0
    
    ! Subdomain-to-domain DOF-mapping. This array defines for every DOF
    ! in the local domain the corresponding DOF on the global domain.
    integer :: h_IglobalDofs = ST_NOHANDLE
  
  end type
  
  !</typeblock>
  
  public :: t_sdfemFeSpace

  !<typeblock>
  
  ! A FE-space structure that represents a set of FE spaces on a subdomain.
  type t_sdfemBlkFeSpace
  
    ! Reference to the underlying discretisation of the full domain.
    ! this may point to NULL if no discretisation structure for the full
    ! domain is available.
    type(t_blockDiscretisation), pointer :: p_rdiscrGlobal => null()
    
    ! Number of blocks
    integer :: nblocks = 0
    
    ! Array of nblocks t_sdfemFeSpace structures for each block
    type(t_sdfemFeSpace), dimension(:), pointer :: p_RsdfemSpaces => null()
    
    ! Handle to the element list of all the handles defining the subdomain.
    integer :: h_Ielements = ST_NOHANDLE
  
    ! Reference to the underlying triangulation that represents the subdomain.
    type(t_triangulation), pointer :: p_rtriaSubdomain => null()
  
    ! Reference to the underlying triangulation of the domain
    type(t_triangulation), pointer :: p_rtriaDomain => null()

  end type
  
  !</typeblock>
  
  public :: t_sdfemBlkFeSpace

  !<typeblock>
  
  ! A hierarchy of FE-space block structures.
  type t_sdfemBlkFeHierarchy
  
    ! An underlying mesh hierarchy.
    type(t_meshHierarchy) :: rmeshHierarchy

    ! Number of levels available in this structure.
    integer :: nlevels = 0

    ! Level information.
    type(t_sdfemBlkFeSpace), dimension(:), pointer :: p_RfeBlkSpaces => null()


  end type
  
  !</typeblock>
  
  public :: t_sdfemBlkFeHierarchy

!</types>
  
  public :: sdfem_createSubdFem
  public :: sdfem_releaseSubdFem
  public :: sdfem_extractData
  public :: sdfem_linearCombToGlobal
  public :: sdfem_linearCombToLocal
  
  public :: sdfem_createSubdBlkFem
  public :: sdfem_releaseSubdBlkFem
  public :: sdfem_subdBlkFemNewBlock
  public :: spdiscr_createCompDiscrManif1D
  
  public :: sdfem_createBlkHierarchy
  public :: sdfem_releaseBlkHierarchy
  
  interface sdfem_createSubdFem
    module procedure sdfem_createSubdFemFromElList
    module procedure sdfem_createSubdFemFrom2DBdReg
  end interface

  interface sdfem_createSubdBlkFem
    module procedure sdfem_createSubdBlkFemElList
    module procedure sdfem_createSubdBlkFem2DBdReg
  end interface
  
  interface sdfem_createBlkHierarchy
    module procedure sdfem_createBlkHierarchyMeshH
  end interface
  
contains

  ! ***************************************************************************

!<subrotine>
  
  subroutine sdfem_createSubdFemFromElList(rsdfemFeSpace,rdiscretisation,Ielements,&
      rtriaSubdomain)

!<description>
  ! Creates a subdomain-FEM-structure for a set of elements. Ielements is a
  !´list of elements in the triangulation connected to rdiscretisation.
  ! The routine extracts these elements, creates a new triangulation and
  ! discretisation and saves all this in rsdfemFeSpace.
!</description>

!<input>
  ! Underlying discretisation of the full domain.
  type(t_spatialDiscretisation), intent(in), target :: rdiscretisation
  
  ! List of elements that form the subdomain.
  integer, dimension(:), intent(in) :: Ielements
  
  ! OPTIONAL: Triangulation of the subdomain.
  type(t_triangulation), target, optional :: rtriaSubdomain
!</input>

!<output>
  ! Structure that represents the FEM space on the subdomain.
  type(t_sdfemFeSpace), intent(out) :: rsdfemFeSpace
!</output>

!</subroutine>
    
    ! local variables
    integer :: i,ieldist,ielold,ielnew,ndof,idof
    type(t_spatialDiscretisation), pointer :: p_rdiscrSubdomain
    integer, dimension(:), pointer :: p_IelementDistrOld,p_IelementCounterOld
    integer, dimension(:), pointer :: p_IelementDistrNew,p_IelementCounterNew
    integer, dimension(:), pointer :: p_IelementListNew
    integer, dimension(:), pointer :: p_IglobalDofs
    integer, dimension(:), allocatable :: p_IdofsOld,p_IdofsNew

    if (size(Ielements) .eq. 0) return

    ! Save references to the new structure and create nonexisting
    ! substructures.
    rsdfemFeSpace%p_rdiscrGlobal => rdiscretisation
    rsdfemFeSpace%p_rboundary => rdiscretisation%p_rboundary
    allocate(rsdfemFeSpace%p_rdiscrSubdomain)
    
    ! Create a triangulation for the subdomain.
    if (present(rtriaSubdomain)) then
      rsdfemFeSpace%p_rtriaSubdomain => rtriaSubdomain
      
      ! Prevent the structure from being deallocated
      rsdfemFeSpace%cflags = ior (rsdfemFeSpace%cflags,SDFEM_SHAREDTRIA)
    else
      ! Create the subdomain triangulation
      allocate(rsdfemFeSpace%p_rtriaSubdomain)
      if (associated(rsdfemFeSpace%p_rboundary)) then
        call tria_generateSubdomain(rdiscretisation%p_rtriangulation, Ielements, &
            rsdfemFeSpace%p_rtriaSubdomain, rsdfemFeSpace%p_rboundary)
      else
        call tria_generateSubdomain(rdiscretisation%p_rtriangulation, Ielements, &
            rsdfemFeSpace%p_rtriaSubdomain)
      end if
    end if

    ! -----------------------------------------------------
    ! Duplicate the discretisation, set up the new one
    ! -----------------------------------------------------
    call spdiscr_duplicateDiscrSc2 (&
        rsdfemFeSpace%p_rdiscrGlobal, rsdfemFeSpace%p_rdiscrSubdomain, .false.)
        
    ! Rebuild the element sets
    p_rdiscrSubdomain => rsdfemFeSpace%p_rdiscrSubdomain
    p_rdiscrSubdomain%RelementDistr(:)%NEL = 0
    p_rdiscrSubdomain%p_rtriangulation => rsdfemFeSpace%p_rtriaSubdomain
    
    if (p_rdiscrSubdomain%h_IelementDistr .ne. ST_NOHANDLE) then
      call storage_getbase_int (p_rdiscrSubdomain%h_IelementDistr,p_IelementDistrOld)
      call storage_getbase_int (p_rdiscrSubdomain%h_IelementDistr,p_IelementDistrNew)
    else
      ! Uniform discretisation
      nullify(p_IelementDistrOld)
      nullify(p_IelementDistrNew)
    end if
    
    if (p_rdiscrSubdomain%h_IelementCounter .ne. ST_NOHANDLE) then
      call storage_getbase_int (p_rdiscrSubdomain%h_IelementCounter,p_IelementCounterOld)
      call storage_getbase_int (p_rdiscrSubdomain%h_IelementCounter,p_IelementCounterNew)
    else
      nullify(p_IelementCounterOld)
      nullify(p_IelementCounterNew)
    end if

    ! Loop over the element distributions.
    do ieldist = 1,p_rdiscrSubdomain%inumFESpaces
    
      ! On every element distribution, loop over all elements and
      ! insert the element if it belongs to this distribution.
      ! this will be faster than finding out the corresponding element
      ! distribution and obtaining a pointer to the corresponding list
      if (p_rdiscrSubdomain%RelementDistr(ieldist)%h_IelementList .eq. ST_NOHANDLE) then
        cycle
      end if
      
      call storage_getbase_int (&
          p_rdiscrSubdomain%RelementDistr(ieldist)%h_IelementList,p_IelementListNew)
          
      ielnew = 0
      
      do i=1,size(Ielements)
      
        ! Get the old element
        ielold = Ielements(i)
        
        if (associated(p_IelementDistrOld)) then
          ! Insert it to its corresponding element distribution
          if (p_IelementDistrOld(ielold) .eq. ieldist) then
          
            ! Global data.
            !
            ! New element in this list.
            ielnew = ielnew + 1
                
            p_IelementDistrNew(i) = ieldist
            
            if (associated(p_IelementCounterNew)) then
              p_IelementCounterNew (i) = ielnew
            end if
            
            ! local data
            p_IelementListNew(ielnew) = i
          
          end if
        else
          ! Uniform discretisation. Always insert.

          ! Global data.
          !
          ! New element in this list.
          ielnew = ielnew + 1
              
          ! local data
          p_IelementListNew(ielnew) = i
        end if
        
      end do ! i
      
      ! Are there any elements? If not, release the memory.
      if (ielnew .gt. 0) then
        p_rdiscrSubdomain%RelementDistr(ieldist)%NEL = ielnew
        call storage_realloc ("sdfem_createSubdFemFromElList",ielnew,&
            p_rdiscrSubdomain%RelementDistr(ieldist)%h_IelementList,&
            ST_NEWBLOCK_NOINIT)
      else
        call storage_free (p_rdiscrSubdomain%RelementDistr(ieldist)%h_IelementList)
      end if
    
    end do ! iellist
    
    ! -----------------------------------------------------
    ! Build up the DOF mapping structures
    ! -----------------------------------------------------
    
    ! At first, figure out how many DOFs we have on the subdomain.
    ndof = dof_igetNDofGlob(p_rdiscrSubdomain)
    rsdfemFeSpace%ndof = ndof
    
    call storage_new ("sdfem_createSubdFemFromElList", &
        "h_IglobalDofs", ndof, ST_INT, rsdfemFeSpace%h_IglobalDofs, &
        ST_NEWBLOCK_NOINIT)
    call storage_getbase_int (rsdfemFeSpace%h_IglobalDofs,p_IglobalDofs)
    
    ! Get the maximum number of local DOFs
    ndof = 0
    do ieldist = 1,p_rdiscrSubdomain%inumFESpaces
      ndof = max(ndof,&
                 elem_igetNDofLoc(p_rdiscrSubdomain%RelementDistr(ieldist)%celement))
    end do
    
    ! For every element, determine the global DOFs in the old and new
    ! discretisation. Save a reference to the old DOFs for the new DOFs.
    allocate(p_IdofsOld(ndof))
    allocate(p_IdofsNew(ndof))
    ieldist = 1
    do i=1,size(Ielements)
    
      ! Original element
      ielold = Ielements(i)
      
      ! Corresponding element distribution
      if (associated(p_IelementDistrNew)) then
        ieldist = p_IelementDistrNew(i)
      end if

      ! Number of local DOFs
      ndof = elem_igetNDofLoc(p_rdiscrSubdomain%RelementDistr(ieldist)%celement)
      
      ! Get the old and new global DOFs
      call dof_locGlobMapping(rdiscretisation, ielold, p_IdofsOld)
      call dof_locGlobMapping(p_rdiscrSubdomain, i, p_IdofsNew)
    
      ! Save the link between them.
      ! Overwrite any previously calculated value.
      ! (May be the same value is calculated more than once, but we do not
      ! care, it is always the same.)
      do idof = 1,ndof
        p_IglobalDofs(p_IdofsNew(idof)) = p_IdofsOld(idof)
      end do
    
    end do

    deallocate(p_IdofsOld)
    deallocate(p_IdofsNew)
    
  end subroutine

  ! ***************************************************************************

!<subrotine>
  
  subroutine sdfem_createSubdFemFrom2DBdReg(rsdfemFeSpace,rdiscretisation,rbdRegion,&
      rtriaSubdomain)

!<description>
  ! Creates a subdomain-FEM-structure from all elements that touch
  ! a 2d boundary region.
!</description>

!<input>
  ! Underlying discretisation of the full domain.
  type(t_spatialDiscretisation), intent(in), target :: rdiscretisation
  
  ! Boundary region. All elements touching this are extracted.
  type(t_boundaryRegion), intent(in) :: rbdRegion
  
  ! OPTIONAL: Triangulation of the subdomain.
  type(t_triangulation), target, optional :: rtriaSubdomain
!</input>

!<output>
  ! Structure that represents the FEM space on the subdomain.
  type(t_sdfemFeSpace), intent(out) :: rsdfemFeSpace
!</output>

!</subroutine>

    ! lcoal variables
    integer :: ncount,i,j
    integer, dimension(:), allocatable :: p_Ielements
    
    ! Get all the elements in the boundary refion
    call bcasm_getElementsInBdRegion (rdiscretisation%p_rtriangulation,&
        rbdRegion,ncount)

    if (ncount .eq. 0) return ! Nothing to do
    
    allocate (p_Ielements(ncount))
    call bcasm_getElementsInBdRegion (rdiscretisation%p_rtriangulation,&
        rbdRegion,ncount,IelList=p_Ielements)
    
    ! Remove duplicates
    i = 1
    j = 1
    do i = 2,ncount
      if (p_Ielements(j) .ne. p_Ielements(i)) then
        j = j + 1
        p_Ielements(j) = p_Ielements(i)
      end if
    end do
    
    ! Extract the subdomain containing these elements.
    call sdfem_createSubdFem(rsdfemFeSpace,rdiscretisation,p_Ielements(1:j))
    deallocate(p_Ielements)
    
  end subroutine

  ! ***************************************************************************

!<subrotine>
  
  subroutine sdfem_createSubdFemManif1D(rsdfemFeSpace,rdiscretisation,&
      rtriaSubdomain)

!<description>
  ! Creates a subdomain-FEM-structure for the manifold created by the
  ! points on the boundary of a domain.
!</description>

!<input>
  ! Underlying discretisation of the full domain.
  type(t_spatialDiscretisation), intent(in), target :: rdiscretisation
  
  ! OPTIONAL: Triangulation of the subdomain.
  ! Describes the manifold of all boundary components of the domain.
  type(t_triangulation), target, optional :: rtriaSubdomain
!</input>

!<output>
  ! Structure that represents the FEM space on the subdomain.
  type(t_sdfemFeSpace), intent(out) :: rsdfemFeSpace
!</output>

!</subroutine>

    ! Save references to the new structure and create nonexisting
    ! substructures.
    rsdfemFeSpace%p_rdiscrGlobal => rdiscretisation
    rsdfemFeSpace%p_rboundary => rdiscretisation%p_rboundary
    allocate(rsdfemFeSpace%p_rdiscrSubdomain)
    
    ! Create a triangulation for the subdomain.
    if (present(rtriaSubdomain)) then
      rsdfemFeSpace%p_rtriaSubdomain => rtriaSubdomain
      
      ! Prevent the structure from being deallocated
      rsdfemFeSpace%cflags = ior (rsdfemFeSpace%cflags,SDFEM_SHAREDTRIA)
    else
      ! Create the subdomain triangulation
      allocate(rsdfemFeSpace%p_rtriaSubdomain)
      call tria_createRawBdryTria2D(&
          rsdfemFeSpace%p_rtriaSubdomain, rdiscretisation%p_rtriangulation)
    end if

    ! -----------------------------------------------------
    ! Create the discretisation structure.
    ! -----------------------------------------------------
    call spdiscr_createCompDiscrManif1D(rsdfemFeSpace%p_rdiscrSubdomain,&
        rdiscretisation,rsdfemFeSpace%p_rtriaSubdomain)

    ! -----------------------------------------------------
    ! Build up the DOF mapping structures
    ! -----------------------------------------------------
    
    ! Get the DOFs on the boundary. They are exactly in the
    ! order which is necessary for the 1D discretisation.
    call bcasm_getDOFsOnBoundary (rdiscretisation, &
        rsdfemFeSpace%h_IglobalDofs,rsdfemFeSpace%ndof)
    
  end subroutine

  ! ***************************************************************************

!<subrotine>
  
  subroutine spdiscr_createCompDiscrManif1D(rdiscrManifold,rdiscretisation,&
      rtriaManifold)

!<description>
  ! Creates a "compatible" discretisation structure for the manifold created 
  ! by the points on the boundary of a domain. The new discretisation 
  ! is created compatible to rdiscretisation on the boudary.
!</description>

!<input>
  ! Underlying discretisation of the full domain.
  type(t_spatialDiscretisation), intent(in), target :: rdiscretisation
  
  ! Triangulation of the manifold.
  ! Describes the manifold of all boundary components of the domain.
  type(t_triangulation), target, optional :: rtriaManifold
!</input>

!<output>
  ! Discretisation structure for the discretisation of the buondary.
  type(t_spatialDiscretisation), intent(inout), target :: rdiscrManifold
!</output>

!</subroutine>
    
    ! local variables
    integer :: i,ieldist
    integer, dimension(:), pointer :: p_IelementDistrOld,p_IelementCounterOld
    integer, dimension(:), pointer :: p_IelementDistrNew,p_IelementCounterNew
    integer, dimension(:), pointer :: p_IelementListNew

    ! -----------------------------------------------------
    ! Duplicate the discretisation, set up the new one
    ! -----------------------------------------------------
    call spdiscr_duplicateDiscrSc2 (&
        rdiscretisation, rdiscrManifold, .false.)
        
    ! Change the element IDs to compatible 1D elements IDs.
    do i=1,rdiscrManifold%inumFESpaces
      
      select case (rdiscrManifold%RelementDistr(i)%celement)
      
      case (EL_P0_2D,   EL_Q0_2D)
        rdiscrManifold%RelementDistr(i)%celement = EL_P0_1D
        
      case (EL_P1_2D,   EL_Q1_2D)
        rdiscrManifold%RelementDistr(i)%celement = EL_P1_1D

      case (EL_P2_2D,   EL_Q2_2D)
        rdiscrManifold%RelementDistr(i)%celement = EL_P2_1D

!      case (EL_P1T,     EL_QP1)
!        rdiscrManifold%RelementDistr(i)%celement = EL_DG_T1_1D
!        
!      case (EL_E031_2D, EL_EM31_2D,&
!            EL_QP1_2D,  EL_E030_2D,     EL_EM30_2D)
!        rdiscrManifold%RelementDistr(i)%celement = EL_DG_T2_1D
      
      case default
        call output_line("Incompatible element",&
            OU_CLASS_ERROR,OU_MODE_STD,"sdfem_createSubdFemManif1D")
        call sys_halt()
        
      end select
      
    end do
    
    ! Rebuild the element sets
    rdiscrManifold%RelementDistr(:)%NEL = 0
    rdiscrManifold%p_rtriangulation => rtriaManifold
    
    if (rdiscrManifold%h_IelementDistr .ne. ST_NOHANDLE) then
      call storage_getbase_int (rdiscrManifold%h_IelementDistr,p_IelementDistrOld)
      call storage_getbase_int (rdiscrManifold%h_IelementDistr,p_IelementDistrNew)
    else
      ! Uniform discretisation
      nullify(p_IelementDistrOld)
      nullify(p_IelementDistrNew)
    end if
    
    if (rdiscrManifold%h_IelementCounter .ne. ST_NOHANDLE) then
      call storage_getbase_int (rdiscrManifold%h_IelementCounter,p_IelementCounterOld)
      call storage_getbase_int (rdiscrManifold%h_IelementCounter,p_IelementCounterNew)
    else
      nullify(p_IelementCounterOld)
      nullify(p_IelementCounterNew)
    end if
    
    if (rdiscrManifold%inumFESpaces .eq. 1) then
    
      ! Uniform discretisation.
      ! All elements of the new triangulation are put to the one and only
      ! element distribution.
      
      if (rdiscrManifold%p_rtriangulation%NEL .gt. 0) then

        rdiscrManifold%RelementDistr(1)%NEL = rdiscrManifold%p_rtriangulation%NEL
        call storage_realloc ("sdfem_createSubdFemManif1D",&
            rdiscrManifold%p_rtriangulation%NEL,&
            rdiscrManifold%RelementDistr(1)%h_IelementList,&
            ST_NEWBLOCK_NOINIT)
            
        call storage_getbase_int (&
            rdiscrManifold%RelementDistr(1)%h_IelementList,p_IelementListNew)
        
        do i=1,rdiscrManifold%p_rtriangulation%NEL
          p_IelementListNew(i) = i
        end do
      
      else

        call storage_free (rdiscrManifold%RelementDistr(ieldist)%h_IelementList)

      end if
      
    else
    
      ! Loop over the elements on the boundary. Determine their
      ! element distribution and plug them in the corresponding container.
      !
      ! ... to be implementsd

      call output_line("Mixed discretisations not yet supported.",&
          OU_CLASS_ERROR,OU_MODE_STD,"sdfem_createSubdFemManif1D")
      call sys_halt()
    
    end if

  end subroutine

  ! ***************************************************************************

!<subrotine>
  
  subroutine sdfem_releaseSubdFem(rsdfemFeSpace)

!<description>
  ! Releases a subdomain-FEM structure.
!</description>

!<inputoutput>
  ! Structure to release
  type(t_sdfemFeSpace), intent(inout) :: rsdfemFeSpace
!</inputoutput>

!</subroutine>

    ! Reset all data, release memory
    
    if (associated(rsdfemFeSpace%p_rdiscrSubdomain)) then
      call spdiscr_releaseDiscr(rsdfemFeSpace%p_rdiscrSubdomain)
      deallocate(rsdfemFeSpace%p_rdiscrSubdomain)
    end if
    
    if (iand(rsdfemFeSpace%cflags,SDFEM_SHAREDTRIA) .eq. 0) then
      if (associated(rsdfemFeSpace%p_rtriaSubdomain)) then
        call tria_done(rsdfemFeSpace%p_rtriaSubdomain)
        deallocate(rsdfemFeSpace%p_rtriaSubdomain)
      end if
    else
      nullify(rsdfemFeSpace%p_rtriaSubdomain)
    end if
    
    if (rsdfemFeSpace%h_IglobalDofs .ne. ST_NOHANDLE) then
      call storage_free (rsdfemFeSpace%h_IglobalDofs)
    end if
    
    rsdfemFeSpace%cflags = 0
    nullify(rsdfemFeSpace%p_rdiscrGlobal)
    nullify(rsdfemFeSpace%p_rboundary)
    rsdfemFeSpace%ndof = 0

  end subroutine

  ! ***************************************************************************

!<subrotine>
  
  subroutine sdfem_extractData(rsdfemFeSpace,rvectorSource,rvectorDest)

!<description>
  ! Extracts the DOFs of the subdomain to a vector.
!</description>

!<input>
  ! Subdomain-FEM-structure that defines the mapping between
  ! the global and the sub-domain
  type(t_sdfemFeSpace), intent(in) :: rsdfemFeSpace
  
  ! Source vector, defined on the global domain
  type(t_vectorScalar), intent(in) :: rvectorSource
!</input>

!<inputoutput>
  ! Destination vector that receives the data from the subdomain
  type(t_vectorScalar), intent(inout) :: rvectorDest
!</inputoutput>

!</subroutine>

    ! local data
    integer :: i
    real(DP), dimension(:), pointer :: p_Ddata1,p_Ddata2
    integer, dimension(:), pointer :: p_IglobalDofs

    ! Cancel if nothing to do or error.
    if (rsdfemFeSpace%h_IglobalDofs .eq. ST_NOHANDLE) return
    
    if ((rvectorSource%cdataType .ne. ST_DOUBLE) .or. &
        (rvectorDest%cdataType .ne. ST_DOUBLE) .or. &
        (rvectorSource%cdataType .ne. rvectorDest%cdataType)) then
      call output_line("Unsupported precision",&
          OU_CLASS_ERROR,OU_MODE_STD,"sdfem_extractData")
      call sys_halt()
    end if

    if (rvectorSource%bisSorted .or. rvectorDest%bisSorted) then
      call output_line("Sorted vectors not supported",&
          OU_CLASS_ERROR,OU_MODE_STD,"sdfem_extractData")
      call sys_halt()
    end if
    
    ! Get the data arrays
    call lsyssc_getbase_double (rvectorSource,p_Ddata1)
    call lsyssc_getbase_double (rvectorDest,p_Ddata2)
    call storage_getbase_int (rsdfemFeSpace%h_IglobalDofs,p_IglobalDofs)
    
    ! Transfer the data
    do i=1,rsdfemFeSpace%ndof
      p_Ddata2(i) = p_Ddata1(p_IglobalDofs(i))
    end do

  end subroutine

  ! ***************************************************************************

!<subrotine>
  
  subroutine sdfem_linearCombToGlobal(rsdfemFeSpace,rvectorSource,rvectorDest,dcx,dcy)

!<description>
  ! Linear combination of a subdomain vector into a global vector:
  !
  !    rvectorDest = cx * rvectorSource +  cy * rvectorDest
!</description>

!<input>
  ! Subdomain-FEM-structure that defines the mapping between
  ! the global and the sub-domain
  type(t_sdfemFeSpace), intent(in) :: rsdfemFeSpace
  
  ! Source vector, defined on the subdomain
  type(t_vectorScalar), intent(in) :: rvectorSource
  
  ! Weight for rvectorSource
  real(DP), intent(in) :: dcx

  ! Weight for rvectorDest
  real(DP), intent(in) :: dcy
!</input>

!<inputoutput>
  ! Destination vector.
  type(t_vectorScalar), intent(inout) :: rvectorDest
!</inputoutput>

!</subroutine>

    ! local data
    integer :: i
    real(DP), dimension(:), pointer :: p_Ddata1,p_Ddata2
    integer, dimension(:), pointer :: p_IglobalDofs

    ! Cancel if nothing to do or error.
    if (rsdfemFeSpace%h_IglobalDofs .eq. ST_NOHANDLE) return
    
    if ((rvectorSource%cdataType .ne. ST_DOUBLE) .or. &
        (rvectorDest%cdataType .ne. ST_DOUBLE) .or. &
        (rvectorSource%cdataType .ne. rvectorDest%cdataType)) then
      call output_line("Unsupported precision",&
          OU_CLASS_ERROR,OU_MODE_STD,"sdfem_extractData")
      call sys_halt()
    end if

    if (rvectorSource%bisSorted .or. rvectorDest%bisSorted) then
      call output_line("Sorted vectors not supported",&
          OU_CLASS_ERROR,OU_MODE_STD,"sdfem_extractData")
      call sys_halt()
    end if
    
    ! Get the data arrays
    call lsyssc_getbase_double (rvectorSource,p_Ddata1)
    call lsyssc_getbase_double (rvectorDest,p_Ddata2)
    call storage_getbase_int (rsdfemFeSpace%h_IglobalDofs,p_IglobalDofs)
    
    ! Linear combination
    do i=1,rsdfemFeSpace%ndof
      p_Ddata2(p_IglobalDofs(i)) = dcx * p_Ddata1(i) + dcy * p_Ddata2(p_IglobalDofs(i))
    end do

  end subroutine

  ! ***************************************************************************

!<subrotine>
  
  subroutine sdfem_linearCombToLocal(rsdfemFeSpace,rvectorSource,rvectorDest,dcx,dcy)

!<description>
  ! Linear combination of a global vector into a subdomain vector:
  !
  !    rvectorDest = cx * rvectorSource +  cy * rvectorDest
!</description>

!<input>
  ! Subdomain-FEM-structure that defines the mapping between
  ! the global and the sub-domain
  type(t_sdfemFeSpace), intent(in) :: rsdfemFeSpace
  
  ! Source vector, defined on the global domain
  type(t_vectorScalar), intent(in) :: rvectorSource
  
  ! Weight for rvectorSource
  real(DP), intent(in) :: dcx

  ! Weight for rvectorDest
  real(DP), intent(in) :: dcy
!</input>

!<inputoutput>
  ! Destination vector, defined on the subdomain
  type(t_vectorScalar), intent(inout) :: rvectorDest
!</inputoutput>

!</subroutine>

    ! local data
    integer :: i
    real(DP), dimension(:), pointer :: p_Ddata1,p_Ddata2
    integer, dimension(:), pointer :: p_IglobalDofs

    ! Cancel if nothing to do or error.
    if (rsdfemFeSpace%h_IglobalDofs .eq. ST_NOHANDLE) return
    
    if ((rvectorSource%cdataType .ne. ST_DOUBLE) .or. &
        (rvectorDest%cdataType .ne. ST_DOUBLE) .or. &
        (rvectorSource%cdataType .ne. rvectorDest%cdataType)) then
      call output_line("Unsupported precision",&
          OU_CLASS_ERROR,OU_MODE_STD,"sdfem_extractData")
      call sys_halt()
    end if

    if (rvectorSource%bisSorted .or. rvectorDest%bisSorted) then
      call output_line("Sorted vectors not supported",&
          OU_CLASS_ERROR,OU_MODE_STD,"sdfem_extractData")
      call sys_halt()
    end if
    
    ! Get the data arrays
    call lsyssc_getbase_double (rvectorSource,p_Ddata1)
    call lsyssc_getbase_double (rvectorDest,p_Ddata2)
    call storage_getbase_int (rsdfemFeSpace%h_IglobalDofs,p_IglobalDofs)
    
    ! Linear combination
    do i=1,rsdfemFeSpace%ndof
      p_Ddata2(i) = dcx * p_Ddata1(p_IglobalDofs(i)) + dcy * p_Ddata2(i)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdiscr_duplicateDiscrSc2 (rsourceDiscr, rdestDiscr, bshare)

!<description>
  ! This routine creates a copy of the discretisation structure rsourceDiscr.
  ! Depending on bshare, the destination structure rdestDiscr will either
  ! obtain a "simple" copy (i.e. sharing all handles and all information
  ! with rsourceDiscr) or a separate copy (which costs memory for all the
  ! element information!).
!</description>

!<input>
  ! A source discretisation structure that should be used as template
  type(t_spatialDiscretisation), intent(in) :: rsourceDiscr

  ! OPTIONAL: Whether the new discretisation structure should share its information
  ! with rsourceDiscr.
  ! =FALSE: Create a complete copy of rsourceDiscr which is independent
  !  of rsourceDiscr.
  ! =TRUE: The new discretisation will not be a complete new structure, but a
  !  "derived" structure, i.e. it uses the same dynamic information
  !  (handles and therefore element lists) as rsourceDiscr.
  ! If not specified, TRUE is assumed.
  logical, intent(in), optional :: bshare
!</input>

!<output>
  ! The new discretisation structure. Any old existing information in rdestDiscr
  ! is released if necessary.
  type(t_spatialDiscretisation), intent(inout), target :: rdestDiscr
!</output>

!</subroutine>

    logical :: bshr
    integer :: i

    bshr = .true.
    if (present(bshare)) bshr = bshare

    ! Release old information if present
    call spdiscr_releaseDiscr(rdestDiscr)

    ! Currently, this routine supports only bshare=TRUE!
    if (bshr) then

      ! Copy all information
      rdestDiscr = rsourceDiscr

      ! Duplicate the element distribution structure
      allocate(rdestDiscr%RelementDistr(rdestDiscr%inumFESpaces))
      rdestDiscr%RelementDistr = rsourceDiscr%RelementDistr

      ! Mark the new discretisation structure as "copy", to prevent
      ! the dynamic information to be released.
      ! The dynamic information "belongs" to rdiscrSource and not to the
      ! newly created rdiscrDest!
      rdestDiscr%bisCopy = .true.

    else

      ! Copy all information
      rdestDiscr = rsourceDiscr
      
      ! This is not a copy
      rdestDiscr%bisCopy = .false.
      
      ! Duplicate the handles and the memory behind.
      if (rdestDiscr%h_IelementDistr .ne. ST_NOHANDLE) then
        rdestDiscr%h_IelementDistr = ST_NOHANDLE
        call storage_copy (rsourceDiscr%h_IelementDistr,rdestDiscr%h_IelementDistr)
        
        ! otherwise: uniform discretisation.
      end if

      rdestDiscr%h_IelementCounter = ST_NOHANDLE
      if (rsourceDiscr%h_IelementCounter .ne. ST_NOHANDLE) then
        call storage_copy (rsourceDiscr%h_IelementCounter,rdestDiscr%h_IelementCounter)
      end if

      ! Do not calculate any precompiled DOF mapping.
      ! May be implemented later.
      rdestDiscr%h_IelementDofs = ST_NOHANDLE
      rdestDiscr%h_IelementDofIdx = ST_NOHANDLE
      rdestDiscr%bprecompiledDofMapping = .false.
      rdestDiscr%ndof = 0

      ! Duplicate the element distribution structure
      allocate(rdestDiscr%RelementDistr(rdestDiscr%inumFESpaces))
      rdestDiscr%RelementDistr = rsourceDiscr%RelementDistr

      do i=1,rdestDiscr%inumFESpaces
        rdestDiscr%RelementDistr(i)%h_IelementList = ST_NOHANDLE
        call storage_copy (rsourceDiscr%RelementDistr(i)%h_IelementList,&
            rdestDiscr%RelementDistr(i)%h_IelementList)
      end do

    end if

  end subroutine

  ! ***************************************************************************

!<subrotine>
  
  subroutine sdfem_createSubdBlkFemElList(rsdfemBlkFeSpace,rtriangulation,Ielements,&
      rblockDiscr,rboundary)

!<description>
  ! Creates a block-subdomain-FEM-structure for a set of elements. Ielements is a
  !´list of elements in the triangulation connected to rtriangulation.
  ! The routine extracts these elements, creates a new triangulatio and 
  ! saves all this in rsdfemFeSpace.
  ! Specifying rblockDiscr additionally sets up all the sub-discretisations
  ! in rsdfemBlkFeSpace. If rblockDiscr is not specified, the caller has
  ! to add all blocks one after the other using sdfem_subdBlkFemNewBlock.
!</description>

!<input>
  ! Underlying triangulation
  type(t_triangulation), target, optional :: rtriangulation
  
  ! List of elements that form the subdomain.
  integer, dimension(:), intent(in) :: Ielements
  
  ! OPTIONAL: Block discretisation.
  type(t_blockDiscretisation), intent(in), target, optional :: rblockDiscr
  
  ! OPTIONAL: Definition of the domain. This can be specified if
  ! rblockDiscr is not specified.
  type(t_boundary), intent(in), target, optional :: rboundary
!</input>

!<output>
  ! Structure that represents the FEM space on the subdomain.
  type(t_sdfemBlkFeSpace), intent(out) :: rsdfemBlkFeSpace
!</output>

!</subroutine>

    ! local variables
    integer, dimension(:), pointer :: p_Ielements
    type(t_boundary), pointer :: p_rboundary
    integer :: i

    if (size(Ielements) .eq. 0) return

    ! Remember the element list
    call storage_new ("sdfem_createSubdBlkFemElList", &
        "h_Ielements", size(Ielements), ST_INT, rsdfemBlkFeSpace%h_Ielements, &
        ST_NEWBLOCK_NOINIT)
    call storage_getbase_int (rsdfemBlkFeSpace%h_Ielements,p_Ielements)
    call lalg_copyVector (Ielements,p_Ielements)
    
    ! Is a domain or boundary given?
    nullify(p_rboundary)
    if (present(rboundary)) p_rboundary => rboundary
    if (present(rblockDiscr)) p_rboundary => rblockDiscr%p_rboundary
    
    ! Remember the triangulation
    rsdfemBlkFeSpace%p_rtriaDomain => rtriangulation
    
    ! Create a triangulation for the subdomain.
    allocate(rsdfemBlkFeSpace%p_rtriaSubdomain)
    if (associated(p_rboundary)) then
      call tria_generateSubdomain(rtriangulation, Ielements, &
          rsdfemBlkFeSpace%p_rtriaSubdomain, p_rboundary)
    else
      call tria_generateSubdomain(rtriangulation, Ielements, &
          rsdfemBlkFeSpace%p_rtriaSubdomain)
    end if
    
    ! Discretisation given?
    rsdfemBlkFeSpace%nblocks = 0
    if (present(rblockDiscr)) then
    
      ! Set up the blocks
      allocate(rsdfemBlkFeSpace%p_RsdfemSpaces(rblockDiscr%ncomponents))
      
      ! Add the spatial discretisation structures
      do i=1,rblockDiscr%ncomponents
        call sdfem_subdBlkFemNewBlock(rsdfemBlkFeSpace,rblockDiscr%RspatialDiscr(i))
      end do
    
    else
    
      ! Allocate some space in-advance
      allocate(rsdfemBlkFeSpace%p_RsdfemSpaces(8))
    
    end if

  end subroutine

  ! ***************************************************************************

!<subrotine>
  
  subroutine sdfem_createSubdBlkFem2DBdReg(rsdfemBlkFeSpace,rtriangulation,rbdRegion,&
      rblockDiscr,rboundary)

!<description>
  ! Creates a subdomain-block-FEM-structure from all elements that touch
  ! a 2d boundary region.
!</description>

!<input>
  ! Underlying triangulation
  type(t_triangulation), target, optional :: rtriangulation
  
  ! Boundary region. All elements touching this are extracted.
  type(t_boundaryRegion), intent(in) :: rbdRegion
  
  ! OPTIONAL: Block discretisation.
  type(t_blockDiscretisation), intent(in), target, optional :: rblockDiscr
  
  ! OPTIONAL: Definition of the domain. This can be specified if
  ! rblockDiscr is not specified.
  type(t_boundary), intent(in), target, optional :: rboundary
!</input>

!<output>
  ! Structure that represents the FEM space on the subdomain.
  type(t_sdfemBlkFeSpace), intent(out) :: rsdfemBlkFeSpace
!</output>

!</subroutine>

    ! lcoal variables
    integer :: ncount,i,j
    integer, dimension(:), allocatable :: p_Ielements
    
    ! Get all the elements in the boundary refion
    call bcasm_getElementsInBdRegion (rblockDiscr%p_rtriangulation,&
        rbdRegion,ncount)

    if (ncount .eq. 0) return ! Nothing to do
    
    allocate (p_Ielements(ncount))
    call bcasm_getElementsInBdRegion (rblockDiscr%p_rtriangulation,&
        rbdRegion,ncount,IelList=p_Ielements)
    
    ! Remove duplicates
    i = 1
    j = 1
    do i = 2,ncount
      if (p_Ielements(j) .ne. p_Ielements(i)) then
        j = j + 1
        p_Ielements(j) = p_Ielements(i)
      end if
    end do
    
    ! Extract the subdomain containing these elements.
    call sdfem_createSubdBlkFem(rsdfemBlkFeSpace,rtriangulation,p_Ielements(1:j),&
        rblockDiscr,rboundary)
    deallocate(p_Ielements)
    
  end subroutine

  ! ***************************************************************************

!<subrotine>
  
  subroutine sdfem_releaseSubdBlkFem(rsdfemBlkFeSpace)

!<description>
  ! Releases a subdomain-block-FEM-structure.
!</description>

!<inputoutput>
  ! Structure to be released
  type(t_sdfemBlkFeSpace), intent(inout) :: rsdfemBlkFeSpace
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i

    ! Release the substructures.
    nullify(rsdfemBlkFeSpace%p_rdiscrGlobal)
    
    do i=1,rsdfemBlkFeSpace%nblocks
      call sdfem_releaseSubdFem(rsdfemBlkFeSpace%p_RsdfemSpaces(i))
    end do
    
    rsdfemBlkFeSpace%nblocks = 0
    deallocate(rsdfemBlkFeSpace%p_RsdfemSpaces)
    
    if (rsdfemBlkFeSpace%h_Ielements .ne. ST_NOHANDLE) then
      call storage_free(rsdfemBlkFeSpace%h_Ielements)
    end if  
    
    if (associated(rsdfemBlkFeSpace%p_rtriaSubdomain)) then
      call tria_done(rsdfemBlkFeSpace%p_rtriaSubdomain)
      deallocate(rsdfemBlkFeSpace%p_rtriaSubdomain)
    end if
    
    nullify(rsdfemBlkFeSpace%p_rtriaDomain)
    
  end subroutine

  ! ***************************************************************************

!<subrotine>
  
  subroutine sdfem_subdBlkFemNewBlock(rsdfemBlkFeSpace,rdiscretisation)

!<description>
  ! Appends a discretisation structure to a subdomain-block-FE structure.
  ! This creates a new "local" discretisation structure for the subdomain
  ! defined by rsdfemBlkFeSpace.
  !
  ! After sdfem_createSubdBlkFem is called, this can be used to append
  ! a new discretisation structure to rsdfemBlkFeSpace, resulting in a new
  ! block.
!</description>

!<input>
  ! Space discretisation structure to be added as new block to rsdfemBlkFeSpace.
  type(t_spatialDiscretisation), intent(in) :: rdiscretisation
!</input>

!<inputoutput>
  ! Structure that represents the FEM space on the subdomain.
  type(t_sdfemBlkFeSpace), intent(inout) :: rsdfemBlkFeSpace
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_sdfemFeSpace), dimension(:), pointer :: p_RsdfemSpaces
    integer, dimension(:), pointer :: p_Ielements

    ! Reallocate the blocks if necessary.
    p_RsdfemSpaces => rsdfemBlkFeSpace%p_RsdfemSpaces
    if (rsdfemBlkFeSpace%nblocks .ge. size(p_RsdfemSpaces)) then
      allocate (p_RsdfemSpaces(rsdfemBlkFeSpace%nblocks+8))
      p_RsdfemSpaces(1:rsdfemBlkFeSpace%nblocks) = &
          rsdfemBlkFeSpace%p_RsdfemSpaces(1:rsdfemBlkFeSpace%nblocks)
      deallocate(rsdfemBlkFeSpace%p_RsdfemSpaces)
      rsdfemBlkFeSpace%p_RsdfemSpaces => p_RsdfemSpaces
    end if
    
    ! Add a block.
    rsdfemBlkFeSpace%nblocks = rsdfemBlkFeSpace%nblocks + 1

    ! Set up the new block
    call storage_getbase_int (rsdfemBlkFeSpace%h_Ielements,p_Ielements)
    call sdfem_createSubdFemFromElList(p_RsdfemSpaces(rsdfemBlkFeSpace%nblocks),&
        rdiscretisation,p_Ielements,rsdfemBlkFeSpace%p_rtriaSubdomain)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sdfem_createBlkHierarchyMeshH (rsdfemBlkFeHierarchy,rmeshHierarchy)

!<description>
  ! Creates a FE space subdomain hierarchy based on a mesh hierarchy.
  !
  ! The structures of the different levels are NOT initialised.
  ! The caller has to initialise them according to the submeshes which
  ! are needed.
!</description>

!<input>
  ! An underlying mesh hierarchy.
  type(t_meshHierarchy), intent(in) :: rmeshHierarchy
!</input>

!<output>
  ! Structure that represents a hierarchy of block FE spaces on
  ! subdomains.
  type(t_sdfemBlkFeHierarchy), intent(out) :: rsdfemBlkFeHierarchy
!</output>

!</subroutine>

    if (rmeshHierarchy%nlevels .le. 0) return

    ! Initialise the structure
    rsdfemBlkFeHierarchy%nlevels = rmeshHierarchy%nlevels
    allocate(rsdfemBlkFeHierarchy%p_RfeBlkSpaces(rsdfemBlkFeHierarchy%nlevels))

    ! Duplicate the mesh hierarchy.
    call mshh_initHierarchy(rmeshHierarchy,rsdfemBlkFeHierarchy%rmeshHierarchy)

    ! Subdomain structures on the levels are not initialised.

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sdfem_releaseBlkHierarchy (rsdfemBlkFeHierarchy)

!<description>
  ! Releases a FE space subdomain hierarchy.
!</description>

!<inputoutput>
  ! Structure to be released.
  type(t_sdfemBlkFeHierarchy), intent(out) :: rsdfemBlkFeHierarchy
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i

    if (rsdfemBlkFeHierarchy%nlevels .le. 0) return

    ! Release substructures
    do i=1,rsdfemBlkFeHierarchy%nlevels
      call sdfem_releaseSubdBlkFem(rsdfemBlkFeHierarchy%p_RfeBlkSpaces(i))
    end do
    deallocate(rsdfemBlkFeHierarchy%p_RfeBlkSpaces)
    
    call mshh_releaseHierarchy(rsdfemBlkFeHierarchy%rmeshHierarchy)
    rsdfemBlkFeHierarchy%nlevels = 0

  end subroutine

end module
