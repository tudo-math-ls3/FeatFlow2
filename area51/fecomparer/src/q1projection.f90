module q1projection

  ! Include basic Feat-2 modules
  use fsystem
  use genoutput
  use storage
  
  use cubature
  use triangulation
  use spatialdiscretisation
  
  use meshhierarchy
  use fespacehierarchybase
  use fespacehierarchy
  
  use linearsystemscalar
  use linearsystemblock
  
  use blockmatassemblybase
  use blockmatassembly
  use feevaluation
  use feevaluation2
  
  use collection
  
  implicit none
  
  public :: sol_projectToDGQ1
  
contains

  subroutine fcalcLocalIntegral (dintvalue,rassemblyData,rintegralAssembly,&
    npointsPerElement,nelements,revalVectors,rcollection)

  use collection
  use blockmatassemblybase

  ! Calculates the value of the integral for a set of elements.
  
  ! Returns the value of the integral
  real(DP), intent(out) :: dintvalue

  ! Data necessary for the assembly. Contains determinants and
  ! cubature weights for the cubature,...
  type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

  ! Structure with all data about the assembly
  type(t_bmaIntegralAssembly), intent(in) :: rintegralAssembly

  ! Number of points per element
  integer, intent(in) :: npointsPerElement

  ! Number of elements
  integer, intent(in) :: nelements

  ! Values of FEM functions automatically evaluated in the
  ! cubature points.
  type(t_fev2Vectors), intent(in) :: revalVectors

  ! User defined collection structure
  type(t_collection), intent(inout), target, optional :: rcollection
  
    ! local variables
    real(DP), dimension(:), pointer :: p_DdestDofs
    integer, dimension(:), pointer :: p_IdofIndex
    real(DP), dimension(:,:), pointer :: p_DsourceDofs
    integer :: isourcecomp,idestcomp,i,j,iel
    type(t_vectorScalar), pointer :: p_rvectorSource
    type(t_vectorScalar), pointer :: p_rvectorTarget
    
    ! Get the DOF index array
    call storage_getbase_int (rcollection%IquickAccess(1),p_IdofIndex)
    
    call output_line (trim(sys_siL(rassemblyData%p_IelementList(1),10)))
    
    ! Loop over the vector components
    isourcecomp = 1
    idestcomp = 1
    do
    
      ! Stop if maximum number of components is reached.
      if (isourcecomp .gt. rcollection%p_rvectorQuickAccess1%nblocks) exit
      if (idestcomp .gt. rcollection%p_rvectorQuickAccess2%nblocks) exit
    
      ! Source function
      p_rvectorSource => rcollection%p_rvectorQuickAccess1%RvectorBlock(isourcecomp)
      p_rvectorTarget => rcollection%p_rvectorQuickAccess2%RvectorBlock(idestcomp)

      ! Get the target component.
      call lsyssc_getbase_double (p_rvectorTarget,p_DdestDofs)
    
      ! Get the values in the corners.
      ! This is a temporary array we have to fill with data.
      
      p_DsourceDofs => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC)
      
      ! Loop over the elements
      do iel = 1,nelements
        ! On each element, evaluate the FE function in the points
        call fevl_evaluate (DER_FUNC, p_DsourceDofs(:,iel), p_rvectorSource, &
            rassemblyData%revalElementSet%p_DpointsReal(:,:,iel))
      end do
      
      ! Copy the corner values to the destination array.
      do i=1,nelements
        ! Current element, serves as index
        iel = rassemblyData%p_IelementList(i)
        do j=1,npointsPerElement
          p_DdestDofs(p_IdofIndex(iel)+j-1) = p_DsourceDofs(j,i)
        end do
      end do
      
      ! Next component
      isourcecomp = isourcecomp + 1
      idestcomp = idestcomp + 1
      
    end do
    
    ! Dummy return value
    dintvalue = 0.0_DP

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine sol_projectToDGQ1 (rsource,rtarget)

!<description>
  ! Projects a solution into the DG Q1 space.
!</description>

!<input>
  ! Source solution
  type(t_vectorBlock), intent(in), target :: rsource
!</input>

!<inputoutput>
  ! Target solution, must be DG-Q1.
  type(t_vectorBlock), intent(inout), target :: rtarget
!</inputoutput>

!</subroutine>

    type(t_triangulation), pointer :: p_rtriangulation
    type(t_fev2Vectors) :: revalVectors
    type(t_collection) :: rcollection
    type(t_scalarCubatureInfo) :: rcubatureInfo
    integer :: h_IdofCounter
    integer, dimension(:), pointer :: p_IdofCounter
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer :: i,j
    real(DP) :: dintvalue
   
    ! Count the number of vertices per cell, create an index array
    ! to the target vector.
    p_rtriangulation => rtarget%p_rblockDiscr%p_rtriangulation
    
    h_IdofCounter = ST_NOHANDLE
    call storage_new ("sol_projectToDGQ1", "index", p_rtriangulation%NEL, ST_INT, &
        h_IdofCounter, ST_NEWBLOCK_ZERO)
    
    call storage_getbase_int (h_IdofCounter,p_IdofCounter)
    call storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
    
    p_IdofCounter(1) = 1
    do i=2,p_rtriangulation%NEL
      do j=ubound(p_IverticesAtElement,1),1,-1
        if (p_IverticesAtElement(j,i) .ne. 0) then
          p_IdofCounter(i) = p_IdofCounter(i-1) + j
          exit
        end if
      end do
    end do
    
    ! We use the trapezoidal rule to get the corner values.
    call spdiscr_createDefCubStructure (rtarget%p_rblockDiscr%RspatialDiscr(1),&
        rcubatureInfo,CUB_GEN_AUTO_TRZ)
        
    rcollection%p_rvectorQuickAccess1 => rsource
    rcollection%p_rvectorQuickAccess2 => rtarget
    rcollection%IquickAccess(1) = h_IdofCounter

    ! Dummy-call to the integration. Directly work on our underlying mesh.
    ! Provide some temporary memory for intermediate calculations.
    call fev2_addDummyVectorToEvalList(revalVectors,1)
    call fev2_addVectorToEvalList(revalVectors,rtarget%RvectorBlock(1),0)
    
    call bma_buildIntegral (dintvalue,BMA_CALC_STANDARD,&
        fcalcLocalIntegral,rtarget%p_rblockDiscr%p_rtriangulation,&
        rcollection=rcollection,&
        revalVectors=revalVectors,rcubatureInfo=rcubatureInfo)
        
    call storage_free (h_IdofCounter)
    call fev2_releaseVectorList(revalVectors)
    call spdiscr_releaseCubStructure (rcubatureInfo)

  end subroutine

end module
