!##############################################################################
!# ****************************************************************************
!# <name> cchrzlumping </name>
!# ****************************************************************************
!#
!# <purpose>
!# ... 
!#
!# It contains the following set of routines:
!#
!# 1.) bilf_buildMatrixScalar_hrz
!#     -> ...
!#
!# 2.) bilf_assembleSubmeshMatrix9
!#     -> ...
!#
!# </purpose>
!##############################################################################

  module cchrzlumping
  
    use bilinearformevaluation
    use basicgeometry
    use boundary
    use boundaryaux
    use collection, only: t_collection
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
    use mprimitives
    use scalarpde
    use spatialdiscretisation
    use storage
    use transformation
    use triangulation
    use extstdassemblyinfo
    use perfconfig
  
  implicit none
  
  public :: bilf_buildMatrixScalar_hrz
  public :: bilf_assembleSubmeshMatrix9_hrz
  
  !************************************************************************

  ! global performance configuration
  type(t_perfconfig), target, save :: bilf_perfconfig

  !************************************************************************
  
  contains

!**************************************************************************************************
!
  subroutine bilf_buildMatrixScalar_hrz (rform,bclear,rmatrix,&
      rcubatureInfo, fcoeff_buildMatrixSc_sim, rcollection, rperfconfig)
!
!**************************************************************************************************
!<description>
! ....
!</description>

!<input>
  ! The bilinear form specifying the underlying PDE of the discretisation.
  type(t_bilinearForm), intent(in) :: rform
  
  ! Whether to clear the matrix before calculating the entries.
  ! If .FALSE., the new matrix entries are added to the existing entries.
  logical, intent(in) :: bclear
  
  ! (OPTINOAL:) A scalar cubature information structure that specifies the cubature
  ! formula(s) to use. If not specified, default settings are used.
  type(t_scalarCubatureInfo), intent(in), optional, target :: rcubatureInfo

  ! OPTIONAL: A collection structure. This structure is given to the
  ! callback function for nonconstant coefficients to provide additional
  ! information. 
  type(t_collection), intent(inout), target, optional :: rcollection
  
  ! OPTIONAL: A callback routine for nonconstant coefficient matrices.
  ! Must be present if the matrix has nonconstant coefficients!
  include '../../../kernel/DOFMaintenance/intf_coefficientMatrixSc.inc'
  optional :: fcoeff_buildMatrixSc_sim

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig
!</input>

!<inputoutput>
  ! The FE matrix. Calculated matrix entries are imposed to this matrix.
  type(t_matrixScalar), intent(inout) :: rmatrix
!</inputoutput>

!</subroutine>

  ! local variables
  type(t_bilfMatrixAssembly) :: rmatrixAssembly
  integer :: ielementDistr,icubatureBlock,NEL
  integer, dimension(:), pointer :: p_IelementList
!  type(t_scalarCubatureInfo), target :: rtempCubatureInfo
  type(t_scalarCubatureInfo), pointer :: p_rcubatureInfo
  integer(I32) :: celementTrial,celementTest,ccubature
  type(t_scalarCubatureInfo), target :: rtempCubatureInfo
  
  ! Pointer to the performance configuration
  type(t_perfconfig), pointer :: p_rperfconfig

  if (present(rperfconfig)) then
    p_rperfconfig => rperfconfig
  else
    p_rperfconfig => bilf_perfconfig
  end if

  ! The matrix must be unsorted, otherwise we can not set up the matrix.
  ! Note that we cannot switch off the sorting as easy as in the case
  ! of a vector, since there is a structure behind the matrix! So the caller
  ! has to make sure, the matrix is unsorted when this routine is called.
  if (rmatrix%isortStrategy .gt. 0) then
    call output_line ('Matrix-structure must be unsorted!', &
        OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalar')
    call sys_halt()
  end if

  if ((.not. associated(rmatrix%p_rspatialDiscrTest)) .or. &
      (.not. associated(rmatrix%p_rspatialDiscrTrial))) then
    call output_line ('No discretisation associated!', &
        OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalar')
    call sys_halt()
  end if

  ! If we do not have it, create a cubature info structure that
  ! defines how to do the assembly.
  if (.not. present(rcubatureInfo)) then
    call spdiscr_createDefCubStructure(rmatrix%p_rspatialDiscrTrial,&
        rtempCubatureInfo,CUB_GEN_DEPR_BILFORM)
    p_rcubatureInfo => rtempCubatureInfo
  else
    p_rcubatureInfo => rcubatureInfo
  end if

  ! Do we have a uniform triangulation? Would simplify a lot...
  select case (rmatrix%p_rspatialDiscrTest%ccomplexity)
  case (SPDISC_UNIFORM,SPDISC_CONFORMAL) 
    ! Uniform and conformal discretisations
    select case (rmatrix%cdataType)
    case (ST_DOUBLE) 
      ! Which matrix structure do we have?
      select case (rmatrix%cmatrixFormat) 
      case (LSYSSC_MATRIX9)
      
        ! Probably allocate/clear the matrix
        if (rmatrix%h_DA .eq. ST_NOHANDLE) then
          call lsyssc_allocEmptyMatrix(rmatrix,LSYSSC_SETM_ZERO)
        else
          if (bclear) call lsyssc_clearMatrix (rmatrix)
        end if
      
        ! Loop over the cubature blocks to discretise
        do icubatureBlock = 1,p_rcubatureInfo%ninfoBlockCount
        
          ! Get information about that block.
          call spdiscr_getStdDiscrInfo(icubatureBlock,p_rcubatureInfo,&
              rmatrix%p_rspatialDiscrTest,ielementDistr,celementTest,ccubature,NEL,p_IelementList)
          
          call spdiscr_getStdDiscrInfo(icubatureBlock,p_rcubatureInfo,&
              rmatrix%p_rspatialDiscrTrial,celement=celementTrial)
        
          ! Check if element distribution is empty
          if (NEL .le. 0 ) cycle

          ! Initialise a matrix assembly structure for that element distribution
          call bilf_initAssembly(rmatrixAssembly,rform,&
              celementTest,celementTrial,ccubature,min(p_rperfconfig%NELEMSIM,NEL),rperfconfig)
              
          ! Assemble the data for all elements in this element distribution
          call bilf_assembleSubmeshMatrix9_hrz (rmatrixAssembly,rmatrix,&
              p_IelementList,fcoeff_buildMatrixSc_sim,rcollection,rperfconfig)
          
          ! Release the assembly structure.
          call bilf_doneAssembly(rmatrixAssembly)
        end do
                                       
      case default
        call output_line ('Not supported matrix structure!', &
            OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalar')
        call sys_halt()
      end select
      
    case default
      call output_line ('Single precision matrices currently not supported!', &
          OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalar')
      call sys_halt()
    end select
    
  case default
    call output_line ('General discretisation not implemented!', &
        OU_CLASS_ERROR,OU_MODE_STD,'bilf_buildMatrixScalar')
    call sys_halt()
  end select
  
  ! Release the assembly structure if necessary.
  if (.not. present(rcubatureInfo)) then
    call spdiscr_releaseCubStructure(rtempCubatureInfo)
  end if
end subroutine
  
!**************************************************************************************************
!  
  subroutine bilf_assembleSubmeshMatrix9_hrz (rmatrixAssembly, rmatrix, IelementList,&
      fcoeff_buildMatrixSc_sim, rcollection, rperfconfig)
!      
!************************************************************************************************** 
!<description>

! Assembles the matrix entries for a submesh by integrating over the domain.

!</description>
 
!<input>
  
  ! List of elements where to assemble the bilinear form.
  integer, dimension(:), intent(in), target :: IelementList
  
  ! OPTIONAL: A callback routine for nonconstant coefficient matrices.
  ! Must be present if the matrix has nonconstant coefficients!
  include '../../../kernel/DOFMaintenance/intf_coefficientMatrixSc.inc'
  optional :: fcoeff_buildMatrixSc_sim
  
  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), target, optional :: rperfconfig

!</input>

!<inputoutput>
  
  ! A matrix assembly structure prepared with bilf_initAssembly.
  type(t_bilfMatrixAssembly), intent(inout), target :: rmatrixAssembly
  
  ! A matrix where to assemble the contributions to.
  type(t_matrixScalar), intent(inout) :: rmatrix
  
  ! OPTIONAL: A pointer to a collection structure. This structure is given to the
  ! callback function for nonconstant coefficients to provide additional
  ! information. 
  type(t_collection), intent(inout), target, optional :: rcollection
  
!</inputoutput>
  
!</subroutine>
  
    ! local variables, used by all processors
    real(DP), dimension(:), pointer :: p_DA
    integer :: indofTest,indofTrial,ncubp
    
    ! local data of every processor when using OpenMP
    integer :: IELset,IELmax
    integer :: iel,icubp,ialbet,ia,ib,idofe,jdofe
    real(DP) :: domega,daux,db,dalpha,ddgmass,dtemass
    integer(I32) :: cevaluationTag
    type(t_bilfMatrixAssembly), target :: rlocalMatrixAssembly
    type(t_domainIntSubset) :: rintSubset
    integer, dimension(:,:,:), pointer :: p_Kentry
    real(DP), dimension(:,:,:), pointer :: p_Dentry
    real(DP), dimension(:,:), pointer :: p_Ddetj
    real(DP), dimension(:), pointer :: p_Domega
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTest
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial
    real(DP), dimension(:,:,:), pointer :: p_Dcoefficients
    real(DP), dimension(:), pointer :: p_DcoefficientsBilf
    integer, dimension(:,:), pointer :: p_IdofsTest
    integer, dimension(:,:), pointer :: p_IdofsTrial
    type(t_evalElementSet), pointer :: p_revalElementSet
    integer, dimension(:,:),pointer :: p_Idescriptors
    type(t_perfconfig), pointer :: p_rperfconfig
  
    ! Get some pointers for faster access
    call lsyssc_getbase_double (rmatrix,p_DA)
    indofTest = rmatrixAssembly%indofTest
    indofTrial = rmatrixAssembly%indofTrial
    ncubp = rmatrixAssembly%ncubp

    if (present(rperfconfig)) then
      p_rperfconfig => rperfconfig
    else
      p_rperfconfig => bilf_perfconfig
    end if
    ! Copy the matrix assembly data to the local matrix assembly data,
    ! where we can allocate memory.
    ! For single processor machines, this is actually boring and nonsense.
    ! But using OpenMP, here we get a local copy of the matrix
    ! assembly structure to where we can add some local data which
    ! is released upon return without changing the original matrix assembly
    ! stucture or disturbing the data of the other processors.
    rlocalMatrixAssembly = rmatrixAssembly
    call bilf_allocAssemblyData(rlocalMatrixAssembly)
    
    ! Get some more pointers to local data.
    p_Kentry            => rlocalMatrixAssembly%p_Kentry
    p_Dentry            => rlocalMatrixAssembly%p_Dentry
    p_Domega            => rlocalMatrixAssembly%p_Domega
    p_DbasTest          => rlocalMatrixAssembly%p_DbasTest
    p_DbasTrial         => rlocalMatrixAssembly%p_DbasTrial
    p_Dcoefficients     => rlocalMatrixAssembly%p_Dcoefficients
    p_Idescriptors      => rlocalMatrixAssembly%rform%Idescriptors
    p_IdofsTest         => rlocalMatrixAssembly%p_IdofsTest
    p_IdofsTrial        => rlocalMatrixAssembly%p_IdofsTrial
    p_revalElementSet   => rlocalMatrixAssembly%revalElementSet
    p_DcoefficientsBilf => rlocalMatrixAssembly%rform%Dcoefficients
        
    ! Loop over the elements - blockwise.
    !
    ! Open-MP-Extension: Each loop cycle is executed in a different thread,
    ! so nelementsPerBlock local matrices are simultaneously calculated in the
    ! inner loop(s).
    ! The blocks have all the same size, so we can use static scheduling.
    !
    !%OMP do schedule(static,1)
    do IELset = 1, size(IelementList), rlocalMatrixAssembly%nelementsPerBlock
    
      ! We always handle nelementsPerBlock elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most BILF_NELEMSIM
      ! elements simultaneously.
      
      IELmax = min(size(IelementList),IELset-1+rlocalMatrixAssembly%nelementsPerBlock)
    
      ! --------------------- DOF SEARCH PHASE ------------------------
    
      ! The outstanding feature with finite elements is: A basis
      ! function for a DOF on one element has common support only
      ! with the DOF`s on the same element! E.g. for Q1:
      !
      !        #. . .#. . .#. . .#
      !        .     .     .     .
      !        .  *  .  *  .  *  .
      !        #-----O-----O. . .#
      !        |     |     |     .
      !        |     | iel |  *  .
      !        #-----X-----O. . .#
      !        |     |     |     .
      !        |     |     |  *  .
      !        #-----#-----#. . .#
      !
      ! --> On element iel, the basis function at "X" only interacts
      !     with the basis functions in "O". Elements in the 
      !     neighbourhood ("*") have no support, therefore we only have
      !     to collect all "O" DOF`s.
      !
      ! Calculate the global DOF`s into IdofsTrial / IdofsTest.
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF`s of our NELEMSIM elements simultaneously.
      call dof_locGlobMapping_mult(rmatrix%p_rspatialDiscrTest, &
          IelementList(IELset:IELmax), p_IdofsTest)
                                   
      ! If the DOF`s for the test functions are different, calculate them, too.
      if (.not. rlocalMatrixAssembly%bIdenticalTrialAndTest) then
        call dof_locGlobMapping_mult(rmatrix%p_rspatialDiscrTrial, &
            IelementList(IELset:IELmax), p_IdofsTrial)
      end if
      
      ! ------------------- LOCAL MATRIX SETUP PHASE -----------------------
      
      ! For the assembly of the global matrix, we use a "local"
      ! approach. At first we build a "local" system matrix according
      ! to the current element. This contains all additive
      ! contributions of element iel, which are later added at the
      ! right positions to the elements in the global system matrix.
      !
      ! We have indofTrial trial DOF`s per element and
      ! indofTest test DOF`s per element. Therefore there are
      ! indofTrial*indofTest tupel of basis-/testfunctions (phi_i,psi_j) 
      ! "active" (i.e. have common support) on our current element, each 
      ! giving an additive contribution to the system matrix.
      !
      ! We build a quadratic indofTrial*indofTest local matrix:
      ! Kentry(1..indofTrial,1..indofTest) receives the position 
      ! in the global system matrix, where the corresponding value 
      ! has to be added to.
      ! (The corresponding contributions can be saved separately, 
      ! but we directly add them to the global matrix in this 
      ! approach.)
      !
      ! We build local matrices for all our elements 
      ! in the set simultaneously. Get the positions of the local matrices
      ! in the global matrix.
      call bilf_getLocalMatrixIndices (rmatrix,p_IdofsTest,p_IdofsTrial,p_Kentry,&
          ubound(p_IdofsTest,1),ubound(p_IdofsTrial,1),IELmax-IELset+1)
      
      ! -------------------- ELEMENT EVALUATION PHASE ----------------------
      
      ! Ok, we found the positions of the local matrix entries
      ! that we have to change.
      ! To calculate the matrix contributions, we have to evaluate
      ! the elements to give us the values of the basis functions
      ! in all the DOF`s in all the elements in our set.

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag. 
      cevaluationTag = rlocalMatrixAssembly%cevaluationTag
      
      ! In the first loop, calculate the coordinates on the reference element.
      ! In all later loops, use the precalculated information.
      !
      ! If the cubature points are already initialised, do not do it again.
      ! We check this by taking a look to iinitialisedElements which
      ! gives the current maximum of initialised elements.
      if (IELmax .gt. rlocalMatrixAssembly%iinitialisedElements) then

        ! (Re-)initialise!
        cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REFPOINTS)

        ! Remember the new number of initialised elements
        rlocalMatrixAssembly%iinitialisedElements = IELmax

      else
        ! No need.
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
      end if

      ! Calculate all information that is necessary to evaluate the finite element
      ! on all cells of our subset. This includes the coordinates of the points
      ! on the cells.
      call elprep_prepareSetForEvaluation (p_revalElementSet,&
          cevaluationTag, rmatrix%p_rspatialDiscrTest%p_rtriangulation, &
          IelementList(IELset:IELmax), rlocalMatrixAssembly%ctrafoType, &
          rlocalMatrixAssembly%p_DcubPtsRef(:,1:ncubp), rperfconfig=rperfconfig)
      p_Ddetj => p_revalElementSet%p_Ddetj
      
      ! If the matrix has nonconstant coefficients, calculate the coefficients now.
      if (.not. rlocalMatrixAssembly%rform%ballCoeffConstant) then
        if (present(fcoeff_buildMatrixSc_sim)) then
          call domint_initIntegrationByEvalSet (p_revalElementSet,rintSubset)
          !rintSubset%ielementDistribution =  0
          rintSubset%ielementStartIdx     =  IELset
          rintSubset%p_Ielements          => IelementList(IELset:IELmax)
          rintSubset%p_IdofsTrial         => p_IdofsTrial
          rintSubset%celement             =  rlocalMatrixAssembly%celementTrial
          call fcoeff_buildMatrixSc_sim (rmatrix%p_rspatialDiscrTest,&
              rmatrix%p_rspatialDiscrTrial,&
              rlocalMatrixAssembly%rform, IELmax-IELset+1, ncubp,&
              p_revalElementSet%p_DpointsReal(:,:,1:IELmax-IELset+1),&
              p_IdofsTrial, p_IdofsTest, rintSubset, &
              p_Dcoefficients(:,:,1:IELmax-IELset+1), rcollection)
          call domint_doneIntegration (rintSubset)
        else
          p_Dcoefficients(:,:,1:IELmax-IELset+1) = 1.0_DP
        end if
      end if
      
      ! Calculate the values of the basis functions.
      call elem_generic_sim2 (rlocalMatrixAssembly%celementTest, &
          p_revalElementSet, rlocalMatrixAssembly%BderTest, &
          rlocalMatrixAssembly%p_DbasTest)
      
      ! Omit the calculation of the trial function values if they
      ! are identical to the test function values.
      if (.not. rlocalMatrixAssembly%bidenticalTrialAndTest) then
        call elem_generic_sim2 (rlocalMatrixAssembly%celementTrial, &
            p_revalElementSet, rlocalMatrixAssembly%BderTrial, &
            rlocalMatrixAssembly%p_DbasTrial)
      end if
      
      ! --------------------- DOF COMBINATION PHASE ------------------------
      
      ! Values of all basis functions calculated. Now we can start
      ! to integrate!

      ! Clear the local matrices
      p_Dentry(:,:,1:IELmax-IELset+1) = 0.0_DP
      
      ! We have two different versions for the integration - one
      ! with constant coefficients and one with nonconstant coefficients.
      !
      ! Check the bilinear form which one to use:
      
      if (rlocalMatrixAssembly%rform%ballCoeffConstant) then
      
        ! Constant coefficients. The coefficients are to be found in
        ! the Dcoefficients variable of the form.
        !
        ! Loop over the elements in the current set.

        do iel = 1,IELmax-IELset+1
          
          ! Loop over all cubature points on the current element
          do icubp = 1, ncubp

            ! calculate the current weighting factor in the cubature formula
            ! in that cubature point.
            !
            ! Take the absolut value of the determinant of the mapping.
            ! In 2D, the determinant is always positive, whereas in 3D,
            ! the determinant might be negative -- that is normal!

            domega = p_Domega(icubp)*abs(p_Ddetj(icubp,iel))

            ! Loop over the additive factors in the bilinear form.
            do ialbet = 1,rlocalMatrixAssembly%rform%itermcount
            
              ! Get from Idescriptors the type of the derivatives for the 
              ! test and trial functions. The summand we calculate
              ! here will be added to the matrix entry:
              !
              ! a_ij  =  int_... ( psi_j )_ib  *  ( phi_i )_ia
              !
              ! -> Ix=0: function value, 
              !      =1: first derivative, ...
              !    as defined in the module 'derivative'.
              
              ia = p_Idescriptors(1,ialbet)
              ib = p_Idescriptors(2,ialbet)
              
              ! Multiply domega with the coefficient of the form.
              ! This gives the actual value to multiply the
              ! function value with before summing up to the integral.
              daux = domega * p_DcoefficientsBilf(ialbet)
            
              ! Now loop through all possible combinations of DOF`s
              ! in the current cubature point. The outer loop
              ! loops through the "O"`s in the above picture,
              ! the test functions:

              do idofe = 1,indofTest
              
                ! Get the value of the (test) basis function 
                ! phi_i (our "O") in the cubature point:
                db = p_DbasTest(idofe,ib,icubp,iel)
                
                ! Perform an inner loop through the other DOF`s
                ! (the "X"). 

                do jdofe = 1,indofTrial
                
                  ! Get the value of the basis function 
                  ! psi_j (our "X") in the cubature point. 
                  ! Them multiply:
                  !    db * dbas(..) * daux
                  ! ~= phi_i * psi_j * coefficient * cub.weight
                  ! Summing this up gives the integral, so the contribution
                  ! to the global matrix. 
                  !
                  ! Simply summing up db * dbas(..) * daux would give
                  ! the coefficient of the local matrix. We save this
                  ! contribution in the local matrix.

                  !JCOLB = Kentry(jdofe,idofe,iel)
                  !p_DA(JCOLB) = p_DA(JCOLB) + db*p_DbasTrial(jdofe,ia,icubp,iel)*daux
                  p_Dentry(jdofe,idofe,iel) = p_Dentry(jdofe,idofe,iel) + &
                                        db*p_DbasTrial(jdofe,ia,icubp,iel)*daux
                
                end do ! jdofe
              
              end do ! idofe
              
            end do ! ialbet

          end do ! icubp 
          
        end do ! iel
        
      else
      
        ! Nonconstant coefficients. The coefficients are to be found in
        ! the Dcoefficients variable as computed above.
        !
        ! Loop over the elements in the current set.

        do iel = 1,IELmax-IELset+1
          
          ! Loop over all cubature points on the current element
          do icubp = 1, ncubp

            ! calculate the current weighting factor in the cubature formula
            ! in that cubature point.
            !
            ! Take the absolut value of the determinant of the mapping.
            ! In 2D, the determinant is always positive, whereas in 3D,
            ! the determinant might be negative -- that is normal!

            domega = p_Domega(icubp)*abs(p_Ddetj(icubp,iel))

            ! Loop over the additive factors in the bilinear form.
            do ialbet = 1,rlocalMatrixAssembly%rform%itermcount
            
              ! Get from Idescriptors the type of the derivatives for the 
              ! test and trial functions. The summand we calculate
              ! here will be added to the matrix entry:
              !
              ! a_ij  =  int_... ( psi_j )_ia  *  ( phi_i )_ib
              !
              ! -> Ix=0: function value, 
              !      =1: first derivative, ...
              !    as defined in the module 'derivative'.
              
              ia = rlocalMatrixAssembly%rform%Idescriptors(1,ialbet)
              ib = rlocalMatrixAssembly%rform%Idescriptors(2,ialbet)
              
              ! Multiply domega with the coefficient of the form.
              ! This gives the actual value to multiply the
              ! function value with before summing up to the integral.
              ! Get the precalculated coefficient from the coefficient array.
              daux = domega * p_Dcoefficients(ialbet,icubp,iel)
            
              ! Now loop through all possible combinations of DOF`s
              ! in the current cubature point. The outer loop
              ! loops through the "O" in the above picture,
              ! the test functions:

              do idofe = 1,indofTest
                
                ! Get the value of the (test) basis function 
                ! phi_i (our "O") in the cubature point:
                db = p_DbasTest(idofe,ib,icubp,iel)
                
                ! Perform an inner loop through the other DOF`s
                ! (the "X"). 

                do jdofe = 1,indofTrial
              
                  ! Get the value of the basis function 
                  ! psi_j (our "X") in the cubature point. 
                  ! Them multiply:
                  !    db * dbas(..) * daux
                  ! ~= phi_i * psi_j * coefficient * cub.weight
                  ! Summing this up gives the integral, so the contribution
                  ! to the global matrix. 
                  !
                  ! Simply summing up db * dbas(..) * daux would give
                  ! the coefficient of the local matrix. We save this
                  ! contribution in the local matrix of element iel.

                  !JCOLB = Kentry(jdofe,idofe,iel)
                  !p_DA(JCOLB) = p_DA(JCOLB) + db*p_DbasTrial(jdofe,ia,icubp,iel)*daux
                  p_Dentry(jdofe,idofe,iel) = &
                      p_Dentry(jdofe,idofe,iel)+db*p_DbasTrial(jdofe,ia,icubp,iel)*daux
                
                end do
              
              end do ! jdofe
              
            end do ! ialbet

          end do ! icubp 
          
        end do ! iel

      end if ! rform%ballCoeffConstant

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
      
      ddgmass = 0.0_DP
      dtemass = 0.0_DP      
      do iel = 1,IELmax-IELset+1          
      
!         do idofe = 1,indofTest
!           do jdofe = 1,indofTrial
!             p_DA(p_Kentry(jdofe,idofe,iel)) = &
!                 p_DA(p_Kentry(jdofe,idofe,iel)) + p_Dentry(jdofe,idofe,iel)
!           end do
!         end do

        do idofe = 1,indofTest
          ddgmass = ddgmass + p_Dentry(idofe,idofe,iel)
          do jdofe = 1,indofTrial
            dtemass = dtemass + p_Dentry(jdofe,idofe,iel)
!            p_DA(p_Kentry(jdofe,idofe,iel)) = &
!                p_DA(p_Kentry(jdofe,idofe,iel)) + p_Dentry(jdofe,idofe,iel)
          end do
        end do

        if (ddgmass .eq. 0.0_DP) then
          dalpha = 1.0_DP
        else
          dalpha = dtemass/ddgmass
        end if
   
        do idofe = 1,indoftest
          p_DA(p_Kentry(idofe,idofe,iel)) = &
               p_DA(p_Kentry(idofe,idofe,iel)) + dalpha*p_Dentry(idofe,idofe,iel)
        end do       

      end do ! iel
      !%OMP END CRITICAL

    end do ! IELset
    
    ! Release the local matrix assembly structure
    call bilf_releaseAssemblyData(rlocalMatrixAssembly)
  
  end subroutine

end module