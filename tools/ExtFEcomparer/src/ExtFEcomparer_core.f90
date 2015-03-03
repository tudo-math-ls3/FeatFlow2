! This module contains the core of the compare-application
! It starts at the point where you have 2 FEM Solutions,
! created the according discretisations and parametrisations

module ExtFEcomparer_core

  use fsystem
  use storage
  use genoutput

  use element
  use cubature
  use linearsystemblock
  use linearsystemscalar

  use derivatives
  use collection
  use feevaluation2
  use feevaluation
  use triasearch

  use blockmatassembly
  use blockmatassemblybase
  use blockmatassemblystdop
  use spdiscprojection

  use ExtFEcomparer_typedefs

  use fparser
  use paramlist


implicit none

    private

    public :: ExtFEcomparer_calculate

contains


! This routine is the one that is called from outside
! In it, we prepare everything for the calculations
! that we want to do and then execute them
subroutine ExtFEcomparer_calculate(rproblem_1,rproblem_2, rpostprocessing)
!<input>
! 2 Problem structures that contain everything
  type(t_problem), intent(inout), target :: rproblem_1, rproblem_2
! </input>

!<inputoutput>
  ! The postprocessing structure that gets all results
  type(t_postprocessing), intent(inout) :: rpostprocessing
!</inputoutput>


  ! Do the Integral-Calculations if we want to do some
  ! I prefer the if here instead of an if-statement in
  ! the routine
  if (rpostprocessing%nIntCalculations .gt. 0 ) then
    call ExtFE_calc_Integral(rproblem_1,rproblem_2,rpostprocessing)
  end if

  ! Do the L1-Calculations if we want to do some
  ! I prefer the if here instead of an if-statement in
  ! the routine
  if (rpostprocessing%nL1Calculations .gt. 0 ) then
    call ExtFE_calc_L1(rproblem_1,rproblem_2,rpostprocessing)
  end if

  ! Do the L2-Calculations if we want to do some
  ! I prefer the if here instead of an if-statement in
  ! the routine
  if (rpostprocessing%nL2Calculations .gt. 0 ) then
    call ExtFE_calc_L2(rproblem_1,rproblem_2,rpostprocessing)
  end if

  ! Do the point-evaluations and calculations
  ! I prefer the if here instead of an if-statement in
  ! the routine
  if (rpostprocessing%nPointCalculations .gt. 0) then
    call ExtFE_calc_pointvalues(rproblem_1,rproblem_2,rpostprocessing)
  end if


  ! Do the computations neccesary for the UCD-Output

  call ExtFE_calc_UCD(rproblem_1,rproblem_2,rpostprocessing)


  ! Do the neccesary stuff for the vector output
  ! Is basically just setting some pointers
  call ExtFE_calc_OutOrigVec(rproblem_1,rproblem_2,rpostprocessing)


end subroutine


subroutine ExtFE_calc_Integral(rproblem_1,rproblem_2, rpostprocessing)

!<input>
! 2 Problem structures that contain everything
  type(t_problem), intent(inout), target :: rproblem_1, rproblem_2
! </input>

!<inputoutput>
! The postprocessing structure that receives all results
  type(t_postprocessing), intent(inout), target :: rpostprocessing
!</output>

  ! For the integral calculation we need a collection structure
  type(t_collection) :: collection

  ! Local variables

  ! We need to know how many calculations we want to do
  integer :: nIntCalculations

  ! pointer to the storage of the postprocessing structure
  ! where we will find what we want to calculate
  real(DP), dimension(:),  pointer :: p_Intresults =>NULL()
  integer, dimension (:,:), pointer :: p_IntComp => NULL()
  integer(I32), dimension(:), pointer :: p_IntCubRule => NULL()
  ! For the collection we need to convert the scalar
  ! vectors to block vectors
  type(t_vectorBlock), target :: blockSecond

  ! We need cubature information
  type(t_scalarCubatureInfo) :: RcubatureInformation

  ! We need a FE-Vector that is evaluated automatically
  ! during the integral calls
  type(t_fev2Vectors) :: fefunction

  ! We use this as a storage for the results
  real(DP) :: Dresult

  ! To figure out if we jump somewhere or not
  logical :: lOneFunction, lBothFunctions
  ! and pointers that we set so that we calculate the
  ! right things
  type(t_vectorScalar), pointer :: rCalcThis => NULL()
  type(t_vectorScalar), pointer :: rFirst => NULL()
  type(t_vectorScalar), pointer :: rSecond => NULL()


  ! We also need loop-variables
  integer :: i,k
  ! and strings for some local output so the user sees that
  ! calculations are going on
  character(LEN=ExtFE_STRLEN) :: soutput, stmp

  ! Find out how many L1-compuations to do
  nIntCalculations = rpostprocessing%nIntCalculations


  call output_lbrk()
  call output_line('Calculate the requested Integrals ...')
  call output_lbrk()

  ! Get the arrays
  call storage_getbase_double(rpostprocessing%h_IntResults,p_Intresults)
  call storage_getbase_int2D(rpostprocessing%h_IntCompFunc,p_IntComp)
  call storage_getbase_int(rpostprocessing%h_IntCubRule,p_IntCubRule)

  ! Attatch the parser to the collection
  ! It has already parsed the region of interest
  collection%p_rfparserQuickAccess1=>rpostprocessing%pIntChiOmegaParser

  ! We calculate the L1-Difference between L1comp(1,i) and L1comp(2,i)
  ! on L1ChiOmega(i) and use cubature rule L1CubRule(i)
  ! If one of the components is set to -1, we calculate the L1-Norm of
  ! the other one

   do i=1,nIntCalculations

      ! Init everything with false and empty in every loop
      ! so that the branching works
      lOneFunction = .FALSE.
      lBothFunctions = .FALSE.
      stmp = ''

      ! Find out what to calculate and set the pointers
      ! Pointers rFirst/rSecond points directly to the right component
      ! of the blockvectors if we calculate a difference
      ! rCalcThis points directly to the component of the first
      ! or second function if we want to calculate only a norm and not a difference
      ! First option: both are > 0 - so we want
      ! to calculate the difference
      if( (p_IntComp(1,i) .gt. 0) .and. &
          (p_Intcomp(2,i) .gt. 0)) then
          lOneFunction = .FALSE.
          lBothFunctions = .TRUE.
          rFirst => rproblem_1%coeffVector%RvectorBlock(p_IntComp(1,i))
          rSecond=> rproblem_2%coeffVector%RvectorBlock(p_IntComp(2,i))
          !write(stmp,'(A7,I2,A8,I2,A13,I2,A4)') '||f_1(', p_L1comp(1,i) ,') - f_2(', p_L1comp(2,i) ,')||_L1(Omega_', i, ') = '
          write(stmp,'(A11,I2,A8,I2,A11,I2,A3)') '\int{ f_1(', p_IntComp(1,i) , &
                                ') - f_2(', p_IntComp(2,i) ,') }_Omega_', i, ' = '

      ! Second option: id1 >0, id2 <=0
      ! => calculate ||id1||
      else if((p_IntComp(1,i) .gt. 0) .and. &
              (p_IntComp(2,i) .le. 0)) then
          lOneFunction = .TRUE.
          lBothFunctions = .FALSE.
          rCalcThis => rproblem_1%coeffVector%RvectorBlock(p_IntComp(1,i))
          !write(stmp,'(A7,I2,A13,I2,A4)') '||f_1(', p_L1comp(1,i) ,')||_L1(Omega_', i, ') = '
          write(stmp,'(A11,I2,A11,I2,A3)') '\int{ f_1(', p_IntComp(1,i) , &
                                ') }_Omega_', i, ' = '

      ! Third option: id1<=0, id2 >0
      ! => calculate ||id2||
      else if((p_IntComp(1,i) .le. 0) .and. &
              (p_IntComp(2,i) .gt. 0)) then
          lOneFunction = .TRUE.
          lBothFunctions = .FALSE.
          rCalcThis => rproblem_2%coeffVector%RvectorBlock(p_IntComp(2,i))
          !write(stmp,'(A7,I2,A13,I2,A4)') '||f_2(', p_L1comp(2,i) ,')||_L1(Omega_', i, ') = '
          write(stmp,'(A11,I2,A11,I2,A3)') '\int{ f_2(', p_IntComp(2,i) , &
                                ') }_Omega_', i, ' = '

      else
          call output_line('You need always one component id that &
                   &is greater than 0', &
                   OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_calculate_L1")
          call sys_halt()
      end if ! Which component

      ! Now pick up which region of interest
      ! We know where we stored this information:
      ! It is in the parser of the postprocessing structure
      ! We attatched that one to the collection already, we just
      ! need to tell it which component to use.
      ! In component i we have parsed the region of interest for
      ! calculation i
      collection%IquickAccess(1) = i

      ! Here we trigger the real computation
      ! Do we calculate a norm or a difference?
      ! Norm:
      if((lOneFunction .eqv. .TRUE.) .AND. &
         (lBothFunctions .eqv. .FALSE.)) then
         ! Create the cubature information
          call spdiscr_createDefCubStructure(rCalcThis%p_rspatialDiscr,&
                    RcubatureInformation, p_IntCubRule(i))
          ! mark the function to be evaluated - only function values
          call fev2_addVectorToEvalList(fefunction,rCalcThis,0)

          ! We need different computation
          ! routines for each dimension due to the evaluation
          ! of L1ChiOmega - thats the only difference.
          ! The difference is hidden in the parser but still is there
          ! Therefore we have to branch here again
          ! Store the result in dresult
          select case (rproblem_1%iDimension)
              case(ExtFE_NDIM1)
                call bma_buildIntegral(Dresult,BMA_CALC_STANDARD,calc_integral_onefunc_1D, &
                    rcollection=collection,revalVectors=fefunction,rcubatureInfo=RcubatureInformation)
              case(ExtFE_NDIM2)
                call bma_buildIntegral(Dresult,BMA_CALC_STANDARD,calc_integral_onefunc_2D, &
                    rcollection=collection,revalVectors=fefunction,rcubatureInfo=RcubatureInformation)
              case default
                call output_line('Dimension of your problem is &
                      &not supported', &
                      OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_calculate")
                call sys_halt()
          end select

          ! Store the result
          p_Intresults(i) = Dresult

          ! Release everything
          call fev2_releaseVectorList(fefunction)
          call spdiscr_releaseCubStructure(RcubatureInformation)
          Dresult = -1.0_DP
          rCalcThis  => NULL()
          rFirst => NULL()
          rSecond => NULL()
        ! Difference:
        else if ( (lOneFunction .eqv. .FALSE.) .AND. &
                 (lBothFunctions .eqv. .TRUE.)) then
          ! Create a cubature information
          call spdiscr_createDefCubStructure(rFirst%p_rspatialDiscr,&
                  RcubatureInformation, p_IntCubRule(i))

          ! mark the function to be evaluated
          call fev2_addVectorToEvalList(fefunction,rFirst,0)

          ! Now the tricky part: We cannot pass the second
          ! function directly since it lives on another mesh
          ! so we pass it via the collection
          ! The collection supports only block-vector as quick-access,
          ! so we have to convert it
          call lsysbl_createVecFromScalar(rSecond,blockSecond)
           collection%p_rvectorQuickAccess1 => blockSecond

          ! We need to make it dimension specific because of
          ! the evaluation of the second function + the evaluation
          ! of the region of interest
          select case (rproblem_1%iDimension)
             case(ExtFE_NDIM1)
               call bma_buildIntegral(Dresult,BMA_CALC_STANDARD,calc_L1error_twofunc_1D, &
                   rcollection=collection,revalVectors=fefunction,rcubatureInfo=RcubatureInformation)
             case(ExtFE_NDIM2)
               call bma_buildIntegral(Dresult,BMA_CALC_STANDARD,calc_L1error_twofunc_2D, &
                   rcollection=collection,revalVectors=fefunction,rcubatureInfo=RcubatureInformation)
             case default
               call output_line('Dimension of your problem is &
                     &not supported', &
                     OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_calculate")
               call sys_halt()
          end select


          ! Store the result
          p_Intresults(i) = Dresult
          ! Release everything
          call fev2_releaseVectorList(fefunction)
          call spdiscr_releaseCubStructure(RcubatureInformation)
          Dresult = -1.0_DP
          collection%p_rvectorQuickAccess1 => NULL()
          call lsysbl_releaseVector(blockSecond)
          rFirst => NULL()
          rSecond=>NULL()
          rCalcThis=>NULL()
        end if ! triggering L1 calculations

    ! Now generate some output on the terminal so that the user
    ! sees something is happening and can read the result
    ! Here, we do not write out any output in a file!
    ! That is done in the postprocessing!
    write(soutput,'(A1,E16.10)') ' ', p_Intresults(i)

    call output_line(trim(trim(stmp) // soutput), &
                    OU_CLASS_MSG,OU_MODE_STD+OU_MODE_BENCHLOG)

  end do ! Loop over all L1-Calculations

end subroutine


!---------------------------------------------------------------------

subroutine ExtFE_calc_L1(rproblem_1,rproblem_2, rpostprocessing)

!<input>
! 2 Problem structures that contain everything
  type(t_problem), intent(inout), target :: rproblem_1, rproblem_2
! </input>

!<inputoutput>
! The postprocessing structure that receives all results
  type(t_postprocessing), intent(inout), target :: rpostprocessing
!</output>

  ! For the integral calculation we need a collection structure
  type(t_collection) :: collection

  ! Local variables

  ! We need to know how many calculations we want to do
  integer :: nL1Calculations

  ! pointer to the storage of the postprocessing structure
  ! where we will find what we want to calculate
  real(DP), dimension(:),  pointer :: p_L1results =>NULL()
  integer, dimension (:,:), pointer :: p_L1Comp => NULL()
  integer(I32), dimension(:), pointer :: p_L1CubRule => NULL()
  ! For the collection we need to convert the scalar
  ! vectors to block vectors
  type(t_vectorBlock), target :: blockSecond

  ! We need cubature information
  type(t_scalarCubatureInfo) :: RcubatureInformation

  ! We need a FE-Vector that is evaluated automatically
  ! during the integral calls
  type(t_fev2Vectors) :: fefunction

  ! We use this as a storage for the results
  real(DP) :: Dresult

  ! To figure out if we jump somewhere or not
  logical :: lOneFunction, lBothFunctions
  ! and pointers that we set so that we calculate the
  ! right things
  type(t_vectorScalar), pointer :: rCalcThis => NULL()
  type(t_vectorScalar), pointer :: rFirst => NULL()
  type(t_vectorScalar), pointer :: rSecond => NULL()


  ! We also need loop-variables
  integer :: i,k
  ! and strings for some local output so the user sees that
  ! calculations are going on
  character(LEN=ExtFE_STRLEN) :: soutput, stmp

  ! Find out how many L1-compuations to do
  nL1Calculations = rpostprocessing%nL1Calculations


  call output_lbrk()
  call output_line('Calculate the requested L1-norms ...')
  call output_lbrk()

  ! Get the arrays
  call storage_getbase_double(rpostprocessing%h_L1Results,p_L1results)
  call storage_getbase_int2D(rpostprocessing%h_L1CompFunc,p_L1Comp)
  call storage_getbase_int(rpostprocessing%h_L1CubRule,p_L1CubRule)

  ! Attatch the parser to the collection
  ! It has already parsed the region of interest
  collection%p_rfparserQuickAccess1=>rpostprocessing%pL1ChiOmegaParser

  ! We calculate the L1-Difference between L1comp(1,i) and L1comp(2,i)
  ! on L1ChiOmega(i) and use cubature rule L1CubRule(i)
  ! If one of the components is set to -1, we calculate the L1-Norm of
  ! the other one

   do i=1,nL1Calculations

      ! Init everything with false and empty in every loop
      ! so that the branching works
      lOneFunction = .FALSE.
      lBothFunctions = .FALSE.
      stmp = ''

      ! Find out what to calculate and set the pointers
      ! Pointers rFirst/rSecond points directly to the right component
      ! of the blockvectors if we calculate a difference
      ! rCalcThis points directly to the component of the first
      ! or second function if we want to calculate only a norm and not a difference
      ! First option: both are > 0 - so we want
      ! to calculate the difference
      if( (p_L1Comp(1,i) .gt. 0) .and. &
          (p_L1comp(2,i) .gt. 0)) then
          lOneFunction = .FALSE.
          lBothFunctions = .TRUE.
          rFirst => rproblem_1%coeffVector%RvectorBlock(p_L1Comp(1,i))
          rSecond=> rproblem_2%coeffVector%RvectorBlock(p_L1Comp(2,i))
          write(stmp,'(A7,I2,A8,I2,A13,I2,A4)') '||f_1(', p_L1comp(1,i) ,') - f_2(', p_L1comp(2,i) ,')||_L1(Omega_', i, ') = '

      ! Second option: id1 >0, id2 <=0
      ! => calculate ||id1||
      else if((p_L1Comp(1,i) .gt. 0) .and. &
              (p_L1Comp(2,i) .le. 0)) then
          lOneFunction = .TRUE.
          lBothFunctions = .FALSE.
          rCalcThis => rproblem_1%coeffVector%RvectorBlock(p_L1Comp(1,i))
          write(stmp,'(A7,I2,A13,I2,A4)') '||f_1(', p_L1comp(1,i) ,')||_L1(Omega_', i, ') = '


      ! Third option: id1<=0, id2 >0
      ! => calculate ||id2||
      else if((p_L1Comp(1,i) .le. 0) .and. &
              (p_L1Comp(2,i) .gt. 0)) then
          lOneFunction = .TRUE.
          lBothFunctions = .FALSE.
          rCalcThis => rproblem_2%coeffVector%RvectorBlock(p_L1Comp(2,i))
          write(stmp,'(A7,I2,A13,I2,A4)') '||f_2(', p_L1comp(2,i) ,')||_L1(Omega_', i, ') = '

      else
          call output_line('You need always one component id that &
                   &is greater than 0', &
                   OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_calculate_L1")
          call sys_halt()
      end if ! Which component

      ! Now pick up which region of interest
      ! We know where we stored this information:
      ! It is in the parser of the postprocessing structure
      ! We attatched that one to the collection already, we just
      ! need to tell it which component to use.
      ! In component i we have parsed the region of interest for
      ! calculation i
      collection%IquickAccess(1) = i

      ! Here we trigger the real computation
      ! Do we calculate a norm or a difference?
      ! Norm:
      if((lOneFunction .eqv. .TRUE.) .AND. &
         (lBothFunctions .eqv. .FALSE.)) then
         ! Create the cubature information
          call spdiscr_createDefCubStructure(rCalcThis%p_rspatialDiscr,&
                    RcubatureInformation, p_L1CubRule(i))
          ! mark the function to be evaluated - only function values
          call fev2_addVectorToEvalList(fefunction,rCalcThis,0)

          ! We need different computation
          ! routines for each dimension due to the evaluation
          ! of L1ChiOmega - thats the only difference.
          ! The difference is hidden in the parser but still is there
          ! Therefore we have to branch here again
          ! Store the result in dresult
          select case (rproblem_1%iDimension)
              case(ExtFE_NDIM1)
                call bma_buildIntegral(Dresult,BMA_CALC_STANDARD,calc_L1norm_onefunc_1D, &
                    rcollection=collection,revalVectors=fefunction,rcubatureInfo=RcubatureInformation)
              case(ExtFE_NDIM2)
                call bma_buildIntegral(Dresult,BMA_CALC_STANDARD,calc_L1norm_onefunc_2D, &
                    rcollection=collection,revalVectors=fefunction,rcubatureInfo=RcubatureInformation)
              case default
                call output_line('Dimension of your problem is &
                      &not supported', &
                      OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_calculate")
                call sys_halt()
          end select

          ! Store the result
          p_L1results(i) = Dresult

          ! Release everything
          call fev2_releaseVectorList(fefunction)
          call spdiscr_releaseCubStructure(RcubatureInformation)
          Dresult = -1.0_DP
          rCalcThis  => NULL()
          rFirst => NULL()
          rSecond => NULL()
        ! Difference:
        else if ( (lOneFunction .eqv. .FALSE.) .AND. &
                 (lBothFunctions .eqv. .TRUE.)) then
          ! Create a cubature information
          call spdiscr_createDefCubStructure(rFirst%p_rspatialDiscr,&
                  RcubatureInformation, p_L1CubRule(i))

          ! mark the function to be evaluated
          call fev2_addVectorToEvalList(fefunction,rFirst,0)

          ! Now the tricky part: We cannot pass the second
          ! function directly since it lives on another mesh
          ! so we pass it via the collection
          ! The collection supports only block-vector as quick-access,
          ! so we have to convert it
          call lsysbl_createVecFromScalar(rSecond,blockSecond)
           collection%p_rvectorQuickAccess1 => blockSecond

          ! We need to make it dimension specific because of
          ! the evaluation of the second function + the evaluation
          ! of the region of interest
          select case (rproblem_1%iDimension)
             case(ExtFE_NDIM1)
               call bma_buildIntegral(Dresult,BMA_CALC_STANDARD,calc_L1error_twofunc_1D, &
                   rcollection=collection,revalVectors=fefunction,rcubatureInfo=RcubatureInformation)
             case(ExtFE_NDIM2)
               call bma_buildIntegral(Dresult,BMA_CALC_STANDARD,calc_L1error_twofunc_2D, &
                   rcollection=collection,revalVectors=fefunction,rcubatureInfo=RcubatureInformation)
             case default
               call output_line('Dimension of your problem is &
                     &not supported', &
                     OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_calculate")
               call sys_halt()
          end select


          ! Store the result
          p_L1results(i) = Dresult
          ! Release everything
          call fev2_releaseVectorList(fefunction)
          call spdiscr_releaseCubStructure(RcubatureInformation)
          Dresult = -1.0_DP
          collection%p_rvectorQuickAccess1 => NULL()
          call lsysbl_releaseVector(blockSecond)
          rFirst => NULL()
          rSecond=>NULL()
          rCalcThis=>NULL()
        end if ! triggering L1 calculations

    ! Now generate some output on the terminal so that the user
    ! sees something is happening and can read the result
    ! Here, we do not write out any output in a file!
    ! That is done in the postprocessing!
    write(soutput,'(A1,E16.10)') ' ', p_L1results(i)

    call output_line(trim(trim(stmp) // soutput), &
                    OU_CLASS_MSG,OU_MODE_STD+OU_MODE_BENCHLOG)

  end do ! Loop over all L1-Calculations

end subroutine


! ------L2 - calculations
subroutine ExtFE_calc_L2(rproblem_1,rproblem_2, rpostprocessing)

!<input>
! 2 Problem structures that contain everything
  type(t_problem), intent(inout), target :: rproblem_1, rproblem_2
! </input>

!<inputoutput>
! The postprocessing structure that receives all results
  type(t_postprocessing), intent(inout), target :: rpostprocessing
!</output>

  ! For the integral calculation we need a collection structure
  type(t_collection) :: collection

  ! Local variables

  ! We need to know how many calculations we want to do
  integer :: nL2Calculations

  ! pointer to the storage of the postprocessing structure
  ! where we will find what we want to calculate
  real(DP), dimension(:),  pointer :: p_L2results =>NULL()
  integer, dimension (:,:), pointer :: p_L2Comp => NULL()
  integer(I32), dimension(:), pointer :: p_L2CubRule => NULL()
  ! For the collection we need to convert the scalar
  ! vectors to block vectors
  type(t_vectorBlock), target :: blockSecond

  ! We need cubature information
  type(t_scalarCubatureInfo) :: RcubatureInformation

  ! We need a FE-Vector that is evaluated automatically
  ! during the integral calls
  type(t_fev2Vectors) :: fefunction

  ! We use this as a storage for the results
  real(DP) :: Dresult

  ! To figure out if we jump somewhere or not
  logical :: lOneFunction, lBothFunctions
  ! and pointers that we set so that we calculate the
  ! right things
  type(t_vectorScalar), pointer :: rCalcThis => NULL()
  type(t_vectorScalar), pointer :: rFirst => NULL()
  type(t_vectorScalar), pointer :: rSecond => NULL()


  ! We also need loop-variables
  integer :: i,k
  ! and strings for some local output so the user sees that
  ! calculations are going on
  character(LEN=ExtFE_STRLEN) :: soutput, stmp

  ! Find out how many L2-compuations to do
  nL2Calculations = rpostprocessing%nL2Calculations


  call output_lbrk()
  call output_line('Calculate the requested L2-norms ...')
  call output_lbrk()

  ! Get the arrays
  call storage_getbase_double(rpostprocessing%h_L2Results,p_L2results)
  call storage_getbase_int2D(rpostprocessing%h_L2CompFunc,p_L2Comp)
  call storage_getbase_int(rpostprocessing%h_L2CubRule,p_L2CubRule)

  ! Attatch the parser to the collection
  ! It has already parsed the region of interest
  collection%p_rfparserQuickAccess1=>rpostprocessing%pL2ChiOmegaParser

  ! We calculate the L2-Difference between L2comp(1,i) and L2comp(2,i)
  ! on L2ChiOmega(i) and use cubature rule L2CubRule(i)
  ! If one of the components is set to -1, we calculate the L2-Norm of
  ! the other one

   do i=1,nL2Calculations

      ! Init everything with false and empty in every loop
      ! so that the branching works
      lOneFunction = .FALSE.
      lBothFunctions = .FALSE.
      stmp = ''

      ! Find out what to calculate and set the pointers
      ! Pointers rFirst/rSecond points directly to the right component
      ! of the blockvectors if we calculate a difference
      ! rCalcThis points directly to the component of the first
      ! or second function if we want to calculate only a norm and not a difference
      ! First option: both are > 0 - so we want
      ! to calculate the difference
      if( (p_L2Comp(1,i) .gt. 0) .and. &
          (p_L2comp(2,i) .gt. 0)) then
          lOneFunction = .FALSE.
          lBothFunctions = .TRUE.
          rFirst => rproblem_1%coeffVector%RvectorBlock(p_L2Comp(1,i))
          rSecond=> rproblem_2%coeffVector%RvectorBlock(p_L2Comp(2,i))
          write(stmp,'(A7,I2,A8,I2,A13,I2,A4)') '||f_1(', p_L2comp(1,i) ,') - f_2(', p_L2comp(2,i) ,')||_L2(Omega_', i, ') = '

      ! Second option: id1 >0, id2 <=0
      ! => calculate ||id1||
      else if((p_L2Comp(1,i) .gt. 0) .and. &
              (p_L2Comp(2,i) .le. 0)) then
          lOneFunction = .TRUE.
          lBothFunctions = .FALSE.
          rCalcThis => rproblem_1%coeffVector%RvectorBlock(p_L2Comp(1,i))
          write(stmp,'(A7,I2,A13,I2,A4)') '||f_1(', p_L2comp(1,i) ,')||_L2(Omega_', i, ') = '


      ! Third option: id1<=0, id2 >0
      ! => calculate ||id2||
      else if((p_L2Comp(1,i) .le. 0) .and. &
              (p_L2Comp(2,i) .gt. 0)) then
          lOneFunction = .TRUE.
          lBothFunctions = .FALSE.
          rCalcThis => rproblem_2%coeffVector%RvectorBlock(p_L2Comp(2,i))
          write(stmp,'(A7,I2,A13,I2,A4)') '||f_2(', p_L2comp(2,i) ,')||_L2(Omega_', i, ') = '

      else
          call output_line('You need always one component id that &
                   &is greater than 0', &
                   OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_calculate_L2")
          call sys_halt()
      end if ! Which component

      ! Now pick up which region of interest
      ! We know where we stored this information:
      ! It is in the parser of the postprocessing structure
      ! We attatched that one to the collection already, we just
      ! need to tell it which component to use.
      ! In component i we have parsed the region of interest for
      ! calculation i
      collection%IquickAccess(1) = i

      ! Here we trigger the real computation
      ! Do we calculate a norm or a difference?
      ! Norm:
      if((lOneFunction .eqv. .TRUE.) .AND. &
         (lBothFunctions .eqv. .FALSE.)) then
         ! Create the cubature information
          call spdiscr_createDefCubStructure(rCalcThis%p_rspatialDiscr,&
                    RcubatureInformation, p_L2CubRule(i))
          ! mark the function to be evaluated - only function values
          call fev2_addVectorToEvalList(fefunction,rCalcThis,0)

          ! We need different computation
          ! routines for each dimension due to the evaluation
          ! of L2ChiOmega - thats the only difference.
          ! The difference is hidden in the parser but still is there
          ! Therefore we have to branch here again
          ! Store the result in dresult
          select case (rproblem_1%iDimension)
              case(ExtFE_NDIM1)
                call bma_buildIntegral(Dresult,BMA_CALC_STANDARD,calc_L2norm_onefunc_1D, &
                    rcollection=collection,revalVectors=fefunction,rcubatureInfo=RcubatureInformation)
              case(ExtFE_NDIM2)
                call bma_buildIntegral(Dresult,BMA_CALC_STANDARD,calc_L2norm_onefunc_2D, &
                    rcollection=collection,revalVectors=fefunction,rcubatureInfo=RcubatureInformation)
              case(ExtFE_NDIM3)
                call bma_buildIntegral(Dresult,BMA_CALC_STANDARD,calc_L2norm_onefunc_3D, &
                    rcollection=collection,revalVectors=fefunction,rcubatureInfo=RcubatureInformation)
              case default
                call output_line('Dimension of your problem is &
                      &not supported', &
                      OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_calculate")
                call sys_halt()
          end select

          ! The routine always returns the squared value, so
          ! take the square-root and store it
          p_L2results(i) = sqrt(Dresult)

          ! Release everything
          call fev2_releaseVectorList(fefunction)
          call spdiscr_releaseCubStructure(RcubatureInformation)
          Dresult = -1.0_DP
          rCalcThis  => NULL()
          rFirst => NULL()
          rSecond => NULL()
        ! Difference:
        else if ( (lOneFunction .eqv. .FALSE.) .AND. &
                 (lBothFunctions .eqv. .TRUE.)) then
          ! Create a cubature information
          call spdiscr_createDefCubStructure(rFirst%p_rspatialDiscr,&
                  RcubatureInformation, p_L2CubRule(i))

          ! mark the function to be evaluated
          call fev2_addVectorToEvalList(fefunction,rFirst,0)

          ! Now the tricky part: We cannot pass the second
          ! function directly since it lives on another mesh
          ! so we pass it via the collection
          ! The collection supports only block-vector as quick-access,
          ! so we have to convert it
          call lsysbl_createVecFromScalar(rSecond,blockSecond)
           collection%p_rvectorQuickAccess1 => blockSecond

          ! We need to make it dimension specific because of
          ! the evaluation of the second function + the evaluation
          ! of the region of interest
          select case (rproblem_1%iDimension)
             case(ExtFE_NDIM1)
               call bma_buildIntegral(Dresult,BMA_CALC_STANDARD,calc_L2error_twofunc_1D, &
                   rcollection=collection,revalVectors=fefunction,rcubatureInfo=RcubatureInformation)
             case(ExtFE_NDIM2)
               call bma_buildIntegral(Dresult,BMA_CALC_STANDARD,calc_L2error_twofunc_2D, &
                   rcollection=collection,revalVectors=fefunction,rcubatureInfo=RcubatureInformation)
             case default
               call output_line('Dimension of your problem is &
                     &not supported', &
                     OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_calculate")
               call sys_halt()
          end select


          ! The routine always returns the squared value, so
          ! take the square-root and store it
          p_L2results(i) = sqrt(Dresult)
          ! Release everything
          call fev2_releaseVectorList(fefunction)
          call spdiscr_releaseCubStructure(RcubatureInformation)
          Dresult = -1.0_DP
          collection%p_rvectorQuickAccess1 => NULL()
          call lsysbl_releaseVector(blockSecond)
          rFirst => NULL()
          rSecond=>NULL()
          rCalcThis=>NULL()
        end if ! triggering L2 calculations

    ! Now generate some output on the terminal so that the user
    ! sees something is happening and can read the result
    ! Here, we do not write out any output in a file!
    ! That is done in the postprocessing!
    write(soutput,'(A1,E16.10)') ' ', p_L2results(i)

    call output_line(trim(trim(stmp) // soutput), &
                    OU_CLASS_MSG,OU_MODE_STD+OU_MODE_BENCHLOG)

  end do ! Loop over all L2-Calculations

end subroutine

!---------------------------------------------------------------------

subroutine ExtFE_calc_pointvalues(rproblem_1,rproblem_2, rpostprocessing)

!<input>
! 2 Problem structures that contain everything
  type(t_problem), intent(inout) :: rproblem_1, rproblem_2
! </input>


!<inputoutput>
! The postprocessing structure that receives all results
  type(t_postprocessing), intent(inout) :: rpostprocessing
!</output>

! Local variables

  ! Pointer to the arrays
  real(DP), dimension(:), pointer :: p_PointResults => NULL()
  real(DP), dimension(:,:), pointer :: p_PointCoords => NULL()
  integer, dimension(:,:), pointer :: p_PointFuncComp => NULL()

  integer :: nPointEvaluations,iDeriv
  real(DP), dimension(:), allocatable :: DEvalPoint

  ! To figure out if we jump somewhere or not
  logical :: lOneFunction, lBothFunctions
  ! Pointers to point to the function + comp. we want to calc
  type(t_vectorScalar), pointer :: rCalcThis => NULL()
  type(t_vectorScalar), pointer :: rFirst => NULL()
  type(t_vectorScalar), pointer :: rSecond => NULL()

  real(DP) :: Dresult
  ! We also need loop-variables
  integer :: i,k

  ! and strings for local output
  character(LEN=ExtFE_STRLEN) :: sparam
  character(LEN=ExtFE_STRLEN) :: soutput, stmp, stmp2

  ! For the output we need the names of the derivatives
  character(LEN=6), dimension(ExtFE_NumDerivs) :: sderivnames
  ! We start to count the derivatives from 0, so
  ! we need to shift by 1
  sderivnames(ExtFE_DER_FUNC_3D+1) = '      '
  sderivnames(ExtFE_DER_DERIV_3D_X+1) = '(d/dx)'
  sderivnames(ExtFE_DER_DERIV_3D_Y+1) = '(d/dy)'
  sderivnames(ExtFE_DER_DERIV_3D_Z+1) = '(d/dz)'

  nPointEvaluations = rpostprocessing%nPointCalculations

  call output_lbrk()
  call output_line('Calculate the requested pointvalues ...')
  call output_lbrk()

  ! get the arrays with the coordinates, the results and the components
  ! as we set it up in the init of the postprocessing
  call storage_getbase_double2D(rpostprocessing%h_PointCoordinates,p_PointCoords)
  call storage_getbase_double(rpostprocessing%h_PointResults,p_PointResults)
  call storage_getbase_int2D(rpostprocessing%h_PointFuncComponents,p_PointFuncComp)

  ! Allocate the evaluationpoint
  ! Needs to be done as we do some wrapping and call several
  ! routines. This way it is easier to read
  allocate(DEvalPoint(rproblem_1%iDimension))

  ! Now prepare and then trigger the
  ! evaluations
  do i=1,nPointEvaluations
      lOneFunction = .false.
      lBothFunctions = .false.
      iDeriv = -1

      ! We calculate a difference, so set both pointers
      ! and prepare local output string: (deriv (f_1(comp1) - deriv f_2(comp2))
      ! We cannot store the derivative ID in this case since we would need 2
      if( (p_PointFuncComp(1,i) .gt. 0) .and. &
          (p_PointFuncComp(3,i) .gt. 0)) then
          lOneFunction = .FALSE.
          lBothFunctions = .TRUE.
          rFirst => rproblem_1%coeffVector%RvectorBlock(p_PointFuncComp(1,i))
          rSecond=> rproblem_2%coeffVector%RvectorBlock(p_PointFuncComp(3,i))
          write(stmp,'(A3,I2,A3)') 'f1_', p_PointFuncComp(1,i), ' - '
          write(stmp2,'(A3,I2)') 'f2_', p_PointFuncComp(3,i)
          stmp = trim(sderivnames(p_PointFuncComp(2,i)+1)) // trim(stmp)
          stmp = '(' // trim(stmp)
          stmp2 = trim(sderivnames(p_PointFuncComp(4,i)+1)) // trim(stmp2)
          stmp2 = trim(stmp2) // ')'
          stmp = trim(stmp) // trim(stmp2)

      ! We evaluate only 1 function. Set the rCalcThis pointer
      ! and prepare local output string: (deriv (f_1(comp1)))
      ! Store the local derivative Identifier in iDeriv
      else if((p_PointFuncComp(1,i) .gt. 0) .and. &
              (p_PointFuncComp(3,i) .le. 0)) then
          lOneFunction = .TRUE.
          lBothFunctions = .FALSE.
          rCalcThis => rproblem_1%coeffVector%RvectorBlock(p_PointFuncComp(1,i))
          iDeriv = p_PointFuncComp(2,i)
          write(stmp,'(A3,I2)') 'f1_', p_PointFuncComp(1,i)
          stmp = (sderivnames(p_PointFuncComp(2,i)+1)) // trim(stmp)
          stmp = '(' // trim(stmp)
          stmp = trim(stmp) // ')'
      ! We evaluate only 1 function. Set the rCalcThis pointer
      ! and prepare local output string: (deriv (f_2(comp2)))
      ! Store the local derivative ID in iDeriv
      else if((p_PointFuncComp(1,i) .le. 0) .and. &
              (p_PointFuncComp(3,i) .gt. 0)) then
          lOneFunction = .TRUE.
          lBothFunctions = .FALSE.
          rCalcThis => rproblem_2%coeffVector%RvectorBlock(p_PointFuncComp(3,i))
          iDeriv = p_PointFuncComp(4,i)
          write(stmp,'(A3,I2)') 'f2_', p_PointFuncComp(3,i)
          stmp = (sderivnames(p_PointFuncComp(4,i)+1)) // trim(stmp)
          stmp = '(' // trim(stmp)
          stmp = trim(stmp) // ')'
      else
          call output_line('You need always one component id that &
                  &is greater than 0', &
                  OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_calculate_pointvalues")
          call sys_halt()
      end if ! Which component


      ! Now we can do the calculations
      ! Put the point in the right shape
      do k=1,rproblem_1%iDimension
           DEvalPoint(k) = p_PointCoords(k,i)
      end do

      ! and trigger the computation
      ! Only 1 function:
      if((lOneFunction .eqv. .true.) .AND. &
         (lBothFunctions .eqv. .false.)) then

         ! Call the wrapper that translates the ExtFE-Variables to
         ! global featvariables and then triggers the real computation
         ! We are in the evaluation case, so we have the derivative we want
         ! in iDeriv
         call ExtFE_eval_function(rCalcThis,DEvalPoint,iDeriv,Dresult)

         ! Save the result in the array
         p_PointResults(i) = Dresult

      ! Both functions -> calc f1 - f2
      else if((lOneFunction .eqv. .false.) .AND. &
              (lBothFunctions .eqv. .true.)) then
          ! Store the derivative ID of the first function in iDeriv
          iDeriv = p_PointFuncComp(2,i)
          ! Call the wrapper with the first function
          call ExtFE_eval_function(rFirst,DEvalPoint,iDeriv,Dresult)
          ! And save it in the array
          p_PointResults(i) = Dresult

          ! Store the derivative ID of the second function in iDeriv
          iDeriv = p_PointFuncComp(4,i)
          ! call the wrapper with the second function
          call ExtFE_eval_function(rSecond,DEvalPoint,iDeriv,Dresult)
          ! and subtract the result from p_PointResults(i)
          p_PointResults(i) = p_PointResults(i) - Dresult

      end if

      ! We still have to add the point coordinates and the result
      ! to our local output string:
      stmp = trim(stmp) // '('
      write(stmp2,'(F8.4)') p_PointCoords(1,i)
      stmp = trim(stmp) // trim(stmp2)
      do k=2,rproblem_1%iDimension
          write(stmp2,'(A1,F8.4)') ',', p_PointCoords(k,i)
          stmp = trim(stmp) // trim(stmp2)
          stmp = stmp // ' '
      end do
      stmp = trim(stmp) // ' ) '
      write(stmp2,'(F20.10)') p_PointResults(i)
      soutput = trim(stmp) // ' = ' // trim(adjustl(stmp2))
      ! Make the output
      call output_line(soutput,OU_CLASS_MSG,OU_MODE_STD+OU_MODE_BENCHLOG)


    end do ! all point calculations

    ! clean up
    deallocate(DEvalPoint)

end subroutine



! Do the computations/preparations necessary for the UCD-Output
! This includes: All projections that are necessary and setting
! of pointers

subroutine ExtFE_calc_UCD(rproblem_1,rproblem_2, rpostprocessing)

!<input>
! 2 Problem structures that contain everything
  type(t_problem), intent(inout), target :: rproblem_1, rproblem_2
! </input>

!<inputoutput>
! The Postprocessing structure that receives all results
  type(t_postprocessing), intent(inout) :: rpostprocessing
!</output>

  ! Local variables
  integer :: nDim, nVar,i
  integer, dimension(:), pointer :: p_IntPointerElemProject
  integer(I32) :: ID_ProjectConst, ID_ProjectLin
  real(DP), dimension(:), pointer :: p_VecEntries


  !----------------------------------------------!
  ! Prepare the UCD-Output: Calculate everything !
  ! and store it in the postprocessing structure !
  !----------------------------------------------!
  nDim = rproblem_1%rdiscretisation%ndimension

  ! The element on that we want to project.
  ! The id depends on the dimension
  select case(nDim)
    case(ExtFE_NDIM1)
        ID_ProjectConst = EL_P0_1D
        ID_ProjectLin = EL_P1_1D
    case(ExtFE_NDIM2)
        ID_ProjectConst = EL_Q0_2D
        ID_ProjectLin = EL_Q1_2D
    case(ExtFE_NDIM3)
        ID_ProjectConst = EL_Q0_3D
        ID_ProjectLin = EL_Q1_3D
    end select

  ! If we want to write out the meshes, we need to set the pointers in the postprocessing
  ! structure to the triangulations of the functions
  if(rpostprocessing%ucd_OUT_meshes .eqv. .true.) then
    rpostprocessing%UCD_MeshOnePointer => rproblem_1%coeffVector%p_rblockDiscr%p_rtriangulation
    rpostprocessing%UCD_MeshTwoPointer => rproblem_2%coeffVector%p_rblockDiscr%p_rtriangulation
  end if

  ! Write out first function?
  if(rpostprocessing%ucd_OUT_orig_functions_one .eqv. .true. ) then

    nVar = rproblem_1%coeffVector%nblocks
    ! We need a block-discretisation so that we can do the projection
    allocate(rpostprocessing%UCDBlockDiscrFirst)
    call spdiscr_initBlockDiscr(rpostprocessing%UCDBlockDiscrFirst,nVar,&
                    rproblem_1%coeffVector%p_rblockDiscr%p_rtriangulation)

    ! Now lets see how which variable we want to project in which way and
    ! create an according discretisation in the according block
    call storage_getbase_int(rpostprocessing%h_UCD_AddElemProjectFirst,p_IntPointerElemProject)

    do i=1,nVar
        select case(p_IntPointerElemProject(i))

        case(ExtFE_UCD_POLY_CONST)
            call spdiscr_initDiscr_simple(rpostprocessing%UCDBlockDiscrFirst%RspatialDiscr(i),&
                                    ID_ProjectConst,&
                                    rpostprocessing%UCDBlockDiscrFirst%p_rtriangulation)
        case(ExtFE_UCD_POLY_LINEAR)
            call spdiscr_initDiscr_simple(rpostprocessing%UCDBlockDiscrFirst%RspatialDiscr(i),&
                                    ID_ProjectLin,&
                                    rpostprocessing%UCDBlockDiscrFirst%p_rtriangulation)
        end select
    end do

    ! Now init the projection vector by the discretisation
    allocate(rpostprocessing%UCD_feFunction_first_orig)
    call lsysbl_createVector(rpostprocessing%UCDBlockDiscrFirst,rpostprocessing%UCD_feFunction_first_orig)

    ! Project the solution to the UCD_OUT_Vector
    do i=1,nVar
        call lsyssc_getbase_double(rpostprocessing%UCD_feFunction_first_orig%RvectorBlock(i),p_VecEntries)

        select case(p_IntPointerElemProject(i))

        case(ExtFE_UCD_POLY_CONST)
            call spdp_projectToCells(rproblem_1%coeffVector%RvectorBlock(i),p_VecEntries)
        case(ExtFE_UCD_POLY_LINEAR)
            call spdp_projectToVertices(rproblem_1%coeffVector%RvectorBlock(i),p_VecEntries)
        end select
        p_VecEntries => NULL()
    end do
    p_IntPointerElemProject => NULL()
  end if ! UCD_OUT_First_Function = TRUE


  ! Write out second function?
  if(rpostprocessing%ucd_OUT_orig_functions_two .eqv. .true. ) then
    nVar = rproblem_2%coeffVector%nblocks
    ! We need a block-discretisation so that we can do the projection
    allocate(rpostprocessing%UCDBlockDiscrSecond)
    call spdiscr_initBlockDiscr(rpostprocessing%UCDBlockDiscrSecond,nVar,&
                    rproblem_2%coeffVector%p_rblockDiscr%p_rtriangulation)

    ! Now lets see how which variable we want to project in which way and
    ! create an according discretisation in the according block

    call storage_getbase_int(rpostprocessing%h_UCD_AddElemProjectSecond,p_IntPointerElemProject)

    do i=1,nVar
        select case(p_IntPointerElemProject(i))

        case(ExtFE_UCD_POLY_CONST)
            call spdiscr_initDiscr_simple(rpostprocessing%UCDBlockDiscrSecond%RspatialDiscr(i),&
                                    ID_ProjectConst,&
                                    rpostprocessing%UCDBlockDiscrSecond%p_rtriangulation)
        case(ExtFE_UCD_POLY_LINEAR)
            call spdiscr_initDiscr_simple(rpostprocessing%UCDBlockDiscrSecond%RspatialDiscr(i),&
                                    ID_ProjectLin,&
                                    rpostprocessing%UCDBlockDiscrSecond%p_rtriangulation)
        end select
    end do

    ! Now init the projection vector by the discretisation
    allocate(rpostprocessing%UCD_feFunction_second_orig)
    call lsysbl_createVector(rpostprocessing%UCDBlockDiscrSecond,rpostprocessing%UCD_feFunction_second_orig)


    ! Project the solution to the UCD_OUT_Vector
    do i=1,nVar
        call lsyssc_getbase_double(rpostprocessing%UCD_feFunction_second_orig%RvectorBlock(i),p_VecEntries)

        select case(p_IntPointerElemProject(i))

        case(ExtFE_UCD_POLY_CONST)
            call spdp_projectToCells(rproblem_2%coeffVector%RvectorBlock(i),p_VecEntries)
        case(ExtFE_UCD_POLY_LINEAR)
            call spdp_projectToVertices(rproblem_2%coeffVector%RvectorBlock(i),p_VecEntries)
        end select
        p_VecEntries => NULL()
    end do

    p_IntPointerElemProject => NULL()

  end if ! UCD_OUT_Second_Function = TRUE


end subroutine


! Do the preparations neccesary for the Vector-Output
! Sounds weird,  but thats the structure of the code.
! Here we set up everything and store it in the postprocessing structure

subroutine ExtFE_calc_OutOrigVec(rproblem_1,rproblem_2, rpostprocessing)

!<input>
! 2 Problem structures that contain everything
  type(t_problem), intent(inout), target :: rproblem_1, rproblem_2
! </input>

!<inputoutput>
! The postprocessing structure that receives everything we do
  type(t_postprocessing), intent(inout) :: rpostprocessing
!</output>

  ! All we need to do is set some pointers if we want to write them out

  if(rpostprocessing%writeOutOrigVector1 .eqv. .true.) then
    rpostprocessing%OrigVecFirst => rproblem_1%coeffVector
  end if

  if(rpostprocessing%writeOutOrigVector2 .eqv. .true.) then
    rpostprocessing%OrigVecSecond => rproblem_2%coeffVector
  end if


end subroutine

!-----------------------------------------!
! Here, the real working routines follow: !
!-----------------------------------------!

!---------------------------------------------------------------------------!
!                                                                           !
!                   1D integral routines come here                          !
!                                                                           !
!---------------------------------------------------------------------------!

! Calculates the L1-Norm of a 1D function
! on a Domain specified in the *.dat-files
subroutine calc_integral_onefunc_1D(Dintvalue,rassemblyData,rintAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)
    !<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
  !</input>

!<output>
    ! Returns the value of the integral(s)
    real(DP), dimension(:), intent(out) :: Dintvalue
!</output>

    ! Local variables
    real(DP) :: dx,dval1, work


    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight => NULL()
    real(DP), dimension(:,:), pointer :: p_Dfunc1 => NULL()
    real(DP), dimension(:,:,:), pointer :: p_Dpoints => NULL()
    integer, dimension(:), pointer :: p_Ielements => NULL()

    ! Get cubature weights
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Element numbers of the elements in process
    p_Ielements => rassemblyData%p_IelementList

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal

    ! Get a pointer to the values
    p_Dfunc1 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC)

    ! We should have everything prepared, so let us start
    ! Initial value of the integral
    Dintvalue(1) = 0.0_DP

    ! Loop over all elements
    do iel = 1,nelements

        ! Loop over all cubature points
        do icubp = 1,npointsPerElement

          ! Get the coordinates of the cubature point
          dx = p_Dpoints(1,icubp,iel)

          ! Now we check if we need to work or not
          ! If we need to work, the domainDescription
          ! will return a 1, if not it will return a 0
          ! We saved this one in the collection structure already
          ! and parsed it

           call fparser_evalFunction(rcollection%p_rfparserQuickAccess1,&
                        rcollection%IquickAccess(1),(/dx/),work)

           ! Check if we need to work.
           ! If "work" is 1, we need to work. Since "work"
           ! can only be 0 or 1, this should run. However, you could
           ! hack this by defining functions that do not return 0 or 1
           ! as a function value, but. i.e. 0.5.
           if (abs(work - 1.0_DP) .le. 0.01_DP) then
                ! We need to work
                dval1 = p_Dfunc1(icubp,iel)
                Dintvalue(1) = Dintvalue(1) + dval1*p_DcubWeight(icubp,iel)
          end if


        end do ! icubp

    end do ! iel

end subroutine

!*********************************************************************************************
! Calculates the integral (f1 - f2) dx of 2 1D functions
! on a Domain specified in the *.dat-files

subroutine calc_integralDifference_twofunc_1D(Dintvalue,rassemblyData,rintAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)
    !<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
  !</input>

!<output>
    ! Returns the value of the integral(s)
    real(DP), dimension(:), intent(out) :: Dintvalue
!</output>

    ! Local variables
    real(DP) :: dx,dval1, work


    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight => NULL()
    real(DP), dimension(:,:), pointer :: p_Dfunc1 => NULL()
    real(DP), dimension(:,:,:), pointer :: p_Dpoints => NULL()
    integer, dimension(:), pointer :: p_Ielements => NULL()

    type(t_vectorScalar), pointer :: secondFunction => NULL()
    real(DP), dimension(ExtFE_NDIM1,1) :: evaluationPoint
    real(DP), dimension(1) :: dFuncValue
    integer :: ielEval
    integer, dimension(1) :: iElement
    real(DP) , dimension (1) :: dpoint

    ! Get cubature weights
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Element numbers of the elements in process
    p_Ielements => rassemblyData%p_IelementList

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal

    ! Get a pointer to the values
    p_Dfunc1 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC1D)

    ! Get the second function
    secondFunction => rcollection%p_rvectorQuickAccess1%RvectorBlock(1)

    ! We should have everything prepared, so let us start
    ! Initial value of the integral
    Dintvalue(1) = 0.0_DP

    ! Loop over all elements
    do iel = 1,nelements

        ! Loop over all cubature points
        do icubp = 1,npointsPerElement

          ! Get the coordinates of the cubature point
          dx = p_Dpoints(1,icubp,iel)

          ! Now we check if we need to work or not
          ! If we need to work, the domainDescription
          ! will return a 1, if not it will return a 0

           call fparser_evalFunction(rcollection%p_rfparserQuickAccess1, &
                    rcollection%IquickAccess(1),(/dx/),work)

           ! Check if we need to work.
           ! If "work" is 1, we need to work. Since "work"
           ! can only be 0 or 1, this should run. However, you could
           ! hack this by defining functions that do not return 0 or 1
           ! as a function value, but. i.e. 0.5.
           if (abs(work - 1.0_DP) .le. 0.01_DP) then
                ! We need to work
                ! This one is getting really ugly
                ! Make the evaluation point as an array
                evaluationPoint(1,1) = dx
                ! Find the element nr.
                ! Usually the evaluation routines do this on their
                ! own - but they support only 2D-search algorithms
                ! So we search the element nr. on our own - because
                ! if we support the element nr. by ourselves, the
                ! evaluation routine does not branch into the searching part
                ! Therefore we need some preparation
                dpoint(1) = dx
                call tsrch_getElem_BruteForce(dpoint,&
                    secondFunction%p_rspatialDiscr%p_rtriangulation,ielEval)
                iElement(1) = ielEval

                ! Now we have the element number so we can call the
                ! evaluation routine
                call fevl_evaluate(DER_FUNC1D,dFuncValue,secondFunction, &
                    evaluationPoint,Ielements=iElement)
                ! Calculate the integral value
                dval1 = p_Dfunc1(icubp,iel) - dFuncValue(1)
                Dintvalue(1) = Dintvalue(1) + dval1*p_DcubWeight(icubp,iel)
          end if


        end do ! icubp

    end do ! iel


end subroutine


! Calculates the L1-Norm of a 1D function
! on a Domain specified in the *.dat-files
subroutine calc_L1norm_onefunc_1D(Dintvalue,rassemblyData,rintAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)
    !<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
  !</input>

!<output>
    ! Returns the value of the integral(s)
    real(DP), dimension(:), intent(out) :: Dintvalue
!</output>

    ! Local variables
    real(DP) :: dx,dval1, work


    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight => NULL()
    real(DP), dimension(:,:), pointer :: p_Dfunc1 => NULL()
    real(DP), dimension(:,:,:), pointer :: p_Dpoints => NULL()
    integer, dimension(:), pointer :: p_Ielements => NULL()

    ! Get cubature weights
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Element numbers of the elements in process
    p_Ielements => rassemblyData%p_IelementList

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal

    ! Get a pointer to the values
    p_Dfunc1 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC)

    ! We should have everything prepared, so let us start
    ! Initial value of the integral
    Dintvalue(1) = 0.0_DP

    ! Loop over all elements
    do iel = 1,nelements

        ! Loop over all cubature points
        do icubp = 1,npointsPerElement

          ! Get the coordinates of the cubature point
          dx = p_Dpoints(1,icubp,iel)

          ! Now we check if we need to work or not
          ! If we need to work, the domainDescription
          ! will return a 1, if not it will return a 0
          ! We saved this one in the collection structure already
          ! and parsed it

           call fparser_evalFunction(rcollection%p_rfparserQuickAccess1,&
                        rcollection%IquickAccess(1),(/dx/),work)

           ! Check if we need to work.
           ! If "work" is 1, we need to work. Since "work"
           ! can only be 0 or 1, this should run. However, you could
           ! hack this by defining functions that do not return 0 or 1
           ! as a function value, but. i.e. 0.5.
           if (abs(work - 1.0_DP) .le. 0.01_DP) then
                ! We need to work
                dval1 = abs(p_Dfunc1(icubp,iel))
                Dintvalue(1) = Dintvalue(1) + dval1*p_DcubWeight(icubp,iel)
          end if


        end do ! icubp

    end do ! iel

end subroutine

!*********************************************************************************************
! Calculates the L1-distance of 2 1D functions
! on a Domain specified in the *.dat-files

subroutine calc_L1error_twofunc_1D(Dintvalue,rassemblyData,rintAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)
    !<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
  !</input>

!<output>
    ! Returns the value of the integral(s)
    real(DP), dimension(:), intent(out) :: Dintvalue
!</output>

    ! Local variables
    real(DP) :: dx,dval1, work


    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight => NULL()
    real(DP), dimension(:,:), pointer :: p_Dfunc1 => NULL()
    real(DP), dimension(:,:,:), pointer :: p_Dpoints => NULL()
    integer, dimension(:), pointer :: p_Ielements => NULL()

    type(t_vectorScalar), pointer :: secondFunction => NULL()
    real(DP), dimension(ExtFE_NDIM1,1) :: evaluationPoint
    real(DP), dimension(1) :: dFuncValue
    integer :: ielEval
    integer, dimension(1) :: iElement
    real(DP) , dimension (1) :: dpoint

    ! Get cubature weights
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Element numbers of the elements in process
    p_Ielements => rassemblyData%p_IelementList

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal

    ! Get a pointer to the values
    p_Dfunc1 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC1D)

    ! Get the second function
    secondFunction => rcollection%p_rvectorQuickAccess1%RvectorBlock(1)

    ! We should have everything prepared, so let us start
    ! Initial value of the integral
    Dintvalue(1) = 0.0_DP

    ! Loop over all elements
    do iel = 1,nelements

        ! Loop over all cubature points
        do icubp = 1,npointsPerElement

          ! Get the coordinates of the cubature point
          dx = p_Dpoints(1,icubp,iel)

          ! Now we check if we need to work or not
          ! If we need to work, the domainDescription
          ! will return a 1, if not it will return a 0

           call fparser_evalFunction(rcollection%p_rfparserQuickAccess1, &
                    rcollection%IquickAccess(1),(/dx/),work)

           ! Check if we need to work.
           ! If "work" is 1, we need to work. Since "work"
           ! can only be 0 or 1, this should run. However, you could
           ! hack this by defining functions that do not return 0 or 1
           ! as a function value, but. i.e. 0.5.
           if (abs(work - 1.0_DP) .le. 0.01_DP) then
                ! We need to work
                ! This one is getting really ugly
                ! Make the evaluation point as an array
                evaluationPoint(1,1) = dx
                ! Find the element nr.
                ! Usually the evaluation routines do this on their
                ! own - but they support only 2D-search algorithms
                ! So we search the element nr. on our own - because
                ! if we support the element nr. by ourselves, the
                ! evaluation routine does not branch into the searching part
                ! Therefore we need some preparation
                dpoint(1) = dx
                call tsrch_getElem_BruteForce(dpoint,&
                    secondFunction%p_rspatialDiscr%p_rtriangulation,ielEval)
                iElement(1) = ielEval

                ! Now we have the element number so we can call the
                ! evaluation routine
                call fevl_evaluate(DER_FUNC1D,dFuncValue,secondFunction, &
                    evaluationPoint,Ielements=iElement)
                ! Calculate the integral value
                dval1 = abs(p_Dfunc1(icubp,iel) - dFuncValue(1))
                Dintvalue(1) = Dintvalue(1) + dval1*p_DcubWeight(icubp,iel)
          end if


        end do ! icubp

    end do ! iel


end subroutine

! Calculates the L2-Norm of a 1D function
! on a Domain specified in the *.dat-files
subroutine calc_L2norm_onefunc_1D(Dintvalue,rassemblyData,rintAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)
    !<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
  !</input>

!<output>
    ! Returns the value of the integral(s)
    real(DP), dimension(:), intent(out) :: Dintvalue
!</output>

    ! Local variables
    real(DP) :: dx,dval1, work


    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight => NULL()
    real(DP), dimension(:,:), pointer :: p_Dfunc1 => NULL()
    real(DP), dimension(:,:,:), pointer :: p_Dpoints => NULL()
    integer, dimension(:), pointer :: p_Ielements => NULL()

    ! Get cubature weights
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Element numbers of the elements in process
    p_Ielements => rassemblyData%p_IelementList

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal

    ! Get a pointer to the values
    p_Dfunc1 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC)

    ! We should have everything prepared, so let us start
    ! Initial value of the integral
    Dintvalue(1) = 0.0_DP

    ! Loop over all elements
    do iel = 1,nelements

        ! Loop over all cubature points
        do icubp = 1,npointsPerElement

          ! Get the coordinates of the cubature point
          dx = p_Dpoints(1,icubp,iel)

          ! Now we check if we need to work or not
          ! If we need to work, the domainDescription
          ! will return a 1, if not it will return a 0
          ! We saved this one in the collection structure already
          ! and parsed it

           call fparser_evalFunction(rcollection%p_rfparserQuickAccess1,&
                        rcollection%IquickAccess(1),(/dx/),work)

           ! Check if we need to work.
           ! If "work" is 1, we need to work. Since "work"
           ! can only be 0 or 1, this should run. However, you could
           ! hack this by defining functions that do not return 0 or 1
           ! as a function value, but. i.e. 0.5.
           if (abs(work - 1.0_DP) .le. 0.01_DP) then
                ! We need to work
                dval1 = p_Dfunc1(icubp,iel)**2
                Dintvalue(1) = Dintvalue(1) + dval1*p_DcubWeight(icubp,iel)
          end if


        end do ! icubp

    end do ! iel

end subroutine

!*********************************************************************************************
! Calculates the L2-distance of 2 1D functions
! on a Domain specified in the *.dat-files

subroutine calc_L2error_twofunc_1D(Dintvalue,rassemblyData,rintAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)
    !<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
  !</input>

!<output>
    ! Returns the value of the integral(s)
    real(DP), dimension(:), intent(out) :: Dintvalue
!</output>

    ! Local variables
    real(DP) :: dx,dval1, work


    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight => NULL()
    real(DP), dimension(:,:), pointer :: p_Dfunc1 => NULL()
    real(DP), dimension(:,:,:), pointer :: p_Dpoints => NULL()
    integer, dimension(:), pointer :: p_Ielements => NULL()

    type(t_vectorScalar), pointer :: secondFunction => NULL()
    real(DP), dimension(ExtFE_NDIM1,1) :: evaluationPoint
    real(DP), dimension(1) :: dFuncValue
    integer :: ielEval
    integer, dimension(1) :: iElement
    real(DP) , dimension (1) :: dpoint

    ! Get cubature weights
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Element numbers of the elements in process
    p_Ielements => rassemblyData%p_IelementList

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal

    ! Get a pointer to the values
    p_Dfunc1 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC1D)

    ! Get the second function
    secondFunction => rcollection%p_rvectorQuickAccess1%RvectorBlock(1)

    ! We should have everything prepared, so let us start
    ! Initial value of the integral
    Dintvalue(1) = 0.0_DP

    ! Loop over all elements
    do iel = 1,nelements

        ! Loop over all cubature points
        do icubp = 1,npointsPerElement

          ! Get the coordinates of the cubature point
          dx = p_Dpoints(1,icubp,iel)

          ! Now we check if we need to work or not
          ! If we need to work, the domainDescription
          ! will return a 1, if not it will return a 0

           call fparser_evalFunction(rcollection%p_rfparserQuickAccess1, &
                    rcollection%IquickAccess(1),(/dx/),work)

           ! Check if we need to work.
           ! If "work" is 1, we need to work. Since "work"
           ! can only be 0 or 1, this should run. However, you could
           ! hack this by defining functions that do not return 0 or 1
           ! as a function value, but. i.e. 0.5.
           if (abs(work - 1.0_DP) .le. 0.01_DP) then
                ! We need to work
                ! This one is getting really ugly
                ! Make the evaluation point as an array
                evaluationPoint(1,1) = dx
                ! Find the element nr.
                ! Usually the evaluation routines do this on their
                ! own - but they support only 2D-search algorithms
                ! So we search the element nr. on our own - because
                ! if we support the element nr. by ourselves, the
                ! evaluation routine does not branch into the searching part
                ! Therefore we need some preparation
                dpoint(1) = dx
                call tsrch_getElem_BruteForce(dpoint,&
                    secondFunction%p_rspatialDiscr%p_rtriangulation,ielEval)
                iElement(1) = ielEval

                ! Now we have the element number so we can call the
                ! evaluation routine
                call fevl_evaluate(DER_FUNC1D,dFuncValue,secondFunction, &
                    evaluationPoint,Ielements=iElement)
                ! Calculate the integral value
                dval1 = (p_Dfunc1(icubp,iel) - dFuncValue(1))**2
                Dintvalue(1) = Dintvalue(1) + dval1*p_DcubWeight(icubp,iel)
          end if


        end do ! icubp

    end do ! iel


end subroutine


!---------------------------------------------------------------------------!
!                                                                           !
!                   2D integral routines come here                          !
!                                                                           !
!---------------------------------------------------------------------------!


!******************************************************************************


! Calculates the integral of a 2D function
! on a Domain specified in the *.dat-files

subroutine calc_integral_onefunc_2D(Dintvalue,rassemblyData,rintAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)
    !<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
  !</input>

!<output>
    ! Returns the value of the integral(s)
    real(DP), dimension(:), intent(out) :: Dintvalue
!</output>

    ! Local variables
    real(DP) :: dx,dy,dval1, work


    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight => NULL()
    real(DP), dimension(:,:), pointer :: p_Dfunc1 => NULL()
    real(DP), dimension(:,:,:), pointer :: p_Dpoints => NULL()
    integer, dimension(:), pointer :: p_Ielements => NULL()

    ! Get cubature weights
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Element numbers of the elements in process
    p_Ielements => rassemblyData%p_IelementList

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal

    ! Get a pointer to the values
    p_Dfunc1 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC)

    ! We should have everything prepared, so let us start
    ! Initial value of the integral
    Dintvalue(1) = 0.0_DP

    ! Loop over all elements
    do iel = 1,nelements

        ! Loop over all cubature points
        do icubp = 1,npointsPerElement

          ! Get the coordinates of the cubature point
          dx = p_Dpoints(1,icubp,iel)
          dy = p_Dpoints(2,icubp,iel)

          ! Now we check if we need to work or not
          ! If we need to work, the domainDescription
          ! will return a 1, if not it will return a 0
          ! We already parsed everything

           call fparser_evalFunction(rcollection%p_rfparserQuickAccess1, &
                    rcollection%IquickAccess(1),(/dx,dy/),work)

           ! Check if we need to work.
           ! If "work" is 1, we need to work. Since "work"
           ! can only be 0 or 1, this should run. However, you could
           ! hack this by defining functions that do not return 0 or 1
           ! as a function value, but. i.e. 0.5.
           if (abs(work - 1.0_DP) .le. 0.01_DP) then
                ! We need to work
                dval1 = (p_Dfunc1(icubp,iel))
                Dintvalue(1) = Dintvalue(1) + dval1*p_DcubWeight(icubp,iel)
          end if


        end do ! icubp

    end do ! iel


end subroutine


! Calculates the integral (f1 - f2)dx with 2 2D functions
! on a Domain specified in the *.dat-files

subroutine calc_integralDifference_twofunc_2D(Dintvalue,rassemblyData,rintAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)
    !<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
  !</input>

!<output>
    ! Returns the value of the integral(s)
    real(DP), dimension(:), intent(out) :: Dintvalue
!</output>

    ! Local variables
    real(DP) :: dx,dy,dval1, work


    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight => NULL()
    real(DP), dimension(:,:), pointer :: p_Dfunc1 => NULL()
    real(DP), dimension(:,:,:), pointer :: p_Dpoints => NULL()
    integer, dimension(:), pointer :: p_Ielements => NULL()

    type(t_vectorScalar), pointer :: secondFunction => NULL()
    real(DP), dimension(1) :: dfuncValue
    real(DP), dimension(ExtFE_NDIM2,1) :: evaluationPoint


    ! Get cubature weights
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Element numbers of the elements in process
    p_Ielements => rassemblyData%p_IelementList

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal

    ! Get a pointer to the values
    p_Dfunc1 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC2D)
    ! Get the second function
    secondFunction => rcollection%p_rvectorQuickAccess1%RvectorBlock(1)

    ! We should have everything prepared, so let us start
    ! Initial value of the integral
    Dintvalue(1) = 0.0_DP

    ! Loop over all elements
    do iel = 1,nelements

        ! Loop over all cubature points
        do icubp = 1,npointsPerElement

          ! Get the coordinates of the cubature point
          dx = p_Dpoints(1,icubp,iel)
          dy = p_Dpoints(2,icubp,iel)

          ! Now we check if we need to work or not
          ! If we need to work, the domainDescription
          ! will return a 1, if not it will return a 0
          ! We already parsed everything

           call fparser_evalFunction(rcollection%p_rfparserQuickAccess1,&
                    rcollection%IquickAccess(1),(/dx,dy/),work)

           ! Check if we need to work.
           ! If "work" is 1, we need to work. Since "work"
           ! can only be 0 or 1, this should run. However, you could
           ! hack this by defining functions that do not return 0 or 1
           ! as a function value, but. i.e. 0.5.
           if (abs(work - 1.0_DP) .le. 0.01_DP) then
                ! We need to work
                ! Evaluate the second function
                evaluationPoint(1,1) = dx
                evaluationPoint(2,1) = dy
                call fevl_evaluate(DER_FUNC2D,dFuncValue,secondFunction, &
                    evaluationPoint)
                dval1 = (p_Dfunc1(icubp,iel) -  dfuncValue(1))
                Dintvalue(1) = Dintvalue(1) + dval1*p_DcubWeight(icubp,iel)
          end if


        end do ! icubp

    end do ! iel

end subroutine


! Calculates the L1-Norm of a 2D function
! on a Domain specified in the *.dat-files

subroutine calc_L1norm_onefunc_2D(Dintvalue,rassemblyData,rintAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)
    !<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
  !</input>

!<output>
    ! Returns the value of the integral(s)
    real(DP), dimension(:), intent(out) :: Dintvalue
!</output>

    ! Local variables
    real(DP) :: dx,dy,dval1, work


    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight => NULL()
    real(DP), dimension(:,:), pointer :: p_Dfunc1 => NULL()
    real(DP), dimension(:,:,:), pointer :: p_Dpoints => NULL()
    integer, dimension(:), pointer :: p_Ielements => NULL()

    ! Get cubature weights
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Element numbers of the elements in process
    p_Ielements => rassemblyData%p_IelementList

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal

    ! Get a pointer to the values
    p_Dfunc1 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC)

    ! We should have everything prepared, so let us start
    ! Initial value of the integral
    Dintvalue(1) = 0.0_DP

    ! Loop over all elements
    do iel = 1,nelements

        ! Loop over all cubature points
        do icubp = 1,npointsPerElement

          ! Get the coordinates of the cubature point
          dx = p_Dpoints(1,icubp,iel)
          dy = p_Dpoints(2,icubp,iel)

          ! Now we check if we need to work or not
          ! If we need to work, the domainDescription
          ! will return a 1, if not it will return a 0
          ! We already parsed everything

           call fparser_evalFunction(rcollection%p_rfparserQuickAccess1, &
                    rcollection%IquickAccess(1),(/dx,dy/),work)

           ! Check if we need to work.
           ! If "work" is 1, we need to work. Since "work"
           ! can only be 0 or 1, this should run. However, you could
           ! hack this by defining functions that do not return 0 or 1
           ! as a function value, but. i.e. 0.5.
           if (abs(work - 1.0_DP) .le. 0.01_DP) then
                ! We need to work
                dval1 = abs(p_Dfunc1(icubp,iel))
                Dintvalue(1) = Dintvalue(1) + dval1*p_DcubWeight(icubp,iel)
          end if


        end do ! icubp

    end do ! iel


end subroutine


! Calculates the L1-Distance between 2 2D functions
! on a Domain specified in the *.dat-files

subroutine calc_L1error_twofunc_2D(Dintvalue,rassemblyData,rintAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)
    !<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
  !</input>

!<output>
    ! Returns the value of the integral(s)
    real(DP), dimension(:), intent(out) :: Dintvalue
!</output>

    ! Local variables
    real(DP) :: dx,dy,dval1, work


    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight => NULL()
    real(DP), dimension(:,:), pointer :: p_Dfunc1 => NULL()
    real(DP), dimension(:,:,:), pointer :: p_Dpoints => NULL()
    integer, dimension(:), pointer :: p_Ielements => NULL()

    type(t_vectorScalar), pointer :: secondFunction => NULL()
    real(DP), dimension(1) :: dfuncValue
    real(DP), dimension(ExtFE_NDIM2,1) :: evaluationPoint


    ! Get cubature weights
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Element numbers of the elements in process
    p_Ielements => rassemblyData%p_IelementList

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal

    ! Get a pointer to the values
    p_Dfunc1 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC2D)
    ! Get the second function
    secondFunction => rcollection%p_rvectorQuickAccess1%RvectorBlock(1)

    ! We should have everything prepared, so let us start
    ! Initial value of the integral
    Dintvalue(1) = 0.0_DP

    ! Loop over all elements
    do iel = 1,nelements

        ! Loop over all cubature points
        do icubp = 1,npointsPerElement

          ! Get the coordinates of the cubature point
          dx = p_Dpoints(1,icubp,iel)
          dy = p_Dpoints(2,icubp,iel)

          ! Now we check if we need to work or not
          ! If we need to work, the domainDescription
          ! will return a 1, if not it will return a 0
          ! We already parsed everything

           call fparser_evalFunction(rcollection%p_rfparserQuickAccess1,&
                    rcollection%IquickAccess(1),(/dx,dy/),work)

           ! Check if we need to work.
           ! If "work" is 1, we need to work. Since "work"
           ! can only be 0 or 1, this should run. However, you could
           ! hack this by defining functions that do not return 0 or 1
           ! as a function value, but. i.e. 0.5.
           if (abs(work - 1.0_DP) .le. 0.01_DP) then
                ! We need to work
                ! Evaluate the second function
                evaluationPoint(1,1) = dx
                evaluationPoint(2,1) = dy
                call fevl_evaluate(DER_FUNC2D,dFuncValue,secondFunction, &
                    evaluationPoint)
                dval1 = abs(p_Dfunc1(icubp,iel) -  dfuncValue(1))
                Dintvalue(1) = Dintvalue(1) + dval1*p_DcubWeight(icubp,iel)
          end if


        end do ! icubp

    end do ! iel

end subroutine



!******************************************************************************

! Calculates the L2-Norm of a 2D function
! on a Domain specified in the *.dat-files

subroutine calc_L2norm_onefunc_2D(Dintvalue,rassemblyData,rintAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)
    !<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
  !</input>

!<output>
    ! Returns the value of the integral(s)
    real(DP), dimension(:), intent(out) :: Dintvalue
!</output>

    ! Local variables
    real(DP) :: dx,dy,dval1, work


    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight => NULL()
    real(DP), dimension(:,:), pointer :: p_Dfunc1 => NULL()
    real(DP), dimension(:,:,:), pointer :: p_Dpoints => NULL()
    integer, dimension(:), pointer :: p_Ielements => NULL()

    ! Get cubature weights
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Element numbers of the elements in process
    p_Ielements => rassemblyData%p_IelementList

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal

    ! Get a pointer to the values
    p_Dfunc1 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC)

    ! We should have everything prepared, so let us start
    ! Initial value of the integral
    Dintvalue(1) = 0.0_DP

    ! Loop over all elements
    do iel = 1,nelements

        ! Loop over all cubature points
        do icubp = 1,npointsPerElement

          ! Get the coordinates of the cubature point
          dx = p_Dpoints(1,icubp,iel)
          dy = p_Dpoints(2,icubp,iel)

          ! Now we check if we need to work or not
          ! If we need to work, the domainDescription
          ! will return a 1, if not it will return a 0
          ! We already parsed everything

           call fparser_evalFunction(rcollection%p_rfparserQuickAccess1, &
                    rcollection%IquickAccess(1),(/dx,dy/),work)

           ! Check if we need to work.
           ! If "work" is 1, we need to work. Since "work"
           ! can only be 0 or 1, this should run. However, you could
           ! hack this by defining functions that do not return 0 or 1
           ! as a function value, but. i.e. 0.5.
           if (abs(work - 1.0_DP) .le. 0.01_DP) then
                ! We need to work
                dval1 = p_Dfunc1(icubp,iel)**2
                Dintvalue(1) = Dintvalue(1) + dval1*p_DcubWeight(icubp,iel)
          end if


        end do ! icubp

    end do ! iel


end subroutine


! Calculates the L2-Distance between 2 2D functions
! on a Domain specified in the *.dat-files

subroutine calc_L2error_twofunc_2D(Dintvalue,rassemblyData,rintAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)
    !<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
  !</input>

!<output>
    ! Returns the value of the integral(s)
    real(DP), dimension(:), intent(out) :: Dintvalue
!</output>

    ! Local variables
    real(DP) :: dx,dy,dval1, work


    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight => NULL()
    real(DP), dimension(:,:), pointer :: p_Dfunc1 => NULL()
    real(DP), dimension(:,:,:), pointer :: p_Dpoints => NULL()
    integer, dimension(:), pointer :: p_Ielements => NULL()

    type(t_vectorScalar), pointer :: secondFunction => NULL()
    real(DP), dimension(1) :: dfuncValue
    real(DP), dimension(ExtFE_NDIM2,1) :: evaluationPoint


    ! Get cubature weights
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Element numbers of the elements in process
    p_Ielements => rassemblyData%p_IelementList

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal

    ! Get a pointer to the values
    p_Dfunc1 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC2D)
    ! Get the second function
    secondFunction => rcollection%p_rvectorQuickAccess1%RvectorBlock(1)

    ! We should have everything prepared, so let us start
    ! Initial value of the integral
    Dintvalue(1) = 0.0_DP

    ! Loop over all elements
    do iel = 1,nelements

        ! Loop over all cubature points
        do icubp = 1,npointsPerElement

          ! Get the coordinates of the cubature point
          dx = p_Dpoints(1,icubp,iel)
          dy = p_Dpoints(2,icubp,iel)

          ! Now we check if we need to work or not
          ! If we need to work, the domainDescription
          ! will return a 1, if not it will return a 0
          ! We already parsed everything

           call fparser_evalFunction(rcollection%p_rfparserQuickAccess1,&
                    rcollection%IquickAccess(1),(/dx,dy/),work)

           ! Check if we need to work.
           ! If "work" is 1, we need to work. Since "work"
           ! can only be 0 or 1, this should run. However, you could
           ! hack this by defining functions that do not return 0 or 1
           ! as a function value, but. i.e. 0.5.
           if (abs(work - 1.0_DP) .le. 0.01_DP) then
                ! We need to work
                ! Evaluate the second function
                evaluationPoint(1,1) = dx
                evaluationPoint(2,1) = dy
                call fevl_evaluate(DER_FUNC2D,dFuncValue,secondFunction, &
                    evaluationPoint)
                dval1 = (p_Dfunc1(icubp,iel) -  dfuncValue(1))**2
                Dintvalue(1) = Dintvalue(1) + dval1*p_DcubWeight(icubp,iel)
          end if


        end do ! icubp

    end do ! iel

end subroutine

!---------------------------------------------------------------------------!
!                                                                           !
!                   3D integral routines come here                          !
!                                                                           !
!---------------------------------------------------------------------------!

! Calculates the L2-Norm of a 3D function
! on a Domain specified in the *.dat-files

subroutine calc_L2norm_onefunc_3D(Dintvalue,rassemblyData,rintAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)
    !<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
  !</input>

!<output>
    ! Returns the value of the integral(s)
    real(DP), dimension(:), intent(out) :: Dintvalue
!</output>

    ! Local variables
    real(DP) :: dx,dy,dz,dval1, work


    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight => NULL()
    real(DP), dimension(:,:), pointer :: p_Dfunc1 => NULL()
    real(DP), dimension(:,:,:), pointer :: p_Dpoints => NULL()
    integer, dimension(:), pointer :: p_Ielements => NULL()

    ! Get cubature weights
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Element numbers of the elements in process
    p_Ielements => rassemblyData%p_IelementList

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal

    ! Get a pointer to the values
    p_Dfunc1 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC)

    ! We should have everything prepared, so let us start
    ! Initial value of the integral
    Dintvalue(1) = 0.0_DP

    ! Loop over all elements
    do iel = 1,nelements

        ! Loop over all cubature points
        do icubp = 1,npointsPerElement

          ! Get the coordinates of the cubature point
          dx = p_Dpoints(1,icubp,iel)
          dy = p_Dpoints(2,icubp,iel)
          dz = p_Dpoints(3,icubp,iel)

          ! Now we check if we need to work or not
          ! If we need to work, the domainDescription
          ! will return a 1, if not it will return a 0
          ! We already parsed everything

           call fparser_evalFunction(rcollection%p_rfparserQuickAccess1,&
                        rcollection%IquickAccess(1),(/dx,dy,dz/),work)

           ! Check if we need to work.
           ! If "work" is 1, we need to work. Since "work"
           ! can only be 0 or 1, this should run. However, you could
           ! hack this by defining functions that do not return 0 or 1
           ! as a function value, but. i.e. 0.5.
           if (abs(work - 1.0_DP) .le. 0.01_DP) then
                ! We need to work
                dval1 = p_Dfunc1(icubp,iel)**2
                Dintvalue(1) = Dintvalue(1) + dval1*p_DcubWeight(icubp,iel)
          end if


        end do ! icubp

    end do ! iel


end subroutine


!-----------------------------------------------------------------------------------------------------------------------------------!

! This function is just a wrapper for the functions provided
! by the kernel. We need to branch for dimension etc,so it is
! easier to do it like this

subroutine ExtFE_eval_function(fefunction,Dpoint,CDeriv_ExtFE,Dvalue)
    !<intput>
    ! The function that we want to evaluate
    type(t_vectorScalar), intent(inout) :: fefunction
    ! The point where we want to evaluate
    real(DP), dimension(:), intent(inout)  :: Dpoint
    ! Which derivative
    integer, intent(in)  :: CDeriv_ExtFE
    !</input>

    !<output>
    real(DP), intent(out) :: Dvalue
    !</output>


    ! local variables
    ! We have to convert the coordinates
    real(DP), dimension(:,:),allocatable :: evalPoint
    integer, dimension(1) :: evalElement
    integer :: nDim, iElement
    integer :: CDeriv
    real(DP), dimension(1) :: locValue
    nDim = fefunction%p_rspatialDiscr%ndimension

    ! We have to do several branches here:
    ! We need to find out the element number
    ! and set if we want to evaluate the function
    ! or derivative
    iElement = 0

    select case (nDim)

      ! 1D
      case(ExtFE_NDIM1)
        ! We have to find out the element Nr.
        call tsrch_getElem_BruteForce(Dpoint, &
           fefunction%p_rspatialDiscr%p_rtriangulation,iElement)
        ! Now translate the ExtFE_derivative names to global Feat2
        select case(CDeriv_ExtFE)
            case(ExtFE_DER_FUNC_1D)
                CDeriv = DER_FUNC1D
            case(ExtFE_DER_DERIV_1D_X)
                CDeriv = DER_DERIV1D_X
            case default
                 call output_line('Derivative not supported', &
                    OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_eval_function")
                    call sys_halt()
        end select
        ! Save the evaluation point
        allocate(evalPoint(ExtFE_NDIM1,1))
        evalPoint(1,1) = Dpoint(1)

      case(ExtFE_NDIM2)
        ! Find the element
        call tsrch_getElem_raytrace2D(Dpoint, &
            fefunction%p_rspatialDiscr%p_rtriangulation,iElement)
        ! if we did not find it, go to bruteforce
        if (iElement .eq. 0) then
            call tsrch_getElem_BruteForce(Dpoint, &
                fefunction%p_rspatialDiscr%p_rtriangulation,iElement)
        end if

        ! Translate the local ExtFE-Names to global Feat2
        select case(CDeriv_ExtFE)
            case(ExtFE_DER_FUNC_2D)
                CDeriv = DER_FUNC2D
            case(ExtFE_DER_DERIV_2D_X)
                CDeriv = DER_DERIV2D_X
            case(ExtFE_DER_DERIV_2D_Y)
                CDeriv = DER_DERIV2D_Y
            case default
                 call output_line('Derivative not supported', &
                    OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_eval_function")
                    call sys_halt()
        end select

        allocate(evalPoint(ExtFE_NDIM2,1))
        evalPoint(1,1) = Dpoint(1)
        evalPoint(2,1) = Dpoint(2)

        case(ExtFE_NDIM3)
          call tsrch_getElem_raytrace3D(Dpoint, &
                fefunction%p_rspatialDiscr%p_rtriangulation,iElement)
          if (iElement .eq. 0 ) then
              call tsrch_getElem_BruteForce(Dpoint, &
                  fefunction%p_rspatialDiscr%p_rtriangulation,iElement)
          end if

        select case(CDeriv_ExtFE)
            case(ExtFE_DER_FUNC_3D)
                CDeriv = DER_FUNC3D
            case(ExtFE_DER_DERIV_3D_X)
                CDeriv = DER_DERIV3D_X
            case(ExtFE_DER_DERIV_3D_Y)
                CDeriv = DER_DERIV3D_Y
            case(ExtFE_DER_DERIV_3D_Z)
                CDeriv = DER_DERIV3D_Z
            case default
                 call output_line('Derivative not supported', &
                    OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_eval_function")
                    call sys_halt()
        end select

        allocate(evalPoint(ExtFE_NDIM3,1))
        evalPoint(1,1) = Dpoint(1)
        evalPoint(2,1) = Dpoint(2)
        evalPoint(3,1) = Dpoint(3)

    end select

    ! Okay, we have all we need: The global constants
    ! and the element number. So now we can call the real evaluation
    ! routine provided by the kernel
    evalElement(1) = iElement
    call fevl_evaluate(CDeriv,locValue,fefunction,evalPoint,Ielements=evalElement)

    ! Set the output
    Dvalue = locValue(1)

    ! Clean up
    deallocate(evalPoint)

end subroutine


end module
