! This module contains the core of the compare-application
! It starts at the point where you have 2 FEM Solutions,
! created the according discretisations and parametrisations
! and now want to calculate the differences.

module ExtFEcomparer_core

  use fsystem
  use storage
  use genoutput

  use element
  use cubature
  use spatialdiscretisation
  use linearsystemblock
  use linearsystemscalar

  use bilinearformevaluation
  use derivatives
  use collection
  use feevaluation2
  use feevaluation

  use blockmatassembly
  use blockmatassemblybase
  use blockmatassemblystdop

  use ExtFEcomparer_typedefs

  use fparser
  use paramlist

  use vectorio
  use io

implicit none

contains

  ! This routine shall prepare everything for the calculation
  ! It is the one that is called from outside. The other routines
  ! (that do the actual calculation) are called by this routine.

  subroutine calculate_difference_2D(rproblem_1,rproblem_2)

  type(t_problem), intent(inout), target :: rproblem_1, rproblem_2

  ! Level where we want to evaluate
  integer :: lvl1,lvl2

  ! We also need to create 2 block vectors so we can attatch
  ! the discretisation to it
  type(t_vectorBlock), pointer :: rx1, rx2

  ! We want to calculate integrals, so we need the cubature info
  type(t_scalarCubatureInfo), target:: rcubatureInfoVelocityFirst, rcubatureInfoVelocitySecond
  type(t_scalarCubatureInfo), target:: rcubatureInfoPressureFirst, rcubatureInfoPressureSecond


  ! We need to evaluate our FEM-function in the cubature points
  ! for that purpose we need a vector
  type(t_fev2Vectors) :: rcoeffVectors

 ! We need the collection to pass the second vector to the
 ! actual computation
 type(t_collection) :: rcollection

  ! Local variables
  real(DP) :: dintvalue
  integer :: calcL2DifferenceVelocity, calcL2NormVelocityFirst, calcL2NormVelocitySecond
  integer :: calcL2DifferencePressure, calcL2NormPressureFirst, calcL2NormPressureSecond
  integer :: calcPointDifference, calcPointFirst, calcPointSecond

  ! For the output
  real(DP), dimension(6) :: L2NormResults
  integer :: writeOutL2Norm, writeOutPointValues
  character(LEN=SYS_STRLEN):: outFileNameL2Norm, outFileNamePointValues

  ! for the point values we need an array of points
  real(DP), dimension(:,:), allocatable :: evaluationCoords
  integer, dimension(:), allocatable :: Cder
  integer, dimension(:), allocatable :: IType
  integer :: nPoints,i
  real(DP), dimension(:), allocatable :: pointvaluesFirst, pointvaluesSecond, pointValueDifference
  character(LEN=SYS_STRLEN):: sparam

  ! To write out the point values
  integer :: iunit
  logical :: bfileexists
  character(LEN=SYS_STRLEN):: outputString
  character(LEN=10), dimension(3,3), parameter :: Sfctnames = reshape (&
      (/ "       u1 ","     u1_x ","     u1_y " , &
         "       u2 ","     u2_x ","     u2_y " , &
         "        p ","      p_x ","      p_y " /) ,&
       (/ 3,3 /) )

  !=============================================================!
  ! First of all, we read in what the user wants to calculate.  !
  !=============================================================!

  ! Which L2-Norms?
  call parlst_getvalue_int(rproblem_1%rparamlist,"GENERALSETTINGS","calcL2DifferenceVelocity",calcL2DifferenceVelocity,0)
  call parlst_getvalue_int(rproblem_1%rparamlist,"GENERALSETTINGS","calcL2NormVelocityFirst",calcL2NormVelocityFirst,0)
  call parlst_getvalue_int(rproblem_1%rparamlist,"GENERALSETTINGS","calcL2NormVelocitySecond",calcL2NormVelocitySecond,0)
  call parlst_getvalue_int(rproblem_1%rparamlist,"GENERALSETTINGS","calcL2DifferencePressure",calcL2DifferencePressure,0)
  call parlst_getvalue_int(rproblem_1%rparamlist,"GENERALSETTINGS","calcL2NormPressureFirst",calcL2NormPressureFirst,0)
  call parlst_getvalue_int(rproblem_1%rparamlist,"GENERALSETTINGS","calcL2NormPressureSecond",calcL2NormPressureSecond,0)

  ! Which type of point values?
  call parlst_getvalue_int(rproblem_1%rparamlist,"GENERALSETTINGS","calcPointDifference",calcPointDifference,0)
  call parlst_getvalue_int(rproblem_1%rparamlist,"GENERALSETTINGS","calcPointFirst",calcPointFirst,0)
  call parlst_getvalue_int(rproblem_1%rparamlist,"GENERALSETTINGS","calcPointSecond",calcPointSecond,0)

  ! Postprocessing settings
  ! L2-Norms
  call parlst_getvalue_int(rproblem_1%rparamlist,"POSTPROCESSING","writeOutL2Norm",writeOutL2Norm,0)
  call parlst_getvalue_string(rproblem_1%rparamList,"POSTPROCESSING", "sfilenameL2Norm",outFileNameL2Norm,bdequote=.true.)

  ! Pointvalues
  call parlst_getvalue_int(rproblem_1%rparamlist,"POSTPROCESSING","writeOutPointValues",writeOutPointValues,0)
  call parlst_getvalue_string(rproblem_1%rparamList,"POSTPROCESSING", "sfilenamePointValues",outFileNamePointValues,bdequote=.true.)

  !====================!
  ! Final preparations !
  !====================!

  lvl1 = rproblem_1%NLMAX
  lvl2 = rproblem_2%NLMAX


  ! This part is actually neccesary because we want to have the
  ! situation that each block vector has its blockdiscretisation
  ! attatched to it.

  rx1 => rproblem_1%coeffVector
  rx2 => rproblem_2%coeffVector

  ! Attatch the Block-Discritisation to rx1 and rx2
  rx1%p_rblockDiscr => rproblem_1%RlevelInfo(lvl1)%rdiscretisation
  rx2%p_rblockDiscr => rproblem_2%RlevelInfo(lvl2)%rdiscretisation

  ! Now we attatch the spatial discretisations to rx1 and rx2
  rx1%RvectorBlock(1)%p_rspatialDiscr => rproblem_1%RlevelInfo(lvl1)%rdiscretisation%RspatialDiscr(1)
  rx1%RvectorBlock(2)%p_rspatialDiscr => rproblem_1%RlevelInfo(lvl1)%rdiscretisation%RspatialDiscr(2)
  rx1%RvectorBlock(3)%p_rspatialDiscr => rproblem_1%RlevelInfo(lvl1)%rdiscretisation%RspatialDiscr(3)

  rx2%RvectorBlock(1)%p_rspatialDiscr => rproblem_2%RlevelInfo(lvl2)%rdiscretisation%RspatialDiscr(1)
  rx2%RvectorBlock(2)%p_rspatialDiscr => rproblem_2%RlevelInfo(lvl2)%rdiscretisation%RspatialDiscr(2)
  rx2%RvectorBlock(3)%p_rspatialDiscr => rproblem_2%RlevelInfo(lvl2)%rdiscretisation%RspatialDiscr(3)

  ! In the subroutines we need the parameter list
  ! The lists in both problem structures are the same, so
  ! we can pick any
  rcollection%p_rparlistQuickAccess1 => rproblem_1%rparamlist


  !======================================!
  ! We might have everything, so we start!
  !======================================!

  !--------------------------------------!
  ! We start with the velocity integrals !
  !--------------------------------------!

  ! ##################################
  ! # FIRST Integral: ||u1 - u2||_L2 #
  ! ##################################
  ! Init the result with 0
  L2NormResults(1) = 0.0_DP

  ! Check if the user wants to calculate this:
    if (abs(calcL2DifferenceVelocity - 1.0_DP) .le. 0.1_DP) then

        ! We need the cubature informations
        ! First step: only for Integral ||u1-u2|| on common part, we can
        ! actually use the cubature info from any discretisation
        ! So we use the discrectisation from vector 1

        call spdiscr_createDefCubStructure(rx1%RvectorBlock(1)%p_rspatialDiscr, &
               rcubatureInfoVelocityFirst,rproblem_1%I_Cubature_Formula)

        ! first the u component
        call fev2_addVectorToEvalList(rcoeffVectors,rx1%RvectorBlock(1),0)
        ! then the v component
        call fev2_addVectorToEvalList(rcoeffVectors,rx1%RvectorBlock(2),0)

        ! we cannot pass the second vector via this way because it lives
        ! on a different mesh. Therefore we pass it via the collection
        ! The evaluation of this one will be expensive
        rcollection%p_rvectorQuickAccess1=> rx2


        ! We need space in the vector for the second function
        call fev2_addDummyVectorToEvalList(rcoeffVectors)
        call fev2_addDummyVectorToEvalList(rcoeffVectors)

        ! Calculate the L2 difference
        call bma_buildIntegral(dintvalue,BMA_CALC_STANDARD,calc_L2error_velocity_2D_common, &
            rcollection=rcollection,revalVectors=rcoeffVectors,rcubatureInfo=rcubatureInfoVelocityFirst)

        dintvalue = sqrt(dintvalue)
        call output_line ("L2-difference of the velocity on the common part= "//trim(sys_sdEL(dintvalue,10)))

        ! Store it in the result array
        L2NormResults(1) = dintvalue

        ! Release everything
        call fev2_releaseVectorList(rcoeffVectors)
        call spdiscr_releaseCubStructure(rcubatureInfoVelocityFirst)
        rcollection%p_rvectorQuickAccess1 => NULL()
        dintvalue = 0.0_DP

    end if



    !##############################!
    !# Second integral: ||u1||_L2 #!
    !##############################!
    ! init the result with 0
    L2NormResults(2) = 0.0_DP

    ! We need the name of the variable in the dat-file in the
    ! callback routine
    rcollection%SquickAccess(1) = 'onlyFirst'

    ! Check if the user wants to calculate this

    if (abs(calcL2NormVelocityFirst - 1.0_DP) .le. 0.1_DP) then
        ! Pass function 1 to bma_integral
        call fev2_addVectorToEvalList(rcoeffVectors,rx1%RvectorBlock(1),0)
        call fev2_addVectorToEvalList(rcoeffVectors,rx1%RvectorBlock(2),0)

        ! Create the according cubature structure
        call spdiscr_createDefCubStructure(rx1%RvectorBlock(1)%p_rspatialDiscr, &
               rcubatureInfoVelocityFirst,rproblem_1%I_Cubature_Formula)

        ! Calculate the L2 error
        call bma_buildIntegral(dintvalue,BMA_CALC_STANDARD,calc_L2norm_velocity_2D, &
            rcollection=rcollection,revalVectors=rcoeffVectors,rcubatureInfo=rcubatureInfoVelocityFirst)

        dintvalue = sqrt(dintvalue)
        call output_line ("L2-norm of the velocity of the first function= "//trim(sys_sdEL(dintvalue,10)))

        ! Store it in the result array
        L2NormResults(2) = dintvalue

        ! Release everything
        call fev2_releaseVectorList(rcoeffVectors)
        call spdiscr_releaseCubStructure(rcubatureInfoVelocityFirst)
        dintvalue = 0.0_DP

    end if

    !#############################!
    !# Third integral: ||u2||_L2 #!
    !#############################!
    ! init the result with 0
    L2NormResults(3) = 0.0_DP

    ! We need the name of the variable in the dat-file in the
    ! callback routine
    rcollection%SquickAccess(1) = 'onlySecond'

    ! Check if the user wants to calculate this

    if (abs( calcL2NormVelocitySecond- 1.0_DP) .le. 0.1_DP) then
        ! Pass function 2 to bma_build_integral
        call fev2_addVectorToEvalList(rcoeffVectors,rx2%RvectorBlock(1),0)
        call fev2_addVectorToEvalList(rcoeffVectors,rx2%RvectorBlock(2),0)

        ! Create the according cubature structure
        call spdiscr_createDefCubStructure(rx2%RvectorBlock(1)%p_rspatialDiscr, &
                rcubatureInfoVelocitySecond,rproblem_2%I_Cubature_Formula)

        ! Calculate the L2 error
        call bma_buildIntegral(dintvalue,BMA_CALC_STANDARD,calc_L2norm_velocity_2D, &
            rcollection=rcollection,revalVectors=rcoeffVectors,rcubatureInfo=rcubatureInfoVelocitySecond)

        dintvalue = sqrt(dintvalue)
        call output_line ("L2-norm of the velocity of the second function= "//trim(sys_sdEL(dintvalue,10)))

        ! Write it out in the results array
        L2NormResults(3) = dintvalue

        ! Release everything
        call fev2_releaseVectorList(rcoeffVectors)
        call spdiscr_releaseCubStructure(rcubatureInfoVelocitySecond)
        dintvalue = 0.0_DP

    end if

    !----------------------------------!
    ! Now we do the pressure integrals !
    !----------------------------------!

    !###################################!
    !# Fourth integral: ||p1 - p2||_L2 #!
    !###################################!

      ! Check if the user wants to calculate this:
    if (abs(calcL2DifferencePressure - 1.0_DP) .le. 0.1_DP) then

        ! We need the cubature informations
        ! First step: only for Integral ||p1-p2|| on common part, we can
        ! actually use the cubature info from any discretisation
        ! So we use the discrectisation from vector 1

        call spdiscr_createDefCubStructure(rx1%RvectorBlock(3)%p_rspatialDiscr, &
               rcubatureInfoPressureFirst,rproblem_1%I_Cubature_Formula)

        ! Add it to the evaluation list
        call fev2_addVectorToEvalList(rcoeffVectors,rx1%RvectorBlock(3),0)

        ! we cannot pass the second vector via this way because it lives
        ! on a different mesh. Therefore we pass it via the collection
        ! The evaluation of this one will be expensive
        rcollection%p_rvectorQuickAccess1=> rx2


        ! We need space in the vector for the second function
        call fev2_addDummyVectorToEvalList(rcoeffVectors)

        ! Calculate the L2 error
        call bma_buildIntegral(dintvalue,BMA_CALC_STANDARD,calc_L2error_pressure_2D_common, &
            rcollection=rcollection,revalVectors=rcoeffVectors,rcubatureInfo=rcubatureInfoPressureFirst)

        dintvalue = sqrt(dintvalue)
        call output_line ("L2-difference of the pressure on the common part= "//trim(sys_sdEL(dintvalue,10)))

        ! Store it in the result array
        L2NormResults(4) = dintvalue

        ! Release everything
        call fev2_releaseVectorList(rcoeffVectors)
        call spdiscr_releaseCubStructure(rcubatureInfoPressureFirst)
        rcollection%p_rvectorQuickAccess1 => NULL()
        dintvalue = 0.0_DP

    end if

    !#############################!
    !# Fifth integral: ||p1||_L2 #!
    !#############################!

    ! init the result with 0
    L2NormResults(5) = 0.0_DP

    ! We need the name of the variable in the dat-file in the
    ! callback routine
    rcollection%SquickAccess(1) = 'onlyFirst'

    ! Check if the user wants to calculate this

    if (abs(calcL2NormPressureFirst - 1.0_DP) .le. 0.1_DP) then
        ! Pass function 1 to bma_integral
        call fev2_addVectorToEvalList(rcoeffVectors,rx1%RvectorBlock(3),0)

        ! Create the according cubature structure
        call spdiscr_createDefCubStructure(rx1%RvectorBlock(3)%p_rspatialDiscr, &
               rcubatureInfoPressureFirst,rproblem_1%I_Cubature_Formula)

        ! Calculate the L2 error
        call bma_buildIntegral(dintvalue,BMA_CALC_STANDARD,calc_L2norm_pressure_2D, &
            rcollection=rcollection,revalVectors=rcoeffVectors,rcubatureInfo=rcubatureInfoPressureFirst)

        dintvalue = sqrt(dintvalue)
        call output_line ("L2-norm of the pressure of the first function= "//trim(sys_sdEL(dintvalue,10)))

        ! Store it in the result array
        L2NormResults(5) = dintvalue

        ! Release everything
        call fev2_releaseVectorList(rcoeffVectors)
        call spdiscr_releaseCubStructure(rcubatureInfoPressureFirst)
        dintvalue = 0.0_DP

    end if


    !#############################!
    !# Sixth integral: ||p2||_L2 #!
    !#############################!
    ! init the result with 0
    L2NormResults(6) = 0.0_DP

    ! We need the name of the variable in the dat-file in the
    ! callback routine
    rcollection%SquickAccess(1) = 'onlySecond'

    ! Check if the user wants to calculate this

    if (abs( calcL2NormPressureSecond- 1.0_DP) .le. 0.1_DP) then
        ! Pass function 2 to bma_build_integral
        call fev2_addVectorToEvalList(rcoeffVectors,rx2%RvectorBlock(3),0)

        ! Create the according cubature structure
        call spdiscr_createDefCubStructure(rx2%RvectorBlock(3)%p_rspatialDiscr, &
                rcubatureInfoPressureSecond,rproblem_2%I_Cubature_Formula)

        ! Calculate the L2 error
        call bma_buildIntegral(dintvalue,BMA_CALC_STANDARD,calc_L2norm_pressure_2D, &
            rcollection=rcollection,revalVectors=rcoeffVectors,rcubatureInfo=rcubatureInfoPressureSecond)

        dintvalue = sqrt(dintvalue)
        call output_line ("L2-norm of the pressure of the second function= "//trim(sys_sdEL(dintvalue,10)))

        ! Write it out in the results array
        L2NormResults(6) = dintvalue

        ! Release everything
        call fev2_releaseVectorList(rcoeffVectors)
        call spdiscr_releaseCubStructure(rcubatureInfoPressureSecond)
        dintvalue = 0.0_DP

    end if

    ! Make some space on the Terminal
    call output_line("")

    !--------------------------------------!
    ! Now we take care of the point values !
    !--------------------------------------!

    ! Read the list of points
    ! First: find out how many points
    npoints = parlst_querysubstrings(rproblem_1%rparamlist,"POINTVALUESETTING", &
       "evaluationpoints")

    if (npoints .ge. 0) then
        ! Allocate everything we need
        allocate(pointValueDifference(npoints))
        allocate(pointvaluesFirst(npoints))
        allocate(pointvaluesSecond(npoints))
        allocate(Cder(npoints))
        allocate(IType(npoints))
        allocate(evaluationCoords(2,nPoints))

        ! and now read in the data
        do i=1,npoints
            call parlst_getvalue_string(rproblem_1%rparamlist,"POINTVALUESETTING","evaluationpoints", sparam,"",i)
            read (sparam,*) evaluationCoords(1,i), evaluationCoords(2,i), IType(i), Cder(i)
        end do

        ! If we want to calculate the difference we need
        ! to calculate both values.
        if (calcPointDifference .eq. 1) then
            calcPointFirst = 1
            calcPointSecond = 1
        end if

        ! We have everything, so lets go
        if (abs(calcPointFirst - 1.0_DP) .le. 0.1_DP) then

            ! Evaluate the function. This is done in another routine
            call extComparer_evaluate_function(rx1, evaluationCoords, Cder, IType, pointvaluesFirst)

            ! Print the results to the terminal
            call output_line("The point values of the first function are")
            call output_line("------------------------------------------")
            do i=1,npoints
                write (outputString,"(A10,A,F9.4,A,F9.4,A,E16.10)") Sfctnames(1+CDer(i),IType(i)),&
                        "(",evaluationCoords(1,i),",",evaluationCoords(2,i),") = ",pointvaluesFirst(i)
                call output_line(trim(outputString),coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
            end do

            !make some space on the terminal
            call output_line("")

        else
            ! If we did not calculate anything, we need a 0 in the array.
            ! This is necessary to have a know value in the array.
            ! If we leave this empty, the write-out routine will produce unpredictable
            ! results in the file, if we write a 0 we know what we end up with
            do i=1,nPoints
                pointvaluesFirst(i) = 0.0_DP
            end do
        end if

        if (abs(calcPointSecond - 1.0_DP) .le. 0.1_DP) then

            ! Evaluate the function. This is done in another routine.
            call extComparer_evaluate_function(rx2, evaluationCoords, Cder, IType, pointvaluesSecond)

            ! Print the results to the terminal
            call output_line("The point values of the second function are")
            call output_line("-------------------------------------------")
            do i=1,npoints
                write (outputString,"(A10,A,F9.4,A,F9.4,A,E16.10)") Sfctnames(1+CDer(i),IType(i)),&
                        "(",evaluationCoords(1,i),",",evaluationCoords(2,i),") = ",pointvaluesSecond(i)
                call output_line(trim(outputString),coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
            end do

            ! make some space on the terminal
            call output_line("")

        else
            ! If we did not calculate anything, we need a 0 in the array.
            ! This is necessary to have a know value in the array.
            ! If we leave this empty, the write-out routine will produce unpredictable
            ! results in the file, if we write a 0 we know what we end up with
            do i=1,nPoints
                pointvaluesSecond(i) = 0.0_DP
            end do
        end if

        ! and calc the difference
        if ( abs(calcPointDifference - 1.0_DP) .le. 0.1_DP ) then

            ! We have everything what we need - we just need to calculate the difference
            do i=1,nPoints
                pointValueDifference(i) = abs(pointvaluesFirst(i) - pointvaluesSecond(i))
            end do

            ! and print it out to the terminal
            call output_line("The absolute of the difference in the point values is")
            call output_line("-----------------------------------------------------")
            do i=1,npoints
                write (outputString,"(A10,A,F9.4,A,F9.4,A,E16.10)") Sfctnames(1+CDer(i),IType(i)),&
                        "(",evaluationCoords(1,i),",",evaluationCoords(2,i),") = ",pointValueDifference(i)
                call output_line(trim(outputString),coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
            end do

            ! make some space on the terminal
            call output_line("")

        else
            do i=1,nPoints
                pointValueDifference(i) = 0.0_DP
            end do
        end if

    end if !npoints

    !##################################################
    !# Postprocessing: write the results out in a file#
    !##################################################

    ! Write out the L2-Norm
    if( abs(writeOutL2Norm - 1.0_DP) .le. 0.1_DP ) then
        call vecio_writeArray(L2NormResults,0,outFileNameL2Norm,'(E20.10)')
    end if

    ! Write out point values

    if (npoints .ge. 0) then

        if (abs(writeOutPointValues - 1.0_DP) .le. 0.1_DP) then
            ! Write it out
           call io_openFileForWriting(outFileNamePointValues, iunit, &
           SYS_REPLACE, bfileExists,.true.)
           ! Write a Headline
           write (iunit,"(A)") "First function"
           do i=1,nPoints
            write (outputString,"(A10,A,F9.4,A,F9.4,A,E16.10)") Sfctnames(1+CDer(i),IType(i)),&
                 "(",evaluationCoords(1,i),",",evaluationCoords(2,i),") = ",pointvaluesFirst(i)
            write(iunit,*) outputString
           end do

           write(iunit,*) ""
           write (iunit,"(A)") "Second function"
           do i=1,nPoints
            write (outputString,"(A10,A,F9.4,A,F9.4,A,E16.10)") Sfctnames(1+CDer(i),IType(i)),&
                 "(",evaluationCoords(1,i),",",evaluationCoords(2,i),") = ",pointvaluesSecond(i)
            write(iunit,*) outputString
           end do

           write(iunit,*) ""
           write (iunit,"(A)") "The difference of the function values"
           do i=1,nPoints
            write (outputString,"(A10,A,F9.4,A,F9.4,A,E16.10)") Sfctnames(1+CDer(i),IType(i)),&
                 "(",evaluationCoords(1,i),",",evaluationCoords(2,i),") = ",pointValueDifference(i)
            write(iunit,*) outputString
           end do

          close(iunit)

        end if

        ! Point values written out, clean up
        deallocate(pointValueDifference)
        deallocate(pointvaluesFirst)
        deallocate(pointvaluesSecond)
        deallocate(Cder)
        deallocate(IType)
        deallocate(evaluationCoords)
    end if

  end subroutine



subroutine calc_L2error_velocity_2D_common(Dintvalue,rassemblyData,rintAssembly,&
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
    real(DP) :: dx,dy,dval1,dval2, work
    real(DP), dimension(1) :: dval_u21, dval_u22
    real(DP), dimension(2,1) :: evaluationPoint

    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:), pointer :: p_Dfunc11, p_Dfunc12
    real(DP), dimension(:,:,:), pointer :: p_Dpoints
    type(t_vectorBlock), pointer :: p_rx
    integer, dimension(:), pointer :: p_Ielements

    ! We need a parser to read in the values
    type(t_fparser) :: rparser
    ! storage for the parser
    character(256) :: domainDescription
    ! And we need the parameter list
    type(t_parlist) :: rparamlist


    ! Read out the area description in the string
    call parlst_getvalue_string(rcollection%p_rparlistQuickAccess1, 'DOMAININFO','common',domainDescription)

    call fparser_create(rparser,1)
    call fparser_parseFunction(rparser,1,domainDescription,(/'x','y'/))

    ! Now we can evaluate the area description with the parser-tools

    ! Get cubature weights
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Element numbers of the elements in process
    p_Ielements => rassemblyData%p_IelementList

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal

    ! Get a pointer to the values we already have
    p_Dfunc11 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC)
    p_Dfunc12 => revalVectors%p_RvectorData(2)%p_Ddata(:,:,DER_FUNC)

    ! and a pointer to the second function
    p_rx => rcollection%p_rvectorQuickAccess1

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

           call fparser_evalFunction(rparser,1,(/dx,dy/),work)

           ! Check if we need to work.
           ! If "work" is 1, we need to work. Since "work"
           ! can only be 0 or 1, this should run. However, you could
           ! hack this by defining functions that do not return 0 or 1
           ! as a function value, but. i.e. 0.5.
           if (abs(work - 1.0_DP) .le. 0.01_DP) then
                ! We need to work
                evaluationPoint(1,1) = dx
                evaluationPoint(2,1) = dy

                ! Evaluate the functions and store results in dval_u21 and dval_u22
                ! This will throw out an error if the domain-description is bad,
                ! to be precise if the domain description would force us to evaluate u2
                ! at non-mesh points!
                call fevl_evaluate(DER_FUNC, dval_u21, p_rx%RvectorBlock(1), evaluationPoint)
                call fevl_evaluate(DER_FUNC, dval_u22, p_rx%RvectorBlock(2), evaluationPoint)

                ! We did all evaluation
                dval1 = (p_Dfunc11(icubp,iel) - dval_u21(1) )**2
                dval2 = (p_Dfunc12(icubp,iel) - dval_u22(1) )**2

                Dintvalue(1) = Dintvalue(1) + (dval1 + dval2)*p_DcubWeight(icubp,iel)

          end if


        end do ! icubp

    end do ! iel


end subroutine

subroutine calc_L2norm_velocity_2D(Dintvalue,rassemblyData,rintAssembly,&
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
    real(DP) :: dx,dy,dval1,dval2, work


    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:), pointer :: p_Dfunc1, p_Dfunc2
    real(DP), dimension(:,:,:), pointer :: p_Dpoints
    integer, dimension(:), pointer :: p_Ielements

    ! We need a parser to read in the values
    type(t_fparser) :: rparser
    ! storage for the parser
    character(256) :: domainDescription
    ! And we need the parameter list
    type(t_parlist) :: rparamlist


    ! Read out the area description in the string
    call parlst_getvalue_string(rcollection%p_rparlistQuickAccess1, 'DOMAININFO', &
            rcollection%SquickAccess(1), domainDescription)

    call fparser_create(rparser,1)
    call fparser_parseFunction(rparser,1,domainDescription,(/'x','y'/))

    ! Now we can evaluate the area description with the parser-tools

    ! Get cubature weights
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Element numbers of the elements in process
    p_Ielements => rassemblyData%p_IelementList

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal

    ! Get a pointer to the values
    p_Dfunc1 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC)
    p_Dfunc2 => revalVectors%p_RvectorData(2)%p_Ddata(:,:,DER_FUNC)

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

           call fparser_evalFunction(rparser,1,(/dx,dy/),work)

           ! Check if we need to work.
           ! If "work" is 1, we need to work. Since "work"
           ! can only be 0 or 1, this should run. However, you could
           ! hack this by defining functions that do not return 0 or 1
           ! as a function value, but. i.e. 0.5.
           if (abs(work - 1.0_DP) .le. 0.01_DP) then
                ! We need to work
                dval1 = p_Dfunc1(icubp,iel)**2
                dval2 = p_Dfunc2(icubp,iel)**2
                Dintvalue(1) = Dintvalue(1) + (dval1 + dval2)*p_DcubWeight(icubp,iel)
          end if


        end do ! icubp

    end do ! iel


end subroutine


subroutine calc_L2error_pressure_2D_common(Dintvalue,rassemblyData,rintAssembly,&
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
    real(DP) :: dx,dy,dval1,dval2, work
    real(DP), dimension(1) :: dval_p21
    real(DP), dimension(2,1) :: evaluationPoint

    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:), pointer :: p_Dfunc11
    real(DP), dimension(:,:,:), pointer :: p_Dpoints
    type(t_vectorBlock), pointer :: p_rx
    integer, dimension(:), pointer :: p_Ielements

    ! We need a parser to read in the values
    type(t_fparser) :: rparser
    ! storage for the parser
    character(256) :: domainDescription
    ! And we need the parameter list
    type(t_parlist) :: rparamlist


    ! Read out the area description in the string
    call parlst_getvalue_string(rcollection%p_rparlistQuickAccess1, 'DOMAININFO','common',domainDescription)

    call fparser_create(rparser,1)
    call fparser_parseFunction(rparser,1,domainDescription,(/'x','y'/))

    ! Now we can evaluate the area description with the parser-tools

    ! Get cubature weights
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Element numbers of the elements in process
    p_Ielements => rassemblyData%p_IelementList

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal

    ! Get a pointer to the values we already have
    p_Dfunc11 => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC)

    ! and a pointer to the second function
    p_rx => rcollection%p_rvectorQuickAccess1

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

           call fparser_evalFunction(rparser,1,(/dx,dy/),work)

           ! Check if we need to work.
           ! If "work" is 1, we need to work. Since "work"
           ! can only be 0 or 1, this should run. However, you could
           ! hack this by defining functions that do not return 0 or 1
           ! as a function value, but. i.e. 0.5.
           if (abs(work - 1.0_DP) .le. 0.01_DP) then
                ! We need to work
                evaluationPoint(1,1) = dx
                evaluationPoint(2,1) = dy

                ! Evaluate the functions and store results in dval_u21 and dval_u22
                ! This will throw out an error if the domain-description is bad,
                ! to be precise if the domain description would force us to evaluate u2
                ! at non-mesh points!
                call fevl_evaluate(DER_FUNC, dval_p21, p_rx%RvectorBlock(3), evaluationPoint)

                ! We did all evaluation
                dval1 = (p_Dfunc11(icubp,iel) - dval_p21(1) )**2

                Dintvalue(1) = Dintvalue(1) + dval1*p_DcubWeight(icubp,iel)

          end if


        end do ! icubp

    end do ! iel


end subroutine

subroutine calc_L2norm_pressure_2D(Dintvalue,rassemblyData,rintAssembly,&
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
    real(DP) :: dx,dy,dval1,dval2, work


    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:), pointer :: p_Dfunc1
    real(DP), dimension(:,:,:), pointer :: p_Dpoints
    integer, dimension(:), pointer :: p_Ielements

    ! We need a parser to read in the values
    type(t_fparser) :: rparser
    ! storage for the parser
    character(256) :: domainDescription
    ! And we need the parameter list
    type(t_parlist) :: rparamlist


    ! Read out the area description in the string
    call parlst_getvalue_string(rcollection%p_rparlistQuickAccess1, 'DOMAININFO', &
            rcollection%SquickAccess(1), domainDescription)

    call fparser_create(rparser,1)
    call fparser_parseFunction(rparser,1,domainDescription,(/'x','y'/))

    ! Now we can evaluate the area description with the parser-tools

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

           call fparser_evalFunction(rparser,1,(/dx,dy/),work)

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

subroutine extComparer_evaluate_function(fefunction,Dpoints,Cder,IType,Dvalues )
        ! This is just a wrapper for the evaluation tools provided in the kernel.
        ! We need special settings for them, and it is more handy to do this
        ! with a wrapper like this.
        !<input>
        ! The function we want to evaluate
        type(t_vectorBlock), intent(in) :: fefunction

        ! List of points where we want to evaluate the function
        ! Dimension: (2,Number of Points)
        real(DP), dimension(:,:), intent(in) :: Dpoints

        ! For each point we specify if we want to evaluate
        ! the function or a derivative
        ! Size: Number of points
        integer, dimension(:), intent(in) :: Cder
        integer :: iderType

        ! We need to specify if we want to evaluate
        ! u or v or p
        integer, dimension(:), intent(in) :: IType

        !</input>

        !<output>
        ! The according function values
        real(DP), dimension(:), intent(out) :: Dvalues
        !</output>

        ! Local variables
        integer :: npoints
        integer :: i

        ! We somehow need the number of points
        npoints = ubound(Dpoints,2)

        ! Now we find out which derivative we want to evaluate
        do i=1,npoints
            select case (Cder(i))
                case (0)
                    iderType = DER_FUNC2D
                case (1)
                    iderType = DER_DERIV2D_X
                case (2)
                    iderType = DER_DERIV2D_Y
                case default
                    call output_line (&
                        "Unknown function derivative = "//sys_siL(Cder(i),10), &
                        OU_CLASS_ERROR,OU_MODE_STD,"extComparer_evaluate_function")
                        call sys_halt()
            end select

            ! And now we evaluate the function.
            ! since we evaluate both functions at the same points, it might happen that
            ! this leads the the case that we need a function value at a point where the
            ! function has no value or the code will crash. Therefore we assume 0 in that case.
            call fevl_evaluate(iderType,Dvalues(i:i),fefunction%RvectorBlock(IType(i)), &
                Dpoints(:,i:i), cnonmeshPoints=FEVL_NONMESHPTS_ZERO)
        end do

end subroutine

end module
