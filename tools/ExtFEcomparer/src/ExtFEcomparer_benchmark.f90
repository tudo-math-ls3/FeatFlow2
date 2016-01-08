! <description>
 ! This module has only one purpose:
 ! For debugging/regression-benchmark/... it is not too
 ! good to depend on data that has to be read in, i.e. if you
 ! just don't have a Q1-Solution for a grid, you cannot test the
 ! routines for Q1-vectors. For these tests the entries of the
 ! vector usually do not matter too much, we set up some vectors
 ! here, set up a postprocessing structure and then run the usual
 ! routines.
 !
 ! The first digit of the testID is the dimension!
! </description>

module ExtFEcomparer_benchmark

  ! What do we need?
  ! Include basic Feat-2 modules
  use genoutput
  use paramlist
  use collection
  use storage

  use ExtFEcomparer_typedefs

  use ExtFEcomparer_parameters
  use ExtFEcomparer_discretisation
  use ExtFEcomparer_paramtriang
  use ExtFEcomparer_core
  use ExtFEcomparer_vector
  use ExtFEcomparer_postprocessing

  use triangulation
  use linearsystemblock
  use linearsystemscalar
  use element
  use meshgeneration

  use ucd

  implicit none

  private

  public :: ExtFEcomparer_benchMod

  contains


  subroutine ExtFEcomparer_benchMod(modus)
    !<input>
    ! If modus == 1: no VTK output
    ! If modus == 2: VTK-output, good for debugging
    integer, intent(in) :: modus
    !</input>
    ! we need 2 problem structures + a posptrocessing structure
    type(t_problem), pointer :: p_rproblem1, p_rproblem2

    type(t_postprocessing), pointer :: p_rpostprocessing

    ! Tell the user we are in the benchmark modus
    call output_separator(OU_SEP_AT)
    call output_separator(OU_SEP_DOLLAR)
    call output_line("$ Starting in the bechmark modus")
    call output_line("$ All input functions are ignored")
    call output_line("$ Instead, we use the benchmark functions")
    call output_line("$ and perform the benchmark calculations")
    call output_line("$ for the test that was chosen on the commandline")
    call output_separator(OU_SEP_DOLLAR)
    call output_separator(OU_SEP_AT)

    allocate(p_rproblem1)
    allocate(p_rproblem2)
    allocate(p_rpostprocessing)

    ! We need the parlists
    call parlst_init(p_rproblem1%rparamlist)
    call parlst_init(p_rproblem2%rparamlist)


    ! Set up the tests
    call ExtFE_setup_test(p_rproblem1,&
                          p_rproblem2,&
                          p_rpostprocessing,&
                          modus)

    ! We have everything, so lets calculate
    call ExtFEcomparer_calculate(p_rproblem1,&
                                 p_rproblem2,&
                                 p_rpostprocessing)

    ! We even might want to do some output - ie when debugging
    call ExtFEcomparer_postprocess(p_rpostprocessing)

    ! General clean-up
    call parlst_done(p_rproblem1%rparamlist)
    call parlst_done(p_rproblem2%rparamlist)

    call ExtFEcomparer_doneDiscretisation(p_rproblem1)
    call ExtFEcomparer_doneDiscretisation(p_rproblem2)

    call ExtFEcomparer_doneParamTriang(p_rproblem1)
    call ExtFEcomparer_doneParamTriang(p_rproblem2)

    call ExtFEcomparer_done_postprocessing(p_rpostprocessing)

    call lsysbl_releaseVector(p_rproblem1%coeffVector)
    call lsysbl_releaseVector(p_rproblem2%coeffVector)

    if(associated(p_rproblem1%iElemList)) then
        deallocate(p_rproblem1%iElemList)
    end if
    if(associated(p_rproblem2%iElemList)) then
        deallocate(p_rproblem2%iElemList)
    end if
    deallocate(p_rproblem1)
    deallocate(p_rproblem2)
    deallocate(p_rpostprocessing)

  end subroutine




subroutine ExtFE_setup_test(rproblem1,rproblem2,rpostprocessing,modus)
    ! <inout>
    type(t_problem) , intent(inout):: rproblem1
    type(t_problem) , intent(inout):: rproblem2
    type(t_postprocessing), intent(inout) :: rpostprocessing
    integer, intent(in) :: modus
    ! </inout>


    ! local variables
    logical :: bexists
    character(LEN=ExtFE_STRLEN) :: smaster
    integer :: testID,i
    real(DP), dimension(:), pointer :: p_VecEntries => NULL()

    ! We might want to do some UCD-Output - depends on the modus
    integer, dimension(:), pointer :: p_FuncAddType1
    integer, dimension(:), pointer :: p_FuncAddType2
    integer, dimension(:), pointer :: p_FuncProjType1
    integer, dimension(:), pointer :: p_FuncProjType2
    integer, dimension(:), pointer :: p_FuncScalarVars1
    integer, dimension(:), pointer :: p_FuncScalarVars2


    ! At first, we parse all arguments we have to the parlist
    ! As we want to keep the benchmark modus "hidden" we will not
    ! have a standard place to search, so all arguments are on
    ! the command line. However, there can be a dat-file.
    ! We call it master-dat here
    smaster = ''
    call ExtFEcomparer_getCmdlMasterDat(smaster)
    ! Read the file "master.dat" if it exists
    inquire(file=smaster, exist=bexists)

    if (bexists) then
      ! Read the master file. That either one contains all parameters or
      ! contains references to subfiles with data.
      call parlst_readfromfile (rproblem1%rparamlist, smaster)
      call parlst_readfromfile (rproblem2%rparamlist, smaster)
    end if

    ! Now we add the parameters from the command-line
    ! Thus, command-lines overrules file!
    call ExtFEcomparer_parseCmdlArguments(rproblem1%rparamlist)
    call ExtFEcomparer_parseCmdlArguments(rproblem2%rparamlist)

    ! Print the actual parameters on the screen
    call parlst_info(rproblem1%rparamlist)

    ! Figure out the testID. It is in the Section
    ! ExtFE-Benchmark
    call parlst_getvalue_int(rproblem1%rparamlist,"ExtFE-Benchmark", &
                            "testID",testID)

    ! now, according to the testID, we set up different
    ! geometries and different elements, but the rest
    ! stays the same
    select case(testID)

        case(101)
            call ExtFE_setup_problem_test101(rproblem1)
            call ExtFE_setup_problem_test101(rproblem2)
            call ExtFE_setup_calc_test101(rproblem1,rproblem2,rpostprocessing)
        case(201)
            call ExtFE_setup_problem_test201(rproblem1)
            call ExtFE_setup_problem_test201(rproblem2)
            call ExtFE_setup_calc_test201(rproblem1,rproblem2,rpostprocessing)

        case default
            call output_line("Invalid TestID: "//trim(adjustl(sys_siL(testID,5))) , &
            OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_benchmark")
            call sys_halt()

    end select

    ! The tests may differ in Number of elements/mesh type,
    ! but have 1 thing in common: Output is always
    ! a problem structure with following valid fields:
    ! rproblem%ielementSetting
    ! rproblem%ielemlist/ielempairID/ielemType
    ! rproblem%triangulation
    ! rproblem%NVAR
    ! rproblem%iDimension

    ! With this information, our discretisation routine should work
    call ExtFEcomparer_init_discretisation(rproblem1)
    call ExtFEcomparer_init_discretisation(rproblem2)

    ! After the part above us, we now should have a valid
    ! block-discretisation
    ! so we create a vector according to the block-discretisation
    call lsysbl_createVecBlockByDiscr(rproblem1%rdiscretisation,&
                rproblem1%coeffVector,bclear=.TRUE.)
    call lsysbl_createVecBlockByDiscr(rproblem2%rdiscretisation,&
                rproblem2%coeffVector,bclear=.TRUE.)


    ! Now we have to fill the vector
    select case(testID)
        case(101,201)
            call ExtFE_fill_vec_ones(rproblem1)
            call ExtFE_fill_vec_ones(rproblem2)
        case default
            call output_line("No strategy to fill the vec for TestID: "//trim(adjustl(sys_siL(testID,5))) , &
            OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_benchmark")
            call sys_halt()
    end select

    ! Depending on the modus we do ucd-output
    if(modus .eq. 2) then
        ! To debug we might want to do some UCD-Output
        rpostprocessing%ucd_OUT_meshes = .TRUE.
        rpostprocessing%UCD_Format = ExtFE_UCD_VTK
        rpostprocessing%UCD_Style = UCD_FLAG_STANDARD
        rpostprocessing%UCD_meshOneOutPath = "./benchout/benchmesh1.vtk"
        rpostprocessing%UCD_meshTwoOutPath = "./benchout/benchmesh2.vtk"
        rpostprocessing%UCD_MeshOnePointer => rproblem1%coeffVector%p_rblockDiscr%p_rtriangulation
        rpostprocessing%UCD_MeshTwoPointer => rproblem2%coeffVector%p_rblockDiscr%p_rtriangulation
        rpostprocessing%ucd_OUT_orig_functions_one = .TRUE.
        rpostprocessing%ucd_OUT_orig_functions_two = .TRUE.
        rpostprocessing%UCD_FEfunctionOneOrigOutPath = "./benchout/benchfunction1.vtk"
        rpostprocessing%UCD_FEfunctionTwoOrigOutPath = "./benchout/benchfunction2.vtk"
        ! The add-type (element-based or vertex-based)?
        ! We only add scalar vars - no need for vectors
        call storage_new("ExtFE_benchmark","FuncOneAddType",&
                        rproblem1%NVAR,ST_INT,rpostprocessing%h_UCD_AddTypeOrigFirst,&
                        ST_NEWBLOCK_ZERO)
        ! Project them linear or constant?
        call storage_new("ExtFE_benchmark","FuncOneProjections",&
                        rproblem1%NVAR,ST_INT,rpostprocessing%h_UCD_AddElemProjectFirst,&
                        ST_NEWBLOCK_ZERO)
        call storage_new("ExtFE_benchmark","FuncOneScalars",&
                    rproblem1%NVAR,ST_INT,rpostprocessing%h_UCD_ScalarFirstOrig,&
                    ST_NEWBLOCK_ZERO)
        ! Second function
        call storage_new("ExtFE_benchmark","FuncTwoAddType",&
                        rproblem2%NVAR,ST_INT,rpostprocessing%h_UCD_AddTypeOrigSecond,&
                        ST_NEWBLOCK_ZERO)
        ! Project them linear or constant?
        call storage_new("ExtFE_benchmark","FuncTwoProjections",&
                        rproblem2%NVAR,ST_INT,rpostprocessing%h_UCD_AddElemProjectSecond,&
                        ST_NEWBLOCK_ZERO)
        call storage_new("ExtFE_benchmark","FuncTwoScalars",&
                    rproblem2%NVAR,ST_INT,rpostprocessing%h_UCD_ScalarSecondOrig,&
                    ST_NEWBLOCK_ZERO)
        call storage_getbase_int(rpostprocessing%h_UCD_AddElemProjectFirst,p_FuncProjType1)
        call storage_getbase_int(rpostprocessing%h_UCD_AddElemProjectSecond,p_FuncProjType2)
        call storage_getbase_int(rpostprocessing%h_UCD_AddTypeOrigFirst,p_FuncAddType1)
        call storage_getbase_int(rpostprocessing%h_UCD_AddTypeOrigSecond,p_FuncAddType2)
        call storage_getbase_int(rpostprocessing%h_UCD_ScalarFirstOrig,p_FuncScalarVars1)
        call storage_getbase_int(rpostprocessing%h_UCD_ScalarSecondOrig,p_FuncScalarVars2)

        do i=1,rproblem1%NVAR
            p_FuncAddType1(i) = ExtFE_UCD_ELEM_BASED
            p_FuncProjType1(i) = ExtFE_UCD_POLY_CONST
            p_FuncScalarVars1(i) = i

            p_FuncAddType2(i) = ExtFE_UCD_VERT_BASED
            p_FuncProjType2(i) = ExtFE_UCD_POLY_LINEAR
            p_FuncScalarVars2(i) = i
        end do
    end if


    ! all things are set up - so we can start to calculate!
end subroutine


!-----------------------------------------------------
! Here, the routines that set up a test follow!


! In this subroutine we set up TestID 101:
! configuration: 1D, mesh: [0,1], different element types
subroutine ExtFE_setup_problem_test101(rproblem)
    !<inputoutput>
    type(t_problem), intent(inout) :: rproblem
    !</inputoutput>

    ! How many elements?
    integer :: nElemTest

    ! How often will we refine the mesh
    integer :: nrefinements

    ! 1D!
    rproblem%iDimension = ExtFE_NDIM1

    ! Create a mesh: [0,1]
    call tria_createRawTria1D(rproblem%rtriangulation,0.0_DP,1.0_DP,4)

    ! and refine it
    nrefinements = 6
    call tria_quickRefine2LevelOrdering(nrefinements,rproblem%rtriangulation)

    ! Init a standard mesh
    call tria_initStandardMeshFromRaw(rproblem%rtriangulation)

    ! For the different elements we will provide a list
    rproblem%elementSetting = ExtFE_ElementList

    ! We need to count manually how many elements we test
    nElemTest = 7

    ! For each element 1 component
    ! Do not change this line! Only change nElemTest!
    rproblem%NVAR = nElemTest

    ! Since the list is allocatable, we have to allocate
    ! it first. Therefore you have to count first
    allocate(rproblem%iElemList(nElemTest))

    rproblem%iElemList(1) = EL_P0_1D
    rproblem%iElemList(2) = EL_P1_1D
    rproblem%iElemList(3) = EL_P2_1D
    rproblem%iElemList(4) = EL_S31_1D
    rproblem%iElemList(5) = EL_DG_T0_1D
    rproblem%iElemList(6) = EL_DG_T1_1D
    rproblem%iElemList(7) = EL_DG_T2_1D

end subroutine


! Here, we set up the calculations for test 101:
! At first, we calculate ||f1 - f2|| for every component on the whole
! domain. Then, we calculate ||f_1|| on the subset with x<0.5
! and ||f_1|| on the whole domain, then ||f_2|| on x<0.5 and ||f_2|| on
! the whole domain
subroutine ExtFE_setup_calc_test101(rproblem1,rproblem2,rpostprocessing)
    ! <input>
    type(t_problem) , intent(inout):: rproblem1
    type(t_problem) , intent(inout):: rproblem2
    ! <input>
    ! <inout>
    type(t_postprocessing), intent(inout) :: rpostprocessing
    ! </inout>

    ! Local variables:
    integer :: nL2calculations,nPointCalcs,nL1calculations,nIntCalculations
    integer :: i
    character(LEN=ExtFE_STRLEN), dimension(:), allocatable :: sParseVariables
    character(LEN=ExtFE_STRLEN) :: L2RegionOfInterest1,L2RegionOfInterest2
    character(LEN=ExtFE_STRLEN) :: L1RegionOfInterest1, L1RegionOfInterest2
    character(LEN=ExtFE_STRLEN) :: IntRegionOfInterest1, IntRegionOfInterest2

    ! To get the arrays we need these pointers
    integer, dimension(:,:), pointer :: p_L2FuncComp, p_L1FuncComp, p_IntFuncComp
    integer(I32), dimension(:), pointer :: p_L2CubRule, p_L1CubRule, p_IntCubRule

    ! a pointer to the component + derivative array
    integer, dimension(:,:), pointer :: p_FuncComp
    ! a pointer to the array for the evaluation point
    real(DP), dimension(:,:), pointer :: p_EvalPoint


    ! We test the Integral-Calculations
    ! How many tests do we do?
    ! first Integral-Series: f_1(i) - f_2(i) on the entire domain
    ! second Series: f_1(i) on half of the domain
    ! third series: f_1(i) on the full domain
    ! 4th Series: f_2(i) on x<0.5
    ! 5th series: f_2(i) on the whole domain
    ! The first test should always return something close to 0
    ! the second and third one together verify that parsing
    ! still works. So we have 5*NVAR tests
    nIntcalculations = 5*rproblem1%NVAR
    rpostprocessing%nIntCalculations = nIntcalculations

    ! Init the arrays
    ! Array for components
    call storage_new("ExtFE_setup_tests", &
            "IntCompFunc", (/2,nIntCalculations/),ST_INT, &
             rpostprocessing%h_IntCompFunc,ST_NEWBLOCK_ZERO)

    ! Array for the results
    call storage_new("ExtFE_setup_tests", &
            "IntResults", nIntCalculations,ST_DOUBLE, &
             rpostprocessing%h_IntResults,ST_NEWBLOCK_ZERO)

    ! Array for the cubature rules
    call storage_new("ExtFE_setup_tests", &
            "IntCubRule",nIntCalculations,ST_INT32, &
            rpostprocessing%h_IntCubRule, ST_NEWBLOCK_ZERO)


    ! We will actually parse the region of interest here in a parser
    ! In component i is the region of interest for calculation i
    call fparser_create(rpostprocessing%pIntChiOmegaParser,nIntCalculations)

    ! Now get the arrays and fill them with data
    call storage_getbase_int2D(rpostprocessing%h_IntCompFunc,p_IntFuncComp)
    call storage_getbase_int32(rpostprocessing%h_IntCubRule,p_IntCubRule)


    IntRegionOfInterest1 = '1'
    IntRegionOfInterest2 = 'IF(x<=0.5,1,0)'
    do i=1,rproblem1%NVAR
        ! First series of test
        ! cub rule
        p_IntCubRule(5*i-4) = CUB_G3_1D
        ! components
        p_IntFuncComp(1,5*i-4) = i
        p_IntFuncComp(2,5*i-4) = i
        ! Parse the region of interest
        call fparser_parseFunction(rpostprocessing%pIntChiOmegaParser,5*i-4,&
                    IntRegionOfInterest1,(/'x'/))

        ! Second series of tests
        p_IntCubRule(5*i -3) = CUB_G3_1D
        p_IntFuncComp(1,5*i -3) = i
        p_IntFuncComp(2,5*i -3) = -1
        call fparser_parseFunction(rpostprocessing%pIntChiOmegaParser,5*i-3,&
                    IntRegionOfInterest2,(/'x'/))

        ! Third series of tests
        p_IntCubRule(5*i-2) = CUB_G3_1D
        p_IntFuncComp(1,5*i-2) = i
        p_IntFuncComp(2,5*i-2) = -1
        call fparser_parseFunction(rpostprocessing%pIntChiOmegaParser,5*i-2,&
                    IntRegionOfInterest1,(/'x'/))


        ! 4th series of tests
        p_IntCubRule(5*i-1) = CUB_G3_1D
        p_IntFuncComp(1,5*i-1) = -1
        p_IntFuncComp(2,5*i-1) = i
        call fparser_parseFunction(rpostprocessing%pIntChiOmegaParser,5*i-1,&
                    IntRegionOfInterest2,(/'x'/))

        ! 5th series of tests
        p_IntCubRule(5*i) = CUB_G3_1D
        p_IntFuncComp(1,5*i) = -1
        p_IntFuncComp(2,5*i) = i
        call fparser_parseFunction(rpostprocessing%pIntChiOmegaParser,5*i,&
                    IntRegionOfInterest1,(/'x'/))

    end do




    ! We test the L1-Calculations
    ! How many tests do we do?
    ! first L1-Series: ||f_1(i) - f_2(i)||_L2 on the entire domain
    ! second L1-Series: ||f_1(i)|| on half of the domain
    ! third L1-series: ||f_1(i)|| on the full domain
    ! 4th L1-Series: ||f_2(i)|| on x<0.5
    ! 5th L1-series: ||f_2(i)|| on the whole domain
    ! The first test should always return something close to 0
    ! the second and third one together verify that parsing
    ! still works. So we have 5*NVAR tests
    nL1calculations = 5*rproblem1%NVAR
    rpostprocessing%nL1Calculations = nL1calculations

    ! Init the arrays
    ! Array for components
    call storage_new("ExtFE_setup_tests", &
            "L1CompFunc", (/2,nL1Calculations/),ST_INT, &
             rpostprocessing%h_L1CompFunc,ST_NEWBLOCK_ZERO)

    ! Array for the results
    call storage_new("ExtFE_setup_tests", &
            "L1Results", nL1Calculations,ST_DOUBLE, &
             rpostprocessing%h_L1Results,ST_NEWBLOCK_ZERO)

    ! Array for the cubature rules
    call storage_new("ExtFE_setup_tests", &
            "L1CubRule",nL1Calculations,ST_INT32, &
            rpostprocessing%h_L1CubRule, ST_NEWBLOCK_ZERO)


    ! We will actually parse the L1RegionOfInterest here in a parser
    ! In component i is the region of interest for calculation i
    call fparser_create(rpostprocessing%pL1ChiOmegaParser,nL1Calculations)

    ! Now get the arrays and fill them with data
    call storage_getbase_int2D(rpostprocessing%h_L1CompFunc,p_L1FuncComp)
    call storage_getbase_int32(rpostprocessing%h_L1CubRule,p_L1CubRule)


    L1RegionOfInterest1 = '1'
    L1RegionOfInterest2 = 'IF(x<=0.5,1,0)'
    do i=1,rproblem1%NVAR
        ! First series of test
        ! cub rule
        p_L1CubRule(5*i-4) = CUB_G3_1D
        ! components
        p_L1FuncComp(1,5*i-4) = i
        p_L1FuncComp(2,5*i-4) = i
        ! Parse the region of interest
        call fparser_parseFunction(rpostprocessing%pL1ChiOmegaParser,5*i-4,&
                    L1RegionOfInterest1,(/'x'/))

        ! Second series of tests
        p_L1CubRule(5*i -3) = CUB_G3_1D
        p_L1FuncComp(1,5*i -3) = i
        p_L1FuncComp(2,5*i -3) = -1
        call fparser_parseFunction(rpostprocessing%pL1ChiOmegaParser,5*i-3,&
                    L1RegionOfInterest2,(/'x'/))

        ! Third series of tests
        p_L1CubRule(5*i-2) = CUB_G3_1D
        p_L1FuncComp(1,5*i-2) = i
        p_L1FuncComp(2,5*i-2) = -1
        call fparser_parseFunction(rpostprocessing%pL1ChiOmegaParser,5*i-2,&
                    L1RegionOfInterest1,(/'x'/))


        ! 4th series of tests
        p_L1CubRule(5*i-1) = CUB_G3_1D
        p_L1FuncComp(1,5*i-1) = -1
        p_L1FuncComp(2,5*i-1) = i
        call fparser_parseFunction(rpostprocessing%pL1ChiOmegaParser,5*i-1,&
                    L1RegionOfInterest2,(/'x'/))

        ! 5th series of tests
        p_L1CubRule(5*i) = CUB_G3_1D
        p_L1FuncComp(1,5*i) = -1
        p_L1FuncComp(2,5*i) = i
        call fparser_parseFunction(rpostprocessing%pL1ChiOmegaParser,5*i,&
                    L1RegionOfInterest1,(/'x'/))

    end do



    ! We test the L2-Calculations
    ! How many tests do we do?
    ! first L2-Series: ||f_1(i) - f_2(i)||_L2 on the entire domain
    ! second L2-Series: ||f_1(i)|| on half of the domain
    ! third L2-series: ||f_1(i)|| on the full domain
    ! 4th L2-Series: ||f_2(i)|| on x<0.5
    ! 5th L2-series: ||f_2(i)|| on the whole domain
    ! The first test should always return something close to 0
    ! the second and third one together verify that parsing
    ! still works. So we have 5*NVAR tests
    nL2calculations = 5*rproblem1%NVAR
    rpostprocessing%nL2Calculations = nL2calculations

    ! Init the arrays
    ! Array for components
    call storage_new("ExtFE_setup_tests", &
            "L2CompFunc", (/2,nL2Calculations/),ST_INT, &
             rpostprocessing%h_L2CompFunc,ST_NEWBLOCK_ZERO)

    ! Array for the results
    call storage_new("ExtFE_setup_tests", &
            "L2Results", nL2Calculations,ST_DOUBLE, &
             rpostprocessing%h_L2Results,ST_NEWBLOCK_ZERO)

    ! Array for the cubature rules
    call storage_new("ExtFE_setup_tests", &
            "L2CubRule",nL2Calculations,ST_INT32, &
            rpostprocessing%h_L2CubRule, ST_NEWBLOCK_ZERO)


    ! We will actually parse the L2RegionOfInterest here in a parser
    ! In component i is the region of interest for calculation i
    call fparser_create(rpostprocessing%pL2ChiOmegaParser,nL2Calculations)

    ! Now get the arrays and fill them with data
    call storage_getbase_int2D(rpostprocessing%h_L2CompFunc,p_L2FuncComp)
    call storage_getbase_int32(rpostprocessing%h_L2CubRule,p_L2CubRule)


    L2RegionOfInterest1 = '1'
    L2RegionOfInterest2 = 'IF(x<=0.5,1,0)'
    do i=1,rproblem1%NVAR
        ! First series of test
        ! cub rule
        p_L2CubRule(5*i-4) = CUB_G3_1D
        ! components
        p_L2FuncComp(1,5*i-4) = i
        p_L2FuncComp(2,5*i-4) = i
        ! Parse the region of interest
        call fparser_parseFunction(rpostprocessing%pL2ChiOmegaParser,5*i-4,&
                    L2RegionOfInterest1,(/'x'/))

        ! Second series of tests
        p_L2CubRule(5*i -3) = CUB_G3_1D
        p_L2FuncComp(1,5*i -3) = i
        p_L2FuncComp(2,5*i -3) = -1
        call fparser_parseFunction(rpostprocessing%pL2ChiOmegaParser,5*i-3,&
                    L2RegionOfInterest2,(/'x'/))

        ! Third series of tests
        p_L2CubRule(5*i-2) = CUB_G3_1D
        p_L2FuncComp(1,5*i-2) = i
        p_L2FuncComp(2,5*i-2) = -1
        call fparser_parseFunction(rpostprocessing%pL2ChiOmegaParser,5*i-2,&
                    L2RegionOfInterest1,(/'x'/))


        ! 4th series of tests
        p_L2CubRule(5*i-1) = CUB_G3_1D
        p_L2FuncComp(1,5*i-1) = -1
        p_L2FuncComp(2,5*i-1) = i
        call fparser_parseFunction(rpostprocessing%pL2ChiOmegaParser,5*i-1,&
                    L2RegionOfInterest2,(/'x'/))

        ! 5th series of tests
        p_L2CubRule(5*i) = CUB_G3_1D
        p_L2FuncComp(1,5*i) = -1
        p_L2FuncComp(2,5*i) = i
        call fparser_parseFunction(rpostprocessing%pL2ChiOmegaParser,5*i,&
                    L2RegionOfInterest1,(/'x'/))



    end do

    ! We also want to test the point-evaluation on their own!

    ! How many tests do we want to do?
    ! We want to calculate f-g, (d/dx)f - (d/dx)g, f, (d/dx)f
    nPointCalcs =4*rproblem1%NVAR

    rpostprocessing%nPointCalculations = nPointCalcs

    call storage_new("ExtFEcomparer_benchmark", &
            "PointvalueResults",nPointCalcs,ST_DOUBLE,&
                rpostprocessing%h_PointResults,ST_NEWBLOCK_ZERO)

     ! PointComponent(1,i) = component of function 1 in calculation i
     ! PointComponent(2,i) = derivative of component of function 1 in calculation i
     ! PointComponent(3,i) = component of function 2 in calculation i
     ! PointComponent(4,i) = derivative of component of function 2 in calculation i
     call storage_new("ExtFEcomparer_benchmark", &
                "PointFuncComponents",(/4,nPointCalcs/),ST_INT, &
                 rpostprocessing%h_PointFuncComponents,ST_NEWBLOCK_ZERO)

     ! PointCoordinates(:,i) = coordinates in calculation i
     ! We are in a 1D-Routine!
     call storage_new("ExtFEcomparer_init_postprocessing", &
                "PointCoordinates",(/ExtFE_NDIM1,nPointCalcs/),&
                ST_DOUBLE,rpostprocessing%h_PointCoordinates, &
                ST_NEWBLOCK_ZERO)

    call storage_getbase_int2D(rpostprocessing%h_PointFuncComponents,p_FuncComp)
    call storage_getbase_double2D(rpostprocessing%h_PointCoordinates,p_EvalPoint)

    ! We want to calculate f-g, (d/dx)f - (d/dx)g, f, (d/dx)f
    do i=1,rproblem1%NVAR
        ! First one - f-g
        p_FuncComp(1,4*i-3) = i
        p_FuncComp(2,4*i-3) = ExtFE_DER_FUNC_1D
        p_FuncComp(3,4*i-3) = i
        p_FuncComp(4,4*i-3) = ExtFE_DER_FUNC_1D
        p_EvalPoint(1,4*i-3) = 0.3_DP

        ! derivatives
        p_FuncComp(1,4*i-2) = i
        p_FuncComp(2,4*i-2) = ExtFE_DER_DERIV_1D_X
        p_FuncComp(3,4*i-2) = i
        p_FuncComp(4,4*i-2) = ExtFE_DER_DERIV_1D_X
        p_EvalPoint(1,4*i-2) = 0.3_DP

        ! only first
        p_FuncComp(1,4*i-1) = i
        p_FuncComp(2,4*i-1) = ExtFE_DER_FUNC_1D
        p_FuncComp(3,4*i-1) = -1
        p_FuncComp(4,4*i-1) = ExtFE_DER_FUNC_1D
        p_EvalPoint(1,4*i-1) = 0.3_DP

        ! only first - derivative
        p_FuncComp(1,4*i) = i
        p_FuncComp(2,4*i) = ExtFE_DER_DERIV_1D_X
        p_FuncComp(3,4*i) = -1
        p_FuncComp(4,4*i) = ExtFE_DER_FUNC_1D
        p_EvalPoint(1,4*i) = 0.3_DP

    end do


end subroutine


!----------------------------------------------------------------------


! In this subroutine we set up TestID 201:
! configuration: 2D, mesh with quads [0,1]^2, different element types
subroutine ExtFE_setup_problem_test201(rproblem)
    !<inputoutput>
    type(t_problem), intent(inout) :: rproblem
    !</inputoutput>

    ! How many elements?
    integer :: nElemTest

    ! How often will we refine the mesh
    integer :: nrefinements

    ! 2D!
    rproblem%iDimension = ExtFE_NDIM2

    ! Create a mesh: [0,1]^2
    call meshgen_rectangular2DQuadMesh (rproblem%rtriangulation, 0.0_DP, 1.0_DP, 0.0_DP, 1.0_DP, 4, 4)

    ! and refine it
    nrefinements = 4
    call tria_quickRefine2LevelOrdering(nrefinements,rproblem%rtriangulation)

    ! Init a standard mesh
    call tria_initStandardMeshFromRaw(rproblem%rtriangulation)

    ! For the different elements we will provide a list
    rproblem%elementSetting = ExtFE_ElementList

    ! We need to count manually how many elements we test
    nElemTest = 14

    ! For each element 1 component
    ! Do not change this line! Only change nElemTest!
    rproblem%NVAR = nElemTest

    ! Since the list is allocatable, we have to allocate
    ! it first. Therefore you have to count first
    allocate(rproblem%iElemList(nElemTest))

    rproblem%iElemList(1) = EL_Q0_2D
    rproblem%iElemList(2) = EL_Q1_2D
    rproblem%iElemList(3) = EL_Q2_2D
    rproblem%iElemList(4) = EL_Q3_2D
    rproblem%iElemList(5) = EL_QPW4P0_2D
    rproblem%iElemList(6) = EL_QPW4P1_2D
    rproblem%iElemList(7) = EL_QPW4P2_2D
    rproblem%iElemList(8) = EL_QPW4P1T_2D
    rproblem%iElemList(9) = EL_QP1_2D
    rproblem%iElemList(10) = EL_QP1NPD_2D
    rproblem%iElemList(11) = EL_Q1T_2D
    rproblem%iElemList(12) = EL_EM31_2D
    rproblem%iElemList(13) = EL_Q1NP_2D
    rproblem%iElemList(14) = EL_EM30_2D

end subroutine


subroutine ExtFE_setup_calc_test201(rproblem1,rproblem2,rpostprocessing)
    ! <input>
    type(t_problem) , intent(inout):: rproblem1
    type(t_problem) , intent(inout):: rproblem2
    ! <input>
    ! <inout>
    type(t_postprocessing), intent(inout) :: rpostprocessing
    ! </inout>


    ! Local variables:
    integer :: nL2calculations, nIntCalculations, nL1calculations
    integer :: i
    character(LEN=ExtFE_STRLEN), dimension(:), allocatable :: sParseVariables
    character(LEN=ExtFE_STRLEN) :: L2RegionOfInterest1, L2RegionOfInterest2
    character(LEN=ExtFE_STRLEN) :: L1RegionOfInterest1, L1RegionOfInterest2
    character(LEN=ExtFE_STRLEN) :: IntRegionOfInterest1, IntRegionOfInterest2

    ! To get the arrays we need these pointers
    integer, dimension(:,:), pointer :: p_L2FuncComp, p_L1FuncComp, p_IntFuncComp
    integer(I32), dimension(:), pointer :: p_L2CubRule, p_L1CubRule, p_IntCubRule


    ! Local variables:
    integer :: nPointCalcs
    ! a pointer to the component + derivative array
    integer, dimension(:,:), pointer :: p_FuncComp
    ! a pointer to the array for the evaluation point
    real(DP), dimension(:,:), pointer :: p_EvalPoint



    ! We test the Integral-Calculations
    ! How many tests do we do?
    ! first Integral-Series: f_1(i) - f_2(i) on the entire domain
    ! second Series: f_1(i) on half of the domain
    ! third series: f_1(i) on the full domain
    ! 4th Series: f_2(i) on x<0.5
    ! 5th series: f_2(i) on the whole domain
    ! The first test should always return something close to 0
    ! the second and third one together verify that parsing
    ! still works. So we have 5*NVAR tests
    nIntcalculations = 5*rproblem1%NVAR
    rpostprocessing%nIntCalculations = nIntcalculations

    ! Init the arrays
    ! Array for components
    call storage_new("ExtFE_setup_tests", &
            "IntCompFunc", (/2,nIntCalculations/),ST_INT, &
             rpostprocessing%h_IntCompFunc,ST_NEWBLOCK_ZERO)

    ! Array for the results
    call storage_new("ExtFE_setup_tests", &
            "IntResults", nIntCalculations,ST_DOUBLE, &
             rpostprocessing%h_IntResults,ST_NEWBLOCK_ZERO)

    ! Array for the cubature rules
    call storage_new("ExtFE_setup_tests", &
            "IntCubRule",nIntCalculations,ST_INT32, &
            rpostprocessing%h_IntCubRule, ST_NEWBLOCK_ZERO)


    ! We will actually parse the region of interest here in a parser
    ! In component i is the region of interest for calculation i
    call fparser_create(rpostprocessing%pIntChiOmegaParser,nIntCalculations)

    ! Now get the arrays and fill them with data
    call storage_getbase_int2D(rpostprocessing%h_IntCompFunc,p_IntFuncComp)
    call storage_getbase_int32(rpostprocessing%h_IntCubRule,p_IntCubRule)


    IntRegionOfInterest1 = '1'
    IntRegionOfInterest2 = 'IF(x<=0.5,1,0)'
    do i=1,rproblem1%NVAR
        ! First series of test
        ! cub rule
        p_IntCubRule(5*i-4) = CUB_G3_2D
        ! components
        p_IntFuncComp(1,5*i-4) = i
        p_IntFuncComp(2,5*i-4) = i
        ! Parse the region of interest
        call fparser_parseFunction(rpostprocessing%pIntChiOmegaParser,5*i-4,&
                    IntRegionOfInterest1,(/'x','y'/))

        ! Second series of tests
        p_IntCubRule(5*i -3) = CUB_G3_2D
        p_IntFuncComp(1,5*i -3) = i
        p_IntFuncComp(2,5*i -3) = -1
        call fparser_parseFunction(rpostprocessing%pIntChiOmegaParser,5*i-3,&
                    IntRegionOfInterest2,(/'x','y'/))

        ! Third series of tests
        p_IntCubRule(5*i-2) = CUB_G3_2D
        p_IntFuncComp(1,5*i-2) = i
        p_IntFuncComp(2,5*i-2) = -1
        call fparser_parseFunction(rpostprocessing%pIntChiOmegaParser,5*i-2,&
                    IntRegionOfInterest1,(/'x','y'/))


        ! 4th series of tests
        p_IntCubRule(5*i-1) = CUB_G3_2D
        p_IntFuncComp(1,5*i-1) = -1
        p_IntFuncComp(2,5*i-1) = i
        call fparser_parseFunction(rpostprocessing%pIntChiOmegaParser,5*i-1,&
                    IntRegionOfInterest2,(/'x','y'/))

        ! 5th series of tests
        p_IntCubRule(5*i) = CUB_G3_2D
        p_IntFuncComp(1,5*i) = -1
        p_IntFuncComp(2,5*i) = i
        call fparser_parseFunction(rpostprocessing%pIntChiOmegaParser,5*i,&
                    IntRegionOfInterest1,(/'x','y'/))

    end do




    ! We test the L1-Calculations
    ! How many tests do we do?
    ! first L1-Series: ||f_1(i) - f_2(i)||_L2 on the entire domain
    ! second L1-Series: ||f_1(i)|| on half of the domain
    ! third L1-series: ||f_1(i)|| on the full domain
    ! 4th L1-Series: ||f_2(i)|| on x<0.5
    ! 5th L1-series: ||f_2(i)|| on the whole domain
    ! The first test should always return something close to 0
    ! the second and third one together verify that parsing
    ! still works. So we have 5*NVAR tests
    nL1calculations = 5*rproblem1%NVAR
    rpostprocessing%nL1Calculations = nL1calculations

    ! Init the arrays
    ! Array for components
    call storage_new("ExtFE_setup_tests", &
            "L1CompFunc", (/2,nL1Calculations/),ST_INT, &
             rpostprocessing%h_L1CompFunc,ST_NEWBLOCK_ZERO)

    ! Array for the results
    call storage_new("ExtFE_setup_tests", &
            "L1Results", nL1Calculations,ST_DOUBLE, &
             rpostprocessing%h_L1Results,ST_NEWBLOCK_ZERO)

    ! Array for the cubature rules
    call storage_new("ExtFE_setup_tests", &
            "L1CubRule",nL1Calculations,ST_INT32, &
            rpostprocessing%h_L1CubRule, ST_NEWBLOCK_ZERO)


    ! We will actually parse the L1RegionOfInterest here in a parser
    ! In component i is the region of interest for calculation i
    call fparser_create(rpostprocessing%pL1ChiOmegaParser,nL1Calculations)

    ! Now get the arrays and fill them with data
    call storage_getbase_int2D(rpostprocessing%h_L1CompFunc,p_L1FuncComp)
    call storage_getbase_int32(rpostprocessing%h_L1CubRule,p_L1CubRule)


    L1RegionOfInterest1 = '1'
    L1RegionOfInterest2 = 'IF(x<=0.5,1,0)'
    do i=1,rproblem1%NVAR
        ! First series of test
        ! cub rule
        p_L1CubRule(5*i-4) = CUB_G3_2D
        ! components
        p_L1FuncComp(1,5*i-4) = i
        p_L1FuncComp(2,5*i-4) = i
        ! Parse the region of interest
        call fparser_parseFunction(rpostprocessing%pL1ChiOmegaParser,5*i-4,&
                    L1RegionOfInterest1,(/'x','y'/))

        ! Second series of tests
        p_L1CubRule(5*i -3) = CUB_G3_2D
        p_L1FuncComp(1,5*i -3) = i
        p_L1FuncComp(2,5*i -3) = -1
        call fparser_parseFunction(rpostprocessing%pL1ChiOmegaParser,5*i-3,&
                    L1RegionOfInterest2,(/'x','y'/))

        ! Third series of tests
        p_L1CubRule(5*i-2) = CUB_G3_2D
        p_L1FuncComp(1,5*i-2) = i
        p_L1FuncComp(2,5*i-2) = -1
        call fparser_parseFunction(rpostprocessing%pL1ChiOmegaParser,5*i-2,&
                    L1RegionOfInterest1,(/'x','y'/))


        ! 4th series of tests
        p_L1CubRule(5*i-1) = CUB_G3_2D
        p_L1FuncComp(1,5*i-1) = -1
        p_L1FuncComp(2,5*i-1) = i
        call fparser_parseFunction(rpostprocessing%pL1ChiOmegaParser,5*i-1,&
                    L1RegionOfInterest2,(/'x','y'/))

        ! 5th series of tests
        p_L1CubRule(5*i) = CUB_G3_2D
        p_L1FuncComp(1,5*i) = -1
        p_L1FuncComp(2,5*i) = i
        call fparser_parseFunction(rpostprocessing%pL1ChiOmegaParser,5*i,&
                    L1RegionOfInterest1,(/'x','y'/))

    end do


    ! We test the L2-Calculations
    ! How many tests do we do?
    ! first L2-Series: ||f_1(i) - f_2(i)||_L2 on the entire domain
    ! second L2-Series: ||f_1(i)|| on half of the domain
    ! third L2-series: ||f_1(i)|| on the full domain
    ! 4th Series: ||f_2(i)|| on half of the domain
    ! 5th Series: ||f_2(i)|| on the full domain
    ! The first test should always return something close to 0
    ! the second and third one together verify that parsing
    ! still works. So we have 5*NVAR tests
    nL2calculations = 5*rproblem1%NVAR
    rpostprocessing%nL2Calculations = nL2calculations


    ! Init the arrays
    ! Array for components
    call storage_new("ExtFE_setup_tests", &
            "L2CompFunc", (/2,nL2Calculations/),ST_INT, &
             rpostprocessing%h_L2CompFunc,ST_NEWBLOCK_ZERO)

    ! Array for the results
    call storage_new("ExtFE_setup_tests", &
            "L2Results", nL2Calculations,ST_DOUBLE, &
             rpostprocessing%h_L2Results,ST_NEWBLOCK_ZERO)

    ! Array for the cubature rules
    call storage_new("ExtFE_setup_tests", &
            "L2CubRule",nL2Calculations,ST_INT32, &
            rpostprocessing%h_L2CubRule, ST_NEWBLOCK_ZERO)


    ! We will actually parse the L2RegionOfInterest here in a parser
    ! In component i is the region of interest for calculation i
    call fparser_create(rpostprocessing%pL2ChiOmegaParser,nL2Calculations)

    ! Now get the arrays and fill them with data
    call storage_getbase_int2D(rpostprocessing%h_L2CompFunc,p_L2FuncComp)
    call storage_getbase_int32(rpostprocessing%h_L2CubRule,p_L2CubRule)


    L2RegionOfInterest1 = '1'
    L2RegionOfInterest2 = 'IF(x<=0.5,1,0)'
    do i=1,rproblem1%NVAR
        ! First series of test
        ! cub rule
        p_L2CubRule(5*i-4) = CUB_G3_2D
        ! components
        p_L2FuncComp(1,5*i-4) = i
        p_L2FuncComp(2,5*i-4) = i
        ! Parse the region of interest
        call fparser_parseFunction(rpostprocessing%pL2ChiOmegaParser,5*i-4,&
                    L2RegionOfInterest1,(/'x','y'/))

        ! Second series of tests
        p_L2CubRule(5*i -3) = CUB_G3_2D
        p_L2FuncComp(1,5*i -3) = i
        p_L2FuncComp(2,5*i -3) = -1
        call fparser_parseFunction(rpostprocessing%pL2ChiOmegaParser,5*i-3,&
                    L2RegionOfInterest2,(/'x','y'/))

        ! Third series of tests
        p_L2CubRule(5*i-2) = CUB_G3_2D
        p_L2FuncComp(1,5*i-2) = i
        p_L2FuncComp(2,5*i-2) = -1
        call fparser_parseFunction(rpostprocessing%pL2ChiOmegaParser,5*i-2,&
                    L2RegionOfInterest1,(/'x','y'/))

        ! 4th series of tests
        p_L2CubRule(5*i-1) = CUB_G3_2D
        p_L2FuncComp(1,5*i-1) = -1
        p_L2FuncComp(2,5*i-1) = i
        call fparser_parseFunction(rpostprocessing%pL2ChiOmegaParser,5*i-1,&
                    L2RegionOfInterest2,(/'x','y'/))

        ! 5th series of tests
        p_L2CubRule(5*i) = CUB_G3_2D
        p_L2FuncComp(1,5*i) = -1
        p_L2FuncComp(2,5*i) = i
        call fparser_parseFunction(rpostprocessing%pL2ChiOmegaParser,5*i,&
                    L2RegionOfInterest1,(/'x','y'/))
    end do

    ! Point calculations
    ! How many tests do we want to do?
    ! We want to calculate f-g, (d/dx)f - (d/dx)g, (d/dy)f - (d/dy)g
    ! f, (d/dx)f, (d/dy)f
    nPointCalcs =6*rproblem1%NVAR

    rpostprocessing%nPointCalculations = nPointCalcs

    call storage_new("ExtFEcomparer_benchmark", &
            "PointvalueResults",nPointCalcs,ST_DOUBLE,&
                rpostprocessing%h_PointResults,ST_NEWBLOCK_ZERO)

     ! PointComponent(1,i) = component of function 1 in calculation i
     ! PointComponent(2,i) = derivative of component of function 1 in calculation i
     ! PointComponent(3,i) = component of function 2 in calculation i
     ! PointComponent(4,i) = derivative of component of function 2 in calculation i
     call storage_new("ExtFEcomparer_benchmark", &
                "PointFuncComponents",(/4,nPointCalcs/),ST_INT, &
                 rpostprocessing%h_PointFuncComponents,ST_NEWBLOCK_ZERO)

     ! PointCoordinates(:,i) = coordinates in calculation i
     ! We are in a 2D-Routine!
     call storage_new("ExtFEcomparer_init_postprocessing", &
                "PointCoordinates",(/ExtFE_NDIM2,nPointCalcs/),&
                ST_DOUBLE,rpostprocessing%h_PointCoordinates, &
                ST_NEWBLOCK_ZERO)

    call storage_getbase_int2D(rpostprocessing%h_PointFuncComponents,p_FuncComp)
    call storage_getbase_double2D(rpostprocessing%h_PointCoordinates,p_EvalPoint)

    ! We want to calculate f-g, (d/dx)f-(d/dx)g, (d/dy)f-(d/dy)g,
    !f, (d/dx)f, (d/dy)f
    do i=1,rproblem1%NVAR
            ! First series
            p_EvalPoint(1,6*i-5) = 0.4_DP
            p_EvalPoint(2,6*i-5) = 0.3_DP
            p_FuncComp(1,6*i-5) = i
            p_FuncComp(2,6*i-5) = ExtFE_DER_FUNC_2D
            p_FuncComp(3,6*i-5) = i
            p_FuncComp(4,6*i-5) = ExtFE_DER_FUNC_2D

            ! Second series
            p_EvalPoint(1,6*i-4) = 0.4_DP
            p_EvalPoint(2,6*i-4) = 0.3_DP
            p_FuncComp(1,6*i-4) = i
            p_FuncComp(2,6*i-4) = ExtFE_DER_DERIV_2D_X
            p_FuncComp(3,6*i-4) = i
            p_FuncComp(4,6*i-4) = ExtFE_DER_DERIV_2D_X

            ! Third series
            p_EvalPoint(1,6*i-3) = 0.4_DP
            p_EvalPoint(2,6*i-3) = 0.3_DP
            p_FuncComp(1,6*i-3) = i
            p_FuncComp(2,6*i-3) = ExtFE_DER_DERIV_2D_Y
            p_FuncComp(3,6*i-3) = i
            p_FuncComp(4,6*i-3) = ExtFE_DER_DERIV_2D_Y

            ! Fourth series
            p_EvalPoint(1,6*i-2) = 0.4_DP
            p_EvalPoint(2,6*i-2) = 0.3_DP
            p_FuncComp(1,6*i-2) = i
            p_FuncComp(2,6*i-2) = ExtFE_DER_FUNC_2D
            p_FuncComp(3,6*i-2) = -1
            p_FuncComp(4,6*i-2) = ExtFE_DER_FUNC_2D

            ! Fifth series
            p_EvalPoint(1,6*i-1) = 0.4_DP
            p_EvalPoint(2,6*i-1) = 0.3_DP
            p_FuncComp(1,6*i-1) = i
            p_FuncComp(2,6*i-1) = ExtFE_DER_DERIV_2D_X
            p_FuncComp(3,6*i-1) = -1
            p_FuncComp(4,6*i-1) = ExtFE_DER_FUNC_2D

            ! Sixth series
            p_EvalPoint(1,6*i) = 0.4_DP
            p_EvalPoint(2,6*i) = 0.3_DP
            p_FuncComp(1,6*i) = i
            p_FuncComp(2,6*i) = ExtFE_DER_DERIV_2D_Y
            p_FuncComp(3,6*i) = -1
            p_FuncComp(4,6*i) = ExtFE_DER_FUNC_2D
    end do

end subroutine


! We have some routines for filling the vector, too
! In this routine we fill the vector with ones - so every entry is 1
! no matter to which element/variable it belongs
subroutine ExtFE_fill_vec_ones(rproblem)
    type(t_problem) :: rproblem
    integer :: ifuncComp, ientry
    real(DP), dimension(:), pointer :: p_VecEntries

    do ifuncComp=1,rproblem%coeffVector%nblocks
        ! Get the data array
        call lsyssc_getbase_double(rproblem%coeffVector%RvectorBlock(ifuncComp),p_VecEntries)
        do ientry=1,rproblem%coeffVector%RvectorBlock(ifuncComp)%NEQ
            p_VecEntries(ientry) = 1.0_DP
        end do
    end do ! Loop over all blocks

end subroutine

end module
