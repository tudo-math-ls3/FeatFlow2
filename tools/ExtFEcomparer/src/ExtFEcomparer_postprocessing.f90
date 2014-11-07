! <description>
! This module contains everything related to the
! output of files and postprocessing the results
! from the calculations.
! The postprocessing structure serves in 2 ways:
! 1) We store all results from the calculations in there
! 2) In it, we also set up what we want to calculate.
! The trick behind this? We want only 1 place where we store
! the information what we calculate, and for the output we
! need to know what we did calculate.
! Instead of having to syncronise 2 structures, we just
! use one.
!</description>

module ExtFEcomparer_postprocessing

! Include basic Feat-2 modules
  use fsystem
  use storage
  use genoutput
  use fparser
  use paramlist
  use collection


  use vectorio
  use io
  use ucd
  use cubature
  use fparser
  use linearsystemblock
  use linearsystemscalar

  use ExtFEcomparer_typedefs

  implicit none

  private

  public :: ExtFEcomparer_init_postprocessing
  public :: ExtFEcomparer_done_postprocessing
  public :: ExtFEcomparer_postprocess

contains

! This routine just calls the according initialisations.
subroutine ExtFEcomparer_init_postprocessing(rpostprocessing,rparlist)

    !<input>
    type(t_parlist), intent(inout) :: rparlist
    !</input>

    !<inputoutput>
    type(t_postprocessing) , intent(inout):: rpostprocessing
    !</inputoutput>

    ! General parameters
    call parlst_getvalue_int(rparlist,"ExtFE-DOMAININFO", &
                    "dim",rpostprocessing%nDim)
    call parlst_getvalue_int(rparlist,"ExtFE-FIRST", &
                    "iNVAR",rpostprocessing%nVarFirst)
    call parlst_getvalue_int(rparlist,"ExtFE-SECOND", &
                    "iNVAR",rpostprocessing%nVarSecond)

    ! L2 - Postprocessing
    call ExtFE_init_postprocessing_L2(rpostprocessing,rparlist)

    ! Pointvalue - Postprocessing
    call ExtFE_init_postprocessing_PointValues(rpostprocessing,rparlist)

    ! Ucd-Postprocessing
    call ExtFE_init_postprocessing_UCD(rpostprocessing,rparlist)

    ! Output of original vectors
    call ExtFE_init_postprocessing_OutOrigVec(rpostprocessing,rparlist)


end subroutine

! This routine calls the real postprocessing routines
subroutine ExtFEcomparer_postprocess(rpostprocessing)

    !<inputoutput>
    ! the postprocessing structure
    type(t_postprocessing), intent(inout) :: rpostprocessing
    !</inputoutput>

    call ExtFE_postprocess_L2(rpostprocessing)

    call ExtFE_postprocess_PointValues(rpostprocessing)

    call ExtFE_postprocess_UCD(rpostprocessing)

    call ExtFE_postprocess_OutOrigVec(rpostprocessing)

end subroutine


! This routine calls those routines that clear out everything
subroutine ExtFEcomparer_done_postprocessing(rpostprocessing)

!<input>
   type(t_postprocessing) , intent(inout):: rpostprocessing
!</input>

    call ExtFE_done_postprocessing_L2(rpostprocessing)

    call ExtFE_done_postprocessing_PointValues(rpostprocessing)

    call ExtFE_done_postprocessing_UCD(rpostprocessing)

    call ExtFE_done_postprocessing_OutOrigVec(rpostprocessing)

end subroutine





! Here, the real routines follow

!-----------------------------------------------------!
! Init of L2-Postprocessing

subroutine ExtFE_init_postprocessing_L2(rpostprocessing,rparlist)

    !<input>
    type(t_parlist), intent(in) :: rparlist
    !</input>

    !<inputoutput>
    type(t_postprocessing) , intent(inout):: rpostprocessing
    !</inputoutput>

    integer :: nL2Calculations, nL2ChiOmega, nL2CubRule
    character(LEN=ExtFE_STRLEN) :: L2filepath
    integer :: writeOutL2,i,k,nDIM

    ! pointer to the storage of the postprocessing structure
    real(DP), dimension(:),  pointer :: p_L2results =>NULL()
    integer, dimension (:,:), pointer :: p_L2Comp => NULL()
    character, dimension(:), pointer :: p_L2ChiOmega => NULL()
    character, dimension(:), pointer :: p_L2TriFile => NULL()
    integer(I32), dimension (:), pointer :: p_L2CubRule => NULL()
    character(LEN=ExtFE_STRLEN) :: sparam
    character(LEN=ExtFE_STRLEN) ::sTriFileFirst, sTriFileSecond




    !------------------------------------------!
    ! Set up the L2-Calculation-Postprocessing !
    !------------------------------------------!

    ! Find out how many L2-compuations to do
    nL2Calculations = parlst_querysubstrings(rparlist, &
                    "ExtFE-CALCULATIONS","L2calculations")
    nL2ChiOmega = parlst_querysubstrings(rparlist, &
                    "ExtFE-CALCULATIONS","L2RegionOfInterest")
    nL2CubRule = parlst_querysubstrings(rparlist, &
                    "ExtFE-CALCULATIONS","L2CubatureRule")

  ! Validate that all numbers are the same
  if((nL2Calculations .ne. nL2ChiOmega ).or. &
     (nL2Calculations .ne. nL2CubRule) .or. &
     (nL2ChiOmega .ne. nL2CubRule)) then
     call output_line('Input Error: The number of calculations &
     &must match the number of cubature rules and regions of interest', &
     OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_init_postprocessing")
     call sys_halt()
  end if

  ! Get the dimension
  nDIM = rpostprocessing%nDim

  ! Easy one: Init arrays for results and which component
  ! and store the filepath
  if (nL2Calculations .gt. 0) then

    rpostprocessing%nL2Calculations = nL2Calculations
    call storage_new("ExtFEcomparer_init_postprocessing", &
            "L2CompFunc", (/2,nL2Calculations/),ST_INT, &
             rpostprocessing%h_L2CompFunc,ST_NEWBLOCK_ZERO)
    call storage_new("ExtFEcomparer_init_postprocessing", &
            "L2Results", nL2Calculations,ST_DOUBLE, &
             rpostprocessing%h_L2Results,ST_NEWBLOCK_ZERO)
    call storage_new("ExtFEcomparer_init_postprocessing", &
            "L2CubRule",nL2Calculations,ST_INT32, &
            rpostprocessing%h_L2CubRule, ST_NEWBLOCK_ZERO)
    call storage_new("ExtFEcomparer_init_postprocessing", &
            "L2ChiOmega",ExtFE_STRLEN*nL2Calculations,ST_CHAR, &
            rpostprocessing%h_L2ChiOmega,ST_NEWBLOCK_ZERO)

    call parlst_getvalue_int(rparlist,"ExtFE-POSTPROCESSING", &
                    "writeOutL2Calc",writeOutL2)
    if(writeOutL2 .eq. ExtFE_DO) then

        call storage_new("ExtFEcomparer_init_postprocessing", &
            "L2TriFile",ExtFE_STRLEN*nL2Calculations,ST_CHAR, &
            rpostprocessing%h_L2TriFile,ST_NEWBLOCK_ZERO)
        rpostprocessing%writeOutL2results = .true.
        call parlst_getvalue_string(rparlist,"ExtFE-POSTPROCESSING", &
                        "sL2FilePath",L2filepath,bdequote=.TRUE.)
        if (L2filepath .eq. '' ) then
             call output_line('You want to write out the &
                    &L2-results but have not specified a path', &
                    OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_init_postprocessing")
                    call sys_halt()
        end if

        rpostprocessing%L2filepath = L2filepath
    end if

    ! We will actually parse the L2RegionOfInterest here in a parser
    call fparser_create(rpostprocessing%pL2ChiOmegaParser,nL2Calculations)

    ! Now get the arrays and fill them with data
    call storage_getbase_double(rpostprocessing%h_L2Results,p_L2results)
    call storage_getbase_int2D(rpostprocessing%h_L2CompFunc,p_L2Comp)
    call storage_getbase_char(rpostprocessing%h_L2ChiOmega,p_L2ChiOmega)
    call storage_getbase_char(rpostprocessing%h_L2TriFile,p_L2TriFile)
    call storage_getbase_int(rpostprocessing%h_L2CubRule,p_L2CubRule)

    ! We will store which mesh is used
    call parlst_getvalue_string (rparlist,"ExtFE-FIRST",&
                                 "sMesh",sTriFileFirst,bdequote=.true.)
    call parlst_getvalue_string (rparlist,"ExtFE-SECOND",&
                                 "sMesh",sTriFileFirst,bdequote=.true.)

    ! Now read in the data
    do i=1,nL2Calculations
        ! fetch the whole line
        call parlst_getvalue_string(rparlist, &
                        "ExtFE-CALCULATIONS", "L2calculations", &
                        sparam,sdefault="",isubstring=i)
        ! Now read the parameters from the string
        read(sparam,*) p_L2Comp(1,i), p_L2Comp(2,i)

        call parlst_getvalue_string(rparlist, &
                        "ExtFE-CALCULATIONS", "L2RegionOfInterest", &
                        sparam,sdefault="",isubstring=i)
        ! Copy it to the postprocessing structure for
        ! some output. It is nasty but i did not find
        ! any other way
         do k=1,ExtFE_STRLEN
            p_L2ChiOmega((i-1)*ExtFE_STRLEN+k:(i-1)*ExtFE_STRLEN+k)=sparam(k:k)
         end do

         ! Parse the string in the parser. Unfortunately we have to branch here :(
        select case (nDIM)
             case(ExtFE_NDIM1)
               call fparser_parseFunction(rpostprocessing%pL2ChiOmegaParser,i,sparam,(/'x'/))
             case(ExtFE_NDIM2)
               call fparser_parseFunction(rpostprocessing%pL2ChiOmegaParser,i,sparam,(/'x','y'/))
             case default
               call output_line('Dimension of your problem is &
                     &not supported', &
                     OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_init_postprocessing")
               call sys_halt()
        end select

        call parlst_getvalue_string(rparlist, &
                        "ExtFE-CALCULATIONS", "L2CubatureRule", &
                        sparam,sdefault="",isubstring=i)
        p_L2CubRule(i) = cub_igetID(sparam)

        ! Store which mesh is to be used
        ! ||f-g|| => take first mesh
        if( (p_L2Comp(1,i) .gt. 0) .and. &
            (p_L2comp(2,i) .gt. 0)) then
          do k=1,ExtFE_STRLEN
            p_L2TriFile((i-1)*ExtFE_STRLEN+k:(i-1)*ExtFE_STRLEN+k)=sTriFileFirst(k:k)
          end do

        ! Second option: id1 >0, id2 <=0
        ! => calculate ||id1|| => mesh of first function
        else if((p_L2Comp(1,i) .gt. 0) .and. &
                (p_L2Comp(2,i) .le. 0)) then
          do k=1,ExtFE_STRLEN
            p_L2TriFile((i-1)*ExtFE_STRLEN+k:(i-1)*ExtFE_STRLEN+k)=sTriFileFirst(k:k)
          end do

        ! Third option: id1<=0, id2 >0
        ! => calculate ||id2|| => mesh of second function
        else if((p_L2Comp(1,i) .le. 0) .and. &
                (p_L2Comp(2,i) .gt. 0)) then
          do k=1,ExtFE_STRLEN
            p_L2TriFile((i-1)*ExtFE_STRLEN+k:(i-1)*ExtFE_STRLEN+k)=sTriFileSecond(k:k)
          end do

        else
            call output_line('You need always one component id that &
                    &is greater than 0', &
                    OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_init_postprocessing")
            call sys_halt()
        end if ! Which component


    end do ! Read in of the data

  end if


end subroutine

! Postprocess the L2-Results

subroutine ExtFE_postprocess_L2(rpostprocessing)

    !<inputoutput>
    ! the postprocessing structure
    type(t_postprocessing), intent(inout) :: rpostprocessing
    !</inputoutput>

    ! For L2-Calculations
    integer :: nL2calc
    integer, dimension(:,:), pointer :: iL2FuncComp => NULL()
    real(DP), dimension(:), pointer :: dL2Results => NULL()
    character, dimension(:), pointer :: cL2ChiOmega => NULL()
    character, dimension(:), pointer :: cL2TriFile => NULL()

    ! General variables
    integer :: ifile,i,k
    character(LEN=ExtFE_STRLEN) :: soutString, stmpString, stmpString2

    ! Write out L2-Calculation-Results?
    if (rpostprocessing%writeOutL2results .eqv. .TRUE.) then
        !get pointer do the data arrays
        call storage_getbase_double(rpostprocessing%h_L2Results,dL2Results)
        call storage_getbase_int2D(rpostprocessing%h_L2CompFunc,iL2FuncComp)
        call storage_getbase_char(rpostprocessing%h_L2ChiOmega,cL2ChiOmega)
        call storage_getbase_char(rpostprocessing%h_L2TriFile,cL2TriFile)

        ! How many calculations?
        nL2calc = ubound(dL2Results,1)

        ! Open a file for writing out the results
        call io_openFileForWriting(rpostprocessing%L2filepath, &
                    ifile,SYS_REPLACE,bformatted=.TRUE. )

        ! Ok, file is open, so we write
        do i=1,nL2calc
            ! Empty the string
            soutString = ''
            stmpString = ''
            if((iL2FuncComp(1,i) .gt. 0) .AND. &
                (iL2FuncComp(2,i) .gt. 0)) then
                 write(stmpString,'(A7,I2,A8,I2,A13,I2,A4)') '||f_1(', iL2FuncComp(1,i) , &
                                ') - f_2(', iL2FuncComp(2,i) ,')||_L2(Omega_', i, ') = '
            else if((iL2FuncComp(1,i) .gt. 0 ) .AND. &
                  (iL2FuncComp(2,i) .le. 0) ) then
                 write(stmpString,'(A7,I2,A13,I2,A4)') '||f_1(', iL2FuncComp(1,i) , &
                                ')||_L2(Omega_', i, ') = '
            else
                ! No sanity checks needed - done already during the calculations
                 write(stmpString,'(A7,I2,A13,I2,A4)') '||f_2(', iL2FuncComp(2,i) , &
                                ')||_L2(Omega_', i, ') = '
            end if

            ! Add the result
            write(soutString,'(A1,E16.10)') ' ', dL2Results(i)

            write(ifile,*) trim( trim(stmpString) // trim(soutString) )
        end do

        ! Now we write out what Omega_i actually is
        ! first a linebreak
        write(ifile,*) ''
        do i=1,nL2calc
            ! Empty the strings
            stmpString = ''
            stmpString2 = ''
            write(stmpString,'(A11,I2,A14)') 'with Omega_', i ,' the subset of'

            do k = 1,ExtFE_STRLEN
                stmpString2(k:k) = cL2TriFile((i-1)*ExtFE_STRLEN+k)
            end do !k

            stmpString = trim(stmpString) // ' ' //trim(stmpString2)
            write(stmpString2,'(A24)') ' for that the expression'
            stmpString = trim(stmpString) // trim(stmpString2)

            do k = 1,ExtFE_STRLEN
                stmpString2(k:k) = cL2ChiOmega((i-1)*ExtFE_STRLEN+k)
            end do !k

            write(ifile,'(A,A1,A,A10)') trim(stmpString), ' ', trim(stmpString2), ' returns 1'
        end do ! all L2Calculations


        ! Clean up
        close(ifile)
        dL2Results => NULL()
        iL2FuncComp => NULL()
        cL2ChiOmega => NULL()
    end if


end subroutine

! Release of L2-Postprocessing

subroutine ExtFE_done_postprocessing_L2(rpostprocessing)

!<input>
   type(t_postprocessing) , intent(inout):: rpostprocessing
!</input>

    ! L2-Calculations
    if (rpostprocessing%h_L2CompFunc .gt. ST_NOHANDLE) then
        call storage_free(rpostprocessing%h_L2CompFunc)
    end if
    if (rpostprocessing%h_L2ChiOmega .gt. ST_NOHANDLE) then
        call storage_free(rpostprocessing%h_L2ChiOmega)
    end if
    if (rpostprocessing%h_L2Results .gt. ST_NOHANDLE) then
        call storage_free(rpostprocessing%h_L2Results)
    end if
    if (rpostprocessing%h_L2TriFile .gt. ST_NOHANDLE) then
        call storage_free(rpostprocessing%h_L2TriFile)
    end if
    if(rpostprocessing%h_L2CubRule .gt. ST_NOHANDLE) then
        call storage_free(rpostprocessing%h_L2CubRule)
    end if

    call fparser_release(rpostprocessing%pL2ChiOmegaParser)
end subroutine


! Init of Pointvalue-Postprocessing
subroutine ExtFE_init_postprocessing_PointValues(rpostprocessing,rparlist)
    !<input>
    type(t_parlist), intent(in) :: rparlist
    !</input>

    !<inputoutput>
    type(t_postprocessing) , intent(inout):: rpostprocessing
    !</inputoutput>

    integer :: nPointEvaluations,nDim,i
    character(LEN=ExtFE_STRLEN) :: PointFilepath
    integer :: writeOutPoint
    integer, dimension(2) :: iSizeArray

    real(DP), dimension(:), pointer :: p_PointResults => NULL()
    real(DP), dimension(:,:), pointer :: p_PointCoords => NULL()
    integer, dimension(:,:), pointer :: p_PointFuncComp => NULL()
    character(LEN=ExtFE_STRLEN) :: sparam


    ! Find out the dimension of the problem
    nDim = rpostprocessing%nDim

        nPointEvaluations = parlst_querysubstrings(rparlist, &
                    "ExtFE-CALCULATIONS","evaluationPoints")

    !--------------------------------------!
    ! Set up the pointvalue postprocessing !
    !--------------------------------------!

    if(nPointEvaluations .gt. 0) then
        rpostprocessing%nPointCalculations = nPointEvaluations

        iSizeArray = (/4,nPointEvaluations/)
        ! Create the storage for the information
        ! Result(i) = Result of calculation i
        call storage_new("ExtFEcomparer_init_postprocessing", &
            "PointvalueResults",nPointEvaluations,ST_DOUBLE,&
                rpostprocessing%h_PointResults,ST_NEWBLOCK_ZERO)
        ! PointComponent(1,i) = component of function 1 in calculation i
        ! PointComponent(2,i) = derivative of component of function 1 in calculation i
        ! PointComponent(3,i) = component of function 2 in calculation i
        ! PointComponent(4,i) = derivative of component of function 2 in calculation i

        call storage_new("ExtFEcomparer_init_postprocessing", &
                "PointFuncComponents",iSizeArray,ST_INT, &
                 rpostprocessing%h_PointFuncComponents,ST_NEWBLOCK_ZERO)
        ! PointCoordinates(:,i) = coordinates in calculation i
        call storage_new("ExtFEcomparer_init_postprocessing", &
                "PointCoordinates",(/nDim,nPointEvaluations/),&
                ST_DOUBLE,rpostprocessing%h_PointCoordinates, &
                ST_NEWBLOCK_ZERO)

        call parlst_getvalue_int(rparlist,"ExtFE-POSTPROCESSING", &
                    "writeOutPointValues",writeOutPoint)

        if(writeOutPoint .eq. ExtFE_DO) then
            rpostprocessing%writeOutPointCalucations = .true.
            call parlst_getvalue_string(rparlist,"ExtFE-POSTPROCESSING", &
                        "sPointValuesFilePath",PointFilepath,bdequote=.TRUE.)
            if (PointFilepath .eq. '' ) then
                call output_line('You want to write out the &
                        &results of the point evaluations but have not specified a path', &
                        OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_init_ppostprocessing")
                call sys_halt()
            end if

        rpostprocessing%PointFilepath = PointFilepath

        end if ! writeOUtPoint

        ! Write in the structure what we want to calculate
        ! First step: get the arrays
        call storage_getbase_double2D(rpostprocessing%h_PointCoordinates,p_PointCoords)
        call storage_getbase_double(rpostprocessing%h_PointResults,p_PointResults)
        call storage_getbase_int2D(rpostprocessing%h_PointFuncComponents,p_PointFuncComp)

        do i=1,nPointEvaluations
            ! fetch the whole line
            call parlst_getvalue_string(rparlist, &
                          "ExtFE-CALCULATIONS", "evaluationpoints", &
                            sparam,sdefault="",isubstring=i)
            !and read it in. Which function comp and deriv?
            select case(nDim)

            case(ExtFE_NDIM1)
                read(sparam,*) p_PointFuncComp(1,i), p_PointFuncComp(2,i), &
                            p_PointFuncComp(3,i), p_PointFuncComp(4,i), &
                            p_PointCoords(1,i)
            case(ExtFE_NDIM2)
                read(sparam,*) p_PointFuncComp(1,i), p_PointFuncComp(2,i), &
                            p_PointFuncComp(3,i), p_PointFuncComp(4,i), &
                            p_PointCoords(1,i)  , p_PointCoords(2,i)
            case(ExtFE_NDIM3)
                read(sparam,*) p_PointFuncComp(1,i), p_PointFuncComp(2,i), &
                            p_PointFuncComp(3,i), p_PointFuncComp(4,i), &
                            p_PointCoords(1,i)  , p_PointCoords(2,i), &
                            p_PointCoords(3,i)
            case default
                call output_line('Dimension not implemented', &
                        OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_init_postprocessing")
                call sys_halt()

            end select

        end do ! Read in the Data

    end if

end subroutine


! Postprocess the Pointvalues
subroutine ExtFE_postprocess_PointValues(rpostprocessing)
    !<inputoutput>
    ! the postprocessing structure
    type(t_postprocessing), intent(inout) :: rpostprocessing
    !</inputoutput>

    integer :: nPointEvals
    integer, dimension(:,:), pointer :: iPointFuncComp => NULL()
    real(DP), dimension(:), pointer :: dPointValues => NULL()
    real(DP), dimension(:,:), pointer :: dPointCoordinates => NULL()

    character(LEN=6), dimension(ExtFE_NumDerivs) :: sderivnames

    ! General variables
    integer :: ifile,i,k
    character(LEN=ExtFE_STRLEN) :: soutString, stmpString, stmpString2

    ! We start to count the derivatives from 0, so
    ! we need to shift by 1
    sderivnames(ExtFE_DER_FUNC_3D+1) = '      '
    sderivnames(ExtFE_DER_DERIV_3D_X+1) = '(d/dx)'
    sderivnames(ExtFE_DER_DERIV_3D_Y+1) = '(d/dy)'
    sderivnames(ExtFE_DER_DERIV_3D_Z+1) = '(d/dz)'


    ! Write out Pointvalue-Calculation-Results?
    if(rpostprocessing%writeOutPointCalucations .eqv. .TRUE.) then
        ! Get the arrays
        call storage_getbase_double2D(rpostprocessing%h_PointCoordinates,dPointCoordinates)
        call storage_getbase_double(rpostprocessing%h_PointResults,dPointValues)
        call storage_getbase_int2D(rpostprocessing%h_PointFuncComponents,iPointFuncComp)

        ! find out how many calculations we have
        nPointEvals = ubound(dPointValues,1)

        ! Open a file for writing out the results
        call io_openFileForWriting(rpostprocessing%PointFilepath, &
                    ifile,SYS_REPLACE,bformatted=.TRUE. )

        do i=1,nPointEvals

            ! Empty the string - better save than sorry
            stmpString = ''
            stmpString2 = ''
            soutString = ''

            ! Find out what we have
            ! Both functions
            if((iPointFuncComp(1,i) .gt. 0) .AND. &
               (iPointFuncComp(3,i) .gt. 0)) then
                  write(stmpString,'(A3,I2,A3)') 'f1_', iPointFuncComp(1,i), ' - '
                  write(stmpString2,'(A3,I2)') 'f2_', iPointFuncComp(3,i)
                  stmpString = trim(sderivnames(iPointFuncComp(2,i)+1)) // trim(stmpString)
                  stmpString = '(' // trim(stmpString)
                  stmpString2 = trim(sderivnames(iPointFuncComp(4,i)+1)) // trim(stmpString2)
                  stmpString2 = trim(stmpString2) // ')'
                  stmpString = trim(stmpString) // trim(stmpString2)
            ! Only first
            else if ((iPointFuncComp(1,i) .gt. 0) .AND. &
                    (iPointFuncComp(3,i) .le. 0)) then
                  write(stmpString,'(A3,I2)') 'f1_', iPointFuncComp(1,i)
                  stmpString = (sderivnames(iPointFuncComp(2,i)+1)) // trim(stmpString)
                  stmpString = '(' // trim(stmpString)
                  stmpString = trim(stmpString) // ')'
            ! No sanity check needed - already done during the calculation
            ! it can only be the case: only the second
            else
                  write(stmpString,'(A3,I2)') 'f2_', iPointFuncComp(3,i)
                  stmpString = (sderivnames(iPointFuncComp(4,i)+1)) // trim(stmpString)
                  stmpString = '(' // trim(stmpString)
                  stmpString = trim(stmpString) // ')'
            end if

            ! Now we are still not done, we have to add the function value
            stmpString = trim(stmpString) // '('
            write(stmpString2,'(F8.4)') dPointCoordinates(1,i)
            stmpString = trim(stmpString) // trim(stmpString2)
            do k=2,ubound(dPointCoordinates,1)
                write(stmpString2,'(A1,F8.4)') ',', dPointCoordinates(k,i)
                stmpString = trim(stmpString) // trim(stmpString2)
                stmpString = stmpString // ' '
            end do
            stmpString = trim(stmpString) // ' ) '
            write(stmpString2 ,'(F8.5)') dPointValues(i)
            soutString = trim(stmpString) // ' = ' // trim(stmpString2)

            ! Now we have everything - write it out in a file
            write(ifile,*) trim(soutString)
        end do


        ! Clean up
        close(ifile)
        dPointCoordinates => NULL()
        dPointValues=> NULL()
        iPointFuncComp => NULL()

    end if

end subroutine

! Release of Pointvalue-Postprocessing
subroutine ExtFE_done_postprocessing_PointValues(rpostprocessing)
    !<input>
   type(t_postprocessing) , intent(inout):: rpostprocessing
!</input>

    ! Pointvalue-Calculations
    if (rpostprocessing%h_PointCoordinates .gt. ST_NOHANDLE) then
        call storage_free(rpostprocessing%h_PointCoordinates)
    end if
    if (rpostprocessing%h_PointFuncComponents .gt. ST_NOHANDLE) then
        call storage_free(rpostprocessing%h_PointFuncComponents)
    end if
    if (rpostprocessing%h_PointResults .gt. ST_NOHANDLE) then
        call storage_free(rpostprocessing%h_PointResults)
    end if

end subroutine


! Init of UCD-Postprocessing
subroutine ExtFE_init_postprocessing_UCD(rpostprocessing,rparlist)
    !<input>
    type(t_parlist), intent(inout) :: rparlist
    !</input>

    !<inputoutput>
    type(t_postprocessing) , intent(inout):: rpostprocessing
    !</inputoutput>

    integer :: ucdType,ucdFormat
    integer :: nDim, nVars
    integer :: writeOutMeshes, writeOutFirstFunction, writeOutSecondFunction
    character(LEN=ExtFE_STRLEN) :: sMeshPathFirst, sMeshPathSecond
    character(LEN=ExtFE_STRLEN) :: sFirstFunctionOutput, sSecondFunctionOutput
    character(LEN=ExtFE_STRLEN) :: sUCD_AddElemProjection
    character(LEN=ExtFE_STRLEN) :: sUCD_AddTypeOrig
    character(LEN=ExtFE_STRLEN) :: sUCD_VecComp
    character(LEN=ExtFE_STRLEN) :: sUCD_ScalarComp
    integer :: h_UCD_AddTypeOrig
    integer :: h_UCD_AddElemProject
    integer :: h_UCD_VecVars
    integer :: h_UCD_ScalarVars

    nDim = rpostprocessing%nDim

    !--------------------------------------!
    ! Set up the UCD-Postprocessing        !
    !--------------------------------------!

    ! Write out meshes?
    call parlst_getvalue_int(rparlist,"ExtFE-POSTPROCESSING", &
                    "writeOutMeshesUCD",writeOutMeshes)
    if(writeOutMeshes .eq. ExtFE_DO) then
        call parlst_getvalue_string(rparlist,"ExtFE-POSTPROCESSING", &
                        "sUCDFirstMeshOutput",sMeshPathFirst,bdequote=.TRUE.)
        call parlst_getvalue_string(rparlist,"ExtFE-POSTPROCESSING", &
                        "sUCDSecondMeshOutput",sMeshPathSecond,bdequote=.TRUE.)

        if ( (sMeshPathFirst .eq. '') .or. &
             (sMeshPathSecond .eq. '')  ) then
             call output_line('You want to write out the &
                    &meshes but have not specified a path', &
                    OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_init_postprocessing")
                    call sys_halt()
        end if

        ! Save it in the structure
        rpostprocessing%UCD_meshOneOutPath = sMeshPathFirst
        rpostprocessing%UCD_meshTwoOutPath = sMeshPathSecond
        rpostprocessing%ucd_OUT_meshes = .TRUE.

    end if

    ! Write out first original FE-Function?
    call parlst_getvalue_int(rparlist,"ExtFE-POSTPROCESSING", &
                    "writeOutFirstFunctionUCD",writeOutFirstFunction)
    if(writeOutFirstFunction .eq. ExtFE_DO) then
        call parlst_getvalue_string(rparlist,"ExtFE-POSTPROCESSING", &
                        "sUCDFirstFunctionOutput",sFirstFunctionOutput,bdequote=.TRUE.)

        if (sFirstFunctionOutput .eq. '') then
             call output_line('You want to write out the &
                    &first function but have not specified a path', &
                    OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_init_postprocessing")
                    call sys_halt()
        end if

        ! Save it in the structure
        rpostprocessing%UCD_FEfunctionOneOrigOutPath = sFirstFunctionOutput
        rpostprocessing%ucd_OUT_orig_functions_one = .TRUE.

        sUCD_AddElemProjection = "sUCDElemProjectFirstFunction"
        sUCD_AddTypeOrig = "sUCDaddTypeFirstFunction"
        sUCD_VecComp = "sUCDVecCompFirstFunction"
        sUCD_ScalarComp = "UCD_ScalarFirstOrig"
        nVars = rpostprocessing%nVarFirst

        call ExtFE_init_pp_UCD_fefuncout(h_UCD_AddElemProject,h_UCD_AddTypeOrig,&
                                    h_UCD_VecVars,h_UCD_ScalarVars,&
                            sUCD_AddElemProjection,sUCD_AddTypeOrig,&
                                        sUCD_VecComp,sUCD_ScalarComp,&
                                        nDim,&
                                        nVars,&
                                        rparlist)

        rpostprocessing%h_UCD_AddTypeOrigFirst = h_UCD_AddTypeOrig
        rpostprocessing%h_UCD_AddElemProjectFirst = h_UCD_AddElemProject
        rpostprocessing%h_UCD_VecsFirstOrig = h_UCD_VecVars
        rpostprocessing%h_UCD_ScalarFirstOrig = h_UCD_ScalarVars


    end if

    ! Write out second original FE-Function?
    call parlst_getvalue_int(rparlist,"ExtFE-POSTPROCESSING", &
                    "writeOutSecondFunctionUCD",writeOutSecondFunction)
    if(writeOutSecondFunction .eq. ExtFE_DO) then
        call parlst_getvalue_string(rparlist,"ExtFE-POSTPROCESSING", &
                        "sUCDSecondFunctionOutput",sSecondFunctionOutput,bdequote=.TRUE.)

        if (sSecondFunctionOutput .eq. '') then
             call output_line('You want to write out the &
                    &second function but have not specified a path', &
                    OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_init_postprocessing")
                    call sys_halt()
        end if

        ! Save it in the structure
        rpostprocessing%UCD_FEfunctionTwoOrigOutPath = sSecondFunctionOutput
        rpostprocessing%ucd_OUT_orig_functions_two = .TRUE.
        sUCD_AddElemProjection = "sUCDElemProjectSecondFunction"
        sUCD_AddTypeOrig = "sUCDaddTypeSecondFunction"
        sUCD_VecComp = "sUCDVecCompSecondFunction"
        sUCD_ScalarComp = "UCD_ScalarSecondOrig"
        nVars = rpostprocessing%nVarSecond

        call ExtFE_init_pp_UCD_fefuncout(h_UCD_AddElemProject,h_UCD_AddTypeOrig,&
                                    h_UCD_VecVars,h_UCD_ScalarVars,&
                            sUCD_AddElemProjection,sUCD_AddTypeOrig,&
                                        sUCD_VecComp,sUCD_ScalarComp,&
                                        nDim,&
                                        nVars,&
                                        rparlist)

        rpostprocessing%h_UCD_AddTypeOrigSecond = h_UCD_AddTypeOrig
        rpostprocessing%h_UCD_AddElemProjectSecond = h_UCD_AddElemProject
        rpostprocessing%h_UCD_VecsSecondOrig = h_UCD_VecVars
        rpostprocessing%h_UCD_ScalarSecondOrig = h_UCD_ScalarVars


    end if



    if((rpostprocessing%ucd_OUT_meshes .eqv. .TRUE.) .or. &
        (rpostprocessing%ucd_OUT_orig_functions_one .eqv. .TRUE.) .or. &
        (rpostprocessing%ucd_OUT_orig_functions_two .eqv. .TRUE.))then
        call parlst_getvalue_int(rparlist,"ExtFE-POSTPROCESSING", &
                        "ucdFormat",ucdFormat)
        rpostprocessing%UCD_Format = ucdFormat

        call parlst_getvalue_int(rparlist,"ExtFE-POSTPROCESSING", &
                        "ucdType",ucdType)
        select case(ucdType)
            case(ExtFE_UCD_OUT_STANDARD)
                rpostprocessing%UCD_Style = UCD_FLAG_STANDARD
            case default
                call output_line('Wrong Input for the UCD-type',&
                   OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_init_postprocessing")
                call sys_halt()
            end select
    end if


contains

 !This one serves as a wrapper as we do 2 times the same
 subroutine ExtFE_init_pp_UCD_fefuncout(h_UCD_AddElemProject,h_UCD_AddTypeOrig,&
                                    h_UCD_VecVars,h_UCD_ScalarVars,&
                            sUCD_AddElemProjection,sUCD_AddTypeOrig,&
                                        sUCD_VecComp,sUCD_ScalarComp,&
                                        nDimension, nVarsFunc,rparlist)

    ! Input
    type(t_parlist), intent(inout) :: rparlist
    ! Strings for where to search
    character(LEN=ExtFE_STRLEN), intent(in) :: sUCD_AddTypeOrig
    character(LEN=ExtFE_STRLEN), intent(in) :: sUCD_AddElemProjection
    character(LEN=ExtFE_STRLEN), intent(in) :: sUCD_VecComp
    character(LEN=ExtFE_STRLEN), intent(in) :: sUCD_ScalarComp
    integer, intent(in) :: nDimension
    integer, intent(in) :: nVarsFunc

    ! InOut: Handles to the storage of the arrays
    integer, intent(inout) :: h_UCD_AddTypeOrig
    integer, intent(inout) :: h_UCD_AddElemProject
    integer, intent(inout) :: h_UCD_VecVars
    integer, intent(inout) :: h_UCD_ScalarVars

    ! Local variables
    integer :: i,k, iVarsScalar, nVecs
    integer, dimension (:,:), pointer :: p_IntPointerVecs => NULL()
    integer, dimension(:), pointer :: p_IntPointerScalars => NULL()
    integer, dimension(:), pointer :: p_IntPointerAddElemProject => NULL()
    integer, dimension(:), pointer :: p_IntPointerAddType => NULL()
    integer, dimension(:), allocatable :: VarsNotVector

    character(LEN=ExtFE_STRLEN) :: sparam

    ! Allocate the arrays for the UCD-Output
    ! Therefore we need the number of variables in the function
    ! and how many vectors are in there

    call storage_new("ExtFE_init_postprocessing",sUCD_AddTypeOrig,&
                      nVarsFunc,ST_INT,h_UCD_AddTypeOrig,&
                      ST_NEWBLOCK_ZERO)
    call storage_new("ExtFE_init_postprocessing",sUCD_AddElemProjection,&
                      nVarsFunc,ST_INT,h_UCD_AddElemProject,&
                      ST_NEWBLOCK_ZERO)

    nVecs = parlst_querysubstrings(rparlist, &
                  "ExtFE-POSTPROCESSING",sUCD_VecComp)
    if(nVecs .gt. 0) then
        call storage_new("ExtFE_init_postprocessing",sUCD_VecComp,&
                        (/nVecs,nDimension/),ST_INT,h_UCD_VecVars,&
                        ST_NEWBLOCK_ZERO)

        ! Fetch the whole line with the flags and then translate them
        ! Type of variable
        call storage_getbase_int2D(h_UCD_VecVars,p_IntPointerVecs)
        do i=1,nVecs
            ! fetch the whole line
            call parlst_getvalue_string(rparlist, &
                          "ExtFE-Postprocessing", sUCD_VecComp, &
                            sparam,sdefault="",isubstring=i)
            read(sparam,*) p_IntPointerVecs(i,:)

        end do ! Read in the Data
        ! Done with that one

        ! Assumtion: Every variable not in a vector shall be written out as scalar variable
        ! We know we have maximum nVar variables
        allocate(VarsNotVector(nVarsFunc))
        VarsNotVector(:) = 1
        do i=1,nVecs
            do k=1,nDimension
                VarsNotVector(p_IntPointerVecs(i,k)) = 0
            end do
        end do
        ! Now all remaining 1 tell us which variable are scalar ones, count them
        iVarsScalar = 0
        do i=1,nVarsFunc
            iVarsScalar = iVarsScalar + VarsNotVector(i)
        end do


        if(iVarsScalar .gt. 0) then
            call storage_new("ExtFE_init_postprocessing",sUCD_ScalarComp,&
                   iVarsScalar,ST_INT,h_UCD_ScalarVars,&
                    ST_NEWBLOCK_ZERO)
            call storage_getbase_int(h_UCD_ScalarVars,p_IntPointerScalars)
            k = 0
            do i=1,nVarsFunc
                if (VarsNotVector(i) .eq. 1) then
                    k = k+1
                    p_IntPointerScalars(k) = i
                end if
            end do
        end if
        deallocate(VarsNotVector)
        ! Another case is not possible: if nVecs = 0 and also no variable is a scalar one there is no
        ! component in the function
    else
        iVarsScalar = nVarsFunc
        call storage_new("ExtFE_init_postprocessing",sUCD_ScalarComp,&
                  iVarsScalar,ST_INT,h_UCD_ScalarVars,&
                  ST_NEWBLOCK_ZERO)
        call storage_getbase_int(h_UCD_ScalarVars,p_IntPointerScalars)
        do i=1,nVarsFunc
            p_IntPointerScalars(i) = i
        end do

    end if ! nVecs > 0

    ! Now how to add them
    call storage_getbase_int(h_UCD_AddTypeOrig,p_IntPointerAddType)
    call parlst_getvalue_string(rparlist,"ExtFE-Postprocessing",sUCD_AddTypeOrig,sparam)

    read(sparam,*) p_IntPointerAddType(:)
    ! Done with that
    ! Sanity check
    do i=1,nVarsFunc
        select case(p_IntPointerAddType(i))

        case(ExtFE_UCD_ELEM_BASED,ExtFE_UCD_VERT_BASED)
!            write(sparam,*) p_IntPointerAddType(i)
!            call output_line('Function one: add Type= ' //trim(adjustl(sparam)))
        case default
             write(sparam,*) i
             call output_line('wrong add type for a function, variable '//trim(adjustl(sparam)), &
                      OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_init_postprocessing")
             call sys_halt()
         end select
    end do

    ! Poly degree?
    call storage_getbase_int(h_UCD_AddElemProject,p_IntPointerAddElemProject)
    call parlst_getvalue_string(rparlist,"ExtFE-Postprocessing",sUCD_AddElemProjection,sparam)

    read(sparam,*) p_IntPointerAddElemProject(:)
    ! Done with that
    ! Sanity check
    do i=1,nVarsFunc
         select case(p_IntPointerAddElemProject(i))

         case(ExtFE_UCD_POLY_CONST,ExtFE_UCD_POLY_LINEAR)
           case default
             write(sparam,*) i
             call output_line('not supported poly degree for variable '//trim(adjustl(sparam)), &
                     OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_init_postprocessing")
              call sys_halt()
         end select
    end do

    ! Clean up
    p_IntPointerAddElemProject => NULL()
    p_IntPointerAddType => NULL()
    p_IntPointerScalars => NULL()
    p_IntPointerVecs => NULL()
 end subroutine ! ExtFE_init_pp_UCD_fefuncout

end subroutine ! ExtFE_init_UCD

! Make the desired UCD-Output
subroutine ExtFE_postprocess_UCD(rpostprocessing)

    !<inputoutput>
    ! the postprocessing structure
    type(t_postprocessing), intent(inout) :: rpostprocessing
    !</inputoutput>

    type(t_ucdExport) :: rexport_first_mesh,rexport_second_mesh

    ! Write out the meshes with UCD-Output?
    if(rpostprocessing%ucd_OUT_meshes .eqv. .TRUE.) then
        select case(rpostprocessing%UCD_Format)
            case(ExtFE_UCD_VTK)
                ! Open / write / close
                call ucd_startVTK (rexport_first_mesh,UCD_FLAG_STANDARD,&
                        rpostprocessing%UCD_MeshOnePointer,&
                        rpostprocessing%UCD_meshOneOutPath)
                call ucd_startVTK (rexport_second_mesh,UCD_FLAG_STANDARD,&
                        rpostprocessing%UCD_MeshTwoPointer,&
                        rpostprocessing%UCD_meshTwoOutPath)

            case default
                call output_line('Type of UCD-output &
                     &not supported', &
                     OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_postprocessing")
                call sys_halt()
        end select

        call ucd_write (rexport_first_mesh)
        call ucd_write (rexport_second_mesh)
        call ucd_release (rexport_first_mesh)
        call ucd_release (rexport_second_mesh)

    end if

    ! Write out the first function?
    if(rpostprocessing%ucd_OUT_orig_functions_one .eqv. .TRUE. ) then

        call ExtFE_write_UCD(rpostprocessing%UCD_feFunction_first_orig,&
                            rpostprocessing%UCD_Format,&
                            rpostprocessing%UCD_Style,&
                            rpostprocessing%UCD_FEfunctionOneOrigOutPath,&
                            rpostprocessing%h_UCD_AddTypeOrigFirst,&
                            rpostprocessing%h_UCD_ScalarFirstOrig,&
                            rpostprocessing%h_UCD_VecsFirstOrig)
    end if

    if(rpostprocessing%ucd_OUT_orig_functions_two .eqv. .TRUE. ) then

        call ExtFE_write_UCD(rpostprocessing%UCD_feFunction_second_orig,&
                            rpostprocessing%UCD_Format,&
                            rpostprocessing%UCD_Style,&
                            rpostprocessing%UCD_FEfunctionTwoOrigOutPath,&
                            rpostprocessing%h_UCD_AddTypeOrigSecond,&
                            rpostprocessing%h_UCD_ScalarSecondOrig,&
                            rpostprocessing%h_UCD_VecsSecondOrig)
    end if


    ! We create a wrapper here since we basically do every time the same,
    ! just some different input, ie another function/other flags

  contains

  subroutine ExtFE_write_UCD(fefunction,ucd_fmt,ucd_initType,outpath,h_addType,h_scalarValues,h_VectorValues)
  ! Input
  type(t_vectorBlock), intent(inout) :: fefunction
  integer, intent(in) :: ucd_fmt
  integer(I32), intent(in) :: ucd_initType
  ! outpath
  character(LEN=ExtFE_STRLEN), intent(in) :: outpath
  ! add them vertex-based, midpoint-based,...
  integer, intent(in) :: h_addType
  ! scalar values
  integer, intent(in) :: h_scalarValues
  integer, intent(in) :: h_VectorValues

  ! Local variables
  integer :: i, nVecs, nScalars, nDim
  type(t_ucdExport) :: rexportFunction
  character(LEN=ExtFE_STRLEN) :: svarname
  integer, dimension(:), pointer :: p_addType => NULL()
  integer, dimension(:), pointer :: p_scalarValues => NULL()
  integer, dimension(:,:), pointer :: p_VectorValues => NULL()
  real(DP), dimension(:), pointer :: p_DataX => NULL()
  real(DP), dimension(:), pointer :: p_DataY => NULL()
  real(DP), dimension(:), pointer :: p_DataZ => NULL()


 ! Select ucd-type
 select case(ucd_fmt)

   case(ExtFE_UCD_VTK)
   ! Open
     call ucd_startVTK (rexportFunction,ucd_initType,&
                    fefunction%p_rblockDiscr%p_rtriangulation,&
                    outpath)

   case default

     call output_line('Type of UCD-output &
                    &not supported', &
           OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_postprocessing")
     call sys_halt()
   end select

   ! How to add?
   call storage_getbase_int(h_addType,p_addType)

   ! First we go through the scalar ones
   if (h_scalarValues .gt. ST_NOHANDLE) then
     call storage_getbase_int(h_scalarValues,p_scalarValues)
     nScalars = ubound(p_scalarValues,1)
     do i=1,nScalars
       write(svarname,*) i
        call lsyssc_getbase_double(fefunction%RvectorBlock(p_scalarValues(i)),p_DataX)
         select case(p_addType(p_scalarValues(i)))

          case(ExtFE_UCD_VERT_BASED)
             !call ucd_addVarVertBasedVec(rexportFunction,"SclarVar_"//trim(adjustl(svarname)),p_DataX)
             call ucd_addVariableVertexBased(rexportFunction,"ScalarVar_"//trim(adjustl(svarname)),p_DataX)
          case(ExtFE_UCD_ELEM_BASED)
             !call ucd_addVarElemBasedVec(rexportFunction,"ScalarVar_"//trim(adjustl(svarname)),p_DataX)
             call ucd_addVariableElementBased(rexportFunction,"ScalarVar_"//trim(adjustl(svarname)),p_DataX)
        end select

     end do
   end if

   ! Now the vector valued
   ! Assumtion: All Variables in the vector are added the same way!
   ! first on vertex-based => All vertex based etc
   ! We don't need sanity checks as we did them during the init
   if(h_VectorValues .gt. ST_NOHANDLE) then
        call storage_getbase_INT2D(h_VectorValues,p_VectorValues)
        nVecs = ubound(p_VectorValues,1)
        nDim = fefunction%p_rblockDiscr%ndimension
        svarname = ''
        do i=1,nVecs
            write(svarname,*) i
            svarname = "VectorVar_" //trim(adjustl(svarname))
            select case(nDim)
            case(ExtFE_NDIM1)

                call lsyssc_getbase_double(fefunction%RvectorBlock(p_VectorValues(i,1)),p_DataX)

                select case(p_addType(p_VectorValues(i,1)))
                case(ExtFE_UCD_VERT_BASED)
                    call ucd_addVarVertBasedVec(rexportFunction,svarname,p_DataX)
                case (ExtFE_UCD_ELEM_BASED)
                    call ucd_addVarElemBasedVec(rexportFunction,svarname,p_DataX)
                end select

            case(ExtFE_NDIM2)
                call lsyssc_getbase_double(fefunction%RvectorBlock(p_VectorValues(i,1)),p_DataX)
                call lsyssc_getbase_double(fefunction%RvectorBlock(p_VectorValues(i,2)),p_DataY)

                select case(p_addType(p_VectorValues(i,1)))
                case(ExtFE_UCD_VERT_BASED)
                    call ucd_addVarVertBasedVec(rexportFunction,svarname,p_DataX,p_DataY)
                case (ExtFE_UCD_ELEM_BASED)
                    call ucd_addVarElemBasedVec(rexportFunction,svarname,p_DataX,p_DataY)
                end select

            case(ExtFE_NDIM3)
                call lsyssc_getbase_double(fefunction%RvectorBlock(p_VectorValues(i,1)),p_DataX)
                call lsyssc_getbase_double(fefunction%RvectorBlock(p_VectorValues(i,2)),p_DataY)
                call lsyssc_getbase_double(fefunction%RvectorBlock(p_VectorValues(i,3)),p_DataZ)

                select case(p_addType(p_VectorValues(i,1)))
                case(ExtFE_UCD_VERT_BASED)
                    call ucd_addVarVertBasedVec(rexportFunction,svarname,p_DataX,p_DataY,p_DataZ)
                case (ExtFE_UCD_ELEM_BASED)
                    call ucd_addVarElemBasedVec(rexportFunction,svarname,p_DataX,p_DataY,p_DataZ)
                end select

            end select

      end do


   end if


  ! Write and release
  call ucd_write(rexportFunction)
  call ucd_release(rexportFunction)

 end subroutine !ExtFE_write_UCD

end subroutine !postprocess_UCD

! Release of UCD-Postprocessing
subroutine ExtFE_done_postprocessing_UCD(rpostprocessing)
!<input>
   type(t_postprocessing) , intent(inout):: rpostprocessing
!</input>
    ! Meshes
    if(rpostprocessing%ucd_OUT_meshes .eqv. .true. ) then
        rpostprocessing%UCD_MeshOnePointer => NULL()
        rpostprocessing%UCD_MeshTwoPointer => NULL()
    end if

    ! UCD Out of first function
    if(rpostprocessing%ucd_OUT_orig_functions_one .eqv. .true.) then
        call lsysbl_releaseVector(rpostprocessing%UCD_feFunction_first_orig)
        deallocate(rpostprocessing%UCD_feFunction_first_orig)
        call spdiscr_releaseBlockDiscr(rpostprocessing%UCDBlockDiscrFirst)
        deallocate(rpostprocessing%UCDBlockDiscrFirst)
    end if
    if(rpostprocessing%h_UCD_AddElemProjectFirst .gt. ST_NOHANDLE) then
        call storage_free(rpostprocessing%h_UCD_AddElemProjectFirst)
    end if

    if(rpostprocessing%h_UCD_VecsFirstOrig .gt. ST_NOHANDLE) then
        call storage_free(rpostprocessing%h_UCD_VecsFirstOrig)
    end if

    if(rpostprocessing%h_UCD_AddTypeOrigFirst .gt. ST_NOHANDLE) then
        call storage_free(rpostprocessing%h_UCD_AddTypeOrigFirst)
    end if
    if(rpostprocessing%h_UCD_ScalarFirstOrig .gt. ST_NOHANDLE) then
        call storage_free(rpostprocessing%h_UCD_ScalarFirstOrig)
    end if

    ! UCD-Out of second function
    if(rpostprocessing%ucd_OUT_orig_functions_two .eqv. .true.) then
        call lsysbl_releaseVector(rpostprocessing%UCD_feFunction_second_orig)
        deallocate(rpostprocessing%UCD_feFunction_second_orig)
        call spdiscr_releaseBlockDiscr(rpostprocessing%UCDBlockDiscrSecond)
        deallocate(rpostprocessing%UCDBlockDiscrSecond)
    end if
    if(rpostprocessing%h_UCD_AddElemProjectSecond .gt. ST_NOHANDLE) then
        call storage_free(rpostprocessing%h_UCD_AddElemProjectSecond)
    end if

    if(rpostprocessing%h_UCD_VecsSecondOrig .gt. ST_NOHANDLE) then
        call storage_free(rpostprocessing%h_UCD_VecsSecondOrig)
    end if

    if(rpostprocessing%h_UCD_AddTypeOrigSecond .gt. ST_NOHANDLE) then
        call storage_free(rpostprocessing%h_UCD_AddTypeOrigSecond)
    end if
    if(rpostprocessing%h_UCD_ScalarSecondOrig .gt. ST_NOHANDLE) then
        call storage_free(rpostprocessing%h_UCD_ScalarSecondOrig)
    end if



end subroutine


! Init the output of the vectors
subroutine ExtFE_init_postprocessing_OutOrigVec(rpostprocessing,rparlist)

    !<input>
    type(t_parlist), intent(inout) :: rparlist
    !</input>

    !<inputoutput>
    type(t_postprocessing) , intent(inout):: rpostprocessing
    !</inputoutput>

    integer :: writeOutOrigVector1, writeOutOrigVector2
    character(LEN=ExtFE_STRLEN) :: sPathFirst, sPathSecond, sFMT1, sFMT2

    call parlst_getvalue_int(rparlist,"ExtFE-POSTPROCESSING", &
                    "writeOutOrigVector1",writeOutOrigVector1)

    if(writeOutOrigVector1 .eq. ExtFE_DO) then
        call parlst_getvalue_string(rparlist,"ExtFE-POSTPROCESSING", &
                        "sOutPathFirstOrigVec",sPathFirst,bdequote=.TRUE.)
        if (sPathFirst .eq. '') then
             call output_line('You want to convert the &
                    &first vector but have not specified a path to save it', &
                    OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_init_postprocessing")
                    call sys_halt()
        end if

        call parlst_getvalue_string(rparlist,"ExtFE-POSTPROCESSING", &
                        "sOutVec1OrigFMT",sFMT1,bdequote=.TRUE.)
        if(sFMT1 .eq. '') then
            sFMT1 = '(E22.15)'
        end if

        ! Save it in the postprocessing structure
        rpostprocessing%sOrigVec1OutFMT = sFMT1
        rpostprocessing%writeOutOrigVector1 = .TRUE.
        rpostprocessing%sOrigVecPathOutFirst = sPathFirst
    end if


    call parlst_getvalue_int(rparlist,"ExtFE-POSTPROCESSING", &
                    "writeOutOrigVector2",writeOutOrigVector2)

    if(writeOutOrigVector2 .eq. ExtFE_DO) then
        call parlst_getvalue_string(rparlist,"ExtFE-POSTPROCESSING", &
                        "sOutPathSecondOrigVec", sPathSecond ,bdequote=.TRUE.)

        if (sPathSecond .eq. '') then
             call output_line('You want to convert the &
                    &second vector but have not specified a path to save it', &
                    OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_init_postprocessing")
                    call sys_halt()
        end if

        call parlst_getvalue_string(rparlist,"ExtFE-POSTPROCESSING", &
                        "sOutVec2OrigFMT",sFMT2,bdequote=.TRUE.)
        if(sFMT2 .eq. '') then
            sFMT2 = '(E22.15)'
        end if

        ! Save it in the postprocessing structure
        rpostprocessing%sOrigVec2OutFMT = sFMT2
        rpostprocessing%writeOutOrigVector2 = .TRUE.
        rpostprocessing%sOrigVecPathOutSecond = sPathSecond
    end if

end subroutine

! Write out the vectors
subroutine ExtFE_postprocess_OutOrigVec(rpostprocessing)

    !<inputoutput>
    ! the postprocessing structure
    type(t_postprocessing), intent(inout) :: rpostprocessing
    !</inputoutput>

    logical :: bunsorted = .FALSE.
    character(LEN=ExtFE_STRLEN) :: comment
    comment = 'Converted with the ExtFEcomparer'

    if(rpostprocessing%writeOutOrigVector1 .eqv. .TRUE.) then
        call vecio_writeBlockVectorHR(rpostprocessing%OrigVecFirst,'SOLUTION', &
                bunsorted,0,rpostprocessing%sOrigVecPathOutFirst,&
                sformat=trim(rpostprocessing%sOrigVec1OutFMT),&
                scomment=comment)
    end if

    if(rpostprocessing%writeOutOrigVector2 .eqv. .TRUE.) then
        call vecio_writeBlockVectorHR(rpostprocessing%OrigVecSecond,'SOLUTION', &
                bunsorted,0,rpostprocessing%sOrigVecPathOutSecond,&
                sformat=trim(rpostprocessing%sOrigVec2OutFMT),&
                scomment=comment)
    end if


end subroutine


subroutine ExtFE_done_postprocessing_OutOrigVec(rpostprocessing)

!<input>
   type(t_postprocessing) , intent(inout):: rpostprocessing
!</input>

    rpostprocessing%OrigVecFirst => NULL()
    rpostprocessing%OrigVecSecond => NULL()

end subroutine

end module
