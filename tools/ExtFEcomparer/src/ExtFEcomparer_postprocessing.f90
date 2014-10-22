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
  use ExtFEcomparer_typedefs

  implicit none

  private

  public :: ExtFEcomparer_init_postprocessing
  public :: ExtFEcomparer_done_postprocessing
  public :: ExtFEcomparer_postprocess

contains


subroutine ExtFEcomparer_init_postprocessing(rpostprocessing,rparlist)
    !<description> Init of the postprocessing structure:
    ! Save some space for the data
    !</description>

    !<input>
    type(t_parlist), intent(in) :: rparlist
    !</input>

    !<inputoutput>
    type(t_postprocessing) , intent(inout):: rpostprocessing
    !</inputoutput>

    integer :: nL2Calculations, nL2ChiOmega, nL2CubRule
    character(LEN=ExtFE_STRLEN) :: L2filepath
    integer :: writeOutL2

    integer :: nPointEvaluations,nDim
    character(LEN=ExtFE_STRLEN) :: PointFilepath
    integer :: writeOutPoint
    integer, dimension(2) :: iSizeArray


    ! Find out the dimension of the problem
    call parlst_getvalue_int(rparlist,"ExtFE-DOMAININFO","dim",nDim)
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

  ! Easy one: Init arrays for results and which component
  ! and store the filepath
  if (nL2Calculations .gt. 0) then

    call storage_new("ExtFEcomparer_init_postprocessing", &
            "L2CompFunc", (/2,nL2Calculations/),ST_INT, &
             rpostprocessing%h_L2CompFunc,ST_NEWBLOCK_ZERO)
    call storage_new("ExtFEcomparer_init_postprocessing", &
            "L2Results", nL2Calculations,ST_DOUBLE, &
             rpostprocessing%h_L2Results,ST_NEWBLOCK_ZERO)

    call parlst_getvalue_int(rparlist,"ExtFE-POSTPROCESSING", &
                    "writeOutL2Calc",writeOutL2)
    if(writeOutL2 .eq. ExtFE_DO) then

        call storage_new("ExtFEcomparer_init_postprocessing", &
            "L2ChiOmega",ExtFE_STRLEN*nL2Calculations,ST_CHAR, &
            rpostprocessing%h_L2ChiOmega,ST_NEWBLOCK_ZERO)
        call storage_new("ExtFEcomparer_init_postprocessing", &
            "L2ChiOmega",ExtFE_STRLEN*nL2Calculations,ST_CHAR, &
            rpostprocessing%h_L2TriFile,ST_NEWBLOCK_ZERO)
        rpostprocessing%writeOutL2results = .true.
        call parlst_getvalue_string(rparlist,"ExtFE-POSTPROCESSING", &
                        "sL2FilePath",L2filepath,bdequote=.TRUE.)
        if (L2filepath .eq. '' ) then
             call output_line('You want to write out the &
                    &L2-results but have not specified a path', &
                    OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_init_ppostprocessing")
                    call sys_halt()
        end if

        rpostprocessing%L2filepath = L2filepath
    end if

  end if


    !--------------------------------------!
    ! Set up the pointvalue postprocessing !
    !--------------------------------------!

    nPointEvaluations = parlst_querysubstrings(rparlist, &
                    "ExtFE-CALCULATIONS","evaluationPoints")
    if(nPointEvaluations .gt. 0) then
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
        end if
        if (PointFilepath .eq. '' ) then
             call output_line('You want to write out the &
                    &results of the point evaluations but have not specified a path', &
                    OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_init_ppostprocessing")
                    call sys_halt()
        end if

        rpostprocessing%PointFilepath = PointFilepath

    end if


end subroutine

subroutine ExtFEcomparer_done_postprocessing(rpostprocessing)

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


subroutine ExtFEcomparer_postprocess(rpostprocessing)

    !<inputoutput>
    ! the postprocessing structure
    type(t_postprocessing), intent(inout) :: rpostprocessing
    !</inputoutput>

    ! General variables
    integer :: ifile,i,k
    character(LEN=ExtFE_STRLEN) :: soutString, stmpString, stmpString2


    ! For L2-Calculations
    integer :: nL2calc
    integer, dimension(:,:), pointer :: iL2FuncComp => NULL()
    real(DP), dimension(:), pointer :: dL2Results => NULL()
    character, dimension(:), pointer :: cL2ChiOmega => NULL()
    character, dimension(:), pointer :: cL2TriFile => NULL()

    ! For Pointvalues
    integer :: nPointEvals
    integer, dimension(:,:), pointer :: iPointFuncComp => NULL()
    real(DP), dimension(:), pointer :: dPointValues => NULL()
    real(DP), dimension(:,:), pointer :: dPointCoordinates => NULL()

    character(LEN=6), dimension(ExtFE_NumDerivs) :: sderivnames
    ! We start to count the derivatives from 0, so
    ! we need to shift by 1
    sderivnames(ExtFE_DER_FUNC_3D+1) = '      '
    sderivnames(ExtFE_DER_DERIV_3D_X+1) = '(d/dx)'
    sderivnames(ExtFE_DER_DERIV_3D_Y+1) = '(d/dy)'
    sderivnames(ExtFE_DER_DERIV_3D_Z+1) = '(d/dz)'

    ! Let's see what we have to do
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


end module
