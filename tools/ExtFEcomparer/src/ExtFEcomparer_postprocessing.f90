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
! In this program, we only know 2 types of operations
! with a FE-Function: calculations and postprocessing.
! Everything regarding file i/o is considered to be
! postprocessing.
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

    ! Integral - Postprocessing
    call ExtFE_init_postprocessing_Integral(rpostprocessing,rparlist)

    ! L1 - Postprocessing
    call ExtFE_init_postprocessing_L1(rpostprocessing,rparlist)

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

    call ExtFE_postprocess_Integral(rpostprocessing)

    call ExtFE_postprocess_L1(rpostprocessing)

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

    call ExtFE_done_postprocessing_Integral(rpostprocessing)

    call ExtFE_done_postprocessing_L1(rpostprocessing)

    call ExtFE_done_postprocessing_L2(rpostprocessing)

    call ExtFE_done_postprocessing_PointValues(rpostprocessing)

    call ExtFE_done_postprocessing_UCD(rpostprocessing)

    call ExtFE_done_postprocessing_OutOrigVec(rpostprocessing)

end subroutine





! Here, the real routines follow


!-----------------------------------------------------!
! Init of the Integral-Postprocessing

subroutine ExtFE_init_postprocessing_Integral(rpostprocessing,rparlist)

    !<input>
    type(t_parlist), intent(in) :: rparlist
    !</input>

    !<inputoutput>
    type(t_postprocessing) , intent(inout):: rpostprocessing
    !</inputoutput>

    integer :: nIntCalculations, nIntChiOmega, nIntCubRule
    character(LEN=ExtFE_STRLEN) :: Intfilepath
    integer :: writeOutInt,i,k,nDIM

    ! pointer to the storage of the postprocessing structure
    real(DP), dimension(:),  pointer :: p_Intresults =>NULL()
    integer, dimension (:,:), pointer :: p_IntComp => NULL()
    character, dimension(:), pointer :: p_IntChiOmega => NULL()
    character, dimension(:), pointer :: p_IntTriFile => NULL()
    integer(I32), dimension (:), pointer :: p_IntCubRule => NULL()
    character(LEN=ExtFE_STRLEN) :: sparam
    character(LEN=ExtFE_STRLEN) ::sTriFileFirst, sTriFileSecond
    integer :: iaux




    !------------------------------------------------!
    ! Set up the Integral-Calculation-Postprocessing !
    !------------------------------------------------!

    ! Find out how many Integral-compuations to do
    nIntCalculations = parlst_querysubstrings(rparlist, &
                    "ExtFE-CALCULATIONS","IntCalculations")
    ! How many Regions of interest are there?
    nIntChiOmega = parlst_querysubstrings(rparlist, &
                    "ExtFE-CALCULATIONS","IntRegionOfInterest")
    ! How many cubature rules are there?
    nIntCubRule = parlst_querysubstrings(rparlist, &
                    "ExtFE-CALCULATIONS","IntCubatureRule")

  ! Validate that all numbers are the same
  if((nIntCalculations .ne. nIntChiOmega ).or. &
     (nIntCalculations .ne. nIntCubRule) .or. &
     (nIntChiOmega .ne. nIntCubRule)) then
     call output_line('Input Error: The number of calculations &
     &must match the number of cubature rules and regions of interest', &
     OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_init_postprocessing")
     call sys_halt()
  end if

  ! Get the dimension
  nDIM = rpostprocessing%nDim

  ! Init arrays for results and which component is used
  ! and store the filepath
  if (nIntCalculations .gt. 0) then

    ! Store the number in the structure
    rpostprocessing%nIntCalculations = nIntCalculations

    ! Array for components
    call storage_new("ExtFEcomparer_init_postprocessing", &
            "IntCompFunc", (/2,nIntCalculations/),ST_INT, &
             rpostprocessing%h_IntCompFunc,ST_NEWBLOCK_ZERO)

    ! Array for the results
    call storage_new("ExtFEcomparer_init_postprocessing", &
            "IntResults", nIntCalculations,ST_DOUBLE, &
             rpostprocessing%h_IntResults,ST_NEWBLOCK_ZERO)

    ! Array for the cubature rules
    iaux = ST_INT32
    call storage_new("ExtFEcomparer_init_postprocessing", &
            "IntCubRule",nIntCalculations, iaux, &
            rpostprocessing%h_IntCubRule, ST_NEWBLOCK_ZERO)

    ! Array to store the region of interest for some
    ! postprocessing
    ! Only of interest if we write out the results in a file,
    ! but else we have some branches during the init and I
    ! don't want to do them right now.
    call storage_new("ExtFEcomparer_init_postprocessing", &
            "IntChiOmega",ExtFE_STRLEN*nIntCalculations,ST_CHAR, &
            rpostprocessing%h_IntChiOmega,ST_NEWBLOCK_ZERO)

    ! We want to write in the file which mesh was used
    ! Get an array for that
    call storage_new("ExtFEcomparer_init_postprocessing", &
          "IntTriFile",ExtFE_STRLEN*nIntCalculations,ST_CHAR, &
           rpostprocessing%h_IntTriFile,ST_NEWBLOCK_ZERO)

    ! Do we want to save the results in a file?
    call parlst_getvalue_int(rparlist,"ExtFE-POSTPROCESSING", &
                    "writeOutIntCalc",writeOutInt)

    if(writeOutInt .eq. ExtFE_DO) then

        ! Tell the postprocessing that we want to write a file
        rpostprocessing%writeOutIntresults = .true.

        ! We need a filepath
        call parlst_getvalue_string(rparlist,"ExtFE-POSTPROCESSING", &
                        "sIntFilePath",Intfilepath,bdequote=.TRUE.)
        if (Intfilepath .eq. '' ) then
             call output_line('You want to write out the &
                    &L1-results but have not specified a path', &
                    OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_init_postprocessing")
                    call sys_halt()
        end if

        ! Store the filepath
        rpostprocessing%Intfilepath = Intfilepath

    end if

    ! We will actually parse the L1RegionOfInterest here in a parser
    ! In component i is the region of interest for calculation i
    call fparser_create(rpostprocessing%pIntChiOmegaParser,nIntCalculations)

    ! Now get the arrays and fill them with data
    call storage_getbase_double(rpostprocessing%h_IntResults,p_Intresults)
    call storage_getbase_int2D(rpostprocessing%h_IntCompFunc,p_IntComp)
    call storage_getbase_char(rpostprocessing%h_IntChiOmega,p_IntChiOmega)
    call storage_getbase_char(rpostprocessing%h_IntTriFile,p_IntTriFile)
    call storage_getbase_int32(rpostprocessing%h_IntCubRule,p_IntCubRule)

    ! We will store which mesh is used
    call parlst_getvalue_string (rparlist,"ExtFE-FIRST",&
                                 "sMesh",sTriFileFirst,bdequote=.true.)
    call parlst_getvalue_string (rparlist,"ExtFE-SECOND",&
                                 "sMesh",sTriFileSecond,bdequote=.true.)

    ! Now read in the data
    do i=1,nIntCalculations

        ! Which components?
        ! fetch the whole line
        call parlst_getvalue_string(rparlist, &
                        "ExtFE-CALCULATIONS", "IntCalculations", &
                        sparam,sdefault="",isubstring=i)
        ! Now read the components from the string
        read(sparam,*) p_IntComp(1,i), p_IntComp(2,i)

        ! Which region of interest?
        call parlst_getvalue_string(rparlist, &
                        "ExtFE-CALCULATIONS", "IntRegionOfInterest", &
                        sparam,sdefault="",isubstring=i)
        ! Copy it to the postprocessing structure for
        ! some output. It is nasty but i did not find
        ! any other way
         do k=1,ExtFE_STRLEN
            p_IntChiOmega((i-1)*ExtFE_STRLEN+k:(i-1)*ExtFE_STRLEN+k)=sparam(k:k)
         end do

         ! Parse the string in the parser. Unfortunately we have to branch here
         ! due to the evaluation
         ! could be solved with a branch somewhere else but the code is easier
         ! to understand like this
        select case (nDIM)
             case(ExtFE_NDIM1)
               call fparser_parseFunction(rpostprocessing%pIntChiOmegaParser,i,sparam,(/'x'/))
             case(ExtFE_NDIM2)
               call fparser_parseFunction(rpostprocessing%pIntChiOmegaParser,i,sparam,(/'x','y'/))
             case default
               call output_line('Dimension of your problem is &
                     &not supported', &
                     OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_init_postprocessing")
               call sys_halt()
        end select

        ! Now we need a cubature rule. Fetch the string...
        call parlst_getvalue_string(rparlist, &
                        "ExtFE-CALCULATIONS", "IntCubatureRule", &
                        sparam,sdefault="",isubstring=i)
        ! and convert it to an ID
        p_IntCubRule(i) = cub_igetID(sparam)

        ! Store which mesh is to be used
        ! ||f-g|| => take first mesh
        if( (p_IntComp(1,i) .gt. 0) .and. &
            (p_Intcomp(2,i) .gt. 0)) then
          do k=1,ExtFE_STRLEN
            p_IntTriFile((i-1)*ExtFE_STRLEN+k:(i-1)*ExtFE_STRLEN+k)=sTriFileFirst(k:k)
          end do

        ! Second option: id1 >0, id2 <=0
        ! => calculate ||id1|| => mesh of first function
        else if((p_IntComp(1,i) .gt. 0) .and. &
                (p_IntComp(2,i) .le. 0)) then
          do k=1,ExtFE_STRLEN
            p_IntTriFile((i-1)*ExtFE_STRLEN+k:(i-1)*ExtFE_STRLEN+k)=sTriFileFirst(k:k)
          end do

        ! Third option: id1<=0, id2 >0
        ! => calculate ||id2|| => mesh of second function
        else if((p_IntComp(1,i) .le. 0) .and. &
                (p_IntComp(2,i) .gt. 0)) then
          do k=1,ExtFE_STRLEN
            p_IntTriFile((i-1)*ExtFE_STRLEN+k:(i-1)*ExtFE_STRLEN+k)=sTriFileSecond(k:k)
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

! Postprocess the L1-Results

subroutine ExtFE_postprocess_Integral(rpostprocessing)

    !<inputoutput>
    ! the postprocessing structure
    type(t_postprocessing), intent(inout) :: rpostprocessing
    !</inputoutput>

    ! For L1-Calculations
    integer :: nIntcalc
    integer, dimension(:,:), pointer :: iIntFuncComp => NULL()
    real(DP), dimension(:), pointer :: dIntResults => NULL()
    character, dimension(:), pointer :: cIntChiOmega => NULL()
    character, dimension(:), pointer :: cIntTriFile => NULL()

    ! General variables
    integer :: ifile,i,k
    character(LEN=ExtFE_STRLEN) :: soutString, stmpString, stmpString2

    ! Write out Integral-Calculation-Results?
    ! If not we are already done
    if (rpostprocessing%writeOutIntresults .eqv. .TRUE.) then
        !get pointer do the data arrays
        call storage_getbase_double(rpostprocessing%h_IntResults,dIntResults)
        call storage_getbase_int2D(rpostprocessing%h_IntCompFunc,iIntFuncComp)
        call storage_getbase_char(rpostprocessing%h_IntChiOmega,cIntChiOmega)
        call storage_getbase_char(rpostprocessing%h_IntTriFile,cIntTriFile)

        ! How many calculations?
        nIntcalc = rpostprocessing%nIntCalculations

        ! Open a file for writing out the results
        call io_openFileForWriting(rpostprocessing%Intfilepath, &
                    ifile,SYS_REPLACE,bformatted=.TRUE. )

        ! Ok, file is open, so we write
        do i=1,nIntcalc
            ! Empty the string - better save than sorry
            soutString = ''
            stmpString = ''
            ! Both components were used - so we need to
            ! write ||f_1(comp1) - f_2(comp2)||_L2(Omega_i) =
            ! into a string
            if((iIntFuncComp(1,i) .gt. 0) .AND. &
                (iIntFuncComp(2,i) .gt. 0)) then
                 write(stmpString,'(A11,I2,A8,I2,A11,I2,A3)') '\int{ f_1(', iIntFuncComp(1,i) , &
                                ') - f_2(', iIntFuncComp(2,i) ,') }_Omega_', i, ' = '
            ! Only first component was used - write
            ! ||f_1(comp1)||_L2(Omega_i) =  in the string
            else if((iIntFuncComp(1,i) .gt. 0 ) .AND. &
                  (iIntFuncComp(2,i) .le. 0) ) then
                 write(stmpString,'(A11,I2,A11,I2,A3)') '\int{ f_1(', iIntFuncComp(1,i) , &
                                ') }_Omega_', i, ' = '
            ! Only second component was used - write
            ! ||f_2(comp1)||_L2(Omega_i) =  in the string
            ! No sanity checks needed - done already during the init
            ! it can only be this case
            else
                 write(stmpString,'(A11,I2,A11,I2,A3)') '\int{ f_2(', iIntFuncComp(2,i) , &
                                ') }_Omega_', i, ' = '
            end if

            ! Add the result to the string
            write(soutString,'(A1,E16.10)') ' ', dIntResults(i)

            ! and write it in the file
            write(ifile,*) trim( trim(stmpString) // trim(soutString) )
        end do

        ! Now we write out what Omega_i actually is
        ! At the moment it should work  like this - if something is
        ! cut off I will rewrite this one.
        ! first a linebreak
        write(ifile,*) ''
        do i=1,nIntcalc
            ! Empty the strings
            stmpString = ''
            stmpString2 = ''
            write(stmpString,'(A11,I2,A14)') 'with Omega_', i ,' the subset of'

            do k = 1,ExtFE_STRLEN
                stmpString2(k:k) = cIntTriFile((i-1)*ExtFE_STRLEN+k)
            end do !k

            stmpString = trim(stmpString) // ' ' //trim(stmpString2)
            write(stmpString2,'(A24)') ' for that the expression'
            stmpString = trim(stmpString) // trim(stmpString2)

            do k = 1,ExtFE_STRLEN
                stmpString2(k:k) = cIntChiOmega((i-1)*ExtFE_STRLEN+k)
            end do !k

            write(ifile,'(A,A1,A,A10)') trim(stmpString), ' ', trim(stmpString2), ' returns 1'
        end do ! all L2Calculations


        ! Clean up
        close(ifile)
        dIntResults => NULL()
        iIntFuncComp => NULL()
        cIntChiOmega => NULL()
    end if


end subroutine

! Release of L1-Postprocessing

subroutine ExtFE_done_postprocessing_Integral(rpostprocessing)

!<input>
   type(t_postprocessing) , intent(inout):: rpostprocessing
!</input>

    ! Release all arrays from the storage that were allocated
    ! during the init. Without the If-condition it leads to an
    ! error since then the storage wants to deallocate handles
    ! that were not allocated
    if (rpostprocessing%h_IntCompFunc .gt. ST_NOHANDLE) then
        call storage_free(rpostprocessing%h_IntCompFunc)
    end if
    if (rpostprocessing%h_IntChiOmega .gt. ST_NOHANDLE) then
        call storage_free(rpostprocessing%h_IntChiOmega)
    end if
    if (rpostprocessing%h_IntResults .gt. ST_NOHANDLE) then
        call storage_free(rpostprocessing%h_IntResults)
    end if
    if (rpostprocessing%h_IntTriFile .gt. ST_NOHANDLE) then
        call storage_free(rpostprocessing%h_IntTriFile)
    end if
    if(rpostprocessing%h_IntCubRule .gt. ST_NOHANDLE) then
        call storage_free(rpostprocessing%h_IntCubRule)
    end if

    call fparser_release(rpostprocessing%pIntChiOmegaParser)
end subroutine


!-----------------------------------------------------!
! Init of the L1-Postprocessing

subroutine ExtFE_init_postprocessing_L1(rpostprocessing,rparlist)

    !<input>
    type(t_parlist), intent(in) :: rparlist
    !</input>

    !<inputoutput>
    type(t_postprocessing) , intent(inout):: rpostprocessing
    !</inputoutput>

    integer :: nL1Calculations, nL1ChiOmega, nL1CubRule
    character(LEN=ExtFE_STRLEN) :: L1filepath
    integer :: writeOutL1,i,k,nDIM

    ! pointer to the storage of the postprocessing structure
    real(DP), dimension(:),  pointer :: p_L1results =>NULL()
    integer, dimension (:,:), pointer :: p_L1Comp => NULL()
    character, dimension(:), pointer :: p_L1ChiOmega => NULL()
    character, dimension(:), pointer :: p_L1TriFile => NULL()
    integer(I32), dimension (:), pointer :: p_L1CubRule => NULL()
    character(LEN=ExtFE_STRLEN) :: sparam
    character(LEN=ExtFE_STRLEN) ::sTriFileFirst, sTriFileSecond




    !------------------------------------------!
    ! Set up the L1-Calculation-Postprocessing !
    !------------------------------------------!

    ! Find out how many L1-compuations to do
    nL1Calculations = parlst_querysubstrings(rparlist, &
                    "ExtFE-CALCULATIONS","L1calculations")
    ! How many Regions of interest are there?
    nL1ChiOmega = parlst_querysubstrings(rparlist, &
                    "ExtFE-CALCULATIONS","L1RegionOfInterest")
    ! How many cubature rules are there?
    nL1CubRule = parlst_querysubstrings(rparlist, &
                    "ExtFE-CALCULATIONS","L1CubatureRule")

  ! Validate that all numbers are the same
  if((nL1Calculations .ne. nL1ChiOmega ).or. &
     (nL1Calculations .ne. nL1CubRule) .or. &
     (nL1ChiOmega .ne. nL1CubRule)) then
     call output_line('Input Error: The number of calculations &
     &must match the number of cubature rules and regions of interest', &
     OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_init_postprocessing")
     call sys_halt()
  end if

  ! Get the dimension
  nDIM = rpostprocessing%nDim

  ! Init arrays for results and which component is used
  ! and store the filepath
  if (nL1Calculations .gt. 0) then

    ! Store the number in the structure
    rpostprocessing%nL1Calculations = nL1Calculations

    ! Array for components
    call storage_new("ExtFEcomparer_init_postprocessing", &
            "L1CompFunc", (/2,nL1Calculations/),ST_INT, &
             rpostprocessing%h_L1CompFunc,ST_NEWBLOCK_ZERO)

    ! Array for the results
    call storage_new("ExtFEcomparer_init_postprocessing", &
            "L1Results", nL1Calculations,ST_DOUBLE, &
             rpostprocessing%h_L1Results,ST_NEWBLOCK_ZERO)

    ! Array for the cubature rules
    call storage_new("ExtFEcomparer_init_postprocessing", &
            "L1CubRule",nL1Calculations,ST_INT32, &
            rpostprocessing%h_L1CubRule, ST_NEWBLOCK_ZERO)

    ! Array to store the region of interest for some
    ! postprocessing
    ! Only of interest if we write out the results in a file,
    ! but else we have some branches during the init and I
    ! don't want to do them right now.
    call storage_new("ExtFEcomparer_init_postprocessing", &
            "L1ChiOmega",ExtFE_STRLEN*nL1Calculations,ST_CHAR, &
            rpostprocessing%h_L1ChiOmega,ST_NEWBLOCK_ZERO)

    ! We want to write in the file which mesh was used
    ! Get an array for that
    call storage_new("ExtFEcomparer_init_postprocessing", &
          "L1TriFile",ExtFE_STRLEN*nL1Calculations,ST_CHAR, &
           rpostprocessing%h_L1TriFile,ST_NEWBLOCK_ZERO)

    ! Do we want to save the results in a file?
    call parlst_getvalue_int(rparlist,"ExtFE-POSTPROCESSING", &
                    "writeOutL1Calc",writeOutL1)

    if(writeOutL1 .eq. ExtFE_DO) then

        ! Tell the postprocessing that we want to write a file
        rpostprocessing%writeOutL1results = .true.

        ! We need a filepath
        call parlst_getvalue_string(rparlist,"ExtFE-POSTPROCESSING", &
                        "sL1FilePath",L1filepath,bdequote=.TRUE.)
        if (L1filepath .eq. '' ) then
             call output_line('You want to write out the &
                    &L1-results but have not specified a path', &
                    OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_init_postprocessing")
                    call sys_halt()
        end if

        ! Store the filepath
        rpostprocessing%L1filepath = L1filepath

    end if

    ! We will actually parse the L1RegionOfInterest here in a parser
    ! In component i is the region of interest for calculation i
    call fparser_create(rpostprocessing%pL1ChiOmegaParser,nL1Calculations)

    ! Now get the arrays and fill them with data
    call storage_getbase_double(rpostprocessing%h_L1Results,p_L1results)
    call storage_getbase_int2D(rpostprocessing%h_L1CompFunc,p_L1Comp)
    call storage_getbase_char(rpostprocessing%h_L1ChiOmega,p_L1ChiOmega)
    call storage_getbase_char(rpostprocessing%h_L1TriFile,p_L1TriFile)
    call storage_getbase_int32(rpostprocessing%h_L1CubRule,p_L1CubRule)

    ! We will store which mesh is used
    call parlst_getvalue_string (rparlist,"ExtFE-FIRST",&
                                 "sMesh",sTriFileFirst,bdequote=.true.)
    call parlst_getvalue_string (rparlist,"ExtFE-SECOND",&
                                 "sMesh",sTriFileSecond,bdequote=.true.)

    ! Now read in the data
    do i=1,nL1Calculations

        ! Which components?
        ! fetch the whole line
        call parlst_getvalue_string(rparlist, &
                        "ExtFE-CALCULATIONS", "L1calculations", &
                        sparam,sdefault="",isubstring=i)
        ! Now read the components from the string
        read(sparam,*) p_L1Comp(1,i), p_L1Comp(2,i)

        ! Which region of interest?
        call parlst_getvalue_string(rparlist, &
                        "ExtFE-CALCULATIONS", "L1RegionOfInterest", &
                        sparam,sdefault="",isubstring=i)
        ! Copy it to the postprocessing structure for
        ! some output. It is nasty but i did not find
        ! any other way
         do k=1,ExtFE_STRLEN
            p_L1ChiOmega((i-1)*ExtFE_STRLEN+k:(i-1)*ExtFE_STRLEN+k)=sparam(k:k)
         end do

         ! Parse the string in the parser. Unfortunately we have to branch here
         ! due to the evaluation
         ! could be solved with a branch somewhere else but the code is easier
         ! to understand like this
        select case (nDIM)
             case(ExtFE_NDIM1)
               call fparser_parseFunction(rpostprocessing%pL1ChiOmegaParser,i,sparam,(/'x'/))
             case(ExtFE_NDIM2)
               call fparser_parseFunction(rpostprocessing%pL1ChiOmegaParser,i,sparam,(/'x','y'/))
             case default
               call output_line('Dimension of your problem is &
                     &not supported', &
                     OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_init_postprocessing")
               call sys_halt()
        end select

        ! Now we need a cubature rule. Fetch the string...
        call parlst_getvalue_string(rparlist, &
                        "ExtFE-CALCULATIONS", "L1CubatureRule", &
                        sparam,sdefault="",isubstring=i)
        ! and convert it to an ID
        p_L1CubRule(i) = cub_igetID(sparam)

        ! Store which mesh is to be used
        ! ||f-g|| => take first mesh
        if( (p_L1Comp(1,i) .gt. 0) .and. &
            (p_L1comp(2,i) .gt. 0)) then
          do k=1,ExtFE_STRLEN
            p_L1TriFile((i-1)*ExtFE_STRLEN+k:(i-1)*ExtFE_STRLEN+k)=sTriFileFirst(k:k)
          end do

        ! Second option: id1 >0, id2 <=0
        ! => calculate ||id1|| => mesh of first function
        else if((p_L1Comp(1,i) .gt. 0) .and. &
                (p_L1Comp(2,i) .le. 0)) then
          do k=1,ExtFE_STRLEN
            p_L1TriFile((i-1)*ExtFE_STRLEN+k:(i-1)*ExtFE_STRLEN+k)=sTriFileFirst(k:k)
          end do

        ! Third option: id1<=0, id2 >0
        ! => calculate ||id2|| => mesh of second function
        else if((p_L1Comp(1,i) .le. 0) .and. &
                (p_L1Comp(2,i) .gt. 0)) then
          do k=1,ExtFE_STRLEN
            p_L1TriFile((i-1)*ExtFE_STRLEN+k:(i-1)*ExtFE_STRLEN+k)=sTriFileSecond(k:k)
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

! Postprocess the L1-Results

subroutine ExtFE_postprocess_L1(rpostprocessing)

    !<inputoutput>
    ! the postprocessing structure
    type(t_postprocessing), intent(inout) :: rpostprocessing
    !</inputoutput>

    ! For L1-Calculations
    integer :: nL1calc
    integer, dimension(:,:), pointer :: iL1FuncComp => NULL()
    real(DP), dimension(:), pointer :: dL1Results => NULL()
    character, dimension(:), pointer :: cL1ChiOmega => NULL()
    character, dimension(:), pointer :: cL1TriFile => NULL()

    ! General variables
    integer :: ifile,i,k
    character(LEN=ExtFE_STRLEN) :: soutString, stmpString, stmpString2

    ! Write out L1-Calculation-Results?
    ! If not we are already done
    if (rpostprocessing%writeOutL1results .eqv. .TRUE.) then
        !get pointer do the data arrays
        call storage_getbase_double(rpostprocessing%h_L1Results,dL1Results)
        call storage_getbase_int2D(rpostprocessing%h_L1CompFunc,iL1FuncComp)
        call storage_getbase_char(rpostprocessing%h_L1ChiOmega,cL1ChiOmega)
        call storage_getbase_char(rpostprocessing%h_L1TriFile,cL1TriFile)

        ! How many calculations?
        nL1calc = rpostprocessing%nL1Calculations

        ! Open a file for writing out the results
        call io_openFileForWriting(rpostprocessing%L1filepath, &
                    ifile,SYS_REPLACE,bformatted=.TRUE. )

        ! Ok, file is open, so we write
        do i=1,nL1calc
            ! Empty the string - better save than sorry
            soutString = ''
            stmpString = ''
            ! Both components were used - so we need to
            ! write ||f_1(comp1) - f_2(comp2)||_L2(Omega_i) =
            ! into a string
            if((iL1FuncComp(1,i) .gt. 0) .AND. &
                (iL1FuncComp(2,i) .gt. 0)) then
                 write(stmpString,'(A7,I2,A8,I2,A13,I2,A4)') '||f_1(', iL1FuncComp(1,i) , &
                                ') - f_2(', iL1FuncComp(2,i) ,')||_L1(Omega_', i, ') = '
            ! Only first component was used - write
            ! ||f_1(comp1)||_L2(Omega_i) =  in the string
            else if((iL1FuncComp(1,i) .gt. 0 ) .AND. &
                  (iL1FuncComp(2,i) .le. 0) ) then
                 write(stmpString,'(A7,I2,A13,I2,A4)') '||f_1(', iL1FuncComp(1,i) , &
                                ')||_L1(Omega_', i, ') = '
            ! Only second component was used - write
            ! ||f_2(comp1)||_L2(Omega_i) =  in the string
            ! No sanity checks needed - done already during the init
            ! it can only be this case
            else
                 write(stmpString,'(A7,I2,A13,I2,A4)') '||f_2(', iL1FuncComp(2,i) , &
                                ')||_L1(Omega_', i, ') = '
            end if

            ! Add the result to the string
            write(soutString,'(A1,E16.10)') ' ', dL1Results(i)

            ! and write it in the file
            write(ifile,*) trim( trim(stmpString) // trim(soutString) )
        end do

        ! Now we write out what Omega_i actually is
        ! At the moment it should work  like this - if something is
        ! cut off I will rewrite this one.
        ! first a linebreak
        write(ifile,*) ''
        do i=1,nL1calc
            ! Empty the strings
            stmpString = ''
            stmpString2 = ''
            write(stmpString,'(A11,I2,A14)') 'with Omega_', i ,' the subset of'

            do k = 1,ExtFE_STRLEN
                stmpString2(k:k) = cL1TriFile((i-1)*ExtFE_STRLEN+k)
            end do !k

            stmpString = trim(stmpString) // ' ' //trim(stmpString2)
            write(stmpString2,'(A24)') ' for that the expression'
            stmpString = trim(stmpString) // trim(stmpString2)

            do k = 1,ExtFE_STRLEN
                stmpString2(k:k) = cL1ChiOmega((i-1)*ExtFE_STRLEN+k)
            end do !k

            write(ifile,'(A,A1,A,A10)') trim(stmpString), ' ', trim(stmpString2), ' returns 1'
        end do ! all L2Calculations


        ! Clean up
        close(ifile)
        dL1Results => NULL()
        iL1FuncComp => NULL()
        cL1ChiOmega => NULL()
    end if


end subroutine

! Release of L1-Postprocessing

subroutine ExtFE_done_postprocessing_L1(rpostprocessing)

!<input>
   type(t_postprocessing) , intent(inout):: rpostprocessing
!</input>

    ! Release all arrays from the storage that were allocated
    ! during the init. Without the If-condition it leads to an
    ! error since then the storage wants to deallocate handles
    ! that were not allocated
    if (rpostprocessing%h_L1CompFunc .gt. ST_NOHANDLE) then
        call storage_free(rpostprocessing%h_L1CompFunc)
    end if
    if (rpostprocessing%h_L1ChiOmega .gt. ST_NOHANDLE) then
        call storage_free(rpostprocessing%h_L1ChiOmega)
    end if
    if (rpostprocessing%h_L1Results .gt. ST_NOHANDLE) then
        call storage_free(rpostprocessing%h_L1Results)
    end if
    if (rpostprocessing%h_L1TriFile .gt. ST_NOHANDLE) then
        call storage_free(rpostprocessing%h_L1TriFile)
    end if
    if(rpostprocessing%h_L1CubRule .gt. ST_NOHANDLE) then
        call storage_free(rpostprocessing%h_L1CubRule)
    end if

    call fparser_release(rpostprocessing%pL1ChiOmegaParser)
end subroutine


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
    ! How many Regions of interest are there?
    nL2ChiOmega = parlst_querysubstrings(rparlist, &
                    "ExtFE-CALCULATIONS","L2RegionOfInterest")
    ! How many cubature rules are there?
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

  ! Init arrays for results and which component is used
  ! and store the filepath
  if (nL2Calculations .gt. 0) then

    ! Store the number in the structure
    rpostprocessing%nL2Calculations = nL2Calculations

    ! Array for components
    call storage_new("ExtFEcomparer_init_postprocessing", &
            "L2CompFunc", (/2,nL2Calculations/),ST_INT, &
             rpostprocessing%h_L2CompFunc,ST_NEWBLOCK_ZERO)

    ! Array for the results
    call storage_new("ExtFEcomparer_init_postprocessing", &
            "L2Results", nL2Calculations,ST_DOUBLE, &
             rpostprocessing%h_L2Results,ST_NEWBLOCK_ZERO)

    ! Array for the cubature rules
    call storage_new("ExtFEcomparer_init_postprocessing", &
            "L2CubRule",nL2Calculations,ST_INT32, &
            rpostprocessing%h_L2CubRule, ST_NEWBLOCK_ZERO)

    ! Array to store the region of interest for some
    ! postprocessing
    ! Only of interest if we write out the results in a file,
    ! but else we have some branches during the init and I
    ! don't want to do them right now.
    call storage_new("ExtFEcomparer_init_postprocessing", &
            "L2ChiOmega",ExtFE_STRLEN*nL2Calculations,ST_CHAR, &
            rpostprocessing%h_L2ChiOmega,ST_NEWBLOCK_ZERO)

    ! We want to write in the file which mesh was used
    ! Get an array for that
    call storage_new("ExtFEcomparer_init_postprocessing", &
          "L2TriFile",ExtFE_STRLEN*nL2Calculations,ST_CHAR, &
           rpostprocessing%h_L2TriFile,ST_NEWBLOCK_ZERO)

    ! Do we want to save the results in a file?
    call parlst_getvalue_int(rparlist,"ExtFE-POSTPROCESSING", &
                    "writeOutL2Calc",writeOutL2)

    if(writeOutL2 .eq. ExtFE_DO) then

        ! Tell the postprocessing that we want to write a file
        rpostprocessing%writeOutL2results = .true.

        ! We need a filepath
        call parlst_getvalue_string(rparlist,"ExtFE-POSTPROCESSING", &
                        "sL2FilePath",L2filepath,bdequote=.TRUE.)
        if (L2filepath .eq. '' ) then
             call output_line('You want to write out the &
                    &L2-results but have not specified a path', &
                    OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_init_postprocessing")
                    call sys_halt()
        end if

        ! Store the filepath
        rpostprocessing%L2filepath = L2filepath

    end if

    ! We will actually parse the L2RegionOfInterest here in a parser
    ! In component i is the region of interest for calculation i
    call fparser_create(rpostprocessing%pL2ChiOmegaParser,nL2Calculations)

    ! Now get the arrays and fill them with data
    call storage_getbase_double(rpostprocessing%h_L2Results,p_L2results)
    call storage_getbase_int2D(rpostprocessing%h_L2CompFunc,p_L2Comp)
    call storage_getbase_char(rpostprocessing%h_L2ChiOmega,p_L2ChiOmega)
    call storage_getbase_char(rpostprocessing%h_L2TriFile,p_L2TriFile)
    call storage_getbase_int32(rpostprocessing%h_L2CubRule,p_L2CubRule)

    ! We will store which mesh is used
    call parlst_getvalue_string (rparlist,"ExtFE-FIRST",&
                                 "sMesh",sTriFileFirst,bdequote=.true.)
    call parlst_getvalue_string (rparlist,"ExtFE-SECOND",&
                                 "sMesh",sTriFileSecond,bdequote=.true.)

    ! Now read in the data
    do i=1,nL2Calculations

        ! Which components?
        ! fetch the whole line
        call parlst_getvalue_string(rparlist, &
                        "ExtFE-CALCULATIONS", "L2calculations", &
                        sparam,sdefault="",isubstring=i)
        ! Now read the components from the string
        read(sparam,*) p_L2Comp(1,i), p_L2Comp(2,i)

        ! Which region of interest?
        call parlst_getvalue_string(rparlist, &
                        "ExtFE-CALCULATIONS", "L2RegionOfInterest", &
                        sparam,sdefault="",isubstring=i)
        ! Copy it to the postprocessing structure for
        ! some output. It is nasty but i did not find
        ! any other way
         do k=1,ExtFE_STRLEN
            p_L2ChiOmega((i-1)*ExtFE_STRLEN+k:(i-1)*ExtFE_STRLEN+k)=sparam(k:k)
         end do

         ! Parse the string in the parser. Unfortunately we have to branch here
         ! due to the evaluation
         ! could be solved with a branch somewhere else but the code is easier
         ! to understand like this
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

        ! Now we need a cubature rule. Fetch the string...
        call parlst_getvalue_string(rparlist, &
                        "ExtFE-CALCULATIONS", "L2CubatureRule", &
                        sparam,sdefault="",isubstring=i)
        ! and convert it to an ID
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
    ! If not we are already done
    if (rpostprocessing%writeOutL2results .eqv. .TRUE.) then
        !get pointer do the data arrays
        call storage_getbase_double(rpostprocessing%h_L2Results,dL2Results)
        call storage_getbase_int2D(rpostprocessing%h_L2CompFunc,iL2FuncComp)
        call storage_getbase_char(rpostprocessing%h_L2ChiOmega,cL2ChiOmega)
        call storage_getbase_char(rpostprocessing%h_L2TriFile,cL2TriFile)

        ! How many calculations?
        nL2calc = rpostprocessing%nL2Calculations

        ! Open a file for writing out the results
        call io_openFileForWriting(rpostprocessing%L2filepath, &
                    ifile,SYS_REPLACE,bformatted=.TRUE. )

        ! Ok, file is open, so we write
        do i=1,nL2calc
            ! Empty the string - better save than sorry
            soutString = ''
            stmpString = ''
            ! Both components were used - so we need to
            ! write ||f_1(comp1) - f_2(comp2)||_L2(Omega_i) =
            ! into a string
            if((iL2FuncComp(1,i) .gt. 0) .AND. &
                (iL2FuncComp(2,i) .gt. 0)) then
                 write(stmpString,'(A7,I2,A8,I2,A13,I2,A4)') '||f_1(', iL2FuncComp(1,i) , &
                                ') - f_2(', iL2FuncComp(2,i) ,')||_L2(Omega_', i, ') = '
            ! Only first component was used - write
            ! ||f_1(comp1)||_L2(Omega_i) =  in the string
            else if((iL2FuncComp(1,i) .gt. 0 ) .AND. &
                  (iL2FuncComp(2,i) .le. 0) ) then
                 write(stmpString,'(A7,I2,A13,I2,A4)') '||f_1(', iL2FuncComp(1,i) , &
                                ')||_L2(Omega_', i, ') = '
            ! Only second component was used - write
            ! ||f_2(comp1)||_L2(Omega_i) =  in the string
            ! No sanity checks needed - done already during the init
            ! it can only be this case
            else
                 write(stmpString,'(A7,I2,A13,I2,A4)') '||f_2(', iL2FuncComp(2,i) , &
                                ')||_L2(Omega_', i, ') = '
            end if

            ! Add the result to the string
            write(soutString,'(A1,E16.10)') ' ', dL2Results(i)

            ! and write it in the file
            write(ifile,*) trim( trim(stmpString) // trim(soutString) )
        end do

        ! Now we write out what Omega_i actually is
        ! At the moment it should work  like this - if something is
        ! cut off I will rewrite this one.
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

    ! Release all arrays from the storage that were allocated
    ! during the init. Without the If-condition it leads to an
    ! error since then the storage wants to deallocate handles
    ! that were not allocated
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

!-------------------------------------------------------------------------

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

        ! Create the storage for the information
        ! Result(i) = Result of calculation i
        call storage_new("ExtFEcomparer_init_postprocessing", &
            "PointvalueResults",nPointEvaluations,ST_DOUBLE,&
                rpostprocessing%h_PointResults,ST_NEWBLOCK_ZERO)

        ! To make the next allocation easier - write up in the array
        ! how much space we need
        iSizeArray = (/4,nPointEvaluations/)

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

        ! Do we save the results in a file?
        call parlst_getvalue_int(rparlist,"ExtFE-POSTPROCESSING", &
                    "writeOutPointValues",writeOutPoint)

        if(writeOutPoint .eq. ExtFE_DO) then
            ! tell the postprocessing structure that we want to save
            rpostprocessing%writeOutPointCalucations = .true.

            ! get the filepath for the file
            call parlst_getvalue_string(rparlist,"ExtFE-POSTPROCESSING", &
                        "sPointValuesFilePath",PointFilepath,bdequote=.TRUE.)
            if (PointFilepath .eq. '' ) then
                call output_line('You want to write out the &
                        &results of the point evaluations but have not specified a path', &
                        OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_init_ppostprocessing")
                call sys_halt()
            end if

        rpostprocessing%PointFilepath = PointFilepath

        end if ! writeOutPoint

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
            !and read it in. Which function comp and deriv and which point?
            ! We have to branch since in ie 1D there is no y-component
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

    ! We want to make the output good-looking so we write what
    ! derivative was evaluated
    ! We start to count the derivatives from 0, so
    ! we need to shift by 1
    sderivnames(ExtFE_DER_FUNC_3D+1) = '      '
    sderivnames(ExtFE_DER_DERIV_3D_X+1) = '(d/dx)'
    sderivnames(ExtFE_DER_DERIV_3D_Y+1) = '(d/dy)'
    sderivnames(ExtFE_DER_DERIV_3D_Z+1) = '(d/dz)'


    ! Write out Pointvalue-Calculation-Results?
    ! If not we are done here
    if(rpostprocessing%writeOutPointCalucations .eqv. .TRUE.) then
        ! Get the arrays
        call storage_getbase_double2D(rpostprocessing%h_PointCoordinates,dPointCoordinates)
        call storage_getbase_double(rpostprocessing%h_PointResults,dPointValues)
        call storage_getbase_int2D(rpostprocessing%h_PointFuncComponents,iPointFuncComp)

        ! find out how many calculations we have
        nPointEvals = rpostprocessing%nPointCalculations

        ! Open a file for writing out the results
        call io_openFileForWriting(rpostprocessing%PointFilepath, &
                    ifile,SYS_REPLACE,bformatted=.TRUE. )

        do i=1,nPointEvals

            ! Empty the string - better save than sorry
            stmpString = ''
            stmpString2 = ''
            soutString = ''

            ! Find out what we have
            ! Both functions => write (deriv f_1(comp1) - deriv f_2(comp2))
            ! in the string
            if((iPointFuncComp(1,i) .gt. 0) .AND. &
               (iPointFuncComp(3,i) .gt. 0)) then
                  write(stmpString,'(A3,I2,A3)') 'f1_', iPointFuncComp(1,i), ' - '
                  write(stmpString2,'(A3,I2)') 'f2_', iPointFuncComp(3,i)
                  stmpString = trim(sderivnames(iPointFuncComp(2,i)+1)) // trim(stmpString)
                  stmpString = '(' // trim(stmpString)
                  stmpString2 = trim(sderivnames(iPointFuncComp(4,i)+1)) // trim(stmpString2)
                  stmpString2 = trim(stmpString2) // ')'
                  stmpString = trim(stmpString) // trim(stmpString2)
            ! Only first => write (deriv f_1(comp1)) in the string
            else if ((iPointFuncComp(1,i) .gt. 0) .AND. &
                    (iPointFuncComp(3,i) .le. 0)) then
                  write(stmpString,'(A3,I2)') 'f1_', iPointFuncComp(1,i)
                  stmpString = (sderivnames(iPointFuncComp(2,i)+1)) // trim(stmpString)
                  stmpString = '(' // trim(stmpString)
                  stmpString = trim(stmpString) // ')'
            ! No sanity check needed - already done during the calculation
            ! it can only be the case: only the second
            ! => write (deriv f_2(comp2)) in the string
            else
                  write(stmpString,'(A3,I2)') 'f2_', iPointFuncComp(3,i)
                  stmpString = (sderivnames(iPointFuncComp(4,i)+1)) // trim(stmpString)
                  stmpString = '(' // trim(stmpString)
                  stmpString = trim(stmpString) // ')'
            end if

            ! Now we are still not done, we have to add the
            ! evaluation point and the value. The evaluation point
            ! depends on the dimension
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

            ! Now we have everything - write it out in the file
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

    ! deallocate all arrays that we allocated during the init.
    ! Without if we get an error because the storage then tries
    ! to deallocate arrays that were not allocated
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
    integer :: handle_UCD_AddTypeOrig
    integer :: handle_UCD_AddElemProject
    integer :: handle_UCD_VecVars
    integer :: handle_UCD_ScalarVars


    ! We need the dim everywhere - save it
    nDim = rpostprocessing%nDim

    !--------------------------------------!
    ! Set up the UCD-Postprocessing        !
    !--------------------------------------!

    ! Write out meshes?
    call parlst_getvalue_int(rparlist,"ExtFE-POSTPROCESSING", &
                    "writeOutMeshesUCD",writeOutMeshes)
    if(writeOutMeshes .eq. ExtFE_DO) then
        ! We need a filepath
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
        ! We need a filepath
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

        ! Now we set these temporary values. What do they do/mean?
        ! We have wrapper for the init.
        ! With the strings we tell where in the parlist we have to search.
        ! The handles are set to NOHANDLE, the wrapper will return us the handle
        ! If this is not done and we i.e. have no scalars, we will not do anything
        ! with the variable handle_UCD_ScalarVars. After the wrapper we just save
        ! everything in the structure. So it might happen that if the compiler initially
        ! sets the handle_UCD_ScalarVars to 2, we save in the structure that in handle
        ! 2 all scalar variables are. this will lead to unpredictable behaviour
        ! For the init, the wrapper needs to know how many variables are in the vector

        sUCD_AddElemProjection = "sUCDElemProjectFirstFunction"
        sUCD_AddTypeOrig = "sUCDaddTypeFirstFunction"
        sUCD_VecComp = "sUCDVecCompFirstFunction"
        sUCD_ScalarComp = "UCD_ScalarFirstOrig"
        nVars = rpostprocessing%nVarFirst
        handle_UCD_AddElemProject = ST_NOHANDLE
        handle_UCD_AddTypeOrig = ST_NOHANDLE
        handle_UCD_VecVars = ST_NOHANDLE
        handle_UCD_ScalarVars = ST_NOHANDLE

        call ExtFE_init_pp_UCD_fefuncout(handle_UCD_AddElemProject,handle_UCD_AddTypeOrig,&
                                    handle_UCD_VecVars,handle_UCD_ScalarVars,&
                            sUCD_AddElemProjection,sUCD_AddTypeOrig,&
                                        sUCD_VecComp,sUCD_ScalarComp,&
                                        nDim,&
                                        nVars,&
                                        rparlist)
        ! Now save the handles in the structure.
        rpostprocessing%h_UCD_AddTypeOrigFirst = handle_UCD_AddTypeOrig
        rpostprocessing%h_UCD_AddElemProjectFirst = handle_UCD_AddElemProject
        rpostprocessing%h_UCD_VecsFirstOrig = handle_UCD_VecVars
        rpostprocessing%h_UCD_ScalarFirstOrig = handle_UCD_ScalarVars


    end if

    ! Write out second original FE-Function?
    call parlst_getvalue_int(rparlist,"ExtFE-POSTPROCESSING", &
                    "writeOutSecondFunctionUCD",writeOutSecondFunction)
    if(writeOutSecondFunction .eq. ExtFE_DO) then
        ! If yes, we need a filepath
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

        ! Now we set these temporary values. What do they do/mean?
        ! We have wrapper for the init.
        ! With the strings we tell where in the parlist we have to search.
        ! The handles are set to NOHANDLE, the wrapper will return us the handle
        ! If this is not done and we i.e. have no scalars, we will not do anything
        ! with the variable handle_UCD_ScalarVars. After the wrapper we just save
        ! everything in the structure. So it might happen that if the compiler initially
        ! sets the handle_UCD_ScalarVars to 2, we save in the structure that in handle
        ! 2 all scalar variables are. this will lead to unpredictable behaviour
        ! For the init, the wrapper needs to know how many variables are in the vector
        sUCD_AddElemProjection = "sUCDElemProjectSecondFunction"
        sUCD_AddTypeOrig = "sUCDaddTypeSecondFunction"
        sUCD_VecComp = "sUCDVecCompSecondFunction"
        sUCD_ScalarComp = "UCD_ScalarSecondOrig"
        nVars = rpostprocessing%nVarSecond
        handle_UCD_AddElemProject = ST_NOHANDLE
        handle_UCD_AddTypeOrig = ST_NOHANDLE
        handle_UCD_VecVars = ST_NOHANDLE
        handle_UCD_ScalarVars = ST_NOHANDLE

        call ExtFE_init_pp_UCD_fefuncout(handle_UCD_AddElemProject,handle_UCD_AddTypeOrig,&
                                    handle_UCD_VecVars,handle_UCD_ScalarVars,&
                            sUCD_AddElemProjection,sUCD_AddTypeOrig,&
                                        sUCD_VecComp,sUCD_ScalarComp,&
                                        nDim,&
                                        nVars,&
                                        rparlist)

        ! Save the handles in the structure
        rpostprocessing%h_UCD_AddTypeOrigSecond = handle_UCD_AddTypeOrig
        rpostprocessing%h_UCD_AddElemProjectSecond = handle_UCD_AddElemProject
        rpostprocessing%h_UCD_VecsSecondOrig = handle_UCD_VecVars
        rpostprocessing%h_UCD_ScalarSecondOrig = handle_UCD_ScalarVars


    end if


    ! IF we want to do some ucd output, we read in which format (VTK,...)
    ! and which type (Standard, ...)
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
    ! Which Dimension do we have
    integer, intent(in) :: nDimension
    ! How many variables/components are in the function?
    integer, intent(in) :: nVarsFunc

    ! InOut: Handles to the storage of the arrays
    integer, intent(inout) :: h_UCD_AddTypeOrig
    integer, intent(inout) :: h_UCD_AddElemProject
    integer, intent(inout) :: h_UCD_VecVars
    integer, intent(inout) :: h_UCD_ScalarVars

    ! Local variables
    integer :: i,k, iVarsScalar, nVecs
    ! After the init, we want to do something with the arrays, so we
    ! need some pointers for them
    integer, dimension (:,:), pointer :: p_IntPointerVecs => NULL()
    integer, dimension(:), pointer :: p_IntPointerScalars => NULL()
    integer, dimension(:), pointer :: p_IntPointerAddElemProject => NULL()
    integer, dimension(:), pointer :: p_IntPointerAddType => NULL()
    integer, dimension(:), allocatable :: VarsNotVector

    character(LEN=ExtFE_STRLEN) :: sparam

    ! Allocate the arrays for the UCD-Output
    ! Therefore we need the number of variables in the function
    ! and how many vectors are in there

    ! The add-type (element-based or vertex-based)?
    call storage_new("ExtFE_init_postprocessing",sUCD_AddTypeOrig,&
                      nVarsFunc,ST_INT,h_UCD_AddTypeOrig,&
                      ST_NEWBLOCK_ZERO)
    ! Project them linear or constant?
    call storage_new("ExtFE_init_postprocessing",sUCD_AddElemProjection,&
                      nVarsFunc,ST_INT,h_UCD_AddElemProject,&
                      ST_NEWBLOCK_ZERO)
    ! How many vectors do we have?
    nVecs = parlst_querysubstrings(rparlist, &
                  "ExtFE-POSTPROCESSING",sUCD_VecComp)

    ! If we have vectors we have to take care of them
    if(nVecs .gt. 0) then
        ! An array to store which components are in which vector
        call storage_new("ExtFE_init_postprocessing",sUCD_VecComp,&
                        (/nVecs,nDimension/),ST_INT,h_UCD_VecVars,&
                        ST_NEWBLOCK_ZERO)

        ! Get a pointer to the array
        call storage_getbase_int2D(h_UCD_VecVars,p_IntPointerVecs)
        do i=1,nVecs
            ! fetch the whole line and then read it in
            call parlst_getvalue_string(rparlist, &
                          "ExtFE-Postprocessing", sUCD_VecComp, &
                            sparam,sdefault="",isubstring=i)
            read(sparam,*) p_IntPointerVecs(i,:)

        end do ! Read in the Data
        ! Done with that one

        ! Now figure out which scalar variables are left
        ! Assumtion: Every variable not in a vector shall be written out as scalar variable
        ! We know we have maximum nVar variables
        ! Idea: Init an array with 1, then set each entry to 0 where we have a vector component
        ! Then count what is left
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

        ! Do we have scalar variables?
        if(iVarsScalar .gt. 0) then
            ! We need an array for them
            call storage_new("ExtFE_init_postprocessing",sUCD_ScalarComp,&
                   iVarsScalar,ST_INT,h_UCD_ScalarVars,&
                    ST_NEWBLOCK_ZERO)
            call storage_getbase_int(h_UCD_ScalarVars,p_IntPointerScalars)
            ! Now we want to write in the array. Problem: We have to do
            ! some indexing. Easy solution: loop over all components and count with
            ! another variable
            k = 0
            do i=1,nVarsFunc
                if (VarsNotVector(i) .eq. 1) then
                    k = k+1
                    p_IntPointerScalars(k) = i
                end if
            end do

        end if
        deallocate(VarsNotVector)

    ! Else we have only scalar variables
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

    ! Now how to add them (element/vertex-based)
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

    ! Which projection?
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

        ! We made a wrapper. give all handles to the arrays + the filepath
        ! + the UCD-Format + the UCD-Style to the wrapper - thats it
        call ExtFE_write_UCD(rpostprocessing%UCD_feFunction_first_orig,&
                            rpostprocessing%UCD_Format,&
                            rpostprocessing%UCD_Style,&
                            rpostprocessing%UCD_FEfunctionOneOrigOutPath,&
                            rpostprocessing%h_UCD_AddTypeOrigFirst,&
                            rpostprocessing%h_UCD_ScalarFirstOrig,&
                            rpostprocessing%h_UCD_VecsFirstOrig)
    end if

    if(rpostprocessing%ucd_OUT_orig_functions_two .eqv. .TRUE. ) then

        ! We made a wrapper. give all handles to the arrays + the filepath
        ! + the UCD-Format + the UCD-Style to the wrapper - thats it
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
    ! Get the array
     call storage_getbase_int(h_scalarValues,p_scalarValues)
     nScalars = ubound(p_scalarValues,1)
     do i=1,nScalars
        call lsyssc_getbase_double(fefunction%RvectorBlock(p_scalarValues(i)),p_DataX)
         select case(p_addType(p_scalarValues(i)))

          case(ExtFE_UCD_VERT_BASED)
             !call ucd_addVarVertBasedVec(rexportFunction,"SclarVar_"//trim(adjustl(svarname)),p_DataX)
             call ucd_addVariableVertexBased(rexportFunction,"ScalarVar_"//trim(adjustl(sys_siL(i,3))),p_DataX)
          case(ExtFE_UCD_ELEM_BASED)
             !call ucd_addVarElemBasedVec(rexportFunction,"ScalarVar_"//trim(adjustl(svarname)),p_DataX)
             call ucd_addVariableElementBased(rexportFunction,"ScalarVar_"//trim(adjustl(sys_siL(i,3))),p_DataX)
        end select

     end do
   end if

   ! Now the vector valued
   ! Assumtion: All Variables in the vector are added the same way!
   ! first on vertex-based => All vertex based etc
   ! We don't need sanity checks as we did them during the init
   ! Everytime the same: Get how to add them, then call the ucd_addVarxyz with
   ! the structure and variable name and data_x,data_y, data_z (depending on the dimension)
   if(h_VectorValues .gt. ST_NOHANDLE) then
        call storage_getbase_INT2D(h_VectorValues,p_VectorValues)
        nVecs = ubound(p_VectorValues,1)
        nDim = fefunction%p_rblockDiscr%ndimension
        svarname = ''
        do i=1,nVecs
            svarname = "VectorVar_" //trim(adjustl(sys_siL(i,3)))
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
    if(associated(rpostprocessing%UCDBlockDiscrFirst)) then
        if(rpostprocessing%UCDBlockDiscrFirst%ndimension .gt. 0) then
           call spdiscr_releaseBlockDiscr(rpostprocessing%UCDBlockDiscrFirst)
        end if
        deallocate(rpostprocessing%UCDBlockDiscrFirst)
    end if



    if(associated(rpostprocessing%UCD_feFunction_first_orig)) then
        if(rpostprocessing%UCD_feFunction_first_orig%h_Ddata .gt. ST_NOHANDLE) then
            call lsysbl_releaseVector(rpostprocessing%UCD_feFunction_first_orig)
        end if
        deallocate(rpostprocessing%UCD_feFunction_first_orig)
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
    if(associated(rpostprocessing%UCDBlockDiscrSecond)) then
        if(rpostprocessing%UCDBlockDiscrSecond%ndimension .gt. 0) then
            call spdiscr_releaseBlockDiscr(rpostprocessing%UCDBlockDiscrSecond)
        end if
        deallocate(rpostprocessing%UCDBlockDiscrSecond)
     end if

    if(associated(rpostprocessing%UCD_feFunction_second_orig)) then
        if(rpostprocessing%UCD_feFunction_second_orig%h_Ddata .gt. ST_NOHANDLE) then
            call lsysbl_releaseVector(rpostprocessing%UCD_feFunction_second_orig)
        end if
        deallocate(rpostprocessing%UCD_feFunction_second_orig)
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
        ! If we want to do that, we need a path where to save
        call parlst_getvalue_string(rparlist,"ExtFE-POSTPROCESSING", &
                        "sOutPathFirstOrigVec",sPathFirst,bdequote=.TRUE.)
        if (sPathFirst .eq. '') then
             call output_line('You want to convert the &
                    &first vector but have not specified a path to save it', &
                    OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_init_postprocessing")
                    call sys_halt()
        end if

        ! and we need a format for writing.
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
        ! IF we want to do that we need a path
        call parlst_getvalue_string(rparlist,"ExtFE-POSTPROCESSING", &
                        "sOutPathSecondOrigVec", sPathSecond ,bdequote=.TRUE.)

        if (sPathSecond .eq. '') then
             call output_line('You want to convert the &
                    &second vector but have not specified a path to save it', &
                    OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_init_postprocessing")
                    call sys_halt()
        end if

        ! And a format for writing
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
        ! Call vecio_writeBlockVectorHR with the according arguments
        call vecio_writeBlockVectorHR(rpostprocessing%OrigVecFirst,'SOLUTION', &
                bunsorted,0,rpostprocessing%sOrigVecPathOutFirst,&
                sformat=trim(rpostprocessing%sOrigVec1OutFMT),&
                scomment=trim(adjustl(comment)))
    end if

    if(rpostprocessing%writeOutOrigVector2 .eqv. .TRUE.) then
        call vecio_writeBlockVectorHR(rpostprocessing%OrigVecSecond,'SOLUTION', &
                bunsorted,0,rpostprocessing%sOrigVecPathOutSecond,&
                sformat=trim(rpostprocessing%sOrigVec2OutFMT),&
                scomment=trim(adjustl(comment)))
    end if


end subroutine


subroutine ExtFE_done_postprocessing_OutOrigVec(rpostprocessing)

!<input>
   type(t_postprocessing) , intent(inout):: rpostprocessing
!</input>
    ! Just some pointers to set to NULL()
    rpostprocessing%OrigVecFirst => NULL()
    rpostprocessing%OrigVecSecond => NULL()

end subroutine

end module
