module ExtFEcomparer_vector

  !<description>
  ! In this module we take care of loading the vector
  ! in the problem structure and combine it with the
  ! discretisation so that it becomes a FE-function
  !</description>

  ! Include basic Feat-2 modules
  use fsystem
  use storage
  use genoutput
  use vectorio
  use linearsystemblock
  use linearsystemscalar

  use io
  use ExtFEcomparer_typedefs

  implicit none

  private

  public :: ExtFEcomparer_load_vector


contains

subroutine ExtFEcomparer_load_vector(rproblem)

  ! We could just start with the load-command, but we don't
  ! The reason is simple: While it is very easy to load a block-vector,
  ! some output we want to read in might not be a block-vector
  ! So we first create the vector according to the discretisation

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout) :: rproblem
!</inputoutput>

!<local variables>
  integer :: NLMAX
  integer :: N, iunit,NVAR,i,j, NEst
  real(DP), dimension(:), pointer :: p_Ddata => NULL()
  character(LEN=ExtFE_STRLEN) :: sVectorName
  character :: sVarname
  real(DP) :: dval

  logical :: lFormatted, lUnsorted
!</local variables>

  NLMAX = rproblem%NLMAX

  ! Create the vector according to the discretisation
  ! The benefit of this approach is that all flags in the
  ! vectors are set correct. (especially attaching the
  ! discretisation to the vector can sometimes become tricky).
  ! Another big plus is that for the second input format we
  ! just can read the coefficents in the data array
  ! and are fine then
  call lsysbl_createVecBlockByDiscr(rproblem%rdiscretisation,&
            rproblem%coeffVector,bclear=.TRUE.,cdataType=ST_DOUBLE )

  ! Now we have a vector with the right size and a block discretisation
  ! so if we fill it with data now we have a FE-Function
  ! If the vector-file actually is a block-vector it is easy, else
  ! it is a bit more tricky.


  ! What type of Vector? Real Block-Vector or a postprocessed?

  select case (rproblem%vectorType)

    ! First of all, the easy case: a real block vector
    case (ExtFE_isSolutionBlockvector)
        ! first we need to know what we get in the vector file
        ! Formatted, unformatted, sorted, unsorted?
        select case (rproblem%vectorFileFormat)
            case (ExtFE_formatted_unsorted)
            lFormatted = .TRUE.
            lUnsorted = .FALSE.
        case (ExtFE_unformatted_unsorted)
            lFormatted = .FALSE.
            lUnsorted = .FALSE.
        case (ExtFE_formatted_sorted)
            lFormatted = .TRUE.
            lUnsorted = .TRUE.
        case (ExtFE_unformatted_sorted)
            lFormatted = .FALSE.
            lUnsorted = .TRUE.
        case default
            call output_line("Your vector file format is not supported",&
                OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_load_vector: vector format")
            call sys_halt()
        end select ! file format

        ! Now we can load it directly
        call vecio_readBlockVectorHR(rproblem%coeffVector,sVectorName,lUnsorted,0, &
                                     rproblem%sVectorFile,lFormatted)

    case(ExtFE_isPostprocessed)
        ! We need to know if it is formatted or unformatted
        ! At the moment we only support unformatted
        select case (rproblem%vectorFileFormat)
            case (ExtFE_unformatted)
                ! We have an vector written out in an unformatted
                ! (that "non-human-readable") file.
                ! The format specification is:
                ! First an Integer N, then a character and then N doubles
                ! The integer N is not used for allocation since we know
                ! from the discretisation how many entries the vector must have
                lFormatted = .FALSE.

                ! Read out how many variables are in the vector
                ! 1 variable is something like velocity_x
                NVAR = rproblem%NVAR

                ! Open the file
                call io_openFileForReading(rproblem%sVectorFile,iunit,lFormatted)
                ! Loop over all variables
                do i=1,NVAR
                    NEst = rproblem%coeffVector%RvectorBlock(i)%NEQ
                    ! We don't really need the name of the variable
                    read(iunit) N, sVarname
                    if(NEst .lt. N) then
                        call output_line("According to the discretisation your vector should have length "&
                         //trim(sys_siL(NEst,20)) ,OU_CLASS_WARNING,OU_MODE_STD)
                        call output_line(" but according to you file it has length "//trim(sys_siL(N,20)),&
                            OU_CLASS_WARNING,OU_MODE_STD)
                        call output_line("To avoid a segmentation fault we are just going to read in the first "&
                            //trim(sys_siL(NEst,20)), OU_CLASS_WARNING,OU_MODE_STD)
                        call output_line("values of this vector. However we continue since this is left open as small",&
                            OU_CLASS_WARNING,OU_MODE_STD)
                        call output_line("hack if you want to analyze only one function",OU_CLASS_WARNING,OU_MODE_STD)
                        call output_line("Pay attention to your results as now some strings are interpreted as doubles",&
                            OU_CLASS_WARNING,OU_MODE_STD)
                    elseif(NEst .gt. N) then
                        call output_line("According to the discretisation you vector should have length "&
                         //trim(sys_siL(NEst,20)) ,OU_CLASS_ERROR,OU_MODE_STD)
                        call output_line(" but according to you file it has length "//trim(sys_siL(N,20)),&
                            OU_CLASS_ERROR,OU_MODE_STD)
                        call output_line(" Since we cannot manipulate to avoid segmentation faults we quit here", &
                            OU_CLASS_ERROR,OU_MODE_STD)
                        !Close the file
                        close(iunit)
                        call sys_halt()
                    end if
                    call lsyssc_getbase_double(rproblem%coeffVector%RvectorBlock(i),p_Ddata)
                       do j=1,NEst
                          read(iunit) dval
                          p_Ddata(j) = dval
                       end do
                end do

                close(iunit)

            case default
                call output_line("Wrong input for vector format",&
                    OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_load_vector: vector format")
                call sys_halt()
        end select

    case default
        call output_line("Wrong input for the vector type",&
             OU_CLASS_ERROR,OU_MODE_STD,"ExtFEcomparer_load_vector: vector type")
        call sys_halt()
  end select


end subroutine




end module
