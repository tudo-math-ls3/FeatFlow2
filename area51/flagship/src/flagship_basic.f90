!##############################################################################
!# ****************************************************************************
!# <name> flagship_basic </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the basic routines for the main program
!#
!# The following routines are available:
!#
!# 1.) flagship_readParserFromFile
!#     -> Reads expressions from file and initialise the function parser
!#
!# 2.) flagship_initUCDexport
!#     -> Initialises the UCD exporter for file output
!#
!# 3.) flagship_updateSolverMatrix
!#     -> Updates the matrix in the solver structure
!#
!# 4.) flagship_initPerfConfig
!#     -> Initialises performance configurations from file
!# </purpose>
!##############################################################################

module flagship_basic

#include "flagship.h"

!$use omp_lib
  use fsystem
  use fparser
  use genoutput
  use io
  use perfconfig
  use problem
  use solveraux
  use triangulation
  use ucd

  use linearsystemscalar, only     : lsyssc_initPerfConfig
  use pprocerror, only             : pperr_initPerfConfig
  use pprocgradients, only         : ppgrd_initPerfConfig
  use afcstabbase, only            : afcstab_initPerfConfig
  use afcstabscalar, only          : afcsc_initPerfConfig
  use afcstabsystem, only          : afcsys_initPerfConfig
  use groupfemscalar, only         : gfsc_initPerfConfig
  use groupfemsystem, only         : gfsys_initPerfConfig
  use elementpreprocessing, only   : elprep_initPerfConfig
  use transformation, only         : trafo_initPerfConfig
  use linearalgebra, only          : lalg_initPerfConfig
  use bilinearformevaluation, only : bilf_initPerfConfig
  use linearformevaluation, only   : linf_initPerfConfig

  implicit none

  private
  public :: flagship_readParserFromFile
  public :: flagship_initUCDexport
  public :: flagship_updateSolverMatrix
  public :: flagship_initPerfConfig

!<constants>

!<constantblock description="Flags for matrix update specification bitfield">

  ! Update the matrix only for the current solver
  integer(I32), parameter, public :: UPDMAT_NORECURSIVE            = 2**0

  ! Update the matrix for coarse grid solver
  integer(I32), parameter, public :: UPDMAT_LINEAR_SOLVER          = 2**1

  ! Update the preconditioner for the coarse grid solver
  integer(I32), parameter, public :: UPDMAT_LINEAR_PRECOND         = 2**2

  ! Update the solver of the linear multigrid solver
  integer(I32), parameter, public :: UPDMAT_LINEARMG_SOLVER        = 2**3

  ! Update the smoothers of the linear multigrid solver
  integer(I32), parameter, public :: UPDMAT_LINEARMG_SMOOTHER      = 2**4

  ! Update the coarsegrid solver of the linear multigrid solver
  integer(I32), parameter, public :: UPDMAT_LINEARMG_COARSEGRID    = 2**5

!</constantblock>

!<constantblock description="Predified bitfields for matrix update">

  ! Update everything
  integer(I32), parameter, public :: UPDMAT_ALL = UPDMAT_LINEAR_SOLVER +&
                                                  UPDMAT_LINEAR_PRECOND +&
                                                  UPDMAT_LINEARMG_SOLVER +&
                                                  UPDMAT_LINEARMG_SMOOTHER +&
                                                  UPDMAT_LINEARMG_COARSEGRID

  ! Update Jacobian matrix for steady-state flows
  integer(I32), parameter, public :: UPDMAT_JAC_STEADY = UPDMAT_LINEAR_SOLVER +&
                                                         UPDMAT_LINEARMG_SOLVER +&
                                                         UPDMAT_LINEARMG_SMOOTHER

  ! Update Jacobian matrix for transient flows
  integer(I32), parameter, public :: UPDMAT_JAC_TRANSIENT = UPDMAT_LINEAR_SOLVER +&
                                                            UPDMAT_LINEAR_PRECOND +&
                                                            UPDMAT_LINEARMG_SOLVER +&
                                                            UPDMAT_LINEARMG_SMOOTHER

!</constantblock>


!<constantblock description="Global type of system structure">

  ! time-dependent flow
  integer, parameter, public :: SYSTEM_INTERLEAVEFORMAT        = 0

  ! steady-state flow
  integer, parameter, public :: SYSTEM_BLOCKFORMAT             = 1

!</constantblock>

!</constants>


contains

  ! ***************************************************************************

!<subroutine>

  subroutine flagship_readParserFromFile(sfilename, ssectionname,&
                                         cvariables, rparser)

!<description>
    ! This subroutine initialises a vector profile from the parameters
    ! specified in the parameter file by calling a function parser w.r.t.
    ! the specified variable names
!</description>

!<input>
    ! name of parameter file
    character(LEN=*), intent(in) :: sfilename

    ! name of the parameter section
    character(LEN=*), intent(in) :: ssectionname

    ! symbolic variable names
    character(LEN=*), dimension(:), intent(in) :: cvariables
!</input>

!<output>
    ! function parser
    type(t_fparser), intent(out) :: rparser
!</output>
!</subroutine>

    ! local variables
    character(SYS_STRLEN) :: skeyword
    character(1024) :: sdata,sexpression
    integer :: iunit,ios,ipos,idatalen,icomp,ncomp

    ! Try to open the file
    call io_openFileForReading(sfilename, iunit, .true.)

    ! Oops...
    if (iunit .eq. -1) then
      call output_line('Unable to open input file!',&
                       OU_CLASS_WARNING,OU_MODE_STD,'flagship_readParserFromFile')
      call sys_halt()
    end if

    ! Read through the input file until the given keyword is reached.
    ios = 0
    do while(ios .eq. 0)

      ! Read next line in file
      call io_readlinefromfile(iunit, sdata, idatalen, ios)
      if (ios .ne. 0) then
        call output_line('Unable to read KEYWORD from input file!',&
                         OU_CLASS_WARNING,OU_MODE_STD,'flagship_readParserFromFile')
        call sys_halt()
      end if

      ! Check for keyword
      call sys_tolower(sdata(1:idatalen), skeyword)
      if (trim(adjustl(skeyword)) .eq. trim(adjustl(ssectionname))) exit
    end do

    ! We found the keyword. String NCOMP must be the next line to read.
    call io_readlinefromfile(iunit, sdata, idatalen, ios)
    if (ios .ne. 0) then
      call output_line('Unable to read data from input file!',&
                       OU_CLASS_WARNING,OU_MODE_STD,'flagship_readParserFromFile')
      call sys_halt()
    end if

    ! Check for keyword NCOMP
    call sys_tolower(sdata(1:idatalen), skeyword)
    if (trim(adjustl(skeyword)) .ne. 'ncomp') then
      call output_line('Syntax error in input file!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'flagship_readParserFromFile')
      call sys_halt()
    end if

    ! Read value for NCOMP
    read(iunit,*,IOSTAT=ios) ncomp
    if (ios .ne. 0) then
      call output_line('Unable to read value of NCOMP from input file!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'flagship_readParserFromFile')
      call sys_halt()
    end if

    ! Create parser structure
    call fparser_create(rparser, ncomp)

    ! Read expressions into parser, relations and concatinations
    do icomp = 1, ncomp

      ! Read until expression is finished
      ios  = 0
      ipos = 1
      sexpression = " "

      do while(ios .eq. 0)

        ! Read next line in file
        call io_readlinefromfile(iunit, sdata, idatalen, ios)
        if (ios .ne. 0) then
          call output_line('Syntax error in input file!',&
                           OU_CLASS_ERROR,OU_MODE_STD,'flagship_readParserFromFile')
          call sys_halt()
        end if

        ! Append line to expression
        sexpression(ipos:) = sdata(1:idatalen)
        ipos = len_trim(sexpression)

        ! Check if expression is continued in the following line
        if (sexpression(max(1,ipos-2):ipos) .eq. '...') then
          ipos = ipos-2
        else
          exit
        end if
      end do

      call fparser_parseFunction(rparser, icomp, sexpression(1:ipos), cvariables)
    end do

    ! We are done. Close the unit.
    close(iunit)

  end subroutine flagship_readParserFromFile

  !*****************************************************************************

!<subroutine>

  subroutine flagship_initUCDexport(rproblemLevel, sfilename, ioutputUCD,&
                                    rexport, ifilenumber, rtriangulation)

!<description>
    ! This subroutine initialises the UCD exporter structure. If the
    ! optional parameter ifilenumber is given, the outputfile is named
    ! 'sfilename'.<ifilenumber>.'ext' where 'ext' is the file
    ! extension that corresponds to the UCD format.
!</description>

!<input>
    ! multigrid structure
    type(t_problemLevel), intent(in), target :: rproblemLevel

    ! name of the output file
    character(LEN=*), intent(in) :: sfilename

    ! type of UCD output
    integer, intent(in) :: ioutputUCD

    ! OPTIONAL: number of the output file
    integer, intent(in), optional :: ifilenumber

    ! OPTIONAL: triangulation which overwrites the triangulation
    ! associated with the problem level structure
    type(t_triangulation), target, optional :: rtriangulation
!</input>

!<output>
    ! UCD export structure
    type(t_ucdExport), intent(out) :: rexport
!</output>
!</subroutine>

    ! local variables
    type(t_triangulation), pointer :: p_rtriangulation

    ! Set pointer to triangulation
    if (present(rtriangulation)) then
      p_rtriangulation => rtriangulation
    else
      p_rtriangulation => rproblemLevel%rtriangulation
    end if

    select case(ioutputUCD)
    case (UCD_FORMAT_GMV)
      if (present(ifilenumber)) then
        call ucd_startGMV(rexport, UCD_FLAG_STANDARD, p_rtriangulation,&
            trim(adjustl(sfilename))//'.'//trim(sys_si0(ifilenumber,5))//'.gmv')
      else
        call ucd_startGMV(rexport, UCD_FLAG_STANDARD,&
            p_rtriangulation, trim(adjustl(sfilename))//'.gmv')
      end if

    case (UCD_FORMAT_BGMV)
      if (present(ifilenumber)) then
        call ucd_startBGMV(rexport, UCD_FLAG_STANDARD, p_rtriangulation,&
            trim(adjustl(sfilename))//'.'//trim(sys_si0(ifilenumber,5))//'.gmv')
      else
        call ucd_startBGMV(rexport, UCD_FLAG_STANDARD,&
            p_rtriangulation, trim(adjustl(sfilename))//'.gmv')
      end if

    case (UCD_FORMAT_AVS)
      if (present(ifilenumber)) then
        call ucd_startAVS(rexport, UCD_FLAG_STANDARD, p_rtriangulation,&
            trim(adjustl(sfilename))//'.'//trim(sys_si0(ifilenumber,5))//'.avs')
      else
        call ucd_startAVS(rexport, UCD_FLAG_STANDARD,&
            p_rtriangulation, trim(adjustl(sfilename))//'.avs')
      end if

    case (UCD_FORMAT_VTK)
      if (present(ifilenumber)) then
        call ucd_startVTK(rexport, UCD_FLAG_STANDARD, p_rtriangulation,&
            trim(adjustl(sfilename))//'.'//trim(sys_si0(ifilenumber,5))//'.vtu')
      else
        call ucd_startVTK(rexport, UCD_FLAG_STANDARD,&
            p_rtriangulation, trim(adjustl(sfilename))//'.vtu')
      end if

    case default
      call output_line('Invalid UCD output type!', &
                       OU_CLASS_ERROR,OU_MODE_STD,'flagship_initUCDexport')
      call sys_halt()
    end select

  end subroutine flagship_initUCDexport

  !*****************************************************************************

!<subroutine>

  subroutine flagship_updateSolverMatrix(rproblemLevel, rsolver, imatrix,&
                                         isystemformat, iupdflag, nlminOpt, nlmaxOpt)

!<description>
    ! This subroutine updates the solver structure by setting the matrices.
    ! If the optional parameters NLMINOPT and NLMAXOPT are not given, then
    ! only the current level of the given multigrid structure is processed.
!</description>

!<input>
    ! multigrid level to start with
    type(t_problemLevel), intent(in), target :: rproblemLevel

    ! number of the matrix in the problem level structure
    integer, intent(in) :: imatrix

    ! type of the system format
    integer, intent(in) :: isystemformat

    ! flags which is used to specify the matrices to be updated
    integer(I32), intent(in) :: iupdflag

    ! OPTIONAL: minimal multigrid level
    integer, intent(in), optional :: nlminOpt

    ! OPTIONAL: maximal multigrid level
    integer, intent(in), optional :: nlmaxOpt
!</input>

!<inputoutput>
    ! solver structure
    type(t_solver), intent(inout) :: rsolver
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: nlmin,nlmax

    ! Set minimal level
    if (present(nlminOpt)) then
      nlmin = nlminOpt
    else
      nlmin = solver_getMinimumMultigridlevel(rsolver, rproblemLevel%ilev)
    end if

    ! Set maximum level
    if (present(nlmaxOpt)) then
      nlmax = nlmaxOpt
    else
      nlmax = solver_getMaximumMultigridlevel(rsolver, rproblemLevel%ilev)
    end if

    ! Ok, let us update the solver (recursively?)
    call updateMatrix(rproblemLevel, rsolver, imatrix,&
                      isystemformat, iupdflag, nlmin, nlmax)

  contains

    ! Here, the real working routines follow.

    !**************************************************************
    ! Set the solver matrices recursively. The routine must be called
    ! with a valid solver structure RSOLVER which serves as top-level
    ! solver. Depending on the parameter IUPDFLAG, the matrices
    ! are only set for the given solver structure or recursively for
    ! all subnodes. The parameter IMATRIX denotes the matrix number.
    ! The parameters NLMIN/NLMAX denote the minimum/maximum level to
    ! be considered.

    recursive subroutine updateMatrix(rproblemLevel, rsolver, imatrix,&
                                      isystemformat, iupdflag, nlmin, nlmax)

      type(t_problemLevel), intent(in), target :: rproblemLevel
      type(t_solver), intent(inout) :: rsolver
      integer(I32), intent(in) :: iupdflag
      integer, intent(in) :: imatrix,isystemformat,nlmin,nlmax


      ! local variables
      type(t_problemLevel), pointer :: p_rproblemLevelTmp
      type(t_problemLevel), pointer :: p_rproblemLevelCoarse
      integer :: i

      ! What kind of solver are we?
      select case(rsolver%csolverType)
      case (SV_FMG)

        ! There are no matrices for the full multigrid solver.
        if (iand(iupdflag, UPDMAT_NORECURSIVE) .eq. UPDMAT_NORECURSIVE) return

        ! Directly proceed to the coarsegrid solver which serves as
        ! nonlinear solver on each level
        call updatematrix(rproblemLevel,&
                          rsolver%p_solverMultigrid%p_solverCoarsegrid,&
                          imatrix, isystemformat, iupdflag, nlmin, nlmax)


      case (SV_NONLINEARMG)

        ! There are no matrices for the nonlinear multigrid solver.
        if (iand(iupdflag, UPDMAT_NORECURSIVE) .eq. UPDMAT_NORECURSIVE) return

        ! Directly proceed to the coarsegrid solver
        call updateMatrix(rproblemLevel,&
                          rsolver%p_solverMultigrid%p_solverCoarsegrid,&
                          imatrix, isystemformat, iupdflag,&
                          rsolver%p_solverMultigrid%nlmin,&
                          rsolver%p_solverMultigrid%nlmin)

        ! Proceed to the smoothers
        if (associated(rsolver%p_solverMultigrid%p_smoother)) then
          do i = lbound(rsolver%p_solverMultigrid%p_smoother,1),&
                 ubound(rsolver%p_solverMultigrid%p_smoother,1)

            call updateMatrix(rproblemLevel,&
                              rsolver%p_solverMultigrid%p_smoother(i),&
                              imatrix, isystemformat, iupdflag,&
                              rsolver%p_solverMultigrid%nlmin,&
                              rsolver%p_solverMultigrid%nlmin)
          end do
        end if


      case (SV_NONLINEAR)

        ! There are no matrices for the nonlinear single-grid solver.
        if (iand(iupdflag, UPDMAT_NORECURSIVE) .eq. UPDMAT_NORECURSIVE) return

        ! Directly proceed to the linear solver subnode.
        call updateMatrix(rproblemLevel, rsolver%p_solverSubnode(1),&
                          imatrix, isystemformat, iupdflag, nlmin, nlmax)


      case (SV_LINEARMG)

        ! Are there multiple levels?
        if (rsolver%p_solverMultigrid%nlmin .eq.&
            rsolver%p_solverMultigrid%nlmax) then

          ! Proceed to single grid solver
          call updateMatrix(rproblemLevel,&
                            rsolver%p_solverMultigrid%p_solverCoarsegrid,&
                            imatrix, isystemformat, iupdflag,&
                            rsolver%p_solverMultigrid%nlmin,&
                            rsolver%p_solverMultigrid%nlmin)

        else

          ! We are on the level of the linear multigrid solver.
          p_rproblemLevelTmp    => rproblemLevel
          p_rproblemLevelCoarse => rproblemLevel
          do while(associated(p_rproblemLevelTmp))

            ! Do we have to set matrices for this level?
            if (p_rproblemLevelTmp%ilev > nlmax) then
              p_rproblemLevelTmp => p_rproblemLevelTmp%p_rproblemLevelCoarse
              cycle
            elseif(p_rproblemLevelTmp%ilev < nlmin) then
              exit
            end if

            ! What type of matrix format are we
            select case(isystemFormat)
            case (SYSTEM_INTERLEAVEFORMAT)

              if (iand(iupdflag, UPDMAT_LINEARMG_SOLVER) .eq.&
                                 UPDMAT_LINEARMG_SOLVER) then
                ! Set the system matrix for the linear solver
                call solver_setSolverMatrix(rsolver,&
                    p_rproblemLevelTmp%Rmatrix(imatrix),&
                    p_rproblemLevelTmp%ilev)
              end if

              if (iand(iupdflag, UPDMAT_LINEARMG_SMOOTHER) .eq.&
                                 UPDMAT_LINEARMG_SMOOTHER) then
                ! Set the system matrix for the linear smoother
                ! Note that the smoother is not required in the coarsest level
                if (p_rproblemLevelTmp%ilev > nlmin) then
                  call solver_setSmootherMatrix(rsolver,&
                      p_rproblemLevelTmp%Rmatrix(imatrix),&
                      p_rproblemLevelTmp%ilev)
                end if
              end if

            case (SYSTEM_BLOCKFORMAT)

              if (iand(iupdflag, UPDMAT_LINEARMG_SOLVER) .eq.&
                                 UPDMAT_LINEARMG_SOLVER) then
                ! Set the system matrix for the linear solver
                call solver_setSolverMatrix(rsolver,&
                    p_rproblemLevelTmp%RmatrixBlock(imatrix),&
                    p_rproblemLevelTmp%ilev)
              end if

              if (iand(iupdflag, UPDMAT_LINEARMG_SMOOTHER) .eq.&
                                 UPDMAT_LINEARMG_SMOOTHER) then
                ! Set the system matrix for the linear smoother
                ! Note that the smoother is not required in the coarsest level
                if (p_rproblemLevelTmp%ilev > nlmin) then
                  call solver_setSmootherMatrix(rsolver,&
                      p_rproblemLevelTmp%RmatrixBlock(imatrix),&
                      p_rproblemLevelTmp%ilev)
                end if
              end if

            case default
              call output_line('Unsupported system format!',&
                               OU_CLASS_ERROR,OU_MODE_STD,'updateMatrix')
              call sys_halt()
            end select

            ! Switch to next coarser level
            p_rproblemLevelCoarse => p_rproblemLevelTmp
            p_rproblemLevelTmp    => p_rproblemLevelTmp%p_rproblemLevelCoarse
          end do

          if (iand(iupdflag, UPDMAT_LINEARMG_COARSEGRID) .eq.&
                             UPDMAT_LINEARMG_COARSEGRID) then
            ! Set the system matrix for the linear coarse grid solver
            call updateMatrix(p_rproblemLevelCoarse,&
                              rsolver%p_solverMultigrid%p_solverCoarsegrid,&
                              imatrix, isystemformat, iupdflag,&
                              rsolver%p_solverMultigrid%nlmin,&
                              rsolver%p_solverMultigrid%nlmin)
          end if
        end if


      case (SV_LINEAR)

        ! The solver matrix and preconditioner matrix are only updated if
        ! the current level satisfies nlmin <= ilev <= nlmax
        if (nlmin .eq. rproblemLevel%ilev) then

          ! What type of matrix format are we
          select case(isystemFormat)
          case (SYSTEM_INTERLEAVEFORMAT)

            if (iand(iupdflag, UPDMAT_LINEAR_SOLVER) .eq.&
                               UPDMAT_LINEAR_SOLVER) then
              ! Set the system matrix of the single-grid solver
              call solver_setSolverMatrix(rsolver,&
                  rproblemLevel%Rmatrix(imatrix))
            end if

            if (iand(iupdflag, UPDMAT_LINEAR_PRECOND) .eq.&
                               UPDMAT_LINEAR_PRECOND) then
              ! Set the system matrix of the preconditioner
              call solver_setPrecondMatrix(rsolver,&
                  rproblemLevel%Rmatrix(imatrix))
            end if

          case (SYSTEM_BLOCKFORMAT)

            if (iand(iupdflag, UPDMAT_LINEAR_SOLVER) .eq.&
                               UPDMAT_LINEAR_SOLVER) then
              ! Set the system matrix of the single-grid solver
              call solver_setSolverMatrix(rsolver,&
                  rproblemLevel%RmatrixBlock(imatrix))
            end if

            if (iand(iupdflag, UPDMAT_LINEAR_PRECOND) .eq.&
                               UPDMAT_LINEAR_PRECOND) then
              ! Set the system matrix of the preconditioner
              call solver_setPrecondMatrix(rsolver,&
                  rproblemLevel%RmatrixBlock(imatrix))
            end if

          case default
            call output_line('Unsupported system format!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'updateMatrix')
            call sys_halt()
          end select

        else

          p_rproblemLevelTmp => rproblemLevel
          do while(associated(p_rproblemLevelTmp))

            ! Are we on the coarse grid level?
            if (p_rproblemLevelTmp%ilev > nlmin) then
              p_rproblemLevelTmp => p_rproblemLevelTmp%p_rproblemLevelCoarse
              cycle
            elseif(p_rproblemLevelTmp%ilev .eq. nlmin) then
              exit
            else
              call output_line('Invalid multigrid level!',&
                               OU_CLASS_ERROR,OU_MODE_STD,'updateMatrix')
              call sys_halt()
            end if
          end do

          ! What type of matrix format are we
          select case(isystemFormat)
          case (SYSTEM_INTERLEAVEFORMAT)

            if (iand(iupdflag, UPDMAT_LINEAR_SOLVER) .eq.&
                               UPDMAT_LINEAR_SOLVER) then
              ! Set the system matrix of the single-grid solver
              call solver_setSolverMatrix(rsolver,&
                  p_rproblemLevelTmp%Rmatrix(imatrix))
            end if

            if (iand(iupdflag, UPDMAT_LINEAR_PRECOND) .eq.&
                               UPDMAT_LINEAR_PRECOND) then
              ! Set the system matrix of the preconditioner
              call solver_setPrecondMatrix(rsolver,&
                  p_rproblemLevelTmp%Rmatrix(imatrix))
            end if

          case (SYSTEM_BLOCKFORMAT)

            if (iand(iupdflag, UPDMAT_LINEAR_SOLVER) .eq.&
                               UPDMAT_LINEAR_SOLVER) then
              ! Set the system matrix of the single-grid solver
              call solver_setSolverMatrix(rsolver,&
                  p_rproblemLevelTmp%RmatrixBlock(imatrix))
            end if

            if (iand(iupdflag, UPDMAT_LINEAR_PRECOND) .eq.&
                               UPDMAT_LINEAR_PRECOND) then
              ! Set the system matrix of the preconditioner
              call solver_setPrecondMatrix(rsolver,&
                  p_rproblemLevelTmp%RmatrixBlock(imatrix))
            end if

          case default
            call output_line('Unsupported system format!',&
                             OU_CLASS_ERROR,OU_MODE_STD,'updateMatrix')
            call sys_halt()
          end select
        end if

      case default
        call output_line('Unsupported solver type!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'updateMatrix')
        call sys_halt()
      end select
    end subroutine updateMatrix
  end subroutine flagship_updateSolverMatrix

  !*****************************************************************************

!<subroutine>

  subroutine flagship_initPerfConfig(sfilename, rperfconfig, ssectionName)

!<description>
    ! This subroutine initialises a performance configuration by the
    ! values stored in an INI-file. If the optional parameters
    ! rperfconfig and ssectionName are not present, then the global
    ! performance configurations of the kernerl modules are
    ! initialised. Otherwise, the given performance configuration
    ! rperfconfig is initialised by the values given in section
    ! ssectionName.
!</description>

!<input>
    ! file name
    character(len=*), intent(in) :: sfilename

    ! OPTIONAL: name of the section
    character(len=*), intent(in), optional :: ssectionName
!</input>

!<output>
    ! OPTIONAL: performance configuration
    type(t_perfconfig), intent(out), optional :: rperfconfig
!</output>
!</subroutine>

    type(t_perfconfig) :: rperfconfigTmp

    if (present(rperfconfig)) then
      
      ! Initialise the given performance configuration
      if (present(ssectionName)) then
        call pcfg_readPerfConfig(sfilename, ssectionName, rperfconfig)
      else
        call pcfg_readPerfConfig(sfilename, '', rperfconfig)
      end if

    else

      ! Initialise all global performance configurations
      call pcfg_readPerfConfig(sfilename, 'LSYSSC', rperfconfigTmp)
      call lsyssc_initPerfConfig(rperfConfigTmp)

      call pcfg_readPerfConfig(sfilename, 'PPERR', rperfconfigTmp)
      call pperr_initPerfConfig(rperfConfigTmp)

      call pcfg_readPerfConfig(sfilename, 'PPGRD', rperfconfigTmp)
      call ppgrd_initPerfConfig(rperfConfigTmp)

      call pcfg_readPerfConfig(sfilename, 'AFCSTAB', rperfconfigTmp)
      call afcstab_initPerfConfig(rperfConfigTmp)

      call pcfg_readPerfConfig(sfilename, 'AFCSC', rperfconfigTmp)
      call afcsc_initPerfConfig(rperfConfigTmp)

      call pcfg_readPerfConfig(sfilename, 'AFCSYS', rperfconfigTmp)
      call afcsys_initPerfConfig(rperfConfigTmp)

      call pcfg_readPerfConfig(sfilename, 'GFSC', rperfconfigTmp)
      call gfsc_initPerfConfig(rperfConfigTmp)

      call pcfg_readPerfConfig(sfilename, 'GFSYS', rperfconfigTmp)
      call gfsys_initPerfConfig(rperfConfigTmp)

      call pcfg_readPerfConfig(sfilename, 'EL', rperfconfigTmp)
      call elprep_initPerfConfig(rperfConfigTmp)

      call pcfg_readPerfConfig(sfilename, 'TRAFO', rperfconfigTmp)
      call trafo_initPerfConfig(rperfConfigTmp)

      call pcfg_readPerfConfig(sfilename, 'LALG', rperfconfigTmp)
      call lalg_initPerfConfig(rperfConfigTmp)

      call pcfg_readPerfConfig(sfilename, 'BILF', rperfconfigTmp)
      call bilf_initPerfConfig(rperfConfigTmp)

      call pcfg_readPerfConfig(sfilename, 'LINF', rperfconfigTmp)
      call linf_initPerfConfig(rperfConfigTmp)

      call pcfg_readPerfConfig(sfilename, 'FPAR', rperfconfigTmp)
      call fparser_initPerfConfig(rperfConfigTmp)

    end if

  end subroutine flagship_initPerfConfig

end module flagship_basic
