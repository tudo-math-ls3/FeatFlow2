module meshadaptbase

  use boundary
  use fsystem
  use genoutput
  use hadaptaux
  use hadaptivity
  use io
  use linearsystemscalar
  use signals
  use storage
  use triangulation

  implicit none
  
  private
  public :: perform_adaptation
  
  ! global variables
  type(t_boundary), save, public :: rbnd
  type(t_hadapt), save, public :: rhadapt
  type(t_triangulation), save, public :: rtria

  character(len=256), save, public :: smesh,serror,sdata
  integer, save, public :: ilinelen
  logical, save, public :: bbnd

contains
  
  function perform_adaptation(isignum) result(iresult)
    
    integer, intent(in) :: isignum
    integer :: iresult

    ! local variables
    type(t_vectorScalar) :: rindicator
    real(DP), dimension(:), pointer :: p_Dindicator
    integer :: iel,ios,iunit
       
    select case(isignum)
    case (SIGINT)
      ! Read indicator vector
      call output_line("Reading indicator field from './"//trim(serror)//"'...")
      call io_openFileForReading('./'//trim(serror), iunit, .true.)
      
      ! Read first line from file
      read(iunit, fmt=*) iel
      if (iel .ne. rtria%NEL) then
        call output_line ("Mismatch in number of elements!", &
                          OU_CLASS_ERROR,OU_MODE_STD,"meshadapt")
        call sys_halt()
      end if
      
      ! Create indicator
      call lsyssc_createVector(rindicator, rtria%NEL, .true., ST_DOUBLE)
      call lsyssc_getbase_double(rindicator, p_Dindicator)
      
      do iel=1,rtria%NEL
        call io_readlinefromfile(iunit, sdata, ilinelen, ios)
        p_Dindicator(iel) = sys_str2Double(sdata, "(F20.10)")
      end do
      close(iunit)
      
      ! Perform mesh adaptation
      call hadapt_refreshAdaptation(rhadapt, rtria)
      call hadapt_performAdaptation(rhadapt, rindicator)
      
      ! Release indicator
      call lsyssc_releaseVector(rindicator)
      
      ! Update triangulation structure
      call hadapt_generateRawMesh(rhadapt, rtria)
      
      ! Initialise standard mesh
      if(bbnd) then
        call tria_initStandardMeshFromRaw (rtria, rbnd)
      else
        call tria_initStandardMeshFromRaw (rtria)
      end if
      
      ! Export triangulation structure
      call output_line("Exporting triangulation to './"//trim(smesh)//"_ref.tri'...")
      if (bbnd) then
        call tria_exportTriFile(rtria, './'//trim(smesh)//'_ref.tri', TRI_FMT_STANDARD)
      else
        call tria_exportTriFile(rtria, './'//trim(smesh)//'_ref.tri', TRI_FMT_NOPARAMETRISATION)
      end if
      iresult = 0

    case (SIGQUIT)

      ! Clean up
      call hadapt_releaseAdaptation(rhadapt)
      call tria_done(rtria)
      if(bbnd) call boundary_release(rbnd)
      
      ! Clean up the storage management, finish
      call storage_done()
      iresult = 0

      stop
      
    case default
       iresult = 0

    end select

  end function perform_adaptation

end module meshadaptbase
