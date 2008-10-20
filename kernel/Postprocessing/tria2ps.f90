!##############################################################################
!# ****************************************************************************
!# <name> tria2ps </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains some routines which export a 2D triangulation into
!# a postscript file.
!#
!# </purpose>
!#
!##############################################################################
module tria2ps

use fsystem
use genoutput
use storage
use io
use triangulation

implicit none

!<constants>

!<constantblock description="Flags for the PostScript exporter">
  ! If set, then the aspect ratio of the mesh is kept. If not set, then
  ! the mesh is stretched such that it fills the draw box.
  integer(I32), parameter :: TRIA2PS_FLAG_KEEP_ASPECT_RATIO = 2**0

!</constantblock>

!<constantblock description="Alignment identifiers">

  ! If neither TRIA2PS_ALIGN_TOP nor TRIA2PS_ALIGN_BOTTOM is set, then
  ! the mesh is centered vertically in the drawing area - dito for
  ! horizontally centered.

  ! The mesh is aligned to the top border of the drawing area.
  !integer(I32), parameter :: TRIA2PS_ALIGN_TOP              = 2**0

  ! The mesh is aligned to the bottom border of the drawing area.
  !integer(I32), parameter :: TRIA2PS_ALIGN_BOTTOM           = 2**1
  
  ! The mesh is aligned to the left border of the drawing area.
  ! Cannot be used together with TRIA2PS_ALIGN_RIGHT.
  !integer(I32), parameter :: TRIA2PS_ALIGN_LEFT             = 2**2
  
  ! The mesh is aligned to the right border of the drawing area.
  ! Cannot be used together with TRIA2PS_ALIGN_LEFT.
  !integer(I32), parameter :: TRIA2PS_ALIGN_RIGHT            = 2**3

!</constantblock>

!<types>

!<typeblock description="Configuration structure for the tria2ps routines">

  type t_psconfig
  
    ! Flags for the exporter. A combination of TRIA2PS_FLAG_XXXX constants
    ! defined above.
    integer(I32) :: cflags = TRIA2PS_FLAG_KEEP_ASPECT_RATIO
    
    ! Specifies the alignment of the mesh in the drawing area.
    ! A valid combination of TRIA2PS_ALIGN_XXXX constants defined above.
    ! By default, the mesh is centered in the drawing area.
    !integer(I32) :: calign = 0
  
    ! Size of the PostScript drawing area in millimeters.
    real(DP) :: dwidth  = 100.0_DP
    real(DP) :: dheight = 100.0_DP
    
    ! Transformation matrix. Each vertice is multiplied by this matrix before
    ! it is rendered into the postscript file. If you e.g. want to draw a
    ! rotated mesh, then Dtrafo would be the corresponding rotation matrix.
    real(DP), dimension(2,2) :: Dtrafo = (/(/1_DP,0_DP/),(/0_DP,1_DP/)/)
    
    ! Line width that is to be used for drawing the edges.
    real(DP) :: dedgeWidth = 0.1_DP

  end type

!</typeblock>

!</types>


contains

  ! ***************************************************************************

!<subroutine>

  subroutine tria2ps_writePostScript(rtria, sfilename, rconfig)

!<description>
  ! Exports a mesh into a PostScript file.
!</description>

!<input>
  ! The mesh that is to be exported. Must be a valid 2D mesh.
  type(t_triangulation), intent(IN) :: rtria
  
  ! The filename of the PostScript file that is to be generated.
  character(len=*), intent(IN) :: sfilename
  
  ! OPTIONAL: A configuation structure which specifies more information
  ! for the exporter.
  type(t_psconfig), optional, intent(IN) :: rconfig
!</input>

!</subroutine>

  ! Local variables
  type(t_psconfig) :: rcfg
  real(DP) :: dtmp
  real(DP), dimension(2) :: Dv
  real(DP), dimension(:,:), pointer :: p_Dcoords
  integer, dimension(:,:), pointer :: p_Iedges
  integer :: i,j,NVT, NMT
  integer :: iunit
  
  ! Bounding box of mesh
  real(DP) :: dbboxMinX, dbboxMinY, dbboxMaxX, dbboxMaxY,dbboxWidth,dbboxHeight
  
  ! Scaling parameters
  real(DP) :: dscaleX, dscaleY
  
    ! Did the user specify a config structure?
    if(present(rconfig)) then
      ! Okay, so let's check whether the configuration is valid.
      
      ! First of all, let's make sure the drawing area is not empty.
      if((rconfig%dwidth .le. SYS_EPSREAL) .or. &
         (rconfig%dheight .le. SYS_EPSREAL)) then
        call output_line('Drawing area is invalid!', OU_CLASS_ERROR,&
                         OU_MODE_STD, 'tria2ps_writePostScript')
        call sys_halt()
      end if
      
      ! And let's make sure the transformation matrix is regular - otherwise
      ! we would divide by zero later!
      dtmp = rconfig%Dtrafo(1,1)*rconfig%Dtrafo(2,2) &
           - rconfig%Dtrafo(1,2)*rconfig%Dtrafo(2,1)
      if(abs(dtmp) .le. SYS_EPSREAL) then
        call output_line('Transformation matrix is singular!', OU_CLASS_ERROR,&
                         OU_MODE_STD, 'tria2ps_writePostScript')
        call sys_halt()
      end if
      
      ! And make sure the line width is positive.
      if(rconfig%dedgeWidth .le. SYS_EPSREAL) then
        call output_line('Edge width must be positive!', OU_CLASS_ERROR,&
                         OU_MODE_STD, 'tria2ps_writePostScript')
        call sys_halt()
      end if
      
      ! Okay, the configuration seems to be correct.
      rcfg = rconfig
      
    end if
    
    ! Now make sure the triangulation is a 2D mesh.
    if(rtria%ndim .ne. NDIM2D) then
      call output_line('Only 2D triangulations supported!', OU_CLASS_ERROR,&
                       OU_MODE_STD, 'tria2ps_writePostScript')
      call sys_halt()
    end if
    
    ! Okay, get the necessary information from the mesh.
    NVT = rtria%NVT
    NMT = rtria%NMT
    call storage_getbase_double2D(rtria%h_DvertexCoords, p_Dcoords)
    call storage_getbase_int2D(rtria%h_IverticesAtEdge, p_Iedges)
    
    ! We now need to calculate the bounding box of the mesh.
    ! So get the first vertice and initialise the bounding box with it.
    Dv = matmul(rcfg%Dtrafo, p_Dcoords(1:2,1))
    dbboxMinX = Dv(1)
    dbboxMaxX = Dv(1)
    dbboxMinY = Dv(2)
    dbboxMaxY = Dv(2)
    
    ! And loop through the rest of the vertices.
    do i = 2, NVT
      
      ! Calculate transformed vertice
      Dv = matmul(rcfg%Dtrafo, p_Dcoords(1:2,i))
      
      ! And update the bounding box.
      dbboxMinX = min(dbboxMinX, Dv(1))
      dbboxMaxX = max(dbboxMaxX, Dv(1))
      dbboxMinY = min(dbboxMinY, Dv(2))
      dbboxMaxY = max(dbboxMaxY, Dv(2))
      
    end do
    
    ! Calculate the dimensions of the bounding box
    dbboxWidth  = dbboxMaxX - dbboxMinX
    dbboxHeight = dbboxMaxY - dbboxMinY
    
    ! Now comes the interesting part: We need to calculate the transformation
    ! from mesh vertice coordinates into the coordinates of the PostScript file.
    
    ! Calculate the scaling parameters:
    dscaleX = rcfg%dwidth / dbboxWidth
    dscaleY = rcfg%dheight / dbboxHeight
    
    ! Do we have to keep the aspect ratio?
    if(iand(rcfg%cflags,TRIA2PS_FLAG_KEEP_ASPECT_RATIO) .ne. 0) then
      
      ! Yes, so choose the minimum scaling factor.
      dscaleX = min(dscaleX, dscaleY)
      dscaleY = dscaleX
      
    end if
    
    ! Okay, open a file for writing
    call io_openFileForWriting(sfilename,iunit,SYS_REPLACE,bformatted=.true.)
    
    ! Fail?
    if(iunit .le. 0) then
      call output_line('Failed to open file for writing!', OU_CLASS_ERROR,&
                       OU_MODE_STD, 'tria2ps_writePostScript')
      call sys_halt()
    end if
    
    ! Okay, write PostScript header
    write(iunit,'(A)') '%!!PS-Adobe-3.0 EPSF-3.0'
    
    ! Write the line width
    write(iunit,'(F12.6)') rcfg%dedgeWidth
    
    ! Now go through all edges
    do i = 1, NMT
    
      ! Get the first vertice
      Dv = matmul(rcfg%Dtrafo,p_Dcoords(1:2,p_Iedges(1,i)))
      
      ! Transform the coordinates
      Dv(1) = (Dv(1) - dbboxMinX) * dscaleX
      Dv(2) = (Dv(2) - dbboxMinY) * dscaleY
      
      ! Write first vertice
      write(iunit,'("newpath ",F12.6,F12.6," moveto")') Dv(1), Dv(2)
      
      ! Get the second vertice
      Dv = matmul(rcfg%Dtrafo,p_Dcoords(1:2,p_Iedges(2,i)))
      
      ! Transform the coordinates
      Dv(1) = (Dv(1) - dbboxMinX) * dscaleX
      Dv(2) = (Dv(2) - dbboxMinY) * dscaleY
      
      write(iunit,'(F12.6,F12.6," lineto")') Dv(1), Dv(2)
      write(iunit,'(A)') 'stroke'
    
    end do
    
    ! Close the file
    close(iunit)
  
  end subroutine

end module