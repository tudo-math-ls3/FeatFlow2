!##############################################################################
!# ****************************************************************************
!# <name> meshadapt </name>
!# ****************************************************************************
!#
!# <purpose>
!# This application performs local mesh adaptation based on a given
!# indicator function and prescribed refinement/recoarsening tolerances.
!#
!# Calling example:
!#
!# meshadapt -read2d mymesh -indicator myindicator.dat -refmax 3 \
!#           -reftol 0.8 -crstol 0.2
!#
!# The initial 2D mesh is read from the TRI/PRM file mymesh.tri/prm
!# and refined based on the error element-wise indicator given in file
!# myerror.dat. Cells with error larger than 0.8 are refined and cells
!# with error smalled than 0.2 are re-coarsened. The maximum
!# refinement level is 3, that is, an element of the initial grid can
!# be refined 3 times at most.
!#
!# </purpose>
!##############################################################################

program meshadapt

  use fsystem
  use genoutput
  use signals
  use storage
  
  use meshadaptbase

  implicit none
  
  ! Initialise system-wide settings
  call sys_init()
  
  ! Initialise the output system
  call output_init()

  ! Initialise the FEAT 2.0 storage management
  call storage_init(100, 100)
  
  if (madapt_signalhandler(SIGUSR1) .eq. 0) then
    
    ! We are in daemon mode, hence, register signal handler
    call fsignal(SIGUSR2, madapt_signalhandler) ! export mesh
    call fsignal(SIGINT,  madapt_signalhandler) ! single step
    call fsignal(SIGQUIT, madapt_signalhandler) ! finalise
    call fsignal(SIGHUP,  madapt_signalhandler) ! finalise
    call fsignal(SIGTERM, madapt_signalhandler) ! finalise

    daemon: do
      ! Perform mesh adaptation step?
      call fsignal(SIGINT, madapt_signalhandler)

      ! Export mesh to file?
      call fsignal(SIGUSR2, madapt_signalhandler)

      ! Finalise and exit?
      call fsignal(SIGQUIT, madapt_signalhandler)
      call fsignal(SIGTERM, madapt_signalhandler)
      call fsignal(SIGHUP,  madapt_signalhandler)
    end do daemon

  else
    
    ! We are not in daemon mode, hence, perform a single mesh
    ! adaptation step, export the mesh to file and exit
    if (madapt_signalhandler(SIGINT) .ne. 0) then
      call output_line("An error occured during mesh adaptation!",&
          OU_CLASS_ERROR,OU_MODE_STD,"meshadapt")
      call sys_halt()
    end if

    if (madapt_signalhandler(SIGUSR2) .ne. 0) then
      call output_line("An error occured while exporting mesh to file!",&
          OU_CLASS_ERROR,OU_MODE_STD,"meshadapt")
      call sys_halt()
    end if

    if (madapt_signalhandler(SIGQUIT) .ne. 0) then
      call output_line("An error occured during finalisation!",&
          OU_CLASS_ERROR,OU_MODE_STD,"meshadapt")
      call sys_halt()
    end if
  end if
  
  ! Clean up the storage management, finish
1 call storage_done()

contains

  subroutine meshadaptQuick(ndim,smesh,nref,nrefmax,dreftol,dcrstol)

    integer, intent(in) :: ndim,nref,nrefmax
    character(len=*), intent(in) :: smesh
    real(DP), intent(in) :: dreftol,dcrstol

    ! local variables
    type(t_meshAdapt) :: rmeshAdapt
    real(DP), dimension(:,:), allocatable :: Dcoords
    real(DP), dimension(:), allocatable :: Dind
    integer, dimension(:,:), allocatable :: IverticesAtElement
    real(DP) :: xc,yc,phi,W
    integer :: nel,nvt,nnve,iref,iel,ive,istatus
    
    ! Initialisation
    call madapt_init(rmeshAdapt,ndim,smesh)

    ! Get data from adaptation structure
    nel = madapt_getnel(rmeshAdapt)
    nvt = madapt_getnvt(rmeshAdapt)
    nnve = madapt_getnnve(rmeshAdapt)

    ! Get mesh from adaptation structure
    allocate(Dcoords(ndim,nvt), IverticesAtElement(nnve,nel))
    call madapt_getVertexCoords(rmeshAdapt, Dcoords)
    call madapt_getVerticesAtElement(rmeshAdapt, IverticesAtElement)
    
    W = 0.0_DP
    do iel=1,nel
      do ive=1,nnve
        W = W + sqrt((Dcoords(1,IverticesAtElement(ive,iel))-&
                      Dcoords(1,IverticesAtElement(mod(ive,nnve)+1,iel)))**2&
              +      (Dcoords(2,IverticesAtElement(ive,iel))-&
                      Dcoords(2,IverticesAtElement(mod(ive,nnve)+1,iel)))**2)       
      end do
    end do
    W = W/real(nel*nnve)

    ! Deallocate local data
    deallocate(Dcoords,IverticesAtElement)

    ! Perform mesh refinement
    do iref=1,nref
      
      ! Get mesh from adaptation structure
      allocate(Dcoords(ndim,nvt), IverticesAtElement(nnve,nel),Dind(nel))
      call madapt_getVertexCoords(rmeshAdapt, Dcoords)
      call madapt_getVerticesAtElement(rmeshAdapt, IverticesAtElement)
    
      do iel=1,nel
        ! Create circular indicator function
        xc = sum(Dcoords(1,IverticesAtElement(:,iel)))/real(nnve,DP)
        yc = sum(Dcoords(2,IverticesAtElement(:,iel)))/real(nnve,DP)
        phi = 0.2_DP - sqrt( (xc-0.5_DP)**2 + (yc-0.5_DP)**2 )

        ! Define indicator array
        if (abs(phi) .le. W) then
          Dind(iel) = 2.0_DP
        else
          Dind(iel) = 1.0_DP
        end if
      end do

      ! Perform one step of mesh adaptation
      istatus = madapt_step(rmeshAdapt,nel,Dind,nrefmax,dreftol,dcrstol)

      ! Deallocate local data
      deallocate(Dcoords,IverticesAtElement,Dind)

      ! Get data from adaptation structure
      nel = madapt_getnel(rmeshAdapt)
      nvt = madapt_getnvt(rmeshAdapt)
      nnve = madapt_getnnve(rmeshAdapt)

      print *, "REFINE:",nel,nvt

      W = W/2.0_DP
    end do
    
    ! Perform mesh re-coarsening
    do iref=1,nref

      ! Define indicator array
      allocate(Dind(nel)); Dind = 0.1
      
      ! Perform one step of mesh adaptation
      istatus = madapt_step(rmeshAdapt,nel,Dind,nrefmax,dreftol,dcrstol)
      
      ! Deallocate local data
      deallocate(Dind)

      ! Get data from adaptation structure
      nel = madapt_getnel(rmeshAdapt)
      nvt = madapt_getnvt(rmeshAdapt)

      print *, "COARSEN:",nel,nvt

      ! Get mesh from adaptation structure
      allocate(Dcoords(ndim,nvt), IverticesAtElement(nnve,nel))

      call madapt_getVertexCoords(rmeshAdapt, Dcoords)
      call madapt_getVerticesAtElement(rmeshAdapt, IverticesAtElement)


      ! Deallocate local data
      deallocate(Dcoords,IverticesAtElement)
    end do

    
    
    ! Finalisation
    call madapt_done(rmeshAdapt)

  end subroutine meshadaptQuick

end program meshadapt
