module domeasures

  use fsystem
  use triangulation

  implicit none

contains

!************************************************************************

!<subroutine>

  subroutine GetMeasure(Nel,Done,Dmeasure)

  !<description>
  ! This subroutine is used in the grid adaption for measuring the  
  ! domain area
  !</description>

  !<input>
  ! Number of elements
  integer, intent(IN) :: Nel
  ! Array of element areas
  real(DP), dimension(:), intent (IN) :: Done
  !</input>
  
  !<output>
  ! area of the domain
  real(DP), intent(OUT) :: Dmeasure
  !</output>

!</subroutine>

  ! local variables
  integer(I32) :: i

    Dmeasure=0D0
    do i=1,Nel
      Dmeasure=Dmeasure+Done(i)
    end do
      
  end subroutine GetMeasure


!************************************************************************

!<subroutine>

  subroutine GetAvgSizeElems(rtriangulation,Kcount,Dsize,Dsizem,Dratio)

!<description>
  ! Calculates the average size of elements surrounding a node and 
  ! normalize it such that 1d0 is the max. over the whole vector dsize
!</description>

 include 'cout.inc'

!<input>
    type (t_triangulation), intent(IN) :: rtriangulation
!</input>

!<output>
  ! area of the domain
  real(DP), intent(OUT), dimension(:) :: Dratio,Dsizem,Dsize
  integer(I32), intent(OUT), dimension(:):: Kcount
!</output>

!</subroutine>

  ! local variables
  real(DP), dimension(:,:), pointer :: p_Dcorvg
  integer(I32), dimension(:,:), pointer :: p_Kvert
  integer(I32), dimension(:,:), pointer :: p_Kadj
  integer(I32) :: iel, ivt, ive, jel
  real(DP) :: dax, day, dnx, dny, dnrmn, dpx, dpy, dist, dist2
  real(DP) :: dmaxs, dmins, dquot
  
    call storage_getbase_int2D(rtriangulation%h_IverticesAtElement,p_Kvert)
    call storage_getbase_double2D(rtriangulation%h_DcornerCoordinates,p_Dcorvg)
    call storage_getbase_int2D(rtriangulation%h_IneighboursAtElement,p_Kadj)

    ! Initialise the arrays with 0
    call lalg_vectorClearDble (Dsize)
    call lalg_vectorClearInt (Kcount)

    do iel=1,rtriangulation%NEL

      ! Get first point in quad-element
      dax=p_Dcorvg(1,p_Kvert(1,iel))
      day=p_Dcorvg(2,p_Kvert(1,iel))

      ! Get 3rd point, build normal vector (DNX,DNY) to the diagonal of the
      ! element, normalise it
      dnx= p_Dcorvg(2,p_Kvert(3,iel))-day
      dny=-p_Dcorvg(1,p_Kvert(3,iel))+dax
      dnrmn=DSQRT(dnx*dnx+dny*dny)
      dnx=dnx/dnrmn
      dny=dny/dnrmn
        
      ! triangle 1,2,3
      dpx=p_dcorvg(1,p_kvert(2,iel))
      dpy=p_dcorvg(2,p_kvert(2,iel))
      dist=abs(dnx*dpx+dny*dpy-dnx*dax-dny*day)

      ! triangle 1,3,4
      dpx=p_dcorvg(1,p_kvert(4,iel))
      dpy=p_dcorvg(2,p_kvert(4,iel))
      dist2=abs(dnx*dpx+dny*dpy-dnx*dax-dny*day)
        
      ! Calculate size of the element
      Dsizem(iel)=0.5d0*dnrmn*(dist+dist2)
        
      ! Add the area information to each vertex, building an array with
      ! the sum of areas meeting in a vertex

      Dsize(p_Kvert(1,iel))=Dsize(p_Kvert(1,iel))+0.5d0*dnrmn*(dist+dist2)
      Dsize(p_Kvert(2,iel))=Dsize(p_Kvert(2,iel))+0.5d0*dnrmn*(dist+dist2)
      Dsize(p_Kvert(3,iel))=Dsize(p_Kvert(3,iel))+0.5d0*dnrmn*(dist+dist2)
      Dsize(p_Kvert(4,iel))=Dsize(p_Kvert(4,iel))+0.5d0*dnrmn*(dist+dist2)
      Kcount(p_Kvert(1,iel))=Kcount(p_Kvert(1,iel))+1
      Kcount(p_Kvert(2,iel))=Kcount(p_Kvert(2,iel))+1
      Kcount(p_Kvert(3,iel))=Kcount(p_Kvert(3,iel))+1
      Kcount(p_Kvert(4,iel))=Kcount(p_Kvert(4,iel))+1

    end do      

    ! Build the average size of the elements around each vertex:

    ! Divide the summed area by the number of elements meeting in a vertex
    ! and normalise it with the maximum area.

    dmaxs=0d0
    dmins=1d0
    do ivt=1,rtriangulation%NVT
      Dsize(ivt)=Dsize(ivt)/dble(Kcount(ivt))
      dmaxs=max(dmaxs,Dsize(ivt))
      dmins=min(dmins,Dsize(ivt))
    end do

    do ivt=1,rtriangulation%NVT
      Dsize(ivt)=Dsize(ivt)/dmaxs
    end do

    ! Build the vector with the ratio between the area of each element
    ! and its surrounding elements

    do iel=1,rtriangulation%NEL
      Dratio(iel)=0d0
      do ive=1,4
        jel=p_Kadj(ive,iel)
        if (jel.ne.0) then
          ! ratio between area of current and its neighbour element
          dquot=Dsizem(iel)/Dsizem(jel)
          if (dquot.lt.1d0) dquot=1d0/dquot
          ! find the maximum ratio
          if (dquot.gt.Dratio(iel)) Dratio(iel)=dquot
        end if
      end do
    end do

  end subroutine GetAvgSizeElems

end module 
    

