MODULE domeasures

  USE fsystem
  USE triangulation

  IMPLICIT NONE

CONTAINS

!************************************************************************

!<subroutine>

  SUBROUTINE GetMeasure(Nel,Done,Dmeasure)

  !<description>
  ! This subroutine is used in the grid adaption for measuring the  
  ! domain area
  !</description>

  !<input>
  ! Number of elements
  INTEGER, INTENT(IN) :: Nel
  ! Array of element areas
  REAL(DP), DIMENSION(:), INTENT (IN) :: Done
  !</input>
  
  !<output>
  ! area of the domain
  REAL(DP), INTENT(OUT) :: Dmeasure
  !</output>

!</subroutine>

  ! local variables
  INTEGER(I32) :: i

    Dmeasure=0D0
    DO i=1,Nel
      Dmeasure=Dmeasure+Done(i)
    END DO
      
  END SUBROUTINE GetMeasure


!************************************************************************

!<subroutine>

  SUBROUTINE GetAvgSizeElems(rtriangulation,Kcount,Dsize,Dsizem,Dratio)

!<description>
  ! Calculates the average size of elements surrounding a node and 
  ! normalize it such that 1d0 is the max. over the whole vector dsize
!</description>

 INCLUDE 'cout.inc'

!<input>
    TYPE (t_triangulation), INTENT(IN) :: rtriangulation
!</input>

!<output>
  ! area of the domain
  REAL(DP), INTENT(OUT), DIMENSION(:) :: Dratio,Dsizem,Dsize
  INTEGER(I32), INTENT(OUT), DIMENSION(:):: Kcount
!</output>

!</subroutine>

  ! local variables
  REAL(DP), DIMENSION(:,:), POINTER :: p_Dcorvg
  INTEGER(I32), DIMENSION(:,:), POINTER :: p_Kvert
  INTEGER(I32), DIMENSION(:,:), POINTER :: p_Kadj
  INTEGER(I32) :: iel, ivt, ive, jel
  REAL(DP) :: dax, day, dnx, dny, dnrmn, dpx, dpy, dist, dist2
  REAL(DP) :: dmaxs, dmins, dquot
  
    CALL storage_getbase_int2D(rtriangulation%h_IverticesAtElement,p_Kvert)
    CALL storage_getbase_double2D(rtriangulation%h_DcornerCoordinates,p_Dcorvg)
    CALL storage_getbase_int2D(rtriangulation%h_IneighboursAtElement,p_Kadj)

    ! Initialise the arrays with 0
    CALL lalg_vectorClearDble (Dsize)
    CALL lalg_vectorClearInt (Kcount)

    DO iel=1,rtriangulation%NEL

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
      dist=ABS(dnx*dpx+dny*dpy-dnx*dax-dny*day)

      ! triangle 1,3,4
      dpx=p_dcorvg(1,p_kvert(4,iel))
      dpy=p_dcorvg(2,p_kvert(4,iel))
      dist2=ABS(dnx*dpx+dny*dpy-dnx*dax-dny*day)
        
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

    END DO      

    ! Build the average size of the elements around each vertex:

    ! Divide the summed area by the number of elements meeting in a vertex
    ! and normalise it with the maximum area.

    dmaxs=0d0
    dmins=1d0
    DO ivt=1,rtriangulation%NVT
      Dsize(ivt)=Dsize(ivt)/DBLE(Kcount(ivt))
      dmaxs=max(dmaxs,Dsize(ivt))
      dmins=min(dmins,Dsize(ivt))
    END DO

    DO ivt=1,rtriangulation%NVT
      Dsize(ivt)=Dsize(ivt)/dmaxs
    END DO

    ! Build the vector with the ratio between the area of each element
    ! and its surrounding elements

    DO iel=1,rtriangulation%NEL
      Dratio(iel)=0d0
      DO ive=1,4
        jel=p_Kadj(ive,iel)
        IF (jel.ne.0) THEN
          ! ratio between area of current and its neighbour element
          dquot=Dsizem(iel)/Dsizem(jel)
          IF (dquot.lt.1d0) dquot=1d0/dquot
          ! find the maximum ratio
          IF (dquot.gt.Dratio(iel)) Dratio(iel)=dquot
        END IF
      END DO
    END DO

  END SUBROUTINE GetAvgSizeElems

END MODULE 
    

