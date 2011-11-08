!##############################################################################
!# ****************************************************************************
!# <name> disto3d_aux </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains a hand full of auxiliary routines for the
!# elemdbg3d_testX debugger modules.
!#
!# The following routines can be found in this module:
!#
!# .) disto3d_calcDomainVolume
!#    -> Calculates the total volume of the domain represented by a
!#       triangulation by integration of the 1-function.
!#
!# .) disto3d_calcFaceNormalDeviation
!#    -> Calculates the maximal and mean face normal deviation (in respect to
!#       the face's midpoint normal) over all faces of a given triangulation.
!#
!# .) disto3d_calcFaceNLDeviation
!#    -> Calculates the maximal and mean non-linear deviation of all faces
!#       of a given triangulation.
!#
!# .) disto3d_distortCubeCoords
!#    -> Distorts the vertices of a CUBE mesh by coordinate-based stochastical
!#       distortion.
!#
!# .) disto3d_distortCubePlaneXY
!#    -> Distorts the vertices of a CUBE mesh by XY-plane-wise index-based
!#       stochastical distortion.
!#
!# .) disto3d_distortCubePlaneXZ
!#    -> Distorts the vertices of a CUBE mesh by XZ-plane-wise index-based
!#       stochastical distortion.
!#
!# .) disto3d_distortCubePlaneYZ
!#    -> Distorts the vertices of a CUBE mesh by YZ-plane-wise index-based
!#       stochastical distortion.
!#
!# </purpose>
!##############################################################################

module disto3d_aux

use fsystem
use genoutput
use storage
use cubature
use basicgeometry
use transformation
use triangulation

implicit none

contains

  ! ***************************************************************************

!<subroutine>

  pure subroutine disto3d_calcFaceTrafo(Dv, Dtrafo)

!<description>
  ! AUXILIARY ROUTINE:
  ! Calculates the transformation coefficients for a bilinear face trafo from
  ! the reference quadrilateral onto a real face.
!</description>

!<input>
  ! The coordinates of the four corner vertices of the face.
  real(DP), dimension(3,4), intent(IN) :: Dv
!</input>

!<output>
  ! The transformation coefficients for the bilinear trafo.
  real(DP), dimension(4,3), intent(OUT) :: Dtrafo
!</output>

!</subroutine>
  
  integer :: i

    do i = 1, 3
      Dtrafo(1,i) = 0.25_DP*( Dv(i,1)+Dv(i,2)+Dv(i,3)+Dv(i,4))
      Dtrafo(2,i) = 0.25_DP*(-Dv(i,1)+Dv(i,2)+Dv(i,3)-Dv(i,4))
      Dtrafo(3,i) = 0.25_DP*(-Dv(i,1)-Dv(i,2)+Dv(i,3)+Dv(i,4))
      Dtrafo(4,i) = 0.25_DP*( Dv(i,1)-Dv(i,2)+Dv(i,3)-Dv(i,4))
    end do
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine disto3d_calcFaceNormal(Dtrafo, dx, dy, Dnormal)

!<description>
  ! AUXILIARY ROUTINE:
  ! Calculates the normal vector on a face.
!</description>

!<input>
  ! The transformation coefficients of the face, as returned by the
  ! disto3d_calcFaceTrafo routine.
  real(DP), dimension(4,3), intent(IN) :: Dtrafo
  
  ! A point on the reference quadrilateral [-1,1]x[-1,1] for which the normal
  ! vector is to be calculated.
  real(DP), intent(IN) :: dx, dy
!</input>

!<output>
  ! The normal vector in the point Dx.
  real(DP), dimension(3), intent(OUT) :: Dnormal
!</output>

!</subroutine>

  ! jacobian matrix
  real(DP), dimension(3,2) :: Djac
  real(DP) :: daux
  integer :: i
  
    ! Calculate the jacobian matrix in the point Dx.
    do i = 1, 3
      Djac(i,1) = Dtrafo(2,i) + Dtrafo(4,i)*dy
      Djac(i,2) = Dtrafo(3,i) + Dtrafo(4,i)*dx
    end do
    
    ! Calculate the normal vector as the cross-product of the jacobian matrix's
    ! columns.
    Dnormal(1) = Djac(2,1)*Djac(3,2) - Djac(3,1)*Djac(2,2)
    Dnormal(2) = Djac(3,1)*Djac(1,2) - Djac(1,1)*Djac(3,2)
    Dnormal(3) = Djac(1,1)*Djac(2,2) - Djac(2,1)*Djac(1,2)
    
    ! Calculate the length of the normal vector
    daux = sqrt(Dnormal(1)**2 + Dnormal(2)**2 + Dnormal(3)**2)
    if(daux .gt. 0.0_DP) then
      ! Normalise the vector
      Dnormal = Dnormal * (1.0_DP / daux)
    end if

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine disto3d_calcDomainVolume(rtria, dvolume)

!<description>
  ! Calculates the volume of the domain represented by the triangulation.
!</description>

!<input>
  ! The triangulation whose total volume is to be determined.
  type(t_triangulation), intent(IN) :: rtria
!</input>

!<output>
  ! The total volume of the domain.
  real(DP), intent(OUT) :: dvolume
!</output>

!</subroutine>

  ! Parameter: Cubature formula
  integer, parameter :: ccubature = CUB_G2_3D

  ! Arrays for the cubature formula
  integer :: ncubp
  real(DP), dimension(:,:), allocatable :: DcubPts
  real(DP), dimension(:), allocatable :: Domega
  
  ! Arrays for the triangulation
  real(DP), dimension(:,:), pointer :: p_Dvtx
  integer, dimension(:,:), pointer :: p_Ielem
  
  ! Array for the corner vertice coordinates
  real(DP), dimension(3,8) :: Dcorners
  
  ! Array for the transformation
  real(DP), dimension(TRAFO_NAUXJAC3D) :: DjacPrep
  
  ! Dummy array for jacobian matrix
  real(DP), dimension(9) :: Djac
  
  ! Other local variables
  integer :: iel, i
  real(DP) :: ddet

    ! Make sure the triangulation is 3D.
    if(rtria%ndim .ne. NDIM3D) then
      call output_line ('Mesh must be a 3D mesh!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'disto3d_calcDomainVolume')
      call sys_halt()
    end if

    ! Make sure the mesh consists only of hexahedra.
    if(rtria%NEL .ne. rtria%InelOfType(TRIA_NVEHEXA3D)) then
      call output_line ('Mesh must contain only hexahedra!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'disto3d_calcDomainVolume')
      call sys_halt()
    end if
    
    ! Get the arrays from the triangulation
    call storage_getbase_double2D(rtria%h_DvertexCoords, p_Dvtx)
    call storage_getbase_int2D(rtria%h_IverticesAtElement, p_Ielem)
    
    ! Init the cubature formula
    ncubp = cub_igetNumPts(ccubature)
    allocate(DcubPts(3,ncubp))
    allocate(Domega(ncubp))
    call cub_getCubature(ccubature, DcubPts, Domega)
    
    ! Loop over all elements
    dvolume = 0.0_DP
    !$omp parallel do private(i,Dcorners,DjacPrep,Djac,ddet) reduction(+:dvolume)
    do iel = 1, rtria%NEL
    
      ! Get the element's corner vertice coordinates
      do i = 1, 8
        Dcorners(:,i) = p_Dvtx(:,p_Ielem(i,iel))
      end do
      
      ! Calculate trafo
      call trafo_prepJac_hexa3D(Dcorners, DjacPrep)
      
      ! Calculate determinants for all points and perform the integration
      do i = 1, ncubp
        call trafo_calcJac_hexa3D (DjacPrep, Djac, ddet, &
                     DcubPts(1,i), DcubPts(2,i), DcubPts(3,i))
                     
        dvolume = dvolume + abs(ddet)*Domega(i)
      end do
    
    end do ! iel
    !$omp end parallel do
    
    deallocate(Domega)
    deallocate(DcubPts)
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine disto3d_calcFaceNormalDeviation(rtria, dmaxDev, dmeanDev, ccubature)

!<description>
  ! Calculates the mean and maximum normal vector deviation for the faces of
  ! a given triangulation.
  !
  ! Let T:[-1,1]^2 -> R^3 denote the bilinear transformation from the
  ! reference quadrilateral onto a real face and let N:[-1,1]^2 -> R^3 denote
  ! the normal vector field of T.
  !
  ! Then the normal deviation of the transformation is defined as:
  !
  !                          k
  !        Normal_Dev(T) := max || N(x_i,y_i) - N(0,0) ||_2
  !                         i=1
  !
  ! for some set { (x_1,y_1),...,(x_k,y_k) } of points on the reference quad.
  !
  ! Remark:
  ! If the face is planar, then the normal deviation of the face is always 0.
!</description>

!<input>
  ! The triangulation whose face normal deviation is to be determined.
  type(t_triangulation), intent(IN) :: rtria
  
  ! OPTIONAL: A cubature formula identifier which specifies in which points
  ! the normals are to be calculated. If given, must be a valid quadrilateral
  ! cubatuture rule. If not specified, the trapezoidal rule is used.
  integer, optional, intent(IN) :: ccubature
!</input>

!<output>
  ! The maximum of all face normal deviations.
  real(DP), intent(OUT) :: dmaxDev
  
  ! The mean of all face normal deviations.
  real(DP), intent(OUT) :: dmeanDev
!</output>

!</subroutine>

  ! Arrays for the cubature formula
  integer :: ccub,ncubp
  real(DP), dimension(:,:), allocatable :: DcubPts
  real(DP), dimension(:), allocatable :: Domega
  
  ! Arrays for the triangulation
  real(DP), dimension(:,:), pointer :: p_Dvtx
  integer, dimension(:,:), pointer :: p_Iface
  real(DP), dimension(3,4) :: Dcorners
  
  ! Array for the bilinear trafo
  real(DP), dimension(4,3) :: Dtrafo
  
  ! Arrays for the normal vectors
  real(DP), dimension(3) :: Dnormal1, Dnormal2
  
  ! Other local variables
  integer :: iat, i
  real(DP) :: daux, daux2
  
    ! Make sure the triangulation is 3D.
    if(rtria%ndim .ne. NDIM3D) then
      call output_line ('Mesh must be a 3D mesh!', &
          OU_CLASS_ERROR,OU_MODE_STD,'disto3d_calcFaceNormalDeviation')
      call sys_halt()
    end if

    ! Make sure the mesh consists only of hexahedra.
    if(rtria%NEL .ne. rtria%InelOfType(TRIA_NVEHEXA3D)) then
      call output_line ('Mesh must contain only hexahedra!', &
          OU_CLASS_ERROR,OU_MODE_STD,'disto3d_calcFaceNormalDeviation')
      call sys_halt()
    end if
    
    ! Get the arrays from the triangulation
    call storage_getbase_double2D(rtria%h_DvertexCoords, p_Dvtx)
    call storage_getbase_int2D(rtria%h_IverticesAtFace, p_Iface)
    
    ! Init the cubature formula
    ccub = CUB_TRZ
    if(present(ccubature)) then
      
      ! Check the shape
      if(cub_igetShape(ccubature) .ne. BGEOM_SHAPE_QUAD) then
        call output_line ('Cubature rule is not a valid quadrilateral rule!', &
            OU_CLASS_ERROR,OU_MODE_STD,'disto3d_calcFaceNormalDeviation')
        call sys_halt()
      else
        ccub = ccubature
      end if
    end if

    ncubp = cub_igetNumPts(ccub)
    allocate(DcubPts(2,ncubp))
    allocate(Domega(ncubp))
    call cub_getCubature(ccub, DcubPts, Domega)
    
    ! Deallocate Domega, as we don't need it
    deallocate(Domega)
    
    ! Reset the deviations
    dmaxDev = 0.0_DP
    dmeanDev = 0.0_DP

    ! Loop over all faces
    do iat = 1, rtria%NAT
    
      ! Get the corner vertice coordinates
      do i = 1, 4
        Dcorners(:,i) = p_Dvtx(:,p_Iface(i,iat))
      end do
    
      ! Calculate the trafo coefficients
      call disto3d_calcFaceTrafo(Dcorners, Dtrafo)
      
      ! Calculate the normal vector in the face midpoint
      call disto3d_calcFaceNormal(Dtrafo, 0.0_DP, 0.0_DP, Dnormal1)
      
      ! Loop over all other normal vectors
      daux = 0.0_DP
      do i = 1, ncubp
        
        ! Calculate normal vector in cubature point
        call disto3d_calcFaceNormal(Dtrafo, DcubPts(1,i), DcubPts(2,i), &
                                    Dnormal2)
        
        ! Subtract mid-point normal vector
        Dnormal2 = Dnormal2 - Dnormal1
        
        ! Calculate euclid norm of deviation vector
        daux2 = sqrt(Dnormal2(1)**2 + Dnormal2(2)**2 + Dnormal2(3)**2)
        
        ! Calculate maximum for the current face
        daux = max(daux, daux2)
        
      end do
      
      ! Incorporate local deviation to global ones
      dmaxDev = max(dmaxDev, daux)
      dmeanDev = dmeanDev + daux
    
    end do ! iat
    
    ! Deallocate cubature points
    deallocate(DcubPts)
    
    ! Calculate mean deviation
    if(rtria%NAT .gt. 0) dmeanDev = dmeanDev / real(rtria%NAT,DP)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine disto3d_calcFaceNLDeviation(rtria, dmaxDev, dmeanDev)

!<description>
  ! Calculates the mean and maximum non-linear deviation for the faces of
  ! a given triangulation.
  !
  ! Let T denote the bilinear transformation from the reference quadrilateral
  ! onto a real face in the triangulation of the form:
  !
  !        T: [-1,1]^2 -> R^3
  !
  !        T(x,y) = a1 + a2*x + a3*y + a4*x*y
  !
  ! with a1,...,a4 in R^3 being the corresponding transformation vectors.
  !
  ! Then the non-linear deviation of the transformation is defined as:
  !
  !        NL_Dev(T) := || a4 ||_2 / || a2 X a3 ||_2
  !
  ! where 'X' denotes the 3D cross-product (outer product).
  !
  ! Remark:
  ! If a face is a parallelogram, then its non-linear deviation is 0.
!</description>

!<input>
  ! The triangulation whose face non-linear deviation is to be determined.
  type(t_triangulation), intent(IN) :: rtria
!</input>

!<output>
  ! The maximum of all face non-linear deviations.
  real(DP), intent(OUT) :: dmaxDev
  
  ! The mean of all face non-linear deviations.
  real(DP), intent(OUT) :: dmeanDev
!</output>

!</subroutine>
  
  ! Arrays for the triangulation
  real(DP), dimension(:,:), pointer :: p_Dvtx
  integer, dimension(:,:), pointer :: p_Iface
  real(DP), dimension(3,4) :: Dcorners
  
  ! Arrays for the bilinear trafo
  real(DP), dimension(4,3) :: Dtrafo
  
  ! Other local variables
  integer :: iat, i
  real(DP) :: daux, daux2
  
    ! Make sure the triangulation is 3D.
    if(rtria%ndim .ne. NDIM3D) then
      call output_line ('Mesh must be a 3D mesh!', &
          OU_CLASS_ERROR,OU_MODE_STD,'disto3d_calcFaceNLDeviation')
      call sys_halt()
    end if

    ! Make sure the mesh consists only of hexahedra.
    if(rtria%NEL .ne. rtria%InelOfType(TRIA_NVEHEXA3D)) then
      call output_line ('Mesh must contain only hexahedra!', &
          OU_CLASS_ERROR,OU_MODE_STD,'disto3d_calcFaceNLDeviation')
      call sys_halt()
    end if
    
    ! Get the arrays from the triangulation
    call storage_getbase_double2D(rtria%h_DvertexCoords, p_Dvtx)
    call storage_getbase_int2D(rtria%h_IverticesAtFace, p_Iface)
    
    ! Reset the deviations
    dmaxDev = 0.0_DP
    dmeanDev = 0.0_DP

    ! Loop over all faces
    do iat = 1, rtria%NAT
    
      ! Get the corner vertice coordinates
      do i = 1, 4
        Dcorners(:,i) = p_Dvtx(:,p_Iface(i,iat))
      end do
    
      ! Calculate the trafo coefficients
      call disto3d_calcFaceTrafo(Dcorners, Dtrafo)
      
      ! Calculate norm of linear trafo factors
      daux2 = (Dtrafo(2,2)*Dtrafo(3,3) - Dtrafo(2,3)*Dtrafo(3,2))**2 &
            + (Dtrafo(2,3)*Dtrafo(3,1) - Dtrafo(2,1)*Dtrafo(3,3))**2 &
            + (Dtrafo(2,1)*Dtrafo(3,2) - Dtrafo(2,2)*Dtrafo(3,1))**2
            
      ! Jump out if the linear factors are zero
      if(daux2 .eq. 0.0_DP) then
        dmaxDev = SYS_INFINITY_DP
        dmeanDev = SYS_INFINITY_DP
        return
      end if
      
      ! Calculate norm of non-linear factors
      daux = Dtrafo(4,1)**2 + Dtrafo(4,2)**2 + Dtrafo(4,3)**2
      
      ! Calculate non-linear deviation
      daux = sqrt(daux / daux2)
    
      ! Incorporate local deviation to global ones
      dmaxDev = max(dmaxDev, daux)
      dmeanDev = dmeanDev + daux
        
    end do ! iat
    
    ! Calculate mean deviation
    if(rtria%NAT .gt. 0) dmeanDev = dmeanDev / real(rtria%NAT,DP)
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine disto3d_distortCubeCoords(rtria, ddist)

!<description>
  ! Distorts the vertices of a CUBE mesh by coordinate-based stochastical
  ! distortion.
!</description>

!<input>
  ! The distortion parameter
  real(DP), intent(IN) :: ddist
!</input>

!<inputoutput>
  ! The CUBE mesh that is to be distorted.
  type(t_triangulation), intent(INOUT) :: rtria
!</inputoutput>

!</subroutine>

  real(DP), parameter :: dtol = 1e-8_DP

  ! local variables
  real(DP), dimension(:,:), pointer :: p_Dvtx
  real(DP) :: dhdist,dh2
  integer :: ivt,iX,iY,iZ

    ! Calculate number of vertices in each dimension
    dh2 = real(rtria%NVT,DP)**(1.0_DP/3.0_DP)
    
    ! Calculate distortion parameter
    dhdist = ddist / (dh2 + 1.0_DP)
    if(dhdist .eq. 0.0_DP) return
    
    ! Get arrays from the triangulation
    call storage_getbase_double2d (rtria%h_DvertexCoords, p_Dvtx)
    
    ! Loop over the vertices
    do ivt = 1, rtria%NVT
    
      ! Calculate the X-, Y- and Z-position of this vertice
      iX = int(p_Dvtx(1,ivt) * dh2) + 1
      iY = int(p_Dvtx(2,ivt) * dh2) + 1
      iZ = int(p_Dvtx(3,ivt) * dh2) + 1
      
      ! Distort the vertice's coordiantes.
      ! We also need to distort the
      if((p_Dvtx(1,ivt) .gt. dtol) .and. (p_Dvtx(1,ivt)+dtol .lt. 1.0_DP)) &
        p_Dvtx(1,ivt) = p_Dvtx(1,ivt) + real((-1)**mod(3*iX,17),DP)*dhdist
      if((p_Dvtx(2,ivt) .gt. dtol) .and. (p_Dvtx(2,ivt)+dtol .lt. 1.0_DP)) &
        p_Dvtx(2,ivt) = p_Dvtx(2,ivt) + real((-1)**mod(5*iY,17),DP)*dhdist
      if((p_Dvtx(3,ivt) .gt. dtol) .and. (p_Dvtx(3,ivt)+dtol .lt. 1.0_DP)) &
        p_Dvtx(3,ivt) = p_Dvtx(3,ivt) + real((-1)**mod(7*iZ,17),DP)*dhdist
        
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine disto3d_distortCubePlaneXY(rtria, ddist, dsecDist)

!<description>
  ! Distorts the vertices of a CUBE mesh by XY-plane-wise index-based
  ! stochastical distortion.
!</description>

!<input>
  ! The distortion parameters
  real(DP), intent(IN) :: ddist, dsecDist
!</input>

!<inputoutput>
  ! The CUBE mesh that is to be distorted.
  type(t_triangulation), intent(INOUT) :: rtria
!</inputoutput>

!</subroutine>

  real(DP), parameter :: dtol = 1e-8_DP

  ! local variables
  real(DP), dimension(:,:), pointer :: p_Dvtx
  real(DP) :: dhdist,dhdist2,dh2
  integer :: ivt,iX,iY,iZ,idx,invt

    ! Calculate number of vertices in each dimension
    dh2 = real(rtria%NVT,DP)**(1.0_DP/3.0_DP)
    
    ! Get number of vertices
    invt = int(dh2) + 10
    
    ! Calculate distortion parameters
    dhdist = ddist / (dh2 + 1.0_DP)
    dhdist2 = dsecDist / (dh2 + 1.0_DP)
    if((dhdist .eq. 0.0_DP) .and. (dhdist2 .eq. 0.0_DP)) return
    
    ! Get arrays from the triangulation
    call storage_getbase_double2d (rtria%h_DvertexCoords, p_Dvtx)
    
    ! Loop over the vertices
    do ivt = 1, rtria%NVT
    
      ! Calculate the X-, Y- and Z-position of this vertice
      iX = int(p_Dvtx(1,ivt) * dh2) + 1
      iY = int(p_Dvtx(2,ivt) * dh2) + 1
      iZ = int(p_Dvtx(3,ivt) * dh2) + 1
      
      ! Calculate plane index
      idx = iX*invt + iY
      
      ! Distort the vertice's coordiantes.
      ! We also need to distort the
      if((p_Dvtx(1,ivt) .gt. dtol) .and. (p_Dvtx(1,ivt)+dtol .lt. 1.0_DP)) &
        p_Dvtx(1,ivt) = p_Dvtx(1,ivt) + real((-1)**mod(idx,17),DP)*dhdist
      if((p_Dvtx(2,ivt) .gt. dtol) .and. (p_Dvtx(2,ivt)+dtol .lt. 1.0_DP)) &
        p_Dvtx(2,ivt) = p_Dvtx(2,ivt) + real((-1)**mod(idx+7,17),DP)*dhdist
      if((p_Dvtx(3,ivt) .gt. dtol) .and. (p_Dvtx(3,ivt)+dtol .lt. 1.0_DP)) &
        p_Dvtx(3,ivt) = p_Dvtx(3,ivt) + real((-1)**mod(ivt,17),DP)*dhdist2
        
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine disto3d_distortCubePlaneXZ(rtria, ddist, dsecDist)

!<description>
  ! Distorts the vertices of a CUBE mesh by XZ-plane-wise index-based
  ! stochastical distortion.
!</description>

!<input>
  ! The distortion parameters
  real(DP), intent(IN) :: ddist, dsecDist
!</input>

!<inputoutput>
  ! The CUBE mesh that is to be distorted.
  type(t_triangulation), intent(INOUT) :: rtria
!</inputoutput>

!</subroutine>

  real(DP), parameter :: dtol = 1e-8_DP

  ! local variables
  real(DP), dimension(:,:), pointer :: p_Dvtx
  real(DP) :: dhdist,dhdist2,dh2
  integer :: ivt,iX,iY,iZ,idx,invt

    ! Calculate number of vertices in each dimension
    dh2 = real(rtria%NVT,DP)**(1.0_DP/3.0_DP)
    
    ! Get number of vertices
    invt = int(dh2) + 10
    
    ! Calculate distortion parameters
    dhdist = ddist / (dh2 + 1.0_DP)
    dhdist2 = dsecDist / (dh2 + 1.0_DP)
    if((dhdist .eq. 0.0_DP) .and. (dhdist2 .eq. 0.0_DP)) return
    
    ! Get arrays from the triangulation
    call storage_getbase_double2d (rtria%h_DvertexCoords, p_Dvtx)
    
    ! Loop over the vertices
    do ivt = 1, rtria%NVT
    
      ! Calculate the X-, Y- and Z-position of this vertice
      iX = int(p_Dvtx(1,ivt) * dh2) + 1
      iY = int(p_Dvtx(2,ivt) * dh2) + 1
      iZ = int(p_Dvtx(3,ivt) * dh2) + 1
      
      ! Calculate plane index
      idx = iX*invt + iZ
      
      ! Distort the vertice's coordiantes.
      ! We also need to distort the
      if((p_Dvtx(1,ivt) .gt. dtol) .and. (p_Dvtx(1,ivt)+dtol .lt. 1.0_DP)) &
        p_Dvtx(1,ivt) = p_Dvtx(1,ivt) + real((-1)**mod(idx,17),DP)*dhdist
      if((p_Dvtx(2,ivt) .gt. dtol) .and. (p_Dvtx(2,ivt)+dtol .lt. 1.0_DP)) &
        p_Dvtx(2,ivt) = p_Dvtx(2,ivt) + real((-1)**mod(ivt,17),DP)*dhdist2
      if((p_Dvtx(3,ivt) .gt. dtol) .and. (p_Dvtx(3,ivt)+dtol .lt. 1.0_DP)) &
        p_Dvtx(3,ivt) = p_Dvtx(3,ivt) + real((-1)**mod(idx+7,17),DP)*dhdist
        
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine disto3d_distortCubePlaneYZ(rtria, ddist, dsecDist)

!<description>
  ! Distorts the vertices of a CUBE mesh by YZ-plane-wise index-based
  ! stochastical distortion.
!</description>

!<input>
  ! The distortion parameters
  real(DP), intent(IN) :: ddist, dsecDist
!</input>

!<inputoutput>
  ! The CUBE mesh that is to be distorted.
  type(t_triangulation), intent(INOUT) :: rtria
!</inputoutput>

!</subroutine>

  real(DP), parameter :: dtol = 1e-8_DP

  ! local variables
  real(DP), dimension(:,:), pointer :: p_Dvtx
  real(DP) :: dhdist,dhdist2,dh2
  integer :: ivt,iX,iY,iZ,idx,invt

    ! Calculate number of vertices in each dimension
    dh2 = real(rtria%NVT,DP)**(1.0_DP/3.0_DP)
    
    ! Get number of vertices
    invt = int(dh2) + 10
    
    ! Calculate distortion parameters
    dhdist = ddist / (dh2 + 1.0_DP)
    dhdist2 = dsecDist / (dh2 + 1.0_DP)
    if((dhdist .eq. 0.0_DP) .and. (dhdist2 .eq. 0.0_DP)) return
    
    ! Get arrays from the triangulation
    call storage_getbase_double2d (rtria%h_DvertexCoords, p_Dvtx)
    
    ! Loop over the vertices
    do ivt = 1, rtria%NVT
    
      ! Calculate the X-, Y- and Z-position of this vertice
      iX = int(p_Dvtx(1,ivt) * dh2) + 1
      iY = int(p_Dvtx(2,ivt) * dh2) + 1
      iZ = int(p_Dvtx(3,ivt) * dh2) + 1
      
      ! Calculate plane index
      idx = iY*invt + iZ
      
      ! Distort the vertice's coordiantes.
      ! We also need to distort the
      if((p_Dvtx(1,ivt) .gt. dtol) .and. (p_Dvtx(1,ivt)+dtol .lt. 1.0_DP)) &
        p_Dvtx(1,ivt) = p_Dvtx(1,ivt) + real((-1)**mod(ivt,17),DP)*dhdist2
      if((p_Dvtx(2,ivt) .gt. dtol) .and. (p_Dvtx(2,ivt)+dtol .lt. 1.0_DP)) &
        p_Dvtx(2,ivt) = p_Dvtx(2,ivt) + real((-1)**mod(idx,17),DP)*dhdist
      if((p_Dvtx(3,ivt) .gt. dtol) .and. (p_Dvtx(3,ivt)+dtol .lt. 1.0_DP)) &
        p_Dvtx(3,ivt) = p_Dvtx(3,ivt) + real((-1)**mod(idx+7,17),DP)*dhdist
        
    end do

  end subroutine

end module