!#########################################################################
!# ***********************************************************************
!# <name> pprocindicator </name>
!# ***********************************************************************
!#
!# <purpose>
!# This module contains various routines for calculating error indicators.
!#
!# The following routines can be found in this module:
!#
!# 1.) ppind_SecondDifference
!#     -> Calculates the second-difference indicator by R. Loehner
!# </purpose>
!#########################################################################

module pprocindicator

!$use omp_lib
  use basicgeometry
  use element
  use fsystem
  use genoutput
  use linearalgebra
  use linearsystemscalar
  use spatialdiscretisation
  use storage
  use triangulation

  implicit none

  private

  public :: ppind_secondDifference

contains

  !*****************************************************************************

!<subroutine>

  subroutine ppind_secondDifference(rvectorScalar, dnoiseFilter,&
                                    dabsFilter, rindicator)

!<description>
    ! This subroutine computes the second-difference indicator proposed
    ! by R. Loehner. It is restricted to linear and multilinear finite
    ! elements in arbitrary space dimensions.
!</description>

!<input>
    ! FE solution vector
    type(t_vectorScalar), intent(in) :: rvectorScalar

    ! parameter for noise filtering
    real(DP), intent(in) :: dnoiseFilter

    ! parameter for absolute filtering
    real(DP), intent(in) :: dabsFilter
!</input>

!<output>
    ! indicator vector
    type(t_vectorScalar), intent(out) :: rindicator
!</output>
!</subroutine>

    ! Pointer to the spatial discretization
    type(t_spatialDiscretisation), pointer :: p_rspatialDiscr

    ! Pointer to the triangulation
    type(t_triangulation), pointer :: p_rtriangulation

    ! Pointer to element indicator
    real(DP), dimension(:), pointer :: p_Dindicator

    ! Pointer to scalar error vector
    real(DP), dimension(:), pointer :: p_Ddata

    ! Pointer to vertex coordinates
    real(DP), dimension(:,:), pointer :: p_DvertexCoords

    ! Pointer to vertices at element
    integer, dimension(:,:), pointer :: p_IverticesAtElement

    ! Pointer to neighbours at element
    integer, dimension(:,:), pointer :: p_IneighboursAtElement

    ! Pointer to nodal property
    integer, dimension(:), pointer :: p_InodalProperty

    ! Handle for nodal coefficients
    integer :: h_Dcoefficients

    ! Pointer to p_Dcoefficients
    real(DP), dimension(:,:), pointer :: p_Dcoefficients

    ! Handle for nodal indicator
    integer :: h_DnodalIndicator

    ! Pointer to p_DnodalIndicator
    real(DP), dimension(:), pointer :: p_DnodalIndicator

    ! local variable
    integer, dimension(2) :: Isize
    integer :: i

    ! Set pointer to the underlying spatial discretization
    p_rspatialDiscr => rvectorScalar%p_rspatialDiscr
    if (.not.(associated(p_rspatialDiscr))) then
      call output_line('Vector does not provide discretization structure!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'errest_calcSecondDifferenceIndicator')
      call sys_halt()
    end if

    ! Loop over all FE spaces an check that linear or multi-linear elements are used
    do i = 1, p_rspatialDiscr%inumFESpaces
      select case (p_rspatialDiscr%RelementDistr(i)%celement)

      case (EL_P1_1D,&
            EL_P1_2D, EL_Q1_2D,&
            EL_P1_3D, EL_Q1_3D)
        ! These element types are valid

      case DEFAULT
        call output_line('Unsupported type of finite elements!',&
                         OU_CLASS_ERROR,OU_MODE_STD,'errest_calcSecondDifferenceIndicator')
        call sys_halt()
      end select
    end do

    ! Set pointer to the underlying triangulation
    p_rtriangulation => p_rspatialDiscr%p_rtriangulation
    if (.not.(associated(p_rtriangulation))) then
      call output_line('Discretization structure does not provide triangulation!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'errest_calcSecondDifferenceIndicator')
      call sys_halt()
    end if

    ! Create scalar vector for element indicator
    call lsyssc_createVector(rindicator, p_rtriangulation%NEL, .false.)

    ! Create temporal memory
    h_DnodalIndicator = ST_NOHANDLE
    call storage_new('ppind_secondDifference',' DnodalIndicator',&
                     p_rtriangulation%NVT, ST_DOUBLE, h_DnodalIndicator, ST_NEWBLOCK_NOINIT)
    call storage_getbase_double(h_DnodalIndicator, p_DnodalIndicator)

    ! Set pointers
    call storage_getbase_double2D(p_rtriangulation%h_DvertexCoords,&
                                  p_DvertexCoords)
    call storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement,&
                               p_IverticesAtElement)
    call storage_getbase_int2D(p_rtriangulation%h_IneighboursAtElement,&
                               p_IneighboursAtElement)
    call storage_getbase_int(p_rtriangulation%h_InodalProperty,&
                             p_InodalProperty)
    call lsyssc_getbase_double(rvectorScalar, p_Ddata)
    call lsyssc_getbase_double(rindicator, p_Dindicator)

    ! How many spatial dimensions are we?
    select case (p_rtriangulation%ndim)

    case (NDIM2D)
      ! Create temporal memory
      h_Dcoefficients = ST_NOHANDLE
      Isize = (/12, p_rtriangulation%NVT/)
      call storage_new('ppind_secondDifference',' Dcoefficients',&
                       Isize, ST_DOUBLE, h_Dcoefficients, ST_NEWBLOCK_NOINIT)
      call storage_getbase_double2d(h_Dcoefficients, p_Dcoefficients)

      ! Compute the second-difference indicator in 2D
      call doSecondDiffIndicator2D(p_Ddata, p_DvertexCoords, p_IverticesAtElement,&
                                   p_IneighboursAtElement, p_InodalProperty,&
                                   dnoiseFilter, dabsFilter, p_rtriangulation%NEL,&
                                   p_rtriangulation%NVT, p_Dcoefficients,&
                                   p_DnodalIndicator, p_Dindicator)

      ! Release temporal memory
      call storage_free(h_Dcoefficients)


    case DEFAULT
      call output_line('Unsupported spatial dimension!',&
                        OU_CLASS_ERROR,OU_MODE_STD,'ppind_secondDifference')
      call sys_halt()
    end select

    ! Release temporal memory
    call storage_free(h_DnodalIndicator)

  contains

    ! Here, the real working routines follow.

    !**************************************************************
    ! Second-difference indicator in 2D

    subroutine doSecondDiffIndicator2D(Ddata, DvertexCoords, IverticesAtElement,&
                                       IneighboursAtElement, InodalProperty, &
                                       dweight, dfilter, NEL, NVT, Dcoefficients,&
                                       DnodalIndicator, Dindicator)

      real(DP), dimension(:), intent(in) :: Ddata
      real(DP), dimension(:,:), intent(in) :: DvertexCoords
      integer, dimension(:,:), intent(in) :: IverticesAtElement
      integer, dimension(:,:), intent(in) :: IneighboursAtElement
      integer, dimension(:), intent(in) :: InodalProperty
      real(DP), intent(in) :: dweight,dfilter
      integer, intent(in) :: NEL,NVT

      real(DP), dimension(:,:), intent(out) :: Dcoefficients
      real(DP), dimension(:), intent(out) :: DnodalIndicator
      real(DP), dimension(:), intent(out) :: Dindicator

      ! local variables
      real(DP), dimension(2) :: DabsDeriv,Dderiv,Dfunc
      real(DP), dimension(2,3) :: Dbas,dabsBas
      real(DP) :: darea,ddet,dlength,dbdrInt
      real(DP) :: daux1,daux2,daux3
      real(DP) :: u1,u2,u3,u4,x1,x2,x3,x4,y1,y2,y3,y4
      integer :: iel,ivt,i1,i2,i3,i4,ive,ncontributions



      ! Clear array for coefficients
      call lalg_clearVectorDble2D(Dcoefficients)

      ! Loop over all elements in triangulation
      do iel = 1, NEL

        ! What type of element are we
        select case(tria_getNVE(IverticesAtElement, iel))

        case(TRIA_NVETRI2D)

          ! Determine global degrees of freedom
          i1 = IverticesAtElement(1, iel)
          i2 = IverticesAtElement(2, iel)
          i3 = IverticesAtElement(3, iel)

          ! Determine vertex coordinates
          x1 = DvertexCoords(1, i1)
          y1 = DvertexCoords(2, i1)
          x2 = DvertexCoords(1, i2)
          y2 = DvertexCoords(2, i2)
          x3 = DvertexCoords(1, i3)
          y3 = DvertexCoords(2, i3)

          ! Determine solution values at vertices
          u1 = Ddata(i1)
          u2 = Ddata(i2)
          u3 = Ddata(i3)


          ! Calculate determinant and area for triangle
          ddet  = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)
          darea = 0.5_DP*abs(ddet)

          ! Calculate averages of solution on triangle
          Dfunc(1)  = (abs(u1*(y2-y3)) + abs(u2*(y3-y1)) + abs(u3*(y1-y2))) / abs(ddet)
          Dfunc(2)  = (abs(u1*(x3-x2)) + abs(u2*(x1-x3)) + abs(u3*(x2-x1))) / abs(ddet)

          ! Calculate derivatives of solution on triangle
          Dderiv(1) = (u1*(y2-y3) + u2*(y3-y1) + u3*(y1-y2)) / ddet
          Dderiv(2) = (u1*(x3-x2) + u2*(x1-x3) + u3*(x2-x1)) / ddet

          ! Calculate absolute values of derivatives
          DabsDeriv = abs(Dderiv)

          ! Calculate derivatives of basis functions on triangle
          Dbas(1,1) = (y2-y3) / ddet
          Dbas(2,1) = (x3-x2) / ddet
          Dbas(1,2) = (y3-y1) / ddet
          Dbas(2,2) = (x1-x3) / ddet
          Dbas(1,3) = (y1-y2) / ddet
          Dbas(2,3) = (x2-x1) / ddet

          ! Calculate absolute values of basis functions
          DabsBas = abs(Dbas)

          ! Update nodal coefficient vector for node I1
          Dcoefficients(1,i1) = Dcoefficients(1,i1) + Dbas(1,1) * Dderiv(1) * darea
          Dcoefficients(2,i1) = Dcoefficients(2,i1) + Dbas(1,1) * Dderiv(2) * darea
          Dcoefficients(3,i1) = Dcoefficients(3,i1) + Dbas(2,1) * Dderiv(1) * darea
          Dcoefficients(4,i1) = Dcoefficients(4,i1) + Dbas(2,1) * Dderiv(2) * darea

          Dcoefficients(5,i1) = Dcoefficients(5,i1) + DabsBas(1,1) * DabsDeriv(1) * darea
          Dcoefficients(6,i1) = Dcoefficients(6,i1) + DabsBas(1,1) * DabsDeriv(2) * darea
          Dcoefficients(7,i1) = Dcoefficients(7,i1) + DabsBas(2,1) * DabsDeriv(1) * darea
          Dcoefficients(8,i1) = Dcoefficients(8,i1) + DabsBas(2,1) * DabsDeriv(2) * darea

          Dcoefficients(9, i1) = Dcoefficients(9, i1) + DabsBas(1,1) * dweight * Dfunc(1) * darea
          Dcoefficients(10,i1) = Dcoefficients(10,i1) + DabsBas(1,1) * dweight * Dfunc(2) * darea
          Dcoefficients(11,i1) = Dcoefficients(11,i1) + DabsBas(2,1) * dweight * Dfunc(1) * darea
          Dcoefficients(12,i1) = Dcoefficients(12,i1) + DabsBas(2,1) * dweight * Dfunc(2) * darea


          ! Update nodal coefficient vector for node I2
          Dcoefficients(1,i2) = Dcoefficients(1,i2) + Dbas(1,2) * Dderiv(1) * darea
          Dcoefficients(2,i2) = Dcoefficients(2,i2) + Dbas(1,2) * Dderiv(2) * darea
          Dcoefficients(3,i2) = Dcoefficients(3,i2) + Dbas(2,2) * Dderiv(1) * darea
          Dcoefficients(4,i2) = Dcoefficients(4,i2) + Dbas(2,2) * Dderiv(2) * darea

          Dcoefficients(5,i2) = Dcoefficients(5,i2) + DabsBas(1,2) * DabsDeriv(1) * darea
          Dcoefficients(6,i2) = Dcoefficients(6,i2) + DabsBas(1,2) * DabsDeriv(2) * darea
          Dcoefficients(7,i2) = Dcoefficients(7,i2) + DabsBas(2,2) * DabsDeriv(1) * darea
          Dcoefficients(8,i2) = Dcoefficients(8,i2) + DabsBas(2,2) * DabsDeriv(2) * darea

          Dcoefficients(9, i2) = Dcoefficients(9, i2) + DabsBas(1,2) * dweight * Dfunc(1) * darea
          Dcoefficients(10,i2) = Dcoefficients(10,i2) + DabsBas(1,2) * dweight * Dfunc(2) * darea
          Dcoefficients(11,i2) = Dcoefficients(11,i2) + DabsBas(2,2) * dweight * Dfunc(1) * darea
          Dcoefficients(12,i2) = Dcoefficients(12,i2) + DabsBas(2,2) * dweight * Dfunc(2) * darea


          ! Update nodal coefficient vector for node I3
          Dcoefficients(1,i3) = Dcoefficients(1,i3) + Dbas(1,3) * Dderiv(1) * darea
          Dcoefficients(2,i3) = Dcoefficients(2,i3) + Dbas(1,3) * Dderiv(2) * darea
          Dcoefficients(3,i3) = Dcoefficients(3,i3) + Dbas(2,3) * Dderiv(1) * darea
          Dcoefficients(4,i3) = Dcoefficients(4,i3) + Dbas(2,3) * Dderiv(2) * darea

          Dcoefficients(5,i3) = Dcoefficients(5,i3) + DabsBas(1,3) * DabsDeriv(1) * darea
          Dcoefficients(6,i3) = Dcoefficients(6,i3) + DabsBas(1,3) * DabsDeriv(2) * darea
          Dcoefficients(7,i3) = Dcoefficients(7,i3) + DabsBas(2,3) * DabsDeriv(1) * darea
          Dcoefficients(8,i3) = Dcoefficients(8,i3) + DabsBas(2,3) * DabsDeriv(2) * darea

          Dcoefficients(9, i3) = Dcoefficients(9, i3) + DabsBas(1,3) * dweight * Dfunc(1) * darea
          Dcoefficients(10,i3) = Dcoefficients(10,i3) + DabsBas(1,3) * dweight * Dfunc(2) * darea
          Dcoefficients(11,i3) = Dcoefficients(11,i3) + DabsBas(2,3) * dweight * Dfunc(1) * darea
          Dcoefficients(12,i3) = Dcoefficients(12,i3) + DabsBas(2,3) * dweight * Dfunc(2) * darea


        case (TRIA_NVEQUAD2D)

          ! Each quadrilateral 1-2-3-4 is subdivided into two triangles 1-2-3 and 1-3-4

          ! Determine global degrees of freedom
          i1 = IverticesAtElement(1, iel)
          i2 = IverticesAtElement(2, iel)
          i3 = IverticesAtElement(3, iel)
          i4 = IverticesAtElement(4, iel)

          ! Determine vertex coordinates
          x1 = DvertexCoords(1, i1)
          y1 = DvertexCoords(2, i1)
          x2 = DvertexCoords(1, i2)
          y2 = DvertexCoords(2, i2)
          x3 = DvertexCoords(1, i3)
          y3 = DvertexCoords(2, i3)
          x4 = DvertexCoords(1, i4)
          y4 = DvertexCoords(2, i4)

          ! Determine solution values at vertices
          u1 = Ddata(i1)
          u2 = Ddata(i2)
          u3 = Ddata(i3)
          u4 = Ddata(i4)


          ! Calculate determinant and area for triangle (1-2-3)
          ddet  = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)
          darea = 0.5_DP*abs(ddet)

          ! Calculate averages of solution on triangle
          Dfunc(1)  = (abs(u1*(y2-y3)) + abs(u2*(y3-y1)) + abs(u3*(y1-y2))) / abs(ddet)
          Dfunc(2)  = (abs(u1*(x3-x2)) + abs(u2*(x1-x3)) + abs(u3*(x2-x1))) / abs(ddet)

          ! Calculate derivatives of solution on triangle
          Dderiv(1) = (u1*(y2-y3) + u2*(y3-y1) + u3*(y1-y2)) / ddet
          Dderiv(2) = (u1*(x3-x2) + u2*(x1-x3) + u3*(x2-x1)) / ddet

          ! Calculate absolute values of derivatives
          DabsDeriv = abs(Dderiv)

          ! Calculate derivatives of basis functions on triangle
          Dbas(1,1) = (y2-y3) / ddet
          Dbas(2,1) = (x3-x2) / ddet
          Dbas(1,2) = (y3-y1) / ddet
          Dbas(2,2) = (x1-x3) / ddet
          Dbas(1,3) = (y1-y2) / ddet
          Dbas(2,3) = (x2-x1) / ddet

          ! Calculate absolute values of basis functions
          DabsBas = abs(Dbas)


          ! Update nodal coefficient vector for node I1
          Dcoefficients(1,i1) = Dcoefficients(1,i1) + Dbas(1,1) * Dderiv(1) * darea
          Dcoefficients(2,i1) = Dcoefficients(2,i1) + Dbas(1,1) * Dderiv(2) * darea
          Dcoefficients(3,i1) = Dcoefficients(3,i1) + Dbas(2,1) * Dderiv(1) * darea
          Dcoefficients(4,i1) = Dcoefficients(4,i1) + Dbas(2,1) * Dderiv(2) * darea

          Dcoefficients(5,i1) = Dcoefficients(5,i1) + DabsBas(1,1) * DabsDeriv(1) * darea
          Dcoefficients(6,i1) = Dcoefficients(6,i1) + DabsBas(1,1) * DabsDeriv(2) * darea
          Dcoefficients(7,i1) = Dcoefficients(7,i1) + DabsBas(2,1) * DabsDeriv(1) * darea
          Dcoefficients(8,i1) = Dcoefficients(8,i1) + DabsBas(2,1) * DabsDeriv(2) * darea

          Dcoefficients(9, i1) = Dcoefficients(9, i1) + DabsBas(1,1) * dweight * Dfunc(1) * darea
          Dcoefficients(10,i1) = Dcoefficients(10,i1) + DabsBas(1,1) * dweight * Dfunc(2) * darea
          Dcoefficients(11,i1) = Dcoefficients(11,i1) + DabsBas(2,1) * dweight * Dfunc(1) * darea
          Dcoefficients(12,i1) = Dcoefficients(12,i1) + DabsBas(2,1) * dweight * Dfunc(2) * darea


          ! Update nodal coefficient vector for node I2
          Dcoefficients(1,i2) = Dcoefficients(1,i2) + Dbas(1,2) * Dderiv(1) * darea
          Dcoefficients(2,i2) = Dcoefficients(2,i2) + Dbas(1,2) * Dderiv(2) * darea
          Dcoefficients(3,i2) = Dcoefficients(3,i2) + Dbas(2,2) * Dderiv(1) * darea
          Dcoefficients(4,i2) = Dcoefficients(4,i2) + Dbas(2,2) * Dderiv(2) * darea

          Dcoefficients(5,i2) = Dcoefficients(5,i2) + DabsBas(1,2) * DabsDeriv(1) * darea
          Dcoefficients(6,i2) = Dcoefficients(6,i2) + DabsBas(1,2) * DabsDeriv(2) * darea
          Dcoefficients(7,i2) = Dcoefficients(7,i2) + DabsBas(2,2) * DabsDeriv(1) * darea
          Dcoefficients(8,i2) = Dcoefficients(8,i2) + DabsBas(2,2) * DabsDeriv(2) * darea

          Dcoefficients(9, i2) = Dcoefficients(9, i2) + DabsBas(1,2) * dweight * Dfunc(1) * darea
          Dcoefficients(10,i2) = Dcoefficients(10,i2) + DabsBas(1,2) * dweight * Dfunc(2) * darea
          Dcoefficients(11,i2) = Dcoefficients(11,i2) + DabsBas(2,2) * dweight * Dfunc(1) * darea
          Dcoefficients(12,i2) = Dcoefficients(12,i2) + DabsBas(2,2) * dweight * Dfunc(2) * darea


          ! Update nodal coefficient vector for node I3
          Dcoefficients(1,i3) = Dcoefficients(1,i3) + Dbas(1,3) * Dderiv(1) * darea
          Dcoefficients(2,i3) = Dcoefficients(2,i3) + Dbas(1,3) * Dderiv(2) * darea
          Dcoefficients(3,i3) = Dcoefficients(3,i3) + Dbas(2,3) * Dderiv(1) * darea
          Dcoefficients(4,i3) = Dcoefficients(4,i3) + Dbas(2,3) * Dderiv(2) * darea

          Dcoefficients(5,i3) = Dcoefficients(5,i3) + DabsBas(1,3) * DabsDeriv(1) * darea
          Dcoefficients(6,i3) = Dcoefficients(6,i3) + DabsBas(1,3) * DabsDeriv(2) * darea
          Dcoefficients(7,i3) = Dcoefficients(7,i3) + DabsBas(2,3) * DabsDeriv(1) * darea
          Dcoefficients(8,i3) = Dcoefficients(8,i3) + DabsBas(2,3) * DabsDeriv(2) * darea

          Dcoefficients(9, i3) = Dcoefficients(9, i3) + DabsBas(1,3) * dweight * Dfunc(1) * darea
          Dcoefficients(10,i3) = Dcoefficients(10,i3) + DabsBas(1,3) * dweight * Dfunc(2) * darea
          Dcoefficients(11,i3) = Dcoefficients(11,i3) + DabsBas(2,3) * dweight * Dfunc(1) * darea
          Dcoefficients(12,i3) = Dcoefficients(12,i3) + DabsBas(2,3) * dweight * Dfunc(2) * darea


          ! Calculate determinant and area for triangle (1-3-4)
          ddet  = (x3-x1)*(y4-y1) - (x4-x1)*(y3-y1)
          darea = 0.5_DP*abs(ddet)

          ! Calculate averages of solution on triangle
          Dfunc(1)  = (abs(u1*(y3-y4)) + abs(u3*(y4-y1)) + abs(u4*(y1-y3))) / abs(ddet)
          Dfunc(2)  = (abs(u1*(x4-x3)) + abs(u3*(x1-x4)) + abs(u4*(x3-x1))) / abs(ddet)

          ! Calculate derivatives of solution on triangle
          Dderiv(1) = (u1*(y3-y4) + u3*(y4-y1) + u4*(y1-y3)) / ddet
          Dderiv(2) = (u1*(x4-x3) + u3*(x1-x4) + u4*(x3-x1)) / ddet

          ! Calculate absolute values of derivatives
          DabsDeriv = abs(Dderiv)

          ! Calculate derivatives of basis functions on triangle
          Dbas(1,1) = (y3-y4) / ddet
          Dbas(2,1) = (x4-x3) / ddet
          Dbas(1,2) = (y4-y1) / ddet
          Dbas(2,2) = (x1-x4) / ddet
          Dbas(1,3) = (y1-y3) / ddet
          Dbas(2,3) = (x3-x1) / ddet

          ! Calculate absolute values of basis functions
          DabsBas = abs(Dbas)

          ! Update nodal coefficient vector for node I1
          Dcoefficients(1,i1) = Dcoefficients(1,i1) + Dbas(1,1) * Dderiv(1) * darea
          Dcoefficients(2,i1) = Dcoefficients(2,i1) + Dbas(1,1) * Dderiv(2) * darea
          Dcoefficients(3,i1) = Dcoefficients(3,i1) + Dbas(2,1) * Dderiv(1) * darea
          Dcoefficients(4,i1) = Dcoefficients(4,i1) + Dbas(2,1) * Dderiv(2) * darea

          Dcoefficients(5,i1) = Dcoefficients(5,i1) + DabsBas(1,1) * DabsDeriv(1) * darea
          Dcoefficients(6,i1) = Dcoefficients(6,i1) + DabsBas(1,1) * DabsDeriv(2) * darea
          Dcoefficients(7,i1) = Dcoefficients(7,i1) + DabsBas(2,1) * DabsDeriv(1) * darea
          Dcoefficients(8,i1) = Dcoefficients(8,i1) + DabsBas(2,1) * DabsDeriv(2) * darea

          Dcoefficients(9, i1) = Dcoefficients(9, i1) + DabsBas(1,1) * dweight * Dfunc(1) * darea
          Dcoefficients(10,i1) = Dcoefficients(10,i1) + DabsBas(1,1) * dweight * Dfunc(2) * darea
          Dcoefficients(11,i1) = Dcoefficients(11,i1) + DabsBas(2,1) * dweight * Dfunc(1) * darea
          Dcoefficients(12,i1) = Dcoefficients(12,i1) + DabsBas(2,1) * dweight * Dfunc(2) * darea


          ! Update nodal coefficient vector for node I3
          Dcoefficients(1,i3) = Dcoefficients(1,i3) + Dbas(1,2) * Dderiv(1) * darea
          Dcoefficients(2,i3) = Dcoefficients(2,i3) + Dbas(1,2) * Dderiv(2) * darea
          Dcoefficients(3,i3) = Dcoefficients(3,i3) + Dbas(2,2) * Dderiv(1) * darea
          Dcoefficients(4,i3) = Dcoefficients(4,i3) + Dbas(2,2) * Dderiv(2) * darea

          Dcoefficients(5,i3) = Dcoefficients(5,i3) + DabsBas(1,2) * DabsDeriv(1) * darea
          Dcoefficients(6,i3) = Dcoefficients(6,i3) + DabsBas(1,2) * DabsDeriv(2) * darea
          Dcoefficients(7,i3) = Dcoefficients(7,i3) + DabsBas(2,2) * DabsDeriv(1) * darea
          Dcoefficients(8,i3) = Dcoefficients(8,i3) + DabsBas(2,2) * DabsDeriv(2) * darea

          Dcoefficients(9, i3) = Dcoefficients(9, i3) + DabsBas(1,2) * dweight * Dfunc(1) * darea
          Dcoefficients(10,i3) = Dcoefficients(10,i3) + DabsBas(1,2) * dweight * Dfunc(2) * darea
          Dcoefficients(11,i3) = Dcoefficients(11,i3) + DabsBas(2,2) * dweight * Dfunc(1) * darea
          Dcoefficients(12,i3) = Dcoefficients(12,i3) + DabsBas(2,2) * dweight * Dfunc(2) * darea


          ! Update nodal coefficient vector for node I4
          Dcoefficients(1,i4) = Dcoefficients(1,i4) + Dbas(1,3) * Dderiv(1) * darea
          Dcoefficients(2,i4) = Dcoefficients(2,i4) + Dbas(1,3) * Dderiv(2) * darea
          Dcoefficients(3,i4) = Dcoefficients(3,i4) + Dbas(2,3) * Dderiv(1) * darea
          Dcoefficients(4,i4) = Dcoefficients(4,i4) + Dbas(2,3) * Dderiv(2) * darea

          Dcoefficients(5,i4) = Dcoefficients(5,i4) + DabsBas(1,3) * DabsDeriv(1) * darea
          Dcoefficients(6,i4) = Dcoefficients(6,i4) + DabsBas(1,3) * DabsDeriv(2) * darea
          Dcoefficients(7,i4) = Dcoefficients(7,i4) + DabsBas(2,3) * DabsDeriv(1) * darea
          Dcoefficients(8,i4) = Dcoefficients(8,i4) + DabsBas(2,3) * DabsDeriv(2) * darea

          Dcoefficients(9, i4) = Dcoefficients(9, i4) + DabsBas(1,3) * dweight * Dfunc(1) * darea
          Dcoefficients(10,i4) = Dcoefficients(10,i4) + DabsBas(1,3) * dweight * Dfunc(2) * darea
          Dcoefficients(11,i4) = Dcoefficients(11,i4) + DabsBas(2,3) * dweight * Dfunc(1) * darea
          Dcoefficients(12,i4) = Dcoefficients(12,i4) + DabsBas(2,3) * dweight * Dfunc(2) * darea


        case DEFAULT
          call output_line('Unsupproted element type!',&
              OU_CLASS_ERROR,OU_MODE_STD,'doSecondDiffIndicator2D')
          call sys_halt()
        end select
      end do

      ! Loop over all vertices in triangulation
      do ivt = 1, NVT
        ! Compute numerator
        daux1 = Dcoefficients(1, ivt)*Dcoefficients(1, ivt) +&
                Dcoefficients(2, ivt)*Dcoefficients(2, ivt) +&
                Dcoefficients(3, ivt)*Dcoefficients(3, ivt) +&
                Dcoefficients(4, ivt)*Dcoefficients(4, ivt)

        ! Compute derivative part of denominator
        daux2 = Dcoefficients(5, ivt)*Dcoefficients(5, ivt) +&
                Dcoefficients(6, ivt)*Dcoefficients(6, ivt) +&
                Dcoefficients(7, ivt)*Dcoefficients(7, ivt) +&
                Dcoefficients(8, ivt)*Dcoefficients(8, ivt)

        ! Compute average part of denominator
        daux3 = Dcoefficients( 9, ivt)*Dcoefficients( 9, ivt) +&
                Dcoefficients(10, ivt)*Dcoefficients(10, ivt) +&
                Dcoefficients(11, ivt)*Dcoefficients(11, ivt) +&
                Dcoefficients(12, ivt)*Dcoefficients(12, ivt)

        ! Compute nodal indicator
        DnodalIndicator(ivt) = sqrt(daux1/(daux2+max(dfilter, daux3)))
      end do


      ! Clear array
      call lalg_clearVectorDble(Dindicator)

      ! Loop over all elements in triangulation
      do iel = 1, NEL

        ! Initialise number of contributions
        ncontributions = 0

        ! Loop over all vertices of the element
        do ive = 1, tria_getNVE(IverticesAtElement, iel)

          ! Determine global degree of freedom
          ivt = IverticesAtElement(ive, iel)

          ! Skip boundary nodes
          if (InodalProperty(ivt) .ne. 0) cycle

          ! Apply nodal contribution to element indicator
          Dindicator(iel) = Dindicator(iel) + DnodalIndicator(ivt)

          ! Increase number of contributions
          ncontributions = ncontributions+1
        end do

        ! Scale element indicator by number of contributions
        if (ncontributions .ne. 0) then
          Dindicator(iel) = Dindicator(iel) / real(ncontributions, DP)
        else
          Dindicator(iel) = 0.0_DP
        end if
      end do

    end subroutine doSecondDiffIndicator2D
  end subroutine ppind_secondDifference

end module pprocindicator
