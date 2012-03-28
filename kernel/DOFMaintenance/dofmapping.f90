!##############################################################################
!# ****************************************************************************
!# <name> dofmapping </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module is a basic module of the discretisation. It contains routines
!# to map local degrees of freedom on one element primitive onto the global
!# degrees of freedom of a solution vector.
!#
!# The following functions provide support for global DOF`s:
!#
!# 1.) dof_igetNDofGlob
!#     -> Get number of global DOF`s described by a given discretisation
!#
!# 2.) dof_igetNDofGlobBlock
!#     -> Get number of global DOF`s described by a given block discretisation
!#
!# 3.) dof_locGlobMapping
!#     -> Map the 'local' degrees of freedom 1..n on one element to the global
!#        degrees of freedom according to a discretisaion
!#
!# 4.) dof_locGlobMapping_mult
!#     -> Map the 'local' degrees of freedom 1..n on a set of elements to
!#        the global degrees of freedom according to a discretisaion.
!#
!# 5.) dof_infoDiscr
!#     -> Prints out information about a discretisation to the terminal
!#
!# 6.) dof_infoDiscrBlock
!#     -> Prints out information about a block discretisation to the terminal
!#
!# 7.) dof_precomputeDofMapping
!#     -> Precompute the DOF-mapping into a discretisation structure
!#
!# </purpose>
!##############################################################################

module dofmapping

!$use omp_lib
  use basicgeometry
  use element
  use fsystem
  use genoutput
  use spatialdiscretisation
  use storage
  use triangulation

  implicit none

  private

!<constants>

  !<constantblock description="Kind values for global DOF`s">

  ! kind value for indexing global DOF`s
  ! !!! DEPRECATED: DO NOT USE THIS CONSTANT ANYMORE !!!
  integer, parameter, public :: PREC_DOFIDX     = I32

  !</constantblock>

!</constants>

  public :: dof_igetNDofGlob
  public :: dof_igetNDofGlobBlock
  public :: dof_locGlobMapping
  public :: dof_locGlobMapping_mult
  public :: dof_infoDiscr
  public :: dof_infoDiscrBlock
  public :: dof_precomputeDofMapping

contains

  ! ***************************************************************************

!<function>

  integer function dof_igetNDofGlob(rdiscretisation)

!<description>
  ! This function returns for a given discretisation the number of global
  ! degrees of freedom in the corresponding scalar DOF vector.
!</description>

!<input>
  ! The discretisation structure that specifies the (scalar) discretisation.
  type(t_spatialDiscretisation), intent(in) :: rdiscretisation
!</input>

!<result>
  ! global number of equations on current grid
!</result>

!</function>

  integer(I32) :: celement
  integer(I32), dimension(2) :: Celements

  dof_igetNDofGlob = 0

  if (rdiscretisation%bprecompiledDofMapping) then
    ! The value is already available.
    dof_igetNDofGlob = rdiscretisation%ndof
    return
  end if

  select case(rdiscretisation%ndimension)
  case (NDIM1D)

      celement = rdiscretisation%RelementDistr(1)%celement
      dof_igetNDofGlob = NDFG_uniform1D (rdiscretisation%p_rtriangulation, celement)

  case (NDIM2D)
    if (rdiscretisation%ccomplexity .eq. SPDISC_UNIFORM) then

      celement = rdiscretisation%RelementDistr(1)%celement
      ! Uniform discretisation - fall back to the old FEAT mapping
      dof_igetNDofGlob = NDFG_uniform2D (rdiscretisation%p_rtriangulation, celement)

    else if (rdiscretisation%ccomplexity .eq. SPDISC_CONFORMAL) then

      ! Conformal discretisation. That is a little bit tricky!
      ! At first, we support only the case where two element types are mixed.
      if (rdiscretisation%inumFESpaces .eq. 2) then
        Celements(1) = rdiscretisation%RelementDistr(1)%celement
        Celements(2) = rdiscretisation%RelementDistr(2)%celement

        dof_igetNDofGlob = NDFG_conformal2D_2el (&
            rdiscretisation%p_rtriangulation, Celements(1:2))

      end if

    end if

  case (NDIM3D)
    ! Currently, only uniform discretisations are supported.
    if (rdiscretisation%ccomplexity .eq. SPDISC_UNIFORM) then

      celement = rdiscretisation%RelementDistr(1)%celement

      ! Uniform discretisation - fall back to the old FEAT mapping
      dof_igetNDofGlob = NDFG_uniform3D (rdiscretisation%p_rtriangulation, celement)

    end if

  case default

    ! Dimension not supported
    call output_line('Invalid discretisation dimension!',&
                     OU_CLASS_ERROR,OU_MODE_STD,'dof_igetNDofGlob')
    call sys_halt()

  end select

  contains

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Internal subroutine: Get global DOF number for uniform discretisation
    ! with only one element type.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    integer function NDFG_uniform1D (rtriangulation, celement)

    ! IN: The underlying triangulation
    type(t_triangulation), intent(in) :: rtriangulation

    ! IN: The element type of the discretisation
    integer(I32), intent(in) :: celement

    ! OUT: number of global DOF`s.

    ! The number of global DOF`s depends on the element type...
    select case (elem_getPrimaryElement(celement))
    case (EL_P0_1D)
      ! DOF`s in the cell midpoints
      NDFG_uniform1D = rtriangulation%NEL
    case (EL_P1_1D)
      ! DOF`s in the vertices
      NDFG_uniform1D = rtriangulation%NVT
    case (EL_P2_1D)
      ! DOF`s in the vertices and cell midpoints
      NDFG_uniform1D = rtriangulation%NVT + rtriangulation%NEL
    case (EL_S31_1D)
      ! DOF`s in the vertices
      NDFG_uniform1D = 2*rtriangulation%NVT
    case (EL_PN_1D)
      ! DOF`s in the vertices + cells
      NDFG_uniform1D = rtriangulation%NVT + rtriangulation%NEL * &
                       iand(ishft(celement,-16),255_I32)
    case (EL_DG_T0_1D)
      ! DOF`s in the cell midpoints
      NDFG_uniform1D = rtriangulation%NEL
    case (EL_DG_T1_1D)
      ! DOF`s in the cell midpoints
      NDFG_uniform1D = 2*rtriangulation%NEL
    case (EL_DG_T2_1D)
      ! DOF`s in the cell midpoints
      NDFG_uniform1D = 3*rtriangulation%NEL
    end select

    end function

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Internal subroutine: Get global DOF number for uniform discretisation
    ! with only one element type.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! This is roughly the NDFG routine of the old FEAT library...

    integer function NDFG_uniform2D (rtriangulation, celement)

    ! IN: The underlying triangulation
    type(t_triangulation), intent(in) :: rtriangulation

    ! IN: The element type of the discretisation
    integer(I32), intent(in) :: celement

    ! OUT: number of global DOF`s.

    ! The number of global DOF`s depends on the element type...
    select case (elem_getPrimaryElement(celement))
    case (EL_P0, EL_Q0)
      ! DOF`s in the cell midpoints
      NDFG_uniform2D = rtriangulation%NEL
    case (EL_P1, EL_Q1)
      ! DOF`s in the vertices
      NDFG_uniform2D = rtriangulation%NVT
    case (EL_DG_P1_2D)
      ! DOF`s in the elements
      NDFG_uniform2D = rtriangulation%NEL*3
    case (EL_DG_Q1_2D)
      ! DOF`s in the elements
      NDFG_uniform2D = rtriangulation%NEL*4
    case (EL_DG_Q2_2D)
      ! Each element has 9 DOFS
      ! in the vertices, edge midpoints and element midpoints
      NDFG_uniform2D = rtriangulation%NEL*9
    case (EL_P2)
      ! DOF`s in the vertices and edge midpoints
      NDFG_uniform2D = rtriangulation%NVT + rtriangulation%NMT
    case (EL_Q2)
      ! DOF`s in the vertices, edge midpoints and element midpoints
      NDFG_uniform2D = rtriangulation%NVT + rtriangulation%NMT + rtriangulation%NEL
    case (EL_P3)
      ! 1 DOF`s per vertices, 2 DOF per edge
      NDFG_uniform2D = rtriangulation%NVT + 2*rtriangulation%NMT
    case (EL_Q3)
      ! 1 DOF`s per vertices, 2 DOF per edge, 4 DOF in the inner
      NDFG_uniform2D = rtriangulation%NVT + 2*rtriangulation%NMT + 4*rtriangulation%NEL
    case (EL_QP1)
      ! 3 DOF`s in the midpoint of the element.
      NDFG_uniform2D = 3*rtriangulation%NEL
    case (EL_QPW4P1_2D)
      ! 1 DOF in each vertex and one in the element midpoint
      NDFG_uniform2D = rtriangulation%NVT + rtriangulation%NEL
    case (EL_QPW4P2_2D)
      ! 1 DOF in each vertex, 1 DOF in each edge and 5 in the element
      NDFG_uniform2D = rtriangulation%NVT + rtriangulation%NMT + 5*rtriangulation%NEL
    case (EL_P1T, EL_Q1T)
      ! 1 DOF per edge
      NDFG_uniform2D = rtriangulation%NMT
    case (EL_Q1TB)
      ! 1 DOF per edge, one per element
      NDFG_uniform2D = rtriangulation%NMT + rtriangulation%NEL
    case (EL_Q2T)
      ! 2 DOF in the edges, one per element
      NDFG_uniform2D = 2*rtriangulation%NMT + rtriangulation%NEL
    case (EL_Q2TB)
      ! E037
      NDFG_uniform2D = 2*rtriangulation%NMT + 2*rtriangulation%NEL
    case (EL_Q3T_2D)
      ! 3 DOF`s in the edges, three per element
      NDFG_uniform2D = 3*(rtriangulation%NMT + rtriangulation%NEL)
    case (EL_DG_T0_2D)
      ! DOF`s in the cell midpoints
      NDFG_uniform2D = rtriangulation%NEL
    case (EL_DG_T1_2D, EL_DCTP1_2D, EL_DCQP1_2D)
      ! DOF`s in the cell midpoints
      NDFG_uniform2D = 3*rtriangulation%NEL
    case (EL_DG_T2_2D, EL_DCTP2_2D, EL_DCQP2_2D)
      ! DOF`s in the cell midpoints
      NDFG_uniform2D = 6*rtriangulation%NEL
    case (EL_QPW4DCP1_2D)
      ! DOF`s in the cell midpoints
      NDFG_uniform2D = 11*rtriangulation%NEL
    end select

    end function

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Internal subroutine: Get global DOF number for conformal discretisation
    ! with two element types.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    integer function NDFG_conformal2D_2el (rtriangulation, Celements)

    ! IN: The underlying triangulation
    type(t_triangulation), intent(in) :: rtriangulation

    ! IN: List of element types in the discretisation. Celements(1) is one element
    ! type identifier, Celements(2) the other one.
    integer(I32), dimension(:), intent(in) :: Celements

    ! OUT: number of global DOF`s.

    ! local variables
    integer(I32), dimension(size(Celements)) :: IelTypesPrimary

    ! Get the primary element number
    IelTypesPrimary = elem_getPrimaryElement(Celements)

    ! The number of global DOF`s depends on the element type...
    select case (IelTypesPrimary(1))
    case (EL_P0, EL_Q0)
      select case (IelTypesPrimary(2))
      case (EL_P0, EL_Q0)
        ! DOF`s in the cell midpoints
        NDFG_conformal2D_2el = rtriangulation%NEL
      end select

    case (EL_P1, EL_Q1)
      select case (IelTypesPrimary(2))
      case (EL_P1, EL_Q1)
        ! DOF`s in the vertices
        NDFG_conformal2D_2el = rtriangulation%NVT
      end select

    case (EL_P2,EL_Q2)
      select case (IelTypesPrimary(2))
      case (EL_Q2)
        ! Number of vertices + Number of edges (edge midpoints) +
        ! Number of quads (quad midpoints)
        NDFG_conformal2D_2el = rtriangulation%NVT + rtriangulation%NMT + &
            rtriangulation%InelOfType(TRIA_NVEQUAD2D)
      end select

    case (EL_P1T,EL_Q1T)
      select case (IelTypesPrimary(2))
      case (EL_P1T,EL_Q1T)
        ! DOF`s in the edge midpoints
        NDFG_conformal2D_2el = rtriangulation%NMT

      end select

    end select

    end function

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Internal subroutine: Get global DOF number for uniform discretisation
    ! with only one element type.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! This is roughly the NDFG routine of the old FEAT library...

    integer function NDFG_uniform3D (rtriangulation, celement)

    ! IN: The underlying triangulation
    type(t_triangulation), intent(in) :: rtriangulation

    ! IN: The element type of the discretisation
    integer(I32), intent(in) :: celement

    ! OUT: number of global DOF`s.

    ! The number of global DOF`s depends on the element type...
    select case (elem_getPrimaryElement(celement))
    case (EL_P0_3D, EL_Q0_3D, EL_Y0_3D, EL_R0_3D)
      ! DOF`s in the cell midpoints
      NDFG_uniform3D = rtriangulation%NEL
    case (EL_P1_3D, EL_Q1_3D, EL_Y1_3D, EL_R1_3D)
      ! DOF`s in the vertices
      NDFG_uniform3D = rtriangulation%NVT
    case (EL_Q2_3D)
      ! DOF`s in everything
      NDFG_uniform3D = rtriangulation%NVT + rtriangulation%NMT + &
                       rtriangulation%NAT + rtriangulation%NEL
    case (EL_QP1_3D)
      ! 4 DOF`s in the midpoint of the element.
      NDFG_uniform3D = 4*rtriangulation%NEL
    case (EL_Q1T_3D)
      ! DOF`s in the face midpoints
      NDFG_uniform3D = rtriangulation%NAT
    case (EL_Q2T_3D)
      ! DOF`s in the face midpoints
      NDFG_uniform3D = 3*rtriangulation%NAT + rtriangulation%NEL
    end select

    end function

  end function

  ! ***************************************************************************

!<function>

  integer function dof_igetNDofGlobBlock(rdiscretisation)

!<description>
  ! This function returns for a given block discretisation the number of global
  ! degrees of freedom in the corresponding block DOF vector.
!</description>

!<input>
  ! The discretisation structure that specifies the (block) discretisation.
  type(t_blockDiscretisation), intent(in) :: rdiscretisation
!</input>

!<result>
  ! global number of equations on current grid
!</result>

!</function>

    ! Sum up the DOF`s of every scalar sub-discretisation structure
    integer :: i
    integer :: icount

    icount = 0
    do i=1,rdiscretisation%ncomponents
      icount = icount + &
          dof_igetNDofGlob(rdiscretisation%RspatialDiscr(i))
    end do

    dof_igetNDofGlobBlock = icount

  end function

  ! ***************************************************************************

!<subroutine>

  subroutine dof_locGlobMapping(rdiscretisation, ielIdx, IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the element ielIdx. It is a wrapper routine
  ! for the corresponding routines for a specific element type.
  !
  ! On each element, there are a number of local DOF`s 1..n (where n can
  ! be obtained using elem_igetNDofLoc). This subroutine returns the
  ! corresponding global degrees of freedom, corresponding to these local
  ! ones.
!</description>

!<input>

  ! The discretisation structure that specifies the (scalar) discretisation.
  type(t_spatialDiscretisation), intent(in) :: rdiscretisation

  ! Element index, where the mapping should be computed.
  integer, intent(in) :: ielIdx

!</input>

!<output>

  ! array of global DOF numbers
  integer, dimension(:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer, dimension(1) :: ielIdx_array
  integer, dimension(size(IdofGlob),1) :: IdofGlob_array

    ! Wrapper to dof_locGlobMapping_mult - better use this directly!
    ielIdx_array(1) = ielidx
    IdofGlob_array(:,1) = IdofGlob(:)
    call dof_locGlobMapping_mult (rdiscretisation, ielIdx_array,  &
                                  IdofGlob_array)
    IdofGlob = IdofGlob_array(:,1)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine dof_locGlobMapping_mult(rdiscretisation, IelIdx, IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements IelIdx. It is a wrapper routine
  ! for the corresponding routines for a specific element type.
  !
  ! On each element, there are a number of local DOF`s 1..n (where n can
  ! be obtained using elem_igetNDofLoc). This subroutine returns the
  ! corresponding global degrees of freedom, corresponding to these local
  ! ones.
  ! The routine allows to calculate the DOF`s to a list of elements
  ! which are all of the same type. IelIdx is this list of elements and
  ! IdofGlob is a 2D array which receives for every element the
  ! corresponding global DOF`s.
!</description>

!<input>

  ! The discretisation structure that specifies the (scalar) discretisation.
  type(t_spatialDiscretisation), intent(in) :: rdiscretisation

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

    ! local variables
    integer, dimension(:,:), pointer :: p_2darray,p_2darray2,p_2darray3
    integer, dimension(:), pointer :: p_IelementDofs,p_IelementDofIdx
    integer, dimension(:), pointer :: p_IelementCounter
    type(t_triangulation), pointer :: p_rtriangulation
    integer(I32) :: celement
    integer(I32), dimension(2) :: Celements
    integer :: i,iel,isize,ipos,ielnr

    if (rdiscretisation%bprecompiledDofMapping) then
      ! Just copy the global DOF`s from the array with the
      ! DOF`s of the element. All elements have the same number of DOF`s.
      call storage_getbase_int (rdiscretisation%h_IelementDofs,p_IelementDofs)
      call storage_getbase_int (rdiscretisation%h_IelementDofIdx,p_IelementDofIdx)

      do iel = 1,size(IelIdx)
        ielnr = IelIdx(iel)
        isize = p_IelementDofIdx(ielnr+1) - p_IelementDofIdx(ielnr)
        ipos = p_IelementDofIdx(ielnr)
        do i=1,isize
          IdofGlob(i,iel) = p_IelementDofs(ipos+i-1)
        end do
      end do

      ! That is it.
      return
    end if

    ! No, we have to compute...
    p_rtriangulation => rdiscretisation%p_rtriangulation

    select case(rdiscretisation%ndimension)
    case (NDIM1D)
      ! Call the right 'multiple-get' routines for global DOF`s.
      ! For this purpose we evaluate the pointers in the discretisation
      ! structure (if necessary) to prevent another call using pointers...
      ! The number of global DOF`s depends on the element type...

      celement = rdiscretisation%RelementDistr(1)%celement

      select case (elem_getPrimaryElement(celement))
      case (EL_P0_1D)
        ! DOF`s for P0
        call dof_locGlobUniMult_P0_1D(IelIdx, IdofGlob)
        return
      case (EL_P1_1D)
        ! DOF`s in the vertices
        call storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
        call dof_locGlobUniMult_P1_1D(p_2darray, IelIdx, IdofGlob)
        return
      case (EL_P2_1D)
        ! DOF`s in the vertices and line midpoints
        call storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
        call dof_locGlobUniMult_P2_1D(p_rtriangulation%NVT, p_2darray, IelIdx,&
                                      IdofGlob)
        return
      case (EL_S31_1D)
        ! DOF`s in the vertices
        call storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
        call dof_locGlobUniMult_S31_1D(p_rtriangulation%NVT, p_2darray, IelIdx,&
                                       IdofGlob)
        return
      case (EL_PN_1D)
        ! DOF`s in the vertices
        call storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
        call dof_locGlobUniMult_PN_1D(iand(int(ishft(celement,-16)),255), &
            p_rtriangulation%NVT, p_rtriangulation%NEL, p_2darray, IelIdx, IdofGlob)
        return
      case (EL_DG_T0_1D)
        ! DOF`s for DG_T0_1D
        call dof_locGlobUniMult_DG_T0_1D(IelIdx, IdofGlob)
        return
      case (EL_DG_T1_1D)
        ! DOF`s for DG_T1_1D
        call dof_locGlobUniMult_DG_T1_1D(IelIdx, IdofGlob)
        return
      case (EL_DG_T2_1D)
        ! DOF`s for DG_T2_1D
        call dof_locGlobUniMult_DG_T2_1D(IelIdx, IdofGlob)
        return
      end select

    case (NDIM2D)
      ! At first we deal only with uniform discretisations
      if (rdiscretisation%ccomplexity .eq. SPDISC_UNIFORM) then

        ! Call the right 'multiple-get' routines for global DOF`s.
        ! For this purpose we evaluate the pointers in the discretisation
        ! structure (if necessary) to prevent another call using pointers...
        ! The number of global DOF`s depends on the element type...

        celement = rdiscretisation%RelementDistr(1)%celement

        select case (elem_getPrimaryElement(celement))
        case (EL_P0, EL_Q0)
          ! DOF`s for Q0
          call dof_locGlobUniMult_P0Q0(IelIdx, IdofGlob)
          return
        case (EL_P1, EL_Q1)
          ! DOF`s in the vertices
          call storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
          call dof_locGlobUniMult_P1Q1(p_2darray, IelIdx, IdofGlob)
          return
        case (EL_DG_P1_2D)
          ! DOF`s in the vertices
          call dof_locGlobUniMult_DGP12D(IelIdx, IdofGlob)
          return
        case (EL_DG_Q1_2D)
          ! DOF`s in the vertices
          call dof_locGlobUniMult_DGQ12D(IelIdx, IdofGlob)
          return
        case (EL_DG_Q2_2D)
          ! DOF`s in the vertices, edge midpoints and element midpoint
          call dof_locGlobUniMult_DGQ22D(IelIdx, IdofGlob)
          return
        case (EL_P2)
          ! DOF`s in the vertices and egde midpoints
          call storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
          call storage_getbase_int2D (p_rtriangulation%h_IedgesAtElement,p_2darray2)
          call dof_locGlobUniMult_P2(p_rtriangulation%NVT,p_2darray, p_2darray2,&
                                     IelIdx, IdofGlob)
          return
        case (EL_Q2)
          ! DOF`s in the vertices, egde midpoints and element midpoints
          call storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
          call storage_getbase_int2D (p_rtriangulation%h_IedgesAtElement,p_2darray2)
          call dof_locGlobUniMult_Q2(p_rtriangulation%NVT,p_rtriangulation%NMT,&
                                     p_2darray, p_2darray2, IelIdx, IdofGlob)
          return
        case (EL_P3)
        case (EL_Q3)

        case (EL_QP1)
          ! DOF`s for Q1
          call dof_locGlobUniMult_QP1(p_rtriangulation%NEL,IelIdx, IdofGlob)
          return
        case (EL_P1T)
          ! DOF`s in the edges
          call storage_getbase_int2D (p_rtriangulation%h_IedgesAtElement,p_2darray)
          call dof_locGlobUniMult_P1T(p_2darray, IelIdx, IdofGlob)
          return
        case (EL_Q1T)
          ! DOF`s in the edges
          call storage_getbase_int2D (p_rtriangulation%h_IedgesAtElement,p_2darray)
          call dof_locGlobUniMult_Q1T(p_2darray, IelIdx, IdofGlob)
          return
        case (EL_QPW4P1_2D)
          ! DOF`s in the vertices and element midpoints
          call storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
          call dof_locGlobUniMult_QPW4P1(p_2darray, IelIdx, IdofGlob, p_rtriangulation%NVT)
          return
        case (EL_QPW4P2_2D)
          ! DOF`s in the vertices and element midpoints
          call storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
          call storage_getbase_int2D (p_rtriangulation%h_IedgesAtElement,p_2darray2)
          call dof_locGlobUniMult_QPW4P2(p_2darray, p_2darray2, IelIdx, IdofGlob, &
              p_rtriangulation%NVT,p_rtriangulation%NMT)
          return
        case (EL_Q1TB)
          ! DOF`s in the edges
          call storage_getbase_int2D (p_rtriangulation%h_IedgesAtElement,p_2darray)
          call dof_locGlobUniMult_Q1TB(p_rtriangulation%NMT,p_2darray, IelIdx, IdofGlob)
          return
        case (EL_Q2T)
          ! DOF`s in the edges and the element center
          call storage_getbase_int2D (p_rtriangulation%h_IedgesAtElement,p_2darray)
          call dof_locGlobUniMult_Q2T(p_rtriangulation%NMT,p_2darray, IelIdx, IdofGlob)
          return
        case (EL_Q2TB)
          ! DOF`s in the edges and the element center
          call storage_getbase_int2D (p_rtriangulation%h_IedgesAtElement,p_2darray)
          call dof_locGlobUniMult_Q2TB(p_rtriangulation%NMT,p_rtriangulation%NEL,&
                                       p_2darray, IelIdx, IdofGlob)
          return
        case (EL_Q3T_2D)
          ! DOF`s in the edges and the element center
          call storage_getbase_int2D (p_rtriangulation%h_IedgesAtElement,p_2darray)
          call dof_locGlobUniMult_Q3T_2D(p_rtriangulation%NMT,p_rtriangulation%NEL,&
                                         p_2darray, IelIdx, IdofGlob)
          return
        case (EL_DG_T0_2D)
          ! DOF`s for DG_T0_1D
          call dof_locGlobUniMult_DG_T0_2D(IelIdx, IdofGlob)
          return
        case (EL_DG_T1_2D, EL_DCTP1_2D, EL_DCQP1_2D)
          ! DOF`s for DG_T1_1D
          call dof_locGlobUniMult_DG_T1_2D(IelIdx, IdofGlob)
          return
        case (EL_DG_T2_2D, EL_DCTP2_2D, EL_DCQP2_2D)
          ! DOF`s for DG_T2_1D
          call dof_locGlobUniMult_DG_T2_2D(IelIdx, IdofGlob)
          return
        case (EL_QPW4DCP1_2D)
          call dof_locGlobUniMult_QPW4DCP1_2D(IelIdx, IdofGlob)
          return
        end select

        return

      else if (rdiscretisation%ccomplexity .eq. SPDISC_CONFORMAL) then

        ! Conformal discretisation. That is a little bit tricky!
        ! At first, we support only the case where two element types are mixed.
        if (rdiscretisation%inumFESpaces .eq. 2) then

          Celements(1) = rdiscretisation%RelementDistr(1)%celement
          Celements(2) = rdiscretisation%RelementDistr(2)%celement

          ! Get the primary element type(s)
          Celements = elem_getPrimaryElement(Celements)

          ! Now the actual mappings...
          select case (Celements(1))
          case (EL_P0, EL_Q0)
            select case (Celements(2))
            case (EL_P0, EL_Q0)
              ! DOF`s in the cell midpoints.
              ! That works like P0 elements.
              call dof_locGlobUniMult_P0Q0(IelIdx, IdofGlob)
              return
            end select

          case (EL_P1, EL_Q1)
            select case (Celements(2))
            case (EL_P1, EL_Q1)
              ! DOF`s in the vertices.
              ! That works like P1 elements.
              call storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
              call dof_locGlobUniMult_P1Q1(p_2darray, IelIdx, IdofGlob)
              return
            end select

          case (EL_P2, EL_Q2)
            select case (Celements(2))
            case (EL_P2, EL_Q2)
              ! DOF`s in the vertices, edges and element mitpoints of the quads.
              ! For this purpose, we need the element counter array that counts
              ! every quad element.
              call storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
              call storage_getbase_int2D (p_rtriangulation%h_IedgesAtElement,p_2darray2)
              call storage_getbase_int (rdiscretisation%h_IelementCounter,p_IelementCounter)

              ! Use p_IverticesAtElement evaluated at the first element in the element
              ! set to determine NVE. It is either 3 or 4 and valid for all elements
              ! in the current element set.
              if (ubound(p_2darray,1) .ge. 4) then
                if (p_2darray(4,IelIdx(1)) .eq. 0) then
                  call dof_locGlobUniMult_P2Q2(p_rtriangulation%NVT,p_rtriangulation%NMT,3,&
                      p_2darray, p_2darray2, p_IelementCounter, IelIdx, IdofGlob)
                else
                  call dof_locGlobUniMult_P2Q2(p_rtriangulation%NVT,p_rtriangulation%NMT,4,&
                      p_2darray, p_2darray2, p_IelementCounter, IelIdx, IdofGlob)
                end if
              else
                ! Pure triangular mesh
                call dof_locGlobUniMult_P2Q2(p_rtriangulation%NVT,p_rtriangulation%NMT,3,&
                    p_2darray, p_2darray2, p_IelementCounter, IelIdx, IdofGlob)
              end if
              return
            end select

          case (EL_P1T, EL_Q1T)
            select case (Celements(2))
            case (EL_P1T, EL_Q1T)
              ! DOF`s in the edges
              ! That works like P1 elements.
              call storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
              call dof_locGlobUniMult_P1Q1(p_2darray, IelIdx, IdofGlob)
              return
            end select

          end select

        end if

      end if

    case (NDIM3D)
      ! At first we deal only with uniform discretisations
      if (rdiscretisation%ccomplexity .eq. SPDISC_UNIFORM) then

        ! Call the right 'multiple-get' routines for global DOF`s.
        ! For this purpose we evaluate the pointers in the discretisation
        ! structure (if necessary) to prevent another call using pointers...
        ! The number of global DOF`s depends on the element type...
        celement = rdiscretisation%RelementDistr(1)%celement

        select case (elem_getPrimaryElement(celement))
        case (EL_P0_3D, EL_Q0_3D, EL_Y0_3D, EL_R0_3D)
          ! DOF`s for Q0
          call dof_locGlobUniMult_P0Q0_3D(IelIdx, IdofGlob)
          return
        case (EL_P1_3D, EL_Q1_3D, EL_Y1_3D, EL_R1_3D)
          ! DOF`s in the vertices
          call storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
          call dof_locGlobUniMult_P1Q1_3D(p_2darray, IelIdx, IdofGlob)
          return
        case (EL_Q2_3D)
          ! DOF`s in everything
          call storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
          call storage_getbase_int2D (p_rtriangulation%h_IedgesAtElement,p_2darray2)
          call storage_getbase_int2D (p_rtriangulation%h_IfacesAtElement,p_2darray3)
          call dof_locGlobUniMult_Q2_3D(p_rtriangulation%NVT,p_rtriangulation%NMT, &
              p_rtriangulation%NAT, p_2darray, p_2darray2, p_2darray3, IelIdx, IdofGlob)
          return
        case (EL_QP1_3D)
          ! DOF`s for QP1
          call dof_locGlobUniMult_QP1_3D(p_rtriangulation%NEL,IelIdx, IdofGlob)
          return
        case (EL_Q1T_3D)
          ! DOF`s in the face midpoints
          call storage_getbase_int2D(p_rtriangulation%h_IfacesAtElement,p_2darray)
          call dof_locGlobUniMult_Q1T_3D(p_2darray, IelIdx, IdofGlob)
          return
        case (EL_Q2T_3D)
          ! DOF`s in the face midpoints
          call storage_getbase_int2D(p_rtriangulation%h_IfacesAtElement,p_2darray)
          call dof_locGlobUniMult_Q2T_3D(p_rtriangulation%NAT,p_2darray, IelIdx, IdofGlob)
          return
        end select

      end if

    case default
      call output_line('Invalid discretisation!',&
                       OU_CLASS_ERROR,OU_MODE_STD,'dof_locGlobMapping_mult')
      call sys_halt()
    end select

    call output_line('Unsupported discretisation!',&
                     OU_CLASS_ERROR,OU_MODE_STD,'dof_locGlobMapping_mult')
    call sys_halt()

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_P0_1D(IelIdx, IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! All elements in the list are assumed to be P0_1D.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Global DOF = number of the element
      IdofGlob(1,i) = IelIdx(i)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_P1_1D(IverticesAtElement, IelIdx, IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! all elements in the list are assumed to be Q1.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! An array with the number of vertices adjacent to each element of the
  ! triangulation.
  integer, dimension(:,:), intent(in) :: IverticesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF`s - which are simply the vertex numbers of the
      ! corners.
      IdofGlob(1:2,i) = IverticesAtElement(1:2,IelIdx(i))
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_P2_1D(NVT, IverticesAtElement, IelIdx,&
                                           IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! all elements in the list are assumed to be Q1.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! Number of corner vertices in the triangulation
  integer, intent(in) :: NVT

  ! An array with the number of vertices adjacent to each element of the
  ! triangulation.
  integer, dimension(:,:), intent(in) :: IverticesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF`s - which are simply the vertex numbers of the
      ! corners and the cell midpoints.
      IdofGlob(1:2,i) = IverticesAtElement(1:2,IelIdx(i))
      IdofGlob(3,i) = NVT + IelIdx(i)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_S31_1D(NVT, IverticesAtElement, IelIdx,&
                                            IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! all elements in the list are assumed to be S31.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! Number of corner vertices in the triangulation
  integer, intent(in) :: NVT

  ! An array with the number of vertices adjacent to each element of the
  ! triangulation.
  integer, dimension(:,:), intent(in) :: IverticesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF`s - which are simply the vertex numbers of the
      ! corners.
      IdofGlob(1:2,i) = IverticesAtElement(1:2,IelIdx(i))
      IdofGlob(3:4,i) = NVT + IdofGlob(1:2,i)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_PN_1D(ndegree,NVT,NEL,IverticesAtElement,&
                                           IelIdx, IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! all elements in the list are assumed to be Q1.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! Specifies the local degree of the PN element
  integer, intent(in) :: ndegree

  ! The number of vertices in the triangulation
  integer, intent(in) :: NVT

  ! The number of elements in the triangulation
  integer, intent(in) :: NEL

  ! An array with the number of vertices adjacent to each element of the
  ! triangulation.
  integer, dimension(:,:), intent(in) :: IverticesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i,j

    ! Loop through the elements to handle
    do i = 1, size(IelIdx)
      IdofGlob(1:2,i) = IverticesAtElement(1:2,IelIdx(i))
      do j = 1, ndegree-1
        IdofGlob(2+j,i) = NVT + (j-1)*NEL + IelIdx(i)
      end do
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_P0Q0(IelIdx, IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! all elements in the list are assumed to be P0 or Q0.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Global DOF = number of the element
      IdofGlob(1,i) = IelIdx(i)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_P1Q1(IverticesAtElement, IelIdx, IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! all elements in the list are assumed to be either P1 or Q1.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! An array with the number of vertices adjacent to each element of the
  ! triangulation.
  integer, dimension(:,:), intent(in) :: IverticesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i,j

    ! Get the number of local DOF`s - usually either 3 or 4, depending on
    ! the element. The first dimension of IdofGlob indicates the number of
    ! DOF`s.
    j = min(ubound(IverticesAtElement,1),ubound(IdofGlob,1))

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF`s - which are simply the vertex numbers of the
      ! corners.
      IdofGlob(1:j,i) = IverticesAtElement(1:j,IelIdx(i))
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_QPW4P1(IverticesAtElement, IelIdx, IdofGlob, NVT)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! all elements in the list are assumed to be quadrilateral, piecewise linear
  ! P1 elements.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! An array with the number of vertices adjacent to each element of the
  ! triangulation.
  integer, dimension(:,:), intent(in) :: IverticesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

  ! Total number of vertices
  integer, intent(in) :: NVT

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF`s - which are simply the vertex numbers of the
      ! corners.
      IdofGlob(1:4,i) = IverticesAtElement(1:4,IelIdx(i))
      IdofGlob(5,i) = NVT + IelIdx(i)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_QPW4P2(IverticesAtElement, IedgesAtElement,&
      IelIdx, IdofGlob, NVT, NMT)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! all elements in the list are assumed to be quadrilateral, piecewise linear
  ! P1 elements.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! An array with the number of vertices adjacent to each element of the
  ! triangulation.
  integer, dimension(:,:), intent(in) :: IverticesAtElement

  ! An array with the number of edges adjacent to each element of the
  ! triangulation.
  integer, dimension(:,:), intent(in) :: IedgesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

  ! Total number of vertices
  integer, intent(in) :: NVT

  ! Total number of edges
  integer, intent(in) :: NMT

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF`s - which are simply the vertex numbers of the
      ! corners...
      IdofGlob(1:4,i) = IverticesAtElement(1:4,IelIdx(i))

      ! Plus the edge numbers...
      IdofGlob(5:8,i) = NVT + IedgesAtElement(1:4,IelIdx(i))

      ! plus 5 times the element, for the midpoint and the four inner edges
      IdofGlob(9,i)  = NVT + NMT + 5*(IelIdx(i)-1) + 1
      IdofGlob(10,i) = NVT + NMT + 5*(IelIdx(i)-1) + 2
      IdofGlob(11,i) = NVT + NMT + 5*(IelIdx(i)-1) + 3
      IdofGlob(12,i) = NVT + NMT + 5*(IelIdx(i)-1) + 4
      IdofGlob(13,i) = NVT + NMT + 5*(IelIdx(i)-1) + 5
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_DGP12D(IelIdx, IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! all elements in the list are assumed to be the DG P1.
  ! A uniform grid is assumed, i.e. a grid completely discretised by the
  ! same element.
!</description>

!<input>

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF`s. Every element gives 3 DOF`s which are
      ! 'internally' associated to the vertices but are not coupled.
      IdofGlob(1,i) = 1+3*(IelIdx(i)-1)
      IdofGlob(2,i) = 2+3*(IelIdx(i)-1)
      IdofGlob(3,i) = 3+3*(IelIdx(i)-1)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_DGQ12D(IelIdx, IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! all elements in the list are assumed to be the DG Q1.
  ! A uniform grid is assumed, i.e. a grid completely discretised by the
  ! same element.
!</description>

!<input>

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF`s. Every element gives 4 DOF`s which are
      ! 'internally' associated to the vertices but are not coupled.
      IdofGlob(1,i) = 1+4*(IelIdx(i)-1)
      IdofGlob(2,i) = 2+4*(IelIdx(i)-1)
      IdofGlob(3,i) = 3+4*(IelIdx(i)-1)
      IdofGlob(4,i) = 4+4*(IelIdx(i)-1)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_DGQ22D(IelIdx, IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! all elements in the list are assumed to be the DG Q2.
  ! A uniform grid is assumed, i.e. a grid completely discretised by the
  ! same element.
!</description>

!<input>

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF`s. Every element gives 9 DOF`s which are
      ! not coupled to any other element
      IdofGlob(1,i) = 1+9*(IelIdx(i)-1)
      IdofGlob(2,i) = 2+9*(IelIdx(i)-1)
      IdofGlob(3,i) = 3+9*(IelIdx(i)-1)
      IdofGlob(4,i) = 4+9*(IelIdx(i)-1)
      IdofGlob(5,i) = 5+9*(IelIdx(i)-1)
      IdofGlob(6,i) = 6+9*(IelIdx(i)-1)
      IdofGlob(7,i) = 7+9*(IelIdx(i)-1)
      IdofGlob(8,i) = 8+9*(IelIdx(i)-1)
      IdofGlob(9,i) = 9+9*(IelIdx(i)-1)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_P2(NVT, IverticesAtElement, &
                                        IedgesAtElement,IelIdx, IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! all elements in the list are assumed to be P2.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>
  integer, intent(in) :: NVT
  ! An array with the number of vertices adjacent to each element of the
  ! triangulation.
  integer, dimension(:,:), intent(in) :: IverticesAtElement

  ! An array with the number of edges adjacent to each element of the
  ! triangulation.
  integer, dimension(:,:), intent(in) :: IedgesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF`s.
      ! The P2 element has global DOF`s in the corners and edge midpoints
      ! of the triangles.
      !
      ! Take the numbers of the corners of the triangles at first.
      IdofGlob(1:3,i) = IverticesAtElement(1:3,IelIdx(i))

      ! Then append the numbers of the edges as midpoint numbers.
      ! Note that the number in this array is NVT+1..NVT+NMT.
      IdofGlob(4:6,i) = NVT + IedgesAtElement(1:3,IelIdx(i))

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_Q2(NVT,NMT,IverticesAtElement, &
                                        IedgesAtElement,IelIdx, IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! all elements in the list are assumed to be Q2.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>
  ! An array with the number of vertices adjacent to each element of the
  ! triangulation.
  integer, dimension(:,:), intent(in) :: IverticesAtElement

  ! An array with the number of edges adjacent to each element of the
  ! triangulation.
  integer, dimension(:,:), intent(in) :: IedgesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

  ! Number of corner vertices in the triangulation
  integer, intent(in) :: NVT

  ! Number of edes in the triangulation
  integer, intent(in) :: NMT
!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF`s.
      ! The Q2 element has global DOF`s in the corners, edge midpoints
      ! and element midpoints of the quads.
      !
      ! Take the numbers of the corners of the triangles at first.
      IdofGlob(1:4,i) = IverticesAtElement(1:4,IelIdx(i))

      ! Then append the numbers of the edges as midpoint numbers.
      ! Note that the number in this array is NVT+1..NVT+NMT.
      IdofGlob(5:8,i) = NVT+IedgesAtElement(1:4,IelIdx(i))

      ! At last append the element number - shifted by NVT+NMT to get
      ! a number behind.
      IdofGlob(9,i) = IelIdx(i)+NVT+NMT

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_P2Q2(NVT,NMT,NVE,IverticesAtElement, &
     IedgesAtElement,IelementCounter,IelIdx,IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! all elements in the list are assumed to be P2 or Q2.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>
  ! An array with the number of vertices adjacent to each element of the
  ! triangulation.
  integer, dimension(:,:), intent(in) :: IverticesAtElement

  ! An array with the number of edges adjacent to each element of the
  ! triangulation.
  integer, dimension(:,:), intent(in) :: IedgesAtElement

  ! Element counter array. This gives every triangle and every quad a
  ! unique running number (1,2,3,...)
  integer, dimension(:), intent(in) :: IelementCounter

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

  ! Number of corner vertices in the triangulation
  integer, intent(in) :: NVT

  ! Number of edes in the triangulation
  integer, intent(in) :: NMT

  ! Element type identifier for which type of elements is currently
  ! under view in IelIdx. All elements in IelIdx are assumed to be of
  ! the same type.
  ! =3: triangular, =4: quad.
  integer, intent(in) :: NVE
!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i

    if (NVE .eq. 3) then
      ! This element set consists of triangular elements.

      ! Loop through the elements to handle
      do i=1,size(IelIdx)
        ! Calculate the global DOF`s.
        ! The P2 element has global DOF`s in the corners and edge midpoints.
        !
        ! Take the numbers of the corners of the triangles at first.
        IdofGlob(1:3,i) = IverticesAtElement(1:3,IelIdx(i))

        ! Then append the numbers of the edges as midpoint numbers.
        ! Note that the number in this array is NVT+1..NVT+NMT.
        IdofGlob(4:6,i) = NVT+IedgesAtElement(1:3,IelIdx(i))

      end do

    else
      ! This element set consists of quad elements.

      ! Loop through the elements to handle
      do i=1,size(IelIdx)
        ! Calculate the global DOF`s.
        ! The Q2 element has global DOF`s in the corners, edge midpoints
        ! and element midpoints of the quads.
        !
        ! Take the numbers of the corners of the triangles at first.
        IdofGlob(1:4,i) = IverticesAtElement(1:4,IelIdx(i))

        ! Then append the numbers of the edges as midpoint numbers.
        ! Note that the number in this array is NVT+1..NVT+NMT.
        IdofGlob(5:8,i) = NVT+IedgesAtElement(1:4,IelIdx(i))

        ! At last append the element number - shifted by NVT+NMT to get
        ! a number behind. Note that we must not specify the actual element
        ! number here, but the element number in the set of quad elements!
        ! This is due to the fact that the element midpoints of triangular
        ! elements do not contribute to DOF`s in a mixed P2/Q2 discretisatrion!
        IdofGlob(9,i) = IelementCounter(IelIdx(i))+NVT+NMT

      end do

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_QP1(NEL, IelIdx, IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! all elements in the list are assumed to be QP1.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! Number of elements in the triangulation
  integer, intent(in) :: NEL

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! 1st Global DOF = number of the element = function value
      IdofGlob(1,i) = IelIdx(i)
      ! 2nd Global DOF = NEL + number of the element = X-derivative
      IdofGlob(2,i) = NEL+IelIdx(i)
      ! 3rd Global DOF = 2*NEL + number of the element = Y-derivative
      IdofGlob(3,i) = 2*NEL+IelIdx(i)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_P1T(IedgesAtElement, IelIdx, IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! All elements in the list are assumed to be E020.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! An array with the number of edges adjacent on each element of the
  ! triangulation.
  integer, dimension(:,:), intent(in) :: IedgesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF`s - which are simply the vertex numbers of the
      ! corners.
      ! We always copy all elements of IedgesAtElement (:,.).
      ! There is no harm and the compiler can optimise better.

      IdofGlob(1:TRIA_NVETRI2D,i) = IedgesAtElement(1:TRIA_NVETRI2D,IelIdx(i))
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_Q1T(IedgesAtElement, IelIdx, IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! All elements in the list are assumed to be E030, E031, EM30 or EM31.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! An array with the number of edges adjacent on each element of the
  ! triangulation.
  integer, dimension(:,:), intent(in) :: IedgesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF`s - which are simply the edge numbers.
      ! We always copy all elements of IedgesAtElement (:,.).
      ! There is no harm and the compiler can optimise better.

      IdofGlob(1:TRIA_NVEQUAD2D,i) = IedgesAtElement(1:TRIA_NVEQUAD2D,IelIdx(i))
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_Q1TB(iNMT, IedgesAtElement, IelIdx, IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! All elements in the list are assumed to be EB30.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! Number of edges in the triangulation.
  integer, intent(in) :: iNMT

  ! An array with the number of edges adjacent on each element of the
  ! triangulation.
  integer, dimension(:,:), intent(in) :: IedgesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF`s - which are simply the numbers of the
      ! edges. The DOF in the element gets the element number.
      ! We always copy all elements of IedgesAtElement (:,.).
      ! There is no harm and the compiler can optimise better.

      IdofGlob(1:TRIA_NVEQUAD2D,i) = IedgesAtElement(1:TRIA_NVEQUAD2D,IelIdx(i))
      IdofGlob(TRIA_NVEQUAD2D+1,i) = iNMT + IelIdx(i)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_Q2T(iNMT, IedgesAtElement, IelIdx, IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! All elements in the list are assumed to be E050.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! Number of edges in the triangulation.
  integer, intent(in) :: iNMT

  ! An array with the number of edges adjacent on each element of the
  ! triangulation.
  integer, dimension(:,:), intent(in) :: IedgesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF`s.
      ! The first 4 DOF`s are the number of the edges.
      ! The next 4 DOF`s are the number of the edges + nmt.
      ! The last DOF is the element number + 2*nmt.
      ! We always copy all elements of IedgesAtElement (:,.).
      ! There is no harm and the compiler can optimise better.

      IdofGlob(1:TRIA_NVEQUAD2D,i) = IedgesAtElement(1:TRIA_NVEQUAD2D,IelIdx(i))
      IdofGlob(TRIA_NVEQUAD2D+1:2*TRIA_NVEQUAD2D,i) = &
          IedgesAtElement(1:TRIA_NVEQUAD2D,IelIdx(i))+iNMT
      IdofGlob(2*TRIA_NVEQUAD2D+1,i) = 2*iNMT + IelIdx(i)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_Q2TB(iNMT, iNEL, &
      IedgesAtElement, IelIdx, IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! All elements in the list are assumed to be EB50.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! Number of edges in the triangulation.
  integer, intent(in) :: iNMT

  ! Number of elements in the triangulation.
  integer, intent(in) :: iNEL

  ! An array with the number of edges adjacent on each element of the
  ! triangulation.
  integer, dimension(:,:), intent(in) :: IedgesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF`s.
      ! The first 4 DOF`s are the number of the edges.
      ! The next 4 DOF`s are the number of the edges + nmt.
      ! The 9th DOF is at 2*nmt + element number
      ! The 10th DOF is at 2*nmt + inel + element number
      ! We always copy all elements of IedgesAtElement (:,.).
      ! There is no harm and the compiler can optimise better.

      IdofGlob(1:TRIA_NVEQUAD2D,i) = IedgesAtElement(1:TRIA_NVEQUAD2D,IelIdx(i))
      IdofGlob(TRIA_NVEQUAD2D+1:2*TRIA_NVEQUAD2D,i) = &
          IedgesAtElement(1:TRIA_NVEQUAD2D,IelIdx(i))+iNMT
      IdofGlob(2*TRIA_NVEQUAD2D+1,i) = 2*iNMT + IelIdx(i)
      IdofGlob(2*TRIA_NVEQUAD2D+2,i) = 2*iNMT + iNEL + IelIdx(i)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_Q3T_2D(iNMT, iNEL, IedgesAtElement, IelIdx, IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! All elements in the list are assumed to be E050.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! Number of edges in the triangulation.
  integer, intent(in) :: iNMT

  ! Number of elements in the triangulation.
  integer, intent(in) :: iNEL

  ! An array with the number of edges adjacent on each element of the
  ! triangulation.
  integer, dimension(:,:), intent(in) :: IedgesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF`s.
      ! The first 4 DOF`s are the number of the edges.
      ! The next 4 DOF`s are the number of the edges + nmt.
      ! The last DOF is the element number + 2*nmt.
      ! We always copy all elements of IedgesAtElement (:,.).
      ! There is no harm and the compiler can optimise better.

      IdofGlob(1:TRIA_NVEQUAD2D,i) = IedgesAtElement(1:TRIA_NVEQUAD2D,IelIdx(i))
      IdofGlob(TRIA_NVEQUAD2D+1:2*TRIA_NVEQUAD2D,i) = &
          IedgesAtElement(1:TRIA_NVEQUAD2D,IelIdx(i)) + iNMT
      IdofGlob(2*TRIA_NVEQUAD2D+1:3*TRIA_NVEQUAD2D,i) = &
          IedgesAtElement(1:TRIA_NVEQUAD2D,IelIdx(i)) + 2*iNMT
      IdofGlob(3*TRIA_NVEQUAD2D+1,i) = 3*iNMT + IelIdx(i)
      IdofGlob(3*TRIA_NVEQUAD2D+2,i) = 3*iNMT + iNEL + IelIdx(i)
      IdofGlob(3*TRIA_NVEQUAD2D+3,i) = 3*iNMT + 2*iNEL + IelIdx(i)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_P0Q0_3D(IelIdx, IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! all elements in the list are assumed to be P0 or Q0.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Global DOF = number of the element
      IdofGlob(1,i) = IelIdx(i)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_P1Q1_3D(IverticesAtElement, IelIdx,&
                                             IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! all elements in the list are assumed to be Q1.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! An array with the number of vertices adjacent to each element of the
  ! triangulation.
  integer, dimension(:,:), intent(in) :: IverticesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i,j

    ! Get the number of local DOF`s - usually either 3 or 4, depending on
    ! the element. The first dimension of IdofGlob indicates the number of
    ! DOF`s.
    j = min(ubound(IverticesAtElement,1),ubound(IdofGlob,1))

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF`s - which are simply the vertex numbers of the
      ! corners.
      IdofGlob(1:j,i) = IverticesAtElement(1:j,IelIdx(i))
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_Q2_3D(NVT,NMT,NAT,IverticesAtElement, &
                           IedgesAtElement,IfacesAtElement,IelIdx, IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! all elements in the list are assumed to be Q2.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>
  ! An array with the number of vertices adjacent to each element of the
  ! triangulation.
  integer, dimension(:,:), intent(in) :: IverticesAtElement

  ! An array with the number of edges adjacent to each element of the
  ! triangulation.
  integer, dimension(:,:), intent(in) :: IedgesAtElement

  ! An array with the number of faces adjacent to each element of the
  ! triangulation.
  integer, dimension(:,:), intent(in) :: IfacesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

  ! Number of corner vertices in the triangulation
  integer, intent(in) :: NVT

  ! Number of edes in the triangulation
  integer, intent(in) :: NMT

  ! Number of faces in the triangulation
  integer, intent(in) :: NAT
!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      IdofGlob(1:8,i) = IverticesAtElement(1:8,IelIdx(i))
      IdofGlob(9:20,i) = NVT + IedgesAtElement(1:12,IelIdx(i))
      IdofGlob(21:26,i) = NVT+NMT + IfacesAtElement(1:6,IelIdx(i))
      IdofGlob(27,i) = IelIdx(i)+NVT+NMT+NAT
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_QP1_3D(NEL, IelIdx, IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! all elements in the list are assumed to be QP1.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! Number of elements in the triangulation
  integer, intent(in) :: NEL

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! 1st Global DOF = number of the element = function value
      IdofGlob(1,i) = IelIdx(i)
      ! 2nd Global DOF = NEL + number of the element = X-derivative
      IdofGlob(2,i) = NEL+IelIdx(i)
      ! 3rd Global DOF = 2*NEL + number of the element = Y-derivative
      IdofGlob(3,i) = 2*NEL+IelIdx(i)
      ! 4th Global DOF = 3*NEL + number of the element = Z-derivative
      IdofGlob(4,i) = 3*NEL+IelIdx(i)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_Q1T_3D(IfacesAtElement, IelIdx, IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! all elements in the list are assumed to be Q1~.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! An array with the number of vertices adjacent to each element of the
  ! triangulation.
  integer, dimension(:,:), intent(in) :: IfacesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF`s - which are simply the face numbers.
      IdofGlob(1:6,i) = IfacesAtElement(1:6,IelIdx(i))
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_Q2T_3D(iNAT,IfacesAtElement, IelIdx, IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! all elements in the list are assumed to be Q2~.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! Number of faces in the triangulation.
  integer, intent(in) :: iNAT

  ! An array with the number of vertices adjacent to each element of the
  ! triangulation.
  integer, dimension(:,:), intent(in) :: IfacesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF`s - which are simply the face numbers.
      IdofGlob(1:6,i) = IfacesAtElement(1:6,IelIdx(i))
      IdofGlob( 7,i) = iNAT + 2*(IfacesAtElement(1,IelIdx(i))-1)+1
      IdofGlob( 8,i) = iNAT + 2*(IfacesAtElement(1,IelIdx(i))-1)+2
      IdofGlob( 9,i) = iNAT + 2*(IfacesAtElement(2,IelIdx(i))-1)+1
      IdofGlob(10,i) = iNAT + 2*(IfacesAtElement(2,IelIdx(i))-1)+2
      IdofGlob(11,i) = iNAT + 2*(IfacesAtElement(3,IelIdx(i))-1)+1
      IdofGlob(12,i) = iNAT + 2*(IfacesAtElement(3,IelIdx(i))-1)+2
      IdofGlob(13,i) = iNAT + 2*(IfacesAtElement(4,IelIdx(i))-1)+1
      IdofGlob(14,i) = iNAT + 2*(IfacesAtElement(4,IelIdx(i))-1)+2
      IdofGlob(15,i) = iNAT + 2*(IfacesAtElement(5,IelIdx(i))-1)+1
      IdofGlob(16,i) = iNAT + 2*(IfacesAtElement(5,IelIdx(i))-1)+2
      IdofGlob(17,i) = iNAT + 2*(IfacesAtElement(6,IelIdx(i))-1)+1
      IdofGlob(18,i) = iNAT + 2*(IfacesAtElement(6,IelIdx(i))-1)+2
      IdofGlob(19,i) = 3*iNAT + IelIdx(i)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_DG_T0_1D(IelIdx, IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! All elements in the list are assumed to be DG_T0_1D.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Global DOF = number of the element
      IdofGlob(1,i) = IelIdx(i)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_DG_T1_1D(IelIdx, IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! All elements in the list are assumed to be DG_T1_1D.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Global DOF = number of the element
      IdofGlob(1,i) = 2*IelIdx(i)-1
      IdofGlob(2,i) = 2*IelIdx(i)
    end do

  end subroutine


  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_DG_T2_1D(IelIdx, IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! All elements in the list are assumed to be DG_T2_1D.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Global DOF = number of the element
      IdofGlob(1,i) = 3*IelIdx(i)-2
      IdofGlob(2,i) = 3*IelIdx(i)-1
      IdofGlob(3,i) = 3*IelIdx(i)
    end do

  end subroutine

    ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_DG_T0_2D(IelIdx, IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! All elements in the list are assumed to be DG_T0_1D.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Global DOF = number of the element
      IdofGlob(1,i) = IelIdx(i)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_DG_T1_2D(IelIdx, IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! All elements in the list are assumed to be DG_T1_1D.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Global DOF = number of the element
      IdofGlob(1,i) = 3*IelIdx(i)-2
      IdofGlob(2,i) = 3*IelIdx(i)-1
      IdofGlob(3,i) = 3*IelIdx(i)
    end do

  end subroutine


  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_DG_T2_2D(IelIdx, IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! All elements in the list are assumed to be DG_T2_1D.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Global DOF = number of the element
      IdofGlob(1,i) = 6*IelIdx(i)-5
      IdofGlob(2,i) = 6*IelIdx(i)-4
      IdofGlob(3,i) = 6*IelIdx(i)-3
      IdofGlob(4,i) = 6*IelIdx(i)-2
      IdofGlob(5,i) = 6*IelIdx(i)-1
      IdofGlob(6,i) = 6*IelIdx(i)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine dof_locGlobUniMult_QPW4DCP1_2D(IelIdx, IdofGlob)

!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! All elements in the list are assumed to be DG_T2_1D.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(in) :: IelIdx

!</input>

!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF`s.
  integer, dimension(:,:), intent(out) :: IdofGlob

!</output>

!</subroutine>

  ! local variables
  integer :: i,j

    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      do j = 1, 11
        IdofGlob(j,i) = 11*(IelIdx(i)-1) + j
      end do
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine dof_infoDiscr (rspatialDiscr)

!<description>
  ! This routine prints out statistical information about a discretisation
  ! to the terminal.
!</description>

!<inputoutput>
  ! The discretisation structure where information should be printed.
  type(t_spatialDiscretisation), intent(in), target :: rspatialDiscr
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i
    type(t_elementDistribution), pointer :: p_relementDistr

    ! General information:
    call output_line ('Dimension:                    '&
        //trim(sys_siL(rspatialDiscr%ndimension,10)))
    call output_line ('Complexity:                   ',bnolinebreak=.true.,&
        bnoTrim=.true.)
    select case (rspatialDiscr%ccomplexity)
    case (SPDISC_UNIFORM)
      call output_line ('uniform',cdateTimeLogPolicy=OU_DTP_NONE)
    case (SPDISC_CONFORMAL)
      call output_line ('conformal',cdateTimeLogPolicy=OU_DTP_NONE)
    case (SPDISC_MIXED)
      call output_line ('mixed',cdateTimeLogPolicy=OU_DTP_NONE)
    case DEFAULT
      call output_line ('undefined',cdateTimeLogPolicy=OU_DTP_NONE)
    end select
    call output_line ('#DOFs:                        '&
        //trim(sys_siL(dof_igetNDofGlob(rspatialDiscr),16)))
    call output_line ('#finite element spaces:       '&
        //trim(sys_siL(rspatialDiscr%inumFESpaces,10)))

    ! Print out detailed information about the FE spaces.
    call output_line ('Discretisation details:')
    call output_line ('FE-space #elements       NVE   trial-element   test-element')

    ! Loop through all element distributions
    do i=1,rspatialDiscr%inumFESpaces

      p_relementDistr => rspatialDiscr%RelementDistr(i)

      call output_line ( ' ' &
        // sys_siL(i,8) &
        // sys_siL(p_relementDistr%NEL,16) &
        // sys_siL(elem_igetNVE(p_relementDistr%celement),6) &
        // sys_siL(int(iand(elem_getPrimaryElement(p_relementDistr%celement),&
                   not(EL_DIMENSION))),16) )

    end do

  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine dof_infoDiscrBlock (rblockDiscr,bdetailed)

!<description>
  ! This routine prints out statistical information about a block
  ! discretisation to the terminal.
!</description>

!<input>
  ! Whether a detailed description is printed to the terminal or not.
  ! FALSE prints out information only about the block discretisation.
  ! TRUE prints out more detailed information about the block
  ! discretisation, the structure of the blocks, the used FE spaces etc.
  logical, intent(in) :: bdetailed

  ! The discretisation structure where information should be printed.
  type(t_blockDiscretisation), intent(in), target :: rblockDiscr
!</input>

!</subroutine>

    integer :: i

    call output_line ('Dimension:                    '&
        //trim(sys_siL(rblockDiscr%ndimension,10)))
    call output_line ('Complexity:                   ',bnolinebreak=.true.,&
        bnoTrim=.true.)
    select case (rblockDiscr%ccomplexity)
    case (SPDISC_UNIFORM)
      call output_line ('uniform',cdateTimeLogPolicy=OU_DTP_NONE)
    case (SPDISC_CONFORMAL)
      call output_line ('conformal',cdateTimeLogPolicy=OU_DTP_NONE)
    case (SPDISC_MIXED)
      call output_line ('mixed',cdateTimeLogPolicy=OU_DTP_NONE)
    case DEFAULT
      call output_line ('undefined',cdateTimeLogPolicy=OU_DTP_NONE)
    end select
    call output_line ('#DOFs:                        '&
        //trim(sys_siL(dof_igetNDofGlobBlock(rblockDiscr),16)))

    call output_line ('Number of components:         '&
        //trim(sys_siL(rblockDiscr%ncomponents,10)))

    if (bdetailed) then
      do i=1,rblockDiscr%ncomponents
        call output_lbrk ()
        call output_line ('Solution component:           '//trim(sys_siL(i,10)))
        call dof_infoDiscr(rblockDiscr%RspatialDiscr(i))
      end do
    end if

  end subroutine

! ***************************************************************************
!<subroutine>

  subroutine dof_precomputeDofMapping(rdiscretisation)

!<description>
  ! Precomputes the DOF-mapping. Simple discretisations may work with
  ! a direct mapping. More complex discretisations (hanging nodes e.g.) must
  ! precompute the DOF-mapping to determine how to map local DOF`s to
  ! global DOF`s.
!</description>

!<inputoutput>
  ! The discretisation structure where to precompute the DOF-mapping.
  type(t_spatialDiscretisation), intent(inout) :: rdiscretisation
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: ndoflocal,i,j,ieldistr,ipos
    integer, dimension(:,:), allocatable :: IelementDofs
    integer, dimension(:), pointer :: p_IelementDofs,p_IelementDofIdx,p_IelementList
    type(t_triangulation), pointer :: p_rtriangulation

    ! Release old arrays if necessary.
    call spdiscr_releaseDofMapping(rdiscretisation)

    ! Is that a simple Dof-mapping that can directly be precomputed?
    rdiscretisation%ndof = dof_igetNDofGlob(rdiscretisation)
    if (rdiscretisation%ndof .ne. 0) then

      ! Allocate lists for precomputing the DOF`s.
      p_rtriangulation => rdiscretisation%p_rtriangulation

      ! Allocate a temp array accepting the local DOF`s.
      !
      ! At first allocate an index array; this has to be build in advance.
      call storage_new ('dof_precomputeDofMapping', 'h_IelementDofIdx', &
          p_rtriangulation%NEL+1, ST_INT, rdiscretisation%h_IelementDofIdx,   &
          ST_NEWBLOCK_ZERO)
      call storage_getbase_int(rdiscretisation%h_IelementDofIdx,p_IelementDofIdx)

      ! Build that array by looping through the element distributions and
      ! counting the number of local DOF`s.
      do ieldistr = 1,rdiscretisation%inumFESpaces
        call storage_getbase_int(rdiscretisation%RelementDistr(ieldistr)%h_IelementList,&
            p_IelementList)
        ndoflocal = elem_igetNDofLoc(rdiscretisation%RelementDistr(ieldistr)%celement)
        do i = 1,rdiscretisation%RelementDistr(ieldistr)%NEL
          p_IelementDofIdx(1+p_IelementList(i)) = ndoflocal
        end do
      end do

      ! Sum the values up to get the actual index array.
      p_IelementDofIdx(1) = 1
      do i = 2,p_rtriangulation%NEL+1
        p_IelementDofIdx(i) = p_IelementDofIdx(i) + p_IelementDofIdx(i-1)
      end do

      ! Now get the actual DOF`s.
      call storage_new ('dof_precomputeDofMapping', 'h_IelementDofs', &
          p_IelementDofIdx(p_rtriangulation%NEL+1)-1, ST_INT, rdiscretisation%h_IelementDofs,&
          ST_NEWBLOCK_ZERO)
      call storage_getbase_int(rdiscretisation%h_IelementDofs,p_IelementDofs)

      do ieldistr = 1,rdiscretisation%inumFESpaces

        call storage_getbase_int(rdiscretisation%RelementDistr(ieldistr)%h_IelementList,&
            p_IelementList)
        ndoflocal = elem_igetNDofLoc(rdiscretisation%RelementDistr(ieldistr)%celement)

        ! Allocate temp memory for the DOF`s. Write them to this array at first
        ! and copy them later to the actual array.
        allocate(IelementDofs(ndoflocal,rdiscretisation%RelementDistr(ieldistr)%NEL))

        call dof_locGlobMapping_mult(rdiscretisation, p_IelementList, IelementDofs)

        do i = 1,rdiscretisation%RelementDistr(ieldistr)%NEL
          ipos = p_IelementDofIdx(p_IelementList(i))
          do j=1,ndoflocal
            p_IelementDofs(ipos+j-1) = IelementDofs(j,i)
          end do
        end do

        ! Release temp memory, that is it.
        deallocate (IelementDofs)

      end do

      ! Remember that the DOF-mapping is precomputed now.
      rdiscretisation%bprecompiledDofMapping = .true.

    else
      call output_line ('Discretisation does not support precomputed DOF''s!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'dof_precomputeDofMapping')
      call sys_halt()
    end if

  end subroutine

end module
