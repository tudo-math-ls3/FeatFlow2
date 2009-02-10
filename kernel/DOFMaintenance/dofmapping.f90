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
!# The following functions provide support for global DOF's:
!#
!# 1.) dof_igetNDofGlob
!#     -> Get number of global DOF's described by a given discretisation
!#
!# 2.) dof_igetNDofGlobBlock
!#     -> Get number of global DOF's described by a given block discretisation
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
!# </purpose>
!##############################################################################

module dofmapping

  use fsystem
  use spatialdiscretisation
  use triangulation
  use element
  
  implicit none

!<constants>

  !<constantblock description="Kind values for global DOF's">

  ! kind value for indexing global DOF's
  ! !!! DEPRECATED: DO NOT USE THIS CONSTANT ANYMORE !!!
  integer, parameter :: PREC_DOFIDX     = I32

  !</constantblock>

!</constants>

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
  type(t_spatialDiscretisation), intent(IN) :: rdiscretisation
!</input>

!<result>
  ! global number of equations on current grid
!</result>

!</function>

  integer(I32) :: ieltyp
  integer(I32), dimension(2) :: IelTypes

  dof_igetNDofGlob = 0

  select case(rdiscretisation%ndimension)
  case (NDIM1D)
  
      ieltyp = rdiscretisation%RelementDistr(1)%celement
      dof_igetNDofGlob = NDFG_uniform1D (rdiscretisation%p_rtriangulation, ieltyp)

  case (NDIM2D)
    if (rdiscretisation%ccomplexity .eq. SPDISC_UNIFORM) then

      ieltyp = rdiscretisation%RelementDistr(1)%celement
      ! Uniform discretisation - fall back to the old FEAT mapping
      dof_igetNDofGlob = NDFG_uniform2D (rdiscretisation%p_rtriangulation, ieltyp)

    else if (rdiscretisation%ccomplexity .eq. SPDISC_CONFORMAL) then

      ! Conformal discretisation. That's a little bit tricky!
      ! At first, we support only the case where two element types are mixed.
      if (rdiscretisation%inumFESpaces .eq. 2) then
        IelTypes(1) = rdiscretisation%RelementDistr(1)%celement
        IelTypes(2) = rdiscretisation%RelementDistr(2)%celement

        dof_igetNDofGlob = NDFG_conformal2D_2el (&
            rdiscretisation%p_rtriangulation, IelTypes(1:2))
        
      end if

    end if
  
  case (NDIM3D)
    ! Currently, only uniform discretisations are supported.
    if (rdiscretisation%ccomplexity .eq. SPDISC_UNIFORM) then

      ieltyp = rdiscretisation%RelementDistr(1)%celement

      ! Uniform discretisation - fall back to the old FEAT mapping
      dof_igetNDofGlob = NDFG_uniform3D (rdiscretisation%p_rtriangulation, ieltyp)
    
    end if

  case default
  
    ! Dimension not supported
    print *,'dof_igetNDofGlob: Invalid discretisation dimension!'
    call sys_halt()
  
  end select
  
  contains
  
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Internal subroutine: Get global DOF number for uniform discretisation
    ! with only one element type.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    integer function NDFG_uniform1D (rtriangulation, ieltype)
    
    ! IN: The underlying triangulation
    type(t_triangulation), intent(IN) :: rtriangulation
    
    ! IN: The element type of the discretisation
    integer(I32), intent(IN) :: ieltype
    
    ! OUT: number of global DOF's.
    
    ! The number of global DOF's depends on the element type...
    select case (elem_getPrimaryElement(ieltype))
    case (EL_P0_1D)
      ! DOF's in the cell midpoints
      NDFG_uniform1D = rtriangulation%NEL
    case (EL_P1_1D)
      ! DOF's in the vertices
      NDFG_uniform1D = rtriangulation%NVT
    case (EL_P2_1D)
      ! DOF's in the vertices and cell midpoints
      NDFG_uniform1D = rtriangulation%NVT + rtriangulation%NEL
    case (EL_S31_1D)
      ! DOF's in the vertices
      NDFG_uniform1D = 2*rtriangulation%NVT
    end select
    
    end function

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Internal subroutine: Get global DOF number for uniform discretisation
    ! with only one element type.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! This is roughly the NDFG routine of the old FEAT library...
    
    integer function NDFG_uniform2D (rtriangulation, ieltype)
    
    ! IN: The underlying triangulation
    type(t_triangulation), intent(IN) :: rtriangulation
    
    ! IN: The element type of the discretisation
    integer(I32), intent(IN) :: ieltype
    
    ! OUT: number of global DOF's.
    
    ! The number of global DOF's depends on the element type...
    select case (elem_getPrimaryElement(ieltype))
    case (EL_P0, EL_Q0)
      ! DOF's in the cell midpoints
      NDFG_uniform2D = rtriangulation%NEL
    case (EL_P1, EL_Q1)
      ! DOF's in the vertices
      NDFG_uniform2D = rtriangulation%NVT
    case (EL_P2)
      ! DOF's in the vertices and edge midpoints
      NDFG_uniform2D = rtriangulation%NVT + rtriangulation%NMT
    case (EL_Q2)
      ! DOF's in the vertices, edge midpoints and element midpoints
      NDFG_uniform2D = rtriangulation%NVT + rtriangulation%NMT + rtriangulation%NEL
    case (EL_P3)
      ! 1 DOF's per vertices, 2 DOF per edge 
      NDFG_uniform2D = rtriangulation%NVT + 2*rtriangulation%NMT
    case (EL_Q3)
      ! 1 DOF's per vertices, 2 DOF per edge, 4 DOF in the inner
      NDFG_uniform2D = rtriangulation%NVT + 2*rtriangulation%NMT + 4*rtriangulation%NEL
    case (EL_QP1)
      ! 3 DOF's in the midpoint of the element.
      NDFG_uniform2D = 3*rtriangulation%NEL
    case (EL_P1T, EL_Q1T)
      ! 1 DOF per edge
      NDFG_uniform2D = rtriangulation%NMT
    case (EL_Q1TB)
      ! 1 DOF per edge, one per element
      NDFG_uniform2D = rtriangulation%NMT + rtriangulation%NEL
    case (EL_Q2T)
      ! DOF's in the vertices, edge midpoints and element midpoints
      NDFG_uniform2D = 2*rtriangulation%NMT + rtriangulation%NEL
    case (EL_Q2TB)
      ! E037
      NDFG_uniform2D = 2*rtriangulation%NMT + 2*rtriangulation%NEL
    end select
    
    end function

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Internal subroutine: Get global DOF number for conformal discretisation
    ! with two element types.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    integer function NDFG_conformal2D_2el (rtriangulation, IelTypes)
    
    ! IN: The underlying triangulation
    type(t_triangulation), intent(IN) :: rtriangulation
    
    ! IN: List of element types in the discretisation. IelTypes(1) is one element
    ! type identifier, IelTypes(2) the other one.
    integer(I32), dimension(:), intent(IN) :: IelTypes
    
    ! OUT: number of global DOF's.
    
    ! local variables
    integer(I32), dimension(size(IelTypes)) :: IelTypesPrimary
    
    ! Get the primary element number
    IelTypesPrimary = elem_getPrimaryElement(IelTypes)
    
    ! The number of global DOF's depends on the element type...
    select case (IelTypesPrimary(1))
    case (EL_P0, EL_Q0)
      select case (IelTypesPrimary(2))
      case (EL_P0, EL_Q0)
        ! DOF's in the cell midpoints
        NDFG_conformal2D_2el = rtriangulation%NEL
      end select
      
    case (EL_P1, EL_Q1)
      select case (IelTypesPrimary(2))
      case (EL_P1, EL_Q1)
        ! DOF's in the vertices
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
        ! DOF's in the edge midpoints
        NDFG_conformal2D_2el = rtriangulation%NMT
      
      end select
      
    end select

    end function

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Internal subroutine: Get global DOF number for uniform discretisation
    ! with only one element type.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! This is roughly the NDFG routine of the old FEAT library...
    
    integer function NDFG_uniform3D (rtriangulation, ieltype)
    
    ! IN: The underlying triangulation
    type(t_triangulation), intent(IN) :: rtriangulation
    
    ! IN: The element type of the discretisation
    integer(I32), intent(IN) :: ieltype
    
    ! OUT: number of global DOF's.
    
    ! The number of global DOF's depends on the element type...
    select case (elem_getPrimaryElement(ieltype))
    case (EL_P0_3D, EL_Q0_3D, EL_Y0_3D, EL_R0_3D)
      ! DOF's in the cell midpoints
      NDFG_uniform3D = rtriangulation%NEL
    case (EL_P1_3D, EL_Q1_3D, EL_Y1_3D, EL_R1_3D)
      ! DOF's in the vertices
      NDFG_uniform3D = rtriangulation%NVT
    case (EL_QP1_3D)
      ! 4 DOF's in the midpoint of the element.
      NDFG_uniform3D = 4*rtriangulation%NEL
    case (EL_Q1T_3D)
      ! DOF's in the face midpoints
      NDFG_uniform3D = rtriangulation%NAT
    case (EL_Q2T_3D)
      ! DOF's in the face midpoints
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
  type(t_blockDiscretisation), intent(IN) :: rdiscretisation
!</input>

!<result>
  ! global number of equations on current grid
!</result>

!</function>

    ! Sum up the DOF's of every scalar sub-discretisation structure
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
  ! On each element, there are a number of local DOF's 1..n (where n can
  ! be obtained using elem_igetNDofLoc). This subroutine returns the
  ! corresponding global degrees of freedom, corresponding to these local
  ! ones.
!</description>

!<input>

  ! The discretisation structure that specifies the (scalar) discretisation.
  type(t_spatialDiscretisation), intent(IN) :: rdiscretisation

  ! Element index, where the mapping should be computed.
  integer, intent(IN) :: ielIdx

!</input>
    
!<output>

  ! array of global DOF numbers
  integer, dimension(:), intent(OUT) :: IdofGlob

!</output>

! </subroutine>

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
  ! On each element, there are a number of local DOF's 1..n (where n can
  ! be obtained using elem_igetNDofLoc). This subroutine returns the
  ! corresponding global degrees of freedom, corresponding to these local
  ! ones.
  ! The routine allows to calculate the DOF's to a list of elements
  ! which are all of the same type. IelIdx is this list of elements and
  ! IdofGlob is a 2D array which receives for every element the
  ! corresponding global DOF's.
!</description>

!<input>

  ! The discretisation structure that specifies the (scalar) discretisation.
  type(t_spatialDiscretisation), intent(IN) :: rdiscretisation

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(IN) :: IelIdx
  
!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  integer, dimension(:,:), intent(OUT) :: IdofGlob

!</output>

! </subroutine>

    ! local variables
    integer, dimension(:,:), pointer :: p_2darray,p_2darray2
    integer, dimension(:), pointer :: p_IelementCounter
    type(t_triangulation), pointer :: p_rtriangulation     
    integer(I32) :: ieltype
    integer(I32), dimension(2) :: IelTypes

    p_rtriangulation => rdiscretisation%p_rtriangulation
    
    select case(rdiscretisation%ndimension)
    case (NDIM1D)
      ! Call the right 'multiple-get' routines for global DOF's.
      ! For this purpose we evaluate the pointers in the discretisation
      ! structure (if necessary) to prevent another call using pointers...
      ! The number of global DOF's depends on the element type...
      
      ieltype = rdiscretisation%RelementDistr(1)%celement
      
      select case (elem_getPrimaryElement(ieltype))
      case (EL_P0_1D)
        ! DOF's for P0
        call dof_locGlobUniMult_P0_1D(IelIdx, IdofGlob)
        return
      case (EL_P1_1D)
        ! DOF's in the vertices
        call storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
        call dof_locGlobUniMult_P1_1D(p_2darray, IelIdx, IdofGlob)
        return
      case (EL_P2_1D)
        ! DOF's in the vertices and line midpoints
        call storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
        call dof_locGlobUniMult_P2_1D(p_rtriangulation%NVT, p_2darray, IelIdx,&
                                      IdofGlob)
        return
      case (EL_S31_1D)
        ! DOF's in the vertices
        call storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
        call dof_locGlobUniMult_S31_1D(p_rtriangulation%NVT, p_2darray, IelIdx,&
                                       IdofGlob)
        return
      end select
        
    case (NDIM2D)
      ! At first we deal only with uniform discretisations
      if (rdiscretisation%ccomplexity .eq. SPDISC_UNIFORM) then
      
        ! Call the right 'multiple-get' routines for global DOF's.
        ! For this purpose we evaluate the pointers in the discretisation
        ! structure (if necessary) to prevent another call using pointers...
        ! The number of global DOF's depends on the element type...
        
        ieltype = rdiscretisation%RelementDistr(1)%celement
        
        select case (elem_getPrimaryElement(ieltype))
        case (EL_P0, EL_Q0)
          ! DOF's for Q0
          call dof_locGlobUniMult_P0Q0(IelIdx, IdofGlob)
          return
        case (EL_P1, EL_Q1)
          ! DOF's in the vertices
          call storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
          call dof_locGlobUniMult_P1Q1(p_2darray, IelIdx, IdofGlob)
          return
        case (EL_P2)
          ! DOF's in the vertices and egde midpoints
          call storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
          call storage_getbase_int2D (p_rtriangulation%h_IedgesAtElement,p_2darray2)
          call dof_locGlobUniMult_P2(p_rtriangulation%NVT,p_2darray, p_2darray2,&
                                     IelIdx, IdofGlob)
          return
        case (EL_Q2)
          ! DOF's in the vertices, egde midpoints and element midpoints
          call storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
          call storage_getbase_int2D (p_rtriangulation%h_IedgesAtElement,p_2darray2)
          call dof_locGlobUniMult_Q2(p_rtriangulation%NVT,p_rtriangulation%NMT,&
                                    p_2darray, p_2darray2, IelIdx, IdofGlob)
          return
        case (EL_P3) 
        case (EL_Q3) 

        case (EL_QP1)
          ! DOF's for Q1
          call dof_locGlobUniMult_QP1(p_rtriangulation%NEL,IelIdx, IdofGlob)
          return
        case (EL_P1T)
          ! DOF's in the edges
          call storage_getbase_int2D (p_rtriangulation%h_IedgesAtElement,p_2darray)
          call dof_locGlobUniMult_P1T(p_2darray, IelIdx, IdofGlob)
          return
        case (EL_Q1T)
          ! DOF's in the edges
          call storage_getbase_int2D (p_rtriangulation%h_IedgesAtElement,p_2darray)
          call dof_locGlobUniMult_Q1T(p_2darray, IelIdx, IdofGlob)
          return
        case (EL_Q1TB)
          ! DOF's in the edges
          call storage_getbase_int2D (p_rtriangulation%h_IedgesAtElement,p_2darray)
          call dof_locGlobUniMult_Q1TB(p_rtriangulation%NMT,p_2darray, IelIdx, IdofGlob)
          return
        case (EL_Q2T)
          ! DOF's in the edges and the element center
          call storage_getbase_int2D (p_rtriangulation%h_IedgesAtElement,p_2darray)
          call dof_locGlobUniMult_Q2T(p_rtriangulation%NMT,p_2darray, IelIdx, IdofGlob)
          return
        case (EL_Q2TB)
          ! DOF's in the edges and the element center
          call storage_getbase_int2D (p_rtriangulation%h_IedgesAtElement,p_2darray)
          call dof_locGlobUniMult_Q2TB(p_rtriangulation%NMT,p_rtriangulation%NEL,&
                                       p_2darray, IelIdx, IdofGlob)
          return
        end select
        
        return
      
      else if (rdiscretisation%ccomplexity .eq. SPDISC_CONFORMAL) then
      
        ! Conformal discretisation. That's a little bit tricky!
        ! At first, we support only the case where two element types are mixed.
        if (rdiscretisation%inumFESpaces .eq. 2) then

          IelTypes(1) = rdiscretisation%RelementDistr(1)%celement
          IelTypes(2) = rdiscretisation%RelementDistr(2)%celement

          ! Get the primary element type(s)
          IelTypes = elem_getPrimaryElement(IelTypes)

          ! Now the actual mappings...
          select case (IelTypes(1))
          case (EL_P0, EL_Q0)
            select case (IelTypes(2))
            case (EL_P0, EL_Q0)
              ! DOF's in the cell midpoints.
              ! That works like P0 elements.
              call dof_locGlobUniMult_P0Q0(IelIdx, IdofGlob)
              return
            end select
            
          case (EL_P1, EL_Q1)
            select case (IelTypes(2))
            case (EL_P1, EL_Q1)
              ! DOF's in the vertices.
              ! That works like P1 elements.
              call storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
              call dof_locGlobUniMult_P1Q1(p_2darray, IelIdx, IdofGlob)
              return
            end select
            
          case (EL_P2, EL_Q2)
            select case (IelTypes(2))
            case (EL_P2, EL_Q2)
              ! DOF's in the vertices, edges and element mitpoints of the quads.
              ! For this purpose, we need the element counter array that counts
              ! every quad element.
              call storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
              call storage_getbase_int2D (p_rtriangulation%h_IedgesAtElement,p_2darray2)
              call storage_getbase_int (rdiscretisation%h_IelementCounter,p_IelementCounter)
              
              ! Use p_IverticesAtElement evaluated at the first element in the element
              ! set to determine NVE. It's either 3 or 4 and valid for all elements
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
            select case (IelTypes(2))
            case (EL_P1T, EL_Q1T)
              ! DOF's in the edges
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
      
        ! Call the right 'multiple-get' routines for global DOF's.
        ! For this purpose we evaluate the pointers in the discretisation
        ! structure (if necessary) to prevent another call using pointers...
        ! The number of global DOF's depends on the element type...
        ieltype = rdiscretisation%RelementDistr(1)%celement
        
        select case (elem_getPrimaryElement(ieltype))
        case (EL_P0_3D, EL_Q0_3D, EL_Y0_3D, EL_R0_3D)
          ! DOF's for Q0
          call dof_locGlobUniMult_P0Q0_3D(IelIdx, IdofGlob)
          return
        case (EL_P1_3D, EL_Q1_3D, EL_Y1_3D, EL_R1_3D)
          ! DOF's in the vertices
          call storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
          call dof_locGlobUniMult_P1Q1_3D(p_2darray, IelIdx, IdofGlob)
          return
        case (EL_QP1)
          ! DOF's for QP1
          call dof_locGlobUniMult_QP1_3D(p_rtriangulation%NEL,IelIdx, IdofGlob)
          return
        case (EL_Q1T_3D)
          ! DOF's in the face midpoints
          call storage_getbase_int2D(p_rtriangulation%h_IfacesAtElement,p_2darray)
          call dof_locGlobUniMult_Q1T_3D(p_2darray, IelIdx, IdofGlob)
          return
        case (EL_Q2T_3D)
          ! DOF's in the face midpoints
          call storage_getbase_int2D(p_rtriangulation%h_IfacesAtElement,p_2darray)
          call dof_locGlobUniMult_Q2T_3D(p_rtriangulation%NAT,p_2darray, IelIdx, IdofGlob)
          return
        end select
        
      end if
    
    case default
      print *,'dof_locGlobMapping_mult: invalid discretisation!'
      call sys_halt()
    end select

    print *,'dof_locGlobMapping_mult: Unsupported discretisation!'
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
  integer, dimension(:), intent(IN) :: IelIdx

!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  integer, dimension(:,:), intent(OUT) :: IdofGlob

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
  integer, dimension(:,:), intent(IN) :: IverticesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(IN) :: IelIdx
  
!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  integer, dimension(:,:), intent(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  integer :: i
  
    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF's - which are simply the vertex numbers of the 
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
  integer, intent(IN) :: NVT

  ! An array with the number of vertices adjacent to each element of the
  ! triangulation.
  integer, dimension(:,:), intent(IN) :: IverticesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(IN) :: IelIdx
  
!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  integer, dimension(:,:), intent(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  integer :: i
  
    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF's - which are simply the vertex numbers of the 
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
  integer, intent(IN) :: NVT

  ! An array with the number of vertices adjacent to each element of the
  ! triangulation.
  integer, dimension(:,:), intent(IN) :: IverticesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(IN) :: IelIdx
  
!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  integer, dimension(:,:), intent(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  integer :: i
  
    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF's - which are simply the vertex numbers of the 
      ! corners.
      IdofGlob(1:2,i) = IverticesAtElement(1:2,IelIdx(i))
      IdofGlob(3:4,i) = NVT + IdofGlob(1:2,i)
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
  integer, dimension(:), intent(IN) :: IelIdx

!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  integer, dimension(:,:), intent(OUT) :: IdofGlob

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
  integer, dimension(:,:), intent(IN) :: IverticesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(IN) :: IelIdx
  
!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  integer, dimension(:,:), intent(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  integer :: i,j
  
    ! Get the number of local DOF's - usually either 3 or 4, depending on
    ! the element. The first dimension of IdofGlob indicates the number of 
    ! DOF's.
    j = min(ubound(IverticesAtElement,1),ubound(IdofGlob,1))
    
    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF's - which are simply the vertex numbers of the 
      ! corners.
      IdofGlob(1:j,i) = IverticesAtElement(1:j,IelIdx(i))
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
  integer, intent(IN) :: NVT
  ! An array with the number of vertices adjacent to each element of the
  ! triangulation.
  integer, dimension(:,:), intent(IN) :: IverticesAtElement

  ! An array with the number of edges adjacent to each element of the
  ! triangulation.
  integer, dimension(:,:), intent(IN) :: IedgesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(IN) :: IelIdx
  
!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  integer, dimension(:,:), intent(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  integer :: i
  
    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF's.
      ! The P2 element has global DOF's in the corners and edge midpoints
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
  integer, dimension(:,:), intent(IN) :: IverticesAtElement

  ! An array with the number of edges adjacent to each element of the
  ! triangulation.
  integer, dimension(:,:), intent(IN) :: IedgesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(IN) :: IelIdx
  
  ! Number of corner vertices in the triangulation
  integer, intent(IN) :: NVT
  
  ! Number of edes in the triangulation
  integer, intent(IN) :: NMT
!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  integer, dimension(:,:), intent(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  integer :: i
  
    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF's.
      ! The Q2 element has global DOF's in the corners, edge midpoints
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
  integer, dimension(:,:), intent(IN) :: IverticesAtElement

  ! An array with the number of edges adjacent to each element of the
  ! triangulation.
  integer, dimension(:,:), intent(IN) :: IedgesAtElement

  ! Element counter array. This gives every triangle and every quad a
  ! unique running number (1,2,3,...)
  integer, dimension(:), intent(IN) :: IelementCounter

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(IN) :: IelIdx
  
  ! Number of corner vertices in the triangulation
  integer, intent(IN) :: NVT
  
  ! Number of edes in the triangulation
  integer, intent(IN) :: NMT
  
  ! Element type identifier for which type of elements is currently
  ! under view in IelIdx. All elements in IelIdx are assumed to be of
  ! the same type.
  ! =3: triangular, =4: quad.
  integer, intent(IN) :: NVE
!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  integer, dimension(:,:), intent(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  integer :: i
  
    if (NVE .eq. 3) then
      ! This element set consists of triangular elements.
       
      ! Loop through the elements to handle
      do i=1,size(IelIdx)
        ! Calculate the global DOF's.
        ! The P2 element has global DOF's in the corners and edge midpoints.
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
        ! Calculate the global DOF's.
        ! The Q2 element has global DOF's in the corners, edge midpoints
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
        ! elements don't contribute to DOF's in a mixed P2/Q2 discretisatrion!
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
  integer, intent(IN) :: NEL

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(IN) :: IelIdx

!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  integer, dimension(:,:), intent(OUT) :: IdofGlob

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
  integer, dimension(:,:), intent(IN) :: IedgesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(IN) :: IelIdx
  
!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  integer, dimension(:,:), intent(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  integer :: i
  
    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF's - which are simply the vertex numbers of the 
      ! corners.
      ! We always copy all elements of IedgesAtElement (:,.).
      ! There's no harm and the compiler can optimise better.
      
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
  integer, dimension(:,:), intent(IN) :: IedgesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(IN) :: IelIdx
  
!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  integer, dimension(:,:), intent(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  integer :: i
  
    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF's - which are simply the edge numbers.
      ! We always copy all elements of IedgesAtElement (:,.).
      ! There's no harm and the compiler can optimise better.
      
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
  integer, intent(IN) :: iNMT

  ! An array with the number of edges adjacent on each element of the
  ! triangulation.
  integer, dimension(:,:), intent(IN) :: IedgesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(IN) :: IelIdx
  
!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  integer, dimension(:,:), intent(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  integer :: i
  
    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF's - which are simply the numbers of the 
      ! edges. The DOF in the element gets the element number.
      ! We always copy all elements of IedgesAtElement (:,.).
      ! There's no harm and the compiler can optimise better.
      
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
  integer, intent(IN) :: iNMT

  ! An array with the number of edges adjacent on each element of the
  ! triangulation.
  integer, dimension(:,:), intent(IN) :: IedgesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(IN) :: IelIdx
  
!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  integer, dimension(:,:), intent(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  integer :: i
  
    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF's.
      ! The first 4 DOF's are the number of the edges.
      ! The next 4 DOF's are the number of the edges + nmt.
      ! The last DOF is the element number + 2*nmt.
      ! We always copy all elements of IedgesAtElement (:,.).
      ! There's no harm and the compiler can optimise better.
      
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
  integer, intent(IN) :: iNMT

  ! Number of elements in the triangulation.
  integer, intent(IN) :: iNEL

  ! An array with the number of edges adjacent on each element of the
  ! triangulation.
  integer, dimension(:,:), intent(IN) :: IedgesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(IN) :: IelIdx
  
!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  integer, dimension(:,:), intent(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  integer :: i
  
    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF's.
      ! The first 4 DOF's are the number of the edges.
      ! The next 4 DOF's are the number of the edges + nmt.
      ! The 9th DOF is at 2*nmt + element number
      ! The 10th DOF is at 2*nmt + inel + element number
      ! We always copy all elements of IedgesAtElement (:,.).
      ! There's no harm and the compiler can optimise better.
      
      IdofGlob(1:TRIA_NVEQUAD2D,i) = IedgesAtElement(1:TRIA_NVEQUAD2D,IelIdx(i))
      IdofGlob(TRIA_NVEQUAD2D+1:2*TRIA_NVEQUAD2D,i) = &
          IedgesAtElement(1:TRIA_NVEQUAD2D,IelIdx(i))+iNMT
      IdofGlob(2*TRIA_NVEQUAD2D+1,i) = 2*iNMT + IelIdx(i)
      IdofGlob(2*TRIA_NVEQUAD2D+2,i) = 2*iNMT + iNEL + IelIdx(i)
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
  integer, dimension(:), intent(IN) :: IelIdx

!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  integer, dimension(:,:), intent(OUT) :: IdofGlob

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
  integer, dimension(:,:), intent(IN) :: IverticesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(IN) :: IelIdx
  
!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  integer, dimension(:,:), intent(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  integer :: i,j
  
    ! Get the number of local DOF's - usually either 3 or 4, depending on
    ! the element. The first dimension of IdofGlob indicates the number of 
    ! DOF's.
    j = min(ubound(IverticesAtElement,1),ubound(IdofGlob,1))
    
    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF's - which are simply the vertex numbers of the 
      ! corners.
      IdofGlob(1:j,i) = IverticesAtElement(1:j,IelIdx(i))
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
  integer, intent(IN) :: NEL

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(IN) :: IelIdx

!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  integer, dimension(:,:), intent(OUT) :: IdofGlob

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
  integer, dimension(:,:), intent(IN) :: IfacesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(IN) :: IelIdx
  
!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  integer, dimension(:,:), intent(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  integer :: i
  
    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF's - which are simply the face numbers.
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
  integer, intent(IN) :: iNAT

  ! An array with the number of vertices adjacent to each element of the
  ! triangulation.
  integer, dimension(:,:), intent(IN) :: IfacesAtElement

  ! Element indices, where the mapping should be computed.
  integer, dimension(:), intent(IN) :: IelIdx
  
!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  integer, dimension(:,:), intent(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  integer :: i
  
    ! Loop through the elements to handle
    do i=1,size(IelIdx)
      ! Calculate the global DOF's - which are simply the face numbers.
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

  subroutine dof_infoDiscr (rspatialDiscr)
  
!<description>
  ! This routine prints out statistical information about a discretisation
  ! to the terminal.
!</description>

!<inputoutput>
  ! The discretisation structure where information should be printed.
  type(t_spatialDiscretisation), intent(IN), target :: rspatialDiscr
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
      call output_line ('uniform')
    case (SPDISC_CONFORMAL)
      call output_line ('conformal')
    case (SPDISC_MIXED)
      call output_line ('mixed')
    case DEFAULT
      call output_line ('undefined')
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
        // sys_siL(iand(elem_getPrimaryElement(p_relementDistr%celement),&
                   not(EL_DIMENSION)),16) )
      
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
  logical, intent(IN) :: bdetailed

  ! The discretisation structure where information should be printed.
  type(t_blockDiscretisation), intent(IN), target :: rblockDiscr
!</input>
  
!</subroutine>

    integer :: i

    call output_line ('Dimension:                    '&
        //trim(sys_siL(rblockDiscr%ndimension,10)))
    call output_line ('Complexity:                   ',bnolinebreak=.true.,&
        bnoTrim=.true.)
    select case (rblockDiscr%ccomplexity)
    case (SPDISC_UNIFORM)
      call output_line ('uniform')
    case (SPDISC_CONFORMAL)
      call output_line ('conformal')
    case (SPDISC_MIXED)
      call output_line ('mixed')
    case DEFAULT
      call output_line ('undefined')
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

end module
