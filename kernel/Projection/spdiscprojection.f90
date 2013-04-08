!##############################################################################
!# ****************************************************************************
!# <name> spdiscprojection </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines to project solution vectors between different
!# finite element spaces/discretisations. This allows e.g. to 'convert'
!# a solution vector to another one by projecting it.
!#
!# The following routines can be found here:
!#
!# 1.) spdp_projectSolutionScalar
!#     -> Convert a scalar solution vector to a new space.
!#
!# 2.) spdp_projectSolution
!#     -> Convert a block solution vector to a new space.
!#
!# 3.) spdp_stdProjectionToP1Q1Scalar
!#     -> Project a scalar solution vector to another one discretised by
!#        <tex>$P_1$</tex> and/or <tex>$Q_1$</tex> elements.
!#
!# 4.) spdp_projectToVertices
!#     -> Project a scalar vector from primal space to the vertices of the
!#        underlying mesh.
!#
!# 5.) spdp_projectToCells
!#     -> Project a scalar vector from primal space to the cells of the
!#        underlying mesh.
!# </purpose>
!##############################################################################

module spdiscprojection

!$use omp_lib
  use basicgeometry
  use derivatives
  use dofmapping
  use element
  use elementpreprocessing
  use fsystem
  use genoutput
  use linearalgebra
  use linearsystemblock
  use linearsystemscalar
  use perfconfig
  use spatialdiscretisation
  use storage
  use triangulation
  use cubature

  implicit none

  private

  public :: spdp_projectSolutionScalar
  public :: spdp_projectSolution
  public :: spdp_stdProjectionToP1Q1Scalar
  public :: spdp_projectToVertices
  public :: spdp_projectToCells
  public :: spdp_aux_getCornerRefCoords

contains

  ! ***************************************************************************

!<subroutine>
  subroutine spdp_projectSolutionScalar (rsourceVector,rdestVector,rperfconfig)

!<description>
  ! This routine 'converts' a given scalar solution vector rsourceVector to
  ! another solution vector rdestVector. The scalar discretisation structure
  ! in rdestVector specifies the new FE spaces, rsourceVector should be
  ! converted to. The new 'projected' solution is build in rdestVector.
  !
  ! Source and destination vector must be unsorted.
!</description>

!<input>
  ! The source vector to be projected.
  type(t_vectorScalar), intent(in) :: rsourceVector

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>

!<output>
  ! An existing scalar vector structure that receives the projected
  ! solution vector. Must provide a scalar discretisation structure
  ! that specifies the destination FE spaces.
  type(t_vectorScalar), intent(inout) :: rdestVector
!</output>

!</subroutine>

    ! local variables
    type(t_spatialDiscretisation), pointer :: p_rsourceDiscr,p_rdestDiscr
    type(t_triangulation), pointer :: p_rtriangulation
    real(DP), dimension(:), pointer :: p_Dsource,p_Ddest
    integer, dimension(:), pointer :: p_IelementsAtVertexIdx
    integer, dimension(:), pointer :: p_IelementsAtVertex
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    integer, dimension(:,:), pointer :: p_IedgesAtElement,&
        p_IfacesAtElement

    ! Up to now, this routine is rather rudimentary.
    ! We only support
    ! - the same triangulation + domain in the source and destination
    !   discretisation structure
    ! - only uniform discretisation structures
    ! - conversion to Q1 only (for GMV output e.g.), except both (source and
    !   destination space) are identical)
    ! - double precision vectors only
    ! - all vectors must be unsorted

    p_rsourceDiscr => rsourceVector%p_rspatialDiscr
    p_rdestDiscr => rdestVector%p_rspatialDiscr

    if (.not. associated(p_rsourceDiscr)) then
      call output_line ('No source discretisation!', &
                        OU_CLASS_ERROR,OU_MODE_STD,&
                        'spdp_projectSolutionScalar')
      call sys_halt()
    end if

    if (.not. associated(p_rdestDiscr)) then
      call output_line ('No destination discretisation!', &
                        OU_CLASS_ERROR,OU_MODE_STD,&
                        'spdp_projectSolutionScalar')
      call sys_halt()
    end if

    if (associated(p_rsourceDiscr%p_rboundary) .and.&
        associated(p_rdestDiscr%p_rboundary)) then

      if (.not. associated(p_rsourceDiscr%p_rboundary,&
                           p_rdestDiscr%p_rboundary)) then
        call output_line ('Different boundaries!', &
                          OU_CLASS_ERROR,OU_MODE_STD,&
                          'spdp_projectSolutionScalar')
        call sys_halt()
      end if
    end if

    if (.not. associated(p_rsourceDiscr%p_rtriangulation,&
                         p_rdestDiscr%p_rtriangulation)) then
      call output_line ('Different triangulations!', &
                        OU_CLASS_ERROR,OU_MODE_STD,&
                        'spdp_projectSolutionScalar')
      call sys_halt()
    end if

    if ((p_rsourceDiscr%ccomplexity .ne. SPDISC_UNIFORM) .or. &
        (p_rdestDiscr%ccomplexity .ne. SPDISC_UNIFORM)) then
      call output_line ('Only uniform discretisations supported!', &
                        OU_CLASS_ERROR,OU_MODE_STD,&
                        'spdp_projectSolutionScalar')
      call sys_halt()
    end if

    if (rsourceVector%bisSorted .or. rdestVector%bisSorted) then
      call output_line ('Vectors must be unsorted for projection!', &
                        OU_CLASS_ERROR,OU_MODE_STD,&
                        'spdp_projectSolutionScalar')
      call sys_halt()
    end if

    ! Ok, now we have a chance that we can convert.

    ! If the spaces are identical, we can simply copy the vector
    if (p_rsourceDiscr%RelementDistr(1)%celement .eq. &
        p_rdestDiscr%RelementDistr(1)%celement) then

      ! Ok, that is easy.
      ! Copy the vector data but prevent structural data from being overwritten.
      ! Let us hope the vectors have the same length :)
      ! (otherwise the copy-routine will quit)
      call lsyssc_duplicateVector (rsourceVector,rdestVector,&
          LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPYOVERWRITE)
      return

    end if

    if ((rsourceVector%cdataType .ne. ST_DOUBLE) .or. &
        (rdestVector%cdataType .ne. ST_DOUBLE)) then
      call output_line ('Only double precision vectors supported!', &
                        OU_CLASS_ERROR,OU_MODE_STD,&
                        'spdp_projectSolutionScalar')
      call sys_halt()
    end if

    ! Clear the destination vector
    call lsyssc_clearVector (rdestVector)

    ! What is the destination space?
    select case (elem_getPrimaryElement(p_rdestDiscr%RelementDistr(1)%&
                                        celement))
    case (EL_Q0, EL_P0)
      select case (elem_getPrimaryElement(p_rsourceDiscr%RelementDistr(1)%&
                                          celement))
      case (EL_Q1, EL_P1)
        ! Get geometric information from the triangulation.
        p_rtriangulation => p_rsourceDiscr%p_rtriangulation
        call storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
                                                     p_IverticesAtElement)

        ! Get the vector data
        call lsyssc_getbase_double (rsourceVector,p_Dsource)
        call lsyssc_getbase_double (rdestVector,p_Ddest)

        ! Call the conversion routine
        call spdp_Q1P1toQ0P0_DP (p_Dsource, p_Ddest, p_rtriangulation%NEL, &
                                 p_IverticesAtElement)
      case default
        ! Fallback to projection into the cells.
        call lsyssc_getbase_double (rdestVector,p_Ddest)
        call spdp_projectToCells (rsourceVector, p_Ddest, rperfconfig=rperfconfig)

      end select

    case (EL_Q1)
      ! So we should convert the source vector into a Q1 destination vector.
      ! Which element is used in the trial space?
      select case (elem_getPrimaryElement(p_rsourceDiscr%RelementDistr(1)%&
                                          celement))
      case (EL_Q0)
        ! Not too hard. Basically, take the mean of all elements adjacent to a vertex.
        !
        ! Get geometric information from the triangulation.
        p_rtriangulation => p_rsourceDiscr%p_rtriangulation
        call storage_getbase_int (p_rtriangulation%h_IelementsAtVertexIdx,&
                                                   p_IelementsAtVertexIdx)
        call storage_getbase_int (p_rtriangulation%h_IelementsAtVertex,&
                                                   p_IelementsAtVertex)

        ! Get the vector data
        call lsyssc_getbase_double (rsourceVector,p_Dsource)
        call lsyssc_getbase_double (rdestVector,p_Ddest)

        ! Call the conversion routine
        call spdp_Q0P0toQ1P1_DP (p_Dsource, p_Ddest, p_rtriangulation%NVT, &
                                 p_IelementsAtVertexIdx,p_IelementsAtVertex)

      case (EL_Q1T,EL_Q1TB,EL_Q2T,EL_Q2TB)
        ! That is a little bit harder. We have to convert an FE space with DOF`s
        ! in the midpoints to Q1. (For simplicity, the integral mean value variant
        ! is treated as if the DOF`s were in the edge midpoints. The error
        ! is negligible.
        !
        ! Get geometric information from the triangulation.
        p_rtriangulation => p_rsourceDiscr%p_rtriangulation
        call storage_getbase_int (p_rtriangulation%h_IelementsAtVertexIdx,&
                                                   p_IelementsAtVertexIdx)
        call storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
                                                     p_IverticesAtElement)
        call storage_getbase_int2d (p_rtriangulation%h_IedgesAtElement,&
                                                     p_IedgesAtElement)

        ! Get the vector data
        call lsyssc_getbase_double (rsourceVector,p_Dsource)
        call lsyssc_getbase_double (rdestVector,p_Ddest)

        ! Call the conversion routine
        call spdp_E030toQ1_DP (p_Dsource, p_Ddest, &
                               p_rtriangulation%NVT, p_rtriangulation%NEL, &
                               p_IverticesAtElement,p_IedgesAtElement,&
                               p_IelementsAtVertexIdx)

      case (EL_Q2)
        ! Rather easy. Take the first NVT elements of the Q2-vector
        ! as values in the corners of Q1.
        call lsyssc_getbase_double (rsourceVector,p_Dsource)
        call lsyssc_getbase_double (rdestVector,p_Ddest)
        call lalg_copyVector (p_Dsource(1:size(p_Ddest)),p_Ddest)

      case (EL_QP1)
        ! Also not completely trivial. Interpolation of the values on the element
        ! midpoints to the corners, neglecting the error.
        ! Get geometric information from the triangulation.
        p_rtriangulation => p_rsourceDiscr%p_rtriangulation
        call storage_getbase_int (p_rtriangulation%h_IelementsAtVertexIdx,&
                                                   p_IelementsAtVertexIdx)
        call storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
                                                     p_IverticesAtElement)

        ! Get the vector data
        call lsyssc_getbase_double (rsourceVector,p_Dsource)
        call lsyssc_getbase_double (rdestVector,p_Ddest)

        ! Call the conversion routine
        call spdp_QP1toQ1_DP (p_Dsource, p_Ddest, &
                              p_rtriangulation%NVT, p_rtriangulation%NEL, &
                              p_IverticesAtElement,&
                              p_IelementsAtVertexIdx)

      case default
        ! Fallback to projection into the vertices.
        call lsyssc_getbase_double (rdestVector,p_Ddest)
        call spdp_projectToVertices (rsourceVector, p_Ddest, rperfconfig=rperfconfig)

      end select

    case (EL_P1)
      ! So we should convert the source vector into a P1 destination vector.
      ! Which element is used in the trial space?
      select case (elem_getPrimaryElement(p_rsourceDiscr%RelementDistr(1)%&
                                          celement))
      case (EL_P0)
        ! Not too hard. Basically, take the mean of all elements adjacent to a vertex.
        !
        ! Get geometric information from the triangulation.
        p_rtriangulation => p_rsourceDiscr%p_rtriangulation
        call storage_getbase_int (p_rtriangulation%h_IelementsAtVertexIdx,&
                                                   p_IelementsAtVertexIdx)
        call storage_getbase_int (p_rtriangulation%h_IelementsAtVertex,&
                                                   p_IelementsAtVertex)

        ! Get the vector data
        call lsyssc_getbase_double (rsourceVector,p_Dsource)
        call lsyssc_getbase_double (rdestVector,p_Ddest)

        ! Call the conversion routine
        call spdp_Q0P0toQ1P1_DP (p_Dsource, p_Ddest, p_rtriangulation%NVT, &
                                 p_IelementsAtVertexIdx,p_IelementsAtVertex)

      case (EL_P1T)
        ! That is a little bit harder. We have to convert an FE space with DOF`s
        ! in the midpoints to P1. (For simplicity, the integral mean value variant
        ! is treated as if the DOF`s were in the edge midpoints. The error
        ! is negligible.
        !
        ! Get geometric information from the triangulation.
        p_rtriangulation => p_rsourceDiscr%p_rtriangulation
        call storage_getbase_int (p_rtriangulation%h_IelementsAtVertexIdx,&
                                                   p_IelementsAtVertexIdx)
        call storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
                                                     p_IverticesAtElement)
        call storage_getbase_int2d (p_rtriangulation%h_IedgesAtElement,&
                                                     p_IedgesAtElement)

        ! Get the vector data
        call lsyssc_getbase_double (rsourceVector,p_Dsource)
        call lsyssc_getbase_double (rdestVector,p_Ddest)

        ! Call the conversion routine
        call spdp_E020toP1_DP (p_Dsource, p_Ddest, &
                               p_rtriangulation%NVT, p_rtriangulation%NEL, &
                               p_IverticesAtElement,p_IedgesAtElement,&
                               p_IelementsAtVertexIdx)

      case default
        ! Fallback to projection into the vertices.
        call lsyssc_getbase_double (rdestVector,p_Ddest)
        call spdp_projectToVertices (rsourceVector, p_Ddest, rperfconfig=rperfconfig)

      end select

    case (EL_Q1_3D)
        ! So we should convert the source vector into a 3D Q1 destination vector.
      ! Which element is used in the trial space?
      select case (elem_getPrimaryElement(p_rsourceDiscr%RelementDistr(1)%&
                                          celement))
      case (EL_Q0_3D)
        ! Not too hard. Basically, take the mean of all elements adjacent to a vertex.
        !
        !
        ! Get geometric information from the triangulation.
        p_rtriangulation => p_rsourceDiscr%p_rtriangulation
        call storage_getbase_int (p_rtriangulation%h_IelementsAtVertexIdx,&
                                                   p_IelementsAtVertexIdx)
        call storage_getbase_int (p_rtriangulation%h_IelementsAtVertex,&
                                                   p_IelementsAtVertex)

        ! Get the vector data
        call lsyssc_getbase_double (rsourceVector,p_Dsource)
        call lsyssc_getbase_double (rdestVector,p_Ddest)

        ! Call the conversion routine - we can use the 2D version here...
        call spdp_Q0P0toQ1P1_DP (p_Dsource, p_Ddest, p_rtriangulation%NVT, &
                                 p_IelementsAtVertexIdx,p_IelementsAtVertex)
      case (EL_Q1T_3D)
        ! That is a little bit harder. We have to convert an FE space with DOF`s
        ! in the midpoints to Q1. (For simplicity, the integral mean value variant
        ! is treated as if the DOF`s were in the edge midpoints. The error
        ! is negligible.
        !
        ! Get geometric information from the triangulation.
        p_rtriangulation => p_rsourceDiscr%p_rtriangulation
        call storage_getbase_int (p_rtriangulation%h_IelementsAtVertexIdx,&
                                                   p_IelementsAtVertexIdx)
        call storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
                                                     p_IverticesAtElement)
        call storage_getbase_int2d (p_rtriangulation%h_IfacesAtElement,&
                                                     p_IfacesAtElement)

        ! Get the vector data
        call lsyssc_getbase_double (rsourceVector,p_Dsource)
        call lsyssc_getbase_double (rdestVector,p_Ddest)

        ! Call the conversion routine
        call spdp_E030toQ1_3D_DP (p_Dsource, p_Ddest, p_rtriangulation%NVT,&
                                  p_rtriangulation%NMT, p_rtriangulation%NEL, &
                                  p_IverticesAtElement,p_IfacesAtElement,&
                                  p_IelementsAtVertexIdx)

        case default
          call output_line ('Unsupported element in source space!', &
                            OU_CLASS_ERROR,OU_MODE_STD,&
                            'spdp_projectSolutionScalar')
          call sys_halt()

      case (EL_Q2_3D)
        ! Rather easy. Take the first NVT elements of the Q2-vector
        ! as values in the corners of Q1.
        call lsyssc_getbase_double (rsourceVector,p_Dsource)
        call lsyssc_getbase_double (rdestVector,p_Ddest)
        call lalg_copyVector (p_Dsource(1:size(p_Ddest)),p_Ddest)

      end select

    case default
      call output_line ('Unsupported element in destination space!', &
                        OU_CLASS_ERROR,OU_MODE_STD,&
                        'spdp_projectSolutionScalar')
      call sys_halt()
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine spdp_Q0P0toQ1P1_DP (Dsource, Ddest, NVT, &
                                 IelementsAtVertexIdx,IelementsAtVertex)

!<description>
  ! AUXILIARY ROUTINE.
  ! Convert a solution vector based on a uniform discretisation with
  !  Q0/P0 to a solution vector based on the Q1/P1 element.
!</description>

!<input>
  ! Source vector to be converted
  real(DP), dimension(:), intent(in) :: Dsource

  ! Number of vertices in the triangulation
  integer, intent(in) :: NVT

  ! IelementsAtVertexIdx array of the triangulation
  integer, dimension(:), intent(in) :: IelementsAtVertexIdx

  ! IelementsAtVertex array of the triangulation
  integer, dimension(:), intent(in) :: IelementsAtVertex
!</input>

!<output>
  ! Destination vector of size NVT; receives the interpolated solution
  real(DP), dimension(:), intent(out) :: Ddest
!</output>

!</subroutine>

    ! local variables
    integer :: iv
    integer :: nadj
    integer :: ielidx

    ! Loop through the vertices
    do iv=1,NVT

      ! On each vertex, loop through the adjacent elements
      do ielidx = IelementsAtVertexIdx(iv),IelementsAtVertexIdx(iv+1)-1
        ! Sum up the element contributions into the vertex
        Ddest(iv) = Ddest(iv) + Dsource(IelementsAtVertex(ielidx))
      end do

      ! Divide by the number of adjacent elements, this results
      ! in the interpolated solution.
      nadj = IelementsAtVertexIdx(iv+1) - IelementsAtVertexIdx(iv)
      Ddest(iv) = Ddest(iv) / real(nadj,DP)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine spdp_Q1P1toQ0P0_DP (Dsource, Ddest, NEL, IverticesAtElement)

!<description>
  ! AUXILIARY ROUTINE.
  ! Convert a solution vector based on a uniform discretisation with
  ! Q1/P1 to a solution vector based on the Q0/P0 element.
!</description>

!<input>
  ! Source vector to be converted
  real(DP), dimension(:), intent(in) :: Dsource

  ! Number of elements in the triangulation
  integer, intent(in) :: NEL

  ! IverticesAtElement array of the triangulation
  integer, dimension(:,:), intent(in) :: IverticesAtElement
!</input>

!<output>
  ! Destination vector of size NVT; receives the interpolated solution
  real(DP), dimension(:), intent(out) :: Ddest
!</output>

!</subroutine>

    ! local variables
    integer :: iel,iv,NVE

    ! Loop through the elements
    do iel=1,NEL

      ! On each element, loop through the vertices to get the mean value on
      ! the element
      Ddest(iel) = 0.0_DP
      NVE = tria_getNVE(IverticesAtElement,iel)
      do iv = 1, NVE
        Ddest(iel) = Ddest(iel) + Dsource(IverticesAtElement(iv,iel))
      end do

      ! Divide by the number of vertices
      Ddest(iel) = Ddest(iel) / real(NVE,DP)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine spdp_E030toQ1_DP (Dsource, Ddest, NVT, NEL, &
                               IverticesAtElement,IedgesAtElement,&
                               IelementsAtVertexIdx)

!<description>
  ! AUXILIARY ROUTINE.
  ! Convert a solution vector based on a uniform discretisation with E030, E031,
  ! EM30, EM31 to a solution vector based on the Q1 element.
!</description>

!<input>
  ! Source vector to be converted
  real(DP), dimension(:), intent(in) :: Dsource

  ! Number of vertices in the triangulation
  integer, intent(in) :: NVT

  ! Number of elements in the triangulation
  integer, intent(in) :: NEL

  ! IelementsAtVertexIdx array of the triangulation
  integer, dimension(:), intent(in) :: IelementsAtVertexIdx

  ! IverticesAtElement array of the triangulation (old KVERT)
  integer, dimension(:,:), intent(in) :: IverticesAtElement

  ! IedgesAtElement array of the triangulation (old KMID)
  integer, dimension(:,:), intent(in) :: IedgesAtElement
!</input>

!<output>
  ! Destination vector of size NVT; receives the interpolated solution
  real(DP), dimension(:), intent(out) :: Ddest
!</output>

!</subroutine>

    ! local variables
    integer :: iv
    integer :: iel
    integer :: IM1,IM2,IM3,IM4
    integer :: IV1,IV2,IV3,IV4
    real(DP) :: DUH1,DUH2,DUH3,DUH4
    integer :: nadj

    ! Clear the output array
    call lalg_clearVector (Ddest)

    ! Loop through the elements
    do iel=1,NEL

      ! Get the global DOF`s on the current element in the E030 space
      IM1 = IedgesAtElement(1,iel)
      IM2 = IedgesAtElement(2,iel)
      IM3 = IedgesAtElement(3,iel)
      IM4 = IedgesAtElement(4,iel)

      ! Get the global DOF`s on the current element in the Q1 space
      IV1 = IverticesAtElement(1,iel)
      IV2 = IverticesAtElement(2,iel)
      IV3 = IverticesAtElement(3,iel)
      IV4 = IverticesAtElement(4,iel)

      ! Get the values of the DOF`s in the E030 space
      DUH1 = Dsource(IM1)
      DUH2 = Dsource(IM2)
      DUH3 = Dsource(IM3)
      DUH4 = Dsource(IM4)

      ! Bilinear interpolation gives what we have to add to the
      ! value in each corner:
      Ddest(IV1) = Ddest(IV1) + 0.75_DP*(DUH1+DUH4) - 0.25_DP*(DUH2+DUH3)
      Ddest(IV2) = Ddest(IV2) + 0.75_DP*(DUH2+DUH1) - 0.25_DP*(DUH3+DUH4)
      Ddest(IV3) = Ddest(IV3) + 0.75_DP*(DUH3+DUH2) - 0.25_DP*(DUH4+DUH1)
      Ddest(IV4) = Ddest(IV4) + 0.75_DP*(DUH4+DUH3) - 0.25_DP*(DUH1+DUH2)

    end do

    ! Divide by the number of adjacent elements, this results
    ! in the interpolated solution.
    do iv=1,NVT
      nadj = IelementsAtVertexIdx(iv+1) - IelementsAtVertexIdx(iv)
      Ddest(iv) = Ddest(iv) / real(nadj,DP)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdp_E037toQ1_DP (Dsource, Ddest, NVT, NEL, &
                               IverticesAtElement,IedgesAtElement,&
                               IelementsAtVertexIdx,ItwistIndexEdges)

!<description>
  ! AUXILIARY ROUTINE.
  ! Convert a solution vector based on a uniform discretisation with E030, E031,
  ! EM30, EM31 to a solution vector based on the Q1 element.
!</description>

!<input>
  ! Source vector to be converted
  real(DP), dimension(:), intent(in) :: Dsource

  ! Number of vertices in the triangulation
  integer, intent(in) :: NVT

  ! Number of elements in the triangulation
  integer, intent(in) :: NEL

  ! IelementsAtVertexIdx array of the triangulation
  integer, dimension(:), intent(in) :: IelementsAtVertexIdx

  ! IverticesAtElement array of the triangulation (old KVERT)
  integer, dimension(:,:), intent(in) :: IverticesAtElement

  ! IedgesAtElement array of the triangulation (old KMID)
  integer, dimension(:,:), intent(in) :: IedgesAtElement

  ! ItwistIndexEdges array of the triangulation
  integer(I32), dimension(:), intent(in) :: ItwistIndexEdges
!</input>

!<output>
  ! Destination vector of size NVT; receives the interpolated solution
  real(DP), dimension(:), intent(out) :: Ddest
!</output>

!</subroutine>

    ! local variables
    integer :: iv
    integer :: iel
    integer :: IM1,IM2,IM3,IM4
    integer :: IV1,IV2,IV3,IV4
    real(DP) :: DUH1,DUH2,DUH3,DUH4
    integer :: nadj

    ! Clear the output array
    call lalg_clearVector (Ddest)

    ! Loop through the elements
    do iel=1,NEL

      ! Get the global DOF`s on the current element in the E030 space
      IM1 = IedgesAtElement(1,iel)
      IM2 = IedgesAtElement(2,iel)
      IM3 = IedgesAtElement(3,iel)
      IM4 = IedgesAtElement(4,iel)

      ! Get the global DOF`s on the current element in the Q1 space
      IV1 = IverticesAtElement(1,iel)
      IV2 = IverticesAtElement(2,iel)
      IV3 = IverticesAtElement(3,iel)
      IV4 = IverticesAtElement(4,iel)

      ! Get the values of the DOF`s in the E030 space
      DUH1 = Dsource(IM1)
      DUH2 = Dsource(IM2)
      DUH3 = Dsource(IM3)
      DUH4 = Dsource(IM4)

      ! Bilinear interpolation gives what we have to add to the
      ! value in each corner:
      Ddest(IV1) = Ddest(IV1) + 0.75_DP*(DUH1+DUH4) - 0.25_DP*(DUH2+DUH3)
      Ddest(IV2) = Ddest(IV2) + 0.75_DP*(DUH2+DUH1) - 0.25_DP*(DUH3+DUH4)
      Ddest(IV3) = Ddest(IV3) + 0.75_DP*(DUH3+DUH2) - 0.25_DP*(DUH4+DUH1)
      Ddest(IV4) = Ddest(IV4) + 0.75_DP*(DUH4+DUH3) - 0.25_DP*(DUH1+DUH2)

    end do

    ! Divide by the number of adjacent elements, this results
    ! in the interpolated solution.
    do iv=1,NVT
      nadj = IelementsAtVertexIdx(iv+1) - IelementsAtVertexIdx(iv)
      Ddest(iv) = Ddest(iv) / real(nadj,DP)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine spdp_QP1toQ1_DP (Dsource, Ddest, NVT, NEL, &
                              IverticesAtElement,&
                              IelementsAtVertexIdx)

!<description>
  ! AUXILIARY ROUTINE.
  ! Convert a solution vector based on a uniform discretisation with QP1,
  ! to a solution vector based on the Q1 element.
  !
  ! Linear interpolation of values in the midpoints of the elements.
!</description>

!<input>
  ! Source vector to be converted
  real(DP), dimension(:), intent(in) :: Dsource

  ! Number of vertices in the triangulation
  integer, intent(in) :: NVT

  ! Number of elements in the triangulation
  integer, intent(in) :: NEL

  ! IelementsAtVertexIdx array of the triangulation
  integer, dimension(:), intent(in) :: IelementsAtVertexIdx

  ! IverticesAtElement array of the triangulation (old KVERT)
  integer, dimension(:,:), intent(in) :: IverticesAtElement

!</input>

!<output>
  ! Destination vector of size NVT; receives the interpolated solution
  real(DP), dimension(:), intent(out) :: Ddest
!</output>

!</subroutine>

    ! local variables
    integer :: iv
    integer :: iel
    integer :: IV1,IV2,IV3,IV4
    integer :: nadj

    ! Clear the output array
    call lalg_clearVector (Ddest)

    ! Loop through the elements
    do iel=1,NEL

      ! Get the global DOF`s on the current element in the Q1 space
      IV1 = IverticesAtElement(1,iel)
      IV2 = IverticesAtElement(2,iel)
      IV3 = IverticesAtElement(3,iel)
      IV4 = IverticesAtElement(4,iel)

      ! Get the value in the midpoint of the current element
      ! and add it to all corners.
      Ddest(IV1) = Ddest(IV1) + Dsource(iel)
      Ddest(IV2) = Ddest(IV2) + Dsource(iel)
      Ddest(IV3) = Ddest(IV3) + Dsource(iel)
      Ddest(IV4) = Ddest(IV4) + Dsource(iel)

    end do

    ! Divide by the number of adjacent elements, this results
    ! in the interpolated solution.
    do iv=1,NVT
      nadj = IelementsAtVertexIdx(iv+1) - IelementsAtVertexIdx(iv)
      Ddest(iv) = Ddest(iv) / real(nadj,DP)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine spdp_E030toQ1_3D_DP (Dsource, Ddest, NVT, NMT, NEL, &
                                  IverticesAtElement,IfacesAtElement,&
                                  IelementsAtVertexIdx)

!<description>
  ! AUXILIARY ROUTINE.
  ! Convert a solution vector based on a uniform discretisation with E030, E031,
  ! EM30, EM31 to a solution vector based on the Q1 element.
!</description>

!<input>
  ! Source vector to be converted
  real(DP), dimension(:), intent(in) :: Dsource

  ! Number of vertices in the triangulation
  integer, intent(in) :: NVT

  ! Number of edges in the triangulation
  integer, intent(in) :: NMT

  ! Number of elements in the triangulation
  integer, intent(in) :: NEL

  ! IelementsAtVertexIdx array of the triangulation
  integer, dimension(:), intent(in) :: IelementsAtVertexIdx

  ! IverticesAtElement array of the triangulation (old KVERT)
  integer, dimension(:,:), intent(in) :: IverticesAtElement

  ! IfacesAtElement array of the triangulation
  integer, dimension(:,:), intent(in) :: IfacesAtElement
!</input>

!<output>
  ! Destination vector of size NVT; receives the interpolated solution
  real(DP), dimension(:), intent(out) :: Ddest
!</output>

!</subroutine>

    ! local variables
    integer :: iv
    integer :: iel
    integer, dimension(6) :: F
    integer, dimension(8) :: V
    real(DP), dimension(6) :: D
    real(DP),parameter :: R13 = 0.333333333333333_DP
    real(DP),parameter :: R23 = 0.666666666666667_DP
    integer :: nadj

    ! Clear the output array
    call lalg_clearVector (Ddest)

    ! Loop through the elements
    do iel=1,NEL

      ! Get the global DOF`s on the current element in the E030 space
      F(1:6) = IfacesAtElement(1:6,iel)

      ! Get the global DOF`s on the current element in the Q1 space
      V(1:8) = IverticesAtElement(1:8,iel)

      ! Get the values of the DOF`s in the E030 space
      D(1:6) = Dsource(F(1:6))

      ! Bilinear interpolation gives what we have to add to the
      ! value in each corner:
      Ddest(V(1)) = Ddest(V(1)) - R13*(D(3)+D(4)+D(6)) + R23*(D(1)+D(2)+D(5))
      Ddest(V(2)) = Ddest(V(2)) - R13*(D(4)+D(5)+D(6)) + R23*(D(1)+D(2)+D(3))
      Ddest(V(3)) = Ddest(V(3)) - R13*(D(2)+D(5)+D(6)) + R23*(D(1)+D(3)+D(4))
      Ddest(V(4)) = Ddest(V(4)) - R13*(D(2)+D(3)+D(6)) + R23*(D(1)+D(4)+D(5))
      Ddest(V(5)) = Ddest(V(5)) - R13*(D(1)+D(3)+D(4)) + R23*(D(2)+D(5)+D(6))
      Ddest(V(6)) = Ddest(V(6)) - R13*(D(1)+D(4)+D(5)) + R23*(D(2)+D(3)+D(6))
      Ddest(V(7)) = Ddest(V(7)) - R13*(D(1)+D(2)+D(5)) + R23*(D(3)+D(4)+D(6))
      Ddest(V(8)) = Ddest(V(8)) - R13*(D(1)+D(2)+D(3)) + R23*(D(4)+D(5)+D(6))

    end do

    ! Divide by the number of adjacent elements, this results
    ! in the interpolated solution.
    do iv=1,NVT
      nadj = IelementsAtVertexIdx(iv+1) - IelementsAtVertexIdx(iv)
      Ddest(iv) = Ddest(iv) / real(nadj,DP)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine spdp_E020toP1_DP (Dsource, Ddest, NVT, NEL, &
                               IverticesAtElement,IedgesAtElement,&
                               IelementsAtVertexIdx)

!<description>
  ! AUXILIARY ROUTINE.
  ! Convert a solution vector based on a uniform discretisation with E020
  ! to a solution vector based on the P1 element.
!</description>

!<input>
  ! Source vector to be converted
  real(DP), dimension(:), intent(in) :: Dsource

  ! Number of vertices in the triangulation
  integer, intent(in) :: NVT

  ! Number of elements in the triangulation
  integer, intent(in) :: NEL

  ! IelementsAtVertexIdx array of the triangulation
  integer, dimension(:), intent(in) :: IelementsAtVertexIdx

  ! IverticesAtElement array of the triangulation (old KVERT)
  integer, dimension(:,:), intent(in) :: IverticesAtElement

  ! IedgesAtElement array of the triangulation (old KMID)
  integer, dimension(:,:), intent(in) :: IedgesAtElement
!</input>

!<output>
  ! Destination vector of size NVT; receives the interpolated solution
  real(DP), dimension(:), intent(out) :: Ddest
!</output>

!</subroutine>

    ! local variables
    integer :: iv
    integer :: iel
    integer :: IM1,IM2,IM3
    integer :: IV1,IV2,IV3
    real(DP) :: DUH1,DUH2,DUH3
    integer :: nadj

    ! Clear the output array
    call lalg_clearVector (Ddest)

    ! Loop through the elements
    do iel=1,NEL

      ! Get the global DOF`s on the current element in the E020 space
      IM1 = IedgesAtElement(1,iel)
      IM2 = IedgesAtElement(2,iel)
      IM3 = IedgesAtElement(3,iel)

      ! Get the global DOF`s on the current element in the P1 space
      IV1 = IverticesAtElement(1,iel)
      IV2 = IverticesAtElement(2,iel)
      IV3 = IverticesAtElement(3,iel)

      ! Get the values of the DOF`s in the E020 space
      DUH1 = Dsource(IM1)
      DUH2 = Dsource(IM2)
      DUH3 = Dsource(IM3)

      ! Linear interpolation gives what we have to add to the
      ! value in each corner:
      Ddest(IV1) = Ddest(IV1) + DUH1 - DUH2 + DUH3
      Ddest(IV2) = Ddest(IV2) + DUH1 + DUH2 - DUH3
      Ddest(IV3) = Ddest(IV3) - DUH1 + DUH2 + DUH3

    end do

    ! Divide by the number of adjacent elements, this results
    ! in the interpolated solution.
    do iv=1,NVT
      nadj = IelementsAtVertexIdx(iv+1) - IelementsAtVertexIdx(iv)
      Ddest(iv) = Ddest(iv) / real(nadj,DP)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine spdp_projectSolution (rsourceVector,rdestVector,rperfconfig)

!<description>
  ! This routine 'converts' a given solution vector rsourceVector to
  ! another solution vector rdestVector. The discretisation structures
  ! in the subvertors of rdestVector specifies the new FE spaces, rsourceVector
  ! should be converted to. The new 'projected' solution is build in
  ! rdestVector.
!</description>

!<input>
  ! The source vector to be projected.
  type(t_vectorBlock), intent(in) :: rsourceVector

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>

!<output>
  ! An existing vector structure that receives the projected
  ! solution vector. Must provide a scalar discretisation structures
  ! that specifies the destination FE spaces.
  type(t_vectorBlock), intent(inout) :: rdestVector
!</output>

!</subroutine>

    ! local variables
    integer :: i

    if (rsourceVector%nblocks .ne. rdestVector%nblocks) then
      call output_line ('Different block structure!', &
                        OU_CLASS_ERROR,OU_MODE_STD,&
                        'spdp_projectSolution')
      call sys_halt()
    end if

    ! Apply spdp_projectSolutionScalar to every subvector, that is all
    do i=1,rsourceVector%nblocks
      call spdp_projectSolutionScalar (rsourceVector%RvectorBlock(i),&
                                       rdestVector%RvectorBlock(i), rperfconfig)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine spdp_stdProjectionToP1Q1Scalar (rsourceVector,&
      rdestVector,rdestDiscretisation,rperfconfig)

!<description>
  ! This routine 'converts' a given scalar solution vector rsourceVector to
  ! another solution vector rdestVector. If necessary, the destination
  ! vector is allocated.
  ! The source vector can be an arbitrary FE solution vector. The destination
  ! vector will be a solution vector in the <tex>$Q_1$</tex> space (for quad
  ! elements) and in the $P_1$ space (for triangular elements), respectively.
  ! (Such vertex based solution vectors are typically written out
  ! to external files like GMV during the postprocessing).
  !
  ! If rdestDiscretisation is defined, it must describe a discretisation
  ! with <tex>$P_1$</tex> and <tex>$Q_1$</tex> elements, respectively.
  ! If rdestDiscretisation is undefined, the routine automatically
  ! creates and returns rdestDiscretisation for the destination vector
  ! rdestVector based on the discretisation of the source vector.
  ! The caller must take care, that this discretisation structure
  ! is destroyed when rdestVector is destroyed to avoid memory leaks!
  !
  ! (-> Can be used to project multiple vectors at once:
  !  The first call creates rdestDiscretisation, the next calls use it
  !  and afterwards the application destroys it together with the vector.)
  !
  ! Source and destination vector must be unsorted.
!</description>

!<input>
  ! The source vector to be projected.
  type(t_vectorScalar), intent(in) :: rsourceVector

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>

!<inputoutput>
  ! A scalar vector structure that receives the projected
  ! solution vector. If undefined, a vector is created.
  type(t_vectorScalar), intent(inout) :: rdestVector

  ! A discretisation structure that defines a discretisation with
  ! $P_1$ and/or <tex>$Q_1$</tex> elements. If undefines, the structure
  ! is automatically created.
  type(t_spatialDiscretisation), intent(inout) :: rdestDiscretisation
!</inputoutput>

!</subroutine>

    ! Check that the source discretisation structure is valid.
    if (rsourceVector%p_rspatialDiscr%ndimension .ne. NDIM2D) then
      call output_line ('Only 2D discretisation supported!', &
                        OU_CLASS_ERROR,OU_MODE_STD,&
                        'spdp_stdProjectionToP1Q1')
      call sys_halt()
    end if

    ! Is the destination discretisation given?
    if (rdestDiscretisation%ndimension .eq. 0) then

      ! Create a discretisation with $P_1$ and/or $Q_1$ elements.
      ! Derive it from the existing discretisation structure in
      ! the vector, so the element information does not use too
      ! much memory...

      call spdiscr_deriveDiscr_triquad(&
          rsourceVector%p_rspatialDiscr,&
          EL_P1,EL_Q1,rdestDiscretisation)

    end if

    ! Is the destination vector given?
    if (rdestVector%NEQ .eq. 0) then
      call lsyssc_createVecByDiscr (&
          rdestDiscretisation,rdestVector,.true.)
    end if

    ! Project the solution to the $P_1$ / $Q_1$ space.
    call spdp_projectSolutionScalar (rsourceVector,rdestVector, rperfconfig)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdp_projectToVertices (rvector, p_Dvalues, ideriv, rperfconfig)

!<description>
  ! This routines projects a vector from primal space to the vertices of the
  ! mesh, i.e. the FE function which is described by the coefficient vector
  ! is evaluated in the vertices of its underlying mesh.
  ! The coefficient vector must be unsorted.
!</description>

!<input>
  ! The coefficient vector to be projected.
  type(t_vectorScalar), intent(in) :: rvector

  ! OPTIONAL: A derivative quantifier specifying which derivative is to be
  ! projected onto the vertices. If not given, DER_FUNC is used.
  integer, optional, intent(in) :: ideriv

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>

!<inputoutput>
  ! An array that recieves the values of the FE function in the vertices of
  ! the mesh. If set to NULL on entry, an array of the correct size is
  ! allocated.
  ! For vector-valued FEM spaces, the array will/must be of size n*NVT
  ! with n the dimension of the FE space. At 1..NVT, the first coordinate
  ! can be found, at NVT+1..2*NVT the 2nd coordinate etc.
  real(DP), dimension(:), pointer :: p_Dvalues
!</inputoutput>

!</subroutine>

  ! A hand full of local variables
  type(t_spatialDiscretisation), pointer :: p_rdiscr
  type(t_triangulation), pointer :: p_rtria
  type(t_elementDistribution), pointer :: p_relemDist
  integer, dimension(:), pointer :: p_IelemList, p_IelemAtVertIdx, p_IcurEL
  integer, dimension(:,:), pointer :: p_IvertAtElem
  real(DP), dimension(4,TRIA_MAXNVE) :: Dcorners
  integer :: NEL,NVT,NVE,NDIM,NBAS,NDER,ied,ivt,iel,i,j,k
  type(t_evalElementSet)  :: reval
  integer(I32) :: cevalTag, celement, ctrafo
  integer :: NELtodo, NELdone, NELpatch
  logical, dimension(EL_MAXNDER) :: Bder
  real(DP), dimension(:,:,:,:), allocatable :: Dbas
  integer, dimension(:,:), allocatable :: Idofs
  real(DP), dimension(:), pointer :: p_Dx
  real(DP) :: dt
  integer :: ider,nadj,idimfe,ndimfe
  logical :: bvertexDofs,ballocated

    ! Which derivative?
    ider = DER_FUNC
    if(present(ideriv)) ider = ideriv

    ! First of all, get the spatial discretisation
    p_rdiscr => rvector%p_rspatialDiscr

    if(.not. associated(p_rdiscr)) then
      call output_line ('Vector does not have a discretisation!', &
                        OU_CLASS_ERROR,OU_MODE_STD, 'spdp_projectToVertices')
      call sys_halt()
    end if

    ! Get the coefficient vector`s data array
    call lsyssc_getbase_double(rvector, p_Dx)

    ! Now get the triangulation
    p_rtria => p_rdiscr%p_rtriangulation

    ! Get the dimension and the number of vertices
    NVT = p_rtria%NVT
    NDIM = p_rtria%ndim

    ! Project function values?
    if(ider .eq. DER_FUNC) then

      ! Let' see if our FE space has DOFs in the vertices
      bvertexDofs = .true.
      do i = 1, p_rdiscr%inumFESpaces

        ! check the elements
        select case(elem_getPrimaryElement(p_rdiscr%RelementDistr(i)%celement))
        case (EL_P1_1D, EL_P2_1D, EL_S31_1D)
        case (EL_P1_2D, EL_P2_2D)
        case (EL_Q1_2D, EL_Q2_2D)
        case (EL_P1_3D)
        case (EL_Q1_3D, EL_Q2_3D)
          ! For all of the previous elements we can simply copy the first NVT entries
          ! from the coefficient vector to obtain the function values in the vertices.

        case default
          ! Some other FE space with no function values in the vertices
          bvertexDofs = .false.

        end select

      end do

      ! Vertex-DOF based FE space?
      if(bvertexDofs) then

        ! Prepare the output array
        if(.not. associated(p_Dvalues)) then
          allocate(p_Dvalues(NVT))
        else if(ubound(p_Dvalues,1) .lt. NVT) then
          deallocate(p_Dvalues)
          allocate(p_Dvalues(NVT))
        end if

        ! Then we only need to copy the first NVT DOFs from our vector
        call lalg_copyVector(p_Dx, p_Dvalues, NVT)
        return

      end if

    end if

    ! Get the vertices-at-element and elements-at-vertex-idx arrays
    call storage_getbase_int2d(p_rtria%h_IverticesAtElement, p_IvertAtElem)
    call storage_getbase_int(p_rtria%h_IelementsAtVertexIdx, p_IelemAtVertIdx)

    ! Set up Bder
    Bder = .false.
    Bder(ider) = .true.
    
    ballocated = .false.
    ndimfe = 0

    ! Okay, now loop through all element distributions
    do ied = 1, p_rdiscr%inumFESpaces

      ! Get the element distribution
      p_relemDist => p_rdiscr%RelementDistr(ied)

      ! Get the element list for this distribution
      call storage_getbase_int(p_relemDist%h_IelementList, p_IelemList)

      ! Get the element, its evaluation tag and trafo type
      celement = p_relemDist%celement
      cevalTag = elem_getEvaluationTag(celement)
      ctrafo = elem_igetTrafoType(celement)

      ! Get the number of vertices per element
      NVE = tria_getNVE(p_rtria, p_IelemList(1))

      ! Get the number of basis functions and derivatives
      NBAS = elem_igetNDofLoc(celement)
      NDER = elem_getMaxDerivative(celement)

      ! Make sure ider is supported
      if(ider .gt. NDER) then
        call output_line ('Derivative not supported by element!', &
            OU_CLASS_ERROR,OU_MODE_STD, 'spdp_projectToVertices')
        call sys_halt()
      end if

      if (.not. ballocated) then

        ! Get the dimension of the FE space
        ndimfe = elem_igetFeDimension(celement)

        ! Prepare the output array
        if(.not. associated(p_Dvalues)) then
          allocate(p_Dvalues(NVT*ndimfe))
        else if(ubound(p_Dvalues,1) .lt. NVT*ndimfe) then
          deallocate(p_Dvalues)
          allocate(p_Dvalues(NVT*ndimfe))
        end if

        ! Clear the output array
        call lalg_clearVector(p_Dvalues,NVT*ndimfe)

        ballocated = .true.
      
      end if

      ! Calculate the corner vertice reference coordinates
      call spdp_aux_getCornerRefCoords(Dcorners, NDIM, NVE)

      ! Get the number of elements in this distribution
      NEL = size(p_IelemList)

      ! Determine the element patch size
      NELpatch = min(1000, NEL)

      ! Allocate evaluation array
      allocate(Dbas(NBAS*ndimfe,NDER,NVE,NELpatch))

      ! Allocate DOF-mapping array
      allocate(Idofs(NBAS,NELpatch))

      ! Loop over the element patches
      NELdone = 0
      do while(NELdone .lt. NEL)

        ! How many elements do we process this time?
        NELtodo = min(NELpatch, NEL-NELdone)

        ! Get a pointer to the current element list
        p_IcurEL => p_IelemList(NELdone+1:NELdone+NELtodo)

        ! Prepare the element for evaluation
        call elprep_prepareSetForEvaluation(reval, cevalTag, p_rtria, &
            p_IcurEL, ctrafo, Dcorners(:,1:NVE), rperfconfig=rperfconfig)

        ! Evaluate the element
        call elem_generic_sim2(celement, reval, Bder, Dbas)

        ! Perform DOF-mapping
        call dof_locGlobMapping_mult(p_rdiscr, p_IcurEL, Idofs(:,1:NELtodo))

        ! Now loop over all elements
        do i = 1, NELtodo

          ! Get the index of the current element
          iel = p_IcurEL(i)

          ! Loop over all vertices of the element
          do j = 1, NVE

            ! Get the index of the vertice
            ivt = p_IvertAtElem(j,iel)
            if(ivt .le. 0) cycle

            ! Loop over the dimension of the FE space
            do idimfe = 0,ndimfe-1
              ! Okay, loop over all basis functions and calculate the value
              dt = 0.0_DP
              do k = 1, NBAS
                dt = dt + Dbas(k+idimfe*NBAS,ider,j,i)*p_Dx(Idofs(k,i))
              end do ! k

              ! Add the value onto the vertice
              p_Dvalues(ivt+idimfe*NVT) = p_Dvalues(ivt+idimfe*NVT) + dt
            end do

          end do ! j

        end do ! i

        ! Release element evaluation set
        call elprep_releaseElementSet(reval)

        ! Go for the next element patch
        NELdone = NELdone + NELtodo

      end do

      ! Deallocate arrays
      deallocate(Idofs)
      deallocate(Dbas)

      ! Go for the next element distribution

    end do ! ied

    do idimfe = 0,ndimfe-1
      ! And loop through all vertices
      do ivt = 1, NVT

        ! Calculate the number of elements adjacent to this vertice
        nadj = p_IelemAtVertIdx(ivt+1) - p_IelemAtVertIdx(ivt)

        ! And divide the value in this vertice by the number of elements adjacent
        ! to this vertice.
        if(nadj .gt. 1) then
          p_Dvalues(ivt+idimfe*NVT) = p_Dvalues(ivt+idimfe*NVT) / real(nadj,DP)
        end if

      end do
    end do

    ! That is it

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine spdp_projectToCells (rvector, p_Dvalues, ideriv, rperfconfig)

!<description>
  ! This routines projects a vector from primal space to the cells of the
  ! mesh, i.e. the FE function which is described by the coefficient vector
  ! is evaluated in the element midpoints of its underlying mesh.
  ! The coefficient vector must be unsorted.
!</description>

!<input>
  ! The coefficient vector to be projected.
  type(t_vectorScalar), intent(in) :: rvector

  ! OPTIONAL: A derivative quantifier specifying which derivative is to be
  ! projected onto the vertices. If not given, DER_FUNC is used.
  integer, optional, intent(in) :: ideriv

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>

!<inputoutput>
  ! An array that recieves the values of the FE function in the cells of
  ! the mesh. If set to NULL on entry, an array of the correct size is
  ! allocated.
  ! For vector-valued FEM spaces, the array will/must be of size n*NEL
  ! with n the dimension of the FE space. At 1..NEL, the first coordinate
  ! can be found, at NEL+1..2*NEL the 2nd coordinate etc.
  real(DP), dimension(:), pointer :: p_Dvalues
!</inputoutput>

!</subroutine>

  ! A hand full of local variables
  type(t_spatialDiscretisation), pointer :: p_rdiscr
  type(t_triangulation), pointer :: p_rtria
  type(t_elementDistribution), pointer :: p_relemDist
  integer, dimension(:), pointer :: p_IelemList, p_IcurEL
  real(DP), dimension(4,1) :: DmidPoint
  integer :: NEL,NVE,NDIM,NBAS,NDER,ied,iel,i,j
  type(t_evalElementSet)  :: reval
  integer(I32) :: cevalTag, celement, ctrafo
  integer :: NELtodo, NELdone, NELpatch,NELinList
  logical, dimension(EL_MAXNDER) :: Bder
  real(DP), dimension(:,:,:,:), allocatable :: Dbas
  integer, dimension(:,:), allocatable :: Idofs
  real(DP), dimension(:), pointer :: p_Dx
  real(DP) :: dt
  integer :: ider,idimfe,ndimfe
  logical :: bcellDofs,ballocated

    ! Which derivative?
    ider = DER_FUNC
    if(present(ideriv)) ider = ideriv

    ! First of all, get the spatial discretisation
    p_rdiscr => rvector%p_rspatialDiscr

    if(.not. associated(p_rdiscr)) then
      call output_line ('Vector does not have a discretisation!', &
                        OU_CLASS_ERROR,OU_MODE_STD, 'spdp_projectToCells')
      call sys_halt()
    end if

    ! Get the coefficient vector`s data array
    call lsyssc_getbase_double(rvector, p_Dx)

    ! Now get the triangulation
    p_rtria => p_rdiscr%p_rtriangulation

    ! Get the dimension
    NDIM = p_rtria%ndim
    NEL = p_rtria%NEL

    ! Project function values?
    if(ider .eq. DER_FUNC) then

      ! Let' see if our FE space has DOFs in the cells
      bcellDofs = .true.
      do i = 1, p_rdiscr%inumFESpaces

        ! check the elements
        select case(elem_getPrimaryElement(p_rdiscr%RelementDistr(i)%celement))
        case (EL_P0_1D)
        case (EL_P0_2D)
        case (EL_Q0_2D, EL_QP1_2D)
        case (EL_P0_3D)
        case (EL_Q0_3D, EL_QP1_3D)
          ! For all of the previous elements we can simply copy the first NEL entries
          ! from the coefficient vector to obtain the function values in the cells.

        case default
          ! Some other FE space with no function values in the cells
          bcellDofs = .false.

        end select

      end do

      ! Cell-DOF based FE space?
      if(bcellDofs) then

        ! Prepare the output array
        if(.not. associated(p_Dvalues)) then
          allocate(p_Dvalues(NEL))
        else if(ubound(p_Dvalues,1) .lt. NEL) then
          deallocate(p_Dvalues)
          allocate(p_Dvalues(NEL))
        end if

        ! Then we only need to copy the first NVT DOFs from our vector
        call lalg_copyVector(p_Dx, p_Dvalues, NEL)
        return

      end if

    end if

    ! Set up Bder
    Bder = .false.
    Bder(ider) = .true.
    
    ballocated = .false.
    ndimfe = 0

    ! Okay, now loop through all element distributions
    do ied = 1, p_rdiscr%inumFESpaces

      ! Get the element distribution
      p_relemDist => p_rdiscr%RelementDistr(ied)

      ! Get the element list for this distribution
      call storage_getbase_int(p_relemDist%h_IelementList, p_IelemList)

      ! Get the element, its evaluation tag and trafo type
      celement = p_relemDist%celement
      cevalTag = elem_getEvaluationTag(celement)
      ctrafo = elem_igetTrafoType(celement)

      ! Get the number of vertices per element
      NVE = tria_getNVE(p_rtria, p_IelemList(1))

      ! Get the number of basis functions and derivatives
      NBAS = elem_igetNDofLoc(celement)
      NDER = elem_getMaxDerivative(celement)
      
      ! Make sure ider is supported
      if(ider .gt. NDER) then
        call output_line ('Derivative not supported by element!', &
            OU_CLASS_ERROR,OU_MODE_STD, 'spdp_projectToCells')
        call sys_halt()
      end if
      
      if (.not. ballocated) then

        ! Get the Dimension of the FE space
        ndimfe = elem_igetFeDimension(celement)

        ! Prepare the output array
        if(.not. associated(p_Dvalues)) then
          allocate(p_Dvalues(NEL*ndimfe))
        else if(ubound(p_Dvalues,1) .lt. NEL*ndimfe) then
          deallocate(p_Dvalues)
          allocate(p_Dvalues(NEL*ndimfe))
        end if

        ! Clear the output array
        call lalg_clearVector(p_Dvalues,NEL*ndimfe)
        
        ballocated = .true.
      
      end if

      ! Calculate the cell midpoint reference coordinates
      call spdp_aux_getMidpointRefCoords(DmidPoint, NDIM, NVE)

      ! Get the number of elements in this distribution
      NELinList = size(p_IelemList)

      ! Determine the element patch size
      NELpatch = min(1000, NELinList)

      ! Allocate evaluation array
      allocate(Dbas(NBAS*ndimfe,NDER,1,NELpatch))

      ! Allocate DOF-mapping array
      allocate(Idofs(NBAS,NELpatch))

      ! Loop over the element patches
      NELdone = 0
      do while(NELdone .lt. NELinList)

        ! How many elements do we process this time?
        NELtodo = min(NELpatch, NELinList-NELdone)

        ! Get a pointer to the current element list
        p_IcurEL => p_IelemList(NELdone+1:NELdone+NELtodo)

        ! Prepare the element for evaluation
        call elprep_prepareSetForEvaluation(reval, cevalTag, p_rtria, &
            p_IcurEL, ctrafo, DmidPoint, rperfconfig=rperfconfig)

        ! Evaluate the element
        call elem_generic_sim2(celement, reval, Bder, Dbas)

        ! Perform DOF-mapping
        call dof_locGlobMapping_mult(p_rdiscr, p_IcurEL, Idofs(:,1:NELtodo))

        ! Now loop over all elements
        do i = 1, NELtodo

          ! Get the index of the current element
          iel = p_IcurEL(i)

          ! Loop through the dimension of the FE space
          do idimfe = 0,ndimfe-1
            ! Okay, loop over all basis functions and calculate the value
            dt = 0.0_DP
            do j = 1, NBAS
              dt = dt + Dbas(j+idimfe*NBAS,ider,1,i)*p_Dx(Idofs(j,i))
            end do ! k

            ! Add the value onto the cell
            p_Dvalues(iel + idimfe*NEL) = p_Dvalues(iel + idimfe*NEL) + dt
          end do

        end do ! i

        ! Release element evaluation set
        call elprep_releaseElementSet(reval)

        ! Go for the next element patch
        NELdone = NELdone + NELtodo

      end do

      ! Deallocate arrays
      deallocate(Idofs)
      deallocate(Dbas)

      ! Go for the next element distribution

    end do ! ied

    ! That is it

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine spdp_aux_getCornerRefCoords(Dcoords, ndim, nve)

!<description>
  ! Auxiliary routine:
  ! Returns the reference coordinates of the corner vertices for a cell type.
!</description>

!<output>
  ! The reference coordinates of the corner vertices.
  real(DP), dimension(:,:), intent(out) :: Dcoords
!</output>

!<input>
  ! The dimension of the element.
  integer, intent(in) :: ndim

  ! The number of corner vertices of the element.
  integer, intent(in) :: nve
!</input>

!</subroutine>

    ! What cell type to we have here?
    select case(ndim)
    case(1)
      ! 1D, so it has to be an edge -> reference coordinates
      Dcoords(1,1) = -1.0_DP
      Dcoords(1,2) =  1.0_DP

    case(2)
      ! 2D
      select case(nve)
      case(3)
        ! Triangle -> barycentric coordinates
        Dcoords(1,1) = 1.0_DP
        Dcoords(2,1) = 0.0_DP
        Dcoords(3,1) = 0.0_DP
        Dcoords(1,2) = 0.0_DP
        Dcoords(2,2) = 1.0_DP
        Dcoords(3,2) = 0.0_DP
        Dcoords(1,3) = 0.0_DP
        Dcoords(2,3) = 0.0_DP
        Dcoords(3,3) = 1.0_DP

      case(4)
        ! Quadrilateral -> reference coordinates
        Dcoords(1,1) = -1.0_DP
        Dcoords(2,1) = -1.0_DP
        Dcoords(1,2) =  1.0_DP
        Dcoords(2,2) = -1.0_DP
        Dcoords(1,3) =  1.0_DP
        Dcoords(2,3) =  1.0_DP
        Dcoords(1,4) = -1.0_DP
        Dcoords(2,4) =  1.0_DP

      end select

    case(3)
      ! 3D
      select case(nve)
      case(4)
        ! Tetrahedron -> barycentric coordinates
        Dcoords(1,1) = 1.0_DP
        Dcoords(2,1) = 0.0_DP
        Dcoords(3,1) = 0.0_DP
        Dcoords(4,1) = 0.0_DP
        Dcoords(1,2) = 0.0_DP
        Dcoords(2,2) = 1.0_DP
        Dcoords(3,2) = 0.0_DP
        Dcoords(4,2) = 0.0_DP
        Dcoords(1,3) = 0.0_DP
        Dcoords(2,3) = 0.0_DP
        Dcoords(3,3) = 1.0_DP
        Dcoords(4,3) = 0.0_DP
        Dcoords(1,4) = 0.0_DP
        Dcoords(2,4) = 0.0_DP
        Dcoords(3,4) = 0.0_DP
        Dcoords(4,4) = 1.0_DP

      case(5)
        ! Pyramid -> reference coordinates
        Dcoords(1,1) = -1.0_DP
        Dcoords(2,1) = -1.0_DP
        Dcoords(3,1) =  0.0_DP
        Dcoords(1,2) =  1.0_DP
        Dcoords(2,2) = -1.0_DP
        Dcoords(3,2) =  0.0_DP
        Dcoords(1,3) =  1.0_DP
        Dcoords(2,3) =  1.0_DP
        Dcoords(3,3) =  0.0_DP
        Dcoords(1,4) = -1.0_DP
        Dcoords(2,4) =  1.0_DP
        Dcoords(3,4) =  0.0_DP
        Dcoords(1,5) =  0.0_DP
        Dcoords(2,5) =  0.0_DP
        Dcoords(3,5) =  1.0_DP

      case(6)
        ! Prism -> reference coordinates
        Dcoords(1,1) =  0.0_DP
        Dcoords(2,1) =  0.0_DP
        Dcoords(3,1) = -1.0_DP
        Dcoords(1,2) =  1.0_DP
        Dcoords(2,2) =  0.0_DP
        Dcoords(3,2) = -1.0_DP
        Dcoords(1,3) =  0.0_DP
        Dcoords(2,3) =  1.0_DP
        Dcoords(3,3) = -1.0_DP
        Dcoords(1,4) =  0.0_DP
        Dcoords(2,4) =  0.0_DP
        Dcoords(3,4) =  1.0_DP
        Dcoords(1,5) =  1.0_DP
        Dcoords(2,5) =  0.0_DP
        Dcoords(3,5) =  1.0_DP
        Dcoords(1,6) =  0.0_DP
        Dcoords(2,6) =  1.0_DP
        Dcoords(3,6) =  1.0_DP

      case(8)
        ! Hexahedron -> reference coordinates
        Dcoords(1,1) = -1.0_DP
        Dcoords(2,1) = -1.0_DP
        Dcoords(3,1) = -1.0_DP
        Dcoords(1,2) =  1.0_DP
        Dcoords(2,2) = -1.0_DP
        Dcoords(3,2) = -1.0_DP
        Dcoords(1,3) =  1.0_DP
        Dcoords(2,3) =  1.0_DP
        Dcoords(3,3) = -1.0_DP
        Dcoords(1,4) = -1.0_DP
        Dcoords(2,4) =  1.0_DP
        Dcoords(3,4) = -1.0_DP
        Dcoords(1,5) = -1.0_DP
        Dcoords(2,5) = -1.0_DP
        Dcoords(3,5) =  1.0_DP
        Dcoords(1,6) =  1.0_DP
        Dcoords(2,6) = -1.0_DP
        Dcoords(3,6) =  1.0_DP
        Dcoords(1,7) =  1.0_DP
        Dcoords(2,7) =  1.0_DP
        Dcoords(3,7) =  1.0_DP
        Dcoords(1,8) = -1.0_DP
        Dcoords(2,8) =  1.0_DP
        Dcoords(3,8) =  1.0_DP

      end select

    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  pure subroutine spdp_aux_getMidpointRefCoords(Dcoords, ndim, nve)

!<description>
  ! Auxiliary routine:
  ! Returns the reference coordinates of the element midpoint for a cell type.
!</description>

!<output>
  ! The reference coordinates of the element midpoint.
  real(DP), dimension(:,:), intent(out) :: Dcoords
!</output>

!<input>
  ! The dimension of the element.
  integer, intent(in) :: ndim

  ! The number of corner vertices of the element.
  integer, intent(in) :: nve
!</input>

!</subroutine>

    ! What cell type to we have here?
    select case(ndim)
    case(1)
      ! 1D, so it has to be an edge -> reference coordinates
      Dcoords = 0.0_DP

    case(2)
      ! 2D
      select case(nve)
      case(3)
        ! Triangle -> barycentric coordinates
        Dcoords = 1.0_DP / 3.0_DP

      case(4)
        ! Quadrilateral -> reference coordinates
        Dcoords = 0.0_DP

      end select

    case(3)
      ! 3D
      select case(nve)
      case(4)
        ! Tetrahedron -> barycentric coordinates
        Dcoords = 0.25_DP

      case(5)
        ! Pyramid -> reference coordinates
        Dcoords(1,1) = 0.0_DP
        Dcoords(2,1) = 0.0_DP
        Dcoords(3,1) = 0.5_DP

      case(6)
        ! Prism -> reference coordinates
        Dcoords(1,1) = 1.0_DP / 3.0_DP
        Dcoords(2,1) = 1.0_DP / 3.0_DP
        Dcoords(3,1) = 0.0_DP

      case(8)
        ! Hexahedron -> reference coordinates
        Dcoords = 0.0_DP

      end select

    end select

  end subroutine

end module
