!##############################################################################
!# ****************************************************************************
!# <name> spdiscprojection </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines to project solution vectors between different
!# finite element spaces/discretisations. This allowes e.g. to 'convert'
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
!#        $P_1$ and/or $Q_1$ elements.
!#
!# 4.) spdp_projectToVertices
!#     -> Project a scalar vector from primal space to the vertices of the
!#        underlying mesh.
!#
!# 5.) spdp_projectToToCells
!#     -> Project a scalar vector from primal space to the cells of the
!#        underlying mesh.
!# </purpose>
!##############################################################################

module spdiscprojection

  use fsystem
  use storage
  use triangulation
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use feevaluation
  use derivatives
  use linearalgebra
  
  implicit none

  private
  
  public :: spdp_projectSolutionScalar
  public :: spdp_projectSolution
  public :: spdp_stdProjectionToP1Q1Scalar
  public :: spdp_projectToVertices
  public :: spdp_projectToCells

contains

  ! ***************************************************************************

!<subroutine>
  subroutine spdp_projectSolutionScalar (rsourceVector,rdestVector)
  
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
  type(t_vectorScalar), intent(IN) :: rsourceVector
!</input>

!<output>
  ! An existing scalar vector structure that receives the projected 
  ! solution vector. Must provide a scalar discretisation structure 
  ! that specifies the destination FE spaces.
  type(t_vectorScalar), intent(INOUT) :: rdestVector
!</output>

!</subroutine>

    ! local variables
    type(t_spatialDiscretisation), pointer :: p_rsourceDiscr,p_rdestDiscr
    type(t_triangulation), pointer :: p_rtriangulation
    real(DP), dimension(:), pointer :: p_Dsource,p_Ddest
    integer(PREC_VERTEXIDX), dimension(:), pointer :: p_IelementsAtVertexIdx 
    integer(PREC_ELEMENTIDX), dimension(:), pointer :: p_IelementsAtVertex
    integer(PREC_VERTEXIDX), dimension(:,:), pointer :: p_IverticesAtElement
    integer(PREC_EDGEIDX), dimension(:,:), pointer :: p_IedgesAtElement,&
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
      print *,'spdp_projectSolutionScalar: No source discretisation!'
      call sys_halt()
    end if

    if (.not. associated(p_rdestDiscr)) then
      print *,'spdp_projectSolutionScalar: No destination discretisation!'
      call sys_halt()
    end if
    
    ! We need to comment this out for 1D/3D projections :(
!    IF (.NOT. ASSOCIATED(p_rsourceDiscr%p_rboundary, p_rdestDiscr%p_rboundary)) THEN
!      PRINT *,'spdp_projectSolutionScalar: Different domains'
!      CALL sys_halt()
!    END IF
    
    if (.not. associated(p_rsourceDiscr%p_rtriangulation,&
                         p_rdestDiscr%p_rtriangulation)) then
      print *,'spdp_projectSolutionScalar: Different triangulations'
      call sys_halt()
    end if
    
    if ((p_rsourceDiscr%ccomplexity .ne. SPDISC_UNIFORM) .or. &
        (p_rdestDiscr%ccomplexity .ne. SPDISC_UNIFORM)) then
      print *,'spdp_projectSolutionScalar: Only uniform discretisations supported!'
      call sys_halt()
    end if
    
    if ((rsourceVector%isortStrategy .gt. 0) .or. &
        (rdestVector%isortStrategy .gt. 0)) then
      print *,'spdp_projectSolutionScalar: Vectors must be unsorted for projection!'
      call sys_halt()
    end if
    
    ! Ok, now we have a chance that we can convert.
    ! If the spaces are identical, we can simply copy the vector
    if (p_rsourceDiscr%RelementDistr(1)%celement .eq. &
        p_rdestDiscr%RelementDistr(1)%celement) then
        
      ! Ok, that's easy.
      ! Copy the vector data but prevent structural data from being overwritten.
      ! Let's hope the vectors have the same length :)
      ! (otherwise the copy-routine will quit)
      call lsyssc_duplicateVector (rsourceVector,rdestVector,&
          LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPYOVERWRITE)
      return
      
    end if
    
    if ((rsourceVector%cdataType .ne. ST_DOUBLE) .or. &
        (rdestVector%cdataType .ne. ST_DOUBLE)) then
      print *,'spdp_projectSolutionScalar: Only double precision vectors supported!'
      call sys_halt()
    end if
    
    ! Clear the destination vector
    call lsyssc_clearVector (rdestVector)
    
    ! What is the destination space?
    ! Remark: Currently only 2D / 3D Q1 supported...
    select case (elem_getPrimaryElement(p_rdestDiscr%RelementDistr(1)%&
                                        celement))
    case (EL_Q1)
      ! So we should convert the source vector into a Q1 destination vector.
      ! Which element is used in the trial space?
      select case (elem_getPrimaryElement(p_rsourceDiscr%RelementDistr(1)%&
                                          celement))
      case (EL_Q0)
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
        
        ! Call the conversion routine
        call spdp_Q0toQ1_dble (p_Dsource, p_Ddest, p_rtriangulation%NVT, &
                               p_IelementsAtVertexIdx,p_IelementsAtVertex)
      case (EL_Q1T,EL_Q1TB,EL_Q2T,EL_Q2TB)
        ! That's a little bit harder. We have to convert an FE space with DOF's
        ! in the midpoints to Q1. (For simplicity, the integral mean value variant
        ! is treated as if the DOF's were in the edge midpoints. The error
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
        call spdp_E030toQ1_dble (p_Dsource, p_Ddest, &
                                 p_rtriangulation%NVT, p_rtriangulation%NEL, &
                                 p_IverticesAtElement,p_IedgesAtElement,&
                                 p_IelementsAtVertexIdx)
                          
      case (EL_Q2)
        ! Rather easy. Take the forst NVT elements of the Q2-vector
        ! as values in the corners of Q1.
        call lsyssc_getbase_double (rsourceVector,p_Dsource)
        call lsyssc_getbase_double (rdestVector,p_Ddest)
        call lalg_copyVectorDble (p_Dsource(1:size(p_Ddest)),p_Ddest)

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
        call spdp_QP1toQ1_dble (p_Dsource, p_Ddest, &
                                 p_rtriangulation%NVT, p_rtriangulation%NEL, &
                                 p_IverticesAtElement,&
                                 p_IelementsAtVertexIdx)
     
      case default
        print *,'spdp_projectSolutionScalar: Unsupported element in source space!'
        call sys_halt()
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
        call spdp_Q0toQ1_dble (p_Dsource, p_Ddest, p_rtriangulation%NVT, &
                                  p_IelementsAtVertexIdx,p_IelementsAtVertex)
      case (EL_Q1T_3D)
        ! That's a little bit harder. We have to convert an FE space with DOF's
        ! in the midpoints to Q1. (For simplicity, the integral mean value variant
        ! is treated as if the DOF's were in the edge midpoints. The error
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
        call spdp_E030toQ1_3D_dble (p_Dsource, p_Ddest, p_rtriangulation%NVT,&
                                    p_rtriangulation%NMT, p_rtriangulation%NEL, &
                                    p_IverticesAtElement,p_IfacesAtElement,&
                                    p_IelementsAtVertexIdx)
                          
        case default
        print *,'spdp_projectSolutionScalar: Unsupported element in source space!'
        call sys_halt()
      end select

    case default
      print *,'spdp_projectSolutionScalar: Unsupported element in destination space!'
      call sys_halt()
    end select
    
  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>
  subroutine spdp_Q0toQ1_dble (Dsource, Ddest, NVT, &
                               IelementsAtVertexIdx,IelementsAtVertex)
  
!<description>
  ! AUXILIARY ROUTINE.
  ! Convert a solution vector based on a uniform discretisation with Q0
  ! to a solution vector based on the Q1 element.
!</description>

!<input>
  ! Source vector to be converted
  real(DP), dimension(:), intent(IN) :: Dsource
  
  ! Number of vertices in the triangulation
  integer(PREC_VERTEXIDX), intent(IN) :: NVT
  
  ! IelementsAtVertexIdx array of the triangulation
  integer(PREC_VERTEXIDX), dimension(:), intent(IN) :: IelementsAtVertexIdx 

  ! IelementsAtVertex array of the triangulation
  integer(PREC_ELEMENTIDX), dimension(:), intent(IN) :: IelementsAtVertex
!</input>
  
!<output>
  ! Destination vector of size NVT; receives the interpolated solution
  real(DP), dimension(:), intent(OUT) :: Ddest
!</output>

!</subroutine>

    ! local variables
    integer(PREC_VERTEXIDX) :: iv
    integer :: nadj
    integer(PREC_ELEMENTIDX) :: ielidx

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
  subroutine spdp_E030toQ1_dble (Dsource, Ddest, NVT, NEL, &
                                 IverticesAtElement,IedgesAtElement,&
                                 IelementsAtVertexIdx)
  
!<description>
  ! AUXILIARY ROUTINE.
  ! Convert a solution vector based on a uniform discretisation with E030, E031,
  ! EM30, EM31 to a solution vector based on the Q1 element.
!</description>

!<input>
  ! Source vector to be converted
  real(DP), dimension(:), intent(IN) :: Dsource
  
  ! Number of vertices in the triangulation
  integer(PREC_VERTEXIDX), intent(IN) :: NVT
  
  ! Number of elements in the triangulation
  integer(PREC_ELEMENTIDX), intent(IN) :: NEL
  
  ! IelementsAtVertexIdx array of the triangulation
  integer(PREC_VERTEXIDX), dimension(:), intent(IN) :: IelementsAtVertexIdx 
  
  ! IverticesAtElement array of the triangulation (old KVERT)
  integer(PREC_VERTEXIDX), dimension(:,:), intent(IN) :: IverticesAtElement

  ! IedgesAtElement array of the triangulation (old KMID)
  integer(PREC_EDGEIDX), dimension(:,:), intent(IN) :: IedgesAtElement
!</input>
  
!<output>
  ! Destination vector of size NVT; receives the interpolated solution
  real(DP), dimension(:), intent(OUT) :: Ddest
!</output>

!</subroutine>

    ! local variables
    integer(PREC_VERTEXIDX) :: iv
    integer(PREC_ELEMENTIDX) :: iel
    integer(PREC_EDGEIDX) :: IM1,IM2,IM3,IM4
    integer(PREC_VERTEXIDX) :: IV1,IV2,IV3,IV4
    real(DP) :: DUH1,DUH2,DUH3,DUH4
    integer :: nadj
    
    ! Clear the output array
    call lalg_clearVectorDble (Ddest)
    
    ! Loop through the elements
    do iel=1,NEL
    
      ! Get the global DOF's on the current element in the E030 space
      IM1 = IedgesAtElement(1,iel)
      IM2 = IedgesAtElement(2,iel)
      IM3 = IedgesAtElement(3,iel)
      IM4 = IedgesAtElement(4,iel)
      
      ! Get the global DOF's on the current element in the Q1 space
      IV1 = IverticesAtElement(1,iel)
      IV2 = IverticesAtElement(2,iel)
      IV3 = IverticesAtElement(3,iel)
      IV4 = IverticesAtElement(4,iel)
      
      ! Get the values of the DOF's in the E030 space
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

  subroutine spdp_E037toQ1_dble (Dsource, Ddest, NVT, NEL, &
                                 IverticesAtElement,IedgesAtElement,&
                                 IelementsAtVertexIdx,ItwistIndexEdges)
  
!<description>
  ! AUXILIARY ROUTINE.
  ! Convert a solution vector based on a uniform discretisation with E030, E031,
  ! EM30, EM31 to a solution vector based on the Q1 element.
!</description>

!<input>
  ! Source vector to be converted
  real(DP), dimension(:), intent(IN) :: Dsource
  
  ! Number of vertices in the triangulation
  integer(PREC_VERTEXIDX), intent(IN) :: NVT
  
  ! Number of elements in the triangulation
  integer(PREC_ELEMENTIDX), intent(IN) :: NEL
  
  ! IelementsAtVertexIdx array of the triangulation
  integer(PREC_VERTEXIDX), dimension(:), intent(IN) :: IelementsAtVertexIdx 
  
  ! IverticesAtElement array of the triangulation (old KVERT)
  integer(PREC_VERTEXIDX), dimension(:,:), intent(IN) :: IverticesAtElement

  ! IedgesAtElement array of the triangulation (old KMID)
  integer(PREC_EDGEIDX), dimension(:,:), intent(IN) :: IedgesAtElement
  
  ! ItwistIndexEdges array of the triangulation
  integer(I32), dimension(:), intent(IN) :: ItwistIndexEdges
!</input>
  
!<output>
  ! Destination vector of size NVT; receives the interpolated solution
  real(DP), dimension(:), intent(OUT) :: Ddest
!</output>

!</subroutine>

    ! local variables
    integer(PREC_VERTEXIDX) :: iv
    integer(PREC_ELEMENTIDX) :: iel
    integer(PREC_EDGEIDX) :: IM1,IM2,IM3,IM4
    integer(PREC_VERTEXIDX) :: IV1,IV2,IV3,IV4
    real(DP) :: DUH1,DUH2,DUH3,DUH4
    integer :: nadj
    
    ! Clear the output array
    call lalg_clearVectorDble (Ddest)
    
    ! Loop through the elements
    do iel=1,NEL
    
      ! Get the global DOF's on the current element in the E030 space
      IM1 = IedgesAtElement(1,iel)
      IM2 = IedgesAtElement(2,iel)
      IM3 = IedgesAtElement(3,iel)
      IM4 = IedgesAtElement(4,iel)
      
      ! Get the global DOF's on the current element in the Q1 space
      IV1 = IverticesAtElement(1,iel)
      IV2 = IverticesAtElement(2,iel)
      IV3 = IverticesAtElement(3,iel)
      IV4 = IverticesAtElement(4,iel)
      
      ! Get the values of the DOF's in the E030 space
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
  subroutine spdp_QP1toQ1_dble (Dsource, Ddest, NVT, NEL, &
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
  real(DP), dimension(:), intent(IN) :: Dsource
  
  ! Number of vertices in the triangulation
  integer(PREC_VERTEXIDX), intent(IN) :: NVT
  
  ! Number of elements in the triangulation
  integer(PREC_ELEMENTIDX), intent(IN) :: NEL
  
  ! IelementsAtVertexIdx array of the triangulation
  integer(PREC_VERTEXIDX), dimension(:), intent(IN) :: IelementsAtVertexIdx 
  
  ! IverticesAtElement array of the triangulation (old KVERT)
  integer(PREC_VERTEXIDX), dimension(:,:), intent(IN) :: IverticesAtElement

!</input>
  
!<output>
  ! Destination vector of size NVT; receives the interpolated solution
  real(DP), dimension(:), intent(OUT) :: Ddest
!</output>

!</subroutine>

    ! local variables
    integer(PREC_VERTEXIDX) :: iv
    integer(PREC_ELEMENTIDX) :: iel
    integer(PREC_VERTEXIDX) :: IV1,IV2,IV3,IV4
    integer :: nadj
    
    ! Clear the output array
    call lalg_clearVectorDble (Ddest)
    
    ! Loop through the elements
    do iel=1,NEL
    
      ! Get the global DOF's on the current element in the Q1 space
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
  subroutine spdp_E030toQ1_3D_dble (Dsource, Ddest, NVT, NMT, NEL, &
                                 IverticesAtElement,IfacesAtElement,&
                                 IelementsAtVertexIdx)
  
!<description>
  ! AUXILIARY ROUTINE.
  ! Convert a solution vector based on a uniform discretisation with E030, E031,
  ! EM30, EM31 to a solution vector based on the Q1 element.
!</description>

!<input>
  ! Source vector to be converted
  real(DP), dimension(:), intent(IN) :: Dsource
  
  ! Number of vertices in the triangulation
  integer(PREC_VERTEXIDX), intent(IN) :: NVT
  
  ! Number of edges in the triangulation
  integer(PREC_EDGEIDX), intent(IN) :: NMT
  
  ! Number of elements in the triangulation
  integer(PREC_ELEMENTIDX), intent(IN) :: NEL
  
  ! IelementsAtVertexIdx array of the triangulation
  integer(PREC_VERTEXIDX), dimension(:), intent(IN) :: IelementsAtVertexIdx 
  
  ! IverticesAtElement array of the triangulation (old KVERT)
  integer(PREC_VERTEXIDX), dimension(:,:), intent(IN) :: IverticesAtElement

  ! IfacesAtElement array of the triangulation
  integer(PREC_EDGEIDX), dimension(:,:), intent(IN) :: IfacesAtElement
!</input>
  
!<output>
  ! Destination vector of size NVT; receives the interpolated solution
  real(DP), dimension(:), intent(OUT) :: Ddest
!</output>

!</subroutine>

    ! local variables
    integer(PREC_VERTEXIDX) :: iv
    integer(PREC_ELEMENTIDX) :: iel
    integer(PREC_EDGEIDX), dimension(6) :: F
    integer(PREC_VERTEXIDX), dimension(8) :: V
    real(DP), dimension(6) :: D
    real(DP),parameter :: R13 = 0.333333333333333_DP
    real(DP),parameter :: R23 = 0.666666666666667_DP
    integer :: nadj
    
    ! Clear the output array
    call lalg_clearVectorDble (Ddest)
    
    ! Loop through the elements
    do iel=1,NEL
    
      ! Get the global DOF's on the current element in the E030 space
      F(1:6) = IfacesAtElement(1:6,iel)
      
      ! Get the global DOF's on the current element in the Q1 space
      V(1:8) = IverticesAtElement(1:8,iel)
      
      ! Get the values of the DOF's in the E030 space
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
  subroutine spdp_projectSolution (rsourceVector,rdestVector)
  
!<description>
  ! This routine 'converts' a given solution vector rsourceVector to
  ! another solution vector rdestVector. The discretisation structures
  ! in the subvertors of rdestVector specifies the new FE spaces, rsourceVector 
  ! should be converted to. The new 'projected' solution is build in 
  ! rdestVector.
!</description>

!<input>
  ! The source vector to be projected.
  type(t_vectorBlock), intent(IN) :: rsourceVector
!</input>

!<output>
  ! An existing vector structure that receives the projected 
  ! solution vector. Must provide a scalar discretisation structures
  ! that specifies the destination FE spaces.
  type(t_vectorBlock), intent(INOUT) :: rdestVector
!</output>

!</subroutine>

    ! local variables
    integer :: i
    
    if (rsourceVector%nblocks .ne. rdestVector%nblocks) then
      print *,'spdp_projectSolution: Different block structure!'
      call sys_halt()
    end if
    
    ! Apply spdp_projectSolutionScalar to every subvector, that's all
    do i=1,rsourceVector%nblocks
      call spdp_projectSolutionScalar (rsourceVector%RvectorBlock(i),&
                                       rdestVector%RvectorBlock(i))
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine spdp_stdProjectionToP1Q1Scalar (rsourceVector,&
      rdestVector,rdestDiscretisation)
  
!<description>
  ! This routine 'converts' a given scalar solution vector rsourceVector to
  ! another solution vector rdestVector. If necessary, the destination
  ! vector is allocated.
  ! The source vector can be an arbitrary FE solution vector. The destination
  ! vector will be a solution vector in the $Q_1$ space (for quad
  ! elements) and in the $P_1$ space (for triangular elements), respectively.
  ! (Such vertex based solution vectors are typically written out
  ! to external files like GMV during the postprocessing).
  !
  ! If rdestDiscretisation is defined, it must describe a discretisation 
  ! with $P_1$ and $Q_1$ elements, respectively. 
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
  type(t_vectorScalar), intent(IN) :: rsourceVector
!</input>

!<inputoutput>
  ! A scalar vector structure that receives the projected 
  ! solution vector. If undefined, a vector is created.
  type(t_vectorScalar), intent(INOUT) :: rdestVector

  ! A discretisation structure that defines a discretisation with
  ! $P_1$ and/or $Q_1$ elements. If undefines, the structure
  ! is automatically created.  
  type(t_spatialDiscretisation), intent(INOUT) :: rdestDiscretisation
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
          EL_P1,EL_Q1,SPDISC_CUB_AUTOMATIC,SPDISC_CUB_AUTOMATIC,&
          rdestDiscretisation)
    
    end if
    
    ! Is the destination vector given?
    if (rdestVector%NEQ .eq. 0) then
      call lsyssc_createVecByDiscr (&
          rdestDiscretisation,rdestVector,.true.)
    end if
    
    ! Project the solution to the $P_1$ / $Q_1$ space.
    call spdp_projectSolutionScalar (rsourceVector,rdestVector)

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine spdp_projectToVertices (rvector, p_Dvalues)
  
!<description>
  ! This routines projects a vector from primal space to the vertices of the
  ! mesh, i.e. the FE function which is described by the coefficient vector
  ! is evaluated in the vertices of its underlying mesh.
  ! The coefficient vector must be unsorted.
!</description>

!<input>
  ! The coefficient vector to be projected.
  type(t_vectorScalar), intent(IN) :: rvector
!</input>

!<inputoutput>
  ! An array that recieves the values of the FE function in the vertices of
  ! the mesh. If set to NULL on entry, an array of the correct size is
  ! allocated.
  real(DP), dimension(:), pointer :: p_Dvalues
!</inputoutput>

!</subroutine>

  ! A hand full of local variables
  type(t_spatialDiscretisation), pointer :: p_rdiscr
  type(t_triangulation), pointer :: p_rtria
  type(t_elementDistribution), pointer :: p_relemDist
  integer, dimension(:), pointer :: p_IelemList, p_IelemAtVertIdx
  integer, dimension(:,:), pointer :: p_IvertAtElem
  real(DP), dimension(4,TRIA_MAXNVE) :: DpointsRef
  real(DP), dimension(TRIA_MAXNVE) :: Deval
  integer :: NVT,NVE,NDIM,ied,ivt,iel,i,j
  
    ! First of all, get the spatial discretisation
    p_rdiscr => rvector%p_rspatialDiscr
    
    if(.not. associated(p_rdiscr)) then
      print *, 'ERROR: spdp_projectToVertices'
      print *, 'Vector does not have a discretisation!'
      call sys_halt()
    end if
    
    ! Now get the triangulation
    p_rtria => p_rdiscr%p_rtriangulation
    
    ! Get the dimension and the number of vertices
    NVT = p_rtria%NVT
    NDIM = p_rtria%ndim
    
    ! Get the vertices-at-element and elements-at-vertex-idx arrays
    call storage_getbase_int2d(p_rtria%h_IverticesAtElement, p_IvertAtElem)
    call storage_getbase_int(p_rtria%h_IelementsAtVertexIdx, p_IelemAtVertIdx)
    
    ! Prepare the values array
    if(.not. associated(p_Dvalues)) then
      allocate(p_Dvalues(NVT))
    else if(ubound(p_Dvalues,1) .lt. NVT) then
      deallocate(p_Dvalues)
      allocate(p_Dvalues(NVT))
    end if
    call lalg_clearVectorDble(p_Dvalues,NVT)
    
    ! Although we call the FE-evaluation routines to evaluate the function,
    ! we still need to loop through all FE spaces, as in a mixed conformal
    ! discretisation, we need to set up the reference coordinates for each
    ! cell type...
    
    ! Okay, now loop through all element distributions
    do ied = 1, p_rdiscr%inumFESpaces
    
      ! Get the element distribution
      p_relemDist => p_rdiscr%RelementDistr(ied)
    
      ! Get the element list for this distribution
      call storage_getbase_int(p_relemDist%h_IelementList, p_IelemList)
      
      ! Get the number of vertice coordinates
      NVE = elem_igetNVE(p_relemDist%celement)
      
      ! Calculate the corner vertice reference coordinates
      call spdp_aux_getCornerRefCoords(DpointsRef, NDIM, NVE)
      
      ! Now go through all elements in this element set
      do i = 1, size(p_IelemList)
      
        ! Get the index of the element
        iel = p_IelemList(i)
        
        ! Call the FE evaluation
        call fevl_evaluate_mult1 (DER_FUNC, Deval, rvector, iel, DpointsRef)
        
        ! Now incorporate the evaluation into the global array
        do j = 1, ubound(p_IvertAtElem,1)
        
          ! Get the index of the global vertice
          ivt = p_IvertAtElem(j,iel)
          
          if(ivt .gt. 0) then
          
            ! Add the local solution onto the global entry
            p_Dvalues(ivt) = p_Dvalues(ivt) + Deval(j)
            
          end if
        
        end do ! j
        
        ! Go for the next element in this distribution
        
      end do ! i
    
      ! Go for the next element distribution
    
    end do ! ied
    
    ! And loop through all vertices
    do ivt = 1, NVT
    
      ! Calculate the number of elements adjacent to this vertice
      j = p_IelemAtVertIdx(ivt+1) - p_IelemAtVertIdx(ivt)
      
      ! And divide the value in this vertice by the number of elements adjacent
      ! to this vertice.
      if(j .gt. 1) then
        p_Dvalues(ivt) = p_Dvalues(ivt) / real(j,DP)
      end if
    
    end do
    
    ! That's it

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine spdp_projectToCells (rvector, p_Dvalues)
  
!<description>
  ! This routines projects a vector from primal space to the cells of the
  ! mesh, i.e. the FE function which is described by the coefficient vector
  ! is evaluated in the element midpoints of its underlying mesh.
  ! The coefficient vector must be unsorted.
!</description>

!<input>
  ! The coefficient vector to be projected.
  type(t_vectorScalar), intent(IN) :: rvector
!</input>

!<inputoutput>
  ! An array that recieves the values of the FE function in the cells of
  ! the mesh. If set to NULL on entry, an array of the correct size is
  ! allocated.
  real(DP), dimension(:), pointer :: p_Dvalues
!</inputoutput>

!</subroutine>

  ! A hand full of local variables
  type(t_spatialDiscretisation), pointer :: p_rdiscr
  type(t_triangulation), pointer :: p_rtria
  type(t_elementDistribution), pointer :: p_relemDist
  integer, dimension(:), pointer :: p_IelemList, p_IelemAtVertIdx
  integer, dimension(:,:), pointer :: p_IvertAtElem
  real(DP), dimension(4,1) :: DpointsRef
  real(DP), dimension(1) :: Deval
  integer :: NEL,NVE,NDIM,ied,iel,i
  
    ! First of all, get the spatial discretisation
    p_rdiscr => rvector%p_rspatialDiscr
    
    if(.not. associated(p_rdiscr)) then
      print *, 'ERROR: spdp_projectToCells'
      print *, 'Vector does not have a discretisation!'
      call sys_halt()
    end if
    
    ! Now get the triangulation
    p_rtria => p_rdiscr%p_rtriangulation
    
    ! Get the dimension and the number of elements
    NEL = p_rtria%NEL
    NDIM = p_rtria%ndim
    
    ! Prepare the values array
    if(.not. associated(p_Dvalues)) then
      allocate(p_Dvalues(NEL))
    else if(ubound(p_Dvalues,1) .lt. NEL) then
      deallocate(p_Dvalues)
      allocate(p_Dvalues(NEL))
    end if
    call lalg_clearVectorDble(p_Dvalues,NEL)
    
    ! Although we call the FE-evaluation routines to evaluate the function,
    ! we still need to loop through all FE spaces, as in a mixed conformal
    ! discretisation, we need to set up the reference coordinates for each
    ! cell type...
    
    ! Okay, now loop through all element distributions
    do ied = 1, p_rdiscr%inumFESpaces
    
      ! Get the element distribution
      p_relemDist => p_rdiscr%RelementDistr(ied)
    
      ! Get the element list for this distribution
      call storage_getbase_int(p_relemDist%h_IelementList, p_IelemList)
      
      ! Get the number of vertice coordinates
      NVE = elem_igetNVE(p_relemDist%celement)
      
      ! Calculate the element midpoint reference coordinates
      call spdp_aux_getMidpointRefCoords(DpointsRef, NDIM, NVE)
      
      ! Now go through all elements in this element set
      do i = 1, size(p_IelemList)
      
        ! Get the index of the element
        iel = p_IelemList(i)
        
        ! Call the FE evaluation
        call fevl_evaluate_mult1 (DER_FUNC, Deval, rvector, iel, DpointsRef)
        
        ! Now incorporate the evaluation into the global array
        p_Dvalues(iel) = Deval(iel)
        
        ! Go for the next element in this distribution
        
      end do ! i
    
      ! Go for the next element distribution
    
    end do ! ied
    
    ! That's it

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
  real(DP), dimension(:,:), intent(OUT) :: Dcoords
!</output>

!<input>
  ! The dimension of the element.
  integer, intent(IN) :: ndim
  
  ! The number of corner vertices of the element.
  integer, intent(IN) :: nve
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
  real(DP), dimension(:,:), intent(OUT) :: Dcoords
!</output>

!<input>
  ! The dimension of the element.
  integer, intent(IN) :: ndim
  
  ! The number of corner vertices of the element.
  integer, intent(IN) :: nve
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
      
      case(8)
        ! Hexahedron -> reference coordinates
        Dcoords = 0.0_DP
      
      end select
    
    end select
  
  end subroutine

end module
