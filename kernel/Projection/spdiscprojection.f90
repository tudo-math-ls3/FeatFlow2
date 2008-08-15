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
!# </purpose>
!##############################################################################

MODULE spdiscprojection

  USE fsystem
  USE triangulation
  USE spatialdiscretisation
  USE linearsystemscalar
  USE linearsystemblock
  
  IMPLICIT NONE

CONTAINS

  ! ***************************************************************************

!<subroutine>
  SUBROUTINE spdp_projectSolutionScalar (rsourceVector,rdestVector)
  
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
  TYPE(t_vectorScalar), INTENT(IN) :: rsourceVector
!</input>

!<output>
  ! An existing scalar vector structure that receives the projected 
  ! solution vector. Must provide a scalar discretisation structure 
  ! that specifies the destination FE spaces.
  TYPE(t_vectorScalar), INTENT(INOUT) :: rdestVector
!</output>

!</subroutine>

    ! local variables
    TYPE(t_spatialDiscretisation), POINTER :: p_rsourceDiscr,p_rdestDiscr
    TYPE(t_triangulation), POINTER :: p_rtriangulation
    REAL(DP), DIMENSION(:), POINTER :: p_Dsource,p_Ddest
    INTEGER(PREC_VERTEXIDX), DIMENSION(:), POINTER :: p_IelementsAtVertexIdx 
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:), POINTER :: p_IelementsAtVertex
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER :: p_IedgesAtElement,&
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
    
    IF (.NOT. ASSOCIATED(p_rsourceDiscr)) THEN
      PRINT *,'spdp_projectSolutionScalar: No source discretisation!'
      CALL sys_halt()
    END IF

    IF (.NOT. ASSOCIATED(p_rdestDiscr)) THEN
      PRINT *,'spdp_projectSolutionScalar: No destination discretisation!'
      CALL sys_halt()
    END IF
    
    ! We need to comment this out for 1D/3D projections :(
!    IF (.NOT. ASSOCIATED(p_rsourceDiscr%p_rboundary, p_rdestDiscr%p_rboundary)) THEN
!      PRINT *,'spdp_projectSolutionScalar: Different domains'
!      CALL sys_halt()
!    END IF
    
    IF (.NOT. ASSOCIATED(p_rsourceDiscr%p_rtriangulation,&
                         p_rdestDiscr%p_rtriangulation)) THEN
      PRINT *,'spdp_projectSolutionScalar: Different triangulations'
      CALL sys_halt()
    END IF
    
    IF ((p_rsourceDiscr%ccomplexity .NE. SPDISC_UNIFORM) .OR. &
        (p_rdestDiscr%ccomplexity .NE. SPDISC_UNIFORM)) THEN
      PRINT *,'spdp_projectSolutionScalar: Only uniform discretisations supported!'
      CALL sys_halt()
    END IF
    
    IF ((rsourceVector%isortStrategy .GT. 0) .OR. &
        (rdestVector%isortStrategy .GT. 0)) THEN
      PRINT *,'spdp_projectSolutionScalar: Vectors must be unsorted for projection!'
      CALL sys_halt()
    END IF
    
    ! Ok, now we have a chance that we can convert.
    ! If the spaces are identical, we can simply copy the vector
    IF (p_rsourceDiscr%RelementDistr(1)%celement .EQ. &
        p_rdestDiscr%RelementDistr(1)%celement) THEN
        
      ! Ok, that's easy.
      ! Copy the vector data but prevent structural data from being overwritten.
      ! Let's hope the vectors have the same length :)
      ! (otherwise the copy-routine will quit)
      CALL lsyssc_duplicateVector (rsourceVector,rdestVector,&
          LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPYOVERWRITE)
      RETURN
      
    END IF
    
    IF ((rsourceVector%cdataType .NE. ST_DOUBLE) .OR. &
        (rdestVector%cdataType .NE. ST_DOUBLE)) THEN
      PRINT *,'spdp_projectSolutionScalar: Only double precision vectors supported!'
      CALL sys_halt()
    END IF
    
    ! Clear the destination vector
    CALL lsyssc_clearVector (rdestVector)
    
    ! What is the destination space?
    ! Remark: Currently only 2D / 3D Q1 supported...
    SELECT CASE (elem_getPrimaryElement(p_rdestDiscr%RelementDistr(1)%&
                                        celement))
    CASE (EL_Q1)
      ! So we should convert the source vector into a Q1 destination vector.
      ! Which element is used in the trial space?
      SELECT CASE (elem_getPrimaryElement(p_rsourceDiscr%RelementDistr(1)%&
                                          celement))
      CASE (EL_Q0)
        ! Not too hard. Basically, take the mean of all elements adjacent to a vertex.
        !
        !
        ! Get geometric information from the triangulation.
        p_rtriangulation => p_rsourceDiscr%p_rtriangulation
        CALL storage_getbase_int (p_rtriangulation%h_IelementsAtVertexIdx,&
                                                   p_IelementsAtVertexIdx)
        CALL storage_getbase_int (p_rtriangulation%h_IelementsAtVertex,&
                                                   p_IelementsAtVertex)
                                                     
        ! Get the vector data
        CALL lsyssc_getbase_double (rsourceVector,p_Dsource)
        CALL lsyssc_getbase_double (rdestVector,p_Ddest)
        
        ! Call the conversion routine
        CALL spdp_Q0toQ1_dble (p_Dsource, p_Ddest, p_rtriangulation%NVT, &
                               p_IelementsAtVertexIdx,p_IelementsAtVertex)
      CASE (EL_Q1T,EL_Q1TB,EL_Q2T,EL_Q2TB)
        ! That's a little bit harder. We have to convert an FE space with DOF's
        ! in the midpoints to Q1. (For simplicity, the integral mean value variant
        ! is treated as if the DOF's were in the edge midpoints. The error
        ! is negligible.
        !
        ! Get geometric information from the triangulation.
        p_rtriangulation => p_rsourceDiscr%p_rtriangulation
        CALL storage_getbase_int (p_rtriangulation%h_IelementsAtVertexIdx,&
                                                   p_IelementsAtVertexIdx)
        CALL storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
                                                     p_IverticesAtElement)
        CALL storage_getbase_int2d (p_rtriangulation%h_IedgesAtElement,&
                                                     p_IedgesAtElement)
                                                     
        ! Get the vector data
        CALL lsyssc_getbase_double (rsourceVector,p_Dsource)
        CALL lsyssc_getbase_double (rdestVector,p_Ddest)
        
        ! Call the conversion routine
        CALL spdp_E030toQ1_dble (p_Dsource, p_Ddest, &
                                 p_rtriangulation%NVT, p_rtriangulation%NEL, &
                                 p_IverticesAtElement,p_IedgesAtElement,&
                                 p_IelementsAtVertexIdx)
                          
      CASE (EL_Q2)
        ! Rather easy. Take the forst NVT elements of the Q2-vector
        ! as values in the corners of Q1.
        CALL lsyssc_getbase_double (rsourceVector,p_Dsource)
        CALL lsyssc_getbase_double (rdestVector,p_Ddest)
        CALL lalg_copyVectorDble (p_Dsource(1:SIZE(p_Ddest)),p_Ddest)

      CASE (EL_QP1)
        ! Also not completely trivial. Interpolation of the values on the element
        ! midpoints to the corners, neglecting the error.
        ! Get geometric information from the triangulation.
        p_rtriangulation => p_rsourceDiscr%p_rtriangulation
        CALL storage_getbase_int (p_rtriangulation%h_IelementsAtVertexIdx,&
                                                   p_IelementsAtVertexIdx)
        CALL storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
                                                     p_IverticesAtElement)
                                                     
        ! Get the vector data
        CALL lsyssc_getbase_double (rsourceVector,p_Dsource)
        CALL lsyssc_getbase_double (rdestVector,p_Ddest)
        
        ! Call the conversion routine
        CALL spdp_QP1toQ1_dble (p_Dsource, p_Ddest, &
                                 p_rtriangulation%NVT, p_rtriangulation%NEL, &
                                 p_IverticesAtElement,&
                                 p_IelementsAtVertexIdx)
     
      CASE DEFAULT
        PRINT *,'spdp_projectSolutionScalar: Unsupported element in source space!'
        CALL sys_halt()
      END SELECT
    
    CASE (EL_Q1_3D)
        ! So we should convert the source vector into a 3D Q1 destination vector.
      ! Which element is used in the trial space?
      SELECT CASE (elem_getPrimaryElement(p_rsourceDiscr%RelementDistr(1)%&
                                          celement))
      CASE (EL_Q0_3D)
        ! Not too hard. Basically, take the mean of all elements adjacent to a vertex.
        !
        !
        ! Get geometric information from the triangulation.
        p_rtriangulation => p_rsourceDiscr%p_rtriangulation
        CALL storage_getbase_int (p_rtriangulation%h_IelementsAtVertexIdx,&
                                                   p_IelementsAtVertexIdx)
        CALL storage_getbase_int (p_rtriangulation%h_IelementsAtVertex,&
                                                   p_IelementsAtVertex)
                                                     
        ! Get the vector data
        CALL lsyssc_getbase_double (rsourceVector,p_Dsource)
        CALL lsyssc_getbase_double (rdestVector,p_Ddest)
        
        ! Call the conversion routine - we can use the 2D version here...
        CALL spdp_Q0toQ1_dble (p_Dsource, p_Ddest, p_rtriangulation%NVT, &
                                  p_IelementsAtVertexIdx,p_IelementsAtVertex)
      CASE (EL_Q1T_3D)
        ! That's a little bit harder. We have to convert an FE space with DOF's
        ! in the midpoints to Q1. (For simplicity, the integral mean value variant
        ! is treated as if the DOF's were in the edge midpoints. The error
        ! is negligible.
        !
        ! Get geometric information from the triangulation.
        p_rtriangulation => p_rsourceDiscr%p_rtriangulation
        CALL storage_getbase_int (p_rtriangulation%h_IelementsAtVertexIdx,&
                                                   p_IelementsAtVertexIdx)
        CALL storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
                                                     p_IverticesAtElement)
        CALL storage_getbase_int2d (p_rtriangulation%h_IfacesAtElement,&
                                                     p_IfacesAtElement)
                                                     
        ! Get the vector data
        CALL lsyssc_getbase_double (rsourceVector,p_Dsource)
        CALL lsyssc_getbase_double (rdestVector,p_Ddest)
        
        ! Call the conversion routine
        CALL spdp_E030toQ1_3D_dble (p_Dsource, p_Ddest, p_rtriangulation%NVT,&
                                    p_rtriangulation%NMT, p_rtriangulation%NEL, &
                                    p_IverticesAtElement,p_IfacesAtElement,&
                                    p_IelementsAtVertexIdx)
                          
        CASE DEFAULT
        PRINT *,'spdp_projectSolutionScalar: Unsupported element in source space!'
        CALL sys_halt()
      END SELECT

    CASE DEFAULT
      PRINT *,'spdp_projectSolutionScalar: Unsupported element in destination space!'
      CALL sys_halt()
    END SELECT
    
  END SUBROUTINE
  
  ! ***************************************************************************
  
!<subroutine>
  SUBROUTINE spdp_Q0toQ1_dble (Dsource, Ddest, NVT, &
                               IelementsAtVertexIdx,IelementsAtVertex)
  
!<description>
  ! AUXILIARY ROUTINE.
  ! Convert a solution vector based on a uniform discretisation with Q0
  ! to a solution vector based on the Q1 element.
!</description>

!<input>
  ! Source vector to be converted
  REAL(DP), DIMENSION(:), INTENT(IN) :: Dsource
  
  ! Number of vertices in the triangulation
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NVT
  
  ! IelementsAtVertexIdx array of the triangulation
  INTEGER(PREC_VERTEXIDX), DIMENSION(:), INTENT(IN) :: IelementsAtVertexIdx 

  ! IelementsAtVertex array of the triangulation
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:), INTENT(IN) :: IelementsAtVertex
!</input>
  
!<output>
  ! Destination vector of size NVT; receives the interpolated solution
  REAL(DP), DIMENSION(:), INTENT(OUT) :: Ddest
!</output>

!</subroutine>

    ! local variables
    INTEGER(PREC_VERTEXIDX) :: iv
    INTEGER :: nadj
    INTEGER(PREC_ELEMENTIDX) :: ielidx

    ! Loop through the vertices   
    DO iv=1,NVT
    
      ! On each vertex, loop through the adjacent elements
      DO ielidx = IelementsAtVertexIdx(iv),IelementsAtVertexIdx(iv+1)-1
        ! Sum up the element contributions into the vertex
        Ddest(iv) = Ddest(iv) + Dsource(IelementsAtVertex(ielidx))
      END DO

      ! Divide by the number of adjacent elements, this results
      ! in the interpolated solution.
      nadj = IelementsAtVertexIdx(iv+1) - IelementsAtVertexIdx(iv)
      Ddest(iv) = Ddest(iv) / REAL(nadj,DP)
    END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>
  SUBROUTINE spdp_E030toQ1_dble (Dsource, Ddest, NVT, NEL, &
                                 IverticesAtElement,IedgesAtElement,&
                                 IelementsAtVertexIdx)
  
!<description>
  ! AUXILIARY ROUTINE.
  ! Convert a solution vector based on a uniform discretisation with E030, E031,
  ! EM30, EM31 to a solution vector based on the Q1 element.
!</description>

!<input>
  ! Source vector to be converted
  REAL(DP), DIMENSION(:), INTENT(IN) :: Dsource
  
  ! Number of vertices in the triangulation
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NVT
  
  ! Number of elements in the triangulation
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NEL
  
  ! IelementsAtVertexIdx array of the triangulation
  INTEGER(PREC_VERTEXIDX), DIMENSION(:), INTENT(IN) :: IelementsAtVertexIdx 
  
  ! IverticesAtElement array of the triangulation (old KVERT)
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElement

  ! IedgesAtElement array of the triangulation (old KMID)
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), INTENT(IN) :: IedgesAtElement
!</input>
  
!<output>
  ! Destination vector of size NVT; receives the interpolated solution
  REAL(DP), DIMENSION(:), INTENT(OUT) :: Ddest
!</output>

!</subroutine>

    ! local variables
    INTEGER(PREC_VERTEXIDX) :: iv
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER(PREC_EDGEIDX) :: IM1,IM2,IM3,IM4
    INTEGER(PREC_VERTEXIDX) :: IV1,IV2,IV3,IV4
    REAL(DP) :: DUH1,DUH2,DUH3,DUH4
    INTEGER :: nadj
    
    ! Clear the output array
    CALL lalg_clearVectorDble (Ddest)
    
    ! Loop through the elements
    DO iel=1,NEL
    
      ! Get the global DOF's on the current element in the E030 space
      IM1 = IedgesAtElement(1,iel)-NVT
      IM2 = IedgesAtElement(2,iel)-NVT
      IM3 = IedgesAtElement(3,iel)-NVT
      IM4 = IedgesAtElement(4,iel)-NVT
      
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
      
    END DO

    ! Divide by the number of adjacent elements, this results
    ! in the interpolated solution.
    DO iv=1,NVT
      nadj = IelementsAtVertexIdx(iv+1) - IelementsAtVertexIdx(iv)
      Ddest(iv) = Ddest(iv) / REAL(nadj,DP)
    END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>
  SUBROUTINE spdp_QP1toQ1_dble (Dsource, Ddest, NVT, NEL, &
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
  REAL(DP), DIMENSION(:), INTENT(IN) :: Dsource
  
  ! Number of vertices in the triangulation
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NVT
  
  ! Number of elements in the triangulation
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NEL
  
  ! IelementsAtVertexIdx array of the triangulation
  INTEGER(PREC_VERTEXIDX), DIMENSION(:), INTENT(IN) :: IelementsAtVertexIdx 
  
  ! IverticesAtElement array of the triangulation (old KVERT)
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElement

!</input>
  
!<output>
  ! Destination vector of size NVT; receives the interpolated solution
  REAL(DP), DIMENSION(:), INTENT(OUT) :: Ddest
!</output>

!</subroutine>

    ! local variables
    INTEGER(PREC_VERTEXIDX) :: iv
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER(PREC_VERTEXIDX) :: IV1,IV2,IV3,IV4
    INTEGER :: nadj
    
    ! Clear the output array
    CALL lalg_clearVectorDble (Ddest)
    
    ! Loop through the elements
    DO iel=1,NEL
    
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
      
    END DO

    ! Divide by the number of adjacent elements, this results
    ! in the interpolated solution.
    DO iv=1,NVT
      nadj = IelementsAtVertexIdx(iv+1) - IelementsAtVertexIdx(iv)
      Ddest(iv) = Ddest(iv) / REAL(nadj,DP)
    END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>
  SUBROUTINE spdp_E030toQ1_3D_dble (Dsource, Ddest, NVT, NMT, NEL, &
                                 IverticesAtElement,IfacesAtElement,&
                                 IelementsAtVertexIdx)
  
!<description>
  ! AUXILIARY ROUTINE.
  ! Convert a solution vector based on a uniform discretisation with E030, E031,
  ! EM30, EM31 to a solution vector based on the Q1 element.
!</description>

!<input>
  ! Source vector to be converted
  REAL(DP), DIMENSION(:), INTENT(IN) :: Dsource
  
  ! Number of vertices in the triangulation
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NVT
  
  ! Number of edges in the triangulation
  INTEGER(PREC_EDGEIDX), INTENT(IN) :: NMT
  
  ! Number of elements in the triangulation
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NEL
  
  ! IelementsAtVertexIdx array of the triangulation
  INTEGER(PREC_VERTEXIDX), DIMENSION(:), INTENT(IN) :: IelementsAtVertexIdx 
  
  ! IverticesAtElement array of the triangulation (old KVERT)
  INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), INTENT(IN) :: IverticesAtElement

  ! IfacesAtElement array of the triangulation
  INTEGER(PREC_EDGEIDX), DIMENSION(:,:), INTENT(IN) :: IfacesAtElement
!</input>
  
!<output>
  ! Destination vector of size NVT; receives the interpolated solution
  REAL(DP), DIMENSION(:), INTENT(OUT) :: Ddest
!</output>

!</subroutine>

    ! local variables
    INTEGER(PREC_VERTEXIDX) :: iv
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER(PREC_EDGEIDX), DIMENSION(6) :: F
    INTEGER(PREC_VERTEXIDX), DIMENSION(8) :: V
    REAL(DP), DIMENSION(6) :: D
    REAL(DP),PARAMETER :: R13 = 0.333333333333333_DP
    REAL(DP),PARAMETER :: R23 = 0.666666666666667_DP
    INTEGER :: nadj
    
    ! Clear the output array
    CALL lalg_clearVectorDble (Ddest)
    
    ! Loop through the elements
    DO iel=1,NEL
    
      ! Get the global DOF's on the current element in the E030 space
      F(1:6) = IfacesAtElement(1:6,iel)-NVT-NMT
      
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
      
    END DO

    ! Divide by the number of adjacent elements, this results
    ! in the interpolated solution.
    DO iv=1,NVT
      nadj = IelementsAtVertexIdx(iv+1) - IelementsAtVertexIdx(iv)
      Ddest(iv) = Ddest(iv) / REAL(nadj,DP)
    END DO

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  SUBROUTINE spdp_projectSolution (rsourceVector,rdestVector)
  
!<description>
  ! This routine 'converts' a given solution vector rsourceVector to
  ! another solution vector rdestVector. The discretisation structures
  ! in the subvertors of rdestVector specifies the new FE spaces, rsourceVector 
  ! should be converted to. The new 'projected' solution is build in 
  ! rdestVector.
!</description>

!<input>
  ! The source vector to be projected.
  TYPE(t_vectorBlock), INTENT(IN) :: rsourceVector
!</input>

!<output>
  ! An existing vector structure that receives the projected 
  ! solution vector. Must provide a scalar discretisation structures
  ! that specifies the destination FE spaces.
  TYPE(t_vectorBlock), INTENT(INOUT) :: rdestVector
!</output>

!</subroutine>

    ! local variables
    INTEGER :: i
    
    IF (rsourceVector%nblocks .NE. rdestVector%nblocks) THEN
      PRINT *,'spdp_projectSolution: Different block structure!'
      CALL sys_halt()
    END IF
    
    ! Apply spdp_projectSolutionScalar to every subvector, that's all
    DO i=1,rsourceVector%nblocks
      CALL spdp_projectSolutionScalar (rsourceVector%RvectorBlock(i),&
                                       rdestVector%RvectorBlock(i))
    END DO

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>
  SUBROUTINE spdp_stdProjectionToP1Q1Scalar (rsourceVector,&
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
  TYPE(t_vectorScalar), INTENT(IN) :: rsourceVector
!</input>

!<inputoutput>
  ! A scalar vector structure that receives the projected 
  ! solution vector. If undefined, a vector is created.
  TYPE(t_vectorScalar), INTENT(INOUT) :: rdestVector

  ! A discretisation structure that defines a discretisation with
  ! $P_1$ and/or $Q_1$ elements. If undefines, the structure
  ! is automatically created.  
  TYPE(t_spatialDiscretisation), INTENT(INOUT) :: rdestDiscretisation
!</inputoutput>

!</subroutine>

    ! Check that the source discretisation structure is valid.
    IF (rsourceVector%p_rspatialDiscr%ndimension .NE. NDIM2D) THEN
      CALL output_line ('Only 2D discretisation supported!', &
                        OU_CLASS_ERROR,OU_MODE_STD,&
                        'spdp_stdProjectionToP1Q1')  
      CALL sys_halt()
    END IF

    ! Is the destination discretisation given?
    IF (rdestDiscretisation%ndimension .EQ. 0) THEN
    
      ! Create a discretisation with $P_1$ and/or $Q_1$ elements.
      ! Derive it from the existing discretisation structure in
      ! the vector, so the element information does not use too
      ! much memory...
      
      CALL spdiscr_deriveDiscr_triquad(&
          rsourceVector%p_rspatialDiscr,&
          EL_P1,EL_Q1,SPDISC_CUB_AUTOMATIC,SPDISC_CUB_AUTOMATIC,&
          rdestDiscretisation)
    
    END IF
    
    ! Is the destination vector given?
    IF (rdestVector%NEQ .EQ. 0) THEN
      CALL lsyssc_createVecByDiscr (&
          rdestDiscretisation,rdestVector,.TRUE.)
    END IF
    
    ! Project the solution to the $P_1$ / $Q_1$ space.
    CALL spdp_projectSolutionScalar (rsourceVector,rdestVector)

  END SUBROUTINE

END MODULE
