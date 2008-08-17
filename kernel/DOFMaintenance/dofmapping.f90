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

MODULE dofmapping

  USE fsystem
  USE spatialdiscretisation
  USE triangulation
  USE element
  
  IMPLICIT NONE

!<constants>

  !<constantblock description="Kind values for global DOF's">

  ! kind value for indexing global DOF's
  INTEGER, PARAMETER :: PREC_DOFIDX     = I32

  !</constantblock>

!</constants>

CONTAINS

  ! ***************************************************************************

!<function>  

  INTEGER(PREC_DOFIDX) FUNCTION dof_igetNDofGlob(rdiscretisation)

!<description>
  ! This function returns for a given discretisation the number of global
  ! degrees of freedom in the corresponding scalar DOF vector.
!</description>

!<input>    
  ! The discretisation structure that specifies the (scalar) discretisation.
  TYPE(t_spatialDiscretisation), INTENT(IN) :: rdiscretisation
!</input>

!<result>
  ! global number of equations on current grid
!</result>

!</function>

  INTEGER(I32) :: ieltyp
  INTEGER(I32), DIMENSION(2) :: IelTypes

  dof_igetNDofGlob = 0

  SELECT CASE(rdiscretisation%ndimension)
  CASE (NDIM1D)
  
      ieltyp = rdiscretisation%RelementDistr(1)%celement
      dof_igetNDofGlob = NDFG_uniform1D (rdiscretisation%p_rtriangulation, ieltyp)

  CASE (NDIM2D)
    IF (rdiscretisation%ccomplexity .EQ. SPDISC_UNIFORM) THEN

      ieltyp = rdiscretisation%RelementDistr(1)%celement
      ! Uniform discretisation - fall back to the old FEAT mapping
      dof_igetNDofGlob = NDFG_uniform2D (rdiscretisation%p_rtriangulation, ieltyp)

    ELSE IF (rdiscretisation%ccomplexity .EQ. SPDISC_CONFORMAL) THEN

      ! Conformal discretisation. That's a little bit tricky!
      ! At first, we support only the case where two element types are mixed.
      IF (rdiscretisation%inumFESpaces .EQ. 2) THEN
        IelTypes(1) = rdiscretisation%RelementDistr(1)%celement
        IelTypes(2) = rdiscretisation%RelementDistr(2)%celement

        dof_igetNDofGlob = NDFG_conformal2D_2el (&
            rdiscretisation%p_rtriangulation, IelTypes(1:2))
        
      END IF

    END IF
  
  CASE (NDIM3D)
    ! Currently, only uniform discretisations are supported.
    IF (rdiscretisation%ccomplexity .EQ. SPDISC_UNIFORM) THEN

      ieltyp = rdiscretisation%RelementDistr(1)%celement

      ! Uniform discretisation - fall back to the old FEAT mapping
      dof_igetNDofGlob = NDFG_uniform3D (rdiscretisation%p_rtriangulation, ieltyp)
    
    END IF

  CASE DEFAULT
  
    ! Dimension not supported
    PRINT *,'dof_igetNDofGlob: Invalid discretisation dimension!'
    CALL sys_halt()
  
  END SELECT
  
  CONTAINS
  
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Internal subroutine: Get global DOF number for uniform discretisation
    ! with only one element type.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    INTEGER(PREC_DOFIDX) FUNCTION NDFG_uniform1D (rtriangulation, ieltype)
    
    ! IN: The underlying triangulation
    TYPE(t_triangulation), INTENT(IN) :: rtriangulation
    
    ! IN: The element type of the discretisation
    INTEGER(I32), INTENT(IN) :: ieltype
    
    ! OUT: number of global DOF's.
    
    ! The number of global DOF's depends on the element type...
    SELECT CASE (elem_getPrimaryElement(ieltype))
    CASE (EL_P0_1D)
      ! DOF's in the cell midpoints
      NDFG_uniform1D = rtriangulation%NEL
    CASE (EL_P1_1D)
      ! DOF's in the vertices
      NDFG_uniform1D = rtriangulation%NVT
    CASE (EL_P2_1D)
      ! DOF's in the vertices and cell midpoints
      NDFG_uniform1D = rtriangulation%NVT + rtriangulation%NEL
    CASE (EL_S31_1D)
      ! DOF's in the vertices
      NDFG_uniform1D = 2*rtriangulation%NVT
    END SELECT
    
    END FUNCTION

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Internal subroutine: Get global DOF number for uniform discretisation
    ! with only one element type.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! This is roughly the NDFG routine of the old FEAT library...
    
    INTEGER(PREC_DOFIDX) FUNCTION NDFG_uniform2D (rtriangulation, ieltype)
    
    ! IN: The underlying triangulation
    TYPE(t_triangulation), INTENT(IN) :: rtriangulation
    
    ! IN: The element type of the discretisation
    INTEGER(I32), INTENT(IN) :: ieltype
    
    ! OUT: number of global DOF's.
    
    ! The number of global DOF's depends on the element type...
    SELECT CASE (elem_getPrimaryElement(ieltype))
    CASE (EL_P0, EL_Q0)
      ! DOF's in the cell midpoints
      NDFG_uniform2D = rtriangulation%NEL
    CASE (EL_P1, EL_Q1)
      ! DOF's in the vertices
      NDFG_uniform2D = rtriangulation%NVT
    CASE (EL_P2)
      ! DOF's in the vertices and edge midpoints
      NDFG_uniform2D = rtriangulation%NVT + rtriangulation%NMT
    CASE (EL_Q2)
      ! DOF's in the vertices, edge midpoints and element midpoints
      NDFG_uniform2D = rtriangulation%NVT + rtriangulation%NMT + rtriangulation%NEL
    CASE (EL_P3)
      ! 1 DOF's per vertices, 2 DOF per edge 
      NDFG_uniform2D = rtriangulation%NVT + 2*rtriangulation%NMT
    CASE (EL_Q3)
      ! 1 DOF's per vertices, 2 DOF per edge, 4 DOF in the inner
      NDFG_uniform2D = rtriangulation%NVT + 2*rtriangulation%NMT + 4*rtriangulation%NEL
    CASE (EL_QP1)
      ! 3 DOF's in the midpoint of the element.
      NDFG_uniform2D = 3*rtriangulation%NEL
    CASE (EL_P1T, EL_Q1T)
      ! 1 DOF per edge
      NDFG_uniform2D = rtriangulation%NMT
    CASE (EL_Q1TB)
      ! 1 DOF per edge, one per element
      NDFG_uniform2D = rtriangulation%NMT + rtriangulation%NEL
    CASE (EL_Q2T)
      ! DOF's in the vertices, edge midpoints and element midpoints
      NDFG_uniform2D = 2*rtriangulation%NMT + rtriangulation%NEL
    CASE (EL_Q2TB)
      ! E037
      NDFG_uniform2D = 2*rtriangulation%NMT + 2*rtriangulation%NEL
    END SELECT
    
    END FUNCTION

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Internal subroutine: Get global DOF number for conformal discretisation
    ! with two element types.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    INTEGER(PREC_DOFIDX) FUNCTION NDFG_conformal2D_2el (rtriangulation, IelTypes)
    
    ! IN: The underlying triangulation
    TYPE(t_triangulation), INTENT(IN) :: rtriangulation
    
    ! IN: List of element types in the discretisation. IelTypes(1) is one element
    ! type identifier, IelTypes(2) the other one.
    INTEGER(I32), DIMENSION(:), INTENT(IN) :: IelTypes
    
    ! OUT: number of global DOF's.
    
    ! local variables
    INTEGER(I32), DIMENSION(SIZE(IelTypes)) :: IelTypesPrimary
    
    ! Get the primary element number
    IelTypesPrimary = elem_getPrimaryElement(IelTypes)
    
    ! The number of global DOF's depends on the element type...
    SELECT CASE (IelTypesPrimary(1))
    CASE (EL_P0, EL_Q0)
      SELECT CASE (IelTypesPrimary(2))
      CASE (EL_P0, EL_Q0)
        ! DOF's in the cell midpoints
        NDFG_conformal2D_2el = rtriangulation%NEL
      END SELECT
      
    CASE (EL_P1, EL_Q1)
      SELECT CASE (IelTypesPrimary(2))
      CASE (EL_P1, EL_Q1)
        ! DOF's in the vertices
        NDFG_conformal2D_2el = rtriangulation%NVT
      END SELECT
    
    CASE (EL_P2,EL_Q2)
      SELECT CASE (IelTypesPrimary(2))
      CASE (EL_Q2)
        ! Number of vertices + Number of edges (edge midpoints) +
        ! Number of quads (quad midpoints)
        NDFG_conformal2D_2el = rtriangulation%NVT + rtriangulation%NMT + &
            rtriangulation%InelOfType(TRIA_NVEQUAD2D)
      END SELECT

    CASE (EL_P1T,EL_Q1T)
      SELECT CASE (IelTypesPrimary(2))
      CASE (EL_P1T,EL_Q1T)
        ! DOF's in the edge midpoints
        NDFG_conformal2D_2el = rtriangulation%NMT
      END SELECT
      
    END SELECT

    END FUNCTION

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Internal subroutine: Get global DOF number for uniform discretisation
    ! with only one element type.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! This is roughly the NDFG routine of the old FEAT library...
    
    INTEGER(PREC_DOFIDX) FUNCTION NDFG_uniform3D (rtriangulation, ieltype)
    
    ! IN: The underlying triangulation
    TYPE(t_triangulation), INTENT(IN) :: rtriangulation
    
    ! IN: The element type of the discretisation
    INTEGER(I32), INTENT(IN) :: ieltype
    
    ! OUT: number of global DOF's.
    
    ! The number of global DOF's depends on the element type...
    SELECT CASE (elem_getPrimaryElement(ieltype))
    CASE (EL_P0_3D, EL_Q0_3D)
      ! DOF's in the cell midpoints
      NDFG_uniform3D = rtriangulation%NEL
    CASE (EL_P1_3D, EL_Q1_3D)
      ! DOF's in the vertices
      NDFG_uniform3D = rtriangulation%NVT
    CASE (EL_Q1T_3D)
      ! DOF's in the face midpoints
      NDFG_uniform3D = rtriangulation%NAT
    END SELECT
    
    END FUNCTION

  END FUNCTION 

  ! ***************************************************************************

!<function>  

  INTEGER(PREC_DOFIDX) FUNCTION dof_igetNDofGlobBlock(rdiscretisation)

!<description>
  ! This function returns for a given block discretisation the number of global
  ! degrees of freedom in the corresponding block DOF vector.
!</description>

!<input>    
  ! The discretisation structure that specifies the (block) discretisation.
  TYPE(t_blockDiscretisation), INTENT(IN) :: rdiscretisation
!</input>

!<result>
  ! global number of equations on current grid
!</result>

!</function>

    ! Sum up the DOF's of every scalar sub-discretisation structure
    INTEGER :: i
    INTEGER(PREC_DOFIDX) :: icount
    
    icount = 0
    DO i=1,rdiscretisation%ncomponents
      icount = icount + &
          dof_igetNDofGlob(rdiscretisation%RspatialDiscr(i))
    END DO
    
    dof_igetNDofGlobBlock = icount

  END FUNCTION

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE dof_locGlobMapping(rdiscretisation, ielIdx, IdofGlob)
  
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
  TYPE(t_spatialDiscretisation), INTENT(IN) :: rdiscretisation

  ! Element index, where the mapping should be computed.
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: ielIdx

!</input>
    
!<output>

  ! array of global DOF numbers
  INTEGER(PREC_DOFIDX), DIMENSION(:), INTENT(OUT) :: IdofGlob

!</output>

! </subroutine>

  ! local variables
  INTEGER(PREC_ELEMENTIDX), DIMENSION(1) :: ielIdx_array
  INTEGER(PREC_DOFIDX), DIMENSION(SIZE(IdofGlob),1) :: IdofGlob_array

    ! Wrapper to dof_locGlobMapping_mult - better use this directly!
    ielIdx_array(1) = ielidx
    IdofGlob_array(:,1) = IdofGlob(:)
    CALL dof_locGlobMapping_mult (rdiscretisation, ielIdx_array,  &
                                  IdofGlob_array)
    IdofGlob = IdofGlob_array(:,1)

  END SUBROUTINE

  ! ***************************************************************************
!<subroutine>

  SUBROUTINE dof_locGlobMapping_mult(rdiscretisation, IelIdx, IdofGlob)
  
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
  TYPE(t_spatialDiscretisation), INTENT(IN) :: rdiscretisation

  ! Element indices, where the mapping should be computed.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:), INTENT(IN) :: IelIdx
  
!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(OUT) :: IdofGlob

!</output>

! </subroutine>

    ! local variables
    INTEGER(I32), DIMENSION(:,:), POINTER :: p_2darray,p_2darray2
    INTEGER(I32), DIMENSION(:), POINTER :: p_IelementCounter
    TYPE(t_triangulation), POINTER :: p_rtriangulation     
    INTEGER(I32) :: ieltype
    INTEGER(I32), DIMENSION(2) :: IelTypes

    p_rtriangulation => rdiscretisation%p_rtriangulation
    
    SELECT CASE(rdiscretisation%ndimension)
    CASE (NDIM1D)
      ! Call the right 'multiple-get' routines for global DOF's.
      ! For this purpose we evaluate the pointers in the discretisation
      ! structure (if necessary) to prevent another call using pointers...
      ! The number of global DOF's depends on the element type...
      
      ieltype = rdiscretisation%RelementDistr(1)%celement
      
      SELECT CASE (elem_getPrimaryElement(ieltype))
      CASE (EL_P0_1D)
        ! DOF's for P0
        CALL dof_locGlobUniMult_P0_1D(IelIdx, IdofGlob)
        RETURN
      CASE (EL_P1_1D)
        ! DOF's in the vertices
        CALL storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
        CALL dof_locGlobUniMult_P1_1D(p_2darray, IelIdx, IdofGlob)
        RETURN
      CASE (EL_P2_1D)
        ! DOF's in the vertices and line midpoints
        CALL storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
        CALL dof_locGlobUniMult_P2_1D(p_rtriangulation%NVT, p_2darray, IelIdx,&
                                      IdofGlob)
        RETURN
      CASE (EL_S31_1D)
        ! DOF's in the vertices
        CALL storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
        CALL dof_locGlobUniMult_S31_1D(p_rtriangulation%NVT, p_2darray, IelIdx,&
                                       IdofGlob)
        RETURN
      END SELECT
        
    CASE (NDIM2D)
      ! At first we deal only with uniform discretisations
      IF (rdiscretisation%ccomplexity .EQ. SPDISC_UNIFORM) THEN
      
        ! Call the right 'multiple-get' routines for global DOF's.
        ! For this purpose we evaluate the pointers in the discretisation
        ! structure (if necessary) to prevent another call using pointers...
        ! The number of global DOF's depends on the element type...
        
        ieltype = rdiscretisation%RelementDistr(1)%celement
        
        SELECT CASE (elem_getPrimaryElement(ieltype))
        CASE (EL_P0, EL_Q0)
          ! DOF's for Q0
          CALL dof_locGlobUniMult_P0Q0(IelIdx, IdofGlob)
          RETURN
        CASE (EL_P1, EL_Q1)
          ! DOF's in the vertices
          CALL storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
          CALL dof_locGlobUniMult_P1Q1(p_2darray, IelIdx, IdofGlob)
          RETURN
        CASE (EL_P2)
          ! DOF's in the vertices and egde midpoints
          CALL storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
          CALL storage_getbase_int2D (p_rtriangulation%h_IedgesAtElement,p_2darray2)
          CALL dof_locGlobUniMult_P2(p_rtriangulation%NVT,p_2darray, p_2darray2,&
                                     IelIdx, IdofGlob)
          RETURN
        CASE (EL_Q2)
          ! DOF's in the vertices, egde midpoints and element midpoints
          CALL storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
          CALL storage_getbase_int2D (p_rtriangulation%h_IedgesAtElement,p_2darray2)
          CALL dof_locGlobUniMult_Q2(p_rtriangulation%NVT,p_rtriangulation%NMT,&
                                    p_2darray, p_2darray2, IelIdx, IdofGlob)
          RETURN
        CASE (EL_P3) 
        CASE (EL_Q3) 

        CASE (EL_QP1)
          ! DOF's for Q1
          CALL dof_locGlobUniMult_QP1(p_rtriangulation%NEL,IelIdx, IdofGlob)
          RETURN
        CASE (EL_P1T)
          ! DOF's in the edges
          CALL storage_getbase_int2D (p_rtriangulation%h_IedgesAtElement,p_2darray)
          CALL dof_locGlobUniMult_E20(p_2darray, IelIdx, IdofGlob)
          RETURN
        CASE (EL_Q1T)
          ! DOF's in the edges
          CALL storage_getbase_int2D (p_rtriangulation%h_IedgesAtElement,p_2darray)
          CALL dof_locGlobUniMult_E30(p_2darray, IelIdx, IdofGlob)
          RETURN
        CASE (EL_Q1TB)
          ! DOF's in the edges
          CALL storage_getbase_int2D (p_rtriangulation%h_IedgesAtElement,p_2darray)
          CALL dof_locGlobUniMult_E30B(p_rtriangulation%NMT,p_2darray, IelIdx, IdofGlob)
          RETURN
        CASE (EL_Q2T)
          ! DOF's in the edges and the element center
          CALL storage_getbase_int2D (p_rtriangulation%h_IedgesAtElement,p_2darray)
          CALL dof_locGlobUniMult_E035(p_rtriangulation%NMT,p_2darray, IelIdx, IdofGlob)
          RETURN
        CASE (EL_Q2TB)
          ! DOF's in the edges and the element center
          CALL storage_getbase_int2D (p_rtriangulation%h_IedgesAtElement,p_2darray)
          CALL dof_locGlobUniMult_E037(p_rtriangulation%NMT,p_rtriangulation%NEL,&
                                       p_2darray, IelIdx, IdofGlob)
          RETURN
        END SELECT
        
        RETURN
      
      ELSE IF (rdiscretisation%ccomplexity .EQ. SPDISC_CONFORMAL) THEN
      
        ! Conformal discretisation. That's a little bit tricky!
        ! At first, we support only the case where two element types are mixed.
        IF (rdiscretisation%inumFESpaces .EQ. 2) THEN

          IelTypes(1) = rdiscretisation%RelementDistr(1)%celement
          IelTypes(2) = rdiscretisation%RelementDistr(2)%celement

          ! Get the primary element type(s)
          IelTypes = elem_getPrimaryElement(IelTypes)

          ! Now the actual mappings...
          SELECT CASE (IelTypes(1))
          CASE (EL_P0, EL_Q0)
            SELECT CASE (IelTypes(2))
            CASE (EL_P0, EL_Q0)
              ! DOF's in the cell midpoints.
              ! That works like P0 elements.
              CALL dof_locGlobUniMult_P0Q0(IelIdx, IdofGlob)
              RETURN
            END SELECT
            
          CASE (EL_P1, EL_Q1)
            SELECT CASE (IelTypes(2))
            CASE (EL_P1, EL_Q1)
              ! DOF's in the vertices.
              ! That works like P1 elements.
              CALL storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
              CALL dof_locGlobUniMult_P1Q1(p_2darray, IelIdx, IdofGlob)
              RETURN
            END SELECT
            
          CASE (EL_P2, EL_Q2)
            SELECT CASE (IelTypes(2))
            CASE (EL_P2, EL_Q2)
              ! DOF's in the vertices, edges and element mitpoints of the quads.
              ! For this purpose, we need the element counter array that counts
              ! every quad element.
              CALL storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
              CALL storage_getbase_int2D (p_rtriangulation%h_IedgesAtElement,p_2darray2)
              CALL storage_getbase_int (rdiscretisation%h_IelementCounter,p_IelementCounter)
              
              ! Use p_IverticesAtElement evaluated at the first element in the element
              ! set to determine NVE. It's either 3 or 4 and valid for all elements
              ! in the current element set.
              IF (UBOUND(p_2darray,1) .GE. 4) THEN
                IF (p_2darray(4,IelIdx(1)) .EQ. 0) THEN
                  CALL dof_locGlobUniMult_P2Q2(p_rtriangulation%NVT,p_rtriangulation%NMT,3,&
                      p_2darray, p_2darray2, p_IelementCounter, IelIdx, IdofGlob)
                ELSE
                  CALL dof_locGlobUniMult_P2Q2(p_rtriangulation%NVT,p_rtriangulation%NMT,4,&
                      p_2darray, p_2darray2, p_IelementCounter, IelIdx, IdofGlob)
                END IF
              ELSE
                ! Pure triangular mesh
                CALL dof_locGlobUniMult_P2Q2(p_rtriangulation%NVT,p_rtriangulation%NMT,3,&
                    p_2darray, p_2darray2, p_IelementCounter, IelIdx, IdofGlob)
              END IF
              RETURN
            END SELECT

          CASE (EL_P1T, EL_Q1T)
            SELECT CASE (IelTypes(2))
            CASE (EL_P1T, EL_Q1T)
              ! DOF's in the edges
              ! That works like P1 elements.
              CALL storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
              CALL dof_locGlobUniMult_P1Q1(p_2darray, IelIdx, IdofGlob)
              RETURN
            END SELECT
            
          END SELECT

        END IF
      
      END IF
    
    CASE (NDIM3D)
      ! At first we deal only with uniform discretisations
      IF (rdiscretisation%ccomplexity .EQ. SPDISC_UNIFORM) THEN
      
        ! Call the right 'multiple-get' routines for global DOF's.
        ! For this purpose we evaluate the pointers in the discretisation
        ! structure (if necessary) to prevent another call using pointers...
        ! The number of global DOF's depends on the element type...
        ieltype = rdiscretisation%RelementDistr(1)%celement
        
        SELECT CASE (elem_getPrimaryElement(ieltype))
        CASE (EL_P0_3D, EL_Q0_3D)
          ! DOF's for Q0
          CALL dof_locGlobUniMult_P0Q0_3D(IelIdx, IdofGlob)
          RETURN
        CASE (EL_P1_3D, EL_Q1_3D)
          ! DOF's in the vertices
          CALL storage_getbase_int2D (p_rtriangulation%h_IverticesAtElement,p_2darray)
          CALL dof_locGlobUniMult_P1Q1_3D(p_2darray, IelIdx, IdofGlob)
          RETURN
        CASE (EL_Q1T_3D)
          ! DOF's in the face midpoints
          CALL storage_getbase_int2D(p_rtriangulation%h_IfacesAtElement,p_2darray)
          CALL dof_locGlobUniMult_Q1T_3D(p_2darray, IelIdx, IdofGlob)
          RETURN
        END SELECT
        
      END IF
    
    CASE DEFAULT
      PRINT *,'dof_locGlobMapping_mult: invalid discretisation!'
      CALL sys_halt()
    END SELECT

    PRINT *,'dof_locGlobMapping_mult: Unsupported discretisation!'
    CALL sys_halt()

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  PURE SUBROUTINE dof_locGlobUniMult_P0_1D(IelIdx, IdofGlob)
  
!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! All elements in the list are assumed to be P0_1D.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! Element indices, where the mapping should be computed.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:), INTENT(IN) :: IelIdx

!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  INTEGER(I32) :: i
  
  ! Loop through the elements to handle
  DO i=1,SIZE(IelIdx)
    ! Global DOF = number of the element
    IdofGlob(1,i) = IelIdx(i)
  END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  PURE SUBROUTINE dof_locGlobUniMult_P1_1D(IverticesAtElement, IelIdx, IdofGlob)
  
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
  INTEGER(I32), DIMENSION(:,:), INTENT(IN) :: IverticesAtElement

  ! Element indices, where the mapping should be computed.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:), INTENT(IN) :: IelIdx
  
!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  INTEGER(I32) :: i
  
  ! Loop through the elements to handle
  DO i=1,SIZE(IelIdx)
    ! Calculate the global DOF's - which are simply the vertex numbers of the 
    ! corners.
    IdofGlob(1:2,i) = IverticesAtElement(1:2,IelIdx(i))
  END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  PURE SUBROUTINE dof_locGlobUniMult_P2_1D(NVT, IverticesAtElement, IelIdx,&
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
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NVT

  ! An array with the number of vertices adjacent to each element of the
  ! triangulation.
  INTEGER(I32), DIMENSION(:,:), INTENT(IN) :: IverticesAtElement

  ! Element indices, where the mapping should be computed.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:), INTENT(IN) :: IelIdx
  
!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  INTEGER(I32) :: i
  
  ! Loop through the elements to handle
  DO i=1,SIZE(IelIdx)
    ! Calculate the global DOF's - which are simply the vertex numbers of the 
    ! corners and the cell midpoints.
    IdofGlob(1:2,i) = IverticesAtElement(1:2,IelIdx(i))
    IdofGlob(3,i) = NVT + IelIdx(i)
  END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  PURE SUBROUTINE dof_locGlobUniMult_S31_1D(NVT, IverticesAtElement, IelIdx,&
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
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NVT

  ! An array with the number of vertices adjacent to each element of the
  ! triangulation.
  INTEGER(I32), DIMENSION(:,:), INTENT(IN) :: IverticesAtElement

  ! Element indices, where the mapping should be computed.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:), INTENT(IN) :: IelIdx
  
!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  INTEGER(I32) :: i
  
  ! Loop through the elements to handle
  DO i=1,SIZE(IelIdx)
    ! Calculate the global DOF's - which are simply the vertex numbers of the 
    ! corners.
    IdofGlob(1:2,i) = IverticesAtElement(1:2,IelIdx(i))
    IdofGlob(3:4,i) = NVT + IdofGlob(1:2,i)
  END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  PURE SUBROUTINE dof_locGlobUniMult_P0Q0(IelIdx, IdofGlob)
  
!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! all elements in the list are assumed to be P0 or Q0.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! Element indices, where the mapping should be computed.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:), INTENT(IN) :: IelIdx

!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  INTEGER(I32) :: i
  
  ! Loop through the elements to handle
  DO i=1,SIZE(IelIdx)
    ! Global DOF = number of the element
    IdofGlob(1,i) = IelIdx(i)
  END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  PURE SUBROUTINE dof_locGlobUniMult_P1Q1(IverticesAtElement, IelIdx, IdofGlob)
  
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
  INTEGER(I32), DIMENSION(:,:), INTENT(IN) :: IverticesAtElement

  ! Element indices, where the mapping should be computed.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:), INTENT(IN) :: IelIdx
  
!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  INTEGER(I32) :: i,j
  
  ! Get the number of local DOF's - usually either 3 or 4, depending on
  ! the element. The first dimension of IdofGlob indicates the number of 
  ! DOF's.
  j = MIN(UBOUND(IverticesAtElement,1),UBOUND(IdofGlob,1))
  
  ! Loop through the elements to handle
  DO i=1,SIZE(IelIdx)
    ! Calculate the global DOF's - which are simply the vertex numbers of the 
    ! corners.
    IdofGlob(1:j,i) = IverticesAtElement(1:j,IelIdx(i))
  END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  PURE SUBROUTINE dof_locGlobUniMult_P2(NVT, IverticesAtElement, &
                                        IedgesAtElement,IelIdx, IdofGlob)
  
!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! all elements in the list are assumed to be P2.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NVT
  ! An array with the number of vertices adjacent to each element of the
  ! triangulation.
  INTEGER(I32), DIMENSION(:,:), INTENT(IN) :: IverticesAtElement

  ! An array with the number of edges adjacent to each element of the
  ! triangulation.
  INTEGER(I32), DIMENSION(:,:), INTENT(IN) :: IedgesAtElement

  ! Element indices, where the mapping should be computed.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:), INTENT(IN) :: IelIdx
  
!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  INTEGER(I32) :: i
  
  ! Loop through the elements to handle
  DO i=1,SIZE(IelIdx)
    ! Calculate the global DOF's.
    ! The P2 element has global DOF's in the corners and edge midpoints
    ! of the triangles. 
    !
    ! Take the numbers of the corners of the triangles at first.
    IdofGlob(1:3,i) = IverticesAtElement(1:3,IelIdx(i))

    ! Then append the numbers of the edges as midpoint numbers.
    ! Note that the number in this array is NVT+1..NVT+NMT.
    IdofGlob(4:6,i) = NVT + IedgesAtElement(1:3,IelIdx(i))
    
  END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  PURE SUBROUTINE dof_locGlobUniMult_Q2(NVT,NMT,IverticesAtElement, &
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
  INTEGER(I32), DIMENSION(:,:), INTENT(IN) :: IverticesAtElement

  ! An array with the number of edges adjacent to each element of the
  ! triangulation.
  INTEGER(I32), DIMENSION(:,:), INTENT(IN) :: IedgesAtElement

  ! Element indices, where the mapping should be computed.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:), INTENT(IN) :: IelIdx
  
  ! Number of corner vertices in the triangulation
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NVT
  
  ! Number of edes in the triangulation
  INTEGER(PREC_EDGEIDX), INTENT(IN) :: NMT
!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  INTEGER(I32) :: i
  
  ! Loop through the elements to handle
  DO i=1,SIZE(IelIdx)
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
    
  END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  PURE SUBROUTINE dof_locGlobUniMult_P2Q2(NVT,NMT,NVE,IverticesAtElement, &
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
  INTEGER(I32), DIMENSION(:,:), INTENT(IN) :: IverticesAtElement

  ! An array with the number of edges adjacent to each element of the
  ! triangulation.
  INTEGER(I32), DIMENSION(:,:), INTENT(IN) :: IedgesAtElement

  ! Element counter array. This gives every triangle and every quad a
  ! unique running number (1,2,3,...)
  INTEGER(I32), DIMENSION(:), INTENT(IN) :: IelementCounter

  ! Element indices, where the mapping should be computed.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:), INTENT(IN) :: IelIdx
  
  ! Number of corner vertices in the triangulation
  INTEGER(PREC_VERTEXIDX), INTENT(IN) :: NVT
  
  ! Number of edes in the triangulation
  INTEGER(PREC_EDGEIDX), INTENT(IN) :: NMT
  
  ! Element type identifier for which type of elements is currently
  ! under view in IelIdx. All elements in IelIdx are assumed to be of
  ! the same type.
  ! =3: triangular, =4: quad.
  INTEGER, INTENT(IN) :: NVE
!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  INTEGER(I32) :: i
  
  IF (NVE .EQ. 3) THEN
    ! This element set consists of triangular elements.
     
    ! Loop through the elements to handle
    DO i=1,SIZE(IelIdx)
      ! Calculate the global DOF's.
      ! The P2 element has global DOF's in the corners and edge midpoints.
      !
      ! Take the numbers of the corners of the triangles at first.
      IdofGlob(1:3,i) = IverticesAtElement(1:3,IelIdx(i))

      ! Then append the numbers of the edges as midpoint numbers.
      ! Note that the number in this array is NVT+1..NVT+NMT.
      IdofGlob(4:6,i) = NVT+IedgesAtElement(1:3,IelIdx(i))
      
    END DO
     
  ELSE
    ! This element set consists of quad elements.

    ! Loop through the elements to handle
    DO i=1,SIZE(IelIdx)
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
      
    END DO
    
  END IF

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  PURE SUBROUTINE dof_locGlobUniMult_QP1(NEL, IelIdx, IdofGlob)
  
!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! all elements in the list are assumed to be QP1.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! Number of elements in the triangulation
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: NEL

  ! Element indices, where the mapping should be computed.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:), INTENT(IN) :: IelIdx

!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  INTEGER(I32) :: i
  
  ! Loop through the elements to handle
  DO i=1,SIZE(IelIdx)
    ! 1st Global DOF = number of the element = function value
    IdofGlob(1,i) = IelIdx(i)
    ! 2nd Global DOF = NEL + number of the element = X-derivative
    IdofGlob(2,i) = NEL+IelIdx(i)
    ! 3rd Global DOF = 2*NEL + number of the element = Y-derivative
    IdofGlob(3,i) = 2*NEL+IelIdx(i)
  END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  PURE SUBROUTINE dof_locGlobUniMult_E20(IedgesAtElement, IelIdx, IdofGlob)
  
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
  INTEGER(I32), DIMENSION(:,:), INTENT(IN) :: IedgesAtElement

  ! Element indices, where the mapping should be computed.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:), INTENT(IN) :: IelIdx
  
!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  INTEGER(I32) :: i
  
  ! Loop through the elements to handle
  DO i=1,SIZE(IelIdx)
    ! Calculate the global DOF's - which are simply the vertex numbers of the 
    ! corners.
    ! We always copy all elements of IedgesAtElement (:,.).
    ! There's no harm and the compiler can optimise better.
    
    IdofGlob(1:TRIA_NVETRI2D,i) = IedgesAtElement(1:TRIA_NVETRI2D,IelIdx(i))
  END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  PURE SUBROUTINE dof_locGlobUniMult_E30(IedgesAtElement, IelIdx, IdofGlob)
  
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
  INTEGER(I32), DIMENSION(:,:), INTENT(IN) :: IedgesAtElement

  ! Element indices, where the mapping should be computed.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:), INTENT(IN) :: IelIdx
  
!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  INTEGER(I32) :: i
  
  ! Loop through the elements to handle
  DO i=1,SIZE(IelIdx)
    ! Calculate the global DOF's - which are simply the edge numbers.
    ! We always copy all elements of IedgesAtElement (:,.).
    ! There's no harm and the compiler can optimise better.
    
    IdofGlob(1:TRIA_NVEQUAD2D,i) = IedgesAtElement(1:TRIA_NVEQUAD2D,IelIdx(i))
  END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  PURE SUBROUTINE dof_locGlobUniMult_E30B(iNMT, IedgesAtElement, IelIdx, IdofGlob)
  
!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! All elements in the list are assumed to be E030B, E031B, EM30B or EM31B
  ! (with bubble).
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! Number of edges in the triangulation.
  INTEGER(I32), INTENT(IN) :: iNMT

  ! An array with the number of edges adjacent on each element of the
  ! triangulation.
  INTEGER(I32), DIMENSION(:,:), INTENT(IN) :: IedgesAtElement

  ! Element indices, where the mapping should be computed.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:), INTENT(IN) :: IelIdx
  
!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  INTEGER(I32) :: i
  
  ! Loop through the elements to handle
  DO i=1,SIZE(IelIdx)
    ! Calculate the global DOF's - which are simply the numbers of the 
    ! edges. The DOF in the element gets the element number.
    ! We always copy all elements of IedgesAtElement (:,.).
    ! There's no harm and the compiler can optimise better.
    
    IdofGlob(1:TRIA_NVEQUAD2D,i) = IedgesAtElement(1:TRIA_NVEQUAD2D,IelIdx(i))
    IdofGlob(TRIA_NVEQUAD2D+1,i) = iNMT + IelIdx(i)
  END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  PURE SUBROUTINE dof_locGlobUniMult_E035(iNMT, IedgesAtElement, IelIdx, IdofGlob)
  
!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! All elements in the list are assumed to be E035.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! Number of edges in the triangulation.
  INTEGER(I32), INTENT(IN) :: iNMT

  ! An array with the number of edges adjacent on each element of the
  ! triangulation.
  INTEGER(I32), DIMENSION(:,:), INTENT(IN) :: IedgesAtElement

  ! Element indices, where the mapping should be computed.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:), INTENT(IN) :: IelIdx
  
!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  INTEGER(I32) :: i
  
  ! Loop through the elements to handle
  DO i=1,SIZE(IelIdx)
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
  END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  PURE SUBROUTINE dof_locGlobUniMult_E037(iNMT, iNEL, &
      IedgesAtElement, IelIdx, IdofGlob)
  
!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! All elements in the list are assumed to be E037.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! Number of edges in the triangulation.
  INTEGER(I32), INTENT(IN) :: iNMT

  ! Number of elements in the triangulation.
  INTEGER(I32), INTENT(IN) :: iNEL

  ! An array with the number of edges adjacent on each element of the
  ! triangulation.
  INTEGER(I32), DIMENSION(:,:), INTENT(IN) :: IedgesAtElement

  ! Element indices, where the mapping should be computed.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:), INTENT(IN) :: IelIdx
  
!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  INTEGER(I32) :: i
  
  ! Loop through the elements to handle
  DO i=1,SIZE(IelIdx)
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
  END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  PURE SUBROUTINE dof_locGlobUniMult_P0Q0_3D(IelIdx, IdofGlob)
  
!<description>
  ! This subroutine calculates the global indices in the array IdofGlob
  ! of the degrees of freedom of the elements in the list IelIdx.
  ! all elements in the list are assumed to be P0 or Q0.
  ! A uniform grid is assumed, i.e. a grid completely discretised the
  ! same element.
!</description>

!<input>

  ! Element indices, where the mapping should be computed.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:), INTENT(IN) :: IelIdx

!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  INTEGER(I32) :: i
  
  ! Loop through the elements to handle
  DO i=1,SIZE(IelIdx)
    ! Global DOF = number of the element
    IdofGlob(1,i) = IelIdx(i)
  END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  PURE SUBROUTINE dof_locGlobUniMult_P1Q1_3D(IverticesAtElement, IelIdx,&
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
  INTEGER(I32), DIMENSION(:,:), INTENT(IN) :: IverticesAtElement

  ! Element indices, where the mapping should be computed.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:), INTENT(IN) :: IelIdx
  
!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  INTEGER(I32) :: i,j
  
  ! Get the number of local DOF's - usually either 3 or 4, depending on
  ! the element. The first dimension of IdofGlob indicates the number of 
  ! DOF's.
  j = MIN(UBOUND(IverticesAtElement,1),UBOUND(IdofGlob,1))
  
  ! Loop through the elements to handle
  DO i=1,SIZE(IelIdx)
    ! Calculate the global DOF's - which are simply the vertex numbers of the 
    ! corners.
    IdofGlob(1:j,i) = IverticesAtElement(1:j,IelIdx(i))
  END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  PURE SUBROUTINE dof_locGlobUniMult_Q1T_3D(IfacesAtElement, IelIdx, IdofGlob)
  
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
  INTEGER(I32), DIMENSION(:,:), INTENT(IN) :: IfacesAtElement

  ! Element indices, where the mapping should be computed.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:), INTENT(IN) :: IelIdx
  
!</input>
    
!<output>

  ! Array of global DOF numbers; for every element in IelIdx there is
  ! a subarray in this list receiving the corresponding global DOF's.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(OUT) :: IdofGlob

!</output>

!</subroutine>

  ! local variables 
  INTEGER(I32) :: i
  
  ! Loop through the elements to handle
  DO i=1,SIZE(IelIdx)
    ! Calculate the global DOF's - which are simply the face numbers.
    IdofGlob(1:6,i) = IfacesAtElement(1:6,IelIdx(i))
  END DO

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE dof_infoDiscr (rspatialDiscr)
  
!<description>
  ! This routine prints out statistical information about a discretisation
  ! to the terminal.
!</description>

!<inputoutput>
  ! The discretisation structure where information should be printed.
  TYPE(t_spatialDiscretisation), INTENT(IN), TARGET :: rspatialDiscr
!</inputoutput>
  
!</subroutine>

    ! local variables
    INTEGER :: i
    TYPE(t_elementDistribution), POINTER :: p_relementDistr

    ! General information:
    CALL output_line ('Dimension:                    '&
        //TRIM(sys_siL(rspatialDiscr%ndimension,10)))
    CALL output_line ('Complexity:                   ',bnolinebreak=.TRUE.,&
        bnoTrim=.TRUE.)
    SELECT CASE (rspatialDiscr%ccomplexity)
    CASE (SPDISC_UNIFORM)
      CALL output_line ('uniform')
    CASE (SPDISC_CONFORMAL)
      CALL output_line ('conformal')
    CASE (SPDISC_MIXED)
      CALL output_line ('mixed')
    CASE DEFAULT
      CALL output_line ('undefined')
    END SELECT
    CALL output_line ('#DOFs:                        '&
        //TRIM(sys_siL(dof_igetNDofGlob(rspatialDiscr),16)))
    CALL output_line ('#finite element spaces:       '&
        //TRIM(sys_siL(rspatialDiscr%inumFESpaces,10)))
        
    ! Print out detailed information about the FE spaces.
    CALL output_line ('Discretisation details:')
    CALL output_line ('FE-space #elements       NVE   trial-element   test-element')
    
    ! Loop through all element distributions
    DO i=1,rspatialDiscr%inumFESpaces
    
      p_relementDistr => rspatialDiscr%RelementDistr(i)
      
      CALL output_line ( ' ' &
        // sys_siL(i,8) &
        // sys_siL(p_relementDistr%NEL,16) &
        // sys_siL(elem_igetNVE(p_relementDistr%celement),6) &
        // sys_siL(IAND(elem_getPrimaryElement(p_relementDistr%celement),&
                   NOT(EL_DIMENSION)),16) )
      
    END DO
    
  END SUBROUTINE  


  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE dof_infoDiscrBlock (rblockDiscr,bdetailed)
  
!<description>
  ! This routine prints out statistical information about a block 
  ! discretisation to the terminal.
!</description>

!<input>
  ! Whether a detailed description is printed to the terminal or not.
  ! FALSE prints out information only about the block discretisation.
  ! TRUE prints out more detailed information about the block 
  ! discretisation, the structure of the blocks, the used FE spaces etc.
  LOGICAL, INTENT(IN) :: bdetailed

  ! The discretisation structure where information should be printed.
  TYPE(t_blockDiscretisation), INTENT(IN), TARGET :: rblockDiscr
!</input>
  
!</subroutine>

    INTEGER :: i

    CALL output_line ('Dimension:                    '&
        //TRIM(sys_siL(rblockDiscr%ndimension,10)))
    CALL output_line ('Complexity:                   ',bnolinebreak=.TRUE.,&
        bnoTrim=.TRUE.)
    SELECT CASE (rblockDiscr%ccomplexity)
    CASE (SPDISC_UNIFORM)
      CALL output_line ('uniform')
    CASE (SPDISC_CONFORMAL)
      CALL output_line ('conformal')
    CASE (SPDISC_MIXED)
      CALL output_line ('mixed')
    CASE DEFAULT
      CALL output_line ('undefined')
    END SELECT
    CALL output_line ('#DOFs:                        '&
        //TRIM(sys_siL(dof_igetNDofGlobBlock(rblockDiscr),16)))

    CALL output_line ('Number of components:         '&
        //TRIM(sys_siL(rblockDiscr%ncomponents,10)))
        
    IF (bdetailed) THEN
      DO i=1,rblockDiscr%ncomponents
        CALL output_lbrk ()
        CALL output_line ('Solution component:           '//TRIM(sys_siL(i,10)))
        CALL dof_infoDiscr(rblockDiscr%RspatialDiscr(i))
      END DO
    END IF

  END SUBROUTINE

END MODULE
