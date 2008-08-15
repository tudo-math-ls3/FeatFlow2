!##############################################################################
!# ****************************************************************************
!# <name> matrixrestriction </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# This module contains routines to realise a matrix restriction, i.e. the
!# construction of coarse grid matrix operators from fine grid matrices by the
!# Galerkin approach (c.f. p. 224ff in Turek's book).
!#
!# In the Galerkin approach a coarse grid matrix $A_{2h}$ 
!# is build from a fine grid matrix $A_h$ with the help of the prolongation
!# operator $P_{2h}^h$ and the restriction operator $R_{h}^{2h}$ using
!#
!#       $$ A_{2h}  =  P_{2h}^h  A_h  R_{h}^{2h}. $$
!#
!# This is equivalent to the standard approach $a_{ij} = a_h(\phi_j,\phi_i)$
!# in the case of conformal finite elements, but yields a different matrix
!# in the case of nonconformal FE spaces.
!#
!# The routines in this file allow for some situations (discretisation
!# with $\tilde Q_1$ e.g.) to rebuild some lines in a given matrix according
!# to this approach. This is a kind of stabilisation in case the mesh elements
!# show high aspect ratios.
!#
!# The following subroutines can be found here:
!#
!# 1.) mrest_matrixRestrictionEX3Y
!#     -> Build coarse grid matrix entries on elements with too high
!#        aspect ratio by using a fine grid matrix.
!#
!# </purpose>
!##############################################################################

MODULE matrixrestriction

  USE fsystem
  USE linearsystemscalar
  USE matrixmodification
  USE geometryaux
  
  IMPLICIT NONE

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE mrest_matrixRestrictionEX3Y (rfineMatrix,rcoarseMatrix,&
                                          iARindicator, dARbound)
  
!<description>
  ! This routine rebuilds some line in the scalar coarse grid matrix 
  ! rcoarseMatrix with the help of the Galerkin approach using the scalar 
  ! fine grid matrix rfineMatrix. Only the lines of those DOF's will be 
  ! rebuild where the elements that belong to the DOF show too high aspect 
  ! ratios.
  !
  ! The routine works only for matrices that are build using a uniform
  ! discretisation with $\tilde Q_1$.
!</description>

!<input>
  ! Fine grid matrix.
  TYPE(t_matrixScalar), INTENT(IN) :: rfineMatrix
  
  ! Aspect-ratio indicator. Configures when to rebuild a line of the
  ! coarse matrix by constant restriction.
  ! <=0: do nothing, don't modify coarse grid matrix
  !  =1: switch depending on aspect ratio of current element
  !  =2: switch depending on aspect ratio of current element
  !      and aspect ratio of neighbour element
  INTEGER, INTENT(IN) :: iARindicator
  
  ! Maximum allowed aspect ratio. Rows in the matrix corresponding
  ! to elements (i.e. a DOF's of an element) with an aspect ratio
  ! larger than dARbound are rebuild with the Galerkin approach
  ! by constant restriction.
  REAL(DP), INTENT(IN) :: dARbound
  
!</input>

!<inputoutput>
  ! Coarse grid matrix to be modified.
  TYPE(t_matrixScalar), INTENT(IN) :: rcoarseMatrix
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER(I32) :: i1,i2
    INTEGER :: iedge,iedge1,iedge2,iedge4,jedge,jedge1,jedge2,jedge4
    INTEGER(PREC_ELEMENTIDX) :: iel,iadj1,iadj2,iadj4,iel1,iel2,jadj2,jadj4,jel1,jel2
    INTEGER(PREC_EDGEIDX) :: imid1,imid2,imid3,imid4,imid5
    INTEGER(PREC_EDGEIDX) :: im1,im2,im3,im4,im5,im6,im7,im8,im9,im10,im11,im12
    INTEGER(PREC_VERTEXIDX) :: nvt1,nvt2
    INTEGER(PREC_MATIDX) :: ild
    TYPE(t_triangulation), POINTER :: p_rtriaCoarse,p_rtriaFine
    REAL(DP), DIMENSION(NDIM2D,TRIA_MAXNVE2D) :: dcoords
    REAL(DP), DIMENSION(0:TRIA_MAXNME2D) :: daspectRatio
    INTEGER(PREC_ELEMENTIDX), DIMENSION(TRIA_MAXNME2D) :: Iiel
    INTEGER(PREC_ELEMENTIDX), DIMENSION(0:TRIA_MAXNME2D) :: ielAdjacent
    REAL(DP) :: dval1,dval2,dval3,dval4,dval5
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IneighboursAtElementCoarse
    INTEGER(PREC_ELEMENTIDX), DIMENSION(:,:), POINTER :: p_IneighboursAtElementFine
    REAL(DP), DIMENSION(:,:), POINTER                 :: p_DvertexCoordsCoarse
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER    :: p_IedgesAtElementCoarse
    INTEGER(PREC_EDGEIDX), DIMENSION(:,:), POINTER    :: p_IedgesAtElementFine
    INTEGER(PREC_VERTEXIDX), DIMENSION(:,:), POINTER   :: p_IverticesAtElementCoarse
    
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER       :: p_KldCoarse,p_KldFine
    INTEGER(PREC_MATIDX), DIMENSION(:), POINTER       :: p_KdiagonalCoarse,p_KdiagonalFine
    INTEGER(PREC_VECIDX), DIMENSION(:), POINTER       :: p_KcolCoarse,p_KcolFine
    REAL(DP), DIMENSION(:), POINTER                   :: p_DaCoarse,p_DaFine
    REAL(DP) :: dv1,dv2,dv3,dv4,dv5,dv6,dv7,dv8,dv9,dv10,dv11,dv12
    
    ! No modification if the parameters are out of bounds!
    IF (iARindicator .LE. 0) RETURN
    IF (dARbound .LT. 0.0_DP) RETURN
    
    ! At first some basic checks if we are able to complete our task at all.
    
    IF ((rfineMatrix%cmatrixFormat .NE. LSYSSC_MATRIX9) .OR. &
       (rcoarseMatrix%cmatrixFormat .NE. LSYSSC_MATRIX9)) THEN
      PRINT *,'mrest_matrixRestrictionEX3Y: Only format 9 matrices supported!'
      CALL sys_halt()
    END IF
    
    IF ((IAND(rfineMatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .NE. 0) .OR. &
        (IAND(rcoarseMatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .NE. 0)) THEN
      PRINT *,'mrest_matrixRestrictionEX3Y: Matrix must not be transposed!'
      CALL sys_halt()
    END IF
    
    IF ((rfineMatrix%p_rspatialDiscrTest%ccomplexity .NE. SPDISC_UNIFORM) .OR. &
        (rcoarseMatrix%p_rspatialDiscrTest%ccomplexity .NE. SPDISC_UNIFORM)) THEN
      PRINT *,'mrest_matrixRestrictionEX3Y: Only uniform discretisation supported!'
      CALL sys_halt()
    END IF

    IF ((rfineMatrix%cdataType .NE. ST_DOUBLE) .OR. &
        (rcoarseMatrix%cdataType .NE. ST_DOUBLE)) THEN
      PRINT *,'mrest_matrixRestrictionEX3Y: Only double precision matrices supported!'
      CALL sys_halt()
    END IF
    
    IF ((rfineMatrix%isortStrategy .GT. 0) .OR. &
        (rcoarseMatrix%isortStrategy .GT. 0)) THEN
      PRINT *,'mrest_matrixRestrictionEX3Y: Sorted matrices not supported!'
      CALL sys_halt()
    END IF

    IF ((rfineMatrix%dscaleFactor .NE. 1.0_DP) .OR. &
        (rcoarseMatrix%dscaleFactor .NE. 1.0_DP)) THEN
      PRINT *,'mrest_matrixRestrictionEX3Y: Scaled matrices not supported!'
      CALL sys_halt()
    END IF

    i1 = rfineMatrix%p_rspatialDiscrTrial%RelementDistr(1)%celement
    i2 = rcoarseMatrix%p_rspatialDiscrTrial%RelementDistr(1)%celement
    IF ((elem_getPrimaryElement(i1) .NE. EL_Q1T) .OR. &
        (elem_getPrimaryElement(i2) .NE. EL_Q1T)) THEN
      PRINT *,'mrest_matrixRestrictionEX3Y: Only Q1~-discretisation supported!'
      CALL sys_halt()
    END IF

    ! Looks good, so let's start.
    !
    ! Get information about the triangulation on the coarse and fine grid.
    p_rtriaCoarse => rcoarseMatrix%p_rspatialDiscrTest%p_rtriangulation
    p_rtriaFine => rfineMatrix%p_rspatialDiscrTest%p_rtriangulation

    ! Fetch all the information we need from the triangulation and the matrices
    ! for easier access:
    
    CALL storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                               p_IneighboursAtElementCoarse)
    CALL storage_getbase_int2d(p_rtriaCoarse%h_IedgesAtElement, &
                               p_IedgesAtElementCoarse)

    CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                               p_IneighboursAtElementFine)
    CALL storage_getbase_int2d(p_rtriaFine%h_IedgesAtElement, &
                               p_IedgesAtElementFine)

    CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                               p_IverticesAtElementCoarse)
                               
    CALL storage_getbase_double2d(p_rtriaCoarse%h_DvertexCoords, &
                                  p_DvertexCoordsCoarse)

    CALL lsyssc_getbase_Kld (rcoarseMatrix,p_KldCoarse)
    CALL lsyssc_getbase_Kcol (rcoarseMatrix,p_KcolCoarse)
    CALL lsyssc_getbase_Kdiagonal (rcoarseMatrix,p_KdiagonalCoarse)
    CALL lsyssc_getbase_double (rcoarseMatrix,p_DaCoarse)
    CALL lsyssc_getbase_Kld (rfineMatrix,p_KldFine)
    CALL lsyssc_getbase_Kcol (rfineMatrix,p_KcolFine)
    CALL lsyssc_getbase_Kdiagonal (rfineMatrix,p_KdiagonalFine)
    CALL lsyssc_getbase_double (rfineMatrix,p_DaFine)
    
    nvt1 = p_rtriaCoarse%NVT
    nvt2 = p_rtriaFine%NVT
    
    ! Loop through all elements of the coarse grid.

    DO iel=1,p_rtriaCoarse%NEL

      ! Get the numbers of the elements that are neighbours to our current coarse
      ! grid element:
      !             +--------+
      !             |        |
      !             | IELA3  |
      !             |        |
      !    +--------4--------3--------+
      !    |        |        |        |
      !    | IELA4  |  iel   | IELA2  |
      !    |        |        |        |
      !    +--------1--------2--------+
      !             |        |
      !             | IELA1  |
      !             |        |
      !             +--------+
      
      ielAdjacent(0) = iel
      ielAdjacent(1) = p_IneighboursAtElementCoarse(1,iel)
      ielAdjacent(2) = p_IneighboursAtElementCoarse(2,iel)
      ielAdjacent(3) = p_IneighboursAtElementCoarse(3,iel)
      ielAdjacent(4) = p_IneighboursAtElementCoarse(4,iel)

      ! Get the aspect ratio of the current element.

      dcoords = p_DvertexCoordsCoarse(:, &
                  p_IverticesAtElementCoarse(:,ielAdjacent(0)))
      daspectRatio(0) = gaux_getAspectRatio_quad2D (dcoords)
      IF (daspectRatio(0) .LT. 1.0_DP) daspectRatio(0) = 1.0_DP/daspectRatio(0)

      ! Get the aspect ratio of all the other neighbour elements - if they
      ! exist!
      DO iedge=1,TRIA_MAXNME2D
        IF (ielAdjacent(iedge) .NE. 0) THEN
          ! Get the aspect ratio of the current coarse grid element;
          ! if necessary, calculate the reciprocal.
          dcoords = p_DvertexCoordsCoarse(:, &
                      p_IverticesAtElementCoarse(:,ielAdjacent(iedge)))
          daspectRatio(iedge) = gaux_getAspectRatio_quad2D (dcoords)
          IF (daspectRatio(iedge) .LT. 1.0_DP) daspectRatio(iedge) = 1.0_DP/daspectRatio(iedge)
        ELSE
          daspectRatio(iedge) = 0.0_DP
        END IF
      END DO

      ! Get the elements on the fine grid that are inside of the
      ! current coarse grid element. We can easily fetch them because
      ! of the two-level ordering:

      Iiel(1) = iel
      Iiel(2) = p_IneighboursAtElementFine(2,Iiel(1))
      Iiel(3) = p_IneighboursAtElementFine(2,Iiel(2))
      Iiel(4) = p_IneighboursAtElementFine(2,Iiel(3))

      ! So we have:
      !
      !     4==============X===============3
      !     |              |               |
      !     |              |               |
      !     |    Iiel(4)   |     Iiel(3)   |
      !     |              |               |
      !     |                              |
      !     X----------- iel  -------------X 
      !     |                              |
      !     |              |               |
      !     |    Iiel(1)   |     Iiel(2)   |
      !     |              |               |
      !     |              |               |
      !     1==============X===============2
      !
      ! Now loop through the four edges/corners of the element
      ! iel on the coarse grid:

      edgeloop: DO iedge=1,4
      
        ! We assign:
        !
        ! iedge1 = local number of current edge
        ! iedge2 = local number of preceding edge in counterclockwise
        !        sense on that element
        ! iedge4 = local number of succeding edge in counterclockwise
        !        sense on that element

        iedge1=iedge
        iedge2=MOD(iedge1,4)+1
        iedge4=MOD(iedge1+2,4)+1

        ! i.e. in the example with iedge=1:
        !
        ! iedge4 O----------- iel  -------------O iedge2=2
        !  =4    |                              |
        !        |              |               |
        !        |    Iiel(1)   |     Iiel(2)   |
        !        |              |               |
        !        |              |               |
        !        ===============O================  
        !                    iedge1=1
        !
        ! Now, assign iadj1, iadj2 and iadj4 the numbers of the
        ! neighbouring elements on these three edges:
        
        iadj1 = ielAdjacent(iedge1)
        iadj2 = ielAdjacent(iedge2)
        iadj4 = ielAdjacent(iedge4)

        ! i.e.
        !
        !  +-------+-------+-------+
        !  |       |       |       |
        !  | iadj4 4  iel  2 iadj2 |
        !  |       |       |       |
        !  +-------+---1---+-------+
        !          |       |
        !          | iadj1 |
        !          |       |
        !          +-------+
        !
        ! iel1 and iel2 now receive the numbers of the elements on the fine
        ! grid that share the same edge iedge=iedge1:

        iel1=Iiel(iedge1)
        iel2=Iiel(iedge2)
      
        ! i.e.:
        !
        !      O------------ iel ------------O     
        !      |                             |
        !      |              |              |
        !      |     iel1     |     iel2     |
        !      |              |              |
        !      |              |              |
        !      1==============O==============2
        !                    iedge1
        !
        ! Now consider neighbouring element iadj1 of our current element iel:
        !
        !  +-------+-------+-------+
        !  |       |       |       |
        !  |       4 iel1  2       |
        !  |       |       |       |
        !  +-------+---1---+-------+
        !          |       |
        !          | iadj1 |
        !          |       |
        !          +-------+
        !
        ! In case iadj1 exists and has a smaller number than our current
        ! element, the matrix on that position is already modified, so we
        ! can skip the computation.

        IF ((iadj1 .LT. iel) .AND. (iadj1 .NE. 0)) CYCLE edgeloop

        ! In case a neighbour element exists and has a larger number,
        ! check the aspect ratio - of our current element as well as
        ! (if configured by the parameters) that of the neighbour element:

        IF (iadj1.NE.0) THEN

          ! In case our current element (and probably our neighbour) has to
          ! big jump in the aspect ratio, the matrix does not have to be modified, 
          ! we skip the modification.
          !
          ! Both aspect ratios in-range?

          IF ((daspectRatio(0) .LT. dARbound) .AND. &
              (daspectRatio (iedge) .LT. dARbound)) CYCLE edgeloop      
        
          ! At least one element is out-of-bounds.
          ! Check if it's the current element and if the neighbour element
          ! is important or not.

          IF ((iARindicator .EQ. 1) .AND. &
              (daspectRatio(0) .LT. darbound)) CYCLE edgeloop

        ELSE
      
          ! No neighbour. Only the aspect ratio of the current element
          ! decides on whether the matrix is modified or not. 

          IF (daspectRatio(0) .LT. dARbound) CYCLE edgeloop

        ENDIF
        
        ! Ok, we are in the case where the matrix must be modified
        ! because of too large anisotropy.
        !
        ! Get the global numbers of the edges imid1..3 corresponding to
        ! the local numbers iedge1,2,4 of the current element iel on the
        ! coarse grid.
        ! Subtract NVT to get the corresponding degree of freedom
        ! in the Q1~ discretisation.

        imid1 =p_IedgesAtElementCoarse(iedge1,iel)-nvt1
        imid2 =p_IedgesAtElementCoarse(iedge2,iel)-nvt1
        imid3 =p_IedgesAtElementCoarse(iedge4,iel)-nvt1

        ! i.e.:
        !
        ! imid3 O------------ iel ------------O imid2=iedge2
        ! =IVE3 |                             |
        !       |              |              |
        !       |     iel1     |     iel2     |
        !       |              |              |
        !       |              |              |
        !       1==============O==============2
        !                 imid1=iedge1
        !
        ! The same way, calculate the DOF's on the fine grid:

        im1 =p_IedgesAtElementFine(1,iel1)-nvt2
        im2 =p_IedgesAtElementFine(4,iel2)-nvt2
        im3 =p_IedgesAtElementFine(2,iel1)-nvt2
        im4 =p_IedgesAtElementFine(4,iel1)-nvt2
        im5 =p_IedgesAtElementFine(1,iel2)-nvt2
        im6 =p_IedgesAtElementFine(3,iel1)-nvt2
        im7 =p_IedgesAtElementFine(2,iel2)-nvt2
      
        ! i.e.:
        !
        !             im6            im7 
        !       +------O----- iel ----O-------+
        !       |4            3|3            2|
        !       |     iel1     |     iel2     |
        !   im4 O              O im3          O im5
        !       |              |              |
        !       |1            2|4            1|
        !       1=======O======+=======O======2
        !              im1            im2      
        !
        ! If there is a neighbour element on the coarse grid,
        ! we have to calculate the DOF's on that neighbour element:


        IF (iadj1.NE.0) THEN
      
          ! Loop through the four edges on the neighbour element
          ! to find the edge adjacent to our current element iel:
      
          DO jedge=1,4
            IF (p_IneighboursAtElementCoarse(jedge,iadj1).EQ.iel) EXIT
          END DO  

          !           4--------3
          !           |        |
          !           |  iel   |
          !           |   ^    |
          !  +--------1---^----2--------+
          !  |        |   ^    |        |
          !  |   ?    | iadj1  |   ?    |
          !  |        |        |        |
          !  +--------+--------+--------+
          !           |        |
          !           |   ?    |
          !           |        |
          !           +--------+
          !
          ! So jedge contains the local number (1..4) of the edge
          ! adjacent to iel and iadj1, looking from iadj1.
          !
          ! As in the case of iel, calculate the DOF's on that element.

          jedge1=jedge
          jedge2=MOD(jedge1,4)+1
          jedge4=MOD(jedge1+2,4)+1

          ! As an example, consider jedge=jedge1=3.
          ! Then, we calculate the left and right edge jedge2 and jedge4
          ! that correspond to the DOF's we have to take into account
          ! in the constant prolongation/restriction, i.e.:
          !
          ! iedge4 O----         iel1         ----O iedge2  
          !        |                              |
          !        |                              |
          !        |                              |
          !        |              |               |
          !        |              |               |
          !        4==============O===============3  
          !        |         jedge=jedge1=3       |
          !        |                              |
          !        |                              |
          !        |                              |
          !        |                              |
          ! jedge2 O                              O jedge4=2
          !   =4   |                              |
          !        |                              |
          !        |                              |
          !        |                              |
          !        |                              |
          !        1==============O===============2  
          !
          ! Get the global DOF's on the coarse grid corresponding to
          ! these two edges:
 
          imid4 =p_IedgesAtElementCoarse(jedge2,iadj1)-nvt1
          imid5 =p_IedgesAtElementCoarse(jedge4,iadj1)-nvt1

          ! i.e.
          !
          !         O----         iel          ----O       
          !         |                              |
          !         |                              |
          !         |                              |
          !         |              |               |
          !         |              |               |
          !         4==============O===============3  
          !         |         jedge=jedge1=3       |
          !         |                              |
          !         |                              |
          !         |                              |
          ! imid4   |                              | imid5
          ! =jedge2 O            iadj1             O =jedge4
          !         |                              |
          !         |                              |
          !         |                              |
          !         |                              |
          !         |                              |
          !         1==============O===============2  
          !
          ! Get the elements adjacent to iadj1 along these two edges:

          jadj2=p_IneighboursAtElementCoarse(jedge2,iadj1)
          jadj4=p_IneighboursAtElementCoarse(jedge4,iadj1)

          ! i.e. we have:
          !
          !           +--------+
          !           |   |    |
          !           |- iel --|
          !           |   |    |
          !  +--------+--------+--------+
          !  |        |        |        |
          !  | jadj2  | iadj1  | jadj4  |
          !  |        |        |        |
          !  +--------+--------+--------+
          !           |        |
          !           |   ?    |
          !           |        |
          !           +--------+
          !
          ! Now switch to the fine grid.
          ! Consider elements iel1 and iel2 in the coarse grid element iel.
          ! Find out the two fine grid elements on the neighbour element iadj1
          ! that touch iel1 and iel2:

          jel1 = p_IneighboursAtElementFine(4,iel2)
          jel2 = p_IneighboursAtElementFine(1,iel1)

          ! i.e.
          !
          !      O------------ iel -------------O       
          !      |                              |
          !      |              |               |
          !      |     iel1     |     iel2      |
          !      |      v       |      v        |
          !      |1     v       |      v       2|
          !      +======v=======O======v========+  
          !      |4     v       |      v       3|
          !      |      v       |      v        |
          !      |     jel2     |     jel1      |
          !      |              |               |
          !      |                              |
          !      O----------- iadj1 ------------O       
          !      |                              |
          !      |              |               |
          !      |              |               |
          !      |              |               |
          !      |1             |              2|
          !      +==============O===============+  C
          !
          ! Also on the fine grid, get the global numbers of the
          ! DOF's of the elements jel1 and jel2 that give a contribution
          ! to the constant prolongation/restriction:

          im8  = p_IedgesAtElementFine(2,jel1)-nvt2
          im9  = p_IedgesAtElementFine(4,jel1)-nvt2
          im10 = p_IedgesAtElementFine(1,jel2)-nvt2
          im11 = p_IedgesAtElementFine(3,jel1)-nvt2
          im12 = p_IedgesAtElementFine(2,jel2)-nvt2

          ! i.e.:
          !
          !       O----         iel          ----O       
          !       |                              |
          !       |                              |
          !       |                              |
          !       |              |               |
          !       |              |               |
          !       4==============+===============3  
          !       |1            4|2             1|
          !       |     jel2     |     jel1      |
          !  im10 O              O im8           O im9
          !       |              |               |
          !       |2            3|3             4|      
          !       +------O---- iadj1 ----O-------+      
          !       |4   im12     3|3    im11     2|
          !       |              |               |
          !       |              |               |
          !       |              |               |
          !       |1            2|4             1|
          !       1==============+===============2  


        ELSE

          ! In case there is no neighbour iadj1, set im8=0 to indicate that.
          im8 =0
          
        END IF
        
        ! Finally, the IMx variables now contain the numbers of the global 
        ! DOF's on the following edges:
        !
        !       |     im6             im7      |
        !       +------O----- iel ----O--------+
        !       |                              |
        !       |              |               |
        !   im4 O              O im3           O im5
        !       |              |               |
        !       |              |               |
        !       +=======O======+=======O=======+
        !       |      im1     |      im2      |
        !       |              |               |
        !  im10 O              O im8           O im9
        !       |              |               |
        !       |                              |      
        !       +------O---- iadj1 ----O-------+      
        !       |    im12            im11      |
        !
        ! Finally initialise some variables for later computations.

        dval1=0.0_DP
        dval2=0.0_DP
        dval3=0.0_DP
        dval4=0.0_DP
        dval5=0.0_DP

        ! Now use the values of the fine grid matrix at the positions IMx
        ! to cvalculate the entries in the coarse grid matrix at 
        ! imid1,...,imid5 with the Galerkin approach.
        !
        ! This corresponds to an application of the prolongation and 
        ! restriction operator to the fine grid matrix:

        ild = getXYindex (im1,im3,p_KcolFine,p_KldFine)
        dv1=p_DaFine(p_KdiagonalFine(im1))+p_DaFine(ild)
   
        IF (im8.NE.0) THEN
          ild = getXYindex (im1,im8,p_KcolFine,p_KldFine)
          dv1=dv1+p_DaFine(ild)
        END IF

        ild = getXYindex (im2,im3,p_KcolFine,p_KldFine)
        dv2=p_DaFine(p_KdiagonalFine(im2))+p_DaFine(ild)

        IF (im8.NE.0) THEN
          ild = getXYindex (im2,im8,p_KcolFine,p_KldFine)
          dv2=dv2+p_DaFine(ild)
        END IF

        ild = getXYindex (im3,im1,p_KcolFine,p_KldFine)
        dv3=p_DaFine(p_KdiagonalFine(im3))+p_DaFine(ild)

        ild = getXYindex (im3,im2,p_KcolFine,p_KldFine)
        dv3=dv3+p_DaFine(ild)

        ild = getXYindex (im4,im1,p_KcolFine,p_KldFine)
        dv4=p_DaFine(ild)

        ild = getXYindex (im4,im3,p_KcolFine,p_KldFine)
        dv4=dv4+p_DaFine(ild)

        ild = getXYindex (im5,im2,p_KcolFine,p_KldFine)
        dv5=p_DaFine(ild)

        ild = getXYindex (im5,im3,p_KcolFine,p_KldFine)
        dv5=dv5+p_DaFine(ild)

        ild = getXYindex (im6,im1,p_KcolFine,p_KldFine)
        dv6=p_DaFine(ild)

        ild = getXYindex (im6,im3,p_KcolFine,p_KldFine)
        dv6=dv6+p_DaFine(ild)

        ild = getXYindex (im7,im2,p_KcolFine,p_KldFine)
        dv7=p_DaFine(ild)

        ild = getXYindex (im7,im3,p_KcolFine,p_KldFine)
        dv7=dv7+p_DaFine(ild)

        IF (im8.NE.0) THEN
          ild = getXYindex (im8,im1,p_KcolFine,p_KldFine)
          dv8=p_DaFine(p_KdiagonalFine(im8))+p_DaFine(ild)

          ild = getXYindex (im8,im2,p_KcolFine,p_KldFine)
          dv8=dv8+p_DaFine(ild)

          ild = getXYindex (im9,im8,p_KcolFine,p_KldFine)
          dv9=p_DaFine(ild)

          ild = getXYindex (im9,im2,p_KcolFine,p_KldFine)
          dv9=dv9+p_DaFine(ild)

          ild = getXYindex (im10,im8,p_KcolFine,p_KldFine)
          dv10=p_DaFine(ild)

          ild = getXYindex (im10,im1,p_KcolFine,p_KldFine)
          dv10=dv10+p_DaFine(ild)

          ild = getXYindex (im11,im8,p_KcolFine,p_KldFine)
          dv11=p_DaFine(ild)

          ild = getXYindex (im11,im2,p_KcolFine,p_KldFine)
          dv11=dv11+p_DaFine(ild)

          ild = getXYindex (im12,im8,p_KcolFine,p_KldFine)
          dv12=p_DaFine(ild)

          ild = getXYindex (im12,im1,p_KcolFine,p_KldFine)
          dv12=dv12+p_DaFine(ild)
        ENDIF

        IF (iadj1.EQ.0) THEN
          dval1=1.0_DP/1.0_DP*(dv1+dv2+dv3)
        ELSE
          dval1=1.0_DP/1.0_DP*(dv1+dv2+dv3+dv8)
        ENDIF

        IF (iadj2.EQ.0) THEN
          dval2=1.0_DP/1.0_DP*(dv5+dv7)
        ELSE
          dval2=1.0_DP/1.0_DP*(dv5+dv7)
        ENDIF

        IF (iadj4.EQ.0) THEN
          dval3=1.0_DP/1.0_DP*(dv4+dv6)
        ELSE
          dval3=1.0_DP/1.0_DP*(dv4+dv6)
        ENDIF

        IF (iadj1.NE.0) THEN

          IF (jadj2.EQ.0) THEN
            dval4=1.0_DP/1.0_DP*(dv10+dv12)
          ELSE
            dval4=1.0_DP/1.0_DP*(dv10+dv12)
          ENDIF

          IF (jadj4.EQ.0) THEN
            dval5=1.0_DP/1.0_DP*(dv9+dv11)
          ELSE
            dval5=1.0_DP/1.0_DP*(dv9+dv11)
          ENDIF

        ENDIF

        ! Calculate the actual entries in the coarse grid matrix.
        ! First, clear the row imid1 in the matrix.
        p_DaCoarse(p_KldCoarse(imid1):p_KldCoarse(imid1+1)-1) = 0.0_DP

        ! Rebuild the entries in that row:

        p_DaCoarse(p_KdiagonalCoarse(imid1))=dval1

        ild = getXYindex (imid2,imid1,p_KcolCoarse,p_KldCoarse)
        p_DaCoarse(ild)=dval2
      
        ild = getXYindex (imid3,imid1,p_KcolCoarse,p_KldCoarse)
        p_DaCoarse(ild)=dval3

        IF (iadj1.NE.0) THEN
          ild = getXYindex (imid4,imid1,p_KcolCoarse,p_KldCoarse)
          p_DaCoarse(ild)=dval4

          ild = getXYindex (imid5,imid1,p_KcolCoarse,p_KldCoarse)
          p_DaCoarse(ild)=dval5
        ENDIF
        
        ! We finished with the current edge imid1. Switch to the
        ! next edge in counterclockwise sense and go on.
        
      END DO edgeloop

      ! Current element finished. Proceed with next element

    END DO

    ! That's it.
    
  CONTAINS

    ! -----------------------------------------------------------------------
    ! Auxiliary routine: Get X/Y-Index
    !
    ! Performs a search for the index ild such that matrix(IX,IY)=KLA(ild) 
    ! holds. I.e. searches in the matrix array for the index belonging
    ! to the position IX/IY, so that the caller can modify this matrix
    ! element directly.
    ! -----------------------------------------------------------------------

    PURE INTEGER FUNCTION getXYindex (IX, IY, Kcol, Kld)

      ! input: column number to search for
      INTEGER(PREC_VECIDX), INTENT(IN) :: IX

      ! input: row number where to search
      INTEGER(PREC_VECIDX), INTENT(IN) :: IY

      ! input: Column structure of the matrix
      INTEGER(PREC_VECIDX), DIMENSION(:), INTENT(IN) :: Kcol
      
      ! input: Row structure of the matrix
      INTEGER(PREC_MATIDX), DIMENSION(:), INTENT(IN) :: Kld
    
      ! result: index of entry (ix,iy) in the matrix array.
      ! =-1, if the entry does not exist.
    
      ! local variables:
      INTEGER(PREC_MATIDX) :: ild
      INTEGER(PREC_VECIDX) :: icol
      
      ! Look through row IY:

      DO ild=Kld(IY),Kld(IY+1)-1
        icol=Kcol(ild)
        ! If there's column IX in this row we can stop here
        IF (icol .EQ. IX) THEN
          getXYindex = ild
          RETURN
        END IF
      END DO
      
      ! Otherwise: error - this element does not exist in our matrix
      ! indicate that by returning -1. This will usually result in an 'array
      ! of bounds exception' in our caller if compiled in DEBUG mode.
      
      getXYindex = -1

    END FUNCTION

  END SUBROUTINE

END MODULE
