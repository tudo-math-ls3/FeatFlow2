!##############################################################################
!# ****************************************************************************
!# <name> matrixrestriction </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# This module contains routines to realise a matrix restriction, i.e. the
!# construction of coarse grid matrix operators from fine grid matrices by the
!# Galerkin approach (c.f. p. 224ff in Turek`s book).
!#
!# In the Galerkin approach a coarse grid matrix <tex>$A_{2h}$</tex>
!# is build from a fine grid matrix $A_h$ with the help of the prolongation
!# operator <tex>$P_{2h}^h$</tex> and the restriction operator
!# <tex>$R_{h}^{2h}$</tex> using
!#
!#       <tex> $$ A_{2h}  =  P_{2h}^h  A_h  R_{h}^{2h}. $$ </tex>
!#
!# This is equivalent to the standard approach
!# <tex>$a_{ij} = a_h(\phi_j,\phi_i)$</tex>
!# in the case of conformal finite elements, but yields a different matrix
!# in the case of nonconformal FE spaces.
!#
!# The routines in this file allow for some situations (discretisation
!# with <tex>$\tilde Q_1$</tex> e.g.) to rebuild some lines in a given matrix
!# according to this approach. This is a kind of stabilisation in case the mesh
!# elements show high aspect ratios.
!#
!# The following subroutines can be found here:
!#
!# 1.) mrest_matrixRestrictionEX3Y
!#     -> Build coarse grid matrix entries on elements with too high
!#        aspect ratio by using a fine grid matrix.
!#
!# </purpose>
!##############################################################################

module matrixrestriction

!$use omp_lib
  use fsystem
  use storage
  use genoutput
  use matrixmodification
  use basicgeometry
  use geometryaux
  use triangulation
  use element
  use spatialdiscretisation
  use linearsystemscalar

  implicit none

  private

  public :: mrest_matrixRestrictionEX3Y

contains

  ! ***************************************************************************

!<subroutine>

  subroutine mrest_matrixRestrictionEX3Y (rfineMatrix,rcoarseMatrix,&
                                          iARindicator, dARbound)

!<description>
  ! This routine rebuilds some line in the scalar coarse grid matrix
  ! rcoarseMatrix with the help of the Galerkin approach using the scalar
  ! fine grid matrix rfineMatrix. Only the lines of those DOF`s will be
  ! rebuild where the elements that belong to the DOF show too high aspect
  ! ratios.
  !
  ! The routine works only for matrices that are build using a uniform
  ! discretisation with <tex>$\tilde Q_1$</tex>.
!</description>

!<input>
  ! Fine grid matrix.
  type(t_matrixScalar), intent(in) :: rfineMatrix

  ! Aspect-ratio indicator. Configures when to rebuild a line of the
  ! coarse matrix by constant restriction.
  ! <=0: do nothing, do not modify coarse grid matrix
  !  =1: switch depending on aspect ratio of current element
  !  =2: switch depending on aspect ratio of current element
  !      and aspect ratio of neighbour element
  integer, intent(in) :: iARindicator

  ! Maximum allowed aspect ratio. Rows in the matrix corresponding
  ! to elements (i.e. a DOF`s of an element) with an aspect ratio
  ! larger than dARbound are rebuild with the Galerkin approach
  ! by constant restriction.
  real(DP), intent(in) :: dARbound

!</input>

!<inputoutput>
  ! Coarse grid matrix to be modified.
  type(t_matrixScalar), intent(in) :: rcoarseMatrix
!</inputoutput>

!</subroutine>

    ! local variables
    integer(I32) :: i1,i2
    integer :: iedge,iedge1,iedge2,iedge4,jedge,jedge1,jedge2,jedge4
    integer :: iel,iadj1,iadj2,iadj4,iel1,iel2,jadj2,jadj4,jel1,jel2
    integer :: imid1,imid2,imid3,imid4,imid5
    integer :: im1,im2,im3,im4,im5,im6,im7,im8,im9,im10,im11,im12
    integer :: nvt1,nvt2
    integer :: ild
    type(t_triangulation), pointer :: p_rtriaCoarse,p_rtriaFine
    real(DP), dimension(NDIM2D,TRIA_MAXNVE2D) :: dcoords
    real(DP), dimension(0:TRIA_MAXNME2D) :: daspectRatio
    integer, dimension(TRIA_MAXNME2D) :: Iiel
    integer, dimension(0:TRIA_MAXNME2D) :: ielAdjacent
    real(DP) :: dval1,dval2,dval3,dval4,dval5
    integer, dimension(:,:), pointer :: p_IneighboursAtElementCoarse
    integer, dimension(:,:), pointer :: p_IneighboursAtElementFine
    real(DP), dimension(:,:), pointer                 :: p_DvertexCoordsCoarse
    integer, dimension(:,:), pointer    :: p_IedgesAtElementCoarse
    integer, dimension(:,:), pointer    :: p_IedgesAtElementFine
    integer, dimension(:,:), pointer   :: p_IverticesAtElementCoarse

    integer, dimension(:), pointer       :: p_KldCoarse,p_KldFine
    integer, dimension(:), pointer       :: p_KdiagonalCoarse,p_KdiagonalFine
    integer, dimension(:), pointer       :: p_KcolCoarse,p_KcolFine
    real(DP), dimension(:), pointer                   :: p_DaCoarse,p_DaFine
    real(DP) :: dv1,dv2,dv3,dv4,dv5,dv6,dv7,dv8,dv9,dv10,dv11,dv12

    ! No modification if the parameters are out of bounds!
    if (iARindicator .le. 0) return
    if (dARbound .lt. 0.0_DP) return

    ! At first some basic checks if we are able to complete our task at all.

    if ((rfineMatrix%cmatrixFormat .ne. LSYSSC_MATRIX9) .or. &
       (rcoarseMatrix%cmatrixFormat .ne. LSYSSC_MATRIX9)) then
      call output_line ("Only format 9 matrices supported!",&
          OU_CLASS_ERROR,OU_MODE_STD,"mrest_matrixRestrictionEX3Y")
      call sys_halt()
    end if

    if ((iand(rfineMatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0) .or. &
        (iand(rcoarseMatrix%imatrixSpec,LSYSSC_MSPEC_TRANSPOSED) .ne. 0)) then
      call output_line ("Matrix must not be transposed!",&
          OU_CLASS_ERROR,OU_MODE_STD,"mrest_matrixRestrictionEX3Y")
      call sys_halt()
    end if

    if ((rfineMatrix%p_rspatialDiscrTest%ccomplexity .ne. SPDISC_UNIFORM) .or. &
        (rcoarseMatrix%p_rspatialDiscrTest%ccomplexity .ne. SPDISC_UNIFORM)) then
      call output_line ("Only uniform discretisation supported!",&
          OU_CLASS_ERROR,OU_MODE_STD,"mrest_matrixRestrictionEX3Y")
      call sys_halt()
    end if

    if ((rfineMatrix%cdataType .ne. ST_DOUBLE) .or. &
        (rcoarseMatrix%cdataType .ne. ST_DOUBLE)) then
      call output_line ("Only double precision matrices supported!",&
          OU_CLASS_ERROR,OU_MODE_STD,"mrest_matrixRestrictionEX3Y")
      call sys_halt()
    end if

    if ((rfineMatrix%bcolumnsSorted .or. rcoarseMatrix%bcolumnsSorted) .or. &
        (rfineMatrix%browsSorted .or. rcoarseMatrix%browsSorted)) then
      call output_line ("Sorted matrices not supported!",&
          OU_CLASS_ERROR,OU_MODE_STD,"mrest_matrixRestrictionEX3Y")
      call sys_halt()
    end if

    if ((rfineMatrix%dscaleFactor .ne. 1.0_DP) .or. &
        (rcoarseMatrix%dscaleFactor .ne. 1.0_DP)) then
      call output_line ("Scaled matrices not supported!",&
          OU_CLASS_ERROR,OU_MODE_STD,"mrest_matrixRestrictionEX3Y")
      call sys_halt()
    end if

    i1 = rfineMatrix%p_rspatialDiscrTrial%RelementDistr(1)%celement
    i2 = rcoarseMatrix%p_rspatialDiscrTrial%RelementDistr(1)%celement
    if ((elem_getPrimaryElement(i1) .ne. EL_Q1T) .or. &
        (elem_getPrimaryElement(i2) .ne. EL_Q1T)) then
      call output_line ("Only Q1~-discretisation supported!",&
          OU_CLASS_ERROR,OU_MODE_STD,"mrest_matrixRestrictionEX3Y")
      call sys_halt()
    end if

    ! Looks good, so let us start.
    !
    ! Get information about the triangulation on the coarse and fine grid.
    p_rtriaCoarse => rcoarseMatrix%p_rspatialDiscrTest%p_rtriangulation
    p_rtriaFine => rfineMatrix%p_rspatialDiscrTest%p_rtriangulation

    ! Fetch all the information we need from the triangulation and the matrices
    ! for easier access:

    call storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                               p_IneighboursAtElementCoarse)
    call storage_getbase_int2d(p_rtriaCoarse%h_IedgesAtElement, &
                               p_IedgesAtElementCoarse)

    call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                               p_IneighboursAtElementFine)
    call storage_getbase_int2d(p_rtriaFine%h_IedgesAtElement, &
                               p_IedgesAtElementFine)

    call storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                               p_IverticesAtElementCoarse)

    call storage_getbase_double2d(p_rtriaCoarse%h_DvertexCoords, &
                                  p_DvertexCoordsCoarse)

    call lsyssc_getbase_Kld (rcoarseMatrix,p_KldCoarse)
    call lsyssc_getbase_Kcol (rcoarseMatrix,p_KcolCoarse)
    call lsyssc_getbase_Kdiagonal (rcoarseMatrix,p_KdiagonalCoarse)
    call lsyssc_getbase_double (rcoarseMatrix,p_DaCoarse)
    call lsyssc_getbase_Kld (rfineMatrix,p_KldFine)
    call lsyssc_getbase_Kcol (rfineMatrix,p_KcolFine)
    call lsyssc_getbase_Kdiagonal (rfineMatrix,p_KdiagonalFine)
    call lsyssc_getbase_double (rfineMatrix,p_DaFine)

    nvt1 = p_rtriaCoarse%NVT
    nvt2 = p_rtriaFine%NVT

    ! Loop through all elements of the coarse grid.

    do iel=1,p_rtriaCoarse%NEL

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
      if (daspectRatio(0) .lt. 1.0_DP) daspectRatio(0) = 1.0_DP/daspectRatio(0)

      ! Get the aspect ratio of all the other neighbour elements - if they
      ! exist!
      do iedge=1,TRIA_MAXNME2D
        if (ielAdjacent(iedge) .ne. 0) then
          ! Get the aspect ratio of the current coarse grid element;
          ! if necessary, calculate the reciprocal.
          dcoords = p_DvertexCoordsCoarse(:, &
                      p_IverticesAtElementCoarse(:,ielAdjacent(iedge)))
          daspectRatio(iedge) = gaux_getAspectRatio_quad2D (dcoords)
          if (daspectRatio(iedge) .lt. 1.0_DP) daspectRatio(iedge) = 1.0_DP/daspectRatio(iedge)
        else
          daspectRatio(iedge) = 0.0_DP
        end if
      end do

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

      edgeloop: do iedge=1,4

        ! We assign:
        !
        ! iedge1 = local number of current edge
        ! iedge2 = local number of preceding edge in counterclockwise
        !        sense on that element
        ! iedge4 = local number of succeding edge in counterclockwise
        !        sense on that element

        iedge1=iedge
        iedge2=mod(iedge1,4)+1
        iedge4=mod(iedge1+2,4)+1

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

        if ((iadj1 .lt. iel) .and. (iadj1 .ne. 0)) cycle edgeloop

        ! In case a neighbour element exists and has a larger number,
        ! check the aspect ratio - of our current element as well as
        ! (if configured by the parameters) that of the neighbour element:

        if (iadj1.ne.0) then

          ! In case our current element (and probably our neighbour) has to
          ! big jump in the aspect ratio, the matrix does not have to be modified,
          ! we skip the modification.
          !
          ! Both aspect ratios in-range?

          if ((daspectRatio(0) .lt. dARbound) .and. &
              (daspectRatio (iedge) .lt. dARbound)) cycle edgeloop

          ! At least one element is out-of-bounds.
          ! Check if it is the current element and if the neighbour element
          ! is important or not.

          if ((iARindicator .eq. 1) .and. &
              (daspectRatio(0) .lt. darbound)) cycle edgeloop

        else

          ! No neighbour. Only the aspect ratio of the current element
          ! decides on whether the matrix is modified or not.

          if (daspectRatio(0) .lt. dARbound) cycle edgeloop

        endif

        ! Ok, we are in the case where the matrix must be modified
        ! because of too large anisotropy.
        !
        ! Get the global numbers of the edges imid1..3 corresponding to
        ! the local numbers iedge1,2,4 of the current element iel on the
        ! coarse grid.
        ! Subtract NVT to get the corresponding degree of freedom
        ! in the Q1~ discretisation.

        imid1 =p_IedgesAtElementCoarse(iedge1,iel)
        imid2 =p_IedgesAtElementCoarse(iedge2,iel)
        imid3 =p_IedgesAtElementCoarse(iedge4,iel)

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
        ! The same way, calculate the DOF`s on the fine grid:

        im1 =p_IedgesAtElementFine(1,iel1)
        im2 =p_IedgesAtElementFine(4,iel2)
        im3 =p_IedgesAtElementFine(2,iel1)
        im4 =p_IedgesAtElementFine(4,iel1)
        im5 =p_IedgesAtElementFine(1,iel2)
        im6 =p_IedgesAtElementFine(3,iel1)
        im7 =p_IedgesAtElementFine(2,iel2)

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
        ! we have to calculate the DOF`s on that neighbour element:


        if (iadj1.ne.0) then

          ! Loop through the four edges on the neighbour element
          ! to find the edge adjacent to our current element iel:

          do jedge=1,4
            if (p_IneighboursAtElementCoarse(jedge,iadj1).eq.iel) exit
          end do

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
          ! As in the case of iel, calculate the DOF`s on that element.

          jedge1=jedge
          jedge2=mod(jedge1,4)+1
          jedge4=mod(jedge1+2,4)+1

          ! As an example, consider jedge=jedge1=3.
          ! Then, we calculate the left and right edge jedge2 and jedge4
          ! that correspond to the DOF`s we have to take into account
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
          ! Get the global DOF`s on the coarse grid corresponding to
          ! these two edges:

          imid4 =p_IedgesAtElementCoarse(jedge2,iadj1)
          imid5 =p_IedgesAtElementCoarse(jedge4,iadj1)

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
          ! DOF`s of the elements jel1 and jel2 that give a contribution
          ! to the constant prolongation/restriction:

          im8  = p_IedgesAtElementFine(2,jel1)
          im9  = p_IedgesAtElementFine(4,jel1)
          im10 = p_IedgesAtElementFine(1,jel2)
          im11 = p_IedgesAtElementFine(3,jel1)
          im12 = p_IedgesAtElementFine(2,jel2)

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


        else

          ! In case there is no neighbour iadj1, set im8=0 to indicate that.
          im8 =0

        end if

        ! Finally, the IMx variables now contain the numbers of the global
        ! DOF`s on the following edges:
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

        if (im8.ne.0) then
          ild = getXYindex (im1,im8,p_KcolFine,p_KldFine)
          dv1=dv1+p_DaFine(ild)
        end if

        ild = getXYindex (im2,im3,p_KcolFine,p_KldFine)
        dv2=p_DaFine(p_KdiagonalFine(im2))+p_DaFine(ild)

        if (im8.ne.0) then
          ild = getXYindex (im2,im8,p_KcolFine,p_KldFine)
          dv2=dv2+p_DaFine(ild)
        end if

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

        if (im8.ne.0) then
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
        endif

        if (iadj1.eq.0) then
          dval1=1.0_DP/1.0_DP*(dv1+dv2+dv3)
        else
          dval1=1.0_DP/1.0_DP*(dv1+dv2+dv3+dv8)
        endif

        if (iadj2.eq.0) then
          dval2=1.0_DP/1.0_DP*(dv5+dv7)
        else
          dval2=1.0_DP/1.0_DP*(dv5+dv7)
        endif

        if (iadj4.eq.0) then
          dval3=1.0_DP/1.0_DP*(dv4+dv6)
        else
          dval3=1.0_DP/1.0_DP*(dv4+dv6)
        endif

        if (iadj1.ne.0) then

          if (jadj2.eq.0) then
            dval4=1.0_DP/1.0_DP*(dv10+dv12)
          else
            dval4=1.0_DP/1.0_DP*(dv10+dv12)
          endif

          if (jadj4.eq.0) then
            dval5=1.0_DP/1.0_DP*(dv9+dv11)
          else
            dval5=1.0_DP/1.0_DP*(dv9+dv11)
          endif

        endif

        ! Calculate the actual entries in the coarse grid matrix.
        ! First, clear the row imid1 in the matrix.
        p_DaCoarse(p_KldCoarse(imid1):p_KldCoarse(imid1+1)-1) = 0.0_DP

        ! Rebuild the entries in that row:

        p_DaCoarse(p_KdiagonalCoarse(imid1))=dval1

        ild = getXYindex (imid2,imid1,p_KcolCoarse,p_KldCoarse)
        p_DaCoarse(ild)=dval2

        ild = getXYindex (imid3,imid1,p_KcolCoarse,p_KldCoarse)
        p_DaCoarse(ild)=dval3

        if (iadj1.ne.0) then
          ild = getXYindex (imid4,imid1,p_KcolCoarse,p_KldCoarse)
          p_DaCoarse(ild)=dval4

          ild = getXYindex (imid5,imid1,p_KcolCoarse,p_KldCoarse)
          p_DaCoarse(ild)=dval5
        endif

        ! We finished with the current edge imid1. Switch to the
        ! next edge in counterclockwise sense and go on.

      end do edgeloop

      ! Current element finished. Proceed with next element

    end do

    ! That is it.

  contains

    ! -----------------------------------------------------------------------
    ! Auxiliary routine: Get X/Y-Index
    !
    ! Performs a search for the index ild such that matrix(IX,IY)=KLA(ild)
    ! holds. I.e. searches in the matrix array for the index belonging
    ! to the position IX/IY, so that the caller can modify this matrix
    ! element directly.
    ! -----------------------------------------------------------------------

    pure integer function getXYindex (IX, IY, Kcol, Kld)

      ! input: column number to search for
      integer, intent(in) :: IX

      ! input: row number where to search
      integer, intent(in) :: IY

      ! input: Column structure of the matrix
      integer, dimension(:), intent(in) :: Kcol

      ! input: Row structure of the matrix
      integer, dimension(:), intent(in) :: Kld

      ! result: index of entry (ix,iy) in the matrix array.
      ! =-1, if the entry does not exist.

      ! local variables:
      integer :: ild
      integer :: icol

      ! Look through row IY:

      do ild=Kld(IY),Kld(IY+1)-1
        icol=Kcol(ild)
        ! If there is column IX in this row we can stop here
        if (icol .eq. IX) then
          getXYindex = ild
          return
        end if
      end do

      ! Otherwise: error - this element does not exist in our matrix
      ! indicate that by returning -1. This will usually result in an 'array
      ! of bounds exception' in our caller if compiled in DEBUG mode.

      getXYindex = -1

    end function

  end subroutine

end module
