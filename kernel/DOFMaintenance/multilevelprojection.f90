!##############################################################################
!# ****************************************************************************
!# <name> multilevelprojection </name>
!# ****************************************************************************
!#
!# <purpose>
!# Contains routines for prolongation, restriction and interpolation
!# of solution vectors between different levels, for scalar systems as well
!# as for block systems.
!#
!# The module contains the following routines:
!#
!# 1.) mlprj_initProjectionDirect /
!#     mlprj_initProjectionDiscr /
!#     mlprj_initProjectionVec /
!#     mlprj_initProjectionMat
!#     -> Initialises a projection structure with values according to a
!#        spatial discretisation.
!#        Uses a discretisation structure, a (block) vector, a (block)
!#        matrix on the fine/coarse grid or a block discretisation structure
!#        as template.
!#
!# 2.) mlprj_doneProjection
!#     -> Cleans up a projection structure.
!#
!# 3.) mlprj_getTempMemoryDirect /
!#     mlprj_getTempMemoryVec /
!#     mlprj_getTempMemoryMat
!#     -> Determine the amount of temporary storage that is necessary
!#        to transport a vector from one level to another
!#        Uses a discretisation structure, a (block) vector or a (block)
!#        matrix on the fine/coarse grid as template.
!#
!# 4.) mlprj_performProlongation
!#     -> Interpolates a solution vector from a coarse grid to a fine grid
!#        (L2-projection in the primal space)
!#
!# 5.) mlprj_performRestriction
!#     -> Restricts a defect vector from a fine grid to a coarse grid
!#        (L2-projection in the dual space)
!#
!# 6.) mlprj_performInterpolation
!#     -> Interpolated a solution from a fine grid to a coarse grid
!#        (L2-projection in the primal space)
!#
!# 7.) mlprj_initL2Projection
!#     -> Sets the matrices for an scalar projection structure which are
!#        needed for L2-projection.
!#
!# 8.) mlprj_initMatrixProjection
!#     -> Directly sets the prolongation (and thus implicitly also the
!#        restriction) matrix that is to be used.
!#
!# 9.) mlprj_initPrjHierarchy
!#     -> Creates a projection hierarchy containing a set of levels.
!#
!# 10.) mlprj_initPrjHierarchyLevel
!#      -> Initialises a level in a projection hierarchy
!#
!# 11.) mlprj_commitPrjHierarchy
!#      -> Prepares a projection hierarchy for use
!#
!# 12.) mlprj_releasePrjHierarchy
!#      -> Releases a projection hierarchy
!#
!# 13.) mlprj_performProlongationHier
!#      -> Performs a prolongation using a projection hierarchy
!#
!# 14.) mlprj_performRestrictionHier
!#      -> Performs a restriction using a projection hierarchy
!#
!# 15.) mlprj_performInterpolationHier
!#      -> Performs an interpolation using a projection hierarchy
!#
!# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- \\
!#  Some notes on the L2-Projection grid transfer  \\
!# ----------------------------------------------- \\
!# The possibility to choose "real" L2-projection for multi-level grid
!# transfer has been added mainly for research purposes, especially with
!# non-conforming elements.
!#
!# In contrast to the "hard-coded" grid transfer this method has 2 advantages:
!#
!# +1. It is a black-box projection, i.e. one can apply the multigrid solver
!#     on discretisations for which the developer has not (yet) written hard-
!#     coded grid transfer operators.
!#
!# +2. It provides true L2-orthogonality also for non-conforming elements.
!#     The "hard-coded" grid transfer for non-conforming elements is just an
!#     approximation of the "true" L2-projection operator, thus usually giving
!#     worse results.
!#
!# But it also has 2 major disadvantages:
!#
!# -1. The user also needs to assemble a mass (on the fine grid) and a
!#     2-level-mass matrix for the grid transfer and pass it to the
!#     mlprj_initL2Projection routine. The mass matrix can be assembled using
!#     the bilinearformevaluation. The 2-level-mass matrix needs to be
!#     assembled by the multileveloperators module which, until now, supports
!#     only uniform discretisations (and no triangular meshes).
!#
!# -2. The grid transfer is (in contrast to the "hard-coded" projection)
!#     very expensive and time consuming, as each prolongation and each
!#     restriction involves solving a linear system with the fine mesh mass
!#     matrix.
!#
!# <verb>
!# ------------
!# - WARNING! -
!# ------------
!# </verb>
!# There is no support for the "mlprj_performInterpolation" routine using
!# L2-projection, and it is NOT planned to implement it! If you want to
!# restrict a solution vector using L2-projection you will need to do it by
!# yourself
!#
!# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- \\
!#  Some notes on the matrix-based grid transfer  \\
!# ---------------------------------------------- \\
!# There also is a possibility to directly specify a matrix for grid transfer.
!# In this case, when prolongating a coarse mesh solution, the prolongation
!# matrix is multiplied by the coarse mesh solution vector to get the fine
!# mesh solution vector. When restricting a defect vector, the transposed
!# prolongation matrix is multiplied by the fine mesh defect vector to get
!# the coarse mesh defect vector.
!#
!# <verb>
!# ------------
!# - WARNING! -
!# ------------
!# </verb>
!# There is no support for the "mlprj_performInterpolation" routine using
!# matrix-based grid transfer, and it is NOT planned to implement it! If you
!# want to restrict a solution vector using a matrix you will need to do it
!# by yourself
!#
!# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- \\
!#  Some notes on the projection hierarchy  \\
!# ---------------------------------------- \\
!#
!# The projection hierarchy can be used to maintain the projection
!# between a set of levels. For its use, one has to prepare the structure
!# with a set of commands like this:
!#
!# <code>
!#   type(t_interlevelProjectionHier) :: rprjHierarchy
!#
!#   call mlprj_initPrjHierarchy(rprjHierarchy,nlmin,nlmax)
!#   do ilevel = nlmin,nlmax
!#     call mlprj_initPrjHierarchyLevel(rprjHierarchy,ilevel,Rdiscretisation(ilevel))
!#   end do
!#   ...
!#   ... do something else with rprjHierarchy if necessary, e.g. preparing
!#   ... rprjHierarchy%p_Rprojection(ilevel) with special parameters
!#   ...
!#   call mlprj_commitPrjHierarchy (rprjHierarchy)
!# </code>
!#
!# The above code initialises the hierarchy for levels nlmin..nlmax,
!# attaches the discretisation structures of these levels and
!# "commits" the structure. Afterwards it can he used with
!# mlprj_performProlongationHier/mlprj_performRestrictionHier/
!# mlprj_performInterpolationHier to projection solution/RHS vectors between
!# adjacent levels. These routines accept a coarse grid/fine grid vector
!# and the number of the level corresponding to the fine grid vector.
!#
!# When the computation is finished, the whole structure must be released
!# with
!#
!# <code>
!#   call mlprj_releasePrjHierarchy (rprjHierarchy)
!# </code>
!#
!# </purpose>
!##############################################################################

module multilevelprojection

!$use omp_lib
  use fsystem
  use storage
  use basicgeometry
  use genoutput
  use spatialdiscretisation
  use linearsystemscalar
  use linearsystemblock
  use triangulation
  use geometryaux
  use element
  use quicksolver
  use linearalgebra
  use dofmapping

  implicit none

  private

  public :: t_interlevelProjectionDiscrPtr
  public :: t_interlevelProjectionScalar
  public :: t_interlevelProjectionBlock
  public :: t_interlevelProjectionHier
  public :: mlprj_initProjectionDirect
  public :: mlprj_initProjectionDiscr
  public :: mlprj_initProjectionVec
  public :: mlprj_initProjectionMat
  public :: mlprj_initL2Projection
  public :: mlprj_initMatrixProjection
  public :: mlprj_doneProjection
  public :: mlprj_getTempMemoryScalar
  public :: mlprj_getTempMemoryDirect
  public :: mlprj_getTempMemoryVec
  public :: mlprj_getTempMemoryMat
  public :: mlprj_performProlongation
  public :: mlprj_performRestriction
  public :: mlprj_performInterpolation
  public :: mlprj_initPrjHierarchy
  public :: mlprj_initPrjHierarchyLevel
  public :: mlprj_commitPrjHierarchy
  public :: mlprj_releasePrjHierarchy
  public :: mlprj_performProlongationHier
  public :: mlprj_performRestrictionHier
  public :: mlprj_performInterpolationHier

!<constants>

!<constantblock description="Projection Types">
  ! Hard-coded projection operators
  integer(I32), parameter, public :: MLP_PROJ_TYPE_HARDCODED = 0_I32

  ! L2-projection operators
  integer(I32), parameter, public :: MLP_PROJ_TYPE_L2_PROJ   = 1_I32

  ! Matrix-based projection operators
  integer(I32), parameter, public :: MLP_PROJ_TYPE_MATRIX    = 2_I32

!</constantblock>

!</constants>

!<types>

!<typeblock>

  ! Describes for a scalar equation which type of prolongation, restriction
  ! or interpolation to be used. The information is (of course) element
  ! dependent.
  ! One t_interlevelProjectionScalar structure is associated to one
  ! element group of a spatial discretisation.
  ! t_interlevelProjectionScalar is not part if the t_spatialDiscretisation
  ! structure, as t_interlevelProjectionScalar resides "between" levels:
  ! If there is one t_spatialDiscretisation structure for level 3 and
  ! one for level 4, the t_interlevelProjectionScalar structure is
  ! logically arranged at level "3.5" to configure the projection
  ! between level 3 and 4. For this reason, there is one
  ! t_interlevelProjectionScalar less than levels: As there is no level 0,
  ! the structure for level "0.5" is missing!

  type t_interlevelProjectionScalar

    ! Specifies the projection type of this projection structure.
    ! One of the MLP_PROJ_TYPE_XXXX constants defined above.
    integer(I32)                :: iprojType = MLP_PROJ_TYPE_HARDCODED

    ! Element type that should be assumed for the prolongation. Must fit
    ! to the element type of the discretisation concerning the DOF`s:
    ! If element "EM30" is used e.g., prolongation for "EM30", "EM31",
    ! "E031" and "E030" is possible, although performance will suffer
    ! in using the "wrong" prolongation. An error will be displayed if
    ! the DOF`s do not fit together at all, e.g. when using <tex>$Q_1$</tex>
    ! prolongation when using <tex>$Q_2$</tex> for the discretisation.
    !
    ! A value of EL_UNDEFINED indicates that the "natural" prolongation
    ! (i.e. that configured by the spatial discretisation) should be used.
    integer(I32)                :: ielementTypeProlongation = EL_UNDEFINED

    ! Order of the prolongation to use.
    ! -1=use default prolongation (e.g. linear for <tex>$Q_1$</tex>, quadratic for
    !    <tex>$Q_2$</tex>,...). Some elements allow to configure the type of the
    !    prolongation, e.g. when <tex>$Q_2$</tex> is used, apart from -1 the following
    !    values are allowed:
    ! 0=constant prolongation
    ! 1=linear prolongation of a once refined mesh
    ! 2=quadratic interpolation
    integer                     :: iprolongationOrder = -1

    ! Prolongation variant for nonconforming elements.
    ! Allows to switch to a special-type prolongation. Whether or not this
    ! has an effect depends on the discretisation.
    ! = 0: Use default prolongation.
    ! Uniform discretisation with E030/E031/EM30/EM31:
    ! = 1: Use standard prolongation, equally weighted (1/2 from left, 1/2 from right),
    ! = 2: Use extended prolongation, equally weighted (1/2 from left, 1/2 from right),
    ! = 3: Use extended prolongation, weighted by element size (L2 projection),
    ! = 4: Use extended prolongation, weighted by element size of neighbour element
    ! To activate extended prolongation, set this to >= 2 after initialising the
    ! interlevel projection structure!
    ! Uniform discretisation with Q1:
    ! = 2: FEAST mirror boundary prolongation with zero on the boundary (experimentally)
    integer                     :: iprolVariant = 0

    ! Configuration parameter for extended prolongation of E030/EM30/E031/EM31
    ! element. Only valid if iprolVariant >= 2.
    ! Aspect-ratio indicator; controls switching to constant prolongation.
    ! <=1: switch depending on aspect ratio of current element (standard),
    !  =2: switch depending on aspect ratio of current element and
    !      neighbour element
    integer                     :: iprolARIndicatorEX3Y = 1

    ! Configuration parameter for extended prolongation of E030/EM30/E031/EM31
    ! element. Only valid if iprolVariant >= 2.
    ! Upper bound aspect ratio; for all elements with higher AR
    ! the prolongation is switched to constant prolongation .
    ! This is set to 20.0 by default according to the analysis in
    ! [Michael Koester, Robuste Mehrgitter-Krylowraum-Techniken fuer FEM-Verfahren,
    !  2004, Diploma-Theses, Chair of Mathematics, University of Dortmund]
    real(DP)                    :: dprolARboundEX3Y = 20.0_DP

    ! Element type that should be assumed for the restriction. Must fit
    ! to the element type of the discretisation concerning the DOF`s:
    ! If element "EM30" is used e.g., restriction for "EM30", "EM31",
    ! "E031" and "E030" is possible, although performance will suffer
    ! in using the "wrong" restriction. An error will be displayed if
    ! the DOF`s do not fit together at all, e.g. when using <tex>$Q_1$</tex>
    ! restriction when using <tex>$Q_2$</tex> for the discretisation.
    !
    ! A value of EL_UNDEFINED indicates that the "natural" restriction
    ! (i.e. that configured by the spatial discretisation) should be used.
    integer(I32)                :: ielementTypeRestriction = EL_UNDEFINED

    ! Order of the restriction to use.
    ! -1=use default restriction (e.g. linear for <tex>$Q_1$</tex>, quadratic for
    !    <tex>$Q_2$</tex>,...). Some elements allow to configure the type of the
    !    restriction, e.g. when <tex>$Q_2$</tex> is used, apart from -1 the following
    !    values are allowed:
    ! 0=constant restriction
    ! 1=linear restriction of a once refined mesh
    ! 2=quadratic restriction
    integer                     :: irestrictionOrder = -1

    ! Restriction variant for nonconforming elements.
    ! Allows to switch to a special-type restriction. Whether or not this
    ! has an effect depends on the discretisation.
    ! = 0: Use default restriction.
    ! Uniform discretisation with E030/E031/EM30/EM31:
    ! = 1: Use standard restriction, equally weighted (1/2 from left, 1/2 from right),
    ! = 2: Use extended restriction, equally weighted (1/2 from left, 1/2 from right),
    ! = 3: Use extended restriction, weighted by element size (L2 projection),
    ! = 4: Use extended restriction, weighted by element size of neighbour element
    ! To activate extended prolongation, set this to >= 2 after initialising the
    ! interlevel projection structure!
    ! Uniform discretisation with Q1:
    ! = 1: Element-wise FEAST mirror boundary restriction (experimentally)
    ! = 2: FEAST mirror boundary restriction with zero on the boundary (experimentally)
    integer                     :: irestVariant = 0

    ! Configuration parameter for extended restriction of E030/EM30/E031/EM31
    ! element. Only valid if irestVariant >= 2.
    ! Aspect-ratio indicator; controls switching to constant prolongation.
    ! <=1: switch depending on aspect ratio of current element (standard),
    !  =2: switch depending on aspect ratio of current element and
    !      neighbour element
    integer                     :: irestARIndicatorEX3Y = 1

    ! Configuration parameter for extended restriction of E030/EM30/E031/EM31
    ! element. Only valid if irestVariant >= 2.
    ! Upper bound aspect ratio; for all elements with higher AR
    ! the prolongation is switched to constant prolongation .
    ! This is set to 20.0 by default according to the analysis in
    ! [Michael Koester, Robuste Mehrgitter-Krylowraum-Techniken fuer FEM-Verfahren,
    !  2004, Diploma-Theses, Chair of Mathematics, University of Dortmund]
    real(DP)                    :: drestARboundEX3Y = 20.0_DP

    ! Element type that should be assumed for the interpolation of a
    ! solution to a lower level. Must fit to the element type of the
    ! discretisation concerning the DOF`s.
    ! An error will be displayed if the DOF`s do not fit together at
    ! all, e.g. when using <tex>$Q_1$</tex> interpolation when using <tex>$Q_2$</tex>
    ! for the discretisation.
    !
    ! A value of EL_UNDEFINED indicates that the "natural" interpolation
    ! (i.e. that configured by the spatial discretisation) should be used.
    integer(I32)                :: ielementTypeInterpolation = EL_UNDEFINED

    ! Order of the interpolation to use when interpolating a solution vector
    ! to a lower level.
    ! -1=use default interpolation (e.g. linear for <tex>$Q_1$</tex>, quadratic for
    !    <tex>$Q_2$</tex>,...). Some elements allow to configure the type of the
    !    interpolation, e.g. when <tex>$Q_2$</tex> is used, apart from -1 the following
    !    values are allowed:
    ! 0=constant interpolation
    ! 1=linear interpolation of a once refined mesh
    ! 2=quadratic interpolation
    integer                     :: iinterpolationOrder = -1

    ! -------------------------------------------------------------------------
    ! L2-Projection Structures
    ! -------------------------------------------------------------------------
    ! Mass matrix of the fine mesh spatial discretisation.
    type(t_matrixScalar)        :: rmatrixMass

    ! 2-Level-Mass matrix between the spatial discretisations.
    type(t_matrixScalar)        :: rmatrix2Lvl

    ! Maximum number of iterations
    integer                     :: nmaxL2Iter = 50

    ! Relaxation parameter for SOR solver
    real(DP)                    :: drelaxL2 = 1.2_DP

    ! Absolute tolerance for L2-projection
    real(DP)                    :: depsL2 = 1E-15_DP

    ! -------------------------------------------------------------------------
    ! Matrix-based Projection Structures
    ! -------------------------------------------------------------------------
    ! The prolongation matrix that is to be used
    type(t_matrixScalar)        :: rmatrixProl

    ! The restriction matrix that is to be used.
    ! May be undefined, in this case, the restriction matrix is automatically
    ! derived from the prolongation matrix.
    type(t_matrixScalar)        :: rmatrixRest

    ! The interpolation matrix that is to be used.
    type(t_matrixScalar)        :: rmatrixInterp

  end type

!</typeblock>

!<typeblock>

  ! Contains a list of t_interlevelProjectionScalar structures that
  ! describes the projection of block vectors between different levels
  ! in the discretisation.
  ! The structure contains a 2-dimensional array: There is one
  ! t_interlevelProjectionScalar structure for every equation and
  ! for every element group in the equation.
  !
  ! The default initialisation initialises this structure with the
  ! usual values for the whole grid transfer process. If there is special
  ! need for "reconfiguring" the grid transfer for a special element
  ! type in a special equation, the application can change the content.

  type t_interlevelProjectionBlock

    ! A list of t_interlevelProjectionScalar structures for every
    ! equation and every element group in the discretisation.
    ! DIMENSION(1..#FE-spaces, 1..#equations)
    type(t_interlevelProjectionScalar), dimension(:,:), pointer :: RscalarProjection => null()

  end type

!</typeblock>

!<typeblock>

  ! A type encapsuling a pointer to a discretisation structure. Private.
  type t_interlevelProjectionDiscrPtr

    private

    ! Pointer to a discretisation
    type(t_blockDiscretisation), pointer :: p_rdiscretisation

  end type

!</typeblock>

!<typeblock>

  ! This structure defines a hierarchy of projection structures
  ! for a hierarchy of refinement levels. It can be used to project
  ! FEM functions given as block vectors from one level to another.
  type t_interlevelProjectionHier

    ! Minimum level in the hierarchy. Corresponds to index 1 in
    ! p_Rprojection and p_RdiscrPointer.
    integer :: nlmin = 0

    ! Maximum level in the hierarchy
    integer :: nlmax = 0

    ! An array of interlevel projection structures for all the levels.
    ! Note 1: This array is allocated as 1:nlmax-nlmin+1.
    ! Note 2: The projection structure for the coarse level p_Rprojection(1)
    !   is not used and not needed!
    type(t_interlevelProjectionBlock), dimension(:), pointer :: p_Rprojection

    ! An array of discretisation structures for all the levels
    type(t_interlevelProjectionDiscrPtr), dimension(:), pointer :: p_RdiscrPointer

    ! Temporary memory for the projection
    type(t_vectorScalar) :: rtempVector

  end type

!</typeblock>


!</types>

contains

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_initProjectionDirect (rprojection,RspatialDiscr)

!<description>
  ! This subroutine initialises an t_interlevelProjectionBlock with default
  ! values for a given PDE system. This allows the interlevel-projection
  ! routines to calculate constants that are necessary during a level
  ! change. (This is used e.g. in the general prolongation/restriction
  ! where the prolongation/restriction matrix must be calculated).
  !
  ! The calculated information is saved to rprojection and released
  ! with mlprj_doneProjection.
!</description>

!<input>
  ! An array of discretisation structures. Each discretisation structure
  ! corresponds to one scalar equation. The projection structure will be
  ! prepared to according to this discretisation list.
  type(t_spatialDiscretisation), dimension(:), intent(in) :: RspatialDiscr
!</input>

!<output>
  ! A t_interlevelProjectionBlock structure that will be filled with data
  ! about the projection of all the equations described by RspatialDiscr.
  type(t_interlevelProjectionBlock), intent(out) :: rprojection
!</output>

!</subroutine>

    integer :: nFEspaces,nequations,i

    ! Get the max. number of FE-spaces and the number of equations
    nequations = size(RspatialDiscr)
    nFEspaces = 0
    do i=1,nequations
      nFEspaces = max(nFEspaces,spdiscr_getNelemGroups(RspatialDiscr(i)))
    end do

    ! Allocate spatial discretisation structures for all the subblocks
    allocate(rprojection%RscalarProjection(nFEspaces,nequations))

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_initProjectionDiscr (rprojection,rdiscretisation)

!<description>
  ! This subroutine initialises an t_interlevelProjectionBlock with default
  ! values for a given PDE system. This allows the interlevel-projection
  ! routines to calculate constants that are necessary during a level
  ! change. (This is used e.g. in the general prolongation/restriction
  ! where the prolongation/restriction matrix must be calculated).
  !
  ! The calculated information is saved to rprojection and released
  ! with mlprj_doneProjection.
  !
  ! The PDE system is specified by the given block discretisation structure.
!</description>

!<input>
  ! A block discretisation structure specifying the discretisation
  type(t_blockDiscretisation), intent(in) :: rdiscretisation
!</input>

!<output>
  ! A t_interlevelProjectionBlock structure that will be filled with data
  ! about the projection of all the equations described by RspatialDiscr.
  type(t_interlevelProjectionBlock), intent(out) :: rprojection
!</output>

!</subroutine>

    if (rdiscretisation%ncomponents .eq. 0) then
      call output_line ("No discretisation!", &
          OU_CLASS_ERROR,OU_MODE_STD,"mlprj_initProjectionDiscr")
      call sys_halt()
    end if

    ! Call the standard initialisation routine
    call mlprj_initProjectionDirect (rprojection,&
         rdiscretisation%RspatialDiscr)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_initProjectionVec (rprojection,rvector)

!<description>
  ! This subroutine initialises an t_interlevelProjectionBlock with default
  ! values for a given PDE system. This allows the interlevel-projection
  ! routines to calculate constants that are necessary during a level
  ! change. (This is used e.g. in the general prolongation/restriction
  ! where the prolongation/restriction matrix must be calculated).
  !
  ! The calculated information is saved to rprojection and released
  ! with mlprj_doneProjection.
  !
  ! The PDE system is specified by the p_rspatialDiscretisation structures
  ! saved in the given vector.
!</description>

!<input>
  ! A vector containing information about the spatial discretisation of
  ! the given PDE.
  type(t_vectorBlock), intent(in) :: rvector
!</input>

!<output>
  ! A t_interlevelProjectionBlock structure that will be filled with data
  ! about the projection of all the equations described by RspatialDiscr.
  type(t_interlevelProjectionBlock), intent(out) :: rprojection
!</output>

!</subroutine>

    ! local variables
    type(t_spatialDiscretisation), dimension(rvector%nblocks) :: Rdiscr
    integer :: i

    if (rvector%nblocks .eq. 0) then
      call output_line ("No discretisation!", &
          OU_CLASS_ERROR,OU_MODE_STD,"mlprj_initProjectionVec")
      call sys_halt()
    end if

    ! Set up an array of discretisation structures for all the equations
    do i=1,rvector%nblocks
      Rdiscr(i) = rvector%RvectorBlock(i)%p_rspatialDiscr
    end do

    ! Call the standard initialisation routine
    call mlprj_initProjectionDirect (rprojection,Rdiscr)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_initProjectionMat (rprojection,rmatrix)

!<description>
  ! This subroutine initialises an t_interlevelProjectionBlock with default
  ! values for a given PDE system. This allows the interlevel-projection
  ! routines to calculate constants that are necessary during a level
  ! change. (This is used e.g. in the general prolongation/restriction
  ! where the prolongation/restriction matrix must be calculated).
  !
  ! The calculated information is saved to rprojection and released
  ! with mlprj_doneProjection.
  !
  ! The PDE system is specified by the p_rspatialDiscretisation structures
  ! saved in the given matrix.
!</description>

!<input>
  ! A matrix containing information about the spatial discretisation of
  ! the given PDE.
  type(t_matrixBlock), intent(in) :: rmatrix
!</input>

!<output>
  ! A t_interlevelProjectionBlock structure that will be filled with data
  ! about the projection of all the equations described by RspatialDiscr.
  type(t_interlevelProjectionBlock), intent(out) :: rprojection
!</output>

!</subroutine>

    ! local variables;
    ! type(t_spatialDiscretisation), dimension(max(rmatrix%nblocksPerRow,1)) :: Rdiscr
    type(t_spatialDiscretisation), dimension(:), pointer :: p_Rdiscr
    integer :: i,j

    if (rmatrix%nblocksPerRow .le. 0) then
      call output_line ("No discretisation!", &
          OU_CLASS_ERROR,OU_MODE_STD,"mlprj_initProjectionMat")
      call sys_halt()
    end if

    allocate(p_Rdiscr(rmatrix%nblocksPerRow))

    ! Set up an array of discretisation structures for all the equations.
    ! In every "column" of the block matrix, search for the first existing
    ! matrix and use its properties for initialisation
    do i=1,rmatrix%nblocksPerRow
      do j=1,rmatrix%nblocksPerCol
        if (rmatrix%RmatrixBlock(j,i)%NEQ .ne. 0) then
          p_Rdiscr(i) = &
            rmatrix%RmatrixBlock(j,i)%p_rspatialDiscrTrial
          exit
        end if
      end do
    end do

    ! Call the standard initialisation routine
    call mlprj_initProjectionDirect (rprojection,p_Rdiscr)

    deallocate(p_Rdiscr)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_doneProjection (rprojection)

!<description>
  ! Cleans up a t_interlevelProjectionBlock structure. All dynamically allocated
  ! memory is released.
!</description>

!<inputoutput>
  ! The t_interlevelProjectionBlock structure which is to be cleaned up.
  type(t_interlevelProjectionBlock), intent(inout) :: rprojection
!</inputoutput>

!</subroutine>
  integer :: i,j
  type(t_interlevelProjectionScalar), pointer :: p_rprj

    ! Release allocated memory
    if (associated(rprojection%RscalarProjection)) then

      ! Go through all scalar projections
      do i = 1, ubound(rprojection%RscalarProjection,1)
        do j = 1, ubound(rprojection%RscalarProjection,2)

          p_rprj => rprojection%RscalarProjection(i,j)

          ! What projection type did we use?
          select case(p_rprj%iprojType)
          case (MLP_PROJ_TYPE_HARDCODED)
            ! Nothing to do here...

          case (MLP_PROJ_TYPE_L2_PROJ)
            ! Release the mass matrices  for L2-projection
            call lsyssc_releaseMatrix(p_rprj%rmatrixMass)
            call lsyssc_releaseMatrix(p_rprj%rmatrix2Lvl)

          case (MLP_PROJ_TYPE_MATRIX)
            ! Release the prolongation matrix
            call lsyssc_releaseMatrix(p_rprj%rmatrixProl)
            if(lsyssc_hasMatrixStructure(p_rprj%rmatrixInterp)) &
              call lsyssc_releaseMatrix(p_rprj%rmatrixInterp)
            if(lsyssc_hasMatrixStructure(p_rprj%rmatrixRest)) &
              call lsyssc_releaseMatrix(p_rprj%rmatrixRest)

          end select

        end do ! j
      end do ! i

      ! Deallocate the scalar projection array
      deallocate(rprojection%RscalarProjection)

    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_getProjectionStrategy (rprojection,&
      rdiscrCoarse, rdiscrFine, ielemGroup, ractProjection)

!<description>
  ! Internal subroutine. This creates a projection structure ractProjection
  ! from a template projection structure rprojection and the discretisation
  ! structures on the coarse and fine grid and checks compatibility.
  !
  ! An error is thrown if either the DOF`s on the coarse- and fine grid do not
  ! fit together (i.e. there is no compatible projection) or if they do not fit
  ! to the projection template (e.g. if ielementTypeProlongation=EL_Q1 when the
  ! actual discretisation structure uses itrialElement=EL_Q2).
  ! Standard values (e.g. ielementTypeProlongation=EL_UNDEFINED) are replaced
  ! by the actual values (e.g. ielementTypeProlongation=EL_Q2,...).
!</description>

!<input>
  ! The t_interlevelProjectionScalar structure which is used as a template
  type(t_interlevelProjectionScalar), intent(in) :: rprojection

  ! Discretisation on the coarse mesh.
  type(t_spatialDiscretisation), intent(in) :: rdiscrCoarse

  ! Discretisation on the fine mesh.
  type(t_spatialDiscretisation), intent(in) :: rdiscrFine
  
  ! Current element group
  integer, intent(in) :: ielemGroup
!</input>

!<output>
  ! The t_interlevelProjectionScalar structure which configures the actual
  ! grid transfer between rdiscrCoarse and rdiscrFine.
  type(t_interlevelProjectionScalar), intent(out) :: ractProjection
!</output>

!</subroutine>
    
    ! local variables
    integer(I32) :: celemCoarse, celemFine

    ! Coarse and fine element
    call spdiscr_getElemGroupInfo (rdiscrCoarse,ielemGroup,celemCoarse)
    call spdiscr_getElemGroupInfo (rdiscrFine,ielemGroup,celemFine)

    ! Check to see if the two discretisation structures fit together
    if (celemCoarse .ne. celemFine) then
      call output_line ("Elements on the coarse and fine grid incompatible!", &
          OU_CLASS_ERROR,OU_MODE_STD,"mlprj_getProjectionStrategy")
      call sys_halt()
    end if

    ! Copy the template to the actual projection structure.
    ractProjection = rprojection

    ! Is any element type to be replaced according to a discretisation?
    if (ractProjection%ielementTypeProlongation .eq. EL_UNDEFINED) then
      ractProjection%ielementTypeProlongation = celemCoarse
    end if

    if (ractProjection%ielementTypeRestriction .eq. EL_UNDEFINED) then
      ractProjection%ielementTypeRestriction = celemCoarse
    end if

    if (ractProjection%ielementTypeInterpolation .eq. EL_UNDEFINED) then
      ractProjection%ielementTypeInterpolation = celemCoarse
    end if

    ! Check to see if the discretisation structures fit to the projection structure
    if (ractProjection%ielementTypeProlongation .ne. celemCoarse) then
      ! Here, we can insert some additional code so that E030 is compatible eith EM30
      ! and so on...
      call output_line ("element group of the grid and interlevel projection incompatible!", &
          OU_CLASS_ERROR,OU_MODE_STD,"mlprj_getProjectionStrategy")
      call sys_halt()
    end if

    if (ractProjection%ielementTypeRestriction .ne. celemCoarse) then
      ! Here, we can insert some additional code so that E030 is compatible eith EM30
      ! and so on...
      call output_line ("element group of the grid and interlevel projection incompatible!", &
          OU_CLASS_ERROR,OU_MODE_STD,"mlprj_getProjectionStrategy")
      call sys_halt()
    end if

    if (ractProjection%ielementTypeInterpolation .ne. celemCoarse) then
      ! Here, we can insert some additional code so that E030 is compatible eith EM30
      ! and so on...
      call output_line ("element group of the grid and interlevel projection incompatible!", &
          OU_CLASS_ERROR,OU_MODE_STD,"mlprj_getProjectionStrategy")
      call sys_halt()
    end if

  end subroutine

  ! ***************************************************************************

!<function>

  integer function mlprj_getTempMemoryScalar (rprojectionScalar,rdiscrCoarse,rdiscrFine)

!<description>
  ! This function returns for a given projection-structure and given
  ! discretisation structures on the coarse and fine grid the amount of temporary
  ! memory that is needed to do prolongation, restriction or interpolation.
  ! A value of 0 of course indicates that no temporary memory is needed.
!</description>

!<input>
  ! The t_interlevelProjectionScalar structure which configures
  ! the projection between rdiscrCoarse and rdiscrFine
  ! for each of the element groups in rdiscrCoarse /rdiscrFine
  type(t_interlevelProjectionScalar), dimension(:), intent(in) :: RprojectionScalar

  ! The element group structure of the equation on the coarse grid.
  type(t_spatialDiscretisation), intent(in) :: rdiscrCoarse

  ! The element group structure of the equation on the fine grid.
  type(t_spatialDiscretisation), intent(in) :: rdiscrFine
!</input>

!<result>
  ! Amount of temporary memory needed for prolongation/restriction/interpolation
  ! between vectors corresponding to the given combination of spatial
  ! discretisations.
  ! =0 if no memory is necessary.
!</result>

!</function>

  ! local variables
  integer :: i, nmem, nmemtotal

    ! Reserve at least as much memory as DOFs on the fine grid.
    ! Necessary for resorting sorted vectors.
    nmemtotal = dof_igetNDofGlob(rdiscrFine)

    ! Now loop through all scalar projections
    do i = 1, ubound(RprojectionScalar,1)

      ! What type of projection do we perform here?
      select case(RprojectionScalar(i)%iprojType)
      case (MLP_PROJ_TYPE_HARDCODED)
        ! In this case, no additional memory is needed.
        nmem = 0

      case (MLP_PROJ_TYPE_L2_PROJ)
        ! In this case, we will definately need additional memory.
        ! The amount of memory we need is two vectors at the size of
        ! the fine mesh mass matrix, as we need them for the internal
        ! solving process.
        nmem = 2*RprojectionScalar(i)%rmatrixMass%NEQ

      case (MLP_PROJ_TYPE_MATRIX)
        ! No temporary memory necessary here.
        nmem = 0

      case default
        ! Unknown - so nothing is needed.
        nmem = 0

      end select

      ! And let us choose the maximum
      nmemtotal = max(nmemtotal, nmem)

    end do

    ! Basically, there is no additional memory needed.
    mlprj_getTempMemoryScalar = nmemtotal

  end function

  ! ***************************************************************************

!<function>

  integer function mlprj_getTempMemoryDirect (rprojection,RdiscrCoarse,RdiscrFine)

!<description>
  ! Returns the amount of temporary memory that is needed by the
  ! interlevel-routines to transfer vectors between two grids.
  ! RdiscrCoarse and RdiscrFine are arrays of discretisation structures.
  ! Each discretisation structure corresponds to one scalar equation.
!</description>

!<input>
  ! Projection structure that configures the grid transfer for all equations
  ! and all element groups in each equation.
  type(t_interlevelProjectionBlock), intent(in) :: rprojection

  ! List of disretisation structures for the equations on the Coarse grid
  type(t_spatialDiscretisation), dimension(:), intent(in) :: RdiscrCoarse

  ! List of disretisation structures for the equations on the Fine grid
  type(t_spatialDiscretisation), dimension(:), intent(in) :: RdiscrFine
!</input>

!<result>
  ! Length of a temporary vector that is necessary for transferring vectors
  ! on the coarse grid to the fine grid and vice versa.
  ! =0 if no memory is necessary.
!</result>

!</function>

  ! local variables
  integer :: i
  integer :: imemmax,imemact

    if (size(RdiscrCoarse) .ne. size(RdiscrFine)) then
      call output_line ("Coarse and fine grid incompatible!", &
          OU_CLASS_ERROR,OU_MODE_STD,"mlprj_getTempMemoryDirect")
      call sys_halt()
    end if

    ! How much memory do we need?
    ! As we perform the grid transfer equation-by-equation and element-
    ! distribution by element-distribution, we need the maximum of all
    ! mlprj_getTempMemoryScalar calls!
    !
    ! Loop through all blocks:

    imemmax = 0
    do i=1,size(RdiscrCoarse)
      ! How much memory needed for that equation?
      imemact = mlprj_getTempMemoryScalar (rprojection%rscalarProjection(:,i),&
                                           RdiscrCoarse(i), RdiscrFine(i))

      imemmax = max(imemmax,imemact)
    end do

    ! Now we know how mich we need:
    mlprj_getTempMemoryDirect = imemmax

  end function

  ! ***************************************************************************

!<function>

  integer function mlprj_getTempMemoryVec (rprojection,rvectorCoarse,rvectorFine)

!<description>
  ! Returns the amount of temporary memory that is needed by the
  ! interlevel-routines to transfer vectors between two grids.
  ! rvectorCoarse and rvectorFine define template vectors on the coarse
  ! and fine grid, respectively; the memory is calculated using
  ! the discretisation structures associated to these.
!</description>

!<input>
  ! Projection structure that configures the grid transfer for all equations
  ! and all element groups in each equation.
  type(t_interlevelProjectionBlock), intent(in) :: rprojection

  ! Coarse grid vector. Must have the discretisation structure of the
  ! coarse grid attached.
  type(t_vectorBlock), intent(in) :: rvectorCoarse

  ! Fine grid vector. Must have the discretisation structure of the
  ! Fine grid attached.
  type(t_vectorBlock), intent(in) :: rvectorFine
!</input>

!<result>
  ! Length of a temporary vector that is necessary for transferring vectors
  ! on the coarse grid to the fine grid and vice versa.
  ! =0 if no memory is necessary.
!</result>

!</function>

  ! local variables
  type(t_spatialDiscretisation), dimension(rvectorCoarse%nblocks) :: RdiscrCoarse
  type(t_spatialDiscretisation), dimension(rvectorFine%nblocks) :: RdiscrFine
  integer :: i

    if ((rvectorCoarse%nblocks .eq. 0) .or. (rvectorFine%nblocks .eq. 0)) then
      call output_line ("No discretisation!", &
          OU_CLASS_ERROR,OU_MODE_STD,"mlprj_getTempMemoryVec")
      call sys_halt()
    end if

    ! Set up an array of discretisation structures for all the equations
    do i=1,rvectorCoarse%nblocks
      RdiscrCoarse(i) = &
        rvectorCoarse%RvectorBlock(i)%p_rspatialDiscr
    end do

    do i=1,rvectorFine%nblocks
      RdiscrFine(i) = &
        rvectorFine%RvectorBlock(i)%p_rspatialDiscr
    end do

    ! Call the standard getTempMemory routine
    mlprj_getTempMemoryVec = &
      mlprj_getTempMemoryDirect (rprojection, RdiscrCoarse,RdiscrFine)

  end function

  ! ***************************************************************************

!<function>

  integer function mlprj_getTempMemoryMat (rprojection,rmatrixCoarse,rmatrixFine)

!<description>
  ! Returns the amount of temporary memory that is needed by the
  ! interlevel-routines to transfer vectors between two grids.
  ! rmatrixCoarse and rmatrixFine define template matrices on the coarse
  ! and fine grid, respectively; the memory is calculated using
  ! the discretisation structures associated to these.
!</description>

!<input>
  ! Projection structure that configures the grid transfer for all equations
  ! and all element groups in each equation.
  type(t_interlevelProjectionBlock), intent(in) :: rprojection

  ! Coarse grid matrix. Must have the discretisation structure of the
  ! coarse grid attached.
  type(t_matrixBlock), intent(in) :: rmatrixCoarse

  ! Fine grid matrix. Must have the discretisation structure of the
  ! Fine grid attached.
  type(t_matrixBlock), intent(in) :: rmatrixFine
!</input>

!<result>
  ! Length of a temporary vector that is necessary for transferring vectors
  ! on the coarse grid to the fine grid and vice versa.
  ! =0 if no memory is necessary.
!</result>

!</function>

  ! local variables;
  ! type(t_spatialDiscretisation), dimension(max(rmatrixCoarse%nblocksPerRow,1)) :: RdiscrCoarse
  ! type(t_spatialDiscretisation), dimension(max(rmatrixFine%nBlocksPerRow,1)) :: RdiscrFine
  type(t_spatialDiscretisation), dimension(:), pointer :: p_RdiscrCoarse
  type(t_spatialDiscretisation), dimension(:), pointer :: p_RdiscrFine
  integer :: i,j

    if ((rmatrixCoarse%nblocksPerRow .le. 0) .or. (rmatrixFine%nblocksPerRow .le. 0)) then
      call output_line ("No discretisation!", &
          OU_CLASS_ERROR,OU_MODE_STD,"mlprj_getTempMemoryMat")
      call sys_halt()
    end if

    allocate(p_RdiscrCoarse(rmatrixCoarse%nblocksPerRow))
    allocate(p_RdiscrFine(rmatrixFine%nblocksPerRow))

    ! Set up an array of discretisation structures for all the equations
    do i=1,rmatrixCoarse%nblocksPerRow
      do j=1,rmatrixCoarse%nblocksPerCol
        if (lsysbl_isSubmatrixPresent(rmatrixCoarse,j,i)) then
          if (.not. &
              associated(rmatrixCoarse%RmatrixBlock(j,i)%p_rspatialDiscrTrial)) then
            call output_line ("No discretisation structure in coarse matrix at " // &
                trim(sys_siL(i,4)) // "," // trim(sys_siL(j,4)), &
                OU_CLASS_ERROR,OU_MODE_STD,"mlprj_getTempMemoryMat")
            call sys_halt()
          end if
          p_RdiscrCoarse(i) = &
            rmatrixCoarse%RmatrixBlock(j,i)%p_rspatialDiscrTrial
          exit
        end if
      end do
    end do

    do i=1,rmatrixFine%nblocksPerRow
      do j=1,rmatrixFine%nblocksPerCol
        if (lsysbl_isSubmatrixPresent(rmatrixFine,j,i)) then
          if (.not. &
              associated(rmatrixFine%RmatrixBlock(j,i)%p_rspatialDiscrTrial)) then
            call output_line ("No discretisation structure in fine matrix at " // &
                trim(sys_siL(i,4)) // "," // trim(sys_siL(j,4)), &
                OU_CLASS_ERROR,OU_MODE_STD,"mlprj_getTempMemoryMat")
            call sys_halt()
          end if
          p_RdiscrFine(i) = &
            rmatrixFine%RmatrixBlock(j,i)%p_rspatialDiscrTrial
        end if
      end do
    end do

    ! Call the standard getTempMemory routine
    mlprj_getTempMemoryMat = &
      mlprj_getTempMemoryDirect (rprojection, p_RdiscrCoarse,p_RdiscrFine)

    deallocate(p_RdiscrCoarse)
    deallocate(p_RdiscrFine)

  end function

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_performProlongation (rprojection,rcoarseVector, &
                                        rfineVector,rtempVector,icomponent)

!<description>
  ! Performs a prolongation for a given block vector (i.e. a projection
  ! in the primal space where the solution lives). The vector
  ! rcoarseVector on a coarser grid is projected to the vector
  ! rfineVector on a finer grid.
  ! rprojection configures how the grid transfer is performed.
  ! This projection structure must be build corresponding to the spatial
  ! discretisation structures in rcoarseVector and rfineVector!
!</description>

!<input>
  ! The t_interlevelProjectionBlock structure that configures the grid transfer
  type(t_interlevelProjectionBlock), intent(in) :: rprojection

  ! Coarse grid vector
  type(t_vectorBlock), intent(inout) :: rcoarseVector

  ! OPTIONAL: Component of the vector, the prolongation should be applied to.
  ! If not present, the prolongation is applied to all components.
  integer, intent(in), optional :: icomponent
!</input>

!<inputoutput>
  ! Temporary vector. The vector must have the same data type as
  ! rfineVector and must be at least as long as indicated by the function
  ! mlprj_getTempMemory. If mlprj_getTempMemory was =0, the vector may
  ! be a dummy vector.
  ! The vector does not have to be connected to a discretisation structure
  ! or something similar; the content is undefined at entry and will be
  ! undefined when leaving this routine.
  type(t_vectorScalar), intent(inout) :: rtempVector
!</inputoutput>

!<output>
  ! Fine grid vector
  type(t_vectorBlock), intent(inout) :: rfineVector
!</output>

!</subroutine>

  ! local variables
  integer :: i,istart,iend
  type(t_spatialDiscretisation), pointer :: p_rdiscrCoarse,p_rdiscrFine
  type(t_triangulation), pointer :: p_rtriaCoarse,p_rtriaFine
  type(t_interlevelProjectionScalar) :: ractProjection
  logical :: bsorted

  ! Pointers into the triangulation
  integer, dimension(:,:), pointer :: p_IverticesAtElementCoarse,&
      p_IverticesAtEdgeCoarse, p_IverticesAtElementFine, p_IverticesAtFaceCoarse
  integer, dimension(:,:), pointer :: p_IneighboursAtElementCoarse,&
      p_IneighboursAtElementFine
  integer, dimension(:,:), pointer :: p_IedgesAtElementCoarse,&
      p_IfacesAtElementCoarse, p_IedgesAtElementFine, p_IfacesAtElementFine
  integer(I32), dimension(:), pointer :: p_ItwistIndexCoarse, &
      p_ItwistIndexFine
  real(DP), dimension(:,:), pointer :: p_DvertexCoordsCoarse
  real(DP), dimension(:), pointer   :: p_DelementAreaCoarse

  ! Data arrays
  real(DP), dimension(:), pointer :: p_DuCoarse, p_DuFine

    ! The vectors must be of data type DOUBLE - we do not support anything
    ! different at the moment...
    if (rcoarseVector%cdataType .ne. ST_DOUBLE) then
      call output_line ("Coarse grid vector has unsupported data type!", &
          OU_CLASS_ERROR,OU_MODE_STD,"mlprj_performProlongation")
      call sys_halt()
    end if

    if (rfineVector%cdataType .ne. ST_DOUBLE) then
      call output_line ("Fine grid vector has unsupported data type!", &
          OU_CLASS_ERROR,OU_MODE_STD,"mlprj_performProlongation")
      call sys_halt()
    end if

    ! Calls the correct projection routine for each block in the
    ! discretisation...

    istart = 1
    iend = rcoarseVector%nblocks
    if (present(icomponent)) then
      istart = icomponent
      iend = icomponent
    end if

    do i=istart,iend

      ! Skip this block if it is empty
      if (rcoarseVector%RvectorBlock(i)%NEQ .le. 0) cycle
      
      ! What type of prolongation do we use here?
      select case(rprojection%RscalarProjection(1,i)%iprojType)
      case (MLP_PROJ_TYPE_L2_PROJ)
        ! Call scalar L2-prolongation
        call mlprj_prolScalarL2(rprojection%RscalarProjection(1,i), &
          rcoarseVector%RvectorBlock(i), rfineVector%RvectorBlock(i),rtempVector)

        ! Continue with next block
        cycle

      case (MLP_PROJ_TYPE_MATRIX)
        ! Matrix-based prolongation
        call lsyssc_matVec(rprojection%RscalarProjection(1,i)%rmatrixProl,&
          rcoarseVector%RvectorBlock(i), rfineVector%RvectorBlock(i), &
          1.0_DP, 0.0_DP, .false.)

        ! Continue with next block
        cycle

      end select

      p_rdiscrCoarse => rcoarseVector%RvectorBlock(i)%p_rspatialDiscr
      p_rdiscrFine => rfineVector%RvectorBlock(i)%p_rspatialDiscr

      ! We need a discretisation:
      if ((.not. associated(p_rdiscrCoarse)) .or. &
          (.not. associated(p_rdiscrFine))) then
        call output_line ("Intergrid transfer: No discretisation!", &
            OU_CLASS_ERROR,OU_MODE_STD,"mlprj_performProlongation")
        call sys_halt()
      end if

      ! Currently, we support only uniform triangulations.
      if ((p_rdiscrCoarse%ccomplexity .ne. SPDISC_UNIFORM) .or. &
          (p_rdiscrCoarse%ccomplexity .ne. SPDISC_UNIFORM)) then
        call output_line ("Intergrid transfer supports currently only uniform discretisations!", &
            OU_CLASS_ERROR,OU_MODE_STD,"mlprj_performProlongation")
        call sys_halt()
      end if

      ! Get the pointers to the vectors
      call lsyssc_getbase_double (rcoarseVector%RvectorBlock(i),p_DuCoarse)
      call lsyssc_getbase_double (rfineVector%RvectorBlock(i),p_DuFine)

      ! Use the first projection structure as template and create
      ! the actual projection structure for our situation.
      ! Remember, in a uniform grid we only have one projection structure
      ! and one element group!
      call mlprj_getProjectionStrategy (rprojection%RscalarProjection(1,i), &
          p_rdiscrCoarse, p_rdiscrFine, 1, ractProjection)

      ! Is the vector sorted?
      bsorted = rcoarseVector%RvectorBlock(i)%bisSorted
      
      ! The following hard-coded variants only work for unsorted vectors.
      if (bsorted) then
        call lsyssc_sortVector (rcoarseVector%RvectorBlock(i),.false.,rtempVector)
        
        ! Fine vector is set to "unsorted". The entries are overwritten,
        ! so we do not need to change them
        rfineVector%RvectorBlock(i)%bisSorted = .false.
      end if

      ! Depending on the element type of the trial functions in the
      ! discretisation, choose the right prolongation and call it.
      p_rtriaCoarse => p_rdiscrCoarse%p_rtriangulation
      p_rtriaFine => p_rdiscrFine%p_rtriangulation
      select case (elem_getPrimaryElement(ractProjection%ielementTypeProlongation))
      case (EL_P0_1D)
        ! P0 prolongation
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call mlprj_prolUniformP0_1D_double (p_DuCoarse,p_DuFine, &
             p_IneighboursAtElementFine, p_rtriaCoarse%NEL)

      case (EL_P1_1D)
        ! P1 prolongation
        call storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                             p_IverticesAtElementCoarse)
        call storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                             p_IverticesAtElementFine)
        call mlprj_prolUniformP1_1D_double (p_DuCoarse,p_DuFine, &
             p_IverticesAtElementCoarse,p_IverticesAtElementFine,&
             p_rtriaCoarse%NEL)

      case (EL_P2_1D)
        ! P2 prolongation
        call storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                             p_IverticesAtElementCoarse)
        call storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                             p_IverticesAtElementFine)
        call mlprj_prolUniformP2_1D_double (p_DuCoarse,p_DuFine, &
             p_IverticesAtElementCoarse, p_rtriaCoarse%NVT, p_rtriaCoarse%NEL)

      case (EL_S31_1D)
        ! S31 prolongation
        call storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                             p_IverticesAtElementCoarse)
        call storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                             p_IverticesAtElementFine)
        call storage_getbase_double2d(p_rtriaCoarse%h_DvertexCoords, &
                             p_DvertexCoordsCoarse)
        call mlprj_prolUniformS31_1D_double (p_DuCoarse,p_DuFine, &
             p_IverticesAtElementCoarse,p_IverticesAtElementFine,&
             p_DvertexCoordsCoarse, p_rtriaCoarse%NVT, &
             p_rtriaFine%NVT, p_rtriaCoarse%NEL)

      case (EL_P0)
        ! P0 prolongation
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call mlprj_prolUniformP0_double (p_DuCoarse,p_DuFine, &
             p_IneighboursAtElementFine, p_rtriaCoarse%NEL)

      case (EL_P1)
        ! P1 prolongation
        call storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                             p_IverticesAtElementCoarse)
        call storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                             p_IverticesAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                             p_IneighboursAtElementCoarse)
        call mlprj_prolUniformP1_double (p_DuCoarse,p_DuFine, &
             p_IverticesAtElementCoarse,p_IverticesAtElementFine,&
             p_IneighboursAtElementCoarse,p_rtriaCoarse%NEL)

      case (EL_P2)
        ! P2 prolongation
        call storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                             p_IverticesAtElementCoarse)
        call storage_getbase_int2d(p_rtriaCoarse%h_IedgesAtElement, &
                             p_IedgesAtElementCoarse)
        call storage_getbase_int2d(p_rtriaFine%h_IedgesAtElement, &
                             p_IedgesAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                             p_IneighboursAtElementCoarse)
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call mlprj_prolUniformP2_double (p_DuCoarse,p_DuFine, &
             p_IverticesAtElementCoarse,&
             p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
             p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
             p_rtriaCoarse%NEL,p_rtriaCoarse%NVT,p_rtriaFine%NVT)

      case (EL_P1T)
        ! P1~ prolongation, DOF`s = edge midpoint based
        call storage_getbase_int2d(p_rtriaFine%h_IedgesAtElement, &
                             p_IedgesAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IedgesAtElement, &
                             p_IedgesAtElementCoarse)
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                             p_IneighboursAtElementCoarse)
        call mlprj_prolUniformEx20_double (p_DuCoarse,p_DuFine, &
                p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
                p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
                p_rtriaCoarse%NEL)

      case (EL_Q0,EL_DG_T0_2d)
        ! Q0 prolongation
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call mlprj_prolUniformQ0_double (p_DuCoarse,p_DuFine, &
             p_IneighboursAtElementFine, p_rtriaCoarse%NEL)

      case (EL_Q1)
        ! Q1 prolongation
        select case (ractProjection%iprolVariant)
        case (:1) ! Standard prolongation
          call storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtEdge, &
                              p_IverticesAtEdgeCoarse)
          call storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                              p_IverticesAtElementCoarse)
          call storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                              p_IverticesAtElementFine)
          call mlprj_prolUniformQ1_double (p_DuCoarse, p_DuFine, &
                    p_IverticesAtEdgeCoarse, p_IverticesAtElementCoarse, &
                    p_rtriaCoarse%NVT, p_rtriaCoarse%NMT, p_rtriaCoarse%NEL)
           ! "old" implementation
!            CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
!                                p_IverticesAtElementCoarse)
!            CALL storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
!                                p_IverticesAtElementFine)
!            CALL storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
!                                p_IneighboursAtElementCoarse)
!            CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
!                                p_IneighboursAtElementFine)
!            CALL mlprj_prolUniformQ1_double (p_DuCoarse,p_DuFine, &
!                p_IverticesAtElementCoarse,p_IverticesAtElementFine,&
!                p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,p_rtriaCoarse%NEL)

        case (2:) !Experimental FEAST MIRROR prolongation with zero boundary
          call storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                              p_IverticesAtElementCoarse)
          call storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                              p_IverticesAtElementFine)
          call storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                              p_IneighboursAtElementCoarse)
          call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                              p_IneighboursAtElementFine)
          call mlprj_prolUnifQ1FMzero_double (p_DuCoarse,p_DuFine, &
              p_IverticesAtElementCoarse,p_IverticesAtElementFine,&
              p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
              p_rtriaCoarse%NEL,p_rtriaFine%NEL)
        end select

      case (EL_Q2)
        ! Q2 prolongation
        call storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                             p_IverticesAtElementFine)
        call storage_getbase_int2d(p_rtriaFine%h_IedgesAtElement, &
                             p_IedgesAtElementFine)
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)

        call mlprj_prolUniformQ2_double (p_DuCoarse,p_DuFine, &
             p_IverticesAtElementFine,p_IedgesAtElementFine,&
             p_IneighboursAtElementFine,&
             p_rtriaFine%NVT, p_rtriaFine%NMT, p_rtriaCoarse%NEL)

      case (EL_QP1)
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call mlprj_prolUniformQP1_double (p_DuCoarse,p_DuFine, &
                             p_IneighboursAtElementFine,&
                             p_rtriaCoarse%NEL, p_rtriaFine%NEL)

      case (EL_Q1T)
        ! Q1~ prolongation
        call storage_getbase_int2d(p_rtriaFine%h_IedgesAtElement, &
                             p_IedgesAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IedgesAtElement, &
                             p_IedgesAtElementCoarse)
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                             p_IneighboursAtElementCoarse)

        ! Type of prolongation? Extended or not?
        select case (ractProjection%iprolVariant)
        case (:1) ! Standard prolongation
          if (iand(ractProjection%ielementTypeProlongation,int(2**16,I32)) .ne. 0) then
            ! DOF`s = integral mean values
            call mlprj_prolUniformEx30_double (p_DuCoarse,p_DuFine, &
                p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
                p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
                p_rtriaCoarse%NEL)
          else
            ! DOF`s = edge midpoint based
            call mlprj_prolUniformEx31_double (p_DuCoarse,p_DuFine, &
                p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
                p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
                p_rtriaCoarse%NEL)
          end if

        case (2:) ! Extended prolongation; modified weights, local switch
                  ! to constant prolongation
          call storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                                     p_IverticesAtElementCoarse)
          call storage_getbase_double2d(p_rtriaCoarse%h_DvertexCoords, &
                                        p_DvertexCoordsCoarse)
          call storage_getbase_double(p_rtriaCoarse%h_DelementVolume, &
                                      p_DelementAreaCoarse)
          ! (what a nasty call...)
          if (iand(ractProjection%ielementTypeProlongation,int(2**16,I32)) .ne. 0) then
            ! DOF`s = integral mean values
            call mlprj_prolUniformEx30ext_double (p_DuCoarse,p_DuFine, &
                    p_DvertexCoordsCoarse,p_IverticesAtElementCoarse, &
                    p_DelementAreaCoarse,&
                    p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
                    p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
                    p_rtriaCoarse%NEL, &
                    min(4,ractProjection%iprolVariant)-2, &
                    ractProjection%dprolARboundEX3Y, &
                    ractProjection%iprolARIndicatorEX3Y)
          else
            ! DOF`s = edge midpoint based
            call mlprj_prolUniformEx31ext_double (p_DuCoarse,p_DuFine, &
                    p_DvertexCoordsCoarse,p_IverticesAtElementCoarse, &
                    p_DelementAreaCoarse,&
                    p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
                    p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
                    p_rtriaCoarse%NEL, &
                    min(4,ractProjection%iprolVariant)-2, &
                    ractProjection%dprolARboundEX3Y, &
                    ractProjection%iprolARIndicatorEX3Y)
          end if
        end select

      case (EL_Q2T)
        ! Q2~ prolongation
        call storage_getbase_int2d(p_rtriaFine%h_IedgesAtElement, &
                             p_IedgesAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IedgesAtElement, &
                             p_IedgesAtElementCoarse)
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                             p_IneighboursAtElementCoarse)
        call storage_getbase_int32(p_rtriaFine%h_ItwistIndex, &
                             p_ItwistIndexFine)
        call storage_getbase_int32(p_rtriaCoarse%h_ItwistIndex, &
                             p_ItwistIndexCoarse)

        call mlprj_prolUniformEx50_double (p_DuCoarse,p_DuFine, &
            p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
            p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
            p_ItwistIndexCoarse,p_ItwistIndexFine,&
            p_rtriaCoarse%NMT,p_rtriaFine%NMT,p_rtriaCoarse%NEL)

      case (EL_Q2TB)
        ! Q2~ with bubble prolongation
        call storage_getbase_int2d(p_rtriaFine%h_IedgesAtElement, &
                             p_IedgesAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IedgesAtElement, &
                             p_IedgesAtElementCoarse)
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                             p_IneighboursAtElementCoarse)
        call storage_getbase_int32(p_rtriaFine%h_ItwistIndex, &
                             p_ItwistIndexFine)
        call storage_getbase_int32(p_rtriaCoarse%h_ItwistIndex, &
                             p_ItwistIndexCoarse)

        call mlprj_prolUniformEB50_double (p_DuCoarse,p_DuFine, &
            p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
            p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
            p_ItwistIndexCoarse,p_ItwistIndexFine,&
            p_rtriaCoarse%NMT,p_rtriaFine%NMT,&
            p_rtriaCoarse%NEL,p_rtriaFine%NEL)

      case (EL_Q0_3D)
        ! Q0 prolongation
        call mlprj_prolUniformQ0_3D_double (p_DuCoarse,p_DuFine, p_rtriaCoarse%NEL)

      case (EL_Q1_3D)
        ! Q1 prolongation
        call storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtEdge, &
                            p_IverticesAtEdgeCoarse)
        call storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtFace, &
                            p_IverticesAtFaceCoarse)
        call storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                            p_IverticesAtElementCoarse)
        call storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                            p_IverticesAtElementFine)
        call mlprj_prolUniformQ1_3D_double (p_DuCoarse, p_DuFine, &
                  p_IverticesAtEdgeCoarse, p_IverticesAtFaceCoarse,&
                  p_IverticesAtElementCoarse, p_rtriaCoarse%NVT, &
                  p_rtriaCoarse%NMT, p_rtriaCoarse%NAT, p_rtriaCoarse%NEL)

      case (EL_Q1T_3D)
        ! Q1~ prolongation, DOF`s = integral mean value
        call storage_getbase_int2d(p_rtriaFine%h_IfacesAtElement, &
                             p_IfacesAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IfacesAtElement, &
                             p_IfacesAtElementCoarse)
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                             p_IneighboursAtElementCoarse)

        ! Type of prolongation? Extended or not?
!          SELECT CASE (ractProjection%iprolVariant)
!          CASE (:1) ! Standard prolongation
!            IF (IAND(ractProjection%ielementTypeProlongation,INT(2**16,I32)) .NE. 0) THEN
!              ! DOF`s = integral mean values
!              CALL mlprj_prolUniformEx30_double (p_DuCoarse,p_DuFine, &
!                  p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
!                  p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
!                  p_rtriaCoarse%NVT,p_rtriaFine%NVT,p_rtriaCoarse%NEL)
!            ELSE
            ! DOF`s = face midpoint based
            call mlprj_prolUniformEx3x_3D_double (p_DuCoarse,p_DuFine, &
                p_IfacesAtElementCoarse,p_IfacesAtElementFine,&
                p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
                p_rtriaCoarse%NEL,ractProjection%ielementTypeProlongation)
!            END IF
!
!          CASE (2:) ! Extended prolongation; modified weights, local switch
!                    ! to constant prolongation
!            CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
!                                       p_IverticesAtElementCoarse)
!            CALL storage_getbase_double2d(p_rtriaCoarse%h_DvertexCoords, &
!                                          p_DvertexCoordsCoarse)
!            CALL storage_getbase_double(p_rtriaCoarse%h_DelementVolume, &
!                                        p_DelementAreaCoarse)
!            ! (what a nasty call...)
!            IF (IAND(ractProjection%ielementTypeProlongation,INT(2**16,I32)) .NE. 0) THEN
!              ! DOF`s = integral mean values
!              CALL mlprj_prolUniformEx30ext_double (p_DuCoarse,p_DuFine, &
!                      p_DvertexCoordsCoarse,p_IverticesAtElementCoarse, &
!                      p_DelementAreaCoarse,&
!                      p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
!                      p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
!                      p_rtriaCoarse%NVT,p_rtriaFine%NVT,p_rtriaCoarse%NEL, &
!                      MIN(4,ractProjection%iprolVariant)-2, &
!                      ractProjection%dprolARboundEX3Y, &
!                      ractProjection%iprolARIndicatorEX3Y)
!            ELSE
!              ! DOF`s = edge midpoint based
!              CALL mlprj_prolUniformEx31ext_double (p_DuCoarse,p_DuFine, &
!                      p_DvertexCoordsCoarse,p_IverticesAtElementCoarse, &
!                      p_DelementAreaCoarse,&
!                      p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
!                      p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
!                      p_rtriaCoarse%NVT,p_rtriaFine%NVT,p_rtriaCoarse%NEL, &
!                      MIN(4,ractProjection%iprolVariant)-2, &
!                      ractProjection%dprolARboundEX3Y, &
!                      ractProjection%iprolARIndicatorEX3Y)
!            END IF
!          END SELECT

      case (EL_DG_T1_2D)
        ! DG_T1 prolongation
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call mlprj_prolUniformDGT1_double (p_DuCoarse,p_DuFine, &
             p_IneighboursAtElementFine, p_rtriaCoarse%NEL)

      case (EL_DG_T2_2D)
        ! DG_T2 prolongation
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call mlprj_prolUniformDGT2_double (p_DuCoarse,p_DuFine, &
             p_IneighboursAtElementFine, p_rtriaCoarse%NEL)

      case (EL_DG_Q1_2D)
        ! DG_Q1 prolongation
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call mlprj_prolUniformDGQ1_double (p_DuCoarse,p_DuFine, &
             p_IneighboursAtElementFine, p_rtriaCoarse%NEL)

      case (EL_DG_Q2_2D)
        ! DG_Q2 prolongation
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call mlprj_prolUniformDGQ2_double (p_DuCoarse,p_DuFine, &
             p_IneighboursAtElementFine, p_rtriaCoarse%NEL)

      case default
        call output_line ("Unknown element!", &
            OU_CLASS_ERROR,OU_MODE_STD,"mlprj_performProlongation")
        call sys_halt()
      end select
      
      ! Switch back sorted vectors.
      if (bsorted) then
        call lsyssc_sortVector (rfineVector%RvectorBlock(i),.true.,rtempVector)
        call lsyssc_sortVector (rcoarseVector%RvectorBlock(i),.true.,rtempVector)
      end if

    end do  ! i

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_performRestriction (rprojection,rcoarseVector, &
                                       rfineVector,rtempVector,icomponent)

!<description>
  ! Performs a restriction for a given block vector (i.e. a projection
  ! in the dual space where the RHS vector lives). The vector
  ! rfineVector on a finer grid is projected to the vector
  ! rcoarseVector on a coarser grid.
  ! rprojection configures how the grid transfer is performed.
  ! This projection structure must be build corresponding to the spatial
  ! discretisation structures in rcoarseVector and rfineVector!
!</description>

!<input>
  ! The t_interlevelProjectionBlock structure that configures the grid transfer
  type(t_interlevelProjectionBlock), intent(in) :: rprojection

  ! Fine grid vector
  type(t_vectorBlock), intent(inout) :: rfineVector

  ! OPTIONAL: Component of the vector, the restriction should be applied to.
  ! If not present, the restriction is applied to all components.
  integer, intent(in), optional :: icomponent
!</input>

!<inputoutput>
  ! Temporary vector. The vector must have the same data type as
  ! rcoarseVector and must be at least as long as indicated by the function
  ! mlprj_getTempMemory. If mlprj_getTempMemory was =0, the vector may
  ! be a dummy vector.
  ! The vector does not have to be connected to a discretisation structure
  ! or something similar; the content is undefined at entry and will be
  ! undefined when leaving this routine.
  type(t_vectorScalar), intent(inout) :: rtempVector

  ! Coarse grid vector
  type(t_vectorBlock), intent(inout) :: rcoarseVector
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i,istart,iend
  logical :: bsorted
  type(t_spatialDiscretisation), pointer :: p_rdiscrCoarse,p_rdiscrFine
  type(t_interlevelProjectionScalar) :: ractProjection
  type(t_triangulation), pointer :: p_rtriaCoarse,p_rtriaFine

  ! Pointers into the triangulation
  integer, dimension(:,:), pointer :: p_IverticesAtElementCoarse,&
      p_IverticesAtElementFine,p_IverticesAtEdgeCoarse,p_IverticesAtFaceCoarse
  integer, dimension(:,:), pointer :: p_IneighboursAtElementCoarse,&
      p_IneighboursAtElementFine
  integer, dimension(:,:), pointer :: p_IedgesAtElementCoarse,&
      p_IedgesAtElementFine, p_IfacesAtElementCoarse, p_IfacesAtElementFine
  integer(I32), dimension(:), pointer :: p_ItwistIndexCoarse, &
      p_ItwistIndexFine
  real(DP), dimension(:,:), pointer :: p_DvertexCoordsCoarse
  real(DP), dimension(:), pointer   :: p_DelementAreaCoarse

  ! Data arrays
  real(DP), dimension(:), pointer :: p_DuCoarse, p_DuFine

    ! The vectors must be of data type DOUBLE - we do not support anything
    ! different at the moment...
    if (rcoarseVector%cdataType .ne. ST_DOUBLE) then
      call output_line ("Coarse grid vector has unsupported data type!", &
          OU_CLASS_ERROR,OU_MODE_STD,"mlprj_performRestriction")
      call sys_halt()
    end if

    if (rfineVector%cdataType .ne. ST_DOUBLE) then
      call output_line ("Fine grid vector has unsupported data type!", &
          OU_CLASS_ERROR,OU_MODE_STD,"mlprj_performRestriction")
      call sys_halt()
    end if

    ! Calls the correct projection routine for each block in the
    ! discretisation...

    istart = 1
    iend = rcoarseVector%nblocks
    if (present(icomponent)) then
      istart = icomponent
      iend = icomponent
    end if

    do i=istart,iend

      ! Skip this block if it is empty
      if (rcoarseVector%RvectorBlock(i)%NEQ .le. 0) cycle

      ! What type of restriction do we use here?
      select case(rprojection%RscalarProjection(1,i)%iprojType)
      case(MLP_PROJ_TYPE_L2_PROJ)

        ! Call scalar L2-restriction
        call mlprj_restScalarL2(rprojection%RscalarProjection(1,i), &
          rcoarseVector%RvectorBlock(i), rfineVector%RvectorBlock(i),rtempVector)

        ! Continue with next block
        cycle

      case(MLP_PROJ_TYPE_MATRIX)

        ! Matrix-based restriction
        if (lsyssc_hasMatrixStructure(rprojection%RscalarProjection(1,i)%rmatrixRest)) then
          ! Uwse the restriction matrix specified in the structure.
          call lsyssc_matVec(rprojection%RscalarProjection(1,i)%rmatrixRest,&
            rfineVector%RvectorBlock(i), rcoarseVector%RvectorBlock(i), &
            1.0_DP, 0.0_DP)
        else
          ! The restriction is the transposed f the prolongation.
          call lsyssc_matVec(rprojection%RscalarProjection(1,i)%rmatrixProl,&
            rfineVector%RvectorBlock(i), rcoarseVector%RvectorBlock(i), &
            1.0_DP, 0.0_DP, .true.)
        end if

        ! Continue with next block
        cycle

      end select

      p_rdiscrCoarse => rcoarseVector%RvectorBlock(i)%p_rspatialDiscr
      p_rdiscrFine => rfineVector%RvectorBlock(i)%p_rspatialDiscr

      ! We need a discretisation:
      if ((.not. associated(p_rdiscrCoarse)) .or. &
          (.not. associated(p_rdiscrFine))) then
        call output_line ("Intergrid transfer: No discretisation!", &
            OU_CLASS_ERROR,OU_MODE_STD,"mlprj_performRestriction")
        call sys_halt()
      end if

      ! Currently, we support only uniform triangulations.
      if ((p_rdiscrCoarse%ccomplexity .ne. SPDISC_UNIFORM) .or. &
          (p_rdiscrFine%ccomplexity .ne. SPDISC_UNIFORM)) then
        call output_line ("Intergrid transfer supports currently only uniform discretisations!", &
            OU_CLASS_ERROR,OU_MODE_STD,"mlprj_performRestriction")
        call sys_halt()
      end if

      ! Get the pointers to the vectors
      call lsyssc_getbase_double (rcoarseVector%RvectorBlock(i),p_DuCoarse)
      call lsyssc_getbase_double (rfineVector%RvectorBlock(i),p_DuFine)

      ! Use the first projection structure as template and create
      ! the actual projection structure for our situation.
      ! Remember, in a uniform grid we only have one projection structure
      ! and one element group!
      call mlprj_getProjectionStrategy (rprojection%RscalarProjection(1,i), &
          p_rdiscrCoarse,p_rdiscrFine,1,ractProjection)

      ! Is the vector sorted?
      bsorted = rfineVector%RvectorBlock(i)%bisSorted
      
      ! The following hard-coded variants only work for unsorted vectors.
      if (bsorted) then
        call lsyssc_sortVector (rfineVector%RvectorBlock(i),.false.,rtempVector)
        
        ! Coarse vector is set to "unsorted". The entries are overwritten,
        ! so we do not need to change them
        rcoarseVector%RvectorBlock(i)%bisSorted = .false.
      end if

      ! Depending on the element type of the trial functions in the
      ! discretisation, choose the right prolongation and call it.
      p_rtriaCoarse => p_rdiscrCoarse%p_rtriangulation
      p_rtriaFine => p_rdiscrFine%p_rtriangulation
      select case (elem_getPrimaryElement(ractProjection%ielementTypeRestriction))
      case (EL_P0_1D)
        ! P0 restriction
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call mlprj_restUniformP0_1D_double (p_DuCoarse,p_DuFine, &
             p_IneighboursAtElementFine, p_rtriaCoarse%NEL)

      case (EL_P1_1D)
        ! P1 restriction
        call storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                             p_IverticesAtElementCoarse)
        call storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                             p_IverticesAtElementFine)
        call mlprj_restUniformP1_1D_double (p_DuCoarse,p_DuFine, &
             p_IverticesAtElementCoarse,p_IverticesAtElementFine, &
             p_rtriaCoarse%NEL)

      case (EL_P2_1D)
        ! P2 restriction
        call storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                             p_IverticesAtElementCoarse)
        call storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                             p_IverticesAtElementFine)
        call mlprj_restUniformP2_1D_double (p_DuCoarse,p_DuFine, &
             p_IverticesAtElementCoarse,p_rtriaCoarse%NVT, p_rtriaCoarse%NEL)

      case (EL_S31_1D)
        ! S31 restriction
        call storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                             p_IverticesAtElementCoarse)
        call storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                             p_IverticesAtElementFine)
        call storage_getbase_double2d(p_rtriaCoarse%h_DvertexCoords, &
                             p_DvertexCoordsCoarse)
        call mlprj_restUniformS31_1D_double (p_DuCoarse,p_DuFine, &
             p_IverticesAtElementCoarse,p_IverticesAtElementFine, &
             p_DvertexCoordsCoarse,p_rtriaCoarse%NVT, &
             p_rtriaFine%NVT, p_rtriaCoarse%NEL)

      case (EL_P0)
        ! P0 restriction
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call mlprj_restUniformP0_double (p_DuCoarse,p_DuFine, &
             p_IneighboursAtElementFine, p_rtriaCoarse%NEL)

      case (EL_P1)
        ! P1 restriction
        call storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                             p_IverticesAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                             p_IneighboursAtElementCoarse)
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call mlprj_restUniformP1_double (p_DuCoarse,p_DuFine, &
             p_IverticesAtElementFine,p_IneighboursAtElementFine, &
             p_rtriaCoarse%NEL,p_rtriaFine%NEL)

      case (EL_P2)
        ! P2 restriction
        call storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                             p_IverticesAtElementCoarse)
        call storage_getbase_int2d(p_rtriaCoarse%h_IedgesAtElement, &
                             p_IedgesAtElementCoarse)
        call storage_getbase_int2d(p_rtriaFine%h_IedgesAtElement, &
                             p_IedgesAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                             p_IneighboursAtElementCoarse)
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call mlprj_restUniformP2_double (p_DuCoarse,p_DuFine, &
             p_IverticesAtElementCoarse,&
             p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
             p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
             p_rtriaCoarse%NEL,p_rtriaCoarse%NVT,p_rtriaFine%NVT)

      case (EL_P1T)
        ! P1~ restriction, DOF`s = edge midpoint based
        call storage_getbase_int2d(p_rtriaFine%h_IedgesAtElement, &
                             p_IedgesAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IedgesAtElement, &
                             p_IedgesAtElementCoarse)
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                             p_IneighboursAtElementCoarse)
        call mlprj_restUniformEx20_double (p_DuCoarse,p_DuFine, &
                p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
                p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
                p_rtriaCoarse%NEL)

      case (EL_Q0,EL_DG_T0_2D)
        ! Q0 restriction
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call mlprj_restUniformQ0_double (p_DuCoarse,p_DuFine, &
             p_IneighboursAtElementFine, p_rtriaCoarse%NEL)

      case (EL_Q1)
        ! Q1 restriction

        ! Type of restriction?
        select case (ractProjection%irestVariant)
        case (:0) ! Standard restriction
          call storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtEdge, &
                              p_IverticesAtEdgeCoarse)
          call storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                              p_IverticesAtElementCoarse)
          call storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                              p_IverticesAtElementFine)
          call mlprj_restUniformQ1_double (p_DuCoarse, p_DuFine, &
                    p_IverticesAtEdgeCoarse, p_IverticesAtElementCoarse, &
                    p_rtriaCoarse%NVT, p_rtriaCoarse%NMT, p_rtriaCoarse%NEL)
           ! "old" implementation
!            CALL storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
!                                p_IverticesAtElementFine)
!            CALL storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
!                                p_IneighboursAtElementFine)
!            CALL mlprj_restUniformQ1_double (p_DuCoarse,p_DuFine, &
!                p_IverticesAtElementFine,p_IneighboursAtElementFine,&
!                p_rtriaFine%NEL)
        case (1) ! All boundaries are FEAST mirror boundaries
          call storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                              p_IverticesAtElementFine)
          call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                              p_IneighboursAtElementFine)
          call mlprj_restUniformQ1FM_double (p_DuCoarse,p_DuFine, &
              p_IverticesAtElementFine,p_IneighboursAtElementFine,&
              p_rtriaFine%NEL)
        case (2:) ! All boundaries are FEAST mirror boundaries with zero boundary
          call storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                              p_IverticesAtElementFine)
          call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                              p_IneighboursAtElementFine)
          call storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                              p_IverticesAtElementCoarse)
          call storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                              p_IneighboursAtElementCoarse)
          call mlprj_restUnifQ1FMzero_double (p_DuCoarse,p_DuFine, &
              p_IverticesAtElementFine,p_IneighboursAtElementFine,&
              p_IverticesAtElementCoarse,p_IneighboursAtElementCoarse,&
              p_rtriaFine%NEL,p_rtriaCoarse%NEL)

        end select

      case (EL_Q2)
        ! Q2 restriction
        call storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                             p_IverticesAtElementCoarse)
        call storage_getbase_int2d(p_rtriaCoarse%h_IedgesAtElement, &
                             p_IedgesAtElementCoarse)
        call storage_getbase_int2d(p_rtriaFine%h_IedgesAtElement, &
                             p_IedgesAtElementFine)
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)

        call mlprj_restUniformQ2_double (p_DuCoarse,p_DuFine, &
             p_IverticesAtElementCoarse, p_IedgesAtElementCoarse, &
             p_IedgesAtElementFine, p_IneighboursAtElementFine,&
             p_rtriaCoarse%NVT,p_rtriaFine%NVT,p_rtriaCoarse%NMT,&
             p_rtriaFine%NMT,p_rtriaCoarse%NEL)

      case (EL_QP1)
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call mlprj_restUniformQP1_double (p_DuCoarse,p_DuFine, &
             p_IneighboursAtElementFine, p_rtriaCoarse%NEL,p_rtriaFine%NEL)

      case (EL_Q1T)
        ! Q1~ restriction
        call storage_getbase_int2d(p_rtriaFine%h_IedgesAtElement, &
                             p_IedgesAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IedgesAtElement, &
                             p_IedgesAtElementCoarse)
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                             p_IneighboursAtElementCoarse)

        ! Type of restriction? Extended or not?
        select case (ractProjection%irestVariant)
        case (:1) ! Standard prolongation
          if (iand(ractProjection%ielementTypeRestriction,int(2**16,I32)) .ne. 0) then
            ! DOF`s = integral mean values
            call mlprj_restUniformEx30_double (p_DuCoarse,p_DuFine, &
                p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
                p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
                p_rtriaCoarse%NEL)
          else
            ! DOF`s = edge midpoint based
            call mlprj_restUniformEx31_double (p_DuCoarse,p_DuFine, &
                p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
                p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
                p_rtriaCoarse%NEL)
          end if

        case (2:) ! Extended prolongation; modified weights, local switch
                  ! to constant prolongation
          call storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                                     p_IverticesAtElementCoarse)
          call storage_getbase_double2d(p_rtriaCoarse%h_DvertexCoords, &
                                        p_DvertexCoordsCoarse)
          call storage_getbase_double(p_rtriaCoarse%h_DelementVolume, &
                                      p_DelementAreaCoarse)
          ! (what a nasty call...)
          if (iand(ractProjection%ielementTypeRestriction,int(2**16,I32)) .ne. 0) then
            ! DOF`s = integral mean values
            call mlprj_restUniformEx30ext_double (p_DuCoarse,p_DuFine, &
                    p_DvertexCoordsCoarse,p_IverticesAtElementCoarse, &
                    p_DelementAreaCoarse,&
                    p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
                    p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
                    p_rtriaCoarse%NEL, &
                    min(4,ractProjection%irestVariant)-2, &
                    ractProjection%dprolARboundEX3Y, &
                    ractProjection%iprolARIndicatorEX3Y)
          else
            ! DOF`s = edge midpoint based
            call mlprj_restUniformEx31ext_double (p_DuCoarse,p_DuFine, &
                    p_DvertexCoordsCoarse,p_IverticesAtElementCoarse, &
                    p_DelementAreaCoarse,&
                    p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
                    p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
                    p_rtriaCoarse%NEL, &
                    min(4,ractProjection%irestVariant)-2, &
                    ractProjection%dprolARboundEX3Y, &
                    ractProjection%iprolARIndicatorEX3Y)
          end if
        end select

      case (EL_Q2T)
        ! Q2~ prolongation
        call storage_getbase_int2d(p_rtriaFine%h_IedgesAtElement, &
                             p_IedgesAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IedgesAtElement, &
                             p_IedgesAtElementCoarse)
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                             p_IneighboursAtElementCoarse)
        call storage_getbase_int32(p_rtriaFine%h_ItwistIndex, &
                             p_ItwistIndexFine)
        call storage_getbase_int32(p_rtriaCoarse%h_ItwistIndex, &
                             p_ItwistIndexCoarse)

        call mlprj_restUniformEx50_double (p_DuCoarse,p_DuFine, &
            p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
            p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
            p_ItwistIndexCoarse,p_ItwistIndexFine,&
            p_rtriaCoarse%NMT,p_rtriaFine%NMT,p_rtriaCoarse%NEL)

      case (EL_Q2TB)
        ! Q2~ with bubble prolongation
        call storage_getbase_int2d(p_rtriaFine%h_IedgesAtElement, &
                             p_IedgesAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IedgesAtElement, &
                             p_IedgesAtElementCoarse)
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                             p_IneighboursAtElementCoarse)
        call storage_getbase_int32(p_rtriaFine%h_ItwistIndex, &
                             p_ItwistIndexFine)
        call storage_getbase_int32(p_rtriaCoarse%h_ItwistIndex, &
                             p_ItwistIndexCoarse)

        call mlprj_restUniformEB50_double (p_DuCoarse,p_DuFine, &
            p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
            p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
            p_ItwistIndexCoarse,p_ItwistIndexFine,&
            p_rtriaCoarse%NMT,p_rtriaFine%NMT,&
            p_rtriaCoarse%NEL,p_rtriaFine%NEL)

      case (EL_Q0_3D)
        ! Q0 restriction
        call mlprj_restUniformQ0_3D_double (p_DuCoarse,p_DuFine, p_rtriaCoarse%NEL)

      case (EL_Q1_3D)
        ! Q1 prolongation
        call storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtEdge, &
                            p_IverticesAtEdgeCoarse)
        call storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtFace, &
                            p_IverticesAtFaceCoarse)
        call storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
                            p_IverticesAtElementCoarse)
        call storage_getbase_int2d(p_rtriaFine%h_IverticesAtElement, &
                            p_IverticesAtElementFine)
        call mlprj_restUniformQ1_3D_double (p_DuCoarse, p_DuFine, &
                  p_IverticesAtEdgeCoarse, p_IverticesAtFaceCoarse,&
                  p_IverticesAtElementCoarse, p_rtriaCoarse%NVT, &
                  p_rtriaCoarse%NMT, p_rtriaCoarse%NAT, p_rtriaCoarse%NEL)

      case (EL_Q1T_3D)
        ! Q1~ restriction
        call storage_getbase_int2d(p_rtriaFine%h_IfacesAtElement, &
                             p_IfacesAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IfacesAtElement, &
                             p_IfacesAtElementCoarse)
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                             p_IneighboursAtElementCoarse)

!          ! Type of restriction? Extended or not?
!          SELECT CASE (ractProjection%irestVariant)
!          CASE (:1) ! Standard prolongation
!            IF (IAND(ractProjection%ielementTypeRestriction,INT(2**16,I32)) .NE. 0) THEN
!              ! DOF`s = integral mean values
!              CALL mlprj_restUniformEx30_double (p_DuCoarse,p_DuFine, &
!                  p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
!                  p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
!                  p_rtriaCoarse%NVT,p_rtriaFine%NVT,p_rtriaCoarse%NEL)
!            ELSE
            ! DOF`s = face midpoint based
            call mlprj_restUniformEx3x_3D_double (p_DuCoarse,p_DuFine, &
                p_IfacesAtElementCoarse,p_IfacesAtElementFine,&
                p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
                p_rtriaCoarse%NEL,ractProjection%ielementTypeRestriction)
!            END IF
!
!          CASE (2:) ! Extended prolongation; modified weights, local switch
!                    ! to constant prolongation
!            CALL storage_getbase_int2d(p_rtriaCoarse%h_IverticesAtElement, &
!                                       p_IverticesAtElementCoarse)
!            CALL storage_getbase_double2d(p_rtriaCoarse%h_DvertexCoords, &
!                                          p_DvertexCoordsCoarse)
!            CALL storage_getbase_double(p_rtriaCoarse%h_DelementVolume, &
!                                        p_DelementAreaCoarse)
!            ! (what a nasty call...)
!            IF (IAND(ractProjection%ielementTypeRestriction,INT(2**16,I32)) .NE. 0) THEN
!              ! DOF`s = integral mean values
!              CALL mlprj_restUniformEx30ext_double (p_DuCoarse,p_DuFine, &
!                      p_DvertexCoordsCoarse,p_IverticesAtElementCoarse, &
!                      p_DelementAreaCoarse,&
!                      p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
!                      p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
!                      p_rtriaCoarse%NVT,p_rtriaFine%NVT,p_rtriaCoarse%NEL, &
!                      MIN(4,ractProjection%irestVariant)-2, &
!                      ractProjection%dprolARboundEX3Y, &
!                      ractProjection%iprolARIndicatorEX3Y)
!            ELSE
!              ! DOF`s = edge midpoint based
!              CALL mlprj_restUniformEx31ext_double (p_DuCoarse,p_DuFine, &
!                      p_DvertexCoordsCoarse,p_IverticesAtElementCoarse, &
!                      p_DelementAreaCoarse,&
!                      p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
!                      p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
!                      p_rtriaCoarse%NVT,p_rtriaFine%NVT,p_rtriaCoarse%NEL, &
!                      MIN(4,ractProjection%irestVariant)-2, &
!                      ractProjection%dprolARboundEX3Y, &
!                      ractProjection%iprolARIndicatorEX3Y)
!            END IF
!          END SELECT

      case (EL_DG_T1_2D)
        ! DG_T1 restriction
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call mlprj_restUniformDGT1_double (p_DuCoarse,p_DuFine, &
             p_IneighboursAtElementFine, p_rtriaCoarse%NEL)

      case (EL_DG_T2_2D)
        ! DG_T2 restriction
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call mlprj_restUniformDGT2_double (p_DuCoarse,p_DuFine, &
             p_IneighboursAtElementFine, p_rtriaCoarse%NEL)

      case (EL_DG_Q1_2D)
        ! DG_Q1 restriction
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call mlprj_restUniformDGQ1_double (p_DuCoarse,p_DuFine, &
             p_IneighboursAtElementFine, p_rtriaCoarse%NEL)

      case (EL_DG_Q2_2D)
        ! DG_Q2 restriction
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call mlprj_restUniformDGQ2_double (p_DuCoarse,p_DuFine, &
             p_IneighboursAtElementFine, p_rtriaCoarse%NEL)

      case default
        call output_line ("Unknown element!", &
            OU_CLASS_ERROR,OU_MODE_STD,"mlprj_performRestriction")
        call sys_halt()
      end select

      ! Switch back sorted vectors.
      if (bsorted) then
        call lsyssc_sortVector (rfineVector%RvectorBlock(i),.true.,rtempVector)
        call lsyssc_sortVector (rcoarseVector%RvectorBlock(i),.true.,rtempVector)
      end if

    end do  ! i

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_performInterpolation (rprojection,rcoarseVector, &
                                         rfineVector,rtempVector,icomponent)

!<description>
  ! Performs an interpolation for a given block vector (i.e. a projection
  ! in the primal space where the solution vector lives). The vector
  ! rfineVector on a finer grid is projected to the vector
  ! rcoarseVector on a coarser grid.
  ! rprojection configures how the grid transfer is performed.
  ! This projection structure must be build corresponding to the spatial
  ! discretisation structures in rcoarseVector and rfineVector!
!</description>

!<input>

  ! The t_interlevelProjectionBlock structure that configures the grid transfer
  type(t_interlevelProjectionBlock), intent(in) :: rprojection

  ! Fine grid vector
  type(t_vectorBlock), intent(inout) :: rfineVector

  ! OPTIONAL: Component of the vector, the interpolation should be applied to.
  ! If not present, the interpolation is applied to all components.
  integer, intent(in), optional :: icomponent
!</input>

!<inputoutput>
  ! Temporary vector. The vector must have the same data type as
  ! rcoarseVector and must be at least as long as indicated by the function
  ! mlprj_getTempMemory. If mlprj_getTempMemory was =0, the vector may
  ! be a dummy vector.
  ! The vector does not have to be connected to a discretisation structure
  ! or something similar; the content is undefined at entry and will be
  ! undefined when leaving this routine.
  type(t_vectorScalar), intent(inout) :: rtempVector

  ! Coarse grid vector
  type(t_vectorBlock), intent(inout) :: rcoarseVector
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i,istart,iend
  logical :: bsorted
  type(t_spatialDiscretisation), pointer :: p_rdiscrCoarse,p_rdiscrFine
  type(t_interlevelProjectionScalar) :: ractProjection
  type(t_triangulation), pointer :: p_rtriaCoarse,p_rtriaFine

  ! Pointers into the triangulation
  integer, dimension(:,:), pointer :: p_IneighboursAtElementCoarse
  integer, dimension(:,:), pointer :: p_IneighboursAtElementFine
  integer, dimension(:,:), pointer :: p_IedgesAtElementCoarse
  integer, dimension(:,:), pointer :: p_IedgesAtElementFine
  integer, dimension(:,:), pointer :: p_IfacesAtElementCoarse
  integer, dimension(:,:), pointer :: p_IfacesAtElementFine
  integer(I32), dimension(:), pointer :: p_ItwistIndexCoarse, &
      p_ItwistIndexFine

  ! Data arrays
  real(DP), dimension(:), pointer :: p_DuCoarse, p_DuFine

    ! The vectors must be of data type DOUBLE - we do not support anything
    ! different at the moment...
    if (rcoarseVector%cdataType .ne. ST_DOUBLE) then
      call output_line ("Coarse grid vector has unsupported data type!", &
          OU_CLASS_ERROR,OU_MODE_STD,"mlprj_performInterpolation")
      call sys_halt()
    end if

    if (rfineVector%cdataType .ne. ST_DOUBLE) then
      call output_line ("Fine grid vector has unsupported data type!", &
          OU_CLASS_ERROR,OU_MODE_STD,"mlprj_performInterpolation")
      call sys_halt()
    end if

    ! Calls the correct projection routine for each block in the
    ! discretisation...

    istart = 1
    iend = rcoarseVector%nblocks
    if (present(icomponent)) then
      istart = icomponent
      iend = icomponent
    end if

    do i=istart,iend

      ! Skip this block if it is empty
      if (rcoarseVector%RvectorBlock(i)%NEQ .le. 0) cycle

      ! What type of interpolation do we use here?
      select case(rprojection%RscalarProjection(1,i)%iprojType)
      case(MLP_PROJ_TYPE_L2_PROJ)
        ! NOT SUPPORTED!
        call output_line ("L2-projection not supported for interpolation!", &
            OU_CLASS_ERROR,OU_MODE_STD,"mlprj_performInterpolation")
        call sys_halt()

      case(MLP_PROJ_TYPE_MATRIX)

        ! Matrix-based interpolation
        ! Remember that the interpolation matrix is stored transposed!
        call lsyssc_matVec(rprojection%RscalarProjection(1,i)%rmatrixInterp,&
          rfineVector%RvectorBlock(i), rcoarseVector%RvectorBlock(i), &
          1.0_DP, 0.0_DP, .true.)

        ! Continue with next block
        cycle

      end select


      p_rdiscrCoarse => rcoarseVector%RvectorBlock(i)%p_rspatialDiscr
      p_rdiscrFine => rfineVector%RvectorBlock(i)%p_rspatialDiscr

      ! We need a discretisation:
      if ((.not. associated(p_rdiscrCoarse)) .or. &
          (.not. associated(p_rdiscrFine))) then
        call output_line ("Intergrid transfer: No discretisation!", &
          OU_CLASS_ERROR,OU_MODE_STD,"mlprj_performInterpolation")
        call sys_halt()
      end if

      ! Currently, we support only uniform triangulations.
      if ((p_rdiscrCoarse%ccomplexity .ne. SPDISC_UNIFORM) .or. &
          (p_rdiscrCoarse%ccomplexity .ne. SPDISC_UNIFORM)) then
        call output_line ("Intergrid transfer supports currently only uniform discretisations!", &
          OU_CLASS_ERROR,OU_MODE_STD,"mlprj_performInterpolation")
        call sys_halt()
      end if

      ! Get the pointers to the vectors
      call lsyssc_getbase_double (rcoarseVector%RvectorBlock(i),p_DuCoarse)
      call lsyssc_getbase_double (rfineVector%RvectorBlock(i),p_DuFine)

      ! Use the first projection structure as template and create
      ! the actual projection structure for our situation.
      ! Remember, in a uniform grid we only have one projection structure
      ! and one element group!
      call mlprj_getProjectionStrategy (rprojection%RscalarProjection(1,i), &
          p_rdiscrCoarse,p_rdiscrFine,1,ractProjection)

      ! Is the vector sorted?
      bsorted = rfineVector%RvectorBlock(i)%bisSorted
      
      ! The following hard-coded variants only work for unsorted vectors.
      if (bsorted) then
        call lsyssc_sortVector (rfineVector%RvectorBlock(i),.false.,rtempVector)
        
        ! Coarse vector is set to "unsorted". The entries are overwritten,
        ! so we do not need to change them
        rcoarseVector%RvectorBlock(i)%bisSorted = .false.
      end if

      ! Depending on the element type of the trial functions in the
      ! discretisation, choose the right prolongation and call it.
      p_rtriaCoarse => p_rdiscrCoarse%p_rtriangulation
      p_rtriaFine => p_rdiscrFine%p_rtriangulation
      select case (elem_getPrimaryElement(ractProjection%ielementTypeProlongation))
      case (EL_P0_1D)
        ! P0 interpolation
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call mlprj_interpUniformP0_1D_double (p_DuCoarse,p_DuFine, &
             p_IneighboursAtElementFine, p_rtriaCoarse%NEL)

      case (EL_P1_1D)
        ! P1 interpolation
        call mlprj_interpUniformP1_1D_double (p_DuCoarse,p_DuFine,p_rtriaCoarse%NVT)

      case (EL_P2_1D)
        ! P2 interpolation
        call mlprj_interpUniformP2_1D_double (p_DuCoarse,p_DuFine,&
             p_rtriaCoarse%NVT, p_rtriaCoarse%NEL)

      case (EL_S31_1D)
        ! S31 interpolation
        call mlprj_interpUniS31_1D_double (p_DuCoarse,p_DuFine,&
             p_rtriaCoarse%NVT, p_rtriaFine%NVT)

      case (EL_P0)
        ! P0 interpolation
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call mlprj_interpUniformP0_double (p_DuCoarse,p_DuFine, &
            p_IneighboursAtElementFine, p_rtriaCoarse%NEL)

      case (EL_P1)
        ! P1 interpolation
        call mlprj_interpUniformP1_double (p_DuCoarse,p_DuFine,p_rtriaCoarse%NVT)

      case (EL_P2)
        ! P2 interpolation
        call mlprj_interpUniformP2_double (p_DuCoarse,p_DuFine, &
                                           p_rtriaCoarse%NVT, p_rtriaCoarse%NMT)

      case (EL_P1T)
        ! P1~ interpolation, DOF`s = edge midpoint based
        call storage_getbase_int2d(p_rtriaFine%h_IedgesAtElement, &
                             p_IedgesAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IedgesAtElement, &
                             p_IedgesAtElementCoarse)
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                             p_IneighboursAtElementCoarse)
        call mlprj_interpUniformEx20_double (p_DuCoarse,p_DuFine, &
             p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
             p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
             p_rtriaCoarse%NEL)

      case (EL_Q0,EL_DG_T0_2D)
        ! Q0 interpolation
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call mlprj_interpUniformQ0_double (p_DuCoarse,p_DuFine, &
             p_IneighboursAtElementFine, p_rtriaCoarse%NEL)

      case (EL_Q1)
        ! Q1 interpolation
        call mlprj_interpUniformQ1_double (p_DuCoarse,p_DuFine, p_rtriaCoarse%NVT)

      case (EL_Q2)
        ! Q2 interpolation
        call mlprj_interpUniformQ2_double (p_DuCoarse,p_DuFine, &
                   p_rtriaCoarse%NVT, p_rtriaCoarse%NMT, p_rtriaCoarse%NEL)

      case (EL_QP1)
        ! QP1 interpolation
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call mlprj_interpUniformQP1_double (p_DuCoarse,p_DuFine, &
                p_IneighboursAtElementFine, p_rtriaCoarse%NEL,p_rtriaFine%NEL)

      case (EL_Q1T)
        ! Q1~ interpolation, DOF`s = integral mean values
        ! We use the same routine also for interpolating Ex31 solutions - there is
        ! not too much difference...
        call storage_getbase_int2d(p_rtriaFine%h_IedgesAtElement, &
                             p_IedgesAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IedgesAtElement, &
                             p_IedgesAtElementCoarse)
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                             p_IneighboursAtElementCoarse)
        call mlprj_interpUniformEx30_double (p_DuCoarse,p_DuFine, &
             p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
             p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
             p_rtriaCoarse%NEL)

      case (EL_Q2T)
        ! Q2~ interpolation
        call storage_getbase_int2d(p_rtriaFine%h_IedgesAtElement, &
                             p_IedgesAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IedgesAtElement, &
                             p_IedgesAtElementCoarse)
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                             p_IneighboursAtElementCoarse)
        call storage_getbase_int32(p_rtriaFine%h_ItwistIndex, &
                             p_ItwistIndexFine)
        call storage_getbase_int32(p_rtriaCoarse%h_ItwistIndex, &
                             p_ItwistIndexCoarse)

        call mlprj_interpUniformEx50_double (p_DuCoarse,p_DuFine, &
            p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
            p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
            p_ItwistIndexCoarse,p_ItwistIndexFine,&
            p_rtriaCoarse%NMT,p_rtriaFine%NMT,p_rtriaCoarse%NEL)

      case (EL_Q2TB)
        ! Q2~ with bubble interpolation
        call storage_getbase_int2d(p_rtriaFine%h_IedgesAtElement, &
                             p_IedgesAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IedgesAtElement, &
                             p_IedgesAtElementCoarse)
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                             p_IneighboursAtElementCoarse)
        call storage_getbase_int32(p_rtriaFine%h_ItwistIndex, &
                             p_ItwistIndexFine)
        call storage_getbase_int32(p_rtriaCoarse%h_ItwistIndex, &
                             p_ItwistIndexCoarse)

        call mlprj_interpUniformEB50_double (p_DuCoarse,p_DuFine, &
            p_IedgesAtElementCoarse,p_IedgesAtElementFine,&
            p_IneighboursAtElementCoarse,p_IneighboursAtElementFine,&
            p_ItwistIndexCoarse,p_ItwistIndexFine,&
            p_rtriaCoarse%NMT,p_rtriaFine%NMT,&
            p_rtriaCoarse%NEL,p_rtriaFine%NEL)

      case (EL_Q0_3D)
        ! Q0 interpolation
        call mlprj_interpUniformQ0_3D_double (p_DuCoarse,p_DuFine,p_rtriaCoarse%NEL)

      case (EL_Q1_3D)
        ! Q1 interpolation
        call mlprj_interpUniformQ1_3D_double (p_DuCoarse,p_DuFine,p_rtriaCoarse%NVT)

      case (EL_Q1T_3D)
        ! Q1~ interpolation, DOF`s = integral mean values
        ! We use the same routine also for interpolating Ex31 solutions - there is
        ! not too much difference...
        call storage_getbase_int2d(p_rtriaFine%h_IfacesAtElement, &
                             p_IfacesAtElementFine)
        call storage_getbase_int2d(p_rtriaCoarse%h_IfacesAtElement, &
                             p_IfacesAtElementCoarse)
        call storage_getbase_int2d(p_rtriaCoarse%h_IneighboursAtElement, &
                             p_IneighboursAtElementCoarse)
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call mlprj_interpUniformEx3x_3D_dbl (p_DuCoarse,p_DuFine, &
                p_IfacesAtElementCoarse,p_IfacesAtElementFine,&
                p_IneighboursAtElementCoarse,p_IneighboursAtElementFine, &
                p_rtriaCoarse%NEL,ractProjection%ielementTypeInterpolation)

      case (EL_DG_T1_2D)
        ! DG_T1 interpolation
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call mlprj_interpUniformDGT1_double (p_DuCoarse,p_DuFine, &
             p_IneighboursAtElementFine, p_rtriaCoarse%NEL)

      case (EL_DG_T2_2D)
        ! DG_T2 interpolation
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call mlprj_interpUniformDGT2_double (p_DuCoarse,p_DuFine, &
             p_IneighboursAtElementFine, p_rtriaCoarse%NEL)

      case (EL_DG_Q1_2D)
        ! DG_Q1 interpolation
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call mlprj_interpUniformDGQ1_double (p_DuCoarse,p_DuFine, &
             p_IneighboursAtElementFine, p_rtriaCoarse%NEL)

      case (EL_DG_Q2_2D)
        ! DG_Q2 interpolation
        call storage_getbase_int2d(p_rtriaFine%h_IneighboursAtElement, &
                             p_IneighboursAtElementFine)
        call mlprj_interpUniformDGQ2_double (p_DuCoarse,p_DuFine, &
             p_IneighboursAtElementFine, p_rtriaCoarse%NEL)

      case default
        call output_line ("Unknown element!", &
          OU_CLASS_ERROR,OU_MODE_STD,"mlprj_performInterpolation")
        call sys_halt()
      end select

      ! Switch back sorted vectors.
      if (bsorted) then
        call lsyssc_sortVector (rfineVector%RvectorBlock(i),.true.,rtempVector)
        call lsyssc_sortVector (rcoarseVector%RvectorBlock(i),.true.,rtempVector)
      end if

    end do  ! i

  end subroutine

  ! ***************************************************************************
  ! Now the actual prolongation/restriction routines follow.
  ! For every type of element there is a set of routines that perform
  ! the actual grid transfer.
  ! We separate for
  ! - double and single precision vectors
  ! - different elements
  ! - uniform / conformal / ... discretisations
  ! ***************************************************************************

  ! ***************************************************************************
  ! Support for 1D P0 element
  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_prolUniformP0_1D_double (DuCoarse,DuFine, &
               IneighboursAtElementFine, NELcoarse)

!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! $P_0$, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Coarse grid vector
  real(DP), dimension(:), intent(in) :: DuCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Fine grid vector
  real(DP), dimension(:), intent(out) :: DuFine
!</output>

!</subroutine>

  ! local variables
  integer :: iel
  integer :: IELH1,IELH2
  real(DP) :: duh

    ! Loop over the elements
    do iel=1,NELCoarse

      ! Get the two child elements of the coarse grid element
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)

      ! Put the value on the coarse grid into all four child
      ! elements
      duh = DuCoarse(iel)
      DuFine(IELH1) = duh
      DuFine(IELH2) = duh
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_restUniformP0_1D_double (DuCoarse,DuFine, &
               IneighboursAtElementFine, NELcoarse)

!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! $P_0$, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  integer :: iel
  integer :: IELH1,IELH2

    ! The information that was "distributed" in the prolongation has to
    ! be "collected".

    ! Loop over the elements to collect the missing additive contributions:
    do iel=1,NELcoarse

      ! Get the elements on the fine grid that are children of the
      ! coarse grid element
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)

      ! Sum up the values in these nodes to get the
      ! value in the coarse grid element
      DuCoarse(iel)= DuFine(IELH1)+DuFine(IELH2)

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_interpUniformP0_1D_double (DuCoarse,DuFine, &
               IneighboursAtElementFine, NELcoarse)

!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! $P_0$, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  real(DP), parameter :: Q2 = 0.5_DP

  integer :: iel
  integer :: IELH1,IELH2

    ! The information that was "distributed" in the prolongation has to
    ! be "collected".

    ! Loop over the elements to collect the missing additive contributions:
    do iel=1,NELcoarse

      ! Get the elements on the fine grid that are children of the
      ! coarse grid element
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)

      ! From all fourtwo values, build the mean and use that as
      ! value of the coarse grid element
      DuCoarse(iel)= Q2*(DuFine(IELH1)+DuFine(IELH2))

    end do

  end subroutine

  ! ***************************************************************************
  ! Support for 1D P1 element
  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_prolUniformP1_1D_double (DuCoarse,DuFine,IvertsAtElemCoarse,&
                                            IvertsAtElemFine,NELcoarse)

!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! $P_1$, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Coarse grid vector
  real(DP), dimension(:), intent(in) :: DuCoarse

  ! IverticesAtElement array (KVERT) on the coarse grid
  integer, dimension(:,:), intent(in) :: IvertsAtElemCoarse

  ! IverticesAtElement array (KVERT) on the fine grid
  integer, dimension(:,:), intent(in) :: IvertsAtElemFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Fine grid vector
  real(DP), dimension(:), intent(out) :: DuFine
!</output>

!</subroutine>

  ! local variables
  real(DP), parameter :: Q2 = 0.5_DP
  integer :: iel
  real(DP) :: duh1,duh2

    ! Copy the first NVT entries - they belong to the coarse grid vertices
    ! that are fine grid vertices at the same time.
    call lalg_copyVector (DuCoarse,DuFine(1:size(DuCoarse)))

    ! Loop over the elements
    do iel=1,NELCoarse

      duh1=DuCoarse(IvertsAtElemCoarse(1,iel))
      duh2=DuCoarse(IvertsAtElemCoarse(2,iel))

      ! If IEL is the current element index of the coarse grid element,
      ! then it was refined into 2 new fine grid elements with indices
      ! IEL and NVT+IEL, where NVT is the number of coarse grid vertices.
      ! The "new" vertice in the fine grid corresponding to the coarse grid
      ! element IEL is the second vertice of the fine grid element IEL
      ! (and the first of the fine grid element NVT+IEL).
      DuFine(IvertsAtElemFine(2,iel)) = Q2 * (duh1 + duh2)

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_restUniformP1_1D_double (DuCoarse,DuFine,IvertsAtElemCoarse,&
               IvertsAtElemFine,NELcoarse)

!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! $P_1$, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IverticesAtElement array (KVERT) on the caorse grid
  integer, dimension(:,:), intent(in) :: IvertsAtElemCoarse

  ! IverticesAtElement array (KVERT) on the fine grid
  integer, dimension(:,:), intent(in) :: IvertsAtElemFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  real(DP), parameter :: Q2 = 0.5_DP
  integer :: iel
  integer :: ifgv, icgv1, icgv2

    ! The information that was "distributed" in the prolongation has to
    ! be "collected".
    !
    ! Copy the first NVT entries - this gives the first additive contribution.
    call lalg_copyVector (DuFine(1:size(DuCoarse)),DuCoarse)

    ! Loop over the elements to collect the missing additive contributions:
    do iel=1,NELcoarse

      ! Get the "new" fine grid vertice
      ifgv = IvertsAtElemFine(2,iel)

      ! Get the "old" coarse grid vertices
      icgv1 = IvertsAtElemCoarse(1,iel)
      icgv2 = IvertsAtElemCoarse(2,iel)

      DuCoarse(icgv1) = DuCoarse(icgv1) + Q2 * DuFine(ifgv)
      DuCoarse(icgv2) = DuCoarse(icgv2) + Q2 * DuFine(ifgv)

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_interpUniformP1_1D_double (DuCoarse,DuFine,NVTcoarse)

!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! $P_1$, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! Number of vertices in the coarse grid
  integer, intent(in) :: NVTcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables

    ! The first coase.NVT entries of the fine grid vector define
    ! the values on the coarse grid - because of the two-level ordering!
    call lalg_copyVector(DUfine(1:NVTcoarse),DUcoarse(1:NVTCoarse))

  end subroutine

  ! ***************************************************************************
  ! Support for 1D P2 element
  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_prolUniformP2_1D_double (DuCoarse,DuFine,IvertsAtElemCoarse,&
                                            NVTcoarse,NELcoarse)

!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! $P_1$, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Coarse grid vector
  real(DP), dimension(:), intent(in) :: DuCoarse

  ! IverticesAtElement array (KVERT) on the coarse grid
  integer, dimension(:,:), intent(in) :: IvertsAtElemCoarse

  ! Number of vertices in the coarse grid
  integer, intent(in) :: NVTcoarse

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Fine grid vector
  real(DP), dimension(:), intent(out) :: DuFine
!</output>

!</subroutine>

  ! local variables
  real(DP), parameter :: P_3_8 = 0.375_DP
  real(DP), parameter :: P_1_8 = 0.125_DP
  real(DP), parameter :: P_3_4 = 0.75_DP
  integer :: iel
  integer :: NVTfine
  real(DP) :: dv1,dv2,del

    ! Copy all DOFs from the coarse grid into the fine grid -
    ! the DOFs belonging to the "new" fine grid vertices get
    ! their values from the edge midpoints in the coarse grid.
    call lalg_copyVector (DuCoarse,DuFine(1:size(DuCoarse)))

    ! Calculate number of vertices in fine grid
    NVTfine = NVTcoarse+NELcoarse

    ! Loop over the elements
    do iel=1,NELCoarse

      ! Get the DOFs from the coarse grid vertices...
      dv1 = DuCoarse(IvertsAtElemCoarse(1,iel))
      dv2 = DuCoarse(IvertsAtElemCoarse(2,iel))

      ! ...and get the DOF from the coarse grid edge midpoint
      del = DuCoarse(NVTcoarse + iel)

      ! Perform quadratic interpolation to calculate the DOFs for the
      ! fine grid edge midpoints.
      DuFine(NVTfine + iel) = P_3_8*dv1 - P_1_8*dv2 + P_3_4*del
      DuFine(NVTfine + NELcoarse + iel) = -P_1_8*dv1 + P_3_8*dv2 + P_3_4*del

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_restUniformP2_1D_double (DuCoarse,DuFine,IvertsAtElemCoarse,&
               NVTcoarse,NELcoarse)

!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! $P_1$, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IverticesAtElement array (KVERT) on the caorse grid
  integer, dimension(:,:), intent(in) :: IvertsAtElemCoarse

  ! Number of vertices in the coarse grid
  integer, intent(in) :: NVTcoarse

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  real(DP), parameter :: P_3_8 = 0.375_DP
  real(DP), parameter :: P_1_8 = 0.125_DP
  real(DP), parameter :: P_3_4 = 0.75_DP
  integer :: iel
  integer :: NVTfine
  real(DP) :: dem1,dem2
  integer :: icgv1, icgv2,icgem

    ! The information that was "distributed" in the prolongation has to
    ! be "collected".
    !
    ! Copy the first NVT entries - this gives the first additive contribution.
    call lalg_copyVector (DuFine(1:size(DuCoarse)),DuCoarse)

    ! Calculate number of vertices in fine grid
    NVTfine = NVTcoarse+NELcoarse

    ! Loop over the elements to collect the missing additive contributions:
    do iel=1,NELcoarse

      ! Get the fine grid edge midpoints
      dem1 = DuFine(NVTfine+iel)
      dem2 = DuFine(NVTfine+NELcoarse+iel)

      icgv1 = IvertsAtElemCoarse(1,iel)
      icgv2 = IvertsAtElemCoarse(2,iel)
      icgem = NVTcoarse+iel

      ! distribute the DOFs
      DuCoarse(icgv1) = DuCoarse(icgv1) + P_3_8*dem1 - P_1_8*dem2
      DuCoarse(icgv2) = DuCoarse(icgv2) - P_1_8*dem1 + P_3_8*dem2
      DuCoarse(icgem) = DuCoarse(icgem) + P_3_4*(dem1 + dem2)

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_interpUniformP2_1D_double (DuCoarse,DuFine,NVTcoarse,NELcoarse)

!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! $P_1$, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! Number of vertices in the coarse grid
  integer, intent(in) :: NVTcoarse

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  integer :: len

    len = NVTcoarse+NELcoarse

    ! The first coase.NVT entries of the fine grid vector define
    ! the values on the coarse grid - because of the two-level ordering!
    call lalg_copyVector(DuFine(1:len),DuCoarse(1:len))

  end subroutine

  ! ***************************************************************************
  ! Support for 1D S31 element
  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_prolUniformS31_1D_double (DuCoarse,DuFine,IvertsAtElemCoarse,&
             IvertsAtElemFine,DvertexCoordsCoarse,NVTcoarse,NVTfine,NELcoarse)

!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! $S_31$, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Coarse grid vector
  real(DP), dimension(:), intent(in) :: DuCoarse

  ! IverticesAtElement array (KVERT) on the coarse grid
  integer, dimension(:,:), intent(in) :: IvertsAtElemCoarse

  ! IverticesAtElement array (KVERT) on the fine grid
  integer, dimension(:,:), intent(in) :: IvertsAtElemFine

  ! DvertexCoords array (DCORVG) on the coarse grid
  real(DP), dimension(:,:), intent(in) :: DvertexCoordsCoarse

  ! Number of vertices in the coarse grid
  integer, intent(in) :: NVTcoarse

  ! Number of vertices in the fine grid
  integer, intent(in) :: NVTfine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Fine grid vector
  real(DP), dimension(:), intent(out) :: DuFine
!</output>

!</subroutine>

  ! local variables
  real(DP), parameter :: Q12 = 0.5_DP
  real(DP), parameter :: Q14 = 0.25_DP
  real(DP), parameter :: Q34 = 0.75_DP
  integer :: iel
  integer :: ivt1, ivt2,ivt
  real(DP) :: dp1,dp2,dq1,dq2,ddet

    ! Copy the first NVTcoarse entries - these are the coefficients of the
    ! basis functions for the function values.
    call lalg_copyVector(DuCoarse(1:NVTcoarse),DuFine(1:NVTcoarse))

    ! Copy NVTcoarse entries beginning at NVTcoarse+1 of the coarse
    ! vector into the fine vector beginning at NVTfine+1 - these are the
    ! coefficients for the function derivative values.
    call lalg_copyVector(DuCoarse(NVTcoarse + 1 : 2*NVTcoarse), &
                         DuFine(NVTfine + 1 : NVTfine + NVTcoarse))

    ! Loop over the elements
    do iel=1,NELCoarse

      ! Vertices of element in coarse vector
      ivt1 = IvertsAtElemCoarse(1,iel)
      ivt2 = IvertsAtElemCoarse(2,iel)

      ! Calculate the determinant of this line
      ddet = 0.5 * (DvertexCoordsCoarse(1,ivt2) - DvertexCoordsCoarse(1,ivt1))

      ! Function value coefficients
      dp1=DuCoarse(ivt1)
      dp2=DuCoarse(ivt2)
      ! Function derivative coefficients
      dq1=DuCoarse(ivt1 + NVTcoarse)
      dq2=DuCoarse(ivt2 + NVTcoarse)

      ! Refined Vertice in fine vector
      ivt = IvertsAtElemFine(2,iel)

      ! Function Value
      DuFine(ivt) = Q12 * (dp1 + dp2) + Q14 * ddet * (dq1 - dq2)
      ! Function Derivative
      DuFine(ivt+NVTfine) = Q34 * (dp2 - dp1) / ddet - Q14 * (dq1 + dq2)

    end do


  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_restUniformS31_1D_double (DuCoarse,DuFine,IvertsAtElemCoarse,&
             IvertsAtElemFine,DvertexCoordsCoarse,NVTcoarse,NVTfine,NELcoarse)

!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! $P_1$, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IverticesAtElement array (KVERT) on the caorse grid
  integer, dimension(:,:), intent(in) :: IvertsAtElemCoarse

  ! IverticesAtElement array (KVERT) on the fine grid
  integer, dimension(:,:), intent(in) :: IvertsAtElemFine

  ! DvertexCoords array (DCORVG) on the coarse grid
  real(DP), dimension(:,:), intent(in) :: DvertexCoordsCoarse

  ! Number of vertices in the coarse grid
  integer, intent(in) :: NVTcoarse

  ! Number of vertices in the fine grid
  integer, intent(in) :: NVTfine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  real(DP), parameter :: Q12 = 0.5_DP
  real(DP), parameter :: Q14 = 0.25_DP
  real(DP), parameter :: Q34 = 0.75_DP
  integer :: iel
  integer :: ivtp1, ivtp2, ivtq1, ivtq2,ivt
  real(DP) :: dpf,dqf,ddet

    ! The information that was "distributed" in the prolongation has to
    ! be "collected".

    ! Copy the first NVTcoarse entries - these are the coefficients of the
    ! basis functions for the function values.
    call lalg_copyVector(DuFine(1:NVTcoarse),DuCoarse(1:NVTcoarse))

    ! Copy NVTcoarse entries beginning at NVTfine+1 of the fine vector into
    ! the coarse vector beginning at NVTfine+1 - these are the
    ! coefficients for the function derivative values.
    call lalg_copyVector(DuFine(NVTfine + 1 : NVTfine + NVTcoarse),&
                         DuCoarse(NVTcoarse + 1 : 2*NVTcoarse))

    ! Loop over the elements to collect the missing additive contributions:
    do iel=1,NELcoarse

      ! Vertices of element in coarse vector
      ivtp1 = IvertsAtElemCoarse(1,iel)
      ivtp2 = IvertsAtElemCoarse(2,iel)
      ivtq1 = ivtp1+NVTcoarse
      ivtq2 = ivtp2+NVTcoarse

      ! Calculate the determinant of this line
      ddet = 0.5 * (DvertexCoordsCoarse(1,ivtp2) - DvertexCoordsCoarse(1,ivtp1))

      ! Refined Vertice in fine vector
      ivt = IvertsAtElemFine(2,iel)
      dpf = DuFine(ivt)
      dqf = DuFine(ivt+NVTfine)

      ! Function values
      DuCoarse(ivtp1) = DuCoarse(ivtp1) + Q12*dpf - Q34*dqf/ddet
      DuCoarse(ivtp2) = DuCoarse(ivtp2) + Q12*dpf + Q34*dqf/ddet

      ! Function derivatives
      DuCoarse(ivtq1) = DuCoarse(ivtq1) + Q14*(dpf*ddet - dqf)
      DuCoarse(ivtq2) = DuCoarse(ivtq2) - Q14*(dpf*ddet + dqf)

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_interpUniS31_1D_double (DuCoarse,DuFine,NVTcoarse,NVTfine)

!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! $P_1$, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! Number of vertices in the coarse grid
  integer, intent(in) :: NVTcoarse

  ! Number of vertices in the fine grid
  integer, intent(in) :: NVTfine
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables

    ! The first NVTcoarse entries of the fine grid vector define
    ! the values on the coarse grid - because of the two-level ordering!
    call lalg_copyVector(DUfine(1:NVTcoarse), DUcoarse(1:NVTCoarse))
    call lalg_copyVector(DUfine(NVTfine + 1 : NVTfine + NVTcoarse),&
                         DUcoarse(NVTCoarse + 1 : 2*NVTCoarse))

  end subroutine

  ! ***************************************************************************
  ! Support for P0 element
  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_prolUniformP0_double (DuCoarse,DuFine, &
               IneighboursAtElementFine, NELcoarse)

!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! $P_0$, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Coarse grid vector
  real(DP), dimension(:), intent(in) :: DuCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Fine grid vector
  real(DP), dimension(:), intent(out) :: DuFine
!</output>

!</subroutine>

  ! local variables
  integer :: iel
  integer :: IELH1,IELH2,IELH3,IELH4
  real(DP) :: duh

    ! Loop over the elements
    do iel=1,NELCoarse

      ! Get the four child elements of the coarse grid element
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(1,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH1)
      IELH4 = IneighboursAtElementFine(3,IELH1)

      ! Put the value on the coarse grid into all four child
      ! elements
      duh = DuCoarse(iel)
      DuFine(IELH1) = duh
      DuFine(IELH2) = duh
      DuFine(IELH3) = duh
      DuFine(IELH4) = duh
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_restUniformP0_double (DuCoarse,DuFine, &
               IneighboursAtElementFine, NELcoarse)

!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! $P_0$, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  integer :: iel
  integer :: IELH1,IELH2,IELH3,IELH4

    ! The information that was "distributed" in the prolongation has to
    ! be "collected".

    ! Loop over the elements to collect the missing additive contributions:
    do iel=1,NELcoarse

      ! Get the elements on the fine grid that are children of the
      ! coarse grid element
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(1,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH1)
      IELH4 = IneighboursAtElementFine(3,IELH1)

      ! Sum up the values in these nodes to get the
      ! value in the coarse grid element
      DuCoarse(iel)= DuFine(IELH1)+DuFine(IELH2)+DuFine(IELH3)+DuFine(IELH4)

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_interpUniformP0_double (DuCoarse,DuFine, &
               IneighboursAtElementFine, NELcoarse)

!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! $P_0$, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  real(DP), parameter :: Q4 = 0.25_DP

  integer :: iel
  integer :: IELH1,IELH2,IELH3,IELH4

    ! The information that was "distributed" in the prolongation has to
    ! be "collected".

    ! Loop over the elements to collect the missing additive contributions:
    do iel=1,NELcoarse

      ! Get the elements on the fine grid that are children of the
      ! coarse grid element
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(1,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(3,IELH3)

      ! From all four values, build the mean and use that as
      ! value of the coarse grid element
      DuCoarse(iel)= Q4*(DuFine(IELH1)+DuFine(IELH2)+DuFine(IELH3)+DuFine(IELH4))

    end do

  end subroutine

  ! ***************************************************************************
  ! Support for P1 element
  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_prolUniformP1_double (DuCoarse,DuFine, &
               IverticesAtElementCoarse,IverticesAtElementFine,&
               IneighboursAtElementCoarse,NELcoarse)

!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! $P_1$, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Coarse grid vector
  real(DP), dimension(:), intent(in) :: DuCoarse

  ! IverticesAtElement array (KVERT) on the coarse grid
  integer, dimension(:,:), intent(in) :: IverticesAtElementCoarse

  ! IverticesAtElement array (KVERT) on the fine grid
  integer, dimension(:,:), intent(in) :: IverticesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementCoarse

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Fine grid vector
  real(DP), dimension(:), intent(out) :: DuFine
!</output>

!</subroutine>

  ! local variables
  real(DP), parameter :: Q2 = 0.5_DP
  integer :: iel
  real(DP) :: duh1,duh2,duh3

    ! Copy the first NVT entries - they belong to the coarse grid vertices
    ! that are fine grid vertices at the same time.
    call lalg_copyVector (DuCoarse,DuFine(1:size(DuCoarse)))

    ! Loop over the elements
    do iel=1,NELCoarse

      duh1=DuCoarse(IverticesAtElementCoarse(1,iel))
      duh2=DuCoarse(IverticesAtElementCoarse(2,iel))
      duh3=DuCoarse(IverticesAtElementCoarse(3,iel))

      ! Now check on every of the edges, if we already computed
      ! the value in the midpoint: Compute only if the neighbour
      ! element has smaller number.
      if (IneighboursAtElementCoarse(1,iel) .lt. iel) &
        DuFine(IverticesAtElementFine(1,iel)) = Q2*(duh1+duh2)

      if (IneighboursAtElementCoarse(2,iel) .lt. iel) &
        DuFine(IverticesAtElementFine(2,iel)) = Q2*(duh2+duh3)

      if (IneighboursAtElementCoarse(3,iel) .lt. iel) &
        DuFine(IverticesAtElementFine(3,iel)) = Q2*(duh3+duh1)

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_restUniformP1_double (DuCoarse,DuFine, &
               IverticesAtElementFine,IneighboursAtElementFine,&
               NELcoarse,NELfine)

!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! $P_1$, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IverticesAtElement array (KVERT) on the fine grid
  integer, dimension(:,:), intent(in) :: IverticesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse

  ! Number of elements in the fine grid
  integer, intent(in) :: NELfine
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  real(DP), parameter :: Q2 = 0.5_DP
  integer :: iel
  integer :: ih1

    ! The information that was "distributed" in the prolongation has to
    ! be "collected".
    !
    ! Copy the first NVT entries - this gives the first additive contribution.
    call lalg_copyVector (DuFine(1:size(DuCoarse)),DuCoarse)

    ! Loop over the elements to collect the missing additive contributions:
    do iel=NELcoarse+1,NELfine
      IH1 = IverticesAtElementFine(1,iel)

      ! Treat every edge only once:

      if (IneighboursAtElementFine(1,iel) .lt. iel) &
        DuCoarse(IH1) = DuCoarse(IH1) + Q2*DuFine(IverticesAtElementFine(2,iel))

      if (IneighboursAtElementFine(3,iel) .lt. iel) &
        DuCoarse(IH1) = DuCoarse(IH1) + Q2*DuFine(IverticesAtElementFine(3,iel))

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_interpUniformP1_double (DuCoarse,DuFine,NVTcoarse)

!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! $P_1$, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! Number of vertices in the coarse grid
  integer, intent(in) :: NVTcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

    ! The first coase.NVT entries of the fine grid vector define
    ! the values on the coarse grid - because of the two-level ordering!
    call lalg_copyVector(DUfine(1:NVTcoarse),DUcoarse(1:NVTCoarse))

  end subroutine

  ! ***************************************************************************
  ! Support for P2 element
  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_prolUniformP2_double (DuCoarse,DuFine, &
               IverticesAtElementCoarse,&
               IedgesAtElementCoarse,IedgesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               NELcoarse,NVTcoarse,NVTfine)

!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! $P_2$, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Coarse grid vector
  real(DP), dimension(:), intent(in) :: DuCoarse

  ! IverticesAtElement array (KVERT) on the coarse grid
  integer, dimension(:,:), intent(in) :: IverticesAtElementCoarse

  ! IedgesAtElement array (KMID) on the coarse grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementCoarse

  ! IedgesAtElement array (KMID) on the fine grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse

  ! Number of vertices in the coarse grid
  integer, intent(in) :: NVTcoarse

  ! Number of vertices in the fine grid
  integer, intent(in) :: NVTfine
!</input>

!<output>
  ! Fine grid vector
  real(DP), dimension(:), intent(out) :: DuFine
!</output>

!</subroutine>

  ! local variables
  integer :: IEL,IEL1,IEL2
  real(DP) :: duc1,duc2,duc3,dum1,dum2,dum3
  real(DP), parameter :: Q8 = 0.125_DP
  real(DP), parameter :: Q4 = 0.25_DP
  real(DP), parameter :: Q2 = 0.5_DP

    ! First, we remember the refinement scheme to clarify the
    ! local numbering of the vertices in the coarse- and fine-
    ! grid triangles.
    ! Let a coarse grid triangle be locally numbered as:
    !
    !   2
    !   |  \
    !   |    \
    !   | IEL  \
    !   |        \
    !   3----------1
    !
    ! Then the refinement process assigns the following numbers:
    !
    !   2
    !   |  \
    !   |     \
    !   |        \
    !   2*-------- 1*
    !   | \   IEL  |   \
    !   |    \     |      \
    !   |       \  |         \
    !   3----------3*----------1
    !
    ! i.e. the element number "IEL" is put into the middle with
    ! the corner vertices numbered according to the local edge
    ! numbers of the coarse grid element.
    !
    ! To access information on the edges of the fine grid element,
    ! we have to work with adjacencies!
    !
    ! First copy DuCoarse to DuFine. This will transfer the corner values
    ! from the coarse grid to the fine grid. More precisely, this
    ! will map:
    !   Coarse grid vertices 1..NVTC  -> Fine grid vertices 1..NVTC
    !   Coarse grid midpoints 1..NMTC -> Fine grid vertices NVTC+1..NVTC+NMTC = NVTF
    ! Afterwards, we only have to create the missing midpoint values!

    call lalg_copyVector (DuCoarse,DuFine(1:size(DuCoarse)))

    ! loop over the elements

    do IEL = 1,nelCoarse

      ! We fetch the function values of the coarse grid element
      ! into variables following the following scheme:
      !
      !     DUC2
      !     |   \
      !     |      \
      !     |         \
      !     DUM2         DUM1
      !     |                \
      !     |                   \
      !     |                      \
      !     DUC3 --------DUM3-------  DUC1

      DUC1 = DuCoarse(IverticesAtElementCoarse(1,IEL))
      DUC2 = DuCoarse(IverticesAtElementCoarse(2,IEL))
      DUC3 = DuCoarse(IverticesAtElementCoarse(3,IEL))
      DUM1 = DuCoarse(IedgesAtElementCoarse(1,IEL)+NVTcoarse)
      DUM2 = DuCoarse(IedgesAtElementCoarse(2,IEL)+NVTcoarse)
      DUM3 = DuCoarse(IedgesAtElementCoarse(3,IEL)+NVTcoarse)

      ! We have to calculate the function values in the new DOF`s on the
      ! fine grid, which are located here:
      !
      !     DUC2
      !     |   \
      !     X     X
      !     |        \
      !     DUM2 --X-- DUM1
      !     |   \      |   \
      !     X     X    X     X
      !     |       \  |        \
      !     DUC3---X---DUM3---X---DUC1
      !
      ! On this trangle, the function is a quadratic polynomial:
      !
      !   P(X,Y) = c1 + c2*x + c3*y + c4*x^2 + c5*x*y + c6*y^2
      !
      ! Solving for the coefficients such that the polynomial takes
      ! the DUxy-values in the corners/midpoints, we obtain the
      ! polynomial as:
      !
      !   P(X,Y) = DUC3 +
      !             (-3*DUC3-DUC1+4*DUM3)*x +
      !             (-3*DUC3-DUC2+4*DUM2)*y +
      !             (4*DUC3-4*DUM3+4*DUM1-4*DUM2)*x*y +
      !             (2*DUC3+2*DUC1-4*DUM3)*x^2 +
      !             (2*DUC3+2*DUC2-4*DUM2)*y^2
      !
      ! This has to be evaluated in the new points, marked as "X"
      ! in the above sketch.
      !
      ! Remember, that the coarse grid element IEL is moved
      ! to the inner fine-grid element. The corners of the coars
      ! grid element are always the first vertices of the triangles
      ! on the fine grid elements.
      !
      ! First calculate the edge mitpoint values of that one!
      !
      !   |        \
      !   DUM2---X---DUM1
      !   | \   IEL  |   \
      !         X     X
      !           \  |
      !           --DUM3--
      !
      ! DUF(IedgesAtElementFine(1,IEL)) = P(1/4,1/2)
      !                   = -1/8*DUC3-1/8*DUC1+1/4*DUM3+1/2*DUM2+1/2*DUM1
      ! DUF(IedgesAtElementFine(2,IEL)) = P(1/4,1/2)
      !                   = -1/8*DUC1+1/2*DUM3-1/8*DUC2+1/2*DUM2+1/4*DUM1
      ! DUF(IedgesAtElementFine(3,IEL)) = P(1/2,1/4)
      !                   = -1/8*DUC3+1/2*DUM3-1/8*DUC2+1/4*DUM2+1/2*DUM1

      DuFine(IedgesAtElementFine(1,IEL)+NVTfine) = -Q8*DUC3-Q8*DUC1+Q4*DUM3+Q2*DUM2+Q2*DUM1
      DuFine(IedgesAtElementFine(2,IEL)+NVTfine) = -Q8*DUC1+Q2*DUM3-Q8*DUC2+Q2*DUM2+Q4*DUM1
      DuFine(IedgesAtElementFine(3,IEL)+NVTfine) = -Q8*DUC3+Q2*DUM3-Q8*DUC2+Q4*DUM2+Q2*DUM1

      ! Now to the values on the other edges.
      ! Only calculate information on edges that are adjacent to
      ! an element with lower element number. This prevents
      ! information from being computed twice.
      !
      ! Check the first edge:

      if (IneighboursAtElementCoarse(1,IEL).lt.IEL) then

        ! Neighbor element has a smaller number than our element.
        ! Calculate the information on the current edge.
        ! This is the edge corresponding to the edge midpoint DUM1!
        ! We have to calculate the values in the new edge midpoints.
        !
        ! DUC2
        !  |   \
        !  |     X
        !  |  IEL1  \
        ! DUM2 ----- DUM1
        !      \ IEL  |   \
        !        \    |     X
        !          \  |  IEL2  \
        !            DUM3-------DUC1
        !
        ! Use adjacencies to get the fine-grid elements IEL1 and IEL2.

        IEL1 = IneighboursAtElementFine(1,IEL)
        IEL2 = IneighboursAtElementFine(3,IEL)

        ! Calculate the new edge midpoints:
        !
        !  DUF(IedgesAtElementFine(3,IEL1)) = P(1/4,3/4)
        !                     = -1/8*DUC1+3/8*DUC2+3/4*DUM1
        !  DUF(IedgesAtElementFine(1,IEL2)) = P(3/4,1/4)
        !                     = 3/8*DUC1-1/8*DUC2+3/4*DUM1

        DuFine(IedgesAtElementFine(3,IEL1)+NVTfine) = -Q8*DUC1+3.0_DP*Q8*DUC2+3.0_DP*Q4*DUM1
        DuFine(IedgesAtElementFine(1,IEL2)+NVTfine) = 3.0_DP*Q8*DUC1-Q8*DUC2+3.0_DP*Q4*DUM1

      end if

      ! Check the next edge:

      if (IneighboursAtElementCoarse(2,IEL).lt.IEL) then

        ! DUC2
        !  |   \
        !  X     \
        !  |   IEL1 \
        ! DUM2 ----- DUM1
        !  |   \      |
        !  X     \IEL |
        !  |  IEL2 \  |
        ! DUC3-------DUM3
        !
        ! Use adjacencies to get the fine-grid elements IEL1 and IEL2.

        IEL1 = IneighboursAtElementFine(1,IEL)
        IEL2 = IneighboursAtElementFine(2,IEL)

        ! Calculate the new edge midpoints:
        !
        !  DUF(IedgesAtElementFine(1,IEL1)) = P(0,3/4)
        !                     = -1/8*DUC3+3/8*DUC2+3/4*DUM2
        !  DUF(IedgesAtElementFine(3,IEL2)) = P(0,1/4)
        !                     = 3/8*DUC3-1/8*DUC2+3/4*DUM2

        DuFine(IedgesAtElementFine(1,IEL1)+NVTfine) = &
            -Q8*DUC3+3.0_DP*Q8*DUC2+3.0_DP*Q4*DUM2
        DuFine(IedgesAtElementFine(3,IEL2)+NVTfine) = &
            3.0_DP*Q8*DUC3-1.0_DP*Q8*DUC2+3.0_DP*Q4*DUM2

      end if

      ! Check the last edge

      if (IneighboursAtElementCoarse(3,IEL).lt.IEL) then

        ! DUM2 ----- DUM1
        !  |   \  IEL |   \
        !  | IEL2\    | IEL1\
        !  |       \  |        \
        ! DUC3---X---DUM3---X---DUC1
        !
        ! Use adjacencies to get the fine-grid elements IEL1 and IEL2.

        IEL1 = IneighboursAtElementFine(3,IEL)
        IEL2 = IneighboursAtElementFine(2,IEL)

        ! Calculate the new edge midpoints:
        !
        !  DUF(IedgesAtElementFine(3,IEL1)) = P(3/4,0)
        !                     = -1/8*DUC3+3/8*DUC1+3/4*DUM3
        !  DUF(IedgesAtElementFine(1,IEL2)) = P(1/4,0)
        !                     = 3/8*DUC3-1/8*DUC1+3/4*DUM3

        DuFine(IedgesAtElementFine(3,IEL1)+NVTfine) = -Q8*DUC3+3.0_DP*Q8*DUC1+3.0_DP*Q4*DUM3
        DuFine(IedgesAtElementFine(1,IEL2)+NVTfine) = 3.0_DP*Q8*DUC3-Q8*DUC1+3.0_DP*Q4*DUM3

      end if

      !  That is it - next element.

    end do ! IEL

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_restUniformP2_double (DuCoarse,DuFine, &
               IverticesAtElementCoarse,&
               IedgesAtElementCoarse,IedgesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               NELcoarse,NVTcoarse,NVTfine)

!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! $P_2$, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Coarse grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IverticesAtElement array (KVERT) on the coarse grid
  integer, dimension(:,:), intent(in) :: IverticesAtElementCoarse

  ! IedgesAtElement array (KMID) on the coarse grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementCoarse

  ! IedgesAtElement array (KMID) on the fine grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse

  ! Number of vertices in the coarse grid
  integer, intent(in) :: NVTcoarse

  ! Number of vertices in the fine grid
  integer, intent(in) :: NVTfine

!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

    ! local variables
    integer :: IEL,IEL2,IEL3,IEL4
    real(DP) :: dn1,dn2,dn3
    real(DP) :: duf11,duf12,duf13,duf21,duf23,duf31,duf33
    real(DP) :: duf41,DUF43
    real(DP), parameter :: Q8 = 0.125_DP
    real(DP), parameter :: Q4 = 0.25_DP
    real(DP), parameter :: Q2 = 0.5_DP

    ! First we remember the refinement scheme to clarify the
    ! local numbering of the vertices in the coarse- and fine-
    ! grid triangles.
    ! Let a coarse grid triangle be locally numbered as:
    !
    !   2
    !   |  \
    !   |    \
    !   | IEL  \
    !   |        \
    !   3----------1
    !
    ! Then the refinement process assigns the following numbers:
    !
    !   2
    !   |  \
    !   |     \
    !   |        \
    !   2*-------- 1*
    !   | \   IEL  |   \
    !   |    \     |      \
    !   |       \  |         \
    !   3----------3*----------1
    !
    ! i.e. the element number "IEL" is put into the middle with
    ! the corner vertices numbered according to the local edge
    ! numbers of the coarse grid element.
    !
    ! To access information on the edges of the fine grid element,
    ! we have to work with adjacencies!
    !
    ! Copy the first NVTC+NMTC values from DUF to DUC. This will
    ! transfer the contribution of the values from the
    ! fine-grid vertices that are coarse-grid vertices or midpoints
    ! as well. More precisely, this will transfer:
    !
    !   Fine grid vertices 1..NVTC -> Coarse grid vertices 1..NVTC
    !   Fine grid vertices NVTC+1..NVTC+NMTC = NVTF -> Coarse grid midpoints 1..NMTC
    !
    ! Afterwards, we have to add only the contribution of the fine grid
    ! edge midpoints to the coarse grid vertices/midpoints.

    call lalg_copyVector (DuFine(1:size(DuCoarse)),DuCoarse)

    ! loop over the elements

    do IEL = 1,nelCoarse

      ! The prolongation created from DUC1-3 and DUM1-3 in the
      ! following sketch the fine grid values DUFxy:
      !
      !    DUC2
      !     |  \
      !     |    \
      !   DUF21  DUF23
      !     |        \
      !     |          \
      !    DUM2--DUF11--DUM1
      !     |  \         |  \
      !     |    \       |    \
      !   DUF33  DUF12 DUF13  DUF41
      !     |        \   |        \
      !     |          \ |          \
      !    DUC3--DUF31--DUM3--DUF43--DUC1
      !
      ! This was done by a weighted averaging. The restriction now
      ! builds the DUC1-3 and DUM1-3 values as a weighted sum
      ! of themselves (DUM1-3, DUC1-3) and the new midpoints
      ! DUFxy using the same weights.
      !
      ! We had
      !   DUCx(fine grid) := 1*DUCx(coarse grid)
      !   DUMy(fine grid) := 1*DUCy(coarse grid)
      ! Therefore:
      !   DUCx(coarse grid) := 1*DUCx(fine grid) + ...
      !   DUCy(coarse grid) := 1*DUMy(fine grid) + ...
      !
      ! This part of the sum is already written to DUC with the
      ! above LCP1 command! now comes the rest of the sum.
      !
      ! The prolongation used the following formulas to calculate
      ! the fine grid vertices:
      !
      !  DUF11 = P(1/4,1/2)
      !        = -1/8*DUC3-1/8*DUC1+1/4*DUM3+1/2*DUM2+1/2*DUM1
      !  DUF12 = P(1/4,1/4)
      !        = -1/8*DUC1+1/2*DUM3-1/8*DUC2+1/2*DUM2+1/4*DUM1
      !  DUF13 = P(1/2,1/4)
      !        = -1/8*DUC3+1/2*DUM3-1/8*DUC2+1/4*DUM2+1/2*DUM1
      !
      !  DUF23 = P(1/4,3/4)
      !        = -1/8*DUC1+3/8*DUC2+3/4*DUM1
      !  DUF41 = P(3/4,1/4)
      !        = 3/8*DUC1-1/8*DUC2+3/4*DUM1
      !
      !  DUF21 = P(0,3/4)
      !        = -1/8*DUC3+3/8*DUC2+3/4*DUM2
      !  DUF33 = P(0,1/4)
      !        = 3/8*DUC3-1/8*DUC2+3/4*DUM2
      !
      !  DUF43 = P(3/4,0)
      !        = -1/8*DUC3+3/8*DUC1+3/4*DUM3
      !  DUF31 = P(1/4,0)
      !        = 3/8*DUC3-1/8*DUC1+3/4*DUM3
      !
      ! This is equivalent to the system
      !
      !  DUF11    [-1/8,    0, -1/8, 1/2,  1/2,  1/4]   DUC1
      !  DUF12    [-1/8, -1/8 ,   0, 1/4,  1/2,  1/2]   DUC2
      !  DUF13    [   0, -1/8, -1/8, 1/2,  1/4,  1/2]   DUC3
      !  DUF21    [   0,  3/8, -1/8,   0,  3/4,    0]   DUM1
      !  DUF23 =  [-1/8,  3/8,    0, 3/4,    0,    0] * DUM2
      !  DUF31    [-1/8,    0,  3/8,   0,    0,  3/4]   DUM3
      !  DUF33    [   0, -1/8,  3/8,   0,  3/4,    0]
      !  DUF41    [ 3/8, -1/8 ,   0, 3/4,    0,    0]
      !  DUF43    [ 3/8,    0, -1/8,   0,    0,  3/4]
      !
      ! Transposing it says what is left to add to DUC1-3/DUM1-3:
      !
      !   DUC1    [-1/8, -1/8,    0,    0, -1/8, -1/8,    0,  3/8,  3/8]   DUF11
      !   DUC2    [   0, -1/8, -1/8,  3/8,  3/8,    0, -1/8, -1/8,    0]   DUF12
      !   DUC3 += [-1/8,    0, -1/8, -1/8,    0,  3/8,  3/8,    0, -1/8] * DUF13
      !   DUM1    [ 1/2,  1/4,  1/2,    0,  3/4,    0,    0,  3/4,    0]   DUF21
      !   DUM2    [ 1/2,  1/2,  1/4, 3/4,     0,    0,  3/4,    0,    0]   DUF23
      !   DUM3    [ 1/4,  1/2,  1/2,    0,    0,  3/4,    0,    0,  3/4]   DUF31
      !                                                                    DUF33
      !                                                                    DUF41
      !                                                                    DUF43
      !
      ! Fetch the fine grid values into local variables:

      IEL2 = IneighboursAtElementFine(1,IEL)
      IEL3 = IneighboursAtElementFine(2,IEL)
      IEL4 = IneighboursAtElementFine(3,IEL)
      DUF11 = DuFine(IedgesAtElementFine(1,IEL)+NVTfine)
      DUF12 = DuFine(IedgesAtElementFine(2,IEL)+NVTfine)
      DUF13 = DuFine(IedgesAtElementFine(3,IEL)+NVTfine)
      DUF21 = DuFine(IedgesAtElementFine(1,IEL2)+NVTfine)
      DUF23 = DuFine(IedgesAtElementFine(3,IEL2)+NVTfine)
      DUF31 = DuFine(IedgesAtElementFine(1,IEL3)+NVTfine)
      DUF33 = DuFine(IedgesAtElementFine(3,IEL3)+NVTfine)
      DUF41 = DuFine(IedgesAtElementFine(1,IEL4)+NVTfine)
      DUF43 = DuFine(IedgesAtElementFine(3,IEL4)+NVTfine)

      ! When we add the information to DUC1-3/DUM1-3 we have to take
      ! into account whether there is a neighbor element or not!
      ! If there is a neighbor, we only add half of the information
      ! of the edge with the neighbor to DUC1-3/DUM1-3.
      ! In the loop over the elements here, we will
      ! later reach the neighbor and add another time half of the
      ! information, which that way completes that edge.
      !
      ! Set Nx=0.5 if there is a neighbor on edge x on the
      ! coarse grid or =1.0 if there is no neighbor

      dn1 = 0.5_DP
      dn2 = 0.5_DP
      dn3 = 0.5_DP

      if (IneighboursAtElementCoarse(1,IEL) .eq. 0) dn1 = 1.0_DP
      if (IneighboursAtElementCoarse(2,IEL) .eq. 0) dn2 = 1.0_DP
      if (IneighboursAtElementCoarse(3,IEL) .eq. 0) dn3 = 1.0_DP

      ! Now sum up the restriction.

      ! DUC1:

      DuCoarse(IverticesAtElementCoarse(1,IEL)) = DuCoarse(IverticesAtElementCoarse(1,IEL)) &
             -Q8*DUF11 -Q8*DUF12                        &
        +dn1*(-Q8*DUF23 +3.0_DP*Q8*DUF41)               &
        +dn3*(-Q8*DUF31 +3.0_DP*Q8*DUF43)

      ! DUC2:

      DuCoarse(IverticesAtElementCoarse(2,IEL)) = DuCoarse(IverticesAtElementCoarse(2,IEL)) &
             -Q8*DUF12 -Q8*DUF13                        &
        +dn2*(3.0_DP*Q8*DUF21-Q8*DUF33)                 &
        +dn1*(3.0_DP*Q8*DUF23-Q8*DUF41)

      ! DUC3:

      DuCoarse(IverticesAtElementCoarse(3,IEL)) = DuCoarse(IverticesAtElementCoarse(3,IEL)) &
             -Q8*DUF11 -Q8*DUF13                        &
        +dn2*(-Q8*DUF21 +3.0_DP*Q8*DUF33)               &
        +dn3*(-Q8*DUF43 +3.0_DP*Q8*DUF31)

      ! DUM1:

      DuCoarse(IedgesAtElementCoarse(1,IEL)+NVTcoarse) = DuCoarse(IedgesAtElementCoarse(1,IEL)+NVTcoarse) &
       +     Q2*DUF11 +Q4*DUF12 +Q2*DUF13                                             &
       + dn1*(3.0_DP*Q4*DUF23 +3*Q4*DUF41)

      ! DUM2:

      DuCoarse(IedgesAtElementCoarse(2,IEL)+NVTcoarse) = DuCoarse(IedgesAtElementCoarse(2,IEL)+NVTcoarse) &
       +     Q2*DUF11 +Q2*DUF12 +Q4*DUF13                                             &
       + dn2*(3.0_DP*Q4*DUF21+3.0_DP*Q4*DUF33)

      ! DUM3:

      DuCoarse(IedgesAtElementCoarse(3,IEL)+NVTcoarse) = DuCoarse(IedgesAtElementCoarse(3,IEL)+NVTcoarse) &
       +     Q4*DUF11 +Q2*DUF12 +Q2*DUF13                                             &
       +dn3*(3.0_DP*Q4*DUF43 +3.0_DP*Q4*DUF31)

      ! That is it - next element.

    end do ! IEL

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_interpUniformP2_double (DuCoarse,DuFine, &
                                           NVTcoarse, NMTcoarse)

!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! $P_2$, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! Number of vertices in the coarse grid
  integer, intent(in) :: NVTcoarse

  ! Number of edges in the coarse grid
  integer, intent(in) :: NMTcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables

    ! The first coase.NVT+NMT entries of the fine grid vector define
    ! the values on the coarse grid - because of the two-level ordering!
    call lalg_copyVector(DUfine(1:NVTcoarse+NMTcoarse),&
                         DUcoarse(1:NVTCoarse+NMTcoarse))

  end subroutine

  ! ***************************************************************************
  ! Support for P1~ element, DOF`s = edge midpoints
  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_prolUniformEx20_double (DuCoarse,DuFine, &
               IedgesAtElementCoarse,IedgesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               NELcoarse)

!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! E020, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Coarse grid vector
  real(DP), dimension(:), intent(in) :: DuCoarse

  ! IedgesAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementCoarse

  ! IedgesAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Fine grid vector
  real(DP), dimension(:), intent(out) :: DuFine
!</output>

!</subroutine>

  ! local variables
  real(DP) :: DUH1,DUH2,DUH3
  integer :: IM1,IM2,IM3, IA,IB,IC
  integer :: iel, IELH1, IELH2, IELH3, IELH4

  ! Weights for the prolongation; all coefficients are
  ! halfed, so dividing by 2 is not necessary in the calculation routines.
  real(DP), parameter :: A1=0.125_DP, A2=0.0625_DP, A3=-0.0625_DP, A4=0.25_DP

    ! Clear the output vector
    call lalg_clearVector(DuFine)

    ! Loop over the coarse grid elements
    do iel=1,NELcoarse

      ! Get the DOF`s of the coarse grid element
      IM1 = IedgesAtElementCoarse(1,iel)
      IM2 = IedgesAtElementCoarse(2,iel)
      IM3 = IedgesAtElementCoarse(3,iel)

      ! Get the values of the corresponding DOF`s
      DUH1 = DuCoarse(IM1)
      DUH2 = DuCoarse(IM2)
      DUH3 = DuCoarse(IM3)

      ! Get the element numbers of the fine-grid elements in the coarse grid;
      ! the numbers are defined by the two-level ordering.
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(1,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH1)
      IELH4 = IneighboursAtElementFine(3,IELH1)

      ! Distribute the value at the edge IM1 to the
      ! corresponding fine inner nodes

      IA=IedgesAtElementFine(3,IELH1)
      IB=IedgesAtElementFine(1,IELH3)

      if (IneighboursAtElementCoarse(1,iel).ne.0) then
        ! There is a neighbour at the edge
        DuFine(IA)=DuFine(IA)+A1*DUH1+A2*DUH2+A3*DUH3
        DuFine(IB)=DuFine(IB)+A1*DUH1+A3*DUH2+A2*DUH3
      else
        ! No neighbour; boundary element
        DuFine(IA)=DuFine(IA)+2.0_DP*(A1*DUH1+A2*DUH2+A3*DUH3)
        DuFine(IB)=DuFine(IB)+2.0_DP*(A1*DUH1+A3*DUH2+A2*DUH3)
      end if

      ! Distribute the value at the edge IM2 to the
      ! corresponding fine inner nodes

      IA=IedgesAtElementFine(3,IELH3)
      IB=IedgesAtElementFine(1,IELH2)

      if (IneighboursAtElementCoarse(2,iel).ne.0) then
        ! There is a neighbour at the edge
        DuFine(IA)=DuFine(IA)+A1*DUH2+A2*DUH3+A3*DUH1
        DuFine(IB)=DuFine(IB)+A1*DUH2+A3*DUH3+A2*DUH1
      else
        ! No neighbour; boundary element
        DuFine(IA)=DuFine(IA)+2.0_DP*(A1*DUH2+A2*DUH3+A3*DUH1)
        DuFine(IB)=DuFine(IB)+2.0_DP*(A1*DUH2+A3*DUH3+A2*DUH1)
      end if

      ! Distribute the value at the edge IM3 to the
      ! corresponding fine inner nodes

      IA=IedgesAtElementFine(3,IELH4)
      IB=IedgesAtElementFine(1,IELH3)

      if (IneighboursAtElementCoarse(3,iel).ne.0) then
        ! There is a neighbour at the edge
        DuFine(IA)=DuFine(IA)+A1*DUH3+A2*DUH1+A3*DUH2
        DuFine(IB)=DuFine(IB)+A1*DUH3+A3*DUH1+A2*DUH2
      else
        ! No neighbour; boundary element
        DuFine(IA)=DuFine(IA)+2.0_DP*(A1*DUH3+A2*DUH1+A3*DUH2)
        DuFine(IB)=DuFine(IB)+2.0_DP*(A1*DUH3+A3*DUH1+A2*DUH2)
      end if

      ! Distribute the value at the interior triangle to the
      ! corresponding fine inner nodes

      IA=IedgesAtElementFine(1,IELH1)
      IB=IedgesAtElementFine(2,IELH1)
      IC=IedgesAtElementFine(3,IELH1)

      DuFine(IA)=DuFine(IA)+A4*(DUH1+DUH2)
      DuFine(IB)=DuFine(IB)+A4*(DUH2+DUH3)
      DuFine(IC)=DuFine(IC)+A4*(DUH1+DUH3)

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_restUniformEx20_double (DuCoarse,DuFine, &
               IedgesAtElementCoarse,IedgesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,NELcoarse)

!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! E020 element, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IedgesAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementCoarse

  ! IedgesAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

    ! local variables
    real(DP) :: DUH1,DUH2,DUH3,DUH4,DUH5,DUH6,DUH7,DUH8,DUH9
    integer :: IM1,IM2,IM3
    integer :: I1,I2,I3,I4,I5,I6,I7,I8,I9
    integer :: iel,IELH1,IELH2,IELH3,IELH4

    ! Weights for the restriction; all coefficients are
    ! halfed, so dividing by 2 is not necessary in the calculation routines.
    real(DP), parameter :: A1=0.125_DP, A2=0.0625_DP, A3=-0.0625_DP

    ! Clear the output vector
    call lalg_clearVector(DuCoarse)

    ! Loop over the coarse grid elements
    do iel=1,NELcoarse

      ! Get the DOF`s of the coarse grid element
      IM1 = IedgesAtElementCoarse(1,iel)
      IM2 = IedgesAtElementCoarse(2,iel)
      IM3 = IedgesAtElementCoarse(3,iel)

      ! Get the element numbers of the fine-grid elements in the coarse grid;
      ! the numbers are defined by the two-level ordering.
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(1,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH1)
      IELH4 = IneighboursAtElementFine(3,IELH1)

      ! Get the DOF`s on the fine grid
      I1=IedgesAtElementFine(3,IELH2)
      I2=IedgesAtElementFine(1,IELH2)
      I3=IedgesAtElementFine(3,IELH3)
      I4=IedgesAtElementFine(1,IELH3)
      I5=IedgesAtElementFine(3,IELH4)
      I6=IedgesAtElementFine(1,IELH4)
      I7=IedgesAtElementFine(1,IELH1)
      I8=IedgesAtElementFine(2,IELH1)
      I9=IedgesAtElementFine(3,IELH1)

      ! Get the values of the DOF`s on the fine grid
      DUH1= DuFine(I1)
      DUH2= DuFine(I2)
      DUH3= DuFine(I3)
      DUH4= DuFine(I4)
      DUH5= DuFine(I5)
      DUH6= DuFine(I6)
      DUH7= DuFine(I7)
      DUH8= DuFine(I8)
      DUH9= DuFine(I9)

      ! Now collect the information in the same way as it was
      ! "distributed" by the prolongation routine.
      ! This realises the adjoint operator of the prolongation.

      ! Calculate the value of the edge IM1

      if (IneighboursAtElementCoarse(1,iel).ne.0) then
        ! inner edge
        DuCoarse(IM1)= DuCoarse(IM1) &
                     + A1*(DUH1+DUH6+DUH7+DUH9) &
                     + A2*(DUH2+DUH5) &
                     + A3*(DUH3+DUH4)
      else
        ! boundary edge
        DuCoarse(IM1)= 2.0_DP*( A1*(DUH1+DUH6+DUH7+DUH9) &
                              + A2*(DUH2+DUH5) &
                              + A3*(DUH3+DUH4) )
      end if

      ! Calculate the value of the edge IM2

      if (IneighboursAtElementCoarse(2,iel).ne.0) then
        ! inner edge
        DuCoarse(IM2)= DuCoarse(IM2) &
                     + A1*(DUH2+DUH3+DUH7+DUH8) &
                     + A2*(DUH1+DUH4) &
                     + A3*(DUH5+DUH6)
      else
        ! boundary edge
        DuCoarse(IM2)= 2.0_DP*( A1*(DUH2+DUH3+DUH7+DUH8) &
                              + A2*(DUH1+DUH4) &
                              + A3*(DUH5+DUH6) )
      end if

      ! Calculate the value of the edge IM3

      if (IneighboursAtElementCoarse(3,iel).ne.0) then
        ! inner edge
        DuCoarse(IM3)= DuCoarse(IM3) &
                     + A1*(DUH4+DUH5+DUH8+DUH9) &
                     + A2*(DUH3+DUH6) &
                     + A3*(DUH1+DUH2)
      else
        ! boundary edge
        DuCoarse(IM3)= 2.0_DP*( A1*(DUH4+DUH5+DUH8+DUH9) &
                              + A2*(DUH3+DUH6) &
                              + A3*(DUH1+DUH2) )
      end if

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_interpUniformEx20_double (DuCoarse,DuFine, &
               IedgesAtElementCoarse,IedgesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               NELcoarse)

!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! Ex20, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IedgesAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementCoarse

  ! IedgesAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

    ! local variables
    real(DP) :: DUH1,DUH2,DUH3,DUH4,DUH5,DUH6,DUH7,DUH8,DUH9
    integer :: IM1,IM2,IM3
    integer :: I1,I2,I3,I4,I5,I6,I7,I8,I9
    integer :: iel,IELH1,IELH2,IELH3,IELH4

    ! Weights for the restriction; all coefficients are
    ! halfed, so dividing by 2 is not necessary in the calculation routines.
    real(DP), parameter :: A1=0.125_DP, A2=0.0625_DP, A3=-0.0625_DP

    ! Clear the output vector
    call lalg_clearVector(DuCoarse)

    ! Loop over the coarse grid elements
    do iel=1,NELcoarse

      ! Get the DOF`s of the coarse grid element
      IM1 = IedgesAtElementCoarse(1,iel)
      IM2 = IedgesAtElementCoarse(2,iel)
      IM3 = IedgesAtElementCoarse(3,iel)

      ! Get the element numbers of the fine-grid elements in the coarse grid;
      ! the numbers are defined by the two-level ordering.
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(1,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH1)
      IELH4 = IneighboursAtElementFine(3,IELH1)

      ! Get the DOF`s on the fine grid
      I1=IedgesAtElementFine(3,IELH2)
      I2=IedgesAtElementFine(1,IELH2)
      I3=IedgesAtElementFine(3,IELH3)
      I4=IedgesAtElementFine(1,IELH3)
      I5=IedgesAtElementFine(3,IELH4)
      I6=IedgesAtElementFine(1,IELH4)
      I7=IedgesAtElementFine(1,IELH1)
      I8=IedgesAtElementFine(2,IELH1)
      I9=IedgesAtElementFine(3,IELH1)

      ! Get the values of the DOF`s on the fine grid
      DUH1= DuFine(I1)
      DUH2= DuFine(I2)
      DUH3= DuFine(I3)
      DUH4= DuFine(I4)
      DUH5= DuFine(I5)
      DUH6= DuFine(I6)
      DUH7= DuFine(I7)
      DUH8= DuFine(I8)
      DUH9= DuFine(I9)

      ! Now collect the information in the same way as it was
      ! "distributed" by the prolongation routine.
      ! This realises the adjoint operator of the prolongation.

      ! Calculate the value of the edge IM1

      if (IneighboursAtElementCoarse(1,iel).ne.0) then
        ! inner edge
        DuCoarse(IM1)= DuCoarse(IM1) &
                     + A1*(DUH1+DUH6+DUH7+DUH9) &
                     + A2*(DUH2+DUH5) &
                     + A3*(DUH3+DUH4)
      else
        ! boundary edge
        DuCoarse(IM1)= 2.0_DP*( A1*(DUH1+DUH6+DUH7+DUH9) &
                              + A2*(DUH2+DUH5) &
                              + A3*(DUH3+DUH4) )
      end if

      ! Calculate the value of the edge IM2

      if (IneighboursAtElementCoarse(2,iel).ne.0) then
        ! inner edge
        DuCoarse(IM2)= DuCoarse(IM2) &
                     + A1*(DUH2+DUH3+DUH7+DUH8) &
                     + A2*(DUH1+DUH4) &
                     + A3*(DUH5+DUH6)
      else
        ! boundary edge
        DuCoarse(IM2)= 2.0_DP*( A1*(DUH2+DUH3+DUH7+DUH8) &
                              + A2*(DUH1+DUH4) &
                              + A3*(DUH5+DUH6) )
      end if

      ! Calculate the value of the edge IM3

      if (IneighboursAtElementCoarse(3,iel).ne.0) then
        ! inner edge
        DuCoarse(IM3)= DuCoarse(IM3) &
                     + A1*(DUH4+DUH5+DUH8+DUH9) &
                     + A2*(DUH3+DUH6) &
                     + A3*(DUH1+DUH2)
      else
        ! boundary edge
        DuCoarse(IM3)= 2.0_DP*( A1*(DUH4+DUH5+DUH8+DUH9) &
                              + A2*(DUH3+DUH6) &
                              + A3*(DUH1+DUH2) )
      end if

    end do

  end subroutine

  ! ***************************************************************************
  ! Support for Q0 element
  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_prolUniformQ0_double (DuCoarse,DuFine, &
               IneighboursAtElementFine, NELcoarse)

!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! $Q_0$, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Coarse grid vector
  real(DP), dimension(:), intent(in) :: DuCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Fine grid vector
  real(DP), dimension(:), intent(out) :: DuFine
!</output>

!</subroutine>

  ! local variables
  integer :: iel
  integer :: IELH1,IELH2,IELH3,IELH4
  real(DP) :: duh

    ! Loop over the elements
    do iel=1,NELCoarse

      ! Get the four child elements of the coarse grid element
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(2,IELH3)

      ! Put the value on the coarse grid into all four child
      ! elements
      duh = DuCoarse(iel)
      DuFine(IELH1) = duh
      DuFine(IELH2) = duh
      DuFine(IELH3) = duh
      DuFine(IELH4) = duh
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_restUniformQ0_double (DuCoarse,DuFine, &
               IneighboursAtElementFine, NELcoarse)

!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! $Q_0$, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  integer :: iel
  integer :: IELH1,IELH2,IELH3,IELH4

    ! The information that was "distributed" in the prolongation has to
    ! be "collected".

    ! Loop over the elements to collect the missing additive contributions:
    do iel=1,NELcoarse

      ! Get the elements on the fine grid that are children of the
      ! coarse grid element
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(2,IELH3)

      ! Sum up the values in these nodes to get the
      ! value in the coarse grid element
      DuCoarse(iel)= DuFine(IELH1)+DuFine(IELH2)+DuFine(IELH3)+DuFine(IELH4)

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_interpUniformQ0_double (DuCoarse,DuFine, &
               IneighboursAtElementFine, NELcoarse)

!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! $Q_0$, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  real(DP), parameter :: Q4 = 0.25_DP

  integer :: iel
  integer :: IELH1,IELH2,IELH3,IELH4

    ! The information that was "distributed" in the prolongation has to
    ! be "collected".

    ! Loop over the elements to collect the missing additive contributions:
    do iel=1,NELcoarse

      ! Get the elements on the fine grid that are children of the
      ! coarse grid element
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(2,IELH3)

      ! From all four values, build the mean and use that as
      ! value of the coarse grid element
      DuCoarse(iel)= Q4*(DuFine(IELH1)+DuFine(IELH2)+DuFine(IELH3)+DuFine(IELH4))

    end do

  end subroutine

  ! ***************************************************************************
  ! Support for Q1 element
  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_prolUniformQ1_double (DuCoarse,DuFine, &
               IverticesAtEdgeCoarse, IverticesAtElementCoarse, &
               NVTcoarse, NMTcoarse, NELcoarse)
! "old" parameter list
!               IverticesAtElementCoarse,IverticesAtElementFine,&
!               IneighboursAtElementCoarse,IneighboursAtElementFine,NELcoarse)

!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! <tex>$Q_1$</tex>, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Coarse grid vector
  real(DP), dimension(:), intent(in) :: DuCoarse

  ! IverticesAtEdge array on the coarse grid.
  integer, dimension(:,:), intent(in) :: IverticesAtEdgeCoarse

  ! IverticesAtElement array (KVERT) on the coarse grid
  integer, dimension(:,:), intent(in) :: IverticesAtElementCoarse

!  ! IverticesAtElement array (KVERT) on the fine grid
!  INTEGER, DIMENSION(:,:), INTENT(in) :: IverticesAtElementFine
!
!  ! IneighboursAtElement array on the coarse grid
!  INTEGER, DIMENSION(:,:), INTENT(in) :: IneighboursAtElementCoarse
!
!  ! IneighboursAtElement array on the fine grid
!  INTEGER, DIMENSION(:,:), INTENT(in) :: IneighboursAtElementFine

  ! Number of vertices in the coarse grid
  integer, intent(in) :: NVTcoarse

  ! Number of edges in the coarse grid
  integer, intent(in) :: NMTcoarse

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Fine grid vector
  real(DP), dimension(:), intent(out) :: DuFine
!</output>

!</subroutine>

  ! local variables
  integer :: iel
  integer :: iedge

    ! Copy the first NVT entries - they belong to the coarse grid vertices
    ! that are fine grid vertices at the same time.
    call lalg_copyVector (DuCoarse,DuFine(1:NVTcoarse))

    ! Loop over the edges
    do iedge = 1, NMTcoarse
      ! Calculate the edge midpoint DOF
      DuFine(NVTcoarse + iedge) = 0.5_DP * (&
          DuCoarse(IverticesAtEdgeCoarse(1,iedge))+&
          DuCoarse(IverticesAtEdgeCoarse(2,iedge)))

    end do

    ! Loop over the elements
    do iel = 1, NELcoarse
      ! Calculate the quad cell midpoint DOF
      DuFine(NVTcoarse + NMTcoarse + iel) = 0.25_DP * (&
          DuCoarse(IverticesAtElementCoarse(1,iel))+&
          DuCoarse(IverticesAtElementCoarse(2,iel))+&
          DuCoarse(IverticesAtElementCoarse(3,iel))+&
          DuCoarse(IverticesAtElementCoarse(4,iel)))
    end do

!  This is the "old" implemention, based on Feat 1.x
!  ! local variables
!  REAL(DP), PARAMETER :: Q2 = 0.5_DP
!  REAL(DP), PARAMETER :: Q4 = 0.25_DP
!
!  INTEGER :: iel,ielh1,ielh2,ielh3,ielh4
!  REAL(DP) :: duh1,duh2,duh3,duh4
!
!    ! Copy the first NVT entries - they belong to the coarse grid vertices
!    ! that are fine grid vertices at the same time.
!    CALL lalg_copyVector (DuCoarse,DuFine(1:SIZE(DuCoarse)))
!
!    ! Loop over the elements
!    DO iel=1,NELCoarse
!
!      duh1=DuCoarse(IverticesAtElementCoarse(1,iel))
!      duh2=DuCoarse(IverticesAtElementCoarse(2,iel))
!      duh3=DuCoarse(IverticesAtElementCoarse(3,iel))
!      duh4=DuCoarse(IverticesAtElementCoarse(4,iel))
!
!      ielh1=iel
!      ielh2=IneighboursAtElementFine(2,ielh1)
!      ielh3=IneighboursAtElementFine(2,ielh2)
!      ielh4=IneighboursAtElementFine(2,ielh3)
!
!      ! Now check on every of the edges, if we already computed
!      ! the value in the midpoint: Compute only if the neighbour
!      ! element has smaller number.
!      IF (IneighboursAtElementCoarse(1,iel) .LT. iel) &
!        DuFine(IverticesAtElementFine(2,ielh1)) = Q2*(duh1+duh2)
!
!      IF (IneighboursAtElementCoarse(2,iel) .LT. iel) &
!        DuFine(IverticesAtElementFine(2,ielh2)) = Q2*(duh2+duh3)
!
!      IF (IneighboursAtElementCoarse(3,iel) .LT. iel) &
!        DuFine(IverticesAtElementFine(2,ielh3)) = Q2*(duh3+duh4)
!
!      IF (IneighboursAtElementCoarse(4,iel) .LT. iel) &
!        DuFine(IverticesAtElementFine(2,ielh4)) = Q2*(duh4+duh1)
!
!      ! Do not forget the DOF in the midpoint of the element
!      DuFine(IverticesAtElementFine(3,iel)) = Q4*(duh1+duh2+duh3+duh4)
!
!    END DO

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_restUniformQ1_double (DuCoarse,DuFine, &
             IverticesAtEdgeCoarse, IverticesAtElementCoarse, &
             NVTcoarse, NMTcoarse, NELcoarse)
!               IverticesAtElementFine,IneighboursAtElementFine,&
!               NELfine)

!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! Q1, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IverticesAtEdge array on the coarse grid.
  integer, dimension(:,:), intent(in) :: IverticesAtEdgeCoarse

  ! IverticesAtElement array (KVERT) on the coarse grid
  integer, dimension(:,:), intent(in) :: IverticesAtElementCoarse

  ! Number of vertices in the coarse grid
  integer, intent(in) :: NVTcoarse

  ! Number of edges in the coarse grid
  integer, intent(in) :: NMTcoarse

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse

! "old" parameters
!  ! IverticesAtElement array (KVERT) on the fine grid
!  INTEGER, DIMENSION(:,:), INTENT(in) :: IverticesAtElementFine
!
!  ! IneighboursAtElement array on the coarse grid
!  INTEGER, DIMENSION(:,:), INTENT(in) :: IneighboursAtElementFine
!
!  ! Number of elements in the fine grid
!  INTEGER, INTENT(in) :: NELfine
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  integer :: iel
  integer :: iedge
  integer :: ivt
  real(DP) :: dx

    ! The information that was "distributed" in the prolongation has to
    ! be "collected".
    !
    ! Copy the first NVT entries - this gives the first additive contribution.
    call lalg_copyVector (DuFine(1:NVTcoarse),DuCoarse)

    ! Loop over the edges
    do iedge = 1, NMTcoarse
      ! get the fine grid DOF
      dx = 0.5_DP * DuFine(NVTcoarse + iedge)

      ! distribute it to the coarse grid DOFs
      ivt = IverticesAtEdgeCoarse(1,iedge)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
      ivt = IverticesAtEdgeCoarse(2,iedge)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
    end do

    ! Loop over the elements
    do iel = 1, NELcoarse
      ! get the fine grid DOF
      dx = 0.25_DP * DuFine(NVTcoarse + NMTcoarse + iel)

      ! distribute it to the coarse grid DOFs
      ivt = IverticesAtElementCoarse(1, iel)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
      ivt = IverticesAtElementCoarse(2, iel)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
      ivt = IverticesAtElementCoarse(3, iel)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
      ivt = IverticesAtElementCoarse(4, iel)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
    end do

! This is the "old" implementation, based on Feat 1.x
!  ! local variables
!  REAL(DP), PARAMETER :: Q2 = 0.5_DP
!  REAL(DP), PARAMETER :: Q4 = 0.25_DP
!
!  INTEGER :: iel
!  INTEGER :: i1,i2,i3,i4
!
!    ! The information that was "distributed" in the prolongation has to
!    ! be "collected".
!    !
!    ! Copy the first NVT entries - this gives the first additive contribution.
!    CALL lalg_copyVector (DuFine(1:SIZE(DuCoarse)),DuCoarse)
!
!    ! Loop over the elements to collect the missing additive contributions:
!    DO iel=1,NELfine
!      i1=IverticesAtElementFine(1,iel)
!      i2=IverticesAtElementFine(2,iel)
!      i3=IverticesAtElementFine(3,iel)
!      i4=IverticesAtElementFine(4,iel)
!
!      ! Additive contribution of the midpoint
!      DuCoarse(i1) = DuCoarse(i1)+Q4*(DuFine(i2)+DuFine(i3)+DuFine(i4))
!
!      ! Additional contribution on the boundary:
!      IF (IneighboursAtElementFine(1,iel) .EQ. 0) &
!        DuCoarse(i1) = DuCoarse(i1)+Q4*DuFine(i2)
!      IF (IneighboursAtElementFine(4,iel) .EQ. 0) &
!        DuCoarse(i1) = DuCoarse(i1)+Q4*DuFine(i4)
!    END DO

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_prolUnifQ1FMzero_double (DuCoarse,DuFine, &
               IverticesAtElementCoarse,IverticesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               NELcoarse,NELfine)

!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! <tex>$Q_1$</tex>, uniform triangulation, double precision vector.
  ! Experimental FEAST MIRROR variant, zeroes the entries on the boundary.
!</description>

!<input>
  ! Coarse grid vector
  real(DP), dimension(:), intent(in) :: DuCoarse

  ! IverticesAtElement array (KVERT) on the coarse grid
  integer, dimension(:,:), intent(in) :: IverticesAtElementCoarse

  ! IverticesAtElement array (KVERT) on the fine grid
  integer, dimension(:,:), intent(in) :: IverticesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse

  ! Number of elements in the fine grid
  integer, intent(in) :: NELfine
!</input>

!<output>
  ! Fine grid vector
  real(DP), dimension(:), intent(out) :: DuFine
!</output>

!</subroutine>

  ! local variables
  real(DP), parameter :: Q2 = 0.5_DP
  real(DP), parameter :: Q4 = 0.25_DP

  integer :: iel,ielh1,ielh2,ielh3,ielh4
  integer :: i1,i2,i3,i4
  real(DP) :: duh1,duh2,duh3,duh4

    ! Copy the first NVT entries - they belong to the coarse grid vertices
    ! that are fine grid vertices at the same time.
    call lalg_copyVector (DuCoarse,DuFine(1:size(DuCoarse)))

    ! Loop over the elements
    do iel=1,NELCoarse

      duh1=DuCoarse(IverticesAtElementCoarse(1,iel))
      duh2=DuCoarse(IverticesAtElementCoarse(2,iel))
      duh3=DuCoarse(IverticesAtElementCoarse(3,iel))
      duh4=DuCoarse(IverticesAtElementCoarse(4,iel))

      ielh1=iel
      ielh2=IneighboursAtElementFine(2,ielh1)
      ielh3=IneighboursAtElementFine(2,ielh2)
      ielh4=IneighboursAtElementFine(2,ielh3)

      ! Now check on every of the edges, if we already computed
      ! the value in the midpoint: Compute only if the neighbour
      ! element has smaller number.
      if (IneighboursAtElementCoarse(1,iel) .lt. iel) &
        DuFine(IverticesAtElementFine(2,ielh1)) = Q2*(duh1+duh2)

      if (IneighboursAtElementCoarse(2,iel) .lt. iel) &
        DuFine(IverticesAtElementFine(2,ielh2)) = Q2*(duh2+duh3)

      if (IneighboursAtElementCoarse(3,iel) .lt. iel) &
        DuFine(IverticesAtElementFine(2,ielh3)) = Q2*(duh3+duh4)

      if (IneighboursAtElementCoarse(4,iel) .lt. iel) &
        DuFine(IverticesAtElementFine(2,ielh4)) = Q2*(duh4+duh1)

      ! Do not forget the DOF in the midpoint of the element
      DuFine(IverticesAtElementFine(3,iel)) = Q4*(duh1+duh2+duh3+duh4)

    end do

    ! DOF`s on the boundary get value 0.0.
    do iel=1,NELfine
      i1=IverticesAtElementFine(1,iel)
      i2=IverticesAtElementFine(2,iel)
      i3=IverticesAtElementFine(3,iel)
      i4=IverticesAtElementFine(4,iel)

      ! Additional contribution on the boundary:
      if (IneighboursAtElementFine(1,iel) .eq. 0) then
        DuFine(i1) = 0.0_DP
        DuFine(i2) = 0.0_DP
      end if
      if (IneighboursAtElementFine(2,iel) .eq. 0) then
        DuFine(i2) = 0.0_DP
        DuFine(i3) = 0.0_DP
      end if
      if (IneighboursAtElementFine(3,iel) .eq. 0) then
        DuFine(i3) = 0.0_DP
        DuFine(i4) = 0.0_DP
      end if
      if (IneighboursAtElementFine(4,iel) .eq. 0) then
        DuFine(i4) = 0.0_DP
        DuFine(i1) = 0.0_DP
      end if
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_restUniformQ1FM_double (DuCoarse,DuFine, &
               IverticesAtElementFine,IneighboursAtElementFine,&
               NELfine)

!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! Q1, uniform triangulation, double precision vector.
  ! Cellwise approach for FEAST mirror boundary.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IverticesAtElement array (KVERT) on the fine grid
  integer, dimension(:,:), intent(in) :: IverticesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the fine grid
  integer, intent(in) :: NELfine
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  real(DP), parameter :: Q2 = 0.5_DP
  real(DP), parameter :: Q4 = 0.25_DP

  integer :: iel
  integer :: i1,i2,i3,i4

    ! The information that was "distributed" in the prolongation has to
    ! be "collected".
    call lalg_clearVector (DuCoarse)

    ! Loop over the elements to collect the missing additive contributions:
    do iel=1,NELfine
      i1=IverticesAtElementFine(1,iel)
      i2=IverticesAtElementFine(2,iel)
      i3=IverticesAtElementFine(3,iel)
      i4=IverticesAtElementFine(4,iel)

      ! Additive contribution of all vertices to the coarse grid vertex
      DuCoarse(i1) = DuCoarse(i1)+Q4*(DUfine(i1)+DuFine(i2)+DuFine(i3)+DuFine(i4))

      ! Additional contribution on the boundary:
      if (IneighboursAtElementFine(1,iel) .eq. 0) then
        DuCoarse(i1) = DuCoarse(i1)+2.0_DP*Q4*DuFine(i2)
        DuCoarse(i1) = DuCoarse(i1)+2.0_DP*Q4*DuFine(i1)
      end if
      if (IneighboursAtElementFine(4,iel) .eq. 0) then
        DuCoarse(i1) = DuCoarse(i1)+2.0_DP*Q4*DuFine(i4)
        DuCoarse(i1) = DuCoarse(i1)+2.0_DP*Q4*DuFine(i1)
      end if
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_restUnifQ1FMzero_double (DuCoarse,DuFine,&
               IverticesAtElementFine,IneighboursAtElementFine,&
               IverticesAtElementCoarse,IneighboursAtElementCoarse,&
               NELfine,NELcoarse)

!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! Q1, uniform triangulation, double precision vector.
  ! Cellwise approach for FEAST mirror boundary.
  ! The DOF`s on the boundary are set to 0.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IverticesAtElement array (KVERT) on the fine grid
  integer, dimension(:,:), intent(in) :: IverticesAtElementFine

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! IverticesAtElement array (KVERT) on the coarse grid
  integer, dimension(:,:), intent(in) :: IverticesAtElementCoarse

  ! IneighboursAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementCoarse

  ! Number of elements in the fine grid
  integer, intent(in) :: NELfine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  real(DP), parameter :: Q2 = 0.5_DP
  real(DP), parameter :: Q4 = 0.25_DP

  integer :: iel
  integer :: i1,i2,i3,i4

    ! The information that was "distributed" in the prolongation has to
    ! be "collected".
    !
    ! Copy the first NVT entries - this gives the first additive contribution.
    call lalg_copyVector (DuFine(1:size(DuCoarse)),DuCoarse)

    ! Loop over the elements to collect the missing additive contributions:
    do iel=1,NELfine
      i1=IverticesAtElementFine(1,iel)
      i2=IverticesAtElementFine(2,iel)
      i3=IverticesAtElementFine(3,iel)
      i4=IverticesAtElementFine(4,iel)

      ! Additive contribution of the midpoint
      DuCoarse(i1) = DuCoarse(i1)+Q4*(DuFine(i2)+DuFine(i3)+DuFine(i4))

      ! Additional contribution on the boundary:
      if (IneighboursAtElementFine(1,iel) .eq. 0) &
        DuCoarse(i1) = DuCoarse(i1)+Q4*DuFine(i2)
      if (IneighboursAtElementFine(4,iel) .eq. 0) &
        DuCoarse(i1) = DuCoarse(i1)+Q4*DuFine(i4)
    end do

    ! DOF`s on the boundary get value 0.0.
    do iel=1,NELcoarse
      i1=IverticesAtElementCoarse(1,iel)
      i2=IverticesAtElementCoarse(2,iel)
      i3=IverticesAtElementCoarse(3,iel)
      i4=IverticesAtElementCoarse(4,iel)

      ! Additional contribution on the boundary:
      if (IneighboursAtElementCoarse(1,iel) .eq. 0) then
        DuCoarse(i1) = 0.0_DP
        DuCoarse(i2) = 0.0_DP
      end if
      if (IneighboursAtElementCoarse(2,iel) .eq. 0) then
        DuCoarse(i2) = 0.0_DP
        DuCoarse(i3) = 0.0_DP
      end if
      if (IneighboursAtElementCoarse(3,iel) .eq. 0) then
        DuCoarse(i3) = 0.0_DP
        DuCoarse(i4) = 0.0_DP
      end if
      if (IneighboursAtElementCoarse(4,iel) .eq. 0) then
        DuCoarse(i4) = 0.0_DP
        DuCoarse(i1) = 0.0_DP
      end if
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_interpUniformQ1_double (DuCoarse,DuFine, NVTcoarse)

!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! Q1, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! Number of vertices in the coarse grid
  integer, intent(in) :: NVTcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

    ! The first coase.NVT entries of the fine grid vector define the values
    ! on the coarse grid - because of the two-level ordering!
    call lalg_copyVector(DUfine(1:NVTcoarse),DUcoarse(1:NVTCoarse))

  end subroutine

  ! ***************************************************************************
  ! Support for Q2 element
  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_prolUniformQ2_double (DuCoarse,DuFine, &
               IverticesAtElementFine,IedgesAtElementFine,&
               IneighboursAtElementFine,&
               NVTfine, NMTfine, NELcoarse)

!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! <tex>$Q_2$</tex>, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Coarse grid vector
  real(DP), dimension(:), intent(in) :: DuCoarse

  ! IverticesAtElement array (KVERT) on the fine grid
  integer, dimension(:,:), intent(in) :: IverticesAtElementFine

  ! IedgesAtElement array (KMID) on the fine grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementFine

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of vertices on the fine grid
  integer, intent(in) :: NVTfine

  ! Number of edges in the fine grid
  integer, intent(in) :: NMTfine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Fine grid vector
  real(DP), dimension(:), intent(out) :: DuFine
!</output>

!</subroutine>

  ! local variables
  integer :: iel
  integer :: i
  integer, dimension(4) :: IelFine

    ! Copy the first NVT+NMT+NEL entries - they belong to the coarse grid
    ! vertices/edge midpoints/element midpoints and
    ! are fine grid vertices at the same time.
    call lalg_copyVector (DuCoarse,DuFine(1:size(DuCoarse)))

    ! Loop over the elements of the coarse grid
    do iel=1,NELcoarse

      ! Obtain the numbers of the fine grid elements inside of the
      ! coarse grid element. According to the regular refinement, the element
      ! number of the coarse grid element is the element number of the
      ! first element inside the coarse grid element on the fine grid.
      IelFine(1)=iel
      IelFine(2)=IneighboursAtElementFine(2,IelFine(1))
      IelFine(3)=IneighboursAtElementFine(2,IelFine(2))
      IelFine(4)=IneighboursAtElementFine(2,IelFine(3))

      ! Loop over the fine grid elements in the coarse grid element.
      ! "Distribute" the information from the edge midpoints and the element
      ! midpoint to the edge midpoints/element midpoints of the
      ! fine grid element.
      do i=1,4
        ! Distribute information on the edges of the coarse grid element
        ! to the edges of the fine grid element i inside of the coarse
        ! grid element.
        DUfine(IedgesAtElementFine(1,IelFine(i))+NVTfine)= &
             +(3.0/8.0)*DUfine(IverticesAtElementFine(1,IelFine(i))) &
             +(3.0/4.0)*DUfine(IverticesAtElementFine(2,IelFine(i))) &
             -(1.0/8.0)*DUfine(IverticesAtElementFine(1,IelFine(mod(i,4)+1)))
        DUfine(IedgesAtElementFine(4,IelFine(i))+NVTfine)= &
             +(3.0/8.0)*DUfine(IverticesAtElementFine(1,IelFine(i))) &
             +(3.0/4.0)*DUfine(IverticesAtElementFine(4,IelFine(i))) &
             -(1.0/8.0)*DUfine(IverticesAtElementFine(1,IelFine(mod(i+2,4)+1)))
        DUfine(IedgesAtElementFine(2,IelFine(i))+NVTfine)= &
             +(3.0/8.0)*DUfine(IverticesAtElementFine(2,IelFine(i))) &
             +(3.0/4.0)*DUfine(IverticesAtElementFine(3,IelFine(i))) &
             -(1.0/8.0)*DUfine(IverticesAtElementFine(4,IelFine(mod(i+2,4)+1)))

        ! Distribute information of the coarse grid midpoint to
        ! the fine grid midpoint of the fine grid element i inside
        ! of the coarse grid element.
        DUfine(NVTfine+NMTfine+IelFine(i))= &
             +(9.0/64.0)*DUfine(IverticesAtElementFine(1,IelFine(i))) &
             +(18.0/64.0)*DUfine(IverticesAtElementFine(2,IelFine(i))) &
             +(36.0/64.0)*DUfine(IverticesAtElementFine(3,IelFine(i))) &
             +(18.0/64.0)*DUfine(IverticesAtElementFine(4,IelFine(i))) &
             -(3.0/64.0)*DUfine(IverticesAtElementFine(1,IelFine(mod(i,4)+1))) &
             -(6.0/64.0)*DUfine(IverticesAtElementFine(2,IelFine(mod(i,4)+1))) &
             -(3.0/64.0)*DUfine(IverticesAtElementFine(1,IelFine(mod(i+2,4)+1))) &
             -(6.0/64.0)*DUfine(IverticesAtElementFine(4,IelFine(mod(i+2,4)+1))) &
             +(1.0/64.0)*DUfine(IverticesAtElementFine(1,IelFine(mod(i+1,4)+1)))
      enddo
    enddo

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_restUniformQ2_double (DuCoarse,DuFine, &
               IverticesAtElementCoarse, IedgesAtElementCoarse, &
               IedgesAtElementFine, IneighboursAtElementFine,&
               NVTcoarse,NVTfine,NMTcoarse,NMTfine,NELcoarse)

!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! <tex>$Q_2$</tex>, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IverticesAtElement array (KVERT) on the coarse grid
  integer, dimension(:,:), intent(in) :: IverticesAtElementCoarse

  ! IedgesAtElement array (KMID) on the coarse grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementCoarse

  ! IedgesAtElement array (KMID) on the fine grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of vertices in the coarse grid
  integer, intent(in) :: NVTcoarse

  ! Number of vertices in the fine grid
  integer, intent(in) :: NVTfine

  ! Number of edges in the coarse grid
  integer, intent(in) :: NMTcoarse

  ! Number of elements in the fine grid
  integer, intent(in) :: NMTfine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  integer :: iel
  integer :: i
  integer, dimension(4) :: IelFine

    ! The information that was "distributed" in the prolongation has to
    ! be "collected".
    !
    ! Copy the first NVT+NMT+NEL (coarse) entries - this gives the first
    ! additive contribution: The values in the corners/edge midpoints/
    ! element midpoints of the coarse grid stem from with the values of the
    ! corners of the fine grid.
    call lalg_copyVector (DuFine(1:size(DuCoarse)),DuCoarse)

    ! Loop over the elements to collect the missing additive contributions:
    do iel=1,NELcoarse

      ! Obtain the numbers of the fine grid elements inside of the
      ! coarse grid element. According to the regular refinement, the element
      ! number of the coarse grid element is the element number of the
      ! first element inside the coarse grid element on the fine grid.
      IelFine(1)=iel
      IelFine(2)=IneighboursAtElementFine(2,IelFine(1))
      IelFine(3)=IneighboursAtElementFine(2,IelFine(2))
      IelFine(4)=IneighboursAtElementFine(2,IelFine(3))

      do i=1,4

        ! Collect information from the corners of the fine grid element
        DUcoarse(IverticesAtElementCoarse(i,iel)) = &
             DUcoarse(IverticesAtElementCoarse(i,iel)) &
             +0.5*(24.0/64.0)*DUfine(IEdgesAtElementFine(1,IelFine(i))+NVTfine) &
             +0.5*(24.0/64.0)*DUfine(IEdgesAtElementFine(4,IelFine(i))+NVTfine) &
             -0.5*(8.0/64.0)*DUfine(IEdgesAtElementFine(4,IelFine(mod(i,4)+1))+NVTfine) &
             -0.5*(8.0/64.0)*DUfine(IEdgesAtElementFine(1,IelFine(mod(i+2,4)+1))+NVTfine) &
             +(9.0/64.0)*DUfine(NVTfine+NMTFine+IelFine(i)) &
             -(3.0/64.0)*DUfine(NVTfine+NMTFine+IelFine(mod(i,4)+1)) &
             +(1.0/64.0)*DUfine(NVTfine+NMTFine+IelFine(mod(i+1,4)+1)) &
             -(3.0/64.0)*DUfine(NVTfine+NMTFine+IelFine(mod(i+2,4)+1))

        if(IneighboursAtElementFine(1,IelFine(i)).eq.0) then
          DUcoarse(IverticesAtElementCoarse(i,iel))= &
              DUcoarse(IverticesAtElementCoarse(i,iel)) &
              +0.5*(24.0/64.0)*DUfine(IEdgesAtElementFine(1,IelFine(i))+NVTfine)
        end if

        if(IneighboursAtElementFine(4,IelFine(i)).eq.0) then
          DUcoarse(IverticesAtElementCoarse(i,iel))= &
              DUcoarse(IverticesAtElementCoarse(i,iel)) &
              +0.5*(24.0/64.0)*DUfine(IEdgesAtElementFine(4,IelFine(i))+NVTfine)
        end if

        if(IneighboursAtElementFine(4,IelFine(mod(i,4)+1)).eq.0) then
          DUcoarse(IverticesAtElementCoarse(i,iel))= &
              DUcoarse(IverticesAtElementCoarse(i,iel)) &
              -0.5*(8.0/64.0)*DUfine(IEdgesAtElementFine(4,IelFine(mod(i,4)+1))+NVTfine)
        end if

        if(IneighboursAtElementFine(1,IelFine(mod(i+2,4)+1)).eq.0) then
          DUcoarse(IverticesAtElementCoarse(i,iel))= &
              DUcoarse(IverticesAtElementCoarse(i,iel)) &
              -0.5*(8.0/64.0)*DUfine(IEdgesAtElementFine(1,IelFine(mod(i+2,4)+1))+NVTfine)
        end if

        ! Collect information from the edge midpoints of the fine grid element
        DUcoarse(IEdgesAtElementCoarse(i,iel)+NVTcoarse)= &
            DUcoarse(IEdgesAtElementCoarse(i,iel)+NVTcoarse) &
             +0.5*(48.0/64.0)*DUfine(IEdgesAtElementFine(1,IelFine(i))+NVTfine) &
             +(24.0/64.0)*DUfine(IEdgesAtElementFine(2,IelFine(i))+NVTfine) &
             +0.5*(48.0/64.0)*DUfine(IEdgesAtElementFine(4,IelFine(mod(i,4)+1))+NVTfine) &
             -(8.0/64.0)*DUfine(IEdgesAtElementFine(2,IelFine(mod(i+1,4)+1))+NVTfine) &
             +(18.0/64.0)*DUfine(NVTfine+NMTFine+IelFine(i)) &
             +(18.0/64.0)*DUfine(NVTfine+NMTFine+IelFine(mod(i,4)+1)) &
             -(6.0/64.0)*DUfine(NVTfine+NMTFine+IelFine(mod(i+1,4)+1)) &
             -(6.0/64.0)*DUfine(NVTfine+NMTFine+IelFine(mod(i+2,4)+1))

        if(IneighboursAtElementFine(1,IelFine(i)).eq.0) then
          DUcoarse(IEdgesAtElementCoarse(i,iel)+NVTcoarse)= &
              DUcoarse(IEdgesAtElementCoarse(i,iel)+NVTcoarse) &
              +0.5*(48.0/64.0)*DUfine(IEdgesAtElementFine(1,IelFine(i))+NVTfine)
        end if

        if(IneighboursAtElementFine(4,IelFine(mod(i,4)+1)).eq.0) then
          DUcoarse(IEdgesAtElementCoarse(i,iel)+NVTcoarse)= &
              DUcoarse(IEdgesAtElementCoarse(i,iel)+NVTcoarse) &
              +0.5*(48.0/64.0)*DUfine(IEdgesAtElementFine(4,IelFine(mod(i,4)+1))+NVTfine)
        end if

      enddo

      ! Collect information from the midpoints of the fine grid elements
      DUcoarse(NVTcoarse+NMTCoarse+iel)= &
          DUcoarse(NVTcoarse+NMTCoarse+iel) &
           +(48.0/64.0)*DUfine(IEdgesAtElementFine(2,IelFine(1))+NVTfine) &
           +(48.0/64.0)*DUfine(IEdgesAtElementFine(2,IelFine(2))+NVTfine) &
           +(48.0/64.0)*DUfine(IEdgesAtElementFine(2,IelFine(3))+NVTfine) &
           +(48.0/64.0)*DUfine(IEdgesAtElementFine(2,IelFine(4))+NVTfine) &
           +(36.0/64.0)*DUfine(NVTfine+NMTfine+IelFine(1)) &
           +(36.0/64.0)*DUfine(NVTfine+NMTfine+IelFine(2)) &
           +(36.0/64.0)*DUfine(NVTfine+NMTfine+IelFine(3)) &
           +(36.0/64.0)*DUfine(NVTfine+NMTfine+IelFine(4))

    enddo

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_interpUniformQ2_double (DuCoarse,DuFine, &
                                           NVTcoarse, NMTcoarse, NELcoarse)

!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! <tex>$Q_2$</tex>, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! Number of vertices in the coarse grid
  integer, intent(in) :: NVTcoarse

  ! Number of edges in the coarse grid
  integer, intent(in) :: NMTcoarse

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables

    ! The first coase.NVT+NMT+NEL entries of the fine grid vector define
    ! the values on the coarse grid - because of the two-level ordering!
    call lalg_copyVector(DUfine(1:NVTcoarse+NMTcoarse+NELcoarse),&
                         DUcoarse(1:NVTCoarse+NMTcoarse+NELcoarse))

  end subroutine

  ! ***************************************************************************
  ! Support for QP1 element
  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_prolUniformQP1_double (DuCoarse,DuFine, &
               IneighboursAtElementFine,NELcoarse,NELfine)

!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! QP1, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Coarse grid vector
  real(DP), dimension(:), intent(in) :: DuCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse

  ! Number of elements in the fine grid
  integer, intent(in) :: NELfine
!</input>

!<output>
  ! Fine grid vector
  real(DP), dimension(:), intent(out) :: DuFine
!</output>

!</subroutine>

  ! In the QP1 element, the DOF`s in the coarse and fine grid are
  ! organised as follows:
  !
  ! +-----------+-----------+
  ! |           |           |
  ! |           |           |
  ! |     X4->      <-X1    |
  ! |     |     ^     |     |
  ! |     v     |     v     |
  ! +--------   O->  -------+
  ! |     ^           ^     |
  ! |     |     |     |     |
  ! |     X1->  |   <-X1    |
  ! |           |           |
  ! |           |           |
  ! +-----------+-----------+
  !
  ! The function value of "O" must be transported by a linear mapping to
  ! all the X. The derivative in "O" is mapped to all "X" as it is constant
  ! in the whole element. This gives the following prolongation matrix
  ! (where O=(f,u,v), Xi = (fi,ui,vi) )
  !
  !   ( 1 -.5 -.5 )  *  ( f )  =  ( f1 )
  !   ( 1  .5 -.5 )     ( u )     ( f2 )
  !   ( 1  .5  .5 )     ( v )     ( f3 )
  !   ( 1 -.5  .5 )               ( f4 )
  !   (     1   0 )               ( u1 )
  !   (     0  -1 )               ( u2 )
  !   (    -1   0 )               ( u3 )
  !   (     1   1 )               ( u4 )
  !   (     0   1 )               ( v1 )
  !   (     1   0 )               ( v2 )
  !   (     0  -1 )               ( v3 )
  !   (    -1   0 )               ( v4 )
  !
  ! The restriction matrix is the transposed of that...

  ! local variables
  real(DP), parameter :: Q2 = 0.5_DP

  integer :: iel,ielh1,ielh2,ielh3,ielh4
  real(DP) :: duh1,duh2,duh3

    ! Loop over the elements
    do iel=1,NELCoarse

      ! Get the values of the basis functions of the coarse grid element
      duh1=DuCoarse(iel)
      duh2=DuCoarse(iel+NELcoarse)
      duh3=DuCoarse(iel+2*NELcoarse)

      ! Get fine grid element numbers
      ielh1=iel
      ielh2=IneighboursAtElementFine(2,ielh1)
      ielh3=IneighboursAtElementFine(2,ielh2)
      ielh4=IneighboursAtElementFine(2,ielh3)

      ! Apply the prolonfation matrix to the coarse grid basis functions
      ! to get the fine grid values.
      DuFine(ielh1) = duh1 - Q2*duh2 - Q2*duh3
      DuFine(ielh2) = duh1 + Q2*duh2 - Q2*duh3
      DuFine(ielh3) = duh1 + Q2*duh2 + Q2*duh3
      DuFine(ielh4) = duh1 - Q2*duh2 + Q2*duh3

      DuFine(NELfine+ielh1) = duh2
      DuFine(NELfine+ielh2) = -duh3
      DuFine(NELfine+ielh3) = -duh2
      DuFine(NELfine+ielh4) = duh3

      DuFine(2*NELfine+ielh1) = duh3
      DuFine(2*NELfine+ielh2) = duh2
      DuFine(2*NELfine+ielh3) = -duh3
      DuFine(2*NELfine+ielh4) = -duh2

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_restUniformQP1_double (DuCoarse,DuFine, &
               IneighboursAtElementFine, NELcoarse,NELfine)

!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! QP1, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IneighboursAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse

  ! Number of elements in the fine grid
  integer, intent(in) :: NELfine
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  real(DP), parameter :: Q2 = 0.5_DP

  integer :: iel
  integer :: ielh1,ielh2,ielh3,ielh4

    ! The information that was "distributed" in the prolongation has to
    ! be "collected". This means, we apply the transposed prolongation
    ! matrix to the RHS vector.

    ! Loop over the elements to collect the additive contributions:
    do iel=1,NELcoarse

      ! Get fine grid element numbers
      ielh1=iel
      ielh2=IneighboursAtElementFine(2,ielh1)
      ielh3=IneighboursAtElementFine(2,ielh2)
      ielh4=IneighboursAtElementFine(2,ielh3)

      ! Collect the distributed values to form the coarse grid RHS:
      DuCoarse(iel) = DuFine(ielh1) + DuFine(ielh2) &
                  + DuFine(ielh3) + DuFine(ielh4)

      DuCoarse(NELcoarse+iel) = &
          -Q2*DuFine(ielh1) + Q2*DuFine(ielh2)      &
            + Q2*DuFine(ielh3) - Q2*DuFine(ielh4)   &
          + DuFine(NELfine+ielh1) - DuFine(NELfine+ielh3)    &
          + DuFine(2*NELfine+ielh2) - DuFine(2*NELfine+ielh4)

      DuCoarse(2*NELcoarse+iel) = &
          -Q2*DuFine(ielh1) - Q2*DuFine(ielh2)      &
            + Q2*DuFine(ielh3) + Q2*DuFine(ielh4)   &
          - DuFine(NELfine+ielh2) + DuFine(NELfine+ielh4)    &
          + DuFine(2*NELfine+ielh1) - DuFine(2*NELfine+ielh3)

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_interpUniformQP1_double (DuCoarse,DuFine, &
                  IneighboursAtElementFine, NELcoarse,NELfine)

!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! QP1, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IneighboursAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse

  ! Number of elements in the fine grid
  integer, intent(in) :: NELfine
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  real(DP), parameter :: Q2 = 0.5_DP
  real(DP), parameter :: Q4 = 0.25_DP

  integer :: iel
  integer :: ielh1,ielh2,ielh3,ielh4

    ! Loop over the elements
    do iel=1,NELcoarse

      ! Get fine grid element numbers
      ielh1=iel
      ielh2=IneighboursAtElementFine(2,ielh1)
      ielh3=IneighboursAtElementFine(2,ielh2)
      ielh4=IneighboursAtElementFine(2,ielh3)

      ! Interpolate the solution on the fine grid to the coarse grid;
      ! take the mean of 4 values each. Take care of the orientation
      ! of the basis functions!

      DuCoarse(iel) = Q4*(DuFine(ielh1) + DuFine(ielh2) &
                         +DuFine(ielh3) + DuFine(ielh4) )

      DuCoarse(NELcoarse+iel) = &
        Q4*(DuFine(ielh1+NELfine) - DuFine(ielh2+2*NELfine) &
           -DuFine(ielh3+NELfine) + DuFine(ielh4+2*NELfine) )

      DuCoarse(2*NELcoarse+iel) = &
        Q4*(DuFine(ielh1+2*NELfine) + DuFine(ielh2+NELfine) &
           -DuFine(ielh3+2*NELfine) - DuFine(ielh4+NELfine) )

    end do

  end subroutine

  ! ***************************************************************************
  ! Support for Q1~ element, DOF`s = integral mean values in the edges
  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_prolUniformEx30_double (DuCoarse,DuFine, &
               IedgesAtElementCoarse,IedgesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               NELcoarse)

!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! E030/EM30, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Coarse grid vector
  real(DP), dimension(:), intent(in) :: DuCoarse

  ! IedgesAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementCoarse

  ! IedgesAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Fine grid vector
  real(DP), dimension(:), intent(out) :: DuFine
!</output>

!</subroutine>

  ! local variables
  real(DP) :: DUH1,DUH2,DUH3,DUH4
  integer :: IM1,IM2, IM3,IM4, IA, IB, IC
  integer :: iel, IELH1, IELH2, IELH3, IELH4

  ! Weights for the restriction; all coefficients are halfed, so dividing
  ! by 2 is not necessary in the calculation routines.
  real(DP), parameter :: A1=0.5_DP, A2=-0.125_DP, A3=0.0_DP, A4=0.125_DP
  real(DP), parameter :: A5=0.625_DP, A6=0.125_DP, A7=0.125_DP, A8=0.125_DP

    ! Clear the output vector
    call lalg_clearVector(DuFine)

    ! Loop over the coarse grid elements
    do iel=1,NELcoarse

      ! Get the DOF`s of the coarse grid element
      IM1 = IedgesAtElementCoarse(1,iel)
      IM2 = IedgesAtElementCoarse(2,iel)
      IM3 = IedgesAtElementCoarse(3,iel)
      IM4 = IedgesAtElementCoarse(4,iel)

      ! Get the values of the corresponding DOF`s
      DUH1 = DuCoarse(IM1)
      DUH2 = DuCoarse(IM2)
      DUH3 = DuCoarse(IM3)
      DUH4 = DuCoarse(IM4)

      ! Get the element numbers of the fine-grid elements in the coarse grid;
      ! the numbers are defined by the two-level ordering.
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(2,IELH3)

      ! Distribute the value at the edge IM1 to the
      ! corresponding fine inner nodes

      if (IneighboursAtElementCoarse(1,iel).ne.0) then
        ! There is a neighbour at the edge
        IA=IedgesAtElementFine(1,IELH1)
        IB=IedgesAtElementFine(4,IELH2)
        DuFine(IA)=DuFine(IA)+   A1*DUH1+A2*DUH2+A3*DUH3+A4*DUH4
        DuFine(IB)=DuFine(IB)+   A1*DUH1+A4*DUH2+A3*DUH3+A2*DUH4
      else
        ! No neighbour; boundary element
        IA=IedgesAtElementFine(1,IELH1)
        IB=IedgesAtElementFine(4,IELH2)
        DuFine(IA)=DuFine(IA)+2.0_DP*(A1*DUH1+A2*DUH2+A3*DUH3+A4*DUH4)
        DuFine(IB)=DuFine(IB)+2.0_DP*(A1*DUH1+A4*DUH2+A3*DUH3+A2*DUH4)
      endif
      IC=IedgesAtElementFine(2,IELH1)
      DuFine(IC)=A5*DUH1+A6*(DUH2+DUH4)+A7*DUH3

      ! Distribute the value at the edge IM2 to the
      ! corresponding fine inner nodes

      if (IneighboursAtElementCoarse(2,iel).ne.0) then
        ! There is a neighbour at the edge
       IA=IedgesAtElementFine(1,IELH2)
       IB=IedgesAtElementFine(4,IELH3)
       DuFine(IA)=DuFine(IA)+   A1*DUH2+A2*DUH3+A3*DUH4+A4*DUH1
       DuFine(IB)=DuFine(IB)+   A1*DUH2+A4*DUH3+A3*DUH4+A2*DUH1
      else
        ! No neighbour; boundary element
       IA=IedgesAtElementFine(1,IELH2)
       IB=IedgesAtElementFine(4,IELH3)
       DuFine(IA)=DuFine(IA)+2.0_DP*(A1*DUH2+A2*DUH3+A3*DUH4+A4*DUH1)
       DuFine(IB)=DuFine(IB)+2.0_DP*(A1*DUH2+A4*DUH3+A3*DUH4+A2*DUH1)
      endif
      IC=IedgesAtElementFine(2,IELH2)
      DuFine(IC)=A5*DUH2+A6*(DUH3+DUH1)+A7*DUH4

      ! Distribute the value at the edge IM3 to the
      ! corresponding fine inner nodes

      if (IneighboursAtElementCoarse(3,iel).ne.0) then
        ! There is a neighbour at the edge
       IA=IedgesAtElementFine(1,IELH3)
       IB=IedgesAtElementFine(4,IELH4)
       DuFine(IA)=DuFine(IA)+   A1*DUH3+A2*DUH4+A3*DUH1+A4*DUH2
       DuFine(IB)=DuFine(IB)+   A1*DUH3+A4*DUH4+A3*DUH1+A2*DUH2
      else
        ! No neighbour; boundary element
       IA=IedgesAtElementFine(1,IELH3)
       IB=IedgesAtElementFine(4,IELH4)
       DuFine(IA)=DuFine(IA)+2.0_DP*(A1*DUH3+A2*DUH4+A3*DUH1+A4*DUH2)
       DuFine(IB)=DuFine(IB)+2.0_DP*(A1*DUH3+A4*DUH4+A3*DUH1+A2*DUH2)
      endif
      IC=IedgesAtElementFine(2,IELH3)
      DuFine(IC)=A5*DUH3+A6*(DUH4+DUH2)+A7*DUH1

      ! Distribute the value at the edge IM4 to the
      ! corresponding fine inner nodes

      if (IneighboursAtElementCoarse(4,iel).ne.0) then
        ! There is a neighbour at the edge
        IA=IedgesAtElementFine(1,IELH4)
        IB=IedgesAtElementFine(4,IELH1)
        DuFine(IA)=DuFine(IA)+   A1*DUH4+A2*DUH1+A3*DUH2+A4*DUH3
        DuFine(IB)=DuFine(IB)+   A1*DUH4+A4*DUH1+A3*DUH2+A2*DUH3
      else
        ! No neighbour; boundary element
        IA=IedgesAtElementFine(1,IELH4)
        IB=IedgesAtElementFine(4,IELH1)
        DuFine(IA)=DuFine(IA)+2.0_DP*(A1*DUH4+A2*DUH1+A3*DUH2+A4*DUH3)
        DuFine(IB)=DuFine(IB)+2.0_DP*(A1*DUH4+A4*DUH1+A3*DUH2+A2*DUH3)
      end if
      IC=IedgesAtElementFine(2,IELH4)
      DuFine(IC)=A5*DUH4+A6*(DUH1+DUH3)+A7*DUH2

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_prolUniformEx30ext_double (DuCoarse,DuFine, &
               DvertexCoordsCoarse,IverticesAtElementCoarse,DelementAreaCoarse,&
               IedgesAtElementCoarse,IedgesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               NELcoarse, iweightingType, daspectRatioBound, iarIndicator)

!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! E030/EM30, uniform triangulation, double precision vector.
  !
  ! Extended version. Switch to constant prolongation if aspect ratio
  ! of an element is too large.
!</description>

!<input>
  ! Coarse grid vector
  real(DP), dimension(:), intent(in) :: DuCoarse

  ! DvertexCoords array on the coarse grid
  real(DP), dimension(:,:), intent(in)                :: DvertexCoordsCoarse

  ! IverticesAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IverticesAtElementCoarse

  ! DelementArea array on the coarse grid
  real(DP), dimension(:), intent(in)                  :: DelementAreaCoarse

  ! IedgesAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementCoarse

  ! IedgesAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse

  ! Type of the averaging on the element edges
  ! <=0: standard averaging of both contributions by 1/2,
  !  =1: weighted averaging of the interpolated function values:
  !      The area of the current coarse grid element determines
  !      the weight. (L2-projection, standard),
  !  =2: weighted averaging of the interpolated function values:
  !      The area of the neightbour element of the coarse grid
  !      the weight.
  integer, intent(in)  :: iweightingType

  ! Upper bound aspect ratio; for all elements with higher AR
  ! the prolongation is switched to constant prolongation
  real(DP), intent(in) :: daspectRatioBound

  ! Aspect-ratio indicator.
  ! Controls switching to constant prolongation.
  ! <=1: switch depending on aspect ratio of current element,
  !  =2: switch depending on aspect ratio of current element and
  !      neighbour element
  integer, intent(in)  :: iarIndicator
!</input>

!<output>
  ! Fine grid vector
  real(DP), dimension(:), intent(out) :: DuFine
!</output>

!</subroutine>

  ! local variables
  real(DP) :: DUH1,DUH2,DUH3,DUH4
  real(DP), dimension(0:TRIA_MAXNME2D) :: daspectRatio,darea
  real(DP), dimension(TRIA_MAXNME2D) :: dweight
  integer :: IM1,IM2, IM3,IM4, IA, IB, IC
  integer :: iel, IELH1, IELH2, IELH3, IELH4
  integer, dimension(0:TRIA_MAXNME2D) :: IELA
  integer :: i
  integer, dimension(TRIA_MAXNME2D) :: idoConstant
  real(DP), dimension(NDIM2D,TRIA_MAXNVE2D) :: dcoords

  ! Weights for the prolongation.
  ! PRWEIG (.,1) gives the constants for the standard prolongation,
  ! PRWEIG (.,2) gives the constants for the constant prolongation.
  real(DP), dimension(8,2), parameter :: prweight = &
      reshape((/1.0_DP, -0.25_DP, 0.0_DP, 0.25_DP, &
                0.625_DP, 0.125_DP, 0.125_DP, 0.125_DP, &
                1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, &
                1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP/),(/8,2/))

    ! Clear the output vector
    call lalg_clearVector(DuFine)

    ! Loop over the coarse grid elements
    do iel=1,NELcoarse

      ! Get the numbers of the elements that are neighbours to our current coarse
      ! grid element:
      !             +--------+
      !             |        |
      !             | IELA3  |
      !             |        |
      !    +--------4--------3--------+
      !    |        |        |        |
      !    | IELA4  |  IEL   | IELA2  |
      !    |        |        |        |
      !    +--------1--------2--------+
      !             |        |
      !             | IELA1  |
      !             |        |
      !             +--------+
      IELA(0) = iel
      IELA(1) = IneighboursAtElementCoarse(1,iel)
      IELA(2) = IneighboursAtElementCoarse(2,iel)
      IELA(3) = IneighboursAtElementCoarse(3,iel)
      IELA(4) = IneighboursAtElementCoarse(4,iel)

      ! For these five elements, determine the aspect ratio and their area.
      !
      ! At first the element in the center, which always exists.
      !
      ! Get the aspect ratio of the current coarse grid element;
      ! if necessary, calculate the reciprocal.
      dcoords = DvertexCoordsCoarse(:,IverticesAtElementCoarse(:,IELA(0)))
      daspectRatio(0) = gaux_getAspectRatio_quad2D (dcoords)
      if (daspectRatio(0) .lt. 1.0_DP) daspectRatio(0) = 1.0_DP/daspectRatio(0)

      ! and the area of that element.
      darea(0) = DelementAreaCoarse(iel)

      ! Then the remaining neighbours.
      do i=1,TRIA_MAXNME2D
        if (IELA(i) .ne. 0) then
          ! Get the aspect ratio of the current coarse grid element;
          ! if necessary, calculate the reciprocal.
          dcoords = DvertexCoordsCoarse(:,IverticesAtElementCoarse(:,IELA(i)))
          daspectRatio(i) = gaux_getAspectRatio_quad2D (dcoords)
          if (daspectRatio(i) .lt. 1.0_DP) daspectRatio(i) = 1.0_DP/daspectRatio(i)

          ! and the area of that element.
          darea(i) = DelementAreaCoarse(IELA(i))
        else
          daspectRatio(i) = 0.0_DP
          darea(i) = 0.0_DP
        end if
      end do

      ! Calculate weighting factors for the interpolation.
      ! The iweightingType parameter describes:
      ! <= 0: simple interpolation, weight both contribtutions by 1/2
      !  = 1: take the weighted mean of the interpolated function values
      !       by weighting with the area of the current coarse grid element
      ! >= 2: take the weighted mean of the interpolated function values
      !       by weighting with the area of the neighboured coarse grid element

      select case (iweightingType)
      case (:0)
        dweight = 0.5_DP
      case (1)
        dweight = darea(0) / (darea(0)+darea(1:TRIA_MAXNME2D))
      case (2:)
        dweight = darea(1:TRIA_MAXNME2D) / (darea(0)+darea(1:TRIA_MAXNME2D))
      end select

      ! Where there is no neighbour, set the weighting factor to 1.0
!      DO i=1,TRIA_MAXNMT
!        IF (IELA(i) .EQ. 0) dweight(i) = 1.0_DP
!      END DO
      where(IELA(1:TRIA_MAXNME2D) .eq. 0) dweight(1:TRIA_MAXNME2D) = 1.0_DP

      ! Now determine on which edge to switch to constant prolongation
      ! By default, we do not use constant prolongation
      idoConstant = 1

      ! ... but if the caller wants us to switch in a special situation...
      if ((iarIndicator .ge. 1) .and. (daspectRatioBound .ge. 0.0_DP)) then

        ! ... switch to constant of our element is too large...
        if (daspectRatio(0) .gt. daspectRatioBound) idoConstant = 2

        ! and if iarIndicator>2, also check the neighbour element
        if (iarIndicator .ge. 2) then
!          DO i=1,TRIA_MAXNME2D
!            IF (daspectRatio(i) .GT. daspectRatioBound) idoConstant(i) = 2
!          END DO
          where (daspectRatio(1:4) .gt. daspectRatioBound) idoConstant = 2
        end if

      end if

      ! Get the DOF`s of the coarse grid element
      IM1 = IedgesAtElementCoarse(1,iel)
      IM2 = IedgesAtElementCoarse(2,iel)
      IM3 = IedgesAtElementCoarse(3,iel)
      IM4 = IedgesAtElementCoarse(4,iel)

      ! Get the values of the corresponding DOF`s
      DUH1 = DuCoarse(IM1)
      DUH2 = DuCoarse(IM2)
      DUH3 = DuCoarse(IM3)
      DUH4 = DuCoarse(IM4)

      ! Get the element numbers of the fine-grid elements in the coarse grid;
      ! the numbers are defined by the two-level ordering.
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(2,IELH3)

      ! Now let us start with the actual prolongation
      ! ---------------------------------------------

      ! Get the DOF`s on the fine grid
      IA=IedgesAtElementFine(1,IELH1)
      IB=IedgesAtElementFine(4,IELH2)
      IC=IedgesAtElementFine(2,IELH1)

      ! Now we have the following situation:

      !   4               IM3                3
      !     ===============X================
      !     |              |               |
      !     |              |               |
      !     |    IELH4     |     IELH3     |
      !     |              |               |
      !     |                              |
      ! IM4 X----------- IEL1 -------------X IM2
      !     |                              |
      !     |              |               |
      !     |    IELH1     o IC  IELH2     |
      !     |              |               |
      !     |              |               |
      !   1 =======o=======X=======o======== 2
      !     |     IA      IM1      IB      |
      !     |                              |
      !     |                              |
      !     |                              |
      !     |                              |
      !     |            IELA1             |
      !     |                              |
      !     |                              |
      !     |                              |
      !     |                              |
      !     |                              |
      !     ================================

      ! Distribute the value at the edge IM1 to the
      ! corresponding fine inner nodes

      DuFine(IA) = DuFine(IA) &
                  + dweight(1)*(prweight(1,idoConstant(1))*DUH1 &
                              +prweight(2,idoConstant(1))*DUH2 &
                              +prweight(3,idoConstant(1))*DUH3 &
                              +prweight(4,idoConstant(1))*DUH4)
      DuFine(IB) = DuFine(IB) &
                  + dweight(1)*(prweight(1,idoConstant(1))*DUH1 &
                              +prweight(4,idoConstant(1))*DUH2 &
                              +prweight(3,idoConstant(1))*DUH3 &
                              +prweight(2,idoConstant(1))*DUH4)
      DuFine(IC) = prweight(5,idoConstant(1))*DUH1 &
                 + prweight(6,idoConstant(1))*(DUH2+DUH4) &
                 + prweight(7,idoConstant(1))*DUH3

      ! Distribute the value at the edge IM2 to the
      ! corresponding fine inner nodes
      IA=IedgesAtElementFine(1,IELH2)
      IB=IedgesAtElementFine(4,IELH3)
      IC=IedgesAtElementFine(2,IELH2)

      DuFine(IA) = DuFine(IA) &
                 + dweight(2)*(prweight(1,idoConstant(2))*DUH2 &
                              +prweight(2,idoConstant(2))*DUH3 &
                              +prweight(3,idoConstant(2))*DUH4 &
                              +prweight(4,idoConstant(2))*DUH1)
      DuFine(IB) = DuFine(IB) &
                 + dweight(2)*(prweight(1,idoConstant(2))*DUH2 &
                              +prweight(4,idoConstant(2))*DUH3 &
                              +prweight(3,idoConstant(2))*DUH4 &
                              +prweight(2,idoConstant(2))*DUH1)
      DuFine(IC) = prweight(5,idoConstant(2))*DUH2 &
                 + prweight(6,idoConstant(2))*(DUH3+DUH1) &
                 + prweight(7,idoConstant(2))*DUH4

      ! Distribute the value at the edge IM3 to the
      ! corresponding fine inner nodes
      IA=IedgesAtElementFine(1,IELH3)
      IB=IedgesAtElementFine(4,IELH4)
      IC=IedgesAtElementFine(2,IELH3)

      DuFine(IA) = DuFine(IA) &
                 + dweight(3)*(prweight(1,idoConstant(3))*DUH3 &
                              +prweight(2,idoConstant(3))*DUH4 &
                              +prweight(3,idoConstant(3))*DUH1 &
                              +prweight(4,idoConstant(3))*DUH2)
      DuFine(IB) = DuFine(IB) &
                 + dweight(3)*(prweight(1,idoConstant(3))*DUH3 &
                              +prweight(4,idoConstant(3))*DUH4 &
                              +prweight(3,idoConstant(3))*DUH1 &
                              +prweight(2,idoConstant(3))*DUH2)
      DuFine(IC) = prweight(5,idoConstant(3))*DUH3 &
                 + prweight(6,idoConstant(3))*(DUH4+DUH2) &
                 + prweight(7,idoConstant(3))*DUH1

      ! Distribute the value at the edge IM4 to the
      ! corresponding fine inner nodes
      IA=IedgesAtElementFine(1,IELH4)
      IB=IedgesAtElementFine(4,IELH1)
      IC=IedgesAtElementFine(2,IELH4)

      DuFine(IA) = DuFine(IA) &
                 + dweight(4)*(prweight(1,idoConstant(4))*DUH4 &
                              +prweight(2,idoConstant(4))*DUH1 &
                              +prweight(3,idoConstant(4))*DUH2 &
                              +prweight(4,idoConstant(4))*DUH3)
      DuFine(IB) = DuFine(IB) &
                 + dweight(4)*(prweight(1,idoConstant(4))*DUH4 &
                              +prweight(4,idoConstant(4))*DUH1 &
                              +prweight(3,idoConstant(4))*DUH2 &
                              +prweight(2,idoConstant(4))*DUH3)
      DuFine(IC) = prweight(5,idoConstant(4))*DUH4 &
                 + prweight(6,idoConstant(4))*(DUH1+DUH3) &
                 + prweight(7,idoConstant(4))*DUH2

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_restUniformEx30_double (DuCoarse,DuFine, &
               IedgesAtElementCoarse,IedgesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               NELcoarse)

!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! E030/EM30 element, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IedgesAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementCoarse

  ! IedgesAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

    ! local variables
    real(DP) :: DUH1,DUH2,DUH3,DUH4,DUH5,DUH6,DUH7,DUH8,DUH9,DUH10,DUH11,DUH12
    integer :: IM1,IM2, IM3,IM4
    integer :: I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12
    integer :: iel, IELH1, IELH2, IELH3, IELH4

    ! Weights for the restriction
    real(DP), parameter :: A1=1.0_DP, A2=-0.125_DP, A3=0.0_DP, A4=0.125_DP
    real(DP), parameter :: A5=0.625_DP, A6=0.125_DP, A7=0.125_DP, A8=0.125_DP

    ! Loop over the coarse grid elements
    do iel=1,NELcoarse

      ! Get the DOF`s of the coarse grid element
      IM1 = IedgesAtElementCoarse(1,iel)
      IM2 = IedgesAtElementCoarse(2,iel)
      IM3 = IedgesAtElementCoarse(3,iel)
      IM4 = IedgesAtElementCoarse(4,iel)

      ! Get the element numbers of the fine-grid elements in the coarse grid;
      ! the numbers are defined by the two-level ordering.
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(2,IELH3)

      ! Get the DOF`s on the fine grid
      I1=IedgesAtElementFine(1,IELH1)
      I2=IedgesAtElementFine(4,IELH2)
      I3=IedgesAtElementFine(1,IELH2)
      I4=IedgesAtElementFine(4,IELH3)
      I5=IedgesAtElementFine(1,IELH3)
      I6=IedgesAtElementFine(4,IELH4)
      I7=IedgesAtElementFine(1,IELH4)
      I8=IedgesAtElementFine(4,IELH1)
      I9=IedgesAtElementFine(2,IELH1)
      I10=IedgesAtElementFine(2,IELH2)
      I11=IedgesAtElementFine(2,IELH3)
      I12=IedgesAtElementFine(2,IELH4)

      ! Get the values of the DOF`s on the fine grid
      DUH1= DuFine(I1)
      DUH2= DuFine(I2)
      DUH3= DuFine(I3)
      DUH4= DuFine(I4)
      DUH5= DuFine(I5)
      DUH6= DuFine(I6)
      DUH7= DuFine(I7)
      DUH8= DuFine(I8)
      DUH9= DuFine(I9)
      DUH10=DuFine(I10)
      DUH11=DuFine(I11)
      DUH12=DuFine(I12)

      ! Now collect the information in the same way as it was
      ! "distributed" by the prolongation routine.
      ! This realises the adjoint operator of the prolongation.
      !
      ! Calculate the value of the edge IM1

      if (IneighboursAtElementCoarse(1,iel).ne.0) then
        ! inner edge
        if (IneighboursAtElementCoarse(1,iel).gt.iel) then
          DuCoarse(IM1)= A1*(DUH1+DUH2)+A2*(DUH4+DUH7) &
                       + A3*(DUH5+DUH6)+A4*(DUH3+DUH8) &
                       + A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
        else
          DuCoarse(IM1)= DuCoarse(IM1)+A2*(DUH4+DUH7) &
                       + A3*(DUH5+DUH6)+A4*(DUH3+DUH8) &
                       + A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
        endif
      else
        ! boundary edge
        DuCoarse(IM1)=     A1*(DUH1+DUH2)+2.0_DP*A2*(DUH4+DUH7) &
                +2.0_DP*A3*(DUH5+DUH6)+2.0_DP*A4*(DUH3+DUH8) &
                +       A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
      endif

      ! Calculate the value of the edge IM2

      if (IneighboursAtElementCoarse(2,iel).ne.0) then
        ! inner edge
        if (IneighboursAtElementCoarse(2,iel).gt.iel) then
           DuCoarse(IM2)= A1*(DUH3+DUH4)+A2*(DUH6+DUH1) &
                         +A3*(DUH7+DUH8)+A4*(DUH5+DUH2) &
                         +A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
        else
          DuCoarse(IM2)=DuCoarse(IM2)+A2*(DUH6+DUH1) &
                              +A3*(DUH7+DUH8)+A4*(DUH5+DUH2) &
                              +A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
        endif
      else
        ! boundary edge
        DuCoarse(IM2)= A1*(DUH3+DUH4)+2.0_DP*A2*(DUH6+DUH1) &
                +2.0_DP*A3*(DUH7+DUH8)+2.0_DP*A4*(DUH5+DUH2) &
                +       A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
      endif

      ! Calculate the value of the edge IM3

      if (IneighboursAtElementCoarse(3,iel).ne.0) then
        ! inner edge
        if (IneighboursAtElementCoarse(3,iel).gt.iel) then
           DuCoarse(IM3)= A1*(DUH5+DUH6)+A2*(DUH8+DUH3) &
                        + A3*(DUH1+DUH2)+A4*(DUH7+DUH4) &
                        + A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
        else
           DuCoarse(IM3)= DuCoarse(IM3)+A2*(DUH8+DUH3) &
                        + A3*(DUH1+DUH2)+A4*(DUH7+DUH4) &
                        + A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
        endif
      else
        ! boundary edge
        DuCoarse(IM3)= A1*(DUH5+DUH6)+2.0_DP*A2*(DUH8+DUH3) &
                +2.0_DP*A3*(DUH1+DUH2)+2.0_DP*A4*(DUH7+DUH4) &
                +       A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
      endif

      ! Calculate the value of the edge IM4

      if (IneighboursAtElementCoarse(4,iel).ne.0) then
        ! inner edge
        if (IneighboursAtElementCoarse(4,iel).gt.iel) then
          DuCoarse(IM4)= A1*(DUH7+DUH8)+A2*(DUH2+DUH5) &
                       + A3*(DUH3+DUH4)+A4*(DUH1+DUH6) &
                       + A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
        else
          DuCoarse(IM4)= DuCoarse(IM4)+A2*(DUH2+DUH5) &
                       + A3*(DUH3+DUH4)+A4*(DUH1+DUH6) &
                       + A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
        endif
      else
        ! boundary edge
        DuCoarse(IM4)= A1*(DUH7+DUH8)+2.0_DP*A2*(DUH2+DUH5) &
                +2.0_DP*A3*(DUH3+DUH4)+2.0_DP*A4*(DUH1+DUH6) &
                +       A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
      endif

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_restUniformEx30ext_double (DuCoarse,DuFine, &
               DvertexCoordsCoarse,IverticesAtElementCoarse,DelementAreaCoarse,&
               IedgesAtElementCoarse,IedgesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               NELcoarse, iweightingType, daspectRatioBound, iarIndicator)

!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! E030/EM30 element, uniform triangulation, double precision vector.
  !
  ! Extended version. Switch to constant restriction if aspect ratio
  ! of an element is too large.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! DvertexCoords array on the coarse grid
  real(DP), dimension(:,:), intent(in)                :: DvertexCoordsCoarse

  ! IverticesAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IverticesAtElementCoarse

  ! DelementArea array on the coarse grid
  real(DP), dimension(:), intent(in)                  :: DelementAreaCoarse

  ! IedgesAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementCoarse

  ! IedgesAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse

  ! Type of the averaging on the element edges
  ! <=0: standard averaging of both contributions by 1/2,
  !  =1: weighted averaging of the interpolated function values:
  !      The area of the current coarse grid element determines
  !      the weight. (L2-projection, standard),
  !  =2: weighted averaging of the interpolated function values:
  !      The area of the neightbour element of the coarse grid
  !      the weight.
  integer, intent(in)  :: iweightingType

  ! Upper bound aspect ratio; for all elements with higher AR
  ! the prolongation is switched to constant prolongation
  real(DP), intent(in) :: daspectRatioBound

  ! Aspect-ratio indicator.
  ! Controls switching to constant prolongation.
  ! <=1: switch depending on aspect ratio of current element,
  !  =2: switch depending on aspect ratio of current element and
  !      neighbour element
  integer, intent(in)  :: iarIndicator
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

    ! local variables
    real(DP) :: DUH1,DUH2,DUH3,DUH4,DUH5,DUH6,DUH7,DUH8,DUH9,DUH10,DUH11,DUH12
    integer :: IM1,IM2, IM3,IM4
    integer :: I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12
    integer :: iel, IELH1, IELH2, IELH3, IELH4
    integer, dimension(0:TRIA_MAXNME2D) :: IELA
    integer :: i
    integer, dimension(TRIA_MAXNME2D) :: idoConstant
    real(DP), dimension(0:TRIA_MAXNME2D) :: daspectRatio, darea
    real(DP), dimension(TRIA_MAXNME2D) :: dweight
    real(DP), dimension(NDIM2D,TRIA_MAXNVE2D) :: dcoords

    ! Weights for the restriction.
    ! PRWEIG (.,1) gives the constants for the standard restriction,
    ! PRWEIG (.,2) gives the constants for the constant restriction.
    real(DP), dimension(8,2), parameter :: prweight = &
        reshape((/1.0_DP, -0.25_DP, 0.0_DP, 0.25_DP, &
                  0.625_DP, 0.125_DP, 0.125_DP, 0.125_DP, &
                  1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, &
                  1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP/),(/8,2/))

    ! Clear the output vector
    call lalg_clearVector(DuCoarse)

    ! Loop over the coarse grid elements
    do iel=1,NELcoarse

      ! Get the numbers of the elements that are neighbours to our current coarse
      ! grid element:
      !             +--------+
      !             |        |
      !             | IELA3  |
      !             |        |
      !    +--------4--------3--------+
      !    |        |        |        |
      !    | IELA4  |  IEL   | IELA2  |
      !    |        |        |        |
      !    +--------1--------2--------+
      !             |        |
      !             | IELA1  |
      !             |        |
      !             +--------+
      IELA(0) = iel
      IELA(1) = IneighboursAtElementCoarse(1,iel)
      IELA(2) = IneighboursAtElementCoarse(2,iel)
      IELA(3) = IneighboursAtElementCoarse(3,iel)
      IELA(4) = IneighboursAtElementCoarse(4,iel)

      ! For these five elements, determine the aspect ratio and their area.
      !
      ! At first the element in the center, which always exists.
      !
      ! Get the aspect ratio of the current coarse grid element;
      ! if necessary, calculate the reciprocal.
      dcoords = DvertexCoordsCoarse(:,IverticesAtElementCoarse(:,IELA(0)))
      daspectRatio(0) = gaux_getAspectRatio_quad2D (dcoords)
      if (daspectRatio(0) .lt. 1.0_DP) daspectRatio(0) = 1.0_DP/daspectRatio(0)

      ! and the area of that element.
      darea(0) = DelementAreaCoarse(iel)

      ! Then the remaining neighbours.
      do i=1,TRIA_MAXNME2D
        if (IELA(i) .ne. 0) then
          ! Get the aspect ratio of the current coarse grid element;
          ! if necessary, calculate the reciprocal.
          dcoords = DvertexCoordsCoarse(:,IverticesAtElementCoarse(:,IELA(i)))
          daspectRatio(i) = gaux_getAspectRatio_quad2D (dcoords)
          if (daspectRatio(i) .lt. 1.0_DP) daspectRatio(i) = 1.0_DP/daspectRatio(i)

          ! and the area of that element.
          darea(i) = DelementAreaCoarse(IELA(i))
        else
          daspectRatio(i) = 0.0_DP
          darea(i) = 0.0_DP
        end if
      end do

      ! Calculate weighting factors for the interpolation.
      ! The iweightingType parameter describes:
      ! <= 0: simple interpolation, weight both contribtutions by 1/2
      !  = 1: take the weighted mean of the interpolated function values
      !       by weighting with the area of the current coarse grid element
      ! >= 2: take the weighted mean of the interpolated function values
      !       by weighting with the area of the neighboured coarse grid element

      select case (iweightingType)
      case (:0)
        dweight = 0.5_DP
      case (1)
        dweight = darea(0) / (darea(0)+darea(1:TRIA_MAXNME2D))
      case (2:)
        dweight = darea(1:TRIA_MAXNME2D) / (darea(0)+darea(1:TRIA_MAXNME2D))
      end select

      ! Where there is no neighbour, set the weighting factor to 1.0
!      DO i=1,TRIA_MAXNMT
!        IF (IELA(i) .EQ. 0) dweight(i) = 1.0_DP
!      END DO
      where(IELA(1:TRIA_MAXNME2D) .eq. 0) dweight(1:TRIA_MAXNME2D) = 1.0_DP

      ! Now determine on which edge to switch to constant prolongation
      ! By default, we do not use constant prolongation
      idoConstant = 1

      ! ... but if the caller wants us to switch in a special situation...
      if ((iarIndicator .ge. 1) .and. (daspectRatioBound .ge. 0.0_DP)) then

        ! ... switch to constant of our element is too large...
        if (daspectRatio(0) .gt. daspectRatioBound) idoConstant = 2

        ! and if iarIndicator>2, also check the neighbour element
        if (iarIndicator .ge. 2) then
!          DO i=1,TRIA_MAXNME2D
!            IF (daspectRatio(i) .GT. daspectRatioBound) idoConstant(i) = 2
!          END DO
          where (daspectRatio(1:4) .gt. daspectRatioBound) idoConstant = 2
        end if

      end if

      ! Now let us strt with the actual restriction
      ! ------------------------------------------

      ! Get the DOF`s of the coarse grid element
      IM1 = IedgesAtElementCoarse(1,iel)
      IM2 = IedgesAtElementCoarse(2,iel)
      IM3 = IedgesAtElementCoarse(3,iel)
      IM4 = IedgesAtElementCoarse(4,iel)

      ! Get the element numbers of the fine-grid elements in the coarse grid;
      ! the numbers are defined by the two-level ordering.
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(2,IELH3)

      ! Get the DOF`s on the fine grid
      I1=IedgesAtElementFine(1,IELH1)
      I2=IedgesAtElementFine(4,IELH2)
      I3=IedgesAtElementFine(1,IELH2)
      I4=IedgesAtElementFine(4,IELH3)
      I5=IedgesAtElementFine(1,IELH3)
      I6=IedgesAtElementFine(4,IELH4)
      I7=IedgesAtElementFine(1,IELH4)
      I8=IedgesAtElementFine(4,IELH1)
      I9=IedgesAtElementFine(2,IELH1)
      I10=IedgesAtElementFine(2,IELH2)
      I11=IedgesAtElementFine(2,IELH3)
      I12=IedgesAtElementFine(2,IELH4)

      ! Get the values of the DOF`s on the fine grid
      DUH1= DuFine(I1)
      DUH2= DuFine(I2)
      DUH3= DuFine(I3)
      DUH4= DuFine(I4)
      DUH5= DuFine(I5)
      DUH6= DuFine(I6)
      DUH7= DuFine(I7)
      DUH8= DuFine(I8)
      DUH9= DuFine(I9)
      DUH10=DuFine(I10)
      DUH11=DuFine(I11)
      DUH12=DuFine(I12)

      ! Now we have the following situation:
      !   4       I6      IM3      I5      3
      !     =======o=======X=======o========
      !     |              |               |
      !     |              |               |
      !  I7 o    IELH4     o I11 IELH3     o I4
      !     |              |               |
      !     |                              |
      ! IM4 X------o---- IEL1 -----o-------X IM2
      !     |    I12               I10     |
      !     |              |               |
      !  I8 o    IELH1     o I9  IELH2     o I3
      !     |              |               |
      !     |              |               |
      !   1 =======o=======X=======o======== 2
      !     |     I1      IM1      I2      |
      !     |                              |
      !     |                              |
      !     |                              |
      !     |                              |
      !     |            IELA1             |
      !     |                              |
      !     |                              |
      !     |                              |
      !     |                              |
      !     |                              |
      !     ================================


      ! Collect the information in the same way as it was
      ! "distributed" by the prolongation routine.
      ! This realises the adjoint operator of the prolongation.
      !
      ! Calculate the value of the edge IM1

      DuCoarse(IM1)= DuCoarse(IM1) &
                    +dweight(1)*(prweight(1,idoConstant(1))*(DUH1+DUH2) &
                                +prweight(2,idoConstant(1))*(DUH4+DUH7) &
                                +prweight(3,idoConstant(1))*(DUH5+DUH6) &
                                +prweight(4,idoConstant(1))*(DUH3+DUH8)) &
                    +prweight(5,idoConstant(1))*DUH9 &
                    +prweight(6,idoConstant(1))*(DUH10+DUH12) &
                    +prweight(7,idoConstant(1))*DUH11

      ! Calculate the value of the edge IM2

      DuCoarse(IM2)= DuCoarse(IM2) &
                    +dweight(2)*(prweight(1,idoConstant(2))*(DUH3+DUH4) &
                                +prweight(2,idoConstant(2))*(DUH6+DUH1) &
                                +prweight(3,idoConstant(2))*(DUH7+DUH8) &
                                +prweight(4,idoConstant(2))*(DUH5+DUH2)) &
                    +prweight(5,idoConstant(2))*DUH10 &
                    +prweight(6,idoConstant(2))*(DUH11+DUH9) &
                    +prweight(7,idoConstant(2))*DUH12

      ! Calculate the value of the edge IM3

      DuCoarse(IM3)= DuCoarse(IM3) &
                    +dweight(3)*(prweight(1,idoConstant(3))*(DUH5+DUH6) &
                                +prweight(2,idoConstant(3))*(DUH8+DUH3) &
                                +prweight(3,idoConstant(3))*(DUH1+DUH2) &
                                +prweight(4,idoConstant(3))*(DUH7+DUH4)) &
                    +prweight(5,idoConstant(3))*DUH11 &
                    +prweight(6,idoConstant(3))*(DUH12+DUH10) &
                    +prweight(7,idoConstant(3))*DUH9

      ! Calculate the value of the edge IM4

      DuCoarse(IM4)= DuCoarse(IM4) &
                    +dweight(4)*(prweight(1,idoConstant(4))*(DUH7+DUH8) &
                                +prweight(2,idoConstant(4))*(DUH2+DUH5) &
                                +prweight(3,idoConstant(4))*(DUH3+DUH4) &
                                +prweight(4,idoConstant(4))*(DUH1+DUH6)) &
                    +prweight(5,idoConstant(4))*DUH12 &
                    +prweight(6,idoConstant(4))*(DUH9+DUH11) &
                    +prweight(7,idoConstant(4))*DUH10

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_interpUniformEx30_double (DuCoarse,DuFine, &
               IedgesAtElementCoarse,IedgesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               NELcoarse)

!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! Ex30, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IedgesAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementCoarse

  ! IedgesAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

    ! local variables
    real(DP) :: DUH1,DUH2,DUH3,DUH4,DUH5,DUH6,DUH7,DUH8,DUH9,DUH10,DUH11,DUH12
    integer :: IM1,IM2, IM3,IM4
    integer :: I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12
    integer :: iel, IELH1, IELH2, IELH3, IELH4

    ! Weights for the restriction
    real(DP), parameter :: A1=0.1875_DP, A2=0.375_DP, A3=-0.0625_DP
    real(DP), parameter :: R1=0.375_DP, R2=0.75_DP, R3=-0.125_DP

    ! Loop over the coarse grid elements
    do iel=1,NELcoarse

      ! Get the DOF`s of the coarse grid element
      IM1 = IedgesAtElementCoarse(1,iel)
      IM2 = IedgesAtElementCoarse(2,iel)
      IM3 = IedgesAtElementCoarse(3,iel)
      IM4 = IedgesAtElementCoarse(4,iel)

      ! Get the element numbers of the fine-grid elements in the coarse grid;
      ! the numbers are defined by the two-level ordering.
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(2,IELH3)

      ! Get the DOF`s on the fine grid
      I1=IedgesAtElementFine(1,IELH1)
      I2=IedgesAtElementFine(4,IELH2)
      I3=IedgesAtElementFine(1,IELH2)
      I4=IedgesAtElementFine(4,IELH3)
      I5=IedgesAtElementFine(1,IELH3)
      I6=IedgesAtElementFine(4,IELH4)
      I7=IedgesAtElementFine(1,IELH4)
      I8=IedgesAtElementFine(4,IELH1)
      I9=IedgesAtElementFine(2,IELH1)
      I10=IedgesAtElementFine(2,IELH2)
      I11=IedgesAtElementFine(2,IELH3)
      I12=IedgesAtElementFine(2,IELH4)

      ! Get the values of the DOF`s on the fine grid
      DUH1= DuFine(I1)
      DUH2= DuFine(I2)
      DUH3= DuFine(I3)
      DUH4= DuFine(I4)
      DUH5= DuFine(I5)
      DUH6= DuFine(I6)
      DUH7= DuFine(I7)
      DUH8= DuFine(I8)
      DUH9= DuFine(I9)
      DUH10=DuFine(I10)
      DUH11=DuFine(I11)
      DUH12=DuFine(I12)

      ! Now interpolate the fine-grid values to the coarse grid nodes.
      !
      ! Calculate the value of the edge IM1

      if (IneighboursAtElementCoarse(1,iel).ne.0) then
        ! inner edge
        if (IneighboursAtElementCoarse(1,iel).gt.iel) then
          DuCoarse(IM1)= A1*(DUH1+DUH2) +A2*DUH9 +A3*(DUH8+DUH3+DUH10+DUH12)
        else
          DuCoarse(IM1)= DuCoarse(IM1) &
                       + A1*(DUH1+DUH2) +A2*DUH9 +A3*(DUH8+DUH3+DUH10+DUH12)
        endif
      else
        ! boundary edge
        DuCoarse(IM1)= R1*(DUH1+DUH2) +R2*DUH9 +R3*(DUH8+DUH3+DUH10+DUH12)
      endif

      ! Calculate the value of the edge IM2

      if (IneighboursAtElementCoarse(2,iel).ne.0) then
        ! inner edge
        if (IneighboursAtElementCoarse(2,iel).gt.iel) then
           DuCoarse(IM2)= A1*(DUH3+DUH4) +A2*DUH10 +A3*(DUH2+DUH5+DUH9 +DUH11)
        else
          DuCoarse(IM2)= DuCoarse(IM2) &
                       + A1*(DUH3+DUH4) +A2*DUH10 +A3*(DUH2+DUH5+DUH9 +DUH11)
        endif
      else
        ! boundary edge
        DuCoarse(IM2)= R1*(DUH3+DUH4) +R2*DUH10 +R3*(DUH2+DUH5+DUH9 +DUH11)
      endif

      ! Calculate the value of the edge IM3

      if (IneighboursAtElementCoarse(3,iel).ne.0) then
        ! inner edge
        if (IneighboursAtElementCoarse(3,iel).gt.iel) then
           DuCoarse(IM3)= A1*(DUH5+DUH6) +A2*DUH11 +A3*(DUH4+DUH7+DUH10+DUH12)
        else
           DuCoarse(IM3)= DuCoarse(IM3) &
                        + A1*(DUH5+DUH6) +A2*DUH11 +A3*(DUH4+DUH7+DUH10+DUH12)
        endif
      else
        ! boundary edge
        DuCoarse(IM3)= R1*(DUH5+DUH6) +R2*DUH11 +R3*(DUH4+DUH7+DUH10+DUH12)
      endif

      ! Calculate the value of the edge IM4

      if (IneighboursAtElementCoarse(4,iel).ne.0) then
        ! inner edge
        if (IneighboursAtElementCoarse(4,iel).gt.iel) then
          DuCoarse(IM4)= A1*(DUH7+DUH8) +A2*DUH12 +A3*(DUH6+DUH1+DUH9 +DUH11)
        else
          DuCoarse(IM4)= DuCoarse(IM4) &
                       + A1*(DUH7+DUH8) +A2*DUH12 +A3*(DUH6+DUH1+DUH9 +DUH11)
        endif
      else
        ! boundary edge
        DuCoarse(IM4)= R1*(DUH7+DUH8) +R2*DUH12 +R3*(DUH6+DUH1+DUH9 +DUH11)
      endif

    end do

  end subroutine

  ! ***************************************************************************
  ! Support for Q1~ element, DOF`s = edge midpoints
  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_prolUniformEx31_double (DuCoarse,DuFine, &
               IedgesAtElementCoarse,IedgesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               NELcoarse)

!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! E031/EM31, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Coarse grid vector
  real(DP), dimension(:), intent(in) :: DuCoarse

  ! IedgesAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementCoarse

  ! IedgesAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Fine grid vector
  real(DP), dimension(:), intent(out) :: DuFine
!</output>

!</subroutine>

  ! local variables
  real(DP) :: DUH1,DUH2,DUH3,DUH4
  integer :: IM1,IM2, IM3,IM4, IA, IB, IC
  integer :: iel, IELH1, IELH2, IELH3, IELH4

  ! Weights for the restriction; all coefficients are halfed, so dividing
  ! by 2 is not necessary in the calculation routines.
  real(DP), parameter :: A1=0.46875_DP, A2=-0.09375_DP, A3=-0.03125_DP, A4=0.15625_DP
  real(DP), parameter :: A5=0.5625_DP, A6=0.1875_DP, A7=0.0625_DP, A8=0.1875_DP

    ! Clear the output vector
    call lalg_clearVector(DuFine)

    ! Loop over the coarse grid elements
    do iel=1,NELcoarse

      ! Get the DOF`s of the coarse grid element
      IM1 = IedgesAtElementCoarse(1,iel)
      IM2 = IedgesAtElementCoarse(2,iel)
      IM3 = IedgesAtElementCoarse(3,iel)
      IM4 = IedgesAtElementCoarse(4,iel)

      ! Get the values of the corresponding DOF`s
      DUH1 = DuCoarse(IM1)
      DUH2 = DuCoarse(IM2)
      DUH3 = DuCoarse(IM3)
      DUH4 = DuCoarse(IM4)

      ! Get the element numbers of the fine-grid elements in the coarse grid;
      ! the numbers are defined by the two-level ordering.
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(2,IELH3)

      ! Distribute the value at the edge IM1 to the
      ! corresponding fine inner nodes

      if (IneighboursAtElementCoarse(1,iel).ne.0) then
        ! There is a neighbour at the edge
        IA=IedgesAtElementFine(1,IELH1)
        IB=IedgesAtElementFine(4,IELH2)
        DuFine(IA)=DuFine(IA)+   A1*DUH1+A2*DUH2+A3*DUH3+A4*DUH4
        DuFine(IB)=DuFine(IB)+   A1*DUH1+A4*DUH2+A3*DUH3+A2*DUH4
      else
        ! No neighbour; boundary element
        IA=IedgesAtElementFine(1,IELH1)
        IB=IedgesAtElementFine(4,IELH2)
        DuFine(IA)=DuFine(IA)+2.0_DP*(A1*DUH1+A2*DUH2+A3*DUH3+A4*DUH4)
        DuFine(IB)=DuFine(IB)+2.0_DP*(A1*DUH1+A4*DUH2+A3*DUH3+A2*DUH4)
      endif
      IC=IedgesAtElementFine(2,IELH1)
      DuFine(IC)=A5*DUH1+A6*(DUH2+DUH4)+A7*DUH3

      ! Distribute the value at the edge IM2 to the
      ! corresponding fine inner nodes

      if (IneighboursAtElementCoarse(2,iel).ne.0) then
        ! There is a neighbour at the edge
       IA=IedgesAtElementFine(1,IELH2)
       IB=IedgesAtElementFine(4,IELH3)
       DuFine(IA)=DuFine(IA)+   A1*DUH2+A2*DUH3+A3*DUH4+A4*DUH1
       DuFine(IB)=DuFine(IB)+   A1*DUH2+A4*DUH3+A3*DUH4+A2*DUH1
      else
        ! No neighbour; boundary element
       IA=IedgesAtElementFine(1,IELH2)
       IB=IedgesAtElementFine(4,IELH3)
       DuFine(IA)=DuFine(IA)+2.0_DP*(A1*DUH2+A2*DUH3+A3*DUH4+A4*DUH1)
       DuFine(IB)=DuFine(IB)+2.0_DP*(A1*DUH2+A4*DUH3+A3*DUH4+A2*DUH1)
      endif
      IC=IedgesAtElementFine(2,IELH2)
      DuFine(IC)=A5*DUH2+A6*(DUH3+DUH1)+A7*DUH4

      ! Distribute the value at the edge IM3 to the
      ! corresponding fine inner nodes

      if (IneighboursAtElementCoarse(3,iel).ne.0) then
        ! There is a neighbour at the edge
       IA=IedgesAtElementFine(1,IELH3)
       IB=IedgesAtElementFine(4,IELH4)
       DuFine(IA)=DuFine(IA)+   A1*DUH3+A2*DUH4+A3*DUH1+A4*DUH2
       DuFine(IB)=DuFine(IB)+   A1*DUH3+A4*DUH4+A3*DUH1+A2*DUH2
      else
        ! No neighbour; boundary element
       IA=IedgesAtElementFine(1,IELH3)
       IB=IedgesAtElementFine(4,IELH4)
       DuFine(IA)=DuFine(IA)+2.0_DP*(A1*DUH3+A2*DUH4+A3*DUH1+A4*DUH2)
       DuFine(IB)=DuFine(IB)+2.0_DP*(A1*DUH3+A4*DUH4+A3*DUH1+A2*DUH2)
      endif
      IC=IedgesAtElementFine(2,IELH3)
      DuFine(IC)=A5*DUH3+A6*(DUH4+DUH2)+A7*DUH1

      ! Distribute the value at the edge IM4 to the
      ! corresponding fine inner nodes

      if (IneighboursAtElementCoarse(4,iel).ne.0) then
        ! There is a neighbour at the edge
        IA=IedgesAtElementFine(1,IELH4)
        IB=IedgesAtElementFine(4,IELH1)
        DuFine(IA)=DuFine(IA)+   A1*DUH4+A2*DUH1+A3*DUH2+A4*DUH3
        DuFine(IB)=DuFine(IB)+   A1*DUH4+A4*DUH1+A3*DUH2+A2*DUH3
      else
        ! No neighbour; boundary element
        IA=IedgesAtElementFine(1,IELH4)
        IB=IedgesAtElementFine(4,IELH1)
        DuFine(IA)=DuFine(IA)+2.0_DP*(A1*DUH4+A2*DUH1+A3*DUH2+A4*DUH3)
        DuFine(IB)=DuFine(IB)+2.0_DP*(A1*DUH4+A4*DUH1+A3*DUH2+A2*DUH3)
      end if
      IC=IedgesAtElementFine(2,IELH4)
      DuFine(IC)=A5*DUH4+A6*(DUH1+DUH3)+A7*DUH2

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_prolUniformEx31ext_double (DuCoarse,DuFine, &
               DvertexCoordsCoarse,IverticesAtElementCoarse,DelementAreaCoarse,&
               IedgesAtElementCoarse,IedgesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               NELcoarse, iweightingType, daspectRatioBound, iarIndicator)

!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! E031/EM31, uniform triangulation, double precision vector.
  !
  ! Extended version. Switch to constant prolongation if aspect ratio
  ! of an element is too large.
!</description>

!<input>
  ! Coarse grid vector
  real(DP), dimension(:), intent(in) :: DuCoarse

  ! DvertexCoords array on the coarse grid
  real(DP), dimension(:,:), intent(in)                :: DvertexCoordsCoarse

  ! IverticesAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IverticesAtElementCoarse

  ! DelementArea array on the coarse grid
  real(DP), dimension(:), intent(in)                  :: DelementAreaCoarse

  ! IedgesAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementCoarse

  ! IedgesAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse

  ! Type of the averaging on the element edges
  ! <=0: standard averaging of both contributions by 1/2,
  !  =1: weighted averaging of the interpolated function values:
  !      The area of the current coarse grid element determines
  !      the weight. (L2-projection, standard),
  !  =2: weighted averaging of the interpolated function values:
  !      The area of the neightbour element of the coarse grid
  !      the weight.
  integer, intent(in)  :: iweightingType

  ! Upper bound aspect ratio; for all elements with higher AR
  ! the prolongation is switched to constant prolongation
  real(DP), intent(in) :: daspectRatioBound

  ! Aspect-ratio indicator.
  ! Controls switching to constant prolongation.
  ! <=1: switch depending on aspect ratio of current element,
  !  =2: switch depending on aspect ratio of current element and
  !      neighbour element
  integer, intent(in)  :: iarIndicator
!</input>

!<output>
  ! Fine grid vector
  real(DP), dimension(:), intent(out) :: DuFine
!</output>

!</subroutine>

  ! local variables
  real(DP) :: DUH1,DUH2,DUH3,DUH4
  real(DP), dimension(0:TRIA_MAXNME2D) :: daspectRatio,darea
  real(DP), dimension(TRIA_MAXNME2D) :: dweight
  integer :: IM1,IM2, IM3,IM4, IA, IB, IC
  integer :: iel, IELH1, IELH2, IELH3, IELH4
  integer, dimension(0:TRIA_MAXNME2D) :: IELA
  integer :: i
  integer, dimension(TRIA_MAXNME2D) :: idoConstant
  real(DP), dimension(NDIM2D,TRIA_MAXNVE2D) :: dcoords

  ! Weights for the prolongation.
  ! PRWEIG (.,1) gives the constants for the standard prolongation,
  ! PRWEIG (.,2) gives the constants for the constant prolongation.
  real(DP), dimension(8,2), parameter :: prweight = &
      reshape((/0.9375_DP, -0.1875_DP, -0.0625_DP, 0.3125_DP, &
                0.5625_DP, 0.1875_DP, 0.0625_DP, 0.1875_DP, &
                1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, &
                1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP/),(/8,2/))

    ! Clear the output vector
    call lalg_clearVector(DuFine)

    ! Loop over the coarse grid elements
    do iel=1,NELcoarse

      ! Get the numbers of the elements that are neighbours to our current coarse
      ! grid element:
      !             +--------+
      !             |        |
      !             | IELA3  |
      !             |        |
      !    +--------4--------3--------+
      !    |        |        |        |
      !    | IELA4  |  IEL   | IELA2  |
      !    |        |        |        |
      !    +--------1--------2--------+
      !             |        |
      !             | IELA1  |
      !             |        |
      !             +--------+
      IELA(0) = iel
      IELA(1) = IneighboursAtElementCoarse(1,iel)
      IELA(2) = IneighboursAtElementCoarse(2,iel)
      IELA(3) = IneighboursAtElementCoarse(3,iel)
      IELA(4) = IneighboursAtElementCoarse(4,iel)

      ! For these five elements, determine the aspect ratio and their area.
      !
      ! At first the element in the center, which always exists.
      !
      ! Get the aspect ratio of the current coarse grid element;
      ! if necessary, calculate the reciprocal.
      dcoords = DvertexCoordsCoarse(:,IverticesAtElementCoarse(:,IELA(0)))
      daspectRatio(0) = gaux_getAspectRatio_quad2D (dcoords)
      if (daspectRatio(0) .lt. 1.0_DP) daspectRatio(0) = 1.0_DP/daspectRatio(0)

      ! and the area of that element.
      darea(0) = DelementAreaCoarse(iel)

      ! Then the remaining neighbours.
      do i=1,TRIA_MAXNME2D
        if (IELA(i) .ne. 0) then
          ! Get the aspect ratio of the current coarse grid element;
          ! if necessary, calculate the reciprocal.
          dcoords = DvertexCoordsCoarse(:,IverticesAtElementCoarse(:,IELA(i)))
          daspectRatio(i) = gaux_getAspectRatio_quad2D (dcoords)
          if (daspectRatio(i) .lt. 1.0_DP) daspectRatio(i) = 1.0_DP/daspectRatio(i)

          ! and the area of that element.
          darea(i) = DelementAreaCoarse(IELA(i))
        else
          daspectRatio(i) = 0.0_DP
          darea(i) = 0.0_DP
        end if
      end do

      ! Calculate weighting factors for the interpolation.
      ! The iweightingType parameter describes:
      ! <= 0: simple interpolation, weight both contribtutions by 1/2
      !  = 1: take the weighted mean of the interpolated function values
      !       by weighting with the area of the current coarse grid element
      ! >= 2: take the weighted mean of the interpolated function values
      !       by weighting with the area of the neighboured coarse grid element

      select case (iweightingType)
      case (:0)
        dweight = 0.5_DP
      case (1)
        dweight = darea(0) / (darea(0)+darea(1:TRIA_MAXNME2D))
      case (2:)
        dweight = darea(1:TRIA_MAXNME2D) / (darea(0)+darea(1:TRIA_MAXNME2D))
      end select

      ! Where there is no neighbour, set the weighting factor to 1.0
!      DO i=1,TRIA_MAXNMT
!        IF (IELA(i) .EQ. 0) dweight(i) = 1.0_DP
!      END DO
      where(IELA(1:TRIA_MAXNME2D) .eq. 0) dweight(1:TRIA_MAXNME2D) = 1.0_DP

      ! Now determine on which edge to switch to constant prolongation
      ! By default, we do not use constant prolongation
      idoConstant = 1

      ! ... but if the caller wants us to switch in a special situation...
      if ((iarIndicator .ge. 1) .and. (daspectRatioBound .ge. 0.0_DP)) then

        ! ... switch to constant of our element is too large...
        if (daspectRatio(0) .gt. daspectRatioBound) idoConstant = 2

        ! and if iarIndicator>2, also check the neighbour element
        if (iarIndicator .ge. 2) then
!          DO i=1,TRIA_MAXNME2D
!            IF (daspectRatio(i) .GT. daspectRatioBound) idoConstant(i) = 2
!          END DO
          where (daspectRatio(1:4) .gt. daspectRatioBound) idoConstant = 2
        end if

      end if

      ! Get the DOF`s of the coarse grid element
      IM1 = IedgesAtElementCoarse(1,iel)
      IM2 = IedgesAtElementCoarse(2,iel)
      IM3 = IedgesAtElementCoarse(3,iel)
      IM4 = IedgesAtElementCoarse(4,iel)

      ! Get the values of the corresponding DOF`s
      DUH1 = DuCoarse(IM1)
      DUH2 = DuCoarse(IM2)
      DUH3 = DuCoarse(IM3)
      DUH4 = DuCoarse(IM4)

      ! Get the element numbers of the fine-grid elements in the coarse grid;
      ! the numbers are defined by the two-level ordering.
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(2,IELH3)

      ! Now let us start with the actual prolongation
      ! ---------------------------------------------

      ! Get the DOF`s on the fine grid
      IA=IedgesAtElementFine(1,IELH1)
      IB=IedgesAtElementFine(4,IELH2)
      IC=IedgesAtElementFine(2,IELH1)

      ! Now we have the following situation:

      !   4               IM3                3
      !     ===============X================
      !     |              |               |
      !     |              |               |
      !     |    IELH4     |     IELH3     |
      !     |              |               |
      !     |                              |
      ! IM4 X----------- IEL1 -------------X IM2
      !     |                              |
      !     |              |               |
      !     |    IELH1     o IC  IELH2     |
      !     |              |               |
      !     |              |               |
      !   1 =======o=======X=======o======== 2
      !     |     IA      IM1      IB      |
      !     |                              |
      !     |                              |
      !     |                              |
      !     |                              |
      !     |            IELA1             |
      !     |                              |
      !     |                              |
      !     |                              |
      !     |                              |
      !     |                              |
      !     ================================

      ! Distribute the value at the edge IM1 to the
      ! corresponding fine inner nodes

      DuFine(IA) = DuFine(IA) &
                  + dweight(1)*(prweight(1,idoConstant(1))*DUH1 &
                              +prweight(2,idoConstant(1))*DUH2 &
                              +prweight(3,idoConstant(1))*DUH3 &
                              +prweight(4,idoConstant(1))*DUH4)
      DuFine(IB) = DuFine(IB) &
                  + dweight(1)*(prweight(1,idoConstant(1))*DUH1 &
                              +prweight(4,idoConstant(1))*DUH2 &
                              +prweight(3,idoConstant(1))*DUH3 &
                              +prweight(2,idoConstant(1))*DUH4)
      DuFine(IC) = prweight(5,idoConstant(1))*DUH1 &
                 + prweight(6,idoConstant(1))*(DUH2+DUH4) &
                 + prweight(7,idoConstant(1))*DUH3

      ! Distribute the value at the edge IM2 to the
      ! corresponding fine inner nodes
      IA=IedgesAtElementFine(1,IELH2)
      IB=IedgesAtElementFine(4,IELH3)
      IC=IedgesAtElementFine(2,IELH2)

      DuFine(IA) = DuFine(IA) &
                 + dweight(2)*(prweight(1,idoConstant(2))*DUH2 &
                              +prweight(2,idoConstant(2))*DUH3 &
                              +prweight(3,idoConstant(2))*DUH4 &
                              +prweight(4,idoConstant(2))*DUH1)
      DuFine(IB) = DuFine(IB) &
                 + dweight(2)*(prweight(1,idoConstant(2))*DUH2 &
                              +prweight(4,idoConstant(2))*DUH3 &
                              +prweight(3,idoConstant(2))*DUH4 &
                              +prweight(2,idoConstant(2))*DUH1)
      DuFine(IC) = prweight(5,idoConstant(2))*DUH2 &
                 + prweight(6,idoConstant(2))*(DUH3+DUH1) &
                 + prweight(7,idoConstant(2))*DUH4

      ! Distribute the value at the edge IM3 to the
      ! corresponding fine inner nodes
      IA=IedgesAtElementFine(1,IELH3)
      IB=IedgesAtElementFine(4,IELH4)
      IC=IedgesAtElementFine(2,IELH3)

      DuFine(IA) = DuFine(IA) &
                 + dweight(3)*(prweight(1,idoConstant(3))*DUH3 &
                              +prweight(2,idoConstant(3))*DUH4 &
                              +prweight(3,idoConstant(3))*DUH1 &
                              +prweight(4,idoConstant(3))*DUH2)
      DuFine(IB) = DuFine(IB) &
                 + dweight(3)*(prweight(1,idoConstant(3))*DUH3 &
                              +prweight(4,idoConstant(3))*DUH4 &
                              +prweight(3,idoConstant(3))*DUH1 &
                              +prweight(2,idoConstant(3))*DUH2)
      DuFine(IC) = prweight(5,idoConstant(3))*DUH3 &
                 + prweight(6,idoConstant(3))*(DUH4+DUH2) &
                 + prweight(7,idoConstant(3))*DUH1

      ! Distribute the value at the edge IM4 to the
      ! corresponding fine inner nodes
      IA=IedgesAtElementFine(1,IELH4)
      IB=IedgesAtElementFine(4,IELH1)
      IC=IedgesAtElementFine(2,IELH4)

      DuFine(IA) = DuFine(IA) &
                 + dweight(4)*(prweight(1,idoConstant(4))*DUH4 &
                              +prweight(2,idoConstant(4))*DUH1 &
                              +prweight(3,idoConstant(4))*DUH2 &
                              +prweight(4,idoConstant(4))*DUH3)
      DuFine(IB) = DuFine(IB) &
                 + dweight(4)*(prweight(1,idoConstant(4))*DUH4 &
                              +prweight(4,idoConstant(4))*DUH1 &
                              +prweight(3,idoConstant(4))*DUH2 &
                              +prweight(2,idoConstant(4))*DUH3)
      DuFine(IC) = prweight(5,idoConstant(4))*DUH4 &
                 + prweight(6,idoConstant(4))*(DUH1+DUH3) &
                 + prweight(7,idoConstant(4))*DUH2

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_restUniformEx31_double (DuCoarse,DuFine, &
               IedgesAtElementCoarse,IedgesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,NELcoarse)

!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! E031/EM31 element, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IedgesAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementCoarse

  ! IedgesAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

    ! local variables
    real(DP) :: DUH1,DUH2,DUH3,DUH4,DUH5,DUH6,DUH7,DUH8,DUH9,DUH10,DUH11,DUH12
    integer :: IM1,IM2,IM3,IM4
    integer :: I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12
    integer :: iel,IELH1,IELH2,IELH3,IELH4

    ! Weights for the restriction
    real(DP), parameter :: A1=0.9375_DP, A2=-0.09375_DP, A3=-0.03125_DP, A4=0.15625_DP
    real(DP), parameter :: A5=0.5625_DP, A6=0.1875_DP, A7=0.0625_DP, A8=0.1875_DP

    ! Loop over the coarse grid elements
    do iel=1,NELcoarse

      ! Get the DOF`s of the coarse grid element
      IM1 = IedgesAtElementCoarse(1,iel)
      IM2 = IedgesAtElementCoarse(2,iel)
      IM3 = IedgesAtElementCoarse(3,iel)
      IM4 = IedgesAtElementCoarse(4,iel)

      ! Get the element numbers of the fine-grid elements in the coarse grid;
      ! the numbers are defined by the two-level ordering.
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(2,IELH3)

      ! Get the DOF`s on the fine grid
      I1=IedgesAtElementFine(1,IELH1)
      I2=IedgesAtElementFine(4,IELH2)
      I3=IedgesAtElementFine(1,IELH2)
      I4=IedgesAtElementFine(4,IELH3)
      I5=IedgesAtElementFine(1,IELH3)
      I6=IedgesAtElementFine(4,IELH4)
      I7=IedgesAtElementFine(1,IELH4)
      I8=IedgesAtElementFine(4,IELH1)
      I9=IedgesAtElementFine(2,IELH1)
      I10=IedgesAtElementFine(2,IELH2)
      I11=IedgesAtElementFine(2,IELH3)
      I12=IedgesAtElementFine(2,IELH4)

      ! Get the values of the DOF`s on the fine grid
      DUH1= DuFine(I1)
      DUH2= DuFine(I2)
      DUH3= DuFine(I3)
      DUH4= DuFine(I4)
      DUH5= DuFine(I5)
      DUH6= DuFine(I6)
      DUH7= DuFine(I7)
      DUH8= DuFine(I8)
      DUH9= DuFine(I9)
      DUH10=DuFine(I10)
      DUH11=DuFine(I11)
      DUH12=DuFine(I12)

      ! Now collect the information in the same way as it was
      ! "distributed" by the prolongation routine.
      ! This realises the adjoint operator of the prolongation.
      !
      ! Calculate the value of the edge IM1

      if (IneighboursAtElementCoarse(1,iel).ne.0) then
        ! inner edge
        if (IneighboursAtElementCoarse(1,iel).gt.iel) then
          DuCoarse(IM1)= A1*(DUH1+DUH2)+A2*(DUH4+DUH7) &
                       + A3*(DUH5+DUH6)+A4*(DUH3+DUH8) &
                       + A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
        else
          DuCoarse(IM1)= DuCoarse(IM1)+A2*(DUH4+DUH7) &
                       + A3*(DUH5+DUH6)+A4*(DUH3+DUH8) &
                       + A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
        endif
      else
        ! boundary edge
        DuCoarse(IM1)=     A1*(DUH1+DUH2)+2.0_DP*A2*(DUH4+DUH7) &
                +2.0_DP*A3*(DUH5+DUH6)+2.0_DP*A4*(DUH3+DUH8) &
                +       A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
      endif

      ! Calculate the value of the edge IM2

      if (IneighboursAtElementCoarse(2,iel).ne.0) then
        ! inner edge
        if (IneighboursAtElementCoarse(2,iel).gt.iel) then
           DuCoarse(IM2)= A1*(DUH3+DUH4)+A2*(DUH6+DUH1) &
                         +A3*(DUH7+DUH8)+A4*(DUH5+DUH2) &
                         +A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
        else
          DuCoarse(IM2)=DuCoarse(IM2)+A2*(DUH6+DUH1) &
                              +A3*(DUH7+DUH8)+A4*(DUH5+DUH2) &
                              +A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
        endif
      else
        ! boundary edge
        DuCoarse(IM2)= A1*(DUH3+DUH4)+2.0_DP*A2*(DUH6+DUH1) &
                +2.0_DP*A3*(DUH7+DUH8)+2.0_DP*A4*(DUH5+DUH2) &
                +       A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
      endif

      ! Calculate the value of the edge IM3

      if (IneighboursAtElementCoarse(3,iel).ne.0) then
        ! inner edge
        if (IneighboursAtElementCoarse(3,iel).gt.iel) then
           DuCoarse(IM3)= A1*(DUH5+DUH6)+A2*(DUH8+DUH3) &
                        + A3*(DUH1+DUH2)+A4*(DUH7+DUH4) &
                        + A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
        else
           DuCoarse(IM3)= DuCoarse(IM3)+A2*(DUH8+DUH3) &
                        + A3*(DUH1+DUH2)+A4*(DUH7+DUH4) &
                        + A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
        endif
      else
        ! boundary edge
        DuCoarse(IM3)= A1*(DUH5+DUH6)+2.0_DP*A2*(DUH8+DUH3) &
                +2.0_DP*A3*(DUH1+DUH2)+2.0_DP*A4*(DUH7+DUH4) &
                +       A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
      endif

      ! Calculate the value of the edge IM4

      if (IneighboursAtElementCoarse(4,iel).ne.0) then
        ! inner edge
        if (IneighboursAtElementCoarse(4,iel).gt.iel) then
          DuCoarse(IM4)= A1*(DUH7+DUH8)+A2*(DUH2+DUH5) &
                       + A3*(DUH3+DUH4)+A4*(DUH1+DUH6) &
                       + A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
        else
          DuCoarse(IM4)= DuCoarse(IM4)+A2*(DUH2+DUH5) &
                       + A3*(DUH3+DUH4)+A4*(DUH1+DUH6) &
                       + A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
        endif
      else
        ! boundary edge
        DuCoarse(IM4)= A1*(DUH7+DUH8)+2.0_DP*A2*(DUH2+DUH5) &
                +2.0_DP*A3*(DUH3+DUH4)+2.0_DP*A4*(DUH1+DUH6) &
                +       A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
      endif

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_restUniformEx31ext_double (DuCoarse,DuFine, &
               DvertexCoordsCoarse,IverticesAtElementCoarse,DelementAreaCoarse,&
               IedgesAtElementCoarse,IedgesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               NELcoarse, iweightingType, daspectRatioBound, iarIndicator)

!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! E031/EM31 element, uniform triangulation, double precision vector.
  !
  ! Extended version. Switch to constant restriction if aspect ratio
  ! of an element is too large.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! DvertexCoords array on the coarse grid
  real(DP), dimension(:,:), intent(in)                :: DvertexCoordsCoarse

  ! IverticesAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IverticesAtElementCoarse

  ! DelementArea array on the coarse grid
  real(DP), dimension(:), intent(in)                  :: DelementAreaCoarse

  ! IedgesAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementCoarse

  ! IedgesAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse

  ! Type of the averaging on the element edges
  ! <=0: standard averaging of both contributions by 1/2,
  !  =1: weighted averaging of the interpolated function values:
  !      The area of the current coarse grid element determines
  !      the weight. (L2-projection, standard),
  !  =2: weighted averaging of the interpolated function values:
  !      The area of the neightbour element of the coarse grid
  !      the weight.
  integer, intent(in)  :: iweightingType

  ! Upper bound aspect ratio; for all elements with higher AR
  ! the prolongation is switched to constant prolongation
  real(DP), intent(in) :: daspectRatioBound

  ! Aspect-ratio indicator.
  ! Controls switching to constant prolongation.
  ! <=1: switch depending on aspect ratio of current element,
  !  =2: switch depending on aspect ratio of current element and
  !      neighbour element
  integer, intent(in)  :: iarIndicator
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

    ! local variables
    real(DP) :: DUH1,DUH2,DUH3,DUH4,DUH5,DUH6,DUH7,DUH8,DUH9,DUH10,DUH11,DUH12
    integer :: IM1,IM2, IM3,IM4
    integer :: I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12
    integer :: iel, IELH1, IELH2, IELH3, IELH4
    integer, dimension(0:TRIA_MAXNME2D) :: IELA
    integer :: i
    integer, dimension(TRIA_MAXNME2D) :: idoConstant
    real(DP), dimension(0:TRIA_MAXNME2D) :: daspectRatio, darea
    real(DP), dimension(TRIA_MAXNME2D) :: dweight
    real(DP), dimension(NDIM2D,TRIA_MAXNVE2D) :: dcoords

    ! Weights for the restriction.
    ! PRWEIG (.,1) gives the constants for the standard restriction,
    ! PRWEIG (.,2) gives the constants for the constant restriction.
    real(DP), dimension(8,2), parameter :: prweight = &
        reshape((/0.9375_DP, -0.1875_DP, -0.0625_DP, 0.3125_DP, &
                  0.5625_DP, 0.1875_DP, 0.0625_DP, 0.1875_DP, &
                  1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, &
                  1.0_DP, 0.0_DP, 0.0_DP, 0.0_DP/),(/8,2/))

    ! Clear the output vector
    call lalg_clearVector(DuCoarse)

    ! Loop over the coarse grid elements
    do iel=1,NELcoarse

      ! Get the numbers of the elements that are neighbours to our current coarse
      ! grid element:
      !             +--------+
      !             |        |
      !             | IELA3  |
      !             |        |
      !    +--------4--------3--------+
      !    |        |        |        |
      !    | IELA4  |  IEL   | IELA2  |
      !    |        |        |        |
      !    +--------1--------2--------+
      !             |        |
      !             | IELA1  |
      !             |        |
      !             +--------+
      IELA(0) = iel
      IELA(1) = IneighboursAtElementCoarse(1,iel)
      IELA(2) = IneighboursAtElementCoarse(2,iel)
      IELA(3) = IneighboursAtElementCoarse(3,iel)
      IELA(4) = IneighboursAtElementCoarse(4,iel)

      ! For these five elements, determine the aspect ratio and their area.
      !
      ! At first the element in the center, which always exists.
      !
      ! Get the aspect ratio of the current coarse grid element;
      ! if necessary, calculate the reciprocal.
      dcoords = DvertexCoordsCoarse(:,IverticesAtElementCoarse(:,IELA(0)))
      daspectRatio(0) = gaux_getAspectRatio_quad2D (dcoords)
      if (daspectRatio(0) .lt. 1.0_DP) daspectRatio(0) = 1.0_DP/daspectRatio(0)

      ! and the area of that element.
      darea(0) = DelementAreaCoarse(iel)

      ! Then the remaining neighbours.
      do i=1,TRIA_MAXNME2D
        if (IELA(i) .ne. 0) then
          ! Get the aspect ratio of the current coarse grid element;
          ! if necessary, calculate the reciprocal.
          dcoords = DvertexCoordsCoarse(:,IverticesAtElementCoarse(:,IELA(i)))
          daspectRatio(i) = gaux_getAspectRatio_quad2D (dcoords)
          if (daspectRatio(i) .lt. 1.0_DP) daspectRatio(i) = 1.0_DP/daspectRatio(i)

          ! and the area of that element.
          darea(i) = DelementAreaCoarse(IELA(i))
        else
          daspectRatio(i) = 0.0_DP
          darea(i) = 0.0_DP
        end if
      end do

      ! Calculate weighting factors for the interpolation.
      ! The iweightingType parameter describes:
      ! <= 0: simple interpolation, weight both contribtutions by 1/2
      !  = 1: take the weighted mean of the interpolated function values
      !       by weighting with the area of the current coarse grid element
      ! >= 2: take the weighted mean of the interpolated function values
      !       by weighting with the area of the neighboured coarse grid element

      select case (iweightingType)
      case (:0)
        dweight = 0.5_DP
      case (1)
        dweight = darea(0) / (darea(0)+darea(1:TRIA_MAXNME2D))
      case (2:)
        dweight = darea(1:TRIA_MAXNME2D) / (darea(0)+darea(1:TRIA_MAXNME2D))
      end select

      ! Where there is no neighbour, set the weighting factor to 1.0
!      DO i=1,TRIA_MAXNMT
!        IF (IELA(i) .EQ. 0) dweight(i) = 1.0_DP
!      END DO
      where(IELA(1:TRIA_MAXNME2D) .eq. 0) dweight(1:TRIA_MAXNME2D) = 1.0_DP

      ! Now determine on which edge to switch to constant prolongation
      ! By default, we do not use constant prolongation
      idoConstant = 1

      ! ... but if the caller wants us to switch in a special situation...
      if ((iarIndicator .ge. 1) .and. (daspectRatioBound .ge. 0.0_DP)) then

        ! ... switch to constant of our element is too large...
        if (daspectRatio(0) .gt. daspectRatioBound) idoConstant = 2

        ! and if iarIndicator>2, also check the neighbour element
        if (iarIndicator .ge. 2) then
!          DO i=1,TRIA_MAXNME2D
!            IF (daspectRatio(i) .GT. daspectRatioBound) idoConstant(i) = 2
!          END DO
          where (daspectRatio(1:4) .gt. daspectRatioBound) idoConstant = 2
        end if

      end if

      ! Now let us strt with the actual restriction
      ! ------------------------------------------

      ! Get the DOF`s of the coarse grid element
      IM1 = IedgesAtElementCoarse(1,iel)
      IM2 = IedgesAtElementCoarse(2,iel)
      IM3 = IedgesAtElementCoarse(3,iel)
      IM4 = IedgesAtElementCoarse(4,iel)

      ! Get the element numbers of the fine-grid elements in the coarse grid;
      ! the numbers are defined by the two-level ordering.
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(2,IELH3)

      ! Get the DOF`s on the fine grid
      I1=IedgesAtElementFine(1,IELH1)
      I2=IedgesAtElementFine(4,IELH2)
      I3=IedgesAtElementFine(1,IELH2)
      I4=IedgesAtElementFine(4,IELH3)
      I5=IedgesAtElementFine(1,IELH3)
      I6=IedgesAtElementFine(4,IELH4)
      I7=IedgesAtElementFine(1,IELH4)
      I8=IedgesAtElementFine(4,IELH1)
      I9=IedgesAtElementFine(2,IELH1)
      I10=IedgesAtElementFine(2,IELH2)
      I11=IedgesAtElementFine(2,IELH3)
      I12=IedgesAtElementFine(2,IELH4)

      ! Get the values of the DOF`s on the fine grid
      DUH1= DuFine(I1)
      DUH2= DuFine(I2)
      DUH3= DuFine(I3)
      DUH4= DuFine(I4)
      DUH5= DuFine(I5)
      DUH6= DuFine(I6)
      DUH7= DuFine(I7)
      DUH8= DuFine(I8)
      DUH9= DuFine(I9)
      DUH10=DuFine(I10)
      DUH11=DuFine(I11)
      DUH12=DuFine(I12)

      ! Now we have the following situation:
      !   4       I6      IM3      I5      3
      !     =======o=======X=======o========
      !     |              |               |
      !     |              |               |
      !  I7 o    IELH4     o I11 IELH3     o I4
      !     |              |               |
      !     |                              |
      ! IM4 X------o---- IEL1 -----o-------X IM2
      !     |    I12               I10     |
      !     |              |               |
      !  I8 o    IELH1     o I9  IELH2     o I3
      !     |              |               |
      !     |              |               |
      !   1 =======o=======X=======o======== 2
      !     |     I1      IM1      I2      |
      !     |                              |
      !     |                              |
      !     |                              |
      !     |                              |
      !     |            IELA1             |
      !     |                              |
      !     |                              |
      !     |                              |
      !     |                              |
      !     |                              |
      !     ================================


      ! Collect the information in the same way as it was
      ! "distributed" by the prolongation routine.
      ! This realises the adjoint operator of the prolongation.
      !
      ! Calculate the value of the edge IM1

      DuCoarse(IM1)= DuCoarse(IM1) &
                    +dweight(1)*(prweight(1,idoConstant(1))*(DUH1+DUH2) &
                                +prweight(2,idoConstant(1))*(DUH4+DUH7) &
                                +prweight(3,idoConstant(1))*(DUH5+DUH6) &
                                +prweight(4,idoConstant(1))*(DUH3+DUH8)) &
                    +prweight(5,idoConstant(1))*DUH9 &
                    +prweight(6,idoConstant(1))*(DUH10+DUH12) &
                    +prweight(7,idoConstant(1))*DUH11

      ! Calculate the value of the edge IM2

      DuCoarse(IM2)= DuCoarse(IM2) &
                    +dweight(2)*(prweight(1,idoConstant(2))*(DUH3+DUH4) &
                                +prweight(2,idoConstant(2))*(DUH6+DUH1) &
                                +prweight(3,idoConstant(2))*(DUH7+DUH8) &
                                +prweight(4,idoConstant(2))*(DUH5+DUH2)) &
                    +prweight(5,idoConstant(2))*DUH10 &
                    +prweight(6,idoConstant(2))*(DUH11+DUH9) &
                    +prweight(7,idoConstant(2))*DUH12

      ! Calculate the value of the edge IM3

      DuCoarse(IM3)= DuCoarse(IM3) &
                    +dweight(3)*(prweight(1,idoConstant(3))*(DUH5+DUH6) &
                                +prweight(2,idoConstant(3))*(DUH8+DUH3) &
                                +prweight(3,idoConstant(3))*(DUH1+DUH2) &
                                +prweight(4,idoConstant(3))*(DUH7+DUH4)) &
                    +prweight(5,idoConstant(3))*DUH11 &
                    +prweight(6,idoConstant(3))*(DUH12+DUH10) &
                    +prweight(7,idoConstant(3))*DUH9

      ! Calculate the value of the edge IM4

      DuCoarse(IM4)= DuCoarse(IM4) &
                    +dweight(4)*(prweight(1,idoConstant(4))*(DUH7+DUH8) &
                                +prweight(2,idoConstant(4))*(DUH2+DUH5) &
                                +prweight(3,idoConstant(4))*(DUH3+DUH4) &
                                +prweight(4,idoConstant(4))*(DUH1+DUH6)) &
                    +prweight(5,idoConstant(4))*DUH12 &
                    +prweight(6,idoConstant(4))*(DUH9+DUH11) &
                    +prweight(7,idoConstant(4))*DUH10

    end do

  end subroutine

  ! ***************************************************************************
  ! Support for 3D Q0 element
  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_prolUniformQ0_3D_double (DuCoarse,DuFine,NELcoarse)

!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! $Q_0$, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Coarse grid vector
  real(DP), dimension(:), intent(in) :: DuCoarse

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Fine grid vector
  real(DP), dimension(:), intent(out) :: DuFine
!</output>

!</subroutine>

  ! local variables
  integer :: iel
  integer, dimension(8) :: ielf
  real(DP) :: duh

    ! Loop over the elements
    do iel = 1, NELCoarse

      ! Get the element numbers of the fine-grid elements in the coarse grid;
      ! the numbers are defined by the two-level ordering.
      ielf(1) = iel
      ielf(2) = NELcoarse + 7*(iel-1) + 1
      ielf(3) = ielf(2) + 1
      ielf(4) = ielf(3) + 1
      ielf(5) = ielf(4) + 1
      ielf(6) = ielf(5) + 1
      ielf(7) = ielf(6) + 1
      ielf(8) = ielf(7) + 1

      ! Put the value on the coarse grid into all four child
      ! elements
      duh = DuCoarse(iel)
      DuFine(ielf(1:8)) = duh
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_restUniformQ0_3D_double (DuCoarse,DuFine,NELcoarse)

!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! $Q_0$, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  integer :: iel
  integer, dimension(8) :: ielf

    ! The information that was "distributed" in the prolongation has to
    ! be "collected".

    ! Loop over the elements to collect the missing additive contributions:
    do iel=1,NELcoarse

      ! Get the element numbers of the fine-grid elements in the coarse grid;
      ! the numbers are defined by the two-level ordering.
      ielf(1) = iel
      ielf(2) = NELcoarse + 7*(iel-1) + 1
      ielf(3) = ielf(2) + 1
      ielf(4) = ielf(3) + 1
      ielf(5) = ielf(4) + 1
      ielf(6) = ielf(5) + 1
      ielf(7) = ielf(6) + 1
      ielf(8) = ielf(7) + 1

      ! Sum up the values in these nodes to get the
      ! value in the coarse grid element
      DuCoarse(iel)= DuFine(ielf(1))+DuFine(ielf(2))+&
                     DuFine(ielf(3))+DuFine(ielf(4))+DuFine(ielf(5))+&
                     DuFine(ielf(6))+DuFine(ielf(7))+DuFine(ielf(8))

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_interpUniformQ0_3D_double (DuCoarse,DuFine,NELcoarse)

!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! $Q_0$, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  integer :: iel
  integer, dimension(8) :: ielf

    ! The information that was "distributed" in the prolongation has to
    ! be "collected".

    ! Loop over the elements to collect the missing additive contributions:
    do iel=1,NELcoarse

      ! Get the element numbers of the fine-grid elements in the coarse grid;
      ! the numbers are defined by the two-level ordering.
      ielf(1) = iel
      ielf(2) = NELcoarse + 7*(iel-1) + 1
      ielf(3) = ielf(2) + 1
      ielf(4) = ielf(3) + 1
      ielf(5) = ielf(4) + 1
      ielf(6) = ielf(5) + 1
      ielf(7) = ielf(6) + 1
      ielf(8) = ielf(7) + 1

      ! Sum up the values in these nodes to get the
      ! value in the coarse grid element
      DuCoarse(iel)= 0.125_DP * (DuFine(ielf(1))+DuFine(ielf(2))+&
                     DuFine(ielf(3))+DuFine(ielf(4))+DuFine(ielf(5))+&
                     DuFine(ielf(6))+DuFine(ielf(7))+DuFine(ielf(8)))

    end do

  end subroutine

  ! ***************************************************************************
  ! Support for 3D Q1 element
  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_prolUniformQ1_3D_double (DuCoarse, DuFine, &
             IverticesAtEdgeCoarse, IverticesAtFaceCoarse, &
             IverticesAtElementCoarse, NVTcoarse, NMTcoarse, NATcoarse,NELcoarse)

!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! <tex>$Q_1$</tex>, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Coarse grid vector
  real(DP), dimension(:), intent(in) :: DuCoarse

  ! IverticesAtEdge array on the coarse grid.
  integer, dimension(:,:), intent(in) :: IverticesAtEdgeCoarse

  ! IverticesAtFace array on the coarse grid.
  integer, dimension(:,:), intent(in) :: IverticesAtFaceCoarse

  ! IverticesAtElement array (KVERT) on the coarse grid
  integer, dimension(:,:), intent(in) :: IverticesAtElementCoarse

  ! Number of vertices in the coarse grid
  integer, intent(in) :: NVTcoarse

  ! Number of edges in the coarse grid
  integer, intent(in) :: NMTcoarse

  ! Number of faces in the coarse grid
  integer, intent(in) :: NATcoarse

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Fine grid vector
  real(DP), dimension(:), intent(out) :: DuFine
!</output>

!</subroutine>

  ! local variables
  integer :: iel
  integer :: iedge, iface

    ! Copy the first NVT entries - they belong to the coarse grid vertices
    ! that are fine grid vertices at the same time.
    call lalg_copyVector (DuCoarse,DuFine(1:NVTcoarse))

    ! Loop over the edges
    do iedge = 1, NMTcoarse
      ! Calculate the edge midpoint DOF
      DuFine(NVTcoarse + iedge) = 0.5_DP * (&
          DuCoarse(IverticesAtEdgeCoarse(1,iedge))+&
          DuCoarse(IverticesAtEdgeCoarse(2,iedge)))

    end do

    ! Loop over the faces
    do iface = 1, NATcoarse
      ! Calculate the face midpoint DOF
      DuFine(NVTCoarse + NMTCoarse + iface) = 0.25_DP * (&
          DuCoarse(IverticesAtFaceCoarse(1,iface))+&
          DuCoarse(IverticesAtFaceCoarse(2,iface))+&
          DuCoarse(IverticesAtFaceCoarse(3,iface))+&
          DuCoarse(IverticesAtFaceCoarse(4,iface)))
    end do

    ! Loop over the elements
    do iel = 1, NELcoarse
      ! Calculate the hexahedron cell midpoint DOF
      DuFine(NVTcoarse + NMTcoarse + NATcoarse + iel) = 0.125_DP * (&
          DuCoarse(IverticesAtElementCoarse(1,iel))+&
          DuCoarse(IverticesAtElementCoarse(2,iel))+&
          DuCoarse(IverticesAtElementCoarse(3,iel))+&
          DuCoarse(IverticesAtElementCoarse(4,iel))+&
          DuCoarse(IverticesAtElementCoarse(5,iel))+&
          DuCoarse(IverticesAtElementCoarse(6,iel))+&
          DuCoarse(IverticesAtElementCoarse(7,iel))+&
          DuCoarse(IverticesAtElementCoarse(8,iel)))
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_restUniformQ1_3D_double (DuCoarse,DuFine, &
             IverticesAtEdgeCoarse, IverticesAtFaceCoarse, &
             IverticesAtElementCoarse, NVTcoarse, NMTcoarse, NATcoarse,NELcoarse)

!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! Q1, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IverticesAtEdge array on the coarse grid.
  integer, dimension(:,:), intent(in) :: IverticesAtEdgeCoarse

  ! IverticesAtFace array on the coarse grid.
  integer, dimension(:,:), intent(in) :: IverticesAtFaceCoarse

  ! IverticesAtElement array (KVERT) on the coarse grid
  integer, dimension(:,:), intent(in) :: IverticesAtElementCoarse

  ! Number of vertices in the coarse grid
  integer, intent(in) :: NVTcoarse

  ! Number of edges in the coarse grid
  integer, intent(in) :: NMTcoarse

  ! Number of faces in the coarse grid
  integer, intent(in) :: NATcoarse

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  integer :: iel
  integer :: iedge, iface
  integer :: ivt
  real(DP) :: dx

    ! The information that was "distributed" in the prolongation has to
    ! be "collected".
    !
    ! Copy the first NVT entries - this gives the first additive contribution.
    call lalg_copyVector (DuFine(1:NVTcoarse),DuCoarse)

    ! Loop over the edges
    do iedge = 1, NMTcoarse
      ! get the fine grid DOF
      dx = 0.5_DP * DuFine(NVTcoarse + iedge)

      ! distribute it to the coarse grid DOFs
      ivt = IverticesAtEdgeCoarse(1,iedge)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
      ivt = IverticesAtEdgeCoarse(2,iedge)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
    end do

    ! Loop over the faces
    do iface = 1, NATcoarse
      ! get the fine grid DOF
      dx = 0.25_DP * DuFine(NVTcoarse + NMTcoarse + iface)

      ! distribute it to the coarse grid DOFs
      ivt = IverticesAtFaceCoarse(1,iface)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
      ivt = IverticesAtFaceCoarse(2,iface)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
      ivt = IverticesAtFaceCoarse(3,iface)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
      ivt = IverticesAtFaceCoarse(4,iface)
    end do

    ! Loop over the elements
    do iel = 1, NELcoarse
      ! get the fine grid DOF
      dx = 0.125_DP * DuFine(NVTcoarse + NMTcoarse + NATcoarse + iel)

      ! distribute it to the coarse grid DOFs
      ivt = IverticesAtElementCoarse(1, iel)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
      ivt = IverticesAtElementCoarse(2, iel)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
      ivt = IverticesAtElementCoarse(3, iel)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
      ivt = IverticesAtElementCoarse(4, iel)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
      ivt = IverticesAtElementCoarse(5, iel)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
      ivt = IverticesAtElementCoarse(6, iel)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
      ivt = IverticesAtElementCoarse(7, iel)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
      ivt = IverticesAtElementCoarse(8, iel)
      DuCoarse(ivt) = DuCoarse(ivt) + dx
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_interpUniformQ1_3D_double (DuCoarse,DuFine, NVTcoarse)

!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! Q1, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! Number of vertices in the coarse grid
  integer, intent(in) :: NVTcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

    ! The first coase.NVT entries of the fine grid vector define the values
    ! on the coarse grid - because of the two-level ordering!
    call lalg_copyVector(DUfine(1:NVTcoarse),DUcoarse(1:NVTCoarse))

  end subroutine


  ! ***************************************************************************
  ! Support for 3D Q1~ element
  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_prolUniformEx3x_3D_double (DuCoarse,DuFine, &
               IfacesAtElementCoarse,IfacesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               NELcoarse,celement)

!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! E030/E031/EM30/EM31, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Coarse grid vector
  real(DP), dimension(:), intent(in) :: DuCoarse

  ! IfacesAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IfacesAtElementCoarse

  ! IfacesAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IfacesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse

  ! Element type to use for prolongation
  integer(I32), intent(in) :: celement
!</input>

!<output>
  ! Fine grid vector
  real(DP), dimension(:), intent(out) :: DuFine
!</output>

!</subroutine>

  ! local variables
  real(DP) :: dw, A1, A2, A3, A4, A5, A6, A7
  real(DP), dimension(6) :: Dx
  integer, dimension(6) :: IM
  integer, dimension(4) :: idx
  integer :: iel
  integer, dimension(8) :: ielf

    ! Weights for the restriction; all coefficients are halfed, so dividing
    ! by 2 is not necessary in the calculation routines.
    if (iand(celement, int(2**16,I32)) .eq. 0) then
      ! Weights for Ex31 element
      A1 =  0.458333333333333_DP   ! = 11/24
      A2 =  0.145833333333333_DP   ! =  7/48
      A3 = -0.104166666666667_DP   ! = -5/48
      A4 = -0.041666666666667_DP   ! = -1/24
      A5 =  0.458333333333333_DP   ! = 11/24
      A6 =  0.083333333333333_DP   ! =  1/12
      A7 = -0.041666666666667_DP   ! = -1/24
    else
      ! Weights for Ex30 element
      A1 = 0.5_DP
      A2 = 0.125_DP
      A3 = -0.125_DP
      A4 = 0.0_DP
      A5 = 0.5_DP
      A6 = 0.0_DP
      A7 = 0.0_DP
    end if

    ! Clear the output vector
    call lalg_clearVector(DuFine)

    ! Loop over the coarse grid elements
    do iel = 1, NELcoarse

      ! Get the DOF`s of the coarse grid element
      IM(1) = IfacesAtElementCoarse(1,iel)
      IM(2) = IfacesAtElementCoarse(2,iel)
      IM(3) = IfacesAtElementCoarse(3,iel)
      IM(4) = IfacesAtElementCoarse(4,iel)
      IM(5) = IfacesAtElementCoarse(5,iel)
      IM(6) = IfacesAtElementCoarse(6,iel)

      ! Get the values of the corresponding DOF`s
      Dx(1) = DuCoarse(IM(1))
      Dx(2) = DuCoarse(IM(2))
      Dx(3) = DuCoarse(IM(3))
      Dx(4) = DuCoarse(IM(4))
      Dx(5) = DuCoarse(IM(5))
      Dx(6) = DuCoarse(IM(6))

      ! Get the element numbers of the fine-grid elements in the coarse grid;
      ! the numbers are defined by the two-level ordering.
      ielf(1) = iel
      ielf(2) = IneighboursAtElementFine(3,ielf(1))
      ielf(3) = IneighboursAtElementFine(3,ielf(2))
      ielf(4) = IneighboursAtElementFine(3,ielf(3))
      ielf(5) = IneighboursAtElementFine(6,ielf(1))
      ielf(6) = IneighboursAtElementFine(4,ielf(5))
      ielf(7) = IneighboursAtElementFine(4,ielf(6))
      ielf(8) = IneighboursAtElementFine(4,ielf(7))

      ! First, we are going to prolongate the DOFs on the boundary of
      ! the coarse grid hexahedron.
      ! Face 1
      ! Calculate scaling factor. If this face of the hexahedron belongs
      ! to the domain`s boundary (i.e. the hexahedron does not have a neighbour
      ! at this face), then the scaling factor is 2, otherwise 1.
      dw = 1.0_DP
      if (IneighboursAtElementCoarse(1,iel) .eq. 0) dw = 2.0_DP
      ! Calculate fine grid DOF indices for this coarse grid face
      idx(1) = IfacesAtElementFine(1, ielf(1))
      idx(2) = IfacesAtElementFine(1, ielf(2))
      idx(3) = IfacesAtElementFine(1, ielf(3))
      idx(4) = IfacesAtElementFine(1, ielf(4))
      ! distribute the DOF
      DuFine(idx(1))=DuFine(idx(1))+&
         dw*(A1*Dx(1)+A2*Dx(2)+A3*Dx(3)+A3*Dx(4)+A2*Dx(5)+A4*Dx(6))
      DuFine(idx(2))=DuFine(idx(2))+&
         dw*(A1*Dx(1)+A2*Dx(2)+A2*Dx(3)+A3*Dx(4)+A3*Dx(5)+A4*Dx(6))
      DuFine(idx(3))=DuFine(idx(3))+&
         dw*(A1*Dx(1)+A3*Dx(2)+A2*Dx(3)+A2*Dx(4)+A3*Dx(5)+A4*Dx(6))
      DuFine(idx(4))=DuFine(idx(4))+&
         dw*(A1*Dx(1)+A3*Dx(2)+A3*Dx(3)+A2*Dx(4)+A2*Dx(5)+A4*Dx(6))

      ! Face 2
      dw = 1.0_DP
      if (IneighboursAtElementCoarse(2,iel) .eq. 0) dw = 2.0_DP
      idx(1) = IfacesAtElementFine(2, ielf(1))
      idx(2) = IfacesAtElementFine(5, ielf(2))
      !idx(3) = IfacesAtElementFine(5, ielf(6))
      !idx(4) = IfacesAtElementFine(2, ielf(5))
      idx(3) = IfacesAtElementFine(2, ielf(6))
      idx(4) = IfacesAtElementFine(5, ielf(5))
      DuFine(idx(1))=DuFine(idx(1))+&
        dw*(A2*Dx(1)+A1*Dx(2)+A3*Dx(3)+A4*Dx(4)+A2*Dx(5)+A3*Dx(6))
      DuFine(idx(2))=DuFine(idx(2))+&
        dw*(A2*Dx(1)+A1*Dx(2)+A2*Dx(3)+A4*Dx(4)+A3*Dx(5)+A3*Dx(6))
      DuFine(idx(3))=DuFine(idx(3))+&
        dw*(A3*Dx(1)+A1*Dx(2)+A2*Dx(3)+A4*Dx(4)+A3*Dx(5)+A2*Dx(6))
      DuFine(idx(4))=DuFine(idx(4))+&
        dw*(A3*Dx(1)+A1*Dx(2)+A3*Dx(3)+A4*Dx(4)+A2*Dx(5)+A2*Dx(6))

      ! Face 3
      dw = 1.0_DP
      if (IneighboursAtElementCoarse(3,iel) .eq. 0) dw = 2.0_DP
      idx(1) = IfacesAtElementFine(2, ielf(2))
      idx(2) = IfacesAtElementFine(5, ielf(3))
      !idx(3) = IfacesAtElementFine(5, ielf(7))
      !idx(4) = IfacesAtElementFine(2, ielf(6))
      idx(3) = IfacesAtElementFine(2, ielf(7))
      idx(4) = IfacesAtElementFine(5, ielf(6))
      DuFine(idx(1))=DuFine(idx(1))+&
        dw*(A2*Dx(1)+A2*Dx(2)+A1*Dx(3)+A3*Dx(4)+A4*Dx(5)+A3*Dx(6))
      DuFine(idx(2))=DuFine(idx(2))+&
        dw*(A2*Dx(1)+A3*Dx(2)+A1*Dx(3)+A2*Dx(4)+A4*Dx(5)+A3*Dx(6))
      DuFine(idx(3))=DuFine(idx(3))+&
        dw*(A3*Dx(1)+A3*Dx(2)+A1*Dx(3)+A2*Dx(4)+A4*Dx(5)+A2*Dx(6))
      DuFine(idx(4))=DuFine(idx(4))+&
        dw*(A3*Dx(1)+A2*Dx(2)+A1*Dx(3)+A3*Dx(4)+A4*Dx(5)+A2*Dx(6))

      ! Face 4
      dw = 1.0_DP
      if (IneighboursAtElementCoarse(4,iel) .eq. 0) dw = 2.0_DP
      idx(1) = IfacesAtElementFine(2, ielf(3))
      idx(2) = IfacesAtElementFine(5, ielf(4))
      !idx(3) = IfacesAtElementFine(5, ielf(8))
      !idx(4) = IfacesAtElementFine(2, ielf(7))
      idx(3) = IfacesAtElementFine(2, ielf(8))
      idx(4) = IfacesAtElementFine(5, ielf(7))
      DuFine(idx(1))=DuFine(idx(1))+&
        dw*(A2*Dx(1)+A4*Dx(2)+A2*Dx(3)+A1*Dx(4)+A3*Dx(5)+A3*Dx(6))
      DuFine(idx(2))=DuFine(idx(2))+&
        dw*(A2*Dx(1)+A4*Dx(2)+A3*Dx(3)+A1*Dx(4)+A2*Dx(5)+A3*Dx(6))
      DuFine(idx(3))=DuFine(idx(3))+&
        dw*(A3*Dx(1)+A4*Dx(2)+A3*Dx(3)+A1*Dx(4)+A2*Dx(5)+A2*Dx(6))
      DuFine(idx(4))=DuFine(idx(4))+&
        dw*(A3*Dx(1)+A4*Dx(2)+A2*Dx(3)+A1*Dx(4)+A3*Dx(5)+A2*Dx(6))

      ! Face 5
      dw = 1.0_DP
      if (IneighboursAtElementCoarse(5,iel) .eq. 0) dw = 2.0_DP
      idx(1) = IfacesAtElementFine(2, ielf(4))
      idx(2) = IfacesAtElementFine(5, ielf(1))
      !idx(3) = IfacesAtElementFine(5, ielf(5))
      !idx(4) = IfacesAtElementFine(2, ielf(8))
      idx(3) = IfacesAtElementFine(2, ielf(5))
      idx(4) = IfacesAtElementFine(5, ielf(8))
      DuFine(idx(1))=DuFine(idx(1))+&
        dw*(A2*Dx(1)+A3*Dx(2)+A4*Dx(3)+A2*Dx(4)+A1*Dx(5)+A3*Dx(6))
      DuFine(idx(2))=DuFine(idx(2))+&
        dw*(A2*Dx(1)+A2*Dx(2)+A4*Dx(3)+A3*Dx(4)+A1*Dx(5)+A3*Dx(6))
      DuFine(idx(3))=DuFine(idx(3))+&
        dw*(A3*Dx(1)+A2*Dx(2)+A4*Dx(3)+A3*Dx(4)+A1*Dx(5)+A2*Dx(6))
      DuFine(idx(4))=DuFine(idx(4))+&
        dw*(A3*Dx(1)+A3*Dx(2)+A4*Dx(3)+A2*Dx(4)+A1*Dx(5)+A2*Dx(6))

      ! Face 6
      dw = 1.0_DP
      if (IneighboursAtElementCoarse(6,iel) .eq. 0) dw = 2.0_DP
      idx(1) = IfacesAtElementFine(1, ielf(5))
      idx(2) = IfacesAtElementFine(1, ielf(6))
      idx(3) = IfacesAtElementFine(1, ielf(7))
      idx(4) = IfacesAtElementFine(1, ielf(8))
      DuFine(idx(1))=DuFine(idx(1))+&
        dw*(A4*Dx(1)+A2*Dx(2)+A3*Dx(3)+A3*Dx(4)+A2*Dx(5)+A1*Dx(6))
      DuFine(idx(2))=DuFine(idx(2))+&
        dw*(A4*Dx(1)+A2*Dx(2)+A2*Dx(3)+A3*Dx(4)+A3*Dx(5)+A1*Dx(6))
      DuFine(idx(3))=DuFine(idx(3))+&
        dw*(A4*Dx(1)+A3*Dx(2)+A2*Dx(3)+A2*Dx(4)+A3*Dx(5)+A1*Dx(6))
      DuFine(idx(4))=DuFine(idx(4))+&
        dw*(A4*Dx(1)+A3*Dx(2)+A3*Dx(3)+A2*Dx(4)+A2*Dx(5)+A1*Dx(6))


      ! Now we need to calculate the DOFs of the fine grid which lie
      ! inside the coarse grid tetrahedron.
      idx(1) = IfacesAtElementFine(3, ielf(1))
      idx(2) = IfacesAtElementFine(3, ielf(2))
      idx(3) = IfacesAtElementFine(3, ielf(3))
      idx(4) = IfacesAtElementFine(3, ielf(4))
      DuFine(idx(1))=A5*Dx(1)+A5*Dx(2)+A6*Dx(3)+A7*Dx(4)+A6*Dx(5)+A7*Dx(6)
      DuFine(idx(2))=A5*Dx(1)+A6*Dx(2)+A5*Dx(3)+A6*Dx(4)+A7*Dx(5)+A7*Dx(6)
      DuFine(idx(3))=A5*Dx(1)+A7*Dx(2)+A6*Dx(3)+A5*Dx(4)+A6*Dx(5)+A7*Dx(6)
      DuFine(idx(4))=A5*Dx(1)+A6*Dx(2)+A7*Dx(3)+A6*Dx(4)+A5*Dx(5)+A7*Dx(6)

      idx(1) = IfacesAtElementFine(6, ielf(1))
      idx(2) = IfacesAtElementFine(6, ielf(2))
      idx(3) = IfacesAtElementFine(6, ielf(3))
      idx(4) = IfacesAtElementFine(6, ielf(4))
      DuFine(idx(1))=A6*Dx(1)+A5*Dx(2)+A7*Dx(3)+A7*Dx(4)+A5*Dx(5)+A6*Dx(6)
      DuFine(idx(2))=A6*Dx(1)+A5*Dx(2)+A5*Dx(3)+A7*Dx(4)+A7*Dx(5)+A6*Dx(6)
      DuFine(idx(3))=A6*Dx(1)+A7*Dx(2)+A5*Dx(3)+A5*Dx(4)+A7*Dx(5)+A6*Dx(6)
      DuFine(idx(4))=A6*Dx(1)+A7*Dx(2)+A7*Dx(3)+A5*Dx(4)+A5*Dx(5)+A6*Dx(6)

      !idx(1) = IfacesAtElementFine(3, ielf(5))
      !idx(2) = IfacesAtElementFine(3, ielf(6))
      !idx(3) = IfacesAtElementFine(3, ielf(7))
      !idx(4) = IfacesAtElementFine(3, ielf(8))
      idx(1) = IfacesAtElementFine(4, ielf(5))
      idx(2) = IfacesAtElementFine(4, ielf(6))
      idx(3) = IfacesAtElementFine(4, ielf(7))
      idx(4) = IfacesAtElementFine(4, ielf(8))
      DuFine(idx(1))=A7*Dx(1)+A5*Dx(2)+A6*Dx(3)+A7*Dx(4)+A6*Dx(5)+A5*Dx(6)
      DuFine(idx(2))=A7*Dx(1)+A6*Dx(2)+A5*Dx(3)+A6*Dx(4)+A7*Dx(5)+A5*Dx(6)
      DuFine(idx(3))=A7*Dx(1)+A7*Dx(2)+A6*Dx(3)+A5*Dx(4)+A6*Dx(5)+A5*Dx(6)
      DuFine(idx(4))=A7*Dx(1)+A6*Dx(2)+A7*Dx(3)+A6*Dx(4)+A5*Dx(5)+A5*Dx(6)

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_restUniformEx3x_3D_double (DuCoarse,DuFine, &
               IfacesAtElementCoarse,IfacesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               NELcoarse,celement)

!<description>
  ! Restrict a defect vector from a fine grid to a coarse grid.
  ! E031/EM31, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IfacesAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IfacesAtElementCoarse

  ! IfacesAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IfacesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse

  ! Element type to use for prolongation
  integer(I32), intent(in) :: celement
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  real(DP) :: dw, A1, A2, A3, A4, A5, A6, A7
  real(DP), dimension(36) :: Dx
  integer, dimension(6) :: IM
  integer :: iel
  integer, dimension(8) :: ielf

    ! Weights for the restriction; all coefficients are halfed, so dividing
    ! by 2 is not necessary in the calculation routines.
    if (iand(celement, int(2**16,I32)) .eq. 0) then
      ! Weights for Ex31 element
      A1 =  0.458333333333333_DP   ! = 11/24
      A2 =  0.145833333333333_DP   ! =  7/48
      A3 = -0.104166666666667_DP   ! = -5/48
      A4 = -0.041666666666667_DP   ! = -1/24
      A5 =  0.458333333333333_DP   ! = 11/24
      A6 =  0.083333333333333_DP   ! =  1/12
      A7 = -0.041666666666667_DP   ! = -1/24
    else
      ! Weights for Ex30 element
      A1 =  0.5_DP
      A2 =  0.125_DP
      A3 = -0.125_DP
      A4 =  0.0_DP
      A5 =  0.5_DP
      A6 =  0.0_DP
      A7 =  0.0_DP
    end if

    ! Clear the output vector
    call lalg_clearVector(DuCoarse)

    ! Loop over the coarse grid elements
    do iel = 1, NELcoarse

      ! Get the DOF`s of the coarse grid element
      IM(1) = IfacesAtElementCoarse(1,iel)
      IM(2) = IfacesAtElementCoarse(2,iel)
      IM(3) = IfacesAtElementCoarse(3,iel)
      IM(4) = IfacesAtElementCoarse(4,iel)
      IM(5) = IfacesAtElementCoarse(5,iel)
      IM(6) = IfacesAtElementCoarse(6,iel)

      ! Get the element numbers of the fine-grid elements in the coarse grid;
      ! the numbers are defined by the two-level ordering.
      ielf(1) = iel
      ielf(2) = IneighboursAtElementFine(3,ielf(1))
      ielf(3) = IneighboursAtElementFine(3,ielf(2))
      ielf(4) = IneighboursAtElementFine(3,ielf(3))
      ielf(5) = IneighboursAtElementFine(6,ielf(1))
      ielf(6) = IneighboursAtElementFine(4,ielf(5))
      ielf(7) = IneighboursAtElementFine(4,ielf(6))
      ielf(8) = IneighboursAtElementFine(4,ielf(7))

      ! Get the fine grid DOFs
      Dx( 1)=DuFine(IfacesAtElementFine(1,ielf(1)))
      Dx( 2)=DuFine(IfacesAtElementFine(1,ielf(2)))
      Dx( 3)=DuFine(IfacesAtElementFine(1,ielf(3)))
      Dx( 4)=DuFine(IfacesAtElementFine(1,ielf(4)))
      Dx( 5)=DuFine(IfacesAtElementFine(2,ielf(1)))
      Dx( 6)=DuFine(IfacesAtElementFine(5,ielf(2)))
      !Dx( 7)=DuFine(IfacesAtElementFine(5,ielf(6)))
      Dx( 7)=DuFine(IfacesAtElementFine(2,ielf(6)))
      !Dx( 8)=DuFine(IfacesAtElementFine(2,ielf(5)))
      Dx( 8)=DuFine(IfacesAtElementFine(5,ielf(5)))
      Dx( 9)=DuFine(IfacesAtElementFine(2,ielf(2)))
      Dx(10)=DuFine(IfacesAtElementFine(5,ielf(3)))
      !Dx(11)=DuFine(IfacesAtElementFine(5,ielf(7)))
      Dx(11)=DuFine(IfacesAtElementFine(2,ielf(7)))
      !Dx(12)=DuFine(IfacesAtElementFine(2,ielf(6)))
      Dx(12)=DuFine(IfacesAtElementFine(5,ielf(6)))
      Dx(13)=DuFine(IfacesAtElementFine(2,ielf(3)))
      Dx(14)=DuFine(IfacesAtElementFine(5,ielf(4)))
      !Dx(15)=DuFine(IfacesAtElementFine(5,ielf(8)))
      Dx(15)=DuFine(IfacesAtElementFine(2,ielf(8)))
      !Dx(16)=DuFine(IfacesAtElementFine(2,ielf(7)))
      Dx(16)=DuFine(IfacesAtElementFine(5,ielf(7)))
      Dx(17)=DuFine(IfacesAtElementFine(2,ielf(4)))
      Dx(18)=DuFine(IfacesAtElementFine(5,ielf(1)))
      !Dx(19)=DuFine(IfacesAtElementFine(5,ielf(5)))
      Dx(19)=DuFine(IfacesAtElementFine(2,ielf(5)))
      !Dx(20)=DuFine(IfacesAtElementFine(2,ielf(8)))
      Dx(20)=DuFine(IfacesAtElementFine(5,ielf(8)))
      Dx(21)=DuFine(IfacesAtElementFine(1,ielf(5)))
      Dx(22)=DuFine(IfacesAtElementFine(1,ielf(6)))
      Dx(23)=DuFine(IfacesAtElementFine(1,ielf(7)))
      Dx(24)=DuFine(IfacesAtElementFine(1,ielf(8)))
      Dx(25)=DuFine(IfacesAtElementFine(3,ielf(1)))
      Dx(26)=DuFine(IfacesAtElementFine(3,ielf(2)))
      Dx(27)=DuFine(IfacesAtElementFine(3,ielf(3)))
      Dx(28)=DuFine(IfacesAtElementFine(3,ielf(4)))
      Dx(29)=DuFine(IfacesAtElementFine(6,ielf(1)))
      Dx(30)=DuFine(IfacesAtElementFine(6,ielf(2)))
      Dx(31)=DuFine(IfacesAtElementFine(6,ielf(3)))
      Dx(32)=DuFine(IfacesAtElementFine(6,ielf(4)))
      !Dx(33)=DuFine(IfacesAtElementFine(3,ielf(5)))
      !Dx(34)=DuFine(IfacesAtElementFine(3,ielf(6)))
      !Dx(35)=DuFine(IfacesAtElementFine(3,ielf(7)))
      !Dx(36)=DuFine(IfacesAtElementFine(3,ielf(8)))
      Dx(33)=DuFine(IfacesAtElementFine(4,ielf(5)))
      Dx(34)=DuFine(IfacesAtElementFine(4,ielf(6)))
      Dx(35)=DuFine(IfacesAtElementFine(4,ielf(7)))
      Dx(36)=DuFine(IfacesAtElementFine(4,ielf(8)))

      ! Face 1
      dw = 1.0_DP
      if (IneighboursAtElementCoarse(1,iel) .eq. 0) dw = 2.0_DP
      DuCoarse(IM(1))= DuCoarse(IM(1)) + dw * (&
        +A1*(Dx(1)+Dx(2)+Dx(3)+Dx(4))&
        +A2*(Dx(5)+Dx(6)+Dx(9)+Dx(10)+Dx(13)+Dx(14)+Dx(17)+Dx(18))&
        +A3*(Dx(7)+Dx(8)+Dx(11)+Dx(12)+Dx(15)+Dx(16)+Dx(19)+Dx(20))&
        +A4*(Dx(21)+Dx(22)+Dx(23)+Dx(24))&
        +A5*(Dx(25)+Dx(26)+Dx(27)+Dx(28))&
        +A6*(Dx(29)+Dx(30)+Dx(31)+Dx(32))&
        +A7*(Dx(33)+Dx(34)+Dx(35)+Dx(36)))

      ! Face 2
      dw = 1.0_DP
      if (IneighboursAtElementCoarse(2,iel) .eq. 0) dw = 2.0_DP
      DuCoarse(IM(2))= DuCoarse(IM(2)) + dw * (&
        +A1*(Dx(5)+Dx(6)+Dx(7)+Dx(8))&
        +A2*(Dx(1)+Dx(2)+Dx(9)+Dx(12)+Dx(21)+Dx(22)+Dx(18)+Dx(19))&
        +A3*(Dx(3)+Dx(4)+Dx(10)+Dx(11)+Dx(23)+Dx(24)+Dx(17)+Dx(20))&
        +A4*(Dx(13)+Dx(14)+Dx(15)+Dx(16))&
        +A5*(Dx(25)+Dx(29)+Dx(30)+Dx(33))&
        +A6*(Dx(26)+Dx(28)+Dx(34)+Dx(36))&
        +A7*(Dx(27)+Dx(31)+Dx(32)+Dx(35)))

      ! Face 3
      dw = 1.0_DP
      if (IneighboursAtElementCoarse(3,iel) .eq. 0) dw = 2.0_DP
      DuCoarse(IM(3))= DuCoarse(IM(3)) + dw * (&
        +A1*(Dx(9)+Dx(10)+Dx(11)+Dx(12))&
        +A2*(Dx(2)+Dx(3)+Dx(6)+Dx(7)+Dx(22)+Dx(23)+Dx(13)+Dx(16))&
        +A3*(Dx(1)+Dx(4)+Dx(5)+Dx(8)+Dx(21)+Dx(24)+Dx(14)+Dx(15))&
        +A4*(Dx(17)+Dx(18)+Dx(19)+Dx(20))&
        +A5*(Dx(26)+Dx(30)+Dx(31)+Dx(34))&
        +A6*(Dx(25)+Dx(27)+Dx(33)+Dx(35))&
        +A7*(Dx(28)+Dx(29)+Dx(32)+Dx(36)))

      ! Face 4
      dw = 1.0_DP
      if (IneighboursAtElementCoarse(4,iel) .eq. 0) dw = 2.0_DP
      DuCoarse(IM(4))= DuCoarse(IM(4)) + dw * (&
        +A1*(Dx(13)+Dx(14)+Dx(15)+Dx(16))&
        +A2*(Dx(3)+Dx(4)+Dx(10)+Dx(11)+Dx(23)+Dx(24)+Dx(17)+Dx(20))&
        +A3*(Dx(1)+Dx(2)+Dx(9)+Dx(12)+Dx(21)+Dx(22)+Dx(18)+Dx(19))&
        +A4*(Dx(5)+Dx(6)+Dx(7)+Dx(8))&
        +A5*(Dx(27)+Dx(31)+Dx(32)+Dx(35))&
        +A6*(Dx(26)+Dx(28)+Dx(34)+Dx(36))&
        +A7*(Dx(25)+Dx(29)+Dx(30)+Dx(33)))

      ! Face 5
      dw = 1.0_DP
      if (IneighboursAtElementCoarse(5,iel) .eq. 0) dw = 2.0_DP
      DuCoarse(IM(5))= DuCoarse(IM(5)) + dw * (&
        +A1*(Dx(17)+Dx(18)+Dx(19)+Dx(20))&
        +A2*(Dx(1)+Dx(4)+Dx(5)+Dx(8)+Dx(21)+Dx(24)+Dx(14)+Dx(15))&
        +A3*(Dx(2)+Dx(3)+Dx(6)+Dx(7)+Dx(22)+Dx(23)+Dx(13)+Dx(16))&
        +A4*(Dx(9)+Dx(10)+Dx(11)+Dx(12))&
        +A5*(Dx(28)+Dx(29)+Dx(32)+Dx(36))&
        +A6*(Dx(25)+Dx(27)+Dx(33)+Dx(35))&
        +A7*(Dx(26)+Dx(30)+Dx(31)+Dx(34)))

      ! Face 6
      dw = 1.0_DP
      if (IneighboursAtElementCoarse(6,iel) .eq. 0) dw = 2.0_DP
      DuCoarse(IM(6))= DuCoarse(IM(6)) + dw * (&
        +A1*(Dx(21)+Dx(22)+Dx(23)+Dx(24))&
        +A2*(Dx(7)+Dx(8)+Dx(11)+Dx(12)+Dx(15)+Dx(16)+Dx(19)+Dx(20))&
        +A3*(Dx(5)+Dx(6)+Dx(9)+Dx(10)+Dx(13)+Dx(14)+Dx(17)+Dx(18))&
        +A4*(Dx(1)+Dx(2)+Dx(3)+Dx(4))&
        +A5*(Dx(33)+Dx(34)+Dx(35)+Dx(36))&
        +A6*(Dx(29)+Dx(30)+Dx(31)+Dx(32))&
        +A7*(Dx(25)+Dx(26)+Dx(27)+Dx(28)))

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_interpUniformEx3x_3D_dbl (DuCoarse,DuFine, &
               IfacesAtElementCoarse,IfacesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               NELcoarse,celement)

!<description>
  ! Restrict a solution vector from a fine grid to a coarse grid.
  ! E031/EM31, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IfacesAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IfacesAtElementCoarse

  ! IfacesAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IfacesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse

  ! Element type to use for prolongation
  integer(I32), intent(in) :: celement
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  real(DP) :: dw, A1, A2, A3
  real(DP), dimension(36) :: Dx
  integer, dimension(6) :: IM
  integer :: iel
  integer, dimension(8) :: ielf

    ! Weights for the restriction; all coefficients are halfed, so dividing
    ! by 2 is not necessary in the calculation routines.
!    IF (IAND(celement, 2**16) .EQ. 0) THEN
!      ! Weights for Ex31 element
!      A1 =  0.458333333333333_DP   ! = 11/24
!      A2 =  0.145833333333333_DP   ! =  7/48
!      A3 = -0.104166666666667_DP   ! = -5/48
!      A4 = -0.041666666666667_DP   ! = -1/24
!      A5 =  0.458333333333333_DP   ! = 11/24
!      A6 =  0.083333333333333_DP   ! =  1/12
!      A7 = -0.041666666666667_DP   ! = -1/24
!    ELSE
      ! Weights for Ex30 element
      A1 =  0.125_DP
      A2 =  0.25_DP
      A3 =  0.0833333333333333_DP
!    END IF

    ! Clear the output vector
    call lalg_clearVector(DuCoarse)

    ! Loop over the coarse grid elements
    do iel = 1, NELcoarse

      ! Get the DOF`s of the coarse grid element
      IM(1) = IfacesAtElementCoarse(1,iel)
      IM(2) = IfacesAtElementCoarse(2,iel)
      IM(3) = IfacesAtElementCoarse(3,iel)
      IM(4) = IfacesAtElementCoarse(4,iel)
      IM(5) = IfacesAtElementCoarse(5,iel)
      IM(6) = IfacesAtElementCoarse(6,iel)

      ! Get the element numbers of the fine-grid elements in the coarse grid;
      ! the numbers are defined by the two-level ordering.
      ielf(1) = iel
      ielf(2) = IneighboursAtElementFine(3,ielf(1))
      ielf(3) = IneighboursAtElementFine(3,ielf(2))
      ielf(4) = IneighboursAtElementFine(3,ielf(3))
      ielf(5) = IneighboursAtElementFine(6,ielf(1))
      ielf(6) = IneighboursAtElementFine(4,ielf(5))
      ielf(7) = IneighboursAtElementFine(4,ielf(6))
      ielf(8) = IneighboursAtElementFine(4,ielf(7))

      ! Get the fine grid DOFs
      Dx(1)=DuFine(IfacesAtElementFine(1,ielf(1)))
      Dx(2)=DuFine(IfacesAtElementFine(2,ielf(1)))
      Dx(3)=DuFine(IfacesAtElementFine(3,ielf(1)))
      Dx(4)=DuFine(IfacesAtElementFine(4,ielf(1)))
      Dx(5)=DuFine(IfacesAtElementFine(5,ielf(1)))
      Dx(6)=DuFine(IfacesAtElementFine(6,ielf(1)))
      Dx(7)=DuFine(IfacesAtElementFine(1,ielf(2)))
      Dx(8)=DuFine(IfacesAtElementFine(2,ielf(2)))
      Dx(9)=DuFine(IfacesAtElementFine(3,ielf(2)))
      Dx(10)=DuFine(IfacesAtElementFine(5,ielf(2)))
      Dx(11)=DuFine(IfacesAtElementFine(6,ielf(2)))
      Dx(12)=DuFine(IfacesAtElementFine(1,ielf(3)))
      Dx(13)=DuFine(IfacesAtElementFine(2,ielf(3)))
      Dx(14)=DuFine(IfacesAtElementFine(3,ielf(3)))
      Dx(15)=DuFine(IfacesAtElementFine(5,ielf(3)))
      Dx(16)=DuFine(IfacesAtElementFine(6,ielf(3)))
      Dx(17)=DuFine(IfacesAtElementFine(1,ielf(4)))
      Dx(18)=DuFine(IfacesAtElementFine(2,ielf(4)))
      Dx(19)=DuFine(IfacesAtElementFine(5,ielf(4)))
      Dx(20)=DuFine(IfacesAtElementFine(6,ielf(4)))
      Dx(21)=DuFine(IfacesAtElementFine(1,ielf(5)))
      !Dx(22)=DuFine(IfacesAtElementFine(2,ielf(5)))
      !Dx(23)=DuFine(IfacesAtElementFine(3,ielf(5)))
      !Dx(24)=DuFine(IfacesAtElementFine(4,ielf(5)))
      !Dx(25)=DuFine(IfacesAtElementFine(5,ielf(5)))
      Dx(22)=DuFine(IfacesAtElementFine(5,ielf(5)))
      Dx(23)=DuFine(IfacesAtElementFine(4,ielf(5)))
      Dx(24)=DuFine(IfacesAtElementFine(3,ielf(5)))
      Dx(25)=DuFine(IfacesAtElementFine(2,ielf(5)))
      Dx(26)=DuFine(IfacesAtElementFine(1,ielf(6)))
      !Dx(27)=DuFine(IfacesAtElementFine(2,ielf(6)))
      !Dx(28)=DuFine(IfacesAtElementFine(3,ielf(6)))
      !Dx(29)=DuFine(IfacesAtElementFine(5,ielf(6)))
      Dx(27)=DuFine(IfacesAtElementFine(5,ielf(6)))
      Dx(28)=DuFine(IfacesAtElementFine(4,ielf(6)))
      Dx(29)=DuFine(IfacesAtElementFine(2,ielf(6)))
      Dx(30)=DuFine(IfacesAtElementFine(1,ielf(7)))
      !Dx(31)=DuFine(IfacesAtElementFine(2,ielf(7)))
      !Dx(32)=DuFine(IfacesAtElementFine(3,ielf(7)))
      !Dx(33)=DuFine(IfacesAtElementFine(5,ielf(7)))
      Dx(31)=DuFine(IfacesAtElementFine(5,ielf(7)))
      Dx(32)=DuFine(IfacesAtElementFine(4,ielf(7)))
      Dx(33)=DuFine(IfacesAtElementFine(2,ielf(7)))
      Dx(34)=DuFine(IfacesAtElementFine(1,ielf(8)))
      !Dx(35)=DuFine(IfacesAtElementFine(2,ielf(8)))
      !Dx(36)=DuFine(IfacesAtElementFine(5,ielf(8)))
      Dx(35)=DuFine(IfacesAtElementFine(5,ielf(8)))
      Dx(36)=DuFine(IfacesAtElementFine(2,ielf(8)))

      ! Face 1
      dw = 1.0_DP
      if (IneighboursAtElementCoarse(1,iel) .eq. 0) dw = 2.0_DP
      DuCoarse(IM(1))= DuCoarse(IM(1)) + dw * (&
         A1*(Dx(1)+Dx(7)+Dx(12)+Dx(17))+&
         A2*(Dx(3)+Dx(14)+Dx(4)+Dx(9))-&
         A3*(Dx(2)+Dx(10)+Dx(8)+Dx(15)+Dx(13)+Dx(19)+&
             Dx(18)+Dx(5)+Dx(6)+Dx(11)+Dx(16)+Dx(20)))

      ! Face 2
      dw = 1.0_DP
      if (IneighboursAtElementCoarse(2,iel) .eq. 0) dw = 2.0_DP
      DuCoarse(IM(2))= DuCoarse(IM(2)) + dw * (&
          A1*(Dx(2)+Dx(10)+Dx(29)+Dx(22))+&
          A2*(Dx(3)+Dx(23)+Dx(6)+Dx(11))-&
          A3*(Dx(1)+Dx(7)+Dx(8)+Dx(27)+Dx(26)+Dx(21)+&
              Dx(25)+Dx(5)+Dx(4)+Dx(9)+Dx(24)+Dx(28)))

      ! Face 3
      dw = 1.0_DP
      if (IneighboursAtElementCoarse(3,iel) .eq. 0) dw = 2.0_DP
      DuCoarse(IM(3))= DuCoarse(IM(3)) + dw * (&
          A1*(Dx(8)+Dx(15)+Dx(33)+Dx(27))+&
          A2*(Dx(11)+Dx(16)+Dx(9)+Dx(28))-&
          A3*(Dx(3)+Dx(14)+Dx(32)+Dx(23)+Dx(7)+Dx(12)+&
              Dx(13)+Dx(31)+Dx(30)+Dx(26)+Dx(29)+Dx(10)))

      ! Face 4
      dw = 1.0_DP
      if (IneighboursAtElementCoarse(4,iel) .eq. 0) dw = 2.0_DP
      DuCoarse(IM(4))= DuCoarse(IM(4)) + dw * (&
        A1*(Dx(13)+Dx(31)+Dx(36)+Dx(19))+&
        A2*(Dx(32)+Dx(14)+Dx(16)+Dx(20))-&
        A3*(Dx(17)+Dx(12)+Dx(15)+Dx(33)+Dx(30)+Dx(34)+&
            Dx(35)+Dx(18)+Dx(4)+Dx(9)+Dx(28)+Dx(24)))

      ! Face 5
      dw = 1.0_DP
      if (IneighboursAtElementCoarse(5,iel) .eq. 0) dw = 2.0_DP
      DuCoarse(IM(5))= DuCoarse(IM(5)) + dw * (&
         A1*(Dx(5)+Dx(18)+Dx(35)+Dx(25))+&
         A2*(Dx(6)+Dx(20)+Dx(4)+Dx(24))-&
         A3*(Dx(1)+Dx(17)+Dx(19)+Dx(36)+Dx(34)+Dx(21)+&
             Dx(22)+Dx(2)+Dx(3)+Dx(14)+Dx(32)+Dx(23)))

      ! Face 6
      dw = 1.0_DP
      if (IneighboursAtElementCoarse(6,iel) .eq. 0) dw = 2.0_DP
      DuCoarse(IM(6))= DuCoarse(IM(6)) + dw *  (&
         A1*(Dx(26)+Dx(30)+Dx(34)+Dx(21))+&
         A2*(Dx(23)+Dx(32)+Dx(24)+Dx(28))-&
         A3*(Dx(22)+Dx(29)+Dx(27)+Dx(33)+Dx(31)+Dx(36)+&
             Dx(21)+Dx(25)+Dx(6)+Dx(11)+Dx(16)+Dx(20)))

    end do

  end subroutine

  ! ***************************************************************************
  ! Support for Q2~ element
  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_prolUniformEx50_double (DuCoarse,DuFine,&
               IedgesAtElementCoarse,IedgesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               ItwistCoarse,ItwistFine,NMTcoarse,NMTfine,NELcoarse)

!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! E050/EM50, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Coarse grid vector
  real(DP), dimension(:), intent(in) :: DuCoarse

  ! IedgesAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementCoarse

  ! IedgesAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! ItwistIndexEdges array on the coarse grid
  integer(I32), dimension(:), intent(in) :: ItwistCoarse

  ! ItwistIndexEdges array on the fine grid
  integer(I32), dimension(:), intent(in) :: ItwistFine

  ! Number of egdes in the coarse grid
  integer, intent(in) :: NMTcoarse

  ! Number of egdes in the fine grid
  integer, intent(in) :: NMTfine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Fine grid vector
  real(DP), dimension(:), intent(out) :: DuFine
!</output>

!</subroutine>

  ! local variables
  integer :: iel,i
  integer, dimension(4) :: Ielf, Imc

  ! Local vectors
  real(DP), dimension(9) :: Dv
  real(DP), dimension(4:6,4) :: Dtf
  real(DP), dimension(4) :: Dtn
  integer, dimension(7,4) :: Idf

  ! "Projection" coefficients
  real(DP), parameter :: R11 =   27.0_DP / 808.0_DP
  real(DP), parameter :: R12 = 1239.0_DP / 808.0_DP

  real(DP), parameter :: R21 =   81.0_DP / 808.0_DP
  real(DP), parameter :: R22 =  222.0_DP / 808.0_DP
  real(DP), parameter :: R23 =  323.0_DP / 808.0_DP
  real(DP), parameter :: R24 =  384.0_DP / 808.0_DP

    ! Clear the output vector
    call lalg_clearVector(DuFine)

    ! Loop over the coarse grid elements
    do iel = 1, NELcoarse

      ! Get the element numbers of the fine grid elements
      Ielf(1) = iel
      Ielf(2) = IneighboursAtElementFine(2,Ielf(1))
      Ielf(3) = IneighboursAtElementFine(2,Ielf(2))
      Ielf(4) = IneighboursAtElementFine(2,Ielf(3))

      ! Get the edge indices of the coarse grid element
      Imc(1) = IedgesAtElementCoarse(1,iel)
      Imc(2) = IedgesAtElementCoarse(2,iel)
      Imc(3) = IedgesAtElementCoarse(3,iel)
      Imc(4) = IedgesAtElementCoarse(4,iel)

      Dtn = 1.0_DP
      do i = 1, 4

        ! Calculate the Edge-Int-Mean DOFs on the fine mesh
        Idf(1,i) = IedgesAtElementFine(1,Ielf(i))
        Idf(2,i) = IedgesAtElementFine(2,Ielf(i))
        Idf(3,i) = IedgesAtElementFine(4,Ielf(i))
        ! Calculate the Edge-Legendre-Int-Mean DOFs on the fine mesh
        Idf(4,i) = Idf(1,i) + NMTfine
        Idf(5,i) = Idf(2,i) + NMTfine
        Idf(6,i) = Idf(3,i) + NMTfine
        ! Calculate the Quad-Int-Mean DOF on the fine mesh
        Idf(7,i) = 2*NMTfine + Ielf(i)

        ! Calculate the twist factors for the Edge-Legendre-Int-Mean DOFs
        ! on the fine mesh
        Dtf(4,i) = real(1-2*iand(      ItwistFine(Ielf(i))    ,1_I32),DP)
        Dtf(5,i) = real(1-2*iand(ishft(ItwistFine(Ielf(i)),-1),1_I32),DP)
        Dtf(6,i) = real(1-2*iand(ishft(ItwistFine(Ielf(i)),-3),1_I32),DP)

        ! If there is a neighbour at the coarse mesh edge, then we
        ! set the corresponding factor to 1/2
        if (IneighboursAtElementCoarse(i,iel) .ne. 0) Dtn(i) = 0.5_DP

      end do

      ! Get the values of the corresponding coarse mesh DOFs
      Dv(1:4) = DuCoarse(Imc(1:4))
      Dv(  5) = DuCoarse(Imc(1) + NMTcoarse) &
              * real(1-2*iand(      ItwistCoarse(iel)    ,1_I32),DP)
      Dv(  6) = DuCoarse(Imc(2) + NMTcoarse) &
              * real(1-2*iand(ishft(ItwistCoarse(iel),-1),1_I32),DP)
      Dv(  7) = DuCoarse(Imc(3) + NMTcoarse) &
              * real(1-2*iand(ishft(ItwistCoarse(iel),-2),1_I32),DP)
      Dv(  8) = DuCoarse(Imc(4) + NMTcoarse) &
              * real(1-2*iand(ishft(ItwistCoarse(iel),-3),1_I32),DP)
      Dv(  9) = DuCoarse(2*NMTcoarse + iel)

      DuFine(Idf(1,1)) = DuFine(Idf(1,1)) + (Dv(1)-R12*Dv(5)+&
                         R11*(-Dv(6)-Dv(7)-Dv(8)))*Dtn(1)
      DuFine(Idf(3,2)) = DuFine(Idf(3,2)) + (Dv(1)+R12*Dv(5)+&
                         R11*(Dv(6)+Dv(7)+Dv(8)))*Dtn(1)
      DuFine(Idf(1,2)) = DuFine(Idf(1,2)) + (Dv(2)-R12*Dv(6)+&
                         R11*(-Dv(5)-Dv(7)-Dv(8)))*Dtn(2)
      DuFine(Idf(3,3)) = DuFine(Idf(3,3)) + (Dv(2)+R12*Dv(6)+&
                         R11*(Dv(5)+Dv(7)+Dv(8)))*Dtn(2)
      DuFine(Idf(1,3)) = DuFine(Idf(1,3)) + (Dv(3)-R12*Dv(7)+&
                         R11*(-Dv(5)-Dv(6)-Dv(8)))*Dtn(3)
      DuFine(Idf(3,4)) = DuFine(Idf(3,4)) + (Dv(3)+R12*Dv(7)+&
                         R11*(Dv(5)+Dv(6)+Dv(8)))*Dtn(3)
      DuFine(Idf(1,4)) = DuFine(Idf(1,4)) + (Dv(4)-R12*Dv(8)+&
                         R11*(-Dv(5)-Dv(6)-Dv(7)))*Dtn(4)
      DuFine(Idf(3,1)) = DuFine(Idf(3,1)) + (Dv(4)+R12*Dv(8)+&
                         R11*(Dv(5)+Dv(6)+Dv(7)))*Dtn(4)
      DuFine(Idf(2,1)) = DuFine(Idf(2,1)) - 0.25_DP*(Dv(2)+Dv(4))+&
                         0.375_DP*(Dv(1)-Dv(3)+Dv(6)-Dv(8))+1.5_DP*Dv(9)
      DuFine(Idf(2,2)) = DuFine(Idf(2,2)) - 0.25_DP*(Dv(1)+Dv(3))+&
                         0.375_DP*(Dv(2)-Dv(4)-Dv(5)+Dv(7))+1.5_DP*Dv(9)
      DuFine(Idf(2,3)) = DuFine(Idf(2,3)) - 0.25_DP*(Dv(2)+Dv(4))+&
                         0.375_DP*(-Dv(1)+Dv(3)-Dv(6)+Dv(8))+1.5_DP*Dv(9)
      DuFine(Idf(2,4)) = DuFine(Idf(2,4)) - 0.25_DP*(Dv(1)+Dv(3))+&
                         0.375_DP*(-Dv(2)+Dv(4)+Dv(5)-Dv(7))+1.5_DP*Dv(9)

      DuFine(Idf(4,1)) = DuFine(Idf(4,1)) + Dtn(1)*Dtf(4,1)*(&
           0.125_DP*(Dv(1)-Dv(2)-Dv(3)-Dv(4))+R23*Dv(5)+R22*Dv(6)&
           -R21*Dv(7)-R24*Dv(8)+0.25_DP*Dv(9))
      DuFine(Idf(6,2)) = DuFine(Idf(6,2)) + Dtn(1)*Dtf(6,2)*(&
           0.125_DP*(-Dv(1)+Dv(2)+Dv(3)+Dv(4))+R23*Dv(5)-R24*Dv(6)&
           -R21*Dv(7)+R22*Dv(8)-0.25_DP*Dv(9))
      DuFine(Idf(4,2)) = DuFine(Idf(4,2)) + Dtn(2)*Dtf(4,2)*(&
           0.125_DP*(-Dv(1)+Dv(2)-Dv(3)-Dv(4))-R24*Dv(5)+R23*Dv(6)&
           +R22*Dv(7)-R21*Dv(8)+0.25_DP*Dv(9))
      DuFine(Idf(6,3)) = DuFine(Idf(6,3)) + Dtn(2)*Dtf(6,3)*(&
           0.125_DP*(Dv(1)-Dv(2)+Dv(3)+Dv(4))+R22*Dv(5)+R23*Dv(6)&
           -R24*Dv(7)-R21*Dv(8)-0.25_DP*Dv(9))
      DuFine(Idf(4,3)) = DuFine(Idf(4,3)) + Dtn(3)*Dtf(4,3)*(&
           0.125_DP*(-Dv(1)-Dv(2)+Dv(3)-Dv(4))-R21*Dv(5)-R24*Dv(6)&
           +R23*Dv(7)+R22*Dv(8)+0.25_DP*Dv(9))
      DuFine(Idf(6,4)) = DuFine(Idf(6,4)) + Dtn(3)*Dtf(6,4)*(&
           0.125_DP*(Dv(1)+Dv(2)-Dv(3)+Dv(4))-R21*Dv(5)+R22*Dv(6)&
           +R23*Dv(7)-R24*Dv(8)-0.25_DP*Dv(9))
      DuFine(Idf(4,4)) = DuFine(Idf(4,4)) + Dtn(4)*Dtf(4,4)*(&
           0.125_DP*(-Dv(1)-Dv(2)-Dv(3)+Dv(4))+R22*Dv(5)-R21*Dv(6)&
           -R24*Dv(7)+R23*Dv(8)+0.25_DP*Dv(9))
      DuFine(Idf(6,1)) = DuFine(Idf(6,1)) + Dtn(4)*Dtf(6,1)*(&
           0.125_DP*(Dv(1)+Dv(2)+Dv(3)-Dv(4))-R24*Dv(5)-R21*Dv(6)&
           +R22*Dv(7)+R23*Dv(8)-0.25_DP*Dv(9))

      DuFine(Idf(5,1)) = DuFine(Idf(5,1)) + Dtf(5,1)*(0.25_DP*(-Dv(1)+Dv(9))&
                         +0.125_DP*(-Dv(6)+Dv(8)))
      DuFine(Idf(5,2)) = DuFine(Idf(5,2)) + Dtf(5,2)*(0.25_DP*(-Dv(2)+Dv(9))&
                         +0.125_DP*(Dv(5)-Dv(7)))
      DuFine(Idf(5,3)) = DuFine(Idf(5,3)) + Dtf(5,3)*(0.25_DP*(-Dv(3)+Dv(9))&
                         +0.125_DP*(Dv(6)-Dv(8)))
      DuFine(Idf(5,4)) = DuFine(Idf(5,4)) + Dtf(5,4)*(0.25_DP*(-Dv(4)+Dv(9))&
                         +0.125_DP*(-Dv(5)+Dv(7)))

      DuFine(Idf(7,1)) = DuFine(Idf(7,1)) + 0.25_DP*(Dv(1)-Dv(2)-Dv(3)+Dv(4))+&
                         0.1875_DP*(-Dv(5)+Dv(6)-Dv(7)+Dv(8))+Dv(9)
      DuFine(Idf(7,2)) = DuFine(Idf(7,2)) + 0.25_DP*(Dv(1)+Dv(2)-Dv(3)-Dv(4))+&
                         0.1875_DP*(Dv(5)-Dv(6)+Dv(7)-Dv(8))+Dv(9)
      DuFine(Idf(7,3)) = DuFine(Idf(7,3)) + 0.25_DP*(-Dv(1)+Dv(2)+Dv(3)-Dv(4))+&
                         0.1875_DP*(-Dv(5)+Dv(6)-Dv(7)+Dv(8))+Dv(9)
      DuFine(Idf(7,4)) = DuFine(Idf(7,4)) + 0.25_DP*(-Dv(1)-Dv(2)+Dv(3)+Dv(4))+&
                         0.1875_DP*(Dv(5)-Dv(6)+Dv(7)-Dv(8))+Dv(9)

      ! Go for the next element

    end do

    ! That is it

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_restUniformEx50_double (DuCoarse,DuFine, &
               IedgesAtElementCoarse,IedgesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               ItwistCoarse,ItwistFine,NMTcoarse,NMTfine,NELcoarse)

!<description>
  ! Restricts a defect vector from a fine grid to a coarse grid.
  ! E050/EM50, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IedgesAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementCoarse

  ! IedgesAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! ItwistIndexEdges array on the coarse grid
  integer(I32), dimension(:), intent(in) :: ItwistCoarse

  ! ItwistIndexEdges array on the fine grid
  integer(I32), dimension(:), intent(in) :: ItwistFine

  ! Number of egdes in the coarse grid
  integer, intent(in) :: NMTcoarse

  ! Number of egdes in the fine grid
  integer, intent(in) :: NMTfine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  integer :: iel,i
  integer, dimension(4) :: Ielf

  ! local vectors
  real(DP), dimension(28) :: dv
  real(DP), dimension(4) :: Dtn,Dtf
  integer, dimension(9) :: Idf

  ! "Projection" coefficients
  real(DP), parameter :: R11 =   27.0_DP / 808.0_DP
  real(DP), parameter :: R12 = 1239.0_DP / 808.0_DP

  real(DP), parameter :: R21 =   81.0_DP / 808.0_DP
  real(DP), parameter :: R22 =  222.0_DP / 808.0_DP
  real(DP), parameter :: R23 =  323.0_DP / 808.0_DP
  real(DP), parameter :: R24 =  384.0_DP / 808.0_DP

    ! Clear the output vector
    call lalg_clearVector(DuCoarse)

    ! Loop over the coarse grid elements
    do iel = 1, NELcoarse

      ! Get the element numbers of the fine grid elements
      Ielf(1) = iel
      Ielf(2) = IneighboursAtElementFine(2,Ielf(1))
      Ielf(3) = IneighboursAtElementFine(2,Ielf(2))
      Ielf(4) = IneighboursAtElementFine(2,Ielf(3))

      ! Calculate the coarse grid DOFs
      Idf(1:4) = IedgesAtElementCoarse(1:4,iel)
      Idf(5:8) = Idf(1:4) + NMTcoarse
      Idf(  9) = 2*NMTcoarse + iel

      ! Calculate twist index factors and neighbour scales for the
      ! coarse grid edges
      Dtn = 1.0_DP
      do i = 1, 4

        ! Twist index factor
        Dtf(i) = real(1-2*iand(ishft(ItwistCoarse(iel),1-i),1_I32),DP)

        ! Neighbour factor
        if (IneighboursAtElementCoarse(i,iel) .ne. 0) Dtn(i) = 0.5_DP

      end do

      ! Get the values of the corresponding fine mesh DOFs
      Dv( 1) = DuFine(IedgesAtElementFine(1,Ielf(1)))*Dtn(1)
      Dv( 2) = DuFine(IedgesAtElementFine(4,Ielf(2)))*Dtn(1)
      Dv( 3) = DuFine(IedgesAtElementFine(1,Ielf(2)))*Dtn(2)
      Dv( 4) = DuFine(IedgesAtElementFine(4,Ielf(3)))*Dtn(2)
      Dv( 5) = DuFine(IedgesAtElementFine(1,Ielf(3)))*Dtn(3)
      Dv( 6) = DuFine(IedgesAtElementFine(4,Ielf(4)))*Dtn(3)
      Dv( 7) = DuFine(IedgesAtElementFine(1,Ielf(4)))*Dtn(4)
      Dv( 8) = DuFine(IedgesAtElementFine(4,Ielf(1)))*Dtn(4)
      Dv( 9) = DuFine(IedgesAtElementFine(2,Ielf(1)))
      Dv(10) = DuFine(IedgesAtElementFine(2,Ielf(2)))
      Dv(11) = DuFine(IedgesAtElementFine(2,Ielf(3)))
      Dv(12) = DuFine(IedgesAtElementFine(2,Ielf(4)))
      Dv(13) = DuFine(IedgesAtElementFine(1,Ielf(1))+NMTfine)&
             * real(1-2*iand(      ItwistFine(Ielf(1))    ,1_I32),DP)*Dtn(1)
      Dv(14) = DuFine(IedgesAtElementFine(4,Ielf(2))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(2)),-3),1_I32),DP)*Dtn(1)
      Dv(15) = DuFine(IedgesAtElementFine(1,Ielf(2))+NMTfine)&
             * real(1-2*iand(      ItwistFine(Ielf(2))    ,1_I32),DP)*Dtn(2)
      Dv(16) = DuFine(IedgesAtElementFine(4,Ielf(3))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(3)),-3),1_I32),DP)*Dtn(2)
      Dv(17) = DuFine(IedgesAtElementFine(1,Ielf(3))+NMTfine)&
             * real(1-2*iand(      ItwistFine(Ielf(3))    ,1_I32),DP)*Dtn(3)
      Dv(18) = DuFine(IedgesAtElementFine(4,Ielf(4))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(4)),-3),1_I32),DP)*Dtn(3)
      Dv(19) = DuFine(IedgesAtElementFine(1,Ielf(4))+NMTfine)&
             * real(1-2*iand(      ItwistFine(Ielf(4))    ,1_I32),DP)*Dtn(4)
      Dv(20) = DuFine(IedgesAtElementFine(4,Ielf(1))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(1)),-3),1_I32),DP)*Dtn(4)
      Dv(21) = DuFine(IedgesAtElementFine(2,Ielf(1))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(1)),-1),1_I32),DP)
      Dv(22) = DuFine(IedgesAtElementFine(2,Ielf(2))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(2)),-1),1_I32),DP)
      Dv(23) = DuFine(IedgesAtElementFine(2,Ielf(3))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(3)),-1),1_I32),DP)
      Dv(24) = DuFine(IedgesAtElementFine(2,Ielf(4))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(4)),-1),1_I32),DP)
      Dv(25) = DuFine(2*NMTfine + Ielf(1))
      Dv(26) = DuFine(2*NMTfine + Ielf(2))
      Dv(27) = DuFine(2*NMTfine + Ielf(3))
      Dv(28) = DuFine(2*NMTfine + Ielf(4))

      DuCoarse(Idf(1)) = DuCoarse(Idf(1)) + 0.125_DP*(+Dv(13)-Dv(14)-Dv(15)+Dv(16)&
                        -Dv(17)+Dv(18)-Dv(19)+Dv(20))+0.25_DP*(-Dv(21)+Dv(25)&
                        +Dv(26)-Dv(27)-Dv(28)-Dv(10)-Dv(12))+Dv(1)+Dv(2)&
                        +0.375_DP*(Dv(9)-Dv(11))
      DuCoarse(Idf(2)) = DuCoarse(Idf(2)) + 0.125_DP*(-Dv(13)+Dv(14)+Dv(15)-Dv(16)&
                        -Dv(17)+Dv(18)-Dv(19)+Dv(20))+0.25_DP*(-Dv(22)-Dv(25)&
                        +Dv(26)+Dv(27)-Dv(28)-Dv(9)-Dv(11))+Dv(3)+Dv(4)&
                        +0.375_DP*(Dv(10)-Dv(12))
      DuCoarse(Idf(3)) = DuCoarse(Idf(3)) + 0.125_DP*(-Dv(13)+Dv(14)-Dv(15)+Dv(16)&
                        +Dv(17)-Dv(18)-Dv(19)+Dv(20))+0.25_DP*(-Dv(23)-Dv(25)&
                        -Dv(26)+Dv(27)+Dv(28)-Dv(10)-Dv(12))+Dv(5)+Dv(6)&
                        +0.375_DP*(-Dv(9)+Dv(11))
      DuCoarse(Idf(4)) = DuCoarse(Idf(4)) + 0.125_DP*(-Dv(13)+Dv(14)-Dv(15)+Dv(16)&
                        -Dv(17)+Dv(18)+Dv(19)-Dv(20))+0.25_DP*(-Dv(24)+Dv(25)&
                        -Dv(26)-Dv(27)+Dv(28)-Dv(9)-Dv(11))+Dv(7)+Dv(8)&
                        +0.375_DP*(-Dv(10)+Dv(12))

      DuCoarse(Idf(5)) = DuCoarse(Idf(5)) + Dtf(1)*(R24*(-Dv(15)-Dv(20))&
                        +R22*(Dv(16)+Dv(19))+0.125_DP*(Dv(22)-Dv(24))&
                        +0.1875_DP*(-Dv(25)+Dv(26)-Dv(27)+Dv(28))+R12*(-Dv(1)+Dv(2))&
                        +R11*(-Dv(3)+Dv(4)-Dv(5)+Dv(6)-Dv(7)+Dv(8))&
                        +R21*(-Dv(17)-Dv(18))&
                        +0.375_DP*(-Dv(10)+Dv(12))+R23*(Dv(13)+Dv(14)))
      DuCoarse(Idf(6)) = DuCoarse(Idf(6)) + Dtf(2)*(R24*(-Dv(14)-Dv(17))&
                        +R22*(Dv(13)+Dv(18))+0.125_DP*(-Dv(21)+Dv(23))&
                        +0.1875_DP*(Dv(25)-Dv(26)+Dv(27)-Dv(28))+R12*(-Dv(3)+Dv(4))&
                        +R11*(-Dv(1)+Dv(2)-Dv(5)+Dv(6)-Dv(7)+Dv(8))&
                        +R21*(-Dv(19)-Dv(20))&
                        +0.375_DP*(Dv(9)-Dv(11))+R23*(Dv(15)+Dv(16)))
      DuCoarse(Idf(7)) = DuCoarse(Idf(7)) + Dtf(3)*(R24*(-Dv(16)-Dv(19))&
                        +R22*(Dv(15)+Dv(20))+0.125_DP*(Dv(24)-Dv(22))&
                        +0.1875_DP*(-Dv(25)+Dv(26)-Dv(27)+Dv(28))+R12*(-Dv(5)+Dv(6))&
                        +R11*(-Dv(1)+Dv(2)-Dv(3)+Dv(4)-Dv(7)+Dv(8))&
                        +R21*(-Dv(13)-Dv(14))&
                        +0.375_DP*(Dv(10)-Dv(12))+R23*(Dv(17)+Dv(18)))
      DuCoarse(Idf(8)) = DuCoarse(Idf(8)) + Dtf(4)*(R24*(-Dv(13)-Dv(18))&
                        +R22*(Dv(14)+Dv(17))+0.125_DP*(Dv(21)-Dv(23))&
                        +0.1875_DP*(Dv(25)-Dv(26)+Dv(27)-Dv(28))+R12*(-Dv(7)+Dv(8))&
                        +R11*(-Dv(1)+Dv(2)-Dv(3)+Dv(4)-Dv(5)+Dv(6))&
                        +R21*(-Dv(15)-Dv(16))&
                        +0.375_DP*(-Dv(9)+Dv(11))+R23*(Dv(19)+Dv(20)))

      DuCoarse(Idf(9)) = DuCoarse(Idf(9)) + 0.25_DP*(Dv(13)-Dv(14)+Dv(15)-Dv(16)&
                        +Dv(17)-Dv(18)+Dv(19)-Dv(20)+Dv(21)+Dv(22)+Dv(23)+Dv(24))&
                        +Dv(25)+Dv(26)+Dv(27)+Dv(28)+1.5_DP*(Dv(9)+Dv(10)+Dv(11)+Dv(12))

    end do

    ! That is it

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_interpUniformEx50_double (DuCoarse,DuFine, &
               IedgesAtElementCoarse,IedgesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               ItwistCoarse,ItwistFine,NMTcoarse,NMTfine,NELcoarse)

!<description>
  ! Restricts a solution vector from a fine grid to a coarse grid.
  ! E050/EM50, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IedgesAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementCoarse

  ! IedgesAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! ItwistIndexEdges array on the coarse grid
  integer(I32), dimension(:), intent(in) :: ItwistCoarse

  ! ItwistIndexEdges array on the fine grid
  integer(I32), dimension(:), intent(in) :: ItwistFine

  ! Number of egdes in the coarse grid
  integer, intent(in) :: NMTcoarse

  ! Number of egdes in the fine grid
  integer, intent(in) :: NMTfine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  integer :: iel,i
  integer, dimension(4) :: Ielf

  ! local vectors
  real(DP), dimension(28) :: dv
  real(DP), dimension(4) :: Dtn,Dtf
  integer, dimension(9) :: Idf

  ! "Projection" coefficients
  real(DP), parameter :: Q11 =  1.0_DP /  8.0_DP    ! 0.125
  real(DP), parameter :: Q12 =  1.0_DP /  4.0_DP    ! 0.25
  real(DP), parameter :: Q13 =  3.0_DP / 16.0_DP    ! 0.1875
  real(DP), parameter :: Q14 =  9.0_DP / 16.0_DP    ! 0.5625

  real(DP), parameter :: Q21 = 19.0_DP /  96.0_DP   ! 0.197916666...
  real(DP), parameter :: Q22 =  1.0_DP /   8.0_DP   ! 0.125
  real(DP), parameter :: Q23 =  5.0_DP / 192.0_DP   ! 0.026041666...
  real(DP), parameter :: Q24 =  1.0_DP /  96.0_DP   ! 0.010416666...
  real(DP), parameter :: Q25 = 11.0_DP / 192.0_DP   ! 0.057291666...
  real(DP), parameter :: Q26 = 13.0_DP / 128.0_DP   ! 0.1015625
  real(DP), parameter :: Q27 =  3.0_DP /  64.0_DP   ! 0.046875
  real(DP), parameter :: Q28 =  3.0_DP / 128.0_DP   ! 0.0234375
  real(DP), parameter :: Q29 =  1.0_DP / 128.0_DP   ! 0.0078125
  real(DP), parameter :: Q30 =  7.0_DP / 128.0_DP   ! 0.0546875
  real(DP), parameter :: Q31 =  9.0_DP /  32.0_DP
  real(DP), parameter :: Q32 =  3.0_DP /  32.0_DP

    ! Clear the output vector
    call lalg_clearVector(DuCoarse)

    ! Loop over the coarse grid elements
    do iel = 1, NELcoarse

      ! Get the element numbers of the fine grid elements
      Ielf(1) = iel
      Ielf(2) = IneighboursAtElementFine(2,Ielf(1))
      Ielf(3) = IneighboursAtElementFine(2,Ielf(2))
      Ielf(4) = IneighboursAtElementFine(2,Ielf(3))

      ! Calculate the coarse grid DOFs
      Idf(1:4) = IedgesAtElementCoarse(1:4,iel)
      Idf(5:8) = Idf(1:4) + NMTcoarse
      Idf(  9) = 2*NMTcoarse + iel

      ! Calculate twist index factors and neighbour scales for the
      ! coarse grid edges
      Dtn = 1.0_DP
      do i = 1, 4

        ! Twist index factor
        Dtf(i) = real(1-2*iand(ishft(ItwistCoarse(iel),1-i),1_I32),DP)

        ! Neighbour factor
        if (IneighboursAtElementCoarse(i,iel) .ne. 0) Dtn(i) = 0.5_DP

      end do

      ! Get the values of the corresponding fine mesh DOFs
      Dv( 1) = DuFine(IedgesAtElementFine(1,Ielf(1)))
      Dv( 2) = DuFine(IedgesAtElementFine(4,Ielf(2)))
      Dv( 3) = DuFine(IedgesAtElementFine(1,Ielf(2)))
      Dv( 4) = DuFine(IedgesAtElementFine(4,Ielf(3)))
      Dv( 5) = DuFine(IedgesAtElementFine(1,Ielf(3)))
      Dv( 6) = DuFine(IedgesAtElementFine(4,Ielf(4)))
      Dv( 7) = DuFine(IedgesAtElementFine(1,Ielf(4)))
      Dv( 8) = DuFine(IedgesAtElementFine(4,Ielf(1)))
      Dv( 9) = DuFine(IedgesAtElementFine(2,Ielf(1)))
      Dv(10) = DuFine(IedgesAtElementFine(2,Ielf(2)))
      Dv(11) = DuFine(IedgesAtElementFine(2,Ielf(3)))
      Dv(12) = DuFine(IedgesAtElementFine(2,Ielf(4)))
      Dv(13) = DuFine(IedgesAtElementFine(1,Ielf(1))+NMTfine)&
             * real(1-2*iand(      ItwistFine(Ielf(1))    ,1_I32),DP)
      Dv(14) = DuFine(IedgesAtElementFine(4,Ielf(2))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(2)),-3),1_I32),DP)
      Dv(15) = DuFine(IedgesAtElementFine(1,Ielf(2))+NMTfine)&
             * real(1-2*iand(      ItwistFine(Ielf(2))    ,1_I32),DP)
      Dv(16) = DuFine(IedgesAtElementFine(4,Ielf(3))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(3)),-3),1_I32),DP)
      Dv(17) = DuFine(IedgesAtElementFine(1,Ielf(3))+NMTfine)&
             * real(1-2*iand(      ItwistFine(Ielf(3))    ,1_I32),DP)
      Dv(18) = DuFine(IedgesAtElementFine(4,Ielf(4))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(4)),-3),1_I32),DP)
      Dv(19) = DuFine(IedgesAtElementFine(1,Ielf(4))+NMTfine)&
             * real(1-2*iand(      ItwistFine(Ielf(4))    ,1_I32),DP)
      Dv(20) = DuFine(IedgesAtElementFine(4,Ielf(1))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(1)),-3),1_I32),DP)
      Dv(21) = DuFine(IedgesAtElementFine(2,Ielf(1))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(1)),-1),1_I32),DP)
      Dv(22) = DuFine(IedgesAtElementFine(2,Ielf(2))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(2)),-1),1_I32),DP)
      Dv(23) = DuFine(IedgesAtElementFine(2,Ielf(3))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(3)),-1),1_I32),DP)
      Dv(24) = DuFine(IedgesAtElementFine(2,Ielf(4))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(4)),-1),1_I32),DP)
      Dv(25) = DuFine(2*NMTfine + Ielf(1))
      Dv(26) = DuFine(2*NMTfine + Ielf(2))
      Dv(27) = DuFine(2*NMTfine + Ielf(3))
      Dv(28) = DuFine(2*NMTfine + Ielf(4))

      DuCoarse(Idf(1)) = DuCoarse(Idf(1)) + Dtn(1)*(Q11*(Dv(5)+Dv(6)) &
                       +Q12*(Dv(1)+Dv(2)-Dv(10)-Dv(12))+Q14*(Dv(25)+Dv(26))&
                       -Q13*(Dv(27)+Dv(28)))
      DuCoarse(Idf(2)) = DuCoarse(Idf(2)) + Dtn(2)*(Q11*(Dv(7)+Dv(8)) &
                       +Q12*(Dv(3)+Dv(4)-Dv(9)-Dv(11))+Q14*(Dv(26)+Dv(27))&
                       -Q13*(Dv(28)+Dv(25)))
      DuCoarse(Idf(3)) = DuCoarse(Idf(3)) + Dtn(3)*(Q11*(Dv(1)+Dv(2)) &
                       +Q12*(Dv(5)+Dv(6)-Dv(10)-Dv(12))+Q14*(Dv(27)+Dv(28))&
                       -Q13*(Dv(25)+Dv(26)))
      DuCoarse(Idf(4)) = DuCoarse(Idf(4)) + Dtn(4)*(Q11*(Dv(3)+Dv(4)) &
                       +Q12*(Dv(7)+Dv(8)-Dv(9)-Dv(11))+Q14*(Dv(28)+Dv(25))&
                       -Q13*(Dv(26)+Dv(27)))

      DuCoarse(Idf(5)) = DuCoarse(Idf(5)) + Dtn(1)*Dtf(1)*(Q21*(Dv(2)-Dv(1)) &
                       +Q22*(Dv(12)-Dv(10))+Q23*(Dv(8)-Dv(3))+Q24*(Dv(6)-Dv(5)) &
                       +Q25*(Dv(4)-Dv(7))+Q26*(Dv(13)+Dv(14))+Q27*(Dv(22)-Dv(24)) &
                       -Q28*(Dv(20)+Dv(15))-Q29*(Dv(17)+Dv(18))+Q30*(Dv(16)+Dv(19)) &
                       +Q31*(Dv(26)-Dv(25))+Q32*(Dv(28)-Dv(27)))
      DuCoarse(Idf(6)) = DuCoarse(Idf(6)) + Dtn(2)*Dtf(2)*(Q21*(Dv(4)-Dv(3)) &
                       +Q22*(Dv(9)-Dv(11))+Q23*(Dv(2)-Dv(5))+Q24*(Dv(8)-Dv(7)) &
                       +Q25*(Dv(6)-Dv(1))+Q26*(Dv(15)+Dv(16))+Q27*(Dv(23)-Dv(21)) &
                       -Q28*(Dv(14)+Dv(17))-Q29*(Dv(19)+Dv(20))+Q30*(Dv(13)+Dv(18)) &
                       +Q31*(Dv(27)-Dv(26))+Q32*(Dv(25)-Dv(28)))
      DuCoarse(Idf(7)) = DuCoarse(Idf(7)) + Dtn(3)*Dtf(3)*(Q21*(Dv(6)-Dv(5)) &
                       +Q22*(Dv(10)-Dv(12))+Q23*(Dv(4)-Dv(7))+Q24*(Dv(2)-Dv(1)) &
                       +Q25*(Dv(8)-Dv(3))+Q26*(Dv(17)+Dv(18))-Q27*(Dv(22)-Dv(24)) &
                       -Q28*(Dv(16)+Dv(19))-Q29*(Dv(13)+Dv(14))+Q30*(Dv(15)+Dv(20)) &
                       +Q31*(Dv(28)-Dv(27))+Q32*(Dv(26)-Dv(25)))
      DuCoarse(Idf(8)) = DuCoarse(Idf(8)) + Dtn(4)*Dtf(4)*(Q21*(Dv(8)-Dv(7)) &
                       +Q22*(Dv(11)-Dv(9))+Q23*(Dv(6)-Dv(1))+Q24*(Dv(4)-Dv(3)) &
                       +Q25*(Dv(2)-Dv(5))+Q26*(Dv(19)+Dv(20))+Q27*(Dv(21)-Dv(23)) &
                       -Q28*(Dv(13)+Dv(18))-Q29*(Dv(15)+Dv(16))+Q30*(Dv(14)+Dv(17)) &
                       +Q31*(Dv(25)-Dv(28))+Q32*(Dv(27)-Dv(26)))

      DuCoarse(Idf(9)) = 0.25_DP*(Dv(25)+Dv(26)+Dv(27)+Dv(28))

    end do

    ! That is it

  end subroutine

  ! ***************************************************************************
  ! Support for Q2~ element with bubble
  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_prolUniformEB50_double (DuCoarse,DuFine,&
               IedgesAtElementCoarse,IedgesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               ItwistCoarse,ItwistFine,NMTcoarse,NMTfine,NELcoarse,NELfine)

!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! EB50, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Coarse grid vector
  real(DP), dimension(:), intent(in) :: DuCoarse

  ! IedgesAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementCoarse

  ! IedgesAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! ItwistIndexEdges array on the coarse grid
  integer(I32), dimension(:), intent(in) :: ItwistCoarse

  ! ItwistIndexEdges array on the fine grid
  integer(I32), dimension(:), intent(in) :: ItwistFine

  ! Number of egdes in the coarse grid
  integer, intent(in) :: NMTcoarse

  ! Number of egdes in the fine grid
  integer, intent(in) :: NMTfine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse

  ! Number of elements in the fine grid
  integer, intent(in) :: NELfine
!</input>

!<output>
  ! Fine grid vector
  real(DP), dimension(:), intent(out) :: DuFine
!</output>

!</subroutine>

  ! local variables
  integer :: iel,i
  integer, dimension(4) :: Ielf, Imc

  ! Local vectors
  real(DP), dimension(10) :: Dv
  real(DP), dimension(4:6,4) :: Dtf
  real(DP), dimension(4) :: Dtn
  integer, dimension(8,4) :: Idf

  ! "Projection" coefficients
  real(DP), parameter :: R11 =   27.0_DP / 808.0_DP
  real(DP), parameter :: R12 = 1239.0_DP / 808.0_DP

  real(DP), parameter :: R21 =   81.0_DP / 808.0_DP
  real(DP), parameter :: R22 =  222.0_DP / 808.0_DP
  real(DP), parameter :: R23 =  323.0_DP / 808.0_DP
  real(DP), parameter :: R24 =  384.0_DP / 808.0_DP

  real(DP), parameter :: R31 = 0.0_DP
  real(DP), parameter :: R32 = 0.0_DP
  real(DP), parameter :: R33 = 0.0_DP
  real(DP), parameter :: R34 = 1 / 16.0_DP

!  ! "Interpolation" coefficients
!  real(DP), parameter :: R11 =  3.0_DP / 32.0_DP  ! 0.09375_DP
!  real(DP), parameter :: R12 = 51.0_DP / 32.0_DP  ! 1.59375_DP
!
!  real(DP), parameter :: R21 =  3.0_DP / 32.0_DP  ! 0.09375_DP
!  real(DP), parameter :: R22 =  9.0_DP / 32.0_DP  ! 0.28125_DP
!  real(DP), parameter :: R23 = 13.0_DP / 32.0_DP  ! 0.40625_DP
!  real(DP), parameter :: R24 = 15.0_DP / 32.0_DP  ! 0.46875_DP
!
!  real(DP), parameter :: R31 =  1.0_DP / 12.0_DP  ! 0.0833333...
!  real(DP), parameter :: R32 =  1.0_DP / 24.0_DP  ! 0.0416666...
!  real(DP), parameter :: R33 =  1.0_DP / 48.0_DP  ! 0.0208333...
!  real(DP), parameter :: R34 = 25.0_DP / 16.0_DP  ! 1.5625_DP

    ! Clear the output vector
    call lalg_clearVector(DuFine)

    ! Loop over the coarse grid elements
    do iel = 1, NELcoarse

      ! Get the element numbers of the fine grid elements
      Ielf(1) = iel
      Ielf(2) = IneighboursAtElementFine(2,Ielf(1))
      Ielf(3) = IneighboursAtElementFine(2,Ielf(2))
      Ielf(4) = IneighboursAtElementFine(2,Ielf(3))

      ! Get the edge indices of the coarse grid element
      Imc(1) = IedgesAtElementCoarse(1,iel)
      Imc(2) = IedgesAtElementCoarse(2,iel)
      Imc(3) = IedgesAtElementCoarse(3,iel)
      Imc(4) = IedgesAtElementCoarse(4,iel)

      Dtn = 1.0_DP
      do i = 1, 4

        ! Calculate the Edge-Int-Mean DOFs on the fine mesh
        Idf(1,i) = IedgesAtElementFine(1,Ielf(i))
        Idf(2,i) = IedgesAtElementFine(2,Ielf(i))
        Idf(3,i) = IedgesAtElementFine(4,Ielf(i))
        ! Calculate the Edge-Legendre-Int-Mean DOFs on the fine mesh
        Idf(4,i) = Idf(1,i) + NMTfine
        Idf(5,i) = Idf(2,i) + NMTfine
        Idf(6,i) = Idf(3,i) + NMTfine
        ! Calculate the Quad-Int-Mean DOF on the fine mesh
        Idf(7,i) = 2*NMTfine + Ielf(i)
        ! Calculate the Quad-Legendre-Int-Mean DOF on the fine mesh
        Idf(8,i) = 2*NMTfine + NELfine + Ielf(i)

        ! Calculate the twist factors for the Edge-Legendre-Int-Mean DOFs
        ! on the fine mesh
        Dtf(4,i) = real(1-2*iand(      ItwistFine(Ielf(i))    ,1_I32),DP)
        Dtf(5,i) = real(1-2*iand(ishft(ItwistFine(Ielf(i)),-1),1_I32),DP)
        Dtf(6,i) = real(1-2*iand(ishft(ItwistFine(Ielf(i)),-3),1_I32),DP)

        ! If there is a neighbour at the coarse mesh edge, then we
        ! set the corresponding factor to 1/2
        if (IneighboursAtElementCoarse(i,iel) .ne. 0) Dtn(i) = 0.5_DP

      end do

      ! Get the values of the corresponding coarse mesh DOFs
      Dv(1:4) = DuCoarse(Imc(1:4))
      Dv(  5) = DuCoarse(Imc(1) + NMTcoarse) &
              * real(1-2*iand(      ItwistCoarse(iel)    ,1_I32),DP)
      Dv(  6) = DuCoarse(Imc(2) + NMTcoarse) &
              * real(1-2*iand(ishft(ItwistCoarse(iel),-1),1_I32),DP)
      Dv(  7) = DuCoarse(Imc(3) + NMTcoarse) &
              * real(1-2*iand(ishft(ItwistCoarse(iel),-2),1_I32),DP)
      Dv(  8) = DuCoarse(Imc(4) + NMTcoarse) &
              * real(1-2*iand(ishft(ItwistCoarse(iel),-3),1_I32),DP)
      Dv(  9) = DuCoarse(2*NMTcoarse + iel)
      Dv( 10) = DuCoarse(2*NMTcoarse + NELcoarse + iel)

      DuFine(Idf(1,1)) = DuFine(Idf(1,1)) + (Dv(1)-R12*Dv(5)+&
                         R11*(-Dv(6)-Dv(7)-Dv(8)))*Dtn(1)
      DuFine(Idf(3,2)) = DuFine(Idf(3,2)) + (Dv(1)+R12*Dv(5)+&
                         R11*(Dv(6)+Dv(7)+Dv(8)))*Dtn(1)
      DuFine(Idf(1,2)) = DuFine(Idf(1,2)) + (Dv(2)-R12*Dv(6)+&
                         R11*(-Dv(5)-Dv(7)-Dv(8)))*Dtn(2)
      DuFine(Idf(3,3)) = DuFine(Idf(3,3)) + (Dv(2)+R12*Dv(6)+&
                         R11*(Dv(5)+Dv(7)+Dv(8)))*Dtn(2)
      DuFine(Idf(1,3)) = DuFine(Idf(1,3)) + (Dv(3)-R12*Dv(7)+&
                         R11*(-Dv(5)-Dv(6)-Dv(8)))*Dtn(3)
      DuFine(Idf(3,4)) = DuFine(Idf(3,4)) + (Dv(3)+R12*Dv(7)+&
                         R11*(Dv(5)+Dv(6)+Dv(8)))*Dtn(3)
      DuFine(Idf(1,4)) = DuFine(Idf(1,4)) + (Dv(4)-R12*Dv(8)+&
                         R11*(-Dv(5)-Dv(6)-Dv(7)))*Dtn(4)
      DuFine(Idf(3,1)) = DuFine(Idf(3,1)) + (Dv(4)+R12*Dv(8)+&
                         R11*(Dv(5)+Dv(6)+Dv(7)))*Dtn(4)
      DuFine(Idf(2,1)) = DuFine(Idf(2,1)) - 0.25_DP*(Dv(2)+Dv(4))+&
                         0.375_DP*(Dv(1)-Dv(3)+Dv(6)-Dv(8))+1.5_DP*Dv(9)
      DuFine(Idf(2,2)) = DuFine(Idf(2,2)) - 0.25_DP*(Dv(1)+Dv(3))+&
                         0.375_DP*(Dv(2)-Dv(4)-Dv(5)+Dv(7))+1.5_DP*Dv(9)
      DuFine(Idf(2,3)) = DuFine(Idf(2,3)) - 0.25_DP*(Dv(2)+Dv(4))+&
                         0.375_DP*(-Dv(1)+Dv(3)-Dv(6)+Dv(8))+1.5_DP*Dv(9)
      DuFine(Idf(2,4)) = DuFine(Idf(2,4)) - 0.25_DP*(Dv(1)+Dv(3))+&
                         0.375_DP*(-Dv(2)+Dv(4)+Dv(5)-Dv(7))+1.5_DP*Dv(9)

      DuFine(Idf(4,1)) = DuFine(Idf(4,1)) + Dtn(1)*Dtf(4,1)*(&
           0.125_DP*(Dv(1)-Dv(2)-Dv(3)-Dv(4))+R23*Dv(5)+R22*Dv(6)&
           -R21*Dv(7)-R24*Dv(8)+0.25_DP*Dv(9)-6.25_DP*Dv(10))
      DuFine(Idf(6,2)) = DuFine(Idf(6,2)) + Dtn(1)*Dtf(6,2)*(&
           0.125_DP*(-Dv(1)+Dv(2)+Dv(3)+Dv(4))+R23*Dv(5)-R24*Dv(6)&
           -R21*Dv(7)+R22*Dv(8)-0.25_DP*Dv(9)+6.25_DP*Dv(10))
      DuFine(Idf(4,2)) = DuFine(Idf(4,2)) + Dtn(2)*Dtf(4,2)*(&
           0.125_DP*(-Dv(1)+Dv(2)-Dv(3)-Dv(4))-R24*Dv(5)+R23*Dv(6)&
           +R22*Dv(7)-R21*Dv(8)+0.25_DP*Dv(9)-6.25_DP*Dv(10))
      DuFine(Idf(6,3)) = DuFine(Idf(6,3)) + Dtn(2)*Dtf(6,3)*(&
           0.125_DP*(Dv(1)-Dv(2)+Dv(3)+Dv(4))+R22*Dv(5)+R23*Dv(6)&
           -R24*Dv(7)-R21*Dv(8)-0.25_DP*Dv(9)+6.25_DP*Dv(10))
      DuFine(Idf(4,3)) = DuFine(Idf(4,3)) + Dtn(3)*Dtf(4,3)*(&
           0.125_DP*(-Dv(1)-Dv(2)+Dv(3)-Dv(4))-R21*Dv(5)-R24*Dv(6)&
           +R23*Dv(7)+R22*Dv(8)+0.25_DP*Dv(9)-6.25_DP*Dv(10))
      DuFine(Idf(6,4)) = DuFine(Idf(6,4)) + Dtn(3)*Dtf(6,4)*(&
           0.125_DP*(Dv(1)+Dv(2)-Dv(3)+Dv(4))-R21*Dv(5)+R22*Dv(6)&
           +R23*Dv(7)-R24*Dv(8)-0.25_DP*Dv(9)+6.25_DP*Dv(10))
      DuFine(Idf(4,4)) = DuFine(Idf(4,4)) + Dtn(4)*Dtf(4,4)*(&
           0.125_DP*(-Dv(1)-Dv(2)-Dv(3)+Dv(4))+R22*Dv(5)-R21*Dv(6)&
           -R24*Dv(7)+R23*Dv(8)+0.25_DP*Dv(9)-6.25_DP*Dv(10))
      DuFine(Idf(6,1)) = DuFine(Idf(6,1)) + Dtn(4)*Dtf(6,1)*(&
           0.125_DP*(Dv(1)+Dv(2)+Dv(3)-Dv(4))-R24*Dv(5)-R21*Dv(6)&
           +R22*Dv(7)+R23*Dv(8)-0.25_DP*Dv(9)+6.25_DP*Dv(10))

      DuFine(Idf(5,1)) = DuFine(Idf(5,1)) + Dtf(5,1)*(0.25_DP*(-Dv(1)+Dv(9))&
                         +0.125_DP*(-Dv(6)+Dv(8))+3.125_DP*Dv(10))
      DuFine(Idf(5,2)) = DuFine(Idf(5,2)) + Dtf(5,2)*(0.25_DP*(-Dv(2)+Dv(9))&
                         +0.125_DP*(Dv(5)-Dv(7))+3.125_DP*Dv(10))
      DuFine(Idf(5,3)) = DuFine(Idf(5,3)) + Dtf(5,3)*(0.25_DP*(-Dv(3)+Dv(9))&
                         +0.125_DP*(Dv(6)-Dv(8))+3.125_DP*Dv(10))
      DuFine(Idf(5,4)) = DuFine(Idf(5,4)) + Dtf(5,4)*(0.25_DP*(-Dv(4)+Dv(9))&
                         +0.125_DP*(-Dv(5)+Dv(7))+3.125_DP*Dv(10))

      DuFine(Idf(7,1)) = DuFine(Idf(7,1)) + 0.25_DP*(Dv(1)-Dv(2)-Dv(3)+Dv(4))+&
                         0.1875_DP*(-Dv(5)+Dv(6)-Dv(7)+Dv(8))+Dv(9)
      DuFine(Idf(7,2)) = DuFine(Idf(7,2)) + 0.25_DP*(Dv(1)+Dv(2)-Dv(3)-Dv(4))+&
                         0.1875_DP*(Dv(5)-Dv(6)+Dv(7)-Dv(8))+Dv(9)
      DuFine(Idf(7,3)) = DuFine(Idf(7,3)) + 0.25_DP*(-Dv(1)+Dv(2)+Dv(3)-Dv(4))+&
                         0.1875_DP*(-Dv(5)+Dv(6)-Dv(7)+Dv(8))+Dv(9)
      DuFine(Idf(7,4)) = DuFine(Idf(7,4)) + 0.25_DP*(-Dv(1)-Dv(2)+Dv(3)+Dv(4))+&
                         0.1875_DP*(Dv(5)-Dv(6)+Dv(7)-Dv(8))+Dv(9)

      DuFine(Idf(8,1)) = DuFine(Idf(8,1)) + R33*(-Dv(1)+Dv(2)+Dv(3)-Dv(4))+&
                         R31*(-Dv(5)+Dv(8))+R32*(-Dv(6)+Dv(7))+R34*Dv(10)
      DuFine(Idf(8,2)) = DuFine(Idf(8,2)) + R33*(-Dv(1)-Dv(2)+Dv(3)+Dv(4))+&
                         R31*(Dv(5)-Dv(6))+R32*(-Dv(7)+Dv(8))+R34*Dv(10)
      DuFine(Idf(8,3)) = DuFine(Idf(8,3)) + R33*(Dv(1)-Dv(2)-Dv(3)+Dv(4))+&
                         R32*(Dv(5)-Dv(8))+R31*(Dv(6)-Dv(7))+R34*Dv(10)
      DuFine(Idf(8,4)) = DuFine(Idf(8,4)) + R33*(Dv(1)+Dv(2)-Dv(3)-Dv(4))+&
                         R32*(-Dv(5)+Dv(6))+R31*(Dv(7)-Dv(8))+R34*Dv(10)

      ! Go for the next element

    end do

    ! That is it

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_restUniformEB50_double (DuCoarse,DuFine, &
               IedgesAtElementCoarse,IedgesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               ItwistCoarse,ItwistFine,NMTcoarse,NMTfine,NELcoarse,NELfine)

!<description>
  ! Restricts a defect vector from a fine grid to a coarse grid.
  ! EB50, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IedgesAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementCoarse

  ! IedgesAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! ItwistIndexEdges array on the coarse grid
  integer(I32), dimension(:), intent(in) :: ItwistCoarse

  ! ItwistIndexEdges array on the fine grid
  integer(I32), dimension(:), intent(in) :: ItwistFine

  ! Number of egdes in the coarse grid
  integer, intent(in) :: NMTcoarse

  ! Number of egdes in the fine grid
  integer, intent(in) :: NMTfine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse

  ! Number of elements in the fine grid
  integer, intent(in) :: NELfine
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  integer :: iel,i
  integer, dimension(4) :: Ielf

  ! local vectors
  real(DP), dimension(32) :: dv
  real(DP), dimension(4) :: Dtn,Dtf
  integer, dimension(10) :: Idf

  ! "Projection" coefficients
  real(DP), parameter :: R11 =   27.0_DP / 808.0_DP
  real(DP), parameter :: R12 = 1239.0_DP / 808.0_DP

  real(DP), parameter :: R21 =   81.0_DP / 808.0_DP
  real(DP), parameter :: R22 =  222.0_DP / 808.0_DP
  real(DP), parameter :: R23 =  323.0_DP / 808.0_DP
  real(DP), parameter :: R24 =  384.0_DP / 808.0_DP

  real(DP), parameter :: R31 = 0.0_DP
  real(DP), parameter :: R32 = 0.0_DP
  real(DP), parameter :: R33 = 0.0_DP
  real(DP), parameter :: R34 = 1 / 16.0_DP

!  ! "Interpolation" coefficients
!  real(DP), parameter :: R11 =  3.0_DP / 32.0_DP  ! 0.09375_DP
!  real(DP), parameter :: R12 = 51.0_DP / 32.0_DP  ! 1.59375_DP
!
!  real(DP), parameter :: R21 =  3.0_DP / 32.0_DP  ! 0.09375_DP
!  real(DP), parameter :: R22 =  9.0_DP / 32.0_DP  ! 0.28125_DP
!  real(DP), parameter :: R23 = 13.0_DP / 32.0_DP  ! 0.40625_DP
!  real(DP), parameter :: R24 = 15.0_DP / 32.0_DP  ! 0.46875_DP
!
!  real(DP), parameter :: R31 =  1.0_DP / 12.0_DP  ! 0.0833333...
!  real(DP), parameter :: R32 =  1.0_DP / 24.0_DP  ! 0.0416666...
!  real(DP), parameter :: R33 =  1.0_DP / 48.0_DP  ! 0.0208333...
!  real(DP), parameter :: R34 = 25.0_DP / 16.0_DP  ! 1.5625_DP

    ! Clear the output vector
    call lalg_clearVector(DuCoarse)

    ! Loop over the coarse grid elements
    do iel = 1, NELcoarse

      ! Get the element numbers of the fine grid elements
      Ielf(1) = iel
      Ielf(2) = IneighboursAtElementFine(2,Ielf(1))
      Ielf(3) = IneighboursAtElementFine(2,Ielf(2))
      Ielf(4) = IneighboursAtElementFine(2,Ielf(3))

      ! Calculate the coarse grid DOFs
      Idf(1:4) = IedgesAtElementCoarse(1:4,iel)
      Idf(5:8) = Idf(1:4) + NMTcoarse
      Idf(  9) = 2*NMTcoarse + iel
      Idf( 10) = 2*NMTcoarse + NELcoarse + iel

      ! Calculate twist index factors and neighbour scales for the
      ! coarse grid edges
      Dtn = 1.0_DP
      do i = 1, 4

        ! Twist index factor
        Dtf(i) = real(1-2*iand(ishft(ItwistCoarse(iel),1-i),1_I32),DP)

        ! Neighbour factor
        if (IneighboursAtElementCoarse(i,iel) .ne. 0) Dtn(i) = 0.5_DP

      end do

      ! Get the values of the corresponding fine mesh DOFs
      Dv( 1) = DuFine(IedgesAtElementFine(1,Ielf(1)))*Dtn(1)
      Dv( 2) = DuFine(IedgesAtElementFine(4,Ielf(2)))*Dtn(1)
      Dv( 3) = DuFine(IedgesAtElementFine(1,Ielf(2)))*Dtn(2)
      Dv( 4) = DuFine(IedgesAtElementFine(4,Ielf(3)))*Dtn(2)
      Dv( 5) = DuFine(IedgesAtElementFine(1,Ielf(3)))*Dtn(3)
      Dv( 6) = DuFine(IedgesAtElementFine(4,Ielf(4)))*Dtn(3)
      Dv( 7) = DuFine(IedgesAtElementFine(1,Ielf(4)))*Dtn(4)
      Dv( 8) = DuFine(IedgesAtElementFine(4,Ielf(1)))*Dtn(4)
      Dv( 9) = DuFine(IedgesAtElementFine(2,Ielf(1)))
      Dv(10) = DuFine(IedgesAtElementFine(2,Ielf(2)))
      Dv(11) = DuFine(IedgesAtElementFine(2,Ielf(3)))
      Dv(12) = DuFine(IedgesAtElementFine(2,Ielf(4)))
      Dv(13) = DuFine(IedgesAtElementFine(1,Ielf(1))+NMTfine)&
             * real(1-2*iand(      ItwistFine(Ielf(1))    ,1_I32),DP)*Dtn(1)
      Dv(14) = DuFine(IedgesAtElementFine(4,Ielf(2))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(2)),-3),1_I32),DP)*Dtn(1)
      Dv(15) = DuFine(IedgesAtElementFine(1,Ielf(2))+NMTfine)&
             * real(1-2*iand(      ItwistFine(Ielf(2))    ,1_I32),DP)*Dtn(2)
      Dv(16) = DuFine(IedgesAtElementFine(4,Ielf(3))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(3)),-3),1_I32),DP)*Dtn(2)
      Dv(17) = DuFine(IedgesAtElementFine(1,Ielf(3))+NMTfine)&
             * real(1-2*iand(      ItwistFine(Ielf(3))    ,1_I32),DP)*Dtn(3)
      Dv(18) = DuFine(IedgesAtElementFine(4,Ielf(4))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(4)),-3),1_I32),DP)*Dtn(3)
      Dv(19) = DuFine(IedgesAtElementFine(1,Ielf(4))+NMTfine)&
             * real(1-2*iand(      ItwistFine(Ielf(4))    ,1_I32),DP)*Dtn(4)
      Dv(20) = DuFine(IedgesAtElementFine(4,Ielf(1))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(1)),-3),1_I32),DP)*Dtn(4)
      Dv(21) = DuFine(IedgesAtElementFine(2,Ielf(1))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(1)),-1),1_I32),DP)
      Dv(22) = DuFine(IedgesAtElementFine(2,Ielf(2))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(2)),-1),1_I32),DP)
      Dv(23) = DuFine(IedgesAtElementFine(2,Ielf(3))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(3)),-1),1_I32),DP)
      Dv(24) = DuFine(IedgesAtElementFine(2,Ielf(4))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(4)),-1),1_I32),DP)
      Dv(25) = DuFine(2*NMTfine + Ielf(1))
      Dv(26) = DuFine(2*NMTfine + Ielf(2))
      Dv(27) = DuFine(2*NMTfine + Ielf(3))
      Dv(28) = DuFine(2*NMTfine + Ielf(4))
      Dv(29) = DuFine(2*NMTfine + NELfine + Ielf(1))
      Dv(30) = DuFine(2*NMTfine + NELfine + Ielf(2))
      Dv(31) = DuFine(2*NMTfine + NELfine + Ielf(3))
      Dv(32) = DuFine(2*NMTfine + NELfine + Ielf(4))

      DuCoarse(Idf(1)) = DuCoarse(Idf(1)) + 0.125_DP*(+Dv(13)-Dv(14)-Dv(15)+Dv(16)&
                        -Dv(17)+Dv(18)-Dv(19)+Dv(20))+0.25_DP*(-Dv(21)+Dv(25)&
                        +Dv(26)-Dv(27)-Dv(28)-Dv(10)-Dv(12))+R33*(-Dv(29)-Dv(30)&
                        +Dv(31)+Dv(32))+Dv(1)+Dv(2)+0.375_DP*(Dv(9)-Dv(11))
      DuCoarse(Idf(2)) = DuCoarse(Idf(2)) + 0.125_DP*(-Dv(13)+Dv(14)+Dv(15)-Dv(16)&
                        -Dv(17)+Dv(18)-Dv(19)+Dv(20))+0.25_DP*(-Dv(22)-Dv(25)&
                        +Dv(26)+Dv(27)-Dv(28)-Dv(9)-Dv(11))+R33*(Dv(29)-Dv(30)&
                        -Dv(31)+Dv(32))+Dv(3)+Dv(4)+0.375_DP*(Dv(10)-Dv(12))
      DuCoarse(Idf(3)) = DuCoarse(Idf(3)) + 0.125_DP*(-Dv(13)+Dv(14)-Dv(15)+Dv(16)&
                        +Dv(17)-Dv(18)-Dv(19)+Dv(20))+0.25_DP*(-Dv(23)-Dv(25)&
                        -Dv(26)+Dv(27)+Dv(28)-Dv(10)-Dv(12))+R33*(Dv(29)+Dv(30)&
                        -Dv(31)-Dv(32))+Dv(5)+Dv(6)+0.375_DP*(-Dv(9)+Dv(11))
      DuCoarse(Idf(4)) = DuCoarse(Idf(4)) + 0.125_DP*(-Dv(13)+Dv(14)-Dv(15)+Dv(16)&
                        -Dv(17)+Dv(18)+Dv(19)-Dv(20))+0.25_DP*(-Dv(24)+Dv(25)&
                        -Dv(26)-Dv(27)+Dv(28)-Dv(9)-Dv(11))+R33*(-Dv(29)+Dv(30)&
                        +Dv(31)-Dv(32))+Dv(7)+Dv(8)+0.375_DP*(-Dv(10)+Dv(12))

      DuCoarse(Idf(5)) = DuCoarse(Idf(5)) + Dtf(1)*(R24*(-Dv(15)-Dv(20))&
                        +R22*(Dv(16)+Dv(19))+0.125_DP*(Dv(22)-Dv(24))&
                        +0.1875_DP*(-Dv(25)+Dv(26)-Dv(27)+Dv(28))+R31*(-Dv(29)+Dv(30))&
                        +R32*(Dv(31)-Dv(32))+R12*(-Dv(1)+Dv(2))&
                        +R11*(-Dv(3)+Dv(4)-Dv(5)+Dv(6)-Dv(7)+Dv(8))&
                        +R21*(-Dv(17)-Dv(18))&
                        +0.375_DP*(-Dv(10)+Dv(12))+R23*(Dv(13)+Dv(14)))
      DuCoarse(Idf(6)) = DuCoarse(Idf(6)) + Dtf(2)*(R24*(-Dv(14)-Dv(17))&
                        +R22*(Dv(13)+Dv(18))+0.125_DP*(-Dv(21)+Dv(23))&
                        +0.1875_DP*(Dv(25)-Dv(26)+Dv(27)-Dv(28))+R31*(-Dv(30)+Dv(31))&
                        +R32*(-Dv(29)+Dv(32))+R12*(-Dv(3)+Dv(4))&
                        +R11*(-Dv(1)+Dv(2)-Dv(5)+Dv(6)-Dv(7)+Dv(8))&
                        +R21*(-Dv(19)-Dv(20))&
                        +0.375_DP*(Dv(9)-Dv(11))+R23*(Dv(15)+Dv(16)))
      DuCoarse(Idf(7)) = DuCoarse(Idf(7)) + Dtf(3)*(R24*(-Dv(16)-Dv(19))&
                        +R22*(Dv(15)+Dv(20))+0.125_DP*(Dv(24)-Dv(22))&
                        +0.1875_DP*(-Dv(25)+Dv(26)-Dv(27)+Dv(28))+R31*(-Dv(31)+Dv(32))&
                        +R32*(Dv(29)-Dv(30))+R12*(-Dv(5)+Dv(6))&
                        +R11*(-Dv(1)+Dv(2)-Dv(3)+Dv(4)-Dv(7)+Dv(8))&
                        +R21*(-Dv(13)-Dv(14))&
                        +0.375_DP*(Dv(10)-Dv(12))+R23*(Dv(17)+Dv(18)))
      DuCoarse(Idf(8)) = DuCoarse(Idf(8)) + Dtf(4)*(R24*(-Dv(13)-Dv(18))&
                        +R22*(Dv(14)+Dv(17))+0.125_DP*(Dv(21)-Dv(23))&
                        +0.1875_DP*(Dv(25)-Dv(26)+Dv(27)-Dv(28))+R31*(Dv(29)-Dv(32))&
                        +R32*(Dv(30)-Dv(31))+R12*(-Dv(7)+Dv(8))&
                        +R11*(-Dv(1)+Dv(2)-Dv(3)+Dv(4)-Dv(5)+Dv(6))&
                        +R21*(-Dv(15)-Dv(16))&
                        +0.375_DP*(-Dv(9)+Dv(11))+R23*(Dv(19)+Dv(20)))

      DuCoarse(Idf(9)) = DuCoarse(Idf(9)) + 0.25_DP*(Dv(13)-Dv(14)+Dv(15)-Dv(16)&
                        +Dv(17)-Dv(18)+Dv(19)-Dv(20)+Dv(21)+Dv(22)+Dv(23)+Dv(24))&
                        +Dv(25)+Dv(26)+Dv(27)+Dv(28)+1.5_DP*(Dv(9)+Dv(10)+Dv(11)+Dv(12))

      DuCoarse(Idf(10)) = DuCoarse(Idf(10)) + 6.25_DP*(-Dv(13)+Dv(14)-Dv(15)&
                        +Dv(16)-Dv(17)+Dv(18)-Dv(19)+Dv(20))+3.125_DP*(Dv(21)+Dv(22)&
                        +Dv(23)+Dv(24))+R34*(Dv(29)+Dv(30)+Dv(31)+Dv(32))

    end do

    ! That is it

  end subroutine


  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_interpUniformEB50_double (DuCoarse,DuFine, &
               IedgesAtElementCoarse,IedgesAtElementFine,&
               IneighboursAtElementCoarse,IneighboursAtElementFine,&
               ItwistCoarse,ItwistFine,NMTcoarse,NMTfine,NELcoarse,NELfine)

!<description>
  ! Restricts a solution vector from a fine grid to a coarse grid.
  ! EB50, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IedgesAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementCoarse

  ! IedgesAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IedgesAtElementFine

  ! IneighboursAtElement array on the coarse grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! ItwistIndexEdges array on the coarse grid
  integer(I32), dimension(:), intent(in) :: ItwistCoarse

  ! ItwistIndexEdges array on the fine grid
  integer(I32), dimension(:), intent(in) :: ItwistFine

  ! Number of egdes in the coarse grid
  integer, intent(in) :: NMTcoarse

  ! Number of egdes in the fine grid
  integer, intent(in) :: NMTfine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse

  ! Number of elements in the fine grid
  integer, intent(in) :: NELfine
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  integer :: iel,i
  integer, dimension(4) :: Ielf

  ! local vectors
  real(DP), dimension(32) :: dv
  real(DP), dimension(4) :: Dtn,Dtf
  integer, dimension(10) :: Idf

  ! "Projection" coefficients
  real(DP), parameter :: Q11 =  1.0_DP /  8.0_DP    ! 0.125
  real(DP), parameter :: Q12 =  1.0_DP /  4.0_DP    ! 0.25
  real(DP), parameter :: Q13 =  3.0_DP / 16.0_DP    ! 0.1875
  real(DP), parameter :: Q14 =  9.0_DP / 16.0_DP    ! 0.5625

  real(DP), parameter :: Q21 = 19.0_DP /  96.0_DP   ! 0.197916666...
  real(DP), parameter :: Q22 =  1.0_DP /   8.0_DP   ! 0.125
  real(DP), parameter :: Q23 =  5.0_DP / 192.0_DP   ! 0.026041666...
  real(DP), parameter :: Q24 =  1.0_DP /  96.0_DP   ! 0.010416666...
  real(DP), parameter :: Q25 = 11.0_DP / 192.0_DP   ! 0.057291666...
  real(DP), parameter :: Q26 = 13.0_DP / 128.0_DP   ! 0.1015625
  real(DP), parameter :: Q27 =  3.0_DP /  64.0_DP   ! 0.046875
  real(DP), parameter :: Q28 =  3.0_DP / 128.0_DP   ! 0.0234375
  real(DP), parameter :: Q29 =  1.0_DP / 128.0_DP   ! 0.0078125
  real(DP), parameter :: Q30 =  7.0_DP / 128.0_DP   ! 0.0546875
  real(DP), parameter :: Q31 =  9.0_DP /  32.0_DP
  real(DP), parameter :: Q32 =  3.0_DP /  32.0_DP

  real(DP), parameter :: Q41 =  1.0_DP /  640.0_DP   ! 0.0015625
  real(DP), parameter :: Q42 =  1.0_DP /  320.0_DP   ! 0.003125
  real(DP), parameter :: Q43 = 21.0_DP / 1280.0_DP   ! 0.01640625
  real(DP), parameter :: Q44 =  9.0_DP /  640.0_DP   ! 0.0140625
  real(DP), parameter :: Q45 =  1.0_DP /   64.0_DP   ! 0.015625


    ! Clear the output vector
    call lalg_clearVector(DuCoarse)

    ! Loop over the coarse grid elements
    do iel = 1, NELcoarse

      ! Get the element numbers of the fine grid elements
      Ielf(1) = iel
      Ielf(2) = IneighboursAtElementFine(2,Ielf(1))
      Ielf(3) = IneighboursAtElementFine(2,Ielf(2))
      Ielf(4) = IneighboursAtElementFine(2,Ielf(3))

      ! Calculate the coarse grid DOFs
      Idf(1:4) = IedgesAtElementCoarse(1:4,iel)
      Idf(5:8) = Idf(1:4) + NMTcoarse
      Idf(  9) = 2*NMTcoarse + iel
      Idf( 10) = 2*NMTcoarse + NELcoarse + iel

      ! Calculate twist index factors and neighbour scales for the
      ! coarse grid edges
      Dtn = 1.0_DP
      do i = 1, 4

        ! Twist index factor
        Dtf(i) = real(1-2*iand(ishft(ItwistCoarse(iel),1-i),1_I32),DP)

        ! Neighbour factor
        if (IneighboursAtElementCoarse(i,iel) .ne. 0) Dtn(i) = 0.5_DP

      end do

      ! Get the values of the corresponding fine mesh DOFs
      Dv( 1) = DuFine(IedgesAtElementFine(1,Ielf(1)))
      Dv( 2) = DuFine(IedgesAtElementFine(4,Ielf(2)))
      Dv( 3) = DuFine(IedgesAtElementFine(1,Ielf(2)))
      Dv( 4) = DuFine(IedgesAtElementFine(4,Ielf(3)))
      Dv( 5) = DuFine(IedgesAtElementFine(1,Ielf(3)))
      Dv( 6) = DuFine(IedgesAtElementFine(4,Ielf(4)))
      Dv( 7) = DuFine(IedgesAtElementFine(1,Ielf(4)))
      Dv( 8) = DuFine(IedgesAtElementFine(4,Ielf(1)))
      Dv( 9) = DuFine(IedgesAtElementFine(2,Ielf(1)))
      Dv(10) = DuFine(IedgesAtElementFine(2,Ielf(2)))
      Dv(11) = DuFine(IedgesAtElementFine(2,Ielf(3)))
      Dv(12) = DuFine(IedgesAtElementFine(2,Ielf(4)))
      Dv(13) = DuFine(IedgesAtElementFine(1,Ielf(1))+NMTfine)&
             * real(1-2*iand(      ItwistFine(Ielf(1))    ,1_I32),DP)
      Dv(14) = DuFine(IedgesAtElementFine(4,Ielf(2))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(2)),-3),1_I32),DP)
      Dv(15) = DuFine(IedgesAtElementFine(1,Ielf(2))+NMTfine)&
             * real(1-2*iand(      ItwistFine(Ielf(2))    ,1_I32),DP)
      Dv(16) = DuFine(IedgesAtElementFine(4,Ielf(3))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(3)),-3),1_I32),DP)
      Dv(17) = DuFine(IedgesAtElementFine(1,Ielf(3))+NMTfine)&
             * real(1-2*iand(      ItwistFine(Ielf(3))    ,1_I32),DP)
      Dv(18) = DuFine(IedgesAtElementFine(4,Ielf(4))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(4)),-3),1_I32),DP)
      Dv(19) = DuFine(IedgesAtElementFine(1,Ielf(4))+NMTfine)&
             * real(1-2*iand(      ItwistFine(Ielf(4))    ,1_I32),DP)
      Dv(20) = DuFine(IedgesAtElementFine(4,Ielf(1))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(1)),-3),1_I32),DP)
      Dv(21) = DuFine(IedgesAtElementFine(2,Ielf(1))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(1)),-1),1_I32),DP)
      Dv(22) = DuFine(IedgesAtElementFine(2,Ielf(2))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(2)),-1),1_I32),DP)
      Dv(23) = DuFine(IedgesAtElementFine(2,Ielf(3))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(3)),-1),1_I32),DP)
      Dv(24) = DuFine(IedgesAtElementFine(2,Ielf(4))+NMTfine)&
             * real(1-2*iand(ishft(ItwistFine(Ielf(4)),-1),1_I32),DP)
      Dv(25) = DuFine(2*NMTfine + Ielf(1))
      Dv(26) = DuFine(2*NMTfine + Ielf(2))
      Dv(27) = DuFine(2*NMTfine + Ielf(3))
      Dv(28) = DuFine(2*NMTfine + Ielf(4))
      Dv(29) = DuFine(2*NMTfine + NELfine + Ielf(1))
      Dv(30) = DuFine(2*NMTfine + NELfine + Ielf(2))
      Dv(31) = DuFine(2*NMTfine + NELfine + Ielf(3))
      Dv(32) = DuFine(2*NMTfine + NELfine + Ielf(4))

      DuCoarse(Idf(1)) = DuCoarse(Idf(1)) + Dtn(1)*(Q11*(Dv(5)+Dv(6)) &
                       +Q12*(Dv(1)+Dv(2)-Dv(10)-Dv(12))+Q14*(Dv(25)+Dv(26))&
                       -Q13*(Dv(27)+Dv(28)))
      DuCoarse(Idf(2)) = DuCoarse(Idf(2)) + Dtn(2)*(Q11*(Dv(7)+Dv(8)) &
                       +Q12*(Dv(3)+Dv(4)-Dv(9)-Dv(11))+Q14*(Dv(26)+Dv(27))&
                       -Q13*(Dv(28)+Dv(25)))
      DuCoarse(Idf(3)) = DuCoarse(Idf(3)) + Dtn(3)*(Q11*(Dv(1)+Dv(2)) &
                       +Q12*(Dv(5)+Dv(6)-Dv(10)-Dv(12))+Q14*(Dv(27)+Dv(28))&
                       -Q13*(Dv(25)+Dv(26)))
      DuCoarse(Idf(4)) = DuCoarse(Idf(4)) + Dtn(4)*(Q11*(Dv(3)+Dv(4)) &
                       +Q12*(Dv(7)+Dv(8)-Dv(9)-Dv(11))+Q14*(Dv(28)+Dv(25))&
                       -Q13*(Dv(26)+Dv(27)))

      DuCoarse(Idf(5)) = DuCoarse(Idf(5)) + Dtn(1)*Dtf(1)*(Q21*(Dv(2)-Dv(1)) &
                       +Q22*(Dv(12)-Dv(10))+Q23*(Dv(8)-Dv(3))+Q24*(Dv(6)-Dv(5)) &
                       +Q25*(Dv(4)-Dv(7))+Q26*(Dv(13)+Dv(14))+Q27*(Dv(22)-Dv(24)) &
                       -Q28*(Dv(20)+Dv(15))-Q29*(Dv(17)+Dv(18))+Q30*(Dv(16)+Dv(19)) &
                       +Q31*(Dv(26)-Dv(25))+Q32*(Dv(28)-Dv(27)))
      DuCoarse(Idf(6)) = DuCoarse(Idf(6)) + Dtn(2)*Dtf(2)*(Q21*(Dv(4)-Dv(3)) &
                       +Q22*(Dv(9)-Dv(11))+Q23*(Dv(2)-Dv(5))+Q24*(Dv(8)-Dv(7)) &
                       +Q25*(Dv(6)-Dv(1))+Q26*(Dv(15)+Dv(16))+Q27*(Dv(23)-Dv(21)) &
                       -Q28*(Dv(14)+Dv(17))-Q29*(Dv(19)+Dv(20))+Q30*(Dv(13)+Dv(18)) &
                       +Q31*(Dv(27)-Dv(26))+Q32*(Dv(25)-Dv(28)))
      DuCoarse(Idf(7)) = DuCoarse(Idf(7)) + Dtn(3)*Dtf(3)*(Q21*(Dv(6)-Dv(5)) &
                       +Q22*(Dv(10)-Dv(12))+Q23*(Dv(4)-Dv(7))+Q24*(Dv(2)-Dv(1)) &
                       +Q25*(Dv(8)-Dv(3))+Q26*(Dv(17)+Dv(18))-Q27*(Dv(22)-Dv(24)) &
                       -Q28*(Dv(16)+Dv(19))-Q29*(Dv(13)+Dv(14))+Q30*(Dv(15)+Dv(20)) &
                       +Q31*(Dv(28)-Dv(27))+Q32*(Dv(26)-Dv(25)))
      DuCoarse(Idf(8)) = DuCoarse(Idf(8)) + Dtn(4)*Dtf(4)*(Q21*(Dv(8)-Dv(7)) &
                       +Q22*(Dv(11)-Dv(9))+Q23*(Dv(6)-Dv(1))+Q24*(Dv(4)-Dv(3)) &
                       +Q25*(Dv(2)-Dv(5))+Q26*(Dv(19)+Dv(20))+Q27*(Dv(21)-Dv(23)) &
                       -Q28*(Dv(13)+Dv(18))-Q29*(Dv(15)+Dv(16))+Q30*(Dv(14)+Dv(17)) &
                       +Q31*(Dv(25)-Dv(28))+Q32*(Dv(27)-Dv(26)))

      DuCoarse(Idf(9)) = 0.25_DP*(Dv(25)+Dv(26)+Dv(27)+Dv(28))

      DuCoarse(Idf(10)) = &
                        -Q41*(Dv(1)+Dv(2)+Dv(3)+Dv(4)+Dv(5)+Dv(6)+Dv(7)+Dv(8)) &
                        +Q42*(Dv(9)+Dv(10)+Dv(11)+Dv(12)) &
                        +Q43*(-Dv(13)+Dv(14)-Dv(15)+Dv(16) &
                              -Dv(17)+Dv(18)-Dv(19)+Dv(20)) &
                        +Q44*(Dv(21)+Dv(22)+Dv(23)+Dv(24)) &
                        +Q45*(Dv(29)+Dv(30)+Dv(31)+Dv(32))

    end do

    ! That is it

  end subroutine

  ! ***************************************************************************
  ! Support for DGT1 element
  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_prolUniformDGT1_double (DuCoarse,DuFine, &
               IneighboursAtElementFine, NELcoarse)

!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! DG_T1, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Coarse grid vector
  real(DP), dimension(:), intent(in) :: DuCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Fine grid vector
  real(DP), dimension(:), intent(out) :: DuFine
!</output>

!</subroutine>

  ! local variables
  integer :: iel
  integer :: IELH1,IELH2,IELH3,IELH4
  real(DP) :: duh, ddx, ddy

    ! Loop over the elements
    do iel=1,NELCoarse

      ! Get the four child elements of the coarse grid element
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(2,IELH3)

      ! Put the value on the coarse grid into all four child
      ! elements
      duh = DuCoarse(3*(iel-1)+1)
      ddx = 0.5_dp*DuCoarse(3*(iel-1)+2)
      ddy = 0.5_dp*DuCoarse(3*(iel-1)+3)

      DuFine(3*(IELH1-1)+1) = duh -ddx-ddy
      DuFine(3*(IELH1-1)+2) = ddx
      DuFine(3*(IELH1-1)+3) = ddy

      DuFine(3*(IELH2-1)+1) = duh +ddx-ddy
      DuFine(3*(IELH2-1)+2) = ddy
      DuFine(3*(IELH2-1)+3) = -ddx

      DuFine(3*(IELH3-1)+1) = duh +ddx+ddy
      DuFine(3*(IELH3-1)+2) = -ddx
      DuFine(3*(IELH3-1)+3) = -ddy

      DuFine(3*(IELH4-1)+1) = duh -ddx+ddy
      DuFine(3*(IELH4-1)+2) = -ddy
      DuFine(3*(IELH4-1)+3) = ddx

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_restUniformDGT1_double (DuCoarse,DuFine, &
               IneighboursAtElementFine, NELcoarse)

!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! DG_T1, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  integer :: iel
  integer :: IELH1,IELH2,IELH3,IELH4

    ! The information that was "distributed" in the prolongation has to
    ! be "collected".

    ! Loop over the elements to collect the missing additive contributions:
    do iel=1,NELcoarse

      ! Get the elements on the fine grid that are children of the
      ! coarse grid element
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(2,IELH3)

      ! Sum up the values in these nodes to get the
      ! value in the coarse grid element
      DuCoarse(3*(iel-1)+1)= DuFine(3*(IELH1-1)+1)+DuFine(3*(IELH2-1)+1)+DuFine(3*(IELH3-1)+1)+DuFine(3*(IELH4-1)+1)
      DuCoarse(3*(iel-1)+2)= 0.5_dp*(-DuFine(3*(IELH1-1)+1)+DuFine(3*(IELH2-1)+1)+DuFine(3*(IELH3-1)+1)&
       -DuFine(3*(IELH4-1)+1) + DuFine(3*(IELH1-1)+2)-DuFine(3*(IELH2-1)+3)-DuFine(3*(IELH3-1)+2)+DuFine(3*(IELH4-1)+3) )
      DuCoarse(3*(iel-1)+3)= 0.5_dp*(-DuFine(3*(IELH1-1)+1)-DuFine(3*(IELH2-1)+1)+DuFine(3*(IELH3-1)+1)&
       +DuFine(3*(IELH4-1)+1) + DuFine(3*(IELH1-1)+3)+DuFine(3*(IELH2-1)+2)-DuFine(3*(IELH3-1)+3)-DuFine(3*(IELH4-1)+2) )

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_interpUniformDGT1_double (DuCoarse,DuFine, &
               IneighboursAtElementFine, NELcoarse)

!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! DG_T1, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  real(DP), parameter :: Q4 = 0.25_DP

  integer :: iel
  integer :: IELH1,IELH2,IELH3,IELH4

    ! The information that was "distributed" in the prolongation has to
    ! be "collected".

    ! Loop over the elements to collect the missing additive contributions:
    do iel=1,NELcoarse

      ! Get the elements on the fine grid that are children of the
      ! coarse grid element
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(3,IELH1)
      IELH3 = IneighboursAtElementFine(3,IELH2)
      IELH4 = IneighboursAtElementFine(3,IELH3)

      ! From all four values, build the mean and use that as
      ! value of the coarse grid element
      DuCoarse(3*(iel-1)+1)= 0.25_dp*(DuFine(3*(IELH1-1)+1)+DuFine(3*(IELH2-1)+1)+DuFine(3*(IELH3-1)+1)+DuFine(3*(IELH4-1)+1))
      DuCoarse(3*(iel-1)+2)= 0.5_dp*(DuFine(3*(IELH1-1)+2)+DuFine(3*(IELH2-1)+3)-DuFine(3*(IELH3-1)+2)-DuFine(3*(IELH4-1)+3))
      DuCoarse(3*(iel-1)+3)= 0.5_dp*(DuFine(3*(IELH1-1)+3)-DuFine(3*(IELH2-1)+2)-DuFine(3*(IELH3-1)+3)+DuFine(3*(IELH4-1)+2))

    end do

  end subroutine

  ! ***************************************************************************
  ! Support for DGT2 element
  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_prolUniformDGT2_double (DuCoarse,DuFine, &
               IneighboursAtElementFine, NELcoarse)

!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! DG_T2, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Coarse grid vector
  real(DP), dimension(:), intent(in) :: DuCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Fine grid vector
  real(DP), dimension(:), intent(out) :: DuFine
!</output>

!</subroutine>

  ! local variables
  integer :: iel
  integer :: IELH1,IELH2,IELH3,IELH4
  real(DP) :: duh, ddx, ddy, ddxx, ddxy, ddyy

    ! Loop over the elements
    do iel=1,NELCoarse

      ! Get the four child elements of the coarse grid element
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(2,IELH3)

      ! Put the value on the coarse grid into all four child
      ! elements
      duh = DuCoarse(6*(iel-1)+1)
      ddx = DuCoarse(6*(iel-1)+2)
      ddy = DuCoarse(6*(iel-1)+3)
      ddxx = DuCoarse(6*(iel-1)+4)
      ddyy = DuCoarse(6*(iel-1)+5)
      ddxy = DuCoarse(6*(iel-1)+6)

      DuFine(6*(IELH1-1)+1) = duh-0.5_dp*ddx-0.5_dp*ddy !-0.5_dp*ddxx-0.5_dp*ddyy+0.25_dp*ddxy
      DuFine(6*(IELH1-1)+2) = 0.5_dp*(ddx-0.5_dp*ddxx-0.5_dp*ddxy)
      DuFine(6*(IELH1-1)+3) = 0.5_dp*(ddy-0.5_dp*ddyy-0.5_dp*ddxy)
      DuFine(6*(IELH1-1)+4) = 0.25_dp*ddxx
      DuFine(6*(IELH1-1)+5) = 0.25_dp*ddyy
      DuFine(6*(IELH1-1)+6) = 0.25_dp*ddxy

      DuFine(6*(IELH2-1)+1) = duh+0.5_dp*ddx-0.5_dp*ddy !-0.5_dp*ddxx-0.5_dp*ddyy-0.25_dp*ddxy
      DuFine(6*(IELH2-1)+2) = 0.5_dp*(ddy-0.5_dp*ddyy+0.5_dp*ddxy)
      DuFine(6*(IELH2-1)+3) = -0.5_dp*(ddx+0.5_dp*ddxx-0.5_dp*ddxy)
      DuFine(6*(IELH2-1)+4) = 0.25_dp*ddyy
      DuFine(6*(IELH2-1)+5) = 0.25_dp*ddxx
      DuFine(6*(IELH2-1)+6) = -0.25_dp*ddxy

      DuFine(6*(IELH3-1)+1) = duh+0.5_dp*ddx+0.5_dp*ddy !-0.5_dp*ddxx-0.5_dp*ddyy+0.25_dp*ddxy
      DuFine(6*(IELH3-1)+2) = -0.5_dp*(ddx+0.5_dp*ddxx+0.5_dp*ddxy)
      DuFine(6*(IELH3-1)+3) = -0.5_dp*(ddy+0.5_dp*ddyy+0.5_dp*ddxy)
      DuFine(6*(IELH3-1)+4) = 0.25_dp*ddxx
      DuFine(6*(IELH3-1)+5) = 0.25_dp*ddyy
      DuFine(6*(IELH3-1)+6) = 0.25_dp*ddxy

      DuFine(6*(IELH4-1)+1) = duh-0.5_dp*ddx+0.5_dp*ddy !-0.5_dp*ddxx-0.5_dp*ddyy-0.25_dp*ddxy
      DuFine(6*(IELH4-1)+2) = -0.5_dp*(ddy+0.5_dp*ddyy-0.5_dp*ddxy)
      DuFine(6*(IELH4-1)+3) = 0.5_dp*(ddx-0.5_dp*ddxx+0.5_dp*ddxy)
      DuFine(6*(IELH4-1)+4) = 0.25_dp*ddyy
      DuFine(6*(IELH4-1)+5) = 0.25_dp*ddxx
      DuFine(6*(IELH4-1)+6) = -0.25_dp*ddxy

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_restUniformDGT2_double (DuCoarse,DuFine, &
               IneighboursAtElementFine, NELcoarse)

!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! DG_T0, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  integer :: iel
  integer :: IELH1,IELH2,IELH3,IELH4
  real(dp) :: dh1, dh2, dh3, dh4

    ! The information that was "distributed" in the prolongation has to
    ! be "collected".

    ! Loop over the elements to collect the missing additive contributions:
    do iel=1,NELcoarse

      ! Get the elements on the fine grid that are children of the
      ! coarse grid element
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(2,IELH3)

      ! Sum up the values in these nodes to get the
      ! value in the coarse grid element
      DuCoarse(6*(iel-1)+1)= DuFine(6*(IELH1-1)+1)+DuFine(6*(IELH2-1)+1)+DuFine(6*(IELH3-1)+1)+DuFine(6*(IELH4-1)+1)

      dh1 = -DuFine(6*(IELH1-1)+1)+DuFine(6*(IELH1-1)+2)
      dh2 = +DuFine(6*(IELH2-1)+1)-DuFine(6*(IELH2-1)+3)
      dh3 = +DuFine(6*(IELH3-1)+1)-DuFine(6*(IELH3-1)+2)
      dh4 = -DuFine(6*(IELH4-1)+1)+DuFine(6*(IELH4-1)+3)
      DuCoarse(6*(iel-1)+2)= 0.5_dp*( dh1+dh2+dh3+dh4 )

      dh1 = -DuFine(6*(IELH1-1)+1)+DuFine(6*(IELH1-1)+3)
      dh2 = -DuFine(6*(IELH2-1)+1)+DuFine(6*(IELH2-1)+2)
      dh3 = +DuFine(6*(IELH3-1)+1)-DuFine(6*(IELH3-1)+3)
      dh4 = +DuFine(6*(IELH4-1)+1)-DuFine(6*(IELH4-1)+2)
      DuCoarse(6*(iel-1)+3)= 0.5_dp*( dh1+dh2+dh3+dh4 )

      dh1 = -DuFine(6*(IELH1-1)+2)+DuFine(6*(IELH1-1)+4)
      dh2 = -DuFine(6*(IELH2-1)+3)+DuFine(6*(IELH2-1)+5)
      dh3 = -DuFine(6*(IELH3-1)+2)+DuFine(6*(IELH3-1)+4)
      dh4 = -DuFine(6*(IELH4-1)+3)+DuFine(6*(IELH4-1)+5)
      DuCoarse(6*(iel-1)+4)= 0.25_dp*( dh1+dh2+dh3+dh4 )

      dh1 = -DuFine(6*(IELH1-1)+3)+DuFine(6*(IELH1-1)+5)
      dh2 = -DuFine(6*(IELH2-1)+2)+DuFine(6*(IELH2-1)+4)
      dh3 = -DuFine(6*(IELH3-1)+3)+DuFine(6*(IELH3-1)+5)
      dh4 = -DuFine(6*(IELH4-1)+2)+DuFine(6*(IELH4-1)+4)
      DuCoarse(6*(iel-1)+5)= 0.25_dp*( dh1+dh2+dh3+dh4 )

      dh1 = -DuFine(6*(IELH1-1)+2)-DuFine(6*(IELH1-1)+3)+DuFine(6*(IELH1-1)+6)
      dh2 = +DuFine(6*(IELH2-1)+2)+DuFine(6*(IELH2-1)+3)-DuFine(6*(IELH2-1)+6)
      dh3 = -DuFine(6*(IELH3-1)+2)-DuFine(6*(IELH3-1)+3)+DuFine(6*(IELH3-1)+6)
      dh4 = +DuFine(6*(IELH4-1)+2)+DuFine(6*(IELH4-1)+3)-DuFine(6*(IELH4-1)+6)
      DuCoarse(6*(iel-1)+6)= 0.25_dp*( dh1+dh2+dh3+dh4 )


    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_interpUniformDGT2_double (DuCoarse,DuFine, &
               IneighboursAtElementFine, NELcoarse)

!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! DG_T2, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  real(DP), parameter :: Q4 = 0.25_DP

  integer :: iel
  integer :: IELH1,IELH2,IELH3,IELH4
  real(dp) :: dh1, dh2, dh3, dh4

    ! The information that was "distributed" in the prolongation has to
    ! be "collected".

    ! Loop over the elements to collect the missing additive contributions:
    do iel=1,NELcoarse

      ! Get the elements on the fine grid that are children of the
      ! coarse grid element
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(2,IELH3)

      ! From all four values, build the mean and use that as
      ! value of the coarse grid element
      DuCoarse(6*(iel-1)+1)= 0.25_dp*(DuFine(6*(IELH1-1)+1)+DuFine(6*(IELH2-1)+1)+DuFine(6*(IELH3-1)+1)+DuFine(6*(IELH4-1)+1))

      dh1 = DuFine(6*(IELH1-1)+2)+DuFine(6*(IELH1-1)+4)+DuFine(6*(IELH1-1)+6)
      dh2 = DuFine(6*(IELH2-1)+3)+DuFine(6*(IELH2-1)+5)+DuFine(6*(IELH2-1)+6)
      dh3 = DuFine(6*(IELH3-1)+2)+DuFine(6*(IELH3-1)+4)+DuFine(6*(IELH3-1)+6)
      dh4 = DuFine(6*(IELH4-1)+3)+DuFine(6*(IELH4-1)+5)+DuFine(6*(IELH4-1)+6)
      DuCoarse(6*(iel-1)+2)= 0.5_dp*( dh1-dh2-dh3+dh4 )

      dh1 = DuFine(6*(IELH1-1)+3)+DuFine(6*(IELH1-1)+5)+DuFine(6*(IELH1-1)+6)
      dh2 = DuFine(6*(IELH2-1)+2)+DuFine(6*(IELH2-1)+4)+DuFine(6*(IELH2-1)+6)
      dh3 = DuFine(6*(IELH3-1)+3)+DuFine(6*(IELH3-1)+5)+DuFine(6*(IELH3-1)+6)
      dh4 = DuFine(6*(IELH4-1)+2)+DuFine(6*(IELH4-1)+4)+DuFine(6*(IELH4-1)+6)
      DuCoarse(6*(iel-1)+3)= 0.5_dp*( dh1+dh2-dh3-dh4 )

      DuCoarse(6*(iel-1)+4) = DuFine(6*(IELH1-1)+4)+DuFine(6*(IELH2-1)+5)+DuFine(6*(IELH3-1)+4)+DuFine(6*(IELH4-1)+5)

      DuCoarse(6*(iel-1)+5) = DuFine(6*(IELH1-1)+5)+DuFine(6*(IELH2-1)+4)+DuFine(6*(IELH3-1)+5)+DuFine(6*(IELH4-1)+4)

      DuCoarse(6*(iel-1)+6) = DuFine(6*(IELH1-1)+6)-DuFine(6*(IELH2-1)+6)+DuFine(6*(IELH3-1)+6)-DuFine(6*(IELH4-1)+6)

    end do

  end subroutine

  ! ***************************************************************************
  ! Support for DGQ1 element
  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_prolUniformDGQ1_double (DuCoarse,DuFine, &
               IneighboursAtElementFine, NELcoarse)

!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! DG_Q1, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Coarse grid vector
  real(DP), dimension(:), intent(in) :: DuCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Fine grid vector
  real(DP), dimension(:), intent(out) :: DuFine
!</output>

!</subroutine>

  ! local variables
  integer :: iel
  integer :: IELH1,IELH2,IELH3,IELH4
  real(DP) :: da, db, dc, dd

    ! Loop over the elements
    do iel=1,NELCoarse

      ! Get the four child elements of the coarse grid element
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(2,IELH3)

      ! Put the value on the coarse grid into all four child
      ! elements
      da = DuCoarse(4*(iel-1)+1)
      db = DuCoarse(4*(iel-1)+2)
      dc = DuCoarse(4*(iel-1)+3)
      dd = DuCoarse(4*(iel-1)+4)

      DuFine(4*(IELH1-1)+1) = da
      DuFine(4*(IELH1-1)+2) = 0.5_dp*(da+db)
      DuFine(4*(IELH1-1)+3) = 0.25_dp*(da+db+dc+dd)
      DuFine(4*(IELH1-1)+4) = 0.5_dp*(da+dd)

      DuFine(4*(IELH2-1)+1) = db
      DuFine(4*(IELH2-1)+2) = 0.5_dp*(db+dc)
      DuFine(4*(IELH2-1)+3) = 0.25_dp*(da+db+dc+dd)
      DuFine(4*(IELH2-1)+4) = 0.5_dp*(da+db)

      DuFine(4*(IELH3-1)+1) = dc
      DuFine(4*(IELH3-1)+2) = 0.5_dp*(dc+dd)
      DuFine(4*(IELH3-1)+3) = 0.25_dp*(da+db+dc+dd)
      DuFine(4*(IELH3-1)+4) = 0.5_dp*(db+dc)

      DuFine(4*(IELH4-1)+1) = dd
      DuFine(4*(IELH4-1)+2) = 0.5_dp*(dd+da)
      DuFine(4*(IELH4-1)+3) = 0.25_dp*(da+db+dc+dd)
      DuFine(4*(IELH4-1)+4) = 0.5_dp*(dc+dd)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_restUniformDGQ1_double (DuCoarse,DuFine, &
               IneighboursAtElementFine, NELcoarse)

!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! DG_Q1, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  integer :: iel
  integer :: IELH1,IELH2,IELH3,IELH4
  real(dp) :: dh1, dh2, dh3, dh4

    ! The information that was "distributed" in the prolongation has to
    ! be "collected".

    ! Loop over the elements to collect the missing additive contributions:
    do iel=1,NELcoarse

      ! Get the elements on the fine grid that are children of the
      ! coarse grid element
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(2,IELH3)

      ! Sum up the values in these nodes to get the
      ! value in the coarse grid element
      dh1 = DuFine(4*(IELH1-1)+1) + 0.5_dp*DuFine(4*(IELH1-1)+2) + 0.25_dp*DuFine(4*(IELH1-1)+3) + 0.5_dp*DuFine(4*(IELH1-1)+4)
      dh2 = 0.25_dp*DuFine(4*(IELH2-1)+3) + 0.5_dp*DuFine(4*(IELH2-1)+4)
      dh3 = 0.25_dp*DuFine(4*(IELH3-1)+3)
      dh4 = 0.5_dp*DuFine(4*(IELH4-1)+2) + 0.25_dp*DuFine(4*(IELH4-1)+3)
      DuCoarse(4*(iel-1)+1) = dh1 + dh2 + dh3 + dh4

      dh1 = DuFine(4*(IELH2-1)+1) + 0.5_dp*DuFine(4*(IELH2-1)+2) + 0.25_dp*DuFine(4*(IELH2-1)+3) + 0.5_dp*DuFine(4*(IELH2-1)+4)
      dh2 = 0.25_dp*DuFine(4*(IELH3-1)+3) + 0.5_dp*DuFine(4*(IELH3-1)+4)
      dh3 = 0.25_dp*DuFine(4*(IELH4-1)+3)
      dh4 = 0.5_dp*DuFine(4*(IELH1-1)+2) + 0.25_dp*DuFine(4*(IELH1-1)+3)
      DuCoarse(4*(iel-1)+2) = dh1 + dh2 + dh3 + dh4

      dh1 = DuFine(4*(IELH3-1)+1) + 0.5_dp*DuFine(4*(IELH3-1)+2) + 0.25_dp*DuFine(4*(IELH3-1)+3) + 0.5_dp*DuFine(4*(IELH3-1)+4)
      dh2 = 0.25_dp*DuFine(4*(IELH4-1)+3) + 0.5_dp*DuFine(4*(IELH4-1)+4)
      dh3 = 0.25_dp*DuFine(4*(IELH1-1)+3)
      dh4 = 0.5_dp*DuFine(4*(IELH2-1)+2) + 0.25_dp*DuFine(4*(IELH2-1)+3)
      DuCoarse(4*(iel-1)+3) = dh1 + dh2 + dh3 + dh4

      dh1 = DuFine(4*(IELH4-1)+1) + 0.5_dp*DuFine(4*(IELH4-1)+2) + 0.25_dp*DuFine(4*(IELH4-1)+3) + 0.5_dp*DuFine(4*(IELH4-1)+4)
      dh2 = 0.25_dp*DuFine(4*(IELH1-1)+3) + 0.5_dp*DuFine(4*(IELH1-1)+4)
      dh3 = 0.25_dp*DuFine(4*(IELH2-1)+3)
      dh4 = 0.5_dp*DuFine(4*(IELH3-1)+2) + 0.25_dp*DuFine(4*(IELH3-1)+3)
      DuCoarse(4*(iel-1)+4) = dh1 + dh2 + dh3 + dh4

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_interpUniformDGQ1_double (DuCoarse,DuFine, &
               IneighboursAtElementFine, NELcoarse)

!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! DG_Q1, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  real(DP), parameter :: Q4 = 0.25_DP

  integer :: iel
  integer :: IELH1,IELH2,IELH3,IELH4

    ! The information that was "distributed" in the prolongation has to
    ! be "collected".

    ! Loop over the elements to collect the missing additive contributions:
    do iel=1,NELcoarse

      ! Get the elements on the fine grid that are children of the
      ! coarse grid element
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(2,IELH3)

      ! From all four values, build the mean and use that as
      ! value of the coarse grid element
      DuCoarse(4*(iel-1)+1) = DuFine(4*(IELH1-1)+1)
      DuCoarse(4*(iel-1)+2) = DuFine(4*(IELH2-1)+1)
      DuCoarse(4*(iel-1)+3) = DuFine(4*(IELH3-1)+1)
      DuCoarse(4*(iel-1)+4) = DuFine(4*(IELH4-1)+1)

    end do

  end subroutine


  ! ***************************************************************************
  ! Support for DGQ2 element
  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_prolUniformDGQ2_double (DuCoarse,DuFine, &
               IneighboursAtElementFine, NELcoarse)

!<description>
  ! Prolongate a solution vector from a coarse grid to a fine grid.
  ! DG_Q2, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Coarse grid vector
  real(DP), dimension(:), intent(in) :: DuCoarse

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Fine grid vector
  real(DP), dimension(:), intent(out) :: DuFine
!</output>

!</subroutine>

  ! local variables
  integer :: iel
  integer :: IELH1,IELH2,IELH3,IELH4
  real(DP) :: da, db, dc, dd, de, df, dg, dh, di

    ! Loop over the elements
    do iel=1,NELCoarse

      ! Get the four child elements of the coarse grid element
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(2,IELH3)

      ! Put the value on the coarse grid into all four child
      ! elements
      da = DuCoarse(9*(iel-1)+1)
      db = DuCoarse(9*(iel-1)+2)
      dc = DuCoarse(9*(iel-1)+3)
      dd = DuCoarse(9*(iel-1)+4)
      de = DuCoarse(9*(iel-1)+5)
      df = DuCoarse(9*(iel-1)+6)
      dg = DuCoarse(9*(iel-1)+7)
      dh = DuCoarse(9*(iel-1)+8)
      di = DuCoarse(9*(iel-1)+9)

      DuFine(9*(IELH1-1)+1) = da
      DuFine(9*(IELH1-1)+2) = de
      DuFine(9*(IELH1-1)+3) = di
      DuFine(9*(IELH1-1)+4) = dh
      DuFine(9*(IELH1-1)+5) = 0.375_dp*da+0.75_dp*de-0.125_dp*db
      DuFine(9*(IELH1-1)+6) = 0.375_dp*de+0.75_dp*di-0.125_dp*dg
      DuFine(9*(IELH1-1)+7) = 0.375_dp*dh+0.75_dp*di-0.125_dp*df
      DuFine(9*(IELH1-1)+8) = 0.375_dp*da+0.75_dp*dh-0.125_dp*dd
      DuFine(9*(IELH1-1)+9) = 9.0_dp/64.0_dp*da + 1.0_dp/64.0_dp*dc &
                              -3.0_dp/64.0_dp*(db+dd) &
                              +9.0_dp/32.0_dp*(de+dh) &
                              -3.0_dp/32.0_dp*(df+dg) &
                              +9.0_dp/16.0_dp*di

      DuFine(9*(IELH2-1)+1) = db
      DuFine(9*(IELH2-1)+2) = df
      DuFine(9*(IELH2-1)+3) = di
      DuFine(9*(IELH2-1)+4) = de
      DuFine(9*(IELH2-1)+5) = 0.375_dp*db+0.75_dp*df-0.125_dp*dc
      DuFine(9*(IELH2-1)+6) = 0.375_dp*df+0.75_dp*di-0.125_dp*dh
      DuFine(9*(IELH2-1)+7) = 0.375_dp*de+0.75_dp*di-0.125_dp*dg
      DuFine(9*(IELH2-1)+8) = 0.375_dp*db+0.75_dp*de-0.125_dp*da
      DuFine(9*(IELH2-1)+9) = 9.0_dp/64.0_dp*db + 1.0_dp/64.0_dp*dd &
                              -3.0_dp/64.0_dp*(da+dc) &
                              +9.0_dp/32.0_dp*(de+df) &
                              -3.0_dp/32.0_dp*(dg+dh) &
                              +9.0_dp/16.0_dp*di

      DuFine(9*(IELH3-1)+1) = dc
      DuFine(9*(IELH3-1)+2) = dg
      DuFine(9*(IELH3-1)+3) = di
      DuFine(9*(IELH3-1)+4) = df
      DuFine(9*(IELH3-1)+5) = 0.375_dp*dc+0.75_dp*dg-0.125_dp*dd
      DuFine(9*(IELH3-1)+6) = 0.375_dp*dg+0.75_dp*di-0.125_dp*de
      DuFine(9*(IELH3-1)+7) = 0.375_dp*df+0.75_dp*di-0.125_dp*dh
      DuFine(9*(IELH3-1)+8) = 0.375_dp*dc+0.75_dp*df-0.125_dp*db
      DuFine(9*(IELH3-1)+9) = 9.0_dp/64.0_dp*dc + 1.0_dp/64.0_dp*da &
                              -3.0_dp/64.0_dp*(db+dd) &
                              +9.0_dp/32.0_dp*(df+dg) &
                              -3.0_dp/32.0_dp*(de+dh) &
                              +9.0_dp/16.0_dp*di

      DuFine(9*(IELH4-1)+1) = dd
      DuFine(9*(IELH4-1)+2) = dh
      DuFine(9*(IELH4-1)+3) = di
      DuFine(9*(IELH4-1)+4) = dg
      DuFine(9*(IELH4-1)+5) = 0.375_dp*dd+0.75_dp*dh-0.125_dp*da
      DuFine(9*(IELH4-1)+6) = 0.375_dp*dh+0.75_dp*di-0.125_dp*df
      DuFine(9*(IELH4-1)+7) = 0.375_dp*dg+0.75_dp*di-0.125_dp*de
      DuFine(9*(IELH4-1)+8) = 0.375_dp*dd+0.75_dp*dg-0.125_dp*dc
      DuFine(9*(IELH4-1)+9) = 9.0_dp/64.0_dp*dd + 1.0_dp/64.0_dp*db &
                              -3.0_dp/64.0_dp*(da+dc) &
                              +9.0_dp/32.0_dp*(dg+dh) &
                              -3.0_dp/32.0_dp*(de+df) &
                              +9.0_dp/16.0_dp*di

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_restUniformDGQ2_double (DuCoarse,DuFine, &
               IneighboursAtElementFine, NELcoarse)

!<description>
  ! Restricts a RHS vector from a fine grid to a coarse grid.
  ! DG_Q2, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  integer :: iel
  integer :: IELH1,IELH2,IELH3,IELH4
  real(dp) :: dh1, dh2, dh3, dh4

    ! The information that was "distributed" in the prolongation has to
    ! be "collected".

    ! Loop over the elements to collect the missing additive contributions:
    do iel=1,NELcoarse

      ! Get the elements on the fine grid that are children of the
      ! coarse grid element
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(2,IELH3)

      ! Sum up the values in these nodes to get the
      ! value in the coarse grid element
      dh1 = DuFine(9*(IELH1-1)+1) + 0.375_dp*DuFine(9*(IELH1-1)+5) &
       + 0.375_dp*DuFine(9*(IELH1-1)+8) + 9.0_dp/64.0_dp*DuFine(9*(IELH1-1)+9)
      dh2 = -0.125_dp*DuFine(9*(IELH2-1)+8) - 3.0_dp/64.0_dp*DuFine(9*(IELH2-1)+9)
      dh3 = 1.0_dp/64.0_dp*DuFine(9*(IELH3-1)+9)
      dh4 = -0.125_dp*DuFine(9*(IELH4-1)+5) - 3.0_dp/64.0_dp*DuFine(9*(IELH4-1)+9)
      DuCoarse(9*(iel-1)+1) = dh1 + dh2 + dh3 + dh4

      dh1 = DuFine(9*(IELH2-1)+1) + 0.375_dp*DuFine(9*(IELH2-1)+5) &
       + 0.375_dp*DuFine(9*(IELH2-1)+8) + 9.0_dp/64.0_dp*DuFine(9*(IELH2-1)+9)
      dh2 = -0.125_dp*DuFine(9*(IELH3-1)+8) - 3.0_dp/64.0_dp*DuFine(9*(IELH3-1)+9)
      dh3 = 1.0_dp/64.0_dp*DuFine(9*(IELH4-1)+9)
      dh4 = -0.125_dp*DuFine(9*(IELH1-1)+5) - 3.0_dp/64.0_dp*DuFine(9*(IELH1-1)+9)
      DuCoarse(9*(iel-1)+2) = dh1 + dh2 + dh3 + dh4

      dh1 = DuFine(9*(IELH3-1)+1) + 0.375_dp*DuFine(9*(IELH3-1)+5) &
       + 0.375_dp*DuFine(9*(IELH3-1)+8) + 9.0_dp/64.0_dp*DuFine(9*(IELH3-1)+9)
      dh2 = -0.125_dp*DuFine(9*(IELH4-1)+8) - 3.0_dp/64.0_dp*DuFine(9*(IELH4-1)+9)
      dh3 = 1.0_dp/64.0_dp*DuFine(9*(IELH1-1)+9)
      dh4 = -0.125_dp*DuFine(9*(IELH2-1)+5) - 3.0_dp/64.0_dp*DuFine(9*(IELH2-1)+9)
      DuCoarse(9*(iel-1)+3) = dh1 + dh2 + dh3 + dh4

      dh1 = DuFine(9*(IELH4-1)+1) + 0.375_dp*DuFine(9*(IELH4-1)+5) &
       + 0.375_dp*DuFine(9*(IELH4-1)+8) + 9.0_dp/64.0_dp*DuFine(9*(IELH4-1)+9)
      dh2 = -0.125_dp*DuFine(9*(IELH1-1)+8) - 3.0_dp/64.0_dp*DuFine(9*(IELH1-1)+9)
      dh3 = 1.0_dp/64.0_dp*DuFine(9*(IELH2-1)+9)
      dh4 = -0.125_dp*DuFine(9*(IELH3-1)+5) - 3.0_dp/64.0_dp*DuFine(9*(IELH3-1)+9)
      DuCoarse(9*(iel-1)+4) = dh1 + dh2 + dh3 + dh4

      dh1 = DuFine(9*(IELH1-1)+2) + 0.75_dp*DuFine(9*(IELH1-1)+5) + 0.375_dp*DuFine(9*(IELH1-1)+6) &
       + 9.0_dp/32.0_dp*DuFine(9*(IELH1-1)+9)
      dh2 = DuFine(9*(IELH2-1)+4) + 0.375_dp*DuFine(9*(IELH2-1)+7) + 0.75_dp*DuFine(9*(IELH2-1)+8) &
       + 9.0_dp/32.0_dp*DuFine(9*(IELH2-1)+9)
      dh3 = -0.125_dp*DuFine(9*(IELH3-1)+6) -3.0_dp/32.0_dp*DuFine(9*(IELH3-1)+9)
      dh4 = -0.125_dp*DuFine(9*(IELH4-1)+7) -3.0_dp/32.0_dp*DuFine(9*(IELH4-1)+9)
      DuCoarse(9*(iel-1)+5) = dh1 + dh2 + dh3 + dh4

      dh1 = DuFine(9*(IELH2-1)+2) + 0.75_dp*DuFine(9*(IELH2-1)+5) + 0.375_dp*DuFine(9*(IELH2-1)+6) &
       + 9.0_dp/32.0_dp*DuFine(9*(IELH2-1)+9)
      dh2 = DuFine(9*(IELH3-1)+4) + 0.375_dp*DuFine(9*(IELH3-1)+7) + 0.75_dp*DuFine(9*(IELH3-1)+8) &
       + 9.0_dp/32.0_dp*DuFine(9*(IELH3-1)+9)
      dh3 = -0.125_dp*DuFine(9*(IELH4-1)+6) -3.0_dp/32.0_dp*DuFine(9*(IELH4-1)+9)
      dh4 = -0.125_dp*DuFine(9*(IELH1-1)+7) -3.0_dp/32.0_dp*DuFine(9*(IELH1-1)+9)
      DuCoarse(9*(iel-1)+6) = dh1 + dh2 + dh3 + dh4

      dh1 = DuFine(9*(IELH3-1)+2) + 0.75_dp*DuFine(9*(IELH3-1)+5) + 0.375_dp*DuFine(9*(IELH3-1)+6) &
       + 9.0_dp/32.0_dp*DuFine(9*(IELH3-1)+9)
      dh2 = DuFine(9*(IELH4-1)+4) + 0.375_dp*DuFine(9*(IELH4-1)+7) + 0.75_dp*DuFine(9*(IELH4-1)+8) &
       + 9.0_dp/32.0_dp*DuFine(9*(IELH4-1)+9)
      dh3 = -0.125_dp*DuFine(9*(IELH1-1)+6) -3.0_dp/32.0_dp*DuFine(9*(IELH1-1)+9)
      dh4 = -0.125_dp*DuFine(9*(IELH2-1)+7) -3.0_dp/32.0_dp*DuFine(9*(IELH2-1)+9)
      DuCoarse(9*(iel-1)+7) = dh1 + dh2 + dh3 + dh4

      dh1 = DuFine(9*(IELH4-1)+2) + 0.75_dp*DuFine(9*(IELH4-1)+5) + 0.375_dp*DuFine(9*(IELH4-1)+6) &
       + 9.0_dp/32.0_dp*DuFine(9*(IELH4-1)+9)
      dh2 = DuFine(9*(IELH1-1)+4) + 0.375_dp*DuFine(9*(IELH1-1)+7) + 0.75_dp*DuFine(9*(IELH1-1)+8) &
       + 9.0_dp/32.0_dp*DuFine(9*(IELH1-1)+9)
      dh3 = -0.125_dp*DuFine(9*(IELH2-1)+6) -3.0_dp/32.0_dp*DuFine(9*(IELH2-1)+9)
      dh4 = -0.125_dp*DuFine(9*(IELH3-1)+7) -3.0_dp/32.0_dp*DuFine(9*(IELH3-1)+9)
      DuCoarse(9*(iel-1)+8) = dh1 + dh2 + dh3 + dh4

      dh1 = DuFine(9*(IELH1-1)+3) + 0.75_dp*DuFine(9*(IELH1-1)+6) + 0.75_dp*DuFine(9*(IELH1-1)+7) &
       + 9.0_dp/16.0_dp*DuFine(9*(IELH1-1)+9)
      dh2 = DuFine(9*(IELH2-1)+3) + 0.75_dp*DuFine(9*(IELH2-1)+6) + 0.75_dp*DuFine(9*(IELH2-1)+7) &
       + 9.0_dp/16.0_dp*DuFine(9*(IELH2-1)+9)
      dh3 = DuFine(9*(IELH3-1)+3) + 0.75_dp*DuFine(9*(IELH3-1)+6) + 0.75_dp*DuFine(9*(IELH3-1)+7) &
       + 9.0_dp/16.0_dp*DuFine(9*(IELH3-1)+9)
      dh4 = DuFine(9*(IELH4-1)+3) + 0.75_dp*DuFine(9*(IELH4-1)+6) + 0.75_dp*DuFine(9*(IELH4-1)+7) &
       + 9.0_dp/16.0_dp*DuFine(9*(IELH4-1)+9)
      DuCoarse(9*(iel-1)+9) = dh1 + dh2 + dh3 + dh4

    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_interpUniformDGQ2_double (DuCoarse,DuFine, &
               IneighboursAtElementFine, NELcoarse)

!<description>
  ! Interpolates a solution vector from a fine grid to a coarse grid.
  ! DG_Q2, uniform triangulation, double precision vector.
!</description>

!<input>
  ! Fine grid vector
  real(DP), dimension(:), intent(in) :: DuFine

  ! IneighboursAtElement array on the fine grid
  integer, dimension(:,:), intent(in) :: IneighboursAtElementFine

  ! Number of elements in the coarse grid
  integer, intent(in) :: NELcoarse
!</input>

!<output>
  ! Coarse grid vector
  real(DP), dimension(:), intent(out) :: DuCoarse
!</output>

!</subroutine>

  ! local variables
  real(DP), parameter :: Q4 = 0.25_DP

  integer :: iel
  integer :: IELH1,IELH2,IELH3,IELH4

    ! The information that was "distributed" in the prolongation has to
    ! be "collected".

    ! Loop over the elements to collect the missing additive contributions:
    do iel=1,NELcoarse

      ! Get the elements on the fine grid that are children of the
      ! coarse grid element
      IELH1 = iel
      IELH2 = IneighboursAtElementFine(2,IELH1)
      IELH3 = IneighboursAtElementFine(2,IELH2)
      IELH4 = IneighboursAtElementFine(2,IELH3)

      ! From all four values, build the mean and use that as
      ! value of the coarse grid element
      DuCoarse(9*(iel-1)+1) = DuFine(9*(IELH1-1)+1)
      DuCoarse(9*(iel-1)+2) = DuFine(9*(IELH2-1)+1)
      DuCoarse(9*(iel-1)+3) = DuFine(9*(IELH3-1)+1)
      DuCoarse(9*(iel-1)+4) = DuFine(9*(IELH4-1)+1)

      DuCoarse(9*(iel-1)+5) = 0.5_dp*( DuFine(9*(IELH1-1)+2) + DuFine(9*(IELH2-1)+4))
      DuCoarse(9*(iel-1)+6) = 0.5_dp*( DuFine(9*(IELH2-1)+2) + DuFine(9*(IELH3-1)+4))
      DuCoarse(9*(iel-1)+7) = 0.5_dp*( DuFine(9*(IELH3-1)+2) + DuFine(9*(IELH4-1)+4))
      DuCoarse(9*(iel-1)+8) = 0.5_dp*( DuFine(9*(IELH4-1)+2) + DuFine(9*(IELH1-1)+4))

      DuCoarse(9*(iel-1)+9) = 0.25_dp*( DuFine(9*(IELH1-1)+3) &
                                       +DuFine(9*(IELH2-1)+3) &
                                       +DuFine(9*(IELH3-1)+3) &
                                       +DuFine(9*(IELH4-1)+3) )
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_initL2Projection (rprojection,r2Lvlmass,rmass)

!<description>
  ! Sets the mass matrices which are needed for L2-projection.
!</description>

!<input>
  ! The 2-Level-Mass matrix
  type(t_matrixScalar), intent(in) :: r2LvlMass

  ! The mass matrix of the fine grid
  type(t_matrixScalar), intent(in) :: rmass
!</input>

!<inputoutput>
  ! The scalar projection structure for which the matrices are to be set.
  type(t_interlevelProjectionScalar), intent(inout) :: rprojection
!</inputoutput>

!</subroutine>

     ! This is an L2-projection
     rprojection%iprojType = MLP_PROJ_TYPE_L2_PROJ

     ! First of all make sure that the mass matrix is a type 7/9 matrix.
!     if(rprojection%drelaxL2 .eq. 0.0_DP) then
!
!       ! Type 7 or 9 is allowed
!       if((rmass%cmatrixFormat .ne. LSYSSC_MATRIX7) .and. &
!          (rmass%cmatrixFormat .ne. LSYSSC_MATRIX9) then
!
!          ! We cannot handle this...
!          call output_line("Mass matrix must be a type 7/9 matrix!", &
!              OU_CLASS_ERROR,OU_MODE_STD,"mlprj_initL2Projection")
!
!          call sys_halt()
!
!       end if
!
!     else

       ! Only type 9 is allowed
       if(rmass%cmatrixFormat .ne. LSYSSC_MATRIX9) then

          ! We cannot handle this...
          call output_line("Mass matrix must be a type 9 matrix!", &
              OU_CLASS_ERROR,OU_MODE_STD,"mlprj_initL2Projection")

          call sys_halt()

       end if

!     end if

     ! Create shared copies of the mass matrices
     call lsyssc_duplicateMatrix(r2LvlMass, rprojection%rmatrix2Lvl, &
                                 LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE)
     call lsyssc_duplicateMatrix(rmass, rprojection%rmatrixMass, &
                                 LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE)

     ! That is it

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_prolScalarL2 (rprojection, rcoarseVector, rfineVector, &
                                 rtempVector)

!<description>
  ! Performs an L2-prolongation of a solution vector on the coarse grid
  ! to a solution vector of the fine grid.
!</description>

!<input>
  ! The t_interlevelProjectionScalar structure that configures the grid transfer
  type(t_interlevelProjectionScalar), intent(in) :: rprojection

  ! Coarse grid vector
  type(t_vectorScalar), intent(inout) :: rcoarseVector

  ! Temporary vector
  type(t_vectorScalar), intent(inout) :: rtempVector
!</input>

!<output>
  ! Fine grid vector
  type(t_vectorScalar), intent(inout) :: rfineVector
!</output>

!</subroutine>

  ! local variables
  real(DP) :: dtol
  integer :: niter, cinfo
  real(DP), dimension(:), pointer :: p_Dtmp

    ! Make sure the temporary vector`s size is sufficient. For the
    ! prolongation, NEQ is sufficient.
    if(rtempVector%NEQ .lt. rfineVector%NEQ) then

      ! The temporary vector is too small...
      call output_line("Temporary vector is too small!", &
          OU_CLASS_ERROR,OU_MODE_STD,"mlprj_prolScalarL2")

      call sys_halt()

    end if

    ! Calculate the product of the 2-Level-Mass matrix and the coarse
    ! mesh solution vector.
    call lsyssc_matVec (rprojection%rmatrix2Lvl, rcoarseVector, &
                              rfineVector, 1.0_DP, 0.0_DP, .false.)

    ! Now get the temporary and fine mesh vector`s data arrays
    call lsyssc_getbase_double(rtempVector, p_Dtmp)

    ! Call the SOR solver
    niter = rprojection%nmaxL2Iter
    dtol = rprojection%depsL2
    call qsol_solveSOR(rprojection%rmatrixMass, rfineVector, p_Dtmp, &
                       cinfo, niter, dtol, rprojection%drelaxL2)

    ! That is it

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_restScalarL2 (rprojection, rcoarseVector, rfineVector, &
                                 rtempVector)

!<description>
  ! Performs an L2-restriction of a defect vector on the fine grid
  ! to a defect vector of the coarse grid.
!</description>

!<input>
  ! The t_interlevelProjectionScalar structure that configures the grid transfer
  type(t_interlevelProjectionScalar), intent(in) :: rprojection

  ! Fine grid vector
  type(t_vectorScalar), intent(inout) :: rfineVector

  ! Temporary vector
  type(t_vectorScalar), intent(inout) :: rtempVector
!</input>

!<output>
  ! Coarse grid vector
  type(t_vectorScalar), intent(inout) :: rcoarseVector
!</output>

!</subroutine>

  ! vector data arrays
  integer :: n, niter, cinfo
  real(DP), dimension(:), pointer :: p_Dtmp, p_Dwork, p_Dfine, p_Dcoarse
  real(DP) :: dtol

  ! variables concerning the mass matrices
  integer :: i,j,k
  real(DP), dimension(:), pointer :: p_Da
  integer, dimension(:), pointer :: p_Kcol, p_Kld, p_Kdiag

    ! Make sure the temporary vector`s size is sufficient. For the
    ! restriction, we need 2*NEQ.
    if(rtempVector%NEQ .lt. 2*rfineVector%NEQ) then

      ! The temporary vector is too small...
      call output_line("Temporary vector is too small!", &
          OU_CLASS_ERROR,OU_MODE_STD,"mlprj_restScalarL2")

      call sys_halt()

    end if

    ! Now get the temporary and fine mesh vector`s data arrays
    call lsyssc_getbase_double(rtempVector, p_Dtmp)
    call lsyssc_getbase_double(rfineVector, p_Dfine)

    ! Get the size of the fine mesh vector
    n = rfineVector%NEQ

    ! Now copy the fine mesh vector into the temporary vector. We will not call
    ! the lsyssc_copyVector routine as the vectors have different length!
    call lsyssc_getbase_double(rfineVector, p_Dfine)
    call lalg_copyVector(p_Dfine,p_Dtmp,n)

    ! Okay, now the fine mesh vector has been copied into p_Dtmp(1:n).
    ! Now we need to assign the work array. As the length of p_Dtmp is
    ! at least 2*n, we can use the rest as a work array:
    p_Dwork => p_Dtmp(n+1:)

    ! Get the arrays from the mass matrix
    call lsyssc_getbase_double(rprojection%rmatrixMass, p_Da)
    call lsyssc_getbase_Kcol(rprojection%rmatrixMass, p_Kcol)
    call lsyssc_getbase_Kld(rprojection%rmatrixMass, p_Kld)
    call lsyssc_getbase_Kdiagonal(rprojection%rmatrixMass, p_Kdiag)

    ! And call the SOR solver
    niter = rprojection%nmaxL2iter
    dtol = rprojection%depsL2
    call qsol_solveSOR(n, p_Kld, p_Kcol, p_Kdiag, p_Da, p_Dtmp, p_Dwork, &
                       cinfo, niter, dtol, rprojection%drelaxL2)

    ! If everything went well, p_Dtmp(1:n) now contains the solution of
    ! the Mass system. The only thing we now need to do is to multiply
    ! the solution by the transposed 2-Level-Mass matrix to get the
    ! restricted coarse mesh vector. Unfortunately, the solution vector`s size
    ! does not match the matrix`s size, so we will perform the matrix-vector
    ! multiplication by hand here:

    ! Get the coarse mesh vector`s data array
    call lsyssc_getbase_double(rcoarseVector, p_Dcoarse)

    ! Get the arrays from the 2-Level-Mass matrix
    call lsyssc_getbase_double(rprojection%rmatrix2Lvl, p_Da)
    call lsyssc_getbase_Kcol(rprojection%rmatrix2Lvl, p_Kcol)
    call lsyssc_getbase_Kld(rprojection%rmatrix2Lvl, p_Kld)

    ! Format the coarse mesh vector
    call lalg_clearVector(p_Dcoarse, rcoarseVector%NEQ)

    ! And perform the transposed matrix-vector multiplication
    do i = 1, rfinevector%NEQ
      do j = p_Kld(i), p_Kld(i+1)-1
        k = p_Kcol(j)
        p_Dcoarse(k) = p_Dcoarse(k) + p_Da(j)*p_Dtmp(i)
      end do
    end do

    ! That is it

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_initMatrixProjection (rprojection,rmatrixProl,rmatrixInterp,&
      rmatrixRest)

!<description>
  ! Sets the prolongation and (optionally) the interpolation matrices which
  ! are needed for matrix-based projection.
!</description>

!<input>
  ! The prolongation matrix
  type(t_matrixScalar), intent(in) :: rmatrixProl

  ! OPTIONAL: The interpolation matrix. If given, the interpolation matrix
  ! is assumed to be stored transposed!
  type(t_matrixScalar), optional, intent(in) :: rmatrixInterp

  ! OPTIONAL: The restriction matrix. If not specified, the restriction is
  ! automatically derived from the prolongation.
  type(t_matrixScalar), optional, intent(in) :: rmatrixRest
!</input>

!<inputoutput>
  ! The scalar projection structure for which the matrix is to be set.
  type(t_interlevelProjectionScalar), intent(inout) :: rprojection
!</inputoutput>

!</subroutine>

    ! This is a matrix-based projection
    rprojection%iprojType = MLP_PROJ_TYPE_MATRIX

    ! Create a shared copy of the prolongation matrix
    call lsyssc_duplicateMatrix(rmatrixProl, rprojection%rmatrixProl, &
                                LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE)

    ! Do we have an interpolation matrix?
    if(present(rmatrixInterp)) &
      call lsyssc_duplicateMatrix(rmatrixInterp, rprojection%rmatrixInterp, &
                                  LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE)

    ! Is there a restriction matrix?
    if(present(rmatrixRest)) &
      call lsyssc_duplicateMatrix(rmatrixRest, rprojection%rmatrixRest, &
                                  LSYSSC_DUP_SHARE, LSYSSC_DUP_SHARE)

    ! That is it

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_initPrjHierarchy (rprjHierarchy,nlmin,nlmax)

!<description>
  ! Initialises a projection hierarchy for projection between levels NLMIN
  ! and NLMAX.
!</description>

!<input>
  ! Minimum level in the hierarchy.
  integer, intent(in) :: nlmin

  ! Maximum level in the hierarchy.
  integer, intent(in) :: nlmax
!</input>

!<output>
  ! Projection hierarchy structure to initialise.
  type(t_interlevelProjectionHier), intent(out) :: rprjHierarchy
!</output>

!</subroutine>

    ! Basic initialisation of the structure.
    rprjHierarchy%nlmin = nlmin
    rprjHierarchy%nlmax = nlmax

    allocate(rprjHierarchy%p_Rprojection(nlmax-nlmin+1))
    allocate(rprjHierarchy%p_RdiscrPointer(nlmax-nlmin+1))

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_initPrjHierarchyLevel (rprjHierarchy,ilevel,rdiscretisation)

!<description>
  ! Initialises level ilevel of the projection hierarchy using the block
  ! discretisation rdiscretisation.
  ! After the initialisation of all levels, the routine
  ! mlprj_commitProjectionHierarchy must be called!
!</description>

!<input>
  ! Number of the level to initialise
  integer, intent(in) :: ilevel

  ! Block discretisation of this level. A pointer to this structure
  ! is written to rprjHierarchy, so the structure must persist until
  ! the hierarchy is released.
  type(t_blockDiscretisation), intent(in),target :: rdiscretisation
!</input>

!<inputoutput>
  ! Projection hierarchy structure.
  type(t_interlevelProjectionHier), intent(inout) :: rprjHierarchy
!</inputoutput>

!</subroutine>

    ! Initialise all levels except for the coarse one, this does not need
    ! initialisation.
    if (ilevel .gt. rprjHierarchy%nlmin) then
      call mlprj_initProjectionDiscr (&
          rprjHierarchy%p_Rprojection(ilevel-rprjHierarchy%nlmin+1),rdiscretisation)
    end if

    ! Remember the discretisation
    rprjHierarchy%p_RdiscrPointer(ilevel-rprjHierarchy%nlmin+1)%p_rdiscretisation => &
        rdiscretisation

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_commitPrjHierarchy (rprjHierarchy)

!<description>
  ! Commits a level hierarchy. This routine must be called after the
  ! initialisation of all levels using mlprj_initHierarchyLevel.
  ! It does some final changes to rprjHierarchy to make it ready to use.
!</description>

!<inputoutput>
  ! Projection hierarchy structure.
  type(t_interlevelProjectionHier), intent(inout), target :: rprjHierarchy
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i,imaxmem
    type(t_blockDiscretisation), pointer :: p_rdiscrCoarse,p_rdiscrFine

    ! Calculate the amount of temp memory needed for the projection
    ! and allocate a temp vector with that size.

    imaxmem = 1
    do i=rprjHierarchy%nlmin+1,rprjHierarchy%nlmax
      ! Pass the system metrices on the coarse/fine grid to
      ! mlprj_getTempMemoryDirect to specify the discretisation structures
      ! of all equations in the PDE there.
      ! This calculates the overall size of the temp vector needed for the
      ! projection between an arbitrary pair of consecutive levels.
      p_rdiscrCoarse => &
          rprjHierarchy%p_RdiscrPointer(i-rprjHierarchy%nlmin)%p_rdiscretisation
      p_rdiscrFine => &
          rprjHierarchy%p_RdiscrPointer(i-rprjHierarchy%nlmin+1)%p_rdiscretisation
      imaxmem = max(imaxmem,mlprj_getTempMemoryDirect (&
          rprjHierarchy%p_Rprojection(i),&
          p_rdiscrCoarse%RspatialDiscr(1:p_rdiscrCoarse%ncomponents),&
          p_rdiscrFine%RspatialDiscr(1:p_rdiscrFine%ncomponents)))
    end do

    ! Set up a scalar temporary vector for the projection.
    if (rprjHierarchy%rtempVector%NEQ .ne. 0) then
      call lsyssc_releaseVector (rprjHierarchy%rtempVector)
    end if
    call lsyssc_createVector (rprjHierarchy%rtempVector,imaxmem,.false.)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_releasePrjHierarchy (rprjHierarchy)

!<description>
  ! Releases a projection hierarchy.
!</description>

!<inputoutput>
  ! Projection hierarchy structure to clean up.
  type(t_interlevelProjectionHier), intent(inout) :: rprjHierarchy
!</inputoutput>

!</subroutine>

    if (associated(rprjHierarchy%p_Rprojection)) then
      call lsyssc_releaseVector(rprjHierarchy%rtempVector)
      deallocate(rprjHierarchy%p_Rprojection)
      nullify(rprjHierarchy%p_Rprojection)
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_performProlongationHier (rprjHierarchy,&
      ifineLevel,rcoarseVector,rfineVector)

!<description>
  ! Performs a prolongation for a given block vector (i.e. a projection
  ! in the primal space where the solution lives) based on a projection
  ! hierarchy.
!</description>

!<input>
  ! Underlying projection hierarchy structure.
  type(t_interlevelProjectionHier), intent(inout) :: rprjHierarchy

  ! Coarse grid vector
  type(t_vectorBlock), intent(inout) :: rcoarseVector

  ! Level of the fine grid vector.
  integer, intent(in) :: ifineLevel
!</input>

!<output>
  ! Fine grid vector
  type(t_vectorBlock), intent(inout) :: rfineVector
!</output>

!</subroutine>

    call mlprj_performProlongation (&
        rprjHierarchy%P_Rprojection(ifineLevel-rprjHierarchy%nlmin+1),&
        rcoarseVector,rfineVector,rprjHierarchy%rtempVector)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_performRestrictionHier (rprjHierarchy,&
      ifineLevel,rcoarseVector,rfineVector)

!<description>
  ! Performs a restriction for a given block vector (i.e. a projection
  ! in the dual space where the RHS vector lives) based on a projection
  ! hierarchy.
!</description>

!<input>
  ! Underlying projection hierarchy structure.
  type(t_interlevelProjectionHier), intent(inout) :: rprjHierarchy

  ! Fine grid vector
  type(t_vectorBlock), intent(inout) :: rfineVector

  ! Level of the fine grid vector.
  integer, intent(in) :: ifineLevel
!</input>

!<output>
  ! Coarse grid vector
  type(t_vectorBlock), intent(inout) :: rcoarseVector
!</output>

!</subroutine>

    call mlprj_performRestriction (&
        rprjHierarchy%P_Rprojection(ifineLevel-rprjHierarchy%nlmin+1),&
        rcoarseVector,rfineVector,rprjHierarchy%rtempVector)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine mlprj_performInterpolationHier (rprjHierarchy,&
      ifineLevel,rcoarseVector,rfineVector)

!<description>
  ! Performs an interpolation for a given block vector (i.e. a projection
  ! in the primal space where the solution vector lives) based on a projection
  ! hierarchy.
!</description>

!<input>
  ! Underlying projection hierarchy structure.
  type(t_interlevelProjectionHier), intent(inout) :: rprjHierarchy

  ! Fine grid vector
  type(t_vectorBlock), intent(inout) :: rfineVector

  ! Level of the fine grid vector.
  integer, intent(in) :: ifineLevel
!</input>

!<output>
  ! Coarse grid vector
  type(t_vectorBlock), intent(inout) :: rcoarseVector
!</output>

!</subroutine>

    call mlprj_performInterpolation (&
        rprjHierarchy%P_Rprojection(ifineLevel-rprjHierarchy%nlmin+1),&
        rcoarseVector,rfineVector,rprjHierarchy%rtempVector)

  end subroutine

end module
