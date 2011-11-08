!##############################################################################
!# ****************************************************************************
!# <name> adaptivecubature </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# This module provides routines for adaptive cubature. The very basic
!# implementation contains routines to support the summed cubature by
!# providing functions that automatically adapt a standard cubature
!# rule until the cubature error drops below a given tolerance.
!#
!# For some reference about adaptive cubature see e.g. here:
!#
!#  [Werner Vogt: Adaptive Verfahren zur numerischen Quadratur und Kubatur,
!#   Preprint, TU Ilmenau, 2006]
!#
!# The following subroutines can be found here:
!#
!# 1.) adcub_determineSummedCubature
!#     -> Analyses a function and determines an id for an adaptive cubature
!#        formula
!#
!# 2.) adcub_integrateFunction
!#     -> Integrate a function on a set of elements with a given cubature rule
!#
!# </purpose>
!##############################################################################

module adaptivecubature

  use basicgeometry
  use collection
  use cubature
  use element
  use elementpreprocessing
  use fsystem
  use genoutput
  use perfconfig
  use triangulation

  implicit none

  private
  
  public :: adcub_determineSummedCubature
  public :: adcub_integrateFunction

contains

  !****************************************************************************

!<subroutine>

  subroutine adcub_determineSummedCubature(ccubType,depsRel,depsAbs,rtriangulation,&
      Ielements,ffunctionRefSimple,rcollection,rperfconfig)

!<description>
  ! Determines an adaptive cubature formula which is accurate enough
  ! to integrate the function ffunctionReference on the elements Ielements
  ! of the mesh rtriangulation.
  !
  ! The cubature rule is determined in such a way, that the relative and
  ! absolute error is lower than depsRel / depsAbs.
!</description>

!<inputoutput>
  ! Cubature type identifier.
  ! On entry: Must be initialised with a basic cubature formula CCUB_xxxx
  ! that defines the underlying cubature rule (e.g. CUB_G2_2D for 2x2-gauss in 2D).
  ! On exit: Is replaced by a new cubature rule identifier that defines the
  ! adaptive cubature formula.
  integer(I32), intent(inout) :: ccubType
!</inputoutput>

!<input>
  ! Bound for the relative error in the integration.
  real(DP), intent(in) :: depsRel

  ! Bound for the absolute error in the integration.
  real(DP), intent(in) :: depsAbs

  ! Underlying mesh where ffunctionReference should be evaluated.
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! List of elements where the integration should be carried out.
  ! All elements must be of the same type!!! (i.e. all must be triangles,
  ! quads, hexas or similar); it is not allowed to mix e.g. triangles and
  ! quads in this list.
  integer, dimension(:), intent(in) :: Ielements
  
  ! Reference function which should be tested for integration.
  include '../Postprocessing/intf_functionScSimple.inc'
  
  ! OPTIONAL: Collection structure to be passed to ffunctionReference.
  type(t_collection), intent(inout), optional :: rcollection

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>

!</subroutine>

    real(DP) :: dvalue, dvalue2, derror, dvalueinit
    integer(I32) :: ccurrentcubtype
    integer :: icubref
    
    ! The inital cubature rule is without refinement.
    icubref = 0
    ccurrentcubtype = cub_getStdCubType(ccubType)
    
    ! Determine the initial value of the integral.
    call adcub_integrateFunction(dvalueinit,ccurrentcubtype,rtriangulation,&
        Ielements,ffunctionRefSimple,rcollection,rperfconfig)
    dvalue = dvalueinit
        
    ! Now refine the integration until the error is small enough.
    do
      icubref = icubref + 1
      ccurrentcubtype = cub_getSummedCubType(ccubType,icubref)
      dvalue2 = dvalue
      call adcub_integrateFunction(dvalue,ccurrentcubtype,rtriangulation,&
          Ielements,ffunctionRefSimple,rcollection,rperfconfig)
      
      derror = (dvalue-dvalue2) / (2**icubref - 1)
      
      ! Stop if the error is small enough. Ok, this formula is actually
      ! slightly wrong, since we stop, if the error on the 'fine' mesh
      ! is small enough while we return the cubature formula of the coarse
      ! mesh, but anyway, this is only a rough approximation.
      ! If we do not do this, there would be always at least one refinement
      ! what prevents the unrefined cubature rule to be used.
      !
      ! And: At most 5 refinements!
      if ((derror .lt. depsAbs) .and. &
          (derror*dvalueinit .lt. depsRel) .or. &
          (icubref .gt. 5)) then
        ccubType = cub_getSummedCubType(ccubType,icubref-1)
        exit
      end if
    
    end do

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine adcub_integrateFunction(dvalue,ccubType,rtriangulation,Ielements,&
      ffunctionRefSimple,rcollection,rperfconfig)

!<description>
  ! Integrates a scalar function ffunctionReference on a set of elements
  ! Ielements of a given mesh rtriangulation with a cubature formula ccubType.
!</description>

!<input>
  ! Cubature type identifier of the cubature formula to use for integration.
  integer(I32), intent(in) :: ccubType

  ! Underlying mesh where ffunctionReference should be evaluated.
  type(t_triangulation), intent(in) :: rtriangulation
  
  ! List of elements where the integration should be carried out.
  ! All elements must be of the same type!!! (i.e. all must be triangles,
  ! quads, hexas or similar); it is not allowed to mix e.g. triangles and
  ! quads in this list.
  integer, dimension(:), intent(in), target :: Ielements
  
  ! Reference function which should be tested for integration.
  include '../Postprocessing/intf_functionScSimple.inc'
  
  ! OPTIONAL: Collection structure to be passed to ffunctionReference.
  type(t_collection), intent(inout), optional :: rcollection

  ! OPTIONAL: local performance configuration. If not given, the
  ! global performance configuration is used.
  type(t_perfconfig), intent(in), optional :: rperfconfig
!</input>

!<output>
  ! The value of the integral.
  real(DP), intent(out) :: dvalue
!</output>

!</subroutine>

    integer :: IELset, IELmax, iel, NEL, idim, icubp
    integer :: ncubp,nelementsPerBlock,nmaxelements
    integer(I32) :: cevaluationtag
    real(DP), dimension(:,:), allocatable :: DpointsRef
    real(DP), dimension(:), allocatable :: DomegaRef
    type(t_evalElementSet) :: revalElementSet
    integer(I32) :: ctrafoType
    real(DP), dimension(:,:), allocatable :: Dcoefficients
    real(DP), dimension(:,:), pointer :: p_Ddetj
    real(DP) :: om

    ! Get the cubature points on the reference element
    ncubp = cub_igetNumPts(ccubType)
    idim = cub_igetCoordDim(ccubType)
    allocate (DpointsRef(idim,ncubp))
    allocate (DomegaRef(ncubp))
    call cub_getCubature(ccubType, DpointsRef, DomegaRef)

    ! Number of elements
    NEL = size(Ielements)
    
    ! Determine the number of elements in each set in such a way,
    ! that the element blocks are not too large and that the
    ! number of simultaneous treated cubature points is not too large.
    ! We allow the treatment of at most 16x16x16 cubature points per element
    ! on a set of 1000 elements, resulting in at most 16x16x16x1000=4.096.000
    ! cubature points. So reduce the number of simultaneouly treated
    ! elements if necessary -- but treat at least one element.
    nmaxelements = max(1,min(ncubp*1000,16*16*16*1000)/ncubp)
    nelementsPerBlock = min(nmaxelements,NEL)
    
    ! Allocate memory for the coefficients
    allocate(Dcoefficients(ncubp,nelementsPerBlock))

    ! Initialisation of the element set.
    call elprep_init(revalElementSet)
    cevaluationtag = EL_EVLTAG_REALPOINTS + EL_EVLTAG_JAC + EL_EVLTAG_DETJ + &
        EL_EVLTAG_REFPOINTS + EL_EVLTAG_COORDS
    
    ! Get the transformation type from the basic P1/Q1 element in the
    ! current dimension
    ! Get the shape of the cubature id
    select case(cub_igetShape(ccubType))
    case(BGEOM_SHAPE_LINE)
      ! 1D line formula -> reference coordinates
      ctrafoType = elem_igetTrafoType(EL_P1_1D)
    case(BGEOM_SHAPE_TRIA)
      ! 2D triangle formula -> barycentric coordinates
      ctrafoType = elem_igetTrafoType(EL_P1_2D)
    case(BGEOM_SHAPE_QUAD)
      ! 2D quadrilateral formula -> reference coordinates
      ctrafoType = elem_igetTrafoType(EL_Q1_2D)
    case default
      ! unknown formula
      call output_line ('Unsupported element.', &
                        OU_CLASS_ERROR,OU_MODE_STD,'adcub_integrateFunction')
      call sys_halt()
    end select
    
    dvalue = 0.0_DP
    
    ! Loop throgh the elements in sets.
    do IELset = 1,NEL,nmaxelements
    
      ! Maximum element of the set
      IELmax = min(NEL, IELset-1+nmaxelements)
    
      ! Prepare the element set to calculate information of the transformation.
      call elprep_prepareSetForEvaluation (revalElementSet,&
          cevaluationTag, rtriangulation, Ielements(IELset:IELmax), &
          ctrafoType, DpointsRef(:,:), rperfconfig=rperfconfig)
      p_Ddetj => revalElementSet%p_Ddetj

      ! In the next loop, we do not have to evaluate the coordinates
      ! on the reference elements anymore.
      cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))

      ! Calculate the values in the cubature points.
      call ffunctionRefSimple (IELmax-IELset+1,ncubp,Ielements(IELset:IELmax),&
          revalElementSet%p_DpointsReal,Dcoefficients(:,1:IELmax-IELset+1),rcollection,&
          revalElementSet%p_DpointsRef,revalElementSet%p_Djac,revalElementSet%p_Ddetj)
      
      ! Sum up to the contribution of the integral.
      do IEL = 1, IELmax-IELset+1
        
        ! Loop over all cubature points on the current element
        do icubp = 1, ncubp
          
          ! calculate the current weighting factor in the cubature formula
          ! in that cubature point.
          !
          ! Take the absolut value of the determinant of the mapping.
          ! In 2D, the determinant is always positive, whereas in 3D,
          ! the determinant might be negative -- that is normal!
          
          OM = DomegaRef(icubp)*abs(p_Ddetj(icubp,IEL))
          
          dvalue = dvalue + OM * Dcoefficients(icubp,IEL)
          
        end do ! ICUBP
        
      end do ! iel
    
    end do
    
    ! Release memory, finish
    call elprep_releaseElementSet(revalElementSet)
    deallocate(DomegaRef)
    deallocate(DpointsRef)

  end subroutine

end module
