!<types>

!<typeblock>

  ! Structure which collects information for function evaluation routines.
  ! Remarks:
  ! a) If the structure describes a number of points where to evaluate 
  !    without an underlying mesh, there is
  !        nelements = 1
  !        npointPerElement = #points
  !        p_Dpoints = coordinates 
  !    and p_rtriangulation => null().
  ! b) If the structure describes only one point where to evaluate, there is
  !        nelements = 1
  !        npointPerElement = 1
  !        p_Dpoints = coordinates.
  type t_paramsSimEval

    ! <!-- ############################################### -->
    ! <!-- Mandatory information, must always be available -->
    ! <!-- ############################################### -->

    ! Number of elements.
    integer :: nelements = 0

    ! Number of points per element.
    integer :: npointsPerElement = 0

    ! Array of all points on all the elements where coefficients are needed. 
    ! DIMENSION(#dimensions,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), pointer :: p_Dpoints => null()

    ! <!-- ###################################### -->
    ! <!-- Optional information, mostly available -->
    ! <!-- ###################################### -->

    ! Pointer to the underlying boundary structure.
    ! If not available, this points to null().
    type(t_boundary), pointer :: p_rboundary => null()

    ! Pointer to the underlying triangulation
    ! If not available, this points to null().
    type(t_triangulation), pointer :: p_rtriangulation => null()
    
    ! List of elements.
    ! If not available, this points to null().
    integer, dimension(:), pointer :: p_Ielements => null()

    ! Coordinates of the evaluation points in the reference
    ! coordinate system. 
    ! If this information is not available, this parameter points to null().
    ! DIMENSION(#dimensions,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), pointer :: p_DpointsRef => null()

    ! The Jacobian matrix of the mapping between the reference and each
    ! real element, for all points on all elements in progress.
    ! If this information is not available, this parameter points to null().
    ! DIMENSION(#dimensions*#dimensions,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), pointer :: p_Djac => null()
    
    ! The Jacobian determinant of the mapping of each point from the
    ! reference element to each real element in progress.
    ! If this information is not available, this parameter points to null().
    ! DIMENSION(npointsPerElement,nelements)
    real(DP), dimension(:,:), pointer :: p_Ddetj => null()
    
    ! <!-- ################################ -->
    ! <!-- Evaluation on the boundary in 2D -->
    ! <!-- ################################ -->
    
    ! Number of the current boundary component, or =0 if not available.
    ! All elements in p_Ielements refer to this boundary component.
    integer :: ibct = 0
    
    ! Parameter values of the points on the boundary
    ! If this information is not available, this parameter points to null().
    real(dp), dimension(:), pointer :: p_DpointsPar => null()
    
    ! Boundary region that is currently being processed.
    ! If this information is not available, this parameter points to null().
    type(t_boundaryRegion), pointer :: rboundaryRegion => null()
    
    ! Mesh region that is currently being processed.
    ! If this information is not available, this parameter points to null().
    type(t_meshRegion), pointer :: p_rmeshRegion => null()
    
    ! <!-- ########################## -->
    ! <!-- Evaluation of linear forms -->
    ! <!-- ########################## -->
    
    ! The bilinear form which is currently being evaluated.
    ! If this information is not available, this parameter points to null().
    type(t_linearForm), pointer :: p_rlinearForm => null()

    ! <!-- ############################ -->
    ! <!-- Evaluation of bilinear forms -->
    ! <!-- ############################ -->

    ! The bilinear form which is currently being evaluated.
    ! If this information is not available, this parameter points to null().
    type(t_bilinearForm), pointer :: p_rbilinearForm => null()

    ! <!-- ############################# -->
    ! <!-- Evaluation of trilinear forms -->
    ! <!-- ############################# -->

    ! The trilinear form which is currently being evaluated.
    ! If this information is not available, this parameter points to null().
    type(t_trilinearForm), pointer :: p_rtrilinearForm => null()
  end type

!</typeblock>

!</types>

!<constants>
!<constantblock description="Task identifiers for callback function evaluation">
  ! Do nothing - this is just a dummy.
  integer, parameter :: FUNC_TASK_NONE                  = 0

  ! Evaluate the function in a given set of points.
  ! Results must be saved to Dvalues(:,:,1).
  integer, parameter :: FUNC_TASK_EVAL_FUNC             = 1

  ! Evaluate one component of a derivative.
  ! The variable itag is a derivative quantifier (DER_DERIV_xxxx) that
  ! defines the component that should be evaluated.
  !   Dvalues(:,:,1)=derivatives identified by itag
  integer, parameter :: FUNC_TASK_EVAL_DERIV            = 2

  ! Evaluate the full 1st order derivative of the function.
  ! Results must be saved to:
  !   Dvalues(:,:,1)=x-derivatives,
  !   Dvalues(:,:,2)=y-derivatives (if exists),
  !   Dvalues(:,:,3)=z-derivatives (if exists).
  integer, parameter :: FUNC_TASK_EVAL_DERIV1FULL       = 3

  ! Evaluate the function in a given set of points on the boundary.
  ! Results must be saved to Dvalues(:,:,1).
  integer, parameter :: FUNC_TASK_EVAL_FUNC_BD          = 4

  ! Evaluate one component of a derivativeon the boundary.
  ! The variable itag is a derivative quantifier (DER_DERIV_xxxx) that
  ! defines the component that should be evaluated.
  !   Dvalues(:,:,1)=derivatives identified by itag
  integer, parameter :: FUNC_TASK_EVAL_DERIV_BD         = 5

  ! Evaluate the 1st derivative of the function on the boundary.
  ! Results must be saved to:
  !   Dvalues(:,:,1)=x-derivatives,
  !   Dvalues(:,:,2)=y-derivatives (if exists),
  !   Dvalues(:,:,3)=z-derivatives (if exists).
  integer, parameter :: FUNC_TASK_EVAL_DERIV1_BD        = 6

  ! Evaluate the coefficients in front of the terms of a linear form.
  ! Results must be saved to Dvalues(:,:,1..nterms) with nterms the number
  ! of terms in the linear form.
  integer, parameter :: FUNC_TASK_EVAL_LINFCOEFF        = 7

  ! Evaluate the coefficients in front of the terms of a bilinear form.
  ! Results must be saved to Dvalues(:,:,1..nterms) with nterms the number
  ! of terms in the bilinear form.
  integer, parameter :: FUNC_TASK_EVAL_BILFCOEFF        = 8

  ! Evaluate the coefficients in front of the terms of a trilinear form.
  ! Results must be saved to Dvalues(:,:,1..nterms) with nterms the number
  ! of terms in the trilinear form.
  integer, parameter :: FUNC_TASK_EVAL_TRILFCOEFF       = 9
!</constantblock>
!</constants>

  interface
  
  !<subroutine>

    subroutine ffunctionSimEval (ctask,itag,rparams,Dvalues,rcollection)
    
  !<description>
    ! This interface defines a general implementation of a callback routine
    ! used for the evaluation of functions in the kernel. It is used for matrix
    ! and vector assembly as well as for the calculation of itegrals.
    !
    ! The purpose of this routine is to calculate the values of a (usually
    ! analytically given) function in a number of points on a set of elements.
    ! The stricture rparams defines all parameters where to evaluate.
    ! The points where to evaluate can e.g. be accessed by:
    !    dx = rparams%p_Dpoints(1,ipoint,ielement) 
    !    dy = rparams%p_Dpoints(2,ipoint,ielement) 
    ! -> x- and y-coordinates of the point ipoint on element ielement.
    ! Depending on the purpose for which this routine is used,
    ! the parameters in rparams may provide additional information
    ! (like a triangulation, boundary, element list etc.).
    !
    ! The values that are calculated in this routine have to be written
    ! to the Dvalues(:,:,:) array. The last dimension of Dvalues defines
    ! here the number of a 'term'; for standard callback functions, 1 has
    ! to be used as an index here, e.g.
    !   Dvalues(:,:,1) = 16.0_DP * rparams%p_Dpoints(2,:,:) * (1.0_DP-rparams%p_Dpoints(2,:,:))
    ! returns a standard parabolic profile.
    !
    ! For more complicated function evaluations like the evaluation of a
    ! bilinear form, this routine may have to calculate multiple terms,
    ! like the coefficients of multiple additive integral terms; In this
    ! case, iterm in Dvalues(:,:,iterm) allows to access each term to
    ! calculate.
  !</description>
    
    use fsystem
    use collection

  !<input>
    ! A task identifier that specifies in which way the function should return
    ! its values. One of the FUNC_TASK_xxxx constants.
    ! The standard value is e.g. FUNC_TASK_EVAL_FUNC which tells the routine just to
    ! return function values.
    integer, intent(in) :: ctask
  !</input>

  !<inputoutput>
    ! A tag variable. The exact meaning depends on the assembly routine.
    ! If not used, 0 is passed here.
    integer, intent(inout) :: itag

    ! Parameter structure defining where to evaluate
    type(t_paramsSimEval), intent(inout) :: rparams

    ! Optional: A collection structure to provide additional 
    ! information to the coefficient routine. 
    type(t_collection), intent(inout), optional :: rcollection
  !</inputoutput>
  
  !<output>
    ! Result values. #terms varies, depending on ctask.
    ! DIMENSION(#points per element, #elements, #terms)
    real(DP), dimension(:,:,:), intent(out) :: Dvalues
  !</output>
    
  !</subroutine>
  
    end subroutine
    
  end interface

