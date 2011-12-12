  ! -*- mode: f90; -*-
  
  ! This file describes the interface for a coefficient function that is called
  ! during the evaluation of a finite element function.

  interface

  !<subroutine>


    subroutine ffunctionCoefficient (iderType, nelements, npointsPerElement, indof,&
                                     Idofs, Dpoints, Dbas, Dvalues, rcollection)

      use collection
      use fsystem

    !<description>
      ! This subroutine is called during the evaluation of a
      ! non-trivial FE-function. It has to evaluate the FE-function in
      ! the given cubature points Dpoints and store the result in
      ! Dvalues. In contrast to the standard feevaluation routines,
      ! this callback function can realise an arbitrary FE-function
      ! which is not available as precomputed vector.
    !</description>

    !<input>
      ! Type of function value to evaluate. One of the DER_xxxx constants,
      ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
      integer, intent(in) :: iderType

      ! Number of elements, where the FE-function must be computed.
      integer, intent(in) :: nelements
    
      ! Number of points per element, where the FE-function must be computed
      integer, intent(in) :: npointsPerElement

      ! Number of local DOF`s for trial space
      integer, intent(in) :: indof

      ! List of global DOF`s for trial space
      ! DIMENSION(1..indofs,1..nelements)
      integer, dimension(:,:), intent(in) :: Idofs

      ! List of points, where the FE-function must be computed.
      ! DIMENSION(1..ndim,1..npoints,1..nelements)
      real(DP), dimension(:,:,:), intent(in) :: Dpoints

      ! Precomputed array of basis functions and their derivatives
      ! DIMENSION(1..indof,1..MAXNBAS,1..npoints,1..nelements)
      real(DP), dimension(:,:,:,:), intent(in) :: Dbas
    !</input>

    !<inputoutput>
      ! Optional: A collection structure to provide additional
      ! information to the coefficient routine.
      type(t_collection), intent(inout), optional :: rcollection
    !</inputoutput>
  
    !<output>
      ! Evaluated FE-function for all given points on all given elements.
      ! DIMENSION(npointsPerElement,nelements)
      real(DP), dimension(:,:), intent(out) :: Dvalues
    !</output>
    
    !</subroutine>
  
    end subroutine
    
  end interface
