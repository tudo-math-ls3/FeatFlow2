!##############################################################################
!# ****************************************************************************
!# <name> poissoncallback </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains callback functions for the poisson problem that are
!# used during the matrix/vector assembly for specifying analytical data.
!# There are three callback functions involved, which may be called depending
!# on the situation. All of them correspond to a specific interface for
!# callback functions, defined in 'intf_xxxx.inc' files.
!#
!# 1.) coeff_Laplace
!#     -> Returns the coefficients for the Laplace matrix. This routine is
!#        only used if the problem to calculate has nonconstant coefficients!
!#        Otherwise the routine is dead.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientMatrixSc.inc'
!#
!# 2.) coeff_RHS
!#     -> Returns analytical values for the right hand side of the Laplace
!#        equation.
!#     -> Corresponds to the interface defined in the file
!#        'intf_coefficientVectorSc.inc'
!#
!# 3.) getBoundaryValues
!#     -> Returns analitical values on the (Dirichlet) boundary of the
!#        problem to solve.
!#     -> Corresponds to the interface defined in the file
!#        'intf_bcassembly.inc'
!#
!# </purpose>
!##############################################################################

MODULE poisson_callback

  USE fsystem
  USE storage
  USE linearsolver
  USE boundary
  USE bilinearformevaluation
  USE linearformevaluation
  USE cubature
  USE matrixfilters
  USE vectorfilters
  USE bcassembly
  
  IMPLICIT NONE

CONTAINS

  ! ***************************************************************************
  !<subroutine>

  SUBROUTINE coeff_Laplace (rdiscretisation,ielementDistribution, rform, &
                ielementStartIdx,nelements,npointsPerElement,Ielements,Dcoords, &
                DcubPtsRef,DcubPtsReal,IdofsTrial,IdofsTest,Djac,Ddetj,p_rcollection, &
                Dcoefficients)
  
!<description>
  ! This subroutine is called during the matrix assembly. It has to compute
  ! the coefficients in front of the terms of the bilinear form.
  !
  ! The routine accepts a set of elements and a set of points on these
  ! elements (cubature points), in reference as well as in real coordinates.
  ! According to the terms in the bilinear form, the routine has to compute
  ! simultaneously for all these points and all the terms in the bilinear form
  ! the corresponding coefficients in front of the terms.
!</description>
  
!<input>
  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  TYPE(t_spatialDiscretisation), INTENT(IN)                   :: rdiscretisation
  
  ! The currently active element distribution in the discretisation.
  ! Allows the routine to get the currently active element type for
  ! trial and test functions.
  INTEGER, INTENT(IN)                                         :: ielementDistribution
  
  ! The bilinear form which is currently to be evaluated:
  TYPE(t_bilinearForm), INTENT(IN)                            :: rform
  
  ! Start index of the current element block Ielements in the current element 
  ! distribution ielementDistribution. If this is =1, the routine is called the 
  ! first time for the current element distribution.
  INTEGER(I32), INTENT(IN)                                    :: ielementStartIdx

  ! Number of elements, where the coefficients must be computed.
  ! This is always a part of the element distribution.
  INTEGER, INTENT(IN)                                         :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  INTEGER, INTENT(IN)                                         :: npointsPerElement
  
  ! A list of elements of length nelements where coefficients must
  ! be computed by this routine.
  INTEGER(I32), DIMENSION(nelements), INTENT(IN)              :: Ielements
  
  ! A list of the corner vertices of all elements in Ielements
  REAL(DP), DIMENSION(2,TRIA_MAXNVE2D,nelements), INTENT(IN)  :: Dcoords
  
  ! A list of points in coordinates on the reference element.
  ! Each set of points corresponds to the corresponding element
  ! in Ielements
  REAL(DP), DIMENSION(NDIM2D,npointsPerElement,nelements), INTENT(IN) :: DcubPtsRef

  ! A list of points, corresponding to DcubPtsRef, in real coordinates.
  ! Each set of points corresponds to the corresponding element
  ! in Ielements.
  REAL(DP), DIMENSION(NDIM2D,npointsPerElement,nelements), INTENT(IN) :: DcubPtsReal
  
  ! An array accepting the DOF's on all elements trial in the trial space.
  ! DIMENSION(#local DOF's in trial space,nelements)
  INTEGER(PREC_DOFIDX), DIMENSION(:,:) :: IdofsTrial
  
  ! An array accepting the DOF's on all elements trial in the trial space.
  ! DIMENSION(#local DOF's in test space,nelements)
  INTEGER(PREC_DOFIDX), DIMENSION(:,:) :: IdofsTest
  
  ! The Jacobian matrix of the mapping between the reference and each
  ! real element, for all points on all elements.
  REAL(DP), DIMENSION(TRAFO_NJACENTRIES,npointsPerElement,nelements), INTENT(IN) :: Djac
  
  ! The Jacobian determinant of the mapping of each point from the
  ! reference element to each real element
  REAL(DP), DIMENSION(npointsPerElement,nelements), INTENT(IN) :: Ddetj
  
  ! A pointer to a collection structure to provide additional 
  ! information to the coefficient routine. May point to NULL() if not defined.
  TYPE(t_collection), POINTER                   :: p_rcollection
  
!</input>

!<output>
  ! A list of all coefficients in front of all terms in the bilinear form -
  ! for all given points on all given elements.
  !   DIMENSION(itermCount,npointsPerElement,nelements)
  ! with itermCount the number of terms in the bilinear form.
  REAL(DP), DIMENSION(:,:,:), INTENT(OUT)                      :: Dcoefficients
!</output>
  
!</subroutine>

  Dcoefficients = 1.0_DP

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE coeff_RHS (rdiscretisation,ielementDistribution, rform, &
                ielementStartIdx,nelements,npointsPerElement,Ielements,Dcoords, &
                DcubPtsRef,DcubPtsReal,IdofsTest,Djac,Ddetj,p_rcollection, &
                Dcoefficients)
  
!<description>
  ! This subroutine is called during the matrix assembly. It has to compute
  ! the coefficients in front of the terms of the bilinear form.
  !
  ! The routine accepts a set of elements and a set of points on these
  ! elements (cubature points), in reference as well as in real coordinates.
  ! According to the terms in the bilinear form, the routine has to compute
  ! simultaneously for all these points and all the terms in the bilinear form
  ! the corresponding coefficients in front of the terms.
!</description>
  
!<input>
  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  TYPE(t_spatialDiscretisation), INTENT(IN)                   :: rdiscretisation
  
  ! The currently active element distribution in the discretisation.
  ! Allows the routine to get the currently active element type for
  ! trial and test functions.
  INTEGER, INTENT(IN)                                         :: ielementDistribution
  
  ! The linear form which is currently to be evaluated:
  TYPE(t_linearForm), INTENT(IN)                              :: rform
  
  ! Start index of the current element block Ielements in the current element 
  ! distribution ielementDistribution. If this is =1, the routine is called the 
  ! first time for the current element distribution.
  INTEGER(I32), INTENT(IN)                                    :: ielementStartIdx

  ! Number of elements, where the coefficients must be computed.
  ! This is always a part of the element distribution.
  INTEGER, INTENT(IN)                                         :: nelements
  
  ! Number of points per element, where the coefficients must be computed
  INTEGER, INTENT(IN)                                         :: npointsPerElement
  
  ! A list of elements of length nelements where coefficients must
  ! be computed by this routine.
  INTEGER(I32), DIMENSION(nelements), INTENT(IN)              :: Ielements
  
  ! A list of the corner vertices of all elements in Ielements
  REAL(DP), DIMENSION(2,TRIA_MAXNVE2D,nelements), INTENT(IN)  :: Dcoords
  
  ! A list of points in coordinates on the reference element.
  ! Each set of points corresponds to the corresponding element
  ! in Ielements
  REAL(DP), DIMENSION(NDIM2D,npointsPerElement,nelements), INTENT(IN) :: DcubPtsRef

  ! A list of points, corresponding to DcubPtsRef, in real coordinates.
  ! Each set of points corresponds to the corresponding element
  ! in Ielements.
  REAL(DP), DIMENSION(NDIM2D,npointsPerElement,nelements), INTENT(IN) :: DcubPtsReal
  
  ! An array accepting the DOF's on all elements trial in the trial space.
  ! DIMENSION(#local DOF's in test space,nelements)
  INTEGER(PREC_DOFIDX), DIMENSION(:,:) :: IdofsTest
  
  ! The Jacobian matrix of the mapping between the reference and each
  ! real element, for all points on all elements.
  REAL(DP), DIMENSION(TRAFO_NJACENTRIES,npointsPerElement,nelements), INTENT(IN) :: Djac
  
  ! The Jacobian determinant of the mapping of each point from the
  ! reference element to each real element
  REAL(DP), DIMENSION(npointsPerElement,nelements), INTENT(IN) :: Ddetj
  
  ! A pointer to a collection structure to provide additional 
  ! information to the coefficient routine. May point to NULL() if not defined.
  TYPE(t_collection), POINTER                   :: p_rcollection
  
!</input>

!<output>
  ! A list of all coefficients in front of all terms in the bilinear form -
  ! for all given points on all given elements.
  !   DIMENSION(itermCount,npointsPerElement,nelements)
  ! with itermCount the number of terms in the bilinear form.
  REAL(DP), DIMENSION(:,:,:), INTENT(OUT)                      :: Dcoefficients
!</output>
  
!</subroutine>

  Dcoefficients = 1.0_DP

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE getBoundaryValues (rdiscretisation,rbcRegion,ielement, &
                                  cinfoNeeded,iwhere,dwhere, p_rcollection, Dvalues)
  
  USE collection
  USE spatialdiscretisation
  USE discretebc
  
!<description>
  ! This subroutine is called during the discretisation of boundary
  ! conditions. It calculates a special quantity on the boundary, which is
  ! then used by the discretisation routines to generate a discrete
  ! 'snapshot' of the (actually analytic) boundary conditions.
!</description>
  
!<input>
  ! The discretisation structure that defines the basic shape of the
  ! triangulation with references to the underlying triangulation,
  ! analytic boundary boundary description etc.
  TYPE(t_spatialDiscretisation), INTENT(IN)                   :: rdiscretisation
  
  ! Boundary condition region that is currently being processed.
  ! (This e.g. defines the type of boundary conditions that are
  !  currently being calculated, as well as information about the current
  !  boundary segment 'where we are at the moment'.)
  TYPE(t_bcRegion), INTENT(IN)                                :: rbcRegion
  
  
  ! The element number on the boundary which is currently being processed
  INTEGER(I32), INTENT(IN)                                    :: ielement
  
  ! The type of information, the routine should calculate. One of the
  ! DISCBC_NEEDxxxx constants. Depending on the constant, the routine has
  ! to return one or multiple information value in the result array.
  INTEGER, INTENT(IN)                                         :: cinfoNeeded
  
  ! A reference to a geometric object where information should be computed.
  ! cinfoNeeded=DISCBC_NEEDFUNC : 
  !   iwhere = number of the point in the triangulation or
  !          = 0, if only the parameter value of the point is known; this
  !               can be found in dwhere,
  ! cinfoNeeded=DISCBC_NEEDDERIV : 
  !   iwhere = number of the point in the triangulation or
  !          = 0, if only the parameter value of the point is known; this
  !               can be found in dwhere,
  ! cinfoNeeded=DISCBC_NEEDINTMEAN : 
  !   iwhere = number of the edge where the value integral mean value
  !            should be computed
  INTEGER, INTENT(IN)                                         :: iwhere

  ! A reference to a geometric object where information should be computed.
  ! cinfoNeeded=DISCBC_NEEDFUNC : 
  !   dwhere = parameter value of the point where the value should be computed,
  ! cinfoNeeded=DISCBC_NEEDDERIV : 
  !   dwhere = parameter value of the point where the value should be computed,
  ! cinfoNeeded=DISCBC_NEEDINTMEAN : 
  !   dwhere = 0 (not used)
  REAL(DP), INTENT(IN)                                        :: dwhere
    
  ! A pointer to a collection structure to provide additional 
  ! information to the coefficient routine. May point to NULL() if not defined.
  TYPE(t_collection), POINTER                  :: p_rcollection

!</input>

!<output>
  ! This array receives the calculated information. If the caller
  ! only needs one value, the computed quantity is put into Dvalues(1). 
  ! If multiple values are needed, they are collected here (e.g. for 
  ! DISCBC_NEEDDERIV: Dvalues(1)=x-derivative, Dvalues(2)=y-derivative,...)
  REAL(DP), DIMENSION(:), INTENT(OUT)                         :: Dvalues
!</output>
  
!</subroutine>

  ! Return zero Dirichlet boundary values for all situations.
  Dvalues(1) = 0.0_DP

  END SUBROUTINE

END MODULE
