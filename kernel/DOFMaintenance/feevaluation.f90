!##############################################################################
!# ****************************************************************************
!# <name> feevaluation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines to evaluate a finite element function u
!# represented by a scalar vector rvector in one or multiple points in the
!# domain.
!# </purpose>
!##############################################################################

MODULE feevaluation

  USE fsystem
  USE linearsystemscalar
  
  IMPLICIT NONE

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE fevl_evaluate_sim (rvectorScalar, Dcoords, Djac, Ddetj, &
                  ieltyp, IdofsTrial, npoints,  nelements, Dpoints, iderType,&
                  Dvalues)
                                      
!<description>
  ! This routine allows to evaluate a finite element solution vector
  ! rvectorScalar simultaneously in multiple points on multiple elements in a 
  ! discretisation.
  ! The caller must provide all necessary information for the evaluation in the
  ! parameters.
  ! This routine is specialised to evaluate in multiple elements. For this
  ! purpose, the caller must make sure, that the same finite element type
  ! is used on all elements where to evaluate!
  ! So, evaluating 'simultaneously' on some $Q_1$ and some $P_1$ elements
  ! is not allowed e.g.. 
!</description>

!<input>
  ! The scalar solution vector that is to be evaluated.
  TYPE(t_vectorScalar), INTENT(IN)              :: rvectorScalar
  
  ! The FE function must be discretised with the same trial functions on all
  ! elements where it should be evaluated here. ieltyp defines the type
  ! of FE trial function that was used for the discretisation on those 
  ! elements that we are concerning here.
  INTEGER, INTENT(IN)                           :: ieltyp

  ! A list of the corner vertices of all elements in progress.
  ! array [1..NDIM2D,1..TRIA_MAXNVE2D,1..Number of elements] of double
  REAL(DP), DIMENSION(:,:,:), INTENT(IN)        :: Dcoords
  
  ! The Jacobian matrix of the mapping between the reference and each
  ! real element, for all points on all elements in progress.
  ! array [1..TRAFO_NJACENTRIES,1..npointsPerElement,1..Number of elements]
  REAL(DP), DIMENSION(:,:,:),INTENT(IN)         :: Djac
  
  ! The Jacobian determinant of the mapping of each point from the
  ! reference element to each real element in progress.
  ! array [1..npointsPerElement,1..Number of elements]
  REAL(DP), DIMENSION(:,:), INTENT(IN)          :: Ddetj
  
  ! An array accepting the DOF's on all elements in the trial space
  ! of the FE function.
  ! DIMENSION(\#local DOF's in trial space,nelements)
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(IN) :: IdofsTrial
  
  ! Number of points on every element where to evalate the function
  INTEGER, INTENT(IN) :: npoints
  
  ! Number of elements, the function is evaluated at
  INTEGER, INTENT(IN)  :: nelements
  
  ! Array with coordinates of the points where to evaluate.
  ! DIMENSION(NDIM2D,npoints,nelements).
  ! The coordinates are expected 
  ! - on the reference element, if ieltyp identifies a parametric element
  ! - on the real element, if ieltyp identifies a nonparametric element
  ! It's assumed that:
  !  Dpoints(1,.)=x-coordinates,
  !  Dpoints(2,.)=y-coordinates.
  ! furthermore:
  !  Dpoints(:,i,.) = Coordinates of point i
  ! furthermore:
  !  Dpoints(:,:,j) = Coordinates of all points on element j
  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: Dpoints

  ! Type of function value to evaluate. One of the DER_xxxx constants,
  ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
  INTEGER, INTENT(IN)                            :: iderType

!</input>

!<output>
  ! Values of the FE function at the points specified by Dpoints.
  ! DIMENSION(npoints,nelements).
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dvalues
!</output>

!</subroutine>

  ! local variables
  LOGICAL, DIMENSION(EL_MAXNDER) :: Bder
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: DbasTrial
  INTEGER :: indofTrial
  REAL(DP) :: dval
  INTEGER :: iel,ipoint,ibas
  REAL(DP), DIMENSION(:), POINTER :: p_Ddata
  REAL(SP), DIMENSION(:), POINTER :: p_Fdata
  
  ! What to evaluate?
  Bder = .FALSE.
  Bder(iderType) = .TRUE.
  
  ! Allocate memory for the basis function values
  indofTrial = elem_igetNDofLoc(ieltyp)
  ALLOCATE(DbasTrial(indofTrial,elem_getMaxDerivative(ieltyp),npoints,nelements))
  
  ! Evaluate the basis functions
  CALL elem_generic_sim (ieltyp, Dcoords, Djac, Ddetj, &
                         Bder, DbasTrial, npoints, nelements, Dpoints)  
  
  IF (rvectorScalar%cdataType .EQ. ST_DOUBLE) THEN
  
    ! Get the data array from the vector
    CALL lsyssc_getbase_double(rvectorScalar,p_Ddata)
    
    ! Now that we have the basis functions, we want to have the function values.
    ! We get them by multiplying the FE-coefficients with the values of the
    ! basis functions and summing up.
    DO iel=1,nelements
      DO ipoint = 1,npoints
        ! Calculate the value in the point
        dval = 0.0_DP
        DO ibas = 1,indofTrial
          dval = dval + &
                 p_Ddata(IdofsTrial(ibas,iel)) * DbasTrial(ibas,iderType,ipoint,iel)
        END DO
        ! Save the value in the point
        Dvalues(ipoint,iel) = dval
      END DO
    END DO
    
  ELSE IF (rvectorScalar%cdataType .EQ. ST_SINGLE) THEN
  
    ! Get the data array from the vector
    CALL lsyssc_getbase_single(rvectorScalar,p_Fdata)
    
    ! Now that we have the basis functions, we want to have the function values.
    ! We get them by multiplying the FE-coefficients with the values of the
    ! basis functions and summing up.
    DO iel=1,nelements
      DO ipoint = 1,npoints
        ! Calculate the value in the point
        dval = 0.0_DP
        DO ibas = 1,indofTrial
          dval = dval + &
                 p_Fdata(IdofsTrial(ibas,iel)) * DbasTrial(ibas,iderType,ipoint,iel)
        END DO
        ! Save the value in the point
        Dvalues(ipoint,iel) = dval
      END DO
    END DO
    
  ELSE
    PRINT *,'fevl_evaluate_sim: Unsupported vector precision'
    STOP
  END IF
  
  ! Release memory, finish
  DEALLOCATE(DbasTrial)

  END SUBROUTINE

END MODULE