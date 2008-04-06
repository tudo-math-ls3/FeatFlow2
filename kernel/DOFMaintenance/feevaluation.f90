!##############################################################################
!# ****************************************************************************
!# <name> feevaluation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines to evaluate finite element functions
!# in given points.
!#
!# The module contains the following routines:
!#
!# 1.) fevl_evaluate
!#     -> Evaluate a FE function in an arbitrary set of points.
!#        Determines automatically the elements that contain the points.
!#
!# 3.) fevl_evaluate_mult
!#     -> Evaluate a FE function simultaneously in multiple points on
!#        one element. The element containing the point must be given.
!#
!# 2.) fevl_evaluate_sim
!#     -> Evaluate a FE function simultaneously in multiple points on
!#        multiple elements. Elements containing the points must be given.
!#
!# </purpose>
!##############################################################################

MODULE feevaluation

  USE fsystem
  USE linearsystemscalar
  USE triasearch
  USE element
  USE elementpreprocessing
  
  IMPLICIT NONE

!<constants>

!<constantblock description="Constants for the cnonmeshPoints parameter.">

  ! Points outside of the domain not allowed.
  INTEGER, PARAMETER :: FEVL_NONMESHPTS_NONE   = 0

  ! Points outside of the domain allowed, try to find an element nearby to evaluate.
  INTEGER, PARAMETER :: FEVL_NONMESHPTS_NEARBY = 1

  ! Points outside of the domain allowed, assume 0 as value there.
  INTEGER, PARAMETER :: FEVL_NONMESHPTS_ZERO   = 2

!</constantblock>

!</constants>

  ! There are two functions fevl_evaluate_mult which do the same -- with
  ! different calling conventions and different complexities.
  INTERFACE fevl_evaluate_mult
    MODULE PROCEDURE fevl_evaluate_mult1
    MODULE PROCEDURE fevl_evaluate_mult2
  END INTERFACE

  ! There are two functions fevl_evaluate_sim which do the same -- with
  ! different calling conventions and different complexities.
  INTERFACE fevl_evaluate_sim
    MODULE PROCEDURE fevl_evaluate_sim1
    MODULE PROCEDURE fevl_evaluate_sim2
    MODULE PROCEDURE fevl_evaluate_sim3
  END INTERFACE

CONTAINS

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE fevl_evaluate (iderType, Dvalues, rvectorScalar, Dpoints, &
      Ielements, IelementsHint, cnonmeshPoints)
                                      
!<description>
  ! This is the most general (and completely slowest) finite element evaluation
  ! routine. It allows to evaluate a general (scalar) FE function specified
  ! by rvectorScalar in a set of points Dpoints. The values of the
  ! FE function are written into Dvalues.
!</description>

!<input>
  ! Type of function value to evaluate. One of the DER_xxxx constants,
  ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
  INTEGER, INTENT(IN)                            :: iderType

  ! The scalar solution vector that is to be evaluated.
  TYPE(t_vectorScalar), INTENT(IN)              :: rvectorScalar
  
  ! A list of points where to evaluate.
  ! DIMENSION(1..ndim,1..npoints)
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: Dpoints
  
  ! OPTIONAL: A list of elements containing the points Dpoints.
  ! If this is not specified, the element numbers containing the points
  ! are determined automatically.
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: Ielements

  ! OPTIONAL: A list of elements that are near the points in Dpoints.
  ! This gives only a hint where to start searching for the actual elements
  ! containing the points. This is ignored if Ielements is specified!
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:), INTENT(IN), OPTIONAL :: IelementsHint
  
  ! OPTIONAL: A FEVL_NONMESHPTS_xxxx constant that defines what happens
  ! if a point is located outside of the domain. May happen e.g. in
  ! nonconvex domains. FEVL_NONMESHPTS_NONE is the default 
  ! parameter if cnonmeshPoints is not specified. 
  INTEGER, INTENT(IN), OPTIONAL :: cnonmeshPoints
  
!</input>

!<output>
  ! Values of the FE function at the points specified by Dpoints.
  REAL(DP), DIMENSION(:), INTENT(OUT) :: Dvalues
!</output>

!</subroutine>

    ! local variables
    INTEGER :: cnonmesh
    INTEGER :: ipoint,ieltype,indof,nve,ibas
    INTEGER(PREC_ELEMENTIDX) :: iel
    INTEGER(I32), DIMENSION(:), POINTER :: p_IelementDistr
    LOGICAL, DIMENSION(EL_MAXNDER) :: Bder
    REAL(DP) :: dval
    
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    REAL(SP), DIMENSION(:), POINTER :: p_Fdata
    
    ! Transformation
    INTEGER(I32) :: ctrafoType
    REAL(DP), DIMENSION(TRAFO_MAXDIMREFCOORD) :: DparPoint
    
    ! Values of basis functions and DOF's
    REAL(DP), DIMENSION(EL_MAXNBAS,EL_MAXNDER) :: Dbas
    INTEGER(PREC_DOFIDX), DIMENSION(EL_MAXNBAS) :: Idofs
    
    ! List of element distributions in the discretisation structure
    TYPE(t_elementDistribution), DIMENSION(:), POINTER :: p_RelementDistribution

    ! Evaluation structure and tag
    TYPE(t_evalElement) :: revalElement
    INTEGER(I32) :: cevaluationTag

    ! Ok, slow but general.
    
    p_RelementDistribution => rvectorScalar%p_rspatialDiscretisation%RelementDistribution
    
    ! For uniform discretisations, we get the element type in advance...
    IF (rvectorScalar%p_rspatialDiscretisation%ccomplexity .EQ. SPDISC_UNIFORM) THEN
      
      ! Element type
      ieltype = rvectorScalar%p_rspatialDiscretisation%&
          RelementDistribution(1)%itrialElement

      ! Get the number of local DOF's for trial and test functions
      indof = elem_igetNDofLoc(ieltype)
      
      ! Number of vertices on the element
      nve = elem_igetNVE(ieltype)
      
      ! Type of transformation from/to the reference element
      ctrafoType = elem_igetTrafoType(ieltype)
      
      ! Get the element evaluation tag; necessary for the preparation of the element
      cevaluationTag = elem_getEvaluationTag(ieltype)
      
      NULLIFY(p_IelementDistr)
    ELSE
      CALL storage_getbase_int (&
          rvectorScalar%p_rspatialDiscretisation%h_IelementDistr,p_IelementDistr)
    END IF
    
    ! Get the data vector
    SELECT CASE (rvectorScalar%cdataType)
    CASE (ST_DOUBLE) 
      CALL lsyssc_getbase_double(rvectorScalar,p_Ddata)
    CASE (ST_SINGLE)
      CALL lsyssc_getbase_single(rvectorScalar,p_Fdata)
    CASE DEFAULT
      CALL output_line ('Unsupported vector precision!',&
          OU_CLASS_ERROR,OU_MODE_STD,'fevl_evaluate')
      CALL sys_halt()
    END SELECT
    
    ! What to evaluate?
    Bder = .FALSE.
    Bder(iderType) = .TRUE.
    
    cnonmesh = FEVL_NONMESHPTS_NONE
    IF (PRESENT(cnonmeshPoints)) cnonmesh = cnonmeshPoints
    
    ! We loop over all points.

    iel = 1

    DO ipoint = 1,UBOUND(Dpoints,2)
    
      ! Get the element number that contains the point. 
      IF (PRESENT(Ielements)) THEN
        ! Either we have that...
        iel = Ielements (ipoint)
      ELSE
        ! Or we have to find the element. Probably, the caller gave us a
        ! hint where we can start the search...
        IF (PRESENT(IelementsHint)) THEN
          iel = IelementsHint (ipoint)
        END IF
        
        ! Otherwise, we use iel from the previous iteration as starting point
        ! to search for the element. Use raytracing search to find the element
        ! containing the point.
        CALL tsrch_getElem_raytrace2D (&
          Dpoints(:,ipoint),rvectorScalar%p_rspatialDiscretisation%p_rtriangulation,iel)

        ! Ok, not found... Brute force search
        IF (iel .EQ. 0) THEN
          CALL tsrch_getElem_BruteForce (Dpoints(:,ipoint), &
            rvectorScalar%p_rspatialDiscretisation%p_rtriangulation,iel)
        END IF
        
        IF (iel .EQ. 0) THEN
        
          ! We really have not luck here... Are nonmesh-points allowed?
          IF (cnonmesh .EQ. FEVL_NONMESHPTS_NEARBY) THEN
          
            ! Yes. Find the closest element!
            CALL tsrch_getNearestElem_BruteForce (Dpoints(:,ipoint), &
              rvectorScalar%p_rspatialDiscretisation%p_rtriangulation,iel)
              
            ! The evaluation routine then evaluates the FE function
            ! outside of the element...
            
          ELSE IF (cnonmesh .EQ. FEVL_NONMESHPTS_ZERO) THEN
            
            ! Assume zero here.
            Dvalues(ipoint) = dval
            
            CYCLE
            
          END IF
          
        END IF

      END IF

      IF (iel .EQ. 0) THEN
        CALL output_line ('Point '//TRIM(sys_siL(ipoint,10))//' not found!', &
                          OU_CLASS_ERROR,OU_MODE_STD,'fevl_evaluate')
        CYCLE
      END IF
    
      ! Get the type of the element iel
      IF (ASSOCIATED(p_IelementDistr)) THEN
        ieltype = p_RelementDistribution(p_IelementDistr(iel))%itrialElement

        ! Get the number of local DOF's for trial and test functions
        indof = elem_igetNDofLoc(ieltype)
        
        ! Number of vertices on the element
        nve = elem_igetNVE(ieltype)
        
        ! Type of transformation from/to the reference element
        ctrafoType = elem_igetTrafoType(ieltype)
        
        ! Get the element evaluation tag; necessary for the preparation of the element
        cevaluationTag = elem_getEvaluationTag(ieltype)
      END IF
        
      ! Calculate the global DOF's on that element into IdofsTest.
      CALL dof_locGlobMapping (rvectorScalar%p_rspatialDiscretisation, &
          iel,.FALSE.,Idofs)
     
      ! Get the element shape information
      CALL elprep_prepareForEvaluation (revalElement, EL_EVLTAG_COORDS, &
          rvectorScalar%p_rspatialDiscretisation%p_rtriangulation, iel, ctrafoType)
          
      ! Calculate the transformation of the point to the reference element
      CALL trafo_calcRefCoords (ctrafoType,revalElement%Dcoords,&
          Dpoints(:,ipoint),DparPoint)

      ! Now calculate everything else what is necessary for the element    
      CALL elprep_prepareForEvaluation (revalElement, &
          IAND(cevaluationTag,NOT(EL_EVLTAG_COORDS)), &
          rvectorScalar%p_rspatialDiscretisation%p_rtriangulation, iel, &
          ctrafoType, DparPoint, Dpoints(:,ipoint))

      ! Call the element to calculate the values of the basis functions
      ! in the point.
      CALL elem_generic2 (ieltype, revalElement, Bder, Dbas)
      
      ! Combine the basis functions to get the function value.
      dval = 0.0_DP
      IF (rvectorScalar%cdataType .EQ. ST_DOUBLE) THEN
      
        ! Now that we have the basis functions, we want to have the function values.
        ! We get them by multiplying the FE-coefficients with the values of the
        ! basis functions and summing up.
        !          
        ! Calculate the value in the point
        DO ibas = 1,indof
          dval = dval + p_Ddata(Idofs(ibas)) * Dbas(ibas,iderType)
        END DO
      
      ELSE IF (rvectorScalar%cdataType .EQ. ST_SINGLE) THEN
      
        ! Now that we have the basis functions, we want to have the function values.
        ! We get them by multiplying the FE-coefficients with the values of the
        ! basis functions and summing up.
        !
        ! Calculate the value in the point
        DO ibas = 1,indof
          dval = dval + p_Fdata(Idofs(ibas)) * Dbas(ibas,iderType)
        END DO
        
      END IF

      ! Save the value in the point
      Dvalues(ipoint) = dval
      
    END DO ! ipoint

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE fevl_evaluate_mult1 (iderType, Dvalues, rvectorScalar, ielement, &
      DpointsRef, Dpoints)
                                      
!<description>
  ! This is a rather general finite element evaluation
  ! routine. It allows to evaluate a general (scalar) FE function specified
  ! by rvectorScalar in a set of points Dpoints on one element ielement. 
  ! The values of the FE function are written into Dvalues.
!</description>

!<input>
  ! Type of function value to evaluate. One of the DER_xxxx constants,
  ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
  INTEGER, INTENT(IN)                            :: iderType

  ! The scalar solution vector that is to be evaluated.
  TYPE(t_vectorScalar), INTENT(IN)              :: rvectorScalar
  
  ! The element number containing the points in Dpoints.
  INTEGER(PREC_ELEMENTIDX), INTENT(IN) :: ielement
  
  ! OPTIONAL: Coordinates of the points on the reference element.
  ! If not specified, the coordinates are automatically calculated.
  ! Either Dpoints or DpointsRef must be specified!
  ! DIMENSION(1..ndim,1..npoints)
  REAL(DP), DIMENSION(:,:), INTENT(IN), OPTIONAL :: DpointsRef
  
  ! OPTIONAL: A list of points where to evaluate. All points must be inside
  ! of element ielement. 
  ! If not specified, the coordinates are automatically calculated.
  ! Either Dpoints or DpointsRef must be specified!
  ! DIMENSION(1..ndim,1..npoints)
  REAL(DP), DIMENSION(:,:), INTENT(IN), OPTIONAL :: Dpoints
!</input>

!<output>
  ! Values of the FE function at the points specified by Dpoints.
  REAL(DP), DIMENSION(:), INTENT(OUT) :: Dvalues
!</output>

!</subroutine>

    ! local variables
    INTEGER :: ipoint,ieltype,indof,nve,ibas,npoints
    INTEGER(I32), DIMENSION(:), POINTER :: p_IelementDistr
    LOGICAL, DIMENSION(EL_MAXNDER) :: Bder
    REAL(DP) :: dval
    
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    REAL(SP), DIMENSION(:), POINTER :: p_Fdata
    
    ! Transformation
    INTEGER(I32) :: ctrafoType
    REAL(DP), DIMENSION(TRAFO_MAXDIMREFCOORD) :: DparPoint
    
    ! Values of basis functions and DOF's
    REAL(DP), DIMENSION(EL_MAXNBAS,EL_MAXNDER) :: Dbas
    INTEGER(PREC_DOFIDX), DIMENSION(EL_MAXNBAS) :: Idofs
    
    ! List of element distributions in the discretisation structure
    TYPE(t_elementDistribution), DIMENSION(:), POINTER :: p_RelementDistribution

    ! Evaluation structure and tag
    TYPE(t_evalElement) :: revalElement
    INTEGER(I32) :: cevaluationTag

    ! Ok, slow but general.
    
    ! Are points given?
    IF ((.NOT. PRESENT(Dpoints)) .AND. (.NOT. PRESENT(DpointsRef))) THEN
      CALL output_line ('Evaluation points not specified!', &
                        OU_CLASS_ERROR,OU_MODE_STD,'fevl_evaluate_mult')  
    END IF
    
    IF (PRESENT(DpointsRef)) THEN
      npoints = UBOUND(DpointsRef,2)
    ELSE
      npoints = UBOUND(Dpoints,2)
    END IF
    
    p_RelementDistribution => rvectorScalar%p_rspatialDiscretisation%RelementDistribution
    
    ! For uniform discretisations, we get the element type in advance...
    IF (rvectorScalar%p_rspatialDiscretisation%ccomplexity .EQ. SPDISC_UNIFORM) THEN
      
      ! Element type
      ieltype = p_RelementDistribution(1)%itrialElement

    ELSE
      CALL storage_getbase_int (rvectorScalar%p_rspatialDiscretisation%h_IelementDistr,&
          p_IelementDistr)

      ieltype = p_RelementDistribution(p_IelementDistr(ielement))%itrialElement
    END IF

    ! Get the number of local DOF's for trial and test functions
    indof = elem_igetNDofLoc(ieltype)
    
    ! Number of vertices on the element
    nve = elem_igetNVE(ieltype)
    
    ! Type of transformation from/to the reference element
    ctrafoType = elem_igetTrafoType(ieltype)
    
    ! Get the element evaluation tag; necessary for the preparation of the element
    cevaluationTag = elem_getEvaluationTag(ieltype)

    ! Get the data vector
    SELECT CASE (rvectorScalar%cdataType)
    CASE (ST_DOUBLE) 
      CALL lsyssc_getbase_double(rvectorScalar,p_Ddata)
    CASE (ST_SINGLE)
      CALL lsyssc_getbase_single(rvectorScalar,p_Fdata)
    CASE DEFAULT
      CALL output_line ('Unsupported vector precision!',&
          OU_CLASS_ERROR,OU_MODE_STD,'fevl_evaluate')
      CALL sys_halt()
    END SELECT
    
    ! What to evaluate?
    Bder = .FALSE.
    Bder(iderType) = .TRUE.
    
    ! Calculate the global DOF's on that element into IdofsTest.
    CALL dof_locGlobMapping (rvectorScalar%p_rspatialDiscretisation, &
        ielement,.FALSE.,Idofs)
        
    ! Get the element shape information
    CALL elprep_prepareForEvaluation (revalElement, EL_EVLTAG_COORDS, &
        rvectorScalar%p_rspatialDiscretisation%p_rtriangulation, ielement, ctrafoType)
  
    ! We loop over all points
    DO ipoint = 1,npoints
    
      IF (PRESENT(DpointsRef)) THEN
        DparPoint(1:UBOUND(DpointsRef,1)) = DpointsRef(:,ipoint)
      ELSE
        ! Calculate the transformation of the point to the reference element
        CALL trafo_calcRefCoords (ctrafoType,revalElement%Dcoords,Dpoints(:,ipoint),DparPoint)
      END IF

      ! Now calculate everything else what is necessary for the element.
      ! Don't calculate the shape of the cell again since we did this in advance above.
      IF (PRESENT(Dpoints)) THEN
        CALL elprep_prepareForEvaluation (revalElement, &
            IAND(cevaluationTag,NOT(EL_EVLTAG_COORDS)), &
            rvectorScalar%p_rspatialDiscretisation%p_rtriangulation, ielement, &
            ctrafoType, DparPoint, Dpoints(:,ipoint))
      ELSE
        CALL elprep_prepareForEvaluation (revalElement, &
            IAND(cevaluationTag,NOT(EL_EVLTAG_COORDS)), &
            rvectorScalar%p_rspatialDiscretisation%p_rtriangulation, ielement, &
            ctrafoType, DparPoint)
      END IF

      ! Call the element to calculate the values of the basis functions
      ! in the point.
      CALL elem_generic2 (ieltype, revalElement, Bder, Dbas)
      
      dval = 0.0_DP
      IF (rvectorScalar%cdataType .EQ. ST_DOUBLE) THEN
      
        ! Now that we have the basis functions, we want to have the function values.
        ! We get them by multiplying the FE-coefficients with the values of the
        ! basis functions and summing up.
        !          
        ! Calculate the value in the point
        DO ibas = 1,indof
          dval = dval + p_Ddata(Idofs(ibas)) * Dbas(ibas,iderType)
        END DO
      
      ELSE IF (rvectorScalar%cdataType .EQ. ST_SINGLE) THEN
      
        ! Now that we have the basis functions, we want to have the function values.
        ! We get them by multiplying the FE-coefficients with the values of the
        ! basis functions and summing up.
        !
        ! Calculate the value in the point
        DO ibas = 1,indof
          dval = dval + p_Fdata(Idofs(ibas)) * Dbas(ibas,iderType)
        END DO
        
      END IF

      ! Save the value in the point
      Dvalues(ipoint) = dval
      
    END DO ! ipoint

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE fevl_evaluate_mult2 (rvectorScalar, Dcoords, Djac, Ddetj, &
                  ieltyp, IdofsTrial, npoints,  Dpoints, iderType,&
                  Dvalues, ItwistIndex)
                                      
!<description>
  ! DEPRECATED!
  ! This routine allows to evaluate a finite element solution vector
  ! rvectorScalar simultaneously in multiple points on one elements in a 
  ! discretisation.
  ! The caller must provide all necessary information for the evaluation in the
  ! parameters.
!</description>

!<input>
  ! The scalar solution vector that is to be evaluated.
  TYPE(t_vectorScalar), INTENT(IN)               :: rvectorScalar
  
  ! The FE function must be discretised with the same trial functions on all
  ! elements where it should be evaluated here. ieltyp defines the type
  ! of FE trial function that was used for the discretisation on those 
  ! elements that we are concerning here.
  INTEGER(I32), INTENT(IN)                       :: ieltyp

  ! A list of the corner vertices of the element.
  ! array [1..NDIM2D,1..TRIA_MAXNVE2D] of double
  REAL(DP), DIMENSION(:,:), INTENT(IN)           :: Dcoords
  
  ! The Jacobian matrix of the mapping between the reference and the
  ! real element, for all points on the element.
  ! array [1..TRAFO_NJACENTRIES,1..npointsPerElement]
  REAL(DP), DIMENSION(:,:),INTENT(IN)            :: Djac
  
  ! The Jacobian determinant of the mapping of each point from the
  ! reference element to the real element.
  ! array [1..npointsPerElement]
  REAL(DP), DIMENSION(:), INTENT(IN)             :: Ddetj
  
  ! An array accepting the DOF's on the element in the trial space
  ! of the FE function.
  ! DIMENSION(\#local DOF's in trial space)
  INTEGER(PREC_DOFIDX), DIMENSION(:), INTENT(IN) :: IdofsTrial
  
  ! Number of points on the element where to evalate the function
  INTEGER, INTENT(IN) :: npoints
  
  ! Array with coordinates of the points where to evaluate.
  ! DIMENSION(NDIM2D,npoints).
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
  REAL(DP), DIMENSION(:,:), INTENT(IN)           :: Dpoints

  ! Type of function value to evaluate. One of the DER_xxxx constants,
  ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
  INTEGER, INTENT(IN)                            :: iderType

  ! OPTIONAL: Twist index bitfield of the element. Defines for every edge its
  ! orientation.
  ! Can be omitted if the element does not need this information.
  INTEGER(I32), INTENT(IN), OPTIONAL :: itwistIndex

!</input>

!<output>
  ! Values of the FE function at the points specified by Dpoints.
  ! DIMENSION(npoints).
  REAL(DP), DIMENSION(:), INTENT(OUT)            :: Dvalues
!</output>

!</subroutine>

  ! local variables
  LOGICAL, DIMENSION(EL_MAXNDER) :: Bder
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: DbasTrial
  INTEGER :: indofTrial
  REAL(DP) :: dval
  INTEGER :: ipoint,ibas
  REAL(DP), DIMENSION(:), POINTER :: p_Ddata
  REAL(SP), DIMENSION(:), POINTER :: p_Fdata
  
  ! What to evaluate?
  Bder = .FALSE.
  Bder(iderType) = .TRUE.
  
  ! Allocate memory for the basis function values
  indofTrial = elem_igetNDofLoc(ieltyp)
  ALLOCATE(DbasTrial(indofTrial,elem_getMaxDerivative(ieltyp),npoints))
  
  ! Evaluate the basis functions
  CALL elem_generic_mult (ieltyp, Dcoords, Djac, Ddetj, &
                         Bder, DbasTrial, npoints, Dpoints, itwistIndex)  
  
  IF (rvectorScalar%cdataType .EQ. ST_DOUBLE) THEN
  
    ! Get the data array from the vector
    CALL lsyssc_getbase_double(rvectorScalar,p_Ddata)
    
    ! Now that we have the basis functions, we want to have the function values.
    ! We get them by multiplying the FE-coefficients with the values of the
    ! basis functions and summing up.
    DO ipoint = 1,npoints
      ! Calculate the value in the point
      dval = 0.0_DP
      DO ibas = 1,indofTrial
        dval = dval + &
            p_Ddata(IdofsTrial(ibas)) * DbasTrial(ibas,iderType,ipoint)
      END DO
      ! Save the value in the point
      Dvalues(ipoint) = dval
    END DO
    
  ELSE IF (rvectorScalar%cdataType .EQ. ST_SINGLE) THEN
  
    ! Get the data array from the vector
    CALL lsyssc_getbase_single(rvectorScalar,p_Fdata)
    
    ! Now that we have the basis functions, we want to have the function values.
    ! We get them by multiplying the FE-coefficients with the values of the
    ! basis functions and summing up.
    DO ipoint = 1,npoints
      ! Calculate the value in the point
      dval = 0.0_DP
      DO ibas = 1,indofTrial
        dval = dval + &
            p_Fdata(IdofsTrial(ibas)) * DbasTrial(ibas,iderType,ipoint)
      END DO
      ! Save the value in the point
      Dvalues(ipoint) = dval
    END DO
    
  ELSE
    PRINT *,'fevl_evaluate_sim: Unsupported vector precision'
    CALL sys_halt()
  END IF
  
  ! Release memory, finish
  DEALLOCATE(DbasTrial)
  
  END SUBROUTINE

!  ! ***************************************************************************
!
!!<subroutine>
!
!  SUBROUTINE fevl_evaluate_sim (iderType, Dvalues, rvectorScalar, Dpoints, &
!      Ielements, DpointsRef)
!                                      
!!<description>
!  ! This is a rather general finite element evaluation
!  ! routine. It allows to evaluate a general (scalar) FE function specified
!  ! by rvectorScalar in a set of points Dpoints on a set of elements Ielements. 
!  ! The values of the FE function are written into Dvalues.
!!</description>
!
!!<input>
!  ! Type of function value to evaluate. One of the DER_xxxx constants,
!  ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
!  INTEGER, INTENT(IN)                            :: iderType
!
!  ! The scalar solution vector that is to be evaluated.
!  TYPE(t_vectorScalar), INTENT(IN)              :: rvectorScalar
!  
!  ! A list of points where to evaluate. All points must be inside
!  ! of element ielement.
!  ! DIMENSION(1..ndim,1..npoints,1..nelements)
!  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: Dpoints
!  
!  ! A list of elements containing the points in Dpoints.
!  ! All elements in this list must be of the same type!!!
!  INTEGER(PREC_ELEMENTIDX), DIMENSION(:), INTENT(IN) :: Ielements
!  
!  ! OPTIONAL: Coordinates of the points on the reference element.
!  ! If not specified, the coordinates are automatically calculated.
!  ! DIMENSION(1..ndim,1..npoints,1..nelements)
!  REAL(DP), DIMENSION(:,:,:), INTENT(IN), TARGET, OPTIONAL :: DpointsRef
!!</input>
!
!!<output>
!  ! Values of the FE function at the points specified by Dpoints.
!  ! DIMENSION(1..npoints,1..nelements)
!  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dvalues
!!</output>
!
!!</subroutine>
!
!    ! local variables
!    LOGICAL :: bnonpar
!    INTEGER :: ipoint,ieltype,indof,nve,ibas,iel,ndim,ntwistsize
!    INTEGER(I32), DIMENSION(:), POINTER :: p_IelementDistr
!    LOGICAL, DIMENSION(EL_MAXNDER) :: Bder
!    REAL(DP) :: dval
!    REAL(DP), DIMENSION(:,:,:), POINTER :: p_DpointsRef
!    
!    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
!    REAL(SP), DIMENSION(:), POINTER :: p_Fdata
!    
!    ! Triangulation information
!    REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoords
!    INTEGER(PREC_POINTIDX), DIMENSION(:,:), POINTER :: p_IverticesAtElement
!    
!    ! Transformation
!    REAL(DP), DIMENSION(:,:,:),ALLOCATABLE :: Djac
!    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Ddetj
!    INTEGER(I32) :: ctrafoType
!    
!    ! Values of basis functions and DOF's
!    REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: Dbas
!    INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE :: Idofs
!    
!    ! Coordinates of the corners of one element
!    REAL(DP), DIMENSION(:,:,:),ALLOCATABLE :: Dcoord
!
!    ! List of element distributions in the discretisation structure
!    TYPE(t_elementDistribution), DIMENSION(:), POINTER :: p_RelementDistribution
!
!    ! Twist index array to define the orientation of faces/edges
!    INTEGER(I32), DIMENSION(:,:), ALLOCATABLE :: ItwistIndex
!
!    ! Ok, slow but general.
!    
!    ! Get triangulation information
!    CALL storage_getbase_double2d (&
!        rvectorScalar%p_rspatialDiscretisation%p_rtriangulation%h_DvertexCoords,&
!        p_DvertexCoords)
!    CALL storage_getbase_int2d (&
!        rvectorScalar%p_rspatialDiscretisation%p_rtriangulation%h_IverticesAtElement,&
!        p_IverticesAtElement)
!    
!    p_RelementDistribution => rvectorScalar%p_rspatialDiscretisation%RelementDistribution
!
!    ! For uniform discretisations, we get the element type in advance...
!    IF (rvectorScalar%p_rspatialDiscretisation%ccomplexity .EQ. SPDISC_UNIFORM) THEN
!      
!      ! Element type
!      ieltype = p_RelementDistribution(1)%itrialElement
!
!      ! Get the number of local DOF's for trial and test functions
!      indof = elem_igetNDofLoc(ieltype)
!      
!      ! Number of vertices on the element
!      nve = elem_igetNVE(ieltype)
!      
!      ! Type of transformation from/to the reference element
!      ctrafoType = elem_igetTrafoType(ieltype)
!      
!      ! Element nonparametric?
!      bnonpar = elem_isNonparametric(ieltype)
!      
!      NULLIFY(p_IelementDistr)
!    ELSE
!      CALL storage_getbase_int (rvectorScalar%p_rspatialDiscretisation%h_IelementDistr,&
!          p_IelementDistr)
!    END IF
!    
!    ! Get the data vector
!    SELECT CASE (rvectorScalar%cdataType)
!    CASE (ST_DOUBLE) 
!      CALL lsyssc_getbase_double(rvectorScalar,p_Ddata)
!    CASE (ST_SINGLE)
!      CALL lsyssc_getbase_single(rvectorScalar,p_Fdata)
!    CASE DEFAULT
!      CALL output_line ('Unsupported vector precision!',&
!          OU_CLASS_ERROR,OU_MODE_STD,'fevl_evaluate')
!      CALL sys_halt()
!    END SELECT
!    
!    ! What to evaluate?
!    Bder = .FALSE.
!    Bder(iderType) = .TRUE.
!    
!    ! Get the type of the element ielement
!    IF (ASSOCIATED(p_IelementDistr)) THEN
!      ! As all elements have the same type, we get the element
!      ! characteristics by checking the first element.
!      ieltype = p_RelementDistribution(p_IelementDistr(Ielements(1)))%itrialElement
!
!      ! Get the number of local DOF's for trial and test functions
!      indof = elem_igetNDofLoc(ieltype)
!      
!      ! Number of vertices on the element
!      nve = elem_igetNVE(ieltype)
!      
!      ! Type of transformation from/to the reference element
!      ctrafoType = elem_igetTrafoType(ieltype)
!      
!      ! Element nonparametric?
!      bnonpar = elem_isNonparametric(ieltype)
!      
!    END IF
!      
!    ! Calculate the global DOF's on that element into IdofsTest.
!    ALLOCATE(Idofs(indof,SIZE(Ielements)))
!    CALL dof_locGlobMapping_mult(rvectorScalar%p_rspatialDiscretisation, &
!        Ielements, .FALSE.,Idofs)
!        
!    ! Get the coordinates forming the elements
!    ALLOCATE(Dcoord(UBOUND(p_DvertexCoords,1),NVE,SIZE(Ielements)))
!    DO iel = 1,SIZE(Ielements)
!      Dcoord(:,1:nve,iel) = &
!        p_DvertexCoords(:,p_IverticesAtElement(:,Ielements(iel)))
!    END DO
!    
!    ! Get the coordinates of all points on the reference element
!    IF (PRESENT(DpointsRef)) THEN
!      p_DpointsRef => DpointsRef
!    ELSE
!      ALLOCATE(p_DpointsRef(UBOUND(Dpoints,1),UBOUND(Dpoints,2),&
!                            UBOUND(Dpoints,3)))
!      ! Calculate the transformation of the point to the reference element
!      DO iel = 1,SIZE(Ielements)
!        DO ipoint = 1,UBOUND(Dpoints,2)
!          CALL trafo_calcRefCoords (ctrafoType,Dcoord(:,:,iel),&
!              Dpoints(:,ipoint,iel),p_DpointsRef(:,ipoint,iel))
!        END DO
!      END DO
!    END IF
!    
!    ! Calculate the transformation of all the points from the reference
!    ! to the real element(s).
!    ndim = UBOUND(Dcoord,1)
!    ALLOCATE(Djac(ndim*ndim,UBOUND(Dpoints,2),UBOUND(Dpoints,3)))
!    ALLOCATE(Ddetj(UBOUND(Dpoints,2),UBOUND(Dpoints,3)))
!    CALL trafo_calctrafo_sim (elem_igetTrafoType(ieltype),SIZE(Ielements),&
!        UBOUND(Dpoints,2),Dcoord,&
!        p_DpointsRef,Djac,Ddetj)    
!  
!    ! Does the element need twist indices?
!    ntwistsize = elem_getTwistIndexSize(ieltype)
!
!    ! If necessary, calculate the twist index array. The element may need it.
!    ! We always allocate so that we have something we can pass to the element...
!    ALLOCATE(ItwistIndex(MAX(1,ntwistSize),SIZE(Ielements)))
!    IF (ntwistSize .NE. 0) THEN
!      CALL trafo_calcTwistIndices(&
!        rvectorScalar%p_rspatialDiscretisation%p_rtriangulation,Ielements,ItwistIndex)
!    END IF
!
!    ! Calculate the values of the basis functions in the given points.
!    ALLOCATE(Dbas(indof,&
!             elem_getMaxDerivative(ieltype),&
!             UBOUND(Dpoints,2), UBOUND(Dpoints,3)))
!    IF (bnonpar) THEN
!      CALL elem_generic_sim (ieltype, Dcoord, Djac, Ddetj, &
!                           Bder, Dbas, UBOUND(Dpoints,2), UBOUND(Dpoints,3), &
!                           Dpoints,ItwistIndex)    
!    ELSE
!      CALL elem_generic_sim (ieltype, Dcoord, Djac, Ddetj, &
!                           Bder, Dbas, UBOUND(Dpoints,2), UBOUND(Dpoints,3), &
!                           p_DpointsRef,ItwistIndex)    
!    END IF
!  
!    ! Calculate the desired values. We loop over all points and all elements
!    IF (rvectorScalar%cdataType .EQ. ST_DOUBLE) THEN
!      DO iel = 1, UBOUND(Dpoints,3)
!        DO ipoint = 1,UBOUND(Dpoints,2)
!      
!          dval = 0.0_DP
!        
!          ! Now that we have the basis functions, we want to have the function values.
!          ! We get them by multiplying the FE-coefficients with the values of the
!          ! basis functions and summing up.
!          !          
!          ! Calculate the value in the point
!          DO ibas = 1,indof
!            dval = dval + p_Ddata(Idofs(ibas,iel)) * Dbas(ibas,iderType,ipoint,iel)
!          END DO
!        
!          ! Save the value in the point
!          Dvalues(ipoint,iel) = dval
!        
!        END DO ! ipoint
!      END DO ! iel
!
!    ELSE IF (rvectorScalar%cdataType .EQ. ST_SINGLE) THEN
!    
!      DO iel = 1, UBOUND(Dpoints,3)
!        DO ipoint = 1,UBOUND(Dpoints,2)
!          
!          ! Now that we have the basis functions, we want to have the function values.
!          ! We get them by multiplying the FE-coefficients with the values of the
!          ! basis functions and summing up.
!          !
!          ! Calculate the value in the point
!          DO ibas = 1,indof
!            dval = dval + p_Fdata(Idofs(ibas,iel)) * Dbas(ibas,iderType,ipoint,iel)
!          END DO
!        
!          ! Save the value in the point
!          Dvalues(ipoint,iel) = dval
!        
!        END DO ! ipoint
!      END DO ! iel
!      
!    END IF
!    
!    ! Release allocated memory
!    DEALLOCATE(ItwistIndex)
!    DEALLOCATE(Dbas)
!    DEALLOCATE(Ddetj)
!    DEALLOCATE(Djac)
!    IF (.NOT. PRESENT(DpointsRef)) THEN
!      DEALLOCATE(p_DpointsRef)
!    END IF
!    DEALLOCATE(Dcoord)
!    DEALLOCATE(Idofs)
!
!  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE fevl_evaluate_sim1 (iderType, Dvalues, rvectorScalar, Dpoints, &
      Ielements, DpointsRef)
                                      
!<description>
  ! This is a rather general finite element evaluation
  ! routine. It allows to evaluate a general (scalar) FE function specified
  ! by rvectorScalar in a set of points Dpoints on a set of elements Ielements. 
  ! The values of the FE function are written into Dvalues.
!</description>

!<input>
  ! Type of function value to evaluate. One of the DER_xxxx constants,
  ! e.g. DER_FUNC for function values, DER_DERIV_X for x-derivatives etc.
  INTEGER, INTENT(IN)                            :: iderType

  ! The scalar solution vector that is to be evaluated.
  TYPE(t_vectorScalar), INTENT(IN)              :: rvectorScalar
  
  ! A list of points where to evaluate. All points must be inside
  ! of element ielement.
  ! DIMENSION(1..ndim,1..npoints,1..nelements)
  REAL(DP), DIMENSION(:,:,:), INTENT(IN) :: Dpoints
  
  ! A list of elements containing the points in Dpoints.
  ! All elements in this list must be of the same type!!!
  INTEGER(PREC_ELEMENTIDX), DIMENSION(:), INTENT(IN) :: Ielements
  
  ! OPTIONAL: Coordinates of the points on the reference element.
  ! If not specified, the coordinates are automatically calculated.
  ! DIMENSION(1..ndim,1..npoints,1..nelements)
  REAL(DP), DIMENSION(:,:,:), INTENT(IN), TARGET, OPTIONAL :: DpointsRef
!</input>

!<output>
  ! Values of the FE function at the points specified by Dpoints.
  ! DIMENSION(1..npoints,1..nelements)
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: Dvalues
!</output>

!</subroutine>

    ! local variables
    LOGICAL :: bnonpar
    INTEGER :: ipoint,ieltype,indof,nve,ibas,iel
    INTEGER(I32), DIMENSION(:), POINTER :: p_IelementDistr
    LOGICAL, DIMENSION(EL_MAXNDER) :: Bder
    REAL(DP) :: dval
    
    REAL(DP), DIMENSION(:,:,:), POINTER :: p_DpointsRef
    
    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
    REAL(SP), DIMENSION(:), POINTER :: p_Fdata
    
    ! The triangulation structure - to shorten some things...
    TYPE(t_triangulation), POINTER :: p_rtriangulation
    
    ! Transformation
    INTEGER(I32) :: ctrafoType
    
    ! Values of basis functions and DOF's
    REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: Dbas
    INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE :: Idofs

    ! Element evaluation set that collects element specific information
    ! during the evaluation
    TYPE(t_evalElementSet)  :: revalElementSet
    INTEGER(I32) :: cevaluationTag
    
    ! List of element distributions in the discretisation structure
    TYPE(t_elementDistribution), DIMENSION(:), POINTER :: p_RelementDistribution

    ! Ok, slow but general.
    
    ! Get triangulation information
    p_rtriangulation => rvectorScalar%p_rspatialDiscretisation%p_rtriangulation
    
    p_RelementDistribution => rvectorScalar%p_rspatialDiscretisation%RelementDistribution

    ! For uniform discretisations, we get the element type in advance...
    IF (rvectorScalar%p_rspatialDiscretisation%ccomplexity .EQ. SPDISC_UNIFORM) THEN
      
      ! Element type
      ieltype = p_RelementDistribution(1)%itrialElement

      ! Get the number of local DOF's for trial and test functions
      indof = elem_igetNDofLoc(ieltype)
      
      ! Number of vertices on the element
      nve = elem_igetNVE(ieltype)
      
      ! Type of transformation from/to the reference element
      ctrafoType = elem_igetTrafoType(ieltype)
      
      ! Element nonparametric?
      bnonpar = elem_isNonparametric(ieltype)
      
      NULLIFY(p_IelementDistr)
    ELSE
      CALL storage_getbase_int (rvectorScalar%p_rspatialDiscretisation%h_IelementDistr,&
          p_IelementDistr)
    END IF
    
    ! Get the data vector
    SELECT CASE (rvectorScalar%cdataType)
    CASE (ST_DOUBLE) 
      CALL lsyssc_getbase_double(rvectorScalar,p_Ddata)
    CASE (ST_SINGLE)
      CALL lsyssc_getbase_single(rvectorScalar,p_Fdata)
    CASE DEFAULT
      CALL output_line ('Unsupported vector precision!',&
          OU_CLASS_ERROR,OU_MODE_STD,'fevl_evaluate')
      CALL sys_halt()
    END SELECT
    
    ! What to evaluate?
    Bder = .FALSE.
    Bder(iderType) = .TRUE.
    
    ! Get the type of the element ielement
    IF (ASSOCIATED(p_IelementDistr)) THEN
      ! As all elements have the same type, we get the element
      ! characteristics by checking the first element.
      ieltype = p_RelementDistribution(p_IelementDistr(Ielements(1)))%itrialElement

      ! Get the number of local DOF's for trial and test functions
      indof = elem_igetNDofLoc(ieltype)
      
      ! Number of vertices on the element
      nve = elem_igetNVE(ieltype)
      
      ! Type of transformation from/to the reference element
      ctrafoType = elem_igetTrafoType(ieltype)
      
      ! Element nonparametric?
      bnonpar = elem_isNonparametric(ieltype)
      
    END IF
      
    ! Calculate the global DOF's on that element into IdofsTest.
    ALLOCATE(Idofs(indof,SIZE(Ielements)))
    CALL dof_locGlobMapping_mult(rvectorScalar%p_rspatialDiscretisation, &
        Ielements, .FALSE.,Idofs)
        
    ! Get the element evaluation tag of all FE spaces. We need it to evaluate
    ! the elements later. All of them can be combined with OR, what will give
    ! a combined evaluation tag. 
    cevaluationTag = elem_getEvaluationTag(ieltype)
    
    ! Don't create coordinates on the reference element; we do this manually!
    cevaluationTag = IAND(cevaluationTag,NOT(EL_EVLTAG_REFPOINTS))
                    
    ! Get the coordinates of all points on the reference element
    IF (.NOT. PRESENT(DpointsRef)) THEN
      p_DpointsRef => DpointsRef
    ELSE
      ALLOCATE(p_DpointsRef(UBOUND(Dpoints,1),UBOUND(Dpoints,2),&
               UBOUND(Dpoints,3)))
      ! Calculate the transformation of the point to the reference element
      DO iel = 1,SIZE(Ielements)
        DO ipoint = 1,UBOUND(Dpoints,2)
          CALL trafo_calcRefCoords (ctrafoType,revalElementSet%p_Dcoords(:,:,iel),&
              Dpoints(:,ipoint,iel),p_DpointsRef(:,ipoint,iel))
        END DO
      END DO
    END IF
    
    ! Calculate all information that is necessary to evaluate the finite element
    ! on all cells of our subset. This includes the coordinates of the points
    ! on the cells.
    CALL elprep_prepareSetForEvaluation (revalElementSet,&
        cevaluationTag, p_rtriangulation, Ielements, ctrafoType,DpointsRef=p_DpointsRef)

    ! Calculate the values of the basis functions in the given points.
    ALLOCATE(Dbas(indof,&
             elem_getMaxDerivative(ieltype),&
             UBOUND(Dpoints,2), UBOUND(Dpoints,3)))
    CALL elem_generic_sim2 (ieltype, revalElementSet, Bder, Dbas)
             
    ! Calculate the desired values. We loop over all points and all elements
    IF (rvectorScalar%cdataType .EQ. ST_DOUBLE) THEN
      DO iel = 1, UBOUND(Dpoints,3)
        DO ipoint = 1,UBOUND(Dpoints,2)
      
          dval = 0.0_DP
        
          ! Now that we have the basis functions, we want to have the function values.
          ! We get them by multiplying the FE-coefficients with the values of the
          ! basis functions and summing up.
          !          
          ! Calculate the value in the point
          DO ibas = 1,indof
            dval = dval + p_Ddata(Idofs(ibas,iel)) * Dbas(ibas,iderType,ipoint,iel)
          END DO
        
          ! Save the value in the point
          Dvalues(ipoint,iel) = dval
        
        END DO ! ipoint
      END DO ! iel

    ELSE IF (rvectorScalar%cdataType .EQ. ST_SINGLE) THEN
    
      DO iel = 1, UBOUND(Dpoints,3)
        DO ipoint = 1,UBOUND(Dpoints,2)
          
          ! Now that we have the basis functions, we want to have the function values.
          ! We get them by multiplying the FE-coefficients with the values of the
          ! basis functions and summing up.
          !
          ! Calculate the value in the point
          DO ibas = 1,indof
            dval = dval + p_Fdata(Idofs(ibas,iel)) * Dbas(ibas,iderType,ipoint,iel)
          END DO
        
          ! Save the value in the point
          Dvalues(ipoint,iel) = dval
        
        END DO ! ipoint
      END DO ! iel
      
    END IF
    
    ! Release allocated memory
    ! Remove the reference to DpointsRef again
    DEALLOCATE(Dbas)

    IF (.NOT. PRESENT(DpointsRef)) THEN
      DEALLOCATE(revalElementSet%p_DpointsRef)
    END IF
  
    CALL elprep_releaseElementSet(revalElementSet)
    
    IF (.NOT. PRESENT(DpointsRef)) THEN
      DEALLOCATE(p_DpointsRef)
    END IF
    DEALLOCATE(Idofs)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE fevl_evaluate_sim2 (rvectorScalar, Dcoords, Djac, Ddetj, &
                  ieltyp, IdofsTrial, npoints,  nelements, Dpoints, iderType,&
                  Dvalues,ItwistIndexEdges)
                                      
!<description>
  ! DEPRECATED!
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
  INTEGER(I32), INTENT(IN)                      :: ieltyp

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

  ! OPTIONAL: List of twist indices. Defines for every element the
  ! orgientation of the edges.
  ! Can be omitted if the element does not need it.
  ! Array with DIMENSION(1:NVE/NVA,nelements)
  INTEGER(I32), DIMENSION(:), INTENT(IN), OPTIONAL :: ItwistIndexEdges
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
                         Bder, DbasTrial, npoints, nelements, Dpoints,ItwistIndexEdges)  
  
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
    CALL sys_halt()
  END IF
  
  ! Release memory, finish
  DEALLOCATE(DbasTrial)

  END SUBROUTINE

  ! ***************************************************************************

!<subroutine>

  SUBROUTINE fevl_evaluate_sim3 (rvectorScalar, revalElementSet,&
                  ieltyp, IdofsTrial, iderType, Dvalues)
                                      
!<description>
  ! This routine allows to evaluate a finite element solution vector
  ! rvectorScalar simultaneously in multiple points on multiple elements in a 
  ! discretisation.
  ! revalElementScalar must specify all information about where and how
  ! to evaluate; e.g. the coordinates of the evaluation points are
  ! to be fo9und here. The routine will then evaluate the element ieltyp
  ! in these points.
  !
  ! This routine is specialised to evaluate in multiple elements. For this
  ! purpose, the caller must make sure, that the same finite element type
  ! is used on all elements where to evaluate!
  ! So, evaluating 'simultaneously' on some $Q_1$ and some $P_1$ elements
  ! is not allowed e.g.. 
!</description>

!<input>
  ! Element evaluation set that contains all information necessary
  ! for the evaluation (coordinates of the points, transformation,...)
  TYPE(t_evalElementSet), INTENT(IN) :: revalElementSet

  ! The scalar solution vector that is to be evaluated.
  TYPE(t_vectorScalar), INTENT(IN)              :: rvectorScalar
  
  ! The FE function must be discretised with the same trial functions on all
  ! elements where it should be evaluated here. ieltyp defines the type
  ! of FE trial function that was used for the discretisation on those 
  ! elements that we are concerning here.
  INTEGER(I32), INTENT(IN)                      :: ieltyp

  ! An array accepting the DOF's on all elements in the trial space
  ! of the FE function.
  ! DIMENSION(\#local DOF's in trial space,nelements)
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), INTENT(IN) :: IdofsTrial
  
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
  INTEGER :: indofTrial,npoints,nelements
  REAL(DP) :: dval
  INTEGER :: iel,ipoint,ibas
  REAL(DP), DIMENSION(:), POINTER :: p_Ddata
  REAL(SP), DIMENSION(:), POINTER :: p_Fdata
  
  npoints = UBOUND(Dvalues,1)
  nelements = UBOUND(Dvalues,2)
  
  ! What to evaluate?
  Bder = .FALSE.
  Bder(iderType) = .TRUE.
  
  ! Allocate memory for the basis function values
  indofTrial = elem_igetNDofLoc(ieltyp)
  ALLOCATE(DbasTrial(indofTrial,elem_getMaxDerivative(ieltyp),npoints,nelements))
  
  ! Evaluate the basis functions
  CALL elem_generic_sim2 (ieltyp, revalElementSet, Bder, DbasTrial)
  
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
    CALL sys_halt()
  END IF
  
  ! Release memory, finish
  DEALLOCATE(DbasTrial)

  END SUBROUTINE

END MODULE
