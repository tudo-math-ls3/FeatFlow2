!##############################################################################
!# ****************************************************************************
!# <name> linearformevaluation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines for the discretisation of linear forms,
!# i.e. the creation of vectors (which usually appear as RHS vectors). 
!# It contains the
!# following set of routines:
!#
!# 1.) linf_buildVectorScalar
!#     -> Assembles the entries of a vector according to a linear form.
!# </purpose>
!##############################################################################

MODULE linearformevaluation

  USE fsystem
  USE linearsystemscalar
  USE spatialdiscretisation
  USE scalarpde
  USE derivatives
  USE cubature
  USE collection
  USE domainintegration
  
  IMPLICIT NONE

!<constants>
!<constantblock description="Constants defining the blocking of the assembly">

  ! Number of elements to handle simultaneously when building vectors
  INTEGER :: LINF_NELEMSIM   = 1000
  
!</constantblock>
!</constants>

CONTAINS

  !****************************************************************************

!<subroutine>

  SUBROUTINE linf_buildVectorScalar (rdiscretisation,rform,bclear,rvectorScalar,&
                                     fcoeff_buildVectorSc_sim,rcollection)
  
!<description>
  ! This routine assembles the entries of a vector according to a linear form
  ! (typically used for assembling RHS vectors).
  !
  ! If bclear=TRUE, the vector is cleared before the assembly and any 
  ! sorting of the entries is switched off - the vector is set up unsorted.
  !
  ! If bclear=FALSE, the vector must be unsorted when this routine is called, 
  ! otherwise an error is thrown.
!</description>

!<input>
  ! The underlying discretisation structure which is to be used to
  ! create the vector.
  TYPE(t_spatialDiscretisation), INTENT(IN), TARGET :: rdiscretisation
  
  ! The linear form specifying the underlying PDE of the discretisation.
  TYPE(t_linearForm), INTENT(IN) :: rform
  
  ! Whether to clear the vector before calculating the entries.
  ! If .FALSE., the new entries are added to the existing entries.
  LOGICAL, INTENT(IN) :: bclear
  
  ! OPTIONAL: A collection structure. This structure is 
  ! given to the callback function for calculating the function
  ! which should be discretised in the linear form.
  TYPE(t_collection), INTENT(INOUT), TARGET, OPTIONAL :: rcollection
  
  ! A callback routine for the function to be discretised.
  INCLUDE 'intf_coefficientVectorSc.inc'
  OPTIONAL :: fcoeff_buildVectorSc_sim
!</input>

!<inputoutput>
  ! The FE vector. Calculated entries are imposed to this vector.
  TYPE(t_vectorScalar), INTENT(INOUT) :: rvectorScalar
!</inputoutput>

!</subroutine>
  
  ! If the vector is not set up as new vector, it has to be unsorted.
  ! If it's a new vector, we switch off the sorting.
  IF (bclear) THEN
    rvectorScalar%isortStrategy = -ABS(rvectorScalar%isortStrategy)
  END IF
  
  ! The vector must be unsorted, otherwise we can't set up the vector.
  IF (rvectorScalar%isortStrategy .GT. 0) THEN
    PRINT *,'linf_buildVectorScalar: Vector must be unsorted!'
    CALL sys_halt()
  END IF

  ! Do we have a uniform triangulation? Would simplify a lot...
  IF (rdiscretisation%ccomplexity .EQ. SPDISC_UNIFORM) THEN 
  
    IF (rvectorScalar%cdataType .EQ. ST_DOUBLE) THEN
  
      CALL linf_buildVectord_conf2 (rdiscretisation,rform,bclear,rVectorScalar,&  
                                   fcoeff_buildVectorSc_sim,rcollection)
    ELSE
      PRINT *,'linf_buildVectorScalar: Single precision vectors currently not supported!'
    END IF
  
  ! Do we have a uniform triangulation? Would simplify a lot...
  ELSE IF (rdiscretisation%ccomplexity .EQ. SPDISC_CONFORMAL) THEN 
  
    IF (rvectorScalar%cdataType .EQ. ST_DOUBLE) THEN
  
      CALL linf_buildVectord_conf2 (rdiscretisation,rform,bclear,rVectorScalar,&  
                                   fcoeff_buildVectorSc_sim,rcollection)
    ELSE
      PRINT *,'linf_buildVectorScalar: Single precision vectors currently not supported!'
    END IF
  
  ELSE
    PRINT *,'linf_buildVectorScalar: General discretisation &
            & not implemented!'
    CALL sys_halt()
  END IF

  END SUBROUTINE
  
  !****************************************************************************

!<subroutine>

  SUBROUTINE linf_buildVectord_conf (rdiscretisation,rform,bclear,rVectorScalar,&
                                     fcoeff_buildVectorSc_sim,rcollection)

!<description>
  ! This routine calculates the entries of a discretised finite element vector.
  ! The discretisation is assumed to be conformal, i.e. the DOF's
  ! of all finite elements must 'match'. 
  ! The linear form is defined by
  !        (f,$phi_i$), i=1..*
  ! with $Phi_i$ being the test functions defined in the discretisation
  ! structure.
  ! In case the array for the vector entries does not exist, the routine
  ! allocates memory in size of the matrix of the heap for the matrix entries
  ! and initialises all necessary variables of the vector according to the
  ! parameters (NEQ, pointer to the discretisation,...)
  !
  ! Double-precision version.
!</description>

!<input>
  ! The underlying discretisation structure which is to be used to
  ! create the vector.
  TYPE(t_spatialDiscretisation), INTENT(IN), TARGET :: rdiscretisation
  
  ! The linear form specifying the underlying PDE of the discretisation.
  TYPE(t_linearForm), INTENT(IN) :: rform
  
  ! Whether to clear the matrix before calculating the entries.
  ! If .FALSE., the new matrix entries are added to the existing entries.
  LOGICAL, INTENT(IN) :: bclear
  
  ! OPTIONAL: A collection structure. This structure is given to the
  ! callback function for nonconstant coefficients to provide additional
  ! information. 
  TYPE(t_collection), INTENT(INOUT), TARGET, OPTIONAL :: rcollection
  
  ! A callback routine which is able to calculate the values of the
  ! function $f$ which is to be discretised.
  INCLUDE 'intf_coefficientVectorSc2.inc'
  OPTIONAL :: fcoeff_buildVectorSc_sim
!</input>

!<inputoutput>
  ! The FE vector. Calculated matrix entries are added to this vector.
  TYPE(t_vectorScalar), INTENT(INOUT) :: rvectorScalar
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i,i1,j,icurrentElementDistr, ICUBP, IALBET, IA, NVE, IDOFE
  LOGICAL :: bnonparTest
  INTEGER(PREC_ELEMENTIDX) :: IEL, IELmax, IELset
  REAL(DP) :: OM,AUX
  
  ! Array to tell the element which derivatives to calculate
  LOGICAL, DIMENSION(EL_MAXNDER) :: Bder
  
  ! Cubature point coordinates on the reference element
  REAL(DP), DIMENSION(CUB_MAXCUBP, NDIM3D) :: Dxi

  ! For every cubature point on the reference element,
  ! the corresponding cubature weight
  REAL(DP), DIMENSION(CUB_MAXCUBP) :: Domega
  
  ! number of cubature points on the reference element
  INTEGER :: ncubp
  
  ! Pointer to the vector entries
  REAL(DP), DIMENSION(:), POINTER :: p_Ddata

  ! An allocateable array accepting the DOF's of a set of elements.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsTest
  !INTEGER(PREC_DOFIDX), DIMENSION(EL_MAXNBAS,BILF_NELEMSIM), TARGET :: IdofsTest, IdofsTrial
  !INTEGER(PREC_DOFIDX), DIMENSION(:,:), POINTER :: p_IdofsTrial
  
  ! Allocateable arrays for the values of the basis functions - 
  ! for test space.
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: DbasTest
  
  ! Number of entries in the vector - for quicker access
  INTEGER(I32) :: NEQ
  
  ! Number of local degees of freedom for test functions
  INTEGER :: indofTest
  
  ! The triangulation structure - to shorten some things...
  TYPE(t_triangulation), POINTER :: p_rtriangulation
  
  ! A pointer to an element-number list
  INTEGER(I32), DIMENSION(:), POINTER :: p_IelementList
  
  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: DcubPtsRef

  ! An array receiving the coordinates of cubature points on
  ! the real element for all elements in a set.
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: DcubPtsReal

  ! Pointer to the point coordinates to pass to the element function.
  ! Point either to DcubPtsRef or to DcubPtsReal, depending on whether
  ! the trial/test element is parametric or not.
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPtsTest
  
  ! Array with coordinates of the corners that form the real element.
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Dcoords
  
  ! A small vector holding only the additive contributions of
  ! one element
  REAL(DP), DIMENSION(EL_MAXNBAS) :: DlocalData
  
  ! Arrays for saving Jacobian determinants and matrices
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Ddetj
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Djac
  
  ! Pointer to KVERT of the triangulation
  INTEGER(I32), DIMENSION(:,:), POINTER :: p_IverticesAtElement
  
  ! Pointer to DCORVG of the triangulation
  REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoords
  
  ! Current element distribution
  TYPE(t_elementDistribution), POINTER :: p_elementDistribution
  
  ! Number of elements in the current element distribution
  INTEGER(PREC_ELEMENTIDX) :: NEL
  
  ! Number of elements in a block. Normally =BILF_NELEMSIM,
  ! except if there are less elements in the discretisation.
  INTEGER :: nelementsPerBlock
  
  ! Pointer to the coefficients that are computed by the callback routine.
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Dcoefficients
  
  !REAL(DP), DIMENSION(11) :: DT
  
  !CHARACTER(LEN=20) :: CFILE
  
  ! Which derivatives of basis functions are needed?
  ! Check the descriptors of the bilinear form and set BDER
  ! according to these.

  !CALL ZTIME(DT(1))

  Bder = .FALSE.
  
  ! Loop through the additive terms
  DO i=1,rform%itermCount
    ! The desriptor Idescriptors gives directly the derivative
    ! which is to be computed!
    I1=rform%Idescriptors(i)
    
    IF ((I1 .LE.0) .OR. (I1 .GT. DER_MAXNDER)) THEN
      PRINT *,'linf_buildVectord_conf: Invalid descriptor'
      CALL sys_halt()
    ENDIF
    
    Bder(I1)=.TRUE.
  END DO
  
  IF (rvectorScalar%h_Ddata .EQ. ST_NOHANDLE) THEN
  
    ! Get the size of the vector and put it to the matrix structure.
    NEQ = dof_igetNDofGlob(rdiscretisation)
    
    ! Initialise the vector parameters
    rvectorScalar%NEQ            = NEQ
    rvectorScalar%iidxFirstEntry = 1
    rvectorScalar%p_rspatialDiscretisation => rdiscretisation
    rvectorScalar%cdataType      = ST_DOUBLE

    ! Clear the entries in the vector - we need to start with zero
    ! when assembling a new vector.
    CALL storage_new1D ('linf_buildVectord_conf', 'vector', &
                        NEQ, ST_DOUBLE, rvectorScalar%h_Ddata, &
                        ST_NEWBLOCK_ZERO)
    CALL storage_getbase_double (rvectorScalar%h_Ddata,p_Ddata)

  ELSE
  
    ! Get information about the vector:
    NEQ = rvectorScalar%NEQ
  
    CALL storage_getbase_double (rvectorScalar%h_Ddata,p_Ddata)
    
    ! Maybe the vector is a partial vector of a larger one.
    ! Let the pointer point to the right position.
    IF (rvectorScalar%iidxFirstEntry .NE. 1) THEN
      p_Ddata => p_Ddata (rvectorScalar%iidxFirstEntry : &
                          rvectorScalar%iidxFirstEntry + rvectorScalar%NEQ - 1)
    END IF

    ! If desired, clear the vector before assembling.
    IF (bclear) THEN
      CALL lalg_clearVectorDble (p_Ddata)
    END IF
    
  END IF
  
  ! Get a pointer to the triangulation - for easier access.
  p_rtriangulation => rdiscretisation%p_rtriangulation
  
  ! For saving some memory in smaller discretisations, we calculate
  ! the number of elements per block. For smaller triangulations,
  ! this is NEL. If there are too many elements, it's at most
  ! BILF_NELEMSIM. This is only used for allocating some arrays.
  nelementsPerBlock = MIN(LINF_NELEMSIM,p_rtriangulation%NEL)
  
  ! Get a pointer to the KVERT and DCORVG array
  CALL storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement, &
                             p_IverticesAtElement)
  CALL storage_getbase_double2D(p_rtriangulation%h_DvertexCoords, &
                             p_DvertexCoords)

  ! Allocate memory for corner coordinates
  ALLOCATE(DCoords(2,TRIA_MAXNVE2D,nelementsPerBlock))
  
  ! Now loop over the different element distributions (=combinations
  ! of trial and test functions) in the discretisation.
  !CALL ZTIME(DT(2))

  DO icurrentElementDistr = 1,rdiscretisation%inumFESpaces
  
    ! Activate the current element distribution
    p_elementDistribution => rdiscretisation%RelementDistribution(icurrentElementDistr)
  
    ! Cancel if this element distribution is empty.
    IF (p_elementDistribution%NEL .EQ. 0) CYCLE

    ! Get the number of corner vertices of the element
    NVE = elem_igetNVE(p_elementDistribution%itrialElement)
    IF (NVE .NE. p_elementDistribution%itestElement) THEN
      PRINT *,'bilf_buildMatrix9d_conf2: element spaces incompatible!'
      CALL sys_halt()
    END IF
    
    ! Get the number of local DOF's for trial and test functions
    indofTest = elem_igetNDofLoc(p_elementDistribution%itestElement)
    
    ! Initialise the cubature formula,
    ! Get cubature weights and point coordinates on the reference element
    CALL cub_getCubPoints(p_elementDistribution%ccubTypeLinForm, ncubp, Dxi, Domega)
    
    ! Allocate arrays accepting cubature point coordinates.
    ! It's at most as large as number of elements or length
    ! of the element set.
    ALLOCATE(DcubPtsRef(p_rtriangulation%ndim,ncubp,nelementsPerBlock))
    ALLOCATE(DcubPtsReal(p_rtriangulation%ndim,ncubp,nelementsPerBlock))
    
    ! Put the cubature point coordinates in the right format to the
    ! cubature-point array.
    ! Initialise all entries in DcubPtsRef with the same coordinates -
    ! as the cubature point coordinates are identical on all elements
    DO j=1,SIZE(DcubPtsRef,3)
      DO i=1,ncubp
        DcubPtsRef(1,i,j) = Dxi(i,1)
        DcubPtsRef(2,i,j) = Dxi(i,2)
      END DO
    END DO
    
    ! Allocate an array saving the coordinates of corner vertices of elements
    ALLOCATE(Djac(4,ncubp,nelementsPerBlock))
    ALLOCATE(Ddetj(ncubp,nelementsPerBlock))
    
    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  ALLOCATE(DbasTest(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    !  ALLOCATE(DbasTrial(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    ! would lead to nonused memory blocks in these arrays during the assembly, 
    ! which reduces the speed by 50%!
    
    ALLOCATE(DbasTest(indofTest,elem_getMaxDerivative(p_elementDistribution%itestElement),&
             ncubp,nelementsPerBlock))

    ! Allocate memory for the DOF's of all the elements.
    ALLOCATE(IdofsTest(indofTest,nelementsPerBlock))

    ! Allocate memory for the coefficients
    ALLOCATE(Dcoefficients(rform%itermCount,ncubp,nelementsPerBlock))
  
    ! Check if one of the trial/test elements is nonparametric
    bnonparTest  = elem_isNonparametric(p_elementDistribution%itestElement)
                    
    ! Let p_DcubPtsTest point either to DcubPtsReal or
    ! DcubPtsRef - depending on whether the space is parametric or not.
    IF (bnonparTest) THEN
      p_DcubPtsTest => DcubPtsReal
    ELSE
      p_DcubPtsTest => DcubPtsRef
    END IF
    
    !CALL ZTIME(DT(3))
    ! p_IelementList must point to our set of elements in the discretisation
    ! with that combination of trial/test functions
    CALL storage_getbase_int (p_elementDistribution%h_IelementList, &
                              p_IelementList)

    ! Get the number of elements there.
    NEL = p_elementDistribution%NEL
                              
    ! Loop over the elements - blockwise.
    DO IELset = 1, NEL, LINF_NELEMSIM
    
      ! We always handle LINF_NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
      ! elements simultaneously.
      
      IELmax = MIN(NEL,IELset-1+LINF_NELEMSIM)
    
      ! Calculate the global DOF's into IdofsTest.
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF's of our LINF_NELEMSIM elements simultaneously.
      CALL dof_locGlobMapping_mult(rdiscretisation, p_IelementList(IELset:IELmax), &
                                  .TRUE.,IdofsTest)
                                   
      !CALL ZTIME(DT(4))
      
      ! We have the coordinates of the cubature points saved in the
      ! coordinate array from above. Unfortunately for nonparametric
      ! elements, we need the real coordinate.
      ! Furthermore, we anyway need the coordinates of the element
      ! corners and the Jacobian determinants corresponding to
      ! all the points.
      !
      ! At first, get the coordinates of the corners of all the
      ! elements in the current set. 
      
!      DO IEL=1,IELmax-IELset+1
!        DCoords(:,:,IEL) = p_DvertexCoords(:, &
!                            p_IverticesAtElement(:,p_IelementList(IELset+IEL-1)))
!      END DO
      DO IEL=1,IELmax-IELset+1
        DO J = 1,NVE
          DO I = 1,p_rtriangulation%ndim
            DCoords(I,J,IEL) = p_DvertexCoords(I, &
                               p_IverticesAtElement(J,p_IelementList(IELset+IEL-1)))
          END DO
        END DO
      END DO
      !CALL ZTIME(DT(6))
      
      ! Depending on the type of transformation, we must now choose
      ! the mapping between the reference and the real element.
      ! In case we use a nonparametric element as test function, we need the 
      ! coordinates of the points on the real element, too.
      ! Unfortunately, we need the real coordinates of the cubature points
      ! anyway for the function - so calculate them all.
      CALL trafo_calctrafo_sim (&
            rdiscretisation%RelementDistribution(icurrentElementDistr)%ctrafoType,&
            IELmax-IELset+1,ncubp,Dcoords,&
            DcubPtsRef,Djac(:,:,1:IELmax-IELset+1),Ddetj(:,1:IELmax-IELset+1),DcubPtsReal)
    
      !CALL ZTIME(DT(7))
      
      ! Now it's time to call our coefficient function to calculate the
      ! function values in the cubature points:
      
      CALL fcoeff_buildVectorSc_sim (rdiscretisation,icurrentElementDistr, rform, &
                IELset,IELmax-IELset+1_I32,ncubp,p_IelementList(IELset:IELmax),Dcoords, &
                DcubPtsRef,DcubPtsReal,IdofsTest,Djac,Ddetj, &
                Dcoefficients(:,:,1:IELmax-IELset+1_I32),rcollection)
      
      !CALL ZTIME(DT(8))                              
      ! Calculate the values of the basis functions.
      ! Pass p_DcubPts as point coordinates, which point either to the
      ! coordinates on the reference element (the same for all elements)
      ! or on the real element - depending on whether this is a 
      ! parametric or nonparametric element.
      CALL elem_generic_sim (p_elementDistribution%itestElement, Dcoords, &
            Djac(:,:,1:IELmax-IELset+1), Ddetj(:,1:IELmax-IELset+1), &
            Bder, DbasTest, ncubp, IELmax-IELset+1, p_DcubPtsTest)
            
      !CALL ZTIME(DT(9))
      ! Values of all basis functions calculated. Now we can start 
      ! to integrate!
      !
      ! Loop through elements in the set and for each element,
      ! loop through the DOF's and cubature points to calculate the
      ! integral:
      
      DO IEL=1,IELmax-IELset+1
        
        ! We make a 'local' approach, i.e. we calculate the values of the
        ! integral into the vector DlocalData and add them later into
        ! the large solution vector.
        
        ! Clear the output vector.
        DlocalData(1:indofTest) = 0.0_DP

        ! Loop over all cubature points on the current element
        DO ICUBP = 1, ncubp
        
          ! calculate the current weighting factor in the cubature formula
          ! in that cubature point.
          !
          ! Take the absolut value of the determinant of the mapping.
          ! In 2D, the determinant is always positive, whereas in 3D,
          ! the determinant might be negative -- that's normal!

          OM = Domega(ICUBP)*ABS(Ddetj(ICUBP,IEL))

          ! Loop over the additive factors in the linear form.
          DO IALBET = 1,rform%itermcount
          
            ! Get from Idescriptors the type of the derivatives for the 
            ! test and trial functions. The summand we calculate
            ! here will be:
            !
            ! int_...  f * ( phi_i )_IA
            !
            ! -> IA=0: function value, 
            !      =1: first derivative, 
            !      =2: 2nd derivative,...
            !    as defined in the module 'derivative'.
            
            IA = rform%Idescriptors(IALBET)
            
            ! Multiply OM with the coefficient of the form.
            ! This gives the actual value to multiply the
            ! function value with before summing up to the integral.
            ! Get the precalculated coefficient from the coefficient array.
            AUX = OM * Dcoefficients(IALBET,ICUBP,IEL)
          
            ! Now loop through all possible combinations of DOF's
            ! in the current cubature point. 

            DO IDOFE=1,indofTest
            
              ! Get the value of the basis function 
              ! phi_o in the cubature point. 
              ! Them multiply:
              !    DBAS(..) * AUX
              ! ~= phi_i * coefficient * cub.weight
              ! Summing this up gives the integral, so the contribution
              ! to the vector. 
              !
              ! Simply summing up DBAS(..) * AUX would give
              ! the additive contribution for the vector. We save this
              ! contribution in the local array.

              DlocalData(IDOFE) = DlocalData(IDOFE)+DbasTest(IDOFE,IA,ICUBP,IEL)*AUX
            
            END DO ! IDOFE
            
          END DO ! IALBET

        END DO ! ICUBP 

        ! Incorporate the local vector into the global one.
        ! The 'local' DOF 1..indofTest is mapped to the global DOF using
        ! the IdofsTest array.
        DO IDOFE=1,indofTest
          p_Ddata(IdofsTest(IDOFE,IEL)) = p_Ddata(IdofsTest(IDOFE,IEL)) + DlocalData(IDOFE)
        END DO

      END DO ! IEL

      !CALL ZTIME(DT(10))
    END DO ! IELset
    
    DEALLOCATE(Dcoefficients)
    DEALLOCATE(IdofsTest)
    DEALLOCATE(DbasTest)
    DEALLOCATE(Ddetj)
    DEALLOCATE(Djac)
    DEALLOCATE(DcubPtsReal)
    DEALLOCATE(DcubPtsRef)

  END DO ! icurrentElementDistr

  ! Clean up memory, finish

  DEALLOCATE(Dcoords)
  !CALL ZTIME(DT(11))
  
  !DO i=2,11
  !  PRINT *,'Time for assembly part ',i,': ',DT(i)-DT(i-1)
  !END DO
  
  END SUBROUTINE
  
  !****************************************************************************

!<subroutine>

  SUBROUTINE linf_buildVectord_conf2 (rdiscretisation,rform,bclear,rVectorScalar,&
                                     fcoeff_buildVectorSc_sim,rcollection)

!<description>
  ! This routine calculates the entries of a discretised finite element vector.
  ! The discretisation is assumed to be conformal, i.e. the DOF's
  ! of all finite elements must 'match'. 
  ! The linear form is defined by
  !        (f,$phi_i$), i=1..*
  ! with $Phi_i$ being the test functions defined in the discretisation
  ! structure.
  ! In case the array for the vector entries does not exist, the routine
  ! allocates memory in size of the matrix of the heap for the matrix entries
  ! and initialises all necessary variables of the vector according to the
  ! parameters (NEQ, pointer to the discretisation,...)
  !
  ! Double-precision version.
!</description>

!<input>
  ! The underlying discretisation structure which is to be used to
  ! create the vector.
  TYPE(t_spatialDiscretisation), INTENT(IN), TARGET :: rdiscretisation
  
  ! The linear form specifying the underlying PDE of the discretisation.
  TYPE(t_linearForm), INTENT(IN) :: rform
  
  ! Whether to clear the matrix before calculating the entries.
  ! If .FALSE., the new matrix entries are added to the existing entries.
  LOGICAL, INTENT(IN) :: bclear
  
  ! OPTIONAL: A pointer to a collection structure. This structure is given to the
  ! callback function for nonconstant coefficients to provide additional
  ! information. 
  TYPE(t_collection), INTENT(INOUT), TARGET, OPTIONAL :: rcollection
  
  ! A callback routine which is able to calculate the values of the
  ! function $f$ which is to be discretised.
  INCLUDE 'intf_coefficientVectorSc.inc'
  OPTIONAL :: fcoeff_buildVectorSc_sim
!</input>

!<inputoutput>
  ! The FE vector. Calculated matrix entries are added to this vector.
  TYPE(t_vectorScalar), INTENT(INOUT) :: rvectorScalar
!</inputoutput>

!</subroutine>

  ! local variables
  INTEGER :: i,i1,j,k,icurrentElementDistr, ICUBP, IALBET, IA, NVE
  LOGICAL :: bnonparTest
  INTEGER(I32) :: IEL, IELmax, IELset, IDOFE
  REAL(DP) :: OM,AUX
  
  ! Array to tell the element which derivatives to calculate
  LOGICAL, DIMENSION(EL_MAXNDER) :: Bder
  
  ! Cubature point coordinates on the reference element
  REAL(DP), DIMENSION(CUB_MAXCUBP, NDIM3D) :: Dxi

  ! For every cubature point on the reference element,
  ! the corresponding cubature weight
  REAL(DP), DIMENSION(CUB_MAXCUBP) :: Domega
  
  ! number of cubature points on the reference element
  INTEGER :: ncubp
  
  ! Pointer to the vector entries
  REAL(DP), DIMENSION(:), POINTER :: p_Ddata

  ! An allocateable array accepting the DOF's of a set of elements.
  INTEGER(PREC_DOFIDX), DIMENSION(:,:), ALLOCATABLE, TARGET :: IdofsTest
  !INTEGER(PREC_DOFIDX), DIMENSION(EL_MAXNBAS,BILF_NELEMSIM), TARGET :: IdofsTest, IdofsTrial
  !INTEGER(PREC_DOFIDX), DIMENSION(:,:), POINTER :: p_IdofsTrial
  
  ! Allocateable arrays for the values of the basis functions - 
  ! for test space.
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: DbasTest
  
  ! Number of entries in the vector - for quicker access
  INTEGER(I32) :: NEQ
  
  ! Number of local degees of freedom for test functions
  INTEGER :: indofTest
  
  ! The triangulation structure - to shorten some things...
  TYPE(t_triangulation), POINTER :: p_rtriangulation
  
  ! A pointer to an element-number list
  INTEGER(I32), DIMENSION(:), POINTER :: p_IelementList
  
  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPtsRef

  ! An array receiving the coordinates of cubature points on
  ! the real element for all elements in a set.
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPtsReal

  ! Pointer to the point coordinates to pass to the element function.
  ! Point either to p_DcubPtsRef or to p_DcubPtsReal, depending on whether
  ! the trial/test element is parametric or not.
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_DcubPtsTest
  
  ! Array with coordinates of the corners that form the real element.
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_Dcoords
  
  ! A small vector holding only the additive controbutions of
  ! one element
  REAL(DP), DIMENSION(EL_MAXNBAS) :: DlocalData
  
  ! Arrays for saving Jacobian determinants and matrices
  REAL(DP), DIMENSION(:,:), POINTER :: p_Ddetj
  REAL(DP), DIMENSION(:,:,:), POINTER :: p_Djac
  
  ! Pointer to KVERT of the triangulation
  INTEGER(I32), DIMENSION(:,:), POINTER :: p_IverticesAtElement
  
  ! Pointer to DCORVG of the triangulation
  REAL(DP), DIMENSION(:,:), POINTER :: p_DvertexCoords
  
  ! Current element distribution
  TYPE(t_elementDistribution), POINTER :: p_elementDistribution
  
  ! Number of elements in the current element distribution
  INTEGER(PREC_ELEMENTIDX) :: NEL

  ! Number of elements in a block. Normally =BILF_NELEMSIM,
  ! except if there are less elements in the discretisation.
  INTEGER :: nelementsPerBlock
  
  ! Pointer to the coefficients that are computed by the callback routine.
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: Dcoefficients
  
  ! A t_domainIntSubset structure that is used for storing information
  ! and passing it to callback routines.
  TYPE(t_domainIntSubset) :: rintSubset

  !REAL(DP), DIMENSION(11) :: DT
  
  !CHARACTER(LEN=20) :: CFILE
  
  ! Which derivatives of basis functions are needed?
  ! Check the descriptors of the bilinear form and set BDER
  ! according to these.

  !CALL ZTIME(DT(1))

  Bder = .FALSE.
  
  ! Loop through the additive terms
  DO i=1,rform%itermCount
    ! The desriptor Idescriptors gives directly the derivative
    ! which is to be computed!
    I1=rform%Idescriptors(i)
    
    IF ((I1 .LE.0) .OR. (I1 .GT. DER_MAXNDER)) THEN
      PRINT *,'linf_buildVectord_conf: Invalid descriptor'
      CALL sys_halt()
    ENDIF
    
    Bder(I1)=.TRUE.
  END DO
  
  IF (rvectorScalar%h_Ddata .EQ. ST_NOHANDLE) THEN
  
    ! Get the size of the vector and put it to the matrix structure.
    NEQ = dof_igetNDofGlob(rdiscretisation)
    
    ! Initialise the vector parameters
    rvectorScalar%NEQ            = NEQ
    rvectorScalar%iidxFirstEntry = 1
    rvectorScalar%p_rspatialDiscretisation => rdiscretisation
    rvectorScalar%cdataType      = ST_DOUBLE

    ! Clear the entries in the vector - we need to start with zero
    ! when assembling a new vector.
    CALL storage_new1D ('linf_buildVectord_conf', 'vector', &
                        NEQ, ST_DOUBLE, rvectorScalar%h_Ddata, &
                        ST_NEWBLOCK_ZERO)
    CALL storage_getbase_double (rvectorScalar%h_Ddata,p_Ddata)

  ELSE
  
    ! Get information about the vector:
    NEQ = rvectorScalar%NEQ
  
    CALL storage_getbase_double (rvectorScalar%h_Ddata,p_Ddata)
    
    ! Maybe the vector is a partial vector of a larger one.
    ! Let the pointer point to the right position.
    IF (rvectorScalar%iidxFirstEntry .NE. 1) THEN
      p_Ddata => p_Ddata (rvectorScalar%iidxFirstEntry : &
                          rvectorScalar%iidxFirstEntry + rvectorScalar%NEQ - 1)
    END IF

    ! If desired, clear the vector before assembling.
    IF (bclear) THEN
      CALL lalg_clearVectorDble (p_Ddata)
    END IF
    
  END IF
  
  ! Get a pointer to the triangulation - for easier access.
  p_rtriangulation => rdiscretisation%p_rtriangulation
  
  ! For saving some memory in smaller discretisations, we calculate
  ! the number of elements per block. For smaller triangulations,
  ! this is NEL. If there are too many elements, it's at most
  ! BILF_NELEMSIM. This is only used for allocating some arrays.
  nelementsPerBlock = MIN(LINF_NELEMSIM,p_rtriangulation%NEL)
  
  ! Get a pointer to the KVERT and DCORVG array
  CALL storage_getbase_int2D(p_rtriangulation%h_IverticesAtElement, &
                             p_IverticesAtElement)
  CALL storage_getbase_double2D(p_rtriangulation%h_DvertexCoords, &
                             p_DvertexCoords)

  ! Now loop over the different element distributions (=combinations
  ! of trial and test functions) in the discretisation.
  !CALL ZTIME(DT(2))

  DO icurrentElementDistr = 1,rdiscretisation%inumFESpaces
  
    ! Activate the current element distribution
    p_elementDistribution => rdiscretisation%RelementDistribution(icurrentElementDistr)
  
    ! Cancel if this element distribution is empty.
    IF (p_elementDistribution%NEL .EQ. 0) CYCLE

    ! Get the number of local DOF's for trial and test functions
    indofTest = elem_igetNDofLoc(p_elementDistribution%itestElement)
    
    ! Get the number of corner vertices of the element
    NVE = elem_igetNVE(p_elementDistribution%itrialElement)
    IF (NVE .NE. elem_igetNVE(p_elementDistribution%itestElement)) THEN
      PRINT *,'linf_buildVectord_conf2: element spaces incompatible!'
      CALL sys_halt()
    END IF
    
    ! Initialise the cubature formula,
    ! Get cubature weights and point coordinates on the reference element
    CALL cub_getCubPoints(p_elementDistribution%ccubTypeLinForm, ncubp, Dxi, Domega)
    
    ! Open-MP-Extension: Open threads here.
    ! "j" is declared as private; shared gave errors with the Intel compiler
    ! in Windows!?!
    ! Each thread will allocate its own local memory...
    !
    !$OMP PARALLEL PRIVATE(rintSubset, p_DcubPtsRef,p_DcubPtsReal, &
    !$OMP   p_Djac,p_Ddetj,p_Dcoords,DbasTest, &
    !$OMP   IdofsTest,bnonparTest,&
    !$OMP   p_DcubPtsTest,Dcoefficients, &
    !$OMP   j, ielmax,IEL, idofe, &
    !$OMP   ICUBP, IALBET,OM,IA,aux)    
    
    ! Get from the trial element space the type of coordinate system
    ! that is used there:
    j = elem_igetCoordSystem(p_elementDistribution%itrialElement)
    
    ! Allocate memory and get local references to it.
    CALL domint_initIntegration (rintSubset,nelementsPerBlock,ncubp,&
        j,p_rtriangulation%ndim,NVE)
    p_DcubPtsRef =>  rintSubset%p_DcubPtsRef
    p_DcubPtsReal => rintSubset%p_DcubPtsReal
    p_Djac =>        rintSubset%p_Djac
    p_Ddetj =>       rintSubset%p_Ddetj
    p_Dcoords =>     rintSubset%p_DCoords

    ! Put the cubature point coordinates in the right format to the
    ! cubature-point array.
    ! Initialise all entries in p_DcubPtsRef with the same coordinates -
    ! as the cubature point coordinates are identical on all elements
    DO j=1,SIZE(p_DcubPtsRef,3)
      DO i=1,ncubp
        DO k=1,SIZE(p_DcubPtsRef,1)
          ! Could be solved using the TRANSPOSE operator - but often is's 
          ! faster this way...
          p_DcubPtsRef(k,i,j) = Dxi(i,k)
        END DO
      END DO
    END DO
    
    ! Quickly check if one of the specified derivatives is out of the allowed range:
    DO IALBET = 1,rform%itermcount
      IA = rform%Idescriptors(IALBET)
      IF ((IA.LT.0) .OR. &
          (IA .GT. elem_getMaxDerivative(p_elementDistribution%itrialElement))) THEN
        PRINT *,'linf_buildVectord_conf2: Specified test-derivative',IA,&
                ' not available'
        CALL sys_halt()
      END IF
    END DO

    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  ALLOCATE(DbasTest(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    !  ALLOCATE(DbasTrial(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    ! would lead to nonused memory blocks in these arrays during the assembly, 
    ! which reduces the speed by 50%!
    ALLOCATE(DbasTest(indofTest,elem_getMaxDerivative(p_elementDistribution%itestElement),&
             ncubp,nelementsPerBlock))

    ! Allocate memory for the DOF's of all the elements.
    ALLOCATE(IdofsTest(indofTest,nelementsPerBlock))

    ! Allocate memory for the coefficients
    ALLOCATE(Dcoefficients(rform%itermCount,ncubp,nelementsPerBlock))
  
    ! Check if one of the trial/test elements is nonparametric
    bnonparTest  = elem_isNonparametric(p_elementDistribution%itestElement)
                    
    ! Let p_DcubPtsTest point either to p_DcubPtsReal or
    ! p_DcubPtsRef - depending on whether the space is parametric or not.
    IF (bnonparTest) THEN
      p_DcubPtsTest => p_DcubPtsReal
    ELSE
      p_DcubPtsTest => p_DcubPtsRef
    END IF
    
    !CALL ZTIME(DT(3))
    ! p_IelementList must point to our set of elements in the discretisation
    ! with that combination of trial/test functions
    CALL storage_getbase_int (p_elementDistribution%h_IelementList, &
                              p_IelementList)
                              
    ! Get the number of elements there.
    NEL = p_elementDistribution%NEL
  
  
    ! Loop over the elements - blockwise.
    !$OMP do schedule(static,1)
    DO IELset = 1, NEL, LINF_NELEMSIM
    
      ! We always handle LINF_NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most LINF_NELEMSIM
      ! elements simultaneously.
      
      IELmax = MIN(NEL,IELset-1+LINF_NELEMSIM)
    
      ! Calculate the global DOF's into IdofsTest.
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF's of our LINF_NELEMSIM elements simultaneously.
      CALL dof_locGlobMapping_mult(rdiscretisation, p_IelementList(IELset:IELmax), &
                                  .TRUE.,IdofsTest)
                                   
      !CALL ZTIME(DT(4))
      
      ! We have the coordinates of the cubature points saved in the
      ! coordinate array from above. Unfortunately for nonparametric
      ! elements, we need the real coordinate.
      ! Furthermore, we anyway need the coordinates of the element
      ! corners and the Jacobian determinants corresponding to
      ! all the points.
      !
      ! At first, get the coordinates of the corners of all the
      ! elements in the current set. 
      
!      DO IEL=1,IELmax-IELset+1
!        p_Dcoords(:,:,IEL) = p_DvertexCoords(:, &
!                            p_IverticesAtElement(:,p_IelementList(IELset+IEL-1)))
!      END DO
!      DO IEL=1,IELmax-IELset+1
!        DO J = 1,NVE
!          DO I = 1,p_rtriangulation%ndim
!            p_Dcoords(I,J,IEL) = p_DvertexCoords(I, &
!                               p_IverticesAtElement(J,p_IelementList(IELset+IEL-1)))
!          END DO
!        END DO
!      END DO

      CALL trafo_getCoords_sim (elem_igetTrafoType(p_elementDistribution%itrialElement),&
          p_rtriangulation,p_IelementList(IELset:IELmax),p_Dcoords)

      !CALL ZTIME(DT(6))
      
      ! Depending on the type of transformation, we must now choose
      ! the mapping between the reference and the real element.
      ! In case we use a nonparametric element as test function, we need the 
      ! coordinates of the points on the real element, too.
      ! Unfortunately, we need the real coordinates of the cubature points
      ! anyway for the function - so calculate them all.
      CALL trafo_calctrafo_sim (&
            p_elementDistribution%ctrafoType,&
            IELmax-IELset+1,ncubp,p_Dcoords,&
            p_DcubPtsRef,p_Djac(:,:,1:IELmax-IELset+1),p_Ddetj(:,1:IELmax-IELset+1),&
            p_DcubPtsReal)
    
      !CALL ZTIME(DT(7))
      
      ! Now it's time to call our coefficient function to calculate the
      ! function values in the cubature points:
      
      rintSubset%ielementDistribution = icurrentElementDistr
      rintSubset%ielementStartIdx = IELset
      rintSubset%p_Ielements => p_IelementList(IELset:IELmax)
      CALL fcoeff_buildVectorSc_sim (rdiscretisation,rform, &
                IELmax-IELset+1_I32,ncubp,p_DcubPtsReal, &
                IdofsTest,rintSubset, &
                Dcoefficients(:,:,1:IELmax-IELset+1_I32),rcollection)
      
      !CALL ZTIME(DT(8))                              
      ! Calculate the values of the basis functions.
      ! Pass p_DcubPts as point coordinates, which point either to the
      ! coordinates on the reference element (the same for all elements)
      ! or on the real element - depending on whether this is a 
      ! parametric or nonparametric element.
      CALL elem_generic_sim (p_elementDistribution%itestElement, p_Dcoords, &
            p_Djac(:,:,1:IELmax-IELset+1), p_Ddetj(:,1:IELmax-IELset+1), &
            Bder, DbasTest, ncubp, IELmax-IELset+1, p_DcubPtsTest)
            
      !CALL ZTIME(DT(9))
      ! Values of all basis functions calculated. Now we can start 
      ! to integrate!
      !
      ! Loop through elements in the set and for each element,
      ! loop through the DOF's and cubature points to calculate the
      ! integral:
      
      DO IEL=1,IELmax-IELset+1
        
        ! We make a 'local' approach, i.e. we calculate the values of the
        ! integral into the vector DlocalData and add them later into
        ! the large solution vector.
        
        ! Clear the output vector.
        DlocalData(1:indofTest) = 0.0_DP

        ! Loop over all cubature points on the current element
        DO ICUBP = 1, ncubp
        
          ! calculate the current weighting factor in the cubature formula
          ! in that cubature point.
          !
          ! Take the absolut value of the determinant of the mapping.
          ! In 2D, the determinant is always positive, whereas in 3D,
          ! the determinant might be negative -- that's normal!

          OM = Domega(ICUBP)*ABS(p_Ddetj(ICUBP,IEL))

          ! Loop over the additive factors in the linear form.
          DO IALBET = 1,rform%itermcount
          
            ! Get from Idescriptors the type of the derivatives for the 
            ! test and trial functions. The summand we calculate
            ! here will be:
            !
            ! int_...  f * ( phi_i )_IA
            !
            ! -> IA=0: function value, 
            !      =1: first derivative, 
            !      =2: 2nd derivative,...
            !    as defined in the module 'derivative'.
            
            IA = rform%Idescriptors(IALBET)
            
            ! Multiply OM with the coefficient of the form.
            ! This gives the actual value to multiply the
            ! function value with before summing up to the integral.
            ! Get the precalculated coefficient from the coefficient array.
            AUX = OM * Dcoefficients(IALBET,ICUBP,IEL)
          
            ! Now loop through all possible combinations of DOF's
            ! in the current cubature point. 

            DO IDOFE=1,indofTest
            
              ! Get the value of the basis function 
              ! phi_o in the cubature point. 
              ! Them multiply:
              !    DBAS(..) * AUX
              ! ~= phi_i * coefficient * cub.weight
              ! Summing this up gives the integral, so the contribution
              ! to the vector. 
              !
              ! Simply summing up DBAS(..) * AUX would give
              ! the additive contribution for the vector. We save this
              ! contribution in the local array.

              DlocalData(IDOFE) = DlocalData(IDOFE)+DbasTest(IDOFE,IA,ICUBP,IEL)*AUX
            
            END DO ! IDOFE
            
          END DO ! IALBET

        END DO ! ICUBP 

        ! Incorporate the local vector into the global one.
        ! The 'local' DOF 1..indofTest is mapped to the global DOF using
        ! the IdofsTest array.
        DO IDOFE=1,indofTest
          p_Ddata(IdofsTest(IDOFE,IEL)) = p_Ddata(IdofsTest(IDOFE,IEL)) + DlocalData(IDOFE)
        END DO

      END DO ! IEL

      !CALL ZTIME(DT(10))
    END DO ! IELset
    !$OMP END DO
    
    ! Release memory
    CALL domint_doneIntegration(rintSubset)

    DEALLOCATE(Dcoefficients)
    DEALLOCATE(IdofsTest)
    DEALLOCATE(DbasTest)
    
    !$OMP END PARALLEL

  END DO ! icurrentElementDistr

  !CALL ZTIME(DT(11))
  
  !DO i=2,11
  !  PRINT *,'Time for assembly part ',i,': ',DT(i)-DT(i-1)
  !END DO
  
  
  
  END SUBROUTINE

END MODULE
