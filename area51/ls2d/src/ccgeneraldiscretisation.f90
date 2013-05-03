!##############################################################################
!# ****************************************************************************
!# <name> ccgeneraldiscretisation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains the basic spatial discretisation related routines for
!# CC2D. Here, matrix and RHS creation routines can be found as well as
!# routines to initialise/clean up discretisation structures and routines
!# to read/write solution vectors.
!#
!# The following routines can be found here:
!#
!# 1.) cc_initDiscretisation
!#     -> Initialise the discretisation structure inside of the problem
!#        structure using the parameters from the INI/DAT files.
!#
!# 2.) cc_allocMatVec
!#     -> Allocates memory for vectors/matrices on all levels.
!#
!# 4.) cc_generateStaticMatrices
!#     -> Assembles matrix entries of static matrices (Stokes, B)
!#        on one level
!#
!# 5.) cc_generateBasicMatrices
!#     -> Assembles the matrix entries of all static matrices on all levels.
!#
!# 6.) cc_generateBasicRHS
!#     -> Generates a general RHS vector without any boundary conditions
!#        implemented
!#
!# 7.) cc_doneMatVec
!#     -> Cleanup of matrices/vectors, release all memory
!#
!# 8.) cc_doneDiscretisation
!#     -> Cleanup of the underlying discretisation structures
!#
!# 9.) cc_initInitialSolution
!#     -> Init solution vector according to parameters in the DAT file
!#
!# 10.) cc_writeSolution
!#      -> Write solution vector as configured in the DAT file.
!#
!# Auxiliary routines
!#
!# 1.) cc_getDiscretisation
!#     -> Initialise a block discretisation structure according to
!#        an element type combination.
!#
!# 2.) cc_deriveDiscretisation
!#     -> Derives a block discretisation structure according to
!#        an element type combination from an existing discretisation
!#        structure.
!# </purpose>
!##############################################################################

module ccgeneraldiscretisation

  use fsystem
  use storage
  use linearsolver
  use boundary
  use cubature
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use element
  use spatialdiscretisation
  use bilinearformevaluation
  use linearformevaluation
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use paramlist
  use scalarpde
  use stdoperators
  
  use collection
  use convection
  use vectorio
    
  use ccbasic
  use cccallback
  use ccnonlinearcoreinit
  
  USE levelset
  
  implicit none
  
  !REAL(DP) :: dnu1 = 0.001_DP, dnu2= 0.001_DP, rho1= 2.0_DP,rho2=1.0_DP
  
contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine cc_initDiscretisation (rproblem)
  
!<description>
  ! This routine initialises the discretisation structure of the underlying
  ! problem and saves it to the problem structure.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: I,j,k,ielementType,iElementTypeStabil,icubtemp
  integer(I32) :: icubA,icubB,icubF, icubM
  character(LEN=SYS_NAMELEN) :: sstr
  
    ! An object for saving the domain:
    type(t_boundary), pointer :: p_rboundary
    
    ! An object for saving the triangulation on the domain
    type(t_triangulation), pointer :: p_rtriangulation
    
    ! An object for the block discretisation on one level
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    type(t_spatialDiscretisation), pointer :: p_rdiscretisationMass
    
    
    ! Which discretisation is to use?
    ! Which cubature formula should be used?
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'iElementType',ielementType,3)

    call parlst_getvalue_string (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'scubStokes',sstr,'')
    if (sstr .eq. '') then
      icubtemp = CUB_G2X2
      call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                                'icubStokes',icubtemp,icubtemp)
      icubA = icubtemp
    else
      icubA = cub_igetID(sstr)
    end if

    call parlst_getvalue_string (rproblem%rparamList,'CC-DISCRETISATION',&
                                'scubB',sstr,'')
    if (sstr .eq. '') then
      icubtemp = CUB_G2X2
      call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                                'icubB',icubtemp,icubtemp)
      icubB = icubtemp
    else
      icubB = cub_igetID(sstr)
    end if

    call parlst_getvalue_string (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'scubF',sstr,'')
    if (sstr .eq. '') then
      icubtemp = CUB_G2X2
      call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                                'icubF',icubtemp,icubtemp)
      icubF = icubtemp
    else
      icubF = cub_igetID(sstr)
    end if
    
    ! Stabilisation parameters.
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'IUPWIND',rproblem%rstabilisation%iupwind,0)
    
    call parlst_getvalue_double (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'DUPSAM',rproblem%rstabilisation%dupsam,0.0_DP)

    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'IELEMENTTYPESTABIL',iElementTypeStabil,0)

    ! Now set up discrezisation structures on all levels:

    do i=rproblem%NLMIN,rproblem%NLMAX
    
      ! Ask the problem structure to give us the boundary and triangulation.
      ! We need it for the discretisation.
      p_rboundary => rproblem%rboundary
      p_rtriangulation => rproblem%RlevelInfo(i)%rtriangulation
      p_rdiscretisation => rproblem%RlevelInfo(i)%rdiscretisation
      
      
      ! -----------------------------------------------------------------------
      ! Initialise discretisation structures for the spatial discretisation
      ! -----------------------------------------------------------------------

      ! Initialise the block discretisation according to the element specifier.
      call cc_getDiscretisation (ielementType,p_rdiscretisation,&
          rproblem%RlevelInfo(i)%rtriangulation, rproblem%rboundary, &
          icubA, icubB, icubF)
      
        
    

      
      ! Probably initialise the element type that is to be used for the jump
      ! stabilisation. A value of -1 means: use the same element(s).
      ! In this case, we have to initialise rdiscretisationUEOStabil.
      if (iElementTypeStabil .ne. -1) then
        call cc_deriveDiscretisation (iElementTypeStabil,p_rdiscretisation,&
            rproblem%RlevelInfo(i)%rdiscretisationStabil)
      else
        call spdiscr_deriveBlockDiscr (p_rdiscretisation,&
            rproblem%RlevelInfo(i)%rdiscretisationStabil)
      end if
      
      ! -----------------------------------------------------------------------
      ! Time-dependent problem
      ! -----------------------------------------------------------------------

      call parlst_getvalue_int (rproblem%rparamList, 'TIME-DISCRETISATION', &
                                       'ITIMEDEPENDENCE', j, 0)
      if (j .ne. 0) then
      
        ! Initialise a discretisation structure for the mass matrix.
        ! Copy the discretisation structure of the first (Stokes) block
        ! and replace the cubature-formula identifier by that which is to be
        ! used for the mass matrix.
        call spdiscr_duplicateDiscrSc (p_rdiscretisation%RspatialDiscr(1), &
            rproblem%RlevelInfo(i)%rdiscretisationMass, .true.)
            
        p_rdiscretisationMass => rproblem%RlevelInfo(i)%rdiscretisationMass
        
        call parlst_getvalue_string (rproblem%rparamList,'CC-DISCRETISATION',&
                                    'scubStokes',sstr,'')
        if (sstr .eq. '') then
          icubtemp = CUB_G2X2
          call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                                    'icubStokes',icubtemp,icubtemp)
          icubM = icubtemp
        else
          icubM = cub_igetID(sstr)
        end if

        ! Initialise the cubature formula appropriately.
        do k = 1,p_rdiscretisationMass%inumFESpaces
          p_rdiscretisationMass%RelementDistr(k)%ccubTypeBilForm = icubM
        end do

        ! Should we do mass lumping?
        call parlst_getvalue_int (rproblem%rparamList, 'CC-DISCRETISATION', &
            'IMASS', j, 0)
                                        
        if (j .eq. 0) then
        
          ! How to do lumping?
          call parlst_getvalue_int (rproblem%rparamList, 'CC-DISCRETISATION', &
              'IMASSLUMPTYPE', j, 0)
                                          
          ! Set cubature formula for lumping. The constant from the DAT file corresponds
          ! to one of the LSYSSC_LUMP_xxxx constants for lsyssc_lumpMatrixScalar.
          ! When to do simple mass lumping, replace the cubature formula by one
          ! that is compatible with the corresponding element to generate
          ! a diagonal mass matrix.
          if (j .eq. LSYSSC_LUMP_STD) then
            
            do k = 1,p_rdiscretisationMass%inumFESpaces
              
              j = spdiscr_getLumpCubature (&
                  p_rdiscretisationMass%RelementDistr(k)%celement)
              if (j .ne. 0) then
                icubM = j
              else
                call output_line (&
                    'Unknown cubature formula for mass lumping!', &
                    OU_CLASS_ERROR,OU_MODE_STD,'cc_initDiscretisation')
                call sys_halt()
              end if
              
              ! Set the cubature formula appropriately
              p_rdiscretisationMass%RelementDistr(k)%ccubTypeBilForm = icubM

            end do
          
          end if
        
        end if
       
      end if
      
    end do
                                   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_getDiscretisation (ielementType,rdiscretisation,&
     rtriangulation, rboundary, icubA, icubB, icubF)

!<description>
  ! Initialises a discretisation structure according to an element combination
  ! identifier.
!</description>

!<input>
  ! Element combination identifier to use. The following identifiers are supported:
  ! 0 = Q1~(E031) / Q1~(E031) / Q0
  ! 1 = Q1~(E030) / Q1~(E030) / Q0
  ! 2 = Q1~(EM31) / Q1~(EM31) / Q0
  ! 3 = Q1~(EM30) / Q1~(EM30) / Q0 = standard
  ! 4 = Q2 (E013) / Q2 (E013) / QP1
  ! 5 = Q1~(EM30) / Q1~(EM30) / Q0 unpivoted (much faster than 3 but less stable)
  ! 6 = Q1~(EM30) / Q1~(EM30) / Q0 unscaled (slightly faster than 3 but less stable)
  integer, intent(in) :: ielementType
  
  ! Boundary structure for the discretisation
  type(t_boundary), intent(in), target :: rboundary
  
  ! Triangulation structure for the discretisation
  type(t_triangulation), intent(in), target :: rtriangulation
  
  ! Cubature formula for the velocity matrices
  integer(I32), intent(in) :: icubA
  
  ! Cubature formula for the gradient/divergence matrices
  integer(I32), intent(in) :: icubB
  
  ! Cubature formula for linear forms (RHS vectors)
  integer(I32), intent(in) :: icubF
!</input>

!<output>
  ! Block discretisation structure that defines the discretisation according
  ! to ielementType.
  type(t_blockDiscretisation), intent(out), target :: rdiscretisation
!</output>

!</subroutine>

  ! Number of equations in our problem. velocity+velocity+pressure = 3
  integer, parameter :: nequations = 3
  
  ! local variables
  integer(I32) :: ieltypeUV, ieltypeP
  
    ! Initialise the element type identifiers according to ielementType
    select case (ielementType)
    case (0)
      ieltypeUV = EL_E031
      ieltypeP = EL_Q0

    case (1)
      ieltypeUV = EL_E030
      ieltypeP = EL_Q0

    case (2)
      ieltypeUV = EL_EM31
      ieltypeP = EL_Q0

    case (3)
      ieltypeUV = EL_EM30
      ieltypeP = EL_Q0

    case (4)
      ieltypeUV = EL_Q2
      ieltypeP = EL_QP1
                  
    case (5)
      ieltypeUV = EL_EM30_UNPIVOTED
      ieltypeP = EL_Q0

    case (6)
      ieltypeUV = EL_EM30_UNSCALED
      ieltypeP = EL_Q0

    case (7)
      ieltypeUV = EL_EB50
      ieltypeP = EL_QP1

    case default
      call output_line (&
          'Unknown discretisation: iElementType = '//sys_siL(ielementType,10), &
          OU_CLASS_ERROR,OU_MODE_STD,'cc_initDiscretisation')
      call sys_halt()
    end select
  
    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies 3 blocks in the
    ! solution vector.
    call spdiscr_initBlockDiscr (&
        rdiscretisation,nequations,rtriangulation,rboundary)

    ! rdiscretisation%RspatialDiscr is a list of scalar
    ! discretisation structures for every component of the solution vector.
    ! We have a solution vector with three components:
    !  Component 1 = X-velocity
    !  Component 2 = Y-velocity
    !  Component 3 = Pressure
    ! For simplicity, we set up one discretisation structure for the
    ! velocity...
    call spdiscr_initDiscr_simple ( &
        rdiscretisation%RspatialDiscr(1), &
        ieltypeUV,icubA,rtriangulation, rboundary)
                
    ! Manually set the cubature formula for the RHS as the above routine
    ! uses the same for matrix and vectors.
    rdiscretisation%RspatialDiscr(1)% &
      RelementDistr(1)%ccubTypeLinForm = icubF
                
    ! ...and copy this structure also to the discretisation structure
    ! of the 2nd component (Y-velocity). This needs no additional memory,
    ! as both structures will share the same dynamic information afterwards.
    call spdiscr_duplicateDiscrSc(rdiscretisation%RspatialDiscr(1),&
        rdiscretisation%RspatialDiscr(2))

    ! For the pressure (3rd component), we set up a separate discretisation
    ! structure, as this uses different finite elements for trial and test
    ! functions.
    call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(1),  &
        ieltypeP, icubB,rdiscretisation%RspatialDiscr(3))
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_deriveDiscretisation (ielementType,rsourceDiscretisation,&
     rdestDiscretisation)

!<description>
  ! Derives a 'compatible' discretisation from an existing discretisation
  ! and changes the element pair.
!</description>

!<input>
  ! Element combination identifier to use for the destination discretisation.
  ! The following identifiers are supported:
  ! 0 = Q1~(E031) / Q1~(E031) / Q0
  ! 1 = Q1~(E030) / Q1~(E030) / Q0
  ! 2 = Q1~(EM31) / Q1~(EM31) / Q0
  ! 3 = Q1~(EM30) / Q1~(EM30) / Q0 = standard
  ! 4 = Q2 (E013) / Q2 (E013) / QP1
  ! 5 = Q1~(EM30) / Q1~(EM30) / Q0 unpivoted (much faster than 3 but less stable)
  ! 6 = Q1~(EM30) / Q1~(EM30) / Q0 unscaled (slightly faster than 3 but less stable)
  integer, intent(in) :: ielementType

  ! Block discretisation structure that defines the template discretisation
  ! which element pair is to be changed.
  type(t_blockDiscretisation), intent(in), target :: rsourceDiscretisation
  
!</input>

!<output>
  ! Destination discretisation structure. Receives a copy of the source
  ! discretisation with the element type changed according to ielementType.
  type(t_blockDiscretisation), intent(out), target :: rdestDiscretisation
!</output>

!</subroutine>

  ! Number of equations in our problem. velocity+velocity+pressure = 3
  integer, parameter :: nequations = 3
  
  ! local variables
  integer(I32) :: ieltypeUV, ieltypeP
  
    ! Initialise the element type identifiers according to ielementType
    select case (ielementType)
    case (0)
      ieltypeUV = EL_E031
      ieltypeP = EL_Q0

    case (1)
      ieltypeUV = EL_E030
      ieltypeP = EL_Q0

    case (2)
      ieltypeUV = EL_EM31
      ieltypeP = EL_Q0

    case (3)
      ieltypeUV = EL_EM30
      ieltypeP = EL_Q0

    case (4)
      ieltypeUV = EL_Q2
      ieltypeP = EL_QP1
                  
    case (5)
      ieltypeUV = EL_EM30_UNPIVOTED
      ieltypeP = EL_Q0

    case (6)
      ieltypeUV = EL_EM30_UNSCALED
      ieltypeP = EL_Q0

    case (7)
      ieltypeUV = EL_EB50
      ieltypeP = EL_QP1

    case default
      call output_line (&
          'Unknown discretisation: iElementType = '//sys_siL(ielementType,10), &
          OU_CLASS_ERROR,OU_MODE_STD,'cc_initDiscretisation')
      call sys_halt()
    end select
  
    ! Copy the source discretisation structure.
    call spdiscr_deriveBlockDiscr (rsourceDiscretisation,rdestDiscretisation)
    
    ! Replace the spatial discretisation structures inside by new ones
    ! with the correct element type.
    call spdiscr_deriveSimpleDiscrSc (rsourceDiscretisation%RspatialDiscr(1),  &
        ieltypeUV, SPDISC_CUB_NOCHANGE,rdestDiscretisation%RspatialDiscr(1))
        
    call spdiscr_deriveSimpleDiscrSc (rsourceDiscretisation%RspatialDiscr(2),  &
        ieltypeUV, SPDISC_CUB_NOCHANGE,rdestDiscretisation%RspatialDiscr(2))
        
    call spdiscr_deriveSimpleDiscrSc (rsourceDiscretisation%RspatialDiscr(3),  &
        ieltypeP, SPDISC_CUB_NOCHANGE,rdestDiscretisation%RspatialDiscr(3))
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_doneDiscretisation (rproblem)
  
!<description>
  ! Releases the discretisation from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i,j
   

    do i=rproblem%NLMAX,rproblem%NLMIN,-1
      
      ! Remove the block discretisation structure and all substructures.
      call spdiscr_releaseBlockDiscr(rproblem%RlevelInfo(i)%rdiscretisation)
       call spdiscr_releaseBlockDiscr(rproblem%RlevelInfo(i)%rdiscretisationLS)
      
      ! -----------------------------------------------------------------------
      ! Time-dependent problem
      ! -----------------------------------------------------------------------

      call parlst_getvalue_int (rproblem%rparamList, 'TIME-DISCRETISATION', &
          'ITIMEDEPENDENCE', j, 0)
      if (j .ne. 0) then
        ! Release the mass matrix discretisation.
        call spdiscr_releaseDiscr (rproblem%RlevelInfo(i)%rdiscretisationMass)
      end if

    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_allocMatVec (rproblem,rvector,rrhs)
  
!<description>
  ! Allocates memory for all matrices and vectors of the problem on the heap
  ! by evaluating the parameters in the problem structure.
  ! Matrices/vectors of global importance are added to the collection
  ! structure of the problem, given in rproblem. Matrix/vector entries
  ! are not calculated, they are left uninitialised.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem
  
  ! A vector structure for the solution vector. The structure is initialised,
  ! memory is allocated for the data entries.
  type(t_vectorBlock), intent(INOUT) :: rvector

  ! A vector structure for the RHS vector. The structure is initialised,
  ! memory is allocated for the data entries.
  type(t_vectorBlock), intent(INOUT) :: rrhs
  
 
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i,cmatBuildType,istrongDerivativeBmatrix
  
    ! A pointer to the system matrix and the RHS/solution vectors.
    type(t_matrixScalar), pointer :: p_rmatrixStokes
    type(t_matrixScalar), pointer :: p_rmatrixTemplateFEM,p_rmatrixTemplateGradient
   
    
    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
  
    call parlst_getvalue_int (rproblem%rparamList, 'CC-DISCRETISATION', &
        'ISTRONGDERIVATIVEBMATRIX', istrongDerivativeBmatrix, 0)
  
    ! When the jump stabilisation is used, we have to create an extended
    ! matrix stencil!
    cmatBuildType = BILF_MATC_ELEMENTBASED
    
    if (rproblem%rstabilisation%iupwind .eq. 2) cmatBuildType = BILF_MATC_EDGEBASED
  
    ! Initialise all levels...
    do i = rproblem%NLMIN, rproblem%NLMAX

      ! -----------------------------------------------------------------------
      ! Basic (Navier-) Stokes problem
      ! -----------------------------------------------------------------------

      ! Ask the problem structure to give us the discretisation structure
      p_rdiscretisation => rproblem%RlevelInfo(i)%rdiscretisation
      
      ! The global system looks as follows:
      !
      !    ( A         B1 )
      !    (      A    B2 )
      !    ( B1^T B2^T    )
      !
      ! with A = L + nonlinear Convection. We compute in advance
      ! a standard Stokes matrix L which can be added later to the
      ! convection matrix, resulting in the nonlinear system matrix.
      !
      ! At first, we create 'template' matrices that define the structure
      ! of each of the submatrices in that global system. These matrices
      ! are later used to 'derive' the actual Laplace/Stokes matrices
      ! by sharing the structure (KCOL/KLD).
      !
      ! Get a pointer to the template FEM matrix. This is used for the
      ! Laplace/Stokes matrix and probably for the mass matrix.
      
      
      p_rmatrixTemplateFEM => rproblem%RlevelInfo(i)%rmatrixTemplateFEM

      ! Create the matrix structure
      call bilf_createMatrixStructure (&
                p_rdiscretisation%RspatialDiscr(1),LSYSSC_MATRIX9,&
                p_rmatrixTemplateFEM,cconstrType=cmatBuildType)

      ! In the global system, there are two gradient matrices B1 and B2.
      ! Create a template matrix that defines their structure.
      p_rmatrixTemplateGradient => rproblem%RlevelInfo(i)%rmatrixTemplateGradient
      
      ! Create the matrices structure of the pressure using the 3rd
      ! spatial discretisation structure in p_rdiscretisation%RspatialDiscr.
      call bilf_createMatrixStructure (&
                p_rdiscretisation%RspatialDiscr(3),LSYSSC_MATRIX9,&
                p_rmatrixTemplateGradient,p_rdiscretisation%RspatialDiscr(1))
      
      ! Ok, now we use the matrices from above to create the actual submatrices
      ! that are used in the global system.
      !
      ! Get a pointer to the (scalar) Stokes matrix:
      p_rmatrixStokes => rproblem%RlevelInfo(i)%rmatrixStokes
      
      ! Connect the Stokes matrix to the template FEM matrix such that they
      ! use the same structure.
      !
      ! Don't create a content array yet, it will be created by
      ! the assembly routines later.
      call lsyssc_duplicateMatrix (p_rmatrixTemplateFEM,&
                  p_rmatrixStokes,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
      
      ! Allocate memory for the entries; don't initialise the memory.
      call lsyssc_allocEmptyMatrix (p_rmatrixStokes,LSYSSC_SETM_UNDEFINED)
      
      ! In the global system, there are two coupling matrices B1 and B2.
      ! Both have the structure of the template gradient matrix.
      ! So connect the two B-matrices to the template gradient matrix
      ! such that they share the same structure.
      ! Create the matrices structure of the pressure using the 3rd
      ! spatial discretisation structure in p_rdiscretisation%RspatialDiscr.
      !
      ! Don't create a content array yet, it will be created by
      ! the assembly routines later.
      ! Allocate memory for the entries; don't initialise the memory.
      
      call lsyssc_duplicateMatrix (p_rmatrixTemplateGradient,&
          rproblem%RlevelInfo(i)%rmatrixB1,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                  
      call lsyssc_duplicateMatrix (p_rmatrixTemplateGradient,&
          rproblem%RlevelInfo(i)%rmatrixB2,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                
      call lsyssc_allocEmptyMatrix (rproblem%RlevelInfo(i)%rmatrixB1,&
                                           LSYSSC_SETM_UNDEFINED)
                                           
      call lsyssc_allocEmptyMatrix (rproblem%RlevelInfo(i)%rmatrixB2,&
                                           LSYSSC_SETM_UNDEFINED)

      ! Set up memory for the divergence matrices B1^T and B2^T.
      ! Virtually transpose B1/B2 and probably allocate memory, depending on what
      ! we need. This gives us the divergence matrices D1 and D2.
      call lsyssc_transposeMatrix (rproblem%RlevelInfo(i)%rmatrixB1,&
          rproblem%RlevelInfo(i)%rmatrixD1,LSYSSC_TR_VIRTUAL)

      call lsyssc_transposeMatrix (rproblem%RlevelInfo(i)%rmatrixB2,&
          rproblem%RlevelInfo(i)%rmatrixD2,LSYSSC_TR_VIRTUAL)

      if (istrongDerivativeBmatrix .eq. 1) then

        ! Allocate new memory for the entries of D1 and D2.
        
        call lsyssc_releaseMatrixContent(rproblem%RlevelInfo(i)%rmatrixD1)
        call lsyssc_releaseMatrixContent(rproblem%RlevelInfo(i)%rmatrixD2)

        call lsyssc_allocEmptyMatrix (rproblem%RlevelInfo(i)%rmatrixD1,&
            LSYSSC_SETM_UNDEFINED)
                                            
        call lsyssc_allocEmptyMatrix (rproblem%RlevelInfo(i)%rmatrixD2,&
            LSYSSC_SETM_UNDEFINED)

      end if

      ! -----------------------------------------------------------------------
      ! Temporary vectors
      !
      ! Now on all levels except for the maximum one, create a temporary
      ! vector on that level, based on the block discretisation structure.
      ! It's used for building the matrices on lower levels.
      if (i .lt. rproblem%NLMAX) then
        call lsysbl_createVecBlockByDiscr (&
            rproblem%RlevelInfo(i)%rdiscretisation,&
            rproblem%RlevelInfo(i)%rtempVector,.true.)
      end if

    end do
    
    ! (Only) on the finest level, we need to have to allocate a RHS vector
    ! and a solution vector.
    !
    ! Although we could manually create the solution/RHS vector,
    ! the easiest way to set up the vector structure is
    ! to create it by using our matrix as template.
    ! Initialise the vectors with 0.
    call lsysbl_createVecBlockByDiscr (&
        rproblem%RlevelInfo(rproblem%NLMAX)%rdiscretisation,rrhs,.true.)
       
        
    call lsysbl_createVecBlockByDiscr (&
        rproblem%RlevelInfo(rproblem%NLMAX)%rdiscretisation,rvector,.true.)

  end subroutine

  ! ***************************************************************************
!<subroutine>
  subroutine cc_generateMatVecLS (rproblem,rvector,rrhs)
  
   type(t_problem), intent(INOUT), target :: rproblem
   type(t_vectorblock), INTENT(INOUT) :: rvector,rrhs
  
   type(t_blockDiscretisation), pointer :: p_rdiscretisation
   type(t_matrixScalar), pointer :: p_rmatrixLS
   
   integer :: i, cmatBuildType
   
   !cmatBuildType = BILF_MATC_ELEMENTBASED
   cmatBuildType = BILF_MATC_EDGEBASED
    
    do i=rproblem%NLMIN,rproblem%NLMAX
      p_rdiscretisation => rproblem%RlevelInfo(i)%rdiscretisationLS
     
      call spdiscr_initBlockDiscr (p_rdiscretisation,1,&
          rproblem%RlevelInfo(i)%rtriangulation,rproblem%rboundary)
      call spdiscr_initDiscr_simple ( &
        p_rdiscretisation%RspatialDiscr(1),EL_EM30,CUB_G4X4, &
           rproblem%RlevelInfo(i)%rtriangulation,rproblem%rboundary)
      p_rdiscretisation%RspatialDiscr(1)% &
      RelementDistr(1)%ccubTypeLinForm = CUB_G4X4
    
     
      p_rdiscretisation => rproblem%RlevelInfo(i)%rdiscretisationLS
      p_rmatrixLS => rproblem%RlevelInfo(i)%rmatrixLS
      call bilf_createMatrixStructure (&
                p_rdiscretisation%RspatialDiscr(1),LSYSSC_MATRIX9,&
                p_rmatrixLS,cconstrType=cmatBuildType)
  
      call lsyssc_allocEmptyMatrix (p_rmatrixLS,LSYSSC_SETM_UNDEFINED)
      call stdop_assembleSimpleMatrix (p_rmatrixLS,DER_FUNC,DER_FUNC)
    end do
    
    call lsysbl_createVecBlockByDiscr (p_rdiscretisation,rrhs,.true.)
        
    call lsysbl_createVecBlockByDiscr (p_rdiscretisation,rvector,.true.)
    
    call lsyssc_clearVector (rrhs%rvectorBlock(1))
  
  end subroutine
  ! ***************************************************************************
!<subroutine>

  subroutine cc_generateStaticMatrices (rproblem,rvector,rlevelInfo)
  
!<description>
  ! Calculates entries of all static matrices (Stokes, B-matrices,...)
  ! in the specified problem structure, i.e. the entries of all matrices
  ! that do not change during the computation or which serve as a template for
  ! generating other matrices.
  !
  ! Memory for those matrices must have been allocated before with
  ! allocMatVec!
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT) :: rproblem
  
  type(t_vectorBlock), intent(in), target :: rvector

  ! A level-info structure. The static matrices in this structure are generated.
  type(t_problem_lvl), intent(INOUT),target :: rlevelInfo
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: j,istrongDerivativeBmatrix

    ! A pointer to the Stokes and mass matrix
    type(t_matrixScalar), pointer :: p_rmatrixStokes,p_rmatrixMass

    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    
    ! Structure for the bilinear form for assembling Stokes,...
    TYPE(t_bilinearForm) :: rform
     
    real(DP):: rho1, rho2, dnu1,dnu2
    
    
    call parlst_getvalue_int (rproblem%rparamList, 'CC-DISCRETISATION', &
        'ISTRONGDERIVATIVEBMATRIX', istrongDerivativeBmatrix, 0)
        
    call parlst_getvalue_double(rproblem%rparamList, 'CC-DISCRETISATION', &
        'dnu_first', dnu1, 0.001_DP)
    call parlst_getvalue_double(rproblem%rparamList, 'CC-DISCRETISATION', &
        'dnu_second', dnu2, 0.001_DP)
    
    call parlst_getvalue_double(rproblem%rparamList, 'CC-DISCRETISATION', &
        'rho_first', rho1, 2.0_DP)
    call parlst_getvalue_double(rproblem%rparamList, 'CC-DISCRETISATION', &
        'rho_second', rho2, 1.0_DP)
        

    ! Initialise the collection for the assembly process with callback routines.
    ! Basically, this stores the simulation time in the collection if the
    ! simulation is nonstationary.
    call cc_initCollectForAssembly (rproblem,rproblem%rcollection)
    rproblem%rcollection%p_rvectorQuickAccess1 => rvector
    rproblem%rcollection%DquickAccess(4) = dnu1
    rproblem%rcollection%DquickAccess(5) = dnu2
    
    ! -----------------------------------------------------------------------
    ! Basic (Navier-) Stokes problem
    ! -----------------------------------------------------------------------
    
    ! Ask the problem structure to give us the discretisation structure
    p_rdiscretisation => rlevelInfo%rdiscretisation
    
    ! Get a pointer to the (scalar) Stokes matrix:
    p_rmatrixStokes => rlevelInfo%rmatrixStokes
    
    
    
    ! The global system looks as follows:
    !
    !    ( A         B1 )   ( A         B1 )
    !    (      A    B2 ) = (      A    B2 )
    !    ( B1^T B2^T    )   ( D1   D2      )
    !
    ! with A = L + nonlinear Convection. We compute in advance
    ! a standard Stokes matrix L which can be added later to the
    ! convection matrix, resulting in the nonlinear system matrix,
    ! as well as both B-matrices.
    
!    ! For assembling of the entries, we need a bilinear form,
!    ! which first has to be set up manually.
!    ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
!    ! scalar system matrix in 2D.
!
    rform%itermCount = 2
    rform%Idescriptors(1,1) = DER_DERIV_X
    rform%Idescriptors(2,1) = DER_DERIV_X
    rform%Idescriptors(1,2) = DER_DERIV_Y
    rform%Idescriptors(2,2) = DER_DERIV_Y

    rform%ballCoeffConstant = .FALSE.
    rform%BconstantCoeff = .FALSE.
    rform%Dcoefficients(1)  = rproblem%dnu
    rform%Dcoefficients(2)  = rproblem%dnu
    
   
!
!    ! Now we can build the matrix entries.
!    ! We specify the callback function coeff_Stokes for the coefficients.
!    ! As long as we use constant coefficients, this routine is not used.
!    ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
!    ! the framework will call the callback routine to get analytical data.
!    !
!    ! We pass our collection structure as well to this routine,
!    ! so the callback routine has access to everything what is
!    ! in the collection.
    
   
!    CALL bilf_buildMatrixScalar (rform,.TRUE.,&
!                                 p_rmatrixStokes,coeff_Stokes,&
!                                 rproblem%rcollection)
    
    ! The following call is a replacement for all the lines commented out
    ! above. It directly sets up the Laplace matrix.
    ! If it's necessary to modify the Laplace matrix, remove this command
    ! and comment in the stuff above.
    call stdop_assembleLaplaceMatrix(p_rmatrixStokes,.true.,rproblem%dnu)
    
    ! In the global system, there are two coupling matrices B1 and B2.
    ! These are build either as "int p grad(phi)" (standard case)
    ! or as "int grad(p) phi" (special case, only for nonconstant pressure).
    
    if (istrongDerivativeBmatrix .eq. 0) then

      ! Standard case.
    
      ! Build the first pressure matrix B1.
      call stdop_assembleSimpleMatrix (rlevelInfo%rmatrixB1,&
          DER_FUNC,DER_DERIV_X,-1.0_DP)

      ! Build the second pressure matrix B2.
      call stdop_assembleSimpleMatrix (rlevelInfo%rmatrixB2,&
          DER_FUNC,DER_DERIV_Y,-1.0_DP)
      
      ! Set up the matrices D1 and D2.
      ! For that purpose, virtually transpose B1/B2.
      call lsyssc_transposeMatrix (rlevelInfo%rmatrixB1,&
          rlevelInfo%rmatrixD1,LSYSSC_TR_VIRTUAL)

      call lsyssc_transposeMatrix (rlevelInfo%rmatrixB2,&
          rlevelInfo%rmatrixD2,LSYSSC_TR_VIRTUAL)
          
    else if (istrongDerivativeBmatrix .eq. 1) then
    
      ! Special case for nonconstant pressure.
    
      ! Build the first pressure matrix B1.
      call stdop_assembleSimpleMatrix (rlevelInfo%rmatrixB1,&
          DER_DERIV_X,DER_FUNC,1.0_DP)

      ! Build the second pressure matrix B2.
      call stdop_assembleSimpleMatrix (rlevelInfo%rmatrixB2,&
          DER_DERIV_Y,DER_FUNC,1.0_DP)
      
      ! Virtually transpose D1 and D2 so that we can set them up.
      ! They are saved virtually transposed at the moment, so doing that again
      ! will give us non-transposed matrices.
      call lsyssc_transposeMatrixInSitu (rlevelInfo%rmatrixD1,&
          LSYSSC_TR_VIRTUAL)

      call lsyssc_transposeMatrixInSitu (rlevelInfo%rmatrixD2,&
          LSYSSC_TR_VIRTUAL)
      
      ! Build the first divergence matrix D1^T
      call stdop_assembleSimpleMatrix (rlevelInfo%rmatrixD1,&
          DER_FUNC,DER_DERIV_X,-1.0_DP)

      ! Build the second divergence matrix D2^T
      call stdop_assembleSimpleMatrix (rlevelInfo%rmatrixD2,&
          DER_FUNC,DER_DERIV_Y,-1.0_DP)
          
      ! Ok, finally return their state to 'virtually transposed'.
      call lsyssc_transposeMatrixInSitu (rlevelInfo%rmatrixD1,LSYSSC_TR_VIRTUAL)

      call lsyssc_transposeMatrixInSitu (rlevelInfo%rmatrixD2,LSYSSC_TR_VIRTUAL)
    
    end if
                                
    ! -----------------------------------------------------------------------
    ! Time-dependent problem
    ! -----------------------------------------------------------------------

    call parlst_getvalue_int (rproblem%rparamList, 'TIME-DISCRETISATION', &
                                     'ITIMEDEPENDENCE', j, 0)
    if (j .ne. 0) then
    
      p_rmatrixMass => rlevelInfo%rmatrixMass

      ! If there is an existing mass matrix, release it.
      call lsyssc_releaseMatrix (rlevelInfo%rmatrixMass)

      ! Generate mass matrix. The matrix has basically the same structure as
      ! our template FEM matrix, so we can take that.
      call lsyssc_duplicateMatrix (rlevelInfo%rmatrixTemplateFEM,&
                  p_rmatrixMass,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                  
      ! Change the discretisation structure of the mass matrix to the
      ! correct one; at the moment it points to the discretisation structure
      ! of the Stokes matrix...
      call lsyssc_assignDiscrDirectMat (p_rmatrixMass,rlevelInfo%rdiscretisationMass)

      ! Call the standard matrix setup routine to build the matrix.
     rproblem%rcollection%p_rvectorQuickAccess1 => rvector
     rproblem%rcollection%DquickAccess(4) = rho1
     rproblem%rcollection%DquickAccess(5) = rho2
       
     rform%itermCount = 1
     rform%Idescriptors(1,1) = DER_FUNC
     rform%Idescriptors(2,1) = DER_FUNC

     rform%ballCoeffConstant = .FALSE.
     rform%BconstantCoeff = .FALSE.
     rform%Dcoefficients(1)  = 1.0_DP
    CALL bilf_buildMatrixScalar (rform,.TRUE.,&
                                 p_rmatrixMass,coeff_Mass,&
                                 rproblem%rcollection)
      
                        
!      call stdop_assembleSimpleMatrix (p_rmatrixMass,DER_FUNC,DER_FUNC)
                  
      ! Should we do mass lumping?
      call parlst_getvalue_int (rproblem%rparamList, 'CC-DISCRETISATION', &
          'IMASS', j, 0)
      
                                    
      if (j .eq. 0) then
      
        ! How to do lumping?
        call parlst_getvalue_int (rproblem%rparamList, 'CC-DISCRETISATION', &
            'IMASSLUMPTYPE', j, 0)
                                        
        ! Lump the mass matrix. The constant from the DAT file corresponds
        ! to one of the LSYSSC_LUMP_xxxx constants for lsyssc_lumpMatrixScalar.
        call lsyssc_lumpMatrixScalar (p_rmatrixMass,j)
      
      end if
      
    end if
   
    ! Clean up the collection (as we are done with the assembly, that's it.
    call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_generateBasicMatrices (rproblem,rvector)
  
!<description>
  ! Calculates the entries of all static matrices (Mass, B,...) on all levels.
  !
  ! Memory for those matrices must have been allocated before with
  ! allocMatVec!
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT) :: rproblem
  
  type(t_vectorblock), intent(IN) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_vectorblock) :: rvectorTemp
    integer :: i
    CALL lsysbl_copyvector(rvector,rvectorTemp)
    
    do i=rproblem%NLMAX,rproblem%NLMIN,-1
      call cc_generateStaticMatrices (&
          rproblem,rvectorTemp,rproblem%RlevelInfo(i))
      IF(i.GT.rproblem%NLMIN) &
       ! CALL prjTrace(rvectorTemp,rproblem%RlevelInfo(i-1))
       CALL prjF2C(rvectorTemp,rproblem%rlevelinfo(i-1)%rdiscretisationLS)
    end do
    call lsysbl_releaseVector(rvectortemp)
   
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_generateBasicRHS (rproblem,rrhs,rvector)
  
!<description>
  ! Calculates the entries of the basic right-hand-side vector on the finest
  ! level. Boundary conditions or similar things are not implemented into
  ! the vector.
  ! Memory for the RHS vector must have been allocated in advance.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem
  
  ! The RHS vector which is to be filled with data.
  type(t_vectorBlock), intent(INOUT) :: rrhs
  
  type(t_vectorBlock),intent(in),target :: rvector
!</inputoutput>

!</subroutine>

  ! local variables
  
    ! A bilinear and linear form describing the analytic problem to solve
    type(t_linearForm) :: rlinform
    
    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    
    real(DP) :: rho1, rho2, gravity1, gravity2

    ! Get a pointer to the RHS on the finest level as well as to the
    ! block discretisation structure:
    p_rdiscretisation => rrhs%p_rblockDiscr
    
    call parlst_getvalue_double(rproblem%rparamList, 'CC-DISCRETISATION', &
        'rho_first', rho1, 2.0_DP)
    call parlst_getvalue_double(rproblem%rparamList, 'CC-DISCRETISATION', &
        'rho_second', rho2, 1.0_DP)
    call parlst_getvalue_double(rproblem%rparamList, 'CC-DISCRETISATION', &
        'gravity_first', gravity1, 0.0_DP)
    call parlst_getvalue_double(rproblem%rparamList, 'CC-DISCRETISATION', &
        'gravity_second', gravity2, 0.0_DP)
     
    
    ! The vector structure is already prepared, but the entries are missing.
    !
    ! At first set up the corresponding linear form (f,Phi_j):
    rlinform%itermCount = 1
    rlinform%Idescriptors(1) = DER_FUNC
    
    ! ... and then discretise the RHS to the first subvector of
    ! the block vector using the discretisation structure of the
    ! first block.
    !
    ! We pass our collection structure as well to this routine,
    ! so the callback routine has access to everything what is
    ! in the collection.
    !
    ! Note that the vector is unsorted after calling this routine!
    !
    ! Initialise the collection for the assembly process with callback routines.
    ! Basically, this stores the simulation time in the collection if the
    ! simulation is nonstationary.
    !call cc_initCollectForAssembly (rproblem,rproblem%rcollection)
    rproblem%rcollection%p_rvectorQuickAccess1 => rvector
    rproblem%rcollection%DquickAccess(1) = gravity1
    rproblem%rcollection%DquickAccess(2) = gravity2
    rproblem%rcollection%DquickAccess(4) = rho1
    rproblem%rcollection%DquickAccess(5) = rho2

    ! Discretise the X-velocity part:
    call linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscr(1),rlinform,.true.,&
              rrhs%RvectorBlock(1),coeff_RHS_x,&
              rproblem%rcollection)

    ! And the Y-velocity part:
    call linf_buildVectorScalar (&
              p_rdiscretisation%RspatialDiscr(2),rlinform,.true.,&
              rrhs%RvectorBlock(2),coeff_RHS_y,&
              rproblem%rcollection)
                                
    ! The third subvector must be zero initially - as it represents the RHS of
    ! the equation "div(u) = 0".
    call lsyssc_clearVector(rrhs%RvectorBlock(3))
                                
    ! Clean up the collection (as we are done with the assembly, that's it.
    
    !call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_initInitialSolution (rproblem,rvector)
  
!<description>
  ! Initialises the initial solution vector into rvector. Depending on the settings
  ! in the DAT file this is either zero or read from a file.
!</description>

!<input>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(IN) :: rproblem
!</input>

!<inputoutput>
  ! The solution vector to be initialised. Must be set up according to the
  ! maximum level NLMAX in rproblem!
  type(t_vectorBlock), intent(INOUT) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: istart
    type(t_vectorBlock) :: rvector1,rvector2
    type(t_vectorScalar) :: rvectorTemp
    character(LEN=SYS_STRLEN) :: sarray,sfile,sfileString
    integer :: ilev
    integer :: NEQ
    type(t_interlevelProjectionBlock) :: rprojection

    ! Get the parameter what to do with rvector
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'isolutionStart',istart,0)
    call parlst_getvalue_string (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'ssolutionStart',sfileString,'')

    ! Create a temp vector at level NLMAX-istart+1.
    ilev = rproblem%NLMAX-abs(istart)+1
    
    if (ilev .lt. rproblem%NLMIN) then
      call output_line (&
          'Level of start vector is < NLMIN! Initialising with zero!', &
          OU_CLASS_WARNING,OU_MODE_STD,'cc_initInitialSolution')
      istart = 0
    end if
    
    ! Init with zero? Read from file?
    if (istart .eq. 0) then
      ! Init with zero
      call lsysbl_clearVector (rvector)
    else
      ! Remove possible ''-characters
      read(sfileString,*) sfile

      call lsysbl_createVectorBlock (&
          rproblem%RlevelInfo(ilev)%rdiscretisation,rvector1,.false.)
      
      ! Read in the vector
      call vecio_readBlockVectorHR (&
          rvector1, sarray, .true., 0, sfile, istart .gt. 0)
          
      ! If the vector is on level < NLMAX, we have to bring it to level NLMAX
      do while (ilev .lt. rproblem%NLMAX)
        
        ! Initialise a vector for the higher level and a prolongation structure.
        call lsysbl_createVectorBlock (&
            rproblem%RlevelInfo(ilev+1)%rdiscretisation,rvector2,.false.)
        
        call mlprj_initProjectionVec (rprojection,rvector2)
        
        ! Prolongate to the next higher level.

        NEQ = mlprj_getTempMemoryVec (rprojection,rvector1,rvector2)
        if (NEQ .ne. 0) call lsyssc_createVector (rvectorTemp,NEQ,.false.)
        call mlprj_performProlongation (rprojection,rvector1, &
                                        rvector2,rvectorTemp)
        if (NEQ .ne. 0) call lsyssc_releaseVector (rvectorTemp)
        
        ! Swap rvector1 and rvector2. Release the coarse grid vector.
        call lsysbl_swapVectors (rvector1,rvector2)
        call lsysbl_releaseVector (rvector2)
        
        call mlprj_doneProjection (rprojection)
        
        ! rvector1 is now on level ilev+1
        ilev = ilev+1
        
      end do
      
      ! Copy the resulting vector rvector1 to the output.
      call lsysbl_copyVector (rvector1,rvector)
      
      ! Release the temp vector
      call lsysbl_releaseVector (rvector1)
      
    end if

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_writeSolution (rproblem,rvector)
  
!<description>
  ! Writes a solution vector rvector to a file as configured in the parameters
  ! in the DAT file.
!</description>

!<input>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(IN) :: rproblem

  ! The solution vector to be written out. Must be set up according to the
  ! maximum level NLMAX in rproblem!
  type(t_vectorBlock), intent(IN) :: rvector
!</input>

!</subroutine>

    ! local variables
    integer :: idestLevel
    type(t_vectorBlock) :: rvector1,rvector2
    type(t_vectorScalar) :: rvectorTemp
    character(LEN=SYS_STRLEN) :: sfile,sfileString
    integer :: ilev
    integer :: NEQ
    type(t_interlevelProjectionBlock) :: rprojection
    logical :: bformatted

    ! Get the parameter what to do with rvector
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'iSolutionWrite',idestLevel,0)
    call parlst_getvalue_string (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'sSolutionWrite',sfileString,'')

    if (idestLevel .eq. 0) return ! nothing to do.
    
    ! Remove possible ''-characters
    read(sfileString,*) sfile
    
    bformatted = idestLevel .gt. 0
    idestLevel = rproblem%NLMAX-abs(idestLevel)+1 ! level where to write out

    if (idestLevel .lt. rproblem%NLMIN) then
      call output_line (&
          'Warning: Level for solution vector is < NLMIN! Writing out at level NLMIN!',&
          OU_CLASS_WARNING,OU_MODE_STD,'cc_initInitialSolution')
      idestLevel = rproblem%NLMIN
    end if
    
    ! Interpolate the solution down to level istart.
    call lsysbl_copyVector (rvector,rvector1)   ! creates new rvector1!

    do ilev = rproblem%NLMAX,idestLevel+1,-1
      
      ! Initialise a vector for the lower level and a prolongation structure.
      call lsysbl_createVectorBlock (&
          rproblem%RlevelInfo(ilev-1)%rdiscretisation,rvector2,.false.)
      
      call mlprj_initProjectionVec (rprojection,rvector2)
      
      ! Interpolate to the next higher level.
      ! (Don't 'restrict'! Restriction would be for the dual space = RHS vectors!)

      NEQ = mlprj_getTempMemoryVec (rprojection,rvector2,rvector1)
      if (NEQ .ne. 0) call lsyssc_createVector (rvectorTemp,NEQ,.false.)
      call mlprj_performInterpolation (rprojection,rvector2,rvector1, &
                                       rvectorTemp)
      if (NEQ .ne. 0) call lsyssc_releaseVector (rvectorTemp)
      
      ! Swap rvector1 and rvector2. Release the fine grid vector.
      call lsysbl_swapVectors (rvector1,rvector2)
      call lsysbl_releaseVector (rvector2)
      
      call mlprj_doneProjection (rprojection)
      
    end do

    ! Write out the solution.
    if (bformatted) then
      call vecio_writeBlockVectorHR (rvector1, 'SOLUTION', .true.,&
                                    0, sfile, '(E22.15)')
    else
      call vecio_writeBlockVectorHR (rvector1, 'SOLUTION', .true.,0, sfile)
    end if

    ! Release temp memory.
    call lsysbl_releaseVector (rvector1)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_doneMatVec (rproblem,rvector,rrhs,rvecLS,rrhsLS)
  
!<description>
  ! Releases system matrix and vectors.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(INOUT), target :: rproblem

  ! A vector structure for the solution vector. The structure is cleaned up,
  ! memory is released.
  type(t_vectorBlock), intent(INOUT) :: rvector,rvecLS

  ! A vector structure for the RHS vector. The structure is cleaned up,
  ! memory is released.
  type(t_vectorBlock), intent(INOUT) :: rrhs,rrhsLS
!</inputoutput>

!</subroutine>

    integer :: i

    ! Release matrices and vectors on all levels
    do i=rproblem%NLMAX,rproblem%NLMIN,-1
      ! If there is an existing mass matrix, release it.
      call lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixMass)

      ! Release Stokes, B1 and B2 matrix
      call lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixB2)
      call lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixB1)
      call lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixStokes)
      call lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixLS)
      
      ! Release the template matrices. This is the point, where the
      ! memory of the matrix structure is released.
      call lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixTemplateGradient)
      call lsyssc_releaseMatrix (rproblem%RlevelInfo(i)%rmatrixTemplateFEM)
      
      ! Remove the temp vector that was used for interpolating the solution
      ! from higher to lower levels in the nonlinear iteration.
      if (i .lt. rproblem%NLMAX) then
        call lsysbl_releaseVector(rproblem%RlevelInfo(i)%rtempVector)
      end if
      
    end do

    ! Delete solution/RHS vector
    
    call lsysbl_releaseVector (rvector)
    call lsysbl_releaseVector (rrhs)
    
    
    call lsysbl_releaseVector (rvecLS)
    call lsysbl_releaseVector (rrhsLS)

  end subroutine
  ! ***************************************************************************
! <subroutine>
  SUBROUTINE initSolutionLS(rproblem,rvector)
  
  TYPE(t_problem), INTENT(INOUT) :: rproblem
  TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
  
  TYPE(t_blockDiscretisation), POINTER :: p_rdiscretisation
 
  INTEGER :: i
  
  i = rproblem%NLMAX
  p_rdiscretisation => rproblem%RlevelInfo(i)%rdiscretisationLS
  call lsysbl_createVecBlockByDiscr (p_rdiscretisation,rvector,.true.)
  
   call initLS_analytical(rproblem,rvector)
              
  END SUBROUTINE

  ! ***************************************************************************
end module
