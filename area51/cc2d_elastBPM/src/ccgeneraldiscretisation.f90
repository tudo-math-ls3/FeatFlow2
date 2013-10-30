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
!# 2.) cc_allocMatVec / cc_generateBasicMat / cc_doneMatVec
!#     -> Allocates/ Generates / Releases memory for vectors/matrices on all levels.
!#
!# 4.) cc_allocTemplateMatrices / cc_generateTemplateMatrices /
!#     cc_releaseTemplateMatrices
!#     -> Allocates/Assembles/Releases matrix entries of template matrices
!#        (Stokes, B) on one level
!#
!# 5.) cc_generateBasicMat
!#     -> Assembles the matrix entries of all template matrices on all levels.
!#
!# 6.) cc_generateBasicRHS
!#     -> Generates a general RHS vector without any boundary conditions
!#        implemented
!#
!# 7.) cc_doneDiscretisation
!#     -> Cleanup of the underlying discretisation structures
!#
!# 8.) cc_initInitialSolution
!#     -> Init solution vector according to parameters in the DAT file
!#
!# 9.) cc_writeSolution
!#     -> Write solution vector as configured in the DAT file.
!#
!# 10.) cc_initDynamicLevelInfo
!#      -> Initialises the dynamic level information structure of one level
!#         with basic information
!#
!# 11.) cc_doneDynamicLevelInfo
!#      -> Releases the dynamic level information structure of one level
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
  use multileveloperators
  use analyticprojection
  use feevaluation
  use matrixio !obaid
  
  use collection
  use convection
  use vectorio
  use matrixio
    
  use ccbasic
  use cccallback
  use ccnonlinearcoreinit
  
  implicit none
  
contains
  
  ! ***************************************************************************

!<subroutine>

  subroutine cc_initDynamicLevelInfo (rdynamicLevelInfo)
  
!<description>
  ! Initialises a dynamic level information structure with basic information.
!</description>

!<output>
  ! A dynamic level information structure to be initialised.
  type(t_dynamicLevelInfo), intent(out), target :: rdynamicLevelInfo
!</output>

!</subroutine>

    ! Initialise the BC/FBC structures.
    call bcasm_initDiscreteBC(rdynamicLevelInfo%rdiscreteBC)
    call bcasm_initDiscreteFBC(rdynamicLevelInfo%rdiscreteFBC)
    
    rdynamicLevelInfo%bhasNeumannBoundary = .false.
    rdynamicLevelInfo%hedgesDirichletBC = ST_NOHANDLE
    rdynamicLevelInfo%nedgesDirichletBC = 0

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_doneDynamicLevelInfo (rdynamicLevelInfo)
  
!<description>
  ! Releases a dynamic level information structure.
!</description>

!<inputoutput>
  ! A dynamic level information structure to be initialised.
  type(t_dynamicLevelInfo), intent(inout), target :: rdynamicLevelInfo
!</inputoutput>

!</subroutine>

    ! Release our discrete version of the boundary conditions
    call bcasm_releaseDiscreteBC (rdynamicLevelInfo%rdiscreteBC)
    
    ! as well as the discrete version of the BC`s for fictitious boundaries
    call bcasm_releaseDiscreteFBC (rdynamicLevelInfo%rdiscreteFBC)

    ! Release the Dirichlet edges.
    if (rdynamicLevelInfo%hedgesDirichletBC .ne. ST_NOHANDLE) then
      call storage_free (rdynamicLevelInfo%hedgesDirichletBC)
    end if
    rdynamicLevelInfo%nedgesDirichletBC = 0

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_initDiscretisation (rproblem)
  
!<description>
  ! This routine initialises the discretisation structure of the underlying
  ! problem and saves it to the problem structure.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i,j,k,ielementType,icubA,icubB,icubF, icubM
  integer :: iElementTypeStabil
  character(LEN=SYS_NAMELEN) :: sstr
  
    ! An object for saving the domain:
    type(t_boundary), pointer :: p_rboundary
    
    ! An object for saving the triangulation on the domain
    type(t_triangulation), pointer :: p_rtriangulation
    
    ! An object for the block discretisation on one level
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    type(t_spatialDiscretisation), pointer :: p_rdiscretisationMass
    type(t_spatialDiscretisation), pointer :: p_rdiscretisationMassPressure
    
    ! Which discretisation is to use?
    ! Which cubature formula should be used?
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'iElementType',ielementType,3)

    call parlst_getvalue_string (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'scubStokes',sstr,'')
    if (sstr .eq. '') then
      call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                                'icubStokes',icubA,int(SPDISC_CUB_AUTOMATIC))
    else
      icubA = cub_igetID(sstr)
    end if
!/***/ what about the icubK11 ... icubK22  ?
    call parlst_getvalue_string (rproblem%rparamList,'CC-DISCRETISATION',&
                                'scubB',sstr,'')
    if (sstr .eq. '') then
      call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                                'icubB',icubB,int(SPDISC_CUB_AUTOMATIC))
    else
      icubB = cub_igetID(sstr)
    end if

    call parlst_getvalue_string (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'scubF',sstr,'')
    if (sstr .eq. '') then
      call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                                'icubF',icubF,int(SPDISC_CUB_AUTOMATIC))
    else
      icubF = cub_igetID(sstr)
    end if
    
    ! Stabilisation parameters.
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'IUPWIND',rproblem%rstabilisation%iupwind,0)
    
    call parlst_getvalue_double (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'DUPSAM',rproblem%rstabilisation%dupsam,0.0_DP)

    call parlst_getvalue_double (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'DUPSAMSTAR',rproblem%rstabilisation%dupsamstar,0.0_DP)

    call parlst_getvalue_double (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'DEOJEDGEEXP',rproblem%rstabilisation%deojEdgeExp,2.0_DP)
    
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'ILOCALH',rproblem%rstabilisation%clocalH,1)

    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'IELEMENTTYPESTABIL',iElementTypeStabil,0)

    call parlst_getvalue_string (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'scubEOJ',sstr,'')
    if (sstr .eq. '') then
      call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                                'icubEOJ',i,int(CUB_G4_1D))
      rproblem%rstabilisation%ccubEOJ = i
    else
      rproblem%rstabilisation%ccubEOJ = cub_igetID(sstr)
    end if

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
          int(icubA,I32), int(icubB,I32), int(icubF,I32))
      
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
!/***/     
      ! Add a reference of the velocity discretisation to the rasmTempl
      ! structure. The actual block structure we keep in the levelInfo-structure.
      call spdiscr_duplicateDiscrSc (&
          rproblem%RlevelInfo(i)%rdiscretisationStabil%RspatialDiscr(1), &
          rproblem%RlevelInfo(i)%rasmTempl%rdiscretisationStabil, .true.)
      
      ! -----------------------------------------------------------------------
      ! Mass matrices. They are used in so many cases, it is better we always
      ! have them available.
      ! -----------------------------------------------------------------------

      ! Initialise a discretisation structure for the mass matrices.
      ! Copy the discretisation structure of the first (Stokes) block
      ! and replace the cubature-formula identifier by that which is to be
      ! used for the mass matrix.
      call spdiscr_duplicateDiscrSc (p_rdiscretisation%RspatialDiscr(1), &
          rproblem%RlevelInfo(i)%rasmTempl%rdiscretisationMass,.true.)
      call spdiscr_duplicateDiscrSc (p_rdiscretisation%RspatialDiscr(7), &
          rproblem%RlevelInfo(i)%rasmTempl%rdiscretisationMassPressure,.true.)
          
      p_rdiscretisationMass => rproblem%RlevelInfo(i)%rasmTempl%rdiscretisationMass
      p_rdiscretisationMassPressure => &
          rproblem%RlevelInfo(i)%rasmTempl%rdiscretisationMassPressure
      
      call parlst_getvalue_string (rproblem%rparamList,'CC-DISCRETISATION',&
                                  'scubMass',sstr,'')
      if (sstr .eq. '') then
        call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                                  'icubM',icubM,int(SPDISC_CUB_AUTOMATIC))
      else
        icubM = cub_igetID(sstr)
      end if
      
      if (icubM .ne. SPDISC_CUB_AUTOMATIC) then
        ! Initialise the cubature formula appropriately.
        do k = 1,p_rdiscretisationMass%inumFESpaces
          p_rdiscretisationMass%RelementDistr(k)%ccubTypeBilForm = icubM
          p_rdiscretisationMassPressure%RelementDistr(k)%ccubTypeBilForm = icubM
        end do
      end if

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
            p_rdiscretisationMassPressure%RelementDistr(k)%ccubTypeBilForm = icubM

          end do
        
        end if
      
      end if
      
      ! Initialise the dynamic level information structure with basic information.
      call cc_initDynamicLevelInfo (rproblem%RlevelInfo(i)%rdynamicInfo)
      
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
  !  0 = Q1~(E031) / Q1~(E031) / Q0
  !  1 = Q1~(E030) / Q1~(E030) / Q0
  !  2 = Q1~(EM31) / Q1~(EM31) / Q0
  !  3 = Q1~(EM30) / Q1~(EM30) / Q0 = standard
  !  4 = Q2 (E013) / Q2 (E013) / QP1
  !  5 = Q1~(EM30) / Q1~(EM30) / Q0 unpivoted (much faster than 3 but less stable)
  !  6 = Q1~(EM30) / Q1~(EM30) / Q0 unscaled (slightly faster than 3 but less stable)
  !  7 = Q1~(EM30) / Q1~(EM30) / Q0 new interface implementation
  !  8 = Q1~(EB30) / Q1~(EB30) / Q0
  ! 10 = Q2~(EB50) / Q2~(EB50) / QP1
  ! 11 = Q2~(EM50) / Q2~(EM50) / QP1
  ! 12 = Q2 (E013) / Q2 (E013) / QP1 (nonparametric)
  ! 20 = Q1        / Q1        / Q1
  ! 30 = P1~(E020) / P1~(E020) / P0
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

  ! Number of equations in our problem. uSx + uSy + vSx + vSy + vFx + vFy + p = 7
  integer, parameter :: nequations = 7
  
  ! local variables
  integer(I32) :: ieltypeUV, ieltypeP
  logical :: btriallowed, bquadallowed
  
    btriallowed  = .false.
    bquadallowed = .false.
  
    ! Initialise the element type identifiers according to ielementType
    select case (ielementType)
    case (0)
      ieltypeUV = EL_E031
      ieltypeP = EL_Q0
      bquadallowed = .true.

    case (1)
      ieltypeUV = EL_E030
      ieltypeP = EL_Q0
      bquadallowed = .true.

    case (2)
      ieltypeUV = EL_EM31
      ieltypeP = EL_Q0
      bquadallowed = .true.

    case (3)
      ieltypeUV = EL_EM30
      ieltypeP = EL_Q0
      bquadallowed = .true.

    case (4)
      ieltypeUV = EL_Q2
      ieltypeP = EL_QP1
      bquadallowed = .true.
                  
    case (5)
      ieltypeUV = EL_EM30_UNPIVOTED
      ieltypeP = EL_Q0
      bquadallowed = .true.

    case (6)
      ieltypeUV = EL_EM30_UNSCALED
      ieltypeP = EL_Q0
      bquadallowed = .true.

    case (7)
      ieltypeUV = EL_EM30_NEW
      ieltypeP = EL_Q0
      bquadallowed = .true.

    case (8)
      ieltypeUV = EL_EB30
      ieltypeP = EL_Q0
      bquadallowed = .true.

    case (10)
      ieltypeUV = EL_EB50
      ieltypeP = EL_QP1
      bquadallowed = .true.

    case (11)
      ieltypeUV = EL_EM50
      ieltypeP = EL_QP1
      bquadallowed = .true.

    case (12)
      ieltypeUV = EL_Q2
      ieltypeP = EL_QP1NP
      bquadallowed = .true.

    case (13)
      ieltypeUV = EL_Q2
      ieltypeP = EL_QP1NPD
      bquadallowed = .true.

    case (14)
      ieltypeUV = EL_Q1TB_2D
      ieltypeP = EL_QP1
      bquadallowed = .true.

    case (15)
      ieltypeUV = EL_Q1TBNP_2D
      ieltypeP = EL_QP1NPD
      bquadallowed = .true.

    case (20)
      ieltypeUV = EL_Q1
      ieltypeP = EL_Q1
      bquadallowed = .true.

    case (30)
      ieltypeUV = EL_P1T
      ieltypeP = EL_P0
      btriallowed  = .true.

    case default
      call output_line (&
          'Unknown discretisation: iElementType = '//sys_siL(ielementType,10), &
          OU_CLASS_ERROR,OU_MODE_STD,'cc_getDiscretisation')
      call sys_halt()
    end select
    
    ! Check the mesh to ensure that the discretisation is compatible.
    if ( ((rtriangulation%InelOfType(TRIA_NVETRI2D) .ne. 0) .and. .not. btriallowed) .or. &
         ((rtriangulation%InelOfType(TRIA_NVEQUAD2D) .ne. 0) .and. .not. bquadallowed)) then
      call output_line (&
          "Discretisation does not support the current mesh!", &
          OU_CLASS_ERROR,OU_MODE_STD,'cc_getDiscretisation')
      call sys_halt()
    end if
  
    ! Now we can start to initialise the discretisation. At first, set up
    ! a block discretisation structure that specifies 7 blocks in the
    ! solution vector.
    call spdiscr_initBlockDiscr (&
        rdiscretisation,nequations,rtriangulation,rboundary)

    ! rdiscretisation%RspatialDiscr is a list of scalar
    ! discretisation structures for every component of the solution vector.
    ! We have a solution vector with three components:
    !  Component 1 = X-displacement for the solid skeleton
    !  Component 2 = Y-displacement for the solid skeleton
    !  Component 3 = x-velocity for the solid skeleton
    !  Component 4 = y-velocity for the solid skeleton
    !  Component 5 = x-velocity for the pore fluid
    !  Component 6 = y-velocity for the pore fluid
    !  Component 7 = Pressure
    ! For simplicity, we set up one discretisation structure for the
    ! velocity...
    call spdiscr_initDiscr_simple ( &
        rdiscretisation%RspatialDiscr(1), &
        ieltypeUV,icubA,rtriangulation, rboundary)
                
    ! Manually set the cubature formula for the RHS as the above routine
    ! uses the same for matrix and vectors.
    if (icubF .ne. SPDISC_CUB_AUTOMATIC) then
      rdiscretisation%RspatialDiscr(1)% &
        RelementDistr(1)%ccubTypeLinForm = icubF
    end if
                
    ! ...and copy this structure also to the discretisation structure
    ! of the 2nd component (Y-displacement). This needs no additional memory,
    ! as both structures will share the same dynamic information afterwards.
    call spdiscr_duplicateDiscrSc(rdiscretisation%RspatialDiscr(1),&
        rdiscretisation%RspatialDiscr(2))
    call spdiscr_duplicateDiscrSc(rdiscretisation%RspatialDiscr(1),&
        rdiscretisation%RspatialDiscr(3))
    call spdiscr_duplicateDiscrSc(rdiscretisation%RspatialDiscr(1),&
        rdiscretisation%RspatialDiscr(4))
    call spdiscr_duplicateDiscrSc(rdiscretisation%RspatialDiscr(1),&
        rdiscretisation%RspatialDiscr(5))
    call spdiscr_duplicateDiscrSc(rdiscretisation%RspatialDiscr(1),&
        rdiscretisation%RspatialDiscr(6))

    ! For the pressure (7th component), we set up a separate discretisation
    ! structure, as this uses different finite elements for trial and test
    ! functions.
    call spdiscr_deriveSimpleDiscrSc (rdiscretisation%RspatialDiscr(1),  &
        ieltypeP, icubB,rdiscretisation%RspatialDiscr(7))
  
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
  !  0 = Q1~(E031) / Q1~(E031) / Q0
  !  1 = Q1~(E030) / Q1~(E030) / Q0
  !  2 = Q1~(EM31) / Q1~(EM31) / Q0
  !  3 = Q1~(EM30) / Q1~(EM30) / Q0 = standard
  !  4 = Q2 (E013) / Q2 (E013) / QP1
  !  5 = Q1~(EM30) / Q1~(EM30) / Q0 unpivoted (much faster than 3 but less stable)
  !  6 = Q1~(EM30) / Q1~(EM30) / Q0 unscaled (slightly faster than 3 but less stable)
  !  7 = Q2~(EB50) / Q2~(EB50) / QP1
  !  8 = Q1~(EB30) / Q1~(EB30) / Q0
  ! 10 = Q2~(EB50) / Q2~(EB50) / QP1
  ! 11 = Q2~(EM50) / Q2~(EM50) / QP1
  ! 12 = Q2 (E013) / Q2 (E013) / QP1 (nonparametric)
  ! 20 = Q1        / Q1        / Q1
  ! 30 = P1~(E020) / P1~(E020) / P0
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

  ! Number of equations in our problem: uSx + uSy + vSx + vSy + vFx + vFy + p = 7
  integer, parameter :: nequations = 7
  
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
      ieltypeUV = EL_EM30_NEW
      ieltypeP = EL_Q0

    case (8)
      ieltypeUV = EL_EB30
      ieltypeP = EL_Q0

    case (10)
      ieltypeUV = EL_EB50
      ieltypeP = EL_QP1

    case (11)
      ieltypeUV = EL_EM50
      ieltypeP = EL_QP1

    case (12)
      ieltypeUV = EL_Q2
      ieltypeP = EL_QP1NP

    case (13)
      ieltypeUV = EL_Q2
      ieltypeP = EL_QP1NPD

    case (14)
      ieltypeUV = EL_Q1TB_2D
      ieltypeP = EL_QP1

    case (15)
      ieltypeUV = EL_Q1TBNP_2D
      ieltypeP = EL_QP1NPD

    case (20)
      ieltypeUV = EL_Q1
      ieltypeP = EL_Q1

    case (30)
      ieltypeUV = EL_P1T
      ieltypeP = EL_P0

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

    call spdiscr_deriveSimpleDiscrSc (rsourceDiscretisation%RspatialDiscr(4),  &
        ieltypeUV, SPDISC_CUB_NOCHANGE,rdestDiscretisation%RspatialDiscr(4))
        
    call spdiscr_deriveSimpleDiscrSc (rsourceDiscretisation%RspatialDiscr(5),  &
        ieltypeUV, SPDISC_CUB_NOCHANGE,rdestDiscretisation%RspatialDiscr(5))
        
    call spdiscr_deriveSimpleDiscrSc (rsourceDiscretisation%RspatialDiscr(6),  &
        ieltypeUV, SPDISC_CUB_NOCHANGE,rdestDiscretisation%RspatialDiscr(6))
! pay attention to the last one, it is not like the six before !!!
    call spdiscr_deriveSimpleDiscrSc (rsourceDiscretisation%RspatialDiscr(7),  &
        ieltypeP, SPDISC_CUB_NOCHANGE,rdestDiscretisation%RspatialDiscr(7))
  
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_doneDiscretisation (rproblem)
  
!<description>
  ! Releases the discretisation from the heap.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: i

    do i=rproblem%NLMAX,rproblem%NLMIN,-1
      
      ! Remove the block discretisation structure and all substructures.
      call spdiscr_releaseBlockDiscr(rproblem%RlevelInfo(i)%rdiscretisation)
      call spdiscr_releaseBlockDiscr(rproblem%RlevelInfo(i)%rdiscretisationStabil)
      call spdiscr_releaseDiscr(rproblem%RlevelInfo(i)%rasmTempl%rdiscretisationStabil)
      
      ! Release the mass matrix discretisation.
      call spdiscr_releaseDiscr (rproblem%RlevelInfo(i)%rasmTempl%rdiscretisationMass)
      call spdiscr_releaseDiscr (rproblem%RlevelInfo(i)%rasmTempl%rdiscretisationMassPressure)

      ! Release dynamic level information
      call cc_doneDynamicLevelInfo (rproblem%RlevelInfo(i)%rdynamicInfo)

    end do
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_allocTemplateMatrices (rproblem,rdiscretisation,rasmTempl,&
      rdiscretisationCoarse)
  
!<description>
  ! Allocates memory and generates the structure of all template matrices
  ! in rasmTempl.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout) :: rproblem
  
  ! Discretisation structure that defines how to discretise the different
  ! operators.
  type(t_blockDiscretisation), intent(in), target :: rdiscretisation

  ! A t_asmTemplates structure. The matrices in this structure are generated.
  type(t_asmTemplates), intent(inout), target :: rasmTempl
  
  ! OPTIONAL: Discretisation structure of the level below the level
  ! identified by rdiscretisation. Must be specified on all levels except
  ! for the coarse mesh.
  type(t_blockDiscretisation), intent(in), optional :: rdiscretisationCoarse
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: cmatBuildType,istrongDerivativeBmatrix
  integer :: iprojTypeVelocity, iprojTypePressure
  
    call parlst_getvalue_int (rproblem%rparamList, 'CC-DISCRETISATION', &
        'ISTRONGDERIVATIVEBMATRIX', istrongDerivativeBmatrix, 0)
    
    ! Get the projection types
    call parlst_getvalue_int (rproblem%rparamList, 'CC-PROLREST', &
        'IPROJTYPEVELOCITY', iprojTypeVelocity, 0)
! attention ! iprojTypeVelocity for non-pressure components
! althought the name is pretty misleadig (velocity)
!the name was inherited from the raw cc2d code and it will be
! changed latter
    call parlst_getvalue_int (rproblem%rparamList, 'CC-PROLREST', &
        'IPROJTYPEPRESSURE', iprojTypePressure, 0)
  
    ! When the jump stabilisation is used, we have to create an extended
    ! matrix stencil!
    cmatBuildType = BILF_MATC_ELEMENTBASED
    
    if ((rproblem%rstabilisation%iupwind .eq. CCMASM_STAB_EDGEORIENTED) .or. &
        (rproblem%rstabilisation%iupwind .eq. CCMASM_STAB_EDGEORIENTED2) .or. &
        (rproblem%rstabilisation%iupwind .eq. CCMASM_STAB_FASTEDGEORIENTED)) then
      cmatBuildType = BILF_MATC_EDGEBASED
    end if
  
    ! -----------------------------------------------------------------------
    ! Basic (Navier-) Stokes problem
    ! -----------------------------------------------------------------------

    ! The global system looks as follows:
    !
      !    ( A11   A12   A13   .     A15   .      BS1  )
      !    ( A21   A22   .     A24   .     A26    BS2  )
      !    ( A31   .     A33   .     .     .      .    )
      !    ( .     A42   .     A44   .     .      .    )
      !    ( .     .     A53   .     A55   .      BF1  )
      !    ( .     .     .     A64   .     A66    BF2  )
      !    ( .     .     DS1   DS2   DF1   DF2    .    )
    !
    ! with A55 =  A66 = L + nonlinear Convection. We compute in advance
    ! a standard Stokes matrix L which can be added later to the
    ! convection matrix, resulting in the nonlinear system matrix.
    !
    ! A11 = (2*mu + lambda) * u1_x * v1_x + mu * u1_y * v1_y
    ! A12 = mu * u2_x * v1_y + lambda * u2_y * v1_x
    ! A21 = mu * u1_y * v2_x + lambda * u1_x * v2_y
    ! A22 = mu * u2_x * v2_x + (2*mu + lambda) * u2_y * v2_y
    ! Note: A21 = A12^T
    !
    !/***/ complet the rest !!!
    ! :
    ! :
    ! At first, we create 'template' matrices that define the structure
    ! of each of the submatrices in that global system. These matrices
    ! are later used to 'derive' the actual Laplace/Stokes matrices
    ! by sharing the structure (KCOL/KLD).
    !
    ! Get a pointer to the template FEM matrix. This is used for the
    ! Laplace/Stokes matrix and probably for the mass matrix.
    
    ! Create the matrix structure
    call bilf_createMatrixStructure (&
              rdiscretisation%RspatialDiscr(1),LSYSSC_MATRIX9,&
              rasmTempl%rmatrixTemplateFEM,cconstrType=cmatBuildType)

    ! Create the matrices structure of the pressure using the 3rd
    ! spatial discretisation structure in rdiscretisation%RspatialDiscr.
    call bilf_createMatrixStructure (&
              rdiscretisation%RspatialDiscr(7),LSYSSC_MATRIX9,&
              rasmTempl%rmatrixTemplateFEMPressure)

    ! Create the matrices structure of the pressure using the 3rd
    ! spatial discretisation structure in rdiscretisation%RspatialDiscr.
    call bilf_createMatrixStructure (&
              rdiscretisation%RspatialDiscr(7),LSYSSC_MATRIX9,&
              rasmTempl%rmatrixTemplateGradient,&
              rdiscretisation%RspatialDiscr(1))
              
    ! Transpose the B-structure to get the matrix template for the
    ! divergence matrices.
    call lsyssc_transposeMatrix (rasmTempl%rmatrixTemplateGradient,&
        rasmTempl%rmatrixTemplateDivergence,LSYSSC_TR_STRUCTURE)
    
    ! Ok, now we use the matrices from above to create the actual submatrices
    ! that are used in the global system.
    !
    ! Connect the Stokes matrix to the template FEM matrix such that they
    ! use the same structure.
    !
    ! Do not create a content array yet, it will be created by
    ! the assembly routines later.
    call lsyssc_duplicateMatrix (rasmTempl%rmatrixTemplateFEM,&
                rasmTempl%rmatrixStokes,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)

    call lsyssc_duplicateMatrix (rasmTempl%rmatrixTemplateFEM,&
                rasmTempl%rmatrixK11,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)

    call lsyssc_duplicateMatrix (rasmTempl%rmatrixTemplateFEM,&
                rasmTempl%rmatrixK12,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
!/***/ transpose, K21 = K12^T will be done latter
    call lsyssc_duplicateMatrix (rasmTempl%rmatrixTemplateFEM,&
                rasmTempl%rmatrixK21,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)

    call lsyssc_duplicateMatrix (rasmTempl%rmatrixTemplateFEM,&
                rasmTempl%rmatrixK22,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
    
    ! Allocate memory for the entries; do not initialise the memory.
    call lsyssc_allocEmptyMatrix (rasmTempl%rmatrixStokes,LSYSSC_SETM_UNDEFINED)
    call lsyssc_allocEmptyMatrix (rasmTempl%rmatrixK11,LSYSSC_SETM_UNDEFINED)
    call lsyssc_allocEmptyMatrix (rasmTempl%rmatrixK12,LSYSSC_SETM_UNDEFINED)
!/***/ transpose, K21 = K12^T will be done latter
    call lsyssc_allocEmptyMatrix (rasmTempl%rmatrixK21,LSYSSC_SETM_UNDEFINED)
    call lsyssc_allocEmptyMatrix (rasmTempl%rmatrixK22,LSYSSC_SETM_UNDEFINED)
    
    ! In the global system, there are two coupling matrices B1 and B2.
    ! Both have the structure of the template gradient matrix.
    ! So connect the two B-matrices to the template gradient matrix
    ! such that they share the same structure.
    ! Create the matrices structure of the pressure using the 3rd
    ! spatial discretisation structure in rdiscretisation%RspatialDiscr.
    !
    ! Do not create a content array yet, it will be created by
    ! the assembly routines later.
    ! Allocate memory for the entries; do not initialise the memory.
    
    call lsyssc_duplicateMatrix (rasmTempl%rmatrixTemplateGradient,&
        rasmTempl%rmatrixBS1,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                
    call lsyssc_duplicateMatrix (rasmTempl%rmatrixTemplateGradient,&
        rasmTempl%rmatrixBS2,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)

    call lsyssc_duplicateMatrix (rasmTempl%rmatrixTemplateGradient,&
        rasmTempl%rmatrixBF1,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                
    call lsyssc_duplicateMatrix (rasmTempl%rmatrixTemplateGradient,&
        rasmTempl%rmatrixBF2,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
              
    call lsyssc_allocEmptyMatrix (rasmTempl%rmatrixBS1,LSYSSC_SETM_UNDEFINED)
    call lsyssc_allocEmptyMatrix (rasmTempl%rmatrixBS2,LSYSSC_SETM_UNDEFINED)
    call lsyssc_allocEmptyMatrix (rasmTempl%rmatrixBF1,LSYSSC_SETM_UNDEFINED)
    call lsyssc_allocEmptyMatrix (rasmTempl%rmatrixBF2,LSYSSC_SETM_UNDEFINED)

    ! Set up memory for the divergence matrices D1 and D2.
    call lsyssc_duplicateMatrix (rasmTempl%rmatrixTemplateDivergence,&
        rasmTempl%rmatrixDS1,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                
    call lsyssc_duplicateMatrix (rasmTempl%rmatrixTemplateDivergence,&
        rasmTempl%rmatrixDS2,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)

    call lsyssc_duplicateMatrix (rasmTempl%rmatrixTemplateDivergence,&
        rasmTempl%rmatrixDF1,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                
    call lsyssc_duplicateMatrix (rasmTempl%rmatrixTemplateDivergence,&
        rasmTempl%rmatrixDF2,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)

    call lsyssc_allocEmptyMatrix (rasmTempl%rmatrixDS1,LSYSSC_SETM_UNDEFINED)
    call lsyssc_allocEmptyMatrix (rasmTempl%rmatrixDS2,LSYSSC_SETM_UNDEFINED)
    call lsyssc_allocEmptyMatrix (rasmTempl%rmatrixDF1,LSYSSC_SETM_UNDEFINED)
    call lsyssc_allocEmptyMatrix (rasmTempl%rmatrixDF2,LSYSSC_SETM_UNDEFINED)
    
    ! The D1^T and D2^T matrices are by default the same as B1 and B2.
    ! These matrices may be different for special VANCA variants if
    ! B1 and B2 is different from D1 and D2 (which is actually a rare case).
    call lsyssc_duplicateMatrix (rasmTempl%rmatrixBS1,&
        rasmTempl%rmatrixDS1T,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
                
    call lsyssc_duplicateMatrix (rasmTempl%rmatrixBS2,&
        rasmTempl%rmatrixDS2T,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

    call lsyssc_duplicateMatrix (rasmTempl%rmatrixBF1,&
        rasmTempl%rmatrixDF1T,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
                
    call lsyssc_duplicateMatrix (rasmTempl%rmatrixBF2,&
        rasmTempl%rmatrixDF2T,LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    ! -----------------------------------------------------------------------
    ! Projection matrices
    !
    ! If we use matrix-based projection for prolongation/restriction, then
    ! we need to allocate the memory for the 2-level-matrices here.
    if(present(rdiscretisationCoarse)) then
!/***/ iprojTypeVelocity is for non-pressure components (2 displacements and 4 velocities)
! althought the name is misleading. It was inherited from the raw cc2d code
!it will be changed latter. the same with rmatrixProlVelocity and rmatrixInterpVelocity
! However, one would need different projection types and projection matrices if the FE spaces
! for non pressure components were not same !
      if(iprojTypeVelocity .eq. 1) then
        
        ! Assemble the matrix structure for velocity projection
        call mlop_create2LvlMatrixStruct (&
            rdiscretisationCoarse%RspatialDiscr(1),&
            rdiscretisation%RspatialDiscr(1),&
            LSYSSC_MATRIX9, rasmTempl%rmatrixProlVelocity)
       
        ! And let the interpolation matrix share the structure
        call lsyssc_duplicateMatrix(&
            rasmTempl%rmatrixProlVelocity, &
            rasmTempl%rmatrixInterpVelocity, &
            LSYSSC_DUP_SHARE, LSYSSC_DUP_REMOVE)
        
      end if
      
      if(iprojTypePressure .eq. 1) then
        
        ! Assemble the matrix structure for pressure projection
        call mlop_create2LvlMatrixStruct (&
            rdiscretisationCoarse%RspatialDiscr(7),&
            rdiscretisation%RspatialDiscr(7),&
            LSYSSC_MATRIX9, rasmTempl%rmatrixProlPressure)

        ! And let the interpolation matrix share the structure
        call lsyssc_duplicateMatrix(&
            rasmTempl%rmatrixProlPressure, &
            rasmTempl%rmatrixInterpPressure, &
            LSYSSC_DUP_SHARE, LSYSSC_DUP_REMOVE)
          
      end if

    end if
      
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_generateTemplateMatrices (rproblem,rdiscretisation,rasmTempl)
  
!<description>
  ! Calculates entries of all template matrices (Stokes, B-matrices,...)
  ! in the specified problem structure, i.e. the entries of all matrices
  ! that do not change during the computation or which serve as a template for
  ! generating other matrices.
  !
  ! Memory for those matrices must have been allocated before with
  ! allocMatVec!
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout) :: rproblem
  
  ! Discretisation structure that defines how to discretise the different
  ! operators.
  type(t_blockDiscretisation), intent(in), target :: rdiscretisation

  ! A t_asmTemplates structure. The template matrices in this structure are generated.
  type(t_asmTemplates), intent(inout), target :: rasmTempl
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: j,istrongDerivativeBmatrix

    ! Structure for a precomputed jump stabilisation matrix
    type(t_jumpStabilisation) :: rjumpStabil
    
    ! Structure for the bilinear form for assembling Stokes,...
     TYPE(t_bilinearForm) :: rform

    call parlst_getvalue_int (rproblem%rparamList, 'CC-DISCRETISATION', &
        'ISTRONGDERIVATIVEBMATRIX', istrongDerivativeBmatrix, 0)

    ! Initialise the collection for the assembly process with callback routines.
    ! Basically, this stores the simulation time in the collection if the
    ! simulation is nonstationary.
    call cc_initCollectForAssembly (rproblem,rproblem%rcollection)
    
    ! -----------------------------------------------------------------------
    ! Basic (Navier-) Stokes problem
    ! -----------------------------------------------------------------------
    
    ! The global system looks as follows:
    !
      !    ( A11   A12   A13   .     A15   .      BS1  )
      !    ( A21   A22   .     A24   .     A26    BS2  )
      !    ( A31   .     A33   .     .     .      .    )
      !    ( .     A42   .     A44   .     .      .    )
      !    ( .     .     A53   .     A55   .      BF1  )
      !    ( .     .     .     A64   .     A66    BF2  )
      !    ( .     .     DS1   DS2   DF1   DF2    .    )
    !\***\ Viel zu schreiben here !
    ! with A = L + nonlinear Convection. We compute in advance
    ! a standard Stokes matrix L which can be added later to the
    ! convection matrix, resulting in the nonlinear system matrix,
    ! as well as both B-matrices.
    
!    ! For assembling of the entries, we need a bilinear form,
!    ! which first has to be set up manually.
!    ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
!    ! scalar system matrix in 2D.
!
!    rform%itermCount = 2
!    rform%Idescriptors(1,1) = DER_DERIV_X
!    rform%Idescriptors(2,1) = DER_DERIV_X
!    rform%Idescriptors(1,2) = DER_DERIV_Y
!    rform%Idescriptors(2,2) = DER_DERIV_Y
!
!    ! In the standard case, we have constant coefficients:
!    rform%ballCoeffConstant = .TRUE.
!    rform%BconstantCoeff = .TRUE.
!    rform%Dcoefficients(1)  = rproblem%rphysics%dnu
!    rform%Dcoefficients(2)  = rproblem%rphysics%dnu
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
!                                 rasmTempl%rmatrixStokes,coeff_Stokes,&
!                                 rproblem%rcollection)

    ! The following call is a replacement for all the lines commented out
    ! above. It directly sets up the Laplace matrix.
    ! If it is necessary to modify the Laplace matrix, remove this command
    ! and comment in the stuff above.
    call stdop_assembleLaplaceMatrix (rasmTempl%rmatrixStokes,.true.,rproblem%rphysics%dnu)
!------------------------------- / K11 / -----------------------------------------
   ! For assembling of the entries of K11, we need a bilinear form,
   ! which first has to be set up manually.
   ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
   ! scalar system matrix in 2D.

   rform%itermCount = 2
   rform%Idescriptors(1,1) = DER_DERIV_X
   rform%Idescriptors(2,1) = DER_DERIV_X
   rform%Idescriptors(1,2) = DER_DERIV_Y
   rform%Idescriptors(2,2) = DER_DERIV_Y

   ! In the standard case, we have constant coefficients:
   rform%ballCoeffConstant = .TRUE.
   rform%BconstantCoeff = .TRUE.
   ! (2*mu + lambda) * u1_x * v1_x + mu * u1_y * v1_y
   rform%Dcoefficients(1)  = 2.0_DP*rproblem%rphysics%dmu + rproblem%rphysics%dlambda
   rform%Dcoefficients(2)  = rproblem%rphysics%dmu

   ! Now we can build the matrix entries.
   ! We specify the callback function coeff_Stokes for the coefficients.
   ! As long as we use constant coefficients, this routine is not used.
   ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
   ! the framework will call the callback routine to get analytical data.
   !
   ! We pass our collection structure as well to this routine,
   ! so the callback routine has access to everything what is
   ! in the collection.
   CALL bilf_buildMatrixScalar (rform,.TRUE.,&
                                rasmTempl%rmatrixK11,coeff_K11,&
                                rproblem%rcollection)
!------------------------------- / K22 / -----------------------------------------    
   ! For assembling of the entries of K22, we need a bilinear form,
   ! which first has to be set up manually.
   ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
   ! scalar system matrix in 2D.

   rform%itermCount = 2
   rform%Idescriptors(1,1) = DER_DERIV_X
   rform%Idescriptors(2,1) = DER_DERIV_X
   rform%Idescriptors(1,2) = DER_DERIV_Y
   rform%Idescriptors(2,2) = DER_DERIV_Y

   ! In the standard case, we have constant coefficients:
   rform%ballCoeffConstant = .TRUE.
   rform%BconstantCoeff = .TRUE.
    ! mu * u2_x * v2_x + (2*mu + lambda) * u2_y * v2_y
   rform%Dcoefficients(1)  = rproblem%rphysics%dmu
   rform%Dcoefficients(2)  = 2.0_DP*rproblem%rphysics%dmu + rproblem%rphysics%dlambda

   ! Now we can build the matrix entries.
   ! We specify the callback function coeff_Stokes for the coefficients.
   ! As long as we use constant coefficients, this routine is not used.
   ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
   ! the framework will call the callback routine to get analytical data.
   !
   ! We pass our collection structure as well to this routine,
   ! so the callback routine has access to everything what is
   ! in the collection.
   CALL bilf_buildMatrixScalar (rform,.TRUE.,&
                                rasmTempl%rmatrixK22,coeff_K22,&
                                rproblem%rcollection)
!------------------------------- / K12 / -----------------------------------------    
   ! For assembling of the entries of K12, we need a bilinear form,
   ! which first has to be set up manually.
   ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
   ! scalar system matrix in 2D.

   rform%itermCount = 2
   rform%Idescriptors(1,1) = DER_DERIV_X
   rform%Idescriptors(2,1) = DER_DERIV_Y
   rform%Idescriptors(1,2) = DER_DERIV_Y
   rform%Idescriptors(2,2) = DER_DERIV_X

   ! In the standard case, we have constant coefficients:
   rform%ballCoeffConstant = .TRUE.
   rform%BconstantCoeff = .TRUE.
   ! mu * u2_x * v1_y + lambda * u2_y * v1_x
   rform%Dcoefficients(1)  = rproblem%rphysics%dmu
   rform%Dcoefficients(2)  = rproblem%rphysics%dlambda

   ! Now we can build the matrix entries.
   ! We specify the callback function coeff_Stokes for the coefficients.
   ! As long as we use constant coefficients, this routine is not used.
   ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
   ! the framework will call the callback routine to get analytical data.
   !
   ! We pass our collection structure as well to this routine,
   ! so the callback routine has access to everything what is
   ! in the collection.
   CALL bilf_buildMatrixScalar (rform,.TRUE.,&
                                rasmTempl%rmatrixK12,coeff_K12,&
                                rproblem%rcollection)
!------------------------------- / K21 / ----------------------------------------- 
!/***/ K21 = K12^T ?
   ! For assembling of the entries of K21, we need a bilinear form,
   ! which first has to be set up manually.
   ! We specify the bilinear form (grad Psi_j, grad Phi_i) for the
   ! scalar system matrix in 2D.

   rform%itermCount = 2
   rform%Idescriptors(1,1) = DER_DERIV_Y
   rform%Idescriptors(2,1) = DER_DERIV_X
   rform%Idescriptors(1,2) = DER_DERIV_X
   rform%Idescriptors(2,2) = DER_DERIV_Y

   ! In the standard case, we have constant coefficients:
   rform%ballCoeffConstant = .TRUE.
   rform%BconstantCoeff = .TRUE.
   ! mu * u1_y * v2_x + lambda * u1_x * v2_y
   rform%Dcoefficients(1)  = rproblem%rphysics%dmu
   rform%Dcoefficients(2)  = rproblem%rphysics%dlambda

   ! Now we can build the matrix entries.
   ! We specify the callback function coeff_Stokes for the coefficients.
   ! As long as we use constant coefficients, this routine is not used.
   ! By specifying ballCoeffConstant = BconstantCoeff = .FALSE. above,
   ! the framework will call the callback routine to get analytical data.
   !
   ! We pass our collection structure as well to this routine,
   ! so the callback routine has access to everything what is
   ! in the collection.
   CALL bilf_buildMatrixScalar (rform,.TRUE.,&
                                rasmTempl%rmatrixK21,coeff_K21,&
                                rproblem%rcollection)
!  call matio_writeMatrixHR(rasmTempl%rmatrixK21,'K_21',.true.,0,'K_21','(E20.10)')
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~` 
    ! In the global system, there are two coupling matrices B1 and B2.
    ! These are built either as "int p grad(phi)" (standard case)
    ! or as "int grad(p) phi" (special case, only for nonconstant pressure).
    
    if (istrongDerivativeBmatrix .eq. 0) then

      ! Standard case.
    
      ! Build the first pressure matrices BS1 and BS2
      call stdop_assembleSimpleMatrix (rasmTempl%rmatrixBS1,&
          DER_FUNC,DER_DERIV_X,-rproblem%rphysics%dnSo)

      call stdop_assembleSimpleMatrix (rasmTempl%rmatrixBF1,&
          DER_FUNC,DER_DERIV_X,-rproblem%rphysics%dnFo)

      ! Build the second pressure matrices BS2 and BF2
      call stdop_assembleSimpleMatrix (rasmTempl%rmatrixBS2,&
          DER_FUNC,DER_DERIV_Y,-rproblem%rphysics%dnSo)
      
      call stdop_assembleSimpleMatrix (rasmTempl%rmatrixBF2,&
          DER_FUNC,DER_DERIV_Y,-rproblem%rphysics%dnFo)

      ! Set up the matrices DS1,DS2, DF1 & DF2 by transposing BS1,BS2,BF1 & BF2
      ! For that purpose, virtually transpose BS1 ... BF2.
      call lsyssc_transposeMatrix (rasmTempl%rmatrixBS1,&
          rasmTempl%rmatrixDS1,LSYSSC_TR_CONTENT)

      call lsyssc_transposeMatrix (rasmTempl%rmatrixBF1,&
          rasmTempl%rmatrixDF1,LSYSSC_TR_CONTENT)

      call lsyssc_transposeMatrix (rasmTempl%rmatrixBS2,&
          rasmTempl%rmatrixDS2,LSYSSC_TR_CONTENT)

      call lsyssc_transposeMatrix (rasmTempl%rmatrixBF2,&
          rasmTempl%rmatrixDF2,LSYSSC_TR_CONTENT)
          
    else if (istrongDerivativeBmatrix .eq. 1) then
    
      ! Special case for nonconstant pressure.
    
      ! Build the first pressure matrices BS1 & BF1.
      call stdop_assembleSimpleMatrix (rasmTempl%rmatrixBS1,&
          DER_DERIV_X,DER_FUNC,rproblem%rphysics%dnSo)

      call stdop_assembleSimpleMatrix (rasmTempl%rmatrixBF1,&
          DER_DERIV_X,DER_FUNC,rproblem%rphysics%dnFo)

      ! Build the second pressure matrix BS2 & BF2.
      call stdop_assembleSimpleMatrix (rasmTempl%rmatrixBS2,&
          DER_DERIV_Y,DER_FUNC,rproblem%rphysics%dnSo)

      call stdop_assembleSimpleMatrix (rasmTempl%rmatrixBF2,&
          DER_DERIV_Y,DER_FUNC,rproblem%rphysics%dnFo)

      
      ! Build the first divergence matrices DS1 & DF1
      call stdop_assembleSimpleMatrix (rasmTempl%rmatrixDS1,&
          DER_FUNC,DER_DERIV_X,-rproblem%rphysics%dnSo)

      call stdop_assembleSimpleMatrix (rasmTempl%rmatrixDF1,&
          DER_FUNC,DER_DERIV_X,-rproblem%rphysics%dnFo)

      ! Build the second divergence matrices DS2 & DF2
      call stdop_assembleSimpleMatrix (rasmTempl%rmatrixDS2,&
          DER_FUNC,DER_DERIV_Y,-rproblem%rphysics%dnSo)

      call stdop_assembleSimpleMatrix (rasmTempl%rmatrixDF2,&
          DER_FUNC,DER_DERIV_Y,-rproblem%rphysics%dnFo)
                    
      ! Set up the D1^T and D2^T matrices as transposed D1 and D2.
      ! We have to allocate separate memory for that.
      call lsyssc_allocEmptyMatrix(rasmTempl%rmatrixDS1T,LSYSSC_SETM_UNDEFINED)
      call lsyssc_allocEmptyMatrix(rasmTempl%rmatrixDF1T,LSYSSC_SETM_UNDEFINED)
      call lsyssc_allocEmptyMatrix(rasmTempl%rmatrixDS2T,LSYSSC_SETM_UNDEFINED)
      call lsyssc_allocEmptyMatrix(rasmTempl%rmatrixDF2T,LSYSSC_SETM_UNDEFINED)
      
      ! Transpose only the content, the structure is the same as B1 and B2.
      call lsyssc_transposeMatrix (rasmTempl%rmatrixDS1,&
          rasmTempl%rmatrixDS1T,LSYSSC_TR_CONTENT)

      call lsyssc_transposeMatrix (rasmTempl%rmatrixDF1,&
          rasmTempl%rmatrixDF1T,LSYSSC_TR_CONTENT)

      call lsyssc_transposeMatrix (rasmTempl%rmatrixDS2,&
          rasmTempl%rmatrixDS2T,LSYSSC_TR_CONTENT)

      call lsyssc_transposeMatrix (rasmTempl%rmatrixDF2,&
          rasmTempl%rmatrixDF2T,LSYSSC_TR_CONTENT)
          
    end if
                              
    ! -----------------------------------------------------------------------
    ! Fast edge-oriented stabilisation
    ! -----------------------------------------------------------------------

    if (rproblem%rstabilisation%iupwind .eq. CCMASM_STAB_FASTEDGEORIENTED) then
    
      ! Precompute the jump stabilisation matrix.
      !
      ! Set up the jump stabilisation structure.
      ! There is not much to do, only initialise the viscosity...
      rjumpStabil%dnu = rproblem%rphysics%dnu
      
      ! Set stabilisation parameter
      rjumpStabil%dgamma = rproblem%rstabilisation%dupsam
      rjumpStabil%dgammastar = rproblem%rstabilisation%dupsamstar
      rjumpStabil%deojEdgeExp = rproblem%rstabilisation%deojEdgeExp
      
      ! Matrix weight
      rjumpStabil%dtheta = 1.0_DP
      
      ! Cubature formula
      rjumpStabil%ccubType = rproblem%rstabilisation%ccubEOJ
      
      ! Allocate an empty matrix
      call lsyssc_duplicateMatrix (rasmTempl%rmatrixTemplateFEM,&
                  rasmTempl%rmatrixStabil,LSYSSC_DUP_SHARE,LSYSSC_DUP_EMPTY)
      call lsyssc_clearMatrix (rasmTempl%rmatrixStabil)

      ! Call the jump stabilisation technique to precalculate the
      ! matrix. The operator is assumed to be linear, so we can just
      ! precompute it.
      call conv_jumpStabilisation2d (&
          rjumpStabil, CONV_MODMATRIX, rasmTempl%rmatrixStabil,&
          rdiscretisation=rasmTempl%rdiscretisationStabil)
          
    end if

    ! -----------------------------------------------------------------------
    ! Mass matrices. They are used in so many cases, it is better we always
    ! have them available.
    ! -----------------------------------------------------------------------

    ! If there is an existing mass matrix, release it.
    call lsyssc_releaseMatrix (rasmTempl%rmatrixMass)
    call lsyssc_releaseMatrix (rasmTempl%rmatrixMassPressure)

    ! Generate mass matrix. The matrix has basically the same structure as
    ! our template FEM matrix, so we can take that.
    call lsyssc_duplicateMatrix (rasmTempl%rmatrixTemplateFEM,&
                rasmTempl%rmatrixMass,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
    call lsyssc_duplicateMatrix (rasmTempl%rmatrixTemplateFEMPressure,&
                rasmTempl%rmatrixMassPressure,LSYSSC_DUP_SHARE,LSYSSC_DUP_REMOVE)
                
    ! Change the discretisation structure of the mass matrix to the
    ! correct one; at the moment it points to the discretisation structure
    ! of the Stokes matrix...
    call lsyssc_assignDiscrDirectMat (rasmTempl%rmatrixMass,&
        rasmTempl%rdiscretisationMass)
    call lsyssc_assignDiscrDirectMat (rasmTempl%rmatrixMassPressure,&
        rasmTempl%rdiscretisationMassPressure)

    ! Call the standard matrix setup routine to build the matrix.
    call stdop_assembleSimpleMatrix (rasmTempl%rmatrixMass,DER_FUNC,DER_FUNC)
    call stdop_assembleSimpleMatrix (rasmTempl%rmatrixMassPressure,DER_FUNC,DER_FUNC)
                
    ! Should we do mass lumping?
    call parlst_getvalue_int (rproblem%rparamList, 'CC-DISCRETISATION', &
                              'IMASS', j, 0)
                                    
    if (j .eq. 0) then
    
      ! How to do lumping?
      call parlst_getvalue_int (rproblem%rparamList, 'CC-DISCRETISATION', &
                                'IMASSLUMPTYPE', j, 0)
                                      
      ! Lump the mass matrix. The constant from the DAT file corresponds
      ! to one of the LSYSSC_LUMP_xxxx constants for lsyssc_lumpMatrixScalar.
      call lsyssc_lumpMatrixScalar (rasmTempl%rmatrixMass,j)
      call lsyssc_lumpMatrixScalar (rasmTempl%rmatrixMassPressure,j)
    
    end if
      
    ! Clean up the collection (as we are done with the assembly, that is it.
    call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_generateProjectionMatrices (&
      rdiscretisationCoarse,rdiscretisationFine,rasmTemplFine)
  
!<description>
  ! Calculates the matrix entries for the projection matrices in rasmTemplFine.
  !
  ! Memory for those matrices must have been allocated before with
  ! allocMatVec!
!</description>

!<inputoutput>
  ! Discretisation structure that defines how to discretise the different
  ! operators on the coarse mesh.
  type(t_blockDiscretisation), intent(in), target :: rdiscretisationCoarse

  ! Discretisation structure that defines how to discretise the different
  ! operators on the fine mesh.
  type(t_blockDiscretisation), intent(in), target :: rdiscretisationFine

  ! A t_asmTemplates structure for the fine mesg. The projection
  ! matrices in this structure are generated.
  type(t_asmTemplates), intent(inout), target :: rasmTemplFine
!</inputoutput>

!</subroutine>

    ! Assemble projection matrix for displacement and velocity
! rmatrixProlVelocity is not only for velocity but also for displacement
! however the name is misleading and it will be changed latter. the current
! name was inherited from the original raw cc2d code.
    if(lsyssc_hasMatrixStructure(rasmTemplFine%rmatrixProlVelocity)) then
    
      ! Assemble its entries
      call mlop_build2LvlProlMatrix(&
          rdiscretisationCoarse%RspatialDiscr(1),&
          rdiscretisationFine%RspatialDiscr(1),&
          .true., rasmTemplFine%rmatrixProlVelocity, MLOP_AVRG_MASS)

      call mlop_build2LvlInterpMatrix(&
          rdiscretisationCoarse%RspatialDiscr(1),&
          rdiscretisationFine%RspatialDiscr(1),&
          .true., rasmTemplFine%rmatrixInterpVelocity, MLOP_AVRG_MASS)

    end if

    ! Assemble projection matrix for pressure?
    if(lsyssc_hasMatrixStructure(rasmTemplFine%rmatrixProlPressure)) then
    
      ! Assemble its entries
      call mlop_build2LvlProlMatrix(&
          rdiscretisationCoarse%RspatialDiscr(7),&
          rdiscretisationFine%RspatialDiscr(7),&
          .true., rasmTemplFine%rmatrixProlPressure, MLOP_AVRG_MASS)

      call mlop_build2LvlInterpMatrix(&
          rdiscretisationCoarse%RspatialDiscr(7),&
          rdiscretisationFine%RspatialDiscr(7),&
          .true., rasmTemplFine%rmatrixInterpPressure, MLOP_AVRG_MASS)

    end if

  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine cc_releaseTemplateMatrices (rasmTempl)
  
!<description>
  ! Releases all template matrices from a t_asmTemplates structure.
!</description>

!<inputoutput>
  ! A t_asmTemplates structure to be cleaned up.
  type(t_asmTemplates), intent(inout), target :: rasmTempl
!</inputoutput>

!</subroutine>

    ! If there is an existing mass matrix, release it.
    call lsyssc_releaseMatrix (rasmTempl%rmatrixMass)
    call lsyssc_releaseMatrix (rasmTempl%rmatrixMassPressure)

    ! Release Stokes, B1, B2,... matrices
    call lsyssc_releaseMatrix (rasmTempl%rmatrixDS2T)
    call lsyssc_releaseMatrix (rasmTempl%rmatrixDS1T)
    call lsyssc_releaseMatrix (rasmTempl%rmatrixDS2)
    call lsyssc_releaseMatrix (rasmTempl%rmatrixDS1)
    call lsyssc_releaseMatrix (rasmTempl%rmatrixBS2)
    call lsyssc_releaseMatrix (rasmTempl%rmatrixBS1)
    call lsyssc_releaseMatrix (rasmTempl%rmatrixDF2T)
    call lsyssc_releaseMatrix (rasmTempl%rmatrixDF1T)
    call lsyssc_releaseMatrix (rasmTempl%rmatrixDF2)
    call lsyssc_releaseMatrix (rasmTempl%rmatrixDF1)
    call lsyssc_releaseMatrix (rasmTempl%rmatrixBF2)
    call lsyssc_releaseMatrix (rasmTempl%rmatrixBF1)
!
    call lsyssc_releaseMatrix (rasmTempl%rmatrixK21)
    call lsyssc_releaseMatrix (rasmTempl%rmatrixK12)
    call lsyssc_releaseMatrix (rasmTempl%rmatrixK22)
    call lsyssc_releaseMatrix (rasmTempl%rmatrixK11)
!
    call lsyssc_releaseMatrix (rasmTempl%rmatrixStokes)
    if (lsyssc_hasMatrixStructure(rasmTempl%rmatrixStabil)) then
      call lsyssc_releaseMatrix (rasmTempl%rmatrixStabil)
    end if
    
    ! Release the template matrices. This is the point, where the
    ! memory of the matrix structure is released.
    call lsyssc_releaseMatrix (rasmTempl%rmatrixTemplateDivergence)
    call lsyssc_releaseMatrix (rasmTempl%rmatrixTemplateGradient)
    call lsyssc_releaseMatrix (rasmTempl%rmatrixTemplateFEM)
    call lsyssc_releaseMatrix (rasmTempl%rmatrixTemplateFEMPressure)
    
    ! Release prolongation/interpolation matrices.
    if(lsyssc_hasMatrixStructure(rasmTempl%rmatrixProlVelocity)) then
      call lsyssc_releaseMatrix(rasmTempl%rmatrixInterpVelocity)
      call lsyssc_releaseMatrix(rasmTempl%rmatrixProlVelocity)
    end if
    if(lsyssc_hasMatrixStructure(rasmTempl%rmatrixProlPressure)) then
      call lsyssc_releaseMatrix(rasmTempl%rmatrixInterpPressure)
      call lsyssc_releaseMatrix(rasmTempl%rmatrixProlPressure)
    end if
      
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
  type(t_problem), intent(inout), target :: rproblem
  
  ! A vector structure for the solution vector. The structure is initialised,
  ! memory is allocated for the data entries.
  type(t_vectorBlock), intent(inout) :: rvector

  ! A vector structure for the RHS vector. The structure is initialised,
  ! memory is allocated for the data entries.
  type(t_vectorBlock), intent(inout) :: rrhs
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i
  
    ! Initialise all levels...
    do i=rproblem%NLMIN,rproblem%NLMAX

      ! Generate the matrices on this level.
      if (i .eq. rproblem%NLMIN) then
        call cc_allocTemplateMatrices (rproblem,&
            rproblem%RlevelInfo(i)%rdiscretisation,&
            rproblem%RlevelInfo(i)%rasmTempl)
      else
        call cc_allocTemplateMatrices (rproblem,&
            rproblem%RlevelInfo(i)%rdiscretisation,&
            rproblem%RlevelInfo(i)%rasmTempl,&
            rproblem%RlevelInfo(i-1)%rdiscretisation)
      end if
      
      ! -----------------------------------------------------------------------
      ! Temporary vectors
      !
      ! Now on all levels except for the maximum one, create a temporary
      ! vector on that level, based on the block discretisation structure.
      ! It is used for building the matrices on lower levels.
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

  subroutine cc_generateBasicMat (rproblem)
  
!<description>
  ! Calculates the entries of all template matrices (Mass, B,...) on all levels.
  !
  ! Memory for those matrices must have been allocated before with
  ! allocMatVec!
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout) :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: i

    ! Assemble template matrices
    do i = rproblem%NLMIN, rproblem%NLMAX
      call cc_generateTemplateMatrices (rproblem,&
          rproblem%RlevelInfo(i)%rdiscretisation,&
          rproblem%RlevelInfo(i)%rasmTempl)
    end do
    
    ! Assemble projection matrices
    do i = rproblem%NLMIN+1, rproblem%NLMAX
      call cc_generateProjectionMatrices(&
          rproblem%RlevelInfo(i-1)%rdiscretisation, &
          rproblem%RlevelInfo(i)%rdiscretisation, &
          rproblem%RlevelInfo(i)%rasmTempl)
    end do

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_doneMatVec (rproblem,rvector,rrhs)
  
!<description>
  ! Releases system matrix and vectors.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem

  ! A vector structure for the solution vector. The structure is cleaned up,
  ! memory is released.
  type(t_vectorBlock), intent(inout) :: rvector

  ! A vector structure for the RHS vector. The structure is cleaned up,
  ! memory is released.
  type(t_vectorBlock), intent(inout) :: rrhs
!</inputoutput>

!</subroutine>

    integer :: i

    ! Release matrices and vectors on all levels
    do i=rproblem%NLMAX,rproblem%NLMIN,-1
      ! Release the template matrices.
      call cc_releaseTemplateMatrices (rproblem%RlevelInfo(i)%rasmTempl)
      
      ! Remove the temp vector that was used for interpolating the solution
      ! from higher to lower levels in the nonlinear iteration.
      if (i .lt. rproblem%NLMAX) then
        call lsysbl_releaseVector(rproblem%RlevelInfo(i)%rtempVector)
      end if
      
    end do

    ! Delete solution/RHS vector
    call lsysbl_releaseVector (rvector)
    call lsysbl_releaseVector (rrhs)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine coeff_RHS_const (rdiscretisation,rform, &
                  nelements,npointsPerElement,Dpoints, &
                  IdofsTest,rdomainIntSubset,&
                  Dcoefficients,rcollection)
    
    use basicgeometry
    use triangulation
    use collection
    use scalarpde
    use domainintegration
    
  !<description>
    ! Routine that returns a constant function for additional
    ! assembly of RHS parts in the moving frame formulation.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(in)                   :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(in)                              :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(in)                                         :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(in)                                         :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(in)  :: Dpoints

    ! An array accepting the DOF`s on all elements trial in the trial space.
    ! DIMENSION(\#local DOF`s in test space,nelements)
    integer, dimension(:,:), intent(in) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(in)              :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(inout), optional      :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(out)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
  
    Dcoefficients(:,:,:) = rcollection%DquickAccess(1)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_generateBasicRHS (rproblem,rasmTemplates,rrhsassembly,rrhs)
  
!<description>
  ! Calculates the entries of the basic right-hand-side vector on the finest
  ! level. Boundary conditions or similar things are not implemented into
  ! the vector.
  ! Memory for the RHS vector must have been allocated in advance.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
  
  ! Assembly templates structure of the level of the RHS.
  type(t_asmTemplates), intent(in) :: rasmTemplates
  
  ! RHS assembly structure that defines the RHS.
  type(t_rhsAssembly), intent(inout) :: rrhsAssembly
  
  ! The RHS vector which is to be filled with data.
  type(t_vectorBlock), intent(inout) :: rrhs
!</inputoutput>

!</subroutine>

    ! local variables
    real(DP) :: dreltime,dtimeweight
    real(DP) :: dnSrhoSR, dnFrhoFR
    integer :: iidx1,iidx2
    character(len=SYS_STRLEN) :: sfilename,sarray
  
    ! A linear form describing the analytic problem to solve
    type(t_linearForm) :: rlinform
    
    ! A pointer to the discretisation structure with the data.
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    
    ! Parameters used for the moving frame formulation
    integer :: imovingFrame, iblock
    real(DP), dimension(NDIM2D) :: Dvelocity,Dacceleration
    type(t_collection) :: rlocalcoll
! ################## added by obaid #########################
    ! variable for selecting a specifig boundary region
    type(t_boundaryRegion) :: rboundaryRegion
    call boundary_createRegion(rproblem%rboundary, 1, 3, rboundaryRegion)
! 1: boundary componet (for holed cylinder, we have 2 components)
! 3: boundary segment
!  now go to line (2095 - 2101)
! ###########################################################

    dnSrhoSR = rproblem%rphysics%dnSo*rproblem%rphysics%drhoSR
    dnFrhoFR = rproblem%rphysics%dnFo*rproblem%rphysics%drhoFR

    ! Get a pointer to the RHS on the finest level as well as to the
    ! block discretisation structure:
    p_rdiscretisation => rrhs%p_rblockDiscr
    
    ! Clear the RHS at first.
    ! The third subvector must be zero initially - as it represents the RHS of
    ! the equation "div(u) = 0".
    call lsysbl_clearVector (rrhs)
    
    ! If the RHS is nonzero, generate it.
    if (rrhsAssembly%ctype .ne. 0) then
    
      ! In a first step, check what kind of RHS we have. Probably we have
      ! to read some files to generate it.
      if (rrhsAssembly%ctype .eq. 3) then
        
        ! The basic RHS can be found in the rrhsassembly structure,
        ! it was read in at the beginning of the program.
        !
        ! Multiply with mass matrices to calculate the actual RHS from the nodal vector.
        call lsyssc_scalarMatVec (rasmTemplates%rmatrixMass, &
            rrhsAssembly%rrhsVector%RvectorBlock(1), &
            rrhs%RVectorBlock(1), rrhsAssembly%dmultiplyX, 0.0_DP, .false.)

        call lsyssc_scalarMatVec (rasmTemplates%rmatrixMass, &
            rrhsAssembly%rrhsVector%RvectorBlock(2), &
            rrhs%RVectorBlock(2), rrhsAssembly%dmultiplyY, 0.0_DP, .false.)
!............................................................................
        call lsyssc_scalarMatVec (rasmTemplates%rmatrixMass, &
            rrhsAssembly%rrhsVector%RvectorBlock(3), &
            rrhs%RVectorBlock(3), rrhsAssembly%dmultiplyX, 0.0_DP, .false.)

        call lsyssc_scalarMatVec (rasmTemplates%rmatrixMass, &
            rrhsAssembly%rrhsVector%RvectorBlock(4), &
            rrhs%RVectorBlock(4), rrhsAssembly%dmultiplyY, 0.0_DP, .false.)
!............................................................................
        call lsyssc_scalarMatVec (rasmTemplates%rmatrixMass, &
            rrhsAssembly%rrhsVector%RvectorBlock(5), &
            rrhs%RVectorBlock(5), rrhsAssembly%dmultiplyX, 0.0_DP, .false.)

        call lsyssc_scalarMatVec (rasmTemplates%rmatrixMass, &
            rrhsAssembly%rrhsVector%RvectorBlock(6), &
            rrhs%RVectorBlock(6), rrhsAssembly%dmultiplyY, 0.0_DP, .false.)
!............................................................................       
      else if (rrhsAssembly%ctype .eq. 4) then
        ! Determine the file before and after the current simulation time.
        dreltime = (rproblem%rtimedependence%dtime - rrhsAssembly%dtimeInit) / &
                   (rrhsAssembly%dtimeMax - rrhsAssembly%dtimeInit) *&
                    real(rrhsAssembly%inumfiles,DP)
        
        ! Get the indices. If dtime refers to the max. time, we have to look
        ! at the last interval.
        iidx1 = max(0,min(int(dreltime),rrhsAssembly%inumfiles)) + rrhsAssembly%ifirstindex
        iidx2 = min(iidx1+1,rrhsAssembly%inumfiles) + rrhsAssembly%ifirstindex
            
        ! There is iidx1 <= dreltime <= iidx2 and iidx2=iidx1+0/1.
        ! By subtracting iidx1 from dreltime we obtain a value in the range [0,1]
        ! which is the interpolation weight.
        dtimeweight = dreltime - real(iidx1-rrhsAssembly%ifirstindex,DP)
        
        ! Did we read the files already? Try to avoid reading files at all costs,
        ! it is expensive!
        if (rrhsAssembly%icurrentRhs .ne. iidx1) then
          if (rrhsAssembly%icurrentRhs2 .eq. iidx1) then
            ! Shift the RHS.
            call lsysbl_copyVector (rrhsAssembly%rrhsVector2,rrhsAssembly%rrhsVector)
            rrhsAssembly%icurrentRhs = iidx1
          else
            ! Read the file.
            sfilename = trim(rrhsAssembly%sfilename)//"."//sys_si0(iidx1,5)
            call output_line ("Reading RHS file """//trim(sfilename)//""".",&
                coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
            call vecio_readBlockVectorHR (rrhsAssembly%rrhsVector, sarray, .true.,&
                0, sfilename, rrhsAssembly%iformatted .ne. 0)
            rrhsAssembly%icurrentRhs = iidx1
          end if
        end if

        if (rrhsAssembly%icurrentRhs2 .ne. iidx2) then
          if (rrhsAssembly%icurrentRhs .eq. iidx2) then
            ! Shift the RHS.
            ! this may actually only happen in the case where we repeat a timestep
            ! as we walk backwards in time...
            call lsysbl_copyVector (rrhsAssembly%rrhsVector,rrhsAssembly%rrhsVector2)
            rrhsAssembly%icurrentRhs = iidx2
          else
            ! Read the file.
            sfilename = trim(rrhsAssembly%sfilename)//"."//sys_si0(iidx2,5)
            call output_line ("Reading RHS file """//trim(sfilename)//""".",&
                coutputMode=OU_MODE_STD+OU_MODE_BENCHLOG )
            call vecio_readBlockVectorHR (rrhsAssembly%rrhsVector2, sarray, .true.,&
                0, sfilename, rrhsAssembly%iformatted .ne. 0)
            rrhsAssembly%icurrentRhs2 = iidx2
          end if
        end if
        
        ! Get the RHS by linear interpolation.
        ! Multiply with mass matrices to calculate the actual RHS from the nodal vector.
!..............................................................................       
        call lsyssc_scalarMatVec (rasmTemplates%rmatrixMass, &
            rrhsAssembly%rrhsVector%RvectorBlock(1), &
            rrhs%RVectorBlock(1), rrhsAssembly%dmultiplyX, 0.0_DP, .false.)

        call lsyssc_scalarMatVec (rasmTemplates%rmatrixMass, &
            rrhsAssembly%rrhsVector%RvectorBlock(2), &
            rrhs%RVectorBlock(2), rrhsAssembly%dmultiplyY, 0.0_DP, .false.)
!..............................................................................
        call lsyssc_scalarMatVec (rasmTemplates%rmatrixMass, &
            rrhsAssembly%rrhsVector%RvectorBlock(3), &
            rrhs%RVectorBlock(3), rrhsAssembly%dmultiplyX, 0.0_DP, .false.)

        call lsyssc_scalarMatVec (rasmTemplates%rmatrixMass, &
            rrhsAssembly%rrhsVector%RvectorBlock(4), &
            rrhs%RVectorBlock(4), rrhsAssembly%dmultiplyY, 0.0_DP, .false.)
!..............................................................................
        call lsyssc_scalarMatVec (rasmTemplates%rmatrixMass, &
            rrhsAssembly%rrhsVector%RvectorBlock(5), &
            rrhs%RVectorBlock(5), rrhsAssembly%dmultiplyX, 0.0_DP, .false.)

        call lsyssc_scalarMatVec (rasmTemplates%rmatrixMass, &
            rrhsAssembly%rrhsVector%RvectorBlock(6), &
            rrhs%RVectorBlock(6), rrhsAssembly%dmultiplyY, 0.0_DP, .false.)
!..............................................................................

        if (dtimeweight .gt. 0.0_DP) then
          call lsyssc_scalarMatVec (rasmTemplates%rmatrixMass, &
              rrhsAssembly%rrhsVector2%RvectorBlock(1), &
              rrhs%RVectorBlock(1), rrhsAssembly%dmultiplyX*dtimeweight, &
              (1.0_DP-dtimeweight), .false.)

          call lsyssc_scalarMatVec (rasmTemplates%rmatrixMass, &
              rrhsAssembly%rrhsVector2%RvectorBlock(2), &
              rrhs%RVectorBlock(2), rrhsAssembly%dmultiplyY*dtimeweight, &
              (1.0_DP-dtimeweight), .false.)
!..............................................................................
          call lsyssc_scalarMatVec (rasmTemplates%rmatrixMass, &
              rrhsAssembly%rrhsVector2%RvectorBlock(3), &
              rrhs%RVectorBlock(3), rrhsAssembly%dmultiplyX*dtimeweight, &
              (1.0_DP-dtimeweight), .false.)

          call lsyssc_scalarMatVec (rasmTemplates%rmatrixMass, &
              rrhsAssembly%rrhsVector2%RvectorBlock(4), &
              rrhs%RVectorBlock(4), rrhsAssembly%dmultiplyY*dtimeweight, &
              (1.0_DP-dtimeweight), .false.)
!..............................................................................
          call lsyssc_scalarMatVec (rasmTemplates%rmatrixMass, &
              rrhsAssembly%rrhsVector2%RvectorBlock(5), &
              rrhs%RVectorBlock(5), rrhsAssembly%dmultiplyX*dtimeweight, &
              (1.0_DP-dtimeweight), .false.)

          call lsyssc_scalarMatVec (rasmTemplates%rmatrixMass, &
              rrhsAssembly%rrhsVector2%RvectorBlock(6), &
              rrhs%RVectorBlock(6), rrhsAssembly%dmultiplyY*dtimeweight, &
              (1.0_DP-dtimeweight), .false.)
!..............................................................................
        end if

        !call lsysbl_copyVector (rrhsAssembly%rrhsVector,rrhs)
!        if (dtimeweight .gt. 0.0_DP) then
!
!          call lsysbl_vectorLinearComb (rrhs,rrhsAssembly%rrhsVector2,&
!              (1.0_DP-dtimeweight),dtimeweight)
!        end if
        
      end if
      
      ! The vector structure is already prepared, but the entries are missing.
      !\***\ additional work for body force and tractions to be included latter !!!
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
      call cc_initCollectForAssembly (rproblem,rproblem%rcollection)

      ! Discretise the f_uS1 part of the solid skeleton:
      call linf_buildVectorScalar (&
                p_rdiscretisation%RspatialDiscr(1),rlinform,.false.,&
                rrhs%RvectorBlock(1),coeff_RHS_uSx,&
                rproblem%rcollection)

      ! The f_uS2 part of the solid skeleton:
      call linf_buildVectorScalar (&
                p_rdiscretisation%RspatialDiscr(2),rlinform,.false.,&
                rrhs%RvectorBlock(2),coeff_RHS_uSy,&
                rproblem%rcollection)

      ! Discretise the X-velocity part of the solid skeleton:
      call linf_buildVectorScalar (&
                p_rdiscretisation%RspatialDiscr(3),rlinform,.false.,&
                rrhs%RvectorBlock(3),coeff_RHS_vSx,&
                rproblem%rcollection)

      ! The Y-velocity part of the solid skeleton:
      call linf_buildVectorScalar (&
                p_rdiscretisation%RspatialDiscr(4),rlinform,.false.,&
                rrhs%RvectorBlock(4),coeff_RHS_vSy,&
                rproblem%rcollection)

      ! Discretise the X-velocity part of the pore fluid:
      call linf_buildVectorScalar (&
                p_rdiscretisation%RspatialDiscr(5),rlinform,.false.,&
                rrhs%RvectorBlock(5),coeff_RHS_vFx,&
                rproblem%rcollection)

      ! The Y-velocity part of the pore fluid:
      call linf_buildVectorScalar (&
                p_rdiscretisation%RspatialDiscr(6),rlinform,.false.,&
                rrhs%RvectorBlock(6),coeff_RHS_vFy,&
                rproblem%rcollection)

      ! The pressure/divergence part:
      call linf_buildVectorScalar (&
                p_rdiscretisation%RspatialDiscr(7),rlinform,.false.,&
                rrhs%RvectorBlock(7),coeff_RHS_p,&
                rproblem%rcollection)

! #####################   added by Obaid   ##############################
      rboundaryRegion%dminParam = 2.0_DP
      rboundaryRegion%dmaxParam = 3.0_DP
      iblock = 6  ! f_vf2,       iblock2 for f_us2
      call linf_buildVectorScalarBdr2d(rlinform, CUB_G2_1D, .false., &
          rrhs%RvectorBlock(iblock), RHS_2D_surf, rboundaryRegion, rproblem%rcollection)
! you may need also to specify the magnification factors to scale the deformation
! in line 1623 & 1624 in ccpostprocessing.f90
! #######################################################################
                             
!       ! Is the moving-frame formulatino active?
!       call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
!           'imovingFrame',imovingFrame,0)
!           
!       if (imovingFrame .ne. 0) then
!       
!         ! Get the velocity and acceleration from the callback routine.
!         call getMovingFrameVelocity (Dvelocity,Dacceleration,rproblem%rcollection)
!               
!         ! Assemble a constant RHS with the returned acceleration using the
!         ! above coeff_RHS_const, this realises the moving frame in the
!         ! inner of the domain. We pass the constant function in DquickAccess(1).
!         
!         ! Discretise the X-velocity part:
!         rlocalColl%DquickAccess(1) = Dacceleration(1)
!         call linf_buildVectorScalar (&
!                   p_rdiscretisation%RspatialDiscr(1),rlinform,.false.,&
!                   rrhs%RvectorBlock(1),coeff_RHS_const,&
!                   rlocalcoll)
! 
!         ! And the Y-velocity part:
!         rlocalColl%DquickAccess(1) = Dacceleration(2)
!         call linf_buildVectorScalar (&
!                   p_rdiscretisation%RspatialDiscr(2),rlinform,.false.,&
!                   rrhs%RvectorBlock(2),coeff_RHS_const,&
!                   rlocalcoll)
! 
! ! .....
!         rlocalColl%DquickAccess(1) = Dacceleration(3)
!         call linf_buildVectorScalar (&
!                   p_rdiscretisation%RspatialDiscr(3),rlinform,.false.,&
!                   rrhs%RvectorBlock(3),coeff_RHS_const,&
!                   rlocalcoll)
! 
!         ! And the Y-velocity part:
!         rlocalColl%DquickAccess(1) = Dacceleration(4)
!         call linf_buildVectorScalar (&
!                   p_rdiscretisation%RspatialDiscr(4),rlinform,.false.,&
!                   rrhs%RvectorBlock(4),coeff_RHS_const,&
!                   rlocalcoll)
! ! ....      
!         rlocalColl%DquickAccess(1) = Dacceleration(5)
!         call linf_buildVectorScalar (&
!                   p_rdiscretisation%RspatialDiscr(5),rlinform,.false.,&
!                   rrhs%RvectorBlock(5),coeff_RHS_const,&
!                   rlocalcoll)
! 
!         ! And the Y-velocity part:
!         rlocalColl%DquickAccess(1) = Dacceleration(6)
!         call linf_buildVectorScalar (&
!                   p_rdiscretisation%RspatialDiscr(6),rlinform,.false.,&
!                   rrhs%RvectorBlock(6),coeff_RHS_const,&
!                   rlocalcoll)
! 
!       end if
                                
      ! Clean up the collection (as we are done with the assembly, that is it.
      call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)
      
    end if

  end subroutine

  ! ***************************************************************************

  !<subroutine>

    subroutine fcoeff_solProjection (rdiscretisation, rform, &
                  nelements, npointsPerElement, Dpoints, &
                  IdofsTest, rdomainIntSubset, &
                  Dcoefficients, rcollection)
    
    use fsystem
    use basicgeometry
    use triangulation
    use scalarpde
    use domainintegration
    use spatialdiscretisation
    use collection
    
  !<description>
    ! Called when the initial solution has to be projected into another FEM
    ! space. Evaluated the current initial solution in cubature points.
  !</description>
    
  !<input>
    ! The discretisation structure that defines the basic shape of the
    ! triangulation with references to the underlying triangulation,
    ! analytic boundary boundary description etc.
    type(t_spatialDiscretisation), intent(IN) :: rdiscretisation
    
    ! The linear form which is currently to be evaluated:
    type(t_linearForm), intent(IN) :: rform
    
    ! Number of elements, where the coefficients must be computed.
    integer, intent(IN) :: nelements
    
    ! Number of points per element, where the coefficients must be computed
    integer, intent(IN) :: npointsPerElement
    
    ! This is an array of all points on all the elements where coefficients
    ! are needed.
    ! Remark: This usually coincides with rdomainSubset%p_DcubPtsReal.
    ! DIMENSION(dimension,npointsPerElement,nelements)
    real(DP), dimension(:,:,:), intent(IN) :: Dpoints

    ! An array accepting the DOF`s on all elements test in the test space.
    ! DIMENSION(\#local DOF`s in test space,Number of elements)
    integer, dimension(:,:), intent(IN) :: IdofsTest

    ! This is a t_domainIntSubset structure specifying more detailed information
    ! about the element set that is currently being integrated.
    ! It is usually used in more complex situations (e.g. nonlinear matrices).
    type(t_domainIntSubset), intent(IN) :: rdomainIntSubset

    ! Optional: A collection structure to provide additional
    ! information to the coefficient routine.
    type(t_collection), intent(INOUT), optional :: rcollection
    
  !</input>
  
  !<output>
    ! A list of all coefficients in front of all terms in the linear form -
    ! for all given points on all given elements.
    !   DIMENSION(itermCount,npointsPerElement,nelements)
    ! with itermCount the number of terms in the linear form.
    real(DP), dimension(:,:,:), intent(OUT)                      :: Dcoefficients
  !</output>
    
  !</subroutine>
  
      integer :: icomponent
      
      ! Evaluate copmponent icomponent
      icomponent = rcollection%IquickAccess(1)
      
      call fevl_evaluate_sim (&
          rcollection%p_rvectorQuickAccess1%RvectorBlock(icomponent), &
          rdomainIntSubset, DER_FUNC, Dcoefficients, 1)
  
    end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_initInitialSolution (rproblem,rvector)
  
!<description>
  ! Initialises the initial solution vector into rvector. Depending on the settings
  ! in the DAT file this is either zero or read from a file.
  !
  ! The routine assumes that basic mass matrices have already been constructed.
!</description>

!<input>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout) :: rproblem
!</input>

!<inputoutput>
  ! The solution vector to be initialised. Must be set up according to the
  ! maximum level NLMAX in rproblem!
  type(t_vectorBlock), intent(inout) :: rvector
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iinitialSolutionLevel,ctypeInitialSolution,ielementTypeInitialSolution
    type(t_vectorBlock), target :: rvector1,rvector2
    type(t_vectorScalar) :: rvectorTemp
    character(LEN=SYS_STRLEN) :: sarray,sfile,sfileString
    integer :: ilev,ierror
    integer :: NEQ
    type(t_interlevelProjectionBlock) :: rprojection
    type(t_linearForm) :: rlinform
    type(t_blockDiscretisation), pointer :: p_rdiscretisation
    type(t_blockDiscretisation) :: rdiscretisationInitSol
    type(t_linsolNode), pointer :: p_rsolverNode,p_rpreconditioner
    type(t_matrixBlock), dimension(1) :: Rmatrices
    type(t_vectorBlock) :: rsingleRHS,rsingleSol
    type(t_collection) :: rcollection

    ! Get the parameter what to do with rvector
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
        'ctypeInitialSolution',ctypeInitialSolution,0)
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
        'iinitialSolutionLevel',iinitialSolutionLevel,0)
    call parlst_getvalue_string (rproblem%rparamList,'CC-DISCRETISATION',&
        'sinitialSolutionFilename',sfileString,'')

    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
        'ielementTypeInitialSolution',ielementTypeInitialSolution,-1)

    ! What is the type of the initial solution?
    select case (ctypeInitialSolution)
    case (0)
      ! Init with zero
      call lsysbl_clearVector (rvector)
      
    case (1,2) !/***/ in our case we do not consider this case
      ! We have to read a file -- formatted or unformatted.
      !
      ! Create a temp vector at level NLMAX-iinitialSolutionLevel+1.
      if (iinitialSolutionLevel .gt. 0) then
        ilev = iinitialSolutionLevel
      else
        ilev = rproblem%NLMAX-abs(iinitialSolutionLevel)
      end if
      
      if (ilev .lt. rproblem%NLMIN) then
        call output_line (&
            'Level of start vector is < NLMIN! Initialising with zero!', &
            OU_CLASS_WARNING,OU_MODE_STD,'cc_initInitialSolution')
        ilev = rproblem%NLMIN
      end if

      if (ilev .gt. rproblem%NLMAX) then
        call output_line (&
            'Level of start vector is > NLMAX! Initialising with zero!', &
            OU_CLASS_WARNING,OU_MODE_STD,'cc_initInitialSolution')
        ilev = rproblem%NLMAX
      end if
      
      ! Remove possible ''-characters
      read(sfileString,*) sfile

      ! Create a basic block vector that takes our solution.
      call lsysbl_createVectorBlock (&
          rproblem%RlevelInfo(ilev)%rdiscretisation,rvector1,.false.)
      
      ! Can we directly read in the solution?
      if (ielementTypeInitialSolution .eq. -1) then
        ! Read in the vector to the 2nd temp vector.
        ! This will automatically create rvector2.
        call vecio_readBlockVectorHR (&
          rvector2, sarray, .true., 0, sfile, ctypeInitialSolution .eq. 1)
          
        ! Copy the first three components; maybe rvector2 has more...
        ! Copy only the data, ignore the discretisation.
        if ((rvector2%RvectorBlock(1)%NEQ .ne. rvector1%RvectorBlock(1)%NEQ) .or.&
            (rvector2%RvectorBlock(2)%NEQ .ne. rvector1%RvectorBlock(2)%NEQ) .or.&
            (rvector2%RvectorBlock(3)%NEQ .ne. rvector1%RvectorBlock(3)%NEQ) .or.&
            (rvector2%RvectorBlock(4)%NEQ .ne. rvector1%RvectorBlock(4)%NEQ) .or.&
            (rvector2%RvectorBlock(5)%NEQ .ne. rvector1%RvectorBlock(5)%NEQ) .or.&
            (rvector2%RvectorBlock(6)%NEQ .ne. rvector1%RvectorBlock(6)%NEQ) .or.&
            (rvector2%RvectorBlock(7)%NEQ .ne. rvector1%RvectorBlock(7)%NEQ)) then
          call output_line (&
              'Start vector has invalid size!', &
              OU_CLASS_WARNING,OU_MODE_STD,'cc_initInitialSolution')
          call sys_halt()
        end if
        
        call lsyssc_duplicateVector (rvector2%RvectorBlock(1),rvector1%RvectorBlock(1),&
            LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPYOVERWRITE)
        call lsyssc_duplicateVector (rvector2%RvectorBlock(2),rvector1%RvectorBlock(2),&
            LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPYOVERWRITE)
        call lsyssc_duplicateVector (rvector2%RvectorBlock(3),rvector1%RvectorBlock(3),&
            LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPYOVERWRITE)
        call lsyssc_duplicateVector (rvector2%RvectorBlock(4),rvector1%RvectorBlock(4),&
            LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPYOVERWRITE)
        call lsyssc_duplicateVector (rvector2%RvectorBlock(5),rvector1%RvectorBlock(5),&
            LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPYOVERWRITE)
        call lsyssc_duplicateVector (rvector2%RvectorBlock(6),rvector1%RvectorBlock(6),&
            LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPYOVERWRITE)
        call lsyssc_duplicateVector (rvector2%RvectorBlock(7),rvector1%RvectorBlock(7),&
            LSYSSC_DUP_IGNORE,LSYSSC_DUP_COPYOVERWRITE)
            
        ! Release the temp vector
        call lsysbl_releaseVector (rvector2)
      else
        ! Create a compatible discretisation for the alternative element type
        call cc_deriveDiscretisation (ielementTypeInitialSolution,&
            rproblem%RlevelInfo(ilev)%rdiscretisation,&
            rdiscretisationInitSol)
       
        ! Create a temp vector and read the solution
        call lsysbl_createVectorBlock (&
            rdiscretisationInitSol,rvector2,.false.)
        call vecio_readBlockVectorHR (&
          rvector2, sarray, .true., 0, sfile, ctypeInitialSolution .eq. 1)
          
        ! Put the vector to the collection, will be needed for the projection.
        rcollection%p_rvectorQuickAccess1 => rvector2
          
        ! Project the solution to out current FEM space, so to say to the
        ! vector vector1.
        !
        ! X-velocity !/***/ for 3 ...6 rmatrixMass 
        ! are not multiplied by factors dnSrhoSR, dnFrhoFR here and it should not
        rcollection%IquickAccess(1) = 1
        call anprj_analytL2projectionByMass (rvector1%RvectorBlock(1),&
            rproblem%RlevelInfo(ilev)%rasmTempl%rmatrixMass,&
            fcoeff_solProjection,rcollection)

        ! Y-velocity
        rcollection%IquickAccess(1) = 2
        call anprj_analytL2projectionByMass (rvector1%RvectorBlock(2),&
            rproblem%RlevelInfo(ilev)%rasmTempl%rmatrixMass,&
            fcoeff_solProjection,rcollection)
!

        rcollection%IquickAccess(1) = 3
        call anprj_analytL2projectionByMass (rvector1%RvectorBlock(3),&
            rproblem%RlevelInfo(ilev)%rasmTempl%rmatrixMass,&
            fcoeff_solProjection,rcollection)

        ! Y-velocity
        rcollection%IquickAccess(1) = 4
        call anprj_analytL2projectionByMass (rvector1%RvectorBlock(4),&
            rproblem%RlevelInfo(ilev)%rasmTempl%rmatrixMass,&
            fcoeff_solProjection,rcollection)

!

        rcollection%IquickAccess(1) = 5
        call anprj_analytL2projectionByMass (rvector1%RvectorBlock(5),&
            rproblem%RlevelInfo(ilev)%rasmTempl%rmatrixMass,&
            fcoeff_solProjection,rcollection)

        ! Y-velocity
        rcollection%IquickAccess(1) = 6
        call anprj_analytL2projectionByMass (rvector1%RvectorBlock(6),&
            rproblem%RlevelInfo(ilev)%rasmTempl%rmatrixMass,&
            fcoeff_solProjection,rcollection)
            
        ! Pressure
        rcollection%IquickAccess(1) = 7
        call anprj_analytL2projectionByMass (rvector1%RvectorBlock(7),&
            rproblem%RlevelInfo(ilev)%rasmTempl%rmatrixMassPressure,&
            fcoeff_solProjection,rcollection)
        
        ! rvector2 is not needed anymore, as well as the temporary discretisation.
        call lsysbl_releaseVector (rvector2)
        call spdiscr_releaseBlockDiscr (rdiscretisationInitSol)
      end if

          
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

    case (3)
    
      ! We have to create the solution by analytical callback functions.
      ! To do this for an arbitrary finite element, we have to do an
      ! L2 projection for X,Y and P. That means we have to solve:
      !     <u,phi> = <u_analytic,phi>
      ! which means in the FE context:
      !     Mu = f
      ! with f a RHS created by analytical point values.

      if (rproblem%MSHOW_Initialisation .ge. 1) &
        call output_line('Preparing L2 projection of analytical initial solution...')

      ! Get a pointer to the RHS on the finest level as well as to the
      ! block discretisation structure:
      p_rdiscretisation => rvector%p_rblockDiscr
      
      ! The vector structure is already prepared, but the entries are missing.
      !
      ! At first set up the corresponding linear form (u_analytical,Phi_j):
      rlinform%itermCount = 1
      rlinform%Idescriptors(1) = DER_FUNC
      
      ! Create a temp vector as RHS.
      call lsysbl_createVectorBlock (&
          rproblem%RlevelInfo(rproblem%NLMAX)%rdiscretisation,rvector1,.false.)
      
      ! Assemble the RHS.
      !
      ! Initialise the collection for the assembly process with callback routines.
      ! Basically, this stores the simulation time in the collection if the
      ! simulation is nonstationary.
      call cc_initCollectForAssembly (rproblem,rproblem%rcollection)

      ! Discretise the X- and y- solid skeleton displacement parts:
      call linf_buildVectorScalar (&
                p_rdiscretisation%RspatialDiscr(1),rlinform,.true.,&
                rvector1%RvectorBlock(1),coeff_AnalyticSolution_uSx,&
                rproblem%rcollection)

      call linf_buildVectorScalar (&
                p_rdiscretisation%RspatialDiscr(2),rlinform,.true.,&
                rvector1%RvectorBlock(2),coeff_AnalyticSolution_uSy,&
                rproblem%rcollection)

      ! Discretise the X- and y- solid skeleton velocity part:
      call linf_buildVectorScalar (&
                p_rdiscretisation%RspatialDiscr(3),rlinform,.true.,&
                rvector1%RvectorBlock(3),coeff_AnalyticSolution_vSx,&
                rproblem%rcollection)

      ! Y-velocity:
      call linf_buildVectorScalar (&
                p_rdiscretisation%RspatialDiscr(4),rlinform,.true.,&
                rvector1%RvectorBlock(4),coeff_AnalyticSolution_vSy,&
                rproblem%rcollection)

      ! Discretise the X- and y- pore fluid velocity part:
      call linf_buildVectorScalar (&
                p_rdiscretisation%RspatialDiscr(5),rlinform,.true.,&
                rvector1%RvectorBlock(5),coeff_AnalyticSolution_vFx,&
                rproblem%rcollection)

      ! Y-velocity:
      call linf_buildVectorScalar (&
                p_rdiscretisation%RspatialDiscr(6),rlinform,.true.,&
                rvector1%RvectorBlock(6),coeff_AnalyticSolution_vFy,&
                rproblem%rcollection)

      ! Pressure:
      call linf_buildVectorScalar (&
                p_rdiscretisation%RspatialDiscr(7),rlinform,.true.,&
                rvector1%RvectorBlock(7),coeff_AnalyticSolution_P,&
                rproblem%rcollection)
                                  
      ! Clean up the collection (as we are done with the assembly.
      call cc_doneCollectForAssembly (rproblem,rproblem%rcollection)
      
      ! Now prepare a linear solver for solving with the mass matrix.
      ! This is just a simple 1-level Jacobi solver as the mass matrix
      ! is usually well conditioned.
      call linsol_initJacobi (p_rpreconditioner)
      call linsol_initDefCorr (p_rsolverNode, p_rpreconditioner)
      
      p_rpreconditioner%domega = 0.8_DP
      p_rsolverNode%depsRel = SYS_EPSREAL_DP * 100.0_DP
      p_rsolverNode%nmaxiterations = 1000
      p_rsolverNode%ioutputLevel = 1

      if (rproblem%MSHOW_Initialisation .ge. 1) &
        call output_line('Solving L2 projection for U1...')
      
      ! We solve separately for the 7 components.
      ! Prepare the solver for the velocity.
      ! Prepare the solver for the X-velocity.
      call lsysbl_createMatFromScalar(&
          rproblem%RlevelInfo(rproblem%NLMAX)%rasmTempl%rmatrixMass,Rmatrices(1))
      call linsol_setMatrices (p_rsolverNode,Rmatrices)
      call linsol_initStructure (p_rsolverNode,ierror)
      call linsol_initData (p_rsolverNode,ierror)

      ! -----
      ! Solve for the X-displacement of the solid skeleton, uSx
      call lsysbl_createVecFromScalar(rvector1%RvectorBlock(1),rsingleRHS)
      call lsysbl_createVecFromScalar(rvector%RvectorBlock(1),rsingleSol)
      call lsysbl_duplicateVector (rsingleSol,rvector2,&
          LSYSSC_DUP_COPY,LSYSSC_DUP_EMPTY)
    
      call linsol_solveAdaptively (p_rsolverNode,rsingleSol,rsingleRHS,rvector2)
      
      if (p_rsolverNode%iresult .ne. 0) then
        ! Cancel, there is something wrong.
        call output_line('Cannot compute L2 projection, solver broke down. Using zero!')
        call lsysbl_clearVector (rvector)

      else
      
        ! -----
        ! Solve for the Y-displacement of the solid skeleton, uSy
        ! Throw away the temporary references to the 1st component in rvector1/rvector
        ! and replace them by references to the 2nd component.
        call lsysbl_releaseVector (rsingleRHS)
        call lsysbl_releaseVector (rsingleSol)

        if (rproblem%MSHOW_Initialisation .ge. 1) &
          call output_line('Solving L2 projection for U2...')

        call lsysbl_createVecFromScalar(rvector1%RvectorBlock(2),rsingleRHS)
        call lsysbl_createVecFromScalar(rvector%RvectorBlock(2),rsingleSol)

        call linsol_solveAdaptively (p_rsolverNode,rsingleSol,rsingleRHS,rvector2)
        
        if (p_rsolverNode%iresult .ne. 0) then
          ! Cancel, there is something wrong.
          call output_line(&
              'Cannot compute L2 projection, solver broke down. Using zero!')
          call lsysbl_clearVector (rvector)
          
        else

	  ! -----
	  ! Solve for the X-velocity of the solid skeleton.
	  ! Throw away the temporary references to the 1st component in rvector1/rvector
	  ! and replace them by references to the 2nd component.
	  call lsysbl_releaseVector (rsingleRHS)
	  call lsysbl_releaseVector (rsingleSol)

	  if (rproblem%MSHOW_Initialisation .ge. 1) &
	    call output_line('Solving L2 projection for U3...')

	  call lsysbl_createVecFromScalar(rvector1%RvectorBlock(3),rsingleRHS)
	  call lsysbl_createVecFromScalar(rvector%RvectorBlock(3),rsingleSol)

	  call linsol_solveAdaptively (p_rsolverNode,rsingleSol,rsingleRHS,rvector2)
        
	  if (p_rsolverNode%iresult .ne. 0) then
	    ! Cancel, there is something wrong.
	    call output_line(&
		'Cannot compute L2 projection, solver broke down. Using zero!')
	    call lsysbl_clearVector (rvector)
          
	  else

	    ! -----
	    ! Solve for the y-velocity of the solid skeleton.
	    ! Throw away the temporary references to the 1st component in rvector1/rvector
	    ! and replace them by references to the 2nd component.
	    call lsysbl_releaseVector (rsingleRHS)
	    call lsysbl_releaseVector (rsingleSol)

	    if (rproblem%MSHOW_Initialisation .ge. 1) &
	      call output_line('Solving L2 projection for U4...')

	    call lsysbl_createVecFromScalar(rvector1%RvectorBlock(4),rsingleRHS)
	    call lsysbl_createVecFromScalar(rvector%RvectorBlock(4),rsingleSol)

	    call linsol_solveAdaptively (p_rsolverNode,rsingleSol,rsingleRHS,rvector2)
        
	    if (p_rsolverNode%iresult .ne. 0) then
	      ! Cancel, there is something wrong.
	      call output_line(&
		  'Cannot compute L2 projection, solver broke down. Using zero!')
	      call lsysbl_clearVector (rvector)
          
	    else

	      ! -----
	      ! Solve for the x-velocity of the pore fluid.
	      ! Throw away the temporary references to the 1st component in rvector1/rvector
	      ! and replace them by references to the 2nd component.
	      call lsysbl_releaseVector (rsingleRHS)
	      call lsysbl_releaseVector (rsingleSol)

	      if (rproblem%MSHOW_Initialisation .ge. 1) &
		call output_line('Solving L2 projection for U5...')

	      call lsysbl_createVecFromScalar(rvector1%RvectorBlock(5),rsingleRHS)
	      call lsysbl_createVecFromScalar(rvector%RvectorBlock(5),rsingleSol)

	      call linsol_solveAdaptively (p_rsolverNode,rsingleSol,rsingleRHS,rvector2)
        
	      if (p_rsolverNode%iresult .ne. 0) then
		! Cancel, there is something wrong.
		call output_line(&
		    'Cannot compute L2 projection, solver broke down. Using zero!')
		call lsysbl_clearVector (rvector)
          
	      else

		! -----
		! Solve for the y-velocity of the pore fluid.
		! Throw away the temporary references to the 1st component in rvector1/rvector
		! and replace them by references to the 2nd component.
		call lsysbl_releaseVector (rsingleRHS)
		call lsysbl_releaseVector (rsingleSol)

		if (rproblem%MSHOW_Initialisation .ge. 1) &
		call output_line('Solving L2 projection for U6...')

		call lsysbl_createVecFromScalar(rvector1%RvectorBlock(6),rsingleRHS)
		call lsysbl_createVecFromScalar(rvector%RvectorBlock(6),rsingleSol)

		call linsol_solveAdaptively (p_rsolverNode,rsingleSol,rsingleRHS,rvector2)
        
		if (p_rsolverNode%iresult .ne. 0) then
		  ! Cancel, there is something wrong.
		  call output_line(&
		      'Cannot compute L2 projection, solver broke down. Using zero!')
		  call lsysbl_clearVector (rvector)
          
		else
		  ! -----
		  ! And the pressure. For that one we have to reinitialise the solver,
		  ! the matrices,...
		  call lsysbl_releaseVector (rsingleRHS)
		  call lsysbl_releaseVector (rsingleSol)
		  call lsysbl_releaseVector (rvector2)
		  call lsysbl_releaseMatrix (Rmatrices(1))

		  if (rproblem%MSHOW_Initialisation .ge. 1) &
		    call output_line('Solving L2 projection for P...')

		  ! Attach a new matrix.
		  call lsysbl_createMatFromScalar(&
		    rproblem%RlevelInfo(rproblem%NLMAX)%rasmTempl%rmatrixMassPressure,Rmatrices(1))

		  call linsol_doneData (p_rsolverNode,ierror)
		  call linsol_doneStructure (p_rsolverNode,ierror)

		  call linsol_setMatrices (p_rsolverNode,Rmatrices)
          
		  call linsol_initStructure (p_rsolverNode,ierror)
		  call linsol_initData (p_rsolverNode,ierror)
          
		  ! Create references to the pressure part and solve.
		  call lsysbl_createVecFromScalar(rvector1%RvectorBlock(7),rsingleRHS)
		  call lsysbl_createVecFromScalar(rvector%RvectorBlock(7),rsingleSol)
		  call lsysbl_duplicateVector (rsingleSol,rvector2,&
		      LSYSSC_DUP_COPY,LSYSSC_DUP_EMPTY)
          
		  call linsol_solveAdaptively (p_rsolverNode,rsingleSol,rsingleRHS,rvector2)
     
		  if (p_rsolverNode%iresult .ne. 0) then
		    ! Cancel, there is something wrong.
		    call output_line(&
			'Cannot compute L2 projection, solver broke down. Using zero!')
		    call lsysbl_clearVector (rvector)
		  end if ! for P
		end if	 ! for vFy
	      end if	 ! for vFx
	    end if	 ! for vSy
          end if	 ! for vSx
        end if		 ! for uSy
      end if		 ! for uSx
      
      ! -----
      ! That is it, cleanup.
      call lsysbl_releaseVector (rsingleRHS)
      call lsysbl_releaseVector (rsingleSol)
      call lsysbl_releaseVector (rvector2)
      call lsysbl_releaseMatrix (Rmatrices(1))
      
      call linsol_doneData (p_rsolverNode,ierror)
      call linsol_doneStructure (p_rsolverNode,ierror)
      call linsol_releaseSolver (p_rsolverNode)
    
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_writeSolution (rproblem,rvector,dtime)
  
!<description>
  ! Writes a solution vector rvector to a file as configured in the parameters
  ! in the DAT file.
!</description>

!<input>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(in) :: rproblem
  
  ! The solution vector to be written out. Must be set up according to the
  ! maximum level NLMAX in rproblem!
  type(t_vectorBlock), intent(in) :: rvector

  ! OPTIONAL: If present: Current time during a nonstationary simulation.
  ! For a stationary simulation and for the final timestep in a nonstationary
  ! simulation, this parameter must not be present.
  real(dp), intent(in), optional :: dtime
!</input>

!</subroutine>

    ! local variables
    integer :: iwriteSolutionLevel,cwriteFinalSolution
    type(t_vectorBlock) :: rvector1,rvector2
    type(t_vectorScalar) :: rvectorTemp
    character(LEN=SYS_STRLEN) :: sfile,sfileString
    integer :: ilev
    integer :: NEQ
    type(t_interlevelProjectionBlock) :: rprojection
    logical :: bformatted

    ! Get the parameter what to do with rvector
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'cwriteFinalSolution',cwriteFinalSolution,0)
    call parlst_getvalue_int (rproblem%rparamList,'CC-DISCRETISATION',&
                              'iwriteSolutionLevel',iwriteSolutionLevel,0)
    call parlst_getvalue_string (rproblem%rparamList,'CC-DISCRETISATION',&
                                 'swriteSolutionFilename',sfileString,'')

    if (cwriteFinalSolution .eq. 0) return ! nothing to do.
    
    ! Remove possible ''-characters
    read(sfileString,*) sfile
    
    bformatted = (cwriteFinalSolution .eq. 1) .or. (cwriteFinalSolution .eq. 3)
    ! level where to write out; correct if negative.
    if (iwriteSolutionLevel .le. 0) then
      iwriteSolutionLevel = rproblem%NLMAX-abs(iwriteSolutionLevel)
    end if

    if (iwriteSolutionLevel .lt. rproblem%NLMIN) then
      call output_line (&
          'Warning: Level for solution vector is < NLMIN! Writing out at level NLMIN!',&
          OU_CLASS_WARNING,OU_MODE_STD,'cc_initInitialSolution')
      iwriteSolutionLevel = rproblem%NLMIN
    end if
    
    ! Interpolate the solution down to level iinitialSolutionLevel.
    call lsysbl_copyVector (rvector,rvector1)   ! creates new rvector1!

    do ilev = rproblem%NLMAX,iwriteSolutionLevel+1,-1
      
      ! Initialise a vector for the lower level and a prolongation structure.
      call lsysbl_createVectorBlock (&
          rproblem%RlevelInfo(ilev-1)%rdiscretisation,rvector2,.false.)
      
      call mlprj_initProjectionVec (rprojection,rvector2)
      
      ! Interpolate to the next higher level.
      ! (Do not 'restrict'! Restriction would be for the dual space = RHS vectors!)

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

    call output_lbrk ()
    call output_line ('Writing RAW solution file: '//sfile)

    ! Write out the solution.
    if (bformatted) then
      if (.not. present(dtime)) then
        call vecio_writeBlockVectorHR (rvector1, 'SOLUTION', .true.,&
            0, sfile, '(E22.15)')
      else
        call vecio_writeBlockVectorHR (rvector1, 'SOLUTION', .true.,&
            0, sfile, '(E22.15)',"time="//trim(sys_sdL(dtime,15)))
      end if
    else
      call vecio_writeBlockVectorHR (rvector1, 'SOLUTION', .true.,0, sfile)
    end if

    ! Release temp memory.
    call lsysbl_releaseVector (rvector1)

  end subroutine

end module