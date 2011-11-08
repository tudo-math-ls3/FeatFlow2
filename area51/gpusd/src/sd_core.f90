module sd_core

  use fsystem
  use storage
  use genoutput
  use linearalgebra
  use linearsystemscalar
  use linearsystemblock
  use cubature
  use vectorio
  use matrixio
  use domainintegration
  use derivatives
  use statistics
  use basicgeometry
  use geometryaux
  use dofmapping
  use jumpstabilisation
  use feevaluation
  use triangulation
  use element
  use elementpreprocessing
  use transformation
  use spatialdiscretisation
  use bilinearformevaluation
  use collection
  
  implicit none

contains

  ! ***************************************************************************

!<subroutine>

  subroutine sd_asm ( u1Xvel,u1Yvel, rmatrix, &
                  dupsam,dnu,dalpha,dbeta,dtheta, ddelta, &
                  clocalh)
!<description>

!</description>

!<input>

  ! Primary X-velocity of $u_1$
  real(DP), dimension(:), intent(in) :: u1Xvel
  
  ! Primary Y-velocity of $u_1$
  real(DP), dimension(:), intent(in) :: u1Yvel

  ! dupsam  - control parameter.
  !          -1: simple upwind,
  !          =0: Samarskji upwind
  real(DP), intent(in) :: dupsam
  
  ! Viscosity parameter $\nu = 1/Re$ if viscosity is constant
  real(DP), intent(in) :: dnu
  
  ! Weighting factor for the mass matrix.
  real(DP), intent(in) :: dalpha

  ! Weighting factor for the Stokes matrix. (Stokes matrix = 1/Re * Laplace)
  real(DP), intent(in) :: dbeta

  ! Weighting factor of the convective operator: $\theta * u*grad(u)$.
  ! For time-dependent problems, this can be set to the step size
  ! in the $\Theta$-scheme.
  real(DP), intent(in) :: dtheta
  
  ! Weighting factor for the nonlinear term
  real(DP), intent(in) :: ddelta
  
  ! Method how to compute the local h
  integer, intent(in) :: clocalh
  
!</input>

!<inputoutput>
  ! The system matrix. Must be format 7 or 9.
  type(t_matrixScalar), intent(inout), target :: rmatrix
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: indof,indofALE,IEQ,IDOFE,JDOFE,icubp
  integer :: JCOL0,JDFG,JCOL
  integer :: IEL,IELset,IELmax
  logical, dimension(EL_MAXNDER) :: Bder,BderALE
  real(DP) :: dumax,dumaxr, du1loc, du2loc, dunorm,db,OM,AH,dre,dny
  real(DP) :: HBASI1,HBASI2,HBASI3,HBASJ1,HBASJ2,HBASJ3,HSUMI,HSUMJ
  integer :: NVE
  
  ! Matrix structure arrays
  integer, dimension(:), pointer :: p_Kcol
  integer, dimension(:), pointer :: p_Kld
  real(DP), dimension(:), pointer :: p_Da
  
  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  real(DP), dimension(:,:), allocatable :: p_DcubPtsRef

  ! The discretisation - for easier access
  type(t_spatialDiscretisation), pointer :: p_rdiscretisation
  
  ! Triangulation
  type(t_triangulation), pointer :: p_rtriangulation
  real(DP), dimension(:,:), pointer :: p_DvertexCoords
  integer, dimension(:,:), pointer :: p_IedgesAtElement,p_IverticesAtElement

  ! Number of elements in a block. Normally =BILF_NELEMSIM,
  ! except if there are less elements in the discretisation.
  integer :: nelementsPerBlock

  ! One and only element distribution
  type(t_elementDistribution), pointer :: p_relementDistribution

  ! For every cubature point on the reference element,
  ! the corresponding cubature weight
  real(DP), dimension(:), allocatable :: Domega
  
  ! number of cubature points on the reference element
  integer :: ncubp

  ! An element evaluation set for evaluating elements.
  type(t_evalElementSet) :: revalElementSet

  ! Arrays for saving Jacobian determinants
  real(DP), dimension(:,:), pointer :: p_Ddetj
  
  ! An allocateable array accepting the DOF`s of a set of elements.
  integer, dimension(:,:), allocatable, target :: Idofs
  
  ! Allocateable arrays for the values of the basis functions -
  ! for test and trial spaces.
  real(DP), dimension(:,:,:,:), allocatable, target :: Dbas

  ! Local matrices, used during the assembly.
  ! Values and positions of values in the global matrix.
  integer, dimension(:,:,:), allocatable :: Kentry
  real(DP), dimension(:,:), allocatable :: Dentry

  ! A pointer to an element-number list
  integer, dimension(:), pointer :: p_IelementList

  ! Pointer to the velocity field in the cubature points.
  real(DP), dimension(:,:,:), allocatable :: Dvelocity
  
  ! An array with local DELTA`s, each DELTA for one element
  real(DP), dimension(:), allocatable :: DlocalDelta

  ! Type of transformation from the reference to the real element
  integer(I32) :: ctrafoType
  
  ! Element evaluation tag; collects some information necessary for evaluating
  ! the elements.
  integer(I32) :: cevaluationTag
  
    ! Initialise the derivative flags
    Bder = .false.
    Bder(DER_FUNC) = .true.
    Bder(DER_DERIV_X) = .true.
    Bder(DER_DERIV_Y) = .true.

    ! For ALE we do not even need so much
    BderALE = .false.
    BderALE(DER_FUNC) = .true.
    
    ! Shortcut to the spatial discretisation
    p_rdiscretisation => rmatrix%p_rspatialDiscrTest
    
    ! Get the element distribution. Here, we can find information about
    ! the cubature formula etc...
    p_relementDistribution => p_rdiscretisation%RelementDistr(1)
    
    ! Get some information about the triangulation
    p_rtriangulation => p_rdiscretisation%p_rtriangulation
    call storage_getbase_double2d (p_rtriangulation%h_DvertexCoords,&
                                   p_DvertexCoords)
    call storage_getbase_int2d (p_rtriangulation%h_IverticesAtElement,&
                                p_IverticesAtElement)
    call storage_getbase_int2d (p_rtriangulation%h_IedgesAtElement,&
                                p_IedgesAtElement)
    
    ! Get the number of local DOF`s for trial/test functions.
    ! We assume trial and test functions to be the same.
    indof = elem_igetNDofLoc(p_relementDistribution%celement)

    ! Get the number of local DOF`s Q1 -- we need them for ALE.
    indofALE = elem_igetNDofLoc(p_relementDistribution%celement)
    
    ! Number of local DOF`s
    NVE = elem_igetNVE(p_relementDistribution%celement)
    
    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it is at most
    ! BILF_NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = min(BILF_NELEMSIM,p_rtriangulation%NEL)
    
    ! For cdef containing CONV_MODDEFECT, we build the defect vector
    !     D = RHS - A*U
    ! In this case, the defect(rhs vectors must be present
    
    
      ! Get matrix arrays
      call lsyssc_getbase_double (rmatrix,p_Da)
    call lsyssc_getbase_Kcol (rmatrix,p_Kcol)
    call lsyssc_getbase_Kld (rmatrix,p_Kld)
    
    ! Get from the trial element space the type of coordinate system
    ! that is used there:
    ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)
    
    ! Get the number of cubature points for the cubature formula
    ncubp = cub_igetNumPts(p_relementDistribution%ccubTypeBilForm)
    
    ! Allocate two arrays for the points and the weights
    allocate(Domega(ncubp))
    allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),ncubp))
    
    ! Get the cubature formula
    call cub_getCubature(p_relementDistribution%ccubTypeBilForm,p_DcubPtsRef, Domega)
    
    ! Allocate an array saving the coordinates of corner vertices of elements
    
    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  allocate(Dbas(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    ! would lead to nonused memory blocks in these arrays during the assembly,
    ! which reduces the speed by 50%!
    allocate(Dbas(indof,elem_getMaxDerivative(p_relementDistribution%celement), &
             ncubp,nelementsPerBlock))

    ! Allocate memory for the DOF`s of all the elements.
    allocate(Idofs(indof,nelementsPerBlock))
    
    ! Allocate memory for array with local DELTA`s
    allocate(DlocalDelta(nelementsPerBlock))

    ! Allocate an array saving the local matrices for all elements
    ! in an element set.
    ! We could also allocate EL_MAXNBAS*EL_MAXNBAS*BILF_NELEMSIM integers
    ! for this local matrix, but this would normally not fit to the cache
    ! anymore! indofTrial*indofTest*BILF_NELEMSIM is normally much smaller!
    allocate(Kentry(indof,indof,nelementsPerBlock))
    allocate(Dentry(indof,indof))
    
    ! Allocate memory for the velocity in the cubature points.
    allocate(Dvelocity(NDIM2D,ncubp,nelementsPerBlock))
    
    ! Initialisation of the element set.
    call elprep_init(revalElementSet)

    ! What is the reciprocal of nu? We need it later.
    if (dnu .ne. 0.0_DP) then
      dre = 1.0_DP/dnu
      
      ! dny gets the actual multiplier for the Laplace matrix.
      ! Remember: dbeta*Stokes = dbeta*dnu*Laplace = dny*Laplace.
      ! This may be =0.0 if the Stokes operator should not be included into
      ! the matrix.
      dny = dbeta*dnu
    else
      print *,'SD: NU=0 not allowed! Set dbeta=0 to prevent Stokes operator'// &
              ' from being build!'
      call sys_halt()
    end if

    ! If ddelta=0, we have to neglect the nonlinearity. In both cases,
    ! set DlocalDelta=0 which disables the nonlinear term in the assembly.
    ! If dupsam=0, we neglect the stabilisation term (central difference like
    ! discretisation), so we set DlocalDelta=0 as well.
    if ((ddelta .eq. 0.0_DP) .or. (dupsam .eq. 0.0_DP)) then
      call lalg_clearVectorDble (DlocalDelta)
    end if
    
    ! Calculate the maximum norm of the actual velocity field
    ! U = A1*U1 + A2*U2 into DUMAX.
    ! Round up the norm to 1D-8 if it is too small...

    dumax=0.0_DP
      do IEQ=1,size(u1Xvel)
        du1loc = u1Xvel(IEQ)
        du2loc = u1Yvel(IEQ)
        dunorm = sqrt(du1loc**2+du2loc**2)
        dumax = max(DUMAX,DUNORM)
      end do

    if (dumax.lt.1E-8_DP) dumax=1E-8_DP
    dumaxr = 1.0_DP/dumax

    ! p_IelementList must point to our set of elements in the discretisation
    ! with that combination of trial/test functions
    call storage_getbase_int (p_relementDistribution%h_IelementList, &
                              p_IelementList)
    
    ! Loop over the elements - blockwise.
    do IELset = 1, size(p_IelementList), BILF_NELEMSIM
    
      ! We always handle BILF_NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most BILF_NELEMSIM
      ! elements simultaneously.
      
      IELmax = min(size(p_IelementList),IELset-1+BILF_NELEMSIM)
    
      ! The outstanding feature with finite elements is: A basis
      ! function for a DOF on one element has common support only
      ! with the DOF`s on the same element! E.g. for Q1:
      !
      !        #. . .#. . .#. . .#
      !        .     .     .     .
      !        .  *  .  *  .  *  .
      !        #-----O-----O. . .#
      !        |     |     |     .
      !        |     | IEL |  *  .
      !        #-----X-----O. . .#
      !        |     |     |     .
      !        |     |     |  *  .
      !        #-----#-----#. . .#
      !
      ! --> On element IEL, the basis function at "X" only interacts
      !     with the basis functions in "O". Elements in the
      !     neighbourhood ("*") have no support, therefore we only have
      !     to collect all "O" DOF`s.
      !
      ! Calculate the global DOF`s into IdofsTrial / IdofsTest.
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF`s of our BILF_NELEMSIM elements simultaneously.
      call dof_locGlobMapping_mult(p_rdiscretisation, p_IelementList(IELset:IELmax), &
                                  Idofs)
                                  
      ! Calculate local DELTA`s for streamline diffusion method.
      ! (cf. p. 121 in Turek`s CFD book).
      ! For every element, we need a local DELTA.
      ! Every local delta is weighted by the global "ddelta".
      ! If ddelta=0, we do not do anything as this disables the
      ! nonlinear term.
      ! If UPSAM=0.0, we have a central-difference like discretisation, which
      ! is one can see as the local stabilisation weight Delta is also = 0.0.
      ! In this case, we even switch of the calculation of the local Delta,
      ! as it is always =0.0, so we save a little bit time.
      if ((ddelta .ne. 0.0_DP) .and. (dupsam .ne. 0.0_DP))then
        call getLocalDeltaQuad (clocalh,&
                      u1Xvel,u1Yvel,&
                      p_IelementList(IELset:IELmax),&
                      duMaxR,DlocalDelta,p_rtriangulation,Idofs,dupsam,dre)
      end if
                                   
      ! For the assembly of the global matrix, we use a "local"
      ! approach. At first we build a "local" system matrix according
      ! to the current element. This contains all additive
      ! contributions of element IEL, which are later added at the
      ! right positions to the elements in the global system matrix.
      !
      ! We have indofTrial trial DOF`s per element and
      ! indofTest test DOF`s per element. Therefore there are
      ! indofTrial*indofTest tupel of basis-/testfunctions (phi_i,psi_j)
      ! "active" (i.e. have common support) on our current element, each
      ! giving an additive contribution to the system matrix.
      !
      ! We build a quadratic indofTrial*indofTest local matrix:
      ! Kentry(1..indofTrial,1..indofTest) receives the position
      !   in the global system matrix, where the corresponding value
      !   has to be added to.
      ! (The corresponding contrbutions can be saved separately,
      !  but we directly add them to the global matrix in this
      !  approach.)
      !
      ! We build local matrices for all our elements
      ! in the set simultaneously.
      ! Loop through elements in the set and for each element,
      ! loop through the local matrices to initialise them:
      do IEL=1,IELmax-IELset+1
      
        ! For building the local matrices, we have first to
        ! loop through the test functions (the "O"`s), as these
        ! define the rows in the matrix.
        do IDOFE=1,indof
        
          ! Row IDOFE of the local matrix corresponds
          ! to row=global DOF KDFG(IDOFE) in the global matrix.
          ! This is one of the the "O"`s in the above picture.
          ! Get the starting position of the corresponding row
          ! to JCOL0:

          JCOL0=p_KLD(Idofs(IDOFE,IEL))
          
          ! Now we loop through the other DOF`s on the current element
          ! (the "O"`s).
          ! All these have common support with our current basis function
          ! and will therefore give an additive value to the global
          ! matrix.
          
          do JDOFE=1,indof
            
            ! Get the global DOF of the "X" which interacts with
            ! our "O".
            
            JDFG=Idofs(JDOFE,IEL)
            
            ! Starting in JCOL0 (which points to the beginning of
            ! the line initially), loop through the elements in
            ! the row to find the position of column IDFG.
            ! Jump out of the do loop if we find the column.
            
            do JCOL=JCOL0,rmatrix%NA
              if (p_KCOL(JCOL) .eq. JDFG) exit
            end do

            ! Because columns in the global matrix are sorted
            ! ascendingly (except for the diagonal element),
            ! the next search can start after the column we just found.
            
            ! JCOL0=JCOL+1
            
            ! Save the position of the matrix entry into the local
            ! matrix.
            ! Note that a column in Kentry corresponds to a row in
            ! the real matrix. We aligned Kentry/DENTRY this way to get
            ! higher speed of the assembly routine, since this leads
            ! to better data locality.
            
            Kentry(JDOFE,IDOFE,IEL)=JCOL
            
          end do ! IDOFE
          
        end do ! JDOFE
        
      end do ! IEL
      
      ! Ok, we found the positions of the local matrix entries
      ! that we have to change.
      ! To calculate the matrix contributions, we have to evaluate
      ! the elements to give us the values of the basis functions
      ! in all the DOF`s in all the elements in our set.
      !
      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag.
      cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)
                      
      ! In the first loop, calculate the coordinates on the reference element.
      ! In all later loops, use the precalculated information.
      if (IELset .eq. 1) then
        cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REFPOINTS)
      else
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
      end if

      ! Calculate all information that is necessary to evaluate the finite element
      ! on all cells of our subset. This includes the coordinates of the points
      ! on the cells.
      call elprep_prepareSetForEvaluation (revalElementSet,&
          cevaluationTag, p_rtriangulation, p_IelementList(IELset:IELmax), &
          ctrafoType, p_DcubPtsRef(:,1:ncubp))
      p_Ddetj => revalElementSet%p_Ddetj

      ! Calculate the values of the basis functions.
      call elem_generic_sim2 (p_relementDistribution%celement, &
          revalElementSet, Bder, Dbas)

      ! We want to set up the nonlinear part of the matrix
      !
      !   n~_h (u_h, u_h, v_h)
      !
      ! = n_h (u_h, u_h, v_h) + sum_T ( delta_T ( u_h*grad u_h, u_h*grad v_h)_T )
      !   ^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      !  standard nonlin. part                  stabilization
      !
      ! More precisely, as we want to assemble the matrix which is
      ! later multiplied with coefficient vectors, we have to insert
      ! basis functions in the above terms instead of u_h and v_h.
      ! Assuming the representation u_h=sum_j(u_j*Phi_j) and
      ! v_h=sum_i(u_i,Phi_i), the above term is evaluated in the
      ! DOF`s as:
      !
      !   n_h (u_h, Phi_j, Phi_i)
      ! + sum_T ( delta_T ( u_h*grad Phi_j, u_h*grad Phi_i )_T )
      !
      ! In nonstationary simulations, the system matrix typically
      ! contains a mass matrix to respect the time derivative.
      ! The matrix has the form
      !
      ! [  dcmass*M*I  +  THWEIG * (-nu * Laplace(.))  ] + THWEIG * u grad(.)
      !
      ! In a first step, we calculate the velocity field in all
      ! cubature points on all elements of the current block.
      ! If we only have a primary velocity field
      ! (dweight2=0), we can calculate that only by summing up the
      ! velocities in U1Lx, otherwise we have to sum up
      ! dweight1*u1vel + dweight2*u2vel
      
        ! Loop over all elements in the current set
        do IEL=1,IELmax-IELset+1
        
          ! Loop over all cubature points on the current element
          do ICUBP = 1, ncubp
          
            du1loc = 0.0_DP
            du2loc = 0.0_DP
          
            ! Perform a loop through the trial DOF`s.
            do JDOFE=1,indof

              ! Get the value of the (test) basis function
              ! phi_i (our "O") in the cubature point:
              db = Dbas(JDOFE,1,ICUBP,IEL)
              
              ! Sum up to the value in the cubature point
              JDFG = Idofs(JDOFE,IEL)
              du1loc = du1loc + u1Xvel(JDFG)*db
              du2loc = du2loc + u1Yvel(JDFG)*db

            end do ! JDOFE
            
            ! Save the computed velocity
            Dvelocity(1,ICUBP,IEL) = du1loc
            Dvelocity(2,ICUBP,IEL) = du2loc
          
          end do ! ICUBP
          
        end do ! IEL
          

      ! Ok, we now use Dvelocity as coefficient array in the assembly
      ! of a biliinear form!
      !
      ! Loop over the elements in the current set.

      do IEL=1,IELmax-IELset+1
        
        ! Clear the local matrix
        Dentry = 0.0_DP
        
        ! Loop over all cubature points on the current element
        do ICUBP = 1, ncubp

          ! Calculate the current weighting factor in the cubature formula
          ! in that cubature point.
          !
          ! Normally, we have to take the absolut value of the determinant
          ! of the mapping here!
          ! In 2D, the determinant is always positive, whereas in 3D,
          ! the determinant might be negative -- that is normal!
          ! But because this routine only works in 2D, we can skip
          ! the ABS here!

          OM = Domega(ICUBP)*p_Ddetj(ICUBP,IEL)

          ! Current velocity in this cubature point:
          du1loc = Dvelocity (1,ICUBP,IEL)
          du2loc = Dvelocity (2,ICUBP,IEL)
          
          ! We take a more detailed look onto the last scalar product
          ! of n~_h (u_h, u_h, v_h) what we want to calculate here.
          !
          ! The vector u_h=(DU1,DU2) contains both velocity components,
          ! for the X as well as for the Y velocity. On the other hand
          ! the system matrix we want to build here will be designed for
          ! one velocity component only! Therefore, Phi_i and Phi_j
          ! are scalar functions, so grad(Phi_i), grad(Phi_j) are vectors
          ! with two components. Therefore, the last scalar product is more
          ! in detail:
          !
          !     ( u_h*grad Phi_j, u_h*grad Phi_i )_T
          !
          ! =   ( < (DU1) , (grad(Phi_j)_1) > , < (DU1) , (grad(Phi_i)_1) > )_T
          !         (DU2) , (grad(Phi_j)_2)       (DU2) , (grad(Phi_i)_2)
          !
          ! =   < (DU1) , (grad(Phi_j)_1) >  *  < (DU1) , (grad(Phi_j)_1) >
          !       (DU2) , (grad(Phi_j)_2)         (DU2) , (grad(Phi_j)_2)
          !
          ! =   HSUMJ * HSUMI
          !
          ! i.e. a product of two scalar values!
          !
          ! Summing up over all pairs of multiindices.
          !
          ! Outer loop over the DOF`s i=1..indof on our current element,
          ! which corresponds to the basis functions Phi_i:

          do IDOFE=1,indof
          
            ! Fetch the contributions of the (test) basis functions Phi_i
            ! (our "O")  for function value and first derivatives for the
            ! current DOF into HBASIy:
          
            HBASI1 = Dbas(IDOFE,1,ICUBP,IEL)
            HBASI2 = Dbas(IDOFE,2,ICUBP,IEL)
            HBASI3 = Dbas(IDOFE,3,ICUBP,IEL)
           
            ! Calculate
            !
            !     U * grad(Phi_i)  =  < grad(Phi_i), U >
            !
            !   = ( grad(Phi_i)_1 , (DU1) )
            !     ( grad(Phi_i)_2   (DU2) )
            !
            ! Remember: DU1MV=DU2MV=0 in this case.
            !
            ! If ALE is active, use v=mesh velocity and calculate
            !
            !     (U-v) * grad(Phi_i)  =  < grad(Phi_i), U-v >
            !
            !   = ( grad(Phi_i)_1 , (DU1-DU1MV) )
            !     ( grad(Phi_i)_2   (DU2-DU2MV) )

            HSUMI = HBASI2*du1loc + HBASI3*du2loc

            ! Inner loop over the DOF`s j=1..indof, which corresponds to
            ! the basis function Phi_j:

            do JDOFE=1,indof
              
             
                ! Fetch the contributions of the (trial) basis function Phi_j
                ! (out "X") for function value and first derivatives for the
                ! current DOF into HBASJy:
              
                HBASJ1 = Dbas(JDOFE,1,ICUBP,IEL)
                HBASJ2 = Dbas(JDOFE,2,ICUBP,IEL)
                HBASJ3 = Dbas(JDOFE,3,ICUBP,IEL)

                ! Calculate
                !
                !     U * grad(Phi_j)  =  < grad(Phi_j), U >
                !
                !   = ( grad(Phi_j)_1 , (DU1) )
                !     ( grad(Phi_j)_2   (DU2) )
                !
                ! Remember: DU1MV=DU2MV=0 in this case.
                !
                ! If ALE is active, use v=mesh velocity and calculate
                !
                !     (U-v) * grad(Phi_j)  =  < grad(Phi_j), U-v >
                !
                !   = ( grad(Phi_j)_1 , (DU1-DU1MV) )
                !     ( grad(Phi_j)_2   (DU2-DU2MV) )
                !
                ! But as v is already incorporated into DVelocity,
                ! we do not have to worry about that.

                HSUMJ = HBASJ2*du1loc+HBASJ3*du2loc
    
                ! Finally calculate the contribution to the system
                ! matrix. Depending on the configuration of DNU,
                ! dalpha,ddelta,... this decomposes into three
                ! different parts:
                !
                ! AH = n~_h(u_h,phi_j,phi_i)        | nonlinear part
                !    + dny*(grad(phi_j,grad(phi_i)) | -dny*Laplace(u) = -dbeta*Stokes
                !    + dalpha*(phi_j*phi_i)         | Mass matrix
                !
                ! The last two parts are probably not added to the
                ! matrix by setting DNY or CT0 to 0, respectively.
                !
                ! For saving some numerical operations, we write:
                !
                !     HSUMJ * (Delta * HSUMI + HBASI1)
                !
                ! =   Delta * HSUMJ * HSUMI
                !   + HSUMJ * HBASI1
                !
                ! =   Delta * ( U*grad(Phi_j), U*grad(Phi_i) )
                !   + (U*grad(Phi_j),Phi_i)
                !
                ! <->   sum_T ( delta_T ( u_h*grad Phi_j, u_h*grad Phi_i) )_T
                !     + n_h (u_h, Phi_j, Phi_i)
                !
                ! plus the terms for the Stokes and Mass matrix,
                ! if their coefficient is <> 0.
                
                AH = ddelta*HSUMJ*(DlocalDelta(IEL)*HSUMI+HBASI1) &
                    + dny*(HBASI2*HBASJ2+HBASI3*HBASJ3) &
                    + dalpha*HBASI1*HBASJ1
    
              ! Weighten the calculated value AH by the cubature
              ! weight OM and add it to the local matrix. After the
              ! loop over all DOF`s is finished, each entry contains
              ! the calculated integral.

              Dentry(JDOFE,IDOFE) = Dentry(JDOFE,IDOFE)+OM*AH
              
            end do ! IDOFE
            
          end do ! JDOFE

        end do ! ICUBP
        
        ! Now we have set up a "local" system matrix. We can either
        ! include it into the real matrix or we can use it to simply
        ! modify the RHS vector to create a defect vector (throwing
        ! away the information about the matrix afterwards, which would
        ! result in a matrix free modification of the RHS vector).
        !
        ! For cdef= containing CONV_MODMATRIX, incorporate our "local" system matrix
        ! into the global matrix. The position of each entry DENTRY(X,Y)
        ! in the global matrix array A was saved in element Kentry(X,Y)
        ! before.
        ! Kentry gives the position of the additive contributions in Dentry.
        ! The entry is weighted by the current dtheta, which is usually
        ! the weighting parameter of the corresponding THETA-scheme of a
        ! nonstationary simulation. For stationary simulations, dtheta is typically
        ! 1.0 which includes the local matrix into the global one directly.)
        
          do IDOFE=1,indof
            do JDOFE=1,indof
              p_DA(Kentry(JDOFE,IDOFE,IEL)) = p_DA(Kentry(JDOFE,IDOFE,IEL)) + &
                 Dentry(JDOFE,IDOFE)
            end do
          end do

      end do ! IEL

    end do ! IELset
    
    ! Release memory
    call elprep_releaseElementSet(revalElementSet)

    deallocate(p_DcubPtsRef)
    deallocate(Domega)
    deallocate(DlocalDelta)
    deallocate(Dvelocity)
    deallocate(Dentry)
    deallocate(Kentry)
    deallocate(Idofs)
    deallocate(Dbas)

  end subroutine

  ! ----------------------------------------------------------------------

  subroutine getLocalDeltaQuad (clocalh,&
                      Du1x,Du1y,Ielements,&
                      duMaxR,Ddelta,rtriangulation,Idofs,dupsam,dnurec)

  ! This routine calculates a local ddelta=DELTA_T for a set of finite
  ! elements Ielements. This can be used by the streamline diffusion
  ! stabilisation technique as a multiplier of the (local) bilinear form.
  !
  ! The effective velocity that is used for calculating the ddelta
  ! is combined by a weighted mean of the two velocity fields U1,U2
  ! by:
  !                   du = da1*Du1 + da2*Du2
  ! The coefficients A1,A2 allow the caller to take influence on which
  ! velocity field to weight more.
  
  ! Method how to compute the local h.
  ! =0: Use the root of the area of the element as local H
  ! =1: Use the length of the way that a particle travels through
  !     the element in direction of the flow
  integer, intent(in) :: clocalH
  
  ! Main velocity field.
  real(DP), dimension(*), intent(in) :: du1x,du1y
  
  ! Reciprocal of the maximum norm of velocity in the domain:
  ! 1/duMaxR = 1/||u||_Omega
  real(DP), intent(in) :: duMaxR
  
  ! Reciprocal value 1/NU of coefficient NU in front of the
  ! Laplacian term of the Navier-Stokes equation
  !   NU * Laplace(u) + u*grad(u) + ...
  real(DP), intent(in) :: dnuRec
  
  ! user defined parameter for configuring the streamline diffusion.
  ! < 0: Simple calculation of ddelta, using
  !      ddelta = |UPSAM| * h_T.
  ! > 0: usually UPSAM = 0.1 .. 2; Samarskji-like calculation of ddelta using:
  !      ddelta = UPSAM * h_t/||u||_T * 2*Re_T/(1+Re_T)
  real(DP), intent(in) :: dupsam
  
  ! List of elements where the Ddelta should be calculated
  integer, dimension(:), intent(in) :: Ielements
  
  ! Array with global degrees of freedom on the elements
  integer, dimension(:,:), intent(in) :: Idofs
  
  ! Triangulation that defines the mesh.
  type(t_triangulation), intent(in) :: rtriangulation

  ! Out: local Ddelta on all elements
  real(DP), dimension(:), intent(out) :: ddelta

  ! local variables
  real(DP) :: dlocalH,du1,du2,dunorm,dreLoc
  integer :: iel,ielidx
  integer :: idof
  integer, dimension(:,:), pointer :: p_IverticesAtElement
  real(DP), dimension(:,:), pointer :: p_DvertexCoords
  real(DP), dimension(:), pointer :: p_DelementVolume

    ! Get some crucial data
    if (clocalh .eq. 0) then
      call storage_getbase_double (rtriangulation%h_DelementVolume,p_DelementVolume)
      
      ! Loop through all elements
      do ielidx = 1,size(Ielements)
      
        iel = Ielements(ielidx)

        ! Loop through the local degrees of freedom on element IEL.
        ! Sum up the velocities on these DOF`s. This will result
        ! in the vector (DU1,DU2) representing the (mean) X/Y-velocity
        ! through element IEL.

        ! For elements whose DOF`s represent directly the velocity, U1/U2
        ! represent the mean velocity
        ! along an egde/on the midpoint of each edge, so U1/U2 is
        ! clearly an approximation to the velocity in element T.

        du1=0.0_DP
        du2=0.0_DP
        do idof=1,ubound(Idofs,1)
          du1=du1+du1x(Idofs(idof,ielidx))
          du2=du2+du1y(Idofs(idof,ielidx))
        end do

        ! Calculate the norm of that local velocity:

        dunorm = sqrt(du1**2+du2**2) / real(ubound(Idofs,1),DP)
        
        ! Now we have:   dunorm = ||u||_T
        ! and:           u_T = a1*u1_T + a2*u2_T

        ! If the norm of the velocity is small, we choose ddelta = 0,
        ! which results in central difference in the streamline diffusion
        ! matrix assembling:

!        if (dunorm .le. 1E-8_DP) then
!
!          Ddelta(ielidx) = 0.0_DP
!
!        else

          ! Calculate the local h from the area of the element
          dlocalH = sqrt(p_DelementVolume(iel))

          ! Calculate ddelta... (cf. p. 121 in Turek`s CFD book)

!          if (dupsam .lt. 0.0_DP) then
!
!            ! For UPSAM<0, we use simple calculation of ddelta:
!
!            Ddelta(ielidx) = abs(dupsam)*dlocalH
!
!          else
          
            ! For UPSAM >= 0, we use standard Samarskji-like calculation
            ! of ddelta. At first calculate the local Reynolds number
            ! RELOC = Re_T = ||u||_T * h_T / NU
            
            dreLoc = dunorm*dlocalH*dnuRec
            
            ! and then the ddelta = UPSAM * h_t/||u|| * 2*Re_T/(1+Re_T)
            
            Ddelta(ielidx) = dupsam * dlocalH*duMaxR * 2.0_DP*(dreLoc/(1.0_DP+dreLoc))
            
!          end if ! (UPSAM.LT.0.0)
          
!        end if ! (dunorm.LE.1D-8)

      end do
      
    else
    
      call storage_getbase_double2d (rtriangulation%h_DvertexCoords,p_DvertexCoords)
      call storage_getbase_int2d (rtriangulation%h_IverticesAtElement,p_IverticesAtElement)

      ! Loop through all elements
      do ielidx = 1,size(Ielements)
      
        iel = Ielements(ielidx)

        ! Loop through the local degrees of freedom on element IEL.
        ! Sum up the velocities on these DOF`s. This will result
        ! in the vector (DU1,DU2) representing the (mean) X/Y-velocity
        ! through element IEL.

        ! For elements whose DOF`s represent directly the velocity, U1/U2
        ! represent the mean velocity
        ! along an egde/on the midpoint of each edge, so U1/U2 is
        ! clearly an approximation to the velocity in element T.

        du1=0.0_DP
        du2=0.0_DP
        do idof=1,ubound(Idofs,1)
          du1=du1+du1x(Idofs(idof,ielidx))
          du2=du2+du1y(Idofs(idof,ielidx))
        end do

        ! Calculate the norm of that local velocity:

        dunorm = sqrt(du1**2+du2**2) / real(ubound(Idofs,1),DP)
        
        ! Now we have:   dunorm = ||u||_T
        ! and:           u_T = a1*u1_T + a2*u2_T

        ! If the norm of the velocity is small, we choose ddelta = 0,
        ! which results in central difference in the streamline diffusion
        ! matrix assembling:

        if (dunorm .le. 1E-8_DP) then
        
          Ddelta(ielidx) = 0.0_DP

        else

          ! u_T defines the "slope" of the velocity through
          ! the element T. At next, calculate the local mesh width
          ! dlocalH = h = h_T on our element T=IEL:

          call getLocalMeshWidthQuad (dlocalH,dunorm, du1, du2, iel, &
              p_IverticesAtElement,p_DvertexCoords)

          ! Calculate ddelta... (cf. p. 121 in Turek`s CFD book)

          if (dupsam .lt. 0.0_DP) then

            ! For UPSAM<0, we use simple calculation of ddelta:
          
            Ddelta(ielidx) = abs(dupsam)*dlocalH
            
          else
          
            ! For UPSAM >= 0, we use standard Samarskji-like calculation
            ! of ddelta. At first calculate the local Reynolds number
            ! RELOC = Re_T = ||u||_T * h_T / NU
            
            dreLoc = dunorm*dlocalH*dnuRec
            
            ! and then the ddelta = UPSAM * h_t/||u|| * 2*Re_T/(1+Re_T)
            
            Ddelta(ielidx) = dupsam * dlocalH*duMaxR * 2.0_DP*(dreLoc/(1.0_DP+dreLoc))
            
          end if ! (UPSAM.LT.0.0)
          
        end if ! (dunorm.LE.1D-8)

      end do

    end if

  end subroutine

  ! ----------------------------------------------------------------------

  pure subroutine getLocalMeshWidthQuad (dlocalH, dunorm,  XBETA1, &
                      XBETA2, JEL,Kvert,Dcorvg)
  
  ! Determine the local mesh width for an element JEL of a
  ! triangulation.
  
  ! Element where the local h should be calculated
  integer, intent(in)               :: JEL
  
  integer, dimension(TRIA_MAXNVE2D,*), intent(in) :: Kvert
  real(DP), dimension(NDIM2D,*), intent(in)          :: Dcorvg
  
  ! norm ||u||_T = mean velocity through element T=JEL
  real(DP), intent(in)  :: dunorm
  
  ! mean velocity u_T = (xbeta1,xbeta2) through element T=JEL
  real(DP), intent(in)  :: XBETA1, XBETA2
  
  ! local mesh width
  real(DP), intent(out) :: dlocalH
  
  ! local variables
  real(DP) :: dlambda
  integer :: NECK1,NECK2,NECK3,NECK4
  real(DP) :: X1,Y1,X2,Y2,X3,Y3,X4,Y4
  real(DP) :: dalphaMax, dalpha

    ! Fetch the numbers of the four corners of element JEL

    neck1=Kvert(1,JEL)
    neck2=Kvert(2,JEL)
    neck3=Kvert(3,JEL)
    neck4=Kvert(4,JEL)

    ! Fetch the coordinates of these corners

    x1=Dcorvg(1,neck1)
    y1=Dcorvg(2,neck1)
    x2=Dcorvg(1,neck2)
    y2=Dcorvg(2,neck2)
    x3=Dcorvg(1,neck3)
    y3=Dcorvg(2,neck3)
    x4=Dcorvg(1,neck4)
    y4=Dcorvg(2,neck4)

    ! Scale: (deactivated)

    !  dsp=max(xbeta1,xbeta2)

    !  xbeta1=xbeta1
    !  xbeta2=xbeta2

    dalphaMax=0.0_DP
    
    ! In the next step, we calculate the `maximum possible mesh with
    ! in direction of the flow`; this is the maximum possible length
    ! that a particle can cross in the current element.
    ! The picture in mind is the following:
    !
    !          G3
    !   +-------------X-------+
    !   |            /        |
    !   |           /         |
    !   |          /          |
    !   |         /           |
    !   |        /            |
    ! G4|       /             | G2
    !   |      ^ (beta1,beta2)|
    !   |     /               |
    !   |    /                |
    !   |   /                 |
    !   |  /                  |
    !   | /                   |
    !   |/                    |
    !   O---------------------+
    !            G1
    !
    ! The vector (beta1,beta2) gives the direction of the flow.
    ! A particle starting in point O and moves at most up to point X.
    ! The length of the line (O,X) is the local mesh with h.
    !
    ! Loop through the four corners of element JEL and check
    ! of a line with slope BETA=(xbeta1,xbeta2) starting in this
    ! corner really intersects with one of the edges of the element.
    ! Remark that we only have to check the two opposite edges
    ! to the current corner!

    ! -----------------------------------------------------------------
    ! Check the first corner:

    call intersectLines2D(X1,Y1,dalpha,XBETA1,XBETA2, &
                X3,Y3,dlambda,X2,Y2)
    dalphaMax=max(dalpha,dalphaMax)

    call intersectLines2D(X1,Y1,dalpha,XBETA1,XBETA2, &
                X3,Y3,dlambda,X4,Y4)
    dalphaMax=max(dalpha,dalphaMax)
    
    ! -----------------------------------------------------------------
    ! The second one...
    
    call intersectLines2D(X2,Y2,dalpha,XBETA1,XBETA2, &
                X4,Y4,dlambda,X1,Y1)
    dalphaMax=max(dalpha,dalphaMax)

    call intersectLines2D(X2,Y2,dalpha,XBETA1,XBETA2, &
                X4,Y4,dlambda,X3,Y3)
    dalphaMax=max(dalpha,dalphaMax)
    
    ! -----------------------------------------------------------------
    ! The third one...
    
    call intersectLines2D(X3,Y3,dalpha,XBETA1,XBETA2, &
                X1,Y1,dlambda,X2,Y2)
    dalphaMax=max(dalpha,dalphaMax)

    call intersectLines2D(X3,Y3,dalpha,XBETA1,XBETA2, &
                X1,Y1,dlambda,X4,Y4)
    dalphaMax=max(dalpha,dalphaMax)
    
    ! -----------------------------------------------------------------
    ! And the fourth=last one...
    
    call intersectLines2D(X4,Y4,dalpha,XBETA1,XBETA2, &
                X2,Y2,dlambda,X1,Y1)
    dalphaMax=max(dalpha,dalphaMax)

    call intersectLines2D(X4,Y4,dalpha,XBETA1,XBETA2, &
                X2,Y2,dlambda,X3,Y3)
    dalphaMax=max(dalpha,dalphaMax)

    ! -----------------------------------------------------------------
    ! finally determine the local h=h_T
    !
    ! dalphaMax is the maximum alpha, normalised as 'parameter value',
    ! i.e. dalphaMax=1.0 corresponds to a vector 1.0*(dbeta1,dbeta2).
    ! We multiply with dunorm=|(dbeta1,dbeta2)| to get the actual length
    ! of the vector which can be placed inside of the element.
    !
    ! Furthermore, we multiply with an additional weight 4. (why ?!?)

    dlocalH=dalphaMax*4.0_DP*dunorm

  end subroutine
  
  ! ----------------------------------------------------------------------

  pure subroutine intersectLines2D (XO,YO,dalpha,BETA1,BETA2, &
                      XA,YA,dlambda,XB,YB)

  ! Intersect two lines in R^2

  ! Origin of line 1
  real(DP), intent(in) :: XO,YO
  
  ! Direction of line 1
  real(DP), intent(in) :: BETA1,BETA2
  
  ! One point on the second line
  real(DP), intent(in) :: XA,YA
  
  ! Another point on the second line
  real(DP), intent(in) :: XB,YB
  
  ! Parameter value of the intersection point on line 1.
  ! =0.0, if there is no intersection point
  real(DP), intent(out) :: dalpha
  
  real(DP), intent(out) :: dlambda
  
  ! local variables
  double precision :: dsp

    ! Scalar product of the line (xa,ya)->(xb,yb) with the
    ! counterclockwise normal n1 of (beta1,beta2)
    dsp=BETA2*(XB-XA)-BETA1*(YB-YA)
    
    if (dsp.eq.0.0_DP) then
    
      ! beta and the vector are parallel
      dalpha=0.0_DP
      
    else

      ! Scalar product of (beta1,beta2) with the (inner) normal vector n2
      ! of the line (xo,yo)->(xa,ya).
      dlambda=(BETA1*(YA-YO)-BETA2*(XA-XO))/dsp

      !                    (xb,yb)
      !   +-----------------+
      !   |                 |
      !   |                 |
      !   ^ n2              |
      !   !                 |
      !   !  (beta1,beta2)  |    (beta1,beta2)
      !   !    ^            |    ^
      !   !   /  ^__ n1     |   /
      !   !  /      \__     |  /
      !   ! /          \__  | /
      !   !/              \_|/
      !   +-----------------+
      ! (xo,yo)            (xa,ya)
      !
      ! (What is this? Documentation incomplete. Has someone a good
      ! reference?)

      ! is the intersection point inside of the element?
      if ((dlambda.ge.-1E-1_DP).and.(dlambda.le.1.11E0_DP)) then
        if (BETA1 .ne. 0.0_DP) then
          dalpha=((XA-XO)+dlambda*(XB-XA))/BETA1
        else
          if (BETA2 .ne. 0.0_DP) then
            dalpha=((YA-YO)+dlambda*(YB-YA))/BETA2
          else
            dalpha=0.0_DP
          end if
        end if
      else
        dalpha=0.0_DP
      end if
      
    end if

  end subroutine

end module