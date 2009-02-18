!##############################################################################
!# ****************************************************************************
!# <name> stokes2d_callback </name>
!# ****************************************************************************
!#
!# <purpose>
!#
!# Contains extended assembly routines for the matrix in the optimal control
!# of distributed fluid flow.
!#
!# </purpose>
!##############################################################################

module optcontrolconvection

  use fsystem
  use storage
  use linearsystemscalar
  use linearsystemblock
  use derivatives
  use bilinearformevaluation
  
  implicit none

  ! Defines the terms in the Optimal-Control operator to compute.

  type t_optcoperator
  
    ! Stabilisation parameter.
    ! Standard value = 1.0_DP
    ! Note: A value of 0.0_DP sets up the convection part without
    ! any stabilisation (central-difference like discretisation).
    real(DP) :: dupsam = 1.0_DP
    
    ! Whether the viscosity is constant.
    logical :: bconstViscosity = .true.
    
    ! Viscosity parameter $\nu = 1/Re$ if viscosity is constant
    real(dp) :: dnu = 0.0_DP
    
    ! Coefficient in front of the Mass matrix in the primal equation
    real(dp) :: dprimalAlpha = 0.0_DP
  
    ! Coefficient in front of the Stokes operator nu*Laplace in the primal equation
    real(dp) :: dprimalBeta = 0.0_DP
    
    ! Coefficient in front of the convection operator \grad(.) y in the
    ! primal equation.
    real(dp) :: dprimalDelta = 0.0_DP

    ! Coefficient in front of the transposed convection operator \grad(.)^t y
    ! in the primal equation.
    real(dp) :: dprimalDeltaTrans = 0.0_DP
  
    ! Coefficient in front of the Newton operator \grad(y)(.) in the
    ! primal equation.
    real(dp) :: dprimalNewton = 0.0_DP
  
    ! Coefficient in front of the Newton operator \grad(y)^t(.) in the
    ! primal equation.
    real(dp) :: dprimalNewtonTrans = 0.0_DP

    ! Coefficient in front of the Mass matrix in the dual equation
    real(dp) :: ddualAlpha = 0.0_DP
  
    ! Coefficient in front of the Stokes operator nu*Laplace in the dual equation
    real(dp) :: ddualBeta = 0.0_DP
    
    ! Coefficient in front of the convection operator \grad(.) y in the
    ! dual equation.
    real(dp) :: ddualDelta = 0.0_DP

    ! Coefficient in front of the transposed convection operator \grad(.)^t y
    ! in the dual equation.
    real(dp) :: ddualDeltaTrans = 0.0_DP
  
    ! Coefficient in front of the Newton operator \grad(y)(.) in the
    ! dual equation.
    real(dp) :: ddualNewton = 0.0_DP
  
    ! Coefficient in front of the Newton operator \grad(y)^t(.) in the
    ! dual equation.
    real(dp) :: ddualNewtonTrans = 0.0_DP
    
    ! Coefficient in front of the transposed reactive convection operator
    ! \lambda\grad(.)^t in the dual equation
    real(dp) :: ddualRDeltaTrans = 0.0_DP
    
    ! Coefficient in front of the reactive Newton operator
    ! \grad(\lambda)(.) in the dual equation
    real(dp) :: ddualRNewton = 0.0_DP
    
    ! Coefficients in front of the reactive mass matrix that couples the dual
    ! to the primal equation
    real(dp) :: ddualRAlpha = 0.0_DP
    
    ! Multiplier to convert the dual velocity to the control.
    ! This is usually = -1/alpha for alpha coming from the optimal control problem.
    real(dp) :: dcontrolMultiplier = 0.0_DP
    
    ! Weight in front of the control u (= -1/alpha lambda).
    real(dp) :: dcontrolWeight = 0.0_DP
    
    ! Type of projection to use when converting the dual velocity to the
    ! control u.
    ! =0: No projection
    ! =1: DOF-based projection
    ! =2: Cubature-point based projection
    integer :: ccontrolProjection = 0
    
    ! Minimum/Maximum bound for the X-control.
    real(dp) :: dmin1 = -1E99_DP
    real(dp) :: dmax1 = 1E99_DP

    ! Minimum/Maximum bound for the Y-control.
    real(dp) :: dmin2 = -1E99_DP
    real(dp) :: dmax2 = 1E99_DP

    ! Calculation of local H.
    ! In 3D, there are 2 different methods to calculate the local H which
    ! is needed for the assembly of the convective part.
    ! =0: Use the cube-root of the volume of the hexahedron as local H
    ! =1: Use the length of the way that a particle travels through
    !     the hexahedron in direction of the flow
    integer :: clocalH = 1

  end type

  ! Assembly structure that saves information about the current
  ! status during the assembly of the local/global matrices.

  type t_optcassemblyinfo
  
    ! Block discretisation structure that defines the discretisation
    type(t_spatialDiscretisation), pointer :: p_rdiscretisation
    
    ! Id of the element distribution to be processed.
    integer :: ielemDistribution
    
    ! Pointer to the current element distribution
    type(t_elementDistribution), pointer :: p_relementDistribution

    ! Derivative specifier for evaluating elements
    logical, dimension(EL_MAXNDER) :: Bder
    
    ! An allocateable array accepting the DOF's of a set of elements.
    integer(PREC_DOFIDX), dimension(:,:), pointer :: Idofs
    
    ! Allocateable arrays for the values of the basis functions - 
    ! for test and trial spaces.
    real(DP), dimension(:,:,:,:), pointer :: Dbas

    ! An array receiving the coordinates of cubature points on
    ! the reference element for all elements in a set.
    real(DP), dimension(:,:), pointer :: p_DcubPtsRef

    ! For every cubature point on the reference element,
    ! the corresponding cubature weight
    real(DP), dimension(CUB_MAXCUBP) :: Domega

    ! Type of transformation from the reference to the real element 
    integer :: ctrafoType

    ! An element evaluation set for evaluating elements.
    type(t_evalElementSet) :: revalElementSet
    
    ! Specifier to signal when cubature points are initialised.
    logical :: bcubPtsInitialised
    
    ! Number of local DOF's
    integer :: indof
    
    ! Number of cubature points on the reference element
    integer :: ncubp
    
    ! Maximum number of simultaneously processed elements per element set
    integer :: nelementsPerBlock
    
    ! Local matrices, used during the assembly.
    ! Values and positions of values in the global matrix.
    integer(PREC_DOFIDX), dimension(:,:,:), pointer :: Kentry
    
    ! Additional contributions for the submatrices A11, A12, A21, A22 stemming from Newton.
    integer(PREC_DOFIDX), dimension(:,:,:), pointer :: Kentry12
    
    ! Maximum velocity magnitude
    real(DP) :: dumax
    
    ! An array with local DELTA's, each DELTA for one element
    real(DP), dimension(:), pointer :: DlocalDelta
    
    ! A pointer to an element-number list
    integer(I32), dimension(:), pointer :: p_IelementList

    ! Pointer to the primal/dual velocity field in the cubature points.
    real(DP), dimension(:,:,:), pointer :: Dpvel,Ddvel,Dcontrol
    
    ! Pointer to the velocity X- and Y-derivative of the primal velocity
    ! in the cubature points
    real(DP), dimension(:,:,:), pointer :: DpvelXderiv,DpvelYderiv
    real(DP), dimension(:,:,:), pointer :: DdvelXderiv,DdvelYderiv
    
    ! Pointer to the primal/dual velocity DOF's on the elements
    real(DP), dimension(:,:,:), pointer :: DpvelDofs,DdvelDofs
    
  end type

contains

! ***************************************************************************

  !<subroutine>

  subroutine computeLocalMatricesDiag (Dbas,Domega,Ddetj,ndof,ncubp,NEL,&
      Dvelocity,DvelocityXderiv,DvelocityYderiv,DlocalDelta,&
      DentryA11,DentryA22,DentryA44,DentryA55,&
      computeForm,dweight,dnu)
      
  !<description>
    ! Computes a local matrix to be incorporated into the global matrix.
    ! Variant for handling the diagonal matrices.
  !</description>
      
  !<input>
      
    ! Basis functions in all cubature points
    real(DP), dimension(:,:,:,:), intent(in)  :: Dbas
    
    ! Cubature weights
    real(DP), dimension(:), intent(in)        :: Domega
    
    ! Jacobian determinants in all cubature points
    real(DP), dimension(:,:), intent(in)      :: Ddetj
    
    ! Number of local DOF's in each element
    integer, intent(in)                       :: ndof
    
    ! Number of cubature points on each element
    integer, intent(in)                       :: ncubp
    
    ! Number of elements
    integer, intent(in)                       :: NEL
    
    ! Velocity vector and its derivatives in the cubature points
    real(DP), dimension(:,:,:), intent(in)    :: Dvelocity
    real(DP), dimension(:,:,:), intent(in)    :: DvelocityXderiv
    real(DP), dimension(:,:,:), intent(in)    :: DvelocityYderiv
    
    ! Local delta of the SD method
    real(DP), dimension(:), intent(in)        :: DlocalDelta
    
    ! Weight for A11-22, A44-55, A41-52, A14-25
    real(DP), dimension(2,2), intent(in) :: Dweight
    
    ! Viscosity parameter
    real(DP), intent(in) :: dnu

    interface
      subroutine computeForm (dbasI,dbasIX,dbasIY,dbasJ,dbasJX,dbasJY,&
          du1,du1x,du1y,du2,du2x,du2y,dnu,dlocalDelta,Dweight,da11,da22,da44,da55)
          
      use fsystem
          
      !<description>  
        ! Compute the form in a cubature point.
      !</description>
        
      !<input>
        ! Basis function phi_i and its derivatives.
        real(DP), intent(in) :: dbasI,dbasIX,dbasIY
        
        ! Basis function phi_j and its derivatives.
        real(DP), intent(in) :: dbasJ,dbasJX,dbasJY
        
        ! X-velocity of the given solution vector and its derivatives.
        real(DP), intent(in) :: du1,du1x,du1y
        
        ! Y-velocity of the given solution vector and its derivatives.
        real(DP), intent(in) :: du2,du2x,du2y
        
        ! Viscosity parameter
        real(DP), intent(in) :: dnu
        
        ! Local delta 
        real(DP), intent(in) :: dlocalDelta
        
        ! Weight for A11-22, A44-55, A41-52, A14-25
        real(DP), dimension(2,2), intent(in) :: Dweight
      !</input>
        
      !<output>
        ! Values of the form
        real(DP), intent(out) :: da11,da22,da44,da55
      !</output>
      
      end subroutine
    end interface   
    
    external :: computeForm
    
  !</input>
  
  !<output>
    
    ! Entries in the matrix
    real(DP), dimension(:,:,:), intent(inout) :: DentryA11
    real(DP), dimension(:,:,:), intent(inout) :: DentryA22
    real(DP), dimension(:,:,:), intent(inout) :: DentryA44
    real(DP), dimension(:,:,:), intent(inout) :: DentryA55
  
  !</output>
    
  !</subroutine>
    
    ! local variables
    integer :: iel, icubp, idofe, jdofe
    real(DP) :: OM,du1,du2,du1locx,du1locy,du2locx,du2locY
    real(DP) :: HBASI1,HBASI2,HBASI3,HBASJ1,HBASJ2,HBASJ3
    real(DP) :: AH11,AH22,AH44,AH55
    
    AH11 = 0.0_DP
    AH22 = 0.0_DP
    AH44 = 0.0_DP
    AH55 = 0.0_DP
    
    ! Loop over the elements in the current set.
    do iel=1,NEL

      ! Loop over all cubature points on the current element
      do icubp = 1, ncubp

        ! Calculate the current weighting factor in the cubature formula
        ! in that cubature point.
        !
        ! Normally, we have to take the absolut value of the determinant 
        ! of the mapping here!
        ! In 2D, the determinant is always positive, whereas in 3D,
        ! the determinant might be negative -- that's normal!
        ! But because this routine only works in 2D, we can skip
        ! the ABS here!

        OM = Domega(icubp)*Ddetj(icubp,iel)

        ! Current velocity in this cubature point:
        du1 = Dvelocity (1,icubp,iel)
        du2 = Dvelocity (2,icubp,iel)
        du1locx = DvelocityXderiv (1,icubp,iel)
        du1locy = DvelocityXderiv (2,icubp,iel)
        du2locx = DvelocityYderiv (1,icubp,iel)
        du2locy = DvelocityYderiv (2,icubp,iel)
        
        ! Outer loop over the DOF's i=1..indof on our current element, 
        ! which corresponds to the basis functions Phi_i:

        do idofe=1,ndof
        
          ! Fetch the contributions of the (test) basis functions Phi_i
          ! (our "O")  for function value and first derivatives for the 
          ! current DOF into HBASIy:
        
          HBASI1 = Dbas(idofe,1,icubp,iel)
          HBASI2 = Dbas(idofe,2,icubp,iel)
          HBASI3 = Dbas(idofe,3,icubp,iel)
          
          ! Inner loop over the DOF's j=1..indof, which corresponds to
          ! the basis function Phi_j:

          do jdofe=1,ndof
            
            ! Fetch the contributions of the (trial) basis function Phi_j
            ! (out "X") for function value and first derivatives for the 
            ! current DOF into HBASJy:
          
            HBASJ1 = Dbas(jdofe,1,icubp,iel)
            HBASJ2 = Dbas(jdofe,2,icubp,iel)
            HBASJ3 = Dbas(jdofe,3,icubp,iel)

            ! Finally calculate the contribution to the system
            ! matrices A11, A12, A21 and A22.
            call computeForm (HBASI1,HBASI2,HBASI3,HBASJ1,HBASJ2,HBASJ3,&
                du1,du2,du1locx,du1locy,du2locx,du2locy,dnu,DlocalDelta(iel),&
                Dweight,AH11,AH22,AH44,AH55)
                
            ! Weighten the calculated value AHxy by the cubature
            ! weight OM and add it to the local matrices. After the
            ! loop over all DOF's is finished, each entry contains
            ! the calculated integral.

            DentryA11(jdofe,idofe,iel) = DentryA11(jdofe,idofe,iel)+OM*AH11
            DentryA22(jdofe,idofe,iel) = DentryA22(jdofe,idofe,iel)+OM*AH22
            DentryA44(jdofe,idofe,iel) = DentryA44(jdofe,idofe,iel)+OM*AH44
            DentryA55(jdofe,idofe,iel) = DentryA55(jdofe,idofe,iel)+OM*AH55
            
          end do ! idofe
          
        end do ! jdofe

      end do ! icubp 
    
    end do ! iel
    
  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine computeLocalMatricesFullDiag (Dbas,Domega,Ddetj,ndof,ncubp,NEL,&
      Dvelocity,DvelocityXderiv,DvelocityYderiv,DlocalDelta,&
      DentryA11,DentryA22,DentryA44,DentryA55,&
      DentryA12,DentryA21,DentryA45,DentryA54,&
      computeForm,dweight,dnu)
      
  !<description>
    ! Computes a local matrix to be incorporated into the global matrix.
    ! Variant for handling the diagonal matrix blocks.
  !</description>
      
  !<input>
      
    ! Basis functions in all cubature points
    real(DP), dimension(:,:,:,:), intent(in)  :: Dbas
    
    ! Cubature weights
    real(DP), dimension(:), intent(in)        :: Domega
    
    ! Jacobian determinants in all cubature points
    real(DP), dimension(:,:), intent(in)      :: Ddetj
    
    ! Number of local DOF's in each element
    integer, intent(in)                       :: ndof
    
    ! Number of cubature points on each element
    integer, intent(in)                       :: ncubp
    
    ! Number of elements
    integer, intent(in)                       :: NEL
    
    ! Velocity vector and its derivatives in the cubature points
    real(DP), dimension(:,:,:), intent(in)    :: Dvelocity
    real(DP), dimension(:,:,:), intent(in)    :: DvelocityXderiv
    real(DP), dimension(:,:,:), intent(in)    :: DvelocityYderiv
    
    ! Local delta of the SD method
    real(DP), dimension(:), intent(in)        :: DlocalDelta

    ! Weight for A11-22, A44-55, A41-52, A14-25
    real(DP), dimension(2,2), intent(in) :: Dweight
    
    ! Viscosity parameter
    real(DP), intent(in) :: dnu

    interface
      subroutine computeForm (dbasI,dbasIX,dbasIY,dbasJ,dbasJX,dbasJY,&
          du1,du1x,du1y,du2,du2x,du2y,dnu,dlocalDelta,Dweight,&
          da11,da22,da44,da55,da12,da21,da45,da54)
          
      use fsystem
          
      !<description>  
        ! Compute the form in a cubature point.
      !</description>
        
      !<input>
        ! Basis function phi_i and its derivatives.
        real(DP), intent(in) :: dbasI,dbasIX,dbasIY
        
        ! Basis function phi_j and its derivatives.
        real(DP), intent(in) :: dbasJ,dbasJX,dbasJY
        
        ! X-velocity of the given solution vector and its derivatives.
        real(DP), intent(in) :: du1,du1x,du1y
        
        ! Y-velocity of the given solution vector and its derivatives.
        real(DP), intent(in) :: du2,du2x,du2y
        
        ! Viscosity parameter
        real(DP), intent(in) :: dnu
        
        ! Local delta 
        real(DP), intent(in) :: dlocalDelta

        ! Weight for A11-22, A44-55, A41-52, A14-25
        real(DP), dimension(2,2), intent(in) :: Dweight
      !</input>
        
      !<output>
        ! Values of the form
        real(DP), intent(out) :: da11,da22,da44,da55,da12,da21,da45,da54
      !</output>
      
      end subroutine
    end interface   
    
    external :: computeForm

  !</input>    
    
  !<output>
    
    ! Entries in the matrix
    real(DP), dimension(:,:,:), intent(inout) :: DentryA11
    real(DP), dimension(:,:,:), intent(inout) :: DentryA22
    real(DP), dimension(:,:,:), intent(inout) :: DentryA44
    real(DP), dimension(:,:,:), intent(inout) :: DentryA55
    real(DP), dimension(:,:,:), intent(inout) :: DentryA12
    real(DP), dimension(:,:,:), intent(inout) :: DentryA21
    real(DP), dimension(:,:,:), intent(inout) :: DentryA45
    real(DP), dimension(:,:,:), intent(inout) :: DentryA54
    
  !</output>
    
  !</subroutine>
    
    ! local variables
    integer :: iel, icubp, idofe, jdofe
    real(DP) :: OM,du1,du2,du1locx,du1locy,du2locx,du2locY
    real(DP) :: HBASI1,HBASI2,HBASI3,HBASJ1,HBASJ2,HBASJ3
    real(DP) :: AH11,AH22,AH44,AH55,AH12,AH21,AH45,AH54
    
    AH11 = 0.0_DP
    AH22 = 0.0_DP
    AH44 = 0.0_DP
    AH55 = 0.0_DP
    AH12 = 0.0_DP
    AH21 = 0.0_DP
    AH45 = 0.0_DP
    AH54 = 0.0_DP
    
    ! Loop over the elements in the current set.
    do iel=1,NEL

      ! Loop over all cubature points on the current element
      do icubp = 1, ncubp

        ! Calculate the current weighting factor in the cubature formula
        ! in that cubature point.
        !
        ! Normally, we have to take the absolut value of the determinant 
        ! of the mapping here!
        ! In 2D, the determinant is always positive, whereas in 3D,
        ! the determinant might be negative -- that's normal!
        ! But because this routine only works in 2D, we can skip
        ! the ABS here!

        OM = Domega(icubp)*Ddetj(icubp,iel)

        ! Current velocity in this cubature point:
        du1 = Dvelocity (1,icubp,iel)
        du2 = Dvelocity (2,icubp,iel)
        du1locx = DvelocityXderiv (1,icubp,iel)
        du1locy = DvelocityXderiv (2,icubp,iel)
        du2locx = DvelocityYderiv (1,icubp,iel)
        du2locy = DvelocityYderiv (2,icubp,iel)
        
        ! Outer loop over the DOF's i=1..indof on our current element, 
        ! which corresponds to the basis functions Phi_i:

        do idofe=1,ndof
        
          ! Fetch the contributions of the (test) basis functions Phi_i
          ! (our "O")  for function value and first derivatives for the 
          ! current DOF into HBASIy:
        
          HBASI1 = Dbas(idofe,1,icubp,iel)
          HBASI2 = Dbas(idofe,2,icubp,iel)
          HBASI3 = Dbas(idofe,3,icubp,iel)
          
          ! Inner loop over the DOF's j=1..indof, which corresponds to
          ! the basis function Phi_j:

          do jdofe=1,ndof
            
            ! Fetch the contributions of the (trial) basis function Phi_j
            ! (out "X") for function value and first derivatives for the 
            ! current DOF into HBASJy:
          
            HBASJ1 = Dbas(jdofe,1,icubp,iel)
            HBASJ2 = Dbas(jdofe,2,icubp,iel)
            HBASJ3 = Dbas(jdofe,3,icubp,iel)

            ! Finally calculate the contribution to the system matrices.
            call computeForm (HBASI1,HBASI2,HBASI3,HBASJ1,HBASJ2,HBASJ3,&
                du1,du2,du1locx,du1locy,du2locx,du2locy,dnu,DlocalDelta(iel),&
                Dweight,AH11,AH22,AH44,AH55,AH12,AH21,AH45,AH54)
                
            ! Weighten the calculated value AHxy by the cubature
            ! weight OM and add it to the local matrices. After the
            ! loop over all DOF's is finished, each entry contains
            ! the calculated integral.

            DentryA11(jdofe,idofe,iel) = DentryA11(jdofe,idofe,iel)+OM*AH11
            DentryA22(jdofe,idofe,iel) = DentryA22(jdofe,idofe,iel)+OM*AH22
            DentryA44(jdofe,idofe,iel) = DentryA44(jdofe,idofe,iel)+OM*AH44
            DentryA55(jdofe,idofe,iel) = DentryA55(jdofe,idofe,iel)+OM*AH55
            
            DentryA12(jdofe,idofe,iel) = DentryA12(jdofe,idofe,iel)+OM*AH12
            DentryA21(jdofe,idofe,iel) = DentryA21(jdofe,idofe,iel)+OM*AH21
            DentryA45(jdofe,idofe,iel) = DentryA45(jdofe,idofe,iel)+OM*AH45
            DentryA54(jdofe,idofe,iel) = DentryA54(jdofe,idofe,iel)+OM*AH54
            
          end do ! idofe
          
        end do ! jdofe

      end do ! icubp 
    
    end do ! iel
    
  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine computeLocalMatricesFull (Dbas,Domega,Ddetj,ndof,ncubp,NEL,&
      Dvelocity,DvelocityXderiv,DvelocityYderiv,DlocalDelta,&
      DentryA11,DentryA22,DentryA44,DentryA55,&
      DentryA12,DentryA21,DentryA45,DentryA54,&
      DentryA41,DentryA52,DentryA42,DentryA51,&
      DentryA14,DentryA25,DentryA24,DentryA15,&
      computeForm,dweight,dnu)
      
  !<description>
    ! Computes a local matrix to be incorporated into the global matrix.
    ! Variant for handling the full matrix.
  !</description>
      
  !<input>      
  
    ! Basis functions in all cubature points
    real(DP), dimension(:,:,:,:), intent(in)  :: Dbas
    
    ! Cubature weights
    real(DP), dimension(:), intent(in)        :: Domega
    
    ! Jacobian determinants in all cubature points
    real(DP), dimension(:,:), intent(in)      :: Ddetj
    
    ! Number of local DOF's in each element
    integer, intent(in)                       :: ndof
    
    ! Number of cubature points on each element
    integer, intent(in)                       :: ncubp
    
    ! Number of elements
    integer, intent(in)                       :: NEL
    
    ! Velocity vector and its derivatives in the cubature points
    real(DP), dimension(:,:,:), intent(in)    :: Dvelocity
    real(DP), dimension(:,:,:), intent(in)    :: DvelocityXderiv
    real(DP), dimension(:,:,:), intent(in)    :: DvelocityYderiv
    
    ! Local delta of the SD method
    real(DP), dimension(:), intent(in)        :: DlocalDelta
    
    ! Weight for A11-22, A44-55, A41-52, A14-25
    real(DP), dimension(2,2), intent(in) :: Dweight
    
    ! Viscosity parameter
    real(DP), intent(in) :: dnu

    interface
      subroutine computeForm (dbasI,dbasIX,dbasIY,dbasJ,dbasJX,dbasJY,&
          du1,du1x,du1y,du2,du2x,du2y,dnu,dlocalDelta,Dweight,&
          da11,da22,da44,da55,da12,da21,da45,da54,&
          da41,da52,da42,da51,da14,da25,da24,da15)
          
      use fsystem
          
      !<description>  
        ! Compute the form in a cubature point.
        ! Variant for handling the full matrix.
      !</description>
        
      !<input>
        ! Basis function phi_i and its derivatives.
        real(DP), intent(in) :: dbasI,dbasIX,dbasIY
        
        ! Basis function phi_j and its derivatives.
        real(DP), intent(in) :: dbasJ,dbasJX,dbasJY
        
        ! X-velocity of the given solution vector and its derivatives.
        real(DP), intent(in) :: du1,du1x,du1y
        
        ! Y-velocity of the given solution vector and its derivatives.
        real(DP), intent(in) :: du2,du2x,du2y
        
        ! Viscosity parameter
        real(DP), intent(in) :: dnu
        
        ! Local delta 
        real(DP), intent(in) :: dlocalDelta

        ! Weight for A11-22, A44-55, A41-52, A14-25
        real(DP), dimension(2,2), intent(in) :: Dweight
      !</input>
        
      !<output>
        ! Values of the form
        real(DP), intent(out) :: da11,da22,da44,da55,da12,da21,da45,da54
        real(DP), intent(out) :: da41,da52,da42,da51,da14,da25,da24,da15
      !</output>
      
      end subroutine
    end interface   
    
    external :: computeForm
    
  !</input>
    
  !<output>
    
    ! Entries in the matrix
    real(DP), dimension(:,:,:), intent(inout) :: DentryA11
    real(DP), dimension(:,:,:), intent(inout) :: DentryA22
    real(DP), dimension(:,:,:), intent(inout) :: DentryA44
    real(DP), dimension(:,:,:), intent(inout) :: DentryA55
    real(DP), dimension(:,:,:), intent(inout) :: DentryA12
    real(DP), dimension(:,:,:), intent(inout) :: DentryA21
    real(DP), dimension(:,:,:), intent(inout) :: DentryA45
    real(DP), dimension(:,:,:), intent(inout) :: DentryA54
    
    real(DP), dimension(:,:,:), intent(inout) :: DentryA41
    real(DP), dimension(:,:,:), intent(inout) :: DentryA52
    real(DP), dimension(:,:,:), intent(inout) :: DentryA42
    real(DP), dimension(:,:,:), intent(inout) :: DentryA51
    real(DP), dimension(:,:,:), intent(inout) :: DentryA14
    real(DP), dimension(:,:,:), intent(inout) :: DentryA25
    real(DP), dimension(:,:,:), intent(inout) :: DentryA24
    real(DP), dimension(:,:,:), intent(inout) :: DentryA15
 
  !</output>   
  
  !<subroutine>
    
    ! local variables
    integer :: iel, icubp, idofe, jdofe
    real(DP) :: OM,du1,du2,du1locx,du1locy,du2locx,du2locY
    real(DP) :: HBASI1,HBASI2,HBASI3,HBASJ1,HBASJ2,HBASJ3
    real(DP) :: AH11,AH22,AH44,AH55,AH12,AH21,AH45,AH54
    real(DP) :: AH41,AH52,AH42,AH51,AH14,AH25,AH24,AH15
    
    AH11 = 0.0_DP
    AH22 = 0.0_DP
    AH44 = 0.0_DP
    AH55 = 0.0_DP
    AH12 = 0.0_DP
    AH21 = 0.0_DP
    AH45 = 0.0_DP
    AH54 = 0.0_DP
    
    AH41 = 0.0_DP
    AH52 = 0.0_DP
    AH42 = 0.0_DP
    AH51 = 0.0_DP
    AH14 = 0.0_DP
    AH25 = 0.0_DP
    AH24 = 0.0_DP
    AH15 = 0.0_DP
    
    ! Loop over the elements in the current set.
    do iel=1,NEL

      ! Loop over all cubature points on the current element
      do icubp = 1, ncubp

        ! Calculate the current weighting factor in the cubature formula
        ! in that cubature point.
        !
        ! Normally, we have to take the absolut value of the determinant 
        ! of the mapping here!
        ! In 2D, the determinant is always positive, whereas in 3D,
        ! the determinant might be negative -- that's normal!
        ! But because this routine only works in 2D, we can skip
        ! the ABS here!

        OM = Domega(icubp)*Ddetj(icubp,iel)

        ! Current velocity in this cubature point:
        du1 = Dvelocity (1,icubp,iel)
        du2 = Dvelocity (2,icubp,iel)
        du1locx = DvelocityXderiv (1,icubp,iel)
        du1locy = DvelocityXderiv (2,icubp,iel)
        du2locx = DvelocityYderiv (1,icubp,iel)
        du2locy = DvelocityYderiv (2,icubp,iel)
        
        ! Outer loop over the DOF's i=1..indof on our current element, 
        ! which corresponds to the basis functions Phi_i:

        do idofe=1,ndof
        
          ! Fetch the contributions of the (test) basis functions Phi_i
          ! (our "O")  for function value and first derivatives for the 
          ! current DOF into HBASIy:
        
          HBASI1 = Dbas(idofe,1,icubp,iel)
          HBASI2 = Dbas(idofe,2,icubp,iel)
          HBASI3 = Dbas(idofe,3,icubp,iel)
          
          ! Inner loop over the DOF's j=1..indof, which corresponds to
          ! the basis function Phi_j:

          do jdofe=1,ndof
            
            ! Fetch the contributions of the (trial) basis function Phi_j
            ! (out "X") for function value and first derivatives for the 
            ! current DOF into HBASJy:
          
            HBASJ1 = Dbas(jdofe,1,icubp,iel)
            HBASJ2 = Dbas(jdofe,2,icubp,iel)
            HBASJ3 = Dbas(jdofe,3,icubp,iel)

            ! Finally calculate the contribution to the system matrices.
            call computeForm (HBASI1,HBASI2,HBASI3,HBASJ1,HBASJ2,HBASJ3,&
                du1,du2,du1locx,du1locy,du2locx,du2locy,dnu,DlocalDelta(iel),&
                Dweight,AH11,AH22,AH44,AH55,AH12,AH21,AH45,AH54,&
                AH41,AH52,AH42,AH51,AH14,AH25,AH24,AH15)
                
            ! Weighten the calculated value AHxy by the cubature
            ! weight OM and add it to the local matrices. After the
            ! loop over all DOF's is finished, each entry contains
            ! the calculated integral.

            DentryA11(jdofe,idofe,iel) = DentryA11(jdofe,idofe,iel)+OM*AH11
            DentryA22(jdofe,idofe,iel) = DentryA22(jdofe,idofe,iel)+OM*AH22
            DentryA44(jdofe,idofe,iel) = DentryA44(jdofe,idofe,iel)+OM*AH44
            DentryA55(jdofe,idofe,iel) = DentryA55(jdofe,idofe,iel)+OM*AH55
            
            DentryA12(jdofe,idofe,iel) = DentryA12(jdofe,idofe,iel)+OM*AH12
            DentryA21(jdofe,idofe,iel) = DentryA21(jdofe,idofe,iel)+OM*AH21
            DentryA45(jdofe,idofe,iel) = DentryA45(jdofe,idofe,iel)+OM*AH45
            DentryA54(jdofe,idofe,iel) = DentryA54(jdofe,idofe,iel)+OM*AH54

            DentryA14(jdofe,idofe,iel) = DentryA14(jdofe,idofe,iel)+OM*AH14
            DentryA25(jdofe,idofe,iel) = DentryA25(jdofe,idofe,iel)+OM*AH25
            DentryA24(jdofe,idofe,iel) = DentryA24(jdofe,idofe,iel)+OM*AH24
            DentryA15(jdofe,idofe,iel) = DentryA15(jdofe,idofe,iel)+OM*AH15
                                                                           
            DentryA41(jdofe,idofe,iel) = DentryA41(jdofe,idofe,iel)+OM*AH41
            DentryA52(jdofe,idofe,iel) = DentryA52(jdofe,idofe,iel)+OM*AH52
            DentryA51(jdofe,idofe,iel) = DentryA51(jdofe,idofe,iel)+OM*AH51
            DentryA42(jdofe,idofe,iel) = DentryA42(jdofe,idofe,iel)+OM*AH42
            
          end do ! idofe
          
        end do ! jdofe

      end do ! icubp 
    
    end do ! iel
    
  end subroutine

  ! ***************************************************************************

!<subroutine>
  subroutine conv_strdiffOptC2dblk1 ( &
                  rvelocityVector,rmatrix,&
                  dupsam,dnu,Dalpha,Dbeta,Dtheta,Ddelta,Dnewton, &
                  DdeltaTransposed,DnewtonTransposed)
!<description>

!</description>

!<input>

  ! Velocity vector for the nonlinearity.
  ! The first blocks 1/2 and 4/5 in this vector define the evaluation
  ! point (primal/dual velocity).
  type(t_vectorBlock), intent(in) :: rvelocityVector
  
  ! dupsam  - control parameter.
  !          -1: simple upwind,
  !          =0: Samarskji upwind
  real(DP), intent(IN) :: dupsam
  
  ! Viscosity parameter $\nu = 1/Re$ if viscosity is constant
  real(DP), intent(IN) :: dnu 
  
  ! Weighting factor for the mass matrix.
  real(DP), dimension(2,2), intent(IN) :: Dalpha

  ! Weighting factor for the Stokes matrix. (Stokes matrix = 1/Re * Laplace)
  real(DP), dimension(2,2), intent(IN) :: Dbeta

  ! Weighting factor of the convective operator: $\theta * u*grad(u)$. 
  ! For time-dependent problems, this can be set to the step size
  ! in the $\Theta$-scheme.
  real(DP), dimension(2,2), intent(IN) :: Dtheta 
  
  ! Weighting factor for the nonlinear term
  real(DP), dimension(2,2), intent(IN) :: Ddelta

  ! Weighting factor of the Newton matrix. A value of 0.0 deactivates the
  ! Newton part. A value != 0.0 activates Newton; in this case the submatrices
  ! A12 and A21 must be present in rmatrix.
  real(DP), dimension(2,2), intent(IN) :: Dnewton

  ! Weighting factor of the transposed convection matrix. A value of 0.0 deactivates
  ! this operator.
  real(DP), dimension(2,2), intent(IN) :: DdeltaTransposed
  
  ! Weighting factor of the transposed Newton matrix. A value of 0.0 deactivates
  ! this operator.
  real(DP), dimension(2,2), intent(IN) :: DnewtonTransposed
      
!</input>

!<inputoutput>
  ! The system matrix. The submatrices for the velocity must be in block
  ! A11, A12, A21 and A22 and must be in matrix format 7 or 9.
  ! A11 and A22 must have the same structure. A12 and A21 must have
  ! the same structure.
  type(t_matrixBlock), intent(INOUT), target :: rmatrix
!</inputoutput>

!</subroutine>

  ! local variables
  integer :: indof,IEQ,I,K,idofe,jdofe,icubp
  integer(PREC_DOFIDX) :: jcol0,jdfg,jcol
  integer(PREC_ELEMENTIDX) :: iel,ielset,ielmax
  logical, dimension(EL_MAXNDER) :: Bder
  real(DP) :: dumax,dumaxr, du1loc, du2loc, dunorm,db,dre
  real(DP) :: du1locx,du2locx,du1locy,du2locy,dbx,dby
  integer :: NVE
  real(dp) :: dalphamax,dbetamax,dthetamax,ddeltamax,dnewtonmax
  real(dp) :: ddeltaTransposedMax,dnewtonTransposedMax
  
  ! Matrix structure arrays
  integer(PREC_VECIDX), dimension(:), pointer :: p_Kcol
  integer(PREC_MATIDX), dimension(:), pointer :: p_Kld
  real(DP), dimension(:), pointer :: p_Da11,p_Da22,p_Da44,p_Da55

  integer(PREC_VECIDX), dimension(:), pointer :: p_Kcol12
  integer(PREC_MATIDX), dimension(:), pointer :: p_Kld12
  real(DP), dimension(:), pointer :: p_Da12,p_Da21,p_Da45,p_Da54
  
  real(DP), dimension(:), pointer :: p_DpvelocityX, p_DpvelocityY
  
  ! An array receiving the coordinates of cubature points on
  ! the reference element for all elements in a set.
  real(DP), dimension(:,:), pointer :: p_DcubPtsRef

  ! The discretisation - for easier access
  type(t_spatialDiscretisation), pointer :: p_rdiscretisation
  
  ! Triangulation
  type(t_triangulation), pointer :: p_rtriangulation
  real(DP), dimension(:,:), pointer :: p_DvertexCoords
  integer(PREC_VERTEXIDX), dimension(:,:), pointer :: p_IedgesAtElement,p_IverticesAtElement

  ! Number of elements in a block. Normally =BILF_NELEMSIM,
  ! except if there are less elements in the discretisation.
  integer :: nelementsPerBlock

  ! One and only element distribution
  type(t_elementDistribution), pointer :: p_relementDistribution

  ! Cubature point coordinates on the reference element
  real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi

  ! For every cubature point on the reference element,
  ! the corresponding cubature weight
  real(DP), dimension(CUB_MAXCUBP) :: Domega
  
  ! number of cubature points on the reference element
  integer :: ncubp

  ! An element evaluation set for evaluating elements.
  type(t_evalElementSet) :: revalElementSet
  logical :: bcubPtsInitialised

  ! Arrays for saving Jacobian determinants and matrices
  real(DP), dimension(:,:), pointer :: p_Ddetj
  
  ! An allocateable array accepting the DOF's of a set of elements.
  integer(PREC_DOFIDX), dimension(:,:), allocatable, target :: Idofs
  
  ! Allocateable arrays for the values of the basis functions - 
  ! for test and trial spaces.
  real(DP), dimension(:,:,:,:), allocatable, target :: Dbas

  ! Local matrices, used during the assembly.
  ! Values and positions of values in the global matrix.
  integer(PREC_DOFIDX), dimension(:,:,:), allocatable :: Kentry
  
  ! Additional contributions for the submatrices A11, A12, A21, A22 stemming from Newton.
  integer(PREC_DOFIDX), dimension(:,:,:), allocatable :: Kentry12
  real(DP), dimension(:,:,:), allocatable :: DentryA11
  real(DP), dimension(:,:,:), allocatable :: DentryA12
  real(DP), dimension(:,:,:), allocatable :: DentryA21
  real(DP), dimension(:,:,:), allocatable :: DentryA22
  
  real(DP), dimension(:,:,:), allocatable :: DentryA44
  real(DP), dimension(:,:,:), allocatable :: DentryA45
  real(DP), dimension(:,:,:), allocatable :: DentryA54
  real(DP), dimension(:,:,:), allocatable :: DentryA55
  
  real(DP), dimension(:,:,:), allocatable :: DentryA41
  real(DP), dimension(:,:,:), allocatable :: DentryA52
  real(DP), dimension(:,:,:), allocatable :: DentryA42
  real(DP), dimension(:,:,:), allocatable :: DentryA51
  
  real(DP), dimension(:,:,:), allocatable :: DentryA14
  real(DP), dimension(:,:,:), allocatable :: DentryA25
  real(DP), dimension(:,:,:), allocatable :: DentryA24
  real(DP), dimension(:,:,:), allocatable :: DentryA15
  
  ! A pointer to an element-number list
  integer(I32), dimension(:), pointer :: p_IelementList

  ! Pointer to the velocity field in the cubature points.
  real(DP), dimension(:,:,:), allocatable :: Dvelocity
  
  ! Pointer to the velocity X- and Y-derivative in the cubature points
  real(DP), dimension(:,:,:), allocatable :: DvelocityUderiv
  real(DP), dimension(:,:,:), allocatable :: DvelocityVderiv
  
  ! An array with local DELTA's, each DELTA for one element
  real(DP), dimension(:), allocatable :: DlocalDelta

  ! Type of transformation from the reference to the real element 
  integer :: ctrafoType
  
  ! Element evaluation tag; collects some information necessary for evaluating
  ! the elements.
  integer(I32) :: cevaluationTag
 
    ! Initialise the derivative flags
    Bder = .false.
    Bder(DER_FUNC) = .true.
    Bder(DER_DERIV_X) = .true.
    Bder(DER_DERIV_Y) = .true.
    
    ! Compute the maximum weighting factors. This allows us to
    ! figure out which operators are active and which not.
    dalphamax = maxval(abs(Dalpha))
    dbetamax = maxval(abs(Dbeta))
    dthetamax = maxval(abs(Dtheta))
    ddeltamax = maxval(abs(Ddelta))
    dnewtonmax = maxval(abs(Dnewton))
    ddeltaTransposedMax = maxval(abs(DdeltaTransposed))
    dnewtonTransposedMax = maxval(abs(DnewtonTransposed))

    ! Shortcut to the spatial discretisation.
    ! We assume the same for all, A11, A12, A21 and A22.
    p_rdiscretisation => rmatrix%RmatrixBlock(1,1)%p_rspatialDiscrTest
    
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
    
    ! Get the number of local DOF's for trial/test functions.
    ! We assume trial and test functions to be the same.
    indof = elem_igetNDofLoc(p_relementDistribution%celement)

    ! Number of local DOF's
    NVE = elem_igetNVE(p_relementDistribution%celement)
    
    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it's at most
    ! BILF_NELEMSIM. This is only used for allocating some arrays.
    nelementsPerBlock = min(BILF_NELEMSIM,p_rtriangulation%NEL)
    
    ! Get pointers to the matrix content (if necessary)
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,1),p_Da11)
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,2),p_Da22)
    
    ! Get pointers to the velocity to be evaluated.
    call lsyssc_getbase_double (rvelocityVector%RvectorBlock(1),p_DpvelocityX)
    call lsyssc_getbase_double (rvelocityVector%RvectorBlock(2),p_DpvelocityY)
    
    if (.not. lsysbl_isSubmatrixPresent(rmatrix,1,2)) then
      print *,"A12/A21 not present!"
      call sys_halt()
    end if
    
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,2),p_Da12)
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,1),p_Da21)
    
    ! Get pointers to the matrix structure(s).
    call lsyssc_getbase_Kcol (rmatrix%RmatrixBlock(1,1),p_Kcol)
    call lsyssc_getbase_Kld (rmatrix%RmatrixBlock(1,1),p_Kld)
    
    call lsyssc_getbase_Kcol (rmatrix%RmatrixBlock(1,2),p_Kcol12)
    call lsyssc_getbase_Kld (rmatrix%RmatrixBlock(1,2),p_Kld12)
   
    ! Get from the trial element space the type of coordinate system
    ! that is used there:
    ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)
    
    ! Allocate some memory to hold the cubature points on the reference element
    allocate(p_DcubPtsRef(trafo_igetReferenceDimension(ctrafoType),CUB_MAXCUBP))

    ! Initialise the cubature formula,
    ! Get cubature weights and point coordinates on the reference element
    call cub_getCubPoints(p_relementDistribution%ccubTypeBilForm, ncubp, Dxi, Domega)
    
    ! Reformat the cubature points; they are in the wrong shape!
    do i=1,ncubp
      do k=1,ubound(p_DcubPtsRef,1)
        p_DcubPtsRef(k,i) = Dxi(i,k)
      end do
    end do
    
    ! Open-MP-Extension: Open threads here.
    ! "csysTrial" is declared as private; shared gave errors with the Intel compiler
    ! in Windows!?!
    ! Each thread will allocate its own local memory...
        
    !%OMP PARALLEL private( &
    !%OMP p_Ddetj, i,k,Dbas,Idofs, &
    !%OMP DlocalDelta,Kentry,Kentry12,Dentry, &
    !%OMP DentryA11,DentryA12,DentryA21,DentryA22,Dvelocity, &
    !%OMP DvelocityUderiv,DvelocityVderiv,dre,iel,db,icubp,& 
    !%OMP idofe,jcol0,jdofe,jdfg,jcol,du1loc,du2loc,dbx,dby, &
    !%OMP du1locx,du1locy,du2locx,du2locy,OM,AH,HBASI1,HBASI2,& 
    !%OMP HBASI3,HBASJ1,HBASJ2,HBASJ3,HSUMI,HSUMJ,AH11,AH12,AH21, &
    !%OMP AH22,ielmax,revalElementSet)

    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  allocate(Dbas(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    ! would lead to nonused memory blocks in these arrays during the assembly, 
    ! which reduces the speed by 50%!
    
    allocate(Dbas(indof,elem_getMaxDerivative(p_relementDistribution%celement), &
             ncubp,nelementsPerBlock))

    ! Allocate memory for the DOF's of all the elements.
    allocate(Idofs(indof,nelementsPerBlock))
    
    ! Allocate memory for array with local DELTA's
    allocate(DlocalDelta(nelementsPerBlock))

    ! Allocate an array saving the local matrices for all elements
    ! in an element set.
    ! We could also allocate EL_MAXNBAS*EL_MAXNBAS*BILF_NELEMSIM integers
    ! for this local matrix, but this would normally not fit to the cache
    ! anymore! indofTrial*indofTest*BILF_NELEMSIM is normally much smaller!
    !
    ! Kentry (:,:,:) defines the positions of the local matrices
    ! in the submatrices A11 and A22.
    ! KentryA12 (:,:,:) defines the positions of the local matrices
    ! in the submatrices A12 and A21.
    allocate(Kentry(indof,indof,nelementsPerBlock))
    
    allocate(Kentry12(indof,indof,nelementsPerBlock))
    
    ! Dentry (:,:,:) fetches the 'main' matrix entries (Laplace, Mass,
    ! Convection).
    ! DentryA11, DentryA12, DentryA21 and DentryA22 fetches additional entries in 
    ! A11, A12, A21 and A22 of the Newton matrix, which is not always calculated
    ! and therefore not always used!
    allocate(DentryA11(indof,indof,nelementsPerBlock))
    allocate(DentryA12(indof,indof,nelementsPerBlock))
    allocate(DentryA21(indof,indof,nelementsPerBlock))
    allocate(DentryA22(indof,indof,nelementsPerBlock))

    allocate(DentryA44(indof,indof,nelementsPerBlock))
    allocate(DentryA45(indof,indof,nelementsPerBlock))
    allocate(DentryA54(indof,indof,nelementsPerBlock))
    allocate(DentryA55(indof,indof,nelementsPerBlock))

    allocate(DentryA41(indof,indof,nelementsPerBlock))
    allocate(DentryA52(indof,indof,nelementsPerBlock))
    allocate(DentryA42(indof,indof,nelementsPerBlock))
    allocate(DentryA51(indof,indof,nelementsPerBlock))
             
    allocate(DentryA14(indof,indof,nelementsPerBlock))
    allocate(DentryA25(indof,indof,nelementsPerBlock))
    allocate(DentryA24(indof,indof,nelementsPerBlock))
    allocate(DentryA15(indof,indof,nelementsPerBlock))
    
    ! Allocate memory for the velocity in the cubature points.
    allocate(Dvelocity(NDIM2D,ncubp,nelementsPerBlock))
    
    ! Initialisation of the element set.
    call elprep_init(revalElementSet)

    ! Indicate that cubature points must still be initialised in the element set.
    bcubPtsInitialised = .false.
    
    allocate(DvelocityUderiv(NDIM2D,ncubp,nelementsPerBlock))
    allocate(DvelocityVderiv(NDIM2D,ncubp,nelementsPerBlock))
    
    ! What is the reciprocal of nu? We need it later.
    if (dnu .ne. 0.0_DP) then
      dre = 1.0_DP/dnu
    else
      print *,'SD: NU=0 not allowed! Set dbeta=0 to prevent Stokes operator'// &
              ' from being build!'
      call sys_halt()
    end if
    
    ! If ddelta=0, we have to neglect the nonlinearity. In both cases,
    ! set DlocalDelta=0 which disables the nonlinear term in the assembly.
    ! If dupsam=0, we neglect the stabilisation term (central difference like
    ! discretisation), so we set DlocalDelta=0 as well.
    if ((ddeltamax .eq. 0.0_DP) .or. (dupsam .eq. 0.0_DP)) then
      call lalg_clearVectorDble (DlocalDelta)
    end if
    
    ! Calculate the maximum norm of the actual velocity field
    ! Round up the norm to 1D-8 if it's too small...
    !%OMP SINGLE
    dumax=0.0_DP
    
    do IEQ=1,size(p_DpvelocityX)
      dunorm = sqrt(p_DpvelocityX(IEQ)**2+p_DpvelocityX(IEQ)**2)
      dumax = max(DUMAX,DUNORM)
    end do
  
    !print *,"dumax: ",dumax
    if (dumax.lt.1E-8_DP) dumax=1E-8_DP
    dumaxr = 1.0_DP/dumax
    !%OMP end SINGLE

    ! p_IelementList must point to our set of elements in the discretisation
    ! with that combination of trial/test functions
    call storage_getbase_int (p_relementDistribution%h_IelementList, &
                              p_IelementList)


    ! Loop over the elements - blockwise.
    !
    ! Open-MP-Extension: Each loop cycle is executed in a different thread,
    ! so BILF_NELEMSIM local matrices are simultaneously calculated in the
    ! inner loop(s).
    ! The blocks have all the same size, so we can use static scheduling.
    !%OMP do SCHEDULE(dynamic,1)
    do ielset = 1, size(p_IelementList), BILF_NELEMSIM

      ! We always handle BILF_NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most BILF_NELEMSIM
      ! elements simultaneously.
      
      ielmax = min(size(p_IelementList),ielset-1+BILF_NELEMSIM)
    
      ! The outstanding feature with finite elements is: A basis
      ! function for a DOF on one element has common support only
      ! with the DOF's on the same element! E.g. for Q1:
      !
      !        #. . .#. . .#. . .#
      !        .     .     .     .
      !        .  *  .  *  .  *  .
      !        #-----O-----O. . .#
      !        |     |     |     .
      !        |     | iel |  *  .
      !        #-----X-----O. . .#
      !        |     |     |     .
      !        |     |     |  *  .
      !        #-----#-----#. . .#
      !
      ! --> On element iel, the basis function at "X" only interacts
      !     with the basis functions in "O". Elements in the 
      !     neighbourhood ("*") have no support, therefore we only have
      !     to collect all "O" DOF's.
      !
      ! Calculate the global DOF's into IdofsTrial / IdofsTest.
      !
      ! More exactly, we call dof_locGlobMapping_mult to calculate all the
      ! global DOF's of our BILF_NELEMSIM elements simultaneously.
      call dof_locGlobMapping_mult(p_rdiscretisation, p_IelementList(ielset:ielmax), &
                                  Idofs)
                                  
      ! Calculate local DELTA's for streamline diffusion method.
      ! (cf. p. 121 in Turek's CFD book).
      ! For every element, we need a local DELTA.
      ! Every local delta is weighted by the global "ddelta".
      ! If ddelta=0, we don't do anything as this disables the
      ! nonlinear term.
      ! If UPSAM=0.0, we have a central-difference like discretisation, which
      ! is one can see as the local stabilisation weight Delta is also = 0.0.
      ! In this case, we even switch of the calculation of the local Delta,
      ! as it is always =0.0, so we save a little bit time.
      if ((ddeltaMax .ne. 0.0_DP) .and. (dupsam .ne. 0.0_DP))then
        do iel=1,ielmax-ielset+1
          call getLocalDeltaQuad (p_DpvelocityX,p_DpvelocityY,p_DpvelocityX,p_DpvelocityY,&
                      1.0_DP,0.0_DP, &
                      int(iel+ielset-1,PREC_ELEMENTIDX),DUMAXR,DlocalDelta(iel), &
                      p_IverticesAtElement,p_DvertexCoords,Idofs(:,iel),indof, &
                      dupsam,dre)
        end do ! iel
      end if
                                   
      ! For the assembly of the global matrix, we use a "local"
      ! approach. At first we build a "local" system matrix according
      ! to the current element. This contains all additive
      ! contributions of element iel, which are later added at the
      ! right positions to the elements in the global system matrix.
      !
      ! We have indofTrial trial DOF's per element and
      ! indofTest test DOF's per element. Therefore there are
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
      do iel=1,ielmax-ielset+1
      
        ! For building the local matrices, we have first to
        ! loop through the test functions (the "O"'s), as these
        ! define the rows in the matrix.
        do idofe=1,indof
        
          ! Row idofe of the local matrix corresponds 
          ! to row=global DOF KDFG(idofe) in the global matrix.
          ! This is one of the the "O"'s in the above picture.
          ! Get the starting position of the corresponding row
          ! to jcol0:

          jcol0=p_KLD(Idofs(idofe,iel))
          
          ! Now we loop through the other DOF's on the current element
          ! (the "O"'s).
          ! All these have common support with our current basis function
          ! and will therefore give an additive value to the global
          ! matrix.
          
          do jdofe=1,indof
            
            ! Get the global DOF of the "X" which interacts with 
            ! our "O".
            
            jdfg=Idofs(jdofe,iel)
            
            ! Starting in jcol0 (which points to the beginning of
            ! the line initially), loop through the elements in
            ! the row to find the position of column IDFG.
            ! Jump out of the do loop if we find the column.
            
            do jcol=jcol0,rmatrix%RmatrixBlock(1,1)%NA
              if (p_KCOL(jcol) .eq. jdfg) exit
            end do
            
            ! Because columns in the global matrix are sorted 
            ! ascendingly (except for the diagonal element),
            ! the next search can start after the column we just found.
            
            ! jcol0=jcol+1
            
            ! Save the position of the matrix entry into the local
            ! matrix.
            ! Note that a column in Kentry corresponds to a row in
            ! the real matrix. We aligned Kentry/DENTRY this way to get
            ! higher speed of the assembly routine, since this leads
            ! to better data locality.
            
            Kentry(jdofe,idofe,iel)=jcol
            
          end do ! idofe
          
        end do ! jdofe
        
      end do ! iel
      
      ! If the Newton part is to be calculated, we also need the matrix positions
      ! in A12 and A21. We can skip this part if the column structure is
      ! exactly the same!
      if (associated(p_Kcol,p_Kcol12)) then
      
        Kentry12(:,:,:) = Kentry(:,:,:)
        
      else

        do iel=1,ielmax-ielset+1
        
          ! For building the local matrices, we have first to
          ! loop through the test functions (the "O"'s), as these
          ! define the rows in the matrix.
          do idofe=1,indof
          
            ! Row idofe of the local matrix corresponds 
            ! to row=global DOF KDFG(idofe) in the global matrix.
            ! This is one of the the "O"'s in the above picture.
            ! Get the starting position of the corresponding row
            ! to jcol0:

            jcol0=p_KLD12(Idofs(idofe,iel))
            
            ! Now we loop through the other DOF's on the current element
            ! (the "O"'s).
            ! All these have common support with our current basis function
            ! and will therefore give an additive value to the global
            ! matrix.
            
            do jdofe=1,indof
              
              ! Get the global DOF of the "X" which interacts with 
              ! our "O".
              
              jdfg=Idofs(jdofe,iel)
              
              ! Starting in jcol0 (which points to the beginning of
              ! the line initially), loop through the elements in
              ! the row to find the position of column IDFG.
              ! Jump out of the do loop if we find the column.
              
              do jcol=jcol0,rmatrix%RmatrixBlock(1,2)%NA
                if (p_KCOL12(jcol) .eq. jdfg) exit
              end do

              ! Because columns in the global matrix are sorted 
              ! ascendingly (except for the diagonal element),
              ! the next search can start after the column we just found.
              
              ! jcol0=jcol+1
              
              ! Save the position of the matrix entry into the local
              ! matrix.
              ! Note that a column in Kentry corresponds to a row in
              ! the real matrix. We aligned Kentry/DENTRY this way to get
              ! higher speed of the assembly routine, since this leads
              ! to better data locality.
              
              Kentry12(jdofe,idofe,iel)=jcol
              
            end do ! idofe
            
          end do ! jdofe
          
        end do ! iel
        
      end if

      
      ! Ok, we found the positions of the local matrix entries
      ! that we have to change.
      ! To calculate the matrix contributions, we have to evaluate
      ! the elements to give us the values of the basis functions
      ! in all the DOF's in all the elements in our set.

      ! Get the element evaluation tag of all FE spaces. We need it to evaluate
      ! the elements later. All of them can be combined with OR, what will give
      ! a combined evaluation tag. 
      cevaluationTag = elem_getEvaluationTag(p_relementDistribution%celement)
      cevaluationTag = ior(cevaluationTag,elem_getEvaluationTag(EL_Q1))
                      
      ! In the first loop, calculate the coordinates on the reference element.
      ! In all later loops, use the precalculated information.
      !
      ! Note: Why not using
      !   if (ielset .EQ. 1) then
      ! here, but this strange concept with the boolean variable?
      ! Because the if-command does not work with OpenMP! bcubPtsInitialised
      ! is a local variable and will therefore ensure that every thread
      ! is initialising its local set of cubature points!
      if (.not. bcubPtsInitialised) then
        bcubPtsInitialised = .true.
        cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REFPOINTS)
      else
        cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
      end if
      
      ! Calculate all information that is necessary to evaluate the finite element
      ! on all cells of our subset. This includes the coordinates of the points
      ! on the cells.
      call elprep_prepareSetForEvaluation (revalElementSet,&
          cevaluationTag, p_rtriangulation, p_IelementList(ielset:ielmax), &
          ctrafoType, p_DcubPtsRef(:,1:ncubp))
      p_Ddetj => revalElementSet%p_Ddetj

      ! Calculate the values of the basis functions.
      ! Pass p_DcubPts as point coordinates, which point either to the
      ! coordinates on the reference element (the same for all elements)
      ! or on the real element - depending on whether this is a 
      ! parametric or nonparametric element.
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
      ! DOF's as:
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
      do iel=1,ielmax-ielset+1
      
        ! Loop over all cubature points on the current element
        do icubp = 1, ncubp
        
          du1loc = 0.0_DP
          du2loc = 0.0_DP
        
          ! Perform a loop through the trial DOF's.
          do jdofe=1,indof

            ! Get the value of the (test) basis function 
            ! phi_i (our "O") in the cubature point:
            
            db = Dbas(jdofe,1,icubp,iel)
            
            ! Sum up to the value in the cubature point
            
            jdfg = Idofs(jdofe,iel)
            du1loc = du1loc +p_DpvelocityX(jdfg)*db
            du2loc = du2loc +p_DpvelocityY(jdfg)*db

          end do ! jdofe
          
          ! Save the computed velocity
          Dvelocity(1,icubp,iel) = du1loc
          Dvelocity(2,icubp,iel) = du2loc
        
        end do ! icubp
        
      end do ! iel
      
      ! Compute X- and Y-derivative of the velocity?
      do iel=1,ielmax-ielset+1
      
        ! Loop over all cubature points on the current element
        do icubp = 1, ncubp
        
          du1locx = 0.0_DP
          du1locy = 0.0_DP
          du2locx = 0.0_DP
          du2locy = 0.0_DP
        
          ! Perform a loop through the trial DOF's.
          do jdofe=1,indof

            ! Get the value of the (trial) basis function 
            ! phi_i in the cubature point:
            dbx = Dbas(jdofe,DER_DERIV_X,icubp,iel)
            dby = Dbas(jdofe,DER_DERIV_Y,icubp,iel)

            ! Sum up to the value in the cubature point
            jdfg = Idofs(jdofe,iel)
            du1locx = du1locx + p_DpvelocityX(jdfg)*dbx
            du1locy = du1locy + p_DpvelocityY(jdfg)*dby
            du2locx = du2locx + p_DpvelocityX(jdfg)*dbx
            du2locy = du2locy + p_DpvelocityY(jdfg)*dby

          end do ! jdofe
          
          
          ! Save the computed velocity derivative
          DvelocityUderiv(1,icubp,iel) = du1locx
          DvelocityUderiv(2,icubp,iel) = du1locy
          DvelocityVderiv(1,icubp,iel) = du2locx
          DvelocityVderiv(2,icubp,iel) = du2locy
        
        end do ! icubp
        
      end do ! iel
      
      ! Ok, we now use Dvelocity as coefficient array in the assembly
      ! of a bilinear form!
      !
      ! Clear the local matrices. If the Newton part is to be calculated,
      ! we must clear everything, otherwise only Dentry.
      DentryA11 = 0.0_DP
      DentryA12 = 0.0_DP
      DentryA21 = 0.0_DP
      DentryA22 = 0.0_DP
      
      DentryA44 = 0.0_DP
      DentryA45 = 0.0_DP
      DentryA54 = 0.0_DP
      DentryA55 = 0.0_DP
      
      DentryA41 = 0.0_DP
      DentryA52 = 0.0_DP
      DentryA42 = 0.0_DP
      DentryA51 = 0.0_DP
      
      DentryA14 = 0.0_DP
      DentryA25 = 0.0_DP
      DentryA24 = 0.0_DP
      DentryA15 = 0.0_DP
      
      if (dalphaMax .ne. 0.0_DP) then
        ! Mass matrix
        call computeLocalMatricesDiag (Dbas,Domega,p_Ddetj,indof,ncubp,ielmax-ielset+1,&
            Dvelocity,DvelocityUderiv,DvelocityVderiv,DlocalDelta,&
            DentryA11,DentryA22,DentryA44,DentryA55,&
            computeFormMass,dalpha,dnu)
      end if
      
      if ((dbetaMax .ne. 0.0_DP) .and. (dnu .ne. 0.0_DP)) then
        ! Laplace matrix
        call computeLocalMatricesDiag (Dbas,Domega,p_Ddetj,indof,ncubp,ielmax-ielset+1,&
            Dvelocity,DvelocityUderiv,DvelocityVderiv,DlocalDelta,&
            DentryA11,DentryA22,DentryA44,DentryA55,&
            computeFormStokes,dbeta,dnu)
      end if
      
      if (ddeltaMax .ne. 0.0_DP) then
        ! Convection matrix u*grad
        call computeLocalMatricesDiag (Dbas,Domega,p_Ddetj,indof,ncubp,ielmax-ielset+1,&
            Dvelocity,DvelocityUderiv,DvelocityVderiv,DlocalDelta,&
            DentryA11,DentryA22,DentryA44,DentryA55,&
            computeFormConvection,ddelta,dnu)
      end if
      
      if (dnewtonMax .ne. 0.0_DP) then
        ! Newton matrix (.)*grad(u)
        call computeLocalMatricesFullDiag (Dbas,Domega,p_Ddetj,indof,ncubp,ielmax-ielset+1,&
            Dvelocity,DvelocityUderiv,DvelocityVderiv,DlocalDelta,&
            DentryA11,DentryA22,DentryA44,DentryA55,&
            DentryA12,DentryA21,DentryA45,DentryA54,&
            computeFormNewton,dnewton,dnu)
      end if
      
      if (ddeltaTransposedMax .ne. 0.0_DP) then
        ! Convection matrix u^t*grad
        call computeLocalMatricesFullDiag (Dbas,Domega,p_Ddetj,indof,ncubp,ielmax-ielset+1,&
            Dvelocity,DvelocityUderiv,DvelocityVderiv,DlocalDelta,&
            DentryA11,DentryA22,DentryA44,DentryA55,&
            DentryA12,DentryA21,DentryA45,DentryA54,&
            computeFormConvectionTransposed,ddeltaTransposed,dnu)
      end if
      
      if (dnewtonTransposedMax .ne. 0.0_DP) then
        ! Newton matrix (.)*grad(u)^t
        call computeLocalMatricesFullDiag (Dbas,Domega,p_Ddetj,indof,ncubp,ielmax-ielset+1,&
            Dvelocity,DvelocityUderiv,DvelocityVderiv,DlocalDelta,&
            DentryA11,DentryA22,DentryA44,DentryA55,&
            DentryA12,DentryA21,DentryA45,DentryA54,&
            computeFormNewtonTransposed,dnewtonTransposed,dnu)
      end if
        
      ! Now we have set up "local" system matrices. We can either    
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

!      if (iand(cdef,CONV_MODMATRIX) .ne. 0) then
!      
!        ! With or without Newton?
!        if ((dnewton .eq. 0.0_DP) .and. (dnewtonTransposed .eq. 0.0_DP)) then
!        
!          ! Include the local matrices into the global system matrix,
!          ! subblock A11 and (if different from A11) also into A22.
!          !%OMP CRITICAL
!          do iel=1,ielmax-ielset+1
!            do idofe=1,indof
!              do jdofe=1,indof
!                p_Da11(Kentry(jdofe,idofe,iel)) = p_Da11(Kentry(jdofe,idofe,iel)) + &
!                    dtheta * Dentry(jdofe,idofe,iel)
!              end do
!            end do
!          end do
!          !%OMP end CRITICAL
!          
!          if (.not. associated(p_Da11,p_Da22)) then
!            !%OMP CRITICAL
!            do iel=1,ielmax-ielset+1
!              do idofe=1,indof
!                do jdofe=1,indof
!                  p_Da22(Kentry(jdofe,idofe,iel)) = &
!                      p_Da22(Kentry(jdofe,idofe,iel)) + &
!                      dtheta * Dentry(jdofe,idofe,iel)
!                end do
!              end do
!            end do
!            !%OMP end CRITICAL
!
!          end if
!
!        else
!
!          ! Include the local matrices into the global system matrix,
!          ! subblock A11 and A22 (both must exist and be independent from
!          ! each other).
!          !%OMP CRITICAL
!          do iel=1,ielmax-ielset+1
!            do idofe=1,indof
!              do jdofe=1,indof
!                ! Kentry (:,:,:) -> positions of local matrix in A11 and A22.
!                !
!                ! DentryA11 (:,:,:) -> Newton part of A11
!                p_Da11(Kentry(jdofe,idofe,iel)) = p_Da11(Kentry(jdofe,idofe,iel)) + &
!                    dtheta * ( Dentry(jdofe,idofe,iel) + &
!                               DentryA11(jdofe,idofe,iel) )
!
!                ! DentryA22 (:,:,:) -> Newton part of A22
!                p_Da22(Kentry(jdofe,idofe,iel)) = p_Da22(Kentry(jdofe,idofe,iel)) + &
!                    dtheta * ( Dentry(jdofe,idofe,iel) + &
!                               DentryA22(jdofe,idofe,iel) )
!              end do
!            end do
!          end do
!          !%OMP end CRITICAL
!          
!          !%OMP CRITICAL
!          ! Include the local Newton matrix parts into A12 and A21.
!          do iel=1,ielmax-ielset+1
!            do idofe=1,indof
!              do jdofe=1,indof
!                ! Kentry12 (:,:,:) -> positions of local matrix in A12 and A21.
!                !
!                ! Dentry (:,:,:) -> Newton part of A12
!                p_Da12(Kentry12(jdofe,idofe,iel)) = p_Da12(Kentry12(jdofe,idofe,iel)) + &
!                    dtheta * DentryA12(jdofe,idofe,iel) 
!
!                ! Dentry (:,:,:) -> Newton part of A21
!                p_Da21(Kentry12(jdofe,idofe,iel)) = p_Da21(Kentry12(jdofe,idofe,iel)) + &
!                    dtheta * DentryA21(jdofe,idofe,iel) 
!              end do
!            end do
!          end do
!          !%OMP end CRITICAL
!
!        end if        
!        
!      end if
!      
!      ! For cdef containing CONV_MODDEFECT, build the defect vector                     
!      !     D = RHS - A*U                                         
!      ! This is done matrix free, only with the help of the local 
!      ! matrix.                                                   
!      ! In this case, D=(D1,D2) is expected to be the RHS on      
!      ! entry and will be updated to be the defect vector when    
!      ! this routine is left.                                     
!
!      if (iand(cdef,CONV_MODDEFECT) .ne. 0) then
!        
!        ! With or without Newton?
!        if ((dnewton .eq. 0.0_DP) .and. (dnewtonTransposed .eq. 0.0_DP)) then
!          !%OMP CRITICAL
!          do iel=1,ielmax-ielset+1
!            do idofe=1,indof
!
!              IDFG=Idofs(idofe,iel)
!
!              do jdofe=1,indof
!
!                denth = dtheta*Dentry(jdofe,idofe,iel)         
!      
!                jdfg=Idofs(jdofe,iel)
!                Ddef1(IDFG)= Ddef1(IDFG) - denth*Du1(jdfg)
!                Ddef2(IDFG)= Ddef2(IDFG) - denth*Du2(jdfg)
!
!              end do
!            end do
!          end do
!          !%OMP end CRITICAL
!        else
!          !%OMP CRITICAL
!          do iel=1,ielmax-ielset+1
!            do idofe=1,indof
!
!              IDFG=Idofs(idofe,iel)
!
!              do jdofe=1,indof
!
!                denth = dtheta*Dentry(jdofe,idofe,iel)         
!      
!                jdfg=Idofs(jdofe,iel)
!                Ddef1(IDFG)= Ddef1(IDFG) - denth*Du1(jdfg)
!                Ddef2(IDFG)= Ddef2(IDFG) - denth*Du2(jdfg)
!                
!                ! Newton part
!                Ddef1(IDFG)= Ddef1(IDFG) &
!                           - dtheta*DentryA11(jdofe,idofe,iel)*Du1(jdfg) &
!                           - dtheta*DentryA12(jdofe,idofe,iel)*Du2(jdfg)
!                Ddef2(IDFG)= Ddef2(IDFG) &
!                           - dtheta*DentryA21(jdofe,idofe,iel)*Du1(jdfg) &
!                           - dtheta*DentryA22(jdofe,idofe,iel)*Du2(jdfg)
!
!              end do
!            end do
!          end do
!          !%OMP end CRITICAL          
!        end if
!
!      end if
            

    end do ! ielset
    !%OMP end do 
    
    ! Release memory
    call elprep_releaseElementSet(revalElementSet)

    deallocate(DlocalDelta)
    
    deallocate(DentryA11)
    deallocate(DentryA12)
    deallocate(DentryA21)
    deallocate(DentryA22)
               
    deallocate(DentryA44)
    deallocate(DentryA45)
    deallocate(DentryA54)
    deallocate(DentryA55)
               
    deallocate(DentryA41)
    deallocate(DentryA52)
    deallocate(DentryA42)
    deallocate(DentryA51)
               
    deallocate(DentryA14)
    deallocate(DentryA25)
    deallocate(DentryA24)
    deallocate(DentryA15)

    deallocate(Kentry12)
    deallocate(DvelocityUderiv)
    deallocate(DvelocityVderiv)
    deallocate(Dvelocity)
    deallocate(Kentry)
    deallocate(Idofs)
    deallocate(Dbas)
    !%OMP end PARALLEL
    deallocate(p_DcubPtsRef)
    
  end subroutine
  
  subroutine computeFormMass (dbasI,dbasIX,dbasIY,dbasJ,dbasJX,dbasJY,&
      du1,du1x,du1y,du2,du2x,du2y,dnu,dlocalDelta,Dweight,da11,da22,da44,da55)
      
  !<description>  
    ! Compute the form of the mass matrix.
  !</description>
    
  !<input>
    ! Basis function phi_i and its derivatives.
    real(DP), intent(in) :: dbasI,dbasIX,dbasIY
    
    ! Basis function phi_j and its derivatives.
    real(DP), intent(in) :: dbasJ,dbasJX,dbasJY
    
    ! X-velocity of the given solution vector and its derivatives.
    real(DP), intent(in) :: du1,du1x,du1y
    
    ! Y-velocity of the given solution vector and its derivatives.
    real(DP), intent(in) :: du2,du2x,du2y
    
    ! Viscosity parameter
    real(DP), intent(in) :: dnu
    
    ! Local delta 
    real(DP), intent(in) :: dlocalDelta
    
    ! Weight for A11-22, A44-55, A41-52, A14-25
    real(DP), dimension(2,2), intent(in) :: Dweight
  !</input>
    
  !<output>
    ! Values of the form
    real(DP), intent(out) :: da11,da22,da44,da55
  !</output>
  
    real(dp) :: dtemp
    
    ! dalpha*HBASI1*HBASJ1
    dtemp = dbasI*dbasJ
  
    da11 = Dweight(1,1)*dtemp
    da22 = Dweight(1,1)*dtemp
    da44 = Dweight(2,2)*dtemp
    da55 = Dweight(2,2)*dtemp
  
  end subroutine

  subroutine computeFormStokes (dbasI,dbasIX,dbasIY,dbasJ,dbasJX,dbasJY,&
      du1,du1x,du1y,du2,du2x,du2y,dnu,dlocalDelta,Dweight,da11,da22,da44,da55)
      
  !<description>  
    ! Compute the form of the Stokes matrix.
  !</description>
    
  !<input>
    ! Basis function phi_i and its derivatives.
    real(DP), intent(in) :: dbasI,dbasIX,dbasIY
    
    ! Basis function phi_j and its derivatives.
    real(DP), intent(in) :: dbasJ,dbasJX,dbasJY
    
    ! X-velocity of the given solution vector and its derivatives.
    real(DP), intent(in) :: du1,du1x,du1y
    
    ! Y-velocity of the given solution vector and its derivatives.
    real(DP), intent(in) :: du2,du2x,du2y
    
    ! Viscosity parameter
    real(DP), intent(in) :: dnu
    
    ! Local delta 
    real(DP), intent(in) :: dlocalDelta
    
    ! Weight for A11-22, A44-55, A41-52, A14-25
    real(DP), dimension(2,2), intent(in) :: Dweight
  !</input>
    
  !<output>
    ! Values of the form
    real(DP), intent(out) :: da11,da22,da44,da55
  !</output>
  
    real(dp) :: dtemp
    
    ! dny*(grad(phi_j,grad(phi_i))
    dtemp = dnu*(dbasIX*dbasIX+dbasIY*dbasIY)
  
    da11 = Dweight(1,1)*dtemp
    da22 = Dweight(1,1)*dtemp
    da44 = Dweight(2,2)*dtemp
    da55 = Dweight(2,2)*dtemp
  
  end subroutine

  subroutine computeFormConvection (dbasI,dbasIX,dbasIY,dbasJ,dbasJX,dbasJY,&
      du1,du1x,du1y,du2,du2x,du2y,dnu,dlocalDelta,Dweight,da11,da22,da44,da55)
      
  !<description>  
    ! Compute the form of the Convection matrix.
  !</description>
    
  !<input>
    ! Basis function phi_i and its derivatives.
    real(DP), intent(in) :: dbasI,dbasIX,dbasIY
    
    ! Basis function phi_j and its derivatives.
    real(DP), intent(in) :: dbasJ,dbasJX,dbasJY
    
    ! X-velocity of the given solution vector and its derivatives.
    real(DP), intent(in) :: du1,du1x,du1y
    
    ! Y-velocity of the given solution vector and its derivatives.
    real(DP), intent(in) :: du2,du2x,du2y
    
    ! Viscosity parameter
    real(DP), intent(in) :: dnu
    
    ! Local delta 
    real(DP), intent(in) :: dlocalDelta
    
    ! Weight for A11-22, A44-55, A41-52, A14-25
    real(DP), dimension(2,2), intent(in) :: Dweight
  !</input>
    
  !<output>
    ! Values of the form
    real(DP), intent(out) :: da11,da22,da44,da55
  !</output>
  
    real(dp) :: HSUMI,HSUMJ,dtemp
    
    ! Delta*(U*grad(Phi_j), U*grad(Phi_i)) + (U*grad(Phi_j),Phi_i)
    
    HSUMI = dbasIX*du1 + dbasIY*du2
    HSUMJ = dbasJX*du1 + dbasJY*du2
    
    dtemp = HSUMJ*(dlocalDelta*HSUMI+dbasI)
  
    da11 = Dweight(1,1)*dtemp
    da22 = Dweight(1,1)*dtemp
    da44 = Dweight(2,2)*dtemp
    da55 = Dweight(2,2)*dtemp
  
  end subroutine
     
  subroutine computeFormNewton (dbasI,dbasIX,dbasIY,dbasJ,dbasJX,dbasJY,&
      du1,du1x,du1y,du2,du2x,du2y,dnu,dlocalDelta,Dweight,&
      da11,da22,da44,da55,da12,da21,da45,da54)
      
  !<description>  
    ! Compute the form in a cubature point.
  !</description>
    
  !<input>
    ! Basis function phi_i and its derivatives.
    real(DP), intent(in) :: dbasI,dbasIX,dbasIY
    
    ! Basis function phi_j and its derivatives.
    real(DP), intent(in) :: dbasJ,dbasJX,dbasJY
    
    ! X-velocity of the given solution vector and its derivatives.
    real(DP), intent(in) :: du1,du1x,du1y
    
    ! Y-velocity of the given solution vector and its derivatives.
    real(DP), intent(in) :: du2,du2x,du2y
    
    ! Viscosity parameter
    real(DP), intent(in) :: dnu
    
    ! Local delta 
    real(DP), intent(in) :: dlocalDelta

    ! Weight for A11-22, A44-55, A41-52, A14-25
    real(DP), dimension(2,2), intent(in) :: Dweight
  !</input>
    
  !<output>
    ! Values of the form
    real(DP), intent(out) :: da11,da22,da44,da55,da12,da21,da45,da54
  !</output>

  
    real(dp) :: dtemp
    
    dtemp = dbasI*dbasJ
  
    da11 = Dweight(1,1) * du1x * dtemp
    da12 = Dweight(1,1) * du1y * dtemp
    da21 = Dweight(1,1) * du2x * dtemp
    da22 = Dweight(1,1) * du2y * dtemp

    da44 = Dweight(2,2) * du1x * dtemp
    da45 = Dweight(2,2) * du1y * dtemp
    da54 = Dweight(2,2) * du2x * dtemp
    da55 = Dweight(2,2) * du2y * dtemp
  
  end subroutine

  subroutine computeFormConvectionTransposed (dbasI,dbasIX,dbasIY,dbasJ,dbasJX,dbasJY,&
      du1,du1x,du1y,du2,du2x,du2y,dnu,dlocalDelta,Dweight,&
      da11,da22,da44,da55,da12,da21,da45,da54)
      
  !<description>  
    ! Compute the form in a cubature point.
  !</description>
    
  !<input>
    ! Basis function phi_i and its derivatives.
    real(DP), intent(in) :: dbasI,dbasIX,dbasIY
    
    ! Basis function phi_j and its derivatives.
    real(DP), intent(in) :: dbasJ,dbasJX,dbasJY
    
    ! X-velocity of the given solution vector and its derivatives.
    real(DP), intent(in) :: du1,du1x,du1y
    
    ! Y-velocity of the given solution vector and its derivatives.
    real(DP), intent(in) :: du2,du2x,du2y
    
    ! Viscosity parameter
    real(DP), intent(in) :: dnu
    
    ! Local delta 
    real(DP), intent(in) :: dlocalDelta

    ! Weight for A11-22, A44-55, A41-52, A14-25
    real(DP), dimension(2,2), intent(in) :: Dweight
  !</input>
    
  !<output>
    ! Values of the form
    real(DP), intent(out) :: da11,da22,da44,da55,da12,da21,da45,da54
  !</output>

  
    real(dp) :: dtemp1,dtemp2
    
    dtemp1 = dbasJX*dbasI
    dtemp2 = dbasJY*dbasI
  
    da11 = Dweight(1,1) * du1 * dtemp1
    da12 = Dweight(1,1) * du2 * dtemp1
    da21 = Dweight(1,1) * du1 * dtemp2
    da22 = Dweight(1,1) * du2 * dtemp2

    da44 = Dweight(2,2) * du1 * dtemp1
    da45 = Dweight(2,2) * du2 * dtemp1
    da54 = Dweight(2,2) * du1 * dtemp2
    da55 = Dweight(2,2) * du2 * dtemp2
  
  end subroutine
     
  subroutine computeFormNewtonTransposed (dbasI,dbasIX,dbasIY,dbasJ,dbasJX,dbasJY,&
      du1,du1x,du1y,du2,du2x,du2y,dnu,dlocalDelta,Dweight,&
      da11,da22,da44,da55,da12,da21,da45,da54)
      
  !<description>  
    ! Compute the form of the transposed Newton operator.
  !</description>
    
  !<input>
    ! Basis function phi_i and its derivatives.
    real(DP), intent(in) :: dbasI,dbasIX,dbasIY
    
    ! Basis function phi_j and its derivatives.
    real(DP), intent(in) :: dbasJ,dbasJX,dbasJY
    
    ! X-velocity of the given solution vector and its derivatives.
    real(DP), intent(in) :: du1,du1x,du1y
    
    ! Y-velocity of the given solution vector and its derivatives.
    real(DP), intent(in) :: du2,du2x,du2y
    
    ! Viscosity parameter
    real(DP), intent(in) :: dnu
    
    ! Local delta 
    real(DP), intent(in) :: dlocalDelta

    ! Weight for A11-22, A44-55, A41-52, A14-25
    real(DP), dimension(2,2), intent(in) :: Dweight
  !</input>
    
  !<output>
    ! Values of the form
    real(DP), intent(out) :: da11,da22,da44,da55,da12,da21,da45,da54
  !</output>

  
    real(dp) :: dtemp
    
    dtemp = dbasJ*dbasI
  
    da11 = Dweight(1,1) * du1x * dtemp
    da12 = Dweight(1,1) * du2x * dtemp
    da21 = Dweight(1,1) * du1y * dtemp
    da22 = Dweight(1,1) * du2y * dtemp

    da44 = Dweight(2,2) * du1x * dtemp
    da45 = Dweight(2,2) * du2x * dtemp
    da54 = Dweight(2,2) * du1y * dtemp
    da55 = Dweight(2,2) * du2y * dtemp
  
  end subroutine

!  es fehlt: 
!  * Ordentliche Behandlung des R-Terms beim Newton in der dualen Gleichung.
!  Zur Zeit ist die Behandlung noch falsch, das übergebene DWeight-array enthält
!  wahrscheinlich noch die falschen Gewichte! Die Behandlung muss hier und in der
!  Konvektion geschenen.
!  
!  * Bei der neuen SD-Methode: Implementation der berechneten Größen in Defekt/Matrix,
!  ist zur Zeit noch auskommentiert.
!  
!  * Umsetzung des diskreten Newtons
  
  ! ***************************************************************************

  !<subroutine>

  subroutine computeLocalMatricesFullExt (Dbas,Domega,Ddetj,ndof,ncubp,NEL,&
      Dpvel,DpvelXderiv,DpvelYderiv,Ddvel,DdvelXderiv,DdvelYderiv,&
      DentryA11,DentryA22,DentryA44,DentryA55,&
      DentryA12,DentryA21,DentryA45,DentryA54,&
      DentryA41,DentryA52,DentryA42,DentryA51,&
      DentryA14,DentryA25,DentryA24,DentryA15,&
      DlocalDelta,roptcoperator)
      
  !<description>
    ! Computes a local matrix to be incorporated into the global matrix.
    ! Variant for handling the full matrix.
  !</description>
      
  !<input>      
  
    ! Basis functions in all cubature points
    real(DP), dimension(:,:,:,:), intent(in)  :: Dbas
    
    ! Cubature weights
    real(DP), dimension(:), intent(in)        :: Domega
    
    ! Jacobian determinants in all cubature points
    real(DP), dimension(:,:), intent(in)      :: Ddetj
    
    ! Number of local DOF's in each element
    integer, intent(in)                       :: ndof
    
    ! Number of cubature points on each element
    integer, intent(in)                       :: ncubp
    
    ! Number of elements
    integer, intent(in)                       :: NEL
    
    ! Primal velocity vector and its derivatives in the cubature points
    real(DP), dimension(:,:,:), intent(in)    :: Dpvel
    real(DP), dimension(:,:,:), intent(in)    :: DpvelXderiv
    real(DP), dimension(:,:,:), intent(in)    :: DpvelYderiv

    real(DP), dimension(:,:,:), intent(in)    :: Ddvel
    real(DP), dimension(:,:,:), intent(in)    :: DdvelXderiv
    real(DP), dimension(:,:,:), intent(in)    :: DdvelYderiv
    
    ! Local delta of the SD method
    real(DP), dimension(:), intent(in)        :: DlocalDelta
    
    ! Configuration of the operator
    type(t_optcoperator), intent(in)          :: roptcoperator
    
  !</input>
    
  !<output>
    
    ! Entries in the matrix
    real(DP), dimension(:,:,:), intent(inout) :: DentryA11
    real(DP), dimension(:,:,:), intent(inout) :: DentryA22
    real(DP), dimension(:,:,:), intent(inout) :: DentryA44
    real(DP), dimension(:,:,:), intent(inout) :: DentryA55
    real(DP), dimension(:,:,:), intent(inout) :: DentryA12
    real(DP), dimension(:,:,:), intent(inout) :: DentryA21
    real(DP), dimension(:,:,:), intent(inout) :: DentryA45
    real(DP), dimension(:,:,:), intent(inout) :: DentryA54
    
    real(DP), dimension(:,:,:), intent(inout) :: DentryA41
    real(DP), dimension(:,:,:), intent(inout) :: DentryA52
    real(DP), dimension(:,:,:), intent(inout) :: DentryA42
    real(DP), dimension(:,:,:), intent(inout) :: DentryA51
    real(DP), dimension(:,:,:), intent(inout) :: DentryA14
    real(DP), dimension(:,:,:), intent(inout) :: DentryA25
    real(DP), dimension(:,:,:), intent(inout) :: DentryA24
    real(DP), dimension(:,:,:), intent(inout) :: DentryA15
 
  !</output>   
  
  !<subroutine>
    
    ! local variables
    integer :: iel, icubp, idofe, jdofe
    real(DP) :: OM
    real(DP) :: du1p,du2p,du1locxp,du1locyp,du2locxp,du2locyp
    real(DP) :: du1d,du2d,du1locxd,du1locyd,du2locxd,du2locyd
    real(DP) :: HBASI1,HBASI2,HBASI3,HBASJ1,HBASJ2,HBASJ3
    real(DP) :: AH11,AH22,AH44,AH55,AH12,AH21,AH45,AH54
    real(DP) :: AH41,AH52,AH42,AH51,AH14,AH25,AH24,AH15
    
    AH11 = 0.0_DP
    AH22 = 0.0_DP
    AH44 = 0.0_DP
    AH55 = 0.0_DP
    AH12 = 0.0_DP
    AH21 = 0.0_DP
    AH45 = 0.0_DP
    AH54 = 0.0_DP
    
    AH41 = 0.0_DP
    AH52 = 0.0_DP
    AH42 = 0.0_DP
    AH51 = 0.0_DP
    AH14 = 0.0_DP
    AH25 = 0.0_DP
    AH24 = 0.0_DP
    AH15 = 0.0_DP
    
    ! Loop over the elements in the current set.
    do iel=1,NEL

      ! Loop over all cubature points on the current element
      do icubp = 1, ncubp

        ! Calculate the current weighting factor in the cubature formula
        ! in that cubature point.
        !
        ! Normally, we have to take the absolut value of the determinant 
        ! of the mapping here!
        ! In 2D, the determinant is always positive, whereas in 3D,
        ! the determinant might be negative -- that's normal!
        ! But because this routine only works in 2D, we can skip
        ! the ABS here!

        OM = Domega(icubp)*Ddetj(icubp,iel)

        ! Current velocity in this cubature point:
        du1p = Dpvel (1,icubp,iel)
        du2p = Dpvel (2,icubp,iel)
        du1locxp = DpvelXderiv (1,icubp,iel)
        du1locyp = DpvelXderiv (2,icubp,iel)
        du2locxp = DpvelYderiv (1,icubp,iel)
        du2locyp = DpvelYderiv (2,icubp,iel)
        
        du1d = Ddvel (1,icubp,iel)
        du2d = Ddvel (2,icubp,iel)
        du1locxd = DdvelXderiv (1,icubp,iel)
        du1locyd = DdvelXderiv (2,icubp,iel)
        du2locxd = DdvelYderiv (1,icubp,iel)
        du2locyd = DdvelYderiv (2,icubp,iel)
        
        ! Outer loop over the DOF's i=1..indof on our current element, 
        ! which corresponds to the basis functions Phi_i:

        do idofe=1,ndof
        
          ! Fetch the contributions of the (test) basis functions Phi_i
          ! (our "O")  for function value and first derivatives for the 
          ! current DOF into HBASIy:
        
          HBASI1 = Dbas(idofe,1,icubp,iel)
          HBASI2 = Dbas(idofe,2,icubp,iel)
          HBASI3 = Dbas(idofe,3,icubp,iel)
          
          ! Inner loop over the DOF's j=1..indof, which corresponds to
          ! the basis function Phi_j:

          do jdofe=1,ndof
            
            ! Fetch the contributions of the (trial) basis function Phi_j
            ! (out "X") for function value and first derivatives for the 
            ! current DOF into HBASJy:
          
            HBASJ1 = Dbas(jdofe,1,icubp,iel)
            HBASJ2 = Dbas(jdofe,2,icubp,iel)
            HBASJ3 = Dbas(jdofe,3,icubp,iel)

            ! Finally calculate the contribution to the system matrices.
            !call computeForm (HBASI1,HBASI2,HBASI3,HBASJ1,HBASJ2,HBASJ3,&
            !    du1,du2,du1locx,du1locy,du2locx,du2locy,dnu,DlocalDelta(iel),&
            !    Dweight,AH11,AH22,AH44,AH55,AH12,AH21,AH45,AH54,&
            !    AH41,AH52,AH42,AH51,AH14,AH25,AH24,AH15)
                
            ! Weighten the calculated value AHxy by the cubature
            ! weight OM and add it to the local matrices. After the
            ! loop over all DOF's is finished, each entry contains
            ! the calculated integral.

            DentryA11(jdofe,idofe,iel) = DentryA11(jdofe,idofe,iel)+OM*AH11
            DentryA22(jdofe,idofe,iel) = DentryA22(jdofe,idofe,iel)+OM*AH22
            DentryA44(jdofe,idofe,iel) = DentryA44(jdofe,idofe,iel)+OM*AH44
            DentryA55(jdofe,idofe,iel) = DentryA55(jdofe,idofe,iel)+OM*AH55
            
            DentryA12(jdofe,idofe,iel) = DentryA12(jdofe,idofe,iel)+OM*AH12
            DentryA21(jdofe,idofe,iel) = DentryA21(jdofe,idofe,iel)+OM*AH21
            DentryA45(jdofe,idofe,iel) = DentryA45(jdofe,idofe,iel)+OM*AH45
            DentryA54(jdofe,idofe,iel) = DentryA54(jdofe,idofe,iel)+OM*AH54

            DentryA14(jdofe,idofe,iel) = DentryA14(jdofe,idofe,iel)+OM*AH14
            DentryA25(jdofe,idofe,iel) = DentryA25(jdofe,idofe,iel)+OM*AH25
            DentryA24(jdofe,idofe,iel) = DentryA24(jdofe,idofe,iel)+OM*AH24
            DentryA15(jdofe,idofe,iel) = DentryA15(jdofe,idofe,iel)+OM*AH15
                                                                           
            DentryA41(jdofe,idofe,iel) = DentryA41(jdofe,idofe,iel)+OM*AH41
            DentryA52(jdofe,idofe,iel) = DentryA52(jdofe,idofe,iel)+OM*AH52
            DentryA51(jdofe,idofe,iel) = DentryA51(jdofe,idofe,iel)+OM*AH51
            DentryA42(jdofe,idofe,iel) = DentryA42(jdofe,idofe,iel)+OM*AH42
            
          end do ! idofe
          
        end do ! jdofe

      end do ! icubp 
    
    end do ! iel
    
  end subroutine

  ! ***************************************************************************

  !<subroutine>

  subroutine computeLocalOptCMatrices (Dbas,Domega,Ddetj,ndof,ncubp,NEL,&
      Dpvel,DpvelXderiv,DpvelYderiv,Ddvel,DdvelXderiv,DdvelYderiv,&
      DpvelDofs,DdvelDofs,&
      DentryA11,DentryA22,DentryA44,DentryA55,&
      DentryA12,DentryA21,DentryA45,DentryA54,&
      DentryA41,DentryA52,DentryA42,DentryA51,&
      DentryA14,DentryA25,DentryA24,DentryA15,&
      DlocalDelta,roptcoperator)
      
  !<description>
    ! Computes a local matrix of the optimal control operator
    ! which is to be incorporated into the global matrix.
  !</description>
      
  !<input>      
  
    ! Basis functions in all cubature points
    real(DP), dimension(:,:,:,:), intent(in)  :: Dbas
    
    ! Cubature weights
    real(DP), dimension(:), intent(in)        :: Domega
    
    ! Jacobian determinants in all cubature points
    real(DP), dimension(:,:), intent(in)      :: Ddetj
    
    ! Number of local DOF's in each element
    integer, intent(in)                       :: ndof
    
    ! Number of cubature points on each element
    integer, intent(in)                       :: ncubp
    
    ! Number of elements
    integer, intent(in)                       :: NEL
    
    ! Primal velocity vector and its derivatives in the cubature points
    real(DP), dimension(:,:,:), intent(in)    :: Dpvel
    real(DP), dimension(:,:,:), intent(in)    :: DpvelXderiv
    real(DP), dimension(:,:,:), intent(in)    :: DpvelYderiv

    real(DP), dimension(:,:,:), intent(in)    :: Ddvel
    real(DP), dimension(:,:,:), intent(in)    :: DdvelXderiv
    real(DP), dimension(:,:,:), intent(in)    :: DdvelYderiv
    
    ! DOF's in the primal/dual velocity
    real(DP), dimension(:,:,:), intent(in)    :: DpvelDofs
    real(DP), dimension(:,:,:), intent(in)    :: DdvelDofs
    
    ! Local delta of the SD method
    real(DP), dimension(:), intent(in)        :: DlocalDelta
    
    ! Configuration of the operator
    type(t_optcoperator), intent(in)          :: roptcoperator
    
  !</input>
    
  !<output>
    
    ! Entries in the matrix.
    ! All local matrices are saved transposed, i.e. we have
    ! DentryIJ(column,row,element). This has to be taken care of
    ! during the assembly and was introduced to gain higher speed
    ! as it avoids extensive jumps in the memory access.
    real(DP), dimension(:,:,:), intent(inout) :: DentryA11
    real(DP), dimension(:,:,:), intent(inout) :: DentryA22
    real(DP), dimension(:,:,:), intent(inout) :: DentryA44
    real(DP), dimension(:,:,:), intent(inout) :: DentryA55
    
    real(DP), dimension(:,:,:), intent(inout) :: DentryA12
    real(DP), dimension(:,:,:), intent(inout) :: DentryA21
    real(DP), dimension(:,:,:), intent(inout) :: DentryA45
    real(DP), dimension(:,:,:), intent(inout) :: DentryA54
    
    real(DP), dimension(:,:,:), intent(inout) :: DentryA41
    real(DP), dimension(:,:,:), intent(inout) :: DentryA52
    real(DP), dimension(:,:,:), intent(inout) :: DentryA42
    real(DP), dimension(:,:,:), intent(inout) :: DentryA51

    real(DP), dimension(:,:,:), intent(inout) :: DentryA14
    real(DP), dimension(:,:,:), intent(inout) :: DentryA25
    real(DP), dimension(:,:,:), intent(inout) :: DentryA24
    real(DP), dimension(:,:,:), intent(inout) :: DentryA15
 
  !</output>   
  
  !<subroutine>
    
    ! local variables
    integer :: iel, icubp, idofe, jdofe
    real(DP) :: OM,dtemp,dsumI,dsumJ
    real(DP) :: du1p,du2p
    real(DP) :: du1locxp,du1locyp,du2locxp,du2locyp
    real(DP) :: du1d,du2d,du1locxd,du1locyd,du2locxd,du2locyd
    real(DP) :: dbasI,dbasIX,dbasIY,dbasJ,dbasJX,dbasJY
    real(DP) :: AH11,AH22,AH44,AH55,AH12,AH21,AH45,AH54
    real(DP) :: AH41,AH52,AH42,AH51,AH14,AH25,AH24,AH15
    real(dp) :: dtemp1,dtemp2
    
    AH11 = 0.0_DP
    AH22 = 0.0_DP
    AH44 = 0.0_DP
    AH55 = 0.0_DP
    
    AH12 = 0.0_DP
    AH21 = 0.0_DP
    AH45 = 0.0_DP
    AH54 = 0.0_DP
    
    AH41 = 0.0_DP
    AH52 = 0.0_DP
    AH42 = 0.0_DP
    AH51 = 0.0_DP
    
    AH14 = 0.0_DP
    AH25 = 0.0_DP
    AH24 = 0.0_DP
    AH15 = 0.0_DP
    
    ! Compute mass matrix part?
    if ((roptcoperator%dprimalAlpha .ne. 0.0_DP) .or. &
        (roptcoperator%ddualAlpha .ne. 0.0_DP) .or. &
        (roptcoperator%ddualRAlpha .ne. 0.0_DP)) then
    
      ! Loop over the elements in the current set.
      do iel=1,NEL

        ! Loop over all cubature points on the current element
        do icubp = 1, ncubp

          ! Calculate the current weighting factor in the cubature formula
          ! in that cubature point.
          !
          ! Normally, we have to take the absolut value of the determinant 
          ! of the mapping here!
          ! In 2D, the determinant is always positive, whereas in 3D,
          ! the determinant might be negative -- that's normal!
          ! But because this routine only works in 2D, we can skip
          ! the ABS here!

          OM = Domega(icubp)*Ddetj(icubp,iel)

          ! Outer loop over the DOF's i=1..indof on our current element, 
          ! which corresponds to the basis functions Phi_i:

          do idofe=1,ndof
          
            ! Fetch the contributions of the (test) basis functions Phi_i
            ! (our "O")  for function value and first derivatives for the 
            ! current DOF into HBASIy:
          
            dbasI = Dbas(idofe,1,icubp,iel)
            
            ! Inner loop over the DOF's j=1..indof, which corresponds to
            ! the basis function Phi_j:

            do jdofe=1,ndof
              
              ! Fetch the contributions of the (trial) basis function Phi_j
              ! (out "X") for function value and first derivatives for the 
              ! current DOF into HBASJy:
            
              dbasJ = Dbas(jdofe,1,icubp,iel)

              ! Contribution to the mass matrix -- weighted by the cubature weight
              dtemp = OM * dbasI * dbasJ

              ! Weighten the calculated value AHxy by the cubature
              ! weight OM and add it to the local matrices. After the
              ! loop over all DOF's is finished, each entry contains
              ! the calculated integral.

              DentryA11(jdofe,idofe,iel) = DentryA11(jdofe,idofe,iel) + &
                                           roptcoperator%dprimalAlpha * dtemp
                                           
              DentryA22(jdofe,idofe,iel) = DentryA22(jdofe,idofe,iel) + &
                                           roptcoperator%dprimalAlpha * dtemp
                                           
              DentryA44(jdofe,idofe,iel) = DentryA44(jdofe,idofe,iel) + &
                                           roptcoperator%ddualAlpha * dtemp
                                           
              DentryA55(jdofe,idofe,iel) = DentryA55(jdofe,idofe,iel) + &
                                           roptcoperator%ddualAlpha * dtemp

              DentryA41(jdofe,idofe,iel) = DentryA41(jdofe,idofe,iel) + &
                                           roptcoperator%ddualRAlpha * dtemp
                                           
              DentryA52(jdofe,idofe,iel) = DentryA52(jdofe,idofe,iel) + &
                                           roptcoperator%ddualRAlpha * dtemp
              
            end do ! idofe
            
          end do ! jdofe

        end do ! icubp 
      
      end do ! iel
    end if
    
    ! Laplace matrix
    if ((roptcoperator%dprimalBeta .ne. 0.0_DP) .or. &
        (roptcoperator%ddualBeta .ne. 0.0_DP)) then
    
      ! Loop over the elements in the current set.
      do iel=1,NEL

        ! Loop over all cubature points on the current element
        do icubp = 1, ncubp

          ! Calculate the current weighting factor in the cubature formula
          ! in that cubature point.
          !
          ! Normally, we have to take the absolut value of the determinant 
          ! of the mapping here!
          ! In 2D, the determinant is always positive, whereas in 3D,
          ! the determinant might be negative -- that's normal!
          ! But because this routine only works in 2D, we can skip
          ! the ABS here!

          OM = Domega(icubp)*Ddetj(icubp,iel)

          ! Outer loop over the DOF's i=1..indof on our current element, 
          ! which corresponds to the basis functions Phi_i:

          do idofe=1,ndof
          
            ! Fetch the contributions of the (test) basis functions Phi_i
            ! (our "O")  for function value and first derivatives for the 
            ! current DOF into HBASIy:
          
            dbasIX = Dbas(idofe,2,icubp,iel)
            dbasIY = Dbas(idofe,3,icubp,iel)
            
            ! Inner loop over the DOF's j=1..indof, which corresponds to
            ! the basis function Phi_j:

            do jdofe=1,ndof
              
              ! Fetch the contributions of the (trial) basis function Phi_j
              ! (out "X") for function value and first derivatives for the 
              ! current DOF into HBASJy:
            
              dbasJX = Dbas(jdofe,2,icubp,iel)
              dbasJY = Dbas(jdofe,3,icubp,iel)

              ! dny*(grad(phi_j,grad(phi_i))
              dtemp = roptcoperator%dnu * (dbasIX*dbasJX + dbasIY*dbasJY)

              ! Weighten the calculated value AHxy by the cubature
              ! weight OM and add it to the local matrices. After the
              ! loop over all DOF's is finished, each entry contains
              ! the calculated integral.

              DentryA11(jdofe,idofe,iel) = DentryA11(jdofe,idofe,iel) + &
                  OM * roptcoperator%dprimalBeta * dtemp
                  
              DentryA22(jdofe,idofe,iel) = DentryA22(jdofe,idofe,iel) + &
                  OM * roptcoperator%dprimalBeta * dtemp
                  
              DentryA44(jdofe,idofe,iel) = DentryA44(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualBeta * dtemp
                  
              DentryA55(jdofe,idofe,iel) = DentryA55(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualBeta * dtemp
              
            end do ! idofe
            
          end do ! jdofe

        end do ! icubp 
      
      end do ! iel
    end if    

    ! Convection operator grad(.)*y.
    ! This operator implies the transposed Newton \grad(y)^t*(.) in the
    ! dual equation.
    if ((roptcoperator%dprimalDelta .ne. 0.0_DP) .or. &
        (roptcoperator%ddualDelta .ne. 0.0_DP)) then
    
      ! Loop over the elements in the current set.
      do iel=1,NEL

        ! Loop over all cubature points on the current element
        do icubp = 1, ncubp

          ! Calculate the current weighting factor in the cubature formula
          ! in that cubature point.
          !
          ! Normally, we have to take the absolut value of the determinant 
          ! of the mapping here!
          ! In 2D, the determinant is always positive, whereas in 3D,
          ! the determinant might be negative -- that's normal!
          ! But because this routine only works in 2D, we can skip
          ! the ABS here!

          OM = Domega(icubp)*Ddetj(icubp,iel)
          
          du1p = Dpvel(1,icubp,iel)
          du2p = Dpvel(2,icubp,iel)
          du1locxp = DpvelXderiv(1,icubp,iel)
          du1locyp = DpvelXderiv(2,icubp,iel)
          du2locxp = DpvelYderiv(1,icubp,iel)
          du2locyp = DpvelYderiv(2,icubp,iel)

          ! Outer loop over the DOF's i=1..indof on our current element, 
          ! which corresponds to the basis functions Phi_i:

          do idofe=1,ndof
          
            ! Fetch the contributions of the (test) basis functions Phi_i
            ! (our "O")  for function value and first derivatives for the 
            ! current DOF into HBASIy:
          
            dbasI = Dbas(idofe,1,icubp,iel)
            dbasIX = Dbas(idofe,2,icubp,iel)
            dbasIY = Dbas(idofe,3,icubp,iel)
            
            dsumI = dbasIX*du1p + dbasIY*du2p
            
            ! Inner loop over the DOF's j=1..indof, which corresponds to
            ! the basis function Phi_j:

            do jdofe=1,ndof
              
              ! Fetch the contributions of the (trial) basis function Phi_j
              ! (out "X") for function value and first derivatives for the 
              ! current DOF into HBASJy:
            
              dbasJ = Dbas(jdofe,1,icubp,iel)
              dbasJX = Dbas(jdofe,2,icubp,iel)
              dbasJY = Dbas(jdofe,3,icubp,iel)

              dsumJ = dbasJX*du1p + dbasJY*du2p

              ! Delta*(U*grad(Phi_j), U*grad(Phi_i)) + (U*grad(Phi_j),Phi_i)
              dtemp = dsumJ * (DlocalDelta(iel)*dsumI + dbasI)

              ! Transposed Newton in the dual equation.            
              dtemp1 = dbasJ*dbasI
            
              ! Weighten the calculated value AHxy by the cubature
              ! weight OM and add it to the local matrices. After the
              ! loop over all DOF's is finished, each entry contains
              ! the calculated integral.

              DentryA11(jdofe,idofe,iel) = DentryA11(jdofe,idofe,iel) + &
                  OM * roptcoperator%dprimalDelta * dtemp
                  
              DentryA22(jdofe,idofe,iel) = DentryA22(jdofe,idofe,iel) + &
                  OM * roptcoperator%dprimalDelta * dtemp
                  
              DentryA44(jdofe,idofe,iel) = DentryA44(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualDelta * dtemp
                  
              DentryA55(jdofe,idofe,iel) = DentryA55(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualDelta * dtemp


              DentryA44(jdofe,idofe,iel) = DentryA44(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualNewtonTrans * du1locxp * dtemp1
                  
              DentryA45(jdofe,idofe,iel) = DentryA45(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualNewtonTrans * du2locxp * dtemp1
                  
              DentryA54(jdofe,idofe,iel) = DentryA54(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualNewtonTrans * du1locyp * dtemp1
                  
              DentryA55(jdofe,idofe,iel) = DentryA55(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualNewtonTrans * du2locyp * dtemp1
              
            end do ! idofe
            
          end do ! jdofe

        end do ! icubp 
      
      end do ! iel
    end if    
    
    ! Newton operator grad(y)*(.).
    ! This operator implies the reactive operator (.)*grad(\lambda)+(\grad(.)^t)\lambda
    ! in the dual equation.
    if (roptcoperator%dprimalNewton .ne. 0.0_DP) then
    
      ! Loop over the elements in the current set.
      do iel=1,NEL

        ! Loop over all cubature points on the current element
        do icubp = 1, ncubp

          ! Calculate the current weighting factor in the cubature formula
          ! in that cubature point.
          !
          ! Normally, we have to take the absolut value of the determinant 
          ! of the mapping here!
          ! In 2D, the determinant is always positive, whereas in 3D,
          ! the determinant might be negative -- that's normal!
          ! But because this routine only works in 2D, we can skip
          ! the ABS here!

          OM = Domega(icubp)*Ddetj(icubp,iel)
          
          du1locxp = DpvelXderiv(1,icubp,iel)
          du1locyp = DpvelXderiv(2,icubp,iel)
          du2locxp = DpvelYderiv(1,icubp,iel)
          du2locyp = DpvelYderiv(2,icubp,iel)
          
          du1locxd = DdvelXderiv(1,icubp,iel)
          du1locyd = DdvelXderiv(2,icubp,iel)
          du2locxd = DdvelYderiv(1,icubp,iel)
          du2locyd = DdvelYderiv(2,icubp,iel)
          
          du1d = Ddvel(1,icubp,iel)
          du2d = Ddvel(2,icubp,iel)

          ! Outer loop over the DOF's i=1..indof on our current element, 
          ! which corresponds to the basis functions Phi_i:

          do idofe=1,ndof
          
            ! Fetch the contributions of the (test) basis functions Phi_i
            ! (our "O")  for function value and first derivatives for the 
            ! current DOF into HBASIy:
          
            dbasI = Dbas(idofe,1,icubp,iel)
            dbasIX = Dbas(idofe,2,icubp,iel)
            dbasIY = Dbas(idofe,3,icubp,iel)
            
            ! Inner loop over the DOF's j=1..indof, which corresponds to
            ! the basis function Phi_j:

            do jdofe=1,ndof
              
              ! Fetch the contributions of the (trial) basis function Phi_j
              ! (out "X") for function value and first derivatives for the 
              ! current DOF into HBASJy:
            
              dbasJ = Dbas(jdofe,1,icubp,iel)
              dbasJX = Dbas(jdofe,2,icubp,iel)
              dbasJY = Dbas(jdofe,3,icubp,iel)

              ! Delta*(U*grad(Phi_j), U*grad(Phi_i)) + (U*grad(Phi_j),Phi_i)
              dtemp = dbasJ*dbasI
              dtemp1 = dbasJX*dbasI
              dtemp2 = dbasJY*dbasI

              ! Weighten the calculated value AHxy by the cubature
              ! weight OM and add it to the local matrices. After the
              ! loop over all DOF's is finished, each entry contains
              ! the calculated integral.

              ! Newton operator in the primal
              DentryA11(jdofe,idofe,iel) = DentryA11(jdofe,idofe,iel) + &
                  OM * roptcoperator%dprimalNewton * du1locxp * dtemp
                  
              DentryA12(jdofe,idofe,iel) = DentryA12(jdofe,idofe,iel) + &
                  OM * roptcoperator%dprimalNewton * du1locyp * dtemp
                  
              DentryA21(jdofe,idofe,iel) = DentryA21(jdofe,idofe,iel) + &
                  OM * roptcoperator%dprimalNewton * du2locxp * dtemp
                  
              DentryA22(jdofe,idofe,iel) = DentryA22(jdofe,idofe,iel) + &
                  OM * roptcoperator%dprimalNewton * du2locyp * dtemp


              ! In the dual equation, we have to set up the reactive
              ! operator.
              !
              ! Newton part \grad(lambda)*(.)
              DentryA41(jdofe,idofe,iel) = DentryA41(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualRNewton * du1locxd * dtemp
                  
              DentryA42(jdofe,idofe,iel) = DentryA42(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualRNewton * du1locyd * dtemp
                  
              DentryA51(jdofe,idofe,iel) = DentryA51(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualRNewton * du2locxd * dtemp
                  
              DentryA52(jdofe,idofe,iel) = DentryA52(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualRNewton * du2locyd * dtemp

              ! Transposed Convective part grad(lambda)^t*(.)
              DentryA41(jdofe,idofe,iel) = DentryA41(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualRDeltaTrans * du1d * dtemp1
                  
              DentryA42(jdofe,idofe,iel) = DentryA42(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualRDeltaTrans * du2d * dtemp1
                  
              DentryA51(jdofe,idofe,iel) = DentryA51(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualRDeltaTrans * du1d * dtemp2
                  
              DentryA52(jdofe,idofe,iel) = DentryA52(jdofe,idofe,iel) + &
                  OM * roptcoperator%ddualRDeltaTrans * du2d * dtemp2

            end do ! idofe
            
          end do ! jdofe

        end do ! icubp 
      
      end do ! iel
    end if    
    
    ! Compute the control mass matrix part?
    if (roptcoperator%dcontrolWeight .ne. 0.0_DP) then
    
      ! Here, we have two cases. In case this is a standard operator,
      ! we just set up a weighted mass matrix that converts the
      ! dual velocity to the control.
      
      if (roptcoperator%dprimalNewton .eq. 0.0_DP) then
    
        ! Loop over the elements in the current set.
        do iel=1,NEL

          ! Loop over all cubature points on the current element
          do icubp = 1, ncubp

            ! Calculate the current weighting factor in the cubature formula
            ! in that cubature point.
            !
            ! Normally, we have to take the absolut value of the determinant 
            ! of the mapping here!
            ! In 2D, the determinant is always positive, whereas in 3D,
            ! the determinant might be negative -- that's normal!
            ! But because this routine only works in 2D, we can skip
            ! the ABS here!

            OM = Domega(icubp)*Ddetj(icubp,iel)
            
            ! Outer loop over the DOF's i=1..indof on our current element, 
            ! which corresponds to the basis functions Phi_i:

            do idofe=1,ndof
            
              ! Fetch the contributions of the (test) basis functions Phi_i
              ! (our "O")  for function value and first derivatives for the 
              ! current DOF into HBASIy:
            
              dbasI = Dbas(idofe,1,icubp,iel)
              
              ! Inner loop over the DOF's j=1..indof, which corresponds to
              ! the basis function Phi_j:

              do jdofe=1,ndof
                
                ! Fetch the contributions of the (trial) basis function Phi_j
                ! (out "X") for function value and first derivatives for the 
                ! current DOF into HBASJy:
              
                dbasJ = Dbas(jdofe,1,icubp,iel)

                ! Contribution to the mass matrix
                dtemp = dbasI * dbasJ

                ! Weighten the calculated value AHxy by the cubature
                ! weight OM and add it to the local matrices. After the
                ! loop over all DOF's is finished, each entry contains
                ! the calculated integral.

                DentryA14(jdofe,idofe,iel) = DentryA14(jdofe,idofe,iel) + &
                     OM * roptcoperator%dcontrolWeight * roptcoperator%dcontrolMultiplier * dtemp
                                             
                DentryA25(jdofe,idofe,iel) = DentryA25(jdofe,idofe,iel) + &
                     OM * roptcoperator%dcontrolWeight * roptcoperator%dcontrolMultiplier * dtemp
                
              end do ! idofe
              
            end do ! jdofe

          end do ! icubp 
        
        end do ! iel
      
      else
      
        ! In the other case, we have to set up a projected mass matrix.
        ! The question is, which projection to use.
        
        select case (roptcoperator%ccontrolProjection)
        case (0)
          ! No projection. Same loop as above.
          !
          ! Loop over the elements in the current set.
          do iel=1,NEL

            ! Loop over all cubature points on the current element
            do icubp = 1, ncubp

              ! Calculate the current weighting factor in the cubature formula
              ! in that cubature point.
              !
              ! Normally, we have to take the absolut value of the determinant 
              ! of the mapping here!
              ! In 2D, the determinant is always positive, whereas in 3D,
              ! the determinant might be negative -- that's normal!
              ! But because this routine only works in 2D, we can skip
              ! the ABS here!

              OM = Domega(icubp)*Ddetj(icubp,iel)

              ! Outer loop over the DOF's i=1..indof on our current element, 
              ! which corresponds to the basis functions Phi_i:

              do idofe=1,ndof
              
                ! Fetch the contributions of the (test) basis functions Phi_i
                ! (our "O")  for function value and first derivatives for the 
                ! current DOF into HBASIy:
              
                dbasI = Dbas(idofe,1,icubp,iel)
                
                ! Inner loop over the DOF's j=1..indof, which corresponds to
                ! the basis function Phi_j:

                do jdofe=1,ndof
                  
                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! (out "X") for function value and first derivatives for the 
                  ! current DOF into HBASJy:
                
                  dbasJ = Dbas(jdofe,1,icubp,iel)

                  ! Contribution to the mass matrix
                  dtemp = roptcoperator%dcontrolWeight * roptcoperator%dcontrolMultiplier * &
                      dbasI * dbasJ

                  ! Weighten the calculated value AHxy by the cubature
                  ! weight OM and add it to the local matrices. After the
                  ! loop over all DOF's is finished, each entry contains
                  ! the calculated integral.

                  DentryA14(jdofe,idofe,iel) = DentryA14(jdofe,idofe,iel) + OM * dtemp
                  DentryA25(jdofe,idofe,iel) = DentryA25(jdofe,idofe,iel) + OM * dtemp
                  
                end do ! idofe
                
              end do ! jdofe

            end do ! icubp 
          
          end do ! iel
         
        case (1)
          ! DOF-based projection. At first, set up the standard matrix as above.
          !
          ! Loop over the elements in the current set.
          do iel=1,NEL

            ! Loop over all cubature points on the current element
            do icubp = 1, ncubp

              ! Calculate the current weighting factor in the cubature formula
              ! in that cubature point.
              !
              ! Normally, we have to take the absolut value of the determinant 
              ! of the mapping here!
              ! In 2D, the determinant is always positive, whereas in 3D,
              ! the determinant might be negative -- that's normal!
              ! But because this routine only works in 2D, we can skip
              ! the ABS here!

              OM = Domega(icubp)*Ddetj(icubp,iel)

              ! Outer loop over the DOF's i=1..indof on our current element, 
              ! which corresponds to the basis functions Phi_i:

              do idofe=1,ndof
              
                ! Fetch the contributions of the (test) basis functions Phi_i
                ! (our "O")  for function value and first derivatives for the 
                ! current DOF into HBASIy:
              
                dbasI = Dbas(idofe,1,icubp,iel)
                
                ! Check if the DOF violates the bounds. If that's the case,
                ! it must be ignored. At first, create the control from the dual velocity
                dtemp1 = DdvelDofs(idofe,iel,1) * roptcoperator%dcontrolMultiplier
                dtemp2 = DdvelDofs(idofe,iel,2) * roptcoperator%dcontrolMultiplier
                
                ! Inner loop over the DOF's j=1..indof, which corresponds to
                ! the basis function Phi_j:
                if ((dtemp1 .ge. roptcoperator%dmin1) .and.  &
                    (dtemp1 .le. roptcoperator%dmax1)) then
                  do jdofe=1,ndof
                    
                    ! Fetch the contributions of the (trial) basis function Phi_j
                    ! (out "X") for function value and first derivatives for the 
                    ! current DOF into HBASJy:
                  
                    dbasJ = Dbas(jdofe,1,icubp,iel)

                    ! Contribution to the mass matrix
                    dtemp = roptcoperator%dcontrolWeight * roptcoperator%dcontrolMultiplier * &
                        dbasI * dbasJ

                    ! Weighten the calculated value AHxy by the cubature
                    ! weight OM and add it to the local matrices. After the
                    ! loop over all DOF's is finished, each entry contains
                    ! the calculated integral.

                    DentryA14(jdofe,idofe,iel) = DentryA14(jdofe,idofe,iel) + OM * dtemp
                    
                  end do ! idofe
                end if

                if ((dtemp2 .ge. roptcoperator%dmin2) .and.  &
                    (dtemp2 .le. roptcoperator%dmax2)) then
                  do jdofe=1,ndof
                    
                    ! Fetch the contributions of the (trial) basis function Phi_j
                    ! (out "X") for function value and first derivatives for the 
                    ! current DOF into HBASJy:
                  
                    dbasJ = Dbas(jdofe,1,icubp,iel)

                    ! Contribution to the mass matrix
                    dtemp = roptcoperator%dcontrolWeight * roptcoperator%dcontrolMultiplier * &
                        dbasI * dbasJ

                    ! Weighten the calculated value AHxy by the cubature
                    ! weight OM and add it to the local matrices. After the
                    ! loop over all DOF's is finished, each entry contains
                    ! the calculated integral.

                    DentryA25(jdofe,idofe,iel) = DentryA25(jdofe,idofe,iel) + OM * dtemp
                    
                  end do ! idofe
                end if
                
              end do ! jdofe

            end do ! icubp 
          
          end do ! iel
          
        case (2)
          ! Cubature-point based projection.
          !
          ! Loop over the elements in the current set.
          do iel=1,NEL

            ! Loop over all cubature points on the current element
            do icubp = 1, ncubp

              ! Calculate the current weighting factor in the cubature formula
              ! in that cubature point.
              !
              ! Normally, we have to take the absolut value of the determinant 
              ! of the mapping here!
              ! In 2D, the determinant is always positive, whereas in 3D,
              ! the determinant might be negative -- that's normal!
              ! But because this routine only works in 2D, we can skip
              ! the ABS here!

              OM = Domega(icubp)*Ddetj(icubp,iel)

              du1d = roptcoperator%dcontrolMultiplier * Ddvel(1,icubp,iel)
              du2d = roptcoperator%dcontrolMultiplier * Ddvel(2,icubp,iel)

              ! Outer loop over the DOF's i=1..indof on our current element, 
              ! which corresponds to the basis functions Phi_i:

              do idofe=1,ndof
              
                ! Fetch the contributions of the (test) basis functions Phi_i
                ! (our "O")  for function value and first derivatives for the 
                ! current DOF into HBASIy:
              
                dbasI = Dbas(idofe,1,icubp,iel)
                
                ! Inner loop over the DOF's j=1..indof, which corresponds to
                ! the basis function Phi_j:

                do jdofe=1,ndof
                  
                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! (out "X") for function value and first derivatives for the 
                  ! current DOF into HBASJy:
                
                  dbasJ = Dbas(jdofe,1,icubp,iel)

                  ! Contribution to the mass matrix
                  dtemp = roptcoperator%dcontrolWeight * roptcoperator%dcontrolMultiplier * &
                      dbasI * dbasJ
                      
                  ! Calculate the control in the cubature point
                  dtemp1 = dtemp
                  dtemp2 = dtemp
                  
                  ! Compute dtemp1/dtemp2. Set it to dtemp if the control is
                  ! in the bounds, otherwise to 0.
                  if ((du1d .lt. roptcoperator%dmin1) .or. (du1d .gt. roptcoperator%dmax1)) then
                    dtemp1 = 0.0_DP
                  end if

                  if ((du2d .lt. roptcoperator%dmin2) .or. (du2d .gt. roptcoperator%dmax2)) then
                    dtemp2 = 0.0_DP
                  end if

                  ! Weighten the calculated value AHxy by the cubature
                  ! weight OM and add it to the local matrices. After the
                  ! loop over all DOF's is finished, each entry contains
                  ! the calculated integral.

                  DentryA14(jdofe,idofe,iel) = DentryA14(jdofe,idofe,iel) + OM * dtemp1
                  DentryA25(jdofe,idofe,iel) = DentryA25(jdofe,idofe,iel) + OM * dtemp2
                  
                end do ! idofe
                
              end do ! jdofe

            end do ! icubp 
          
          end do ! iel
         
        end select
      
      end if
    end if
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine conv_strdiffOptC2dinitasm (rdiscretisation,ielemDistribution,roptcassemblyinfo)
  
!<description>
  ! Initialise the roptcassemblyinfo structure for the assembly of
  ! local matrices.
!</description>

!<input>
  ! Block discretisation structure that defines the discretisation
  type(t_spatialDiscretisation), intent(in), target :: rdiscretisation
  
  ! Id of the element distribution to be processed.
  integer, intent(in) :: ielemDistribution
  
!</input>
  
!<output>
  ! Structure collecting information about the assembly.
  type(t_optcassemblyinfo), intent(out) :: roptcassemblyinfo
!</output>

!</subroutine>

    ! local variables
    type(t_elementDistribution), pointer :: p_relementDistribution
    integer :: i,k
    
    ! Cubature point coordinates on the reference element
    real(DP), dimension(CUB_MAXCUBP, NDIM3D) :: Dxi

    ! Get the element distribution. Here, we can find information about
    ! the cubature formula etc...
    p_relementDistribution => rdiscretisation%RelementDistr(ielemDistribution)
    
    ! Remember the discretisation and the element distribution
    roptcassemblyinfo%p_rdiscretisation => rdiscretisation
    roptcassemblyinfo%p_relementDistribution => p_relementDistribution
    roptcassemblyinfo%ielemDistribution = ielemDistribution
    
    ! Get the number of local DOF's for trial/test functions.
    ! We assume trial and test functions to be the same.
    roptcassemblyinfo%indof = elem_igetNDofLoc(p_relementDistribution%celement)

    ! For saving some memory in smaller discretisations, we calculate
    ! the number of elements per block. For smaller triangulations,
    ! this is NEL. If there are too many elements, it's at most
    ! BILF_NELEMSIM. This is only used for allocating some arrays.
    roptcassemblyinfo%nelementsPerBlock = min(BILF_NELEMSIM,rdiscretisation%p_rtriangulation%NEL)
    
    ! Get from the trial element space the type of coordinate system
    ! that is used there:
    roptcassemblyinfo%ctrafoType = elem_igetTrafoType(p_relementDistribution%celement)
    
    ! Allocate some memory to hold the cubature points on the reference element
    allocate(roptcassemblyinfo%p_DcubPtsRef(&
        trafo_igetReferenceDimension(roptcassemblyinfo%ctrafoType),CUB_MAXCUBP))

    ! Initialise the cubature formula,
    ! Get cubature weights and point coordinates on the reference element
    call cub_getCubPoints(roptcassemblyinfo%p_relementDistribution%ccubTypeBilForm, &
        roptcassemblyinfo%ncubp, Dxi, roptcassemblyinfo%Domega)
    
    ! Reformat the cubature points; they are in the wrong shape!
    do i=1,roptcassemblyinfo%ncubp
      do k=1,ubound(roptcassemblyinfo%p_DcubPtsRef,1)
        roptcassemblyinfo%p_DcubPtsRef(k,i) = Dxi(i,k)
      end do
    end do
    
    ! Allocate arrays for the values of the test- and trial functions.
    ! This is done here in the size we need it. Allocating it in-advance
    ! with something like
    !  allocate(Dbas(EL_MAXNBAS,EL_MAXNDER,ncubp,nelementsPerBlock))
    ! would lead to nonused memory blocks in these arrays during the assembly, 
    ! which reduces the speed by 50%!
    
    allocate(roptcassemblyinfo%Dbas(&
        roptcassemblyinfo%indof,elem_getMaxDerivative(p_relementDistribution%celement), &
        roptcassemblyinfo%ncubp,roptcassemblyinfo%nelementsPerBlock))
        
    ! Allocate memory for the DOF's of all the elements.
    allocate(roptcassemblyinfo%Idofs(roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    
    ! Allocate memory for array with local DELTA's
    allocate(roptcassemblyinfo%DlocalDelta(roptcassemblyinfo%nelementsPerBlock))

    ! Allocate an array saving the local matrices for all elements
    ! in an element set.
    ! We could also allocate EL_MAXNBAS*EL_MAXNBAS*BILF_NELEMSIM integers
    ! for this local matrix, but this would normally not fit to the cache
    ! anymore! indofTrial*indofTest*BILF_NELEMSIM is normally much smaller!
    !
    ! Kentry (:,:,:) defines the positions of the local matrices
    ! in the submatrices A11 and A22.
    ! KentryA12 (:,:,:) defines the positions of the local matrices
    ! in the submatrices A12 and A21.
    allocate(roptcassemblyinfo%Kentry(&
        roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    
    allocate(roptcassemblyinfo%Kentry12(&
        roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    
    ! Allocate memory for the velocites in the cubature points.
    allocate(roptcassemblyinfo%Dpvel(NDIM2D,roptcassemblyinfo%ncubp,&
        roptcassemblyinfo%nelementsPerBlock))
    allocate(roptcassemblyinfo%Ddvel(NDIM2D,roptcassemblyinfo%ncubp,&
        roptcassemblyinfo%nelementsPerBlock))
    allocate(roptcassemblyinfo%Dcontrol(NDIM2D,roptcassemblyinfo%ncubp,&
        roptcassemblyinfo%nelementsPerBlock))
    allocate(roptcassemblyinfo%DpvelXderiv(NDIM2D,roptcassemblyinfo%ncubp,&
        roptcassemblyinfo%nelementsPerBlock))
    allocate(roptcassemblyinfo%DpvelYderiv(NDIM2D,roptcassemblyinfo%ncubp,&
        roptcassemblyinfo%nelementsPerBlock))
    allocate(roptcassemblyinfo%DdvelXderiv(NDIM2D,roptcassemblyinfo%ncubp,&
        roptcassemblyinfo%nelementsPerBlock))
    allocate(roptcassemblyinfo%DdvelYderiv(NDIM2D,roptcassemblyinfo%ncubp,&
        roptcassemblyinfo%nelementsPerBlock))
    
    ! Allocate memory for the DOF's on the elements.
    allocate (roptcassemblyinfo%DpvelDofs(&
        roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock,NDIM2D))
    allocate (roptcassemblyinfo%DdvelDofs(&
        roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock,NDIM2D))
    
    ! Indicate that cubature points must still be initialised in the element set.
    roptcassemblyinfo%bcubPtsInitialised = .false.
    
    ! Initialise the derivative flags
    roptcassemblyinfo%Bder = .false.
    roptcassemblyinfo%Bder(DER_FUNC) = .true.
    roptcassemblyinfo%Bder(DER_DERIV2D_X) = .true.
    roptcassemblyinfo%Bder(DER_DERIV2D_Y) = .true.
    
    ! Initialisation of the element set.
    call elprep_init(roptcassemblyinfo%revalElementSet)

    ! p_IelementList must point to our set of elements in the discretisation
    ! with that combination of trial/test functions
    call storage_getbase_int (p_relementDistribution%h_IelementList, &
                              roptcassemblyinfo%p_IelementList)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine conv_strdiffOptC2dinitelemset (rvelmatrix,rvelmatrixoffdiag,rvelocityVector,&
      roptcoperator,roptcassemblyinfo,istartElement,iendElement)
  
!<description>
  ! Initialise the roptcassemblyinfo structure for the assembly of
  ! local matrices for the elements istartElement..istartElement+ielCountMax-1 (of the
  ! current element set). 
!</description>

!<input>
  ! Template FE matrix specifying the connectivity of the primal/dual
  ! velocity space for diagonal blocks.
  type(t_matrixScalar), intent(in) :: rvelmatrix

  ! Template FE matrix specifying the connectivity of the primal/dual
  ! velocity space for offdiagonal blocks.
  type(t_matrixScalar), intent(in) :: rvelmatrixoffdiag
  
  ! Velocity vector
  type(t_vectorBlock), intent(in) :: rvelocityVector

  ! Structure defining the operator to set up.
  type(t_optcoperator), intent(in) :: roptcoperator

  ! Index of the start element of the current element set in roptcassemblyinfo.
  integer, intent(in) :: istartElement
  
  ! Index of the last element.
  integer, intent(in) :: iendElement
!</input>

!<inputoutput>
  ! Structure collecting information about the assembly.
  type(t_optcassemblyinfo), intent(inout) :: roptcassemblyinfo
!</inputoutput>

!</subroutine>

    ! local variables
    integer :: iel, idofe, jdofe, jcol0, jcol, icubp, jdfg
    real(DP) :: du1locp,du2locp,du1locd,du2locd,dc1loc,dc2loc,dcx,dcy
    real(DP) :: du1locxp,du2locxp,du1locyp,du2locyp,dbx,dby,db
    real(DP) :: du1locxd,du2locxd,du1locyd,du2locyd
    real(DP) :: dumaxr,dre
    integer(PREC_VECIDX), dimension(:), pointer :: p_Kcol
    integer(PREC_MATIDX), dimension(:), pointer :: p_Kld
    integer(PREC_VECIDX), dimension(:), pointer :: p_Kcol12
    integer(PREC_MATIDX), dimension(:), pointer :: p_Kld12

    ! Primal and dual velocity
    real(DP), dimension(:), pointer :: p_DpvelocityX, p_DpvelocityY
    real(DP), dimension(:), pointer :: p_DdvelocityX, p_DdvelocityY
    
    ! Element evaluation tag; collects some information necessary for evaluating
    ! the elements.
    integer(I32) :: cevaluationTag
    
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    integer(PREC_VERTEXIDX), dimension(:,:), pointer :: p_IedgesAtElement,p_IverticesAtElement
    
    ! Get triangulation information
    call storage_getbase_double2d (roptcassemblyinfo%p_rdiscretisation%&
        p_rtriangulation%h_DvertexCoords,p_DvertexCoords)
    call storage_getbase_int2d (roptcassemblyinfo%p_rdiscretisation%&
        p_rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
    call storage_getbase_int2d (roptcassemblyinfo%p_rdiscretisation%&
        p_rtriangulation%h_IedgesAtElement,p_IedgesAtElement)

    ! Get pointers to the matrix structure(s).
    call lsyssc_getbase_Kcol (rvelmatrix,p_Kcol)
    call lsyssc_getbase_Kld (rvelmatrix,p_Kld)
    
    call lsyssc_getbase_Kcol (rvelmatrixoffdiag,p_Kcol12)
    call lsyssc_getbase_Kld (rvelmatrixoffdiag,p_Kld12)
   
    ! The outstanding feature with finite elements is: A basis
    ! function for a DOF on one element has common support only
    ! with the DOF's on the same element! E.g. for Q1:
    !
    !        #. . .#. . .#. . .#
    !        .     .     .     .
    !        .  *  .  *  .  *  .
    !        #-----O-----O. . .#
    !        |     |     |     .
    !        |     | iel |  *  .
    !        #-----X-----O. . .#
    !        |     |     |     .
    !        |     |     |  *  .
    !        #-----#-----#. . .#
    !
    ! --> On element iel, the basis function at "X" only interacts
    !     with the basis functions in "O". Elements in the 
    !     neighbourhood ("*") have no support, therefore we only have
    !     to collect all "O" DOF's.
    !
    ! Calculate the global DOF's into IdofsTrial / IdofsTest.
    !
    ! More exactly, we call dof_locGlobMapping_mult to calculate all the
    ! global DOF's of our BILF_NELEMSIM elements simultaneously.
    call dof_locGlobMapping_mult(roptcassemblyinfo%p_rdiscretisation, &
        roptcassemblyinfo%p_IelementList(istartElement:iendElement), roptcassemblyinfo%Idofs)
                                
    ! For the assembly of the global matrix, we use a "local"
    ! approach. At first we build a "local" system matrix according
    ! to the current element. This contains all additive
    ! contributions of element iel, which are later added at the
    ! right positions to the elements in the global system matrix.
    !
    ! We have indofTrial trial DOF's per element and
    ! indofTest test DOF's per element. Therefore there are
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
    do iel=1,iendelement-istartElement+1
    
      ! For building the local matrices, we have first to
      ! loop through the test functions (the "O"'s), as these
      ! define the rows in the matrix.
      do idofe=1,roptcassemblyinfo%indof
      
        ! Row idofe of the local matrix corresponds 
        ! to row=global DOF KDFG(idofe) in the global matrix.
        ! This is one of the the "O"'s in the above picture.
        ! Get the starting position of the corresponding row
        ! to jcol0:

        jcol0=p_KLD(roptcassemblyinfo%Idofs(idofe,iel))
        
        ! Now we loop through the other DOF's on the current element
        ! (the "O"'s).
        ! All these have common support with our current basis function
        ! and will therefore give an additive value to the global
        ! matrix.
        
        do jdofe=1,roptcassemblyinfo%indof
          
          ! Get the global DOF of the "X" which interacts with 
          ! our "O".
          
          jdfg=roptcassemblyinfo%Idofs(jdofe,iel)
          
          ! Starting in jcol0 (which points to the beginning of
          ! the line initially), loop through the elements in
          ! the row to find the position of column IDFG.
          ! Jump out of the do loop if we find the column.
          
          do jcol=jcol0,rvelmatrix%NA
            if (p_KCOL(jcol) .eq. jdfg) exit
          end do
          
          ! Because columns in the global matrix are sorted 
          ! ascendingly (except for the diagonal element),
          ! the next search can start after the column we just found.
          
          ! jcol0=jcol+1
          
          ! Save the position of the matrix entry into the local
          ! matrix.
          ! Note that a column in Kentry corresponds to a row in
          ! the real matrix. We aligned Kentry/DENTRY this way to get
          ! higher speed of the assembly routine, since this leads
          ! to better data locality.
          
          roptcassemblyinfo%Kentry(jdofe,idofe,iel)=jcol
          
        end do ! idofe
        
      end do ! jdofe
      
    end do ! iel
    
    ! If the Newton part is to be calculated, we also need the matrix positions
    ! in A12 and A21. We can skip this part if the column structure is
    ! exactly the same!
    if (associated(p_Kcol,p_Kcol12)) then
    
      roptcassemblyinfo%Kentry12(:,:,:) = roptcassemblyinfo%Kentry(:,:,:)
      
    else

      do iel=1,iendelement-istartElement+1
      
        ! For building the local matrices, we have first to
        ! loop through the test functions (the "O"'s), as these
        ! define the rows in the matrix.
        do idofe=1,roptcassemblyinfo%indof
        
          ! Row idofe of the local matrix corresponds 
          ! to row=global DOF KDFG(idofe) in the global matrix.
          ! This is one of the the "O"'s in the above picture.
          ! Get the starting position of the corresponding row
          ! to jcol0:

          jcol0=p_KLD12(roptcassemblyinfo%Idofs(idofe,iel))
          
          ! Now we loop through the other DOF's on the current element
          ! (the "O"'s).
          ! All these have common support with our current basis function
          ! and will therefore give an additive value to the global
          ! matrix.
          
          do jdofe=1,roptcassemblyinfo%indof
            
            ! Get the global DOF of the "X" which interacts with 
            ! our "O".
            
            jdfg=roptcassemblyinfo%Idofs(jdofe,iel)
            
            ! Starting in jcol0 (which points to the beginning of
            ! the line initially), loop through the elements in
            ! the row to find the position of column IDFG.
            ! Jump out of the do loop if we find the column.
            
            do jcol=jcol0,rvelmatrixoffdiag%NA
              if (p_KCOL12(jcol) .eq. jdfg) exit
            end do

            ! Because columns in the global matrix are sorted 
            ! ascendingly (except for the diagonal element),
            ! the next search can start after the column we just found.
            
            ! jcol0=jcol+1
            
            ! Save the position of the matrix entry into the local
            ! matrix.
            ! Note that a column in Kentry corresponds to a row in
            ! the real matrix. We aligned Kentry/DENTRY this way to get
            ! higher speed of the assembly routine, since this leads
            ! to better data locality.
            
            roptcassemblyinfo%Kentry12(jdofe,idofe,iel)=jcol
            
          end do ! idofe
          
        end do ! jdofe
        
      end do ! iel
      
    end if
    
    ! Ok, we found the positions of the local matrix entries
    ! that we have to change.
    ! To calculate the matrix contributions, we have to evaluate
    ! the elements to give us the values of the basis functions
    ! in all the DOF's in all the elements in our set.

    ! Get the element evaluation tag of all FE spaces. We need it to evaluate
    ! the elements later. All of them can be combined with OR, what will give
    ! a combined evaluation tag. 
    cevaluationTag = elem_getEvaluationTag(roptcassemblyinfo%p_relementDistribution%celement)
                    
    ! In the first loop, calculate the coordinates on the reference element.
    ! In all later loops, use the precalculated information.
    !
    ! Note: Why not using
    !   if (ielset .EQ. 1) then
    ! here, but this strange concept with the boolean variable?
    ! Because the if-command does not work with OpenMP! bcubPtsInitialised
    ! is a local variable and will therefore ensure that every thread
    ! is initialising its local set of cubature points!
    if (.not. roptcassemblyinfo%bcubPtsInitialised) then
      roptcassemblyinfo%bcubPtsInitialised = .true.
      cevaluationTag = ior(cevaluationTag,EL_EVLTAG_REFPOINTS)
    else
      cevaluationTag = iand(cevaluationTag,not(EL_EVLTAG_REFPOINTS))
    end if
    
    ! Calculate all information that is necessary to evaluate the finite element
    ! on all cells of our subset. This includes the coordinates of the points
    ! on the cells.
    call elprep_prepareSetForEvaluation (roptcassemblyinfo%revalElementSet,&
        cevaluationTag, roptcassemblyinfo%p_rdiscretisation%p_rtriangulation, &
        roptcassemblyinfo%p_IelementList(istartElement:iendelement), &
        roptcassemblyinfo%ctrafoType, roptcassemblyinfo%p_DcubPtsRef(:,1:roptcassemblyinfo%ncubp))

    ! Calculate the values of the basis functions.
    ! Pass p_DcubPts as point coordinates, which point either to the
    ! coordinates on the reference element (the same for all elements)
    ! or on the real element - depending on whether this is a 
    ! parametric or nonparametric element.
    call elem_generic_sim2 (roptcassemblyinfo%p_relementDistribution%celement, &
        roptcassemblyinfo%revalElementSet, roptcassemblyinfo%Bder, roptcassemblyinfo%Dbas)
          
    ! Calculate the primal and dual velocity and its derivatives
    ! in all the cubature points on all the elements.
    ! The calculation routine might need them.
    
    ! Get pointers to the velocity to be evaluated.
    call lsyssc_getbase_double (rvelocityVector%RvectorBlock(1),p_DpvelocityX)
    call lsyssc_getbase_double (rvelocityVector%RvectorBlock(2),p_DpvelocityY)
    call lsyssc_getbase_double (rvelocityVector%RvectorBlock(4),p_DdvelocityX)
    call lsyssc_getbase_double (rvelocityVector%RvectorBlock(5),p_DdvelocityY)
    
    ! Extract the velocity DOF's omn the elements.
    do iel=1,iendelement-istartElement+1
      do idofe = 1,roptcassemblyinfo%indof
        jdofe = roptcassemblyinfo%Idofs(idofe,iel)
        roptcassemblyinfo%DpvelDofs(idofe,iel,1) = p_DpvelocityX(jdofe)
        roptcassemblyinfo%DpvelDofs(idofe,iel,2) = p_DpvelocityY(jdofe)
        roptcassemblyinfo%DdvelDofs(idofe,iel,1) = p_DdvelocityX(jdofe)
        roptcassemblyinfo%DdvelDofs(idofe,iel,2) = p_DdvelocityY(jdofe)
      end do
    end do
    
    ! Calculate local DELTA's for streamline diffusion method.
    ! (cf. p. 121 in Turek's CFD book).
    ! For every element, we need a local DELTA.
    ! Every local delta is weighted by the global "ddelta".
    ! If ddelta=0, we don't do anything as this disables the
    ! nonlinear term.
    ! If UPSAM=0.0, we have a central-difference like discretisation, which
    ! is one can see as the local stabilisation weight Delta is also = 0.0.
    ! In this case, we even switch of the calculation of the local Delta,
    ! as it is always =0.0, so we save a little bit time.
    dumaxr = 1.0_DP / roptcassemblyinfo%dumax
    dre = 1.0_DP/roptcoperator%dnu
    if (roptcoperator%dupsam .ne. 0.0_DP) then
      do iel=1,iendelement-istartElement+1
        call getLocalDeltaQuad (p_DpvelocityX,p_DpvelocityY,p_DpvelocityX,p_DpvelocityY,&
                    1.0_DP,0.0_DP, &
                    int(iel+istartElement-1,PREC_ELEMENTIDX),dumaxr,&
                    roptcassemblyinfo%DlocalDelta(iel), &
                    p_IverticesAtElement,p_DvertexCoords,roptcassemblyinfo%Idofs(:,iel),&
                    roptcassemblyinfo%indof,roptcoperator%dupsam,dre)
      end do ! iel
    end if
    
    ! Loop over all elements in the current set
    do iel=1,iendelement-istartElement+1
    
      ! Loop over all cubature points on the current element
      do icubp = 1, roptcassemblyinfo%ncubp
      
        ! Primal/dual velocity
        du1locp = 0.0_DP
        du2locp = 0.0_DP
        du1locd = 0.0_DP
        du2locd = 0.0_DP
        
        ! X/Y-derivative of the primal/dual velocity
        du1locxp = 0.0_DP
        du1locyp = 0.0_DP
        du2locxp = 0.0_DP
        du2locyp = 0.0_DP

        du1locxd = 0.0_DP
        du1locyd = 0.0_DP
        du2locxd = 0.0_DP
        du2locyd = 0.0_DP
      
        ! Perform a loop through the trial DOF's.
        do jdofe=1,roptcassemblyinfo%indof

          ! Get the value of the (test) basis function 
          ! phi_i (our "O") in the cubature point:
          
          db = roptcassemblyinfo%Dbas(jdofe,1,icubp,iel)
          dbx = roptcassemblyinfo%Dbas(jdofe,DER_DERIV2D_X,icubp,iel)
          dby = roptcassemblyinfo%Dbas(jdofe,DER_DERIV2D_Y,icubp,iel)
          
          ! Sum up to the value in the cubature point
          
          jdfg = roptcassemblyinfo%Idofs(jdofe,iel)
          du1locp = du1locp + p_DpvelocityX(jdfg)*db
          du2locp = du2locp + p_DpvelocityY(jdfg)*db
          
          du1locd = du1locd + p_DdvelocityX(jdfg)*db
          du2locd = du2locd + p_DdvelocityY(jdfg)*db
          
          du1locxp = du1locxp + p_DpvelocityX(jdfg)*dbx
          du1locyp = du1locyp + p_DpvelocityX(jdfg)*dby
          du2locxp = du2locxp + p_DpvelocityY(jdfg)*dbx
          du2locyp = du2locyp + p_DpvelocityY(jdfg)*dby

          du1locxd = du1locxd + p_DdvelocityX(jdfg)*dbx
          du1locyd = du1locyd + p_DdvelocityX(jdfg)*dby
          du2locxd = du2locxd + p_DdvelocityY(jdfg)*dbx
          du2locyd = du2locyd + p_DdvelocityY(jdfg)*dby

        end do ! jdofe
        
        ! Save the computed velocities and derivatives
        roptcassemblyinfo%Dpvel(1,icubp,iel) = du1locp
        roptcassemblyinfo%Dpvel(2,icubp,iel) = du2locp

        roptcassemblyinfo%Ddvel(1,icubp,iel) = du1locd
        roptcassemblyinfo%Ddvel(2,icubp,iel) = du2locd

        roptcassemblyinfo%DpvelXderiv(1,icubp,iel) = du1locxp
        roptcassemblyinfo%DpvelXderiv(2,icubp,iel) = du1locyp
        roptcassemblyinfo%DpvelYderiv(1,icubp,iel) = du2locxp
        roptcassemblyinfo%DpvelYderiv(2,icubp,iel) = du2locyp
      
        roptcassemblyinfo%DdvelXderiv(1,icubp,iel) = du1locxd
        roptcassemblyinfo%DdvelXderiv(2,icubp,iel) = du1locyd
        roptcassemblyinfo%DdvelYderiv(1,icubp,iel) = du2locxd
        roptcassemblyinfo%DdvelYderiv(2,icubp,iel) = du2locyd

      end do ! icubp
      
    end do ! iel
    
    ! Compute the control in the cubature points.
    ! This is a little bit more tricky as it may involve a
    ! projection.
    select case (roptcoperator%ccontrolProjection)
    case (0)
      ! No projection
      !
      ! Loop over all elements in the current set
      do iel=1,iendelement-istartElement+1
      
        ! Loop over all cubature points on the current element
        do icubp = 1, roptcassemblyinfo%ncubp
        
          ! Control
          dc1loc = 0.0_DP
          dc2loc = 0.0_DP
        
          ! Perform a loop through the trial DOF's.
          do jdofe=1,roptcassemblyinfo%indof

            ! Get the value of the (test) basis function 
            ! phi_i (our "O") in the cubature point:
            
            db = roptcassemblyinfo%Dbas(jdofe,1,icubp,iel)

            ! Convert the dual velocity to the control
            dcx = roptcoperator%dcontrolMultiplier * p_DdvelocityX(jdfg)
            dcy = roptcoperator%dcontrolMultiplier * p_DdvelocityY(jdfg)
            
            ! Sum up to the value in the cubature point
            dc1loc = dc1loc + dcx*db
            dc2loc = dc2loc + dcy*db

          end do ! jdofe
          
          ! Save the computed control
          roptcassemblyinfo%Dcontrol(1,icubp,iel) = dc1loc
          roptcassemblyinfo%Dcontrol(2,icubp,iel) = dc2loc

        end do ! icubp
        
      end do ! iel    

    case (1)
      ! DOF-based projection
      !
      ! Loop over all elements in the current set
      do iel=1,iendelement-istartElement+1
      
        ! Loop over all cubature points on the current element
        do icubp = 1, roptcassemblyinfo%ncubp
        
          ! Control
          dc1loc = 0.0_DP
          dc2loc = 0.0_DP
        
          ! Perform a loop through the trial DOF's.
          do jdofe=1,roptcassemblyinfo%indof

            ! Get the value of the (test) basis function 
            ! phi_i (our "O") in the cubature point:
            
            db = roptcassemblyinfo%Dbas(jdofe,1,icubp,iel)
            
            ! Convert the dual velocity to the control
            dcx = min(max(roptcoperator%dcontrolMultiplier * p_DdvelocityX(jdfg),&
                roptcoperator%dmin1),roptcoperator%dmax1)
            dcy = min(max(roptcoperator%dcontrolMultiplier * p_DdvelocityY(jdfg),&
                roptcoperator%dmin1),roptcoperator%dmax1)
            
            ! Sum up to the value in the cubature point
            dc1loc = dc1loc + dcx*db
            dc2loc = dc2loc + dcy*db

          end do ! jdofe
          
          ! Save the computed control
          roptcassemblyinfo%Dcontrol(1,icubp,iel) = dc1loc
          roptcassemblyinfo%Dcontrol(2,icubp,iel) = dc2loc

        end do ! icubp
        
      end do ! iel    
      
    case (2)
      ! Cubature-point-based projection
      !
      ! Loop over all elements in the current set
      do iel=1,iendelement-istartElement+1
      
        ! Loop over all cubature points on the current element
        do icubp = 1, roptcassemblyinfo%ncubp
        
          ! Control
          dc1loc = 0.0_DP
          dc2loc = 0.0_DP
        
          ! Perform a loop through the trial DOF's.
          do jdofe=1,roptcassemblyinfo%indof

            ! Get the value of the (test) basis function 
            ! phi_i (our "O") in the cubature point:
            
            db = roptcassemblyinfo%Dbas(jdofe,1,icubp,iel)
            
            ! Convert the dual velocity to the control
            dcx = roptcoperator%dcontrolMultiplier * p_DdvelocityX(jdfg)
            dcy = roptcoperator%dcontrolMultiplier * p_DdvelocityY(jdfg)
            
            ! Sum up to the value in the cubature point
            dc1loc = dc1loc + dcx*db
            dc2loc = dc2loc + dcy*db

          end do ! jdofe
          
          ! Save the computed control
          roptcassemblyinfo%Dcontrol(1,icubp,iel) = &
            min(max(dc1loc,roptcoperator%dmin1),roptcoperator%dmax1)
          roptcassemblyinfo%Dcontrol(2,icubp,iel) = &
            min(max(dc2loc,roptcoperator%dmin2),roptcoperator%dmax2)

        end do ! icubp
        
      end do ! iel    
    end select      
      
  end subroutine
  
  ! ***************************************************************************

!<subroutine>

  subroutine conv_strdiffOptC2ddoneasm (roptcassemblyinfo)
  
!<description>
  ! Clean up the roptcassemblyinfo structure.
!</description>

!<inputoutput>
  ! Structure collecting information about the assembly.
  type(t_optcassemblyinfo), intent(inout) :: roptcassemblyinfo
!</inputoutput>

!</subroutine>

    deallocate(roptcassemblyinfo%Dbas)
    deallocate(roptcassemblyinfo%Idofs)
    deallocate(roptcassemblyinfo%DlocalDelta)
    deallocate(roptcassemblyinfo%Kentry)
    deallocate(roptcassemblyinfo%Kentry12)
    
    ! Allocate memory for the velocites in the cubature points.
    deallocate(roptcassemblyinfo%Dpvel)
    deallocate(roptcassemblyinfo%Ddvel)
    deallocate(roptcassemblyinfo%Dcontrol)
    deallocate(roptcassemblyinfo%DpvelXderiv)
    deallocate(roptcassemblyinfo%DpvelYderiv)
    deallocate(roptcassemblyinfo%DdvelXderiv)
    deallocate(roptcassemblyinfo%DdvelYderiv)
    
    deallocate(roptcassemblyinfo%p_DcubPtsRef)

    ! Clean up of the element set.
    call elprep_releaseElementSet(roptcassemblyinfo%revalElementSet)

  end subroutine
  
  ! ***************************************************************************

!<subroutine>
  subroutine conv_strdiffOptC2dgetMatrix (rmatrix,roptcoperator,dweight,&
      rvelocityVector)
!<description>
  ! Calculate the matrix of the nonlinear operator:
  !   dweight*A(rvelocityVector)
!</description>

!<input>

  ! Structure defining the operator to set up.
  type(t_optcoperator), intent(in) :: roptcoperator
  
  ! Weight of the operator. Standard = 1.0
  real(dp), intent(in) :: dweight
      
  ! Velocity vector for the nonlinearity.
  ! The first blocks 1/2 and 4/5 in this vector define the evaluation
  ! point (primal/dual velocity).
  type(t_vectorBlock), intent(in) :: rvelocityVector
  
!</input>

!<inputoutput>
  ! The system matrix. The submatrices for the velocity must be in block
  ! A11, A12, A21 and A22 and must be in matrix format 7 or 9.
  ! A11 and A22 must have the same structure. A12 and A21 must have
  ! the same structure.
  type(t_matrixBlock), intent(inout), target :: rmatrix
!</inputoutput>

!</subroutine>

  ! local variables
  integer(PREC_ELEMENTIDX) :: ielset,ielmax
  
  ! Matrix structure arrays
  real(DP), dimension(:), pointer :: p_Da11,p_Da22,p_Da44,p_Da55
  real(DP), dimension(:), pointer :: p_Da12,p_Da21,p_Da45,p_Da54
  real(DP), dimension(:), pointer :: p_Da41,p_Da51,p_Da42,p_Da52
  real(DP), dimension(:), pointer :: p_Da14,p_Da15,p_Da24,p_Da25
  
  ! Local matrices
  real(DP), dimension(:,:,:), allocatable :: DentryA11
  real(DP), dimension(:,:,:), allocatable :: DentryA12
  real(DP), dimension(:,:,:), allocatable :: DentryA21
  real(DP), dimension(:,:,:), allocatable :: DentryA22
  
  real(DP), dimension(:,:,:), allocatable :: DentryA44
  real(DP), dimension(:,:,:), allocatable :: DentryA45
  real(DP), dimension(:,:,:), allocatable :: DentryA54
  real(DP), dimension(:,:,:), allocatable :: DentryA55
  
  real(DP), dimension(:,:,:), allocatable :: DentryA41
  real(DP), dimension(:,:,:), allocatable :: DentryA52
  real(DP), dimension(:,:,:), allocatable :: DentryA42
  real(DP), dimension(:,:,:), allocatable :: DentryA51
  
  real(DP), dimension(:,:,:), allocatable :: DentryA14
  real(DP), dimension(:,:,:), allocatable :: DentryA25
  real(DP), dimension(:,:,:), allocatable :: DentryA24
  real(DP), dimension(:,:,:), allocatable :: DentryA15
  
  ! Temp vector
  type(t_vectorBlock) :: rvectorBlock
  
  ! Assembly status structure
  type(t_optcassemblyinfo) :: roptcassemblyinfo

    if (roptcoperator%dnu .eq. 0.0_DP) then
      print *,'SD: NU=0 not allowed! Set dbeta=0 to prevent Stokes operator'// &
              ' from being build!'
      call sys_halt()
    end if
    
    ! Get pointers to the matrix content 
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,1),p_Da11)
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,2),p_Da22)

    nullify(p_Da12,p_Da21)
    if (lsysbl_isSubmatrixPresent(rmatrix,1,2)) then
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,2),p_Da12)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,1),p_Da21)
    end if

    call lsyssc_getbase_double (rmatrix%RmatrixBlock(4,4),p_Da44)
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(5,5),p_Da55)

    nullify(p_Da45,p_Da54)
    if (lsysbl_isSubmatrixPresent(rmatrix,4,5)) then
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(4,5),p_Da45)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(5,4),p_Da54)
    end if

    call lsyssc_getbase_double (rmatrix%RmatrixBlock(4,1),p_Da41)
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(5,2),p_Da52)
    
    nullify(p_Da42,p_Da51)
    if (lsysbl_isSubmatrixPresent(rmatrix,4,2)) then
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(4,2),p_Da42)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(5,1),p_Da51)
    end if

    call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,4),p_Da14)
    call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,5),p_Da25)

    nullify(p_Da24,p_Da15)
    if (lsysbl_isSubmatrixPresent(rmatrix,2,4)) then
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(2,4),p_Da24)
      call lsyssc_getbase_double (rmatrix%RmatrixBlock(1,5),p_Da15)
    end if
    
    ! Initialise the asembly of the local matrices.
    call conv_strdiffOptC2dinitasm (rvelocityVector%p_rblockDiscr%RspatialDiscr(1),&
        1,roptcassemblyinfo)

    ! Calculate the maximum norm of the actual velocity field
    ! Round up the norm to 1D-8 if it's too small...
    call lsysbl_deriveSubvector(rvelocityVector,rvectorBlock,1,2,.true.)
    call lsysbl_getVectorMagnitude (rvectorBlock,dumax=roptcassemblyinfo%dumax)
    call lsysbl_releaseVector (rvectorBlock)
  
    if (roptcassemblyinfo%dumax .lt. 1E-8_DP) roptcassemblyinfo%dumax=1E-8_DP

    ! Allocate memory for the local matrices
    allocate(DentryA11(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA12(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA21(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA22(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))

    allocate(DentryA44(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA45(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA54(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA55(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))

    allocate(DentryA41(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA52(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA42(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA51(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))

    allocate(DentryA14(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA25(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA24(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA15(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
  
    ! Initialise the array with the local Delta values for the stabilisation
    call lalg_clearVectorDble (roptcassemblyinfo%DlocalDelta)

    ! Loop over the elements - blockwise.
    !
    ! Open-MP-Extension: Each loop cycle is executed in a different thread,
    ! so BILF_NELEMSIM local matrices are simultaneously calculated in the
    ! inner loop(s).
    ! The blocks have all the same size, so we can use static scheduling.
    !%OMP do SCHEDULE(dynamic,1)
    do ielset = 1, size(roptcassemblyinfo%p_IelementList), BILF_NELEMSIM

      ! We always handle BILF_NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most BILF_NELEMSIM
      ! elements simultaneously.
      
      ielmax = min(size(roptcassemblyinfo%p_IelementList),ielset-1+BILF_NELEMSIM)
    
      ! Initialise the element set, compute the basis functions in the
      ! cubature points.
      if (lsysbl_isSubmatrixPresent(rmatrix,1,2)) then
        call conv_strdiffOptC2dinitelemset (rmatrix%RmatrixBlock(1,1),&
            rmatrix%RmatrixBlock(1,2),rvelocityVector,&
            roptcoperator,roptcassemblyinfo,ielset,ielmax)
      else
        call conv_strdiffOptC2dinitelemset (rmatrix%RmatrixBlock(1,1),&
            rmatrix%RmatrixBlock(1,1),rvelocityVector,&
            roptcoperator,roptcassemblyinfo,ielset,ielmax)
      end if

      ! Clear the local matrices. If the Newton part is to be calculated,
      ! we must clear everything, otherwise only Dentry.
      DentryA11 = 0.0_DP
      DentryA12 = 0.0_DP
      DentryA21 = 0.0_DP
      DentryA22 = 0.0_DP
      
      DentryA44 = 0.0_DP
      DentryA45 = 0.0_DP
      DentryA54 = 0.0_DP
      DentryA55 = 0.0_DP
      
      DentryA41 = 0.0_DP
      DentryA52 = 0.0_DP
      DentryA42 = 0.0_DP
      DentryA51 = 0.0_DP
      
      DentryA14 = 0.0_DP
      DentryA25 = 0.0_DP
      DentryA24 = 0.0_DP
      DentryA15 = 0.0_DP
      
      ! Now calculate the local matrices on all the elements.
      call computeLocalOptCMatrices (roptcassemblyinfo%Dbas,&
          roptcassemblyinfo%Domega,roptcassemblyinfo%revalElementSet%p_Ddetj,&
          roptcassemblyinfo%indof,roptcassemblyinfo%ncubp,ielmax-ielset+1,&
          roptcassemblyinfo%Dpvel,roptcassemblyinfo%DpvelXderiv,&
          roptcassemblyinfo%DpvelYderiv,roptcassemblyinfo%Ddvel,&
          roptcassemblyinfo%DdvelXderiv,roptcassemblyinfo%DdvelYderiv,&
          roptcassemblyinfo%DpvelDofs,roptcassemblyinfo%DdvelDofs,&
          DentryA11,DentryA22,DentryA44,DentryA55,&
          DentryA12,DentryA21,DentryA45,DentryA54,&
          DentryA41,DentryA52,DentryA42,DentryA51,&
          DentryA14,DentryA25,DentryA24,DentryA15,&
          roptcassemblyinfo%DlocalDelta,roptcoperator)      
      
      ! Incorporate the computed local matrices into the global matrix.
      !%OMP CRITICAL
      if (associated(p_Da11)) &
        call incorporateLocalMatrix (DentryA11,p_Da11,roptcassemblyinfo%Kentry,dweight)

      if (associated(p_Da12)) &
        call incorporateLocalMatrix (DentryA12,p_Da12,roptcassemblyinfo%Kentry12,dweight)

      if (associated(p_Da21)) &
        call incorporateLocalMatrix (DentryA21,p_Da21,roptcassemblyinfo%Kentry12,dweight)

      if (associated(p_Da22) .and. .not. associated (p_Da11,p_Da22)) &
        call incorporateLocalMatrix (DentryA22,p_Da22,roptcassemblyinfo%Kentry,dweight)


      if (associated(p_Da44)) &
        call incorporateLocalMatrix (DentryA44,p_Da44,roptcassemblyinfo%Kentry,dweight)

      if (associated(p_Da45)) &
        call incorporateLocalMatrix (DentryA45,p_Da45,roptcassemblyinfo%Kentry12,dweight)

      if (associated(p_Da54)) &
        call incorporateLocalMatrix (DentryA54,p_Da54,roptcassemblyinfo%Kentry12,dweight)

      if (associated(p_Da55)) &
        call incorporateLocalMatrix (DentryA55,p_Da55,roptcassemblyinfo%Kentry,dweight)
      
      
      if (associated(p_Da41)) &
        call incorporateLocalMatrix (DentryA41,p_Da41,roptcassemblyinfo%Kentry,dweight)

      if (associated(p_Da42)) &
        call incorporateLocalMatrix (DentryA42,p_Da42,roptcassemblyinfo%Kentry12,dweight)

      if (associated(p_Da51)) &
        call incorporateLocalMatrix (DentryA51,p_Da51,roptcassemblyinfo%Kentry12,dweight)

      if (associated(p_Da52)) &
        call incorporateLocalMatrix (DentryA52,p_Da52,roptcassemblyinfo%Kentry,dweight)
      
      
      if (associated(p_Da14)) &
        call incorporateLocalMatrix (DentryA14,p_Da14,roptcassemblyinfo%Kentry,dweight)
                                                    
      if (associated(p_Da24)) &                     
        call incorporateLocalMatrix (DentryA24,p_Da24,roptcassemblyinfo%Kentry12,dweight)
                                                    
      if (associated(p_Da15)) &                     
        call incorporateLocalMatrix (DentryA15,p_Da15,roptcassemblyinfo%Kentry12,dweight)
                                                    
      if (associated(p_Da25)) &                     
        call incorporateLocalMatrix (DentryA25,p_Da25,roptcassemblyinfo%Kentry,dweight)
      

    end do ! ielset
    !%OMP end do 
    
    ! Release memory
    deallocate(DentryA11)
    deallocate(DentryA22)
    if (allocated(DentryA12)) deallocate(DentryA12)
    if (allocated(DentryA21)) deallocate(DentryA21)
               
    deallocate(DentryA44)
    deallocate(DentryA55)
    if (allocated(DentryA45)) deallocate(DentryA45)
    if (allocated(DentryA54)) deallocate(DentryA54)
               
    deallocate(DentryA41)
    deallocate(DentryA52)
    if (allocated(DentryA42)) deallocate(DentryA42)
    if (allocated(DentryA51)) deallocate(DentryA51)
               
    deallocate(DentryA14)
    deallocate(DentryA25)
    if (allocated(DentryA24)) deallocate(DentryA24)
    if (allocated(DentryA15)) deallocate(DentryA15)
               
  contains
    
    subroutine incorporateLocalMatrix (DaLocal,Da,Kentry,dweight)
    
    ! Incorporate a local matrix into a global one.
    
    ! The local matrix to incorporate.
    real(dp), dimension(:,:,:), intent(in) :: DaLocal
    
    ! The global matrix
    real(dp), dimension(:), intent(inout) :: Da
    
    ! Positions of the local matrix in the global one.
    integer, dimension(:,:,:), intent(in) :: Kentry
    
    ! Weight of the local matrix
    real(dp) :: dweight

      ! local variables
      integer :: iel,idofe,jdofe
    
      do iel=1,ubound(DaLocal,3)
        do idofe=1,ubound(DaLocal,2)
          do jdofe=1,ubound(DaLocal,1)
            Da(Kentry(jdofe,idofe,iel)) = Da(Kentry(jdofe,idofe,iel)) + &
                dweight * DaLocal(jdofe,idofe,iel) 
          end do
        end do
      end do
    
    end subroutine
    
  end subroutine
     
  ! ***************************************************************************

!<subroutine>
  subroutine conv_strdiffOptC2dgetDefect (rvelMatrix,roptcoperator,rvelocityVector,&
      dweight,rx,rd)
      
!<description>
  ! Calculate the defect of the nonlinear operator:
  !   rd = rd - dweight*A(rvelocityVector)rx 
!</description>

!<input>
  ! Template FE matrix specifying the connectivity in the FE space
  ! of the primal/dual velocity.
  type(t_matrixScalar), intent(in) :: rvelMatrix

  ! Structure defining the operator to set up.
  type(t_optcoperator), intent(in) :: roptcoperator
  
  ! Velocity vector for the nonlinearity.
  ! The first blocks 1/2 and 4/5 in this vector define the evaluation
  ! point (primal/dual velocity).
  type(t_vectorBlock), intent(in) :: rvelocityVector
  
  ! Solution vector, to be multiplied with the matrix.
  type(t_vectorBlock), intent(in) :: rx
  
  ! Weight for the solution
  real(DP), intent(in) :: dweight
!</input>
  
!<inputoutput>  
  ! On entry: RHS. On exit: Defect.
  type(t_vectorBlock), intent(inout) :: rd
!</inputoutput>  

!</subroutine>

  ! local variables
  integer(PREC_ELEMENTIDX) :: ielset,ielmax
  
  ! Local matrices
  real(DP), dimension(:,:,:), allocatable :: DentryA11
  real(DP), dimension(:,:,:), allocatable :: DentryA12
  real(DP), dimension(:,:,:), allocatable :: DentryA21
  real(DP), dimension(:,:,:), allocatable :: DentryA22
  
  real(DP), dimension(:,:,:), allocatable :: DentryA44
  real(DP), dimension(:,:,:), allocatable :: DentryA45
  real(DP), dimension(:,:,:), allocatable :: DentryA54
  real(DP), dimension(:,:,:), allocatable :: DentryA55
  
  real(DP), dimension(:,:,:), allocatable :: DentryA41
  real(DP), dimension(:,:,:), allocatable :: DentryA52
  real(DP), dimension(:,:,:), allocatable :: DentryA42
  real(DP), dimension(:,:,:), allocatable :: DentryA51
  
  real(DP), dimension(:,:,:), allocatable :: DentryA14
  real(DP), dimension(:,:,:), allocatable :: DentryA25
  real(DP), dimension(:,:,:), allocatable :: DentryA24
  real(DP), dimension(:,:,:), allocatable :: DentryA15
  
  ! Temp vector
  type(t_vectorBlock) :: rvectorBlock
  
  ! Pointer to the subvectors
  real(dp), dimension(:), pointer :: p_Dx1,p_Dx2,p_Dx4,p_Dx5
  real(dp), dimension(:), pointer :: p_Dd1,p_Dd2,p_Dd4,p_Dd5
  
  ! Assembly status structure
  type(t_optcassemblyinfo) :: roptcassemblyinfo

    if (roptcoperator%dnu .eq. 0.0_DP) then
      print *,'SD: NU=0 not allowed! Set dbeta=0 to prevent Stokes operator'// &
              ' from being build!'
      call sys_halt()
    end if
    
    ! Initialise the asembly of the local matrices.
    call conv_strdiffOptC2dinitasm (rvelocityVector%p_rblockDiscr%RspatialDiscr(1),&
        1,roptcassemblyinfo)

    ! Calculate the maximum norm of the actual velocity field
    ! Round up the norm to 1D-8 if it's too small...
    call lsysbl_deriveSubvector(rvelocityVector,rvectorBlock,1,2,.true.)
    call lsysbl_getVectorMagnitude (rvectorBlock,dumax=roptcassemblyinfo%dumax)
    call lsysbl_releaseVector (rvectorBlock)
  
    if (roptcassemblyinfo%dumax .lt. 1E-8_DP) roptcassemblyinfo%dumax=1E-8_DP

    ! Allocate memory for the local matrices
    allocate(DentryA11(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA12(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA21(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA22(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))

    allocate(DentryA44(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA45(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA54(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA55(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))

    allocate(DentryA41(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA52(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA42(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA51(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))

    allocate(DentryA14(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA25(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA24(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
    allocate(DentryA15(roptcassemblyinfo%indof,roptcassemblyinfo%indof,roptcassemblyinfo%nelementsPerBlock))
  
    ! Initialise the array with the local Delta values for the stabilisation
    call lalg_clearVectorDble (roptcassemblyinfo%DlocalDelta)
    
    ! Get pointers to the subvectors
    call lsyssc_getbase_double (rx%RvectorBlock(1),p_Dx1)
    call lsyssc_getbase_double (rx%RvectorBlock(2),p_Dx2)
    call lsyssc_getbase_double (rx%RvectorBlock(4),p_Dx4)
    call lsyssc_getbase_double (rx%RvectorBlock(5),p_Dx5)

    call lsyssc_getbase_double (rd%RvectorBlock(1),p_Dd1)
    call lsyssc_getbase_double (rd%RvectorBlock(2),p_Dd2)
    call lsyssc_getbase_double (rd%RvectorBlock(4),p_Dd4)
    call lsyssc_getbase_double (rd%RvectorBlock(5),p_Dd5)

    ! Loop over the elements - blockwise.
    !
    ! Open-MP-Extension: Each loop cycle is executed in a different thread,
    ! so BILF_NELEMSIM local matrices are simultaneously calculated in the
    ! inner loop(s).
    ! The blocks have all the same size, so we can use static scheduling.
    !%OMP do SCHEDULE(dynamic,1)
    do ielset = 1, size(roptcassemblyinfo%p_IelementList), BILF_NELEMSIM

      ! We always handle BILF_NELEMSIM elements simultaneously.
      ! How many elements have we actually here?
      ! Get the maximum element number, such that we handle at most BILF_NELEMSIM
      ! elements simultaneously.
      
      ielmax = min(size(roptcassemblyinfo%p_IelementList),ielset-1+BILF_NELEMSIM)
    
      ! Initialise the element set, compute the basis functions in the
      ! cubature points.
      call conv_strdiffOptC2dinitelemset (rvelMatrix,rvelMatrix,&
          rvelocityVector,roptcoperator,roptcassemblyinfo,ielset,ielmax)

      ! Clear the local matrices. If the Newton part is to be calculated,
      ! we must clear everything, otherwise only Dentry.
      DentryA11 = 0.0_DP
      DentryA12 = 0.0_DP
      DentryA21 = 0.0_DP
      DentryA22 = 0.0_DP
      
      DentryA44 = 0.0_DP
      DentryA45 = 0.0_DP
      DentryA54 = 0.0_DP
      DentryA55 = 0.0_DP
      
      DentryA41 = 0.0_DP
      DentryA52 = 0.0_DP
      DentryA42 = 0.0_DP
      DentryA51 = 0.0_DP
      
      DentryA14 = 0.0_DP
      DentryA25 = 0.0_DP
      DentryA24 = 0.0_DP
      DentryA15 = 0.0_DP
      
      ! Now calculate the local matrices on all the elements.
      call computeLocalOptCMatrices (roptcassemblyinfo%Dbas,&
          roptcassemblyinfo%Domega,roptcassemblyinfo%revalElementSet%p_Ddetj,&
          roptcassemblyinfo%indof,roptcassemblyinfo%ncubp,ielmax-ielset+1,&
          roptcassemblyinfo%Dpvel,roptcassemblyinfo%DpvelXderiv,&
          roptcassemblyinfo%DpvelYderiv,roptcassemblyinfo%Ddvel,&
          roptcassemblyinfo%DdvelXderiv,roptcassemblyinfo%DdvelYderiv,&
          roptcassemblyinfo%DpvelDofs,roptcassemblyinfo%DdvelDofs,&
          DentryA11,DentryA22,DentryA44,DentryA55,&
          DentryA12,DentryA21,DentryA45,DentryA54,&
          DentryA41,DentryA52,DentryA42,DentryA51,&
          DentryA14,DentryA25,DentryA24,DentryA15,&
          roptcassemblyinfo%DlocalDelta,roptcoperator)      
      
      ! Perform a local matrix-vector multiplication with the local matrices.
      ! Incorporate the computed local matrices into the global matrix.
      call calcDefect (DentryA11,roptcassemblyinfo%Idofs,p_Dx1,p_Dd1,dweight)
      call calcDefect (DentryA12,roptcassemblyinfo%Idofs,p_Dx2,p_Dd1,dweight)
      call calcDefect (DentryA21,roptcassemblyinfo%Idofs,p_Dx1,p_Dd2,dweight)
      call calcDefect (DentryA22,roptcassemblyinfo%Idofs,p_Dx2,p_Dd2,dweight)
      call calcDefect (DentryA44,roptcassemblyinfo%Idofs,p_Dx4,p_Dd4,dweight)
      call calcDefect (DentryA45,roptcassemblyinfo%Idofs,p_Dx5,p_Dd4,dweight)
      call calcDefect (DentryA54,roptcassemblyinfo%Idofs,p_Dx4,p_Dd5,dweight)
      call calcDefect (DentryA55,roptcassemblyinfo%Idofs,p_Dx5,p_Dd5,dweight)
      call calcDefect (DentryA41,roptcassemblyinfo%Idofs,p_Dx1,p_Dd4,dweight)
      call calcDefect (DentryA42,roptcassemblyinfo%Idofs,p_Dx2,p_Dd4,dweight)
      call calcDefect (DentryA51,roptcassemblyinfo%Idofs,p_Dx1,p_Dd5,dweight)
      call calcDefect (DentryA52,roptcassemblyinfo%Idofs,p_Dx2,p_Dd5,dweight)
      call calcDefect (DentryA14,roptcassemblyinfo%Idofs,p_Dx4,p_Dd1,dweight)
      call calcDefect (DentryA24,roptcassemblyinfo%Idofs,p_Dx4,p_Dd2,dweight)
      call calcDefect (DentryA15,roptcassemblyinfo%Idofs,p_Dx5,p_Dd1,dweight)
      call calcDefect (DentryA25,roptcassemblyinfo%Idofs,p_Dx5,p_Dd2,dweight)
      
    end do ! ielset
    !%OMP end do 
    
    ! Release memory
    deallocate(DentryA11)
    deallocate(DentryA22)
    deallocate(DentryA12)
    deallocate(DentryA21)
               
    deallocate(DentryA44)
    deallocate(DentryA55)
    deallocate(DentryA45)
    deallocate(DentryA54)
               
    deallocate(DentryA41)
    deallocate(DentryA52)
    deallocate(DentryA42)
    deallocate(DentryA51)
               
    deallocate(DentryA14)
    deallocate(DentryA25)
    deallocate(DentryA24)
    deallocate(DentryA15)
               
  contains
    
    subroutine calcDefect (DaLocal,KcolLocal,Dx,Dd,dweight)
    
    ! Calculate the local defect:
    !   Dd = Dd - DaLocal*Dx
    
    ! The local matrix to incorporate.
    real(dp), dimension(:,:,:), intent(in) :: DaLocal
    
    ! Global column numbers of the local matrix columns on each element.
    integer, dimension(:,:), intent(in) :: KcolLocal
    
    ! Solution vector
    real(dp), dimension(:), intent(in) :: Dx
    
    ! RHS vector; receives the local defect.
    real(dp), dimension(:), intent(inout) :: Dd
    
    ! Weight for the solution vectot
    real(dp), intent(in) :: dweight
    
      ! local variables
      integer :: iel,idofe,jdofe
      real(dp) :: dval
    
      do iel=1,ubound(DaLocal,3)
        do idofe=1,ubound(DaLocal,1)
          dval = 0.0_DP
          ! DaLocal*Dx. Note that DaLocal is saved transposed for quicker access,
          ! so jdofe in the first index describes the column!
          do jdofe=1,ubound(DaLocal,2)
            dval = dval + DaLocal(jdofe,idofe,iel) * Dx(KcolLocal(jdofe,iel))
          end do
          Dd(KcolLocal(idofe,iel)) = Dd(KcolLocal(idofe,iel)) - dweight*dval
        end do
      end do
    
    end subroutine
    
  end subroutine
     
  ! ----------------------------------------------------------------------

  pure subroutine getLocalDeltaQuad (U1L1,U1L2,U2L1,U2L2,A1L,A2L,IEL,&
                      duMaxR,ddelta,Kvert,Dcorvg,KDFG,IDFL,UPSAM,NUREC)

  ! This routine calculates a local ddelta=DELTA_T for a finite element
  ! T=IEL. This can be used by the streamline diffusion stabilisation
  ! technique as a multiplier of the (local) bilinear form.
  !
  ! The effective velocity that is used for calculating the ddelta
  ! is combined by a weighted mean of the two velocity fields U1,U2
  ! by:
  !                   Ux = A1*U1Lx + A2*U2Lx
  ! The coefficients A1,A2 allow the caller to take influence on which
  ! velocity field to weight more.
  
  ! Main velocity field.
  real(DP), dimension(*), intent(IN) :: U1L1,U1L2
  
  ! Secondary velocity field. 
  real(DP), dimension(*), intent(IN) :: U2L1,U2L2
  
  ! weighting factor for U1L1/U1L2
  real(DP), intent(IN) :: A1L
  
  ! weighting factor for U2L1/U2L2
  real(DP), intent(IN) :: A2L
  
  ! Reciprocal of the maximum norm of velocity in the domain:
  ! 1/duMaxR = 1/||u||_Omega
  real(DP), intent(IN) :: duMaxR
  
  ! Reciprocal value 1/NU of coefficient NU in front of the
  ! Laplacian term of the Navier-Stokes equation
  !   NU * Laplace(u) + u*grad(u) + ...
  real(DP), intent(IN) :: NUREC
  
  ! user defined parameter for configuring the streamline diffusion.
  ! < 0: Simple calculation of ddelta, using 
  !      ddelta = |UPSAM| * h_T.
  ! > 0: usually UPSAM = 0.1 .. 2; Samarskji-like calculation of ddelta using:
  !      ddelta = UPSAM * h_t/||u||_T * 2*Re_T/(1+Re_T)
  real(DP), intent(IN) :: UPSAM
  
  ! Element where the ddelta should be calculated
  integer(PREC_ELEMENTIDX), intent(IN) :: IEL
  
  ! Number of degrees of freedom on element IEL
  integer, intent(IN) :: IDFL
  
  ! Array with global degrees of freedom, corresponding to
  ! local degrees of freedom 1..IDFL on element IEL.
  integer(PREC_DOFIDX), dimension(*), intent(IN) :: KDFG
  
  integer(PREC_VERTEXIDX), dimension(TRIA_MAXNVE2D,*), intent(IN) :: Kvert
  real(DP), dimension(NDIM2D,*), intent(IN) :: Dcorvg

  ! local ddelta
  real(DP), intent(OUT) :: ddelta

  ! local variables
  real(DP) :: dlocalH,DU1,DU2,dunorm,RELOC
  integer(PREC_DOFIDX) :: idof

    ! Loop through the local degrees of freedom on element IEL.
    ! Sum up the velocities on these DOF's. This will result
    ! in the vector (DU1,DU2) representing the (mean) X/Y-velocity
    ! through element IEL.

    ! For elements whose DOF's represent directly the velocity, U1/U2 
    ! represent the mean velocity
    ! along an egde/on the midpoint of each edge, so U1/U2 is
    ! clearly an approximation to the velocity in element T.

    DU1=0.0_DP
    DU2=0.0_DP
    do idof=1,IDFL
      DU1=DU1+(A1L*U1L1(KDFG(idof))+A2L*U2L1(KDFG(idof)))
      DU2=DU2+(A1L*U1L2(KDFG(idof))+A2L*U2L2(KDFG(idof)))
    end do

    ! Calculate the norm of that local velocity:

    dunorm = sqrt(DU1**2+DU2**2) / dble(IDFL)
    
    ! Now we have:   dunorm = ||u||_T
    ! and:           u_T = a1*u1_T + a2*u2_T

    ! If the norm of the velocity is small, we choose ddelta = 0,
    ! which results in central difference in the streamline diffusion
    ! matrix assembling:

    if (dunorm.le.1D-8) then
    
      ddelta = 0.0_DP

    else

      ! u_T defines the "slope" of the velocity through
      ! the element T. At next, calculate the local mesh width
      ! dlocalH = h = h_T on our element T=IEL:

      call getLocalMeshWidthQuad (dlocalH,dunorm, DU1, DU2, IEL, Kvert,Dcorvg)

      ! Calculate ddelta... (cf. p. 121 in Turek's CFD book)

      if (UPSAM.lt.0.0_DP) then

        ! For UPSAM<0, we use simple calculation of ddelta:        
      
        ddelta = abs(UPSAM)*dlocalH
        
      else
      
        ! For UPSAM >= 0, we use standard Samarskji-like calculation
        ! of ddelta. At first calculate the local Reynolds number
        ! RELOC = Re_T = ||u||_T * h_T / NU
        
        RELOC = dunorm*dlocalH*NUREC
        
        ! and then the ddelta = UPSAM * h_t/||u|| * 2*Re_T/(1+Re_T)
        
        ddelta = UPSAM * dlocalH*duMaxR * 2.0_DP*(RELOC/(1.0_DP+RELOC))
        
      end if ! (UPSAM.LT.0.0)
      
    end if ! (dunorm.LE.1D-8)

  end subroutine

  ! ----------------------------------------------------------------------

  pure subroutine getLocalMeshWidthQuad (dlocalH, dunorm,  XBETA1, &
                      XBETA2, JEL,Kvert,Dcorvg)
  
  ! Determine the local mesh width for an element JEL of a 
  ! triangulation.
  
  ! Element where the local h should be calculated
  integer(PREC_ELEMENTIDX), intent(IN)               :: JEL
  
  integer(PREC_VERTEXIDX), dimension(TRIA_MAXNVE2D,*), intent(IN) :: Kvert
  real(DP), dimension(NDIM2D,*), intent(IN)          :: Dcorvg
  
  ! norm ||u||_T = mean velocity through element T=JEL
  real(DP), intent(IN)  :: dunorm
  
  ! mean velocity u_T = (xbeta1,xbeta2) through element T=JEL
  real(DP), intent(IN)  :: XBETA1, XBETA2
  
  ! local mesh width
  real(DP), intent(OUT) :: dlocalH
  
  ! local variables
  real(DP) :: dlambda
  integer(PREC_VERTEXIDX) :: NECK1,NECK2,NECK3,NECK4
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
    
    ! In the next step, we calculate the 'maximum possible mesh with
    ! in direction of the flow'; this is the maximum possible length
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
  real(DP), intent(IN) :: XO,YO
  
  ! Direction of line 1
  real(DP), intent(IN) :: BETA1,BETA2
  
  ! One point on the second line
  real(DP), intent(IN) :: XA,YA
  
  ! Another point on the second line
  real(DP), intent(IN) :: XB,YB
  
  ! Parameter value of the intersection point on line 1.
  ! =0.0, if there is no intersection point
  real(DP), intent(OUT) :: dalpha
  
  real(DP), intent(OUT) :: dlambda
  
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
