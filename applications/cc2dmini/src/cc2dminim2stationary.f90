!##############################################################################
!# ****************************************************************************
!# <name> cc2dminim2stationary </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module invokes the nonlinear solver to solve the basic CC2D problem.
!#
!# The following routines can be found here:
!#
!# 1.) c2d2_getDefect
!#     -> Callback routine. Calculate nonlinear defect
!#
!# 2.) c2d2_getOptimalDamping
!#     -> Auxiliary routine. Calculate optimal damping parameter
!#
!# 3.) c2d2_precondDefect
!#     -> Callback routine. Preconditioning of nonlinear defect
!#
!# 4.) c2d2_getProlRest
!#     -> Auxiliary routine: Set up iterlevel projection structure
!#        with information from INI/DAT files
!#
!# 5.) c2d2_preparePreconditioner
!#     -> Auxiliary routine: Prepare preconditioner of nonlinear iteration
!#
!# 6.) c2d2_releasePreconditioner
!#     -> Auxiliary routine: Clean up preconditioner of nonlinear iteration
!#
!# 7.) c2d2_getNonlinearSolver
!#     -> Auxiliary routine: Initialise nonlinear solver configuration
!#        with information from INI/DAT files
!#
!# 8.) c2d2_solve
!#     -> Invoke the nonlinear solver
!#
!# </purpose>
!##############################################################################

module cc2dminim2stationary

  use fsystem
  use storage
  use linearsolver
  use boundary
  use bilinearformevaluation
  use linearformevaluation
  use cubature
  use basicgeometry
  use matrixfilters
  use vectorfilters
  use bcassembly
  use triangulation
  use spatialdiscretisation
  use coarsegridcorrection
  use spdiscprojection
  use nonlinearsolver
  use paramlist
  use linearsolverautoinitialise
  use matrixrestriction
  use element
  use filtersupport
  use multilevelprojection
  
  use collection
  use convection
    
  use cc2dminim2basic
  use cc2dmini_callback
  
  implicit none
  
!<types>

!<typeblock>

  ! Preconditioner structure for CCxD. This structure saves the configuration of the
  ! preconditioner that is used during the nonlinear iteration.
  
  type t_ccPreconditioner
  
    ! Type of preconditioner.
    ! =0: Preconditioning with inverse mass matrix (not yet implemented),
    ! =1: Preconditioning by linear solver, solving the linearised system,
    ! =2: Preconditioning by Newton-Iteration (not yet implemented),
    integer :: itypePreconditioning
    
    ! Pointer to linear solver node if a linear solver is the preconditioner
    type(t_linsolNode), pointer :: p_rsolverNode
    
    ! A filter chain to filter the vectors and the matrix during the
    ! solution process.
    type(t_filterChain), dimension(1) :: RfilterChain
    
    ! An interlevel projection structure for changing levels
    type(t_interlevelProjectionBlock) :: rprojection

    ! Temporary scalar vector; used for calculating the nonlinear matrix
    ! on lower levels / projecting the solution from higher to lower levels.
    type(t_vectorScalar) :: rtempVectorSc

    ! Temporary scalar vector; used for calculating the optimal damping
    ! parameter.
    type(t_vectorScalar) :: rtempVectorSc2

  end type

!</typeblock>

!</types>
  
contains
  
  ! ***************************************************************************
  !<subroutine>
  
    subroutine c2d2_getDefect (ite,rx,rb,rd,p_rcollection)
  
    use fsystem
    use linearsystemblock
    use collection
    
  !<description>
    ! FOR NONLINEAR ITERATION:
    ! Defect vector calculation callback routine. Based on the current iteration
    ! vector rx and the right hand side vector rb, this routine has to compute the
    ! defect vector rd. The routine accepts a pointer to a collection structure
    ! p_rcollection, which allows the routine to access information from the
    ! main application (e.g. system matrices).
  !</description>

  !<input>
    ! Number of current iteration. 0=build initial defect
    integer, intent(in)                           :: ite

    ! Current iteration vector
    type(t_vectorBlock), intent(in),target        :: rx

    ! Right hand side vector of the equation.
    type(t_vectorBlock), intent(in), target       :: rb
  !</input>
               
  !<inputoutput>
    ! Pointer to collection structure of the application. Points to NULL()
    ! if there is none.
    type(t_collection), pointer                   :: p_rcollection

    ! Defect vector b-A(x)x. This must be filled by the callback routine
    ! with data.
    type(t_vectorBlock), intent(inout), target    :: rd
  !</inputoutput>
  
  !</subroutine>

    ! local variables
    integer :: ilvmax
    type(t_matrixBlock), pointer :: p_rmatrix
    type(t_matrixScalar), pointer :: p_rmatrixLaplace
    type(t_matrixBlock) :: rmatrixLaplaceBlock
    type(t_convUpwind) :: rupwind

    ! A filter chain to pre-filter the vectors and the matrix.
    type(t_filterChain), dimension(1), target :: RfilterChain

      ! Get minimum/maximum level from the collection
      ilvmax = collct_getvalue_int (p_rcollection,'NLMAX')
      
      ! Get the system and the Laplace matrix on the maximum level
      p_rmatrix => collct_getvalue_mat (p_rcollection,'SYSTEMMAT',ilvmax)
      p_rmatrixLaplace => collct_getvalue_matsca (p_rcollection,'LAPLACE',ilvmax)
      
      ! Build a temporary 3x3 block matrix rmatrixLaplace with Laplace
      ! on the main diagonal:
      !
      ! (  L    0   B1 )
      ! (  0    L   B2 )
      ! ( B1^T B2^T 0  )
      !
      !
      call lsysbl_duplicateMatrix (p_rmatrix,rmatrixLaplaceBlock,&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
          
      call lsyssc_duplicateMatrix (p_rmatrixLaplace,&
          rmatrixLaplaceBlock%RmatrixBlock(1,1),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)

      call lsyssc_duplicateMatrix (p_rmatrixLaplace,&
          rmatrixLaplaceBlock%RmatrixBlock(2,2),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      
      ! Now, in the first step, we build the linear part of the nonlinear defect:
      !     d_lin = rhs - (-nu * Laplace(.))*solution
      call lsysbl_copyVector (rb,rd)
      call lsysbl_matVec (rmatrixLaplaceBlock, rx, rd, -1.0_DP, 1.0_DP)
      
      ! Release the temporary matrix again.
      call lsysbl_releaseMatrix (rmatrixLaplaceBlock)
      
      ! For the final defect
      !
      !     d = rhs - (-nu * Laplace(.))*solution - u*grad(.)*solution
      !       = d_lin -  u*grad(.)*solution
      !
      ! in case of the Navier-Stokes equation, we need the nonlinearity.
      
      ! Set up a filter that modifies the block vectors/matrix
      ! according to boundary conditions.
      ! Initialise the first filter of the filter chain as boundary
      ! implementation filter for defect vectors:
      RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

      ! Apply the filter chain to the defect vector.
      ! As the filter consists only of an implementation filter for
      ! boundary conditions, this implements the boundary conditions
      ! into the defect vector.
      call filter_applyFilterChainVec (rd, RfilterChain)
      
      ! Should we discretise the Navier-Stokes nonlinearity?
      if (collct_getvalue_int (p_rcollection,'ISTOKES') .eq. 0) then
      
        ! Which type of stabilisation/strategy for setting up the nonlinearity
        ! do we use?
        select case (collct_getvalue_int (p_rcollection,'IUPWIND'))
        case (1)
      
          ! Set up the upwind structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rupwind%dnu = collct_getvalue_real (p_rcollection,'NU')
          
          ! Set stabilisation parameter
          rupwind%dupsam = collct_getvalue_real (p_rcollection,'UPSAM')
          
          ! Call the upwind method to calculate the nonlinear defect.
          ! As we calculate only the defect, the matrix is ignored!
          call conv_upwind2d (rx, rx, 1.0_DP, 0.0_DP,&
                              rupwind, CONV_MODDEFECT, &
                              p_rmatrix%RmatrixBlock(1,1), rx, rd)
                  
        case DEFAULT
          print *,'Don''t know how to set up nonlinearity!?!'
          stop
          
        end select
        
        ! Apply the filter chain to the defect vector again - since the
        ! implementation of the nonlinearity usually changes the Dirichlet
        ! nodes in the vector!
        call filter_applyFilterChainVec (rd, RfilterChain)
        
      end if

      ! That's it
      
    end subroutine
    
  ! ***************************************************************************

  !<subroutine>

    subroutine c2d2_getOptimalDamping (rd,rx,rb,rtemp1,rtemp2,domega,p_rcollection)
  
  !<description>
    ! This subroutine is called inside of the nonlinear loop, to be precise,
    ! inside of c2d2_precondDefect. It calculates an optiman damping parameter
    ! for the nonlinear defect correction.
    !
    ! The nonlinear loop reads:
    !
    !     $$ u_(n+1) = u_n + OMEGA * C^{-1}d_n $$
    !
    ! with $d_n$ the nonlinear defect and $C^{-1}$ a preconditioner (usually
    ! the linearised system).
    ! Based on the current solution $u_n$, the defect vector $d_n$, the RHS
    ! vector $f_n$ and the previous parameter OMEGA, a new
    ! OMEGA=domega value is calculated.
    !
    ! The nonlinear system matrix on the finest level in the collection is
    ! overwritten by $A(u_n+domega_{old}*C^{-1}d_n)$.
  !</description>

  !<input>
    ! Current iteration vector
    type(t_vectorBlock), intent(in)               :: rx

    ! Current RHS vector of the nonlinear equation
    type(t_vectorBlock), intent(in)               :: rb

    ! Defect vector b-A(x)x.
    type(t_vectorBlock), intent(in)               :: rd

    ! Pointer to collection structure of the application. Points to NULL()
    ! if there is none.
    type(t_collection), pointer                   :: p_rcollection
  !</input>

  !<inputoutput>
    ! A temporary vector in the structure of rx
    type(t_vectorBlock), intent(inout)            :: rtemp1

    ! A 2nd temporary vector in the structure of rx
    type(t_vectorBlock), intent(inout)            :: rtemp2

    ! Damping parameter. On entry: an initial value given e.g. by the
    ! previous step.
    ! On return: The new damping parameter.
    real(DP), intent(inout)                       :: domega
  !</inputoutput>
  
  !</subroutine>

    ! local variables
    integer :: ilvmax
    real(DP) :: domegaMin, domegaMax,dskv1,dskv2
    type(t_matrixBlock), pointer :: p_rmatrix
    type(t_matrixScalar), pointer :: p_rmatrixLaplace
    type(t_convUpwind) :: rupwind

    ! A filter chain to pre-filter the vectors and the matrix.
    type(t_filterChain), dimension(1), target :: RfilterChain

!    DEBUG!!!:
!    real(dp), dimension(:), pointer :: p_vec,p_def,p_temp1,p_temp2,p_da
!    call lsysbl_getbase_double (rd,p_def)
!    call lsysbl_getbase_double (rx,p_vec)
!    call lsysbl_getbase_double (rtemp1,p_temp1)
!    call lsysbl_getbase_double (rtemp2,p_temp2)
!    ilvmax = collct_getvalue_int (p_rcollection,'NLMAX')
!    p_rmatrix => collct_getvalue_mat (p_rcollection,'SYSTEMMAT',ilvmax)
!    call storage_getbase_double (p_rmatrix%RmatrixBlock(1,1)%h_da,p_da)

      ! Get minimum/maximum level from the collection
      ilvmax = collct_getvalue_int (p_rcollection,'NLMAX')
      
      ! Get the minimum/maximum damping parameter from the collection.
      ! They were put there by the main routine using the parameters
      ! from the INI/DAT files.
      domegaMin = collct_getvalue_real (p_rcollection,'OMEGAMIN')
      domegaMax = collct_getvalue_real (p_rcollection,'OMEGAMAX')
      
      ! Is there anything to do?
      if (domegaMin .ge. domegaMax) then
        ! No - cancel.
        domega = domegaMin
        return
      end if

      ! Get the system and the Laplace matrix on the maximum level
      p_rmatrix => collct_getvalue_mat (p_rcollection,'SYSTEMMAT',ilvmax)
      p_rmatrixLaplace => collct_getvalue_matsca (p_rcollection,'LAPLACE',ilvmax)

      ! Set up a filter that modifies the block vectors/matrix
      ! according to boundary conditions.
      ! Initialise the first filter of the filter chain as boundary
      ! implementation filter for defect vectors:
      RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL
      
      ! We now want to calculate a new OMEGA parameter
      ! with OMGMIN < OMEGA < OMGMAX.
      !
      ! The defect correction for a problem like T(u)u=f has the form
      !
      !       u_(n+1)  =  u_n  +  OMEGA * C * ( f - T(u_n)u_n )
      !                =  u_n  +  OMEGA * d_n
      !
      ! with an appropriate preconditioner C, which we don't care here.
      ! In our case, this iteration system can be written as:
      !
      ! (u1)     (u1)                     ( (f1)   [ A         B1] (u1) )
      ! (u2)  := (u2)  + OMEGA * C^{-1} * ( (f2) - [      A    B2] (u2) )
      ! (p )     (p )                     ( (fp)   [ B1^T B2^T 0 ] (p ) )
      !
      !                                   |------------------------------|
      !                                              = d_n
      !                            |-------------------------------------|
      !                                        = Y = (y1,y2,yp) = rd
      !
      ! with KST1=KST1(u1,u2,p) and Y=rd being the solution from
      ! the Oseen equation with
      !
      !                  [ A         B1 ]
      !    C = T(u_n) =  [      A    B2 ]
      !                  [ B1^T B2^T 0  ]
      !
      ! The parameter OMEGA is calculated as the result of the 1D
      ! minimization problem:
      !
      !   OMEGA = min_omega || T(u^l+omega*Y)*(u^l+omega*Y) - f ||_E
      !
      !           < T(u^l+omegaold*Y)Y , f - T(u^l+omegaold*Y)u^l >
      !        ~= -------------------------------------------------
      !              < T(u^l+omegaold*Y)Y , T(u^l+omegaold*Y)Y >
      !
      ! when choosing omegaold=previous omega, which is a good choice
      ! as one can see by linearization (see p. 170, Turek's book).
      !
      ! Here, ||.||_E denotes the the Euclidian norm to the Euclidian
      ! scalar product <.,.>.
      
      ! ==================================================================
      ! First term of scalar product in the nominator
      !
      ! Calculate the new nonlinear block A at the
      ! point rtemp1 = u_n + omegaold*Y
      ! ==================================================================
      !
      ! At first, calculate the point rtemp1 = u_n+omegaold*Y where
      ! to evaluate the matrix.

      call lsysbl_copyVector(rd,rtemp1)
      call lsysbl_vectorLinearComb (rx,rtemp1,1.0_DP,domega)

      ! Should we discretise the Navier-Stokes nonlinearity?
      ! Would mean to rebuild the system matrix. Otherwise we can reuse
      ! our previous diffusion matrix without rebuilding it.
      
      if (collct_getvalue_int (p_rcollection,'ISTOKES') .eq. 0) then
      
        ! Which type of stabilisation/strategy for setting up the nonlinearity
        ! do we use?
        select case (collct_getvalue_int (p_rcollection,'IUPWIND'))
        case (1)
      
          ! Construct the linear part of the nonlinear matrix on the maximum
          ! level.
          !
          ! The system matrix looks like:
          !   (  A    0   B1 )
          !   (  0    A   B2 )
          !   ( B1^T B2^T 0  )
          !
          ! The A-matrix consists of Laplace+Convection.
          ! We build them separately and add together.
          !
          ! So at first, initialise the A-matrix with the Laplace contribution.
          ! We ignore the structure and simply overwrite the content of the
          ! system submatrices with the Laplace matrix.
          call lsyssc_duplicateMatrix (p_rmatrixLaplace,p_rmatrix%RmatrixBlock(1,1),&
                                      LSYSSC_DUP_IGNORE, LSYSSC_DUP_COPY)

          ! Set up the upwind structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rupwind%dnu = collct_getvalue_real (p_rcollection,'NU')
          
          ! Set stabilisation parameter
          rupwind%dupsam = collct_getvalue_real (p_rcollection,'UPSAM')
          
          ! Call the upwind method to evaluate nonlinearity part of the matrix
          ! in the point rtemp1.
          call conv_upwind2d (rtemp1, rtemp1, 1.0_DP, 0.0_DP,&
                              rupwind, CONV_MODMATRIX, &
                              p_rmatrix%RmatrixBlock(1,1))
                                    
        case DEFAULT
          print *,'Don''t know how to set up nonlinearity!?!'
          stop
          
        end select
        
        ! Apply the filter chain to the matrix.
        ! As the filter consists only of an implementation filter for
        ! boundary conditions, this implements the boundary conditions
        ! into the system matrix.
        call filter_applyFilterChainMat (p_rmatrix, RfilterChain)
          
      end if
        
      ! ==================================================================
      ! Second term of the scalar product in the nominator
      ! Calculate the defect rtemp2 = F-T*u_n.
      ! ==================================================================

      call lsysbl_copyVector (rb,rtemp2)
      call lsysbl_matVec (p_rmatrix, rx, rtemp2, -1.0_DP, 1.0_DP)
      
      ! This is a defect vector - filter it! This e.g. implements boundary
      ! conditions.
      call filter_applyFilterChainVec (rtemp2, RfilterChain)
      
      ! ==================================================================
      ! For all terms in the fraction:
      ! Calculate the value  rtemp1 = T*Y
      ! ==================================================================

      call lsysbl_matVec (p_rmatrix, rd, rtemp1, 1.0_DP, 0.0_DP)
      
      ! This is a defect vector against 0 - filter it! This e.g.
      ! implements boundary conditions.
      call filter_applyFilterChainVec (rtemp1, RfilterChain)
      
      ! ==================================================================
      ! Calculation of the fraction terms.
      ! Calculate nominator:    dskv1:= (T*Y,D)   = (rtemp1,rtemp2)
      ! Calculate denominator:  dskv2:= (T*Y,T*Y) = (rtemp1,rtemp1)
      ! ==================================================================
      
      dskv1 = lsysbl_scalarProduct (rtemp1, rtemp2)
      dskv2 = lsysbl_scalarProduct (rtemp1, rtemp1)
      
      if (dskv2 .lt. 1.0E-40_DP) then
        print *,'Error in c2d2_getOptimalDamping. dskv2 nearly zero.'
        print *,'Optimal damping parameter singular.'
        print *,'Is the triangulation ok??? .tri-file destroyed?'
        stop
      end if
      
      ! Ok, we have the nominator and the denominator. Divide them
      ! by each other to calculate the new OMEGA.
      
      domega = dskv1 / dskv2
      
      ! And make sure it's in the allowed range:
      
      domega = max(domegamin,min(domegamax,domega))
      
      ! That's it, we have our new Omega.
  
    end subroutine

  ! ***************************************************************************

  !<subroutine>

    subroutine c2d2_precondDefect (ite,rd,rx,rb,domega,bsuccess,p_rcollection)
  
    use fsystem
    use linearsystemblock
    use collection
    
  !<description>
    ! FOR NONLINEAR ITERATION:
    ! Defect vector calculation callback routine. Based on the current iteration
    ! vector rx and the right hand side vector rb, this routine has to compute the
    ! defect vector rd. The routine accepts a pointer to a collection structure
    ! p_rcollection, which allows the routine to access information from the
    ! main application (e.g. system matrices).
  !</description>

  !<inputoutput>
    ! Number of current iteration.
    integer, intent(in)                           :: ite

    ! Defect vector b-A(x)x. This must be replaced by J^{-1} rd by a preconditioner.
    type(t_vectorBlock), intent(inout), target    :: rd

    ! Pointer to collection structure of the application. Points to NULL()
    ! if there is none.
    type(t_collection), pointer                   :: p_rcollection
    
    ! Damping parameter. Is set to rsolverNode%domega (usually = 1.0_DP)
    ! on the first call to the callback routine.
    ! The callback routine can modify this parameter according to any suitable
    ! algorithm to calculate an 'optimal damping' parameter. The nonlinear loop
    ! will then use this for adding rd to the solution vector:
    ! $$ x_{n+1} = x_n + domega*rd $$
    ! domega will stay at this value until it's changed again.
    real(DP), intent(inout)                       :: domega

    ! If the preconditioning was a success. Is normally automatically set to
    ! TRUE. If there is an error in the preconditioner, this flag can be
    ! set to FALSE. In this case, the nonlinear solver breaks down with
    ! the error flag set to 'preconditioner broke down'.
    logical, intent(inout)                        :: bsuccess
  !</inputoutput>
  
  !<input>
    ! Current iteration vector
    type(t_vectorBlock), intent(in), target       :: rx

    ! Current right hand side of the nonlinear system
    type(t_vectorBlock), intent(in), target       :: rb
  !</input>
  
  !</subroutine>
  
    ! local variables
    integer :: istokes, iupwind,iadaptivematrix
    real(DP) :: dadmatthreshold
  
    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), pointer :: p_rmatrix,p_rmatrixFine
    type(t_matrixScalar), pointer :: p_rmatrixLaplace
    type(t_vectorScalar), pointer :: p_rvectorTemp,p_rvectorTemp2
    type(t_vectorBlock) :: rtemp1,rtemp2
    integer :: ierror,NLMAX,NLMIN, ilev
    type(t_linsolNode), pointer :: p_rsolverNode
    type(t_convUpwind) :: rupwind
    type(t_vectorBlock), pointer :: p_rvectorFine,p_rvectorCoarse

    ! An interlevel projection structure for changing levels
    type(t_interlevelProjectionBlock), pointer :: p_rprojection

    ! A filter chain to pre-filter the vectors and the matrix.
    type(t_filterChain), dimension(1), target :: RfilterChain

!    real(dp), dimension(:), pointer :: p_vec,p_def,p_da
!    call lsysbl_getbase_double (rd,p_def)
!    call lsysbl_getbase_double (rx,p_vec)
!    NLMAX = collct_getvalue_int (p_rcollection,'NLMAX')
!    p_rmatrix => collct_getvalue_mat (p_rcollection,'SYSTEMMAT',NLMAX)
!    call storage_getbase_double (p_rmatrix%RmatrixBlock(1,1)%h_da,p_da)

      ! Get minimum and maximum level from the collection
      NLMAX = collct_getvalue_int (p_rcollection,'NLMAX')
      NLMIN = collct_getvalue_int (p_rcollection,'NLMIN')
      
      ! Get parameters about adaptive matrix generation from collection
      iadaptivematrix = collct_getvalue_int (p_rcollection,'IADAPTIVEMATRIX')
      dadmatthreshold = collct_getvalue_real (p_rcollection,'dAdMatThreshold')
      
      ! Set up a filter that modifies the block vectors/matrix
      ! according to boundary conditions.
      ! Initialise the first filter of the filter chain as boundary
      ! implementation filter for defect vectors:
      RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL
      
      ! Get the interlevel projection structure and the temporary vector
      ! from the collection.
      ! Our 'parent' prepared there how to interpolate the solution on the
      ! fine grid to coarser grids.
      p_rprojection => collct_getvalue_ilvp(p_rcollection,'ILVPROJECTION')
      p_rvectorTemp => collct_getvalue_vecsca(p_rcollection,'RTEMPSCALAR')

      ! If we discretise Navier-Stokes, we have to set up a stabilisation
      ! structure dor upwind / streamiline diffusion / ...
      
      istokes = collct_getvalue_int (p_rcollection,'ISTOKES')
      
      if (istokes .eq. 0) then
      
        iupwind = collct_getvalue_int (p_rcollection,'IUPWIND')
      
        select case (iupwind)
        case (1)
          ! Set up the upwind structure for the creation of the defect.
          ! There's not much to do, only initialise the viscosity...
          rupwind%dnu = collct_getvalue_real (p_rcollection,'NU')
          
          ! Set stabilisation parameter
          rupwind%dupsam = collct_getvalue_real (p_rcollection,'UPSAM')
          
        case DEFAULT
          print *,'Don''t know how to set up nonlinearity!?!'
          stop
        
        end select
      end if
      
      ! On all levels, we have to set up the nonlinear system matrix,
      ! so that the linear solver can be applied to it.

      nullify(p_rmatrix)
      do ilev=NLMAX,NLMIN,-1
      
        ! Get the system matrix and the Laplace matrix
        p_rmatrixFine => p_rmatrix
        p_rmatrix => collct_getvalue_mat (p_rcollection,'SYSTEMMAT',ilev)
        p_rmatrixLaplace => collct_getvalue_matsca (p_rcollection,'LAPLACE',ilev)
        
        ! On the highest level, we use rx as solution to build the nonlinear
        ! matrix. On lower levels, we have to create a solution
        ! on that level from a fine-grid solution before we can use
        ! it to build the matrix!
        if (ilev .eq. NLMAX) then
          p_rvectorCoarse => rx
        else
          ! Get the temporary vector on level i. Will receive the solution
          ! vector on that level.
          p_rvectorCoarse => collct_getvalue_vec (p_rcollection,'RTEMPVEC',ilev)
          
          ! Get the solution vector on level i+1. This is either the temporary
          ! vector on that level, or the solution vector on the maximum level.
          if (ilev .lt. NLMAX-1) then
            p_rvectorFine => collct_getvalue_vec (p_rcollection,'RTEMPVEC',ilev+1)
          else
            p_rvectorFine => rx
          end if

          ! Interpolate the solution from the finer grid to the coarser grid.
          ! The interpolation is configured in the interlevel projection
          ! structure we got from the collection.
          call mlprj_performInterpolation (p_rprojection,p_rvectorCoarse, &
                                           p_rvectorFine,p_rvectorTemp)

          ! Apply the filter chain to the temp vector.
          ! This implements the boundary conditions that are attached to it.
          call filter_applyFilterChainVec (p_rvectorCoarse, RfilterChain)

        end if
        
        ! Should we discretise the Navier-Stokes nonlinearity?
        if (istokes .eq. 0) then
        
          ! Which type of stabilisation/strategy for setting up the nonlinearity
          ! do we use?
          select case (iupwind)
          case (1)
        
            ! The system matrix looks like:
            !   (  A    0   B1 )
            !   (  0    A   B2 )
            !   ( B1^T B2^T 0  )
            !
            ! The A-matrix consists of Laplace+Convection.
            ! We build them separately and add together.
            !
            ! So at first, initialise the A-matrix with the Laplace contribution.
            ! We ignore the structure and simply overwrite the content of the
            ! system submatrices with the Laplace matrix.
            call lsyssc_duplicateMatrix (p_rmatrixLaplace,p_rmatrix%RmatrixBlock(1,1),&
                                        LSYSSC_DUP_IGNORE, LSYSSC_DUP_COPY)

            ! Call the upwind method to calculate the nonlinear matrix.
            call conv_upwind2d (p_rvectorCoarse, p_rvectorCoarse, 1.0_DP, 0.0_DP,&
                                rupwind, CONV_MODMATRIX, &
                                p_rmatrix%RmatrixBlock(1,1))
                                         
          case DEFAULT
            print *,'Don''t know how to set up nonlinearity!?!'
            stop
            
          end select
        
        else
          ! The system matrix looks like:
          !   (  A    0   B1 )
          !   (  0    A   B2 )
          !   ( B1^T B2^T 0  )
          !
          ! The A-matrix is a simple Laplace.
          !
          ! Copy the Laplace matrix to A and filter it according to the boundary
          ! conditions - this gives then the system matrix.
          call lsyssc_duplicateMatrix (p_rmatrixLaplace,p_rmatrix%RmatrixBlock(1,1),&
                                      LSYSSC_DUP_IGNORE, LSYSSC_DUP_COPY)

        end if

        ! For the construction of matrices on lower levels, call the matrix
        ! restriction. In case we have a uniform discretisation with Q1~,
        ! iadaptivematrix is <> 0 and so this will rebuild some matrix entries
        ! by a Galerkin approach using constant prolongation/restriction.
        ! This helps to stabilise the solver if there are elements in the
        ! mesh with high aspect ratio.
        if (ilev .lt. NLMAX) then
          call mrest_matrixRestrictionEX3Y (p_rmatrixFine%RmatrixBlock(1,1),&
                                            p_rmatrix%RmatrixBlock(1,1),&
                                            iadaptivematrix, dadmatthreshold)
        end if
      
        ! Apply the filter chain to the matrix.
        ! As the filter consists only of an implementation filter for
        ! boundary conditions, this implements the boundary conditions
        ! into the system matrix.
        call filter_applyFilterChainMat (p_rmatrix, RfilterChain)
        
      end do
      
      ! Our 'parent' (the caller of the nonlinear solver) has prepared
      ! a preconditioner node for us (a linear solver with symbolically
      ! factorised matrices). Get this from the collection.
      
      p_rsolverNode => collct_getvalue_linsol(p_rcollection,'LINSOLVER')

      ! Initialise data of the solver. This in fact performs a numeric
      ! factorisation of the matrices in UMFPACK-like solvers.
      call linsol_initData (p_rsolverNode, ierror)
      if (ierror .ne. LINSOL_ERR_NOERROR) stop
      
      
      ! Finally solve the system. As we want to solve Ax=b with
      ! b being the real RHS and x being the real solution vector,
      ! we use linsol_solveAdaptively. If b is a defect
      ! RHS and x a defect update to be added to a solution vector,
      ! we would have to use linsol_precondDefect instead.
      call linsol_precondDefect (p_rsolverNode,rd)

      ! Release the numeric factorisation of the matrix.
      ! We don't release the symbolic factorisation, as we can use them
      ! for the next iteration.
      call linsol_doneData (p_rsolverNode)

      ! Finally calculate a new damping parameter domega.
      !
      ! For this purpose, we need two temporary vectors.
      ! On one hand, we have p_rvectorTemp.
      ! Get the second temporary vector from the collection as it was
      ! prepared by our 'parent' that invoked the nonlinear solver.
      p_rvectorTemp2 => collct_getvalue_vecsca(p_rcollection,'RTEMP2SCALAR')
      
      ! Both temp vectors are scalar, but we need block-vectors in the
      ! structure of rx/rb. Derive block vectors in that structure that
      ! share their memory with the scalar temp vectors. Note that the
      ! temp vectors are created large enough by our parent!
      call lsysbl_createVecFromScalar (p_rvectorTemp,rtemp1)
      call lsysbl_enforceStructure (rb,rtemp1)

      call lsysbl_createVecFromScalar (p_rvectorTemp2,rtemp2)
      call lsysbl_enforceStructure (rb,rtemp2)

      ! Calculate the omega
      call c2d2_getOptimalDamping (rd,rx,rb,rtemp1,rtemp2,&
                                   domega,p_rcollection)

      ! Release the temp block vectors. This only cleans up the structure.
      ! The data is not released from heap as it belongs to the
      ! scalar temp vectors.
      call lsysbl_releaseVector (rtemp2)
      call lsysbl_releaseVector (rtemp1)

    end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine c2d2_getProlRest (rprojection, rparamList, sname)
  
!<description>
  ! Initialises an existing interlevel projection structure rprojection
  ! with parameters from the INI/DAT files. sname is the section in the
  ! parameter list containing parameters about prolongation restriction.
!</description>

!<input>
  ! Parameter list that contains the parameters from the INI/DAT file(s).
  type(t_parlist), intent(in) :: rparamList
  
  ! Name of the section in the parameter list containing the parameters
  ! of the prolongation/restriction.
  character(LEN=*), intent(in) :: sname
!</input>

!<output>
  ! An interlevel projection block structure containing an initial
  ! configuration of prolongation/restriction. The structure is modified
  ! according to the parameters in the INI/DAT file(s).
  type(t_interlevelProjectionBlock), intent(inout) :: rprojection
!</output>

!</subroutine>

    ! local variables
    type(t_parlstSection), pointer :: p_rsection
    integer :: i1
    real(DP) :: d1

    ! Check that there is a section called sname - otherwise we
    ! cannot create anything!
    
    call parlst_querysection(rparamList, sname, p_rsection)

    if (.not. associated(p_rsection)) then
      ! We use the default configuration; stop here.
      return
    end if
    
    ! Now take a look which parameters appear in that section.

    ! Prolongation/restriction order for velocity components
    call parlst_getvalue_int (p_rsection,'iinterpolationOrderVel',i1,-1)
    
    if (i1 .ne. -1) then
      ! Initialise order of prolongation/restriction for velocity components
      rprojection%RscalarProjection(:,1:NDIM2D)%iprolongationOrder  = i1
      rprojection%RscalarProjection(:,1:NDIM2D)%irestrictionOrder   = i1
      rprojection%RscalarProjection(:,1:NDIM2D)%iinterpolationOrder = i1
    end if

    ! Prolongation/restriction order for pressure
    call parlst_getvalue_int (p_rsection,'iinterpolationOrderPress',i1,-1)
    
    if (i1 .ne. -1) then
      ! Initialise order of prolongation/restriction for velocity components
      rprojection%RscalarProjection(:,NDIM2D+1)%iprolongationOrder  = i1
      rprojection%RscalarProjection(:,NDIM2D+1)%irestrictionOrder   = i1
      rprojection%RscalarProjection(:,NDIM2D+1)%iinterpolationOrder = i1
    end if
    
    ! Prolongation/restriction variant for velocity components
    ! in case of Q1~ discretisation
    call parlst_getvalue_int (p_rsection,'iinterpolationVariantVel',i1,0)
    
    if (i1 .ne. -1) then
      rprojection%RscalarProjection(:,1:NDIM2D)%iprolVariant  = i1
      rprojection%RscalarProjection(:,1:NDIM2D)%irestVariant  = i1
    end if
    
    ! Aspect-ratio indicator in case of Q1~ discretisation
    ! with extended prolongation/restriction
    call parlst_getvalue_int (p_rsection,'iintARIndicatorEX3YVel',i1,1)
    
    if (i1 .ne. 1) then
      rprojection%RscalarProjection(:,1:NDIM2D)%iprolARIndicatorEX3Y  = i1
      rprojection%RscalarProjection(:,1:NDIM2D)%irestARIndicatorEX3Y  = i1
    end if

    ! Aspect-ratio bound for switching to constant prolongation/restriction
    ! in case of Q1~ discretisation with extended prolongation/restriction
    call parlst_getvalue_double (p_rsection,'dintARboundEX3YVel',d1,20.0_DP)
    
    if (d1 .ne. 20.0_DP) then
      rprojection%RscalarProjection(:,1:NDIM2D)%dprolARboundEX3Y  = d1
      rprojection%RscalarProjection(:,1:NDIM2D)%drestARboundEX3Y  = d1
    end if

  end subroutine
    
  ! ***************************************************************************

!<subroutine>

  subroutine c2d2_preparePreconditioner (rproblem,rpreconditioner)
  
!<description>
  ! This routine prepares the preconditioner that us used during the
  ! nonlinear iteration. The structure rpreconditioner will be initialised.
  ! Necessary variables will be added to the collection structure in
  ! rproblem\%rcollection to be available in the callback routines.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!<output>
  ! A preconditioner structure for the CCxD problem. Will be initialised
  ! with data.
  type(t_ccPreconditioner), intent(out), target :: rpreconditioner
!</output>

!</subroutine>

    ! local variables
    integer :: NLMIN,NLMAX
    integer :: i
    real(DP) :: d
    integer :: imaxmem
    character(LEN=PARLST_MLDATA) :: ssolverName,sstring
    type(t_spatialDiscretisation), pointer :: p_rdiscr
    integer(I32) :: celement

    ! Error indicator during initialisation of the solver
    integer :: ierror
  
    ! A pointer to the system matrix and the RHS vector as well as
    ! the discretisation
    type(t_matrixBlock), pointer :: p_rmatrix
    type(t_vectorBlock), pointer :: p_rrhs,p_rvector

    ! An array for the system matrix(matrices) during the initialisation of
    ! the linear solver.
    type(t_matrixBlock), dimension(NNLEV) :: Rmatrices
    
    ! At first, ask the parameters in the INI/DAT file which type of
    ! preconditioner is to be used. The data in the preconditioner structure
    ! is to be initialised appropriately!
    call parlst_getvalue_int (rproblem%rparamList, 'CC2D-NONLINEAR', &
                              'itypePreconditioning', &
                              rpreconditioner%itypePreconditioning, 1)
    
    select case (rpreconditioner%itypePreconditioning)
    case (1)
      ! Ok, we have to initialise a linear solver for solving the linearised
      ! problem.
      
      ! Which levels have we to take care of during the solution process?
      NLMIN = rproblem%NLMIN
      NLMAX = rproblem%NLMAX
    
      ! Get our right hand side / solution / matrix on the finest
      ! level from the problem structure.
      p_rrhs    => rproblem%rrhs
      p_rvector => rproblem%rvector
      p_rmatrix => rproblem%RlevelInfo(NLMAX)%rmatrix
      
      ! During the linear solver, the boundary conditions must
      ! frequently be imposed to the vectors. This is done using
      ! a filter chain. As the linear solver does not work with
      ! the actual solution vectors but with defect vectors instead,
      ! a filter for implementing the real boundary conditions
      ! would be wrong.
      ! Therefore, create a filter chain with one filter only,
      ! which implements Dirichlet-conditions into a defect vector.
      rpreconditioner%RfilterChain(1)%ifilterType = FILTER_DISCBCDEFREAL

      ! Figure out the name of the section that contains the information
      ! about the linear subsolver. Ask the parameter list from the INI/DAT file
      ! for the 'slinearSolver' value
      call parlst_getvalue_string (rproblem%rparamList, 'CC2D-NONLINEAR', &
                                  'slinearSolver', sstring, '')
      ssolverName = ''
      if (sstring .ne. '') read (sstring,*) ssolverName
      if (ssolverName .eq. '') then
        print *,'No linear subsolver!'
        stop
      end if
                                    
      ! Initialise a standard interlevel projection structure. We
      ! can use the same structure for all levels. Therefore it's enough
      ! to initialise one structure using the RHS vector on the finest
      ! level to specify the shape of the PDE-discretisation.
      call mlprj_initProjectionVec (rpreconditioner%rprojection,rproblem%rrhs)
      
      ! Initialise the projection structure with data from the INI/DAT
      ! files. This allows to configure prolongation/restriction.
      call c2d2_getProlRest (rpreconditioner%rprojection, rproblem%rparamList, &
                             'CC-PROLREST')
      
      ! Add the interlevel projection structure to the collection; we can
      ! use it later for setting up nonlinear matrices.
      call collct_setvalue_ilvp(rproblem%rcollection,'ILVPROJECTION',&
                                rpreconditioner%rprojection,.true.)
      
      ! Initialise the linear subsolver using the parameters from the INI/DAT
      ! files, the prepared filter chain and the interlevel projection structure.
      ! This gives us the linear solver node rpreconditioner%p_rsolverNode
      ! which identifies the linear solver.
      call linsolinit_initFromFile (rpreconditioner%p_rsolverNode,&
                                    rproblem%rparamList,ssolverName,&
                                    NLMAX-NLMIN+1,rpreconditioner%RfilterChain,&
                                    rpreconditioner%rprojection)
      
      ! How much memory is necessary for performing the level change?
      ! We ourself must build nonlinear matrices on multiple levels and have
      ! to interpolate the solution vector from finer level to coarser ones.
      ! We need temporary memory for this purpose...

      imaxmem = 0
      do i=NLMIN+1,NLMAX
        ! Pass the system metrices on the coarse/fine grid to
        ! mlprj_getTempMemoryMat to specify the discretisation structures
        ! of all equations in the PDE there.
        imaxmem = max(imaxmem,mlprj_getTempMemoryMat (rpreconditioner%rprojection,&
                              rproblem%RlevelInfo(i-1)%rmatrix,&
                              rproblem%RlevelInfo(i)%rmatrix))
      end do
      
      ! Set up a scalar temporary vector that we need for building up nonlinear
      ! matrices. It must be at least as large as MAXMEM and NEQ(finest level),
      ! as we use it for resorting vectors, too.
      call lsyssc_createVector (rpreconditioner%rtempVectorSc,&
                                max(imaxmem,rproblem%rrhs%NEQ),.false.)
      call collct_setvalue_vecsca(rproblem%rcollection,'RTEMPSCALAR',&
                                  rpreconditioner%rtempVectorSc,.true.)
      
      ! Set up a second temporary vector that we need for calculating
      ! the optimal defect correction.
      call lsyssc_createVector (rpreconditioner%rtempVectorSc2,&
                                rproblem%rrhs%NEQ,.false.,ST_DOUBLE)
      call collct_setvalue_vecsca(rproblem%rcollection,'RTEMP2SCALAR',&
                                  rpreconditioner%rtempVectorSc2,.true.)
      
      ! Attach the system matrices to the solver.
      !
      ! We copy our matrices to a big matrix array and transfer that
      ! to the setMatrices routines. This intitialises then the matrices
      ! on all levels according to that array.
      Rmatrices(NLMIN:NLMAX) = rproblem%RlevelInfo(NLMIN:NLMAX)%rmatrix
      call linsol_setMatrices(rpreconditioner%p_rsolverNode,Rmatrices(NLMIN:NLMAX))
      
      ! Initialise structure/data of the solver. This allows the
      ! solver to allocate memory / perform some precalculation
      ! to the problem.
      call linsol_initStructure (rpreconditioner%p_rsolverNode,ierror)
      if (ierror .ne. LINSOL_ERR_NOERROR) stop
      
      ! Put the prepared solver node to the collection for later use.
      ! Remember, it's our preconditioner we need during the nonlinear
      ! iteration!
      call collct_setvalue_linsol(rproblem%rcollection,'LINSOLVER',&
                                  rpreconditioner%p_rsolverNode,.true.)
      
      ! Add information about adaptive matrix generation from INI/DAT files
      ! to the collection.
      call parlst_getvalue_int (rproblem%rparamList, 'CC-DISCRETISATION', &
                                'iAdaptiveMatrix', i, 0)
                                
      ! Switch off adaptive matrix generation if our discretisation is not a
      ! uniform Q1~ discretisation - the matrix restriction does not support
      ! other cases.
      p_rdiscr => rproblem%RlevelInfo(NLMAX)%p_rdiscretisation%RspatialDiscr(1)
      
      call spdiscr_getElemGroupInfo (p_rdiscr,1,celement)
      
      if ((p_rdiscr%ccomplexity .ne. SPDISC_UNIFORM) .or. &
          ((celement .ne. EL_E030) .and. (celement .ne. EL_E031) .and. &
           (celement .ne. EL_EM30) .and. (celement .ne. EL_EM31))) then
        i = 0
      end if
      
      call collct_setvalue_int(rproblem%rcollection,'IADAPTIVEMATRIX',i,.true.)
                                
      call parlst_getvalue_double(rproblem%rparamList, 'CC-DISCRETISATION', &
                                 'dAdMatThreshold', d, 20.0_DP)
      call collct_setvalue_real(rproblem%rcollection,'DADMATTHRESHOLD',d,.true.)
      
    case DEFAULT
      
      ! Unknown preconditioner
      print *,'Unknown preconditioner for nonlinear iteration!'
      stop
      
    end select

  end subroutine

! ***************************************************************************

!<subroutine>

  subroutine c2d2_releasePreconditioner (rproblem,rpreconditioner)
  
!<description>
  ! This routine releases the preconditioner for the nonlinear iteration
  ! which was prepared in c2d2_preparePreconditioner. Memory is released
  ! from heap. All variables that were attached to the collection
  ! structure rproblem\%rcollection for use in the nonlinear iteration
  ! are removed from the collection.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem

  ! The preconditioner structure for the CCxD problem.
  type(t_ccPreconditioner), intent(inout), target :: rpreconditioner
!</inputoutput>

!</subroutine>

    ! Which preconditioner do we have?
    select case (rpreconditioner%itypePreconditioning)
    case (1)
      ! Preconditioner was a linear solver structure.
      !
      ! Release the temporary vector(s)
      call lsyssc_releaseVector (rpreconditioner%rtempVectorSc)
      call lsyssc_releaseVector (rpreconditioner%rtempVectorSc2)
      
      ! Remove the solver node from the collection - not needed anymore there
      call collct_deletevalue(rproblem%rcollection,'LINSOLVER')
      
      ! Remove the temporary vector from the collection
      call collct_deletevalue(rproblem%rcollection,'RTEMPSCALAR')
      call collct_deletevalue(rproblem%rcollection,'RTEMP2SCALAR')
      
      ! Remove the interlevel projection structure
      call collct_deletevalue(rproblem%rcollection,'ILVPROJECTION')
      
      ! Remove parameters for adaptive matrix generation
      call collct_deletevalue(rproblem%rcollection,'DADMATTHRESHOLD')
      call collct_deletevalue(rproblem%rcollection,'IADAPTIVEMATRIX')

      ! Clean up the linear solver, release all memory, remove the solver node
      ! from memory.
      call linsol_releaseSolver (rpreconditioner%p_rsolverNode)
      
      ! Release the ultilevel projection structure.
      call mlprj_doneProjection (rpreconditioner%rprojection)
      
    case DEFAULT
      
      ! Unknown preconditioner
      print *,'Unknown preconditioner for nonlinear iteration!'
      stop
      
    end select

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine c2d2_getNonlinearSolver (rnlSolver, rparamList, sname)
  
!<description>
  ! Creates a nonlinear solver node rnlSolver and initialises it with parameters
  ! from the INI/DAT files given in the parameter list rparamList.
  ! sname is the name of a section in rparamList that configures the
  ! parameter of the nonlinear solver.
!</description>

!<input>
  ! Parameter list that contains the parameters from the INI/DAT file(s).
  type(t_parlist), intent(in) :: rparamList
  
  ! Name of the section in the parameter list containing the parameters
  ! of the nonlinear solver.
  character(LEN=*), intent(in) :: sname
!</input>

!<output>
  ! A t_nlsolNode structure that contains the configuration of the nonlinear
  ! solver. The parameters are initialised according to the information
  ! in the section sname of the parameter list rparamList
  type(t_nlsolNode) :: rnlSolver
!</output>

!</subroutine>

    ! local variables
    type(t_parlstSection), pointer :: p_rsection

    ! Check that there is a section called sname - otherwise we
    ! cannot create anything!
    
    call parlst_querysection(rparamList, sname, p_rsection)

    if (.not. associated(p_rsection)) then
      print *,'Cannot create nonlinear solver; no section '''&
              //trim(sname)//'''!'
      stop
    end if
    
    ! Parse the given parameters now to initialise the solver node.
    ! This is now CCxD-specific!
    
    call parlst_getvalue_int (p_rsection, 'nminIterations', &
                              rnlSolver%nminIterations, rnlSolver%nminIterations)

    call parlst_getvalue_int (p_rsection, 'nmaxIterations', &
                              rnlSolver%nmaxIterations, rnlSolver%nmaxIterations)

    call parlst_getvalue_int (p_rsection, 'ioutputLevel', &
                              rnlSolver%ioutputLevel, rnlSolver%ioutputLevel)

    call parlst_getvalue_double (p_rsection, 'depsUR', &
                                 rnlSolver%DepsRel(1), rnlSolver%DepsRel(1))
    rnlSolver%DepsRel(2) = rnlSolver%DepsRel(1)

    call parlst_getvalue_double (p_rsection, 'depsPR', &
                                 rnlSolver%DepsRel(3), rnlSolver%DepsRel(3))

    call parlst_getvalue_double (p_rsection, 'depsD', &
                                 rnlSolver%DepsAbs(1), rnlSolver%DepsAbs(1))
    rnlSolver%DepsAbs(2) = rnlSolver%DepsAbs(1)

    call parlst_getvalue_double (p_rsection, 'depsDiv', &
                                 rnlSolver%DepsAbs(3), rnlSolver%DepsAbs(3))

    ! Initial damping parameter.
    call parlst_getvalue_double (p_rsection, 'domegaIni', &
                                 rnlSolver%domega, rnlSolver%domega)

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine c2d2_solve (rproblem)
  
!<description>
  ! Solves the given problem by applying a nonlinear solver with a preconditioner
  ! configured in the INI/DAT files.
!</description>

!<inputoutput>
  ! A problem structure saving problem-dependent information.
  type(t_problem), intent(inout), target :: rproblem
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_ccPreconditioner) :: rpreconditioner
    type(t_vectorBlock), target :: rtempBlock
    real(DP) :: d1

    ! The nonlinear solver configuration
    type(t_nlsolNode) :: rnlSol
    
    ! Initialise the preconditioner for the nonlinear iteration
    call c2d2_preparePreconditioner (rproblem,rpreconditioner)
    
    ! Add the minimum/maximum damping parameter to the collection
    call parlst_getvalue_double (rproblem%rparamList, 'CC2D-NONLINEAR', &
                                 'domegaMin', d1, 0.0_DP)
    call collct_setvalue_real(rproblem%rcollection,'OMEGAMIN',d1,.true.)
    
    call parlst_getvalue_double (rproblem%rparamList, 'CC2D-NONLINEAR', &
                                 'domegaMax', d1, 0.0_DP)
    call collct_setvalue_real(rproblem%rcollection,'OMEGAMAX',d1,.true.)
    
    ! Create a temporary vector we need for the nonliner iteration.
    call lsysbl_createVecBlockIndirect (rproblem%rrhs, rtempBlock, .false.)

    ! Initialise the nonlinear solver node rnlSol with parameters from
    ! the INI/DAT files.
    call c2d2_getNonlinearSolver (rnlSol, rproblem%rparamList, 'CC2D-NONLINEAR')

    ! Call the nonlinear solver. For preconditioning and defect calculation,
    ! the solver calls our callback routines.
    call nlsol_performSolve(rnlSol,rproblem%rvector,rproblem%rrhs,rtempBlock,&
                            c2d2_getDefect,c2d2_precondDefect,&
                            rcollection=rproblem%rcollection)
                            
    ! Delete min./max. damping parameters from the collection
    call collct_deletevalue (rproblem%rcollection,'OMEGAMAX')
    call collct_deletevalue (rproblem%rcollection,'OMEGAMIN')

    ! Release the temporary vector
    call lsysbl_releaseVector (rtempBlock)
    
    ! Release the preconditioner
    call c2d2_releasePreconditioner (rproblem,rpreconditioner)
    
    print *
    print *,'Nonlinear solver statistics'
    print *,'---------------------------'
    print *,'Initial defect: ',rnlSol%DinitialDefect(1)
    print *,'Final defect:  ',rnlSol%DfinalDefect(1)
    print *,'#Iterations:   ',rnlSol%iiterations
    
  end subroutine

end module
