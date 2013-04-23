!##############################################################################
!# ****************************************************************************
!# <name> blockmatassemblystdop </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module defines a set of standard operators that can be used
!# with the block matrix/vector assembly. These routines may serve as
!# a template for more complex callback routines.
!#
!# The following set of callback routines realise standard operators.
!# The "bma_fcalc_XXXX" methods can directly be used with bma_buildMatrix, 
!# bma_buildVector and bma_buildIntegral.
!# The "bma_docalc_XXXX" methods can be called in a callback routine for 
!# bma_buildMatrix, bma_buildVector and bma_buildIntegral.
!#
!# 1.) bma_fcalc_mass / bma_fcalc_massDiag / bma_docalc_mass
!#     -> Used with bma_buildMatrix, this calculates mass matrices 
!#        at a specified position / in all diagonal blocks of a block matrix.
!#
!# 2.) bma_fcalc_laplace / bma_fcalc_laplaceDiag / bma_docalc_laplace
!#     -> Used with bma_buildMatrix, this calculates Laplace matrices
!#        at a specified position / in all diagonal blocks of a block matrix.
!#
!# 3.) bma_fcalc_gradientdiv / bma_docalc_gradientdiv
!#     -> Used with bma_buildMatrix, this calculates a gradient matrix
!#        at a specified position / in all diagonal blocks of a block matrix.
!#        The gradient matrix is set up after a partial integration, ( -p, div u ).
!#
!# 4.) bma_fcalc_divergence / bma_docalc_divergence
!#     -> Used with bma_buildMatrix, this calculates a divergence matrix
!#        at a specified position / in all diagonal blocks of a block matrix.
!#
!# 5.) bma_fcalc_convection_ugradvw / bma_docalc_convection_ugradvw
!#     -> Used with bma_buildMatrix, this routine can be used to assemble
!#        the convection operator ( (u grad) v, w) at a specified position in a block
!#        matrix.
!#
!# 6.) bma_fcalc_convection_vugradw / bma_docalc_convection_vugradw
!#     -> Used with bma_buildMatrix, this routine can be used to assemble
!#        the convection operator ( v, (u grad) w) at a specified position in a block
!#        matrix.
!#
!# 7.) bma_fcalc_convection_graduvw / bma_docalc_convection_graduvw
!#     -> Used with bma_buildMatrix, this routine can be used to assemble
!#        the convection operator ( (grad u) v, w) at a specified position in a block
!#        matrix.
!#
!# 8.) bma_fcalc_convection_vgraduw / bma_docalc_convection_vgraduw
!#     -> Used with bma_buildMatrix, this routine can be used to assemble
!#        the convection operator ( v, (grad u) w) at a specified position in a block
!#        matrix.
!#
!# 9.) bma_fcalc_rhsConst / bma_docalc_rhsConst
!#     -> Calculates the RHS vector based on the function f=const(=1).
!#
!# 10.) bma_fcalc_rhsBubble / bma_docalc_rhsBubble
!#      -> Calculates the RHS vector based on the function f=32*y*(1-y)+32*x*(1-x)
!#         which is the RHS for u=16*x*(1-x)*y*(1-y) in the Laplace
!#         equation -Laplace(u)=f.
!#
!# 11.) bma_fcalc_rhsBubblePlusFE / bma_docalc_rhsBubblePlusFE
!#      -> Calculates the RHS vector based on the function 
!#         f=32*y*(1-y)+32*x*(1-x) + v(x,y)
!#         with v(x,y) being a finite element function passed via parameters.
!#
!# 12.) bma_fcalc_rhsFE / bma_docalc_rhsFE
!#     -> Calculates the RHS vector based on the function 
!#        f=v(x,y)
!#        with v(x,y) being a finite element function passed via parameters.
!#
!# 13.) bma_fcalc_integralOne / bma_docalc_integralOne
!#     -> Calculates the integral of the function v=1 (which results in the
!#        size of the domain).
!#
!# 14.) bma_fcalc_integralFE / bma_fcalc_integralFE
!#     -> Calculates the integral of an arbitrary FEM function.
!#
!# 15.) bma_fcalc_bubbleL2error / bma_docalc_bubbleL2error
!#      -> Calculates the squared L2 error of a FEM function to a 
!#         bubble function
!#        
!# 16.) bma_fcalc_bubbleH1error / bma_docalc_bubbleH1error
!#     -> Calculates the squared H1 error of a FEM function to a 
!#        bubble function
!#
!# 17.) bma_fcalc_L2norm / bma_docalc_L2norm
!#     -> Calculates the squared L2 norm of a FEM function
!#        
!# 18.) bma_fcalc_H1norm / bma_docalc_H1norm
!#     -> Calculates the squared H1 (semi-)norm of a FEM function
!#
!# 19.) bma_fcalc_divergenceL2norm / bma_docalc_divergenceL2norm
!#     -> Calculates the squared L2-norm of the divergence of a vector field.
!# </purpose>
!##############################################################################

module blockmatassemblystdop

  use fsystem
  use storage
  use genoutput
  use geometryaux
  use collection, only: t_collection
  use basicgeometry
  use perfconfig

  use feevaluation2
  use blockmatassemblybase

  implicit none

  private

  !************************************************************************

  ! global performance configuration
  type(t_perfconfig), target, save :: bma_perfconfig

  !************************************************************************

  public :: bma_docalc_mass
  public :: bma_docalc_laplace
  public :: bma_docalc_gradientdiv
  public :: bma_docalc_divergence
  public :: bma_docalc_convection_ugradvw
  public :: bma_docalc_convection_vugradw
  public :: bma_docalc_convection_graduvw
  public :: bma_docalc_convection_vgraduw

  public :: bma_docalc_rhsConst
  public :: bma_docalc_rhsBubble
  public :: bma_docalc_rhsBubblePlusFE
  public :: bma_docalc_rhsFE

  public :: bma_fcalc_massDiag
  public :: bma_fcalc_laplaceDiag
  public :: bma_fcalc_mass
  public :: bma_fcalc_laplace
  public :: bma_fcalc_gradientdiv
  public :: bma_fcalc_divergence
  public :: bma_fcalc_convection_ugradvw
  public :: bma_fcalc_convection_vugradw
  public :: bma_fcalc_convection_graduvw
  public :: bma_fcalc_convection_vgraduw
  
  public :: bma_fcalc_rhsConst
  public :: bma_fcalc_rhsBubble
  public :: bma_fcalc_rhsBubblePlusFE
  public :: bma_fcalc_rhsFE
  
  public :: bma_fcalc_integralOne
  public :: bma_fcalc_integralFE
  public :: bma_fcalc_bubbleL2error
  public :: bma_fcalc_bubbleH1error
  public :: bma_fcalc_L2norm
  public :: bma_fcalc_H1norm
  
  public :: bma_fcalc_divergenceL2norm

contains

  !****************************************************************************

!<subroutine>

  subroutine bma_docalc_mass(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,dscale,ix,iy)

!<description>  
    ! Calculates the Mass operator at position (ix,iy) of the
    ! global matrix with a multiplier scale.
!</description>

!<inputoutput>
    ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
    ! have to be filled with data.
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Multiplier in front of the matrix.
    real(DP), intent(in) :: dscale
    
    ! Column/row in the global matrix where to asseble the matrix to.
    integer, intent(in) :: ix,iy
!</input>

!</subroutine>

    ! Local variables
    real(DP) :: dbasI, dbasJ
    integer :: iel, icubp, idofe, jdofe, ivar, nvar
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix
    real(DP), dimension(:,:,:,:), pointer :: p_DlocalMatrixIntl
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial,p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaMatrixData), pointer :: p_rmatrixData

    integer :: ndimfe, idimfe

    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Get local data
    p_rmatrixData => RmatrixData(iy,ix)
    p_DbasTrial => RmatrixData(iy,ix)%p_DbasTrial
    p_DbasTest => RmatrixData(iy,ix)%p_DbasTest
    
    ! FE space dimension
    ndimfe = RmatrixData(iy,ix)%ndimfeTrial
    
    if (ndimfe .ne. RmatrixData(iy,ix)%ndimfeTest) then
      ! This does not make sense.
      call output_line ("Dimension of trial and test FE space different.",&
          OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_mass")
      call sys_halt()
    end if

    ! Interleaved matrix?
    if (.not. p_rmatrixData%bisInterleaved) then

      ! Get the matrix data      
      p_DlocalMatrix => RmatrixData(iy,ix)%p_Dentry

      if (ndimfe .eq. 1) then

        ! -----------------------
        ! Scalar-valued FE space. 
        ! -----------------------
        ! This IF command prevents an inner IF command and speeds 
        ! up the computation for all scalar standard FE spaces.
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Psi_i:
            do idofe=1,p_rmatrixData%ndofTest

              ! Fetch the contributions of the (test) basis functions Psi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the (trial) basis function Phi_j:
              do jdofe=1,p_rmatrixData%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJ = p_DbasTrial(jdofe,DER_FUNC,icubp,iel)

                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * dbasJ*dbasI

              end do ! jdofe

            end do ! idofe

          end do ! icubp

        end do ! iel
        
      else
      
        ! -----------------------
        ! Vector-valued FE space.
        ! -----------------------
        ! This implementation works for all types of finite elements, but it
        ! is slightly slower for scalar FE spaces due to an additional inner loop.

        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Loop over the dimensions of the FE space
            do idimfe = 0,ndimfe-1

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Psi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Psi_i
                ! into dbasI
                dbasI = p_DbasTest(idofe+idimfe*p_rmatrixData%ndofTest,DER_FUNC,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the (trial) basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJ = p_DbasTrial(jdofe+idimfe*p_rmatrixData%ndofTrial,DER_FUNC,icubp,iel)

                  ! Multiply the values of the basis functions
                  ! (1st derivatives) by the cubature weight and sum up
                  ! into the local matrices.
                  p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                      dscale * p_DcubWeight(icubp,iel) * dbasJ*dbasI

                end do ! jdofe

              end do ! idofe
              
            end do ! idimfe

          end do ! icubp

        end do ! iel
        
      end if ! ndimfe = 1

    else

      ! Get the matrix data      
      p_DlocalMatrixIntl => RmatrixData(iy,ix)%p_DentryIntl

      ! Interleave-data
      nvar = RmatrixData(iy,ix)%nvar

      if (ndimfe .eq. 1) then
        
        ! -----------------------
        ! Scalar-valued FE space. 
        ! -----------------------
        ! This IF command prevents an inner IF command and speeds 
        ! up the computation for all scalar standard FE spaces.

        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Psi_i:
            do idofe=1,p_rmatrixData%ndofTest

              ! Fetch the contributions of the (test) basis functions Psi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the (trial) basis function Phi_j:
              do jdofe=1,p_rmatrixData%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJ = p_DbasTrial(jdofe,DER_FUNC,icubp,iel)

                do ivar = 1,nvar

                  ! Multiply the values of the basis functions
                  ! (1st derivatives) by the cubature weight and sum up
                  ! into the local matrices.
                  p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) = &
                      p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) + &
                          dscale * p_DcubWeight(icubp,iel) * dbasJ*dbasI

                end do ! ivar

              end do ! jdofe

            end do ! idofe

          end do ! icubp

        end do ! iel
        
      else

        ! -----------------------
        ! Vector-valued FE space.
        ! -----------------------
        ! This implementation works for all types of finite elements, but it
        ! is slightly slower for scalar FE spaces due to an additional inner loop.

        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Loop over the dimensions of the FE space
            do idimfe = 0,ndimfe-1

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Psi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Psi_i
                ! into dbasI
                dbasI = p_DbasTest(idofe+idimfe*p_rmatrixData%ndofTest,DER_FUNC,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the (trial) basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJ = p_DbasTrial(jdofe+idimfe*p_rmatrixData%ndofTrial,DER_FUNC,icubp,iel)

                  do ivar = 1,nvar

                    ! Multiply the values of the basis functions
                    ! (1st derivatives) by the cubature weight and sum up
                    ! into the local matrices.
                    p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) = &
                        p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) + &
                            dscale * p_DcubWeight(icubp,iel) * dbasJ*dbasI

                  end do ! ivar

                end do ! jdofe

              end do ! idofe
              
            end do ! idimfe

          end do ! icubp

        end do ! iel

      end if

    end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_massDiag(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>
    ! Calculates the Mass operator in all diagonal matrices.
    !
    ! Note: If rcollection is not specified, the matrix is calculated
    ! in all diagonal blocks with a multiplier of 1.
    ! If rcollection is specified, the following parameters are expected:
    ! rcollection%DquickAccess(1) = multiplier in front of the mass matrix.
    ! rcollection%IquickAccess(1) = First diagonal block or =0 to 
    !                               start from block 1.
    ! rcollection%IquickAccess(2) = Last diagonal block or =0 to calculate 
    !                               up to the last.
!</description>

!<inputoutput>
    ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
    ! have to be filled with data.
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!</subroutine>

    integer :: i,istart,iend
    real(DP) :: dscale

    ! First/last block, multiplier
    istart = 1
    iend = min(ubound(RmatrixData,1),ubound(RmatrixData,2))
    dscale = 1.0_DP
    if (present(rcollection)) then
      if (rcollection%IquickAccess(1) .ne. 0) then
        istart = rcollection%IquickAccess(1)
      end if
      if (rcollection%IquickAccess(2) .ne. 0) then
        iend = rcollection%IquickAccess(2)
      end if
      dscale = rcollection%DquickAccess(1)
    end if
    
    ! Loop through all diagonal blocks   
    do i = istart, iend
    
      ! Calculate the Mass matrix at the diagonal position i.
      call bma_docalc_mass(RmatrixData,rassemblyData,rmatrixAssembly,&
          npointsPerElement,nelements,dscale,i,i)

    end do

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_mass(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the Mass operator in all diagonal matrices.
    !
    ! Note: If rcollection is not specified, the matrix is calculated
    ! at position (1,1) with a multiplier of 1.
    ! If rcollection is specified, the following parameters are expected:
    ! rcollection%DquickAccess(1) = multiplier in front of the mass matrix.
    ! rcollection%IquickAccess(1) = x-position in the destination matrix
    ! rcollection%IquickAccess(2) = y-position in the destination matrix
!</description>

!<inputoutput>
    ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
    ! have to be filled with data.
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!</subroutine>

    ! Local variables
    integer :: ix, iy
    real(DP) :: dscale

    ! Get parameters
    dscale = 1.0_DP
    iy = 1
    ix = 1

    if (present(rcollection)) then
      dscale = rcollection%DquickAccess(1)
      ix = rcollection%IquickAccess(1)
      iy = rcollection%IquickAccess(2)

      ! Cancel if nothing to do or parameters wrong
      if (dscale .eq. 0.0_DP) return

      if ((ix .lt. 1) .or. (iy .lt. 1) .or. &
          (ix .gt. ubound(RmatrixData,2)) .or. (iy .gt. ubound(RmatrixData,1))) then
        call output_line ("Parameters wrong.",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_mass")
        call sys_halt()
      end if

    end if

    ! Calculate the Mass matrix at the desired position
    call bma_docalc_mass(RmatrixData,rassemblyData,rmatrixAssembly,&
        npointsPerElement,nelements,dscale,ix,iy)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_docalc_laplace(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,dscale,ix,iy)

!<description>  
    ! Calculates the Laplace operator at position (ix,iy) in a block matrix
    ! with a scaling factor dscale in front.
!</description>

!<inputoutput>
    ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
    ! have to be filled with data.
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Multiplier in front of the matrix.
    real(DP), intent(in) :: dscale
    
    ! Column/row in the global matrix where to asseble the matrix to.
    integer, intent(in) :: ix,iy
!</input>

!</subroutine>

    real(DP) :: dbasIx, dbasJx, dbasIy, dbasJy, dbasIz, dbasJz
    integer :: iel, icubp, idofe, jdofe, ivar, nvar
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix
    real(DP), dimension(:,:,:,:), pointer :: p_DlocalMatrixIntl
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial,p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaMatrixData), pointer :: p_rmatrixData

    integer :: ndimfe, idimfe

    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Get local data
    p_rmatrixData => RmatrixData(iy,ix)
    p_DbasTrial => RmatrixData(iy,ix)%p_DbasTrial
    p_DbasTest => RmatrixData(iy,ix)%p_DbasTest

    ! FE space dimension
    ndimfe = RmatrixData(iy,ix)%ndimfeTrial
    
    if (ndimfe .ne. RmatrixData(iy,ix)%ndimfeTest) then
      ! This does not make sense.
      call output_line ("Dimension of trial and test FE space different.",&
          OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_laplace")
      call sys_halt()
    end if

    ! Interleaved matrix?
    if (.not. p_rmatrixData%bisInterleaved) then

      ! Get the matrix data
      p_DlocalMatrix => RmatrixData(iy,ix)%p_Dentry

      if (ndimfe .eq. 1) then

        ! -----------------------
        ! Scalar-valued FE space. 
        ! -----------------------
        ! This IF command prevents an inner IF command and speeds 
        ! up the computation for all scalar standard FE spaces.

        ! What's the dimension?
        select case (rmatrixAssembly%p_rtriangulation%ndim)

        case (NDIM1D)

          ! *****************
          ! 1D Laplace matrix
          ! *****************

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Psi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Psi_i
                ! into dbasI
                dbasIx = p_DbasTest(idofe,DER_DERIV1D_X,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the (trial) basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJx = p_DbasTrial(jdofe,DER_DERIV1D_X,icubp,iel)

                  ! Multiply the values of the basis functions
                  ! (1st derivatives) by the cubature weight and sum up
                  ! into the local matrices.
                  p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                      dscale * p_DcubWeight(icubp,iel) * dbasJx*dbasIx

                end do ! jdofe

              end do ! idofe

            end do ! icubp

          end do ! iel

        case (NDIM2D)

          ! *****************
          ! 2D Laplace matrix
          ! *****************

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Psi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Psi_i
                ! into dbasI
                dbasIx = p_DbasTest(idofe,DER_DERIV2D_X,icubp,iel)
                dbasIy = p_DbasTest(idofe,DER_DERIV2D_Y,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the (trial) basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJx = p_DbasTrial(jdofe,DER_DERIV2D_X,icubp,iel)
                  dbasJy = p_DbasTrial(jdofe,DER_DERIV2D_Y,icubp,iel)

                  ! Multiply the values of the basis functions
                  ! (1st derivatives) by the cubature weight and sum up
                  ! into the local matrices.
                  p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                      dscale * p_DcubWeight(icubp,iel) * ( dbasJx*dbasIx + dbasJy*dbasIy )

                end do ! jdofe

              end do ! idofe

            end do ! icubp

          end do ! iel

        case (NDIM3D)

          ! *****************
          ! 3D Laplace matrix
          ! *****************

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Psi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Psi_i
                ! into dbasI
                dbasIx = p_DbasTest(idofe,DER_DERIV3D_X,icubp,iel)
                dbasIy = p_DbasTest(idofe,DER_DERIV3D_Y,icubp,iel)
                dbasIz = p_DbasTest(idofe,DER_DERIV3D_Z,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the (trial) basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJx = p_DbasTrial(jdofe,DER_DERIV3D_X,icubp,iel)
                  dbasJy = p_DbasTrial(jdofe,DER_DERIV3D_Y,icubp,iel)
                  dbasJz = p_DbasTrial(jdofe,DER_DERIV3D_Z,icubp,iel)

                  ! Multiply the values of the basis functions
                  ! (1st derivatives) by the cubature weight and sum up
                  ! into the local matrices.
                  p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                      dscale * p_DcubWeight(icubp,iel) * &
                              ( dbasJx*dbasIx + dbasJy*dbasIy + dbasJz*dbasIz )

                end do ! jdofe

              end do ! idofe

            end do ! icubp

          end do ! iel

        end select
        
      else

        ! -----------------------
        ! Vector-valued FE space.
        ! -----------------------
        ! This implementation works for all types of finite elements, but it
        ! is slightly slower for scalar FE spaces due to an additional inner loop.

        ! What's the dimension?
        select case (rmatrixAssembly%p_rtriangulation%ndim)

        case (NDIM1D)

          ! *****************
          ! 1D Laplace matrix
          ! *****************

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Loop over the dimensions of the FE space
              do idimfe = 0,ndimfe-1
                
                ! Outer loop over the DOF's i=1..ndof on our current element,
                ! which corresponds to the (test) basis functions Psi_i:
                do idofe=1,p_rmatrixData%ndofTest

                  ! Fetch the contributions of the (test) basis functions Psi_i
                  ! into dbasI
                  dbasIx = p_DbasTest(&
                      idofe+idimfe*p_rmatrixData%ndofTest,DER_DERIV1D_X,icubp,iel)

                  ! Inner loop over the DOF's j=1..ndof, which corresponds to
                  ! the (trial) basis function Phi_j:
                  do jdofe=1,p_rmatrixData%ndofTrial

                    ! Fetch the contributions of the (trial) basis function Phi_j
                    ! into dbasJ
                    dbasJx = p_DbasTrial(&
                        jdofe+idimfe*p_rmatrixData%ndofTrial,DER_DERIV1D_X,icubp,iel)

                    ! Multiply the values of the basis functions
                    ! (1st derivatives) by the cubature weight and sum up
                    ! into the local matrices.
                    p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                        dscale * p_DcubWeight(icubp,iel) * dbasJx*dbasIx

                  end do ! jdofe

                end do ! idofe
                
              end do ! idimfe

            end do ! icubp

          end do ! iel

        case (NDIM2D)

          ! *****************
          ! 2D Laplace matrix
          ! *****************

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Loop over the dimensions of the FE space
              do idimfe = 0,ndimfe-1
              
                ! Outer loop over the DOF's i=1..ndof on our current element,
                ! which corresponds to the (test) basis functions Psi_i:
                do idofe=1,p_rmatrixData%ndofTest

                  ! Fetch the contributions of the (test) basis functions Psi_i
                  ! into dbasI
                  dbasIx = p_DbasTest(&
                      idofe+idimfe*p_rmatrixData%ndofTest,DER_DERIV2D_X,icubp,iel)
                  dbasIy = p_DbasTest(&
                      idofe+idimfe*p_rmatrixData%ndofTest,DER_DERIV2D_Y,icubp,iel)

                  ! Inner loop over the DOF's j=1..ndof, which corresponds to
                  ! the (trial) basis function Phi_j:
                  do jdofe=1,p_rmatrixData%ndofTrial

                    ! Fetch the contributions of the (trial) basis function Phi_j
                    ! into dbasJ
                    dbasJx = p_DbasTrial(&
                        jdofe+idimfe*p_rmatrixData%ndofTrial,DER_DERIV2D_X,icubp,iel)
                    dbasJy = p_DbasTrial(&
                        jdofe+idimfe*p_rmatrixData%ndofTrial,DER_DERIV2D_Y,icubp,iel)

                    ! Multiply the values of the basis functions
                    ! (1st derivatives) by the cubature weight and sum up
                    ! into the local matrices.
                    p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                        dscale * p_DcubWeight(icubp,iel) * ( dbasJx*dbasIx + dbasJy*dbasIy )

                  end do ! jdofe

                end do ! idofe
                
              end do ! idimfe

            end do ! icubp

          end do ! iel

        case (NDIM3D)

          ! *****************
          ! 3D Laplace matrix
          ! *****************

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Loop over the dimensions of the FE space
              do idimfe = 0,ndimfe-1
              
                ! Outer loop over the DOF's i=1..ndof on our current element,
                ! which corresponds to the (test) basis functions Psi_i:
                do idofe=1,p_rmatrixData%ndofTest

                  ! Fetch the contributions of the (test) basis functions Psi_i
                  ! into dbasI
                  dbasIx = p_DbasTest(&
                      idofe+idimfe*p_rmatrixData%ndofTest,DER_DERIV3D_X,icubp,iel)
                  dbasIy = p_DbasTest(&
                      idofe+idimfe*p_rmatrixData%ndofTest,DER_DERIV3D_Y,icubp,iel)
                  dbasIz = p_DbasTest(&
                      idofe+idimfe*p_rmatrixData%ndofTest,DER_DERIV3D_Z,icubp,iel)

                  ! Inner loop over the DOF's j=1..ndof, which corresponds to
                  ! the (trial) basis function Phi_j:
                  do jdofe=1,p_rmatrixData%ndofTrial

                    ! Fetch the contributions of the (trial) basis function Phi_j
                    ! into dbasJ
                    dbasJx = p_DbasTrial(&
                        jdofe+idimfe*p_rmatrixData%ndofTrial,DER_DERIV3D_X,icubp,iel)
                    dbasJy = p_DbasTrial(&
                        jdofe+idimfe*p_rmatrixData%ndofTrial,DER_DERIV3D_Y,icubp,iel)
                    dbasJz = p_DbasTrial(&
                        jdofe+idimfe*p_rmatrixData%ndofTrial,DER_DERIV3D_Z,icubp,iel)

                    ! Multiply the values of the basis functions
                    ! (1st derivatives) by the cubature weight and sum up
                    ! into the local matrices.
                    p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                        dscale * p_DcubWeight(icubp,iel) * &
                                ( dbasJx*dbasIx + dbasJy*dbasIy + dbasJz*dbasIz )

                  end do ! jdofe

                end do ! idofe
                
              end do ! idimfe

            end do ! icubp

          end do ! iel

        end select

      end if

    else

      ! Get the matrix data      
      p_DlocalMatrixIntl => RmatrixData(iy,ix)%p_DentryIntl

      ! Interleave-data
      nvar = RmatrixData(iy,ix)%nvar

      if (ndimfe .eq. 1) then

        ! -----------------------
        ! Scalar-valued FE space. 
        ! -----------------------
        ! This IF command prevents an inner IF command and speeds 
        ! up the computation for all scalar standard FE spaces.

        ! What's the dimension?
        select case (rmatrixAssembly%p_rtriangulation%ndim)

        case (NDIM1D)

          ! *****************
          ! 1D Laplace matrix
          ! *****************

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Psi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Psi_i
                ! into dbasI
                dbasIx = p_DbasTest(idofe,DER_DERIV1D_X,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the (trial) basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJx = p_DbasTrial(jdofe,DER_DERIV1D_X,icubp,iel)

                  do ivar = 1,nvar

                    ! Multiply the values of the basis functions
                    ! (1st derivatives) by the cubature weight and sum up
                    ! into the local matrices.
                    p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) = &
                        p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) + &
                          dscale * p_DcubWeight(icubp,iel) * dbasJx*dbasIx

                  end do ! ivar

                end do ! jdofe

              end do ! idofe

            end do ! icubp

          end do ! iel

        case (NDIM2D)

          ! *****************
          ! 2D Laplace matrix
          ! *****************

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Psi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Psi_i
                ! into dbasI
                dbasIx = p_DbasTest(idofe,DER_DERIV2D_X,icubp,iel)
                dbasIy = p_DbasTest(idofe,DER_DERIV2D_Y,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the (trial) basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJx = p_DbasTrial(jdofe,DER_DERIV2D_X,icubp,iel)
                  dbasJy = p_DbasTrial(jdofe,DER_DERIV2D_Y,icubp,iel)

                  do ivar = 1,nvar

                    ! Multiply the values of the basis functions
                    ! (1st derivatives) by the cubature weight and sum up
                    ! into the local matrices.
                    p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) = &
                        p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) + &
                            dscale * p_DcubWeight(icubp,iel) * ( dbasJx*dbasIx + dbasJy*dbasIy )

                  end do ! ivar

                end do ! jdofe

              end do ! idofe

            end do ! icubp

          end do ! iel

        case (NDIM3D)

          ! *****************
          ! 3D Laplace matrix
          ! *****************

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Psi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Psi_i
                ! into dbasI
                dbasIx = p_DbasTest(idofe,DER_DERIV3D_X,icubp,iel)
                dbasIy = p_DbasTest(idofe,DER_DERIV3D_Y,icubp,iel)
                dbasIz = p_DbasTest(idofe,DER_DERIV3D_Z,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the (trial) basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJx = p_DbasTrial(jdofe,DER_DERIV3D_X,icubp,iel)
                  dbasJy = p_DbasTrial(jdofe,DER_DERIV3D_Y,icubp,iel)
                  dbasJz = p_DbasTrial(jdofe,DER_DERIV3D_Z,icubp,iel)

                  do ivar = 1,nvar

                    ! Multiply the values of the basis functions
                    ! (1st derivatives) by the cubature weight and sum up
                    ! into the local matrices.
                    p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) = &
                        p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) + &
                            dscale * p_DcubWeight(icubp,iel) * &
                                    ( dbasJx*dbasIx + dbasJy*dbasIy + dbasJz*dbasIz )

                  end do ! ivar

                end do ! jdofe

              end do ! idofe

            end do ! icubp

          end do ! iel

        end select

      else
      
        ! -----------------------
        ! Vector-valued FE space.
        ! -----------------------
        ! This implementation works for all types of finite elements, but it
        ! is slightly slower for scalar FE spaces due to an additional inner loop.
      
        ! What's the dimension?
        select case (rmatrixAssembly%p_rtriangulation%ndim)

        case (NDIM1D)

          ! *****************
          ! 1D Laplace matrix
          ! *****************

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Loop over the dimensions of the FE space
              do idimfe = 0,ndimfe-1

                ! Outer loop over the DOF's i=1..ndof on our current element,
                ! which corresponds to the (test) basis functions Psi_i:
                do idofe=1,p_rmatrixData%ndofTest

                  ! Fetch the contributions of the (test) basis functions Psi_i
                  ! into dbasI
                  dbasIx = p_DbasTest(&
                      idofe+idimfe*p_rmatrixData%ndofTest,DER_DERIV1D_X,icubp,iel)

                  ! Inner loop over the DOF's j=1..ndof, which corresponds to
                  ! the (trial) basis function Phi_j:
                  do jdofe=1,p_rmatrixData%ndofTrial

                    ! Fetch the contributions of the (trial) basis function Phi_j
                    ! into dbasJ
                    dbasJx = p_DbasTrial(&
                        jdofe+idimfe*p_rmatrixData%ndofTrial,DER_DERIV1D_X,icubp,iel)

                    do ivar = 1,nvar

                      ! Multiply the values of the basis functions
                      ! (1st derivatives) by the cubature weight and sum up
                      ! into the local matrices.
                      p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) = &
                          p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) + &
                            dscale * p_DcubWeight(icubp,iel) * dbasJx*dbasIx

                    end do ! ivar

                  end do ! jdofe

                end do ! idofe
                
              end do ! idimfe

            end do ! icubp

          end do ! iel

        case (NDIM2D)

          ! *****************
          ! 2D Laplace matrix
          ! *****************

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Loop over the dimensions of the FE space
              do idimfe = 0,ndimfe-1

                ! Outer loop over the DOF's i=1..ndof on our current element,
                ! which corresponds to the (test) basis functions Psi_i:
                do idofe=1,p_rmatrixData%ndofTest

                  ! Fetch the contributions of the (test) basis functions Psi_i
                  ! into dbasI
                  dbasIx = p_DbasTest(&
                      idofe+idimfe*p_rmatrixData%ndofTest,DER_DERIV2D_X,icubp,iel)
                  dbasIy = p_DbasTest(&
                      idofe+idimfe*p_rmatrixData%ndofTest,DER_DERIV2D_Y,icubp,iel)

                  ! Inner loop over the DOF's j=1..ndof, which corresponds to
                  ! the (trial) basis function Phi_j:
                  do jdofe=1,p_rmatrixData%ndofTrial

                    ! Fetch the contributions of the (trial) basis function Phi_j
                    ! into dbasJ
                    dbasJx = p_DbasTrial(&
                        jdofe+idimfe*p_rmatrixData%ndofTrial,DER_DERIV2D_X,icubp,iel)
                    dbasJy = p_DbasTrial(&
                        jdofe+idimfe*p_rmatrixData%ndofTrial,DER_DERIV2D_Y,icubp,iel)

                    do ivar = 1,nvar

                      ! Multiply the values of the basis functions
                      ! (1st derivatives) by the cubature weight and sum up
                      ! into the local matrices.
                      p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) = &
                          p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) + &
                              dscale * p_DcubWeight(icubp,iel) * ( dbasJx*dbasIx + dbasJy*dbasIy )

                    end do ! ivar

                  end do ! jdofe

                end do ! idofe
                
              end do ! idimfe

            end do ! icubp

          end do ! iel

        case (NDIM3D)

          ! *****************
          ! 3D Laplace matrix
          ! *****************

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Loop over the dimensions of the FE space
              do idimfe = 0,ndimfe-1

                ! Outer loop over the DOF's i=1..ndof on our current element,
                ! which corresponds to the (test) basis functions Psi_i:
                do idofe=1,p_rmatrixData%ndofTest

                  ! Fetch the contributions of the (test) basis functions Psi_i
                  ! into dbasI
                  dbasIx = p_DbasTest(&
                      idofe+idimfe*p_rmatrixData%ndofTest,DER_DERIV3D_X,icubp,iel)
                  dbasIy = p_DbasTest(&
                      idofe+idimfe*p_rmatrixData%ndofTest,DER_DERIV3D_Y,icubp,iel)
                  dbasIz = p_DbasTest(&
                      idofe+idimfe*p_rmatrixData%ndofTest,DER_DERIV3D_Z,icubp,iel)

                  ! Inner loop over the DOF's j=1..ndof, which corresponds to
                  ! the (trial) basis function Phi_j:
                  do jdofe=1,p_rmatrixData%ndofTrial

                    ! Fetch the contributions of the (trial) basis function Phi_j
                    ! into dbasJ
                    dbasJx = p_DbasTrial(&
                        jdofe+idimfe*p_rmatrixData%ndofTrial,DER_DERIV3D_X,icubp,iel)
                    dbasJy = p_DbasTrial(&
                        jdofe+idimfe*p_rmatrixData%ndofTrial,DER_DERIV3D_Y,icubp,iel)
                    dbasJz = p_DbasTrial(&
                        jdofe+idimfe*p_rmatrixData%ndofTrial,DER_DERIV3D_Z,icubp,iel)

                    do ivar = 1,nvar

                      ! Multiply the values of the basis functions
                      ! (1st derivatives) by the cubature weight and sum up
                      ! into the local matrices.
                      p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) = &
                          p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) + &
                              dscale * p_DcubWeight(icubp,iel) * &
                                      ( dbasJx*dbasIx + dbasJy*dbasIy + dbasJz*dbasIz )

                    end do ! ivar

                  end do ! jdofe

                end do ! idofe
                
              end do ! idimfe

            end do ! icubp

          end do ! iel

        end select

      end if      

    end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_laplaceDiag(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the Laplace operator in all diagonal matrices
    !
    ! Note: If rcollection is not specified, the matrix is calculated
    ! in all diagonal blocks with a multiplier of 1.
    ! If rcollection is specified, the following parameters are expected:
    ! rcollection%DquickAccess(1) = multiplier in front of the matrix.
    ! rcollection%IquickAccess(1) = First diagonal block or =0 to 
    !                               start from block 1.
    ! rcollection%IquickAccess(2) = Last diagonal block or =0 to calculate 
    !                               up to the last.
!</description>

!<inputoutput>
    ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
    ! have to be filled with data.
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!</subroutine>

    integer :: i,istart,iend
    real(DP) :: dscale

    ! First/last block, multiplier
    istart = 1
    iend = min(ubound(RmatrixData,1),ubound(RmatrixData,2))
    dscale = 1.0_DP
    if (present(rcollection)) then
      if (rcollection%IquickAccess(1) .ne. 0) then
        istart = rcollection%IquickAccess(1)
      end if
      if (rcollection%IquickAccess(2) .ne. 0) then
        iend = rcollection%IquickAccess(2)
      end if
      dscale = rcollection%DquickAccess(1)
    end if
    
    ! Loop through all diagonal blocks   
    do i = istart, iend
    
      ! Calculate the Laplace at the diagonal position i.
      call bma_docalc_laplace(RmatrixData,rassemblyData,rmatrixAssembly,&
          npointsPerElement,nelements,1.0_DP,i,i)

    end do

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_laplace(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the Laplace operator at position (x,y) in a block matrix.
    !
    ! Note: If rcollection is not specified, the matrix is calculated
    ! at matrix position (1,1) with a multiplier of 1.
    ! If rcollection is specified, the following parameters are expected:
    ! rcollection%DquickAccess(1) = multiplier in front of the matrix.
    ! rcollection%IquickAccess(1) = x-coordinate in the block matrix
    ! rcollection%IquickAccess(2) = y-coordinate in the block matrix
!</description>

!<inputoutput>
    ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
    ! have to be filled with data.
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!</subroutine>

    ! local variables
    integer :: ix, iy
    real(DP) :: dscale

    ! Get parameters
    dscale = 1.0_DP
    iy = 1
    ix = 1

    if (present(rcollection)) then
      dscale = rcollection%DquickAccess(1)
      ix = rcollection%IquickAccess(1)
      iy = rcollection%IquickAccess(2)

      ! Cancel if nothing to do or parameters wrong
      if (dscale .eq. 0.0_DP) return

      if ((ix .lt. 1) .or. (iy .lt. 1) .or. &
          (ix .gt. ubound(RmatrixData,2)) .or. (iy .gt. ubound(RmatrixData,1))) then
        call output_line ("Parameters wrong.",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_laplace")
        call sys_halt()
      end if

    end if

    ! Calculate the Laplace at the desired position.
    call bma_docalc_laplace(RmatrixData,rassemblyData,rmatrixAssembly,&
        npointsPerElement,nelements,1.0_DP,ix,iy)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_docalc_gradientdiv(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,dscale,ix,iy)

!<description>  
    ! Calculates the gradient operator at position (ix,iy) in a block matrix
    ! with a scaling factor dscale in front.
    ! Uses the "transposed divergence approach", i.e., not the actual
    ! gradient is set up but the derivative is put to the test functions.
    !
    !   ( - p, div w )  ( =  ( grad(p), w )  after partial integration )
    !
    ! For scalar-valued v-spaces, the gradient operator is
    ! imposed to the matrices (iy,ix), (iy+1,ix), (iy+2,ix). It is assumed that
    ! the subspaces in v_i are all of the same FE type.
    ! For vector-valued v-spaces, the gradient operator is imposed
    ! to the matrix at position (ix,iy).
!</description>

!<inputoutput>
    ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
    ! have to be filled with data.
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Multiplier in front of the matrix.
    real(DP), intent(in) :: dscale
    
    ! Column/row in the global matrix where to asseble the matrix to.
    integer, intent(in) :: ix,iy
!</input>

!</subroutine>

    real(DP) :: dbasIx, dbasJ, dbasIy, dbasIz
    integer :: iel, icubp, idofe, jdofe, ivar, nvar
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix,p_DlocalMatrix2,p_DlocalMatrix3
    real(DP), dimension(:,:,:,:), pointer :: p_DlocalMatrixIntl,p_DlocalMatrixIntl2,p_DlocalMatrixIntl3
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial,p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaMatrixData), pointer :: p_rmatrixData

    integer :: ndimfe
    real(DP) :: dscaleLocal

    ! Scaling factor is multiplied by -1 as we have the "divergence" formulation.
    dscaleLocal = -dscale

    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Get local data
    p_rmatrixData => RmatrixData(iy,ix)
    p_DbasTrial => RmatrixData(iy,ix)%p_DbasTrial
    p_DbasTest => RmatrixData(iy,ix)%p_DbasTest

    ! FE space dimension
    ndimfe = RmatrixData(iy,ix)%ndimfeTest
    
    if (RmatrixData(iy,ix)%ndimfeTrial .ne. 1) then
      ! This does not make sense.
      call output_line ("Trial-space must be scalar-valued.",&
          OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_gradientdiv")
      call sys_halt()
    end if

    ! Interleaved matrix?
    if (.not. p_rmatrixData%bisInterleaved) then

      if (ndimfe .eq. 1) then

        ! -----------------------
        ! Scalar-valued FE space. 
        ! -----------------------
        ! This IF command prevents an inner IF command and speeds 
        ! up the computation for all scalar standard FE spaces.

        ! What's the dimension?
        select case (rmatrixAssembly%p_rtriangulation%ndim)

        case (NDIM1D)

          ! *****************
          ! 1D gradient matrix
          ! *****************

          ! Get the matrix data
          p_DlocalMatrix => RmatrixData(iy,ix)%p_Dentry

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Psi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Psi_i
                ! into dbasI
                dbasIx = p_DbasTest(idofe,DER_DERIV1D_X,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the (trial) basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJ = p_DbasTrial(jdofe,DER_FUNC1D,icubp,iel)

                  ! Multiply the values of the basis functions
                  ! (1st derivatives) by the cubature weight and sum up
                  ! into the local matrices.
                  p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                      dscaleLocal * p_DcubWeight(icubp,iel) * dbasJ*dbasIx

                end do ! jdofe

              end do ! idofe

            end do ! icubp

          end do ! iel

        case (NDIM2D)

          ! *****************
          ! 2D gradient matrix
          ! *****************

          ! Get the matrix data
          p_DlocalMatrix => RmatrixData(iy,ix)%p_Dentry
          p_DlocalMatrix2 => RmatrixData(iy+1,ix)%p_Dentry

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Psi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Psi_i
                ! into dbasI
                dbasIx = p_DbasTest(idofe,DER_DERIV2D_X,icubp,iel)
                dbasIy = p_DbasTest(idofe,DER_DERIV2D_Y,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the (trial) basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJ = p_DbasTrial(jdofe,DER_FUNC2D,icubp,iel)

                  ! Multiply the values of the basis functions
                  ! (1st derivatives) by the cubature weight and sum up
                  ! into the local matrices.
                  p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                      dscaleLocal * p_DcubWeight(icubp,iel) * dbasJ * dbasIx

                  p_DlocalMatrix2(jdofe,idofe,iel) = p_DlocalMatrix2(jdofe,idofe,iel) + &
                      dscaleLocal * p_DcubWeight(icubp,iel) * dbasJ * dbasIy

                end do ! jdofe

              end do ! idofe

            end do ! icubp

          end do ! iel

        case (NDIM3D)

          ! *****************
          ! 3D gradient matrix
          ! *****************

          ! Get the matrix data
          p_DlocalMatrix => RmatrixData(iy,ix)%p_Dentry
          p_DlocalMatrix2 => RmatrixData(iy+1,ix)%p_Dentry
          p_DlocalMatrix3 => RmatrixData(iy+2,ix)%p_Dentry

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Psi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Psi_i
                ! into dbasI
                dbasIx = p_DbasTest(idofe,DER_DERIV3D_X,icubp,iel)
                dbasIy = p_DbasTest(idofe,DER_DERIV3D_Y,icubp,iel)
                dbasIz = p_DbasTest(idofe,DER_DERIV3D_Z,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the (trial) basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJ = p_DbasTrial(jdofe,DER_FUNC3D,icubp,iel)

                  ! Multiply the values of the basis functions
                  ! (1st derivatives) by the cubature weight and sum up
                  ! into the local matrices.
                  p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                      dscaleLocal * p_DcubWeight(icubp,iel) * dbasJ * dbasIx
                      
                  p_DlocalMatrix2(jdofe,idofe,iel) = p_DlocalMatrix2(jdofe,idofe,iel) + &
                      dscaleLocal * p_DcubWeight(icubp,iel) * dbasJ * dbasIy
                      
                  p_DlocalMatrix3(jdofe,idofe,iel) = p_DlocalMatrix3(jdofe,idofe,iel) + &
                      dscaleLocal * p_DcubWeight(icubp,iel) * dbasJ * dbasIz

                end do ! jdofe

              end do ! idofe

            end do ! icubp

          end do ! iel

        end select
        
      else

        ! -----------------------
        ! Vector-valued FE space.
        ! -----------------------
        ! This implementation works for all types of finite elements, but it
        ! is slightly slower for scalar FE spaces due to an additional inner loop.

        ! Get the matrix data
        p_DlocalMatrix => RmatrixData(iy,ix)%p_Dentry

        ! What's the dimension?
        select case (rmatrixAssembly%p_rtriangulation%ndim)

        case (NDIM1D)

          ! *****************
          ! 1D gradient matrix
          ! *****************

          ! Does not apply. A vector-valued FE-space in 1D is scalar-valued,
          ! so this case will never be called.

        case (NDIM2D)

          ! *****************
          ! 2D gradient matrix
          ! *****************

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Psi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Psi_i
                ! into dbasI
                dbasIx = p_DbasTest(&
                    idofe+0*p_rmatrixData%ndofTest,DER_DERIV2D_X,icubp,iel)
                dbasIy = p_DbasTest(&
                    idofe+1*p_rmatrixData%ndofTest,DER_DERIV2D_Y,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the (trial) basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJ = p_DbasTrial(jdofe,DER_FUNC2D,icubp,iel)

                  ! Multiply the values of the basis functions
                  ! (1st derivatives) by the cubature weight and sum up
                  ! into the local matrices.
                  p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                      dscaleLocal * p_DcubWeight(icubp,iel) * dbasJ * ( dbasIx + dbasIy )

                end do ! jdofe

              end do ! idofe
                
            end do ! icubp

          end do ! iel

        case (NDIM3D)

          ! *****************
          ! 3D gradient matrix
          ! *****************

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Psi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Psi_i
                ! into dbasI
                dbasIx = p_DbasTest(&
                    idofe+0*p_rmatrixData%ndofTest,DER_DERIV3D_X,icubp,iel)
                dbasIy = p_DbasTest(&
                    idofe+1*p_rmatrixData%ndofTest,DER_DERIV3D_Y,icubp,iel)
                dbasIz = p_DbasTest(&
                    idofe+2*p_rmatrixData%ndofTest,DER_DERIV3D_Z,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the (trial) basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJ = p_DbasTrial(jdofe,DER_FUNC3D,icubp,iel)

                  ! Multiply the values of the basis functions
                  ! (1st derivatives) by the cubature weight and sum up
                  ! into the local matrices.
                  p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                      dscaleLocal * p_DcubWeight(icubp,iel) * dbasJ * ( dbasIx + dbasIy + dbasIz )

                end do ! jdofe

              end do ! idofe
                
            end do ! icubp

          end do ! iel

        end select

      end if

    else

      ! Interleave-data
      nvar = RmatrixData(iy,ix)%nvar

      if (ndimfe .eq. 1) then

        ! -----------------------
        ! Scalar-valued FE space. 
        ! -----------------------
        ! This IF command prevents an inner IF command and speeds 
        ! up the computation for all scalar standard FE spaces.

        ! What's the dimension?
        select case (rmatrixAssembly%p_rtriangulation%ndim)

        case (NDIM1D)

          ! *****************
          ! 1D gradient matrix
          ! *****************

          ! Get the matrix data      
          p_DlocalMatrixIntl => RmatrixData(iy,ix)%p_DentryIntl

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Psi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Psi_i
                ! into dbasI
                dbasIx = p_DbasTest(idofe,DER_DERIV1D_X,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the (trial) basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJ = p_DbasTrial(jdofe,DER_FUNC1D,icubp,iel)

                  do ivar = 1,nvar

                    ! Multiply the values of the basis functions
                    ! (1st derivatives) by the cubature weight and sum up
                    ! into the local matrices.
                    p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) = &
                        p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) + &
                          dscaleLocal * p_DcubWeight(icubp,iel) * dbasJ*dbasIx

                  end do ! ivar

                end do ! jdofe

              end do ! idofe

            end do ! icubp

          end do ! iel

        case (NDIM2D)

          ! *****************
          ! 2D gradient matrix
          ! *****************

          ! Get the matrix data      
          p_DlocalMatrixIntl => RmatrixData(iy,ix)%p_DentryIntl
          p_DlocalMatrixIntl2 => RmatrixData(iy+1,ix)%p_DentryIntl

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Psi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Psi_i
                ! into dbasI
                dbasIx = p_DbasTest(idofe,DER_DERIV2D_X,icubp,iel)
                dbasIy = p_DbasTest(idofe,DER_DERIV2D_Y,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the (trial) basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJ = p_DbasTrial(jdofe,DER_FUNC2D,icubp,iel)

                  do ivar = 1,nvar

                    ! Multiply the values of the basis functions
                    ! (1st derivatives) by the cubature weight and sum up
                    ! into the local matrices.
                    p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) = &
                        p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) + &
                            dscaleLocal * p_DcubWeight(icubp,iel) * dbasJ * dbasIx

                    p_DlocalMatrixIntl2(ivar,jdofe,idofe,iel) = &
                        p_DlocalMatrixIntl2(ivar,jdofe,idofe,iel) + &
                            dscaleLocal * p_DcubWeight(icubp,iel) * dbasJ * dbasIy

                  end do ! ivar

                end do ! jdofe

              end do ! idofe

            end do ! icubp

          end do ! iel

        case (NDIM3D)

          ! *****************
          ! 3D gradient matrix
          ! *****************

          ! Get the matrix data      
          p_DlocalMatrixIntl => RmatrixData(iy,ix)%p_DentryIntl
          p_DlocalMatrixIntl2 => RmatrixData(iy+1,ix)%p_DentryIntl
          p_DlocalMatrixIntl3 => RmatrixData(iy+2,ix)%p_DentryIntl

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Psi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Psi_i
                ! into dbasI
                dbasIx = p_DbasTest(idofe,DER_DERIV3D_X,icubp,iel)
                dbasIy = p_DbasTest(idofe,DER_DERIV3D_Y,icubp,iel)
                dbasIz = p_DbasTest(idofe,DER_DERIV3D_Z,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the (trial) basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJ = p_DbasTrial(jdofe,DER_FUNC3D,icubp,iel)

                  do ivar = 1,nvar

                    ! Multiply the values of the basis functions
                    ! (1st derivatives) by the cubature weight and sum up
                    ! into the local matrices.
                    p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) = &
                        p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) + &
                            dscaleLocal * p_DcubWeight(icubp,iel) * dbasJ * dbasIx

                    p_DlocalMatrixIntl2(ivar,jdofe,idofe,iel) = &
                        p_DlocalMatrixIntl2(ivar,jdofe,idofe,iel) + &
                            dscaleLocal * p_DcubWeight(icubp,iel) * dbasJ * dbasIy

                    p_DlocalMatrixIntl3(ivar,jdofe,idofe,iel) = &
                        p_DlocalMatrixIntl3(ivar,jdofe,idofe,iel) + &
                            dscaleLocal * p_DcubWeight(icubp,iel) * dbasJ * dbasIz
                  end do ! ivar

                end do ! jdofe

              end do ! idofe

            end do ! icubp

          end do ! iel

        end select

      else
      
        ! -----------------------
        ! Vector-valued FE space.
        ! -----------------------
        ! This implementation works for all types of finite elements, but it
        ! is slightly slower for scalar FE spaces due to an additional inner loop.
      
        ! Get the matrix data      
        p_DlocalMatrixIntl => RmatrixData(iy,ix)%p_DentryIntl

        ! What's the dimension?
        select case (rmatrixAssembly%p_rtriangulation%ndim)

        case (NDIM1D)
        
          ! *****************
          ! 1D gradient matrix
          ! *****************
          
          ! Does not apply. A vector-valued FE-space in 1D is scalar-valued,
          ! so this case will never be called.

        case (NDIM2D)

          ! *****************
          ! 2D gradient matrix
          ! *****************

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Psi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Psi_i
                ! into dbasI
                dbasIx = p_DbasTest(&
                    idofe+0*p_rmatrixData%ndofTest,DER_DERIV2D_X,icubp,iel)
                dbasIy = p_DbasTest(&
                    idofe+1*p_rmatrixData%ndofTest,DER_DERIV2D_Y,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the (trial) basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJ = p_DbasTrial(jdofe,DER_FUNC2D,icubp,iel)

                  do ivar = 1,nvar

                    ! Multiply the values of the basis functions
                    ! (1st derivatives) by the cubature weight and sum up
                    ! into the local matrices.
                    p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) = &
                        p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) + &
                            dscaleLocal * p_DcubWeight(icubp,iel) * dbasJ * ( dbasIx + dbasIy )

                  end do ! ivar

                end do ! jdofe

              end do ! idofe
                
            end do ! icubp

          end do ! iel

        case (NDIM3D)

          ! *****************
          ! 3D gradient matrix
          ! *****************

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Psi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Psi_i
                ! into dbasI
                dbasIx = p_DbasTest(&
                    idofe+0*p_rmatrixData%ndofTest,DER_DERIV3D_X,icubp,iel)
                dbasIy = p_DbasTest(&
                    idofe+1*p_rmatrixData%ndofTest,DER_DERIV3D_Y,icubp,iel)
                dbasIz = p_DbasTest(&
                    idofe+2*p_rmatrixData%ndofTest,DER_DERIV3D_Z,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the (trial) basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJ = p_DbasTrial(jdofe,DER_FUNC3D,icubp,iel)

                  do ivar = 1,nvar

                    ! Multiply the values of the basis functions
                    ! (1st derivatives) by the cubature weight and sum up
                    ! into the local matrices.
                    p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) = &
                        p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) + &
                            dscaleLocal * p_DcubWeight(icubp,iel) * dbasJ * ( dbasIx + dbasIy + dbasIz )

                  end do ! ivar

                end do ! jdofe

              end do ! idimfe

            end do ! icubp

          end do ! iel

        end select

      end if      

    end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_gradientdiv(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the gradient operator at position (x,y) in a block matrix.
    ! Uses the "transposed divergence approach", i.e., not the actual
    ! gradient is set up but the derivative is put to the test functions.
    !
    !   ( - p, div w )  ( =  ( grad(p), w )  after partial integration )
    !
    ! Note: If rcollection is not specified, the matrix is calculated
    ! at matrix position (1,1) with a multiplier of 1.
    ! If rcollection is specified, the following parameters are expected:
    !   rcollection%DquickAccess(1) = multiplier in front of the matrix.
    !   rcollection%IquickAccess(1) = x-coordinate in the block matrix
    !   rcollection%IquickAccess(2) = y-coordinate in the block matrix
    ! For scalar-valued v-spaces, the gradient operator is
    ! imposed to the matrices (y,x), (y+1,x), (y+2,x). It is assumed that
    ! the subspaces in v_i are all of the same FE type.
    ! For vector-valued v-spaces, the gradient operator is imposed
    ! to the matrix at position (x,y).
!</description>

!<inputoutput>
    ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
    ! have to be filled with data.
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!</subroutine>

    ! local variables
    integer :: ix, iy
    real(DP) :: dscale

    ! Get parameters
    dscale = 1.0_DP
    iy = 1
    ix = 1

    if (present(rcollection)) then

      dscale = rcollection%DquickAccess(1)
      ix = rcollection%IquickAccess(1)
      iy = rcollection%IquickAccess(2)

      ! Cancel if nothing to do or parameters wrong
      if (dscale .eq. 0.0_DP) return

      if ((ix .lt. 1) .or. (iy .lt. 1) .or. &
          (ix .gt. ubound(RmatrixData,2)) .or. (iy .gt. ubound(RmatrixData,1))) then
        call output_line ("Parameters wrong.",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_gradientdiv")
        call sys_halt()
      end if

    end if
    
    call bma_docalc_gradientdiv(RmatrixData,rassemblyData,rmatrixAssembly,&
        npointsPerElement,nelements,dscale,ix,iy)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_docalc_divergence(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,dscale,ix,iy)

!<description>  
    ! Calculates the divergence operator at position (ix,iy) in a block matrix:
    !
    !    ( div v, q )
    !
    ! For scalar-valued v-spaces, the divergence operator is
    ! imposed to the matrices (iy,ix), (iy,ix+1), (iy,ix+2). It is assumed that
    ! the subspaces in v_i are all of the same FE type.
    ! For vector-valued v-spaces, the divergence operator is imposed
    ! to the matrix at position (ix,iy).
!</description>

!<inputoutput>
    ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
    ! have to be filled with data.
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Multiplier in front of the matrix.
    real(DP), intent(in) :: dscale
    
    ! Column/row in the global matrix where to asseble the matrix to.
    integer, intent(in) :: ix,iy
!</input>

!</subroutine>

    real(DP) :: dbasI, dbasJx, dbasJy, dbasJz
    integer :: iel, icubp, idofe, jdofe, ivar, nvar
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix,p_DlocalMatrix2,p_DlocalMatrix3
    real(DP), dimension(:,:,:,:), pointer :: p_DlocalMatrixIntl,p_DlocalMatrixIntl2,p_DlocalMatrixIntl3
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial,p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaMatrixData), pointer :: p_rmatrixData

    integer :: ndimfe

    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Get local data
    p_rmatrixData => RmatrixData(iy,ix)
    p_DbasTrial => RmatrixData(iy,ix)%p_DbasTrial
    p_DbasTest => RmatrixData(iy,ix)%p_DbasTest

    ! FE space dimension
    ndimfe = RmatrixData(iy,ix)%ndimfeTrial
    
    if (RmatrixData(iy,ix)%ndimfeTest .ne. 1) then
      ! This does not make sense.
      call output_line ("Test-space must be scalar-valued.",&
          OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_gradientdiv")
      call sys_halt()
    end if

    ! Interleaved matrix?
    if (.not. p_rmatrixData%bisInterleaved) then

      if (ndimfe .eq. 1) then

        ! -----------------------
        ! Scalar-valued FE space. 
        ! -----------------------
        ! This IF command prevents an inner IF command and speeds 
        ! up the computation for all scalar standard FE spaces.

        ! What's the dimension?
        select case (rmatrixAssembly%p_rtriangulation%ndim)

        case (NDIM1D)

          ! *****************
          ! 1D divergence matrix
          ! *****************

          ! Get the matrix data
          p_DlocalMatrix => RmatrixData(iy,ix)%p_Dentry

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Psi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Psi_i
                ! into dbasI
                dbasI = p_DbasTest(idofe,DER_FUNC1D,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the (trial) basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJx = p_DbasTrial(jdofe,DER_DERIV1D_X,icubp,iel)

                  ! Multiply the values of the basis functions
                  ! (1st derivatives) by the cubature weight and sum up
                  ! into the local matrices.
                  p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                      dscale * p_DcubWeight(icubp,iel) * dbasJx*dbasI

                end do ! jdofe

              end do ! idofe

            end do ! icubp

          end do ! iel

        case (NDIM2D)

          ! *****************
          ! 2D divergence matrix
          ! *****************

          ! Get the matrix data
          p_DlocalMatrix => RmatrixData(iy,ix)%p_Dentry
          p_DlocalMatrix2 => RmatrixData(iy+1,ix)%p_Dentry

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Psi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Psi_i
                ! into dbasI
                dbasI = p_DbasTest(idofe,DER_FUNC2D,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the (trial) basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJx = p_DbasTrial(jdofe,DER_DERIV2D_X,icubp,iel)
                  dbasJy = p_DbasTrial(jdofe,DER_DERIV2D_Y,icubp,iel)

                  ! Multiply the values of the basis functions
                  ! (1st derivatives) by the cubature weight and sum up
                  ! into the local matrices.
                  p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                      dscale * p_DcubWeight(icubp,iel) * dbasJx * dbasI

                  p_DlocalMatrix2(jdofe,idofe,iel) = p_DlocalMatrix2(jdofe,idofe,iel) + &
                      dscale * p_DcubWeight(icubp,iel) * dbasJy * dbasI

                end do ! jdofe

              end do ! idofe

            end do ! icubp

          end do ! iel

        case (NDIM3D)

          ! *****************
          ! 3D divergence matrix
          ! *****************

          ! Get the matrix data
          p_DlocalMatrix => RmatrixData(iy,ix)%p_Dentry
          p_DlocalMatrix2 => RmatrixData(iy+1,ix)%p_Dentry
          p_DlocalMatrix3 => RmatrixData(iy+2,ix)%p_Dentry

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Psi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Psi_i
                ! into dbasI
                dbasI = p_DbasTest(idofe,DER_FUNC3D,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the (trial) basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJx = p_DbasTrial(jdofe,DER_DERIV3D_X,icubp,iel)
                  dbasJy = p_DbasTrial(jdofe,DER_DERIV3D_Y,icubp,iel)
                  dbasJz = p_DbasTrial(jdofe,DER_DERIV3D_Z,icubp,iel)

                  ! Multiply the values of the basis functions
                  ! (1st derivatives) by the cubature weight and sum up
                  ! into the local matrices.
                  p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                      dscale * p_DcubWeight(icubp,iel) * dbasJx * dbasI

                  p_DlocalMatrix2(jdofe,idofe,iel) = p_DlocalMatrix2(jdofe,idofe,iel) + &
                      dscale * p_DcubWeight(icubp,iel) * dbasJy * dbasI

                  p_DlocalMatrix3(jdofe,idofe,iel) = p_DlocalMatrix3(jdofe,idofe,iel) + &
                      dscale * p_DcubWeight(icubp,iel) * dbasJz * dbasI
                end do ! jdofe

              end do ! idofe

            end do ! icubp

          end do ! iel

        end select
        
      else

        ! -----------------------
        ! Vector-valued FE space.
        ! -----------------------
        ! This implementation works for all types of finite elements, but it
        ! is slightly slower for scalar FE spaces due to an additional inner loop.

        ! Get the matrix data
        p_DlocalMatrix => RmatrixData(iy,ix)%p_Dentry

        ! What's the dimension?
        select case (rmatrixAssembly%p_rtriangulation%ndim)

        case (NDIM1D)

          ! *****************
          ! 1D divergence matrix
          ! *****************

          ! Does not apply. A vector-valued FE-space in 1D is scalar-valued,
          ! so this case will never be called.

        case (NDIM2D)

          ! *****************
          ! 2D divergence matrix
          ! *****************

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Psi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Psi_i
                ! into dbasI
                dbasI = p_DbasTest(idofe,DER_FUNC2D,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the (trial) basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJx = p_DbasTrial(&
                      jdofe+0*p_rmatrixData%ndofTrial,DER_DERIV2D_X,icubp,iel)
                  dbasJy = p_DbasTrial(&
                      jdofe+1*p_rmatrixData%ndofTrial,DER_DERIV2D_Y,icubp,iel)

                  ! Multiply the values of the basis functions
                  ! (1st derivatives) by the cubature weight and sum up
                  ! into the local matrices.
                  p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                      dscale * p_DcubWeight(icubp,iel) * ( dbasJx + dbasJy ) * dbasI

                end do ! jdofe

              end do ! idofe
                
            end do ! icubp

          end do ! iel

        case (NDIM3D)

          ! *****************
          ! 3D divergence matrix
          ! *****************

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Psi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Psi_i
                ! into dbasI
                dbasI = p_DbasTest(idofe,DER_FUNC3D,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the (trial) basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJx = p_DbasTrial(&
                      jdofe+0*p_rmatrixData%ndofTrial,DER_DERIV3D_X,icubp,iel)
                  dbasJy = p_DbasTrial(&
                      jdofe+1*p_rmatrixData%ndofTrial,DER_DERIV3D_Y,icubp,iel)
                  dbasJz = p_DbasTrial(&
                      jdofe+2*p_rmatrixData%ndofTrial,DER_DERIV3D_Z,icubp,iel)

                  ! Multiply the values of the basis functions
                  ! (1st derivatives) by the cubature weight and sum up
                  ! into the local matrices.
                  p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                      dscale * p_DcubWeight(icubp,iel) * &
                              ( dbasJx + dbasJy + dbasJz ) * dbasI

                end do ! jdofe

              end do ! idofe

            end do ! icubp

          end do ! iel

        end select

      end if

    else

      ! Interleave-data
      nvar = RmatrixData(iy,ix)%nvar

      if (ndimfe .eq. 1) then

        ! -----------------------
        ! Scalar-valued FE space. 
        ! -----------------------
        ! This IF command prevents an inner IF command and speeds 
        ! up the computation for all scalar standard FE spaces.

        ! What's the dimension?
        select case (rmatrixAssembly%p_rtriangulation%ndim)

        case (NDIM1D)

          ! *****************
          ! 1D divergence matrix
          ! *****************

          ! Get the matrix data      
          p_DlocalMatrixIntl => RmatrixData(iy,ix)%p_DentryIntl

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Psi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Psi_i
                ! into dbasI
                dbasI = p_DbasTest(idofe,DER_FUNC1D,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the (trial) basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJx = p_DbasTrial(jdofe,DER_DERIV1D_X,icubp,iel)

                  do ivar = 1,nvar

                    ! Multiply the values of the basis functions
                    ! (1st derivatives) by the cubature weight and sum up
                    ! into the local matrices.
                    p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) = &
                        p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) + &
                          dscale * p_DcubWeight(icubp,iel) * dbasJx * dbasI

                  end do ! ivar

                end do ! jdofe

              end do ! idofe

            end do ! icubp

          end do ! iel

        case (NDIM2D)

          ! *****************
          ! 2D divergence matrix
          ! *****************

          ! Get the matrix data      
          p_DlocalMatrixIntl => RmatrixData(iy,ix)%p_DentryIntl
          p_DlocalMatrixIntl2 => RmatrixData(iy+1,ix)%p_DentryIntl

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Psi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Psi_i
                ! into dbasI
                dbasI = p_DbasTest(idofe,DER_FUNC2D,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the (trial) basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJx = p_DbasTrial(jdofe,DER_DERIV2D_X,icubp,iel)
                  dbasJy = p_DbasTrial(jdofe,DER_DERIV2D_Y,icubp,iel)

                  do ivar = 1,nvar

                    ! Multiply the values of the basis functions
                    ! (1st derivatives) by the cubature weight and sum up
                    ! into the local matrices.
                    p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) = &
                        p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) + &
                            dscale * p_DcubWeight(icubp,iel) * dbasJx * dbasI

                    p_DlocalMatrixIntl2(ivar,jdofe,idofe,iel) = &
                        p_DlocalMatrixIntl2(ivar,jdofe,idofe,iel) + &
                            dscale * p_DcubWeight(icubp,iel) * dbasJy * dbasI

                  end do ! ivar

                end do ! jdofe

              end do ! idofe

            end do ! icubp

          end do ! iel

        case (NDIM3D)

          ! *****************
          ! 3D divergence matrix
          ! *****************

          ! Get the matrix data      
          p_DlocalMatrixIntl => RmatrixData(iy,ix)%p_DentryIntl
          p_DlocalMatrixIntl2 => RmatrixData(iy+1,ix)%p_DentryIntl
          p_DlocalMatrixIntl3 => RmatrixData(iy+2,ix)%p_DentryIntl

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Psi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Psi_i
                ! into dbasI
                dbasI = p_DbasTest(idofe,DER_FUNC3D,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the (trial) basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJx = p_DbasTrial(jdofe,DER_DERIV3D_X,icubp,iel)
                  dbasJy = p_DbasTrial(jdofe,DER_DERIV3D_Y,icubp,iel)
                  dbasJz = p_DbasTrial(jdofe,DER_DERIV3D_Z,icubp,iel)

                  do ivar = 1,nvar

                    ! Multiply the values of the basis functions
                    ! (1st derivatives) by the cubature weight and sum up
                    ! into the local matrices.
                    p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) = &
                        p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) + &
                            dscale * p_DcubWeight(icubp,iel) * dbasJx * dbasI

                    p_DlocalMatrixIntl2(ivar,jdofe,idofe,iel) = &
                        p_DlocalMatrixIntl2(ivar,jdofe,idofe,iel) + &
                            dscale * p_DcubWeight(icubp,iel) * dbasJy * dbasI

                    p_DlocalMatrixIntl3(ivar,jdofe,idofe,iel) = &
                        p_DlocalMatrixIntl3(ivar,jdofe,idofe,iel) + &
                            dscale * p_DcubWeight(icubp,iel) * dbasJz * dbasI
                  end do ! ivar

                end do ! jdofe

              end do ! idofe

            end do ! icubp

          end do ! iel

        end select

      else
      
        ! -----------------------
        ! Vector-valued FE space.
        ! -----------------------
        ! This implementation works for all types of finite elements, but it
        ! is slightly slower for scalar FE spaces due to an additional inner loop.
      
        ! Get the matrix data      
        p_DlocalMatrixIntl => RmatrixData(iy,ix)%p_DentryIntl

        ! What's the dimension?
        select case (rmatrixAssembly%p_rtriangulation%ndim)

        case (NDIM1D)

          ! *****************
          ! 1D divergence matrix
          ! *****************

          ! Does not apply. A vector-valued FE-space in 1D is scalar-valued,
          ! so this case will never be called.

        case (NDIM2D)

          ! *****************
          ! 2D divergence matrix
          ! *****************

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Psi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Psi_i
                ! into dbasI
                dbasI = p_DbasTest(&
                    idofe,DER_FUNC2D,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the (trial) basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJx = p_DbasTrial(&
                      jdofe+0*p_rmatrixData%ndofTrial,DER_DERIV2D_X,icubp,iel)
                  dbasJy = p_DbasTrial(&
                      jdofe+1*p_rmatrixData%ndofTrial,DER_DERIV2D_Y,icubp,iel)

                  do ivar = 1,nvar

                    ! Multiply the values of the basis functions
                    ! (1st derivatives) by the cubature weight and sum up
                    ! into the local matrices.
                    p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) = &
                        p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) + &
                            dscale * p_DcubWeight(icubp,iel) * ( dbasJx + dbasJy ) * dbasI

                  end do ! ivar

                end do ! jdofe

              end do ! idofe

            end do ! icubp

          end do ! iel

        case (NDIM3D)

          ! *****************
          ! 3D divergence matrix
          ! *****************

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Psi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Psi_i
                ! into dbasI
                dbasI = p_DbasTest(idofe,DER_FUNC3D,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the (trial) basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJx = p_DbasTrial(&
                      jdofe+0*p_rmatrixData%ndofTrial,DER_DERIV3D_X,icubp,iel)
                  dbasJy = p_DbasTrial(&
                      jdofe+1*p_rmatrixData%ndofTrial,DER_DERIV3D_Y,icubp,iel)
                  dbasJz = p_DbasTrial(&
                      jdofe+2*p_rmatrixData%ndofTrial,DER_DERIV3D_Z,icubp,iel)

                  do ivar = 1,nvar

                    ! Multiply the values of the basis functions
                    ! (1st derivatives) by the cubature weight and sum up
                    ! into the local matrices.
                    p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) = &
                        p_DlocalMatrixIntl(ivar,jdofe,idofe,iel) + &
                            dscale * p_DcubWeight(icubp,iel) * &
                            ( dbasJx + dbasJy + dbasJz ) * dbasI

                  end do ! ivar

                end do ! jdofe

              end do ! idofe

            end do ! icubp

          end do ! iel

        end select

      end if      

    end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_divergence(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the divergence operator at position (x,y) in a block matrix:
    !
    !    ( div v, q )
    !
    ! Note: If rcollection is not specified, the matrix is calculated
    ! at matrix position (1,1) with a multiplier of 1.
    ! If rcollection is specified, the following parameters are expected:
    !   rcollection%DquickAccess(1) = multiplier in front of the matrix.
    !   rcollection%IquickAccess(1) = x-coordinate in the block matrix
    !   rcollection%IquickAccess(2) = y-coordinate in the block matrix
    ! For scalar-valued v-spaces, the divergence operator is
    ! imposed to the matrices (y,x), (y,x+1), (y,x+2). It is assumed that
    ! the subspaces in v_i are all of the same FE type.
    ! For vector-valued v-spaces, the divergence operator is imposed
    ! to the matrix at position (x,y).
!</description>

!<inputoutput>
    ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
    ! have to be filled with data.
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!</subroutine>

    integer :: ix, iy
    real(DP) :: dscale

    ! Get parameters
    dscale = 1.0_DP
    iy = 1
    ix = 1

    if (present(rcollection)) then
      dscale = rcollection%DquickAccess(1)
      ix = rcollection%IquickAccess(1)
      iy = rcollection%IquickAccess(2)

      ! Cancel if nothing to do or parameters wrong
      if (dscale .eq. 0.0_DP) return

      if ((ix .lt. 1) .or. (iy .lt. 1) .or. &
          (ix .gt. ubound(RmatrixData,2)) .or. (iy .gt. ubound(RmatrixData,1))) then
        call output_line ("Parameters wrong.",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_divergence")
        call sys_halt()
      end if

    end if
    
    ! Calculate the matrix.
    call bma_docalc_divergence(RmatrixData,rassemblyData,rmatrixAssembly,&
        npointsPerElement,nelements,dscale,ix,iy)

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_docalc_convection_ugradvw(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,dscale,ix,iy,du1,du2,du3,rvectorField)

!<description>  
    ! Calculates a convection operator "( (u grad) v, w)" at position (x,y)
    ! of a block matrix, with a convection u given as a finite element function
    ! and a scaling factor dscale in front.
    ! v is the trial and w the test basis function.
    !
    ! If rvectorField is not specified, a constant velocity du1/du2/du3
    ! is assumed, and du1/du2/du3 must be present.
    ! If rvectorField is specified, it describes the underlying vector field.
!</description>

!</remarks>
    ! Remark 1:
    ! The routine currently assumes that all velocity components are discretised
    ! with the same FEM space.
    !
    ! Remark 2:
    ! The routine currently assumes that all velocity matrices are independent.
    ! Matrices sharing data are not supported. This cound be realised by
    ! taking care of the flags RmatrixData(:,:)%bsharedMatrixData which indicate
    ! which matrix data is shared.
    !
    ! Remark 3:
    ! Interleaved matrices are currently not supported.
!</remarks>

!<inputoutput>
    ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
    ! have to be filled with data.
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements
    
    ! Scaling factor
    real(DP), intent(in) :: dscale
    
    ! Position in the global matrix where to assemble
    integer, intent(in) :: ix,iy
    
    ! OPTIONAL: Constant velocity. Must be specified if rvectorField is not
    ! present. Ignored if rvectorField is present.
    real(DP), intent(in), optional :: du1,du2,du3

    ! OPTIONAL: Evaluation structure that describes an underlying nonconstant
    ! velocity field. Can be omitted if du1/du2/du3 is given.
    type(t_fev2VectorData), intent(in), optional :: rvectorField
!</input>

!</subroutine>

    ! Local variables
    real(DP) :: dbasI, dbasJx, dbasJy, dbasJz
    integer :: iel, icubp, idofe, jdofe
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix11,p_DlocalMatrix22,p_DlocalMatrix33
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial,p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaMatrixData), pointer :: p_rmatrixData11,p_rmatrixData22,p_rmatrixData33
    real(DP), dimension(:,:,:,:), pointer :: p_Du
    real(DP) :: dvelX, dvelY, dvelZ

    integer :: ndim
    logical :: bvelConst

    ! Dimension of the underlying space
    ndim = rmatrixAssembly%p_rtriangulation%ndim
    
    ! Get parameters
    bvelConst = present(rvectorField)
    if (bvelConst) then
      dvelX = 0.0_DP
      dvelY = 0.0_DP
      dvelZ = 0.0_DP
    
      ! Constant velocity
      if (present(du1)) dvelX = du1
      if (present(du2)) dvelY = du2
      if (present(du3)) dvelZ = du3
    else
      ! Get the velocity field from the parameters
      p_Du => rvectorField%p_DdataVec
    end if

    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Get local data
    p_DbasTrial => RmatrixData(iy,ix)%p_DbasTrial
    p_DbasTest => RmatrixData(iy,ix)%p_DbasTest

    ! Set up the local matrix of the convection.
    select case (ndim)
    case (NDIM1D)
      ! Matrices to be set up
      p_rmatrixData11 => RmatrixData(iy,ix)
      
      p_DlocalMatrix11 => RmatrixData(iy,ix)%p_Dentry
      
      ! Currently, interleaved matrices are not supported
      if (p_rmatrixData11%bisInterleaved) then
        call output_line ("Interleaved matrices not supported",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_docalc_convection_ugradvw")
        call sys_halt()
      end if

      if ((p_rmatrixData11%ndimfeTrial .ne. 1) .or. &
          (p_rmatrixData11%ndimfeTest .ne. 1)) then
        ! This does not make sense.
        call output_line ("Only scalar-valued FE spaces supported.",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_docalc_convection_ugradvw")
        call sys_halt()
      end if

      if (bvelConst) then
      
        ! Set up the matrix for constant velocity.
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Psi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Psi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the (trial) basis function Phi_j:
              do jdofe=1,p_rmatrixData11%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJx = p_DbasTrial(jdofe,DER_DERIV1D_X,icubp,iel)
                
                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX * dbasJx * dbasI )            ! ( u1 phi_x , psi )

              end do ! jdofe

            end do ! idofe

          end do ! icubp

        end do ! iel
        
      else

        ! Set up the matrix for nonconstant velocity.
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Velocity field in this cubature point
            dvelX = p_Du(1,icubp,iel,DER_FUNC)
            
            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Psi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Psi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the (trial) basis function Phi_j:
              do jdofe=1,p_rmatrixData11%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJx = p_DbasTrial(jdofe,DER_DERIV1D_X,icubp,iel)
                
                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX * dbasJx * dbasI )            ! ( u1 phi_x , psi )

              end do ! jdofe

            end do ! idofe

          end do ! icubp

        end do ! iel
        
      end if
      
    case (NDIM2D)

      ! Matrices to be set up
      p_rmatrixData11 => RmatrixData(iy,ix)
      p_rmatrixData22 => RmatrixData(iy+1,ix+1)
      
      p_DlocalMatrix11 => RmatrixData(iy,ix)%p_Dentry
      p_DlocalMatrix22 => RmatrixData(iy+1,ix+1)%p_Dentry

      ! Currently, interleaved matrices are not supported
      if (p_rmatrixData11%bisInterleaved .or. p_rmatrixData22%bisInterleaved) then
        call output_line ("Interleaved matrices not supported",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_docalc_convection_ugradvw")
        call sys_halt()
      end if

      if ((p_rmatrixData11%ndimfeTrial .ne. 1) .or. &
          (p_rmatrixData11%ndimfeTest .ne. 1)) then
        ! This does not make sense.
        call output_line ("Only scalar-valued FE spaces supported.",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_docalc_convection_ugradvw")
        call sys_halt()
      end if

      if (bvelConst) then
      
        ! Set up the matrix for constant velocity.
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Psi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Psi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the (trial) basis function Phi_j:
              do jdofe=1,p_rmatrixData11%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJx = p_DbasTrial(jdofe,DER_DERIV2D_X,icubp,iel)
                dbasJy = p_DbasTrial(jdofe,DER_DERIV2D_Y,icubp,iel)
                
                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX * dbasJx * dbasI + &          ! ( u1 phi_x , psi )
                      dvelY * dbasJy * dbasI )            ! ( u2 phi_y , psi )

                p_DlocalMatrix22(jdofe,idofe,iel) = p_DlocalMatrix22(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX * dbasJx * dbasI + &          ! ( u1 phi_x , psi )
                      dvelY * dbasJy * dbasI )            ! ( u2 phi_y , psi )

              end do ! jdofe

            end do ! idofe

          end do ! icubp

        end do ! iel
        
      else

        ! Set up the matrix for nonconstant velocity.
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Velocity field in this cubature point
            dvelX = p_Du(1,icubp,iel,DER_FUNC)
            dvelY = p_Du(2,icubp,iel,DER_FUNC)
            
            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Psi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Psi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the (trial) basis function Phi_j:
              do jdofe=1,p_rmatrixData11%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJx = p_DbasTrial(jdofe,DER_DERIV2D_X,icubp,iel)
                dbasJy = p_DbasTrial(jdofe,DER_DERIV2D_Y,icubp,iel)
                
                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX * dbasJx * dbasI + &          ! ( u1 phi_x , psi )
                      dvelY * dbasJy * dbasI )            ! ( u2 phi_y , psi )

                p_DlocalMatrix22(jdofe,idofe,iel) = p_DlocalMatrix22(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX * dbasJx * dbasI + &          ! ( u1 phi_x , psi )
                      dvelY * dbasJy * dbasI )            ! ( u2 phi_y , psi )

              end do ! jdofe

            end do ! idofe

          end do ! icubp

        end do ! iel
        
      end if

    case (NDIM3D)

      ! Matrices to be set up
      p_rmatrixData11 => RmatrixData(iy,ix)
      p_rmatrixData22 => RmatrixData(iy+1,ix+1)
      p_rmatrixData33 => RmatrixData(iy+2,ix+2)

      p_DlocalMatrix11 => RmatrixData(iy,ix)%p_Dentry
      p_DlocalMatrix22 => RmatrixData(iy+1,ix+1)%p_Dentry
      p_DlocalMatrix33 => RmatrixData(iy+2,ix+2)%p_Dentry

      ! Currently, interleaved matrices are not supported
      if (p_rmatrixData11%bisInterleaved .or. &
          p_rmatrixData22%bisInterleaved .or. &
          p_rmatrixData33%bisInterleaved) then
        call output_line ("Interleaved matrices not supported",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_docalc_convection_ugradvw")
        call sys_halt()
      end if

      if ((p_rmatrixData11%ndimfeTrial .ne. 1) .or. &
          (p_rmatrixData11%ndimfeTest .ne. 1)) then
        ! This does not make sense.
        call output_line ("Only scalar-valued FE spaces supported.",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_docalc_convection_ugradvw")
        call sys_halt()
      end if

      if (bvelConst) then
      
        ! Set up the matrix for constant velocity.
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Psi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Psi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the (trial) basis function Phi_j:
              do jdofe=1,p_rmatrixData11%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJx = p_DbasTrial(jdofe,DER_DERIV3D_X,icubp,iel)
                dbasJy = p_DbasTrial(jdofe,DER_DERIV3D_Y,icubp,iel)
                dbasJz = p_DbasTrial(jdofe,DER_DERIV3D_Z,icubp,iel)
                
                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX * dbasJx * dbasI + &          ! ( u1 phi_x , psi )
                      dvelY * dbasJy * dbasI + &          ! ( u2 phi_y , psi )
                      dvelZ * dbasJz * dbasI )            ! ( u3 phi_z , psi )

                p_DlocalMatrix22(jdofe,idofe,iel) = p_DlocalMatrix22(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX * dbasJx * dbasI + &          ! ( u1 phi_x , psi )
                      dvelY * dbasJy * dbasI + &          ! ( u2 phi_y , psi )
                      dvelZ * dbasJz * dbasI )            ! ( u3 phi_z , psi )

                p_DlocalMatrix33(jdofe,idofe,iel) = p_DlocalMatrix33(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX * dbasJx * dbasI + &          ! ( u1 phi_x , psi )
                      dvelY * dbasJy * dbasI + &          ! ( u2 phi_y , psi )
                      dvelZ * dbasJz * dbasI )            ! ( u3 phi_z , psi )

              end do ! jdofe

            end do ! idofe

          end do ! icubp

        end do ! iel
        
      else

        ! Set up the matrix for nonconstant velocity.
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Velocity field in this cubature point
            dvelX = p_Du(1,icubp,iel,DER_FUNC)
            dvelY = p_Du(2,icubp,iel,DER_FUNC)
            dvelZ = p_Du(3,icubp,iel,DER_FUNC)
            
            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Psi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Psi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the (trial) basis function Phi_j:
              do jdofe=1,p_rmatrixData11%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJx = p_DbasTrial(jdofe,DER_DERIV3D_X,icubp,iel)
                dbasJy = p_DbasTrial(jdofe,DER_DERIV3D_Y,icubp,iel)
                dbasJz = p_DbasTrial(jdofe,DER_DERIV3D_Z,icubp,iel)
                
                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX * dbasJx * dbasI + &          ! ( u1 phi_x , psi )
                      dvelY * dbasJy * dbasI + &          ! ( u2 phi_y , psi )
                      dvelZ * dbasJz * dbasI )            ! ( u3 phi_z , psi )

                p_DlocalMatrix22(jdofe,idofe,iel) = p_DlocalMatrix22(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX * dbasJx * dbasI + &          ! ( u1 phi_x , psi )
                      dvelY * dbasJy * dbasI + &          ! ( u2 phi_y , psi )
                      dvelZ * dbasJz * dbasI )            ! ( u3 phi_z , psi )

                p_DlocalMatrix33(jdofe,idofe,iel) = p_DlocalMatrix33(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX * dbasJx * dbasI + &          ! ( u1 phi_x , psi )
                      dvelY * dbasJy * dbasI + &          ! ( u2 phi_y , psi )
                      dvelZ * dbasJz * dbasI )            ! ( u3 phi_z , psi )

              end do ! jdofe

            end do ! idofe

          end do ! icubp

        end do ! iel
        
      end if

    end select

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_convection_ugradvw(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates a convection operator "( (u grad) v, w)" at position (x,y)
    ! of a block matrix, with a convection u given as a finite element function.
    ! v is the trial and w the test basis function.
    !
    ! Note: If rcollection is not specified, the matrix is calculated
    ! in all diagonal blocks with a multiplier of 1.
    ! If rcollection is specified, the following parameters are expected:
    !
    ! rcollection%DquickAccess(1) = multiplier in front of the operator.
    ! rcollection%IquickAccess(1) = x-position in the matrix where to set up the operator.
    ! rcollection%IquickAccess(2) = y-position in the matrix where to set up the operator.
    ! rcollection%IquickAccess(3) = 0, if the convection is a constant vector field.
    !                                  in this case:
    !                                  1D: rcollection%DquickAccess(2)   = x-velocity
    !                                  2D: rcollection%DquickAccess(2:3) = x/y-velocity
    !                                  3D: rcollection%DquickAccess(2:4) = x/y/z-velocity
    !                             = 1, if the convection is specified by a
    !                                  finite element velocity field. In this case,
    !                                  a finite element velocity field must be specified
    !                                  as parameter revalVectors to the call of 
    !                                  bma_buildMatrix.
!</description>

!<remarks>
    ! Remark 1: 
    ! Using the routines from feevaluation2, it is possible to specify
    ! a nonconstant velocity field. The corresponding code looks as follows:
    !
    ! <verb>
    !     use feevaluation2
    !
    !     ...
    !
    !     type(t_collection) :: rcollection
    !     type(t_scalarCubatureInfo) :: rcubatureInfo   ! Cubature formula
    !     type(t_matrixBlock) :: rmatrix                ! Matrix to be calculated
    !
    !     type(t_vectorBlock) :: rvelocity        ! The velocity field
    !     type(t_fev2Vectors) :: revalVectors     ! Collection of vectors to evaluate
    !
    !     ! Prepare the cubature formula
    !     call spdiscr_createDefCubStructure (..., rcubatureInfo, CUB_GEN_AUTO)
    !
    !     ...
    !
    !     rcollection%IquickAccess(1) = 1          ! x-Position in the matrix
    !     rcollection%IquickAccess(2) = 1          ! y-Position in the matrix
    !     rcollection%IquickAccess(3) = 1          ! Nonconstant viscosity
    !     rcollection%DquickAccess(1) = 1.0_DP     ! Scaling
    !
    !     ! Add the X-, Y- and Z-velocity to revalVectors
    !     call fev2_addVectorFieldToEvalList(revalVectors,0,&
    !         rvelocity%RvectorBlock(1),rvelocity%RvectorBlock(2),rvelocity%RvectorBlock(3))
    !
    !     ! Set up the matrix
    !     call bma_buildMatrix (rmatrix,BMA_CALC_STANDARD, bma_fcalc_convection_ugradvw, &
    !         rcollection, revalVectors=revalVectors,rcubatureInfo=rcubatureInfo)
    !
    !     ! Release the vector structure
    !     call fev2_releaseVectorList(revalVectors)
    !
    !     ...
    !
    !     ! Release the cubature formula
    !     call spdiscr_releaseCubStructure (rcubatureInfo)
    ! </verb>
    !
    ! Remark 2: 
    ! The routine <verb>fev2_addVectorFieldToEvalList</verb> allows to define
    ! the evaluation of derivatives as well. In 3D, e.g., one may apply
    !
    ! <verb>
    !     call fev2_addVectorFieldToEvalList(revalVectors,1,&
    !         rvelocity%RvectorBlock(1),rvelocity%RvectorBlock(2),rvelocity%RvectorBlock(3))
    ! </verb>
    !
    ! which calculates function values as well as 1st derivatives of the complete
    ! vector field (due to the "1" at the end). The calculated values in the
    ! cubature points can then be found in the "p_Ddata" elements of revalVectors:
    !
    ! <verb>
    !   revalVectors%p_RvectorData(1)%p_DdataVec(1,:,:,DER_FUNC)      = u1
    !   revalVectors%p_RvectorData(1)%p_DdataVec(1,:,:,DER_DERIV3D_X) = d/dx u1
    !   revalVectors%p_RvectorData(1)%p_DdataVec(1,:,:,DER_DERIV3D_Y) = d/dy u1
    !   revalVectors%p_RvectorData(1)%p_DdataVec(1,:,:,DER_DERIV3D_Z) = d/dz u1
    !
    !   revalVectors%p_RvectorData(1)%p_DdataVec(2,:,:,DER_FUNC)      = u2
    !   revalVectors%p_RvectorData(1)%p_DdataVec(2,:,:,DER_DERIV3D_X) = d/dx u2
    !   revalVectors%p_RvectorData(1)%p_DdataVec(2,:,:,DER_DERIV3D_Y) = d/dy u2
    !   revalVectors%p_RvectorData(1)%p_DdataVec(2,:,:,DER_DERIV3D_Z) = d/dz u2
    !
    !   revalVectors%p_RvectorData(1)%p_DdataVec(3,:,:,DER_FUNC)      = u3
    !   revalVectors%p_RvectorData(1)%p_DdataVec(3,:,:,DER_DERIV3D_X) = d/dx u3
    !   revalVectors%p_RvectorData(1)%p_DdataVec(3,:,:,DER_DERIV3D_Y) = d/dy u3
    !   revalVectors%p_RvectorData(1)%p_DdataVec(3,:,:,DER_DERIV3D_Z) = d/dz u3
    ! </verb>
    !
    ! in all cubature points on all elements. The vector data in
    ! revalVectors%p_RvectorData appears exactly in the order, the vectors
    ! are added to revalVectors by fev2_addVectorFiueldToEvalList.
    !
    ! Remark 3:
    ! The routine currently assumes that all velocity components are discretised
    ! with the same FEM space.
    !
    ! Remark 4:
    ! The routine currently assumes that all velocity matrices are independent.
    ! Matrices sharing data are not supported. This cound be realised by
    ! taking care of the flags RmatrixData(:,:)%bsharedMatrixData which indicate
    ! which matrix data is shared.
    !
    ! Remark 5:
    ! Interleaved matrices are currently not supported.
!</remarks>

!<inputoutput>
    ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
    ! have to be filled with data.
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!</subroutine>

    ! Local variables
    integer :: ndim,ix,iy
    real(DP) :: dscale
    real(DP) :: dvelX, dvelY, dvelZ
    logical :: bvelConst

    ! Dimension of the underlying space
    ndim = rmatrixAssembly%p_rtriangulation%ndim
    
    ! Get parameters
    dscale = 1.0_DP
    dvelX = 0.0_DP
    dvelY = 0.0_DP
    dvelZ = 0.0_DP
    bvelConst = .true.
    ix = 1
    iy = 1
    
    if (present(rcollection)) then
      dscale = rcollection%DquickAccess(1)
      
      ! Constant velocity?
      bvelConst = (rcollection%IquickAccess(3) .eq. 0)
      
      if (bvelConst) then
        ! Get the constant velocity
        dvelX = rcollection%DquickAccess(2)
        if (ndim .ge. NDIM2D) dvelY = rcollection%DquickAccess(3)
        if (ndim .ge. NDIM3D) dvelZ = rcollection%DquickAccess(4)
      end if
      
      ! Position
      ix = rcollection%IquickAccess(1)
      iy = rcollection%IquickAccess(2)

    end if
    
    if (bvelConst) then
      call bma_docalc_convection_ugradvw(RmatrixData,rassemblyData,rmatrixAssembly,&
          npointsPerElement,nelements,dscale,ix,iy,dvelX,dvelY,dvelZ)
    else
      call bma_docalc_convection_ugradvw(RmatrixData,rassemblyData,rmatrixAssembly,&
          npointsPerElement,nelements,dscale,ix,iy,rvectorField=revalVectors%p_RvectorData(1))
    end if

  end subroutine
    
  !****************************************************************************

!<subroutine>

  subroutine bma_docalc_convection_vugradw(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,dscale,ix,iy,du1,du2,du3,rvectorField)

!<description>  
    ! Calculates a convection operator "( v, (u grad) w )" at position (x,y)
    ! of a block matrix, with a convection u given as a finite element function
    ! and a scaling factor dscale in front.
    ! v is the trial and w the test basis function.
    !
    ! If rvectorField is not specified, a constant velocity du1/du2/du3
    ! is assumed, and du1/du2/du3 must be present.
    ! If rvectorField is specified, it describes the underlying vector field.
!</description>

!</remarks>
    ! Remark 1:
    ! The routine currently assumes that all velocity components are discretised
    ! with the same FEM space.
    !
    ! Remark 2:
    ! The routine currently assumes that all velocity matrices are independent.
    ! Matrices sharing data are not supported. This cound be realised by
    ! taking care of the flags RmatrixData(:,:)%bsharedMatrixData which indicate
    ! which matrix data is shared.
    !
    ! Remark 3:
    ! Interleaved matrices are currently not supported.
!</remarks>

!<inputoutput>
    ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
    ! have to be filled with data.
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements
    
    ! Scaling factor
    real(DP), intent(in) :: dscale
    
    ! Position in the global matrix where to assemble
    integer, intent(in) :: ix,iy
    
    ! OPTIONAL: Constant velocity. Must be specified if rvectorField is not
    ! present. Ignored if rvectorField is present.
    real(DP), intent(in), optional :: du1,du2,du3

    ! OPTIONAL: Evaluation structure that describes an underlying nonconstant
    ! velocity field. Can be omitted if du1/du2/du3 is given.
    type(t_fev2VectorData), intent(in), optional :: rvectorField
!</input>

!</subroutine>

    ! Local variables
    real(DP) :: dbasJ, dbasIx, dbasIy, dbasIz
    integer :: iel, icubp, idofe, jdofe
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix11,p_DlocalMatrix22,p_DlocalMatrix33
    type(t_bmaMatrixData), pointer :: p_rmatrixData11,p_rmatrixData22,p_rmatrixData33
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial,p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:,:,:), pointer :: p_Du

    integer :: ndim
    real(DP) :: dvelX,dvelY,dvelZ
    logical :: bvelConst

    ! Dimension of the underlying space
    ndim = rmatrixAssembly%p_rtriangulation%ndim
    
    ! Get parameters
    bvelConst = present(rvectorField)
    if (bvelConst) then
      dvelX = 0.0_DP
      dvelY = 0.0_DP
      dvelZ = 0.0_DP

      ! Constant velocity
      if (present(du1)) dvelX = du1
      if (present(du2)) dvelY = du2
      if (present(du3)) dvelZ = du3
    else
      ! Get the velocity field from the parameters
      p_Du => rvectorField%p_DdataVec
    end if
    
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Get local data
    p_DbasTrial => RmatrixData(iy,ix)%p_DbasTrial
    p_DbasTest => RmatrixData(iy,ix)%p_DbasTest

    ! Set up the local matrix of the convection.
    select case (ndim)
    case (NDIM1D)
      ! Matrices to be set up
      p_rmatrixData11 => RmatrixData(iy,ix)
      
      p_DlocalMatrix11 => RmatrixData(iy,ix)%p_Dentry
      
      ! Currently, interleaved matrices are not supported
      if (p_rmatrixData11%bisInterleaved) then
        call output_line ("Interleaved matrices not supported",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_docalc_convection_vugradw")
        call sys_halt()
      end if

      if ((p_rmatrixData11%ndimfeTrial .ne. 1) .or. &
          (p_rmatrixData11%ndimfeTest .ne. 1)) then
        ! This does not make sense.
        call output_line ("Only scalar-valued FE spaces supported.",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_docalc_convection_vugradw")
        call sys_halt()
      end if

      if (bvelConst) then
      
        ! Set up the matrix for constant velocity.
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Psi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Psi_i
              ! into dbasI
              dbasIx = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the (trial) basis function Phi_j:
              do jdofe=1,p_rmatrixData11%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJ = p_DbasTrial(jdofe,DER_DERIV1D_X,icubp,iel)
                
                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelX * dbasIx )            ! ( phi , u1 psi_x )

              end do ! jdofe

            end do ! idofe

          end do ! icubp

        end do ! iel
        
      else

        ! Set up the matrix for nonconstant velocity.
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Velocity field in this cubature point
            dvelX = p_Du(1,icubp,iel,DER_FUNC)
            
            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Psi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Psi_i
              ! into dbasI
              dbasIx = p_DbasTest(idofe,DER_DERIV1D_X,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the (trial) basis function Phi_j:
              do jdofe=1,p_rmatrixData11%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJ = p_DbasTrial(jdofe,DER_FUNC,icubp,iel)
                
                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelX * dbasIx )            ! ( phi , u1 psi_x )

              end do ! jdofe

            end do ! idofe

          end do ! icubp

        end do ! iel
        
      end if
      
    case (NDIM2D)

      ! Matrices to be set up
      p_rmatrixData11 => RmatrixData(iy,ix)
      p_rmatrixData22 => RmatrixData(iy+1,ix+1)
      
      p_DlocalMatrix11 => RmatrixData(iy,ix)%p_Dentry
      p_DlocalMatrix22 => RmatrixData(iy+1,ix+1)%p_Dentry

      ! Currently, interleaved matrices are not supported
      if (p_rmatrixData11%bisInterleaved .or. p_rmatrixData22%bisInterleaved) then
        call output_line ("Interleaved matrices not supported",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_docalc_convection_vugradw")
        call sys_halt()
      end if

      if ((p_rmatrixData11%ndimfeTrial .ne. 1) .or. &
          (p_rmatrixData11%ndimfeTest .ne. 1)) then
        ! This does not make sense.
        call output_line ("Only scalar-valued FE spaces supported.",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_docalc_convection_vugradw")
        call sys_halt()
      end if

      if (bvelConst) then
      
        ! Set up the matrix for constant velocity.
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Psi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Psi_i
              ! into dbasI
              dbasIx = p_DbasTest(idofe,DER_DERIV2D_X,icubp,iel)
              dbasIy = p_DbasTest(idofe,DER_DERIV2D_Y,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the (trial) basis function Phi_j:
              do jdofe=1,p_rmatrixData11%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJ = p_DbasTrial(jdofe,DER_FUNC,icubp,iel)
                
                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelX * dbasIx + &          ! ( phi , u1 psi_x )
                      dbasJ * dvelY * dbasIy )            ! ( phi , u2 psi_y )

                p_DlocalMatrix22(jdofe,idofe,iel) = p_DlocalMatrix22(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelX * dbasIx + &          ! ( phi , u1 psi_x )
                      dbasJ * dvelY * dbasIy )            ! ( phi , u2 psi_y )

              end do ! jdofe

            end do ! idofe

          end do ! icubp

        end do ! iel
        
      else

        ! Set up the matrix for nonconstant velocity.
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Velocity field in this cubature point
            dvelX = p_Du(1,icubp,iel,DER_FUNC)
            dvelY = p_Du(2,icubp,iel,DER_FUNC)
            
            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Psi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Psi_i
              ! into dbasI
              dbasIx = p_DbasTest(idofe,DER_DERIV2D_X,icubp,iel)
              dbasIy = p_DbasTest(idofe,DER_DERIV2D_Y,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the (trial) basis function Phi_j:
              do jdofe=1,p_rmatrixData11%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJ = p_DbasTrial(jdofe,DER_FUNC,icubp,iel)
                
                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelX * dbasIx + &          ! ( phi , u1 psi_x )
                      dbasJ * dvelY * dbasIy )            ! ( phi , u2 psi_y )

                p_DlocalMatrix22(jdofe,idofe,iel) = p_DlocalMatrix22(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelX * dbasIx + &          ! ( phi , u1 psi_x )
                      dbasJ * dvelY * dbasIy )            ! ( phi , u2 psi_y )

              end do ! jdofe

            end do ! idofe

          end do ! icubp

        end do ! iel
        
      end if

    case (NDIM3D)

      ! Matrices to be set up
      p_rmatrixData11 => RmatrixData(iy,ix)
      p_rmatrixData22 => RmatrixData(iy+1,ix+1)
      p_rmatrixData33 => RmatrixData(iy+2,ix+2)

      p_DlocalMatrix11 => RmatrixData(iy,ix)%p_Dentry
      p_DlocalMatrix22 => RmatrixData(iy+1,ix+1)%p_Dentry
      p_DlocalMatrix33 => RmatrixData(iy+2,ix+2)%p_Dentry

      ! Currently, interleaved matrices are not supported
      if (p_rmatrixData11%bisInterleaved .or. &
          p_rmatrixData22%bisInterleaved .or. &
          p_rmatrixData33%bisInterleaved) then
        call output_line ("Interleaved matrices not supported",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_docalc_convection_vugradw")
        call sys_halt()
      end if

      if ((p_rmatrixData11%ndimfeTrial .ne. 1) .or. &
          (p_rmatrixData11%ndimfeTest .ne. 1)) then
        ! This does not make sense.
        call output_line ("Only scalar-valued FE spaces supported.",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_docalc_convection_vugradw")
        call sys_halt()
      end if

      if (bvelConst) then
      
        ! Set up the matrix for constant velocity.
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Psi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Psi_i
              ! into dbasI
              dbasIx = p_DbasTest(idofe,DER_DERIV3D_X,icubp,iel)
              dbasIy = p_DbasTest(idofe,DER_DERIV3D_Y,icubp,iel)
              dbasIz = p_DbasTest(idofe,DER_DERIV3D_Z,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the (trial) basis function Phi_j:
              do jdofe=1,p_rmatrixData11%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJ = p_DbasTrial(jdofe,DER_FUNC,icubp,iel)
                
                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelX * dbasIx + &          ! ( phi , u1 psi_x )
                      dbasJ * dvelY * dbasIy + &          ! ( phi , u2 psi_y )
                      dbasJ * dvelZ * dbasIz )            ! ( phi , u3 psi_z )

                p_DlocalMatrix22(jdofe,idofe,iel) = p_DlocalMatrix22(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelX * dbasIx + &          ! ( phi , u1 psi_x )
                      dbasJ * dvelY * dbasIy + &          ! ( phi , u2 psi_y )
                      dbasJ * dvelZ * dbasIz )            ! ( phi , u3 psi_z )

                p_DlocalMatrix33(jdofe,idofe,iel) = p_DlocalMatrix33(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelX * dbasIx + &          ! ( phi , u1 psi_x )
                      dbasJ * dvelY * dbasIy + &          ! ( phi , u2 psi_y )
                      dbasJ * dvelZ * dbasIz )            ! ( phi , u3 psi_z )

              end do ! jdofe

            end do ! idofe

          end do ! icubp

        end do ! iel
        
      else

        ! Set up the matrix for nonconstant velocity.
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Velocity field in this cubature point
            dvelX = p_Du(1,icubp,iel,DER_FUNC)
            dvelY = p_Du(2,icubp,iel,DER_FUNC)
            dvelZ = p_Du(3,icubp,iel,DER_FUNC)
            
            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Psi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Psi_i
              ! into dbasI
              dbasIx = p_DbasTest(idofe,DER_DERIV3D_X,icubp,iel)
              dbasIy = p_DbasTest(idofe,DER_DERIV3D_Y,icubp,iel)
              dbasIz = p_DbasTest(idofe,DER_DERIV3D_Z,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the (trial) basis function Phi_j:
              do jdofe=1,p_rmatrixData11%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJ = p_DbasTrial(jdofe,DER_FUNC,icubp,iel)
                
                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelX * dbasIx + &          ! ( phi , u1 psi_x )
                      dbasJ * dvelY * dbasIy + &          ! ( phi , u2 psi_y )
                      dbasJ * dvelZ * dbasIz )            ! ( phi , u3 psi_z )

                p_DlocalMatrix22(jdofe,idofe,iel) = p_DlocalMatrix22(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelX * dbasIx + &          ! ( phi , u1 psi_x )
                      dbasJ * dvelY * dbasIy + &          ! ( phi , u2 psi_y )
                      dbasJ * dvelZ * dbasIz )            ! ( phi , u3 psi_z )

                p_DlocalMatrix33(jdofe,idofe,iel) = p_DlocalMatrix33(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelX * dbasIx + &          ! ( phi , u1 psi_x )
                      dbasJ * dvelY * dbasIy + &          ! ( phi , u2 psi_y )
                      dbasJ * dvelZ * dbasIz )            ! ( phi , u3 psi_z )

              end do ! jdofe

            end do ! idofe

          end do ! icubp

        end do ! iel
        
      end if

    end select

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_convection_vugradw(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates a convection operator "( v, (u grad) w )" at position (x,y)
    ! of a block matrix, with a convection u given as a finite element function.
    ! v is the trial and w the test basis function.
    !
    ! Note: If rcollection is not specified, the matrix is calculated
    ! in all diagonal blocks with a multiplier of 1.
    ! If rcollection is specified, the following parameters are expected:
    !
    ! rcollection%DquickAccess(1) = multiplier in front of the operator.
    ! rcollection%IquickAccess(1) = x-position in the matrix where to set up the operator.
    ! rcollection%IquickAccess(2) = y-position in the matrix where to set up the operator.
    ! rcollection%IquickAccess(3) = 0, if the convection is a constant vector field.
    !                                  in this case:
    !                                  1D: rcollection%DquickAccess(2)   = x-velocity
    !                                  2D: rcollection%DquickAccess(2:3) = x/y-velocity
    !                                  3D: rcollection%DquickAccess(2:4) = x/y/z-velocity
    !                             = 1, if the convection is specified by a
    !                                  finite element velocity field. In this case,
    !                                  a finite element velocity field must be specified
    !                                  as parameter revalVectors to the call of 
    !                                  bma_buildMatrix.
!</description>

!<remarks>
    ! Remark 1: 
    ! Using the routines from feevaluation2, it is possible to specify
    ! a nonconstant velocity field. The corresponding code looks as follows:
    !
    ! <verb>
    !     use feevaluation2
    !
    !     ...
    !
    !     type(t_collection) :: rcollection
    !     type(t_scalarCubatureInfo) :: rcubatureInfo   ! Cubature formula
    !     type(t_matrixBlock) :: rmatrix                ! Matrix to be calculated
    !
    !     type(t_vectorBlock) :: rvelocity        ! The velocity field
    !     type(t_fev2Vectors) :: revalVectors     ! Collection of vectors to evaluate
    !
    !     ! Prepare the cubature formula
    !     call spdiscr_createDefCubStructure (..., rcubatureInfo, CUB_GEN_AUTO)
    !
    !     ...
    !
    !     rcollection%IquickAccess(1) = 1          ! x-Position in the matrix
    !     rcollection%IquickAccess(2) = 1          ! y-Position in the matrix
    !     rcollection%IquickAccess(3) = 1          ! Nonconstant viscosity
    !     rcollection%DquickAccess(1) = 1.0_DP     ! Scaling
    !
    !     ! Add the X-, Y- and Z-velocity to revalVectors
    !     call fev2_addVectorFieldToEvalList(revalVectors,0,&
    !         rvelocity%RvectorBlock(1),rvelocity%RvectorBlock(2),rvelocity%RvectorBlock(3))
    !
    !     ! Set up the matrix
    !     call bma_buildMatrix (rmatrix,BMA_CALC_STANDARD, bma_fcalc_convection_vugradw, &
    !         rcollection, revalVectors=revalVectors,rcubatureInfo=rcubatureInfo)
    !
    !     ! Release the vector structure
    !     call fev2_releaseVectorList(revalVectors)
    !
    !     ...
    !
    !     ! Release the cubature formula
    !     call spdiscr_releaseCubStructure (rcubatureInfo)
    ! </verb>
    !
    ! Remark 2:
    ! The routine currently assumes that all velocity components are discretised
    ! with the same FEM space.
    !
    ! Remark 3:
    ! The routine currently assumes that all velocity matrices are independent.
    ! Matrices sharing data are not supported. This cound be realised by
    ! taking care of the flags RmatrixData(:,:)%bsharedMatrixData which indicate
    ! which matrix data is shared.
    !
    ! Remark 4:
    ! Interleaved matrices are currently not supported.
!</remarks>

!<inputoutput>
    ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
    ! have to be filled with data.
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!</subroutine>

    ! Local variables
    integer :: ndim,ix,iy
    real(DP) :: dscale
    real(DP) :: dvelX, dvelY, dvelZ
    logical :: bvelConst

    ! Dimension of the underlying space
    ndim = rmatrixAssembly%p_rtriangulation%ndim
    
    ! Get parameters
    dscale = 1.0_DP
    dvelX = 0.0_DP
    dvelY = 0.0_DP
    dvelZ = 0.0_DP
    bvelConst = .true.
    ix = 1
    iy = 1
    
    if (present(rcollection)) then
      dscale = rcollection%DquickAccess(1)
      
      ! Constant velocity?
      bvelConst = (rcollection%IquickAccess(3) .eq. 0)
      
      if (bvelConst) then
        ! Get the constant velocity
        dvelX = rcollection%DquickAccess(2)
        if (ndim .ge. NDIM2D) dvelY = rcollection%DquickAccess(3)
        if (ndim .ge. NDIM3D) dvelZ = rcollection%DquickAccess(4)
      end if
      
      ! Position
      ix = rcollection%IquickAccess(1)
      iy = rcollection%IquickAccess(2)

    end if
    
    if (bvelConst) then
      call bma_docalc_convection_vugradw(RmatrixData,rassemblyData,rmatrixAssembly,&
          npointsPerElement,nelements,dscale,ix,iy,dvelX,dvelY,dvelZ)
    else
      call bma_docalc_convection_vugradw(RmatrixData,rassemblyData,rmatrixAssembly,&
          npointsPerElement,nelements,dscale,ix,iy,rvectorField=revalVectors%p_RvectorData(1))
    end if

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_docalc_convection_graduvw(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,dscale,ix,iy,&
      du1x,du1y,du1z,du2x,du2y,du2z,du3x,du3y,du3z,rvectorField)

!<description>  
    ! Calculates a convection operator "( (grad u) v, w)" at position (x,y)
    ! of a block matrix, with a convection u given as a finite element function
    ! and a scaling factor dscale in front.
    ! v is the trial and w the test basis function.
    !
    ! If rvectorField is not specified, a constant velocity derivative
    ! du1x,du1y,du1z,du2x,du2y,du2z,du3x,du3y,du3z
    ! is assumed, and (duIJ) must be present.
    ! If rvectorField is specified, it describes the underlying vector field.
!</description>

!</remarks>
    ! Remark 1:
    ! The routine currently assumes that all velocity components are discretised
    ! with the same FEM space.
    !
    ! Remark 2:
    ! The routine currently assumes that all velocity matrices are independent.
    ! Matrices sharing data are not supported. This cound be realised by
    ! taking care of the flags RmatrixData(:,:)%bsharedMatrixData which indicate
    ! which matrix data is shared.
    !
    ! Remark 3:
    ! Interleaved matrices are currently not supported.
!</remarks>

!<inputoutput>
    ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
    ! have to be filled with data.
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Scaling factor
    real(DP), intent(in) :: dscale
    
    ! Position in the global matrix where to assemble
    integer, intent(in) :: ix,iy
    
    ! OPTIONAL: Derivative of the velocity. Must be specified if rvectorField is not
    ! present. Ignored if rvectorField is present.
    real(DP), intent(in), optional :: du1x,du1y,du1z,du2x,du2y,du2z,du3x,du3y,du3z

    ! OPTIONAL: Evaluation structure that describes an underlying nonconstant
    ! velocity field. Can be omitted if du1/du2/du3 is given.
    type(t_fev2VectorData), intent(in), optional :: rvectorField
!</input>

!</subroutine>

    ! Local variables
    real(DP) :: dbasI, dbasJ
    integer :: iel, icubp, idofe, jdofe
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix11,p_DlocalMatrix12,p_DlocalMatrix13
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix21,p_DlocalMatrix22,p_DlocalMatrix23
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix31,p_DlocalMatrix32,p_DlocalMatrix33
    type(t_bmaMatrixData), pointer :: p_rmatrixData11,p_rmatrixData12,p_rmatrixData13
    type(t_bmaMatrixData), pointer :: p_rmatrixData21,p_rmatrixData22,p_rmatrixData23
    type(t_bmaMatrixData), pointer :: p_rmatrixData31,p_rmatrixData32,p_rmatrixData33
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial,p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:,:,:), pointer :: p_Du

    integer :: ndim
    real(DP) :: dvelX1, dvelX2, dvelX3
    real(DP) :: dvelY1, dvelY2, dvelY3
    real(DP) :: dvelZ1, dvelZ2, dvelZ3
    logical :: bvelConst

    ! Dimension of the underlying space
    ndim = rmatrixAssembly%p_rtriangulation%ndim
    
    ! Get parameters
    bvelConst = present(rvectorField)
    if (bvelConst) then
      dvelX1 = 0.0_DP
      dvelX2 = 0.0_DP
      dvelX3 = 0.0_DP
      
      dvelY1 = 0.0_DP
      dvelY2 = 0.0_DP
      dvelY3 = 0.0_DP
      
      dvelZ1 = 0.0_DP
      dvelZ2 = 0.0_DP
      dvelZ3 = 0.0_DP

      ! Constant velocity derivative
      if (present(du1x)) dvelX1 = du1x
      if (present(du1y)) dvelX2 = du1y
      if (present(du1z)) dvelX3 = du1z
      if (present(du2x)) dvelY1 = du2x
      if (present(du2y)) dvelY2 = du2y
      if (present(du2z)) dvelY3 = du2z
      if (present(du3x)) dvelZ1 = du3x
      if (present(du3y)) dvelZ2 = du3y
      if (present(du3z)) dvelZ3 = du3z
    else
      ! Get the velocity field from the parameters
      p_Du => rvectorField%p_DdataVec
    end if
    
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Get local data
    p_DbasTrial => RmatrixData(iy,ix)%p_DbasTrial
    p_DbasTest => RmatrixData(iy,ix)%p_DbasTest

    ! Set up the local matrix of the convection.
    select case (ndim)
    case (NDIM1D)
      ! Matrices to be set up
      p_rmatrixData11 => RmatrixData(iy,ix)
      
      p_DlocalMatrix11 => RmatrixData(iy,ix)%p_Dentry
      
      ! Currently, interleaved matrices are not supported
      if (p_rmatrixData11%bisInterleaved) then
        call output_line ("Interleaved matrices not supported",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_docalc_convection_graduvw")
        call sys_halt()
      end if

      if ((p_rmatrixData11%ndimfeTrial .ne. 1) .or. &
          (p_rmatrixData11%ndimfeTest .ne. 1)) then
        ! This does not make sense.
        call output_line ("Only scalar-valued FE spaces supported.",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_docalc_convection_graduvw")
        call sys_halt()
      end if

      if (bvelConst) then
      
        ! Set up the matrix for constant velocity.
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Psi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Psi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the (trial) basis function Phi_j:
              do jdofe=1,p_rmatrixData11%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJ = p_DbasTrial(jdofe,DER_FUNC,icubp,iel)
                
                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX1 * dbasJ * dbasI )            ! ( u1_x phi , psi )

              end do ! jdofe

            end do ! idofe

          end do ! icubp

        end do ! iel
        
      else

        ! Set up the matrix for nonconstant velocity.
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Velocity field in this cubature point
            dvelX1 = p_Du(1,icubp,iel,DER_DERIV1D_X)
            
            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Psi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Psi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the (trial) basis function Phi_j:
              do jdofe=1,p_rmatrixData11%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJ = p_DbasTrial(jdofe,DER_DERIV1D_X,icubp,iel)
                
                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX1 * dbasJ * dbasI )            ! ( u1_x phi , psi )

              end do ! jdofe

            end do ! idofe

          end do ! icubp

        end do ! iel
        
      end if
      
    case (NDIM2D)

      ! Matrices to be set up
      p_rmatrixData11 => RmatrixData(iy,ix)
      p_rmatrixData12 => RmatrixData(iy,ix+1)
      p_rmatrixData21 => RmatrixData(iy+1,ix)
      p_rmatrixData22 => RmatrixData(iy+1,ix+1)
      
      p_DlocalMatrix11 => RmatrixData(iy,ix)%p_Dentry
      p_DlocalMatrix12 => RmatrixData(iy,ix+1)%p_Dentry
      p_DlocalMatrix21 => RmatrixData(iy+1,ix)%p_Dentry
      p_DlocalMatrix22 => RmatrixData(iy+1,ix+1)%p_Dentry

      ! Currently, interleaved matrices are not supported
      if (p_rmatrixData11%bisInterleaved .or. &
          p_rmatrixData12%bisInterleaved .or. &
          p_rmatrixData21%bisInterleaved .or. &
          p_rmatrixData22%bisInterleaved) then
        call output_line ("Interleaved matrices not supported",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_docalc_convection_graduvw")
        call sys_halt()
      end if

      if ((p_rmatrixData11%ndimfeTrial .ne. 1) .or. &
          (p_rmatrixData11%ndimfeTest .ne. 1)) then
        ! This does not make sense.
        call output_line ("Only scalar-valued FE spaces supported.",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_docalc_convection_graduvw")
        call sys_halt()
      end if

      if (bvelConst) then
      
        ! Set up the matrix for constant velocity.
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Psi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Psi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the (trial) basis function Phi_j:
              do jdofe=1,p_rmatrixData11%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJ = p_DbasTrial(jdofe,DER_FUNC,icubp,iel)
                
                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX1 * dbasJ * dbasI )            ! ( u1_x phi , psi )

                p_DlocalMatrix12(jdofe,idofe,iel) = p_DlocalMatrix12(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX2 * dbasJ * dbasI )            ! ( u1_y phi , psi )

                p_DlocalMatrix21(jdofe,idofe,iel) = p_DlocalMatrix21(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelY1 * dbasJ * dbasI )            ! ( u2_x phi , psi )

                p_DlocalMatrix22(jdofe,idofe,iel) = p_DlocalMatrix22(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelY2 * dbasJ * dbasI )            ! ( u2_y phi , psi )

              end do ! jdofe

            end do ! idofe

          end do ! icubp

        end do ! iel
        
      else

        ! Set up the matrix for nonconstant velocity.
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Velocity field in this cubature point
            dvelX1 = p_Du(1,icubp,iel,DER_DERIV2D_X)
            dvelY1 = p_Du(2,icubp,iel,DER_DERIV2D_X)
            dvelX2 = p_Du(1,icubp,iel,DER_DERIV2D_Y)
            dvelY2 = p_Du(2,icubp,iel,DER_DERIV2D_Y)
            
            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Psi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Psi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the (trial) basis function Phi_j:
              do jdofe=1,p_rmatrixData11%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJ = p_DbasTrial(jdofe,DER_FUNC,icubp,iel)
                
                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX1 * dbasJ * dbasI )            ! ( u1_x phi , psi )

                p_DlocalMatrix12(jdofe,idofe,iel) = p_DlocalMatrix12(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX2 * dbasJ * dbasI )            ! ( u1_y phi , psi )

                p_DlocalMatrix21(jdofe,idofe,iel) = p_DlocalMatrix21(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelY1 * dbasJ * dbasI )            ! ( u2_x phi , psi )

                p_DlocalMatrix22(jdofe,idofe,iel) = p_DlocalMatrix22(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelY2 * dbasJ * dbasI )            ! ( u2_y phi , psi )

              end do ! jdofe

            end do ! idofe

          end do ! icubp

        end do ! iel
        
      end if

    case (NDIM3D)

      ! Matrices to be set up
      p_rmatrixData11 => RmatrixData(iy,ix)
      p_rmatrixData12 => RmatrixData(iy,ix+1)
      p_rmatrixData13 => RmatrixData(iy,ix+2)
      p_rmatrixData21 => RmatrixData(iy+1,ix)
      p_rmatrixData22 => RmatrixData(iy+1,ix+1)
      p_rmatrixData23 => RmatrixData(iy+1,ix+2)
      p_rmatrixData31 => RmatrixData(iy+2,ix)
      p_rmatrixData32 => RmatrixData(iy+2,ix+1)
      p_rmatrixData33 => RmatrixData(iy+2,ix+2)
      
      p_DlocalMatrix11 => RmatrixData(iy,ix)%p_Dentry
      p_DlocalMatrix12 => RmatrixData(iy,ix+1)%p_Dentry
      p_DlocalMatrix13 => RmatrixData(iy,ix+3)%p_Dentry
      p_DlocalMatrix21 => RmatrixData(iy+1,ix)%p_Dentry
      p_DlocalMatrix22 => RmatrixData(iy+1,ix+1)%p_Dentry
      p_DlocalMatrix23 => RmatrixData(iy+1,ix+3)%p_Dentry
      p_DlocalMatrix31 => RmatrixData(iy+2,ix)%p_Dentry
      p_DlocalMatrix32 => RmatrixData(iy+2,ix+1)%p_Dentry
      p_DlocalMatrix33 => RmatrixData(iy+2,ix+3)%p_Dentry

      ! Currently, interleaved matrices are not supported
      if ( p_rmatrixData11%bisInterleaved .or. &
           p_rmatrixData12%bisInterleaved .or. &
           p_rmatrixData13%bisInterleaved .or. &
           p_rmatrixData21%bisInterleaved .or. &
           p_rmatrixData22%bisInterleaved .or. &
           p_rmatrixData23%bisInterleaved .or. &
           p_rmatrixData31%bisInterleaved .or. &
           p_rmatrixData32%bisInterleaved .or. &
           p_rmatrixData33%bisInterleaved) then
        call output_line ("Interleaved matrices not supported",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_docalc_convection_graduvw")
        call sys_halt()
      end if

      if ((p_rmatrixData11%ndimfeTrial .ne. 1) .or. &
          (p_rmatrixData11%ndimfeTest .ne. 1)) then
        ! This does not make sense.
        call output_line ("Only scalar-valued FE spaces supported.",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_docalc_convection_graduvw")
        call sys_halt()
      end if

      if (bvelConst) then
      
        ! Set up the matrix for constant velocity.
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Psi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Psi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the (trial) basis function Phi_j:
              do jdofe=1,p_rmatrixData11%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJ = p_DbasTrial(jdofe,DER_DERIV3D_X,icubp,iel)
                
                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX1 * dbasJ * dbasI )            ! ( u1_x phi , psi )

                p_DlocalMatrix12(jdofe,idofe,iel) = p_DlocalMatrix12(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX2 * dbasJ * dbasI )            ! ( u1_y phi , psi )

                p_DlocalMatrix13(jdofe,idofe,iel) = p_DlocalMatrix13(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX3 * dbasJ * dbasI )            ! ( u1_y phi , psi )

                p_DlocalMatrix21(jdofe,idofe,iel) = p_DlocalMatrix21(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelY1 * dbasJ * dbasI )            ! ( u2_x phi , psi )

                p_DlocalMatrix22(jdofe,idofe,iel) = p_DlocalMatrix22(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelY2 * dbasJ * dbasI )            ! ( u2_y phi , psi )

                p_DlocalMatrix23(jdofe,idofe,iel) = p_DlocalMatrix23(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelY3 * dbasJ * dbasI )            ! ( u2_y phi , psi )

                p_DlocalMatrix31(jdofe,idofe,iel) = p_DlocalMatrix31(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelZ1 * dbasJ * dbasI )            ! ( u3_x phi , psi )

                p_DlocalMatrix32(jdofe,idofe,iel) = p_DlocalMatrix32(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelZ2 * dbasJ * dbasI )            ! ( u3_y phi , psi )

                p_DlocalMatrix33(jdofe,idofe,iel) = p_DlocalMatrix33(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelZ3 * dbasJ * dbasI )            ! ( u3_y phi , psi )

              end do ! jdofe

            end do ! idofe

          end do ! icubp

        end do ! iel
        
      else

        ! Set up the matrix for nonconstant velocity.
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Velocity field in this cubature point
            dvelX1 = p_Du(1,icubp,iel,DER_DERIV3D_X)
            dvelY1 = p_Du(2,icubp,iel,DER_DERIV3D_X)
            dvelZ1 = p_Du(3,icubp,iel,DER_DERIV3D_X)

            dvelX2 = p_Du(1,icubp,iel,DER_DERIV3D_Y)
            dvelY2 = p_Du(2,icubp,iel,DER_DERIV3D_Y)
            dvelZ2 = p_Du(3,icubp,iel,DER_DERIV3D_Y)

            dvelX3 = p_Du(1,icubp,iel,DER_DERIV3D_Z)
            dvelY3 = p_Du(2,icubp,iel,DER_DERIV3D_Z)
            dvelZ3 = p_Du(3,icubp,iel,DER_DERIV3D_Z)
            
            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Psi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Psi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the (trial) basis function Phi_j:
              do jdofe=1,p_rmatrixData11%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJ = p_DbasTrial(jdofe,DER_DERIV3D_X,icubp,iel)
                
                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX1 * dbasJ * dbasI )            ! ( u1_x phi , psi )

                p_DlocalMatrix12(jdofe,idofe,iel) = p_DlocalMatrix12(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX2 * dbasJ * dbasI )            ! ( u1_y phi , psi )

                p_DlocalMatrix13(jdofe,idofe,iel) = p_DlocalMatrix13(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelX3 * dbasJ * dbasI )            ! ( u1_y phi , psi )

                p_DlocalMatrix21(jdofe,idofe,iel) = p_DlocalMatrix21(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelY1 * dbasJ * dbasI )            ! ( u2_x phi , psi )

                p_DlocalMatrix22(jdofe,idofe,iel) = p_DlocalMatrix22(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelY2 * dbasJ * dbasI )            ! ( u2_y phi , psi )

                p_DlocalMatrix23(jdofe,idofe,iel) = p_DlocalMatrix23(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelY3 * dbasJ * dbasI )            ! ( u2_y phi , psi )

                p_DlocalMatrix31(jdofe,idofe,iel) = p_DlocalMatrix31(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelZ1 * dbasJ * dbasI )            ! ( u3_x phi , psi )

                p_DlocalMatrix32(jdofe,idofe,iel) = p_DlocalMatrix32(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelZ2 * dbasJ * dbasI )            ! ( u3_y phi , psi )

                p_DlocalMatrix33(jdofe,idofe,iel) = p_DlocalMatrix33(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dvelZ3 * dbasJ * dbasI )            ! ( u3_y phi , psi )

              end do ! jdofe

            end do ! idofe

          end do ! icubp

        end do ! iel
        
      end if

    end select

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_convection_graduvw(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates a convection operator "( (grad u) v, w)" at position (x,y)
    ! of a block matrix, with a convection u given as a finite element function.
    ! v is the trial and w the test basis function.
    !
    ! Note: If rcollection is not specified, the matrix is calculated
    ! in all diagonal blocks with a multiplier of 1.
    ! If rcollection is specified, the following parameters are expected:
    !
    ! rcollection%DquickAccess(1) = multiplier in front of the operator.
    ! rcollection%IquickAccess(1) = x-position in the matrix where to set up the operator.
    ! rcollection%IquickAccess(2) = y-position in the matrix where to set up the operator.
    ! rcollection%IquickAccess(3) = 0, if the convection is a constant vector field.
    !                                  in this case:
    !                                  1D: rcollection%DquickAccess(2)   = d/dx x-velocity
    !                                  2D: rcollection%DquickAccess(2:3) = (d/dx,d/dy) x-velocity
    !                                      rcollection%DquickAccess(4:5) = (d/dx,d/dy) y-velocity
    !                                  3D: rcollection%DquickAccess(2:4) = (d/dx,d/dy,d/dz) x-velocity
    !                                      rcollection%DquickAccess(5:7) = (d/dx,d/dy,d/dz) y-velocity
    !                                      rcollection%DquickAccess(8:10) = (d/dx,d/dy,d/dz) z-velocity
    !                             = 1, if The convection is specified by a
    !                                  finite element velocity field. In this case,
    !                                  a finite element velocity field must be specified
    !                                  as parameter revalVectors to the call of 
    !                                  bma_buildMatrix.
!</description>

!<remarks>
    ! Remark 1: 
    ! Using the routines from feevaluation2, it is possible to specify
    ! a nonconstant velocity field. The corresponding code looks as follows:
    !
    ! <verb>
    !     use feevaluation2
    !
    !     ...
    !
    !     type(t_collection) :: rcollection
    !     type(t_scalarCubatureInfo) :: rcubatureInfo   ! Cubature formula
    !     type(t_matrixBlock) :: rmatrix                ! Matrix to be calculated
    !
    !     type(t_vectorBlock) :: rvelocity        ! The velocity field
    !     type(t_fev2Vectors) :: revalVectors     ! Collection of vectors to evaluate
    !
    !     ! Prepare the cubature formula
    !     call spdiscr_createDefCubStructure (..., rcubatureInfo, CUB_GEN_AUTO)
    !
    !     ...
    !
    !     rcollection%IquickAccess(1) = 1          ! x-Position in the matrix
    !     rcollection%IquickAccess(2) = 1          ! y-Position in the matrix
    !     rcollection%IquickAccess(3) = 1          ! Nonconstant viscosity
    !     rcollection%DquickAccess(1) = 1.0_DP     ! Scaling
    !
    !     ! Add the X-, Y- and Z-velocity to revalVectors
    !     call fev2_addVectorFieldToEvalList(revalVectors,1,&
    !         rvelocity%RvectorBlock(1),rvelocity%RvectorBlock(2),rvelocity%RvectorBlock(3))
    !
    !     ! Set up the matrix
    !     call bma_buildMatrix (rmatrix,BMA_CALC_STANDARD, bma_fcalc_convection_graduvw, &
    !         rcollection, revalVectors=revalVectors,rcubatureInfo=rcubatureInfo)
    !
    !     ! Release the vector structure
    !     call fev2_releaseVectorList(revalVectors)
    !
    !     ...
    !
    !     ! Release the cubature formula
    !     call spdiscr_releaseCubStructure (rcubatureInfo)
    ! </verb>
    !
    ! Remark 2:
    ! The routine currently assumes that all velocity components are discretised
    ! with the same FEM space.
    !
    ! Remark 3:
    ! The routine currently assumes that all velocity matrices are independent.
    ! Matrices sharing data are not supported. This cound be realised by
    ! taking care of the flags RmatrixData(:,:)%bsharedMatrixData which indicate
    ! which matrix data is shared.
    !
    ! Remark 4:
    ! Interleaved matrices are currently not supported.
!</remarks>

!<inputoutput>
    ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
    ! have to be filled with data.
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!</subroutine>

    ! Local variables
    integer :: ndim,ix,iy
    real(DP) :: dscale
    real(DP) :: dvelX1, dvelX2, dvelX3
    real(DP) :: dvelY1, dvelY2, dvelY3
    real(DP) :: dvelZ1, dvelZ2, dvelZ3
    logical :: bvelConst

    ! Dimension of the underlying space
    ndim = rmatrixAssembly%p_rtriangulation%ndim
    
    ! Get parameters
    dscale = 1.0_DP
    dvelX1 = 0.0_DP
    dvelX2 = 0.0_DP
    dvelX3 = 0.0_DP
    
    dvelY1 = 0.0_DP
    dvelY2 = 0.0_DP
    dvelY3 = 0.0_DP
    
    dvelZ1 = 0.0_DP
    dvelZ2 = 0.0_DP
    dvelZ3 = 0.0_DP
    
    bvelConst = .true.
    ix = 1
    iy = 1
    
    if (present(rcollection)) then
      dscale = rcollection%DquickAccess(1)
      
      ! Constant velocity?
      bvelConst = (rcollection%IquickAccess(3) .eq. 0)
      
      if (bvelConst) then
        ! Get the constant velocity
        dvelX1 = rcollection%DquickAccess(2)
        if (ndim .ge. NDIM2D) then
          dvelX2 = rcollection%DquickAccess(3)
          
          dvelY1 = rcollection%DquickAccess(4)
          dvelY2 = rcollection%DquickAccess(5)
        end if
        if (ndim .ge. NDIM3D) then
          dvelX2 = rcollection%DquickAccess(3)
          dvelX3 = rcollection%DquickAccess(4)
          
          dvelY1 = rcollection%DquickAccess(5)
          dvelY2 = rcollection%DquickAccess(6)
          dvelY3 = rcollection%DquickAccess(7)
          
          dvelZ1 = rcollection%DquickAccess(8)
          dvelZ2 = rcollection%DquickAccess(9)
          dvelZ3 = rcollection%DquickAccess(10)
        end if
      end if
      
      ! Position
      ix = rcollection%IquickAccess(1)
      iy = rcollection%IquickAccess(2)

    end if
    
    if (bvelConst) then
      call bma_docalc_convection_graduvw(RmatrixData,rassemblyData,rmatrixAssembly,&
          npointsPerElement,nelements,dscale,ix,iy,&
          dvelX1,dvelX2,dvelX3,dvelY1,dvelY2,dvelY3,dvelZ1,dvelZ2,dvelZ3)
    else
      call bma_docalc_convection_graduvw(RmatrixData,rassemblyData,rmatrixAssembly,&
          npointsPerElement,nelements,dscale,ix,iy,rvectorField=revalVectors%p_RvectorData(1))
    end if


  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_docalc_convection_vgraduw(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,dscale,ix,iy,&
      du1x,du1y,du1z,du2x,du2y,du2z,du3x,du3y,du3z,rvectorField)

!<description>  
    ! Calculates a convection operator "( v, (grad u) w)" at position (x,y)
    ! of a block matrix, with a convection u given as a finite element function
    ! and a scaling factor dscale in front.
    ! v is the trial and w the test basis function.
    !
    ! If rvectorField is not specified, a constant velocity derivative
    ! du1x,du1y,du1z,du2x,du2y,du2z,du3x,du3y,du3z
    ! is assumed, and (duIJ) must be present.
    ! If rvectorField is specified, it describes the underlying vector field.
!</description>

!</remarks>
    ! Remark 1:
    ! The routine currently assumes that all velocity components are discretised
    ! with the same FEM space.
    !
    ! Remark 2:
    ! The routine currently assumes that all velocity matrices are independent.
    ! Matrices sharing data are not supported. This cound be realised by
    ! taking care of the flags RmatrixData(:,:)%bsharedMatrixData which indicate
    ! which matrix data is shared.
    !
    ! Remark 3:
    ! Interleaved matrices are currently not supported.
!</remarks>

!<inputoutput>
    ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
    ! have to be filled with data.
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Scaling factor
    real(DP), intent(in) :: dscale
    
    ! Position in the global matrix where to assemble
    integer, intent(in) :: ix,iy
    
    ! OPTIONAL: Derivative of the velocity. Must be specified if rvectorField is not
    ! present. Ignored if rvectorField is present.
    real(DP), intent(in), optional :: du1x,du1y,du1z,du2x,du2y,du2z,du3x,du3y,du3z

    ! OPTIONAL: Evaluation structure that describes an underlying nonconstant
    ! velocity field. Can be omitted if du1/du2/du3 is given.
    type(t_fev2VectorData), intent(in), optional :: rvectorField
!</input>

!</subroutine>

    ! Local variables
    real(DP) :: dbasI, dbasJ
    integer :: iel, icubp, idofe, jdofe
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix11,p_DlocalMatrix12,p_DlocalMatrix13
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix21,p_DlocalMatrix22,p_DlocalMatrix23
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix31,p_DlocalMatrix32,p_DlocalMatrix33
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial,p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaMatrixData), pointer :: p_rmatrixData11,p_rmatrixData12,p_rmatrixData13
    type(t_bmaMatrixData), pointer :: p_rmatrixData21,p_rmatrixData22,p_rmatrixData23
    type(t_bmaMatrixData), pointer :: p_rmatrixData31,p_rmatrixData32,p_rmatrixData33
    real(DP), dimension(:,:,:,:), pointer :: p_Du
    
    integer :: ndim
    real(DP) :: dvelX1, dvelX2, dvelX3
    real(DP) :: dvelY1, dvelY2, dvelY3
    real(DP) :: dvelZ1, dvelZ2, dvelZ3
    logical :: bvelConst

    ! Dimension of the underlying space
    ndim = rmatrixAssembly%p_rtriangulation%ndim
    
    ! Get parameters
    bvelConst = present(rvectorField)
    if (bvelConst) then
      dvelX1 = 0.0_DP
      dvelX2 = 0.0_DP
      dvelX3 = 0.0_DP
      
      dvelY1 = 0.0_DP
      dvelY2 = 0.0_DP
      dvelY3 = 0.0_DP
      
      dvelZ1 = 0.0_DP
      dvelZ2 = 0.0_DP
      dvelZ3 = 0.0_DP

      ! Constant velocity derivative
      if (present(du1x)) dvelX1 = du1x
      if (present(du1y)) dvelX2 = du1y
      if (present(du1z)) dvelX3 = du1z
      if (present(du2x)) dvelY1 = du2x
      if (present(du2y)) dvelY2 = du2y
      if (present(du2z)) dvelY3 = du2z
      if (present(du3x)) dvelZ1 = du3x
      if (present(du3y)) dvelZ2 = du3y
      if (present(du3z)) dvelZ3 = du3z
    else
      ! Get the velocity field from the parameters
      p_Du => rvectorField%p_DdataVec
    end if
    
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Get local data
    p_DbasTrial => RmatrixData(iy,ix)%p_DbasTrial
    p_DbasTest => RmatrixData(iy,ix)%p_DbasTest

    ! Set up the local matrix of the convection.
    select case (ndim)
    case (NDIM1D)
      ! Matrices to be set up
      p_rmatrixData11 => RmatrixData(iy,ix)
      
      p_DlocalMatrix11 => RmatrixData(iy,ix)%p_Dentry
      
      ! Currently, interleaved matrices are not supported
      if (p_rmatrixData11%bisInterleaved) then
        call output_line ("Interleaved matrices not supported",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_docalc_convection_vgraduw")
        call sys_halt()
      end if

      if ((p_rmatrixData11%ndimfeTrial .ne. 1) .or. &
          (p_rmatrixData11%ndimfeTest .ne. 1)) then
        ! This does not make sense.
        call output_line ("Only scalar-valued FE spaces supported.",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_docalc_convection_vgraduw")
        call sys_halt()
      end if

      if (bvelConst) then
      
        ! Set up the matrix for constant velocity.
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Psi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Psi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the (trial) basis function Phi_j:
              do jdofe=1,p_rmatrixData11%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJ = p_DbasTrial(jdofe,DER_FUNC,icubp,iel)
                
                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelX1 * dbasI )            ! ( phi , u1_x psi )

              end do ! jdofe

            end do ! idofe

          end do ! icubp

        end do ! iel
        
      else

        ! Set up the matrix for nonconstant velocity.
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Velocity field in this cubature point
            dvelX1 = p_Du(1,icubp,iel,DER_DERIV1D_X)
            
            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Psi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Psi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the (trial) basis function Phi_j:
              do jdofe=1,p_rmatrixData11%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJ = p_DbasTrial(jdofe,DER_DERIV1D_X,icubp,iel)
                
                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelX1 * dbasI )            ! ( phi , u1_x psi )

              end do ! jdofe

            end do ! idofe

          end do ! icubp

        end do ! iel
        
      end if
      
    case (NDIM2D)

      ! Matrices to be set up
      p_rmatrixData11 => RmatrixData(iy,ix)
      p_rmatrixData12 => RmatrixData(iy,ix+1)
      p_rmatrixData21 => RmatrixData(iy+1,ix)
      p_rmatrixData22 => RmatrixData(iy+1,ix+1)
      
      p_DlocalMatrix11 => RmatrixData(iy,ix)%p_Dentry
      p_DlocalMatrix12 => RmatrixData(iy,ix+1)%p_Dentry
      p_DlocalMatrix21 => RmatrixData(iy+1,ix)%p_Dentry
      p_DlocalMatrix22 => RmatrixData(iy+1,ix+1)%p_Dentry

      ! Currently, interleaved matrices are not supported
      if (p_rmatrixData11%bisInterleaved .or. &
          p_rmatrixData12%bisInterleaved .or. &
          p_rmatrixData21%bisInterleaved .or. &
          p_rmatrixData22%bisInterleaved) then
        call output_line ("Interleaved matrices not supported",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_docalc_convection_vgraduw")
        call sys_halt()
      end if

      if ((p_rmatrixData11%ndimfeTrial .ne. 1) .or. &
          (p_rmatrixData11%ndimfeTest .ne. 1)) then
        ! This does not make sense.
        call output_line ("Only scalar-valued FE spaces supported.",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_docalc_convection_vgraduw")
        call sys_halt()
      end if

      if (bvelConst) then
      
        ! Set up the matrix for constant velocity.
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Psi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Psi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the (trial) basis function Phi_j:
              do jdofe=1,p_rmatrixData11%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJ = p_DbasTrial(jdofe,DER_FUNC,icubp,iel)
                
                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelX1 * dbasI )            ! ( phi , u1_x psi )

                p_DlocalMatrix12(jdofe,idofe,iel) = p_DlocalMatrix12(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelY1 * dbasI )            ! ( phi , u2_x psi )

                p_DlocalMatrix21(jdofe,idofe,iel) = p_DlocalMatrix21(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelX2 * dbasI )            ! ( phi , u1_y psi )

                p_DlocalMatrix22(jdofe,idofe,iel) = p_DlocalMatrix22(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelY2 * dbasI )            ! ( phi , u2_y psi )

              end do ! jdofe

            end do ! idofe

          end do ! icubp

        end do ! iel
        
      else

        ! Set up the matrix for nonconstant velocity.
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Velocity field in this cubature point
            dvelX1 = p_Du(1,icubp,iel,DER_DERIV2D_X)
            dvelY1 = p_Du(2,icubp,iel,DER_DERIV2D_X)

            dvelX2 = p_Du(1,icubp,iel,DER_DERIV2D_Y)
            dvelY2 = p_Du(2,icubp,iel,DER_DERIV2D_Y)
            
            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Psi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Psi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the (trial) basis function Phi_j:
              do jdofe=1,p_rmatrixData11%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJ = p_DbasTrial(jdofe,DER_FUNC,icubp,iel)
                
                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelX1 * dbasI )            ! ( phi , u1_x psi )

                p_DlocalMatrix12(jdofe,idofe,iel) = p_DlocalMatrix12(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelY1 * dbasI )            ! ( phi , u2_x psi )

                p_DlocalMatrix21(jdofe,idofe,iel) = p_DlocalMatrix21(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelX2 * dbasI )            ! ( phi , u1_y psi )

                p_DlocalMatrix22(jdofe,idofe,iel) = p_DlocalMatrix22(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelY2 * dbasI )            ! ( phi , u2_y psi )

              end do ! jdofe

            end do ! idofe

          end do ! icubp

        end do ! iel
        
      end if

    case (NDIM3D)

      ! Matrices to be set up
      p_rmatrixData11 => RmatrixData(iy,ix)
      p_rmatrixData12 => RmatrixData(iy,ix+1)
      p_rmatrixData13 => RmatrixData(iy,ix+2)
      p_rmatrixData21 => RmatrixData(iy+1,ix)
      p_rmatrixData22 => RmatrixData(iy+1,ix+1)
      p_rmatrixData23 => RmatrixData(iy+1,ix+2)
      p_rmatrixData31 => RmatrixData(iy+2,ix)
      p_rmatrixData32 => RmatrixData(iy+2,ix+1)
      p_rmatrixData33 => RmatrixData(iy+2,ix+2)
      
      p_DlocalMatrix11 => RmatrixData(iy,ix)%p_Dentry
      p_DlocalMatrix12 => RmatrixData(iy,ix+1)%p_Dentry
      p_DlocalMatrix13 => RmatrixData(iy,ix+3)%p_Dentry
      p_DlocalMatrix21 => RmatrixData(iy+1,ix)%p_Dentry
      p_DlocalMatrix22 => RmatrixData(iy+1,ix+1)%p_Dentry
      p_DlocalMatrix23 => RmatrixData(iy+1,ix+3)%p_Dentry
      p_DlocalMatrix31 => RmatrixData(iy+2,ix)%p_Dentry
      p_DlocalMatrix32 => RmatrixData(iy+2,ix+1)%p_Dentry
      p_DlocalMatrix33 => RmatrixData(iy+2,ix+3)%p_Dentry

      ! Currently, interleaved matrices are not supported
      if (p_rmatrixData11%bisInterleaved .or. &
          p_rmatrixData12%bisInterleaved .or. &
          p_rmatrixData13%bisInterleaved .or. &
          p_rmatrixData21%bisInterleaved .or. &
          p_rmatrixData22%bisInterleaved .or. &
          p_rmatrixData23%bisInterleaved .or. &
          p_rmatrixData31%bisInterleaved .or. &
          p_rmatrixData32%bisInterleaved .or. &
          p_rmatrixData33%bisInterleaved) then
        call output_line ("Interleaved matrices not supported",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_docalc_convection_vgraduw")
        call sys_halt()
      end if

      if ((p_rmatrixData11%ndimfeTrial .ne. 1) .or. &
          (p_rmatrixData11%ndimfeTest .ne. 1)) then
        ! This does not make sense.
        call output_line ("Only scalar-valued FE spaces supported.",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_docalc_convection_vgraduw")
        call sys_halt()
      end if

      if (bvelConst) then
      
        ! Set up the matrix for constant velocity.
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Psi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Psi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the (trial) basis function Phi_j:
              do jdofe=1,p_rmatrixData11%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJ = p_DbasTrial(jdofe,DER_DERIV3D_X,icubp,iel)
                
                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelX1 * dbasI )            ! ( phi , u1_x psi )

                p_DlocalMatrix12(jdofe,idofe,iel) = p_DlocalMatrix12(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelY1 * dbasI )            ! ( phi , u2_x psi )

                p_DlocalMatrix13(jdofe,idofe,iel) = p_DlocalMatrix13(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelZ1 * dbasI )            ! ( phi , u3_x psi )

                p_DlocalMatrix21(jdofe,idofe,iel) = p_DlocalMatrix21(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelX2 * dbasI )            ! ( phi , u1_y psi )

                p_DlocalMatrix22(jdofe,idofe,iel) = p_DlocalMatrix22(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelY2 * dbasI )            ! ( phi , u2_y psi )

                p_DlocalMatrix23(jdofe,idofe,iel) = p_DlocalMatrix23(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelZ2 * dbasI )            ! ( phi , u3_y psi )

                p_DlocalMatrix31(jdofe,idofe,iel) = p_DlocalMatrix31(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelX3 * dbasI )            ! ( phi , u1_z psi )

                p_DlocalMatrix32(jdofe,idofe,iel) = p_DlocalMatrix32(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelY3 * dbasI )            ! ( phi , u2_z psi )

                p_DlocalMatrix33(jdofe,idofe,iel) = p_DlocalMatrix33(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelZ3 * dbasI )            ! ( phi , u3_z psi )

              end do ! jdofe

            end do ! idofe

          end do ! icubp

        end do ! iel
        
      else

        ! Set up the matrix for nonconstant velocity.
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Velocity field in this cubature point
            dvelX1 = p_Du(1,icubp,iel,DER_DERIV3D_X)
            dvelX2 = p_Du(1,icubp,iel,DER_DERIV3D_Y)
            dvelZ1 = p_Du(3,icubp,iel,DER_DERIV3D_X)

            dvelY1 = p_Du(2,icubp,iel,DER_DERIV3D_X)
            dvelY2 = p_Du(2,icubp,iel,DER_DERIV3D_Y)
            dvelZ2 = p_Du(3,icubp,iel,DER_DERIV3D_Y)

            dvelX3 = p_Du(1,icubp,iel,DER_DERIV3D_Z)
            dvelY3 = p_Du(2,icubp,iel,DER_DERIV3D_Z)
            dvelZ3 = p_Du(3,icubp,iel,DER_DERIV3D_Z)
            
            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Psi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Psi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the (trial) basis function Phi_j:
              do jdofe=1,p_rmatrixData11%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJ = p_DbasTrial(jdofe,DER_DERIV3D_X,icubp,iel)
                
                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix11(jdofe,idofe,iel) = p_DlocalMatrix11(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelX1 * dbasI )            ! ( phi , u1_x psi )

                p_DlocalMatrix12(jdofe,idofe,iel) = p_DlocalMatrix12(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelY1 * dbasI )            ! ( phi , u2_x psi )

                p_DlocalMatrix13(jdofe,idofe,iel) = p_DlocalMatrix13(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelZ1 * dbasI )            ! ( phi , u3_x psi )

                p_DlocalMatrix21(jdofe,idofe,iel) = p_DlocalMatrix21(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelX2 * dbasI )            ! ( phi , u1_y psi )

                p_DlocalMatrix22(jdofe,idofe,iel) = p_DlocalMatrix22(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelY2 * dbasI )            ! ( phi , u2_y psi )

                p_DlocalMatrix23(jdofe,idofe,iel) = p_DlocalMatrix23(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelZ2 * dbasI )            ! ( phi , u3_y psi )

                p_DlocalMatrix31(jdofe,idofe,iel) = p_DlocalMatrix31(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelX3 * dbasI )            ! ( phi , u1_z psi )

                p_DlocalMatrix32(jdofe,idofe,iel) = p_DlocalMatrix32(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelY3 * dbasI )            ! ( phi , u2_z psi )

                p_DlocalMatrix33(jdofe,idofe,iel) = p_DlocalMatrix33(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * &
                    ( dbasJ * dvelZ3 * dbasI )            ! ( phi , u3_z psi )

              end do ! jdofe

            end do ! idofe

          end do ! icubp

        end do ! iel
        
      end if

    end select

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_convection_vgraduw(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates a convection operator "( v, (grad u) w)" at position (x,y)
    ! of a block matrix, with a convection u given as a finite element function.
    ! v is the trial and w the test basis function.
    !
    ! Note: If rcollection is not specified, the matrix is calculated
    ! in all diagonal blocks with a multiplier of 1.
    ! If rcollection is specified, the following parameters are expected:
    !
    ! rcollection%DquickAccess(1) = multiplier in front of the operator.
    ! rcollection%IquickAccess(1) = x-position in the matrix where to set up the operator.
    ! rcollection%IquickAccess(2) = y-position in the matrix where to set up the operator.
    ! rcollection%IquickAccess(3) = 0, if the convection is a constant vector field.
    !                                  in this case:
    !                                  1D: rcollection%DquickAccess(2)   = d/dx x-velocity
    !                                  2D: rcollection%DquickAccess(2:3) = (d/dx,d/dy) x-velocity
    !                                      rcollection%DquickAccess(4:5) = (d/dx,d/dy) y-velocity
    !                                  3D: rcollection%DquickAccess(2:4) = (d/dx,d/dy,d/dz) x-velocity
    !                                      rcollection%DquickAccess(5:7) = (d/dx,d/dy,d/dz) y-velocity
    !                                      rcollection%DquickAccess(8:10) = (d/dx,d/dy,d/dz) z-velocity
    !                             = 1, if The convection is specified by a
    !                                  finite element velocity field. In this case,
    !                                  a finite element velocity field must be specified
    !                                  as parameter revalVectors to the call of 
    !                                  bma_buildMatrix.
!</description>

!<remarks>
    ! Remark 1: 
    ! Using the routines from feevaluation2, it is possible to specify
    ! a nonconstant velocity field. The corresponding code looks as follows:
    !
    ! <verb>
    !     use feevaluation2
    !
    !     ...
    !
    !     type(t_collection) :: rcollection
    !     type(t_scalarCubatureInfo) :: rcubatureInfo   ! Cubature formula
    !     type(t_matrixBlock) :: rmatrix                ! Matrix to be calculated
    !
    !     type(t_vectorBlock) :: rvelocity        ! The velocity field
    !     type(t_fev2Vectors) :: revalVectors     ! Collection of vectors to evaluate
    !
    !     ! Prepare the cubature formula
    !     call spdiscr_createDefCubStructure (..., rcubatureInfo, CUB_GEN_AUTO)
    !
    !     ...
    !
    !     rcollection%IquickAccess(1) = 1          ! x-Position in the matrix
    !     rcollection%IquickAccess(2) = 1          ! y-Position in the matrix
    !     rcollection%IquickAccess(3) = 1          ! Nonconstant viscosity
    !     rcollection%DquickAccess(1) = 1.0_DP     ! Scaling
    !
    !     ! Add the X-, Y- and Z-velocity to revalVectors
    !     call fev2_addVectorFieldToEvalList(revalVectors,1,&
    !         rvelocity%RvectorBlock(1),rvelocity%RvectorBlock(2),rvelocity%RvectorBlock(3))
    !
    !     ! Set up the matrix
    !     call bma_buildMatrix (rmatrix,BMA_CALC_STANDARD, bma_fcalc_convection_graduvw, &
    !         rcollection, revalVectors=revalVectors,rcubatureInfo=rcubatureInfo)
    !
    !     ! Release the vector structure
    !     call fev2_releaseVectorList(revalVectors)
    !
    !     ...
    !
    !     ! Release the cubature formula
    !     call spdiscr_releaseCubStructure (rcubatureInfo)
    ! </verb>
    !
    ! Remark 2:
    ! The routine currently assumes that all velocity components are discretised
    ! with the same FEM space.
    !
    ! Remark 3:
    ! The routine currently assumes that all velocity matrices are independent.
    ! Matrices sharing data are not supported. This cound be realised by
    ! taking care of the flags RmatrixData(:,:)%bsharedMatrixData which indicate
    ! which matrix data is shared.
    !
    ! Remark 4:
    ! Interleaved matrices are currently not supported.
!</remarks>

!<inputoutput>
    ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
    ! have to be filled with data.
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!</subroutine>

    ! Local variables
    integer :: ndim,ix,iy
    real(DP) :: dscale
    real(DP) :: dvelX1, dvelX2, dvelX3
    real(DP) :: dvelY1, dvelY2, dvelY3
    real(DP) :: dvelZ1, dvelZ2, dvelZ3
    logical :: bvelConst

    ! Dimension of the underlying space
    ndim = rmatrixAssembly%p_rtriangulation%ndim
    
    ! Get parameters
    dscale = 1.0_DP
    dvelX1 = 0.0_DP
    dvelX2 = 0.0_DP
    dvelX3 = 0.0_DP
    
    dvelY1 = 0.0_DP
    dvelY2 = 0.0_DP
    dvelY3 = 0.0_DP
    
    dvelZ1 = 0.0_DP
    dvelZ2 = 0.0_DP
    dvelZ3 = 0.0_DP
    
    bvelConst = .true.
    ix = 1
    iy = 1
    
    if (present(rcollection)) then
      dscale = rcollection%DquickAccess(1)
      
      ! Constant velocity?
      bvelConst = (rcollection%IquickAccess(3) .eq. 0)
      
      if (bvelConst) then
        ! Get the constant velocity
        dvelX1 = rcollection%DquickAccess(2)
        if (ndim .ge. NDIM2D) then
          dvelX2 = rcollection%DquickAccess(3)
          
          dvelY1 = rcollection%DquickAccess(4)
          dvelY2 = rcollection%DquickAccess(5)
        end if
        if (ndim .ge. NDIM3D) then
          dvelX2 = rcollection%DquickAccess(3)
          dvelX3 = rcollection%DquickAccess(4)
          
          dvelY1 = rcollection%DquickAccess(5)
          dvelY2 = rcollection%DquickAccess(6)
          dvelY3 = rcollection%DquickAccess(7)
          
          dvelZ1 = rcollection%DquickAccess(8)
          dvelZ2 = rcollection%DquickAccess(9)
          dvelZ3 = rcollection%DquickAccess(10)
        end if
      end if
      
      ! Position
      ix = rcollection%IquickAccess(1)
      iy = rcollection%IquickAccess(2)

    end if
    
    if (bvelConst) then
      call bma_docalc_convection_vgraduw(RmatrixData,rassemblyData,rmatrixAssembly,&
          npointsPerElement,nelements,dscale,ix,iy,&
          dvelX1,dvelX2,dvelX3,dvelY1,dvelY2,dvelY3,dvelZ1,dvelZ2,dvelZ3)
    else
      call bma_docalc_convection_vgraduw(RmatrixData,rassemblyData,rmatrixAssembly,&
          npointsPerElement,nelements,dscale,ix,iy,rvectorField=revalVectors%p_RvectorData(1))
    end if

  end subroutine

  !****************************************************************************
  !****************************************************************************
  !****************************************************************************

!<subroutine>

  subroutine bma_docalc_rhsConst(rvectorData,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,icomp,dval)

!<description>  
    ! Calculates a right-hand side vector according to the right-hand
    ! side function f=1 into a specific component.
!</description>

!<inputoutput>
    ! Vector data of all subvectors. The arrays p_Dentry of all subvectors
    ! have to be filled with data.
    type(t_bmaVectorData), dimension(:), intent(inout), target :: RvectorData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaVectorAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaVectorAssembly), intent(in) :: rvectorAssembly
    
    ! Number of points per element
    integer, intent(in) :: npointsPerElement
    
    ! Number of elements
    integer, intent(in) :: nelements
    
    ! ID of the component where to assemble to.
    integer, intent(in) :: icomp
    
    ! Value of the function.
    real(DP), intent(in) :: dval
!</input>
    
!</subroutine>

    ! Local variables
    real(DP) :: dbasI
    integer :: iel, icubp, idofe, idimfe, ndimfe
    real(DP), dimension(:,:), pointer :: p_DlocalVector
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaVectorData), pointer :: p_rvectorData
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    
    ! Get the data arrays of the subvector
    p_rvectorData => RvectorData(icomp)
    p_DlocalVector => p_rvectorData%p_Dentry
    p_DbasTest => p_rvectorData%p_DbasTest
  
    ! FE space dimension
    ndimfe = p_rvectorData%ndimfe

    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement

        ! Loop over the dimensions of the FE space
        do idimfe = 0,ndimfe-1
        
          ! Outer loop over the DOF's i=1..ndof on our current element,
          ! which corresponds to the (test) basis functions Psi_i:
          do idofe=1,p_rvectorData%ndofTest
          
            ! Fetch the contributions of the (test) basis functions Psi_i
            ! into dbasI
            dbasI = p_DbasTest(idofe+idimfe*p_rvectorData%ndofTest,DER_FUNC,icubp,iel)
            
            ! Multiply the values of the basis functions
            ! (1st derivatives) by the cubature weight and sum up
            ! into the local vectors.
            p_DlocalVector(idofe,iel) = p_DlocalVector(idofe,iel) + &
                p_DcubWeight(icubp,iel) * dval * dbasI
            
          end do ! jdofe
          
        end do ! idimfe

      end do ! icubp
    
    end do ! iel
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_rhsConst(rvectorData,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates a right-hand side vector according to the right-hand
    ! side function f=constant.
    !
    ! If rcollection is not specified, the rhs is calculated
    ! in all components using the constant value f=constant=1.
    !
    ! If rcollection is specified, the following parameters are expected:
    ! rcollection%IquickAccess(1) = Number of the component that should
    !                               receive the RHS.
    ! rcollection%DquickAccess(1) = Value of the function.
!</description>

!<inputoutput>
    ! Vector data of all subvectors. The arrays p_Dentry of all subvectors
    ! have to be filled with data.
    type(t_bmaVectorData), dimension(:), intent(inout), target :: RvectorData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaVectorAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaVectorAssembly), intent(in) :: rvectorAssembly
    
    ! Number of points per element
    integer, intent(in) :: npointsPerElement
    
    ! Number of elements
    integer, intent(in) :: nelements
    
    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>
    
!</subroutine>

    ! Local variables
    integer :: icomp,istartcomp,iendcomp
    real(DP) :: dval
  
    ! Either compute to all components or to a specified one.
    if (present(rcollection)) then
      istartcomp = rcollection%IquickAccess(1)
      iendcomp = rcollection%IquickAccess(1)
      dval = rcollection%DquickAccess(1)
    else
      istartcomp = 1
      iendcomp = size(rvectorData)
      dval = 1.0_DP
    end if
    
    ! Loop over the components
    do icomp = istartcomp,iendcomp

      call bma_docalc_rhsConst(rvectorData,rassemblyData,rvectorAssembly,&
          npointsPerElement,nelements,icomp,dval)
      
    end do ! icomp
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_docalc_rhsBubble(rvectorData,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,icomp,dscale)

!<description>  
    ! Calculates a right-hand side vector according to the right-hand
    ! side function f=dscale*32*y*(1-y)+32*x*(1-x).
!</description>

!<inputoutput>
    ! Vector data of all subvectors. The arrays p_Dentry of all subvectors
    ! have to be filled with data.
    type(t_bmaVectorData), dimension(:), intent(inout), target :: RvectorData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaVectorAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaVectorAssembly), intent(in) :: rvectorAssembly
    
    ! Number of points per element
    integer, intent(in) :: npointsPerElement
    
    ! Number of elements
    integer, intent(in) :: nelements
    
    ! ID of the component where to assemble to.
    integer, intent(in) :: icomp
    
    ! Scaling factor
    real(DP), intent(in) :: dscale
!</input>
    
!</subroutine>

    ! Local variables
    real(DP) :: dbasI, dval, dx, dy
    integer :: iel, icubp, idofe, idimfe, ndimfe
    real(DP), dimension(:,:), pointer :: p_DlocalVector
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaVectorData), pointer :: p_rvectorData
    real(DP), dimension(:,:,:), pointer :: p_Dpoints
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal
    
    ! Get the data arrays of the subvector
    p_rvectorData => RvectorData(icomp)
    p_DlocalVector => p_rvectorData%p_Dentry
    p_DbasTest => p_rvectorData%p_DbasTest
  
    ! FE space dimension
    ndimfe = p_rvectorData%ndimfe

    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement
      
        ! Get the coordinates of the cubature point.
        dx = p_Dpoints(1,icubp,iel)
        dy = p_Dpoints(2,icubp,iel)

        ! Calculate the values of the RHS using the coordinates
        ! of the cubature points.
        dval = dscale*32.0_DP*dy*(1.0_DP-dy) + 32_DP*dx*(1.0_DP-dx)
        
        ! Loop over the dimensions of the FE space
        do idimfe = 0,ndimfe-1

          ! Outer loop over the DOF's i=1..ndof on our current element,
          ! which corresponds to the (test) basis functions Psi_i:
          do idofe=1,p_rvectorData%ndofTest
          
            ! Fetch the contributions of the (test) basis functions Psi_i
            ! into dbasI
            dbasI = p_DbasTest(idofe+idimfe*p_rvectorData%ndofTest,DER_FUNC,icubp,iel)
            
            ! Multiply the values of the basis functions
            ! (1st derivatives) by the cubature weight and sum up
            ! into the local vectors.
            p_DlocalVector(idofe,iel) = p_DlocalVector(idofe,iel) + &
                p_DcubWeight(icubp,iel) * dval * dbasI
            
          end do ! idofe
          
        end do ! idimfe

      end do ! icubp
    
    end do ! iel
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_rhsBubble(rvectorData,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates a right-hand side vector according to the right-hand
    ! side function f=dscale*32*y*(1-y)+32*x*(1-x).
    !
    ! If rcollection is not specified, the rhs is calculated
    ! in all components.
    ! If rcollection is specified, the following parameters are expected:
    ! rcollection%IquickAccess(1) = Number of the component that should
    !                               receive the RHS.
    ! rcollection%DquickAccess(1) = Scaling factor dscale
!</description>

!<inputoutput>
    ! Vector data of all subvectors. The arrays p_Dentry of all subvectors
    ! have to be filled with data.
    type(t_bmaVectorData), dimension(:), intent(inout), target :: RvectorData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaVectorAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaVectorAssembly), intent(in) :: rvectorAssembly
    
    ! Number of points per element
    integer, intent(in) :: npointsPerElement
    
    ! Number of elements
    integer, intent(in) :: nelements
    
    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>
    
!</subroutine>

    ! Local variables
    real(DP) :: dscale
    integer :: icomp,istartcomp,iendcomp
  
    ! Either compute to all components or to a specified one.
    if (present(rcollection)) then
      istartcomp = rcollection%IquickAccess(1)
      iendcomp = rcollection%IquickAccess(1)
      dscale = rcollection%DquickAccess(1)
    else
      istartcomp = 1
      iendcomp = size(rvectorData)
      dscale = 1.0_DP
    end if
    
    ! Loop over the components
    do icomp = istartcomp,iendcomp
      
      call bma_docalc_rhsBubble(rvectorData,rassemblyData,rvectorAssembly,&
          npointsPerElement,nelements,icomp,dscale)
      
    end do ! icomp
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_docalc_rhsBubblePlusFE(rvectorData,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,icomp,dscale,rcoefficient)

!<description>  
    ! Calculates a right-hand side vector according to the right-hand
    ! side function f = dscale * 32*y*(1-y)+32*x*(1-x) + v(x,y)
    ! with v being a FEM function.
    !
    ! The FEM function v(x,y) must be provided in rcoefficient.
    ! The routine only supports non-interleaved vectors.
!</description>

!<inputoutput>
    ! Vector data of all subvectors. The arrays p_Dentry of all subvectors
    ! have to be filled with data.
    type(t_bmaVectorData), dimension(:), intent(inout), target :: RvectorData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaVectorAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaVectorAssembly), intent(in) :: rvectorAssembly
    
    ! Number of points per element
    integer, intent(in) :: npointsPerElement
    
    ! Number of elements
    integer, intent(in) :: nelements
    
    ! ID of the component where to assemble to.
    integer, intent(in) :: icomp
    
    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Evaluation structure that describes an underlying nonconstant coefficient.
    type(t_fev2VectorData), intent(in), optional :: rcoefficient
!</input>
    
!</subroutine>

    ! Local variables
    real(DP) :: dbasI, dval, dx, dy
    integer :: iel, icubp, idofe, idimfe, ndimfe
    real(DP), dimension(:,:), pointer :: p_DlocalVector
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaVectorData), pointer :: p_rvectorData
    real(DP), dimension(:,:,:), pointer :: p_Dpoints
    real(DP), dimension(:,:), pointer :: p_Dfunc
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal
    
    ! Get the data array with the values of the FEM function
    ! in the cubature points
    p_Dfunc => rcoefficient%p_Ddata(:,:,DER_FUNC2D)
    
    ! Get the data arrays of the subvector
    p_rvectorData => RvectorData(icomp)
    p_DlocalVector => p_rvectorData%p_Dentry
    p_DbasTest => p_rvectorData%p_DbasTest

    ! FE space dimension
    ndimfe = p_rvectorData%ndimfe
  
    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement
      
        ! Get the coordinates of the cubature point.
        dx = p_Dpoints(1,icubp,iel)
        dy = p_Dpoints(2,icubp,iel)

        ! Calculate the values of the RHS in the cubature point:
        !     f = 32*y*(1-y)+32*x*(1-x) + v(x,y)
        dval = dscale * (32.0_DP*dy*(1.0_DP-dy) + 32_DP*dx*(1.0_DP-dx)) + p_Dfunc(icubp,iel)
        
        ! Loop over the dimensions of the FE space
        do idimfe = 0,ndimfe-1

          ! Outer loop over the DOF's i=1..ndof on our current element,
          ! which corresponds to the (test) basis functions Psi_i:
          do idofe=1,p_rvectorData%ndofTest
          
            ! Fetch the contributions of the (test) basis functions Psi_i
            ! into dbasI
            dbasI = p_DbasTest(idofe+idimfe*p_rvectorData%ndofTest,DER_FUNC,icubp,iel)
            
            ! Multiply the values of the basis functions
            ! (1st derivatives) by the cubature weight and sum up
            ! into the local vectors.
            p_DlocalVector(idofe,iel) = p_DlocalVector(idofe,iel) + &
                p_DcubWeight(icubp,iel) * dval * dbasI
            
          end do ! idofe
          
        end do ! idimfe

      end do ! icubp
    
    end do ! iel
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_rhsBubblePlusFE(rvectorData,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates a right-hand side vector according to the right-hand
    ! side function f = dscale * (32*y*(1-y)+32*x*(1-x)) + v(x,y)
    ! with v being a FEM function.
    !
    ! The FEM function v(x,y) must be provided in revalVectors.
    ! The routine only supports non-interleaved vectors.
    !
    ! If rcollection is not specified, the rhs is calculated
    ! in all components.
    ! If rcollection is specified, the following parameters are expected:
    ! rcollection%IquickAccess(1) = Number of the component that should
    !                               receive the RHS.
    ! rcollection%DquickAccess(1) = Scaling factor dscale
!</description>

!<inputoutput>
    ! Vector data of all subvectors. The arrays p_Dentry of all subvectors
    ! have to be filled with data.
    type(t_bmaVectorData), dimension(:), intent(inout), target :: RvectorData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaVectorAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaVectorAssembly), intent(in) :: rvectorAssembly
    
    ! Number of points per element
    integer, intent(in) :: npointsPerElement
    
    ! Number of elements
    integer, intent(in) :: nelements
    
    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>
    
!</subroutine>

    ! Local variables
    real(DP) :: dscale
    integer :: icomp,istartcomp,iendcomp

    if (revalVectors%ncount .eq. 0) then
      call output_line ("FEM function missing.",&
          OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_rhsBubblePlusFE")
      call sys_halt()
    end if
    
    ! Either compute to all components or to a specified one.
    if (present(rcollection)) then
      istartcomp = rcollection%IquickAccess(1)
      iendcomp = rcollection%IquickAccess(1)
      dscale = rcollection%DquickAccess(1)
    else
      istartcomp = 1
      iendcomp = size(rvectorData)
      dscale = 1.0_DP
    end if
    
    ! Loop over the components
    do icomp = istartcomp,iendcomp

      call bma_docalc_rhsBubblePlusFE(rvectorData,rassemblyData,rvectorAssembly,&
          npointsPerElement,nelements,icomp,dscale,&
          revalVectors%p_RvectorData(icomp-istartcomp+1))

    end do ! icomp
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_docalc_rhsFE(rvectorData,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,icomp,dscale,rcoefficient)

!<description>  
    ! Calculates a right-hand side vector according to the right-hand
    ! side function f = v(x,y)
    ! with v being a FEM function.
    !
    ! The FEM function v(x,y) must be provided in revalVectors.
    ! It must have as many components as the right hand side f,
    ! f_i is computed from v_i.
    !
    ! The routine only supports non-interleaved vectors.
    !
    ! This routine is typically used in L2 projections for setting up the
    ! RHS based on a FEM function.
!</description>

!<inputoutput>
    ! Vector data of all subvectors. The arrays p_Dentry of all subvectors
    ! have to be filled with data.
    type(t_bmaVectorData), dimension(:), intent(inout), target :: RvectorData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaVectorAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaVectorAssembly), intent(in) :: rvectorAssembly
    
    ! Number of points per element
    integer, intent(in) :: npointsPerElement
    
    ! Number of elements
    integer, intent(in) :: nelements
    
    ! ID of the component where to assemble to.
    integer, intent(in) :: icomp
    
    ! Scaling factor
    real(DP), intent(in) :: dscale

    ! Evaluation structure that describes an underlying nonconstant coefficient.
    type(t_fev2VectorData), intent(in), optional :: rcoefficient
!</input>
    
!</subroutine>

    ! Local variables
    real(DP) :: dbasI, dval
    integer :: ndimfe, idimfe
    integer :: iel, icubp, idofe
    real(DP), dimension(:,:), pointer :: p_DlocalVector
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaVectorData), pointer :: p_rvectorData
    real(DP), dimension(:,:), pointer :: p_Dfunc
    real(DP), dimension(:,:,:), pointer :: p_DfuncVec
  
    ! Get the data arrays of the subvector
    p_rvectorData => RvectorData(icomp)
    p_DlocalVector => p_rvectorData%p_Dentry
    p_DbasTest => p_rvectorData%p_DbasTest
    
    ! FE space dimension
    ndimfe = p_rvectorData%ndimfe
  
    if (ndimfe .eq. 1) then

      ! -----------------------
      ! Scalar-valued FE space.
      ! -----------------------

      ! Get the data array with the values of the FEM function v_i
      ! in the cubature points.
      p_Dfunc => rcoefficient%p_Ddata(:,:,DER_FUNC2D)
      
      ! Loop over the elements in the current set.
      do iel = 1,nelements

        ! Loop over all cubature points on the current element
        do icubp = 1,npointsPerElement
        
          ! Calculate the values of the RHS in the cubature point:
          !     f_i = dscale * v_i(x,y)
          dval = p_Dfunc(icubp,iel)
          
          ! Outer loop over the DOF's i=1..ndof on our current element,
          ! which corresponds to the (test) basis functions Psi_i:
          do idofe=1,p_rvectorData%ndofTest
          
            ! Fetch the contributions of the (test) basis functions Psi_i
            ! into dbasI
            dbasI = dscale * p_DbasTest(idofe,DER_FUNC,icubp,iel)
            
            ! Multiply the values of the basis functions
            ! (1st derivatives) by the cubature weight and sum up
            ! into the local vectors.
            p_DlocalVector(idofe,iel) = p_DlocalVector(idofe,iel) + &
                p_DcubWeight(icubp,iel) * dval * dbasI
            
          end do ! jdofe

        end do ! icubp
      
      end do ! iel
      
    else

      ! -----------------------
      ! Vector-valued FE space.
      ! -----------------------
    
      ! Get the data array with the values of the FEM function v_i
      ! in the cubature points.
      p_DfuncVec => rcoefficient%p_DdataVec(:,:,:,DER_FUNC2D)
      
      ! Loop over the elements in the current set.
      do iel = 1,nelements

        ! Loop over all cubature points on the current element
        do icubp = 1,npointsPerElement
        
          ! Loop over the dimensions of the FE space
          do idimfe = 0,ndimfe-1

            ! Calculate the values of the RHS in the cubature point
            ! for component (1+idimfe):
            !     f_i = dscale * v_i(x,y)
            dval = dscale * p_DfuncVec(1+idimfe,icubp,iel)
            
            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Psi_i:
            do idofe=1,p_rvectorData%ndofTest
            
              ! Fetch the contributions of the (test) basis functions Psi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe+idimfe*p_rvectorData%ndofTest,DER_FUNC,icubp,iel)
              
              ! Multiply the values of the basis functions
              ! (1st derivatives) by the cubature weight and sum up
              ! into the local vectors.
              p_DlocalVector(idofe,iel) = p_DlocalVector(idofe,iel) + &
                  p_DcubWeight(icubp,iel) * dval * dbasI
              
            end do ! idofe
          
          end do ! idimfe

        end do ! icubp
      
      end do ! iel
    
    end if
      
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_rhsFE(rvectorData,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates a right-hand side vector according to the right-hand
    ! side function f = dscale * v(x,y)
    ! with v being a FEM function with as many components as f.
    !
    ! The FEM function v(x,y) must be provided in revalVectors.
    ! It must have as many components as the right hand side f,
    ! f_i is computed from v_i.
    !
    ! The routine only supports non-interleaved vectors.
    !
    ! This routine is typically used in L2 projections for setting up the
    ! RHS based on a FEM function.
    !
    ! If rcollection is not specified, the rhs is calculated
    ! in all components.
    ! If rcollection is specified, the following parameters are expected:
    ! rcollection%IquickAccess(1) = Number of the component that should
    !                               receive the RHS.
    ! rcollection%DquickAccess(1) = Scaling factor dscale
!</description>

!<inputoutput>
    ! Vector data of all subvectors. The arrays p_Dentry of all subvectors
    ! have to be filled with data.
    type(t_bmaVectorData), dimension(:), intent(inout), target :: RvectorData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaVectorAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaVectorAssembly), intent(in) :: rvectorAssembly
    
    ! Number of points per element
    integer, intent(in) :: npointsPerElement
    
    ! Number of elements
    integer, intent(in) :: nelements
    
    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>
    
!</subroutine>

    ! Local variables
    real(DP) :: dscale
    integer :: icomp,istartcomp,iendcomp
    
    if (revalVectors%ncount .ne. size(rvectorData)) then
      call output_line ("FEM function missing or wrong length.",&
          OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_rhsFE")
      call sys_halt()
    end if
    
    ! Either compute to all components or to a specified one.
    if (present(rcollection)) then
      istartcomp = rcollection%IquickAccess(1)
      iendcomp = rcollection%IquickAccess(1)
      dscale = rcollection%DquickAccess(1)
    else
      istartcomp = 1
      iendcomp = size(rvectorData)
      dscale = 1.0_DP
    end if
    
    ! Loop over the components
    do icomp = istartcomp,iendcomp

      call bma_docalc_rhsFE(rvectorData,rassemblyData,rvectorAssembly,&
          npointsPerElement,nelements,icomp,dscale,&
          revalVectors%p_RvectorData(icomp-istartcomp+1))

    end do ! icomp
    
  end subroutine

  !****************************************************************************
  !****************************************************************************
  !****************************************************************************

!<subroutine>

  subroutine bma_docalc_integralOne(dintvalue,rassemblyData,rintAssembly,&
      npointsPerElement,nelements)

!<description>  
    ! Calculates the integral of the function v=1.
    ! The result is the size of the domain.
!</description>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements
!</input>

!<output>
    ! Returns the value of the integral
    real(DP), intent(out) :: dintvalue
!</output>    

!</subroutine>

    ! Local variables
    real(DP) :: dval
    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    
    dintvalue = 0.0_DP
    
    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement

        ! The value of the function v is just the 1.0-function here
        dval = 1.0_DP
        
        ! Multiply the values by the cubature weight and sum up
        ! into the integral value
        dintvalue = dintvalue + p_DcubWeight(icubp,iel) * dval
          
      end do ! icubp
    
    end do ! iel
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_integralOne(dintvalue,rassemblyData,rintAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the integral of the function v=1.
    ! The result is the size of the domain.
!</description>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!<output>
    ! Returns the value of the integral
    real(DP), intent(out) :: dintvalue
!</output>    

!</subroutine>

    call bma_docalc_integralOne(dintvalue,rassemblyData,rintAssembly,&
        npointsPerElement,nelements)
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_docalc_integralFE(dintvalue,rassemblyData,rintAssembly,&
      npointsPerElement,nelements,rfunction,icomp)

!<description>  
    ! Calculates the integral of an arbitrary finite element function v.
    !
    ! If v is vector valued components, the sum of the integrals of all
    ! components is returned. If v is vector valued and icomp is specified,
    ! only the integral of component icomp is returned.
    !
    ! The FEM function(s) must be provided in rfunction.
    ! The routine only supports non-interleaved vectors.
!</description>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Function to integrate.
    type(t_fev2VectorData), intent(in) :: rfunction
    
    ! OPTIONAL: Number of the component.
    ! Only applies for vector valued FE functions.
    integer, intent(in), optional :: icomp
!</input>

!<output>
    ! Returns the value of the integral
    real(DP), intent(out) :: dintvalue
!</output>    

!</subroutine>

    ! Local variables
    real(DP) :: dval
    integer :: iel, icubp
    integer :: icomponent,istart,iend
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:), pointer :: p_Dfunc
    real(DP), dimension(:,:,:), pointer :: p_DfuncVec
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    
    dintvalue = 0.0_DP

    ! Skip interleaved vectors.
    if (rfunction%bisInterleaved) return
    
    ! Skip vector-valued FE functions
    if (.not. rfunction%bisVectorField) then

      ! Get the data array with the values of the FEM function
      ! in the cubature points
      p_Dfunc => rfunction%p_Ddata(:,:,DER_FUNC2D)

      ! Loop over the elements in the current set.
      do iel = 1,nelements

        ! Loop over all cubature points on the current element
        do icubp = 1,npointsPerElement

          ! Get the value of the FEM function
          dval = p_Dfunc(icubp,iel)
          
          ! Multiply the values by the cubature weight and sum up
          ! into the integral value
          dintvalue = dintvalue + p_DcubWeight(icubp,iel) * dval
            
        end do ! icubp
      
      end do ! iel
    
    else
    
      ! Get the data array with the values of the FEM function
      ! in the cubature points
      p_DfuncVec => rfunction%p_DdataVec(:,:,:,DER_FUNC2D)

      ! Loop through the components
      istart = 1
      iend = ubound(p_DfuncVec,1)
      if (present(icomp)) then
        istart = icomp
        iend = icomp
      end if
      
      ! Loop over the elements in the current set.
      do iel = 1,nelements

        ! Loop over all cubature points on the current element
        do icubp = 1,npointsPerElement

          do icomponent = istart,iend

            ! Get the value of the FEM function
            dval = p_DfuncVec(icomponent,icubp,iel)
            
            ! Multiply the values by the cubature weight and sum up
            ! into the integral value
            dintvalue = dintvalue + p_DcubWeight(icubp,iel) * dval
            
          end do
            
        end do ! icubp
      
      end do ! iel
    
    end if
      
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_integralFE(dintvalue,rassemblyData,rintAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the integral of an arbitrary finite element function v.
    ! If v has is a vector field, the sum of the integrals of all
    ! components is returned.
    !
    ! The FEM function(s) must be provided in revalVectors.
    ! The routine only supports non-interleaved vectors.
!</description>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!<output>
    ! Returns the value of the integral
    real(DP), intent(out) :: dintvalue
!</output>    

!</subroutine>

    ! Local variables
    integer :: ivec
    real(DP) :: dint
  
    dintvalue = 0.0_DP

    ! Loop through all provided FEM functions
    do ivec = 1,revalVectors%ncount
    
      ! Calculate the integral, sum up
      call bma_docalc_integralFE(dint,rassemblyData,rintAssembly,&
          npointsPerElement,nelements,revalVectors%p_RvectorData(ivec))
          
      dintvalue = dintvalue + dint
      
    end do ! ivec
    
  end subroutine

  !****************************************************************************

!<--
! The L2-norm of a vector field "rvector" can be calculated, e.g., as follows:
!
!    ! Declare an evaluation structure for the subvectors
!    type(t_fev2Vectors) :: revalVectors
!    ...
!    real(DP) :: dresult
!
!    ! Add the subvectors
!    call fev2_addVectorToEvalList(revalVectors,rvector%RvectorBlock(1),0)
!    call fev2_addVectorToEvalList(revalVectors,rvector%RvectorBlock(2),0)
!    ...
!
!    ! Evaluate the norm to "dresult"
!    call bma_buildIntegral(dresult,BMA_CALC_STANDARD,bma_fcalc_L2norm,&
!        revalVectors=revalVectors,rcubatureInfo=rcubatureInfo)
!
!    ! Clean up
!    call fev2_releaseVectorList(revalVectors)
!-->

  !****************************************************************************

!<subroutine>

  subroutine bma_docalc_L2norm(dintvalue,rassemblyData,rintAssembly,&
      npointsPerElement,nelements,rfunction)

!<description>  
    ! Calculates the squared L2-norm of an arbitrary finite element 
    ! function: <tex> $ |||v||^2 $ </tex>.
    ! If v is a vector field, the sum of the integrals of all
    ! components is returned.
    !
    ! The FEM function(s) must be provided in rfunction.
    ! The routine only supports non-interleaved vectors.
!</description>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Function to compute the L2 norm from.
    type(t_fev2VectorData), intent(in) :: rfunction
!</input>

!<output>
    ! Returns the value of the integral
    real(DP), intent(out) :: dintvalue
!</output>    

!</subroutine>

    ! Local variables
    real(DP) :: dval
    integer :: iel, icubp, idimfe
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:), pointer :: p_Dfunc
    real(DP), dimension(:,:,:), pointer :: p_DfuncVec
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    
    dintvalue = 0.0_DP

    ! Skip interleaved vectors.
    if (rfunction%bisInterleaved) return

    if (.not. rfunction%bisVectorField) then
    
      ! -----------------------
      ! Scalar-valued FE space. 
      ! -----------------------
    
      ! Get the data array with the values of the FEM function
      ! in the cubature points
      p_Dfunc => rfunction%p_Ddata(:,:,DER_FUNC)

      ! Loop over the elements in the current set.
      do iel = 1,nelements

        ! Loop over all cubature points on the current element
        do icubp = 1,npointsPerElement

          ! Get the value of the FEM function
          dval = p_Dfunc(icubp,iel)
          
          ! Multiply the values by the cubature weight and sum up
          ! into the integral value
          dintvalue = dintvalue + p_DcubWeight(icubp,iel) * dval**2
            
        end do ! icubp
      
      end do ! iel
      
    else
    
      ! -----------------------
      ! Vector-valued FE space.
      ! -----------------------
      ! Get the data array with the values of the FEM function
      ! in the cubature points
      p_DfuncVec => rfunction%p_DdataVec(:,:,:,DER_FUNC)

      ! Loop over the elements in the current set.
      do iel = 1,nelements

        ! Loop over all cubature points on the current element
        do icubp = 1,npointsPerElement

          ! Loop over the dimensions of the FE space
          do idimfe = 1,rfunction%ndimVectorField

            ! Get the value of the FEM function
            dval = p_DfuncVec(idimfe,icubp,iel)
            
            ! Multiply the values by the cubature weight and sum up
            ! into the integral value
            dintvalue = dintvalue + p_DcubWeight(icubp,iel) * dval**2
            
          end do
            
        end do ! icubp
      
      end do ! iel
      
    end if
      
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_L2norm(dintvalue,rassemblyData,rintAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the squared L2-norm of an arbitrary finite element 
    ! function: <tex> $ |||v||^2 $ </tex>.
    ! If v has multiple components, the sum of the integrals of all
    ! components is returned.
    !
    ! The FEM function(s) must be provided in revalVectors.
    ! The routine only supports non-interleaved vectors.
!</description>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!<output>
    ! Returns the value of the integral
    real(DP), intent(out) :: dintvalue
!</output>    

!</subroutine>

    ! Local variables
    integer :: ivec
    real(DP) :: dint
  
    dintvalue = 0.0_DP

    ! Loop through all provided FEM functions
    do ivec = 1,revalVectors%ncount
    
      ! Calculate the integral, sum up
      call bma_docalc_L2norm(dintvalue,rassemblyData,rintAssembly,&
          npointsPerElement,nelements,revalVectors%p_RvectorData(ivec))
          
      dintvalue = dintvalue + dint
      
    end do ! ivec

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_docalc_H1norm(dintvalue,rassemblyData,rintAssembly,&
      npointsPerElement,nelements,rfunction)

!<description>  
    ! Calculates the squared H1 (semi-)norm of an arbitrary finite element 
    ! function <tex> $ ||v||^2_{H^1} $ </tex>.
    ! If v is vector valued, the sum of the integrals of all
    ! components is returned.
    !
    ! The FEM function(s) must be provided in rfunction.
    ! The routine only supports non-interleaved vectors.
!</description>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Function to compute the H1 norm from.
    type(t_fev2VectorData), intent(in) :: rfunction
!</input>

!<output>
    ! Returns the value of the integral
    real(DP), intent(out) :: dintvalue
!</output>    

!</subroutine>

    ! Local variables
    real(DP) :: dderivX,dderivY,dderivZ
    integer :: iel, icubp, ndim, idimfe
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:), pointer :: p_DderivX,p_DderivY,p_DderivZ
    real(DP), dimension(:,:,:), pointer :: p_DderivXVec,p_DderivYVec,p_DderivZVec
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    
    dintvalue = 0.0_DP

    ! Skip interleaved vectors.
    if (rfunction%bisInterleaved) return

    ! Dimension?
    ndim = rintAssembly%p_rtriangulation%ndim

    if (.not. rfunction%bisVectorField) then
    
      ! -----------------------
      ! Scalar-valued FE space. 
      ! -----------------------

      select case (ndim)
      
      ! ----------------
      ! 1D
      ! ----------------
      case (NDIM1D)
      
        ! Get the data array with the values of the FEM function
        ! in the cubature points
        p_DderivX => rfunction%p_Ddata(:,:,DER_DERIV1D_X)

        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Get the value of the FEM function
            dderivX = p_DderivX(icubp,iel)
            
            ! Multiply the values by the cubature weight and sum up
            ! into the integral value
            dintvalue = dintvalue + p_DcubWeight(icubp,iel) * dderivX**2
              
          end do ! icubp
        
        end do ! iel

      ! ----------------
      ! 2D
      ! ----------------
      case (NDIM2D)
      
        ! Get the data array with the values of the FEM function
        ! in the cubature points
        p_DderivX => rfunction%p_Ddata(:,:,DER_DERIV2D_X)
        p_DderivY => rfunction%p_Ddata(:,:,DER_DERIV2D_Y)

        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Get the value of the FEM function
            dderivX = p_DderivX(icubp,iel)
            dderivY = p_DderivY(icubp,iel)
            
            ! Multiply the values by the cubature weight and sum up
            ! into the integral value
            dintvalue = dintvalue + p_DcubWeight(icubp,iel) * &
                ( dderivX**2 + dderivY**2 )
              
          end do ! icubp
        
        end do ! iel

      ! ----------------
      ! 3D
      ! ----------------
      case (NDIM3D)
      
        ! Get the data array with the values of the FEM function
        ! in the cubature points
        p_DderivX => rfunction%p_Ddata(:,:,DER_DERIV3D_X)
        p_DderivY => rfunction%p_Ddata(:,:,DER_DERIV3D_Y)
        p_DderivZ => rfunction%p_Ddata(:,:,DER_DERIV3D_Z)

        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Get the value of the FEM function
            dderivX = p_DderivX(icubp,iel)
            dderivY = p_DderivY(icubp,iel)
            dderivZ = p_DderivZ(icubp,iel)
            
            ! Multiply the values by the cubature weight and sum up
            ! into the integral value
            dintvalue = dintvalue + p_DcubWeight(icubp,iel) * &
                ( dderivX**2 + dderivY**2 + dderivZ**2 )
              
          end do ! icubp
        
        end do ! iel
        
      end select

    else
    
      ! -----------------------
      ! Vector-valued FE space. 
      ! -----------------------

      select case (ndim)
      
      ! ----------------
      ! 1D
      ! ----------------
      case (NDIM1D)
      
        ! Get the data array with the values of the FEM function
        ! in the cubature points
        p_DderivXVec => rfunction%p_DdataVec(:,:,:,DER_DERIV1D_X)

        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Loop over the dimensions of the FE space
            do idimfe = 1,rfunction%ndimVectorField
              
              ! Get the value of the FEM function
              dderivX = p_DderivXVec(idimfe,icubp,iel)
              
              ! Multiply the values by the cubature weight and sum up
              ! into the integral value
              dintvalue = dintvalue + p_DcubWeight(icubp,iel) * dderivX**2
              
            end do ! idimfe
              
          end do ! icubp
        
        end do ! iel

      ! ----------------
      ! 2D
      ! ----------------
      case (NDIM2D)
      
        ! Get the data array with the values of the FEM function
        ! in the cubature points
        p_DderivXVec => rfunction%p_DdataVec(:,:,:,DER_DERIV2D_X)
        p_DderivYVec => rfunction%p_DdataVec(:,:,:,DER_DERIV2D_Y)

        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Loop over the dimensions of the FE space
            do idimfe = 1,rfunction%ndimVectorField
              
              ! Get the value of the FEM function
              dderivX = p_DderivXVec(idimfe,icubp,iel)
              dderivY = p_DderivYVec(idimfe,icubp,iel)
              
              ! Multiply the values by the cubature weight and sum up
              ! into the integral value
              dintvalue = dintvalue + p_DcubWeight(icubp,iel) * &
                  ( dderivX**2 + dderivY**2 )
                  
            end do ! idimfe
              
          end do ! icubp
        
        end do ! iel

      ! ----------------
      ! 3D
      ! ----------------
      case (NDIM3D)
      
        ! Get the data array with the values of the FEM function
        ! in the cubature points
        p_DderivXVec => rfunction%p_DdataVec(:,:,:,DER_DERIV3D_X)
        p_DderivYVec => rfunction%p_DdataVec(:,:,:,DER_DERIV3D_Y)
        p_DderivZVec => rfunction%p_DdataVec(:,:,:,DER_DERIV3D_Z)

        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Loop over the dimensions of the FE space
            do idimfe = 1,rfunction%ndimVectorField
            
              ! Get the value of the FEM function
              dderivX = p_DderivXVec(idimfe,icubp,iel)
              dderivY = p_DderivYVec(idimfe,icubp,iel)
              dderivZ = p_DderivZVec(idimfe,icubp,iel)
              
              ! Multiply the values by the cubature weight and sum up
              ! into the integral value
              dintvalue = dintvalue + p_DcubWeight(icubp,iel) * &
                  ( dderivX**2 + dderivY**2 + dderivZ**2 )
            
            end do
              
          end do ! icubp
        
        end do ! iel
        
      end select

    end if
      
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_H1norm(dintvalue,rassemblyData,rintAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the squared H1 (semi-)norm of an arbitrary finite element 
    ! function <tex> $ ||v||^2_{H^1} $ </tex>.
    ! If v has multiple components, the sum of the integrals of all
    ! components is returned.
    !
    ! The FEM function(s) must be provided in revalVectors.
    ! The routine only supports non-interleaved vectors.
!</description>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!<output>
    ! Returns the value of the integral
    real(DP), intent(out) :: dintvalue
!</output>    

!</subroutine>

    ! Local variables
    integer :: ivec
    real(DP) :: dint
  
    dintvalue = 0.0_DP

    ! Loop through all provided FEM functions
    do ivec = 1,revalVectors%ncount
    
      ! Calculate the integral, sum up
      call bma_docalc_H1norm(dintvalue,rassemblyData,rintAssembly,&
          npointsPerElement,nelements,revalVectors%p_RvectorData(ivec))
          
      dintvalue = dintvalue + dint
      
    end do ! ivec

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_docalc_bubbleL2error(dintvalue,rassemblyData,rintegralAssembly,&
      npointsPerElement,nelements,rfunction)

!<description>  
    ! Calculates the (squared) L2 error of an arbitrary FEM function v
    ! to the bubble function u=16x(1-x)y(1-y):
    !
    !  <tex> $$ || v - u ||^2_{L2} $$ </tex>
    !
    ! The FEM function(s) must be provided in rfunction.
    ! The routine only supports scalar non-interleaved vectors.
!</description>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintegralAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Function to compute the error from.
    type(t_fev2VectorData), intent(in) :: rfunction
!</input>

!<output>
    ! Returns the value of the integral
    real(DP), intent(out) :: dintvalue
!</output>    

!</subroutine>

    ! Local variables
    real(DP) :: dval1,dval2,dx,dy
    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:), pointer :: p_Dfunc
    real(DP), dimension(:,:,:), pointer :: p_Dpoints
  
    ! Skip interleaved vectors.
    if (rfunction%bisInterleaved) return

    ! Skip vector-valued FE functions
    if (rfunction%bisVectorField) return

    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    
    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal
    
    ! Get the data array with the values of the FEM function
    ! in the cubature points
    p_Dfunc => rfunction%p_Ddata(:,:,DER_FUNC2D)

    dintvalue = 0.0_DP
    
    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement

        ! Get the value of the bubble function
        dx = p_Dpoints(1,icubp,iel)
        dy = p_Dpoints(2,icubp,iel)
        
        dval1 = 16.0_DP * dx * (1.0_DP-dx) * dy * (1.0_DP-dy)

        ! Get the error of the FEM function to the bubble function
        dval2 = p_Dfunc(icubp,iel)
        
        ! Multiply the values by the cubature weight and sum up
        ! into the (squared) L2 error:
        dintvalue = dintvalue + &
            p_DcubWeight(icubp,iel) * (dval1 - dval2)**2
          
      end do ! icubp
    
    end do ! iel
      
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_bubbleL2error(dintvalue,rassemblyData,rintegralAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the (squared) L2 error of an arbitrary FEM function v
    ! to the bubble function u=16x(1-x)y(1-y):
    !
    !  <tex> $$ || v - u ||^2_{L2} $$ </tex>
    !
    ! The FEM function must be provided in revalVectors.
    ! The routine only supports non-interleaved vectors.
!</description>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintegralAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!<output>
    ! Returns the value of the integral
    real(DP), intent(out) :: dintvalue
!</output>    

!</subroutine>

    call bma_docalc_bubbleL2error(dintvalue,rassemblyData,rintegralAssembly,&
        npointsPerElement,nelements,revalVectors%p_RvectorData(1))
      
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_docalc_bubbleH1error(dintvalue,rassemblyData,rintegralAssembly,&
      npointsPerElement,nelements,rfunction)

!<description>  
    ! Calculates the (squared) H1 error of an arbitrary FEM function v
    ! to the bubble function u=16x(1-x)y(1-y)
    ! (based on the H1 semi-norm).
    !
    !  <tex> $$ | v - u |^2_{H1} $$ </tex>
    !
    ! The FEM function must be provided in revalVectors.
    ! The routine only supports non-interleaved vectors.
!</description>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintegralAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Function to compute the error from.
    type(t_fev2VectorData), intent(in) :: rfunction
!</input>

!<output>
    ! Returns the value of the integral
    real(DP), intent(out) :: dintvalue
!</output>    

!</subroutine>

    ! Local variables
    real(DP) :: dderivX1,dderivY1,dderivX2,dderivY2,dx,dy
    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:), pointer :: p_DderivX,p_DderivY
    real(DP), dimension(:,:,:), pointer :: p_Dpoints
  
    ! Skip interleaved vectors.
    if (rfunction%bisInterleaved) return

    ! Skip vector-valued FE functions
    if (rfunction%bisVectorField) return

    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    
    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal
    
    ! Get the data array with the values of the FEM function
    ! in the cubature points
    p_DderivX => rfunction%p_Ddata(:,:,DER_DERIV2D_X)
    p_DderivY => rfunction%p_Ddata(:,:,DER_DERIV2D_Y)

    dintvalue = 0.0_DP
    
    ! Loop over the elements in the current set.
    do iel = 1,nelements

      ! Loop over all cubature points on the current element
      do icubp = 1,npointsPerElement

        ! Get the derivatives of the bubble function in the cubature point
        dx = p_Dpoints(1,icubp,iel)
        dy = p_Dpoints(2,icubp,iel)
        
        dderivX1 = 16.0_DP*dy*(dy-1.0_DP)*(2.0_DP*dx-1.0_DP)
        dderivY1 = 16.0_DP*dx*(dx-1.0_DP)*(2.0_DP*dy-1.0_DP)

        ! Get the error of the FEM function derivatives of the bubble function
        ! in the cubature point
        dderivX2 = p_DderivX(icubp,iel)
        dderivY2 = p_DderivY(icubp,iel)
        
        ! Multiply the values by the cubature weight and sum up
        ! into the (squared) H1 error:
        dintvalue = dintvalue + p_DcubWeight(icubp,iel) * &
            ( (dderivX1 - dderivX2)**2 + (dderivY1 - dderivY2)**2 )
          
      end do ! icubp
    
    end do ! iel
      
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_bubbleH1error(dintvalue,rassemblyData,rintegralAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the (squared) H1 error of an arbitrary FEM function v
    ! to the bubble function u=16x(1-x)y(1-y)
    ! (based on the H1 semi-norm).
    !
    !  <tex> $$ | v - u |^2_{H1} $$ </tex>
    !
    ! The FEM function must be provided in revalVectors.
    ! The routine only supports non-interleaved vectors.
!</description>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintegralAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!<output>
    ! Returns the value of the integral
    real(DP), intent(out) :: dintvalue
!</output>    

!</subroutine>

    call bma_docalc_bubbleH1error(dintvalue,rassemblyData,rintegralAssembly,&
        npointsPerElement,nelements,revalVectors%p_RvectorData(1))
      
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_docalc_divergenceL2norm(dintvalue,rassemblyData,rintAssembly,&
      npointsPerElement,nelements,rfunction)

!<description>  
    ! Calculates the (squared) L2-norm of the divergence of an arbitrary finite 
    ! element vector field <tex> $ ||div v||^2_{L_2} $ </tex>.
    !
    ! The vector field including first derivatives must be provided in 
    ! revalVectors. The routine only supports non-interleaved vectors.
    !
    ! If the FEM-space is scalar valued, exactly as many components must be
    ! specified as there are dimensions.
    !
    ! If the FEM-space is vector valued, exactly one function must be
    ! provided which must have exactly as many components as dimensions.
!</description>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Function to compute the divergence from.
    type(t_fev2VectorData), intent(in) :: rfunction
!</input>

!<output>
    ! Returns the value of the integral
    real(DP), intent(out) :: dintvalue
!</output>    

!</subroutine>

    ! Local variables
    real(DP) :: dderivX,dderivY,dderivZ
    integer :: iel, icubp, ndim, idimfe
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:,:), pointer :: p_DderivXVec,p_DderivYVec,p_DderivZVec
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    
    dintvalue = 0.0_DP
    
    ! Dimension of the spaces
    ndim = rintAssembly%p_rtriangulation%ndim

    ! Must be vector-valued.
    if (.not. rfunction%bisVectorField) return

    ! -----------------------
    ! Vector-valued function.
    ! FE space or multiple components.
    ! -----------------------

    ! Cancel if the number of vectors does not match the dimension.
    ! Each vector corresponds to one dimension.
    if (rfunction%ndimVectorField .ne. ndim) then
      call output_line ("Parameters wrong.",&
          OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_divergenceL2norm")
      call sys_halt()
    end if

    select case (ndim)
    
    ! ----------------
    ! 1D
    ! ----------------
    case (NDIM1D)
    
      ! Get the data array with the values of the FEM function
      ! in the cubature points
      p_DderivXVec => rfunction%p_DdataVec(:,:,:,DER_DERIV1D_X)

      ! Loop over the elements in the current set.
      do iel = 1,nelements

        ! Loop over all cubature points on the current element
        do icubp = 1,npointsPerElement

          ! Get the value of the FEM function
          dderivX = p_DderivXVec(1,icubp,iel)
          
          ! Multiply the values by the cubature weight and sum up
          ! into the integral value
          dintvalue = dintvalue + p_DcubWeight(icubp,iel) * dderivX**2
            
        end do ! icubp
      
      end do ! iel

    ! ----------------
    ! 2D
    ! ----------------
    case (NDIM2D)
    
      ! Get the data array with the values of the FEM function
      ! in the cubature points
      p_DderivXVec => rfunction%p_DdataVec(:,:,:,DER_DERIV2D_X)
      p_DderivYVec => rfunction%p_DdataVec(:,:,:,DER_DERIV2D_Y)

      ! Loop over the elements in the current set.
      do iel = 1,nelements

        ! Loop over all cubature points on the current element
        do icubp = 1,npointsPerElement

          ! Loop over the dimensions of the FE space
          do idimfe = 1,rfunction%ndimVectorField
            
            ! Get the value of the FEM function
            dderivX = p_DderivXVec(1,icubp,iel)
            dderivY = p_DderivYVec(2,icubp,iel)
            
            ! Multiply the values by the cubature weight and sum up
            ! into the integral value
            dintvalue = dintvalue + p_DcubWeight(icubp,iel) * &
                ( dderivX + dderivY )**2
                
          end do ! idimfe
            
        end do ! icubp
      
      end do ! iel

    ! ----------------
    ! 3D
    ! ----------------
    case (NDIM3D)
    
      ! Get the data array with the values of the FEM function
      ! in the cubature points
      p_DderivXVec => rfunction%p_DdataVec(:,:,:,DER_DERIV3D_X)
      p_DderivYVec => rfunction%p_DdataVec(:,:,:,DER_DERIV3D_Y)
      p_DderivZVec => rfunction%p_DdataVec(:,:,:,DER_DERIV3D_Z)

      ! Loop over the elements in the current set.
      do iel = 1,nelements

        ! Loop over all cubature points on the current element
        do icubp = 1,npointsPerElement

          ! Loop over the dimensions of the FE space
          do idimfe = 1,rfunction%ndimVectorField
          
            ! Get the value of the FEM function
            dderivX = p_DderivXVec(1,icubp,iel)
            dderivY = p_DderivYVec(2,icubp,iel)
            dderivZ = p_DderivZVec(3,icubp,iel)
            
            ! Multiply the values by the cubature weight and sum up
            ! into the integral value
            dintvalue = dintvalue + p_DcubWeight(icubp,iel) * &
                ( dderivX + dderivY + dderivZ )**2
          
          end do
            
        end do ! icubp
      
      end do ! iel
      
    end select

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_divergenceL2norm(dintvalue,rassemblyData,rintAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the (squared) L2-norm of the divergence of an arbitrary finite 
    ! element vector field <tex> $ ||div v||^2_{L_2} $ </tex>.
    !
    ! The vector field including first derivatives must be provided in 
    ! revalVectors. The routine only supports vector-valued, non-interleaved vectors:
    ! - If the FEM space is scalar-valued, revalVectors must be set up
    !   using fev2_addVectorFieldToEvalList.
    ! - If the FEM space is vector-valued, revalVectors must be set up
    !   using fev2_addVectorToEvalList providing the complete vector field
    !   in one FE function.
!</description>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaIntegralAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaIntegralAssembly), intent(in) :: rintAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>

!<output>
    ! Returns the value of the integral
    real(DP), intent(out) :: dintvalue
!</output>    

!</subroutine>

    call bma_docalc_divergenceL2norm(dintvalue,rassemblyData,rintAssembly,&
      npointsPerElement,nelements,revalVectors%p_RvectorData(1))
      
  end subroutine

  ! ***************************************************************************
  ! ***************************************************************************
  ! ***************************************************************************

!<subroutine>

  subroutine bma_auxGetLocalMeshWidth_area (Ielements,rtriangulation,DmeshWidth)

!<description>
  ! For a set of elements Ielements, this routine calculates a local mesh width
  ! using the area of the cells. Area/Volume-based approach.
  !
  ! Note: All elements must be of the same type (triangular, quad,...).
!</description>

!<input>
  ! List of elements.
  integer, dimension(:), intent(in) :: Ielements
  
  ! Underlying triangulation
  type(t_triangulation), intent(in) :: rtriangulation
!</input>

!<output>
  ! Mesh with of all elements Ielements.
  real(DP), dimension(:), intent(out) :: DmeshWidth
!</output>

!</subroutine>

    ! local variables
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    real(DP), dimension(:), pointer :: p_DelementVolume
    integer :: nve,iel
    
    if (size(Ielements) .eq. 0) return
    
    ! Get some information about the mesh.
    call storage_getbase_int2d(rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
    call storage_getbase_double(rtriangulation%h_DelementVolume,p_DelementVolume)
    
    ! Get the number of vertices on the element.
    do nve = ubound(p_IverticesAtElement,1),1,-1
      if (p_IverticesAtElement(nve,Ielements(1)) .ne. 0) exit
    end do

    ! What is the current dimension
    select case (rtriangulation%ndim)
    case (NDIM1D)
      do iel = 1,size(Ielements)
        ! Length of the interval
        DmeshWidth(iel) = p_DelementVolume(iel)
      end do
      
    case (NDIM2D)
      
      ! Type of the element
      select case (nve)
      
      ! Triangular cell
      case (3)
        
        do iel = 1,size(Ielements)
          ! Calculate the local h from the area of the element.
          ! As the element is a tri, multiply the volume by 2.
          DmeshWidth(iel) = sqrt(2.0_DP*p_DelementVolume(iel))
        end do
      
      ! Quadrilateral cell
      case (4)
      
        do iel = 1,size(Ielements)
          ! Calculate the local h from the area of the element.
          DmeshWidth(iel) = sqrt(p_DelementVolume(iel))
        end do

      end select

    case (NDIM3D)
      
      ! Type of the element
      select case (nve)
      
      ! Tetrahedral cell
      case (6)
        
        do iel = 1,size(Ielements)
          ! Calculate the local h from the area of the element.
          ! As the element is a tetra, multiply the volume by 3.
          DmeshWidth(iel) = sqrt(3.0_DP*p_DelementVolume(iel))
        end do
      
      ! Hexahedral cell
      case (7)
      
        do iel = 1,size(Ielements)
          ! Calculate the local h from the area of the element.
          DmeshWidth(iel) = sqrt(p_DelementVolume(iel))
        end do

      end select

    end select
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine bma_auxGetLocalMeshWidth_ray (&
      Ielements,Ddirection,rtriangulation,DmeshWidth)

!<description>
  ! For a set of elements Ielements, this routine calculates a local mesh width
  ! using the area of the cells. Ray-shooting approach.
  !
  ! Note: All elements must be of the same type (triangular, quad,...).
!</description>

!<input>
  ! List of elements.
  integer, dimension(:), intent(in) :: Ielements
  
  ! Underlying triangulation
  type(t_triangulation), intent(in) :: rtriangulation

  ! For every cubature point and element, the shooting direction 
  ! by which the element size is calculated. The direction does not have
  ! to be normalised.
  real(DP), dimension(:,:), intent(in) :: Ddirection
!</input>

!<output>
  ! Mesh with of all elements Ielements.
  real(DP), dimension(:), intent(out) :: DmeshWidth
!</output>

!</subroutine>

    ! local variables
    integer, dimension(:,:), pointer :: p_IverticesAtElement
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    integer :: nve,iel,i,ielement,ivt
    real(DP) :: dlen
    real(DP), dimension(2,4) :: Dcorners2D
    
    if (size(Ielements) .eq. 0) return
    
    ! Get some information about the mesh.
    call storage_getbase_int2d(rtriangulation%h_IverticesAtElement,p_IverticesAtElement)
    call storage_getbase_double2d(rtriangulation%h_DvertexCoords,p_DvertexCoords)
    
    ! Get the number of vertices on the element.
    do nve = ubound(p_IverticesAtElement,1),1,-1
      if (p_IverticesAtElement(nve,Ielements(1)) .ne. 0) exit
    end do

    ! What is the current dimension
    select case (rtriangulation%ndim)
      
    case (NDIM2D)
      
      ! Type of the element
      select case (nve)
     
      ! Quadrilateral cell
      case (4)
      
        do iel = 1,size(Ielements)
          ! Get the element corners.
          ielement = Ielements(iel)
          do i=1,4
            ivt = p_IverticesAtElement(i,ielement)
            Dcorners2D(1,i) = p_DvertexCoords(1,ivt)
            Dcorners2D(2,i) = p_DvertexCoords(2,ivt)
          end do
          
          ! Length of the shoting direction
          dlen = 1.0_DP/(sqrt(Ddirection(1,iel)**2 + Ddirection(2,iel)**2))
        
          ! Shoot the ray.
          call gaux_getDirectedExtentQuad2D (Dcorners2D,&
              Ddirection(1,iel)*dlen,Ddirection(2,iel)*dlen,DmeshWidth(iel))
        end do
        
      ! Fallback.
      case default
        call bma_auxGetLocalMeshWidth_area (Ielements,rtriangulation,DmeshWidth)

      end select

    case (NDIM3D)
      
      ! Type of the element
      select case (nve)
      
      ! Fallback.
      case default
        call bma_auxGetLocalMeshWidth_area (Ielements,rtriangulation,DmeshWidth)

      end select

    ! Fallback.
    case default
      call bma_auxGetLocalMeshWidth_area (Ielements,rtriangulation,DmeshWidth)

    end select
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine bma_auxGetMeanVector (DmeanVector,&
      DcubWeights,npointsPerElement,nelements,du1,du2,du3,Dvector)

!<description>
  ! Calculates the mean vector in every element.
!</description>

!<input>
  ! For every cubature point and every element, a cubature weight
  ! to calculate integrals on the element.
  real(DP), dimension(:,:), intent(in) :: DcubWeights
  
  ! Number of points per element
  integer, intent(in) :: npointsPerElement

  ! Number of elements
  integer, intent(in) :: nelements
  
  ! OPTIONAL: A constant velocity field. Can be omitted if Dvector is specified.
  real(DP), intent(in), optional :: du1,du2,du3
  
  ! OPTIONAL: Array with the values of the velocity field. Specifies for all points on
  ! all elements velocity.
  ! dimension(#dim, npointsPerElement,nelements)
  real(DP), dimension(:,:,:), intent(in), optional :: Dvector
!</input>

!<output>
  ! Out: The mean velocity in every element.
  real(DP), dimension(:,:), intent(out) :: DmeanVector
!</output>

!</subroutine>

    ! local variables
    integer :: iel,icubp,ndim
    real(DP) :: dvelX,dvelY,dvelZ,dvolume,domega
    
    dvelX = 0.0_DP
    dvelY = 0.0_DP
    dvelZ = 0.0_DP
    if (present(du1)) dvelX = du1
    if (present(du2)) dvelX = du2
    if (present(du3)) dvelX = du3
    
    ! Dimension
    ndim = ubound(Dvector,1)
    
    if (.not. present (Dvector)) then
    
      ! Constant velocity. Nothing to calculate.
    
      select case (ndim)
      case (NDIM1D)
    
        ! Loop over the elements
        do iel=1,nelements

          ! Calculate the mean velocity
          DmeanVector(1,iel) = dvelX

        end do
          
      case (NDIM2D)
    
        ! Loop over the elements
        do iel=1,nelements
        
          DmeanVector(1,iel) = dvelX
          DmeanVector(2,iel) = dvelY
          
        end do

      case (NDIM3D)
    
        ! Loop over the elements
        do iel=1,nelements
        
          DmeanVector(1,iel) = dvelX
          DmeanVector(2,iel) = dvelY
          DmeanVector(3,iel) = dvelZ

        end do
            
      end select
      
    else          

      select case (ndim)
      case (NDIM1D)
    
        ! Loop over the elements
        do iel=1,nelements
        
          dvelX = 0.0_DP
          dvolume = 0.0_DP
          
          ! Loop over the cubature points, calculate the element size and
          ! the integral of the velocity.
          do icubp = 1,npointsPerElement
            domega = DcubWeights(icubp,iel)
            dvolume = dvolume + domega
            dvelX = dvelX + domega*Dvector(1,icubp,iel)
          end do
          
          ! Calculate the mean velocity
          DmeanVector(1,iel) = dvelX/dvolume

        end do
          
      case (NDIM2D)
    
        ! Loop over the elements
        do iel=1,nelements
        
          dvelX = 0.0_DP
          dvelY = 0.0_DP
          dvolume = 0.0_DP
          
          ! Loop over the cubature points, calculate the element size and
          ! the integral of the velocity.
          do icubp = 1,npointsPerElement
            domega = DcubWeights(icubp,iel)
            dvolume = dvolume + domega
            dvelX = dvelX + domega*Dvector(1,icubp,iel)
            dvelY = dvelZ + domega*Dvector(2,icubp,iel)
          end do
          
          ! Calculate the mean velocity
          DmeanVector(1,iel) = dvelX/dvolume
          DmeanVector(2,iel) = dvelY/dvolume

        end do

      case (NDIM3D)
    
        ! Loop over the elements
        do iel=1,nelements
        
          dvelX = 0.0_DP
          dvelY = 0.0_DP
          dvelZ = 0.0_DP
          dvolume = 0.0_DP
          
          ! Loop over the cubature points, calculate the element size and
          ! the integral of the velocity.
          do icubp = 1,npointsPerElement
            domega = DcubWeights(icubp,iel)
            dvolume = dvolume + domega
            dvelX = dvelX + domega*Dvector(1,icubp,iel)
            dvelY = dvelZ + domega*Dvector(2,icubp,iel)
            dvelZ = dvelY + domega*Dvector(3,icubp,iel)
          end do
          
          ! Calculate the mean velocity
          DmeanVector(1,iel) = dvelX/dvolume
          DmeanVector(2,iel) = dvelY/dvolume
          DmeanVector(3,iel) = dvelZ/dvolume

        end do
            
      end select

    end if
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine bma_auxGetMeanScalar (DmeanScalar,&
      DcubWeights,npointsPerElement,nelements,ddata,Dscalar)

!<description>
  ! Calculates the mean data in every element.
!</description>

!<input>
  ! For every cubature point and every element, a cubature weight
  ! to calculate integrals on the element.
  real(DP), dimension(:,:), intent(in) :: DcubWeights
  
  ! Number of points per element
  integer, intent(in) :: npointsPerElement

  ! Number of elements
  integer, intent(in) :: nelements
  
  ! OPTIONAL: A constant data field. Can be omitted if Dscalar is specified.
  real(DP), intent(in), optional :: ddata
  
  ! OPTIONAL: Array with the values of the velocity field. Specifies for all points on
  ! all elements velocity.
  ! dimension(npointsPerElement,nelements)
  real(DP), dimension(:,:), intent(in), optional :: Dscalar
!</input>

!<output>
  ! Out: The mean velocity in every element.
  real(DP), dimension(:), intent(out) :: DmeanScalar
!</output>

!</subroutine>

    ! local variables
    integer :: iel,icubp
    real(DP) :: da,dvolume,domega
    
    if (.not. present (Dscalar)) then
    
      ! Constant data. Nothing to calculate.
    
      ! Loop over the elements
      do iel=1,nelements

        ! Calculate the mean velocity
        DmeanScalar(iel) = ddata

      end do
          
    else          

      ! Loop over the elements
      do iel=1,nelements
      
        da = 0.0_DP
        dvolume = 0.0_DP
        
        ! Loop over the cubature points, calculate the element size and
        ! the integral of the velocity.
        do icubp = 1,npointsPerElement
          domega = DcubWeights(icubp,iel)
          dvolume = dvolume + domega
          da = da + domega*Dscalar(icubp,iel)
        end do
        
        ! Calculate the mean velocity
        DmeanScalar(iel) = da/dvolume

      end do
      
    end if
        
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine bma_auxGetLDelta (Ddelta,&
      cstabiltype,dupsam,duMaxR,DcubWeights,DmeshWidth,npointsPerElement,nelements,&
      DvelocityMean,dnuConst,Dnu,duMin)

!<description>
  ! This routine calculates a local ddelta=DELTA_T for a set of finite
  ! elements. This can be used by the streamline diffusion
  ! stabilisation technique as a multiplier of the (local) bilinear form.
  !
  ! Element independent version which uses the volume of the elements
  ! for the computation.
  !
  ! Used for vector-valued elements.
!</description>

!<input>
  ! Reciprocal of the maximum norm of velocity in the domain:
  ! 1/duMaxR = 1/||u||_max(Omega)
  real(DP), intent(in) :: duMaxR

  ! For every cubature point and every element, a cubature weight
  ! to calculate integrals on the element.
  real(DP), dimension(:,:), intent(in) :: DcubWeights

  ! Local h (mesh width) of each element
  real(DP), dimension(:), intent(in) :: DmeshWidth

  ! Type of SD method to apply.
  ! = 0: Use simple SD stabilisation:
  !      ddelta = dupsam * h_T.
  ! = 1: Use Samarskji SD stabilisation; usually dupsam = 0.1 .. 2.
  !      ddelta = dupsam * h_t/||u||_T * 2*Re_T/(1+Re_T)
  integer, intent(in) :: cstabilType

  ! User defined parameter for configuring the streamline diffusion.
  real(DP), intent(in) :: dupsam

  ! Number of points per element
  integer, intent(in) :: npointsPerElement

  ! Number of elements
  integer, intent(in) :: nelements
  
  ! Array with the values of the velocity field. Specifies for all elements
  ! the mean velocity in that element.
  ! dimension(#dim,#elements)
  real(DP), dimension(:,:), intent(in) :: DvelocityMean

  ! OPTIONAL: Viscosity coefficient. Can be omitted if Dnu is specified.
  real(DP), intent(in), optional :: dnuConst

  ! OPTIONAL: Viscosity coefficient in all cubature points on all selement
  real(DP), dimension(:,:), intent(in), optional :: Dnu

  ! OPTINOAL: Minimum velocity from which stabilisation is applied.
  ! For a velocity smaller than this, no stabilisation is applied
  ! on the corresponding element. Default = 1E-8.
  real(DP), intent(in), optional :: duMin
!</input>

!<output>
  ! Out: local Ddelta on all elements
  real(DP), dimension(:), intent(out) :: Ddelta
!</output>

!</subroutine>

  ! local variables
  real(DP) :: dlocalH,dunorm,dreLoc,dnuRec
  integer :: iel,icubp
  real(DP) :: domega,dminVelo,dvolume
  logical :: bconstNu
  integer :: ndim
  
    ndim = ubound (DvelocityMean,1)
  
    ! Minimum velocity
    dminVelo = 1E-8_DP
    if (present(duMin)) dminVelo = duMin
    
    ! Some preparations if some things are constant.
    bconstNu = .not. present (Dnu)
    
    if (bconstNu) then
    
      ! =================================================
      ! Nonconstant velocity, constant viscosity
      ! =================================================

      ! Calculate the mean viscosity coefficient -- or more precisely,
      ! its reciprocal.
      dnuRec = 1.0_DP/dnuConst

      ! Calculate ddelta... (cf. p. 121 in Turek`s CFD book)
      select case (cstabiltype)
      
      ! =====================================================
      ! Simple calculation of ddelta.
      ! Viscosity not used.
      case (0)
      
        select case (ndim)
        
        case (NDIM1D)

          ! Loop through all elements
          do iel = 1,nelements

            ! Calculate the norm of the mean local velocity.
            dunorm = abs(DvelocityMean(1,iel))

            ! Now we have:   dunorm = ||u||_T
            ! and:           u_T = a1*u1_T + a2*u2_T

            ! If the norm of the velocity is small, we choose ddelta = 0,
            ! which results in central difference in the streamline diffusion
            ! matrix assembling:

            if (dunorm .le. dminVelo) then

              Ddelta(iel) = 0.0_DP

            else

              ! Calculate ddelta... (cf. p. 121 in Turek`s CFD book)
              Ddelta(iel) = dupsam*DmeshWidth(iel)

            end if

          end do

        case (NDIM2D)

          ! Loop through all elements
          do iel = 1,nelements

            ! Calculate the norm of the mean local velocity.
            dunorm = sqrt(DvelocityMean(1,iel)**2 + DvelocityMean(2,iel)**2)

            ! Now we have:   dunorm = ||u||_T
            ! and:           u_T = a1*u1_T + a2*u2_T

            ! If the norm of the velocity is small, we choose ddelta = 0,
            ! which results in central difference in the streamline diffusion
            ! matrix assembling:

            if (dunorm .le. dminVelo) then

              Ddelta(iel) = 0.0_DP

            else

              ! Calculate ddelta... (cf. p. 121 in Turek`s CFD book)
              Ddelta(iel) = dupsam*DmeshWidth(iel)

            end if

          end do
          
        case (NDIM3D)

          ! Loop through all elements
          do iel = 1,nelements

            ! Calculate the norm of the mean local velocity.
            dunorm = sqrt(DvelocityMean(1,iel)**2 + &
                          DvelocityMean(2,iel)**2 + DvelocityMean(3,iel)**2)

            ! Now we have:   dunorm = ||u||_T
            ! and:           u_T = a1*u1_T + a2*u2_T

            ! If the norm of the velocity is small, we choose ddelta = 0,
            ! which results in central difference in the streamline diffusion
            ! matrix assembling:

            if (dunorm .le. dminVelo) then

              Ddelta(iel) = 0.0_DP

            else

              ! Calculate ddelta... (cf. p. 121 in Turek`s CFD book)
              Ddelta(iel) = dupsam*dlocalH

            end if

          end do

        end select

      ! =====================================================
      ! Standard Samarskji-like calculation
      case (1)

        select case (ndim)
        case (NDIM1D)

          ! Loop through all elements
          do iel = 1,nelements

            ! Calculate the norm of the mean local velocity.
            dunorm = abs(DvelocityMean(1,iel))

            ! Now we have:   dunorm = ||u||_T
            ! and:           u_T = a1*u1_T + a2*u2_T

            ! If the norm of the velocity is small, we choose ddelta = 0,
            ! which results in central difference in the streamline diffusion
            ! matrix assembling:

            if (dunorm .le. dminVelo) then

              Ddelta(iel) = 0.0_DP

            else

              ! At first calculate the local Reynolds number
              ! RELOC = Re_T = ||u||_T * h_T / NU

              dreLoc = dunorm*DmeshWidth(iel)*dnuRec

              ! and then the ddelta = UPSAM * h_t/||u|| * 2*Re_T/(1+Re_T)

              Ddelta(iel) = dupsam * dlocalH*duMaxR * 2.0_DP*(dreLoc/(1.0_DP+dreLoc))

            end if

          end do

        case (NDIM2D)

          ! Loop through all elements
          do iel = 1,nelements

            ! Calculate the norm of the mean local velocity.
            dunorm = sqrt(DvelocityMean(1,iel)**2 + DvelocityMean(2,iel)**2)

            ! Now we have:   dunorm = ||u||_T
            ! and:           u_T = a1*u1_T + a2*u2_T

            ! If the norm of the velocity is small, we choose ddelta = 0,
            ! which results in central difference in the streamline diffusion
            ! matrix assembling:

            if (dunorm .le. dminVelo) then

              Ddelta(iel) = 0.0_DP

            else

              ! At first calculate the local Reynolds number
              ! RELOC = Re_T = ||u||_T * h_T / NU

              dreLoc = dunorm*DmeshWidth(iel)*dnuRec

              ! and then the ddelta = UPSAM * h_t/||u|| * 2*Re_T/(1+Re_T)

              Ddelta(iel) = dupsam * dlocalH*duMaxR * 2.0_DP*(dreLoc/(1.0_DP+dreLoc))

            end if

          end do
          
        case (NDIM3D)

          ! Loop through all elements
          do iel = 1,nelements

            ! Calculate the norm of the mean local velocity.
            dunorm = sqrt(DvelocityMean(1,iel)**2 + &
                          DvelocityMean(2,iel)**2 + DvelocityMean(3,iel)**2)

            ! Now we have:   dunorm = ||u||_T
            ! and:           u_T = a1*u1_T + a2*u2_T

            ! If the norm of the velocity is small, we choose ddelta = 0,
            ! which results in central difference in the streamline diffusion
            ! matrix assembling:

            if (dunorm .le. dminVelo) then

              Ddelta(iel) = 0.0_DP

            else

              ! At first calculate the local Reynolds number
              ! RELOC = Re_T = ||u||_T * h_T / NU

              dreLoc = dunorm*DmeshWidth(iel)*dnuRec

              ! and then the ddelta = UPSAM * h_t/||u|| * 2*Re_T/(1+Re_T)

              Ddelta(iel) = dupsam * dlocalH*duMaxR * 2.0_DP*(dreLoc/(1.0_DP+dreLoc))

            end if

          end do
          
        end select

      end select

    else
    
      ! =================================================
      ! Nonconstant velocity, nonconstant viscosity
      ! =================================================

      ! Calculate ddelta... (cf. p. 121 in Turek`s CFD book)
      select case (cstabiltype)
      
      ! =====================================================
      ! Simple calculation of ddelta
      case (0)
      
        select case (ndim)
        
        case (NDIM1D)

          ! Loop through all elements
          do iel = 1,nelements

            ! Calculate the norm of the mean local velocity.
            dunorm = abs(DvelocityMean(1,iel))

            ! If the norm of the velocity is small, we choose ddelta = 0,
            ! which results in central difference in the streamline diffusion
            ! matrix assembling:

            if (dunorm .le. dminVelo) then

              Ddelta(iel) = 0.0_DP

            else

              ! Calculate ddelta... (cf. p. 121 in Turek`s CFD book)
              Ddelta(iel) = dupsam*DmeshWidth(iel)

            end if

          end do

        case (NDIM2D)

          ! Loop through all elements
          do iel = 1,nelements

            ! Calculate the norm of the mean local velocity.
            dunorm = sqrt(DvelocityMean(1,iel)**2 + DvelocityMean(2,iel)**2)

            ! If the norm of the velocity is small, we choose ddelta = 0,
            ! which results in central difference in the streamline diffusion
            ! matrix assembling:

            if (dunorm .le. dminVelo) then

              Ddelta(iel) = 0.0_DP

            else

              ! Calculate ddelta... (cf. p. 121 in Turek`s CFD book)
              Ddelta(iel) = dupsam*DmeshWidth(iel)

            end if

          end do
          
        case (NDIM3D)

          ! Loop through all elements
          do iel = 1,nelements

            ! Calculate the norm of the mean local velocity.
            dunorm = sqrt(DvelocityMean(1,iel)**2 + &
                          DvelocityMean(2,iel)**2 + DvelocityMean(3,iel)**2)

            ! Now we have:   dunorm = ||u||_T
            ! and:           u_T = a1*u1_T + a2*u2_T

            ! If the norm of the velocity is small, we choose ddelta = 0,
            ! which results in central difference in the streamline diffusion
            ! matrix assembling:

            if (dunorm .le. dminVelo) then

              Ddelta(iel) = 0.0_DP

            else

              ! Calculate ddelta... (cf. p. 121 in Turek`s CFD book)
              Ddelta(iel) = dupsam*dlocalH

            end if

          end do

        end select

      ! =====================================================
      ! Standard Samarskji-like calculation
      case (1)

        select case (ndim)
        case (NDIM1D)

          ! Loop through all elements
          do iel = 1,nelements

            ! Loop through the cubature points on the current element
            ! and calculate the mean velocity there.
            dnuRec = 0.0_DP
            do icubp = 1,npointsPerElement
              domega = DcubWeights(icubp,iel)
              dvolume = dvolume + domega
              dnuRec = dnuRec + domega*Dnu(icubp,iel)
            end do

            ! Calculate the mean viscosity coefficient -- or more precisely,
            ! its reciprocal.
            dnuRec = dvolume / dnuRec

            ! Calculate the norm of the mean local velocity.
            dunorm = abs(DvelocityMean(1,iel))

            ! Now we have:   dunorm = ||u||_T
            ! and:           u_T = a1*u1_T + a2*u2_T

            ! If the norm of the velocity is small, we choose ddelta = 0,
            ! which results in central difference in the streamline diffusion
            ! matrix assembling:

            if (dunorm .le. dminVelo) then

              Ddelta(iel) = 0.0_DP

            else

              ! At first calculate the local Reynolds number
              ! RELOC = Re_T = ||u||_T * h_T / NU

              dreLoc = dunorm*DmeshWidth(iel)*dnuRec

              ! and then the ddelta = UPSAM * h_t/||u|| * 2*Re_T/(1+Re_T)

              Ddelta(iel) = dupsam * dlocalH*duMaxR * 2.0_DP*(dreLoc/(1.0_DP+dreLoc))

            end if

          end do

        case (NDIM2D)

          ! Loop through all elements
          do iel = 1,nelements

            ! Loop through the cubature points on the current element
            ! and calculate the mean velocity there.
            dnuRec = 0.0_DP
            do icubp = 1,npointsPerElement
              domega = DcubWeights(icubp,iel)
              dvolume = dvolume + domega
              dnuRec = dnuRec + domega*Dnu(icubp,iel)
            end do

            ! Calculate the mean viscosity coefficient -- or more precisely,
            ! its reciprocal.
            dnuRec = dvolume / dnuRec

            ! Calculate the norm of the mean local velocity.
            dunorm = sqrt(DvelocityMean(1,iel)**2 + DvelocityMean(2,iel)**2)

            ! Now we have:   dunorm = ||u||_T
            ! and:           u_T = a1*u1_T + a2*u2_T

            ! If the norm of the velocity is small, we choose ddelta = 0,
            ! which results in central difference in the streamline diffusion
            ! matrix assembling:

            if (dunorm .le. dminVelo) then

              Ddelta(iel) = 0.0_DP

            else

              ! At first calculate the local Reynolds number
              ! RELOC = Re_T = ||u||_T * h_T / NU

              dreLoc = dunorm*DmeshWidth(iel)*dnuRec

              ! and then the ddelta = UPSAM * h_t/||u|| * 2*Re_T/(1+Re_T)

              Ddelta(iel) = dupsam * dlocalH*duMaxR * 2.0_DP*(dreLoc/(1.0_DP+dreLoc))

            end if

          end do
          
        case (NDIM3D)

          ! Loop through all elements
          do iel = 1,nelements

            ! Loop through the cubature points on the current element
            ! and calculate the mean velocity there.
            dnuRec = 0.0_DP
            do icubp = 1,npointsPerElement
              domega = DcubWeights(icubp,iel)
              dvolume = dvolume + domega
              dnuRec = dnuRec + domega*Dnu(icubp,iel)
            end do

            ! Calculate the mean viscosity coefficient -- or more precisely,
            ! its reciprocal.
            dnuRec = dvolume / dnuRec

            ! Calculate the norm of the mean local velocity.
            dunorm = sqrt(DvelocityMean(1,iel)**2 + &
                          DvelocityMean(2,iel)**2 + DvelocityMean(3,iel)**2)

            ! Now we have:   dunorm = ||u||_T
            ! and:           u_T = a1*u1_T + a2*u2_T

            ! If the norm of the velocity is small, we choose ddelta = 0,
            ! which results in central difference in the streamline diffusion
            ! matrix assembling:

            if (dunorm .le. dminVelo) then

              Ddelta(iel) = 0.0_DP

            else

              ! At first calculate the local Reynolds number
              ! RELOC = Re_T = ||u||_T * h_T / NU

              dreLoc = dunorm*DmeshWidth(iel)*dnuRec

              ! and then the ddelta = UPSAM * h_t/||u|| * 2*Re_T/(1+Re_T)

              Ddelta(iel) = dupsam * dlocalH*duMaxR * 2.0_DP*(dreLoc/(1.0_DP+dreLoc))

            end if

          end do
          
        end select

      end select

    end if
    
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine bma_auxViscosityTensorNorm2D (cviscoTensorType, DviscoTensorNorm, &
      Du1x, Du1y, Du2x, Du2y)

!<description>
  ! Calculates the norm of a viscosity tensor,
  !    z := D(u):D(u) = ||D(u)||^2
  ! The calculation takes place in every cubature point on every element.
  ! ctensorType defines the type of the tensor.
!</description>

!<input>
  ! Type of the viscosity tensor:
  ! =0:  D(u)=grad(u)
  ! =1:  D(u)=1/2 ( grad(u) + grad(u)^T )
  integer, intent(in) :: cviscoTensorType
  
  ! The x/y derivative of the u1 and u2 in every cubature point on every
  ! element.
  real(DP), dimension(:,:), intent(in) :: Du1x
  real(DP), dimension(:,:), intent(in) :: Du1y
  real(DP), dimension(:,:), intent(in) :: Du2x
  real(DP), dimension(:,:), intent(in) :: Du2y
!</input>

!<output>
  ! Output array. Receives for every point on every element the norm of
  ! the tensor.
  real(DP), dimension(:,:), intent(out) :: DviscoTensorNorm
!</output>

!</subroutine>

    ! local variables
    integer :: npointsPerElement, nelements, i, j
    
    npointsPerElement = ubound (DviscoTensorNorm,1)
    nelements = ubound(DviscoTensorNorm,2)
    
    ! Calculate ||D(u)||^2 to p_Ddata(:,:)
    select case (cviscoTensorType)
    case (0)
      ! D(u) = grad(u)
      do i=1,nelements
        do j=1,npointsPerElement
          DviscoTensorNorm(j,i) = Du1x(j,i)**2 + Du1y(j,i)**2 + &
              Du2x(j,i)**2 + Du2y(j,i)**2
        end do
      end do
      
    case (1)
      ! D(u) = 1/2 ( grad(u) + grad(u)^T )
      do i=1,nelements
        do j=1,npointsPerElement
          DviscoTensorNorm(j,i) = (Du1x(j,i)**2 + &
              0.5_DP * (Du1y(j,i) + Du2x(j,i))**2 + &
              Du2y(j,i)**2)
        end do
      end do
      
    case default
    
      DviscoTensorNorm(:,:) = 0.0_DP
      
    end select
        
  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine bma_auxFluidViscosity (cviscoModel, Dviscosity, &
      dnu,dviscoYield,dviscoEps,dviscoExponent,DviscoTensorNorm)

!<description>
  ! Calculates a nonlinear viscosity.
!</description>

!<input>
  ! Type of the viscosity model:
  ! =0:  Constant viscosity:
  !        nu(x) = dnu 
  ! =1:  Power law:
  !         nu(x) = dnu * z^(dviscoexponent/2 - 1),
  !         z = ||D(u)||^2+dviscoEps
  ! =2:  Bingham fluid:
  !         nu(x) = dnu + sqrt(2)/2 * dviscoyield / sqrt(|D(u)||^2+dviscoEps^2),
  ! =3:  General viscoplastic fluid:
  !         nu(x) = dnu + sqrt(2)/2 * dviscoyield * z^(dviscoexponent/2 - 1),
  !         z = ||D(u)||^2 + dviscoEps^2
  integer, intent(in) :: cviscoModel
  
  ! Viscosity constant
  real(DP), intent(in) :: dnu
  
  ! Additional parameters for the viscosity.
  real(DP), intent(in) :: dviscoYield
  real(DP), intent(in) :: dviscoEps
  real(DP), intent(in) :: dviscoExponent
  
  ! Norm of the viscosity tensor. Must be specified if cviscoModel<>0.
  real(DP), dimension(:,:), intent(in) :: DviscoTensorNorm
!</input>

!<output>
  ! Output array. Receives for every point on every element the viscosity.
  real(DP), dimension(:,:), intent(out) :: Dviscosity
!</output>

!</subroutine>

    ! local variables
    integer :: npointsPerElement, nelements, i, j
    
    npointsPerElement = ubound (Dviscosity,1)
    nelements = ubound(Dviscosity,2)

    ! Calculate the viscosity.
    select case (cviscoModel)
    case (1)
      ! Power law:
      !   nu = dnu * z^(dviscoexponent/2 - 1),
      !   z = ||D(u)||^2+dviscoEps
      do i=1,nelements
        do j = 1,npointsPerElement
          Dviscosity(j,i) = &
            dnu * (DviscoTensorNorm(j,i)+dviscoEps)**(0.5_DP*dviscoExponent - 1.0_DP)
        end do
      end do

    case (2)
      ! Bingham fluid:
      !   nu = dnu + sqrt(2)/2 * dviscoyield / sqrt(|D(u)||^2+dviscoEps^2),
      do i=1,nelements
        do j = 1,npointsPerElement
          Dviscosity(j,i) = dnu + 0.5_DP * sqrt(2.0_DP) * &
              dviscoyield / sqrt(DviscoTensorNorm(j,i) + dviscoEps**2)
        end do
      end do
      
    case (3)
      ! General viscoplastic fluid:
      !   nu = dnu + sqrt(2)/2 * dviscoyield * z^(dviscoexponent/2 - 1),
      !   z = ||D(u)||^2 + dviscoEps^2
      do i=1,nelements
        do j = 1,npointsPerElement
          Dviscosity(j,i) = dnu + 0.5_DP * sqrt(2.0_DP) * dviscoyield * &
              ( DviscoTensorNorm(j,i) + dviscoEps**2 )**( 0.5_DP * dviscoExponent - 1.0_DP)
        end do
      end do

    end select

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_docalc_streamlinediff(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,dscale,ix,iy,btensor,cstabiltype,clocalh,dupsam,dumax,&
      relemTempVectors,relemTempVectorsVec,du1,du2,du3,dviscosity,rvectorField,rvectorNu)

!<description>  
    ! Calculates the Streamline diffusion stabilisation at position (ix,iy)
    ! of the global matrix.
!</description>

!<inputoutput>
    ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
    ! have to be filled with data.
    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!</inputoutput>

!<input>
    ! Data necessary for the assembly. Contains determinants and
    ! cubature weights for the cubature,...
    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData

    ! Structure with all data about the assembly
    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly

    ! Number of points per element
    integer, intent(in) :: npointsPerElement

    ! Number of elements
    integer, intent(in) :: nelements

    ! Scaling factor of the SD term
    real(DP), intent(in) :: dscale
    
    ! Position of the matrix where to assemble the operator to
    integer, intent(in) :: ix,iy
    
    ! Defines whether the operator is assembled to the full tensor.
    ! =false: Operator is only assembled to the matrix at position (ix,iy).
    ! =true:  Operator is assembled to a tensor of the same dimension
    !         as the velocity field, starting at position (ix,iy) in the
    !         matrix.
    logical, intent(in) :: btensor
    
    ! Type of SD method to apply.
    ! = 0: Use simple SD stabilisation:
    !      ddelta = dupsam * h_T.
    ! = 1: Use Samarskji SD stabilisation; usually dupsam = 0.1 .. 2.
    !      ddelta = dupsam * h_t/||u||_T * 2*Re_T/(1+Re_T)
    integer, intent(in) :: cstabilType
    
    ! Method how to compute the local mesh width.
    ! = 0: Element volume based.
    ! = 1: Ray shooting method.
    integer, intent(in) :: clocalh
    
    ! Stabilisation parameter
    real(DP), intent(in) :: dupsam
    
    ! Maximum notm of the velocity over the complete domain.
    real(DP), intent(in) :: dumax

    ! Temporary vector structures. This structure must provide 2 temporary
    ! element-based vector buffers.
    type(t_fev2VectorData), intent(in) :: relemTempVectors

    ! Temporary vector structures. This structure must provide 1 temporary
    ! vector-valued element-based buffer. The buffer must have as many
    ! components as the underlying dimension of the space.
    type(t_fev2VectorData), intent(in) :: relemTempVectorsVec

    ! OPTIONAL: Constant velocity. Must be specified if rvectorField is not
    ! present. Ignored if rvectorField is present.
    real(DP), intent(in), optional :: du1,du2,du3

    ! OPTIONAL: Constant viscosity. Must be specified if rvectorNu is not
    ! present. Ignored if rvectorNu is present.
    real(DP), intent(in), optional :: dviscosity

    ! OPTIONAL: Evaluation structure that describes an underlying nonconstant
    ! velocity field. Can be omitted if du1/du2/du3 is given.
    type(t_fev2VectorData), intent(in), optional :: rvectorField

    ! OPTIONAL: Evaluation structure that describes an underlying nonconstant
    ! viscosity. Can be omitted if dviscosity is given.
    type(t_fev2VectorData), intent(in), optional :: rvectorNu
!</input>

!</subroutine>

    ! Local variables
    real(DP) :: dbasIx, dbasIy, dbasIz, dbasJx, dbasJy, dbasJz
    integer :: iel, icubp, idofe, jdofe, ndim
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial,p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaMatrixData), pointer :: p_rmatrixData
    real(DP), dimension(:,:), pointer :: p_Dnu
    real(DP), dimension(:), pointer :: p_DlocalDelta,p_DmeshWidth

    logical :: bconstNu, bconstVel
    real(DP), dimension(:,:,:), pointer :: p_Du
    real(DP), dimension(:,:), pointer :: p_DuMean
    real(DP) :: dnu,dvelX,dvelY,dvelZ, dsumI, dsumJ
    integer :: i,iend

    integer :: ndimfe
    real(DP) :: dumaxR

    ! Dimension of the underlying space.
    ndim = rmatrixAssembly%p_rtriangulation%ndim
    
    ! Get Dnu and the temp array for the computation of the local Delta.
    p_DlocalDelta => relemTempVectors%p_DcellData(:,1)
    p_DmeshWidth => relemTempVectors%p_DcellData(:,2)
    p_DuMean => relemTempVectorsVec%p_DcellDataVec(:,:,1)
    
    dumaxR = 1.0_DP/duMax
    
    bconstNu = present(rvectorNu)
    bconstVel = present(rvectorField)
    
    if (bconstVel) then
      dnu = dviscosity
    else
      p_Dnu => rvectorNu%p_Ddata(:,:,DER_FUNC)
    end if
    
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Get local data
    p_rmatrixData => RmatrixData(iy,ix)
    p_DbasTrial => RmatrixData(iy,ix)%p_DbasTrial
    p_DbasTest => RmatrixData(iy,ix)%p_DbasTest
    
    ! FE space dimension
    ndimfe = RmatrixData(iy,ix)%ndimfeTrial
    
    if (ndimfe .ne. RmatrixData(iy,ix)%ndimfeTest) then
      ! This does not make sense.
      call output_line ("Dimension of trial and test FE space different.",&
          OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_mass")
      call sys_halt()
    end if
    
    ! Interleaved matrix?
    if (.not. p_rmatrixData%bisInterleaved) then

      ! ---------------------------------------------------
      ! Compute the mean velocity in each element
      if (bconstVel) then

        dvelX = 0.0_DP
        dvelY = 0.0_DP
        dvelZ = 0.0_DP

        ! Constant velocity
        if (present(du1)) dvelX = du1
        if (present(du2)) dvelY = du2
        if (present(du3)) dvelZ = du3
        
        ! Calculate the mean velocity in every element
        call bma_auxGetMeanVector (p_DuMean,&
            p_DcubWeight,npointsPerElement,nelements,dvelX,dvelY,dvelZ)
        
      else
        p_Du => rvectorField%p_DdataVec(:,:,:,DER_FUNC)

        ! Calculate the mean velocity in every element
        call bma_auxGetMeanVector (p_DuMean,&
            p_DcubWeight,npointsPerElement,nelements,Dvector=p_Du)
        
      end if
      
      ! ---------------------------------------------------
      ! Calculate the local mesh width.
      select case (clocalh)
      case (0)
        call bma_auxGetLocalMeshWidth_area (&
            rassemblyData%p_IelementList,rmatrixAssembly%p_rtriangulation,p_DmeshWidth)
            
      case (1)
        if (.not. bconstVel) then
          ! Ray shooting currently only supported for nonconstant velocity fields.
          ! We use the mean velocity as shooting direction.
          call bma_auxGetLocalMeshWidth_ray (&
              rassemblyData%p_IelementList,p_DuMean,rmatrixAssembly%p_rtriangulation,p_DmeshWidth)
              
        else
          call bma_auxGetLocalMeshWidth_area (&
              rassemblyData%p_IelementList,rmatrixAssembly%p_rtriangulation,p_DmeshWidth)
        end if
      end select

      ! ---------------------------------------------------
      ! Compute the local delta
      if (bconstNu) then
        call bma_auxGetLDelta (p_DlocalDelta,&
              cstabiltype,dupsam,duMaxR,p_DcubWeight,&
              p_DmeshWidth,npointsPerElement,nelements,&
              p_DuMean,dnu)
      else
        call bma_auxGetLDelta (p_DlocalDelta,&
              cstabiltype,dupsam,duMaxR,p_DcubWeight,&
              p_DmeshWidth,npointsPerElement,nelements,&
              p_DuMean,Dnu=p_Dnu)
      end if

      if (ndimfe .eq. 1) then

        ! -----------------------
        ! Scalar-valued FE space. 
        ! -----------------------
        ! This IF command prevents an inner IF command and speeds 
        ! up the computation for all scalar standard FE spaces.
        
        ! Do we assemble the full tensor?
        if (btensor) then
          iend = ndim-1
        else
          iend = 0
        end if

        ! Set up the operator
        !
        ! n~_h (u_h, phi_i, phi_j)
        !
        ! =   ( u_h*grad Phi_j, u_h*grad Phi_i )_T
        !
        ! =   ( < (DU1) , (grad(Phi_j)_1) > , < (DU1) , (grad(Phi_i)_1) > )_T
        !         (DU2) , (grad(Phi_j)_2)       (DU2) , (grad(Phi_i)_2)
        !
        ! =   < (DU1) , (grad(Phi_j)_1) >  *  < (DU1) , (grad(Phi_j)_1) >
        !       (DU2) , (grad(Phi_j)_2)         (DU2) , (grad(Phi_j)_2)
        !
        ! =   dsumJ * dsumI
        !
        ! weighted with dlocalDelta, which is computed depending
        ! on the viscosity, velocity etc.
        
        select case (ndim)
        case (NDIM1D)

          do i=0,iend
          
            ! Get the matrix data      
            p_DlocalMatrix => RmatrixData(iy+i,ix+i)%p_Dentry

            if (bconstVel) then
            
              ! Loop over the elements in the current set.
              do iel = 1,nelements

                ! Loop over all cubature points on the current element
                do icubp = 1,npointsPerElement

                  ! Outer loop over the DOF's i=1..ndof on our current element,
                  ! which corresponds to the (test) basis functions Psi_i:
                  do idofe=1,p_rmatrixData%ndofTest

                    ! Fetch the contributions of the (test) basis functions Psi_i
                    ! into dbasIx/y
                    dbasIx = p_DbasTest(idofe,DER_DERIV1D_X,icubp,iel)
                    
                    ! Calculate dsumI, multiply with dscale, the local delta
                    ! and the cubature weight. Saves some multiplications.
                    dsumI = dscale * p_DcubWeight(icubp,iel) * p_DlocalDelta(iel) * &
                            (dvelX*dbasIx)

                    ! Inner loop over the DOF's j=1..ndof, which corresponds to
                    ! the (trial) basis function Phi_j:
                    do jdofe=1,p_rmatrixData%ndofTrial

                      ! Fetch the contributions of the (trial) basis function Phi_j
                      ! into dbasJx/y
                      dbasJx = p_DbasTrial(jdofe,DER_DERIV2D_X,icubp,iel)

                      ! Calculate dsumJ
                      dsumJ = dvelX*dbasJx

                      ! Multiply the values of the basis functions
                      ! (1st derivatives) by the cubature weight and sum up
                      ! into the local matrices.
                      p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                          dsumJ*dsumI

                    end do ! jdofe

                  end do ! idofe

                end do ! icubp

              end do ! iel
              
            else
            
              ! Loop over the elements in the current set.
              do iel = 1,nelements

                ! Loop over all cubature points on the current element
                do icubp = 1,npointsPerElement
                
                  ! Get the velocity in that point.
                  dvelX = p_Du(1,icubp,iel)

                  ! Outer loop over the DOF's i=1..ndof on our current element,
                  ! which corresponds to the (test) basis functions Psi_i:
                  do idofe=1,p_rmatrixData%ndofTest

                    ! Fetch the contributions of the (test) basis functions Psi_i
                    ! into dbasIx/y
                    dbasIx = p_DbasTest(idofe,DER_DERIV1D_X,icubp,iel)
                    
                    ! Calculate dsumI, multiply with dscale, the local delta
                    ! and the cubature weight. Saves some multiplications.
                    dsumI = dscale * p_DcubWeight(icubp,iel) * p_DlocalDelta(iel) * &
                            (dvelX*dbasIx)

                    ! Inner loop over the DOF's j=1..ndof, which corresponds to
                    ! the (trial) basis function Phi_j:
                    do jdofe=1,p_rmatrixData%ndofTrial

                      ! Fetch the contributions of the (trial) basis function Phi_j
                      ! into dbasJx/y
                      dbasJx = p_DbasTrial(jdofe,DER_DERIV2D_X,icubp,iel)

                      ! Calculate dsumJ
                      dsumJ = dvelX*dbasJx

                      ! Multiply the values of the basis functions
                      ! (1st derivatives) by the cubature weight and sum up
                      ! into the local matrices.
                      p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                          dsumJ*dsumI

                    end do ! jdofe

                  end do ! idofe

                end do ! icubp

              end do ! iel

            end if

          end do ! i

        case (NDIM2D)

          do i=0,iend
          
            ! Get the matrix data      
            p_DlocalMatrix => RmatrixData(iy+i,ix+i)%p_Dentry

            if (bconstVel) then
            
              ! Loop over the elements in the current set.
              do iel = 1,nelements

                ! Loop over all cubature points on the current element
                do icubp = 1,npointsPerElement

                  ! Outer loop over the DOF's i=1..ndof on our current element,
                  ! which corresponds to the (test) basis functions Psi_i:
                  do idofe=1,p_rmatrixData%ndofTest

                    ! Fetch the contributions of the (test) basis functions Psi_i
                    ! into dbasIx/y
                    dbasIx = p_DbasTest(idofe,DER_DERIV2D_X,icubp,iel)
                    dbasIy = p_DbasTest(idofe,DER_DERIV2D_Y,icubp,iel)
                    
                    ! Calculate dsumI, multiply with dscale, the local delta
                    ! and the cubature weight. Saves some multiplications.
                    dsumI = dscale * p_DcubWeight(icubp,iel) * p_DlocalDelta(iel) * &
                            (dvelX*dbasIx + dvelY*dbasIy)

                    ! Inner loop over the DOF's j=1..ndof, which corresponds to
                    ! the (trial) basis function Phi_j:
                    do jdofe=1,p_rmatrixData%ndofTrial

                      ! Fetch the contributions of the (trial) basis function Phi_j
                      ! into dbasJx/y
                      dbasJx = p_DbasTrial(jdofe,DER_DERIV2D_X,icubp,iel)
                      dbasJy = p_DbasTrial(jdofe,DER_DERIV2D_Y,icubp,iel)

                      ! Calculate dsumJ
                      dsumJ = dvelX*dbasJx + dvelY*dbasJy

                      ! Multiply the values of the basis functions
                      ! (1st derivatives) by the cubature weight and sum up
                      ! into the local matrices.
                      p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                          dsumJ*dsumI

                    end do ! jdofe

                  end do ! idofe

                end do ! icubp

              end do ! iel
              
            else
            
              ! Loop over the elements in the current set.
              do iel = 1,nelements

                ! Loop over all cubature points on the current element
                do icubp = 1,npointsPerElement
                
                  ! Get the velocity in that point.
                  dvelX = p_Du(1,icubp,iel)
                  dvelY = p_Du(2,icubp,iel)

                  ! Outer loop over the DOF's i=1..ndof on our current element,
                  ! which corresponds to the (test) basis functions Psi_i:
                  do idofe=1,p_rmatrixData%ndofTest

                    ! Fetch the contributions of the (test) basis functions Psi_i
                    ! into dbasIx/y
                    dbasIx = p_DbasTest(idofe,DER_DERIV2D_X,icubp,iel)
                    dbasIy = p_DbasTest(idofe,DER_DERIV2D_Y,icubp,iel)
                    
                    ! Calculate dsumI, multiply with dscale, the local delta
                    ! and the cubature weight. Saves some multiplications.
                    dsumI = dscale * p_DcubWeight(icubp,iel) * p_DlocalDelta(iel) * &
                            (dvelX*dbasIx + dvelY*dbasIy)

                    ! Inner loop over the DOF's j=1..ndof, which corresponds to
                    ! the (trial) basis function Phi_j:
                    do jdofe=1,p_rmatrixData%ndofTrial

                      ! Fetch the contributions of the (trial) basis function Phi_j
                      ! into dbasJx/y
                      dbasJx = p_DbasTrial(jdofe,DER_DERIV2D_X,icubp,iel)
                      dbasJy = p_DbasTrial(jdofe,DER_DERIV2D_Y,icubp,iel)

                      ! Calculate dsumJ
                      dsumJ = dvelX*dbasJx + dvelY*dbasJy

                      ! Multiply the values of the basis functions
                      ! (1st derivatives) by the cubature weight and sum up
                      ! into the local matrices.
                      p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                          dsumJ*dsumI

                    end do ! jdofe

                  end do ! idofe

                end do ! icubp

              end do ! iel

            end if

          end do ! i
          
        case (NDIM3D)

          do i=0,iend
          
            ! Get the matrix data      
            p_DlocalMatrix => RmatrixData(iy+i,ix+i)%p_Dentry

            if (bconstVel) then
            
              ! Loop over the elements in the current set.
              do iel = 1,nelements

                ! Loop over all cubature points on the current element
                do icubp = 1,npointsPerElement

                  ! Outer loop over the DOF's i=1..ndof on our current element,
                  ! which corresponds to the (test) basis functions Psi_i:
                  do idofe=1,p_rmatrixData%ndofTest

                    ! Fetch the contributions of the (test) basis functions Psi_i
                    ! into dbasIx/y
                    dbasIx = p_DbasTest(idofe,DER_DERIV3D_X,icubp,iel)
                    dbasIy = p_DbasTest(idofe,DER_DERIV3D_Y,icubp,iel)
                    dbasIz = p_DbasTest(idofe,DER_DERIV3D_Z,icubp,iel)
                    
                    ! Calculate dsumI, multiply with dscale, the local delta
                    ! and the cubature weight. Saves some multiplications.
                    dsumI = dscale * p_DcubWeight(icubp,iel) * p_DlocalDelta(iel) * &
                            (dvelX*dbasIx + dvelY*dbasIy + dvelZ*dbasIz)

                    ! Inner loop over the DOF's j=1..ndof, which corresponds to
                    ! the (trial) basis function Phi_j:
                    do jdofe=1,p_rmatrixData%ndofTrial

                      ! Fetch the contributions of the (trial) basis function Phi_j
                      ! into dbasJx/y
                      dbasJx = p_DbasTrial(jdofe,DER_DERIV3D_X,icubp,iel)
                      dbasJy = p_DbasTrial(jdofe,DER_DERIV3D_Y,icubp,iel)
                      dbasJz = p_DbasTrial(jdofe,DER_DERIV3D_Z,icubp,iel)

                      ! Calculate dsumJ
                      dsumJ = dvelX*dbasJx + dvelY*dbasJy + dvelZ*dbasJz

                      ! Multiply the values of the basis functions
                      ! (1st derivatives) by the cubature weight and sum up
                      ! into the local matrices.
                      p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                          dsumJ*dsumI

                    end do ! jdofe

                  end do ! idofe

                end do ! icubp

              end do ! iel
              
            else
            
              ! Loop over the elements in the current set.
              do iel = 1,nelements

                ! Loop over all cubature points on the current element
                do icubp = 1,npointsPerElement
                
                  ! Get the velocity in that point.
                  dvelX = p_Du(1,icubp,iel)
                  dvelY = p_Du(2,icubp,iel)
                  dvelZ = p_Du(3,icubp,iel)

                  ! Outer loop over the DOF's i=1..ndof on our current element,
                  ! which corresponds to the (test) basis functions Psi_i:
                  do idofe=1,p_rmatrixData%ndofTest

                    ! Fetch the contributions of the (test) basis functions Psi_i
                    ! into dbasIx/y
                    dbasIx = p_DbasTest(idofe,DER_DERIV3D_X,icubp,iel)
                    dbasIy = p_DbasTest(idofe,DER_DERIV3D_Y,icubp,iel)
                    dbasIz = p_DbasTest(idofe,DER_DERIV3D_Z,icubp,iel)
                    
                    ! Calculate dsumI, multiply with dscale, the local delta
                    ! and the cubature weight. Saves some multiplications.
                    dsumI = dscale * p_DcubWeight(icubp,iel) * p_DlocalDelta(iel) * &
                            (dvelX*dbasIx + dvelY*dbasIy + dvelZ*dbasIz)

                    ! Inner loop over the DOF's j=1..ndof, which corresponds to
                    ! the (trial) basis function Phi_j:
                    do jdofe=1,p_rmatrixData%ndofTrial

                      ! Fetch the contributions of the (trial) basis function Phi_j
                      ! into dbasJx/y
                      dbasJx = p_DbasTrial(jdofe,DER_DERIV3D_X,icubp,iel)
                      dbasJy = p_DbasTrial(jdofe,DER_DERIV3D_Y,icubp,iel)
                      dbasJz = p_DbasTrial(jdofe,DER_DERIV3D_Z,icubp,iel)

                      ! Calculate dsumJ
                      dsumJ = dvelX*dbasJx + dvelY*dbasJy + dvelZ*dbasJz

                      ! Multiply the values of the basis functions
                      ! (1st derivatives) by the cubature weight and sum up
                      ! into the local matrices.
                      p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                          dsumJ*dsumI

                    end do ! jdofe

                  end do ! idofe

                end do ! icubp

              end do ! iel

            end if

          end do ! i

        end select
        
      else
      
        ! -----------------------
        ! Vector-valued FE space.
        ! -----------------------
        ! This implementation works for all types of finite elements, but it
        ! is slightly slower for scalar FE spaces due to an additional inner loop.

      end if ! ndimfe = 1

    else

    end if

  end subroutine

!  !****************************************************************************
!
!!<subroutine>
!
!  subroutine bma_fcalc_streamlinediff(RmatrixData,rassemblyData,rmatrixAssembly,&
!      npointsPerElement,nelements,revalVectors,rcollection)
!
!!<description>  
!    ! Calculates the Streamline diffusion stabilisation.
!    !
!    ! If rcollection must be specified, the following parameters are expected:
!    !
!    ! rcollection%DquickAccess(1) = stabilisation parameter dupsam
!    ! rcollection%DquickAccess(2) = maximum norm of the velocity vector
!    ! rcollection%IquickAccess(1) = x-position in the destination matrix
!    ! rcollection%IquickAccess(2) = y-position in the destination matrix
!    ! rcollection%IquickAccess(3) = Type of SD method to apply.
!    !                               = 0: Use simple SD stabilisation:
!    !                                    ddelta = dupsam * h_T.
!    !                               = 1: Use Samarskji SD stabilisation; usually dupsam = 0.1 .. 2.
!    !                                    ddelta = dupsam * h_t/||u||_T * 2*Re_T/(1+Re_T)
!    ! rcollection%IquickAccess(4) = Method how to compute the local mesh width.
!    !                               = 0: Element volume based.
!    !                               = 1: Ray shooting method.
!    !
!    ! Furthermore, the following vectors must be provided in revalVectors:
!    !
!    !   revalVectors(1) = Element-based temporary vector
!    !   revalVectors(2) = Element-based temporary vector
!    !   revalVectors(1) = velocity field, n-dimensional, including 1st derivative
!    !   revalVectors(n+1) = Viscosity coefficient in all cubature points on all element
!!</description>
!
!!<inputoutput>
!    ! Matrix data of all matrices. The arrays p_Dentry of all submatrices
!    ! have to be filled with data.
!    type(t_bmaMatrixData), dimension(:,:), intent(inout), target :: RmatrixData
!!</inputoutput>
!
!!<input>
!    ! Data necessary for the assembly. Contains determinants and
!    ! cubature weights for the cubature,...
!    type(t_bmaMatrixAssemblyData), intent(in) :: rassemblyData
!
!    ! Structure with all data about the assembly
!    type(t_bmaMatrixAssembly), intent(in) :: rmatrixAssembly
!
!    ! Number of points per element
!    integer, intent(in) :: npointsPerElement
!
!    ! Number of elements
!    integer, intent(in) :: nelements
!
!    ! Values of FEM functions automatically evaluated in the
!    ! cubature points.
!    type(t_fev2Vectors), intent(in) :: revalVectors
!
!    ! User defined collection structure
!    type(t_collection), intent(inout), target, optional :: rcollection
!!</input>
!
!!</subroutine>
!
!    ! Local variables
!    real(DP) :: dbasI, dbasJ
!    integer :: iel, icubp, idofe, jdofe, ivar, nvar, ndim
!    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix11,p_DlocalMatrix22
!    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial,p_DbasTest
!    real(DP), dimension(:,:), pointer :: p_DcubWeight
!    type(t_bmaMatrixData), pointer :: p_rmatrixData
!    real(DP), dimension(:,:,:), pointer :: p_Dnu
!    real(DP), dimension(:), pointer :: p_DlocalDelta,p_DmeshWidth
!    real(DP), dimension(:,:,:), pointer :: p_Du
!
!    integer :: ix, iy, ndimfe, idimfe, cstabiltype, imeshwidthmethod
!    real(DP) :: dupsam, dumaxR
!
!    if (.not. present(rcollection)) then
!      call output_line ("Parameters missing.",&
!          OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_streamlinediff")
!      call sys_halt()
!    end if
!
!    ! Get parameters
!    dupsam = rcollection%DquickAccess(1)
!    dumaxR = 1.0_DP / rcollection%DquickAccess(2)
!    
!    ix = rcollection%IquickAccess(1)
!    iy = rcollection%IquickAccess(2)
!    cstabiltype = rcollection%IquickAccess(3)
!    imeshwidthmethod = rcollection%IquickAccess(4)
!
!    ! Cancel if nothing to do or parameters wrong
!    if (dupsam .eq. 0.0_DP) return
!
!    if ((ix .lt. 1) .or. (iy .lt. 1) .or. &
!        (ix .gt. ubound(RmatrixData,2)) .or. (iy .gt. ubound(RmatrixData,1))) then
!      call output_line ("Parameters wrong.",&
!          OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_streamlinediff")
!      call sys_halt()
!    end if
!  
!  end subroutine

end module
