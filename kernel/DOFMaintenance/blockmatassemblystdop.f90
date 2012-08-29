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
!# They can directly be used with bma_buildMatrix, bma_buildVector and 
!# bma_buildIntegral:
!#
!# 1.) bma_fcalc_mass
!#     -> Used with bma_buildMatrix, this calculates mass matrices in all
!#        diagonal blocks of a block matrix.
!#
!# 2.) bma_fcalc_laplace
!#     -> Used with bma_buildMatrix, this calculates Laplace matrices matrices
!#        in all diagonal blocks of a block matrix.
!#
!# 3.) bma_fcalc_convection
!#     -> Used with bma_buildMatrix, this routine can be used to assemble
!#        a standard convection operator in the first diagonal blocks
!#        of a block matrix
!#
!# 4.) bma_fcalc_rhsOne
!#     -> Calculates the RHS vector based on the function f=1.
!#
!# 5.) bma_fcalc_rhsBubble
!#     -> Calculates the RHS vector based on the function f=32*y*(1-y)+32*x*(1-x)
!#        which is the RHS for u=16*x*(1-x)*y*(1-y) in the Laplace
!#        equation -Laplace(u)=f.
!#
!# 6.) bma_fcalc_rhsBubblePlusFE
!#     -> Calculates the RHS vector based on the function 
!#        f=32*y*(1-y)+32*x*(1-x) + v(x,y)
!#        with v(x,y) being a finite element function passed via parameters.
!#
!# 6.) bma_fcalc_rhsFE
!#     -> Calculates the RHS vector based on the function 
!#        f=v(x,y)
!#        with v(x,y) being a finite element function passed via parameters.
!#
!# 8.) bma_fcalc_integralOne
!#     -> Calculates the integral of the function v=1 (which results in the
!#        size of the domain).
!#
!# 9.) bma_fcalc_integralFE
!#     -> Calculates the integral of an arbitrary FEM function.
!#
!# 10.) bma_fcalc_bubbleL2error
!#      -> Calculates the squared L2 error of a FEM function to a 
!#         bubble function
!#        
!# 11.) bma_fcalc_bubbleH1error
!#     -> Calculates the squared H1 error of a FEM function to a 
!#        bubble function
!#
!# 12.) bma_fcalc_L2norm
!#     -> Calculates the squared L2 norm of a FEM function
!#        
!# 13.) bma_fcalc_H1norm
!#     -> Calculates the squared H1 (semi-)norm of a FEM function
!# </purpose>
!##############################################################################

module blockmatassemblystdop

  use fsystem
  use genoutput
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

  public :: bma_fcalc_mass
  public :: bma_fcalc_laplace
  public :: bma_fcalc_convection
  
  public :: bma_fcalc_rhsOne
  public :: bma_fcalc_rhsBubble
  public :: bma_fcalc_rhsBubblePlusFE
  public :: bma_fcalc_rhsFE
  
  public :: bma_fcalc_integralOne
  public :: bma_fcalc_integralFE
  public :: bma_fcalc_bubbleL2error
  public :: bma_fcalc_bubbleH1error
  public :: bma_fcalc_L2norm
  public :: bma_fcalc_H1norm

contains

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_mass(RmatrixData,rassemblyData,rmatrixAssembly,&
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

!<subroutine>

    ! Local variables
    real(DP) :: dbasI, dbasJ
    integer :: iel, icubp, idofe, jdofe, i, ivar, nvar
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix
    real(DP), dimension(:,:,:,:), pointer :: p_DlocalMatrixIntl
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial,p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaMatrixData), pointer :: p_rmatrixData

    integer :: istart, iend
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

    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Loop over the diagonal blocks
    do i = istart,iend

      ! Get local data
      p_rmatrixData => RmatrixData(i,i)
      p_DbasTrial => RmatrixData(i,i)%p_DbasTrial
      p_DbasTest => RmatrixData(i,i)%p_DbasTest

      ! Interleaved matrix?
      if (.not. p_rmatrixData%bisInterleaved) then

        ! Get the matrix data      
        p_DlocalMatrix => RmatrixData(i,i)%p_Dentry

        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Phi_i:
            do idofe=1,p_rmatrixData%ndofTest

              ! Fetch the contributions of the (test) basis functions Phi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the basis function Phi_j:
              do jdofe=1,p_rmatrixData%ndofTrial

                ! Fetch the contributions of the (trial) basis function Phi_j
                ! into dbasJ
                dbasJ = p_DbasTrial(jdofe,DER_FUNC,icubp,iel)

                ! Multiply the values of the basis functions
                ! (1st derivatives) by the cubature weight and sum up
                ! into the local matrices.
                p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                    dscale * p_DcubWeight(icubp,iel) * dbasJ*dbasI

              end do ! idofe

            end do ! jdofe

          end do ! icubp

        end do ! iel

      else

        ! Get the matrix data      
        p_DlocalMatrixIntl => RmatrixData(i,i)%p_DentryIntl

        ! Interleave-data
        nvar = RmatrixData(i,i)%nvar

        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Phi_i:
            do idofe=1,p_rmatrixData%ndofTest

              ! Fetch the contributions of the (test) basis functions Phi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the basis function Phi_j:
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

              end do ! idofe

            end do ! jdofe

          end do ! icubp

        end do ! iel

      end if

    end do

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_laplace(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the Laplace operator in all diagonal matrices
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

!<subroutine>

    real(DP) :: dbasIx, dbasJx, dbasIy, dbasJy, dbasIz, dbasJz
    integer :: iel, icubp, idofe, jdofe, i, ivar, nvar
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix
    real(DP), dimension(:,:,:,:), pointer :: p_DlocalMatrixIntl
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial,p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaMatrixData), pointer :: p_rmatrixData

    integer :: istart, iend
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

    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Loop through all diagonal blocks   
    do i = istart, iend

      ! Get local data
      p_rmatrixData => RmatrixData(i,i)
      p_DbasTrial => RmatrixData(i,i)%p_DbasTrial
      p_DbasTest => RmatrixData(i,i)%p_DbasTest

      ! Interleaved matrix?
      if (.not. p_rmatrixData%bisInterleaved) then

        ! Get the matrix data
        p_DlocalMatrix => RmatrixData(i,i)%p_Dentry

        ! What's the dimension?
        select case (rmatrixAssembly%p_rtriangulation%ndim)

        case (NDIM1D)

          ! 1D Laplace matrix

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Phi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Phi_i
                ! into dbasI
                dbasIx = p_DbasTest(idofe,DER_DERIV1D_X,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the basis function Phi_j:
                do jdofe=1,p_rmatrixData%ndofTrial

                  ! Fetch the contributions of the (trial) basis function Phi_j
                  ! into dbasJ
                  dbasJx = p_DbasTrial(jdofe,DER_DERIV1D_X,icubp,iel)

                  ! Multiply the values of the basis functions
                  ! (1st derivatives) by the cubature weight and sum up
                  ! into the local matrices.
                  p_DlocalMatrix(jdofe,idofe,iel) = p_DlocalMatrix(jdofe,idofe,iel) + &
                      dscale * p_DcubWeight(icubp,iel) * dbasJx*dbasIx

                end do ! idofe

              end do ! jdofe

            end do ! icubp

          end do ! iel

        case (NDIM2D)

          ! 2D Laplace matrix

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Phi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Phi_i
                ! into dbasI
                dbasIx = p_DbasTest(idofe,DER_DERIV2D_X,icubp,iel)
                dbasIy = p_DbasTest(idofe,DER_DERIV2D_Y,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the basis function Phi_j:
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

                end do ! idofe

              end do ! jdofe

            end do ! icubp

          end do ! iel

        case (NDIM3D)

          ! 3D Laplace matrix

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Phi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Phi_i
                ! into dbasI
                dbasIx = p_DbasTest(idofe,DER_DERIV3D_X,icubp,iel)
                dbasIy = p_DbasTest(idofe,DER_DERIV3D_Y,icubp,iel)
                dbasIz = p_DbasTest(idofe,DER_DERIV3D_Z,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the basis function Phi_j:
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

                end do ! idofe

              end do ! jdofe

            end do ! icubp

          end do ! iel

        end select

      else

        ! Get the matrix data      
        p_DlocalMatrixIntl => RmatrixData(i,i)%p_DentryIntl

        ! Interleave-data
        nvar = RmatrixData(i,i)%nvar

        ! What's the dimension?
        select case (rmatrixAssembly%p_rtriangulation%ndim)

        case (NDIM1D)

          ! 1D Laplace matrix

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Phi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Phi_i
                ! into dbasI
                dbasIx = p_DbasTest(idofe,DER_DERIV1D_X,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the basis function Phi_j:
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

                end do ! idofe

              end do ! jdofe

            end do ! icubp

          end do ! iel

        case (NDIM2D)

          ! 2D Laplace matrix

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Phi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Phi_i
                ! into dbasI
                dbasIx = p_DbasTest(idofe,DER_DERIV2D_X,icubp,iel)
                dbasIy = p_DbasTest(idofe,DER_DERIV2D_Y,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the basis function Phi_j:
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

                end do ! idofe

              end do ! jdofe

            end do ! icubp

          end do ! iel

        case (NDIM3D)

          ! 3D Laplace matrix

          ! Loop over the elements in the current set.
          do iel = 1,nelements

            ! Loop over all cubature points on the current element
            do icubp = 1,npointsPerElement

              ! Outer loop over the DOF's i=1..ndof on our current element,
              ! which corresponds to the (test) basis functions Phi_i:
              do idofe=1,p_rmatrixData%ndofTest

                ! Fetch the contributions of the (test) basis functions Phi_i
                ! into dbasI
                dbasIx = p_DbasTest(idofe,DER_DERIV3D_X,icubp,iel)
                dbasIy = p_DbasTest(idofe,DER_DERIV3D_Y,icubp,iel)
                dbasIz = p_DbasTest(idofe,DER_DERIV3D_Z,icubp,iel)

                ! Inner loop over the DOF's j=1..ndof, which corresponds to
                ! the basis function Phi_j:
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

                end do ! idofe

              end do ! jdofe

            end do ! icubp

          end do ! iel

        end select

      end if

    end do

  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_convection(RmatrixData,rassemblyData,rmatrixAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates a convection operator "(u grad) ." into the top/left position
    ! of a block matrix (position 1/1).
    !
    ! Note: If rcollection is not specified, the matrix is calculated
    ! in all diagonal blocks with a multiplier of 1.
    ! If rcollection is specified, the following parameters are expected:
    !
    ! rcollection%DquickAccess(1) = multiplier in front of the operator.
    ! rcollection%IquickAccess(1) = 0, if the convection is a constant vector field.
    !                                  in this case:
    !                                  1D: rcollection%DquickAccess(2)   = x-velocity
    !                                  2D: rcollection%DquickAccess(2:3) = x/y-velocity
    !                                  3D: rcollection%DquickAccess(2:4) = x/y/z-velocity
    !                             = 1, if the convection is specified by a
    !                                  finite element velocity field. In this case,
    !                                  a finite element velocity field must be specified
    !                                  as parameter revalVectors to the call of 
    !                                  bma_buildMatrix. The first vector must be the
    !                                  X-velocity, the 2nd the Y-velocity and 
    !                                  the third the Z-velocity.
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
    !     rcollection%IquickAccess(1) = 1          ! Nonconstant viscosity
    !     rcollection%DquickAccess(1) = 1.0_DP     ! Scaling
    !
    !     ! Add the X-, Y- and Z-velocity to revalVectors
    !     call fev2_addVectorToEvalList(revalVectors,rvelocity%RvectorBlock(1),0)
    !     call fev2_addVectorToEvalList(revalVectors,rvelocity%RvectorBlock(2),0)
    !     call fev2_addVectorToEvalList(revalVectors,rvelocity%RvectorBlock(3),0)
    !
    !     ! Set up the matrix
    !     call bma_buildMatrix (rmatrix,BMA_CALC_STANDARD, bma_fcalc_convection, &
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
    ! The routine <verb>fev2_addVectorToEvalList</verb> allows to define
    ! the evaluation of derivatives as well. In 3D, e.g., one may apply
    !
    ! <verb>
    !     call fev2_addVectorToEvalList(revalVectors,rvelocity%RvectorBlock(1),1)   ! u1
    !     call fev2_addVectorToEvalList(revalVectors,rvelocity%RvectorBlock(2),1)   ! u2
    !     call fev2_addVectorToEvalList(revalVectors,rvelocity%RvectorBlock(3),1)   ! u3
    ! </verb>
    !
    ! which calculates function values as well as 1st derivatives of the complete
    ! vector field (due to the "1" at the end). The calculated values in the
    ! cubature points can then be found in the "p_Ddata" elements of revalVectors:
    !
    ! <verb>
    !   revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC)      = u1
    !   revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_DERIV3D_X) = d/dx u1
    !   revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_DERIV3D_Y) = d/dy u1
    !   revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_DERIV3D_Z) = d/dz u1
    !
    !   revalVectors%p_RvectorData(2)%p_Ddata(:,:,DER_FUNC)      = u2
    !   revalVectors%p_RvectorData(2)%p_Ddata(:,:,DER_DERIV3D_X) = d/dx u2
    !   revalVectors%p_RvectorData(2)%p_Ddata(:,:,DER_DERIV3D_Y) = d/dy u2
    !   revalVectors%p_RvectorData(2)%p_Ddata(:,:,DER_DERIV3D_Z) = d/dz u2
    !
    !   revalVectors%p_RvectorData(3)%p_Ddata(:,:,DER_FUNC)      = u3
    !   revalVectors%p_RvectorData(3)%p_Ddata(:,:,DER_DERIV3D_X) = d/dx u3
    !   revalVectors%p_RvectorData(3)%p_Ddata(:,:,DER_DERIV3D_Y) = d/dy u3
    !   revalVectors%p_RvectorData(3)%p_Ddata(:,:,DER_DERIV3D_Z) = d/dz u3
    ! </verb>
    !
    ! in all cubature points on all elements. The vector data in
    ! revalVectors%p_RvectorData appears exactly in the order, the vectors
    ! are added to revalVectors by fev2_addVectorToEvalList.
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

!<subroutine>

    ! Local variables
    real(DP) :: dbasI, dbasJx, dbasJy, dbasJz
    integer :: iel, icubp, idofe, jdofe
    real(DP), dimension(:,:,:), pointer :: p_DlocalMatrix11,p_DlocalMatrix22,p_DlocalMatrix33
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTrial,p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaMatrixData), pointer :: p_rmatrixData11,p_rmatrixData22,p_rmatrixData33
    real(DP), dimension(:,:,:), pointer :: p_Du1,p_Du2,p_Du3

    integer :: ndim
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
    
    if (present(rcollection)) then
      dscale = rcollection%DquickAccess(1)
      
      ! Constant velocity?
      bvelConst = (rcollection%IquickAccess(1) .eq. 0)
      
      if (bvelConst) then
        ! Get the constant velocity
        dvelX = rcollection%DquickAccess(2)
        if (ndim .ge. 2) dvelY = rcollection%DquickAccess(3)
        if (ndim .ge. 3) dvelZ = rcollection%DquickAccess(4)
      end if
      
    end if
    
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Get local data
    p_DbasTrial => RmatrixData(1,1)%p_DbasTrial
    p_DbasTest => RmatrixData(1,1)%p_DbasTest

    ! Set up the local matrix of the convection.
    select case (ndim)
    case (NDIM1D)
      ! Matrices to be set up
      p_rmatrixData11 => RmatrixData(1,1)
      
      p_DlocalMatrix11 => RmatrixData(1,1)%p_Dentry
      
      ! Currently, interleaved matrices are not supported
      if (p_rmatrixData11%bisInterleaved) then
        call output_line ("Interleaved matrices not supported",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_convection")
        call sys_halt()
      end if

      if (bvelConst) then
      
        ! Set up the matrix for constant velocity.
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Phi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Phi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the basis function Phi_j:
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

              end do ! idofe

            end do ! jdofe

          end do ! icubp

        end do ! iel
        
      else

        if (revalVectors%ncount .lt. 1) then
          call output_line ("FEM function missing.",&
              OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_convection")
          call sys_halt()
        end if

        ! Set up the matrix for nonconstant velocity.
        !
        ! Get the velocity field from the parameters
        p_Du1 => revalVectors%p_RvectorData(1)%p_Ddata
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Velocity field in this cubature point
            dvelX = p_Du1(icubp,iel,DER_FUNC)
            
            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Phi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Phi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the basis function Phi_j:
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

              end do ! idofe

            end do ! jdofe

          end do ! icubp

        end do ! iel
        
      end if
      
    case (NDIM2D)

      ! Matrices to be set up
      p_rmatrixData11 => RmatrixData(1,1)
      p_rmatrixData22 => RmatrixData(2,2)
      
      p_DlocalMatrix11 => RmatrixData(1,1)%p_Dentry
      p_DlocalMatrix22 => RmatrixData(2,2)%p_Dentry

      ! Currently, interleaved matrices are not supported
      if (p_rmatrixData11%bisInterleaved) then
        call output_line ("Interleaved matrices not supported",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_convection")
        call sys_halt()
      end if

      if (bvelConst) then
      
        ! Set up the matrix for constant velocity.
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Phi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Phi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the basis function Phi_j:
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

              end do ! idofe

            end do ! jdofe

          end do ! icubp

        end do ! iel
        
      else

        if (revalVectors%ncount .lt. 2) then
          call output_line ("FEM function missing.",&
              OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_convection")
          call sys_halt()
        end if

        ! Set up the matrix for nonconstant velocity.
        !
        ! Get the velocity field from the parameters
        p_Du1 => revalVectors%p_RvectorData(1)%p_Ddata
        p_Du2 => revalVectors%p_RvectorData(2)%p_Ddata
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Velocity field in this cubature point
            dvelX = p_Du1(icubp,iel,DER_FUNC)
            dvelY = p_Du2(icubp,iel,DER_FUNC)
            
            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Phi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Phi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the basis function Phi_j:
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

              end do ! idofe

            end do ! jdofe

          end do ! icubp

        end do ! iel
        
      end if

    case (NDIM3D)

      ! Matrices to be set up
      p_rmatrixData11 => RmatrixData(1,1)
      p_rmatrixData22 => RmatrixData(2,2)
      p_rmatrixData33 => RmatrixData(3,3)

      p_DlocalMatrix11 => RmatrixData(1,1)%p_Dentry
      p_DlocalMatrix22 => RmatrixData(2,2)%p_Dentry
      p_DlocalMatrix33 => RmatrixData(3,3)%p_Dentry

      ! Currently, interleaved matrices are not supported
      if (p_rmatrixData11%bisInterleaved) then
        call output_line ("Interleaved matrices not supported",&
            OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_convection")
        call sys_halt()
      end if

      if (bvelConst) then
      
        ! Set up the matrix for constant velocity.
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Phi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Phi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the basis function Phi_j:
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

              end do ! idofe

            end do ! jdofe

          end do ! icubp

        end do ! iel
        
      else

        if (revalVectors%ncount .lt. 3) then
          call output_line ("FEM function missing.",&
              OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_bubbleH1error")
          call sys_halt()
        end if

        ! Set up the matrix for nonconstant velocity.
        !
        ! Get the velocity field from the parameters
        p_Du1 => revalVectors%p_RvectorData(1)%p_Ddata
        p_Du2 => revalVectors%p_RvectorData(2)%p_Ddata
        p_Du3 => revalVectors%p_RvectorData(3)%p_Ddata
      
        ! Loop over the elements in the current set.
        do iel = 1,nelements

          ! Loop over all cubature points on the current element
          do icubp = 1,npointsPerElement

            ! Velocity field in this cubature point
            dvelX = p_Du1(icubp,iel,DER_FUNC)
            dvelY = p_Du2(icubp,iel,DER_FUNC)
            dvelZ = p_Du3(icubp,iel,DER_FUNC)
            
            ! Outer loop over the DOF's i=1..ndof on our current element,
            ! which corresponds to the (test) basis functions Phi_i:
            do idofe=1,p_rmatrixData11%ndofTest

              ! Fetch the contributions of the (test) basis functions Phi_i
              ! into dbasI
              dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)

              ! Inner loop over the DOF's j=1..ndof, which corresponds to
              ! the basis function Phi_j:
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

              end do ! idofe

            end do ! jdofe

          end do ! icubp

        end do ! iel
        
      end if

    end select

  end subroutine

  !****************************************************************************
  !****************************************************************************
  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_rhsOne(rvectorData,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates a right-hand side vector according to the right-hand
    ! side function f=1 in all components.
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
    
!<subroutine>

    ! Local variables
    real(DP) :: dbasI, dval
    integer :: icomp
    integer :: iel, icubp, idofe
    real(DP), dimension(:,:), pointer :: p_DlocalVector
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaVectorData), pointer :: p_rvectorData
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    
    ! Loop over the components
    do icomp = 1,size(rvectorData)

      ! Get the data arrays of the subvector
      p_rvectorData => RvectorData(icomp)
      p_DlocalVector => RvectorData(icomp)%p_Dentry
      p_DbasTest => RvectorData(icomp)%p_DbasTest
    
      ! Loop over the elements in the current set.
      do iel = 1,nelements

        ! Loop over all cubature points on the current element
        do icubp = 1,npointsPerElement

          ! The value of the RHS is just the 1.0-function here
          dval = 1.0_DP
          
          ! Outer loop over the DOF's i=1..ndof on our current element,
          ! which corresponds to the (test) basis functions Phi_i:
          do idofe=1,p_rvectorData%ndofTest
          
            ! Fetch the contributions of the (test) basis functions Phi_i
            ! into dbasI
            dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)
            
            ! Multiply the values of the basis functions
            ! (1st derivatives) by the cubature weight and sum up
            ! into the local vectors.
            p_DlocalVector(idofe,iel) = p_DlocalVector(idofe,iel) + &
                p_DcubWeight(icubp,iel) * dval * dbasI
            
          end do ! jdofe

        end do ! icubp
      
      end do ! iel
      
    end do ! icomp
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_rhsBubble(rvectorData,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates a right-hand side vector according to the right-hand
    ! side function f=32*y*(1-y)+32*x*(1-x).
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
    
!<subroutine>

    ! Local variables
    real(DP) :: dbasI, dval, dx, dy
    integer :: icomp
    integer :: iel, icubp, idofe
    real(DP), dimension(:,:), pointer :: p_DlocalVector
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaVectorData), pointer :: p_rvectorData
    real(DP), dimension(:,:,:), pointer :: p_Dpoints
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal
    
    ! Loop over the components
    do icomp = 1,size(rvectorData)

      ! Get the data arrays of the subvector
      p_rvectorData => RvectorData(icomp)
      p_DlocalVector => RvectorData(icomp)%p_Dentry
      p_DbasTest => RvectorData(icomp)%p_DbasTest
    
      ! Loop over the elements in the current set.
      do iel = 1,nelements

        ! Loop over all cubature points on the current element
        do icubp = 1,npointsPerElement
        
          ! Get the coordinates of the cubature point.
          dx = p_Dpoints(1,icubp,iel)
          dy = p_Dpoints(2,icubp,iel)

          ! Calculate the values of the RHS using the coordinates
          ! of the cubature points.
          dval = 32.0_DP*dy*(1.0_DP-dy) + 32_DP*dx*(1.0_DP-dx)
          
          ! Outer loop over the DOF's i=1..ndof on our current element,
          ! which corresponds to the (test) basis functions Phi_i:
          do idofe=1,p_rvectorData%ndofTest
          
            ! Fetch the contributions of the (test) basis functions Phi_i
            ! into dbasI
            dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)
            
            ! Multiply the values of the basis functions
            ! (1st derivatives) by the cubature weight and sum up
            ! into the local vectors.
            p_DlocalVector(idofe,iel) = p_DlocalVector(idofe,iel) + &
                p_DcubWeight(icubp,iel) * dval * dbasI
            
          end do ! jdofe

        end do ! icubp
      
      end do ! iel
      
    end do ! icomp
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_rhsBubblePlusFE(rvectorData,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates a right-hand side vector according to the right-hand
    ! side function f = 32*y*(1-y)+32*x*(1-x) + v(x,y)
    ! with v being a FEM function.
    !
    ! The FEM function v(x,y) must be provided in revalVectors.
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
    
    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>
    
!<subroutine>

    ! Local variables
    real(DP) :: dbasI, dval, dx, dy
    integer :: icomp
    integer :: iel, icubp, idofe
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
    
    if (revalVectors%ncount .eq. 0) then
      call output_line ("FEM function missing.",&
          OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_rhsBubblePlusFE")
      call sys_halt()
    end if
    
    ! Get the data array with the values of the FEM function
    ! in the cubature points
    p_Dfunc => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC2D)
    
    ! Loop over the components
    do icomp = 1,size(rvectorData)

      ! Get the data arrays of the subvector
      p_rvectorData => RvectorData(icomp)
      p_DlocalVector => RvectorData(icomp)%p_Dentry
      p_DbasTest => RvectorData(icomp)%p_DbasTest
    
      ! Loop over the elements in the current set.
      do iel = 1,nelements

        ! Loop over all cubature points on the current element
        do icubp = 1,npointsPerElement
        
          ! Get the coordinates of the cubature point.
          dx = p_Dpoints(1,icubp,iel)
          dy = p_Dpoints(2,icubp,iel)

          ! Calculate the values of the RHS in the cubature point:
          !     f = 32*y*(1-y)+32*x*(1-x) + v(x,y)
          dval = 32.0_DP*dy*(1.0_DP-dy) + 32_DP*dx*(1.0_DP-dx) + p_Dfunc(icubp,iel)
          
          ! Outer loop over the DOF's i=1..ndof on our current element,
          ! which corresponds to the (test) basis functions Phi_i:
          do idofe=1,p_rvectorData%ndofTest
          
            ! Fetch the contributions of the (test) basis functions Phi_i
            ! into dbasI
            dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)
            
            ! Multiply the values of the basis functions
            ! (1st derivatives) by the cubature weight and sum up
            ! into the local vectors.
            p_DlocalVector(idofe,iel) = p_DlocalVector(idofe,iel) + &
                p_DcubWeight(icubp,iel) * dval * dbasI
            
          end do ! jdofe

        end do ! icubp
      
      end do ! iel
      
    end do ! icomp
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_rhsFE(rvectorData,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

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
    
    ! Values of FEM functions automatically evaluated in the
    ! cubature points.
    type(t_fev2Vectors), intent(in) :: revalVectors

    ! User defined collection structure
    type(t_collection), intent(inout), target, optional :: rcollection
!</input>
    
!<subroutine>

    ! Local variables
    real(DP) :: dbasI, dval
    integer :: icomp
    integer :: iel, icubp, idofe
    real(DP), dimension(:,:), pointer :: p_DlocalVector
    real(DP), dimension(:,:,:,:), pointer :: p_DbasTest
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    type(t_bmaVectorData), pointer :: p_rvectorData
    real(DP), dimension(:,:), pointer :: p_Dfunc
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight

    if (revalVectors%ncount .ne. size(rvectorData)) then
      call output_line ("FEM function missing or wrong length.",&
          OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_rhsFE")
      call sys_halt()
    end if
    
    ! Loop over the components
    do icomp = 1,size(rvectorData)

      ! Get the data arrays of the subvector
      p_rvectorData => RvectorData(icomp)
      p_DlocalVector => RvectorData(icomp)%p_Dentry
      p_DbasTest => RvectorData(icomp)%p_DbasTest
    
      ! Get the data array with the values of the FEM function v_i
      ! in the cubature points.
      p_Dfunc => revalVectors%p_RvectorData(icomp)%p_Ddata(:,:,DER_FUNC2D)
      
      ! Loop over the elements in the current set.
      do iel = 1,nelements

        ! Loop over all cubature points on the current element
        do icubp = 1,npointsPerElement
        
          ! Calculate the values of the RHS in the cubature point:
          !     f_i = v_i(x,y)
          dval = p_Dfunc(icubp,iel)
          
          ! Outer loop over the DOF's i=1..ndof on our current element,
          ! which corresponds to the (test) basis functions Phi_i:
          do idofe=1,p_rvectorData%ndofTest
          
            ! Fetch the contributions of the (test) basis functions Phi_i
            ! into dbasI
            dbasI = p_DbasTest(idofe,DER_FUNC,icubp,iel)
            
            ! Multiply the values of the basis functions
            ! (1st derivatives) by the cubature weight and sum up
            ! into the local vectors.
            p_DlocalVector(idofe,iel) = p_DlocalVector(idofe,iel) + &
                p_DcubWeight(icubp,iel) * dval * dbasI
            
          end do ! jdofe

        end do ! icubp
      
      end do ! iel
      
    end do ! icomp
    
  end subroutine

  !****************************************************************************
  !****************************************************************************
  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_integralOne(dintvalue,rassemblyData,rvectorAssembly,&
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
    type(t_bmaIntegralAssembly), intent(in) :: rvectorAssembly

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

!<subroutine>

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

  subroutine bma_fcalc_integralFE(dintvalue,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the integral of an arbitrary finite element function v.
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
    type(t_bmaIntegralAssembly), intent(in) :: rvectorAssembly

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

!<subroutine>

    ! Local variables
    real(DP) :: dval
    integer :: iel, icubp, ivec
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:), pointer :: p_Dfunc
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    
    dintvalue = 0.0_DP

    ! Loop through all provided FEM functions
    do ivec = 1,revalVectors%ncount
    
      ! Skip interleaved vectors.
      if (revalVectors%p_RvectorData(ivec)%bisInterleaved) cycle

      ! Get the data array with the values of the FEM function
      ! in the cubature points
      p_Dfunc => revalVectors%p_RvectorData(ivec)%p_Ddata(:,:,DER_FUNC2D)

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

  subroutine bma_fcalc_L2norm(dintvalue,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the squared L2-norm of an arbitrary finite element 
    ! function: <tex> ||v||^2 </tex>.
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
    type(t_bmaIntegralAssembly), intent(in) :: rvectorAssembly

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

!<subroutine>

    ! Local variables
    real(DP) :: dval
    integer :: iel, icubp, ivec
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:), pointer :: p_Dfunc
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    
    dintvalue = 0.0_DP

    ! Loop through all provided FEM functions
    do ivec = 1,revalVectors%ncount
    
      ! Skip interleaved vectors.
      if (revalVectors%p_RvectorData(ivec)%bisInterleaved) cycle

      ! Get the data array with the values of the FEM function
      ! in the cubature points
      p_Dfunc => revalVectors%p_RvectorData(ivec)%p_Ddata(:,:,DER_FUNC)

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
      
    end do ! ivec
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_H1norm(dintvalue,rassemblyData,rvectorAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the squared H1 (semi-)norm of an arbitrary finite element 
    ! function <tex> ||v||^2_{H^1} </tex>.
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
    type(t_bmaIntegralAssembly), intent(in) :: rvectorAssembly

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

!<subroutine>

    ! Local variables
    real(DP) :: dderivX,dderivY,dderivZ
    integer :: iel, icubp, ivec,ndim
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:), pointer :: p_DderivX,p_DderivY,p_DderivZ
  
    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    
    dintvalue = 0.0_DP

    ! Loop through all provided FEM functions
    do ivec = 1,revalVectors%ncount
    
      ! Skip interleaved vectors.
      if (revalVectors%p_RvectorData(ivec)%bisInterleaved) cycle

      ! Dimension?
      ndim = revalVectors%p_RvectorData(ivec)%p_rvector%p_rspatialDiscr%p_rtriangulation%ndim
      select case (ndim)
      
      ! ----------------
      ! 1D
      ! ----------------
      case (NDIM1D)
      
        ! Get the data array with the values of the FEM function
        ! in the cubature points
        p_DderivX => revalVectors%p_RvectorData(ivec)%p_Ddata(:,:,DER_DERIV1D_X)

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
        p_DderivX => revalVectors%p_RvectorData(ivec)%p_Ddata(:,:,DER_DERIV2D_X)
        p_DderivY => revalVectors%p_RvectorData(ivec)%p_Ddata(:,:,DER_DERIV2D_Y)

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
        p_DderivX => revalVectors%p_RvectorData(ivec)%p_Ddata(:,:,DER_DERIV3D_X)
        p_DderivY => revalVectors%p_RvectorData(ivec)%p_Ddata(:,:,DER_DERIV3D_Y)
        p_DderivZ => revalVectors%p_RvectorData(ivec)%p_Ddata(:,:,DER_DERIV3D_Z)

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
      
    end do ! ivec
    
  end subroutine

  !****************************************************************************

!<subroutine>

  subroutine bma_fcalc_bubbleL2error(dintvalue,rassemblyData,rintegralAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the (squared) L2 error of an arbitrary FEM function v
    ! to the bubble function u=16x(1-x)y(1-y):
    !
    !  <tex>  || v - u ||^2_{L2}  </tex>
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

!<subroutine>

    ! Local variables
    real(DP) :: dval1,dval2,dx,dy
    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:), pointer :: p_Dfunc
    real(DP), dimension(:,:,:), pointer :: p_Dpoints
  
    ! Cancel if no FEM function is given.
    if (revalVectors%ncount .eq. 0) then
      call output_line ("FEM function missing.",&
          OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_bubbleL2error")
      call sys_halt()
    end if

    ! Skip interleaved vectors.
    if (revalVectors%p_RvectorData(1)%bisInterleaved) return

    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    
    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal
    
    ! Get the data array with the values of the FEM function
    ! in the cubature points
    p_Dfunc => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_FUNC2D)

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

  subroutine bma_fcalc_bubbleH1error(dintvalue,rassemblyData,rintegralAssembly,&
      npointsPerElement,nelements,revalVectors,rcollection)

!<description>  
    ! Calculates the (squared) H1 error of an arbitrary FEM function v
    ! to the bubble function u=16x(1-x)y(1-y)
    ! (based on the H1 semi-norm).
    !
    !  <tex>  | v - u |^2_{H1}  </tex>
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

!<subroutine>

    ! Local variables
    real(DP) :: dderivX1,dderivY1,dderivX2,dderivY2,dx,dy
    integer :: iel, icubp
    real(DP), dimension(:,:), pointer :: p_DcubWeight
    real(DP), dimension(:,:), pointer :: p_DderivX,p_DderivY
    real(DP), dimension(:,:,:), pointer :: p_Dpoints
  
    ! Calcel if no FEM function is given.
    if (revalVectors%ncount .eq. 0) then
      call output_line ("FEM function missing.",&
          OU_CLASS_ERROR,OU_MODE_STD,"bma_fcalc_bubbleH1error")
      call sys_halt()
    end if

    ! Skip interleaved vectors.
    if (revalVectors%p_RvectorData(1)%bisInterleaved) return

    ! Get cubature weights data
    p_DcubWeight => rassemblyData%p_DcubWeight
    
    ! Get the coordinates of the cubature points
    p_Dpoints => rassemblyData%revalElementSet%p_DpointsReal
    
    ! Get the data array with the values of the FEM function
    ! in the cubature points
    p_DderivX => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_DERIV2D_X)
    p_DderivY => revalVectors%p_RvectorData(1)%p_Ddata(:,:,DER_DERIV2D_Y)

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

end module
