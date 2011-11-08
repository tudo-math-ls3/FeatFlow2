!##############################################################################
!# ****************************************************************************
!# <name> errest </name>
!# ****************************************************************************
!#
!# <purpose>
!#    This application is used to test error estimation techniques.
!# </purpose>
!##############################################################################

program errest
  
  use element
  use fsystem
  use genoutput
  use linearsystemscalar
  use linearsystemblock
  use pprocerror
  use pprocgradients
  use pprocindicator
  use spatialdiscretisation
  use storage
  use triangulation
  use ucd

  implicit none

  type(t_triangulation) :: rtriangulation
  type(t_ucdExport) :: rexport, rexportOut
  type(t_vectorScalar) :: rdensity, rpressure, rvelocity_x, rvelocity_y, rvelocity, radvect, rindicator
  type(t_blockDiscretisation) :: rdiscretisation

  real(DP), dimension(:), pointer :: p_Ddensity, p_Dpressure, p_Dvelocity_x, p_Dvelocity_y, p_Dvelocity, p_Dadvect, p_Dindicator
  integer :: nlength

  ! The very first thing in every application:
  ! Initialise system-wide settings:
  call system_init()
  
  ! The very second thing in every program:
  ! Initialise the FEAT 2.0 storage management:
  call storage_init(999, 100)


  ! Read GMV file
  call ucd_readGMV ('../flagship/out/zpinch_energy_AMR4_tria.00008.gmv',&
                    rexport, rtriangulation)
  call ucd_infoVariables(rexport)

  ! Generate discretisation
  call spdiscr_initBlockDiscr(rdiscretisation, 1, rtriangulation)
  call spdiscr_initDiscr_triquad(rdiscretisation%RspatialDiscr(1), &
                                 EL_E001, EL_E011, SPDISC_CUB_AUTOMATIC,&
                                 SPDISC_CUB_AUTOMATIC, rtriangulation)
  
  ! Create scalar vector
  call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1), rdensity,  .true.)
  call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1), rpressure, .true.)
  call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1), rvelocity, .true.)
  call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1), rvelocity_x, .true.)
  call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1), rvelocity_y, .true.)
  call lsyssc_createVecByDiscr(rdiscretisation%RspatialDiscr(1), radvect, .true.)
  call lsyssc_createVector(rindicator, rtriangulation%NEL, .true.)
  call lsyssc_getbase_double(rdensity,  p_Ddensity)
  call lsyssc_getbase_double(rpressure, p_Dpressure)
  call lsyssc_getbase_double(rvelocity, p_Dvelocity)
  call lsyssc_getbase_double(rvelocity_x, p_Dvelocity_x)
  call lsyssc_getbase_double(rvelocity_y, p_Dvelocity_y)
  call lsyssc_getbase_double(radvect,   p_Dadvect)
  call lsyssc_getbase_double(rindicator, p_Dindicator)
  
  ! Get density values
  call ucd_getVariable(rexport, 'density',    p_Ddensity,  nlength)
  call ucd_getVariable(rexport, 'pressure',   p_Dpressure, nlength)
  call ucd_getVariable(rexport, 'velocity_x', p_Dvelocity_x, nlength)
  call ucd_getVariable(rexport, 'velocity_y', p_Dvelocity_y, nlength)
  call ucd_getVariable(rexport, 'advect',     p_Dadvect,   nlength)

  p_Dvelocity = sqrt(p_Dvelocity_x**2 + p_Dvelocity_y**2)
  
  ! Compute indicator
  call calcAdaptationIndicator(rtriangulation, p_Ddensity, p_Dpressure,&
                               p_Dvelocity, p_Dadvect, p_Dindicator)

  ! Export gradient values to UCD file
  call ucd_startGMV(rexportOut, UCD_FLAG_STANDARD, rtriangulation, 'indicator.gmv')

  call ucd_addVariableVertexBased(rexportOut, 'density', UCD_VAR_STANDARD, p_Ddensity)
  call ucd_addVariableVertexBased(rexportOut, 'pressure', UCD_VAR_STANDARD, p_Dpressure)
  call ucd_addVariableVertexBased(rexportOut, 'vec_mag', UCD_VAR_STANDARD, p_Dvelocity)
  call ucd_addVariableVertexBased(rexportOut, 'advect', UCD_VAR_STANDARD, p_Dadvect)
  call ucd_addVariableElementBased(rexportOut, 'indicator', UCD_VAR_STANDARD, p_Dindicator)

  call ucd_write(rexportOut)
  
  ! Release data structures
  call ucd_release(rexport)
  call ucd_release(rexportOut)
  call lsyssc_releaseVector(rdensity)
  call lsyssc_releaseVector(rpressure)
  call lsyssc_releaseVector(rvelocity)
  call lsyssc_releaseVector(rvelocity_x)
  call lsyssc_releaseVector(rvelocity_y)
  call lsyssc_releaseVector(radvect)
  call lsyssc_releaseVector(rindicator)
  call tria_done(rtriangulation)
  call spdiscr_releaseBlockDiscr(rdiscretisation, .true.)

  ! Print out heap statistics - just to check if everything
  ! is cleaned up.
  ! This should display 'Handles in use=0' and 'Memory in use=0'!
  call output_lbrk ()
  call storage_info(.true.)
  
  ! Clean up the storage management, finish
  call storage_done()

contains
  
  !*****************************************************************************
  
  subroutine calcAdaptationIndicator(rtriangulation, Ddensity, Dpressure,&
                                     Dvelocity, Dadvect, Dindicator)

    type(t_triangulation), intent(IN) :: rtriangulation
    real(DP), dimension(:), intent(IN)  :: Ddensity, Dpressure, Dvelocity, Dadvect
    real(DP), dimension(:), intent(OUT) :: Dindicator
          
    ! local variables
    real(DP), dimension(:,:), pointer :: p_DvertexCoords
    integer, dimension(:,:), pointer :: p_IverticesAtElement, p_IneighboursAtElement
    logical, dimension(:), pointer :: p_BisActiveElement
    real(DP) :: dgradient, dlambda_i, dlambda_j
    integer :: iel,jel,ive,nve,i,j,iprotectLayer
        
    real(DP), parameter :: dEpsS = 0.2_DP
    real(DP), parameter :: dEpsC = 0.1_DP
    
    ! Set pointers
    call storage_getbase_double2D(rtriangulation%h_DvertexCoords, p_DvertexCoords)
    call storage_getbase_int2D(rtriangulation%h_IverticesAtElement, p_IverticesAtElement)
    call storage_getbase_int2D(rtriangulation%h_IneighboursAtElement, p_IneighboursAtElement)
    
    ! Loop over all elements
    elements: do iel = 1, rtriangulation%NEL

      ! Get number of vertices at element
      nve = tria_getNVE(p_IverticesAtElement, iel)

      ! Initialize gradient value
      dgradient = 0.0_DP

      ! Check element for shocks and contact discontinuities
      vertices: do ive = 1, nve
        
        ! Get global vertex numbers
        i = p_IverticesAtElement(ive, iel)
        j = p_IverticesAtElement(mod(ive,nve)+1, iel)

        ! Compute nodal values of tracer
        dlambda_i = Dadvect(i)/Ddensity(i)
        dlambda_j = Dadvect(j)/Ddensity(j)
        
        ! Update maximum gradient function
        dgradient = max(dgradient,&
            abs(abs(Ddensity(i))  - abs(Ddensity(j)))  / max(abs(Ddensity(i)),  abs(Ddensity(j))),&
            abs(abs(Dpressure(i)) - abs(Dpressure(j))) / max(abs(Dpressure(i)), abs(Dpressure(j))),&
            abs(abs(dlambda_i)    - abs(dlambda_j))    / max(1e-2, min(abs(dlambda_i), abs(dlambda_j))) )
      end do vertices
      
      ! If we end up here, then the maximum gradient function is adopted
      Dindicator(iel) = dgradient
    end do elements

     
!!$    ! Add protection layers
!!$    allocate(p_BisActiveElement(rtriangulation%NEL))
!!$
!!$    do iprotectLayer = 1, 2
!!$
!!$      p_BisActiveElement = .false.
!!$
!!$      do iel = 1, rtriangulation%NEL
!!$
!!$        if (p_BisactiveElement(iel)) cycle
!!$        if (p_Dindicator(iel) .le. 0.8) cycle
!!$
!!$        do ive = 1, tria_getNVE(p_IverticesAtElement, iel)
!!$          jel = p_IneighboursAtElement(ive, iel)
!!$          if (jel .eq. 0) cycle
!!$          if (p_BisactiveElement(jel)) then
!!$            p_Dindicator(jel) = max(p_Dindicator(jel), p_Dindicator(iel))
!!$          else
!!$            if (p_Dindicator(jel) .lt. 0.8) then
!!$              p_Dindicator(jel) = max(p_Dindicator(jel), p_Dindicator(iel))
!!$              p_BisactiveElement(jel) = .true.
!!$            end if
!!$          end if
!!$        end do
!!$      end do
!!$    end do
!!$
!!$    deallocate(p_BisActiveElement)
    
  end subroutine calcAdaptationIndicator

end program errest
