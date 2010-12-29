!##############################################################################
!# ****************************************************************************
!# <name> geneticalgorithm </name>
!# ****************************************************************************
!#
!# <purpose>
!# This modules provides the basic infrastructure for writing genetic
!# algorithms based on crossover operations and mutation.
!#
!# This module is programmed following the "Introduction to Genetic Algorithms"
!# http://www.obitko.com/tutorials/genetic-algorithms/index.php
!#
!# </purpose>
!##############################################################################
module geneticalgorithm

  use fsystem
  use sort

  implicit none
  
  private

  public :: t_chromosome,t_population
  public :: ga_initPopulation
  public :: ga_releasePopulation
  public :: ga_copyPopulation
  public :: ga_outputPopulation
  public :: ga_sortPopulation
  public :: ga_renewPopulation
  public :: ga_selectChromosome
  public :: ga_initChromosome
  public :: ga_releaseChromosome
  public :: ga_copyChromosome
  public :: ga_outputChromosome
  public :: ga_crossoverChromosome
  public :: ga_mutateChromosome

!<types>
!<typeblock>

  ! A chromosome that can be used to store one single experimen

  type t_chromosome

    ! Fitness level: the larger the better
    real(DP) :: dfitness = 0.0_DP

    ! Character string encoding the DNA
    character, dimension(:), pointer :: p_DNA => null()

  end type t_chromosome
!</typeblock>

!<typeblock>

  ! A population of many chromosomes

  type t_population

    ! Crossover probability: recommendation is 80% - 95%.
    !
    ! This rate determines how often crossover will be performed. If
    ! the crossover probability is 0%, then the whole new generation
    ! is made from exact copies of chromosomes from the old
    ! population. If is is 100%, then all offsping is made by
    ! crossover.
    real(DP) :: dcrossoverrate = 0.85_DP
    
    ! Mutation probability: recommendation is 0.5% - 1%.
    !
    ! This rate determines how often parts of chromosomes will be
    ! mutated. If the mutation probability is 0%, nothing is
    ! changed. If it is 100%, then the whole chromosome is mutated.
    real(DP) :: dmutationrate = 0.01_DP

    ! Elitism probability: 
    !
    ! This rate determines how many of the best chromosomes from the
    ! old population will be saved for the new population.
    real(DP) :: delitismrate = 0.05_DP

    ! Number of chromosomes in population
    integer :: nchromosomes = 0

    ! Set of chromosomes in population
    type(t_chromosome), dimension(:), pointer :: p_Rchromosomes => null()

  end type t_population
!</typeblock>
!</types>

contains

  !************************************************************************
  ! Population routines
  !************************************************************************

!<subroutine>
  subroutine ga_initPopulation(rpopulation, nchromosomes, ndnaLength,&
      dcrossoverrate, dmutationrate, delitismrate)

!<description>
    ! This subroutine initialises a population of chromosomes
!</description>

!<input>
    ! Number of chromosomes in population
    integer, intent(in) :: nchromosomes

    ! Length of DNA string for each chromosome
    integer, intent(in) :: ndnaLength

    ! OPTIONAL: rate of crossover
    real(DP), intent(in), optional :: dcrossoverrate

    ! OPTIONAL: rate of mutation
    real(DP), intent(in), optional :: dmutationrate

    ! OPTIONAL: rate of elitism
    real(DP), intent(in), optional :: delitismrate
!</input>

!<output>
    ! Population of chromosomes
    type(t_population), intent(out) :: rpopulation
!</output>
!</subroutine>

    ! local variables
    integer :: ichromosome,isystemclock,i,isize
    integer, dimension(:), allocatable :: Iseed

    ! Initialise pseudo-random generator
    call system_clock(isystemClock)
    call random_seed(size=isize)
    allocate(Iseed(isize))
    Iseed = isystemClock + 37 * (/ (i - 1, i = 1, isize) /)
    call random_seed(put=Iseed)
    deallocate(Iseed)

    ! Initialisation
    if(present(dcrossoverrate)) rpopulation%dcrossoverrate = dcrossoverrate
    if(present(dmutationrate))  rpopulation%dmutationrate  = dmutationrate
    if(present(delitismrate))   rpopulation%delitismrate   = delitismrate

    rpopulation%nchromosomes = nchromosomes
    allocate(rpopulation%p_Rchromosomes(rpopulation%nchromosomes))

    do ichromosome = 1, rpopulation%nchromosomes
      call ga_initChromosome(rpopulation%p_Rchromosomes(ichromosome),ndnaLength)
    end do
    
  end subroutine ga_initPopulation

  !************************************************************************
  
!<subroutine>
  subroutine ga_releasePopulation(rpopulation)

!<description>
    ! This subroutine releases a population of chromosomes
!</description>
   
!<inputoutput>
    ! Population of chromosomes
    type(t_population), intent(inout) :: rpopulation
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: ichromosome
    
    do ichromosome = 1, rpopulation%nchromosomes
      call ga_releaseChromosome(rpopulation%p_Rchromosomes(ichromosome))
    end do

    deallocate(rpopulation%p_Rchromosomes)
    rpopulation%nchromosomes = 0
    rpopulation%p_Rchromosomes => null()

  end subroutine ga_releasePopulation

  !************************************************************************

!<subroutine>
  subroutine ga_copyPopulation(rpopulationSrc, rpopulationDest)

!<description>
    ! This subroutine copies a population to another population
!</description>

!<input>
    ! Source population of chromosomes
    type(t_population), intent(in) :: rpopulationSrc
!</input>

!<inputoutput>
    ! Destination population of chromosomes
    type(t_population), intent(inout) :: rpopulationDest
!</inputoutput>
!</subroutine>

    ! local variables
    integer :: ichromosome

    if (rpopulationDest%nchromosomes .ne. 0)&
        call ga_releasePopulation(rpopulationDest)

    call ga_initPopulation(rpopulationDest, rpopulationSrc%nchromosomes,&
                           size(rpopulationSrc%p_Rchromosomes(1)%p_DNA),&
                           rpopulationSrc%dcrossoverrate,&
                           rpopulationSrc%dmutationrate,&
                           rpopulationSrc%delitismrate)

    do ichromosome = 1, rpopulationSrc%nchromosomes
      call ga_copyChromosome(rpopulationSrc%p_Rchromosomes(ichromosome),&
                             rpopulationDest%p_Rchromosomes(ichromosome))
    end do

  end subroutine ga_copyPopulation

  !************************************************************************

!<subroutine>
  subroutine ga_outputPopulation(rpopulation, iunit)

!<description>
    ! This subroutine writes out a population of cromosomes
!</description>

!<input>
    ! Population of chromosomes
    type(t_population), intent(in) :: rpopulation

    ! OPTIONAL: number of the output unit
    integer, intent(in), optional :: iunit
!</input>
!</subroutine>

    ! local variables
    integer :: ichromosome

    do ichromosome = 1, rpopulation%nchromosomes
      call ga_outputChromosome(rpopulation%p_Rchromosomes(ichromosome), iunit)
    end do

  end subroutine ga_outputPopulation

  !************************************************************************

!<subroutine>
  subroutine ga_sortPopulation(rpopulation)

!<description>
    ! This subroutine sorts the chromosomes of a population according
    ! to their fitness level (the larger the better)
!</description>

!<inputoutpu>
    ! Population of chromosomes
    type(t_population), intent(inout) :: rpopulation
!</inputoutpu>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: Dfitness
    integer, dimension(:), pointer :: Imapping
    integer :: ichromosome
    
    allocate(Dfitness(rpopulation%nchromosomes))
    allocate(Imapping(rpopulation%nchromosomes))

    do ichromosome = 1, rpopulation%nchromosomes
      Dfitness(ichromosome) = rpopulation%p_Rchromosomes(ichromosome)%dfitness
      Imapping(ichromosome) = ichromosome
    end do

    call sort_dp(Dfitness, SORT_QUICK, Imapping)
    Imapping = Imapping(rpopulation%nchromosomes:1:-1)
    rpopulation%p_Rchromosomes = rpopulation%p_Rchromosomes(Imapping)
    
    deallocate(Dfitness, Imapping)

  end subroutine ga_sortPopulation

  !************************************************************************

!<subroutine>
  subroutine ga_renewPopulation(rpopulation)

!<description>
    ! This subroutine generates a new population from the given population
!</description>

!<inputoutput>
    ! Source and destination population of chromosomes
    type(t_population), intent(inout) :: rpopulation
!</inputoutput>
!</subroutine>
    
    ! local variables
    type(t_population) :: rpopulationTmp
    integer :: ichromosome,ioffset,isrc1,isrc2
    
    ! Sort the population according to the fitness levels
    call ga_sortPopulation(rpopulation)

    ! Step 1: Replace lower population half by upper population half (natural selection)
    ioffset = floor(rpopulation%nchromosomes/2.0)
    do ichromosome = 1, ioffset
      call ga_copyChromosome(rpopulation%p_Rchromosomes(ichromosome),&
                             rpopulation%p_Rchromosomes(ioffset+ichromosome))
    end do
    
    ! Make a copy of the population
    call ga_copyPopulation(rpopulation,rpopulationTmp)

    ! Step 2: Randomly blend 1st quarter chromosomes into 4th quarter (crossover)
    do ichromosome = floor(3.0*rpopulation%nchromosomes/4.0)+1,&
                     rpopulation%nchromosomes

      ! Select parent chromosome from first quarter
      isrc1 = ga_selectChromosome(rpopulationTmp,&
                1, floor(rpopulation%nchromosomes/4.0))

      ! Select parent chromosome from fourth quarter
      isrc2 = ga_selectChromosome(rpopulationTmp,&
                floor(3.0*rpopulation%nchromosomes/4.0)+1,&
                rpopulation%nchromosomes)

      ! Generate new chromosome by crossover of parents
      call ga_crossoverChromosome(rpopulationTmp%p_Rchromosomes(isrc1),&
                                  rpopulationTmp%p_Rchromosomes(isrc2),&
                                  rpopulation%p_Rchromosomes(ichromosome))
    end do

    ! Step 3: Randomly blend 2nd quarter chromosomes into 3rd quarter (crossover)
    do ichromosome = floor(rpopulation%nchromosomes/2.0)+1,&
                     floor(3.0*rpopulation%nchromosomes/4.0)

      ! Select parent chromosome from second quarter
      isrc1 = ga_selectChromosome(rpopulationTmp,&
                floor(rpopulation%nchromosomes/4.0)+1,&
                floor(rpopulation%nchromosomes/2.0))

      ! Select parent chromosome from third quarter
      isrc2 = ga_selectChromosome(rpopulationTmp,&
                floor(rpopulation%nchromosomes/2.0)+1,&
                floor(3.0*rpopulation%nchromosomes/4.0))

      ! Generate new chromosome by crossover of parents
      call ga_crossoverChromosome(rpopulationTmp%p_Rchromosomes(isrc1),&
                                  rpopulationTmp%p_Rchromosomes(isrc2),&
                                  rpopulation%p_Rchromosomes(ichromosome))
    end do

    ! Step 4: Mutate lower 3/4 of population
    if (rpopulation%dmutationrate .ne. 0.0_DP) then
      do ichromosome = floor(rpopulation%nchromosomes/4.0),&
                       rpopulation%nchromosomes
        call ga_mutateChromosome(rpopulation%p_Rchromosomes(ichromosome),&
                                 rpopulation%dmutationrate)
      end do
    end if

    ! Release temporal copy of population
    call ga_releasePopulation(rpopulationTmp)

!!$    ! Step 1: Copy the best chromosomes
!!$    nbestchromosomes = min(rpopulationTmp%nchromosomes,&
!!$                           ceiling(rpopulationTmp%nchromosomes*&
!!$                                   rpopulationTmp%delitismrate))
!!$
!!$    do ichromosome = 1, nbestchromosomes
!!$      call ga_copyChromosome(rpopulationTmp%p_Rchromosomes(ichromosome),&
!!$                             rpopulation%p_Rchromosomes(ichromosome))
!!$    end do
!!$    
!!$    ! Step 2: Generate new chromosomes by crossover and mutation
!!$    do ichromosome = nbestchromosomes+1, min(rpopulationTmp%nchromosomes,&
!!$                                             ceiling(rpopulationTmp%nchromosomes*&
!!$                                                     rpopulationTmp%dcrossoverrate))
!!$      ! Select two parent chromosomes
!!$      isrc1 = ga_selectChromosome(rpopulationTmp)
!!$      isrc2 = ga_selectChromosome(rpopulationTmp)
!!$
!!$      ! Generate new chromosome by crossover of parents
!!$      call ga_crossoverChromosome(rpopulationTmp%p_Rchromosomes(isrc1),&
!!$                                  rpopulationTmp%p_Rchromosomes(isrc2),&
!!$                                  rpopulation%p_Rchromosomes(ichromosome))
!!$      ! Mutate new chromosome
!!$      if (rpopulation%dmutationrate .ne. 0.0_DP) &
!!$          call ga_mutateChromosome(rpopulation%p_Rchromosomes(ichromosome),&
!!$                                   rpopulation%dmutationrate)
!!$    end do
!!$
!!$    ! Step 3: Generate new chromosomes randomly
!!$    do ichromosome = ceiling(rpopulationTmp%nchromosomes*&
!!$                             rpopulationTmp%dcrossoverrate), rpopulationTmp%nchromosomes
!!$      call ga_initChromosome(rpopulation%p_Rchromosomes(ichromosome),&
!!$                             size(rpopulationTmp%p_Rchromosomes(1)%p_DNA))
!!$    end do
!!$
!!$    ! Release temporal copy of population
!!$    call ga_releasePopulation(rpopulationTmp)

  end subroutine ga_renewPopulation

  !************************************************************************

!<function>
  function ga_selectChromosome(rpopulation, ifirst, ilast) result(ichromosome)

!<description>
    ! This function selects a chromosome from the population by the
    ! roulette wheel selection method (i.e. parents are selected
    ! according to their fitness. The better the chromosomes are, the
    ! more chances to be selected they have.)
!</description>

!<input>
    ! Population of chromosomes
    type(t_population), intent(in) :: rpopulation

    ! OPTIONAL: First and last chromosome to consider
    integer, intent(in), optional :: ifirst,ilast
!</input>

!<result>
    ! Number of the selected chromosome
    integer :: ichromosome
!</result>
!</function>
    
    ! local variables
    real(DP) :: dtotalFitness, dfitnessBound
    integer :: ifirstChromosome,ilastChromosome

    ! Initialisation
    ifirstChromosome = 1
    ilastChromosome  = rpopulation%nchromosomes
    if (present(ifirst)) ifirstChromosome = ifirst
    if (present(ilast))  ilastChromosome  = ilast

    ! Compute total fitness level
    dtotalFitness = 0.0_DP
    do ichromosome = ifirstChromosome, ilastChromosome
      dtotalFitness = dtotalFitness +&
          rpopulation%p_Rchromosomes(ichromosome)%dfitness
    end do
    
    ! Compute bound for fitness level
    call random_number(dfitnessBound)
    dfitnessBound = dfitnessBound*dtotalFitness

    ! Roulette wheel selection
    dtotalFitness = 0.0_DP
    do ichromosome = ifirstChromosome, ilastChromosome
      dtotalFitness = dtotalFitness +&
          rpopulation%p_Rchromosomes(ichromosome)%dfitness
      if (dtotalFitness .gt. dfitnessBound) return
    end do
    
    ! Return selected chromosome
    ichromosome = ilastChromosome
    
  end function ga_selectChromosome

  !************************************************************************
  ! Chromosome routines
  !************************************************************************

!<subroutine>
  subroutine ga_initChromosome(rchromosome, ndnaLength)

!<description>
    ! This subroutine initialises a chromosome
!</description>

!<input>
    ! Length of the DNA string
    integer, intent(in) :: ndnaLength
!</input>

!<output>
    ! Chromosome
    type(t_chromosome), intent(out) :: rchromosome
!</output>
!</subroutine>

    ! local variables
    real(DP), dimension(:), pointer :: Drandom
    integer :: idna

    rchromosome%dfitness = 0.0_DP
    allocate(rchromosome%p_DNA(ndnaLength))

    ! Generate pseudo-random number
    allocate(Drandom(ndnaLength))
    call random_number(Drandom)

    ! Reformat random number into DNA string
    do idna = 1, ndnaLength
      rchromosome%p_DNA(idna) = merge ('1','0', Drandom(idna) .ge. 0.5_DP)
    end do
    deallocate(Drandom)

  end subroutine ga_initChromosome

  !************************************************************************

!<subroutine>
  subroutine ga_releaseChromosome(rchromosome)

!<description>
    ! This subroutine releases a chromosome
!</description>

!<inputoutput>
    ! Chromosome
    type(t_chromosome), intent(inout) :: rchromosome
!</inputoutput>
!</subroutine>

    rchromosome%dfitness = 0.0_DP
    deallocate(rchromosome%p_DNA)
    rchromosome%p_DNA => null()

  end subroutine ga_releaseChromosome

  !************************************************************************

!<subroutine>
  subroutine ga_copyChromosome(rchromosomeSrc, rchromosomeDest)

!<description>
    ! This subroutine copies a chromosome to another chromosome
!</description>

!<input>
    ! Source chromosome
    type(t_chromosome), intent(in) :: rchromosomeSrc
!</input>

!<inputoutput>
    ! Destination chromosome
    type(t_chromosome), intent(inout) :: rchromosomeDest
!</inputoutput>
!</subroutine>

    if (associated(rchromosomeDest%p_DNA))&
        call ga_releaseChromosome(rchromosomeDest)

    call ga_initChromosome(rchromosomeDest, size(rchromosomeSrc%p_DNA))
    rchromosomeDest%dfitness = rchromosomeSrc%dfitness
    rchromosomeDest%p_DNA    = rchromosomeSrc%p_DNA

  end subroutine ga_copyChromosome

  !************************************************************************

!<subroutine>
  subroutine ga_outputChromosome(rchromosome, iunit)
    
!<description>
    ! This subroutine writes out a chromosome
!</description>

!<input>
    ! Chromosome
    type(t_chromosome), intent(in) :: rchromosome

    ! OPTIONAL: number of the output unit
    integer, intent(in), optional :: iunit
!</input>
!</subroutine>

    ! local variables
    integer :: idna

    if (present(iunit)) then
      write(iunit,fmt='(G12.4)', advance='no') rchromosome%dfitness
      do idna = 1, size(rchromosome%p_DNA)
        write(iunit,fmt='(A)', advance='no') rchromosome%p_DNA(idna)
      end do
      write(iunit,*)
    else
      write(*,fmt='(G12.4)', advance='no') rchromosome%dfitness
      do idna = 1, size(rchromosome%p_DNA)
        write(*,fmt='(A)', advance='no') rchromosome%p_DNA(idna)
      end do
      write(*,*)
    end if
      
  end subroutine ga_outputChromosome

  !************************************************************************

!<subroutine>
  subroutine ga_crossoverChromosome(rchromosomeSrc1, rchromosomeSrc2, rchromosomeDest)

!<description>
    ! This subroutine generates a new chromosome based on two parent
    ! cromosomes by two point crossover. That is, a random position
    ! ipos is determined and DNA(1:ipos) is copied from the first
    ! parent chromosome, whereas DNA(ipos+1:end) is copied from the
    ! second parent chromosome.
!</description>

!<input>
    ! Parent chromosomes
    type(t_chromosome), intent(in) :: rchromosomeSrc1,rchromosomeSrc2
!</input>

!<inputoutput>
    ! Crossover
    type(t_chromosome), intent(inout) :: rchromosomeDest
!</inputoutput>
!</subroutine>

    ! local variables
    real(DP), dimension(2) :: Drandom
    integer, dimension(2) :: Idna
    integer :: ndnaLength

    rchromosomeDest%dfitness = 0.0_DP
    
    ndnaLength = min(size(rchromosomeSrc1%p_DNA),&
                     size(rchromosomeSrc2%p_DNA),&
                     size(rchromosomeDest%p_DNA))
    
    ! Generate random crossover points
    call random_number(Drandom)
    if (Drandom(1) .gt. Drandom(2)) Drandom(1:2) = Drandom(2:1:-1)
    Idna = int(Drandom*(ndnaLength))+1

    ! Crossover parent chromosomes
    rchromosomeDest%p_DNA(1:Idna(1))         = rchromosomeSrc1%p_DNA(1:Idna(1))
    rchromosomeDest%p_DNA(Idna(1)+1:Idna(2)) = rchromosomeSrc2%p_DNA(Idna(1)+1:Idna(2))
    rchromosomeDest%p_DNA(Idna(2)+1:)        = rchromosomeSrc1%p_DNA(Idna(2)+1:)

  end subroutine ga_crossoverChromosome

  !************************************************************************

!<subroutine>
  subroutine ga_mutateChromosome(rchromosome, dmutationrate)

!<description>
    ! This subroutine mutates a given chromosome, i.e. flips bits,
    ! according to the given mutation probability
!</description>

!<input>
    ! Mutation probability
    real(DP), intent(in) :: dmutationrate
!</input>

!<inputoutput>
    ! Chromosome
    type(t_chromosome), intent(inout) :: rchromosome
!</inputoutput>
!</subroutine>    

    real(DP) :: drandom
    integer :: idna,ndnaLength,imutation,nmutations

    rchromosome%dfitness = 0.0_DP

    ndnaLength = size(rchromosome%p_DNA)
    nmutations = max(1,ceiling(ndnaLength*dmutationrate))

    do imutation = 1, nmutations
      call random_number(drandom)
      idna = int(drandom*ndnaLength)+1
      rchromosome%p_DNA(idna:idna) = merge('0','1', rchromosome%p_DNA(idna:idna).eq.'1')
    end do

  end subroutine ga_mutateChromosome

end module geneticalgorithm
