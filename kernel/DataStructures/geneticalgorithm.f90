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
    integer :: ichromosome,isystemclock

    ! Initialise random number generator
    call system_clock(isystemclock)
    call random_seed(isystemclock)

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
  subroutine ga_outputPopulation(rpopulation)

!<description>
    ! This subroutine writes out a population of cromosomes
!</description>

!<input>
    ! Population of chromosomes
    type(t_population), intent(in) :: rpopulation
!</input>
!</subroutine>

    ! local variables
    integer :: ichromosome

    do ichromosome = 1, rpopulation%nchromosomes
      call ga_outputChromosome(rpopulation%p_Rchromosomes(ichromosome))
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
    integer :: ichromosome,nbestchromosomes,isrc1,isrc2
    
    ! Make a copy of the old population
    call ga_copyPopulation(rpopulation,rpopulationTmp)

    ! Sort the old population according to the fitness levels
    call ga_sortPopulation(rpopulationTmp)
    
    ! Step 1: Copy the best chromosomes
    nbestchromosomes = min(rpopulationTmp%nchromosomes,&
                           ceiling(rpopulationTmp%nchromosomes*&
                                   rpopulationTmp%delitismrate))

    do ichromosome = 1, nbestchromosomes
      call ga_copyChromosome(rpopulationTmp%p_Rchromosomes(ichromosome),&
                             rpopulation%p_Rchromosomes(ichromosome))
    end do
    
    ! Step 2: Generate new chromosomes by crossover and mutation
    do ichromosome = nbestchromosomes+1, min(rpopulationTmp%nchromosomes,&
                                             ceiling(rpopulationTmp%nchromosomes*&
                                                     rpopulationTmp%dcrossoverrate))
      ! Select two parent chromosomes
      isrc1 = ga_selectChromosome(rpopulationTmp)
      isrc2 = ga_selectChromosome(rpopulationTmp)

      ! Generate new chromosome by crossover of parents
      call ga_crossoverChromosome(rpopulationTmp%p_Rchromosomes(isrc1),&
                                  rpopulationTmp%p_Rchromosomes(isrc2),&
                                  rpopulation%p_Rchromosomes(ichromosome))
      ! Mutate new chromosome
      if (rpopulation%dmutationrate .ne. 0.0_DP) &
          call ga_mutateChromosome(rpopulation%p_Rchromosomes(ichromosome),&
                                   rpopulation%dmutationrate)
    end do

    ! Step 3: Generate new chromosomes randomly
    do ichromosome = ceiling(rpopulationTmp%nchromosomes*&
                             rpopulationTmp%dcrossoverrate), rpopulationTmp%nchromosomes
      call ga_initChromosome(rpopulation%p_Rchromosomes(ichromosome),&
                             size(rpopulationTmp%p_Rchromosomes(1)%p_DNA))
    end do

    ! Release temporal copy of population
    call ga_releasePopulation(rpopulationTmp)

  end subroutine ga_renewPopulation

  !************************************************************************

!<function>
  function ga_selectChromosome(rpopulation) result(ichromosome)

!<description>
    ! This function selects a chromosome from the population by the
    ! roulette wheel selection method (i.e. parents are selected
    ! according to their fitness. The better the chromosomes are, the
    ! more chances to be selected they have.)
!</description>

!<input>
    ! Population of chromosomes
    type(t_population), intent(in) :: rpopulation
!</input>

!<result>
    ! Number of the selected chromosome
    integer :: ichromosome
!</result>
!</function>
    
    ! local variables
    real(DP) :: dtotalFitness, dfitnessBound

    ! Compute total fitness level
    dtotalFitness = 0.0_DP
    do ichromosome = 1, rpopulation%nchromosomes
      dtotalFitness = dtotalFitness +&
          rpopulation%p_Rchromosomes(ichromosome)%dfitness
    end do
    
    ! Compute bound for fitness level
    call random_number(dfitnessBound)
    if (dfitnessBound .eq. 0.0_DP) dfitnessBound = SYS_EPSREAL
    dfitnessBound = dfitnessBound*dtotalFitness

    ! Roulette wheel selection
    dtotalFitness = 0.0_DP
    do ichromosome = 1, rpopulation%nchromosomes
      dtotalFitness = dtotalFitness +&
          rpopulation%p_Rchromosomes(ichromosome)%dfitness
      if (dtotalFitness .gt. dfitnessBound) return
    end do
    
    ! Return selected chromosome
    ichromosome = rpopulation%nchromosomes
    
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

    allocate(Drandom(ndnaLength))
    call random_number(Drandom)
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
  subroutine ga_outputChromosome(rchromosome)
    
!<description>
    ! This subroutine writes out a chromosome
!</description>

!<input>
    ! Chromosome
    type(t_chromosome), intent(in) :: rchromosome
!</input>
!</subroutine>

    ! local variables
    integer :: idna

    do idna = 1, size(rchromosome%p_DNA)
      write(*,fmt='(A)', advance='no') rchromosome%p_DNA(idna)
    end do
    write(*,fmt='(3X,F5.3)') rchromosome%dfitness

  end subroutine ga_outputChromosome

  !************************************************************************

!<subroutine>
  subroutine ga_crossoverChromosome(rchromosomeSrc1, rchromosomeSrc2, rchromosomeDest)

!<description>
    ! This subroutine generates a new chromosome based on two parent
    ! cromosomes by single point crossover. That is, a random position
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
    real(DP) :: drandom
    integer :: idna,ndnaLength

    rchromosomeDest%dfitness = 0.0_DP
    
    ndnaLength = min(size(rchromosomeSrc1%p_DNA),&
                     size(rchromosomeSrc2%p_DNA),&
                     size(rchromosomeDest%p_DNA))
    call random_number(drandom)
    idna = min(ceiling(ndnaLength/drandom),ndnaLength)

    rchromosomeDest%p_DNA(1:idna)            = rchromosomeSrc1%p_DNA(1:idna)
    rchromosomeDest%p_DNA(idna+1:ndnaLength) = rchromosomeSrc2%p_DNA(idna+1:ndnaLength)
    rchromosomeDest%p_DNA(ndnaLength+1:)     = '0'

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
    nmutations = max(1,ceiling(ndnaLength/100*dmutationrate))

    do imutation = 1, nmutations
      call random_number(drandom)
      idna = min(ceiling(ndnaLength/drandom),ndnaLength)

      rchromosome%p_DNA(idna:idna) = merge('0','1', rchromosome%p_DNA(idna:idna).eq.'1')
    end do

  end subroutine ga_mutateChromosome

end module geneticalgorithm
