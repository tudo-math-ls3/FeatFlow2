
program gpusd

use gpusd_test1

implicit none

  ! initialise FEAT2 system
  call system_init()
  
  ! initialise FEAT2 output
  call output_init('./log/output.txt')
  
  ! initialise FEAT2 storage
  call storage_init(999, 100)

  ! call test #1
  call output_lbrk()
  call output_line('GPU-SD: Test #1')
  call output_line('---------------')
  call gpusd_1()

  ! print storage statistics
  call output_lbrk()
  call storage_info(.true.)
  
  ! release storage
  call storage_done()

end program