program mgrenum

use fsystem
use storage
use genoutput
use paramlist
use mgrenum2d_test1

implicit none

! variables
type(t_parlist) :: rparam
character(LEN=256) :: slogfile

  ! init FEAT2
  call sys_init()
  call storage_init(999, 100)
  
  ! read in parameters
  call parlst_init(rparam)
  call parlst_readFromFile(rparam, './data/mgrenum.dat')

  ! init output system
  call parlst_getvalue_string(rparam, '', 'SLOGFILE', slogfile, './log/output.log')
  call output_init(slogfile)

  !!!
  call mgrenum2d_1(rparam)

  ! release parameter list
  call parlst_done(rparam)

  ! release FEAT2
  call output_lbrk ()
  call storage_info(.true.)
  call storage_done()

end program