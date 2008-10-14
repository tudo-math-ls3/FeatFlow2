program test_ArraySort

  use fsystem
  use ArraySort
  
  integer(I32), parameter :: nindex = 3
  integer(I32), parameter :: nnode = 100
  
  integer(I32), dimension(nindex, nnode) :: Ielem
  integer(I32) :: iindex, inode
  real(DP) :: tmp

  
  call RANDOM_SEED
  
  do inode=1, nnode
    do iindex=1, nindex
      call random_number(tmp)
      Ielem(iindex,inode) = int (tmp*100,I32)
    end do
    print *, Ielem(:,inode)
  end do
  
  do iindex=1, nindex
    print *, ''
    print *, 'Sort by index:',iindex
    call arraySort_sortByIndex(Ielem,iindex,nnode,nindex,SORT_STABLE)
    do inode=1, nnode
      print *, Ielem(:,inode)
    end do
  end do

end program 
