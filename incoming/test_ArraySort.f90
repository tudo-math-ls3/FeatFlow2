PROGRAM test_ArraySort

  USE fsystem
  USE ArraySort
  
  INTEGER(I32), PARAMETER :: nindex = 3
  INTEGER(I32), PARAMETER :: nnode = 100
  
  INTEGER(I32), DIMENSION(nindex, nnode) :: Ielem
  INTEGER(I32) :: iindex, inode
  REAL(DP) :: tmp

  
  CALL RANDOM_SEED
  
  DO inode=1, nnode
    DO iindex=1, nindex
      CALL RANDOM_NUMBER(tmp)
      Ielem(iindex,inode) = INT (tmp*100,I32)
    END DO
    PRINT *, Ielem(:,inode)
  END DO
  
  DO iindex=1, nindex
    PRINT *, ''
    PRINT *, 'Sort by index:',iindex
    CALL arraySort_sortByIndex(Ielem,iindex,nnode,nindex,SORT_STABLE)
    DO inode=1, nnode
      PRINT *, Ielem(:,inode)
    END DO
  END DO

END PROGRAM 
