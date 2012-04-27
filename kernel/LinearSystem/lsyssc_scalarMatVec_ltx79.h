  ! We have to perform matrix*vector + vector.
  ! What we actually calculate here is:
  !    y  =  cx * A^t * x  +  ( cy * y )
  !       =  ( (cx * x)^t * A  +  (cy * y)^t )^t.

  ! Unfortunately, if cy != 1, then we need to scale y now.
  if(cy .eq. __VecZero__) then
    ! Clear y
    call lalg_clearVector (Dy)
  else if(cy .ne. __VecOne__) then
    ! Scale y
    call lalg_scaleVector(Dy, cy)
  end if

#ifdef __NVAR__

  ! REMARK: OpenMP does not work for transposed matrix-vector multiplication
  ! in the way it is realised below!!!

  ! Perform the multiplication.
  if (cx .ne. __MatOne__) then

    !DOES NOT WORK: omp parallel do private(ia,icol,Ddtmp) default(shared) &
    !DOES NOT WORK: omp if(NEQ > rperfconfig%NEQMIN_OMP)
    do irow = 1, NEQ
      do ivar = 1, __NVAR__
        Ddtmp(ivar) = cx*Dx(__NVAR__*(irow-1)+ivar)
      end do
      do ia = Kld(irow), Kld(irow+1)-1
        icol = Kcol(ia)
        do ivar = 1, __NVAR__
          Dy(__NVAR__*(icol-1)+ivar) = Dy(__NVAR__*(icol-1)+ivar)&
                                     + Da(ia)*Ddtmp(ivar)
        end do
      end do
    end do
    !DOES NOT WORK: omp end parallel do

  else   ! cx = 1.0

    !DOES NOT WORK: omp parallel do private(ia,icol,Ddtmp) default(shared) &
    !DOES NOT WORK: omp if(NEQ > rperfconfig%NEQMIN_OMP)
    do irow = 1, NEQ
      do ia = Kld(irow), Kld(irow+1)-1
        icol = Kcol(ia)
        do ivar = 1, __NVAR__
          Dy(__NVAR__*(icol-1)+ivar) = Dy(__NVAR__*(icol-1)+ivar)&
                                     + Da(ia)*Dx(__NVAR__*(irow-1)+ivar)
        end do
      end do
    end do
    !DOES NOT WORK: omp end parallel do

  end if

#else

  ! Perform the multiplication.
  if (cx .ne. __MatOne__) then

    !DOES NOT WORK: omp parallel do private(ia,icol,dtmp) default(shared) &
    !DOES NOT WORK: omp if(NEQ > rperfconfig%NEQMIN_OMP)
    do irow = 1, NEQ
      dtmp = cx*Dx(irow)
      do ia = Kld(irow), Kld(irow+1)-1
        icol = Kcol(ia)
        Dy(icol) = Dy(icol) + Da(ia)*dtmp
      end do
    end do
    !DOES NOT WORK: omp end parallel do

  else   ! cx = 1.0

    !DOES NOT WORK: omp parallel do private(ia,icol,dtmp) default(shared) &
    !DOES NOT WORK: omp if(NEQ > rperfconfig%NEQMIN_OMP)
    do irow = 1, NEQ
      do ia = Kld(irow), Kld(irow+1)-1
        icol = Kcol(ia)
        Dy(icol) = Dy(icol) + Da(ia)*Dx(irow)
      end do
    end do
    !DOES NOT WORK: omp end parallel do

  end if

#endif
