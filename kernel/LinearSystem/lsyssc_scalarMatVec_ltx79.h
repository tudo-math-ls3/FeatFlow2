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
  if (cx .eq. __MatOne__) then

    !$omp parallel do private(ia,icol,ivar,Ddtmp) default(shared) &
    !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
    do irow = 1, NEQ
      do ia = Kld(irow), Kld(irow+1)-1
        icol = Kcol(ia)
        do ivar = 1, __NVAR__
          !$omp atomic
          Dy(__NVAR__*(icol-1)+ivar) = Dy(__NVAR__*(icol-1)+ivar)&
                                     + Da(ia)*Dx(__NVAR__*(irow-1)+ivar)
        end do
      end do
    end do
    !$omp end parallel do

  else if (cx .eq. -__MatOne__) then

    !$omp parallel do private(ia,icol,ivar,Ddtmp) default(shared) &
    !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
    do irow = 1, NEQ
      do ia = Kld(irow), Kld(irow+1)-1
        icol = Kcol(ia)
        do ivar = 1, __NVAR__
          !$omp atomic
          Dy(__NVAR__*(icol-1)+ivar) = Dy(__NVAR__*(icol-1)+ivar)&
                                     - Da(ia)*Dx(__NVAR__*(irow-1)+ivar)
        end do
      end do
    end do
    !$omp end parallel do

  else

    !$omp parallel do private(ia,icol,ivar,Ddtmp) default(shared) &
    !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
    do irow = 1, NEQ
      do ivar = 1, __NVAR__
        Ddtmp(ivar) = cx*Dx(__NVAR__*(irow-1)+ivar)
      end do
      do ia = Kld(irow), Kld(irow+1)-1
        icol = Kcol(ia)
        do ivar = 1, __NVAR__
          !$omp atomic
          Dy(__NVAR__*(icol-1)+ivar) = Dy(__NVAR__*(icol-1)+ivar)&
                                     + Da(ia)*Ddtmp(ivar)
        end do
      end do
    end do
    !$omp end parallel do

  end if

#else

  ! Perform the multiplication.
  if (cx .eq. __MatOne__) then

    !$omp parallel do private(ia,icol) default(shared) &
    !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
    do irow = 1, NEQ
      do ia = Kld(irow), Kld(irow+1)-1
        icol = Kcol(ia)
        !$omp atomic
        Dy(icol) = Dy(icol) + Da(ia)*Dx(irow)
      end do
    end do
    !$omp end parallel do

  else if (cx .eq. -__MatOne__) then

    !$omp parallel do private(ia,icol) default(shared) &
    !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
    do irow = 1, NEQ
      do ia = Kld(irow), Kld(irow+1)-1
        icol = Kcol(ia)
        !$omp atomic
        Dy(icol) = Dy(icol) - Da(ia)*Dx(irow)
      end do
    end do
    !$omp end parallel do

  else

    !$omp parallel do private(ia,icol,dtmp) default(shared) &
    !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
    do irow = 1, NEQ
      dtmp = cx*Dx(irow)
      do ia = Kld(irow), Kld(irow+1)-1
        icol = Kcol(ia)
        !$omp atomic
        Dy(icol) = Dy(icol) + Da(ia)*dtmp
      end do
    end do
    !$omp end parallel do

  end if

#endif
