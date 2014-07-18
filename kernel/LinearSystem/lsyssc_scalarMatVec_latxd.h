#ifdef __NVAR__

  ! Perform the multiplication
  if (cx .eq. __MatOne__) then

    if(cy .eq. __VecZero__) then

      !$omp parallel do default(shared) if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1,NEQ
        do ivar = 1, __NVAR__
          Dy(__NVAR__*(irow-1)+ivar) = Da(irow)*Dx(__NVAR__*(irow-1)+ivar)
        end do
      end do
      !$omp end parallel do

    else if(cy .eq. __VecOne__) then

      !$omp parallel do default(shared) if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1,NEQ
        do ivar = 1, __NVAR__
          Dy(__NVAR__*(irow-1)+ivar) = Dy(__NVAR__*(irow-1)+ivar)&
                                     + Da(irow)*Dx(__NVAR__*(irow-1)+ivar)
        end do
      end do
      !$omp end parallel do

    else if(cy .eq. -__VecOne__) then

      !$omp parallel do default(shared) if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1,NEQ
        do ivar = 1, __NVAR__
          Dy(__NVAR__*(irow-1)+ivar) = -Dy(__NVAR__*(irow-1)+ivar)&
                                     + Da(irow)*Dx(__NVAR__*(irow-1)+ivar)
        end do
      end do
      !$omp end parallel do

    else   ! arbitrary cy value

      !$omp parallel do default(shared) if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1,NEQ
        do ivar = 1, __NVAR__
          Dy(__NVAR__*(irow-1)+ivar) = cy*Dy(__NVAR__*(irow-1)+ivar)&
                                     + Da(irow)*Dx(__NVAR__*(irow-1)+ivar)
        end do
      end do
      !$omp end parallel do

    end if

  else if (cx .eq. -__MatOne__) then

    if(cy .eq. __VecZero__) then

      !$omp parallel do default(shared) if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1,NEQ
        do ivar = 1, __NVAR__
          Dy(__NVAR__*(irow-1)+ivar) = -Da(irow)*Dx(__NVAR__*(irow-1)+ivar)
        end do
      end do
      !$omp end parallel do

    else if(cy .eq. __VecOne__) then

      !$omp parallel do default(shared) if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1,NEQ
        do ivar = 1, __NVAR__
          Dy(__NVAR__*(irow-1)+ivar) = Dy(__NVAR__*(irow-1)+ivar)&
                                     - Da(irow)*Dx(__NVAR__*(irow-1)+ivar)
        end do
      end do
      !$omp end parallel do

    else if(cy .eq. -__VecOne__) then

      !$omp parallel do default(shared) if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1,NEQ
        do ivar = 1, __NVAR__
          Dy(__NVAR__*(irow-1)+ivar) = -Dy(__NVAR__*(irow-1)+ivar)&
                                     - Da(irow)*Dx(__NVAR__*(irow-1)+ivar)
        end do
      end do
      !$omp end parallel do

    else   ! arbitrary cy value

      !$omp parallel do default(shared) if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1,NEQ
        do ivar = 1, __NVAR__
          Dy(__NVAR__*(irow-1)+ivar) = cy*Dy(__NVAR__*(irow-1)+ivar)&
                                     - Da(irow)*Dx(__NVAR__*(irow-1)+ivar)
        end do
      end do
      !$omp end parallel do

    end if

  else   ! arbitrary cx value

    if (cy .eq. __VecZero__) then

      !$omp parallel do default(shared) if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1,NEQ
        do ivar = 1, __NVAR__
          Dy(__NVAR__*(irow-1)+ivar) = cx*Da(irow)*Dx(__NVAR__*(irow-1)+ivar)
        end do
      end do
      !$omp end parallel do

    else if(cy .eq. __VecOne__) then

      !$omp parallel do default(shared) if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1,NEQ
        do ivar = 1, __NVAR__
          Dy(__NVAR__*(irow-1)+ivar) = Dy(__NVAR__*(irow-1)+ivar)&
                                     + cx*Da(irow)*Dx(__NVAR__*(irow-1)+ivar)
        end do
      end do
      !$omp end parallel do

    else if(cy .eq. -__VecOne__) then

      !$omp parallel do default(shared) if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1,NEQ
        do ivar = 1, __NVAR__
          Dy(__NVAR__*(irow-1)+ivar) = -Dy(__NVAR__*(irow-1)+ivar)&
                                     + cx*Da(irow)*Dx(__NVAR__*(irow-1)+ivar)
        end do
      end do
      !$omp end parallel do

    else   ! arbitrary cy value

      !$omp parallel do default(shared) if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1,NEQ
        do ivar = 1, __NVAR__
          Dy(__NVAR__*(irow-1)+ivar) = cy*Dy(__NVAR__*(irow-1)+ivar)&
                                     + cx*Da(irow)*Dx(__NVAR__*(irow-1)+ivar)
        end do
      end do
      !$omp end parallel do

    end if

  end if

#else

  ! Perform the multiplication
  if (cx .eq. __MatOne__) then

    if(cy .eq. __VecZero__) then

      !$omp parallel do default(shared) if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1,NEQ
        Dy(irow) = Da(irow)*Dx(irow)
      end do
      !$omp end parallel do

    else if(cy .eq. __VecOne__) then

      !$omp parallel do default(shared) if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1,NEQ
        Dy(irow) = Dy(irow) + Da(irow)*Dx(irow)
      end do
      !$omp end parallel do

    else if(cy .eq. -__VecOne__) then

      !$omp parallel do default(shared) if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1,NEQ
        Dy(irow) = -Dy(irow) + Da(irow)*Dx(irow)
      end do
      !$omp end parallel do

    else   ! arbitrary cy value

      !$omp parallel do default(shared) if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1,NEQ
        Dy(irow) = cy*Dy(irow) + Da(irow)*Dx(irow)
      end do
      !$omp end parallel do

    end if

  else if (cx .eq. -__MatOne__) then

    if(cy .eq. __VecZero__) then

      !$omp parallel do default(shared) if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1,NEQ
        Dy(irow) = -Da(irow)*Dx(irow)
      end do
      !$omp end parallel do

    else if(cy .eq. __VecOne__) then

      !$omp parallel do default(shared) if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1,NEQ
        Dy(irow) = Dy(irow) - Da(irow)*Dx(irow)
      end do
      !$omp end parallel do

    else if(cy .eq. -__VecOne__) then

      !$omp parallel do default(shared) if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1,NEQ
        Dy(irow) = -Dy(irow) - Da(irow)*Dx(irow)
      end do
      !$omp end parallel do

    else   ! arbitrary cy value

      !$omp parallel do default(shared) if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1,NEQ
        Dy(irow) = cy*Dy(irow) + Da(irow)*Dx(irow)
      end do
      !$omp end parallel do

    end if

  else   ! arbitrary cx value

    if (cy .eq. __VecZero__) then

      !$omp parallel do default(shared) if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1,NEQ
        Dy(irow) = cx*Da(irow)*Dx(irow)
      end do
      !$omp end parallel do

    else if(cy .eq. __VecOne__) then

      !$omp parallel do default(shared) if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1,NEQ
        Dy(irow) = Dy(irow) + cx*Da(irow)*Dx(irow)
      end do
      !$omp end parallel do

    else if(cy .eq. -__VecOne__) then

      !$omp parallel do default(shared) if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1,NEQ
        Dy(irow) = -Dy(irow) + cx*Da(irow)*Dx(irow)
      end do
      !$omp end parallel do

    else   ! arbitrary cy value

      !$omp parallel do default(shared) if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1,NEQ
        Dy(irow) = cy*Dy(irow) + cx*Da(irow)*Dx(irow)
      end do
      !$omp end parallel do

    end if

  end if

#endif
