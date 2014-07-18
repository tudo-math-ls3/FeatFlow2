#ifdef __NVAR__

  ! Perform the multiplication
  if(cx .eq. __MatOne__) then

    if(cy .eq. __VecZero__) then

      !$omp parallel do private(ia,icol,ivar,Ddtmp) default(shared) &
      !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1, NEQ
        Ddtmp = __MatZero__
        do ia = Kld(irow), Kld(irow+1)-1
          icol = Kcol(ia)
          do ivar = 1, __NVAR__
            Ddtmp(ivar) = Ddtmp(ivar) + Da(ia)*Dx(__NVAR__*(icol-1)+ivar)
          end do
        end do
        do ivar = 1, __NVAR__
          Dy(__NVAR__*(irow-1)+ivar) = Ddtmp(ivar)
        end do
      end do
      !$omp end parallel do

    else if(cy .eq. __VecOne__) then

      !$omp parallel do private(ia,icol,ivar,Ddtmp) default(shared) &
      !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1, NEQ
        Ddtmp = __MatZero__
        do ia = Kld(irow), Kld(irow+1)-1
          icol = Kcol(ia)
          do ivar = 1, __NVAR__
            Ddtmp(ivar) = Ddtmp(ivar) + Da(ia)*Dx(__NVAR__*(icol-1)+ivar)
          end do
        end do
        do ivar = 1, __NVAR__
          Dy(__NVAR__*(irow-1)+ivar) = Dy(__NVAR__*(irow-1)+ivar) + Ddtmp(ivar)
        end do
      end do
      !$omp end parallel do

    else if(cy .eq. -__VecOne__) then

      !$omp parallel do private(ia,icol,ivar,Ddtmp) default(shared) &
      !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1, NEQ
        Ddtmp = __MatZero__
        do ia = Kld(irow), Kld(irow+1)-1
          icol = Kcol(ia)
          do ivar = 1, __NVAR__
            Ddtmp(ivar) = Ddtmp(ivar) + Da(ia)*Dx(__NVAR__*(icol-1)+ivar)
          end do
        end do
        do ivar = 1, __NVAR__
          Dy(__NVAR__*(irow-1)+ivar) = -Dy(__NVAR__*(irow-1)+ivar) + Ddtmp(ivar)
        end do
      end do
      !$omp end parallel do

    else   ! arbitrary cy value

      !$omp parallel do private(ia,icol,ivar,Ddtmp) default(shared) &
      !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1, NEQ
        Ddtmp = __MatZero__
        do ia = Kld(irow), Kld(irow+1)-1
          icol = Kcol(ia)
          do ivar = 1, __NVAR__
            Ddtmp(ivar) = Ddtmp(ivar) + Da(ia)*Dx(__NVAR__*(icol-1)+ivar)
          end do
        end do
        do ivar = 1, __NVAR__
          Dy(__NVAR__*(irow-1)+ivar) = cy*Dy(__NVAR__*(irow-1)+ivar) + Ddtmp(ivar)
        end do
      end do
      !$omp end parallel do

    end if

  else if(cx .eq. -__MatOne__) then

    if(cy .eq. __VecZero__) then

      !$omp parallel do private(ia,icol,ivar,Ddtmp) default(shared) &
      !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1, NEQ
        Ddtmp = __MatZero__
        do ia = Kld(irow), Kld(irow+1)-1
          icol = Kcol(ia)
          do ivar = 1, __NVAR__
            Ddtmp(ivar) = Ddtmp(ivar) + Da(ia)*Dx(__NVAR__*(icol-1)+ivar)
          end do
        end do
        do ivar = 1, __NVAR__
          Dy(__NVAR__*(irow-1)+ivar) = -Ddtmp(ivar)
        end do
      end do
      !$omp end parallel do

    else if(cy .eq. __VecOne__) then

      !$omp parallel do private(ia,icol,ivar,Ddtmp) default(shared) &
      !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1, NEQ
        Ddtmp = __MatZero__
        do ia = Kld(irow), Kld(irow+1)-1
          icol = Kcol(ia)
          do ivar = 1, __NVAR__
            Ddtmp(ivar) = Ddtmp(ivar) + Da(ia)*Dx(__NVAR__*(icol-1)+ivar)
          end do
        end do
        do ivar = 1, __NVAR__
          Dy(__NVAR__*(irow-1)+ivar) = Dy(__NVAR__*(irow-1)+ivar) - Ddtmp(ivar)
        end do
      end do
      !$omp end parallel do

    else if(cy .eq. -__VecOne__) then

      !$omp parallel do private(ia,icol,ivar,Ddtmp) default(shared) &
      !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1, NEQ
        Ddtmp = __MatZero__
        do ia = Kld(irow), Kld(irow+1)-1
          icol = Kcol(ia)
          do ivar = 1, __NVAR__
            Ddtmp(ivar) = Ddtmp(ivar) + Da(ia)*Dx(__NVAR__*(icol-1)+ivar)
          end do
        end do
        do ivar = 1, __NVAR__
          Dy(__NVAR__*(irow-1)+ivar) = -Dy(__NVAR__*(irow-1)+ivar) - Ddtmp(ivar)
        end do
      end do
      !$omp end parallel do

    else   ! arbitrary cy value

      !$omp parallel do private(ia,icol,ivar,Ddtmp) default(shared) &
      !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1, NEQ
        Ddtmp = __MatZero__
        do ia = Kld(irow), Kld(irow+1)-1
          icol = Kcol(ia)
          do ivar = 1, __NVAR__
            Ddtmp(ivar) = Ddtmp(ivar) + Da(ia)*Dx(__NVAR__*(icol-1)+ivar)
          end do
        end do
        do ivar = 1, __NVAR__
          Dy(__NVAR__*(irow-1)+ivar) = cy*Dy(__NVAR__*(irow-1)+ivar) - Ddtmp(ivar)
        end do
      end do
      !$omp end parallel do

    end if

  else   ! arbitrary cx value
    
    if(cy .eq. __VecZero__) then

      !$omp parallel do private(ia,icol,ivar,Ddtmp) default(shared) &
      !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1, NEQ
        Ddtmp = __MatZero__
        do ia = Kld(irow), Kld(irow+1)-1
          icol = Kcol(ia)
          do ivar = 1, __NVAR__
            Ddtmp(ivar) = Ddtmp(ivar) + Da(ia)*Dx(__NVAR__*(icol-1)+ivar)
          end do
        end do
        do ivar = 1, __NVAR__
          Dy(__NVAR__*(irow-1)+ivar) = cx*Ddtmp(ivar)
        end do
      end do
      !$omp end parallel do

    else if(cy .eq. __VecOne__) then

      !$omp parallel do private(ia,icol,ivar,Ddtmp) default(shared) &
      !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1, NEQ
        Ddtmp = __MatZero__
        do ia = Kld(irow), Kld(irow+1)-1
          icol = Kcol(ia)
          do ivar = 1, __NVAR__
            Ddtmp(ivar) = Ddtmp(ivar) + Da(ia)*Dx(__NVAR__*(icol-1)+ivar)
          end do
        end do
        do ivar = 1, __NVAR__
          Dy(__NVAR__*(irow-1)+ivar) = Dy(__NVAR__*(irow-1)+ivar) + cx*Ddtmp(ivar)
        end do
      end do
      !$omp end parallel do

    else if(cy .eq. -__VecOne__) then

      !$omp parallel do private(ia,icol,ivar,Ddtmp) default(shared) &
      !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1, NEQ
        Ddtmp = __MatZero__
        do ia = Kld(irow), Kld(irow+1)-1
          icol = Kcol(ia)
          do ivar = 1, __NVAR__
            Ddtmp(ivar) = Ddtmp(ivar) + Da(ia)*Dx(__NVAR__*(icol-1)+ivar)
          end do
        end do
        do ivar = 1, __NVAR__
          Dy(__NVAR__*(irow-1)+ivar) = -Dy(__NVAR__*(irow-1)+ivar) + cx*Ddtmp(ivar)
        end do
      end do
      !$omp end parallel do

    else   ! arbitrary cy value

      !$omp parallel do private(ia,icol,ivar,Ddtmp) default(shared) &
      !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1, NEQ
        Ddtmp = __MatZero__
        do ia = Kld(irow), Kld(irow+1)-1
          icol = Kcol(ia)
          do ivar = 1, __NVAR__
            Ddtmp(ivar) = Ddtmp(ivar) + Da(ia)*Dx(__NVAR__*(icol-1)+ivar)
          end do
        end do
        do ivar = 1, __NVAR__
          Dy(__NVAR__*(irow-1)+ivar) = cy*Dy(__NVAR__*(irow-1)+ivar) + cx*Ddtmp(ivar)
        end do
      end do
      !$omp end parallel do

    end if

  end if

#else

  ! Perform the multiplication
  if(cx .eq. __MatOne__) then

    if(cy .eq. __VecZero__) then

      !$omp parallel do private(ia,dtmp) default(shared) &
      !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1, NEQ
        dtmp = __MatZero__
        do ia = Kld(irow), Kld(irow+1)-1
          dtmp = dtmp + Da(ia)*Dx(Kcol(ia))
        end do
        Dy(irow) = dtmp
      end do
      !$omp end parallel do

    else if(cy .eq. __VecOne__) then

      !$omp parallel do private(ia,dtmp) default(shared) &
      !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1, NEQ
        dtmp = __MatZero__
        do ia = Kld(irow), Kld(irow+1)-1
          dtmp = dtmp + Da(ia)*Dx(Kcol(ia))
        end do
        Dy(irow) = Dy(irow) + dtmp
      end do
      !$omp end parallel do

    else if(cy .eq. -__VecOne__) then

      !$omp parallel do private(ia,dtmp) default(shared) &
      !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1, NEQ
        dtmp = __MatZero__
        do ia = Kld(irow), Kld(irow+1)-1
          dtmp = dtmp + Da(ia)*Dx(Kcol(ia))
        end do
        Dy(irow) = -Dy(irow) + dtmp
      end do
      !$omp end parallel do

    else   ! arbitrary cy value

      !$omp parallel do private(ia,dtmp) default(shared) &
      !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1, NEQ
        dtmp = __MatZero__
        do ia = Kld(irow), Kld(irow+1)-1
          dtmp = dtmp + Da(ia)*Dx(Kcol(ia))
        end do
        Dy(irow) = cy*Dy(irow) + dtmp
      end do
      !$omp end parallel do

    end if

  else if(cx .eq. -__MatOne__) then

    if(cy .eq. __VecZero__) then

      !$omp parallel do private(ia,dtmp) default(shared) &
      !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1, NEQ
        dtmp = __MatZero__
        do ia = Kld(irow), Kld(irow+1)-1
          dtmp = dtmp + Da(ia)*Dx(Kcol(ia))
        end do
        Dy(irow) = -dtmp
      end do
      !$omp end parallel do

    else if(cy .eq. __VecOne__) then

      !$omp parallel do private(ia,dtmp) default(shared) &
      !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1, NEQ
        dtmp = __MatZero__
        do ia = Kld(irow), Kld(irow+1)-1
          dtmp = dtmp + Da(ia)*Dx(Kcol(ia))
        end do
        Dy(irow) = Dy(irow) - dtmp
      end do
      !$omp end parallel do

    else if(cy .eq. -__VecOne__) then

      !$omp parallel do private(ia,dtmp) default(shared) &
      !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1, NEQ
        dtmp = __MatZero__
        do ia = Kld(irow), Kld(irow+1)-1
          dtmp = dtmp + Da(ia)*Dx(Kcol(ia))
        end do
        Dy(irow) = -Dy(irow) - dtmp
      end do
      !$omp end parallel do

    else   ! arbitrary cy value

      !$omp parallel do private(ia,dtmp) default(shared) &
      !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1, NEQ
        dtmp = __MatZero__
        do ia = Kld(irow), Kld(irow+1)-1
          dtmp = dtmp + Da(ia)*Dx(Kcol(ia))
        end do
        Dy(irow) = cy*Dy(irow) - dtmp
      end do
      !$omp end parallel do

    end if

  else   ! arbitrary cx value

    if(cy .eq. __VecZero__) then

      !$omp parallel do private(ia,dtmp) default(shared) &
      !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1, NEQ
        dtmp = __MatZero__
        do ia = Kld(irow), Kld(irow+1)-1
          dtmp = dtmp + Da(ia)*Dx(Kcol(ia))
        end do
        Dy(irow) = cx*dtmp
      end do
      !$omp end parallel do

    else if(cy .eq. __VecOne__) then

      !$omp parallel do private(ia,dtmp) default(shared) &
      !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1, NEQ
        dtmp = __MatZero__
        do ia = Kld(irow), Kld(irow+1)-1
          dtmp = dtmp + Da(ia)*Dx(Kcol(ia))
        end do
        Dy(irow) = Dy(irow) + cx*dtmp
      end do
      !$omp end parallel do

    else if(cy .eq. -__VecOne__) then

      !$omp parallel do private(ia,dtmp) default(shared) &
      !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1, NEQ
        dtmp = __MatZero__
        do ia = Kld(irow), Kld(irow+1)-1
          dtmp = dtmp + Da(ia)*Dx(Kcol(ia))
        end do
        Dy(irow) = -Dy(irow) + cx*dtmp
      end do
      !$omp end parallel do

    else   ! arbitrary cy value

      !$omp parallel do private(ia,dtmp) default(shared) &
      !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1, NEQ
        dtmp = __MatZero__
        do ia = Kld(irow), Kld(irow+1)-1
          dtmp = dtmp + Da(ia)*Dx(Kcol(ia))
        end do
        Dy(irow) = cy*Dy(irow) + cx*dtmp
      end do
      !$omp end parallel do

    end if   

  end if

#endif
