!-*- mode: f90; -*-

  ! Perform the multiplication
  if (cx .ne. __MatOne__) then
    
    if (cy .eq. __VecZero__) then
      
      !$omp parallel do default(shared) private(ia,icol,ivar,jvar,Ddtmp) &
      !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1,NEQ
        Ddtmp = __MatZero__
        do ia = Kld(irow), Kld(irow+1)-1
          icol = Kcol(ia)
          do ivar = 1,__NVAR__
            do jvar = 1,__NVAR__
              Ddtmp(ivar) = Ddtmp(ivar) + Dx(__NVAR__*(icol-1)+jvar)&
                          * Da(__NVAR__*__NVAR__*(ia-1)+__NVAR__*(jvar-1)+ivar)
            end do
          end do
        end do
        do ivar = 1,__NVAR__
          Dy(__NVAR__*(irow-1)+ivar) = cx*Ddtmp(ivar)
        end do
      end do
      !$omp end parallel do
      
    else if(cy .eq. __VecOne__) then
      
      !$omp parallel do default(shared) private(ia,icol,ivar,jvar,Ddtmp) &
      !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1,NEQ
        Ddtmp = __MatZero__
        do ia = Kld(irow), Kld(irow+1)-1
          icol = Kcol(ia)
          do ivar = 1,__NVAR__
            do jvar = 1,__NVAR__
              Ddtmp(ivar) = Ddtmp(ivar) + Dx(__NVAR__*(icol-1)+jvar)&
                          * Da(__NVAR__*__NVAR__*(ia-1)+__NVAR__*(jvar-1)+ivar)
            end do
          end do
        end do
        do ivar = 1,__NVAR__
          Dy(__NVAR__*(irow-1)+ivar) = Dy(__NVAR__*(irow-1)+ivar) + cx*Ddtmp(ivar)
        end do
      end do
      !$omp end parallel do
      
    else   ! arbitrary cy value
      
      !$omp parallel do default(shared) private(ia,icol,ivar,jvar,Ddtmp) &
      !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1,NEQ
        Ddtmp = __MatZero__
        do ia = Kld(irow), Kld(irow+1)-1
          icol = Kcol(ia)
          do ivar = 1,__NVAR__
            do jvar = 1,__NVAR__
              Ddtmp(ivar) = Ddtmp(ivar) + Dx(__NVAR__*(icol-1)+jvar)&
                          * Da(__NVAR__*__NVAR__*(ia-1)+__NVAR__*(jvar-1)+ivar)
            end do
          end do
        end do
        do ivar = 1,__NVAR__
          Dy(__NVAR__*(irow-1)+ivar) = cy*Dy(__NVAR__*(irow-1)+ivar) + cx*Ddtmp(ivar)
        end do
      end do
      !$omp end parallel do
      
    end if
    
  else   ! arbitrary cx value
    
    if(cy .eq. __VecZero__) then
      
      !$omp parallel do default(shared) private(ia,icol,ivar,jvar,Ddtmp) &
      !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1,NEQ
        Ddtmp = __MatZero__
        do ia = Kld(irow), Kld(irow+1)-1
          icol = Kcol(ia)
          do ivar = 1,__NVAR__
            do jvar = 1,__NVAR__
              Ddtmp(ivar) = Ddtmp(ivar) + Dx(__NVAR__*(icol-1)+jvar)&
                          * Da(__NVAR__*__NVAR__*(ia-1)+__NVAR__*(jvar-1)+ivar)
            end do
          end do
        end do
        do ivar = 1,__NVAR__
          Dy(__NVAR__*(irow-1)+ivar) = Ddtmp(ivar)
        end do
      end do
      !$omp end parallel do
      
    else if(cy .eq. __VecOne__) then
      
      !$omp parallel do default(shared) private(ia,icol,ivar,jvar,Ddtmp) &
      !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1,NEQ
        Ddtmp = __MatZero__
        do ia = Kld(irow), Kld(irow+1)-1
          icol = Kcol(ia)
          do ivar = 1,__NVAR__
            do jvar = 1,__NVAR__
              Ddtmp(ivar) = Ddtmp(ivar) + Dx(__NVAR__*(icol-1)+jvar)&
                          * Da(__NVAR__*__NVAR__*(ia-1)+__NVAR__*(jvar-1)+ivar)
            end do
          end do
        end do
        do ivar = 1,__NVAR__
          Dy(__NVAR__*(irow-1)+ivar) = Dy(__NVAR__*(irow-1)+ivar) + Ddtmp(ivar)
        end do
      end do
      !$omp end parallel do
      
    else   ! arbitrary cy value
      
      !$omp parallel do default(shared) private(ia,icol,ivar,jvar,Ddtmp) &
      !$omp if(NEQ > rperfconfig%NEQMIN_OMP)
      do irow = 1,NEQ
        Ddtmp = __MatZero__
        do ia = Kld(irow), Kld(irow+1)-1
          icol = Kcol(ia)
          do ivar = 1,__NVAR__
            do jvar = 1,__NVAR__
              Ddtmp(ivar) = Ddtmp(ivar) + Dx(__NVAR__*(icol-1)+jvar)&
                          * Da(__NVAR__*__NVAR__*(ia-1)+__NVAR__*(jvar-1)+ivar)
            end do
          end do
        end do
        do ivar = 1,__NVAR__
          Dy(__NVAR__*(irow-1)+ivar) = cy*Dy(__NVAR__*(irow-1)+ivar) + Ddtmp(ivar)
        end do
      end do
      !$omp end parallel do
      
    end if
    
  end if
