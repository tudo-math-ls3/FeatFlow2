module inv
  implicit none

  integer, parameter :: DP=8

contains

!<subroutine>

subroutine invert(Da,Df,Dx,ndim,ipar)
  
  !<description>
  
  ! This subroutine performs the direct inversion of a NxN system.
  ! 
  ! If the parameter ipar=0, then only factorization of matrix A is performed. 
  ! For ipar=1, the vector x is calculated using the factorized matrix A. 
  ! For ipar=2, LAPACK routine DGESV is used to solve the dense linear system Ax=f. 
  ! In addition, for NDIM=2,3,4 explicit formulas are employed to replace the
  ! more expensive LAPACK routine.
  ! 

  !</description>


  !<input>

  real(DP), dimension(ndim), intent(IN) :: Df
  integer, intent(IN) :: ndim,ipar

  !</input>


  !<inputoutput>
  
  real(DP), dimension(ndim,ndim), intent(INOUT) :: Da

  !</inputoutput>


  !<output>
  
  real(DP), dimension(ndim), intent(OUT) :: Dx
  
  !</output>

!</subroutine>


  ! local variables
  real(DP), dimension(ndim,ndim) :: Db
  real(DP), dimension(ndim) :: Dpiv
  integer, dimension(ndim) :: Kindx,Kindy
  
  real(DP) :: dpivot,daux
  integer :: idim1,idim2,idim3,ix,iy,indx,indy,info
    
  select case (ipar)
  case (0)
     ! Perform factorization of matrix Da
     
     ! Initialization
     Kindx=0;  Kindy=0
     
     do idim1=1,ndim
        
        ! Determine pivotal element
        dpivot=0
        
        do iy=1,ndim
           if (Kindy(iy) /= 0) cycle
           
           do ix=1,ndim
              if (Kindx(ix) /= 0) cycle
              
              if (abs(Da(ix,iy)) .le. abs(dpivot)) cycle
              dpivot=Da(ix,iy);  indx=ix;  indy=iy
           end do
        end do
        
        ! Return if pivotal element is zero
        if (abs(dpivot) .le. 0._DP) return
        
        Kindx(indx)=indy;  Kindy(indy)=indx;  Da(indx,indy)=1._DP/dpivot
        
        do idim2=1,ndim
           if (idim2 == indy) cycle 
           Da(1:indx-1,   idim2)=Da(1:indx-1,   idim2)-Da(1:indx-1,   indy)*Da(indx,idim2)/dpivot
           Da(indx+1:ndim,idim2)=Da(indx+1:ndim,idim2)-Da(indx+1:ndim,indy)*Da(indx,idim2)/dpivot
        end do
        
        do ix=1,ndim
           if (ix /= indx) Da(ix,indy)=Da(ix,indy)/dpivot
        end do
        
        do iy=1,ndim
           if (iy /= indy) Da(indx,iy)=-Da(indx,iy)/dpivot
        end do
     end do
     
     do ix=1,ndim
        if (Kindx(ix) == ix) cycle
        
        do iy=1,ndim
           if (Kindx(iy) == ix) exit
        end do
        
        do idim1=1,ndim
           daux=Da(ix,idim1);  Da(ix,idim1)=Da(iy,idim1);  Da(iy,idim1)=daux
        end do
        
        Kindx(iy)=Kindx(ix);  Kindx(ix)=ix
     end do
        
     do ix=1,ndim
        if (Kindy(ix) == ix) cycle
        
        do iy=1,ndim
           if (Kindy(iy) == ix) exit
        end do
        
        do idim1=1,ndim
           daux=Da(idim1,ix);  Da(idim1,ix)=Da(idim1,iy);  Da(idim1,iy)=daux
        end do
        
        Kindy(iy)=Kindy(ix);  Kindy(ix)=ix
     end do
         

  case (1)    
     ! Perform inversion of Da to solve the Da * Dx = Df
     do idim1=1,ndim
        Dx(idim1)=0
        do idim2=1,ndim
           Dx(idim1)=Dx(idim1)+Da(idim1,idim2)*Df(idim2)
        end do
     end do

  case (2)
     ! Solve the dense linear system Ax=f calling LAPACK routine
     
     select case(ndim)
     case (2)
        ! Explicit formula for 2x2 system
        Db(1,1)= Da(2,2)
        Db(2,1)=-Da(2,1)
        Db(1,2)=-Da(1,2)
        Db(2,2)= Da(1,1)
        daux=Da(1,1)*Da(2,2)-Da(1,2)*Da(2,1)
        Dx=matmul(Db,Df)/daux
        
     case (3)
        ! Explicit formula for 3x3 system
        Db(1,1)=Da(2,2)*Da(3,3)-Da(2,3)*Da(3,2)
        Db(2,1)=Da(2,3)*Da(3,1)-Da(2,1)*Da(3,3)
        Db(3,1)=Da(2,1)*Da(3,2)-Da(2,2)*Da(3,1)
        Db(1,2)=Da(1,3)*Da(3,2)-Da(1,2)*Da(3,3)
        Db(2,2)=Da(1,1)*Da(3,3)-Da(1,3)*Da(3,1)
        Db(3,2)=Da(1,2)*Da(3,1)-Da(1,1)*Da(3,2)
        Db(1,3)=Da(1,2)*Da(2,3)-Da(1,3)*Da(2,2)
        Db(2,3)=Da(1,3)*Da(2,1)-Da(1,1)*Da(2,3)
        Db(3,3)=Da(1,1)*Da(2,2)-Da(1,2)*Da(2,1)
        daux=Da(1,1)*Da(2,2)*Da(3,3)+Da(2,1)*Da(3,2)*Da(1,3)+ &
             Da(3,1)*Da(1,2)*Da(2,3)-Da(1,1)*Da(3,2)*Da(2,3)- &
             Da(3,1)*Da(2,2)*Da(1,3)-Da(2,1)*Da(1,2)*Da(3,3)
        Dx=matmul(Db,Df)/daux

     case (4)
        ! Explicit formula for 4x4 system
        Db(1,1)=Da(2,2)*Da(3,3)*Da(4,4)+Da(2,3)*Da(3,4)*Da(4,2)+Da(2,4)*Da(3,2)*Da(4,3)- &
             Da(2,2)*Da(3,4)*Da(4,3)-Da(2,3)*Da(3,2)*Da(4,4)-Da(2,4)*Da(3,3)*Da(4,2)
        Db(2,1)=Da(2,1)*Da(3,4)*Da(4,3)+Da(2,3)*Da(3,1)*Da(4,4)+Da(2,4)*Da(3,3)*Da(4,1)- &
             Da(2,1)*Da(3,3)*Da(4,4)-Da(2,3)*Da(3,4)*Da(4,1)-Da(2,4)*Da(3,1)*Da(4,3)
        Db(3,1)=Da(2,1)*Da(3,2)*Da(4,4)+Da(2,2)*Da(3,4)*Da(4,1)+Da(2,4)*Da(3,1)*Da(4,2)- &
             Da(2,1)*Da(3,4)*Da(4,2)-Da(2,2)*Da(3,1)*Da(4,4)-Da(2,4)*Da(3,2)*Da(4,1)
        Db(4,1)=Da(2,1)*Da(3,3)*Da(4,2)+Da(2,2)*Da(3,1)*Da(4,3)+Da(2,3)*Da(3,2)*Da(4,1)- &
             Da(2,1)*Da(3,2)*Da(4,3)-Da(2,2)*Da(3,3)*Da(4,1)-Da(2,3)*Da(3,1)*Da(4,2)
        Db(1,2)=Da(1,2)*Da(3,4)*Da(4,3)+Da(1,3)*Da(3,2)*Da(4,4)+Da(1,4)*Da(3,3)*Da(4,2)- &
             Da(1,2)*Da(3,3)*Da(4,4)-Da(1,3)*Da(3,4)*Da(4,2)-Da(1,4)*Da(3,2)*Da(4,3)
        Db(2,2)=Da(1,1)*Da(3,3)*Da(4,4)+Da(1,3)*Da(3,4)*Da(4,1)+Da(1,4)*Da(3,1)*Da(4,3)- &
             Da(1,1)*Da(3,4)*Da(4,3)-Da(1,3)*Da(3,1)*Da(4,4)-Da(1,4)*Da(3,3)*Da(4,1)
        Db(3,2)=Da(1,1)*Da(3,4)*Da(4,2)+Da(1,2)*Da(3,1)*Da(4,4)+Da(1,4)*Da(3,2)*Da(4,1)- &
             Da(1,1)*Da(3,2)*Da(4,4)-Da(1,2)*Da(3,4)*Da(4,1)-Da(1,4)*Da(3,1)*Da(4,2)
        Db(4,2)=Da(1,1)*Da(3,2)*Da(4,3)+Da(1,2)*Da(3,3)*Da(4,1)+Da(1,3)*Da(3,1)*Da(4,2)- &
             Da(1,1)*Da(3,3)*Da(4,2)-Da(1,2)*Da(3,1)*Da(4,3)-Da(1,3)*Da(3,2)*Da(4,1)
        Db(1,3)=Da(1,2)*Da(2,3)*Da(4,4)+Da(1,3)*Da(2,4)*Da(4,2)+Da(1,4)*Da(2,2)*Da(4,3)- &
             Da(1,2)*Da(2,4)*Da(4,3)-Da(1,3)*Da(2,2)*Da(4,4)-Da(1,4)*Da(2,3)*Da(4,2)
        Db(2,3)=Da(1,1)*Da(2,4)*Da(4,3)+Da(1,3)*Da(2,1)*Da(4,4)+Da(1,4)*Da(2,3)*Da(4,1)- &
             Da(1,1)*Da(2,3)*Da(4,4)-Da(1,3)*Da(2,4)*Da(4,1)-Da(1,4)*Da(2,1)*Da(4,3)
        Db(3,3)=Da(1,1)*Da(2,2)*Da(4,4)+Da(1,2)*Da(2,4)*Da(4,1)+Da(1,4)*Da(2,1)*Da(4,2)- &
             Da(1,1)*Da(2,4)*Da(4,2)-Da(1,2)*Da(2,1)*Da(4,4)-Da(1,4)*Da(2,2)*Da(4,1)
        Db(4,3)=Da(1,1)*Da(2,3)*Da(4,2)+Da(1,2)*Da(2,1)*Da(4,3)+Da(1,3)*Da(2,2)*Da(4,1)- &
             Da(1,1)*Da(2,2)*Da(4,3)-Da(1,2)*Da(2,3)*Da(4,1)-Da(1,3)*Da(2,1)*Da(4,2)
        Db(1,4)=Da(1,2)*Da(2,4)*Da(3,3)+Da(1,3)*Da(2,2)*Da(3,4)+Da(1,4)*Da(2,3)*Da(3,2)- &
             Da(1,2)*Da(2,3)*Da(3,4)-Da(1,3)*Da(2,4)*Da(3,2)-Da(1,4)*Da(2,2)*Da(3,3)
        Db(2,4)=Da(1,1)*Da(2,3)*Da(3,4)+Da(1,3)*Da(2,4)*Da(3,1)+Da(1,4)*Da(2,1)*Da(3,3)- &
             Da(1,1)*Da(2,4)*Da(3,3)-Da(1,3)*Da(2,1)*Da(3,4)-Da(1,4)*Da(2,3)*Da(3,1)
        Db(3,4)=Da(1,1)*Da(2,4)*Da(3,2)+Da(1,2)*Da(2,1)*Da(3,4)+Da(1,4)*Da(2,2)*Da(3,1)- &
             Da(1,1)*Da(2,2)*Da(3,4)-Da(1,2)*Da(2,4)*Da(3,1)-Da(1,4)*Da(2,1)*Da(3,2)
        Db(4,4)=Da(1,1)*Da(2,2)*Da(3,3)+Da(1,2)*Da(2,3)*Da(3,1)+Da(1,3)*Da(2,1)*Da(3,2)- &
             Da(1,1)*Da(2,3)*Da(3,2)-Da(1,2)*Da(2,1)*Da(3,3)-Da(1,3)*Da(2,2)*Da(3,1)
        daux=Da(1,1)*Da(2,2)*Da(3,3)*Da(4,4)+Da(1,1)*Da(2,3)*Da(3,4)*Da(4,2)+Da(1,1)*Da(2,4)*Da(3,2)*Da(4,3)+ &
             Da(1,2)*Da(2,1)*Da(3,4)*Da(4,3)+Da(1,2)*Da(2,3)*Da(3,1)*Da(4,4)+Da(1,2)*Da(2,4)*Da(3,3)*Da(4,1)+ &
             Da(1,3)*Da(2,1)*Da(3,2)*Da(4,4)+Da(1,3)*Da(2,2)*Da(3,4)*Da(4,1)+Da(1,3)*Da(2,4)*Da(3,1)*Da(4,2)+ &
             Da(1,4)*Da(2,1)*Da(3,3)*Da(4,2)+Da(1,4)*Da(2,2)*Da(3,1)*Da(4,3)+Da(1,4)*Da(2,3)*Da(3,2)*Da(4,1)- &
             Da(1,1)*Da(2,2)*Da(3,4)*Da(4,3)-Da(1,1)*Da(2,3)*Da(3,2)*Da(4,4)-Da(1,1)*Da(2,4)*Da(3,3)*Da(4,2)- &
             Da(1,2)*Da(2,1)*Da(3,3)*Da(4,4)-Da(1,2)*Da(2,3)*Da(3,4)*Da(4,1)-Da(1,2)*Da(2,4)*Da(3,1)*Da(4,3)- &
             Da(1,3)*Da(2,1)*Da(3,4)*Da(4,2)-Da(1,3)*Da(2,2)*Da(3,1)*Da(4,4)-Da(1,3)*Da(2,4)*Da(3,2)*Da(4,1)- &
             Da(1,4)*Da(2,1)*Da(3,2)*Da(4,3)-Da(1,4)*Da(2,2)*Da(3,3)*Da(4,1)-Da(1,4)*Da(2,3)*Da(3,1)*Da(4,2)
        Dx=matmul(Db,Df)/daux
                 
     case DEFAULT
        ! Use LAPACK routine for general NxN system, where N>4
        Dpiv=0; Dx=Df
        call DGESV(ndim,1,Da,ndim,Dpiv,Dx,ndim,info)

     end select
  end select
end subroutine invert

end module inv
