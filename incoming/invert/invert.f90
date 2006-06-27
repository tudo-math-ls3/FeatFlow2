MODULE inv
  IMPLICIT NONE

  INTEGER, PARAMETER :: DP=8

CONTAINS

!<subroutine>

SUBROUTINE invert(Da,Df,Dx,ndim,ipar)
  
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

  REAL(DP), DIMENSION(ndim), INTENT(IN) :: Df
  INTEGER, INTENT(IN) :: ndim,ipar

  !</input>


  !<inputoutput>
  
  REAL(DP), DIMENSION(ndim,ndim), INTENT(INOUT) :: Da

  !</inputoutput>


  !<output>
  
  REAL(DP), DIMENSION(ndim), INTENT(OUT) :: Dx
  
  !</output>

!</subroutine>


  ! local variables
  REAL(DP), DIMENSION(ndim,ndim) :: Db
  REAL(DP), DIMENSION(ndim) :: Dpiv
  INTEGER, DIMENSION(ndim) :: Kindx,Kindy
  
  REAL(DP) :: dpivot,daux
  INTEGER :: idim1,idim2,idim3,ix,iy,indx,indy,info
    
  SELECT CASE (ipar)
  CASE (0)
     ! Perform factorization of matrix Da
     
     ! Initialization
     Kindx=0;  Kindy=0
     
     DO idim1=1,ndim
        
        ! Determine pivotal element
        dpivot=0
        
        DO iy=1,ndim
           IF (Kindy(iy) /= 0) CYCLE
           
           DO ix=1,ndim
              IF (Kindx(ix) /= 0) CYCLE
              
              IF (ABS(Da(ix,iy)) .LE. ABS(dpivot)) CYCLE
              dpivot=Da(ix,iy);  indx=ix;  indy=iy
           END DO
        END DO
        
        ! Return if pivotal element is zero
        IF (ABS(dpivot) .LE. 0._DP) RETURN
        
        Kindx(indx)=indy;  Kindy(indy)=indx;  Da(indx,indy)=1._DP/dpivot
        
        DO idim2=1,ndim
           IF (idim2 == indy) CYCLE 
           Da(1:indx-1,   idim2)=Da(1:indx-1,   idim2)-Da(1:indx-1,   indy)*Da(indx,idim2)/dpivot
           Da(indx+1:ndim,idim2)=Da(indx+1:ndim,idim2)-Da(indx+1:ndim,indy)*Da(indx,idim2)/dpivot
        END DO
        
        DO ix=1,ndim
           IF (ix /= indx) Da(ix,indy)=Da(ix,indy)/dpivot
        END DO
        
        DO iy=1,ndim
           IF (iy /= indy) Da(indx,iy)=-Da(indx,iy)/dpivot
        END DO
     END DO
     
     DO ix=1,ndim
        IF (Kindx(ix) == ix) CYCLE
        
        DO iy=1,ndim
           IF (Kindx(iy) == ix) EXIT
        END DO
        
        DO idim1=1,ndim
           daux=Da(ix,idim1);  Da(ix,idim1)=Da(iy,idim1);  Da(iy,idim1)=daux
        END DO
        
        Kindx(iy)=Kindx(ix);  Kindx(ix)=ix
     END DO
        
     DO ix=1,ndim
        IF (Kindy(ix) == ix) CYCLE
        
        DO iy=1,ndim
           IF (Kindy(iy) == ix) EXIT
        END DO
        
        DO idim1=1,ndim
           daux=Da(idim1,ix);  Da(idim1,ix)=Da(idim1,iy);  Da(idim1,iy)=daux
        END DO
        
        Kindy(iy)=Kindy(ix);  Kindy(ix)=ix
     END DO
         

  CASE (1)    
     ! Perform inversion of Da to solve the Da * Dx = Df
     DO idim1=1,ndim
        Dx(idim1)=0
        DO idim2=1,ndim
           Dx(idim1)=Dx(idim1)+Da(idim1,idim2)*Df(idim2)
        END DO
     END DO

  CASE (2)
     ! Solve the dense linear system Ax=f calling LAPACK routine
     
     SELECT CASE(ndim)
     CASE (2)
        ! Explicit formula for 2x2 system
        Db(1,1)= Da(2,2)
        Db(2,1)=-Da(2,1)
        Db(1,2)=-Da(1,2)
        Db(2,2)= Da(1,1)
        daux=Da(1,1)*Da(2,2)-Da(1,2)*Da(2,1)
        Dx=MATMUL(Db,Df)/daux
        
     CASE (3)
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
        Dx=MATMUL(Db,Df)/daux

     CASE (4)
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
        Dx=MATMUL(Db,Df)/daux
                 
     CASE DEFAULT
        ! Use LAPACK routine for general NxN system, where N>4
        Dpiv=0; Dx=Df
        CALL DGESV(ndim,1,Da,ndim,Dpiv,Dx,ndim,info)

     END SELECT
  END SELECT
END SUBROUTINE invert

END MODULE inv
