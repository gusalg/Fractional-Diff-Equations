PROGRAM FDE_linear_ly_full_memo
!
! This code was written by Gustavo Salgado in Fortran 95 to 
! Solve a Single Variable Linear Caputo Fractional Differential
! Equations
!
IMPLICIT NONE
 !--------------- ----------------------------------------
 INTEGER, PARAMETER :: p_r=8, p_i=8
 REAL(KIND=p_r), DIMENSION(:,:), ALLOCATABLE :: X, ci
 REAL(KIND=p_r), DIMENSION(2,1) :: X_temp, aux
 REAL(KIND=p_r) :: a, b, y_0
 REAL(KIND=p_r) :: q, h, memo
 REAL :: time, timei, timef
 INTEGER(KIND=p_i) :: i, j, n, cont
 !--------------------------------------------------------
 !
 WRITE(*,*) " "
 WRITE(*,*) "System: ^C D^q_0 y = y0 lambda y"
 WRITE(*,*) " "
 !
 WRITE(*,*) "Derivation Order, q, in interval (0,2):"
 READ(*,*) q
 WRITE(*,*) "lambda coordinates (complex - a+ib):"
 READ(*,*) a, b
 WRITE(*,*) "Inicial condition:"
 READ(*,*) y_0
 WRITE(*,*) "Integration step size, h:"
 READ(*,*) h
 WRITE(*,*) "Number of iterations, n:"
 READ(*,*) n
 !
 CALL CPU_TIME(timei)
 ALLOCATE(X(2,n),ci(2,n))
 !
 !
  X(:,1)=(/ y_0 , 0.0_8 /)
  !--------------------------------------------------------
  ci(:,1)=(/ 1.0, 1.0 /)
  memo=-1.0_8 
  !--------------------------------------------------------
  cont=1 
  !--------------------------------------------------------
  DO i=2,n   
    !------------------------------------------------------- 
    ci(:,i)=(1-(1+q)/(i-1))*ci(:,i-1)
    memo=memo-ci(1,i)
    X_temp=s(i)
    aux(:,1)=(/ DF(X(:,i-1)) /)
    X(:,i)= aux(:,1)*h**q-X_temp(:,1)
    !-------------------------------------------------------
    !WRITE(*,*) cont
    cont=cont+1
    END DO
  !-------------------------------------------------------   
 CALL CPU_TIME(timef) 
 !-------------------------------------------------------------------
 OPEN(UNIT=100, FILE="FDE_linear_ly_full_memo.dat", STATUS="REPLACE")
 DO j=1,n
   WRITE(100,*) h*(j-1), X(1,j), X(2,j)
 END DO
 CLOSE(100)
 !-------------------------------------------------------------------
 DEALLOCATE(X,ci)
 time=timef-timei
 !--------------------------------------------------------
 WRITE(*,*) "Time ellipsed:", time, "segundos"
 WRITE(*,*) "Fractional Order:", q
 WRITE(*,*) "Parameters a=", a, "e b=", b
 WRITE(*,*) "Step h:", h
 WRITE(*,*) "Final Time:", h*cont
!-------------------------------------------------
!
CONTAINS
!-------------------------------------------------  
 FUNCTION DF(y)
  REAL(KIND=p_r), DIMENSION(2,1) :: DF
  REAL(KIND=p_r), DIMENSION(2,1), INTENT(IN) :: y
  !------------Fractional Diff Equation-------
   DF(1,1)=a*y(1,1)-b*y(2,1)     
   !------------------------------------------
   DF(2,1)=b*y(1,1)+a*y(2,1)
   !------------------------------------------
!-------------------------------------------------
  RETURN
 END FUNCTION DF
!-------------------------------------------------  
 FUNCTION s(m)
  REAL(KIND=p_r), DIMENSION(2,1) :: s
  INTEGER(KIND=p_i), INTENT(IN) :: m
  INTEGER(KIND=p_i) :: l
  s(:,1)=(/ 0, 0/)
   DO l=2,m
    s(:,1)=s(:,1)+ci(:,l)*X(:,m-l+1)
   END DO
    s(:,1)=s(:,1)+memo*X(:,1)
  RETURN
 END FUNCTION s
!-------------------------------------------------  
! 
!
END PROGRAM FDE_linear_ly_full_memo
