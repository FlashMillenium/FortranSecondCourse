real function fun(x)
real x
fun= tan(x)/x
end function fun

program Lab1 !find data with lagrange and spline function and integral value

use Spline !use Spline.f95 as module

EXTERNAL fun

INTEGER NOFUN
INTEGER, PARAMETER :: N = 7, M=19
REAL fun, A, B, RELERR, ABSERR,RESULT,ERREST, FLAG ! var for quanc8
DATA A/1.0/, B/2.0/, RELERR/1.E-06/,ABSERR/0.0/
REAL :: XAR(N) = (/ -1., -0.96, -0.86, -0.79,  0.22,  0.5,  0.93/)   ! our data
REAL :: FAR(N) = (/ -1., -0.151, 0.894, 0.986, 0.895, 0.5, -0.306 /) ! x and f(x)
INTEGER            :: Out = 0, i
CHARACTER(10)      :: form
REAL Bs(N), Cs(N), Ds(N), SEVAL !var for SPLINE and SEVAL
REAL XiRES(M), LRES(M), SRES(M) ! output array

call SPLIN(N, XAR, FAR, Bs, Cs, Ds)

!create output data from spline and lagrange function
do i=1,M
   XiRES(i) = -1 + 0.1*i ! xk= -1 + 0.1*k (k=1,2,...,19)
   SRES(i) = SEVAL(N,XiRES(i),XAR,FAR,Bs,Cs,Ds)
   LRES(i) = LgrF(XiRES(i))
end do

!integral calculation

CALL QUANC8(fun, A, B, ABSERR, RELERR, RESULT, ERREST, NOFUN, FLAG)
print *, 'Integral value from Quanc8:', RESULT
print *, 'ERREST: ', ERREST, ' NOFUN: ', NOFUN, ' FLAG: ', FLAG


open(file="output.txt", encoding="UTF-8", newunit=Out)
     write(form,'(a, i0, a)') '(', M, 'f12.8)'
     write(Out, *) "Xk: "
     write(Out, form) XiRES
     write(Out, *) "Spline f(x): "
     write(Out, form) SRES
     write(Out, *) "Lagrange f(x): "
     write(Out, form) LRES
close(Out)

contains
! Lagrange function
real function LgrF(x)
real x, numerator, denominator
integer j,i
i=0
LgrF=0
do i=1, N
   numerator=1
   denominator=1
   do j=1, N
      if(i /= j) then
        numerator=numerator * (x - XAR(j))
        ! (x-x1)(x-x2)...(x-x(i-1))(x-x(i+1))...(x-xN)
        denominator=denominator * (XAR(i) - XAR(j))
        ! (xi-x1)(xi-x2)...(xi-x(i-1))(xi-x(i+1))...(xi-xN)
      end if
   end do
   LgrF=LgrF+FAR(i)*numerator/denominator !
end do
end function LgrF

end program Lab1
