program lab3

implicit none



call rkf(0.02)
print *," "
call RungePower4(0.01, 0.02)
print *," "
call RungePower4(0.005,0.02)

contains

subroutine rkf(Tprint)
   integer :: NEQN=2 !equation count
   real    :: Y(2) = (/0.,1./) ! first value 
   real    :: T = 0., TOUT=0. ! t= independ value. tout = point where search solution
   real    :: RELERR = 0.1E-03, ABSERR=0.0 
   real    :: WORK(27), Tprint, Tfinal=0.4
   integer :: IWORK(5), IFLAG=1, i, endloop

   endloop = tfinal/tprint
   print *, "Solution with RKF45: "
   do i = 0, endloop-1
      TOUT= TOUT+tprint

      call RKF45(fun, NEQN, Y, T, TOUT, RELERR, ABSERR, IFLAG, WORK, IWORK)
      print '(a, f7.4, 2x, a, f10.6, 2x, a, f10.6)', 'T= ', T, 'Y1= ', Y(1), 'Y2= ', Y(2)
      select case (IFLAG)
  
         case(3)
         print *, "function RKF45 autochange RELERR"
        
         case(4)
         print *, "count of steps is very large"

         case(5)
         print *, "program have to change ABSERR"
         ABSERR=0.1E-07

         case(6)
         print *, "our error accuracy reachless and RELERR being change"
         RELERR=RELERR*10.0
         IFLAG=2
    
         case(7)
         print *, "Not effective"
         IFLAG=2
    
         case(8)
         print *, "wrong input parameters"

         case default
        
         end select
   end do
end subroutine rkf

subroutine RungePower4(h, Tprint)
   real     :: k1(2), k2(2), k3(2), k4(2)
   real     :: Y(2) = (/0.,1./), tempY(2), yp(2) ! first value end temporary value
   real     :: h, t=0.0
   real     :: Tprint
   integer  :: i, endloop, pcounter, counter=0
   
   y(1)=0
   y(2)=1
   t=0
   counter=0
   print '(a, f7.4)', "Runge-Kutta with power 4 and step is: ", h

   endloop=0.4/h
   pcounter=Tprint/h
   print '(a, f7.4, 2x, a, f10.6, 2x, a, f10.6)', 'T= ', T, 'Y1= ', Y(1), 'Y2= ', Y(2)
   do i= 0, endloop
   tempY=Y
   call fun(t, tempY, yp)
   k1 = h*yp
   call fun(t + h/3, Y + k1/3, yp)
   k2=h*yp
   call fun(t + 2*h/3, Y - k1/3 + k2, yp)
   k3=h*yp
   call fun(t+h,Y+k1-k2+k3, yp)
   k4=h*yp
   Y=Y+(k1+3*k2+3*k3+k4)/8
      if(pcounter .eq. counter) then
         counter = 0
         print '(a, f7.4, 2x, a, f10.6, 2x, a, f10.6)', 'T= ', T, 'Y1= ', Y(1), 'Y2= ', Y(2)
      end if
   t=h+t
   counter = counter + 1 
   end do

end subroutine RungePower4


subroutine fun(t, y, yp)
  real :: t, y(2), yp(2)
  yp(1)= -310*y(1)-3000*y(2)+1/(10*t*t+1)
  yp(2)=y(1)+exp(-2*t)
  return
end subroutine fun

end program lab3
