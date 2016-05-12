program Lab2
  implicit none
  integer, parameter :: N=8, M=5, NDIM = 8
  integer            :: i, j, z
  character(10)      :: form = "(8f15.10)"
  real               :: orignMatrix(N,N), orignb(N), p(M), A(N,N), b(N) 
  real               :: TransA(N,N), Transb(N), delta
  real               :: COND, WORK(N) ! for call subroutine
  integer            :: IPVT(N)       ! 
  orignMatrix = transpose(reshape((/ -3.,  -4.,  -4.,   7.,   2.,   3.,   8.,   7.,&
                                  &   0., -15.,  -1.,   5.,  -3.,   6.,   6.,  -6.,&
                                  &  -4.,   2., -16.,   7.,   0.,   8.,  -7.,   6.,&
                                  &   0.,   8.,  -5., -11.,   1.,   0.,   4.,   5.,&
                                  &   8.,   6.,  -8.,   4.,  27.,  -7.,  -1.,   5.,&
                                  &  -4.,  -2.,   1.,   2.,  -8.,  10.,   7.,   0.,&
                                  &   0.,  -1.,   5.,   2.,  -8.,   2.,  -2.,   0.,&
                                  &   0.,  -8.,  -7.,   3.,  -7.,  -4.,  -8.,   5./), shape(orignMatrix)))
  
  orignb = (/54., -72., -33., -15., 180., -5., -14., -131./);

  p = (/1., 0.1, 0.01, 0.0001, 0.000001/)
  
  do z=1, M
    print '(a, f8.6)', "p is ", p(z)
    
    A=orignMatrix
    b=orignb

    A(1,1)=A(1,1)+p(z)
    b(1)  =b(1)  +2*p(z)
    TransA = matmul(transpose(A),A)  ! left Gauss
    Transb = matmul(transpose(A), b) ! transformation
    
    print *, "Original A:"
    print form, ((A(i,j),j=1,N),i=1,N)
    print *, "Original b:"
    print form, b       
    print *, "Transp A:"
    print form, ((TransA(i,j),j=1,N),i=1,N)
    print *, "Transp b:"
    print "(8f15.7)", Transb

  
    call DECOMP(NDIM,N,A,COND,IPVT,WORK)
    print '(/, a, es14.6, /)', "Conditional number of A: ", COND
    
    call SOLVE(NDIM,N,A,b,IPVT) !in b array subroutine write result (our x1)
    
    print * , "x from normal matrix is: "
    print form, b
    
    call DECOMP(NDIM,N,TransA,COND,IPVT,WORK)
    print '(/, a, es14.6, /)', "Conditional number of Gauss A: ", COND
    
    call SOLVE(NDIM,N,TransA, Transb, IPVT) !in Transb our result (x2 i mean)
    
    print * , "x from Gauss left transpose  matrix is: "
    print form, Transb
    
    delta = norm2(b-Transb)/norm2(b)
    print '(/, a, f14.6,/)', "Delta - ||(x1-x2)||/||x1|| = ", delta
    read(*,*) !getch analog in fortran
  end do
!  print form, ((orignMatrix(i,j),j=1,N),i=1,N)
!  print form, orignb

end program Lab2
