program MetodoNewton

    use funcao_test_sincos
    !use funcao_test_poly
    !use funcao_test1
    implicit none
    
    !integer :: deg
    !real(kind=8), allocatable :: coef(:)
    !dimension (deg+1) :: coef
    
    real(kind=8) :: p1,p2,c1,c2,n1,n2
    
    real(kind=8) :: x0,e1,e2,xmin,fmin
    integer :: maxiter
    
    maxiter = 500
    e1 = 1.0d-9
    e2 = 1.0d-5
    x0 = -2.0d0
!     
    c1 = 15.0d0
    c2 = 23.0d0
    p1 = 0.5d0
    p2 = 127.0d0
    n1 = 128.0d0
    n2 = -3.41d0
    
!     deg = 5
!     allocate(coef(deg))
!     coef = (/-1.0d0,1.0d0,2.0d0,3.0d0,4.0d0,5.0d0/)
    
    xmin = newton(x0,maxiter,e1,e2,n1,n2,p1,p2,c1,c2)
    fmin = f(xmin,n1,n2,p1,p2,c1,c2)
    
!     xmin = newton(x0,maxiter,e1,e2,coef,deg)
!     fmin = f(xmin,coef,deg)
    
!     xmin = newton(x0,maxiter,e1,e2)
!     fmin = f(xmin)
    
    print*, xmin
    print*, fmin

    !----------------------------------------------------------------------------------------------------------
    !------------------------------------------ Método de Newton -------------------------------------
    !----------------------------------------------------------------------------------------------------------
    
    contains
    
    function newton(x0,maxiter,e1,e2,n1,n2,p1,p2,c1,c2)
    !function newton(x0,maxiter,e1,e2,coef,deg)
!     function newton(x0,maxiter,e1,e2)
        implicit none

!         integer :: deg
!         real(kind=8), dimension (deg+1) :: coef
        real(kind=8) :: p1,p2,c1,c2,n1,n2
        real(kind=8) :: e1,e2
        real(kind=8) :: newton,x0,xn,fn,fnder
        integer :: maxiter,numiter
        
        xn = x0
        
        !fn = f(xn)
    
        !fn = f(xn,coef,deg)
    
        fn = f(xn,n1,n2,p1,p2,c1,c2)
    
        numiter = 0d0
    
        do while ( (abs(fn) .GE. e1) .AND. (numiter .LE. maxiter) )
            
            !fnder = fder(xn)
            !xn = xn - fn/fnder
            !fn = f(xn)
            
            !fnder=fder(xn,coef,deg)
            !xn = xn - fn/fnder
            !fn=f(xn,coef,deg)
            
            fnder = fder(xn,n1,n2,p1,p2,c1,c2)
            xn = xn - fn/fnder
            fn = f(xn,n1,n2,p1,p2,c1,c2)
        
            numiter=numiter+1
        end do
    
        print*, numiter
        newton = xn
    
    end function
    
    !----------------------------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------------------
    
end program
