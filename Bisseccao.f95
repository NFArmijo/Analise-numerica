
program MetodoBisseccao

    use funcao_test_sincos
    !use funcao_test_poly
    !use funcao_test1
    implicit none
    
    !integer :: deg
    !real(kind=8), allocatable :: coef(:)
    !dimension (deg+1) :: coef
    
    real(kind=8) :: p1,p2,c1,c2,n1,n2
    
    real(kind=8) :: a,b,e1,e2,xmin,fmin
    integer :: maxiter
    
    maxiter = 500
    e1 = 1.0d-9
    e2 = 1.0d-5
    a = -1.0d0
    b = 1.0d0
    
    c1 = 15.0d0
    c2 = 23.0d0
    p1 = 0.5d0
    p2 = 127.0d0
    n1 = 128.0d0
    n2 = -3.41d0
    
    !deg = 5
    !allocate(coef(deg))
    !coef = (/-1.0,1.0,2.0,3.0,4.0,5.0/)
    
    xmin = bisseccao(a,b,maxiter,e1,e2,n1,n2,p1,p2,c1,c2)
    fmin = f(xmin,n1,n2,p1,p2,c1,c2)
    
    !xmin = bisseccao(a,b,maxiter,e1,e2,coef,deg)
    !fmin = f(xmin,coef,deg)
    
    !xmin = bisseccao(a,b,maxiter,e1,e2)
    !fmin = f(xmin)
    
    print*, xmin
    print*, fmin

    !----------------------------------------------------------------------------------------------------------
    !------------------------------------------ Método de Bissecção -------------------------------------
    !----------------------------------------------------------------------------------------------------------
    
    contains
    
    function bisseccao(a,b,maxiter,e1,e2,n1,n2,p1,p2,c1,c2)
    !function bisseccao(a,b,maxiter,e1,e2,coef,deg)
    !function bisseccao(a,b,maxiter,e1,e2)
        implicit none

        !integer :: deg
        !real(kind=8), dimension (deg+1) :: coef
        real(kind=8) :: p1,p2,c1,c2,n1,n2
        real(kind=8) :: p,a,b,e1,e2,fb,fp
        real*8 :: bisseccao
        integer :: maxiter,numiter
    
        p = (a+b)*0.5d0
        
        !fb = f(b)
        !fp = f(p)
    
        !fb = f(b,coef,deg)
        !fp = f(p,coef,deg)
    
        fb = f(b,n1,n2,p1,p2,c1,c2)
        fp = f(p,n1,n2,p1,p2,c1,c2)
    
        numiter = 0
    
        do while ( ((abs(fp) .GE. e1) .OR. (abs(b-a) .GE. e2)) .AND. (numiter .LE. maxiter) )
            ! não precisa-se abs(b-a) pois a < b.
        
            if (fb*fp .GT. 0.d0) then
                b = p
            else
                a = p
            endif
            
            p = (a+b)*0.5d0
            
            !fb = f(b)
            !fp = f(p)
            
            !fb=f(b,coef,deg)
            !fp=f(p,coef,deg)
            
            fb = f(b,n1,n2,p1,p2,c1,c2)
            fp = f(p,n1,n2,p1,p2,c1,c2)
        
            numiter=numiter+1
        end do
    
        print*, numiter
        bisseccao = p
    
    end function
    
    !----------------------------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------------------
    
end program

