program MetodoSecante

    use funcao_test_raiz
    !use funcao_test_poly
    !use funcao_test_sincos
    implicit none
    
    real(kind=8) :: m !raiz desejada
    
    !integer :: deg
    !real(kind=8), allocatable :: coef(:)
    !dimension (deg+1) :: coef
    
    !real(kind=8) :: p1,p2,c1,c2,n1,n2
    
    real(kind=8) :: x0,x1,e1,e2,xmin,fmin
    integer :: maxiter
    
    maxiter = 500
    e1 = 1.0d-9
    e2 = 1.0d-5
    x0 = -1.5d0
    x1 = 2.0d0
    
    m = 23.0d0 !raiz desejada
    
    !deg = 5
    !allocate(coef(deg))
    !coef = (/-1.0d0,1.0d0,2.0d0,3.0d0,4.0d0,5.0d0/)
    
    !c1 = 15.0d0
    !c2 = 23.0d0
    !p1 = 0.5d0
    !p2 = 127.0d0
    !n1 = 128.0d0
    !n2 = -3.41d0
    
    xmin = secante(m,x0,x1,maxiter,e1,e2)
    fmin = f(xmin,m)
    
    !xmin = secante(x0,x1,maxiter,e1,e2,coef,deg)
    !fmin = f(xmin,coef,deg)
    
    !xmin = secante(x0,x1,maxiter,e1,e2,n1,n2,p1,p2,c1,c2)
    !fmin = f(xmin,n1,n2,p1,p2,c1,c2)

    print*,"Pontos iniciais: "
    print*,"x0 = ",x0
    print*,"x1 = ",x1
    print*,"Solucao: ",xmin
    print*,"Valor da funcao: ",fmin

    !----------------------------------------------------------------------------------------------------------
    !------------------------------------------ Método da Secante -------------------------------------
    !----------------------------------------------------------------------------------------------------------
    
    contains
    
    function secante(m,x0,x1,maxiter,e1,e2)
    !function secante(x0,x1,maxiter,e1,e2,coef,deg)
    !function secante(x0,x1,maxiter,e1,e2,n1,n2,p1,p2,c1,c2)
    
        implicit none

        real(kind=8) :: m
        
        !integer :: deg
        !real(kind=8), dimension (deg+1) :: coef
        
        !real(kind=8) :: p1,p2,c1,c2,n1,n2
        
        real(kind=8) :: e1,e2,div,xaux
        real(kind=8) :: secante,x0,x1,xn,xm,fn,fm
        integer :: maxiter,numiter
        
        xm = x0
        xn = x1
        
        fm = f(xm,m)
        fn = f(xn,m)
        
        !fn = f(xn,coef,deg)
        !fm = f(xm,coef,deg)
    
        !fn = f(xn,n1,n2,p1,p2,c1,c2)
        !fm = f(xm,n1,n2,p1,p2,c1,c2)
    
        numiter = 0
    
        do while ( (abs(fn) .GE. e1) .AND. (numiter .LE. maxiter) )
            
            div = 1.0d0/(fn-fm)
            xaux = xn
            xn = xn-fn*(xn-xm)*div
            xm = xaux
            fn = f(xn,m)
            fm = f(xm,m)
            
            !div = 1.d0/(fn-fm)
            !xaux = xn
            !xn = xn-fn*(xn-xm)*div
            !xm = xaux
            !fn = f(xn,coef,deg)
            !fm = f(xm,coef,deg)
            
            !div = 1.d0/(fn-fm)
            !xaux = xn
            !xn = xn-fn*(xn-xm)*div
            !xm = xaux
            !fn = f(xn,n1,n2,p1,p2,c1,c2)
            !fm = f(xm,n1,n2,p1,p2,c1,c2)
        
            print*,"erro: ",abs(fn)
        
            numiter=numiter+1
        end do
    
        print*,"Numero de iteracoes: ",numiter
        secante = xn
    
    end function
    
    !----------------------------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------------------
    
end program
