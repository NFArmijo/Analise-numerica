program Romberg

    implicit none
    
    real(kind=8) :: a,b,tol,sol
    integer(kind=8) :: i, maxIter
    
    !---------------------------------- Test --------------------------------------
    
    tol = 0.000001d0
    maxIter = 500
    a = 0.0d0
    !b = 1.0d0
    b = 0.995d0
    
    sol = romb(a,b,tol,maxIter)
    
    
    !----------------------------------------------------------------------------------
    
    contains
    
    !------------------------------------------------------------------------------------
    !----------------------------------------- Método de Romberg ---------------------------
    !---------------------------------------------------------------------------------------
    
    
    function romb(a,b,tol,maxIter)
        implicit none
        
        real(kind=8) :: a,b,h,tol,romb,erro
        real(kind=8), allocatable :: R1(:),R2(:)
        integer(kind=8) :: i,j, maxIter
        
        h = b-a
        
        allocate(R1(1))
        allocate(R2(2))
        
        R1(1) = trapez(h,a,b)
        
        h = h/2
        
        R2(1) = trapez(h,a,b)
        R2(2) = R2(1) + ((R2(1)-R1(1)))/3
        
        erro = abs(R2(2)-R2(1))
        
        deallocate(R1)
        allocate(R1(2))
            
        R1 = R2
        
        i=2
        
        !do while (i .LE. maxIter .AND. erro .GT. tol*R2(i))
        do while (i .LE. maxIter .AND. erro .GT. tol)
        
            i=i+1
        
            deallocate(R2)
            allocate(R2(i))
            
            h = h/2
            R2(1) = trapez(h,a,b)
            
            do j=2,i
                R2(j) = R2(j-1) + ((R2(j-1)-R1(j-1))/((4**(real(j)-1))-1))
            end do
            
            !erro = abs(R2(i)-R2(i-1))
            erro = abs(R2(i)-R1(i-1))
            
            deallocate(R1)
            allocate(R1(i))
            
            R1 = R2
        
        end do
        
        romb = R2(i)
        
        print*, ' '
        print*, 'Aproximacao Romberg: ', R2(i)
        print*, ' '
        print*, 'Numero de iteracoes: ', i
        print*, ' '
        print*, 'Erro aproximado: ', erro
        print*, ' '
        
    end function
    
    !------------------------------------------------------------------------------------------
    !------------------------------------------ Trapezio ------------------------------------
    !------------------------------------------------------------------------------------------
    
    
     function trapez(h,a,b)
        implicit none
        
        real(kind=8) :: a,b,h,trapez,aux
        integer(kind=8) :: n,i
        
        n = int((b-a)/h)
        
        trapez = f(a)
        
        do i=1,n-1
            aux = a+real(i)*h
            trapez = trapez + 2*f(aux)
        end do
        
        trapez = (h/2)*(trapez + f(b))
        
    end

    !------------------------------------------------------------------------------------------
    !------------------------------------------ Funções Test ------------------------------------
    !------------------------------------------------------------------------------------------

     !------------------------------------------------------------------------
!     
    function f(x)
        implicit none
        
        real(kind=8) :: x,f
        
        f = 1/(1-x)
    
    end function
!     
!     !------------------------------------------------------------------------
!     
!     function f(x)
!         implicit none
!         
!         real(kind=8) :: x,f
!         
!         f = 1/(1+25*x)
!     
!     end function
!  
    !-----------------------------------------------------------------------------
!
!     function f(x)
!         implicit none
!         
!         real(kind=8) :: x,f
!         
!         f = 1
!     
!     end function
!
!     !-------------------------------------------------------------------------
!         
!     function f(x)
!         implicit none
!         
!         real(kind=8) :: x,f
!         
!         f = x
!     
!     end function
!     
!
    !-------------------------------------------------------------------------
!
!     function f(x)
!         implicit none
!         
!         real(kind=8) :: x,f
!         
!         f = x**9
!     
!     end function
    
end program
