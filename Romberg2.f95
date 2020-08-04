program Romberg2

    implicit none
    
    real(kind=8) :: a,b,tol,sol
    integer(kind=8) :: i, maxIter
    
    !---------------------------------- Test --------------------------------------
    
    tol = 0.000001d0
    maxIter = 10000
    
    a = 0.0d0 ! Integral I_1
    b = 1.0d0
    
!     a = 0.0d0  ! Integral I_2
!     b = 0.995d0

!     a = -1.0d0  ! Integral I_3
!     b = 1.0d0
     
!     a = -5.0d0  ! Integral I_4
!     b = 5.0d0

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
        
        R1(1) = (h/2)*(f(a)+f(b))
        
        h = h/2
        
        i=2
        
        R2(1) = trapez(h,R1(1),a,b,i-1)
        R2(2) = R2(1) + ((R2(1)-R1(1)))/3
        
        erro = abs(R2(2)-R2(1))
        
        deallocate(R1)
        allocate(R1(2))
            
        R1 = R2
        
        !do while (i .LE. maxIter .AND. erro .GT. tol*R2(i)) !CRITERIO 2
        do while (i .LE. maxIter .AND. erro .GT. tol) !CRITERIO 1
        
            i=i+1
        
            deallocate(R2)
            allocate(R2(i))
            
            h = h/2
            R2(1) = trapez(h,R1(1),a,b,i-1)
            
            do j=2,i
                R2(j) = R2(j-1) + ((R2(j-1)-R1(j-1))/((4**(real(j)-1))-1))
            end do
            
            !erro = abs(R2(i)-R2(i-1)) !CRITERIO 2
            erro = abs(R2(i)-R1(i-1)) !CRITERIO 1
            
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
    
    
     function trapez(h,th,a,b,i)
        implicit none
        
        real(kind=8) :: a,b,h,trapez,aux,th
        integer(kind=8) :: i,j
        
        trapez = th/2
        
        do j=1,2**(i-1)
            aux = a+(2*real(j)-1)*h
            trapez = trapez + h*f(aux)
        end do
        
    end

    !------------------------------------------------------------------------------------------
    !------------------------------------------ Funções Test ------------------------------------
    !------------------------------------------------------------------------------------------

    !--------------------- Integral I_1 ------------------------------------------------
!
!     function f(x)
!         implicit none
!         
!         real(kind=8) :: x,f
!         
!         f = x**2
!     
!     end function
    
     !--------------------- Integral I_2 ---------------------------------------------------
!     
!     function f(x)
!         implicit none
!         
!         real(kind=8) :: x,f
!         
!         f = 1/(1-x)
!     
!     end function
!     
!     !--------------------- Integral I_3 ---------------------------------------------------
!     
!     function f(x)
!         implicit none
!         
!         real(kind=8) :: x,f
!         
!         f = 1/(1+25*(x**2))
!     
!     end function
!  
    !----------------------- Integral I_4 ------------------------------------------------------
!
!     function f(x)
!         implicit none
!         
!         real(kind=8) :: x,f
!         
!         f = 1/(1+(x**2))
!     
!     end function
!
    !--------------------------------------------------------------------------------
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

    
end program
