program Rayleigh_Ritz

    implicit none
    
    real(kind=8), allocatable :: a(:),b(:),d(:) !Matriz e vetor independente
    real(kind=8), allocatable :: c(:) !Solução do sistema
    real(kind=8), allocatable :: u(:) !Solução aproximada
    real(kind=8), allocatable :: ureal(:) !Solução real da EDO
    real(kind=8), allocatable :: udif(:) !Diferença das soluções
    real(kind=8) :: erro,pi
    real(kind=8) :: aux,aux0,aux1,aux2
    integer(kind=8) :: n,i,j,k
    
    !-------------------------- Exercício 1 ---------------------------------
    
!     do j=4,8
!         n = 2**j-1
!         
!         allocate(a(n))      !Diagonal da matriz
!         allocate(b(n-1))    !Subdiagonal da matriz
!         allocate(d(n))      !Vetor independente
!         allocate(c(n))      !Solução do sistema
!         allocate(u(n))      !Aproximação da função solução da EDO
!         allocate(ureal(n))  !Função solução da EDO
!         allocate(udif(n))   !Diferença das soluções
!     
!         do i=1,n-1
!             d(i) = (1/(real(n)+1))*(((real(i)+1)**2)/2+real(i)**2+((real(i)-1)**2)/2-(real(i)-1)*real(i)-real(i)*(real(i)+1))
!             a(i) = 2*(real(n)+1)
!             b(i) = -real(n)-1
!         end do
!         a(n) = 2*(real(n)+1)
!         d(n) = (1/(real(n)+1))*(((real(n)+1)**2)/2+real(n)**2+((real(n)-1)**2)/2-(real(n)-1)*real(n)-real(n)*(real(n)+1))
!     
!         c = thomas3diag(a,b,d,n)
!     
!         do i=1,n
!             u(i) = c(i)
!             ureal(i) = 0.5d0*(real(i)/(real(n)+1))*(1-(real(i)/(real(n)+1)))
!             udif(i) = u(i)-ureal(i)
!         end do
!     
!         erro = norm_inf(udif,n)
!         print*, 'Error of the approximation of the sol of the EDO: ', erro
!         print*, ' '
!         !print*, 'Approximated solution of the EDO: ', u
!         print*, ' '
!         
!         deallocate(a)
!         deallocate(b)
!         deallocate(d)
!         deallocate(c)
!         deallocate(u)
!         deallocate(ureal)
!         deallocate(udif)
! 
!     end do
    
    !--------------------------- Exercício 2 --------------------------
        
!     !open(1, file = 'erroExer2v2.txt', status = 'new')
!         
!     do j=4,8
! 
!         !n = j
!         n = 2**j-1
!         
!         allocate(a(n))      !Diagonal da matriz
!         allocate(b(n-1))    !Subdiagonal da matriz
!         allocate(d(n))      !Vetor independente
!         allocate(c(n))      !Solução do sistema
!         allocate(u(10*n+1))      !Aproximação da função solução da EDO
!         allocate(ureal(10*n+1))  !Função solução da EDO
!         allocate(udif(10*n+1))   !Diferença das soluções
!     
!         aux = real(n)+1
!         
!         do i=1,n-1
!             aux0 = (real(i)-1)/aux
!             aux1 = real(i)/aux
!             aux2 = (real(i)+1)/aux
!             
!             d(i) = aux*( -(aux0**4+6*(aux1**4)+aux2**4)+2*(aux0**3)+4*(aux1**3)*(2+aux0+aux2) &
!             +2*(aux2**3)-aux0**2-2*(aux1**2)*(1+3*(aux0+aux2))-aux2**2+2*aux1*(aux0+aux2) )
!             
!             a(i) = 2*aux
!             b(i) = -aux
!             
!         end do
!         
!         aux0 = (real(n)-1)/aux
!         aux1 = real(n)/aux
!         aux2 = 1.0d0
!         
!         a(n) = 2*aux
!         
!         d(n) = aux*( -(aux0**4+6*(aux1**4)+aux2**4)+2*(aux0**3)+4*(aux1**3)*(2+aux0+aux2) &
!         +2*(aux2**3)-aux0**2-2*(aux1**2)*(1+3*(aux0+aux2))-aux2**2+2*aux1*(aux0+aux2) )
!     
!         c = thomas3diag(a,b,d,n)
!         
!         do i=1,10*n
!             k = int(((real(i)-1)*(real(n)+1))/(10*real(n)))+1
!             if (k .NE. 1 .AND. k .NE. n+1) then
!                 u(i) = c(k)*aux*(((real(i)-1)/(10*real(n)))-((real(k)-1)/aux)) + c(k-1)*aux*((real(k)/aux) &
!                 -((real(i)-1)/(10*real(n))))
!             elseif (k .EQ. 1) then
!                 u(i) = c(k)*aux*((real(i)-1)/(10*real(n)))
!             else
!                 u(i) = c(k-1)*aux*(1-((real(i)-1)/(10*real(n))))
!             end if
!             ureal(i) = (((real(i)-1)/(10*real(n)))**2)*((1-((real(i)-1)/(10*real(n))))**2)
!             udif(i) = u(i)-ureal(i)
!         end do
! 
!         
!         u(10*n+1) = 0.0d0
!         ureal(10*n+1) = 0.0d0
!         udif(10*n+1) = u(10*n+1)-ureal(10*n+1)
!     
!         erro = norm_inf(udif,10*n+1)
!         
!         print*, 'Error of the approximation of the sol of the EDO: ', erro
!         print*, ' '
!         print*, 'Approximated solution of the EDO: ', u
!         print*, ' '
!         
!         !write(1,*) erro
!         
!         deallocate(a)
!         deallocate(b)
!         deallocate(d)
!         deallocate(c)
!         deallocate(u)
!         deallocate(ureal)
!         deallocate(udif)
!         
!     end do
!     
!     !close(1)
    
    !----------------------------- Exercício 3 -----------------------------------
    
    !open(1, file = 'erroExer3.txt', status = 'new')
    
    pi = 3.14159265d0
    
    do j=4,8
        
        n = 2**j-1
        
        allocate(a(n))      !Diagonal da matriz
        allocate(b(n-1))    !Subdiagonal da matriz
        allocate(d(n))      !Vetor independente
        allocate(c(n))      !Solução do sistema. Coeficientes da combinação linear splines
        allocate(u(10*n+1))      !Aproximação da função solução da EDO
        allocate(ureal(10*n+1))  !Função solução da EDO
        allocate(udif(10*n+1))   !Diferença das soluções
    
        aux = real(n)+1
        
        do i=1,n-1
            aux0 = pi*((real(i)-1)/aux)
            aux1 = pi*(real(i)/aux)
            aux2 = pi*((real(i)+1)/aux)
            
            d(i) = 2*aux*( 2*sin(aux1)-sin(aux0)-sin(aux2)+(cos(aux1))*(aux0-2*aux1+aux2) )
            
            a(i) = (2*aux) + ((2*(pi**2))/(3*aux))
            b(i) = -aux+((pi**2)/(6*aux))
            
        end do
        
        aux0 = pi*((real(n)-1)/aux)
        aux1 = pi*(real(n)/aux)
        aux2 = pi
        
        a(n) = (2*aux) + ((2*(pi**2))/(3*aux))
        d(n) = 2*aux*( 2*sin(aux1)-sin(aux0)-sin(aux2)+(cos(aux1))*(aux0-2*aux1+aux2) )
    
        c = thomas3diag(a,b,d,n)
        
        do i=1,10*n
            k = int(((real(i)-1)*(real(n)+1))/(10*real(n)))+1
            if (k .NE. 1 .AND. k .NE. n+1) then
                u(i) = c(k)*aux*(((real(i)-1)/(10*real(n)))-((real(k)-1)/aux)) + c(k-1)*aux*((real(k)/aux) &
                -((real(i)-1)/(10*real(n))))
            elseif (k .EQ. 1) then
                u(i) = c(k)*aux*((real(i)-1)/(10*real(n)))
            else
                u(i) = c(k-1)*aux*(1-((real(i)-1)/(10*real(n))))
            end if
            ureal(i) = sin(pi*((real(i)-1)/(10*real(n))))
            udif(i) = u(i)-ureal(i)
        end do
        
        u(10*n+1) = 0.0d0
        ureal(10*n+1) = 0.0d0
        udif(10*n+1) = u(10*n+1)-ureal(10*n+1)
    
        erro = norm_inf(udif,10*n+1)
        
        print*, 'Error of the approximation of the sol of the EDO: ', erro
        print*, ' '
        print*, 'Approximated solution of the EDO: ', u
        print*, ' '
        
        !write(1,*) erro
        
        deallocate(a)
        deallocate(b)
        deallocate(d)
        deallocate(c)
        deallocate(u)
        deallocate(ureal)
        deallocate(udif)
        
        !close(1)
        
    end do
    
    !------------------------------------------------------------------------
    
    contains
    
    subroutine ldl3diag(a,b,n,d,l)
        implicit none
        integer(kind=8), intent(in) :: n
        real(kind=8), dimension(n), intent(in) :: a
        real(kind=8), dimension(n-1), intent(in) :: b
        real(kind=8), dimension(n), intent(out) :: d
        real(kind=8), dimension(n-1), intent(out) :: l
        
        integer :: i
        
        d(1) = a(1)
        l(1) = b(1)/d(1)
        
        do i=2,n-1
            d(i) = a(i)-b(i-1)*l(i-1)
            l(i) = b(i)/d(i)
        end do
        
        d(n) = a(n)-b(n-1)*l(n-1)
        
    end subroutine ldl3diag
    
    !----------------------------------------------------------------------------------------------------------
    !------------------------------------------ Método de Thomas ---------------------------------------------
    !----------------------------------------------------------------------------------------------------------
    
    function thomas3diag(a,b,c,n)
        implicit none
        integer(kind=8) :: n
        real(kind=8), dimension(n) :: a,c
        real(kind=8), dimension(n-1) :: b
        integer(kind=8) :: i
        real(kind=8), dimension(n) :: d
        real(kind=8), dimension(n-1) :: l
        real(kind=8), dimension(n) :: thomas3diag
        real(kind=8), dimension(n) :: aux
        
        call ldl3diag(a,b,n,d,l)
    
        aux(1) = c(1)
        do i=1,n
            aux(i) = c(i)-l(i-1)*aux(i-1)
        end do
        
        do i=1,n
            aux(i) = aux(i)/d(i)
        end do
        
        thomas3diag(n) = aux(n)
        do i=n,1,-1
            thomas3diag(i) = aux(i)-l(i)*thomas3diag(i+1)
        end do
        
        !print*, ' '
        !print*, 'Matrix A, diagonal: ', a
        !print*, ' '
        !print*, 'Matrix A, subdiagonal : ', b
        !print*, ' '
        !print*, 'Independent vector b : ', c
        !print*, ' '
        !print*, 'L matrix : ', l
        !print*, ' '
        !print*, 'D matrix : ', d
        !print*, ' '
        !print*, 'The solution of the system Ac=d is: ', thomas3diag !solSist
        !print*, ' '
        
    end function
    
    !------------------------------------------------------------------------------------------
    !------------------------------------------Norma infinito ------------------------------------
    !------------------------------------------------------------------------------------------
    
     function norm_inf(M,d)
        implicit none
        
        integer(kind=8) :: d,i
        real(kind=8), dimension(d) :: M
        real(kind=8) :: norm_inf,aux
        
        norm_inf = abs(M(1))
        
        do i=2,d
            aux= abs(M(i))
            if (aux > norm_inf) then
                norm_inf = aux
            end if
        end do
        
    end
    
end program
