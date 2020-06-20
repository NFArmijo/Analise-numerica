program Rayleigh_Ritz

    implicit none
    
    real(kind=8), allocatable :: a(:),b(:),d(:) !Matriz e vetor independente
    real(kind=8), allocatable :: c(:) !Solução do sistema
    real(kind=8), allocatable :: u(:) !Solução aproximada
    real(kind=8), allocatable :: ureal(:) !Solução real da EDO
    real(kind=8), allocatable :: udif(:) !Diferença das soluções
    real(kind=8) :: erro
    integer :: n,i,j
    
    !----------- Exercício 1 ---------------------
    
    do j=4,8
        n = 2**j-1
        
        allocate(a(n))      !Diagonal da matriz
        allocate(b(n-1))    !Subdiagonal da matriz
        allocate(d(n))      !Vetor independente
        allocate(c(n))      !Solução do sistema
        allocate(u(n))      !Aproximação da função solução da EDO
        allocate(ureal(n))  !Função solução da EDO
        allocate(udif(n))   !Diferença das soluções
    
        do i=1,n-1
            d(i) = (1/(real(n)+1))*(((real(i)+1)**2)/2+real(i)**2+((real(i)-1)**2)/2-(real(i)-1)*real(i)-real(i)*(real(i)+1))
            a(i) = 2*(real(n)+1)
            b(i) = -real(n)-1
        end do
        a(n) = 2*(real(n)+1)
        d(n) = (1/(real(n)+1))*(((real(n)+1)**2)/2+real(n)**2+((real(n)-1)**2)/2-(real(n)-1)*real(n)-real(n)*(real(n)+1))
    
        c = thomas3diag(a,b,d,n)
    
        do i=1,n
            u(i) = c(i)
            ureal(i) = 0.5d0*(real(i)/(real(n)+1))*(1-(real(i)/(real(n)+1)))
            udif(i) = u(i)-ureal(i)
        end do
    
        erro = norm_inf(udif,n)
        print*, 'Error of the approximation of the sol of the EDO: ', erro
        print*, ' '
        print*, 'Approximated solution of the EDO: ', u
        print*, ' '
        print*, ' '
        
        deallocate(a)
        deallocate(b)
        deallocate(d)
        deallocate(c)
        deallocate(u)
        deallocate(ureal)
        deallocate(udif)
        
    end do
    
    !----------- Exercício 2 ---------------------
    
    !n = 4
    !allocate(a(n))
    !allocate(c(n))
    !allocate(b(n-1))
    !llocate(solSist(n))

    !solSist = thomas3diag(a,b,c,n)
    
    !----------- Exercício 3 ---------------------
    
    !n = 4
    !allocate(a(n))
    !allocate(c(n))
    !allocate(b(n-1))
    !allocate(solSist(n))
    
    !solSist = thomas3diag(a,b,c,n)

    !----------------------------------------------------------------------------------------------------------
    !------------------------------------------ Decomposição LDL Tridiagonal-----------------------------------
    !----------------------------------------------------------------------------------------------------------
    
    contains
    
    subroutine ldl3diag(a,b,n,d,l)
        implicit none
        integer, intent(in) :: n
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
        integer :: n
        real(kind=8), dimension(n) :: a,c
        real(kind=8), dimension(n-1) :: b
        integer :: i
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
        !print*, 'Matrix A, subdiagonal : ', b
        !print*, ' '
        !print*, 'Independent vector b : ', c
        !print*, ' '
        !print*, 'L matrix : ', l
        !print*, 'D matrix : ', d
        !print*, ' '
        print*, 'The solution of the system Ac=d is: ', thomas3diag !solSist
        print*, ' '
        
    end function
    
    !------------------------------------------------------------------------------------------
    !------------------------------------------Norma infinito ------------------------------------
    !------------------------------------------------------------------------------------------
    
     function norm_inf(M,d)
        implicit none
        
        integer :: d,i
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
