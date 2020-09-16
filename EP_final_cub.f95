program EP_Final_cub

    implicit none
    
    real(kind=8), allocatable :: a(:),b(:),c(:),d(:),ind(:) !Matriz e vetor independente
    real(kind=8), allocatable :: coef(:) !Solução do sistema
    real(kind=8), allocatable :: u(:) !Solução aproximada
    real(kind=8), allocatable :: ureal(:) !Solução real da EDO
    real(kind=8), allocatable :: udif(:) !Diferença das soluções
    real(kind=8) :: erro, x0, xn,tol, h, xi, y, zero, um, aux, aux2, aux10
    integer(kind=8) :: n,i,j,k,maxIter,j1,j2,j0,j10
    
    !------------------------------------------------------------------------------------
    !----------------------------------------- Exercício 1 -----------------------------------
    !---------------------------------------------------------------------------------------
        
!     open(1, file = 'sol_EP_final_cub.txt', status = 'new')
    
    tol = 0.0000001d0
    maxIter = 1000
    
    j1 = 1
    j2 = 2
    j0 = 0
    j10 = 10
    
    zero = 0
    um = 1
    aux2 = 2
    aux10 = 10
    
    do j=4,8

        n = 2**j-1
        
        allocate(a(n+2))      !Diagonal da matriz
        allocate(b(n+1))    !Subdiagonal da matriz
        allocate(c(n))    !Subdiagonal da matriz
        allocate(d(n-1))    !Subdiagonal da matriz
        allocate(ind(n+2))  !Vetor independiente do sistema linear
        allocate(coef(n+1))      !Solução do sistema
        allocate(u(10*n+1))      !Aproximação da função solução da EDO
        allocate(ureal(10*n+1))  !Função solução da EDO
        allocate(udif(10*n+1))   !Diferença das soluções
    
        aux = real(n)+1
        h = um/aux
        
        do i=j0,n-j2
            xi = real(i)*h
            
            x0 = -funmin(zero,-(xi-aux2*h))
            xn = funmin(um,xi+aux2*h)
            ind(i+j1) = romb(x0,xn,tol,h,xi,n,maxIter,int(1))
            a(i+j1) = romb(x0,xn,tol,h,xi,n,maxIter,int(2))
            
            x0 = -funmin(zero,-(xi-h))
            xn = funmin(um,xi+aux2*h)
            b(i+j1) = romb(x0,xn,tol,h,xi,n,maxIter,int(3))
            
            x0 = -funmin(zero,-xi)
            xn = funmin(um,xi+aux2*h)
            c(i+j1) = romb(x0,xn,tol,h,xi,n,maxIter,int(4))
            
            x0 = -funmin(zero,-(xi+h))
            xn = funmin(um,xi+aux2*h)
            d(i+j1) = romb(x0,xn,tol,h,xi,n,maxIter,int(5))
            
        end do
        
        xi = (real(n)-um)*h
        x0 = -funmin(zero,-(xi-aux2*h))
        xn = funmin(um,xi+aux2*h)
        ind(n) = romb(x0,xn,tol,h,xi,n,maxIter,int(1)) !i = n-1
        a(n) = romb(x0,xn,tol,h,xi,n,maxIter,int(2))
        
        x0 = -funmin(zero,-(xi-h))
        xn = funmin(um,xi+aux2*h)
        b(n) = romb(x0,xn,tol,h,xi,n,maxIter,int(3))
        
        x0 = -funmin(zero,-xi)
        xn = funmin(um,xi+aux2*h)
        c(n) = romb(x0,xn,tol,h,xi,n,maxIter,int(4))
        
        xi = real(n)*h
        x0 = -funmin(zero,-(xi-aux2*h))
        xn = funmin(um,xi+aux2*h)
        ind(n+j1) = romb(x0,xn,tol,h,xi,n,maxIter,int(1)) !i = n
        a(n+j1) = romb(x0,xn,tol,h,xi,n,maxIter,int(2))
        
        x0 = -funmin(zero,-(xi-h))
        xn = funmin(um,xi+aux2*h)
        b(n+j1) = romb(x0,xn,tol,h,xi,n,maxIter,int(3))
        
        xi = (real(n)+um)*h
        x0 = -funmin(zero,-(xi-aux2*h))
        xn = funmin(um,xi+aux2*h)
        ind(n+j2) = romb(x0,xn,tol,h,xi,n,maxIter,int(1)) !i = n+1
        a(n+j2) = romb(x0,xn,tol,h,xi,n,maxIter,int(2))
        
        coef = solLu7diag(a,b,c,d,ind,n+2)
        
        do i=j1,j10*n
            y = (real(i)-um)/(aux10*real(n))
            k = int(y/h,8)+j1
            xi = real(k)*h
            if (k .NE. j1 .AND. k .NE. n+j1) then
                u(i) = coef(k-j1)*phi(y,h,xi-aux2*h,n) + coef(k)*phi(y,h,xi-h,n) + coef(k+j1)*phi(y,h,xi,n) &
                + coef(k+j2)*phi(y,h,xi+h,n)
            elseif (k .EQ. j1) then
                u(i) = coef(k)*phi(y,h,xi-h,n) + coef(k+j1)*phi(y,h,xi,n) + coef(k+j2)*phi(y,h,xi+h,n)
            else
                u(i) = coef(k-j2)*phi(y,h,xi-aux2*h,n) + coef(k-j1)*phi(y,h,xi-h,n) + coef(k)*phi(y,h,xi,n)
            endif
            ureal(i) = (y**2)*(((um-y))**2)
            udif(i) = u(i)-ureal(i)
        end do

        u(j10*n+j1) = zero
        ureal(j10*n+j1) = zero
        udif(j10*n+j1) = u(j10*n+j1)-ureal(j10*n+j1)
    
        erro = norm_inf(udif,j10*n+j1)
        
        print*, ' '
        print*, 'Error of the approximation of the sol of the EDO: ', erro
        print*, ' '
        
!         write(1,*) u

        deallocate(a)
        deallocate(b)
        deallocate(c)
        deallocate(d)
        deallocate(ind)
        deallocate(coef)
        deallocate(u)
        deallocate(ureal)
        deallocate(udif)
        
    end do
    
!     close(1)
    
    contains
    
    !------------------------------------------------------------------------------------
    !----------------------------------------- Método de Romberg ---------------------------
    !---------------------------------------------------------------------------------------
    
    
    function romb(a,b,tol,h,xi,n,maxIter,tipo)
        implicit none
        
        real(kind=8) :: a,b,step,h,tol,romb,erro,xi
        real(kind=8), allocatable :: R1(:),R2(:)
        integer(kind=8) :: i,j, maxIter,n
        integer :: tipo
        
        step = b-a
        
        allocate(R1(1))
        allocate(R2(2))
        
        if (tipo .EQ. int(1)) then
        
            R1(1) = (step/2)*(f(a)*phi(a,h,xi,n)+f(b)*phi(a,h,xi,n))
            step = step/2
            i=2
            R2(1) = trapez(step,R1(1),a,b,i-1,xi,h,n,tipo)
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
                
                step = step/2
                R2(1) = trapez(step,R1(1),a,b,i-1,xi,h,n,tipo)
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
            
        else if (tipo .EQ. int(2)) then
            R1(1) = (step/2)*(funk(a)*((phid(a,h,xi,n))**2) &
            +funq(a)*((phi(a,h,xi,n))**2)+funk(b)*((phid(b,h,xi,n))**2) &
            +funq(b)*((phi(a,h,xi,n))**2))
            step = step/2
            i=2
            R2(1) = trapez(step,R1(1),a,b,i-1,xi,h,n,tipo)
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
                
                step = step/2
                R2(1) = trapez(step,R1(1),a,b,i-1,xi,h,n,tipo)
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
            
        else if (tipo .EQ. int(3)) then
            R1(1) = (step/2)*(funk(a)*phid(a,h,xi,n)*phid(a,h,xi+h,n) &
            +funq(a)*phi(a,h,xi,n)*phi(a,h,xi+h,n)+funk(b)*phid(b,h,xi,n)*phid(b,h,xi+h,n) &
            +funq(b)*phi(b,h,xi,n)*phi(b,h,xi+h,n))
            step = step/2
            i=2
            R2(1) = trapez(step,R1(1),a,b,i-1,xi,h,n,tipo)
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
                
                step = step/2
                R2(1) = trapez(step,R1(1),a,b,i-1,xi,h,n,tipo)
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
            
        else if (tipo .EQ. int(4)) then
            R1(1) = (step/2)*(funk(a)*phid(a,h,xi,n)*phid(a,h,xi+2*h,n) &
            +funq(a)*phi(a,h,xi,n)*phi(a,h,xi+2*h,n)+funk(b)*phid(b,h,xi,n)*phid(b,h,xi+2*h,n) &
            +funq(b)*phi(b,h,xi,n)*phi(b,h,xi+2*h,n))
            step = step/2
            i=2
            R2(1) = trapez(step,R1(1),a,b,i-1,xi,h,n,tipo)
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
                
                step = step/2
                R2(1) = trapez(step,R1(1),a,b,i-1,xi,h,n,tipo)
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
                
        else
            R1(1) = (step/2)*(funk(a)*phid(a,h,xi,n)*phid(a,h,xi+3*h,n) &
            +funq(a)*phi(a,h,xi,n)*phi(a,h,xi+3*h,n)+funk(b)*phid(b,h,xi,n)*phid(b,h,xi+3*h,n) &
            +funq(b)*phi(b,h,xi,n)*phi(b,h,xi+3*h,n))
            step = step/2
            i=2
            R2(1) = trapez(step,R1(1),a,b,i-1,xi,h,n,tipo)
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
                
                step = step/2
                R2(1) = trapez(step,R1(1),a,b,i-1,xi,h,n,tipo)
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
            
        end if
    
    end function
    
    !------------------------------------------------------------------------------------------
    !------------------------------------------ Trapezio ------------------------------------
    !------------------------------------------------------------------------------------------
    
     function trapez(step,th,a,b,i,xi,h,n,tipo)
        implicit none
        
        real(kind=8) :: a,b,step,trapez,aux,th,h,xi,aux2,aux3,um
        integer(kind=8) :: i,j,n,j1,j2
        integer :: tipo
        
        um = 1
        aux2 = 2
        aux3 = 3
        j1 = 1
        j2 = 2
        
        trapez = th/aux2
        
        if (tipo .EQ. int(1)) then
            do j=j1,j2**(i-j1)
                aux = a+(aux2*real(j)-um)*step
                trapez = trapez + step*f(aux)*phi(aux,h,xi,n)
            end do
        elseif (tipo .EQ. int(2)) then
            do j=j1,j2**(i-j1)
                aux = a+(aux2*real(j)-um)*step
                trapez = trapez + step*(funk(aux)*phid(aux,h,xi,n)*phid(aux,h,xi,n) &
                +funq(aux)*phi(aux,h,xi,n)*phi(aux,h,xi,n))
            end do
        elseif (tipo .EQ. int(3)) then
            do j=j1,j2**(i-j1)
                aux = a+(aux2*real(j)-um)*step
                trapez = trapez + step*(funk(aux)*phid(aux,h,xi,n)*phid(aux,h,xi+h,n) &
                +funq(aux)*phi(aux,h,xi,n)*phi(aux,h,xi+h,n))
            end do
        elseif (tipo .EQ. int(4)) then
            do j=j1,j2**(i-j1)
                aux = a+(aux2*real(j)-um)*step
                trapez = trapez + step*(funk(aux)*phid(aux,h,xi,n)*phid(aux,h,xi+aux2*h,n) &
                +funq(aux)*phi(aux,h,xi,n)*phi(aux,h,xi+aux2*h,n))
            end do
        else
            do j=j1,j2**(i-j1)
                aux = a+(aux2*real(j)-um)*step
                trapez = trapez + step*(funk(aux)*phid(aux,h,xi,n)*phid(aux,h,xi+aux3*h,n) &
                +funq(aux)*phi(aux,h,xi,n)*phi(aux,h,xi+aux3*h,n))
            end do
        end if
        
    end function
    
    !--------------------------------------------------------------------------------------------------
    !------------------------------------------ LU heptadiagonal ----------------------------------
    !---------------------------------------------------------------------------------------------
    
    subroutine lu7diag(a,b,c,d,alfa,beta,gama,theta,delta,omega,n)
        implicit none
        integer(kind=8), intent(in) :: n
        real(kind=8), dimension(n), intent(in) :: a
        real(kind=8), dimension(n-1), intent(in) :: b
        real(kind=8), dimension(n-2), intent(in) :: c
        real(kind=8), dimension(n-3), intent(in) :: d
        real(kind=8), dimension(n), intent(out) :: alfa
        real(kind=8), dimension(n-1), intent(out) :: beta
        real(kind=8), dimension(n-2), intent(out) :: gama
        real(kind=8), dimension(n-1), intent(out) :: theta
        real(kind=8), dimension(n-2), intent(out) :: omega
        real(kind=8), dimension(n-3), intent(out) :: delta
        
        integer(kind=8) :: i,j1,j2,j3,j4
        
        j1 = 1
        j2 = 2
        j3 = 3
        j4 = 4
        
        alfa(j1) = a(j1)
        beta(j1) = b(j1) ! i=1
        gama(j1) = c(j1)
        
        theta(j1) = b(j1)/alfa(j1)
        alfa(j2) = a(j2) - theta(j1)*beta(j1)  ! i=2
        beta(j2) = b(j2) - theta(j1)*gama(j1)
        gama(j2) = c(j2) - theta(j1)*d(j1)
        
        omega(j1) = c(j1)/alfa(j1)
        theta(j2) = (b(j2)-omega(j1)*beta(j1))/alfa(j2)
        alfa(j3) = a(j3) - omega(j1)*gama(j1) - theta(j2)*beta(j2)  ! i=3
        beta(j3) = b(j3) - omega(j1)*d(j1) - theta(j2)*gama(j2)
        gama(j3) = c(j3) - theta(j2)*d(j2)
        
        do i=j4,n-j2
            delta(i-j3) = d(i-j3)/alfa(i-j3)
            omega(i-j2) = (c(i-j2)-delta(i-j3)*beta(i-j3))/alfa(i-j2)   ! i=4,...,n-2
            theta(i-j1) = (b(i-j1)-delta(i-j3)*gama(i-j3)-omega(i-j2)*beta(i-j2))/alfa(i-j1)
            alfa(i) = a(i) - delta(i-j3)*d(i-j3) - omega(i-j2)*gama(i-j2) - theta(i-j1)*beta(i-j1)
            beta(i) = b(i) - omega(i-j2)*d(i-j2) - theta(i-j1)*gama(i-j1)
            gama(i) = c(i) - theta(i-j1)*d(i-j1)
        end do
        
        delta(n-j4) = d(n-j4)/alfa(n-j4)
        omega(n-j3) = (c(n-j3)-delta(n-j4)*beta(n-j4))/alfa(n-j3)
        theta(n-j2) = (b(n-j2)-delta(n-j4)*gama(n-j4)-omega(n-j3)*beta(n-j3))/alfa(n-j2)   ! i=n-1
        alfa(n-j1) = a(n-j1) - delta(n-j4)*d(n-j4) - omega(n-j3)*gama(n-j3) - theta(n-j2)*beta(n-j2)
        beta(n-j1) = b(n-j1) - omega(n-j3)*d(n-j3) - theta(n-j2)*gama(n-j2)
        
        delta(n-j3) = d(n-j3)/alfa(n-j3)
        omega(n-j2) = (c(n-j2)-delta(n-j3)*beta(n-j3))/alfa(n-j2)
        theta(n-j1) = (b(n-j1)-delta(n-j3)*gama(n-j3)-omega(n-j2)*beta(n-j2))/alfa(n-j1)   ! i=n
        alfa(n) = a(n) - delta(n-j3)*d(n-j3) - omega(n-j2)*gama(n-j2) - theta(n-j1)*beta(n-j1)
        
    end subroutine lu7diag
    
    !----------------------------------------------------------------------------------------------------------
    !------------------------------- Solução do Sistema usando LU hepta --------------------------------------
    !----------------------------------------------------------------------------------------------------------
    
    function solLu7diag(a,b,c,d,ind,n)
        implicit none
        
        real(kind=8), dimension(n) :: a,ind
        real(kind=8), dimension(n-1) :: b
        real(kind=8), dimension(n-2) :: c
        real(kind=8), dimension(n-3) :: d
        real(kind=8), dimension(n) :: alfa
        real(kind=8), dimension(n-1) :: beta
        real(kind=8), dimension(n-2) :: gama
        real(kind=8), dimension(n-1) :: theta
        real(kind=8), dimension(n-2) :: omega
        real(kind=8), dimension(n-3) :: delta
        integer(kind=8) :: i,n,j1,j2,j3,j4
        real(kind=8), dimension(n) :: solLu7diag
        real(kind=8), dimension(n) :: aux
        
        j1 = 1
        j2 = 2
        j3 = 3
        j4 = 4
        
        call lu7diag(a,b,c,d,alfa,beta,gama,theta,delta,omega,n)
    
        aux(j1) = ind(j1)/alfa(j1)
        aux(j2) = (ind(j2) - beta(j1)*aux(j1))/alfa(j2)
        aux(j3) = (ind(j3) - gama(j1)*aux(j1) - beta(j2)*aux(j2))/alfa(j3)
        
        do i=j4,n
            aux(i) = (ind(i) - d(i-j3)*aux(i-j3) - gama(i-j2)*aux(i-j2) - beta(i-j1)*aux(i-j1))/alfa(i)
        end do
        
        solLu7diag(n) = aux(n)
        solLu7diag(n-j1) = aux(n-j1) - theta(n-j1)*solLu7diag(n)
        solLu7diag(n-j2) = aux(n-j2) - theta(n-j2)*solLu7diag(n-j1) - omega(n-j2)*solLu7diag(n)
        
        do i=j3,n-j1
            solLu7diag(n-i) = aux(n-i) - theta(n-i)*solLu7diag(n-i+j1) - omega(n-i)*solLu7diag(n-i+j2) &
            - delta(n-i)*solLu7diag(n-i+j3)
        end do
        
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
        
    end function
    
    !------------------------------------------------------------------------------------------
    !------------------------------------------ Função min ------------------------------------
    !------------------------------------------------------------------------------------------
    
     function funmin(a,b)
        implicit none
        
        real(kind=8) :: funmin,a,b
        
        if (a .LE. b) then
            funmin = a
        else
            funmin = b
        endif
        
    end function
    
    !------------------------------------------------------------------------------------------------
    !----------------------------- Função B-Spline Cúbico --------------------------------------------
    !------------------------------------------------------------------------------------------------

     function Bcub(x)
         implicit none
         
         real(kind=8) :: x,Bcub,aux2,aux4,aux1,aux0
         
         aux4 = 4
         aux2 = 2
         aux1 = 1
         aux0 = 0
         
         if (x .LT. -aux2) then
            Bcub = aux0
         elseif (x .LT. -aux1) then
            Bcub = (aux1/aux4)*((aux2+x)**3)
         elseif (x .LT. aux0) then
            Bcub = (aux1/aux4)*((aux2+x)**3 - aux4*((aux1+x)**3))
         elseif (x .LT. aux1) then
            Bcub = (aux1/aux4)*((aux2-x)**3 - aux4*((aux1-x)**3))
         elseif (x .LT. aux2) then
            Bcub = (aux1/aux4)*((aux2-x)**3)
         else
            Bcub = aux0
         endif

     end function
     
    !------------------------------------------------------------------------------------------------
    !--------------------- Função B-Spline Cúbico ------------------------------------------------
    !------------------------------------------------------------------------------------------------

     function Bcubj(x,xi,h)
         implicit none
         
         real(kind=8) :: x,Bcubj,xi,h
         
         Bcubj = Bcub((x-xi)/h)

     end function
     
    !------------------------------------------------------------------------------------------------
    !--------------------- Função Spline Cúbico (Chapéu) --------------------------------------------
    !------------------------------------------------------------------------------------------------

     function phi(x,h,xi,n)
        implicit none
         
        real(kind=8) :: x,phi,h,xi,aux2,aux4
        integer(kind=8) :: n
         
        aux4 = 4
        aux2 = 2
         
        if (int(xi/h) .EQ. int(0)) then
            phi = Bcubj(x,xi,h)-aux4*Bcubj(x,-h,h)
        elseif (int(xi/h) .EQ. int(1)) then
            phi = Bcubj(x,xi,h)-Bcubj(x,xi-aux2*h,h)
        elseif (int(xi/h) .GE. int(2) .AND. int(xi/h) .LE. n-1) then
            phi = Bcubj(x,xi,h)
        elseif (int(xi/h) .EQ. n) then
            phi = Bcubj(x,xi,h)-Bcubj(x,xi+aux2*h,h)
        else
            phi = Bcubj(x,xi,h)-aux4*Bcubj(x,xi+h,h)
        endif
        
     end function

    !------------------------------------------------------------------------------------------------
    !--------------------- Função B-Spline Cúbico Derivada ------------------------------------------------
    !------------------------------------------------------------------------------------------------

     function Bcubd(x)
         implicit none
         
         real(kind=8) :: x,Bcubd, aux0, aux1, aux2, aux3, aux4
         
         aux0 = 0
         aux1 = 1
         aux2 = 2
         aux3 = 3
         aux4 = 4
         
         if (x .LT. -aux2) then
            Bcubd = aux0
         elseif (x .LT. -aux1) then
            Bcubd = (aux3/aux4)*((aux2+x)**2)
         elseif (x .LT. aux0) then
            Bcubd = (aux3/aux4)*((aux2+x)**2 - aux4*((aux1+x)**2))
         elseif (x .LT. aux1) then
            Bcubd = (aux3/aux4)*(-((aux2-x)**2) + aux4*((aux1-x)**2))
         elseif (x .LT. aux2) then
            Bcubd = -(aux3/aux4)*((aux2-x)**2)
         else
            Bcubd = aux0
         endif

     end function
     
    !------------------------------------------------------------------------------------------------
    !--------------------- Função B-Spline Cúbico Derivada ------------------------------------------------
    !------------------------------------------------------------------------------------------------

     function Bcubjd(x,xi,h)
         implicit none
         
         real(kind=8) :: x,Bcubjd,xi,h,um
         
         um = 1
         
         Bcubjd = (um/h)*Bcubd((x-xi)/h)

     end function
     
    !------------------------------------------------------------------------------------------------
    !--------------------- Função Spline Cúbico Derivada ------------------------------------------------
    !------------------------------------------------------------------------------------------------

     function phid(x,h,xi,n)
         implicit none
         
         real(kind=8) :: x,phid,h,xi,aux2,aux4
         integer(kind=8) :: j,n,j1
         
         j1 = 1
         aux2 = 2
         aux4 = 4
         
         if (int(xi/h) .EQ. int(0)) then
            phid = Bcubjd(x,xi,h)-aux4*Bcubjd(x,-h,h)
         elseif (int(xi/h) .EQ. int(1)) then
            phid = Bcubjd(x,xi,h)-Bcubjd(x,xi-aux2*h,h)
        elseif (int(xi/h) .GE. int(2) .AND. int(xi/h) .LE. n-j1) then
            phid = Bcubjd(x,xi,h)
        elseif (int(xi/h) .EQ. n) then
            phid = Bcubjd(x,xi,h)-Bcubjd(x,xi+aux2*h,h)
        else
            phid = Bcubjd(x,xi,h)-aux4*Bcubjd(x,xi+h,h)
        endif
        
     end function
     
    !------------------------------------------------------------------------------------------------
    !--------------------- Função q(x) Test ------------------------------------------------
    !------------------------------------------------------------------------------------------------

     function funq(x)
         implicit none
         
         real(kind=8) :: x,funq
         
         funq = 0
        
     end function
     
     !------------------------------------------------------------------------------------------------
     !--------------------- Função k(x) Test ------------------------------------------------
     !------------------------------------------------------------------------------------------------

     function funk(x)
         implicit none
         
         real(kind=8) :: x,funk
         
         funk = 1
        
     end function
     
    !------------------------------------------------------------------------------------------------
    !--------------------- Função f(x) Test ------------------------------------------------
    !------------------------------------------------------------------------------------------------

     function f(x)
         implicit none
         
         real(kind=8) :: x,f,aux1,aux2,aux12
         
         aux1 = 1
         aux2 = 2
         aux12 = 12
         
         f = aux12*x*(aux1-x)-aux2
        
     end function
    
end program
