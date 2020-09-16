program EP_Final_lin

    implicit none
    
    real(kind=8), allocatable :: a(:),b(:),d(:) !Matriz e vetor independente
    real(kind=8), allocatable :: c(:) !Solução do sistema
    real(kind=8), allocatable :: u(:) !Solução aproximada
    real(kind=8), allocatable :: ureal(:) !Solução real da EDO
    real(kind=8), allocatable :: udif(:) !Diferença das soluções
    real(kind=8) :: erro, x0, xn,tol, h
    real(kind=8) :: aux,aux0,aux1,aux2
    integer(kind=8) :: n,i,j,k,maxIter
    
    !------------------------------------------------------------------------------------
    !----------------------------------------- Exercício 1 -----------------------------------
    !---------------------------------------------------------------------------------------
        
    open(1, file = 'solEpfinal16exp.txt', status = 'new')
    
    tol = 0.000001d0
    maxIter = 1000
    
    do j=4,4

        n = 2**j-1
        
        allocate(a(n))      !Diagonal da matriz
        allocate(b(n-1))    !Subdiagonal da matriz
        allocate(d(n))      !Vetor independente
        allocate(c(n))      !Solução do sistema
        allocate(u(10*n+1))      !Aproximação da função solução da EDO
        allocate(ureal(10*n+1))  !Função solução da EDO
        allocate(udif(10*n+1))   !Diferença das soluções
    
        aux = real(n)+1
        h = 1/aux
        
        do i=1,n-1
            x0 = (real(i)-1)*h
            xn = (real(i)+1)*h
            d(i) = romb(x0,xn,tol,h,real(i)*h,maxIter,int(1))
            a(i) = romb(x0,xn,tol,h,real(i)*h,maxIter,int(2))
            x0 = real(i)*h
            b(i) = romb(x0,xn,tol,h,real(i)*h,maxIter,int(3))
        end do
        
        x0 = (real(n)-1)*h
        xn = (real(n)+1)*h
        d(n) = romb(x0,xn,tol,h,real(n)*h,maxIter,int(1))
        a(n) = romb(x0,xn,tol,h,real(n)*h,maxIter,int(2))
    
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
            ureal(i) = (((real(i)-1)/(10*real(n)))**2)*((1-((real(i)-1)/(10*real(n))))**2)
            udif(i) = u(i)-ureal(i)
        end do
        
        u(10*n+1) = 0.0d0
        ureal(10*n+1) = 0.0d0
        udif(10*n+1) = u(10*n+1)-ureal(10*n+1)
    
        erro = norm_inf(udif,10*n+1)
        
        print*, 'Error of the approximation of the sol of the EDO: ', erro
        print*, ' '
        !print*, 'Approximated solution of the EDO: ', u
        !print*, ' '
        
        !write(1,*) erro, (aux**2)*erro
        write(1,*) u
        
        deallocate(a)
        deallocate(b)
        deallocate(d)
        deallocate(c)
        deallocate(u)
        deallocate(ureal)
        deallocate(udif)
        
    end do
    
    close(1)
    
    contains
    
    !------------------------------------------------------------------------------------
    !----------------------------------------- Método de Romberg ---------------------------
    !---------------------------------------------------------------------------------------
    
    
    function romb(a,b,tol,h,xi,maxIter,tipo)
        implicit none
        
        real(kind=8) :: a,b,step,h,tol,romb,erro,xi
        real(kind=8), allocatable :: R1(:),R2(:)
        integer(kind=8) :: i,j, maxIter
        integer :: tipo 
        
        step = b-a
        
        allocate(R1(1))
        allocate(R2(2))
        
        if (tipo .EQ. int(1)) then
        
            R1(1) = (step/2)*(f(a)*phi(a,h,xi)+f(b)*phi(a,h,xi))
            
            step = step/2
            
            i=2
            
            R2(1) = trapez(step,R1(1),a,b,i-1,xi,h,tipo)
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
                R2(1) = trapez(step,R1(1),a,b,i-1,xi,h,tipo)
                
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
                
            R1(1) = (step/2)*(funk(a)*phid(a,h,xi)*phid(a,h,xi+h) &
            +funq(a)*phi(a,h,xi)*phi(a,h,xi+h)+funk(b)*phid(b,h,xi)*phid(b,h,xi+h) &
            +funq(b)*phi(b,h,xi)*phi(b,h,xi+h))
            
            step = step/2
            
            i=2
            
            R2(1) = trapez(step,R1(1),a,b,i-1,xi,h,tipo)
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
                R2(1) = trapez(step,R1(1),a,b,i-1,xi,h,tipo)
                
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
                
            R1(1) = (step/2)*(funk(a)*phid(a,h,xi)*phid(a,h,xi) &
            +funq(a)*phi(a,h,xi)*phi(a,h,xi)+funk(b)*phid(b,h,xi)*phid(b,h,xi) &
            +funq(b)*phi(b,h,xi)*phi(b,h,xi))
            
            step = step/2
            
            i=2
            
            R2(1) = trapez(step,R1(1),a,b,i-1,xi,h,tipo)
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
                R2(1) = trapez(step,R1(1),a,b,i-1,xi,h,tipo)
                
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
    
     function trapez(step,th,a,b,i,xi,h,tipo)
        implicit none
        
        real(kind=8) :: a,b,step,trapez,aux,th,h,xi
        integer(kind=8) :: i,j
        integer :: tipo
        
        trapez = th/2
        
        if (tipo .EQ. int(1)) then
            do j=1,2**(i-1)
                aux = a+(2*real(j)-1)*step
                trapez = trapez + step*f(aux)*phi(aux,h,xi)
            end do
        else if (tipo .EQ. int(3)) then
            do j=1,2**(i-1)
                aux = a+(2*real(j)-1)*step
                trapez = trapez + step*(funk(aux)*phid(aux,h,xi)*phid(aux,h,xi+h) &
                +funq(aux)*phi(aux,h,xi)*phi(aux,h,xi+h))
            end do
        else
            do j=1,2**(i-1)
                aux = a+(2*real(j)-1)*step
                trapez = trapez + step*(funk(aux)*phid(aux,h,xi)*phid(aux,h,xi) &
                +funq(aux)*phi(aux,h,xi)*phi(aux,h,xi))
            end do
        end if
        
    end function
    
    !--------------------------------------------------------------------------------------------------
    !------------------------------------------ LDL^T tridiagonal ----------------------------------
    !---------------------------------------------------------------------------------------------
    
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
    
    !--------------------- Função Spline Linear (Chapéu) ------------------------------------------------

     function phi(x,h,xi)
         implicit none
         
         real(kind=8) :: x,phi,h,xi
         
         if (x .LE. xi .AND. x .GE. xi-h) then
            phi = (x-(xi-h))/h
         elseif (x .LE. xi+h .AND. x .GE. xi) then
            phi = ((xi+h)-x)/h
        endif
        
     end function
     
    !--------------------- Função Spline Linear Derivada (Chapéu) ------------------------------------------------

     function phid(x,h,xi)
         implicit none
         
         real(kind=8) :: x,phid,h,xi
         
         if (x .LE. xi .AND. x .GE. xi-h) then
            phid = 1/h
         elseif (x .LE. xi+h .AND. x .GE. xi) then
            phid = -1/h
         endif
        
     end function
    
        !--------------------- Função q(x) Test ------------------------------------------------

     function funq(x)
         implicit none
         
         real(kind=8) :: x,funq
         
         funq = exp(x)
        
     end function
     
         !--------------------- Função k(x) Test ------------------------------------------------

     function funk(x)
         implicit none
         
         real(kind=8) :: x,funk
         
         funk = x**3
        
     end function
     
              !--------------------- Função f(x) Test ------------------------------------------------

     function f(x)
         implicit none
         
         real(kind=8) :: x,f
         
         f = 12*x*(1-x)-2
        
     end function
    
end program
