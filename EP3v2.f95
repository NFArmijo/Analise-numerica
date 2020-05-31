program EP3

    implicit none
    
    real(kind=8), allocatable :: solSist(:,:),X0(:,:)
    integer :: n,max_iter !tamanho do sistema e numero maximo de iteracoes
    real(kind=8) :: w,tol !parámetro do método (1,2) e tolerancia para aproximacao
    integer(kind=8), dimension(101) :: iter_vec
    integer(kind=8) :: min_iter !indice tal que numero de iteraciones é menor
    integer(kind=8) :: num_iter
    integer :: i
    
    !---------- EXEMPLO ---------
    
    tol = 0.000001d0 !10^-6
    max_iter = 10000
    n = 256 !completar
    
    allocate(solSist(n-1,n-1),X0(n-1,n-1))
    
    ! PONTO INICIAL 
    
    X0 = 0.0d0
    !CALL RANDOM_NUMBER(X0)
    
    ! CICLO PARA ENCONTRAR W
    
    do i=1,101
        w = 0.99d0+(real(i)/100)
        call sor_method(X0,n,w,max_iter,tol,0,solSist,num_iter)
        iter_vec(i) = num_iter
        print*, 100*w-100, '%'
    end do
    
    call elem_min(iter_vec,101,min_iter)
    print*, 'Iteracao final....'
    call sor_method(X0,n,0.99d0+(real(min_iter)/100),max_iter,tol,1,solSist,num_iter)
    print*, 'Finalizado'
    
    contains
    
    !------------------------------------------------------------------------------------------
    !------------------------------------------Norma infinito ------------------------------------
    !------------------------------------------------------------------------------------------
    
     function norm_inf(M,d)
        implicit none
        
        integer :: d,i,j
        real(kind=8), dimension(d,d) :: M
        real(kind=8) :: norm_inf,aux
        
        norm_inf = abs(M(1,1))
        
        do i=1,d
            do j=1,d
                aux = abs(M(i,j))
                if (aux > norm_inf) then
                    norm_inf = aux
                end if
            end do
        end do
    end
    
    !------------------------------------------------------------------------------------------
    !------------------------------------------ Minimo do vetor ------------------------------------
    !------------------------------------------------------------------------------------------
    
    subroutine elem_min(v,d,m)
        implicit none
        
        integer(kind=8), dimension(d), intent(in) :: v
        integer, intent(in) :: d
        integer(kind=8), intent(out) :: m
        real(kind=8) :: aux
        integer :: i
        
        aux = v(1)
        
        do i=1,d
            if (v(i) < aux) then
                aux = v(i)
                m = i
            end if
        end do
        
    end subroutine
    
    !----------------------------------------------------------------------------------------------------------
    !------------------------------------- Método SOR(X0,n,w,max_iter,tol) ---------------------------------------------
    !----------------------------------------------------------------------------------------------------------
    
    subroutine sor_method(X0,n,w,max_iter,tol,wrt,sor,num_iter)
        implicit none

        real(kind=8), dimension(n-1,n-1), intent(in) :: X0
        integer, intent(in) :: n,max_iter,wrt
        real(kind=8), intent(in) :: tol,w
        
        real(kind=8), dimension(n-1,n-1), intent(out) :: sor
        integer(kind=8), intent(out) :: num_iter
        
        integer :: i,j
        real(kind=8), dimension(n-1,n-1) :: aux
        real(kind=8) :: erro
        
        sor = X0
        erro = 10*tol
        
        num_iter = 0
        
        do while (num_iter .LE. max_iter .AND. erro .GE. tol)
            
            aux = sor
            
            ! ------------------- iteração para i = 1 -----------------------
            
            sor(1,1) = (1-w)*sor(1,1) + (w/4)*(-3+sor(1,2)-3+sor(2,1))
        
            do j=2,n-2
                sor(1,j) = (1-w)*sor(1,j) + (w/4)*(sor(1,j-1)+sor(1,j+1)-3+sor(2,j))
            end do
            
            sor(1,n-1) = (1-w)*sor(1,n-1) + (w/4)*(sor(1,n-2)+((6/n)-3)-3+sor(2,n-1))
            
            ! ------------------- iteração para i = 2,..,n-2 -----------------------
            
            do i=2,n-2
            
                sor(i,1) = (1-w)*sor(i,1) + (w/4)*(-3+sor(i,2)+sor(i-1,1)+sor(i+1,1))
            
                do j=2,n-2
                    sor(i,j) = (1-w)*sor(i,j) + (w/4)*(sor(i,j-1)+sor(i,j+1)+sor(i-1,j)+sor(i+1,j))
                end do
                
                sor(i,n-1) = (1-w)*sor(i,n-1) + (w/4)*(sor(i,n-2)+(6*(i/n)-3)+sor(i-1,n-1)+sor(i+1,n-1))
                
            end do
            
            ! --------------------- iteração para i = n-1 -----------------------
            
            sor(n-1,1) = (1-w)*sor(n-1,1) + (w/4)*( -3+sor(n-1,2)+sor(n-2,1)+((6/n)-3) )
        
            do j=2,n-2
                sor(n-1,j) = (1-w)*sor(n-1,j) + (w/4)*( sor(n-1,j-1)+sor(n-1,j+1)+sor(n-2,j)+(6*(j/n)-3) )
            end do
            
            sor(n-1,n-1) = (1-w)*sor(n-1,n-1) + (w/4)*( sor(n-1,n-2)+(3-(6/n))+sor(n-2,n-1)+(3-(6/n)) )
            
            !--------------------------------------------------------------------------------
            
            erro = norm_inf(sor-aux,n-1) !erro aproximado
            
            num_iter=num_iter+1 ! número de iterações cresce em 1
            
        end do
        
        
        if (wrt .EQ. 1) then
            open(1, file = 'resultadosEP3.txt', status = 'new')
            write(1,*) 'Parametro w otimo: ', w
            write(1,*) 'Numero de iteracoes w otimo: ', num_iter
            write(1,*) 'Erro da apriximacao: ', erro
            write(1,*) 'Solucao aproximada: ', sor
            close(1)
        end if
        
    end subroutine
    
    !----------------------------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------------------
    
end program
