module funcao_test_poly

    contains
    function f(x,coef,deg)
        implicit none
        integer :: deg,i
        real(kind=8), dimension (deg+1) :: coef
        real(kind=8) :: f,x,aux
        
        f = coef(1)
        
        aux = 1.0d0
        do i =2,deg+1
            aux = aux*x
            f = f + aux*coef(i)
        end do
        
    end function
    
    function fder(x,coef,deg)
        implicit none
        integer :: deg,i
        real(kind=8), dimension (deg+1) :: coef
        real(kind=8) :: fder,x,aux
        
        fder = coef(2)
        
        aux = 1.0d0
        do i =2,deg
            aux = aux*x
            fder = fder + i*aux*coef(i+1)
        end do
        
    end function

end module
