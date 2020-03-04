module funcao_test_poly

    contains
    function f(x,coef,deg)
        implicit none
        integer :: deg,i
        real(kind=8), dimension (deg+1) :: coef
        real(kind=8) :: f,x,aux
        
        f = coef(1)
        
        aux = 1.0
        do i =2,deg+1
            aux = aux*x
            f = f + aux*coef(i)
        end do
        
    end function

end module
