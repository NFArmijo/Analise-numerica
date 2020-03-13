
module funcao_test_raiz

    contains
    function f(x,a)
        implicit none
        real(kind=8) :: f,x,a
        f = x**2-a
    end function
    
    function fder(x)
        implicit none
        real(kind=8) :: fder,x
        fder = 2*x
    end function

end module
