
module funcao_test1

    contains
    function f(x)
        implicit none
        real*8 :: f,x
        f = x**2-1
    end function

end module
