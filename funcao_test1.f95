
module funcao_test1

    contains
    function f(x)
        implicit none
        real*8 :: f,x
        f = x**2-1
    end function
    
    function fder(x)
        implicit none
        real*8 :: fder,x
        fder = 2*x
    end function

end module
