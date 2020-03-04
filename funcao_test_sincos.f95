module funcao_test_sincos

    contains
    function f(x,n1,n2,p1,p2,c1,c2)
        implicit none
        
        real(kind=8) :: x,n1,n2,p1,p2,c1,c2,f
        
        f = c1*sin(n1*x+p1) + c2*cos(n2*x+p2)
        
    end function

end module
