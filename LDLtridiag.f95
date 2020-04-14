program LDLtridiag

    implicit none
    
    real(kind=8), allocatable :: a(:),b(:),c(:)
    real(kind=8), allocatable :: solSist(:)
    integer :: n
    
    n = 4
    allocate(a(n))
    allocate(c(n))
    allocate(b(n-1))
    allocate(solSist(n))
    
    a =  (/10,20,30,40/)
    b =  (/1,2,2/)
    c = (/1,1,1,1/)
    
    solSist = thomas3diag(a,b,c,n)
    
    print*, 'The solution of the system Ax=b is: ', solSist

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
        
        integer :: i,j
        
        d(1) = a(1)
        l(1) = b(1)/d(1)
        
        do i=2,n-1
            d(i) = a(i)-b(i-1)*l(i-1)
            l(i) = b(i)/d(i)
        end do
        
        d(n) = a(n)-b(n-1)*l(n-1)
        
        print*, 'L : ', l
        print*, 'D : ', d
        
    end subroutine ldl3diag
    
    !----------------------------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------------------
    
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
        
    end function
    
    !----------------------------------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------------------------------
    
end program
