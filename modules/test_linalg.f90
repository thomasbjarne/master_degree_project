program test_linalg

    use iso_fortran_env, only: real32, real64
    use mod_linalg
    use mod_io
    implicit none
    
    ! test forwards_subst

    real(real64), dimension(:,:), allocatable :: L
    real(real64), dimension(:), allocatable :: b, x
    integer :: i, n

    n = 5
    allocate( L(n,n), b(n), x(n))

    L = 0
    do i = 1, n
        L(i,i) = 2
    end do
    L(2,1) = -53
    L(3,1) = 31
    L(3,2) = -885
    L(4,1) = 8
    L(4,2) = 79
    L(4,3) = -44

    b = 4
    b(1) = -15
    b(2) = 42
    x = forwards_subst(L, b)

    call write_to_file(L, 1)
    call write_to_file(matmul(L, x), 2)
    call write_to_file(b, 3)


end program test_linalg