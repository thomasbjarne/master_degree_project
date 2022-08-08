program test_linalg

    use mod_linalg
    use mod_io
    implicit none

    ! test upper_hessenberg

    real, dimension(10,10) :: A, H
    integer :: i

    A = -192
    do i = 1, size(A,2)
        A(i,i) = 4
    end do
    A(1,2) = 2
    A(9,2) = 5

    H = upper_hessenberg(A)

    call write_to_file(H,1)




end program test_linalg