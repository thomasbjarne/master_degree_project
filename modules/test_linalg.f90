program test_linalg

    use iso_fortran_env, only: real32, real64
    use mod_linalg
    use mod_data_manip
    use mod_io
    implicit none
    
! test qr_decomposition

    real(real64), dimension(:,:), allocatable :: A, Q, R
    integer :: i, n = 3

    allocate( A(n+1,n), Q(n+1,n), R(n,n))

    A = 0
    do i = 1, n
        A(i,i) = 4
    end do

    call qr_decomposition(A, Q, R)

    call write_to_file(Q, 1)
    call write_to_file(R, 2)
    call write_to_file(matmul(Q, R), 3)
    call write_to_file(A, 4)


end program test_linalg