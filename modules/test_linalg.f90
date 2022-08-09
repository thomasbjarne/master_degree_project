program test_linalg

    use mod_linalg
    use mod_io
    implicit none
    
    real, dimension(:,:), allocatable :: A, A_schur, lambda, R, test
    integer :: i, n = 10

    allocate( A(n,n), A_schur(n,n), lambda(n,n), R(n,n), test(n,n))

    A = 9
    do i = 1, n
        A(i,i) = 4
    end do

    !A(2, 4) = 2
    !A(1, 2) = -5
    !A(3, 1) = 29

    A_schur = qr_algorithm(A)

    call write_to_file(A_schur, 1)

    lambda = 0
    do concurrent (i = 1:n)
        lambda(i,i) = A_schur(i,i)
    end do

    R = inverse_iteration(A,lambda)

    call write_to_file(R,2)

    !Test R(:,i) * lambda(i,i) = A * R(:,i) (verify eigenvec and eigenval)
    do concurrent (i=1:n)
        test(:,i) = abs( R(:,i) * lambda(i,i) - matmul(A, R(:,i) ) )
    end do

    call write_to_file(test,3)

end program test_linalg