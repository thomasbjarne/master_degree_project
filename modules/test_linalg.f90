program test_linalg

    use iso_fortran_env, only: real32, real64
    use mod_linalg
    use mod_data_manip
    use mod_io
    implicit none
    
! test inverse_iteration

    real(real64), dimension(:,:), allocatable :: A, A_schur, lambda, R, AR, lambdaR
    integer :: i, n = 3

    allocate( A(n,n), A_schur(n,n), lambda(n,n), R(n,n), AR(n,n), lambdaR(n,n))

    !EIGENVALUES SHOULD BE 4, 2, 1 -- We find the correct eigenvalues
    A(1,:) = [0, 11, -5]
    A(2,:) = [-2, 17, -7]
    A(3,:) = [-4, 26, -10]
 
    A_schur = qr_algorithm(A)

    call write_to_file(A_schur, 1)

    lambda = 0
    do i = 1, n
        lambda(i,i) = A_schur(i,i)
    end do

    lambda = high_to_low_diag(lambda)
    call write_to_file(lambda, 2)

    !something goes wrong after this point

    R = inverse_iteration(A,lambda)

    call write_to_file(R,3)

    !Test R(:,i) * lambda(i,i) = A * R(:,i) (verify eigenvec and eigenval)
    do i = 1, n
        AR = matmul(A, R)
        lambdaR(:,i) = lambda(i,i) * R(:,i)
    end do

    call write_to_file(AR, 4)
    call write_to_file(lambdaR, 5)


end program test_linalg