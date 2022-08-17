program test_linalg

    use iso_fortran_env, only: real32, real64
    use mod_linalg
    use mod_data_manip
    use mod_io
    implicit none
    
! test qr_algorithm

    real(real64), dimension(:,:), allocatable :: A, A_schur
    integer :: i, n = 3

    allocate( A(n,n), A_schur(n,n))

    A(1,:) = [0, 11, -5]
    A(2,:) = [-2, 17, -7]
    A(3,:) = [-4, 26, -10]

    A_schur = qr_algorithm(A, shift='w')

    call write_to_file(A_schur, 1)


end program test_linalg