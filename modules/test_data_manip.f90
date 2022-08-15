program test_data_manip

    use iso_fortran_env, only: real32, real64
    use mod_data_manip
    use mod_linalg
    use mod_io
    implicit none

    !test high_to_low_diag and low_to_high_diag
    
    real(real64), dimension(20,20) :: A, A_high_to_low, A_low_to_high
    integer :: i

    A = 0
    do i = 1, size(A,1)
        A(i,i) = 73*1e-14
    end do
    A(1,1) = -59
    A(2,2) = 393
    A(3,3) = 2
    A(4,4) = -1

    A_high_to_low = high_to_low_diag(A)
    A_low_to_high = low_to_high_diag(A)
    call write_to_file(A,1)
    call write_to_file(A_high_to_low, 2)
    call write_to_file(A_low_to_high, 3)

end program test_data_manip