program test

    use iso_fortran_env, only: real32, real64
    use mod_data_manip
	use mod_derived_types
	use mod_grid_generation
	use mod_interpolation
    use mod_io
	use mod_linalg
	use mod_parallel
	use mod_runge_kutta
	use mod_finite_volume
    implicit none
    
    ! We attempt to solve q_t + Aq_x = 0, A = diag(1,-1)
    ! with (q_i)_0 = sin(2pi*x) and periodic boundary,
    ! using a simple finite volume method.
    
	integer, parameter :: n = 500
	integer :: i, j, k, m
	real, dimension(0:n) :: x, t
	real, dimension(:), allocatable :: vol_avg_q_0
	real, dimension(:,:), allocatable :: vol_avg_q, solution1
	real, dimension(2, 2) :: A
	real, parameter :: pi = 3.14159265358979323846
	real :: h
	h = 1./(n)
	m = size(A,1)
	allocate(vol_avg_q_0(m*n), solution1(n, n))
	do i = 0, n 
		x(i) = i*h
		t(i) = i*h
	end do
	A = 0
	A(1,1) = 1
	A(1,2) = 3
	A(2,2) = -1
	do i = 1, n
			vol_avg_q_0(i) = (1./(2*pi*h)) * (cos(2*pi*x(i-1)) - cos(2*pi*x(i)))
	end do
	do j = 2, m
		vol_avg_q_0(((j-1)*n)+1:j*n) = vol_avg_q_0(1:n)
	end do
	
	vol_avg_q = fully_discrete_second_order_LW_scheme_1dsys(h, A, 0., 1., vol_avg_q_0)

	print*, 'Maxval of Q is', maxval(vol_avg_q)

	call write_to_file(vol_avg_q(:, 1:n), 1)

	solution1(1,:) = vol_avg_q(1,1:n)

	do i = 2, n
		do j = 1, n
			solution1(i,j) = sin(2*pi*x(j) - t(i)) + 3*sin(2*pi*x(j) + t(i))
		end do
	end do

	call write_to_file(solution1(:,:), 2)

	!Solution matches vol_avg_q!

end program test
