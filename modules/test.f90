program test

    use iso_fortran_env, only: real32, real64
    use mod_data_manip
	use mod_derived_types
	use mod_finite_diff
	use mod_grid_generation
	use mod_interpolation
    use mod_io
	use mod_linalg
	use mod_parallel
	use mod_runge_kutta
	use mod_finite_volume
    implicit none
    
    ! We attempt to solve q_t + Aq_x = 0,
    ! with (q_i)_0 = sin((2i-1)pi*x) and periodic boundary,
    ! using a simple finite volume method.
    
	integer, parameter :: n = 500
	integer :: i, j, k
	real, dimension(0:n) :: x
	real, dimension(1:n, 2) :: vol_avg_q_0
	real, dimension(:,:, :), allocatable :: vol_avg_q
	real, dimension(2, 2) :: A
	real, parameter :: pi = 3.14159265358979323846
	real :: h = 1./(n)
	do i = 0, n 
		x(i) = i*h
	end do
	A = 0
	A(1,1) = 1
	A(2,2) = -1
	do i = 1, n
		do j = 1, 2
			vol_avg_q_0(i, j) = (1./(2*pi*h)) * (cos((2*j -1)*pi*x(i-1)) - cos((2*j-1)*pi*x(i)))
		end do
	end do
	
	vol_avg_q = fully_discrete_second_order_LW_scheme_1dsys(h=h, a=A, tmin=0., tmax=1., q_0=vol_avg_q_0)

	print*, 'Maxval of Q is', maxval(vol_avg_q)

	call write_to_file(vol_avg_q(1,:,:), 1)

end program test
