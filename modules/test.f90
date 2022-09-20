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
    
    ! We attempt to solve q_t + aq_x = 0,
    ! with u_0 = sin(2pi*x) and periodic boundary,
    ! using a simple finite volume method.
    
	integer, parameter :: n = 500
	integer :: fileunit, i, j, k
	real, dimension(0:n) :: x, time
	real, dimension(1:n) :: vol_avg_q_0
	real, dimension(:,:), allocatable :: vol_avg_q
	real, parameter :: pi = 3.14159265358979323846
	real :: h = 1./(n)
	do i = 0, n 
		x(i) = i*h
	end do
	
	do i = 1, n
		vol_avg_q_0(i) = (1./(2*pi*h)) * (cos(2*pi*x(i-1)) - cos(2*pi*x(i)))
	end do
	
	vol_avg_q = fully_discrete_Richtmyer_LW_scheme(h=h, speed=1., tmin=0., tmax=1., q_0=vol_avg_q_0)

	print*, 'Maxval of Q is', maxval(vol_avg_q)

	call write_to_file(vol_avg_q, 1)

end program test
