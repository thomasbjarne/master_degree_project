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
    implicit none
    
    ! We attempt to solve q_t + aq_x = 0,
    ! with u_0 = sin(2pi*x) and periodic boundary,
    ! using a simple finite volume method.
    
	integer, parameter :: n = 500
	integer :: fileunit, i, j, k
	real, dimension(0:n) :: x, t
	real, dimension(1:n) :: vol_avg_q_0
	real, dimension(:,:), allocatable :: vol_avg_q
	real, parameter :: pi = 3.14159265358979323846
	real :: h, a

	a = 1.
	h = 1./(n)
	do i = 0, n 
		x(i) = i*h
	end do
	
	do i = 1, n
		vol_avg_q_0(i) = (1./(2*pi*h)) * (cos(2*pi*x(i-1)) - cos(2*pi*x(i)))
	end do
	
	call rk4_1D(Lax_Friedrichs_scheme, vol_avg_q_0, h, 0., 1., vol_avg_q)

	print*, 'Maxval of Q is', maxval(vol_avg_q)

	call write_to_file(vol_avg_q, 1)

contains

	pure function upwind_scheme(t, y) result(res)
		
		real, intent(in), dimension(:) :: y
		real, intent(in), optional :: t
		real, dimension(size(y)) :: res
		integer :: i

		!If a > 0:
		res(1) = -(a/h)*(y(1) - y(size(y)))
		do i = 2, size(y)
			res(i) = -(a/h)*(y(i)-y(i-1))
		end do
		

		!If a < 0:
		!do i = 1, size(y)-1
		!	res(i) = -(a/h)*(y(i+1) - y(i))
		!end do
		!res(size(y)) = -(a/h)*(y(1) - y(size(y)))

	end function upwind_scheme

	pure function naive_central_scheme(t,y) result(res)

		real, intent(in), dimension(:) :: y
		real, intent(in), optional :: t
		real, dimension(size(y)) :: res
		integer :: i

		res(1) = -(a/(2*h))*(y(2)-y(size(y)))
		do i = 2, size(y)-1
			res(i) = -(a/(2*h))*(y(i+1)-y(i-1))
		end do
		res(size(y)) = -(a/(2*h))*(y(1)-y(size(y)-1))

	end function naive_central_scheme

	pure function Lax_Friedrichs_scheme(t,y) result(res)

		real, intent(in), dimension(:) :: y
		real, intent(in), optional :: t
		real, dimension(size(y)) :: res
		integer :: i

		res(1) = -(a/(2*h))*(y(2)-y(size(y))) + (0.5)*(y(2) + y(size(y)))
		do i = 2, size(y)-1
			res(i) = -(a/(2*h))*(y(i+1)-y(i-1)) + (0.5)*(y(i+1)+y(i-1))
		end do
		res(size(y)) = -(a/(2*h))*(y(1)-y(size(y)-1)) + (0.5)*(y(1) + y(size(y)-1))

	end function Lax_Friedrichs_scheme
	

end program test
