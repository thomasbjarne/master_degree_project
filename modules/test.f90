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
    
	real, dimension(3,2) :: points
	real :: result

	points(1,:) = [0,1]
	points(2,:) = [2,2]
	points(3,:) = [3,4]
!	points(4,:) = [3,-1]
	result = newton_div_diff_interpolate(points, 3.)
	print*, result

end program test