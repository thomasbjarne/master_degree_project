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
    
	!We test our bowyer_watson implementation

	real, dimension(5, 2) :: pointset
	type(triangle), dimension(:), allocatable :: triangulation
	integer :: i, j, fileunit

	pointset(1,1) = 0
	pointset(1,2) = 0
	pointset(2,1) = 1
	pointset(2,2) = 0
	pointset(3,1) = 1
	pointset(3,2) = 1
	pointset(4,1) = 0
	pointset(4,2) = 1
	pointset(5,1) = 0.5
	pointset(5,2) = 0.5
	triangulation = bowyer_watson(pointset)

	open(newunit=fileunit, file='data/datafile1.txt')

	do i = 1, size(triangulation)
		do j = 1, 3
			write(unit=fileunit, fmt=*) triangulation(i) % vertices(j,:)
		end do
	end do

	close(unit=fileunit)


end program test
