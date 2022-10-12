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

	real, dimension(20, 2) :: pointset
	type(triangle), dimension(size(pointset,1)*10) :: triangulation
	integer :: i, j, fileunit

	do i = 1, 10
		pointset(i,1) = i
		pointset(i+10, 1) = i
	end do
	pointset(1:10,2) = 1
	pointset(11:20,2) = 2

	triangulation = bowyer_watson(pointset)

	open(newunit=fileunit, file='data/datafile1.txt')

	do i = 1, size(triangulation)
		do j = 1, 3
			write(unit=fileunit, fmt=*) triangulation(i) % vertices(j,:)
		end do
	end do

	close(unit=fileunit)

end program test
