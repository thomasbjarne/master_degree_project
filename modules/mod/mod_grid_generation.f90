module mod_grid_generation

	use mod_derived_types
	use iso_fortran_env, only: real32, real64
	implicit none

contains

	pure function delaunay_triangulation(P) result(triangles)

		real, intent(in), dimension(:,:) :: P
		type(triangle), dimension(size(P,1)/3) :: triangles
		type(triangle) :: supertriangle
		type(circle) :: circumcircle
		integer :: i, j, k, l, n



	end function delaunay_triangulation

	pure subroutine cell_centered_quadrilateral_grid(xl, xu, yl, yu, n, x, y, squares)

		real, intent(in) :: xl, xu, yl, yu
		integer, intent(in) :: n
		type(quadrilateral), intent(out), dimension(n,n) :: squares
		real, intent(out), dimension(0:n) :: x, y
		real, dimension(4,2) :: vertices
		real :: h1, h2
		integer :: i, j

		h1 = (xu-xl)/n
		h2 = (yu-yl)/n

		do i = 0,n
			x(i) = i * h1
			y(i) = i * h2
		end do

		do i = 1, n
			do j = 1, n
				vertices(1:2,1) = x(i-1)
				vertices(3:4,1) = x(i)
				vertices(1,2) = y(j-1)
				vertices(2,2) = y(j)
				vertices(3,2) = y(j)
				vertices(4,2) = y(j-1)
				squares(i,j) = quadrilateral_constructor([i,j], vertices)
			end do
		end do

	end subroutine cell_centered_quadrilateral_grid

end module mod_grid_generation