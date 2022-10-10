module mod_grid_generation

	use mod_derived_types
	use iso_fortran_env, only: real32, real64
	implicit none

contains

	pure function bowyer_watson(P) result(triangulation)

		real, intent(in), dimension(:,:) :: P
		type(triangle), dimension(size(P,1)) :: triangulation
		type(triangle), dimension(size(P,1)*2) :: bad_triangles
		type(triangle), dimension(size(P,1)**2) :: all_triangles
		type(triangle) :: supertriangle
		type(circle), dimension(size(all_triangles)) :: all_circumcircles
		real, dimension(3, 2) :: supertriangle_vertices, triangle_vertices
		integer :: i, j, k, l, n, g

		g = size(P,1)
		triangulation = triangle(0,0)
		bad_triangles = triangle(0,0)
		all_triangles = triangle(0,0)
		supertriangle_vertices(1,:) = [minval(P(:,1))-maxval(P(:,1)), minval(P(:,2))-1.]
		supertriangle_vertices(2,:) = [2*maxval(P(:,1)), minval(P(:,2))-1.]
		supertriangle_vertices(3,:) = [(maxval(P(:,1))+minval(P(:,1)))/2., 2*maxval(P(:,2))+1.]
		supertriangle = triangle(index=g+1, vertices=supertriangle_vertices)

		! Create all possible triangles from pointset and supertriangle_vertices
		do i = 1, g
			do j = 1, 3
				triangle_vertices(1,:) = P(i,:)
				if (j == 3) then
					triangle_vertices(2,:) = supertriangle_vertices(3,:)
					triangle_vertices(3,:) = supertriangle_vertices(1,:)
				else
					triangle_vertices(2,:) = supertriangle_vertices(j,:)
					triangle_vertices(3,:) = supertriangle_vertices(j+1,:)
				end if
				all_triangles(j+i) = triangle(index=j+i, &
				vertices=triangle_vertices)
				if (i > 1) then
					do k = i-1, 1, -1
						triangle_vertices(1,:) = P(i,:)
						triangle_vertices(2,:) = P(i-k,:)
						triangle_vertices(3,:) = supertriangle_vertices(j,:)
						all_triangles(j+i+k) = triangle(index=j+i+k, &
						vertices = triangle_vertices)
					end do
				end if
			end do
		end do

		! Create all circumcircles for all_triangles
		do i = 1, size(all_triangles)
			all_circumcircles(i) = circle( &
				radius=norm2(all_triangles(i) % vertices(1,:) - all_triangles(i) % circumcenter(1,:)), &
				center=all_triangles(i) % circumcenter)
		end do

		! Add triangles with vertices inside some circumcircle to bad_triangles
		do i = 1, size(all_triangles)
			do j = 1, 3
				do k = 1, size(all_circumcircles)
					if (norm2(all_triangles(i) % vertices(j,:) - all_circumcircles(k) % center(1,:)) &
					< all_circumcircles(k) % radius) then
						bad_triangles(i) = all_triangles(i)
					end if
				end do
			end do
		end do

		triangulation = all_triangles(1:size(triangulation))

	end function bowyer_watson

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