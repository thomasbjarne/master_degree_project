module mod_grid_generation

	use mod_derived_types
	use mod_data_manip
	use iso_fortran_env, only: real32, real64
	implicit none

contains

	function bowyer_watson(P) result(triangulation)

		real, intent(in), dimension(:,:) :: P
		type(triangle), dimension(:), allocatable :: triangulation
		type(triangle), dimension(:), allocatable :: bad_triangles
		type(triangle) :: supertriangle
		real, dimension(3, 2) :: supertriangle_vertices, triangle_vertices
		type(line), dimension(:), allocatable :: new_edge_collection, bt_edge_collection
		type(circle), dimension(:), allocatable :: circumcircles
		type(circle) :: supercircle, B
		integer :: i, j, k, l, n, g

		allocate(triangulation(0), circumcircles(0))

		g = size(P,1)
		supertriangle_vertices(1,:) = &
		[-100*(minval(P(:,1))-maxval(P(:,1))), minval(P(:,2))-1.]
		supertriangle_vertices(2,:) = &
		[100*maxval(P(:,1)), minval(P(:,2))-1.]
		supertriangle_vertices(3,:) = &
		[(maxval(P(:,1))+minval(P(:,1)))/2., 100*maxval(P(:,2))+1.]
		supertriangle = triangle(vertices=supertriangle_vertices)
		supercircle = circle(&
		radius=norm2(supertriangle_vertices(1,:)-supertriangle%circumcenter(1,:)),&
		center=supertriangle%circumcenter)
		call union(supertriangle, triangulation)
		call union(supercircle, circumcircles)

		! ref: pseudocode on wikipedia article of bowyer-watson
		do i = 1, g
			allocate( bad_triangles(0), new_edge_collection(0), bt_edge_collection(0))
			do n = 1, size(triangulation)
				B = circle(&
				radius=norm2(triangulation(n)%vertices(1,:)-triangulation(n)%circumcenter(1,:)),&
				center=triangulation(n)%circumcenter)
				if (point_in_circle(P(i,:), B)) then
                    call union(triangulation(n), bad_triangles)
                end if
			end do
			do j = 1, size(bad_triangles)
				do k = 1, 3
					call union(bad_triangles(j)%edges(k), bt_edge_collection)
				end do
			end do
			do j = 1, size(bt_edge_collection)
				if (count_appearance(bt_edge_collection(j),bt_edge_collection) == 1) then
					call union(bt_edge_collection(j), new_edge_collection)
				end if
			end do
			do j = 1, size(bad_triangles)
				!remove bad triangles from triangulation
				do n = 1, size(triangulation)
					if (equals_triangle(bad_triangles(j), triangulation(n))) then
						triangulation(n) = triangle(vertices=0)
					end if
				end do
			end do
			do k = 1, size(new_edge_collection)
				triangle_vertices(1,:) = P(i,:)
				triangle_vertices(2,:) = new_edge_collection(k) % p1(1,:)
				triangle_vertices(3,:) = new_edge_collection(k) % p2(1,:)
				call union(triangle(vertices=triangle_vertices), triangulation)
			end do
			deallocate(bad_triangles, new_edge_collection, bt_edge_collection)
		end do

		!remove all triangles with one vertex equal to a vertex of supertriangle
		do n = 1, size(triangulation)
			do j = 1, 3
				do l = 1, 3
					if (norm2(triangulation(n)%vertices(j,:) - supertriangle_vertices(l,:))==0)then
						triangulation(n) = triangle(0,0)
					end if
				end do
			end do
		end do

	end function bowyer_watson

	pure subroutine cell_centered_quadrilateral(xl, xu, yl, yu, n, x, y, squares)

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

	end subroutine cell_centered_quadrilateral

end module mod_grid_generation
