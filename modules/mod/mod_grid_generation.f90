module mod_grid_generation

	use mod_derived_types
	use iso_fortran_env, only: real32, real64
	implicit none

contains

	pure function bowyer_watson(P) result(triangulation)

		real, intent(in), dimension(:,:) :: P
		type(triangle), dimension(size(P,1)*10) :: triangulation
		type(triangle), dimension(size(P,1)*10) :: bad_triangles
		type(triangle) :: supertriangle
		real, dimension(3, 2) :: supertriangle_vertices, triangle_vertices
		type(line), dimension(size(bad_triangles)) :: edge_collection
		integer :: i, j, k, l, n, g

		g = size(P,1)
		supertriangle_vertices(1,:) = [minval(P(:,1))-maxval(P(:,1)), minval(P(:,2))-1.]
		supertriangle_vertices(2,:) = [2*maxval(P(:,1)), minval(P(:,2))-1.]
		supertriangle_vertices(3,:) = [(maxval(P(:,1))+minval(P(:,1)))/2., 2*maxval(P(:,2))+1.]
		supertriangle = triangle(index=0, vertices=supertriangle_vertices)
		triangulation(1) = supertriangle

		! ref: pseudocode on wikipedia article of bowyer-watson
		do i = 1, g
			bad_triangles = triangle(index=0,vertices=0)
			do n = 1, size(triangulation)
				if (&
				norm2(P(i,:) - triangulation(n) % circumcenter(1,:)) &
				< norm2(triangulation(n) % vertices(1,:) - triangulation(n) % circumcenter(1,:))) then
                    bad_triangles(n) = triangulation(n)
                end if
			end do
			edge_collection = line(p1=0,p2=0)
			do j = 1, size(bad_triangles)
				do k = 1, 3
					do l = 1, size(bad_triangles)
						if (norm2((bad_triangles(j)%edges(k)%p1) - (bad_triangles(l)%edges(k)%p1)) &
						+ norm2((bad_triangles(j)%edges(k)%p2) - (bad_triangles(l)%edges(k)%p2))>0) then
							edge_collection(j) = bad_triangles(j) % edges(k)
						end if
					end do
				end do

				!remove bad triangles from triangulation
				do n = 2, size(triangulation)
					if (norm2(bad_triangles(j) % vertices - triangulation(n) % vertices)==0) then
						triangulation(n) = triangle(0,0)
					end if
				end do
			end do
			do k = 1, size(edge_collection)
				if (norm2(edge_collection(k) % p1 + edge_collection(k) % p2) > 0) then
					triangle_vertices(1,:) = P(i,:)
					triangle_vertices(2,:) = edge_collection(k) % p1(1,:)
					triangle_vertices(3,:) = edge_collection(k) % p2(1,:)
					triangulation(i + 3 + k) = triangle(index=i+3+k, vertices=triangle_vertices)
				end if
			end do
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

	end subroutine cell_centered_quadrilateral_grid

end module mod_grid_generation
