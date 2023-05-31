module mod_derived_types

	use iso_fortran_env, only: real32, real64
	implicit none

	type :: triangulation
		real(real64), dimension(:,:), allocatable 	:: points, barycenters, &
														distancematrix
		integer, dimension(:,:), allocatable 		:: connectivity_list,&
														edge_list, neighbours
		real(real64), dimension(:), allocatable		:: volumes, dual_volumes
		real(real64)								:: h
	end type triangulation
	
	type :: partition1d
		real(real64), dimension(:), allocatable 	:: points, midpoints,volumes
		integer, dimension(:,:), allocatable		:: intervals, neighbours
		real(real64), dimension(:,:), allocatable	:: distancematrix
		real(real64)								:: h
	end type partition1d
	
contains

!partition1d procedures

pure function partition1d_midpoints(points, intervals) result(centers)

	real(real64), intent(in), dimension(:)	:: points
	integer, intent(in), dimension(:,:)		:: intervals
	real(real64), dimension(:), allocatable :: centers
	integer									:: i, N
	
	N = size(intervals,1)
	allocate ( centers(N))
	do i = 1, N
		centers(i) = (points(intervals(i,1)) + points(intervals(i,2)))/2.
	end do

end function partition1d_midpoints

pure function partition1d_volumes(points, intervals) result(V)

	real(real64), intent(in), dimension(:)	:: points
	integer, intent(in), dimension(:,:)		:: intervals
	real(real64), dimension(:), allocatable :: V
	integer									:: i, N
	
	N = size(intervals,1)
	allocate(V(N))
	do i = 1, N
		V(i) = dabs(points(intervals(i,1)) - points(intervals(i,2)))
	end do

end function partition1d_volumes

pure function partition1d_char_length(volumes) result(h)

	real(real64), intent(in), dimension(:)  :: volumes
	real(real64)							:: h
	
	h = maxval(volumes)

end function partition1d_char_length

pure function partition1d_neighbours(intervals) result(res)
	
	integer, intent(in), dimension(:,:)		:: intervals
	integer, dimension(:,:), allocatable	:: res
	logical, dimension(:), allocatable		:: maskvec
	integer									:: N, i, j, l
	
	N = size(intervals, 1)
	allocate (res(N, 2), maskvec(N))
	res 	= 0
	maskvec = .true.
	do i = 1, N
		maskvec(i) = .false.
		do j = 1, 2
			do l = 1, N
				if ((intervals(i, j) == intervals(l, 1) .or. &
					intervals(i, j) == intervals(l, 2)) .and. &
					maskvec(l) .eqv. .true.) then
					res(i,j) = l
				end if
			end do
		end do
		maskvec = .true.
	end do

end function partition1d_neighbours

pure function partition1d_distances(midpoints) result (dmatrix)

	real(real64), intent(in), dimension(:)		:: midpoints
	real(real64), dimension(:,:), allocatable	:: dmatrix
	integer										:: i, j, N
	
	N = size(midpoints,1)
	allocate (dmatrix(N,N))
	do i = 1, N
		do j = 1, N
			dmatrix(i,j) = abs(midpoints(i) - midpoints(j))
		end do
	end do	

end function partition1d_distances

!triangulation procedures

subroutine triangle_fix_orientation(points, connectivity_list)

	real(real64), intent(in), dimension(:,:)		:: points
	integer, intent(inout), dimension(:,:)			:: connectivity_list
	real(real64)									:: cross_product
	integer											:: i, N, temp
	
	N = size(connectivity_list, 1)
	do i = 1, N
		cross_product = (points(connectivity_list(i,2),1)-&
						points(connectivity_list(i,1),1))*&
						(points(connectivity_list(i,3),2)-&
						points(connectivity_list(i,1),2)) -&
						(points(connectivity_list(i,2),2)-&
						points(connectivity_list(i,1),2))*&
						(points(connectivity_list(i,3),1)-&
						points(connectivity_list(i,1),1))
		if (cross_product < 0) then
			temp				   = connectivity_list(i,3)
			connectivity_list(i,3) = connectivity_list(i,2)
			connectivity_list(i,2) = temp
		end if
	end do

end subroutine

pure function triangle_char_length(points, connectivity_list) result(res)

	real(real64), intent(in), dimension(:,:) 		:: points
	integer, intent(in), dimension(:,:) 			:: connectivity_list
	real(real64), dimension(:), allocatable			:: h_vec
	real(real64)									:: res
	integer											:: N, i
	
	N = size(connectivity_list,1)
	allocate(h_vec(N))
	h_vec = 0
	do i = 1, N
		h_vec(i) = maxval([&
			norm2(points(connectivity_list(i,1),:) - &
				  points(connectivity_list(i,2),:)),&
			norm2(points(connectivity_list(i,2),:) - &
				  points(connectivity_list(i,3),:)),&
			norm2(points(connectivity_list(i,3),:) - &
				  points(connectivity_list(i,1),:))])
	end do
	res = maxval(h_vec)

end function triangle_char_length

pure function triangle_barycenters(points, connectivity_list) result(centers)

	real(real64), intent(in), dimension(:,:) 		:: points
	integer, intent(in), dimension(:,:) 			:: connectivity_list
	real(real64), dimension(:,:), allocatable 		:: centers
	
	allocate(centers(size(points,1),size(points,2)))
	centers = (points(connectivity_list(:,1),:) + &
				points(connectivity_list(:,2),:) + &
				points(connectivity_list(:,3),:))/3

end function triangle_barycenters

pure function triangle_volumes(points, connectivity_list) result(volumes)

	real(real64), intent(in), dimension(:,:)	:: points
	integer, intent(in), dimension(:,:) 		:: connectivity_list
	real(real64), dimension(:), allocatable		:: volumes
	integer										:: N, i
	
	N = size(connectivity_list,1)
	allocate (volumes(N))
	do i = 1, N
		volumes(i) = (1./2.)*abs(points(connectivity_list(i,1),1)*&
						(points(connectivity_list(i,2),2)-&
						points(connectivity_list(i,3),2))+&
						points(connectivity_list(i,2),1)*&
						(points(connectivity_list(i,3),2)-&
						points(connectivity_list(i,1),2))+&
						points(connectivity_list(i,3),1)*&
						(points(connectivity_list(i,1),2)-&
						points(connectivity_list(i,2),2)))
	end do

end function triangle_volumes


pure function triangle_neighbours(connectivity_list) result(neighbours)

	integer, intent(in), dimension(:,:)		:: connectivity_list
	integer, dimension(:,:), allocatable	:: neighbours
	logical, dimension(:), allocatable		:: maskvec
	integer									:: N, i, j, l
	
	N = size(connectivity_list,1)
	allocate(neighbours(N,3), maskvec(N))
	neighbours	= 0
	maskvec		= .true.
	do i = 1, N
		maskvec(i) = .false.
		do j = 1, 2
			do l = 1, N
				if ((connectivity_list(i,j)== &
					connectivity_list(l,1) .or. &
					connectivity_list(i,j)== &
					connectivity_list(l,2) .or. &
					connectivity_list(i,j)== &
					connectivity_list(l,3)) .and. &
					((connectivity_list(i,j+1)== &
					connectivity_list(l,1) .or. &
					connectivity_list(i,j+1)== &
					connectivity_list(l,2) .or. &
					connectivity_list(i,j+1)== &
					connectivity_list(l,3)) .and. &
					maskvec(l).eqv. .true.)) then
					neighbours(i,j) = l
				end if
			end do
		end do
		do l = 1, N
				if ((connectivity_list(i,3)== &
					connectivity_list(l,1) .or. &
					connectivity_list(i,3)== &
					connectivity_list(l,2) .or. &
					connectivity_list(i,3)== &
					connectivity_list(l,3)) .and. &
					((connectivity_list(i,1)== &
					connectivity_list(l,1) .or. &
					connectivity_list(i,1)== &
					connectivity_list(l,2) .or. &
					connectivity_list(i,1)== &
					connectivity_list(l,3)) .and. &
					maskvec(l).eqv. .true.)) then
					neighbours(i,3) = l
				end if
			end do
		maskvec = .true.
	end do

end function triangle_neighbours

pure function triangle_distances(barycenters) result(distancematrix)

	real(real64), intent(in), dimension(:,:)		:: barycenters
	real(real64), dimension(:,:), allocatable		:: distancematrix
	integer											:: i, j, N
	
	N = size(barycenters, 1)
	allocate(distancematrix(N,N))
	do i = 1, N
		do j = 1, N
			distancematrix(i,j) = &
					norm2(barycenters(i,:)-barycenters(j,:))
		end do
	end do

end function triangle_distances

pure function node_volumes(points, connectivity_list, volumes)&
														 result(dual_volumes)

	real(real64), intent(in), dimension(:,:)	:: points
	integer, intent(in), dimension(:,:) 		:: connectivity_list
	real(real64), intent(in), dimension(:)		:: volumes
	real(real64), dimension(:), allocatable		:: dual_volumes
	integer										:: N, i, j
	
	N = size(points,1)
	allocate (dual_volumes(N))
	dual_volumes = 0
	do i = 1, N
		do j = 1, size(connectivity_list,1)
			if (i == connectivity_list(j,1) .or. &
				i == connectivity_list(j,2) .or. &
				i == connectivity_list(j,3)) then
				dual_volumes(i) = dual_volumes(i) + volumes(j)/3
			end if
		end do
	end do

end function node_volumes

end module mod_derived_types
