module mod_reconstruction

	use iso_fortran_env,	only: real32, real64
	use mod_linalg
	use mod_moments
	use mod_derived_types
	implicit none
	
	interface recon_matrix_2d
		module procedure :: recon_matrix_2d_kind32, recon_matrix_2d_kind64
	end interface recon_matrix_2d
	
	interface recon_coeffs
		module procedure :: recon_coeffs_kind32, recon_coeffs_kind64
	end interface recon_coeffs
	
contains

!1D procedures:

pure function stencil_select1d(dof, partition, indexnum) result(stencil)

	type(partition1d), intent(in)			:: partition
	integer, intent(in)						:: dof, indexnum
	integer, dimension(:), allocatable		:: stencil
	real(real64), dimension(:), allocatable :: distances
	integer									:: j, N
	logical, dimension(:), allocatable		:: mask_vec
	
	N = size(partition%intervals, 1)
	allocate (stencil(dof+4), distances(N), mask_vec(N))
	
	distances = partition%distancematrix(indexnum,:)
	
	!choose the dof closest intervals as reconstruction stencil
	mask_vec = .true.
	do j = 1, (dof + 4)
		stencil(j) = minloc(distances, dim=1, mask=mask_vec)
		mask_vec(minloc(distances, mask=mask_vec)) = .false.
	end do
	

end function stencil_select1d

pure function recon_matrix_1d(stencil,moments1,moments2,partition,indexnum)&
																	 result(Mi)

	type(partition1d), intent(in)				:: partition
	integer, intent(in), dimension(:)			:: stencil
	integer, intent(in)							:: indexnum
	real(real64), intent(in), dimension(:)		:: moments1
	real(real64), intent(in), dimension(:,:)	:: moments2
	real(real64), dimension(:,:), allocatable	:: Mi
	real(real64), dimension(:), allocatable 	:: distances
	integer										:: i, j, N
	
	N = size(stencil,1)
	allocate (distances(N),Mi(N, N-4))
	distances = partition%distancematrix(indexnum,:)
	
	do j = 1, N-4
		Mi(1, j) = moments1(j)
	end do
	do i = 2, N
		do j = 1, N-4
			if (distances(stencil(i)) /=0) then
				Mi(i,j) = moments2(stencil(i),j)/(abs(distances(stencil(i))))
			else
				Mi(i,j) = moments2(stencil(i), j)
			end if
		end do
	end do
	

end function recon_matrix_1d

pure function recon_coeffs_1d(stencil, vec_u, Mi, partition, indexnum)&
																 result(alpha)

	type(partition1d), intent(in)					:: partition
	integer, intent(in), dimension(:)				:: stencil
	integer, intent(in)								:: indexnum
	real(real64), intent(in), dimension(:)			:: vec_u
	real(real64), intent(in), dimension(:,:)		:: Mi
	real(real64), dimension(:), allocatable			:: alpha, b, distances
	integer											:: sigma, j
	
	sigma = size(stencil,1)
	allocate(alpha(sigma-4), b(sigma), distances(sigma))
	
	distances = partition%distancematrix(indexnum,:)
	
	b(1) = vec_u(stencil(1))
	do j = 2, sigma
		b(j) = vec_u(stencil(j))/(abs(distances(stencil(j))))
	end do
	alpha = constrained_least_squares(Mi, b)

end function recon_coeffs_1d


!2D procedures:

pure function stencil_select(dof, grid, indexnum) result(stencil)

	type(triangulation), intent(in)				:: grid
	integer, intent(in)							:: dof, indexnum
	integer, dimension(:), allocatable			:: stencil
	real(real64), dimension(:), allocatable		:: distances
	integer										:: k, N
	logical, dimension(:), allocatable 			:: mask_vec
	
	N = size(grid%connectivity_list,1)
	allocate (stencil(dof+4), distances(N), mask_vec(N))
	
	distances = grid%distancematrix(indexnum,:)
	
	! choose the dof closest triangles as reconstruction stencil, including
	! triangle(indexnum).
	mask_vec = .true.
	do k = 1, (dof + 4)
		stencil(k) = minloc(distances, dim=1, mask=mask_vec)
		mask_vec(minloc(distances, mask=mask_vec)) = .false.
	end do
	
end function stencil_select

pure function recon_matrix_2d_kind32(stencil,moments1,moments2,grid,indexnum)&
																	result(Mi)

	type(triangulation), intent(in)					:: grid
	integer, intent(in), dimension(:)				:: stencil
	integer, intent(in)								:: indexnum
	real(real32), intent(in), dimension(:,:)		:: moments2
	real(real32), intent(in), dimension(:)			:: moments1
	real(real64), dimension(:), allocatable			:: distances
	real(real32), dimension(:,:), allocatable		:: Mi
	integer											:: sigma, i, j, N
	
	sigma = size(stencil,1)
	N	  = size(grid%connectivity_list,1)
	allocate ( Mi(sigma, sigma-4), distances(N))
	distances = grid%distancematrix(indexnum,:)
	do j = 1, sigma-4
			Mi(1,j) = moments1(j)
		end do
	do i = 2, sigma
		do j = 1, sigma-4
			Mi(i,j) = moments2(stencil(i), j)/(abs(distances(stencil(i))))
		end do
	end do

end function recon_matrix_2d_kind32

pure function recon_matrix_2d_kind64(stencil, moments1, moments2,grid,indexnum)&
																	 result(Mi)

	type(triangulation), intent(in)					:: grid
	integer, intent(in), dimension(:)				:: stencil
	integer, intent(in)								:: indexnum
	real(real64), intent(in), dimension(:,:)		:: moments2
	real(real64), intent(in), dimension(:)			:: moments1
	real(real64), dimension(:), allocatable			:: distances
	real(real64), dimension(:,:), allocatable		:: Mi
	integer											:: sigma, i, j, N
	
	sigma = size(stencil,1)
	N	  = size(grid%connectivity_list,1)
	allocate ( Mi(sigma, sigma-4), distances(N))
	distances = grid%distancematrix(indexnum,:)
	do j = 1, sigma-4
			Mi(1,j) = moments1(j)
		end do
	do i = 2, sigma
		do j = 1, sigma-4
			Mi(i,j) = moments2(stencil(i), j)/(abs(distances(stencil(i))))
		end do
	end do

end function recon_matrix_2d_kind64

!procedures to use each timestep

pure function recon_coeffs_kind32(stencil, vec_u, Mi, grid, indexnum)&
																 result(alpha)

	type(triangulation), intent(in)					:: grid
	integer, intent(in), dimension(:)				:: stencil
	integer, intent(in)								:: indexnum
	real(real32), intent(in), dimension(:)			:: vec_u
	real(real32), intent(in), dimension(:,:)		:: Mi
	real(real32), dimension(:), allocatable			:: alpha, b, distances
	integer											:: sigma, j, N
	
	N 	  = size(grid%connectivity_list,1)
	sigma = size(stencil,1)
	allocate(alpha(sigma-4), b(sigma), distances(N))
	
	distances = grid%distancematrix(indexnum,:)
	b(1) = vec_u(stencil(1))
	do j = 2, sigma
		b(j) = vec_u(stencil(j))/(abs(distances(stencil(j))))
	end do
	alpha = constrained_least_squares(Mi, b)

end function recon_coeffs_kind32

pure function recon_coeffs_kind64(stencil, vec_u, Mi, grid, indexnum)&
																 result(alpha)

	type(triangulation), intent(in)					:: grid
	integer, intent(in), dimension(:)				:: stencil
	integer, intent(in)								:: indexnum
	real(real64), intent(in), dimension(:)			:: vec_u
	real(real64), intent(in), dimension(:,:)		:: Mi
	real(real64), dimension(:), allocatable			:: alpha, b, distances
	integer											:: sigma, j, N
	
	N 	  = size(grid%connectivity_list,1)
	sigma = size(stencil,1)
	allocate(alpha(sigma-4), b(sigma), distances(N))
	
	distances = grid%distancematrix(indexnum,:)
	b(1) = vec_u(stencil(1))
	do j = 2, sigma
		b(j) = vec_u(stencil(j))/(abs(distances(stencil(j))))
	end do
	alpha = constrained_least_squares(Mi, b)

end function recon_coeffs_kind64

end module mod_reconstruction
