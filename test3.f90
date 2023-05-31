program test3

	use iso_fortran_env, only: real32, real64
	use mod_quadrature
	use mod_io
	use mod_linalg
	use mod_moments
	use mod_reconstruction
	implicit none
	
	type(partition1d)							:: partition
	integer										:: N, kexact, rdof, i, j, &
													timesteps
	real(real64)								:: dt, h, pi, &
													t, tmin, tmax, a, b,L1err,&
													L2err
	real(real64), dimension(:), allocatable		:: u, x
	real(real64), dimension(:,:), allocatable	:: primarymoments, V
	real(real64), dimension(:,:), allocatable	:: R_coeffs, Vinv
	real(real64), dimension(:,:,:), allocatable	:: rmatrices, secondarymoments
	integer, dimension(:,:), allocatable		:: rstencils
	integer, dimension(:), allocatable			:: degvec
	
	!set up grid
	N = 100
	allocate (x(N+1), partition%points(N+1), partition%intervals(N,2))
	do j = 1, N+1
		x(j) = (j-1)*(1./N)
	end do
	partition%points = x
	do j = 1, N
		partition%intervals(j,1) = j
		partition%intervals(j,2) = j+1
	end do
	partition%midpoints 	 = partition1d_midpoints(partition%points, &
												partition%intervals)
	partition%volumes 		 = partition1d_volumes(partition%points, &
											  partition%intervals)
	partition%h 			 = partition1d_char_length(partition%volumes)
	partition%neighbours	 = partition1d_neighbours(partition%intervals)
	partition%distancematrix = partition1d_distances(partition%midpoints)
	call write_to_file(partition%intervals, 1)
	call write_to_file(partition%points, 2)
	call write_to_file(partition%volumes, 3)
	
	!set up variables
	kexact	 	= 6
	rdof		= kexact + 1
	N 			= size(partition%intervals, 1)
	h			= partition%h
	allocate (u(N),rstencils(N,rdof+4),rmatrices(N,rdof+4,rdof))
	allocate (primarymoments(N,11),degvec(11),secondarymoments(N,N,11),V(N,N))
	pi 			= 4.D0*DATAN(1.D0)
	tmax 		= 1
	tmin 		= 0
	dt 			= h
	timesteps 	= floor((tmax-tmin)/dt)
	dt			= (tmax-tmin)/timesteps
	degvec 		= [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
	
	!prepare reconstruction tools
	do j = 1, 11
		primarymoments(:,j) = primarymoments1d(degvec(j), partition)
	end do
	do i = 1, N
		do j = 1, 11
			secondarymoments(i,:,j) = secondarymoments1d(degvec(j),partition,i)
		end do
	end do
	do i = 1, N
		rstencils(i,:) 	 = stencil_select1d(rdof, partition, i)
		rmatrices(i,:,:) = recon_matrix_1d(rstencils(i,:),&
								 primarymoments(i,:), secondarymoments(i, :,:),&
								 partition, i)
		V(i,i) 		 	 = partition%volumes(i)
	end do
	call write_to_file(rmatrices(2,:,:), 4)
	

	!set initial data
	do j = 1, N
		a 		= partition%points(partition%intervals(j,1))
		b 		= partition%points(partition%intervals(j,2))
		u(j) 	= transformed_gauss_legendre(a, b, 5, initialfunc)/V(j,j)
	end do
	call write_to_file(u, 5)
	
	!test reconstruction
	allocate(R_coeffs(N, rdof))
	do j = 1, N
		R_coeffs(j,:) = recon_coeffs_1d(rstencils(j,:), u, rmatrices(j,:,:),&
										partition, j)
	end do
	call write_to_file(R_coeffs, 6)
	L1err = 0
	do i = 1, N
		L1err = L1err + &
				abs(initialfunc(partition%points(i)) -&
					 Ri(partition%points(i)))
	end do
	print*, 'L1 error = ', L1err
	L2err = 0
	
	
	
contains

pure function initialfunc(x) result(res)

	real(real64), intent(in)	:: x
	real(real64)				:: res
	
	res = x
	
end function initialfunc

pure function Ri(x) result (res)

	real(real64), intent(in)	:: x
	real(real64)				:: res
	integer						:: m, sigma
	
	sigma 	= size(R_coeffs, 2)
	res = 0
	do m = 1, sigma
		res = res + R_coeffs(i,m)*&
					((x-partition%midpoints(i))**&
					degvec(m))
	end do

end function Ri


end program test3
