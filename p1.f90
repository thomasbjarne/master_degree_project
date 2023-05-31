program p1
	
	use iso_fortran_env, only: real32, real64
	use mod_quadrature
	use mod_io
	use mod_moments
	use mod_derived_types
	use mod_reconstruction
	use mod_scheme
	implicit none
	
	!u_t + u_x 	= 0
	!u(x, 0) 	= sin(2pi x)
	!u(0, t)	= g(t)	= sin(-2pi t)
	
	type(partition1d)							:: partition
	integer										:: N, kexact, rdof, i, j, &
													timesteps
	real(real64)								:: error, dt, h, pi, &
													t, tmin, tmax, a, b
	real(real64), dimension(:), allocatable		:: u, exactsol, a_1, a_2, a_3,&
													a_4, x, xi
	real(real64), dimension(:,:), allocatable	:: primarymoments, V
	real(real64), dimension(:,:,:), allocatable	:: rmatrices, secondarymoments
	integer, dimension(:,:), allocatable		:: rstencils
	integer, dimension(:), allocatable			:: degvec
	
	!set up grid
	N = 400
	allocate (x(N+1), xi(N+1), partition%points(N+1), partition%intervals(N,2))
	
	!structured regular grid
	do j = 1, N+1
		x(j) = (j-1)*(1./N)
	end do
	
	!structured irregular grid
	!do j = 1, N+1
	!	xi(j) 	= (j-1)*(1./N)
	!	x(j)	= (exp(xi(j))-1)/(exp(1.)-1)
	!end do
	
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
	
	!write grid to datafiles for debugging
	call write_to_file(partition%points, 			1)
	call write_to_file(partition%intervals, 		2)
	call write_to_file(partition%midpoints, 		3)
	call write_to_file(partition%volumes, 			4)
	print*, 'h = ', partition%h
	call write_to_file(partition%neighbours, 		5)
	call write_to_file(partition%distancematrix, 	6)
	
	!set up variables
	kexact	 	= 3
	rdof		= kexact + 1
	N 			= size(partition%intervals, 1)
	h			= partition%h
	allocate (u(N),rstencils(N,rdof+4),rmatrices(N,rdof+4,rdof))
	allocate (primarymoments(N,11),degvec(11),secondarymoments(N,N,11),V(N,N))
	pi 			= 4.D0*DATAN(1.D0)
	tmax 		= 1
	tmin 		= 0
	dt 			= h/2.
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

	!write reconstruction tools to datafiles for debugging
	call write_to_file(rstencils, 7)
	call write_to_file(primarymoments, 8)
	call write_to_file(secondarymoments(1,:,:), 9)
	call write_to_file(rmatrices(1,:,:), 10)

	!set initial data
	do j = 1, N
		a 		= partition%points(partition%intervals(j,1))
		b 		= partition%points(partition%intervals(j,2))
		u(j) 	= transformed_gauss_legendre(a, b, 5, initialfunc)/V(j,j)
	end do
	allocate(exactsol(N))
	exactsol = u
	call write_to_file(exactsol,11)
	
	allocate(a_1(N), a_2(N), a_3(N), a_4(N))
	
	!RK4 time stepping scheme
	do i = 1, timesteps
		t = (i-1)*dt
		
		a_1 = scheme_1d(u, 			   t, 		  partition,rstencils,rmatrices)
		a_2 = scheme_1d(u+a_1*(dt/2.), t+(dt/2.), partition,rstencils,rmatrices)
		a_3 = scheme_1d(u+a_2*(dt/2.), t+(dt/2.), partition,rstencils,rmatrices)
		a_4 = scheme_1d(u+a_3*dt, 	   t + dt,    partition,rstencils,rmatrices)
		
		u = u + (dt/6.)*(a_1 + 2*a_2 + 2*a_3 + a_4) ! u at t + dt
	end do
	
	
	!third order SSP-RK time stepping scheme	
	!do i = 1, timesteps
	!	t = (i-1)*dt
	!	
	!	a_1 = u + dt*scheme_1d(u, t, partition, rstencils, rmatrices)
	!	a_2 = (3./4.)*u + (1./4.)*a_1 + &
	!			(1./4.)*dt*scheme_1d(a_1, t, partition,rstencils,rmatrices)
	!	a_3 = (1./3.)*u + (2./3.)*a_2 + &
	!			(2./3.)*dt*scheme_1d(a_2, t, partition,rstencils,rmatrices)
	!			
	!	u = a_3
	!end do
	
	!convergence analysis
	call write_to_file(u, 12)
	error = sqrt(dot_product(exactsol-u, matmul(V, exactsol-u)))
	print*, 'Error at t = ', t+dt, 'is', error
	
contains

pure function initialfunc(x) result(res)

	real(real64), intent(in)	:: x
	real(real64)				:: res
	
	res = sin(2*pi*x)
	
end function initialfunc


end program
