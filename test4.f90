program test4

	use iso_fortran_env, only: real32, real64
	use mod_quadrature
	use mod_io
	use mod_moments
	use mod_derived_types
	use mod_reconstruction
	use mod_scheme
	implicit none
	
	!u_t + u_x + u_y = 0
	!u(x,y,0) = f(x,y) 		= sin(2pi(x/2 + y/2))
	!u(0,y,t) = g_1(y,t)	= sin(2pi(y/2 - t))
	!u(x,0,t) = g_2(x,t)	= sin(2pi(x/2 - t))
	type(triangulation)							:: grid
	integer										:: N, kexact, rdof, i, j, l, &
													timesteps
	real(real64)								:: error, dt, h, pi, &
													t, tmin, tmax, xi, yi
	real(real64), dimension(:), allocatable		:: u, exactsol, a_1, a_2, a_3,&
													a_4
	real(real64), dimension(:,:), allocatable	:: u_0_coeffs, primarymoments, &
													V, R_coeffs
	real(real64), dimension(:,:,:), allocatable	:: rmatrices, secondarymoments
	integer, dimension(:,:), allocatable		:: rstencils, degvec
	
	!set up grid
	call load_grid(grid)
	
	!write grid to datafiles for debugging
	call write_to_file(grid%points, 			1)
	call write_to_file(grid%connectivity_list, 	2)
	call write_to_file(grid%barycenters, 		3)
	call write_to_file(grid%volumes, 			4)
	print*, 'h = ', grid%h
	call write_to_file(grid%neighbours, 		5)
	call write_to_file(grid%distancematrix, 	6)
	
	!set up variables
	kexact	 	= 0
	rdof		= (kexact+2)*(kexact+1)/2.
	N 			= size(grid%connectivity_list, 1)
	h			= grid%h
	allocate(u(N),rstencils(N,rdof+4),rmatrices(N,rdof+4,rdof))
	allocate(u_0_coeffs(N,15))
	allocate(primarymoments(N,15),degvec(15,2),secondarymoments(N,N,15),V(N,N))
	pi 			= 4.D0*DATAN(1.D0)
	tmax 		= 1
	tmin 		= 0
	dt 			= h/4.
	timesteps 	= floor((tmax-tmin)/dt)
	dt			= (tmax-tmin)/timesteps
	degvec(:,1) = [0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0]
	degvec(:,2) = [0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4]
	
	!prepare reconstruction tools
	do j = 1, 15
		primarymoments(:,j) = calculate_moments(degvec(j,1),degvec(j,2),grid)
	end do
	do i = 1, N
		do j = 1, 15
			secondarymoments(i,:,j) = derived_moments(degvec(j,1),degvec(j,2),&
									grid, i, primarymoments(:,:))
		end do
	end do
	do i = 1, N
		rstencils(i,:) 	 = stencil_select(rdof, grid, i)
		rmatrices(i,:,:) = recon_matrix_2d(rstencils(i,:),&
								 primarymoments(i,:), secondarymoments(i, :,:),&
								 grid, i)
		V(i,i) 		 	 = grid%volumes(i)
	end do
	
	!write reconstruction tools to datafiles for debugging
	call write_to_file(rstencils, 7)
	call write_to_file(primarymoments, 8)
	call write_to_file(secondarymoments(1,:,:), 9)
	call write_to_file(rmatrices(1,:,:), 10)
	
	!set initial data
	do i = 1, N
		xi = grid%barycenters(i,1)
		yi = grid%barycenters(i,2)
		u_0_coeffs(i, :) = [f(xi,yi),f10(xi,yi),f01(xi,yi),f20(xi,yi)/2.,&
							f11(xi,yi),f02(xi,yi)/2.,f30(xi,yi)/6.,&
							f21(xi,yi)/2.,f12(xi,yi)/2.,f03(xi,yi)/6.,&
							f40(xi,yi)/24.,f31(xi,yi)/6.,f22(xi,yi)/4.,&
							f13(xi,yi)/6.,f04(xi,yi)/24.]
		u(i) = dot_product(primarymoments(i,1:rdof), u_0_coeffs(i,1:rdof))
	end do
	allocate (exactsol(N))
	exactsol = u
	call write_to_file(exactsol, 11)
	

	
	!test R_coeffs at t = 0
	allocate (R_coeffs(N, rdof))
	do j = 1, N
		R_coeffs(j,:) = recon_coeffs(rstencils(j,:), u, rmatrices(j,:,:),&
									grid, j)
	end do
	call write_to_file(R_coeffs, 13)
	
	!test L1 error
	error = 0
	do i = 1, N
		error = error + &
				abs(f(grid%points(grid%connectivity_list(i,1),1), &
					grid%points(grid%connectivity_list(i,1),2)) -&
					 Ri(grid%points(grid%connectivity_list(i,1),1), &
					 grid%points(grid%connectivity_list(i,1),2))) &
				+ &
				abs(f(grid%points(grid%connectivity_list(i,2),1), &
					grid%points(grid%connectivity_list(i,2),2)) -&
					 Ri(grid%points(grid%connectivity_list(i,2),1), &
					 grid%points(grid%connectivity_list(i,2),2))) &
				+ & 
				abs(f(grid%points(grid%connectivity_list(i,3),1), &
					grid%points(grid%connectivity_list(i,3),2)) -&
					 Ri(grid%points(grid%connectivity_list(i,3),1), &
					 grid%points(grid%connectivity_list(i,3),2)))
	end do
	print*, 'L1 error = ', error
	
	!RK4 time stepping scheme
	!allocate(a_1(N), a_2(N), a_3(N), a_4(N))
	!do i = 1, timesteps
	!	t = (i-1)*dt
	!	
	!	a_1 = scheme(u, 			  t, 		   grid, rstencils, rmatrices)
	!	a_2 = scheme(u + a_1*(dt/2.), t + (dt/2.), grid, rstencils, rmatrices)
	!	a_3 = scheme(u + a_2*(dt/2.), t + (dt/2.), grid, rstencils, rmatrices)
	!	a_4 = scheme(u + a_3*dt, 	  t + dt, 	   grid, rstencils, rmatrices)
	!	
	!	u = u + (dt/6.)*(a_1 + 2*a_2 + 2*a_3 + a_4) ! u at t +dt
	!end do
	!call write_to_file(u, 12)
	
	!convergence analysis
	!error = sqrt(dot_product(exactsol-u, matmul(V, exactsol-u)))
	!print*, 'Error at t = ', t + dt, 'is', error
		
contains

pure function Ri(x, y) result(res)

	real(real64), intent(in)	:: x, y
	real(real64)				:: res
	integer						:: m, sigma
	
	sigma 	= size(R_coeffs, 2)
	res = 0
	do m = 1, sigma
		res = res + &
				R_coeffs(i,m)*((x-grid%barycenters(i,1))**degvec(m,1))&
				*((y-grid%barycenters(i,2))**degvec(m,2))
	end do

end function Ri

pure function f(x,y) result(iv)
	
	real(real64), intent(in)	:: x, y
	real(real64)				:: iv
	
	iv = x**3 + y

end function f

pure function f10(x,y) result(res)

	real(real64), intent(in)	:: x, y
	real(real64)				:: res
	
	res = 3*(x**2)

end function f10

pure function f01(x,y) result(res)

	real(real64), intent(in)	:: x, y
	real(real64)				:: res
	
	res = 1

end function f01

pure function f20(x,y) result(res)

	real(real64), intent(in)	:: x, y
	real(real64)				:: res
	
	res = 6*x

end function f20

pure function f11(x,y) result(res)

	real(real64), intent(in)	:: x, y
	real(real64)				:: res
	
	res = 0

end function f11

pure function f02(x,y) result(res)

	real(real64), intent(in)	:: x, y
	real(real64)				:: res
	
	res = 0

end function f02

pure function f30(x,y) result(res)

	real(real64), intent(in)	:: x, y
	real(real64)				:: res
	
	res = 6

end function f30

pure function f21(x,y) result(res)

	real(real64), intent(in)	:: x, y
	real(real64)				:: res
	
	res = 0

end function f21

pure function f12(x,y) result(res)

	real(real64), intent(in)	:: x, y
	real(real64)				:: res
	
	res = 0

end function f12

pure function f03(x,y) result(res)

	real(real64), intent(in)	:: x, y
	real(real64)				:: res
	
	res = 0

end function f03

pure function f40(x,y) result(res)

	real(real64), intent(in)	:: x, y
	real(real64)				:: res
	
	res = 0

end function f40

pure function f31(x,y) result(res)

	real(real64), intent(in)	:: x, y
	real(real64)				:: res
	
	res = 0

end function f31

pure function f22(x,y) result(res)

	real(real64), intent(in)	:: x, y
	real(real64)				:: res
	
	res = 0

end function f22

pure function f13(x,y) result(res)

	real(real64), intent(in)	:: x, y
	real(real64)				:: res
	
	res = 0

end function f13

pure function f04(x,y) result(res)

	real(real64), intent(in)	:: x, y
	real(real64)				:: res
	
	res = 0

end function f04
	
end program test4
