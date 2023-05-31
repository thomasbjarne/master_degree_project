module mod_scheme

	use iso_fortran_env, only: real32, real64
	use mod_derived_types
	use mod_quadrature
	use mod_reconstruction
	implicit none
	
contains

!1D scheme

pure function scheme_1d(u, t, partition, stencils, matrices) result(res)

	type(partition1d), intent(in)				:: partition
	real(real64), intent(in), dimension(:,:,:)	:: matrices
	real(real64), intent(in), dimension(:)		:: u
	real(real64), intent(in)					:: t
	integer, intent(in), dimension(:,:)			:: stencils
	integer, dimension(:), allocatable			:: degvec
	real(real64), dimension(:,:), allocatable	:: R_coeffs, V, Vinv
	real(real64), dimension(:), allocatable 	:: res
	real(real64)								:: pi
	integer										:: i, j, N, rdof, fluxtype
	
	pi 		= 4.D0*DATAN(1.D0)
	!set up reconstructions
	allocate (degvec(11))
	degvec = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
	N 	 = size(u,1)
	rdof = size(stencils, 2)-4
	allocate(R_coeffs(N, rdof), V(N,N), Vinv(N,N))
	V 	 = 0
	Vinv = 0
	do j = 1, N
		R_coeffs(j,:) = recon_coeffs_1d(stencils(j,:), u, matrices(j,:,:),&
										partition, j)
		V(j,j) 		= partition%volumes(j)
		Vinv(j,j) 	= 1./partition%volumes(j)
	end do
	
	
	!choose numerical flux, i.e. approx of RP(a,b).
	!fluxtype == 1 -> RP(a,b) = a
	!fluxtype == 2 -> RP(a,b) = b
	!fluxtype == 3 -> RP(a,b) = (a+b)/2
	!fluxtype == 4 -> RP(a,b) = (a+b)/2 - lambda(b-a)
	
	fluxtype = 3
	allocate (res(N))
	res = 0
	if (fluxtype==1) then
		do i = 1, N
			do j = 1, 2
				if (j == 1) then
					res(i) = &
						upwindflux(partition%points(&
									partition%intervals(i,j)))
				else
					res(i) = res(i) - &
						upwindflux(partition%points(&
									partition%intervals(i,j)))
				end if
			end do
			if (partition%points(partition%intervals(i,2)) < &
				partition%points(partition%intervals(i,1))) then
				res(i) = -res(i)
			end if
		end do
	else if (fluxtype == 3) then
		do i = 1, N
			do j = 1, 2
				if (j == 1) then
					res(i) = &
						centralflux(partition%points(&
									partition%intervals(i,j)))
				else
					res(i) = res(i) - &
						centralflux(partition%points(&
									partition%intervals(i,j)))
				end if
			end do
			if (partition%points(partition%intervals(i,2)) < &
				partition%points(partition%intervals(i,1))) then
				res(i) = -res(i)
			end if
		end do
	else if (fluxtype == 4) then
		do i = 1, N
			if (partition%points(partition%intervals(i,2)) > &
				partition%points(partition%intervals(i,1))) then
				do j = 1, 2
					if (j == 1) then
						res(i) = &
							LFflux(partition%points(&
										partition%intervals(i,j)))
					else
						res(i) = res(i) - &
							LFflux(partition%points(&
										partition%intervals(i,j)))
					end if
				end do
			else
			 	do j = 1, 2
					if (j == 1) then
						res(i) = &
							-LFflux(partition%points(&
										partition%intervals(i,j)))
					else
						res(i) = res(i) + &
							LFflux(partition%points(&
										partition%intervals(i,j)))
					end if
				end do
			end if
		end do
	end if

	res = matmul(Vinv, res)
	
	contains
	
	pure function centralflux(x) result(fres)

		real(real64), intent(in)	:: x
		integer						:: sigma, m
		real(real64)				:: fres, Ri, Rj
		
		sigma  = size(R_coeffs,2)
		Ri = 0
		Rj = 0
		do m = 1, sigma
			Ri = Ri + R_coeffs(i,m)*((x-partition%midpoints(i))**degvec(m))
		end do
		if (partition%neighbours(i,j) == 0) then
			Rj = Ri
		else
			do m = 1, sigma
				Rj = Rj + &
					R_coeffs(partition%neighbours(i,j),m)*&
					((x-partition%midpoints(partition%neighbours(i,j)))**&
					degvec(m))
			end do
		end if
		fres = (Ri+Rj)/2.
		if (i == 1 .and. j == 1) fres = sin(-2*pi*t)

	end function centralflux
	
	pure function upwindflux(x) result(fres)

		real(real64), intent(in)	:: x
		integer						:: sigma, m
		real(real64)				:: fres, Ri, Rj
		
		sigma  = size(R_coeffs,2)
		Ri = 0
		Rj = 0
		do m = 1, sigma
			Ri = Ri + R_coeffs(i,m)*((x-partition%midpoints(i))**degvec(m))
		end do
		if (partition%neighbours(i,j) == 0) then
			Rj = Ri
		else
			do m = 1, sigma
				Rj = Rj + &
					R_coeffs(partition%neighbours(i,j),m)*&
					((x-partition%midpoints(partition%neighbours(i,j)))**&
					degvec(m))
			end do
		end if
		if (j == 1) fres = Rj
		if (j == 2) fres = Ri
		if (i == 1 .and. j == 1) fres = sin(-2*pi*t)

	end function upwindflux
	
	pure function LFflux(x) result(fres)

		real(real64), intent(in)	:: x
		integer						:: sigma, m
		real(real64)				:: fres, Ri, Rj
		
		sigma  = size(R_coeffs,2)
		Ri = 0
		Rj = 0
		do m = 1, sigma
			Ri = Ri + R_coeffs(i,m)*((x-partition%midpoints(i))**degvec(m))
		end do
		if (partition%neighbours(i,j) == 0) then
			Rj = Ri
		else
			do m = 1, sigma
				Rj = Rj + &
					R_coeffs(partition%neighbours(i,j),m)*&
					((x-partition%midpoints(partition%neighbours(i,j)))**&
					degvec(m))
			end do
		end if
		fres = (Ri+Rj)/2. - maxval([abs(Ri), abs(Rj)])*(Rj-Ri)/2.
		if (i == 1 .and. j == 1) fres = sin(-2*pi*t)

	end function LFflux

end function scheme_1d

!2D scheme

pure function scheme(u, t, grid, stencils, matrices) result(res)

	type(triangulation), intent(in)				:: grid
	real(real64), intent(in), dimension(:,:,:) 	:: matrices
	real(real64), intent(in), dimension(:)		:: u
	real(real64), intent(in)					:: t
	integer, intent(in), dimension(:,:) 		:: stencils
	integer, dimension(:,:), allocatable		:: degvec
	real(real64), dimension(:,:), allocatable	:: R_coeffs, Vinv
	real(real64), dimension(:), allocatable 	:: res
	real(real64)								:: pi
	integer										:: i, j, N, rdof
	
	allocate(degvec(15,2))
	pi 			= 4.D0*DATAN(1.D0)
	degvec(:,1) = [0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0]
	degvec(:,2) = [0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4]
	N 	 		= size(u,1)
	rdof 		= size(stencils, 2)-4
	allocate(R_coeffs(N, rdof), Vinv(N,N))
	do j = 1, N
			R_coeffs(j,:) = recon_coeffs(stencils(j,:), u, matrices(j,:,:),&
										grid, j)
			Vinv = 1./grid%volumes(j)
	end do
	
	allocate (res(N))
	res = 0
	do i = 1, N
		do j = 1, 3
			res(i) = res(i) -&
					 gauss_legendre_2d(&
					 grid%points(grid%connectivity_list(i,j),1), &
					 grid%points(grid%connectivity_list(i,j+1),1), &
					 grid%points(grid%connectivity_list(i,j),2), &
					 grid%points(grid%connectivity_list(i,j+1),2), &
					 2, 1, centralflux) + &
					 gauss_legendre_2d(&
					 grid%points(grid%connectivity_list(i,j),1), &
					 grid%points(grid%connectivity_list(i,j+1),1), &
					 grid%points(grid%connectivity_list(i,j),2), &
					 grid%points(grid%connectivity_list(i,j+1),2), &
					 1, 1, centralflux)
			if (j == 3) then
				res(i) = res(i) -&
					 gauss_legendre_2d(&
					 grid%points(grid%connectivity_list(i,j),1), &
					 grid%points(grid%connectivity_list(i,j-2),1), &
					 grid%points(grid%connectivity_list(i,j),2), &
					 grid%points(grid%connectivity_list(i,j-2),2), &
					 2, 1, centralflux) + &
					 gauss_legendre_2d(&
					 grid%points(grid%connectivity_list(i,j),1), &
					 grid%points(grid%connectivity_list(i,j-2),1), &
					 grid%points(grid%connectivity_list(i,j),2), &
					 grid%points(grid%connectivity_list(i,j-2),2), &
					 1, 1, centralflux)
			end if
		end do
	end do
	
	res = matmul(Vinv, res)
	
	contains
	
	pure function centralflux(x, y) result(res)
	
		real(real64), intent(in)	:: x, y
		integer						:: N, m
		real(real64)				:: res, Ri, Rj
		
		N 	= size(R_coeffs,2)
		Ri = 0
		Rj = 0
		do m = 1, N
			Ri = Ri + &
				R_coeffs(i,m)*((x-grid%barycenters(i,1))**degvec(m,1))&
				*((y-grid%barycenters(i,2))**degvec(m,2))
		end do
		if (grid%neighbours(i,j) == 0) then
				Rj  = Ri
		else
			do m = 1, N
				Rj = Rj + &
					R_coeffs(grid%neighbours(i,j),m)&
					*((x-grid%barycenters(grid%neighbours(i,j),1))&
					**degvec(m,1))&
					*((y-grid%barycenters(grid%neighbours(i,j),2))**degvec(m,2))
			end do
		end if
		res = (Ri+Rj)/2.
		if(grid%points(grid%connectivity_list(i,j),1) == 0) then
			res = 0
		end if
		if(grid%points(grid%connectivity_list(i,j),2) == 0) then
			res = 0
		end if
		
	end function centralflux

end function scheme

end module mod_scheme
