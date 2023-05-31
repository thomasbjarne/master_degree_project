module mod_moments

	use iso_fortran_env, only: real32, real64
	use mod_derived_types
	use mod_quadrature
	implicit none
	
contains

!1D procedures

pure function primarymoments1d(xdegree, partition) result(momentvector)
!works
	type(partition1d), intent(in)			:: partition
	integer, intent(in)						:: xdegree
	real(real64), dimension(:), allocatable	:: momentvector
	real(real64)							:: xl, xu
	integer									:: i, N, quadorder
	
	N = size(partition%intervals,1)
	allocate(momentvector(N))
	quadorder = 1
	if (xdegree == 2 .or. xdegree == 3) quadorder = 2
	if (xdegree == 4 .or. xdegree == 5) quadorder = 3
	if (xdegree == 6 .or. xdegree == 7) quadorder = 4
	if (xdegree > 7) quadorder = 5
	do i = 1, N
		xl = minval([partition%points(partition%intervals(i,1)),&
					 partition%points(partition%intervals(i,2))])
		xu = maxval([partition%points(partition%intervals(i,1)),&
					 partition%points(partition%intervals(i,2))])
		momentvector(i) = transformed_gauss_legendre(xl, xu, quadorder, func)/&
							partition%volumes(i)
	end do
	
	contains
	
	pure function func(x) result(res)
		
		real(real64), intent(in) :: x
		real(real64)			 :: res
		
		res = (x-partition%midpoints(i))**(xdegree)	
		
	end function func

end function primarymoments1d

pure function secondarymoments1d(xdegree, partition, indexnum) result(vec)
!works
	type(partition1d), intent(in)			:: partition
	integer, intent(in)						:: xdegree, indexnum
	real(real64), dimension(:), allocatable :: vec
	real(real64)							:: xl, xu
	integer									:: j, N, quadorder
	
	N = size(partition%intervals,1)
	allocate(vec(N))
	quadorder = 1
	if (xdegree == 2 .or. xdegree == 3) quadorder = 2
	if (xdegree == 4 .or. xdegree == 5) quadorder = 3
	if (xdegree == 6 .or. xdegree == 7) quadorder = 4
	if (xdegree > 7) quadorder = 5
	do j = 1, N
		xl = minval([partition%points(partition%intervals(j,1)),&
					 partition%points(partition%intervals(j,2))])
		xu = maxval([partition%points(partition%intervals(j,1)),&
					 partition%points(partition%intervals(j,2))])
		vec(j) = transformed_gauss_legendre(xl, xu, quadorder, func)/&
					partition%volumes(j)
	end do

	contains
	
	pure function func(x) result(res)
		
		real(real64), intent(in) :: x
		real(real64)			 :: res
		
		res = (x-partition%midpoints(indexnum))**xdegree	
		
	end function func

end function secondarymoments1d

!2D procedures
pure function calculate_moments(xdegree, ydegree, grid) result(momentvector)
!works
	type(triangulation), intent(in)					:: grid
	integer, intent(in)								:: xdegree, ydegree
	real(real64), dimension(:), allocatable			:: momentvector
	real(real64), dimension(2)						:: xb, yb
	integer											:: N, j, m, quadorder
	
	N = size(grid%connectivity_list,1)
	allocate (momentvector(N))
	momentvector = 0
	quadorder = 1
	if (xdegree == 2 .or. xdegree == 3 .or. ydegree == 2 .or. ydegree == 3) then
		quadorder = 3
	end if
	if (xdegree == 4 .or. xdegree == 5 .or. ydegree == 4 .or. ydegree == 5) then
		quadorder = 4
	end if
	do j = 1, N
		xb = 0
		yb = 0
		do m = 1, 2
			xb = [grid%points(grid%connectivity_list(j,m),1), &
				  grid%points(grid%connectivity_list(j,m+1),1)]
			yb = [grid%points(grid%connectivity_list(j,m),2), &
				  grid%points(grid%connectivity_list(j,m+1),2)]
			momentvector(j) = momentvector(j) +&
				gauss_legendre_2d(xb(1), xb(2), yb(1), yb(2), 2,quadorder,func)
		end do
		xb = [grid%points(grid%connectivity_list(j,3),1), &
			  grid%points(grid%connectivity_list(j,1),1)]
		yb = [grid%points(grid%connectivity_list(j,3),2), &
			  grid%points(grid%connectivity_list(j,1),2)]
		momentvector(j) = momentvector(j) +&
			gauss_legendre_2d(xb(1), xb(2), yb(1), yb(2), 2, quadorder, func)
		momentvector(j) = momentvector(j)/(grid%volumes(j)*(xdegree+1))
	end do
	
	momentvector = abs(momentvector)
	contains
	
	pure function func(x, y) result(res)
	
		real(real64), intent(in) 	:: x, y
		real(real64)				:: res
		
		res = ((x-grid%barycenters(j,1))**(xdegree+1)) * &
					((y-grid%barycenters(j,2))**ydegree)
	
	end function func
															
end function calculate_moments

pure function derived_moments(xdegree, ydegree, grid, indexnum, moments) &
														result(momentvector2)
!needs to be fixed.
	type(triangulation), intent(in)				:: grid
	real(real64), intent(in), dimension(:,:) 	:: moments
	integer, intent(in)							:: indexnum, xdegree, ydegree
	real(real64), dimension(:), allocatable		:: momentvector2, innerval
	integer										:: N, j, k, l
	
	N = size(grid%connectivity_list,1)
	allocate (momentvector2(N), innerval(N))
	!moments : 1, x, y, x^2, xy, y^2, x^3, 
	!			x^2y, xy^2, y^3, x^4, x^3y, x^2y^2, xy^3, y^4
	
	momentvector2	= 0
	innerval 		= 0
	do j = 1, N
		do l = 0, ydegree
			do k = 0, xdegree
				innerval(j) = innerval(j) + &
					binom(xdegree, k) *&
					((grid%barycenters(j,1)-grid%barycenters(indexnum,1))**k)*&
					primary(l,k)
			end do
			momentvector2(j) = momentvector2(j) + binom(ydegree, l)*&
					((grid%barycenters(j,2)-grid%barycenters(indexnum,2))**l)* &
					innerval(j)
			innerval = 0
		end do
	end do
	
	!momentvector2 = abs(momentvector2)
	contains
	
	pure function primary(l,k) result (res)
		
		integer,intent(in)	:: l, k
		integer				:: correctarg
		real(real64)		:: res
		
		if (ydegree-l == 0) then
			if (xdegree-k == 0) correctarg = 1
			if (xdegree-k == 1) correctarg = 2
			if (xdegree-k == 2) correctarg = 4
			if (xdegree-k == 3) correctarg = 7
			if (xdegree-k == 4) correctarg = 11
		else if (ydegree-l == 1) then
			if (xdegree-k == 0) correctarg = 3
			if (xdegree-k == 1) correctarg = 5
			if (xdegree-k == 2) correctarg = 8
			if (xdegree-k == 3) correctarg = 12
		else if (ydegree-l == 2) then
			if (xdegree-k == 0) correctarg = 6
			if (xdegree-k == 1) correctarg = 9
			if (xdegree-k == 2) correctarg = 13
		else if (ydegree-l == 3) then
			if (xdegree-k == 0) correctarg = 10
			if (xdegree-k == 1) correctarg = 14
		else if (ydegree-l == 4) then
			correctarg = 15
		end if 				
		
		res = moments(j,correctarg)
	
	end function primary
	
end function derived_moments

pure function factorial(n) result(nfac)

	integer, intent(in) :: n
	integer				:: nfac, i
	
	nfac = 1
	do i = 2, n
		nfac = i*nfac
	end do

end function factorial

pure function binom(n,k) result(binomcoeff)

	integer, intent(in)		:: n, k
	integer					:: binomcoeff
	
	binomcoeff = factorial(n)/(factorial(k)*factorial(n-k))

end function binom

end module mod_moments
