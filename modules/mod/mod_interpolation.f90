module mod_interpolation

	use iso_fortran_env, only: real32, real64
	implicit none
	
	interface lagrange_interpolate
		module procedure :: lagrange_interpolate_kind32, &
			lagrange_interpolate_kind64
	end interface lagrange_interpolate
	
	interface newton_div_diff_interpolate
		module procedure :: newton_div_diff_interpolate_kind32, &
			newton_div_diff_interpolate_kind64
	end interface newton_div_diff_interpolate

contains

	pure function lagrange_interpolate_kind32(pointset, x) result(Lx)

		real(real32), intent(in) :: x
		real(real32), intent(in), dimension(:,:) :: pointset
		real(real32), dimension(size(pointset,1)) :: basis_val
		real(real32) :: Lx
		integer :: i, j, n
		n = size(pointset, 1)
		Lx = 0
		basis_val = 1
		do j = 1, n
			do i = 1, n
				if (i == j) cycle
				basis_val(j) = basis_val(j) * ( (x-pointset(i,1))/(pointset(j,1)-pointset(i,1)) )
			end do
		end do

		Lx = sum(pointset(1:n,2)*basis_val(1:n))

	end function lagrange_interpolate_kind32
	
	pure function lagrange_interpolate_kind64(pointset, x) result(Lx)

		real(real64), intent(in) :: x
		real(real64), intent(in), dimension(:,:) :: pointset
		real(real64), dimension(size(pointset,1)) :: basis_val
		real(real64) :: Lx
		integer :: i, j, n
		n = size(pointset, 1)
		Lx = 0
		basis_val = 1
		do j = 1, n
			do i = 1, n
				if (i == j) cycle
				basis_val(j) = basis_val(j) * ( (x-pointset(i,1))/(pointset(j,1)-pointset(i,1)) )
			end do
		end do

		Lx = sum(pointset(1:n,2)*basis_val(1:n))

	end function lagrange_interpolate_kind64

	pure function newton_div_diff_interpolate_kind32(pointset, x) result(Px)

		real(real32), intent(in) :: x
		real(real32), intent(in), dimension(:,:) :: pointset
		real(real32) :: Px
		real(real32), dimension(size(pointset,1),size(pointset,1)) :: a
		real(real32), dimension(size(pointset,1)) :: cofactors
		integer :: i, j, n

		a = 0
		n = size(pointset,1)
		do j = 1, n
			a(j,1) = pointset(j,2)
		end do

		do i = 2, n
			do j = 1, n+1-i
				a(j,i) = ( a(j+1,i-1) - a(j,i-1) )/( pointset(j+i-1,1)-pointset(j,1) )
			end do
		end do

		Px = 0
		cofactors(1) = 1
		do i = 2, n
			cofactors(i) = cofactors(i-1)*(x-pointset(i-1,1))
		end do
		do i = 1, n
			Px = Px + (a(1,i)*cofactors(i))
		end do

	end function newton_div_diff_interpolate_kind32
	
	pure function newton_div_diff_interpolate_kind64(pointset, x) result(Px)

		real(real64), intent(in) :: x
		real(real64), intent(in), dimension(:,:) :: pointset
		real(real64) :: Px
		real(real64), dimension(size(pointset,1),size(pointset,1)) :: a
		real(real64), dimension(size(pointset,1)) :: cofactors
		integer :: i, j, n

		a = 0
		n = size(pointset,1)
		do j = 1, n
			a(j,1) = pointset(j,2)
		end do

		do i = 2, n
			do j = 1, n+1-i
				a(j,i) = ( a(j+1,i-1) - a(j,i-1) )/( pointset(j+i-1,1)-pointset(j,1) )
			end do
		end do

		Px = 0
		cofactors(1) = 1
		do i = 2, n
			cofactors(i) = cofactors(i-1)*(x-pointset(i-1,1))
		end do
		do i = 1, n
			Px = Px + (a(1,i)*cofactors(i))
		end do

	end function newton_div_diff_interpolate_kind64

	pure function chebyshev_roots(a, b, n) result(roots)

		real, intent(in) :: a, b
		integer, intent(in) :: n
		real, dimension(n) :: roots
		real, parameter :: pi = 3.14159265358979323846
		integer :: i

		do i = 1, n
			roots(i) = (b-a)/2 + ((b-a)/2)*cos((2*i-1)*pi)/2*n
		end do

	end function chebyshev_roots

	pure function natural_cubic_spline(pointset, x) result(Sx)

		real, intent(in), dimension(:,:) :: pointset
		real, intent(in) :: x
		real, dimension(size(pointset,1)) :: a, c 
		real, dimension(size(a)-1) :: h1, h2, d, b
		real, dimension(size(c), size(c)) :: delta_matrix
		real :: Sx
		integer :: i, j, n

		n = size(pointset,1)
		do i = 1, n
			a(i) = pointset(i,2)
		end do

		do i = 1, n-1
			h1(i) = pointset(i+1, 1) - pointset(i, 1)
			h2(i) = pointset(i+1, 2) - pointset(i, 2)
		end do

		




	end function natural_cubic_spline

end module mod_interpolation
