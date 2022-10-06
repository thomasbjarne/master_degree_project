module mod_quadrature

	use mod_interpolation
	implicit none

contains

	pure function integrate_polynomial(low, high, P_coefficients) result(res)

		real, intent(in) :: low, high
		real, intent(in), dimension(:) :: P_coefficients
		real, dimension(size(P_coefficients)) :: res, F_high, F_low
		integer :: i, n

		n = size(P_coefficients)
		res = 0
		do i = 1, n
			F_high = F_high + P_coefficients(i)*(high**i)/(i+1)
			F_low = F_low + P_coefficients(i)*(low**i)/(i+1)
		end do

		res = F_high - F_low

	end function integrate_polynomial

	pure function gauss_quadrature(a, b, fx) result(res)

		real, intent(in) :: a, b
		real, intent(in), dimension(:) :: fx
		real, dimension(size(fx)) :: legendre_coefficients
		real :: res
		integer :: i

		!legendre_coefficients(i) = integrate_polynomial(lagrange_interpolate(pointset, x))

		res = 0
		do i = 1, size(fx)
			res = res + legendre_coefficients(i)*fx(i)
		end do
		

	end function gauss_quadrature

end module mod_quadrature