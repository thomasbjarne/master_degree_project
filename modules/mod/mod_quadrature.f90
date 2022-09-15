module mod_quadrature

	implicit none

contains

	pure function gauss_quadrature(fx) result(res)

		real, intent(in), dimension(:) :: fx
		real :: res
		integer :: i, n

		n = size(fx)
		

	end function gauss_quadrature

end module mod_quadrature