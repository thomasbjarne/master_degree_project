module mod_quadrature

	use iso_fortran_env, only: real32, real64
	implicit none
		
	interface gauss_legendre
		module procedure :: gauss_legendre_kind32, gauss_legendre_kind64
	end interface gauss_legendre
	
	interface transformed_gauss_legendre
		module procedure :: transformed_gauss_legendre_kind32,&
							transformed_gauss_legendre_kind64
	end interface transformed_gauss_legendre
	
	interface gauss_legendre_2d
		module procedure :: gauss_legendre_2d_kind32,&
							gauss_legendre_2d_kind64
	end interface gauss_legendre_2d
	
contains

pure function gauss_legendre_kind32(degree, func) result(value)

	integer, intent(in)			:: degree
	real(real32)				:: value
	real(real32), dimension(2) 	:: deg_2_GL_quadcoeff, deg_2_legendre_roots
	real(real32), dimension(3) 	:: deg_3_GL_quadcoeff, deg_3_legendre_roots
	real(real32), dimension(4) 	:: deg_4_GL_quadcoeff, deg_4_legendre_roots
	real(real32), dimension(5) 	:: deg_5_GL_quadcoeff, deg_5_legendre_roots
	integer						:: i, n
	interface
		pure function func(z)
			import 					 :: real32
			real(real32) 			 :: func
			real(real32), intent(in) :: z
		end function func
	end interface
	
	deg_2_legendre_roots = [-sqrt(1./3.), sqrt(1./3.)]
	deg_3_legendre_roots = [-sqrt(3./5.),0.,sqrt(3./5.)]
	deg_4_legendre_roots = [-sqrt((15. + 2.*sqrt(30.))/35.),&
							-sqrt((15. - 2.*sqrt(30.))/35.),&
							sqrt((15. - 2.*sqrt(30.))/35.) ,&
							sqrt((15. + 2.*sqrt(30.))/35.)]
	deg_5_legendre_roots = [-(1./3.)*sqrt(5.+2.*sqrt(10./7.)),&
							-(1./3.)*sqrt(5.-2.*sqrt(10./7.)),&
							0.,&
							(1./3.)*sqrt(5.-2.*sqrt(10./7.)),&
							(1./3.)*sqrt(5.+2.*sqrt(10./7.))]
	deg_2_GL_quadcoeff = [1., 1.]
	deg_3_GL_quadcoeff = [5./9., 8./9., 5./9.]
	deg_4_GL_quadcoeff = [(90.-5.*sqrt(30.))/180., &
						  (90.+5.*sqrt(30.))/180., &
						  (90.+5.*sqrt(30.))/180., &
						  (90.-5.*sqrt(30.))/180.]
	deg_5_GL_quadcoeff = [(322.-13.*sqrt(70.))/900.,&
						(322.+13.*sqrt(70.))/900.,&
						128./225.,&
						(322.+13.*sqrt(70.))/900.,&
						(322.-13.*sqrt(70.))/900.]
	
	value = 0
	if (degree == 1) then
		value = 2.*func(0.)
	else if (degree == 2) then
		n = 2
		do i = 1, n
			value = value + deg_2_GL_quadcoeff(i)*func(deg_2_legendre_roots(i))
		end do
	else if (degree == 3) then
		n = 3
		do i = 1, n
			value = value + deg_3_GL_quadcoeff(i)*func(deg_3_legendre_roots(i))
		end do
	else if (degree == 4) then
		n = 4
		do i = 1, n
			value = value + deg_4_GL_quadcoeff(i)*func(deg_4_legendre_roots(i))
		end do
	else if (degree == 5) then
		n = 5
		do i = 1, n
			value = value + deg_5_GL_quadcoeff(i)*func(deg_5_legendre_roots(i))
		end do
	end if

end function gauss_legendre_kind32

pure function gauss_legendre_kind64(degree, func) result(value)

	integer, intent(in)			:: degree
	real(real64)				:: value
	real(real64), dimension(2) 	:: deg_2_GL_quadcoeff, deg_2_legendre_roots
	real(real64), dimension(3) 	:: deg_3_GL_quadcoeff, deg_3_legendre_roots
	real(real64), dimension(4) 	:: deg_4_GL_quadcoeff, deg_4_legendre_roots
	real(real64), dimension(5) 	:: deg_5_GL_quadcoeff, deg_5_legendre_roots
	integer						:: i, n
	interface
		pure function func(z)
			import 					 :: real64
			real(real64) 			 :: func
			real(real64), intent(in) :: z
		end function func
	end interface
	
	deg_2_legendre_roots = [-sqrt(1./3.), sqrt(1./3.)]
	deg_3_legendre_roots = [-sqrt(3./5.),0.,sqrt(3./5.)]
	deg_4_legendre_roots = [-sqrt((15. + 2.*sqrt(30.))/35.),&
							-sqrt((15. - 2.*sqrt(30.))/35.),&
							sqrt((15. - 2.*sqrt(30.))/35.) ,&
							sqrt((15. + 2.*sqrt(30.))/35.)]
	deg_5_legendre_roots = [-(1./3.)*sqrt(5.+ 2.*sqrt(10./7.)),&
							-(1./3.)*sqrt(5.- 2.*sqrt(10./7.)),&
							0.,&
							(1./3.)*sqrt(5.-2.*sqrt(10./7.)),&
							(1./3.)*sqrt(5.+2.*sqrt(10./7.))]
	deg_2_GL_quadcoeff = [1., 1.]
	deg_3_GL_quadcoeff = [5./9., 8./9., 5./9.]
	deg_4_GL_quadcoeff = [(90.-5.*sqrt(30.))/180., &
						  (90.+5.*sqrt(30.))/180., &
						  (90.+5.*sqrt(30.))/180., &
						  (90.-5.*sqrt(30.))/180.]
	deg_5_GL_quadcoeff = [(322.-13.*sqrt(70.))/900.,&
						(322.+13.*sqrt(70.))/900.,&
						128./225.,&
						(322.+13.*sqrt(70.))/900.,&
						(322.-13.*sqrt(70.))/900.]
	
	value = 0
	if (degree == 1) then
		value = 2.*func(0.D0)
	else if (degree == 2) then
		n = 2
		do i = 1, n
			value = value + deg_2_GL_quadcoeff(i)*func(deg_2_legendre_roots(i))
		end do
	else if (degree == 3) then
		n = 3
		do i = 1, n
			value = value + deg_3_GL_quadcoeff(i)*func(deg_3_legendre_roots(i))
		end do
	else if (degree == 4) then
		n = 4
		do i = 1, n
			value = value + deg_4_GL_quadcoeff(i)*func(deg_4_legendre_roots(i))
		end do
	else if (degree == 5) then
		n = 5
		do i = 1, n
			value = value + deg_5_GL_quadcoeff(i)*func(deg_5_legendre_roots(i))
		end do
	end if

end function gauss_legendre_kind64

pure function transformed_gauss_legendre_kind32(a,b, degree, func) result(res)

	real(real32), intent(in) 	:: a, b		!a = lower integral limit, b = upper
	integer, intent(in) 		:: degree
	real(real32)				:: res
	interface
		pure function func(z)
			import						:: real32
			real(real32)				:: func
			real(real32), intent (in) 	:: z
		end function func
	end interface
	
	res = (b-a)/2. * gauss_legendre(degree, transformedfunc)
	
	contains
	
	pure function transformedfunc(xi) result(tres)
		
		real(real32), intent(in) :: xi
		real(real32)			 :: tres
		
		tres = func((b-a)*xi/2. + (b+a)/2.)
		
	end function transformedfunc

end function transformed_gauss_legendre_kind32

pure function transformed_gauss_legendre_kind64(a,b, degree, func) result(res)

	real(real64), intent(in) 	:: a, b		!a = lower integral limit, b = upper
	integer, intent(in) 		:: degree
	real(real64)				:: res
	interface
		pure function func(z)
			import						:: real64
			real(real64)				:: func
			real(real64), intent (in) 	:: z
		end function func
	end interface
	
	res = (b-a)/2. * gauss_legendre(degree, transformedfunc)
	
	contains
	
	pure function transformedfunc(xi) result(tres)
		
		real(real64), intent(in) :: xi
		real(real64)			 :: tres
		
		tres = func((b-a)*xi/2. + (b+a)/2.)
		
	end function transformedfunc

end function transformed_gauss_legendre_kind64

pure function gauss_legendre_2d_kind32(x0, x1, y0, y1, intarg, degree, func)&
																	result(sol)

	real(real32), intent(in) 	:: x0, x1, y0, y1
	integer, intent(in)		 	:: degree, intarg
	real(real32)				:: sol
	
	interface
		pure function func(x, y)
			import						:: real32
			real(real32) 				:: func
			real(real32), intent(in) 	:: x, y
		end function func
	end interface
	
	sol = 0
	if (intarg == 1) then
		sol = (x1-x0)/2. * gauss_legendre(degree, transformedfunc2)
	else if (intarg == 2) then
		sol = (y1-y0)/2. * gauss_legendre(degree, transformedfunc2)
	end if
	
	contains
	
	pure function transformedfunc2(xi) result (tres)
		
		real(real32), intent(in)	:: xi
		real(real32)				:: tres
		
		tres = func((x1-x0)*xi/2. + (x1+x0)/2., &
					(y1-y0)*xi/2. + (y1+y0)/2.)
	
	end function transformedfunc2

end function gauss_legendre_2d_kind32

pure function gauss_legendre_2d_kind64(x0, x1, y0, y1, intarg, degree, func)&
																	result(sol)

	real(real64), intent(in) 	:: x0, x1, y0, y1
	integer, intent(in)		 	:: degree, intarg
	real(real64)				:: sol
	
	interface
		pure function func(x, y)
			import						:: real64
			real(real64) 				:: func
			real(real64), intent(in) 	:: x, y
		end function func
	end interface
	
	sol = 0
	if (intarg == 1) then
		sol = (x1-x0)/2. * gauss_legendre(degree, transformedfunc2)
	else if (intarg == 2) then
		sol = (y1-y0)/2. * gauss_legendre(degree, transformedfunc2)
	end if
	
	contains
	
	pure function transformedfunc2(xi) result (tres)
		
		real(real64), intent(in)	:: xi
		real(real64)				:: tres
		
		tres = func((x1-x0)*xi/2. + (x1+x0)/2., &
					(y1-y0)*xi/2. + (y1+y0)/2.)
	
	end function transformedfunc2

end function gauss_legendre_2d_kind64

end module mod_quadrature
