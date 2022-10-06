module mod_derived_types

	implicit none

	type :: triangle
		integer :: index = 1
		real, dimension(3, 2) :: vertices
		real, dimension(1, 2) :: center
		real :: area
	end type triangle

	type :: circle
		real :: radius, area, circumference
		real, dimension(1,2) :: center = 0
	end type circle

	type :: quadrilateral
		integer, dimension(1,2) :: index = 1
		real, dimension(4, 2) :: vertices
		real, dimension(1, 2) :: center
		real :: area
	end type quadrilateral
	
	interface triangle
		module procedure :: triangle_constructor
	end interface triangle

	interface circle
		module procedure :: circle_constructor
	end interface circle
	
	interface quadrilateral
		module procedure :: quadrilateral_constructor
	end interface quadrilateral

contains

	pure type(triangle) function triangle_constructor(index, vertices) &
	result(res)
		
		integer, intent(in) :: index
		real, intent(in), dimension(3,2) :: vertices
		real, dimension(1, 2) :: center
		real :: area
		res % index = index
		res % vertices = vertices
		res % center(1,1) = sum(vertices(1:3, 1))/3
		res % center(1,2) = sum(vertices(1:3, 2))/3

		area = 0.5* abs( (vertices(1,1) - vertices(3,1)) &
		* (vertices(2,2)-vertices(1,2))& 
		- (vertices(1,1) - vertices(2,1)*(vertices(3,2) - vertices(1,2)) ) &
		)

		res % area = area 

	end function triangle_constructor

	pure type(circle) function circle_constructor(radius, center) result(res)

		real, intent(in) :: radius
		real, intent(in), dimension(1,2) :: center
		real, parameter :: pi = 3.14159265358979323846
		res % radius = radius
		res % center = center
		res % circumference = 2 * pi * radius
		res % area = pi * radius ** 2

	end function circle_constructor

	pure type(quadrilateral) function quadrilateral_constructor(index, vertices)&
	result(res)

		integer, intent(in), dimension(1,2) :: index
		real, intent(in), dimension(4,2) :: vertices
		real, dimension(1, 2) :: center
		real :: area
		res % index = index
		res % vertices = vertices
		res % center(1, 1) = abs(maxval(vertices(:,1)) + minval(vertices(:,1)))/2
		res % center(1, 2) = abs(maxval(vertices(:,2)) + minval(vertices(:,2)))/2

		area = 0.5 * abs ( &
		(vertices(1,2) + vertices(2,2))*(vertices(1,1) &
		- vertices(2,1)) + (vertices(2,2)+vertices(3,2)) *&
		(vertices(2,1) - vertices(3,1)) + (vertices(3,2) + vertices(4,2))&
		* (vertices(3,1)-vertices(4,1)) + (vertices(4,2) + vertices(1,2))&
		* (vertices(4,1)-vertices(1,1)) &
		)

		res % area = area

	end function quadrilateral_constructor

end module mod_derived_types
