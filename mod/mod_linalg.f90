module mod_linalg

	use iso_fortran_env, only: real32, real64
	implicit none

    interface qr_decomposition
        module procedure :: qr_decomposition_kind32, qr_decomposition_kind64
    end interface qr_decomposition

    interface vec_outer_product
        module procedure :: vec_outer_product_kind32, vec_outer_product_kind64
    end interface vec_outer_product

    interface backwards_subst
        module procedure :: backwards_subst_kind32, backwards_subst_kind64
    end interface backwards_subst

    interface forwards_subst
        module procedure :: forwards_subst_kind32, forwards_subst_kind64
    end interface forwards_subst
    
    interface least_squares_qr
    	module procedure :: least_squares_qr_kind32, least_squares_qr_kind64
    end interface least_squares_qr
    
    interface constrained_least_squares
    	module procedure :: constrained_least_squares_kind32, &
    						constrained_least_squares_kind64
    end interface constrained_least_squares
    
contains

pure function identity_matrix(m,n) result(I)

	integer, intent(in)		:: m, n
	real, dimension(m,n)	:: I
	integer 				:: k

	I = 0
	if (m > n) then
		do k = 1, n
			I(k,k) = 1
		end do
	else
		do k = 1, m
			I(k,k) = 1
		end do
	end if
	
end function identity_matrix

pure function vec_outer_product_kind32(u, v) result(A)

	real(real32), intent(in), dimension(:) 		:: u, v
	real(real32), dimension(size(u),size(v)) 	:: A
	integer 									:: i, j

	do i = 1, size(u)
		do j = 1, size(v)
		    A(i, j) = u(i) * v(j)
		end do
	end do

end function vec_outer_product_kind32

pure function vec_outer_product_kind64(u, v) result(A)

	real(real64), intent(in), dimension(:) 		:: u, v
	real(real64), dimension(size(u),size(v)) 	:: A
	integer 									:: i, j

	do i = 1, size(u)
		do j = 1, size(v)
		    A(i, j) = u(i) * v(j)
		end do
	end do

end function vec_outer_product_kind64

pure subroutine qr_decomposition_kind32(A, Q, R)

	real(real32), intent(in), dimension(:,:) 					:: A
	real(real32), intent(out), dimension(size(A,1), size(A,1)) 	:: Q
	real(real32), intent(out), dimension(size(A,1), size(A,2)) 	:: R
	real(real32), dimension(size(A,1), size(A,1)) 				:: ID, Q_i
	real(real32), dimension(:), allocatable 					:: u
	real(real32) 												:: alpha
	integer 													:: i, n, m

	m	= size(A,1)
	n 	= size(A,2)
	ID 	= identity_matrix(m,m)
	Q 	= ID
	R 	= A
	do i = 1, n
		allocate (u(size(A,1) - i + 1))
		alpha 					= norm2(R(i:m,i))
		if (R(i,i) < 0) alpha 	= - alpha
		u 						= R(i:m, i) + alpha * ID(i:m, i) 
		if (norm2(u) /= 0) u 	= u / norm2(u)
		Q_i 					= ID
		Q_i(i:, i:) 			= Q_i(i:, i:) - 2 * vec_outer_product(u,u)
		Q 						= matmul(Q, transpose(Q_i))
		R(i:m, i:n) 			= R(i:m, i:n) &
							- 2 * vec_outer_product(u, matmul(u, R(i:m, i:n)))
		deallocate(u)
	end do

end subroutine qr_decomposition_kind32

pure subroutine qr_decomposition_kind64(A, Q, R)

	real(real64), intent(in), dimension(:,:) 					:: A
	real(real64), intent(out), dimension(size(A,1), size(A,1)) 	:: Q
	real(real64), intent(out), dimension(size(A,1), size(A,2)) 	:: R
	real(real64), dimension(size(A,1), size(A,1)) 				:: ID, Q_i
	real(real64), dimension(:), allocatable 					:: u
	real(real64) 												:: alpha
	integer 													:: i, n, m

	m	= size(A,1)
	n 	= size(A,2)
	ID 	= identity_matrix(m,m)
	Q 	= ID
	R 	= A
	do i = 1, n
		allocate(u(size(A,1) - i + 1))
		alpha 					= norm2(R(i:m,i))
		if (R(i,i) < 0) alpha 	= - alpha
		u 						= R(i:m, i) + alpha * ID(i:m, i) 
		if (norm2(u) /= 0) u 	= u / norm2(u)
		Q_i 					= ID
		Q_i(i:, i:) 			= Q_i(i:, i:) - 2 * vec_outer_product(u,u)
		Q 						= matmul(Q, transpose(Q_i))
		R(i:m, i:n) 			= R(i:m, i:n) &
							- 2 * vec_outer_product(u, matmul(u, R(i:m, i:n)))
		deallocate(u)
	end do

end subroutine qr_decomposition_kind64

pure function backwards_subst_kind32(R, b) result(x)

	real(real32), intent(in), dimension(:,:) 		:: R
	real(real32), intent(in), dimension(size(R,2)) 	:: b
	real(real32), dimension(size(R,2)) 				:: x
	integer 										:: i, n

	n = size(x)
	do i = n, 1, -1
		if (i == n) then
		    x(i) = ( b(i) ) / R(i,i)
		else  
		    x(i) = ( b(i)  - dot_product( x(i+1:n), R(i,i+1:n) ) ) / R(i,i)
		end if
	end do

end function backwards_subst_kind32

pure function backwards_subst_kind64(R, b) result(x)

	real(real64), intent(in), dimension(:,:) 		:: R
	real(real64), intent(in), dimension(size(R,2)) 	:: b
	real(real64), dimension(size(R,2)) 				:: x
	integer 										:: i, n

	n = size(x)
	do i = n, 1, -1
		if (i == n) then
		    x(i) = ( b(i) ) / R(i,i)
		else  
		    x(i) = ( b(i)  - dot_product( x(i+1:n), R(i,i+1:n) ) ) / R(i,i)
		end if
	end do

end function backwards_subst_kind64

pure function forwards_subst_kind32(L, b) result(x)

	real(real32), intent(in), dimension(:,:) 		:: L
	real(real32), intent(in), dimension(size(L,2)) 	:: b
	real(real32), dimension(size(L,2)) 				:: x
	integer 										:: i, n

	n = size(L,1)
	x = 0
	do i = 1, n
		if (i == 1) then
		    x(i) = ( b(i) ) / L(i,i)
		else
		    x(i) = ( b(i) - dot_product( x(1:i-1), L(i, 1:i-1) ) ) / L(i,i)
		end if
	end do

end function forwards_subst_kind32

pure function forwards_subst_kind64(L, b) result(x)

	real(real64), intent(in), dimension(:,:) 		:: L
	real(real64), intent(in), dimension(size(L,2)) 	:: b
	real(real64), dimension(size(L,2)) 				:: x
	integer 										:: i, n

	n = size(x)
	do i = 1, n
		if (i == 1) then
		    x(i) = ( b(i) ) / L(i,i)
		else
		    x(i) = ( b(i) - dot_product( x(1:i-1), L(i, 1:i-1) ) ) / L(i,i)
		end if
	end do

end function forwards_subst_kind64

pure function least_squares_qr_kind32(A, b) result(x)

	real(real32), intent(in), dimension(:,:) 	:: A
	real(real32), intent(in), dimension(:)		:: b
	real(real32), dimension(:), allocatable		:: x
	real(real32), dimension(:,:), allocatable 	:: Q, R
	integer										:: m, n
	
	m = size(A,1)
	n = size(A,2)
	allocate (x(n), Q(m,m), R(m,n))
	call qr_decomposition(A, Q, R)
	x = backwards_subst(R(1:n,1:n), matmul(transpose(Q),b))

end function least_squares_qr_kind32

pure function least_squares_qr_kind64(A, b) result(x)

	real(real64), intent(in), dimension(:,:) 	:: A
	real(real64), intent(in), dimension(:)		:: b
	real(real64), dimension(:), allocatable		:: x
	real(real64), dimension(:,:), allocatable 	:: Q, R
	integer										:: m, n
	
	m = size(A,1)
	n = size(A,2)
	allocate (x(n), Q(m,m), R(m,n))
	call qr_decomposition(A, Q, R)
	x = backwards_subst(R(1:n,1:n), matmul(transpose(Q),b))

end function least_squares_qr_kind64

pure function constrained_least_squares_kind32(A, b) result(x)
	
	!constraint is A(1,:)*x = b(1)
	real(real32), intent(in), dimension(:,:) 	:: A
	real(real32), intent(in), dimension(:)		:: b
	real(real32), dimension(:), allocatable		:: x, v
	real(real32), dimension(:,:), allocatable	:: U, Aug
	integer										:: m, n, j
	
	m = size(A,1)
	n = size(A,2)
	allocate (x(n), U(m,n), v(m), Aug(m,n+1))
	U = A
	v = b
	Aug(:,1:n) = U
	Aug(:,n+1) = v
	do j = 2, m
		Aug(j,:) = Aug(j,:) - (Aug(j,1)/Aug(1,1))*Aug(1,:)
	end do
	x 	= least_squares_qr(Aug(:,1:n),Aug(:,n+1))	
	
end function constrained_least_squares_kind32

pure function constrained_least_squares_kind64(A, b) result(x)
	
	!constraint is A(1,:)*x = b(1)
	real(real64), intent(in), dimension(:,:) 	:: A
	real(real64), intent(in), dimension(:)		:: b
	real(real64), dimension(:), allocatable		:: x, v
	real(real64), dimension(:,:), allocatable	:: U, Aug
	integer										:: m, n, j
	
	m = size(A,1)
	n = size(A,2)
	allocate (x(n), U(m,n), v(m), Aug(m,n+1))
	U = A
	v = b
	Aug(:,1:n) = U
	Aug(:,n+1) = v
	do j = 2, m
		Aug(j,:) = Aug(j,:) - (Aug(j,1)/Aug(1,1))*Aug(1,:)
	end do
	x 	= least_squares_qr(Aug(:,1:n),Aug(:,n+1))	

end function constrained_least_squares_kind64


end module mod_linalg
