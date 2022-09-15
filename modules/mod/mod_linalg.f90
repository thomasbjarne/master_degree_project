module mod_linalg

    use iso_fortran_env, only: real32, real64
    implicit none

    interface qr_decomposition
        module procedure :: qr_decomposition_kind32, qr_decomposition_kind64
    end interface qr_decomposition

    interface vec_outer_product
        module procedure :: vec_outer_product_kind32, vec_outer_product_kind64
    end interface vec_outer_product

    interface qr_algorithm
        module procedure :: qr_algorithm_kind32, qr_algorithm_kind64
    end interface qr_algorithm

    interface upper_hessenberg
        module procedure :: upper_hessenberg_kind32, upper_hessenberg_kind64
    end interface upper_hessenberg

    interface backwards_subst
        module procedure :: backwards_subst_kind32, backwards_subst_kind64
    end interface backwards_subst

    interface forwards_subst
        module procedure :: forwards_subst_kind32, forwards_subst_kind64
    end interface forwards_subst

    interface gaussian_elim
        module procedure :: gaussian_elim_kind32, gaussian_elim_kind64
    end interface gaussian_elim

    interface shifted_iteration
        module procedure :: shifted_iteration_kind32, shifted_iteration_kind64
    end interface shifted_iteration

    interface arnoldi_iteration
        module procedure :: arnoldi_iteration_kind32, arnoldi_iteration_kind64
    end interface arnoldi_iteration

contains
! TO FIX : shifted_iteration, gaussian_elim(?)

    pure function identity_matrix(n) result(I)

        integer, intent(in) :: n
        real, dimension(n,n) :: I
        integer :: k

        I = 0
        do k = 1, n
            I(k,k) = 1
        end do

    end function identity_matrix

    pure subroutine arnoldi_iteration_kind32(A, n, Q, H)

        real(real32), intent(in), dimension(:,:) :: A
        integer, intent(in) :: n
        real(real32), intent(out), dimension(size(A,1), n+1) :: Q
        real(real32), intent(out), dimension(n+1, n) :: H
        real(real32), dimension(size(A,1)) :: b
        real(real32) :: eps
        integer :: i, j

        b = 1
        Q = 0
        H = 0
        eps = 1e-10
        Q(:, 1) = b/norm2(b)
        do j = 1, n
            Q(:,j+1) = matmul(A, Q(:, j))
            do i = 1, j
                H(i, j) = dot_product(Q(:,j+1),Q(:,i))
                Q(:,j+1) = Q(:,j+1) - H(i, j) * Q(:,i)
            end do
            H(j+1, j) = norm2(Q(:,j+1))
            Q(:,j+1) = Q(:,j+1)/H(j+1,j)
            if (H(j+1,j) < eps) exit
        end do

    end subroutine arnoldi_iteration_kind32

    pure subroutine arnoldi_iteration_kind64(A, n, Q, H)

        real(real64), intent(in), dimension(:,:) :: A
        integer, intent(in) :: n
        real(real64), intent(out), dimension(size(A,1), n+1) :: Q
        real(real64), intent(out), dimension(n+1, n) :: H
        real(real64), dimension(size(A,1)) :: b
        real(real64) :: eps
        integer :: i, j

        b = 1
        Q = 0
        H = 0
        eps = 1e-10
        Q(:, 1) = b/norm2(b)
        do j = 1, n
            Q(:,j+1) = matmul(A, Q(:, j))
            do i = 1, j
                H(i, j) = dot_product(Q(:,j+1),Q(:,i))
                Q(:,j+1) = Q(:,j+1) - H(i, j) * Q(:,i)
            end do
            H(j+1, j) = norm2(Q(:,j+1))
            Q(:,j+1) = Q(:,j+1)/H(j+1,j)
            if (H(j+1,j) < eps) exit
        end do

    end subroutine arnoldi_iteration_kind64

    pure subroutine qr_decomposition_kind32(A, Q, R)

        real(real32), intent(in), dimension(:,:) :: A
        real(real32), intent(out), dimension(size(A,1), size(A,2)) :: Q
        real(real32), intent(out), dimension(size(A,2), size(A,2)) :: R
        real(real32), dimension(size(A,1), size(A,2)) :: ID, Q_i
        real(real32), dimension(:), allocatable :: u
        real(real32) :: alpha
        integer :: i, n

        n = size(A,2)

        ID = identity_matrix(n)
        Q = ID
        R = A(1:n,1:n)
        do i = 1, n
            allocate( u(size(A,1) - i + 1) )
            alpha = -norm2( R(i:n,i) )
            if (R(i,i) < 0) alpha = - alpha
            u = R(i:n, i) + alpha * ID(i:n, i) 
            if (norm2(u) /= 0) u = u / norm2( u )
            Q_i = ID
            Q_i(i:, i:) = Q_i(i:, i:) - 2 * vec_outer_product(u,u)
            Q = matmul(Q, transpose(Q_i))
            R(i:n, i:n) = R(i:n, i:n) - 2 * vec_outer_product(u, matmul(u, R(i:n, i:n)))
            deallocate(u)
        end do

    end subroutine qr_decomposition_kind32

    pure subroutine qr_decomposition_kind64(A, Q, R)

        real(real64), intent(in), dimension(:,:) :: A
        real(real64), intent(out), dimension(size(A,1), size(A,2)) :: Q
        real(real64), intent(out), dimension(size(A,2), size(A,2)) :: R
        real(real64), dimension(size(A,1), size(A,2)) :: ID, Q_i
        real(real64), dimension(:), allocatable :: u
        real(real64) :: alpha
        integer :: i, n

        n = size(A,2)
        ID = identity_matrix(n)
        Q = ID
        R = A(1:n,1:n)
        do i = 1, n
            allocate( u(size(A,1) - i + 1) )
            alpha = -norm2( R(i:n,i) )
            if (R(i,i) < 0) alpha = - alpha
            u = R(i:n, i) + alpha * ID(i:n, i) 
            if (norm2(u) /= 0) u = u / norm2( u )
            Q_i = ID
            Q_i(i:, i:) = Q_i(i:, i:) - 2 * vec_outer_product(u,u)
            Q = matmul(Q, transpose(Q_i))
            R(i:n, i:n) = R(i:n, i:n) - 2 * vec_outer_product(u, matmul(u, R(i:n, i:n)))
            deallocate(u)
        end do

    end subroutine qr_decomposition_kind64

    pure function qr_algorithm_kind32(A) result(A_schur)

		real(real32), intent(in), dimension(:,:) :: A
		real(real32), dimension(:,:), allocatable :: A_schur, A_k, Q_k, R_k, H
		real(real32), parameter :: tol = 1e-14
		real(real32) :: mu
		integer :: m, n, count

		count = 0
		m = size(A,1)
		allocate(A_schur(m,m), H(m,m))
		n = m
		H = upper_hessenberg(A)
		do while (n > 1)
			allocate (A_k(1:n,1:n), Q_k(1:n,1:n), R_k(1:n,1:n))
			A_k(1:n,1:n) = H(1:n,1:n)
			do while (maxval(abs(A_k(n, 1:n-1))) > tol .and. count < 1500)
				count = count + 1
				mu = A_k(n,n)
				A_k = A_k - mu*identity_matrix(n)
				call qr_decomposition(A_k, Q_k, R_k)
				A_k = matmul(R_k, Q_k) + mu*identity_matrix(n)
			end do
			A_schur(1:n,1:n) = A_k
			n = n - 1
			deallocate(A_k, Q_k, R_k)
		end do

    end function qr_algorithm_kind32

    pure function qr_algorithm_kind64(A) result(A_schur)

		real(real64), intent(in), dimension(:,:) :: A
		real(real64), dimension(:,:), allocatable :: A_schur, A_k, Q_k, R_k, H
		real(real64), parameter :: tol = 1e-14
		real(real64) :: mu
		integer :: m, n, count

		count = 0
		m = size(A,1)
		allocate(A_schur(m,m), H(m,m))
		n = m
		H = upper_hessenberg(A)
		do while (n > 1)
			allocate (A_k(1:n,1:n), Q_k(1:n,1:n), R_k(1:n,1:n))
			A_k(1:n,1:n) = H(1:n,1:n)
			do while (maxval(abs(A_k(n, 1:n-1))) > tol .and. count < 1500)
				count = count + 1
				mu = A_k(n,n)
				A_k = A_k - mu*identity_matrix(n)
				call qr_decomposition(A_k, Q_k, R_k)
				A_k = matmul(R_k, Q_k) + mu*identity_matrix(n)
			end do
			A_schur(1:n,1:n) = A_k
			n = n - 1
			deallocate(A_k, Q_k, R_k)
		end do
		
    end function qr_algorithm_kind64

    pure function vec_outer_product_kind32(u, v) result(A)

        real(real32), intent(in), dimension(:) :: u, v
        real(real32), dimension(size(u),size(v)) :: A
        integer :: i, j

        do i = 1, size(u)
            do j = 1, size(v)
                A(i, j) = u(i) * v(j)
            end do
        end do

    end function vec_outer_product_kind32

    pure function vec_outer_product_kind64(u, v) result(A)

        real(real64), intent(in), dimension(:) :: u, v
        real(real64), dimension(size(u),size(v)) :: A
        integer :: i, j

        do i = 1, size(u)
            do j = 1, size(v)
                A(i, j) = u(i) * v(j)
            end do
        end do

    end function vec_outer_product_kind64

    pure function backwards_subst_kind32(R, b) result(x)

        real(real32), intent(in), dimension(:,:) :: R
        real(real32), intent(in), dimension(size(R,2)) :: b
        real(real32), dimension(size(R,2)) :: x
        integer :: i, n

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

        real(real64), intent(in), dimension(:,:) :: R
        real(real64), intent(in), dimension(size(R,2)) :: b
        real(real64), dimension(size(R,2)) :: x
        integer :: i, n

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

        real(real32), intent(in), dimension(:,:) :: L
        real(real32), intent(in), dimension(size(L,2)) :: b
        real(real32), dimension(size(L,2)) :: x
        integer :: i, n

        n = size(x)
        do i = 1, n
            if (i == 1) then
                x(i) = ( b(i) ) / L(i,i)
            else
                x(i) = ( b(i) - dot_product( x(1:i-1), L(i, 1:i-1) ) ) / L(i,i)
            end if
        end do

    end function forwards_subst_kind32

    pure function forwards_subst_kind64(L, b) result(x)

        real(real64), intent(in), dimension(:,:) :: L
        real(real64), intent(in), dimension(size(L,2)) :: b
        real(real64), dimension(size(L,2)) :: x
        integer :: i, n

        n = size(x)
        do i = 1, n
            if (i == 1) then
                x(i) = ( b(i) ) / L(i,i)
            else
                x(i) = ( b(i) - dot_product( x(1:i-1), L(i, 1:i-1) ) ) / L(i,i)
            end if
        end do

    end function forwards_subst_kind64

    pure function gaussian_elim_kind32(A, b) result(x)

        real(real32), intent(in), dimension(:,:) :: A
        real(real32), intent(in), dimension(:) :: b
        real(real32), dimension(size(A,1),size(A,2)) :: L, U, P, temp_u, temp_l, temp_p
        real(real32), dimension(size(b)) :: x, y
        integer :: k, i, j, n

        n = size(A,1)
        U = A
        temp_u = U
        L = identity_matrix(n)
        temp_l = L
        P = identity_matrix(n)
        temp_p = P
        
        do k = 1, n-1
            i = maxloc( U(:, k), dim=1 )
            U(k, k:n) = U(i, k:n)
            U(i, k:n) = temp_u(k, k:n)
            L(k, 1:k-1) = L(i, 1:k-1)
            L(i, 1:k-1) = temp_L(k, 1:k-1)
            P(k, :) = P(i, :)
            P(i, :) = temp_P(k, :)
            do j = k+1, n
                L(j,k) = U(j,k) / U(k,k)
                U(j, k:n) = U(j, k:n) - L(j,k) * U(k, k:n)
            end do
        end do

        y = forwards_subst(L, matmul(P, b))
        x = backwards_subst(U, y)

    end function gaussian_elim_kind32

    pure function gaussian_elim_kind64(A, b) result(x)

        real(real64), intent(in), dimension(:,:) :: A
        real(real64), intent(in), dimension(:) :: b
        real(real64), dimension(size(A,1),size(A,2)) :: L, U, P, temp_u, temp_l, temp_p
        real(real64), dimension(size(b)) :: x, y
        integer :: k, i, j, n

        n = size(A,1)
        U = A
        temp_u = U
        L = identity_matrix(n)
        temp_l = L
        P = identity_matrix(n)
        temp_p = P
        
        do k = 1, n-1
            i = maxloc( U(:, k), dim=1 )
            U(k, k:n) = U(i, k:n)
            U(i, k:n) = temp_u(k, k:n)
            L(k, 1:k-1) = L(i, 1:k-1)
            L(i, 1:k-1) = temp_L(k, 1:k-1)
            P(k, :) = P(i, :)
            P(i, :) = temp_P(k, :)
            do j = k+1, n
                L(j,k) = U(j,k) / U(k,k)
                U(j, k:n) = U(j, k:n) - L(j,k) * U(k, k:n)
            end do
            temp_U = U
            temp_L = L
            temp_P = P
        end do

        y = forwards_subst(L, matmul(P,b))
        x = backwards_subst(U, y)

    end function gaussian_elim_kind64

    pure function upper_hessenberg_kind32(A) result(H)

        real(real32), intent(in), dimension(:,:) :: A
        real(real32), dimension(size(A,1), size(A,2)) :: Q, H, ID
        real(real32), dimension(:), allocatable :: u
        real(real32) :: alpha
        integer :: i, n, j

        n = size(A,1)
        ID = identity_matrix(n)
        H = A
        do i = 1, n-2
            allocate( u(n - i) )
            alpha = norm2( H(i+1:n,i) )
            if (H(i+1,i) < 0) alpha = - alpha
            u = H(i+1:n, i) + alpha * ID(i+1:n, i+1) 
            if (norm2(u) /= 0) u = u / norm2( u )
            Q = ID
            Q(i+1:n, i+1:n) = Q(i+1:n, i+1:n) - 2 * vec_outer_product(u, u)
            H(i+1:n, i:n) = matmul(Q(i+1:n, i+1:n), H(i+1:n, i:n))
            H(1:n, i+1:n) = matmul(H(1:n, i+1:n), transpose(Q(i+1:n, i+1:n)) )
            deallocate(u)
        end do

    end function upper_hessenberg_kind32

    pure function upper_hessenberg_kind64(A) result(H)

        real(real64), intent(in), dimension(:,:) :: A
        real(real64), dimension(size(A,1), size(A,2)) :: Q, H, ID
        real(real64), dimension(:), allocatable :: u
        real(real64) :: alpha
        integer :: i, n, j

        n = size(A,1)
        ID = identity_matrix(n)
        H = A
        do i = 1, n-2
            allocate( u(n - i) )
            alpha = norm2( H(i+1:n,i) )
            if (H(i+1,i) < 0) alpha = - alpha
            u = H(i+1:n, i) + alpha * ID(i+1:n, i+1) 
            if (norm2(u) /= 0) u = u / norm2( u )
            Q = ID
            Q(i+1:n, i+1:n) = Q(i+1:n, i+1:n) - 2 * vec_outer_product(u, u)
            H(i+1:n, i:n) = matmul(Q(i+1:n, i+1:n), H(i+1:n, i:n))
            H(1:n, i+1:n) = matmul(H(1:n, i+1:n), transpose(Q(i+1:n, i+1:n)) )
            deallocate(u)
        end do

    end function upper_hessenberg_kind64

    pure function shifted_iteration_kind32(A, lambda) result(R)

        real(real32), intent(in), dimension(:,:) :: A, lambda
        real(real32), dimension(size(lambda,1), size(lambda,2)) :: R, ID
        real(real32), dimension(size(lambda,1)) :: alpha, eig_vec, w
        integer :: i, k, n

        n = size(lambda,1)
        ID = identity_matrix(n)
        do i = 1, n
            alpha(i) = lambda(i,i) + lambda(i,i)*0.02
        end do

        do i = 1, n
            eig_vec = 1
            do k = 1, n**2
                w = gaussian_elim((A - (alpha(i)) * ID), eig_vec)
                eig_vec = w / norm2(w)
            end do
            R(:,i) = eig_vec
        end do

    end function shifted_iteration_kind32

    pure function shifted_iteration_kind64(A, lambda) result(R)

        real(real64), intent(in), dimension(:,:) :: A, lambda
        real(real64), dimension(size(lambda,1), size(lambda,2)) :: R, ID
        real(real64), dimension(size(lambda,1)) :: alpha, eig_vec, w
        integer :: i, k, n

        n = size(lambda,1)
        ID = identity_matrix(n)
        do i = 1, n
            alpha(i) = lambda(i,i) + lambda(i,i)*0.05
        end do

        do i = 1, n
            eig_vec(i) = 1
            do k = 1, 2
                w = gaussian_elim((A - alpha(i) * ID), eig_vec)
                eig_vec = w / norm2(w)
            end do
            R(:,i) = eig_vec / norm2(eig_vec)
        end do

    end function shifted_iteration_kind64

! the following are not tested and finished yet:

    pure function GMRES(A, b) result(x)

        real(real64), intent(in), dimension(:,:) :: A
        real(real64), intent(in), dimension(:) :: b
        real(real64), dimension(:,:), allocatable :: Q, H, temp1, temp2, ID
        real(real64), dimension(:), allocatable :: x, y, u
        integer :: i, n

        n = size(b)
        ID = identity_matrix(n)

        allocate(u(n), x(n), y(n))

        u = b/norm2(b)

        do i = 1, 100
            allocate( Q(n, i), H(n, i))

            call arnoldi_iteration(A, i, Q, H)
            call qr_decomposition(H, temp1, temp2)
            y = least_squares_qr(H, b*ID(i:n, i))
            x = matmul(Q, y)

            deallocate(Q, H)
        end do

    end function GMRES

    pure function least_squares_qr(A, b) result(x)

        real(real64), intent(in), dimension(:,:) :: A
        real(real64), intent(in), dimension(:) :: b
        real(real64), dimension(:,:), allocatable :: Q, R
        real(real64), dimension(:), allocatable :: x
        integer :: i, m, n, k

        m = size(A,1)
        n = size(A,2)

        allocate( x(n), Q(n, m), R(n,n))

        call qr_decomposition(A, Q, R)
        x = gaussian_elim(transpose(Q), b)
        x = matmul(R, x)

    end function least_squares_qr


end module mod_linalg