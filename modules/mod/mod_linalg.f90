module mod_linalg

    implicit none

contains

    pure subroutine arnoldi_iteration(A, n, Q, H)

        real, intent(in), dimension(:,:) :: A
        integer, intent(in) :: n
        real, dimension(size(A,1)) :: b

        real, intent(out), dimension(size(A,1), n+1) :: Q
        real, intent(out), dimension(n+1, n) :: H
        real :: eps
        integer :: i, j

        b(:) = 1
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
            if (H(j+1,j) < eps) then
                exit
            else
                continue
            end if
        end do

    end subroutine arnoldi_iteration

    pure subroutine qr_decomposition(A, Q, R)

        real, intent(in), dimension(:,:) :: A
        real, intent(out), dimension(size(A,1), size(A,2)) :: Q, R
        real, dimension(size(A, 1), size(A,2)) :: ID, A_prime, Q_i
        real, dimension(:), allocatable :: u
        integer :: i, n, m

        n = size(A,1)

        A_prime = A
        ID = 0
        do i = 1, n
            ID(i,i) = 1
        end do

        do i = 1, n-1
            allocate( u(n - i + 1) )
            u = A_prime(i:n, i) - norm2( A_prime(i:n,i) ) * ID(i:n, i) 
            u = u / norm2( u )

            Q_i = ID
            Q_i(i:n, i:n) = Q_i(i:n, i:n) - 2 * vec_outer_product(u, u)

            A_prime = matmul(Q_i, A_prime)

            if (i == 1) then
                Q = transpose(Q_i)
                R = matmul(Q_i, A)
            else
                Q = matmul(Q, transpose(Q_i))
                R = matmul(Q_i, R)
            end if
            deallocate(u)
        end do

    end subroutine qr_decomposition

    pure function qr_algorithm(A) result(A_schur)

        real, intent(in), dimension(:,:) :: A
        real, dimension(size(A,1), size(A,2)) :: A_schur

        real, dimension(size(A,1), size(A,2)) :: A_k, Q_k, R_k
        integer :: i

        A_k = upper_hessenberg(A)

        do i = 1, 1000
            call qr_decomposition(A_k, Q_k, R_k)
            A_k = matmul(R_k, Q_k)
        end do

        A_schur = matmul( matmul(transpose(Q_k) , A_k), Q_k )

    end function qr_algorithm

    pure function vec_outer_product(u, v) result(A)

        real, intent(in), dimension(:) :: u, v
        real, dimension(size(u),size(v)) :: A

        integer :: i, j

        do i = 1, size(u)
            do j = 1, size(v)
                A(i, j) = u(i) * v(j)
            end do
        end do

    end function vec_outer_product

    pure function backwards_subst(R, b) result(x)

        real, intent(in), dimension(:,:) :: R
        real, intent(in), dimension(size(R,2)) :: b
        real, dimension(size(R,2)) :: x, row_Rx
        integer :: i

        do i = size(x), 1, -1
            if (i == size(x)) then
                x(i) = ( b(i) ) / R(i,i)
            else
                row_Rx(i+1:size(x)) = R(i, i+1:size(x)) * x(i+1:size(x)) 
                x(i) = ( b(i)  - sum(row_Rx) ) / R(i,i)
            end if
        end do

    end function backwards_subst

    pure function linear_system_solver(A, b) result(x)

        real, intent(in), dimension(:,:) :: A
        real, intent(in), dimension(:) :: b
        real, dimension(size(A,1),size(A,2)) :: Q, R
        real, dimension(size(b)) :: x, y
        integer :: i

        call qr_decomposition(A, Q, R)
        y = backwards_subst(R, b)

        do i = 1, size(x)
            x(i) = dot_product(Q(i,:), y)
        end do

    end function linear_system_solver

    pure function rayleigh_ritz(A) result(eig_pair)

        real, intent(in), dimension(:,:) :: A
        real, dimension(:,:), allocatable :: eig_pair, V, V_star

        allocate(V(size(A,1), size(A,2)-1), V_star(size(A,1)-1, size(A,2)))
        allocate(eig_pair(size(A,1), 2))


    end function rayleigh_ritz

    pure function upper_hessenberg(A) result(H)

        real, intent(in), dimension(:,:) :: A
        real, dimension(size(A,1),size(A,2)) :: H
        real, dimension(:), allocatable :: x, v, e
        real, dimension(:,:), allocatable :: v_transpose
        integer :: i, n

        n = size(A,1)
        H = A

        do i = 1, n-2

            allocate(x(n-i), v(n-i), e(n-i), v_transpose(1,n-i))

            e = 0
            e(1) = 1

            x = A(i+1:n, i)
            if (x(1) >= 0) then
                v = norm2(x) * e + x
            else if (x(1) < 0) then
                v = -norm2(x) * e + x
            end if
            v = v / norm2(v)
            v_transpose = vector_transpose(v)

            H(i+1:n, i:n) = H(i+1:n, i:n) - 2 * vec_outer_product(v, matmul(v, H(i+1:n, i:n)))
            H(1:n, i+1:n) = H(1:n, i+1:n) - 2 * vec_outer_product(matmul(H(1:n, i+1:n), v), v)

            deallocate(x, v, e, v_transpose)

        end do

    end function upper_hessenberg

    pure function vector_transpose(v) result(v_t)
        
        real, intent(in), dimension(:) :: v
        real, dimension(1,size(v)) :: v_t

        v_t(1,:) = v

    end function vector_transpose

end module mod_linalg