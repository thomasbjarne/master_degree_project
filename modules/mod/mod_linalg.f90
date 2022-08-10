module mod_linalg

    implicit none

contains

    pure function identity_matrix(n) result(I)

        integer, intent(in) :: n
        real, dimension(n,n) :: I
        integer :: k

        I = 0
        do k = 1, n
            I(k,k) = 1
        end do

    end function identity_matrix

    pure subroutine arnoldi_iteration(A, n, Q, H)

        real,    intent(in), dimension(:,:) :: A
        integer, intent(in) :: n
        real, intent(out), dimension(size(A,1), n+1) :: Q
        real, intent(out), dimension(n+1, n) :: H
        real, dimension(size(A,1)) :: b
        real :: eps
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
        real, dimension(size(A, 1), size(A,2)) :: ID, Q_i
        real, dimension(:), allocatable :: u
        real :: alpha
        integer :: i, n

        n = size(A,1)

        ID = identity_matrix(n)
        Q = ID
        R = A
        do i = 1, n
            allocate( u(n - i + 1) )

            alpha = -norm2( R(i:n,i) )
            if (R(i,i) < 0) alpha = - alpha
            
            u = R(i:n, i) + alpha * ID(i:n, i) 
            if (norm2(u) /= 0) u = u / norm2( u )

            Q_i = ID
            Q_i(i:n, i:n) = Q_i(i:n, i:n) - 2 * vec_outer_product(u,u)
            
            Q = matmul(Q, transpose(Q_i))
            R(i:n, i:n) = R(i:n, i:n) - 2 * vec_outer_product(u, matmul(u, R(i:n, i:n)))

            deallocate(u)
        end do

    end subroutine qr_decomposition

    pure function qr_algorithm(A) result(A_schur) !Fix this

        real, intent(in), dimension(:,:) :: A
        real, dimension(size(A,1), size(A,2)) :: A_schur, A_k, Q_k, R_k
        real, parameter :: eps = 1e-10
        integer :: i, j, k, n

        n = size(A_k, 1)

        A_k = upper_hessenberg(A)

        do k = 1, 10000
            A_k = A_k - A_k(n, n) * identity_matrix(n)
            call qr_decomposition(A_k, Q_k, R_k)
            A_k = matmul(R_k, Q_k) + A_k(n, n) * identity_matrix(n)
            do i = 2, n-1
                do j = i+1, n
                    if ( A_k(i,j) <= eps ) then
                        A_k(i, j) = 0
                        A_k(j, i) = 0
                    end if
                end do
            end do
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

        x = matmul(Q, y)

    end function linear_system_solver

    pure function upper_hessenberg(A) result(H)

        real, intent(in), dimension(:,:) :: A
        real, dimension(size(A,1), size(A,2)) :: Q, H, ID
        real, dimension(:), allocatable :: u
        real :: alpha
        integer :: i, n

        n = size(A,1)

        ID = identity_matrix(n)
        H = A

        do i = 1, n-2

            allocate( u(n - i) )

            alpha = norm2( A(i+1:n,i) )
            if (A(i+1,i) < 0) alpha = - alpha

            u = A(i+1:n, i) + alpha * ID(i+1:n, i+1) 
            u = u / norm2( u )

            Q = ID
            Q(i+1:n, i+1:n) = Q(i+1:n, i+1:n) - 2 * vec_outer_product(u, u)
            
            H(i+1:n, i:n) = matmul(Q(i+1:n, i+1:n), A(i+1:n, i:n))
            H(1:n, i+1:n) = matmul(A(1:n, i+1:n), transpose(Q(i+1:n, i+1:n)) )

            deallocate(u)

        end do

    end function upper_hessenberg

    pure function inverse_iteration(A, lambda) result(R)

        real, intent(in), dimension(:,:) :: A, lambda
        real, dimension(size(lambda,1), size(lambda,2)) :: R, temp, ID
        real, dimension(size(lambda,1)) :: eig_val_list, eig_vec, w
        real, parameter :: delta = 1e-10
        integer :: i, k, n

        n = size(lambda,1)
        ID = 0
        do i = 1, n
            eig_val_list(i) = lambda(i,i)
            ID(i,i) = 1
        end do

        eig_vec = 1
        eig_vec = eig_vec / norm2(eig_vec)
        do i = 1, n
            temp = A - (eig_val_list(i) + delta) * ID
            do k = 1, 1000
                w = linear_system_solver(temp, eig_vec)
                eig_vec = w / norm2(w)
            end do
            R(:,i) = eig_vec
        end do

    end function inverse_iteration

end module mod_linalg