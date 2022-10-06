module mod_finite_volume

    implicit none

contains
!SEMI DISCRETE SCHEMES: Args must be t and y to pass the schemes to rk4.

    pure function upwind_scheme(t, y) result(res)
            
        real, intent(in), dimension(:) :: y
        real, intent(in), optional :: t
        real, dimension(size(y)) :: res
        integer :: i
        real, parameter :: a = -1.
        real, parameter :: h = 1./500

        if (a > 0) then
            res(1) = -(a/h)*(y(1) - y(size(y)))
            do i = 2, size(y)
                res(i) = -(a/h)*(y(i)-y(i-1))
            end do
        else if (a < 0) then
            do i = 1, size(y)-1
                res(i) = -(a/h)*(y(i+1) - y(i))
            end do
            res(size(y)) = -(a/h)*(y(1) - y(size(y)))
        end if

    end function upwind_scheme

!FULLY DISCRETE SCHEMES

    pure function fully_discrete_upwind_scheme(h, speed, tmin, tmax, Q_0) result(Q)

        real, intent(in) :: h, speed, tmin, tmax
        real, intent(in), dimension(:) :: Q_0
        real, dimension(:,:), allocatable :: Q
        real :: k
        integer :: i, n
        n = size(Q_0)
        k = abs(tmax - tmin)/n
        allocate( Q(n, n) )

        Q(1,:) = Q_0 
        if (speed > 0) then
            do i = 1, n-1
                Q(i+1, 1) = Q(i, 1) - speed*(k/h)*(Q(i, 1) - Q(i, n))
                Q(i+1, 2:n) = Q(i,2:n) - speed*(k/h)*(Q(i, 2:n) - Q(i, 1:n-1))
            end do
        else
            do i = 1, n-1
                Q(i+1, 1:n-1) = Q(i,1:n-1) - speed*(k/h)*(Q(i, 2:n) - Q(i, 1:n-1))
                Q(i+1, n) = Q(i, n) - speed*(k/h)*(Q(i, 1) - Q(i, n))
            end do
        end if

    end function fully_discrete_upwind_scheme

    pure function fully_discrete_naive_central_scheme(h, speed, tmin, tmax, Q_0) result(Q)

        real, intent(in) :: h, speed, tmin, tmax
        real, intent(in), dimension(:) :: Q_0
        real, dimension(:,:), allocatable :: Q
        real :: k
        integer :: i, n
        n = size(Q_0)
        k = abs(tmax - tmin)/n
        allocate( Q(n, n) )

        Q(1,:) = Q_0
        do i = 1, n-1
            Q(i+1,1) = Q(i,1) - (0.5)*((speed*k)/h)*(Q(i,2)-Q(i,n))
            Q(i+1,2:n-1) = Q(i,2:n-1) - (0.5)*((speed*k)/(h))*(Q(i,3:n)-Q(i,1:n-2))
            Q(i+1,n) = Q(i,n) - (0.5)*((speed*k)/(h))*(Q(i,1)-Q(i,n-1))
        end do

    end function fully_discrete_naive_central_scheme

    pure function fully_discrete_Lax_Friedrichs_scheme(h, speed, tmin, tmax, Q_0) result(Q)

        real, intent(in) :: h, speed, tmin, tmax
        real, intent(in), dimension(:) :: Q_0
        real, dimension(:,:), allocatable :: Q
        real :: k
        integer :: i, n
        n = size(Q_0)
        k = abs(tmax - tmin)/n
        allocate( Q(n, n) )

        Q(1,:) = Q_0
        do i = 1, n-1
            Q(i+1,1) = (0.5)*(Q(i,2)+Q(i,n)) - ((speed*k)/(2*h))*(Q(i,2)-Q(i,n))
            Q(i+1,2:n-1) = (0.5)*(Q(i,3:n)+Q(i,1:n-2)) - ((speed*k)/(2*h))*(Q(i,3:n)-Q(i,1:n-2))
            Q(i+1,n) = (0.5)*(Q(i,1)+Q(i,n-1)) - ((speed*k)/(2*h))*(Q(i,1)-Q(i,n-1))
        end do

    end function fully_discrete_Lax_Friedrichs_scheme

    pure function fully_discrete_Richtmyer_LW_scheme(h, speed, tmin, tmax, Q_0) result(Q)

        real, intent(in) :: h, speed, tmin, tmax
        real, intent(in), dimension(:) :: Q_0
        real, dimension(:,:), allocatable :: Q, Q_mid
        real :: k, k_2, h_2
        integer :: i, n
        n = size(Q_0)
        k = abs(tmax - tmin)/n
        h_2 = h/2
        k_2 = k/2
        allocate( Q(n, n), Q_mid(n,n))

        Q(1,:) = Q_0
        Q_mid(1,:) = Q_0
        do i = 1, n-1
            Q_mid(i+1,1) = (0.5)*(Q_mid(i,2)+Q_mid(i,n)) - ((speed*k_2)/(2*h_2))*(Q_mid(i,2)-Q_mid(i,n))
            Q_mid(i+1,2:n-1) = (0.5)*(Q_mid(i,3:n)+Q_mid(i,1:n-2)) - ((speed*k_2)/(2*h_2))*(Q_mid(i,3:n)-Q_mid(i,1:n-2))
            Q_mid(i+1,n) = (0.5)*(Q_mid(i,1)+Q_mid(i,n-1)) - ((speed*k_2)/(2*h_2))*(Q_mid(i,1)-Q_mid(i,n-1))
        end do
        do i = 1, n-1
            Q(i+1,1) = (0.5)*(Q(i,2)+Q(i,n)) - ((speed*k)/(2*h))*(Q(i,2)-Q_mid(i,n))
            Q(i+1,2:n-1) = (0.5)*(Q(i,3:n)+Q(i,1:n-2)) - ((speed*k)/(2*h))*(Q(i,3:n)-Q_mid(i,1:n-2))
            Q(i+1,n) = (0.5)*(Q(i,1)+Q(i,n-1)) - ((speed*k)/(2*h))*(Q(i,1)-Q_mid(i,n-1))
        end do

    end function fully_discrete_Richtmyer_LW_scheme

    pure function fully_discrete_second_order_LW_scheme(h, speed, tmin, tmax, Q_0) result(Q)

        real, intent(in) :: h, speed, tmin, tmax
        real, intent(in), dimension(:) :: Q_0
        real, dimension(:,:), allocatable :: Q
        real :: k
        integer :: i, n
        n = size(Q_0)
        k = abs(tmax - tmin)/n
        allocate( Q(n, n) )

        Q(1,:) = Q_0
        do i = 1, n-1
            Q(i+1, 1) = Q(i, 1) - (0.5)*((speed*k)/(h))*(Q(i, 2) - Q(i, n)) &
            + (0.5)*((speed**2)*(k**2)/(h**2))*(Q(i,2) - 2*Q(i, 1) + Q(i, n))
            Q(i+1, 2:n-1) = Q(i, 2:n-1) - (0.5)*((speed*k)/(h))*(Q(i, 3:n) - Q(i, 1:n-2)) &
            + (0.5)*((speed**2)*(k**2)/(h**2))*(Q(i,3:n) - 2*Q(i, 2:n-1) + Q(i, 1:n-2))
            Q(i+1, n) = Q(i, n) - (0.5)*((speed*k)/(h))*(Q(i, 1) - Q(i, n-1)) &
            + (0.5)*((speed**2)*(k**2)/(h**2))*(Q(i,1) - 2*Q(i, n) + Q(i, n-1))
        end do

    end function fully_discrete_second_order_LW_scheme

    pure function fully_discrete_second_order_LW_scheme_1dsys(h, A, tmin, tmax, Q_0) result(Q)

        real, intent(in) :: h, tmin, tmax
        real, intent(in), dimension(:,:) :: A
        real, intent(in), dimension(:) :: Q_0
        real, dimension(:,:), allocatable :: Q
        real :: k
        integer :: i, n, j, m, r
        n = size(Q_0)
        m = size(A,1)
        r = n/m
        k = abs(tmax - tmin)/n
        allocate( Q(n, n) )

        Q(1, :) = Q_0
        do i = 1, n-1
            Q(i+1, 1:n:r) = Q(i, 1:n:r) - matmul(A, (0.5)*(k/h)*(Q(i, 2:n:r) - Q(i, r:n:r))) &
            + (0.5)*matmul(matmul(A,A), ((k**2)/(h**2))*(Q(i,2:n:r) - 2*Q(i, 1:n:r) + Q(i, r:n:r)))
            do j = 2, r-1
                Q(i+1, j:n:r) = Q(i, j:n:r) - matmul(A, (0.5)*(k/h)*(Q(i, j+1:n:r) - Q(i, (j-1):n:r))) &
                + matmul(matmul(A, A), (0.5)*((k**2)/(h**2))*(Q(i, (j+1):n:r) - 2*Q(i, j:n:r) + Q(i, (j-1):n:r)))
            end do
            Q(i+1, r:n:r) = Q(i, r:n:r) - matmul(A,(0.5)*(k/h)*(Q(i, 1:n:r) - Q(i, ((n-1)/m):n:r))) &
            + matmul(matmul(A,A), (0.5)*((k**2)/(h**2))*(Q(i,1:n:r) - 2*Q(i, r:n:r) + Q(i, ((n-1)/m):n:r)))
        end do

    end function fully_discrete_second_order_LW_scheme_1dsys

end module mod_finite_volume