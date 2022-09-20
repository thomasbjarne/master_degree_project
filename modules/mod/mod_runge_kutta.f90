module mod_runge_kutta

    implicit none

    interface rk4
        module procedure :: rk4_1D, rk4_2D
    end interface rk4

contains

    pure subroutine rk4_1D(dydt, u_0, k, tmin, tmax, u)

        real, intent(in) :: tmin, tmax, k
        real, intent(in), dimension(:) :: u_0
        real, dimension(:), allocatable :: a_1, a_2, a_3, a_4, time
        real, intent(out), dimension(:, :), allocatable :: u
        integer :: i, n

		interface
			pure function dydt(t,y) result(res)
				real, intent(in), dimension(:) :: y
				real, intent(in), optional :: t
                real, dimension(size(y)) :: res
			end function dydt
		end interface

        n = nint((tmax - tmin) / k)
        
        allocate(a_1(n), a_2(n), a_3(n), a_4(n), time(n))
        allocate(u(n, size(u_0)))
        u(1, :) = u_0

        do i = 1, n - 1
            a_1 = dydt(time(i), u(i, :))
            a_2 = dydt(time(i), u(i, :) + k*a_1/2)
            a_3 = dydt(time(i), u(i, :) + k*a_2/2)
            a_4 = dydt(time(i), u(i, :) + k*a_3)

            u(i+1, :) = u(i, :) + k*(a_1 + 2*a_2 + 2*a_3 + a_4) / 6
        end do

        deallocate( a_1, a_2, a_3, a_4, time )

    end subroutine rk4_1D

    pure subroutine rk4_2D(dydt, u_0, k, tmin, tmax, u)

        real, intent(in) :: tmin, tmax, k
        real, intent(in), dimension(:,:) :: u_0
        real, dimension(:), allocatable :: time
        real, dimension(:,:), allocatable :: a_1, a_2, a_3, a_4
        real, intent(out), dimension(:, :, :), allocatable :: u
        integer :: i, n

		interface
			pure function dydt(t,y) result(res)
				real, intent(in), dimension(:,:) :: y
				real, intent(in), optional :: t
                real, dimension(size(y), size(y)) :: res
			end function dydt
		end interface

        n = nint((tmax - tmin) / k)
        
        allocate(a_1(n,n), a_2(n,n), a_3(n,n), a_4(n,n), time(n))
        allocate(u(n, size(u_0,1), size(u_0,2)))
        u(1, :, :) = u_0

        do i = 1, n - 1
            a_1 = dydt(time(i), u(i, :, :))
            a_2 = dydt(time(i), u(i, :, :) + k*a_1/2)
            a_3 = dydt(time(i), u(i, :, :) + k*a_2/2)
            a_4 = dydt(time(i), u(i, :, :) + k*a_3)

            u(i+1, :, :) = u(i, :, :) + k*(a_1 + 2*a_2 + 2*a_3 + a_4) / 6
        end do

        deallocate( a_1, a_2, a_3, a_4, time )

    end subroutine rk4_2D


end module mod_runge_kutta
