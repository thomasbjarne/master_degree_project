module mod_runge_kutta

    use mod_finite_diff
    implicit none

    interface rk4
        module procedure :: rk4_1D
    end interface rk4

contains

    subroutine rk4_1D(y_0, k, tmin, tmax, y)

        real, intent(in) :: tmin, tmax, k
        real, intent(in), dimension(:) :: y_0
        real, dimension(:), allocatable :: a_1, a_2, a_3, a_4, t
        real, intent(out), dimension(:, :), allocatable :: y
        integer :: i, n

        n = nint( (tmax - tmin) * k )
        
        allocate(a_1(n), a_2(n), a_3(n), a_4(n))
        allocate(y(n, size(y_0)))
        y(1, :) = y_0

        do i = 1, n - 1
            a_1 = dydt(t(i), y(i, :))
            a_2 = dydt(t(i), y(i, :) + k*a_1/2)
            a_3 = dydt(t(i), y(i, :) + k*a_2/2)
            a_4 = dydt(t(i), y(i, :) + k*a_3)

            y(i+1, :) = y(i, :) + k*(a_1 + 2*a_2 + 2*a_3 + a_4) / 6
        end do

    end subroutine rk4_1D

    pure function dydt(t, u) result(f)

        real, intent(in), dimension(:) :: u
        real, intent(in), optional :: t
        real, dimension(size(u)) :: f

        f = central_diff(u, periodic=.true.)

    end function dydt

end module mod_runge_kutta