module mod_finite_diff

    implicit none

    interface central_diff
        module procedure :: central_diff_1d, central_diff_2d
    end interface central_diff

contains

    pure function central_diff_1d(v, periodic, h) result(dv)

        logical, intent(in) :: periodic
        real, dimension(:), intent(in) :: v
        real, dimension(size(v)) :: dv
        real, optional, intent(in) :: h
		
        if (periodic .eqv. .true.) then
            dv(1)           = v(2) - v(size(v))
            dv(2:size(v)-1) = v(3:size(v)) - v(1:size(v) - 2)
            dv(size(v))     = v(1) - v(size(v) - 1)
        else if (periodic .eqv. .false.) then !ignore boundary points..
            dv(1)           = 0
            dv(2:size(v)-1) = v(3:size(v)) - v(1:size(v) - 2)
            dv(size(v))     = 0
        end if

        dv = dv * 0.5

        if (present(h) .eqv. .true.) dv = dv / h

    end function central_diff_1d

    pure function central_diff_2d(v, periodic, h) result(dv)

        logical, intent(in) :: periodic
        real, dimension(:,:), intent(in) :: v
        real, dimension(size(v, 1), size(v, 2)) :: dx, dy
        real, dimension(2, size(v, 1), size(v, 2)) :: dv
        real, optional, intent(in) :: h

        if (periodic .eqv. .true.) then
            dx(:, 1)                = v(:, 2) - v(:, size(v, 2))
            dx(:, 2:size(v, 2) - 1) = v(:, 3:size(v, 2)) &
                                            - v(:, 1:size(v, 2) - 2)
            dx(:, size(v, 2))       = v(:, 1) &
                                            - v(:, size(v, 2) - 1)

            dy(1, :)                = v(2, :) - v(size(v, 1), :)
            dy(2:size(v, 1) - 1, :) = v(3:size(v, 1)-1, :) &
                                            - v(1:size(v, 1) - 2, :)
            dy(size(v, 1), :)       = v(1, :) - v(size(v, 1) - 1, :)
        else if (periodic .eqv. .false.) then !ignore all boundary points
            dx(:, 1)                    = 0
            dx(:, 2:size(v, 2) - 1)     = v(:, 3:size(v, 2)) &
                                            - v(:, 1:size(v, 2) - 2)
            dx(:, size(v, 2))           = 0

            dy(1, :)                    = 0
            dy(2:size(v, 1) - 1, :)     = v(3:size(v, 1)-1, :) &
                                            - v(1:size(v, 1) - 2, :)
            dy(size(v, 1), :)           = 0
        end if

        dv(1, :, :) = dx
        dv(2, :, :) = dy

        dv = dv * 0.5

        if (present(h) .eqv. .true.) dv = dv / h

    end function central_diff_2d

end module mod_finite_diff