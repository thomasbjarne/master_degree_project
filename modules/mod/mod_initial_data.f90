module mod_initial_data

    implicit none

    interface initial_func
        module procedure :: initial_func_1D, initial_func_2D
    end interface initial_func

contains

    pure real function initial_func_1D(x) result(func)

        real, intent(in) :: x

        func = sin(x)

    end function initial_func_1D

    pure real function initial_func_2D(x, y) result(func)

        real, intent(in) :: x, y

        func = sin(x + y)

    end function initial_func_2D

end module mod_initial_data