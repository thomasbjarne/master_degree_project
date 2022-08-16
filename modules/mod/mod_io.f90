module mod_io

    use iso_fortran_env, only: real32, real64, int8, int16
    implicit none

    interface write_to_file
        module procedure :: write_to_file_1d_kind32_real, write_to_file_1d_kind64_real, &
             write_to_file_2d_kind32_real, write_to_file_2d_kind64_real, &
             write_to_file_1d_kind16_int, write_to_file_1d_kind8_int, &
             write_to_file_2d_kind8_int, write_to_file_2d_kind16_int
    end interface write_to_file

contains

!write to file 2D

    subroutine write_to_file_2d_kind8_int(data, data_slot)
            
        integer(int8), dimension(:, :), intent(in) :: data
        integer, intent(in) :: data_slot
        character(len=20) :: filename
        character(len=12) :: str
        integer :: i, fileunit

        write (str, *) data_slot
        str = adjustl(str)
        filename = 'data/datafile' // trim(str) // '.txt'
        filename = trim(filename)

        open(newunit=fileunit, file=filename, position='rewind')

        do i = 1, size(data, 1)
            write(unit=fileunit, fmt=*) data(i, :)
        end do

        close(unit=fileunit)

    end subroutine write_to_file_2d_kind8_int

    subroutine write_to_file_2d_kind16_int(data, data_slot)
            
        integer(int16), dimension(:, :), intent(in) :: data
        integer, intent(in) :: data_slot
        character(len=20) :: filename
        character(len=12) :: str
        integer :: i, fileunit

        write (str, *) data_slot
        str = adjustl(str)
        filename = 'data/datafile' // trim(str) // '.txt'
        filename = trim(filename)

        open(newunit=fileunit, file=filename, position='rewind')

        do i = 1, size(data, 1)
            write(unit=fileunit, fmt=*) data(i, :)
        end do

        close(unit=fileunit)

    end subroutine write_to_file_2d_kind16_int

    subroutine write_to_file_2d_kind32_real(data, data_slot)
        
        real(real32), dimension(:, :), intent(in) :: data
        integer, intent(in) :: data_slot
        character(len=20) :: filename
        character(len=12) :: str
        integer :: i, fileunit

        write (str, *) data_slot
        str = adjustl(str)
        filename = 'data/datafile' // trim(str) // '.txt'
        filename = trim(filename)

        open(newunit=fileunit, file=filename, position='rewind')

        do i = 1, size(data, 1)
            write(unit=fileunit, fmt=*) data(i, :)
        end do

        close(unit=fileunit)

    end subroutine write_to_file_2d_kind32_real

    subroutine write_to_file_2d_kind64_real(data, data_slot)
        
        real(real64), dimension(:, :), intent(in) :: data
        integer, intent(in) :: data_slot
        character(len=20) :: filename
        character(len=12) :: str
        integer :: i, fileunit

        write (str, *) data_slot
        str = adjustl(str)
        filename = 'data/datafile' // trim(str) // '.txt'
        filename = trim(filename)

        open(newunit=fileunit, file=filename, position='rewind')

        do i = 1, size(data, 1)
            write(unit=fileunit, fmt=*) data(i, :)
        end do

        close(unit=fileunit)

    end subroutine write_to_file_2d_kind64_real

!write to file 1D

    subroutine write_to_file_1d_kind8_int(data, data_slot)
            
        integer(int8), dimension(:), intent(in) :: data
        integer, intent(in) :: data_slot
        character(len=20) :: filename
        character(len=12) :: str
        integer :: i, fileunit

        write (str, *) data_slot
        str = adjustl(str)
        filename = 'data/datafile' // trim(str) // '.txt'
        filename = trim(filename)

        open(newunit=fileunit, file=filename, position='rewind')

        do i = 1, size(data, 1)
            write(unit=fileunit, fmt=*) data(i)
        end do

        close(unit=fileunit)

    end subroutine write_to_file_1d_kind8_int

    subroutine write_to_file_1d_kind16_int(data, data_slot)
            
        integer(int16), dimension(:), intent(in) :: data
        integer, intent(in) :: data_slot
        character(len=20) :: filename
        character(len=12) :: str
        integer :: i, fileunit

        write (str, *) data_slot
        str = adjustl(str)
        filename = 'data/datafile' // trim(str) // '.txt'
        filename = trim(filename)

        open(newunit=fileunit, file=filename, position='rewind')

        do i = 1, size(data, 1)
            write(unit=fileunit, fmt=*) data(i)
        end do

        close(unit=fileunit)

    end subroutine write_to_file_1d_kind16_int

    subroutine write_to_file_1d_kind32_real(data, data_slot)
        
        real(real32), dimension(:), intent(in) :: data
        integer, intent(in) :: data_slot
        character(len=20) :: filename
        character(len=12) :: str
        integer :: i, fileunit

        write (str, *) data_slot
        str = adjustl(str)
        filename = 'data/datafile' // trim(str) // '.txt'
        filename = trim(filename)

        open(newunit=fileunit, file=filename, position='rewind')

        do i = 1, size(data, 1)
            write(unit=fileunit, fmt=*) data(i)
        end do

        close(unit=fileunit)

    end subroutine write_to_file_1d_kind32_real

    subroutine write_to_file_1d_kind64_real(data, data_slot)
        
        real(real64), dimension(:), intent(in) :: data
        integer, intent(in) :: data_slot
        character(len=20) :: filename
        character(len=12) :: str
        integer :: i, fileunit

        write (str, *) data_slot
        str = adjustl(str)
        filename = 'data/datafile' // trim(str) // '.txt'
        filename = trim(filename)

        open(newunit=fileunit, file=filename, position='rewind')

        do i = 1, size(data, 1)
            write(unit=fileunit, fmt=*) data(i)
        end do

        close(unit=fileunit)

    end subroutine write_to_file_1d_kind64_real

end module mod_io