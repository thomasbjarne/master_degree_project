module mod_io

    use iso_fortran_env, only: real32, real64
    implicit none

    interface write_to_file
        module procedure :: write_to_file_1d_kind32, write_to_file_1d_kind64, &
             write_to_file_2d_kind32, write_to_file_2d_kind64
    end interface write_to_file

contains

    subroutine write_to_file_2d_kind32(data, data_slot)
        
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

    end subroutine write_to_file_2d_kind32

    subroutine write_to_file_2d_kind64(data, data_slot)
        
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

    end subroutine write_to_file_2d_kind64

    subroutine write_to_file_1d_kind32(data, data_slot)
        
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

    end subroutine write_to_file_1d_kind32

    subroutine write_to_file_1d_kind64(data, data_slot)
        
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

    end subroutine write_to_file_1d_kind64

end module mod_io