module mod_io

    implicit none

    interface write_to_file
        module procedure :: write_to_file_1d, write_to_file_2d
    end interface write_to_file

contains

    subroutine write_to_file_2d(data, data_slot)
        
        real, dimension(:, :), intent(in) :: data
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

    end subroutine write_to_file_2d

    subroutine write_to_file_1d(data, data_slot)
        
        real, dimension(:), intent(in) :: data
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

    end subroutine write_to_file_1d

end module mod_io