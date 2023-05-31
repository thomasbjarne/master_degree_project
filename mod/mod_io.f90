module mod_io

    use iso_fortran_env, only: real32, real64, int8, int16, int32
    use mod_derived_types
    implicit none

	interface load_to_var
		module procedure :: load_to_var_1d_kind8_int, &
							load_to_var_1d_kind16_int, &
							load_to_var_1d_kind32_int, &
							load_to_var_1d_kind32_real, &
							load_to_var_1d_kind64_real, &
							load_to_var_2d_kind8_int, &
							load_to_var_2d_kind16_int, &
							load_to_var_2d_kind32_int, &
							load_to_var_2d_kind32_real, &
							load_to_var_2d_kind64_real
	end interface load_to_var

    interface write_to_file
        module procedure :: write_to_file_1d_logical, &
        					write_to_file_1d_kind8_int, &
        					write_to_file_1d_kind16_int, &
        					write_to_file_1d_kind32_int, &
        					write_to_file_1d_kind32_real, &
        					write_to_file_1d_kind64_real, &
        					write_to_file_2d_kind8_int, &
        					write_to_file_2d_kind16_int, &
        					write_to_file_2d_kind32_int, &
             				write_to_file_2d_kind32_real, &
             				write_to_file_2d_kind64_real
    end interface write_to_file

contains

!load grid

subroutine load_grid(grid)

	type(triangulation) :: grid
	
	call load_to_var('point_list.txt', grid%points, 2)
	call load_to_var('connectivity_list.txt', grid%connectivity_list, 3)
	call load_to_var('edge_list.txt', grid%edge_list, 2)
	call triangle_fix_orientation(grid%points, grid%connectivity_list)
	grid%barycenters  = triangle_barycenters(grid%points, &
											 grid%connectivity_list)
	grid%volumes 	  = triangle_volumes(grid%points, grid%connectivity_list)
	grid%dual_volumes = node_volumes(grid%points, grid%connectivity_list, &
									 grid%volumes)
	grid%neighbours	  = triangle_neighbours(grid%connectivity_list)
	grid%h			  = triangle_char_length(grid%points,&
											 grid%connectivity_list)
	grid%distancematrix = triangle_distances(grid%barycenters)
	
	
end subroutine load_grid

!read row count of file

function num_rows(filename) result(n)

	character(len=*), intent(in) 	:: filename
	integer 						:: n, io, fileunit

	n = 0
	open(newunit=fileunit, file=filename, position='rewind', action='read')
	do
		read(fileunit, *, iostat=io)
		if (io/=0) exit
		n = n + 1
	end do
	close(fileunit)

end function num_rows

!load to var 1d

subroutine load_to_var_1d_kind8_int(filename, var)

	character(len=*), intent(in) 							:: filename
	integer(int8), dimension(:), allocatable, intent(inout) :: var
	integer 												:: i, n, fileunit
	character(len=200)										:: str
	logical 												:: file_exists
	
	str = 'data/' // filename
	str = trim(str)
	inquire(file=str, exist=file_exists)
	
	if (file_exists) then
		n = num_rows(str)
		allocate(var(n))
		open(newunit=fileunit, file=str, position= 'rewind', action='read')
		do i = 1, size(var, 1)
			read(unit=fileunit, fmt=*, err=100) var(i) 
		end do
	else if (file_exists .eqv. .false.) then
		stop
	end if

	100 close(fileunit)
	
end subroutine load_to_var_1d_kind8_int

subroutine load_to_var_1d_kind16_int(filename, var)

	character(len=*), intent(in) 							 :: filename
	integer(int16), dimension(:), allocatable, intent(inout) :: var
	integer 												 :: i, n, fileunit
	character(len=200)										 :: str
	logical 												 :: file_exists
	
	str = 'data/' // filename
	str = trim(str)
	inquire(file=str, exist=file_exists)
	
	if (file_exists) then
		n = num_rows(str)
		allocate(var(n))
		open(newunit=fileunit, file=str, position= 'rewind', action='read')
		do i = 1, size(var, 1)
			read(unit=fileunit, fmt=*, err=100) var(i) 
		end do
	else if (file_exists .eqv. .false.) then
		stop
	end if

	100 close(fileunit)
	
end subroutine load_to_var_1d_kind16_int

subroutine load_to_var_1d_kind32_int(filename, var)

	character(len=*), intent(in) 							 :: filename
	integer(int32), dimension(:), allocatable, intent(inout) :: var
	integer 												 :: i, n, fileunit
	character(len=200)										 :: str
	logical 												 :: file_exists
	
	str = 'data/' // filename
	str = trim(str)
	inquire(file=str, exist=file_exists)
	
	if (file_exists) then
		n = num_rows(str)
		allocate(var(n))
		open(newunit=fileunit, file=str, position= 'rewind', action='read')
		do i = 1, size(var, 1)
			read(unit=fileunit, fmt=*, err=100) var(i) 
		end do
	else if (file_exists .eqv. .false.) then
		stop
	end if

	100 close(fileunit)
	
end subroutine load_to_var_1d_kind32_int

subroutine load_to_var_1d_kind32_real(filename, var)

	character(len=*), intent(in) 							:: filename
	real(real32), dimension(:), allocatable, intent(inout)  :: var
	integer 												:: i, n, fileunit
	character(len=200)										:: str
	logical 												:: file_exists
	
	str = 'data/' // filename
	str = trim(str)
	inquire(file=str, exist=file_exists)
	
	if (file_exists) then
		n = num_rows(str)
		allocate(var(n))
		open(newunit=fileunit, file=str, position= 'rewind', action='read')
		do i = 1, size(var, 1)
			read(unit=fileunit, fmt=*, err=100) var(i) 
		end do
	else if (file_exists .eqv. .false.) then
		stop
	end if

	100 close(fileunit)
	
end subroutine load_to_var_1d_kind32_real

subroutine load_to_var_1d_kind64_real(filename, var)

	character(len=*), intent(in) 							:: filename
	real(real64), dimension(:), allocatable, intent(inout)  :: var
	integer 												:: i, n, fileunit
	character(len=200)										:: str
	logical 												:: file_exists
	
	str = 'data/' // filename
	str = trim(str)
	inquire(file=str, exist=file_exists)
	
	if (file_exists) then
		n = num_rows(str)
		allocate(var(n))
		open(newunit=fileunit, file=str, position= 'rewind', action='read')
		do i = 1, size(var, 1)
			read(unit=fileunit, fmt=*, err=100) var(i) 
		end do
	else if (file_exists .eqv. .false.) then
		stop
	end if

	100 close(fileunit)
	
end subroutine load_to_var_1d_kind64_real

!load to var 2d

subroutine load_to_var_2d_kind8_int(filename, var, columns)

	character(len=*), intent(in) 								:: filename
	integer(int8), dimension(:,:), allocatable, intent(inout)  	:: var
	integer, intent(in)											:: columns
	integer 													:: i, n, &
																   fileunit
	character(len=200)											:: str
	logical 													:: file_exists
	
	str = 'data/' // filename
	str = trim(str)
	inquire(file=str, exist=file_exists)
	
	if (file_exists) then
		n = num_rows(str)
		allocate(var(n, columns))
		open(newunit=fileunit, file=str, position= 'rewind', action='read')
		do i = 1, size(var, 1)
			read(unit=fileunit, fmt=*, err=100) var(i, :) 
		end do
	else if (file_exists .eqv. .false.) then
		stop
	end if

	100 close(fileunit)
	
end subroutine load_to_var_2d_kind8_int

subroutine load_to_var_2d_kind16_int(filename, var, columns)

	character(len=*), intent(in) 								:: filename
	integer(int16), dimension(:,:), allocatable, intent(inout)  :: var
	integer, intent(in)											:: columns
	integer 													:: i, n, &
																   fileunit
	character(len=200)											:: str
	logical 													:: file_exists
	
	str = 'data/' // filename
	str = trim(str)
	inquire(file=str, exist=file_exists)
	
	if (file_exists) then
		n = num_rows(str)
		allocate(var(n, columns))
		open(newunit=fileunit, file=str, position= 'rewind', action='read')
		do i = 1, size(var, 1)
			read(unit=fileunit, fmt=*, err=100) var(i, :) 
		end do
	else if (file_exists .eqv. .false.) then
		stop
	end if

	100 close(fileunit)
	
end subroutine load_to_var_2d_kind16_int

subroutine load_to_var_2d_kind32_int(filename, var, columns)

	character(len=*), intent(in) 								:: filename
	integer(int32), dimension(:,:), allocatable, intent(inout)  :: var
	integer, intent(in)											:: columns
	integer 													:: i, n, &
																   fileunit
	character(len=200)											:: str
	logical 													:: file_exists
	
	str = 'data/' // filename
	str = trim(str)
	inquire(file=str, exist=file_exists)
	
	if (file_exists) then
		n = num_rows(str)
		allocate(var(n, columns))
		open(newunit=fileunit, file=str, position= 'rewind', action='read')
		do i = 1, size(var, 1)
			read(unit=fileunit, fmt=*, err=100) var(i, :) 
		end do
	else if (file_exists .eqv. .false.) then
		stop
	end if

	100 close(fileunit)
	
end subroutine load_to_var_2d_kind32_int


subroutine load_to_var_2d_kind32_real(filename, var, columns)

	character(len=*), intent(in) 								:: filename
	real(real32), dimension(:,:), allocatable, intent(inout)  	:: var
	integer, intent(in)											:: columns
	integer 													:: i, n, &
																   fileunit
	character(len=200)											:: str
	logical 													:: file_exists
	
	str = 'data/' // filename
	str = trim(str)
	inquire(file=str, exist=file_exists)
	
	if (file_exists) then
		n = num_rows(str)
		allocate(var(n, columns))
		open(newunit=fileunit, file=str, position= 'rewind', action='read')
		do i = 1, size(var, 1)
			read(unit=fileunit, fmt=*, err=100) var(i, :) 
		end do
	else if (file_exists .eqv. .false.) then
		stop
	end if

	100 close(fileunit)
	
end subroutine load_to_var_2d_kind32_real

subroutine load_to_var_2d_kind64_real(filename, var, columns)

	character(len=*), intent(in) 								:: filename
	real(real64), dimension(:,:), allocatable, intent(inout)  	:: var
	integer, intent(in)											:: columns
	integer 													:: i, n, &
																   fileunit
	character(len=200)											:: str
	logical 													:: file_exists
	
	str = 'data/' // filename
	str = trim(str)
	inquire(file=str, exist=file_exists)
	
	if (file_exists) then
		n = num_rows(str)
		allocate(var(n, columns))
		open(newunit=fileunit, file=str, position= 'rewind', action='read')
		do i = 1, size(var, 1)
			read(unit=fileunit, fmt=*, err=100) var(i, :) 
		end do
	else if (file_exists .eqv. .false.) then
		stop
	end if

	100 close(fileunit)
	
end subroutine load_to_var_2d_kind64_real

!write to file 1D

subroutine write_to_file_1d_logical(var, data_slot)
        
    logical, dimension(:), intent(in) 		:: var
    integer, intent(in) 					:: data_slot
    character(len=20) 						:: filename
    character(len=12) 						:: str
    integer 								:: i, fileunit

    write (str, *) data_slot
    str 		= adjustl(str)
    filename 	= 'data/datafile' // trim(str) // '.txt'
    filename 	= trim(filename)

    open(newunit=fileunit, file=filename, position='rewind', action='write')

    do i = 1, size(var, 1)
        write(unit=fileunit, fmt=*) var(i)
    end do

    close(unit=fileunit)

end subroutine write_to_file_1d_logical

subroutine write_to_file_1d_kind8_int(var, data_slot)
        
    integer(int8), dimension(:), intent(in) :: var
    integer, intent(in) 					:: data_slot
    character(len=20) 						:: filename
    character(len=12) 						:: str
    integer 								:: i, fileunit

    write (str, *) data_slot
    str 		= adjustl(str)
    filename 	= 'data/datafile' // trim(str) // '.txt'
    filename 	= trim(filename)

    open(newunit=fileunit, file=filename, position='rewind', action='write')

    do i = 1, size(var, 1)
        write(unit=fileunit, fmt=*) var(i)
    end do

    close(unit=fileunit)

end subroutine write_to_file_1d_kind8_int

subroutine write_to_file_1d_kind16_int(var, data_slot)
        
    integer(int16), dimension(:), intent(in) 	:: var
    integer, intent(in) 					 	:: data_slot
    character(len=20) 							:: filename
    character(len=12) 							:: str
    integer 									:: i, fileunit

    write (str, *) data_slot
    str 		= adjustl(str)
    filename 	= 'data/datafile' // trim(str) // '.txt'
    filename 	= trim(filename)

    open(newunit=fileunit, file=filename, position='rewind', action='write')

    do i = 1, size(var, 1)
        write(unit=fileunit, fmt=*) var(i)
    end do

    close(unit=fileunit)

end subroutine write_to_file_1d_kind16_int

subroutine write_to_file_1d_kind32_int(var, data_slot)
        
    integer(int32), dimension(:), intent(in) 	:: var
    integer, intent(in) 					 	:: data_slot
    character(len=20) 							:: filename
    character(len=12) 							:: str
    integer 									:: i, fileunit

    write (str, *) data_slot
    str 		= adjustl(str)
    filename 	= 'data/datafile' // trim(str) // '.txt'
    filename 	= trim(filename)

    open(newunit=fileunit, file=filename, position='rewind', action='write')

    do i = 1, size(var, 1)
        write(unit=fileunit, fmt=*) var(i)
    end do

    close(unit=fileunit)

end subroutine write_to_file_1d_kind32_int

subroutine write_to_file_1d_kind32_real(var, data_slot)
    
    real(real32), dimension(:), intent(in) 	:: var
    integer, intent(in) 					:: data_slot
    character(len=20) 						:: filename
    character(len=12) 						:: str
    integer 								:: i, fileunit

    write (str, *) data_slot
    str 		= adjustl(str)
    filename 	= 'data/datafile' // trim(str) // '.txt'
    filename 	= trim(filename)

    open(newunit=fileunit, file=filename, position='rewind', action='write')

    do i = 1, size(var, 1)
        write(unit=fileunit, fmt=*) var(i)
    end do

    close(unit=fileunit)

end subroutine write_to_file_1d_kind32_real

subroutine write_to_file_1d_kind64_real(var, data_slot)
    
    real(real64), dimension(:), intent(in) 	:: var
    integer, intent(in) 					:: data_slot
    character(len=20) 						:: filename
    character(len=12) 						:: str
    integer 								:: i, fileunit

    write (str, *) data_slot
    str 		= adjustl(str)
    filename 	= 'data/datafile' // trim(str) // '.txt'
    filename 	= trim(filename)

    open(newunit=fileunit, file=filename, position='rewind', action='write')

    do i = 1, size(var, 1)
        write(unit=fileunit, fmt=*) var(i)
    end do

    close(unit=fileunit)

end subroutine write_to_file_1d_kind64_real

!write to file 2D

subroutine write_to_file_2d_kind8_int(var, data_slot)
        
    integer(int8), dimension(:, :), intent(in) 	:: var
    integer, intent(in) 						:: data_slot
    character(len=20) 							:: filename
    character(len=12) 							:: str
    integer 									:: i, fileunit

    write (str, *) data_slot
    str 		= adjustl(str)
    filename 	= 'data/datafile' // trim(str) // '.txt'
    filename 	= trim(filename)

    open(newunit=fileunit, file=filename, position='rewind', action='write')

    do i = 1, size(var, 1)
        write(unit=fileunit, fmt=*) var(i, :)
    end do

    close(unit=fileunit)

end subroutine write_to_file_2d_kind8_int

subroutine write_to_file_2d_kind16_int(var, data_slot)
        
    integer(int16), dimension(:, :), intent(in) :: var
    integer, intent(in) 						:: data_slot
    character(len=20) 							:: filename
    character(len=12) 							:: str
    integer 									:: i, fileunit

    write (str, *) data_slot
    str 		= adjustl(str)
    filename	= 'data/datafile' // trim(str) // '.txt'
    filename 	= trim(filename)

    open(newunit=fileunit, file=filename, position='rewind', action='write')

    do i = 1, size(var, 1)
        write(unit=fileunit, fmt=*) var(i, :)
    end do

    close(unit=fileunit)

end subroutine write_to_file_2d_kind16_int

subroutine write_to_file_2d_kind32_int(var, data_slot)
        
    integer(int32), dimension(:, :), intent(in) :: var
    integer, intent(in) 						:: data_slot
    character(len=20) 							:: filename
    character(len=12) 							:: str
    integer 									:: i, fileunit

    write (str, *) data_slot
    str 		= adjustl(str)
    filename	= 'data/datafile' // trim(str) // '.txt'
    filename 	= trim(filename)

    open(newunit=fileunit, file=filename, position='rewind', action='write')

    do i = 1, size(var, 1)
        write(unit=fileunit, fmt=*) var(i, :)
    end do

    close(unit=fileunit)

end subroutine write_to_file_2d_kind32_int


subroutine write_to_file_2d_kind32_real(var, data_slot)
    
    real(real32), dimension(:, :), intent(in) 	:: var
    integer, intent(in) 						:: data_slot
    character(len=20) 							:: filename
    character(len=12) 							:: str
    integer 									:: i, fileunit

    write (str, *) data_slot
    str 		= adjustl(str)
    filename 	= 'data/datafile' // trim(str) // '.txt'
    filename 	= trim(filename)

    open(newunit=fileunit, file=filename, position='rewind', action='write')

    do i = 1, size(var, 1)
        write(unit=fileunit, fmt=*) var(i, :)
    end do

    close(unit=fileunit)

end subroutine write_to_file_2d_kind32_real

subroutine write_to_file_2d_kind64_real(var, data_slot)
    
    real(real64), dimension(:, :), intent(in) 	:: var
    integer, intent(in) 						:: data_slot
    character(len=20) 							:: filename
    character(len=12) 							:: str
    integer 									:: i, fileunit

    write (str, *) data_slot
    str 		= adjustl(str)
    filename 	= 'data/datafile' // trim(str) // '.txt'
    filename 	= trim(filename)

    open(newunit=fileunit, file=filename, position='rewind', action='write')

    do i = 1, size(var, 1)
        write(unit=fileunit, fmt=*) var(i, :)
    end do

    close(unit=fileunit)

end subroutine write_to_file_2d_kind64_real

end module mod_io
