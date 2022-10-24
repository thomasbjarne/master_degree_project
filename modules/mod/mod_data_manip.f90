module mod_data_manip

    use iso_fortran_env, only: real32, real64
    use mod_derived_types
    implicit none

    interface high_to_low_diag
        module procedure :: high_to_low_diag_kind32, high_to_low_diag_kind64
    end interface high_to_low_diag

    interface low_to_high_diag
        module procedure :: low_to_high_diag_kind32, low_to_high_diag_kind64
    end interface low_to_high_diag

    interface union
        module procedure :: union_real, &
            union_real2darray, &
            union_line, &
            union_triangle, &
            union_circle
    end interface union

    interface in
        module procedure :: in_line, &
            in_real, &
            in_triangle
    end interface in

    interface count_appearance
        module procedure :: count_appearance_line
    end interface count_appearance

contains

!SORTING ALGORITHMS
    pure function high_to_low_diag_kind32(A) result(A_sorted)

        real(real32), intent(in), dimension(:,:) :: A
        real(real32), dimension(:,:), allocatable :: A_sorted
        real(real32), dimension(:), allocatable :: diag_vec, diag_vec_sorted
        logical, dimension(:), allocatable :: mask_vec
        integer :: i, n

        n = size(A,1)
        allocate( A_sorted(n,n), diag_vec(n), diag_vec_sorted(n), mask_vec(n) )
        A_sorted=0

        do i = 1, n
            diag_vec(i) = A(i,i)
        end do

        mask_vec = .true.
        do i = 1, n
            diag_vec_sorted(i) = maxval(diag_vec, mask=mask_vec)
            mask_vec(maxloc(diag_vec, mask=mask_vec)) = .false.
        end do

        do i = 1, n
            A_sorted(i,i) = diag_vec_sorted(i)
        end do

		deallocate( diag_vec, diag_vec_sorted, mask_vec )

    end function high_to_low_diag_kind32

    pure function high_to_low_diag_kind64(A) result(A_sorted)

        real(real64), intent(in), dimension(:,:) :: A
        real(real64), dimension(:,:), allocatable :: A_sorted
        real(real64), dimension(:), allocatable :: diag_vec, diag_vec_sorted
        logical, dimension(:), allocatable :: mask_vec
        integer :: i, n

        n = size(A,1)
        allocate( A_sorted(n,n), diag_vec(n), diag_vec_sorted(n), mask_vec(n) )
        A_sorted=0

        do i = 1, n
            diag_vec(i) = A(i,i)
        end do

        mask_vec = .true.
        do i = 1, n
            diag_vec_sorted(i) = maxval(diag_vec, mask=mask_vec)
            mask_vec(maxloc(diag_vec, mask=mask_vec)) = .false.
        end do

        do i = 1, n
            A_sorted(i,i) = diag_vec_sorted(i)
        end do

		deallocate( diag_vec, diag_vec_sorted, mask_vec )

    end function high_to_low_diag_kind64

    pure function low_to_high_diag_kind32(A) result(A_sorted)

        real(real32), intent(in), dimension(:,:) :: A
        real(real32), dimension(:,:), allocatable :: A_sorted
        real(real32), dimension(:), allocatable :: diag_vec, diag_vec_sorted
        logical, dimension(:), allocatable :: mask_vec
        integer :: i, n

        n = size(A,1)
        allocate( A_sorted(n,n), diag_vec(n), diag_vec_sorted(n), mask_vec(n) )
        A_sorted=0

        do i = 1, n
            diag_vec(i) = A(i,i)
        end do

        mask_vec = .true.
        do i = 1, n
            diag_vec_sorted(i) = minval(diag_vec, mask=mask_vec)
            mask_vec(minloc(diag_vec, mask=mask_vec)) = .false.
        end do

        do i = 1, n
            A_sorted(i,i) = diag_vec_sorted(i)
        end do

		deallocate( diag_vec, diag_vec_sorted, mask_vec )

    end function low_to_high_diag_kind32

    pure function low_to_high_diag_kind64(A) result(A_sorted)

        real(real64), intent(in), dimension(:,:) :: A
        real(real64), dimension(:,:), allocatable :: A_sorted
        real(real64), dimension(:), allocatable :: diag_vec, diag_vec_sorted
        logical, dimension(:), allocatable :: mask_vec
        integer :: i, n

        n = size(A,1)
        allocate( A_sorted(n,n), diag_vec(n), diag_vec_sorted(n), mask_vec(n) )
        A_sorted=0

        do i = 1, n
            diag_vec(i) = A(i,i)
        end do

        mask_vec = .true.
        do i = 1, n
            diag_vec_sorted(i) = minval(diag_vec, mask=mask_vec)
            mask_vec(minloc(diag_vec, mask=mask_vec)) = .false.
        end do

        do i = 1, n
            A_sorted(i,i) = diag_vec_sorted(i)
        end do

		deallocate( diag_vec, diag_vec_sorted, mask_vec )
		
    end function low_to_high_diag_kind64

!ARRAY MANIPULATION
    subroutine union_real(element, set) 

        real, intent(in) :: element
        real, intent(in out), dimension(:),allocatable :: set
        real, dimension(:), allocatable :: temp

        allocate(temp(size(set) +1 ))

        temp(1:size(set)) = set
        temp(size(temp)) = element
        deallocate(set)
        allocate(set(size(temp)))
        set = temp

    end subroutine union_real

    subroutine union_real2darray(element, set) 

        real, intent(in), dimension(:,:) :: element
        real, intent(in out), dimension(:,:),allocatable :: set
        real, dimension(:,:), allocatable :: temp
        integer :: n, m

        n = size(set, 1)
        m = size(set, 2)

        allocate(temp(n + size(element,1), m + size(element, 2)))

        temp(1:n, 1:m) = set
        temp(n+1:n+size(element,1), m+1:m+size(element,2)) = element
        deallocate(set)
        allocate(set(size(temp,1), size(temp,2)))
        set = temp

    end subroutine union_real2darray

    subroutine union_line(element, set) 

        type(line), intent(in) :: element
        type(line), intent(in out), dimension(:),allocatable :: set
        type(line), dimension(:), allocatable :: temp

        allocate(temp(size(set) +1 ))

        temp(1:size(set)) = set
        temp(size(temp)) = element
        deallocate(set)
        allocate(set(size(temp)))
        set = temp

    end subroutine union_line

    subroutine union_triangle(element, set) 

        type(triangle), intent(in) :: element
        type(triangle), intent(in out), dimension(:),allocatable :: set
        type(triangle), dimension(:), allocatable :: temp

        allocate(temp(size(set) +1 ))

        temp(1:size(set)) = set
        temp(size(temp)) = element
        deallocate(set)
        allocate(set(size(temp)))
        set = temp

    end subroutine union_triangle

    subroutine union_circle(element, set) 

        type(circle), intent(in) :: element
        type(circle), intent(in out), dimension(:),allocatable :: set
        type(circle), dimension(:), allocatable :: temp

        allocate(temp(size(set) +1 ))

        temp(1:size(set)) = set
        temp(size(temp)) = element
        deallocate(set)
        allocate(set(size(temp)))
        set = temp

    end subroutine union_circle
    
! COMPARE ARRAY

    pure function in_real(element, set) result(bool)

        real, intent(in) :: element
        real, intent(in), dimension(:) :: set
        logical :: bool
        integer :: i, n
        
        bool = .false.
        n = size(set)
        do i = 1, n
            if (element == set(i)) then
                bool = .true.
                exit
            end if
        end do
        
    end function in_real

    pure function in_line(element, set) result(bool)

        type(line), intent(in) :: element
        type(line), intent(in), dimension(:) :: set
        logical :: bool
        integer :: i, n
        
        bool = .false.
        n = size(set)
        do i = 1, n
            if (element%p1(1,1) == set(i)%p1(1,1) &
            .and. element%p1(1,2) == set(i)%p1(1,2) &
            .and. element%p2(1,1) == set(i)%p2(1,1) &
            .and. element%p2(1,2) == set(i)%p2(1,2)) then
                bool = .true.
                exit
            end if
        end do
        
    end function in_line

    pure function in_triangle(element, set) result(bool)

        type(triangle), intent(in) :: element
        type(triangle), intent(in), dimension(:) :: set
        logical :: bool
        integer :: i, n

        bool = .false.
        n = size(set)
        do i = 1, n
            if(in_line(element%edges(1), set(i)%edges) &
            .and. in_line(element%edges(2), set(i)%edges) &
            .and. in_line(element%edges(3), set(i)%edges)) then
                bool = .true.
            end if
        end do

    end function in_triangle

    pure function point_in_circle(point, B) result(bool)

        real, intent(in), dimension(1,2) :: point
        type(circle), intent(in) :: B
        logical :: bool

        bool = .false.
        if (norm2(point - B%center) < B%radius) then
            bool = .true.
        end if

    end function point_in_circle

    pure function count_appearance_line(element, set) result(count)

        type(line), intent(in) :: element
        type(line), intent(in), dimension(:) :: set
        integer :: i, n, count

        count = 0
        do i = 1, n
            if (element%p1(1,1) == set(i)%p1(1,1) &
            .and. element%p1(1,2) == set(i)%p1(1,2) &
            .and. element%p2(1,1) == set(i)%p2(1,1) &
            .and. element%p2(1,2) == set(i)%p2(1,2)) then
                count = count + 1
                exit
            end if
        end do

    end function count_appearance_line

    pure function equals_triangle(t1, t2) result(bool)

        type(triangle), intent(in) :: t1, t2
        logical :: bool

        bool = .false.
        if(in_line(t1%edges(1), t2%edges) &
            .and. in_line(t1%edges(2), t2%edges) &
            .and. in_line(t1%edges(3), t2%edges)) then
                bool = .true.
        end if

    end function equals_triangle

end module mod_data_manip
