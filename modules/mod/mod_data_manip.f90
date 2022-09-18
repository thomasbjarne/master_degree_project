module mod_data_manip

    use iso_fortran_env, only: real32, real64
    implicit none

    interface high_to_low_diag
        module procedure :: high_to_low_diag_kind32, high_to_low_diag_kind64
    end interface high_to_low_diag

    interface low_to_high_diag
        module procedure :: low_to_high_diag_kind32, low_to_high_diag_kind64
    end interface low_to_high_diag

contains

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


end module mod_data_manip
