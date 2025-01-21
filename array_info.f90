module array_info
!  use allocation

  
  !    implicit none

!    save
     integer*4 :: ip,jp,kp

contains

 subroutine find_unique_sorted(temp_array, ik, affected, occurrences)
    implicit none

    ! Arguments
    integer, intent(in) :: temp_array(:)    ! Input array
    integer, intent(in) :: ik            ! Number of elements to consider
    integer, allocatable, intent(out) :: affected(:) ! Output unique sorted array
    integer, allocatable, intent(out) :: occurrences(:) ! Output occurrences array

    ! Local variables
    real, allocatable :: temp(:)
    integer :: n, i, count, current_count

    ! Check if ik is valid
    if (ik <= 0 .or. ik > size(temp_array)) then
      allocate(affected(0), occurrences(0))
      return
    end if

    n = ik
    allocate(temp(n))

    ! Step 1: Copy and sort the relevant portion of temp_array
    temp = temp_array(1:n)
    call quicksort_recursive(temp, 1, n)

    ! Step 2: Find unique elements and their occurrences
    count = 1
    current_count = 1
    do i = 2, n
      if (temp(i) /= temp(i - 1)) then
        count = count + 1
      end if
    end do

    allocate(affected(count), occurrences(count))

    affected(1) = temp(1)
    current_count = 1
    count = 1
    do i = 2, n
      if (temp(i) /= temp(i - 1)) then
        occurrences(count) = current_count
        count = count + 1
        affected(count) = temp(i)
        current_count = 1
      else
        current_count = current_count + 1
      end if
    end do

    occurrences(count) = current_count

    ! Deallocate temporary array
    deallocate(temp)

  end subroutine find_unique_sorted

  recursive subroutine quicksort_recursive(arr, left, right)
    implicit none

    ! Arguments
    real, intent(inout) :: arr(:)
    integer, intent(in) :: left, right

    ! Local variables
    integer :: i, j
    real :: pivot, temp

    if (left >= right) return

    pivot = arr(left)
    i = left
    j = right

    do while (i < j)
      do while (arr(j) >= pivot .and. i < j)
        j = j - 1
      end do

      if (i < j) then
        arr(i) = arr(j)
        i = i + 1
      end if

      do while (arr(i) <= pivot .and. i < j)
        i = i + 1
      end do

      if (i < j) then
        arr(j) = arr(i)
        j = j - 1
      end if
    end do

    arr(i) = pivot
    call quicksort_recursive(arr, left, i - 1)
    call quicksort_recursive(arr, i + 1, right)

  end subroutine quicksort_recursive





end module array_info
