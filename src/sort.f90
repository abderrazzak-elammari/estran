module sort
  use kinds
  implicit none
  private
  public :: argsort

contains

  ! Stable bottom-up merge sort returning the permutation. Stability matters:
  ! every rank sorts the same arrays and must get exactly the same order, ties
  ! included, or the partitions and edge numbering would differ between ranks.
  subroutine argsort(key, idx)
    real(dp), intent(in) :: key(:)
    integer, intent(out) :: idx(:)
    integer, allocatable :: tmp(:)
    integer :: n, w, lo, mid, hi, i, j, k

    n = size(key)
    idx = [(i, i = 1, n)]
    allocate(tmp(n))
    w = 1
    do while (w < n)
      do lo = 1, n, 2*w
        mid = min(lo + w - 1, n)
        hi = min(lo + 2*w - 1, n)
        if (mid >= hi) cycle
        i = lo; j = mid + 1; k = lo
        do while (i <= mid .and. j <= hi)
          if (key(idx(j)) < key(idx(i))) then
            tmp(k) = idx(j); j = j + 1
          else
            tmp(k) = idx(i); i = i + 1
          end if
          k = k + 1
        end do
        tmp(k:k+mid-i) = idx(i:mid)
        k = k + mid - i + 1
        tmp(k:k+hi-j) = idx(j:hi)
        idx(lo:hi) = tmp(lo:hi)
      end do
      w = 2*w
    end do
  end subroutine argsort

end module sort
