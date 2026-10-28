module partition
  use kinds
  use sort
  implicit none
  private
  public :: rcb_partition

contains

  ! Recursive coordinate bisection on the cell centroids. Cuts are made by
  ! cell count, not by length, so refined regions end up in small parts and
  ! every rank gets the same amount of work. Works for any number of parts:
  ! np parts are split into np/2 and np - np/2 with proportional counts.
  subroutine rcb_partition(cx, cy, np, part)
    real(dp), intent(in) :: cx(:), cy(:)
    integer, intent(in) :: np
    integer, intent(out) :: part(:)
    integer, allocatable :: ids(:)
    integer :: i
    ids = [(i, i = 1, size(cx))]
    call split(ids, np, 0)

  contains

    recursive subroutine split(ids, np, p0)
      integer, intent(inout) :: ids(:)
      integer, intent(in) :: np, p0
      integer, allocatable :: o(:)
      integer :: n, nl, cut
      real(dp) :: wx, wy

      n = size(ids)
      if (np == 1) then
        part(ids) = p0
        return
      end if
      wx = maxval(cx(ids)) - minval(cx(ids))
      wy = maxval(cy(ids)) - minval(cy(ids))
      allocate(o(n))
      if (wx >= wy) then
        call argsort(cx(ids), o)
      else
        call argsort(cy(ids), o)
      end if
      ids = ids(o)
      nl = np/2
      cut = int(int(n, i8)*nl/np)
      call split(ids(1:cut), nl, p0)
      call split(ids(cut+1:n), np - nl, p0 + nl)
    end subroutine split

  end subroutine rcb_partition

end module partition
