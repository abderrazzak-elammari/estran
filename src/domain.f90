module domain
  use kinds
  use mpi_f08
  use mesh_mod
  use partition
  implicit none
  private
  public :: dom_t, build_domain, exchange, gsum, gmin, gmax

  ! What one rank sees: its own cells (1:nown), one layer of ghost cells
  ! (nown+1:nloc), the interior edges touching an own cell and the boundary
  ! edges of own cells. Local numbering follows the global numbering, and the
  ! edge lists follow the global edge order, so every cell accumulates its
  ! fluxes in the same order whatever the number of ranks. The solution is
  ! then bitwise independent of the partition.
  type dom_t
    integer :: rank = 0, nproc = 1
    integer :: nown = 0, nloc = 0
    integer, allocatable :: gid(:), g2l(:), part(:)
    real(dp), allocatable :: area(:), c(:,:)
    integer, allocatable :: nb(:,:)          ! (3, nown) face neighbours, 0 on the boundary
    real(dp), allocatable :: r(:,:,:)        ! (2, 3, nown) centroid to edge midpoints
    real(dp), allocatable :: w(:,:,:)        ! (2, 3, nown) least-squares weights

    integer :: ne = 0
    integer, allocatable :: e(:,:)           ! (2, ne) left and right local cell
    real(dp), allocatable :: en(:,:), elen(:), em(:,:)

    integer :: nbe = 0
    integer, allocatable :: bcell(:), btype(:)
    real(dp), allocatable :: bn(:,:), blen(:), bm(:,:)

    integer :: nnb = 0
    integer, allocatable :: nbr_rank(:), soff(:), roff(:), slist(:), rlist(:)

    integer :: nvtx = 0
    integer, allocatable :: vtx(:), conn(:,:)

    real(dp) :: t_comm = 0
  end type dom_t

  real(dp), allocatable, asynchronous :: sbuf(:), rbuf(:)
  type(MPI_Request), allocatable :: req(:)

contains

  subroutine build_domain(m, d)
    type(mesh_t), intent(in) :: m
    type(dom_t), intent(out) :: d
    real(dp), allocatable :: cx(:), cy(:)
    logical, allocatable :: ghost(:), used(:)
    integer, allocatable :: owner(:), cnt(:), nmap(:)
    real(dp) :: dd(2), a(2,2), det, xa(2), xb(2)
    integer :: t, k, i, j, l, q, n, ie, t1, k1, t2, nv

    call MPI_Comm_rank(MPI_COMM_WORLD, d%rank)
    call MPI_Comm_size(MPI_COMM_WORLD, d%nproc)

    allocate(cx(m%nt), cy(m%nt), d%part(m%nt))
    cx = (m%xy(1, m%tri(1,:)) + m%xy(1, m%tri(2,:)) + m%xy(1, m%tri(3,:)))/3
    cy = (m%xy(2, m%tri(1,:)) + m%xy(2, m%tri(2,:)) + m%xy(2, m%tri(3,:)))/3
    if (d%nproc > 1) then
      call rcb_partition(cx, cy, d%nproc, d%part)
    else
      d%part = 0
    end if

    allocate(ghost(m%nt), source=.false.)
    do t = 1, m%nt
      if (d%part(t) /= d%rank) cycle
      do k = 1, 3
        if (m%nbr(k, t) > 0) then
          if (d%part(m%nbr(k, t)) /= d%rank) ghost(m%nbr(k, t)) = .true.
        end if
      end do
    end do
    d%nown = count(d%part == d%rank)
    d%nloc = d%nown + count(ghost)
    allocate(d%gid(d%nloc), d%g2l(m%nt))
    d%g2l = 0
    i = 0
    do t = 1, m%nt
      if (d%part(t) == d%rank) then
        i = i + 1; d%gid(i) = t; d%g2l(t) = i
      end if
    end do
    do t = 1, m%nt
      if (ghost(t)) then
        i = i + 1; d%gid(i) = t; d%g2l(t) = i
      end if
    end do

    allocate(d%area(d%nloc), d%c(2, d%nloc))
    do i = 1, d%nloc
      t = d%gid(i)
      xa = m%xy(:, m%tri(2, t)) - m%xy(:, m%tri(1, t))
      xb = m%xy(:, m%tri(3, t)) - m%xy(:, m%tri(1, t))
      d%area(i) = 0.5_dp*(xa(1)*xb(2) - xa(2)*xb(1))
      d%c(:, i) = [cx(t), cy(t)]
    end do

    allocate(d%nb(3, d%nown), d%r(2, 3, d%nown), d%w(2, 3, d%nown))
    d%w = 0
    do i = 1, d%nown
      t = d%gid(i)
      a = 0
      do k = 1, 3
        xa = m%xy(:, m%tri(k, t)); xb = m%xy(:, m%tri(mod(k, 3) + 1, t))
        d%r(:, k, i) = 0.5_dp*(xa + xb) - d%c(:, i)
        d%nb(k, i) = 0
        if (m%nbr(k, t) > 0) d%nb(k, i) = d%g2l(m%nbr(k, t))
        if (d%nb(k, i) > 0) then
          dd = d%c(:, d%nb(k, i)) - d%c(:, i)
          a(1,1) = a(1,1) + dd(1)**2
          a(1,2) = a(1,2) + dd(1)*dd(2)
          a(2,2) = a(2,2) + dd(2)**2
        end if
      end do
      a(2,1) = a(1,2)
      det = a(1,1)*a(2,2) - a(1,2)**2
      ! a corner cell with a single neighbour has no usable gradient: it stays
      ! first order
      if (det <= 1.0e-10_dp*(a(1,1) + a(2,2))**2) cycle
      do k = 1, 3
        if (d%nb(k, i) == 0) cycle
        dd = d%c(:, d%nb(k, i)) - d%c(:, i)
        d%w(1, k, i) = ( a(2,2)*dd(1) - a(1,2)*dd(2))/det
        d%w(2, k, i) = (-a(1,2)*dd(1) + a(1,1)*dd(2))/det
      end do
    end do

    ! edges, in global order
    n = 0; l = 0
    do ie = 1, m%ne
      t1 = m%edge(1, ie); t2 = m%edge(3, ie)
      if (t2 > 0) then
        if (d%part(t1) == d%rank .or. d%part(t2) == d%rank) n = n + 1
      else if (d%part(t1) == d%rank) then
        l = l + 1
      end if
    end do
    d%ne = n; d%nbe = l
    allocate(d%e(2, n), d%en(2, n), d%elen(n), d%em(2, n))
    allocate(d%bcell(l), d%btype(l), d%bn(2, l), d%blen(l), d%bm(2, l))
    n = 0; l = 0
    do ie = 1, m%ne
      t1 = m%edge(1, ie); k1 = m%edge(2, ie); t2 = m%edge(3, ie)
      xa = m%xy(:, m%tri(k1, t1)); xb = m%xy(:, m%tri(mod(k1, 3) + 1, t1))
      dd = xb - xa
      if (t2 > 0) then
        if (d%part(t1) /= d%rank .and. d%part(t2) /= d%rank) cycle
        n = n + 1
        d%e(:, n) = [d%g2l(t1), d%g2l(t2)]
        d%elen(n) = norm2(dd)
        d%en(:, n) = [dd(2), -dd(1)]/d%elen(n)
        d%em(:, n) = 0.5_dp*(xa + xb)
      else if (d%part(t1) == d%rank) then
        l = l + 1
        d%bcell(l) = d%g2l(t1)
        d%btype(l) = m%ebc(ie)
        d%blen(l) = norm2(dd)
        d%bn(:, l) = [dd(2), -dd(1)]/d%blen(l)
        d%bm(:, l) = 0.5_dp*(xa + xb)
      end if
    end do

    ! halo lists. Both sides list the shared cells by increasing global id,
    ! so no index has to be communicated.
    allocate(owner(d%nloc - d%nown))
    owner = d%part(d%gid(d%nown+1:d%nloc))
    allocate(cnt(0:d%nproc-1), source=0)
    do i = 1, size(owner)
      cnt(owner(i)) = cnt(owner(i)) + 1
    end do
    d%nnb = count(cnt > 0)
    allocate(d%nbr_rank(d%nnb), d%soff(d%nnb+1), d%roff(d%nnb+1))
    d%nbr_rank = pack([(q, q = 0, d%nproc-1)], cnt > 0)
    allocate(d%rlist(size(owner)))
    d%roff(1) = 0
    do j = 1, d%nnb
      n = d%roff(j)
      do i = 1, size(owner)
        if (owner(i) == d%nbr_rank(j)) then
          n = n + 1; d%rlist(n) = d%nown + i
        end if
      end do
      d%roff(j+1) = n
    end do
    allocate(d%slist(d%nown*min(3, max(d%nnb, 1))))
    d%soff(1) = 0
    do j = 1, d%nnb
      n = d%soff(j)
      do i = 1, d%nown
        if (any(d%nb(:, i) > d%nown)) then
          if (any(owner(pack(d%nb(:, i), d%nb(:, i) > d%nown) - d%nown) == d%nbr_rank(j))) then
            n = n + 1; d%slist(n) = i
          end if
        end if
      end do
      d%soff(j+1) = n
    end do
    nv = 8
    allocate(sbuf(nv*max(d%soff(d%nnb+1), 1)), rbuf(nv*max(d%roff(d%nnb+1), 1)), req(2*d%nnb))

    ! vertices of own cells, for the VTK pieces
    allocate(used(m%nn), source=.false.)
    do i = 1, d%nown
      used(m%tri(:, d%gid(i))) = .true.
    end do
    d%nvtx = count(used)
    d%vtx = pack([(k, k = 1, m%nn)], used)
    allocate(nmap(m%nn), source=0)
    nmap(d%vtx) = [(k, k = 0, d%nvtx - 1)]
    allocate(d%conn(3, d%nown))
    do i = 1, d%nown
      d%conn(:, i) = nmap(m%tri(:, d%gid(i)))
    end do
  end subroutine build_domain

  ! Refresh the ghost values of a(nv, nloc), nv <= 8.
  subroutine exchange(d, a)
    type(dom_t), intent(inout) :: d
    real(dp), intent(inout) :: a(:,:)
    integer :: nv, j, i, o, n
    real(dp) :: t0

    if (d%nnb == 0) return
    t0 = MPI_Wtime()
    nv = size(a, 1)
    do j = 1, d%nnb
      o = nv*d%roff(j); n = nv*(d%roff(j+1) - d%roff(j))
      call MPI_Irecv(rbuf(o+1:o+n), n, MPI_DOUBLE_PRECISION, d%nbr_rank(j), 0, MPI_COMM_WORLD, req(j))
    end do
    do j = 1, d%nnb
      do i = d%soff(j) + 1, d%soff(j+1)
        sbuf(nv*(i-1)+1:nv*i) = a(:, d%slist(i))
      end do
      o = nv*d%soff(j); n = nv*(d%soff(j+1) - d%soff(j))
      call MPI_Isend(sbuf(o+1:o+n), n, MPI_DOUBLE_PRECISION, d%nbr_rank(j), 0, MPI_COMM_WORLD, req(d%nnb+j))
    end do
    call MPI_Waitall(2*d%nnb, req, MPI_STATUSES_IGNORE)
    do i = 1, d%roff(d%nnb+1)
      a(:, d%rlist(i)) = rbuf(nv*(i-1)+1:nv*i)
    end do
    d%t_comm = d%t_comm + MPI_Wtime() - t0
  end subroutine exchange

  real(dp) function gsum(x)
    real(dp), intent(in) :: x
    call MPI_Allreduce(x, gsum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD)
  end function gsum

  real(dp) function gmin(x)
    real(dp), intent(in) :: x
    call MPI_Allreduce(x, gmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD)
  end function gmin

  real(dp) function gmax(x)
    real(dp), intent(in) :: x
    call MPI_Allreduce(x, gmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD)
  end function gmax

end module domain
