module mesh_mod
  use kinds
  use sort
  implicit none
  private
  public :: mesh_t, read_msh, BC_WALL, BC_OPEN, BC_WAVE, bc_name

  integer, parameter :: BC_WALL = 1, BC_OPEN = 2, BC_WAVE = 3

  ! Global triangulation, identical on every rank. Edge k of triangle t joins
  ! tri(k,t) and tri(mod(k,3)+1,t); triangles are counter-clockwise, so the
  ! edge normal (dy, -dx) points out of t.
  type mesh_t
    integer :: nn = 0, nt = 0, ne = 0
    real(dp), allocatable :: xy(:,:)    ! (2, nn)
    integer, allocatable :: tri(:,:)    ! (3, nt)
    integer, allocatable :: nbr(:,:)    ! (3, nt), 0 across a boundary edge
    integer, allocatable :: edge(:,:)   ! (3, ne): t1, k in t1, t2 (0 on the boundary)
    integer, allocatable :: ebc(:)      ! (ne) boundary type, 0 inside
  end type mesh_t

contains

  pure function bc_name(i) result(s)
    integer, intent(in) :: i
    character(len=4) :: s
    select case (i)
    case (BC_WALL); s = 'wall'
    case (BC_OPEN); s = 'open'
    case (BC_WAVE); s = 'wave'
    case default;   s = '?'
    end select
  end function bc_name

  ! Gmsh 2.2 ASCII (gmsh -2 -format msh22). Boundary lines must belong to a
  ! physical curve named wall, open or wave.
  subroutine read_msh(fname, m)
    character(len=*), intent(in) :: fname
    type(mesh_t), intent(out) :: m
    character(len=512) :: line
    character(len=64) :: pname
    integer :: u, ios, n, i, id, typ, ntag, nv, dim, tag, nl, ntri
    integer :: ptag(64), pbc(64), np
    integer :: buf(64)
    real(dp) :: ver, z
    integer, allocatable :: lnodes(:,:), lbc(:), tnodes(:,:)

    open(newunit=u, file=fname, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'cannot open mesh file '//trim(fname)

    np = 0
    do
      read(u, '(a)', iostat=ios) line
      if (ios /= 0) exit
      select case (trim(line))
      case ('$MeshFormat')
        read(u, *) ver, i
        if (ver < 2.0_dp .or. ver >= 3.0_dp .or. i /= 0) &
          error stop 'mesh must be Gmsh 2.2 ASCII: mesh with gmsh -2 -format msh22'
      case ('$PhysicalNames')
        read(u, *) n
        do i = 1, n
          read(u, *) dim, tag, pname
          if (dim /= 1) cycle
          np = np + 1
          ptag(np) = tag
          select case (trim(pname))
          case ('wall'); pbc(np) = BC_WALL
          case ('open'); pbc(np) = BC_OPEN
          case ('wave'); pbc(np) = BC_WAVE
          case default
            error stop 'unknown boundary name '//trim(pname)//' (expected wall, open or wave)'
          end select
        end do
      case ('$Nodes')
        read(u, *) m%nn
        allocate(m%xy(2, m%nn))
        do i = 1, m%nn
          read(u, *) id, m%xy(1, i), m%xy(2, i), z
          if (id /= i) error stop 'non-contiguous node numbering in mesh'
        end do
      case ('$Elements')
        read(u, *) n
        allocate(lnodes(2, n), lbc(n), tnodes(3, n))
        nl = 0; ntri = 0
        do i = 1, n
          read(u, '(a)') line
          read(line, *) id, typ, ntag
          select case (typ)
          case (1);  nv = 2
          case (2);  nv = 3
          case (15); cycle
          case default
            error stop 'only points, lines and 3-node triangles are supported'
          end select
          read(line, *) buf(1:3+ntag+nv)
          if (typ == 1) then
            nl = nl + 1
            lnodes(:, nl) = buf(4+ntag:5+ntag)
            lbc(nl) = lookup(buf(4))
          else
            ntri = ntri + 1
            tnodes(:, ntri) = buf(4+ntag:6+ntag)
          end if
        end do
      end select
    end do
    close(u)

    if (m%nn == 0 .or. .not. allocated(tnodes)) error stop 'no nodes or elements in '//trim(fname)
    m%nt = ntri
    m%tri = tnodes(:, 1:ntri)
    call orient(m)
    call renumber(m)
    call connect(m, lnodes(:, 1:nl), lbc(1:nl))

  contains

    integer function lookup(t)
      integer, intent(in) :: t
      integer :: j
      do j = 1, np
        if (ptag(j) == t) then
          lookup = pbc(j)
          return
        end if
      end do
      error stop 'boundary line without a physical name'
    end function lookup

  end subroutine read_msh

  subroutine orient(m)
    type(mesh_t), intent(inout) :: m
    integer :: t, a
    real(dp) :: d(2), e(2)
    do t = 1, m%nt
      d = m%xy(:, m%tri(2, t)) - m%xy(:, m%tri(1, t))
      e = m%xy(:, m%tri(3, t)) - m%xy(:, m%tri(1, t))
      if (d(1)*e(2) - d(2)*e(1) < 0) then
        a = m%tri(2, t); m%tri(2, t) = m%tri(3, t); m%tri(3, t) = a
      end if
    end do
  end subroutine orient

  ! Gmsh numbers triangles in the order the mesher created them, which jumps
  ! around in space. Sorting them along a Z-order (Morton) curve of their
  ! centroids puts neighbours close in memory; on the finest dam-break mesh
  ! (60 000 cells) the run time drops by a third.
  subroutine renumber(m)
    type(mesh_t), intent(inout) :: m
    real(dp), allocatable :: c(:,:), key(:)
    integer, allocatable :: p(:)
    integer(i8) :: ix, iy, k
    real(dp) :: lo(2), w
    integer :: t, bit

    allocate(c(2, m%nt), key(m%nt), p(m%nt))
    c = (m%xy(:, m%tri(1,:)) + m%xy(:, m%tri(2,:)) + m%xy(:, m%tri(3,:)))/3
    lo = minval(c, dim=2)
    w = maxval(maxval(c, dim=2) - lo)
    do t = 1, m%nt
      ix = int((c(1, t) - lo(1))/w*(2**20 - 1), i8)
      iy = int((c(2, t) - lo(2))/w*(2**20 - 1), i8)
      k = 0
      do bit = 0, 19
        k = ior(k, ishft(ibits(ix, bit, 1), 2*bit))
        k = ior(k, ishft(ibits(iy, bit, 1), 2*bit + 1))
      end do
      key(t) = real(k, dp)
    end do
    call argsort(key, p)
    m%tri = m%tri(:, p)
  end subroutine renumber

  ! Edges are matched by sorting the node pair keys a*(nn+1)+b. The keys are
  ! below 2**53 for any mesh we can hold in memory, so a real(dp) sort is exact.
  subroutine connect(m, lnodes, lbc)
    type(mesh_t), intent(inout) :: m
    integer, intent(in) :: lnodes(:,:), lbc(:)
    real(dp), allocatable :: key(:), lkey(:)
    integer, allocatable :: p(:), lp(:), e(:,:)
    integer :: t, k, i, j, a, b, t1, k1, t2, k2, lo, hi, mid
    real(dp) :: s

    s = real(m%nn + 1, dp)
    allocate(key(3*m%nt), p(3*m%nt))
    do t = 1, m%nt
      do k = 1, 3
        a = m%tri(k, t); b = m%tri(mod(k, 3) + 1, t)
        key(3*(t-1)+k) = real(min(a, b), dp)*s + max(a, b)
      end do
    end do
    call argsort(key, p)

    allocate(lkey(size(lbc)), lp(size(lbc)))
    do i = 1, size(lbc)
      lkey(i) = real(minval(lnodes(:, i)), dp)*s + maxval(lnodes(:, i))
    end do
    call argsort(lkey, lp)

    allocate(m%nbr(3, m%nt), e(3, 3*m%nt), m%ebc(3*m%nt))
    m%nbr = 0
    m%ne = 0
    i = 1
    do while (i <= 3*m%nt)
      t1 = (p(i) - 1)/3 + 1; k1 = mod(p(i) - 1, 3) + 1
      m%ne = m%ne + 1
      if (i < 3*m%nt .and. key(p(min(i+1, 3*m%nt))) == key(p(i))) then
        if (i + 2 <= 3*m%nt) then
          if (key(p(i+2)) == key(p(i))) error stop 'non-manifold edge in mesh'
        end if
        t2 = (p(i+1) - 1)/3 + 1; k2 = mod(p(i+1) - 1, 3) + 1
        m%nbr(k1, t1) = t2
        m%nbr(k2, t2) = t1
        e(:, m%ne) = [t1, k1, t2]
        m%ebc(m%ne) = 0
        i = i + 2
      else
        e(:, m%ne) = [t1, k1, 0]
        ! binary search for the boundary line carrying this edge
        lo = 1; hi = size(lkey)
        j = 0
        do while (lo <= hi)
          mid = (lo + hi)/2
          if (lkey(lp(mid)) < key(p(i))) then
            lo = mid + 1
          else if (lkey(lp(mid)) > key(p(i))) then
            hi = mid - 1
          else
            j = lp(mid); exit
          end if
        end do
        if (j == 0) error stop 'boundary edge not covered by a physical curve'
        m%ebc(m%ne) = lbc(j)
        i = i + 1
      end if
    end do
    ! edge order follows the triangles, so the flux loop walks memory forward
    deallocate(key, p)
    allocate(key(m%ne), p(m%ne))
    key = real(e(1, 1:m%ne), dp)*(m%nt + 1) + e(3, 1:m%ne)
    call argsort(key, p)
    m%edge = e(:, p)
    m%ebc = m%ebc(p)
  end subroutine connect

end module mesh_mod
