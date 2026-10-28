module output
  use, intrinsic :: iso_fortran_env, only: int8, int32
  use kinds
  use mpi_f08
  use config
  use mesh_mod
  use domain
  use solver, only: q, pr, b, hmax
  use cases, only: has_exact, exact_h
  implicit none
  private
  public :: open_outputs, write_fields, sample_gauges, write_errors, close_outputs

  integer :: u_diag = 0, u_gauge = 0, u_err = 0
  integer :: nfiles = 0
  real(dp), allocatable :: tfile(:)
  integer, allocatable :: gcell(:)

contains

  recursive function itoa(i, w) result(s)
    integer, intent(in) :: i
    integer, intent(in), optional :: w
    character(len=:), allocatable :: s
    character(len=32) :: buf
    if (present(w)) then
      write(buf, '(i0.'//itoa(w)//')') i
    else
      write(buf, '(i0)') i
    end if
    s = trim(buf)
  end function itoa

  subroutine open_outputs(d, m)
    type(dom_t), intent(in) :: d
    type(mesh_t), intent(in) :: m
    integer :: k, i, t, owner
    real(dp) :: p(2), a(2), bb(2), c(2), l1, l2, l3

    if (d%rank == 0) call execute_command_line('mkdir -p '//trim(out_dir))
    call MPI_Barrier(MPI_COMM_WORLD)
    allocate(tfile(0))

    ! a gauge belongs to the lowest rank that owns a triangle containing it
    allocate(gcell(ng), source=0)
    do k = 1, ng
      p = [gx(k), gy(k)]
      do i = 1, d%nown
        t = d%gid(i)
        a = m%xy(:, m%tri(1, t)); bb = m%xy(:, m%tri(2, t)); c = m%xy(:, m%tri(3, t))
        l1 = cross(bb - p, c - p); l2 = cross(c - p, a - p); l3 = cross(a - p, bb - p)
        if (min(l1, l2, l3) >= -1.0e-12_dp) then
          gcell(k) = i
          exit
        end if
      end do
      owner = merge(d%rank, d%nproc, gcell(k) > 0)
      call MPI_Allreduce(MPI_IN_PLACE, owner, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD)
      if (owner == d%nproc) then
        if (d%rank == 0) write(stderr, '(a, i0, a)') 'gauge ', k, ' is outside the mesh'
        error stop
      end if
      if (owner /= d%rank) gcell(k) = 0
    end do

    if (d%rank /= 0) return
    open(newunit=u_diag, file=trim(out_dir)//'/diag.csv', status='replace')
    write(u_diag, '(a)') 'step,t,dt,volume,h_max'
    if (ng > 0) then
      open(newunit=u_gauge, file=trim(out_dir)//'/gauges.csv', status='replace')
      write(u_gauge, '(a)', advance='no') 't'
      do k = 1, ng
        write(u_gauge, '(a)', advance='no') ',g'//itoa(k)
      end do
      write(u_gauge, *)
    end if
    if (has_exact()) then
      open(newunit=u_err, file=trim(out_dir)//'/errors.csv', status='replace')
      write(u_err, '(a)') 't,l1,linf'
    end if

  contains

    pure real(dp) function cross(x, y)
      real(dp), intent(in) :: x(2), y(2)
      cross = x(1)*y(2) - x(2)*y(1)
    end function cross

  end subroutine open_outputs

  subroutine sample_gauges(d, t)
    type(dom_t), intent(in) :: d
    real(dp), intent(in) :: t
    real(dp) :: v(max(ng, 1)), t_dummy(max(ng, 1))
    integer :: k
    if (ng == 0) return
    v = 0
    do k = 1, ng
      if (gcell(k) > 0) v(k) = pr(4, gcell(k))
    end do
    if (d%rank == 0) then
      call MPI_Reduce(MPI_IN_PLACE, v, ng, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD)
      write(u_gauge, '(es16.8, *(:, ",", es16.8))') t, v(1:ng)
    else
      call MPI_Reduce(v, t_dummy, ng, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD)
    end if
  end subroutine sample_gauges

  subroutine write_errors(d, t)
    type(dom_t), intent(in) :: d
    real(dp), intent(in) :: t
    real(dp) :: e1, einf, area, e
    integer :: i
    if (.not. has_exact()) return
    e1 = 0; einf = 0
    do i = 1, d%nown
      e = abs(q(1, i) - exact_h(d%c(1, i), d%c(2, i), t))
      e1 = e1 + e*d%area(i)
      einf = max(einf, e)
    end do
    area = gsum(sum(d%area(1:d%nown)))
    e1 = gsum(e1)/area
    einf = gmax(einf)
    if (d%rank == 0) write(u_err, '(es16.8, 2(",", es16.8))') t, e1, einf
  end subroutine write_errors

  subroutine log_step(step, t, dt, vol, hm)
    integer, intent(in) :: step
    real(dp), intent(in) :: t, dt, vol, hm
    if (u_diag /= 0) write(u_diag, '(i0, 4(",", es20.12))') step, t, dt, vol, hm
  end subroutine log_step

  ! One raw-binary .vtu piece per rank, a .pvtu tying them together and a
  ! .pvd series so ParaView can play the run.
  subroutine write_fields(d, m, step, t, dt, vol)
    type(dom_t), intent(in) :: d
    type(mesh_t), intent(in) :: m
    integer, intent(in) :: step
    real(dp), intent(in) :: t, dt, vol
    character(len=:), allocatable :: base
    real(dp) :: hm
    integer :: ip

    hm = gmax(maxval(q(1, 1:d%nown)))
    if (d%rank == 0) call log_step(step, t, dt, vol, hm)
    if (dt_out <= 0) return

    base = trim(prefix)//'_'//itoa(nfiles, 5)
    call write_piece(trim(out_dir)//'/'//base//'_'//itoa(d%rank, 4)//'.vtu', d, m)
    if (d%rank /= 0) then
      nfiles = nfiles + 1
      return
    end if
    tfile = [tfile, t]
    nfiles = nfiles + 1
    block
      integer :: u
      open(newunit=u, file=trim(out_dir)//'/'//base//'.pvtu', status='replace')
      write(u, '(a)') '<?xml version="1.0"?>'
      write(u, '(a)') '<VTKFile type="PUnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">'
      write(u, '(a)') '  <PUnstructuredGrid GhostLevel="0">'
      write(u, '(a)') '    <PPoints><PDataArray type="Float64" NumberOfComponents="3"/></PPoints>'
      write(u, '(a)') '    <PCellData Scalars="h">'
      write(u, '(a)') '      <PDataArray type="Float64" Name="h"/>'
      write(u, '(a)') '      <PDataArray type="Float64" Name="b"/>'
      write(u, '(a)') '      <PDataArray type="Float64" Name="eta"/>'
      write(u, '(a)') '      <PDataArray type="Float64" Name="u"/>'
      write(u, '(a)') '      <PDataArray type="Float64" Name="v"/>'
      write(u, '(a)') '      <PDataArray type="Float64" Name="hmax"/>'
      write(u, '(a)') '      <PDataArray type="Int32" Name="rank"/>'
      write(u, '(a)') '    </PCellData>'
      do ip = 0, d%nproc - 1
        write(u, '(a)') '    <Piece Source="'//base//'_'//itoa(ip, 4)//'.vtu"/>'
      end do
      write(u, '(a)') '  </PUnstructuredGrid>'
      write(u, '(a)') '</VTKFile>'
      close(u)

      open(newunit=u, file=trim(out_dir)//'/'//trim(prefix)//'.pvd', status='replace')
      write(u, '(a)') '<?xml version="1.0"?>'
      write(u, '(a)') '<VTKFile type="Collection" version="0.1">'
      write(u, '(a)') '  <Collection>'
      do ip = 1, nfiles
        write(u, '(a, es16.8, a)') '    <DataSet timestep="', tfile(ip), '" file="'// &
          trim(prefix)//'_'//itoa(ip - 1, 5)//'.pvtu"/>'
      end do
      write(u, '(a)') '  </Collection>'
      write(u, '(a)') '</VTKFile>'
      close(u)
    end block
  end subroutine write_fields

  ! Inline base64 rather than raw appended binary: ParaView builds linked
  ! against expat >= 2.6 without the VTK fix (Ubuntu 24.04's ParaView 5.11,
  ! for one) cannot read raw appended data at all. Base64 costs a third more
  ! disk and works everywhere.
  subroutine write_piece(fname, d, m)
    character(len=*), intent(in) :: fname
    type(dom_t), intent(in) :: d
    type(mesh_t), intent(in) :: m
    character(len=1), parameter :: nl = achar(10)
    integer :: u, k, nc, nv
    real(dp), allocatable :: pts(:,:)

    nc = d%nown; nv = d%nvtx
    allocate(pts(3, nv))
    pts(1:2, :) = m%xy(:, d%vtx)
    pts(3, :) = 0

    open(newunit=u, file=fname, access='stream', form='unformatted', status='replace')
    write(u) '<?xml version="1.0"?>'//nl// &
      '<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">'//nl// &
      '  <UnstructuredGrid>'//nl// &
      '    <Piece NumberOfPoints="'//itoa(nv)//'" NumberOfCells="'//itoa(nc)//'">'//nl// &
      '      <Points>'//nl
    call put(u, 'Float64', '', 3, transfer(pts, [0_int8]))
    write(u) '      </Points>'//nl//'      <Cells>'//nl
    call put(u, 'Int32', 'connectivity', 1, transfer(int(d%conn, int32), [0_int8]))
    call put(u, 'Int32', 'offsets', 1, transfer([(int(3*k, int32), k = 1, nc)], [0_int8]))
    call put(u, 'UInt8', 'types', 1, [(5_int8, k = 1, nc)])
    write(u) '      </Cells>'//nl//'      <CellData Scalars="h">'//nl
    call put(u, 'Float64', 'h', 1, transfer(q(1, 1:nc), [0_int8]))
    call put(u, 'Float64', 'b', 1, transfer(b(1:nc), [0_int8]))
    call put(u, 'Float64', 'eta', 1, transfer(pr(4, 1:nc), [0_int8]))
    call put(u, 'Float64', 'u', 1, transfer(pr(2, 1:nc), [0_int8]))
    call put(u, 'Float64', 'v', 1, transfer(pr(3, 1:nc), [0_int8]))
    call put(u, 'Float64', 'hmax', 1, transfer(hmax, [0_int8]))
    call put(u, 'Int32', 'rank', 1, transfer([(int(d%rank, int32), k = 1, nc)], [0_int8]))
    write(u) '      </CellData>'//nl//'    </Piece>'//nl//'  </UnstructuredGrid>'//nl//'</VTKFile>'//nl
    close(u)

  contains

    ! the UInt64 byte count and the data form one base64 stream
    subroutine put(u, typ, name, ncomp, bytes)
      integer, intent(in) :: u, ncomp
      character(len=*), intent(in) :: typ, name
      integer(int8), intent(in) :: bytes(:)
      character(len=:), allocatable :: s
      s = '        <DataArray type="'//typ//'"'
      if (len(name) > 0) s = s//' Name="'//name//'"'
      if (ncomp > 1) s = s//' NumberOfComponents="'//itoa(ncomp)//'"'
      write(u) s//' format="binary">'//nl//'          '
      write(u) base64([transfer(int(size(bytes), i8), [0_int8]), bytes])
      write(u) nl//'        </DataArray>'//nl
    end subroutine put

  end subroutine write_piece

  pure function base64(x) result(s)
    integer(int8), intent(in) :: x(:)
    character(len=:), allocatable :: s
    character(len=64), parameter :: tab = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/'
    integer :: n, i, j, b0, b1, b2, k
    n = size(x)
    allocate(character(len=4*((n + 2)/3)) :: s)
    j = 0
    do i = 1, n, 3
      b0 = iand(int(x(i)), 255)
      b1 = 0; b2 = 0
      if (i + 1 <= n) b1 = iand(int(x(i+1)), 255)
      if (i + 2 <= n) b2 = iand(int(x(i+2)), 255)
      k = ishft(b0, -2)
      s(j+1:j+1) = tab(k+1:k+1)
      k = ior(ishft(iand(b0, 3), 4), ishft(b1, -4))
      s(j+2:j+2) = tab(k+1:k+1)
      k = ior(ishft(iand(b1, 15), 2), ishft(b2, -6))
      s(j+3:j+3) = tab(k+1:k+1)
      k = iand(b2, 63)
      s(j+4:j+4) = tab(k+1:k+1)
      if (i + 1 > n) s(j+3:j+3) = '='
      if (i + 2 > n) s(j+4:j+4) = '='
      j = j + 4
    end do
  end function base64

  subroutine close_outputs()
    if (u_diag /= 0) close(u_diag)
    if (u_gauge /= 0) close(u_gauge)
    if (u_err /= 0) close(u_err)
  end subroutine close_outputs

end module output
