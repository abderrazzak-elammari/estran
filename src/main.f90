program estran
  use mpi_f08
  use kinds
  use config
  use mesh_mod
  use domain
  use solver
  use output
  implicit none

  type(mesh_t) :: m
  type(dom_t) :: d
  character(len=256) :: arg
  real(dp) :: t, dt, t_out, t_gauge, t_stop, vol0, vol, wall0, wall
  integer :: step, ncells, u
  real(dp), parameter :: tiny_t = 1.0e-12_dp

  call MPI_Init()
  call get_command_argument(1, arg)
  if (len_trim(arg) == 0) then
    write(stderr, '(a)') 'usage: mpirun -np N estran case.nml'
    error stop
  end if
  call read_config(arg)

  wall0 = MPI_Wtime()
  call read_msh(mesh_file, m)
  call build_domain(m, d)
  call init_state(d)
  call open_outputs(d, m)
  ncells = m%nt

  if (d%rank == 0) write(stdout, '(a, a, i0, a, i0, a, i0, a)') trim(case_name), ': ', m%nt, ' triangles, ', &
    d%nproc, ' ranks, ', d%nnb, ' neighbours on rank 0'

  t = 0; step = 0
  vol0 = volume(d)
  call write_fields(d, m, step, t, 0.0_dp, vol0)
  call write_errors(d, t)
  call sample_gauges(d, t)
  t_out = merge(dt_out, huge(t), dt_out > 0)
  t_gauge = merge(dt_gauge, huge(t), dt_gauge > 0 .and. ng > 0)
  wall0 = MPI_Wtime()
  d%t_comm = 0

  do while (t < t_end - tiny_t)
    t_stop = min(t_end, t_out, t_gauge)
    ! outputs that add up to t_end within round-off must not leave a 1e-12 s step
    if (t_end - t_stop < 1.0e-9_dp*t_end) t_stop = t_end
    call advance(d, t, t_stop - t, dt)
    step = step + 1
    if (t_stop - (t + dt) < tiny_t) then
      t = t_stop
    else
      t = t + dt
    end if

    if (abs(t - t_gauge) < tiny_t) then
      call sample_gauges(d, t)
      t_gauge = t_gauge + dt_gauge
    end if
    if (abs(t - t_out) < tiny_t .or. abs(t - t_end) < tiny_t) then
      vol = volume(d)
      call write_fields(d, m, step, t, dt, vol)
      call write_errors(d, t)
      if (d%rank == 0) write(stdout, '(a, f10.4, a, i7, a, es10.3, a, es10.3)') &
        't =', t, '  step', step, '  dt', dt, '  dV/V0', (vol - vol0)/vol0
      if (abs(t - t_out) < tiny_t) t_out = t_out + dt_out
    end if
  end do

  wall = MPI_Wtime() - wall0
  if (d%rank == 0) then
    write(stdout, '(a, f9.2, a, f5.1, a, f7.2, a)') 'wall ', wall, ' s, communication ', &
      100*d%t_comm/max(wall, 1.0e-30_dp), ' % on rank 0, ', &
      real(ncells, dp)*step/wall/1.0e6_dp, ' M cell updates per second'
    open(newunit=u, file=trim(out_dir)//'/timing.csv', status='replace')
    write(u, '(a)') 'nproc,cells,steps,wall,comm_rank0'
    write(u, '(i0, ",", i0, ",", i0, ",", es14.6, ",", es14.6)') d%nproc, ncells, step, wall, d%t_comm
    close(u)
  end if
  call close_outputs()
  call MPI_Finalize()
end program estran
