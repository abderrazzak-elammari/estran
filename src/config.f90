module config
  use kinds
  implicit none

  integer, parameter :: max_gauges = 32

  character(len=32)  :: case_name = 'bump'
  character(len=256) :: mesh_file = ''
  character(len=256) :: out_dir = 'out'
  character(len=64)  :: prefix = 'sol'

  real(dp) :: t_end = 1.0_dp
  real(dp) :: dt_out = 0.1_dp       ! 0 disables field output
  real(dp) :: dt_gauge = 0.0_dp     ! 0 disables gauges
  real(dp) :: cfl = 0.45_dp
  integer  :: order = 2
  real(dp) :: h_dry = 1.0e-3_dp     ! below this a cell and its neighbours fall back to first order
  real(dp) :: h_eps = 1.0e-6_dp     ! below this the velocity is set to zero
  real(dp) :: manning = 0.0_dp

  ! case parameters
  real(dp) :: perturb = 0.0_dp      ! bump: height of the initial pulse
  real(dp) :: h_l = 1.0_dp, h_r = 0.1_dp, x_dam = 25.0_dp
  real(dp) :: amp = 0.0144_dp       ! island: solitary wave height
  real(dp) :: depth = 0.32_dp

  integer  :: ng = 0
  real(dp) :: gx(max_gauges) = 0, gy(max_gauges) = 0

  namelist /run/ case_name, mesh_file, out_dir, prefix, t_end, dt_out, dt_gauge, cfl, order, &
                 h_dry, h_eps, manning, perturb, h_l, h_r, x_dam, amp, depth, ng, gx, gy

contains

  subroutine read_config(fname)
    character(len=*), intent(in) :: fname
    integer :: u, ios
    character(len=256) :: msg
    open(newunit=u, file=fname, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'cannot open '//trim(fname)
    read(u, nml=run, iostat=ios, iomsg=msg)
    close(u)
    if (ios /= 0) error stop 'bad namelist in '//trim(fname)//': '//trim(msg)
    if (ng > max_gauges) error stop 'too many gauges'
    if (order /= 1 .and. order /= 2) error stop 'order must be 1 or 2'
  end subroutine read_config

end module config
