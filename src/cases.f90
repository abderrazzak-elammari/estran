module cases
  use kinds
  use config
  implicit none
  private
  public :: bed, initial, has_exact, exact_h, wave_bc

  ! Thacker (1981), planar surface in a paraboloid, as in SWASHES
  real(dp), parameter :: th_a = 1.0_dp, th_h0 = 0.1_dp, th_eta = 0.5_dp, th_l = 4.0_dp

  ! Briggs et al. (1995) conical island, basin coordinates of the experiment
  real(dp), parameter :: xi = 12.96_dp, yi = 13.80_dp
  real(dp), parameter :: r_toe = 3.6_dp, r_top = 1.1_dp, z_top = 0.625_dp

contains

  elemental function bed(x, y) result(b)
    real(dp), intent(in) :: x, y
    real(dp) :: b, r
    select case (trim(case_name))
    case ('bump')
      b = 0.8_dp*exp(-5*(x - 0.9_dp)**2 - 50*(y - 0.5_dp)**2)
    case ('thacker')
      b = -th_h0*(1 - ((x - th_l/2)**2 + (y - th_l/2)**2)/th_a**2)
    case ('island')
      r = hypot(x - xi, y - yi)
      b = -depth + z_top*min(1.0_dp, max(0.0_dp, (r_toe - r)/(r_toe - r_top)))
    case default
      b = 0
    end select
  end function bed

  elemental subroutine initial(x, y, h, u, v)
    real(dp), intent(in) :: x, y
    real(dp), intent(out) :: h, u, v
    u = 0; v = 0
    select case (trim(case_name))
    case ('bump')
      h = 1 - bed(x, y)
      if (x > 0.05_dp .and. x < 0.15_dp) h = h + perturb
    case ('dambreak')
      h = merge(h_l, h_r, x < x_dam)
    case ('thacker')
      call thacker(x, y, 0.0_dp, h, u, v)
    case ('island')
      h = max(0.0_dp, -bed(x, y))
    case default
      h = 0
    end select
  end subroutine initial

  logical function has_exact()
    has_exact = any(trim(case_name) == ['dambreak', 'thacker '])
  end function has_exact

  function exact_h(x, y, t) result(h)
    real(dp), intent(in) :: x, y, t
    real(dp) :: h, u, v
    select case (trim(case_name))
    case ('dambreak')
      h = dambreak(x - x_dam, t)
    case ('thacker')
      call thacker(x, y, t, h, u, v)
    case default
      h = 0
    end select
  end function exact_h

  elemental subroutine thacker(x, y, t, h, u, v)
    real(dp), intent(in) :: x, y, t
    real(dp), intent(out) :: h, u, v
    real(dp) :: w
    w = sqrt(2*g*th_h0)/th_a
    h = th_eta*th_h0/th_a**2*(2*(x - th_l/2)*cos(w*t) + 2*(y - th_l/2)*sin(w*t) - th_eta) - bed(x, y)
    if (h > 0) then
      u = -th_eta*w*sin(w*t)
      v = th_eta*w*cos(w*t)
    else
      h = 0; u = 0; v = 0
    end if
  end subroutine thacker

  ! Stoker (wet bed) or Ritter (h_r = 0). The middle state comes from equating
  ! the rarefaction and shock velocities, solved by bisection on [h_r, h_l].
  function dambreak(xs, t) result(h)
    real(dp), intent(in) :: xs, t
    real(dp) :: h, cl, cm, um, hm, s, lo, hi
    integer :: it
    cl = sqrt(g*h_l)
    if (t <= 0) then
      h = merge(h_l, h_r, xs < 0)
      return
    end if
    if (h_r <= 0) then
      if (xs < -cl*t) then
        h = h_l
      else if (xs < 2*cl*t) then
        h = (2*cl - xs/t)**2/(9*g)
      else
        h = 0
      end if
      return
    end if
    lo = h_r; hi = h_l
    do it = 1, 200
      hm = 0.5_dp*(lo + hi)
      if (2*(cl - sqrt(g*hm)) - (hm - h_r)*sqrt(0.5_dp*g*(hm + h_r)/(hm*h_r)) > 0) then
        lo = hm
      else
        hi = hm
      end if
    end do
    cm = sqrt(g*hm)
    um = 2*(cl - cm)
    s = hm*um/(hm - h_r)
    if (xs < -cl*t) then
      h = h_l
    else if (xs < (um - cm)*t) then
      h = (2*cl - xs/t)**2/(9*g)
    else if (xs < s*t) then
      h = hm
    else
      h = h_r
    end if
  end function dambreak

  ! Solitary wave entering through x = 0. Its crest reaches the boundary at
  ! t0 = 4/(k c), when the leading tail is still below 0.2 % of the height.
  subroutine wave_bc(t, h, un)
    real(dp), intent(in) :: t
    real(dp), intent(out) :: h, un
    real(dp) :: k, c, eta
    k = sqrt(3*amp/(4*depth**3))
    c = sqrt(g*(depth + amp))
    eta = amp/cosh(k*c*(t - 4/(k*c)))**2
    h = depth + eta
    un = -c*eta/(depth + eta)     ! outward normal is -x
  end subroutine wave_bc

end module cases
