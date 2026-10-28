module solver
  use kinds
  use config
  use mesh_mod, only: BC_WALL, BC_OPEN, BC_WAVE
  use domain
  use cases, only: bed, initial, wave_bc
  implicit none
  private
  public :: q, pr, b, hmax, init_state, advance, primitives, volume

  ! q = (h, hu, hv) conserved, pr = (h, u, v, eta) reconstructed, gp = limited
  ! gradients of pr, two components each.
  real(dp), allocatable :: q(:,:), pr(:,:), gp(:,:), b(:), hmax(:)
  real(dp), allocatable :: res(:,:), sa(:), q0(:,:)

contains

  subroutine init_state(d)
    type(dom_t), intent(in) :: d
    allocate(q(3, d%nloc), pr(4, d%nloc), gp(8, d%nloc), b(d%nloc), hmax(d%nown))
    allocate(res(3, d%nown), sa(d%nown), q0(3, d%nown))
    b = bed(d%c(1,:), d%c(2,:))
    call initial(d%c(1,:), d%c(2,:), q(1,:), q(2,:), q(3,:))
    q(2,:) = q(1,:)*q(2,:)
    q(3,:) = q(1,:)*q(3,:)
    gp = 0
    call primitives(d)
    hmax = q(1, 1:d%nown)
  end subroutine init_state

  subroutine primitives(d)
    type(dom_t), intent(in) :: d
    integer :: i
    do i = 1, d%nloc
      pr(1, i) = q(1, i)
      if (q(1, i) > h_eps) then
        pr(2:3, i) = q(2:3, i)/q(1, i)
      else
        pr(2:3, i) = 0
      end if
      pr(4, i) = q(1, i) + b(i)
    end do
  end subroutine primitives

  ! Least-squares gradients limited with Barth-Jespersen, so face values stay
  ! within the range of the neighbouring cell values. For h this is what keeps
  ! the reconstructed depth non-negative. Next to a dry cell the reconstruction
  ! of eta is meaningless (eta = b there), so those cells drop to first order.
  subroutine gradients(d)
    type(dom_t), intent(inout) :: d
    integer :: i, k, v, j
    real(dp) :: gr(2), wmin, wmax, dw, alpha

    gp(:, 1:d%nown) = 0
    if (order == 1) return
    do i = 1, d%nown
      if (pr(1, i) < h_dry) cycle
      if (any(d%nb(:, i) > 0 .and. pr(1, max(d%nb(:, i), 1)) < h_dry)) cycle
      do v = 1, 4
        gr = 0
        wmin = pr(v, i); wmax = pr(v, i)
        do k = 1, 3
          j = d%nb(k, i)
          if (j == 0) cycle
          gr = gr + d%w(:, k, i)*(pr(v, j) - pr(v, i))
          wmin = min(wmin, pr(v, j)); wmax = max(wmax, pr(v, j))
        end do
        alpha = 1
        do k = 1, 3
          dw = gr(1)*d%r(1, k, i) + gr(2)*d%r(2, k, i)
          if (dw > 1.0e-14_dp) then
            alpha = min(alpha, (wmax - pr(v, i))/dw)
          else if (dw < -1.0e-14_dp) then
            alpha = min(alpha, (wmin - pr(v, i))/dw)
          end if
        end do
        gp(2*v-1:2*v, i) = alpha*gr
      end do
    end do
    call exchange(d, gp)
  end subroutine gradients

  pure subroutine face(i, x, h, u, v, z)
    integer, intent(in) :: i
    real(dp), intent(in) :: x(2)
    real(dp), intent(out) :: h, u, v, z
    real(dp) :: eta
    h   = max(0.0_dp, pr(1, i) + gp(1, i)*x(1) + gp(2, i)*x(2))
    u   = pr(2, i) + gp(3, i)*x(1) + gp(4, i)*x(2)
    v   = pr(3, i) + gp(5, i)*x(1) + gp(6, i)*x(2)
    eta = pr(4, i) + gp(7, i)*x(1) + gp(8, i)*x(2)
    z = eta - h
  end subroutine face

  ! HLL in the frame of the edge, all three components. Wave speeds from
  ! Toro's two-rarefaction estimate, with the dry-bed limits. Upwinding the
  ! tangential velocity on the contact (HLLC) diffuses shear less, but on
  ! these triangulations it left spurious plane waves behind a passing pulse,
  ! of a few tenths of its height; plain HLL does not.
  pure subroutine hll(hl, ul, vl, hr, ur, vr, n, f, smax)
    real(dp), intent(in) :: hl, ul, vl, hr, ur, vr, n(2)
    real(dp), intent(out) :: f(3), smax
    real(dp) :: unl, utl, unr, utr, cl, cr, sl, sr, cs, us, fl(3), fr(3), fh(3)

    f = 0; smax = 0
    if (hl <= 0 .and. hr <= 0) return
    unl = ul*n(1) + vl*n(2); utl = -ul*n(2) + vl*n(1)
    unr = ur*n(1) + vr*n(2); utr = -ur*n(2) + vr*n(1)
    cl = sqrt(g*hl); cr = sqrt(g*hr)
    if (hl <= 0) then
      sl = unr - 2*cr; sr = unr + cr
    else if (hr <= 0) then
      sl = unl - cl; sr = unl + 2*cl
    else
      us = 0.5_dp*(unl + unr) + cl - cr
      cs = 0.5_dp*(cl + cr) + 0.25_dp*(unl - unr)
      sl = min(unl - cl, us - cs)
      sr = max(unr + cr, us + cs)
    end if
    fl = [hl*unl, hl*unl**2 + 0.5_dp*g*hl**2, hl*unl*utl]
    fr = [hr*unr, hr*unr**2 + 0.5_dp*g*hr**2, hr*unr*utr]
    if (sl >= 0) then
      fh = fl
    else if (sr <= 0) then
      fh = fr
    else
      fh = (sr*fl - sl*fr + sl*sr*([hr, hr*unr, hr*utr] - [hl, hl*unl, hl*utl]))/(sr - sl)
    end if
    f = [fh(1), fh(2)*n(1) - fh(3)*n(2), fh(2)*n(2) + fh(3)*n(1)]
    smax = max(abs(sl), abs(sr))
  end subroutine hll

  ! Hydrostatic reconstruction (Audusse et al. 2004). The bed source term is
  ! split into the interface pressure corrections, folded into the edge loop
  ! as g/2 h*^2 n, and a centred part -g h grad(eta) per cell. At rest, eta is
  ! flat, h* is the same on both sides and every term cancels exactly: the
  ! lake stays at rest to round-off.
  subroutine residual(d, t)
    type(dom_t), intent(inout) :: d
    real(dp), intent(in) :: t
    integer :: ie, l, r, i
    real(dp) :: hl, ul, vl, zl, hr, ur, vr, zr, zs, hls, hrs, f(3), s, un, n(2), len

    call exchange(d, q)
    call primitives(d)
    call gradients(d)

    res = 0; sa = 0
    do ie = 1, d%ne
      l = d%e(1, ie); r = d%e(2, ie)
      n = d%en(:, ie); len = d%elen(ie)
      call face(l, d%em(:, ie) - d%c(:, l), hl, ul, vl, zl)
      call face(r, d%em(:, ie) - d%c(:, r), hr, ur, vr, zr)
      zs = max(zl, zr)
      hls = max(0.0_dp, hl + zl - zs)
      hrs = max(0.0_dp, hr + zr - zs)
      call hll(hls, ul, vl, hrs, ur, vr, n, f, s)
      if (l <= d%nown) then
        res(:, l) = res(:, l) - len*f
        res(2:3, l) = res(2:3, l) + len*0.5_dp*g*hls**2*n
        sa(l) = sa(l) + len*s
      end if
      if (r <= d%nown) then
        res(:, r) = res(:, r) + len*f
        res(2:3, r) = res(2:3, r) - len*0.5_dp*g*hrs**2*n
        sa(r) = sa(r) + len*s
      end if
    end do

    do ie = 1, d%nbe
      i = d%bcell(ie)
      n = d%bn(:, ie); len = d%blen(ie)
      call face(i, d%bm(:, ie) - d%c(:, i), hl, ul, vl, zl)
      select case (d%btype(ie))
      case (BC_WALL)
        un = ul*n(1) + vl*n(2)
        hr = hl; ur = ul - 2*un*n(1); vr = vl - 2*un*n(2)
      case (BC_WAVE)
        call wave_bc(t, hr, un)
        ur = un*n(1); vr = un*n(2)
      case default
        hr = hl; ur = ul; vr = vl
      end select
      call hll(hl, ul, vl, hr, ur, vr, n, f, s)
      res(:, i) = res(:, i) - len*f
      res(2:3, i) = res(2:3, i) + len*0.5_dp*g*hl**2*n
      sa(i) = sa(i) + len*s
    end do

    do i = 1, d%nown
      res(2:3, i) = res(2:3, i) - g*pr(1, i)*gp(7:8, i)*d%area(i)
    end do
  end subroutine residual

  ! One SSP-RK2 step. dt comes from the first stage: area / sum(len * s) is
  ! the positivity limit of the first-order hydrostatic scheme, and cfl < 0.5
  ! keeps it for the MUSCL stages.
  subroutine advance(d, t, dt_cap, dt)
    type(dom_t), intent(inout) :: d
    real(dp), intent(in) :: t, dt_cap
    real(dp), intent(out) :: dt
    integer :: i
    real(dp) :: dtl

    q0 = q(:, 1:d%nown)
    call residual(d, t)
    dtl = huge(1.0_dp)
    do i = 1, d%nown
      if (sa(i) > 0) dtl = min(dtl, d%area(i)/sa(i))
    end do
    dt = min(cfl*gmin(dtl), dt_cap)

    do i = 1, d%nown
      q(:, i) = q0(:, i) + dt*res(:, i)/d%area(i)
    end do
    call fix(d)
    call residual(d, t + dt)
    do i = 1, d%nown
      q(:, i) = 0.5_dp*(q0(:, i) + q(:, i) + dt*res(:, i)/d%area(i))
    end do
    call fix(d)
    if (manning > 0) call friction(d, dt)
    hmax = max(hmax, q(1, 1:d%nown))
  end subroutine advance

  subroutine fix(d)
    type(dom_t), intent(in) :: d
    integer :: i
    do i = 1, d%nown
      if (q(1, i) < h_eps) then
        q(1, i) = max(q(1, i), 0.0_dp)
        q(2:3, i) = 0
      end if
    end do
  end subroutine fix

  ! Manning friction, implicit in the velocity so it cannot reverse the flow
  subroutine friction(d, dt)
    type(dom_t), intent(in) :: d
    real(dp), intent(in) :: dt
    integer :: i
    real(dp) :: h, vel
    do i = 1, d%nown
      h = q(1, i)
      if (h <= h_eps) cycle
      vel = hypot(q(2, i), q(3, i))/h
      q(2:3, i) = q(2:3, i)/(1 + dt*g*manning**2*vel/h**(4.0_dp/3))
    end do
  end subroutine friction

  real(dp) function volume(d)
    type(dom_t), intent(in) :: d
    volume = gsum(sum(q(1, 1:d%nown)*d%area(1:d%nown)))
  end function volume

end module solver
