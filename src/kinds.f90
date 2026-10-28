module kinds
  use, intrinsic :: iso_fortran_env, only: dp => real64, i8 => int64, stdout => output_unit, stderr => error_unit
  implicit none
  real(dp), parameter :: g = 9.81_dp
  real(dp), parameter :: pi = acos(-1.0_dp)
end module kinds
