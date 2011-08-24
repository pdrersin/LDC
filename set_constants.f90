module set_constants

  use set_precision, only : dp

  implicit none

  real(dp), parameter :: zero   = 0.0_dp
  real(dp), parameter :: sixth  = 1.0_dp/6.0_dp
  real(dp), parameter :: fifth  = 0.2_dp
  real(dp), parameter :: fourth = 0.25_dp
  real(dp), parameter :: third  = 1.0_dp/3.0_dp
  real(dp), parameter :: half   = 0.5_dp
  real(dp), parameter :: one    = 1.0_dp
  real(dp), parameter :: two    = 2.0_dp
  real(dp), parameter :: three  = 3.0_dp
  real(dp), parameter :: four   = 4.0_dp
  real(dp), parameter :: big    = huge(zero)
  real(dp), parameter :: small  = tiny(zero)

end module set_constants
