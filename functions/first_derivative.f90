!============================= first_derivative ==============================80
!
! calculates a 2nd order first derivative
!
!=============================================================================80
pure function first_derivative(dl, var_l, var_r)

  use set_precision, only : dp
  use set_constants, only : half

  real(dp), intent(in) :: dl, var_l, var_r
  real(dp)             :: first_derivative

  continue

  first_derivative = half*(var_r-var_l)/dl

end function first_derivative
