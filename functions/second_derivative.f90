!========================== second_derivative ================================80
!
! calculates a 2nd order second derivative
!
!=============================================================================80
pure function second_derivative(dl, var_l, var_c, var_r)

  use set_precision, only : dp
  use set_constants, only : two

  real(dp), intent(in) :: dl, var_l, var_c, var_r
  real(dp)             :: second_derivative

  continue

  second_derivative = (var_r-two*var_c+var_l)/(dl*dl)

end function second_derivative
