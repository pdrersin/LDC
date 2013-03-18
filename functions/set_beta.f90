!============================= set_beta ======================================80
!
! sets Chorin's artificial compressibility factor
!
!=============================================================================80
pure function set_beta(velocity, k, u_lid)

  use set_precision, only : dp

  real(dp), intent(in) :: velocity, k, u_lid
  real(dp)             :: set_beta

  continue

  set_beta = sqrt( max(velocity**2, k*u_lid**2) )

end function set_beta
