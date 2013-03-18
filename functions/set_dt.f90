!================================ set_dt =====================================80
!
! sets the local timestep
!
!=============================================================================80
pure function set_dt( dx, dtd, cfl, u_vel, v_vel, beta )

  use set_precision,  only : dp
  use set_constants,  only : one, two, four

  real(dp), intent(in) :: dx, dtd, cfl, u_vel, v_vel, beta
  real(dp)             :: dtc, q
  real(dp)             :: set_dt

  continue

  q = sqrt(u_vel**2 + v_vel**2)

  dtc = dx/(two*q+beta)

  set_dt = cfl*(one/((one/dtc) + (one/dtd)))

end function set_dt
