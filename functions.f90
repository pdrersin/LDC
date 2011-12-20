module functions

  implicit none

  contains

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

!================================= det_3x3 ===================================80
!
! computes determinant of a 3x3 matrix
!
!=============================================================================80

    pure function det_3x3( mat )

      use set_precision,  only : dp

      real(dp), dimension(3,3), intent(in) :: mat
      real(dp)                             :: det_3x3

      continue

      det_3x3 = mat(1,1)*(mat(2,2)*mat(3,3)-mat(2,3)*mat(3,2)) &
              - mat(1,2)*(mat(2,1)*mat(3,3)-mat(2,3)*mat(3,1)) &
              + mat(1,3)*(mat(2,1)*mat(3,2)-mat(2,2)*mat(3,1))

    end function det_3x3

!=========================== find_available_unit =============================80
!
! finds an available I/O unit
!
!=============================================================================80

    integer function find_available_unit()

      integer, parameter :: min_unit = 7
      integer, parameter :: max_unit = 99

      logical :: exists, already_used

      continue

      do find_available_unit = min_unit, max_unit
        inquire(unit=find_available_unit, exist=exists, opened=already_used)
        if (exists .and. .not. already_used) return
      end do

    end function find_available_unit

end module functions
