!============================== set_precision ================================80
!
! Provides IEEE 754 compliant kinds via the iso_c_binding module and,
! by extension, C interoperability.
! The r* naming convention means real and rd is real_default
! The c* naming convention means complex and cd is complex_default
! The i* naming convention means integer and id is integer_default
! * == 4 is single precision, * == 8 is double, and * == 16 is quad
!
! FIXME: use a #ifdef to decide between SINGLE, DOUBLE, and QUAD defaults
!
!=============================================================================80
module set_precision

  use iso_c_binding, only : c_float, c_double, c_long_double, c_float_complex, &
                            c_double_complex, c_long_double_complex,           &
                            c_int, c_long, c_long_long

  implicit none

  private

  public :: r4, r8, r16, rd, dp
  public :: c4, c8, c16, cd
  public :: i4, i8, i16, id

! Real
  integer, parameter :: r4  = c_float
  integer, parameter :: r8  = c_double
  integer, parameter :: r16 = c_long_double
  integer, parameter :: dp  = r8
  integer, parameter :: rd  = r8

! Complex, even though c4 should use 8 bytes,
! we keep the naming consistent with the reals
  integer, parameter :: c4  = c_float_complex
  integer, parameter :: c8  = c_double_complex
  integer, parameter :: c16 = c_long_double_complex
  integer, parameter :: cd  = c8

! Integer, even though i4 should use 2 bytes,
! we keep the naming consistent with the reals
  integer, parameter :: i4  = c_int
  integer, parameter :: i8  = c_long
  integer, parameter :: i16 = c_long_long
  integer, parameter :: id  = i8

! For reference, here are how non F2003 compliant compilers would
! prescribe the kinds :
!  integer, parameter :: sngl = selected_real_kind( 6, 37)
!  integer, parameter :: dbl  = selected_real_kind(15, 307)
!  integer, parameter :: quad = selected_real_kind(33, 4931)

end module set_precision
