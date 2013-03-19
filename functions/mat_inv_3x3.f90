!=============================== mat_inv_3x3 =================================80
!
! Calculates inverse of 3x3 matrix for speed!!!!
!
!=============================================================================80
pure function mat_inv_3x3( mat ) result( inv )

  use set_precision, only : dp

  real(dp), dimension(3,3), intent(in) :: mat
  real(dp), dimension(3,3)             :: inv

  continue

  inv(1,1) = mat(2,2)*mat(3,3)-mat(2,3)*mat(3,2)
  inv(2,1) = mat(2,3)*mat(3,1)-mat(2,1)*mat(3,3)
  inv(3,1) = mat(2,1)*mat(3,2)-mat(2,2)*mat(3,1)

  inv(1,2) = mat(1,3)*mat(3,2)-mat(1,2)*mat(3,3)
  inv(2,2) = mat(1,1)*mat(3,3)-mat(1,3)*mat(3,1)
  inv(3,2) = mat(1,2)*mat(3,1)-mat(1,1)*mat(3,2)

  inv(1,3) = mat(1,2)*mat(2,3)-mat(1,3)*mat(2,2)
  inv(2,3) = mat(1,3)*mat(2,1)-mat(1,1)*mat(2,3)
  inv(3,3) = mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)

  inv = inv/det_3x3(mat)

end function mat_inv_3x3
