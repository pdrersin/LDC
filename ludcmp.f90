!================================= ludcmp ====================================80
!
! This subroutine performs LU decomposition without pivoting,
! most ( CFD related ) matrices passed to it will be diagonally dominant
!
!=============================================================================80

subroutine ludcmp(neq, matrix, lower, upper)

  use set_precision, only : dp

  implicit none

  integer,                      intent(in)  :: neq
  real(dp), dimension(neq,neq), intent(in)  :: matrix
  real(dp), dimension(neq,neq), intent(out) :: upper, lower

  integer  :: i, j, k
  real(dp) :: factor

  real(dp), parameter :: zero = 0.0_dp
  real(dp), parameter :: one  = 1.0_dp

  continue

  lower = zero
  upper = zero

! set up diagonal of lower matrix
  do i = 1, neq
    lower(i,i) = one
  end do

! set first row of upper and first column of lower matrices
  upper(1,:) = matrix(1,:)
  lower(:,1) = matrix(:,1)/upper(1,1)

  do i = 2, neq-1
! set diagonal of upper matrix
    factor = zero
    do k = 1, i-1
      factor = factor + lower(i,k)*upper(k,i)
    end do
    upper(i,i) = matrix(i,i) - factor
!  end do

!  do i = 2, neq-1
! off diagonals of upper matrix
    do j = i+1, neq
      factor = zero
      do k = 1, i-1
        factor = factor + lower(i,k)*upper(k,j)
      end do
      upper(i,j) = matrix(i,j) - factor

! off diagonals of lower matrix
      factor = zero
      do k = 1, i-1
        factor = factor + lower(j,k)*upper(k,i)
      end do
      lower(j,i) = (matrix(j,i) - factor)/upper(i,i)
    end do

  end do

  factor = zero
  do k = 1, neq-1
    factor = factor + lower(neq,k)*upper(k,neq)
  end do

  upper(neq,neq) = matrix(neq,neq) - factor

end subroutine ludcmp

!================================= luvkslv ===================================80
!
! Performs the back substitution for the tridiagonal solver
!
!=============================================================================80

subroutine lubkslv(neq, lower, upper, b, x)

  use set_precision, only : dp

  implicit none

  integer,                      intent(in)  :: neq
  real(dp), dimension(neq,neq), intent(in)  :: lower, upper
  real(dp), dimension(neq),     intent(in)  :: b
  real(dp), dimension(neq),     intent(out) :: x

  integer  :: i, j
  real(dp) :: factor
  real(dp), dimension(neq) :: y

  real(dp), parameter :: zero = 0.0_dp

  continue

! forward substitution
  y(1) = b(1)
  do i = 2, neq
    factor = zero
    do j = 1, i-1
      factor = factor + lower(i,j)*y(j)
    end do
    y(i) = b(i) - factor
  end do

! back substitution
  x(neq) = y(neq)/upper(neq,neq)
  do i = neq-1,1,-1
    factor = zero
    do j = i+1, neq
      factor = factor + upper(i,j)*x(j)
    end do
    x(i) = (y(i) - factor)/upper(i,i)
  end do

end subroutine lubkslv

!=============================== matrix_inv ==================================80
!
! Uses ludcmp and lubkslv to find a matrix inverse
!
!=============================================================================80

subroutine matrix_inv(neq, matrix, inv)

  use set_precision, only : dp

  implicit none

  integer,                      intent(in)  :: neq
  real(dp), dimension(neq,neq), intent(in)  :: matrix
  real(dp), dimension(neq,neq), intent(out) :: inv

  integer :: i
  real(dp), dimension(neq) :: b, soln
  real(dp), dimension(neq,neq) :: l,u

  real(dp), parameter :: zero = 0.0_dp
  real(dp), parameter :: one  = 1.0_dp

  continue

! perform LU decomposition

  call ludcmp(neq, matrix, l, u)

  do i = 1, neq
    b    = zero
    b(i) = one
    call lubkslv(neq, l, u, b, soln)
    inv(:,i) = soln
  end do

end subroutine matrix_inv

!=============================== mat_inv_3x3 =================================80
!
! Calculates inverse of 3x3 matrix for speed!!!!
!
!=============================================================================80

subroutine mat_inv_3x3(mat, inv)

  use set_precision, only : dp
  use functions,     only : det_3x3

  implicit none

  real(dp), dimension(3,3), intent(in)  :: mat
  real(dp), dimension(3,3), intent(out) :: inv

  continue

  inv(1,1) = mat(2,2)*mat(3,3)-mat(2,3)*mat(3,2)
  inv(1,2) = mat(1,3)*mat(3,2)-mat(1,2)*mat(3,3)
  inv(1,3) = mat(1,2)*mat(2,3)-mat(1,3)*mat(2,2)

  inv(2,1) = mat(2,3)*mat(3,1)-mat(2,1)*mat(3,3)
  inv(2,2) = mat(1,1)*mat(3,3)-mat(1,3)*mat(3,1)
  inv(2,3) = mat(1,3)*mat(2,1)-mat(1,1)*mat(2,3)

  inv(3,1) = mat(2,1)*mat(3,2)-mat(2,2)*mat(3,1)
  inv(3,2) = mat(1,2)*mat(3,1)-mat(1,1)*mat(3,2)
  inv(3,3) = mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)

  inv = inv/det_3x3(mat)

end subroutine mat_inv_3x3



