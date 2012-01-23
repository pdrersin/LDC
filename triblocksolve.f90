!============================== triblocksolve ================================80
!
! Uses the Thomas algorithm to solve a generic block tridiagonal system
! Destroys the DD, UD, RHS matrices
!
!=============================================================================80
subroutine triblocksolve(neq, dof, LD, DD, UD, RHS, soln)

  use set_precision, only : dp

  implicit none

  integer,                          intent(in)    :: neq, dof
  real(dp), dimension(neq,neq,dof), intent(inout) :: LD, DD, UD
  real(dp), dimension(neq,dof),     intent(inout) :: RHS
  real(dp), dimension(neq,dof),     intent(out)   :: soln

  integer                      :: i
  real(dp), dimension(neq,neq) :: temp

  continue

! Normalize the first row...
!  call matrix_inv(neq, DD(:,:,1), temp)
  call mat_inv_3x3(DD(:,:,1), temp)

  DD(:,:,1) = matmul(temp,DD(:,:,1))
  UD(:,:,1) = matmul(temp,UD(:,:,1))
  RHS(:,1)  = matmul(temp,RHS(:,1))

! Loop to eliminate lower diagonal and then normalize the row
  do i = 2, dof
    DD(:,:,i) = DD(:,:,i) - matmul(LD(:,:,i), UD(:,:,i-1))
    RHS(:,i)  = RHS(:,i)  - matmul(LD(:,:,i), RHS(:,i-1))

!    call matrix_inv(neq, DD(:,:,i), temp)
    call mat_inv_3x3(DD(:,:,i), temp)

    DD(:,:,i) = matmul(temp,DD(:,:,i))
    UD(:,:,i) = matmul(temp,UD(:,:,i))
    RHS(:,i)  = matmul(temp,RHS(:,i))
  end do

! Back solve... since the diagonal is an identity matrix this is easy
  soln(:,dof) = RHS(:,dof)

  do i = dof-1,1,-1
    soln(:,i) = RHS(:,i) - matmul(UD(:,:,i), soln(:,i+1))
  end do

end subroutine triblocksolve
