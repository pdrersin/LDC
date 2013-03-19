module setup

  use set_precision, only : dp

  implicit none

  private

  public :: ldc_allocate, ldc_deallocate

  public :: dt, beta, soln, soln_new

  public :: resid_out

! grid variables
  public :: x_nodes
  public :: y_nodes
  public :: length
  public :: xmin
  public :: xmax
  public :: ymin
  public :: ymax
  public :: dx      !derived
  public :: dy      !derived

! physical variables
  public :: re
  public :: rho
  public :: u_lid
  public :: p_guage
  public :: mu
  public :: nu      !derived

! solver variables
  public :: solver
  public :: max_iter
  public :: cfl
  public :: k
  public :: c2
  public :: conv_toler
  public :: visc_eps

  public :: dtd     !derived

  integer  :: x_nodes, y_nodes
  real(dp) :: xmin, xmax, ymin, ymax, dx, dy, length

  real(dp) :: re, rho, u_lid, p_guage, mu, nu, dtd

  character(len=8) :: solver

  integer  :: max_iter, resid_out
  real(dp) :: cfl, k, c2, conv_toler, visc_eps

! arrays

  real(dp), allocatable, dimension(:,:)   :: dt, beta
  real(dp), allocatable, dimension(:,:,:) :: soln, soln_new

contains

!=============================== ldc_allocate ================================80
!
! Allocates and initializes common arrays and variables
!
!=============================================================================80
  subroutine ldc_allocate

    use set_precision,  only : dp
    use set_constants,  only : zero, four, big

    continue

    allocate(soln(3,x_nodes,y_nodes), soln_new(3,x_nodes,y_nodes))
    allocate(dt(x_nodes,y_nodes), beta(x_nodes,y_nodes))

! Set BC's and initialize matrices

    soln(1,:,:) = p_guage
    soln(2,:,:) = zero
    soln(3,:,:) = zero

    soln(2,:,y_nodes) = u_lid

    soln_new = soln
    dt       = big
    beta     = big

! Set constant variables

    mu = rho*u_lid*length/re  !Dynamic viscosity, Pa*s
    nu = mu/rho               !Kinematic viscosity, m^2/s
    dx = (xmax-xmin)/real(x_nodes-1,dp)
    dy = (ymax-ymin)/real(y_nodes-1,dp)
    dtd = dx**2/(four*nu)

  end subroutine ldc_allocate

!============================== ldc_deallocate ===============================80
!
! Deallocates common arrays
!
!=============================================================================80
  subroutine ldc_deallocate

    deallocate(dt, beta, soln, soln_new)

  end subroutine ldc_deallocate

end module setup
