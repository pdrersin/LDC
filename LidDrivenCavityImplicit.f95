!  This code calculates the flow in a lid driven cavity
module constants

  use set_precision, only : dp

  implicit none

  real(dp), parameter :: zero   = 0.0_dp
  real(dp), parameter :: sixth  = 1.0_dp/6.0_dp
  real(dp), parameter :: fifth  = 0.2_dp
  real(dp), parameter :: fourth = 0.25_dp
  real(dp), parameter :: third  = 1.0_dp/3.0_dp
  real(dp), parameter :: half   = 0.5_dp
  real(dp), parameter :: one    = 1.0_dp
  real(dp), parameter :: two    = 2.0_dp
  real(dp), parameter :: three  = 3.0_dp
  real(dp), parameter :: four   = 4.0_dp
  real(dp), parameter :: big    = huge(zero)
  real(dp), parameter :: small  = tiny(zero)

end module constants

module functions

  implicit none

  contains

function first_derivative(dl, var_l, var_r)

  use set_precision, only : dp
  use constants,     only : half

  real(dp), intent(in) :: dl, var_l, var_r
  real(dp)             :: first_derivative

continue

  first_derivative = half*(var_r-var_l)/dl

end function first_derivative

function second_derivative(dl, var_l, var_c, var_r)

  use set_precision, only : dp
  use constants,     only : two

  real(dp), intent(in) :: dl, var_l, var_c, var_r
  real(dp)             :: second_derivative

continue

  second_derivative = (var_r-two*var_c+var_l)/dl**2

end function second_derivative
end module functions

program LidDrivenCavityImplicit

  use set_precision, only : dp
  use constants,     only : zero, half, one, two, four
  use functions,     only : first_derivative, second_derivative

  implicit none

!User defined constants
  integer :: Re, imax, jmax, maxiter, i, j, var, equ, n

!Derived constants
  real(dp) :: mu, nu, dx, dy, dtd, dtc

!Main loop variables
  real(dp) :: dudx, dudy, dvdx, dvdy, dpdx, dpdy
  real(dp) :: d2udx2, d2udy2, d2vdx2, d2vdy2
  real(dp) :: curvaturepx, averagepx, curvaturepy, averagepy

!Residual loop variables
  real(dp) :: Pweightfactor

!Constants defined by the problem statement
  real(dp), parameter :: L    = 0.05_dp  !Characteristic length, m
  real(dp), parameter :: rho  = 1.0_dp   !Density, kg/m^3
  real(dp), parameter :: Ulid = 1.0_dp   !Lid velocity, m/s
  real(dp), parameter :: Po   = 0.0_dp   !Gauge value of pressure
  real(dp), parameter :: xmin = 0.0_dp   !Min x domain value
  real(dp), parameter :: xmax = 0.05_dp  !Max x domain value
  real(dp), parameter :: ymin = 0.0_dp   !Min y domain value
  real(dp), parameter :: ymax = 0.05_dp  !Max y domain value
  real(dp), parameter :: k    = 0.01_dp  !Scale for Beta
  real(dp), parameter :: C2   = 0.01_dp  !Artificial viscosity weighting term
  real(dp), parameter :: CFL  = 1.0_dp   !CFL number
  real(dp), parameter :: conv    = 1e-9  !Convergence criteria
  real(dp), parameter :: epsilon = 1e-7  !Small term to avoid division by zero
                                         !in artifical viscosity

  real(dp), dimension(3) :: R, L1, L2, Linf
  real(dp), dimension(3,3) :: Ident
  real(dp), allocatable, dimension(:,:,:) :: q, qn, Low, Diag, Up
  real(dp), allocatable, dimension(:,:)   :: dt, B, RHS

  continue

  print *,'What Reynolds Number? '
  read *,Re
  print *,'How many X nodes? '
  read *,imax
  print *,'How many Y nodes? '
  read *,jmax
  print *,'How many iterations? '
  read *,maxiter

  mu = rho*Ulid*L/Re     !Dynamic viscosity, Pa*s
  nu = mu/rho            !Kinematic viscosity, m^2/s
  dx = (xmax-xmin)/real(imax-1,dp)
  dy = (ymax-ymin)/real(jmax-1,dp)
  dtd = dx**2/(four*nu)

  Ident = reshape((/one, zero, zero, zero, one, zero, zero, zero, one/),(/3,3/))

  allocate(q(3,imax,jmax))
  allocate(qn(3,imax,jmax))
  allocate(dt(imax,jmax), B(imax,jmax))

  allocate(Low(3,3,jmax), Diag(3,3,jmax), Up(3,3,jmax), RHS(3,jmax))

!Sets Boundary Conditions and initialize matrices

  q(2,:,:) = zero
  q(3,:,:) = zero
  q(1,:,:) = Po

  qn = zero
  dt = zero

  q(2,:,jmax) = Ulid

  do n = 1,maxiter
    R(:)    = zero
    L1(:)   = zero
    L2(:)   = zero
    Linf(:) = zero

!$OMP DO
!write(*,*) n
    do j = 2,jmax-1
      do i = 2,imax-1

! Set artificial compressibility factor
        B(i,j) = sqrt( max(q(2,i,j)**2, k*Ulid**2) )

! Set time_step
        dtc = dx /                                                             &
          max((abs(q(2,i,j) + sqrt(q(2,i,j)**2 + four*B(i,j)**2)))/two,        &
              (abs(q(3,i,j) + sqrt(q(3,i,j)**2 + four*B(i,j)**2)))/two)
        dt(i,j) = CFL*(one/((one/dtc) + (one/dtd)))

      end do
    end do

! Loop over each line for y implicit
    do i = 2,imax-1

! Set bottom wall
      Low(:,:,1)  = zero
      Up(:,:,1)   = zero
      Up(1,1,1)   = -two
      Diag(:,:,1) = Ident

      RHS(1,1) = -q(1,i,3)
      RHS(2,1) = zero
      RHS(3,1) = zero

! Set top wall
      Up(:,:,jmax)   = zero
      Low(:,:,jmax)  = zero
      Low(1,1,jmax)  = -two
      Diag(:,:,jmax) = Ident

      RHS(1,jmax) = -q(1,i,jmax-2)
      RHS(2,jmax) = Ulid
      RHS(3,jmax) = zero

      do j = 2,jmax-1
! Form tridiagonal system
! X direction, fast
        dpdx   = first_derivative(dx, q(1,i-1,j), q(1,i+1,j))
        curvaturepx = abs(q(1,i+1,j) - two*q(1,i,j) + q(1,i-1,j))
        averagepx   = abs(q(1,i+1,j) + two*q(1,i,j) + q(1,i-1,j))

        dudx   = first_derivative(dx, q(2,i-1,j), q(2,i+1,j))
        d2udx2 = (q(2,i+1,j) - two*q(2,i,j) + q(2,i-1,j))/(dx**2)
        dvdx   = first_derivative(dx, q(3,i-1,j), q(3,i+1,j))
        d2vdx2 = (q(3,i+1,j) - two*q(3,i,j) + q(3,i-1,j))/(dx**2)

! Y direction, slow
        curvaturepy = abs(q(1,i,j+1) - two*q(1,i,j) + q(1,i,j-1))
        averagepy   = abs(q(1,i,j+1) + two*q(1,i,j) + q(1,i,j-1))

        Low(:,:,j) = zero
        Low(1,3,j) = -half*rho/dy - (rho*C2/dy)*curvaturepy/(averagepy+epsilon)
        Low(2,2,j) = -half*rho*q(3,i,j)/dy - mu/dy**2
        Low(3,1,j) = -half/dy
        Low(3,3,j) = -half*rho*q(3,i,j)/dy - mu/dy**2

        Diag(:,:,j) = zero
        Diag(1,1,j) = one/(dt(i,j)*B(i,j)**2)
        Diag(1,3,j) = two*rho*C2*curvaturepy/(dy*(averagepy+epsilon))
        Diag(2,2,j) = rho/dt(i,j) + two*mu*(one/dy**2 + one/dx**2)
        Diag(3,3,j) = rho/dt(i,j) + two*mu*(one/dy**2 + one/dx**2)

        Up(:,:,j) = zero
        Up(1,3,j) = half*rho/dy - (rho*C2/dy)*curvaturepy/(averagepy+epsilon)
        Up(2,2,j) = half*rho*q(3,i,j)/dy - mu/dy**2
        Up(3,1,j) = half/dy
        Up(3,3,j) = half*rho*q(3,i,j)/dy - mu/dy**2

        RHS(1,j) = q(1,i,j)/(dt(i,j)*B(i,j)**2) - rho*dudx                     &
                 + (rho*C2*dx*d2udx2*curvaturepx)/(averagepx+epsilon)
        RHS(2,j) = rho*q(2,i,j)/dt(i,j) - rho*q(2,i,j)*dudx - dpdx             &
                 + mu*(q(2,i+1,j)+q(2,i-1,j))/dx**2
        RHS(3,j) = rho*q(3,i,j)/dt(i,j) - rho*q(2,i,j)*dvdx                    &
                 + mu*(q(3,i+1,j)+q(3,i-1,j))/dx**2
      end do
! Solve the line implicit system
      call triblocksolve(3, jmax, Low, Diag, Up, RHS, qn(:,i,:))
    end do
!$OMP END DO

!    L1(:) =      L1(:)  / real((imax-2)*(jmax-2),dp)
!    L2(:) = sqrt(L2(:)) / real((imax-2)*(jmax-2),dp)

!Calculate side wall pressures
    do j = 1,jmax
      qn(1,1,j)    = two*qn(1,2,j)      - qn(1,3,j)
      qn(1,imax,j) = two*qn(1,imax-1,j) - qn(1,imax-2,j)
    end do

!Pressure rescaling at the center point of the bottom floor
    Pweightfactor = qn(1,(imax-imax/2),1) - Po
    qn(1,:,:) = qn(1,:,:) - Pweightfactor
  
    q = qn

!Residual Calculations
    if (floor(n/1000.0_dp) == n/1000.0_dp) then
      write(*,300) n, L2(1), L2(2), L2(3)
300   format(1X,i8,2(e15.6),3(e15.6),4(e15.6))
    end if
  end do

  open(20, FILE='Lid_Driven_Cavity_Implicit.tec', status='unknown')

  write(20,*) 'TITLE="Lid Driven Cavity Implicit"'
  write(20,*) 'VARIABLES = "X", "Y", "Pressure", "U", "V","S1","S2","S3"'
  write(20,*) 'ZONE DATAPACKING=BLOCK, I=', IMAX, ', J=', JMAX
  write(20,*) 

  do j = 1, jmax
    do i = 1, imax
      write(20,*) (real(i-1,dp)/real(imax-1))
    end do
  end do

  do j = 1, jmax
    do i = 1, imax
      write(20,*) (real(j-1,dp)/real(jmax-1))
    end do
  end do

  do var = 1,3
    do j = 1, jmax
      do i = 1, imax
        write(20,*) q(var, i, j)
      end do
    end do
  end do

  do j = 1, jmax
    do i = 1, imax
      write(20,*) abs(2.5_dp*( 1.4_dp - log(q(1,i,j)+101325._dp) ) &
                - 0.5_dp*(q(2,i,j)**2+q(3,i,j)**2)/(q(1,i,j)+101325._dp))
    end do
  end do

  do j = 1, jmax
    do i = 1, imax
      write(20,*) abs(q(2,i,j) / (q(1,i,j)+101325._dp))
    end do
  end do

  do j = 1, jmax
    do i = 1, imax
      write(20,*) abs(q(3,i,j) / (q(1,i,j)+101325._dp))
    end do
  end do

  close(20)

end program LidDrivenCavityImplicit
