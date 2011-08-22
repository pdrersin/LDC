!  This code calculates the flow in a lid driven cavity
module set_precision

  implicit none
  
  integer, parameter :: Sngl = Selected_Real_Kind(6 , 37)
  integer, parameter :: Dbl  = Selected_Real_Kind(15, 307)
  integer, parameter :: Quad = Selected_Real_Kind(33, 4931)
  integer, parameter :: dp = Dbl

end module

module set_constants

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

end module set_constants

module functions

  implicit none

  contains

function first_derivative(dl, var_l, var_r)

  use set_precision, only : dp
  use set_constants, only : half

  real(dp), intent(in) :: dl, var_l, var_r
  real(dp)             :: first_derivative

continue

  first_derivative = half*(var_r-var_l)/dl

end function first_derivative

function second_derivative(dl, var_l, var_c, var_r)

  use set_precision, only : dp
  use set_constants, only : two

  real(dp), intent(in) :: dl, var_l, var_c, var_r
  real(dp)             :: second_derivative

continue

  second_derivative = (var_r-two*var_c+var_l)/dl**2

end function second_derivative
end module functions

program LidDrivenCavityExplicit

  use set_precision, only : dp
  use set_constants, only : zero, one, two, four
  use functions,     only : first_derivative, second_derivative

  implicit none

!User defined constants
  integer :: Re, imax, jmax, maxiter, i, j, equ, n, var

!Derived constants
  real(dp) :: mu, nu, dx, dy, dtd, dtc

!Main loop variables
  real(dp) :: dudx, dudy, dvdx, dvdy, dpdx, dpdy
  real(dp) :: d2udx2, d2udy2, d2vdx2, d2vdy2
  real(dp) :: curvaturepx, averagepx, curvaturepy, averagepy, B, S

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
  real(dp), parameter :: CFL  = 0.5_dp   !CFL number
  real(dp), parameter :: conv    = 1e-9  !Convergence criteria
  real(dp), parameter :: epsilon = 1e-7  !Small term to avoid division by zero
                                         !in artifical viscosity

  real(dp), dimension(3) :: R, L1, L2, Linf
  real(dp), allocatable, dimension(:,:,:) :: q, qn
  real(dp), allocatable, dimension(:,:)   :: dt

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

  allocate(q(3,imax,jmax))
  allocate(qn(3,imax,jmax))
  allocate(dt(imax,jmax))

!Sets Boundary Conditions and initialize matrices

  q(2,:,:) = zero
  q(3,:,:) = zero
  q(1,:,:) = Po
  q(2,:,jmax) = Ulid
  qn = q

  dt = zero

  do n = 1,maxiter
    R(:)    = zero
    L1(:)   = zero
    L2(:)   = zero
    Linf(:) = zero

!$OMP DO
    do j = 2,jmax-1
      do i = 2,imax-1

! Set artificial compressibility factor
        B = sqrt( max(q(2,i,j)**2, k*Ulid**2) )

! Set time_step
        dtc = two*dx /                                                          &
          max((abs(q(2,i,j) + sqrt(q(2,i,j)**2 + four*B**2))),                  &
              (abs(q(3,i,j) + sqrt(q(3,i,j)**2 + four*B**2))))
        dt(i,j) = CFL*(one/((one/dtc) + (one/dtd)))

! create_residual(dt(i,j), q, R(:))
! X direction, fast
        dpdx   = first_derivative(dx, q(1,i-1,j), q(1,i+1,j))
        curvaturepx = abs(q(1,i+1,j) - two*q(1,i,j) + q(1,i-1,j))
        averagepx   = abs(q(1,i+1,j) + two*q(1,i,j) + q(1,i-1,j))

        dudx   = first_derivative( dx, q(2,i-1,j),           q(2,i+1,j))
        d2udx2 = second_derivative(dx, q(2,i-1,j), q(2,i,j), q(2,i+1,j))
        dvdx   = first_derivative( dx, q(3,i-1,j),           q(3,i+1,j))
        d2vdx2 = second_derivative(dx, q(3,i-1,j), q(3,i,j), q(3,i+1,j))

! Y direction, slow
        dpdy   = first_derivative(dy, q(1,i,j-1), q(1,i,j+1))
        curvaturepy = abs(q(1,i,j+1) - two*q(1,i,j) + q(1,i,j-1))
        averagepy   = abs(q(1,i,j+1) + two*q(1,i,j) + q(1,i,j-1))

        dudy   = first_derivative( dy, q(2,i,j-1),           q(2,i,j+1))
        d2udy2 = second_derivative(dy, q(2,i,j-1), q(2,i,j), q(2,i,j+1))
        dvdy   = first_derivative( dy, q(3,i,j-1),           q(3,i,j+1))
        d2vdy2 = second_derivative(dy, q(3,i,j-1), q(3,i,j), q(3,i,j+1))
      
        S = rho*C2*(dx*(curvaturepx/(epsilon+averagepx))*d2udx2 +              &
                    dy*(curvaturepy/(epsilon+averagepy))*d2vdy2)

        R(1) = (rho*(dudx+dvdy) - S)*B**2      
        R(2) = q(2,i,j)*dudx + q(3,i,j)*dudy + dpdx/rho - nu*(d2udx2 + d2udy2)
        R(3) = q(2,i,j)*dvdx + q(3,i,j)*dvdy + dpdy/rho - nu*(d2vdx2 + d2vdy2)
      
        do equ = 1,3
          L1(equ)   = L1(equ) + R(equ)
          L2(equ)   = L2(equ) + R(equ)**2
          Linf(equ) = max(Linf(equ), R(equ))
        end do

!        call euler_explicit_iteration(dt(i,j), R(:), q(:,i,j), qn(:,i,j))
        qn(1,i,j) = q(1,i,j) - dt(i,j)*R(1)
        qn(2,i,j) = q(2,i,j) - dt(i,j)*R(2)
        qn(3,i,j) = q(3,i,j) - dt(i,j)*R(3)
      end do
    end do
!$OMP END DO

    L1(:) =      L1(:)  / real((imax-2)*(jmax-2),dp)
    L2(:) = sqrt(L2(:)) / real((imax-2)*(jmax-2),dp)

!Calculate side wall pressures
    do j = 1,jmax
      qn(1,1,j)    = two*qn(1,2,j)      - qn(1,3,j)
      qn(1,imax,j) = two*qn(1,imax-1,j) - qn(1,imax-2,j)
    end do

    do i = 1,imax
      qn(1,i,1)    = two*qn(1,i,2)      - qn(1,i,3)
      qn(1,i,jmax) = two*qn(1,i,jmax-1) - qn(1,i,jmax-2)
    end do

!Pressure rescaling at the center point of the bottom floor
    Pweightfactor = qn(1,(imax-imax/2),1) - Po
    qn(1,:,:) = qn(1,:,:) - Pweightfactor
  
    q = qn

!Residual Calculations
    if (mod(n,1000) == 0) then
      write(*,300) n, L2(1), L2(2), L2(3)
300   format(1X,i8,2(e15.6),3(e15.6),4(e15.6))
    end if
  end do

  open(20, FILE='Lid_Driven_Cavity_Explicit.tec', status='unknown')

  write(20,*) 'TITLE="Lid Driven Cavity Explicit"'
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

  deallocate(q, qn, dt)

end program LidDrivenCavityExplicit
