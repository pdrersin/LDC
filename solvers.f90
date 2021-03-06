module solvers

  implicit none

  private

  public :: ldc_explicit
  public :: ldc_implicit
  public :: ldc_sgs

contains

!=============================== ldc_explicit ================================80
!
! The main routine for the explicit solve
!
!=============================================================================80
  subroutine ldc_explicit( x_nodes, y_nodes, dx, dy, soln, soln_new )

    use set_precision, only : dp
    use set_constants, only : zero, two
    use setup,         only : max_iter, dtd, cfl, k, u_lid, p_guage,           &
                              conv_toler, resid_out

    integer,                                  intent(in)    :: x_nodes, y_nodes
    real(dp),                                 intent(in)    :: dx, dy
    real(dp), dimension(3, x_nodes, y_nodes), intent(inout) :: soln, soln_new

    real(dp), dimension(3) :: R, L1, L2, Linf

    integer  :: iter, i, j, eq
    real(dp) :: Pweightfactor, internal_nodes

    real(dp), dimension(x_nodes, y_nodes) :: dt, beta

    continue

    internal_nodes = real((x_nodes-2)*(y_nodes-2),dp)

    iter_loop : do iter = 1, max_iter
      R    = zero
      L1   = zero
      L2   = zero
      Linf = zero

! Calculate artifical compressibility terms
      do j = 2, y_nodes-1
        do i = 2, x_nodes-1
          beta(i,j) = set_beta(soln(2,i,j), k, u_lid)
        end do
      end do

! Calculate local timestep, separate loops seems to be better for cache
      do j = 2, y_nodes-1
        do i = 2, x_nodes-1
          dt(i,j) = set_dt(dx, dtd, cfl, soln(2,i,j), soln(3,i,j), beta(i,j))
        end do
      end do

! Residual loop

!$omp parallel &
!$omp private(i, j, eq, R) &
!$omp reduction(+ : L1, L2) &
!$omp reduction(max : Linf)
  !$omp do
      do j = 2, y_nodes-1
        do i = 2, x_nodes-1

          R = create_residual(i, j, x_nodes, y_nodes, dx, dy, beta(i,j), soln )

! Update residual norms
          L1   = L1 + R
          L2   = L2 + R**2
          Linf = max(Linf, R)

! Euler explicit solve
          do eq = 1,3
            soln_new(eq,i,j) = soln(eq,i,j) - dt(i,j)*R(eq)
          end do

        end do
      end do
  !$omp end do
!$omp end parallel

! Update L1 and L2 residuals
      L1 =      L1  / internal_nodes
      L2 = sqrt(L2) / internal_nodes

! Update soln
      soln = soln_new

!Calculate side wall pressures
      do j = 1,y_nodes
        soln(1,1,j)       = two*soln(1,2,j)         - soln(1,3,j)
        soln(1,x_nodes,j) = two*soln(1,x_nodes-1,j) - soln(1,x_nodes-2,j)
      end do

      do i = 1,x_nodes
        soln(1,i,1)       = two*soln(1,i,2)         - soln(1,i,3)
        soln(1,i,y_nodes) = two*soln(1,i,y_nodes-1) - soln(1,i,y_nodes-2)
      end do

!Pressure rescaling at the center point of the bottom floor
      Pweightfactor = soln(1,(x_nodes-x_nodes/2),1) - p_guage
      soln(1,:,:)   = soln(1,:,:) - Pweightfactor

!Residual Calculations
      if ( mod(iter,resid_out) == 0 ) then
        write(*,300) iter, L2(1), L2(2), L2(3)
300     format(1X,i8,3e15.6)

        if ( L2(2) <= conv_toler .and. L2(3) <= conv_toler ) then
          write(*,*) "Solution has converged"
          return
        end if
      end if

    end do iter_loop

  end subroutine ldc_explicit

!================================= ldc_sgs ===================================80
!
! The main routine for the symmetric Gauss-Seidel solve
!
!=============================================================================80
  subroutine ldc_sgs( x_nodes, y_nodes, dx, dy, soln )

    use set_precision, only : dp
    use set_constants, only : zero, two
    use setup,         only : max_iter, dtd, cfl, k, u_lid, p_guage,           &
                              conv_toler, resid_out

    integer,                                  intent(in)    :: x_nodes, y_nodes
    real(dp),                                 intent(in)    :: dx, dy
    real(dp), dimension(3, x_nodes, y_nodes), intent(inout) :: soln

    real(dp), dimension(3) :: R, L1, L2, Linf

    integer  :: iter, i, j, eq
    real(dp) :: Pweightfactor, internal_nodes

    real(dp), dimension(x_nodes, y_nodes) :: dt, beta

    continue

    internal_nodes = real((x_nodes-2)*(y_nodes-2),dp)

    iter_loop : do iter = 1, max_iter
      R    = zero
      L1   = zero
      L2   = zero
      Linf = zero

! Forward sweep !

! Calculate artifical compressibility terms
      do j = 2, y_nodes-1
        do i = 2, x_nodes-1
          beta(i,j) = set_beta(soln(2,i,j), k, u_lid)
        end do
      end do

! Calculate local timestep, separate loops seems to be better for cache
      do j = 2, y_nodes-1
        do i = 2, x_nodes-1
          dt(i,j) = set_dt(dx, dtd, cfl, soln(2,i,j), soln(3,i,j), beta(i,j))
        end do
      end do

      f_sgs : do j = 2, y_nodes-1
        do i = 2, x_nodes-1
          R = create_residual(i, j, x_nodes, y_nodes, dx, dy, beta(i,j), soln )
          soln(:,i,j) = soln(:,i,j) - dt(i,j)*R(:)
        end do
      end do f_sgs

! Calculate side wall pressures
      do j = 1,y_nodes
        soln(1,1,j)       = two*soln(1,2,j)         - soln(1,3,j)
        soln(1,x_nodes,j) = two*soln(1,x_nodes-1,j) - soln(1,x_nodes-2,j)
      end do

      do i = 1,x_nodes
        soln(1,i,1)       = two*soln(1,i,2)         - soln(1,i,3)
        soln(1,i,y_nodes) = two*soln(1,i,y_nodes-1) - soln(1,i,y_nodes-2)
      end do

! Pressure rescaling at the center point of the bottom floor
      Pweightfactor = soln(1,(x_nodes-x_nodes/2),1) - p_guage
      soln(1,:,:)   = soln(1,:,:) - Pweightfactor

! Backward sweep !

! Calculate artifical compressibility terms
      do j = 2, y_nodes-1
        do i = 2, x_nodes-1
          beta(i,j) = set_beta(soln(2,i,j), k, u_lid)
        end do
      end do

! Calculate local timestep, separate loops seems to be better for cache
      do j = 2, y_nodes-1
        do i = 2, x_nodes-1
          dt(i,j) = set_dt(dx, dtd, cfl, soln(2,i,j), soln(3,i,j), beta(i,j))
        end do
      end do

      b_sgs : do i = 2, x_nodes-1
        do j = 2, y_nodes-1
          R = create_residual(i, j, x_nodes, y_nodes, dx, dy, beta(i,j), soln )
! Update residual norms
          L1   = L1 + R
          L2   = L2 + R**2
          Linf = max(Linf, R)

! Update solution
          soln(:,i,j) = soln(:,i,j) - dt(i,j)*R(:)
        end do
      end do b_sgs

! Update L1 and L2 residuals
      L1(:) =      L1(:)  / internal_nodes
      L2(:) = sqrt(L2(:)) / internal_nodes

!Calculate side wall pressures
      do j = 1,y_nodes
        soln(1,1,j)       = two*soln(1,2,j)         - soln(1,3,j)
        soln(1,x_nodes,j) = two*soln(1,x_nodes-1,j) - soln(1,x_nodes-2,j)
      end do

      do i = 1,x_nodes
        soln(1,i,1)       = two*soln(1,i,2)         - soln(1,i,3)
        soln(1,i,y_nodes) = two*soln(1,i,y_nodes-1) - soln(1,i,y_nodes-2)
      end do
!Pressure rescaling at the center point of the bottom floor
      Pweightfactor = soln(1,(x_nodes-x_nodes/2),1) - p_guage
      soln(1,:,:)   = soln(1,:,:) - Pweightfactor

!Residual Calculations
      if ( mod(iter,resid_out) == 0 ) then
        write(*,300) iter, L2(1), L2(2), L2(3)
300     format(1X,i8,3e15.6)

        if ( L2(2) <= conv_toler .and. L2(3) <= conv_toler ) then
          write(*,*) "Solution has converged"
          return
        end if
      end if

    end do iter_loop

  end subroutine ldc_sgs

!=============================== create_residual =============================80
!
! Forms the ldc residual, it is a pure, inlineable function
!
!=============================================================================80
  pure function create_residual( i, j, x_nodes, y_nodes, dx, dy, beta, soln )

    use set_precision, only : dp
    use set_constants, only : two
    use setup,         only : rho, nu, visc_eps, c2

    integer,                                intent(in) :: i, j, x_nodes, y_nodes
    real(dp),                               intent(in) :: dx, dy, beta
    real(dp), dimension(3,x_nodes,y_nodes), intent(in) :: soln

    real(dp), dimension(3)                             :: create_residual

    real(dp) :: S, dpdx, crvpx, avgpx, dpdy, crvpy, avgpy
    real(dp) :: dudx, d2udx2, dvdx, d2vdx2, dudy, d2udy2, dvdy, d2vdy2

    continue

! X direction, fast
    dpdx  = first_derivative(dx, soln(1,i-1,j), soln(1,i+1,j))
    crvpx = abs(soln(1,i+1,j) - two*soln(1,i,j) + soln(1,i-1,j))
    avgpx = abs(soln(1,i+1,j) + two*soln(1,i,j) + soln(1,i-1,j))

    dudx   = first_derivative( dx, soln(2,i-1,j),              soln(2,i+1,j))
    d2udx2 = second_derivative(dx, soln(2,i-1,j), soln(2,i,j), soln(2,i+1,j))
    dvdx   = first_derivative( dx, soln(3,i-1,j),              soln(3,i+1,j))
    d2vdx2 = second_derivative(dx, soln(3,i-1,j), soln(3,i,j), soln(3,i+1,j))

! Y direction, slow
    dpdy  = first_derivative(dy, soln(1,i,j-1), soln(1,i,j+1))
    crvpy = abs(soln(1,i,j+1) - two*soln(1,i,j) + soln(1,i,j-1))
    avgpy = abs(soln(1,i,j+1) + two*soln(1,i,j) + soln(1,i,j-1))

    dudy   = first_derivative( dy, soln(2,i,j-1),              soln(2,i,j+1))
    d2udy2 = second_derivative(dy, soln(2,i,j-1), soln(2,i,j), soln(2,i,j+1))
    dvdy   = first_derivative( dy, soln(3,i,j-1),              soln(3,i,j+1))
    d2vdy2 = second_derivative(dy, soln(3,i,j-1), soln(3,i,j), soln(3,i,j+1))

! Artificial viscosity term
    S = rho*c2*(dx*(crvpx/(visc_eps+avgpx))*d2udx2 +                           &
                dy*(crvpy/(visc_eps+avgpy))*d2vdy2)

! Residual
    create_residual(1) = (rho*(dudx+dvdy) - S)*beta**2
    create_residual(2) = soln(2,i,j)*dudx + soln(3,i,j)*dudy                   &
                       + dpdx/rho - nu*(d2udx2 + d2udy2)
    create_residual(3) = soln(2,i,j)*dvdx + soln(3,i,j)*dvdy                   &
                       + dpdy/rho - nu*(d2vdx2 + d2vdy2)

  end function create_residual

!============================= ldc_full_implicit =============================80
!
! The main routine for the full implicit solve, requires SPARSKIT2
!
!=============================================================================80
  subroutine ldc_full_implicit( x_nodes, y_nodes, dx, dy, soln, soln_new )

    use set_precision, only : dp
    use set_constants, only : zero, two
    use setup,         only : max_iter, dtd, cfl, k, u_lid, p_guage,           &
                              conv_toler, resid_out

    integer,                                  intent(in)    :: x_nodes, y_nodes
    real(dp),                                 intent(in)    :: dx, dy
    real(dp), dimension(3, x_nodes, y_nodes), intent(inout) :: soln, soln_new

    real(dp), dimension(3) :: R, L1, L2, Linf

    integer  :: iter, i, j, eq, n_nz, n_col, n_row
    real(dp) :: Pweightfactor, internal_nodes


    integer,  dimension(3*x_nodes*y_nodes) :: ia
    integer,  allocatable, dimension(:)    :: ja
    real(dp), dimension(3*x_nodes*y_nodes) :: rhs, dsoln
    real(dp), allocatable, dimension(:)    :: lhs
    real(dp), dimension(x_nodes, y_nodes)  :: dt, beta

    continue

    n_row = 3*x_nodes*y_nodes
    n_col = n_row
    n_nz = 9*(5*(x_nodes-2)*(y_nodes-2) + 6*(x_nodes + y_nodes - 4) + 8)

    allocate( ja(n_nz), lhs(n_nz) )

    deallocate(lhs, ja)

  end subroutine ldc_full_implicit

!=============================== ldc_implicit ================================80
!
! The main routine for the implicit solve
!
!=============================================================================80
  subroutine ldc_implicit( x_nodes, y_nodes, dx, dy, soln, soln_new )

    use set_precision, only : dp
    use set_constants, only : zero, two
    use setup,         only : max_iter, dtd, cfl, k, u_lid, p_guage,           &
                              conv_toler, resid_out
    use matrix_manip,  only : triblocksolve

    integer,                                  intent(in)    :: x_nodes, y_nodes
    real(dp),                                 intent(in)    :: dx, dy
    real(dp), dimension(3, x_nodes, y_nodes), intent(inout) :: soln, soln_new

    real(dp), dimension(3)                :: R, L2
    real(dp), dimension(3, y_nodes)       :: RHS
    real(dp), dimension(x_nodes, y_nodes) :: dt, beta
    real(dp), dimension(3, 3, y_nodes)    :: Low, Diag, Up

    integer  :: iter, i, j, eq
    real(dp) :: Pweightfactor, internal_nodes

    continue

    internal_nodes = real((x_nodes-2)*(y_nodes-2),dp)

    Low  = zero
    Diag = zero
    Up   = zero
    RHS  = zero

    iter_loop : do iter = 1, max_iter
      R    = zero
      L2   = zero

! Calculate artifical compressibility terms
      do j = 2, y_nodes-1
        do i = 2, x_nodes-1
          beta(i,j) = set_beta(soln(2,i,j), k, u_lid)
        end do
      end do

! Calculate local timestep, separate loops seems to be better for cache
      do j = 2, y_nodes-1
        do i = 2, x_nodes-1
          dt(i,j) = set_dt(dx, dtd, cfl, soln(2,i,j), soln(3,i,j), beta(i,j))
        end do
      end do

! Residual loop

!$omp parallel &
!$omp private(i, Low, Diag, Up, RHS)
  !$omp do
      do i = 2, x_nodes-1
! Form LHS
        call lhs_y_implicit(i, x_nodes, y_nodes, dx, dy, dt(i,:), beta(i,:),   &
                            soln, Low, Diag, Up)
! Form RHS
        call rhs_y_implicit(i, x_nodes, y_nodes, dx, dt(i,:), beta(i,:),       &
                            soln, RHS)
! Modify structure for BC
        call bc_mod_y_implicit(y_nodes, Low, Diag, Up, RHS)

! Solve the line implicit system
        call triblocksolve(3, y_nodes, Low, Diag, Up, RHS, soln_new(:,i,:))

      end do
  !$omp end do
!$omp end parallel

! Update L1 and L2 residuals

      do i = 2, x_nodes-1
        do j = 2, y_nodes-1
          L2(1) = L2(1) + ((soln_new(1,i,j)-soln(1,i,j))/(beta(i,j)*dt(i,j)))**2
          L2(2) = L2(2) + ((soln_new(2,i,j)-soln(2,i,j))/dt(i,j))**2
          L2(3) = L2(3) + ((soln_new(3,i,j)-soln(3,i,j))/dt(i,j))**2
        end do
      end do
      L2 = sqrt(L2) / internal_nodes

! Update soln
      soln = soln_new

!Calculate side wall pressures
      do j = 1,y_nodes
        soln(1,1,j)       = two*soln(1,2,j)         - soln(1,3,j)
        soln(1,x_nodes,j) = two*soln(1,x_nodes-1,j) - soln(1,x_nodes-2,j)
      end do

!Pressure rescaling at the center point of the bottom floor
      Pweightfactor = soln(1,(x_nodes-x_nodes/2),1) - p_guage
      soln(1,:,:)   = soln(1,:,:) - Pweightfactor

!Residual Calculations
      if ( mod(iter,resid_out) == 0 ) then
        write(*,300) iter, L2(1), L2(2), L2(3)
300     format(1X,i8,3e15.6)

        if ( L2(2) <= conv_toler .and. L2(3) <= conv_toler ) then
          write(*,*) "Solution has converged"
          return
        end if
      end if

    end do iter_loop

  end subroutine ldc_implicit

!=============================================================================80
!
!
!
!=============================================================================80
  subroutine rhs_y_implicit( i, x_nodes, y_nodes, dx, dt, beta, soln, RHS )

    use set_precision, only : dp
    use set_constants, only : zero, two
    use setup,         only : rho, mu, u_lid, visc_eps, c2

    integer,                                intent(in)  :: i, x_nodes, y_nodes
    real(dp),                               intent(in)  :: dx
    real(dp), dimension(3,x_nodes,y_nodes), intent(in)  :: soln
    real(dp), dimension(1,y_nodes),         intent(in)  :: dt
    real(dp), dimension(1,y_nodes),         intent(in)  :: beta
    real(dp), dimension(3,y_nodes),         intent(out) :: RHS

    integer  :: j
    real(dp) :: dpdx, crvpx, avgpx, dudx, d2udx2, dvdx

    continue

    RHS(1,1) = zero
    RHS(2,1) = zero
    RHS(3,1) = zero

    do j = 2, y_nodes-1

      dpdx  = first_derivative(dx, soln(1,i-1,j), soln(1,i+1,j))
      crvpx = abs(soln(1,i+1,j) - two*soln(1,i,j) + soln(1,i-1,j))
      avgpx = abs(soln(1,i+1,j) + two*soln(1,i,j) + soln(1,i-1,j))

      dudx   = first_derivative( dx, soln(2,i-1,j),              soln(2,i+1,j))
      d2udx2 = second_derivative(dx, soln(2,i-1,j), soln(2,i,j), soln(2,i+1,j))

      dvdx   = first_derivative(dx, soln(3,i-1,j), soln(3,i+1,j))

      RHS(1,j) = soln(1,i,j)/(dt(1,j)*beta(1,j)**2) - rho*dudx                 &
               + (rho*C2*dx*d2udx2*crvpx)/(avgpx+visc_eps)

      RHS(2,j) = rho*soln(2,i,j)/dt(1,j) - rho*soln(2,i,j)*dudx - dpdx         &
               + mu*(soln(2,i+1,j)+soln(2,i-1,j))/dx**2

      RHS(3,j) = rho*soln(3,i,j)/dt(1,j) - rho*soln(2,i,j)*dvdx                &
               + mu*(soln(3,i+1,j)+soln(3,i-1,j))/dx**2

    end do

    RHS(1,y_nodes) = zero
    RHS(2,y_nodes) = u_lid
    RHS(3,y_nodes) = zero

  end subroutine rhs_y_implicit

!=============================================================================80
!
!
!
!=============================================================================80
  subroutine lhs_y_implicit( i, x_nodes, y_nodes, dx, dy, dt, beta,            &
                             soln, Low, Diag, Up )

    use set_precision, only : dp
    use set_constants, only : zero, half, one, two
    use setup,         only : rho, mu, visc_eps, c2

    integer,                                  intent(in)  :: i, x_nodes, y_nodes
    real(dp),                                 intent(in)  :: dx, dy
    real(dp), dimension(1,y_nodes),           intent(in)  :: dt, beta
    real(dp), dimension(3, x_nodes, y_nodes), intent(in)  :: soln
    real(dp), dimension(3, 3, y_nodes),       intent(out) :: Low, Diag, Up

    integer  :: j
    real(dp) :: crvpy, avgpy
    real(dp), dimension(3,3) :: Ident

    continue

    Ident = reshape( [one, zero, zero, zero, one, zero, zero, zero, one],      &
                     [3,3] )

    Low(:,:,1)  = zero
    Diag(:,:,1) = Ident
    Up(:,:,1)   = zero
    Up(1,1,1)   = -two

! Set interior
    do j = 2,y_nodes-1

      crvpy = abs(soln(1,i,j+1) - two*soln(1,i,j) + soln(1,i,j-1))
      avgpy = abs(soln(1,i,j+1) + two*soln(1,i,j) + soln(1,i,j-1))

      Low(1,1,j) = zero
      Low(2,1,j) = zero
      Low(3,1,j) = -half/dy
      Low(1,2,j) = zero
      Low(2,2,j) = -half*rho*soln(3,i,j)/dy - mu/dy**2
      Low(3,2,j) = zero
      Low(1,3,j) = -half*rho/dy - (rho*C2/dy)*crvpy/(avgpy+visc_eps)
      Low(2,3,j) = zero
      Low(3,3,j) = -half*rho*soln(3,i,j)/dy - mu/dy**2

      Diag(1,1,j) = one/(dt(1,j)*beta(1,j)**2)
      Diag(2,1,j) = zero
      Diag(3,1,j) = zero
      Diag(1,2,j) = zero
      Diag(2,2,j) = rho/dt(1,j) + two*mu*(one/dy**2 + one/dx**2)
      Diag(3,2,j) = zero
      Diag(1,3,j) = two*rho*C2*crvpy/(dy*(avgpy+visc_eps))
      Diag(2,3,j) = zero
      Diag(3,3,j) = rho/dt(1,j) + two*mu*(one/dy**2 + one/dx**2)

      Up(1,1,j) = zero
      Up(2,1,j) = zero
      Up(3,1,j) = half/dy
      Up(1,2,j) = zero
      Up(2,2,j) = half*rho*soln(3,i,j)/dy - mu/dy**2
      Up(3,2,j) = zero
      Up(1,3,j) = half*rho/dy - (rho*C2/dy)*crvpy/(avgpy+visc_eps)
      Up(2,3,j) = zero
      Up(3,3,j) = half*rho*soln(3,i,j)/dy - mu/dy**2

    end do
! Set top wall
    Low(:,:,y_nodes)  = zero
    Low(1,1,y_nodes)  = -two
    Diag(:,:,y_nodes) = Ident
    Up(:,:,y_nodes)   = zero

  end subroutine lhs_y_implicit

!============================= bc_mod_y_implicit =============================80
!
! This routine modifies the extrapolated BC eq's to stay block tridiagonal
!
!=============================================================================80
! How it works:
! Starting with...
! |D1 U1 A  0    0    0    0   |       |  RHS1  | Need to eliminate A
! |L2 D2 U2 0    0    0    0   |       |  RHS2  |
! |0  L3 D3 U3   0    0    0   |       |  RHS3  |
! |0  0  L. D.   U.   0    0   |*{Q} = |  RHS.  |
! |0  0  0  LN-2 DN-2 UN-2 0   |       | RHSN-2 |
! |0  0  0  0    LN-1 DN-1 UN-1|       | RHSN-1 |
! |0  0  0  0    B    LN   DN  |       | RHSN   | Need to eliminate B
!
! Deal with A first by temporarily modifying 2nd row
!
! |  D      U    A 0 0 0 0| ... |  RHS1  |
! |AU2'L2 AU2'D2 A 0 0 0 0| ... |AU2'RHS2|
!
! and subtracting to get
!
! |D1-AU2'L2 U1-AU2'D2 0  0 0 0 0| ... |RHS1-AU2'RHS2|
! |    L2        D2    U2 0 0 0 0| ... |     RHS2    |
!
! Performing similar work for the bottom two rows...
!
! |0 0 0 0 LN-1     DN-1           UN-1     | ... |     RHSN-1      |
! |0 0 0 0 0    LN-BLN-1'DN-1  DN-BLN-1'UN-1| ... |RHSN-BLN-1'RHSN-1|
!
! And this preserves the block tridiagonal system

  subroutine bc_mod_y_implicit( y_nodes, L, D, U, RHS )

    use set_precision, only : dp
    use set_constants, only : zero, half, one, two

    integer,                          intent(in)    :: y_nodes
    real(dp), dimension(3,3,y_nodes), intent(inout) :: L, D, U
    real(dp), dimension(3,y_nodes),   intent(inout) :: RHS

    real(dp), dimension(3,3) :: bc, inv

    continue

    bc = zero
    bc(1,1) = one

!    call matrix_inv(3, U(:,:,2), inv)
    inv = mat_inv_3x3( U(:,:,2) )

    inv = matmul(bc, inv)

    D(:,:,1) = D(:,:,1) - matmul(inv, L(:,:,2))
    U(:,:,1) = U(:,:,1) - matmul(inv, D(:,:,2))
    RHS(:,1) = RHS(:,1) - matmul(inv, RHS(:,2))

!    call matrix_inv(3, L(:,:,y_nodes-1), inv)
    inv = mat_inv_3x3( L(:,:,y_nodes-1) )

    inv = matmul(bc, inv)

    L(:,:,y_nodes) = L(:,:,y_nodes) - matmul(inv, D(:,:,y_nodes-1))
    D(:,:,y_nodes) = D(:,:,y_nodes) - matmul(inv, U(:,:,y_nodes-1))
    RHS(:,y_nodes) = RHS(:,y_nodes) - matmul(inv, RHS(:,y_nodes-1))

  end subroutine bc_mod_y_implicit

  include 'set_beta.f90'
  include 'set_dt.f90'
  include 'first_derivative.f90'
  include 'second_derivative.f90'
  include 'mat_inv_3x3.f90'
  include 'det_3x3.f90'

end module solvers
