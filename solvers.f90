module solvers

  implicit none

contains

  subroutine ldc_explicit(x_nodes, y_nodes, dt, beta, soln, soln_new)

    implicit none

    integer,                                  intent(in)    :: x_nodes, y_nodes
    real(dp), dimension(x_nodes, y_nodes),    intent(inout) :: dt, beta
    real(dp), dimension(3, x_nodes, y_nodes), intent(inout) :: soln, soln_new

    real(dp), dimension(3) R, L1, L2, Linf

    continue

    iter_loop : do iter = 1, max_iter
      R(:)    = zero
      L1(:)   = zero
      L2(:)   = zero
      Linf(:) = zero

      do j = 2, y_nodes-1
        do i = 2, x_nodes-1
          beta(i,j) = set_beta(soln(2,i,j), k, u_lid)
          dt(i,j)   = set_dt(dx, dtd, cfl, soln(2,i,j), soln(3,i,j), beta(i,j))
        end do
      end do

      do j = 2, y_nodes-1
        do i = 2, x_nodes-1

          R = create_residual(i, j, x_nodes, y_nodes, dx, dy, beta(i,j), soln )

          do eq = 1,3
            L1(eq)   = L1(eq) + R(eq)
            L2(eq)   = L2(eq) + R(eq)**2
            Linf(eq) = max(Linf(eq), R(eq))
          end do

          do eq = 1,3
            soln_new(eq,i,j) = soln(eq,i,j) - dt(i,j)*R(eq)
          end do

        end do
      end do

! Update L1 and L2 residuals
      L1(:) =      L1(:)  / real((imax-2)*(jmax-2),dp)
      L2(:) = sqrt(L2(:)) / real((imax-2)*(jmax-2),dp)

! Update soln
      soln = soln_new

!Calculate side wall pressures
      do j = 1,jmax
        soln(1,1,j)    = two*soln(1,2,j)      - soln(1,3,j)
        soln(1,imax,j) = two*soln(1,imax-1,j) - soln(1,imax-2,j)
      end do

      do i = 1,imax
        soln(1,i,1)    = two*soln(1,i,2)      - soln(1,i,3)
        soln(1,i,jmax) = two*soln(1,i,jmax-1) - soln(1,i,jmax-2)
      end do

!Pressure rescaling at the center point of the bottom floor
      Pweightfactor = soln(1,(x_nodes-x_nodes/2),1) - p_guage
      soln(1,:,:)   = soln(1,:,:) - Pweightfactor
  
!Residual Calculations
      if (mod(n,1000) == 0) then
        write(*,300) n, L2(1), L2(2), L2(3)
300     format(1X,i8,2(e15.6),3(e15.6),4(e15.6))
      end if

    end do iter_loop

  end subroutine ldc_explicit


  pure function create_residual(i, j, x_nodes, y_nodes, dx, dy, beta, soln)

    use set_precision, only : dp
    use setup,         only : rho, nu, visc_eps, c2

    implicit none

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

end module solvers
