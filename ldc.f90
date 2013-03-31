!================================== ldc ======================================80
!
! Solves the lid driven cavity problem for square, cartesian domains
!
!=============================================================================80

program ldc

  use fileio,  only : read_input, write_output
  use setup,   only : ldc_allocate, ldc_deallocate, solver,                    &
                      x_nodes, y_nodes, dx, dy, soln, soln_new
  use solvers, only : ldc_explicit, ldc_implicit, ldc_sgs

  implicit none

  continue

  write(*,*) 'Reading input file LDC.nml...'
  call read_input()

  write(*,*) 'Allocating and initializing arrays...'
  call ldc_allocate()

  select case(solver)
  case('explicit')
    write(*,*) 'Beginning explict solve...'

    call ldc_explicit(x_nodes, y_nodes, dx, dy, soln, soln_new)

  case('implicit')
    write(*,*) 'Beginning implict solve...'

    call ldc_implicit(x_nodes, y_nodes, dx, dy, soln, soln_new)

  case('sgs')
    write(*,*) 'Beginning SGS solve...'

    call ldc_sgs(x_nodes, y_nodes, dx, dy, soln)

  case default
    write(*,*) 'Unknown solver type... stopping'
    call ldc_deallocate
    stop
  end select

  call write_output(x_nodes, y_nodes, soln)

  call ldc_deallocate()

end program ldc

