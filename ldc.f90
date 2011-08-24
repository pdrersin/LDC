!=============================================================================80
!
!
!
!=============================================================================80

program LDC

  use fileio, only : read_input, write_output
  use setup,  only : ldc_allocate, ldc_deallocate, solver, &
                     x_nodes, y_nodes, soln

  implicit none

  continue

  write(*,*) 'Reading input file LDC.nml...'
  call read_input()

  write(*,*) 'Allocating and initializing arrays...'
  call ldc_allocate()

  select case(solver)
  case('explicit')
    write(*,*) 'Beginning explict solve...'

!    call ldc_explicit

  case('implicit')
    write(*,*) 'Beginning implict solve...'

!    call ldc_implicit

  case default
    write(*,*) 'Unknown solver type... stopping'
    call ldc_deallocate
    stop
  end select

  call write_output(x_nodes, y_nodes, soln)

  call ldc_deallocate()

end program LDC

