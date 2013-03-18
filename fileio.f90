module fileio

  implicit none

contains

!================================= read_input ================================80
!
! Reads the namelists within ldc.nml
!
!=============================================================================80

  subroutine read_input

    use setup,         only : x_nodes, y_nodes, length, xmin, xmax, ymin, ymax,&
                              re, rho, u_lid, p_guage, solver, max_iter, cfl,  &
                              k, c2, conv_toler, visc_eps
    use functions,     only : find_available_unit

    integer :: nml_unit

    namelist /domain/ x_nodes, y_nodes, length, xmin, xmax, ymin, ymax

    namelist /physical_properties/ re, rho, u_lid, p_guage

    namelist /solver_properties/ solver, max_iter, cfl, k, c2,                 &
                                 conv_toler, visc_eps

    continue

    nml_unit = find_available_unit()

    open(nml_unit, file='ldc.nml', status='old')

    rewind(nml_unit)
    read(nml_unit, nml=domain)

    rewind(nml_unit)
    read(nml_unit, nml=physical_properties)

    rewind(nml_unit)
    read(nml_unit, nml=solver_properties)

  end subroutine read_input

!============================= write_output ==================================80
!
! Writes the solution to a tecplot file readable by Tecplot or Paraview
!
!=============================================================================80

  subroutine write_output(x_nodes, y_nodes, soln)

    use set_precision, only : dp
    use functions,     only : find_available_unit

    integer,                                  intent(in) :: x_nodes, y_nodes
    real(dp), dimension(3, x_nodes, y_nodes), intent(in) :: soln

    integer :: out_unit
    integer :: i, j, eq

    continue

    out_unit = find_available_unit()

    open(out_unit, FILE='ldc.dat', status='unknown')

    write(out_unit,*) 'TITLE="Lid Driven Cavity Solution"'
    write(out_unit,*) 'VARIABLES = "X", "Y", "Pressure", "U", "V"'
    write(out_unit,*) 'ZONE DATAPACKING=BLOCK, I=', x_nodes, ', J=', y_nodes
    write(out_unit,*)

! write the grid
    do j = 1, y_nodes
      do i = 1, x_nodes
        write(out_unit,*) (real(i-1,dp)/real(x_nodes-1))
      end do
    end do

    do j = 1, y_nodes
      do i = 1, x_nodes
        write(out_unit,*) (real(j-1,dp)/real(y_nodes-1))
      end do
    end do

! write the solution
    do eq = 1,3
      do j = 1, y_nodes
        do i = 1, x_nodes
          write(out_unit,*) soln(eq, i, j)
        end do
      end do
    end do

  close(out_unit)

  end subroutine write_output

end module fileio
