subroutine indat

    use shared
    implicit none


    ! Read namelist parameters declared in the module.
    
    namelist /PAR_inlet/ mach,u_inf,v_inf,rho_inf,f_gamma
    namelist /PAR_numerical/ scheme,n_iter,dt,print_step
    namelist /PAR_geometrical/ n_ipoints,n_jpoints,x_max,y_max


    open(1,file='config.in')

    read(1,PAR_inlet)
    read(1,PAR_numerical)
    read(1,PAR_geometrical)

    close(1)

end subroutine indat

subroutine do_mesh

    use shared
    implicit none

    real(kind=8) :: dx, dy
    integer(kind=4) :: i, j

    real(kind=8), dimension(n_ipoints) :: x
    real(kind=8), dimension(n_jpoints) :: y

    write(*,'(A)') 
    write(*,'(A)') "---------------------------------------------------------------"
    write(*,'(A,I7,A,I7)') " + Generating mesh: n_ipoints=",n_ipoints," n_jpoints=",n_jpoints
    write(*,'(A)') "---------------------------------------------------------------"

    dx = x_max/real(n_ipoints, 8)
    dy = y_max/real(n_jpoints, 8)

    do i = 1, n_ipoints
        x(i) = i*dx
    end do

    do j = 1, n_jpoints
        y(j) = j*dy
    end do

    do i = 1, n_ipoints
        do j = 1, n_jpoints
            x_coord(i,j) = x(i)
            y_coord(i,j) = y(j)
        end do
    end do


end subroutine do_mesh

subroutine initial

    use shared
    implicit none

    integer(kind=4) :: i, j
    real(kind=8) :: u, v, p, rho

    write(*,'(A)') 
    write(*,'(A)') "---------------------------------------------------------------"
    write(*,'(A,F10.4)') " + Generating initial condition for Mach:", mach
    write(*,'(A)') "---------------------------------------------------------------"


    ! Initial condition of the primitive variables.

    do i = 1, n_ipoints
        do j = 1, n_jpoints

            u   = mach
            v   = 0.0d0
            p   = 1.0d0/f_gamma
            rho = 1.0d0

            q_prim(1,i,j) = rho
            q_prim(2,i,j) = rho*u
            q_prim(3,i,j) = rho*v
            q_prim(4,i,j) = p/(f_gamma-1.0d0) + 0.5d0*rho*(u**2.0d0 + v**2.0d0)

        end do
    end do


end subroutine initial
