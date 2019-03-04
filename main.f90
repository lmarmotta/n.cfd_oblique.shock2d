program simple_shock

    use shared
    implicit none

    real(kind=8) :: rate
    integer(kind=4) :: crr, c11, c22
    integer(kind=4) :: t

    open(5,file='residue.dat')

    call indat


    ! Main coordnate arrays of the mesh.
    allocate(x_coord(n_ipoints,n_jpoints))
    allocate(y_coord(n_ipoints,n_jpoints))

    ! Primitive properties of the mesh.
    allocate(q_prim(4,n_ipoints,n_jpoints))

    ! Allocate the fluxes.
    allocate(f(4,n_ipoints,n_jpoints))
    allocate(g(4,n_ipoints,n_jpoints))

    allocate(f_harten(4,n_ipoints,n_jpoints))
    allocate(g_harten(4,n_ipoints,n_jpoints))

    ! Allocate the residue vector.
    allocate(res(4,n_ipoints,n_jpoints))

    ! Roe avgs.
    allocate(roe_u(n_ipoints,n_jpoints))
    allocate(roe_v(n_ipoints,n_jpoints))
    allocate(roe_h(n_ipoints,n_jpoints))
    allocate(roe_a(n_ipoints,n_jpoints))


    ! Lets prepare the initial mesh and initial conditions.

    call do_mesh

    call initial

    call boundary

    time_per_iter = 0.0d0

    call system_clock(count_rate=crr)

    rate = real(crr,kind=8)

    write(*,*) ""
    write(*,'(A)') "---------------------------------------------------------------"
    if (scheme == 1) then
        write(*,'(A)') " + Running Harten Scheme [1st Order] "
    end if
    write(*,'(A)') "---------------------------------------------------------------"
    write(*,*)


    do t = 1, n_iter

        if (scheme == 1) then

            call system_clock(c11)

            call fluxes          
            call roe_x
            call mod_fluxes_x
            call roe_y
            call mod_fluxes_y
            call residue
            call time_step
            call boundary

            call system_clock(c22)

            time_per_iter = real(c22 - c11,kind=8)/rate

        end if

        call print_iter(t)

    end do

    call tecplot
    call out_midline


    deallocate(x_coord)
    deallocate(y_coord)
    deallocate(q_prim)

    deallocate(res)

    deallocate(roe_u)
    deallocate(roe_v)
    deallocate(roe_h)
    deallocate(roe_a)

    deallocate(f)
    deallocate(f_harten)

    deallocate(g)
    deallocate(g_harten)

    close(5)


end program simple_shock
