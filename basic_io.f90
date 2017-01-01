subroutine tecplot

    use shared
    implicit none

    integer(4) :: i,j
    real(kind=8) :: p

    open(4,file='output.dat')


    write(4,*)'TITLE = "Shock Reflection "'
    write(4,*)'VARIABLES = "X","Y","rho","U","V","e","p"'
    write(4,*)'ZONE T ="Zone-one", I=',n_ipoints,',J=',n_jpoints,',F=POINT'

    do j = 1,n_jpoints
        do i = 1,n_ipoints

            p = (f_gamma-1.0d0) * (q_prim(4,i,j) - (0.5d0*(q_prim(2,i,j)**2.0d0) + q_prim(3,i,j)**2.0d0)/q_prim(1,i,j))

            write(4,*) x_coord(i,j),y_coord(i,j),q_prim(1,i,j),q_prim(2,i,j)/q_prim(1,i,j),q_prim(3,i,j)/q_prim(1,i,j),q_prim(4,i,j),p

        end do
    end do

    close(4)

end subroutine tecplot

subroutine print_iter(actual)

    use shared
    implicit none

    integer(kind=4) :: actual


    if (actual == 1) then
        pux = print_step
    end if
    
    if (actual == 1) then


    else if (actual == print_step) then

        write(*,'(A,I7,A,F10.4,A,F15.8,A)') " | Iteration: ", actual, " | Time/iter: ", time_per_iter, " [s] | Rho residue: ", log10(res_rho)," |"
        write(5,'(I7,F10.4)') actual,log10(res_rho)

        print_step = print_step + pux

    end if

end subroutine print_iter

subroutine out_midline

    use shared
    implicit none

    integer(kind=4) :: half, i
    real(kind=8) :: max_p = 0.0d0
    real(kind=8), dimension(n_ipoints) :: p

    open(7,file='mid_pressure.dat')

    half = n_jpoints/2

    do i = 1, n_ipoints
        p(i) = (f_gamma-1.0d0) * (q_prim(4,i,half) - (0.5d0*(q_prim(2,i,half)**2.0d0) + q_prim(3,i,half)**2.0d0)/q_prim(1,i,half))
    end do

    max_p = maxval(p)

    do i = 1, n_ipoints
        write(7,'(3F10.4)')  x_coord(i,half), p(i)
    end do

    close(7)

end subroutine out_midline

