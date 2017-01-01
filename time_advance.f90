subroutine residue

    use shared
    implicit none

    integer(kind=4) :: i, j
    real(kind=8) :: dx, dy


    max_res = 0.0d0

    do j = 2, n_jpoints - 1
        do i = 2, n_ipoints - 1

            dx = x_coord(i+1,j) - x_coord(i,j)
            dy = y_coord(i,j+1) - y_coord(i,j)

            res(1,i,j) = -(1.0d0/dx)*(f_harten(1,i,j) - f_harten(1,i-1,j)) - (1.0d0/dy)*(g_harten(1,i,j) - g_harten(1,i,j-1))
            res(2,i,j) = -(1.0d0/dx)*(f_harten(2,i,j) - f_harten(2,i-1,j)) - (1.0d0/dy)*(g_harten(2,i,j) - g_harten(2,i,j-1))
            res(3,i,j) = -(1.0d0/dx)*(f_harten(3,i,j) - f_harten(3,i-1,j)) - (1.0d0/dy)*(g_harten(3,i,j) - g_harten(3,i,j-1))
            res(4,i,j) = -(1.0d0/dx)*(f_harten(4,i,j) - f_harten(4,i-1,j)) - (1.0d0/dy)*(g_harten(4,i,j) - g_harten(4,i,j-1))


            ! Get the max residue.

            if (dabs(res(1,i,j)) > max_res) then 
                max_res = res(1,i,j)
            end if

            res_rho = max_res

        end do
    end do


end subroutine residue

subroutine time_step

    use shared
    implicit none

    integer(kind=4) :: i, j


    do i = 2, n_ipoints - 1
        do j = 2, n_jpoints - 1

            q_prim(1,i,j) = q_prim(1,i,j) + dt*res(1,i,j)
            q_prim(2,i,j) = q_prim(2,i,j) + dt*res(2,i,j)
            q_prim(3,i,j) = q_prim(3,i,j) + dt*res(3,i,j)
            q_prim(4,i,j) = q_prim(4,i,j) + dt*res(4,i,j)

        end do
    end do

    res = 0.0d0


end subroutine time_step

