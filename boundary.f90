subroutine boundary

    use shared
    implicit none

    integer(kind=4) :: i, j
    real(kind=8) :: u, v, rho, p

    do j = 1, n_jpoints

        i = 1

        u   = mach
        v   = 0.0d0
        p   = 1.0d0/f_gamma
        rho = 1.0d0

        q_prim(1,i,j) = rho
        q_prim(2,i,j) = rho*u
        q_prim(3,i,j) = rho*v
        q_prim(4,i,j) = p/(f_gamma-1.0d0) + 0.5d0*rho*(u**2.0d0 + v**2.0d0)

    end do

    do i = 1, n_ipoints

        j = 1

        v = 0.0d0
        u = q_prim(2,i,j+1)/q_prim(1,i,j+1)
        rho = q_prim(1,i,j+1)
        p = (f_gamma-1.0d0) * (q_prim(4,i,j+1) - 0.5d0*(q_prim(2,i,j+1)**2.0d0 + q_prim(3,i,j+1)**2.0d0)/rho)

        q_prim(1,i,j) = rho
        q_prim(2,i,j) = rho*u
        q_prim(3,i,j) = rho*v
        q_prim(4,i,j) = (p/(f_gamma-1.0d0)) + 0.5d0*rho*(u**2.0d0 + v**2.0d0)

    end do

    do j = 1, n_jpoints

        i = n_ipoints

        q_prim(1,i,j) = q_prim(1,i-1,j)
        q_prim(2,i,j) = q_prim(2,i-1,j)
        q_prim(3,i,j) = q_prim(3,i-1,j)
        q_prim(4,i,j) = q_prim(4,i-1,j)

    end do

    do i = 1, n_ipoints
        
        j = n_jpoints

        rho =  1.69996629d0
        u   =  2.3780719d0 * 0.981826459d0
        v   = -2.3780719d0 * 0.189780933d0
        p   =  2.13947107d0 / f_gamma


        q_prim(1,i,j) = rho
        q_prim(2,i,j) = rho*u
        q_prim(3,i,j) = rho*v
        q_prim(4,i,j) = (p/(f_gamma-1.0d0)) + 0.5d0*rho*(u**2.0d0 + v**2.0d0)

    end do

end subroutine boundary
