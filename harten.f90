subroutine roe_x

    use shared
    implicit none

    integer(kind=4) :: i, j
    
    real(kind=8) :: q1_l, q2_l, q3_l, q4_l
    real(kind=8) :: q1_r, q2_r, q3_r, q4_r
    real(kind=8) :: rho_l, u_l, v_l, h_l
    real(kind=8) :: rho_r, u_r, v_r, h_r

    do j = 1, n_jpoints
        do i = 1, n_ipoints - 1

            q1_l = q_prim(1,i,j)
            q2_l = q_prim(2,i,j)
            q3_l = q_prim(3,i,j)
            q4_l = q_prim(4,i,j)

            q1_r = q_prim(1,i+1,j)
            q2_r = q_prim(2,i+1,j)
            q3_r = q_prim(3,i+1,j)
            q4_r = q_prim(4,i+1,j)


            rho_l = q1_l
            rho_r = q1_r

            u_l = q2_l/q1_l
            u_r = q2_r/q1_r

            v_l = q3_l/q1_l
            v_r = q3_r/q1_r

            h_l = ( q4_l + (f_gamma-1.0d0) * (q4_l - 0.5d0 * ( q2_l**2.0d0 + q3_l**2.00d0 ) / q1_l ) ) / q1_l
            h_r = ( q4_r + (f_gamma-1.0d0) * (q4_r - 0.5d0 * ( q2_r**2.0d0 + q3_r**2.00d0 ) / q1_r ) ) / q1_r


            roe_u(i,j) = ( dsqrt(rho_l) * u_l + dsqrt(rho_r) * u_r) / ( dsqrt(rho_l) + dsqrt(rho_r) )
            roe_v(i,j) = ( dsqrt(rho_l) * v_l + dsqrt(rho_r) * v_r) / ( dsqrt(rho_l) + dsqrt(rho_r) )
            roe_h(i,j) = ( dsqrt(rho_l) * h_l + dsqrt(rho_r) * h_r) / ( dsqrt(rho_l) + dsqrt(rho_r) )
            roe_a(i,j) = dsqrt( (f_gamma-1.0d0)*(roe_h(i,j) - 0.5d0*(roe_u(i,j)**2.0d0 + roe_v(i,j)**2.0d0)) )

        end do
    end do

end subroutine roe_x

subroutine roe_y

    use shared
    implicit none

    integer(kind=4) :: i, j
    
    real(kind=8) :: q1_l, q2_l, q3_l, q4_l
    real(kind=8) :: q1_r, q2_r, q3_r, q4_r
    real(kind=8) :: rho_l, u_l, v_l, h_l
    real(kind=8) :: rho_r, u_r, v_r, h_r

    do j = 1, n_jpoints - 1
        do i = 1, n_ipoints

            q1_l = q_prim(1,i,j)
            q2_l = q_prim(2,i,j)
            q3_l = q_prim(3,i,j)
            q4_l = q_prim(4,i,j)

            q1_r = q_prim(1,i,j+1)
            q2_r = q_prim(2,i,j+1)
            q3_r = q_prim(3,i,j+1)
            q4_r = q_prim(4,i,j+1)


            rho_l = q1_l
            rho_r = q1_r

            u_l = q2_l/q1_l
            u_r = q2_r/q1_r

            v_l = q3_l/q1_l
            v_r = q3_r/q1_r

            h_l = ( q4_l + (f_gamma-1.0d0) * (q4_l - 0.5d0 * ( q2_l**2.0d0 + q3_l**2.00d0 ) / q1_l ) ) / q1_l
            h_r = ( q4_r + (f_gamma-1.0d0) * (q4_r - 0.5d0 * ( q2_r**2.0d0 + q3_r**2.00d0 ) / q1_r ) ) / q1_r


            roe_u(i,j) = ( dsqrt(rho_l) * u_l + dsqrt(rho_r) * u_r) / ( dsqrt(rho_l) + dsqrt(rho_r) )
            roe_v(i,j) = ( dsqrt(rho_l) * v_l + dsqrt(rho_r) * v_r) / ( dsqrt(rho_l) + dsqrt(rho_r) )
            roe_h(i,j) = ( dsqrt(rho_l) * h_l + dsqrt(rho_r) * h_r) / ( dsqrt(rho_l) + dsqrt(rho_r) )
            roe_a(i,j) = dsqrt( (f_gamma-1.0d0)*(roe_h(i,j) - 0.5d0*(roe_u(i,j)**2.0d0 + roe_v(i,j)**2.0d0)) )

        end do
    end do

end subroutine roe_y

subroutine mod_fluxes_x

    use shared
    implicit none

    integer(kind=4) :: i, j

    real(kind=8) :: ax_1, ax_2, ax_3, ax_4
    real(kind=8) :: u, v, h, a
    real(kind=8) :: d_rho, d_rhou, d_rhov, d_e
    real(kind=8) :: aa, bb, cc
    real(kind=8) :: eps1, eps2, eps3, eps4


    real(kind=8), dimension(4) :: r1
    real(kind=8), dimension(4) :: r2
    real(kind=8), dimension(4) :: r3
    real(kind=8), dimension(4) :: r4

    real(kind=8), dimension(4) :: a_k

    real(kind=8), dimension(4) :: psi

    eps1 = 0.01d0
    eps2 = 0.0d0
    eps3 = 0.0d0
    eps4 = 0.0d0
    


    do i = 1, n_ipoints - 1
        do j = 1, n_jpoints

            u = roe_u(i,j)
            v = roe_v(i,j)
            h = roe_h(i,j)
            a = roe_a(i,j)

            ax_1 = u - a
            ax_2 = u
            ax_3 = u + a
            ax_4 = u


            r1(1) = 1.0d0
            r1(2) = u - a
            r1(3) = v
            r1(4) = h - (u*a)

            r2(1) = 1.0d0
            r2(2) = u
            r2(3) = v
            r2(4) = (u**2.0d0 + v**2.0d0)/2.0d0

            r3(1) = 1.0d0
            r3(2) = u + a
            r3(3) = v
            r3(4) = h + (u*a)

            r4(1) = 0.0d0
            r4(2) = 0.0d0
            r4(3) = 1.0d0
            r4(4) = v


            d_rho  = q_prim(1,i+1,j) - q_prim(1,i,j)
            d_rhou = q_prim(2,i+1,j) - q_prim(2,i,j)
            d_rhov = q_prim(3,i+1,j) - q_prim(3,i,j)
            d_e    = q_prim(4,i+1,j) - q_prim(4,i,j)


            aa = (f_gamma - 1.0d0) * (d_e + 0.5d0*(u**2.0d0 + v**2.0d0) * d_rho - u*d_rhou - v*d_rhov) / a**2.0d0
            bb = (d_rhou - u*d_rho)/a
            cc =  d_rhov - v*d_rho


            a_k(1) = (aa - bb) / 2.0d0
            a_k(2) =  d_rho - aa
            a_k(3) = (aa + bb) / 2.0d0
            a_k(4) =  cc

            psi(1) = q_z(ax_1,eps1)
            psi(2) = q_z(ax_2,eps2)
            psi(3) = q_z(ax_3,eps3)
            psi(4) = q_z(ax_4,eps4)

            f_harten(1,i,j) = 0.5d0*(f(1,i,j) + f(1,i+1,j)) - 0.5d0*( psi(1)*a_k(1)*r1(1) + psi(2)*a_k(2)*r2(1) + psi(3)*a_k(3)*r3(1) + psi(4)*a_k(4)*r4(1))
            f_harten(2,i,j) = 0.5d0*(f(2,i,j) + f(2,i+1,j)) - 0.5d0*( psi(1)*a_k(1)*r1(2) + psi(2)*a_k(2)*r2(2) + psi(3)*a_k(3)*r3(2) + psi(4)*a_k(4)*r4(2))
            f_harten(3,i,j) = 0.5d0*(f(3,i,j) + f(3,i+1,j)) - 0.5d0*( psi(1)*a_k(1)*r1(3) + psi(2)*a_k(2)*r2(3) + psi(3)*a_k(3)*r3(3) + psi(4)*a_k(4)*r4(3))
            f_harten(4,i,j) = 0.5d0*(f(4,i,j) + f(4,i+1,j)) - 0.5d0*( psi(1)*a_k(1)*r1(4) + psi(2)*a_k(2)*r2(4) + psi(3)*a_k(3)*r3(4) + psi(4)*a_k(4)*r4(4))

        end do
    end do

end subroutine mod_fluxes_x

subroutine mod_fluxes_y

    use shared
    implicit none

    integer(kind=4) :: i, j

    real(kind=8) :: ax_1, ax_2, ax_3, ax_4
    real(kind=8) :: u, v, h, a
    real(kind=8) :: d_rho, d_rhou, d_rhov, d_e
    real(kind=8) :: aa, bb, cc
    real(kind=8) :: eps1, eps2, eps3, eps4


    real(kind=8), dimension(4) :: r1
    real(kind=8), dimension(4) :: r2
    real(kind=8), dimension(4) :: r3
    real(kind=8), dimension(4) :: r4

    real(kind=8), dimension(4) :: a_k

    real(kind=8), dimension(4) :: psi

    eps1 = 0.01d0
    eps2 = 0.0d0
    eps3 = 0.0d0
    eps4 = 0.0d0
    


    do i = 1, n_ipoints
        do j = 1, n_jpoints - 1

            u = roe_u(i,j)
            v = roe_v(i,j)
            h = roe_h(i,j)
            a = roe_a(i,j)

            ax_1 = v - a
            ax_2 = v
            ax_3 = v + a
            ax_4 = v


            r1(1) = 1.0d0
            r1(2) = u 
            r1(3) = v - a
            r1(4) = h - v*a

            r2(1) = 1.0d0
            r2(2) = u
            r2(3) = v
            r2(4) = (u**2.0d0 + v**2.0d0)/2.0d0

            r3(1) = 1.0d0
            r3(2) = u 
            r3(3) = v + a
            r3(4) = h + v*a

            r4(1) = 0.0d0
            r4(2) = 1.0d0
            r4(3) = 0.0d0
            r4(4) = u


            d_rho  = q_prim(1,i,j+1) - q_prim(1,i,j)
            d_rhou = q_prim(2,i,j+1) - q_prim(2,i,j)
            d_rhov = q_prim(3,i,j+1) - q_prim(3,i,j)
            d_e    = q_prim(4,i,j+1) - q_prim(4,i,j)

            aa = (f_gamma - 1.0d0) * (d_e + (u**2.0d0 + v**2.0d0)/2.0d0 * d_rho - u*d_rhou - v*d_rhov) / a**2.0d0
            bb = (d_rhov - v*d_rho)/a
            cc = d_rhou - u*d_rho


            a_k(1) = (aa - bb) / 2.0d0
            a_k(2) = d_rho - aa
            a_k(3) = (aa + bb) / 2.0d0
            a_k(4) = cc
            

            psi(1) = q_z(ax_1,eps1)
            psi(2) = q_z(ax_2,eps2)
            psi(3) = q_z(ax_3,eps3)
            psi(4) = q_z(ax_4,eps4)

            g_harten(1,i,j) = 0.5d0*(g(1,i,j) + g(1,i,j+1)) - 0.5d0*( psi(1)*a_k(1)*r1(1) + psi(2)*a_k(2)*r2(1) + psi(3)*a_k(3)*r3(1) + psi(4)*a_k(4)*r4(1))
            g_harten(2,i,j) = 0.5d0*(g(2,i,j) + g(2,i,j+1)) - 0.5d0*( psi(1)*a_k(1)*r1(2) + psi(2)*a_k(2)*r2(2) + psi(3)*a_k(3)*r3(2) + psi(4)*a_k(4)*r4(2))
            g_harten(3,i,j) = 0.5d0*(g(3,i,j) + g(3,i,j+1)) - 0.5d0*( psi(1)*a_k(1)*r1(3) + psi(2)*a_k(2)*r2(3) + psi(3)*a_k(3)*r3(3) + psi(4)*a_k(4)*r4(3))
            g_harten(4,i,j) = 0.5d0*(g(4,i,j) + g(4,i,j+1)) - 0.5d0*( psi(1)*a_k(1)*r1(4) + psi(2)*a_k(2)*r2(4) + psi(3)*a_k(3)*r3(4) + psi(4)*a_k(4)*r4(4))

        end do
    end do

end subroutine mod_fluxes_y

