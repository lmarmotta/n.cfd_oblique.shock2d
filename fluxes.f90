subroutine fluxes

    use shared
    implicit none

    integer(kind=4) :: i, j
    real(kind=8) :: u, v, e, rho, p


    ! Initial conditions of the fluxes.

    do j = 1, n_jpoints
        do i = 1, n_ipoints


            ! Initial physical malabarisms.

            rho = q_prim(1,i,j)
            e   = q_prim(4,i,j)
            u   = q_prim(2,i,j) / q_prim(1,i,j)
            v   = q_prim(3,i,j) / q_prim(1,i,j)
            p   = (f_gamma - 1.0d0) * ( e -  0.5d0*rho*(u**2.0d0 + v**2.0d0) )

            f(1,i,j) = rho * u
            f(2,i,j) = rho * u**2.0d0 + p
            f(3,i,j) = rho * u * v
            f(4,i,j) = ( e + p ) * u

            g(1,i,j) = rho * v
            g(2,i,j) = rho * u * v
            g(3,i,j) = rho * v**2.0d0 + p
            g(4,i,j) = ( e + p ) * v
            
        end do
    end do

end subroutine fluxes


