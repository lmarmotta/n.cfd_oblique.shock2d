module shared
    implicit none


    ! 
    ! +++ NAMELIST DECLARATION +++
    ! 

    integer(kind=4) :: scheme

    real(kind=8) :: mach 
    real(kind=8) :: u_inf
    real(kind=8) :: v_inf
    real(kind=8) :: rho_inf
    real(kind=8) :: p_inf
    real(kind=8) :: f_gamma 

    integer(kind=4) :: n_iter 
    real(kind=8) :: dt

    integer(kind=4) :: print_step 

    ! Number of points in the mesh.
    integer(kind=4) :: n_ipoints 
    integer(kind=4) :: n_jpoints

    ! Mesh dimension.
    real(kind=8) :: x_max 
    real(kind=8) :: y_max

    ! Residues.
    real(kind=8) :: res_rho
    real(kind=8) :: max_res = 10.0e-12


    !
    ! +++ GENERAL CODE VARIABLES +++
    !

    ! Arrays to store mesh coordnates.
    real(kind=8), allocatable, dimension(:,:) :: x_coord
    real(kind=8), allocatable, dimension(:,:) :: y_coord

    ! Properties arrays.
    real(kind=8), allocatable, dimension(:,:,:) :: q_prim
    real(kind=8), allocatable, dimension(:,:,:) :: f
    real(kind=8), allocatable, dimension(:,:,:) :: g

    real(kind=8), allocatable, dimension(:,:,:) :: f_harten
    real(kind=8), allocatable, dimension(:,:,:) :: g_harten

    ! Roe avg's
    real(kind=8), allocatable, dimension(:,:) :: roe_u
    real(kind=8), allocatable, dimension(:,:) :: roe_v
    real(kind=8), allocatable, dimension(:,:) :: roe_h
    real(kind=8), allocatable, dimension(:,:) :: roe_a

    ! Residue vector.
    real(kind=8), allocatable, dimension(:,:,:) :: res

    ! Print reference variable.
    integer(kind=4) :: pux
    real(kind=8) :: time_per_iter


    contains 

        real(kind=8) function q_z(z,eps)

            implicit none

            real(kind=8) :: z, eps

            q_z = dabs(z)

            if (dabs(z) < eps) then
                q_z = 0.5d0*( (z**2.0d0/eps) + eps)
            end if

        end function q_z


end module shared


