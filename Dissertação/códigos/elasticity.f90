!===============================================================================
! MÓDULO DE ELASTICIDADE
!===============================================================================
module elasticity_mod
    use math_geometry_mod, only: dp
    implicit none
    private

    public :: compute_stress_2d, compute_von_mises_2d, get_lame_parameters
    public :: build_C_plane_stress, build_C_plane_strain
    public :: split_C_vol_dev, split_C_shear_normal
    public :: compute_vol_dev_energies_2d

contains

    pure subroutine compute_stress_2d(C_mat, eps_elem, sigma_elem)
        real(dp), intent(in)  :: C_mat(3,3), eps_elem(3)
        real(dp), intent(out) :: sigma_elem(3)
        sigma_elem = matmul(C_mat, eps_elem)
    end subroutine

    pure function compute_von_mises_2d(sigma_elem) result(von_mises)
        real(dp), intent(in) :: sigma_elem(3)
        real(dp) :: von_mises
        von_mises = sqrt(sigma_elem(1)**2 - sigma_elem(1)*sigma_elem(2) + &
                         sigma_elem(2)**2 + 3.0_dp*sigma_elem(3)**2)
    end function

    pure subroutine get_lame_parameters(E, nu, lambda, mu)
        real(dp), intent(in)  :: E, nu
        real(dp), intent(out) :: lambda, mu
        mu     = E / (2.0_dp*(1.0_dp + nu))
        lambda = (E*nu) / ((1.0_dp + nu)*(1.0_dp - 2.0_dp*nu))
    end subroutine

    pure subroutine build_C_plane_stress(E_mod, nu_val, C_mat)
        real(dp), intent(in)  :: E_mod, nu_val
        real(dp), intent(out) :: C_mat(3,3)
        C_mat = 0.0_dp
        C_mat(1,1) = E_mod/(1.0_dp - nu_val**2)
        C_mat(2,2) = C_mat(1,1)
        C_mat(1,2) = C_mat(1,1)*nu_val
        C_mat(2,1) = C_mat(1,2)
        C_mat(3,3) = E_mod/(2.0_dp*(1.0_dp + nu_val))
    end subroutine

    pure subroutine build_C_plane_strain(E_mod, nu_val, C_mat)
        real(dp), intent(in)  :: E_mod, nu_val
        real(dp), intent(out) :: C_mat(3,3)
        real(dp) :: lambda, mu
        mu     = E_mod/(2.0_dp*(1.0_dp + nu_val))
        lambda = (E_mod*nu_val)/((1.0_dp + nu_val)*(1.0_dp - 2.0_dp*nu_val))
        C_mat = 0.0_dp
        C_mat(1,1) = lambda + 2.0_dp*mu
        C_mat(2,2) = C_mat(1,1)
        C_mat(1,2) = lambda
        C_mat(2,1) = lambda
        C_mat(3,3) = mu
    end subroutine

    !--- Volumétrico / desviador (locking volumétrico) -------------------------
    pure subroutine split_C_vol_dev(E, nu, is_plane_strain, C_vol, C_dev)
        real(dp), intent(in)  :: E, nu
        logical,  intent(in)  :: is_plane_strain
        real(dp), intent(out) :: C_vol(3,3), C_dev(3,3)
        real(dp) :: lambda, mu, kappa, const
        call get_lame_parameters(E, nu, lambda, mu)
        C_vol = 0.0_dp; C_dev = 0.0_dp
        if (is_plane_strain) then
            kappa = lambda + mu
            C_vol(1,1) = kappa; C_vol(1,2) = kappa
            C_vol(2,1) = kappa; C_vol(2,2) = kappa
            C_dev(1,1) = (lambda + 2.0_dp*mu) - kappa
            C_dev(2,2) = C_dev(1,1)
            C_dev(1,2) = lambda - kappa
            C_dev(2,1) = C_dev(1,2)
            C_dev(3,3) = mu
        else
            const = E / (2.0_dp*(1.0_dp - nu))
            C_vol(1,1) = const; C_vol(1,2) = const
            C_vol(2,1) = const; C_vol(2,2) = const
            C_dev(1,1) = E/(1.0_dp - nu**2) - const
            C_dev(2,2) = C_dev(1,1)
            C_dev(1,2) = (E*nu)/(1.0_dp - nu**2) - const
            C_dev(2,1) = C_dev(1,2)
            C_dev(3,3) = E/(2.0_dp*(1.0_dp + nu))
        end if
    end subroutine

    !--- Cisalhamento / normal (shear locking) --------------------------------
    ! C_shear: apenas a parcela de engenharia γ_xy (entrada 3,3)
    ! C_normal: bloco (εx, εy)
    pure subroutine split_C_shear_normal(E, nu, is_plane_strain, C_shear, C_normal)
        real(dp), intent(in)  :: E, nu
        logical,  intent(in)  :: is_plane_strain
        real(dp), intent(out) :: C_shear(3,3), C_normal(3,3)
        real(dp) :: C_full(3,3)
        if (is_plane_strain) then
            call build_C_plane_strain(E, nu, C_full)
        else
            call build_C_plane_stress(E, nu, C_full)
        end if
        C_shear  = 0.0_dp
        C_normal = 0.0_dp
        C_shear(3,3) = C_full(3,3)
        C_normal(1:2,1:2) = C_full(1:2,1:2)
    end subroutine

    pure subroutine compute_vol_dev_energies_2d(eps, area, C_vol, C_dev, U_vol, U_dev)
        real(dp), intent(in)  :: eps(3), area, C_vol(3,3), C_dev(3,3)
        real(dp), intent(out) :: U_vol, U_dev
        real(dp) :: sv(3), sd(3)
        sv = matmul(C_vol, eps); sd = matmul(C_dev, eps)
        U_vol = 0.5_dp*(sv(1)*eps(1)+sv(2)*eps(2)+sv(3)*eps(3))*area
        U_dev = 0.5_dp*(sd(1)*eps(1)+sd(2)*eps(2)+sd(3)*eps(3))*area
    end subroutine

end module elasticity_mod
