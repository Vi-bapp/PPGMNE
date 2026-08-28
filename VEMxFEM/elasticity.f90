!===============================================================================
! MÓDULO DE ELASTICIDADE (Agnóstico ao Método Numérico)
!===============================================================================
module elasticity_mod
    use math_geometry_mod, only: dp
    implicit none
    private

    public :: compute_stress_2d
    public :: compute_von_mises_2d
    public :: get_lame_parameters
    public :: build_C_plane_stress, build_C_plane_strain

    contains

    ! Aplica a Lei de Hooke generalizada
    pure subroutine compute_stress_2d(C_mat, eps_elem, sigma_elem)
        real(dp), intent(in)  :: C_mat(3,3)
        real(dp), intent(in)  :: eps_elem(3)
        real(dp), intent(out) :: sigma_elem(3)

        sigma_elem = matmul(C_mat, eps_elem)
    end subroutine compute_stress_2d

    ! Calcula a Tensão Equivalente de Von Mises (Estado Plano)
    pure function compute_von_mises_2d(sigma_elem) result(von_mises)
        real(dp), intent(in) :: sigma_elem(3)
        real(dp) :: von_mises

        ! sigma_vM = sqrt(sigma_xx^2 - sigma_xx*sigma_yy + sigma_yy^2 + 3*tau_xy^2)
        von_mises = sqrt(sigma_elem(1)**2 - sigma_elem(1)*sigma_elem(2) + &
                         sigma_elem(2)**2 + 3.0_dp * sigma_elem(3)**2)
    end function compute_von_mises_2d

    pure subroutine get_lame_parameters(E, nu, lambda, mu)
        real(dp), intent(in)  :: E, nu
        real(dp), intent(out) :: lambda, mu

        mu = E / (2.0_dp * (1.0_dp + nu))
        lambda = (E * nu) / ((1.0_dp + nu) * (1.0_dp - 2.0_dp * nu))
    end subroutine get_lame_parameters

    ! Constrói a matriz para Estado Plano de Tensões
    pure subroutine build_C_plane_stress(E, nu, C_mat)
        real(dp), intent(in)  :: E, nu
        real(dp), intent(out) :: C_mat(3,3)
        real(dp) :: const

        C_mat = 0.0_dp

        C_mat(1,1) = 1.0_dp
        C_mat(2,2) = 1.0_dp
        C_mat(1,2) = nu
        C_mat(2,1) = nu
        C_mat(3,3) = E / (2.0_dp * (1.0_dp + nu)) ! Este é o próprio módulo de cisalhamento (mu)

        C_mat = C_mat * E / (1.0_dp - nu**2)
    end subroutine build_C_plane_stress

    ! Constrói a matriz para Estado Plano de Tensões
    pure subroutine build_C_plane_strain(E, nu, C_mat)
        real(dp), intent(in)  :: E, nu
        real(dp), intent(out) :: C_mat(3,3)
        real(dp) :: const

        C_mat = 0.0_dp

        ! Estado Plano de Deformação
        C_mat(1,1) = 1.0_dp - nu
        C_mat(2,2) = 1.0_dp - nu
        C_mat(1,2) = nu
        C_mat(2,1) = nu
        C_mat(3,3) = (1.0_dp - 2.0_dp * nu) / 2.0_dp

        C_mat = C_mat * (E / ((1.0_dp + nu) * (1.0_dp - 2.0_dp * nu)))
    end subroutine build_C_plane_strain

end module elasticity_mod