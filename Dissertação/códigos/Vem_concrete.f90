!===============================================================================
! MÓDULO DE OPERADORES CONCRETOS DE BORDO E DOMÍNIO
!===============================================================================
module vem_concrete_operators_mod
    use math_geometry_mod
    use vem_operators_mod
    implicit none
    private

    public :: op_identity, op_grad, op_lap, op_eps_l2_boundary, op_eps_boundary
    public :: op_coeff_lap, op_coeff_div, op_coeff_grad, eval_monomial_derivatives

contains

    subroutine op_identity(x_pt, y_pt, xc, yc, h_E, p_exp, q_exp, normal, dof_dir, op_val)
        real(dp), intent(in) :: x_pt, y_pt, xc, yc, h_E
        integer, intent(in)  :: p_exp, q_exp, dof_dir
        real(dp), intent(in) :: normal(2)
        real(dp), intent(out):: op_val(:)
        real(dp) :: m_val

        call eval_monomial(x_pt, y_pt, xc, yc, h_E, p_exp, q_exp, m_val)
        op_val = 0.0_dp
        if (dof_dir >= 1 .and. dof_dir <= size(op_val)) then
            op_val(dof_dir) = m_val
        end if
    end subroutine op_identity

    subroutine op_grad(x_pt, y_pt, xc, yc, h_E, p_exp, q_exp, normal, dof_dir, op_val)
        real(dp), intent(in) :: x_pt, y_pt, xc, yc, h_E
        integer, intent(in)  :: p_exp, q_exp, dof_dir
        real(dp), intent(in) :: normal(2)
        real(dp), intent(out):: op_val(:)
        real(dp) :: g(2)

        call eval_monomial_grad(x_pt, y_pt, xc, yc, h_E, p_exp, q_exp, g)
        op_val = 0.0_dp
        if (size(op_val) >= 2) op_val(1:2) = g
    end subroutine op_grad

    subroutine op_lap(x_pt, y_pt, xc, yc, h_E, p_exp, q_exp, normal, dof_dir, op_val)
        real(dp), intent(in) :: x_pt, y_pt, xc, yc, h_E
        integer, intent(in)  :: p_exp, q_exp, dof_dir
        real(dp), intent(in) :: normal(2)
        real(dp), intent(out):: op_val(:)
        real(dp) :: l_val

        call eval_monomial_lap(x_pt, y_pt, xc, yc, h_E, p_exp, q_exp, l_val)
        op_val = 0.0_dp
        if (size(op_val) >= 1) op_val(1) = l_val
    end subroutine op_lap

    subroutine op_eps_l2_boundary(x_pt, y_pt, xc, yc, h_E, p_exp, q_exp, normal, dof_dir, op_val)
        real(dp), intent(in) :: x_pt, y_pt, xc, yc, h_E
        integer, intent(in)  :: p_exp, q_exp, dof_dir
        real(dp), intent(in) :: normal(2)
        real(dp), intent(out):: op_val(:)
        real(dp) :: m_val, nx, ny

        nx = normal(1); ny = normal(2)
        call eval_monomial(x_pt, y_pt, xc, yc, h_E, p_exp, q_exp, m_val)

        op_val = 0.0_dp
        if (size(op_val) >= 3) then
            select case (dof_dir)
            case (1)
                op_val(1) = m_val * nx
                op_val(3) = m_val * ny
            case (2)
                op_val(2) = m_val * ny
                op_val(3) = m_val * nx
            end select
        end if
    end subroutine op_eps_l2_boundary

    subroutine op_eps_boundary(x_pt, y_pt, xc, yc, h_E, p_exp, q_exp, normal, dof_dir, op_val)
        real(dp), intent(in) :: x_pt, y_pt, xc, yc, h_E
        integer, intent(in)  :: p_exp, q_exp, dof_dir
        real(dp), intent(in) :: normal(2)
        real(dp), intent(out):: op_val(:)
        real(dp) :: nx, ny, m_val 

        nx = normal(1); ny = normal(2)
        op_val = 0.0_dp

        call eval_monomial(x_pt, y_pt, xc, yc, h_E, p_exp, q_exp, m_val)

        if (size(op_val) >= 3) then
            select case (dof_dir)
            case (1)
                op_val(1) = m_val * nx
                op_val(2) = 0.0_dp
                op_val(3) = m_val * ny
            case (2)
                op_val(1) = 0.0_dp
                op_val(2) = m_val * ny
                op_val(3) = m_val * nx
            end select
        end if
    end subroutine op_eps_boundary

    subroutine eval_monomial_derivatives(x, y, xc, yc, h_E, p, q, dx_m, dy_m)
        real(dp), intent(in) :: x, y, xc, yc, h_E
        integer, intent(in)  :: p, q
        real(dp), intent(out):: dx_m, dy_m

        real(dp) :: rx, ry

        rx = (x - xc) / h_E
        ry = (y - yc) / h_E

        dx_m = merge((real(p, dp) / h_E) * (rx**max(0, p - 1)) * (ry**q), 0.0_dp, p > 0)
        dy_m = merge((real(q, dp) / h_E) * (rx**p) * (ry**max(0, q - 1)), 0.0_dp, q > 0)
    end subroutine eval_monomial_derivatives

    subroutine op_coeff_lap(p_in, q_in, h_E, dof_dir, n_terms, p_out, q_out, coeffs)
        integer, intent(in)  :: p_in, q_in, dof_dir
        real(dp), intent(in) :: h_E
        integer, intent(out) :: n_terms
        integer, intent(out) :: p_out(:), q_out(:)
        real(dp), intent(out):: coeffs(:)
        integer :: p_temp(2), q_temp(2)
        real(dp) :: c_temp(2)

        call eval_lap_coeffs(p_in, q_in, h_E, n_terms, p_temp, q_temp, c_temp)
        p_out(1:n_terms)  = p_temp(1:n_terms)
        q_out(1:n_terms)  = q_temp(1:n_terms)
        coeffs(1:n_terms) = c_temp(1:n_terms)
    end subroutine op_coeff_lap

    subroutine op_coeff_div(p_in, q_in, h_E, dof_dir, n_terms, p_out, q_out, coeffs)
        integer, intent(in)  :: p_in, q_in, dof_dir
        real(dp), intent(in) :: h_E
        integer, intent(out) :: n_terms
        integer, intent(out) :: p_out(:), q_out(:)
        real(dp), intent(out):: coeffs(:)
        integer :: p_temp(1), q_temp(1)
        real(dp) :: c_temp(1)

        call eval_div_coeffs(p_in, q_in, h_E, dof_dir, n_terms, p_temp, q_temp, c_temp)
        p_out(1:n_terms)  = p_temp(1:n_terms)
        q_out(1:n_terms)  = q_temp(1:n_terms)
        coeffs(1:n_terms) = c_temp(1:n_terms)
    end subroutine op_coeff_div

    subroutine op_coeff_grad(p_in, q_in, h_E, dof_dir, n_terms, p_out, q_out, coeffs)
        integer, intent(in)  :: p_in, q_in, dof_dir
        real(dp), intent(in) :: h_E
        integer, intent(out) :: n_terms
        integer, intent(out) :: p_out(:), q_out(:)
        real(dp), intent(out):: coeffs(:)
        integer :: p_temp(1), q_temp(1)
        real(dp) :: c_temp(1)

        call eval_grad_coeffs(p_in, q_in, h_E, dof_dir, n_terms, p_temp, q_temp, c_temp)
        p_out(1:n_terms)  = p_temp(1:n_terms)
        q_out(1:n_terms)  = q_temp(1:n_terms)
        coeffs(1:n_terms) = c_temp(1:n_terms)
    end subroutine op_coeff_grad

end module vem_concrete_operators_mod