!===============================================================================
! MÓDULO DE ÁLGEBRA LINEAR, GEOMETRIA E OPERADORES MONOMIAIS
!===============================================================================
module math_geometry_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: dp
    public :: combination, eye, inverse_matrix, polygon_geometry, polygon_moment
    public :: get_monomial_exponents, find_monomial_index
    public :: get_gauss_lobatto 
    public :: solve_linear_system
    public :: eval_monomial, eval_monomial_grad, eval_monomial_div, eval_monomial_lap, eval_monomial_curl
    public :: eval_lap_coeffs, eval_div_coeffs, eval_grad_coeffs
    public :: compute_monomial_domain_integral

contains

    pure function compute_monomial_domain_integral(p, q, xc, yc, h_E, polygon_coords, n_verts) result(res)
        integer, intent(in) :: p, q, n_verts
        real(dp), intent(in) :: xc, yc, h_E
        real(dp), intent(in) :: polygon_coords(n_verts, 2)
        real(dp) :: res
        real(dp) :: x_shift(n_verts), y_shift(n_verts), mom

        x_shift = polygon_coords(:, 1) - xc
        y_shift = polygon_coords(:, 2) - yc

        call polygon_moment(x_shift, y_shift, n_verts, p, q, mom)
        res = mom / (h_E**(p + q))
    end function compute_monomial_domain_integral

    pure subroutine get_gauss_lobatto(n, pts, w)
        integer, intent(in) :: n
        real(dp), intent(out) :: pts(n), w(n)
        
        select case(n)
            case(1)
                pts(1) = 0.5_dp; w(1) = 1.0_dp
            case(2)
                pts(1) = 0.0_dp; w(1) = 0.5_dp
                pts(2) = 1.0_dp; w(2) = 0.5_dp
            case(3)
                pts(1) = 0.0_dp;        w(1) = 1.0_dp / 6.0_dp
                pts(2) = 0.5_dp;        w(2) = 4.0_dp / 6.0_dp
                pts(3) = 1.0_dp;        w(3) = 1.0_dp / 6.0_dp
            case(4)
                pts(1) = 0.0_dp;        w(1) = 1.0_dp / 12.0_dp
                pts(2) = 0.276393202250021_dp; w(2) = 5.0_dp / 12.0_dp
                pts(3) = 0.723606797749979_dp; w(3) = 5.0_dp / 12.0_dp
                pts(4) = 1.0_dp;        w(4) = 1.0_dp / 12.0_dp
            case default
                pts = 0.0_dp; w = 0.0_dp
        end select
    end subroutine get_gauss_lobatto

    subroutine solve_linear_system(A, b, x, n)
        integer, intent(in) :: n
        real(dp), intent(inout) :: A(n,n)
        real(dp), intent(in) :: b(n)
        real(dp), intent(out) :: x(n)
        
        integer :: ipiv(n), info
        
        x = b
        call dgesv(n, 1, A, n, ipiv, x, n, info)
        if (info /= 0) error stop 'Erro na subrotina DGESV do LAPACK ao resolver sistema.'
    end subroutine solve_linear_system

    pure function combination(n, k) result(res)
        integer, intent(in) :: n, k
        real(dp) :: res
        integer :: i, m

        if (k < 0 .or. k > n) then
            res = 0.0_dp; return
        end if
        if (k == 0 .or. k == n) then
            res = 1.0_dp; return
        end if
        m = min(k, n - k)
        res = 1.0_dp
        do i = 1, m
            res = res * real(n - i + 1, dp) / real(i, dp)
        end do
    end function combination

    pure function eye(n) result(Imat)
        integer, intent(in) :: n
        real(dp) :: Imat(n, n)
        integer :: idx
        Imat = 0.0_dp
        do concurrent (idx = 1:n)
            Imat(idx, idx) = 1.0_dp
        end do
    end function eye

    function inverse_matrix(A) result(Ainv)
        real(dp), intent(in) :: A(:,:)
        real(dp) :: Ainv(size(A,1), size(A,2))
        real(dp), allocatable :: work(:)
        integer, allocatable :: ipiv(:)
        integer :: n, info, lwork

        n = size(A, 1)
        Ainv = A
        allocate(ipiv(n), work(1))

        call dgetrf(n, n, Ainv, n, ipiv, info)
        if (info == 0) then
            call dgetri(n, Ainv, n, ipiv, work, -1, info)
            lwork = nint(work(1))
            deallocate(work); allocate(work(lwork))
            call dgetri(n, Ainv, n, ipiv, work, lwork, info)
        else
            error stop 'Erro na fatoracao LU do LAPACK em inverse_matrix.'
        end if
    end function inverse_matrix

pure subroutine polygon_geometry(x, y, n, area, xc, yc, h_E)
        integer, intent(in) :: n
        real(dp), intent(in) :: x(n), y(n)
        real(dp), intent(out) :: area, xc, yc, h_E
        integer :: i, inext
        real(dp) :: cross_sum, x_sum, y_sum, max_dist_sq, dist_sq
        real(dp) :: xi, yi, xnext, ynext, cross

        cross_sum = 0.0_dp; x_sum = 0.0_dp; y_sum = 0.0_dp

        ! Laço sequencial padrão compatível com gfortran
        do i = 1, n
            inext = mod(i, n) + 1
            xi = x(i);         yi = y(i)
            xnext = x(inext); ynext = y(inext)

            cross = xi * ynext - xnext * yi
            cross_sum = cross_sum + cross
            x_sum     = x_sum + (xi + xnext) * cross
            y_sum     = y_sum + (yi + ynext) * cross
        end do

        area = 0.5_dp * abs(cross_sum)
        if (area > 0.0_dp) then
            xc = x_sum / (6.0_dp * area)
            yc = y_sum / (6.0_dp * area)
        else
            xc = 0.0_dp; yc = 0.0_dp
        end if

        max_dist_sq = 0.0_dp
        do i = 1, n - 1
            dist_sq = maxval((x(i+1:n) - x(i))**2 + (y(i+1:n) - y(i))**2)
            if (dist_sq > max_dist_sq) max_dist_sq = dist_sq
        end do
        h_E = sqrt(max_dist_sq)
    end subroutine polygon_geometry

    pure subroutine polygon_moment(x, y, n, p, q, nu)
        integer, intent(in) :: n, p, q
        real(dp), intent(in) :: x(n), y(n)
        real(dp), intent(out) :: nu
        real(dp) :: coef_denom, poly_sum
        integer :: i, inext, k, l
        real(dp) :: xi, yi, xnext, ynext, cross_i, inner_sum

        coef_denom = real(p + q + 2, dp) * real(p + q + 1, dp) * combination(p + q, p)
        poly_sum = 0.0_dp

        ! Laço sequencial padrão compatível com gfortran
        do i = 1, n
            inext = mod(i, n) + 1
            xi = x(i);         yi = y(i)
            xnext = x(inext); ynext = y(inext)
            cross_i = xi * ynext - xnext * yi

            inner_sum = 0.0_dp
            do k = 0, p
                do l = 0, q
                    inner_sum = inner_sum + combination(k + l, l) * &
                                            combination((p + q) - (k + l), q - l) * &
                                            (xnext**k) * (xi**(p - k)) * &
                                            (ynext**l) * (yi**(q - l))
                end do
            end do

            poly_sum = poly_sum + cross_i * inner_sum
        end do

        nu = poly_sum / coef_denom
    end subroutine polygon_moment

    pure subroutine get_monomial_exponents(deg_k, n_mon, p_exp, q_exp)
        integer, intent(in)  :: deg_k, n_mon
        integer, intent(out) :: p_exp(n_mon), q_exp(n_mon)
        integer :: deg, p, q, idx

        idx = 0
        do deg = 0, deg_k
            do p = deg, 0, -1
                q = deg - p
                idx = idx + 1
                if (idx <= n_mon) then
                    p_exp(idx) = p
                    q_exp(idx) = q
                end if
            end do
        end do
    end subroutine get_monomial_exponents

    pure function find_monomial_index(p_target, q_target, p_arr, q_arr, n_arr) result(idx)
        integer, intent(in) :: p_target, q_target, n_arr
        integer, intent(in) :: p_arr(n_arr), q_arr(n_arr)
        integer :: idx, i

        idx = 0
        do i = 1, n_arr
            if (p_arr(i) == p_target .and. q_arr(i) == q_target) then
                idx = i
                exit
            end if
        end do
    end function find_monomial_index

    pure subroutine eval_monomial(x, y, xc, yc, h_E, p, q, val)
        real(dp), intent(in) :: x, y, xc, yc, h_E
        integer, intent(in)  :: p, q
        real(dp), intent(out):: val
        real(dp) :: rx, ry

        rx = (x - xc) / h_E
        ry = (y - yc) / h_E
        val = (rx**p) * (ry**q)
    end subroutine eval_monomial

    pure subroutine eval_monomial_grad(x, y, xc, yc, h_E, p, q, grad)
        real(dp), intent(in) :: x, y, xc, yc, h_E
        integer, intent(in)  :: p, q
        real(dp), intent(out):: grad(2)
        real(dp) :: rx, ry

        rx = (x - xc) / h_E
        ry = (y - yc) / h_E

        grad(1) = merge((real(p, dp) / h_E) * (rx**max(0, p - 1)) * (ry**q), 0.0_dp, p > 0)
        grad(2) = merge((real(q, dp) / h_E) * (rx**p) * (ry**max(0, q - 1)), 0.0_dp, q > 0)
    end subroutine eval_monomial_grad

    pure subroutine eval_monomial_div(x, y, xc, yc, h_E, p1, q1, p2, q2, div_val)
        real(dp), intent(in) :: x, y, xc, yc, h_E
        integer, intent(in)  :: p1, q1, p2, q2
        real(dp), intent(out):: div_val
        real(dp) :: g1(2), g2(2)

        call eval_monomial_grad(x, y, xc, yc, h_E, p1, q1, g1)
        call eval_monomial_grad(x, y, xc, yc, h_E, p2, q2, g2)
        div_val = g1(1) + g2(2)
    end subroutine eval_monomial_div

    pure subroutine eval_monomial_lap(x, y, xc, yc, h_E, p, q, lap)
        real(dp), intent(in) :: x, y, xc, yc, h_E
        integer, intent(in)  :: p, q
        real(dp), intent(out):: lap
        real(dp) :: rx, ry, d2x, d2y

        rx = (x - xc) / h_E
        ry = (y - yc) / h_E

        if (p >= 2) then
            d2x = (real(p * (p - 1), dp) / (h_E**2)) * (rx**(p - 2)) * (ry**q)
        else
            d2x = 0.0_dp
        end if

        if (q >= 2) then
            d2y = (real(q * (q - 1), dp) / (h_E**2)) * (rx**p) * (ry**(q - 2))
        else
            d2y = 0.0_dp
        end if

        lap = d2x + d2y
    end subroutine eval_monomial_lap

    pure subroutine eval_monomial_curl(x, y, xc, yc, h_E, p1, q1, p2, q2, curl_val)
        real(dp), intent(in) :: x, y, xc, yc, h_E
        integer, intent(in)  :: p1, q1, p2, q2
        real(dp), intent(out):: curl_val
        real(dp) :: g1(2), g2(2)

        call eval_monomial_grad(x, y, xc, yc, h_E, p1, q1, g1)
        call eval_monomial_grad(x, y, xc, yc, h_E, p2, q2, g2)
        curl_val = g2(1) - g1(2)
    end subroutine eval_monomial_curl

    pure subroutine eval_lap_coeffs(p, q, h_E, n_terms, p_out, q_out, coeffs)
        integer, intent(in)  :: p, q
        real(dp), intent(in) :: h_E
        integer, intent(out) :: n_terms
        integer, intent(out) :: p_out(2), q_out(2)
        real(dp), intent(out):: coeffs(2)

        n_terms = 0
        if (p >= 2) then
            n_terms = n_terms + 1
            p_out(n_terms)  = p - 2
            q_out(n_terms)  = q
            coeffs(n_terms) = real(p * (p - 1), dp) / (h_E**2)
        end if
        if (q >= 2) then
            n_terms = n_terms + 1
            p_out(n_terms)  = p
            q_out(n_terms)  = q - 2
            coeffs(n_terms) = real(q * (q - 1), dp) / (h_E**2)
        end if
    end subroutine eval_lap_coeffs

    pure subroutine eval_div_coeffs(p, q, h_E, dof_dir, n_terms, p_out, q_out, coeffs)
        integer, intent(in)  :: p, q, dof_dir
        real(dp), intent(in) :: h_E
        integer, intent(out) :: n_terms
        integer, intent(out) :: p_out(1), q_out(1)
        real(dp), intent(out):: coeffs(1)

        n_terms = 0
        select case (dof_dir)
        case (1)
            if (p >= 1) then
                n_terms = 1
                p_out(1)  = p - 1
                q_out(1)  = q
                coeffs(1) = real(p, dp) / h_E
            end if
        case (2)
            if (q >= 1) then
                n_terms = 1
                p_out(1)  = p
                q_out(1)  = q - 1
                coeffs(1) = real(q, dp) / h_E
            end if
        end select
    end subroutine eval_div_coeffs

    pure subroutine eval_grad_coeffs(p, q, h_E, dof_dir, n_terms, p_out, q_out, coeffs)
        integer, intent(in)  :: p, q, dof_dir
        real(dp), intent(in) :: h_E
        integer, intent(out) :: n_terms
        integer, intent(out) :: p_out(1), q_out(1)
        real(dp), intent(out):: coeffs(1)

        call eval_div_coeffs(p, q, h_E, dof_dir, n_terms, p_out, q_out, coeffs)
    end subroutine eval_grad_coeffs

end module math_geometry_mod