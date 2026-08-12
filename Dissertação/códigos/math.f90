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
    public :: integrate_by_reduction, one_point_quadrature, integrate_by_interpolation
    public :: get_num_monomials

    ! ---------------------------------------------------------
    ! Definição do Tipo Polinômio
    ! ---------------------------------------------------------
    type :: Polynomial
        integer :: degree
        real(8), allocatable :: coeffs(:)
    end type Polynomial

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
        integer :: i, iter, m, k
        real(dp) :: x, x_prev, p, dp_val, d2p, theta, pi
        real(dp) :: p_km1, p_k, p_kp1

        if (n <= 0) return

        select case (n)
        case (1)
            pts(1) = 0.5_dp
            w(1)   = 1.0_dp
        case (2)
            pts(1) = 0.0_dp; w(1) = 0.5_dp
            pts(2) = 1.0_dp; w(2) = 0.5_dp
        case default
            pi = 3.14159265358979323846_dp
            m = n - 1

            pts(1) = -1.0_dp
            pts(n) =  1.0_dp

            do i = 2, n - 1
                theta = pi * real(i - 1, dp) / real(n - 1, dp)
                x = -cos(theta)

                do iter = 1, 25
                    p_km1 = 1.0_dp
                    p_k   = x
                    do k = 1, m - 1
                        p_kp1 = (real(2*k + 1, dp) * x * p_k - real(k, dp) * p_km1) / real(k + 1, dp)
                        p_km1 = p_k
                        p_k   = p_kp1
                    end do
                    p = p_k
                    dp_val = real(m, dp) * (x * p_k - p_km1) / (x**2 - 1.0_dp)
                    d2p    = (2.0_dp * x * dp_val - real(m*(m + 1), dp) * p) / (1.0_dp - x**2)

                    if (abs(d2p) < 1.0e-30_dp) exit
                    x_prev = x
                    x = x - dp_val / d2p
                    if (abs(x - x_prev) < 1.0e-15_dp) exit
                end do
                pts(i) = x
            end do

            do i = 1, n
                x = pts(i)
                p_km1 = 1.0_dp
                p_k   = x
                do k = 1, m - 1
                    p_kp1 = (real(2*k + 1, dp) * x * p_k - real(k, dp) * p_km1) / real(k + 1, dp)
                    p_km1 = p_k
                    p_k   = p_kp1
                end do
                w(i) = 2.0_dp / (real(n*(n - 1), dp) * (p_k**2))
                pts(i) = (pts(i) + 1.0_dp) / 2.0_dp
                w(i)   = w(i) / 2.0_dp
            end do
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
        if (cross_sum /= 0.0_dp) then
            ! Correção: Usar 3.0 * cross_sum para manter a consistência do sinal com x_sum
            xc = x_sum / (3.0_dp * cross_sum)
            yc = y_sum / (3.0_dp * cross_sum)
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

    pure subroutine eval_monomial_grad(x, y, xc, yc, h_E, p, q, g)
        real(dp), intent(in) :: x, y, xc, yc, h_E
        integer, intent(in)  :: p, q
        real(dp), intent(out):: g(2)
        real(dp) :: rx, ry

        rx = (x - xc) / h_E
        ry = (y - yc) / h_E

        g(1) = merge((real(p, dp) / h_E) * (rx**max(0, p - 1)) * (ry**q), 0.0_dp, p > 0)
        g(2) = merge((real(q, dp) / h_E) * (rx**p) * (ry**max(0, q - 1)), 0.0_dp, q > 0)
    end subroutine eval_monomial_grad

    pure subroutine eval_monomial_div(x, y, xc, yc, h_E, p, q, dof_dir, div_val)
        real(dp), intent(in) :: x, y, xc, yc, h_E
        integer, intent(in)  :: p, q, dof_dir
        real(dp), intent(out):: div_val
        real(dp) :: g(2)

        call eval_monomial_grad(x, y, xc, yc, h_E, p, q, g)
        if (dof_dir == 1) then
            div_val = g(1)
        else if (dof_dir == 2) then
            div_val = g(2)
        else
            div_val = 0.0_dp
        end if
    end subroutine eval_monomial_div

    pure subroutine eval_monomial_lap(x, y, xc, yc, h_E, p, q, l_val)
        real(dp), intent(in) :: x, y, xc, yc, h_E
        integer, intent(in)  :: p, q
        real(dp), intent(out):: l_val
        real(dp) :: rx, ry, d2x, d2y

        rx = (x - xc) / h_E
        ry = (y - yc) / h_E

        d2x = merge((real(p * (p - 1), dp) / (h_E**2)) * (rx**max(0, p - 2)) * (ry**q), 0.0_dp, p >= 2)
        d2y = merge((real(q * (q - 1), dp) / (h_E**2)) * (rx**p) * (ry**max(0, q - 2)), 0.0_dp, q >= 2)
        l_val = d2x + d2y
    end subroutine eval_monomial_lap

    pure subroutine eval_monomial_curl(x, y, xc, yc, h_E, p, q, curl_val)
        real(dp), intent(in) :: x, y, xc, yc, h_E
        integer, intent(in)  :: p, q
        real(dp), intent(out):: curl_val(2)
        real(dp) :: g(2)

        call eval_monomial_grad(x, y, xc, yc, h_E, p, q, g)
        curl_val(1) =  g(2)
        curl_val(2) = -g(1)
    end subroutine eval_monomial_curl

    pure subroutine eval_lap_coeffs(p_in, q_in, h_E, n_terms, p_out, q_out, coeffs)
        integer, intent(in)  :: p_in, q_in
        real(dp), intent(in) :: h_E
        integer, intent(out) :: n_terms
        integer, intent(out) :: p_out(2), q_out(2)
        real(dp), intent(out):: coeffs(2)

        n_terms = 0
        p_out = 0; q_out = 0; coeffs = 0.0_dp

        if (p_in >= 2) then
            n_terms = n_terms + 1
            p_out(n_terms)  = p_in - 2
            q_out(n_terms)  = q_in
            coeffs(n_terms) = real(p_in * (p_in - 1), dp) / (h_E**2)
        end if

        if (q_in >= 2) then
            n_terms = n_terms + 1
            p_out(n_terms)  = p_in
            q_out(n_terms)  = q_in - 2
            coeffs(n_terms) = real(q_in * (q_in - 1), dp) / (h_E**2)
        end if
    end subroutine eval_lap_coeffs

    pure subroutine eval_div_coeffs(p_in, q_in, h_E, dof_dir, n_terms, p_out, q_out, coeffs)
        integer, intent(in)  :: p_in, q_in, dof_dir
        real(dp), intent(in) :: h_E
        integer, intent(out) :: n_terms
        integer, intent(out) :: p_out(1), q_out(1)
        real(dp), intent(out):: coeffs(1)

        n_terms = 0
        p_out = 0; q_out = 0; coeffs = 0.0_dp

        if (dof_dir == 1 .and. p_in >= 1) then
            n_terms = 1
            p_out(1)  = p_in - 1
            q_out(1)  = q_in
            coeffs(1) = real(p_in, dp) / h_E
        else if (dof_dir == 2 .and. q_in >= 1) then
            n_terms = 1
            p_out(1)  = p_in
            q_out(1)  = q_in - 1
            coeffs(1) = real(q_in, dp) / h_E
        end if
    end subroutine eval_div_coeffs

    pure subroutine eval_grad_coeffs(p_in, q_in, h_E, dof_dir, n_terms, p_out, q_out, coeffs)
        integer, intent(in)  :: p_in, q_in, dof_dir
        real(dp), intent(in) :: h_E
        integer, intent(out) :: n_terms
        integer, intent(out) :: p_out(1), q_out(1)
        real(dp), intent(out):: coeffs(1)

        call eval_div_coeffs(p_in, q_in, h_E, dof_dir, n_terms, p_out, q_out, coeffs)
    end subroutine eval_grad_coeffs

    !---------------------------------------------------------
    ! Integrações numéricas
    !---------------------------------------------------------

    subroutine integrate_by_reduction(poly_in, ortho_polys, domain_area, integral_result)
        type(Polynomial), intent(in)               :: poly_in
        type(Polynomial), dimension(:), intent(in) :: ortho_polys
        real(dp), intent(in)                       :: domain_area
        real(dp), intent(out)                      :: integral_result
        
        type(Polynomial) :: poly_remainder
        real(dp)         :: remainder_integral
        
        call get_polynomial_remainder(poly_in, ortho_polys, poly_remainder)
        call integrate_low_degree(poly_remainder, domain_area, remainder_integral)
        
        integral_result = remainder_integral
    end subroutine integrate_by_reduction

    subroutine get_polynomial_remainder(poly_in, ortho_polys, poly_remainder)
        type(Polynomial), intent(in)               :: poly_in
        type(Polynomial), dimension(:), intent(in) :: ortho_polys
        type(Polynomial), intent(out)              :: poly_remainder
        
        poly_remainder%degree = poly_in%degree
        if (allocated(poly_in%coeffs)) then
            poly_remainder%coeffs = poly_in%coeffs
        end if
    end subroutine get_polynomial_remainder

    subroutine integrate_low_degree(poly_remainder, domain_area, remainder_integral)
        type(Polynomial), intent(in) :: poly_remainder
        real(dp), intent(in)         :: domain_area
        real(dp), intent(out)        :: remainder_integral
        
        if (allocated(poly_remainder%coeffs) .and. size(poly_remainder%coeffs) > 0) then
            remainder_integral = poly_remainder%coeffs(1) * domain_area
        else
            remainder_integral = 0.0_dp
        end if
    end subroutine integrate_low_degree

    subroutine one_point_quadrature(func, x0, y0, domain_area, integral_result)
        ! Interface explicitando o uso do dp do módulo pai
        interface
            function func(x, y)
                import :: dp
                implicit none
                real(dp) :: func
                real(dp), intent(in) :: x, y
            end function func
        end interface
        
        real(dp), intent(in)  :: x0, y0
        real(dp), intent(in)  :: domain_area
        real(dp), intent(out) :: integral_result
        
        integral_result = func(x0, y0) * domain_area
        
    end subroutine one_point_quadrature

    subroutine integrate_by_interpolation(poly_matrix, f_vals, n_pts, domain_area, integral_result)
        integer, intent(in)   :: n_pts
        ! poly_matrix(i, j) contém o polinômio 'j' avaliado no ponto 'i'
        real(dp), intent(in)  :: poly_matrix(n_pts, n_pts) 
        ! f_vals(i) contém o valor da função no ponto 'i'
        real(dp), intent(in)  :: f_vals(n_pts)             
        real(dp), intent(in)  :: domain_area
        real(dp), intent(out) :: integral_result
        
        real(dp) :: A_copy(n_pts, n_pts)
        real(dp) :: coeffs(n_pts)
        
        ! O solver do LAPACK (dgesv) sobrescreve a matriz original com os 
        ! fatores LU. Por isso, precisamos trabalhar com uma cópia local.
        A_copy = poly_matrix
        
        ! Resolve o sistema: A_copy * coeffs = f_vals
        ! A rotina 'solve_linear_system' foi reutilizada do seu módulo
        call solve_linear_system(A_copy, f_vals, coeffs, n_pts)
        
        ! A integral é o coeficiente do termo constante multiplicado pela área.
        integral_result = coeffs(1) * domain_area
        
    end subroutine integrate_by_interpolation

    !---------------------------------------------------------
    !Matemática semi-simbólica para integração de polinômios
    !---------------------------------------------------------

    pure function get_num_monomials(k) result(n)
        integer, intent(in) :: k
        integer :: n
        n = ((k + 1) * (k + 2)) / 2
    end function get_num_monomials

    subroutine evaluate_monomial_val(x, y, xc, yc, h_E, p, q, dof_dir, op_val)
        real(dp), intent(in) :: x, y, xc, yc, h_E
        integer, intent(in)  :: p, q, dof_dir
        real(dp), intent(out):: op_val(:)
        
        real(dp) :: rx, ry
        rx = (x - xc) / h_E
        ry = (y - yc) / h_E
        
        op_val = 0.0_dp
        if (dof_dir >= 1 .and. dof_dir <= size(op_val)) then
            op_val(dof_dir) = (rx**p) * (ry**q)
        end if
    end subroutine evaluate_monomial_val

    subroutine evaluate_monomial_grad(x, y, xc, yc, h_E, p, q, dof_dir, op_val)
        real(dp), intent(in) :: x, y, xc, yc, h_E
        integer, intent(in)  :: p, q, dof_dir
        real(dp), intent(out):: op_val(:)
        
        real(dp) :: rx, ry, dx, dy
        rx = (x - xc) / h_E
        ry = (y - yc) / h_E
        
        dx = merge((real(p, dp) / h_E) * (rx**max(0, p - 1)) * (ry**q), 0.0_dp, p > 0)
        dy = merge((real(q, dp) / h_E) * (rx**p) * (ry**max(0, q - 1)), 0.0_dp, q > 0)
        
        op_val = 0.0_dp
        if (size(op_val) >= 2) then
            op_val(1) = dx
            op_val(2) = dy
        end if
    end subroutine evaluate_monomial_grad

    subroutine evaluate_monomial_div(x, y, xc, yc, h_E, p, q, dof_dir, op_val)
        real(dp), intent(in) :: x, y, xc, yc, h_E
        integer, intent(in)  :: p, q, dof_dir
        real(dp), intent(out):: op_val(:)
        
        real(dp) :: rx, ry
        rx = (x - xc) / h_E
        ry = (y - yc) / h_E
        
        op_val = 0.0_dp
        if (size(op_val) >= 1) then
            if (dof_dir == 1) then
                ! Componente X: derivada em relação a x
                op_val(1) = merge((real(p, dp) / h_E) * (rx**max(0, p - 1)) * (ry**q), 0.0_dp, p > 0)
            else if (dof_dir == 2) then
                ! Componente Y: derivada em relação a y
                op_val(1) = merge((real(q, dp) / h_E) * (rx**p) * (ry**max(0, q - 1)), 0.0_dp, q > 0)
            end if
        end if
    end subroutine evaluate_monomial_div

    subroutine evaluate_monomial_curl(x, y, xc, yc, h_E, p, q, dof_dir, op_val)
        real(dp), intent(in) :: x, y, xc, yc, h_E
        integer, intent(in)  :: p, q, dof_dir
        real(dp), intent(out):: op_val(:)
        
        real(dp) :: rx, ry
        rx = (x - xc) / h_E
        ry = (y - yc) / h_E
        
        op_val = 0.0_dp
        if (size(op_val) >= 1) then
            if (dof_dir == 1) then
                ! Base (m, 0): Rotacional = - dm/dy
                op_val(1) = -merge((real(q, dp) / h_E) * (rx**p) * (ry**max(0, q - 1)), 0.0_dp, q > 0)
            else if (dof_dir == 2) then
                ! Base (0, m): Rotacional = + dm/dx
                op_val(1) = merge((real(p, dp) / h_E) * (rx**max(0, p - 1)) * (ry**q), 0.0_dp, p > 0)
            end if
        end if
    end subroutine evaluate_monomial_curl

    subroutine evaluate_monomial_lap(x, y, xc, yc, h_E, p, q, dof_dir, op_val)
        real(dp), intent(in) :: x, y, xc, yc, h_E
        integer, intent(in)  :: p, q, dof_dir
        real(dp), intent(out):: op_val(:)
        
        real(dp) :: rx, ry, dxx, dyy
        rx = (x - xc) / h_E
        ry = (y - yc) / h_E
        
        dxx = 0.0_dp
        if (p >= 2) dxx = (real(p, dp) * real(p - 1, dp) / (h_E**2)) * (rx**(p - 2)) * (ry**q)
        
        dyy = 0.0_dp
        if (q >= 2) dyy = (real(q, dp) * real(q - 1, dp) / (h_E**2)) * (rx**p) * (ry**(q - 2))
        
        op_val = 0.0_dp
        if (size(op_val) >= 1) op_val(1) = dxx + dyy
    end subroutine evaluate_monomial_lap

    subroutine evaluate_monomial_hessian(x, y, xc, yc, h_E, p, q, dof_dir, op_val)
        real(dp), intent(in) :: x, y, xc, yc, h_E
        integer, intent(in)  :: p, q, dof_dir
        real(dp), intent(out):: op_val(:)
        
        real(dp) :: rx, ry, dxx, dyy, dxy
        rx = (x - xc) / h_E
        ry = (y - yc) / h_E
        
        dxx = 0.0_dp
        if (p >= 2) dxx = (real(p, dp) * real(p - 1, dp) / (h_E**2)) * (rx**(p - 2)) * (ry**q)
        
        dyy = 0.0_dp
        if (q >= 2) dyy = (real(q, dp) * real(q - 1, dp) / (h_E**2)) * (rx**p) * (ry**(q - 2))
        
        dxy = 0.0_dp
        if (p >= 1 .and. q >= 1) dxy = (real(p * q, dp) / (h_E**2)) * (rx**(p - 1)) * (ry**(q - 1))
        
        op_val = 0.0_dp
        if (size(op_val) >= 3) then
            op_val(1) = dxx ! Derivada parcial dupla em x
            op_val(2) = dxy ! Derivada mista em xy
            op_val(3) = dyy ! Derivada parcial dupla em y
        end if
    end subroutine evaluate_monomial_hessian

end module math_geometry_mod