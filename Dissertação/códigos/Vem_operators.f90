!===============================================================================
! MÓDULO DE OPERADORES GENÉRICOS VEM
!===============================================================================
module vem_operators_mod
    use math_geometry_mod
    implicit none
    private

    public :: compute_matrix_B_boundary, compute_matrix_B_domain, compute_matrix_B, invert_and_project
    public :: compute_D, compute_D_eps
    public :: compute_matrix_G, apply_P0
    public :: operator_interface, coeff_operator_interface

    abstract interface
        subroutine operator_interface(x_pt, y_pt, xc, yc, h_E, p_exp, q_exp, &
                                     normal, dof_dir, op_val)
            import :: dp
            real(dp), intent(in) :: x_pt, y_pt, xc, yc, h_E
            integer, intent(in)  :: p_exp, q_exp, dof_dir
            real(dp), intent(in) :: normal(2)
            real(dp), intent(out):: op_val(:)
        end subroutine operator_interface
    end interface

    abstract interface
        subroutine coeff_operator_interface(p_in, q_in, h_E, dof_dir, n_terms, p_out, q_out, coeffs)
            import :: dp
            integer, intent(in)  :: p_in, q_in, dof_dir
            real(dp), intent(in) :: h_E
            integer, intent(out) :: n_terms
            integer, intent(out) :: p_out(:), q_out(:)
            real(dp), intent(out):: coeffs(:)
        end subroutine coeff_operator_interface
    end interface

contains

    subroutine compute_matrix_B(B, n_mon, n_dofs, n_verts, deg_k, dim_op, &
                                edge_dof_map, edge_length, normals, &
                                gauss_pts, gauss_w, n_gauss, &
                                node_coords, xc, yc, h_E, op_func, &
                                area, n_internal, internal_dof_map, &
                                p_internal, q_internal, op_coeff_func)
        integer, intent(in)  :: n_mon, n_dofs, n_verts, deg_k, dim_op, n_gauss
        integer, intent(in)  :: edge_dof_map(n_gauss, 2, n_verts)
        real(dp), intent(in) :: edge_length(n_verts), normals(2, n_verts)
        real(dp), intent(in) :: gauss_pts(n_gauss), gauss_w(n_gauss)
        real(dp), intent(in) :: node_coords(:,:), xc, yc, h_E
        real(dp), intent(in), optional :: area
        integer, intent(in), optional  :: n_internal
        integer, intent(in), optional  :: internal_dof_map(:,:)
        integer, intent(in), optional  :: p_internal(:), q_internal(:)
        procedure(coeff_operator_interface), optional :: op_coeff_func
        procedure(operator_interface) :: op_func
        real(dp), intent(out):: B(dim_op * n_mon, n_dofs)

        B = 0.0_dp
        call compute_matrix_B_boundary(B, n_mon, n_dofs, n_verts, deg_k, dim_op, &
                                       edge_dof_map, edge_length, normals, &
                                       gauss_pts, gauss_w, n_gauss, &
                                       node_coords, xc, yc, h_E, op_func)
        if (dim_op >= 2 .and. present(op_coeff_func) .and. present(area)) then
            call compute_matrix_B_domain(B, n_mon, n_dofs, dim_op, deg_k, &
                                         area, h_E, n_internal, internal_dof_map, &
                                         p_internal, q_internal, op_coeff_func)
        end if
    end subroutine compute_matrix_B

    subroutine compute_matrix_B_domain(B, n_mon, n_dofs, dim_op, deg_k, area, h_E, &
                                       n_internal, internal_dof_map, &
                                       p_internal, q_internal, op_coeff_func, sign_factor)
        integer, intent(in)     :: n_mon, n_dofs, dim_op, deg_k
        integer, intent(in), optional :: n_internal
        real(dp), intent(in)    :: area, h_E
        integer, intent(in), optional :: internal_dof_map(:,:)
        integer, intent(in), optional :: p_internal(:), q_internal(:)
        procedure(coeff_operator_interface) :: op_coeff_func
        real(dp), intent(in), optional :: sign_factor
        real(dp), intent(inout) :: B(dim_op * n_mon, n_dofs)

        integer :: m, d, t, n_terms, match_idx, col_idx, row_idx, n_int
        real(dp) :: s_factor
        integer :: p_out(2), q_out(2)
        real(dp) :: coeffs(2)
        integer, allocatable :: p_mon(:), q_mon(:)

        if (.not. present(n_internal) .or. .not. present(internal_dof_map)) return

        n_int = n_internal
        s_factor = 1.0_dp
        if (present(sign_factor)) s_factor = sign_factor

        allocate(p_mon(n_mon), q_mon(n_mon))
        call get_monomial_exponents(deg_k, n_mon, p_mon, q_mon)

        do m = 1, n_mon
            do d = 1, 2
                call op_coeff_func(p_mon(m), q_mon(m), h_E, d, n_terms, p_out, q_out, coeffs)

                do t = 1, n_terms
                    match_idx = find_monomial_index(p_out(t), q_out(t), p_internal, q_internal, n_int)

                    if (match_idx > 0) then
                        col_idx = internal_dof_map(match_idx, d)
                        if (col_idx > 0 .and. col_idx <= n_dofs) then
                            row_idx = dim_op * (m - 1) + 1
                            B(row_idx, col_idx) = B(row_idx, col_idx) + s_factor * coeffs(t) * area
                        end if
                    end if
                end do
            end do
        end do
        deallocate(p_mon, q_mon)
    end subroutine compute_matrix_B_domain

    subroutine compute_matrix_B_boundary(B, n_mon, n_dofs, n_verts, deg_k, dim_op, &
                                   edge_dof_map, edge_length, normals, &
                                   gauss_pts, gauss_w, n_gauss, &
                                   node_coords, xc, yc, h_E, op_func)
        integer, intent(in)  :: n_mon, n_dofs, n_verts, deg_k, dim_op, n_gauss
        integer, intent(in)  :: edge_dof_map(n_gauss, 2, n_verts)
        real(dp), intent(in) :: edge_length(n_verts), normals(2, n_verts)
        real(dp), intent(in) :: gauss_pts(n_gauss), gauss_w(n_gauss)
        real(dp), intent(in) :: node_coords(:,:), xc, yc, h_E
        procedure(operator_interface) :: op_func
        real(dp), intent(inout) :: B(dim_op * n_mon, n_dofs)

        integer :: e, g, d, m, r, col_idx, row_idx, v1, v2
        real(dp) :: t, w, wt, x_pt, y_pt, x1, y1, x2, y2
        real(dp) :: op_val(dim_op)
        integer, allocatable :: p_mon(:), q_mon(:)

        allocate(p_mon(n_mon), q_mon(n_mon))
        call get_monomial_exponents(deg_k, n_mon, p_mon, q_mon)

        do e = 1, n_verts
            v1 = e
            v2 = mod(e, n_verts) + 1
            x1 = node_coords(v1, 1); y1 = node_coords(v1, 2)
            x2 = node_coords(v2, 1); y2 = node_coords(v2, 2)

            do g = 1, n_gauss
                t = gauss_pts(g)
                w = gauss_w(g)
                x_pt = x1 + t * (x2 - x1)
                y_pt = y1 + t * (y2 - y1)
                wt = w * edge_length(e)

                do d = 1, 2
                    col_idx = edge_dof_map(g, d, e)
                    if (col_idx > 0 .and. col_idx <= n_dofs) then
                        do m = 1, n_mon
                            call op_func(x_pt, y_pt, xc, yc, h_E, p_mon(m), q_mon(m), &
                                         normals(:, e), d, op_val)
                            do r = 1, dim_op
                                row_idx = dim_op * (m - 1) + r
                                B(row_idx, col_idx) = B(row_idx, col_idx) + wt * op_val(r)
                            end do
                        end do
                    end if
                end do
            end do
        end do

        deallocate(p_mon, q_mon)
    end subroutine compute_matrix_B_boundary

    !===========================================================================
    ! COMPUTE_D_GENERIC
    ! Calcula a matriz D = dof_i(p_alpha) para qualquer ordem k >= 1
    !===========================================================================
    subroutine compute_D(D_mat, boundary_pts, n_boundary_pts, k, xc, yc, h_E, &
                                ndof, area, polygon_coords, n_verts)
        real(dp), intent(out) :: D_mat(:,:)        ! Dimensão: (n_dofs x (ndof * n_mon))
        real(dp), intent(in)  :: boundary_pts(:,:) ! Coordenadas dos nós de borda (n_boundary_pts x 2)
        integer, intent(in)   :: n_boundary_pts    ! Número total de nós na borda
        integer, intent(in)   :: k                 ! Ordem polinomial do VEM
        real(dp), intent(in)  :: xc, yc, h_E       ! Geometria do elemento
        integer, intent(in)   :: ndof              ! Graus de liberdade por nó (1 para escalar, 2 para 2D)
        real(dp), intent(in), optional :: area     ! Área do elemento |E| (necessário se k >= 2)
        real(dp), intent(in), optional :: polygon_coords(:,:) ! Vértices do polígono (necessário se k >= 2)
        integer, intent(in), optional  :: n_verts  ! Número de vértices (necessário se k >= 2)

        integer :: n_mon, n_mon_int, v, m, m_int, d, col_base, row_idx, k_int
        integer, allocatable :: p_exp(:), q_exp(:)         ! Expoentes da base P_k
        integer, allocatable :: p_int(:), q_int(:)         ! Expoentes da base P_{k-2}
        real(dp) :: rx, ry, m_val, int_mon_val

        ! 1. Número de monômios da base P_k
        n_mon = ((k + 1) * (k + 2)) / 2

        allocate(p_exp(n_mon), q_exp(n_mon))
        call get_monomial_exponents(k, n_mon, p_exp, q_exp)

        D_mat = 0.0_dp

        !-----------------------------------------------------------------------
        ! PARTE 1: Graus de Liberdade de Borda (Avaliação Pontual dos Nós)
        !-----------------------------------------------------------------------
        do v = 1, n_boundary_pts
            rx = (boundary_pts(v,1) - xc) / h_E
            ry = (boundary_pts(v,2) - yc) / h_E

            do m = 1, n_mon
                m_val = (rx**p_exp(m)) * (ry**q_exp(m))
                
                ! Mapeamento bloco a bloco segundo ndof
                do d = 1, ndof
                    row_idx  = ndof * (v - 1) + d
                    col_base = ndof * (m - 1) + d
                    D_mat(row_idx, col_base) = m_val
                end do
            end do
        end do

        !-----------------------------------------------------------------------
        ! PARTE 2: Graus de Liberdade Internos / Momento de Domínio (para k >= 2)
        !-----------------------------------------------------------------------
        k_int = k - 2
        if (k_int >= 0) then
            if (.not. present(area) .or. .not. present(polygon_coords) .or. .not. present(n_verts)) then
                error stop "compute_D_generic: Area e dados do poligono sao obrigatorios para k >= 2"
            end if

            n_mon_int = ((k_int + 1) * (k_int + 2)) / 2
            allocate(p_int(n_mon_int), q_int(n_mon_int))
            call get_monomial_exponents(k_int, n_mon_int, p_int, q_int)

            do m_int = 1, n_mon_int
                do m = 1, n_mon
                    ! Integral do produto dos monômios: \int_E (p_alpha * m_m) dE / |E|
                    int_mon_val = compute_monomial_domain_integral( &
                                    p_exp(m) + p_int(m_int), q_exp(m) + q_int(m_int), &
                                    xc, yc, h_E, polygon_coords, n_verts) / area

                    do d = 1, ndof
                        ! Linha do DoF interno correspondente no vetor global do elemento
                        row_idx  = ndof * n_boundary_pts + ndof * (m_int - 1) + d
                        col_base = ndof * (m - 1) + d
                        D_mat(row_idx, col_base) = int_mon_val
                    end do
                end do
            end do

            deallocate(p_int, q_int)
        end if

        deallocate(p_exp, q_exp)
    end subroutine compute_D

    !===========================================================================
    ! COMPUTE_D_EPS_GENERIC
    ! Calcula a matriz D_eps para a projeção de deformações para qualquer k >= 1
    !===========================================================================
    subroutine compute_D_eps(D_mat, boundary_pts, n_boundary_pts, k_eps, &
                                    xc, yc, h_E, ndof, area, polygon_coords, n_verts)
        real(dp), intent(out) :: D_mat(:,:)        ! Dimensão: (n_dofs x (3 * n_mon_eps))
        real(dp), intent(in)  :: boundary_pts(:,:) ! Coordenadas dos nós de borda
        integer, intent(in)   :: n_boundary_pts    ! Número de nós na borda
        integer, intent(in)   :: k_eps             ! Grau da deformacao (k_eps = k - 1)
        real(dp), intent(in)  :: xc, yc, h_E       ! Centroide e diâmetro
        integer, intent(in)   :: ndof              ! Graus de liberdade por nó (ndof = 2)
        real(dp), intent(in), optional :: area     ! Área do elemento
        real(dp), intent(in), optional :: polygon_coords(:,:)
        integer, intent(in), optional  :: n_verts

        integer :: n_mon_eps, n_mon_int, v, m, m_int, dof_x, dof_y, col_base, k_int, row_x, row_y
        integer, allocatable :: p_exp(:), q_exp(:)
        integer, allocatable :: p_int(:), q_int(:)
        real(dp) :: rx, ry, m_val, int_mon_val

        n_mon_eps = ((k_eps + 1) * (k_eps + 2)) / 2

        allocate(p_exp(n_mon_eps), q_exp(n_mon_eps))
        call get_monomial_exponents(k_eps, n_mon_eps, p_exp, q_exp)

        D_mat = 0.0_dp

        !-----------------------------------------------------------------------
        ! 1. Nós de Borda: Mapeamento dos DoFs de Deslocamento Nodal
        !-----------------------------------------------------------------------
        do v = 1, n_boundary_pts
            rx = (boundary_pts(v,1) - xc) / h_E
            ry = (boundary_pts(v,2) - yc) / h_E

            dof_x = ndof * (v - 1) + 1
            dof_y = ndof * (v - 1) + 2

            do m = 1, n_mon_eps
                m_val = (rx**p_exp(m)) * (ry**q_exp(m))
                col_base = 3 * (m - 1)

                ! Componente u_x afeta eps_xx (col 1) e u_y afeta eps_yy (col 2)
                D_mat(dof_x, col_base + 1) = m_val
                D_mat(dof_y, col_base + 2) = m_val
            end do
        end do

        !-----------------------------------------------------------------------
        ! 2. DoFs de Domínio (Se k_eps >= 1, isto é, k >= 2)
        !-----------------------------------------------------------------------
        k_int = k_eps - 1  ! Grau dos momentos internos para deformação = k - 2
        if (k_int >= 0) then
            if (.not. present(area) .or. .not. present(polygon_coords) .or. .not. present(n_verts)) then
                error stop "compute_D_eps_generic: Dados geometricos ausentes para k >= 2"
            end if

            n_mon_int = ((k_int + 1) * (k_int + 2)) / 2
            allocate(p_int(n_mon_int), q_int(n_mon_int))
            call get_monomial_exponents(k_int, n_mon_int, p_int, q_int)

            do m_int = 1, n_mon_int
                do m = 1, n_mon_eps
                    int_mon_val = compute_monomial_domain_integral( &
                                    p_exp(m) + p_int(m_int), q_exp(m) + q_int(m_int), &
                                    xc, yc, h_E, polygon_coords, n_verts) / area

                    row_x = ndof * n_boundary_pts + ndof * (m_int - 1) + 1
                    row_y = ndof * n_boundary_pts + ndof * (m_int - 1) + 2
                    col_base = 3 * (m - 1)

                    D_mat(row_x, col_base + 1) = int_mon_val
                    D_mat(row_y, col_base + 2) = int_mon_val
                end do
            end do

            deallocate(p_int, q_int)
        end if

        deallocate(p_exp, q_exp)
    end subroutine compute_D_eps

    subroutine compute_matrix_G(G, k_eps, xc, yc, h_E, area, polygon_coords, n_verts)
        ! G matrix size: (3 * n_mon) x (3 * n_mon)
        ! k_eps = k - 1 (degree of strain polynomials)
        integer, intent(in)  :: k_eps, n_verts
        real(dp), intent(in) :: xc, yc, h_E, area, polygon_coords(n_verts, 2)
        real(dp), intent(out):: G(:,:)

        integer :: n_mon, m, j, r_base, c_base
        real(dp) :: int_qq ! integral of q_m * q_j over element E
        integer, allocatable :: p_exp(:), q_exp(:)

        n_mon = ((k_eps + 1) * (k_eps + 2)) / 2
        allocate(p_exp(n_mon), q_exp(n_mon))
        call get_monomial_exponents(k_eps, n_mon, p_exp, q_exp)

        G = 0.0_dp

        do m = 1, n_mon
            do j = 1, n_mon
                ! Compute \int_E q_m * q_j dE using polygonal domain integration
                int_qq = compute_monomial_domain_integral(p_exp(m)+p_exp(j), q_exp(m)+q_exp(j), &
                                                       xc, yc, h_E, polygon_coords, n_verts)
            
                ! Assemble 3x3 block: int_qq * I_3
                r_base = 3 * (m - 1)
                c_base = 3 * (j - 1)

                G(r_base + 1, c_base + 1) = int_qq ! eps_xx
                G(r_base + 2, c_base + 2) = int_qq ! eps_yy
                G(r_base + 3, c_base + 3) = int_qq ! gamma_xy
         end do
        end do

        deallocate(p_exp, q_exp)
    end subroutine compute_matrix_G

    subroutine apply_P0(B, k_order, n_verts, n_monomials, loc_dof, ndof)
        implicit none
        integer, intent(in) :: k_order, n_verts, n_monomials, loc_dof, ndof
        real(dp), intent(inout) :: B(ndof*n_monomials, loc_dof)

        integer :: d, v, row_idx, col_idx, boundary_dofs, int_dof_idx

        boundary_dofs = n_verts * k_order * ndof

        do d = 1, ndof
            ! Row index for constant monomial m_1 = 1 for component d
            row_idx = d
            B(row_idx, :) = 0.0_dp

            if (k_order == 1) then
                ! For k = 1: P_0(phi_i) = 1 / N_V for vertex DOFs matching component d
                do v = 1, n_verts
                    col_idx = ndof * (v - 1) + d
                    B(row_idx, col_idx) = 1.0_dp / real(n_verts, dp)
                end do
            else
                ! For k >= 2: P_0(phi_i) = 1 for the 1st internal moment DOF (m_1 = 1)
                int_dof_idx = boundary_dofs + d
                if (int_dof_idx <= loc_dof) then
                    B(row_idx, int_dof_idx) = 1.0_dp
                end if
            end if
        end do
    end subroutine apply_P0

    subroutine invert_and_project(G, B, Pi)
        real(dp), intent(in)  :: G(:,:), B(:,:)
        real(dp), intent(out) :: Pi(:,:)

        integer :: n_rows, n_cols, info
        integer, allocatable :: ipiv(:)
        real(dp), allocatable :: G_lu(:,:)

        n_rows = size(G, 1)
        n_cols = size(B, 2)

        allocate(G_lu(n_rows, n_rows), ipiv(n_rows))
        G_lu = G
        Pi = B

        call dgetrf(n_rows, n_rows, G_lu, n_rows, ipiv, info)
        if (info /= 0) error stop "Erro na fatoracao LU de G em invert_and_project"

        call dgetrs('N', n_rows, n_cols, G_lu, n_rows, ipiv, Pi, n_rows, info)
        if (info /= 0) error stop "Erro na resolucao do sistema G * Pi = B"

        deallocate(G_lu, ipiv)
    end subroutine invert_and_project

end module vem_operators_mod