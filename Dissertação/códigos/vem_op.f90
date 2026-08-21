!===============================================================================
! MÓDULO DE OPERADORES GENÉRICOS VEM (ATUALIZADO E CONSOLIDADO PARA k >= 1)
!===============================================================================
module vem_operators_mod
    use math_geometry_mod
    implicit none
    private

    public :: compute_matrix_B_boundary, compute_matrix_B_domain, compute_matrix_B, invert_and_project
    public :: compute_matrix_D, compute_matrix_Q, compute_matrix_H0_exact
    public :: compute_matrix_G, apply_P0
    public :: operator_interface, coeff_operator_interface

    interface
        subroutine operator_interface(x_pt, y_pt, xc, yc, h_E, p_exp, q_exp, &
                                     normal, dof_dir, op_val)
            import :: dp
            real(dp), intent(in) :: x_pt, y_pt, xc, yc, h_E
            integer, intent(in)  :: p_exp, q_exp, dof_dir
            real(dp), intent(in) :: normal(2)
            real(dp), intent(out):: op_val(:)
        end subroutine operator_interface
    end interface

    interface
        subroutine coeff_operator_interface(p_in, q_in, h_E, dof_dir, n_terms, p_out, q_out, coeffs)
            import :: dp
            integer, intent(in)  :: p_in, q_in, dof_dir
            real(dp), intent(in) :: h_E
            integer, intent(out) :: n_terms
            integer, intent(out) :: p_out(:), q_out(:)
            real(dp), intent(out):: coeffs(:)
        end subroutine coeff_operator_interface
    end interface

    abstract interface
        subroutine op_coeff_interface(p_in, q_in, h_E, row_op, dof_dir, n_terms, p_out, q_out, coeffs)
            import :: dp               
            implicit none
            integer, intent(in)  :: p_in, q_in
            real(dp), intent(in) :: h_E
            integer, intent(in)  :: row_op
            integer, intent(in)  :: dof_dir
            integer, intent(out) :: n_terms
            integer, intent(out) :: p_out(:), q_out(:)
            real(dp), intent(out):: coeffs(:)
        end subroutine op_coeff_interface
    end interface
    
    contains

    subroutine compute_matrix_B(B, n_mon, n_dofs, n_verts, deg_k, dim_op, &
                            edge_dof_map, edge_length, normals, &
                            gauss_pts, gauss_w, n_gauss, &
                            node_coords, xc, yc, h_E, op_func, &
                            area, n_internal, internal_dof_map, &
                            p_internal, q_internal, op_coeff_func, sign_factor)
        implicit none
        integer, intent(in)  :: n_mon, n_dofs, n_verts, deg_k, dim_op, n_gauss
        integer, intent(in)  :: edge_dof_map(n_gauss, 2, n_verts)
        real(dp), intent(in) :: edge_length(n_verts), normals(2, n_verts)
        real(dp), intent(in) :: gauss_pts(n_gauss), gauss_w(n_gauss)
        real(dp), intent(in) :: node_coords(:,:), xc, yc, h_E
        real(dp), intent(in), optional :: area
        integer, intent(in), optional  :: n_internal
        integer, intent(in), optional  :: internal_dof_map(:,:)
        integer, intent(in), optional  :: p_internal(:), q_internal(:)
        procedure(op_coeff_interface), optional :: op_coeff_func
        procedure(operator_interface) :: op_func
        real(dp), intent(in), optional :: sign_factor
        real(dp), intent(out):: B(dim_op * n_mon, n_dofs)

        B = 0.0_dp

        ! 1. Contribuição de Contorno
        call compute_matrix_B_boundary(B, n_mon, n_dofs, n_verts, deg_k, dim_op, &
                                   edge_dof_map, edge_length, normals, &
                                   gauss_pts, gauss_w, n_gauss, &
                                   node_coords, xc, yc, h_E, op_func)

        ! 2. Contribuição de Domínio (k >= 2)
        if (present(op_coeff_func) .and. present(area) .and. present(n_internal) .and. &
            present(internal_dof_map) .and. present(p_internal) .and. present(q_internal)) then
        
            call compute_matrix_B_domain(B, n_mon, n_dofs, dim_op, deg_k, area, h_E, &
                                   op_coeff_func, n_internal, internal_dof_map, &
                                   p_internal, q_internal, sign_factor=sign_factor)
        end if
    end subroutine compute_matrix_B

    subroutine compute_matrix_B_domain(B, n_mon, n_dofs, dim_op, deg_k, area, h_E, &
                                   op_coeff_func, n_internal, internal_dof_map, &
                                   p_internal, q_internal, sign_factor)
    
        implicit none
    
        integer, intent(in)           :: n_mon, n_dofs, dim_op, deg_k
        real(dp), intent(in)          :: area, h_E
        procedure(op_coeff_interface) :: op_coeff_func
        integer, intent(in), optional :: n_internal
        integer, intent(in), optional :: internal_dof_map(:,:)
        integer, intent(in), optional :: p_internal(:), q_internal(:)
        real(dp), intent(in), optional :: sign_factor
        real(dp), intent(inout)       :: B(dim_op * n_mon, n_dofs)

        integer  :: m, d, dof_dir, t, n_terms, match_idx, col_idx, row_idx
        integer  :: n_int, max_dof_dim
        real(dp) :: s_factor
        integer  :: p_out(2), q_out(2)
        real(dp) :: coeffs(2)
        integer, allocatable :: p_mon(:), q_mon(:)

        if (.not. present(n_internal) .or. .not. present(internal_dof_map) .or. &
            .not. present(p_internal) .or. .not. present(q_internal)) return

        n_int = n_internal
        if (n_int <= 0) return

        s_factor = 1.0_dp
        if (present(sign_factor)) s_factor = sign_factor

        max_dof_dim = size(internal_dof_map, 2)

        allocate(p_mon(n_mon), q_mon(n_mon))
        call get_monomial_exponents(deg_k, n_mon, p_mon, q_mon)

        do m = 1, n_mon
            do d = 1, dim_op
                row_idx = dim_op * (m - 1) + d
            
                do dof_dir = 1, max_dof_dim
                    call op_coeff_func(p_mon(m), q_mon(m), h_E, d, dof_dir, n_terms, p_out, q_out, coeffs)
                
                    do t = 1, n_terms
                        match_idx = find_monomial_index(p_out(t), q_out(t), p_internal, q_internal, n_int)
                        if (match_idx > 0) then
                            col_idx = internal_dof_map(match_idx, dof_dir)
                            if (col_idx > 0 .and. col_idx <= n_dofs) then
                                B(row_idx, col_idx) = B(row_idx, col_idx) + s_factor * coeffs(t) * area
                            end if
                        end if
                    end do
                end do
            
            end do
        end do

        deallocate(p_mon, q_mon)
    end subroutine compute_matrix_B_domain

    subroutine compute_matrix_B_boundary(B, n_mon, n_dofs, n_verts, deg_k, dim_op, &
                               edge_dof_map, edge_length, normals, &
                               gauss_pts, gauss_w, n_gauss, &
                               node_coords, xc, yc, h_E, op_func)
                               
        implicit none
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

    subroutine compute_matrix_D(D_mat, boundary_pts, n_boundary_pts, k, xc, yc, h_E, &
                     ndof, loc_ndof, area, polygon_coords, n_verts, &
                     boundary_dof_map, internal_dof_map)
        implicit none
        
        real(dp), intent(out) :: D_mat(:,:)        
        real(dp), intent(in)  :: boundary_pts(:,:) 
        integer, intent(in)   :: n_boundary_pts    
        integer, intent(in)   :: k                 
        real(dp), intent(in)  :: xc, yc, h_E       
        integer, intent(in)   :: ndof              
        integer, intent(in)   :: loc_ndof          
        real(dp), intent(in), optional :: area     
        real(dp), intent(in), optional :: polygon_coords(:,:) 
        integer, intent(in), optional  :: n_verts  
        integer, intent(in), optional  :: boundary_dof_map(:,:) 
        integer, intent(in), optional  :: internal_dof_map(:,:) 

        integer :: n_mon, n_mon_int, v, m, m_int, d, col_base, row_idx, k_int
        integer, allocatable :: p_exp(:), q_exp(:)         
        integer, allocatable :: p_int(:), q_int(:)         
        real(dp) :: rx, ry, m_val, int_mon_val

        n_mon = ((k + 1) * (k + 2)) / 2
        allocate(p_exp(n_mon), q_exp(n_mon))
        call get_monomial_exponents(k, n_mon, p_exp, q_exp)

        D_mat = 0.0_dp

        ! 1. DOFs de contorno: Avaliação pontual escalada
        do v = 1, n_boundary_pts
            rx = (boundary_pts(v,1) - xc) / h_E
            ry = (boundary_pts(v,2) - yc) / h_E

            do m = 1, n_mon
                m_val = (rx**p_exp(m)) * (ry**q_exp(m))
                do d = 1, ndof
                    if (present(boundary_dof_map)) then
                        row_idx = boundary_dof_map(v, d)
                    else
                        row_idx = ndof * (v - 1) + d
                    end if
                    
                    col_base = ndof * (m - 1) + d
                    
                    if (row_idx > 0 .and. row_idx <= loc_ndof) then
                        D_mat(row_idx, col_base) = m_val
                    end if
                end do
            end do
        end do

        ! 2. DOFs internos (k >= 2): Momentos normalizados divididos pela área |E|
        k_int = k - 2
        if (k_int >= 0) then
            if (.not. present(area) .or. .not. present(polygon_coords) .or. .not. present(n_verts)) then
                error stop "compute_D: Area e coordenadas do poligono sao obrigatorias para k >= 2"
            end if

            n_mon_int = ((k_int + 1) * (k_int + 2)) / 2
            allocate(p_int(n_mon_int), q_int(n_mon_int))
            call get_monomial_exponents(k_int, n_mon_int, p_int, q_int)

            do m_int = 1, n_mon_int
                do m = 1, n_mon
                    ! Integral dividida exatamente pela Área |E| (funcional do DOF interno)
                    int_mon_val = compute_monomial_domain_integral( &
                                p_exp(m) + p_int(m_int), q_exp(m) + q_int(m_int), &
                                xc, yc, h_E, polygon_coords, n_verts) / area

                    do d = 1, ndof
                        if (present(internal_dof_map)) then
                            row_idx = internal_dof_map(m_int, d)
                        else
                            row_idx = ndof * n_boundary_pts + ndof * (m_int - 1) + d
                        end if
                        
                        col_base = ndof * (m - 1) + d
                        
                        if (row_idx > 0 .and. row_idx <= loc_ndof) then
                            D_mat(row_idx, col_base) = int_mon_val
                        end if
                    end do
                end do
            end do

            deallocate(p_int, q_int)
        end if

        deallocate(p_exp, q_exp)
    end subroutine compute_matrix_D

    subroutine compute_matrix_G(G, k_eps, xc, yc, h_E, area, polygon_coords, n_verts)
        implicit none
        integer, intent(in)  :: k_eps, n_verts
        real(dp), intent(in) :: xc, yc, h_E, area, polygon_coords(n_verts, 2)
        real(dp), intent(out):: G(:,:)

        integer :: n_mon, m, j, r_base, c_base
        real(dp) :: int_qq 
        integer, allocatable :: p_exp(:), q_exp(:)

        n_mon = ((k_eps + 1) * (k_eps + 2)) / 2
        allocate(p_exp(n_mon), q_exp(n_mon))
        call get_monomial_exponents(k_eps, n_mon, p_exp, q_exp)

        G = 0.0_dp

        do m = 1, n_mon
            do j = 1, n_mon
                int_qq = compute_monomial_domain_integral(p_exp(m)+p_exp(j), q_exp(m)+q_exp(j), &
                                                       xc, yc, h_E, polygon_coords, n_verts)
            
                r_base = 3 * (m - 1)
                c_base = 3 * (j - 1)

                G(r_base + 1, c_base + 1) = int_qq ! eps_xx
                G(r_base + 2, c_base + 2) = int_qq ! eps_yy
                G(r_base + 3, c_base + 3) = int_qq ! gamma_xy
            end do
        end do

        deallocate(p_exp, q_exp)
    end subroutine compute_matrix_G

    subroutine apply_P0(B, k_order, n_verts, n_monomials, loc_dof, ndof, &
                    internal_dof_map, area)
        implicit none

        integer, intent(in) :: k_order, n_verts, n_monomials, loc_dof, ndof
        integer, intent(in), optional :: internal_dof_map(:,:)
        real(dp), intent(in), optional :: area
        real(dp), intent(inout) :: B(ndof*n_monomials, loc_dof)

        integer :: d, v, row_idx, col_idx

        do d = 1, ndof
            row_idx = d
            B(row_idx, :) = 0.0_dp

            if (k_order == 1) then
                ! k = 1: Média aritmética nos vértices
                do v = 1, n_verts
                    col_idx = ndof * (v - 1) + d
                    if (col_idx <= loc_dof) then
                        B(row_idx, col_idx) = 1.0_dp / real(n_verts, dp)
                    end if
                end do
            else
                ! k >= 2: O DOF interno de grau zero representa a média exata do deslocamento
                if (.not. present(internal_dof_map)) then
                    error stop "apply_P0: internal_dof_map eh obrigatorio para k >= 2"
                end if

                col_idx = internal_dof_map(1, d)
                if (col_idx > 0 .and. col_idx <= loc_dof) then
                    B(row_idx, col_idx) = 1.0_dp
                end if
            end if
        end do

    end subroutine apply_P0

    subroutine compute_matrix_Q(H0, Pi_nabla, k_order, n_monomials, ndof, loc_dof, &
                           n_internal, internal_dof_map, area, Q_mat)
        implicit none

        real(dp), intent(in)  :: H0(:,:), Pi_nabla(:,:)
        integer, intent(in)   :: k_order, n_monomials, ndof, loc_dof, n_internal
        integer, intent(in), optional :: internal_dof_map(:,:)
        real(dp), intent(in)  :: area
        real(dp), intent(out) :: Q_mat(ndof * n_monomials, loc_dof)

        integer :: i, d, row_idx, col_idx, n_mon_int

        ! 1. Graus superiores (alpha > n_{k-2}): utiliza a projeção cinemática
        Q_mat = matmul(H0, Pi_nabla)

        ! 2. Bloco Superior (alpha <= n_{k-2}): atribuição direta via integral do momento interno
        if (k_order >= 2 .and. n_internal > 0) then
            n_mon_int = ((k_order - 2 + 1) * (k_order - 2 + 2)) / 2
            do i = 1, n_mon_int
                do d = 1, ndof
                    row_idx = ndof * (i - 1) + d
                
                    Q_mat(row_idx, :) = 0.0_dp
                
                    if (present(internal_dof_map)) then
                        col_idx = internal_dof_map(i, d)
                        if (col_idx > 0 .and. col_idx <= loc_dof) then
                            Q_mat(row_idx, col_idx) = area
                        end if
                    end if
                end do
            end do
        end if

    end subroutine compute_matrix_Q

    subroutine compute_matrix_H0_exact(x, y, n, k, n_monomials, ndof, xc, yc, h_E, H0)
        implicit none
        integer, intent(in) :: n, k, n_monomials, ndof
        real(dp), intent(in) :: x(n), y(n), xc, yc, h_E
        real(dp), intent(out) :: H0(ndof*n_monomials, ndof*n_monomials)
        integer :: i, j, p, q, d, idx1, idx2
        real(dp) :: raw_moment
        integer, allocatable :: p_exp(:), q_exp(:)
        real(dp) :: x_shift(n), y_shift(n)

        H0 = 0.0_dp
        x_shift = x - xc
        y_shift = y - yc

        allocate(p_exp(n_monomials), q_exp(n_monomials))
        call get_monomial_exponents(k, n_monomials, p_exp, q_exp)

        do i = 1, n_monomials
            do j = 1, n_monomials
                p = p_exp(i) + p_exp(j)
                q = q_exp(i) + q_exp(j)
                call polygon_moment(x_shift, y_shift, n, p, q, raw_moment)
                raw_moment = raw_moment / (h_E**(p + q))

                do d = 1, ndof
                    idx1 = ndof * (i - 1) + d
                    idx2 = ndof * (j - 1) + d
                    H0(idx1, idx2) = raw_moment
                end do
            end do
        end do
        deallocate(p_exp, q_exp)
    end subroutine compute_matrix_H0_exact

    subroutine invert_and_project(G, B, Pi)
        implicit none
        real(dp), intent(in)  :: G(:,:), B(:,:)
        real(dp), intent(out) :: Pi(:,:)

        integer :: n_rows, n_cols, info
        integer, allocatable :: ipiv(:)
        real(dp), allocatable :: G_lu(:,:)
        real(dp) :: norm_G

        n_rows = size(G, 1)
        n_cols = size(B, 2)

        if (size(B, 1) /= n_rows) then
            error stop "Erro em invert_and_project: Numero de linhas de G e B sao incompativeis!"
        end if

        norm_G = maxval(abs(G))
        if (norm_G < 1.0e-14_dp) then
            error stop "Erro em invert_and_project: Matriz G eh nula ou degenerada (norma ~ 0)."
        end if

        allocate(G_lu(n_rows, n_rows), ipiv(n_rows))
        G_lu = G
        Pi = B

        call dgetrf(n_rows, n_rows, G_lu, n_rows, ipiv, info)
        if (info /= 0) then
            write(*, '(A, I0)') "Erro LAPACK DGETRF em invert_and_project. Info = ", info
            error stop "Erro na fatoracao LU de G em invert_and_project"
        end if

        call dgetrs('N', n_rows, n_cols, G_lu, n_rows, ipiv, Pi, n_rows, info)
        if (info /= 0) error stop "Erro na resolucao do sistema G * Pi = B"

        deallocate(G_lu, ipiv)
    end subroutine invert_and_project

end module vem_operators_mod