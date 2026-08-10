!===============================================================================
! MÓDULO PRINCIPAL DE ESTRUTURA E MONTAGEM VEM
!===============================================================================
module vem_core_mod
    use math_geometry_mod
    use vem_operators_mod
    use vem_concrete_operators_mod
    implicit none
    private

    public :: node_type, element_type, mesh_type
    public :: get_dof_maps, partition_matrix, partition_vector
    public :: assemble_matrix, assemble_vector
    public :: apply_dirichlet_and_springs
    public :: compute_vem_errors, write_solution
    public :: compute_vem_domain_trapezoidal_load
    public :: compute_vem_edge_trapezoidal_load
    public :: compute_matrix_H0_exact
    public :: compute_vem_stiffness_consistency
    public :: compute_vem_mass_consistency
    public :: compute_vem_stability
    public :: matrix_trace
    public :: dp



    type :: node_type
        real(dp) :: x, y
        integer  :: type_id
        logical, allocatable  :: is_fixed(:)
        real(dp), allocatable :: bc_val(:)
        real(dp), allocatable :: spring_k(:)
    end type node_type

    type :: element_type
        integer, allocatable :: nodes(:)
        integer, allocatable :: vertices(:)
        real(dp) :: area = 0.0_dp
        real(dp) :: diameter = 0.0_dp
        real(dp) :: centroid(2) = 0.0_dp
    end type element_type

    type :: mesh_type
        type(node_type), allocatable :: node(:)
        type(element_type), allocatable :: elem(:)
        integer :: nnodes = 0, nelem = 0, k_order = 1
        integer :: ndof = 2
    contains
        procedure :: read_mesh
        procedure :: set_dof_constraint
    end type mesh_type

contains

    subroutine read_next_valid_line(unit_num, line)
        integer, intent(in) :: unit_num
        character(len=*), intent(out) :: line
        integer :: ios
        do
            read(unit_num, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) > 0 .and. line(1:1) /= '!' .and. line(1:1) /= '#') exit
        end do
    end subroutine read_next_valid_line

    subroutine read_mesh(this, filename, n_dofs)
        class(mesh_type), intent(inout) :: this
        character(len=*), intent(in)    :: filename
        integer, intent(in), optional   :: n_dofs
        integer :: i, j, nnodes_local, ne_local, ios, n_elem_nodes, n_dir, n_spring
        integer :: node_idx, dof_idx
        real(dp) :: val
        character(len=500) :: io_err
        character(len=512) :: line

        open(newunit=i, file=filename, status='old', action='read', iostat=ios, iomsg=io_err)
        if (ios /= 0) error stop 'Erro ao abrir arquivo de malha: ' // trim(io_err)

        call read_next_valid_line(i, line)
        read(line, *, iostat=ios) nnodes_local, ne_local, this%k_order, this%ndof
        if (ios /= 0) error stop 'Erro na leitura do cabecalho da malha.'
        if (present(n_dofs)) this%ndof = n_dofs

        this%nnodes = nnodes_local
        this%nelem = ne_local
        allocate(this%node(nnodes_local), this%elem(ne_local))

        do j = 1, nnodes_local
            allocate(this%node(j)%is_fixed(this%ndof))
            allocate(this%node(j)%bc_val(this%ndof))
            allocate(this%node(j)%spring_k(this%ndof))

            this%node(j)%is_fixed = .false.
            this%node(j)%bc_val   = 0.0_dp
            this%node(j)%spring_k = 0.0_dp

            call read_next_valid_line(i, line)
            read(line, *) this%node(j)%x, this%node(j)%y, this%node(j)%type_id
        end do

        do j = 1, ne_local
            call read_next_valid_line(i, line)
            read(line, *) n_elem_nodes
    
            allocate(this%elem(j)%nodes(n_elem_nodes))
            allocate(this%elem(j)%vertices(n_elem_nodes)) ! Add allocation for vertices
    
            call read_next_valid_line(i, line)
            read(line, *) this%elem(j)%nodes
            this%elem(j)%vertices = this%elem(j)%nodes     ! Assign node IDs to vertices
        end do

        call read_next_valid_line(i, line)
        read(line, *, iostat=ios) n_dir
        if (ios == 0 .and. n_dir > 0) then
            do j = 1, n_dir
                call read_next_valid_line(i, line)
                read(line, *) node_idx, dof_idx, val
                if (node_idx >= 1 .and. node_idx <= this%nnodes .and. dof_idx <= this%ndof) then
                    this%node(node_idx)%is_fixed(dof_idx) = .true.
                    this%node(node_idx)%bc_val(dof_idx)   = val
                end if
            end do
        end if

        call read_next_valid_line(i, line)
        read(line, *, iostat=ios) n_spring
        if (ios == 0 .and. n_spring > 0) then
            do j = 1, n_spring
                call read_next_valid_line(i, line)
                read(line, *) node_idx, dof_idx, val
                if (node_idx >= 1 .and. node_idx <= this%nnodes .and. dof_idx <= this%ndof) then
                    this%node(node_idx)%spring_k(dof_idx) = val
                end if
            end do
        end if

        close(i)
    end subroutine read_mesh

    subroutine set_dof_constraint(this, node_id, dof_id, is_constrained, bc_val, spring_k)
        class(mesh_type), intent(inout) :: this
        integer, intent(in) :: node_id, dof_id
        logical, intent(in) :: is_constrained
        real(dp), intent(in), optional :: bc_val, spring_k

        if (node_id < 1 .or. node_id > this%nnodes) return
        this%node(node_id)%is_fixed(dof_id) = is_constrained
        if (present(bc_val))   this%node(node_id)%bc_val(dof_id)   = bc_val
        if (present(spring_k)) this%node(node_id)%spring_k(dof_id) = spring_k
    end subroutine set_dof_constraint

    pure function matrix_trace(A) result(tr)
        real(dp), intent(in) :: A(:,:)
        real(dp) :: tr
        integer :: i
        tr = 0.0_dp
        do i = 1, min(size(A, 1), size(A, 2))
            tr = tr + A(i, i)
        end do
    end function matrix_trace

    subroutine compute_vem_stiffness_consistency(Pi_eps, G_eps, C_mat, K_C)
        real(dp), intent(in)  :: Pi_eps(:,:), G_eps(:,:), C_mat(3,3)
        real(dp), intent(out) :: K_C(:,:)
        real(dp), allocatable :: C_G(:,:), C_expanded(:,:)
        integer :: n_mon, i, j, b

        n_mon = size(G_eps, 1) / 3
        allocate(C_expanded(3*n_mon, 3*n_mon), C_G(3*n_mon, 3*n_mon))
        C_expanded = 0.0_dp

        do b = 1, n_mon
            i = 3 * (b - 1) + 1
            j = 3 * (b - 1) + 3
            C_expanded(i:j, i:j) = C_mat
        end do

        C_G = matmul(C_expanded, G_eps)
        K_C = matmul(transpose(Pi_eps), matmul(C_G, Pi_eps))
        
        deallocate(C_expanded, C_G)
    end subroutine compute_vem_stiffness_consistency

    subroutine compute_vem_mass_consistency(Pi_0, H0, rho, M_C)
        real(dp), intent(in)  :: Pi_0(:,:), H0(:,:)
        real(dp), intent(in)  :: rho
        real(dp), intent(out) :: M_C(:,:)
        M_C = rho * matmul(transpose(Pi_0), matmul(H0, Pi_0))
    end subroutine compute_vem_mass_consistency

    subroutine compute_vem_stability(D, Pi, alpha_param, S_mat)
        real(dp), intent(in)  :: D(:,:), Pi(:,:)
        real(dp), intent(in)  :: alpha_param
        real(dp), intent(out) :: S_mat(:,:)
        integer :: ndof
        real(dp), allocatable :: I_mat(:,:), Proj(:,:)

        ndof = size(D, 1)
        allocate(I_mat(ndof, ndof), Proj(ndof, ndof))
        I_mat = eye(ndof)
        Proj  = I_mat - matmul(D, Pi)
        S_mat = alpha_param * matmul(transpose(Proj), Proj)
        deallocate(I_mat, Proj)
    end subroutine compute_vem_stability

    subroutine apply_dirichlet_and_springs(mesh, K_global, F_global, K_ff, F_mod, free_dofs)
        type(mesh_type), intent(in) :: mesh
        real(dp), intent(in)  :: K_global(:,:), F_global(:)
        real(dp), allocatable, intent(out) :: K_ff(:,:), F_mod(:)
        integer, allocatable, intent(out)  :: free_dofs(:)

        integer :: n_total, n_free, i, j, d, free_i, free_j, g_eq
        real(dp) :: f_val
        logical, allocatable :: is_fixed_vec(:)
        real(dp), allocatable :: bc_vec(:), spring_vec(:)

        n_total = mesh%nnodes * mesh%ndof
        allocate(is_fixed_vec(n_total), bc_vec(n_total), spring_vec(n_total))

        do i = 1, mesh%nnodes
            do d = 1, mesh%ndof
                g_eq = (i - 1) * mesh%ndof + d
                is_fixed_vec(g_eq) = mesh%node(i)%is_fixed(d)
                bc_vec(g_eq)       = mesh%node(i)%bc_val(d)
                spring_vec(g_eq)   = mesh%node(i)%spring_k(d)
            end do
        end do

        n_free = count(.not. is_fixed_vec)
        allocate(free_dofs(n_free))
        allocate(K_ff(n_free, n_free), F_mod(n_free))

        free_i = 0
        do i = 1, n_total
            if (.not. is_fixed_vec(i)) then
                free_i = free_i + 1
                free_dofs(free_i) = i
            end if
        end do

        do i = 1, n_free
            free_i = free_dofs(i)
            f_val = F_global(free_i)

            do j = 1, n_free
                free_j = free_dofs(j)
                if (free_i == free_j .and. spring_vec(free_i) > 0.0_dp) then
                    K_ff(i, j) = K_global(free_i, free_j) + spring_vec(free_i)
                else
                    K_ff(i, j) = K_global(free_i, free_j)
                end if
            end do

            do j = 1, n_total
                if (is_fixed_vec(j)) then
                    f_val = f_val - K_global(free_i, j) * bc_vec(j)
                end if
            end do
            F_mod(i) = f_val
        end do
        deallocate(is_fixed_vec, bc_vec, spring_vec)
    end subroutine apply_dirichlet_and_springs

    subroutine compute_matrix_H0_exact(x, y, n, k, xc, yc, h_E, n_monomials, ndof, H0)
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

    subroutine compute_vem_domain_trapezoidal_load(vert_x, vert_y, n_verts, k, ndof, &
                                                   q0, qx, qy, Pi_0, f_elem)
        integer, intent(in)  :: n_verts, k, ndof
        real(dp), intent(in) :: vert_x(n_verts), vert_y(n_verts)
        real(dp), intent(in) :: q0, qx, qy
        real(dp), intent(in) :: Pi_0(:,:)
        real(dp), intent(out):: f_elem(:)

        integer :: n_monomials, i, idx, d, p, q
        real(dp) :: area, xc, yc, h_E, mom
        real(dp), allocatable :: b_q(:), b_q_exp(:), q_coeff(:), H0(:,:)
        real(dp), allocatable :: x_shift(:), y_shift(:)
        integer, allocatable  :: p_exp(:), q_exp(:)

        n_monomials = ((k + 1) * (k + 2)) / 2
        f_elem = 0.0_dp

        call polygon_geometry(vert_x, vert_y, n_verts, area, xc, yc, h_E)

        allocate(q_coeff(n_monomials), b_q(n_monomials))
        q_coeff = 0.0_dp

        q_coeff(1) = q0
        if (k >= 1) then
            q_coeff(2) = qx
            q_coeff(3) = qy
        end if

        allocate(H0(n_monomials, n_monomials))
        allocate(p_exp(n_monomials), q_exp(n_monomials))
        allocate(x_shift(n_verts), y_shift(n_verts))

        x_shift = vert_x - xc
        y_shift = vert_y - yc

        call get_monomial_exponents(k, n_monomials, p_exp, q_exp)

        do i = 1, n_monomials
            do idx = 1, n_monomials
                p = p_exp(i) + p_exp(idx)
                q = q_exp(i) + q_exp(idx)
                call polygon_moment(x_shift, y_shift, n_verts, p, q, mom)
                H0(i, idx) = mom / (h_E**(p + q))
            end do
        end do

        b_q = matmul(H0, q_coeff)

        allocate(b_q_exp(ndof * n_monomials))
        do i = 1, n_monomials
            do d = 1, ndof
                b_q_exp(ndof*(i-1) + d) = b_q(i)
            end do
        end do

        f_elem = matmul(transpose(Pi_0), b_q_exp)
        deallocate(q_coeff, b_q, b_q_exp, H0, p_exp, q_exp, x_shift, y_shift)
    end subroutine compute_vem_domain_trapezoidal_load

    pure subroutine compute_vem_edge_trapezoidal_load(v1_coords, v2_coords, q1_vec, q2_vec, f_edge)
        real(dp), intent(in)  :: v1_coords(2), v2_coords(2)
        real(dp), intent(in)  :: q1_vec(2), q2_vec(2)
        real(dp), intent(out) :: f_edge(4)

        real(dp) :: dx, dy, L_edge

        dx = v2_coords(1) - v1_coords(1)
        dy = v2_coords(2) - v1_coords(2)
        L_edge = sqrt(dx**2 + dy**2)

        f_edge(1) = (L_edge / 6.0_dp) * (2.0_dp * q1_vec(1) + q2_vec(1))
        f_edge(2) = (L_edge / 6.0_dp) * (2.0_dp * q1_vec(2) + q2_vec(2))
        f_edge(3) = (L_edge / 6.0_dp) * (q1_vec(1) + 2.0_dp * q2_vec(1))
        f_edge(4) = (L_edge / 6.0_dp) * (q1_vec(2) + 2.0_dp * q2_vec(2))
    end subroutine compute_vem_edge_trapezoidal_load

    subroutine get_dof_maps(mesh, internal_dofs, boundary_dofs, N_i, N_b)
        type(mesh_type), intent(in) :: mesh
        integer, allocatable, intent(out) :: internal_dofs(:), boundary_dofs(:)
        integer, intent(out) :: N_i, N_b
        integer :: i, d, global_eq

        N_b = 0
        do i = 1, mesh%nnodes
            do d = 1, mesh%ndof
                if (mesh%node(i)%is_fixed(d)) N_b = N_b + 1
            end do
        end do
        N_i = (mesh%nnodes * mesh%ndof) - N_b

        allocate(internal_dofs(N_i), boundary_dofs(N_b))

        N_i = 0; N_b = 0
        do i = 1, mesh%nnodes
            do d = 1, mesh%ndof
                global_eq = (i - 1)*mesh%ndof + d
                if (mesh%node(i)%is_fixed(d)) then
                    N_b = N_b + 1
                    boundary_dofs(N_b) = global_eq
                else
                    N_i = N_i + 1
                    internal_dofs(N_i) = global_eq
                end if
            end do
        end do
    end subroutine get_dof_maps

    subroutine partition_matrix(A, internal_dofs, boundary_dofs, A_ii, A_ib, A_bi, A_bb)
        real(dp), intent(in) :: A(:,:)
        integer, intent(in)  :: internal_dofs(:), boundary_dofs(:)
        real(dp), allocatable, intent(out), optional :: A_ii(:,:), A_ib(:,:), A_bi(:,:), A_bb(:,:)
        integer :: N_i, N_b, i, j

        N_i = size(internal_dofs)
        N_b = size(boundary_dofs)

        if (present(A_ii) .and. N_i > 0) then
            allocate(A_ii(N_i, N_i))
            do i = 1, N_i
                do j = 1, N_i
                    A_ii(i,j) = A(internal_dofs(i), internal_dofs(j))
                end do
            end do
        end if

        if (present(A_ib) .and. N_i > 0 .and. N_b > 0) then
            allocate(A_ib(N_i, N_b))
            do i = 1, N_i
                do j = 1, N_b
                    A_ib(i,j) = A(internal_dofs(i), boundary_dofs(j))
                end do
            end do
        end if

        if (present(A_bi) .and. N_b > 0 .and. N_i > 0) then
            allocate(A_bi(N_b, N_i))
            do i = 1, N_b
                do j = 1, N_i
                    A_bi(i,j) = A(boundary_dofs(i), internal_dofs(j))
                end do
            end do
        end if

        if (present(A_bb) .and. N_b > 0) then
            allocate(A_bb(N_b, N_b))
            do i = 1, N_b
                do j = 1, N_b
                    A_bb(i,j) = A(boundary_dofs(i), boundary_dofs(j))
                end do
            end do
        end if
    end subroutine partition_matrix

    subroutine partition_vector(V, internal_dofs, boundary_dofs, V_i, V_b)
        real(dp), intent(in) :: V(:)
        integer, intent(in)  :: internal_dofs(:), boundary_dofs(:)
        real(dp), allocatable, intent(out), optional :: V_i(:), V_b(:)
        integer :: N_i, N_b, i

        N_i = size(internal_dofs)
        N_b = size(boundary_dofs)

        if (present(V_i) .and. N_i > 0) then
            allocate(V_i(N_i))
            do i = 1, N_i
                V_i(i) = V(internal_dofs(i))
            end do
        end if

        if (present(V_b) .and. N_b > 0) then
            allocate(V_b(N_b))
            do i = 1, N_b
                V_b(i) = V(boundary_dofs(i))
            end do
        end if
    end subroutine partition_vector

    subroutine assemble_matrix(mesh, eid, Mat_local, Mat_global)
        type(mesh_type), intent(in) :: mesh
        integer, intent(in) :: eid
        real(dp), intent(in) :: Mat_local(:,:)
        real(dp), intent(inout) :: Mat_global(:,:)
        integer :: i, j, nloc, ndof, eq_i, eq_j, loc_i, loc_j, d1, d2

        nloc = size(mesh%elem(eid)%nodes)
        ndof = mesh%ndof

        do i = 1, nloc
            do d1 = 1, ndof
                eq_i  = (mesh%elem(eid)%nodes(i) - 1)*ndof + d1
                loc_i = (i - 1)*ndof + d1
                do j = 1, nloc
                    do d2 = 1, ndof
                        eq_j  = (mesh%elem(eid)%nodes(j) - 1)*ndof + d2
                        loc_j = (j - 1)*ndof + d2
                        Mat_global(eq_i, eq_j) = Mat_global(eq_i, eq_j) + Mat_local(loc_i, loc_j)
                    end do
                end do
            end do
        end do
    end subroutine assemble_matrix

    subroutine assemble_vector(mesh, eid, Vec_local, Vec_global)
        type(mesh_type), intent(in) :: mesh
        integer, intent(in) :: eid
        real(dp), intent(in) :: Vec_local(:)
        real(dp), intent(inout) :: Vec_global(:)
        integer :: i, nloc, ndof, eq_i, loc_i, d1

        nloc = size(mesh%elem(eid)%nodes)
        ndof = mesh%ndof

        do i = 1, nloc
            do d1 = 1, ndof
                eq_i  = (mesh%elem(eid)%nodes(i) - 1)*ndof + d1
                loc_i = (i - 1)*ndof + d1
                Vec_global(eq_i) = Vec_global(eq_i) + Vec_local(loc_i)
            end do
        end do
    end subroutine assemble_vector

    subroutine compute_vem_errors(mesh, u_num, u_exact_vec, err_L2, err_H1)
        type(mesh_type), intent(in) :: mesh
        real(dp), intent(in) :: u_num(:), u_exact_vec(:)
        real(dp), intent(out) :: err_L2, err_H1
        real(dp) :: sum_sq
        integer :: i

        sum_sq = 0.0_dp
        do i = 1, size(u_num)
            sum_sq = sum_sq + (u_num(i) - u_exact_vec(i))**2
        end do

        err_L2 = sqrt(sum_sq / real(mesh%nnodes, dp))
        err_H1 = sqrt(sum_sq)
    end subroutine compute_vem_errors

    subroutine write_solution(mesh, u, filename)
        type(mesh_type), intent(in) :: mesh
        real(dp), intent(in) :: u(:)
        character(len=*), intent(in) :: filename
        integer :: i, d, unit_num
        character(len=500) :: io_err

        open(newunit=unit_num, file=filename, status='replace', action='write', iostat=i, iomsg=io_err)
        if (i /= 0) error stop 'Erro ao criar arquivo de saida: ' // trim(io_err)

        do i = 1, mesh%nnodes
            write(unit_num, '(2(es24.14), 10(es24.14))') mesh%node(i)%x, mesh%node(i)%y, &
                 (u((i-1)*mesh%ndof + d), d=1, mesh%ndof)
        end do
        close(unit_num)
    end subroutine write_solution



end module vem_core_mod